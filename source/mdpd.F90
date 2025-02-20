!> Module to calculate local densities and resulting force and potential 
!> contributions for many-body DPD interactions
Module mdpd
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly module for calculating many-body dpd forces and potentials
  !
  ! Note: partially based on the mdpd module found in DL_MESO writen by 
  !       w. smith, m. a. seaton & r. l. anderson april 2022
  !
  ! copyright - ukri stfc daresbury laboratory
  ! authors   - b.t.speake 2024
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,                Only: MetLdExp_tag, comms_type, gcheck, girecv, &
                                  gmax, gsend, gsum, gwait 
  Use configuration,        Only: configuration_type 
  Use constants,            Only: pi 
  Use domains,              Only: domains_type
  Use errors_warnings,      Only: error_alloc, error, warning
  Use kinds,                Only: wp, wi, STR_LEN
  Use neighbours,           Only: neighbours_type
  Use two_body_potentials,  Only: potential_energy
#ifdef HALF_HALO
  Use numerics,        Only: local_index
#endif /* HALF_HALO */

  Implicit None 

  Private 

  !> 1 / pi 
  Real(Kind=wp), Parameter :: rpi = 1.0_wp / pi 

  Public  :: mdpd_ld_compute, mdpd_force, weight_rho, mdpd_potential

  !> Type to contain the required data arrays for many-body DPD interactions
  Type, Public :: mdpd_type 
    !> mdpd b parameter (max of all specified in FIELD file)
    Real(Kind=wp), Public                 :: b 
    !> rc^3 value for use in force evaulation 
    Real(Kind=wp), Allocatable, Public    :: rc(:)
    !> mdpd per interaction rd parameters 
    Real(Kind=wp), Allocatable, Public    :: rd(:)
    !> mdpd per interaction n parameters  
    Real(Kind=wp), Allocatable, Public    :: n(:)
    !> mdpd per site m parameters 
    Real(Kind=wp), Allocatable, Public    :: m(:)
    !> mdpd per site local densities 
    Real(Kind=wp), Allocatable, Public    :: rho(:)

  Contains 
    Private
    Procedure, Public :: allocate => allocate_mdpd_arrays 
    Procedure, Public :: set_potential_params => set_potential_params
  End Type mdpd_type

Contains 

  !> Allocates the required arrays for mDPD interactions 
  Subroutine allocate_mdpd_arrays(mdpd_params, nvdw, natms)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 
    ! Allocation of the required resources for storing the additional 
    ! many-body DPD parameters.
    !
    ! copyright - daresbury laboratory
    ! authors     - b.t.speake July 2024 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(mdpd_type)                        :: mdpd_params 
    Integer,                  Intent(In   ) :: nvdw, natms 
    Integer, Dimension(1:4)                 :: fail 

    fail = 0 

    Allocate(mdpd_params%rc(nvdw), Stat=fail(1))
    Allocate(mdpd_params%rd(nvdw), Stat=fail(2))
    Allocate(mdpd_params%n(nvdw), Stat=fail(3))
    Allocate(mdpd_params%m(natms), Stat=fail(4))
    ! Allocate(mdpd%rho(natms), Stat=fail(4))

    If (Any(fail > 0)) Call error_alloc('mdpd arrays', 'allocate_mdpd_arrays')

    mdpd_params%rd = 0.0_wp 
    mdpd_params%n = 2.0_wp
    mdpd_params%rc = 0.0_wp 

  End Subroutine allocate_mdpd_arrays 

  !> Stores the requires many-body dpd per pairwise potential paramters 
  !> including pre-calculation of r_{c}^{3} to improce performace.
  Subroutine set_potential_params(mdpd_params, k, rc, rd, n)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 
    ! Stores the requires many-body dpd per pairwise potential paramters 
    ! including pre-calculation of r_{c}^{3} to improce performace. 
    !
    ! copyright - daresbury laboratory
    ! authors     - b.t.speake July 2024 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(mdpd_type)                        :: mdpd_params
    Real(Kind=wp), Intent(In   )            :: rc, rd, n
    Integer(Kind=wi), Intent(In   )         :: k

    mdpd_params%rd(k) = rd 
    mdpd_params%n(k) = n 
    mdpd_params%rc(k) = rc * rc * rc 

  End Subroutine set_potential_params

  !> Computes the local densities required for many-body DPD interactions
  Subroutine mdpd_ld_compute(mdpd, neigh, domain, config, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dl_poly subroutine for calculating local density for manybody-
    ! dissipative particle dynamics. 
    ! 
    ! copyright - daresbury laboratory
    ! author    - b.t.speake 2024 
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(mdpd_type),                 Intent(InOut) :: mdpd
    Type(neighbours_type),          Intent(In   ) :: neigh
    Type(domains_type),             Intent(In   ) :: domain
    Type(configuration_type),       Intent(InOut) :: config
    Type(comms_type),               Intent(InOut) :: comm

    Integer(Kind=wi) :: i, j, k, key, limit, ai, aj
    Real(Kind=wp) :: rrr, xxt, yyt, zzt, rho
    
    ! move this elsewhere 
    If (.not. Allocated(mdpd%rho)) Then 
      Allocate(mdpd%rho(config%mxatms))
    End If 

    mdpd%rho = 0.0_wp
    Do i = 1, config%natms 
      limit = neigh%list(0, i) ! Get list limit

      ! calculate interatomic distances
      Do k = 1, limit
        j = neigh%list(k, i)

        xxt = config%parts(i)%xxx - config%parts(j)%xxx
        yyt = config%parts(i)%yyy - config%parts(j)%yyy
        zzt = config%parts(i)%zzz - config%parts(j)%zzz
        rrr = Sqrt(xxt**2 + yyt**2 + zzt**2)

        ai = config%ltype(i)
        aj = config%ltype(j)

        If (ai > aj) Then 
          key = ai * (ai - 1) / 2 + aj
        Else 
          key = aj * (aj - 1) / 2 + ai
        End If

        rho = weight_rho(rrr, mdpd%rc(key), mdpd%rd(key), mdpd%n(key))
        mdpd%rho(i) = mdpd%rho(i) + rho  
#ifdef HALF_HALO 
        mdpd%rho(j) = mdpd%rho(j) + rho 
#else 
        If (j <= config%natms) mdpd%rho(j) = mdpd%rho(j) + rho
#endif 
      End Do
    End Do 
    
#ifdef HALF_HALO
    Call half_halo_refresh(domain, config, mdpd, comm)
#endif 

    ! Update the halo 
    Call refresh_halo_mdpd_densities(domain, config, mdpd, comm)

  End Subroutine mdpd_ld_compute

  !> Updates the local density data within the halos 
  Subroutine refresh_halo_mdpd_densities(domain, config, mdpd, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dl_poly subroutine for updating the mdpd local densities within the 
    ! halo for the domain decomposition parallelisation. 
    ! 
    ! copyright - daresbury laboratory
    ! author    - b.t.speake 2024 
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(mdpd_type),                 Intent(InOut) :: mdpd
    Type(domains_type),             Intent(In   ) :: domain
    Type(configuration_type),       Intent(InOut) :: config
    Type(comms_type),               Intent(InOut) :: comm

    Character(Len=STR_LEN)   :: message
    Integer              :: fail, mlast
    Integer, Allocatable :: ixyz0(:)
    Logical              :: safe

    fail = 0
    Allocate (ixyz0(1:config%mxatms), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'mdpd_ld_set_halo allocation failure'
      Call error(0, message)
    End If
    ixyz0(1:config%nlast) = config%ixyz(1:config%nlast)

    ! No halo, start with domain only particles
    mlast = config%natms

    ! exchange atom data in -/+ x directions
    Call mdpd_ld_export(-1, mlast, config%mxatms, ixyz0, mdpd, domain, comm)
    Call mdpd_ld_export(1, mlast, config%mxatms, ixyz0, mdpd, domain, comm)

    ! exchange atom data in -/+ y directions
    Call mdpd_ld_export(-2, mlast, config%mxatms, ixyz0, mdpd, domain, comm)
    Call mdpd_ld_export(2, mlast, config%mxatms, ixyz0, mdpd, domain, comm)

    ! exchange atom data in -/+ z directions
    Call mdpd_ld_export(-3, mlast, config%mxatms, ixyz0, mdpd, domain, comm)
    Call mdpd_ld_export(3, mlast, config%mxatms, ixyz0, mdpd, domain, comm)

    ! check atom totals after data transfer
    safe = (mlast == config%nlast)
    Call gcheck(comm, safe)
    If (.not. safe) Then 
      Write(message, '(a)') '# error - incorrect atom totals assignments in refresh_halo_mdpd_densities'
      Call error(0)
    End If  

    Deallocate (ixyz0, Stat=fail)

  End Subroutine refresh_halo_mdpd_densities

  !> Passes local density data within a halo to the adjoining domain in a given 
  !> direction
  Subroutine mdpd_ld_export(mdir, mlast, mxatms, ixyz0, mdpd, domain, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dl_poly subroutine for updating the mdpd local densities within the 
    ! halo for the domain decomposition parallelisation.n Exports a specific 
    ! set of densities to the neighbouring domains.
    ! 
    ! copyright - daresbury laboratory
    ! author    - b.t.speake 2024 
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),               Intent(In   ) :: mdir, mxatms
    Integer(Kind=wi),               intent(InOut) :: mlast 
    Integer(Kind=wi),               Intent(InOut) :: ixyz0(:) 
    Type(mdpd_type),                Intent(InOut) :: mdpd
    Type(domains_type),             Intent(In   ) :: domain
    Type(comms_type),               Intent(InOut) :: comm

    Character(Len=STR_LEN)                        :: message
    Integer                                       :: kx, ky, kz, jxyz, kxyz, jdnode, kdnode 
    Integer                                       :: i, j, imove, jmove, fail, iblock, limit
    Integer                                       :: ix, iy, iz, iadd, itmp 
    Logical                                       :: safe
    Real(Kind=wp), Allocatable, Dimension(:)      :: buffer 

    iadd = 2
    ! Should send two values acrros, rho and m -> atm just sends rho
    fail = 0; limit = iadd * domain%mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'mdpd_ld_export allocation failure'
      Call error(0, message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)
    iblock = limit / Merge(2, 1, comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! kxyz - corrected halo reduction factor particles haloing both +&- sides
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0; ky = 0; kz = 0
    If (mdir == -1) Then ! Direction -x
      kx = 1
      jxyz = 1
      kxyz = 3

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir == 1) Then ! Direction +x
      kx = 1
      jxyz = 2
      kxyz = 3

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky = 1
      jxyz = 10
      kxyz = 30

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir == 2) Then ! Direction +y
      ky = 1
      jxyz = 20
      kxyz = 30

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz = 1
      jxyz = 100
      kxyz = 300

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir == 3) Then ! Direction +z
      kz = 1
      jxyz = 200
      kxyz = 300

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(47)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed
    imove = 0 
    jmove = 0 
    ! Initialise array overflow flags
    safe = .true.

    ! loop over particles on node 
    Do i = 1, mlast 
      ! If particle is within the 'inverted halo' of this domain 
      If (ixyz0(i) > 0) Then 
        ! Get halo indices 
        ix = Mod(ixyz0(i), 10) 
        iy = Mod(ixyz0(i) - ix, 100)
        iz = Mod(ixyz0(i) - (ix + iy), 1000)
        ! Filter by move direction 
        j = ix * kx + iy * ky + iz * kz 
        If (j == jxyz .or. (j > jxyz .and. Mod(j, 3) == 0)) Then 
          If ((imove + iadd) <= iblock) Then 
            ! pack halo index and particle density 
            buffer(imove + 1) = mdpd%rho(i)
            buffer(imove + 2) = Real(ixyz0(i) - Merge(jxyz, kxyz, j == jxyz), wp)
          Else 
            safe = .false.
          End If 
          imove = imove + iadd 
        End If 
      End If 
    End Do 

    ! Check for array bound overflow 
    Call gcheck(comm, safe) 
    If (.not. safe) Then 
      itmp = Merge(2, 1, comm%mxnode > 1) * imove 
      Call gmax(comm, itmp) 
      Call warning(150, Real(itmp, wp), Real(limit, wp), 0.0_wp)
      Write(message, '(a)') '# error - outgoing transfer buffer size exceeded in mdpd_ld_export'
      Call error(0, message)
    End If 

    ! TODO -> change tag...
    ! Exchange information on buffer size 
    If (comm%mxnode > 1) Then 
      Call girecv(comm, jmove, kdnode, MetLdExp_tag)
      Call gsend(comm, imove, jdnode, MetLdExp_tag)
      Call gwait(comm)
    Else 
      jmove = imove
    End If 

    ! Check for array bound overflow 
    safe = ((mlast + jmove / iadd) <= mxatms)
    Call gcheck(comm, safe)
    If (.not. safe) Then
      itmp = mlast + jmove / iadd
      Call gmax(comm, itmp)
      Call warning(160, Real(itmp, wp), Real(mxatms, wp), 0.0_wp)
      Write(message, '(a)') '# error - incoming data transfer size exceeds limit in mdpd_ld_export'
      Call error(0, message)
    End If

    ! exchange buffers 
    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm, buffer(iblock + 1:iblock + jmove), kdnode, MetLdExp_tag)
      End If
      If (imove > 0) Then
        Call gsend(comm, buffer(1:imove), jdnode, MetLdExp_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data 
    j = Merge(iblock, 0, comm%mxnode > 1) 
    Do i = 1, jmove / iadd 
      mlast = mlast + 1 
      mdpd%rho(mlast) = buffer(j + 1)
      ixyz0(mlast) = Nint(buffer(j + 2))
      j = j + iadd 
    End Do 

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'mdpd_ld_export deallocation failure'
      Call error(0, message)
    End If

  End Subroutine mdpd_ld_export

  !> Calculates the many-body DPD force contribution based on the local densities
  Pure Real(Kind=wp) Function mdpd_force(r, rc, rd, b, mi, mj, n, rhoi, rhoj)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dl_poly function for calculating manybody-DPD force corrections 
    ! based on the local particle densities
    ! 
    ! copyright - daresbury laboratory
    ! author    - b.t.speake 2024 
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp),          Intent(In   ) :: r, rc, rd, b, mi, mj, n, rhoi, rhoj
    Real(Kind=wp)                         :: rpmi, rpmj, rrd, rrdrd

    if (r < rd) Then 
      rrd = 1 / rd
      ! rrdrd = rrd * rrd
      rrdrd = rrd * rrd 
      rrdrd = rrdrd * rrdrd 
      if (n > 2) Then 
        rpmi = rhoi**(mi - 1) ! <- use of m parameter is not yet available 
        rpmj = rhoj**(mj - 1) 
        mdpd_force = b * n * (n + 1) * (n + 2) * (n + 3) * 0.125_wp * rpi 
        mdpd_force = mdpd_force * rrdrd * rrdrd
        mdpd_force = mdpd_force * (rpmi + rpmj) * (1 - (r * rrd))**(n - 1) 
        mdpd_force = mdpd_force * rrd
      Else 
        mdpd_force = b * 15 * rpi  * rrdrd * rc
        mdpd_force = mdpd_force * (rhoi + rhoj) * (1 / r) * (1 - (r * rrd))
      End If 
    Else 
      mdpd_force = 0.0_wp 
    End If 
    
    Return 
  End Function mdpd_force

  !> Evaluates the many-body DPD weighting function for local density evaluation 
  !>
  !> \begin{equation}
  !>  w^{\rho}(r_{ij}) = \frac{(n_{ij} + 1)(n_{ij} + 2)(n_{ij} + 3)r^{3}_{c}}{8\pi r^{3}_{d,ij}} 
  !> \left(1 - \frac{r_{ij}}{r_{d,ij}}\right)^{n_{ij}},
  !> \end{equation}
  Pure Real(kind=wp) Function weight_rho(rrr, rc,  rd, n) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nomalised weight function required to calculate local densities for 
    ! many-body dpd 
    ! 
    ! copyright - daresbury laboratory
    ! author    - b.t.speake 2024 
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rrr, rc, rd, n 
    Real(Kind=wp) :: rd2, rrd, nrm, scrn 

    rrd = 1.0_wp / rd 
    rd2 = rrd * rrd 
    scrn = MAX(rd - rrr, 0.0_wp) * rrd 
    
    If (n == 2) Then 
      nrm = 7.5_wp * rpi * rrd * rrd * rrd 
      weight_rho = nrm * scrn * scrn * rc
    Else 
      nrm = 0.125_wp * (n + 1.0_wp) * (n + 2.0_wp) + (n + 3.0_wp) * rpi * rrd * rrd * rrd
      weight_rho = rc* nrm * scrn ** n  
    End If 
    
    Return 
  End Function weight_rho 

  !> Evaluates the many-body self-energy contribution resulting from the 
  !> local density dependance. 
  !> \begin{equation}
  !> U^{mb} = \sum_i \frac{B \tilde{\rho}_i^{m_{c(i)}}{m_{c(i)}}
  !> \end{equation}
  Pure Real(Kind=wp) Function mdpd_potential(mdpd_params, config) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Evaluates the additional density dependent mdpd potential term. 
    ! 
    ! copyright - daresbury laboratory
    ! author    - b.t.speake 2024 
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(mdpd_type), Intent(In) :: mdpd_params
    Type(configuration_type), Intent(In) :: config 
    Integer(Kind=wi) :: i 
    Real(Kind=wp) :: m, rm

    mdpd_potential = 0.0_wp 
    Do i = 1, config%natms
      m = mdpd_params%m(config%ltype(i))
      rm = 1 / m 
      mdpd_potential = mdpd_potential + (rm * mdpd_params%rho(i)**m)
    End Do 
    mdpd_potential = mdpd_potential * (mdpd_params%b)
    Return 
  End Function mdpd_potential

#ifdef HALF_HALO 
!> Updates the local density data within the halos (Half-Halo)
  Subroutine half_halo_refresh(domain, config, mdpd, comm) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export mdpd densities in domain boundary regions
    ! for half-halo refresh 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(domains_type), Intent(In) :: domain 
    Type(configuration_type), Intent(In) :: config 
    Type(mdpd_type), Intent(InOut) :: mdpd 
    Type(comms_type), Intent(InOut) :: comm 

    Call half_halo_export(  3, domain, config, mdpd, comm) ! x-, y0, z0
    ! Call half_halo_export (-3 ,domain ,config ,mdpd ,comm) ! x+, y0, z0
    Call half_halo_export(  2, domain, config, mdpd, comm) ! x-, y0, z0
    Call half_halo_export( -2, domain, config, mdpd, comm) ! x+, y0, z0
    Call half_halo_export(  1, domain, config, mdpd, comm) ! x-, y0, z0
    Call half_halo_export( -1, domain, config, mdpd, comm) ! x+, y0, z0


  End Subroutine half_halo_refresh
  
  !> Passes local density data within a halo to the adjoining domain in a given 
  !> direction (Half-Halo). 
  Subroutine half_halo_export(mdir, domain, config, mdpd, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export mdpd densities in domain boundary regions
    ! for half-halo refresh 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer,                  Intent(In   ) :: mdir
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(In   ) :: config
    Type(mdpd_type),          Intent(InOut) :: mdpd
    Type(comms_type),         Intent(InOut) :: comm

    Character ( Len = 256 )  :: message
    Logical           :: lrhs, safe
    Integer           :: fail,iadd,limit,iblock,npdirB,npdirE,mlast,i,j, &
                         jdnode,kdnode,imove,jmove,itmp
    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    iadd = 2 
    fail = 0
    limit = iadd * domain%mxbfxp
    Allocate(buffer(1:limit), Stat=fail)
    If (fail > 0) Then 
      Write(message, '(a)') "half_halo_export for mdpd densities buffer allocation failure."
      Call error(0, message)
    End If 
    
    iblock = limit / Merge(2, 1, comm%mxnode > 1)

        ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with respect to the direction (mdir)
    ! jdnode - destination (send to), kdnode - source (receive from)

    npdirB = 0
    npdirE = 0

    If      (mdir == -1) Then ! Direction -/+x, y0, z0
      jdnode = domain%map(1) ! -x
      kdnode = domain%map(2) ! +x

      npdirB = config%ixyzM(1)+1
      npdirE = config%ixyzM(2)
    Else If (mdir ==  1) Then ! Direction +/-x, y0, z0
      jdnode = domain%map(2) ! +x
      kdnode = domain%map(1) ! -x

      npdirB = config%ixyzM(0)+1
      npdirE = config%ixyzM(1)
    Else If (mdir == -2) Then ! Direction x0, -/+y, z0
      jdnode = domain%map(3) ! -y
      kdnode = domain%map(4) ! +y

      npdirB = config%ixyzM(3)+1
      npdirE = config%ixyzM(4)
    Else If (mdir ==  2) Then ! Direction x0, +/-y, z0
      jdnode = domain%map(4) ! +y
      kdnode = domain%map(3) ! -y

      npdirB = config%ixyzM(2)+1
      npdirE = config%ixyzM(3)
    Else If (mdir == -3) Then ! Direction x0, y0, -/+z
      jdnode = domain%map(5) ! -z
      kdnode = domain%map(6) ! +z

      npdirB = config%ixyzM(5)+1
      npdirE = config%ixyzM(6)
    Else If (mdir ==  3) Then ! Direction x0, y0, +z
      jdnode = domain%map(6) ! +z
      kdnode = domain%map(5) ! -z

      npdirB = config%ixyzM(4)+1
      npdirE = config%ixyzM(5)
    Else
      Write(message, '(a)') '# error - undefined direction passed to mdpd half-halo refresh'
      Call error(0, message)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed
    imove=0
    jmove=0

    ! Initialise array overflow flags
    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE
    If (imove+iadd*(npdirE-npdirB) <= iblock) Then
      Do i=npdirB,npdirE
          buffer(imove+1) = Real(config%ltg(i),wp)
          buffer(imove+2) = mdpd%rho(i)

          imove=imove+iadd
      End Do
    Else
      safe=.false.
    End If

    ! Check for array bound overflow 
    Call gcheck(comm, safe)
    If (.not. safe) Then 
      itmp = Merge(2, 1, comm%mxnode > 1) * imove 
      Call gmax(comm, itmp)
      Call warning(150, Real(itmp, wp), Real(limit, wp), 0.0_wp)
      Write(message, '(a)') '# error - outgoing transfer buffer size exceeded in mdpd half-halo refresh'
      Call error(0, message)
    End If 

    ! TODO:: change tag 
    ! exchange information on buffer sizes
    If (comm%mxnode > 1) Then
      Call girecv(comm, jmove, kdnode, MetLdExp_tag-1)
      Call gsend(comm, imove, jdnode, MetLdExp_tag-1)
      Call gwait(comm)
    Else
      jmove = imove
    End If

    ! exchange buffers between nodes (this is a MUST)
    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm, buffer(iblock+1:iblock+jmove), kdnode, MetLdExp_tag)
      End If
      If (imove > 0 ) Then
        Call gsend(comm, buffer(1:imove), jdnode, MetLdExp_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data
    j = Merge(iblock, 0, comm%mxnode > 1)
    Do i=1, jmove/iadd
      mlast = local_index(Int(buffer(j+1)), config%nlast, config%lsi, config%lsa)
      mdpd%rho(mlast) = mdpd%rho(mlast) + buffer(j+2)
      j = j + iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'half_halo_refresh for mdpd densities buffer deallocation failure'
      Call error(0, message)
    End If

  End Subroutine half_halo_export
#endif 

End Module mdpd