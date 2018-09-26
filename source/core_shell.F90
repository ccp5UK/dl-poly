Module core_shell

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for defining global core-shell interaction variables
  ! and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov december 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use comms,           Only : comms_type,gsync,gsum,gcheck,gmax
  Use configuration,   Only : configuration_type,freeze_atoms
  Use particle,        Only : corePart
  Use setup,           Only : nrite,boltz,engunit,output,mxatms,mxatdm,zero_plus,&
    mxtmls,mxlshp
  Use parse,           Only : strip_blanks,lower_case
  Use shared_units,    Only : update_shared_units, SHARED_UNIT_UPDATE_FORCES
  Use numerics,        Only : local_index,images
  Use errors_warnings, Only : error,warning,info
  Use statistics, Only : stats_type
  Use domains, Only : domains_type

  Implicit None
  Private
  Integer, Parameter, Public :: SHELL_ADIABATIC = 1
  InTeger, Parameter, Public :: SHELL_RELAXED = 2
  Integer, Parameter :: NO_SEARCH = 0
  Integer, Parameter :: LINE_SEARCH = 1
  Integer, Parameter :: CONJUGATE_SEARCH = 2

  Type, Public :: core_shell_type
    Private
    Integer, Public :: mxshl,mxtshl,mxfshl,megshl,mxlshp
    Logical,                        Public :: lshmv_shl = .false.

    Integer,                        Public :: ntshl  = 0 , &
      ntshl1 = 0 , &
      ntshl2 = 0 , keyshl

    Real( Kind = wp ),              Public :: smax = 0.0_wp


    Integer,           Allocatable, Public :: numshl(:)
    Integer,           Allocatable, Public :: lstshl(:,:),listshl(:,:),legshl(:,:)
    Integer,           Allocatable, Public :: lishp_shl(:),lashp_shl(:)

    Real( Kind = wp ), Allocatable, Public :: prmshl(:,:)
    Logical            :: newjob = .true. , l_rdf
    Integer            :: keyopt
    Real( Kind = wp )  :: grad_tol,eng_tol,dist_tol(1:2),   &
      step,eng,eng0,eng1,eng2,          &
      grad,grad0,grad1,grad2,onorm,sgn, &
      stride,gamma,x(1),y(1),z(1),fff(0:3)
    Real( Kind = wp )  :: grad_pass
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)
  Contains
    Private
    Procedure, Public :: init => allocate_core_shell_arrays
    Procedure, Public :: deallocate_core_shell_tmp_arrays
    Final :: deallocate_core_shell_arrays
  End Type core_shell_type

  Public :: core_shell_forces
  Public :: core_shell_kinetic
  Public :: core_shell_on_top
  Public :: core_shell_quench
  Public :: core_shell_relax

Contains

  Subroutine allocate_core_shell_arrays(T,mxatdm,mxtmls,mxlshp,neighbours)
   Class( core_shell_type ) :: T
     Integer, Intent(In) :: mxatdm,mxtmls,mxlshp,neighbours

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (T%numshl(1:mxtmls),                        Stat = fail(1))
    Allocate (T%lstshl(1:2,1:T%mxtshl),                    Stat = fail(2))
    Allocate (T%listshl(0:2,1:T%mxshl),                    Stat = fail(3))
    Allocate (T%legshl(0:T%mxfshl,1:mxatdm),               Stat = fail(4))
    Allocate (T%lishp_shl(1:mxlshp),T%lashp_shl(1:neighbours), Stat = fail(5))
    Allocate (T%prmshl(1:2,1:T%mxtshl),                    Stat = fail(6))

    If (Any(fail > 0)) Call error(1005)

    T%numshl  = 0
    T%lstshl  = 0
    T%listshl = 0
    T%legshl  = 0

    T%lishp_shl = 0 ; T%lashp_shl = 0

    T%prmshl  = 0.0_wp

  End Subroutine allocate_core_shell_arrays

  Subroutine deallocate_core_shell_tmp_arrays(T)
   Class( core_shell_type ) :: T

    Integer :: fail(2)

    fail = 0

    If (Allocated(T%numshl)) Deallocate (T%numshl, Stat = fail(1))
    If (Allocated(T%lstshl)) Deallocate (T%lstshl, Stat = fail(2))

    If (Any(fail > 0)) Call error(1030)

  End Subroutine deallocate_core_shell_tmp_arrays

  Subroutine deallocate_core_shell_arrays(T)
    Type( core_shell_type ) :: T

    Integer :: fail(7)

    fail = 0

    If (Allocated(T%numshl)) Deallocate (T%numshl, Stat = fail(1))
    If (Allocated(T%lstshl)) Deallocate (T%lstshl, Stat = fail(2))
    If (Allocated(T%listshl)) Deallocate (T%listshl, Stat = fail(3))
    If (Allocated(T%legshl)) Deallocate (T%legshl, Stat = fail(4))
    If (Allocated(T%lishp_shl)) Deallocate (T%lishp_shl, Stat = fail(5))
    If (Allocated(T%lashp_shl)) Deallocate (T%lashp_shl, Stat = fail(6))
    If (Allocated(T%prmshl)) Deallocate (T%prmshl, Stat = fail(7))
    If (Any(fail > 0)) Call error(1030)
  End Subroutine deallocate_core_shell_arrays


  Subroutine core_shell_forces(cshell,stat,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating core-shell model spring energy
    ! and force terms
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( stats_type ), Intent( InOut ) :: stat
    Type( comms_type ),                  Intent( InOut ) :: comm
    Type( configuration_type ),          Intent( InOut ) :: config

    Logical           :: safe
    Integer           :: fail(1:2),i,j,ia,ib,kk
    Real( Kind = wp ) :: rabsq,fx,fy,fz,gamma,omega,r_4_fac, &
      strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

    Logical,           Allocatable :: lunsafe(:)
    Integer,           Allocatable :: lstopt(:,:)
    Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
    Character( Len = 256 ) :: message

    fail=0
    Allocate (lunsafe(1:cshell%mxshl),lstopt(0:2,1:cshell%mxshl),      Stat=fail(1))
    Allocate (xdab(1:cshell%mxshl),ydab(1:cshell%mxshl),zdab(1:cshell%mxshl), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_forces allocation failure'
      Call error(0,message)
    End If

    r_4_fac = 1.0_wp/24.0_wp ! aharmonic shell coefficient = 1/(4!)

    ! calculate core-shell separation vectors

    Do i=1,cshell%ntshl
      lunsafe(i)=.false.

      ! indices of atoms in a core-shell

      ia=local_index(cshell%listshl(1,i),config%nlast,config%lsi,config%lsa) ; lstopt(1,i)=ia
      ib=local_index(cshell%listshl(2,i),config%nlast,config%lsi,config%lsa) ; lstopt(2,i)=ib

      lstopt(0,i)=0
      If (ia > 0 .and. ib > 0) Then ! Tag
        If (ia <= config%natms .or. ib <= config%natms) Then
          lstopt(0,i)=1
        End If
      Else                          ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= config%natms) .or.   &
          (ib > 0 .and. ib <= config%natms)) .and. &
          (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
      End If

      ! components of bond vector

      If (lstopt(0,i) > 0) Then
        xdab(i)=config%parts(ia)%xxx-config%parts(ib)%xxx
        ydab(i)=config%parts(ia)%yyy-config%parts(ib)%yyy
        zdab(i)=config%parts(ia)%zzz-config%parts(ib)%zzz
      Else ! (DEBUG)
        xdab(i)=0.0_wp
        ydab(i)=0.0_wp
        zdab(i)=0.0_wp
      End If
    End Do

    ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:cshell%ntshl))
    Call gcheck(comm,safe)
    If (.not.safe) Then
      Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
          Do i=1,cshell%ntshl
            If (lunsafe(i)) Then
              Write(message,'(a,2(i10,a))')     &
                'global unit number', cshell%listshl(0,i), &
                ' , with a head particle number', cshell%listshl(1,i),   &
                ' contributes towards next error'
              Call warning(message)
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(100)
    End If

    ! periodic boundary condition

    Call images(config%imcon,config%cell,cshell%ntshl,xdab,ydab,zdab)

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! zero core-shell energy and virial accumulators

    stat%engshl=0.0_wp
    stat%virshl=0.0_wp

    ! loop over all specified core-shell units

    Do i=1,cshell%ntshl
      If (lstopt(0,i) > 0) Then

        ! indices of atoms in a core-shell

        ia=lstopt(1,i)
        ib=lstopt(2,i)

        ! define components of bond vector

        rabsq = xdab(i)**2+ydab(i)**2+zdab(i)**2

        ! index of potential function parameters

        kk=cshell%listshl(0,i)

        ! calculate scalar constant terms using spring potential function
        ! and the parameters in array prmshl

        omega=(0.5_wp*cshell%prmshl(1,kk)+r_4_fac*cshell%prmshl(2,kk)*rabsq)*rabsq
        gamma=cshell%prmshl(1,kk)+cshell%prmshl(2,kk)*rabsq

        ! calculate forces

        fx = -gamma*xdab(i)
        fy = -gamma*ydab(i)
        fz = -gamma*zdab(i)

        If (ia <= config%natms) Then

          config%parts(ia)%fxx=config%parts(ia)%fxx+fx
          config%parts(ia)%fyy=config%parts(ia)%fyy+fy
          config%parts(ia)%fzz=config%parts(ia)%fzz+fz

          ! calculate core-shell unit energy

          stat%engshl=stat%engshl+omega
          stat%virshl=stat%virshl+gamma*rabsq

          ! calculate stress tensor

          strs1 = strs1 + xdab(i)*fx
          strs2 = strs2 + xdab(i)*fy
          strs3 = strs3 + xdab(i)*fz
          strs5 = strs5 + ydab(i)*fy
          strs6 = strs6 + ydab(i)*fz
          strs9 = strs9 + zdab(i)*fz

        End If

        If (ib <= config%natms) Then

          config%parts(ib)%fxx=config%parts(ib)%fxx-fx
          config%parts(ib)%fyy=config%parts(ib)%fyy-fy
          config%parts(ib)%fzz=config%parts(ib)%fzz-fz

        End If

      End If
    End Do

    ! complete stress tensor

    stat%stress(1) = stat%stress(1) + strs1
    stat%stress(2) = stat%stress(2) + strs2
    stat%stress(3) = stat%stress(3) + strs3
    stat%stress(4) = stat%stress(4) + strs2
    stat%stress(5) = stat%stress(5) + strs5
    stat%stress(6) = stat%stress(6) + strs6
    stat%stress(7) = stat%stress(7) + strs3
    stat%stress(8) = stat%stress(8) + strs6
    stat%stress(9) = stat%stress(9) + strs9

    ! sum contributions to potential and virial

    If (comm%mxnode > 1) Then
      buffer(1)=stat%engshl
      buffer(2)=stat%virshl
      Call gsum(comm,buffer(1:2))
      stat%engshl=buffer(1)
      stat%virshl=buffer(2)
    End If

    Deallocate (lunsafe,lstopt, Stat=fail(1))
    Deallocate (xdab,ydab,zdab, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_forces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine core_shell_forces

  Subroutine core_shell_kinetic(config,shlke,cshell,domain,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the internal kinetic energy of
    ! core-shell units in the shell polarisation model
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2008
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent(   Out ) :: shlke
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Integer           :: i,j,k
    Real( Kind = wp ) :: rmu,rvx,rvy,rvz

    ! gather velocities of shared particles

    If (cshell%lshmv_shl) Then
      Call update_shared_units(config,cshell%lishp_shl,cshell%lashp_shl,config%vxx,&
                               config%vyy,config%vzz,domain,comm)
    End If

    ! initialise energy

    shlke=0.0_wp

    ! loop over all specified core-shell pairs

    Do k=1,cshell%ntshl

      ! indices of atoms involved

      i=local_index(cshell%listshl(1,k),config%nlast,config%lsi,config%lsa)
      j=local_index(cshell%listshl(2,k),config%nlast,config%lsi,config%lsa)

      ! for all native and natively shared core-shell units

      If ((i > 0 .and. j > 0) .and. (i <= config%natms .or. j <= config%natms)) Then

        ! calculate reduced mass

        rmu=(config%weight(i)*config%weight(j))/(config%weight(i)+config%weight(j))

        ! frozen particles' velocities are zero
        ! calculate constraint vector normal

        rvx=config%vxx(j)-config%vxx(i)
        rvy=config%vyy(j)-config%vyy(i)
        rvz=config%vzz(j)-config%vzz(i)

        ! calculate core-shell internal kinetic energy

        If (i <= config%natms .and. j <= config%natms) Then

          ! for native core-shell units - full contribution

          shlke=shlke+0.5_wp*rmu*(rvx*rvx+rvy*rvy+rvz*rvz)

        Else

          ! If ( (i <= config%natms .and. j > config%natms) .or. (j <= config%natms .and. i > config%natms) ) Then
          ! for shared core-shell units - halved contribution

          shlke=shlke+0.25_wp*rmu*(rvx*rvx+rvy*rvy+rvz*rvz)

        End If

      End If

    End Do

    ! global sum of core-shell internal kinetic energy

    Call gsum(comm,shlke)

  End Subroutine core_shell_kinetic

  Subroutine core_shell_on_top(cshell,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for positioning shells on top of their cores
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2010
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type( comms_type ),  Intent( InOut ) :: comm
    Type( core_shell_type ) , Intent( InOut ) :: cshell
    Type( configuration_type ), Intent( InOut ) :: config
    Logical :: safe
    Integer :: fail,i,j,ia,ib

    Logical, Allocatable :: lunsafe(:)
    Character( Len = 256 ) :: message

    fail=0
    Allocate (lunsafe(1:cshell%mxshl), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'core_shell_on_top allocation failure'
      Call error(0,message)
    End If


    ! Coincide shells with their cores

    Do i=1,cshell%ntshl
      lunsafe(i)=.false.

      ! indices of atoms in a core-shell

      ia=local_index(cshell%listshl(1,i),config%nlast,config%lsi,config%lsa) ! This is a core
      ib=local_index(cshell%listshl(2,i),config%nlast,config%lsi,config%lsa) ! This is a shell

      ! For every shell in the domain get the coordinates of its
      ! corresponding core (which must be in the domain+hello
      ! area by construction, if not go to a controlled termination)

      If (ib > 0 .and. ib <= config%natms .and. ia > 0) Then
        config%parts(ib)%xxx=config%parts(ia)%xxx
        config%parts(ib)%yyy=config%parts(ia)%yyy
        config%parts(ib)%zzz=config%parts(ia)%zzz
      End If

      ! Detect uncompressed unit

      If ( ((ia > 0 .and. ia <= config%natms) .or.   &
        (ib > 0 .and. ib <= config%natms)) .and. &
        (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
    End Do

    ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:cshell%ntshl))
    Call gcheck(comm,safe)
    If (.not.safe) Then
      Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
          Do i=1,cshell%ntshl
            If (lunsafe(i)) Then
              Write(message,'(a,2(i10,a))')     &
                'global unit number', cshell%listshl(0,i), &
                ' , with a head particle number', cshell%listshl(1,i),   &
                ' contributes towards next error'
              Call warning(message)
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(100)
    End If

    Deallocate (lunsafe, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'core_shell_on_top deallocation failure'
      Call error(0,message)
    End If

  End Subroutine core_shell_on_top

  Subroutine core_shell_quench(config,safe,temp,cshell,domain,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for quenching the internal bond energies of ions
    ! defined by shell model
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1999
    ! amended   - i.t.todorov november 2012
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,           Intent(   Out ) :: safe
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Real( Kind = wp ), Intent( In    ) :: temp
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type), Intent( InOut ) :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Logical           :: safek
    Integer           :: ia,ib,k
    Real( Kind = wp ) :: dvx,dvy,dvz,pke,rmu,scl,tke,tmx,tmy,tmz

    ! Initialise safe flag

    safe=.true.

    ! gather velocities of shared particles

    If (cshell%lshmv_shl) Then
      Call update_shared_units(config,cshell%lishp_shl,cshell%lashp_shl,config%vxx,config%vyy,config%vzz,domain,comm)
    End If

    ! permitted core-shell internal kinetic energy

    pke=boltz*temp*1.0e-4_wp

    ! amend core and shell velocities preserving total momentum

    Do k=1,cshell%ntshl
      ia=local_index(cshell%listshl(1,k),config%nlast,config%lsi,config%lsa)
      ib=local_index(cshell%listshl(2,k),config%nlast,config%lsi,config%lsa)

      If ((ia > 0 .and. ib > 0) .and. (ia <= config%natms .or. ib <= config%natms)) Then
        rmu=(config%weight(ia)*config%weight(ib))/(config%weight(ia)+config%weight(ib))

        ! frozen atoms have zero velocity (no worries)

        dvx=config%vxx(ib)-config%vxx(ia)
        dvy=config%vyy(ib)-config%vyy(ia)
        dvz=config%vzz(ib)-config%vzz(ia)

        tke=rmu*(dvx*dvx+dvy*dvy+dvz*dvz)

        safek = .not.( tke > pke .and. (tke-pke)/pke > 1.0e-3_wp )
        If (.not.safek) Then
          scl=Sqrt(pke/tke)

          tmx=config%weight(ia)*config%vxx(ia)+config%weight(ib)*config%vxx(ib)
          tmy=config%weight(ia)*config%vyy(ia)+config%weight(ib)*config%vyy(ib)
          tmz=config%weight(ia)*config%vzz(ia)+config%weight(ib)*config%vzz(ib)

          ! no corrections for frozen cores

          If (config%lfrzn(ia) == 0) Then
            config%vxx(ia)=tmx/(config%weight(ia)+config%weight(ib))-scl*rmu*dvx/config%weight(ia)
            config%vyy(ia)=tmy/(config%weight(ia)+config%weight(ib))-scl*rmu*dvy/config%weight(ia)
            config%vzz(ia)=tmz/(config%weight(ia)+config%weight(ib))-scl*rmu*dvz/config%weight(ia)

            config%vxx(ib)=tmx/(config%weight(ia)+config%weight(ib))+scl*rmu*dvx/config%weight(ib)
            config%vyy(ib)=tmy/(config%weight(ia)+config%weight(ib))+scl*rmu*dvy/config%weight(ib)
            config%vzz(ib)=tmz/(config%weight(ia)+config%weight(ib))+scl*rmu*dvz/config%weight(ib)
          Else
            config%vxx(ib)=tmx/(config%weight(ia)+config%weight(ib))+scl*rmu*dvx/config%weight(ib)
            config%vyy(ib)=tmy/(config%weight(ia)+config%weight(ib))+scl*rmu*dvy/config%weight(ib)
            config%vzz(ib)=tmz/(config%weight(ia)+config%weight(ib))+scl*rmu*dvz/config%weight(ib)
          End If

          safe=.false.
        End If
      End If
    End Do

    Call gcheck(comm,safe)

  End Subroutine core_shell_quench

  Subroutine core_shell_relax(l_str,relaxed,rdf_collect,rlx_tol,stpcfg,cshell, &
      stat,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for relaxing shells to zero force using conjugate
    ! gradient method
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & w.smith august 2014
    ! contrib   - a.m.elena february 2017
    ! contrib   - i.t.todorov february 2017
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( In    ) :: l_str
    Logical,            Intent( InOut ) :: relaxed,rdf_collect
    Real( Kind = wp ),  Intent( In    ) :: rlx_tol(1:2),stpcfg
    Type(core_shell_type), Intent( InOut ) :: cshell
    Type( stats_type ), Intent( InOut ) :: stat
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Integer                 :: fail(1:2),i,ia,ib,jshl

    ! OUTPUT existence

    Logical               :: l_out
    Character( Len = 10 ) :: c_out

    ! Optimisation iteration and convergence limits

    Integer,      Parameter :: mxpass = 100


    Integer,           Allocatable       :: lstopt(:,:),lst_sh(:)
    Real( Kind = wp ), Allocatable       :: fxt(:),fyt(:),fzt(:)
    Character( Len = 256 ) :: message

    fail=0
    Allocate (lstopt(1:2,1:cshell%mxshl),lst_sh(1:mxatms),      Stat=fail(1))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_relax allocation failure'
      Call error(0,message)
    End If

    If (cshell%newjob) Then
      cshell%newjob = .false.

      ! At start no optimisation has been attempted yet

      cshell%keyopt = NO_SEARCH

      ! Passage accumulators are initialised in core_shell
      ! stat%passshl(1) - cycles counter
      ! stat%passshl(2) - access counter
      ! stat%passshl(3) - average cycles
      ! stat%passshl(4) - minimum cycles
      ! stat%passshl(5) - maximum cycles

      Allocate (cshell%oxt(1:cshell%mxshl),cshell%oyt(1:cshell%mxshl),cshell%ozt(1:cshell%mxshl), Stat=fail(1))
      If (fail(1) > 0) Then
        Write(message,'(a)') 'core_shell_relax allocation failure SAVE'
        Call error(0,message)
      End If
    End If

    ! Step length for relaxation

    If (rlx_tol(2) > zero_plus) Then

      ! Optionally specified

      cshell%step=rlx_tol(2)

    Else

      ! default if unspecified

      cshell%step=0.5_wp/cshell%smax

    End If

    If (cshell%keyopt == NO_SEARCH) Then

      ! No relaxation is yet attempted

      relaxed=.false.

      ! Minimum needed for a pass for this minimisation cycle

      cshell%grad_pass = Huge(1.0_wp)

      ! Avoid rdf calculation redundancy

      cshell%l_rdf=rdf_collect
      If (rdf_collect) rdf_collect=.false.

      ! Print header

      If (l_str) Then
        Write(message,'(a,3x,a,6x,a,11x,a,8x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
          'Relaxing shells to cores:','pass','eng_tot','grad_tol','dis_tol','dcs_max','tol=',rlx_tol(1),'step=',cshell%step
        Call info(message,.true.)
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)
      End If
    End If

    ! gather new forces on shared shells

    If (cshell%lshmv_shl) Then
      Call update_shared_units(config,cshell%lishp_shl,cshell%lashp_shl,&
                               SHARED_UNIT_UPDATE_FORCES,domain,comm)
    End If

    ! Load shell forces on cores (cores don't move during the shell relaxation)

    fxt(1:cshell%mxshl)=0.0_wp ; fyt(1:cshell%mxshl)=0.0_wp ; fzt(1:cshell%mxshl)=0.0_wp
    jshl=0
    Do i=1,cshell%ntshl
      ia=local_index(cshell%listshl(1,i),config%nlast,config%lsi,config%lsa)
      ib=local_index(cshell%listshl(2,i),config%nlast,config%lsi,config%lsa)
      If (ia > 0 .and. ia <= config%natms) Then ! THERE IS AN ib>0 FOR SURE
        jshl=jshl+1
        fxt(jshl)=config%parts(ib)%fxx
        fyt(jshl)=config%parts(ib)%fyy
        fzt(jshl)=config%parts(ib)%fzz
     End If
     lstopt(1,i)=ia
     lstopt(2,i)=ib
  End Do

    ! Current configuration energy

    cshell%eng=stpcfg

    ! Initialise/get eng_tol & verify relaxed condition

    cshell%eng_tol=0.0_wp
    If (cshell%keyopt > NO_SEARCH) Then
      cshell%eng_tol=Abs(1.0_wp-cshell%eng2/cshell%eng)
    End If
    ! Current gradient (modulus of the total force on shells)

    cshell%grad=0.0_wp
    Do i=1,jshl
      cshell%grad=cshell%grad+fxt(i)**2+fyt(i)**2+fzt(i)**2
    End Do
    Call gsum(comm,cshell%grad)
    cshell%grad=Sqrt(cshell%grad)

    ! Get grad_tol & verify relaxed condition

    cshell%grad_tol=cshell%grad/Real(cshell%megshl,wp)
    relaxed=(cshell%grad_tol < rlx_tol(1))

    ! Initialise cshell%dist_tol

    cshell%dist_tol=0.0_wp

    ! CHECK FOR CONVERGENCE

    If (.not.relaxed) Then

      ! Increment main passage counter

      stat%passshl(1)=stat%passshl(1)+1.0_wp

      ! Minimum for passing

      cshell%grad_pass = Min(cshell%grad_pass,cshell%grad_tol)

      ! If in mxpass iterations we are not there, give up but
      ! allow for ten-fold boost in iteration cycle length
      ! for the very first MD step

      If (Nint(stat%passshl(2)) == 0) Then
        If (Nint(stat%passshl(1)) >= 10*mxpass) Then
          Call warning(330,rlx_tol(1),cshell%grad_pass,0.0_wp)
          Call error(474)
        End If
      Else
        If (Nint(stat%passshl(1)) >= mxpass) Then
          Call warning(330,rlx_tol(1),cshell%grad_pass,0.0_wp)
          Call error(474)
        End If
      End If

    Else

      Go To 100

    End If

    If      (cshell%keyopt == NO_SEARCH) Then

      ! Original configuration energy

      cshell%eng0=cshell%eng
      cshell%eng2=cshell%eng

      ! Original gradient (modulus of the total force on shells)

      cshell%onorm=cshell%grad
      cshell%grad0=cshell%grad
      cshell%grad2=cshell%grad

      ! Set original search direction

      cshell%oxt=0.0_wp ; cshell%oyt=0.0_wp ; cshell%ozt=0.0_wp
      Do i=1,jshl
        cshell%oxt(i)=fxt(i)
        cshell%oyt(i)=fyt(i)
        cshell%ozt(i)=fzt(i)
      End Do

      cshell%keyopt=LINE_SEARCH
      cshell%sgn=1.0_wp
      cshell%stride=cshell%sgn*cshell%step

    Else If (cshell%keyopt == LINE_SEARCH) Then

      ! Line search along chosen direction

      cshell%eng1=cshell%eng0
      cshell%eng2=cshell%eng

      cshell%grad1=cshell%grad2
      cshell%grad2=0.0_wp
      Do i=1,jshl
        cshell%grad2=cshell%grad2+cshell%oxt(i)*fxt(i)+cshell%oyt(i)*fyt(i)+cshell%ozt(i)*fzt(i)
      End Do
      Call gsum(comm,cshell%grad2)
      cshell%grad2=cshell%sgn*cshell%grad2/cshell%onorm

      ! Linear extrapolation to minimum

      If (cshell%grad2 < 0.0_wp) Then ! BACK UP FROM THIS DIRECTION
        cshell%keyopt=CONJUGATE_SEARCH
        cshell%stride=cshell%sgn*cshell%step*cshell%grad2/(cshell%grad1-cshell%grad2)
      Else                     ! CARRY ON IN THIS DIRECTION
        cshell%stride=cshell%sgn*cshell%step
      End If

    Else If (cshell%keyopt == CONJUGATE_SEARCH) Then

      ! Construct conjugate search vector

      cshell%eng1=cshell%eng2
      cshell%eng2=cshell%eng

      cshell%gamma=(cshell%grad/cshell%grad0)**2
      cshell%grad0=cshell%grad
      cshell%grad2=0.0_wp
      cshell%onorm=0.0_wp
      Do i=1,jshl
        cshell%oxt(i)=fxt(i)+cshell%gamma*cshell%oxt(i)
        cshell%oyt(i)=fyt(i)+cshell%gamma*cshell%oyt(i)
        cshell%ozt(i)=fzt(i)+cshell%gamma*cshell%ozt(i)

        cshell%onorm=cshell%onorm+cshell%oxt(i)**2+cshell%oyt(i)**2+cshell%ozt(i)**2
        cshell%grad2=cshell%grad2+cshell%oxt(i)*fxt(i)+cshell%oyt(i)*fyt(i)+cshell%ozt(i)*fzt(i)
      End Do
      Call gsum(comm,cshell%onorm)
      cshell%onorm=Sqrt(cshell%onorm)
      Call gsum(comm,cshell%grad2)
      cshell%grad2=cshell%grad2/cshell%onorm
      cshell%sgn=Sign(1.0_wp,cshell%grad2)
      cshell%grad2=cshell%sgn*cshell%grad2

      cshell%keyopt=LINE_SEARCH
      cshell%stride=cshell%sgn*cshell%step

    End If

    ! Load original shell forces on their cores in DD representation

    fxt(1:mxatdm)=0.0_wp ; fyt(1:mxatdm)=0.0_wp ; fzt(1:mxatdm)=0.0_wp
    jshl=0
    Do i=1,cshell%ntshl
      ia=lstopt(1,i)
      If (ia > 0 .and. ia <= config%natms) Then
        jshl=jshl+1
        fxt(ia)=cshell%oxt(jshl)
        fyt(ia)=cshell%oyt(jshl)
        fzt(ia)=cshell%ozt(jshl)
      End If
    End Do

    ! Exchange original shell forces on shared cores across domains

    If (cshell%lshmv_shl) Then
      Call update_shared_units(config,cshell%lishp_shl,cshell%lashp_shl,fxt,fyt,fzt,domain,comm)
    End If

    ! Move shells accordingly to their new positions

   Do i=1,cshell%ntshl
     ia=lstopt(1,i)
     ib=lstopt(2,i)
     If (ia > 0 .and. (ib > 0 .and. ib <= config%natms)) Then
        config%parts(ib)%xxx=config%parts(ib)%xxx+cshell%stride*fxt(ia)
        config%parts(ib)%yyy=config%parts(ib)%yyy+cshell%stride*fyt(ia)
        config%parts(ib)%zzz=config%parts(ib)%zzz+cshell%stride*fzt(ia)
        cshell%dist_tol(1)=Max(cshell%dist_tol(1),fxt(ia)**2+fyt(ia)**2+fzt(ia)**2) ! - shell move
        cshell%x(1)=config%parts(ib)%xxx-config%parts(ia)%xxx 
        cshell%y(1)=config%parts(ib)%yyy-config%parts(ia)%yyy
        cshell%z(1)=config%parts(ib)%zzz-config%parts(ia)%zzz
        Call images(config%imcon,config%cell,1,cshell%x,cshell%y,cshell%z)
        cshell%dist_tol(2)=Max(cshell%dist_tol(2),cshell%x(1)**2+cshell%y(1)**2+cshell%z(1)**2) ! - core-shell separation
      End If
    End Do
    cshell%dist_tol=Sqrt(cshell%dist_tol)
    cshell%dist_tol(1)=cshell%dist_tol(1)*Abs(cshell%stride)
    Call gmax(comm,cshell%dist_tol)

    ! Fit headers in and Close and Open OUTPUT at every 25th print-out

    i=Nint(stat%passshl(1))
    If (l_str) Then
      Write(message,'(1x,i31,1x,1p,2e18.8,4x,f7.4,4x,f7.4,12x,e18.8)') i-1,stpcfg/engunit,cshell%grad_tol,&
        cshell%dist_tol(1),cshell%dist_tol(2),cshell%eng_tol
      Call info(message,.true.)
      If (Mod(i,25) == 0) Then
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)
        Write(message,'(1x,a,3x,a,6x,a,11x,a,9x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
          'Relaxing shells to cores:','pass','eng_tot','grad_tol','ds_tol','dcs_max','tol=',rlx_tol(1),'step=',cshell%step
        Call info(message,.true.)
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)

        If (comm%idnode == 0) Then
          Inquire(File=Trim(output), Exist=l_out, Position=c_out)
          Call strip_blanks(c_out)
          Call lower_case(c_out)
          If (l_out .and. c_out(1:6) == 'append') Then
            Close(Unit=nrite)
            Open(Unit=nrite, File=Trim(output), Position='append')
          End If
        End If
      End If
    End If

    100 Continue

    If (relaxed) Then

      ! Final printout

      i=Nint(stat%passshl(1))
      If (.not.l_str) Then
        Write(message,'(a,4x,a,6x,a,11x,a,8x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
          'Relaxed shells to cores:','pass','eng_tot','grad_tol','dis_tol','dcs_max','tol=',rlx_tol(1),'step=',cshell%step
        Call info(message,.true.)
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)
      End If
      Write(message,'(1x,i31,1x,1p,2e18.8,4x,f7.4,4x,f7.4,12x,e18.8)') &
        i-1,stpcfg/engunit,cshell%grad_tol,cshell%dist_tol(1),cshell%dist_tol(2),cshell%eng_tol
      Call info(message,.true.)
      Write(message,"(1x,130('-'))")
      Call info(message,.true.)

      ! Collect passage statistics

      stat%passshl(3)=stat%passshl(2)*stat%passshl(3)
      stat%passshl(2)=stat%passshl(2)+1.0_wp
      stat%passshl(3)=stat%passshl(3)/stat%passshl(2)+stat%passshl(1)/stat%passshl(2)
      stat%passshl(4)=Min(stat%passshl(1),stat%passshl(4))
      stat%passshl(5)=Max(stat%passshl(1),stat%passshl(5))

      ! Rewind keyopt and main passage counter

      cshell%keyopt = NO_SEARCH
      stat%passshl(1)=0.0_wp

      ! Resume rdf calculations

      If (cshell%l_rdf) rdf_collect=cshell%l_rdf

      ! Zero shells' velocities and forces and redistribute
      ! the residual force to the rest of the system to prevent
      ! COM force generation

      lst_sh(1:config%natms)=0
      cshell%fff(0)=Real(config%natms,wp)
      cshell%fff(1:3)=0.0_wp
      Do i=1,cshell%ntshl
        ib=lstopt(2,i)
        If (ib > 0) Then
           lst_sh(ib)=1
           cshell%fff(0)=cshell%fff(0)-1.0_wp
           cshell%fff(1)=cshell%fff(1)+config%parts(ib)%fxx ; config%parts(ib)%fxx=0.0_wp 
           config%vxx(ib)=0.0_wp
           cshell%fff(2)=cshell%fff(2)+config%parts(ib)%fyy ; config%parts(ib)%fyy=0.0_wp 
           config%vyy(ib)=0.0_wp
           cshell%fff(3)=cshell%fff(3)+config%parts(ib)%fzz ; config%parts(ib)%fzz=0.0_wp 
           config%vzz(ib)=0.0_wp
        End If
      End Do
      Call gsum(comm,cshell%fff)
      cshell%fff(1:3)=cshell%fff(1:3)/cshell%fff(0)
      Do i=1,config%natms
        If (lst_sh(i) == 0) Then
           config%parts(i)%fxx=config%parts(i)%fxx+cshell%fff(1)
           config%parts(i)%fyy=config%parts(i)%fyy+cshell%fff(2)
           config%parts(i)%fzz=config%parts(i)%fzz+cshell%fff(3)
        End If
      End Do

      ! Frozen atoms option

      Call freeze_atoms(config)
    End If

    Deallocate (lstopt,lst_sh, Stat=fail(1))
    Deallocate (fxt,fyt,fzt,   Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_relax deallocation failure'
      Call error(0,message)
    End If

  End Subroutine core_shell_relax
End module core_shell
