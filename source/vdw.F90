Module vdw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global VdW interaction variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use comms,  Only : comms_type,gsum,gbcast
  Use setup
  Use site, Only : site_type
  Use configuration, Only : imcon,volm,natms,ltype,lfrzn, &
                            ltg,fxx,fyy,fzz
  Use mm3lrc
  Use zbl_pots,         Only : ab, intRadZBL, intdRadZBL, &
                           zbl,zbls,zblb

  Use site, Only : site_type
  Use parse, Only : get_line,get_word,word_2_real
  Use neighbours, Only : neighbours_type
  Use errors_warnings, Only : error,warning,info
  Implicit None

  Private

  ! Mixing rule parameters
  !> Null
  Integer( Kind = wi), Parameter, Public :: MIX_NULL = 0
  !> Lorentz-Berthelot: $e_{ij}=(e_i*e_j)^{1/2} \quad s_{ij}=(s_i+s_j)/2$
  Integer( Kind = wi ), Parameter, Public :: MIX_LORENTZ_BERTHELOT = 1
  !> Fender-Hasley: $e_{ij}=(2*e_i*e_j)/(e_i+e_j) \quad s_{ij}=(s_i+s_j)/2$
  Integer( Kind = wi ), Parameter, Public :: MIX_FENDER_HASLEY = 2
  !> Hogervorst Good-Hope: $e_{ij}=(e_i*e_j)^{1/2} \quad s_{ij}=(s_i*s_j)^{1/2}$
  Integer( Kind = wi ), Parameter, Public :: MIX_HOGERVORST = 3
  !> Halgren HHG: $e_{ij}=(4*e_i*e_j)/(e_i^{1/2}+e_j^{1/2})^2 \quad s_{ij}=(s_i^3+s_j^3)/(s_i^2+s_j^2)$
  Integer( Kind = wi ), Parameter, Public :: MIX_HALGREN = 4
  !> Waldman–Hagler: $e_{ij}=2*(e_i*e_j)^{1/2}*(s_i*s_j)^3/(s_i^6+s_j^6) \quad s_{ij}=[(s_i^6+s_j^6)/2]^{1/6}$
  Integer( Kind = wi ), Parameter, Public :: MIX_WALDMAN_HAGLER = 5
  !> Tang-Toennies: $e_{ij}=[(e_i*s_i^6)*(e_j*s_j^6)] / ([(e_i*s_i^12)^{1/13}+(e_j*s_j^12)^{1/13}]/2)^13$
  !>                $s_{ij}=(1/3) \sum_{L=0}^2 [(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^{1/(6-2L)}$
  Integer( Kind = wi ), Parameter, Public :: MIX_TANG_TOENNIES = 6
  !> Functional: $e_{ij}=3*(e_i*e_j)^{1/2} * (s_i*s_j)^3 / \sum_{L=0}^2 [(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(6/(6-2L))$
  !>             $s_ij=(1/3) \sum_{L=0}^2 [(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))$
  Integer( Kind = wi ), Parameter, Public :: MIX_FUNCTIONAL = 7

  !> Type containing Van der Waals data
  Type, Public :: vdw_type
    Private

    !> Flag for any tabulated potential
    Logical, Public :: l_tab = .false.
    !> Direct calculation flag
    Logical, Public :: l_direct = .false.
    !> Force shifting flag
    Logical, Public :: l_force_shift = .false.

    !> Number of two body interactoins
    Integer( Kind = wi ), Public :: n_vdw = 0
    !> Mixing type
    Integer( Kind = wi ), Public :: mixing = MIX_NULL

    Integer( Kind = wi ), Allocatable, Public :: list(:)
    Integer( Kind = wi ), Allocatable, Public :: ltp(:)

    !> VdW parameters
    Real( Kind = wp ), Allocatable, Public :: param(:,:)
    !> VdW cut off
    Real( Kind = wp ), Public :: cutoff

    Real( Kind = wp ), Allocatable, Public :: sigeps(:,:)

    !> Energy long range correction
    Real( Kind = wp ), Public :: elrc
    !> Virial long range correction
    Real( Kind = wp ), Public :: vlrc

    ! Possible tabulated calculation arrays
    !> Tabulated potential
    Real( Kind = wp ), Allocatable, Public :: tab_potential(:,:)
    !> Tabulated force
    Real( Kind = wp ), Allocatable, Public :: tab_force(:,:)
    !> Maximum number of grid points
    Integer( Kind = wi ), Public :: max_grid

    ! Possible force-shifting arrays
    Real( Kind = wp ), Allocatable, Public :: afs(:)
    Real( Kind = wp ), Allocatable, Public :: bfs(:)

    !> Maximum number of VdW interations
    Integer( Kind = wi ), Public :: max_vdw
    !> Maximum number of VdW parameters
    Integer( Kind = wi ), Public :: max_param
  Contains
    Private

    Procedure, Public :: init => allocate_vdw_arrays
    Procedure, Public :: init_table => allocate_vdw_table_arrays
    Procedure, Public :: init_direct => allocate_vdw_direct_fs_arrays
    Final :: cleanup
  End Type vdw_type

  Public :: vdw_forces,vdw_generate,vdw_table_read,vdw_lrc

Contains

  Subroutine allocate_vdw_arrays(T)
    Class( vdw_type ) :: T

    Integer, Dimension( 1:4 ) :: fail

    fail = 0

    Allocate (T%list(1:T%max_vdw),          Stat = fail(1))
    Allocate (T%ltp(1:T%max_vdw),          Stat = fail(2))
    Allocate (T%param(1:T%max_param,1:T%max_vdw), Stat = fail(3))
    Allocate (T%sigeps(1:2,1:T%max_vdw),      Stat = fail(4))

    If (Any(fail > 0)) Call error(1022)

    T%list = 0
    T%ltp = 0

    T%param = 0.0_wp
    T%sigeps = 0.0_wp
  End Subroutine allocate_vdw_arrays

  Subroutine allocate_vdw_table_arrays(T)
    Class( vdw_type ) :: T

    Integer, Dimension( 1:2 ) :: fail

    fail = 0

    Allocate (T%tab_potential(0:T%max_grid,1:T%max_vdw),   Stat = fail(1))
    Allocate (T%tab_force(0:T%max_grid,1:T%max_vdw),   Stat = fail(2))

    If (Any(fail > 0)) Call error(1063)

    T%tab_potential = 0.0_wp
    T%tab_force = 0.0_wp
  End Subroutine allocate_vdw_table_arrays

  Subroutine allocate_vdw_direct_fs_arrays(T)
    Class( vdw_type ) :: T

    Integer :: fail

    fail = 0

    Allocate (T%afs(1:T%max_vdw),T%bfs(1:T%max_vdw), Stat = fail)

    If (fail > 0) Call error(1066)

    T%afs  = 0.0_wp
    T%bfs  = 0.0_wp
  End Subroutine allocate_vdw_direct_fs_arrays

  Subroutine vdw_lrc(site,vdw,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to evaluate vdw long-range corrections to
! pressure and energy in a 3D periodic system
!
! copyright - daresbury laboratory
! author    - t.forester may 1993
! amended   - i.t.todorov september 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (rydberg)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( site_type ), Intent( In    ) :: site
  Type( vdw_type ), Intent( InOut ) :: vdw
  Type( comms_type ), Intent( InOut ) :: comm

  Integer           :: fail,i,j,k,ivdw,keypot,n,m
  Real( Kind = wp ) :: a,b,c,d,e0,nr,mr,r0,r,eps,sig, &
                       eadd,padd,denprd,plrc,t,kk,s9, &
                       z1,z2,rm,al

  Real( Kind = wp ), Dimension( : ), Allocatable :: numfrz
  Character( Len = 256 ) :: message,messages(3)

  fail=0
  Allocate (numfrz(mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'vdw_lrc allocation failure'
     Call error(0,message)
  End If

! initialise long-range corrections to energy and pressure

  plrc = 0.0_wp
  vdw%elrc = 0.0_wp

  If (vdw%l_force_shift) Go To 10 ! force-shifting

! initialise counter arrays and evaluate number density in system

  numfrz = 0.0_wp
  Do i=1,natms
     k = ltype(i)
     If (lfrzn(i) /= 0) numfrz(k)=numfrz(k)+1.0_wp
  End Do
  Call gsum(comm,numfrz(1:site%ntype_atom))

! Evaluate only for 3D periodic systems

  If (imcon /= 0 .and. imcon /= 6) Then
     ivdw = 0

     Do i=1,site%ntype_atom
        Do j=1,i

           eadd = 0.0_wp
           padd = 0.0_wp

           ivdw = ivdw + 1
           k = vdw%list(ivdw)

           keypot=vdw%ltp(k)
           If      (keypot ==  0) Then

! tabulated energy and pressure lrc

              eadd = vdw%param(1,k)
              padd =-vdw%param(2,k)

           Else If (keypot ==  1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=vdw%param(1,k)
              b=vdw%param(2,k)
              r=vdw%cutoff

              eadd = a/(9.0_wp*r**9) - b/(3.0_wp*r**3)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 6.0_wp*b/(3.0_wp*r**3)

           Else If (keypot ==  2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)
              r  =vdw%cutoff

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - sig**6/(3.0_wp*r**3))
              padd = 8.0_wp*eps*(6.0_wp*sig**12/(9.0_wp*r**9) - sig**6/(r**3))

           Else If (keypot ==  3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

              e0=vdw%param(1,k)
              n =Nint(vdw%param(2,k)) ; nr=Real(n,wp)
              m =Nint(vdw%param(3,k)) ; mr=Real(m,wp)
              r0=vdw%param(4,k)
              r =vdw%cutoff

              eadd = e0/(nr-mr)*( mr*r0**n/((nr-3.0_wp)*r**(n-3)) - nr*r0**m/((mr-3.0_wp)*r**(m-3)) )
              padd = e0/(nr-mr)*nr*mr*( r0**n/((nr-3.0_wp)*r**(n-3)) - r0**m/((mr-3.0_wp)*r**(m-3)) )

           Else If (keypot ==  4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              c=vdw%param(3,k)
              r=vdw%cutoff

              eadd = -c/(3.0_wp*r**3)
              padd = -2.0_wp*c/(r**3)

           Else If (keypot ==  5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              c=vdw%param(4,k)
              d=vdw%param(5,k)
              r=vdw%cutoff

              eadd = -c/(3.0_wp*r**3) - d/(5.0_wp*r**5)
              padd = -2.0_wp*c/(r**3) - 8.0_wp*d/(5.0_wp*r**5)

           Else If (keypot ==  6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=vdw%param(1,k)
              b=vdw%param(2,k)
              r=vdw%cutoff

              eadd = a/(9.0_wp*r**9) - b/(7.0_wp*r**7)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 10.0_wp*b/(7.0_wp*r**7)

           Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

              e0=vdw%param(1,k)
              r0=vdw%param(2,k)
              kk=vdw%param(3,k)
              If (kk > Tiny(kk)) Then
                 t = Exp(-kk*(vdw%cutoff - r0))

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*vdw%cutoff+1)**2 + 1) + &
                    e0*t*t/(4.0_wp*kk*kk*kk)*((kk*vdw%cutoff+1)**2 + kk*kk*vdw%cutoff*vdw%cutoff)
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*vdw%cutoff**3 + &
                      3*kk**2*vdw%cutoff**2 +6*kk*vdw%cutoff + 6) + &
                      e0*t*t/(4.0_wp*kk*kk*kk)* & 
                      (4.0_wp*kk**3*vdw%cutoff**3 + 6*kk**2*vdw%cutoff**2 + 6*kk*vdw%cutoff + 3)
              End If

           Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)

              a =0.07_wp
              b =0.12_wp
              e0=1.0e-12_wp

              eadd = intRadMM3(sig,a,b,eps,vdw%cutoff,e0)
              padd = -intRaddMM3(sig,a,b,eps,vdw%cutoff,e0)

           Else If (keypot ==  12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)
              c  =vdw%param(3,k)
              r  =vdw%cutoff

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - c*sig**6/(3.0_wp*r**3))
              padd = 8.0_wp*eps*(6.0_wp*sig**12/(9.0_wp*r**9) - c*sig**6/(r**3))

           Else If (keypot == 13) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

              e0 = vdw%param(1,k)
              r0 = vdw%param(2,k)
              kk = vdw%param(3,k)
               c = vdw%param(4,k)

              If (kk>Tiny(kk)) Then

                 t = Exp(-kk*(vdw%cutoff - r0))
                 s9 = c/(9.0_wp*vdw%cutoff**9)

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*vdw%cutoff+1)**2 + 1) + &
                     e0*t*t/(4.0_wp*kk*kk*kk)*((kk*vdw%cutoff+1)**2 + & 
                     kk*kk*vdw%cutoff*vdw%cutoff) + s9
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*vdw%cutoff**3 + & 
                       3*kk**2*vdw%cutoff**2 + 6*kk*vdw%cutoff + 6) + & 
                       e0*t*t/(4.0_wp*kk*kk*kk)* (4.0_wp*kk**3*vdw%cutoff**3 + & 
                       6*kk**2*vdw%cutoff**2 + 6*kk*vdw%cutoff + 3) + 12.0_wp*s9
              End If

           Else If (keypot == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)

              a = vdw%param(1,k)
              b = vdw%param(2,k)
              c = vdw%param(3,k)
              t = exp(-vdw%cutoff/c)

              eadd = (b*c*vdw%cutoff**3+(3*b*c**2+a*c)*vdw%cutoff**2+(6*b*c**3+2*a*c**2)*vdw%cutoff&
                +6*b*c**4+2*a*c**3)*t
              padd = (b*vdw%cutoff**4+(3*b*c+a)*vdw%cutoff**3+(9*b*c**2+3*a*c)*vdw%cutoff**2+& 
                (18*b*c**3+6*a*c**2)*vdw%cutoff+18*b*c**4+6*a*c**3)*t

           Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

              z1 = vdw%param(1,k)
              z2 = vdw%param(2,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0
              eadd = intRadZBL(kk,a,vdw%cutoff,1e-12_wp)
              padd = intdRadZBL(kk,a,vdw%cutoff,1e-12_wp)

           Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

              e0 = vdw%param(5,k)
              r0 = vdw%param(6,k)
              kk = vdw%param(7,k)

              If (kk > Tiny(kk)) Then
                 t = Exp(-kk*(vdw%cutoff - r0))

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*vdw%cutoff+1)**2 + 1) + &
                    e0*t*t/(4.0_wp*kk*kk*kk)*((kk*vdw%cutoff+1)**2 + kk*kk*vdw%cutoff*vdw%cutoff)
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*vdw%cutoff**3 + &
                      3*kk**2*vdw%cutoff**2 +6*kk*vdw%cutoff + 6) + &
                      e0*t*t/(4.0_wp*kk*kk*kk)* & 
                      (4.0_wp*kk**3*vdw%cutoff**3 + 6*kk**2*vdw%cutoff**2 + 6*kk*vdw%cutoff + 3)
              End If

           Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

              A = vdw%param(5,k)
              r0 = vdw%param(6,k)
              c = vdw%param(7,k)

              t=A*Exp(-vdw%cutoff/r0)

              eadd = (vdw%cutoff**2+2*r0*vdw%cutoff+2*r0**2)*t*r0-c/(3.0_wp*vdw%cutoff**3)
              padd = (vdw%cutoff**3+3*r0*vdw%cutoff**2+6*r0**2*vdw%cutoff+6*r0**3)*t -2.0_wp*c/(vdw%cutoff**3)

           End If

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

           If (i /= j) Then
              eadd = eadd*2.0_wp
              padd = padd*2.0_wp
           End If

           denprd=twopi * (site%num_type(i)*site%num_type(j) - numfrz(i)*numfrz(j)) / volm**2

           vdw%elrc = vdw%elrc + volm*denprd*eadd
           plrc = plrc + denprd*padd/3.0_wp

        End Do
     End Do

  End If

10 Continue

  Write(messages(1),'(a)') 'long-range correction for:'
  Write(messages(2),'(2x,a,e15.6)') 'vdw energy ',vdw%elrc/engunit
  Write(messages(3),'(2x,a,e15.6)') 'vdw pressure ',plrc*prsunt
  Call info(messages,3,.true.)

! convert plrc to a viral term

  vdw%vlrc = plrc*(-3.0_wp*volm)

  Deallocate (numfrz, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'vdw_lrc deallocation failure'
     Call error(0,message)
  End If

End Subroutine vdw_lrc

Subroutine vdw_direct_fs_generate(vdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating force-shifted constant arrays for
! direct vdw evaluation
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (ryd)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Type( vdw_type ), Intent( InOut ) :: vdw

  Integer           :: ivdw,keypot,n,m
  Real( Kind = wp ) :: r0,r0rn,r0rm,r_6,sor6,   &
                       rho,a,b,c,d,e0,kk,nr,mr, &
                       sig,eps,t1,t2,t3,z,dz,   &
                       z1,z2,rm,ic,k

! allocate arrays for force-shifted corrections

  Call vdw%init_direct()

! construct arrays for all types of vdw potential

  Do ivdw=1,vdw%n_vdw

     keypot=vdw%ltp(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=vdw%param(1,ivdw)
        b=vdw%param(2,ivdw)

        r_6=vdw%cutoff**(-6)

        vdw%afs(ivdw) = 6.0_wp*r_6*(2.0_wp*a*r_6-b)
        vdw%bfs(ivdw) =-r_6*(a*r_6-b) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)

        sor6=(sig/vdw%cutoff)**6

        vdw%afs(ivdw) = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        vdw%bfs(ivdw) =-4.0_wp*eps*sor6*(sor6-1.0_wp) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=vdw%param(1,ivdw)
        n =Nint(vdw%param(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(vdw%param(3,ivdw)) ; mr=Real(m,wp)
        r0=vdw%param(4,ivdw)

        a=r0/vdw%cutoff
        b=1.0_wp/Real(n-m,wp)
        r0rn=a**n
        r0rm=a**m

        vdw%afs(ivdw) = e0*mr*nr*(r0rn-r0rm)*b
        vdw%bfs(ivdw) =-e0*(mr*r0rn-nr*r0rm)*b - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a  =vdw%param(1,ivdw)
        rho=vdw%param(2,ivdw)
        c  =vdw%param(3,ivdw)

        If (Abs(rho) <= zero_plus) Then
           If (Abs(a) <= zero_plus) Then
              rho=1.0_wp
           Else
              Call error(467)
           End If
        End If

        b=vdw%cutoff/rho
        t1=a*Exp(-b)
        t2=-c/vdw%cutoff**6

        vdw%afs(ivdw) = (t1*b+6.0_wp*t2)
        vdw%bfs(ivdw) =-(t1+t2) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =vdw%param(1,ivdw)
        b  =vdw%param(2,ivdw)
        sig=vdw%param(3,ivdw)
        c  =vdw%param(4,ivdw)
        d  =vdw%param(5,ivdw)

        t1=a*Exp(b*(sig-vdw%cutoff))
        t2=-c/vdw%cutoff**6
        t3=-d/vdw%cutoff**8

        vdw%afs(ivdw) = (t1*vdw%cutoff*b+6.0_wp*t2+8.0_wp*t3)
        vdw%bfs(ivdw) =-(t1+t2+t3) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=vdw%param(1,ivdw)
        b=vdw%param(2,ivdw)

        t1=a/vdw%cutoff**12
        t2=-b/vdw%cutoff**10

        vdw%afs(ivdw) = (12.0_wp*t1+10.0_wp*t2)
        vdw%bfs(ivdw) =-(t1+t2) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=vdw%param(1,ivdw)
        r0=vdw%param(2,ivdw)
        kk=vdw%param(3,ivdw)

        t1=Exp(-kk*(vdw%cutoff-r0))

        vdw%afs(ivdw) =-2.0_wp*e0*kk*t1*(1.0_wp-t1)*vdw%cutoff
        vdw%bfs(ivdw) =-e0*t1*(t1-2.0_wp) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)
        d  =vdw%param(3,ivdw)

        sor6=(sig/(vdw%cutoff-d))**6

        vdw%afs(ivdw) = (24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(vdw%cutoff-d))*vdw%cutoff
        vdw%bfs(ivdw) =-(4.0_wp*eps*sor6*(sor6-1.0_wp)+eps) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 10) Then ! all zeroed in vdw

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

!       vdw%afs(ivdw) = 0.0_wp !initialised in vdw
!       vdw%bfs(ivdw) = 0.0_wp !initialised in vdw

     Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)

        rho=sig/vdw%cutoff
        t1=1.0_wp/(0.07_wp+rho)
        t2=1.0_wp/(0.12_wp+rho**7)
        t3=eps*(1.07_wp/t1**7)

        vdw%afs(ivdw) =-7.0_wp*t3*rho*(((1.12_wp/t2)-2.0_wp)/t1 + (1.12_wp/t2**2)*rho**6)
        vdw%bfs(ivdw) =-t3*((1.12_wp/t2)-2.0_wp) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

      Else If (keypot == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)
        c  =vdw%param(3,ivdw)

        sor6=(sig/vdw%cutoff)**6

        vdw%afs(ivdw) = 24.0_wp*eps*sor6*(2.0_wp*sor6-c)
        vdw%bfs(ivdw) =-4.0_wp*eps*sor6*(sor6-c) - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 13) Then

! Morse potential with twelve term:: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0=vdw%param(1,ivdw)
        r0=vdw%param(2,ivdw)
        kk=vdw%param(3,ivdw)
        c=vdw%param(4,ivdw)

        t1=Exp(-kk*(vdw%cutoff-r0))
        sor6 = c/vdw%cutoff**12

        vdw%afs(ivdw) =-2.0_wp*e0*kk*t1*(1.0_wp-t1)*vdw%cutoff + 12.0_wp*sor6
        vdw%bfs(ivdw) =-e0*t1*(t1-2.0_wp) + sor6 - vdw%afs(ivdw)
        vdw%afs(ivdw) = vdw%afs(ivdw)/vdw%cutoff

     Else If (keypot == 14) Then

! Morse potential with twelve term:: u=(a+b*r)Exp(-r/c)

        a = vdw%param(1,ivdw)
        b = vdw%param(2,ivdw)
        c = vdw%param(3,ivdw)

        kk=1.0_wp/c
        t1=Exp(-vdw%cutoff*kk)
        vdw%afs(ivdw) = (a+b*vdw%cutoff)*kk*t1-b*t1
        vdw%bfs(ivdw) = -(a*c+a*vdw%cutoff+b*vdw%cutoff*vdw%cutoff)*kk*t1

     Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

        z1 = vdw%param(1,ivdw)
        z2 = vdw%param(2,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        call zbl(vdw%cutoff,kk,a,z,dz)
        vdw%afs(ivdw) = dz/vdw%cutoff
        vdw%bfs(ivdw) = -z-dz

     Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = vdw%param(1,ivdw)
        z2 = vdw%param(2,ivdw)
        rm = vdw%param(3,ivdw)
        ic = 1.0_wp/vdw%param(4,ivdw)
        e0 = vdw%param(5,ivdw)
        r0 = vdw%param(6,ivdw)
        k = vdw%param(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0
        Call zbls(vdw%cutoff,kk,a,rm,ic,e0,k,r0,z,dz)
        vdw%afs(ivdw) = dz/vdw%cutoff
        vdw%bfs(ivdw) = -z-dz

     Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = vdw%param(1,ivdw)
        z2 = vdw%param(2,ivdw)
        rm = vdw%param(3,ivdw)
        ic = 1.0_wp/vdw%param(4,ivdw)
        e0 = vdw%param(5,ivdw)
        r0 = vdw%param(6,ivdw)
        k = vdw%param(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0
        Call zblb(vdw%cutoff,kk,a,rm,ic,e0,r0,k,z,dz)
        vdw%afs(ivdw) = dz/vdw%cutoff
        vdw%bfs(ivdw) = -z-dz

     Else

        Call error(150)

     End If

  End Do

End Subroutine vdw_direct_fs_generate


Subroutine vdw_table_read(vdw,site,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABLE file (for van der waals forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 1994
! amended   - i.t.todorov april 2016
! amended   - a.m.elena january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( vdw_type ), Intent( InOut ) :: vdw
  Type( site_type ), Intent( In    ) :: site
  Type( comms_type ), Intent( InOut ) :: comm

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail,ngrid,katom1,katom2,ivdw,jtpatm,keyvdw,i,j,l
  Real( Kind = wp )      :: delpot,cutpot,dlrpot,rdr,rrr,ppp,vk,vk1,vk2,t,t1,t2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message,messages(4)

  If (comm%idnode == 0) Open(Unit=ntable, File='TABLE')

! skip header record

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

! read mesh resolution

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

  Call get_word(record,word)
  delpot = word_2_real(word)

  Call get_word(record,word)
  cutpot = word_2_real(word)

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  dlrpot = vdw%cutoff/Real(vdw%max_grid-4,wp)

! check grid spacing

  safe=.false.
  If (Abs(delpot-dlrpot) <= 1.0e-8_wp) Then
     safe=.true.
     delpot=dlrpot
  End If
  If (delpot > delr_max .and. (.not.safe)) Then
     Write(messages(1),'(a,1p,e15.7)') 'expected (maximum) radial increment: ',delr_max
     Write(messages(2),'(a,1p,e15.7)') 'TABLE  file actual radial increment: ',delpot
     Write(messages(3),'(a,i10)') 'expected (minimum) number of grid points: ',vdw%max_grid
     Write(messages(4),'(a,i10)') 'TABLE  file actual number of grid points: ',ngrid
     Call info(messages,4,.true.)
     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     Write(message,'(a,i10)') 'TABLE arrays resized for mxgrid = ',vdw%max_grid-4
     Call info(message,.true.)
  End If

! compare grids dimensions

  If (ngrid < vdw%max_grid-4) Then
     Call warning(270,Real(ngrid,wp),Real(vdw%max_grid-4,wp),0.0_wp)
     Call error(48)
  End If

  If (cutpot < vdw%cutoff) Call error(504)

  fail=0
  Allocate (buffer(0:ngrid), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'vdw_table_read allocation failure'
     Call error(0,message)
  End If

! read potential arrays for all pairs

  Do ivdw=1,vdw%n_vdw

! read potential arrays if potential not already defined

     If (vdw%ltp(ivdw) == 0) Then

! read pair potential labels and long-range corrections

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Go To 100

        Call get_word(record,atom1)
        Call get_word(record,atom2)

        Call get_word(record,word)
        vdw%param(1,ivdw)=word_2_real(word)*engunit

        Call get_word(record,word)
        vdw%param(2,ivdw)=word_2_real(word)*engunit

        katom1=0
        katom2=0

        Do jtpatm=1,site%ntype_atom
           If (atom1 == site%unique_atom(jtpatm)) katom1=jtpatm
           If (atom2 == site%unique_atom(jtpatm)) katom2=jtpatm
        End Do

        If (katom1 == 0 .or. katom2 == 0) Then
           Write(message,'(a,i0,a,i0,a)') '****',atom1,'***',atom2,'**** entry in TABLE'
           Call error(81,message,.true.)
        End If

        keyvdw=(Max(katom1,katom2)*(Max(katom1,katom2)-1))/2 + Min(katom1,katom2)

! Only one vdw potential per pair is allowed
! (FIELD AND TABLE potentials overlapping)

        If (vdw%list(keyvdw) /= ivdw) Call error(23)

! read in potential arrays

        Do i=1,(ngrid+3)/4
           j=Min(4,ngrid-(i-1)*4)
           If (comm%idnode == 0) Then
              Read(Unit=ntable, Fmt=*, End=100) buffer((i-1)*4+1:(i-1)*4+j)
           Else
              buffer((i-1)*4+1:(i-1)*4+j)=0.0_wp
           End If
        End Do
        Call gbcast(comm,buffer,0)
! linear extrapolation for grid point 0 at distances close to 0

        vdw%tab_potential(0,ivdw) = 2.0_wp*buffer(1)-buffer(2)

! reconstruct arrays using 3pt interpolation

        If (remake) Then
           Do i=1,vdw%max_grid-4
              rrr = Real(i,wp)*dlrpot
              l   = Int(rrr*rdr)
              ppp=rrr*rdr-Real(l,wp)

              vk  = buffer(l)

! linear extrapolation for the grid points just beyond the cutoff

              If (l+2 > ngrid) Then
                 If (l+1 > ngrid) Then
                    vk1 = 2.0_wp*buffer(l)-buffer(l-1)
                    vk2 = 2.0_wp*vk1-buffer(l)
                 Else
                    vk1 = buffer(l+1)
                    vk2 = 2.0_wp*buffer(l+1)-buffer(l)
                 End If
              Else
                 vk1 = buffer(l+1)
                 vk2 = buffer(l+2)
              End If

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
              vdw%tab_potential(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,vdw%max_grid-4
              vdw%tab_potential(i,ivdw) = buffer(i)
           End Do

! linear extrapolation for the grid point just beyond the cutoff

           vdw%tab_potential(vdw%max_grid-3,ivdw) = 2.0_wp*vdw%tab_potential(vdw%max_grid-4,ivdw) - &
             vdw%tab_potential(vdw%max_grid-5,ivdw)
        End If

! linear extrapolation for the grid point at vdw%max_grid-2

        vdw%tab_potential(vdw%max_grid-2,ivdw) = 2.0_wp*vdw%tab_potential(vdw%max_grid-3,ivdw) - &
          vdw%tab_potential(vdw%max_grid-4,ivdw)

! read in force arrays

        Do i=1,(ngrid+3)/4
           j=Min(4,ngrid-(i-1)*4)
           If (comm%idnode == 0) Then
              Read(Unit=ntable, Fmt=*, End=100) buffer((i-1)*4+1:(i-1)*4+j)
           Else
              buffer((i-1)*4+1:(i-1)*4+j)=0.0_wp
           End If
        End Do
        Call gbcast(comm,buffer,0)
! linear extrapolation for grid point 0 at distances close to 0

        vdw%tab_force(0,ivdw) = (2.0_wp*buffer(1)-0.5_wp*buffer(2))/delpot

! reconstruct arrays using 3pt interpolation

        If (remake) Then
           Do i=1,vdw%max_grid-4
              rrr = Real(i,wp)*dlrpot
              l   = Int(rrr*rdr)
              ppp=rrr*rdr-Real(l,wp)

              vk  = buffer(l)

! linear extrapolation for the grid points just beyond the cutoff

              If (l+2 > ngrid) Then
                 If (l+1 > ngrid) Then
                    vk1 = 2.0_wp*buffer(l)-buffer(l-1)
                    vk2 = 2.0_wp*vk1-buffer(l)
                 Else
                    vk1 = buffer(l+1)
                    vk2 = 2.0_wp*buffer(l+1)-buffer(l)
                 End If
              Else
                 vk1 = buffer(l+1)
                 vk2 = buffer(l+2)
              End If

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              vdw%tab_force(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,vdw%max_grid-4
              vdw%tab_force(i,ivdw) = buffer(i)
           End Do

! linear extrapolation for the grid point just beyond the cutoff

           vdw%tab_force(vdw%max_grid-3,ivdw) = 2.0_wp*vdw%tab_force(vdw%max_grid-4,ivdw) - &
             vdw%tab_force(vdw%max_grid-5,ivdw)
        End If

! linear extrapolation for the grid point at vdw%max_grid-2

        vdw%tab_force(vdw%max_grid-2,ivdw) = 2.0_wp*vdw%tab_force(vdw%max_grid-3,ivdw) - &
          vdw%tab_force(vdw%max_grid-4,ivdw)

! We must distinguish that something has been defined

        If (Abs(vdw%tab_potential(0,ivdw)) <= zero_plus) Then
          vdw%tab_potential(0,ivdw) = Sign(Tiny(vdw%tab_potential(0,ivdw)),vdw%tab_potential(0,ivdw))
        End If

     End If

  End Do

  Call info('potential tables read from TABLE file',.true.)
  If (comm%idnode == 0) Then
     Close(Unit=ntable)
  End If

! convert to internal units

  Do ivdw=1,vdw%n_vdw
     If (vdw%ltp(ivdw) == 0) Then

! Sigma-epsilon initialisation

        vdw%sigeps(1,ivdw)=-1.0_wp
        vdw%sigeps(2,ivdw)= 0.0_wp

        Do i=0,vdw%max_grid
           vdw%tab_potential(i,ivdw)=vdw%tab_potential(i,ivdw)*engunit
           vdw%tab_force(i,ivdw)=vdw%tab_force(i,ivdw)*engunit

! Sigma-epsilon search

           If ((.not.vdw%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vdw%tab_potential(i-1,ivdw))) == -Nint(Sign(1.0_wp,vdw%tab_potential(i,ivdw)))) &
                    vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vdw%tab_potential(i-2,ivdw) >= vdw%tab_potential(i-1,ivdw) .and.  &
                       vdw%tab_potential(i-1,ivdw) <= vdw%tab_potential(i  ,ivdw)) .and. &
                      (vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i-1,ivdw) .or.   &
                       vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i  ,ivdw) .or.   &
                       vdw%tab_potential(i-1,ivdw) /= vdw%tab_potential(i  ,ivdw)) )     &
                    vdw%sigeps(2,ivdw)=-vdw%tab_potential(i-1,ivdw)
              End If
           End If
        End Do
     End If
  End Do

  If (vdw%l_force_shift) Then
     Do ivdw=1,vdw%n_vdw
        If (vdw%ltp(ivdw) == 0) Then

! Sigma-epsilon initialisation

           vdw%sigeps(1,ivdw)=-1.0_wp
           vdw%sigeps(2,ivdw)= 0.0_wp

! Sigma-epsilon search

           Do i=1,vdw%max_grid-4
              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 t  = vdw%tab_potential(i  ,ivdw) + vdw%tab_force(vdw%max_grid-4,ivdw) * &
                      (Real(i  ,wp)*dlrpot/vdw%cutoff-1.0_wp) - vdw%tab_potential(vdw%max_grid-4,ivdw)
                 t1 = vdw%tab_potential(i-1,ivdw) + vdw%tab_force(vdw%max_grid-4,ivdw) * &
                      (Real(i-1,wp)*dlrpot/vdw%cutoff-1.0_wp) - vdw%tab_potential(vdw%max_grid-4,ivdw)
                 If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,t1)) == -Nint(Sign(1.0_wp,t))) &
                       vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    t2 = vdw%tab_potential(i-2,ivdw) + vdw%tab_force(vdw%max_grid-4,ivdw) * &
                         (Real(i-2,wp)*dlrpot/vdw%cutoff-1.0_wp) - vdw%tab_potential(vdw%max_grid-4,ivdw)

                    If ( (t2 >= t1 .and. t1 <= t) .and.         &
                         (t2 /= t1 .or. t2 /= t .or. t1 /= t) ) &
                       vdw%sigeps(2,ivdw)=-t1
                 End If
              End If
           End Do
           vdw%tab_potential(vdw%max_grid-3,ivdw) = 0.0_wp ; vdw%tab_potential(vdw%max_grid-2,ivdw) = 0.0_wp
           vdw%tab_force(vdw%max_grid-3,ivdw) = 0.0_wp ; vdw%tab_force(vdw%max_grid-2,ivdw) = 0.0_wp
        End If
     End Do
  End If

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'vdw_table_read deallocation failure'
     Call error(0,message)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine vdw_table_read

Subroutine vdw_generate(vdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating potential energy and force arrays
! for van der waals forces only
!
! copyright - daresbury laboratory
! author    - w.smith may 1992
! amended   - i.t.todorov march 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (rydberg)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Type( vdw_type ), Intent( InOut ) :: vdw

  Integer           :: i,ivdw,keypot,n,m
  Real( Kind = wp ) :: dlrpot,r,r0,r0rn,r0rm,r_6,sor6,  &
                       rho,a,b,c,d,e0,kk,nr,mr,rc,sig,eps, &
                       alpha,beta,t1,t2,t3,t,z1,z2,dphi,phi, &
                       rm,k

! allocate arrays for tabulating

  Call vdw%init_table()

! define grid resolution for potential arrays

  dlrpot=vdw%cutoff/Real(vdw%max_grid-4,wp)

! construct arrays for all types of vdw potential

  Do ivdw=1,vdw%n_vdw

     keypot=vdw%ltp(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=vdw%param(1,ivdw)
        b=vdw%param(2,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           r_6=r**(-6)

           vdw%tab_potential(i,ivdw)=r_6*(a*r_6-b)
           vdw%tab_force(i,ivdw)=6.0_wp*r_6*(2.0_wp*a*r_6-b)
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           If (a*b > zero_plus) Then
              vdw%sigeps(1,ivdw)=(a/b)**(1.0_wp/6.0_wp)
              vdw%sigeps(2,ivdw)=b**2/(4.0_wp*a)
           End If ! else leave undetermined
        End If

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vdw%tab_potential(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)
           vdw%tab_force(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=sig
           vdw%sigeps(2,ivdw)=eps
        End If

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=vdw%param(1,ivdw)
        n =Nint(vdw%param(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(vdw%param(3,ivdw)) ; mr=Real(m,wp)
        r0=vdw%param(4,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           a=r0/r
           b=1.0_wp/(nr-mr)
           r0rn=(a)**n
           r0rm=(a)**m

           vdw%tab_potential(i,ivdw)=e0*(mr*r0rn-nr*r0rm)*b
           vdw%tab_force(i,ivdw)=e0*mr*nr*(r0rn-r0rm)*b
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=r0*(nr/mr)**(1.0_wp/(mr-nr))
           vdw%sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a  =vdw%param(1,ivdw)
        rho=vdw%param(2,ivdw)
        c  =vdw%param(3,ivdw)

        If (Abs(rho) <= zero_plus) Then
           If (Abs(a) <= zero_plus) Then
              rho=1.0_wp
           Else
              Call error(467)
           End If
        End If

! Sigma-epsilon initialisation

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=-1.0_wp
           vdw%sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           b=r/rho
           t1=a*Exp(-b)
           t2=-c/r**6

           vdw%tab_potential(i,ivdw)=t1+t2
           vdw%tab_force(i,ivdw)=t1*b+6.0_wp*t2

! Sigma-epsilon search

           If ((.not.vdw%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vdw%tab_potential(i-1,ivdw))) == -Nint(Sign(1.0_wp,vdw%tab_potential(i,ivdw)))) &
                    vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vdw%tab_potential(i-2,ivdw) >= vdw%tab_potential(i-1,ivdw) .and.  &
                       vdw%tab_potential(i-1,ivdw) <= vdw%tab_potential(i  ,ivdw)) .and. &
                      (vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i-1,ivdw) .or.   &
                       vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i  ,ivdw) .or.   &
                       vdw%tab_potential(i-1,ivdw) /= vdw%tab_potential(i  ,ivdw)) )     &
                    vdw%sigeps(2,ivdw)=-vdw%tab_potential(i-1,ivdw)
              End If
           End If
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =vdw%param(1,ivdw)
        b  =vdw%param(2,ivdw)
        sig=vdw%param(3,ivdw)
        c  =vdw%param(4,ivdw)
        d  =vdw%param(5,ivdw)

! Sigma-epsilon initialisation

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=-1.0_wp
           vdw%sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           t1=a*Exp(b*(sig-r))
           t2=-c/r**6
           t3=-d/r**8

           vdw%tab_potential(i,ivdw)=t1+t2+t3
           vdw%tab_force(i,ivdw)=t1*r*b+6.0_wp*t2+8.0_wp*t3

! Sigma-epsilon search

           If ((.not.vdw%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vdw%tab_potential(i-1,ivdw))) == -Nint(Sign(1.0_wp,vdw%tab_potential(i,ivdw)))) &
                    vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vdw%tab_potential(i-2,ivdw) >= vdw%tab_potential(i-1,ivdw) .and.  &
                       vdw%tab_potential(i-1,ivdw) <= vdw%tab_potential(i  ,ivdw)) .and. &
                      (vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i-1,ivdw) .or.   &
                       vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i  ,ivdw) .or.   &
                       vdw%tab_potential(i-1,ivdw) /= vdw%tab_potential(i  ,ivdw)) )     &
                    vdw%sigeps(2,ivdw)=-vdw%tab_potential(i-1,ivdw)
              End If
           End If
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=vdw%param(1,ivdw)
        b=vdw%param(2,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           t1=a/r**12
           t2=-b/r**10

           vdw%tab_potential(i,ivdw)=t1+t2
           vdw%tab_force(i,ivdw)=12.0_wp*t1+10.0_wp*t2
        End Do

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=Sqrt(a/b)
           vdw%sigeps(2,ivdw)=((b/6.0_wp)**6)*((5.0_wp/a)**5)
        End If
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

        e0=vdw%param(1,ivdw)
        n =Nint(vdw%param(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(vdw%param(3,ivdw)) ; mr=Real(m,wp)
        r0=vdw%param(4,ivdw)
        rc=vdw%param(5,ivdw) ; If (rc < 1.0e-6_wp) rc=vdw%cutoff

        If (n <= m) Call error(470)

! Sigma-epsilon initialisation

        vdw%sigeps(1,ivdw)=-1.0_wp
        vdw%sigeps(2,ivdw)= 0.0_wp

        t=Real(n-m,wp)

        b=1.0_wp/t
        c = rc/r0 ; If (c < 1.0_wp) Call error(468)

        beta = c*( (c**(m+1)-1.0_wp) / (c**(n+1)-1.0_wp) )**b
        alpha= -t / (  mr*(beta**n)*(1.0_wp+(nr/c-nr-1.0_wp)/c**n) &
                      -nr*(beta**m)*(1.0_wp+(mr/c-mr-1.0_wp)/c**m) )
        e0 = e0*alpha

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot
           If (r <= rc) Then
              a=r0/r

              vdw%tab_potential(i,ivdw)=e0*(  mr*(beta**n)*(a**n-(1.0_wp/c)**n) &
                                -nr*(beta**m)*(a**m-(1.0_wp/c)**m) &
                                +nr*mr*((r/rc-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
              vdw%tab_force(i,ivdw)=e0*mr*nr*( (beta**n)*a**n-(beta**m)*a**m &
                                    -r/rc*((beta/c)**n-(beta/c)**m) )*b

! Sigma-epsilon search

              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,vdw%tab_potential(i-1,ivdw))) == -Nint(Sign(1.0_wp,vdw%tab_potential(i,ivdw)))) &
                       vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    If ( (vdw%tab_potential(i-2,ivdw) >= vdw%tab_potential(i-1,ivdw) .and.  &
                          vdw%tab_potential(i-1,ivdw) <= vdw%tab_potential(i  ,ivdw)) .and. &
                         (vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i-1,ivdw) .or.   &
                          vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i  ,ivdw) .or.   &
                          vdw%tab_potential(i-1,ivdw) /= vdw%tab_potential(i  ,ivdw)) )     &
                       vdw%sigeps(2,ivdw)=-vdw%tab_potential(i-1,ivdw)
                 End If
              End If
           End If ! The else condition is satisfied by the vdw initialisation
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=vdw%param(1,ivdw)
        r0=vdw%param(2,ivdw)
        kk=vdw%param(3,ivdw)

        Do i=0,vdw%max_grid
           r=Real(i,wp)*dlrpot

           t1=Exp(-kk*(r-r0))

           vdw%tab_potential(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)
           vdw%tab_force(i,ivdw)=-2.0_wp*r*e0*kk*(1.0_wp-t1)*t1
        End Do
        t1=Exp(+kk*r0)
        vdw%tab_force(0,ivdw)=-2.0_wp*e0*kk*(1.0_wp-t1)*t1

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=r0-log(2.0_wp)/kk
           vdw%sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)
        d  =vdw%param(3,ivdw)

! Sigma-epsilon initialisation

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=-1.0_wp
           vdw%sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           If (r < vdw%param(4,ivdw) .or. Abs(r-d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6=(sig/(r-d))**6

              vdw%tab_potential(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
              vdw%tab_force(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*r/(r-d)

! Sigma-epsilon search

              If ((.not.vdw%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
                 If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,vdw%tab_potential(i-1,ivdw))) == -Nint(Sign(1.0_wp,vdw%tab_potential(i,ivdw)))) &
                       vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    If ( (vdw%tab_potential(i-2,ivdw) >= vdw%tab_potential(i-1,ivdw) .and.  &
                          vdw%tab_potential(i-1,ivdw) <= vdw%tab_potential(i  ,ivdw)) .and. &
                         (vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i-1,ivdw) .or.   &
                          vdw%tab_potential(i-2,ivdw) /= vdw%tab_potential(i  ,ivdw) .or.   &
                          vdw%tab_potential(i-1,ivdw) /= vdw%tab_potential(i  ,ivdw)) )     &
                       vdw%sigeps(2,ivdw)=-vdw%tab_potential(i-1,ivdw)
                 End If
              End If
           End If
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

     Else If (keypot == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.rc.(1-r/rc)^2

        a =vdw%param(1,ivdw)
        rc=vdw%param(2,ivdw)

        Do i=0,vdw%max_grid
           r=Real(i,wp)*dlrpot

           If (r < rc) Then
              t1=0.5_wp*a*rc
              t2=1.0_wp-r/rc

              vdw%tab_potential(i,ivdw)=t1*t2**2
              vdw%tab_force(i,ivdw)=a*t2*r
           End If
        End Do
        vdw%tab_force(0,ivdw)=a

        vdw%sigeps(1,ivdw)=rc
        vdw%sigeps(2,ivdw)=a

     Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           rho=r/sig
           t1=1.0_wp/(0.07_wp+rho)
           t2=1.0_wp/(0.12_wp+rho**7)
           t3=eps*(1.07_wp*t1)**7

           t=t3*((1.12_wp*t2) - 2.0_wp)

           vdw%tab_potential(i,ivdw)=t
           vdw%tab_force(i,ivdw)=7.0_wp*(t1*t + 1.12_wp*t3*t2**2*rho**6)*rho
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=sig
           vdw%sigeps(2,ivdw)=eps
        End If

      Else If (keypot == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps=vdw%param(1,ivdw)
        sig=vdw%param(2,ivdw)
        c  =vdw%param(3,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vdw%tab_potential(i,ivdw)=4.0_wp*eps*sor6*(sor6-c)
           vdw%tab_force(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-c)
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=sig
           vdw%sigeps(2,ivdw)=eps
        End If

     Else If (keypot == 13) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0 = vdw%param(1,ivdw)
        r0 = vdw%param(2,ivdw)
        kk = vdw%param(3,ivdw)
        c  = vdw%param(4,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           t1=Exp(-kk*(r - r0))
           sor6=c/r**12

           vdw%tab_potential(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)+sor6
           vdw%tab_force(i,ivdw)=-2.0_wp*r*e0*kk*(1.0_wp-t1)*t1+12.0_wp*sor6
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then !???
           vdw%sigeps(1,ivdw)=r0-log(2.0_wp)/kk
           vdw%sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)
        
        a = vdw%param(1,ivdw)
        b = vdw%param(2,ivdw)
        c = vdw%param(3,ivdw)

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           kk = r/c
           t1=Exp(-kk)           

           vdw%tab_potential(i,ivdw) = (a+b*r)*t1
           vdw%tab_force(i,ivdw) = t1*kk*(a-b*c+b*r)
        End Do
        vdw%tab_potential(0,ivdw)= a
        vdw%tab_force(0,ivdw)= 0

        If (.not.vdw%l_force_shift) Then !???
           vdw%sigeps(1,ivdw)=1.0_wp
           vdw%sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

        z1 = vdw%param(1,ivdw)
        z2 = vdw%param(2,ivdw)
        
        ! this is in fact inverse a
        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           call zbl(r,kk,a,phi,dphi)

           vdw%tab_potential(i,ivdw) = phi
           vdw%tab_force(i,ivdw) = dphi

        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=0.0_wp
           vdw%sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = vdw%param(1,ivdw)
        z2 = vdw%param(2,ivdw)
        rm = vdw%param(3,ivdw)
        c = 1.0_wp/vdw%param(4,ivdw)
        e0 = vdw%param(5,ivdw)
        r0 = vdw%param(6,ivdw)
        k = vdw%param(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           Call zbls(r,kk,a,rm,c,e0,k,r0,phi,dphi)
           vdw%tab_potential(i,ivdw) = phi
           vdw%tab_force(i,ivdw) = dphi
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=0.0_wp
           vdw%sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = vdw%param(1,ivdw)
        z2 = vdw%param(2,ivdw)
        rm = vdw%param(3,ivdw)
        c = 1.0_wp/vdw%param(4,ivdw)
        e0 = vdw%param(5,ivdw)
        r0 = vdw%param(6,ivdw)
        k = vdw%param(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,vdw%max_grid
           r=Real(i,wp)*dlrpot

           Call zblb(r,kk,a,rm,c,e0,r0,k,phi,dphi)
           vdw%tab_potential(i,ivdw) = phi
           vdw%tab_force(i,ivdw) = dphi
        End Do
        vdw%tab_potential(0,ivdw)=Huge(vdw%tab_potential(1,ivdw))
        vdw%tab_force(0,ivdw)=Huge(vdw%tab_force(1,ivdw))

        If (.not.vdw%l_force_shift) Then
           vdw%sigeps(1,ivdw)=0.0_wp
           vdw%sigeps(2,ivdw)=0.0_wp
        End If

     Else

        If (.not.vdw%l_tab) Call error(150)

     End If

     If (vdw%l_force_shift .and. (keypot /= 7 .and. keypot /= 10)) Then ! no shifting to shifted n-m and DPD

        vdw%sigeps(1,ivdw)=-1.0_wp
        vdw%sigeps(2,ivdw)= 0.0_wp

        Do i=1,vdw%max_grid-4
          t  = vdw%tab_potential(i  ,ivdw) + &
            vdw%tab_force(vdw%max_grid-4,ivdw)*(Real(i  ,wp)*dlrpot/vdw%cutoff-1.0_wp) - &
            vdw%tab_potential(vdw%max_grid-4,ivdw)
          t1 = vdw%tab_potential(i-1,ivdw) + &
            vdw%tab_force(vdw%max_grid-4,ivdw)*(Real(i-1,wp)*dlrpot/vdw%cutoff-1.0_wp) - &
            vdw%tab_potential(vdw%max_grid-4,ivdw)

          ! Sigma-epsilon search
          If (i > 20) Then ! Assumes some safety against numeric black holes!!!
            If (Sign(1.0_wp,vdw%sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
              If (Nint(Sign(1.0_wp,t1)) == -Nint(Sign(1.0_wp,t))) &
                vdw%sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
            Else                                           ! find epsilon
              t2 = vdw%tab_potential(i-2,ivdw) + &
                vdw%tab_force(vdw%max_grid-4,ivdw)*(Real(i-2,wp)*dlrpot/vdw%cutoff-1.0_wp) - &
                vdw%tab_potential(vdw%max_grid-4,ivdw)
              If ( (t2 >= t1 .and. t1 <= t) .and.         &
                (t2 /= t1 .or. t2 /= t .or. t1 /= t) ) &
                vdw%sigeps(2,ivdw)=-t1
            End If
          End If
        End Do
     End If

! Needed to distinguish that something has been defined

     If (Abs(vdw%tab_potential(0,ivdw)) <= zero_plus) Then
       vdw%tab_potential(0,ivdw) = Sign(Tiny(vdw%tab_potential(0,ivdw)),vdw%tab_potential(0,ivdw))
     End If

  End Do

End Subroutine vdw_generate


Subroutine vdw_forces &
           (iatm,xxt,yyt,zzt,rrt,engvdw,virvdw,stress,neigh,vdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating vdw energy and force terms using
! verlet neighbour list
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov march 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (rydberg)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,                                  Intent( In    ) :: iatm
  Type( neighbours_type ), Intent( In    ) :: neigh
  Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engvdw,virvdw
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Type( vdw_type ), Intent( InOut ) :: vdw

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: dlrpot,rdr

  Integer           :: mm,idi,ai,aj,jatm,key,k,l,ityp,n,m
  Real( Kind = wp ) :: rrr,ppp,gamma,eng,                 &
                       rsq,r_rrr,r_rsq,r_rvdw,r_rrv,rscl, &
                       r0,r0rn,r0rm,r_6,sor6,             &
                       rho,a,b,c,d,e0,kk,                 &
                       nr,mr,rc,sig,eps,alpha,beta,       &
                       fix,fiy,fiz,fx,fy,fz,              &
                       gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3,t,  &
                       strs1,strs2,strs3,strs5,strs6,     &
                       strs9,z1,z2,rm

! define grid resolution for potential arrays and interpolation spacing

  If (newjob) Then
     newjob = .false.

     dlrpot = vdw%cutoff/Real(vdw%max_grid-4,wp)
     rdr    = 1.0_wp/dlrpot
  End If

! initialise potential energy and virial

  engvdw=0.0_wp
  virvdw=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity and type of iatm

  idi=ltg(iatm)
  ai=ltype(iatm)

! load forces

  fix=fxx(iatm)
  fiy=fyy(iatm)
  fiz=fzz(iatm)

! start of primary loop for forces evaluation

  Do mm=1,neigh%list(0,iatm)

! atomic and potential function indices

     jatm=neigh%list(mm,iatm)
     aj=ltype(jatm)

     If (ai > aj) Then
        key=ai*(ai-1)/2 + aj
     Else
        key=aj*(aj-1)/2 + ai
     End If

     k=vdw%list(key)

! interatomic distance

     rrr = rrt(mm)

! validity and truncation of potential

     ityp=vdw%ltp(k)
     If (ityp >= 0 .and. rrr < vdw%cutoff) Then

! Distance derivatives

        r_rrr = 1.0_wp/rrr
        r_rvdw= 1.0_wp/vdw%cutoff
        rsq   = rrr**2
        r_rsq = 1.0_wp/rsq
        r_rrv = r_rrr*r_rvdw
        rscl  = rrr*r_rvdw

! Zero energy and force components

        eng   = 0.0_wp
        gamma = 0.0_wp

        If (vdw%l_direct) Then ! direct calculation

           If      (ityp == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=vdw%param(1,k)
              b=vdw%param(2,k)

              r_6=rrr**(-6)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = r_6*(a*r_6-b)
              gamma = 6.0_wp*r_6*(2.0_wp*a*r_6-b)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)

              sor6=(sig*r_rrr)**6

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

              e0=vdw%param(1,k)
              n =Nint(vdw%param(2,k)) ; nr=Real(n,wp)
              m =Nint(vdw%param(3,k)) ; mr=Real(m,wp)
              r0=vdw%param(4,k)

              a=r0*r_rrr
              b=1.0_wp/(nr-mr)
              r0rn=a**n
              r0rm=a**m

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*(mr*r0rn-nr*r0rm)*b
              gamma = e0*mr*nr*(r0rn-r0rm)*b*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              a  =vdw%param(1,k)
              rho=vdw%param(2,k)
              c  =vdw%param(3,k)

              If (Abs(rho) <= zero_plus) Then
                 If (Abs(a) <= zero_plus) Then
                    rho=1.0_wp
                 Else
                    Call error(467)
                 End If
              End If

              b=rrr/rho
              t1=a*Exp(-b)
              t2=-c*r_rrr**6

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2
              gamma = (t1*b+6.0_wp*t2)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              a  =vdw%param(1,k)
              b  =vdw%param(2,k)
              sig=vdw%param(3,k)
              c  =vdw%param(4,k)
              d  =vdw%param(5,k)

              t1=a*Exp(b*(sig-rrr))
              t2=-c*r_rrr**6
              t3=-d*r_rrr**8

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2+t3
              gamma = (t1*rrr*b+6.0_wp*t2+8.0_wp*t3)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=vdw%param(1,k)
              b=vdw%param(2,k)

              t1= a*r_rrr**12
              t2=-b*r_rrr**10

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2
              gamma = (12.0_wp*t1+10.0_wp*t2)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

              e0=vdw%param(1,k)
              n =Nint(vdw%param(2,k)) ; nr=Real(n,wp)
              m =Nint(vdw%param(3,k)) ; mr=Real(m,wp)
              r0=vdw%param(4,k)
              rc=vdw%param(5,k) ; If (rc < 1.0e-6_wp) rc=vdw%cutoff

              If (n <= m) Call error(470)

              t=Real(n-m,wp)

              b=1.0_wp/t
              c=rc/r0 ; If (c < 1.0_wp) Call error(468)

              beta = c*( (c**(m+1)-1.0_wp) / (c**(n+1)-1.0_wp) )**b
              alpha= -t / (  mr*(beta**n)*(1.0_wp+(nr/c-nr-1.0_wp)/c**n) &
                            -nr*(beta**m)*(1.0_wp+(mr/c-mr-1.0_wp)/c**m) )
              e0 = e0*alpha

              If (rrr <= rc) Then
                 a=r0*r_rrr

                 If (jatm <= natms .or. idi < ltg(jatm))            &
                    eng   = e0*(  mr*(beta**n)*(a**n-(1.0_wp/c)**n) &
                                 -nr*(beta**m)*(a**m-(1.0_wp/c)**m) &
                                 +nr*mr*((rrr/rc-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
                 gamma = e0*Real(m*n,wp)*(  (beta**n)*a**n-(beta**m)*a**m &
                                           -rrr/rc*((beta/c)**n-(beta/c)**m) )*b*r_rsq
              End If

           Else If (ityp == 8) Then

! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}

              e0=vdw%param(1,k)
              r0=vdw%param(2,k)
              kk=vdw%param(3,k)

              t1=Exp(-kk*(rrr-r0))

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*t1*(t1-2.0_wp)
              gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)*r_rrr

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (ityp == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)
              d  =vdw%param(3,k)

              If (rrr < vdw%param(4,k) .or. Abs(rrr-d) < 1.0e-10_wp) Then ! Else leave them zeros
                 sor6=(sig/(rrr-d))**6

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
                 gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rrr*(rrr-d))

                 If (vdw%l_force_shift) Then ! force-shifting
                    If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                    gamma = gamma - vdw%afs(k)*r_rrr
                 End If
              End If

           Else If (ityp == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

              a =vdw%param(1,k)
              rc=vdw%param(2,k)

              If (rrr < rc) Then ! Else leave them zeros
                 t2=rrr/rc
                 t1=0.5_wp*a*rrr*(1.0_wp-t2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = t1*(1.0_wp-t2)
                 gamma = t1*(3.0_wp*t2-1.0_wp)*r_rsq
              End If

           Else If (ityp == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)

              rho=rrr/sig
              t1=1.0_wp/(0.07_wp+rho)
              t2=1.0_wp/(0.12_wp+rho**7)
              t3=eps*(1.07_wp*t1)**7

              t=t3*((1.12_wp*t2) - 2.0_wp)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t
              gamma = 7.0_wp*(t1*t + 1.12_wp*t3*t2**2*rho**6)*rho*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

            Else If (ityp == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

              eps=vdw%param(1,k)
              sig=vdw%param(2,k)
              c  =vdw%param(3,k)

              sor6=(sig*r_rrr)**6

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = 4.0_wp*eps*sor6*(sor6-c)
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-c)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

            Else If (ityp == 13) Then

! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}+c/r^12

              e0=vdw%param(1,k)
              r0=vdw%param(2,k)
              kk=vdw%param(3,k)
              c=vdw%param(4,k)

              t1=Exp(-kk*(rrr-r0))
              sor6 = c*r_rrr**12
              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*t1*(t1-2.0_wp)+sor6
              gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)*r_rrr-12.0_wp*sor6*r_rrr

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

            Else If (ityp == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)

              a = vdw%param(1,k)
              b = vdw%param(2,k)
              c = vdw%param(3,k)

              kk = rrr/c
              t1 = Exp(-kk)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = (a+b*rrr)*t1
              gamma = kk*t1*(a-b*c+b*rrr)*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

            Else If (ityp == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

              z1 = vdw%param(1,k)
              z2 = vdw%param(2,k)
        
        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0

              Call zbl(rrr,kk,a,t1,gamma)
              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng = t1
              gamma = gamma*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

            Else If (ityp == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

              z1 = vdw%param(1,k)
              z2 = vdw%param(2,k)
              rm = vdw%param(3,k)
              c = 1.0_wp/vdw%param(4,k)
              e0 = vdw%param(5,k)
              r0 = vdw%param(6,k)
              t2 = vdw%param(7,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0

              Call zbls(rrr,kk,a,rm,c,e0,t2,r0,t1,gamma)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng = t1
              gamma = gamma*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

            Else If (ityp == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

              z1 = vdw%param(1,k)
              z2 = vdw%param(2,k)
              rm = vdw%param(3,k)
              c = 1.0_wp/vdw%param(4,k)
              e0 = vdw%param(5,k)
              r0 = vdw%param(6,k)
              t2 = vdw%param(7,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0

              Call zblb(rrr,kk,a,rm,c,e0,r0,t2,t1,gamma)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng = t1
              gamma = gamma*r_rsq

              If (vdw%l_force_shift) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + vdw%afs(k)*rrr + vdw%bfs(k)
                 gamma = gamma - vdw%afs(k)*r_rrr
              End If

           Else If (Abs(vdw%tab_potential(0,k)) > zero_plus) Then ! potential read from TABLE - (ityp == 0)

              l   = Int(rrr*rdr)
              ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

              If (jatm <= natms .or. idi < ltg(jatm)) Then
                 vk  = vdw%tab_potential(l,k)
                 vk1 = vdw%tab_potential(l+1,k)
                 vk2 = vdw%tab_potential(l+2,k)

                 t1 = vk  + (vk1 - vk )*ppp
                 t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                 eng = t1 + (t2-t1)*ppp*0.5_wp
                 ! force-shifting
                 If (vdw%l_force_shift) Then
                   eng = eng + vdw%tab_force(vdw%max_grid-4,k)*(rscl-1.0_wp) - &
                     vdw%tab_potential(vdw%max_grid-4,k)
                 End If
              End If

! calculate forces using 3-point interpolation

              gk  = vdw%tab_force(l,k) ; If (l == 0) gk = gk*rrr
              gk1 = vdw%tab_force(l+1,k)
              gk2 = vdw%tab_force(l+2,k)

              t1 = gk  + (gk1 - gk )*ppp
              t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

              gamma = (t1 + (t2-t1)*ppp*0.5_wp)*r_rsq
              If (vdw%l_force_shift) gamma = gamma - vdw%tab_force(vdw%max_grid-4,k)*r_rrv ! force-shifting

           End If

        Else If (Abs(vdw%tab_potential(0,k)) > zero_plus) Then ! no direct = fully tabulated calculation

           l   = Int(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

           If (jatm <= natms .or. idi < ltg(jatm)) Then
              vk  = vdw%tab_potential(l,k)
              vk1 = vdw%tab_potential(l+1,k)
              vk2 = vdw%tab_potential(l+2,k)

              t1 = vk  + (vk1 - vk )*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              eng = t1 + (t2-t1)*ppp*0.5_wp
              ! force-shifting
              If (vdw%l_force_shift) Then
                eng = eng + vdw%tab_force(vdw%max_grid-4,k)*(rscl-1.0_wp) - &
                  vdw%tab_potential(vdw%max_grid-4,k)
              End If
           End If

! calculate forces using 3-point interpolation

           gk  = vdw%tab_force(l,k) ; If (l == 0) gk = gk*rrr
           gk1 = vdw%tab_force(l+1,k)
           gk2 = vdw%tab_force(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)*r_rsq
           If (vdw%l_force_shift) gamma = gamma - vdw%tab_force(vdw%max_grid-4,k)*r_rrv ! force-shifting

        End If

! calculate forces

        fx = gamma*xxt(mm)
        fy = gamma*yyt(mm)
        fz = gamma*zzt(mm)

        fix=fix+fx
        fiy=fiy+fy
        fiz=fiz+fz

        If (jatm <= natms) Then

           fxx(jatm)=fxx(jatm)-fx
           fyy(jatm)=fyy(jatm)-fy
           fzz(jatm)=fzz(jatm)-fz

        End If

        If (jatm <= natms .or. idi < ltg(jatm)) Then

! add interaction energy

           engvdw = engvdw + eng

! add virial

           virvdw = virvdw - gamma*rsq

! add stress tensor

           strs1 = strs1 + xxt(mm)*fx
           strs2 = strs2 + xxt(mm)*fy
           strs3 = strs3 + xxt(mm)*fz
           strs5 = strs5 + yyt(mm)*fy
           strs6 = strs6 + yyt(mm)*fz
           strs9 = strs9 + zzt(mm)*fz

        End If

     End If

  End Do

! load back forces

  fxx(iatm)=fix
  fyy(iatm)=fiy
  fzz(iatm)=fiz

! complete stress tensor

  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9

End Subroutine vdw_forces

  Subroutine cleanup(T)
    Type( vdw_type ) :: T

    If (Allocated(T%list)) Then
      Deallocate(T%list)
    End If
    If (Allocated(T%ltp)) Then
      Deallocate(T%ltp)
    End If

    If (Allocated(T%param)) Then
      Deallocate(T%param)
    End If

    If (Allocated(T%sigeps)) Then
      Deallocate(T%sigeps)
    End If

    If (Allocated(T%tab_potential)) Then
      Deallocate(T%tab_potential)
    End If
    If (Allocated(T%tab_force)) Then
      Deallocate(T%tab_force)
    End If

    If (Allocated(T%afs)) Then
      Deallocate(T%afs)
    End If
    If (Allocated(T%bfs)) Then
      Deallocate(T%bfs)
    End If
  End Subroutine cleanup
End Module vdw
