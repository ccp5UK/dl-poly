Module vdw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global vdw interaction variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,  Only : comms_type,gsum
  Use setup_module
  Use site_module,   Only : ntpatm,numtyp
  Use configuration, Only : imcon,volm,natms,ltype,lfrzn, &
                            ltg,list,fxx,fyy,fzz
  Use mm3lrc
  Use m_zbl,         Only : ab, intRadZBL, intdRadZBL, &
                            zbl,zbls,zblb

  Use site_module,  Only : ntpatm,unqatm
  Use parse_module, Only : get_line,get_word,word_2_real
  Implicit None

  Logical,                        Save :: lt_vdw = .false., & ! no tabulated potentials are present
                                          ld_vdw = .false., & ! no direct calculations are opted
                                          ls_vdw = .false.    ! no force-shifting is opted

  Integer,                        Save :: ntpvdw = 0, &       ! number of 2 body interactions
                                          mxtvdw = 0          ! type of mixing


  Integer,           Allocatable, Save :: lstvdw(:),ltpvdw(:)

  Real( Kind = wp ), Allocatable, Save :: prmvdw(:,:),sigeps(:,:)

  Real( Kind = wp ),              Save :: elrc   = 0.0_wp, &
                                          virlrc = 0.0_wp

! Possible tabulated calculation arrays

  Real( Kind = wp ), Allocatable, Save :: vvdw(:,:),gvdw(:,:)

! Possible force-shifting arrays

  Real( Kind = wp ), Allocatable, Save :: afs(:),bfs(:)

  Public :: allocate_vdw_arrays, allocate_vdw_table_arrays, &
            allocate_vdw_direct_fs_arrays

Contains

  Subroutine allocate_vdw_arrays()

    Use setup_module, Only : mxvdw,mxpvdw

    Implicit None

    Integer, Dimension( 1:4 ) :: fail

    fail = 0

    Allocate (lstvdw(1:mxvdw),          Stat = fail(1))
    Allocate (ltpvdw(1:mxvdw),          Stat = fail(2))
    Allocate (prmvdw(1:mxpvdw,1:mxvdw), Stat = fail(3))
    Allocate (sigeps(1:2,1:mxvdw),      Stat = fail(4))

    If (Any(fail > 0)) Call error(1022)

    lstvdw = 0
    ltpvdw = 0

    prmvdw = 0.0_wp
    sigeps = 0.0_wp

  End Subroutine allocate_vdw_arrays

  Subroutine allocate_vdw_table_arrays()

    Use setup_module, Only : mxvdw,mxgvdw

    Implicit None

    Integer, Dimension( 1:2 ) :: fail

    fail = 0

    Allocate (vvdw(0:mxgvdw,1:mxvdw),   Stat = fail(1))
    Allocate (gvdw(0:mxgvdw,1:mxvdw),   Stat = fail(2))

    If (Any(fail > 0)) Call error(1063)

    vvdw = 0.0_wp
    gvdw = 0.0_wp

  End Subroutine allocate_vdw_table_arrays

  Subroutine allocate_vdw_direct_fs_arrays()

    Use setup_module, Only : mxvdw

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (afs(1:mxvdw),bfs(1:mxvdw), Stat = fail)

    If (fail > 0) Call error(1066)

    afs  = 0.0_wp
    bfs  = 0.0_wp

  End Subroutine allocate_vdw_direct_fs_arrays

  Subroutine vdw_lrc(rvdw,elrc,virlrc,comm)

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



  Real( Kind = wp ), Intent( In    ) :: rvdw
  Real( Kind = wp ), Intent(   Out ) :: elrc,virlrc
  Type( comms_type ), Intent( inOut ) :: comm

  Integer           :: fail,i,j,k,ivdw,keypot,n,m
  Real( Kind = wp ) :: a,b,c,d,e0,nr,mr,r0,r,eps,sig, &
                       eadd,padd,denprd,plrc,t,kk,s9, &
                       z1,z2,rm,al

  Real( Kind = wp ), Dimension( : ), Allocatable :: numfrz

  fail=0
  Allocate (numfrz(mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_lrc allocation failure, node: ', comm%idnode
     Call error(0)
  End If

! initialise long-range corrections to energy and pressure

  plrc = 0.0_wp
  elrc = 0.0_wp

  If (ls_vdw) Go To 10 ! force-shifting

! initialise counter arrays and evaluate number density in system

  numfrz = 0.0_wp
  Do i=1,natms
     k = ltype(i)
     If (lfrzn(i) /= 0) numfrz(k)=numfrz(k)+1.0_wp
  End Do
  Call gsum(comm,numfrz(1:ntpatm))

! Evaluate only for 3D periodic systems

  If (imcon /= 0 .and. imcon /= 6) Then
     ivdw = 0

     Do i=1,ntpatm
        Do j=1,i

           eadd = 0.0_wp
           padd = 0.0_wp

           ivdw = ivdw + 1
           k = lstvdw(ivdw)

           keypot=ltpvdw(k)
           If      (keypot ==  0) Then

! tabulated energy and pressure lrc

              eadd = prmvdw(1,k)
              padd =-prmvdw(2,k)

           Else If (keypot ==  1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=prmvdw(1,k)
              b=prmvdw(2,k)
              r=rvdw

              eadd = a/(9.0_wp*r**9) - b/(3.0_wp*r**3)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 6.0_wp*b/(3.0_wp*r**3)

           Else If (keypot ==  2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              r  =rvdw

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - sig**6/(3.0_wp*r**3))
              padd = 8.0_wp*eps*(6.0_wp*sig**12/(9.0_wp*r**9) - sig**6/(r**3))

           Else If (keypot ==  3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

              e0=prmvdw(1,k)
              n =Nint(prmvdw(2,k)) ; nr=Real(n,wp)
              m =Nint(prmvdw(3,k)) ; mr=Real(m,wp)
              r0=prmvdw(4,k)
              r =rvdw

              eadd = e0/(nr-mr)*( mr*r0**n/((nr-3.0_wp)*r**(n-3)) - nr*r0**m/((mr-3.0_wp)*r**(m-3)) )
              padd = e0/(nr-mr)*nr*mr*( r0**n/((nr-3.0_wp)*r**(n-3)) - r0**m/((mr-3.0_wp)*r**(m-3)) )

           Else If (keypot ==  4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              c=prmvdw(3,k)
              r=rvdw

              eadd = -c/(3.0_wp*r**3)
              padd = -2.0_wp*c/(r**3)

           Else If (keypot ==  5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              c=prmvdw(4,k)
              d=prmvdw(5,k)
              r=rvdw

              eadd = -c/(3.0_wp*r**3) - d/(5.0_wp*r**5)
              padd = -2.0_wp*c/(r**3) - 8.0_wp*d/(5.0_wp*r**5)

           Else If (keypot ==  6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=prmvdw(1,k)
              b=prmvdw(2,k)
              r=rvdw

              eadd = a/(9.0_wp*r**9) - b/(7.0_wp*r**7)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 10.0_wp*b/(7.0_wp*r**7)

           Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

              e0=prmvdw(1,k)
              r0=prmvdw(2,k)
              kk=prmvdw(3,k)
              If (kk > Tiny(kk)) Then
                 t = Exp(-kk*(rvdw - r0))

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*rvdw+1)**2 + 1) + &
                    e0*t*t/(4.0_wp*kk*kk*kk)*((kk*rvdw+1)**2 + kk*kk*rvdw*rvdw)
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*rvdw**3 + &
                      3*kk**2*rvdw**2 +6*kk*rvdw + 6) + &
                      e0*t*t/(4.0_wp*kk*kk*kk)* & 
                      (4.0_wp*kk**3*rvdw**3 + 6*kk**2*rvdw**2 + 6*kk*rvdw + 3)
              End If

           Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)

              a =0.07_wp
              b =0.12_wp
              e0=1.0e-12_wp

              eadd = intRadMM3(sig,a,b,eps,rvdw,e0)
              padd = -intRaddMM3(sig,a,b,eps,rvdw,e0)

           Else If (keypot ==  12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              c  =prmvdw(3,k)
              r  =rvdw

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - c*sig**6/(3.0_wp*r**3))
              padd = 8.0_wp*eps*(6.0_wp*sig**12/(9.0_wp*r**9) - c*sig**6/(r**3))

           Else If (keypot == 13) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

              e0 = prmvdw(1,k)
              r0 = prmvdw(2,k)
              kk = prmvdw(3,k)
               c = prmvdw(4,k)

              If (kk>Tiny(kk)) Then

                 t = Exp(-kk*(rvdw - r0))
                 s9 = c/(9.0_wp*rvdw**9)

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*rvdw+1)**2 + 1) + &
                     e0*t*t/(4.0_wp*kk*kk*kk)*((kk*rvdw+1)**2 + & 
                     kk*kk*rvdw*rvdw) + s9
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*rvdw**3 + & 
                       3*kk**2*rvdw**2 + 6*kk*rvdw + 6) + & 
                       e0*t*t/(4.0_wp*kk*kk*kk)* (4.0_wp*kk**3*rvdw**3 + & 
                       6*kk**2*rvdw**2 + 6*kk*rvdw + 3) + 12.0_wp*s9
              End If

           Else If (keypot == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)

              a = prmvdw(1,k)
              b = prmvdw(2,k)
              c = prmvdw(3,k)
              t = exp(-rvdw/c)

              eadd = (b*c*rvdw**3+(3*b*c**2+a*c)*rvdw**2+(6*b*c**3+2*a*c**2)*rvdw&
                +6*b*c**4+2*a*c**3)*t
              padd = (b*rvdw**4+(3*b*c+a)*rvdw**3+(9*b*c**2+3*a*c)*rvdw**2+& 
                (18*b*c**3+6*a*c**2)*rvdw+18*b*c**4+6*a*c**3)*t

           Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

              z1 = prmvdw(1,k)
              z2 = prmvdw(2,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0
              eadd = intRadZBL(kk,a,rvdw,1e-12_wp)
              padd = intdRadZBL(kk,a,rvdw,1e-12_wp)

           Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

              e0 = prmvdw(5,k)
              r0 = prmvdw(6,k)
              kk = prmvdw(7,k)

              If (kk > Tiny(kk)) Then
                 t = Exp(-kk*(rvdw - r0))

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*rvdw+1)**2 + 1) + &
                    e0*t*t/(4.0_wp*kk*kk*kk)*((kk*rvdw+1)**2 + kk*kk*rvdw*rvdw)
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*rvdw**3 + &
                      3*kk**2*rvdw**2 +6*kk*rvdw + 6) + &
                      e0*t*t/(4.0_wp*kk*kk*kk)* & 
                      (4.0_wp*kk**3*rvdw**3 + 6*kk**2*rvdw**2 + 6*kk*rvdw + 3)
              End If

           Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

              A = prmvdw(5,k)
              r0 = prmvdw(6,k)
              c = prmvdw(7,k)

              t=A*Exp(-rvdw/r0)

              eadd = (rvdw**2+2*r0*rvdw+2*r0**2)*t*r0-c/(3.0_wp*rvdw**3)
              padd = (rvdw**3+3*r0*rvdw**2+6*r0**2*rvdw+6*r0**3)*t -2.0_wp*c/(rvdw**3)

           End If

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

           If (i /= j) Then
              eadd = eadd*2.0_wp
              padd = padd*2.0_wp
           End If

           denprd=twopi * (numtyp(i)*numtyp(j) - numfrz(i)*numfrz(j)) / volm**2

           elrc = elrc + volm*denprd*eadd
           plrc = plrc + denprd*padd/3.0_wp

        End Do
     End Do

  End If

10 Continue

  If (comm%idnode == 0) Write(nrite,"(/,/,1x, &
     & 'long-range correction for: vdw energy  ',e15.6,/,26x, &
     & ': vdw pressure',e15.6)") elrc/engunit,plrc*prsunt

! convert plrc to a viral term

  virlrc = plrc*(-3.0_wp*volm)

  Deallocate (numfrz, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_lrc deallocation failure, node: ', comm%idnode
     Call error(0)
  End If

End Subroutine vdw_lrc

Subroutine vdw_direct_fs_generate(rvdw)

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


  Real( Kind = wp ), Intent( In    ) :: rvdw

  Integer           :: ivdw,keypot,n,m
  Real( Kind = wp ) :: r0,r0rn,r0rm,r_6,sor6,   &
                       rho,a,b,c,d,e0,kk,nr,mr, &
                       sig,eps,t1,t2,t3,z,dz,   &
                       z1,z2,rm,ic,k

! allocate arrays for force-shifted corrections

  Call allocate_vdw_direct_fs_arrays()

! construct arrays for all types of vdw potential

  Do ivdw=1,ntpvdw

     keypot=ltpvdw(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        r_6=rvdw**(-6)

        afs(ivdw) = 6.0_wp*r_6*(2.0_wp*a*r_6-b)
        bfs(ivdw) =-r_6*(a*r_6-b) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        sor6=(sig/rvdw)**6

        afs(ivdw) = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        bfs(ivdw) =-4.0_wp*eps*sor6*(sor6-1.0_wp) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=prmvdw(1,ivdw)
        n =Nint(prmvdw(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(prmvdw(3,ivdw)) ; mr=Real(m,wp)
        r0=prmvdw(4,ivdw)

        a=r0/rvdw
        b=1.0_wp/Real(n-m,wp)
        r0rn=a**n
        r0rm=a**m

        afs(ivdw) = e0*mr*nr*(r0rn-r0rm)*b
        bfs(ivdw) =-e0*(mr*r0rn-nr*r0rm)*b - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a  =prmvdw(1,ivdw)
        rho=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        If (Abs(rho) <= zero_plus) Then
           If (Abs(a) <= zero_plus) Then
              rho=1.0_wp
           Else
              Call error(467)
           End If
        End If

        b=rvdw/rho
        t1=a*Exp(-b)
        t2=-c/rvdw**6

        afs(ivdw) = (t1*b+6.0_wp*t2)
        bfs(ivdw) =-(t1+t2) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =prmvdw(1,ivdw)
        b  =prmvdw(2,ivdw)
        sig=prmvdw(3,ivdw)
        c  =prmvdw(4,ivdw)
        d  =prmvdw(5,ivdw)

        t1=a*Exp(b*(sig-rvdw))
        t2=-c/rvdw**6
        t3=-d/rvdw**8

        afs(ivdw) = (t1*rvdw*b+6.0_wp*t2+8.0_wp*t3)
        bfs(ivdw) =-(t1+t2+t3) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        t1=a/rvdw**12
        t2=-b/rvdw**10

        afs(ivdw) = (12.0_wp*t1+10.0_wp*t2)
        bfs(ivdw) =-(t1+t2) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        kk=prmvdw(3,ivdw)

        t1=Exp(-kk*(rvdw-r0))

        afs(ivdw) =-2.0_wp*e0*kk*t1*(1.0_wp-t1)*rvdw
        bfs(ivdw) =-e0*t1*(t1-2.0_wp) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        d  =prmvdw(3,ivdw)

        sor6=(sig/(rvdw-d))**6

        afs(ivdw) = (24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rvdw-d))*rvdw
        bfs(ivdw) =-(4.0_wp*eps*sor6*(sor6-1.0_wp)+eps) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 10) Then ! all zeroed in vdw_module

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

!       afs(ivdw) = 0.0_wp !initialised in vdw_module
!       bfs(ivdw) = 0.0_wp !initialised in vdw_module

     Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        rho=sig/rvdw
        t1=1.0_wp/(0.07_wp+rho)
        t2=1.0_wp/(0.12_wp+rho**7)
        t3=eps*(1.07_wp/t1**7)

        afs(ivdw) =-7.0_wp*t3*rho*(((1.12_wp/t2)-2.0_wp)/t1 + (1.12_wp/t2**2)*rho**6)
        bfs(ivdw) =-t3*((1.12_wp/t2)-2.0_wp) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

      Else If (keypot == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        sor6=(sig/rvdw)**6

        afs(ivdw) = 24.0_wp*eps*sor6*(2.0_wp*sor6-c)
        bfs(ivdw) =-4.0_wp*eps*sor6*(sor6-c) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 13) Then

! Morse potential with twelve term:: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        kk=prmvdw(3,ivdw)
        c=prmvdw(4,ivdw)

        t1=Exp(-kk*(rvdw-r0))
        sor6 = c/rvdw**12

        afs(ivdw) =-2.0_wp*e0*kk*t1*(1.0_wp-t1)*rvdw + 12.0_wp*sor6
        bfs(ivdw) =-e0*t1*(t1-2.0_wp) + sor6 - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 14) Then

! Morse potential with twelve term:: u=(a+b*r)Exp(-r/c)

        a = prmvdw(1,ivdw)
        b = prmvdw(2,ivdw)
        c = prmvdw(3,ivdw)

        kk=1.0_wp/c
        t1=Exp(-rvdw*kk)
        afs(ivdw) = (a+b*rvdw)*kk*t1-b*t1
        bfs(ivdw) = -(a*c+a*rvdw+b*rvdw*rvdw)*kk*t1

     Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        call zbl(rvdw,kk,a,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        ic = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0
        Call zbls(rvdw,kk,a,rm,ic,e0,k,r0,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        ic = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0
        Call zblb(rvdw,kk,a,rm,ic,e0,r0,k,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else

        Call error(150)

     End If

  End Do

End Subroutine vdw_direct_fs_generate


Subroutine vdw_table_read(rvdw,comm)

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


  Real( Kind = wp ), Intent( In    ) :: rvdw

  Type( comms_type ), Intent( InOut ) :: comm
  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail,ngrid,katom1,katom2,ivdw,jtpatm,keyvdw,i,j,l
  Real( Kind = wp )      :: delpot,cutpot,dlrpot,rdr,rrr,ppp,vk,vk1,vk2,t,t1,t2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer


  If (comm%idnode == 0) Open(Unit=ntable, File='TABLE')

! skip header record

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

! read mesh resolution

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

  Call get_word(record,word)
  delpot = word_2_real(word,comm)

  Call get_word(record,word)
  cutpot = word_2_real(word,comm)

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word,comm))

  dlrpot = rvdw/Real(mxgvdw-4,wp)

! check grid spacing

  safe=.false.
  If (Abs(delpot-dlrpot) <= 1.0e-8_wp) Then
     safe=.true.
     delpot=dlrpot
  End If
  If (delpot > delr_max .and. (.not.safe)) Then
     If (comm%idnode == 0) Then
        Write(nrite,"(/,                                             &
             & ' expected (maximum) radial increment : ',1p,e15.7,/, &
             & ' TABLE  file actual radial increment : ',1p,e15.7)") &
             delr_max, delpot
        Write(nrite,"(/,                                                &
             & ' expected (minimum) number of grid points : ',0p,i10,/, &
             & ' TABLE  file actual number of grid points : ',0p,i10)") &
             mxgvdw, ngrid
     End If
     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     If (comm%idnode == 0) Write(nrite,"(/,' TABLE arrays resized for mxgrid = ',i10)") mxgvdw-4
  End If

! compare grids dimensions

  If (ngrid < mxgvdw-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgvdw-4,wp),0.0_wp)
     Call error(48)
  End If

  If (cutpot < rvdw) Call error(504)

  fail=0
  Allocate (buffer(0:ngrid), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_table_read allocation failure, node: ', comm%idnode
     Call error(0)
  End If

! read potential arrays for all pairs

  Do ivdw=1,ntpvdw

! read potential arrays if potential not already defined

     If (ltpvdw(ivdw) == 0) Then

! read pair potential labels and long-range corrections

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Go To 100

        Call get_word(record,atom1)
        Call get_word(record,atom2)

        Call get_word(record,word)
        prmvdw(1,ivdw)=word_2_real(word,comm)*engunit

        Call get_word(record,word)
        prmvdw(2,ivdw)=word_2_real(word,comm)*engunit

        katom1=0
        katom2=0

        Do jtpatm=1,ntpatm
           If (atom1 == unqatm(jtpatm)) katom1=jtpatm
           If (atom2 == unqatm(jtpatm)) katom2=jtpatm
        End Do

        If (katom1 == 0 .or. katom2 == 0) Then
           If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'**** entry in TABLE'
           Call error(81)
        End If

        keyvdw=(Max(katom1,katom2)*(Max(katom1,katom2)-1))/2 + Min(katom1,katom2)

! Only one vdw potential per pair is allowed
! (FIELD AND TABLE potentials overlapping)

        If (lstvdw(keyvdw) /= ivdw) Call error(23)

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

        vvdw(0,ivdw) = 2.0_wp*buffer(1)-buffer(2)

! reconstruct arrays using 3pt interpolation

        If (remake) Then
           Do i=1,mxgvdw-4
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
              vvdw(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,mxgvdw-4
              vvdw(i,ivdw) = buffer(i)
           End Do

! linear extrapolation for the grid point just beyond the cutoff

           vvdw(mxgvdw-3,ivdw) = 2.0_wp*vvdw(mxgvdw-4,ivdw) - vvdw(mxgvdw-5,ivdw)
        End If

! linear extrapolation for the grid point at mxgvdw-2

        vvdw(mxgvdw-2,ivdw) = 2.0_wp*vvdw(mxgvdw-3,ivdw) - vvdw(mxgvdw-4,ivdw)

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

        gvdw(0,ivdw) = (2.0_wp*buffer(1)-0.5_wp*buffer(2))/delpot

! reconstruct arrays using 3pt interpolation

        If (remake) Then
           Do i=1,mxgvdw-4
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

              gvdw(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,mxgvdw-4
              gvdw(i,ivdw) = buffer(i)
           End Do

! linear extrapolation for the grid point just beyond the cutoff

           gvdw(mxgvdw-3,ivdw) = 2.0_wp*gvdw(mxgvdw-4,ivdw) - gvdw(mxgvdw-5,ivdw)
        End If

! linear extrapolation for the grid point at mxgvdw-2

        gvdw(mxgvdw-2,ivdw) = 2.0_wp*gvdw(mxgvdw-3,ivdw) - gvdw(mxgvdw-4,ivdw)

! We must distinguish that something has been defined

        If (Abs(vvdw(0,ivdw)) <= zero_plus) vvdw(0,ivdw) = Sign(Tiny(vvdw(0,ivdw)),vvdw(0,ivdw))

     End If

  End Do

  If (comm%idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABLE file'
  End If

! convert to internal units

  Do ivdw=1,ntpvdw
     If (ltpvdw(ivdw) == 0) Then

! Sigma-epsilon initialisation

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        Do i=0,mxgvdw
           vvdw(i,ivdw)=vvdw(i,ivdw)*engunit
           gvdw(i,ivdw)=gvdw(i,ivdw)*engunit

! Sigma-epsilon search

           If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                       vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                      (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                       vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                       vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                    sigeps(2,ivdw)=-vvdw(i-1,ivdw)
              End If
           End If
        End Do
     End If
  End Do

  If (ls_vdw) Then
     Do ivdw=1,ntpvdw
        If (ltpvdw(ivdw) == 0) Then

! Sigma-epsilon initialisation

           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp

! Sigma-epsilon search

           Do i=1,mxgvdw-4
              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 t  = vvdw(i  ,ivdw) + gvdw(mxgvdw-4,ivdw) * &
                      (Real(i  ,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
                 t1 = vvdw(i-1,ivdw) + gvdw(mxgvdw-4,ivdw) * &
                      (Real(i-1,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,t1)) == -Nint(Sign(1.0_wp,t))) &
                       sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    t2 = vvdw(i-2,ivdw) + gvdw(mxgvdw-4,ivdw) * &
                         (Real(i-2,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)

                    If ( (t2 >= t1 .and. t1 <= t) .and.         &
                         (t2 /= t1 .or. t2 /= t .or. t1 /= t) ) &
                       sigeps(2,ivdw)=-t1
                 End If
              End If
           End Do
           vvdw(mxgvdw-3,ivdw) = 0.0_wp ; vvdw(mxgvdw-2,ivdw) = 0.0_wp
           gvdw(mxgvdw-3,ivdw) = 0.0_wp ; gvdw(mxgvdw-2,ivdw) = 0.0_wp
        End If
     End Do
  End If

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_table_read deallocation failure, node: ', comm%idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine vdw_table_read

Subroutine vdw_generate(rvdw)

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

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Integer           :: i,ivdw,keypot,n,m
  Real( Kind = wp ) :: dlrpot,r,r0,r0rn,r0rm,r_6,sor6,  &
                       rho,a,b,c,d,e0,kk,nr,mr,rc,sig,eps, &
                       alpha,beta,t1,t2,t3,t,z1,z2,dphi,phi, &
                       rm,k

! allocate arrays for tabulating

  Call allocate_vdw_table_arrays()

! define grid resolution for potential arrays

  dlrpot=rvdw/Real(mxgvdw-4,wp)

! construct arrays for all types of vdw potential

  Do ivdw=1,ntpvdw

     keypot=ltpvdw(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           r_6=r**(-6)

           vvdw(i,ivdw)=r_6*(a*r_6-b)
           gvdw(i,ivdw)=6.0_wp*r_6*(2.0_wp*a*r_6-b)
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           If (a*b > zero_plus) Then
              sigeps(1,ivdw)=(a/b)**(1.0_wp/6.0_wp)
              sigeps(2,ivdw)=b**2/(4.0_wp*a)
           End If ! else leave undetermined
        End If

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)
           gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=sig
           sigeps(2,ivdw)=eps
        End If

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=prmvdw(1,ivdw)
        n =Nint(prmvdw(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(prmvdw(3,ivdw)) ; mr=Real(m,wp)
        r0=prmvdw(4,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           a=r0/r
           b=1.0_wp/(nr-mr)
           r0rn=(a)**n
           r0rm=(a)**m

           vvdw(i,ivdw)=e0*(mr*r0rn-nr*r0rm)*b
           gvdw(i,ivdw)=e0*mr*nr*(r0rn-r0rm)*b
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=r0*(nr/mr)**(1.0_wp/(mr-nr))
           sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a  =prmvdw(1,ivdw)
        rho=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        If (Abs(rho) <= zero_plus) Then
           If (Abs(a) <= zero_plus) Then
              rho=1.0_wp
           Else
              Call error(467)
           End If
        End If

! Sigma-epsilon initialisation

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           b=r/rho
           t1=a*Exp(-b)
           t2=-c/r**6

           vvdw(i,ivdw)=t1+t2
           gvdw(i,ivdw)=t1*b+6.0_wp*t2

! Sigma-epsilon search

           If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                       vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                      (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                       vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                       vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                    sigeps(2,ivdw)=-vvdw(i-1,ivdw)
              End If
           End If
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =prmvdw(1,ivdw)
        b  =prmvdw(2,ivdw)
        sig=prmvdw(3,ivdw)
        c  =prmvdw(4,ivdw)
        d  =prmvdw(5,ivdw)

! Sigma-epsilon initialisation

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=a*Exp(b*(sig-r))
           t2=-c/r**6
           t3=-d/r**8

           vvdw(i,ivdw)=t1+t2+t3
           gvdw(i,ivdw)=t1*r*b+6.0_wp*t2+8.0_wp*t3

! Sigma-epsilon search

           If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                       vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                      (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                       vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                       vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                    sigeps(2,ivdw)=-vvdw(i-1,ivdw)
              End If
           End If
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=a/r**12
           t2=-b/r**10

           vvdw(i,ivdw)=t1+t2
           gvdw(i,ivdw)=12.0_wp*t1+10.0_wp*t2
        End Do

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=Sqrt(a/b)
           sigeps(2,ivdw)=((b/6.0_wp)**6)*((5.0_wp/a)**5)
        End If
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

        e0=prmvdw(1,ivdw)
        n =Nint(prmvdw(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(prmvdw(3,ivdw)) ; mr=Real(m,wp)
        r0=prmvdw(4,ivdw)
        rc=prmvdw(5,ivdw) ; If (rc < 1.0e-6_wp) rc=rvdw

        If (n <= m) Call error(470)

! Sigma-epsilon initialisation

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        t=Real(n-m,wp)

        b=1.0_wp/t
        c = rc/r0 ; If (c < 1.0_wp) Call error(468)

        beta = c*( (c**(m+1)-1.0_wp) / (c**(n+1)-1.0_wp) )**b
        alpha= -t / (  mr*(beta**n)*(1.0_wp+(nr/c-nr-1.0_wp)/c**n) &
                      -nr*(beta**m)*(1.0_wp+(mr/c-mr-1.0_wp)/c**m) )
        e0 = e0*alpha

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot
           If (r <= rc) Then
              a=r0/r

              vvdw(i,ivdw)=e0*(  mr*(beta**n)*(a**n-(1.0_wp/c)**n) &
                                -nr*(beta**m)*(a**m-(1.0_wp/c)**m) &
                                +nr*mr*((r/rc-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
              gvdw(i,ivdw)=e0*mr*nr*( (beta**n)*a**n-(beta**m)*a**m &
                                    -r/rc*((beta/c)**n-(beta/c)**m) )*b

! Sigma-epsilon search

              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                       sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                          vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                         (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                          vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                          vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                       sigeps(2,ivdw)=-vvdw(i-1,ivdw)
                 End If
              End If
           End If ! The else condition is satisfied by the vdw_module initialisation
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        kk=prmvdw(3,ivdw)

        Do i=0,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=Exp(-kk*(r-r0))

           vvdw(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)
           gvdw(i,ivdw)=-2.0_wp*r*e0*kk*(1.0_wp-t1)*t1
        End Do
        t1=Exp(+kk*r0)
        gvdw(0,ivdw)=-2.0_wp*e0*kk*(1.0_wp-t1)*t1

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=r0-log(2.0_wp)/kk
           sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        d  =prmvdw(3,ivdw)

! Sigma-epsilon initialisation

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           If (r < prmvdw(4,ivdw) .or. Abs(r-d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6=(sig/(r-d))**6

              vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
              gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*r/(r-d)

! Sigma-epsilon search

              If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                       sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                          vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                         (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                          vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                          vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                       sigeps(2,ivdw)=-vvdw(i-1,ivdw)
                 End If
              End If
           End If
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.rc.(1-r/rc)^2

        a =prmvdw(1,ivdw)
        rc=prmvdw(2,ivdw)

        Do i=0,mxgvdw
           r=Real(i,wp)*dlrpot

           If (r < rc) Then
              t1=0.5_wp*a*rc
              t2=1.0_wp-r/rc

              vvdw(i,ivdw)=t1*t2**2
              gvdw(i,ivdw)=a*t2*r
           End If
        End Do
        gvdw(0,ivdw)=a

        sigeps(1,ivdw)=rc
        sigeps(2,ivdw)=a

     Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           rho=r/sig
           t1=1.0_wp/(0.07_wp+rho)
           t2=1.0_wp/(0.12_wp+rho**7)
           t3=eps*(1.07_wp*t1)**7

           t=t3*((1.12_wp*t2) - 2.0_wp)

           vvdw(i,ivdw)=t
           gvdw(i,ivdw)=7.0_wp*(t1*t + 1.12_wp*t3*t2**2*rho**6)*rho
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=sig
           sigeps(2,ivdw)=eps
        End If

      Else If (keypot == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-c)
           gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-c)
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=sig
           sigeps(2,ivdw)=eps
        End If

     Else If (keypot == 13) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0 = prmvdw(1,ivdw)
        r0 = prmvdw(2,ivdw)
        kk = prmvdw(3,ivdw)
        c  = prmvdw(4,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=Exp(-kk*(r - r0))
           sor6=c/r**12

           vvdw(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)+sor6
           gvdw(i,ivdw)=-2.0_wp*r*e0*kk*(1.0_wp-t1)*t1+12.0_wp*sor6
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then !???
           sigeps(1,ivdw)=r0-log(2.0_wp)/kk
           sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)
        
        a = prmvdw(1,ivdw)
        b = prmvdw(2,ivdw)
        c = prmvdw(3,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           kk = r/c
           t1=Exp(-kk)           

           vvdw(i,ivdw) = (a+b*r)*t1
           gvdw(i,ivdw) = t1*kk*(a-b*c+b*r)
        End Do
        vvdw(0,ivdw)= a
        gvdw(0,ivdw)= 0

        If (.not.ls_vdw) Then !???
           sigeps(1,ivdw)=1.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        
        ! this is in fact inverse a
        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           call zbl(r,kk,a,phi,dphi)

           vvdw(i,ivdw) = phi
           gvdw(i,ivdw) = dphi

        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=0.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        c = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           Call zbls(r,kk,a,rm,c,e0,k,r0,phi,dphi)
           vvdw(i,ivdw) = phi
           gvdw(i,ivdw) = dphi
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=0.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        c = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           Call zblb(r,kk,a,rm,c,e0,r0,k,phi,dphi)
           vvdw(i,ivdw) = phi
           gvdw(i,ivdw) = dphi
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=0.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else

        If (.not.lt_vdw) Call error(150)

     End If

     If (ls_vdw .and. (keypot /= 7 .and. keypot /= 10)) Then ! no shifting to shifted n-m and DPD

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        Do i=1,mxgvdw-4
           t  = vvdw(i  ,ivdw) + gvdw(mxgvdw-4,ivdw)*(Real(i  ,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
           t1 = vvdw(i-1,ivdw) + gvdw(mxgvdw-4,ivdw)*(Real(i-1,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)

! Sigma-epsilon search

           If (i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,t1)) == -Nint(Sign(1.0_wp,t))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 t2 = vvdw(i-2,ivdw) + gvdw(mxgvdw-4,ivdw)*(Real(i-2,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
                 If ( (t2 >= t1 .and. t1 <= t) .and.         &
                      (t2 /= t1 .or. t2 /= t .or. t1 /= t) ) &
                    sigeps(2,ivdw)=-t1
              End If
           End If
        End Do
     End If

! Needed to distinguish that something has been defined

     If (Abs(vvdw(0,ivdw)) <= zero_plus) vvdw(0,ivdw) = Sign(Tiny(vvdw(0,ivdw)),vvdw(0,ivdw))

  End Do

End Subroutine vdw_generate


Subroutine vdw_forces &
           (iatm,rvdw,xxt,yyt,zzt,rrt,engvdw,virvdw,stress)

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
  Real( Kind = wp ),                        Intent( In    ) :: rvdw
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engvdw,virvdw
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

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

     dlrpot = rvdw/Real(mxgvdw-4,wp)
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

  Do mm=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(mm,iatm)
     aj=ltype(jatm)

     If (ai > aj) Then
        key=ai*(ai-1)/2 + aj
     Else
        key=aj*(aj-1)/2 + ai
     End If

     k=lstvdw(key)

! interatomic distance

     rrr = rrt(mm)

! validity and truncation of potential

     ityp=ltpvdw(k)
     If (ityp >= 0 .and. rrr < rvdw) Then

! Distance derivatives

        r_rrr = 1.0_wp/rrr
        r_rvdw= 1.0_wp/rvdw
        rsq   = rrr**2
        r_rsq = 1.0_wp/rsq
        r_rrv = r_rrr*r_rvdw
        rscl  = rrr*r_rvdw

! Zero energy and force components

        eng   = 0.0_wp
        gamma = 0.0_wp

        If (ld_vdw) Then ! direct calculation

           If      (ityp == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=prmvdw(1,k)
              b=prmvdw(2,k)

              r_6=rrr**(-6)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = r_6*(a*r_6-b)
              gamma = 6.0_wp*r_6*(2.0_wp*a*r_6-b)*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)

              sor6=(sig*r_rrr)**6

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

              e0=prmvdw(1,k)
              n =Nint(prmvdw(2,k)) ; nr=Real(n,wp)
              m =Nint(prmvdw(3,k)) ; mr=Real(m,wp)
              r0=prmvdw(4,k)

              a=r0*r_rrr
              b=1.0_wp/(nr-mr)
              r0rn=a**n
              r0rm=a**m

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*(mr*r0rn-nr*r0rm)*b
              gamma = e0*mr*nr*(r0rn-r0rm)*b*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              a  =prmvdw(1,k)
              rho=prmvdw(2,k)
              c  =prmvdw(3,k)

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

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              a  =prmvdw(1,k)
              b  =prmvdw(2,k)
              sig=prmvdw(3,k)
              c  =prmvdw(4,k)
              d  =prmvdw(5,k)

              t1=a*Exp(b*(sig-rrr))
              t2=-c*r_rrr**6
              t3=-d*r_rrr**8

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2+t3
              gamma = (t1*rrr*b+6.0_wp*t2+8.0_wp*t3)*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=prmvdw(1,k)
              b=prmvdw(2,k)

              t1= a*r_rrr**12
              t2=-b*r_rrr**10

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2
              gamma = (12.0_wp*t1+10.0_wp*t2)*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

              e0=prmvdw(1,k)
              n =Nint(prmvdw(2,k)) ; nr=Real(n,wp)
              m =Nint(prmvdw(3,k)) ; mr=Real(m,wp)
              r0=prmvdw(4,k)
              rc=prmvdw(5,k) ; If (rc < 1.0e-6_wp) rc=rvdw

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

              e0=prmvdw(1,k)
              r0=prmvdw(2,k)
              kk=prmvdw(3,k)

              t1=Exp(-kk*(rrr-r0))

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*t1*(t1-2.0_wp)
              gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)*r_rrr

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (ityp == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              d  =prmvdw(3,k)

              If (rrr < prmvdw(4,k) .or. Abs(rrr-d) < 1.0e-10_wp) Then ! Else leave them zeros
                 sor6=(sig/(rrr-d))**6

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
                 gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rrr*(rrr-d))

                 If (ls_vdw) Then ! force-shifting
                    If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng   = eng + afs(k)*rrr + bfs(k)
                    gamma = gamma - afs(k)*r_rrr
                 End If
              End If

           Else If (ityp == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

              a =prmvdw(1,k)
              rc=prmvdw(2,k)

              If (rrr < rc) Then ! Else leave them zeros
                 t2=rrr/rc
                 t1=0.5_wp*a*rrr*(1.0_wp-t2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = t1*(1.0_wp-t2)
                 gamma = t1*(3.0_wp*t2-1.0_wp)*r_rsq
              End If

           Else If (ityp == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)

              rho=rrr/sig
              t1=1.0_wp/(0.07_wp+rho)
              t2=1.0_wp/(0.12_wp+rho**7)
              t3=eps*(1.07_wp*t1)**7

              t=t3*((1.12_wp*t2) - 2.0_wp)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t
              gamma = 7.0_wp*(t1*t + 1.12_wp*t3*t2**2*rho**6)*rho*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

            Else If (ityp == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              c  =prmvdw(3,k)

              sor6=(sig*r_rrr)**6

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = 4.0_wp*eps*sor6*(sor6-c)
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-c)*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

            Else If (ityp == 13) Then

! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}+c/r^12

              e0=prmvdw(1,k)
              r0=prmvdw(2,k)
              kk=prmvdw(3,k)
              c=prmvdw(4,k)

              t1=Exp(-kk*(rrr-r0))
              sor6 = c*r_rrr**12
              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*t1*(t1-2.0_wp)+sor6
              gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)*r_rrr-12.0_wp*sor6*r_rrr

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

            Else If (ityp == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)

              a = prmvdw(1,k)
              b = prmvdw(2,k)
              c = prmvdw(3,k)

              kk = rrr/c
              t1 = Exp(-kk)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = (a+b*rrr)*t1
              gamma = kk*t1*(a-b*c+b*rrr)*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

            Else If (ityp == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

              z1 = prmvdw(1,k)
              z2 = prmvdw(2,k)
        
        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0

              Call zbl(rrr,kk,a,t1,gamma)
              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng = t1
              gamma = gamma*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

            Else If (ityp == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

              z1 = prmvdw(1,k)
              z2 = prmvdw(2,k)
              rm = prmvdw(3,k)
              c = 1.0_wp/prmvdw(4,k)
              e0 = prmvdw(5,k)
              r0 = prmvdw(6,k)
              t2 = prmvdw(7,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0

              Call zbls(rrr,kk,a,rm,c,e0,t2,r0,t1,gamma)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng = t1
              gamma = gamma*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

            Else If (ityp == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

              z1 = prmvdw(1,k)
              z2 = prmvdw(2,k)
              rm = prmvdw(3,k)
              c = 1.0_wp/prmvdw(4,k)
              e0 = prmvdw(5,k)
              r0 = prmvdw(6,k)
              t2 = prmvdw(7,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0

              Call zblb(rrr,kk,a,rm,c,e0,r0,t2,t1,gamma)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng = t1
              gamma = gamma*r_rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)*r_rrr
              End If

           Else If (Abs(vvdw(0,k)) > zero_plus) Then ! potential read from TABLE - (ityp == 0)

              l   = Int(rrr*rdr)
              ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

              If (jatm <= natms .or. idi < ltg(jatm)) Then
                 vk  = vvdw(l,k)
                 vk1 = vvdw(l+1,k)
                 vk2 = vvdw(l+2,k)

                 t1 = vk  + (vk1 - vk )*ppp
                 t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                 eng = t1 + (t2-t1)*ppp*0.5_wp
                 If (ls_vdw) eng = eng + gvdw(mxgvdw-4,k)*(rscl-1.0_wp) - vvdw(mxgvdw-4,k) ! force-shifting
              End If

! calculate forces using 3-point interpolation

              gk  = gvdw(l,k) ; If (l == 0) gk = gk*rrr
              gk1 = gvdw(l+1,k)
              gk2 = gvdw(l+2,k)

              t1 = gk  + (gk1 - gk )*ppp
              t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

              gamma = (t1 + (t2-t1)*ppp*0.5_wp)*r_rsq
              If (ls_vdw) gamma = gamma - gvdw(mxgvdw-4,k)*r_rrv ! force-shifting

           End If

        Else If (Abs(vvdw(0,k)) > zero_plus) Then ! no direct = fully tabulated calculation

           l   = Int(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

           If (jatm <= natms .or. idi < ltg(jatm)) Then
              vk  = vvdw(l,k)
              vk1 = vvdw(l+1,k)
              vk2 = vvdw(l+2,k)

              t1 = vk  + (vk1 - vk )*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              eng = t1 + (t2-t1)*ppp*0.5_wp
              If (ls_vdw) eng = eng + gvdw(mxgvdw-4,k)*(rscl-1.0_wp) - vvdw(mxgvdw-4,k) ! force-shifting
           End If

! calculate forces using 3-point interpolation

           gk  = gvdw(l,k) ; If (l == 0) gk = gk*rrr
           gk1 = gvdw(l+1,k)
           gk2 = gvdw(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)*r_rsq
           If (ls_vdw) gamma = gamma - gvdw(mxgvdw-4,k)*r_rrv ! force-shifting

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


End Module vdw
