Module mpoles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring configuration variables and arrays for
! multipoles
!
! copyright - daresbury laboratory
! author    - h.a.boateng & i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,  Only : comms_type
  Use setup_module, Only : mxsite,mxexcl,mxspl,mxompl,mximpl,mxatdm,mxatms, &
                           mxlist, r4pie0, zero_plus, mxgele, nrite, sqrpi
  Use configuration, Only : imcon,natms,ltg,fxx,fyy,fzz,xxx,yyy,zzz,cell,list

  Implicit None

  Private

! A derived data type for defining the rotation matrix for multipolar sites

  Type rot_mat
     Real( Kind = wp ), Dimension(1:9) :: mtrxa
     Real( Kind = wp ), Dimension(1:3) :: p1,p2
     Integer                           :: flag,mbnd(1:2)
  End Type rot_mat

! Type of inducible (self-polarisation) scheme

  Integer,           Save :: keyind = 0 ! 0 - default :: unscreened & undamped - iAMOEBA like
                                        ! 1 - CHARMM  :: q_shell == -Sign(q_core) * Sqrt(alpha * k) ; 
                                        !                k_CHARMM = 1000 kcal*mol^−1*Å^−2

  Real( Kind = wp ), Save :: thole  = 1.3_wp  ! default thole dumping for CHARMM representation

! variables for multipolar interactions

  Logical,           Save :: induce=.false. , &
                             gear,aspc,lstsq
  Integer,           Save :: numcof,politer
  Real( Kind = wp ), Save :: convcrit,enepol

  Integer,           Allocatable, Save :: mplmap(:,:,:),mplltg(:) ! mappings from three indices multipole to a one index multipole
  Integer,           Allocatable, Save :: mplflg(:)               ! rotation counter flag
  Integer,           Allocatable, Save :: ltpatm(:,:)             ! bonded connectivity
  Integer,           Allocatable, Save :: lchatm(:,:)             ! CHARMM core-shell screened electrostatics induction list

  Real( Kind = wp ), Allocatable, Save :: mpllfr(:,:),mplgfr(:,:) ! local/lab(site) and global(atom) frames
  Real( Kind = wp ), Allocatable, Save :: plrsit(:),plratm(:)     ! induced dipole polarisation for sites 
                                                                  ! and atoms (inversed if non-zero)
  Real( Kind = wp ), Allocatable, Save :: dmpsit(:),dmpatm(:)     ! sites' and atoms' (thole) dumping 
                                                                  ! coefficient/factor (for self-polarisation)

  Type( rot_mat ),   Allocatable, Save :: mprotm(:)                           ! rotation matrices

  Real( Kind = wp ), Allocatable, Save :: mprotx(:,:),mproty(:,:),mprotz(:,:) ! infinitesimal rotations
  Real( Kind = wp ), Allocatable, Save :: mptrqx(:),mptrqy(:),mptrqz(:)       ! torques due to infinitesimal rotations

  Real( Kind = wp ), Allocatable, Save :: ncombk(:,:)                         ! n combination k values for usage 
                                                                              ! in computing the reciprocal space Ewald sum

  Real( Kind = wp ), Allocatable, Save :: mpfldx(:),mpfldy(:),mpfldz(:)       ! field

  Real( Kind = wp ), Allocatable, Save :: muindx(:),muindy(:),muindz(:)
  Real( Kind = wp ), Allocatable, Save :: indipx(:),indipy(:),indipz(:)

  Real( Kind = wp ), Allocatable, Save :: upidpx(:,:),upidpy(:,:),upidpz(:,:)
  Real( Kind = wp ), Allocatable, Save :: rsdx(:),rsdy(:),rsdz(:)
  Real( Kind = wp ), Allocatable, Save :: polcof(:)


  Public :: keyind, lchatm, ltpatm, mplmap, mpllfr, mprotm, mplltg, mplflg, &
            mplgfr, mprotx, mproty, mprotz, plratm, dmpatm
  Public :: allocate_mpoles_arrays

Contains

  Subroutine allocate_mpoles_arrays()

    Integer           :: n,k,om1,numpl,fail(1:9)
    Real( Kind = wp ) :: Factorial,gearp(1:7),aspcp(1:7)

    If (mximpl < 1) Return ! no MPOLES file read <= no multipoles directive in FIELD

    om1 = mxompl + 1
    numpl = (3**om1 - 1)/2

    fail = 0

    Allocate (mplmap(0:mxompl,0:mxompl,0:mxompl),mplltg(1:numpl),  Stat = fail(1))
    Allocate (mplflg(1:mxatdm),ltpatm(0:mxexcl,1:mxatdm),          Stat = fail(2))
    If (keyind == 1) &
    Allocate (lchatm(0:mxexcl,1:mxatdm),                           Stat = fail(3))
    Allocate (mpllfr(1:mximpl,1:mxsite),mplgfr(1:mximpl,1:mxatms), Stat = fail(4))
    Allocate (plrsit(1:mxsite),plratm(1:mxatms),                   Stat = fail(5))
    Allocate (dmpsit(1:mxsite),dmpatm(1:mxatms),                   Stat = fail(6))
    Allocate (mprotm(1:mxatdm),ncombk(0:mxspl,0:mxspl),            Stat = fail(7))
    Allocate (mptrqx(1:mxatdm),mptrqy(1:mxatdm),mptrqz(1:mxatdm),  Stat = fail(8))
    Allocate (mprotx(1:mximpl,1:mxatms), &
              mproty(1:mximpl,1:mxatms), &
              mprotz(1:mximpl,1:mxatms),                           Stat = fail(9))

    If (Any(fail > 0)) Call error(1025)

    mplflg = 0 ; ltpatm = 0
    If (keyind == 1) lchatm = 0

    mpllfr = 0.0_wp ; mplgfr = 0.0_wp
    plrsit = 0.0_wp ; plratm = 0.0_wp
    dmpsit = 0.0_wp ; dmpatm = 0.0_wp

    Do n=1,mxatdm
       mprotm(n)%mtrxa = 0.0_wp
       mprotm(n)%p1    = 0.0_wp
       mprotm(n)%p2    = 0.0_wp
       mprotm(n)%flag  = 0
       mprotm(n)%mbnd  = 0
    End Do

    mptrqx = 0.0_wp ; mptrqy = 0.0_wp ; mptrqz = 0.0_wp

    mprotx = 0.0_wp ; mproty = 0.0_wp ; mprotz = 0.0_wp

! Build the multipole map (polymap) and compute the constants ncombk
! Also build the map (mplltg) that converts between index of a local
! multipole to the index of the corresponding global multipole

    mplmap(0,0,0)=1

    mplltg(1)=1

    If (mxompl >= 1) Then

       mplmap(1,0,0)=2 ; mplmap(0,1,0)=3 ; mplmap(0,0,1)=4

       mplltg(2)=2 ; mplltg(3)=3 ; mplltg(4)=4

    End If

    If (mxompl >= 2) Then

       mplmap(2,0,0)=5 ; mplmap(1,1,0)=6 ; mplmap(1,0,1)=7
       mplmap(0,2,0)=8 ; mplmap(0,1,1)=9 ; mplmap(0,0,2)=10

       mplltg(5) =5 ; mplltg(6) =6 ; mplltg(7) =7
       mplltg(8) =6 ; mplltg(9) =8 ; mplltg(10)=9
       mplltg(11)=7 ; mplltg(12)=9 ; mplltg(13)=10

    End If

    If (mxompl >= 3) Then

       mplmap(3,0,0)=11 ; mplmap(2,1,0)=12 ; mplmap(2,0,1)=13
       mplmap(1,2,0)=14 ; mplmap(1,1,1)=15 ; mplmap(1,0,2)=16
       mplmap(0,3,0)=17 ; mplmap(0,2,1)=18 ; mplmap(0,1,2)=19
       mplmap(0,0,3)=20

       mplltg(14)=11 ; mplltg(15)=12 ; mplltg(16)=13
       mplltg(17)=12 ; mplltg(18)=14 ; mplltg(19)=15
       mplltg(20)=13 ; mplltg(21)=15 ; mplltg(22)=16
       mplltg(23)=12 ; mplltg(24)=14 ; mplltg(25)=15
       mplltg(26)=14 ; mplltg(27)=17 ; mplltg(28)=18
       mplltg(29)=15 ; mplltg(30)=18 ; mplltg(31)=19
       mplltg(32)=13 ; mplltg(33)=15 ; mplltg(34)=16
       mplltg(35)=15 ; mplltg(36)=18 ; mplltg(37)=19
       mplltg(38)=16 ; mplltg(39)=19 ; mplltg(40)=20

    End If

    If (mxompl >=4) Then

       mplmap(4,0,0)=21 ; mplmap(3,1,0)=22 ; mplmap(3,0,1)=23
       mplmap(2,2,0)=24 ; mplmap(2,1,1)=25 ; mplmap(2,0,2)=26
       mplmap(1,3,0)=27 ; mplmap(1,2,1)=28 ; mplmap(1,1,2)=29
       mplmap(1,0,3)=30 ; mplmap(0,4,0)=31 ; mplmap(0,3,1)=32
       mplmap(0,2,2)=33 ; mplmap(0,1,3)=34 ; mplmap(0,0,4)=35

       mplltg(41)=21 ; mplltg(42)=22 ; mplltg(43)=23
       mplltg(44)=22 ; mplltg(45)=24 ; mplltg(46)=25
       mplltg(47)=23 ; mplltg(48)=25 ; mplltg(49)=26
       mplltg(50)=22 ; mplltg(51)=24 ; mplltg(52)=25
       mplltg(53)=24 ; mplltg(54)=27 ; mplltg(55)=28
       mplltg(56)=25 ; mplltg(57)=28 ; mplltg(58)=29
       mplltg(59)=23 ; mplltg(60)=25 ; mplltg(61)=26
       mplltg(62)=25 ; mplltg(63)=28 ; mplltg(64)=29
       mplltg(65)=26 ; mplltg(66)=29 ; mplltg(67)=30
       mplltg(68)=22 ; mplltg(69)=24 ; mplltg(70)=25
       mplltg(71)=24 ; mplltg(72)=27 ; mplltg(73)=28
       mplltg(74)=25 ; mplltg(75)=28 ; mplltg(76)=29
       mplltg(77)=24 ; mplltg(78)=27 ; mplltg(79)=28
       mplltg(80)=27 ; mplltg(81)=31 ; mplltg(82)=32
       mplltg(83)=28 ; mplltg(84)=32 ; mplltg(85)=33
       mplltg(86)=25 ; mplltg(87)=28 ; mplltg(88)=29
       mplltg(89)=28 ; mplltg(90)=32 ; mplltg(91)=33
       mplltg(92)=29 ; mplltg(93)=33 ; mplltg(94)=34
       mplltg(95)=23 ; mplltg(96)=25 ; mplltg(97)=26
       mplltg(98)=25 ; mplltg(99)=28 ; mplltg(100)=29

       mplltg(101)=26 ; mplltg(102)=29 ; mplltg(103)=30
       mplltg(104)=25 ; mplltg(105)=28 ; mplltg(106)=29
       mplltg(107)=28 ; mplltg(108)=32 ; mplltg(109)=33
       mplltg(110)=29 ; mplltg(111)=33 ; mplltg(112)=34
       mplltg(113)=26 ; mplltg(114)=29 ; mplltg(115)=30
       mplltg(116)=29 ; mplltg(117)=33 ; mplltg(118)=34
       mplltg(119)=30 ; mplltg(120)=34 ; mplltg(121)=35

    End If

    If (mxompl < 0 .or. mxompl >=5) Call error(2071)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute n choose k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Do k = 0, mxspl
       Do n = 0, mxspl
          ncombk(n,k) = Exp(Factorial(n)-Factorial(n-k)-Factorial(k))
       End Do
    End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    If (induce) Then

       Allocate (mpfldx(1:mxatms),mpfldy(1:mxatms),mpfldz(1:mxatms),            Stat = fail(1))
       Allocate (muindx(1:mxatms),muindy(1:mxatms),muindz(1:mxatms),            Stat = fail(2))
       Allocate (indipx(1:mxatms),indipy(1:mxatms),indipz(1:mxatms),            Stat = fail(3))
       Allocate (polcof(numcof),                                                Stat = fail(4))
       Allocate (upidpx(1:numcof,1:mxatms),upidpy(1:numcof,1:mxatms),upidpz(1:numcof,1:mxatms),&
                                                                                Stat = fail(5))
       Allocate (rsdx(1:mxatms),rsdy(1:mxatms),rsdz(1:mxatms),                  Stat = fail(6))

       If (Any(fail > 0)) Call error(1025)

       mpfldx = 0.0_wp; mpfldy = 0.0_wp; mpfldz = 0.0_wp

       muindx = 0.0_wp; muindy = 0.0_wp; muindz = 0.0_wp

       indipx = 0.0_wp; indipy = 0.0_wp; indipz = 0.0_wp

       upidpx = 0.0_wp; upidpy = 0.0_wp; upidpz = 0.0_wp

       rsdx   = 0.0_wp; rsdy   = 0.0_wp; rsdz   = 0.0_wp  ! arrays to store
                                                          ! residuals for
                                                          ! conjugate gradient
                                                          ! minimization of
                                                          ! polarization energy

! Coefficients for polynomial predictor

       polcof = 0.0_wp

! Gear predictor-corrector

       gearp = 0.0_wp

       gearp(1) =   6.0_wp ; gearp(2) = -15.0_wp ; gearp(3) =  20.0_wp
       gearp(4) = -15.0_wp ; gearp(5) =   6.0_wp ; gearp(6) =  -1.0_wp


! Always stable predictor-corrector (aspc)

       aspcp = 0.0_wp

       aspcp(1) =  22.0_wp/ 7.0_wp ; aspcp(2) = -55.0_wp/14.0_wp; aspcp(3) =  55.0_wp/21.0_wp
       aspcp(4) = -22.0_wp/21.0_wp ; aspcp(5) =   5.0_wp/21.0_wp; aspcp(6) =  -1.0_wp/42.0_wp

       If (gear) Then

          Do n = 1, numcof
             polcof(n) = gearp(n)
          End Do

       Else If (aspc) Then

          Do n = 1, numcof
             polcof(n) = aspcp(n)
          End Do

       End If

    End If

  End Subroutine allocate_mpoles_arrays

  Subroutine coul_fscp_mforces &
             (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles assuming a force shifted
  ! coulomb potential kernel
  !
  ! U is proportional to ( 1/r + aa*r  + bb ) such that dU(rcut)/dr = 0
  ! therefore aa = 1/(rcut)**2 and U(rcut) = 0 therefore bb = -2/(rcut)
  !
  ! Note: FS potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use mpoles_container, Only : explicit_fscp_rfp_loops, explicit_ewald_real_loops

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( comms_type ),                       Intent( In    ) :: comm

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                               rdrewd = 0.0_wp , &
                               aa     = 0.0_wp , &
                               bb     = 0.0_wp

    Integer           :: fail,idi,jatm,k1,k2,k3,s1,s2,s3,m,n, &
                         k,ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,rrr,alphan,engmpl,fix,fiy,fiz,fx,fy,fz, &
                         strs1,strs2,strs3,strs5,strs6,strs9,        &
                         ppp,vk0,vk1,vk2,t1,t2,kx,ky,kz,             &
                         txyz,erfcr,tix,tiy,tiz,tjx,tjy,             &
                         tjz,tmp,tmpi,tmpj,sx,sy,sz

    Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

    If (newjob) Then
       newjob = .false.

       If (alpha > zero_plus) Then
          damp = .true.
       Else
          damp = .false.
       End If

       If (damp) Then

  ! interpolation interval

          drewd = rcut/Real(mxgele-4,wp)

  ! reciprocal of interpolation interval

          rdrewd = 1.0_wp/drewd

          fail=0
          Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'coul_fscp_mforces allocation failure, idnode: ', comm%idnode
             Call error(0)
          End If

  ! generate error function complement tables for ewald sum

          Call erfcgen(rcut,alpha,mxgele,erc,fer)

  ! set force and potential shifting parameters (screened terms)

          aa =   fer(mxgele-4)*rcut
          bb = -(erc(mxgele-4)+aa*rcut)

       Else

  ! set force and potential shifting parameters (screened terms)

          aa =  1.0_wp/rcut**2
          bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)

       End If
    End If

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

  ! start of primary loop for forces evaluation

       Do m=1,list(0,iatm)

  ! atomic index

          jatm=list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! interatomic distance

          rrr = rrt(m)

  ! truncation of potential

          If (Maxval(Abs(jmp)) > zero_plus .and. rrr < rcut) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             If (damp) Then

  ! compute derivatives of 'r'

                Call coul_deriv(-1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r'

                a1 = aa*a1

  ! get the value of the ewald real space kernel using 3pt interpolation

                k   = Int(rrr*rdrewd)
                ppp = rrr*rdrewd - Real(k,wp)

                vk0 = erc(k)
                vk1 = erc(k+1)
                vk2 = erc(k+2)

                t1 = vk0 + (vk1 - vk0)*ppp
                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/alpha

  ! compute derivatives of the ewald real space kernel

                Call ewald_deriv(-2,2*mxompl+1,1,erfcr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

  ! scale the derivatives into the right form

                d1 = 2.0_wp*alpha*d1/sqrpi

  ! calculate forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then
                              Call explicit_fscp_rfp_loops &
                               (2*mxompl+1, k1,k2,k3, alpha, d1,a1,               &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n      = ks1+ks2+ks3
                                        alphan = alpha**n

                                        ii     = mplmap(s1,s2,s3)

                                        tmp    = alphan*d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                                        tmpi   = txyz       * tmp
                                        tmpj   = sx*imp(ii) * tmp

                                        t2     = txyz*imp(ii)
                                        t1     = alphan*t2

  ! energy

                                        engmpl = engmpl + t1*d1(ks1,ks2,ks3) + t2*a1(ks1,ks2,ks3)

  ! force
                                        t1     = t1*alpha

                                        fx     = fx     - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

                                        sx = -sx
                                     End Do

                                     sy = -sy
                                  End Do

                                  sz = -sz
                               End Do

                            End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                End If

             Else

  ! compute derivatives of '1/r'

                Call coul_deriv( 1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

  ! compute derivatives of 'r'

                Call coul_deriv(-1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r' and add to d1

                d1 = d1 + a1/rcut**2

  ! calculate potential forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus)  Then
                              Call explicit_ewald_real_loops &
                               (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj=mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n    = ks1+ks2+ks3

                                        ii   = mplmap(s1,s2,s3)

                                        tmp  = d1(ks1,ks2,ks3)

                                        tmpi = txyz       * tmp
                                        tmpj = sx*imp(ii) * tmp

                                        t1   = txyz*imp(ii)

  ! energy

                                        engmpl = engmpl + t1*tmp

  ! force

                                        fx     = fx     - t1*d1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

                                        sx = -sx
                                     End Do

                                     sy = -sy
                                  End Do

                                  sz = -sz
                               End Do

                            End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                End If

             End If

  ! shift potential

             tmp     = aa*rrr + bb
             engmpl  = engmpl + tmp*imp(1)*jmp(1)

  ! shift torque

             tmpi    = tmp*jmp(1)
             tix     = tix    + impx(1)*tmpi
             tiy     = tiy    + impy(1)*tmpi
             tiz     = tiz    + impz(1)*tmpi

             tmpj    = tmp*imp(1)
             tjx     = tjx    + jmpx(1)*tmpj
             tjy     = tjy    + jmpy(1)*tmpj
             tjz     = tjz    + jmpz(1)*tmpj

             fix=fix+fx
             fiy=fiy+fy
             fiz=fiz+fz

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

  ! calculate virial

                vircpe = vircpe - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

  ! calculate stress tensor

                strs1 = strs1 + xxt(m)*fx
                strs2 = strs2 + xxt(m)*fy
                strs3 = strs3 + xxt(m)*fz
                strs5 = strs5 + yyt(m)*fy
                strs6 = strs6 + yyt(m)*fz
                strs9 = strs9 + zzt(m)*fz

             End If

          End If

       End Do

  ! load back forces

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

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

    End If

  End Subroutine coul_fscp_mforces

  Subroutine coul_rfp_mforces &
             (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles assuming a reaction field
  ! potential (corrected for the existence of a dipole moment outside rcut)
  !
  ! Note: RF potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  ! R2: M Neumann, J Chem Phys, 82 (12), 5663, (1985)
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use mpoles_container, Only : explicit_fscp_rfp_loops, explicit_ewald_real_loops

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( comms_type ),                       Intent( In    ) :: comm

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                               rdrewd = 0.0_wp , &
                               aa     = 0.0_wp , &
                               bb     = 0.0_wp , &
                               b0     = 0.0_wp , &
                               rfld0  = 0.0_wp , &
                               rfld1  = 0.0_wp , &
                               rfld2  = 0.0_wp

    Integer           :: fail,idi,jatm,k1,k2,k3,s1,s2,s3,m,n, &
                         k,ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,rrr,alphan,engmpl,fix,fiy,fiz,fx,fy,fz, &
                         strs1,strs2,strs3,strs5,strs6,strs9,        &
                         ppp,vk0,vk1,vk2,t1,t2,kx,ky,kz,             &
                         txyz,erfcr,tmp,tmpi,tmpj,tix,tiy,tiz,       &
                         tjx,tjy,tjz,sx,sy,sz

    Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

    If (newjob) Then
       newjob = .false.

       If (alpha > zero_plus) Then
          damp = .true.
       Else
          damp = .false.
       End If

  ! reaction field terms

       b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
       rfld0 = b0/rcut**3
       rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
       rfld2 = 0.5_wp*rfld0

       If (damp) Then

  ! interpolation interval

          drewd = rcut/Real(mxgele-4,wp)

  ! reciprocal of interpolation interval

          rdrewd = 1.0_wp/drewd

          fail=0
          Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'coul_rfp_mforces allocation failure, idnode: ', comm%idnode
             Call error(0)
          End If

  ! generate error function complement tables for ewald sum

          Call erfcgen(rcut,alpha,mxgele,erc,fer)

  ! set force and potential shifting parameters (screened terms)

          aa =   fer(mxgele-4)*rcut
          bb = -(erc(mxgele-4)+aa*rcut)

       Else

  ! set force and potential shifting parameters (screened terms)

          aa =  1.0_wp/rcut**2
          bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)

       End If
    End If

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

       Do m=1,list(0,iatm)

  ! atomic index

          jatm=list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! interatomic distance

          rrr = rrt(m)

  ! truncation of potential

          If (Maxval(Abs(jmp)) > zero_plus .and. rrr < rcut) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             If (damp) Then

  ! compute derivatives of 'r^2'

                Call coul_deriv(-2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r^2'

                a1 = rfld2*a1

  ! compute derivatives of 'r'

                Call coul_deriv(-1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

  ! scale the derivatives of 'r' and add to a1

                a1 = a1 + aa*d1

  ! get the value of the ewald real space kernel using 3pt interpolation

                k   = Int(rrr*rdrewd)
                ppp = rrr*rdrewd - Real(k,wp)

                vk0 = erc(k)
                vk1 = erc(k+1)
                vk2 = erc(k+2)

                t1 = vk0 + (vk1 - vk0)*ppp
                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/alpha

  ! compute derivatives of the ewald real space kernel

                Call ewald_deriv(-2,2*mxompl+1,1,erfcr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

  ! scale the derivatives into the right form

                d1 = 2.0_wp*alpha*d1/sqrpi

  ! calculate forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then
                              Call explicit_fscp_rfp_loops &
                               (2*mxompl+1, k1,k2,k3, alpha, d1,a1,               &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n      = ks1+ks2+ks3
                                        alphan = alpha**n

                                        ii     = mplmap(s1,s2,s3)

                                        tmp    = alphan*d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                                        tmpi   = txyz       * tmp
                                        tmpj   = sx*imp(ii) * tmp

                                        t2     = txyz*imp(ii)
                                        t1     = alphan*t2

  ! energy

                                        engmpl = engmpl + t1*d1(ks1,ks2,ks3) + t2*a1(ks1,ks2,ks3)

  ! force
                                        t1     = t1*alpha

                                        fx     = fx     - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

                                        sx = -sx
                                     End Do

                                     sy = -sy
                                  End Do

                                  sz = -sz
                               End Do

                            End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                End If

  ! shift potential

                tmp    = bb-0.5_wp*b0/rcut
                engmpl = engmpl + tmp*imp(1)*jmp(1)

  ! shift torque

                tmpi   = tmp*jmp(1)
                tix    = tix    + impx(1)*tmpi
                tiy    = tiy    + impy(1)*tmpi
                tiz    = tiz    + impz(1)*tmpi

                tmpj   = tmp*imp(1)
                tjx    = tjx    + jmpx(1)*tmpj
                tjy    = tjy    + jmpy(1)*tmpj
                tjz    = tjz    + jmpz(1)*tmpj

             Else

  ! compute derivatives of '1/r'

                Call coul_deriv( 1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

  ! compute derivatives of 'r^2'

                Call coul_deriv(-2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r' and add to d1

                d1 = d1 + rfld2*a1

  ! calculate potential forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then
                              Call explicit_ewald_real_loops &
                               (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n    = ks1+ks2+ks3

                                        ii   = mplmap(s1,s2,s3)

                                        tmp  = d1(ks1,ks2,ks3)

                                        tmpi = txyz       * tmp
                                        tmpj = sx*imp(ii) * tmp

                                        t1   = txyz*imp(ii)

  ! energy

                                        engmpl = engmpl + t1*tmp

  ! force

                                        fx     = fx     - t1*d1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

                                        sx = -sx
                                     End Do

                                     sy = -sy
                                  End Do

                                  sz = -sz
                               End Do

                            End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                End If

  ! shift potential

                engmpl = engmpl - rfld1*imp(1)*jmp(1)

  ! shift torque

                tmpi   = -rfld1*jmp(1)
                tix    = tix    + impx(1)*tmpi
                tiy    = tiy    + impy(1)*tmpi
                tiz    = tiz    + impz(1)*tmpi

                tmpj   = -rfld1*imp(1)
                tjx    = tjx    + jmpx(1)*tmpj
                tjy    = tjy    + jmpy(1)*tmpj
                tjz    = tjz    + jmpz(1)*tmpj

             End If

             fix=fix+fx
             fiy=fiy+fy
             fiz=fiz+fz

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

  ! calculate virial

                vircpe = vircpe - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

  ! calculate stress tensor

                strs1 = strs1 + xxt(m)*fx
                strs2 = strs2 + xxt(m)*fy
                strs3 = strs3 + xxt(m)*fz
                strs5 = strs5 + yyt(m)*fy
                strs6 = strs6 + yyt(m)*fz
                strs9 = strs9 + zzt(m)*fz

             End If

          End If

       End Do

  ! load back forces

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

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

    End If

  End Subroutine coul_rfp_mforces

  Subroutine coul_cp_mforces &
             (iatm,rcut,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles with 1/r kernel with no
  ! truncation or damping
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use mpoles_container, Only : explicit_ewald_real_loops

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: idi,jatm,k1,k2,k3,s1,s2,s3,m, &
                         ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,engmpl,fix,fiy,fiz,fx,fy,fz,     &
                         strs1,strs2,strs3,strs5,strs6,strs9, &
                         t1,kx,ky,kz,txyz,tix,tiy,tiz,tjx,    &
                         tjy,tjz,tmp,tmpi,tmpj,sx,sy,sz

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

  ! start of primary loop for forces evaluation

       Do m=1,list(0,iatm)

  ! atomic index

          jatm=list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! truncation of potential - rrt(m) is the interatomic distance

          If (Maxval(Abs(jmp)) > zero_plus .and. rrt(m) < rcut) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             Call coul_deriv(1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrt(m),d1)

  ! calculate forces

             engmpl = 0.0_wp
             fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
             tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

             If (mxompl < 5) Then

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then
                           Call explicit_ewald_real_loops &
                             (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                             imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                             kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                             engmpl,fx,fy,fz)
                         End If

                         kx = -kx

                      End Do

                      ky = -ky

                   End Do

                   kz = -kz

                End Do

             Else

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then

                            txyz=kx*jmp(jj)

                            sz = 1.0_wp
                            Do s3=0,mxompl
                               ks3=k3+s3; ks31=ks3+1

                               sy = sz
                               Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                  sx = sy
                                  Do s1=0,mxompl-s3-s2
                                     ks1=k1+s1; ks11=ks1+1

                                     ii   = mplmap(s1,s2,s3)

                                     tmp  = d1(ks1,ks2,ks3)

                                     tmpi = txyz       * tmp
                                     tmpj = sx*imp(ii) * tmp

                                     t1   = txyz*imp(ii)

  ! energy

                                     engmpl = engmpl + t1*tmp

  ! force

                                     fx     = fx     - t1*d1(ks11,ks2,ks3)
                                     fy     = fy     - t1*d1(ks1,ks21,ks3)
                                     fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                     tix    = tix    + impx(ii)*tmpi
                                     tiy    = tiy    + impy(ii)*tmpi
                                     tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                     tjx    = tjx    + jmpx(jj)*tmpj
                                     tjy    = tjy    + jmpy(jj)*tmpj
                                     tjz    = tjz    + jmpz(jj)*tmpj

                                     sx = -sx
                                  End Do

                                  sy = -sy
                               End Do

                               sz = -sz
                            End Do

                         End If

                         kx = -kx

                      End Do

                      ky = -ky

                   End Do

                   kz = -kz

                End Do

             End If

             fix=fix+fx
             fiy=fiy+fy
             fiz=fiz+fz

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

  ! calculate stress tensor

                strs1 = strs1 + xxt(m)*fx
                strs2 = strs2 + xxt(m)*fy
                strs3 = strs3 + xxt(m)*fz
                strs5 = strs5 + yyt(m)*fy
                strs6 = strs6 + yyt(m)*fz
                strs9 = strs9 + zzt(m)*fz

             End If

          End If

       End Do

  ! load back forces

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

  ! virial

       vircpe = -engcpe

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

    End If

  End Subroutine coul_cp_mforces

  Subroutine coul_dddp_mforces &
             (iatm,rcut,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles with 1/r kernel assuming a
  ! distance dependent dielectric 'constant'
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Use mpoles_container, Only : explicit_ewald_real_loops

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: idi,jatm,k1,k2,k3,s1,s2,s3,m, &
                         ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,engmpl,fix,fiy,fiz,fx,fy,fz,     &
                         strs1,strs2,strs3,strs5,strs6,strs9, &
                         t1,kx,ky,kz,txyz,tix,tiy,tiz,tjx,    &
                         tjy,tjz,tmp,tmpi,tmpj,sx,sy,sz

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

  ! start of primary loop for forces evaluation

       Do m=1,list(0,iatm)

  ! atomic index

          jatm=list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! truncation of potential - rrt(m) is the interatomic distance

          If (Maxval(Abs(jmp)) > zero_plus .and. rrt(m) < rcut) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             Call coul_deriv(2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrt(m),d1)

  ! calculate forces

             engmpl = 0.0_wp
             fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
             tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

             If (mxompl < 5) Then

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then
                           Call explicit_ewald_real_loops &
                             (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                             imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                             kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                             engmpl,fx,fy,fz)
                         End If

                         kx = -kx

                      End Do

                      ky = -ky

                   End Do

                   kz = -kz

                End Do

             Else

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then

                            txyz=kx*jmp(jj)

                            sz = 1.0_wp
                            Do s3=0,mxompl
                               ks3=k3+s3; ks31=ks3+1

                               sy = sz
                               Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                  sx = sy
                                  Do s1=0,mxompl-s3-s2
                                     ks1=k1+s1; ks11=ks1+1

                                     ii   = mplmap(s1,s2,s3)

                                     tmp  = d1(ks1,ks2,ks3)

                                     tmpi = txyz       * tmp
                                     tmpj = sx*imp(ii) * tmp

                                     t1   = txyz*imp(ii)

  ! energy

                                     engmpl = engmpl  + t1*tmp

  ! force

                                     fx     = fx     - t1*d1(ks11,ks2,ks3)
                                     fy     = fy     - t1*d1(ks1,ks21,ks3)
                                     fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                     tix    = tix    + impx(ii)*tmpi
                                     tiy    = tiy    + impy(ii)*tmpi
                                     tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                     tjx    = tjx    + jmpx(jj)*tmpj
                                     tjy    = tjy    + jmpy(jj)*tmpj
                                     tjz    = tjz    + jmpz(jj)*tmpj

                                     sx = -sx
                                  End Do

                                  sy = -sy
                               End Do

                               sz = -sz
                            End Do

                         End If

                         kx = -kx

                      End Do

                      ky = -ky

                   End Do

                   kz = -kz

                End Do

             End If

             fix=fix+fx
             fiy=fiy+fy
             fiz=fiz+fz

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

  ! calculate stress tensor

                strs1 = strs1 + xxt(m)*fx
                strs2 = strs2 + xxt(m)*fy
                strs3 = strs3 + xxt(m)*fz
                strs5 = strs5 + yyt(m)*fy
                strs6 = strs6 + yyt(m)*fz
                strs9 = strs9 + zzt(m)*fz

             End If

          End If

       End Do

  ! load back forces

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

  ! virial

       vircpe = -engcpe

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

    End If

  End Subroutine coul_dddp_mforces

  Subroutine d_ene_trq_mpoles(vircpe_dt,stress)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating the change in energy produced by
  ! an infinitesimal rotation of multipoles
  !
  ! Reference : Sagui, Pedersen, Darden, J. Chem. Phys. 120, 73 (2004)
  !             doi: 10.1063/1.1630791
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng april 2015
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Real( Kind = wp ),                   Intent(   Out ) :: vircpe_dt
    Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

    Integer           :: idi,j,iatm,jatm

    Real( Kind = wp ) :: rsq,rrr,p1(1:3),p2(1:3),u(1:3),v(1:3),w(1:3),magp1,magp2, &
                         p2u(1:3),p1perp2(1:3),p2perp1(1:3),magp1perp2,magp2perp1, &
                         tx,ty,tz,dedv,dedw,stx(1:2),sty(1:2),stz(1:2),            &
                         dux_dx,dux_dy,dux_dz,duy_dy,duy_dz,duz_dz,rmag,rmag3,     &
                         fx,fy,fz,fix,fiy,fiz,dedux,deduy,deduz,                   &
                         tmptx,tmpty,tmptz,xdf,ydf,zdf,                            &
                         strs1,strs2,strs3,strs4,strs5,strs6,strs7,strs8,strs9

  ! initialise virial

    vircpe_dt=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs4=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs7=0.0_wp
    strs8=0.0_wp
    strs9=0.0_wp

    Do iatm = 1, natms

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! global identity of iatm

       idi=ltg(iatm)

       If (mprotm(iatm)%flag == 1) Then

  ! p1 and p2 define the local frame

          p1 = mprotm(iatm)%p1
          p2 = mprotm(iatm)%p2

          magp1 = Sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
          magp2 = Sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))

          p2u = p2/magp2

  ! standard basis for coordinate system

          u(1)  = mprotm(iatm)%mtrxa(1)
          u(2)  = mprotm(iatm)%mtrxa(4)
          u(3)  = mprotm(iatm)%mtrxa(7)

          v(1)  = mprotm(iatm)%mtrxa(2)
          v(2)  = mprotm(iatm)%mtrxa(5)
          v(3)  = mprotm(iatm)%mtrxa(8)

          w(1)  = mprotm(iatm)%mtrxa(3)
          w(2)  = mprotm(iatm)%mtrxa(6)
          w(3)  = mprotm(iatm)%mtrxa(9)

  ! change in energy (E) due to infinitesimal rotation (torque => \tau) dE_{\omega} = -\tau * d\omega
  ! we omit the negative here and introduce it in the final force computation, i.e., because force is
  ! the negative of the change in energy with respect to energy, we'll add the magnitude of the force
  ! computed instead of subtracting the magnitude

          tmptx = mptrqx(iatm)
          tmpty = mptrqy(iatm)
          tmptz = mptrqz(iatm)

          tx = tmptx*u(1) + tmpty*u(2) + tmptz*u(3)
          ty = (tmptx*p2(1) + tmpty*p2(2) + tmptz*p2(3)) / magp2
          tz = tmptx*w(1) + tmpty*w(2) + tmptz*w(3)

  ! component of p2 perpendicular to p1

          p2perp1    = p2 - (u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3))*u
          magp2perp1 = Sqrt(p2perp1(1)*p2perp1(1)+p2perp1(2)*p2perp1(2)+p2perp1(3)*p2perp1(3))

  ! component of p1 perpendicular to p2

          p1perp2    = p1 - (p2u(1)*p1(1) + p2u(2)*p1(2) + p2u(3)*p1(3))*p2u
          magp1perp2 = Sqrt(p1perp2(1)*p1perp2(1)+p1perp2(2)*p1perp2(2)+p1perp2(3)*p1perp2(3))

  ! For p1

  ! dedu = 0.0 since movement du along the u axis does not rotate the frame

          dedv   = tz/magp1
          dedw   =-ty/magp1perp2

          stx(1) = dedv*v(1) + dedw*w(1)
          sty(1) = dedv*v(2) + dedw*w(2)
          stz(1) = dedv*v(3) + dedw*w(3)

  ! for p2

          dedw   = tx/magp2perp1

          stx(2) = dedw*w(1)
          sty(2) = dedw*w(2)
          stz(2) = dedw*w(3)

  ! now compute forces and virial

          fx = 0.0_wp ; fy = 0.0_wp ; fz = 0.0_wp

          Do j = 1, 2

             jatm = mprotm(iatm)%mbnd(j)

             If (jatm > 0) Then
                xdf = xxx(jatm) - xxx(iatm)
                ydf = yyy(jatm) - yyy(iatm)
                zdf = zzz(jatm) - zzz(iatm)

                Call images_s(imcon,cell,xdf,ydf,zdf)

                rsq   = xdf*xdf + ydf*ydf + zdf*zdf
                rrr   = sqrt(rsq)
                rmag  = 1.0_wp/rrr
                rmag3 = 1.0_wp/(rsq*rrr)

                dux_dx = rmag - xdf * xdf * rmag3
                dux_dy =      - xdf * ydf * rmag3       ! duy_dx = dux_dy
                dux_dz =      - xdf * zdf * rmag3       ! duz_dx = dux_dz
                duy_dy = rmag - ydf * ydf * rmag3
                duy_dz =      - ydf * zdf * rmag3       ! duz_dy = duy_dz
                duz_dz = rmag - zdf * zdf * rmag3

                dedux  = stx(j)
                deduy  = sty(j)
                deduz  = stz(j)

  ! Now to find the forces (derivatives of energy with respect to cartesian positions),
  ! i.e. fx=dedx, fy=dedy, fz=dedz

                fx = dedux * dux_dx + deduy * dux_dy + deduz * dux_dz
                fy = dedux * dux_dy + deduy * duy_dy + deduz * duy_dz
                fz = dedux * dux_dz + deduy * duy_dz + deduz * duz_dz

                fix = fix + fx
                fiy = fiy + fy
                fiz = fiz + fz

                If (jatm <= natms) Then

                   fxx(jatm)=fxx(jatm)-fx
                   fyy(jatm)=fyy(jatm)-fy
                   fzz(jatm)=fzz(jatm)-fz

                End If

                If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! calculate virial

                   vircpe_dt = vircpe_dt - (fx*xdf + fy*ydf + fz*zdf)

  ! calculate stress tensor

                   strs1 = strs1 + xdf*fx
                   strs2 = strs2 + xdf*fy
                   strs3 = strs3 + xdf*fz
                   strs4 = strs4 + ydf*fx
                   strs5 = strs5 + ydf*fy
                   strs6 = strs6 + ydf*fz
                   strs7 = strs7 + zdf*fx
                   strs8 = strs8 + zdf*fy
                   strs9 = strs9 + zdf*fz

                End If

             End If

          End Do

  ! load back forces

          fxx(iatm)=fix
          fyy(iatm)=fiy
          fzz(iatm)=fiz

       End If

  ! complete stress tensor (and symmetrise)

       stress(1) = stress(1) + strs1
       stress(2) = stress(2) + 0.5_wp * (strs2 + strs4)
       stress(3) = stress(3) + 0.5_wp * (strs3 + strs7)
       stress(4) = stress(4) + 0.5_wp * (strs2 + strs4)
       stress(5) = stress(5) + strs5
       stress(6) = stress(6) + 0.5_wp * (strs6 + strs8)
       stress(7) = stress(7) + 0.5_wp * (strs3 + strs7)
       stress(8) = stress(8) + 0.5_wp * (strs6 + strs8)
       stress(9) = stress(9) + strs9

    End Do

  End Subroutine d_ene_trq_mpoles

End Module mpoles
