Module mpole
  Use kinds, Only : wp
  Use setup_module, Only : mxatdm, mxatms, mxexcl, mximpl, mxompl, mxsite, mxspl, &
                           sqrpi,r4pie0,zero_plus,nrite,nmpldt
  Use configuration,Only : natms
  Use site_module
  Use core_shell, Only : numshl,lstshl
  Use parse_module

  Implicit None

  !Type, Public mpole_type
  !
  !End Type mpole_type

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
  
  Subroutine intra_mcoul(keyfce,rcut,alpha,epsq,iatm,jatm,scale, &
                      rrr,xdf,ydf,zdf,coul,virele,fx,fy,fz,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bond's or 1-4 dihedral
! electrostatics: adjusted by a weighting factor
!
! copyright - daresbury laboratory
! amended   - i.t.todorov & h.a.boateng february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: keyfce,iatm,jatm
  Real( Kind = wp ), Intent( In    ) :: rcut,alpha,epsq,scale
  Real( Kind = wp ), Intent( In    ) :: xdf,ydf,zdf,rrr
  Real( Kind = wp ), Intent(   Out ) :: coul,virele,fx,fy,fz
  Logical,           Intent( InOut ) :: safe

  Integer                 :: k1,k2,k3,s1,s2,s3,n
  Integer                 :: ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

  Logical,           Save :: newjob = .true. , damp
  Real( Kind = wp ), Save :: aa     = 0.0_wp , &
                             bb     = 0.0_wp , &
                             rfld0  = 0.0_wp , &
                             rfld1  = 0.0_wp , &
                             rfld2  = 0.0_wp

  Real( Kind = wp ) :: exp1,tt,erc,fer,b0,      &
                       tix,tiy,tiz,tjx,tjy,tjz, &
                       talpha,alphan,tmp,tmpi,tmpj,t1,t2,sx,sy,sz,kx,ky,kz,txyz

  Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: b1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
  Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
  Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

  If (newjob) Then
     newjob = .false.

! Check for damped force-shifted coulombic and reaction field interactions
! and set force and potential shifting parameters dependingly

     damp=.false.
     If (alpha > zero_plus) Then
        damp=.true.

        exp1= Exp(-(alpha*rcut)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rcut)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rcut
        fer = (erc + 2.0_wp*(alpha/sqrpi)*exp1)/rcut**2

        aa  = fer*rcut
        bb  = -(erc + aa*rcut)
     Else If (keyfce == 8) Then
        aa =  1.0_wp/rcut**2
        bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)
     End If

! set reaction field terms for RFC

     If (keyfce == 10) Then
        b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
        rfld0 = b0/rcut**3
        rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
        rfld2 = 0.5_wp*rfld0
     End If
  End If

! initialise defaults for coulombic energy and force contributions

   coul =0.0_wp ; virele=0.0_wp
   fx =0.0_wp ; fy =0.0_wp ; fz =0.0_wp
   tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp
   tjx=0.0_wp ; tjy=0.0_wp ; tjz=0.0_wp

! get the multipoles for sites i and j

  imp=mplgfr(:,iatm)
  jmp=mplgfr(:,jatm)

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) <= zero_plus .or. Maxval(Abs(jmp)) <= zero_plus) Return

  If(mxompl > 0 .and. induce) Then

     imp(2)=imp(2)+indipx(iatm)
     imp(3)=imp(3)+indipy(iatm)
     imp(4)=imp(4)+indipz(iatm)

     jmp(2)=jmp(2)+indipx(jatm)
     jmp(3)=jmp(3)+indipy(jatm)
     jmp(4)=jmp(4)+indipz(jatm)

  End If

! get the components for site i and j infinitesimal rotations

  impx=mprotx(:,iatm)
  impy=mproty(:,iatm)
  impz=mprotz(:,iatm)

  jmpx=mprotx(:,jatm)
  jmpy=mproty(:,jatm)
  jmpz=mprotz(:,jatm)

! default convergence factor and derivative of 'r'

  talpha = 1.0_wp
  a1     = 0.0_wp

! scale imp multipoles

  imp = scale*imp*r4pie0/epsq

! Electrostatics by ewald sum = direct coulombic

  If      (keyfce ==  2 .or. keyfce ==  6) Then

! compute derivatives of 1/r kernel

     Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! distance dependent dielectric

  Else If (keyfce ==  4) Then

! Compute derivatives of 1/r^2 kernel

     Call coul_deriv(2,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! force shifted coulombic and reaction field

  Else If (keyfce ==  8 .or. keyfce == 10) Then

     If (damp) Then ! calculate damping contributions

! compute derivatives of 'r'

        Call coul_deriv(-1,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

! scale the derivatives of 'r'

        a1 = aa*a1

        exp1= Exp(-(alpha*rrr)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rrr)

        fer = erc/alpha

! compute derivatives of the ewald real space kernel

        Call ewald_deriv(-2,2*mxompl+1,1,fer,alpha*xdf,alpha*ydf,alpha*zdf,alpha*rrr,d1)

! scale the derivatives into the right form

        d1 = 2.0_wp*alpha*d1/sqrpi

     End If

     If      (keyfce ==  8) Then ! force shifted coulombic
        If (.not.damp) Then ! pure

! compute derivatives of '1/r'

           Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! compute derivatives of 'r'

           Call coul_deriv(-1,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

! scale the derivatives of 'r' and add to d1

           d1 = d1 + aa*a1
           a1 = 0.0_wp

        Else                ! damped

           talpha = alpha

        End If
     Else If (keyfce == 10) Then ! reaction field
        If (.not.damp) Then ! pure

! compute derivatives of '1/r'

           Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! compute derivatives of 'r^2'

           Call coul_deriv(-2,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

! scale the derivatives of 'r' and add to d1

           d1 = d1 + rfld2*a1
           a1 = 0.0_wp

        Else

! compute derivatives of 'r^2'

           Call coul_deriv(-2,2*mxompl+1,xdf,ydf,zdf,rrr,b1)

! scale the derivatives of 'r^2' and add to scaled derivatives of 'r'

           a1 = a1 + rfld2*b1

           talpha = alpha

        End If
     End If

  Else

     safe = .false.

  End If

  If (safe) Then

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
                          alphan = talpha**n

                          ii     = mplmap(s1,s2,s3)

                          tmp    = alphan * d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                          tmpi   = txyz       * tmp
                          tmpj   = sx*imp(ii) * tmp

                          t2     = txyz*imp(ii)
                          t1     = alphan*t2

!  energy

                          coul    = coul + t1*d1(ks1,ks2,ks3)  + t2*a1(ks1,ks2,ks3)

!  force
                          t1      = t1*talpha

                          fx      = fx   - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                          fy      = fy   - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                          fz      = fz   - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

!  torque on iatm

                          tix     = tix  + impx(ii)*tmpi
                          tiy     = tiy  + impy(ii)*tmpi
                          tiz     = tiz  + impz(ii)*tmpi

!  torque on jatm

                          tjx     = tjx  + jmpx(jj)*tmpj
                          tjy     = tjy  + jmpy(jj)*tmpj
                          tjz     = tjz  + jmpz(jj)*tmpj

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

     If (keyfce == 2 .or. keyfce == 4 .or. keyfce == 6) Then

        virele = -coul

     Else If (keyfce ==  8) Then ! force shifted coulombic

        virele = -(fx*xdf + fy*ydf + fz*zdf)

! shift potential

        tmp    = aa*rrr + bb
        coul   = coul   + tmp*imp(1)*jmp(1)

! shift torque

        tmpi   = tmp*jmp(1)
        tix    = tix    + impx(1)*tmpi
        tiy    = tiy    + impy(1)*tmpi
        tiz    = tiz    + impz(1)*tmpi

        tmpj   = tmp*imp(1)
        tjx    = tjx    + jmpx(1)*tmpj
        tjy    = tjy    + jmpy(1)*tmpj
        tjz    = tjz    + jmpz(1)*tmpj

     Else If (keyfce == 10) Then ! reaction field

        virele = -(fx*xdf + fy*ydf + fz*zdf)

! shift potential

        coul   = coul   - rfld1*imp(1)*jmp(1)

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

     tix = tix*r4pie0/epsq
     tiy = tiy*r4pie0/epsq
     tiz = tiz*r4pie0/epsq

     If (iatm <= natms) Then

        mptrqx(iatm)=mptrqx(iatm)+tix
        mptrqy(iatm)=mptrqy(iatm)+tiy
        mptrqz(iatm)=mptrqz(iatm)+tiz

     End If

     If (jatm <= natms) Then

        mptrqx(jatm)=mptrqx(jatm)+tjx
        mptrqy(jatm)=mptrqy(jatm)+tjy
        mptrqz(jatm)=mptrqz(jatm)+tjz

     End If

  End If

End Subroutine intra_mcoul

Subroutine read_mpoles(l_top,sumchg,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the molecular mulitpole
! specifications of the system to be simulated
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent ( In    ) :: l_top
  Real( Kind = wp ), Intent ( InOut ) :: sumchg
  Type( comms_type ), Intent ( InOut ) :: comm
  Logical                :: safe,l_rsh,l_ord=.false.

  Character( Len = 200 ) :: record,record1,record2
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: atom

  Integer                :: itmols,nrept,i,j,k,l,                 &
                            isite,jsite,ksite,lsite,nsite,sitmpl, &
                            ishls,nshels,kshels,isite2,           &
                            ordmpl,ordmpl_start,ordmpl_next,      &
                            ordmpl_min,ordmpl_max,                &
                            indmpl,indmpl_start,indmpl_final

  Real( Kind = wp )      :: Factorial,charge,scl,polarity,dumping

! open MPOLES data file

  If (comm%idnode == 0) Then
     Open(Unit=nmpldt, File = 'MPOLES', Status = 'old')
     Write(nrite,"(/,/,1x,'ELECTROSTATICS MULTIPOLES SPECIFICATION')")
     If (.not.l_top) Write(nrite,"(/,1x,'detailed specification opted out')")
  End If

  Call get_line(safe,nmpldt,record,comm)
  If (.not.safe) Go To 2000

! omit first line

  nsite  = 0
  nshels = 0
  sumchg = 0.0_wp

  ordmpl_min = 4
  ordmpl_max = 0

! read and process directives from mpols/field file

  Do

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nmpldt,record,comm)
        If (.not.safe) Go To 2000
        Call get_word(record,word) ; Call lower_case(word)
     End Do

! specify molecular species

     If (word(1:7) == 'molecul') Then

        Call get_word(record,word) ; Call lower_case(word)
        If (word(1:4) == 'type') Call get_word(record,word)

        If (ntpmls == Nint(word_2_real(word,comm))) Then
           If (comm%idnode == 0) Write(nrite,"(/,/,1x,'number of molecular types',6x,i10)") ntpmls
        Else
           If (comm%idnode == 0) Write(nrite,'(/,1x,a,2(/,1x,a,i0))')                        &
  "*** warning - number of molecular types mistmatch between FIELD and MPOLES !!! ***", &
  "***           FIELD  reports: ", ntpmls,                                             &
  "***           MPOLES reports: ", Nint(word_2_real(word,comm))

           Call error(623)
        End If

! read in molecular characteristics for every molecule

        Do itmols=1,ntpmls

           If (comm%idnode == 0 .and. l_top) Write(nrite,"(/,/,1x,'molecular species type',9x,i10)") itmols

! name of molecular species

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nmpldt,record,comm)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do
           Call strip_blanks(record)
           record1=word(1:Len_Trim(word)+1)//record ; Call lower_case(record1)
           record2=molnam(itmols) ;                   Call lower_case(record2)

           If (record1 == record2) Then
              If (comm%idnode == 0 .and. l_top) Write(nrite,"(/,1x,'name of species:',13x,a40)") molnam(itmols)
           Else
              If (comm%idnode == 0) Write(nrite,'(/,1x,a,i0)') &
  "*** warning - molecular names mistmatch between FIELD and MPOLES for type !!! *** ", itmols

              Call error(623)
           End If

! read molecular data

           Do

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nmpldt,record,comm)
                 If (.not.safe) Go To 2000
                 Call get_word(record,word) ; Call lower_case(word)
              End Do

! number of molecules of this type

              If (word(1:6) == 'nummol') Then

                 Call get_word(record,word)

                 If (nummols(itmols) == Nint(word_2_real(word,comm))) Then
                    If (comm%idnode == 0 .and. l_top) Write(nrite,"(/,1x,'number of molecules  ',10x,i10)") nummols(itmols)
                 Else
                    If (comm%idnode == 0) Write(nrite,'(/,1x,a,2(/,1x,a,i0))')         &
  "*** warning - number of molecules mistmatch between FIELD and MPOLES !!! ***", &
  "***           FIELD  reports: ", nummols(itmols),                              &
  "***           MPOLES reports: ", Nint(word_2_real(word,comm))

                    Call error(623)
                 End If

! read in atomic details

              Else If (word(1:5) == 'atoms') Then

                 Call get_word(record,word)

                 If (numsit(itmols) == Nint(word_2_real(word,comm))) Then
                    If (comm%idnode == 0 .and. l_top) Then
  Write(nrite,"(/,1x,'number of atoms/sites',10x,i10)") numsit(itmols)
  Write(nrite,"(/,1x,'atomic characteristics:', &
       & /,/,15x,'site',4x,'name',2x,'multipolar order',2x,'repeat'/)")
                    End If
                 Else
                    If (comm%idnode == 0) Write(nrite,'(/,1x,a,2(/,1x,a,i0))')                        &
  "*** warning - number of atoms/sites per molecule mistmatch between FIELD and MPOLES !!! ***", &
  "***           FIELD  reports: ", numsit(itmols),                                              &
  "***           MPOLES reports: ", Nint(word_2_real(word,comm))

                    Call error(623)
                 End If

! for every molecule of this type get site and atom description

                 ksite=0 ! reference point
                 Do isite=1,numsit(itmols)
                    If (ksite < numsit(itmols)) Then

! read atom name, highest pole order supplied, repeat

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nmpldt,record,comm)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       atom=word(1:8)

! read supplied pole order

                       Call get_word(record,word)
                       ordmpl=Abs(Nint(word_2_real(word,comm)))
                       indmpl=(ordmpl+3)*(ordmpl+2)*(ordmpl+1)/6

! read supplied repetition

                       Call get_word(record,word)
                       nrept=Abs(Nint(word_2_real(word,comm)))
                       If (nrept == 0) nrept=1

                       jsite=nsite+1
                       lsite=jsite+nrept-1

                       Do i=jsite,lsite
                          If (sitnam(i) /= atom) Then ! detect mish-mash
                             If (comm%idnode == 0) Write(nrite,'(/,1x,a,i0,a)') &
  "*** warning - site names mistmatch between FIELD and MPOLES for site ", ksite+1+i-jsite, " !!! ***"

                             Call error(623)
                          End If
                       End Do

! read supplied site polarisation and dumping factor

                       Call get_word(record,word) ; polarity=Abs(word_2_real(word,comm,0.0_wp))
                       Call get_word(record,word) ; dumping =Abs(word_2_real(word,comm,0.0_wp))

                       l_rsh=.true. ! regular or no shelling (Drude)
                       kshels=nshels
                       Do ishls=1,numshl(itmols) ! detect beyond charge shelling
                          kshels=kshels+1

                          isite2=nsite+lstshl(2,kshels)
                          If ((isite2 >= jsite .and. isite2 <= lsite)) Then
                             l_rsh=.false.
                             If (comm%idnode == 0) Then
                                If (ordmpl > 0) Write(nrite,'(/,1x,a)') &
  "*** warning - a shell (of a polarisable multipolar ion) can only bear a charge to emulate a self-iduced dipole !!! ***"
                                If (polarity > zero_plus) Write(nrite,'(/,1x,a)') &
  "*** warning - a shell (of a polarisable multipolar ion) cannot have its own associated polarisability !!! ***"
                                If (dumping  > zero_plus) Write(nrite,'(/,1x,a)') &
  "*** warning - a shell (of a polarisable multipolar ion) cannot have its own associated dumping factor !!! ***"
                             End If
                          End If
                       End Do

! get the min and max order defined for cores/nucleus, ignore irregular shells

                       If (l_rsh) Then
                          ordmpl_min=Min(ordmpl_min,ordmpl)
                          ordmpl_max=Max(ordmpl_max,ordmpl)
                       End If

                       If (comm%idnode == 0 .and. l_top) Then
                          If (l_rsh) Then
  Write(nrite,"(9x,i10,4x,a8,4x,i2,5x,i10,2f7.3)") ksite+1,atom,ordmpl,nrept,polarity,dumping
                          Else
  Write(nrite,"(9x,i10,4x,a8,4x,i2,5x,i10,2a)") ksite+1,atom,1,nrept,                                 &
                                             Merge(' *ignored* ','           ',polarity > zero_plus), &
                                             Merge(' *ignored* ','           ',dumping  > zero_plus)
                          End If
                       End If

! monopole=charge

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nmpldt,record,comm)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       sitmpl = 1

                       charge=word_2_real(word,comm)
                       chgsit(jsite:lsite)=charge
                       mpllfr(sitmpl,jsite:lsite)=charge
                       If (l_rsh) Then
                          plrsit(jsite:lsite)=polarity
                          dmpsit(jsite:lsite)=dumping
!                      Else ! initilised to zero in mpoles_module
                       End If

! sum absolute charges

                       sumchg=sumchg+Abs(charge)

! report

                       If (comm%idnode == 0 .and. l_top) &
  Write(nrite,"(3x,a12,3x,f10.5)") 'charge',charge

! higher poles counters

                       ordmpl_start = 0
                       ordmpl_next  = ordmpl_start+1
                       indmpl_start = sitmpl+1
                       indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                       Do While (ordmpl_next <= ordmpl)

! read line per pole order

                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nmpldt,record,comm)
                             If (.not.safe) Go To 2000
                             Call get_word(record,word)
                          End Do

! Only assign what FIELD says is needed or it is a shell

                          If (ordmpl_next <= Merge(mxompl,1,l_rsh)) Then

                             Do i=indmpl_start,indmpl_final
                                sitmpl = sitmpl+1
                                mpllfr(sitmpl,jsite:lsite)=word_2_real(word,comm)
                                Call get_word(record,word)
                             End Do

! report

                             If (comm%idnode == 0 .and. l_top) Then
                                If      (ordmpl_next == 1) Then
  Write(nrite,"(3x,a12,3x, 3f10.5)") 'dipole',       mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 2) Then
  Write(nrite,"(3x,a12,3x, 6f10.5)") 'quadrupole',   mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 3) Then
  Write(nrite,"(3x,a12,3x,10f10.5)") 'octupole',     mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 4) Then
  Write(nrite,"(3x,a12,3x,15f10.5)") 'hexadecapole', mpllfr(indmpl_start:indmpl_final,jsite)
                                End If
                             End If

! rescale poles values by their degeneracy

                             If (ordmpl_next > 1) Then
                                sitmpl=sitmpl-(indmpl_final-indmpl_start+1) ! rewind
                                Do i=ordmpl_next,0,-1
                                   l=ordmpl_next-i
                                   Do j=l,0,-1
                                      k=l-j

                                      scl=Exp(Factorial(ordmpl_next)-Factorial(k)-Factorial(j)-Factorial(i))
!                                      Write(*,*) i,j,k,Nint(scl)

                                      sitmpl = sitmpl+1 ! forward and apply scaling if degeneracy exists
                                      If (Nint(scl) /= 1) mpllfr(sitmpl,jsite:lsite)=mpllfr(sitmpl,jsite:lsite)/scl
                                   End Do
                                End Do
                             End If

                          Else

                             l_ord=.true.

! update actual order marker

                             sitmpl=indmpl_final

! report

                             If (comm%idnode == 0 .and. l_top) Then
                                If (l_rsh) Then
                                   If      (ordmpl_next == 1) Then
  Write(nrite,"(3x,a12,1x,a)") 'dipole',                 '     *** supplied but not required ***'
                                   Else If (ordmpl_next == 2) Then
  Write(nrite,"(3x,a12,1x,a)") 'quadrupole',             '     *** supplied but not required ***'
                                   Else If (ordmpl_next == 3) Then
  Write(nrite,"(3x,a12,1x,a)") 'octupole',               '     *** supplied but not required ***'
                                   Else If (ordmpl_next == 4) Then
  Write(nrite,"(3x,a12,1x,a)") 'hexadecapole',           '     *** supplied but not required ***'
                                   Else
  Write(nrite,"(3x,a12,i0,a)") 'pole order ',ordmpl_next,'     *** supplied but not required ***'
                                   End If
                                Else
                                   If      (ordmpl_next == 1) Then
  Write(nrite,"(3x,a12,1x,a)") 'dipole',                 '     *** supplied but ignored as invalid ***'
                                   Else If (ordmpl_next == 2) Then
  Write(nrite,"(3x,a12,1x,a)") 'quadrupole',             '     *** supplied but ignored as invalid ***'
                                   Else If (ordmpl_next == 3) Then
  Write(nrite,"(3x,a12,1x,a)") 'octupole',               '     *** supplied but ignored as invalid ***'
                                   Else If (ordmpl_next == 4) Then
  Write(nrite,"(3x,a12,1x,a)") 'hexadecapole',           '     *** supplied but ignored as invalid ***'
                                   Else
  Write(nrite,"(3x,a12,i0,a)") 'pole order ',ordmpl_next,'     *** supplied but ignored as invalid ***'
                                   End If
                                End If
                             End If

                          End If

! update poles counters

                          ordmpl_next  = ordmpl_next+1
                          indmpl_start = sitmpl+1
                          indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                       End Do

                       nsite=nsite+nrept
                       ksite=ksite+nrept

                       If (ksite == numsit(itmols)) nshels=kshels

                    End If
                 End Do

! finish of data for one molecular type

              Else If (word(1:6) == 'finish') Then

                 If (comm%idnode == 0) Then
                    Write(nrite,'(/,1x,3(a,i0),a)') &
  "*** warning - multipolar electrostatics requested up to order ", &
  mxompl, " with specified interactions up order ",                 &
  ordmpl_max," and least order ", ordmpl_min," !!! ***"
                    If (ordmpl_max*mxompl == 0) Write(nrite,'(1x,2a)') &
  "*** warning - multipolar electrostatics machinery to be used for ", &
  "monopoles only electrostatic interactions (point charges only) !!! ***"
                    If (ordmpl_max > 4) Write(nrite,'(1x,2a)')     &
  "*** warning - electrostatic interactions beyond hexadecapole ", &
  "order can not be considered and are thus ignored !!! ***"
                 End If

                 Go To 1000

              Else

! error exit for unidentified directive in molecular data

                 Call strip_blanks(record)
                 If (comm%idnode == 0) Write(nrite,'(/,1x,2a)') word(1:Len_Trim(word)+1),record
                 Call error(12)

              End If

           End Do

! just finished with this type molecule data

1000       Continue

        End Do

! close MPOLES data file

     Else If (word(1:5) == 'close') Then

        If (comm%idnode == 0) Close(Unit=nmpldt)

! EXIT IF ALL IS OK

        Return

     Else

! error exit for unidentified directive

        If (comm%idnode == 0) Write(nrite,'(/,1x,a)') word(1:Len_Trim(word))
        Call error(4)

     End If

  End Do

  Return

2000 Continue

  If (comm%idnode == 0) Close(Unit=nmpldt)
  Call error(52)

End Subroutine read_mpoles



End Module mpole
