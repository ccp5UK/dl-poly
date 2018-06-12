Module mpole
  Use kinds, Only : wp
  Use setup, Only : mxatdm, mxatms, mxexcl, mximpl, mxompl, site_data%max_site, mxspl, &
                           sqrpi,r4pie0,zero_plus,nrite,nmpldt
  Use configuration,Only : natms
  Use site, Only : site_type
  Use core_shell, Only : numshl,lstshl
  Use parse
  Use comms, Only : comms_type
  Use numerics, Only : factorial
  Use errors_warnings, Only : error,warning,info

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
  Integer,           Allocatable, Save :: lchatm(:,:)             ! CHARMM core-shell screened electrostatics induction neigh%list

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
    Real( Kind = wp ) :: gearp(1:7),aspcp(1:7)

    If (mximpl < 1) Return ! no MPOLES file read <= no multipoles directive in FIELD

    om1 = mxompl + 1
    numpl = (3**om1 - 1)/2

    fail = 0

    Allocate (mplmap(0:mxompl,0:mxompl,0:mxompl),mplltg(1:numpl),  Stat = fail(1))
    Allocate (mplflg(1:mxatdm),ltpatm(0:mxexcl,1:mxatdm),          Stat = fail(2))
    If (keyind == 1) &
    Allocate (lchatm(0:mxexcl,1:mxatdm),                           Stat = fail(3))
    Allocate (mpllfr(1:mximpl,1:site_data%max_site),mplgfr(1:mximpl,1:mxatms), Stat = fail(4))
    Allocate (plrsit(1:site_data%max_site),plratm(1:mxatms),                   Stat = fail(5))
    Allocate (dmpsit(1:site_data%max_site),dmpatm(1:mxatms),                   Stat = fail(6))
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

    Real( Kind = wp )      :: charge,scl,polarity,dumping
    Character( Len = 256 ) :: message

  ! open MPOLES data file

  If (comm%idnode == 0) Then
    Open(Unit=nmpldt, File = 'MPOLES', Status = 'old')
  End If
  Write(message,"('ELECTROSTATICS MULTIPOLES SPECIFICATION')")
  Call info(message,.true.)
  If (.not.l_top) Then
    Write(message,"('detailed specification opted out')")
    Call info(message,.true.)
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

          If (site_data%ntype_mol == Nint(word_2_real(word))) Then
             Write(message,"('number of molecular types',6x,i10)") site_data%ntype_mol
             Call info(message,.true.)
          Else
            Call warning("number of molecular types mistmatch between FIELD and MPOLES",.true.)
            Write(message,'(2(a,i0))') &
            "FIELD  reports: ", site_data%ntype_mol, ", MPOLES reports: ", Nint(word_2_real(word))
            Call info(message,.true.)
            Call error(623,master_only=.true.)
          End If

  ! read in molecular characteristics for every molecule

          Do itmols=1,site_data%ntype_mol

            If (l_top) Then
              Write(message,"('molecular species type',9x,i10)") itmols
              Call info(message,.true.)
            End If

  ! name of molecular species

             word(1:1)='#'
             Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                Call get_line(safe,nmpldt,record,comm)
                If (.not.safe) Go To 2000
                Call get_word(record,word)
             End Do
             Call strip_blanks(record)
             record1=word(1:Len_Trim(word)+1)//record ; Call lower_case(record1)
             record2=site_data%mol_name(itmols) ;                   Call lower_case(record2)

             If (record1 == record2) Then
               If (l_top) Then
                 Write(message,"('name of species:',13x,a40)") site_data%mol_name(itmols)
                 Call info(message,.true.)
               End If
             Else
               Write(message,'(a)') 'molecular names mistmatch between FIELD and MPOLES for type'
               Call error(623,message,.true.)
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

                   If (site_data%num_mols(itmols) == Nint(word_2_real(word))) Then
                     If (l_top) Then
                       Write(message,"('number of molecules  ',10x,i10)") site_data%num_mols(itmols)
                       Call info(message,.true.)
                     End If
                   Else
                     Call warning("number of molecules mistmatch between FIELD and MPOLES",.true.)
                     Write(message,'(2(a,i0))') &
                     "FIELD  reports: ", site_data%num_mols(itmols), ", MPOLES reports: ", Nint(word_2_real(word))
                     Call info(message,.true.)
                     Call error(623,master_only=.true.)
                   End If

  ! read in atomic details

                Else If (word(1:5) == 'atoms') Then

                   Call get_word(record,word)

                   If (site_data%num_site(itmols) == Nint(word_2_real(word))) Then
                      If (l_top) Then
                        Write(message,"('number of atoms/sites',10x,i10)") site_data%num_site(itmols)
                        Call info(message,.true.)
                        Write(message,"('atomic characteristics:', &
                          & /,/,15x,'site',4x,'name',2x,'multipolar order',2x,'repeat'/)")
                        Call info(message,.true.)
                      End If
                   Else
                     Call warning("number of molecules mistmatch between FIELD and MPOLES",.true.)
                     Write(message,'(2(a,i0))') &
                     "FIELD  reports: ", site_data%num_mols(itmols), ", MPOLES reports: ", Nint(word_2_real(word))
                     Call info(message,.true.)
                     Call error(623,master_only=.true.)
                   End If

  ! for every molecule of this type get site and atom description

                   ksite=0 ! reference point
                   Do isite=1,site_data%num_site(itmols)
                      If (ksite < site_data%num_site(itmols)) Then

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
                         ordmpl=Abs(Nint(word_2_real(word)))
                         indmpl=(ordmpl+3)*(ordmpl+2)*(ordmpl+1)/6

  ! read supplied repetition

                         Call get_word(record,word)
                         nrept=Abs(Nint(word_2_real(word)))
                         If (nrept == 0) nrept=1

                         jsite=nsite+1
                         lsite=jsite+nrept-1

                         Do i=jsite,lsite
                           If (site_data%site_name(i) /= atom) Then ! detect mish-mash
                             Write(message,'(a,i0)') &
                               "site names mistmatch between FIELD and MPOLES for site ", ksite+1+i-jsite
                             Call error(623,message,.true.)
                           End If
                         End Do

  ! read supplied site polarisation and dumping factor

                         Call get_word(record,word) ; polarity=Abs(word_2_real(word,0.0_wp))
                         Call get_word(record,word) ; dumping =Abs(word_2_real(word,0.0_wp))

                         l_rsh=.true. ! regular or no shelling (Drude)
                         kshels=nshels
                         Do ishls=1,numshl(itmols) ! detect beyond charge shelling
                            kshels=kshels+1

                            isite2=nsite+lstshl(2,kshels)
                            If ((isite2 >= jsite .and. isite2 <= lsite)) Then
                               l_rsh=.false.
                               If (ordmpl > 0) Write(message,'(a)') &
                                 "a shell (of a polarisable multipolar ion) can only bear a charge to emulate a self-iduced dipole"
                               Call warning(message,.true.)
                               If (polarity > zero_plus) Write(message,'(a)') &
                                 "a shell (of a polarisable multipolar ion) cannot have its own associated polarisability"
                               Call warning(message,.true.)
                               If (dumping  > zero_plus) Write(message,'(a)') &
                                 "a shell (of a polarisable multipolar ion) cannot have its own associated dumping factor"
                               Call warning(message,.true.)
                            End If
                         End Do

  ! get the min and max order defined for cores/nucleus, ignore irregular shells

                         If (l_rsh) Then
                            ordmpl_min=Min(ordmpl_min,ordmpl)
                            ordmpl_max=Max(ordmpl_max,ordmpl)
                         End If

                         If (l_top) Then
                           If (l_rsh) Then
                             Write(message,"(9x,i10,4x,a8,4x,i2,5x,i10,2f7.3)") ksite+1,atom,ordmpl,nrept,polarity,dumping
                             Call info(message,.true.)
                           Else
                             Write(message,"(9x,i10,4x,a8,4x,i2,5x,i10,2a)") ksite+1,atom,1,nrept,                                 &
                               Merge(' *ignored* ','           ',polarity > zero_plus), &
                               Merge(' *ignored* ','           ',dumping  > zero_plus)
                             Call info(message,.true.)
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

                         charge=word_2_real(word)
                         site_data%charge_site(jsite:lsite)=charge
                         mpllfr(sitmpl,jsite:lsite)=charge
                         If (l_rsh) Then
                            plrsit(jsite:lsite)=polarity
                            dmpsit(jsite:lsite)=dumping
  !                      Else ! initilised to zero in mpoles_module
                         End If

  ! sum absolute charges

                         sumchg=sumchg+Abs(charge)

  ! report

                          If (l_top) Then
                            Write(message,"(3x,a12,3x,f10.5)") 'charge',charge
                            Call info(message,.true.)
                          End If

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
                                  mpllfr(sitmpl,jsite:lsite)=word_2_real(word)
                                  Call get_word(record,word)
                               End Do

  ! report

                               If (l_top) Then
                                 If      (ordmpl_next == 1) Then
                                   Write(message,"(3x,a12,3x, 3f10.5)") 'dipole',       mpllfr(indmpl_start:indmpl_final,jsite)
                                   Call info(message,.true.)
                                 Else If (ordmpl_next == 2) Then
                                   Write(message,"(3x,a12,3x, 6f10.5)") 'quadrupole',   mpllfr(indmpl_start:indmpl_final,jsite)
                                   Call info(message,.true.)
                                 Else If (ordmpl_next == 3) Then
                                   Write(message,"(3x,a12,3x,10f10.5)") 'octupole',     mpllfr(indmpl_start:indmpl_final,jsite)
                                   Call info(message,.true.)
                                 Else If (ordmpl_next == 4) Then
                                   Write(message,"(3x,a12,3x,15f10.5)") 'hexadecapole', mpllfr(indmpl_start:indmpl_final,jsite)
                                   Call info(message,.true.)
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

                               If (l_top) Then
                                 If (l_rsh) Then
                                   If      (ordmpl_next == 1) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'dipole', '     *** supplied but not required ***'
                                     Call info(message,.true.)
                                   Else If (ordmpl_next == 2) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'quadrupole', '     *** supplied but not required ***'
                                     Call info(message,.true.)
                                   Else If (ordmpl_next == 3) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'octupole', '     *** supplied but not required ***'
                                     Call info(message,.true.)
                                   Else If (ordmpl_next == 4) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'hexadecapole', '     *** supplied but not required ***'
                                     Call info(message,.true.)
                                   Else
                                     Write(message,"(3x,a12,i0,a)") &
                                       'pole order ',ordmpl_next,'     *** supplied but not required ***'
                                     Call info(message,.true.)
                                   End If
                                 Else
                                   If      (ordmpl_next == 1) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'dipole', '     *** supplied but ignored as invalid ***'
                                     Call info(message,.true.)
                                   Else If (ordmpl_next == 2) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'quadrupole', '     *** supplied but ignored as invalid ***'
                                     Call info(message,.true.)
                                   Else If (ordmpl_next == 3) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'octupole', '     *** supplied but ignored as invalid ***'
                                     Call info(message,.true.)
                                   Else If (ordmpl_next == 4) Then
                                     Write(message,"(3x,a12,1x,a)") &
                                       'hexadecapole', '     *** supplied but ignored as invalid ***'
                                     Call info(message,.true.)
                                   Else
                                     Write(message,"(3x,a12,i0,a)") &
                                       'pole order ',ordmpl_next,'     *** supplied but ignored as invalid ***'
                                     Call info(message,.true.)
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

                         If (ksite == site_data%num_site(itmols)) nshels=kshels

                      End If
                   End Do

  ! finish of data for one molecular type

                Else If (word(1:6) == 'finish') Then

                  Write(message,'(3(a,i0))') &
                    "multipolar electrostatics requested up to order ", &
                    mxompl, " with specified interactions up order ", &
                    ordmpl_max," and least order ", ordmpl_min
                  Call warning(message,.true.)
                  If (ordmpl_max*mxompl == 0) Then
                    Write(message,'(1x,2a)') &
                      "multipolar electrostatics machinery to be used for ", &
                      "monopoles only electrostatic interactions (point charges only)"
                    Call warning(message,.true.)
                  End If
                  If (ordmpl_max > 4) Then
                    Write(message,'(1x,2a)')     &
                      "electrostatic interactions beyond hexadecapole ", &
                      "order can not be considered and are thus ignored"
                    Call warning(message,.true.)
                  End If

                   Go To 1000

                Else

  ! error exit for unidentified directive in molecular data

                   Call strip_blanks(record)
                   Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                   Call error(12,message,.true.)

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

          Write(message,'(a)') word(1:Len_Trim(word))
          Call error(4,message,.true.)

       End If

    End Do

    Return

  2000 Continue

    If (comm%idnode == 0) Close(Unit=nmpldt)
    Call error(52)

  End Subroutine read_mpoles
End Module mpole
