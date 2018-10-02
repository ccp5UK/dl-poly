Module mpole
  Use kinds, Only : wp,wi
  Use setup, Only : mxatdm,mxatms,sqrpi,r4pie0,zero_plus,nmpldt
  Use site, Only : site_type
  Use core_shell, Only : core_shell_type
  Use parse
  Use comms, Only : comms_type
  Use numerics, Only : factorial
  Use errors_warnings, Only : error,warning,info

  Implicit None

  Private

  !> A derived data type for defining the rotation matrix for multipolar sites
  Type, Public :: rot_mat
     Real( Kind = wp ), Dimension(1:9) :: mtrxa
     Real( Kind = wp ), Dimension(1:3) :: p1,p2
     Integer                           :: flag,mbnd(1:2)
  End Type rot_mat

  ! Polarisation keys
  !> Default polarisation, unscreened & undamped - iAMOEBA like
  Integer( Kind = wi ), Parameter, Public :: POLARISATION_DEFAULT = 0
  !> CHARMM polarisation:
  !>
  !> - $q_\text{shell} = -\lvert q_core \rvert * \sqrt(\alpha * k)
  !> - $k_\text{CHARMM} = 1000 \, \mathrm{kcal} \mathrm{mol}^−1 \mathrm{Å}^−2
  Integer( Kind = wi ), Parameter, Public :: POLARISATION_CHARMM = 1

  !> Type containing multipole data
  Type, Public :: mpole_type
    Private

    !> Type of inducible (self-polarisation) scheme
    Integer( Kind = wi ), Public :: key = POLARISATION_DEFAULT

    !> Default thole dumping for CHARMM representation
    Real( Kind = wp ), Public :: thole  = 1.3_wp

    ! Variables for multipolar interactions
    !> Mappings from three indices multipole to a one index multipole
    Integer( Kind = wi ), Allocatable, Public :: map(:,:,:)
    Integer( Kind = wi ), Allocatable, Public :: ltg(:)
    !> Rotation counter flag
    Integer( Kind = wi ), Allocatable, Public :: flg(:)
    !> Bonded connectivity
    Integer( Kind = wi ), Allocatable, Public :: ltp(:,:)
    !> CHARMM core-shell screened electrostatics induction list
    Integer( Kind = wi ), Allocatable, Public :: charmm(:,:)

    !> Local frame
    Real( Kind = wp ), Allocatable, Public :: local_frame(:,:)
    !> Global frame
    Real( Kind = wp ), Allocatable, Public :: global_frame(:,:)

    !> Induced dipole polarisation for sites (inversed if non-zero)
    Real( Kind = wp ), Allocatable, Public :: polarisation_site(:)
    !> Induced dipole polarisation for atoms (inversed if non-zero)
    Real( Kind = wp ), Allocatable, Public :: polarisation_atom(:)

    !> thole dumping coefficient/factor (for self-polarisation) for sites
    Real( Kind = wp ), Allocatable, Public :: dump_site(:)
    !> thole dumping coefficient/factor (for self-polarisation) for atoms
    Real( Kind = wp ), Allocatable, Public :: dump_atom(:)

    !> Rotation matrices
    Type( rot_mat ),   Allocatable, Public :: rotation(:)

    !> Infinitesimal rotations about x
    Real( Kind = wp ), Allocatable, Public :: rotation_x(:,:)
    !> Infinitesimal rotations about y
    Real( Kind = wp ), Allocatable, Public :: rotation_y(:,:)
    !> Infinitesimal rotations about z
    Real( Kind = wp ), Allocatable, Public :: rotation_z(:,:)

    !> Torques due to infinitesimal rotations about x
    Real( Kind = wp ), Allocatable, Public :: torque_x(:)
    !> Torques due to infinitesimal rotations about y
    Real( Kind = wp ), Allocatable, Public :: torque_y(:)
    !> Torques due to infinitesimal rotations about z
    Real( Kind = wp ), Allocatable, Public :: torque_z(:)

    !> n choose k values, used in computing the reciprocal space Ewald sum field
    Real ( Kind = wp ), Allocatable, Public :: n_choose_k(:,:)

    !> Maxmimum number of multipoles
    Integer( Kind = wi ), Public :: max_mpoles
    !> Maxmimum multipolar order
    Integer( Kind = wi ), Public :: max_order

  Contains
    Private

    Procedure, Public :: init => allocate_mpoles_arrays
    Final :: cleanup
  End Type mpole_type

  Public :: read_mpoles

Contains

  Subroutine allocate_mpoles_arrays(T,max_site,max_exclude,mxatdm,bspline,mxatms)
  Class( mpole_type ) :: T
    Integer( Kind = wi ), Intent( In    ) :: max_site
    Integer( Kind = wi ), Intent( In    ) :: max_exclude
    Integer( Kind = wi ), Intent( In    ) :: mxatdm
    Integer( Kind = wi ), Intent( In    ) :: bspline
    Integer( Kind = wi ), Intent( In    ) :: mxatms

    Integer           :: n,k,om1,numpl,fail(1:9)

    ! No MPOLES file read <= no multipoles directive in FIELD
    If (T%max_mpoles < 1) Return

    om1 = T%max_order + 1
    numpl = (3**om1 - 1)/2

    fail = 0

    Allocate (T%map(0:T%max_order,0:T%max_order,0:T%max_order),T%ltg(1:numpl), stat=fail(1))
    Allocate (T%flg(1:mxatdm),T%ltp(0:max_exclude,1:mxatdm), stat=fail(2))
    If (T%key == 1) Then
      Allocate (T%charmm(0:max_exclude,1:mxatdm), stat=fail(3))
    End If
    Allocate (T%local_frame(1:T%max_mpoles,1:max_site),T%global_frame(1:T%max_mpoles,1:mxatms), stat=fail(4))
    Allocate (T%polarisation_site(1:max_site),T%polarisation_atom(1:mxatms), stat=fail(5))
    Allocate (T%dump_site(1:max_site),T%dump_atom(1:mxatms), stat=fail(6))
    Allocate (T%rotation(1:mxatdm),T%n_choose_k(0:bspline,0:bspline), stat=fail(7))
    Allocate (T%torque_x(1:mxatdm),T%torque_y(1:mxatdm),T%torque_z(1:mxatdm), stat=fail(8))
    Allocate (T%rotation_x(1:T%max_mpoles,1:mxatms), &
      T%rotation_y(1:T%max_mpoles,1:mxatms), &
      T%rotation_z(1:T%max_mpoles,1:mxatms), stat=fail(9))

    If (Any(fail > 0)) Call error(1025)

    T%flg = 0 ; T%ltp = 0
    If (T%key == 1) Then
      T%charmm = 0
    End If

    T%local_frame = 0.0_wp ; T%global_frame = 0.0_wp
    T%polarisation_site = 0.0_wp ; T%polarisation_atom = 0.0_wp
    T%dump_site = 0.0_wp ; T%dump_atom = 0.0_wp

    Do n=1,mxatdm
      T%rotation(n)%mtrxa = 0.0_wp
      T%rotation(n)%p1    = 0.0_wp
      T%rotation(n)%p2    = 0.0_wp
      T%rotation(n)%flag  = 0
      T%rotation(n)%mbnd  = 0
    End Do

    T%torque_x = 0.0_wp ; T%torque_y = 0.0_wp ; T%torque_z = 0.0_wp

    T%rotation_x = 0.0_wp ; T%rotation_y = 0.0_wp ; T%rotation_z = 0.0_wp

    ! Build the multipole map (polymap) and compute the constants T%n_choose_k
    ! Also build the map (T%ltg) that converts between index of a local
    ! multipole to the index of the corresponding global multipole
    T%map(0,0,0)=1

    T%ltg(1)=1

    If (T%max_order >= 1) Then
      T%map(1,0,0)=2 ; T%map(0,1,0)=3 ; T%map(0,0,1)=4

      T%ltg(2)=2 ; T%ltg(3)=3 ; T%ltg(4)=4
    End If

    If (T%max_order >= 2) Then
      T%map(2,0,0)=5 ; T%map(1,1,0)=6 ; T%map(1,0,1)=7
      T%map(0,2,0)=8 ; T%map(0,1,1)=9 ; T%map(0,0,2)=10

      T%ltg(5) =5 ; T%ltg(6) =6 ; T%ltg(7) =7
      T%ltg(8) =6 ; T%ltg(9) =8 ; T%ltg(10)=9
      T%ltg(11)=7 ; T%ltg(12)=9 ; T%ltg(13)=10
    End If

    If (T%max_order >= 3) Then
      T%map(3,0,0)=11 ; T%map(2,1,0)=12 ; T%map(2,0,1)=13
      T%map(1,2,0)=14 ; T%map(1,1,1)=15 ; T%map(1,0,2)=16
      T%map(0,3,0)=17 ; T%map(0,2,1)=18 ; T%map(0,1,2)=19
      T%map(0,0,3)=20

      T%ltg(14)=11 ; T%ltg(15)=12 ; T%ltg(16)=13
      T%ltg(17)=12 ; T%ltg(18)=14 ; T%ltg(19)=15
      T%ltg(20)=13 ; T%ltg(21)=15 ; T%ltg(22)=16
      T%ltg(23)=12 ; T%ltg(24)=14 ; T%ltg(25)=15
      T%ltg(26)=14 ; T%ltg(27)=17 ; T%ltg(28)=18
      T%ltg(29)=15 ; T%ltg(30)=18 ; T%ltg(31)=19
      T%ltg(32)=13 ; T%ltg(33)=15 ; T%ltg(34)=16
      T%ltg(35)=15 ; T%ltg(36)=18 ; T%ltg(37)=19
      T%ltg(38)=16 ; T%ltg(39)=19 ; T%ltg(40)=20
    End If

    If (T%max_order >=4) Then
      T%map(4,0,0)=21 ; T%map(3,1,0)=22 ; T%map(3,0,1)=23
      T%map(2,2,0)=24 ; T%map(2,1,1)=25 ; T%map(2,0,2)=26
      T%map(1,3,0)=27 ; T%map(1,2,1)=28 ; T%map(1,1,2)=29
      T%map(1,0,3)=30 ; T%map(0,4,0)=31 ; T%map(0,3,1)=32
      T%map(0,2,2)=33 ; T%map(0,1,3)=34 ; T%map(0,0,4)=35

      T%ltg(41)=21 ; T%ltg(42)=22 ; T%ltg(43)=23
      T%ltg(44)=22 ; T%ltg(45)=24 ; T%ltg(46)=25
      T%ltg(47)=23 ; T%ltg(48)=25 ; T%ltg(49)=26
      T%ltg(50)=22 ; T%ltg(51)=24 ; T%ltg(52)=25
      T%ltg(53)=24 ; T%ltg(54)=27 ; T%ltg(55)=28
      T%ltg(56)=25 ; T%ltg(57)=28 ; T%ltg(58)=29
      T%ltg(59)=23 ; T%ltg(60)=25 ; T%ltg(61)=26
      T%ltg(62)=25 ; T%ltg(63)=28 ; T%ltg(64)=29
      T%ltg(65)=26 ; T%ltg(66)=29 ; T%ltg(67)=30
      T%ltg(68)=22 ; T%ltg(69)=24 ; T%ltg(70)=25
      T%ltg(71)=24 ; T%ltg(72)=27 ; T%ltg(73)=28
      T%ltg(74)=25 ; T%ltg(75)=28 ; T%ltg(76)=29
      T%ltg(77)=24 ; T%ltg(78)=27 ; T%ltg(79)=28
      T%ltg(80)=27 ; T%ltg(81)=31 ; T%ltg(82)=32
      T%ltg(83)=28 ; T%ltg(84)=32 ; T%ltg(85)=33
      T%ltg(86)=25 ; T%ltg(87)=28 ; T%ltg(88)=29
      T%ltg(89)=28 ; T%ltg(90)=32 ; T%ltg(91)=33
      T%ltg(92)=29 ; T%ltg(93)=33 ; T%ltg(94)=34
      T%ltg(95)=23 ; T%ltg(96)=25 ; T%ltg(97)=26
      T%ltg(98)=25 ; T%ltg(99)=28 ; T%ltg(100)=29

      T%ltg(101)=26 ; T%ltg(102)=29 ; T%ltg(103)=30
      T%ltg(104)=25 ; T%ltg(105)=28 ; T%ltg(106)=29
      T%ltg(107)=28 ; T%ltg(108)=32 ; T%ltg(109)=33
      T%ltg(110)=29 ; T%ltg(111)=33 ; T%ltg(112)=34
      T%ltg(113)=26 ; T%ltg(114)=29 ; T%ltg(115)=30
      T%ltg(116)=29 ; T%ltg(117)=33 ; T%ltg(118)=34
      T%ltg(119)=30 ; T%ltg(120)=34 ; T%ltg(121)=35
    End If

    If (T%max_order < 0 .or. T%max_order >=5) Call error(2071)

    ! compute n choose k
    Do k = 0, bspline
      Do n = 0, bspline
        T%n_choose_k(n,k) = Exp(Factorial(n)-Factorial(n-k)-Factorial(k))
      End Do
    End Do
  End Subroutine allocate_mpoles_arrays

  Subroutine cleanup(T)
    Type( mpole_type ) :: T

    If (Allocated(T%map)) Then
      Deallocate(T%map)
    End If
    If (Allocated(T%ltg)) Then
      Deallocate(T%ltg)
    End If
    If (Allocated(T%flg)) Then
      Deallocate(T%flg)
    End If
    If (Allocated(T%ltp)) Then
      Deallocate(T%ltp)
    End If
    If (Allocated(T%charmm)) Then
      Deallocate(T%charmm)
    End If

    If (Allocated(T%local_frame)) Then
      Deallocate(T%local_frame)
    End If
    If (Allocated(T%global_frame)) Then
      Deallocate(T%global_frame)
    End If

    If (Allocated(T%polarisation_site)) Then
      Deallocate(T%polarisation_site)
    End If
    If (Allocated(T%polarisation_atom)) Then
      Deallocate(T%polarisation_atom)
    End If

    If (Allocated(T%dump_site)) Then
      Deallocate(T%dump_site)
    End If
    If (Allocated(T%dump_atom)) Then
      Deallocate(T%dump_atom)
    End If

    If (Allocated(T%rotation)) Then
      Deallocate(T%rotation)
    End If

    If (Allocated(T%rotation_x)) Then
      Deallocate(T%rotation_x)
    End If
    If (Allocated(T%rotation_y)) Then
      Deallocate(T%rotation_y)
    End If
    If (Allocated(T%rotation_z)) Then
      Deallocate(T%rotation_z)
    End If

    If (Allocated(T%torque_x)) Then
      Deallocate(T%torque_x)
    End If
    If (Allocated(T%torque_y)) Then
      Deallocate(T%torque_y)
    End If
    If (Allocated(T%torque_z)) Then
      Deallocate(T%torque_z)
    End If

    If (Allocated(T%n_choose_k)) Then
      Deallocate(T%n_choose_k)
    End If
  End Subroutine cleanup

  Subroutine read_mpoles(l_top,sumchg,cshell,sites,mpoles,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for reading in the molecular mulitpole
  ! specifications of the system to be simulated
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( In    ) :: l_top
    Real( Kind = wp ),  Intent( InOut ) :: sumchg
    Type( site_type ), Intent( InOut ) :: sites
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( comms_type ), Intent( InOut ) :: comm

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

    Character( Len = 256 ) :: message,messages(3)

  ! open MPOLES data file

    If (comm%idnode == 0) Then
       Open(Unit=nmpldt, File = 'MPOLES', Status = 'old')
    End If
    Call info('electrostatics multipoles specification',.true.)
    If (.not.l_top) Then
      Call info('detailed specification opted out',.true.)
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

          If (sites%ntype_mol == Nint(word_2_real(word))) Then
            Write(message,'(a,i10)') 'number of molecular types ',sites%ntype_mol
            Call info(message,.true.)
          Else
            Write(message,'(2(a,i0),a)') &
              'number of molecular types mismatch between FIELD(',sites%ntype_mol, &
              ') and MPOLES(',Nint(word_2_real(word)),')'
            Call warning(message,.true.)
            Call error(623)
          End If

  ! read in molecular characteristics for every molecule

          Do itmols=1,sites%ntype_mol

             If (l_top) Then
               Write(message,'(a,i10)') 'molecular species type ',itmols
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
             record2=sites%mol_name(itmols) ;                   Call lower_case(record2)

             If (record1 == record2) Then
               If (l_top) Then
                 Write(message,'(2a)') 'name of species: ',Trim(sites%mol_name(itmols))
                 Call info(message,.true.)
               End If
             Else
               Call warning('molecular names mismatch between FIELD and MPOLES for type',.true.)
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

                   If (sites%num_mols(itmols) == Nint(word_2_real(word))) Then
                     If (l_top) Then
                       Write(message,'(a,i10)') 'number of molecules ',sites%num_mols(itmols)
                       Call info(message,.true.)
                     End If
                   Else
                     Write(message,'(2(a,i0),a)') &
                       'number of molecular types mismatch between FIELD(',sites%ntype_mol, &
                       ') and MPOLES(',Nint(word_2_real(word)),')'
                     Call warning(message,.true.)
                     Call error(623)
                   End If

  ! read in atomic details

                Else If (word(1:5) == 'atoms') Then

                   Call get_word(record,word)

                   If (sites%num_site(itmols) == Nint(word_2_real(word))) Then
                     If (l_top) Then
                       Write(messages(1),'(a,i10)') 'number of atoms/sites ',sites%num_site(itmols)
                       Write(messages(2),'(a)') 'atomic characteristics:'
                       Write(messages(3),'(8x,a4,4x,a4,2x,a16,2x,a6)') &
                         'site','name','multipolar order','repeat'
                       Call info(messages,3,.true.)
                     End If
                   Else
                     Write(message,'(2(a,i0),a)') &
                       'number of molecular types mismatch between FIELD(',sites%ntype_mol, &
                       ') and MPOLES(',Nint(word_2_real(word)),')'
                     Call warning(message,.true.)
                     Call error(623)
                   End If

  ! for every molecule of this type get site and atom description

                   ksite=0 ! reference point
                   Do isite=1,sites%num_site(itmols)
                      If (ksite < sites%num_site(itmols)) Then

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
                           If (sites%site_name(i) /= atom) Then ! detect mish-mash
                             Write(message,'(a,i0)') &
                               'site names mismatch between FIELD and MPOLES for site ',ksite+1+i-jsite
                             Call warning(message,.true.)
                             Call error(623)
                           End If
                         End Do

  ! read supplied site polarisation and dumping factor

                         Call get_word(record,word) ; polarity=Abs(word_2_real(word,0.0_wp))
                         Call get_word(record,word) ; dumping =Abs(word_2_real(word,0.0_wp))

                         l_rsh=.true. ! regular or no shelling (Drude)
                         kshels=nshels
                         Do ishls=1,cshell%numshl(itmols) ! detect beyond charge shelling
                            kshels=kshels+1

                            isite2=nsite+cshell%lstshl(2,kshels)
                            If ((isite2 >= jsite .and. isite2 <= lsite)) Then
                              l_rsh=.false.
                              If (comm%idnode == 0) Then
                                If (ordmpl > 0) Then
                                  Call warning( &
                                    'a shell (of a polarisable multipolar ion)' &
                                    //'can only bear a charge to emulate a' &
                                    //'self-iduced dipole',.true.)
                                End If
                                If (polarity > zero_plus) Then
                                  Call warning( &
                                    'a shell (of a polarisable multipolar ion)' &
                                    //'cannot have its own associated polarisability', &
                                    .true.)
                                End If
                                If (dumping  > zero_plus) Then
                                  Call warning( &
                                    'a shell (of a polarisable multipolar ion)' &
                                    //'cannot have its own associated dumping factor', &
                                    .true.)
                                End IF
                              End If
                            End If
                         End Do

  ! get the min and max order defined for cores/nucleus, ignore irregular shells

                         If (l_rsh) Then
                            ordmpl_min=Min(ordmpl_min,ordmpl)
                            ordmpl_max=Max(ordmpl_max,ordmpl)
                         End If

                         If (l_top) Then
                           If (l_rsh) Then
                             Write(message,'(2x,i10,4x,a8,4x,i2,5x,i10,2f7.3)') &
                               ksite+1,atom,ordmpl,nrept,polarity,dumping
                           Else
                             Write(message,'(2x,i10,4x,a8,4x,i2,5x,i10,2a)') &
                               ksite+1,atom,ordmpl,nrept,polarity,dumping, &
                               Merge(' *ignored* ','           ',polarity>zero_plus), &
                               Merge(' *ignored* ','           ',dumping >zero_plus)
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

                         sites%charge_site(jsite:lsite)=charge
                         mpoles%local_frame(sitmpl,jsite:lsite)=charge
                         If (l_rsh) Then
                            mpoles%polarisation_site(jsite:lsite)=polarity
                            mpoles%dump_site(jsite:lsite)=dumping
  !                      Else ! initilised to zero in mpoles_module
                         End If

  ! sum absolute charges

                         sumchg=sumchg+Abs(charge)

  ! report

                         If (l_top) Then
                           Write(message,'(2x,a,f10.5)') 'charge ',charge
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

                            If (ordmpl_next <= Merge(mpoles%max_order,1,l_rsh)) Then

                              Do i=indmpl_start,indmpl_final
                                sitmpl = sitmpl+1
                                mpoles%local_frame(sitmpl,jsite:lsite)=word_2_real(word)
                                Call get_word(record,word)
                              End Do

  ! report

                              If (l_top) Then
                                If      (ordmpl_next == 1) Then
                                  Write(message,'(2x,a12,1x,3f10.5)') 'dipole', &
                                    mpoles%local_frame(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 2) Then
                                  Write(message,'(2x,a12,1x,6f10.5)') 'quadrupole', &
                                    mpoles%local_frame(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 3) Then
                                  Write(message,'(2x,a12,1x,10f10.5)') 'octupole', &
                                    mpoles%local_frame(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 4) Then
                                  Write(message,'(2x,a12,1x,15f10.5)') 'hexadecapole', &
                                    mpoles%local_frame(indmpl_start:indmpl_final,jsite)
                                End If
                                Call info(message,.true.)
                              End If

  ! rescale poles values by their degeneracy

                               If (ordmpl_next > 1) Then
                                  sitmpl=sitmpl-(indmpl_final-indmpl_start+1) ! rewind
                                  Do i=ordmpl_next,0,-1
                                     l=ordmpl_next-i
                                     Do j=l,0,-1
                                        k=l-j

                                        scl=Exp(factorial(ordmpl_next)-factorial(k)-factorial(j)-factorial(i))
                                        sitmpl = sitmpl+1 ! forward and apply scaling if degeneracy exists
                                        If (Nint(scl) /= 1) Then
                                          mpoles%local_frame(sitmpl,jsite:lsite)= &
                                            mpoles%local_frame(sitmpl,jsite:lsite)/scl
                                        End If
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
                                     Write(message,'(2x,a12,1x,a)') &
                                       'dipole',' supplied but not required'
                                   Else If (ordmpl_next == 2) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'quadrupole',' supplied but not required'
                                   Else If (ordmpl_next == 3) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'octupole',' supplied but not required'
                                   Else If (ordmpl_next == 4) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'hexadecapole',' supplied but not required'
                                   Else
                                     Write(message,'(2x,a12,i0,a)') &
                                       'pole order ',ordmpl_next,' supplied but not required'
                                   End If
                                 Else
                                   If      (ordmpl_next == 1) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'dipole',' supplied but ignored as invalid'
                                   Else If (ordmpl_next == 2) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'quadrupole',' supplied but ignored as invalid'
                                   Else If (ordmpl_next == 3) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'octupole',' supplied but ignored as invalid'
                                   Else If (ordmpl_next == 4) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'hexadecapole',' supplied but ignored as invalid'
                                   Else
                                     Write(message,'(2x,a11,i0,a)') &
                                       'pole order ',ordmpl_next,' supplied but ignored as invalid'
                                   End If
                                 End If
                                 Call info(message,.true.)
                               End If

                            End If

  ! update poles counters

                            ordmpl_next  = ordmpl_next+1
                            indmpl_start = sitmpl+1
                            indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                         End Do

                         nsite=nsite+nrept
                         ksite=ksite+nrept

                         If (ksite == sites%num_site(itmols)) nshels=kshels

                      End If
                   End Do

  ! finish of data for one molecular type

                Else If (word(1:6) == 'finish') Then

                  Write(message,'(3(a,i0))') &
                    'multipolar electrostatics requested up to order ', &
                    mpoles%max_order, ' with specified interactions up order ',  &
                    ordmpl_max,' and least order ', ordmpl_min
                  Call warning(message,.true.)

                  If (ordmpl_max*mpoles%max_order == 0) Then
                    Call warning( &
                      'multipolar electrostatics machinery to be used for monompoles ' &
                      //'only electrostatic interactions (point charges only)', &
                      .true.)
                  End If
                  If (ordmpl_max > 4) Then
                    Call warning( &
                      'electrostatic interactions beyond hexadecapole order can ' &
                      //'not be considered and are thus ignored',.true.)
                  End If

                  Go To 1000

                Else

  ! error exit for unidentified directive in molecular data

                   Call strip_blanks(record)
                   Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                   Call info(message,.true.)
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

          Call info(word(1:Len_Trim(word)),.true.)
          Call error(4)

       End If

    End Do

    Return

  2000 Continue

    If (comm%idnode == 0) Close(Unit=nmpldt)
    Call error(52)

  End Subroutine read_mpoles

End Module mpole
