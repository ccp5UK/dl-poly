Module development

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 development module
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov june 2013
  ! contrib   - i.j.bush november 2008
  ! contrib   - a.m.elena march 2016
  ! contrib   - a.m.elena february 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only: wp
  Use parse,           Only: get_line, get_word, lower_case, clean_string, word_2_real
  Use filename,        Only: file_type, FILE_CONTROL
  Use errors_warnings, Only: info, error , set_print_level
  Use iso_fortran_env, ONly: compiler_version
#ifdef OLDMPI
  Use comms,           Only: mpi_ver, mpi_subver, comms_type, gcheck
#else
  Use comms,           Only: mpi_ver, mpi_subver, lib_version, comms_type, gcheck
#endif

  Implicit None

  Private

  !> Logicals indicating whether tests should be run for
  !> corresponding module. Add as required
  Type, Public :: testing_type
     Logical :: configuration
     Logical :: dftb_library
   Contains
     Procedure :: all => set_all_tests_true
     !Procedure :: set => set_tests
  End Type testing_type

  !> Type containing development module variables
  Type, Public :: development_type
    Private
    !> OUTPUT redirection to the default output (screen)
    Logical, Public       :: l_scr = .false.
    !> avoid global safety checks (no elegant parallel failures)
    Logical, Public       :: l_fast = .false.
    !> OUTPUT inclusion of an extra last line with E_tot
    Logical, Public       :: l_eng = .false.
    !> REVIVE writing in ASCII (default is binary)
    Logical, Public       :: l_rout = .false.
    !> REVOLD reading in ASCII (default is binary)
    Logical, Public       :: l_rin = .false.
    !> translate CONFIG along a vector to CFGORG
    Logical, Public       :: l_org = .false.
    !> CONFIG rescaling to CFGSCL after reading with termination
    Logical, Public       :: l_scl = .false.
    !> HISTORY generation after reading with termination
    Logical, Public       :: l_his = .false.
    !> termination flag
    Logical, Public       :: l_trm = .false.
    !> detailed timing
    Logical               :: l_tim = .false.
    !> no production of REVCON & REVIVE
    Logical, Public       :: l_tor = .false.
    !> check on minimum separation distance between VNL pairs at re/start
    Logical, Public       :: l_dis = .false.

    !> See whether the new control file format should be used
    Logical, Public        :: new_control = .false.

    !> CFGORG levcfg
    Integer, Public       :: lvcforg = -1
    !> reorigin vector
    Real(Kind=wp), Public :: xorg = 0.0_wp, yorg = 0.0_wp, zorg = 0.0_wp
    !> CFGSCL levcfg
    Integer, Public       :: lvcfscl = -1
    !> CFGSCL lattice parameters
    Real(Kind=wp), Public :: cels(1:9) = 0.0_wp
    !> l_dis default check condition
    Real(Kind=wp), Public :: r_dis = 0.5_wp
    !> Devel start time
    Real(Kind=wp), Public :: t_zero
    !> Unit testing
    Logical, Public :: run_unit_tests = .false.
    Type(testing_type), Public :: unit_test
    !> App testing
    Logical, Public :: run_app_tests = .false.
    Type(testing_type), Public :: app_test
 End Type development_type

  Public :: scan_development
  Public :: build_info

Contains

  Subroutine scan_development(devel, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 development subroutine for raw scanning the contents of the
    ! control file
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(development_type), Intent(InOut) :: devel
    Type(file_type),        Intent(InOut) :: files(:)
    Type(comms_type),       Intent(InOut) :: comm
    Integer :: print_level
    Character(Len=200) :: record
    Character(Len=40)  :: word
    Logical            :: carry, safe

    ! Set safe flag

    safe = .true.

    ! Open the simulation input file

    If (comm%idnode == 0) Inquire (File=files(FILE_CONTROL)%filename, Exist=safe)
    Call gcheck(comm, safe, "enforce")
    If (.not. safe) Then
      Return
    Else
      If (comm%idnode == 0) Then
        Open (Newunit=files(FILE_CONTROL)%unit_no, File=files(FILE_CONTROL)%filename, Status='old')
      End If
    End If

    Call get_line(safe, files(FILE_CONTROL)%unit_no, record, comm)
    If (safe) Then

      carry = .true.
      Do While (carry)

        Call get_line(safe, files(FILE_CONTROL)%unit_no, record, comm)
        If (.not. safe) Exit

        Call lower_case(record)
        Call get_word(record, word)

        ! read DEVELOPMENT option: OUTPUT to screen

        If (word(1:5) == 'l_scr') Then
          devel%l_scr = .true.
!          Call info('%%% OUTPUT redirected to the default output (screen) !!! %%%',.true.)
        Else If (word(1:6) == 'l_fast') Then
          devel%l_fast = .true.
          Call info('%%% speed up by avoiding global safety checks !!! %%%', .true.)
        Else If (word(1:5) == 'l_tim') Then
           devel%l_tim = .true.
        Else If (word(1:7) == 'l_print') Then
           Call get_word(record, word)
           print_level = nint(abs(word_2_real(word)))
           Call set_print_level(print_level)
        Else If (word(1:6) == 'finish') Then
          carry = .false.
        End If

      End Do
    End If
    If (comm%idnode == 0) Call files(FILE_CONTROL)%close ()

  End Subroutine scan_development

#ifdef HOST
#define __HOSTNAME__ HOST
#else
#define __HOSTNAME__ 'unknown'
#endif

#ifdef BUILDER
#define __BUILDER__ BUILDER
#else
#define __BUILDER__  'dr faustroll'
#endif

#ifdef __GFORTRAN__
#define __COMPILER__ 'gfortran'
#elif __INTEL__
#define __COMPILER__ 'ifort'
#elif CRAY
#define __COMPILER__ 'ftn'
#else
#define __COMPILER__ 'noidea'
#endif

#ifndef __VERSION__
#define __VERSION__ 'XYZ'
#endif

  Subroutine build_info()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 development subroutine for ending timing
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena & i.t.todorov april 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=10) :: time
    Character(Len=47) :: aux
    Character(Len=256) :: aux2
    Character(Len=5)  :: zone
    Character(Len=66) :: message
    Character(Len=8)  :: date
    Integer           :: i, value(1:8)

    Call info('', .true.)
    Call info(Repeat("*", 66), .true.)
    If (Len_trim(__DATE__//"  @  "//__TIME__) > 47) Then
      Write (aux, '(a47)') __DATE__//"  @  "//__TIME__
    Else
      Write (aux, *) __DATE__//"  @  "//__TIME__
    End If
    Call clean_string(aux)
    Write (message, '(a4,1x,a9,1x,a46,1x,a4)') "****", "birthday:", aux, "****"
    Call info(message, .true.)

    If (Len_trim(__HOSTNAME__) > 47) Then
      Write (aux, '(a47)') __HOSTNAME__
    Else
      Write (aux, *) __HOSTNAME__
    End If
    Call clean_string(aux)
    Write (message, '(a4,1x,a9,1x,a46,1x,a4)') "****", " machine:", aux, "****"
    Call info(message, .true.)

    If (Len_trim(__BUILDER__) > 47) Then
      Write (aux, '(a47)') __BUILDER__
    Else
      Write (aux, *) __BUILDER__
    End If
    Call clean_string(aux)
    Write (message, '(a4,1x,a9,1x,a46,1x,a4)') "****", " builder:", aux, "****"
    Call info(message, .true.)

    aux2 = compiler_version()
    Call clean_string(aux2)
    Do i = 1, Len_trim(aux2), 46
      aux = aux2(i:Min(i + 45, Len_trim(aux2)))
      Write (message, '(a4,1x,a9,1x,a,1x,a4)') "****", "compiler:", Trim(aux), "****"
      Call info(message, .true.)
    End Do

    If (mpi_ver > 0) Then
      Write (aux, '(a1,i0,a1,i0)') "v", mpi_ver, ".", mpi_subver
      Write (message, '(a4,1x,a9,1x,a46,1x,a4)') "****", "     MPI:", aux, "****"
      Call info(message, .true.)
#ifndef OLDMPI
      Call clean_string(lib_version)
      Do i = 1, Len_trim(lib_version), 46
        aux = lib_version(i:Min(i + 45, Len_trim(lib_version)))
        Write (message, '(a4,1x,a9,1x,a46,1x,a4)') "****", "MPI libs:", aux, "****"
        Call info(message, .true.)
      End Do
#endif
    Else If (mpi_ver < 0) Then
      Write (aux, *) "MPI Library too old.  Please update!!!"
      Write (message, '(a4,1x,a9,1x,a46,1x,a4)') "****", "MPI libs:", aux, "****"
      Call info(message, .true.)
    End If

    Call Date_and_time(date, time, zone, value)
    Write (aux, *) date(1:4), "-", date(5:6), "-", date(7:8), "  @  ", &
      time(1:2), ":", time(3:4), ":", time(5:10), "  (GMT", &
      zone(1:3), ":", zone(4:5), ")"
    Write (message, '(a4,1x,a9,a47,1x,a4)') "****", "executed:", aux, "****"
    Call info(message, .true.)
    Call info(Repeat("*", 66), .true.)
    Call info('', .true.)

  End Subroutine build_info


  Subroutine set_all_tests_true(this)
    Class(testing_type), Intent(inout) :: this
    this%configuration = .true.
    this%dftb_library = .true.
  End Subroutine set_all_tests_true

End Module development
