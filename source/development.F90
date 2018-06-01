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
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup, Only : nread, control
#ifdef OLDMPI
  Use comms, Only : mpi_ver,mpi_subver, comms_type,gcheck, gtime, gsync
#else
  Use comms, Only : mpi_ver,mpi_subver,lib_version,comms_type,gcheck, &
    gtime, gsync
#endif
  Use parse, Only : get_line,get_word,lower_case,clean_string
  Use errors_warnings, Only : info

  Implicit None

  Private

  !> Type containing development module variables
  Type, Public :: development_type
    Private

    !> OUTPUT redirection to the default output (screen)
    Logical, Public :: l_scr  = .false.
    !> avoid global safety checks (no elegant parallel failures)
    Logical, Public :: l_fast = .false.
    !> OUTPUT inclusion of an extra last line with E_tot
    Logical, Public :: l_eng  = .false.
    !> REVIVE writing in ASCII (default is binary)
    Logical, Public :: l_rout = .false.
    !> REVOLD reading in ASCII (default is binary)
    Logical, Public :: l_rin  = .false.
    !> translate CONFIG along a vector to CFGORG
    Logical, Public :: l_org  = .false.
    !> CONFIG rescaling to CFGSCL after reading with termination
    Logical, Public :: l_scl  = .false.
    !> HISTORY generation after reading with termination
    Logical, Public :: l_his  = .false.
    !> termination flag
    Logical, Public :: l_trm  = .false.
    !> detailed timing
    Logical         :: l_tim  = .false.
    !> no production of REVCON & REVIVE
    Logical, Public :: l_tor  = .false.
    !> check on minimum separation distance between VNL pairs at re/start
    Logical, Public :: l_dis  = .false.


    !> CFGORG levcfg
    Integer, Public           :: lvcforg = -1
    !> reorigin vector
    Real( Kind = wp ), Public :: xorg = 0.0_wp, yorg = 0.0_wp, zorg = 0.0_wp

    !> CFGSCL levcfg
    Integer, Public           :: lvcfscl = -1
    !> CFGSCL lattice parameters
    Real( Kind = wp ), Public :: cels(1:9) = 0.0_wp

    !> l_dis default check condition
    Real( Kind = wp ), Public :: r_dis = 0.5_wp

    !> Devel start time
    Real( Kind = wp ), Public :: t_zero
  End Type development_type

  Public :: scan_development
  Public :: build_info

Contains

  Subroutine scan_development(devel,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 development subroutine for raw scanning the contents of the
    ! control file
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( development_type ), Intent( InOut ) :: devel
    Type( comms_type ), Intent( InOut ) :: comm
    Logical                :: carry,safe
    Character( Len = 200 ) :: record
    Character( Len = 40  ) :: word

    ! Set safe flag

    safe=.true.

    ! Open the simulation input file

    If (comm%idnode == 0) Inquire(File=Trim(control), Exist=safe)
    If (comm%mxnode > 1) Call gcheck(comm,safe,"enforce")
    If (.not.safe) Then
      Return
    Else
      If (comm%idnode == 0) Open(Unit=nread, File=Trim(control), Status='old')
    End If

    Call get_line(safe,nread,record,comm)
    If (safe) Then

      carry = .true.
      Do While (carry)

        Call get_line(safe,nread,record,comm)
        If (.not.safe) Exit

        Call lower_case(record)
        Call get_word(record,word)

        ! read DEVELOPMENT option: OUTPUT to screen

        If      (word(1:5) == 'l_scr') Then

          devel%l_scr = .true.

        Else If (word(1:6) == 'l_fast') Then

          devel%l_fast=.true.

        Else If (word(1:5) == 'l_tim') Then

          devel%l_tim=.true.

        Else If (word(1:6) == 'finish') Then

          carry=.false.

        End If

      End Do
    End If
    If (comm%idnode == 0) Close(Unit=nread)

  End Subroutine scan_development

  Subroutine start_devel_time(name,devel,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 development subroutine for starting timing
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush november 2009
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character( Len = * ), Intent( In    ) :: name
    Type( development_type ), Intent( InOut ) :: devel
    Type( comms_type ), Intent( InOut ) :: comm

    If (devel%l_tim) Then
      Call gsync(comm)
      Call gtime(devel%t_zero)
    End If

  End Subroutine start_devel_time

  Subroutine end_devel_time(name,devel,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 development subroutine for ending timing
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush november 2009
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character( Len = * ), Intent( In    ) :: name
    Type( development_type ), Intent( InOut ) :: devel
    Type( comms_type ), Intent( InOut ) :: comm

    Real( Kind = wp ) :: t
    Character( Len = 256 ) :: message

    If (devel%l_tim) Then
      Call gsync(comm)
      Call gtime(t)

      Write(message,'(2(a,3x),f0.3)') 'DEVEL TIME: Time in', name, t-devel%t_zero
      Call info(message,.true.)

    End If

  End Subroutine end_devel_time

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

  Subroutine build_info(devel)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 development subroutine for ending timing
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena & i.t.todorov april 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( development_type ), Intent( In    ) :: devel

    Character( Len =  8 ) :: date
    Character( Len = 10 ) :: time
    Character( Len =  5 ) :: zone
    Integer               :: value(1:8)

    Character( Len = 47 ) :: aux
    Character( Len = 66 ) :: message
    Integer               :: i,l

    Call info('',.true.)
    Call info(Repeat("*",66),.true.)
    If (Len_Trim( __DATE__//"  @  "//__TIME__) > 47) Then
      Write(aux,'(a47)') __DATE__//"  @  "//__TIME__
    Else
      Write(aux,*) __DATE__//"  @  "//__TIME__
    End If
    Call clean_string(aux)
    Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", "birthday:", aux, "****"
    Call info(message,.true.)

    If (Len_Trim(__HOSTNAME__) > 47) Then
      Write(aux,'(a47)') __HOSTNAME__
    Else
      Write(aux,*) __HOSTNAME__
    End If
    Call clean_string(aux)
    Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", " machine:", aux, "****"
    Call info(message,.true.)

    If (Len_Trim(__BUILDER__) > 47) Then
      Write(aux,'(a47)') __BUILDER__
    Else
      Write(aux,*) __BUILDER__
    End If
    Call clean_string(aux)
    Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", " builder:", aux, "****"
    Call info(message,.true.)

    If      (mpi_ver == 0) Then
      If (Len_Trim(__COMPILER__//" v"//__VERSION__//" (serial build)") > 47) Then
        Write(aux,'(a47)') __COMPILER__//" v"//__VERSION__//" (serial build)"
      Else
        Write(aux,*) __COMPILER__//" v"//__VERSION__//" (serial build)"
      End If
    Else If (mpi_ver >  0) Then
      If (Len_Trim(__COMPILER__//" v"//__VERSION__) > 47) Then
        Write(aux,'(a47)') __COMPILER__//" v"//__VERSION__
      Else
        Write(aux,*) __COMPILER__//" v"//__VERSION__
      End If
    End If
    Call clean_string(aux)
    Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", "compiler:", aux, "****"
    Call info(message,.true.)

    If (mpi_ver > 0) Then
      Write(aux,'(a1,i0,a1,i0)') "v",mpi_ver,".",mpi_subver
      Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", "     MPI:", aux, "****"
      Call info(message,.true.)
#ifndef OLDMPI
      Call clean_string(lib_version)
      Do i=1,Len_Trim(lib_version),46
        aux=lib_version(i:Min(i+45,Len_Trim(lib_version)))
        Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", "MPI libs:", aux, "****"
        Call info(message,.true.)
      End Do
#endif
    Else If (mpi_ver < 0) Then
      Write(aux,*) "MPI Library too old.  Please update!!!"
      Write(message,'(a4,1x,a9,1x,a46,1x,a4)') "****", "MPI libs:", aux, "****"
      Call info(message,.true.)
    End If

    Call date_and_time(date,time,zone,value)
    Write(aux,*) date(1:4),"-",date(5:6),"-",date(7:8),"  @  ",   &
      time(1:2),":",time(3:4),":",time(5:10),"  (GMT", &
      zone(1:3),":",zone(4:5),")"
    Write(message,'(a4,1x,a9,a47,1x,a4)') "****", "executed:", aux, "****"
    Call info(message,.true.)
    Call info(Repeat("*",66),.true.)
    Call info('',.true.)

  End Subroutine build_info

End Module development