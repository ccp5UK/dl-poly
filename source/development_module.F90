Module development_module

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
  Use setup_module, Only : nrite, nread, control
#ifdef OLDMPI
  Use comms_module, Only : mpi_ver,mpi_subver
#else
  Use comms_module, Only : mpi_ver,mpi_subver,lib_version
#endif

  Implicit None

  Private

  Logical, Save, Public :: l_scr  = .false. ! OUTPUT redirection to the default output (screen)
  Logical, Save, Public :: l_fast = .false. ! avoid global safety checks (no elegant parallel failures)
  Logical, Save, Public :: l_eng  = .false. ! OUTPUT inclusion of an extra last line with E_tot
  Logical, Save, Public :: l_rout = .false. ! REVIVE writing in ASCII (default is binary)
  Logical, Save, Public :: l_rin  = .false. ! REVOLD reading in ASCII (default is binary)
  Logical, Save, Public :: l_org  = .false. ! translate CONFIG along a vector to CFGORG
  Logical, Save, Public :: l_scl  = .false. ! CONFIG rescaling to CFGSCL after reading with termination
  Logical, Save, Public :: l_his  = .false. ! HISTORY generation after reading with termination
  Logical, Save, Public :: l_trm  = .false. ! termination flag
  Logical, Save         :: l_tim  = .false. ! detailed timing
  Logical, Save, Public :: l_tor  = .false. ! no production of REVCON & REVIVE
  Logical, Save, Public :: l_dis  = .false. ! check on minimum separation distance between VNL pairs at re/start


  Integer, Save, Public           :: lvcforg = -1                                ! CFGORG levcfg
  Real( Kind = wp ), Save, Public :: xorg = 0.0_wp, yorg = 0.0_wp, zorg = 0.0_wp ! reorigin vector

  Integer, Save, Public           :: lvcfscl = -1       ! CFGSCL levcfg
  Real( Kind = wp ), Save, Public :: cels(1:9) = 0.0_wp ! CFGSCL lattice parameters

  Real( Kind = wp ), Save, Public :: r_dis = 0.5_wp ! l_dis default check condition

  Real( Kind = wp ), Save, Public :: t_zero

  Public :: scan_development
  Public :: build_info

Contains

  Subroutine scan_development()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development subroutine for raw scanning the contents of the
! control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use comms_module, Only : idnode,mxnode,gcheck
    Use parse_module, Only : get_line,get_word,lower_case

    Implicit None

    Logical                :: carry,safe
    Character( Len = 200 ) :: record
    Character( Len = 40  ) :: word

! Set safe flag

    safe=.true.

! Open the simulation input file

    If (idnode == 0) Inquire(File=Trim(control), Exist=safe)
    If (mxnode > 1) Call gcheck(safe,"enforce")
    If (.not.safe) Then
       Return
    Else
       If (idnode == 0) Open(Unit=nread, File=Trim(control), Status='old')
    End If

    Call get_line(safe,nread,record)
    If (.not.safe) Go To 10

    carry = .true.
    Do While (carry)

       Call get_line(safe,nread,record)
       If (.not.safe) Go To 10

       Call lower_case(record)
       Call get_word(record,word)

! read DEVELOPMENT option: OUTPUT to screen

       If      (word(1:5) == 'l_scr') Then

          l_scr = .true.

       Else If (word(1:6) == 'l_fast') Then

          l_fast=.true.

       Else If (word(1:5) == 'l_tim') Then

          l_tim=.true.

       Else If (word(1:6) == 'finish') Then

          carry=.false.

       End If

    End Do

10  Continue
    If (idnode == 0) Close(Unit=nread)

  End Subroutine scan_development

  Subroutine start_devel_time(name)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development subroutine for starting timing
!
! copyright - daresbury laboratory
! author    - i.j.bush november 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use comms_module, Only : gtime, gsync

    Character( Len = * ), Intent( In    ) :: name

    If (l_tim) Then
       Call gsync()
       Call gtime(t_zero)
    End If

  End Subroutine start_devel_time

  Subroutine end_devel_time(name)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development subroutine for ending timing
!
! copyright - daresbury laboratory
! author    - i.j.bush november 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use comms_module, Only : gtime, gsync, idnode

    Character( Len = * ), Intent( In    ) :: name

    Real( Kind = wp ) :: t

    If (l_tim) Then
       Call gsync()
       Call gtime(t)

       If (idnode == 0) Then
          Write(nrite,'(1x,2(a,3x),f0.3)') 'DEVEL TIME: Time in', name, t-t_zero
       End If

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

  Subroutine build_info()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development subroutine for ending timing
!
! copyright - daresbury laboratory
! author    - a.m.elena & i.t.todorov april 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use parse_module, Only : clean_string

    Implicit None

    Character( Len =  8 ) :: date
    Character( Len = 10 ) :: time
    Character( Len =  5 ) :: zone
    Integer               :: value(1:8)

    Character( Len = 47 ) :: aux
    Integer               :: i,l

    Write(nrite,'(1x,a66)') Repeat("*",66)
    If (Len_Trim( __DATE__//"  @  "//__TIME__) > 47) Then
      Write(aux,'(a47)') __DATE__//"  @  "//__TIME__
    Else
      Write(aux,*) __DATE__//"  @  "//__TIME__
    End If
    Call clean_string(aux)
    Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", "birthday:", aux, "****"
    If (Len_Trim(__HOSTNAME__) > 47) Then
      Write(aux,'(a47)') __HOSTNAME__
    Else
      Write(aux,*) __HOSTNAME__
    End If
    Call clean_string(aux)
    Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", " machine:", aux, "****"
    If (Len_Trim(__BUILDER__) > 47) Then
      Write(aux,'(a47)') __BUILDER__
    Else
      Write(aux,*) __BUILDER__
    End If
    Call clean_string(aux)
    Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", " builder:", aux, "****"
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
    Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", "compiler:", aux, "****"
    If (mpi_ver > 0) Then
       Write(aux,'(a1,i0,a1,i0)') "v",mpi_ver,".",mpi_subver
       Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", "     MPI:", aux, "****"
#ifndef OLDMPI
       Call clean_string(lib_version)
       Do i=1,Len_Trim(lib_version),46
          aux=lib_version(i:Min(i+45,Len_Trim(lib_version)))
          Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", "MPI libs:", aux, "****"
       End Do
#endif
    Else If (mpi_ver < 0) Then
       Write(aux,*) "MPI Library too old.  Please update!!!"
       Write(nrite,'(1x,a4,1x,a9,1x,a46,1x,a4)') "****", "MPI libs:", aux, "****"
    End If
    Call date_and_time(date,time,zone,value)
    Write(aux,*) date(1:4),"-",date(5:6),"-",date(7:8),"  @  ",   &
                 time(1:2),":",time(3:4),":",time(5:10),"  (GMT", &
                 zone(1:3),":",zone(4:5),")"
    Write(nrite,'(1x,a4,1x,a9,a47,1x,a4)') "****", "executed:", aux, "****"
    Write(nrite,'(1x,a66,/)') Repeat("*",66)

  End Subroutine build_info

End Module development_module
