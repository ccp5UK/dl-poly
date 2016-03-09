Module development_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development module
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
! contrib   - i.j.bush november 2008
! contrib   - a.m.elena march 2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : nrite, nread
  Use comms_module, Only : mpi_ver,mpi_subver,lib_version

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
  Logical, Save :: l_tim  = .false. ! detailed timing
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

    If (idnode == 0) Inquire(File='CONTROL', Exist=safe)
    If (mxnode > 1) Call gcheck(safe,"enforce")
    If (.not.safe) Then
       Return
    Else
       If (idnode == 0) Open(Unit=nread, File='CONTROL', Status='old')
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

  Subroutine build_info()
    Character(len=48) :: aux

    Integer :: i

    Write(nrite,'(1x,a66)')Repeat("*",66)
    Write(aux,*)__DATE__//"@"//__TIME__
    Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***","birthday: ",aux,"****"
    Write(aux,*)__HOSTNAME__
    Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***","machine: ",aux,"****"
    Write(aux,*)__BUILDER__
    Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***","builder: ",aux,"****"
    Write(aux,*)__COMPILER__
    Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***","compiler: ",aux,"****"
    Write(aux,*)__VERSION__
    Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***","version: ",aux,"****"
    If (mpi_ver > 0) Then
      Write(nrite,'(1x,a3,1x,a10,1x,i1,a1,i1,a48)')"***","MPI: ",mpi_ver,".",mpi_subver,Repeat(" ",44)//"****"
      Do i=1,Len_trim(lib_version),47
        aux=lib_version(i:Min(i+46,len_trim(lib_version)))
        Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***","library: ",aux,"****"
      End Do
    Else IF (mpi_ver < 0) Then
      Write(aux,*)"MPI Library too old. Update!"
      Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***"," ",aux,"****"
    Else
      Write(aux,*)"Serial mode selected"
      Write(nrite,'(1x,a3,1x,a10,a48,a4)')"***"," ",aux,"****"
    End If
    Write(nrite,'(1x,a66,/)')Repeat("*",66)
  End Subroutine build_info
End Module development_module
