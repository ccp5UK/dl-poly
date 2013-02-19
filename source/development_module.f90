Module development_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development module
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical, Save :: l_scr  = .false. ! OUTPUT redirection to the default output (screen)
  Logical, Save :: l_eng  = .false. ! OUTPUT inclusion of an extra last line with E_tot
  Logical, Save :: l_rout = .false. ! REVIVE writing in ASCII (default is binary)
  Logical, Save :: l_rin  = .false. ! REVOLD reading in ASCII (default is binary)
  Logical, Save :: l_org  = .false. ! translate CONFIG along a vector to CFGORG
  Logical, Save :: l_scl  = .false. ! CONFIG rescaling to CFGSCL after reading with termination
  Logical, Save :: l_his  = .false. ! HISTORY generation after reading with termination
  Logical, Save :: l_trm  = .false. ! termination flag
  Logical, Save :: l_tim  = .false. ! detailed timing
  Logical, Save :: l_tor  = .false. ! no production of REVCON & REVIVE
  Logical, Save :: l_dis  = .false. ! check on minimum separation distance between VNL pairs at re/start


  Integer, Save           :: lvcforg = -1                                ! CFGORG levcfg
  Real( Kind = wp ), Save :: xorg = 0.0_wp, yorg = 0.0_wp, zorg = 0.0_wp ! reorigin vector

  Integer, Save           :: lvcfscl = -1       ! CFGSCL levcfg
  Real( Kind = wp ), Save :: cels(1:9) = 0.0_wp ! CFGSCL lattice parameters

  Real( Kind = wp ), Save :: r_dis = 0.5_wp ! l_dis default check condition

  Real( Kind = wp ), Save :: t_zero

Contains

  Subroutine scan_development()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 development subroutine for raw scanning the contents of the
! control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use comms_module, Only : idnode,mxnode,gcheck
    Use setup_module, Only : nread
    Use parse_module, Only : get_line,get_word,lower_case

    Implicit None

    Logical                :: carry,safe
    Character( Len = 200 ) :: record
    Character( Len = 40  ) :: word

! Set safe flag

    safe=.true.

! Open the simulation input file

    If (idnode == 0) Inquire(File='CONTROL', Exist=safe)
    If (mxnode > 1) Call gcheck(safe)
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
    Use setup_module, Only : nrite

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

End Module development_module
