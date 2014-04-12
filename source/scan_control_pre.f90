Subroutine scan_control_pre(imc_n,dvar)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning the imc_n & dvar options in the
! control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module,  Only : nread
  Use parse_module,  Only : get_line,get_word,lower_case,word_2_real
  Implicit None

  Integer,           Intent(   Out ) :: imc_n
  Real( Kind = wp ), Intent(   Out ) :: dvar

  Logical                :: carry,safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word

  safe   = .true.  ! all is safe

! density variation parameter default

  dvar = 1.0_wp

  If (idnode == 0) Inquire(File='CONTROL', Exist=safe)
  If (mxnode > 1) Call gcheck(safe,"enforce")
  If (.not.safe) Then
     Go To 10
  Else
     If (idnode == 0) Open(Unit=nread, File='CONTROL', Status='old')
  End If

! Read TITLE record

  Call get_line(safe,nread,record)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! read density variation option
! this is really a pre-scan in order to get the MD box dimensions
! from scan_config giving it failure estimates for when reading

     If      (word(1:7) == 'densvar') Then

        Call get_word(record,word)
        dvar = Abs(word_2_real(word))
        dvar = 1.0_wp + Abs(dvar)/100.0_wp

! read slab option
! limiting DD slicing in z direction to 2 for load balancing purposes
! this is really a pre-scan in order to get the MD box dimensions
! from scan_config before the option is read again in scan_control

     Else If (word(1:4) == 'slab') Then

        imc_n=6

! io options

! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     End If

  End Do

  If (idnode == 0) Close(Unit=nread)

  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control_pre
