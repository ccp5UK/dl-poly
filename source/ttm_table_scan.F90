Subroutine ttm_table_scan(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning specific heat capacity,
! thermal conductivity and electron-phonon coupling
! constant table files to determine numbers of data points
!
! copyright - daresbury laboratory
! author    - m.a.seaton may 2012
! contrib   - g.khara    may 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms, Only : comms_type
  Use setup, Only : ntable,nrite
  Use parse, Only : get_line,get_word,word_2_real
  Use ttm

  Implicit None

   Type(comms_type), Intent(InOut) :: comm
  Logical                :: safe,lexist
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: fail
  Real( Kind = wp )      :: vk1,vk2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  If (l_ttm) Then

    fail=0
    Allocate (buffer(1:mxbuff), Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', comm%idnode
       Call error(0)
    End If

! check existence of thermal conductivity table file

    If (KeType == 3) Then

      Inquire (File='Ke.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 100
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')
      End If

! determine number of lines of data to read

      kel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 5
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word,comm)
          Call get_word(record,word)
          vk2 = word_2_real(word,comm)
          If (vk1>=zero_plus) kel=kel+1
        End If

      End Do
5  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (kel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(675)
      Else
        Allocate (ketable(1:kel,2), Stat=fail)
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', comm%idnode
          Call error(0)
        End If
        ketable(:,:) = 0.0_wp
      End If

    End If

! check existence of specific heat capacity table file

    If (CeType == 3) Then

      Inquire (File='Ce.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 200
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')
      End If

! determine number of lines of data to read

      cel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 10
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word,comm)
          Call get_word(record,word)
          vk2 = word_2_real(word,comm)
          If (vk1>=zero_plus) cel=cel+1
        End If

      End Do
10  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (cel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(677)
      Else
        Allocate (cetable(1:cel,2), Stat=fail)
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', comm%idnode
          Call error(0)
        End If
        cetable(:,:) = 0.0_wp
      End If

    End If

! check existence of thermal diffusivity table file

    If (DeType == 3) Then

      Inquire (File='De.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 300
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')
      End If

! determine number of lines of data to read

      del = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 15
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word,comm)
          Call get_word(record,word)
          vk2 = word_2_real(word,comm)
          If (vk1>=zero_plus) del=del+1
        End If

      End Do
15  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (del>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(679)
      Else
        Allocate (detable(1:del,2), Stat=fail)
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', comm%idnode
          Call error(0)
        End If
        detable(:,:) = 0.0_wp
      End If

    End If

! check existence of coupling constant table file

    If (gvar>0) Then

      Inquire (File='g.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 400
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')
      End If

! determine number of lines of data to read

      gel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 20
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word,comm)
          Call get_word(record,word)
          vk2 = word_2_real(word,comm)
          If (vk1>=zero_plus) gel=gel+1
        End If

      End Do
20  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (gel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(681)
      Else
        Allocate (gtable(1:gel,2), Stat=fail) ! [GK] array length corrected
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', comm%idnode
          Call error(0)
        End If
        gtable(:,:) = 0.0_wp
      End If

    End If

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

  End If

  Return

! end of Ke.dat file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(674)

! end of Ce.dat file error exit

200 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(676)

! end of g.dat file error exit

300 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(678)

400 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(680)

End Subroutine ttm_table_scan

