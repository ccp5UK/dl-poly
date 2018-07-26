Subroutine ttm_table_scan()

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

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : ntable,nrite
  Use parse_module, Only : get_line,get_word,word_2_real
  Use ttm_module

  Implicit None

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
       Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', idnode
       Call error(0)
    End If

! check existence of thermal conductivity table file

    If (KeType == 3) Then

      Inquire (File='Ke.dat', Exist=lexist)
      If (mxnode > 1) Call gcheck(lexist)

      If (.not.lexist) Then
        Go To 100
      Else
        If (idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')
      End If

! determine number of lines of data to read

      kel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record)
        If (.not.safe) Then
          Go To 5
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) kel=kel+1
        End If

      End Do
5  Continue

      If (idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (kel>0)
      If (mxnode > 1) Call gcheck(safe)
      If (.not. safe) Then
        Call error(685)
      Else
        Allocate (ketable(1:kel,2), Stat=fail)
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', idnode
          Call error(0)
        End If
        ketable(:,:) = 0.0_wp
      End If

    End If

! check existence of specific heat capacity table file

    If (CeType == 3) Then

      Inquire (File='Ce.dat', Exist=lexist)
      If (mxnode > 1) Call gcheck(lexist)

      If (.not.lexist) Then
        Go To 200
      Else
        If (idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')
      End If

! determine number of lines of data to read

      cel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record)
        If (.not.safe) Then
          Go To 10
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) cel=cel+1
        End If

      End Do
10  Continue

      If (idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (cel>0)
      If (mxnode > 1) Call gcheck(safe)
      If (.not. safe) Then
        Call error(687)
      Else
        Allocate (cetable(1:cel,2), Stat=fail)
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', idnode
          Call error(0)
        End If
        cetable(:,:) = 0.0_wp
      End If

    End If

! check existence of thermal diffusivity table file

    If (DeType == 3) Then

      Inquire (File='De.dat', Exist=lexist)
      If (mxnode > 1) Call gcheck(lexist)

      If (.not.lexist) Then
        Go To 300
      Else
        If (idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')
      End If

! determine number of lines of data to read

      del = 0
      Do While(.true.)

        Call get_line(safe,ntable,record)
        If (.not.safe) Then
          Go To 15
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) del=del+1
        End If

      End Do
15  Continue

      If (idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (del>0)
      If (mxnode > 1) Call gcheck(safe)
      If (.not. safe) Then
        Call error(689)
      Else
        Allocate (detable(1:del,2), Stat=fail)
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', idnode
          Call error(0)
        End If
        detable(:,:) = 0.0_wp
      End If

    End If

! check existence of coupling constant table file

    If (gvar>0) Then

      Inquire (File='g.dat', Exist=lexist)
      If (mxnode > 1) Call gcheck(lexist)

      If (.not.lexist) Then
        Go To 400
      Else
        If (idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')
      End If

! determine number of lines of data to read

      gel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record)
        If (.not.safe) Then
          Go To 20
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) gel=gel+1
        End If

      End Do
20  Continue

      If (idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (gel>0)
      If (mxnode > 1) Call gcheck(safe)
      If (.not. safe) Then
        Call error(691)
      Else
        Allocate (gtable(1:gel,2), Stat=fail) ! [GK] array length corrected
        If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan allocation failure, node: ', idnode
          Call error(0)
        End If
        gtable(:,:) = 0.0_wp
      End If

    End If

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'ttm_table_scan deallocation failure, node: ', idnode
       Call error(0)
    End If

  End If

  Return

! end of Ke.dat file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(684)

! end of Ce.dat file error exit

200 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(686)

! end of g.dat file error exit

300 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(688)

400 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(690)

End Subroutine ttm_table_scan

