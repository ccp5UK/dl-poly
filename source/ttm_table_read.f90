Subroutine ttm_table_read()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading specific heat capacity and coupling
! constant table files
!
! copyright - daresbury laboratory
! author    - m.a.seaton may 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : ntable,nrite
  Use parse_module, Only : get_line,get_word,word_2_real
  Use ttm_module

  Implicit None

  Logical                :: safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: i
  Real( Kind = wp )      :: vk1,vk2

! read thermal conductivity data

  If (KeType == 3) Then

    If (idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')

    i = 0
    Do While(i<kel)

      Call get_line(safe,ntable,record)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          ketable(i,1) = vk1
          ketable(i,2) = vk2
        End If
      End If

    End Do

    If (idnode==0) Then
      Close(Unit=ntable)
      Write(nrite,'(/,1x,a)') 'thermal conductivity table read from Ke.dat file for two-temperature model'
      Write(nrite,'(1x,"minimum temperature            (K) = ",ES12.4,&
                 &/,1x,"maximum temperature            (K) = ",ES12.4,&
                 &/,1x,"minimum t.c. value   (W m^-1 K^-1) = ",ES12.4,&
                 &/,1x,"maximum t.c. value   (W m^-1 K^-1) = ",ES12.4)') &
                 Minval(ketable(:,1)),Maxval(ketable(:,1)),Minval(ketable(:,2)),Maxval(ketable(:,2))
    End If

! convert thermal conductivity values from W m^-1 K^-1 to kB A^-1 ps^-1

    ketable(1:kel,2) = ketable(1:kel,2)*JKms_to_kBAps

  End If

! read volumetric heat capacity data

  If (CeType == 3) Then

    If (idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')

    i = 0
    Do While(i<cel)

      Call get_line(safe,ntable,record)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          cetable(i,1) = vk1
          cetable(i,2) = vk2
        End If
      End If

    End Do

    If (idnode==0) Then
      Close(Unit=ntable)
      Write(nrite,'(/,1x,a)') 'electronic volumetric heat capacity table read from Ce.dat file for two-temperature model'
      Write(nrite,'(1x,"minimum temperature            (K) = ",ES12.4,&
                 &/,1x,"maximum temperature            (K) = ",ES12.4,&
                 &/,1x,"minimum v.h.c. value (J m^-3 K^-1) = ",ES12.4,&
                 &/,1x,"maximum v.h.c. value (J m^-3 K^-1) = ",ES12.4)') &
                 Minval(cetable(:,1)),Maxval(cetable(:,1)),Minval(cetable(:,2)),Maxval(cetable(:,2))
    End If

! convert volumetric heat capacity values from J m^-3 K^-1 to kB A^-3

    cetable(1:cel,2) = cetable(1:cel,2)*Jm3K_to_kBA3

  End If

! read thermal diffusivity data

  If (DeType == 3) Then

    If (idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')

    i = 0
    Do While(i<del)

      Call get_line(safe,ntable,record)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          detable(i,1) = vk1
          detable(i,2) = vk2
        End If
      End If

    End Do

    If (idnode==0) Then
      Close(Unit=ntable)
      Write(nrite,'(/,1x,a)') 'thermal diffusivity table read from De.dat file for two-temperature model'
      Write(nrite,'(1x,"minimum temperature            (K) = ",ES12.4,&
                 &/,1x,"maximum temperature            (K) = ",ES12.4,&
                 &/,1x,"minimum diffusivity value  (m^2/s) = ",ES12.4,&
                 &/,1x,"maximum diffusivity value  (m^2/s) = ",ES12.4)') &
                 Minval(detable(:,1)),Maxval(detable(:,1)),Minval(detable(:,2)),Maxval(detable(:,2))
    End If

! convert thermal diffusivity values from m^2 s^-1 to A^2 ps^-1

    detable(1:del,2) = detable(1:del,2)*1e8_wp

  End If

! read coupling constant data

  If (gvar>0) Then

    If (idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')

    i = 0
    Do While(i<gel)

      Call get_line(safe,ntable,record)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          gtable(i,1) = vk1
          gtable(i,2) = vk2
        End If
      End If

    End Do

    If (idnode==0) Then
      Close(Unit=ntable)
      Write(nrite,'(/,1x,a)') 'electron-phonon coupling table read from g.dat file for two-temperature model'
      Write(nrite,'(1x,"minimum temperature            (K) = ",ES12.4,&
                 &/,1x,"maximum temperature            (K) = ",ES12.4,&
                 &/,1x,"minimum e-p value    (W m^-3 K^-1) = ",ES12.4,&
                 &/,1x,"maximum e-p value    (W m^-3 K^-1) = ",ES12.4)') &
                 Minval(gtable(:,1)),Maxval(gtable(:,1)),Minval(gtable(:,2)),Maxval(gtable(:,2))
    End If

! convert electron-phonon coupling values from W m^-3 K^-1 to ps^-1

    gtable(1:gel,2) = gtable(1:gel,2)*epc_to_chi

  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(682)

End Subroutine ttm_table_read

