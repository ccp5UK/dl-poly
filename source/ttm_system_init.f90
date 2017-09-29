Subroutine ttm_system_init(nstep,nsteql,keyres,dumpfile,time,temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing electronic temperature restart files
! at job termination or selected intervals in simulation
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use setup_module
  Use ttm_module
  Use ttm_utils
  Use ttm_track_module
  Use comms_module
  Use parse_module,   Only : tabs_2_blanks, get_line, get_word, &
                             strip_blanks, word_2_real

  Implicit None

  Integer,             Intent ( In ) :: keyres,nstep,nsteql
  Real ( Kind = wp ),  Intent ( In ) :: temp,time
  Character (Len = *), Intent ( In ) :: dumpfile

  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Logical                :: safe,l_tmp = .true.
  Integer                :: nxx,nyy,nzz,nstp,i,ix,iy,iz,ii,jj,kk,ijk,ipos(3)
  Real ( Kind = wp )     :: eltmp,tme,lat_sum,lat_max,lat_min
  Integer                :: iounit = 225

! check existence of readable restart file (DUMP_E)

  If (idnode == 0) Inquire(File=dumpfile, Exist=l_tmp)
  If (mxnode > 1) Call gcheck(l_tmp)
  If ((.not. l_tmp) .and. keyres==keyres0) Call error(684)

! if restarting simulation, read restart file

  If (l_tmp .and. keyres==keyres0) Then

    If (idnode==0) Open (Unit=iounit, File=dumpfile)
    Call get_line(safe,iounit,record); If (.not.safe) Goto 100
    Call get_word(record,word) ; nxx=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; nyy=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; nzz=Nint(word_2_real(word,0.0_wp))
    Call get_line(safe,iounit,record); If (.not.safe) Goto 100
    Call get_word(record,word) ; nstp=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; tme=word_2_real(word,0.0_wp)
    Call get_word(record,word) ; depostart=word_2_real(word,0.0_wp)
    Call get_word(record,word) ; depoend=word_2_real(word,0.0_wp)
    ! check size of electronic temperature grid matches with size given in CONTROL file
    If (nxx/=eltsys(1) .or. nyy/=eltsys(2) .or. nzz/=eltsys(3)) Call error(685)
    ! check restart file is at same timestep as restart
    ! (can proceed if not, but need to warn user)
    If (nstp/=nstep .or. Abs(tme-time)>zero_plus) Call warning(520,0.0_wp,0.0_wp,0.0_wp)
    ! read in each line, find appropriate grid cell and assign
    ! electronic temperature if processor has that cell
    Do i=1,eltsys(1)*eltsys(2)*eltsys(3)
      Call get_line(safe,iounit,record); If (.not.safe) Goto 100
      Call get_word(record,word) ; ipos(1)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; ipos(2)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; ipos(3)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; eltmp=word_2_real(word,0.0_wp)
      ix = ipos(1) + midI(1) - 1
      iy = ipos(2) + midI(2) - 1
      iz = ipos(3) + midI(3) - 1
      ii = Floor(Real(ix,Kind=wp)/Real(ntsys(1),Kind=wp))
      jj = Floor(Real(iy,Kind=wp)/Real(ntsys(2),Kind=wp))
      kk = Floor(Real(iz,Kind=wp)/Real(ntsys(3),Kind=wp))
      ix = Mod(ix+ntsys(1)*eltcell(1),ntsys(1)) + 1 - ntcelloff(1)
      iy = Mod(iy+ntsys(2)*eltcell(2),ntsys(2)) + 1 - ntcelloff(2)
      iz = Mod(iz+ntsys(3)*eltcell(3),ntsys(3)) + 1 - ntcelloff(3)
      If (ix>0 .and. ix<=ntcell(1) .and. iy>0 .and. iy<=ntcell(2) .and. iz>0 .and. iz<=ntcell(3)) Then
        ijk = 1 + ix + (ntcell(1)+2) * (iy + (ntcell(2)+2) * iz)
        eltemp(ijk,ii,jj,kk) = eltmp
      End If
    End Do
    ! fill boundary halo values and deal with required boundary conditions
    Call boundaryHalo ()
    Call boundaryCond (bcTypeE, temp)
    ! check whether or not energy deposition has happened yet
    If(nstep>nsteql .and. time>=depostart .and. time<depoend) Then
      Call depoinit(time)
    Else If (time>=depoend) Then
      findepo = .true.
    End If

    ! report successful reading and minimum, maximum and sums of
    ! electronic temperatures
    Call eltemp_sum (lat_sum)
    Call eltemp_max (lat_max)
    Call eltemp_min (lat_min)
    If (idnode==0) Then
      Write(nrite,'(/,1x,a)') 'electronic temperatures read from DUMP_E file for two-temperature model'
      Write(nrite,'(1x,"minimum temperature (K) = ",ES8.4,&
                 &/,1x,"maximum temperature (K) = ",ES8.4,&
                 &/,1x,"sum of temperatures (K) = ",ES8.4)') &
                 lat_min, lat_max, lat_sum
      Close (iounit)
    End If

  Else

! if not restarting simulation, set electronic temperature grid
! to system temperature

    eltemp = temp

  End If

  Return

! Abnormal exit from electronic temperature dump file read

100 Continue

  If (idnode == 0) Write(nrite,"(/,1x,a)") dumpfile, ' data mishmash detected'
  Call error(686)
  Return

End Subroutine ttm_system_init

