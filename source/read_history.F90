Subroutine read_history(l_str,fname,megatm,levcfg,dvar,nstep,tstep,time,exout)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading the trajectory data file
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, li
  Use comms_module
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz,nprx_r,npry_r,nprz_r
  Use site_module
  Use config_module
  Use parse_module,   Only : tabs_2_blanks, get_line, get_word, &
                             strip_blanks, word_2_real
  Use io_module,      Only : io_set_parameters,     &
                             io_get_parameters,     &
                             io_init, io_open,      &
                             io_read_batch,         &
                             io_nc_get_dim,         &
                             io_nc_get_var,         &
                             io_close, io_finalize, &
                             IO_READ_MPIIO,         &
                             IO_READ_DIRECT,        &
                             IO_READ_NETCDF,        &
                             IO_READ_MASTER

  Implicit None

  Character( Len = * ), Intent( In    ) :: fname
  Logical,              Intent( In    ) :: l_str
  Integer,              Intent( In    ) :: megatm

  Integer,              Intent( InOut ) :: levcfg,nstep
  Real( Kind = wp ),    Intent( In    ) :: dvar
  Real( Kind = wp ),    Intent( InOut ) :: tstep,time
  Integer,              Intent(   Out ) :: exout

  Logical,               Save :: newjob = .true.  , &
                                 l_ind  = .true.  , &
                                 l_his  = .true.  , &
                                 l_xtr  = .false. , &
                                 fast   = .true.
  Integer(Kind=li),      Save :: rec  = 0_li , &
                                 frm  = 0_li , &
                                 frm1 = 0_li

  Logical                :: safe, lexist
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: fail(1:5),i,k,l,m,idm,ipx,ipy,ipz, &
                            indatm,nattot,totatm
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz,xhi,yhi,zhi, &
                            cell_vecs(1:3,1:3)

! Some parameters and variables needed by io_module interfaces

  Integer,                           Save :: recsz ! record size
  Integer,                           Save :: fh, io_read
  Integer( Kind = MPI_OFFSET_KIND ), Save :: top_skip
  Character( Len = 1),  Allocatable, Save :: buffer(:,:)

  Character( Len = 8 ), Allocatable :: chbuf(:)
  Integer,              Allocatable :: iwrk(:)
  Real( Kind = wp ),    Allocatable :: axx(:),ayy(:),azz(:), &
                                       bxx(:),byy(:),bzz(:), &
                                       cxx(:),cyy(:),czz(:)

  If (newjob) Then
     newjob = .false.

! Get type of I/O

     Call io_get_parameters( user_method_read = io_read )

! ASCII read

     If (io_read /= IO_READ_NETCDF) Then

! Define the default record size

        recsz = 73

! Does HISTORY exist

        lexist=.true.
        If (idnode == 0) Inquire(File=fname, Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist,"enforce")
        If (.not.lexist) Go To 400

! Open HISTORY

        If (idnode == 0) Open(Unit=nconf, File=fname)

! read the HISTORY file header

        Call get_line(safe,nconf,record); If (.not.safe) Go To 300

        Call get_line(safe,nconf,record); If (.not.safe) Go To 300
        Call get_word(record,word) ; levcfg=Nint(word_2_real(word,0.0_wp))
        Call get_word(record,word)
        Call get_word(record,word) ; If (Nint(word_2_real(word)) /= megatm) Go To 300
        Call get_word(record,word) ; frm=Nint(word_2_real(word,0.0_wp),li)
        Call get_word(record,word) ; rec=Nint(word_2_real(word,0.0_wp),li)

! Change fast if no database records exists or the file is new

        If (frm == Int(0,li) .and. rec == Int(0,li)) fast = .false.

        If (levcfg /= 3) Then

! Detect DL_POLY_2 HISTORY and amend io_read and recsz

           If ((.not.fast) .and. io_read == IO_READ_MPIIO) Then
              io_read = IO_READ_DIRECT
              recsz   = 200
           End If

        Else

! Detect compressed HISTORY and make amendments

           l_ind = .false.
           If (fast .and. io_read /= IO_READ_MASTER) Then
              io_read = IO_READ_MPIIO
              recsz   = 35
           Else
              io_read = IO_READ_DIRECT
           End If

        End If

        If (io_read /= IO_READ_MASTER) Then
           top_skip = Int(2,MPI_OFFSET_KIND)

           fail(1) = 0
           Allocate (buffer(1:recsz,1:4), Stat=fail(1))
           If (fail(1) > 0) Then
              Write(nrite,'(/,1x,a,i0)') 'read_history allocation failure 1, node: ', idnode
              Call error(0)
           End If

           If (io_read == IO_READ_MPIIO) Then
              Close(Unit=nconf)

              Call io_set_parameters( user_comm = dlp_comm_world )
              Call io_init( recsz )
              Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )
           End If
        End If

     Else ! netCDF read

! Does HISTORY exist

        lexist=.true.
        If (idnode == 0) Inquire(File=fname, Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist,"enforce")
        If (.not.lexist) Go To 400

! fast and rec are irrelevant for netCDF (initialised at declaration)

        fast = .true.
        rec  = Int(0,li)

        Call io_set_parameters( user_comm = dlp_comm_world )
        Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )
        Call io_nc_get_dim( 'frame', fh, i )

        If (i > 0) Then
           frm=Int(i,li)
        Else
           Go To 300
        End If

     End If
  Else
     If (io_read == IO_READ_MPIIO .or. io_read == IO_READ_NETCDF) &
        Call io_set_parameters( user_comm = dlp_comm_world )
     If (io_read == IO_READ_MPIIO) Call io_init( recsz )
  End If

! Reinitialise local-to-global counters and all levels of information at every bloody read
! Ilian's & Alin's fix to clear incoherencies with re-emptied/re-filled domains

  lsi=0 ; lsa=0 ; ltg=0                     ! A must unfortunately

! Necessary to not initialise keep the last frame info after hitting EOF
!  xxx=0.0_wp ; yyy = 0.0_wp ; zzz = 0.0_wp ! unfortunate but
!  vxx=0.0_wp ; vyy = 0.0_wp ; vzz = 0.0_wp ! unfortunate but
!  fxx=0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp ! unfortunate but

! MASTER READ

  If      (io_read == IO_READ_MASTER) Then

     fail=0
     Allocate (chbuf(1:mxatms),                           Stat=fail(1))
     Allocate (iwrk(1:mxatms),                            Stat=fail(2))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(3))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(4))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_history allocation failure, node: ', idnode
        Call error(0)
     End If

! read timestep and time

     Call get_line(safe,nconf,record); If (.not.safe) Go To 200

     xxx=0.0_wp ; yyy = 0.0_wp ; zzz = 0.0_wp
     vxx=0.0_wp ; vyy = 0.0_wp ; vzz = 0.0_wp
     fxx=0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp

     Call get_word(record,word) ! timestep
     Call get_word(record,word) ; nstep = Nint(word_2_real(word))
     Call get_word(record,word) ; If (Nint(word_2_real(word)) /= megatm) Go To 300
     Call get_word(record,word) ; levcfg = Nint(word_2_real(word))
     Call get_word(record,word) ; imcon = Nint(word_2_real(word))

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

     Call get_word(record,word) ; tstep = word_2_real(word)
     Call get_word(record,word) ; time = word_2_real(word)

     If (idnode == 0) Write(nrite,"(/,1x,'HISTORY step',i10,' (',f10.3,' ps) is being read')") nstep,time

! read cell vectors

     Call get_line(safe,nconf,record); If (.not.safe) Go To 300
     Call get_word(record,word); cell(1)=word_2_real(word)
     Call get_word(record,word); cell(2)=word_2_real(word)
     Call get_word(record,word); cell(3)=word_2_real(word)

     Call get_line(safe,nconf,record); If (.not.safe) Go To 300
     Call get_word(record,word); cell(4)=word_2_real(word)
     Call get_word(record,word); cell(5)=word_2_real(word)
     Call get_word(record,word); cell(6)=word_2_real(word)

     Call get_line(safe,nconf,record); If (.not.safe) Go To 300
     Call get_word(record,word); cell(7)=word_2_real(word)
     Call get_word(record,word); cell(8)=word_2_real(word)
     Call get_word(record,word); cell(9)=word_2_real(word)

     Call invert(cell,rcell,det)

! Initialise domain localised atom counter (config_module)

     natms=0

! Initialise dispatched atom counter

     indatm=0

! Initialise total number of atoms counter

     nattot=0

     Do k=1,ntpmls
        Do l=1,nummols(k)
           Do m=1,numsit(k)

! Increase counters

              indatm=indatm+1
              nattot=nattot+1

! Initialise transmission arrays

              chbuf(indatm)=' '
              iwrk(indatm)=0

              axx(indatm)=0.0_wp
              ayy(indatm)=0.0_wp
              azz(indatm)=0.0_wp

! Read in transmission arrays

              If (idnode == 0 .and. safe) Then
                 record=' '; Read(Unit=nconf, Fmt='(a)', End=30) record
                 Call tabs_2_blanks(record) ; Call strip_blanks(record)
                 Call get_word(record,word) ; chbuf(indatm)=word(1:8)
                 If (l_ind) Then
                    Call get_word(record,word)
                    iwrk(indatm)=Nint(word_2_real(word,0.0_wp,l_str))
                    If (iwrk(indatm) /= 0) Then
                       iwrk(indatm)=Abs(iwrk(indatm))
                    Else
                       iwrk(indatm)=nattot
                    End If
                 Else
                    iwrk(indatm)=nattot
                 End If

                 If (levcfg /= 3) Then
                    Read(Unit=nconf, Fmt=*, End=30) axx(indatm),ayy(indatm),azz(indatm)

                    If (levcfg > 0) Then
                       Read(Unit=nconf, Fmt=*, End=30) bxx(indatm),byy(indatm),bzz(indatm)
                       If (levcfg > 1) Read(Unit=nconf, Fmt=*, End=30) cxx(indatm),cyy(indatm),czz(indatm)
                    End If
                 Else
                    Call get_word(record,word) ; axx(indatm)=word_2_real(word)
                    Call get_word(record,word) ; bxx(indatm)=word_2_real(word)
                    Call get_word(record,word) ; cxx(indatm)=word_2_real(word)
                 End If
                 Go To 40

30               Continue
                 safe=.false. ! catch error

40               Continue
              End If

! Circulate configuration data to all nodes when transmission arrays are filled up

              If (indatm == mxatms .or. nattot == megatm) Then

! Check if batch was read fine

                 If (mxnode > 1) Call gcheck(safe)
                 If (.not.safe) Go To 300 !Call error(25)

! Briadcaste all atoms

                 If (mxnode > 1) Then
                    Call MPI_BCAST(chbuf,indatm*8,MPI_CHARACTER,0,dlp_comm_world,ierr)
                    Call MPI_BCAST(iwrk,indatm,MPI_INTEGER,0,dlp_comm_world,ierr)

                    If (levcfg /= 3) Then
                       Call MPI_BCAST(axx,indatm,wp_mpi,0,dlp_comm_world,ierr)
                       Call MPI_BCAST(ayy,indatm,wp_mpi,0,dlp_comm_world,ierr)
                       Call MPI_BCAST(azz,indatm,wp_mpi,0,dlp_comm_world,ierr)
                       If (levcfg > 0) Then
                          Call MPI_BCAST(bxx,indatm,wp_mpi,0,dlp_comm_world,ierr)
                          Call MPI_BCAST(byy,indatm,wp_mpi,0,dlp_comm_world,ierr)
                          Call MPI_BCAST(bzz,indatm,wp_mpi,0,dlp_comm_world,ierr)
                          If (levcfg > 1) Then
                             Call MPI_BCAST(bxx,indatm,wp_mpi,0,dlp_comm_world,ierr)
                             Call MPI_BCAST(byy,indatm,wp_mpi,0,dlp_comm_world,ierr)
                             Call MPI_BCAST(bzz,indatm,wp_mpi,0,dlp_comm_world,ierr)
                          End If
                       End If
                    Else
                       Call MPI_BCAST(axx,indatm,wp_mpi,0,dlp_comm_world,ierr)
                       Call MPI_BCAST(ayy,indatm,wp_mpi,0,dlp_comm_world,ierr)
                       Call MPI_BCAST(azz,indatm,wp_mpi,0,dlp_comm_world,ierr)
                    End If
                 End If

! Assign atoms to correct domains (DD bound)

                 Do i=1,indatm
                    sxx=rcell(1)*axx(i)+rcell(4)*ayy(i)+rcell(7)*azz(i)
                    syy=rcell(2)*axx(i)+rcell(5)*ayy(i)+rcell(8)*azz(i)
                    szz=rcell(3)*axx(i)+rcell(6)*ayy(i)+rcell(9)*azz(i)

! sxx,syy,szz are in [-0.5,0.5) interval and values as 0.4(9) may pose a problem

                    sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
                    syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
                    szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

! fold back coordinates

                    axx(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
                    ayy(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
                    azz(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz

! assign domain coordinates (call for errors)

                    ipx=Int((sxx+0.5_wp)*nprx_r)
                    ipy=Int((syy+0.5_wp)*npry_r)
                    ipz=Int((szz+0.5_wp)*nprz_r)

                    idm=ipx+nprx*(ipy+npry*ipz)
                    If      (idm < 0 .or. idm > (mxnode-1)) Then
                       Call error(513)
                    Else If (idm == idnode)                 Then
                       natms=natms+1

                       If (natms < mxatms) Then
                          atmnam(natms)=chbuf(i)
                          ltg(natms)=iwrk(i)

                          If (levcfg /= 3) Then
                             xxx(natms)=axx(i)
                             yyy(natms)=ayy(i)
                             zzz(natms)=azz(i)
                             If (levcfg > 0) Then
                                vxx(natms)=bxx(i)
                                vyy(natms)=byy(i)
                                vzz(natms)=bzz(i)
                                If (levcfg > 1) Then
                                   fxx(natms)=cxx(i)
                                   fyy(natms)=cyy(i)
                                   fzz(natms)=czz(i)
                                End If
                             End If
                          Else
                             xxx(natms)=axx(i)
                             yyy(natms)=ayy(i)
                             zzz(natms)=azz(i)
                          End If
                       Else
                          safe=.false.
                       End If
                    End If
                 End Do

! Check if all is dispatched fine

                 If (mxnode > 1) Call gcheck(safe)
                 If (.not.safe) Call error(45)

! Nullify dispatch counter

                 indatm=0

              End If
           End Do
        End Do
     End Do

     Deallocate (chbuf,       Stat=fail(1))
     Deallocate (iwrk,        Stat=fail(2))
     Deallocate (axx,ayy,azz, Stat=fail(3))
     Deallocate (bxx,byy,bzz, Stat=fail(4))
     Deallocate (cxx,cyy,czz, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_history deallocation failure, node: ', idnode
        Call error(0)
     End If

! PROPER ASCII read

  Else If (io_read /= IO_READ_NETCDF) Then

     Call io_read_batch( fh, top_skip, 4, buffer, ierr )
     If (ierr < 0) Go To 300
     top_skip = top_skip + Int(4,MPI_OFFSET_KIND)

     record = ' '
     Do i = 1, Min( Size( buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = buffer( i, 1 )
     End Do
     Call get_word(record,word) ! timestep
     Call get_word(record,word) ; nstep = Nint(word_2_real(word))
     Call get_word(record,word) ; If (Nint(word_2_real(word)) /= megatm) Go To 300
     Call get_word(record,word) ; levcfg = Nint(word_2_real(word))
     Call get_word(record,word) ; imcon = Nint(word_2_real(word))

! image conditions not compliant with DD and link-cell

    If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

     Call get_word(record,word) ; tstep = word_2_real(word)
     Call get_word(record,word) ; time = word_2_real(word)

     If (idnode == 0) Write(nrite,"(/,1x,'HISTORY step',i10,' (',f10.3,' ps) is being read')") nstep,time

! read cell vectors

     record = ' '
     Do i = 1, Min( Size( buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = buffer( i, 2 )
     End Do
     Call get_word(record,word); cell(1)=word_2_real(word)
     Call get_word(record,word); cell(2)=word_2_real(word)
     Call get_word(record,word); cell(3)=word_2_real(word)

     record = ' '
     Do i = 1, Min( Size( buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = buffer( i, 3 )
     End Do
     Call get_word(record,word); cell(4)=word_2_real(word)
     Call get_word(record,word); cell(5)=word_2_real(word)
     Call get_word(record,word); cell(6)=word_2_real(word)

     record = ' '
     Do i = 1, Min( Size( buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = buffer( i, 4 )
     End Do
     Call get_word(record,word); cell(7)=word_2_real(word)
     Call get_word(record,word); cell(8)=word_2_real(word)
     Call get_word(record,word); cell(9)=word_2_real(word)

     Call read_config_parallel                  &
           (levcfg, dvar, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi)

     If (fast) Then
        If (levcfg /= 3) Then
           i=levcfg+2
        Else
           i=levcfg-2
        End If

        top_skip = top_skip + Int(i,MPI_OFFSET_KIND)*Int(megatm,MPI_OFFSET_KIND)
     Else
        top_skip = Int(0,MPI_OFFSET_KIND)
     End If

! Update current frame and exit gracefully

     frm1 = frm1+Int(1,li)
     If (frm1 == frm) Go To 200

  Else ! netCDF read

! Update current frame

     frm1 = frm1+Int(1,li)
     i = Int( frm1 )

     Call io_nc_get_var( 'time'           , fh,   time, i, 1 )
     Call io_nc_get_var( 'datalevel'      , fh, levcfg, i, 1 )
     Call io_nc_get_var( 'imageconvention', fh,  imcon, i, 1 )

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

     Call io_nc_get_var( 'timestep'       , fh,  tstep, i, 1 )
     Call io_nc_get_var( 'step'           , fh,  nstep, i, 1 )

     If (idnode == 0) Write(nrite,"(/,1x,'HISTORY step',i10,' (',f10.3,' ps) is being read')") nstep,time

! Note that in netCDF the frames are not long integers - Int( frm1 )

     Call io_nc_get_var( 'cell'           , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
     cell = Reshape( cell_vecs, (/ Size( cell ) /) )

     Call read_config_parallel                  &
           (levcfg, dvar, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, Int( i, Kind( top_skip ) ), xhi, yhi, zhi)

     If (frm1 == frm) Go To 200

  End If

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  nlast=natms

! Does the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?
! Check number of atoms in system (CONFIG = FIELD)

  totatm=natms
  If (mxnode > 1) Call gsum(totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (config_module)

  Do i=1,natms
     lsi(i)=i
     lsa(i)=ltg(i)
  End Do
  Call shellsort2(natms,lsi,lsa)

  exout = 0 ! more to read indicator

  Return

! Normal exit from HISTORY file read

200 Continue

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  nlast=natms

! Does the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?
! Check number of atoms in system (CONFIG = FIELD)

  totatm=natms
  If (mxnode > 1) Call gsum(totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (config_module)

  Do i=1,natms
     lsi(i)=i
     lsa(i)=ltg(i)
  End Do
  Call shellsort2(natms,lsi,lsa)

  exout = 1 ! It's an indicator of the end of reading.

  If (idnode == 0) Write(nrite,"(1x,a)") 'HISTORY end of file reached'

  If (io_read == IO_READ_MASTER) Then
     If (imcon == 0) Close(Unit=nconf)
     Deallocate (chbuf,       Stat=fail(1))
     Deallocate (iwrk,        Stat=fail(2))
     Deallocate (axx,ayy,azz, Stat=fail(3))
     Deallocate (bxx,byy,bzz, Stat=fail(4))
     Deallocate (cxx,cyy,czz, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_history deallocation failure, node: ', idnode
        Call error(0)
     End If
  Else
     Call io_close( fh )
     Call io_finalize
  End If

  Return

! Abnormal exit from HISTORY file read

300 Continue

  If (idnode == 0) Write(nrite,"(/,1x,a)") 'HISTORY data mishmash detected'
  exout = -1 ! It's an indicator of the end of reading.
  If (io_read == IO_READ_MASTER) Then
     If (imcon == 0) Close(Unit=nconf)
     Deallocate (chbuf,       Stat=fail(1))
     Deallocate (iwrk,        Stat=fail(2))
     Deallocate (axx,ayy,azz, Stat=fail(3))
     Deallocate (bxx,byy,bzz, Stat=fail(4))
     Deallocate (cxx,cyy,czz, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_history deallocation failure, node: ', idnode
        Call error(0)
     End If
  Else
     Call io_close( fh )
     Call io_finalize
  End If

  Return

! HISTORY not found

400 Continue

  Call error(585)

End Subroutine read_history
