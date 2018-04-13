Module trajectory
  Use kinds,         Only : wp, li
  Use comms,         Only : comms_type,Traject_tag,gsync,wp_mpi,gbcast,gcheck, &
                            gsum,gsend,grecv,offset_kind,mode_wronly, &
                            mode_rdonly, comm_self
  Use domains,       Only : nprx,npry,nprz,nprx_r,npry_r,nprz_r
  Use site
  Use setup
  Use parse,         Only : tabs_2_blanks, get_line, get_word, &
                            strip_blanks, word_2_real
  Use configuration, Only : cfgname,imcon,cell,natms, &
                            ltg,atmnam,chge,weight,   &
                            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,&
                            lsa,lsi,ltg,nlast
  Use statistics,    Only : rsd
  Use io,            Only : io_set_parameters,             &
                            io_get_parameters,             &
                            io_init, io_nc_create,         &
                            io_open, io_write_record,      &
                            io_write_batch,io_read_batch,  &
                            io_nc_get_dim,                 &
                            io_nc_put_var,io_nc_get_var,   &
                            io_write_sorted_file,          &
                            io_close, io_finalize,         &
                            io_nc_get_real_precision,      &
                            io_nc_get_file_real_precision, &
                            IO_HISTORY,IO_HISTORD,         &
                            IO_BASE_COMM_NOT_SET,          &
                            IO_ALLOCATION_ERROR,           &
                            IO_UNKNOWN_WRITE_OPTION,       &
                            IO_UNKNOWN_WRITE_LEVEL,        &
                            IO_WRITE_UNSORTED_MPIIO,       &
                            IO_WRITE_UNSORTED_DIRECT,      &
                            IO_WRITE_UNSORTED_MASTER,      &
                            IO_WRITE_SORTED_MPIIO,         &
                            IO_WRITE_SORTED_DIRECT,        &
                            IO_WRITE_SORTED_NETCDF,        &
                            IO_WRITE_SORTED_MASTER,        &
                            IO_READ_MPIIO,         &
                            IO_READ_DIRECT,        &
                            IO_READ_NETCDF,        &
                            IO_READ_MASTER
  Use numerics,        Only : dcell, invert, shellsort2
  Use configuration,   Only : read_config_parallel
  Use errors_warnings, Only : error
  Implicit None

  Private
  Public :: read_history
  Public :: trajectory_write
Contains


Subroutine read_history(l_str,fname,megatm,levcfg,dvar,nstep,tstep,time,exout,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading the trajectory data file
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Character( Len = * ), Intent( In    ) :: fname
  Logical,              Intent( In    ) :: l_str
  Integer,              Intent( In    ) :: megatm

  Integer,              Intent( InOut ) :: levcfg,nstep
  Real( Kind = wp ),    Intent( In    ) :: dvar
  Real( Kind = wp ),    Intent( InOut ) :: tstep,time
  Integer,              Intent(   Out ) :: exout
  Type( comms_type),    Intent( InOut ) :: comm

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

! Some parameters and variables needed by io interfaces

  Integer,                           Save :: recsz ! record size
  Integer,                           Save :: fh, io_read
  Integer( Kind = offset_kind ), Save :: top_skip
  Character( Len = 1),  Allocatable, Save :: buffer(:,:)

  Character( Len = 8 ), Allocatable :: chbuf(:)
  Integer,              Allocatable :: iwrk(:)
  Real( Kind = wp ),    Allocatable :: axx(:),ayy(:),azz(:), &
                                       bxx(:),byy(:),bzz(:), &
                                       cxx(:),cyy(:),czz(:)
  Integer  :: ierr
  Character ( Len = 256 )  :: message

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
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
        If (.not.lexist) Go To 400

! Open HISTORY

        If (comm%idnode == 0) Open(Unit=nconf, File=fname)

! read the HISTORY file header

        Call get_line(safe,nconf,record,comm); If (.not.safe) Go To 300

        Call get_line(safe,nconf,record,comm); If (.not.safe) Go To 300
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
           top_skip = Int(2,offset_kind)

           fail(1) = 0
           Allocate (buffer(1:recsz,1:4), Stat=fail(1))
           If (fail(1) > 0) Then
              Write(message,'(a)') 'read_history allocation failure 1'
              Call error(0,message)
           End If

           If (io_read == IO_READ_MPIIO) Then
              Close(Unit=nconf)

              Call io_set_parameters( user_comm = comm%comm )
              Call io_init( recsz )
              Call io_open( io_read, comm%comm, fname, mode_rdonly, fh )
           End If
        End If

     Else ! netCDF read

! Does HISTORY exist

        lexist=.true.
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
        If (.not.lexist) Go To 400

! fast and rec are irrelevant for netCDF (initialised at declaration)

        fast = .true.
        rec  = Int(0,li)

        Call io_set_parameters( user_comm = comm%comm )
        Call io_open( io_read, comm%comm, fname, mode_rdonly, fh )
        Call io_nc_get_dim( 'frame', fh, i )

        If (i > 0) Then
           frm=Int(i,li)
        Else
           Go To 300
        End If

     End If
  Else
     If (io_read == IO_READ_MPIIO .or. io_read == IO_READ_NETCDF) &
        Call io_set_parameters( user_comm = comm%comm )
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
        Write(message,'(a)') 'read_history allocation failure'
        Call error(0,message)
     End If

! read timestep and time

     Call get_line(safe,nconf,record,comm); If (.not.safe) Go To 200

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

     If (comm%idnode == 0) Write(nrite,"(/,1x,'HISTORY step',i10,' (',f10.3,' ps) is being read')") nstep,time

! read cell vectors

     Call get_line(safe,nconf,record,comm); If (.not.safe) Go To 300
     Call get_word(record,word); cell(1)=word_2_real(word)
     Call get_word(record,word); cell(2)=word_2_real(word)
     Call get_word(record,word); cell(3)=word_2_real(word)

     Call get_line(safe,nconf,record,comm); If (.not.safe) Go To 300
     Call get_word(record,word); cell(4)=word_2_real(word)
     Call get_word(record,word); cell(5)=word_2_real(word)
     Call get_word(record,word); cell(6)=word_2_real(word)

     Call get_line(safe,nconf,record,comm); If (.not.safe) Go To 300
     Call get_word(record,word); cell(7)=word_2_real(word)
     Call get_word(record,word); cell(8)=word_2_real(word)
     Call get_word(record,word); cell(9)=word_2_real(word)

     Call invert(cell,rcell,det)

! Initialise domain localised atom counter (configuration)

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

              If (comm%idnode == 0 .and. safe) Then
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

                 Call gcheck(comm,safe)
                 If (.not.safe) Go To 300 !Call error(25)

! Briadcaste all atoms

                    Call gbcast(comm,chbuf,0)
                    Call gbcast(comm,iwrk,0)

                    If (levcfg /= 3) Then

                       Call gbcast(comm,axx,0)
                       Call gbcast(comm,ayy,0)
                       Call gbcast(comm,azz,0)
                       If (levcfg > 0) Then
                         Call gbcast(comm,bxx,0)
                         Call gbcast(comm,byy,0)
                         Call gbcast(comm,bzz,0)
                          If (levcfg > 1) Then
                             Call gbcast(comm,cxx,0)
                             Call gbcast(comm,cyy,0)
                             Call gbcast(comm,czz,0)
                          End If
                       End If
                    Else
                       Call gbcast(comm,axx,0)
                       Call gbcast(comm,ayy,0)
                       Call gbcast(comm,azz,0)
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

! fold back coordinatesc

                    axx(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
                    ayy(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
                    azz(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz

! assign domain coordinates (call for errors)

                    ipx=Int((sxx+0.5_wp)*nprx_r)
                    ipy=Int((syy+0.5_wp)*npry_r)
                    ipz=Int((szz+0.5_wp)*nprz_r)

                    idm=ipx+nprx*(ipy+npry*ipz)
                    If      (idm < 0 .or. idm > (comm%mxnode-1)) Then
                       Call error(513)
                    Else If (idm == comm%idnode)                 Then
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

                 Call gcheck(comm,safe)
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
        Write(message,'(a)') 'read_history deallocation failure'
        Call error(0,message)
     End If

! PROPER ASCII read

  Else If (io_read /= IO_READ_NETCDF) Then

     Call io_read_batch( fh, top_skip, 4, buffer, ierr )
     If (ierr < 0) Go To 300
     top_skip = top_skip + Int(4,offset_kind)

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

     If (comm%idnode == 0) Write(nrite,"(/,1x,'HISTORY step',i10,' (',f10.3,' ps) is being read')") nstep,time

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
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi,comm)

     If (fast) Then
        If (levcfg /= 3) Then
           i=levcfg+2
        Else
           i=levcfg-2
        End If

        top_skip = top_skip + Int(i,offset_kind)*Int(megatm,offset_kind)
     Else
        top_skip = Int(0,offset_kind)
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

     If (comm%idnode == 0) Write(nrite,"(/,1x,'HISTORY step',i10,' (',f10.3,' ps) is being read')") nstep,time

! Note that in netCDF the frames are not long integers - Int( frm1 )

     Call io_nc_get_var( 'cell'           , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
     cell = Reshape( cell_vecs, (/ Size( cell ) /) )

     Call read_config_parallel                  &
           (levcfg, dvar, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, Int( i, Kind( top_skip ) ), xhi, yhi, zhi,comm)

     If (frm1 == frm) Go To 200

  End If

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  nlast=natms

! Does the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?
! Check number of atoms in system (CONFIG = FIELD)

  totatm=natms
  Call gsum(comm,totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (configuration)

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
  Call gsum(comm,totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (configuration)

  Do i=1,natms
     lsi(i)=i
     lsa(i)=ltg(i)
  End Do
  Call shellsort2(natms,lsi,lsa)

  exout = 1 ! It's an indicator of the end of reading.

  If (comm%idnode == 0) Write(nrite,"(1x,a)") 'HISTORY end of file reached'

  If (io_read == IO_READ_MASTER) Then
     If (imcon == 0) Close(Unit=nconf)
     Deallocate (chbuf,       Stat=fail(1))
     Deallocate (iwrk,        Stat=fail(2))
     Deallocate (axx,ayy,azz, Stat=fail(3))
     Deallocate (bxx,byy,bzz, Stat=fail(4))
     Deallocate (cxx,cyy,czz, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'read_history deallocation failure'
        Call error(0,message)
     End If
  Else
     Call io_close( fh )
     Call io_finalize
  End If

  Return

! Abnormal exit from HISTORY file read

300 Continue

  If (comm%idnode == 0) Write(nrite,"(/,1x,a)") 'HISTORY data mishmash detected'
  exout = -1 ! It's an indicator of the end of reading.
  If (io_read == IO_READ_MASTER) Then
     If (imcon == 0) Close(Unit=nconf)
     Deallocate (chbuf,       Stat=fail(1))
     Deallocate (iwrk,        Stat=fail(2))
     Deallocate (axx,ayy,azz, Stat=fail(3))
     Deallocate (bxx,byy,bzz, Stat=fail(4))
     Deallocate (cxx,cyy,czz, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'read_history deallocation failure'
        Call error(0,message)
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

Subroutine trajectory_write(keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing history file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
! contrib   - w.smith, i.j.bush
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: keyres,               &
                                        nstraj,istraj,keytrj, &
                                        megatm,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep,time
  Type( comms_type ), Intent( InOut ) :: comm


  Logical,               Save :: newjob = .true. , &
                                 fast   = .true.
  Character( Len = 40 ), Save :: fname
  Integer,               Save :: recsz  = 73 ! default record size
  Integer(Kind=li),      Save :: rec    = 0_li , &
                                 frm    = 0_li

  Logical                :: lexist,safe,ready
  Character( Len = 40 )  :: word
  Integer                :: fail(1:5),i,jj,k,jdnode,jatms
  Integer(Kind=li)       :: rec1
  Real( Kind = wp )      :: buffer(1:2)
  Real( Kind = wp )      :: celprp(1:10),cell_vecs(1:3,1:3)
  Real( Kind = wp )      :: lengths(1:3), angles(1:3)

! Some parameters and variables needed by io interfaces

  Integer                           :: fh, io_write, batsz
  Integer( Kind = offset_kind ) :: rec_mpi_io,jj_io
  Character                         :: lf
  Character( Len = 200 )            :: record

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: cxx,cyy,czz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: ddd,eee,fff

! netCDF check

  Integer :: file_p, file_r
  Integer :: io_p, io_r,ierr
 
  Character ( Len = 256 )  ::  message

  If (.not.(nstep >= nstraj .and. Mod(nstep-nstraj,istraj) == 0)) Return

! Get write method, buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

  If (keytrj == 3) Then
     recsz = 35 ! record size for compressed HISTORY file
     Go To 100
  End If

  If (newjob) Then
     newjob = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = Trim(history)
     Else
        fname = Trim(history) // '.nc'
     End If

! If keyres=1, is HISTORY old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (comm%idnode == 0) Then
              Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Status='replace', Recl=recsz)
              Write(Unit=nhist, Fmt='(a72,a1)',       Rec=Int(1,li)) cfgname(1:72),lf
              Write(Unit=nhist, Fmt='(3i10,2i21,a1)', Rec=Int(2,li)) keytrj,imcon,megatm,frm,rec,lf
              Close(Unit=nhist)
           End If
           rec=Int(2,li)
           frm=Int(0,li)
        Else
           If (comm%idnode == 0) Then
              Call io_set_parameters( user_comm = comm_self )
              Call io_nc_create( comm_self, fname, cfgname, megatm )
           End If
        End If

! Get some sense out of it

     Else

        safe=.true.
        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           If (comm%idnode == 0) Then

              Open(Unit=nhist, File=fname, Form='formatted')

              Do

                 record(1:recsz)=' '

! Assume new style of HISTORY with bookkeeping.

                 If (fast) Then

                    Read(Unit=nhist, Fmt=*, End=20)                     ! title record
                    rec=rec+Int(1,li)
                    Read(Unit=nhist, Fmt='(a)', End=20) record(1:recsz) ! bookkeeping record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    rec=rec+Int(1,li)

                    Call get_word(record(1:recsz),word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record(1:recsz),word) ; Call get_word(record(1:recsz),word)
                       Call get_word(record(1:recsz),word) ; frm=Nint(word_2_real(word,0.0_wp),li)
                       Call get_word(record(1:recsz),word) ; rec=Nint(word_2_real(word,0.0_wp),li)
                       If (frm /= Int(0,li) .and. rec > Int(2,li)) Then
                          Go To 20 ! New style
                       Else
                          fast=.false. ! TOUGH, old style
                          rec=Int(2,li)
                          frm=Int(0,li)
                       End If
                    Else
                       safe=.false. ! Overwrite the file, it's junk to me
                       Go To 20
                    End If

! TOUGH, it needs scanning through

                 Else

                    Read(Unit=nhist, Fmt='(a)', End=20) record(1:recsz) ! timestep record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    rec=rec+Int(1,li)

                    Call get_word(record(1:recsz),word) ; Call get_word(record(1:recsz),word)
                    Call get_word(record(1:recsz),word) ; jj=Nint(word_2_real(word))
                    Call get_word(record(1:recsz),word) ; k=Nint(word_2_real(word))

                    word=' '
                    i = 3 + (2+k)*jj ! total number of lines to read
                    Write(word,'( "(", i0, "( / ) )" )') i-1
                    Read(Unit=nhist, Fmt=word, End=20)
                    rec=rec+Int(i,li)
                    frm=frm+Int(1,li)

                 End If

              End Do

20            Continue
              Close(Unit=nhist)

           End If

           Call gcheck(comm,safe,"enforce")
           If (.not.safe) Then
              lexist=.false.

              rec=Int(0,li)
              frm=Int(0,li)

              Go To 10
           Else If (comm%mxnode > 1) Then
              buffer(1)=Real(frm,wp)
              buffer(2)=Real(rec,wp)

              Call gbcast(comm,buffer,0)

              frm=Nint(buffer(1),li)
              rec=Nint(buffer(2),li)
           End If

        Else ! netCDF read

           If (comm%idnode == 0) Then
              Call io_set_parameters( user_comm = comm_self )
              Call io_open( io_write, comm_self, fname, mode_rdonly, fh )

! Get the precision that the history file was written in
! and check it matches the requested precision

              Call io_nc_get_file_real_precision( fh, file_p, file_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(comm,safe)
           If (.not.safe) Then
              If (comm%idnode == 0) Write(message,'(a)') &
  "Can not determine precision in an exisiting HISTORY.nc file in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync(comm)
              Call error(0,message)
           End If

           If (comm%idnode == 0) Then
              Call io_nc_get_real_precision( io_p, io_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(comm,safe)
           If (.not.safe) Then
              If (comm%idnode == 0) Write(message,'(a)') &
  "Can not determine the desired writing precision in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync(comm)
              Call error(0,message)
           End If

           If (comm%idnode == 0) safe = (io_p == file_p .and. io_r == file_r)
           Call gcheck(comm,safe)
           If (.not.safe) Then
              If (comm%idnode == 0) Then
                 Write(message,'(a)') &
  "Requested writing precision inconsistent with that in an existing HISTORY.nc"
                 Write(nrite, Fmt='(1x,a)', Advance='No') "Precision requested:"
                 Select Case( Selected_real_kind( io_p, io_r ) )
                 Case( Kind( 1.0 ) )
                    Write(message,'(1x,a)') "Single"
                 Case( Kind( 1.0d0 ) )
                    Write(message,'(1x,a)') "Double"
                 End Select
                 Write(nrite, Fmt='(1x,a)', Advance='No') "Precision in file  :"
                 Select Case( Selected_real_kind( file_p, file_r ) )
                 Case( Kind( 1.0 ) )
                    Write(message,'(1x,a)') "Single"
                 Case( Kind( 1.0d0 ) )
                    Write(message,'(1x,a)') "Double"
                 End Select
              End If

              Call gsync(comm)
              Call error(0,message)
           End If

! Get the frame number to check
! For netCDF this is the "frame number" which is not a long integer!

           If (comm%idnode == 0) Call io_nc_get_dim( 'frame', fh, jj )
           Call gbcast(comm,jj,0)

           If (jj > 0) Then
              frm=Int(jj,li)

              If (comm%idnode == 0) Call io_close( fh )
           Else ! Overwrite the file, it's junk to me
              lexist=.false.

              rec=Int(0,li)
              frm=Int(0,li)

              Go To 10
           End If
        End If

     End If
  End If

! Get offsets and define batch

  fail=0
  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Allocate (n_atm(0:comm%mxnode),        Stat=fail(1))
     Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure 0'
        Call error(0,message)
     End If

     chbat=' '
     n_atm=0 ; n_atm(comm%idnode+1)=natms
     Call gsum(comm,n_atm)
     n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = offset_kind)

! Update frame

  frm=frm+Int(1,li)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,offset_kind)
     jj=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_open( io_write, comm_self, fname, mode_wronly, fh )

        Write(record(1:recsz), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record(1:recsz), Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and cell information

        Call io_write_batch( fh, rec_mpi_io, jj, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        jj=jj+4

     End If
     Call gsync(comm)

! Start of file

     rec_mpi_io=Int(rec,offset_kind)+Int(jj,offset_kind)+Int(n_atm(0),offset_kind)*Int(keytrj+2,offset_kind)
     jj=0

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, fname, mode_wronly, fh )

     Do i=1,natms
        Write(record(1:recsz), Fmt='(a8,i10,3f12.6,a18,a1)') atmnam(i),ltg(i),weight(i),chge(i),rsd(i),Repeat(' ',18),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') xxx(i),yyy(i),zzz(i),Repeat(' ',12),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        If (keytrj > 0) Then
           Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') vxx(i),vyy(i),vzz(i),Repeat(' ',12),lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End If

        If (keytrj > 1) Then
           Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') fxx(i),fyy(i),fzz(i),Repeat(' ',12),lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End If

! Dump batch and update start of file

        If (jj + keytrj + 2 >= batsz .or. i == natms) Then
           Call io_write_batch( fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,offset_kind)
           jj=0
        End If
     End Do

! Update and save offset pointer

     rec=rec+Int(4,li)+Int(megatm,li)*Int(keytrj+2,li)
     If (comm%idnode == 0) Then
        Write(record(1:recsz), Fmt='(3i10,2i21,a1)') keytrj,imcon,megatm,frm,rec,lf
        Call io_write_record( fh, Int(1,offset_kind), record(1:recsz) )
     End If

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     Allocate (ddd(1:mxatms),eee(1:mxatms),fff(1:mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     jj=0
     If (comm%idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

! Accumulate header

        Write(record(1:recsz), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record(1:recsz), Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=nhist, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
        rec=rec+Int(jj,li)
        jj=0

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           chbuf(i)=atmnam(i)
           ddd(i)=chge(i)
           eee(i)=weight(i)
           fff(i)=rsd(i)
        End Do

        If (keytrj >= 1) Then
           Do i=1,natms
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)
           End Do
        End If

        If (keytrj >= 2) Then
           Do i=1,natms
              cxx(i)=fxx(i)
              cyy(i)=fyy(i)
              czz(i)=fzz(i)
           End Do
        End If

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,Traject_tag)

              Call grecv(comm,jatms,jdnode,Traject_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,ddd(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,eee(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,fff(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,axx(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,ayy(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,azz(1:jatms),jdnode,Traject_tag)

                 If (keytrj > 0) Then
                    Call grecv(comm,bxx(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,byy(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,bzz(1:jatms),jdnode,Traject_tag)
                 End If

                 If (keytrj > 1) Then
                    Call grecv(comm,cxx(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,cyy(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,czz(1:jatms),jdnode,Traject_tag)
                 End If
              End If
           End If

           Do i=1,jatms
              Write(record(1:recsz), Fmt='(a8,i10,3f12.6,a18,a1)') chbuf(i),iwrk(i),eee(i),ddd(i),fff(i),Repeat(' ',18),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              If (keytrj >= 1) Then
                 Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,recsz
                    chbat(k,jj) = record(k:k)
                 End Do
              End If

              If (keytrj >= 2) Then
                 Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,recsz
                    chbat(k,jj) = record(k:k)
                 End Do
              End If

! Dump batch and update start of file

              If (jj + keytrj + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nhist, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
                 rec=rec+Int(jj,li)
                 jj=0
              End If
           End Do
        End Do

! Update main header

        Write(Unit=nhist, Fmt='(3i10,2i21,a1)', Rec=Int(2,li)) keytrj,imcon,megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Call grecv(comm,ready,0,Traject_tag)

        Call gsend(comm,natms,0,Traject_tag)
        If (natms > 0) Then
           Call gsend(comm,atmnam(:),0,Traject_tag)
           Call gsend(comm,ltg(:),0,Traject_tag)

           Call gsend(comm,chge(:),0,Traject_tag)
           Call gsend(comm,weight(:),0,Traject_tag)
           Call gsend(comm,rsd(:),0,Traject_tag)

           Call gsend(comm,xxx(:),0,Traject_tag)
           Call gsend(comm,yyy(:),0,Traject_tag)
           Call gsend(comm,zzz(:),0,Traject_tag)

           If (keytrj > 0) Then
              Call gsend(comm,vxx(:),0,Traject_tag)
              Call gsend(comm,vyy(:),0,Traject_tag)
              Call gsend(comm,vzz(:),0,Traject_tag)
           End If

           If (keytrj > 1) Then
              Call gsend(comm,fxx(:),0,Traject_tag)
              Call gsend(comm,fyy(:),0,Traject_tag)
              Call gsend(comm,fzz(:),0,Traject_tag)
           End If
        End If

! Save offset pointer

        rec=rec+Int(4,li)+Int(megatm,li)*Int(keytrj+2,li)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     Deallocate (ddd,eee,fff, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write deallocation failure'
        Call error(0,message)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,offset_kind)
     jj_io=rec_mpi_io
     jj=0 ! netCDF current frame
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_open( io_write, comm_self, fname, mode_wronly, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and cell information

           Write(record(1:recsz), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf
           Call io_write_record( fh, jj_io, record(1:recsz) )
           jj_io=jj_io + Int(1,offset_kind)

           Do i = 0, 2
              Write(record(1:recsz), Fmt='(3f20.10,a12,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              Call io_write_record( fh, jj_io, record(1:recsz) )
              jj_io=jj_io+Int(1,offset_kind)
           End Do

        Else ! netCDF write

! Get the current and new frame numbers

           Call io_nc_get_dim( 'frame', fh, jj )
           jj = jj + 1

           Call io_nc_put_var( 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var( 'datalevel'      , fh, keytrj, jj, 1 )
           Call io_nc_put_var( 'imageconvention', fh,  imcon, jj, 1 )
           Call io_nc_put_var( 'timestep '      , fh,  tstep, jj, 1 )

           Call dcell(cell,celprp) ! get cell properties

           cell_vecs = Reshape( cell, (/ 3, 3 /) )

           lengths( 1 ) = celprp( 1 )
           lengths( 2 ) = celprp( 2 )
           lengths( 3 ) = celprp( 3 )

           angles ( 1 ) = Acos( celprp( 5 ) )
           angles ( 2 ) = Acos( celprp( 6 ) )
           angles ( 3 ) = Acos( celprp( 4 ) )
           angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

! Print

           Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, jj /), (/ 3, 3, 1 /) )
           Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, jj /), (/    3, 1 /) )
           Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, jj /), (/    3, 1 /) )

        End If

        Call io_close( fh )
        Call io_finalize

     End If
     Call gsync(comm)

! Start of file

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec_mpi_io=Int(rec,offset_kind)+Int(4,offset_kind)
     Else ! netCDF write
        Call gbcast(comm,jj,0)
        rec_mpi_io = Int(jj,offset_kind)
     End If

! Write the rest

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, fname, mode_wronly, fh )

     Call io_write_sorted_file( fh, keytrj, IO_HISTORY, rec_mpi_io, natms, &
          ltg, atmnam, weight, chge, rsd, xxx, yyy, zzz,                   &
          vxx, vyy, vzz, fxx, fyy, fzz, ierr )

     If ( ierr /= 0 ) Then
        Select Case( ierr )
        Case( IO_BASE_COMM_NOT_SET )
           Call error( 1050 )
        Case( IO_ALLOCATION_ERROR )
           Call error( 1053 )
        Case( IO_UNKNOWN_WRITE_OPTION )
           Call error( 1056 )
        Case( IO_UNKNOWN_WRITE_LEVEL )
           Call error( 1059 )
        End Select
     End If

! Update and save offset pointer

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec=rec+Int(4,li)+Int(megatm,li)*Int(keytrj+2,li)
        If (comm%idnode == 0) Then
           Write(record(1:recsz), Fmt='(3i10,2i21,a1)') keytrj,imcon,megatm,frm,rec,lf
           Call io_write_record( fh, Int(1,offset_kind), record(1:recsz) )
        End If
     End If

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     Allocate (ddd(1:mxatms),eee(1:mxatms),fff(1:mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O

     If (comm%idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,li)
        Write(Unit=nhist, Fmt='(a8,2i10,2i2,2f20.6,a1)', Rec=rec) 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf

        Do i = 0, 2
           rec=rec+Int(1,li)
           Write(Unit=nhist, Fmt='(3f20.10,a12,a1)', Rec=rec) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
        End Do

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           chbuf(i)=atmnam(i)
           ddd(i)=chge(i)
           eee(i)=weight(i)
           fff(i)=rsd(i)
        End Do

        If (keytrj >= 1) Then
           Do i=1,natms
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)
           End Do
        End If

        If (keytrj >= 2) Then
           Do i=1,natms
              cxx(i)=fxx(i)
              cyy(i)=fyy(i)
              czz(i)=fzz(i)
           End Do
        End If

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,Traject_tag)

              Call grecv(comm,jatms,jdnode,Traject_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,ddd(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,eee(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,fff(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,axx(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,ayy(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,azz(1:jatms),jdnode,Traject_tag)

                 If (keytrj > 0) Then
                    Call grecv(comm,bxx(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,byy(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,bzz(1:jatms),jdnode,Traject_tag)
                 End If

                 If (keytrj > 1) Then
                    Call grecv(comm,cxx(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,cyy(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,czz(1:jatms),jdnode,Traject_tag)
                 End If
              End If
           End If

           Do i=1,jatms
              rec1=rec+Int(iwrk(i)-1,li)*Int(keytrj+2,li)+Int(1,li)
              Write(Unit=nhist, Fmt='(a8,i10,3f12.6,a18,a1)', Rec=rec1) chbuf(i),iwrk(i),eee(i),ddd(i),fff(i),Repeat(' ',18),lf

              rec1=rec1+Int(1,li)
              Write(Unit=nhist, Fmt='(3g20.10,a12,a1)', Rec=rec1) axx(i),ayy(i),azz(i),Repeat(' ',12),lf

              If (keytrj >= 1) Then
                 rec1=rec1+Int(1,li)
                 Write(Unit=nhist, Fmt='(3g20.10,a12,a1)', Rec=rec1) bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
              End If

              If (keytrj >= 2) Then
                 rec1=rec1+Int(1,li)
                 Write(Unit=nhist, Fmt='(3g20.10,a12,a1)', Rec=rec1) cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
              End If
           End Do
        End Do

! Update main header

        rec=rec+Int(megatm,li)*Int(keytrj+2,li)
        Write(Unit=nhist, Fmt='(3i10,2i21,a1)', Rec=Int(2,li)) keytrj,imcon,megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Call grecv(comm,ready,0,Traject_tag)

        Call gsend(comm,natms,0,Traject_tag)
        If (natms > 0) Then
           Call gsend(comm,atmnam(:),0,Traject_tag)
           Call gsend(comm,ltg(:),0,Traject_tag)

           Call gsend(comm,chge(:),0,Traject_tag)
           Call gsend(comm,weight(:),0,Traject_tag)
           Call gsend(comm,rsd(:),0,Traject_tag)

           Call gsend(comm,xxx(:),0,Traject_tag)
           Call gsend(comm,yyy(:),0,Traject_tag)
           Call gsend(comm,zzz(:),0,Traject_tag)

           If (keytrj > 0) Then
              Call gsend(comm,vxx(:),0,Traject_tag)
              Call gsend(comm,vyy(:),0,Traject_tag)
              Call gsend(comm,vzz(:),0,Traject_tag)
           End If

           If (keytrj > 1) Then
              Call gsend(comm,fxx(:),0,Traject_tag)
              Call gsend(comm,fyy(:),0,Traject_tag)
              Call gsend(comm,fzz(:),0,Traject_tag)
           End If
        End If

! Save offset pointer

        rec=rec+Int(4,li)+Int(megatm,li)*Int(keytrj+2,li)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     Deallocate (ddd,eee,fff, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write deallocation failure'
        Call error(0,message)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write deallocation failure 0'
        Call error(0,message)
     End If
  End If

  Call gsync(comm)

  Return

100 Continue

  If (newjob) Then
     newjob = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = Trim(history)
     Else
        fname = Trim(history) // '.nc'
     End If

! If keyres=1, is HISTORY old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

110  Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (comm%idnode == 0) Then
              Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Status='replace', Recl=recsz)
              Write(Unit=nhist, Fmt='(a34,a1)',      Rec=Int(1,li)) cfgname(1:34),lf
              Write(Unit=nhist, Fmt='(2i2,3i10,a1)', Rec=Int(2,li)) keytrj,imcon,megatm,frm,rec,lf
              Close(Unit=nhist)
           End If
           rec=Int(2,li)
           frm=Int(0,li)
        Else
           If (comm%idnode == 0) Then
              Call io_set_parameters( user_comm = comm_self )
              Call io_nc_create( comm_self, fname, cfgname, megatm )
           End If
        End If

! Get some sense of it

     Else

        safe=.true.
        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           If (comm%idnode == 0) Then

              Open(Unit=nhist, File=fname, Form='formatted')

              Do

                 record(1:recsz)=' '

! Assume new style of HISTORY with bookkeeping.

                 If (fast) Then

                    Read(Unit=nhist, Fmt=*, End=120)                     ! title record
                    rec=rec+Int(1,li)
                    Read(Unit=nhist, Fmt='(a)', End=120) record(1:recsz) ! bookkeeping record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    rec=rec+Int(1,li)

                    Call get_word(record(1:recsz),word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record(1:recsz),word) ; Call get_word(record(1:recsz),word)
                       Call get_word(record(1:recsz),word) ; frm=Nint(word_2_real(word,0.0_wp),li)
                       Call get_word(record(1:recsz),word) ; rec=Nint(word_2_real(word,0.0_wp),li)
                       If (frm /= Int(0,li) .and. rec > Int(2,li)) Then
                          Go To 120 ! New style
                       Else
                          fast=.false. ! TOUGH, old style
                          rec=Int(2,li)
                          frm=Int(0,li)
                       End If
                    Else
                       safe=.false. ! Overwrite the file, it's junk to me
                       Go To 120
                    End If

! TOUGH, it needs scanning through

                 Else

                    Read(Unit=nhist, Fmt=*, End=120)                 ! timestep record
                    rec=rec+Int(1,li)

                    word=' '
                    i = 3 + megatm ! total number of lines to read
                    Write(word,'( "(", i0, "( / ) )" )') i-1
                    Read(Unit=nhist, Fmt=word, End=120)
                    rec=rec+Int(i,li)
                    frm=frm+Int(1,li)

                 End If

              End Do

120         Continue
            Close(Unit=nhist)

           End If

           Call gcheck(comm,safe,"enforce")
           If (.not.safe) Then
              lexist=.false.

              rec=Int(0,li)
              frm=Int(0,li)

              Go To 110
           Else If (comm%mxnode > 1) Then
              buffer(1)=Real(frm,wp)
              buffer(2)=Real(rec,wp)

              Call gbcast(comm,buffer,0)

              frm=Nint(buffer(1),li)
              rec=Nint(buffer(2),li)
           End If

        Else ! netCDF read

           If (comm%idnode == 0) Then
              Call io_set_parameters( user_comm = comm_self )
              Call io_open( io_write, comm_self, fname, mode_rdonly, fh )

! Get the precision that the history file was written in
! and check it matches the requested precision

              Call io_nc_get_file_real_precision( fh, file_p, file_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(comm,safe)
           If (.not.safe) Then
              If (comm%idnode == 0) Write(message,'(a)') &
  "Can not determine precision in an exisiting HISTORY.nc file in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync(comm)
              Call error(0,message)
           End If

           If (comm%idnode == 0) Then
              Call io_nc_get_real_precision( io_p, io_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(comm,safe)
           If (.not.safe) Then
              If (comm%idnode == 0) Write(message,'(a)') &
  "Can not determine the desired writing precision in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync(comm)
              Call error(0,message)
           End If

           If (comm%idnode == 0) safe = (io_p == file_p .and. io_r == file_r)
           Call gcheck(comm,safe)
           If (.not.safe) Then
              If (comm%idnode == 0) Then
                 Write(message,'(a)') &
  "Requested writing precision inconsistent with that in an existing HISTORY.nc"
                 Write(nrite, Fmt='(1x,a)', Advance='No') "Precision requested:"
                 Select Case( Selected_real_kind( io_p, io_r ) )
                 Case( Kind( 1.0 ) )
                    Write(message,'(1x,a)') "Single"
                 Case( Kind( 1.0d0 ) )
                    Write(message,'(1x,a)') "Double"
                 End Select
                 Write(nrite, Fmt='(1x,a)', Advance='No') "Precision in file  :"
                 Select Case( Selected_real_kind( file_p, file_r ) )
                 Case( Kind( 1.0 ) )
                    Write(message,'(1x,a)') "Single"
                 Case( Kind( 1.0d0 ) )
                    Write(message,'(1x,a)') "Double"
                 End Select
              End If

              Call gsync(comm)
              Call error(0,message)
           End If

! Get the frame number to check
! For netCDF this is the "frame number" which is not a long integer!

           jj=0
           If (comm%idnode == 0) Call io_nc_get_dim( 'frame', fh, jj )
           Call gbcast(comm,jj,0)

           If (jj > 0) Then
              frm=Int(jj,li)

           If (comm%idnode == 0) Call io_close( fh )
           Else ! Overwrite the file, it's junk to me
              lexist=.false.

              rec=Int(0,li)
              frm=Int(0,li)

              Go To 110
           End If
        End If

     End If
  End If

! Get offsets and define batch

  fail=0
  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Allocate (n_atm(0:comm%mxnode),        Stat=fail(1))
     Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure 0'
        Call error(0,message)
     End If

     chbat=' '
     n_atm=0 ; n_atm(comm%idnode+1)=natms
     Call gsum(comm,n_atm)
     n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Update frame

  frm=frm+Int(1,li)

! NO UNSORTED WRITING AS NO INDICES ARE WRITTEN
! In other words even if the user asked for unsorted I/O they get it sorted.

! SORTED MPI-I/O or Parallel Direct Access FORTRAN

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT   .or. &
      io_write == IO_WRITE_SORTED_NETCDF) Then

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, fname, mode_wronly, fh )

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,offset_kind)
     jj_io=rec_mpi_io
     jj=0 ! netCDF current frame
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_open( io_write, comm_self, fname, mode_wronly, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and cell information

           Write(record(1:recsz), Fmt='(a8,i6,f8.5,f12.5,a1)') 'timestep',nstep,tstep,time,lf
           Call io_write_record( fh, jj_io, record(1:recsz) )
           jj_io=jj_io + Int(1,offset_kind)

           Do i = 0, 2
              Write(record(1:recsz), Fmt='(3f10.3,a4,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 4 ), lf
              Call io_write_record( fh, jj_io, record(1:recsz) )
              jj_io=jj_io+Int(1,offset_kind)
           End Do

        Else ! netCDF write

! Get the current and new frame numbers

           Call io_nc_get_dim( 'frame', fh, jj )
           jj = jj + 1

           Call io_nc_put_var( 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var( 'timestep '      , fh,  tstep, jj, 1 )

           Call dcell(cell,celprp) ! get cell properties

           cell_vecs = Reshape( cell, (/ 3, 3 /) )

           lengths( 1 ) = celprp( 1 )
           lengths( 2 ) = celprp( 2 )
           lengths( 3 ) = celprp( 3 )

           angles ( 1 ) = Acos( celprp( 5 ) )
           angles ( 2 ) = Acos( celprp( 6 ) )
           angles ( 3 ) = Acos( celprp( 4 ) )
           angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

! Print

           Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, jj /), (/ 3, 3, 1 /) )
           Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, jj /), (/    3, 1 /) )
           Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, jj /), (/    3, 1 /) )

        End If

        Call io_close( fh )
        Call io_finalize

     End If
     Call gsync(comm)

! Start of file

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec_mpi_io=Int(rec,offset_kind)+Int(4,offset_kind)
     Else ! netCDF write
        Call gbcast(comm,jj,0)
        rec_mpi_io = Int(jj,offset_kind)
     End If

! Write the rest

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, fname, mode_wronly, fh )

     Call io_write_sorted_file( fh, 0*keytrj, IO_HISTORD, rec_mpi_io, natms, &
          ltg, atmnam,  (/ 0.0_wp /),  (/ 0.0_wp /), rsd, xxx, yyy, zzz,     &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                        &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),  ierr )

     If ( ierr /= 0 ) Then
        Select Case( ierr )
        Case( IO_BASE_COMM_NOT_SET )
           Call error( 1050 )
        Case( IO_ALLOCATION_ERROR )
           Call error( 1053 )
        Case( IO_UNKNOWN_WRITE_OPTION )
           Call error( 1056 )
        Case( IO_UNKNOWN_WRITE_LEVEL )
           Call error( 1059 )
        End Select
     End If

! Update and save offset pointer

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec=rec+Int(4,li)+Int(megatm,li)
        If (comm%idnode == 0) Then
           Write(record(1:recsz), Fmt='(2i2,3i10,a1)') keytrj,imcon,megatm,frm,rec,lf
           Call io_write_record( fh, Int(1,offset_kind), record(1:recsz) )
        End If
     End If

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (fff(1:mxatms),                             Stat=fail(3))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O

     If (comm%idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,li)
        Write(Unit=nhist, Fmt='(a8,i6,f8.5,f12.5,a1)', Rec=rec) 'timestep',nstep,tstep,time,lf

        Do i = 0, 2
           rec=rec+Int(1,li)
           Write(Unit=nhist, Fmt='(3f10.3,a4,a1)', Rec=rec) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 4 ), lf
        End Do

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           chbuf(i)=atmnam(i)
           fff(i)=rsd(i)
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,Traject_tag)

              Call grecv(comm,jatms,jdnode,Traject_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,fff(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,axx(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,ayy(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,azz(1:jatms),jdnode,Traject_tag)
              End If
           End If

           Do i=1,jatms
              rec1=rec+Int(iwrk(i),li)
              Write(Unit=nhist, Fmt='(a6,4f7.1,a1)', Rec=rec1) chbuf(i),axx(i),ayy(i),azz(i),fff(i),lf
           End Do
        End Do

! Update main header

        rec=rec+Int(megatm,li)
        Write(Unit=nhist, Fmt='(2i2,3i10,a1)', Rec=Int(2,li)) keytrj,imcon,megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Call grecv(comm,ready,0,Traject_tag)

        Call gsend(comm,natms,0,Traject_tag)
        If (natms > 0) Then
           Call gsend(comm,atmnam(:),0,Traject_tag)
           Call gsend(comm,ltg(:),0,Traject_tag)

           Call gsend(comm,rsd(:),0,Traject_tag)

           Call gsend(comm,xxx(:),0,Traject_tag)
           Call gsend(comm,yyy(:),0,Traject_tag)
           Call gsend(comm,zzz(:),0,Traject_tag)
        End If

! Save offset pointer

        rec=rec+Int(4,li)+Int(megatm,li)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (fff,         Stat=fail(3))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write deallocation failure'
        Call error(0,message)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write deallocation failure 0'
        Call error(0,message)
     End If
  End If

  Call gsync(comm)

End Subroutine trajectory_write

End Module trajectory
