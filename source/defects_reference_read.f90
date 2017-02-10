Subroutine defects_reference_read(name,nstep,celr,nrefs,namr,indr,xr,yr,zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading particles data from REFERENCE file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,   Only : nrite,nrefdt,mxatms,half_minus,config
  Use site_module
  Use config_module,  Only : imcon,cell,natms
  Use domains_module, Only : nprx,npry,nprz,nprx_r,npry_r,nprz_r
  Use parse_module,   Only : tabs_2_blanks, get_line, get_word, word_2_real, strip_blanks
  Use io_module,      Only : io_set_parameters, io_get_parameters, &
                             io_init, io_open,                     &
                             io_nc_get_var,                        &
                             io_close, io_finalize,                &
                             IO_READ_MASTER, IO_READ_NETCDF

  Implicit None

  Character( Len = * ), Intent( In    ) :: name
  Integer,              Intent( In    ) :: nstep
  Real( Kind = wp ),    Intent(   Out ) :: celr(1:9)
  Character( Len = 8 ), Intent(   Out ) :: namr(1:mxatms)
  Integer,              Intent(   Out ) :: indr(1:mxatms),nrefs
  Real( Kind = wp ),    Intent(   Out ) :: xr(1:mxatms),yr(1:mxatms),zr(1:mxatms)

  Logical                :: l_ind = .true.  , &
                            l_str = .false. , &
                            lexist,fast,safe,loop,match
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,fname
  Integer                :: fail(1:2),megref,i,j,lvcfgr,imconr, &
                            indatm,idm,ipx,ipy,ipz,             &
                            itmols,isite,nsite,msite,fsite,ifrz,nrept
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz,cell_vecs(1:3,1:3)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: recsz = 73 ! default record size
  Integer                           :: fh, io_read
  Integer( Kind = MPI_OFFSET_KIND ) :: top_skip

  Character( Len = 8 ), Dimension( : ), Allocatable :: chbuf
  Integer,              Dimension( : ), Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ), Allocatable :: axx,ayy,azz


5 Continue

! Get type of I/O for reading

  Call io_get_parameters( user_method_read = io_read )

! Define filename ASCII or netCDF

  If (io_read /= IO_READ_NETCDF) Then
     fname=name
  Else
     fname=name // '.nc'
  End If

! Check if we have a REFERENCE and default megref

  lexist=.true. ; megref=0
  If (idnode == 0) Inquire(File=fname, Exist=lexist)
  If (mxnode > 1) Call gcheck(lexist,"enforce")
  If (.not.lexist) Then
     If (nstep > 1) Then ! That is a problem.  Abort, we must have restarted
        Call error(551)
     Else                ! Use data from CONFIG
        Call warning(320,0.0_wp,0.0_wp,0.0_wp)
        If (io_read /= IO_READ_NETCDF) Then
           fname=trim(config)
        Else
           fname=trim(config)//'.nc'
        End If
        megref=natms
        If (mxnode > 1) Call gsum(megref)
        If (imcon == 0) Call error(552) ! Lattice parameters are a must
     End If
  End If

! Define/Detect the FAST reading status

  If      (io_read == IO_READ_MASTER) Then

     fast = .false.

  Else If (io_read == IO_READ_NETCDF) Then

     fast = .true.

  Else

! Check if the system input file is a new style CONFIG:
! (i)  all lines are 72 ASCII characters long with
!      a UNIX carriage return as end of line;
! (ii) LINE2 has the particles total value
!      after values of levcfg and imcon.
! No fall back if users have mangled with further lines

     safe = .true.
     fast = .true.
     If (idnode == 0) Then

! Open REFERENCE

        Open(Unit=nrefdt, File=fname)

! Read the REFERENCE file header (TITLE record)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
           i = i + 1
           safe = .false.
           Read(Unit=nrefdt, Fmt='(a1)', Advance='No', IOStat=j, End=10) record(i:i)
           safe = .true.
           If (j < 0) Go To 10
        End Do
10      Continue
        fast = (fast .and. i == recsz)

! Read configuration level and image condition (RECORD2)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
           i = i + 1
           safe = .false.
           Read(Unit=nrefdt, Fmt='(a1)', Advance='No', IOStat=j, End=20) record(i:i)
           safe = .true.
           If (j < 0) Go To 20
        End Do
20      Continue
        fast = (fast .and. i == recsz)

! Read particles total value

        Call get_word(record,word) ; Call get_word(record,word)
        Call get_word(record,word) ; i=Nint(word_2_real(word,0.0_wp,l_str))
        fast = (fast .and. i /= 0)

     End If
     If (mxnode > 1) Then
        Call gsync()
        Call gcheck(safe,"enforce")
        Call gcheck(fast,"enforce")
     End If
     If (.not.safe) Go To 100

! Close REFERENCE

     If (idnode == 0) Close(Unit=nrefdt)

     If      (fast       ) Then
        If (megref == 0) Then          ! Define megref
           If (idnode == 0) megref = i ! It's already been read on master
           If (mxnode > 1) Call gsum(megref)
        End If
     Else If (megref == 0) Then
        io_read = IO_READ_MASTER       ! Abort parallel reading
     End If

  End If

!!! SCAN HEADER

  If (io_read /= IO_READ_NETCDF) Then ! ASCII read

! Open file

     If (idnode == 0) Open(Unit=nrefdt, File=fname)

! Read TITLE record (file header)

     Call get_line(safe,nrefdt,record)
     If (.not.safe) Go To 100

! Read configuration level and image condition

     Call get_line(safe,nrefdt,record)
     If (.not.safe) Go To 100

     Call get_word(record,word)
     lvcfgr=Nint(word_2_real(word))

     Call get_word(record,word)
     imconr=Nint(word_2_real(word))

     If (imconr /= 0) Then
        Call get_line(safe,nrefdt,record)
        If (.not.safe) Go To 100
        Call get_word(record,word)
        celr(1)=word_2_real(word)
        Call get_word(record,word)
        celr(2)=word_2_real(word)
        Call get_word(record,word)
        celr(3)=word_2_real(word)

        Call get_line(safe,nrefdt,record)
        If (.not.safe) Go To 100
        Call get_word(record,word)
        celr(4)=word_2_real(word)
        Call get_word(record,word)
        celr(5)=word_2_real(word)
        Call get_word(record,word)
        celr(6)=word_2_real(word)

        Call get_line(safe,nrefdt,record)
        If (.not.safe) Go To 100
        Call get_word(record,word)
        celr(7)=word_2_real(word)
        Call get_word(record,word)
        celr(8)=word_2_real(word)
        Call get_word(record,word)
        celr(9)=word_2_real(word)
     Else
        Call error(552) ! Lattice parameters are a must
     End If

! image conditions not compliant with DD and link-cell

     If (imconr == 4 .or. imconr == 5 .or. imconr == 7) Call error(300)

! Close REFERENCE

     If (idnode == 0) Close(Unit=nrefdt)
     If (mxnode > 1) Call gsync()

  Else ! netCDF read

! Open file

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )

     i=1 ! For config there is only one frame

     Call io_nc_get_var( 'datalevel'      , fh, lvcfgr, i, 1  )

     Call io_nc_get_var( 'imageconvention', fh, imconr, i, 1  )

! image conditions not compliant with DD and link-cell

     If (imconr == 4 .or. imconr == 5 .or. imconr == 7) Call error(300)

     Call io_nc_get_var( 'cell'           , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
     celr = Reshape( cell_vecs, (/ Size( celr ) /) )

! Close REFERENCE

     Call io_close( fh )

  End If

! REFERENCE to CONFIG cell match

  match=.true.
  If (.not.lexist) Then
     Do i=1,9
        match = match .and. Abs(cell(i)-celr(i)) < 1.0e-6_wp
     End Do
  End If

! Get lattice invert

  Call invert(celr,rcell,det)

! If MASTER read

  If (io_read == IO_READ_MASTER) Then

     fail=0
     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_read allocation failure, node: ', idnode
        Call error(0)
     End If

! Open file and skip header

     If (idnode == 0) Then
        Open(Unit=nrefdt, File=fname)

        Read(Unit=nrefdt, Fmt=*) ! REFERENCE file header (TITLE record)
        Read(Unit=nrefdt, Fmt=*) ! configuration level and image condition

        Read(Unit=nrefdt, Fmt=*) ! cell vectors
        Read(Unit=nrefdt, Fmt=*)
        Read(Unit=nrefdt, Fmt=*)
     End If

! Initialise domain localised referent atoms counter

     nrefs=0

! Initialise dispatched atom counter

     indatm=0

! Initialise total number of atoms and index counter

     megref=0
     j=0

     safe=.true.
     loop=.true.
     Do While (loop)

! Read in transmission arrays

        If (idnode == 0 .and. safe) Then
           record=' '
           Read(Unit=nrefdt, Fmt='(a)', End=40) record
           Call tabs_2_blanks(record) ; Call strip_blanks(record)
           Call get_word(record,word) ; chbuf(indatm+1)=word(1:8)
           Call get_word(record,word) ; iwrk(indatm+1) =Nint(word_2_real(word))
           If (iwrk(indatm+1) == 0) iwrk(indatm+1)=megref+1

           j=Max(j,iwrk(indatm+1))

           Read(Unit=nrefdt, Fmt=*, End=30) axx(indatm+1),ayy(indatm+1),azz(indatm+1)
           If (lvcfgr > 0) Read(Unit=nrefdt, Fmt=*, End=30)
           If (lvcfgr > 1) Read(Unit=nrefdt, Fmt=*, End=30)

! Escape End of File loop closure when all is fine

           Go To 50

! When things have gone wrong then

30         Continue
           safe=.false.

! When End Of File is encountered then

40         Continue
           loop=.false.

50         Continue
        End If

! checks

        If (mxnode > 1) Then
            Call gcheck(safe,"enforce")
            Call gcheck(loop,"enforce")
        End If
        If (.not.safe) Go To 100 ! catch error
        If (loop) Then

! increase counters

           indatm=indatm+1
           megref=megref+1

        Else

! Close file

           If (idnode == 0) Close(Unit=nrefdt)

! Check for inconsistencies in REFERENCE

           If (mxnode > 1) Call gmax(j)
           If (j /= megref) Call error(553)

        End If

! Circulate configuration data to all nodes when transmission arrays
! are filled up or this is the last looping

        If (indatm == mxatms .or. ((.not.loop) .and. indatm > 0)) Then

! Ensure all atoms are in prescribed simulation cell (DD bound) and broadcast them
!
!           Call pbcshift(imconr,celr,indatm,axx,ayy,azz)

           If (mxnode > 1) Then
              Call MPI_BCAST(chbuf,indatm*8,MPI_CHARACTER,0,dlp_comm_world,ierr)
              Call MPI_BCAST(iwrk,indatm,MPI_INTEGER,0,dlp_comm_world,ierr)

              Call MPI_BCAST(axx,indatm,wp_mpi,0,dlp_comm_world,ierr)
              Call MPI_BCAST(ayy,indatm,wp_mpi,0,dlp_comm_world,ierr)
              Call MPI_BCAST(azz,indatm,wp_mpi,0,dlp_comm_world,ierr)
           End If

! Assign atoms positions in fractional coordinates to the correct domains
! (DD bounding)

           Do i=1,indatm
              sxx=rcell(1)*axx(i)+rcell(4)*ayy(i)+rcell(7)*azz(i)
              syy=rcell(2)*axx(i)+rcell(5)*ayy(i)+rcell(8)*azz(i)
              szz=rcell(3)*axx(i)+rcell(6)*ayy(i)+rcell(9)*azz(i)

! sxx,syy,szz are in [-0.5,0.5) interval as values as 0.4(9) may pose a problem

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
                 Call error(555)
              Else If (idm == idnode)                 Then
                 nrefs=nrefs+1

                 If (nrefs < mxatms) Then
                    namr(nrefs)=chbuf(i)
                    indr(nrefs)=iwrk(i)

                    xr(nrefs)=sxx
                    yr(nrefs)=syy
                    zr(nrefs)=szz
                 Else
                    safe=.false.
                 End If
              End If
           End Do

! Check if all is dispatched fine

           If (mxnode > 1) Call gcheck(safe,"enforce")
           If (.not.safe) Call error(556)

! Nullify dispatch counter

           indatm=0

        End If

     End Do

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_read deallocation failure, node: ', idnode
        Call error(0)
     End If

! If PROPER read

  Else

! Open file

     If (fast) Then
        Call io_set_parameters( user_comm = dlp_comm_world )
        Call io_init( recsz )
        Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )
     Else
        Open(Unit=nrefdt, File=fname)
     End If

! top_skip is header size

     If (io_read /= IO_READ_NETCDF) Then
        top_skip = Int(5,MPI_OFFSET_KIND) ! imcon is a must
     Else
        top_skip = Int(1,MPI_OFFSET_KIND) ! This is now the frame = 1
     End If

     Call defects_reference_read_parallel       &
           (lvcfgr, celr, l_ind, l_str, megref, &
            fast, fh, top_skip, nrefs, namr, indr, xr, yr, zr)

! Close REFERENCE

     If (fast) Then
        Call io_close( fh )
        Call io_finalize
     Else
        Close(Unit=nrefdt)
     End If

  End If

! Remove frozen sites so they don't come up as vacancies
! only when dealing with CONFIG

  If (trim(fname) /= trim(config)) Then

     nsite=0
     msite=0
     fsite=0
     Do itmols=1,ntpmls
        Do isite=1,numsit(itmols)
           nsite=nsite+1

           If (frzsit(nsite) /= 0) Then
              Do nrept=1,nummols(itmols)
                 ifrz=nsite+msite+(nrept-1)*numsit(itmols)

                 Do i=1,nrefs
                    If (indr(i) == ifrz) Then
                       namr(i)=namr(nrefs) ; namr(nrefs)=' '
                       indr(i)=indr(nrefs) ; indr(nrefs)=0

                       xr(i)=xr(nrefs) ; xr(nrefs)=0.0_wp
                       yr(i)=yr(nrefs) ; yr(nrefs)=0.0_wp
                       zr(i)=zr(nrefs) ; zr(nrefs)=0.0_wp

                       nrefs=nrefs-1
                    End If
                 End Do
              End Do
           End If
        End Do

        msite=msite+(nummols(itmols)-1)*numsit(itmols)
        fsite=fsite+nummols(itmols)*numfrz(itmols)
     End Do

     If (fsite > 0) Then
        nsite=nrefs
        If (mxnode > 1) Call gsum(nsite)
        megref=nsite
     End If
  End If

! MATCH glitch fix

  If (.not.match) Then
     Call defects_reference_write(fname,megref,nrefs,namr,indr,xr,yr,zr)
     Go To 5
  End If

  Return

! REFERENCE format failure

100 Continue
  If (idnode == 0) Close(Unit=nrefdt)
  Call error(554)

End Subroutine defects_reference_read
