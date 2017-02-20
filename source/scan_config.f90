Subroutine scan_config(megatm,imc_n,dvar,cfgname,levcfg,imcon,cell,xhi,yhi,zhi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of configuration file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
! contrib   - i.j.bush april 2010
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : nconf,config
  Use parse_module, Only : get_line, get_word, strip_blanks, word_2_real
  Use io_module,    Only : io_set_parameters,     &
                           io_get_parameters,     &
                           io_init, io_open,      &
                           io_nc_get_att,         &
                           io_nc_get_var,         &
                           io_close, io_finalize, &
                           IO_READ_MASTER,        &
                           IO_READ_NETCDF

  Implicit None

  Integer,              Intent( In    ) :: megatm,imc_n
  Real( Kind = wp ),    Intent( In    ) :: dvar
  Character( Len = * ), Intent(   Out ) :: cfgname
  Integer,              Intent(   Out ) :: levcfg,imcon
  Real( Kind = wp ),    Intent(   Out ) :: cell(1:9),xhi,yhi,zhi

  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,fname
  Logical                :: safe  = .true.  , &
                            l_ind = .false. , &
                            l_str = .false. , &
                            l_his = .false. , &
                            l_xtr = .true.  , &
                            fast
  Integer                :: totatm,i,j
  Real( Kind = wp )      :: xxx,yyy,zzz,buffer(1:4),cell_vecs(1:3,1:3)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: recsz = 73 ! default record size
  Integer                           :: fh, io_read
  Integer( Kind = MPI_OFFSET_KIND ) :: top_skip


! Get type of I/O for reading

  Call io_get_parameters( user_method_read = io_read )

! Define filename ASCII or netCDF

  If (io_read /= IO_READ_NETCDF) Then
     fname=Trim(config)
  Else
     fname=Trim(config)//'nc'
  End If

! Check if we have a CONFIG

  If (idnode == 0) Inquire(File=fname, Exist=safe)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(124)

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

     fast = .true.
     If (idnode == 0) Then

! Open CONFIG

        Open(Unit=nconf, File=fname)

! Read the CONFIG file header (TITLE record)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
           i = i + 1
           safe = .false.
           Read(Unit=nconf, Fmt='(a1)', Advance='No', IOStat=j, End=10) record(i:i)
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
           Read(Unit=nconf, Fmt='(a1)', Advance='No', IOStat=j, End=20) record(i:i)
           safe = .true.
           If (j < 0) Go To 20
        End Do
20      Continue
        fast = (fast .and. i == recsz)

! Read particles total value

        Call get_word(record,word) ; Call get_word(record,word)
        Call get_word(record,word) ; i=Nint(word_2_real(word,0.0_wp,l_str))
        fast = (fast .and. i == megatm)

     End If
     If (mxnode > 1) Then
        Call gsync()
        Call gcheck(safe,"enforce")
        Call gcheck(fast,"enforce")
     End If
     If (.not.safe) Go To 50

! Close CONFIG

     If (idnode == 0) Close(Unit=nconf)

  End If

!!! SCAN HEADER

  If (io_read /= IO_READ_NETCDF) Then ! ASCII read

! Open CONFIG

     If (idnode == 0) Open(Unit=nconf, File=fname)

! Read TITLE record (file header)

     Call get_line(safe,nconf,record)
     If (.not.safe) Go To 50

     Call strip_blanks(record)
     cfgname=record

! Read configuration level and image condition

     Call get_line(safe,nconf,record)
     If (.not.safe) Go To 50

     Call get_word(record,word)
     levcfg=Nint(word_2_real(word))

! halt execution if configuration level is unsupported

     If (levcfg < 0 .or. levcfg > 2) Call error(517)

     Call get_word(record,word)
     imcon=Nint(word_2_real(word))

! halt execution if image conventions is unsupported

     If (imcon < 0 .or. imcon > 7) Call error(514)

! specify MD cell (not defined for imcon=0)

     If (imcon /= 0) Then
        Call get_line(safe,nconf,record)
        If (.not.safe) Go To 50
        Call get_word(record,word)
        cell(1)=word_2_real(word)
        Call get_word(record,word)
        cell(2)=word_2_real(word)
        Call get_word(record,word)
        cell(3)=word_2_real(word)

        Call get_line(safe,nconf,record)
        If (.not.safe) Go To 50
        Call get_word(record,word)
        cell(4)=word_2_real(word)
        Call get_word(record,word)
        cell(5)=word_2_real(word)
        Call get_word(record,word)
        cell(6)=word_2_real(word)

        Call get_line(safe,nconf,record)
        If (.not.safe) Go To 50
        Call get_word(record,word)
        cell(7)=word_2_real(word)
        Call get_word(record,word)
        cell(8)=word_2_real(word)
        Call get_word(record,word)
        cell(9)=word_2_real(word)
     End If

! Close CONFIG

     If (idnode == 0) Close(Unit=nconf)
     If (mxnode > 1) Call gsync()

  Else ! netCDF read

! Open CONFIG

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )

     i=1 ! For config there is only one frame

     Call io_nc_get_att( 'title'          , fh, cfgname )

     Call io_nc_get_var( 'datalevel'      , fh, levcfg, i, 1  )
     If (levcfg < 0 .or. levcfg > 2) Call error(517)

     Call io_nc_get_var( 'imageconvention', fh,  imcon, i, 1  )
     If (imcon < 0 .or. imcon > 7) Call error(514)

     Call io_nc_get_var( 'cell'           , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
     cell = Reshape( cell_vecs, (/ Size( cell ) /) )

! Close CONFIG

     Call io_close( fh )

  End If

! ELABORATE SCAN

  totatm = 0
  xhi = 0.0_wp
  yhi = 0.0_wp
  zhi = 0.0_wp

  If (imcon == 0 .or. imcon == 6 .or. imc_n == 6) Then

! If MASTER read

     If (io_read == IO_READ_MASTER) Then

        If (idnode == 0) Then

! Open CONFIG

           Open(Unit=nconf, File=fname)

! Skip the header (we know exists from the basic scan)

           Read(Unit=nconf, Fmt=*)
           Read(Unit=nconf, Fmt=*)
           If (imcon /= 0) Then
              Read(Unit=nconf, Fmt=*)
              Read(Unit=nconf, Fmt=*)
              Read(Unit=nconf, Fmt=*)
           End If

! Find the extreme dimensions for the system

           Do
              Read(Unit=nconf, Fmt=*, End=40)
              Read(Unit=nconf, Fmt=*, End=30) xxx,yyy,zzz

              If (levcfg > 0) Then
                 Read(Unit=nconf, Fmt=*, End=30)
                 If (levcfg > 1) Read(Unit=nconf, Fmt=*, End=30)
              End If

              totatm=totatm+1
              xhi=Max(xhi,Abs(xxx))
              yhi=Max(yhi,Abs(yyy))
              zhi=Max(zhi,Abs(zzz))
           End Do

30         Continue ! catch error
           safe = .false.

        End If

40      Continue ! catch EoF
        If (mxnode > 1) Then
           Call gsync()
           Call gcheck(safe,"enforce")
        End If
        If (.not.safe) Go To 50

! Close CONFIG

        If (idnode == 0) Close(Unit=nconf)

        If (mxnode > 1) Then
           buffer(1)=xhi
           buffer(2)=yhi
           buffer(3)=zhi
           buffer(4)=Real(totatm,wp)
           Call gsum(buffer)
           xhi=buffer(1)
           yhi=buffer(2)
           zhi=buffer(3)
           totatm=Nint(buffer(4))
           If (totatm /= megatm) Call error(58)
        End If

! If PROPER read

     Else

! Open CONFIG

        If (fast) Then
           Call io_set_parameters( user_comm = dlp_comm_world )
           Call io_init( recsz )
           Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )
        Else
           Open(Unit=nconf, File=fname)
        End If

! top_skip is header size

        If (io_read /= IO_READ_NETCDF) Then
           If (imcon == 0) Then
              top_skip = Int(2,MPI_OFFSET_KIND)
           Else
              top_skip = Int(5,MPI_OFFSET_KIND)
           End If
        Else
           top_skip = Int(1,MPI_OFFSET_KIND) ! This is now the frame = 1
        End If

        Call read_config_parallel               &
           (levcfg, dvar, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi)

! Close CONFIG

        If (fast) Then
           Call io_close( fh )
           Call io_finalize
        Else
           Close(Unit=nconf)
        End If

     End If

  End If

  Return

! error exit for CONFIG file read

50 Continue
  If (idnode == 0) Close(Unit=nconf)
  Call error(55)

End Subroutine scan_config
