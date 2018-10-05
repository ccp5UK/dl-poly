Module trajectory
  Use kinds,         Only : wi,wp,li
  Use comms,         Only : comms_type,Traject_tag,gsync,wp_mpi,gbcast,gcheck, &
                            gsum,gsend,grecv,offset_kind,mode_wronly, &
                            mode_rdonly, comm_self
  Use domains,       Only : domains_type
  Use site, Only : site_type
  Use constants
  Use parse,         Only : tabs_2_blanks, get_line, get_word, &
                            strip_blanks, word_2_real
  Use configuration, Only : configuration_type
  Use netcdf_wrap,   Only : netcdf_param
  Use io,            Only : io_set_parameters,io_type,             &
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
                            IO_READ_MPIIO,                 &
                            IO_READ_DIRECT,                &
                            IO_READ_NETCDF,                &
                            IO_READ_MASTER,                &
                            IO_SUBSET_POSITIONS
  Use numerics,        Only : dcell, invert, shellsort2
  Use configuration,   Only : read_config_parallel
  Use errors_warnings, Only : error,warning,info
  Use particle,        Only : corePart
  Use filename, Only : file_type,FILE_CONFIG,FILE_HISTORY
  Implicit None

  Private

  !> Trajectory data
  Type, Public :: trajectory_type
    Private
    !> Level of detail in trajectory
    Integer( Kind = wi ), Public :: key
    !> Frequency to write trajectory
    Integer( Kind = wi ), Public :: freq
    !> Step to start writing trajectory
    Integer( Kind = wi ), Public :: start

    !> Contribution to record size, replaces egregious old method of calculating
    !> traj%key+2 repeatedly
    Integer( Kind = wi ) :: record_size
    Logical               :: newjob_read = .true.  , &
      l_ind_read  = .true.  , &
      l_his_read  = .true.  , &
      l_xtr_read  = .false. , &
      fast_read   = .true.
    Integer(Kind=li)     :: rec_read  = 0_li , &
      frm_read  = 0_li , &
      frm1_read = 0_li
    Integer                           :: recsz_read ! record size
    Integer                           :: fh_read, io_read
    Integer( Kind = offset_kind )     :: top_skip_read
    Character( Len = 1),  Allocatable :: buffer(:,:)
    Logical :: newjob_write = .true. , &
      fast_write   = .true.
    Character( Len = 40 ) :: fname
    Integer          :: recsz_write  = 73 ! default record size
    Integer(Kind=li) :: rec_write    = 0_li , &
      frm_write    = 0_li
    Logical, Public :: ltraj
    Logical, Public :: restart = .true.
  Contains
    !> Initialise the trajectory type
    Procedure, Public :: init => init_trajectory_type
    !> Return the number used in trajectory files corresponding to the level of
    !> detail, as these do not necessarily need to be identical
    Procedure, Public :: file_key => trajectory_file_key
    Final :: cleanup 
  End Type trajectory_type

  ! Trajectory detail level keys
  !> Coordinates only
  Integer( Kind = wi ), Parameter :: TRAJ_KEY_COORD = 1
  !> Coordinates and velocities
  Integer( Kind = wi ), Parameter :: TRAJ_KEY_COORD_VEL = 2
  !> Coordinates, velocities and forces
  Integer( Kind = wi ), Parameter :: TRAJ_KEY_COORD_VEL_FORCE = 3
  !> Compressed history file
  Integer( Kind = wi ), Parameter :: TRAJ_KEY_COMPRESSED = 4

  Public :: read_history
  Public :: trajectory_write
Contains
  !> cleanup the allocated arrays for trajectory variables
  Subroutine cleanup(T)
    Type(trajectory_type) :: T

    If (Allocated(T%buffer)) Then
      Deallocate(T%buffer)
    End If

  End Subroutine cleanup

  !> Initialise a trajectory type
  Subroutine init_trajectory_type(T,key,freq,start)
    Class(trajectory_type) :: T
    Integer( Kind = wi ), Intent( In ) :: key
    Integer( Kind = wi ), Intent( In ) :: freq
    Integer( Kind = wi ), Intent( In ) :: start

    ! Convert external trajectory key integers to internal parameters and
    ! determine the record size
    If (key == 0) Then
      T%key = TRAJ_KEY_COORD
      T%record_size = 2
    Else If (key == 1) Then
      T%key = TRAJ_KEY_COORD_VEL
      T%record_size = 3
    Else If (key == 2) Then
      T%key = TRAJ_KEY_COORD_VEL_FORCE
      T%record_size = 4
    Else If (key == 3) Then
      T%key = TRAJ_KEY_COMPRESSED
    End If

    T%freq = freq
    T%start = start
  End Subroutine init_trajectory_type

  Function trajectory_file_key(T) Result(key)
    Class(trajectory_type) :: T
    Integer( Kind = wi ) :: key

    If (T%key == TRAJ_KEY_COORD) Then
      key = 0
    Else If (T%key == TRAJ_KEY_COORD_VEL) Then
      key = 1
    Else If (T%key == TRAJ_KEY_COORD_VEL_FORCE) Then
      key = 2
    Else If (T%key == TRAJ_KEY_COMPRESSED) Then
      key = 3
    End If
  End Function trajectory_file_key

Subroutine read_history(l_str,fname,megatm,levcfg,dvar,nstep,tstep,time,exout, &
    io,traj,sites,domain,config,files,comm)

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

  Type( io_type ), Intent( InOut ) :: io
  Integer,              Intent( InOut ) :: levcfg,nstep
  Real( Kind = wp ),    Intent( In    ) :: dvar
  Real( Kind = wp ),    Intent( InOut ) :: tstep,time
  Integer,              Intent(   Out ) :: exout
  Type( trajectory_type), Intent( InOut ) :: traj
  Type( site_type ), Intent( In    ) :: sites
  Type( domains_type ), Intent( In    ) :: domain
  Type( configuration_type ),     Intent( InOut ) :: config
  Type( file_type ), Intent( InOut ) :: files(:)
  Type( comms_type),    Intent( InOut ) :: comm


  Logical                :: safe, lexist
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: fail(1:5),i,k,l,m,idm,ipx,ipy,ipz, &
                            indatm,nattot,totatm
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz,xhi,yhi,zhi, &
                            cell_vecs(1:3,1:3)

! Some parameters and variables needed by io interfaces


  Character( Len = 8 ), Allocatable :: chbuf(:)
  Integer,              Allocatable :: iwrk(:)
  Real( Kind = wp ),    Allocatable :: axx(:),ayy(:),azz(:), &
                                       bxx(:),byy(:),bzz(:), &
                                       cxx(:),cyy(:),czz(:)
  Integer  :: ierr
  Character ( Len = 256 )  :: message

  If (traj%newjob_read) Then
     traj%newjob_read = .false.

! Get type of I/O

     Call io_get_parameters(io, user_method_read = traj%io_read )

! ASCII read

     If (traj%io_read /= IO_READ_NETCDF) Then

! Define the default record size

        traj%recsz_read = 73

! Does HISTORY exist

        lexist=.true.
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
        If (.not.lexist) Go To 400

! Open HISTORY

        If (comm%idnode == 0) Open(Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)

! read the HISTORY file header

        Call get_line(safe,files(FILE_CONFIG)%unit_no,record,comm); If (.not.safe) Go To 300

        Call get_line(safe,files(FILE_CONFIG)%unit_no,record,comm); If (.not.safe) Go To 300
        Call get_word(record,word) ; levcfg=Nint(word_2_real(word,0.0_wp))
        Call get_word(record,word)
        Call get_word(record,word) ; If (Nint(word_2_real(word)) /= megatm) Go To 300
        Call get_word(record,word) ; traj%frm_read=Nint(word_2_real(word,0.0_wp),li)
        Call get_word(record,word) ; traj%rec_read=Nint(word_2_real(word,0.0_wp),li)

! Change traj%fast_read if no database records exists or the file is new

        If (traj%frm_read == Int(0,li) .and. traj%rec_read == Int(0,li)) traj%fast_read = .false.

        If (levcfg /= 3) Then

! Detect DL_POLY_2 HISTORY and amend traj%io_read and traj%recsz_read

           If ((.not.traj%fast_read) .and. traj%io_read == IO_READ_MPIIO) Then
              traj%io_read = IO_READ_DIRECT
              traj%recsz_read   = 200
           End If

        Else

! Detect compressed HISTORY and make amendments

           traj%l_ind_read = .false.
           If (traj%fast_read .and. traj%io_read /= IO_READ_MASTER) Then
              traj%io_read = IO_READ_MPIIO
              traj%recsz_read   = 35
           Else
              traj%io_read = IO_READ_DIRECT
           End If

        End If

        If (traj%io_read /= IO_READ_MASTER) Then
           traj%top_skip_read = Int(2,offset_kind)

           fail(1) = 0
           Allocate (traj%buffer(1:traj%recsz_read,1:4), Stat=fail(1))
           If (fail(1) > 0) Then
              Write(message,'(a)') 'read_history allocation failure 1'
              Call error(0,message)
           End If

           If (traj%io_read == IO_READ_MPIIO) Then
              Close(Unit=files(FILE_CONFIG)%unit_no)

              Call io_set_parameters(io, user_comm = comm%comm )
              Call io_init(io, traj%recsz_read )
              Call io_open(io, traj%io_read, comm%comm, fname, mode_rdonly, traj%fh_read )
           End If
        End If

     Else ! netCDF read

! Does HISTORY exist

        lexist=.true.
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
        If (.not.lexist) Go To 400

! fast and rec are irrelevant for netCDF (initialised at declaration)

        traj%fast_read = .true.
        traj%rec_read  = Int(0,li)

        Call io_set_parameters(io, user_comm = comm%comm )
        Call io_open(io, traj%io_read, comm%comm, fname, mode_rdonly, traj%fh_read )
        Call io_nc_get_dim(io, 'frame', traj%fh_read, i )

        If (i > 0) Then
          traj%frm_read=Int(i,li)
        Else
           Go To 300
        End If

     End If
  Else
     If (traj%io_read == IO_READ_MPIIO .or. traj%io_read == IO_READ_NETCDF) &
        Call io_set_parameters(io, user_comm = comm%comm )
     If (traj%io_read == IO_READ_MPIIO) Call io_init(io, traj%recsz_read )
  End If

! Reinitialise local-to-global counters and all levels of information at every bloody read
! Ilian's & Alin's fix to clear incoherencies with re-emptied/re-filled domains

  config%lsi=0 ; config%lsa=0 ; config%ltg=0                     ! A must unfortunately

! Necessary to not initialise keep the last frame info after hitting EOF
!  xxx=0.0_wp ; yyy = 0.0_wp ; zzz = 0.0_wp ! unfortunate but
!  vxx=0.0_wp ; vyy = 0.0_wp ; vzz = 0.0_wp ! unfortunate but
!  fxx=0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp ! unfortunate but

! MASTER READ

  If      (traj%io_read == IO_READ_MASTER) Then

     fail=0
     Allocate (chbuf(1:config%mxatms),                           Stat=fail(1))
     Allocate (iwrk(1:config%mxatms),                            Stat=fail(2))
     Allocate (axx(1:config%mxatms),ayy(1:config%mxatms),azz(1:config%mxatms), Stat=fail(3))
     Allocate (bxx(1:config%mxatms),byy(1:config%mxatms),bzz(1:config%mxatms), Stat=fail(4))
     Allocate (cxx(1:config%mxatms),cyy(1:config%mxatms),czz(1:config%mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'read_history allocation failure'
        Call error(0,message)
     End If

! read timestep and time

     Call get_line(safe,files(FILE_CONFIG)%unit_no,record,comm); If (.not.safe) Go To 200

     Do i=1,config%mxatms
       config%parts(i)%xxx=0.0_wp
       config%parts(i)%yyy=0.0_wp
       config%parts(i)%zzz=0.0_wp
       config%parts(i)%fxx=0.0_wp
       config%parts(i)%fyy=0.0_wp
       config%parts(i)%fzz=0.0_wp
     End Do
     config%vxx=0.0_wp ; config%vyy = 0.0_wp ; config%vzz = 0.0_wp

     Call get_word(record,word) ! timestep
     Call get_word(record,word) ; nstep = Nint(word_2_real(word))
     Call get_word(record,word) ; If (Nint(word_2_real(word)) /= megatm) Go To 300
     Call get_word(record,word) ; levcfg = Nint(word_2_real(word))
     Call get_word(record,word) ; config%imcon = Nint(word_2_real(word))

! image conditions not compliant with DD and link-cell

     If (config%imcon == 4 .or. config%imcon == 5 .or. config%imcon == 7) Call error(300)

     Call get_word(record,word) ; tstep = word_2_real(word)
     Call get_word(record,word) ; time = word_2_real(word)

     Write(message,'(a,i10,a,f10.3,a)') 'HISTORY step ',nstep,' (',time,' ps) is being read'
     Call info(message,.true.)

! read cell vectors

     Call get_line(safe,files(FILE_CONFIG)%unit_no,record,comm); If (.not.safe) Go To 300
     Call get_word(record,word); config%cell(1)=word_2_real(word)
     Call get_word(record,word); config%cell(2)=word_2_real(word)
     Call get_word(record,word); config%cell(3)=word_2_real(word)

     Call get_line(safe,files(FILE_CONFIG)%unit_no,record,comm); If (.not.safe) Go To 300
     Call get_word(record,word); config%cell(4)=word_2_real(word)
     Call get_word(record,word); config%cell(5)=word_2_real(word)
     Call get_word(record,word); config%cell(6)=word_2_real(word)

     Call get_line(safe,files(FILE_CONFIG)%unit_no,record,comm); If (.not.safe) Go To 300
     Call get_word(record,word); config%cell(7)=word_2_real(word)
     Call get_word(record,word); config%cell(8)=word_2_real(word)
     Call get_word(record,word); config%cell(9)=word_2_real(word)

     Call invert(config%cell,rcell,det)

! Initialise domain localised atom counter (configuration)

     config%natms=0

! Initialise dispatched atom counter

     indatm=0

! Initialise total number of atoms counter

     nattot=0

     Do k=1,sites%ntype_mol
        Do l=1,sites%num_mols(k)
           Do m=1,sites%num_site(k)

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
                 record=' '; Read(Unit=files(FILE_CONFIG)%unit_no, Fmt='(a)', End=30) record
                 Call tabs_2_blanks(record) ; Call strip_blanks(record)
                 Call get_word(record,word) ; chbuf(indatm)=word(1:8)
                 If (traj%l_ind_read) Then
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
                    Read(Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) axx(indatm),ayy(indatm),azz(indatm)

                    If (levcfg > 0) Then
                       Read(Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) bxx(indatm),byy(indatm),bzz(indatm)
                       If (levcfg > 1) Read(Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) cxx(indatm),cyy(indatm),czz(indatm)
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

              If (indatm == config%mxatms .or. nattot == megatm) Then

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

! Assign atoms to cortraj%rec_readt domains (DD bound)

                 Do i=1,indatm
                    sxx=rcell(1)*axx(i)+rcell(4)*ayy(i)+rcell(7)*azz(i)
                    syy=rcell(2)*axx(i)+rcell(5)*ayy(i)+rcell(8)*azz(i)
                    szz=rcell(3)*axx(i)+rcell(6)*ayy(i)+rcell(9)*azz(i)

! sxx,syy,szz are in [-0.5,0.5) interval and values as 0.4(9) may pose a problem

                    sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
                    syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
                    szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

! fold back coordinatesc

                    axx(i)=config%cell(1)*sxx+config%cell(4)*syy+config%cell(7)*szz
                    ayy(i)=config%cell(2)*sxx+config%cell(5)*syy+config%cell(8)*szz
                    azz(i)=config%cell(3)*sxx+config%cell(6)*syy+config%cell(9)*szz

! assign domain coordinates (call for errors)

                    ipx=Int((sxx+0.5_wp)*domain%nx_real)
                    ipy=Int((syy+0.5_wp)*domain%ny_real)
                    ipz=Int((szz+0.5_wp)*domain%nz_real)

                    idm=ipx+domain%nx*(ipy+domain%ny*ipz)
                    If      (idm < 0 .or. idm > (comm%mxnode-1)) Then
                       Call error(513)
                    Else If (idm == comm%idnode)                 Then
                       config%natms=config%natms+1

                       If (config%natms < config%mxatms) Then
                          config%atmnam(config%natms)=chbuf(i)
                          config%ltg(config%natms)=iwrk(i)

                          If (levcfg /= 3) Then
                             config%parts(config%natms)%xxx=axx(i)
                             config%parts(config%natms)%yyy=ayy(i)
                             config%parts(config%natms)%zzz=azz(i)
                             If (levcfg > 0) Then
                                config%vxx(config%natms)=bxx(i)
                                config%vyy(config%natms)=byy(i)
                                config%vzz(config%natms)=bzz(i)
                                If (levcfg > 1) Then
                                   config%parts(config%natms)%fxx=cxx(i)
                                   config%parts(config%natms)%fyy=cyy(i)
                                   config%parts(config%natms)%fzz=czz(i)
                                End If
                             End If
                          Else
                             config%parts(config%natms)%xxx=axx(i)
                             config%parts(config%natms)%yyy=ayy(i)
                             config%parts(config%natms)%zzz=azz(i)
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

  Else If (traj%io_read /= IO_READ_NETCDF) Then

     Call io_read_batch(io, traj%fh_read, traj%top_skip_read, 4, traj%buffer, ierr )
     If (ierr < 0) Go To 300
     traj%top_skip_read = traj%top_skip_read + Int(4,offset_kind)

     record = ' '
     Do i = 1, Min( Size( traj%buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = traj%buffer( i, 1 )
     End Do
     Call get_word(record,word) ! timestep
     Call get_word(record,word) ; nstep = Nint(word_2_real(word))
     Call get_word(record,word) ; If (Nint(word_2_real(word)) /= megatm) Go To 300
     Call get_word(record,word) ; levcfg = Nint(word_2_real(word))
     Call get_word(record,word) ; config%imcon = Nint(word_2_real(word))

! image conditions not compliant with DD and link-cell

    If (config%imcon == 4 .or. config%imcon == 5 .or. config%imcon == 7) Call error(300)

     Call get_word(record,word) ; tstep = word_2_real(word)
     Call get_word(record,word) ; time = word_2_real(word)

     Write(message,'(a,i10,a,f10.3,a)') 'HISTORY step ',nstep,' (',time,' ps) is being read'
     Call info(message,.true.)

! read cell vectors

     record = ' '
     Do i = 1, Min( Size( traj%buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = traj%buffer( i, 2 )
     End Do
     Call get_word(record,word); config%cell(1)=word_2_real(word)
     Call get_word(record,word); config%cell(2)=word_2_real(word)
     Call get_word(record,word); config%cell(3)=word_2_real(word)

     record = ' '
     Do i = 1, Min( Size( traj%buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = traj%buffer( i, 3 )
     End Do
     Call get_word(record,word); config%cell(4)=word_2_real(word)
     Call get_word(record,word); config%cell(5)=word_2_real(word)
     Call get_word(record,word); config%cell(6)=word_2_real(word)

     record = ' '
     Do i = 1, Min( Size( traj%buffer, Dim = 1 ) - 1, Len( record ) )
        record( i:i ) = traj%buffer( i, 4 )
     End Do
     Call get_word(record,word); config%cell(7)=word_2_real(word)
     Call get_word(record,word); config%cell(8)=word_2_real(word)
     Call get_word(record,word); config%cell(9)=word_2_real(word)

     Call read_config_parallel(config,levcfg,dvar,traj%l_ind_read,l_str,megatm, &
       traj%l_his_read,traj%l_xtr_read,traj%fast_read,traj%fh_read, &
       traj%top_skip_read,xhi,yhi,zhi,io,domain,files,comm)

     If (traj%fast_read) Then
        If (levcfg /= 3) Then
           i=levcfg+2
        Else
           i=levcfg-2
        End If

        traj%top_skip_read = traj%top_skip_read + Int(i,offset_kind)*Int(megatm,offset_kind)
     Else
        traj%top_skip_read = Int(0,offset_kind)
     End If

! Update current frame and exit gracefully

     traj%frm1_read = traj%frm1_read+Int(1,li)
     If (traj%frm1_read == traj%frm_read) Go To 200

  Else ! netCDF read

! Update current frame

     traj%frm1_read = traj%frm1_read+Int(1,li)
     i = Int( traj%frm1_read )

     Call io_nc_get_var(io, 'time'           , traj%fh_read,   time, i, 1 )
     Call io_nc_get_var(io, 'datalevel'      , traj%fh_read, levcfg, i, 1 )
     Call io_nc_get_var(io, 'imageconvention', traj%fh_read,  config%imcon, i, 1 )

! image conditions not compliant with DD and link-cell

     If (config%imcon == 4 .or. config%imcon == 5 .or. config%imcon == 7) Call error(300)

     Call io_nc_get_var(io, 'timestep'       , traj%fh_read,  tstep, i, 1 )
     Call io_nc_get_var(io, 'step'           , traj%fh_read,  nstep, i, 1 )

     Write(message,'(a,i10,a,f10.3,a)') 'HISTORY step ',nstep,' (',time,' ps) is being read'
     Call info(message,.true.)

! Note that in netCDF the frames are not long integers - Int( traj%frm1_read )

     Call io_nc_get_var(io, 'cell'           , traj%fh_read, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
     config%cell = Reshape( cell_vecs, (/ Size( config%cell ) /) )

     Call read_config_parallel(config,levcfg,dvar,traj%l_ind_read,l_str,megatm, &
       traj%l_his_read,traj%l_xtr_read,traj%fast_read,traj%fh_read,  &
       Int(i,Kind(traj%top_skip_read)),xhi,yhi,zhi,io,domain,files,comm)

     If (traj%frm1_read == traj%frm_read) Go To 200

  End If

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  config%nlast=config%natms

! Does the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?
! Check number of atoms in system (CONFIG = FIELD)

  totatm=config%natms
  Call gsum(comm,totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (configuration)

  Do i=1,config%natms
     config%lsi(i)=i
     config%lsa(i)=config%ltg(i)
  End Do
  Call shellsort2(config%natms,config%lsi,config%lsa)

  exout = 0 ! more to read indicator

  Return

! Normal exit from HISTORY file read

200 Continue

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  config%nlast=config%natms

! Does the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?
! Check number of atoms in system (CONFIG = FIELD)

  totatm=config%natms
  Call gsum(comm,totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (configuration)

  Do i=1,config%natms
     config%lsi(i)=i
     config%lsa(i)=config%ltg(i)
  End Do
  Call shellsort2(config%natms,config%lsi,config%lsa)

  exout = 1 ! It's an indicator of the end of reading.

  Call info('HISTORY end of file reached',.true.)

  If (traj%io_read == IO_READ_MASTER) Then
     If (config%imcon == 0) Close(Unit=files(FILE_CONFIG)%unit_no)
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
     Call io_close(io, traj%fh_read )
     Call io_finalize(io)
  End If

  Return

! Abnormal exit from HISTORY file read

300 Continue

  Call info('HISTORY data mishmash detected',.true.)
  exout = -1 ! It's an indicator of the end of reading.
  If (traj%io_read == IO_READ_MASTER) Then
     If (config%imcon == 0) Close(Unit=files(FILE_CONFIG)%unit_no)
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
     Call io_close(io, traj%fh_read )
     Call io_finalize(io)
  End If

  Return

! HISTORY not found

400 Continue

  Call error(585)

End Subroutine read_history

Subroutine trajectory_write(keyres,nstep,tstep,time,io,rsd,netcdf,config, &
  traj,files,comm)

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


  Integer,           Intent( In    ) :: keyres,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep,time
  Real( Kind = wp ), Intent( In    ) :: rsd(:)
  Type( io_type ), Intent( InOut ) :: io
  Type( netcdf_param ), Intent( In    ) :: netcdf
  Type( configuration_type ),     Intent( InOut ) :: config
  Type( trajectory_type ), Intent( InOut ) :: traj
  Type( file_type ), Intent( InOut ) :: files(:)
  Type( comms_type ),   Intent( InOut ) :: comm

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
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: eee,fff
  Type( corePart  ),    Dimension( : ),    Allocatable :: temp_parts

! netCDF check

  Integer :: file_p, file_r
  Integer :: io_p, io_r,ierr
 
  Character ( Len = 256 )  ::  message

  If (.not.(nstep >= traj%start .and. Mod(nstep-traj%start,traj%freq) == 0)) Return

! Get write method, buffer size and line feed character

  Call io_get_parameters(io, user_method_write      = io_write )
  Call io_get_parameters(io, user_buffer_size_write = batsz    )
  Call io_get_parameters(io, user_line_feed         = lf       )

  If (traj%key == TRAJ_KEY_COMPRESSED) Then
     traj%recsz_write = 35 ! record size for compressed HISTORY file
     Go To 100
  End If

  If (traj%newjob_write) Then
     traj%newjob_write = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        traj%fname = files(FILE_HISTORY)%filename
     Else
        traj%fname = Trim(files(FILE_HISTORY)%filename) // '.nc'
     End If

! If keyres=1, is HISTORY old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File=traj%fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (comm%idnode == 0) Then
              Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, &
                Form='formatted', Access='direct', Status='replace', Recl=traj%recsz_write)
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a72,a1)',       Rec=Int(1,li)) config%cfgname(1:72),lf
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3i10,2i21,a1)', Rec=Int(2,li)) &
              traj%file_key(),config%imcon,config%megatm,traj%frm_write,traj%rec_write,lf
              Close(Unit=files(FILE_HISTORY)%unit_no)
           End If
           traj%rec_write=Int(2,li)
           traj%frm_write=Int(0,li)
        Else
           If (comm%idnode == 0) Then
              Call io_set_parameters(io, user_comm = comm_self )
              Call io_nc_create( netcdf, comm_self, traj%fname, config%cfgname, config%megatm )
           End If
        End If

! Get some sense out of it

     Else

        safe=.true.
        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           If (comm%idnode == 0) Then

              Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, Form='formatted')

              Do

                 record(1:traj%recsz_write)=' '

! Assume new style of HISTORY with bookkeeping.

                 If (traj%fast_write) Then

                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt=*, End=20)                     ! title record
                    traj%rec_write=traj%rec_write+Int(1,li)
                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a)', End=20) record(1:traj%recsz_write) ! bookkeeping record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    traj%rec_write=traj%rec_write+Int(1,li)

                    Call get_word(record(1:traj%recsz_write),word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record(1:traj%recsz_write),word) ; Call get_word(record(1:traj%recsz_write),word)
                       Call get_word(record(1:traj%recsz_write),word) ; traj%frm_write=Nint(word_2_real(word,0.0_wp),li)
                       Call get_word(record(1:traj%recsz_write),word) ; traj%rec_write=Nint(word_2_real(word,0.0_wp),li)
                       If (traj%frm_write /= Int(0,li) .and. traj%rec_write > Int(2,li)) Then
                          Go To 20 ! New style
                       Else
                          traj%fast_write=.false. ! TOUGH, old style
                          traj%rec_write=Int(2,li)
                          traj%frm_write=Int(0,li)
                       End If
                    Else
                       safe=.false. ! Overwrite the file, it's junk to me
                       Go To 20
                    End If

! TOUGH, it needs scanning through

                 Else

                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a)', End=20) record(1:traj%recsz_write) ! timestep record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    traj%rec_write=traj%rec_write+Int(1,li)

                    Call get_word(record(1:traj%recsz_write),word) ; Call get_word(record(1:traj%recsz_write),word)
                    Call get_word(record(1:traj%recsz_write),word) ; jj=Nint(word_2_real(word))
                    Call get_word(record(1:traj%recsz_write),word) ; k=Nint(word_2_real(word))

                    word=' '
                    i = 3 + (2+k)*jj ! total number of lines to read
                    Write(word,'( "(", i0, "( / ) )" )') i-1
                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt=word, End=20)
                    traj%rec_write=traj%rec_write+Int(i,li)
                    traj%frm_write=traj%frm_write+Int(1,li)

                 End If

              End Do

20            Continue
              Close(Unit=files(FILE_HISTORY)%unit_no)

           End If

           Call gcheck(comm,safe,"enforce")
           If (.not.safe) Then
              lexist=.false.

              traj%rec_write=Int(0,li)
              traj%frm_write=Int(0,li)

              Go To 10
           Else If (comm%mxnode > 1) Then
              buffer(1)=Real(traj%frm_write,wp)
              buffer(2)=Real(traj%rec_write,wp)

              Call gbcast(comm,buffer,0)

              traj%frm_write=Nint(buffer(1),li)
              traj%rec_write=Nint(buffer(2),li)
           End If

        Else ! netCDF read

           If (comm%idnode == 0) Then
              Call io_set_parameters(io, user_comm = comm_self )
              Call io_open(io, io_write, comm_self, traj%fname, mode_rdonly, fh )

! Get the precision that the history file was written in
! and check it matches the requested precision

              Call io_nc_get_file_real_precision(io, fh, file_p, file_r, ierr )
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
              Call io_nc_get_real_precision( netcdf, io_p, io_r, ierr )
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
              Select Case( Selected_real_kind( io_p, io_r ) )
              Case( Kind( 1.0 ) )
                 Write(message,'(a)') 'Precision requested: Single'
              Case( Kind( 1.0d0 ) )
                 Write(message,'(a)') 'Precision requested: Double'
              End Select
              Call info(message,.true.)
              Select Case( Selected_real_kind( file_p, file_r ) )
              Case( Kind( 1.0 ) )
                 Write(message,'(a)') "Precision in file: Single"
              Case( Kind( 1.0d0 ) )
                 Write(message,'(a)') "Precision in file: Double"
              End Select
              Call info(message,.true.)

              Call gsync(comm)
              Call error(0,'Requested writing precision inconsistent with that in an existing HISTORY.nc')
           End If

! Get the frame number to check
! For netCDF this is the "frame number" which is not a long integer!

           If (comm%idnode == 0) Call io_nc_get_dim(io, 'frame', fh, jj )
           Call gbcast(comm,jj,0)

           If (jj > 0) Then
              traj%frm_write=Int(jj,li)

              If (comm%idnode == 0) Call io_close(io, fh )
           Else ! Overwrite the file, it's junk to me
              lexist=.false.

              traj%rec_write=Int(0,li)
              traj%frm_write=Int(0,li)

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
     Allocate (chbat(1:traj%recsz_write,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure 0'
        Call error(0,message)
     End If

     chbat=' '
     n_atm=0 ; n_atm(comm%idnode+1)=config%natms
     Call gsum(comm,n_atm)
     n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = offset_kind)

! Update frame

  traj%frm_write=traj%frm_write+Int(1,li)

! UNSORTED MPI-I/O or Parallel Ditraj%rec_writet Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(traj%rec_write,offset_kind)
     jj=0
     If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm = comm_self )
        Call io_init(io, traj%recsz_write )
        Call io_open(io, io_write, comm_self, traj%fname, mode_wronly, fh )

        Write(record(1:traj%recsz_write), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,config%megatm,traj%file_key(),&
         config%imcon,tstep,time,lf
        jj=jj+1
        Do k=1,traj%recsz_write
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record(1:traj%recsz_write), Fmt='(3f20.10,a12,a1)') &
                config%cell( 1 + i * 3 ), config%cell( 2 + i * 3 ), config%cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,traj%recsz_write
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and cell information

        Call io_write_batch(io, fh, rec_mpi_io, jj, chbat )

        Call io_close(io, fh )
        Call io_finalize(io)

     Else

        jj=jj+4

     End If
     Call gsync(comm)

! Start of file

     rec_mpi_io=Int(traj%rec_write,offset_kind)+Int(jj,offset_kind)+Int(n_atm(0),offset_kind)*Int(traj%record_size,offset_kind)
     jj=0

     Call io_set_parameters(io, user_comm = comm%comm )
     Call io_init(io, traj%recsz_write )
     Call io_open(io, io_write, comm%comm, traj%fname, mode_wronly, fh )

     Do i=1,config%natms
        Write(record(1:traj%recsz_write), Fmt='(a8,i10,3f12.6,a18,a1)') config%atmnam(i),config%ltg(i),config%weight(i),&
                                                             config%parts(i)%chge,rsd(i),Repeat(' ',18),lf
        jj=jj+1
        Do k=1,traj%recsz_write
           chbat(k,jj) = record(k:k)
        End Do

        Write(record(1:traj%recsz_write), Fmt='(3g20.10,a12,a1)') config%parts(i)%xxx,config%parts(i)%yyy,&
                                                       config%parts(i)%zzz,Repeat(' ',12),lf
        jj=jj+1
        Do k=1,traj%recsz_write
           chbat(k,jj) = record(k:k)
        End Do

        If (traj%key /= TRAJ_KEY_COORD) Then
           Write(record(1:traj%recsz_write), Fmt='(3g20.10,a12,a1)') config%vxx(i),config%vyy(i),config%vzz(i),&
                                                          Repeat(' ',12),lf
           jj=jj+1
           Do k=1,traj%recsz_write
              chbat(k,jj) = record(k:k)
           End Do
        End If

        If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
           Write(record(1:traj%recsz_write), Fmt='(3g20.10,a12,a1)') config%parts(i)%fxx,config%parts(i)%fyy,&
                                                          config%parts(i)%fzz,Repeat(' ',12),lf
           jj=jj+1
           Do k=1,traj%recsz_write
              chbat(k,jj) = record(k:k)
           End Do
        End If

! Dump batch and update start of file

        If (jj + traj%record_size >= batsz .or. i == config%natms) Then
           Call io_write_batch(io, fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,offset_kind)
           jj=0
        End If
     End Do

! Update and save offset pointer

     traj%rec_write=traj%rec_write+Int(4,li)+Int(config%megatm,li)*Int(traj%record_size,li)
     If (comm%idnode == 0) Then
        Write(record(1:traj%recsz_write), Fmt='(3i10,2i21,a1)') traj%file_key(),config%imcon,config%megatm,traj%frm_write,&
           traj%rec_write,lf
        Call io_write_record(io, fh, Int(1,offset_kind), record(1:traj%recsz_write) )
     End If

     Call io_close(io, fh )
     Call io_finalize(io)

! UNSORTED Serial Ditraj%rec_writet Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:config%mxatms),iwrk(1:config%mxatms),                   Stat=fail(1))
     Allocate (bxx(1:config%mxatms),byy(1:config%mxatms),bzz(1:config%mxatms),        Stat=fail(3))
     Allocate (temp_parts(1:config%mxatms),eee(1:config%mxatms),fff(1:config%mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     jj=0
     If (comm%idnode == 0) Then

        Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, Form='formatted', Access='direct', Recl=traj%recsz_write)

! Accumulate header

        Write(record(1:traj%recsz_write), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,config%megatm,traj%file_key(),&
          config%imcon,tstep,time,lf
        jj=jj+1
        Do k=1,traj%recsz_write
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record(1:traj%recsz_write), Fmt='(3f20.10,a12,a1)') &
                config%cell( 1 + i * 3 ), config%cell( 2 + i * 3 ), config%cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,traj%recsz_write
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(73a)', Rec=traj%rec_write+Int(1,li)) (chbat(:,k), k=1,jj)
        traj%rec_write=traj%rec_write+Int(jj,li)
        jj=0

        Do i=1,config%natms
           iwrk(i)=config%ltg(i)

           temp_parts(i)%xxx=config%parts(i)%xxx
           temp_parts(i)%yyy=config%parts(i)%yyy
           temp_parts(i)%zzz=config%parts(i)%zzz

           chbuf(i)=config%atmnam(i)
           temp_parts(i)%chge=config%parts(i)%chge
           eee(i)=config%weight(i)
           fff(i)=rsd(i)
        End Do

        If (Any(traj%key == [TRAJ_KEY_COORD_VEL,TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
           Do i=1,config%natms
              bxx(i)=config%vxx(i)
              byy(i)=config%vyy(i)
              bzz(i)=config%vzz(i)
           End Do
        End If

        If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
           Do i=1,config%natms
              temp_parts(i)%fxx=config%parts(i)%fxx
              temp_parts(i)%fyy=config%parts(i)%fyy
              temp_parts(i)%fzz=config%parts(i)%fzz
           End Do
        End If

        jatms=config%natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,Traject_tag)

              Call grecv(comm,jatms,jdnode,Traject_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,temp_parts(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,eee(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,fff(1:jatms),jdnode,Traject_tag)

                 If (traj%key /= TRAJ_KEY_COORD) Then
                    Call grecv(comm,bxx(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,byy(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,bzz(1:jatms),jdnode,Traject_tag)
                 End If

!                 If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
!                    Call grecv(comm,cxx(1:jatms),jdnode,Traject_tag)
!                    Call grecv(comm,cyy(1:jatms),jdnode,Traject_tag)
!                    Call grecv(comm,czz(1:jatms),jdnode,Traject_tag)
!                 End If
              End If
           End If

           Do i=1,jatms
              Write(record(1:traj%recsz_write), Fmt='(a8,i10,3f12.6,a18,a1)') chbuf(i),iwrk(i),eee(i),temp_parts(i)%chge,&
                                                                   fff(i),Repeat(' ',18),lf
              jj=jj+1
              Do k=1,traj%recsz_write
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record(1:traj%recsz_write), Fmt='(3g20.10,a12,a1)') temp_parts(i)%xxx,temp_parts(i)%yyy,temp_parts(i)%zzz,&
                                                             Repeat(' ',12),lf
              jj=jj+1
              Do k=1,traj%recsz_write
                 chbat(k,jj) = record(k:k)
              End Do

              If (Any(traj%key == [TRAJ_KEY_COORD_VEL,TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
                 Write(record(1:traj%recsz_write), Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,traj%recsz_write
                    chbat(k,jj) = record(k:k)
                 End Do
              End If

              If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
                 Write(record(1:traj%recsz_write), Fmt='(3g20.10,a12,a1)') temp_parts(i)%fxx,temp_parts(i)%fyy,temp_parts(i)%fzz,&
                                                                Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,traj%recsz_write
                    chbat(k,jj) = record(k:k)
                 End Do
              End If

! Dump batch and update start of file

              If (jj + traj%record_size >= batsz .or. i == jatms) Then
                 Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(73a)', Rec=traj%rec_write+Int(1,li)) (chbat(:,k), k=1,jj)
                 traj%rec_write=traj%rec_write+Int(jj,li)
                 jj=0
              End If
           End Do
        End Do

! Update main header

        Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3i10,2i21,a1)', Rec=Int(2,li)) &
        traj%file_key(),config%imcon,config%megatm,traj%frm_write,traj%rec_write,lf

        Close(Unit=files(FILE_HISTORY)%unit_no)

     Else

        Call grecv(comm,ready,0,Traject_tag)

        Call gsend(comm,config%natms,0,Traject_tag)
        If (config%natms > 0) Then
           Call gsend(comm,config%atmnam(:),0,Traject_tag)
           Call gsend(comm,config%ltg(:),0,Traject_tag)

           Call gsend(comm,config%parts(:),0,Traject_tag)
           Call gsend(comm,config%weight(:),0,Traject_tag)
           Call gsend(comm,rsd(:),0,Traject_tag)


           If (traj%key /= TRAJ_KEY_COORD) Then
              Call gsend(comm,config%vxx(:),0,Traject_tag)
              Call gsend(comm,config%vyy(:),0,Traject_tag)
              Call gsend(comm,config%vzz(:),0,Traject_tag)
           End If

        End If

! Save offset pointer

        traj%rec_write=traj%rec_write+Int(4,li)+Int(config%megatm,li)*Int(traj%record_size,li)

     End If

     Deallocate (chbuf,iwrk,         Stat=fail(1))
     Deallocate (bxx,byy,bzz,        Stat=fail(3))
     Deallocate (temp_parts,eee,fff, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write deallocation failure'
        Call error(0,message)
     End If

! SORTED MPI-I/O or Parallel Ditraj%rec_writet Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(traj%rec_write,offset_kind)
     jj_io=rec_mpi_io
     jj=0 ! netCDF current frame
     If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm = comm_self )
        Call io_init(io, traj%recsz_write )
        Call io_open(io, io_write, comm_self, traj%fname, mode_wronly, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and cell information

          Write(record(1:traj%recsz_write), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,config%megatm,traj%file_key(),&
           config%imcon,tstep,time,lf
           Call io_write_record(io, fh, jj_io, record(1:traj%recsz_write) )
           jj_io=jj_io + Int(1,offset_kind)

           Do i = 0, 2
              Write(record(1:traj%recsz_write), Fmt='(3f20.10,a12,a1)') &
                   config%cell( 1 + i * 3 ), config%cell( 2 + i * 3 ), config%cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              Call io_write_record(io, fh, jj_io, record(1:traj%recsz_write) )
              jj_io=jj_io+Int(1,offset_kind)
           End Do

        Else ! netCDF write

! Get the current and new frame numbers

           Call io_nc_get_dim(io, 'frame', fh, jj )
           jj = jj + 1

           Call io_nc_put_var(io, 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var(io, 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var(io, 'datalevel'      , fh, traj%file_key(), jj, 1 )
           Call io_nc_put_var(io, 'imageconvention', fh,  config%imcon, jj, 1 )
           Call io_nc_put_var(io, 'timestep '      , fh,  tstep, jj, 1 )

           Call dcell(config%cell,celprp) ! get config%cell properties

           cell_vecs = Reshape( config%cell, (/ 3, 3 /) )

           lengths( 1 ) = celprp( 1 )
           lengths( 2 ) = celprp( 2 )
           lengths( 3 ) = celprp( 3 )

           angles ( 1 ) = Acos( celprp( 5 ) )
           angles ( 2 ) = Acos( celprp( 6 ) )
           angles ( 3 ) = Acos( celprp( 4 ) )
           angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

! Print

           Call io_nc_put_var(io, 'config%cell'        , fh, cell_vecs, (/ 1, 1, jj /), (/ 3, 3, 1 /) )
           Call io_nc_put_var(io, 'config%cell_lengths', fh, lengths  , (/    1, jj /), (/    3, 1 /) )
           Call io_nc_put_var(io, 'config%cell_angles' , fh, angles   , (/    1, jj /), (/    3, 1 /) )

        End If

        Call io_close(io, fh )
        Call io_finalize(io)

     End If
     Call gsync(comm)

! Start of file

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec_mpi_io=Int(traj%rec_write,offset_kind)+Int(4,offset_kind)
     Else ! netCDF write
        Call gbcast(comm,jj,0)
        rec_mpi_io = Int(jj,offset_kind)
     End If

! Write the rest

     Call io_set_parameters(io, user_comm = comm%comm )
     Call io_init(io, traj%recsz_write )
     Call io_open(io, io_write, comm%comm, traj%fname, mode_wronly, fh )

     Call io_write_sorted_file(io, fh, traj%file_key(), IO_HISTORY, rec_mpi_io, config%natms, &
          config%ltg, config%atmnam, config%weight, rsd, config%parts,                   &
          config%vxx, config%vyy, config%vzz,  ierr )

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
        traj%rec_write=traj%rec_write+Int(4,li)+Int(config%megatm,li)*Int(traj%record_size,li)
        If (comm%idnode == 0) Then
           Write(record(1:traj%recsz_write), Fmt='(3i10,2i21,a1)') traj%file_key(),config%imcon,config%megatm,traj%frm_write,&
              traj%rec_write,lf
           Call io_write_record(io, fh, Int(1,offset_kind), record(1:traj%recsz_write) )
        End If
     End If

     Call io_close(io, fh )
     Call io_finalize(io)

! SORTED Serial Ditraj%rec_writet Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:config%mxatms),iwrk(1:config%mxatms),                   Stat=fail(1))
     Allocate (bxx(1:config%mxatms),byy(1:config%mxatms),bzz(1:config%mxatms),        Stat=fail(3))
     Allocate (temp_parts(1:config%mxatms),eee(1:config%mxatms),fff(1:config%mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O

     If (comm%idnode == 0) Then

        Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, Form='formatted', Access='direct', Recl=traj%recsz_write)

        traj%rec_write=traj%rec_write+Int(1,li)
        Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a8,2i10,2i2,2f20.6,a1)', Rec=traj%rec_write) &
          'timestep',nstep,config%megatm,traj%file_key(),config%imcon,tstep,time,lf

        Do i = 0, 2
           traj%rec_write=traj%rec_write+Int(1,li)
           Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3f20.10,a12,a1)', Rec=traj%rec_write) &
             config%cell(1+i*3),config%cell(2+i*3),config%cell(3+i*3),Repeat(' ',12),lf
        End Do

        Do i=1,config%natms
           iwrk(i)=config%ltg(i)

           temp_parts(i)%xxx=config%parts(i)%xxx
           temp_parts(i)%yyy=config%parts(i)%yyy
           temp_parts(i)%zzz=config%parts(i)%zzz

           chbuf(i)=config%atmnam(i)
           temp_parts(i)%chge=config%parts(i)%chge
           eee(i)=config%weight(i)
           fff(i)=rsd(i)
        End Do

        If (Any(traj%key == [TRAJ_KEY_COORD_VEL,TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
           Do i=1,config%natms
              bxx(i)=config%vxx(i)
              byy(i)=config%vyy(i)
              bzz(i)=config%vzz(i)
           End Do
        End If

        If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
           Do i=1,config%natms
              temp_parts(i)%fxx=config%parts(i)%fxx
              temp_parts(i)%fyy=config%parts(i)%fyy
              temp_parts(i)%fzz=config%parts(i)%fzz
           End Do
        End If

        jatms=config%natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,Traject_tag)

              Call grecv(comm,jatms,jdnode,Traject_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,temp_parts(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,eee(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,fff(1:jatms),jdnode,Traject_tag)

                 If (traj%key /= TRAJ_KEY_COORD) Then
                    Call grecv(comm,bxx(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,byy(1:jatms),jdnode,Traject_tag)
                    Call grecv(comm,bzz(1:jatms),jdnode,Traject_tag)
                 End If

                 !If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
                 !   Call grecv(comm,cxx(1:jatms),jdnode,Traject_tag)
                 !   Call grecv(comm,cyy(1:jatms),jdnode,Traject_tag)
                 !   Call grecv(comm,czz(1:jatms),jdnode,Traject_tag)
                 !End If
              End If
           End If

           Do i=1,jatms
              rec1=traj%rec_write+Int(iwrk(i)-1,li)*Int(traj%record_size,li)+Int(1,li)
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a8,i10,3f12.6,a18,a1)', Rec=rec1) &
                chbuf(i),iwrk(i),eee(i),temp_parts(i)%chge,fff(i),Repeat(' ',18),lf

              rec1=rec1+Int(1,li)
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3g20.10,a12,a1)', Rec=rec1) &
                temp_parts(i)%xxx,temp_parts(i)%yyy,temp_parts(i)%zzz,Repeat(' ',12),lf

              If (Any(traj%key == [TRAJ_KEY_COORD_VEL,TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
                 rec1=rec1+Int(1,li)
                 Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3g20.10,a12,a1)', Rec=rec1) &
                   bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
              End If

              If (Any(traj%key == [TRAJ_KEY_COORD_VEL_FORCE,TRAJ_KEY_COMPRESSED])) Then
                 rec1=rec1+Int(1,li)
                 Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3g20.10,a12,a1)', Rec=rec1) &
                   temp_parts(i)%xxx,temp_parts(i)%yyy,temp_parts(i)%zzz,Repeat(' ',12),lf
              End If
           End Do
        End Do

! Update main header

        traj%rec_write=traj%rec_write+Int(config%megatm,li)*Int(traj%record_size,li)
        Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3i10,2i21,a1)', Rec=Int(2,li)) &
          traj%file_key(),config%imcon,config%megatm,traj%frm_write,traj%rec_write,lf

        Close(Unit=files(FILE_HISTORY)%unit_no)

     Else

        Call grecv(comm,ready,0,Traject_tag)

        Call gsend(comm,config%natms,0,Traject_tag)
        If (config%natms > 0) Then
           Call gsend(comm,config%atmnam(:),0,Traject_tag)
           Call gsend(comm,config%ltg(:),0,Traject_tag)

           Call gsend(comm,config%parts(:),0,Traject_tag)
           Call gsend(comm,config%weight(:),0,Traject_tag)
           Call gsend(comm,rsd(:),0,Traject_tag)


           If (traj%key /= TRAJ_KEY_COORD) Then
              Call gsend(comm,config%vxx(:),0,Traject_tag)
              Call gsend(comm,config%vyy(:),0,Traject_tag)
              Call gsend(comm,config%vzz(:),0,Traject_tag)
           End If

        End If

! Save offset pointer

        traj%rec_write=traj%rec_write+Int(4,li)+Int(config%megatm,li)*Int(traj%record_size,li)

     End If

     Deallocate (chbuf,iwrk,         Stat=fail(1))
     Deallocate (bxx,byy,bzz,        Stat=fail(3))
     Deallocate (temp_parts,eee,fff, Stat=fail(5))
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

  If (traj%newjob_write) Then
     traj%newjob_write = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        traj%fname = files(FILE_HISTORY)%filename
     Else
        traj%fname = Trim(files(FILE_HISTORY)%filename) // '.nc'
     End If

! If keyres=1, is HISTORY old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File=traj%fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

110  Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (comm%idnode == 0) Then
              Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, &
                Form='formatted', Access='direct', Status='replace', Recl=traj%recsz_write)
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a34,a1)',      Rec=Int(1,li)) config%cfgname(1:34),lf
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(2i2,3i10,a1)', Rec=Int(2,li)) &
                traj%file_key(),config%imcon,config%megatm,traj%frm_write,traj%rec_write,lf
              Close(Unit=files(FILE_HISTORY)%unit_no)
           End If
           traj%rec_write=Int(2,li)
           traj%frm_write=Int(0,li)
        Else
           If (comm%idnode == 0) Then
              Call io_set_parameters(io, user_comm = comm_self )
              Call io_nc_create( netcdf, comm_self, traj%fname, config%cfgname, config%megatm )
           End If
        End If

! Get some sense of it

     Else

        safe=.true.
        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           If (comm%idnode == 0) Then

              Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, Form='formatted')

              Do

                 record(1:traj%recsz_write)=' '

! Assume new style of HISTORY with bookkeeping.

                 If (traj%fast_write) Then

                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt=*, End=120)                     ! title record
                    traj%rec_write=traj%rec_write+Int(1,li)
                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a)', End=120) record(1:traj%recsz_write) ! bookkeeping record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    traj%rec_write=traj%rec_write+Int(1,li)

                    Call get_word(record(1:traj%recsz_write),word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record(1:traj%recsz_write),word) ; Call get_word(record(1:traj%recsz_write),word)
                       Call get_word(record(1:traj%recsz_write),word) ; traj%frm_write=Nint(word_2_real(word,0.0_wp),li)
                       Call get_word(record(1:traj%recsz_write),word) ; traj%rec_write=Nint(word_2_real(word,0.0_wp),li)
                       If (traj%frm_write /= Int(0,li) .and. traj%rec_write > Int(2,li)) Then
                          Go To 120 ! New style
                       Else
                          traj%fast_write=.false. ! TOUGH, old style
                          traj%rec_write=Int(2,li)
                          traj%frm_write=Int(0,li)
                       End If
                    Else
                       safe=.false. ! Overwrite the file, it's junk to me
                       Go To 120
                    End If

! TOUGH, it needs scanning through

                 Else

                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt=*, End=120)                 ! timestep record
                    traj%rec_write=traj%rec_write+Int(1,li)

                    word=' '
                    i = 3 + config%megatm ! total number of lines to read
                    Write(word,'( "(", i0, "( / ) )" )') i-1
                    Read(Unit=files(FILE_HISTORY)%unit_no, Fmt=word, End=120)
                    traj%rec_write=traj%rec_write+Int(i,li)
                    traj%frm_write=traj%frm_write+Int(1,li)

                 End If

              End Do

120         Continue
            Close(Unit=files(FILE_HISTORY)%unit_no)

           End If

           Call gcheck(comm,safe,"enforce")
           If (.not.safe) Then
              lexist=.false.

              traj%rec_write=Int(0,li)
              traj%frm_write=Int(0,li)

              Go To 110
           Else If (comm%mxnode > 1) Then
              buffer(1)=Real(traj%frm_write,wp)
              buffer(2)=Real(traj%rec_write,wp)

              Call gbcast(comm,buffer,0)

              traj%frm_write=Nint(buffer(1),li)
              traj%rec_write=Nint(buffer(2),li)
           End If

        Else ! netCDF read

           If (comm%idnode == 0) Then
              Call io_set_parameters(io, user_comm = comm_self )
              Call io_open(io, io_write, comm_self, traj%fname, mode_rdonly, fh )

! Get the precision that the history file was written in
! and check it matches the requested precision

              Call io_nc_get_file_real_precision(io, fh, file_p, file_r, ierr )
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
              Call io_nc_get_real_precision( netcdf, io_p, io_r, ierr )
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
              Select Case( Selected_real_kind( io_p, io_r ) )
              Case( Kind( 1.0 ) )
                 Write(message,'(a)') 'Precision requested: Single'
              Case( Kind( 1.0d0 ) )
                 Write(message,'(a)') 'Precision requested: Double'
              End Select
              Call info(message,.true.)
              Select Case( Selected_real_kind( file_p, file_r ) )
              Case( Kind( 1.0 ) )
                 Write(message,'(a)') "Precision in file: Single"
              Case( Kind( 1.0d0 ) )
                 Write(message,'(a)') "Precision in file: Double"
              End Select
              Call info(message,.true.)

              Call gsync(comm)
              Call error(0,'Requested writing precision inconsistent with that in an existing HISTORY.nc')
           End If

! Get the frame number to check
! For netCDF this is the "frame number" which is not a long integer!

           jj=0
           If (comm%idnode == 0) Call io_nc_get_dim(io, 'frame', fh, jj )
           Call gbcast(comm,jj,0)

           If (jj > 0) Then
              traj%frm_write=Int(jj,li)

           If (comm%idnode == 0) Call io_close(io, fh )
           Else ! Overwrite the file, it's junk to me
              lexist=.false.

              traj%rec_write=Int(0,li)
              traj%frm_write=Int(0,li)

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
     Allocate (chbat(1:traj%recsz_write,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure 0'
        Call error(0,message)
     End If

     chbat=' '
     n_atm=0 ; n_atm(comm%idnode+1)=config%natms
     Call gsum(comm,n_atm)
     n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Update frame

  traj%frm_write=traj%frm_write+Int(1,li)

! NO UNSORTED WRITING AS NO INDICES ARE WRITTEN
! In other words even if the user asked for unsorted I/O they get it sorted.

! SORTED MPI-I/O or Parallel Ditraj%rec_writet Access FORTRAN

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT   .or. &
      io_write == IO_WRITE_SORTED_NETCDF) Then

     Call io_set_parameters(io, user_comm = comm%comm )
     Call io_init(io, traj%recsz_write )
     Call io_open(io, io_write, comm%comm, traj%fname, mode_wronly, fh )

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(traj%rec_write,offset_kind)
     jj_io=rec_mpi_io
     jj=0 ! netCDF current frame
     If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm = comm_self )
        Call io_init(io, traj%recsz_write )
        Call io_open(io, io_write, comm_self, traj%fname, mode_wronly, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and cell information

           Write(record(1:traj%recsz_write), Fmt='(a8,i6,f8.5,f12.5,a1)') 'timestep',nstep,tstep,time,lf
           Call io_write_record(io, fh, jj_io, record(1:traj%recsz_write) )
           jj_io=jj_io + Int(1,offset_kind)

           Do i = 0, 2
              Write(record(1:traj%recsz_write), Fmt='(3f10.3,a4,a1)') &
                   config%cell( 1 + i * 3 ), config%cell( 2 + i * 3 ), config%cell( 3 + i * 3 ), Repeat( ' ', 4 ), lf
              Call io_write_record(io, fh, jj_io, record(1:traj%recsz_write) )
              jj_io=jj_io+Int(1,offset_kind)
           End Do

        Else ! netCDF write

! Get the current and new frame numbers

           Call io_nc_get_dim(io, 'frame', fh, jj )
           jj = jj + 1

           Call io_nc_put_var(io, 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var(io, 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var(io, 'timestep '      , fh,  tstep, jj, 1 )

           Call dcell(config%cell,celprp) ! get config%cell properties

           cell_vecs = Reshape( config%cell, (/ 3, 3 /) )

           lengths( 1 ) = celprp( 1 )
           lengths( 2 ) = celprp( 2 )
           lengths( 3 ) = celprp( 3 )

           angles ( 1 ) = Acos( celprp( 5 ) )
           angles ( 2 ) = Acos( celprp( 6 ) )
           angles ( 3 ) = Acos( celprp( 4 ) )
           angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

! Print

           Call io_nc_put_var(io, 'config%cell'        , fh, cell_vecs, (/ 1, 1, jj /), (/ 3, 3, 1 /) )
           Call io_nc_put_var(io, 'config%cell_lengths', fh, lengths  , (/    1, jj /), (/    3, 1 /) )
           Call io_nc_put_var(io, 'config%cell_angles' , fh, angles   , (/    1, jj /), (/    3, 1 /) )

        End If

        Call io_close(io, fh )
        Call io_finalize(io)

     End If
     Call gsync(comm)

! Start of file

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec_mpi_io=Int(traj%rec_write,offset_kind)+Int(4,offset_kind)
     Else ! netCDF write
        Call gbcast(comm,jj,0)
        rec_mpi_io = Int(jj,offset_kind)
     End If

! Write the rest

     Call io_set_parameters(io, user_comm = comm%comm )
     Call io_init(io, traj%recsz_write )
     Call io_open(io, io_write, comm%comm, traj%fname, mode_wronly, fh )

     Call io_write_sorted_file(io, fh, 0*traj%file_key(), IO_HISTORD, rec_mpi_io, config%natms, &
          config%ltg, config%atmnam, (/ 0.0_wp /),  rsd, config%parts,     &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                        &
          IO_SUBSET_POSITIONS,  ierr )

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
        traj%rec_write=traj%rec_write+Int(4,li)+Int(config%megatm,li)
        If (comm%idnode == 0) Then
           Write(record(1:traj%recsz_write), Fmt='(2i2,3i10,a1)') traj%file_key(),config%imcon,config%megatm,traj%frm_write,&
             traj%rec_write,lf
           Call io_write_record(io, fh, Int(1,offset_kind), record(1:traj%recsz_write) )
        End If
     End If

     Call io_close(io, fh )
     Call io_finalize(io)

! SORTED Serial Ditraj%rec_writet Access FORTRAN

  Else

     Allocate (chbuf(1:config%mxatms),iwrk(1:config%mxatms),            Stat=fail(1))
     Allocate (temp_parts(1:config%mxatms),                      Stat=fail(2))
     Allocate (fff(1:config%mxatms),                             Stat=fail(3))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'trajectory_write allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O

     If (comm%idnode == 0) Then

        Open(Newunit=files(FILE_HISTORY)%unit_no, File=traj%fname, Form='formatted', Access='direct', Recl=traj%recsz_write)

        traj%rec_write=traj%rec_write+Int(1,li)
        Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a8,i6,f8.5,f12.5,a1)', Rec=traj%rec_write) 'timestep',nstep,tstep,time,lf

        Do i = 0, 2
           traj%rec_write=traj%rec_write+Int(1,li)
           Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(3f10.3,a4,a1)', Rec=traj%rec_write) &
             config%cell(1+i*3),config%cell(2+i*3),config%cell(3+i*3),Repeat(' ',4),lf
        End Do

        Do i=1,config%natms
           iwrk(i)=config%ltg(i)

           temp_parts(i)%xxx=config%parts(i)%xxx
           temp_parts(i)%yyy=config%parts(i)%yyy
           temp_parts(i)%zzz=config%parts(i)%zzz

           chbuf(i)=config%atmnam(i)
           fff(i)=rsd(i)
        End Do

        jatms=config%natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,Traject_tag)

              Call grecv(comm,jatms,jdnode,Traject_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,Traject_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,fff(1:jatms),jdnode,Traject_tag)

                 Call grecv(comm,temp_parts(1:jatms),jdnode,Traject_tag)
              End If
           End If

           Do i=1,jatms
              rec1=traj%rec_write+Int(iwrk(i),li)
              Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(a6,4f7.1,a1)', Rec=rec1) &
                chbuf(i),temp_parts(i)%xxx,temp_parts(i)%yyy,temp_parts(i)%zzz,fff(i),lf
           End Do
        End Do

! Update main header

        traj%rec_write=traj%rec_write+Int(config%megatm,li)
        Write(Unit=files(FILE_HISTORY)%unit_no, Fmt='(2i2,3i10,a1)', Rec=Int(2,li)) &
          traj%file_key(),config%imcon,config%megatm,traj%frm_write,traj%rec_write,lf

        Close(Unit=files(FILE_HISTORY)%unit_no)

     Else

        Call grecv(comm,ready,0,Traject_tag)

        Call gsend(comm,config%natms,0,Traject_tag)
        If (config%natms > 0) Then
           Call gsend(comm,config%atmnam(:),0,Traject_tag)
           Call gsend(comm,config%ltg(:),0,Traject_tag)

           Call gsend(comm,rsd(:),0,Traject_tag)

           Call gsend(comm,config%parts(:),0,Traject_tag)
        End If

! Save offset pointer

        traj%rec_write=traj%rec_write+Int(4,li)+Int(config%megatm,li)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (temp_parts, Stat=fail(2))
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
