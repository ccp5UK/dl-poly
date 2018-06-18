Module defects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global defects variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,             Only : wp,li, wi
  Use setup,             Only : mxatms,mxbfxp,ndefdt, &
                                nrefdt,config,half_minus, zero_plus
  Use comms,             Only : comms_type, DefWrite_tag, wp_mpi, DefExport_tag, &
                                DefRWrite_tag,gsum,gcheck,gsync,gmax,gbcast, &
                                gsend,grecv,gwait,girecv,gscatter,gscatterv, &
                                gscatter_columns,offset_kind,mode_wronly, &
                                comm_self,mode_create,mode_rdonly
  Use configuration,     Only : cfgname,imcon,cell,natms,nlast, &
                                atmnam,ltg,lfrzn,xxx,yyy,zzz
  Use parse,             Only : tabs_2_blanks,get_word,word_2_real,get_line,strip_blanks
  Use io,                Only : io_set_parameters,        &
                                io_get_parameters,        &
                                io_init, io_open,         &
                                io_write_record,          &
                                io_write_batch,           &
                                io_close, io_finalize,    &
                                io_get_var,               &
                                io_nc_get_var,            &
                                io_read_batch,            &
                                io_nc_create,             &
                                io_nc_put_var,            &
                                io_write_sorted_file,     &
                                io_delete,                &
                                IO_WRITE_UNSORTED_MPIIO,  &
                                IO_WRITE_UNSORTED_DIRECT, &
                                IO_WRITE_UNSORTED_MASTER, &
                                IO_WRITE_SORTED_MPIIO,    &
                                IO_WRITE_SORTED_DIRECT,   &
                                IO_WRITE_SORTED_NETCDF,   &
                                IO_WRITE_SORTED_MASTER,   &
                                IO_READ_MASTER,           &
                                IO_READ_NETCDF,           &
                                IO_RESTART,               &
                                IO_BASE_COMM_NOT_SET,     &
                                IO_ALLOCATION_ERROR,      &
                                IO_UNKNOWN_WRITE_OPTION,  &
                                IO_UNKNOWN_WRITE_LEVEL
  Use site, Only : site_type
  Use domains,           Only : nprx,npry,nprz,            &
                                nprx_r,npry_r,nprz_r,map,  &
                                idx,idy,idz,r_nprx,r_npry, &
                                r_nprz
  Use numerics,          Only : pbcshift,invert,dcell, shellsort2
  Use errors_warnings,   Only : error,warning
  Use neighbours,        Only : neighbours_type,defects_link_cells
  Use core_shell,        Only : core_shell_type

  Implicit None

  Private
!> REFERENCE data

!> Defect Type  
   Type, Public :: defects_type
     Logical                             :: ldef, newjob
     Integer(Kind =  wi)                 :: isdef, nsdef
     Real(Kind =  wp)                    :: rdef   
     Character( Len = 12 )               :: reffile
     Character( Len = 12 )               :: deffile

     Integer(Kind =  wi)                 :: nrefs = 0 , nlrefs = 0 
     Integer(Kind=li)                    :: rec=0_li,frm = 0_li
     Integer(Kind =  wi)                 :: mxlcdef

     Real( Kind = wp )                   :: celr(1:9)  = 0.0_wp
     Real( Kind = wp )                   :: rdefsq,rcell(1:9),cwx,cwy,cwz, &
                                            dxl,dxr,dyl,dyr,dzl,dzr,cutdef

     Integer(Kind =  wi)   , Allocatable :: lri(:) ,lra(:) ,indr(:)
     Real( Kind = wp )     , Allocatable :: xr(:)  ,yr(:)  ,zr(:)
     Character( Len = 8 )  , Allocatable :: namr(:)
   End Type

  Public :: defects_write

Contains


!> allocate_defects_arrays 
  Subroutine allocate_defects_arrays(dfcts)

    Integer(Kind =  wi), Dimension( 1:3 )  :: fail
    Type(defects_type) , Intent( InOut )   :: dfcts 

    fail = 0
    Allocate (dfcts%namr(1:mxatms)                                 , Stat = fail(1))
    Allocate (dfcts%lri(1:mxatms) ,dfcts%lra(1:mxatms),dfcts%indr(1:mxatms), Stat = fail(2))
    Allocate (dfcts%xr(1:mxatms)  ,dfcts%yr(1:mxatms) ,dfcts%zr(1:mxatms)  , Stat = fail(3))

    If (Any(fail > 0)) Call error(1035)

    dfcts%namr = ' '
    dfcts%lri  = 0 
    dfcts%lra  = 0 
    dfcts%indr = 0

    dfcts%xr   = 0.0_wp 
    dfcts%yr   = 0.0_wp 
    dfcts%zr   = 0.0_wp
 
  End Subroutine allocate_defects_arrays


!> deallocate_defects_arrays 
  Subroutine deallocate_defects_arrays(dfcts)

    Integer(Kind =  wi), Dimension( 1:3 )  :: fail
    Type(defects_type) , Intent( InOut )   :: dfcts 

    fail = 0
    Deallocate (dfcts%namr                     , Stat = fail(1))
    Deallocate (dfcts%lri ,dfcts%lra,dfcts%indr, Stat = fail(2))
    Deallocate (dfcts%xr  ,dfcts%yr ,dfcts%zr  , Stat = fail(3))

    If (Any(fail > 0)) Call error(1035)
 
  End Subroutine deallocate_defects_arrays

!> defects_reference_export   
  Subroutine defects_reference_export(mdir,ixyz,dfcts,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic REFERENCE data in domain boundary
! regions for halo formation
!
! all particle coordinates are in reduced space with origin localised
! onto this node (idnode)
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,                Intent( In    ) :: mdir
  Integer,                Intent( InOut ) :: ixyz(1:mxatms)
  Type( defects_type ),   Intent( InOut ) :: dfcts
  Type( comms_type )  ,   Intent( InOut ) :: comm

  Logical           :: safe,lsx,lsy,lsz,lex,ley,lez
  Integer           :: fail,iadd,limit,iblock,          &
                       i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message

! Number of transported quantities per particle

  iadd=13

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'defects_reference_export allocation failure'
     Call error(0,message)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/Merge(2,1,comm%mxnode > 1)


! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
! kxyz - corrected halo reduction factor particles haloing both +&- sides
! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
! jdnode - destination (send to), kdnode - source (receive from)

  kx = 0 ; ky = 0 ; kz = 0
  lsx = .false. ; lex = .false.
  lsy = .false. ; ley = .false.
  lsz = .false. ; lez = .false.
  If      (mdir == -1) Then ! Direction -x
     kx  = 1
     jxyz= 1
     kxyz= 3
     lsx = (idx == 0)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= 3
     lex = (idx == nprx-1)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= 30
     lsy = (idy == 0)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= 30
     ley = (idy == npry-1)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= 300
     lsz = (idz == 0)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= 300
     lez = (idz == nprz-1)

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(557)
  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,dfcts%nlrefs

! If the particle is within the remaining 'inverted halo' of this domain

     If (ixyz(i) > 0) Then

! Get the necessary halo indices

        ix=Mod(ixyz(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

! If safe to proceed

           If ((imove+iadd) <= iblock) Then

! pack positions and apply possible wrap-around corrections for receiver

              buffer(imove+ 1)=dfcts%xr(i)
              If (lsx) buffer(imove+1)=buffer(imove+1)+1.0_wp
              If (lex) buffer(imove+1)=buffer(imove+1)-1.0_wp
              buffer(imove+ 2)=dfcts%yr(i)
              If (lsy) buffer(imove+2)=buffer(imove+2)+1.0_wp
              If (ley) buffer(imove+2)=buffer(imove+2)-1.0_wp
              buffer(imove+ 3)=dfcts%zr(i)
              If (lsz) buffer(imove+3)=buffer(imove+3)+1.0_wp
              If (lez) buffer(imove+3)=buffer(imove+3)-1.0_wp

! pack config indexing, site name and remaining halo indexing arrays

              buffer(imove+ 4)=Real(dfcts%indr(i),wp)
              buffer(imove+ 5)=Real(Ichar(dfcts%namr(i)(1:1)),wp)
              buffer(imove+ 6)=Real(Ichar(dfcts%namr(i)(2:2)),wp)
              buffer(imove+ 7)=Real(Ichar(dfcts%namr(i)(3:3)),wp)
              buffer(imove+ 8)=Real(Ichar(dfcts%namr(i)(4:4)),wp)
              buffer(imove+ 9)=Real(Ichar(dfcts%namr(i)(5:5)),wp)
              buffer(imove+10)=Real(Ichar(dfcts%namr(i)(6:6)),wp)
              buffer(imove+11)=Real(Ichar(dfcts%namr(i)(7:7)),wp)
              buffer(imove+12)=Real(Ichar(dfcts%namr(i)(8:8)),wp)

! Use the corrected halo reduction factor when the particle is halo to both +&- sides

              buffer(imove+13)=Real(ixyz(i)-Merge(jxyz,kxyz,j == jxyz),wp)

           Else

              safe=.false.

           End If
           imove=imove+iadd

        End If

     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  Call gcheck(comm,safe)
  If (.not.safe) Then
     itmp=Merge(2,1,comm%mxnode > 1)*imove
     Call gmax(comm,itmp)
     Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
     Call error(558)
  End If

! exchange information on buffer sizes

  If (comm%mxnode > 1) Then
     Call girecv(comm,jmove,kdnode,DefExport_tag)
     Call gsend(comm,imove,jdnode,DefExport_tag)
     Call gwait(comm)
  Else
     jmove=imove
  End If

! Check for array bound overflow (can arrays cope with incoming data)

  safe=((dfcts%nlrefs+jmove/iadd) <= mxatms)
  Call gcheck(comm,safe)
  If (.not.safe) Then
     itmp=dfcts%nlrefs+jmove/iadd
     Call gmax(comm,itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
     Call error(559)
  End If

! exchange buffers between nodes (this is a MUST)

  If (comm%mxnode > 1) Then
    If (jmove > 0) Then
      Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,DefExport_tag)
    End If
    If (imove > 0) Then
      Call gsend(comm,buffer(1:imove),jdnode,DefExport_tag)
    End If
    If (jmove > 0) Then
      Call gwait(comm)
    End If
  End If

! load transferred data

  j=Merge(iblock,0,comm%mxnode > 1)
  Do i=1,jmove/iadd
     dfcts%nlrefs=dfcts%nlrefs+1

! unpack positions

     dfcts%xr(dfcts%nlrefs)=buffer(j+1)
     dfcts%yr(dfcts%nlrefs)=buffer(j+2)
     dfcts%zr(dfcts%nlrefs)=buffer(j+3)

! unpack config indexing, site name halo indexing arrays

     dfcts%indr(dfcts%nlrefs)=Nint(buffer(j+4))
     dfcts%namr(dfcts%nlrefs)(1:1)=Char(Nint(buffer(j+ 5)))
     dfcts%namr(dfcts%nlrefs)(2:2)=Char(Nint(buffer(j+ 6)))
     dfcts%namr(dfcts%nlrefs)(3:3)=Char(Nint(buffer(j+ 7)))
     dfcts%namr(dfcts%nlrefs)(4:4)=Char(Nint(buffer(j+ 8)))
     dfcts%namr(dfcts%nlrefs)(5:5)=Char(Nint(buffer(j+ 9)))
     dfcts%namr(dfcts%nlrefs)(6:6)=Char(Nint(buffer(j+10)))
     dfcts%namr(dfcts%nlrefs)(7:7)=Char(Nint(buffer(j+11)))
     dfcts%namr(dfcts%nlrefs)(8:8)=Char(Nint(buffer(j+12)))
     ixyz(dfcts%nlrefs)=Nint(buffer(j+13))

     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'defects_reference_export deallocation failure'
     Call error(0,message)
  End If

End Subroutine defects_reference_export

!> defects_reference_read

  Subroutine defects_reference_read(nstep,dfcts,site,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading particles data from REFERENCE file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  Integer,              Intent( In    ) :: nstep
  Type( defects_type ), Intent( InOut ) :: dfcts
  Type( site_type ), Intent( In    ) :: site
  Type( comms_type ),   Intent( InOut ) :: comm

  Logical                :: l_ind = .true.  , &
                            l_str = .false. , &
                            lexist,fast,safe,loop,match
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,fname
  Integer                :: fail(1:2),megref,i,j,lvcfgr,imconr, &
                            indatm,idm,ipx,ipy,ipz,             &
                            itmols,isite,nsite,msite,fsite,ifrz,nrept
  Real( Kind = wp )      :: det,sxx,syy,szz,cell_vecs(1:3,1:3)

! Some parameters and variables needed by io interfaces

  Integer                           :: recsz = 73 ! default record size
  Integer                           :: fh, io_read
  Integer( Kind = offset_kind ) :: top_skip

  Character( Len = 8 ), Dimension( : ), Allocatable :: chbuf
  Integer,              Dimension( : ), Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ), Allocatable :: axx,ayy,azz
  Character( Len = 256 ) :: message


5 Continue

! Get type of I/O for reading

  Call io_get_parameters( user_method_read = io_read )

! Define filename ASCII or netCDF

  If (io_read /= IO_READ_NETCDF) Then
     fname=trim(dfcts%reffile)
  Else
     fname=trim(dfcts%reffile) // '.nc'
  End If
 
! Check if we have a REFERENCE and default megref

  lexist=.true. ; megref=0
  If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
  Call gcheck(comm,lexist,"enforce")
  If (.not.lexist) Then
     If (nstep > 1) Then ! That is a problem.  Abort, we must have restarted
        Call error(551)
     Else                ! Use data from CONFIG
        Call warning(320,0.0_wp,0.0_wp,0.0_wp)
        If (io_read /= IO_READ_NETCDF) Then
           fname=Trim(config)
        Else
           fname=Trim(config)//'.nc'
        End If
        megref=natms
        Call gsum(comm,megref)
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
     If (comm%idnode == 0) Then

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
     Call gsync(comm)
     Call gcheck(comm,safe,"enforce")
     Call gcheck(comm,fast,"enforce")
     If (.not.safe) Go To 100

! Close REFERENCE

     If (comm%idnode == 0) Close(Unit=nrefdt)

     If      (fast       ) Then
        If (megref == 0) Then          ! Define megref
           If (comm%idnode == 0) megref = i ! It's already been read on master
           Call gsum(comm,megref)
        End If
     Else If (megref == 0) Then
        io_read = IO_READ_MASTER       ! Abort parallel reading
     End If

  End If

!!! SCAN HEADER

  If (io_read /= IO_READ_NETCDF) Then ! ASCII read

! Open file

     If (comm%idnode == 0) Open(Unit=nrefdt, File=fname)

! Read TITLE record (file header)

     Call get_line(safe,nrefdt,record,comm)
     If (.not.safe) Go To 100

! Read configuration level and image condition

     Call get_line(safe,nrefdt,record,comm)
     If (.not.safe) Go To 100

     Call get_word(record,word)
     lvcfgr=Nint(word_2_real(word))

     Call get_word(record,word)
     imconr=Nint(word_2_real(word))

     If (imconr /= 0) Then
        Call get_line(safe,nrefdt,record,comm)
        If (.not.safe) Go To 100
        Call get_word(record,word)
        dfcts%celr(1)=word_2_real(word)
        Call get_word(record,word)
        dfcts%celr(2)=word_2_real(word)
        Call get_word(record,word)
        dfcts%celr(3)=word_2_real(word)

        Call get_line(safe,nrefdt,record,comm)
        If (.not.safe) Go To 100
        Call get_word(record,word)
        dfcts%celr(4)=word_2_real(word)
        Call get_word(record,word)
        dfcts%celr(5)=word_2_real(word)
        Call get_word(record,word)
        dfcts%celr(6)=word_2_real(word)

        Call get_line(safe,nrefdt,record,comm)
        If (.not.safe) Go To 100
        Call get_word(record,word)
        dfcts%celr(7)=word_2_real(word)
        Call get_word(record,word)
        dfcts%celr(8)=word_2_real(word)
        Call get_word(record,word)
        dfcts%celr(9)=word_2_real(word)
     Else
        Call error(552) ! Lattice parameters are a must
     End If

! image conditions not compliant with DD and link-cell

     If (imconr == 4 .or. imconr == 5 .or. imconr == 7) Call error(300)

! Close REFERENCE

     If (comm%idnode == 0) Close(Unit=nrefdt)
     Call gsync(comm)

  Else ! netCDF read

! Open file

     Call io_set_parameters( user_comm = comm%comm )
     Call io_open( io_read, comm%comm, fname, mode_rdonly, fh )

     i=1 ! For config there is only one frame

     Call io_nc_get_var( 'datalevel'      , fh, lvcfgr, i, 1  )

     Call io_nc_get_var( 'imageconvention', fh, imconr, i, 1  )

! image conditions not compliant with DD and link-cell

     If (imconr == 4 .or. imconr == 5 .or. imconr == 7) Call error(300)

     Call io_nc_get_var( 'cell'           , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
     dfcts%celr = Reshape( cell_vecs, (/ Size( dfcts%celr ) /) )

! Close REFERENCE

     Call io_close( fh )

  End If

! REFERENCE to CONFIG cell match

  match=.true.
  If (.not.lexist) Then
     Do i=1,9
        match = match .and. Abs(cell(i)-dfcts%celr(i)) < 1.0e-6_wp
     End Do
  End If

! Get lattice invert

  Call invert(dfcts%celr,dfcts%rcell,det)

! If MASTER read

  If (io_read == IO_READ_MASTER) Then
     fail=0
     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'defects_reference_read allocation failure'
        Call error(0,message)
     End If

! Open file and skip header

     If (comm%idnode == 0) Then
        Open(Unit=nrefdt, File=fname)

        Read(Unit=nrefdt, Fmt=*) ! REFERENCE file header (TITLE record)
        Read(Unit=nrefdt, Fmt=*) ! configuration level and image condition

        Read(Unit=nrefdt, Fmt=*) ! cell vectors
        Read(Unit=nrefdt, Fmt=*)
        Read(Unit=nrefdt, Fmt=*)
     End If

! Initialise domain localised referent atoms counter

     dfcts%nrefs=0

! Initialise dispatched atom counter

     indatm=0

! Initialise total number of atoms and index counter

     megref=0
     j=0

     safe=.true.
     loop=.true.
     Do While (loop)

! Read in transmission arrays

        If (comm%idnode == 0 .and. safe) Then
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

        Call gcheck(comm,safe,"enforce")
        Call gcheck(comm,loop,"enforce")
        
        If (.not.safe) Go To 100 ! catch error
        If (loop) Then

! increase counters

           indatm=indatm+1
           megref=megref+1

        Else

! Close file

           If (comm%idnode == 0) Close(Unit=nrefdt)

! Check for inconsistencies in REFERENCE

           Call gmax(comm,j)
           If (j /= megref) Call error(553)

        End If

! Circulate configuration data to all nodes when transmission arrays
! are filled up or this is the last looping

        If (indatm == mxatms .or. ((.not.loop) .and. indatm > 0)) Then

! Ensure all atoms are in prescribed simulation cell (DD bound) and broadcast them
!
!           Call pbcshift(imconr,dfcts%celr,indatm,axx,ayy,azz)

           Call gbcast(comm,chbuf,0)
           Call gbcast(comm,iwrk,0)

           Call gbcast(comm,axx,0)
           Call gbcast(comm,ayy,0)
           Call gbcast(comm,azz,0)

! Assign atoms positions in fractional coordinates to the correct domains
! (DD bounding)

           Do i=1,indatm
              sxx=dfcts%rcell(1)*axx(i)+dfcts%rcell(4)*ayy(i)+dfcts%rcell(7)*azz(i)
              syy=dfcts%rcell(2)*axx(i)+dfcts%rcell(5)*ayy(i)+dfcts%rcell(8)*azz(i)
              szz=dfcts%rcell(3)*axx(i)+dfcts%rcell(6)*ayy(i)+dfcts%rcell(9)*azz(i)

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
              If      (idm < 0 .or. idm > (comm%mxnode-1)) Then
                 Call error(555)
              Else If (idm == comm%idnode)                 Then
                 dfcts%nrefs=dfcts%nrefs+1

                 If (dfcts%nrefs < mxatms) Then
                    dfcts%namr(dfcts%nrefs)=chbuf(i)
                    dfcts%indr(dfcts%nrefs)=iwrk(i)

                    dfcts%xr(dfcts%nrefs)=sxx
                    dfcts%yr(dfcts%nrefs)=syy
                    dfcts%zr(dfcts%nrefs)=szz
                 Else
                    safe=.false.
                 End If
              End If
           End Do

! Check if all is dispatched fine

           Call gcheck(comm,safe,"enforce")
           If (.not.safe) Call error(556)

! Nullify dispatch counter

           indatm=0

        End If

     End Do

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'defects_reference_read deallocation failure'
        Call error(0,message)
     End If

! If PROPER read

  Else

! Open file

     If (fast) Then
        Call io_set_parameters( user_comm = comm%comm )
        Call io_init( recsz )
        Call io_open( io_read, comm%comm, fname, mode_rdonly, fh )
     Else
        Open(Unit=nrefdt, File=fname)
     End If

! top_skip is header size

     If (io_read /= IO_READ_NETCDF) Then
        top_skip = Int(5,offset_kind) ! imcon is a must
     Else
        top_skip = Int(1,offset_kind) ! This is now the frame = 1
     End If

     Call defects_reference_read_parallel       &
           (lvcfgr,l_ind,l_str,megref,fast,fh,top_skip,dfcts,comm)

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

  If (Trim(fname) /= Trim(config)) Then

     nsite=0
     msite=0
     fsite=0
     Do itmols=1,site%ntype_mol
        Do isite=1,site%num_site(itmols)
           nsite=nsite+1

           If (site%freeze_site(nsite) /= 0) Then
              Do nrept=1,site%num_mols(itmols)
                 ifrz=nsite+msite+(nrept-1)*site%num_site(itmols)

                 Do i=1,dfcts%nrefs
                    If (dfcts%indr(i) == ifrz) Then
                       dfcts%namr(i)=dfcts%namr(dfcts%nrefs) ; dfcts%namr(dfcts%nrefs)=' '
                       dfcts%indr(i)=dfcts%indr(dfcts%nrefs) ; dfcts%indr(dfcts%nrefs)=0

                       dfcts%xr(i)=dfcts%xr(dfcts%nrefs) ; dfcts%xr(dfcts%nrefs)=0.0_wp
                       dfcts%yr(i)=dfcts%yr(dfcts%nrefs) ; dfcts%yr(dfcts%nrefs)=0.0_wp
                       dfcts%zr(i)=dfcts%zr(dfcts%nrefs) ; dfcts%zr(dfcts%nrefs)=0.0_wp

                       dfcts%nrefs=dfcts%nrefs-1
                    End If
                 End Do
              End Do
           End If
        End Do

        msite=msite+(site%num_mols(itmols)-1)*site%num_site(itmols)
        fsite=fsite+site%num_mols(itmols)*site%num_freeze(itmols)
     End Do

     If (fsite > 0) Then
        nsite=dfcts%nrefs
        Call gsum(comm,nsite)
        megref=nsite
     End If
  End If

! MATCH glitch fix

  If (.not.match) Then
     Call defects_reference_write(fname,megref,dfcts,comm)
     Go To 5
  End If

  Return

! REFERENCE format failure

100 Continue
  If (comm%idnode == 0) Close(Unit=nrefdt)
  Call error(554)

End Subroutine defects_reference_read

Subroutine defects_reference_read_parallel      &
           (lvcfgr,l_ind,l_str,megref,fast,fh,top_skip,dfcts,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the REFERENCE data file
! in parallel
!
! copyright - daresbury laboratory
! author    - i.j.bush & i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Logical,                           Intent( In    ) :: l_ind,l_str,fast
  Integer,                           Intent( In    ) :: lvcfgr,megref,fh
  Integer( Kind = offset_kind ),     Intent( In    ) :: top_skip
  Type( defects_type ),              Intent( InOut ) :: dfcts
  Type( comms_type )  ,              Intent( InOut ) :: comm

  Logical                :: safe,do_read
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,forma
  Integer                :: fail(1:6),i,j,k,                &
                            idm,ipx,ipy,ipz,indatm,         &
                            n_read_procs_use,per_read_proc, &
                            my_read_proc_num,ats_per_proc,  &
                            recs_per_at,recs_per_proc,      &
                            wp_vals_per_at,n_loc,           &
                            to_read,which_read_proc,this_base_proc
  Integer( Kind = li )   :: n_sk,n_ii,n_jj
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz

! Some parameters and variables needed by io interfaces

  Integer                           :: io_read,ierr
  Integer                           :: recsz, batsz
  Integer( Kind = offset_kind )     :: rec_mpi_io, n_skip
  Integer                           :: this_rec_buff, recs_to_read

! netCDF

  Integer :: frame, start(1:3), count(1:3)

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx_read,ayy_read,azz_read

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf_read,chbuf_scat
  Integer,              Dimension( : ),    Allocatable :: iwrk_read,iwrk_scat

  Integer,              Dimension( : ),    Allocatable :: first_at,orig_first_at
  Integer,              Dimension( : ),    Allocatable :: n_held,where_buff
  Integer,              Dimension( : ),    Allocatable :: owner_read

  Real( Kind = wp ),    Dimension( :, : ), Allocatable :: scatter_buffer_read
  Real( Kind = wp ),    Dimension( :, : ), Allocatable :: scatter_buffer

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: rec_buff
  Character( Len = 256 ) :: message


! Get reading method, total number of I/O heads and buffer size

  Call io_get_parameters( user_method_read      = io_read          )
  Call io_get_parameters( user_n_io_procs_read  = n_read_procs_use )
  Call io_get_parameters( user_buffer_size_read = batsz            )

  fail = 0 ! fail initialisation

  wp_vals_per_at = 3              ! Scatter buffer sizes
  recs_per_at    = 2 + lvcfgr     ! Scatter buffer sizes

! Note: make 'first_at' and 'orig_first_at' 1 element bigger than strictly
! required to make checking at the end of reading much easier and clearer

  Allocate (first_at(0:n_read_procs_use),orig_first_at(0:n_read_procs_use), Stat=fail(1))
  Allocate (chbuf(1:batsz),iwrk(1:batsz),                                   Stat=fail(2))
  Allocate (scatter_buffer(1:wp_vals_per_at,1:batsz),                       Stat=fail(3))
  If (Any(fail(1:3) > 0)) Then
     Write(message,'(a)') 'defects_reference_read_parallel allocation failure 1'
     Call error(0,message)
  End If

! define basic quantities for the parallel ASCII reading

  per_read_proc = comm%mxnode / n_read_procs_use
  do_read = (Mod( comm%idnode, per_read_proc ) == 0 .and. comm%idnode < per_read_proc * n_read_procs_use)
  my_read_proc_num = comm%idnode / per_read_proc

! Note 'first_at' and 'orig_first_at' have one more element
! in the array than strictly required - makes it easier to
! check that reading by the last I/O processor has finished

  ats_per_proc = megref / n_read_procs_use
  Do i=0,n_read_procs_use
     first_at(i) = i*ats_per_proc + Min(i,megref-ats_per_proc*n_read_procs_use)
  End Do
  orig_first_at = first_at
  ats_per_proc = Max(1,ats_per_proc) ! Fix it if 0
  recs_per_proc = ats_per_proc * recs_per_at

! Catch the case where the first atom belonging to
! a read processor does not actually exists - i.e.
! I/O procs count > megref, and limit reading by do_read

  If (my_read_proc_num < n_read_procs_use) &
     do_read = (do_read .and. first_at(my_read_proc_num) < megref)

! Skip to the point of reading

  If (do_read) Then

     n_skip = Int(recs_per_at,offset_kind) * Int(first_at(my_read_proc_num),offset_kind) + &
              top_skip-Int(1,offset_kind)

     If (.not.fast) Then
        n_sk=Int(n_skip,li)
        n_jj=73*batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
           Do n_ii=1_li,n_sk/n_jj
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_jj
              Read(Unit=nrefdt, Fmt=forma, End=100)
           End Do
           n_ii=Mod(Int(n_skip,li),n_jj)
           If (n_ii > 0_li) Then
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_ii
              Read(Unit=nrefdt, Fmt=forma, End=100)
           End If
        Else
           forma=' '
           Write(forma,'( "(", i0, "/)" )') n_sk
           Read(Unit=nrefdt, Fmt=forma, End=100)
        End If

        recsz=200
        forma=' '
        Write(forma,'( "(", i0, "a1)" )') recsz
     Else
        rec_mpi_io = n_skip + Int(1,offset_kind)
        recsz=73
     End If

! Allocate record buffer, reading buffers, scatter buffers and indexing arrays

     If (io_read /= IO_READ_NETCDF) Then
        Allocate (rec_buff(1:recsz,1:batsz),                                  Stat=fail(1))
     Else
        Allocate (rec_buff(1:Len( chbuf_read ),1:batsz),                      Stat=fail(1))
     End If
     Allocate (chbuf_read(1:batsz),iwrk_read(1:batsz),                        Stat=fail(2))
     Allocate (axx_read(1:batsz),ayy_read(1:batsz),azz_read(1:batsz),         Stat=fail(3))
     Allocate (scatter_buffer_read(1:wp_vals_per_at,1:batsz),                 Stat=fail(4))
     Allocate (chbuf_scat(1:batsz),iwrk_scat(1:batsz),                        Stat=fail(5))
     Allocate (n_held(0:comm%mxnode-1),where_buff(0:comm%mxnode-1),owner_read(1:batsz), Stat=fail(6))
     If (Any(fail(1:6) > 0)) Then
        Write(message,'(a)') 'defects_reference_read_parallel allocation failure 2'
        Call error(0,message)
     End If

  Else

! It is Illegal to pass unallocated allocatable arrays to routines.
! Therefore for arrays that are used by the mpi_scatterv calls
! below allocate them to zero size if they are not used on this core

     Allocate (scatter_buffer_read(1:0,1:0),   Stat=fail(1))
     Allocate (chbuf_scat(1:0),iwrk_scat(1:0), Stat=fail(2))
     Allocate (n_held(0:-1),where_buff(0:-1),  Stat=fail(3))
     If (Any(fail(1:3) > 0)) Then
        Write(message,'(a)') 'defects_reference_read_parallel allocation failure 3'
        Call error(0,message)
     End If

  End If

  Call invert(dfcts%celr,dfcts%rcell,det)

! Initialise domain localised atom counter (defects._module),
! dispatched atom counter and safe dispatch flag

  dfcts%nrefs =0
  indatm=0
  safe  =.true.

  Do k=1,megref

! Read in transmission arrays

     Readers_only: If (do_read .and. indatm == 0) Then
        to_read = Min(batsz,orig_first_at(my_read_proc_num+1)-first_at(my_read_proc_num))

        No_netCDF: If (io_read /= IO_READ_NETCDF) Then

           this_rec_buff = 0
           recs_to_read  = 0
           Do i=1,to_read
              If (this_rec_buff == 0) Then
                 recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at )
                 If (.not.fast) Then
                    Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, comm%ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,offset_kind)
                 End If
              End If

! Atom details

              this_rec_buff = this_rec_buff + 1
              Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                 record( j:j ) = rec_buff( j, this_rec_buff )
              End Do
              Call strip_blanks(record)
              Call get_word(record,word) ; chbuf_read(i)=word(1:8)
              If (l_ind) Then
                 Call get_word(record,word)
                 iwrk_read(i)=Nint(word_2_real(word,0.0_wp,l_str))
                 If (iwrk_read(i) /= 0) Then
                    iwrk_read(i)=Abs(iwrk_read(i))
                 Else
                    iwrk_read(i)=first_at(my_read_proc_num)+i
                 End If
              Else
                 iwrk_read(i)=first_at(my_read_proc_num)+i
              End If

              If (this_rec_buff == recs_to_read) Then
                 this_rec_buff = 0
                 recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 1 )
                 If (.not.fast) Then
                    Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, comm%ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,offset_kind)
                 End If
              End If

! Positions

              this_rec_buff = this_rec_buff + 1
              Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                 record( j:j ) = rec_buff( j, this_rec_buff )
              End Do
              Read(record, Fmt=*, End=100) axx_read(i),ayy_read(i),azz_read(i)
              If (this_rec_buff == recs_to_read) Then
                 this_rec_buff = 0
                 If (lvcfgr > 0) Then
                    recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 2 )
                    If (.not.fast) Then
                       Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                    Else
                       Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, comm%ierr )
                       rec_mpi_io = rec_mpi_io + Int(recs_to_read,offset_kind)
                    End If
                 End If
              End If

! Velocities

              If (lvcfgr > 0) Then
                 this_rec_buff = this_rec_buff + 1
                 If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
                    If (lvcfgr > 1) Then
                       recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 3 )
                       If (.not.fast) Then
                          Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                       Else
                          Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, comm%ierr )
                          rec_mpi_io = rec_mpi_io + Int(recs_to_read,offset_kind)
                       End If
                    End If
                 End If
              End If

! Forces

              If (lvcfgr > 1) Then
                 this_rec_buff = this_rec_buff + 1

                 If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
                 End If
              End If
           End Do

        Else

           If (to_read /= 0) Then
              frame = Int(top_skip,Kind(frame))

              Call io_nc_get_var( 'atomnames', fh, rec_buff, (/ first_at( my_read_proc_num ) + 1, frame /), (/ 8, to_read, 1 /) )
              Do i = 1, to_read
                 Do j = 1, Min( Len( chbuf_read ), Size( rec_buff, Dim = 1 ) )
                    chbuf_read( i )( j:j ) = rec_buff( j, i )
                 End Do
              End Do
              If (l_ind) Then
                 Call io_nc_get_var( 'indices', fh, iwrk_read , (/ first_at( my_read_proc_num ) + 1, frame /), (/ to_read, 1 /) )
              End If

              start = (/ 1, first_at( my_read_proc_num ) + 1, frame /)
              count = (/ 3, to_read, 1 /)

              Call io_get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
           End If

        End If No_netCDF

! Assign atoms positions in fractional coordinates to the correct domains
! (DD bounding) and broadcast them

        n_held=0
        Do i=1,to_read
           sxx=dfcts%rcell(1)*axx_read(i)+dfcts%rcell(4)*ayy_read(i)+dfcts%rcell(7)*azz_read(i)
           syy=dfcts%rcell(2)*axx_read(i)+dfcts%rcell(5)*ayy_read(i)+dfcts%rcell(8)*azz_read(i)
           szz=dfcts%rcell(3)*axx_read(i)+dfcts%rcell(6)*ayy_read(i)+dfcts%rcell(9)*azz_read(i)

! sxx,syy,szz are in [-0.5,0.5) interval and values as 0.4(9) may pose a problem

           sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
           syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
           szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

           axx_read(i)=sxx
           ayy_read(i)=syy
           azz_read(i)=szz

! assign domain coordinates (call for errors)

           ipx=Int((sxx+0.5_wp)*nprx_r)
           ipy=Int((syy+0.5_wp)*npry_r)
           ipz=Int((szz+0.5_wp)*nprz_r)

           idm=ipx+nprx*(ipy+npry*ipz)
           If (idm < 0 .or. idm > (comm%mxnode-1)) Call error(555)
           owner_read(i) = idm
           n_held(idm) = n_held(idm)+1
        End Do

        where_buff(0)=0
        Do i=1,comm%mxnode-1
           where_buff(i) = where_buff(i-1) + n_held(i-1)
        End Do

        Do i=1,to_read
           idm = where_buff(owner_read(i))
           idm = idm+1
           where_buff(owner_read(i)) = idm

           chbuf_scat(idm) = chbuf_read(i)
           iwrk_scat(idm)  = iwrk_read(i)

           scatter_buffer_read(1,idm) = axx_read(i)
           scatter_buffer_read(2,idm) = ayy_read(i)
           scatter_buffer_read(3,idm) = azz_read(i)
        End Do

     End If Readers_only

! Increase buffer counter and update first_at for
! the readers that have something left to read

     indatm = indatm+1
     If (do_read) Then
        If (first_at(my_read_proc_num) < first_at(my_read_proc_num+1)) &
             first_at(my_read_proc_num) = first_at(my_read_proc_num)+1
     End If

! Circulate configuration data to all nodes when transmission arrays are filled up
! Check against megref since at low processors counts (i.e. 1) batsz can be > megref

     Reorganize_buffer: If (indatm == batsz .or. (indatm > 0 .and. k == megref)) Then

        Do which_read_proc = 0 , n_read_procs_use-1
           If (orig_first_at(which_read_proc) >= megref) Exit ! for non-reading readers

           this_base_proc = which_read_proc * per_read_proc
           If (comm%idnode == this_base_proc) Then
              where_buff(0) = 0
              Do i=1,comm%mxnode-1
                 where_buff(i) = where_buff(i-1) + n_held(i-1)
              End Do
           End If

           Call gscatter(comm,n_held(:),n_loc,this_base_proc)

           Call gscatterv(comm,chbuf_scat(:),n_held(:),where_buff(:), &
                          chbuf(1:n_loc),this_base_proc)
           Call gscatterv(comm,iwrk_scat(:),n_held(:),where_buff(:), &
                          iwrk(1:n_loc),this_base_proc)
           Call gscatter_columns(comm,scatter_buffer_read(:,:),n_held(:), &
                                 where_buff(:), &
                                 scatter_buffer(1:wp_vals_per_at,1:n_loc), &
                                 this_base_proc)

! Assign atoms to correct domains

Dispatch:  Do i=1,n_loc
              dfcts%nrefs=dfcts%nrefs+1

! Check safety by the upper bound of: namr,indr,xr,yr,zr &

              If (dfcts%nrefs > Size( dfcts%xr )) Then
                 safe=.false.
                 Exit Dispatch
              End If

              dfcts%namr(dfcts%nrefs)=chbuf(i)
              dfcts%indr(dfcts%nrefs)=iwrk(i)

              dfcts%xr(dfcts%nrefs)=scatter_buffer(1,i)
              dfcts%yr(dfcts%nrefs)=scatter_buffer(2,i)
              dfcts%zr(dfcts%nrefs)=scatter_buffer(3,i)
           End Do Dispatch
        End Do

! Check if all is dispatched fine

        Call gcheck(comm,safe)
        If (.not.safe) Call error(556)

! Nullify dispatch counter

        indatm=0

     End If Reorganize_buffer

  End Do

! The last reader to check for EoFile in REFERENCE
! and if none is hit to call error to abort

  If (do_read) Then
     If (first_at(my_read_proc_num) == megref) Then
        recs_to_read = 1
        If (.not.fast) Then
           Read(Unit=nrefdt, Fmt=forma, Iostat=ierr) rec_buff( :, 1:recs_to_read )
        Else
           Call io_read_batch( fh, rec_mpi_io, 1, rec_buff, ierr )
        End If
        safe = (ierr /= 0)
     End If
  End If
  Call gcheck(comm,safe)
  If (.not.safe) Call error(58)

  If (do_read) Then
     Deallocate (rec_buff,                   Stat=fail(1))
     Deallocate (chbuf_read,iwrk_read,       Stat=fail(2))
     Deallocate (axx_read,ayy_read,azz_read, Stat=fail(3))
     Deallocate (owner_read,                 Stat=fail(4))
     If (Any(fail(1:4) > 0)) Then
        Write(message,'(a)') 'defects_reference_read_parallel deallocation failure 2'
        Call error(0,message)
     End If
  End If

  Deallocate (first_at,orig_first_at, Stat=fail(1))
  Deallocate (n_held,where_buff,      Stat=fail(2))
  Deallocate (chbuf,chbuf_scat,       Stat=fail(3))
  Deallocate (iwrk,iwrk_scat,         Stat=fail(4))
  Deallocate (scatter_buffer_read,    Stat=fail(5))
  Deallocate (scatter_buffer,         Stat=fail(6))
  If (Any(fail(1:6) > 0)) Then
     Write(message,'(a)') 'defects_reference_read_parallel deallocation failure 1'
     Call error(0,message)
  End If

  Return

! error exit for REFERNCE file read

100 Continue
  Call error(554)
End Subroutine defects_reference_read_parallel


Subroutine defects_reference_set_halo(cut,dfcts,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes for REFERENCE
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ),    Intent( In    ) :: cut
  Type( defects_type ), Intent( InOut ) :: dfcts
  Type( comms_type ),   Intent( InOut ) :: comm

  Integer           :: fail,nlx,nly,nlz,i,j,ia,ib
  Real( Kind = wp ) :: celprp(1:10),xdc,ydc,zdc

  Integer, Allocatable :: ixyz(:)
  Character( Len = 256 ) :: message

  fail=0
  Allocate (ixyz(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'defects_reference_set_halo allocation failure'
     Call error(0,message)
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! calculate link cell dimensions per node

  nlx=Int(celprp(7)/(cut*nprx_r))
  nly=Int(celprp(8)/(cut*npry_r))
  nlz=Int(celprp(9)/(cut*nprz_r))

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! link-cell widths in reduced space

  dfcts%cwx=1.0_wp/xdc
  dfcts%cwy=1.0_wp/ydc
  dfcts%cwz=1.0_wp/zdc

! Distance from the - edge of this domain

  dfcts%dxl=Nearest( (-0.5_wp+dfcts%cwx)+Real(idx,wp)*r_nprx , +1.0_wp)+zero_plus
  dfcts%dyl=Nearest( (-0.5_wp+dfcts%cwy)+Real(idy,wp)*r_npry , +1.0_wp)+zero_plus
  dfcts%dzl=Nearest( (-0.5_wp+dfcts%cwz)+Real(idz,wp)*r_nprz , +1.0_wp)+zero_plus

! Distance from the + edge of this domain

  dfcts%dxr=Nearest( (-0.5_wp-dfcts%cwx)+Real(idx+1,wp)*r_nprx , -1.0_wp)-zero_plus
  dfcts%dyr=Nearest( (-0.5_wp-dfcts%cwy)+Real(idy+1,wp)*r_npry , -1.0_wp)-zero_plus
  dfcts%dzr=Nearest( (-0.5_wp-dfcts%cwz)+Real(idz+1,wp)*r_nprz , -1.0_wp)-zero_plus

! Convert atomic positions from MD cell centred
! Cartesian coordinates to reduced space ones
! Populate the halo indicator array

  dfcts%nlrefs=dfcts%nrefs     ! No halo exists yet
  ixyz(1:dfcts%nlrefs)=0 ! Initialise halo indicator
  Do i=1,dfcts%nlrefs
     If (dfcts%xr(i) <= dfcts%dxl) ixyz(i)=ixyz(i)+1
     If (dfcts%xr(i) >= dfcts%dxr) ixyz(i)=ixyz(i)+2

     If (dfcts%yr(i) <= dfcts%dyl) ixyz(i)=ixyz(i)+10
     If (dfcts%yr(i) >= dfcts%dyr) ixyz(i)=ixyz(i)+20

     If (dfcts%zr(i) <= dfcts%dzl) ixyz(i)=ixyz(i)+100
     If (dfcts%zr(i) >= dfcts%dzr) ixyz(i)=ixyz(i)+200
  End Do

! exchange atom data in -/+ x directions

  Call defects_reference_export(-1,ixyz,dfcts,comm)
  Call defects_reference_export( 1,ixyz,dfcts,comm)

! exchange atom data in -/+ y directions

  Call defects_reference_export(-2,ixyz,dfcts,comm)
  Call defects_reference_export( 2,ixyz,dfcts,comm)

! exchange atom data in -/+ z directions

  Call defects_reference_export(-3,ixyz,dfcts,comm)
  Call defects_reference_export( 3,ixyz,dfcts,comm)

  Do i=1,dfcts%nlrefs
     dfcts%lri(i)=i
     dfcts%lra(i)=dfcts%indr(i)
  End Do
  Call shellsort2(dfcts%nlrefs,dfcts%lri,dfcts%lra)

  Do i=1,dfcts%nlrefs-1
     j=1
     Do While ((i+j) <= dfcts%nlrefs)
        If (dfcts%lra(i) == dfcts%lra(i+j)) Then
           ia=Min(dfcts%lri(i),dfcts%lri(i+j))
           ib=Max(dfcts%lri(i),dfcts%lri(i+j))
           dfcts%lri(i)=ia
           dfcts%lri(i+j)=ib
           j=j+1
        Else
           Exit
        End If
     End Do
  End Do

! Get domain halo limits in reduced space

! Distance from the - edge of this domain

  dfcts%dxl=Nearest( (-0.5_wp-dfcts%cwx)+Real(idx,wp)*r_nprx , -1.0_wp)-zero_plus
  dfcts%dyl=Nearest( (-0.5_wp-dfcts%cwy)+Real(idy,wp)*r_npry , -1.0_wp)-zero_plus
  dfcts%dzl=Nearest( (-0.5_wp-dfcts%cwz)+Real(idz,wp)*r_nprz , -1.0_wp)-zero_plus

! Distance from the + edge of this domain

  dfcts%dxr=Nearest( (-0.5_wp+dfcts%cwx)+Real(idx+1,wp)*r_nprx , +1.0_wp)+zero_plus
  dfcts%dyr=Nearest( (-0.5_wp+dfcts%cwy)+Real(idy+1,wp)*r_npry , +1.0_wp)+zero_plus
  dfcts%dzr=Nearest( (-0.5_wp+dfcts%cwz)+Real(idz+1,wp)*r_nprz , +1.0_wp)+zero_plus

  Deallocate (ixyz, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'defects_reference_set_halo deallocation failure'
     Call error(0,message)
  End If

End Subroutine defects_reference_set_halo


Subroutine defects_reference_write(name,megref,dfcts,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for re-writing the reference file needed for
!           defects detection
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Character( Len = * ), Intent( In    ) :: name
  Integer,              Intent( In    ) :: megref
  Type( defects_type),  Intent( InOut ) :: dfcts
  Type( comms_type),    Intent( InOut ) :: comm
  Integer, Parameter :: recsz = 73 ! default record size

  Logical               :: ready
  Integer(Kind=li)      :: recr,recr1 ! record line

  Integer               :: fail(1:4),i,k,jj,jdnode,jatms

  Real( Kind = wp )     :: celprp(1:10),x,y,z,cell_vecs(1:3,1:3)
  Real( Kind = wp )     :: lengths(1:3), angles(1:3)

! Some parameters and variables needed by io interfaces

  Integer                           :: fh
  Integer                           :: io_write,io_read,batsz,ierr
  Integer( Kind = offset_kind )     :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Character( Len = 256 ) :: message


  fail=0
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(1))
  If (fail(1) > 0) Then
     Write(message,'(a)') 'defects_reference_write allocation failure'
     Call error(0,message)
  End If

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_method_read       = io_read  )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! Nasty fix here as this is a REWRITE!!!

  If (io_read == IO_READ_NETCDF) io_write = IO_WRITE_SORTED_NETCDF

! Get offsets and define batch

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Allocate (n_atm(0:comm%mxnode),        Stat=fail(1))
     Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'write_config allocation failure 0'
        Call error(0,message)
     End If

     chbat=' '
     n_atm=0 ; n_atm(comm%idnode+1)=dfcts%nrefs
     Call gsum(comm,n_atm)
     n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = offset_kind)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(1-1,offset_kind)
     jj=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_delete( name ,comm) ! Sort existence issues
        Call io_open( io_write, comm_self, name, mode_wronly + mode_create, fh )

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(3i10,a42,a1)') 0,imcon,megref,Repeat(' ',42),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record, Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header

        Call io_write_batch( fh, rec_mpi_io, jj, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        jj=jj+5

     End If
     Call gsync(comm)! Start of file (updated)

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, name, mode_wronly, fh )

! Get to real space

     Do i=1,dfcts%nrefs
        axx(i)=cell(1)*dfcts%xr(i)+cell(4)*dfcts%yr(i)+cell(7)*dfcts%zr(i)
        ayy(i)=cell(2)*dfcts%xr(i)+cell(5)*dfcts%yr(i)+cell(8)*dfcts%zr(i)
        azz(i)=cell(3)*dfcts%xr(i)+cell(6)*dfcts%yr(i)+cell(9)*dfcts%zr(i)
     End Do

! DD bound

     Call pbcshift(imcon,cell,dfcts%nrefs,axx,ayy,azz)

! Start of file (updated)

     rec_mpi_io=Int(jj,offset_kind)+Int(n_atm(0),offset_kind)*Int(2,offset_kind)
     jj=0
     Do i=1,dfcts%nrefs
        Write(record, Fmt='(a8,i10,a54,a1)') dfcts%namr(i),dfcts%indr(i),Repeat(' ',54),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

! Dump batch and update start of file

        If (jj + 2 >= batsz .or. i == dfcts%nrefs) Then
           Call io_write_batch( fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,offset_kind)
           jj=0
        End If
     End Do

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'defects_reference_write allocation failure 1'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     recr=Int(0,li)
     jj=0
     If (comm%idnode == 0) Then

! Obliterate old REFERENCE if not needed and print header

        Open(Unit=nrefdt, File=name, Form='formatted', Access='direct', Recl=recsz, Status='replace')

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(3i10,a42,a1)') 0,imcon,megref,Repeat(' ',42),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record, Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=nrefdt, Fmt='(73a)', Rec=Int(1,li)) (chbat(:,k), k=1,jj)
        recr=Int(jj,li)
        jj=0

! Get to real space

        Do i=1,dfcts%nrefs
           chbuf(i)=dfcts%namr(i)
           iwrk(i)=dfcts%indr(i)

           axx(i)=cell(1)*dfcts%xr(i)+cell(4)*dfcts%yr(i)+cell(7)*dfcts%zr(i)
           ayy(i)=cell(2)*dfcts%xr(i)+cell(5)*dfcts%yr(i)+cell(8)*dfcts%zr(i)
           azz(i)=cell(3)*dfcts%xr(i)+cell(6)*dfcts%yr(i)+cell(9)*dfcts%zr(i)
        End Do

        jatms=dfcts%nrefs
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,DefRWrite_tag)

              Call grecv(comm,jatms,jdnode,DefRWrite_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,DefRWrite_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,DefRWrite_tag)

                 Call grecv(comm,axx(1:jatms),jdnode,DefRWrite_tag)
                 Call grecv(comm,ayy(1:jatms),jdnode,DefRWrite_tag)
                 Call grecv(comm,azz(1:jatms),jdnode,DefRWrite_tag)
              End If

! Get to real space

              Do i=1,jatms
                 x=axx(i) ; y=ayy(i) ; z=azz(i)
                 axx(i)=cell(1)*x+cell(4)*y+cell(7)*z
                 ayy(i)=cell(2)*x+cell(5)*y+cell(8)*z
                 azz(i)=cell(3)*x+cell(6)*y+cell(9)*z
              End Do
           End If

! DD bound

           Call pbcshift(imcon,cell,jatms,axx,ayy,azz)

           Do i=1,jatms
              Write(record, Fmt='(a8,i10,a54,a1)') chbuf(i),iwrk(i),Repeat(' ',54),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

! Dump batch and update start of file

              If (jj + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nrefdt, Fmt='(73a)', Rec=recr+Int(1,li)) (chbat(:,k), k=1,jj)
                 recr=recr+Int(jj,li)
                 jj=0
              End If
           End Do
        End Do

! Close REFERENCE

        Close(Unit=nrefdt)

     Else

        Call grecv(comm,ready,0,DefRWrite_tag)

        Call gsend(comm,dfcts%nrefs,0,DefRWrite_tag)
        If (dfcts%nrefs > 0) Then
           Call gsend(comm,dfcts%namr(1:dfcts%nrefs),0,DefRWrite_tag)
           Call gsend(comm,dfcts%indr(1:dfcts%nrefs),0,DefRWrite_tag)

           Call gsend(comm,dfcts%xr(1:dfcts%nrefs),0,DefRWrite_tag)
           Call gsend(comm,dfcts%yr(1:dfcts%nrefs),0,DefRWrite_tag)
           Call gsend(comm,dfcts%zr(1:dfcts%nrefs),0,DefRWrite_tag)
        End If

     End If

     Deallocate (chbuf,iwrk, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'defects_reference_write deallocation failure 1'
        Call error(0,message)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(1-1,offset_kind)
     jj=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_delete( name, comm ) ! Sort existence issues
        If (io_write == IO_WRITE_SORTED_NETCDF) &
        Call io_nc_create( comm_self, name, cfgname, megref )

        Call io_open( io_write, comm_self, name, mode_wronly + mode_create, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header

           Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
           Call io_write_record( fh, Int(jj,offset_kind), record )
           jj=jj+1

           Write(record, Fmt='(3i10,a42,a1)') 0,imcon,megref,Repeat(' ',42),lf
           Call io_write_record( fh, Int(jj,offset_kind), record )
           jj=jj+1

           Do i = 0, 2
              Write( record, '( 3f20.10, a12, a1 )' ) cell( 1 + i * 3: 3 + i * 3 ), Repeat( ' ', 12 ), lf
              Call io_write_record( fh, Int(jj,offset_kind), record )
              jj=jj+1
           End Do

        Else ! netCDF write

           jj=1 ! For config there is only one frame

           Call io_nc_put_var( 'time'           , fh, 0.0_wp, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,      0, jj, 1 )
           Call io_nc_put_var( 'datalevel'      , fh,      0, jj, 1 )
           Call io_nc_put_var( 'imageconvention', fh,  imcon, jj, 1 )
           Call io_nc_put_var( 'timestep'       , fh, 0.0_wp, jj, 1 )

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

     Else

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           jj=jj+5
        Else
           jj=1
        End If

     End If
     Call gsync(comm)

! Write the rest

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, name, mode_wronly, fh )

! Get to real space

     Do i=1,dfcts%nrefs
        axx(i)=cell(1)*dfcts%xr(i)+cell(4)*dfcts%yr(i)+cell(7)*dfcts%zr(i)
        ayy(i)=cell(2)*dfcts%xr(i)+cell(5)*dfcts%yr(i)+cell(8)*dfcts%zr(i)
        azz(i)=cell(3)*dfcts%xr(i)+cell(6)*dfcts%yr(i)+cell(9)*dfcts%zr(i)
     End Do

! DD bound

     Call pbcshift(imcon,cell,dfcts%nrefs,axx,ayy,azz)

! Write the rest

     rec_mpi_io=rec_mpi_io+Int(jj,offset_kind)
     Call io_write_sorted_file( fh, 0, IO_RESTART, rec_mpi_io, dfcts%nrefs,          &
          dfcts%indr, dfcts%namr, (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), axx, ayy, azz, &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /),                            &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), ierr )

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

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'defects_reference_write allocation failure 1'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     recr=Int(0,li)
     If (comm%idnode == 0) Then

! Obliterate old REFERENCE if not needed and print header

        Open(Unit=nrefdt, File=name, Form='formatted', Access='direct', Recl=73, Status='replace')
        recr=recr+Int(1,li)
        Write(Unit=nrefdt, Fmt='(a72,a1)',         Rec=recr) cfgname(1:72),lf
        recr=recr+Int(1,li)
        Write(Unit=nrefdt, Fmt='(3i10,a42,a1)',    Rec=recr) 0,imcon,megref,Repeat(' ',42),lf
        Do i = 0, 2
           recr=recr+Int(1,li)
           Write(Unit=nrefdt, Fmt='(3f20.10,a12,a1)', Rec=recr) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
        End Do

! Get to real space

        Do i=1,dfcts%nrefs
           chbuf(i)=dfcts%namr(i)
           iwrk(i)=dfcts%indr(i)

           axx(i)=cell(1)*dfcts%xr(i)+cell(4)*dfcts%yr(i)+cell(7)*dfcts%zr(i)
           ayy(i)=cell(2)*dfcts%xr(i)+cell(5)*dfcts%yr(i)+cell(8)*dfcts%zr(i)
           azz(i)=cell(3)*dfcts%xr(i)+cell(6)*dfcts%yr(i)+cell(9)*dfcts%zr(i)
        End Do

        jatms=dfcts%nrefs
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,DefRWrite_tag)

              Call grecv(comm,jatms,jdnode,DefRWrite_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,DefRWrite_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,DefRWrite_tag)

                 Call grecv(comm,axx(1:jatms),jdnode,DefRWrite_tag)
                 Call grecv(comm,ayy(1:jatms),jdnode,DefRWrite_tag)
                 Call grecv(comm,azz(1:jatms),jdnode,DefRWrite_tag)
              End If

! Get to real space

              Do i=1,jatms
                 x=axx(i) ; y=ayy(i) ; z=azz(i)
                 axx(i)=cell(1)*x+cell(4)*y+cell(7)*z
                 ayy(i)=cell(2)*x+cell(5)*y+cell(8)*z
                 azz(i)=cell(3)*x+cell(6)*y+cell(9)*z
              End Do
           End If

! DD bound

           Call pbcshift(imcon,cell,jatms,axx,ayy,azz)

           Do i=1,jatms
              recr1=recr+Int(iwrk(i)-1,li)*Int(2)+Int(1,li)
              Write(Unit=nrefdt, Fmt='(a8,i10,a54,a1)', Rec=recr1) chbuf(i),iwrk(i),Repeat(' ',54),lf

              recr1=recr1+Int(1,li)
              Write(Unit=nrefdt, Fmt='(3g20.10,a12,a1)', Rec=recr1) axx(i),ayy(i),azz(i),Repeat(' ',12),lf
           End Do
        End Do

! Close REFERENCE

        Close(Unit=nrefdt)

     Else

        Call grecv(comm,ready,0,DefRWrite_tag)

        Call gsend(comm,dfcts%nrefs,0,DefRWrite_tag)
        If (dfcts%nrefs > 0) Then
           Call gsend(comm,dfcts%namr(1:dfcts%nrefs),0,DefRWrite_tag)
           Call gsend(comm,dfcts%indr(1:dfcts%nrefs),0,DefRWrite_tag)

           Call gsend(comm,dfcts%xr(1:dfcts%nrefs),0,DefRWrite_tag)
           Call gsend(comm,dfcts%yr(1:dfcts%nrefs),0,DefRWrite_tag)
           Call gsend(comm,dfcts%zr(1:dfcts%nrefs),0,DefRWrite_tag)
         End If

     End If

     Deallocate (chbuf,iwrk, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'defects_reference_write deallocation failure 1'
        Call error(0,message)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a,i0)') 'write_config deallocation failure 0'
        Call error(0,message)
     End If
  End If

  Deallocate (axx,ayy,azz, Stat=fail(1))
  If (fail(1) > 0) Then
     Write(message,'(a)') 'defects_reference_write deallocation failure'
     Call error(0,message)
  End If

  Call gsync(comm)

End Subroutine defects_reference_write

!> defects_write
  Subroutine defects_write(keyres,ensemble,nstep,tstep,time,cshell,dfcts,neigh,site,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing DEFECTS file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,              Intent( In    ) :: keyres,ensemble,nstep
  Real( Kind = wp )   , Intent( In    ) :: tstep,time
  Type( defects_type ), Intent( InOut ) :: dfcts
  Type( neighbours_type ), Intent( In    ) :: neigh
  Type( site_type ), Intent( In    ) :: site
  Type( core_shell_type ), Intent( In    ) :: cshell
  Type( comms_type)   , Intent( InOut ) :: comm
  
  Integer, Parameter :: recsz = 73 ! default record size

  Logical                 :: safe,lexist,l_tmp,ready
  Character( Len = 40 )   :: word
  Integer                 :: fail(1:7),i,j,k,nlx,nly,nlz,   &
                             ic,jc,kk,ix,iy,iz,jx,jy,jz,    &
                             taken,ni,megni,nv,megnv,jdnode,jatms
  Real( Kind = wp )       :: cut,x,y,z,xs,ys,zs,buffer(1:2)

! Some parameters and variables needed by io interfaces

  Integer                           :: fh, io_write, batsz
  Integer( Kind = offset_kind )     :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

! Number of neighbouring cells to look around for counting defects

  Integer, Parameter :: nsbcll = 27
  Character( Len = 256 ) :: message

! Direction arrays for jumping around in link-cell space

  Integer, Dimension( 1:nsbcll ), Parameter :: &
  nix = (/ 0,  -1,-1,-1, 0, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 0, 1, 1, 1 /) , &
  niy = (/ 0,  -1, 0, 1,-1, 0, 1,-1, 0, 1, -1, 0, 1,-1, 1,-1, 0, 1, -1, 0, 1,-1, 0, 1,-1, 0, 1 /) , &
  niz = (/ 0,  -1,-1,-1,-1,-1,-1,-1,-1,-1,  0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1, 1 /)

  Real( Kind = wp ),    Dimension( : ),    Allocatable :: dr
  Character( Len = 8 ), Dimension( : ),    Allocatable :: namv,nami
  Integer,              Dimension( : ),    Allocatable :: indv,indi,interstitial,occupies
  Integer,              Dimension( : ),    Allocatable :: linkr,lctr,link,lct

  Integer,              Dimension( : ),    Allocatable :: ni_n,nv_n

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: cxx,cyy,czz

  If (.not.(nstep >= dfcts%nsdef .and. Mod(nstep-dfcts%nsdef,dfcts%isdef) == 0)) Return

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! netCDF not implemented for DEFECTS.  Switch to DEFAULT temporarily.

  If (io_write == IO_WRITE_SORTED_NETCDF) io_write = IO_WRITE_SORTED_MPIIO

    
   If (dfcts%newjob) Then
     dfcts%newjob = .false.

! rdef squared

     dfcts%rdefsq=dfcts%rdef**2

! Build lattice sites list from REFERENCE
     Call allocate_defects_arrays(dfcts)
     Call defects_reference_read(nstep,dfcts,site,comm)

! Assume that the MD cell will not change much in size and shape from
! the one provided in REFERENCE, a smaller halo(cutoff(rdef)) is to be set

     cut=dfcts%rdef+0.15_wp
     Call defects_reference_set_halo(cut,dfcts,comm)

! If the keyres=1, is DEFECTS old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File=trim(dfcts%deffile), Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (comm%idnode == 0) Then
           Open (Unit=ndefdt, File=trim(dfcts%deffile), Form='formatted', Access='direct', Status='replace', Recl=recsz)
           Write(Unit=ndefdt, Fmt='(a72,a1)',           Rec=1) cfgname(1:72),lf
           Write(Unit=ndefdt, Fmt='(f7.3,a23,2i21,a1)', Rec=2) dfcts%rdef,Repeat(' ',23),dfcts%frm,dfcts%rec,lf
           Close(Unit=ndefdt)
        End If
        dfcts%rec=Int(2,li)
        dfcts%frm=Int(0,li)

! Get some sense of it

     Else

        safe=.true.
        If (comm%idnode == 0) Then

           Open(Unit=ndefdt, File=trim(dfcts%deffile), Form='formatted')

           l_tmp =.true.
           Do

              record=' '
              If (l_tmp) Then

                 Read(Unit=ndefdt, Fmt=*, End=20)            ! title record
                 dfcts%rec=dfcts%rec+Int(1,li)
                 Read(Unit=ndefdt, Fmt='(a)', End=20) record ! bookkeeping record
                 dfcts%rec=dfcts%rec+Int(1,li)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 If (word(1:Len_Trim(word)) /= 'timestep') Then
                    Call get_word(record,word) ; Call get_word(record,word)
                    Call get_word(record,word) ; dfcts%frm=Nint(word_2_real(word,0.0_wp),li)
                    Call get_word(record,word) ; dfcts%rec=Nint(word_2_real(word,0.0_wp),li)
                    If (dfcts%frm /= Int(0,li) .and. dfcts%rec > Int(2,li)) Then
                       Go To 20 ! New style
                    Else
                       l_tmp=.false. ! TOUGH, old style
                       dfcts%rec=Int(2,li)
                       dfcts%frm=Int(0,li)
                    End If
                 Else
                    safe=.false. ! Overwrite the file, it's junk to me
                    Go To 20
                 End If

              Else

                 Read(Unit=ndefdt, Fmt=*, End=20)            ! timestep record
                 dfcts%rec=dfcts%rec+Int(1,li)

                 Read(Unit=ndefdt, Fmt='(a)', End=20) record ! defects record
                 dfcts%rec=dfcts%rec+Int(1,li)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 Call get_word(record,word) ; j=Nint(word_2_real(word))

                 Do i=1,3+2*j ! 3 lines for cell parameters and 2*j entries for defects
                    Read(Unit=ndefdt, Fmt=*, End=20)
                    dfcts%rec=dfcts%rec+Int(1,li)
                 End Do
                 dfcts%frm=dfcts%frm+Int(1,li)

              End If

           End Do

20         Continue
           Close(Unit=ndefdt)

        End If

        Call gcheck(comm,safe,"enforce")
        If (.not.safe) Then
           lexist=.false.

           dfcts%rec=Int(0,li)
           dfcts%frm=Int(0,li)

           Go To 10
        Else If (comm%mxnode > 1) Then
           buffer(1)=Real(dfcts%frm,wp)
           buffer(2)=Real(dfcts%rec,wp)

           Call gbcast(comm,buffer(1:2),0)

           dfcts%frm=Nint(buffer(1),li)
           dfcts%rec=Nint(buffer(2),li)
        End If

     End If

! Get rcell

     Call invert(cell,dfcts%rcell,cut)

! New real space cutoff and expansion for defects link-cells
     dfcts%cutdef=Min(neigh%cutoff/3.0_wp,2.0_wp*dfcts%rdef)
     dfcts%mxlcdef=Nint(((neigh%cutoff/dfcts%cutdef)**3+0.15_wp)*Real(neigh%max_cell,wp))
  End If

! Update rcell

  If (ensemble >= 20) Call invert(cell,dfcts%rcell,cut)

  fail=0
  Allocate (dr(1:mxatms),                                                Stat=fail(1))
  Allocate (namv(1:mxatms),indv(1:mxatms),nami(1:mxatms),indi(1:mxatms), Stat=fail(2))
  Allocate (interstitial(1:mxatms),occupies(1:mxatms),                   Stat=fail(3))
  Allocate (linkr(1:mxatms),link(1:mxatms),                              Stat=fail(4))
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms),                   Stat=fail(5))
  Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms),                   Stat=fail(6))
  Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms),                   Stat=fail(7))
  If (Any(fail > 0)) Then
     Write(message,'(a)') Trim(dfcts%deffile)//'_write allocation failure'
     Call error(0,message)
  End If

! Build bookkeeping lists: interstitial, occupies

  Do i=1,nlast

! Consider all particles at this point:j=1 - consider in,j=0 - consider out

     j=1

! Get all domain+halo particles' coordinates in MD cell centred reduced space

     cxx(i)=dfcts%rcell(1)*xxx(i)+dfcts%rcell(4)*yyy(i)+dfcts%rcell(7)*zzz(i)
     cyy(i)=dfcts%rcell(2)*xxx(i)+dfcts%rcell(5)*yyy(i)+dfcts%rcell(8)*zzz(i)
     czz(i)=dfcts%rcell(3)*xxx(i)+dfcts%rcell(6)*yyy(i)+dfcts%rcell(9)*zzz(i)

! Exclude particles in the domain's halo farther than cutoff(rdef)
! smaller halo with a width defined by cutoff(rdef)

     If ( i > natms .and.                        &
          ((cxx(i) < dfcts%dxl .or. cxx(i) > dfcts%dxr) .or. &
           (cyy(i) < dfcts%dyl .or. cyy(i) > dfcts%dyr) .or. &
           (czz(i) < dfcts%dzl .or. czz(i) > dfcts%dzr)) ) j=0

! Exclude frozen and shell particles from consideration

     If ( j == 1 .and. (lfrzn(i) /= 0 .or. Any(cshell%listshl(2,1:cshell%ntshl) == ltg(i))) ) j=0

! Assume that every considered particles (1) is an interstitial
! and (2) does not occupy a site yet

     If (j == 1) Then
        interstitial(i)=-1 ! considered particle is assumed to be an interstitial
        occupies(i)    = 0 ! and to not occupy a site yet (but may occupy one later)
     Else
        interstitial(i)= 0 ! excluded particle can neither ever be an interstitial
        occupies(i)    =-1 ! nor can it ever occupy a site
     End If
  End Do

! Partition sites and atoms in link-celled space with same imcon, cell and cutoff!!!

  Allocate (lctr(0:dfcts%mxlcdef),lct(0:dfcts%mxlcdef), Stat=fail(1))
  If (fail(1) > 0) Then
     Write(message,'(a)') Trim(dfcts%deffile)//'_write allocation failure 1'
     Call error(0,message)
  EndIf
  Call defects_link_cells &
           (dfcts%cutdef,dfcts%mxlcdef,dfcts%nrefs,dfcts%nlrefs,dfcts%xr,dfcts%yr,dfcts%zr,nlx,nly,nlz,linkr,lctr)
  Call defects_link_cells &
           (dfcts%cutdef,dfcts%mxlcdef,natms,nlast,cxx,cyy,czz,nlx,nly,nlz,link,lct)

  safe = .true.          ! Initialise safety flag to all safe
  nv = 0                 ! Assume no vacancies (all sites are occupied)
  ni = 0                 ! Assume no interstitials (all particles are occupying a site)
  dr = (dfcts%rdef+0.15_wp)**2 ! Assume all particles are > rdef distance away from the nearest to them site

! Actual defect detection *********************************************

! primary loop over domain subcells of sites

  Do iz=0,nlz+1
     Do iy=0,nly+1
        Do ix=0,nlx+1

! index of primary reference link-cell

           ic=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! bypass reference subcell if empty

           If (lctr(ic) > 0) Then

! head of chain of i-th subcell (reference)

              i=lctr(ic)

! loop over primary cell contents

100           Continue

! Bypass if the site is a shell

              If (Any(cshell%listshl(2,1:cshell%ntshl) == ltg(i))) Go To 400

! Assume the site is vacant

              taken=0 ! Assume the site is vacant

! secondary loop over subcells of particles

              Do kk=1,nsbcll

                 jx=ix+nix(kk)
                 jy=iy+niy(kk)
                 jz=iz+niz(kk)

! SHIFT BACK to the LEFT
! Disregard cells outside the look-up scope

  If ( (jx >= 0) .and. (jx <= nlx+1) .and. &
       (jy >= 0) .and. (jy <= nly+1) .and. &
       (jz >= 0) .and. (jz <= nlz+1) ) Then

! index of neighbouring cell

     jc=1+jx+(nlx+2)*(jy+(nly+2)*jz)

! bypass real subcell if empty

     If (lct(jc) > 0) Then

! head of chain of j-th subcell (real)

        j=lct(jc)

! loop over secondary cell contents

200     Continue

! Bypass if the atom is excluded or has already claimed (occupied) a site

        If (occupies(j) /= 0) Go To 300

! Get vectors particle-site (in reduced space)

        xs=dfcts%xr(i)-cxx(j) ; If (Abs(xs) > dfcts%cwx) Go To 300
        ys=dfcts%yr(i)-cyy(j) ; If (Abs(ys) > dfcts%cwy) Go To 300
        zs=dfcts%zr(i)-czz(j) ; If (Abs(zs) > dfcts%cwz) Go To 300

! Get in real space

        x=cell(1)*xs+cell(4)*ys+cell(7)*zs
        y=cell(2)*xs+cell(5)*ys+cell(8)*zs
        z=cell(3)*xs+cell(6)*ys+cell(9)*zs

! Get particle-site squared distance

        cut=x**2+y**2+z**2

! An atom qualifies to claim a site (occupies > 0)
! when the site is vacant - not occupied/taken (taken=0).
! If the site is not claimed yet then the atom is marked
! as a claimee of this site (occupies > 0) and not an
! interstitial (interstitial = 0).  If the site is
! already taken by another atom then the atom is marked
! as an interstitial of this site (interstitial > 0).

        If (cut < dfcts%rdefsq) Then
           dr(j)=cut             ! Save the distance (site-atom)

! Only general interstitial atoms are allowed to claim.  Otherwise, we are in a big trouble

           If (interstitial(j) >= 0) safe=.false.

           If (taken == 0) Then
              taken=1            ! This site is taken (no longer vacant)
              occupies(j)=i      ! This atom claims a site (this site)
              interstitial(j)=0  ! This atom is not an interstitial any longer
           Else
              interstitial(j)=i  ! This atom is not a general interstitial any longer
           End If
        End If

300     Continue

        j=link(j)
        If (j /= 0) Go To 200

! end of loop over the real subcell contents jc

     End If

! End If of disregard cell outside the look-up scope

  End If

! SHIFT BACK to the RIGHT
! End Do secondary loop around real subcells

              End Do

! Vacancies: taken=0 means we've found a vacancy

              If (taken == 0 .and. i <= dfcts%nrefs) Then
                 nv=nv+1
                 namv(nv)=dfcts%namr(i)
                 indv(nv)=dfcts%indr(i)
                 axx(nv)=cell(1)*dfcts%xr(i)+cell(4)*dfcts%yr(i)+cell(7)*dfcts%zr(i)
                 ayy(nv)=cell(2)*dfcts%xr(i)+cell(5)*dfcts%yr(i)+cell(8)*dfcts%zr(i)
                 azz(nv)=cell(3)*dfcts%xr(i)+cell(6)*dfcts%yr(i)+cell(9)*dfcts%zr(i)
              End If

400           Continue
              i=linkr(i)
              If (i /= 0) Go To 100

! end of bypass of empty subcell ic

           End If
        End Do
     End Do
  End Do

! Check safety on rdef length

  Call gcheck(comm,safe)
  If (.not.safe) Call error(560)

! Define real interstitials and real claimees of a site by swapping
! interstitial/occupation statuses between particles based on
! site-atom distances of these.  Amongst all atoms close to sites -
! interstitials of the site and the site claimee -
! the real claimee is the atom closest to the site.

! primary loop over domain subcells of particles

  Do iz=0,nlz+1
     Do iy=0,nly+1
        Do ix=0,nlx+1

! index of primary link-cell

           ic=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! bypass subcell if empty

           If (lct(ic) > 0) Then

! head of chain of i-th subcell

              i=lct(ic)

! loop over primary cell contents

500           Continue

! Bypass if the atom is excluded or is not a claimee

              If (occupies(i) <= 0) Go To 700

! Circle around domain+halo particles, k is the index of the
! closest to a site atom (real claimee) assuming that the first
! claimee is the true one.  Check for an interstitial of the
! same site, that is closer to it than the claimee.

              k=i
              cut=dr(i)

! secondary loop over subcells of particles

              Do kk=1,nsbcll

                 jx=ix+nix(kk)
                 jy=iy+niy(kk)
                 jz=iz+niz(kk)

! SHIFT BACK to the LEFT
! Disregard cells outside the look-up scope

  If ( (jx >= 0) .and. (jx <= nlx+1) .and. &
       (jy >= 0) .and. (jy <= nly+1) .and. &
       (jz >= 0) .and. (jz <= nlz+1) ) Then

! index of neighbouring cell

     jc=1+jx+(nlx+2)*(jy+(nly+2)*jz)

! bypass real subcell if empty

     If (lct(jc) > 0) Then

! head of chain of j-th subcell (real)

        j=lct(jc)

! loop over secondary cell contents

600     Continue

        If (interstitial(j) == occupies(i) .and. dr(j) < cut) Then
           cut=dr(j)
           k=j
        End If

        j=link(j)
        If (j /= 0) Go To 600

! end of loop over the subcell contents jc

     End If

! End If of disregard cell outside the look-up scope

  End If

! SHIFT BACK to the RIGHT
! End Do secondary loop around subcells

              End Do

! If there is a new true claimee change statuses

              If (k /= i) Then
                 interstitial(i)=interstitial(k) ; interstitial(k)=0
                 occupies(k)=occupies(i)         ; occupies(i)=0
              End If

700           Continue

              i=link(i)
              If (i /= 0) Go To 500

! end of bypass of empty subcell ic

           End If
        End Do
     End Do
  End Do

! Include the following if you want occupant/site type mismatches printed out
! printing from many nodes can be messy
!
!     taken = 0
!     Do i=1,natms
!        If (occupies(i) > 0 .and. atmnam(i) /= namr(occupies(i))) Then
!!           Write(message,'(a,2(1x,a))') 'occupant/site type mismatch:', atmnam(i), namr(occupies(i))
!            Call info(message)
!           taken = taken + 1
!        End If
!     End Do
!     Call gsum(taken)
!     Write(message,'(a,i10,2(1x,a,i10))') 'occupant/site type mismatches', taken
!     Call info(message,.true.)

! Interstitials: i <= natms & interstitial(i) /= 0 means we've found an interstitial

  Do i=1,natms
     If (interstitial(i) /= 0) Then
        ni=ni+1
        nami(ni)=atmnam(i)
        indi(ni)=ltg(i)
        bxx(ni)=xxx(i)
        byy(ni)=yyy(i)
        bzz(ni)=zzz(i)
     End If
  End Do

! Sum global defects values

  megni = ni
  Call gsum(comm,megni)

  megnv = nv
  Call gsum(comm,megnv)

! Get relative offsets (int,vac) for parallel printing

  Allocate (ni_n(0:comm%mxnode),nv_n(0:comm%mxnode), Stat=fail(1))
  Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') Trim(dfcts%deffile)//'_write allocation failure 2'
     Call error(0,message)
  End If

  ni_n=0 ; ni_n(comm%idnode+1)=ni
  Call gsum(comm,ni_n)
  ni_n(0)=Sum(ni_n(0:comm%idnode))

  nv_n=0 ; nv_n(comm%idnode+1)=nv
  Call gsum(comm,nv_n)
  nv_n(0)=Sum(nv_n(0:comm%idnode))

  chbat=' '

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = offset_kind)

! Update frame

  dfcts%frm=dfcts%frm+Int(1,li)

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(dfcts%rec,offset_kind)
     j=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_open( io_write, comm_self, trim(dfcts%deffile), mode_wronly, fh )

        Write(record, Fmt='(a8,i10,2f20.6,i5,f7.3,a2,a1)') &
           'timestep',nstep,tstep,time,imcon,dfcts%rdef,Repeat(' ',2),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(a8,i10,a15,i10,a11,i10,a8,a1)') &
           'defects ',megni+megnv, ' interstitials ',megni, ' vacancies ',megnv,Repeat(' ',8),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Do i = 0, 2
           Write( record, '( 3f20.10, a12, a1 )' ) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           j=j+1
           Do k=1,recsz
              chbat(k,j) = record(k:k)
           End Do
        End Do

! Dump header and cell information

        Call io_write_batch( fh, rec_mpi_io, j, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        j=j+5

     End If
     Call gsync(comm)

! Start of file

     dfcts%rec=dfcts%rec+Int(j,li)

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, trim(dfcts%deffile), mode_wronly, fh )

! Start of file

     rec_mpi_io=Int(dfcts%rec,offset_kind)+Int(2,offset_kind)*Int(ni_n(0),offset_kind)
     j=0
     Do i=1,ni
        Write(record, Fmt='(a2,a8,i10,a52,a1)') 'i_',nami(i),indi(i),Repeat(' ',52),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

! Dump batch and update start of file

        If (j + 2 >= batsz .or. i == ni) Then
           Call io_write_batch( fh, rec_mpi_io, j, chbat )
           rec_mpi_io=rec_mpi_io+Int(j,offset_kind)
           j=0
        End If
     End Do

! Start of file

     dfcts%rec=dfcts%rec+Int(2,li)*Int(megni,li)

! Start of file

     rec_mpi_io=Int(dfcts%rec,offset_kind)+Int(2,offset_kind)*Int(nv_n(0),offset_kind)
     Do i=1,nv
        Write(record, Fmt='(a2,a8,i10,a52,a1)') 'v_',namv(i),indv(i),Repeat(' ',52),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

! Dump batch and update start of file

        If (j + 2 >= batsz .or. i == nv) Then
           Call io_write_batch( fh, rec_mpi_io, j, chbat )
           rec_mpi_io=rec_mpi_io+Int(j,offset_kind)
           j=0
        End If
     End Do
     dfcts%rec=dfcts%rec+Int(2,li)*Int(megnv,li)

! Update main header

     If (comm%idnode == 0) Then
        Write(record, Fmt='(f7.3,a23,2i21,a1)') dfcts%rdef,Repeat(' ',23),dfcts%frm,dfcts%rec,lf
        Call io_write_record( fh, Int(1,offset_kind), record )
     End If

     Call io_close( fh )
     Call io_finalize

  Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
           io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') Trim(dfcts%deffile)//'_write allocation failure 3'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     j=0
     If (comm%idnode == 0) Then
        Open(Unit=ndefdt, File=trim(dfcts%deffile), Form='formatted', Access='direct', Recl=recsz)

! Accumulate header

        Write(record, Fmt='(a8,i10,2f20.6,i5,f7.3,a2,a1)') &
           'timestep',nstep,tstep,time,imcon,dfcts%rdef,Repeat(' ',2),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(a8,i10,a15,i10,a11,i10,a8,a1)') &
           'defects ',megni+megnv, ' interstitials ',megni, ' vacancies ',megnv,Repeat(' ',8),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record, Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           j=j+1
           Do k=1,recsz
              chbat(k,j) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=ndefdt, Fmt='(73a)', Rec=dfcts%rec+Int(1,li)) (chbat(:,k), k=1,j)
        dfcts%rec=dfcts%rec+Int(j,li)
        j=0

        Do i=1,ni
           chbuf(i)=nami(i)
           iwrk(i)=indi(i)

           cxx(i)=bxx(i)
           cyy(i)=byy(i)
           czz(i)=bzz(i)
        End Do

        jatms=ni
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,DefWrite_tag)

              Call grecv(comm,jatms,jdnode,DefWrite_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,DefWrite_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,DefWrite_tag)

                 Call grecv(comm,cxx(1:jatms),jdnode,DefWrite_tag)
                 Call grecv(comm,cyy(1:jatms),jdnode,DefWrite_tag)
                 Call grecv(comm,czz(1:jatms),jdnode,DefWrite_tag)
              End If
           End If

           Do i=1,jatms
              Write(record, Fmt='(a2,a8,i10,a52,a1)') 'i_',chbuf(i),iwrk(i),Repeat(' ',52),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

! Dump batch and update start of file

              If (j + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=ndefdt, Fmt='(73a)', Rec=dfcts%rec+Int(1,li)) (chbat(:,k), k=1,j)
                 dfcts%rec=dfcts%rec+Int(j,li)
                 j=0
              End If
           End Do
        End Do
     Else
        Call grecv(comm,ready,0,DefWrite_tag)

        Call gsend(comm,ni,0,DefWrite_tag)
        If (ni > 0) Then
           Call gsend(comm,nami(1:ni),0,DefWrite_tag)
           Call gsend(comm,indi(1:ni),0,DefWrite_tag)

           Call gsend(comm,bxx(1:ni),0,DefWrite_tag)
           Call gsend(comm,byy(1:ni),0,DefWrite_tag)
           Call gsend(comm,bzz(1:ni),0,DefWrite_tag)
        End If
     End If

     If (comm%idnode == 0) Then
        Do i=1,nv
           chbuf(i)=namv(i)
           iwrk(i)=indv(i)

           cxx(i)=axx(i)
           cyy(i)=ayy(i)
           czz(i)=azz(i)
        End Do

        jatms=nv
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,DefWrite_tag)

              Call grecv(comm,jatms,jdnode,DefWrite_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,DefWrite_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,DefWrite_tag)

                 Call grecv(comm,cxx(1:jatms),jdnode,DefWrite_tag)
                 Call grecv(comm,cyy(1:jatms),jdnode,DefWrite_tag)
                 Call grecv(comm,czz(1:jatms),jdnode,DefWrite_tag)
              End If
           End If

           Do i=1,jatms
              Write(record, Fmt='(a2,a8,i10,a52,a1)') 'v_',chbuf(i),iwrk(i),Repeat(' ',52),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

! Dump batch and update start of file

              If (j + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=ndefdt, Fmt='(73a)', Rec=dfcts%rec+Int(1,li)) (chbat(:,k), k=1,j)
                 dfcts%rec=dfcts%rec+Int(j,li)
                 j=0
              End If
           End Do
        End Do
     Else
        Call grecv(comm,ready,0,DefWrite_tag)

        Call gsend(comm,nv,0,DefWrite_tag)
        If (nv > 0) Then
           Call gsend(comm,namv(1:nv),0,DefWrite_tag)
           Call gsend(comm,indv(1:nv),0,DefWrite_tag)

           Call gsend(comm,axx(1:nv),0,DefWrite_tag)
           Call gsend(comm,ayy(1:nv),0,DefWrite_tag)
           Call gsend(comm,azz(1:nv),0,DefWrite_tag)
        End If
     End If

! Update main header

     If (comm%idnode == 0) Then
        Write(Unit=ndefdt, Fmt='(f7.3,a23,2i21,a1)', Rec=2) dfcts%rdef,Repeat(' ',23),dfcts%frm,dfcts%rec,lf
        Close(Unit=ndefdt)
     Else
        dfcts%rec=dfcts%rec+Int(5,li)+Int(2,li)*Int(megni+megnv,li)
     End If

     Deallocate (chbuf,iwrk, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') Trim(dfcts%deffile)//'_write deallocation failure 3'
        Call error(0,message)
     End If

  End If

  Call gsync(comm)

  Deallocate (ni_n,nv_n, Stat=fail(1))
  Deallocate (chbat,     Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') Trim(dfcts%deffile)//'_write deallocation failure 2'
     Call error(0,message)
  End If

  Deallocate (dr,                    Stat=fail(1))
  Deallocate (namv,indv,nami,indi,   Stat=fail(2))
  Deallocate (interstitial,occupies, Stat=fail(3))
  Deallocate (linkr,lctr,link,lct,   Stat=fail(4))
  Deallocate (axx,ayy,azz,           Stat=fail(5))
  Deallocate (bxx,byy,bzz,           Stat=fail(6))
  Deallocate (cxx,cyy,czz,           Stat=fail(7))
  If (Any(fail > 0)) Then
     Write(message,'(a)') Trim(dfcts%deffile)//'_write deallocation failure'
     Call error(0,message)
  End If


End Subroutine defects_write


End Module defects
