Module configuration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global configuration variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
! contrib   - i.j.bush march 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,li
  Use comms, Only : comms_type,wp_mpi,gbcast,WriteConf_tag,gcheck,gsync,gsum,&
                    gmax,gmin
  Use site

  Use setup,   Only : nconf,nrite,config,mxatms,half_minus,mxrgd,zero_plus, &
                      mxatdm,mxexcl,mxlist
  Use parse,   Only : tabs_2_blanks, &
                             strip_blanks, get_word, word_2_real,get_line
  Use domains, Only : nprx,npry,nprz,nprx_r,npry_r,nprz_r,idx,idy,idz
  Use development, Only : lvcfscl,cels,lvcforg,xorg,yorg,zorg

  Use io,     Only : io_set_parameters,         &
                            io_get_parameters,         &
                            io_init, io_nc_create,     &
                            io_open, io_write_record,  &
                            io_write_batch,            &
                            io_nc_put_var,             &
                            io_write_sorted_file,      &
                            io_delete,                 &
                            io_close, io_finalize,     &
                            io_read_batch,             &
                            io_nc_get_dim,             &
                            io_nc_get_var,             &
                            io_nc_get_att,             &
                            IO_READ_MASTER,            &
                            IO_READ_NETCDF,            &
                            IO_RESTART,                &
                            IO_BASE_COMM_NOT_SET,      &
                            IO_ALLOCATION_ERROR,       &
                            IO_UNKNOWN_WRITE_OPTION,   &
                            IO_UNKNOWN_WRITE_LEVEL,    &
                            IO_WRITE_UNSORTED_MPIIO,   &
                            IO_WRITE_UNSORTED_DIRECT,  &
                            IO_WRITE_UNSORTED_MASTER,  &
                            IO_WRITE_SORTED_MPIIO,     &
                            IO_WRITE_SORTED_DIRECT,    &
                            IO_WRITE_SORTED_NETCDF,    &
                            IO_WRITE_SORTED_MASTER

  Use errors_warnings, Only : error,warning,info
  use numerics, Only : shellsort2,invert,dcell,images,shellsort,pbcshift

#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif
  Implicit None

  Character( Len = 72 ), Save :: cfgname = ' ' , &
                                 sysname = ' '

  Integer,               Save :: imcon =-1 , &
                                 imc_n =-1 , &
                                 natms = 0 , &
                                 nlast = 0 , &
                                 nfree = 0

  Real( Kind = wp ),     Save :: cell(1:9) = 0.0_wp , &
                                 volm      = 0.0_wp , &
                                 sumchg    = 0.0_wp


  Character( Len = 8 ), Allocatable, Save :: atmnam(:)

  Integer,              Allocatable, Save :: lsite(:),ltype(:)
  Integer,              Allocatable, Save :: lfrzn(:),lfree(:)
  Integer,              Allocatable, Save :: lsi(:),lsa(:),ltg(:)
  Integer,              Allocatable, Save :: lexatm(:,:)
  Integer,              Allocatable, Save :: ixyz(:),list(:,:)
  Integer,              Allocatable, Save :: lstfre(:)

  Real( Kind = wp ),    Allocatable, Save :: weight(:),chge(:)
  Real( Kind = wp ),    Allocatable, Save :: xxx(:),yyy(:),zzz(:)
  Real( Kind = wp ),    Allocatable, Save :: vxx(:),vyy(:),vzz(:)
  Real( Kind = wp ),    Allocatable, Save :: fxx(:),fyy(:),fzz(:)

  Public :: reallocate, allocate_config_arrays_read, allocate_config_arrays
  Public :: check_config
  Public :: read_config_parallel
  Public :: scan_config

  Interface reallocate
     Module Procedure reallocate_chr_v
     Module Procedure reallocate_int_v
     Module Procedure reallocate_rwp_v
  End Interface

  Private :: reallocate_chr_v, reallocate_int_v, reallocate_rwp_v

Contains

  Subroutine reallocate_chr_v( delta, a, stat )

    Integer,                           Intent( In    ) :: delta
    Character( Len = * ), Allocatable, Intent( InOut ) :: a(:)
    Integer,                           Intent(   Out ) :: stat

    Integer :: size_old, size_new, size_crs

    Character( Len = Len( a ) ), Allocatable :: tmp(:)

    stat = 0
    If ( delta == 0 ) Return

    size_old = Size( a )
    size_new = size_old + delta
    size_crs = Min( size_old, size_new )

    Allocate ( tmp( 1:size_old ), Stat = stat )
    If ( stat /= 0 ) Return

    tmp = a

    Deallocate ( a, Stat = stat )
    If ( stat /= 0 ) Return

    Allocate ( a( 1:size_new ), Stat = stat )
    If ( stat /= 0 ) Return

    If ( size_crs > 0 ) Then
       a( 1:size_crs ) = tmp( 1:size_crs )
       If ( delta > 0 ) a( size_crs+1:size_new ) = ' '
    End If

    Deallocate ( tmp, Stat = stat )

  End Subroutine reallocate_chr_v

  Subroutine reallocate_int_v( delta, a, stat )

    Integer,              Intent( In    ) :: delta
    Integer, Allocatable, Intent( InOut ) :: a(:)
    Integer,              Intent(   Out ) :: stat

    Integer :: size_old, size_new, size_crs

    Integer, Allocatable :: tmp(:)

    stat = 0
    If ( delta == 0 ) Return

    size_old = Size( a )
    size_new = size_old + delta
    size_crs = Min( size_old, size_new )

    Allocate ( tmp( 1:size_old ), Stat = stat )
    If ( stat /= 0 ) Return

    tmp = a

    Deallocate ( a, Stat = stat )
    If ( stat /= 0 ) Return

    Allocate ( a( 1:size_new ), Stat = stat )
    If ( stat /= 0 ) Return

    If ( size_crs > 0 ) Then
       a( 1:size_crs ) = tmp( 1:size_crs )
       If ( delta > 0 ) a( size_crs+1:size_new ) = 0
    End If

    Deallocate ( tmp, Stat = stat )

  End Subroutine reallocate_int_v

  Subroutine reallocate_rwp_v( delta, a, stat )

    Integer,                        Intent( In    ) :: delta
    Real( Kind = wp ), Allocatable, Intent( InOut ) :: a(:)
    Integer,                        Intent(   Out ) :: stat

    Integer :: size_old, size_new, size_crs

    Real( Kind = wp ), Allocatable :: tmp(:)

    stat = 0
    If ( delta == 0 ) Return

    size_old = Size( a )
    size_new = size_old + delta
    size_crs = Min( size_old, size_new )

    Allocate ( tmp( 1:size_old ), Stat = stat )
    If ( stat /= 0 ) Return

    tmp = a

    Deallocate ( a, Stat = stat )
    If ( stat /= 0 ) Return

    Allocate ( a( 1:size_new ), Stat = stat )
    If ( stat /= 0 ) Return

    If ( size_crs > 0 ) Then
       a( 1:size_crs ) = tmp( 1:size_crs )
       If ( delta > 0 ) a( size_crs+1:size_new ) = 0.0_wp
    End If

    Deallocate ( tmp, Stat = stat )

  End Subroutine reallocate_rwp_v


  Subroutine allocate_config_arrays_read( isize )

    Integer, Intent( In    ) :: isize

    Integer :: fail(1:5)

    fail = 0

    Allocate (atmnam(1:isize),                        Stat = fail(1))
    Allocate (lsi(1:isize),lsa(1:isize),ltg(1:isize), Stat = fail(2))
    Allocate (xxx(1:isize),yyy(1:isize),zzz(1:isize), Stat = fail(3))
    Allocate (vxx(1:isize),vyy(1:isize),vzz(1:isize), Stat = fail(4))
    Allocate (fxx(1:isize),fyy(1:isize),fzz(1:isize), Stat = fail(5))

    If (Any(fail > 0)) Call error(1025)

    atmnam = ' '
    lsi = 0 ; lsa = 0 ; ltg = 0

    xxx = 0.0_wp ; yyy = 0.0_wp ; zzz = 0.0_wp
    vxx = 0.0_wp ; vyy = 0.0_wp ; vzz = 0.0_wp
    fxx = 0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp

  End Subroutine allocate_config_arrays_read

  Subroutine allocate_config_arrays()


    Integer           :: fail(1:5),stat(1:13)

    fail = 0

    Allocate (lsite(1:mxatms),ltype(1:mxatms),           Stat = fail(1))
    Allocate (lfrzn(1:mxatms),lfree(1:mxatms),           Stat = fail(2))
    Allocate (lexatm(0:mxexcl,1:mxatdm),ixyz(1:mxatms),  Stat = fail(3))
    Allocate (list(-3:mxlist,1:mxatdm),lstfre(1:mxatdm), Stat = fail(4))
    Allocate (weight(1:mxatms),chge(1:mxatms),           Stat = fail(5))

    If (Any(fail > 0)) Call error(1025)

    lsite = 0 ; ltype = 0
    lfrzn = 0 ; lfree = 0

    lexatm = 0 ; ixyz = 0
    list = 0 ; lstfre = 0

    weight = 0.0_wp ; chge = 0.0_wp

! Resize the arrays in allocate_config_arrays_read

    stat = 0

    Call reallocate( mxatms - Size( atmnam ), atmnam, stat( 1) )
    Call reallocate( mxatms - Size( lsi    ), lsi,    stat( 2) )
    Call reallocate( mxatms - Size( lsa    ), lsa,    stat( 3) )
    Call reallocate( mxatms - Size( ltg    ), ltg,    stat( 4) )
    Call reallocate( mxatms - Size( xxx    ), xxx,    stat( 5) )
    Call reallocate( mxatms - Size( yyy    ), yyy,    stat( 6) )
    Call reallocate( mxatms - Size( zzz    ), zzz,    stat( 7) )
    Call reallocate( mxatms - Size( vxx    ), vxx,    stat( 8) )
    Call reallocate( mxatms - Size( vyy    ), vyy,    stat( 9) )
    Call reallocate( mxatms - Size( vzz    ), vzz,    stat(10) )
    Call reallocate( mxatms - Size( fxx    ), fxx,    stat(11) )
    Call reallocate( mxatms - Size( fyy    ), fyy,    stat(12) )
    Call reallocate( mxatms - Size( fzz    ), fzz,    stat(13) )

    If ( Any(stat /= 0 )) Call error(1025)

  End Subroutine allocate_config_arrays

  Subroutine check_config(levcfg,l_str,lpse,keyens,iso,keyfce,keyres,megatm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reporting and checking the configuration data
! in CONFIG for: (i) unsuitable simulation options in CONTROL and
! (ii) connectivity to FIELD; before connecting the crystallographic
! data (positions+) to the topology (sites+), i.e. CONFIG to FIELD
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical, Intent( In    ) :: l_str,lpse
  Integer, Intent( In    ) :: levcfg,keyens,iso,keyfce,keyres,megatm
  Type( comms_type ), Intent( InOut ) :: comm

  Logical, Save     :: newjob = .true.
  Logical           :: safe
  Integer           :: fail,k,l,m, &
                       indatm,totatm,mol_sit,loc_ind
  Real( Kind = wp ) :: rcell(1:9),det

  Integer, Allocatable :: iwrk(:)
  Character( Len = 256 ) :: message

  fail=0
  If (l_str) Then
     Allocate (iwrk(1:mxatms), Stat=fail)
     If (fail > 0) Then
        Write(message,'(a)') 'check_config allocation failure'
        Call error(0,message)
    End If
  End If


  If (newjob) Then
     Write(message,"('configuration file name: ',10x,a)") cfgname
     Call info(message,.true.)
     Write(message,"('selected image convention',6x,i10)") imcon
     Call info(message,.true.)
  End If

! Check things for non-periodic systems

  If (imcon == 0 .or. imcon == 6) Then
     If (keyfce == 2) Then
        Call warning(220,0.0_wp,0.0_wp,0.0_wp)
     Else If (keyfce > 0) Then
        Call warning(30,0.0_wp,0.0_wp,0.0_wp)
     End If

     Call warning(260,0.0_wp,0.0_wp,0.0_wp)

     If (keyens >= 20) Call error(390)
  End If

! Check image conditions for nst ensembles

  If (keyens >= 30) Then
     If (iso == 0) Then
        If (imcon == 1 .or. imcon == 2) Then
           Call warning(110,Real(imcon,wp),3.0_wp,0.0_wp)
           imcon = 3
        End If
     Else ! iso > 0
        If (imcon == 1) Then
           Call warning(110,Real(imcon,wp),3.0_wp,0.0_wp)
           imcon = 2
        End If
     End If
  End If

! Check image condition for pseudo

  If (lpse .and. (imcon == 0 .or. imcon == 6)) Call error(540)

  Call invert(cell,rcell,det)

! Specify molecular dynamics simulation cell

  If (newjob) Then
     Write(message,"('simulation cell vectors')")
     Call Info(message,.true.)
     Write(message,"(3f20.10)") cell(1:3)
     Call Info(message,.true.)
     Write(message,"(3f20.10)") cell(4:6)
     Call Info(message,.true.)
     Write(message,"(3f20.10)") cell(7:9)
     Call Info(message,.true.)
     Write(message,"('system volume     ',2x,1p,g22.12)") det
     Call Info(message,.true.)
  End If

! Check on validity of config file contents

  If (keyres > 0 .and. levcfg < 1) Call error(85)

  If (l_str) iwrk(1:natms) = 0 ! initialise

! Safe flag

  safe=.true.

! Global atom counter

  totatm=0

! Local atom counter

  indatm=1

! Site counter for molecules

  mol_sit=0

! Check atom names and assign atomic characteristics
! Connecting FIELD to CONFIG on local basis

  Do k=1,ntpmls

     Do l=1,nummols(k)

        Do m=1,numsit(k)

! Increase global atom counter

           totatm=totatm+1

! If a local atom has a global index totatm

           If (lsa(indatm) == totatm) Then

! Get the local index. mol_sit+m is the global site

              loc_ind=lsi(indatm)

! Second check: Do particle identities and their order from the topology
! (FIELD) match the one found in the crystallographic data (CONFIG)?
! Check for unidentified atoms in CONFIG by their existence in FIELD

              If (atmnam(loc_ind) /= sitnam(mol_sit+m)) Then
                 Write(message,"( 'unidentified atom label :',a8,': atom number ',i5)") atmnam(loc_ind),loc_ind
                 Call info(message)
                 safe=.false.
              End If

! Assign global site, type, weight, charge & frozen status to localised atoms

              lsite(loc_ind)=mol_sit+m
              ltype(loc_ind)=typsit(mol_sit+m)
              weight(loc_ind)=wgtsit(mol_sit+m)
              chge(loc_ind)=chgsit(mol_sit+m)
              lfrzn(loc_ind)=frzsit(mol_sit+m)
              lfree(loc_ind)=fresit(mol_sit+m)

! Print global indices for a later check on ordering (mixed indexing)

              If (l_str) iwrk(indatm) = totatm ! Populate

! Increase local atom counter

              indatm=indatm+1

           End If

        End Do

     End Do

! Increase site counter per molecule

     mol_sit=mol_sit+numsit(k)

  End Do
  indatm=indatm-1 ! Correct presence number

! Global check on successful connection between CONFIG and FIELD

  Call gcheck(comm,safe)
  If (.not.safe) Call error(25)

! Check CONFIG indices eligibility

  If (l_str) Then
     Call all_inds_present( iwrk, indatm, megatm, safe )
     Call gcheck(comm,safe)
     If (.not.safe) Call error(28)

     Deallocate (iwrk, Stat=fail)
     If (fail > 0) Then
        Write(message,'(a)') 'check_config deallocation failure'
        Call error(0,message)
     End If
  End If

! For subsequent checks

  If (newjob) newjob=.false.

Contains

  Subroutine all_inds_present( ind, n_loc, n, all_present )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to check if the array IND which is distributed
! accross the communicator contains all indices from 1 to N once,
! and only once. N_LOC is the number of indices local to this processor
!
! copyright - daresbury laboratory
! author    - i.j.bush march 2009
! adapted   - i.t.todorov june 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Dimension( : ), Intent( In    ) :: ind
    Integer                , Intent( In    ) :: n_loc
    Integer                , Intent( In    ) :: n
    Logical                , Intent(   Out ) :: all_present

    Integer, Dimension( : ), Allocatable :: all_n_loc
    Integer, Dimension( : ), Allocatable :: loc_start
    Integer, Dimension( : ), Allocatable :: local_ind
    Integer, Dimension( : ), Allocatable :: to_send
    Integer, Dimension( : ), Allocatable :: to_recv
    Integer, Dimension( : ), Allocatable :: displs_send
    Integer, Dimension( : ), Allocatable :: displs_recv
    Integer, Dimension( : ), Allocatable :: reorg_ind

    Integer :: me, nproc
    Integer :: where_send
    Integer :: i,fail

    Logical :: loc_present

    me    = comm%idnode
    nproc = comm%mxnode

    ! Ever the optimist !
    all_present = .true.

    ! If the number on each proc do not add up to N there is no
    ! way things can be correct. will need all the local totals so
    ! gather and sum locally rather than just a global sum as
    ! saves a comm.
    Allocate ( all_n_loc( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    ! No check on mpi error as all is pointless because
    ! of the way mpi error handlers are usually dealt with
!
!    all_n_loc = 0
!    all_n_loc( comm%idnode ) = n_loc
!    Call gsum( all_n_loc( 0:nrpocs - 1 ) )
!
    Call MPI_ALLGATHER(     n_loc, 1, MPI_INTEGER, &
                        all_n_loc, 1, MPI_INTEGER, comm%comm, comm%ierr )
    all_present = ( Sum( all_n_loc ) == n )
    If ( .not. all_present ) Return

    ! Work out the first index on each proc
    Allocate ( loc_start( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    loc_start( 0 ) = 1
    Do i = 1, nproc - 1
       loc_start( i ) = loc_start( i - 1 ) + all_n_loc( i - 1 )
    End Do

    ! Array to work with local indices
    Allocate ( local_ind( 1:n_loc ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    local_ind = ind( 1:n_loc )

    ! Sort the local data
    Call shellsort( n_loc, local_ind )

    ! Work out how much data to send to each processor
    Allocate ( to_send( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    to_send = 0
    where_send = 0
    Do i = 1, n_loc
       Do While ( local_ind( i ) > loc_start( where_send ) + all_n_loc( where_send ) - 1 )
          where_send = where_send + 1
       End Do
       to_send( where_send ) = to_send( where_send ) + 1
    End Do

    ! How much node i sends to me is how much I recv
    Allocate ( to_recv( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    Call MPI_ALLTOALL( to_send, 1, MPI_INTEGER, &
                       to_recv, 1, MPI_INTEGER, &
                       comm%comm, comm%ierr )

    ! Work out the displacements in the sending and receiving arrays
    Allocate ( displs_send( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    Allocate ( displs_recv( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    displs_send( 0 ) = 0
    Do i = 1, nproc - 1
       displs_send( i ) = displs_send( i - 1 ) + to_send( i - 1 )
    End Do

    displs_recv( 0 ) = 0
    Do i = 1, nproc - 1
       displs_recv( i ) = displs_recv( i - 1 ) + to_recv( i - 1 )
    End Do

    ! Put the index on the proc that should own it
    Allocate ( reorg_ind( 1:n_loc ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    Call MPI_ALLTOALLV( local_ind, to_send, displs_send, MPI_INTEGER, &
                        reorg_ind, to_recv, displs_recv, MPI_INTEGER, &
                        comm%comm, comm%ierr )

    ! Sort the reorganized data
    Call shellsort( n_loc, reorg_ind )

    ! Check it starts and ends at the right place for any owned particle
    ! i.e. some domaines my not own any (domain on vacuum)
    If ( n_loc > 0 ) &
       all_present = ( reorg_ind( 1     ) == loc_start( me ) .and. &
                       reorg_ind( n_loc ) == loc_start( me ) + n_loc - 1 )

    If ( all_present ) Then
       ! Check it contains all the numbers in between
       Do i = 2, n_loc
          all_present = ( reorg_ind( i ) - reorg_ind( i - 1 ) == 1 )
          If ( .not. all_present ) Then
             Exit
          End If
       End Do
    End If

    ! Is everybody happy?
    loc_present = all_present
    Call MPI_ALLREDUCE( loc_present, all_present, 1, MPI_LOGICAL, MPI_LAND, comm%comm, comm%ierr )

    Deallocate ( reorg_ind   , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( displs_recv , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( displs_send , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( to_recv     , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( to_send     , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( local_ind   , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( loc_start   , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate ( all_n_loc   , Stat = fail )
    If ( fail /= 0 ) Go To 100

    Return

100 Continue
    If (fail > 0) Then
       Write(message,'(a)') 'all_inds_present allocation/deallocation failure'
       Call error(0,message)
    End If

  End Subroutine all_inds_present

End Subroutine check_config




Subroutine read_config(megatm,levcfg,l_ind,l_str,rcut,dvar,xhi,yhi,zhi,dens0,dens,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading CONFIG and getting the average
! particle density
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: megatm,levcfg
  Logical,           Intent( In    ) :: l_ind,l_str
  Real( Kind = wp ), Intent( In    ) :: rcut,dvar
  Real( Kind = wp ), Intent( InOut ) :: xhi,yhi,zhi
  Real( Kind = wp ), Intent(   Out ) :: dens0,dens
  Type( comms_type), Intent( InOut ) :: comm

  Real( Kind = wp ) :: cut

  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,fname
  Logical                :: safe  = .true.  , &
                            l_his = .false. , &
                            l_xtr = .false. , &
                            fast
  Integer                :: fail(1:4),i,j,idm,max_fail,min_fail, &
                            icell,ncells,                        &
                            indatm,nattot,totatm,                &
                            ipx,ipy,ipz,nlx,nly,nlz,             &
                            ix,iy,iz,jx,jy,jz
  Real( Kind = wp )      :: celprp(1:10),rcell(1:9),celh(1:9),det, &
                            volm,vcell,                            &
                            sxx,syy,szz,xdc,ydc,zdc,               &
                            pda_max,pda_min,pda_ave,               &
                            pda_dom_max,pda_dom_min

! Some parameters and variables needed by io interfaces

  Integer                           :: recsz = 73 ! default record size
  Integer                           :: fh, io_read
  Integer( Kind = MPI_OFFSET_KIND ) :: top_skip

  Real( Kind = wp ),    Dimension( : ), Allocatable :: pda

  Character( Len = 8 ), Dimension( : ), Allocatable :: chbuf
  Integer,              Dimension( : ), Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ), Allocatable :: axx,ayy,azz, &
                                                       bxx,byy,bzz, &
                                                       cxx,cyy,czz

  Character( Len = 256) :: message
! image conditions not compliant with DD and link-cell

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Real space cutoff shortened by 50% but not < 1 Angstrom
!(or ==rcut_def in scan_control)

  cut=Max(0.5_wp*rcut,1.0_wp)+1.0e-6_wp

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)
  volm = celprp(10)

! Calculate the number of link-cells per domain in every direction

  nlx=Int(celprp(7)/(cut*nprx_r))
  nly=Int(celprp(8)/(cut*npry_r))
  nlz=Int(celprp(9)/(cut*nprz_r))

  ncells=nlx*nly*nlz

! Check for link cell algorithm violations

  If (ncells == 0) Call error(307)

! Amend volume of density cell if cluster, slab or bulk slab
! cell dimensional properties overwritten but not needed anyway

  If (imcon == 0 .or. imcon == 6 .or. imc_n == 6) Then
     celh=cell

     If (imcon == 0) Then
        celh(1) = Max(1.0_wp,xhi)
        celh(5) = Max(1.0_wp,yhi)
        celh(9) = Max(1.0_wp,zhi)
     Else If (imcon == 6) Then
        celh(9) = Max(1.0_wp,zhi)
     End If

     Call dcell(celh,celprp)
     volm = celprp(10)
  End If

  vcell = volm / (Real(ncells,wp) * Real(comm%mxnode,wp))

! Approximate density and mxatms

  dens = Real(megatm,wp) / volm
  mxatms = Max(1 , Nint( (dvar**1.7_wp) * dens*vcell * Real((nlx+3)*(nly+3)*(nlz+3),wp)))

! Allocate necessary arrays to read CONFIG

  Call allocate_config_arrays_read(mxatms)

! Get type of I/O for reading

  Call io_get_parameters( user_method_read = io_read )

! Define filename ASCII or netCDF

  If (io_read /= IO_READ_NETCDF) Then
     fname=Trim(config)
  Else
     fname=Trim(config) // '.nc'
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

     fast = .true.
     If (comm%idnode == 0) Then

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
        Call gsync(comm)
        Call gcheck(comm,safe,"enforce")
        Call gcheck(comm,fast,"enforce")

     If (.not.safe) Go To 50

! Close CONFIG

     If (comm%idnode == 0) Close(Unit=nconf)

  End If

  fail = 0

! If MASTER read

  If (io_read == IO_READ_MASTER) Then

     Call invert(cell,rcell,det)

! Open CONFIG and skip the header

     If (comm%idnode == 0) Then
        Open(Unit=nconf, File=fname)

        Read(Unit=nconf, Fmt=*)    ! CONFIG file header (TITLE record)
        Read(Unit=nconf, Fmt=*)    ! configuration level and image condition

        If (imcon /= 0) Then
           Read(Unit=nconf, Fmt=*) ! cell vectors (not defined for imcon=0) but cell
           Read(Unit=nconf, Fmt=*) ! is modified in set_bounds for imcon 0 and 6!!!
           Read(Unit=nconf, Fmt=*)
        End If
     End If

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'read_config allocation failure'
        Call error(0,message)
     End If

! Initialise domain localised atom counter (configuration)
! and dispatched atom counter

     natms =0
     indatm=0
     Do nattot=1,megatm
        indatm=indatm+1

! Initialise transmission arrays

        chbuf(indatm)=' '
        iwrk(indatm)=0

        axx(indatm)=0.0_wp
        ayy(indatm)=0.0_wp
        azz(indatm)=0.0_wp

        If (levcfg > 0) Then
           bxx(indatm)=0.0_wp
           byy(indatm)=0.0_wp
           bzz(indatm)=0.0_wp

           If (levcfg > 1) Then
              cxx(indatm)=0.0_wp
              cyy(indatm)=0.0_wp
              czz(indatm)=0.0_wp
           End If
        End If

! Read in transmission arrays

        If (comm%idnode == 0 .and. safe) Then
           record=' '
           Read(Unit=nconf, Fmt='(a)', End=30) record
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

           Read(Unit=nconf, Fmt=*, End=30) axx(indatm),ayy(indatm),azz(indatm)

           If (levcfg > 0) Then
              Read(Unit=nconf, Fmt=*, End=30) bxx(indatm),byy(indatm),bzz(indatm)
              If (levcfg > 1) Read(Unit=nconf, Fmt=*, End=30) cxx(indatm),cyy(indatm),czz(indatm)
           End If
           Go To 40

30         Continue
           safe=.false. ! catch error

40         Continue
        End If

! Circulate configuration data to all nodes when transmission arrays are filled up

        If (indatm == mxatms .or. nattot == megatm) Then

! Check if batch was read fine

           Call gcheck(comm,safe)
           If (.not.safe) Go To 50

! Ensure all atoms are in prescribed simulation cell (DD bound) and broadcast them
!
!           Call pbcshift(imcon,cell,indatm,axx,ayy,azz)

              Call gbcast(comm,chbuf,0)
              Call gbcast(comm,iwrk,0)

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

! Assign atoms to correct domains

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
              If      (idm < 0 .or. idm > (comm%mxnode-1)) Then
                 Call error(513)
               Else If (idm == comm%idnode)                 Then
                 natms=natms+1

                 If (natms < mxatms) Then
                    atmnam(natms)=chbuf(i)
                    ltg(natms)=iwrk(i)

                    xxx(natms)=axx(i)
                    yyy(natms)=ayy(i)
                    zzz(natms)=azz(i)

                    If (levcfg > 0) Then
                       vxx(natms)=bxx(i)
                       vyy(natms)=byy(i)
                       vzz(natms)=bzz(i)
                    Else
                       vxx(natms)=0.0_wp
                       vyy(natms)=0.0_wp
                       vzz(natms)=0.0_wp
                    End If

                    If (levcfg > 1) Then
                       fxx(natms)=cxx(i)
                       fyy(natms)=cyy(i)
                       fzz(natms)=czz(i)
                    Else
                       fxx(natms)=0.0_wp
                       fyy(natms)=0.0_wp
                       fzz(natms)=0.0_wp
                    End If
                 Else
                    safe=.false.
                 End If
              End If
           End Do

! Check if all is dispatched fine

           max_fail=natms
           min_fail=natms
              Call gcheck(comm,safe)
              Call gmax(comm,max_fail)
              Call gmin(comm,min_fail)

           If (.not.safe) Then
              If (comm%idnode == 0) Then
  Write(nrite,'(/,1x,a,i0)')  '*** warning - next error due to maximum number of atoms per domain set to : ', mxatms
  Write(nrite,'(1x,2(a,i0))') '***           but maximum & minumum numbers of atoms per domain asked for : ', &
       max_fail, ' & ', min_fail
  Write(nrite,'(1x,a,i0)')    '***           estimated densvar value for passing this stage safely is : ', &
       Ceiling((dvar*(Real(max_fail,wp)/Real(mxatms,wp))**(1.0_wp/1.7_wp)-1.0_wp)*100.0_wp)
              End If
              Call error(45)
           End If

! Nullify dispatch counter

           indatm=0

        End If
     End Do

! Close CONFIG

     If (comm%idnode == 0) Close(Unit=nconf)
     Call gsync(comm)

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'read_config deallocation failure'
        Call error(0,message)
     End If

! If PROPER read

  Else

! Open CONFIG

     If (fast) Then
        Call io_set_parameters( user_comm = comm%comm )
        Call io_init( recsz )
        Call io_open( io_read, comm%comm, fname, MPI_MODE_RDONLY, fh )
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

     Call read_config_parallel                  &
           (levcfg, dvar, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi,comm)

! Close CONFIG

     If (fast) Then
        Call io_close( fh )
        Call io_finalize
     Else
        Close(Unit=nconf)
     End If

  End If

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  nlast=natms

! Check if the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?

  totatm=natms
  Call gsum(comm,totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (configuration)

  Do i=1,natms
     lsi(i)=i
     lsa(i)=ltg(i)
  End Do
  Call shellsort2(natms,lsi,lsa)

  If (io_read /= IO_READ_MASTER) Then

! This section is not strictly necessary.  However, the new read in method
! means the atoms are not necessarily in the same order in memory as the
! older, slower, method would put them.  This bit makes sure that the order
! is so that 'ltg' is strictly monotonically increasing.  This captures the
! common case where CONFIG has the 'ltg' values all in order (or not
! specified), but there is no easy way for the general case of arbitrary
! ordering of the 'ltg' values in CONFIG.  Of course, this makes no
! difference to the science and to restarts.  However, for initial runs it
! means the initial velocities will not be the same as the old method for
! the arbitrary ordering case.

     atmnam( 1:natms ) = atmnam( lsi( 1:natms ) )
     ltg( 1:natms ) = ltg( lsi( 1:natms ) )

     xxx( 1:natms ) = xxx( lsi( 1:natms ) )
     yyy( 1:natms ) = yyy( lsi( 1:natms ) )
     zzz( 1:natms ) = zzz( lsi( 1:natms ) )

     If (levcfg > 0) Then
        vxx( 1:natms ) = vxx( lsi( 1:natms ) )
        vyy( 1:natms ) = vyy( lsi( 1:natms ) )
        vzz( 1:natms ) = vzz( lsi( 1:natms ) )

        If (levcfg > 1) Then
           fxx( 1:natms ) = fxx( lsi( 1:natms ) )
           fyy( 1:natms ) = fyy( lsi( 1:natms ) )
           fzz( 1:natms ) = fzz( lsi( 1:natms ) )
        End If
     End If
     Do i=1,natms
        lsi(i)=i
        lsa(i)=ltg(i)
     End Do

  End If

! READ CONFIG END

! PARTICLE DENSITY START
! Allocate and initialise particle density array

  Allocate (pda(1:ncells), Stat=fail(1))
  If (fail(1) > 0) Then
     Write(message,'(a)') 'read_config allocation failure'
     Call error(0,message)
  End If
  pda=0.0_wp

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! Shifts from global to local link-cell space:
! (0,0,0) left-most link-cell on the domain (halo)
! (nlx+2*nlp-1,nly+2*nlp-1,nly+2*nlp-1) right-most
! link-cell on the domain (halo)

  jx=1-nlx*idx
  jy=1-nly*idy
  jz=1-nlz*idz

! Get the inverse cell matrix

  Call invert(cell,rcell,celprp(10))

  Do i=1,natms
     sxx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     syy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     szz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! Get cell coordinates accordingly

     ix = Int(xdc*(sxx+0.5_wp)) + jx
     iy = Int(ydc*(syy+0.5_wp)) + jy
     iz = Int(zdc*(szz+0.5_wp)) + jz

! Put all particles in bounded link-cell space: lower and upper
! bounds as 1 <= i_coordinate <= nl_coordinate

     ix = Max( Min( ix , nlx) , 1)
     iy = Max( Min( iy , nly) , 1)
     iz = Max( Min( iz , nlz) , 1)

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/

     icell=1+(ix-1)+nlx*((iy-1)+nly*(iz-1))

     pda(icell)=pda(icell)+1.0_wp
  End Do

  pda_max=  0.0_wp
  pda_min=100.0_wp
  pda_ave=  0.0_wp
  Do icell=1,ncells
     pda_max=Max(pda(icell),pda_max)
     pda_min=Min(pda(icell),pda_min)
     pda_ave=pda_ave+pda(icell)
  End Do
  pda_ave=pda_ave/Real(ncells,wp)

  pda_dom_max=pda_ave
  pda_dom_min=pda_ave
     Call gmax(comm,pda_dom_max)
     Call gmin(comm,pda_dom_min)

     Call gmax(comm,pda_max)
     Call gmin(comm,pda_min)

     Call gsum(comm,pda_ave)
     pda_ave=pda_ave/Real(comm%mxnode,wp)

! Approximation for maximum global density by
! the inter-domain imbalance of domain density

  dens0=pda_max/vcell
  If (comm%mxnode > 1) Then
     If (Nint(pda_dom_min) == 0) Then
        dens = dens0 ! domain(s) matched on vacuum (take no risk)
     Else If (1.15_wp*pda_dom_min > pda_dom_max) Then
        dens = pda_ave/vcell
     Else If (1.25_wp*pda_dom_min > pda_dom_max) Then
        dens = pda_dom_max/vcell
     Else
        dens = dens0 ! too big an imbalance (take no risk)
     End If
  Else
     dens = 1.25_wp*pda_ave/vcell ! allow 25% imbalance
  End If

  Deallocate (pda, Stat=fail(1))
  If (fail(1) > 0) Then
     Write(message,'(a)') 'read_config deallocation failure'
     Call error(0,message)
  End If

! PARTICLE DENSITY END

  Return

! error exit for CONFIG file read

50 Continue
  If (comm%idnode == 0) Close(Unit=nconf)
  Call error(55)

End Subroutine read_config

Subroutine read_config_parallel                 &
           (levcfg, dvar, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the CONFIG data file in parallel
!
! copyright - daresbury laboratory
! author    - i.j.bush & i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Logical,                           Intent( In    ) :: l_ind,l_str,l_his,fast,l_xtr
  Integer,                           Intent( In    ) :: levcfg,megatm,fh
  Integer( Kind = MPI_OFFSET_KIND ), Intent( In    ) :: top_skip
  Real( Kind = wp ),                 Intent( In    ) :: dvar
  Real( Kind = wp ),                 Intent(   Out ) :: xhi,yhi,zhi
  Type( comms_type ),                Intent( InOut ) :: comm

  Logical                :: safe,do_read
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,forma
  Integer                :: fail(1:8),i,j,k,max_fail,min_fail, &
                            idm,ipx,ipy,ipz,indatm,            &
                            n_read_procs_use,per_read_proc,    &
                            my_read_proc_num,ats_per_proc,     &
                            recs_per_at,recs_per_proc,         &
                            wp_vals_per_at,n_loc,              &
                            to_read,which_read_proc,this_base_proc
  Integer( Kind = li )   :: n_sk,n_ii,n_jj
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz

! Some parameters and variables needed by io interfaces

  Integer                           :: io_read
  Integer                           :: recsz, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io, n_skip
  Integer                           :: this_rec_buff, recs_to_read
  Integer                           :: n_ats_in_file

! netCDF

  Integer :: frame, start(1:3), count(1:3)

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx_read,ayy_read,azz_read, &
                                                          bxx_read,byy_read,bzz_read, &
                                                          cxx_read,cyy_read,czz_read

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf_read,chbuf_scat
  Integer,              Dimension( : ),    Allocatable :: iwrk_read,iwrk_scat

  Integer,              Dimension( : ),    Allocatable :: first_at,orig_first_at
  Integer,              Dimension( : ),    Allocatable :: n_held,where_buff
  Integer,              Dimension( : ),    Allocatable :: owner_read

  Real( Kind = wp ),    Dimension( :, : ), Allocatable :: scatter_buffer_read
  Real( Kind = wp ),    Dimension( :, : ), Allocatable :: scatter_buffer

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: rec_buff
  Integer :: ierr

  Character( Len =  256 )  ::  message
 

! Get reading method, total number of I/O heads and buffer size

  Call io_get_parameters( user_method_read      = io_read          )
  Call io_get_parameters( user_n_io_procs_read  = n_read_procs_use )
  Call io_get_parameters( user_buffer_size_read = batsz            )

  fail = 0 ! fail initialisation

  If (levcfg /= 3) Then
     wp_vals_per_at = 3 * (1+levcfg) ! Scatter buffer sizes
     recs_per_at    = 2 + levcfg     ! Scatter buffer sizes
  Else
     wp_vals_per_at = 3 ! Scatter buffer sizes
     recs_per_at    = 1 ! Scatter buffer sizes
  End If

! Note: make 'first_at' and 'orig_first_at' 1 element bigger than strictly
! required to make checking at the end of reading much easier and clearer

  Allocate (first_at(0:n_read_procs_use),orig_first_at(0:n_read_procs_use), Stat=fail(1))
  Allocate (chbuf(1:batsz),iwrk(1:batsz),                                   Stat=fail(2))
  Allocate (scatter_buffer(1:wp_vals_per_at,1:batsz),                       Stat=fail(3))
  If (Any(fail(1:3) > 0)) Then
     Write(message,'(a)') 'read_config_parallel allocation failure 1'
     Call error(0,message)
  End If

! define basic quantities for the parallel ASCII reading

  per_read_proc = comm%mxnode / n_read_procs_use
  do_read = (Mod( comm%idnode, per_read_proc ) == 0 .and. comm%idnode < per_read_proc * n_read_procs_use)
  my_read_proc_num = comm%idnode / per_read_proc

! Note 'first_at' and 'orig_first_at' have one more element
! in the array than strictly required - makes it easier to
! check that reading by the last I/O processor has finished

  ats_per_proc = megatm / n_read_procs_use
  Do i=0,n_read_procs_use
     first_at(i) = i*ats_per_proc + Min(i,megatm-ats_per_proc*n_read_procs_use)
  End Do
  orig_first_at = first_at
  ats_per_proc = Max(1,ats_per_proc) ! Fix it if 0
  recs_per_proc = ats_per_proc * recs_per_at

! Catch the case where the first atom belonging to
! a read processor does not actually exists - i.e.
! I/O procs count > megatm, and limit reading by do_read

  If (my_read_proc_num < n_read_procs_use) &
     do_read = (do_read .and. first_at(my_read_proc_num) < megatm)

! Skip to the point of reading

  If (do_read) Then

     n_skip = Int(recs_per_at,MPI_OFFSET_KIND) * Int(first_at(my_read_proc_num),MPI_OFFSET_KIND) + &
              top_skip-Int(1,MPI_OFFSET_KIND)
     If (.not.fast) Then
        n_sk=Int(n_skip,li)
        n_jj=73*batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
           Do n_ii=1_li,n_sk/n_jj
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_jj
              Read(Unit=nconf, Fmt=forma, End=100)
           End Do
           n_ii=Mod(n_sk,n_jj)-1_li
           If (n_ii > 0_li) Then
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_ii
              Read(Unit=nconf, Fmt=forma, End=100)
           End If
        Else
           forma=' '
           Write(forma,'( "(", i0, "/)" )') n_sk
           Read(Unit=nconf, Fmt=forma, End=100)
        End If

        recsz=200
        forma=' '
        Write(forma,'( "(", i0, "a1)" )') recsz
     Else
        rec_mpi_io = n_skip + Int(1,MPI_OFFSET_KIND)
        recsz=73
        If (levcfg == 3) recsz = 35
     End If

! Allocate record buffer, reading buffers, scatter buffers and indexing arrays

     If (io_read /= IO_READ_NETCDF) Then
        Allocate (rec_buff(1:recsz,1:batsz),                                  Stat=fail(1))
     Else
        Allocate (rec_buff(1:Len( chbuf_read ),1:batsz),                      Stat=fail(1))
     End If
     Allocate (chbuf_read(1:batsz),iwrk_read(1:batsz),                        Stat=fail(2))
     Allocate (axx_read(1:batsz),ayy_read(1:batsz),azz_read(1:batsz),         Stat=fail(3))
     Allocate (bxx_read(1:batsz),byy_read(1:batsz),bzz_read(1:batsz),         Stat=fail(4))
     Allocate (cxx_read(1:batsz),cyy_read(1:batsz),czz_read(1:batsz),         Stat=fail(5))
     Allocate (scatter_buffer_read(1:wp_vals_per_at,1:batsz),                 Stat=fail(6))
     Allocate (chbuf_scat(1:batsz),iwrk_scat(1:batsz),                        Stat=fail(7))
     Allocate (n_held(0:comm%mxnode-1),where_buff(0:comm%mxnode-1),owner_read(1:batsz), Stat=fail(8))
     If (Any(fail(1:8) > 0)) Then
        Write(message,'(a)') 'read_config_parallel allocation failure 2'
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
        Write(message,'(a)') 'read_config_parallel allocation failure 3'
        Call error(0,message)
     End If

  End If

! Initialise extreme box dimensions

  xhi = 0.0_wp
  yhi = 0.0_wp
  zhi = 0.0_wp

  If (.not.l_xtr) Call invert(cell,rcell,det)

! Initialise domain localised atom counter (configuration),
! dispatched atom counter and safe dispatch flag

  natms =0
  indatm=0
  safe  =.true.

  Do k=1,megatm

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
                    Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
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

              If (levcfg == 3) Read(record, Fmt=*, End=100) axx_read(i),ayy_read(i),azz_read(i)

              If (this_rec_buff == recs_to_read) Then
                 this_rec_buff = 0
                 recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 1 )
                 If (.not.fast) Then
                    Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                 End If
              End If

              If (levcfg /= 3) Then

! Positions

                 this_rec_buff = this_rec_buff + 1
                 Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                    record( j:j ) = rec_buff( j, this_rec_buff )
                 End Do
                 Read(record, Fmt=*, End=100) axx_read(i),ayy_read(i),azz_read(i)
                 If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
                    If (levcfg > 0) Then
                       recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 2 )
                       If (.not.fast) Then
                          Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                       Else
                          Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                          rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                       End If
                    End If
                 End If

! Velocities

                 If (levcfg > 0) Then
                    this_rec_buff = this_rec_buff + 1
                    Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                       record( j:j ) = rec_buff( j, this_rec_buff )
                    End Do
                    Read(record, Fmt=*, End=100) bxx_read(i),byy_read(i),bzz_read(i)
                    If (this_rec_buff == recs_to_read) Then
                       this_rec_buff = 0
                       If (levcfg > 1) Then
                          recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 3 )
                          If (.not.fast) Then
                             Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                          Else
                             Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                             rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                          End If
                       End If
                    End If

! Forces

                    If (levcfg > 1) Then
                       this_rec_buff = this_rec_buff + 1
                       Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                          record( j:j ) = rec_buff( j, this_rec_buff )
                       End Do
                       Read(record, Fmt=*, End=100) cxx_read(i),cyy_read(i),czz_read(i)
                       If (this_rec_buff == recs_to_read) Then
                          this_rec_buff = 0
                       End If
                    End If
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

              Select Case( levcfg )
              Case( 0, 3 )
                 Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
              Case( 1 )
                 Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
                 Call get_var( 'velocities' , fh, start, count, bxx_read, byy_read, bzz_read )
              Case( 2 )
                 Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
                 Call get_var( 'velocities' , fh, start, count, bxx_read, byy_read, bzz_read )
                 Call get_var( 'forces'     , fh, start, count, cxx_read, cyy_read, czz_read )
              End Select
           End If

        End If No_netCDF

        If (.not.l_xtr) Then

! Ensure all atoms are in prescribed simulation cell (DD bound)
!
           n_held=0
           Do i=1,to_read
              sxx=rcell(1)*axx_read(i)+rcell(4)*ayy_read(i)+rcell(7)*azz_read(i)
              syy=rcell(2)*axx_read(i)+rcell(5)*ayy_read(i)+rcell(8)*azz_read(i)
              szz=rcell(3)*axx_read(i)+rcell(6)*ayy_read(i)+rcell(9)*azz_read(i)

! sxx,syy,szz are in [-0.5,0.5) interval as values as 0.4(9) may pose a problem

              sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
              syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
              szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

! fold back coordinates

              axx_read(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
              ayy_read(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
              azz_read(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz

! assign domain coordinates (call for errors)

              ipx=Int((sxx+0.5_wp)*nprx_r)
              ipy=Int((syy+0.5_wp)*npry_r)
              ipz=Int((szz+0.5_wp)*nprz_r)

              idm=ipx+nprx*(ipy+npry*ipz)
              If (idm < 0 .or. idm > (comm%mxnode-1)) Call error(513)
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

              If (levcfg /= 3) Then
                 If (levcfg > 0) Then
                    scatter_buffer_read(4,idm) = bxx_read(i)
                    scatter_buffer_read(5,idm) = byy_read(i)
                    scatter_buffer_read(6,idm) = bzz_read(i)

                    If (levcfg > 1) Then
                       scatter_buffer_read(7,idm) = cxx_read(i)
                       scatter_buffer_read(8,idm) = cyy_read(i)
                       scatter_buffer_read(9,idm) = czz_read(i)
                    End If
                 End If
              End If
           End Do

! If only detecting box dimensions for imcon == 0 or 6 or imc_n == 6

        Else

! Get extremes

           xhi = Max( xhi, Maxval( Abs( axx_read( 1:to_read ) ) ) )
           yhi = Max( yhi, Maxval( Abs( ayy_read( 1:to_read ) ) ) )
           zhi = Max( zhi, Maxval( Abs( azz_read( 1:to_read ) ) ) )

        End If
     End If Readers_only

! Increase buffer counter and update first_at for
! the readers that have something left to read

     indatm = indatm+1
     If (do_read) Then
        If (first_at(my_read_proc_num) < first_at(my_read_proc_num+1)) &
             first_at(my_read_proc_num) = first_at(my_read_proc_num)+1
     End If

! Circulate configuration data to all nodes when transmission arrays are filled up
! Check against megatm since at low processors counts (i.e. 1) batsz can be > megatm

     Reorganize_buffer: If (indatm == batsz .or. (indatm > 0 .and. k == megatm)) Then

        Extent_2: If (.not.l_xtr) Then

           Do which_read_proc = 0 , n_read_procs_use-1
              If (orig_first_at(which_read_proc) >= megatm) Exit ! for non-reading readers

              this_base_proc = which_read_proc * per_read_proc
              If (comm%idnode == this_base_proc) Then
                 where_buff(0) = 0
                 Do i=1,comm%mxnode-1
                    where_buff(i) = where_buff(i-1) + n_held(i-1)
                 End Do
              End If

              Call MPI_SCATTER( n_held, 1, MPI_INTEGER, n_loc, 1, MPI_INTEGER, this_base_proc, &
                   comm%comm, comm%ierr )

              Call MPI_SCATTERV( chbuf_scat, 8 * n_held, 8 * where_buff, MPI_CHARACTER, &
                                 chbuf     , 8 * n_loc ,                 MPI_CHARACTER, &
                                 this_base_proc, comm%comm, comm%ierr )

              Call MPI_SCATTERV( iwrk_scat ,     n_held,     where_buff, MPI_INTEGER, &
                                 iwrk      ,     n_loc ,                 MPI_INTEGER, &
                                 this_base_proc, comm%comm, comm%ierr )

              Call MPI_SCATTERV( scatter_buffer_read, wp_vals_per_at * n_held, wp_vals_per_at * where_buff, wp_mpi, &
                                 scatter_buffer     , wp_vals_per_at * n_loc ,                              wp_mpi, &
                                 this_base_proc, comm%comm, comm%ierr )

! Assign atoms to correct domains

              Do i=1,n_loc
                 natms=natms+1

! Check safety by the upper bound of: atmnam,ltg,xxx,yyy,zzz &
! possibly vxx,vyy,vzz & possibly fxx,fyy,fzz as guided by xxx

                 If (natms <= mxatms) Then
                    atmnam(natms)=chbuf(i)
                    ltg(natms)=iwrk(i)

                    xxx(natms)=scatter_buffer(1,i)
                    yyy(natms)=scatter_buffer(2,i)
                    zzz(natms)=scatter_buffer(3,i)

                    If (levcfg /=3 ) Then
                       If (levcfg > 0) Then
                          vxx(natms)=scatter_buffer(4,i)
                          vyy(natms)=scatter_buffer(5,i)
                          vzz(natms)=scatter_buffer(6,i)
                       Else
                          vxx(natms)=0.0_wp
                          vyy(natms)=0.0_wp
                          vzz(natms)=0.0_wp
                       End If

                       If (levcfg > 1) Then
                          fxx(natms)=scatter_buffer(7,i)
                          fyy(natms)=scatter_buffer(8,i)
                          fzz(natms)=scatter_buffer(9,i)
                       Else
                          fxx(natms)=0.0_wp
                          fyy(natms)=0.0_wp
                          fzz(natms)=0.0_wp
                       End If
                    End If
                  Else
                    safe=.false.
                 End If
              End Do
           End Do

! Check if all is dispatched fine

           max_fail=natms
           min_fail=natms
              Call gcheck(comm,safe)
              Call gmax(comm,max_fail)
              Call gmin(comm,min_fail)

           If (.not.safe) Then
              If (comm%idnode == 0) Then
  Write(nrite,'(/,1x,a,i0)')  '*** warning - next error due to maximum number of atoms per domain set to : ', mxatms
  Write(nrite,'(1x,2(a,i0))') '***           but maximum & minumum numbers of atoms per domain asked for : ', &
       max_fail, ' & ', min_fail
  Write(nrite,'(1x,a,i0)')    '***           estimated densvar value for passing this stage safely is : ', &
       Ceiling((dvar*(Real(max_fail,wp)/Real(mxatms,wp))**(1.0_wp/1.7_wp)-1.0_wp)*100.0_wp)
              End If
              Call error(45)
           End If

        End If Extent_2

! Nullify dispatch counter

        indatm=0

     End If Reorganize_buffer

  End Do

! If only detecting box dimensions for imcon == 0 or 6 or imc_n == 6

  If (l_xtr) Then
     Call gmax(comm,xhi)
     Call gmax(comm,yhi)
     Call gmax(comm,zhi)
  End If

  If (l_his) Then

! Skip to the EoFrame of HISTORY when not fast

     If (do_read .and. (.not.fast)) Then
        n_skip = Int(recs_per_at,MPI_OFFSET_KIND) * Int(megatm-first_at(my_read_proc_num),MPI_OFFSET_KIND)

        n_sk=Int(n_skip,li)
        n_jj=73*batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
           Do n_ii=1_li,n_sk/n_jj
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_jj
              Read(Unit=nconf, Fmt=forma, End=100)
           End Do
           n_ii=Mod(Int(n_skip,li),n_jj)
           If (n_ii > 0_li) Then
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_ii
              Read(Unit=nconf, Fmt=forma, End=100)
           End If
        Else
           forma=' '
           Write(forma,'( "(", i0, "/)" )') n_sk
           Read(Unit=nconf, Fmt=forma, End=100)
        End If
     End If

  Else

     If (do_read) Then

        If ( io_read /= IO_READ_NETCDF) Then

! The last reader to check for EoFile in CONFIG
! and if none is hit to call error to abort

           If (first_at(my_read_proc_num) == megatm) Then
              recs_to_read = 1
              If (.not.fast) Then
                 Read(Unit=nconf, Fmt=forma, Iostat=ierr) rec_buff( :, 1:recs_to_read )
              Else
                 Call io_read_batch( fh, rec_mpi_io, 1, rec_buff, ierr )
              End If
              safe = (ierr /= 0)
           End If

        Else

! As netCDF files have no real concept of line numbers,
! instead check the arrays are the correct size

           Call io_nc_get_dim( 'atom', fh, n_ats_in_file )
           safe = n_ats_in_file == megatm

        End If

     End If
     Call gcheck(comm,safe)
     If (.not.safe) Call error(58)

  End If

  If (do_read) Then
     Deallocate (rec_buff,                   Stat=fail(1))
     Deallocate (chbuf_read,iwrk_read,       Stat=fail(2))
     Deallocate (axx_read,ayy_read,azz_read, Stat=fail(3))
     Deallocate (bxx_read,byy_read,bzz_read, Stat=fail(4))
     Deallocate (cxx_read,cyy_read,czz_read, Stat=fail(5))
     Deallocate (owner_read,                 Stat=fail(6))
     If (Any(fail(1:6) > 0)) Then
        Write(message,'(a)') 'read_config_parallel deallocation failure 2'
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
     Write(message,'(a)') 'read_config_parallel deallocation failure 1'
     Call error(0,message)
  End If

  Return

! error exit for CONFIG file read

100 Continue
  Call error(55)

Contains

  Subroutine get_var( what, fh, start, count, x, y, z )

    Character( Len = * )             , Intent( In    ) :: what
    Integer                          , Intent( In    ) :: fh
    Integer   ,        Dimension( : ), Intent( In    ) :: start
    Integer   ,        Dimension( : ), Intent( In    ) :: count
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: x
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: y
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: z

    Real( Kind = wp ), Dimension( :, : ), Allocatable :: buff

    Integer :: to_read
    Integer :: fail
    Integer :: i

    to_read = count( 2 )

    Allocate (buff( 1:3, 1:to_read ), Stat=fail)
    If (fail /= 0) Then
       Write( nrite, '(/,1x,a,i0)') 'read_config_parallel allocation failure 4, node: ', comm%idnode
       Call error( 0 )
    End If

    Call io_nc_get_var( what, fh, buff, start, count )

    Do i = 1, to_read
       x( i ) = buff( 1, i )
       y( i ) = buff( 2, i )
       z( i ) = buff( 3, i )
    End Do

    Deallocate (buff, Stat=fail)
    If (fail /= 0) Then
       Write( nrite, '(/,1x,a,i0)') 'read_config_parallel allocation failure 4, node: ', comm%idnode
       Call error( 0 )
    End If

  End Subroutine get_var

End Subroutine read_config_parallel

Subroutine scan_config(megatm,imc_n,dvar,cfgname,levcfg,imcon,cell,xhi,yhi,zhi,comm)

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

  Integer,              Intent( In    ) :: megatm,imc_n
  Real( Kind = wp ),    Intent( In    ) :: dvar
  Character( Len = * ), Intent(   Out ) :: cfgname
  Integer,              Intent(   Out ) :: levcfg,imcon
  Real( Kind = wp ),    Intent(   Out ) :: cell(1:9),xhi,yhi,zhi
  Type( comms_type ),   Intent( InOut ) :: comm

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

! Some parameters and variables needed by io interfaces

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

  If (comm%idnode == 0) Inquire(File=fname, Exist=safe)
  Call gcheck(comm,safe)
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
     If (comm%idnode == 0) Then

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

        Call gsync(comm)
        Call gcheck(comm,safe,"enforce")
        Call gcheck(comm,fast,"enforce")

     If (.not.safe) Go To 50

! Close CONFIG

     If (comm%idnode == 0) Close(Unit=nconf)

  End If

!!! SCAN HEADER

  If (io_read /= IO_READ_NETCDF) Then ! ASCII read

! Open CONFIG

     If (comm%idnode == 0) Open(Unit=nconf, File=fname)

! Read TITLE record (file header)

     Call get_line(safe,nconf,record,comm)
     If (.not.safe) Go To 50

     Call strip_blanks(record)
     cfgname=record

! Read configuration level and image condition

     Call get_line(safe,nconf,record,comm)
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
        Call get_line(safe,nconf,record,comm)
        If (.not.safe) Go To 50
        Call get_word(record,word)
        cell(1)=word_2_real(word)
        Call get_word(record,word)
        cell(2)=word_2_real(word)
        Call get_word(record,word)
        cell(3)=word_2_real(word)

        Call get_line(safe,nconf,record,comm)
        If (.not.safe) Go To 50
        Call get_word(record,word)
        cell(4)=word_2_real(word)
        Call get_word(record,word)
        cell(5)=word_2_real(word)
        Call get_word(record,word)
        cell(6)=word_2_real(word)

        Call get_line(safe,nconf,record,comm)
        If (.not.safe) Go To 50
        Call get_word(record,word)
        cell(7)=word_2_real(word)
        Call get_word(record,word)
        cell(8)=word_2_real(word)
        Call get_word(record,word)
        cell(9)=word_2_real(word)
     End If

! Close CONFIG

     If (comm%idnode == 0) Close(Unit=nconf)
     Call gsync(comm)

  Else ! netCDF read

! Open CONFIG

     Call io_set_parameters( user_comm = comm%comm )
     Call io_open( io_read, comm%comm, fname, MPI_MODE_RDONLY, fh )

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

        If (comm%idnode == 0) Then

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

        Call gsync(comm)
        Call gcheck(comm,safe,"enforce")

        If (.not.safe) Go To 50

! Close CONFIG

        If (comm%idnode == 0) Close(Unit=nconf)


           buffer(1)=xhi
           buffer(2)=yhi
           buffer(3)=zhi
           buffer(4)=Real(totatm,wp)
           Call gsum(comm,buffer)
           xhi=buffer(1)
           yhi=buffer(2)
           zhi=buffer(3)
           totatm=Nint(buffer(4))
           If (totatm /= megatm) Call error(58)


! If PROPER read

     Else

! Open CONFIG

        If (fast) Then
           Call io_set_parameters( user_comm = comm%comm )
           Call io_init( recsz )
           Call io_open( io_read, comm%comm, fname, MPI_MODE_RDONLY, fh )
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
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi,comm)

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
  If (comm%idnode == 0) Close(Unit=nconf)
  Call error(55)

End Subroutine scan_config

Subroutine scale_config(megatm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for rescaling the crystallographic information
! from CONFIG to new lattice parameters and saving it in CFGSCL
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer, Intent( In    ) :: megatm

  Character ( Len = 6 ) :: name
  Integer               :: i,nstep
  Real( Kind = wp )     :: rcell(1:9),det,uuu,vvv,www,tstep,time
  Type( comms_type )    :: comm

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Rescale

  Do i=1,natms
     uuu=xxx(i)
     vvv=yyy(i)
     www=zzz(i)

     xxx(i)=rcell(1)*uuu+rcell(4)*vvv+rcell(7)*www
     yyy(i)=rcell(2)*uuu+rcell(5)*vvv+rcell(8)*www
     zzz(i)=rcell(3)*uuu+rcell(6)*vvv+rcell(9)*www

     uuu=xxx(i)
     vvv=yyy(i)
     www=zzz(i)

     xxx(i)=cels(1)*uuu+cels(4)*vvv+cels(7)*www
     yyy(i)=cels(2)*uuu+cels(5)*vvv+cels(8)*www
     zzz(i)=cels(3)*uuu+cels(6)*vvv+cels(9)*www
  End Do

! Write REVCON

  name   = 'CFGSCL' ! file name
  nstep  = 0        ! no steps done
  tstep  = 0.0_wp   ! no step exists
  time   = 0.0_wp   ! time is not relevant

  rcell = cell ; cell = cels
  Call write_config(name,lvcfscl,megatm,nstep,tstep,time,comm)
  cell = rcell

End Subroutine scale_config


Subroutine write_config(name,levcfg,megatm,nstep,tstep,time,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing configuration file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Character( Len = * ), Intent( In    ) :: name
  Integer,              Intent( In    ) :: levcfg,megatm,nstep
  Real( Kind = wp ),    Intent( In    ) :: tstep,time
  Type( comms_type ),   Intent( InOut ) :: comm

  Integer, Parameter :: recsz    =   73 ! default record size

  Logical               :: ready
  Character( Len = 40 ) :: fname
  Integer(Kind=li)      :: rec,rec1     ! record line

  Integer               :: fail(1:4),i,k,jj,jdnode,jatms

  Real( Kind = wp )     :: celprp(1:10),cell_vecs(1:3,1:3)
  Real( Kind = wp )     :: lengths(1:3), angles(1:3)

! Some parameters and variables needed by io interfaces

  Integer                           :: fh
  Integer                           :: io_write,batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: cxx,cyy,czz

  Integer :: ierr
  Character ( Len = 256 )  ::  message

! Get write method buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! Get offsets and define batch

  fail=0
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
     n_atm=0 ; n_atm(comm%idnode+1)=natms
     Call gsum(comm,n_atm)
     n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = MPI_OFFSET_KIND)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(1-1,MPI_OFFSET_KIND)
     jj=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_delete( name,comm ) ! Sort existence issues
        Call io_open( io_write, MPI_COMM_SELF, name, MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(4i10,1p,2e16.7,a1)') levcfg,imcon,megatm,nstep,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

! Accumulate header - optional cell information (if present)

        If (imcon > 0) Then
           Do i = 0, 2
              Write(record, Fmt='(3f20.10,a12,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do
           End Do
        End If

! Dump header

        Call io_write_batch( fh, rec_mpi_io, jj, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        jj=jj+2
        If (imcon > 0) jj=jj+3

     End If
     Call gsync(comm)

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_delete( name ,comm)
     Call io_open( io_write, comm%comm, name, MPI_MODE_WRONLY, fh )

! Start of file (updated)

     rec_mpi_io=Int(jj,MPI_OFFSET_KIND)+Int(n_atm(0),MPI_OFFSET_KIND)*Int(levcfg+2,MPI_OFFSET_KIND)
     jj=0
     Do i=1,natms
        Write(record, Fmt='(a8,i10,a54,a1)') atmnam(i),ltg(i),Repeat(' ',54),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') xxx(i),yyy(i),zzz(i),Repeat(' ',12),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        If (levcfg > 0) Then
           Write(record, Fmt='(3g20.10,a12,a1)') vxx(i),vyy(i),vzz(i),Repeat(' ',12),lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do

           If (levcfg > 1) Then
              Write(record, Fmt='(3g20.10,a12,a1)') fxx(i),fyy(i),fzz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do
           End If
        End If

! Dump batch and update start of file

        If (jj + levcfg + 2 >= batsz .or. i == natms) Then
           Call io_write_batch( fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,MPI_OFFSET_KIND)
           jj=0
        End If
     End Do

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'write_config allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     rec=Int(0,li)
     jj=0
     If (comm%idnode == 0) Then

! Write configuration data to new configuration file

        Open(Unit=nconf, File=name, Form='formatted', Access='direct', Recl=recsz, Status='replace')

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(4i10,1p,2e16.7,a1)') levcfg,imcon,megatm,nstep,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

! Accumulate header - optional cell information (if present)

        If (imcon > 0) Then
           Do i = 0, 2
              Write(record, Fmt='(3f20.10,a12,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do
           End Do
        End If

! Dump header and update start of file

        Write(Unit=nconf, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
        rec=Int(jj,li)
        jj=0

        Do i=1,natms
           iwrk(i)=ltg(i)
           chbuf(i)=atmnam(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           If (levcfg > 0) Then
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)

              If (levcfg > 1) Then
                 cxx(i)=fxx(i)
                 cyy(i)=fyy(i)
                 czz(i)=fzz(i)
              End If
           End If
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,WriteConf_tag,comm%comm,comm%ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)

                 If (levcfg > 0) Then
                    Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                    Call MPI_RECV(byy,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                    Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)

                    If (levcfg > 1) Then
                       Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                       Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                       Call MPI_RECV(czz,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                    End If
                 End If
              End If
           End If

           jj=0
           Do i=1,jatms
              Write(record, Fmt='(a8,i10,a54,a1)') atmnam(i),iwrk(i),Repeat(' ',54),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              If (levcfg > 0) Then
                 Write(record, Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,recsz
                    chbat(k,jj) = record(k:k)
                 End Do

                 If (levcfg > 1) Then
                    Write(record, Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
                    jj=jj+1
                    Do k=1,recsz
                       chbat(k,jj) = record(k:k)
                    End Do
                 End If
              End If

! Dump batch and update start of file

              If (jj + levcfg + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nconf, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
                 rec=rec+Int(jj,li)
                 jj=0
              End If
           End Do
        End Do

        Close(Unit=nconf)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,WriteConf_tag,comm%comm,comm%status,comm%ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,WriteConf_tag,comm%comm,comm%ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,WriteConf_tag,comm%comm,comm%ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,WriteConf_tag,comm%comm,comm%ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)

           If (levcfg > 0) Then
              Call MPI_SEND(vxx,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
              Call MPI_SEND(vyy,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
              Call MPI_SEND(vzz,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)

              If (levcfg > 1) Then
                 Call MPI_SEND(fxx,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
                 Call MPI_SEND(fyy,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
                 Call MPI_SEND(fzz,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
              End If
           End If
        End If

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'write_config deallocation failure'
        Call error(0,message)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = name
     Else
        fname = name(1:Len_Trim(name)) // '.nc'
     End If

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(1-1,MPI_OFFSET_KIND)
     jj=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_delete( fname,comm ) ! Sort existence issues
        If (io_write == IO_WRITE_SORTED_NETCDF) Call io_nc_create( MPI_COMM_SELF, fname, cfgname, megatm )
        Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header

           Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
           Call io_write_record( fh, Int(jj,MPI_OFFSET_KIND), record )
           jj=jj+1

           Write(record, Fmt='(4i10,1p,2e16.7,a1)') levcfg,imcon,megatm,nstep,tstep,time,lf
           Call io_write_record( fh, Int(jj,MPI_OFFSET_KIND), record )
           jj=jj+1

! Write optional cell information (if present)

           If (imcon > 0) Then
              Do i = 0, 2
                 Write( record, '( 3f20.10, a12, a1 )' ) &
                      cell( 1 + i * 3: 3 + i * 3 ), Repeat( ' ', 12 ), lf
                 Call io_write_record( fh, Int(jj,MPI_OFFSET_KIND), record )
                 jj=jj+1
              End Do
           End If

        Else ! netCDF write

           jj=1 ! For config there is only one frame

           Call io_nc_put_var( 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var( 'datalevel'      , fh, levcfg, jj, 1 )
           Call io_nc_put_var( 'imageconvention', fh,  imcon, jj, 1 )
           Call io_nc_put_var( 'timestep'       , fh,  tstep, jj, 1 )

           If (imcon > 0) Then
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

        End If

        Call io_close( fh )
        Call io_finalize

     Else

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           jj=jj+2
           If (imcon > 0) jj=jj+3
        Else
           jj=1
        End If

     End If
     Call gsync(comm)

! Write the rest

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, fname, MPI_MODE_WRONLY, fh )

     rec_mpi_io=rec_mpi_io+Int(jj,MPI_OFFSET_KIND)
     Call io_write_sorted_file( fh, levcfg, IO_RESTART, rec_mpi_io, natms,      &
          ltg, atmnam, (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), xxx, yyy, zzz, &
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

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'write_config allocation failure'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     rec=Int(0,li)
     If (comm%idnode == 0) Then

! Write configuration data to new configuration file

        Open(Unit=nconf, File=name, Form='formatted', Access='direct', Recl=recsz, Status='replace')

! Write header

        rec=rec+Int(1,li)
        Write(Unit=nconf, Fmt='(a72,a1)',            Rec=rec) cfgname(1:72),lf
        rec=rec+Int(1,li)
        Write(Unit=nconf, Fmt='(4i10,1p,2e16.7,a1)', Rec=rec) levcfg,imcon,megatm,nstep,tstep,time,lf

! Write optional cell information (if present)

        If (imcon > 0) Then
           Do i = 0, 2
              rec=rec+Int(1,li)
              Write(Unit=nconf, Fmt='(3f20.10,a12,a1)', Rec=rec) &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           End Do
        End If

        Do i=1,natms
           iwrk(i)=ltg(i)
           chbuf(i)=atmnam(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           If (levcfg > 0) Then
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)

              If (levcfg > 1) Then
                 cxx(i)=fxx(i)
                 cyy(i)=fyy(i)
                 czz(i)=fzz(i)
              End If
           End If
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,WriteConf_tag,comm%comm,comm%ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)

                 If (levcfg > 0) Then
                    Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                    Call MPI_RECV(byy,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                    Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)

                    If (levcfg > 1) Then
                       Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                       Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                       Call MPI_RECV(czz,jatms,wp_mpi,jdnode,WriteConf_tag,comm%comm,comm%status,comm%ierr)
                    End If
                 End If
              End If
           End If

           Do i=1,jatms
              rec1=rec+Int(iwrk(i)-1,li)*Int(levcfg+2)+Int(1,li)
              Write(Unit=nconf, Fmt='(a8,i10,a54,a1)',     Rec=rec1) chbuf(i),iwrk(i),Repeat(' ',54),lf
              rec1=rec1+Int(1,li)
              Write(Unit=nconf, Fmt='(3g20.10,a12,a1)',    Rec=rec1) axx(i),ayy(i),azz(i),Repeat(' ',12),lf

              If (levcfg > 0) Then
                 rec1=rec1+Int(1,li)
                 Write(Unit=nconf, Fmt='(3g20.10,a12,a1)', Rec=rec1) bxx(i),byy(i),bzz(i),Repeat(' ',12),lf

                 If (levcfg > 1) Then
                    rec1=rec1+Int(1,li)
                    Write(Unit=nconf, Fmt='(3g20.10,a12,a1)', Rec=rec1) cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
                 End If
              End If
           End Do
        End Do

        Close(Unit=nconf)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,WriteConf_tag,comm%comm,comm%status,comm%ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,WriteConf_tag,comm%comm,comm%ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,WriteConf_tag,comm%comm,comm%ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,WriteConf_tag,comm%comm,comm%ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)

           If (levcfg > 0) Then
              Call MPI_SEND(vxx,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
              Call MPI_SEND(vyy,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
              Call MPI_SEND(vzz,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)

              If (levcfg > 1) Then
                 Call MPI_SEND(fxx,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
                 Call MPI_SEND(fyy,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
                 Call MPI_SEND(fzz,natms,wp_mpi,0,WriteConf_tag,comm%comm,comm%ierr)
              End If
           End If
        End If

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'write_config deallocation failure'
        Call error(0,message)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'write_config deallocation failure 0'
        Call error(0,message)
     End If
  End If

  Call gsync(comm)

End Subroutine write_config


Subroutine getcom(xxx,yyy,zzz,com,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate system centre of mass position
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: xxx,yyy,zzz
    Real( Kind = wp ), Dimension( 1:3 ), Intent(   Out ) :: com
    Type(comms_type), Intent ( InOut )                   :: comm

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: totmas

    Integer                 :: i

! total system mass

    If (newjob) Then
       newjob = .false.

       totmas = 0.0_wp
       Do i=1,natms
          If (lfrzn(i) == 0) totmas = totmas + weight(i)
       End Do

       Call gsum(comm,totmas)
    End If

    com = 0.0_wp

    Do i=1,natms
       If (lfrzn(i) == 0) Then
          com(1) = com(1) + weight(i)*xxx(i)
          com(2) = com(2) + weight(i)*yyy(i)
          com(3) = com(3) + weight(i)*zzz(i)
       End If
    End Do

    Call gsum(comm,com)
    If (totmas >= zero_plus) com = com/totmas

  End Subroutine getcom

  Subroutine getcom_mol(istart,ifinish,cmm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate a molecule's mass and COM
!
! istart  - the global index of the first atom of the molecule
! ifinish - the global index of the last atom of the molecule
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Integer,           Intent( In    ) :: istart,ifinish

    Real( Kind = wp ), Intent(   Out ) :: cmm(0:3)
    Type(comms_type), Intent ( InOut ) :: comm

    Integer           :: fail,i,j,k
    Real( Kind = wp ) :: mass,r(1:3)

    Real( Kind = wp ), Allocatable :: mol(:,:)

    Character ( Len = 256 )   ::  message

    fail = 0
    Allocate (mol(1:(ifinish-istart+1),0:3), Stat = fail)
    If (fail > 0) Then
       Write(message,'(/,1x,a,i0)') 'getcom_mol allocation failure'
       Call error(0,message)
    End If

! Initialise

    mass = 0.0_wp
    cmm  = 0.0_wp

    mol = 0.0_wp
    Do i=1,natms
       j=ltg(i)
       If (j >= istart .and. j <= ifinish) Then
          k=j-istart+1

          mol(k,0) = weight(i)
          mol(k,1) = xxx(i)
          mol(k,2) = yyy(i)
          mol(k,3) = zzz(i)
       End If
    End Do

    Call gsum(comm,mol)

    r(1) = mol(1,1)
    r(2) = mol(1,2)
    r(3) = mol(1,3)

    mol(:,1) = mol(:,1)-r(1)
    mol(:,2) = mol(:,2)-r(2)
    mol(:,3) = mol(:,3)-r(3)

    k=ifinish-istart+1
    Call images(imcon,cell,k,mol(:,1),mol(:,2),mol(:,3))

    mol(:,1) = mol(:,1)+r(1)
    mol(:,2) = mol(:,2)+r(2)
    mol(:,3) = mol(:,3)+r(3)

    Do i=1,k
       mass   = mass   + mol(i,0)
       cmm(0) = cmm(0) + mol(i,0)*Real(1-lfrzn(i),wp)
       cmm(1) = cmm(1) + mol(i,0)*mol(i,1)
       cmm(2) = cmm(2) + mol(i,0)*mol(i,2)
       cmm(3) = cmm(3) + mol(i,0)*mol(i,3)
    End Do

    If (cmm(0) >= zero_plus) cmm(1:3) = cmm(1:3) / mass

    fail = 0
    Deallocate (mol, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'getcom_mol deallocation failure'
       Call error(0,message)
    End If

  End Subroutine getcom_mol


  Subroutine freeze_atoms()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to quench forces and velocities on 'frozen' atoms
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Integer :: i

    Do i=1,natms
       If (lfrzn(i) /= 0) Then
           vxx(i) = 0.0_wp ; vyy(i) = 0.0_wp ; vzz(i) = 0.0_wp
           fxx(i) = 0.0_wp ; fyy(i) = 0.0_wp ; fzz(i) = 0.0_wp
       End If
    End Do

  End Subroutine freeze_atoms

  Subroutine origin_config(megatm,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for translating the origin of the MD box as
  ! defined in CONFIG by the (xorg,yorg,zorg) vector and saving it in
  ! CFGORG
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: megatm
    Type( comms_type ), Intent( InOut ) :: comm

    Character ( Len = 6 ) :: name
    Integer               :: i,nstep
    Real( Kind = wp )     :: tstep,time

  ! Translate

    Do i=1,natms
       xxx(i)=xxx(i)+xorg
       yyy(i)=yyy(i)+yorg
       zzz(i)=zzz(i)+zorg
    End Do

  ! Restore periodic boundaries

    Call pbcshift(imcon,cell,natms,xxx,yyy,zzz)

  ! Write REVCON

    name   = 'CFGORG' ! file name
    nstep  = 0        ! no steps done
    tstep  = 0.0_wp   ! no step exists
    time   = 0.0_wp   ! time is not relevant

    Call write_config(name,lvcforg,megatm,nstep,tstep,time,comm)

  End Subroutine origin_config
End Module configuration
