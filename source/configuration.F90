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

  Use kinds, Only : wp
  Use comms, Only : comms_type
  Use setup_module
  Use site_module
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

  Interface reallocate
     Module Procedure reallocate_chr_v
     Module Procedure reallocate_int_v
     Module Procedure reallocate_rwp_v
  End Interface

  Private :: reallocate_chr_v, reallocate_int_v, reallocate_rwp_v

Contains

  Subroutine reallocate_chr_v( delta, a, stat )

    Implicit None

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

    Implicit None

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

    Implicit None

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

    Implicit None

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

    Use setup_module

    Implicit None

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

  fail=0
  If (l_str) Then
     Allocate (iwrk(1:mxatms), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'check_config allocation failure, node: ', comm%idnode
        Call error(0)
    End If
  End If


  If (comm%idnode == 0 .and. newjob) Then
     Write(nrite,"(/,1x,'configuration file name: ',/,/,10x,a)") cfgname
     Write(nrite,"(/,/,1x,'selected image convention',6x,i10)") imcon
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

  If (comm%idnode == 0 .and. newjob) Then
     Write(nrite,"(/,/,1x,'simulation cell vectors'/)")
     Write(nrite,"(3f20.10)") cell
     Write(nrite,"(/,/,1x,'system volume     ',2x,1p,g22.12)") det
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
                 Write(nrite,"(/,/,1x, 'unidentified atom label :',a8,': atom number ',i5)") atmnam(loc_ind),loc_ind
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
        Write(nrite,'(/,1x,a,i0)') 'check_config deallocation failure, node: ', comm%idnode
        Call error(0)
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

    Implicit None

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
       Write(nrite,'(/,1x,a,i0)') 'all_inds_present allocation/deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

  End Subroutine all_inds_present

End Subroutine check_config


End Module configuration
