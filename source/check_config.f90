Subroutine check_config &
           (levcfg,imcon,l_str,lpse,keyens,iso,keyfce,keyres,megatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reporting and checking the configuration data
! in CONFIG for: (i) unsuitable simulation options in CONTROL and
! (ii) connectivity to FIELD; before connecting the crystallographic
! data (positions+) to the topology (sites+), i.e. CONFIG to FIELD
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use site_module
  Use config_module

  Implicit None

  Logical, Intent( In    ) :: l_str,lpse
  Integer, Intent( In    ) :: levcfg,keyens,iso,keyfce,keyres,megatm
  Integer, Intent( InOut ) :: imcon

  Logical                :: safe
  Integer                :: fail,k,l,m, &
                            indatm,nattot,mol_sit,loc_ind
  Real( Kind = wp )      :: rcell(1:9),det

  Integer, Allocatable :: iwrk(:)

  fail=0
  Allocate (iwrk(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'check_config allocation failure, node: ', idnode
     Call error(0)
  End If


  If (idnode == 0) Then
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

  If (idnode == 0) Then
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

  nattot=0

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

           nattot=nattot+1

! If a local atom has a global index nattot

           If (lsa(indatm) == nattot) Then

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

              If (l_str) iwrk(indatm) = nattot ! Populate

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

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(25)

! Check CONFIG indices eligibility

  If (l_str) Then
     Call all_inds_present( iwrk, indatm, megatm, safe )
     If (mxnode > 1) Call gcheck(safe)
     If (.not.safe) Call error(28)
  End If

  Deallocate (iwrk, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'check_config deallocation failure, node: ', idnode
     Call error(0)
  End If

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

    me    = idnode
    nproc = mxnode

    ! Ever the optimist !
    all_present = .true.

    ! If the number on each proc do not add up to N there is no
    ! way things can be correct. will need all the local totals so
    ! gather and sum locally rather than just a global sum as
    ! saves a comm.
    Allocate( all_n_loc( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    ! No check on mpi error as all is pointless because
    ! of the way mpi error handlers are usually dealt with
!
!    all_n_loc = 0
!    all_n_loc( idnode ) = n_loc
!    Call gsum( all_n_loc( 0:nrpocs - 1 ) )
!
    Call MPI_ALLGATHER(     n_loc, 1, MPI_INTEGER, &
                        all_n_loc, 1, MPI_INTEGER, dlp_comm_world, ierr )
    all_present = ( Sum( all_n_loc ) == n )
    If ( .not. all_present ) Return

    ! Work out the first index on each proc
    Allocate( loc_start( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    loc_start( 0 ) = 1
    Do i = 1, nproc - 1
       loc_start( i ) = loc_start( i - 1 ) + all_n_loc( i - 1 )
    End Do

    ! Array to work with local indices
    Allocate( local_ind( 1:n_loc ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    local_ind = ind( 1:n_loc )

    ! Sort the local data
    Call shellsort( n_loc, local_ind )

    ! Work out how much data to send to each processor
    Allocate( to_send( 0:nproc - 1 ), Stat = fail )
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
    Allocate( to_recv( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    Call MPI_ALLTOALL( to_send, 1, MPI_INTEGER, &
                       to_recv, 1, MPI_INTEGER, &
                       dlp_comm_world, ierr )

    ! Work out the displacements in the sending and receiving arrays
    Allocate( displs_send( 0:nproc - 1 ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    Allocate( displs_recv( 0:nproc - 1 ), Stat = fail )
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
    Allocate( reorg_ind( 1:n_loc ), Stat = fail )
    If ( fail /= 0 ) Go To 100

    Call MPI_ALLTOALLV( local_ind, to_send, displs_send, MPI_INTEGER, &
                        reorg_ind, to_recv, displs_recv, MPI_INTEGER, &
                        dlp_comm_world, ierr )

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
    Call MPI_ALLREDUCE( loc_present, all_present, 1, MPI_LOGICAL, MPI_LAND, dlp_comm_world, ierr )

    Deallocate( reorg_ind   , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( displs_recv , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( displs_send , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( to_recv     , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( to_send     , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( local_ind   , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( loc_start   , Stat = fail )
    If ( fail /= 0 ) Go To 100
    Deallocate( all_n_loc   , Stat = fail )
    If ( fail /= 0 ) Go To 100

    Return

100 Continue
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'all_inds_present allocation/deallocation failure, node: ', idnode
       Call error(0)
    End If

  End Subroutine all_inds_present

End Subroutine check_config
