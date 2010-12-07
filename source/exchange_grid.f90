Subroutine exchange_grid( ixb , ixt, iyb, iyt, izb, izt,          &
                          ixdb, iydb, izdb, qqc_local, qqc_domain )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for excahnging grid data post DaFT
!
! Get all the data required to calculate the reciprocal space contribution
! to the forces on the atoms held by this processor, and send what is required
! by the other procs to them.
! QQC_LOCAL( IXB:IXT, IYB:IYT, IZB:IZT ) holds all the data I have currently.
! QQC_LOCAL( IXDB:IXT, IYDB:IYT, IZDB:IZT ) will hold all the data I need
! MAP tells me about my neighbouring processors ( see the DOMAINS_MODULE )
! MXSPL is the order of the spline
!
! Note that because of the way the splines work I only require to recv
! data from procs in the negative octant relative to me, and hence only
! need to send data to procs in the positive octant.
!
! copyright - daresbury laboratory
! author    - i.j.bush march 2002
! amended   - i.t.todorov november 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module, Only : map

  Implicit None

  Integer,           Intent( In    ) :: ixb , iyb , izb
  Integer,           Intent( In    ) :: ixt , iyt , izt
  Integer,           Intent( In    ) :: ixdb, iydb, izdb
  Real( Kind = wp ), Intent( In    ) :: qqc_local(  ixb :ixt, iyb :iyt, izb :izt )
  Real( Kind = wp ), Intent(   Out ) :: qqc_domain( ixdb:ixt, iydb:iyt, izdb:izt )

  Integer :: lx, ly, lz
  Integer :: me

  ! What's my name ?
  me = idnode

  ! Find length of sides of the domain
  lx = ixt - ixb + 1
  ly = iyt - iyb + 1
  lz = izt - izb + 1

  ! Copy over our local data
  qqc_domain( ixb:ixt, iyb:iyt, izb:izt ) = qqc_local

  ! Now have to recv the `halo' data due to the finite length of the splines
  ! with the ( up to ) 7 processors in the negative octant relative to this
  ! proc, whilst send to the positive octant the halo that they require.

  ! The messages split into 3 groups:
  ! a) 3 `big' messages in the X, Y, Z directions
  ! b) 3 `small' messages in the XY, XZ, YZ
  ! c) 1 `tiny' message in the XYZ
  ! Actually big here is relative. For the X transfer if
  ! the proc grid is Px*Py*Pz, the FFT grid is Nx*Ny*Nz
  ! and the number of points in the spline is is the
  ! size of the transfer is s * Nx * Ny / ( Px * Py )
  ! and analgously for the Y and Z. For the XY transfer
  ! the size is s * s * Nx / Px, and for the XYZ it is
  ! s * s * s. As S is usually small the size of the
  ! messages drops rapidly.

  ! X direction
  Call exchange_grid_halo( map( 1 ), map( 2 ), mxspl, ly, lz, &
       ixt - mxspl + 1, ixt    , iyb, iyt, izb, izt, &
       ixdb           , ixb - 1, iyb, iyt, izb, izt )

  ! Y direction
  Call exchange_grid_halo( map( 3 ), map( 4 ), lx, mxspl, lz, &
       ixb, ixt, iyt - mxspl + 1, iyt    , izb, izt, &
       ixb, ixt, iydb           , iyb - 1, izb, izt )

  ! Z direction
  Call exchange_grid_halo( map( 5 ), map( 6 ), lx, ly, mxspl, &
       ixb, ixt, iyb, iyt, izt - mxspl + 1, izt,              &
       ixb, ixt, iyb, iyt, izdb           , izb - 1 )

  ! XY Direction
  Call exchange_grid_halo( map( 9 ), map( 10 ), mxspl, mxspl, lz,    &
       ixt - mxspl + 1, ixt    , iyt - mxspl + 1, iyt    , izb, izt, &
       ixdb           , ixb - 1, iydb           , iyb - 1, izb, izt )

  ! XZ Direction
  Call exchange_grid_halo( map( 13 ), map( 14 ), mxspl, ly, mxspl, &
       ixt - mxspl + 1, ixt    , iyb, iyt, izt - mxspl + 1, izt,   &
       ixdb           , ixb - 1, iyb, iyt, izdb           , izb - 1 )

  ! YZ Direction
  Call exchange_grid_halo( map( 17 ), map( 18 ), lx, mxspl, mxspl, &
       ixb, ixt, iyt - mxspl + 1, iyt    , izt - mxspl + 1, izt,   &
       ixb, ixt, iydb           , iyb - 1, izdb           , izb - 1 )

  ! XYZ Direction
  Call exchange_grid_halo( map( 19 ), map( 20 ), mxspl, mxspl, mxspl,    &
       ixt - mxspl + 1, ixt, iyt - mxspl + 1, iyt, izt - mxspl + 1, izt, &
       ixdb, ixb - 1, iydb, iyb - 1, izdb, izb - 1 )

Contains

  Subroutine exchange_grid_halo( from, to, lx, ly, lz,         &
                                 xlb, xlt, ylb, ylt, zlb, zlt, &
                                 xdb, xdt, ydb, ydt, zdb, zdt )

    ! Recvs data from proc FROM and sticks it in QQC_DOMAIN( XDB:XDT, YDB:YDT, ZDB:ZDT )
    ! Sends data to   proc TO   from QQC_LOCAL( XLB:XLT, YLB:YLT, ZLB:ZLT )
    ! N.B. Amount of data sent MUST be the same as that recv'd.
    ! If properly written LX, LY, LZ should not be there, but I find this kind of
    ! thing confusing so I included them to make my life easy. They are the dimensions
    ! of bits of the arrays, i.e. LX = XDT - XDB + 1 = XLT - XLB + 1, etc. for LY, LZ.

    Implicit None

    Integer, Intent( In    ) :: from, to
    Integer, Intent( In    ) :: lx, ly, lz
    Integer, Intent( In    ) :: xlb, ylb, zlb
    Integer, Intent( In    ) :: xlt, ylt, zlt
    Integer, Intent( In    ) :: xdb, ydb, zdb
    Integer, Intent( In    ) :: xdt, ydt, zdt

    Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: send_buffer
    Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: recv_buffer

    Integer :: length
    Integer :: recv_tag
    Integer :: fail(1:2) = 0

    ! If the proc to recv from is actually me it means there is only one proc
    ! along this axis ( so the proc to send to is also me ) and so no
    ! message passing need be done. However there is a catch .... The
    ! domain decomposed spme_for routine expects data that would be `seen'
    ! through the periodic boundary conditions to be actually copied
    ! from the high positive indices to negative ones. The Else clause
    ! catches this.
    If ( from /= me ) Then

       ! Length of message to send is the same as that to recv
       length = lx * ly * lz

       ! Copy the data to be sent into a buffer so can send all as one message
       Allocate ( send_buffer( xlb:xlt, ylb:ylt, zlb:zlt ) , stat = fail(1) )
       send_buffer = qqc_local( xlb:xlt, ylb:ylt, zlb:zlt )

       ! Now set up the recv buffer
       Allocate ( recv_buffer( xdb:xdt, ydb:ydt, zdb:zdt ) , stat = fail(2) )

       If (Any(fail > 0)) Then
          Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo allocation failure, node: ', idnode
          Call error(0)
       End If

       ! Exchange the data
       Call MPI_IRECV( recv_buffer, length, wp_mpi, from, 10, dlp_comm_world, recv_tag, ierr )
       Call MPI_SEND(  send_buffer, length, wp_mpi, to  , 10, dlp_comm_world, ierr )
       Call MPI_WAIT(  recv_tag, status, ierr )

       ! And copy the recv'd data into the domain
       qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = recv_buffer

       ! And, as my mum told me, leave things as we found them.
       Deallocate ( recv_buffer , stat = fail(1) )
       Deallocate ( send_buffer , stat = fail(2) )

       If (Any(fail > 0)) Then
          Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo deallocation failure, node: ', idnode
          Call error(0)
       End If

    Else

       ! Simple on node copy
       qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = qqc_local( xlb:xlt, ylb:ylt, zlb:zlt )

    End If

  End Subroutine exchange_grid_halo

End Subroutine exchange_grid
