Subroutine exchange_grid( ixb , ixt , iyb , iyt , izb , izt , qqc_local , &
                          ixdb, iydb, izdb, ixdt, iydt, izdt, qqc_domain  )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for exchanging grid data post DaFT
!
! Get all the data required to calculate the reciprocal space
! contribution to the forces on the atoms held by this processor, and
! send what is required by the other processors to them.
! QQC_LOCAL( IXB:IXT, IYB:IYT, IZB:IZT ) holds all the data I have currently.
! QQC_DOMAIN( IXDB:IXDT, IYDB:IYDT, IZDB:IZDT ) will hold all the data I need.
! MAP tells me about my neighbouring processors (see the DOMAINS_MODULE).
! MXSPL is the order of the spline.
! MXSPL1 is the extended spline should the conditional VNL is at work.
! MXSPL2 covers for MXSPL1 and if particles are not domain bound.
!
! copyright - daresbury laboratory
! author    - i.j.bush & i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module
  Use setup_module,   Only : nrite,mxspl,mxspl2
  Use domains_module, Only : map

  Implicit None

  Integer,           Intent( In    ) :: ixb , iyb , izb
  Integer,           Intent( In    ) :: ixt , iyt , izt
  Integer,           Intent( In    ) :: ixdb, iydb, izdb
  Integer,           Intent( In    ) :: ixdt, iydt, izdt
  Real( Kind = wp ), Intent( In    ) :: qqc_local(   ixb:ixt ,  iyb:iyt ,  izb:izt )
  Real( Kind = wp ), Intent(   Out ) :: qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt )

  Integer :: lx, ly, lz, delspl
  Integer :: me

  delspl=mxspl2-mxspl

! What's my name?

  me = idnode

! Find length of sides of the domain

  lx = ixt - ixb + 1
  ly = iyt - iyb + 1
  lz = izt - izb + 1

! Copy over our local data

  qqc_domain( ixb:ixt, iyb:iyt, izb:izt ) = qqc_local

  If (delspl == 0) Then

! Note that because of the way the splines work when particles don't
! blur off domains (mxspl1==mxspl), i.e. no conditional VNL updates,
! and are bound (mxspl2==mxspl) I only require to receive data from
! processors in the negative octant relative to me, and hence only
! need to send data to processors in the positive octant.

! +X direction face - negative halo

     Call exchange_grid_halo( map(1),                     map(2), &
          mxspl2,                      ly,                          lz, &
          ixt-mxspl2+1, ixt  ,         iyb, iyt,                    izb, izt, &
          ixdb        , ixb-1,         iyb, iyt,                    izb, izt )

! +Y direction face (including the +X face extension) - negative halo

     Call exchange_grid_halo( map(3),                     map(4), &
          lx+mxspl2,                   mxspl2,                      lz, &
          ixdb, ixt,                   iyt-mxspl2+1, iyt  ,         izb, izt, &
          ixdb, ixt,                   iydb        , iyb-1,         izb, izt )

! +Z direction face (including the +Y+X faces extensions) - negative halo

     Call exchange_grid_halo( map(5),                     map(6), &
          lx+mxspl2,                   ly+mxspl2,                   mxspl2, &
          ixdb, ixt,                   iydb, iyt,                   izt-mxspl2+1, izt, &
          ixdb, ixt,                   iydb, iyt,                   izdb        , izb-1 )

  Else

! +X direction face - negative halo

     Call exchange_grid_halo( map(1),                     map(2), &
          mxspl2,                      ly,                          lz, &
          ixt-mxspl2+1, ixt  ,         iyb, iyt,                    izb, izt, &
          ixdb        , ixb-1,         iyb, iyt,                    izb, izt )
!          (ixt)-(ixt-mxspl2+1)+1=mxspl2
!          (ixb-1)-(ixdb)+1=mxspl2

! -X direction face - positive halo

     Call exchange_grid_halo( map(2),                     map(1), &
          delspl,                      ly,                          lz, &
          ixb          , ixb+delspl-1, iyb, iyt,                    izb, izt, &
          ixdt-delspl+1, ixdt        , iyb, iyt,                    izb, izt )
!          (ixb+delspl-1)-(ixb)+1=delspl
!          (ixdt)-(ixdt-delspl+1)+1=delspl


! +Y direction face (including the +&-X faces extensions) - negative halo

     Call exchange_grid_halo( map(3),                     map(4), &
          lx+mxspl2+delspl,            mxspl2,                      lz, &
          ixdb, ixdt,                  iyt-mxspl2+1, iyt  ,         izb, izt, &
          ixdb, ixdt,                  iydb        , iyb-1,         izb, izt )

! -Y direction face (including the +&-X faces extensions) - positive halo

     Call exchange_grid_halo( map(4),                     map(3), &
          lx+mxspl2+delspl,            delspl,                      lz, &
          ixdb, ixdt,                  iyb          , iyb+delspl-1, izb, izt, &
          ixdb, ixdt,                  iydt-delspl+1, iydt        , izb, izt )

! +Z direction face (including the +&-Y+&-X faces extensions) - negative halo

     Call exchange_grid_halo( map(5),                     map(6), &
          lx+mxspl2+delspl,            ly+mxspl2+delspl,            mxspl2, &
          ixdb, ixdt,                  iydb, iydt,                  izt-mxspl2+1, izt, &
          ixdb, ixdt,                  iydb, iydt,                  izdb        , izb-1 )

! -Z direction face (including the +&-Y+&-X faces extensions) - positive halo

     Call exchange_grid_halo( map(6),                     map(5), &
          lx+mxspl2+delspl,            ly+mxspl2+delspl,            delspl, &
          ixdb, ixdt,                  iydb, iydt,                  izb          , izb+delspl-1, &
          ixdb, ixdt,                  iydb, iydt,                  izdt-delspl+1, izdt         )

  End If

Contains

  Subroutine exchange_grid_halo(     from,       to,           &
                                 lx,       ly,       lz,       &
                                 xlb, xlt, ylb, ylt, zlb, zlt, &
                                 xdb, xdt, ydb, ydt, zdb, zdt )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for exchanging grid data post DaFT
!
! Receives data from a processor FROM and sticks it in
!                                QQC( XDB:XDT, YDB:YDT, ZDB:ZDT )
! Sends data to a processor      TO from
!                                QQC( XLB:XLT, YLB:YLT, ZLB:ZLT )
!
! Note: Amount of data sent M-U-S-T be the same as that received!!!
!
! If properly written LX, LY, LZ should not be there, but Ian finds
! this kind of thing confusing so they are included to make his life
! easy.  They are the dimensions of bits of the arrays, i.e.
! LX = XDT - XDB + 1 = XLT - XLB + 1, etc. for LY, LZ.

! copyright - daresbury laboratory
! author    - i.j.bush & i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    Integer :: fail(1:2)

! If the processor to receive FROM is actually ME it means there is
! only one processor along this axis (so the processor to send TO is
! also ME) and so no message passing need be done.  However, there
! is a catch.  The domain decomposed spme_forces routine expects data
! that would be `seen' through the periodic boundary conditions to be
! actually copied from the high positive indices to negative ones and
! vice-versa.  The Else clause catches this.

    If ( from /= me ) Then

! Allocate send and receive buffers (of the same size!!!)
! so all can be sent and received as one message!!!

       fail=0
       Allocate ( send_buffer( xlb:xlt, ylb:ylt, zlb:zlt ) , Stat = fail(1) )
       Allocate ( recv_buffer( xdb:xdt, ydb:ydt, zdb:zdt ) , Stat = fail(2) )
       If (Any(fail > 0)) Then
          Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo allocation failure, node: ', idnode
          Call error(0)
       End If

! Copy the data to be sent

       send_buffer = qqc_domain( xlb:xlt, ylb:ylt, zlb:zlt )

! Length of message to send is the same as that to receive

       length = lx * ly * lz

! Exchange the data

       Call MPI_IRECV( recv_buffer, length, wp_mpi, from, ExchgGrid_tag, dlp_comm_world, request, ierr )
       Call MPI_SEND(  send_buffer, length, wp_mpi, to  , ExchgGrid_tag, dlp_comm_world, ierr )
       Call MPI_WAIT(  request, status, ierr )

! Copy the received data into the domain halo

       qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = recv_buffer

! And, as my mum told me, leave things as we found them

       Deallocate ( recv_buffer , Stat = fail(1) )
       Deallocate ( send_buffer , Stat = fail(2) )
       If (Any(fail > 0)) Then
          Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo deallocation failure, node: ', idnode
          Call error(0)
       End If

    Else

! Simple on node copy - as sizes are the same as 3D shapes

       qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = qqc_domain( xlb:xlt, ylb:ylt, zlb:zlt )

    End If

  End Subroutine exchange_grid_halo

End Subroutine exchange_grid
