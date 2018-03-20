Module ewald

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring ewald routines arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp

  Implicit None

  Private

  Type, Public :: ewald_type
    Private
    Logical, Public :: l_fce  = .true. , l_cp  = .false. , &
                       lf_fce = .true. , lf_cp = .false.

    Real( Kind = wp ), Public :: engsic = 0.0_wp ,                             &
                          e_rc = 0.0_wp , v_rc = 0.0_wp , s_rc(1:9) = 0.0_wp , &
                          e_fr = 0.0_wp , v_fr = 0.0_wp , s_fr(1:9) = 0.0_wp , &
                          ef_fr= 0.0_wp , vf_fr= 0.0_wp , sf_fr(1:9)= 0.0_wp

    Real( Kind = wp ), Public, Allocatable :: fcx(:),fcy(:),fcz(:)
    Real( Kind = wp ), Public, Allocatable :: ffx(:),ffy(:),ffz(:)

    Logical :: newjob_kall = .true., &
               newjob_kfrz = .true.

    Contains
      Private
      Procedure, Public :: check => ewald_check
      Procedure, Public :: refresh => ewald_refresh
      Final :: ewald_deallocate
  End Type ewald_type

Contains

  Subroutine ewald_allocate_kall_arrays(T)
    Use setup_module, Only : mxatms

    Class(ewald_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%fcx(1:mxatms),T%fcy(1:mxatms),T%fcz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1040)

    T%fcx = 0.0_wp
    T%fcy = 0.0_wp
    T%fcz = 0.0_wp
  End Subroutine ewald_allocate_kall_arrays

  Subroutine ewald_allocate_kfrz_arrays(T)
    Use setup_module, Only : mxatms

    Class(ewald_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%ffx(1:mxatms),T%ffy(1:mxatms),T%ffz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1070)

    T%ffx = 0.0_wp
    T%ffy = 0.0_wp
    T%ffz = 0.0_wp
  End Subroutine ewald_allocate_kfrz_arrays

  Subroutine ewald_check(T,keyens,megfrz,nsteql,nstfce,nstep)
    Class(ewald_type) :: T

    Integer, Intent( In    ) :: keyens,megfrz, &
                                nsteql,nstfce,nstep

! Full frozen-frozen evaluation is TRUE by default
! otherwise a "refresh" is applied.

    If (keyens < 20 .and. megfrz > 0) Then
       T%lf_fce = .false.

       If (T%newjob_kfrz) Then
          T%newjob_kfrz = .false.

! At restart allocate the "refresh" k-space frozen-frozen SPME arrays

          Call ewald_allocate_kfrz_arrays(T)

! Allow copying into these arrays

          T%lf_cp=.true. ! == (keyens < 20 .and. megfrz > 0)

! Force the full frozen force evaluation

          T%lf_fce=.true.
       End If

! Reinitialise

       If (T%lf_fce) Then
         T%ffx = 0.0_wp
         T%ffy = 0.0_wp
         T%ffz = 0.0_wp

         T%ef_fr = 0.0_wp
         T%vf_fr = 0.0_wp
         T%sf_fr = 0.0_wp
       End If
    End If

! Full evaluation is TRUE by default and at any 'nstfce' timestep
! otherwise a "refresh" is applied.

    If (nstfce > 1) Then
       T%l_fce = ( (nstep <  nsteql .and. Mod(nstep,nstfce) == 0) .or. &
                 (nstep >= nsteql .and. Mod(nstep-nsteql,nstfce) == 0) )

       If (T%newjob_kall) Then
          T%newjob_kall = .false.

! At restart allocate the "refresh" k-space all SPME arrays

          Call ewald_allocate_kall_arrays(T)

! Allow copying into these arrays

          T%l_cp=.true. ! == (nstfce > 1)

! Force the full force evaluation

          T%l_fce=.true.
       End If

! Reinitialise

       If (T%l_fce) Then
         T%fcx = 0.0_wp
         T%fcy = 0.0_wp
         T%fcz = 0.0_wp

         T%e_rc = 0.0_wp
         T%v_rc = 0.0_wp
         T%s_rc = 0.0_wp

         T%e_fr = 0.0_wp
         T%v_fr = 0.0_wp
         T%s_fr = 0.0_wp
       End If
    End If

  End Subroutine ewald_check

  Subroutine ewald_refresh(T,engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stress)
    Use configuration, Only : natms,fxx,fyy,fzz

    Class( ewald_type )                :: T
    Real( Kind = wp ), Intent( InOut ) :: engcpe_rc,vircpe_rc, &
                                          engcpe_fr,vircpe_fr, &
                                          stress(1:9)

    Integer :: i

    Do i=1,natms
       fxx(i)=fxx(i)+T%fcx(i)
       fyy(i)=fyy(i)+T%fcy(i)
       fzz(i)=fzz(i)+T%fcz(i)
    End Do

    engcpe_rc=T%e_rc
    vircpe_rc=T%v_rc
    stress=stress+T%s_rc

    engcpe_fr=T%e_fr
    vircpe_fr=T%v_fr
    stress=stress+T%s_fr
  End Subroutine ewald_refresh

  Subroutine ewald_deallocate(T)
  !> Deallocate ewald_type arrays
    Type(ewald_type) :: T

    If (Allocated(T%fcx)) Then
      Deallocate(T%fcx)
    End If
    If (Allocated(T%fcy)) Then
      Deallocate(T%fcy)
    End If
    If (Allocated(T%fcz)) Then
      Deallocate(T%fcz)
    End If

    If (Allocated(T%ffx)) Then
      Deallocate(T%ffx)
    End If
    If (Allocated(T%ffy)) Then
      Deallocate(T%ffy)
    End If
    If (Allocated(T%ffz)) Then
      Deallocate(T%ffz)
    End If
  End Subroutine ewald_deallocate

  Subroutine exchange_grid( ixb , ixt , iyb , iyt , izb , izt , qqc_local , &
                            ixdb, iydb, izdb, ixdt, iydt, izdt, qqc_domain, &
                            comm)

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

    Use comms, Only : comms_type
    Use setup_module,   Only : nrite,mxspl,mxspl2
    Use domains_module, Only : map

    Integer,            Intent( In    ) :: ixb , iyb , izb
    Integer,            Intent( In    ) :: ixt , iyt , izt
    Integer,            Intent( In    ) :: ixdb, iydb, izdb
    Integer,            Intent( In    ) :: ixdt, iydt, izdt
    Real( Kind = wp ),  Intent( In    ) :: qqc_local(   ixb:ixt ,  iyb:iyt ,  izb:izt )
    Real( Kind = wp ),  Intent(   Out ) :: qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt )
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: lx, ly, lz, delspl
    Integer :: me

    delspl=mxspl2-mxspl

  ! What's my name?

    me = comm%idnode

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
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

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
            Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo allocation failure, node: ', comm%idnode
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
            Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo deallocation failure, node: ', comm%idnode
            Call error(0)
         End If

      Else

  ! Simple on node copy - as sizes are the same as 3D shapes

         qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = qqc_domain( xlb:xlt, ylb:ylt, zlb:zlt )

      End If

    End Subroutine exchange_grid_halo

  End Subroutine exchange_grid
End Module ewald
