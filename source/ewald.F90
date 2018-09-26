Module ewald

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring ewald routines arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp,wi
  Use comms,           Only : ExchgGrid_tag,comms_type,wp_mpi,gsend,gwait, &
                              girecv
  Use setup,           Only : mxatms,nrite,twopi
  Use configuration,   Only : configuration_type
  Use particle,           Only : corePart
  Use domains,         Only : domains_type
  Use errors_warnings, Only : error
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

    !> SPME FFT B-spline order
    Integer( Kind = wi ), Public :: bspline
    !> SPME FFT B-spline order when padding radius > 0
    Integer( Kind = wi ), Public :: bspline1
    !> And another one
    Integer( Kind = wi ), Public :: bspline2
    !> SPME FFT array dimensions
    Integer( Kind = wi ), Public :: fft_dim_a, fft_dim_b, fft_dim_c
    !> Original(?) SPME FFT array dimensions
    Integer( Kind = wi ), Public :: fft_dim_a1, fft_dim_b1, fft_dim_c1

  Contains
    Private
    Procedure :: ewald_allocate_kall_arrays
    Procedure :: ewald_allocate_kfrz_arrays
    Procedure, Public :: check => ewald_check
    Procedure, Public :: refresh => ewald_refresh
    Final :: ewald_deallocate
  End Type ewald_type

  Public :: bspcoe, bspgen, bspgen_mpoles, dtpbsp, spl_cexp,  exchange_grid

Contains

  Subroutine ewald_allocate_kall_arrays(T)

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

    Class(ewald_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%ffx(1:mxatms),T%ffy(1:mxatms),T%ffz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1070)

    T%ffx = 0.0_wp
    T%ffy = 0.0_wp
    T%ffz = 0.0_wp
  End Subroutine ewald_allocate_kfrz_arrays

  Subroutine ewald_check(T,ensemble,megfrz,nsteql,nstfce,nstep)
    Class(ewald_type) :: T

    Integer, Intent( In    ) :: ensemble,megfrz, &
                                nsteql,nstfce,nstep

! Full frozen-frozen evaluation is TRUE by default
! otherwise a "refresh" is applied.

    If (ensemble < 20 .and. megfrz > 0) Then
       T%lf_fce = .false.

       If (T%newjob_kfrz) Then
          T%newjob_kfrz = .false.

! At restart allocate the "refresh" k-space frozen-frozen SPME arrays

          Call T%ewald_allocate_kfrz_arrays()

! Allow copying into these arrays

          T%lf_cp=.true. ! == (ensemble < 20 .and. megfrz > 0)

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

  Subroutine ewald_refresh(T,engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stress,config)

    Class( ewald_type ), Intent( InOut ) :: T
    Real( Kind = wp ),   Intent( InOut ) :: engcpe_rc,vircpe_rc, &
                                          engcpe_fr,vircpe_fr, &
                                          stress(1:9)
    Type( configuration_type ),    Intent( InOut ) :: config

    Integer :: i

    Do i=1,config%natms
       config%parts(i)%fxx=config%parts(i)%fxx+T%fcx(i)
       config%parts(i)%fyy=config%parts(i)%fyy+T%fcy(i)
       config%parts(i)%fzz=config%parts(i)%fzz+T%fcz(i)
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
                            domain, ewld, comm)

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

    Integer,            Intent( In    ) :: ixb , iyb , izb
    Integer,            Intent( In    ) :: ixt , iyt , izt
    Integer,            Intent( In    ) :: ixdb, iydb, izdb
    Integer,            Intent( In    ) :: ixdt, iydt, izdt
    Real( Kind = wp ),  Intent( In    ) :: qqc_local(   ixb:ixt ,  iyb:iyt ,  izb:izt )
    Real( Kind = wp ),  Intent(   Out ) :: qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt )
    Type( domains_type ), Intent( In    ) :: domain
    Type( ewald_type ), Intent( In    ) :: ewld
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: lx, ly, lz, delspl
    Integer :: me

    delspl=ewld%bspline2-ewld%bspline

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
  ! blur off domains (ewld%bspline1==ewld%bspline), i.e. no conditional VNL updates,
  ! and are bound (ewld%bspline2==ewld%bspline) I only require to receive data from
  ! processors in the negative octant relative to me, and hence only
  ! need to send data to processors in the positive octant.

  ! +X direction face - negative halo

       Call exchange_grid_halo( domain%map(1),                     domain%map(2), &
            ixt-ewld%bspline2+1, ixt  ,         iyb, iyt,                    izb, izt, &
            ixdb        , ixb-1,         iyb, iyt,                    izb, izt )

  ! +Y direction face (including the +X face extension) - negative halo

       Call exchange_grid_halo( domain%map(3),                     domain%map(4), &
            ixdb, ixt,                   iyt-ewld%bspline2+1, iyt  ,         izb, izt, &
            ixdb, ixt,                   iydb        , iyb-1,         izb, izt )

  ! +Z direction face (including the +Y+X faces extensions) - negative halo

       Call exchange_grid_halo( domain%map(5),                     domain%map(6), &
            ixdb, ixt,                   iydb, iyt,                   izt-ewld%bspline2+1, izt, &
            ixdb, ixt,                   iydb, iyt,                   izdb        , izb-1 )

    Else

  ! +X direction face - negative halo

       Call exchange_grid_halo( domain%map(1),                     domain%map(2), &
            ixt-ewld%bspline2+1, ixt  ,         iyb, iyt,                    izb, izt, &
            ixdb        , ixb-1,         iyb, iyt,                    izb, izt )
  !          (ixt)-(ixt-ewld%bspline2+1)+1=ewld%bspline2
  !          (ixb-1)-(ixdb)+1=ewld%bspline2

  ! -X direction face - positive halo

       Call exchange_grid_halo( domain%map(2),                     domain%map(1), &
            ixb          , ixb+delspl-1, iyb, iyt,                    izb, izt, &
            ixdt-delspl+1, ixdt        , iyb, iyt,                    izb, izt )
  !          (ixb+delspl-1)-(ixb)+1=delspl
  !          (ixdt)-(ixdt-delspl+1)+1=delspl


  ! +Y direction face (including the +&-X faces extensions) - negative halo

       Call exchange_grid_halo( domain%map(3),                     domain%map(4), &
            ixdb, ixdt,                  iyt-ewld%bspline2+1, iyt  ,         izb, izt, &
            ixdb, ixdt,                  iydb        , iyb-1,         izb, izt )

  ! -Y direction face (including the +&-X faces extensions) - positive halo

       Call exchange_grid_halo( domain%map(4),                     domain%map(3), &
            ixdb, ixdt,                  iyb          , iyb+delspl-1, izb, izt, &
            ixdb, ixdt,                  iydt-delspl+1, iydt        , izb, izt )

  ! +Z direction face (including the +&-Y+&-X faces extensions) - negative halo

       Call exchange_grid_halo( domain%map(5),                     domain%map(6), &
            ixdb, ixdt,                  iydb, iydt,                  izt-ewld%bspline2+1, izt, &
            ixdb, ixdt,                  iydb, iydt,                  izdb        , izb-1 )

  ! -Z direction face (including the +&-Y+&-X faces extensions) - positive halo

       Call exchange_grid_halo( domain%map(6),                     domain%map(5), &
            ixdb, ixdt,                  iydb, iydt,                  izb          , izb+delspl-1, &
            ixdb, ixdt,                  iydb, iydt,                  izdt-delspl+1, izdt         )

    End If

  Contains

    Subroutine exchange_grid_halo(     from,       to,           &
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
  ! copyright - daresbury laboratory
  ! author    - i.j.bush & i.t.todorov june 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Integer, Intent( In    ) :: from, to
      Integer, Intent( In    ) :: xlb, ylb, zlb
      Integer, Intent( In    ) :: xlt, ylt, zlt
      Integer, Intent( In    ) :: xdb, ydb, zdb
      Integer, Intent( In    ) :: xdt, ydt, zdt

      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: send_buffer
      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: recv_buffer

      Integer :: fail(1:2)

      Character ( Len = 256 )  ::  message

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
            Write(message,'(a)') 'exchange_grid_halo allocation failure'
            Call error(0,message)
         End If

  ! Copy the data to be sent

         send_buffer = qqc_domain( xlb:xlt, ylb:ylt, zlb:zlt )

  ! Exchange the data

         Call girecv(comm,recv_buffer(:,:,:),from,ExchgGrid_tag)
         Call gsend(comm,send_buffer(:,:,:),to,ExchgGrid_tag)
         Call gwait(comm)

  ! Copy the received data into the domain halo

         qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = recv_buffer

  ! And, as my mum told me, leave things as we found them

         Deallocate ( recv_buffer , Stat = fail(1) )
         Deallocate ( send_buffer , Stat = fail(2) )
         If (Any(fail > 0)) Then
            Write(message,'(a)') 'exchange_grid_halo deallocation failure'
            Call error(0,message)
         End If

      Else

  ! Simple on node copy - as sizes are the same as 3D shapes

         qqc_domain( xdb:xdt, ydb:ydt, zdb:zdt ) = qqc_domain( xlb:xlt, ylb:ylt, zlb:zlt )

      End If

    End Subroutine exchange_grid_halo

  End Subroutine exchange_grid

  Subroutine bspcoe(ewld,csp,bscx,bscy,bscz,ww1,ww2,ww3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to calculate B-spline coefficients for
  ! Euler exponential splines
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith july 1998
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( ewald_type ), Intent( In    ) :: ewld
    Real( Kind = wp ),       Intent(   Out ) :: csp(1:ewld%bspline)
    Complex( Kind = wp ),    Intent( In    ) :: ww1(1:ewld%fft_dim_a),ww2(1:ewld%fft_dim_b),ww3(1:ewld%fft_dim_c)
    Complex( Kind = wp ),    Intent(   Out ) :: bscx(1:ewld%fft_dim_a),bscy(1:ewld%fft_dim_b),bscz(1:ewld%fft_dim_c)

    Integer              :: i,j,k
    Complex( Kind = wp ) :: ccc

  ! calculate B-splines at knots

    csp(1)=0.0_wp
    csp(2)=1.0_wp

    Do k=3,ewld%bspline
       csp(k)=0.0_wp

       Do j=k,2,-1
          csp(j)=(Real(j-1,wp)*csp(j)+Real(k-j+1,wp)*csp(j-1))/Real(k-1,wp)
       End Do
    End Do

  ! calculate B-spline coefficients

    Do i=0,ewld%fft_dim_a-1
       ccc=(0.0_wp,0.0_wp)

       Do k=0,ewld%bspline-2
          ccc=ccc+csp(k+2)*ww1(Mod(i*k,ewld%fft_dim_a)+1)
       End Do

       bscx(i+1)=ww1(Mod(i*(ewld%bspline-1),ewld%fft_dim_a)+1)/ccc
    End Do

    Do i=0,ewld%fft_dim_b-1
       ccc=(0.0_wp,0.0_wp)

       Do k=0,ewld%bspline-2
          ccc=ccc+csp(k+2)*ww2(Mod(i*k,ewld%fft_dim_b)+1)
       End Do

       bscy(i+1)=ww2(Mod(i*(ewld%bspline-1),ewld%fft_dim_b)+1)/ccc
    End Do

    Do i=0,ewld%fft_dim_c-1
       ccc=(0.0_wp,0.0_wp)

       Do k=0,ewld%bspline-2
          ccc=ccc+csp(k+2)*ww3(Mod(i*k,ewld%fft_dim_c)+1)
       End Do

       bscz(i+1)=ww3(Mod(i*(ewld%bspline-1),ewld%fft_dim_c)+1)/ccc
    End Do

  End Subroutine bspcoe

  Subroutine bspgen(config,nospl,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,ewld,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to calculate B-splines for SPME method
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith july 1998
  ! amended   - i.t.todorov april 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( ewald_type ), Intent( In    ) :: ewld
    Integer, Intent( In    ) :: nospl
    Type( configuration_type ),               Intent( In    ) :: config
    Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: txx,tyy,tzz

    Real( Kind = wp ), Dimension( 1:ewld%bspline , 1:mxatms ), Intent(   Out ) :: &
                                               bsdx,bsdy,bsdz,bspx,bspy,bspz
    Type( comms_type ),                       Intent( In    ) :: comm

    Integer           :: fail,i,j,k
    Real( Kind = wp ) :: aaa,bbb,ccc, rix0,riy0,riz0, jm1_r,k_r,km1_rr

    Real( Kind = wp ), Dimension( : ), Allocatable :: real_no, inv_no

    Character ( Len = 256 )   ::  message

    fail=0
    Allocate (real_no(1:ewld%bspline),inv_no(1:ewld%bspline), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'bspgen allocation failure'
       Call error(0,message)
    End If

  ! construct B-splines

    Do i=1,ewld%bspline
       real_no(i) = Real(i,wp)
       inv_no(i)  = 1.0_wp / real_no(i)
    End Do

    Do i=1,nospl

  ! initializing 2nd order B-spline
  ! for u where (0<u<1) and (1<u<2)

       bspx(1,i)=txx(i)-Aint(txx(i),wp)
       bspy(1,i)=tyy(i)-Aint(tyy(i),wp)
       bspz(1,i)=tzz(i)-Aint(tzz(i),wp)

       bspx(2,i)=1.0_wp-bspx(1,i)
       bspy(2,i)=1.0_wp-bspy(1,i)
       bspz(2,i)=1.0_wp-bspz(1,i)

  ! Now on to calculate order k B-spline values at k
  ! points where (0<u<k)

       rix0=bspx(1,i)
       riy0=bspy(1,i)
       riz0=bspz(1,i)

       bsdx(1,i)= 1.0_wp
       bsdy(1,i)= 1.0_wp
       bsdz(1,i)= 1.0_wp

       bsdx(2,i)=-1.0_wp
       bsdy(2,i)=-1.0_wp
       bsdz(2,i)=-1.0_wp

       Do k=3,ewld%bspline-1 ! Order of B-spline

          bspx(k,i)=0.0_wp
          bspy(k,i)=0.0_wp
          bspz(k,i)=0.0_wp

          k_r   =real_no(k)
          km1_rr=inv_no(k-1)

          Do j=k,2,-1 ! Compute order k B-spline at points {k,k-1,...,1}

             jm1_r=real_no(j-1)

             aaa=rix0+jm1_r
             bbb=riy0+jm1_r
             ccc=riz0+jm1_r

             bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
             bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
             bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

          End Do

          bspx(1,i)=bspx(1,i)*rix0*km1_rr
          bspy(1,i)=bspy(1,i)*riy0*km1_rr
          bspz(1,i)=bspz(1,i)*riz0*km1_rr

       End Do

  ! Now compute B-splines for order ewld%bspline at k points where
  ! (0<u<ewld%bspline)

       k=ewld%bspline

       bspx(k,i)=0.0_wp
       bspy(k,i)=0.0_wp
       bspz(k,i)=0.0_wp

       k_r   =real_no(k)
       km1_rr=inv_no(k-1)

       Do j=k,2,-1

  ! Derivatives of B-splines with order ewld%bspline at k-1 points

          bsdx(j,i)=bspx(j,i)-bspx(j-1,i)
          bsdy(j,i)=bspy(j,i)-bspy(j-1,i)
          bsdz(j,i)=bspz(j,i)-bspz(j-1,i)

          jm1_r=real_no(j-1)

          aaa=rix0+jm1_r
          bbb=riy0+jm1_r
          ccc=riz0+jm1_r

          bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
          bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
          bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr
       End Do

       bsdx(1,i)=bspx(1,i)
       bsdy(1,i)=bspy(1,i)
       bsdz(1,i)=bspz(1,i)

       bspx(1,i)=bspx(1,i)*rix0*km1_rr
       bspy(1,i)=bspy(1,i)*riy0*km1_rr
       bspz(1,i)=bspz(1,i)*riz0*km1_rr

    End Do

    Deallocate (real_no,inv_no, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'bspgen allocation failure'
       Call error(0,message)
    End If

  End Subroutine bspgen

  Subroutine bspgen_mpoles(nospl,xxx,yyy,zzz,bspx,bspy,bspz, &
      bsddx,bsddy,bsddz,n_choose_k,config,ewld,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to calculate B-splines for SPME method for
  ! multipolar interactions
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov april 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( ewald_type ), Intent( In    ) :: ewld
    Integer,                                                      Intent( In    ) :: nospl
    Real( Kind = wp ), Dimension( 1:mxatms ),                     Intent( In    ) :: xxx,yyy,zzz

    Real( Kind = wp ), Dimension( 1:ewld%bspline , 1:mxatms ),           Intent(   Out ) :: bspx,bspy,bspz
    Real( Kind = wp ), Dimension( 0:ewld%bspline , 1:ewld%bspline , 1:mxatms ), Intent(   Out ) :: bsddx,bsddy,bsddz
    Real( Kind = wp ), Dimension(1:,1:) :: n_choose_k
    Type( configuration_type ),                                   Intent( InOut ) :: config
    Type( comms_type),                                            Intent( In    ) :: comm

    Integer           :: fail,i,j,k,m,n,p,r,s
    Real( Kind = wp ) :: aaa,bbb,ccc, rix0,riy0,riz0, jm1_r,k_r,km1_rr
    Real( Kind = wp ) :: tempx,tempy,tempz,pcombr

    Real( Kind = wp ), Dimension( : ), Allocatable :: real_no, inv_no, pmo_no

    Character  ( Len = 256 )  ::  message

    fail=0
    Allocate (real_no(1:ewld%bspline),inv_no(1:ewld%bspline),pmo_no(0:ewld%bspline), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'bspgen_mpoles allocation failure'
       Call error(0,message)
    End If

  ! initialize derivatives

    bsddx(:,:,:)=0.0_wp
    bsddy(:,:,:)=0.0_wp
    bsddz(:,:,:)=0.0_wp

  ! construct B-splines

    pmo_no(0) = 1.0_wp
    Do i=1,ewld%bspline
       real_no(i) = Real(i,wp)
       inv_no(i)  = 1.0_wp / real_no(i)
       pmo_no(i)  = Real(-1**i,wp)
    End Do

    Do i=1,nospl

  ! initializing 2nd order B-spline
  ! for u where (0<u<1) and (1<u<2)

       bspx(1,i)=config%parts(i)%xxx-Aint(config%parts(i)%xxx,wp)
       bspy(1,i)=config%parts(i)%yyy-Aint(config%parts(i)%yyy,wp)
       bspz(1,i)=config%parts(i)%zzz-Aint(config%parts(i)%zzz,wp)

       bspx(2,i)=1.0_wp-bspx(1,i)
       bspy(2,i)=1.0_wp-bspy(1,i)
       bspz(2,i)=1.0_wp-bspz(1,i)

  ! compute the (ewld%bspline-2)nd derivatives

       k=2; p=ewld%bspline-2

       Do j=1,ewld%bspline

          m=Max(0,j-k)
          n=Min(p,j-1)

          tempx=0.0_wp ; tempy=0.0_wp ; tempz=0.0_wp

          Do r=m,n
             s     = j - r
             pcombr= pmo_no(r)*n_choose_k(p,r)

             tempx = tempx + pcombr*bspx(s,i)
             tempy = tempy + pcombr*bspy(s,i)
             tempz = tempz + pcombr*bspz(s,i)
          End Do

          bsddx(p,j,i)=tempx
          bsddy(p,j,i)=tempy
          bsddz(p,j,i)=tempz

       End Do

  ! Now on to calculate order k B-spline values at k
  ! points where (0<u<k)

       rix0=bspx(1,i)
       riy0=bspy(1,i)
       riz0=bspz(1,i)

       Do k=3,ewld%bspline-1 ! Order of B-spline

          bspx(k,i)=0.0_wp
          bspy(k,i)=0.0_wp
          bspz(k,i)=0.0_wp

          k_r   =real_no(k)
          km1_rr=inv_no(k-1)

          Do j=k,2,-1 ! Compute order k B-spline at points {k,k-1,...,1}

             jm1_r=real_no(j-1)

             aaa=rix0+jm1_r
             bbb=riy0+jm1_r
             ccc=riz0+jm1_r

             bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
             bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
             bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

          End Do

          bspx(1,i)=bspx(1,i)*rix0*km1_rr
          bspy(1,i)=bspy(1,i)*riy0*km1_rr
          bspz(1,i)=bspz(1,i)*riz0*km1_rr

  ! compute the (ewld%bspline-3)rd to 1st derivatives

          p = ewld%bspline-k

          Do j=1,ewld%bspline

             m=Max(0,j-k)
             n=Min(p,j-1)

             tempx=0.0_wp ; tempy=0.0_wp ; tempz=0.0_wp

             Do r=m,n
                s     = j - r
                pcombr= pmo_no(r)*n_choose_k(p,r)

                tempx = tempx + pcombr*bspx(s,i)
                tempy = tempy + pcombr*bspy(s,i)
                tempz = tempz + pcombr*bspz(s,i)
             End Do

             bsddx(p,j,i)=tempx
             bsddy(p,j,i)=tempy
             bsddz(p,j,i)=tempz

          End Do

       End Do

  ! Now compute B-splines for order ewld%bspline at k points where
  ! (0<u<ewld%bspline)

       k=ewld%bspline

       bspx(k,i)=0.0_wp
       bspy(k,i)=0.0_wp
       bspz(k,i)=0.0_wp

       k_r   =real_no(k)
       km1_rr=inv_no(k-1)

       Do j=k,2,-1

  ! B-splines with order ewld%bspline at k-1 points

          jm1_r=real_no(j-1)

          aaa=rix0+jm1_r
          bbb=riy0+jm1_r
          ccc=riz0+jm1_r

          bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
          bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
          bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

       End Do

       bspx(1,i)=bspx(1,i)*rix0*km1_rr
       bspy(1,i)=bspy(1,i)*riy0*km1_rr
       bspz(1,i)=bspz(1,i)*riz0*km1_rr

  ! Now the zeroth derivatives

       bsddx(0,:,i)=bspx(:,i)
       bsddy(0,:,i)=bspy(:,i)
       bsddz(0,:,i)=bspz(:,i)

    End Do

    Deallocate (real_no,inv_no,pmo_no, Stat=fail )
    If (fail > 0) Then
       Write(message,'(a)') 'bspgen_moples deallocation failure'
       Call error(0,message)
    End If

  End Subroutine bspgen_mpoles

  Function Dtpbsp(s1,s2,s3,rcell,bsddx,bsddy,bsddz,n_choose_k,config,ewld)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Function to compute arbitrary derivatives of the product of three
  ! b-splines for use with multipolear interactions
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng april 2014
  ! amended   - i.t.todorov march 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ) :: Dtpbsp

    Type( ewald_type ), Intent( In    ) :: ewld
    Integer,                                 Intent( In   ) :: s1,s2,s3
    Real( Kind = wp ),                       Intent( In   ) :: rcell(9)
    Real( Kind = wp ), Dimension( 0:ewld%bspline ), Intent( In   ) :: bsddx,bsddy,bsddz
    Real( Kind = wp ), Intent( In    ) :: n_choose_k(1:,1:)
    Type( configuration_type ),              Intent( In    ) :: config

    Real( Kind = wp ) :: tx,ty,tz,sx,sy,sz
    Real( Kind = wp ) :: ka11,ka12,ka13,kb21,kb22,kb23,kc31,kc32,kc33
    Integer           :: j1,j2,j3,k1,k2,k3,jj,kk,sk,sk3,sk2,sk1

    Dtpbsp = 0.0_wp

    ka11 = Real(ewld%fft_dim_a,wp)*rcell(1)
    ka12 = Real(ewld%fft_dim_a,wp)*rcell(4)
    ka13 = Real(ewld%fft_dim_a,wp)*rcell(7)
    kb21 = Real(ewld%fft_dim_b,wp)*rcell(2)
    kb22 = Real(ewld%fft_dim_b,wp)*rcell(5)
    kb23 = Real(ewld%fft_dim_b,wp)*rcell(8)
    kc31 = Real(ewld%fft_dim_c,wp)*rcell(3)
    kc32 = Real(ewld%fft_dim_c,wp)*rcell(6)
    kc33 = Real(ewld%fft_dim_c,wp)*rcell(9)

  ! Typically, the box is orthogonal => only diagonals-ka11,kb22,kc33-are non-zero

    If (config%imcon /= 3) Then

       Dtpbsp = ka11**s1 * kb22**s2 * kc33**s3 * bsddx(s1) * bsddy(s2) * bsddz(s3)

    Else

       tz = 1.0_wp
       Do k3 = 0, s3
          ty = tz * n_choose_k(s3,k3); sk3=s3-k3

          Do k2 = 0, s2
             tx = ty * n_choose_k(s2,k2); sk2=s2-k2

             Do k1 = 0, s1
                kk = k1+k2+k3; sk1=s1-k1; sk=sk1+sk2+sk3

                sz = tx * n_choose_k(s1,k1)*bsddx(kk)

                Do j3 = 0, sk3
                   sy = sz * n_choose_k(sk3,j3)*kc33**(sk3-j3)

                   Do j2 = 0, sk2
                      sx = sy * n_choose_k(sk2,j2)*kc32**(sk2-j2)

                      Do j1 = 0, sk1
                         jj = j1+j2+j3

                         Dtpbsp = Dtpbsp + sx * kc31**(sk1-j1) * n_choose_k(sk1,j1)*bsddy(jj)*bsddz(sk-jj)

                         sx=sx*kb21
                      End Do

                      sy=sy*kb22
                   End Do

                   sz=sz*kb23
                End Do

                tx=tx*ka11
             End Do

             ty=ty*ka12
          End Do

          tz=tz*ka13
       End Do

    End If

  End Function Dtpbsp

  Subroutine spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to create complex exponential arrays for b-splines
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith october 1998
  ! amended   - i.t.todorov october 2006
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,              Intent( In    ) :: ndiv1,ndiv2,ndiv3
    Complex( Kind = wp ), Intent(   Out ) :: ww1(1:ndiv1),ww2(1:ndiv2),ww3(1:ndiv3)

    Integer           :: i
    Real( Kind = wp ) :: arg

  ! initialise complex exponential factors

    ww1(1)=(1.0_wp,0.0_wp)

    Do i=1,ndiv1/2
      arg=(twopi/Real(ndiv1,wp))*Real(i,wp)

      ww1(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
      ww1(ndiv1+1-i)=Conjg(ww1(i+1))
    End Do

    ww2(1)=(1.0_wp,0.0_wp)

    Do i=1,ndiv2/2
      arg=(twopi/Real(ndiv2,wp))*Real(i,wp)

      ww2(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
      ww2(ndiv2+1-i)=Conjg(ww2(i+1))
    End Do

    ww3(1)=(1.0_wp,0.0_wp)

    Do i=1,ndiv3/2
      arg=(twopi/Real(ndiv3,wp))*Real(i,wp)

      ww3(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
      ww3(ndiv3+1-i)=Conjg(ww3(i+1))
    End Do

  End Subroutine spl_cexp

End Module ewald
