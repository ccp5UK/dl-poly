Module ewald

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring ewald routines arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use comms,           Only : ExchgGrid_tag,comms_type,wp_mpi,gsend,gwait, &
                              girecv
  Use setup,           Only : mxatms,nrite,mxspl,mxspl2,twopi,kmaxa,kmaxb,kmaxc
  Use configuration,   Only : natms,fxx,fyy,fzz,imcon
  Use domains,         Only : map
  Use mpole,           Only : ncombk

  Use errors_warnings, Only : error

#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

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
      Procedure :: ewald_allocate_kall_arrays
      Procedure :: ewald_allocate_kfrz_arrays
      Procedure, Public :: check => ewald_check
      Procedure, Public :: refresh => ewald_refresh
      Final :: ewald_deallocate
  End Type ewald_type

  Public :: bspcoe, bspgen, bspgen_mpoles, dtpbsp, spl_cexp, dlpfft3, exchange_grid

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

          Call T%ewald_allocate_kfrz_arrays()

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
            ixt-mxspl2+1, ixt  ,         iyb, iyt,                    izb, izt, &
            ixdb        , ixb-1,         iyb, iyt,                    izb, izt )

  ! +Y direction face (including the +X face extension) - negative halo

       Call exchange_grid_halo( map(3),                     map(4), &
            ixdb, ixt,                   iyt-mxspl2+1, iyt  ,         izb, izt, &
            ixdb, ixt,                   iydb        , iyb-1,         izb, izt )

  ! +Z direction face (including the +Y+X faces extensions) - negative halo

       Call exchange_grid_halo( map(5),                     map(6), &
            ixdb, ixt,                   iydb, iyt,                   izt-mxspl2+1, izt, &
            ixdb, ixt,                   iydb, iyt,                   izdb        , izb-1 )

    Else

  ! +X direction face - negative halo

       Call exchange_grid_halo( map(1),                     map(2), &
            ixt-mxspl2+1, ixt  ,         iyb, iyt,                    izb, izt, &
            ixdb        , ixb-1,         iyb, iyt,                    izb, izt )
  !          (ixt)-(ixt-mxspl2+1)+1=mxspl2
  !          (ixb-1)-(ixdb)+1=mxspl2

  ! -X direction face - positive halo

       Call exchange_grid_halo( map(2),                     map(1), &
            ixb          , ixb+delspl-1, iyb, iyt,                    izb, izt, &
            ixdt-delspl+1, ixdt        , iyb, iyt,                    izb, izt )
  !          (ixb+delspl-1)-(ixb)+1=delspl
  !          (ixdt)-(ixdt-delspl+1)+1=delspl


  ! +Y direction face (including the +&-X faces extensions) - negative halo

       Call exchange_grid_halo( map(3),                     map(4), &
            ixdb, ixdt,                  iyt-mxspl2+1, iyt  ,         izb, izt, &
            ixdb, ixdt,                  iydb        , iyb-1,         izb, izt )

  ! -Y direction face (including the +&-X faces extensions) - positive halo

       Call exchange_grid_halo( map(4),                     map(3), &
            ixdb, ixdt,                  iyb          , iyb+delspl-1, izb, izt, &
            ixdb, ixdt,                  iydt-delspl+1, iydt        , izb, izt )

  ! +Z direction face (including the +&-Y+&-X faces extensions) - negative halo

       Call exchange_grid_halo( map(5),                     map(6), &
            ixdb, ixdt,                  iydb, iydt,                  izt-mxspl2+1, izt, &
            ixdb, ixdt,                  iydb, iydt,                  izdb        , izb-1 )

  ! -Z direction face (including the +&-Y+&-X faces extensions) - positive halo

       Call exchange_grid_halo( map(6),                     map(5), &
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

  Subroutine bspcoe(nospl,kmax1,kmax2,kmax3,csp,bscx,bscy,bscz,ww1,ww2,ww3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to calculate B-spline coefficients for
  ! Euler exponential splines
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith july 1998
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                 Intent( In    ) :: nospl,kmax1,kmax2,kmax3
    Real( Kind = wp ),       Intent(   Out ) :: csp(1:mxspl)
    Complex( Kind = wp ),    Intent( In    ) :: ww1(1:kmaxa),ww2(1:kmaxb),ww3(1:kmaxc)
    Complex( Kind = wp ),    Intent(   Out ) :: bscx(1:kmaxa),bscy(1:kmaxb),bscz(1:kmaxc)

    Integer              :: i,j,k
    Complex( Kind = wp ) :: ccc

  ! calculate B-splines at knots

    csp(1)=0.0_wp
    csp(2)=1.0_wp

    Do k=3,nospl
       csp(k)=0.0_wp

       Do j=k,2,-1
          csp(j)=(Real(j-1,wp)*csp(j)+Real(k-j+1,wp)*csp(j-1))/Real(k-1,wp)
       End Do
    End Do

  ! calculate B-spline coefficients

    Do i=0,kmax1-1
       ccc=(0.0_wp,0.0_wp)

       Do k=0,nospl-2
          ccc=ccc+csp(k+2)*ww1(Mod(i*k,kmax1)+1)
       End Do

       bscx(i+1)=ww1(Mod(i*(nospl-1),kmax1)+1)/ccc
    End Do

    Do i=0,kmax2-1
       ccc=(0.0_wp,0.0_wp)

       Do k=0,nospl-2
          ccc=ccc+csp(k+2)*ww2(Mod(i*k,kmax2)+1)
       End Do

       bscy(i+1)=ww2(Mod(i*(nospl-1),kmax2)+1)/ccc
    End Do

    Do i=0,kmax3-1
       ccc=(0.0_wp,0.0_wp)

       Do k=0,nospl-2
          ccc=ccc+csp(k+2)*ww3(Mod(i*k,kmax3)+1)
       End Do

       bscz(i+1)=ww3(Mod(i*(nospl-1),kmax3)+1)/ccc
    End Do

  End Subroutine bspcoe

  Subroutine bspgen(natms,nospl,xxx,yyy,zzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to calculate B-splines for SPME method
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith july 1998
  ! amended   - i.t.todorov april 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: natms,nospl
    Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: xxx,yyy,zzz

    Real( Kind = wp ), Dimension( 1:mxspl , 1:mxatms ), Intent(   Out ) :: &
                                               bsdx,bsdy,bsdz,bspx,bspy,bspz
    Type( comms_type ),                       Intent( In    ) :: comm

    Integer           :: fail,i,j,k
    Real( Kind = wp ) :: aaa,bbb,ccc, rix0,riy0,riz0, jm1_r,k_r,km1_rr

    Real( Kind = wp ), Dimension( : ), Allocatable :: real_no, inv_no

    Character ( Len = 256 )   ::  message

    fail=0
    Allocate (real_no(1:nospl),inv_no(1:nospl), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'bspgen allocation failure'
       Call error(0,message)
    End If

  ! construct B-splines

    Do i=1,nospl
       real_no(i) = Real(i,wp)
       inv_no(i)  = 1.0_wp / real_no(i)
    End Do

    Do i=1,natms

  ! initializing 2nd order B-spline
  ! for u where (0<u<1) and (1<u<2)

       bspx(1,i)=xxx(i)-Aint(xxx(i),wp)
       bspy(1,i)=yyy(i)-Aint(yyy(i),wp)
       bspz(1,i)=zzz(i)-Aint(zzz(i),wp)

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

       Do k=3,nospl-1 ! Order of B-spline

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

  ! Now compute B-splines for order nospl at k points where
  ! (0<u<nospl)

       k=nospl

       bspx(k,i)=0.0_wp
       bspy(k,i)=0.0_wp
       bspz(k,i)=0.0_wp

       k_r   =real_no(k)
       km1_rr=inv_no(k-1)

       Do j=k,2,-1

  ! Derivatives of B-splines with order nospl at k-1 points

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

  Subroutine bspgen_mpoles(natms,nospl,xxx,yyy,zzz,bspx,bspy,bspz,bsddx,bsddy,bsddz,comm)

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

    Integer,                                                      Intent( In    ) :: natms,nospl
    Real( Kind = wp ), Dimension( 1:mxatms ),                     Intent( In    ) :: xxx,yyy,zzz

    Real( Kind = wp ), Dimension( 1:mxspl , 1:mxatms ),           Intent(   Out ) :: bspx,bspy,bspz
    Real( Kind = wp ), Dimension( 0:mxspl , 1:mxspl , 1:mxatms ), Intent(   Out ) :: bsddx,bsddy,bsddz
    Type( comms_type),                                            Intent( In    ) :: comm

    Integer           :: fail,i,j,k,m,n,p,r,s
    Real( Kind = wp ) :: aaa,bbb,ccc, rix0,riy0,riz0, jm1_r,k_r,km1_rr
    Real( Kind = wp ) :: tempx,tempy,tempz,pcombr

    Real( Kind = wp ), Dimension( : ), Allocatable :: real_no, inv_no, pmo_no

    Character  ( Len = 256 )  ::  message

    fail=0
    Allocate (real_no(1:nospl),inv_no(1:nospl),pmo_no(0:nospl), Stat=fail)
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
    Do i=1,nospl
       real_no(i) = Real(i,wp)
       inv_no(i)  = 1.0_wp / real_no(i)
       pmo_no(i)  = Real(-1**i,wp)
    End Do

    Do i=1,natms

  ! initializing 2nd order B-spline
  ! for u where (0<u<1) and (1<u<2)

       bspx(1,i)=xxx(i)-Aint(xxx(i),wp)
       bspy(1,i)=yyy(i)-Aint(yyy(i),wp)
       bspz(1,i)=zzz(i)-Aint(zzz(i),wp)

       bspx(2,i)=1.0_wp-bspx(1,i)
       bspy(2,i)=1.0_wp-bspy(1,i)
       bspz(2,i)=1.0_wp-bspz(1,i)

  ! compute the (nospl-2)nd derivatives

       k=2; p=nospl-2

       Do j=1,nospl

          m=Max(0,j-k)
          n=Min(p,j-1)

          tempx=0.0_wp ; tempy=0.0_wp ; tempz=0.0_wp

          Do r=m,n
             s     = j - r
             pcombr= pmo_no(r)*ncombk(p,r)

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

       Do k=3,nospl-1 ! Order of B-spline

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

  ! compute the (nospl-3)rd to 1st derivatives

          p = nospl-k

          Do j=1,nospl

             m=Max(0,j-k)
             n=Min(p,j-1)

             tempx=0.0_wp ; tempy=0.0_wp ; tempz=0.0_wp

             Do r=m,n
                s     = j - r
                pcombr= pmo_no(r)*ncombk(p,r)

                tempx = tempx + pcombr*bspx(s,i)
                tempy = tempy + pcombr*bspy(s,i)
                tempz = tempz + pcombr*bspz(s,i)
             End Do

             bsddx(p,j,i)=tempx
             bsddy(p,j,i)=tempy
             bsddz(p,j,i)=tempz

          End Do

       End Do

  ! Now compute B-splines for order nospl at k points where
  ! (0<u<nospl)

       k=nospl

       bspx(k,i)=0.0_wp
       bspy(k,i)=0.0_wp
       bspz(k,i)=0.0_wp

       k_r   =real_no(k)
       km1_rr=inv_no(k-1)

       Do j=k,2,-1

  ! B-splines with order nospl at k-1 points

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

  Function Dtpbsp(s1,s2,s3,rcell,bsddx,bsddy,bsddz)

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

    Integer,                                 Intent( In   ) :: s1,s2,s3
    Real( Kind = wp ),                       Intent( In   ) :: rcell(9)
    Real( Kind = wp ), Dimension( 0:mxspl ), Intent( In   ) :: bsddx,bsddy,bsddz

    Real( Kind = wp ) :: tx,ty,tz,sx,sy,sz
    Real( Kind = wp ) :: ka11,ka12,ka13,kb21,kb22,kb23,kc31,kc32,kc33
    Integer           :: j1,j2,j3,k1,k2,k3,jj,kk,sk,sk3,sk2,sk1

    Dtpbsp = 0.0_wp

    ka11 = Real(kmaxa,wp)*rcell(1)
    ka12 = Real(kmaxa,wp)*rcell(4)
    ka13 = Real(kmaxa,wp)*rcell(7)
    kb21 = Real(kmaxb,wp)*rcell(2)
    kb22 = Real(kmaxb,wp)*rcell(5)
    kb23 = Real(kmaxb,wp)*rcell(8)
    kc31 = Real(kmaxc,wp)*rcell(3)
    kc32 = Real(kmaxc,wp)*rcell(6)
    kc33 = Real(kmaxc,wp)*rcell(9)

  ! Typically, the box is orthogonal => only diagonals-ka11,kb22,kc33-are non-zero

    If (imcon /= 3) Then

       Dtpbsp = ka11**s1 * kb22**s2 * kc33**s3 * bsddx(s1) * bsddy(s2) * bsddz(s3)

    Else

       tz = 1.0_wp
       Do k3 = 0, s3
          ty = tz * ncombk(s3,k3); sk3=s3-k3

          Do k2 = 0, s2
             tx = ty * ncombk(s2,k2); sk2=s2-k2

             Do k1 = 0, s1
                kk = k1+k2+k3; sk1=s1-k1; sk=sk1+sk2+sk3

                sz = tx * ncombk(s1,k1)*bsddx(kk)

                Do j3 = 0, sk3
                   sy = sz * ncombk(sk3,j3)*kc33**(sk3-j3)

                   Do j2 = 0, sk2
                      sx = sy * ncombk(sk2,j2)*kc32**(sk2-j2)

                      Do j1 = 0, sk1
                         jj = j1+j2+j3

                         Dtpbsp = Dtpbsp + sx * kc31**(sk1-j1) * ncombk(sk1,j1)*bsddy(jj)*bsddz(sk-jj)

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

  Subroutine dlpfft3(isw,ndiv1,ndiv2,ndiv3,ww1,ww2,ww3,aaa,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 3D fast fourier transform routine (dependent upon spl_cexp)
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith july 1998
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                 Intent( In    ) :: isw,ndiv1,ndiv2,ndiv3
    Complex( Kind = wp ),    Intent( InOut ) :: aaa(1:ndiv1,1:ndiv2,1:ndiv3), &
                                                ww1(1:ndiv1),ww2(1:ndiv2),ww3(1:ndiv3)
    Type( comms_type ),      Intent( In    ) :: comm

    Logical, Save :: newjob = .true.
    Integer, Save :: nu1,nu2,nu3

    Integer              :: fail,iii,jjj,kkk,i,j,k,l,jj2,kk1,k12,num
    Real( Kind = wp )    :: tmp
    Complex( Kind = wp ) :: ttt

    Integer, Dimension( : ), Allocatable, Save :: key1,key2,key3

    Character ( Len = 256 ) :: message   

    If (newjob) Then
       newjob = .false.

  ! check FFT array dimensions

       nu1 = Nint(Log(Real(ndiv1,wp))/Log(2.0_wp))
       nu2 = Nint(Log(Real(ndiv2,wp))/Log(2.0_wp))
       nu3 = Nint(Log(Real(ndiv3,wp))/Log(2.0_wp))

       tmp = Log(Real(ndiv1,wp))*Log(Real(ndiv2,wp))*Log(Real(ndiv3,wp)) / &
             (Log(2.0_wp)**3)

       If (Abs(tmp-Real(nu1*nu2*nu3,wp)) > 1.0e-6_wp) Then
          Write(message,'(a,20i2)') 'error - FFT array dimension not 2^N ',ndiv1,ndiv2,ndiv3
          Call error(0,message)
       End If

  ! allocate reverse bit address arrays

       Allocate (key1(1:ndiv1),key2(1:ndiv2),key3(1:ndiv3), Stat=fail)
       If (fail > 0) Then
          Write(message,'(a)') 'dlpfft3 (FFT reverse bit address arrays) allocation failure'
          Call error(0,message)
       End If

  ! set reverse bit address arrays

       Do kkk=1,ndiv1
          iii=0
          jjj=kkk-1

          Do j=1,nu1
             jj2=jjj/2
             iii=2*(iii-jj2)+jjj
             jjj=jj2
          End Do

          key1(kkk)=iii+1
       End Do

       Do kkk=1,ndiv2
          iii=0
          jjj=kkk-1

          Do j=1,nu2
             jj2=jjj/2
             iii=2*(iii-jj2)+jjj
             jjj=jj2
          End Do

          key2(kkk)=iii+1
       End Do

       Do kkk=1,ndiv3
          iii=0
          jjj=kkk-1

          Do j=1,nu3
             jj2=jjj/2
             iii=2*(iii-jj2)+jjj
             jjj=jj2
          End Do

          key3(kkk)=iii+1
       End Do

  ! initialise complex exponential factors
  ! Call spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

    End If

  ! take conjugate of exponentials if required

    If (isw < 0) Then
       Do i=1,ndiv1
          ww1(i)=Conjg(ww1(i))
       End Do

       Do i=1,ndiv2
          ww2(i)=Conjg(ww2(i))
       End Do

       Do i=1,ndiv3
          ww3(i)=Conjg(ww3(i))
       End Do
    End If

  ! perform fourier transform in X direction

    kkk=0
    num=ndiv1/2
    Do l=1,nu1
       Do While (kkk < ndiv1)
          Do i=1,num
             iii=key1(kkk/num+1)
             kk1=kkk+1
             k12=kk1+num

             Do j=1,ndiv2
                Do k=1,ndiv3
                   ttt=aaa(k12,j,k)*ww1(iii)
                   aaa(k12,j,k)=aaa(kk1,j,k)-ttt
                   aaa(kk1,j,k)=aaa(kk1,j,k)+ttt
                End Do
             End Do

             kkk=kkk+1
          End Do

          kkk=kkk+num
       End Do

       kkk=0
       num=num/2
    End Do

  ! unscramble the fft using bit address array

    Do kkk=1,ndiv1
       iii=key1(kkk)

       If (iii > kkk) Then
          Do j=1,ndiv2
             Do k=1,ndiv3
                ttt=aaa(kkk,j,k)
                aaa(kkk,j,k)=aaa(iii,j,k)
                aaa(iii,j,k)=ttt
             End Do
          End Do
       End If
    End Do

  ! perform fourier transform in Y direction

    kkk=0
    num=ndiv2/2
    Do l=1,nu2
       Do While (kkk < ndiv2)
          Do i=1,num
             iii=key2(kkk/num+1)
             kk1=kkk+1
             k12=kk1+num

             Do j=1,ndiv1
                Do k=1,ndiv3
                   ttt=aaa(j,k12,k)*ww2(iii)
                   aaa(j,k12,k)=aaa(j,kk1,k)-ttt
                   aaa(j,kk1,k)=aaa(j,kk1,k)+ttt
                End Do
             End Do

             kkk=kkk+1
          End Do

          kkk=kkk+num
       End Do

       kkk=0
       num=num/2
    End Do

  ! unscramble the fft using bit address array

    Do kkk=1,ndiv2
       iii=key2(kkk)

       If (iii > kkk) Then
          Do j=1,ndiv1
             Do k=1,ndiv3
                ttt=aaa(j,kkk,k)
                aaa(j,kkk,k)=aaa(j,iii,k)
                aaa(j,iii,k)=ttt
             End Do
          End Do
       End If
    End Do

  ! perform fourier transform in Z direction

     kkk=0
     num=ndiv3/2
     Do l=1,nu3

        Do While (kkk < ndiv3)
           Do i=1,num
              iii=key3(kkk/num+1)
              kk1=kkk+1
              k12=kk1+num

              Do j=1,ndiv1
                 Do k=1,ndiv2
                    ttt=aaa(j,k,k12)*ww3(iii)
                    aaa(j,k,k12)=aaa(j,k,kk1)-ttt
                    aaa(j,k,kk1)=aaa(j,k,kk1)+ttt
                 End Do
              End Do

              kkk=kkk+1
           End Do

           kkk=kkk+num
        End Do

        kkk=0
        num=num/2
     End Do

  ! unscramble the fft using bit address array

     Do kkk=1,ndiv3
        iii=key3(kkk)

        If (iii > kkk) Then
           Do j=1,ndiv1
              Do k=1,ndiv2
                 ttt=aaa(j,k,kkk)
                 aaa(j,k,kkk)=aaa(j,k,iii)
                 aaa(j,k,iii)=ttt
             End Do
          End Do
       End If
    End Do

  ! restore exponentials to unconjugated values if necessary

    If (isw < 0) Then
       Do i=1,ndiv1
          ww1(i)=Conjg(ww1(i))
       End Do

       Do i=1,ndiv2
          ww2(i)=Conjg(ww2(i))
       End Do

       Do i=1,ndiv3
          ww3(i)=Conjg(ww3(i))
       End Do
    End If

  End Subroutine dlpfft3
End Module ewald
