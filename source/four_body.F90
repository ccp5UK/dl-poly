Module four_body

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global four-body potential variables and
  ! arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov june 2008
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,          Only : wp,wi
  Use comms,          Only : comms_type,gsum,gcheck
  Use domains, Only : domains_type
  Use configuration,  Only : configuration_type
  Use particle, Only : corePart
  Use constants, Only : zero_plus
  Use errors_warnings, Only : error, warning
  Use numerics, Only : invert, dcell
  Use statistics, Only : stats_type
  Use neighbours, Only : neighbours_type
  Implicit None

  Private

  Type, Public :: four_body_type
    Private

    !> Number of four-body potentials
    Integer( Kind = wi ), Public :: n_potential = 0

    Logical, Allocatable, Public :: lfr(:)

    Integer( Kind = wi ), Allocatable, Public :: list(:)
    Integer( Kind = wi ), Allocatable, Public :: ltp(:)

    !> Four-body potential parameters
    Real( Kind = wp ), Allocatable, Public :: param(:,:)
    Real( Kind = wp ), Allocatable, Public :: rct(:)

    !> Maximum number of four-body interactions
    Integer( Kind = wi ), Public :: max_four_body
    !> Maximum number of four-body parameters
    Integer( Kind = wi ), Public :: max_param
    Integer( Kind = wi ), Public :: mx3fbp

    !> Four-body potential cutoff
    Real( Kind = wp ), Public :: cutoff
  Contains
    Private

    Procedure, Public :: init => allocate_four_body_arrays
    Final :: cleanup
  End Type

  Public :: four_body_forces

Contains

  Subroutine allocate_four_body_arrays(T,max_site)
  Class( four_body_type ) :: T
    Integer( Kind = wi ), Intent( In    ) :: max_site

    Integer, Dimension(1:5) :: fail

    fail = 0

    Allocate (T%lfr(1:Merge(max_site,0,T%max_four_body > 0)), stat = fail(1))
    Allocate (T%list(1:T%max_four_body), stat = fail(2))
    Allocate (T%ltp(1:T%max_four_body), stat = fail(3))
    Allocate (T%param(1:T%max_param,1:T%max_four_body), stat = fail(4))
    Allocate (T%rct(1:T%max_four_body), stat = fail(5))

    If (Any(fail > 0)) Call error(1024)

    T%lfr = .false.

    T%list = 0
    T%ltp = 0

    T%param = 0.0_wp
    T%rct = 0.0_wp
  End Subroutine allocate_four_body_arrays

  Subroutine four_body_forces(fourbody,stats,neigh,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating four-body forces (energy and
    ! virial) arising from the inversion angle between three atoms around a
    ! nominated central atom
    !
    ! Note: coordinates are converted to reduced units to avoid a call to
    ! images. The link cell algorithm used here necessitates a
    ! parallelepiped cell geometry.
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1996
    ! amended   - i.t.todorov february 2015
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( four_body_type ),  Intent( InOut ) :: fourbody
    Type( stats_type ),      Intent( InOut ) :: stats
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ),      Intent( InOut ) :: comm

    Logical           :: safe,lx0,lx1,ly0,ly1,lz0,lz1

    Integer           :: fail(1:2),                       &
      i,j,k, ii,jj,kk,ll, ia,ib,ic,id, &
      nbx,nby,nbz, ncells,             &
      ix,iy,iz,icell, jx,jy,jz,jcell,  &
      ifbp,jfbp,kfbp,lfbp,             &
      jklbd,kkfbp,limit, ktyp

    Real( Kind = wp ) :: dispx,dispy,dispz, xdc,ydc,zdc,                     &
      rcell(1:9),celprp(1:10),det,                        &
      xab,yab,zab,rab2,rrab, sxab,syab,szab,              &
      xac,yac,zac,rac2,rrac, sxac,syac,szac,              &
      xad,yad,zad,rad2,rrad, sxad,syad,szad, rbc,rcd,rdb, &
      ubx,uby,ubz,ubn,rub, vbx,vby,vbz,vbn,rvb,wwb,       &
      ucx,ucy,ucz,ucn,ruc, vcx,vcy,vcz,vcn,rvc,wwc,       &
      udx,udy,udz,udn,rud, vdx,vdy,vdz,vdn,rvd,wwd,       &
      cosb,cosc,cosd,                                     &
      rubc,rubd,rucd,rucb,rudb,rudc,                      &
      rvbc,rvbd,rvcd,rvcb,rvdb,rvdc,                      &
      thb,thc,thd, k0,th0,cos0, potfbp,gamb,gamc,gamd,    &
      fax,fay,faz, fbx,fby,fbz,                           &
      fcx,fcy,fcz, fdx,fdy,fdz,                           &
      strs1,strs2,strs3,strs5,strs6,strs9

    ! Number of neighbouring cells to look around for counting
    ! many-body interactions

    Integer, Parameter :: nsbcll = 27

    ! Direction arrays for jumping around in link-cell space

    Integer, Dimension( 1:nsbcll ), Parameter :: &
      nix = (/ 0,  -1,-1,-1, 0, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 0, 1, 1, 1 /) , &
      niy = (/ 0,  -1, 0, 1,-1, 0, 1,-1, 0, 1, -1, 0, 1,-1, 1,-1, 0, 1, -1, 0, 1,-1, 0, 1,-1, 0, 1 /) , &
      niz = (/ 0,  -1,-1,-1,-1,-1,-1,-1,-1,-1,  0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1, 1 /)

    Integer,           Dimension( : ), Allocatable :: link,lct,lst,listin
    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt

    Character( Len = 256 ) :: message

    ! Get the dimensional properties of the MD cell

    Call dcell(config%cell,celprp)

    ! Calculate the number of link-cells per domain in every direction

    nbx=Int(domain%nx_recip*celprp(7)/(fourbody%cutoff+1.0e-6_wp))
    nby=Int(domain%ny_recip*celprp(8)/(fourbody%cutoff+1.0e-6_wp))
    nbz=Int(domain%nz_recip*celprp(9)/(fourbody%cutoff+1.0e-6_wp))

    ! check for link cell algorithm violations

    If (nbx < 3 .or. nby < 3 .or. nbz < 3) Call error(305)

    ncells=(nbx+4)*(nby+4)*(nbz+4)
    If (ncells > neigh%max_cell) Then
      Call warning(90,Real(ncells,wp),Real(neigh%max_cell,wp),2.0_wp)
      neigh%max_cell = Nint(1.25_wp*Real(ncells,wp))
      If (ncells > config%mxatms) Call error(69)
    End If

    fail=0
    Allocate (link(1:config%mxatms),listin(1:config%mxatms),lct(1:ncells),lst(1:ncells), Stat=fail(1))
    Allocate (xxt(1:config%mxatms),yyt(1:config%mxatms),zzt(1:config%mxatms),                   Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'four_body_forces allocation failure'
      Call error(0,message)
    End If

    ! Calculate the displacements from the origin of the MD cell
    ! to the bottom left corner of the left-most link-cell (the halo
    ! has two link-cell width)

    ! First term (0.5_wp) = move to the bottom left corner of MD cell
    ! Second term, first term (side) = scale by the number of domains
    ! in the given direction
    ! Second term, second term, first term (id) = move to the bottom
    ! left corner of this domain in the given direction
    ! Second term, second term, second term (2.0_wp/Real(nb,wp)) =
    ! move to the bottom left corner of the left-most link-cell
    ! (the one that constructs the outer layer of the halo)

    dispx=0.5_wp-domain%nx_recip*(Real(domain%idx,wp)-2.0_wp/Real(nbx,wp))
    dispy=0.5_wp-domain%ny_recip*(Real(domain%idy,wp)-2.0_wp/Real(nby,wp))
    dispz=0.5_wp-domain%nz_recip*(Real(domain%idz,wp)-2.0_wp/Real(nbz,wp))

    ! Get the inverse cell matrix

    Call invert(config%cell,rcell,det)

    ! Convert atomic positions (ALL - halo included) from centred
    ! Cartesian coordinates to reduced space coordinates of
    ! the left-most link-cell

    Do i=1,config%nlast
      If (fourbody%lfr(config%ltype(i))) Then
        xxt(i)=rcell(1)*config%parts(i)%xxx+rcell(4)*config%parts(i)%yyy+&
          rcell(7)*config%parts(i)%zzz+dispx
        yyt(i)=rcell(2)*config%parts(i)%xxx+rcell(5)*config%parts(i)%yyy+&
          rcell(8)*config%parts(i)%zzz+dispy
        zzt(i)=rcell(3)*config%parts(i)%xxx+rcell(6)*config%parts(i)%yyy+&
          rcell(9)*config%parts(i)%zzz+dispz
      End If
    End Do

    ! Form linked neigh%list
    ! Initialise link arrays

    link=0
    lct=0
    lst=0

    ! Get the total number of link-cells in MD cell per direction

    xdc=Real(nbx*domain%nx,wp)
    ydc=Real(nby*domain%ny,wp)
    zdc=Real(nbz*domain%nz,wp)

    ! Move ALL particles in link-cell space:
    ! (0,0,0) left-most link-cell on the domain (double link-cell layer)
    ! (nbx+3,nby+3,nby+3) right-most link-cell on the domain
    !***************************************************************
    ! Note: Due to numerical inaccuracy it is possible that some
    ! domain particles (1,natms) may have link-cell space
    ! coordinates in the inner-most layer of the double link-cell
    ! halo / at least one coordinate as shown (nbx+2,nby+2,nbz+2)^
    ! (1,1,1) / as well as particles (natms+1,nlast) from the inner-
    ! most layer of the double link-cell halo may have link-cell
    ! coordinates in the domain region / all coordinates as shown
    ! (nbx,nby,nbz)^(2,2,2) /.  No EDGE EFFECT exists since all
    ! particles (natms+1,nlast) that are outside the double link-
    ! cell layer are not considered.
    !***************************************************************

    Do i=1,config%nlast
      If (fourbody%lfr(config%ltype(i))) Then

        ! Push cell coordinates accordingly

        If (xxt(i) > -zero_plus) Then
          ix=Int(xdc*xxt(i))
        Else
          ix=Int(xdc*(xxt(i)-1.0_wp))
        End If
        If (yyt(i) > -zero_plus) Then
          iy=Int(ydc*yyt(i))
        Else
          iy=Int(ydc*(yyt(i)-1.0_wp))
        End If
        If (zzt(i) > -zero_plus) Then
          iz=Int(zdc*zzt(i))
        Else
          iz=Int(zdc*(zzt(i)-1.0_wp))
        End If

        ! The story has become more complicated with cutoff padding and the
        ! conditional updates of the VNL and thus the halo as now a domain
        ! (1:natms) particle can enter the halo and vice versa.  So LC
        ! bounding is unsafe!!!
        !
        ! Correction for particles (1,natms) belonging to this domain
        ! (idnode) but due to some tiny numerical inaccuracy kicked into
        ! the halo link-cell space and vice-versa particles from the
        ! halo (natms+1,nlast) kicked into the domain link-cell space
        !
        !        If (i <= natms) Then
        !           If (ix < 2)     ix=2
        !           If (ix > nbx+1) ix=nbx+1
        !
        !           If (iy < 2)     iy=2
        !           If (iy > nby+1) iy=nby+1
        !
        !           If (iz < 2)     iz=2
        !           If (iz > nbz+1) iz=nbz+1
        !        Else
        !           lx0=(ix == 2)
        !           lx1=(ix == nbx)
        !           ly0=(iy == 2)
        !           ly1=(iy == nby)
        !           lz0=(iz == 2)
        !           lz1=(iz == nbz)
        !           If ((lx0 .or. lx1) .and. (ly0 .or. ly1) .and. (lz0 .or. lz1)) Then
        !              If      (lx0) Then
        !                 ix=1
        !              Else If (lx1) Then
        !                 ix=nbx+1
        !              Else If (ly0) Then
        !                 iy=1
        !              Else If (ly1) Then
        !                 iy=nby+1
        !              Else If (lz0) Then
        !                 iz=1
        !              Else If (lz1) Then
        !                 iz=nbz+1
        !              End If
        !           End If
        !        End If

        ! Only for particles onto the domain and in the double link-cell
        ! width layer around it. Discard the rest(-:

        If ( (ix >= 0 .and. ix <= nbx+3) .and. &
          (iy >= 0 .and. iy <= nby+3) .and. &
          (iz >= 0 .and. iz <= nbz+3) ) Then

          ! Hypercube function transformation (counting starts from one
          ! rather than zero /map_domains/ and four more link-cells per
          ! dimension are accounted /coming from the assumption that a full
          ! many-body configuration is in 3x3x3 link-cell cub/)

          icell=1+ix+(nbx+4)*(iy+(nby+4)*iz)

          ! at the end of the do-loop lst will give the length of the chain
          ! for this link-cell

          lst(icell)=lst(icell)+1

          ! link points to the next in chain or zero if end of chain occurs
          ! this is the old lct(icell)

          link(i)=lct(icell)

          ! at the end of the do-loop lct will point to the head of chain
          ! for this link-cell (update of lct(icell))

          lct(icell)=i
        End If

      End If
    End Do

    ! flag for undefined potentials

    safe=.true.

    ! initialise potential energy and virial

    stats%engfbp=0.0_wp
    stats%virfbp=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! loop over central atoms of inversion

    Do ix=1,nbx+2
      Do iy=1,nby+2
        Do iz=1,nbz+2

          ! index of primary cell

          icell=1+ix+(nbx+4)*(iy+(nby+4)*iz)

          ! bypass cell if empty

          If (lct(icell) > 0) Then

            ! Initialise extended head of chain array (for all subcells
            ! around icell and icell itself at a very first instance)
            ! and its length (mini-neigh%list of neighbour cell contents)

            k=0
            listin=0

            ! secondary loop over subcells

            Do kk=1,nsbcll

              jx=ix+nix(kk)
              jy=iy+niy(kk)
              jz=iz+niz(kk)

              jcell=1+jx+(nbx+4)*(jy+(nby+4)*jz)

              ! bypass cell if empty

              If (lct(jcell) > 0) Then
                j=lct(jcell)

                Do ii=1,lst(jcell)
                  k=k+1
                  listin(k)=j
                  j=link(j)
                End Do
              End If

            End Do
            limit=k

            ! loop over the primary cell contents

            Do ii=1,lst(icell)

              ! index of the primary atom (and table type)

              ia=listin(ii)
              ifbp=fourbody%mx3fbp*(config%ltype(ia)-1)

              ! bypass if primary atom type is not involved in interaction

              If (fourbody%list(ifbp+1) >= 0) Then

                ! get all possible permutations of triplets and their indices

                Do jj=1,limit-2
                  ib=listin(jj)

                  Do kk=jj+1,limit-1
                    ic=listin(kk)

                    Do ll=kk+1,limit
                      id=listin(ll)

                      ! (FIRST SHIFT TO LEFT)
                      ! if neither of the secondary atoms coincides with the head one

                      If (ib /= ia .and. ic /= ia .and. id /= ia) Then

                        ! get types of atoms

                        jfbp=Max(config%ltype(ib),config%ltype(ic),config%ltype(id))
                        lfbp=Min(config%ltype(ib),config%ltype(ic),config%ltype(id))

                        ! get index of interaction

                        kfbp=config%ltype(ib)+config%ltype(ic)+config%ltype(id)-jfbp-lfbp
                        jklbd=ifbp+lfbp+(kfbp*(kfbp-1))/2+(jfbp*(jfbp**2-1))/6
                        kkfbp=fourbody%list(jklbd)

                        ! check existence of interaction in system

                        If (kkfbp > 0) Then

                          ! only for non-frozen triplets

                          If (config%lfrzn(ia)*config%lfrzn(ib)*config%lfrzn(ic)*config%lfrzn(id) == 0) Then

                            sxab = xxt(ib)-xxt(ia)
                            syab = yyt(ib)-yyt(ia)
                            szab = zzt(ib)-zzt(ia)

                            xab=config%cell(1)*sxab+config%cell(4)*syab+config%cell(7)*szab
                            If (Abs(xab) < fourbody%cutoff) Then

                              yab=config%cell(2)*sxab+config%cell(5)*syab+config%cell(8)*szab
                              If (Abs(yab) < fourbody%cutoff) Then

                                zab=config%cell(3)*sxab+config%cell(6)*syab+config%cell(9)*szab
                                If (Abs(zab) < fourbody%cutoff) Then

                                  rab2=xab*xab+yab*yab+zab*zab

                                  sxac = xxt(ic)-xxt(ia)
                                  syac = yyt(ic)-yyt(ia)
                                  szac = zzt(ic)-zzt(ia)

                                  xac=config%cell(1)*sxac+config%cell(4)*syac+config%cell(7)*szac
                                  If (Abs(xac) < fourbody%cutoff) Then

                                    yac=config%cell(2)*sxac+config%cell(5)*syac+config%cell(8)*szac
                                    If (Abs(yac) < fourbody%cutoff) Then

                                      zac=config%cell(3)*sxac+config%cell(6)*syac+config%cell(9)*szac
                                      If (Abs(zac) < fourbody%cutoff) Then

                                        rac2=xac*xac+yac*yac+zac*zac

                                        sxad = xxt(id)-xxt(ia)
                                        syad = yyt(id)-yyt(ia)
                                        szad = zzt(id)-zzt(ia)

                                        xad=config%cell(1)*sxad+config%cell(4)*syad+config%cell(7)*szad
                                        If (Abs(xad) < fourbody%cutoff) Then

                                          yad=config%cell(2)*sxad+config%cell(5)*syad+config%cell(8)*szad
                                          If (Abs(yad) < fourbody%cutoff) Then

                                            zad=config%cell(3)*sxad+config%cell(6)*syad+config%cell(9)*szad
                                            If (Abs(zad) < fourbody%cutoff) Then

                                              rad2=xad*xad+yad*yad+zad*zad

                                              ! (SECOND SHIFT TO LEFT)

                                              If (Max(rab2,rac2,rad2) <= fourbody%rct(kkfbp)**2) Then

                                                rrab=1.0_wp/Sqrt(rab2)
                                                rrac=1.0_wp/Sqrt(rac2)
                                                rrad=1.0_wp/Sqrt(rad2)

                                                rbc=xab*xac+yab*yac+zab*zac
                                                rcd=xac*xad+yac*yad+zac*zad
                                                rdb=xad*xab+yad*yab+zad*zab

                                                ! calculate bond-angle-plane vectors

                                                ubx=xac*rrac+xad*rrad
                                                uby=yac*rrac+yad*rrad
                                                ubz=zac*rrac+zad*rrad
                                                ubn=1.0_wp/Sqrt(ubx**2+uby**2+ubz**2)
                                                ubx=ubn*ubx
                                                uby=ubn*uby
                                                ubz=ubn*ubz
                                                rub=xab*ubx+yab*uby+zab*ubz

                                                vbx=xac*rrac-xad*rrad
                                                vby=yac*rrac-yad*rrad
                                                vbz=zac*rrac-zad*rrad
                                                vbn=1.0_wp/Sqrt(vbx**2+vby**2+vbz**2)
                                                vbx=vbn*vbx
                                                vby=vbn*vby
                                                vbz=vbn*vbz
                                                rvb=xab*vbx+yab*vby+zab*vbz
                                                wwb=Sqrt(rub**2+rvb**2)

                                                ucx=xad*rrad+xab*rrab
                                                ucy=yad*rrad+yab*rrab
                                                ucz=zad*rrad+zab*rrab
                                                ucn=1.0_wp/Sqrt(ucx**2+ucy**2+ucz**2)
                                                ucx=ucn*ucx
                                                ucy=ucn*ucy
                                                ucz=ucn*ucz
                                                ruc=xac*ucx+yac*ucy+zac*ucz

                                                vcx=xad*rrad-xab*rrab
                                                vcy=yad*rrad-yab*rrab
                                                vcz=zad*rrad-zab*rrab
                                                vcn=1.0_wp/Sqrt(vcx**2+vcy**2+vcz**2)
                                                vcx=vcn*vcx
                                                vcy=vcn*vcy
                                                vcz=vcn*vcz
                                                rvc=xac*vcx+yac*vcy+zac*vcz
                                                wwc=Sqrt(ruc**2+rvc**2)

                                                udx=xab*rrab+xac*rrac
                                                udy=yab*rrab+yac*rrac
                                                udz=zab*rrab+zac*rrac
                                                udn=1.0_wp/Sqrt(udx**2+udy**2+udz**2)
                                                udx=udn*udx
                                                udy=udn*udy
                                                udz=udn*udz
                                                rud=xad*udx+yad*udy+zad*udz

                                                vdx=xab*rrab-xac*rrac
                                                vdy=yab*rrab-yac*rrac
                                                vdz=zab*rrab-zac*rrac
                                                vdn=1.0_wp/Sqrt(vdx**2+vdy**2+vdz**2)
                                                vdx=vdn*vdx
                                                vdy=vdn*vdy
                                                vdz=vdn*vdz
                                                rvd=xad*vdx+yad*vdy+zad*vdz
                                                wwd=Sqrt(rud**2+rvd**2)

                                                ! calculate inversion angle cosines

                                                cosb=wwb*rrab
                                                cosc=wwc*rrac
                                                cosd=wwd*rrad

                                                ! select potential energy function type

                                                ktyp=fourbody%ltp(kkfbp)

                                                ! calculate potential energy and scalar force term

                                                If      (ktyp == 1) Then

                                                  ! harmonic inversion potential

                                                  k0 =fourbody%param(1,kkfbp)
                                                  th0=fourbody%param(2,kkfbp)

                                                  thb=Acos(cosb)
                                                  thc=Acos(cosc)
                                                  thd=Acos(cosd)

                                                  potfbp=k0*((thb-th0)**2+(thc-th0)**2+(thd-th0)**2)/6.0_wp
                                                  gamb=0.0_wp
                                                  If (Abs(thb) > 1.0e-10_wp) gamb=k0*(thb-th0)/(3.0_wp*Sin(thb))
                                                  gamc=0.0_wp
                                                  If (Abs(thc) > 1.0e-10_wp) gamc=k0*(thc-th0)/(3.0_wp*Sin(thc))
                                                  gamd=0.0_wp
                                                  If (Abs(thd) > 1.0e-10_wp) gamd=k0*(thd-th0)/(3.0_wp*Sin(thd))

                                                Else If (ktyp == 2) Then

                                                  ! harmonic cosine inversion potential

                                                  k0  =fourbody%param(1,kkfbp)
                                                  cos0=fourbody%param(2,kkfbp)

                                                  potfbp=k0*((cosb-cos0)**2+(cosc-cos0)**2+(cosd-cos0)**2)/6.0_wp
                                                  gamb=-k0*(cosb-cos0)/3.0_wp
                                                  gamc=-k0*(cosc-cos0)/3.0_wp
                                                  gamd=-k0*(cosd-cos0)/3.0_wp

                                                Else If (ktyp == 3) Then

                                                  ! planar inversion potentials

                                                  k0=fourbody%param(1,kkfbp)

                                                  potfbp=k0*(1.0_wp-(cosb+cosc+cosd)/3.0_wp)
                                                  gamb=k0/3.0_wp
                                                  gamc=k0/3.0_wp
                                                  gamd=k0/3.0_wp

                                                Else

                                                  ! undefined potential

                                                  safe=.false.
                                                  potfbp=0.0_wp
                                                  gamb=0.0_wp
                                                  gamc=0.0_wp
                                                  gamd=0.0_wp

                                                End If

                                                ! calculate bond and u,v scalar products

                                                rubc=xab*ucx+yab*ucy+zab*ucz
                                                rubd=xab*udx+yab*udy+zab*udz
                                                rucd=xac*udx+yac*udy+zac*udz
                                                rucb=xac*ubx+yac*uby+zac*ubz
                                                rudb=xad*ubx+yad*uby+zad*ubz
                                                rudc=xad*ucx+yad*ucy+zad*ucz

                                                rvbc=xab*vcx+yab*vcy+zab*vcz
                                                rvbd=xab*vdx+yab*vdy+zab*vdz
                                                rvcd=xac*vdx+yac*vdy+zac*vdz
                                                rvcb=xac*vbx+yac*vby+zac*vbz
                                                rvdb=xad*vbx+yad*vby+zad*vbz
                                                rvdc=xad*vcx+yad*vcy+zad*vcz

                                                ! calculate atomic forces

                                                fbx = gamb*(-cosb*xab*rrab**2+rrab*(rub*ubx+rvb*vbx)/wwb) +       &
                                                  ( ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2) -   &
                                                  rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2) ) * &
                                                  gamc*rrac/wwc +                                             &
                                                  ( rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2) +   &
                                                  rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2) ) * &
                                                  gamd*rrad/wwd

                                                fby = gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb) +       &
                                                  ( ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2) -   &
                                                  rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2) ) * &
                                                  gamc*rrac/wwc +                                             &
                                                  ( rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2) +   &
                                                  rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2) ) * &
                                                  gamd*rrad/wwd

                                                fbz = gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb) +       &
                                                  ( ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2) -   &
                                                  rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2) ) * &
                                                  gamc*rrac/wwc +                                             &
                                                  ( rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2) +   &
                                                  rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2) ) * &
                                                  gamd*rrad/wwd

                                                fcx = gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc) +       &
                                                  ( rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2) -   &
                                                  rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2) ) * &
                                                  gamd*rrad/wwd +                                             &
                                                  ( rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2) +   &
                                                  rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2) ) * &
                                                  gamb*rrab/wwb

                                                fcy = gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc) +       &
                                                  ( rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2) -   &
                                                  rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2) ) * &
                                                  gamd*rrad/wwd +                                             &
                                                  ( rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2) +   &
                                                  rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2) ) * &
                                                  gamb*rrab/wwb

                                                fcz = gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc) +       &
                                                  ( rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2) -   &
                                                  rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2) ) * &
                                                  gamd*rrad/wwd +                                             &
                                                  ( rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2) +   &
                                                  rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2) ) * &
                                                  gamb*rrab/wwb

                                                fdx = gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd) +       &
                                                  ( rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2) -   &
                                                  rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2) ) * &
                                                  gamb*rrab/wwb +                                             &
                                                  ( ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2) +   &
                                                  rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2) ) * &
                                                  gamc*rrac/wwc

                                                fdy = gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd) +       &
                                                  ( rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2) -   &
                                                  rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2) ) * &
                                                  gamb*rrab/wwb +                                             &
                                                  ( ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2) +   &
                                                  rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2) ) * &
                                                  gamc*rrac/wwc

                                                fdz = gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd) +       &
                                                  ( rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2) -   &
                                                  rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2) ) * &
                                                  gamb*rrab/wwb +                                             &
                                                  ( ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2) +   &
                                                  rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2) ) * &
                                                  gamc*rrac/wwc

                                                fax = -(fbx+fcx+fdx)
                                                fay = -(fby+fcy+fdy)
                                                faz = -(fbz+fcz+fdz)

                                                If (ia <= config%natms) Then

                                                  ! sum inversion energy

                                                  stats%engfbp=stats%engfbp+potfbp

                                                  ! stress tensor calculation for inversion terms

                                                  strs1 = strs1 + xab*fbx + xac*fcx + xad*fdx
                                                  strs2 = strs2 + yab*fbx + yac*fcx + yad*fdx
                                                  strs3 = strs3 + zab*fbx + zac*fcx + zad*fdx
                                                  strs5 = strs5 + yab*fby + yac*fcy + yad*fdy
                                                  strs6 = strs6 + yab*fbz + yac*fcz + yad*fdz
                                                  strs9 = strs9 + zab*fbz + zac*fcz + zad*fdz

                                                  config%parts(ia)%fxx=config%parts(ia)%fxx+fax
                                                  config%parts(ia)%fyy=config%parts(ia)%fyy+fay
                                                  config%parts(ia)%fzz=config%parts(ia)%fzz+faz

                                                End If

                                                If (ib <= config%natms) Then

                                                  config%parts(ib)%fxx=config%parts(ib)%fxx+fbx
                                                  config%parts(ib)%fyy=config%parts(ib)%fyy+fby
                                                  config%parts(ib)%fzz=config%parts(ib)%fzz+fbz

                                                End If

                                                If (ic <= config%natms) Then

                                                  config%parts(ic)%fxx=config%parts(ic)%fxx+fcx
                                                  config%parts(ic)%fyy=config%parts(ic)%fyy+fcy
                                                  config%parts(ic)%fzz=config%parts(ic)%fzz+fcz

                                                End If

                                                If (id <= config%natms) Then

                                                  config%parts(id)%fxx=config%parts(id)%fxx+fdx
                                                  config%parts(id)%fyy=config%parts(id)%fyy+fdy
                                                  config%parts(id)%fzz=config%parts(id)%fzz+fdz

                                                End If

                                              End If

                                              ! (BACK - SECOND SHIFT TO LEFT)

                                            End If
                                          End If
                                        End If
                                      End If
                                    End If
                                  End If
                                End If
                              End If
                            End If
                          End If
                        End If
                      End If

                      ! (BACK - FIRST SHIFT TO LEFT)

                    End Do
                  End Do
                End Do
              End If
            End Do
          End If
        End Do
      End Do
    End Do

    ! check for undefined potentials

    Call gcheck(comm,safe)
    If (.not.safe) Call error(443)

    ! global sum of four-body potential: virial is zero!!!

    Call gsum(comm,stats%engfbp)

    ! complete stress tensor

    stats%stress(1) = stats%stress(1) + strs1
    stats%stress(2) = stats%stress(2) + strs2
    stats%stress(3) = stats%stress(3) + strs3
    stats%stress(4) = stats%stress(4) + strs2
    stats%stress(5) = stats%stress(5) + strs5
    stats%stress(6) = stats%stress(6) + strs6
    stats%stress(7) = stats%stress(7) + strs3
    stats%stress(8) = stats%stress(8) + strs6
    stats%stress(9) = stats%stress(9) + strs9

    Deallocate (link,listin,lct,lst, Stat=fail(1))
    Deallocate (xxt,yyt,zzt,         Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'four_body_forces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine four_body_forces

  Subroutine cleanup(T)
    Type( four_body_type ) :: T

    If (Allocated(T%lfr)) Then
      Deallocate(T%lfr)
    End If

    If (Allocated(T%list)) Then
      Deallocate(T%list)
    End If
    If (Allocated(T%ltp)) Then
      Deallocate(T%ltp)
    End If

    If (Allocated(T%param)) Then
      Deallocate(T%param)
    End If
    If (Allocated(T%rct)) Then
      Deallocate(T%rct)
    End If
  End Subroutine cleanup
End Module four_body
