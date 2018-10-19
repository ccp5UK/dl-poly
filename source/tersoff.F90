Module tersoff

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global tersoff interaction variables and
  ! arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov september 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use comms,   Only : comms_type,gsum
  Use constants, Only : zero_plus, pi
  Use domains, Only : domains_type
  Use configuration,  Only : configuration_type
  Use particle,       Only : corePart
  Use errors_warnings, Only : error,warning
  Use numerics, Only : dcell, invert
  Use statistics, Only : stats_type
  Use neighbours, Only : neighbours_type
  Implicit None

  Private

  !> Type to hold Tersoff potential data
  Type, Public :: tersoff_type
    Private

    !> Number of Tersoff potentials
    Integer( Kind = wi ), Public :: n_potential
    !> Type of potential
    !>
    !> - 1 Tersoff
    !> - 2 Kumagi-Izumi-Hara-Sakai
    Integer( Kind = wi ), Public :: key_pot

    !> Global Tersoff potential cutoff
    Real( Kind = wp ), Public :: cutoff

    !> Tersoff potential switch per atom?
    Logical, Allocatable, Public :: lfr(:)

    Integer( Kind = wi ), Allocatable, Public :: list(:)
    Integer( Kind = wi ), Allocatable, Public :: ltp(:)

    !> Tersoff potential parameters
    Real( Kind = wp ), Allocatable, Public :: param(:,:)
    Real( Kind = wp ), Allocatable, Public :: param2(:,:)

    !> Tersoff screening function? elements(-1,*,1) were cutoffs...
    Real( Kind = wp ), Allocatable :: vmbp(:,:,:)
    Real( Kind = wp ), Allocatable :: gmbp(:,:,:)

    !> Potential cutoff for particular interactions, replaces tersoffs%vmbp(-1,*,1)
    Real( Kind = wp ), Allocatable :: cut(:)

    !> Maximum number of Tersoff interactions
    Integer( Kind = wi ), Public :: max_ter
    !> Maximum number of Tersoff paramters
    Integer( Kind = wi ), Public :: max_param
    !> Maximum number of grid points
    Integer( Kind = wi ), Public :: max_grid
  Contains
    Private

    Procedure, Public :: init => allocate_tersoff_arrays
    Final :: cleanup
  End Type tersoff_type

  Public :: tersoff_forces,tersoff_generate

Contains

  Subroutine allocate_tersoff_arrays(T,max_site)
  Class( tersoff_type ) :: T
    Integer( Kind = wi ), Intent( In    ) :: max_site

    Integer( Kind = wi ) :: nprter
    Integer, Dimension(8) :: fail

    nprter = (T%max_ter*(T%max_ter+1))/2

    fail = 0

    Allocate (T%lfr(1:Merge(max_site,0,T%max_ter > 0)), stat=fail(1))
    Allocate (T%list(1:T%max_ter), stat=fail(2))
    Allocate (T%ltp(1:T%max_ter), stat=fail(3))
    Allocate (T%param(1:T%max_param,1:T%max_ter), stat=fail(4))
    If (T%key_pot == 1) Then
      Allocate (T%param2(1:nprter,1:2), stat=fail(5))
    End If
    Allocate (T%vmbp(0:T%max_grid,1:nprter,1:3), stat=fail(6))
    Allocate (T%gmbp(0:T%max_grid,1:nprter,1:3), stat=fail(7))
    Allocate (T%cut(1:nprter), stat=fail(8))

    If (Any(fail > 0)) Call error(1027)

    T%lfr = .false.

    T%list = 0
    T%ltp = 0

    T%param = 0.0_wp
    If (T%key_pot == 1) Then
      T%param2 = 0.0_wp
    End If
    T%vmbp = 0.0_wp
    T%gmbp = 0.0_wp
    T%cut = 0.0_wp

  End Subroutine allocate_tersoff_arrays

  Subroutine tersoff_forces(tersoffs,stats,neigh,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating tersoff forces (energy and
    ! virial) arising from the tersoff potential as defined by
    !
    ! 1) J. Tersoff, Phys. Rev. B 39 (1989) 5566
    !
    ! 2) T. Kumagai, S. Izumi, S. Hara, and S. Sakai, Computational
    !    Materials Science 39, 457 (2007), ISSN 09270256
    !
    ! Note: coordinates are converted to reduced units to avoid a call to
    ! images.  The link cell algorithm used here necessitates a
    ! parallelepiped cell geometry.
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith october 2004
    ! amended   - i.t.todorov march 2016
    ! amended   - g.khara september 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( tersoff_type ), Intent( InOut )  :: tersoffs
    Type( stats_type ), Intent( InOut )  :: stats
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ),   Intent( InOut )  :: config
    Type( comms_type ), Intent( InOut )  :: comm

    ! flag for undefined potentials NOT NEEDED HERE YET
    Logical           :: lx0,lx1,ly0,ly1,lz0,lz1,flag3 !,safe

    Integer           :: fail(1:7),                      &
      i,j,k, ii,jj,kk,ll,             &
      nbx,nby,nbz, ncells,            &
      ix,iy,iz,icell, jx,jy,jz,jcell, &
      iatm,jatm,katm, iter,jter,kter, &
      ijter,ikter,limit

    Real( Kind = wp ) :: dispx,dispy,dispz, xdc,ydc,zdc,                 &
      rdr,rcell(1:9),celprp(1:10),det,                &
      sxij,syij,szij,                                 &
      gk0,gk1,gk2,vk0,vk1,vk2,                        &
      t1,t2,ppp,bi,ei,ci,di,c1i,c2i,c3i,c4i,c5i,hi,   &
      ak,bk,gtheta,cost,hmct2,c4exp,                  &
      eterm,vterm,gterm,gamma,gam_ij,                 &
      gam_dg,gam_df,gam_dw, fxj,fyj,fzj, fxk,fyk,fzk, &
      strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

    ! Number of neighbouring cells to look around for counting tersoff
    ! interactions

    Integer, Parameter :: nsbcll = 27

    ! Direction arrays for jumping around in link-cell space

    Integer, Dimension( 1:nsbcll ), Parameter :: &
      nix = (/ 0,  -1,-1,-1, 0, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 0, 1, 1, 1 /) , &
      niy = (/ 0,  -1, 0, 1,-1, 0, 1,-1, 0, 1, -1, 0, 1,-1, 1,-1, 0, 1, -1, 0, 1,-1, 0, 1,-1, 0, 1 /) , &
      niz = (/ 0,  -1,-1,-1,-1,-1,-1,-1,-1,-1,  0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1, 1 /)

    Integer,           Dimension( : ), Allocatable :: link,lct,lst,listin
    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt
    Real( Kind = wp ), Dimension( : ), Allocatable :: xtf,ytf,ztf,rtf
    Real( Kind = wp ), Dimension( : ), Allocatable :: ert,eat,grt,gat
    Real( Kind = wp ), Dimension( : ), Allocatable :: scr,gcr,gam,gvr,cst,rkj,wkj


    Character( Len = 256 ) :: message
    ! Get reciprocal of interpolation interval

    rdr=Real(tersoffs%max_grid-4,wp)/tersoffs%cutoff

    ! Get the dimensional properties of the MD cell

    Call dcell(config%cell,celprp)

    ! Calculate the number of link-cells per domain in every direction

    nbx=Int(domain%nx_recip*celprp(7)/(tersoffs%cutoff+1.0e-6_wp))
    nby=Int(domain%ny_recip*celprp(8)/(tersoffs%cutoff+1.0e-6_wp))
    nbz=Int(domain%nz_recip*celprp(9)/(tersoffs%cutoff+1.0e-6_wp))

    ! check for link cell algorithm violations

    If (nbx < 3 .or. nby < 3 .or. nbz < 3) Call error(305)

    ncells=(nbx+4)*(nby+4)*(nbz+4)
    If (ncells > neigh%max_cell) Then
      Call warning(90,Real(ncells,wp),Real(neigh%max_cell,wp),3.0_wp)
      neigh%max_cell = Nint(1.25_wp*Real(ncells,wp))
      If (ncells > config%mxatms) Call error(69)
    End If

    fail=0
    Allocate (link(1:config%mxatms),listin(1:config%mxatms),lct(1:ncells),lst(1:ncells), Stat=fail(1))
    Allocate (xxt(1:config%mxatms),yyt(1:config%mxatms),zzt(1:config%mxatms),                   Stat=fail(2))
    Allocate (xtf(1:neigh%max_list),ytf(1:neigh%max_list),ztf(1:neigh%max_list),rtf(1:neigh%max_list),     Stat=fail(3))
    Allocate (ert(1:neigh%max_list),eat(1:neigh%max_list),grt(1:neigh%max_list),gat(1:neigh%max_list),     Stat=fail(4))
    Allocate (scr(1:neigh%max_list),gcr(1:neigh%max_list),                                 Stat=fail(5))
    Allocate (cst(1:neigh%max_list),gam(1:neigh%max_list),gvr(1:neigh%max_list),                   Stat=fail(6))
    If (tersoffs%key_pot == 2) Allocate (rkj(1:neigh%max_list),wkj(1:neigh%max_list),                Stat=fail(7))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'tersoff_forces allocation failure'
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
      If (tersoffs%lfr(config%ltype(i))) Then
        xxt(i)=rcell(1)*config%parts(i)%xxx+rcell(4)*config%parts(i)%yyy+rcell(7)*&
          config%parts(i)%zzz+dispx
        yyt(i)=rcell(2)*config%parts(i)%xxx+rcell(5)*config%parts(i)%yyy+rcell(8)*&
          config%parts(i)%zzz+dispy
        zzt(i)=rcell(3)*config%parts(i)%xxx+rcell(6)*config%parts(i)%yyy+rcell(9)*&
          config%parts(i)%zzz+dispz
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
      If (tersoffs%lfr(config%ltype(i))) Then

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

        ! Correction for particles (1,natms) belonging to this domain
        ! (idnode) but due to some tiny numerical inaccuracy kicked into
        ! the halo link-cell space and vice-versa particles from the
        ! halo (natms+1,nlast) kicked into the domain link-cell space

        If (i <= config%natms) Then
          If (ix < 2)     ix=2
          If (ix > nbx+1) ix=nbx+1

          If (iy < 2)     iy=2
          If (iy > nby+1) iy=nby+1

          If (iz < 2)     iz=2
          If (iz > nbz+1) iz=nbz+1
        Else
          lx0=(ix == 2)
          lx1=(ix == nbx)
          ly0=(iy == 2)
          ly1=(iy == nby)
          lz0=(iz == 2)
          lz1=(iz == nbz)
          If ((lx0 .or. lx1) .and. (ly0 .or. ly1) .and. (lz0 .or. lz1)) Then
            If      (lx0) Then
              ix=1
            Else If (lx1) Then
              ix=nbx+1
            Else If (ly0) Then
              iy=1
            Else If (ly1) Then
              iy=nby+1
            Else If (lz0) Then
              iz=1
            Else If (lz1) Then
              iz=nbz+1
            End If
          End If
        End If

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

    ! flag for undefined potentials NOT NEEDED HERE YET
    !  safe=.true.

    ! initialise potential energy and virial

    stats%engter=0.0_wp
    stats%virter=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! loop over central atoms of angles

    Do ix=1,nbx+2
      Do iy=1,nby+2
        Do iz=1,nbz+2

          ! index of primary cell

          icell=1+ix+(nbx+4)*(iy+(nby+4)*iz)

          ! bypass cell if empty

          If (lct(icell) > 0) Then

            ! Initialise extended head of chain array (for all subcells
            ! around icell and icell itself at the very first instance)
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

              iatm=listin(ii)
              iter=tersoffs%list(config%ltype(iatm))

              ! bypass if primary atom type is not involved in interaction

              If (iter > 0) Then

                ! construct working arrays by using interpolation

                Do jj=1,limit

                  ! initialise working arrays

                  xtf(jj)=0.0_wp
                  ytf(jj)=0.0_wp
                  ztf(jj)=0.0_wp
                  rtf(jj)=0.0_wp

                  ert(jj)=0.0_wp
                  eat(jj)=0.0_wp
                  grt(jj)=0.0_wp

                  gat(jj)=0.0_wp
                  scr(jj)=0.0_wp
                  gcr(jj)=0.0_wp

                  ! index of the secondary atom

                  jatm=listin(jj)
                  jter=tersoffs%list(config%ltype(jatm))

                  ! bypass if secondary atom type is not involved in interaction

                  If (jter > 0 .and. jatm /= iatm) Then

                    sxij = xxt(jatm)-xxt(iatm)
                    syij = yyt(jatm)-yyt(iatm)
                    szij = zzt(jatm)-zzt(iatm)

                    xtf(jj)=config%cell(1)*sxij+config%cell(4)*syij+config%cell(7)*szij
                    ytf(jj)=config%cell(2)*sxij+config%cell(5)*syij+config%cell(8)*szij
                    ztf(jj)=config%cell(3)*sxij+config%cell(6)*syij+config%cell(9)*szij

                    rtf(jj)=Sqrt(xtf(jj)**2+ytf(jj)**2+ztf(jj)**2)

                    xtf(jj)=xtf(jj)/rtf(jj)
                    ytf(jj)=ytf(jj)/rtf(jj)
                    ztf(jj)=ztf(jj)/rtf(jj)

                    ! interaction index

                    ijter=(Max(iter,jter)*(Max(iter,jter)-1))/2+Min(iter,jter)

                    ! if pair is within cutoff

                    If (rtf(jj) <= tersoffs%cut(ijter)) Then

                      ! if pair is not frozen

                      If (config%lfrzn(iatm)*config%lfrzn(jatm) == 0) Then

                        ll  = Int(rdr*rtf(jj))
                        ppp = rtf(jj)*rdr - Real(ll,wp)

                        ! interpolate screening function

                        vk0 = tersoffs%vmbp(ll,  ijter,1)
                        vk1 = tersoffs%vmbp(ll+1,ijter,1)
                        vk2 = tersoffs%vmbp(ll+2,ijter,1)

                        t1 = vk0 + (vk1 - vk0)*ppp
                        t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                        scr(jj) = t1 + (t2-t1)*ppp*0.5_wp

                        ! interpolate derivative of screening function

                        gk0 = tersoffs%gmbp(ll,  ijter,1)
                        gk1 = tersoffs%gmbp(ll+1,ijter,1)
                        gk2 = tersoffs%gmbp(ll+2,ijter,1)

                        t1 = gk0 + (gk1 - gk0)*ppp
                        t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                        gcr(jj) = -(t1 + (t2 - t1)*ppp*0.5_wp) / rtf(jj)

                        ! interpolate repulsive component of energy

                        vk0 = tersoffs%vmbp(ll,  ijter,2)
                        vk1 = tersoffs%vmbp(ll+1,ijter,2)
                        vk2 = tersoffs%vmbp(ll+2,ijter,2)

                        t1 = vk0 + (vk1 - vk0)*ppp
                        t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                        ert(jj) = t1 + (t2 - t1)*ppp*0.5_wp

                        ! interpolate derivative of repulsive function

                        gk0 = tersoffs%gmbp(ll,  ijter,2)
                        gk1 = tersoffs%gmbp(ll+1,ijter,2)
                        gk2 = tersoffs%gmbp(ll+2,ijter,2)

                        t1 = gk0 + (gk1 - gk0)*ppp
                        t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                        grt(jj) = -(t1 + (t2 - t1)*ppp*0.5_wp) / rtf(jj)

                        ! interpolate attractive component of energy

                        vk0 = tersoffs%vmbp(ll,  ijter,3)
                        vk1 = tersoffs%vmbp(ll+1,ijter,3)
                        vk2 = tersoffs%vmbp(ll+2,ijter,3)

                        t1 = vk0 + (vk1 - vk0)*ppp
                        t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                        eat(jj) = t1 + (t2 - t1)*ppp*0.5_wp

                        ! interpolate derivative of attractive function

                        gk0 = tersoffs%gmbp(ll,  ijter,3)
                        gk1 = tersoffs%gmbp(ll+1,ijter,3)
                        gk2 = tersoffs%gmbp(ll+2,ijter,3)

                        t1 = gk0 + (gk1 - gk0)*ppp
                        t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                        gat(jj) = -(t1 + (t2 - t1)*ppp*0.5_wp) / rtf(jj)

                      End If

                    End If

                  End If

                End Do

                ! calculate three body (attractive) terms

                ! Get parameters for iatm

                If      (tersoffs%key_pot == 1) Then ! TERS
                  bi=tersoffs%param(7, iter)
                  ei=tersoffs%param(8, iter)
                  ci=tersoffs%param(9, iter)
                  di=tersoffs%param(10,iter)
                  hi=tersoffs%param(11,iter)
                Else If (tersoffs%key_pot == 2) Then ! KIHS
                  ei =tersoffs%param(7, iter)
                  di =tersoffs%param(8, iter)
                  c1i=tersoffs%param(9, iter)
                  c2i=tersoffs%param(10,iter)
                  c3i=tersoffs%param(11,iter)
                  c4i=tersoffs%param(12,iter)
                  c5i=tersoffs%param(13,iter)
                  hi =tersoffs%param(14,iter)
                  ak =tersoffs%param(15,iter)
                  bk =tersoffs%param(16,iter)
                End If

                ! bond-angle detection

                Do jj=1,limit

                  ! index of the first secondary atom

                  jatm=listin(jj)
                  jter=tersoffs%list(config%ltype(jatm))

                  ! bypass if secondary atom type is not involved in interaction

                  If (jter > 0 .and. jatm /= iatm) Then

                    ! interaction index

                    ijter=(Max(iter,jter)*(Max(iter,jter)-1))/2+Min(iter,jter)

                    ! if pair is within cutoff

                    If (rtf(jj) <= tersoffs%cut(ijter)) Then

                      ! triplet occurrence flag

                      flag3=.false.

                      ! potential energy and virial terms

                      vterm=0.0_wp
                      eterm=0.0_wp

                      ! initialise work arrays

                      Do kk=1,limit
                        cst(kk)=0.0_wp
                        gam(kk)=0.0_wp
                        gvr(kk)=0.0_wp
                      End Do
                      If (tersoffs%key_pot == 2) Then ! KIHS
                        Do kk=1,limit
                          rkj(kk)=0.0_wp
                          wkj(kk)=0.0_wp
                        End Do
                      End If

                      ! (SHIFT TO LEFT)
                      ! calculate bond factor

                      Do kk=1,limit

                        ! index of the second secondary atom

                        katm=listin(kk)
                        kter=tersoffs%list(config%ltype(katm))

                        ! bypass if secondary atom type is not involved in interaction

                        If (kter > 0 .and. katm /= iatm .and. katm /= jatm) Then

                          ! interaction index

                          ikter=(Max(iter,kter)*(Max(iter,kter)-1))/2+Min(iter,kter)

                          ! if pair is within cutoff

                          If (rtf(kk) <= tersoffs%cut(ikter)) Then

                            ! only for not fully frozen triplets

                            If (config%lfrzn(iatm)*config%lfrzn(jatm)*config%lfrzn(katm) == 0) Then

                              flag3 = .true.

                              ! bond-angle is detected, force and virial parameters calculation

                              cost = xtf(jj)*xtf(kk)+ytf(jj)*ytf(kk)+ztf(jj)*ztf(kk)
                              If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)
                              cst(kk) = cost

                              If      (tersoffs%key_pot == 1) Then ! TERS
                                gtheta= 1.0_wp + (ci/di)**2 - ci**2 / (di**2 + (hi-cost)**2)
                                eterm = eterm + gtheta*tersoffs%param2(ikter,2)*scr(kk) ! L_{ij}
                                vterm = vterm + gtheta*tersoffs%param2(ikter,2)*gcr(kk)*rtf(kk)
                                ! d/dr_k of L_{ij} - angular part as it is used in the virial

                                gam(kk) = gtheta
                                gvr(kk) = 2.0_wp * ci**2 * (hi-cost) / (di**2 + (hi-cost)**2)**2 ! d(gtheta)/sint*d(theta)
                              Else If (tersoffs%key_pot == 2) Then ! KIHS
                                rkj(kk)=rtf(jj)-rtf(kk)
                                wkj(kk)=Exp(ak * rkj(kk)**bk)

                                hmct2=(hi-cost)**2
                                c4exp=c4i*Exp(-c5i*hmct2)

                                gtheta= c1i + c2i*hmct2*(1.0_wp+c4exp)/(c3i+hmct2)
                                eterm = eterm + gtheta*wkj(kk)*scr(kk) ! L_{ij}
                                vterm = vterm + gtheta*wkj(kk)*(gcr(kk)*rtf(kk) + scr(kk)*ak*bk*rkj(kk)**bk)
                                ! d/dr_k of L_{ij} - angular part as it is used in the virial

                                gam(kk) = gtheta
                                gvr(kk) = 2.0_wp * (c2i/(c3i+hmct2)) * (hi-cost) * & ! d(gtheta)/sint*d(theta)
                                  ((c3i/(c3i+hmct2))*(1.0_wp+c4exp) - hmct2*c5i*c4exp)
                              End If

                            End If

                          End If

                        End If

                      End Do

                      ! calculate contribution to energy, virial, two-body stress and forces
                      ! (all associated with the head atom)

                      If      (tersoffs%key_pot == 1) Then ! TERS
                        gam_ij=tersoffs%param2(iter,1)
                        gamma=0.0_wp
                        If (flag3) Then
                          gam_ij = tersoffs%param2(iter,1)*(1.0_wp+(bi*eterm)**ei)**(-0.5_wp/ei) ! gamma_{ij}
                          gamma  = eat(jj) * tersoffs%param2(iter,1) * bi*(bi*eterm)**(ei-1.0_wp) * &
                            0.5_wp*(1.0_wp+(bi*eterm)**ei)**(-0.5_wp/ei - 1.0_wp) ! -FcFa[d/dr gamma_{ij}]/[d/dr Lij]
                        End If
                      Else If (tersoffs%key_pot == 2) Then ! KIHS
                        gam_ij=1.0_wp
                        gamma=0.0_wp
                        If (flag3) Then
                          gam_ij = (1.0_wp+eterm**ei)**(-di) ! gamma_{ij}
                          gamma  = eat(jj) * ei * eterm**(ei-1.0_wp) * &
                            di * (1.0_wp+eterm**ei)**(-di-1.0_wp) ! -FcFa[d/dr gamma_{ij}]/[d/dr Lij]
                        End If
                      End If

                      gterm=(grt(jj)-gam_ij*gat(jj)) ! force_ij (no three body part)

                      ! Counteract the double counting with the 0.5_wp factor

                      If (iatm <= config%natms) Then

                        stats%engter = stats%engter + 0.5_wp*(ert(jj) - gam_ij*eat(jj))    ! energy_ij
                        stats%virter = stats%virter + 0.5_wp*(gamma*vterm + gterm*rtf(jj)) ! virial_ij

                        config%parts(iatm)%fxx=config%parts(iatm)%fxx+0.5_wp*gterm*xtf(jj)
                        config%parts(iatm)%fyy=config%parts(iatm)%fyy+0.5_wp*gterm*ytf(jj)
                        config%parts(iatm)%fzz=config%parts(iatm)%fzz+0.5_wp*gterm*ztf(jj)

                        strs1 = strs1 - 0.5_wp*gterm*rtf(jj)*xtf(jj)*xtf(jj)
                        strs2 = strs2 - 0.5_wp*gterm*rtf(jj)*xtf(jj)*ytf(jj)
                        strs3 = strs3 - 0.5_wp*gterm*rtf(jj)*xtf(jj)*ztf(jj)
                        strs5 = strs5 - 0.5_wp*gterm*rtf(jj)*ytf(jj)*ytf(jj)
                        strs6 = strs6 - 0.5_wp*gterm*rtf(jj)*ytf(jj)*ztf(jj)
                        strs9 = strs9 - 0.5_wp*gterm*rtf(jj)*ztf(jj)*ztf(jj)

                      End If

                      If (jatm <= config%natms) Then

                        config%parts(jatm)%fxx=config%parts(jatm)%fxx-0.5_wp*gterm*xtf(jj)
                        config%parts(jatm)%fyy=config%parts(jatm)%fyy-0.5_wp*gterm*ytf(jj)
                        config%parts(jatm)%fzz=config%parts(jatm)%fzz-0.5_wp*gterm*ztf(jj)

                      End If

                      ! calculate three-body forces contributions

                      Do kk=1,limit

                        ! index of the second secondary atom

                        katm=listin(kk)
                        kter=tersoffs%list(config%ltype(katm))

                        ! bypass if secondary atom type is not involved in interaction

                        If (kter > 0 .and. katm /= iatm .and. katm /= jatm) Then

                          ! interaction index

                          ikter=(Max(iter,kter)*(Max(iter,kter)-1))/2+Min(iter,kter)

                          ! if pair is within cutoff

                          If (rtf(kk) <= tersoffs%cut(ikter)) Then

                            ! only for not fully frozen triplets

                            If (config%lfrzn(iatm)*config%lfrzn(jatm)*config%lfrzn(katm) == 0) Then

                              ! term in d/dr_k L_{ij} no derivative of gamma

                              cost=cst(kk)

                              ! Counteract the double counting with the 0.5_wp factor

                              If      (tersoffs%key_pot == 1) Then ! TERS
                                gam_dg = 0.5_wp*gamma*tersoffs%param2(ikter,2)*scr(kk)*gvr(kk)
                                gam_df = 0.5_wp*gamma*tersoffs%param2(ikter,2)*gcr(kk)*gam(kk)

                                ! calculate contribution to atomic forces

                                fxj = gam_dg*(xtf(kk)-xtf(jj)*cost)/rtf(jj) ! contributions to j
                                fyj = gam_dg*(ytf(kk)-ytf(jj)*cost)/rtf(jj)
                                fzj = gam_dg*(ztf(kk)-ztf(jj)*cost)/rtf(jj)

                                fxk = gam_dg*(xtf(jj)-xtf(kk)*cost)/rtf(kk) - gam_df*xtf(kk) ! contributions to k
                                fyk = gam_dg*(ytf(jj)-ytf(kk)*cost)/rtf(kk) - gam_df*ytf(kk)
                                fzk = gam_dg*(ztf(jj)-ztf(kk)*cost)/rtf(kk) - gam_df*ztf(kk)
                              Else If (tersoffs%key_pot == 2) Then ! KIHS
                                gam_dg = 0.5_wp*gamma*scr(kk)*wkj(kk)*gvr(kk)
                                gam_df = 0.5_wp*gamma*gam(kk)*wkj(kk)*gcr(kk)
                                gam_dw = 0.5_wp*gamma*scr(kk)*gam(kk)*wkj(kk)*(ak*bk*rkj(kk)**(bk-1.0_wp))

                                ! calculate contribution to atomic forces

                                fxj = gam_dg*(xtf(kk)-xtf(jj)*cost)/rtf(jj) - gam_dw*xtf(jj) ! contributions to j
                                fyj = gam_dg*(ytf(kk)-ytf(jj)*cost)/rtf(jj) - gam_dw*ytf(jj)
                                fzj = gam_dg*(ztf(kk)-ztf(jj)*cost)/rtf(jj) - gam_dw*ztf(jj)

                                fxk = gam_dg*(xtf(jj)-xtf(kk)*cost)/rtf(kk) - gam_df*xtf(kk) + gam_dw*xtf(kk)! contributions to k
                                fyk = gam_dg*(ytf(jj)-ytf(kk)*cost)/rtf(kk) - gam_df*ytf(kk) + gam_dw*ytf(kk)
                                fzk = gam_dg*(ztf(jj)-ztf(kk)*cost)/rtf(kk) - gam_df*ztf(kk) + gam_dw*ztf(kk)
                              End If

                              If (iatm <= config%natms) Then

                                config%parts(iatm)%fxx=config%parts(iatm)%fxx-(fxj+fxk)
                                config%parts(iatm)%fyy=config%parts(iatm)%fyy-(fyj+fyk)
                                config%parts(iatm)%fzz=config%parts(iatm)%fzz-(fzj+fzk)

                                ! calculate contribution to stress tensor (associated to the head atom)

                                strs1 = strs1 + (fxj*xtf(jj)*rtf(jj) + fxk*xtf(kk)*rtf(kk))
                                strs2 = strs2 + (fxj*ytf(jj)*rtf(jj) + fxk*ytf(kk)*rtf(kk))
                                strs3 = strs3 + (fxj*ztf(jj)*rtf(jj) + fxk*ztf(kk)*rtf(kk))
                                strs5 = strs5 + (fyj*ytf(jj)*rtf(jj) + fyk*ytf(kk)*rtf(kk))
                                strs6 = strs6 + (fyj*ztf(jj)*rtf(jj) + fyk*ztf(kk)*rtf(kk))
                                strs9 = strs9 + (fzj*ztf(jj)*rtf(jj) + fzk*ztf(kk)*rtf(kk))

                              End If

                              If (jatm <= config%natms) Then

                                config%parts(jatm)%fxx=config%parts(jatm)%fxx+fxj
                                config%parts(jatm)%fyy=config%parts(jatm)%fyy+fyj
                                config%parts(jatm)%fzz=config%parts(jatm)%fzz+fzj

                              End If

                              If (katm <= config%natms) Then

                                config%parts(katm)%fxx=config%parts(katm)%fxx+fxk
                                config%parts(katm)%fyy=config%parts(katm)%fyy+fyk
                                config%parts(katm)%fzz=config%parts(katm)%fzz+fzk

                              End If

                            End If

                          End If

                        End If

                      End Do

                      ! (BACK - SHIFT TO LEFT)

                    End If

                  End If

                End Do

              End If

            End Do

          End If

        End Do
      End Do
    End Do

    ! global sum of tersoff potential and virial


    buffer(1)=stats%engter
    buffer(2)=stats%virter
    Call gsum(comm,buffer(1:2))
    stats%engter=buffer(1)
    stats%virter=buffer(2)


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

    Deallocate (link,listin,lct,lst,      Stat=fail(1))
    Deallocate (xxt,yyt,zzt,              Stat=fail(2))
    Deallocate (xtf,ytf,ztf,rtf,          Stat=fail(3))
    Deallocate (ert,eat,grt,gat,          Stat=fail(4))
    Deallocate (scr,gcr,                  Stat=fail(5))
    Deallocate (cst,gam,gvr,              Stat=fail(6))
    If (tersoffs%key_pot == 2) Deallocate (rkj,wkj, Stat=fail(7))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'tersoff_forces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine tersoff_forces

  Subroutine tersoff_generate(tersoffs)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for generating potential energy and force arrays
    ! for tersoff forces only, based on potential forms as defined by
    ! 1) J. Tersoff, Phys. Rev. B 39 (1989) 5566 and
    ! 2) T. Kumagai, S. Izumi, S. Hara, and S. Sakai, Computational
    !    Materials Science 39, 457 (2007), ISSN 09270256
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith  october 2004
    ! amended   - i.t.todorov september 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( tersoff_type ), Intent( InOut )  :: tersoffs

    Integer           :: i,katom1,katom2,ipt,jpt,kpt
    Real( Kind = wp ) :: dlrpot,baij,saij,bbij,sbij,rij,sij,rrr,arg,rep,att

    ! define grid resolution for potential arrays

    dlrpot=tersoffs%cutoff/Real(tersoffs%max_grid-4,wp)

    ! construct arrays for all types of tersoff potential

    Do katom1=1,tersoffs%n_potential
      Do katom2=1,katom1

        If (tersoffs%ltp(katom1)*tersoffs%ltp(katom2) /= 0) Then

          ipt=tersoffs%list(katom1)
          jpt=tersoffs%list(katom2)
          kpt=(Max(ipt,jpt)*(Max(ipt,jpt)-1))/2+Min(ipt,jpt)

          ! define tersoff parameters

          baij =    Sqrt(tersoffs%param(1,ipt)*tersoffs%param(1,jpt))
          saij = 0.5_wp*(tersoffs%param(2,ipt)+tersoffs%param(2,jpt))
          bbij =    Sqrt(tersoffs%param(3,ipt)*tersoffs%param(3,jpt))
          sbij = 0.5_wp*(tersoffs%param(4,ipt)+tersoffs%param(4,jpt))
          rij  =    Sqrt(tersoffs%param(5,ipt)*tersoffs%param(5,jpt))
          sij  =    Sqrt(tersoffs%param(6,ipt)*tersoffs%param(6,jpt))

          ! store potential cutoff

          tersoffs%cut(kpt) = sij

          ! calculate screening generic function

          Do i=0,tersoffs%max_grid
            rrr=Real(i,wp)*dlrpot

            If      (rrr <= rij) Then
              tersoffs%vmbp(i,kpt,1)=1.0_wp
              tersoffs%gmbp(i,kpt,1)=0.0_wp
            Else
              If (rrr <= sij) Then
                arg=pi*(rrr-rij)/(sij-rij)

                If      (tersoffs%ltp(katom1)*tersoffs%ltp(katom2) == 1) Then

                  tersoffs%vmbp(i,kpt,1)=0.5_wp*(1.0_wp+Cos(arg))
                  tersoffs%gmbp(i,kpt,1)=0.5_wp*pi*rrr*Sin(arg)/(sij-rij)

                Else If (tersoffs%ltp(katom1)*tersoffs%ltp(katom2) == 4) Then

                  ! Murty's correction to screening function (tersoffs%vmbp)
                  ! M.V.R. Murty, H.A. Atwater, Phys. Rev. B 51 (1995) 4889-4993

                  tersoffs%vmbp(i,kpt,1)=0.5_wp+9.0_wp/16.0_wp*Cos(arg)-1.0_wp/16.0_wp*Cos(3.0_wp*arg)
                  tersoffs%gmbp(i,kpt,1)=0.75_wp*pi*rrr*(Sin(arg))**3/(sij-rij)

                End If
              End If
            End If
          End Do

          ! calculate screening repulsion & attraction functions

          Do i=0,tersoffs%max_grid
            rrr=Real(i,wp)*dlrpot

            ! repulsion

            rep=baij*Exp(-saij*rrr)

            tersoffs%vmbp(i,kpt,2)=rep*  tersoffs%vmbp(i,kpt,1)
            tersoffs%gmbp(i,kpt,2)=rep*( tersoffs%gmbp(i,kpt,1) + saij*rrr*tersoffs%vmbp(i,kpt,1) )

            ! attraction

            att=bbij*Exp(-sbij*rrr)

            tersoffs%vmbp(i,kpt,3)=att*  tersoffs%vmbp(i,kpt,1)
            tersoffs%gmbp(i,kpt,3)=att*( tersoffs%gmbp(i,kpt,1) + sbij*rrr*tersoffs%vmbp(i,kpt,1) )
          End Do

        End If

      End Do
    End Do

  End Subroutine tersoff_generate

  Subroutine cleanup(T)
    Type( tersoff_type ) :: T

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
    If (Allocated(T%param2)) Then
      Deallocate(T%param2)
    End If

    If (Allocated(T%vmbp)) Then
      Deallocate(T%vmbp)
    End If
    If (Allocated(T%gmbp)) Then
      Deallocate(T%gmbp)
    End If

    If (Allocated(T%cut)) Then
      Deallocate(T%cut)
    End If
  End Subroutine cleanup
End Module tersoff
