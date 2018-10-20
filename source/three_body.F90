Module three_body

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global three-body potential variables and
  ! arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov june 2008
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, wi

  Use comms,     Only : comms_type,gsum,gcheck
  Use constants, Only : zero_plus, twopi, pi
  Use domains, Only : domains_type
  Use configuration,  Only : configuration_type
  Use particle,        Only : corePart
  Use errors_warnings, Only : error,warning
  Use numerics, Only : dcell, invert
  Use statistics, Only : stats_type
  Use neighbours, Only : neighbours_type

  Implicit None

  Private

  Type, Public ::  threebody_type
    Integer( Kind = wi )               :: ntptbp = 0
    Integer( Kind = wi )               :: mxptbp, mxtbp, mx2tbp
    Real( Kind = wp )                  :: cutoff
    Logical,              Allocatable  :: lfrtbp(:)
    Integer( Kind = wi ), Allocatable  :: lsttbp(:),ltptbp(:)
    Real( Kind = wp ),    Allocatable  :: prmtbp(:,:),rcttbp(:)
  Contains
    Private
    Procedure, Public :: init => allocate_three_body_arrays
    Final ::  deallocate_three_body_arrays
  End Type

  Public ::  three_body_forces

Contains

  Subroutine allocate_three_body_arrays(T,max_site)
    Class(threebody_type) , Intent( InOut )   :: T
    Integer( Kind = wi ), Intent( In    ) :: max_site

    Integer :: fail(5)

    fail=0

    Allocate (T%lfrtbp(1:Merge(max_site,0,T%mxtbp > 0)), Stat = fail(1))
    Allocate (T%lsttbp(1:T%mxtbp),                     Stat = fail(2))
    Allocate (T%ltptbp(1:T%mxtbp),                     Stat = fail(3))
    Allocate (T%prmtbp(1:T%mxptbp,1:T%mxtbp),  Stat = fail(4))
    Allocate (T%rcttbp(1:T%mxtbp),                     Stat = fail(5))

    If (Any(fail > 0)) Call error(1024)

    T%lfrtbp = .false.
    T%lsttbp = 0
    T%ltptbp = 0
    T%prmtbp = 0.0_wp
    T%rcttbp = 0.0_wp
  End Subroutine allocate_three_body_arrays

  Subroutine deallocate_three_body_arrays(T)
    Type( threebody_type ), Intent( InOut ) :: T

    Integer :: fail(5)

    fail = 0 

    Deallocate (T%lfrtbp, Stat = fail(1))
    Deallocate (T%lsttbp,                     Stat = fail(2))
    Deallocate (T%ltptbp,                     Stat = fail(3))
    Deallocate (T%prmtbp,  Stat = fail(4))
    Deallocate (T%rcttbp,                     Stat = fail(5))

    If (Any(fail > 0)) Call error(1024)

  end subroutine deallocate_three_body_arrays

  Subroutine three_body_forces(stats,threebody,neigh,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating three-body forces (energy and
    ! virial) arising from the included angle between three atoms at a
    ! nominated central atom
    !
    ! Note: coordinates are converted to reduced units to avoid a call to
    ! images. The link cell algorithm used here necessitates a
    ! parallelepiped cell geometry.
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith march 1994
    ! amended   - i.t.todorov february 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type(stats_type), Intent( InOut ) :: stats
    Type(threebody_type), Intent( InOut ) :: threebody
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ),     Intent( InOut ) :: config
    Type(comms_type), Intent( InOut ) :: comm

    Logical           :: safe

    Integer           :: fail(1:2),                      &
      i,j,k, ii,jj,kk,jk, ia,ib,ic,   &
      nbx,nby,nbz, ncells,            &
      ix,iy,iz,icell, jx,jy,jz,jcell, &
      itbp,jtbp,ktbp,                 &
      jktbp,kktbp,limit,last, ktyp

    Real( Kind = wp ) :: dispx,dispy,dispz, xdc,ydc,zdc,             &
      rcell(1:9),celprp(1:10),det,                &
      xab,yab,zab,rab,rrab, sxab,syab,szab,       &
      xac,yac,zac,rac,rrac,                       &
      xbc,ybc,zbc,rbc,rrbc, sxbc,sybc,szbc,       &
      theta,cost,sint,rsint,                      &
      fxa,fxc,fya, fyc,fza,fzc,                   &
      k0,theta0,dtheta,dthpi,dth0pi,dth,          &
      a,rho,rho1,rho2,switch,dhb,rhb,term,        &
      pterm,vterm,gamma,gamsa,gamsb,gamsc,        &
      strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

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

    nbx=Int(domain%nx_recip*celprp(7)/(threebody%cutoff+1.0e-6_wp))
    nby=Int(domain%ny_recip*celprp(8)/(threebody%cutoff+1.0e-6_wp))
    nbz=Int(domain%nz_recip*celprp(9)/(threebody%cutoff+1.0e-6_wp))

    ! check for link cell algorithm violations

    If (nbx < 3 .or. nby < 3 .or. nbz < 3) Call error(305)

    ncells=(nbx+4)*(nby+4)*(nbz+4)
    If (ncells > neigh%max_cell) Then
      Call warning(90,Real(ncells,wp),Real(neigh%max_cell,wp),1.0_wp)
      neigh%max_cell = Nint(1.25_wp*Real(ncells,wp))
      If (ncells > config%mxatms) Call error(69)
    End If

    fail=0
    Allocate (link(1:config%mxatms),listin(1:config%mxatms),lct(1:ncells),lst(1:ncells), Stat=fail(1))
    Allocate (xxt(1:config%mxatms),yyt(1:config%mxatms),zzt(1:config%mxatms),                   Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'three_body_forces allocation failure'
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
      If (threebody%lfrtbp(config%ltype(i))) Then
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
      If (threebody%lfrtbp(config%ltype(i))) Then

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
        ! (comm%idnode) but due to some tiny numerical inaccuracy kicked into
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

    stats%engtbp=0.0_wp
    stats%virtbp=0.0_wp

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

              ib=listin(ii)
              itbp=threebody%mx2tbp*(config%ltype(ib)-1)

              ! bypass if primary atom type is not involved in interaction

              If (threebody%lsttbp(itbp+1) >= 0) Then

                last=limit

                ! loop over half the length of listin
                ! (to get the second secondary atom)

                Do kk=1,limit/2

                  ! correction to the loop over the full length of listin
                  ! when the length is an even number

                  If (kk > (limit-1)/2) last=limit/2

                  ! loop over the full length of listin
                  ! (to get the first secondary atom)

                  Do jj=1,last

                    ! index of the first secondary atom

                    ia=listin(jj)

                    jk=jj+kk
                    If (jk > limit) jk=jk-limit

                    ! index of the second secondary atom
                    ! (definitely different form the first one)

                    ic=listin(jk)

                    ! if neither of the secondary atoms coincides with the head one

                    If (ia /= ib .and. ic /= ib) Then

                      ! get types of atoms

                      jtbp=Max(config%ltype(ia),config%ltype(ic))
                      ktbp=Min(config%ltype(ia),config%ltype(ic))

                      ! get index of interaction

                      jktbp=itbp+(jtbp*(jtbp-1))/2+ktbp
                      kktbp=threebody%lsttbp(jktbp)

                      ! check existence of interaction in system

                      If (kktbp > 0) Then

                        ! only for non-frozen triplets

                        If (config%lfrzn(ia)*config%lfrzn(ib)*config%lfrzn(ic) == 0) Then

                          sxab = xxt(ia)-xxt(ib)
                          syab = yyt(ia)-yyt(ib)
                          szab = zzt(ia)-zzt(ib)

                          xab=config%cell(1)*sxab+config%cell(4)*syab+config%cell(7)*szab
                          If (Abs(xab) < threebody%cutoff) Then

                            yab=config%cell(2)*sxab+config%cell(5)*syab+config%cell(8)*szab
                            If (Abs(yab) < threebody%cutoff) Then

                              zab=config%cell(3)*sxab+config%cell(6)*syab+config%cell(9)*szab
                              If (Abs(zab) < threebody%cutoff) Then

                                rab=Sqrt(xab*xab+yab*yab+zab*zab)

                                sxbc = xxt(ic)-xxt(ib)
                                sybc = yyt(ic)-yyt(ib)
                                szbc = zzt(ic)-zzt(ib)

                                xbc=config%cell(1)*sxbc+config%cell(4)*sybc+config%cell(7)*szbc
                                If (Abs(xbc) < threebody%cutoff) Then

                                  ybc=config%cell(2)*sxbc+config%cell(5)*sybc+config%cell(8)*szbc
                                  If (Abs(ybc) < threebody%cutoff) Then

                                    zbc=config%cell(3)*sxbc+config%cell(6)*sybc+config%cell(9)*szbc
                                    If (Abs(zbc) < threebody%cutoff) Then

                                      rbc=Sqrt(xbc*xbc+ybc*ybc+zbc*zbc)

                                      If (Max(rab,rbc) <= threebody%rcttbp(kktbp)) Then

                                        xac = xab - xbc
                                        yac = yab - ybc
                                        zac = zab - zbc
                                        rac=Sqrt(xac*xac+yac*yac+zac*zac)

                                        rrab = 1.0_wp/rab
                                        rrbc = 1.0_wp/rbc
                                        rrac = 1.0_wp/rac

                                        ! normalise direction vectors

                                        xab = xab*rrab
                                        yab = yab*rrab
                                        zab = zab*rrab

                                        xbc = xbc*rrbc
                                        ybc = ybc*rrbc
                                        zbc = zbc*rrbc

                                        xac = xac*rrac
                                        yac = yac*rrac
                                        zac = zac*rrac

                                        cost=(xab*xbc+yab*ybc+zab*zbc)
                                        If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)

                                        ! select potential energy function type

                                        ktyp=threebody%ltptbp(kktbp)

                                        If (ktyp /= 6) Then

                                          ! for potentials different from dreiding/CHARMM hydrogen bond

                                          sint=Max(1.0e-10_wp,Sqrt(1.0_wp-cost**2))
                                          rsint=1.0_wp/sint
                                          theta=Acos(cost)

                                        End If

                                        If      (ktyp == 1) Then

                                          ! harmonic valence angle potential

                                          k0    =threebody%prmtbp(1,kktbp)
                                          theta0=threebody%prmtbp(2,kktbp)
                                          dtheta=theta-theta0

                                          pterm=k0*dtheta

                                          gamma=pterm*rsint
                                          pterm=pterm*0.5_wp*dtheta
                                          vterm=0.0_wp
                                          gamsa=0.0_wp
                                          gamsc=0.0_wp
                                          gamsb=0.0_wp

                                        Else If (ktyp == 2) Then

                                          ! truncated harmonic valence angle potential

                                          k0    =threebody%prmtbp(1,kktbp)
                                          theta0=threebody%prmtbp(2,kktbp)
                                          dtheta=theta-theta0
                                          rho   =threebody%prmtbp(3,kktbp)
                                          switch=-(rab**8+rbc**8)/rho**8

                                          pterm=k0*dtheta*Exp(switch)

                                          gamma=pterm*rsint
                                          pterm=pterm*0.5_wp*dtheta
                                          vterm=pterm*8.0_wp*switch
                                          gamsa=pterm*8.0_wp*rab**7/rho**8
                                          gamsc=pterm*8.0_wp*rbc**7/rho**8
                                          gamsb=0.0_wp

                                        Else If (ktyp == 3) Then

                                          ! screened harmonic valence angle potential

                                          k0    =threebody%prmtbp(1,kktbp)
                                          theta0=threebody%prmtbp(2,kktbp)
                                          dtheta=theta-theta0
                                          rho1  =threebody%prmtbp(3,kktbp)
                                          rho2  =threebody%prmtbp(4,kktbp)
                                          switch=-(rab/rho1+rbc/rho2)

                                          pterm=k0*dtheta*Exp(switch)

                                          gamma=pterm*rsint
                                          pterm=pterm*0.5_wp*dtheta
                                          vterm=pterm*switch
                                          gamsa=pterm/rho1
                                          gamsc=pterm/rho2

                                          gamsa=(pterm/threebody%prmtbp(3,kktbp))
                                          gamsc=(pterm/threebody%prmtbp(4,kktbp))
                                          gamsb=0.0_wp

                                        Else If (ktyp == 4) Then

                                          ! screened Vessal potential type 1

                                          k0    =threebody%prmtbp(1,kktbp)
                                          theta0=threebody%prmtbp(2,kktbp)
                                          dth0pi=theta0-pi
                                          dthpi =theta -pi
                                          dth   =dth0pi**2-dthpi**2
                                          rho1  =threebody%prmtbp(3,kktbp)
                                          rho2  =threebody%prmtbp(4,kktbp)
                                          switch=-(rab/rho1+rbc/rho2)

                                          pterm=(k0*dth/(2.0_wp*dth0pi**2)) * Exp(switch)

                                          gamma=pterm*dthpi*rsint
                                          pterm=pterm*0.25_wp*dth
                                          vterm=pterm*switch
                                          gamsa=pterm/rho1
                                          gamsc=pterm/rho2
                                          gamsb=0.0_wp

                                        Else If (ktyp == 5) Then

                                          ! truncated Vessal potential type 2

                                          k0    =threebody%prmtbp(1,kktbp)
                                          theta0=threebody%prmtbp(2,kktbp)
                                          dtheta=theta-theta0
                                          dth0pi=theta0-pi
                                          dthpi =theta -pi
                                          a     =threebody%prmtbp(3,kktbp)
                                          rho   =threebody%prmtbp(4,kktbp)
                                          switch=-(rab**8+rbc**8)/rho**8

                                          pterm=k0*dtheta*Exp(switch)

                                          gamma=pterm * (theta**(a-1.0_wp) * (dthpi+dth0pi) * &
                                            ((a+4.0_wp)*theta**2 - twopi*(a+2.0_wp)*theta - a*theta0*(dth0pi-pi)) + &
                                            a*pi**(a-1.0_wp) * dth0pi**3) * rsint

                                          pterm=pterm * dtheta * &
                                            (theta**a * (dthpi+dth0pi)**2 + 0.5_wp*a * pi**(a-1.0_wp) * dth0pi**3)
                                          vterm=pterm*8.0_wp*switch
                                          gamsa=pterm*8.0_wp*rab**7/rho**8
                                          gamsc=pterm*8.0_wp*rbc**7/rho**8
                                          gamsb=0.0_wp

                                        Else If (ktyp == 6) Then

                                          ! dreiding/CHARMM hydrogen bond

                                          If (Min(rab,rbc) < 1.5_wp .and. rac < threebody%rcttbp(kktbp)) Then
                                            dhb   =threebody%prmtbp(1,kktbp)
                                            rhb   =threebody%prmtbp(2,kktbp)
                                            switch=(5.0_wp*(rhb/rac)**2-6.0_wp)*(rhb/rac)**10
                                            term  =60.0_wp*(1.0_wp-(rhb/rac)**2)*(rhb/rac)**10

                                            pterm=dhb*cost**3

                                            gamma=pterm*(-4.0_wp)*switch
                                            pterm=pterm*cost
                                            vterm=pterm*term
                                            pterm=pterm*switch
                                            gamsa=0.0_wp
                                            gamsc=0.0_wp
                                            gamsb=-vterm/rac
                                          Else
                                            gamma=0.0_wp
                                            pterm=0.0_wp
                                            vterm=0.0_wp
                                            gamsa=0.0_wp
                                            gamsc=0.0_wp
                                            gamsb=0.0_wp
                                          End If

                                        Else

                                          safe=.false.
                                          gamma=0.0_wp
                                          pterm=0.0_wp
                                          vterm=0.0_wp
                                          gamsa=0.0_wp
                                          gamsc=0.0_wp
                                          gamsb=0.0_wp

                                        End If

                                        ! calculate atomic forces

                                        fxa = gamma*(xbc-xab*cost)*rrab+gamsa*xab+gamsb*xac
                                        fya = gamma*(ybc-yab*cost)*rrab+gamsa*yab+gamsb*yac
                                        fza = gamma*(zbc-zab*cost)*rrab+gamsa*zab+gamsb*zac

                                        fxc = gamma*(xab-xbc*cost)*rrbc+gamsc*xbc-gamsb*xac
                                        fyc = gamma*(yab-ybc*cost)*rrbc+gamsc*ybc-gamsb*yac
                                        fzc = gamma*(zab-zbc*cost)*rrbc+gamsc*zbc-gamsb*zac

                                        If (ia <= config%natms) Then

                                          config%parts(ia)%fxx=config%parts(ia)%fxx+fxa
                                          config%parts(ia)%fyy=config%parts(ia)%fyy+fya
                                          config%parts(ia)%fzz=config%parts(ia)%fzz+fza

                                        End If

                                        If (ib <= config%natms) Then

                                          ! energy and virial (associated to the head atom)

                                          stats%engtbp=stats%engtbp+pterm
                                          stats%virtbp=stats%virtbp+vterm

                                          ! calculate stress tensor (associated to the head atom)

                                          strs1 = strs1 + rab*xab*fxa + rbc*xbc*fxc
                                          strs2 = strs2 + rab*xab*fya + rbc*xbc*fyc
                                          strs3 = strs3 + rab*xab*fza + rbc*xbc*fzc
                                          strs5 = strs5 + rab*yab*fya + rbc*ybc*fyc
                                          strs6 = strs6 + rab*yab*fza + rbc*ybc*fzc
                                          strs9 = strs9 + rab*zab*fza + rbc*zbc*fzc

                                          config%parts(ib)%fxx=config%parts(ib)%fxx-(fxa+fxc)
                                          config%parts(ib)%fyy=config%parts(ib)%fyy-(fya+fyc)
                                          config%parts(ib)%fzz=config%parts(ib)%fzz-(fza+fzc)

                                        End If

                                        If (ic <= config%natms) Then

                                          config%parts(ic)%fxx=config%parts(ic)%fxx+fxc
                                          config%parts(ic)%fyy=config%parts(ic)%fyy+fyc
                                          config%parts(ic)%fzz=config%parts(ic)%fzz+fzc

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
    If (.not.safe) Call error(442)

    ! global sum of three-body potential and virial


    buffer(1)=stats%engtbp
    buffer(2)=stats%virtbp
    Call gsum(comm,buffer(1:2))
    stats%engtbp=buffer(1)
    stats%virtbp=buffer(2)

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
      Write(message,'(a)') 'three_body_forces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine three_body_forces


End Module three_body
