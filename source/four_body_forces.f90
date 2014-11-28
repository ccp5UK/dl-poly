Subroutine four_body_forces(imcon,rcfbp,engfbp,virfbp,stress)

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
! amended   - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum,gcheck
  Use setup_module
  Use domains_module, Only : idx,idy,idz, nprx,npry,nprz, &
                             r_nprx,r_npry,r_nprz
  Use config_module,  Only : cell,natms,nlast,lfrzn,ltype, &
                             xxx,yyy,zzz,fxx,fyy,fzz
  Use four_body_module

  Implicit None

  Integer,                             Intent( In    ) :: imcon
  Real( Kind = wp ),                   Intent( In    ) :: rcfbp
  Real( Kind = wp ),                   Intent(   Out ) :: engfbp,virfbp
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

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


! image conditions not compliant with DD and link-cell

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Calculate the number of link-cells per domain in every direction

  nbx=Int(r_nprx*celprp(7)/(rcfbp+1.0e-6_wp))
  nby=Int(r_npry*celprp(8)/(rcfbp+1.0e-6_wp))
  nbz=Int(r_nprz*celprp(9)/(rcfbp+1.0e-6_wp))

! check for link cell algorithm violations

  If (nbx < 3 .or. nby < 3 .or. nbz < 3) Call error(305)

  ncells=(nbx+4)*(nby+4)*(nbz+4)
  If (ncells > mxcell) Then
     Call warning(90,Real(ncells,wp),Real(mxcell,wp),2.0_wp)
     mxcell = Nint(1.25_wp*Real(ncells,wp))
  End If

  fail=0
  Allocate (link(1:mxatms),listin(1:mxatms),lct(1:ncells),lst(1:ncells), Stat=fail(1))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),                   Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'four_body_forces allocation failure, node: ', idnode
     Call error(0)
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

  dispx=0.5_wp-r_nprx*(Real(idx,wp)-2.0_wp/Real(nbx,wp))
  dispy=0.5_wp-r_npry*(Real(idy,wp)-2.0_wp/Real(nby,wp))
  dispz=0.5_wp-r_nprz*(Real(idz,wp)-2.0_wp/Real(nbz,wp))

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions (ALL - halo included) from centred
! Cartesian coordinates to reduced space coordinates of
! the left-most link-cell

  Do i=1,nlast
     If (lfrfbp(ltype(i))) Then
        xxt(i)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)+dispx
        yyt(i)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)+dispy
        zzt(i)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)+dispz
     End If
  End Do

! Form linked list
! Initialise link arrays

  link=0
  lct=0
  lst=0

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nbx*nprx,wp)
  ydc=Real(nby*npry,wp)
  zdc=Real(nbz*nprz,wp)

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

  Do i=1,nlast
     If (lfrfbp(ltype(i))) Then

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

  engfbp=0.0_wp
  virfbp=0.0_wp

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
! and its length (mini-list of neighbour cell contents)

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
              ifbp=mx3fbp*(ltype(ia)-1)

! bypass if primary atom type is not involved in interaction

              If (lstfbp(ifbp+1) >= 0) Then

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

     jfbp=Max(ltype(ib),ltype(ic),ltype(id))
     lfbp=Min(ltype(ib),ltype(ic),ltype(id))

! get index of interaction

     kfbp=ltype(ib)+ltype(ic)+ltype(id)-jfbp-lfbp
     jklbd=ifbp+lfbp+(kfbp*(kfbp-1))/2+(jfbp*(jfbp**2-1))/6
     kkfbp=lstfbp(jklbd)

! check existence of interaction in system

     If (kkfbp > 0) Then

! only for non-frozen triplets

        If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic)*lfrzn(id) == 0) Then

           sxab = xxt(ib)-xxt(ia)
           syab = yyt(ib)-yyt(ia)
           szab = zzt(ib)-zzt(ia)

           xab=cell(1)*sxab+cell(4)*syab+cell(7)*szab
           If (Abs(xab) < rcfbp) Then

              yab=cell(2)*sxab+cell(5)*syab+cell(8)*szab
              If (Abs(yab) < rcfbp) Then

                 zab=cell(3)*sxab+cell(6)*syab+cell(9)*szab
                 If (Abs(zab) < rcfbp) Then

                    rab2=xab*xab+yab*yab+zab*zab

                    sxac = xxt(ic)-xxt(ia)
                    syac = yyt(ic)-yyt(ia)
                    szac = zzt(ic)-zzt(ia)

                    xac=cell(1)*sxac+cell(4)*syac+cell(7)*szac
                    If (Abs(xac) < rcfbp) Then

                       yac=cell(2)*sxac+cell(5)*syac+cell(8)*szac
                       If (Abs(yac) < rcfbp) Then

                          zac=cell(3)*sxac+cell(6)*syac+cell(9)*szac
                          If (Abs(zac) < rcfbp) Then

                             rac2=xac*xac+yac*yac+zac*zac

                             sxad = xxt(id)-xxt(ia)
                             syad = yyt(id)-yyt(ia)
                             szad = zzt(id)-zzt(ia)

                             xad=cell(1)*sxad+cell(4)*syad+cell(7)*szad
                             If (Abs(xad) < rcfbp) Then

                                yad=cell(2)*sxad+cell(5)*syad+cell(8)*szad
                                If (Abs(yad) < rcfbp) Then

                                   zad=cell(3)*sxad+cell(6)*syad+cell(9)*szad
                                   If (Abs(zad) < rcfbp) Then

                                      rad2=xad*xad+yad*yad+zad*zad

! (SECOND SHIFT TO LEFT)

  If (Max(rab2,rac2,rad2) <= rctfbp(kkfbp)**2) Then

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

     ktyp=ltpfbp(kkfbp)

! calculate potential energy and scalar force term

     If      (ktyp == 1) Then

! harmonic inversion potential

        k0 =prmfbp(1,kkfbp)
        th0=prmfbp(2,kkfbp)

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

        k0  =prmfbp(1,kkfbp)
        cos0=prmfbp(2,kkfbp)

        potfbp=k0*((cosb-cos0)**2+(cosc-cos0)**2+(cosd-cos0)**2)/6.0_wp
        gamb=-k0*(cosb-cos0)/3.0_wp
        gamc=-k0*(cosc-cos0)/3.0_wp
        gamd=-k0*(cosd-cos0)/3.0_wp

     Else If (ktyp == 3) Then

! planar inversion potentials

        k0=prmfbp(1,kkfbp)

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

     If (ia <= natms) Then

! sum inversion energy

        engfbp=engfbp+potfbp

! stress tensor calculation for inversion terms

        strs1 = strs1 + xab*fbx + xac*fcx + xad*fdx
        strs2 = strs2 + yab*fbx + yac*fcx + yad*fdx
        strs3 = strs3 + zab*fbx + zac*fcx + zad*fdx
        strs5 = strs5 + yab*fby + yac*fcy + yad*fdy
        strs6 = strs6 + yab*fbz + yac*fcz + yad*fdz
        strs9 = strs9 + zab*fbz + zac*fcz + zad*fdz

        fxx(ia)=fxx(ia)+fax
        fyy(ia)=fyy(ia)+fay
        fzz(ia)=fzz(ia)+faz

     End If

     If (ib <= natms) Then

        fxx(ib)=fxx(ib)+fbx
        fyy(ib)=fyy(ib)+fby
        fzz(ib)=fzz(ib)+fbz

     End If

     If (ic <= natms) Then

        fxx(ic)=fxx(ic)+fcx
        fyy(ic)=fyy(ic)+fcy
        fzz(ic)=fzz(ic)+fcz

     End If

     If (id <= natms) Then

        fxx(id)=fxx(id)+fdx
        fyy(id)=fyy(id)+fdy
        fzz(id)=fzz(id)+fdz

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

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(443)

! global sum of four-body potential: virial is zero!!!

  If (mxnode > 1) Call gsum(engfbp)

! complete stress tensor

  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9

  Deallocate (link,listin,lct,lst, Stat=fail(1))
  Deallocate (xxt,yyt,zzt,         Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'four_body_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine four_body_forces
