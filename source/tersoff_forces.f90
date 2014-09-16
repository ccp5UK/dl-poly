Subroutine tersoff_forces(imcon,rcter,engter,virter,stress)

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
! author    - w.smith  october 2004
! amended   - i.t.todorov september 2014
! amended   - k.galvin september 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz
  Use config_module,  Only : cell,natms,nlast,lfrzn,ltype, &
                             xxx,yyy,zzz,fxx,fyy,fzz
  Use tersoff_module

  Implicit None

  Integer,                             Intent( In    ) :: imcon
  Real( Kind = wp ),                   Intent( In    ) :: rcter
  Real( Kind = wp ),                   Intent(   Out ) :: engter,virter
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: sidex,sidey,sidez,rdr

! flag for undefined potentials NOT NEEDED HERE YET
  Logical           :: lx0,lx1,ly0,ly1,lz0,lz1,flag3 !,safe

  Integer           :: fail(1:7),                      &
                       i,j,k, ii,jj,kk,ll,             &
                       nbx,nby,nbz, ncells,            &
                       ix,iy,iz,icell, jx,jy,jz,jcell, &
                       iatm,jatm,katm, iter,jter,kter, &
                       ijter,ikter,limit

  Real( Kind = wp ) :: dispx,dispy,dispz, xdc,ydc,zdc,                 &
                       rcell(1:9),celprp(1:10),det,                    &
                       sxij,syij,szij,                                 &
                       gk0,gk1,gk2,vk0,vk1,vk2,                        &
                       t1,t2,ppp,bi,ei,ci,di,c1i,c2i,c3i,c4i,c5i,hi,   &
                       ak,bk,xkj,ykj,zkj,gtheta,cost,hmct2,c4exp,      &
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

  fail=0
  Allocate (link(1:mxatms),listin(1:mxatms),lct(1:mxcell),lst(1:mxcell), Stat=fail(1))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),                   Stat=fail(2))
  Allocate (xtf(1:mxlist),ytf(1:mxlist),ztf(1:mxlist),rtf(1:mxlist),     Stat=fail(3))
  Allocate (ert(1:mxlist),eat(1:mxlist),grt(1:mxlist),gat(1:mxlist),     Stat=fail(4))
  Allocate (scr(1:mxlist),gcr(1:mxlist),                                 Stat=fail(5))
  Allocate (cst(1:mxlist),gam(1:mxlist),gvr(1:mxlist),                   Stat=fail(6))
  If (potter == 2) Allocate (rkj(1:mxlist),wkj(1:mxlist),                Stat=fail(7))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'tersoff_forces allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Get this node's (domain's) coordinates

     idz=idnode/(nprx*npry)
     idy=idnode/nprx-idz*npry
     idx=Mod(idnode,nprx)

! Get the domains' dimensions in reduced space
! (domains are geometrically equivalent)

     sidex=1.0_wp/Real(nprx,wp)
     sidey=1.0_wp/Real(npry,wp)
     sidez=1.0_wp/Real(nprz,wp)

! Get reciprocal of interpolation interval

     rdr=Real(mxgter-4,wp)/rcter
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Calculate the number of link-cells per domain in every direction

  nbx=Int(sidex*celprp(7)/(rcter+1.0e-6_wp))
  nby=Int(sidey*celprp(8)/(rcter+1.0e-6_wp))
  nbz=Int(sidez*celprp(9)/(rcter+1.0e-6_wp))

! check for link cell algorithm violations

  If (nbx < 3 .or. nby < 3 .or. nbz < 3) Call error(305)

  ncells=(nbx+4)*(nby+4)*(nbz+4)
  If (ncells > mxcell) Then
     Call warning(90,Real(ncells,wp),Real(mxcell,wp),0.0_wp)
     Call error(78)
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

  dispx=0.5_wp-sidex*(Real(idx,wp)-2.0_wp/Real(nbx,wp))
  dispy=0.5_wp-sidey*(Real(idy,wp)-2.0_wp/Real(nby,wp))
  dispz=0.5_wp-sidez*(Real(idz,wp)-2.0_wp/Real(nbz,wp))

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions (ALL - halo included) from centred
! Cartesian coordinates to reduced space coordinates of
! the left-most link-cell

  Do i=1,nlast
     If (lfrter(ltype(i))) Then
        xxt(i)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)+dispx
        yyt(i)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)+dispy
        zzt(i)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)+dispz
     End If
  End Do

! Form linked list
! Initialise link arrays

  link=0
  Do i=1,ncells
     lct(i)=0
     lst(i)=0
  End Do

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
     If (lfrter(ltype(i))) Then

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

        If (i <= natms) Then
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

  engter=0.0_wp
  virter=0.0_wp

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

                 iatm=listin(ii)
                 iter=lstter(ltype(iatm))

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
                       jter=lstter(ltype(jatm))

! bypass if secondary atom type is not involved in interaction

                       If (jter > 0 .and. jatm /= iatm) Then

                          sxij = xxt(jatm)-xxt(iatm)
                          syij = yyt(jatm)-yyt(iatm)
                          szij = zzt(jatm)-zzt(iatm)

                          xtf(jj)=cell(1)*sxij+cell(4)*syij+cell(7)*szij
                          ytf(jj)=cell(2)*sxij+cell(5)*syij+cell(8)*szij
                          ztf(jj)=cell(3)*sxij+cell(6)*syij+cell(9)*szij

                          rtf(jj)=Sqrt(xtf(jj)**2+ytf(jj)**2+ztf(jj)**2)

                          xtf(jj)=xtf(jj)/rtf(jj)
                          ytf(jj)=ytf(jj)/rtf(jj)
                          ztf(jj)=ztf(jj)/rtf(jj)

! interaction index

                          ijter=(Max(iter,jter)*(Max(iter,jter)-1))/2+Min(iter,jter)

! if pair is within cutoff

                          If (rtf(jj) <= vmbp(1,ijter,1)) Then

! if pair is not frozen

                             If (lfrzn(iatm)*lfrzn(jatm) == 0) Then

                                ll  = Int(rdr*rtf(jj))
                                ppp = rtf(jj)*rdr - Real(ll,wp)

! interpolate screening function

                                vk0 = vmbp(ll,  ijter,1)
                                vk1 = vmbp(ll+1,ijter,1)
                                vk2 = vmbp(ll+2,ijter,1)

                                t1 = vk0 + (vk1 - vk0)*ppp
                                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                                scr(jj) = t1 + (t2-t1)*ppp*0.5_wp

! interpolate derivative of screening function

                                gk0 = gmbp(ll,  ijter,1)
                                gk1 = gmbp(ll+1,ijter,1)
                                gk2 = gmbp(ll+2,ijter,1)

                                t1 = gk0 + (gk1 - gk0)*ppp
                                t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                                gcr(jj) = -(t1 + (t2 - t1)*ppp*0.5_wp) / rtf(jj)

! interpolate repulsive component of energy

                                vk0 = vmbp(ll,  ijter,2)
                                vk1 = vmbp(ll+1,ijter,2)
                                vk2 = vmbp(ll+2,ijter,2)

                                t1 = vk0 + (vk1 - vk0)*ppp
                                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                                ert(jj) = t1 + (t2 - t1)*ppp*0.5_wp

! interpolate derivative of repulsive function

                                gk0 = gmbp(ll,  ijter,2)
                                gk1 = gmbp(ll+1,ijter,2)
                                gk2 = gmbp(ll+2,ijter,2)

                                t1 = gk0 + (gk1 - gk0)*ppp
                                t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                                grt(jj) = -(t1 + (t2 - t1)*ppp*0.5_wp) / rtf(jj)

! interpolate attractive component of energy

                                vk0 = vmbp(ll,  ijter,3)
                                vk1 = vmbp(ll+1,ijter,3)
                                vk2 = vmbp(ll+2,ijter,3)

                                t1 = vk0 + (vk1 - vk0)*ppp
                                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                                eat(jj) = t1 + (t2 - t1)*ppp*0.5_wp

! interpolate derivative of attractive function

                                gk0 = gmbp(ll,  ijter,3)
                                gk1 = gmbp(ll+1,ijter,3)
                                gk2 = gmbp(ll+2,ijter,3)

                                t1 = gk0 + (gk1 - gk0)*ppp
                                t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                                gat(jj) = -(t1 + (t2 - t1)*ppp*0.5_wp) / rtf(jj)

                             End If

                          End If

                       End If

                    End Do

! calculate three body (attractive) terms

! Get parameters for iatm

                    If      (potter == 1) Then ! TERS
                       bi=prmter(7, iter)
                       ei=prmter(8, iter)
                       ci=prmter(9, iter)
                       di=prmter(10,iter)
                       hi=prmter(11,iter)
                    Else If (potter == 2) Then ! KIHS
                       ei =prmter(7, iter)
                       di =prmter(8, iter)
                       c1i=prmter(9, iter)
                       c2i=prmter(10,iter)
                       c3i=prmter(11,iter)
                       c4i=prmter(12,iter)
                       c5i=prmter(13,iter)
                       hi =prmter(14,iter)
                       ak =prmter(15,iter)
                       bk =prmter(16,iter)
                    End If

! bond-angle detection

                    Do jj=1,limit

! index of the first secondary atom

                       jatm=listin(jj)
                       jter=lstter(ltype(jatm))

! bypass if secondary atom type is not involved in interaction

                       If (jter > 0 .and. jatm /= iatm) Then

! interaction index

                          ijter=(Max(iter,jter)*(Max(iter,jter)-1))/2+Min(iter,jter)

! if pair is within cutoff

                          If (rtf(jj) <= vmbp(1,ijter,1)) Then

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
                             If (potter == 2) Then ! KIHS
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
     kter=lstter(ltype(katm))

! bypass if secondary atom type is not involved in interaction

     If (kter > 0 .and. katm /= iatm .and. katm /= jatm) Then

! interaction index

        ikter=(Max(iter,kter)*(Max(iter,kter)-1))/2+Min(iter,kter)

! if pair is within cutoff

        If (rtf(kk) <= vmbp(1,ikter,1)) Then

! only for not fully frozen triplets

           If (lfrzn(iatm)*lfrzn(jatm)*lfrzn(katm) == 0) Then

              flag3 = .true.

! bond-angle is detected, force and virial parameters calculation

              cost = xtf(jj)*xtf(kk)+ytf(jj)*ytf(kk)+ztf(jj)*ztf(kk)
              If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)
              cst(kk) = cost

              If      (potter == 1) Then ! TERS
                 gtheta= 1.0_wp + (ci/di)**2 - ci**2 / (di**2 + (hi-cost)**2)
                 eterm = eterm + gtheta*prmter2(ikter,2)*scr(kk) ! L_{ij}
                 vterm = vterm + gtheta*prmter2(ikter,2)*gcr(kk)*rtf(kk)
! d/dr_k of L_{ij} - angular part as it is used in the virial

                 gam(kk) = gtheta
                 gvr(kk) = 2.0_wp * ci**2 * (hi-cost) / (di**2 + (hi-cost)**2)**2 ! d(gtheta)/sint*d(theta)
              Else If (potter == 2) Then ! KIHS
                 rkj(kk)=rtf(jj)-rtf(kk)
                 wkj(kk)=Exp(ak * rkj(kk)**bk)

                 hmct2=(hi-cost)**2
                 c4exp=c4i*Exp(-c5i*hmct2)

                 gtheta= c1i + c2i*hmct2*(1.0_wp+c4exp)/(c3i+hmct2)
                 eterm = eterm + gtheta*wkj(kk)*scr(kk) ! L_{ij}
                 vterm = vterm + gtheta*wkj(kk)*(gcr(kk) - scr(kk)*ak*bk*rkj(kk)**(bk-1))*rtf(kk)
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

  If      (potter == 1) Then ! TERS
     gam_ij=prmter2(iter,1)
     gamma=0.0_wp
     If (flag3) Then
        gam_ij = prmter2(iter,1)*(1.0_wp+(bi*eterm)**ei)**(-0.5_wp/ei) ! gamma_{ij}
        gamma  = eat(jj) * prmter2(iter,1) * bi*(bi*eterm)**(ei-1.0_wp) * &
              0.5_wp*(1.0_wp+(bi*eterm)**ei)**(-0.5_wp/ei - 1.0_wp) ! -FcFa[d/dr gamma_{ij}]/[d/dr Lij]
     End If
  Else If (potter == 2) Then ! KIHS
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

  If (iatm <= natms) Then

     engter = engter + 0.5_wp*(ert(jj) - gam_ij*eat(jj))    ! energy_ij
     virter = virter + 0.5_wp*(gamma*vterm + gterm*rtf(jj)) ! virial_ij

     fxx(iatm)=fxx(iatm)+0.5_wp*gterm*xtf(jj)
     fyy(iatm)=fyy(iatm)+0.5_wp*gterm*ytf(jj)
     fzz(iatm)=fzz(iatm)+0.5_wp*gterm*ztf(jj)

     strs1 = strs1 - 0.5_wp*gterm*rtf(jj)*xtf(jj)*xtf(jj)
     strs2 = strs2 - 0.5_wp*gterm*rtf(jj)*xtf(jj)*ytf(jj)
     strs3 = strs3 - 0.5_wp*gterm*rtf(jj)*xtf(jj)*ztf(jj)
     strs5 = strs5 - 0.5_wp*gterm*rtf(jj)*ytf(jj)*ytf(jj)
     strs6 = strs6 - 0.5_wp*gterm*rtf(jj)*ytf(jj)*ztf(jj)
     strs9 = strs9 - 0.5_wp*gterm*rtf(jj)*ztf(jj)*ztf(jj)

  End If

  If (jatm <= natms) Then

     fxx(jatm)=fxx(jatm)-0.5_wp*gterm*xtf(jj)
     fyy(jatm)=fyy(jatm)-0.5_wp*gterm*ytf(jj)
     fzz(jatm)=fzz(jatm)-0.5_wp*gterm*ztf(jj)

  End If

! calculate three-body forces contributions

  Do kk=1,limit

! index of the second secondary atom

     katm=listin(kk)
     kter=lstter(ltype(katm))

! bypass if secondary atom type is not involved in interaction

     If (kter > 0 .and. katm /= iatm .and. katm /= jatm) Then

! interaction index

        ikter=(Max(iter,kter)*(Max(iter,kter)-1))/2+Min(iter,kter)

! if pair is within cutoff

        If (rtf(kk) <= vmbp(1,ikter,1)) Then

! only for not fully frozen triplets

           If (lfrzn(iatm)*lfrzn(jatm)*lfrzn(katm) == 0) Then

! term in d/dr_k L_{ij} no derivative of gamma

              cost=cst(kk)

! Counteract the double counting with the 0.5_wp factor

              If      (potter == 1) Then ! TERS
                 gam_dg = 0.5_wp*gamma*prmter2(ikter,2)*scr(kk)*gvr(kk)
                 gam_df = 0.5_wp*gamma*prmter2(ikter,2)*gcr(kk)*gam(kk)

! calculate contribution to atomic forces

                 fxj = gam_dg*(xtf(kk)-xtf(jj)*cost)/rtf(jj) ! contributions to j
                 fyj = gam_dg*(ytf(kk)-ytf(jj)*cost)/rtf(jj)
                 fzj = gam_dg*(ztf(kk)-ztf(jj)*cost)/rtf(jj)

                 fxk = gam_dg*(xtf(jj)-xtf(kk)*cost)/rtf(kk) - gam_df*xtf(kk) ! contributions to k
                 fyk = gam_dg*(ytf(jj)-ytf(kk)*cost)/rtf(kk) - gam_df*ytf(kk)
                 fzk = gam_dg*(ztf(jj)-ztf(kk)*cost)/rtf(kk) - gam_df*ztf(kk)
              Else If (potter == 2) Then ! KIHS
                 gam_dg = 0.5_wp*gamma*scr(kk)*wkj(kk)*gvr(kk)
                 gam_df = 0.5_wp*gamma*gam(kk)*wkj(kk)*gcr(kk)
                 gam_dw = 0.5_wp*gamma*scr(kk)*gam(kk)*wkj(kk)*(ak*bk*rkj(kk)**(bk-1))

! calculate contribution to atomic forces

                 fxj = gam_dg*(xtf(kk)-xtf(jj)*cost)/rtf(jj) - gam_dw*xtf(jj) ! contributions to j
                 fyj = gam_dg*(ytf(kk)-ytf(jj)*cost)/rtf(jj) - gam_dw*ytf(jj)
                 fzj = gam_dg*(ztf(kk)-ztf(jj)*cost)/rtf(jj) - gam_dw*ztf(jj)

                 fxk = gam_dg*(xtf(jj)-xtf(kk)*cost)/rtf(kk) - gam_df*xtf(kk) + gam_dw*xtf(kk)! contributions to k
                 fyk = gam_dg*(ytf(jj)-ytf(kk)*cost)/rtf(kk) - gam_df*ytf(kk) + gam_dw*ytf(kk)
                 fzk = gam_dg*(ztf(jj)-ztf(kk)*cost)/rtf(kk) - gam_df*ztf(kk) + gam_dw*ztf(kk)
              End If

              If (iatm <= natms) Then

                 fxx(iatm)=fxx(iatm)-(fxj+fxk)
                 fyy(iatm)=fyy(iatm)-(fyj+fyk)
                 fzz(iatm)=fzz(iatm)-(fzj+fzk)

! calculate contribution to stress tensor (associated to the head atom)

                 strs1 = strs1 + (fxj*xtf(jj)*rtf(jj) + fxk*xtf(kk)*rtf(kk))
                 strs2 = strs2 + (fxj*ytf(jj)*rtf(jj) + fxk*ytf(kk)*rtf(kk))
                 strs3 = strs3 + (fxj*ztf(jj)*rtf(jj) + fxk*ztf(kk)*rtf(kk))
                 strs5 = strs5 + (fyj*ytf(jj)*rtf(jj) + fyk*ytf(kk)*rtf(kk))
                 strs6 = strs6 + (fyj*ztf(jj)*rtf(jj) + fyk*ztf(kk)*rtf(kk))
                 strs9 = strs9 + (fzj*ztf(jj)*rtf(jj) + fzk*ztf(kk)*rtf(kk))

              End If

              If (jatm <= natms) Then

                 fxx(jatm)=fxx(jatm)+fxj
                 fyy(jatm)=fyy(jatm)+fyj
                 fzz(jatm)=fzz(jatm)+fzj

              End If

              If (katm <= natms) Then

                 fxx(katm)=fxx(katm)+fxk
                 fyy(katm)=fyy(katm)+fyk
                 fzz(katm)=fzz(katm)+fzk

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

  If (mxnode > 1) Then
     buffer(1)=engter
     buffer(2)=virter
     Call gsum(buffer(1:2))
     engter=buffer(1)
     virter=buffer(2)
  End If

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

  Deallocate (link,listin,lct,lst,      Stat=fail(1))
  Deallocate (xxt,yyt,zzt,              Stat=fail(2))
  Deallocate (xtf,ytf,ztf,rtf,          Stat=fail(3))
  Deallocate (ert,eat,grt,gat,          Stat=fail(4))
  Deallocate (scr,gcr,                  Stat=fail(5))
  Deallocate (cst,gam,gvr,              Stat=fail(6))
  If (potter == 2) Deallocate (rkj,wkj, Stat=fail(7))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'tersoff_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine tersoff_forces
