Subroutine defects_link_cells &
           (imcon,cell,cut,mxcldef,na,nl,xxx,yyy,zzz,nlx,nly,nlz,link,lct)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for the verlet neighbour list based on link-cell
! method.
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz

  Implicit None

  Integer,            Intent( In    ) :: imcon,mxcldef,na,nl
  Real( Kind = wp ) , Intent( In    ) :: cell(1:9),cut, &
                                         xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)

  Integer,            Intent(   Out ) :: nlx,nly,nlz,link(1:mxatms),lct(1:mxcldef)

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: sidex,sidey,sidez

  Logical           :: lx0,lx1,ly0,ly1,lz0,lz1

  Integer           :: fail(1:2),icell,ncells,i,ix,iy,iz

  Real( Kind = wp ) :: celprp(1:10),dispx,dispy,dispz, xdc,ydc,zdc

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt

  fail=0
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_link_cells allocation failure, node: ', idnode
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
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Calculate the number of link-cells per domain in every direction

  nlx=Int(sidex*celprp(7)/cut)
  nly=Int(sidey*celprp(8)/cut)
  nlz=Int(sidez*celprp(9)/cut)

! check for link cell algorithm violations

  If (nlx*nly*nlz == 0) Call error(307)

  ncells=(nlx+2)*(nly+2)*(nlz+2)
  If (ncells > mxcldef) Then
     Call warning(90,Real(ncells,wp),Real(mxcldef,wp),0.0_wp)
     Call error(392)
  End If

! Calculate the displacements from the origin of the MD cell
! to the bottom left corner of the left-most halo link-cell

! First term (0.5_wp) = move to the bottom left corner of MD cell
! Second term, first term (side) = scale by the number of domains
! in the given direction
! Second term, second term, first term (id) = move to the bottom
! left corner of this domain in the given direction
! Second term, second term, second term (1.0_wp/Real(nl,wp)) =
! move to the bottom left corner of the left-most link-cell
! (the one that constructs the halo)

  dispx=0.5_wp-sidex*(Real(idx,wp)-1.0_wp/Real(nlx,wp))
  dispy=0.5_wp-sidey*(Real(idy,wp)-1.0_wp/Real(nly,wp))
  dispz=0.5_wp-sidez*(Real(idz,wp)-1.0_wp/Real(nlz,wp))

! Convert atomic positions (ALL - halo included) from centred
! Cartesian coordinates to reduced space coordinates of
! the left-most link-cell

  Do i=1,nl
     xxt(i)=xxx(i)+dispx
     yyt(i)=yyy(i)+dispy
     zzt(i)=zzz(i)+dispz
  End Do

! Form linked list
! Initialise link arrays

  link=0
  Do i=1,ncells
     lct(i)=0
  End Do

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! Move ALL particles in link-cell space:
! (0,0,0) left-most link-cell on the domain (halo)
! (nlx+1,nly+1,nly+1) right-most link-cell on the domain (halo)
!***************************************************************
! Note: Due to numerical inaccuracy it is possible that some
! domain particles (1,na) may have like-cell space
! coordinates in the halo / at least one coordinate as shown
! (nlx+1,nly+1,nlz+1)^(0,0,0) / as well as halo particles
! (na+1,nl) may have link-cell coordinates in the domain
! / all coordinate as shown (nlx,nly,nlz)^(1,1,1) /.  No
! EDGE EFFECT exists since all particles (na+1,nl) that are
! outside the link-cell layer are not considered.
!***************************************************************

  Do i=1,nl

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

! Correction for particles (1,na) belonging to this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the halo link-cell space and vice-versa particles from
! the halo (na+1,nl) kicked into the domain link-cell space

     If (i <= na) Then
        If (ix == 0)     ix=1
        If (ix == nlx+1) ix=nlx

        If (iy == 0)     iy=1
        If (iy == nly+1) iy=nly

        If (iz == 0)     iz=1
        If (iz == nlz+1) iz=nlz
     Else
        lx0=(ix == 1)
        lx1=(ix == nlx)
        ly0=(iy == 1)
        ly1=(iy == nly)
        lz0=(iz == 1)
        lz1=(iz == nlz)
        If ((lx0 .or. lx1) .and. (ly0 .or. ly1) .and. (lz0 .or. lz1)) Then
           If      (lx0) Then
              ix=0
           Else If (lx1) Then
              ix=nlx
           Else If (ly0) Then
              iy=0
           Else If (ly1) Then
              iy=nly
           Else If (lz0) Then
              iz=0
           Else If (lz1) Then
              iz=nlz
           End If
        End If
     End If

! Only for particles onto the domain and in the link-cell width layer
! around it (since we're in unbounded space naw).  Discard the rest(-:

     If ( (ix >= 0 .and. ix <= nlx+1) .and. &
          (iy >= 0 .and. iy <= nly+1) .and. &
          (iz >= 0 .and. iz <= nlz+1) ) Then

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and two more link-cells per
! dimension are accounted /coming from the halo/)

        icell=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! link points to the next in chain or zero if end of chain occurs
! this is the old lct(icell)

        link(i)=lct(icell)

! at the end of the do-loop lct will point the head of chain
! for this link-cell (update of lct(icell))

        lct(icell)=i
     End If
  End Do

  Deallocate (xxt,yyt,zzt, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_link_cells deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine defects_link_cells
