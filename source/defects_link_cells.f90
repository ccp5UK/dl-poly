Subroutine defects_link_cells &
           (imcon,cell,cut,mxcldef,na,nl,xxt,yyt,zzt,nlx,nly,nlz,link,lct)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for the verlet neighbour list based on link-cell
! method.
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use domains_module

  Implicit None

  Integer,            Intent( In    ) :: imcon,mxcldef,na,nl
  Real( Kind = wp ) , Intent( In    ) :: cell(1:9),cut, &
                                         xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms)

  Integer,            Intent(   Out ) :: nlx,nly,nlz,link(1:mxatms),lct(0:mxcldef)

  Logical, Save :: newjob = .true.

  Logical           :: lx0,lx1,ly0,ly1,lz0,lz1
  Integer           :: icell,ncells,i,ix,iy,iz,jx,jy,jz
  Real( Kind = wp ) :: celprp(1:10), dispx,dispy,dispz, xdc,ydc,zdc

  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Calculate the number of link-cells per domain in every direction

  dispx=celprp(7)/(cut*nprx_r)
  dispy=celprp(8)/(cut*npry_r)
  dispz=celprp(9)/(cut*nprz_r)

  nlx=Int(dispx)
  nly=Int(dispy)
  nlz=Int(dispz)

! check for link cell algorithm violations

  If (nlx*nly*nlz == 0) Call error(307)

  ncells=(nlx+2)*(nly+2)*(nlz+2)
  If (ncells > mxcldef) Then
     Call warning(90,Real(ncells,wp),Real(mxcldef,wp),0.0_wp)
     Call error(392)
  End If

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! Shifts from global to local link-cell space:
! (0,0,0) left-most link-cell on the domain (halo)
! (nlx+1,nly+1,nly+1) right-most
! link-cell on the domain (halo)

  jx=1-nlx*idx
  jy=1-nly*idy
  jz=1-nlz*idz

!***************************************************************
! Note(1): Due to numerical inaccuracy it is possible that some
! domain particles (1,na) may have link-cell space
! coordinates in the halo / at least one coordinate as shown
! (nlx+1,nly+1,nlz+1)^(0,0,0) / as well as halo particles
! (na+1,nl) may have link-cell coordinates in the domain
! / all coordinate as shown (nlx,nly,nlz)^(1,1,1) / or even
! outside the standard one link-cell width, domain-surrounding
! halo / at least one coordinate as shown
! (>nlx+1,>nly+1,>nlz+1)^(<0,<0,<0) /.
!
! Note(2): In SPME, at high accuracy of the ewald summation
! the b-splines may need more positive halo width than the
! standard one of one link-cell (in set_halo_particles it is
! insured that such is supplied).  Such large positive halo may
! lead to EDGE EFFECTs - link-cells constituting the positive
! halo may have larger dimensions than the domain link-cells.
!***************************************************************

! Form linked list
! Initialise link arrays

  Do i=1,nl
     link(i)=0
  End Do
  Do icell=0,ncells
     lct(icell)=0
  End Do

  Do i=nl,na+1,-1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

! Get cell coordinates accordingly

     If (xxt(i) > -half_plus) Then
        ix = Int(xdc*(xxt(i)+0.5_wp)) + jx
     Else
        ix =-Int(xdc*Abs(xxt(i)+0.5_wp)) + jx - 1
     End If
     If (yyt(i) > -half_plus) Then
        iy = Int(ydc*(yyt(i)+0.5_wp)) + jy
     Else
        iy =-Int(ydc*Abs(yyt(i)+0.5_wp)) + jy - 1
     End If
     If (zzt(i) > -half_plus) Then
        iz = Int(zdc*(zzt(i)+0.5_wp)) + jz
     Else
        iz =-Int(zdc*Abs(zzt(i)+0.5_wp)) + jz - 1
     End If

! Correction for halo particles (na+1,nl) of this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the domain only link-cell space

     lx0=(ix == 1)
     lx1=(ix == nlx)
     ly0=(iy == 1)
     ly1=(iy == nly)
     lz0=(iz == 1)
     lz1=(iz == nlz)
     If ( (lx0 .or. lx1) .and. &
          (ly0 .or. ly1) .and. &
          (lz0 .or. lz1) ) Then ! 8 corners of the domain's cube in RS
        If      (lx0) Then
           ix=0
        Else If (lx1) Then
           ix=nlx+1
        Else If (ly0) Then
           iy=0
        Else If (ly1) Then
           iy=nly+1
        Else If (lz0) Then
           iz=0
        Else If (lz1) Then
           iz=nlz+1
        End If
     End If

! Check for residual halo

     lx0=(ix < 0)
     lx1=(ix > nlx+1)
     ly0=(iy < 0)
     ly1=(iy > nly+1)
     lz0=(iz < 0)
     lz1=(iz > nlz+1)
     If ( .not. &
          (lx0 .or. lx1 .or. &
           ly0 .or. ly1 .or. &
           lz0 .or. lz1) ) Then

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and two more link-cells per
! dimension are accounted /coming from the halo/)

        icell = 1 + ix + (nlx + 2)*(iy + (nly + 2)*iz)

     Else

! Put possible residual halo in cell=0

        icell = 0

     End If

! link points to the next in chain or zero if end of chain occurs
! this is the old lct(icell)

     link(i) = lct(icell)

! at the end of the do-loop lct will point to the head of chain
! for this link-cell (update of lct(icell))

     lct(icell) = i

  End Do
  Do i=na,1,-1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

! Get cell coordinates accordingly

     ix = Int(xdc*(xxt(i)+0.5_wp)) + jx
     iy = Int(ydc*(yyt(i)+0.5_wp)) + jy
     iz = Int(zdc*(zzt(i)+0.5_wp)) + jz

! Correction for domain (idnode) only particles (1,na) but due to
! some tiny numerical inaccuracy kicked into its halo link-cell space
! Put all particles in bounded link-cell space: lower and upper
! bounds as 1 <= i_coordinate <= nl_coordinate

     ix = Max( Min( ix , nlx) , 1)
     iy = Max( Min( iy , nly) , 1)
     iz = Max( Min( iz , nlz) , 1)

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and two more link-cells per
! dimension are accounted /coming from the halo/)

     icell = 1 + ix + (nlx + 2)*(iy + (nly + 2)*iz)

! link points to the next in chain or zero if end of chain occurs
! this is the old lct(icell)

     link(i) = lct(icell)

! at the end of the do-loop lct will point to the head of chain
! for this link-cell (update of lct(icell))

     lct(icell) = i

  End Do

End Subroutine defects_link_cells
