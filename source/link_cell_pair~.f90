Subroutine link_cell_pairs(imcon,rcut,lbook,megfrz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for the verlet neighbour list based on link-cell
! method.
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov january 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gcheck,gmax
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz
  Use config_module,  Only : cell,natms,nlast,ltg,lfrzn, &
                             xxx,yyy,zzz,lexatm,list

  Implicit None

  Logical,            Intent( In    ) :: lbook
  Integer,            Intent( In    ) :: imcon,megfrz
  Real( Kind = wp ) , Intent( In    ) :: rcut

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: sidex,sidey,sidez
  Real( Kind = wp ), Save :: cut,rcsq

  Logical           :: linc,safe,lx0,lx1,ly0,ly1,lz0,lz1

  Integer           :: fail(1:2),                &
                       icell,ncells,ipass,kk,ll, &
                       ibig,i,j, nlx,nly,nlz,    &
                       ix,iy,iz,ic, jx,jy,jz,jc

  Real( Kind = wp ) :: rsq,det,rcell(1:9),celprp(1:10), &
                       dispx,dispy,dispz, xdc,ydc,zdc

! Number of neighbouring cells to look around for single counting
! of two body interactions

  Integer, Parameter :: nsbcll = 14

! Direction arrays for jumping around in link-cell space so that
! ALL two body interactions are single counted within a domain
! (but not over-all as double counting still occurs at the
! border for neghbouring link-cells onto neghbouring domains).

  Integer, Dimension( 1:nsbcll ), Parameter :: &
  nix = (/ 0, 1, -1, 0, 1,  -1, 0, 1, -1, 0, 1, -1, 0, 1 /) , &
  niy = (/ 0, 0,  1, 1, 1,  -1,-1,-1,  0, 0, 0,  1, 1, 1 /) , &
  niz = (/ 0, 0,  0, 0, 0,   1, 1, 1,  1, 1, 1,  1, 1, 1 /)

  Integer,           Dimension( : ), Allocatable :: link,lct
  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt

  fail=0
  Allocate (link(1:mxatms),lct(1:mxcell),              Stat=fail(1))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Real space cutoff and squared r.s.c.

     cut=rcut+1.0e-6_wp
     rcsq=rcut**2

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

! halt program if potential cutoff exceeds the minimum half-cell width

  det=Min(celprp(7),celprp(8),celprp(9))
  If (rcut > det/2.0_wp) Then
     Call warning(3,rcut,det/2.0_wp,0.0_wp)
     Call error(95)
  End If

! Calculate the number of link-cells per domain in every direction

  nlx=Int(sidex*celprp(7)/cut)
  nly=Int(sidey*celprp(8)/cut)
  nlz=Int(sidez*celprp(9)/cut)

! check for link cell algorithm violations

  If (nlx*nly*nlz == 0) Call error(307)

  ncells=(nlx+2)*(nly+2)*(nlz+2)
  If (ncells > mxcell) Then
     Call warning(90,Real(ncells,wp),Real(mxcell,wp),0.0_wp)
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

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions (ALL - halo included) from centred
! Cartesian coordinates to reduced space coordinates of
! the left-most link-cell

  Do i=1,nlast
     xxt(i)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)+dispx
     yyt(i)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)+dispy
     zzt(i)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)+dispz
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
! Note(1): Due to numerical inaccuracy it is possible that some
! domain particles (1,natms) may have link-cell space
! coordinates in the halo / at least one coordinate as shown
! (nlx+1,nly+1,nlz+1)^(0,0,0) / as well as halo particles
! (natms+1,nlast) may have link-cell coordinates in the domain
! / all coordinate as shown (nlx,nly,nlz)^(1,1,1) / or even
! outside the standard one link-cell width, domain-surrounding
! halo / at least one coordinate as shown
! (>nlx+1,>nly+1,>nlz+1)^(<0,<0,<0) /.
!
! Note(2): In SPME, at high accuracy of the ewald summation
! the b-splines may need more positive halo width than the
! standard one of one link-cell (in set_halo_particles it is
! insured that such is supplied). Such large positive halo may
! lead to EDGE EFFECTs - link-cells constituting the positive
! halo may have larger dimensions than the domain link-cells.
!***************************************************************

  Do i=1,nlast

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

! Put all particles in bounded link-cell space: lower and upper
! bounds as 0 <= i_coordinate <= nl_coordinate+1

     ix = Max( Min( ix , nlx+1) , 0)
     iy = Max( Min( iy , nly+1) , 0)
     iz = Max( Min( iz , nlz+1) , 0)

! Correction for particles (1,natms) belonging to this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the halo link-cell space and vice-versa particles from the
! halo (natms+1,nlast) kicked into the domain link-cell space

     If (i <= natms) Then
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

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and two more link-cells per
! dimension are accounted /coming from the halo/)

     icell=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! link points to the next in chain or zero if end of chain occurs
! this is the old lct(icell)

     link(i)=lct(icell)

! at the end of the do-loop lct will point to the head of chain
! for this link-cell (update of lct(icell))

     lct(icell)=i

  End Do

! initialise verlet neighbourlist arrays

!  list=0         ! (DEBUG)
  list(-2:0,:)=0 ! (COUNTING DIMENSIONS ONLY)

! initial values of control variables

  ibig=0
  safe=.true.

! loop over the domain's link-cells only (ipass=1) and
! over the domain's border link-cells only (ipass=2)

  Do ipass=1,2

! primary loop over domain subcells

    Do iz=1,nlz
       Do iy=1,nly
          Do ix=1,nlx

! When ipass = 2 be on the domain's border link-cells

             If ( (ipass == 1) .or. (ix == 1) .or. (ix == nlx) &
                               .or. (iy == 1) .or. (iy == nly) &
                               .or. (iz == 1) .or. (iz == nlz) ) Then

! index of primary cell

                ic=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! bypass cell if empty

                If (lct(ic) > 0) Then

! secondary loop over subcells

                   Do kk=ipass,nsbcll

                      If (ipass == 1) Then
                         jx=ix+nix(kk)
                         jy=iy+niy(kk)
                         jz=iz+niz(kk)
                      Else
                         jx=ix-nix(kk)
                         jy=iy-niy(kk)
                         jz=iz-niz(kk)
                      End If

! When ipass = 2 be on the halo link-cells

                      If ( (ipass == 1) .or. (jx == 0) .or. (jx == nlx+1) &
                                        .or. (jy == 0) .or. (jy == nly+1) &
                                        .or. (jz == 0) .or. (jz == nlz+1) ) Then

! index of neighbouring cell

                         jc=1+jx+(nlx+2)*(jy+(nly+2)*jz)

! bypass cell if empty

                         If (lct(jc) > 0) Then

! head of chain of i-th subcell

                            i=lct(ic)

! loop over primary cell contents

100                         Continue

! bypass if primary cell particle > natms

                            If (i <= natms) Then


! get the paired particle to compare i to

                               If (jc /= ic) Then

! if j-th subcell is not the i-th subcell
! get head of chain of j-th subcell

                                  j=lct(jc)

                               Else

! if j-th subcell is the same as the i-th subcell
! get next in line in link

                                  j=link(i)

                               End If

! bypass cell if empty

                               If (j > 0) Then

! loop over secondary cell contents

200                               Continue

! test for frozen atom pairs (frozen atom pairs DO NOT interact)

                                  If (lfrzn(i)*lfrzn(j) == 0) Then

! distance in real space

                                     rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2

! check cutoff criterion (atom pairs MUST BE within the cutoff)

                                     If (rsq <= rcsq) Then

! test for excluded atom (by construction a domain particle has index i <= natms)

                                        linc=.true.
                                        If (lbook) Then
                                           Do ll=1,lexatm(0,i)
                                              If (lexatm(ll,i) == ltg(j)) linc=.false.
                                           End Do
                                        End If

! check for overfloat and add an entry

                                        ll=Max(list(0,i),list(-1,i))+1
                                        If (ll > mxlist) Then
                                           ibig=Max(ibig,ll)
                                           safe=.false.
                                        Else
                                           If (linc) Then
                                              list(0,i)=list(0,i)+1
                                              If (list(-1,i) > 0) Then ! roll
                                                 list(-1,i)=list(-1,i)+1
                                                 list(list(-1,i),i)=list(list(0,i),i)
                                              End If
                                              list(list(0,i),i)=j
                                           Else
                                              If      (list(-1,i) == 0) Then
                                                 list(-1,i)=list(0,i)+1
                                                 list(list(-1,i),i)=j
                                              Else If (list(-1,i) >  0) Then
                                                 list(-1,i)=list(-1,i)+1
                                                 list(list(-1,i),i)=j
                                              Else !If (list(-1,i) <  0) Then
                                                 safe=.false.
                                              End If
                                           End If
                                        End If

! end of if-block for cutoff criterion

                                     End If

! end of if-block on non-frozen atoms

                                  End If

! new secondary cell particle

                                  j=link(j)
                                  If (j /= 0) Go To 200

! end of loop over secondary cell contents

                               End If

! end of bypass if primary cell particle > natms

                            End If

! new primary cell particle

                            i=link(i)
                            If (i /= 0) Go To 100

! end of bypass of empty subcell jc

                         End If

! end of inner if-block on cells and borders

                      End If

! end of secondary loop over subcells

                   End Do

! end of bypass of empty subcell ic

                End If

! end of outer if-block on cells and borders

             End If

! end of loops over ix,iy,iz

          End Do
       End Do
    End Do

! end of loop over passes

  End Do

  If (megfrz > 1) Then
     Do i=1,natms
        If (list(-2,i) == 0) list(-2,i)=list(0,i) ! End of full non-repeatable half VNL (FNRH VNL)
     End Do
  Else
     Do i=1,natms
        list(-2,i)=list(0,i) ! End of full non-repeatable half VNL (FNRH VNL)
     End Do
  End If

  If (lbook) Then
     Do i=1,natms
        If (list(-1,i) == 0) list(-1,i)=list(0,i) ! End of NFP FNRH VNL
     End Do
  Else
     Do i=1,natms
        list(-1,i)=list(0,i) ! End of NFP FNRH VNL
     End Do
  End If

! terminate job if neighbour list array exceeded

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Call gmax(ibig)
     Call warning(290,Real(ibig,wp),Real(mxlist,wp),0.0_wp)
     Call error(106)
  End If

  Deallocate (link,lct,    Stat=fail(1))
  Deallocate (xxt,yyt,zzt, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine link_cell_pairs
