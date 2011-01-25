! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Centre for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

Subroutine link_cell_pairs_helper(       &
             ibegin,iend,nix,niy,niz,    &
             nlx,nly, nlx0e, nly0e,      &
             nlz0e, nlx1s, nly1s, nlz1s, &
             lct_start,at_list,          &
             ncells,nlp,nlp3,nsbcll,megfrz,rcsq)
  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gcheck,gmax
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz
  Use config_module,  Only : cell,natms,nlast,ltg,lfrzn, &
                             xxx,yyy,zzz,lexatm,list

  Implicit None

  Integer, Intent (In) :: nlx,nly, nlx0e, nly0e, nlz0e, nlx1s, &
                          nly1s, nlz1s, ibegin, iend, ncells,  &
                          nlp, nsbcll, nlp3, megfrz

  Integer, Dimension(1:mxatms)  , Intent (In) :: at_list
  Integer, Dimension(1:ncells+1), Intent(In) :: lct_start
  Integer, Dimension(1:nlp3)    , Intent (In) :: nix,niy,niz;

  Real (Kind=wp), Intent (In) :: rcsq

  Logical           :: safe
  Integer :: ix,iy,iz,ic, ix1,ix2,iy1,iy2,iz1,iz2, &
             jx,jy,jz,jc,ipass,nlx0s,nly0s,nlz0s,  &
             nlx1e,nly1e,nlz1e,ibig,               &
             i,ii,j,jj,kk,ll,j_start
  Real (Kind=wp) :: rsq

  ibig=0
  safe=.true.


!$OMP PARALLEL PRIVATE(ix,iy,iz,ic, ix1,ix2,iy1,iy2,iz1,iz2,jx,jy,jz,jc,ipass,&
!$OMP nlx0s,nly0s,nlz0s,nlx1e,nly1e,nlz1e,ibig,i,ii,j,jj,kk,ll,j_start,safe,rsq)
  Do ipass=1,2

! primary loop over domain subcells

!$OMP DO
    Do iz=nlz0e+1,nlz1s-1
       iz1=iz-nlz0e
       iz2=iz-nlz1s
       Do iy=nly0e+1,nly1s-1
          iy1=iy-nly0e
          iy2=iy-nly1s

          Do ix=nlx0e+1,nlx1s-1
             ix1=ix-nlx0e
             ix2=ix-nlx1s

! When ipass = 2 be on the domain's border link-cells

             If ( (ipass == 1) .or.                  &
                  (ix1 >=  1 .and. ix1 <=  nlp) .or. &
                  (ix2 <= -1 .and. ix2 >= -nlp) .or. &
                  (iy1 >=  1 .and. iy1 <=  nlp) .or. &
                  (iy2 <= -1 .and. iy2 >= -nlp) .or. &
                  (iz1 >=  1 .and. iz1 <=  nlp) .or. &
                  (iz2 <= -1 .and. iz2 >= -nlp) ) Then

! index of primary cell

                ic=1+ix+(nlx+2*nlp)*(iy+(nly+2*nlp)*iz)

! loop over primary cell contents

                Do ii=lct_start(ic),lct_start(ic+1)-1

! get domain local particle index

                   i=at_list(ii)

! bypass if primary cell particle > natms

                   If (i <= natms .and. i>=ibegin .and. i<=iend) Then


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

                         If ( (ipass == 1) .or.                     &
                              (jx <= nlx0e) .or. (jx >= nlx1s) .or. &
                              (jy <= nly0e) .or. (jy >= nly1s) .or. &
                              (jz <= nlz0e) .or. (jz >= nlz1s) ) Then

! index of neighbouring cell
                            jc=1+jx+(nlx+2*nlp)*(jy+(nly+2*nlp)*jz)

! get the secondary list particle index in linked-list local description

                            If (jc /= ic) Then

! if j-th subcell is not the i-th subcell
! get head of chain of j-th subcell

                               j_start=lct_start(jc)

                            Else

! if j-th subcell is the same as the i-th subcell
! get next in line in the linked list

                               j_start=ii+1

                            End If

! loop over secondary cell contents

                            Do jj=j_start,lct_start(jc+1)-1

! get domain local particle index

                               j=at_list(jj)

! test for frozen atom pairs (frozen atom pairs DO NOT interact)

                               If (megfrz <= 1 .and. lfrzn(i)*lfrzn(j) == 0) Then

! distance in real space

                                  rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2

! check cutoff criterion (atom pairs MUST BE within the cutoff)
                                  If (rsq <= rcsq) Then

! check for overfloat and add an entry

                                     ll=list(0,i)+1
                                     If (ll > mxlist) Then
                                        ibig=Max(ibig,ll)
                                        safe=.false.
                                     Else
                                        ibig=Max(ibig,ll)
                                        list(0,i)=ll
                                        list(ll,i)=j
                                     End If

! end of if-block for cutoff criterion

                                  End If

! end of if-block on non-frozen atoms

                               End If

! end of loop over secondary cell contents

                            End Do

! end of inner if-block on cells and borders

                         End If

! end of secondary loop over subcells

                      End Do

! end of bypass if primary cell particle > natms

                   End If

! end of loop over primary cell contents
                End Do

! end of outer if-block on cells and borders

             End If

! end of loops over ix,iy,iz

          End Do
       End Do
    End Do

! end of loop over passes

!$OMP END DO
  End Do
!$OMP END PARALLEL
End Subroutine link_cell_pairs_helper

Subroutine link_cell_pairs_remove_exclusions_helper(&
     ibegin,iend)
  Use setup_module
  Use config_module,  Only : natms,ltg,lexatm,list

  Implicit None
  Integer, Intent(In) :: ibegin,iend
  Integer             :: i, ii, kk, ll, jj
  Logical             :: match

!$OMP PARALLEL DO PRIVATE(i,ii,kk,ll,jj)
  Do i=ibegin,iend
     ii=lexatm(0,i)
     If (ii > 0) Then
        kk=1
        ll=list(0,i)
        Do While (kk <= ll)
           jj=ltg(list(kk,i))
           If (match(jj,ii,lexatm(1:ii,i))) Then
              If (kk < ll) list(kk,i)=list(ll,i)
              list(ll,i)=0
              ll=ll-1
              list(0,i)=ll
           Else
              kk=kk+1
           End If
        End Do
     End If
  End Do
!$OMP END PARALLEL DO
End Subroutine link_cell_pairs_remove_exclusions_helper

Subroutine link_cell_pairs(imcon,rcut,lbook,megfrz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for the verlet neighbour list based on link-cell
! method.
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gcheck,gmax
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz
  Use config_module,  Only : cell,natms,nlast,ltg,lfrzn, &
                             xxx,yyy,zzz,lexatm,list

#ifdef COMPILE_CUDA
  Use dl_poly_cuda_module
#endif

  Implicit None

  Logical,            Intent( In    ) :: lbook
  Integer,            Intent( In    ) :: imcon,megfrz
  Real( Kind = wp ) , Intent( In    ) :: rcut

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: sidex,sidey,sidez
  Real( Kind = wp ), Save :: cut,rcsq

  Logical           :: safe,lx0,lx1,ly0,ly1,lz0,lz1,match

  Integer           :: fail(1:4),                            &
                       icell,ncells,ipass,kk,ll,             &
                       ibig,i,ii,j,jj, j_start,              &
                       nlx,nly,nlz,nlp,nlp3,nsbcll,          &
                       nlx0s,nly0s,nlz0s, nlx0e,nly0e,nlz0e, &
                       nlx1s,nly1s,nlz1s, nlx1e,nly1e,nlz1e, &
                       ix,iy,iz,ic, ix1,ix2,iy1,iy2,iz1,iz2, &
                       jx,jy,jz,jc

  Real( Kind = wp ) :: rsq,det,rcell(1:9),celprp(1:10), &
                       dispx,dispy,dispz, xdc,ydc,zdc,nlp2

  Integer,           Dimension( : ), Allocatable :: nix,niy,niz,         &
                                                    lct_count,lct_start, &
                                                    lct_where,which_cell,at_list
  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt

#ifdef COMPILE_CUDA
Call start_timing_link_cell_pairs()
#endif

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
! (domains are geometrically equvalent)

     sidex=1.0_wp/Real(nprx,wp)
     sidey=1.0_wp/Real(npry,wp)
     sidez=1.0_wp/Real(nprz,wp)
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Calculate the number of link-cells per domain in every direction

  dispx=sidex*celprp(7)/cut
  dispy=sidey*celprp(8)/cut
  dispz=sidez*celprp(9)/cut

  nlx=Int(dispx)
  nly=Int(dispy)
  nlz=Int(dispz)

! check for link cell algorithm violations

  If (nlx*nly*nlz == 0) Call error(307)

  ncells=(nlx+2)*(nly+2)*(nlz+2)
  If (ncells > mxcell) Then
     Call warning(90,Real(ncells,wp),Real(mxcell,wp),0.0_wp)
     Call error(392)
  End If

! subcelling and new link-cell parameters

  nlp=1
  nlp2=Real(natms,wp)
  det=nlp2/Real(nlx*nly*nlz,wp)
  Do While (det > 130.0_wp)
     nlp=nlp+1
     rsq=Real(nlp,wp)
     nlx=Int(dispx*rsq)
     nly=Int(dispy*rsq)
     nlz=Int(dispz*rsq)
     det=nlp2/Real(nlx*nly*nlz,wp)
  End Do
  ncells=(nlx+2*nlp)*(nly+2*nlp)*(nlz+2*nlp)
  nlp3=(1+(1+2*nlp)**3)/2

  fail=0
  Allocate (nix(1:nlp3),niy(1:nlp3),niz(1:nlp3),                              Stat=fail(1))
  Allocate (which_cell(1:mxatms),at_list(1:mxatms),                           Stat=fail(2))
  Allocate (lct_count(1:ncells),lct_start(1:ncells+1 ),lct_where(1:ncells+1), Stat=fail(3))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),                        Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs allocation failure, node: ', idnode
     Call error(0)
  End If

! Create jump-around arrays in link-cell space mapping a
! discrete 3D BALL so that ALL two body interactions are
! single counted within a domain (but not over-all as
! double counting still occures globally for all shared
! inter-domain/semi-hello/cross-domain pairs.

  nix=0 ; niy=0 ; niz=0
  nlp2=Real(nlp**2,wp)
  nsbcll=0
  Do iz=0,nlp
     If (iz > 0) Then
        dispz=Real(iz-1,wp)**2
     Else
        dispz=0.0_wp
     End If

     Do iy=-nlp,nlp
        If (iz == 0 .and. iy < 0) Go To 20

        ibig=Abs(iy)
        If (ibig > 0) Then
           dispy=Real(ibig-1,wp)**2
        Else
           dispy=0.0_wp
        End If

        det=dispz+dispy
        If (det > nlp2) Go To 20

        Do ix=-nlp,nlp
           If (iz == 0 .and. iy == 0 .and. ix < 0) Go To 10

           ibig=Abs(ix)
           If (ibig > 0) Then
              dispx=Real(ibig-1,wp)**2
           Else
              dispx=0.0_wp
           End If

           rsq=det+dispx
           If (rsq > nlp2) Go To 10
           nsbcll=nsbcll+1
           nix(nsbcll)=ix
           niy(nsbcll)=iy
           niz(nsbcll)=iz
10         Continue
        End Do
20      Continue
     End Do
  End Do

! Calculate the displacements from the origin of the MD cell
! to the bottom left corner of the left-most halo link-cell

! First term (0.5_wp) = move to the bottom left corner of MD cell
! Second term, first term (side) = scale by the number of domains
! in the given direction
! Second term, second term, first term (id) = move to the bottom
! left corner of this domain in the given direction
! Second term, second term, second term
! (Real(nlp,wp)/Real(nl_coordinate,wp)) =
! move to the bottom left corner of the left-most link-cell
! (the one in the halo)

  dispx=0.5_wp-sidex*(Real(idx,wp)-Real(nlp,wp)/Real(nlx,wp))
  dispy=0.5_wp-sidey*(Real(idy,wp)-Real(nlp,wp)/Real(nly,wp))
  dispz=0.5_wp-sidez*(Real(idz,wp)-Real(nlp,wp)/Real(nlz,wp))

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
! Initilise cell contents counter

  lct_count = 0

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! Move ALL particles in link-cell space:
! (0,0,0) left-most link-cell on the domain (halo)
! (nlx+2*nlp-1,nly+2*nlp-1,nly+2*nlp-1) right-most
! link-cell on the domain (halo)
!***************************************************************
! Note(1): Due to numerical inaccuracy it is possible that some
! domain particles (1,natms) may have like-cell space
! coordinates in the halo / at least one coordinate as shown
! (nlx+1,nly+1,nlz+1)^(nlp-1,nlp-1,nlp-1) /
! as well as halo particles (natms+1,nlast) may have link-cell
! coordinates in the domain / all coordinate as shown
! (nlx,nly,nlz)^(nlp,nlp,nlp) / or even outside the
! standard one link-cell width, domain-surrounding halo
! / at least one coordinate as shown
! (>nlx+2*nlp-1,>nly+2*nlp-1,>nlz+2*nlp-1)^(<0,<0,<0) /.
!
! Note(2): In SPME, at high accuracy of the ewald summation
! the b-splines may need more positive halo width than the
! standard one of one link-cell (in set_halo_particles it is
! insured that such is supplied). Such large positive halo may
! lead to EDGE EFFECTs - link-cells constituting the positive
! halo may have larger dimensions than the domain link-cells.
!***************************************************************

! LC limits
! halo -start

  nlx0s=0
  nly0s=0
  nlz0s=0

! halo -end

  nlx0e=nlp-1
  nly0e=nlp-1
  nlz0e=nlp-1

! halo +start

  nlx1s=nlx+nlp
  nly1s=nly+nlp
  nlz1s=nlz+nlp

! halo +end

  nlx1e=nlx+2*nlp-1
  nly1e=nly+2*nlp-1
  nlz1e=nlz+2*nlp-1

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

! Put all particles in bounded link-cell space: lower and upper bound

     ix = Max( Min( ix , nlx1e) , nlx0s)
     iy = Max( Min( iy , nly1e) , nly0s)
     iz = Max( Min( iz , nlz1e) , nlz0s)

! Correction for particles (1,natms) belonging to this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the halo link-cell space and vice-versa particles from the
! halo (natms+1,nlast) kicked into the domain link-cell space

     If (i <= natms) Then
        If (ix <= nlx0e) ix=nlx0e+1
        If (ix >= nlx1s) ix=nlx1s-1

        If (iy <= nly0e) iy=nly0e+1
        If (iy >= nly1s) iy=nly1s-1

        If (iz <= nlz0e) iz=nlz0e+1
        If (iz >= nlz1s) iz=nlz1s-1
     Else
        lx0=(ix == nlx0e+1)
        lx1=(ix == nlx1s-1)
        ly0=(iy == nly0e+1)
        ly1=(iy == nly1s-1)
        lz0=(iz == nlz0e+1)
        lz1=(iz == nlz1s-1)
        If ((lx0 .or. lx1) .and. (ly0 .or. ly1) .and. (lz0 .or. lz1)) Then
           If      (lx0) Then
              ix=nlx0e
           Else If (lx1) Then
              ix=nlx1s
           Else If (ly0) Then
              iy=nly0e
           Else If (ly1) Then
              iy=nly1s
           Else If (lz0) Then
              iz=nlz0e
           Else If (lz1) Then
              iz=nlz1s
           End If
        End If
     End If

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and 2*nlp more link-cells per
! dimension are accounted /coming from the halo/)

     icell = 1 + ix + (nlx + 2*nlp)*(iy + (nly + 2*nlp)*iz)

! count cell content

     lct_count(icell) = lct_count(icell) + 1

! backwards relationship

     which_cell(i) = icell

  End Do

! break down local list to list of linked-cell lists

  lct_start(1) = 1
  Do i=2,ncells+1
     lct_start(i) = lct_start(i-1) + lct_count(i-1)
  End Do

! domain local to linked-lists local mapping

  lct_where = lct_start
  Do i=1,nlast
     at_list( lct_where( which_cell( i ) ) ) = i
     lct_where( which_cell( i ) ) = lct_where( which_cell( i ) ) + 1
  End Do

! initialise verley neighbourlist arrays

!  list=0      ! (DEBUG)
  list(0,:)=0 !  (FIRST DIMENSION ONLY)

! initial values of control variables

  ibig=0
  safe=.true.


#ifdef COMPILE_CUDA
  if (dl_poly_cuda_offload_link_cell_pairs() .and. dl_poly_cuda_is_cuda_capable()) Then
     Call link_cell_pairs_cuda_initialise(                   &
       natms,mxatms,ncells,mxlist,mxatdm,nsbcll,nlp,nlx,nly, &
       nlx0e, nly0e, nlz0e, nlx1s, nly1s, nlz1s,             &
       xxx,yyy,zzz,lfrzn,at_list,lct_start,list,             &
       nlp3,nix,niy,niz,rcsq,lbook,mxexcl,lexatm,ltg,megfrz)
     Call link_cell_pairs_cuda_invoke()
! 20100129/ck: do not finalise here -- this will keep the list online, i.e.
!   on the card for the rest of the kernels to access without copying it
!   over from the host again.
! 20100218/ck: when needed to finalise, mind if lbook==.true. as the list
!   will be stil lneeded for the exclusions part.
     If ((.not.dl_poly_cuda_offload_tbforces()) .and. (.not.lbook)) Then
        Call link_cell_pairs_cuda_finalise()
     End If

! 20100226/ck: uncomment this to study the effect of sorted atom lists
!   to the CUDA kernels of two_body_forces.
!     If (dl_poly_cuda_offload_tbforces() .and. (.not. lbook)) Then
!        Call link_cell_pairs_cuda_invoke_sort_atoms()
!     End If
  Else
#endif !COMPILE_CUDA

! primary loop over domain subcells

    Do iz=nlz0e+1,nlz1s-1
       iz1=iz-nlz0e
       iz2=iz-nlz1s
       Do iy=nly0e+1,nly1s-1
          iy1=iy-nly0e
          iy2=iy-nly1s

          Do ix=nlx0e+1,nlx1s-1
             ix1=ix-nlx0e
             ix2=ix-nlx1s

! loop over the domain's link-cells only (ipass=1) and
! over the domain's border link-cells only (ipass=2)

           Do ipass=1,2

              If ( (ipass == 1) .or.                  &
                   (ix1 >=  1 .and. ix1 <=  nlp) .or. &
                   (ix2 <= -1 .and. ix2 >= -nlp) .or. &
                   (iy1 >=  1 .and. iy1 <=  nlp) .or. &
                   (iy2 <= -1 .and. iy2 >= -nlp) .or. &
                   (iz1 >=  1 .and. iz1 <=  nlp) .or. &
                   (iz2 <= -1 .and. iz2 >= -nlp) ) Then

! index of primary cell

                 ic=1+ix+(nlx+2*nlp)*(iy+(nly+2*nlp)*iz)

! loop over primary cell contents

                Do ii=lct_start(ic),lct_start(ic+1)-1

! get domain local particle index

                   i=at_list(ii)

! bypass if primary cell particle > natms

                   If (i <= natms) Then

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

! When ipass = 2 be on the halo link-cells only

                         If ( (ipass == 1) .or.                     &
                              (jx <= nlx0e) .or. (jx >= nlx1s) .or. &
                              (jy <= nly0e) .or. (jy >= nly1s) .or. &
                              (jz <= nlz0e) .or. (jz >= nlz1s) ) Then

! index of neighbouring cell

                            jc=1+jx+(nlx+2*nlp)*(jy+(nly+2*nlp)*jz)

! get the secondary list particle index in linked-list local description

                            If (jc /= ic) Then

! if j-th subcell is not the i-th subcell
! get head of chain of j-th subcell

                               j_start=lct_start(jc)

                            Else

! if j-th subcell is the same as the i-th subcell
! get next in line in the linked list

                               j_start=ii+1

                            End If

! loop over secondary cell contents

                            Do jj=j_start,lct_start(jc+1)-1

! get domain local particle index

                               j=at_list(jj)

! test for frozen atom pairs (frozen atom pairs DO NOT interact) is now put outside
!
!                               If (lfrzn(i)*lfrzn(j) == 0) Then

! distance in real space

                                  rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2

! check cutoff criterion (atom pairs MUST BE within the cutoff)

                                  If (rsq <= rcsq) Then

! check for overfloat and add an entry

                                     ll=list(0,i)+1
                                     If (ll > mxlist) Then
                                        ibig=Max(ibig,ll)
                                        safe=.false.
                                     Else
                                        ibig=Max(ibig,ll)
                                        list(0,i)=ll
                                        list(ll,i)=j
                                     End If

! end of if-block for cutoff criterion

                                  End If

! end of if-block on non-frozen atoms
!
!                               End If

! end of loop over secondary cell contents

                            End Do

! end of inner if-block on cells and borders

                         End If

! end of secondary loop over subcells

                      End Do

! end of bypass if primary cell particle > natms

                   End If

! end of loop over primary cell contents

                End Do

! end of outer if-block on cells and borders

             End If

! end of loops over ix,iy,iz

          End Do

! end of loop over passes

       End Do
    End Do
  End Do

#ifdef COMPILE_CUDA
  End If
#endif

! terminate job if neighbour list array exceeded

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Call gmax(ibig)
     Call warning(290,Real(ibig,wp),Real(mxlist,wp),0.0_wp)
     Call error(106)
  End If


! Remove frozen pairs
#ifdef COMPILE_CUDA
  ! In the CUDA port, frozen atom pairs are not included in the
  !  list() to begin with, so they don't need to be removed here
  If (.not.dl_poly_cuda_offload_link_cell_pairs() .and. &
       dl_poly_cuda_is_cuda_capable()) Then
#endif
     If (megfrz > 1) Then
        Do i=1,natms
           ii=lfrzn(i)
           If (ii > 0) Then
              kk=1
              ll=list(0,i)
              Do While (kk <= ll)
                 jj=lfrzn(list(kk,i))
                 If (jj > 0) Then
                    If (kk < ll) list(kk,i)=list(ll,i)
                    list(ll,i)=0
                    ll=ll-1
                    list(0,i)=ll
                 Else
                    kk=kk+1
                 End If
              End Do
           End If
        End Do
     End If
#ifdef COMPILE_CUDA
  End If
#endif


! Remove exluded interactions from the verlet neighbour list


  If (lbook) Then

#ifdef COMPILE_CUDA
     Call start_timing_link_cell_pairs_cuda_remove_excluded()
     If (dl_poly_cuda_offload_link_cell_pairs() .and. &
         dl_poly_cuda_is_cuda_capable()) Then
        Call link_cell_pairs_cuda_invoke_remove_exclusions()
        If (.not.dl_poly_cuda_offload_tbforces()) Then
           Call link_cell_pairs_cuda_finalise()
        End If
     Else
#endif

     Do i=1,natms
        ii=lexatm(0,i)
        If (ii > 0) Then
           kk=1
           ll=list(0,i)
           Do While (kk <= ll)
              jj=ltg(list(kk,i))
              If (match(jj,ii,lexatm(1:ii,i))) Then
                 If (kk < ll) list(kk,i)=list(ll,i)
                 list(ll,i)=0
                 ll=ll-1
                 list(0,i)=ll
              Else
                 kk=kk+1
              End If
           End Do
        End If
     End Do

#ifdef COMPILE_CUDA
     End If
  Call stop_timing_link_cell_pairs_cuda_remove_excluded()
#endif

  End If

  Deallocate (nix,niy,niz,                   Stat=fail(1))
  Deallocate (which_cell,at_list,            Stat=fail(2))
  Deallocate (lct_count,lct_start,lct_where, Stat=fail(3))
  Deallocate (xxt,yyt,zzt,                   Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs deallocation failure, node: ', idnode
     Call error(0)
  End If

#ifdef COMPILE_CUDA
  Call stop_timing_link_cell_pairs()
#endif

End Subroutine link_cell_pairs
