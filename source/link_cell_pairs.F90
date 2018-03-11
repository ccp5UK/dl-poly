Subroutine link_cell_pairs(rcut,rlnk,rvdw,rmet,pdplnc,lbook,megfrz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for the verlet neighbour list based on link-cell
! method.
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2017
! contrib   - i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, li
  Use comms_module,       Only : idnode,mxnode,gcheck,gmax,gsum
  Use setup_module
  Use domains_module,     Only : idx,idy,idz, nprx,npry,nprz, &
                                 r_nprx,r_npry,r_nprz
  Use configuration,      Only : cell,natms,nlast,ltg,lfrzn, &
                                 xxx,yyy,zzz,lexatm,list
  Use core_shell,  Only : listshl,legshl
  Use mpoles_module,      Only : keyind,lchatm
  Use development_module, Only : l_dis,r_dis

  Implicit None

  Logical,            Intent( In    ) :: lbook
  Integer,            Intent( In    ) :: megfrz
  Real( Kind = wp ) , Intent( In    ) :: rcut,rlnk,rvdw,rmet,pdplnc

  Logical           :: safe,lx0,lx1,ly0,ly1,lz0,lz1,match

  Integer           :: fail(1:5),l_end,m_end,                 &
                       icell,ncells,ipass,                    &
                       kk,ll, ibig,i,ii,j,jj, j_start,        &
                       nlx,nly,nlz,nlp,nlp2,nlp3,nlp4,nsbcll, &
                       nlx0s,nly0s,nlz0s, nlx0e,nly0e,nlz0e,  &
                       nlx1s,nly1s,nlz1s, nlx1e,nly1e,nlz1e,  &
                       ix,iy,iz,ic, ix1,ix2,iy1,iy2,iz1,iz2,  &
                       jx,jy,jz,jc

  Real( Kind = wp ) :: cut,rcsq,rsq,det,rcell(1:9),celprp(1:10),cnt(0:4), &
                       x,y,z, x1,y1,z1, dispx,dispy,dispz, xdc,ydc,zdc, nlr2

  Logical,           Dimension( : ), Allocatable :: nir
  Integer,           Dimension( : ), Allocatable :: nix,niy,niz,          &
                                                    lct_count,lct_start,  &
                                                    lct_where,which_cell, &
                                                    cell_dom,cell_bor,at_list
  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt


! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! halt program if potential cutoff exceeds the minimum half-cell width

  det=Min(celprp(7),celprp(8),celprp(9))
  If (rlnk > det/2.0_wp) Then
     Call warning(3,rlnk,det/2.0_wp,0.0_wp)
     Call error(95)
  End If

! Real space cutoff and squared r.s.c.

  cut=rlnk+1.0e-6_wp
  rcsq=rlnk**2

! Calculate the number of link-cells per domain in every direction

  dispx=r_nprx*celprp(7)/cut
  dispy=r_npry*celprp(8)/cut
  dispz=r_nprz*celprp(9)/cut

  nlx=Int(dispx)
  nly=Int(dispy)
  nlz=Int(dispz)

! check for link cell algorithm violations

  If (nlx*nly*nlz == 0) Call error(307)

  ncells=(nlx+2)*(nly+2)*(nlz+2)
  If (ncells > mxcell) Then
     Call warning(90,Real(ncells,wp),Real(mxcell,wp),0.0_wp)
     mxcell = Nint(1.25_wp*Real(ncells,wp))
     If (ncells > mxatms) Call error(69)
  End If

! subcelling and new link-cell parameters

  nlp=1
  nlr2=Real(natms,wp)
  det=nlr2/Real(nlx*nly*nlz,wp)
  Do While (det > pdplnc)
     nlp=nlp+1
     rsq=Real(nlp,wp)
     nlx=Int(dispx*rsq)
     nly=Int(dispy*rsq)
     nlz=Int(dispz*rsq)
     det=nlr2/Real(nlx*nly*nlz,wp)
  End Do
  ncells=(nlx+2*nlp)*(nly+2*nlp)*(nlz+2*nlp)
  nlp2=(1+(1+2*nlp)**3)/2                            ! semi-ball size
  nlp3=nlx*nly*nlz                                   ! domain size
  nlp4=nlp3 -                                      & ! border size
              Merge(nlx-2*nlp, 0, nlx-2*nlp > 0) * &
              Merge(nly-2*nlp, 0, nly-2*nlp > 0) * &
              Merge(nlz-2*nlp, 0, nlz-2*nlp > 0)

  fail=0
  Allocate (nix(1:nlp2),niy(1:nlp2),niz(1:nlp2),nir(1:nlp2),                 Stat=fail(1))
  Allocate (which_cell(1:mxatms),at_list(1:mxatms),                          Stat=fail(2))
  Allocate (lct_count(0:ncells),lct_start(0:ncells+1),lct_where(0:ncells+1), Stat=fail(3))
  Allocate (cell_dom(0:nlp3),cell_bor(0:nlp4),                               Stat=fail(4))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),                       Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs allocation failure, node: ', idnode
     Call error(0)
  End If
  cell_dom(0)=nlp3 ! save array's limit
  cell_bor(0)=nlp4 ! save array's limit

! Create jump-around arrays in link-cell space mapping a
! discrete 3D BALL so that ALL two body interactions are
! single counted within a domain (but not over-all as
! double counting still occurs globally for all shared
! inter-domain/semi-hello/cross-domain pairs.


  nix=0 ; niy=0 ; niz=0 ; nir=.false.
  nlp2=nlp**2
  nlp3=(nlp-1)**2
  nsbcll=0
  Do iz=0,nlp
     If (iz > 0) Then
        iz1=(iz-1)**2
     Else
        iz1=0
     End If

     jz=iz**2

     Do iy=-nlp,nlp
        If (iz == 0 .and. iy < 0) Cycle

        nlp4=Abs(iy)
        If (nlp4 > 0) Then
           iy1=(nlp4-1)**2
        Else
           iy1=0
        End If

        ll=iz1+iy1
        If (ll > nlp2) Cycle

        jy=jz+iy**2

        Do ix=-nlp,nlp
           If (iz == 0 .and. iy == 0 .and. ix < 0) Cycle

           nlp4=Abs(ix)
           If (nlp4 > 0) Then
              ix1=(nlp4-1)**2
           Else
              ix1=0
           End If

           If (ll+ix1 > nlp2) Cycle

           jx=jy+ix**2

           nsbcll=nsbcll+1

           nix(nsbcll)=ix
           niy(nsbcll)=iy
           niz(nsbcll)=iz
           nir(nsbcll)=(jx < nlp3)
        End Do
     End Do
  End Do
!  Write(*,*) 'NLP',nlp,nsbcll,nlx,nly,nlz

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! Shifts from global to local link-cell space:
! (0,0,0) left-most link-cell on the domain (halo)
! (nlx+2*nlp-1,nly+2*nlp-1,nly+2*nlp-1) right-most
! link-cell on the domain (halo)

  jx=nlp-nlx*idx
  jy=nlp-nly*idy
  jz=nlp-nlz*idz

!***************************************************************
! Note(1): Due to numerical inaccuracy it is possible that some
! domain particles (1,natms) may have link-cell space
! coordinates in the halo / at least one coordinate as shown
! (nlx+nlp,nly+nlp,nlz+nlp)^(nlp-1,nlp-1,nlp-1) /
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
! ensured that such is supplied).  Such large positive halo may
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

! Form linked list
! Initialise cell contents counter

  lct_count = 0

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

  Do i=1,natms

! Convert atomic positions from MD cell centred
! Cartesian coordinates to reduced space coordinates

     x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! Get cell coordinates accordingly

     ix = Int(xdc*(x+0.5_wp)) + jx
     iy = Int(ydc*(y+0.5_wp)) + jy
     iz = Int(zdc*(z+0.5_wp)) + jz

! Correction for domain (idnode) only particles (1,natms) but due to
! some tiny numerical inaccuracy kicked into its halo link-cell space
! Put all particles in a bounded link-cell space: lower and upper bounds
! as follows nl_coordinate_0e+1 <= i_coordinate <= nl_coordinate_1s-1 !

     ix = Max( Min( ix , nlx1s-1) , nlx0e+1)
     iy = Max( Min( iy , nly1s-1) , nly0e+1)
     iz = Max( Min( iz , nlz1s-1) , nlz0e+1)

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and 2*nlp more link-cells per
! dimension are accounted /coming from the halo/)

     icell = 1 + ix + (nlx + 2*nlp)*(iy + (nly + 2*nlp)*iz)

! count cell content

     lct_count(icell) = lct_count(icell) + 1

! backwards relationship

     which_cell(i) = icell

  End Do
  Do i=natms+1,nlast

! Convert atomic positions from MD cell centred
! Cartesian coordinates to reduced space coordinates

     x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! Get cell coordinates accordingly

     If (x > -half_plus) Then
        dispx=xdc*(x+0.5_wp)
        ix = Int(dispx) + jx
     Else
        dispx=xdc*Abs(x+0.5_wp)
        ix =-Int(dispx) + jx - 1
     End If
     If (y > -half_plus) Then
        dispy=ydc*(y+0.5_wp)
        iy = Int(dispy) + jy
     Else
        dispy=ydc*Abs(y+0.5_wp)
        iy =-Int(dispy) + jy - 1
     End If
     If (z > -half_plus) Then
        dispz=zdc*(z+0.5_wp)
        iz = Int(dispz) + jz
     Else
        dispz=zdc*Abs(z+0.5_wp)
        iz =-Int(dispz) + jz - 1
     End If

! Exclude any negatively bound residual halo

     If (ix >= nlx0s .and. iy >= nly0s .and. iz >= nlz0s) Then

! Correction for halo particles (natms+1,nlast) of this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the domain only link-cell space

        lx0=(ix > nlx0e)
        lx1=(ix < nlx1s)
        ly0=(iy > nly0e)
        ly1=(iy < nly1s)
        lz0=(iz > nlz0e)
        lz1=(iz < nlz1s)
        If ( (lx0 .and. lx1) .and. &
             (ly0 .and. ly1) .and. &
             (lz0 .and. lz1) ) Then

! Put the closest to the halo coordinate in the halo

           x1=Abs(x-0.5_wp*Sign(1.0_wp,x))
           y1=Abs(y-0.5_wp*Sign(1.0_wp,y))
           z1=Abs(z-0.5_wp*Sign(1.0_wp,z))
           If      (x1 <= y1 .and. x1 <= z1) Then
              If (x < 0.0_wp) Then
                 ix=nlx0e
              Else
                 ix=nlx1s
              End If
           Else If (y1 <= x1 .and. y1 <= z1) Then
              If (y < 0.0_wp) Then
                 iy=nly0e
              Else
                 iy=nly1s
              End If
           Else
              If (z < 0.0_wp) Then
                 iz=nlz0e
              Else
                 iz=nlz1s
              End If
           End If
        End If

! Check for positively bound residual halo

        lx0=(ix < nlx0s)
        lx1=(ix > nlx1e)
        ly0=(iy < nly0s)
        ly1=(iy > nly1e)
        lz0=(iz < nlz0s)
        lz1=(iz > nlz1e)
        If ( .not. &
             (lx0 .or. lx1 .or. &
              ly0 .or. ly1 .or. &
              lz0 .or. lz1) ) Then

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/ and 2*nlp more link-cells per
! dimension are accounted /coming from the halo/)

           icell = 1 + ix + (nlx + 2*nlp)*(iy + (nly + 2*nlp)*iz)

        Else

! Put possible residual halo in cell=0

           icell = 0

        End If

     Else

! Put possible residual halo in cell=0

        icell = 0

     End If

! count cell content

     lct_count(icell) = lct_count(icell) + 1

! backwards relationship

     which_cell(i) = icell

  End Do

! break down local list to list of linked-cell lists

  lct_start(0) = 1
  Do icell=1,ncells+1
     lct_start(icell) = lct_start(icell-1) + lct_count(icell-1)
  End Do

! domain local to linked-lists local mapping

  lct_where = lct_start
  Do i=1,nlast
!     at_list( lct_where( which_cell( i ) ) ) = i
! create a reordered coordinates arrays in the
! same manner to speeds up performance later

     j = lct_where( which_cell( i ) )
     at_list( j ) = i

     xxt(j) = xxx(i)
     yyt(j) = yyy(i)
     zzt(j) = zzz(i)

     lct_where( which_cell( i ) ) = lct_where( which_cell( i ) ) + 1
  End Do

  nlp3=0
  nlp4=0
  Do iz=nlz0e+1,nlz1s-1
     iz1=iz-nlz0e
     iz2=iz-nlz1s
     Do iy=nly0e+1,nly1s-1
        iy1=iy-nly0e
        iy2=iy-nly1s
        Do ix=nlx0e+1,nlx1s-1
           ix1=ix-nlx0e
           ix2=ix-nlx1s

! index of the primary cell

           ic=1+ix+(nlx+2*nlp)*(iy+(nly+2*nlp)*iz)

           nlp3=nlp3+1
           cell_dom(nlp3)=ic

! loop over the domain's border cells only - ipass==2

           If ( (ix1 >=  1 .and. ix1 <=  nlp) .or. &
                (ix2 <= -1 .and. ix2 >= -nlp) .or. &
                (iy1 >=  1 .and. iy1 <=  nlp) .or. &
                (iy2 <= -1 .and. iy2 >= -nlp) .or. &
                (iz1 >=  1 .and. iz1 <=  nlp) .or. &
                (iz2 <= -1 .and. iz2 >= -nlp) ) Then
              nlp4=nlp4+1
              cell_bor(nlp4)=ic
           End If
         End Do
      End Do
   End Do

! initialise verlet neighbourlist (VNL) arrays

!  list=0               ! (DEBUG)
  list(-2:0,1:natms)=0 ! (COUNTING DIMENSIONS ONLY)

! initial values of control variables

  ibig=0
  safe=.true.

! primary loop over domain subcells
! loop over the domain's cells only

  ipass=1
  Do icell=1,cell_dom(0)

! index of the primary cell and its coordinates

     ic=cell_dom(icell)
     ix=Mod(ic-1,nlx + 2*nlp)
     iz=(ic-1)/((nlx + 2*nlp)*(nly + 2*nlp))
     iy=(ic-1)/(nlx + 2*nlp) - (nly + 2*nlp)*iz

! loop over the primary cell contents

     Do ii=lct_start(ic),lct_start(ic+1)-1
        i=at_list(ii) ! get the particle index, by construction [1,natms]

! secondary loop over neighbouring cells
! ipass == 1 - non-negative non-repeatable semi-ball

        Do kk=ipass,nsbcll

! be on domain + possible positive halo cells only - ipass==1

           jx=ix+nix(kk)
           jy=iy+niy(kk)
           jz=iz+niz(kk)

! index of the secondary cell

           jc=1+jx+(nlx+2*nlp)*(jy+(nly+2*nlp)*jz)

! get the secondary list's starting particle index

           If (jc /= ic) Then

! if the secondary cell is different from the primary cell
! get the head of chain of the secondary cell

              j_start=lct_start(jc)

           Else ! only when ipass==1

! if the secondary cell is same as the primary cell
! get the next in line from the primary cell running index

              j_start=ii+1

           End If

! bypass on real space cutoff check for safe cells when
! atom pairs' distances are guaranteed to be within the cutoff

           If (nir(kk)) Then

! check for overflow

              ll=list(0,i)+lct_start(jc+1)-j_start
              If (ll <= mxlist) Then

! loop over the secondary cell contents

                 Do jj=j_start,lct_start(jc+1)-1

! get the particle index

                    j=at_list(jj)

! add an entry
                    ll=list(0,i)+1
                    list(ll,i)=j
                    list(0,i)=ll

! end of loop over the secondary cell contents

                 End Do

! overflow is to occur

              Else

! loop over the secondary cell contents

                 Do jj=j_start,lct_start(jc+1)-1

! get the particle index

                    j=at_list(jj)

! check for overflow and add an entry

                    ll=list(0,i)+1
                    If (ll <= mxlist) Then
                       list(ll,i)=j
                    Else
                       ibig=Max(ibig,ll)
                       safe=.false.
                    End If
                    list(0,i)=ll

! end of loop over the secondary cell contents

                 End Do

              End If

! no bypass on real space cutoff check for safe cells when
! atom pairs' distances are not guaranteed to be within the cutoff
! distances in real space are needed for checking the cutoff criterion

           Else

! loop over the secondary cell contents

              Do jj=j_start,lct_start(jc+1)-1

! get the particle index

                 j=at_list(jj)

! check cutoff criterion (all atom pairs MUST BE within the cutoff)
!                 rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2
! use reordered coordinate list in order to get "ordered" rather than "random" access

                 rsq=(xxt(jj)-xxx(i))**2+(yyt(jj)-yyy(i))**2+(zzt(jj)-zzz(i))**2
                 If (rsq <= rcsq) Then

! check for overflow and add an entry

                    ll=list(0,i)+1
                    If (ll <= mxlist) Then
                       list(ll,i)=j
                    Else
                       safe=.false.
                       ibig=Max(ibig,ll)
                    End If
                    list(0,i)=ll

! end of cutoff criterion check

                 End If

! end of loop over the secondary cell contents

              End Do

! end of bypass on real space cutoff check for safe cells

           End If

! secondary loop over neighbouring cells

        End Do

! end of loop over the primary cell contents

     End Do

! end of loops over domain cells

  End Do

! secondary loop over domain's border subcells only

  ipass=2
  Do icell=1,cell_bor(0)

! index of the primary cell and its coordinates

     ic=cell_bor(icell)
     ix=Mod(ic-1,nlx + 2*nlp)
     iz=(ic-1)/((nlx + 2*nlp)*(nly + 2*nlp))
     iy=(ic-1)/(nlx + 2*nlp) - (nly + 2*nlp)*iz

! loop over the primary cell contents

     Do ii=lct_start(ic),lct_start(ic+1)-1

        i=at_list(ii) ! get the particle index, by construction [1,natms]

! secondary loop over neighbouring cells,
! when ipass==2 exclude self-self (i.e. domain - kk=1)

        Do kk=ipass,nsbcll

! describe a negative non-repeatable semi-ball - ipass==2

           jx=ix-nix(kk)
           jy=iy-niy(kk)
           jz=iz-niz(kk)

! be on halo cells only - ipass==2

           If ( (jx <= nlx0e) .or. (jx >= nlx1s) .or. &
                (jy <= nly0e) .or. (jy >= nly1s) .or. &
                (jz <= nlz0e) .or. (jz >= nlz1s) ) Then

! index of the secondary cell

              jc=1+jx+(nlx+2*nlp)*(jy+(nly+2*nlp)*jz)

! get the secondary list's starting particle index as
! ipass==2 the secondary cell is always different from the
! primary cell get the head of chain of the secondary cell

              If (jc /= ic) j_start=lct_start(jc)

! bypass on real space cutoff check for safe cells when
! atom pairs' distances are guaranteed to be within the cutoff

              If (nir(kk)) Then

! check for overflow

                 ll=list(0,i)+lct_start(jc+1)-j_start
                 If (ll <= mxlist) Then

! loop over the secondary cell contents

                    Do jj=j_start,lct_start(jc+1)-1

! get the particle index

                       j=at_list(jj)

! add an entry
                       ll=list(0,i)+1
                       list(ll,i)=j
                       list(0,i)=ll

! end of loop over the secondary cell contents

                    End Do

! overflow is to occur

                 Else

! loop over the secondary cell contents

                    Do jj=j_start,lct_start(jc+1)-1

! get the particle index

                       j=at_list(jj)

! check for overflow and add an entry

                       ll=list(0,i)+1
                       If (ll <= mxlist) Then
                          list(ll,i)=j
                       Else
                          ibig=Max(ibig,ll)
                          safe=.false.
                       End If
                       list(0,i)=ll

! end of loop over the secondary cell contents

                    End Do

                 End If

! no bypass on real space cutoff check for safe cells when
! atom pairs' distances are not guaranteed to be within the cutoff
! distances in real space are needed for checking the cutoff criterion

              Else

! loop over the secondary cell contents

                 Do jj=j_start,lct_start(jc+1)-1

! get the particle index

                    j=at_list(jj)

! check cutoff criterion (all atom pairs MUST BE within the cutoff)
!                    rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2
! use reordered coordinate list in order to get "ordered" rather than "random" access

                    rsq=(xxt(jj)-xxx(i))**2+(yyt(jj)-yyy(i))**2+(zzt(jj)-zzz(i))**2
                    If (rsq <= rcsq) Then

! check for overflow and add an entry

                       ll=list(0,i)+1
                       If (ll <= mxlist) Then
                          list(ll,i)=j
                       Else
                          safe=.false.
                          ibig=Max(ibig,ll)
                       End If
                       list(0,i)=ll

! end of cutoff criterion check

                    End If

! end of loop over the secondary cell contents

                 End Do

! end of bypass on real space cutoff check for safe cells

              End If

! end of inner if-block on cells and borders

           End If

! secondary loop over neighbouring cells

        End Do

! end of loop over the primary cell contents

     End Do

! end of loops over icell

  End Do

! terminate job if neighbour list array exceeded

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Call gmax(ibig)
     Call warning(290,Real(ibig,wp),Real(mxlist,wp),0.0_wp)
     Call error(106)
  End If

! Rear down frozen pairs

  If (megfrz > 1) Then
     Do i=1,natms
        l_end=list(0,i)
        m_end=l_end

        ii=lfrzn(i)
        If (ii > 0) Then
           Do kk=l_end,1,-1
              j =list(kk,i)
              jj=lfrzn(j)
              If (jj > 0) Then
                 If (kk < m_end) Then
                    list(kk,i)=list(m_end,i)
                    list(m_end,i)=j
                 End If
                 m_end=m_end-1
              End If
           End Do
        End If

        list(-2,i)=list(0,i) ! End of full non-repeatable half VNL (FNRH VNL)
        list( 0,i)=m_end     ! End of new list with no frozen pairs (NFP)
     End Do
  Else
     Do i=1,natms
        list(-2,i)=list(0,i) ! End of full non-repeatable half VNL (FNRH VNL)
     End Do
  End If

! Rear down excluded pairs on top of the frozen ones

  If (lbook) Then
     Do i=1,natms
        l_end=list(0,i)
        m_end=l_end

        ii=lexatm(0,i)
        If (ii > 0) Then
           Do kk=l_end,1,-1
              j =list(kk,i)
              jj=ltg(j)
              If (match(jj,ii,lexatm(1:ii,i))) Then
                 If (kk < m_end) Then
                    list(kk,i)=list(m_end,i)
                    list(m_end,i)=j
                 End If
                 m_end=m_end-1
              End If
           End Do
        End If

        list(-1,i)=list(0,i) ! End of NFP FNRH VNL
        list( 0,i)=m_end     ! End of new list with no excluded interactions (NXI)
     End Do

! CHARMM core-shell screened electrostatic induction interactions
! Push up CHARMM pairs at the top of the bonded part of the list

     If (keyind == 1) Then
        Do i=1,natms
           l_end=list(-1,i) ! search marker to move up
           m_end=list( 0,i) ! CHARMM marker to move down

           ii=lchatm(0,i)
           If (ii > 0) Then ! find what the local sublist CHARMM marker is
outside:      Do While (l_end > m_end+1)    ! Only when space for swap exists

! Check for space at the top

                 j =list(m_end+1,i)
                 jj=ltg(j)
                 If (match(jj,ii,lchatm(1:ii,i))) Then
                    m_end=m_end+1    ! move down CHARMM marker
                    Cycle outside
                 End If

! if a swap can be made at the top space m_end+1
! check for a qualifier (l_end) at the bottom

inside:          Do While (l_end > m_end+1) ! Only when space for swap exists
                    j =list(l_end,i)
                    jj=ltg(j)
                    If (match(jj,ii,lchatm(1:ii,i))) Then
                       ibig            = list(m_end+1,i)
                       list(m_end+1,i) = list(l_end,i)
                       list(l_end,i)   = ibig

                       l_end=l_end-1 ! move up search marker
                       m_end=m_end+1 ! move down CHARMM marker
                       Cycle outside
                    End If
                    l_end=l_end-1    ! move up search marker
                    Cycle inside
                 End Do inside
              End Do outside
           End If
           list(-3,i)=m_end     ! CHARMM end within NFP FNRH VNL, offset wrt NXI end
        End Do
     Else
        Do i=1,natms
           list(-3,i)=list(0,i) ! CHARMM end coincides with NXI end
        End Do
     End If
  Else
     Do i=1,natms
        list(-1,i)=list(0,i)    ! End of NFP FNRH VNL
        list(-3,i)=list(0,i)    ! CHARMM end coincides with NXI end
     End Do
  End If

! check on minimum separation distance between VNL pairs at re/start

  If (l_dis) Then
!     l_dis=.false. ! at re/start ONLY

     cnt=0.0_wp
     Do i=1,natms
        ii=ltg(i)
        ll=legshl(0,i)

!        iz=(which_cell(i)-1)/((nlx + 2*nlp)*(nlx + 2*nlp))
!        iy=(which_cell(i)-1)/(nlx + 2*nlp) - (nly + 2*nlp)*iz
!        ix=Mod(which_cell(i)-1,nlx + 2*nlp)

        Do kk=1,list(-2,i)
           j =list(kk,i)
           jj=ltg(j)

!           jz=(which_cell(j)-1)/((nlx + 2*nlp)*(nlx + 2*nlp))
!           jy=(which_cell(j)-1)/(nlx + 2*nlp) - (nly + 2*nlp)*jz
!           jx=Mod(which_cell(j)-1,nlx + 2*nlp)

! Exclude core-shell units' pairs from the check

           If (ll /= 0) Then ! the primary particle is part of a unit
              If (j <= natms) Then ! can check directly if the pair is part of the same unit
                 If (legshl(0,j) /= 0) Then ! the secondary particle is part of a unit
                    If (legshl(1,i) == legshl(1,j)) Cycle ! both are part of the same unit
                 End If
              Else                 ! cannot check directly
                 If (ll > 0) Then ! the primary particle is a core
                    If (listshl(2,legshl(1,i)) == jj) Cycle ! jj is the shell of that core
                 Else               ! the primary particle is a shell
                    If (listshl(1,legshl(1,i)) == jj) Cycle ! jj is the core of that shell
                 End If
              End If
           End If

           If (j <= natms .or. ii < jj) Then
              cnt(1)=cnt(1)+1.0_wp ! sum up all pairs (rlnk=rcut+rpad)

              det=Sqrt((xxx(i)-xxx(j))**2+(yyy(i)-yyy(j))**2+(zzz(i)-zzz(j))**2)

              If (det < r_dis) Then
                 safe=.false.
                 Write(nrite,'(/,1x,a,2(i10,a),f5.3,a)')                      &
                      '*** warning - the pair with global indeces: '      ,   &
                      ii,'  &',jj,'  violates minimum separation distance (', &
                      det,' Angs) ***'
                 cnt(0)=cnt(0)+1.0_wp ! sum up violators
              End If

              If (kk <= list(0,i)) Then
                 If (det <  rcut) cnt(2)=cnt(2)+1.0_wp ! sum up all pairs (rcut, electrostatics)
                 If (det <  rvdw) cnt(3)=cnt(3)+1.0_wp ! sum up all pairs (rvdw, vdw)
                 If (det <= rmet) cnt(4)=cnt(4)+1.0_wp ! sum up all pairs (rmet, metal)
              End If
           End If
        End Do
     End Do

     If (mxnode > 1) Then
        Call gcheck(safe,"enforce")
        Call gsum(cnt)
     End If

     If (idnode == 0) Then
        If (.not.safe) Write(nrite,'(/,1x,a,i20,2a,f7.3,a,/)')                &
        '*** warning - ', Int(cnt(0),li), ' pair(s) of particles in CONFIG ', &
        'violate(s) the minimum separation distance of ',r_dis,' Angs ***'

        Write(nrite,'(/,1x,a)') &
        'Pair totals of short range interactions over cutoffs (in Angstroms):'
        If (Abs(rlnk-rcut) > 1.0e-6_wp) Write(nrite,'(6x,a,i20,a,f7.3)') &
        'extended       -  ', Int(cnt(1),li), '  within rlnk = ', rlnk
        Write(nrite,'(6x,a,i20,a,f7.3)') &
        'electrostatics -  ', Int(cnt(2),li), '  within rcut = ', rcut

        Write(nrite,'(6x,a,i20,a,f7.3)') &
        'van der Waals  -  ', Int(cnt(3),li), '  within rvdw = ', rvdw

        Write(nrite,'(6x,a,i20,a,f7.3,/)') &
        'metal          -  ', Int(cnt(4),li), '  within rmet = ', rmet
     End If
  End If

  Deallocate (nix,niy,niz,                   Stat=fail(1))
  Deallocate (which_cell,at_list,            Stat=fail(2))
  Deallocate (lct_count,lct_start,lct_where, Stat=fail(3))
  Deallocate (cell_dom,cell_bor,             Stat=fail(4))
  Deallocate (xxt,yyt,zzt,                   Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine link_cell_pairs
