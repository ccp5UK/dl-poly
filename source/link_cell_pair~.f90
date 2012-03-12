Subroutine link_cell_pairs(imcon,rcut,lbook,megfrz)

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
  Use comms_module,       Only : idnode,mxnode,gcheck,gmax,gsum
  Use setup_module
  Use domains_module,     Only : idx,idy,idz, nprx,npry,nprz, &
                                 nprx_r,npry_r,nprz_r
  Use config_module,      Only : cell,natms,nlast,ltg,lfrzn, &
                                 xxx,yyy,zzz,lexatm,list
  Use development_module, Only : l_dis,r_dis

  Implicit None

  Logical,            Intent( In    ) :: lbook
  Integer,            Intent( In    ) :: imcon,megfrz
  Real( Kind = wp ) , Intent( In    ) :: rcut

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: cut,rcsq

  Logical           :: safe,lx0,lx1,ly0,ly1,lz0,lz1,match

  Integer           :: fail(1:2),l_end,m_end, &
                       icell,ncells,ipass,    &
                       kk,ll, ibig,i,ii,j,jj, &
                       nlx,nly,nlz,           &
                       ix,iy,iz,ic,           &
                       jx,jy,jz,jc

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
  Allocate (link(1:mxatms),lct(0:mxcell),              Stat=fail(1))
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

  dispx=celprp(7)/(cut*nprx_r)
  dispy=celprp(8)/(cut*npry_r)
  dispz=celprp(9)/(cut*nprz_r)

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

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions (ALL - halo included) from centred
! Cartesian coordinates to reduced space coordinates

  Do i=1,nlast
     xxt(i)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     yyt(i)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     zzt(i)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
  End Do

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
! ensured that such is supplied).  Such large positive halo may
! lead to EDGE EFFECTs - link-cells constituting the positive
! halo may have larger dimensions than the domain link-cells.
!***************************************************************

! Form linked list
! Initialise link arrays

  Do i=1,nlast
     link(i)=0
  End Do
  Do icell=0,ncells
     lct(icell)=0
  End Do

  Do i=nlast,natms+1,-1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

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

! Correction for halo particles (natms+1,nlast) of this domain
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
  Do i=natms,1,-1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

! Get cell coordinates accordingly

     ix = Int(xdc*(xxt(i)+0.5_wp)) + jx
     iy = Int(ydc*(yyt(i)+0.5_wp)) + jy
     iz = Int(zdc*(zzt(i)+0.5_wp)) + jz

! Correction for domain (idnode) only particles (1,natms) but due to
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

! initialise verlet neighbourlist (VNL) arrays

!  list=0               ! (DEBUG)
  list(-2:0,1:natms)=0 ! (COUNTING DIMENSIONS ONLY)

! initial values of control variables

  ibig=0
  safe=.true.

! loop over the domain's cells only (ipass=1) and
! over the domain's border cells only (ipass=2)

  Do ipass=1,2

! primary loop over domain only cells

     Do iz=1,nlz
        Do iy=1,nly
           Do ix=1,nlx

! When ipass=2 be on the domain's border link-cells only

              If ( (ipass == 1) .or. &
                   (ix == 1) .or. (ix == nlx) .or. &
                   (iy == 1) .or. (iy == nly) .or. &
                   (iz == 1) .or. (iz == nlz) ) Then

! index of the primary cell

                 ic=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! head of chain of the primary cell
! by construction it may contain [0,natms]

                 i=lct(ic)

! bypass the primary cell if empty or emptied

                 If (i > 0) Then

! loop over the primary cell contents

100                 Continue

! secondary loop over neighbouring cells, when ipass=2
! exclude self-self (i.e. domain - kk=1)

                    Do kk=ipass,nsbcll

                       If (ipass == 1) Then ! non-negative non-repeatable semi-ball
                          jx=ix+nix(kk)
                          jy=iy+niy(kk)
                          jz=iz+niz(kk)
                       Else ! If (ipass = 2) Then ! negative non-repeatable semi-ball
                          jx=ix-nix(kk)
                          jy=iy-niy(kk)
                          jz=iz-niz(kk)
                       End If

! be on domain + possible positive halo cells only - ipass=1
! be on halo cells only - ipass=2

                       If ( (ipass == 1) .or. &
                            (jx == 0) .or. (jx == nlx+1) .or. &
                            (jy == 0) .or. (jy == nly+1) .or. &
                            (jz == 0) .or. (jz == nlz+1) ) Then

! index of the secondary cell

                          jc=1+jx+(nlx+2)*(jy+(nly+2)*jz)

! get the secodary cell particle
! by construction it may contain [0,natms], where 0=empty or emptied

! if the secondary cell is different from the primary cell

                          If (jc /= ic) Then

! get the head of chain of the secondary cell

                             j=lct(jc)

! if the secondary cell is same as the primary cell

                          Else ! only when ipass=1

! get the next in line from the primary cell running index

                             j=link(i)

                          End If

! bypass the secondary cell if empty or emptied

                          If (j > 0) Then

! loop over the secondary cell contents

200                          Continue

! check cutoff criterion (all atom pairs MUST BE within the cutoff)


                             rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2
                             If (rsq <= rcsq) Then

! check for overfloat and add an entry

                                ll=list(0,i)+1
                                If (ll <= mxlist) Then
                                   list(ll,i)=j
                                Else
                                   ibig=Max(ibig,ll)
                                   safe=.false.
                                End If
                                list(0,i)=ll

                             End If

! get a new secondary cell particle

                             j=link(j)
                             If (j > 0) Go To 200

! end of bypass of empty or emptied secondary cell

                          End If

! end of inner if-block on cells and borders

                       End If

! end of secondary loop over cells

                    End Do

! get a new primary cell particle

                    i=link(i) ! by construction [0,natms]
                    If (i > 0) Go To 100

! end of bypass of empty or emptied primary cell

                 End If

! end of outer if-block on domain or halo

              End If

! end of loops over ix,iy,iz

           End Do
        End Do
     End Do

! end of loop over passes

  End Do

! terminate job if neighbour list array exceeded

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Call gmax(ibig)
     Call warning(290,Real(ibig,wp),Real(mxlist,wp),0.0_wp)
     Call error(106)
  End If

! check on minimum separation distance between VNL pairs at re/start

  If (l_dis) Then
     l_dis=.false. ! at re/start ONLY
     r_dis=r_dis**2

     det=0.0_wp
     Do i=1,natms
        ii=ltg(i)
        Do kk=1,list(0,i)
           j =list(kk,i)
           jj=ltg(j)
           If (j <= natms .or. ii < jj) Then
              rsq=(xxx(i)-xxx(j))**2+(yyy(i)-yyy(j))**2+(zzz(i)-zzz(j))**2
              If (rsq < r_dis) Then
                 Write(nrite,'(/,1x,a,2i10,a)')                   &
                      '*** warning - pair with global indeces: ', &
                      ii,jj,' violates minimum separation distance'
                 safe=.false.
                 det=det+1.0_wp
              End If
           End If
        End Do
     End Do
     r_dis=Sqrt(r_dis)

     If (mxnode > 1) Then
         Call gcheck(safe)
         Call gsum(det)
     End If

     If ((.not.safe) .and. idnode == 0) Write(nrite,'(/,1x,a,i0,2a)') &
        '*** warning - ', Int(det,ip), ' number of pairs violated ',  &
        'the minimum separation distance ***'
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
  Else
     Do i=1,natms
        list(-1,i)=list(0,i) ! End of NFP FNRH VNL
     End Do
  End If

  Deallocate (link,lct,    Stat=fail(1))
  Deallocate (xxt,yyt,zzt, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine link_cell_pairs
