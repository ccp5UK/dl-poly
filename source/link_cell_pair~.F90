Subroutine link_cell_pairs(rcut,rlnk,rvdw,rmet,pdplnc,lbook,megfrz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for the verlet neighbour list based on link-cell
! method.
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,       Only : idnode,mxnode,gcheck,gmax,gsum
  Use setup_module
  Use domains_module,     Only : idx,idy,idz, nprx,npry,nprz, &
                                 r_nprx,r_npry,r_nprz
  Use config_module,      Only : cell,natms,nlast,ltg,lfrzn, &
                                 xxx,yyy,zzz,lexatm,list
  Use core_shell_module,  Only : listshl,legshl
  Use mpoles_module,      Only : keyind,lchatm
  Use development_module, Only : l_dis,r_dis

  Implicit None

  Logical,            Intent( In    ) :: lbook
  Integer,            Intent( In    ) :: megfrz
  Real( Kind = wp ) , Intent( In    ) :: rcut,rlnk,rvdw,rmet,pdplnc

  Logical           :: safe,lx0,lx1,ly0,ly1,lz0,lz1,match

  Integer           :: fail,l_end,m_end,      &
                       icell,ncells,ipass,    &
                       kk,ll, ibig,i,ii,j,jj, &
                       nlx,nly,nlz,           &
                       ix,iy,iz,ic,           &
                       jx,jy,jz,jc

  Real( Kind = wp ) :: cut,rcsq,rsq,det,rcell(1:9),celprp(1:10),cnt(0:4), &
                       x,y,z, x1,y1,z1, dispx,dispy,dispz, xdc,ydc,zdc

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

  fail=0
  Allocate (link(1:mxatms),lct(0:ncells), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs allocation failure, node: ', idnode
     Call error(0)
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

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

  Do i=nlast,natms+1,-1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

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

     If (ix >= 0 .and. iy >= 0 .and. iz >= 0) Then

! Correction for halo particles (natms+1,nlast) of this domain
! (idnode) but due to some tiny numerical inaccuracy kicked into
! the domain only link-cell space

        lx0=(ix >     0)
        lx1=(ix < nlx+1)
        ly0=(iy >     0)
        ly1=(iy < nly+1)
        lz0=(iz >     0)
        lz1=(iz < nlz+1)
        If ( (lx0 .and. lx1) .and. &
             (ly0 .and. ly1) .and. &
             (lz0 .and. lz1) ) Then

! Put the closest to the halo coordinate in the halo

           x1=Abs(x-0.5_wp*Sign(1.0_wp,x))
           y1=Abs(y-0.5_wp*Sign(1.0_wp,y))
           z1=Abs(z-0.5_wp*Sign(1.0_wp,z))
           If      (x1 <= y1 .and. x1 <= z1) Then
              If (x < 0.0_wp) Then
                 ix=0
              Else
                 ix=nlx+1
              End If
           Else If (y1 <= x1 .and. y1 <= z1) Then
              If (y < 0.0_wp) Then
                 iy=0
              Else
                 iy=nly+1
              End If
           Else
              If (z < 0.0_wp) Then
                 iz=0
              Else
                 iz=nlz+1
              End If
           End If
        End If

! Check for positively bound residual halo

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

     ix = Int(xdc*(x+0.5_wp)) + jx
     iy = Int(ydc*(y+0.5_wp)) + jy
     iz = Int(zdc*(z+0.5_wp)) + jz

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

! loop over the domain's cells only (ipass==1) and
! over the domain's border cells only (ipass==2)

  Do ipass=1,2

! primary loop over domain only cells

     Do iz=1,nlz
        Do iy=1,nly
           Do ix=1,nlx

! When ipass==2 be on the domain's border link-cells only

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

! secondary loop over neighbouring cells, when ipass==2
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

! be on domain + possible positive halo cells only - ipass==1
! be on halo cells only - ipass==2

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

                          Else ! only when ipass==1

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

! check for overflow and add an entry

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

  Deallocate (link,lct, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'link_cell_pairs deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine link_cell_pairs
