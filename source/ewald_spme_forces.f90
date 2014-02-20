Subroutine ewald_spme_forces(alpha,epsq,engcpe_rc,vircpe_rc,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using the smooth particle mesh ewald method
! by Essmann et al. J. Chem. Phys. 103 (1995) 8577
!
! Note: (fourier) reciprocal space terms
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith & i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum,dlp_comm_world
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz,idx,idy,idz
  Use config_module,  Only : cell,volm,natms,nlast,chge,xxx,yyy,zzz
  Use ewald_module
  Use parallel_fft

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: alpha,epsq
  Real( Kind = wp ), Intent(   Out ) :: engcpe_rc,vircpe_rc
  Real( Kind = wp ), Intent( InOut ) :: stress(1:9)

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ixb,iyb,izb, ixt,iyt,izt
  Real( Kind = wp ), Save :: twopi,ixbm1_r,iybm1_r,izbm1_r, &
                                   ixtm0_r,iytm0_r,iztm0_r, &
                                   kmaxa_r,kmaxb_r,kmaxc_r,engsic

  Integer              :: fail(1:4), i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt

  Real( Kind = wp )    :: det,rcell(1:9),celprp(1:10),ralph,rvolm,scale,   &
                          rcpcut,rcpct2,strs(1:9),eng,akv,tmp,bb1,bb2,bb3, &
                          rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3

  Complex( Kind = wp ) :: vterm

! uni is the diagonal unit matrix

  Real( Kind = wp )    :: &
     uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

! B-spline coefficients

  Complex( Kind = wp ), Dimension( : ),   Allocatable, Save :: bscx,bscy,bscz
  Complex( Kind = wp ), Dimension( : ),   Allocatable, Save :: ww1,ww2,ww3

  Real( Kind = wp ),    Dimension( : ),   Allocatable       :: csp
  Real( Kind = wp ),    Dimension( : ),   Allocatable       :: txx,tyy,tzz
  Integer,              Dimension( : ),   Allocatable       :: ixx,iyy,izz,it
  Real( Kind = wp ),    Dimension( :,: ), Allocatable       :: bsdx,bsdy,bsdz
  Real( Kind = wp ),    Dimension( :,: ), Allocatable       :: bspx,bspy,bspz

! ijb context for parallel fft

  Integer, Save :: context

! ijb blocking factors for fft

  Integer, Save :: block_x
  Integer, Save :: block_y
  Integer, Save :: block_z

! ijb indexing arrays for x, y and z as used in parallel fft

  Integer, Dimension( : ), Allocatable, Save :: index_x
  Integer, Dimension( : ), Allocatable, Save :: index_y
  Integer, Dimension( : ), Allocatable, Save :: index_z

! temporary qqc

  Real( Kind = wp ) :: qqc_tmp

! ijb temporary workspace for parallel fft

  Complex( Kind = wp ), Dimension( :, :, : ), Allocatable, Save :: qqq_local
  Complex( Kind = wp ), Dimension( :, :, : ), Allocatable, Save :: pfft_work

  Real( Kind = wp ),    Dimension( :, :, : ), Allocatable, Save :: qqc_local

  Integer :: j_local, k_local, l_local

  fail=0
  Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),                            Stat = fail(1))
  Allocate (ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms),it(1:mxatms),               Stat = fail(2))
  Allocate (bsdx(1:mxspl,1:mxatms),bsdy(1:mxspl,1:mxatms),bsdz(1:mxspl,1:mxatms), Stat = fail(3))
  Allocate (bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms), Stat = fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_forces allocation failure, node: ', idnode
     Call error(0)
  End If

  If (newjob) Then
     newjob = .false.

! derivative of pi

     twopi=2.0_wp*pi

! 3D charge array construction (bottom and top) indices

     ixb=idx*(kmaxa/nprx)+1
     ixt=(idx+1)*(kmaxa/nprx)
     iyb=idy*(kmaxb/npry)+1
     iyt=(idy+1)*(kmaxb/npry)
     izb=idz*(kmaxc/nprz)+1
     izt=(idz+1)*(kmaxc/nprz)

     ixbm1_r=Real(ixb-1,wp)
     ixtm0_r=Nearest( Real(ixt,wp) , -1.0_wp )
     iybm1_r=Real(iyb-1,wp)
     iytm0_r=Nearest( Real(iyt,wp) , -1.0_wp )
     izbm1_r=Real(izb-1,wp)
     iztm0_r=Nearest( Real(izt,wp) , -1.0_wp )

! Real values of kmax vectors

     kmaxa_r=Real(kmaxa,wp)
     kmaxb_r=Real(kmaxb,wp)
     kmaxc_r=Real(kmaxc,wp)

! calculate self-interaction correction (per node)

     engsic=0.0_wp
     Do i=1,natms
        engsic=engsic+chge(i)**2
     End Do

     If (mxnode > 1) Call gsum(engsic)

     engsic=-r4pie0/epsq * alpha*engsic/sqrpi / Real(mxnode,wp)

! allocate the complex exponential arrays (NOT deallocated manually)

     Allocate (ww1(1:kmaxa),ww2(1:kmaxb),ww3(1:kmaxc), Stat = fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ww arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! initialise the complex exponential arrays

     Call spl_cexp(kmaxa,kmaxb,kmaxc,ww1,ww2,ww3)

! ijb set up the parallel fft and useful related quantities

     block_x = kmaxa / nprx
     block_y = kmaxb / npry
     block_z = kmaxc / nprz

     Call initialize_fft( 3, (/ kmaxa, kmaxb, kmaxc /), &
         (/ nprx, npry, nprz /), (/ idx, idy, idz /),   &
         (/ block_x, block_y, block_z /),               &
         dlp_comm_world, context )

! ijb set up the indexing arrays for each dimension (NOT deallocated manually)

     Allocate ( index_x( 1:block_x ), Stat = fail(1) )
     Allocate ( index_y( 1:block_y ), Stat = fail(2) )
     Allocate ( index_z( 1:block_z ), Stat = fail(3) )
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'SPME index arrays allocation failure, node: ', idnode
        Call error(0)
     End If

     Call pfft_indices( kmaxa, block_x, idx, nprx, index_x )
     Call pfft_indices( kmaxb, block_y, idy, npry, index_y )
     Call pfft_indices( kmaxc, block_z, idz, nprz, index_z )

! ijb workspace arrays for DaFT

     Allocate ( qqq_local( 1:block_x, 1:block_y, 1:block_z ), Stat = fail(1) )
     Allocate ( qqc_local( 1:block_x, 1:block_y, 1:block_z ), Stat = fail(2) )
     Allocate ( pfft_work( 1:block_x, 1:block_y, 1:block_z ), Stat = fail(3) )

     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'SPME DaFT workspace arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! calculate B-spline coefficients

     Allocate (bscx(1:kmaxa),bscy(1:kmaxb),bscz(1:kmaxc), Stat = fail(1))
     Allocate (csp(1:mxspl),                              Stat = fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'bsc and cse arrays allocation failure, node: ', idnode
        Call error(0)
     End If

     Call bspcoe(mxspl,kmaxa,kmaxb,kmaxc,csp,bscx,bscy,bscz,ww1,ww2,ww3)

     Deallocate (csp, Stat = fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'cse array deallocation failure, node: ', idnode
        Call error(0)
     End If
  End If

! initialise coulombic potential energy and virial

  engcpe_rc = 0.0_wp
  vircpe_rc = 0.0_wp

! set working parameters

  rvolm=twopi/volm
  ralph=-0.25_wp/alpha**2

! set scaling constant

  scale=rvolm*r4pie0/epsq

! Convert cell coordinates to fractional coordinates intervalled [0,1)
! (bottom left corner of MD cell) and stretch over kmaxs in different
! directions.  Only the halo (natms,nlast] has fractional coordinates
! outside the [0,1) interval.  In the worst case scenario of one
! "effective" link-cell per domain and one domain in the MD cell only,
! the halo will have fractional coordinates intervalled as
! [n,0)u[1,2), where -1 <= n < 0.  Only the positive halo is needed by
! the B-splines since they distribute/spread charge density in
! negative direction with length the length of the spline.
!
! The story has become more complicated with cutoff padding and the
! conditional updates of the VNL and thus the halo as now a domain
! (1:natms) particle can enter the halo and vice versa.  So DD
! bounding is unsafe!!!

  Call invert(cell,rcell,det)
  If (Abs(det) < 1.0e-6_wp) Call error(120)

  Do i=1,nlast
     txx(i)=kmaxa_r*(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)+0.5_wp)
     tyy(i)=kmaxb_r*(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)+0.5_wp)
     tzz(i)=kmaxc_r*(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)+0.5_wp)

! THIS IS NOW UNSAFE
! Get in DD bounds in kmax grid space in case tiny inaccuracies created edge effects
!
!     If (i <= natms) Then
!        If      (txx(i) < ixbm1_r) Then
!           txx(i)=ixbm1_r
!        Else If (txx(i) > ixtm0_r) Then
!           txx(i)=ixtm0_r
!        End If
!
!        If      (tyy(i) < iybm1_r) Then
!           tyy(i)=iybm1_r
!        Else If (tyy(i) > iytm0_r) Then
!           tyy(i)=iytm0_r
!        End If
!
!        If      (tzz(i) < izbm1_r) Then
!           tzz(i)=izbm1_r
!        Else If (tzz(i) > iztm0_r) Then
!           tzz(i)=iztm0_r
!        End If
!     End If

     ixx(i)=Int(txx(i))
     iyy(i)=Int(tyy(i))
     izz(i)=Int(tzz(i))

! Detect if a particle is charged and in the MD cell or in its positive halo
! (t(i) >= -zero_plus) as the B-splines are negative directionally by propagation

     If (tzz(i) >= -zero_plus .and. &
         tyy(i) >= -zero_plus .and. &
         txx(i) >= -zero_plus .and. &
         Abs(chge(i)) > zero_plus) Then
        it(i)=1
     Else
        it(i)=0
     End If
  End Do

! construct B-splines for atoms

  Call bspgen(nlast,mxspl,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz)

  Deallocate (txx,tyy,tzz, Stat = fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_forces allocation failure, node: ', idnode
     Call error(0)
  End If

! zero 3D charge array
! ijb DaFT version - only need set local bit to zero

  qqc_local = 0.0_wp

! construct 3D charge array
! ijb DaFT version - use array that holds only the local data

  Do i=1,nlast

! If a particle is charged and in the MD cell or in its positive halo
! (t(i) >= 0) as the B-splines are negative directionally by propagation

     If (it(i) == 1) Then
        bb3=chge(i)

!        Do l=1,mxspl
!           ll=izz(i)-l+2
!
!! If a particle's B-spline is entering this domain (originating from its
!! positive halo), i.e. <= i.t, and not just to start exiting it, i.e. >= i.b
!! In the limit of one domain in the MD cell (npr.=1, id.=0) i.t=kmax. and i.b=1
!
!           If (ll >= izb .and. ll <= izt) Then
!              l_local = ll - izb + 1
!              bb2=bb3*bspz(l,i)
!
!              Do k=1,mxspl
!                 kk=iyy(i)-k+2
!
!                 If (kk >= iyb .and. kk <= iyt) Then
!                    k_local = kk - iyb + 1
!                    bb1=bb2*bspy(k,i)
!
!                    Do j=1,mxspl
!                       jj=ixx(i)-j+2
!
!                       If (jj >= ixb .and. jj <= ixt) Then
!                          j_local = jj - ixb + 1
!                          det=bb1*bspx(j,i)
!
!                          qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
!                       End If
!                    End Do
!                 End If
!              End Do
!           End If
!        End Do

        llb = Max( izb, izz(i) - mxspl + 2 )
        llt = Min( izt, izz(i) + 1 )

        kkb = Max( iyb, iyy(i) - mxspl + 2 )
        kkt = Min( iyt, iyy(i) + 1 )

        jjb = Max( ixb, ixx(i) - mxspl + 2 )
        jjt = Min( ixt, ixx(i) + 1 )

        Select Case( jjt - jjb + 1 )

        Case Default

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 Do jj = jjb, jjt
                    j = ixx(i) - jj + 2

                    j_local = jj - ixb + 1
                    det=bb1*bspx(j,i)

                    qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
                 End Do
              End Do
           End Do

        Case( 0 )

        Case( 1 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 2 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !2
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 3 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !2
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !3
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 4 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !2
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !3
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !4
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 5 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !2
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !3
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !4
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !5
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 6 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !2
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !3
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !4
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !5
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !6
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 7 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 !1
                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !2
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !3
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !4
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !5
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !6
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 !7
                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        Case( 8 )

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1
              bb2=bb3*bspz(l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1
                 bb1=bb2*bspy(k,i)

                 jj = jjb

                 j = ixx(i) - jj + 2
                 j_local = jj - ixb + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det

                 j = j - 1
                 j_local = j_local + 1
                 det=bb1*bspx(j,i)

                 qqc_local(j_local,k_local,l_local)=qqc_local(j_local,k_local,l_local)+det
              End Do
           End Do

        End Select
     End If

  End Do

! load charge array into complex array for FFT

  qqq_local=Cmplx(qqc_local , Kind = wp)

! calculate inverse 3D FFT of charge array (in place)

  Call pfft(qqq_local,pfft_work,context,1)

! set reciprocal space cutoff

  Call dcell(rcell,celprp)

  rcpcut=0.5_wp*Min(kmaxa_r*celprp(7),kmaxb_r*celprp(8),kmaxc_r*celprp(9))
  rcpcut=rcpcut*1.05_wp*twopi
  rcpct2=rcpcut**2

! initialise temporary stress tensor

  strs = 0.0_wp

! calculate convolution of charge array with gaussian function
! ijb DaFT Version - only loop over the local stuff

  Do l_local=1,block_z
     l=index_z(l_local)

     ll=l-1
     If (l > kmaxc/2) ll=l-kmaxc-1
     tmp=twopi*Real(ll,wp)

     rkx1=tmp*rcell(3)
     rky1=tmp*rcell(6)
     rkz1=tmp*rcell(9)

     bb3=Real( bscz(l)*Conjg(bscz(l)),wp )

     Do k_local=1,block_y
        k=index_y(k_local)

        kk=k-1
        If (k > kmaxb/2) kk=k-kmaxb-1
        tmp=twopi*Real(kk,wp)

        rkx2=rkx1+tmp*rcell(2)
        rky2=rky1+tmp*rcell(5)
        rkz2=rkz1+tmp*rcell(8)

        bb2=bb3*Real( bscy(k)*Conjg(bscy(k)),wp )

        Do j_local=1,block_x
           j=index_x(j_local)

           jj=j-1
           If (j > kmaxa/2) jj=j-kmaxa-1
           tmp=twopi*Real(jj,wp)

           rkx3=rkx2+tmp*rcell(1)
           rky3=rky2+tmp*rcell(4)
           rkz3=rkz2+tmp*rcell(7)

           bb1=bb2*Real( bscx(j)*Conjg(bscx(j)),wp )

           rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

           If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

              vterm=bb1*Exp(ralph*rksq)/rksq*qqq_local(j_local,k_local,l_local)
              akv=2.0_wp*(1.0_wp/rksq-ralph)*Real( vterm*Conjg(qqq_local(j_local,k_local,l_local)),wp )
              strs(1)=strs(1)-rkx3*rkx3*akv
              strs(5)=strs(5)-rky3*rky3*akv
              strs(9)=strs(9)-rkz3*rkz3*akv
              strs(2)=strs(2)-rkx3*rky3*akv
              strs(3)=strs(3)-rkx3*rkz3*akv
              strs(6)=strs(6)-rky3*rkz3*akv
              qqq_local(j_local,k_local,l_local)=vterm

           Else

              qqq_local(j_local,k_local,l_local)=(0.0_wp,0.0_wp)

           End If
        End Do
     End Do
  End Do

! complete strs

  strs(4) = strs(2)
  strs(7) = strs(3)
  strs(8) = strs(6)

! ijb as only looped over local stuff, we need to gsum strs

  If (mxnode > 1) Call gsum(strs)

! scale strs and distribute per node

  strs = strs * scale / Real(mxnode,wp)

  Call pfft(qqq_local,pfft_work,context,-1)

! calculate atomic energy

  eng = 0.0_wp
  Do l=1,block_z
     Do k=1,block_y
        Do j=1,block_x
           qqc_tmp=Real(qqq_local(j,k,l),wp)
           eng=eng+qqc_local(j,k,l)*qqc_tmp
           qqc_local(j,k,l)=qqc_tmp
        End Do
     End Do
  End Do

! ijb as only looped over local stuff, we need to gsum the eng

  If (mxnode > 1) Call gsum(eng)

! scale eng and distribute per node

  eng = eng * scale / Real(mxnode,wp)

! calculate atomic forces

  Call spme_forces(rcell,scale, ixx,iyy,izz, bspx,bspy,bspz, bsdx,bsdy,bsdz, qqc_local, ixb,ixt, iyb,iyt, izb,izt)

! calculate stress tensor (symmetrical, per node)

  strs   = strs + eng*uni
  stress = stress + strs

! distribute energy and virial terms (per node)

  engcpe_rc = eng + engsic
  vircpe_rc = -(strs(1)+strs(5)+strs(9))

! infrequent calculations copying

  If (l_cp) Then
     e_rc=engcpe_rc
     v_rc=vircpe_rc
     s_rc=strs
  End If

  Deallocate (ixx,iyy,izz,it, Stat = fail(1))
  Deallocate (bsdx,bsdy,bsdz, Stat = fail(2))
  Deallocate (bspx,bspy,bspz, Stat = fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

Contains

  Subroutine spme_forces(rcell,scale, ixx,iyy,izz, bspx,bspy,bspz, bsdx,bsdy,bsdz, qqc_local, ixb,ixt, iyb,iyt, izb,izt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic forces in a periodic
! system using smooth particle mesh ewald method (fourier part)
!
! Note: qqc_local is shifted from its definition from above
!       and therefore there is no need for periodic images (!!)
!
! copyright - daresbury laboratory
! author    - w.smith october 1998
! amended   - i.t.todorov october 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Use kinds_f90
    Use comms_module,  Only : idnode,mxnode,gsum
!    Use setup_module
    Use config_module, Only : natms,chge,fxx,fyy,fzz
!    Use ewald_module

    Implicit None

    Integer,           Intent( In    ) :: ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms), &
                                          ixb,ixt, iyb,iyt, izb,izt
    Real( Kind = wp ), Intent( In    ) :: scale,rcell(1:9),                   &
        bsdx(1:mxspl,1:mxatms),bsdy(1:mxspl,1:mxatms),bsdz(1:mxspl,1:mxatms), &
        bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms), &
        qqc_local( ixb:ixt, iyb:iyt, izb:izt )

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: kmaxa_r,kmaxb_r,kmaxc_r

    Integer           :: fail, i,j,k,l, jj,kk,ll, ixdb,iydb,izdb
    Real( Kind = wp ) :: tmp,facx,facy,facz,fff(0:3),fx,fy,fz,fix,fiy,fiz,qsum, &
                         bdxl,bdyl,bdzl,bdxk,bdyk,bdzk,bdxj,bdyj,bdzj

    Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: qqc_domain

! Real values of kmax vectors

    If (newjob) Then
       newjob = .false.

       kmaxa_r=Real(kmaxa,wp)
       kmaxb_r=Real(kmaxb,wp)
       kmaxc_r=Real(kmaxc,wp)
    End If

    ixdb = ixb - mxspl
    iydb = iyb - mxspl
    izdb = izb - mxspl

    fail=0
    Allocate (qqc_domain( ixdb:ixt, iydb:iyt, izdb:izt ), Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_for_domain allocation failure, node: ', idnode
    End If

    Call exchange_grid( ixb, ixt, iyb, iyt, izb, izt,           &
                      ixdb, iydb, izdb, qqc_local, qqc_domain )

    tmp=-2.0_wp*scale
    facx=tmp*kmaxa_r
    facy=tmp*kmaxb_r
    facz=tmp*kmaxc_r

    fff=0.0_wp
    Do i=1,natms
       tmp=chge(i)

       If (Abs(tmp) > zero_plus) Then

! initialise forces

          fix=0.0_wp
          fiy=0.0_wp
          fiz=0.0_wp

          Do l=1,mxspl
             ll=izz(i)-l+2

             bdxl=tmp*facx*bspz(l,i)
             bdyl=tmp*facy*bspz(l,i)
             bdzl=tmp*facz*bsdz(l,i)

             Do k=1,mxspl
                kk=iyy(i)-k+2

                bdxk=bdxl*bspy(k,i)
                bdyk=bdyl*bsdy(k,i)
                bdzk=bdzl*bspy(k,i)

                Do j=1,mxspl
                   jj=ixx(i)-j+2

                   qsum=qqc_domain(jj,kk,ll)
                   bdxj=qsum*bdxk*bsdx(j,i)
                   bdyj=qsum*bdyk*bspx(j,i)
                   bdzj=qsum*bdzk*bspx(j,i)

                   fix=fix+bdxj
                   fiy=fiy+bdyj
                   fiz=fiz+bdzj
                End Do
             End Do
          End Do

          fx=fix*rcell(1)+fiy*rcell(2)+fiz*rcell(3)
          fy=fix*rcell(4)+fiy*rcell(5)+fiz*rcell(6)
          fz=fix*rcell(7)+fiy*rcell(8)+fiz*rcell(9)

! accumulate forces

          fff(0)=fff(0)+1.0_wp
          fff(1)=fff(1)+fx
          fff(2)=fff(2)+fy
          fff(3)=fff(3)+fz

! load forces

          fxx(i)=fxx(i)+fx
          fyy(i)=fyy(i)+fy
          fzz(i)=fzz(i)+fz

! infrequent calculations copying

          If (l_cp) Then
             fcx(i)=fcx(i)+fx
             fcy(i)=fcy(i)+fy
             fcz(i)=fcz(i)+fz
          End If

       End If
    End Do

! remove COM drift arising from SPME approximations

    If (mxnode > 1) Call gsum(fff)
    If (fff(0) > zero_plus) Then
       fff(1:3)=fff(1:3)/fff(0)

       Do i=1,natms
          If (Abs(chge(i)) > zero_plus) Then

             fxx(i)=fxx(i)-fff(1)
             fyy(i)=fyy(i)-fff(2)
             fzz(i)=fzz(i)-fff(3)

! infrequent calculations copying

             If (l_cp) Then
                fcx(i)=fcx(i)-fff(1)
                fcy(i)=fcy(i)-fff(2)
                fcz(i)=fcz(i)-fff(3)
             End If
          End If
       End Do
    End If

    Deallocate (qqc_domain, Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_for_domain dealocation failure, node: ', idnode
    End If

  End Subroutine spme_forces

End Subroutine ewald_spme_forces

Subroutine adjust_kmax( kmax, P )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to adjust a k-vector length
! with what DaFT can handle
!
! copyright - daresbury laboratory
! author    - i.j.bush august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use parallel_fft, Only : pfft_length_ok

  Implicit None

  Integer, Intent( InOut ) :: kmax
  Integer, Intent( In    ) :: P

  ! First make sure kmax is a multiple of P, and is at least as big as the input value
  If ( Mod( kmax, P ) /= 0 ) kmax = ( kmax / P + 1 ) * P

  ! Now check it has suitable factors
  Do While ( .not. pfft_length_ok( kmax / P ) )
     kmax = kmax + P
  End Do

End Subroutine adjust_kmax
