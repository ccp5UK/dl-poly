Subroutine ewald_spme_mforces(alpha,epsq,engcpe_rc,vircpe_rc,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! due to multipolar interactions in a periodic system using the smooth
! particle mesh ewald method for multipoles
!
! This version allows for extension to arbitrary order
!
! Note: (fourier) reciprocal space terms by global fft sumation;
!       (0) It is crucial to have large enough "limit", at least >=
!           kmaxa*kmaxb*kmaxc size needed for the Reshape operation!!!
!       (1) Dcft3, a FFT function, is needed.  It can be found in the
!           standard ESSL library (add -lessl in LDFLAGS in Makefile).
!       (2) ewald_spme_mforce~_d.f90 is not dependent on exchange_grid.f90
!           and parallel_fft.90, which is dependent on gpfaf90.
!           Therefore, all references in Makefile to these three files
!           are not needed.
!
! copyright - daresbury laboratory
! author    - h.a.boateng february 2016
! amended   - i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz,idx,idy,idz
  Use configuration,  Only : cell,volm,natms,nlast,xxx,yyy,zzz
  Use mpoles_module
  Use ewald_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: alpha,epsq
  Real( Kind = wp ), Intent(   Out ) :: engcpe_rc,vircpe_rc
  Real( Kind = wp ), Intent( InOut ) :: stress(1:9)

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ixb,iyb,izb, ixt,iyt,izt
  Real( Kind = wp ), Save :: ixbm1_r,iybm1_r,izbm1_r, &
                             ixtm0_r,iytm0_r,iztm0_r, &
                             kmaxa_r,kmaxb_r,kmaxc_r

  Integer              :: fail(1:4),limit, i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt, inc2,inc3

  Real( Kind = wp )    :: Dtpbsp,det,rcell(1:9),celprp(1:10),ralph,rvolm,scale,exclcoef, &
                          rcpcut,rcpct2,strs(1:9),eng,akv,tmp,bb1,bb2,bb3,               &
                          rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3,           &
                          rrkxy,rrkxz,rrkyx,rrkyz,rrkzx,rrkzy, sq1,sq2,sq3,              &
                          dtp,tq1,tq2,tq3,                                               &
                          imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10,            &
                          ka11,kb22,kc33,ka11sq,kb22sq,kc33sq,kakb,kakc,kbkc

  Real( Kind = wp )    :: imp(1:mximpl)
  Real( Kind = wp )    :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)

  Complex( Kind = wp ) :: vterm,pterm

! uni is the diagonal unit matrix

  Real( Kind = wp )    :: &
     uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

! B-spline coefficients

  Complex( Kind = wp ), Dimension( : ),       Allocatable, Save :: bscx,bscy,bscz
  Complex( Kind = wp ), Dimension( : ),       Allocatable       :: ww1,ww2,ww3

  Real( Kind = wp ),    Dimension( : ),       Allocatable       :: csp,buffer
  Real( Kind = wp ),    Dimension( : ),       Allocatable       :: txx,tyy,tzz
  Integer,              Dimension( : ),       Allocatable       :: ixx,iyy,izz,it
  Real( Kind = wp ),    Dimension( :,:,: ),   Allocatable       :: bsddx,bsddy,bsddz
  Real( Kind = wp ),    Dimension( :,: ),     Allocatable       :: bspx,bspy,bspz
  Real( Kind = wp ),    Dimension( : ),       Allocatable       :: bdx,bdy,bdz

! FFT workspace arrays

  Real( Kind = wp ),    Dimension( :,:,: ),   Allocatable, Save :: qqc
  Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Save :: qqq

  Real( Kind = wp ),    Dimension( :,:,:,: ), Allocatable, Save :: qtc
  Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Save :: qt1,qt2,qt3


  fail=0
  If (newjob) Then
     newjob = .false.

!!! BEGIN DD SPME VARIABLES
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

!!! END DD SPME VARIABLES

!!! BEGIN CARDINAL B-SPLINES SET-UP
! allocate the complex exponential arrays

     Allocate (ww1(1:kmaxa),ww2(1:kmaxb),ww3(1:kmaxc), Stat = fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ww arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! initialise the complex exponential arrays

     Call spl_cexp(kmaxa,kmaxb,kmaxc,ww1,ww2,ww3)

! allocate the global B-spline coefficients and the helper array

     Allocate (bscx(1:kmaxa),bscy(1:kmaxb),bscz(1:kmaxc), Stat = fail(1))
     Allocate (csp(1:mxspl),                              Stat = fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'bsc and cse arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! calculate the global B-spline coefficients

     Call bspcoe(mxspl,kmaxa,kmaxb,kmaxc,csp,bscx,bscy,bscz,ww1,ww2,ww3)

! deallocate the helper array and complex exponential arrays

     Deallocate (csp,         Stat = fail(1))
     Deallocate (ww1,ww2,ww3, Stat = fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'cse and ww arrays deallocation failure, node: ', idnode
        Call error(0)
     End If

!!! END CARDINAL B-SPLINES SET-UP

! allocate the global charge workspace arrays (NOT deallocated manually)

     Allocate (qqc(1:kmaxa,1:kmaxb,1:kmaxc),     Stat = fail(1))
     Allocate (qqq(1:kmaxa,1:kmaxb,1:kmaxc),     Stat = fail(2))
     Allocate (qtc(1:3,1:kmaxa,1:kmaxb,1:kmaxc), Stat = fail(3))
     Allocate (qt1(1:kmaxa,1:kmaxb,1:kmaxc), &
               qt2(1:kmaxa,1:kmaxb,1:kmaxc), &
               qt3(1:kmaxa,1:kmaxb,1:kmaxc),     Stat = fail(4) )
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'global charge workspace arrays allocation failure, node: ', idnode
        Call error(0)
     End If
  End If

  limit=Max(mxbuff,3*kmaxa*kmaxb*kmaxc) ! Try 3*kmaxa*kmaxb*kmaxc
  Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),buffer(1:limit),            Stat = fail(1))
  Allocate (ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms),it(1:mxatms),               Stat = fail(2))
  Allocate (bdx(0:mxspl),bsddx(0:mxspl,1:mxspl,1:mxatms), &
            bdy(0:mxspl),bsddy(0:mxspl,1:mxspl,1:mxatms), &
            bdz(0:mxspl),bsddz(0:mxspl,1:mxspl,1:mxatms),                         Stat = fail(3))
  Allocate (bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms), Stat = fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mforces allocation failure, node: ', idnode
     Call error(0)
  End If

! initialise coulombic potential energy and virial

  engsic    = 0.0_wp
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

     ixx(i)=Int(txx(i))
     iyy(i)=Int(tyy(i))
     izz(i)=Int(tzz(i))

! Detect if a particle is charged and in the MD cell or in its positive halo
! (t(i) >= -zero_plus) as the B-splines are negative directionally by propagation


     If (tzz(i) >= -zero_plus .and. &
         tyy(i) >= -zero_plus .and. &
         txx(i) >= -zero_plus .and. &
         Maxval(Abs(mplgfr(:,i))) > zero_plus)  Then
        it(i)=1
     Else
        it(i)=0
     End If
  End Do

! construct B-splines for atoms

  Call bspgen_mpoles(nlast,mxspl,txx,tyy,tzz,bspx,bspy,bspz,bsddx,bsddy,bsddz)

  Deallocate (txx,tyy,tzz,    Stat = fail(1))
  Deallocate (bspx,bspy,bspz, Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mforces allocation failure, node: ', idnode
     Call error(0)
  End If

! zero 3D charge array

  qtc = 0.0_wp; qqc = 0.0_wp

  ka11 = kmaxa_r*rcell(1)
  kb22 = kmaxb_r*rcell(5)
  kc33 = kmaxc_r*rcell(9)

  ka11sq = ka11*ka11; kb22sq = kb22*kb22; kc33sq = kc33*kc33
  kakb = ka11*kb22; kakc = ka11*kc33; kbkc = kb22*kc33

! construct 3D charge array

  Do i=1,nlast

! If a particle is charged and in the MD cell or in its positive halo
! (t(i) >= 0) as the B-splines are negative directionally by propagation

     If (it(i) == 1) Then

! get the multipoles for site i

        imp=mplgfr(:,i)

        If (mxompl > 0 .and. induce) Then

           imp(2)=imp(2)+indipx(i)
           imp(3)=imp(3)+indipy(i)
           imp(4)=imp(4)+indipz(i)

        End If

        llb = Max( izb, izz(i) - mxspl + 2 )
        llt = Min( izt, izz(i) + 1 )

        kkb = Max( iyb, iyy(i) - mxspl + 2 )
        kkt = Min( iyt, iyy(i) + 1 )

        jjb = Max( ixb, ixx(i) - mxspl + 2 )
        jjt = Min( ixt, ixx(i) + 1 )

        Do ll = llb, llt
           l = izz(i) - ll + 2

           bdz=bsddz(:,l,i)

           Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              bdy=bsddy(:,k,i)

              Do jj = jjb, jjt
                 j = ixx(i) - jj + 2

                 bdx=bsddx(:,j,i)

                 dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

                 dtp = imp1*bdx(0)*bdy(0)*bdz(0)

                 If (mxompl >= 1 ) Then
                    tmp = imp2*ka11*bdx(1)*bdy(0)*bdz(0); dtp=dtp+tmp; tq1=tq1+tmp
                    tmp = imp3*kb22*bdx(0)*bdy(1)*bdz(0); dtp=dtp+tmp; tq2=tq2+tmp
                    tmp = imp4*kc33*bdx(0)*bdy(0)*bdz(1); dtp=dtp+tmp; tq3=tq3+tmp
                 End If

                 If (mxompl == 2) Then
                    tmp = imp5*ka11sq  *bdx(2)*bdy(0)*bdz(0); dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                    tmp = imp8*kb22sq  *bdx(0)*bdy(2)*bdz(0); dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                    tmp = imp10*kc33sq *bdx(0)*bdy(0)*bdz(2); dtp=dtp+tmp; tq3=tq3+3.0_wp*tmp

                    tmp = imp6*kakb*bdx(1)*bdy(1)*bdz(0); dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                    tmp = imp7*kakc*bdx(1)*bdy(0)*bdz(1); dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                    tmp = imp9*kbkc*bdx(0)*bdy(1)*bdz(1); dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
                 End If

                 qqc(jj,kk,ll)   = qqc(jj,kk,ll)   + dtp
                 qtc(1,jj,kk,ll) = qtc(1,jj,kk,ll) + tq1
                 qtc(2,jj,kk,ll) = qtc(2,jj,kk,ll) + tq2
                 qtc(3,jj,kk,ll) = qtc(3,jj,kk,ll) + tq3
              End Do
           End Do
        End Do

     End If

  End Do

! Sum up qqc and qtc

  If (mxnode > 1) Then
     i=kmaxa*kmaxb*kmaxc
     buffer(1:i) = Reshape( qqc, (/i/) )

     Call gsum(buffer(1:i))

     qqc = Reshape( buffer(1:kmaxa*kmaxb*kmaxc), (/kmaxa,kmaxb,kmaxc/) )

     i=3*kmaxa*kmaxb*kmaxc
     buffer(1:i) = Reshape( qtc, (/i/) )

     Call gsum(buffer(1:i))

     qtc = Reshape( buffer(1:3*kmaxa*kmaxb*kmaxc), (/3,kmaxa,kmaxb,kmaxc/) )
  End If

! load charge array into complex array for FFT

  qqq=Cmplx(qqc , Kind = wp)

  qt1=Cmplx(qtc(1,:,:,:), Kind = wp)
  qt2=Cmplx(qtc(2,:,:,:), Kind = wp)
  qt3=Cmplx(qtc(3,:,:,:), Kind = wp)

! calculating inverse 3D FFT of generalized multipolar array (in place)

  inc2=kmaxa
  inc3=kmaxa*kmaxb

  Call Dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmaxa,kmaxb,kmaxc,-1,1.0_wp,buffer,limit)

  If (mxompl > 0) Then

     Call Dcft3(qt1,inc2,inc3,qt1,inc2,inc3,kmaxa,kmaxb,kmaxc,-1,1.0_wp,buffer,limit)
     Call Dcft3(qt2,inc2,inc3,qt2,inc2,inc3,kmaxa,kmaxb,kmaxc,-1,1.0_wp,buffer,limit)
     Call Dcft3(qt3,inc2,inc3,qt3,inc2,inc3,kmaxa,kmaxb,kmaxc,-1,1.0_wp,buffer,limit)

  End If

! set reciprocal space cutoff

  Call dcell(rcell,celprp)

  rcpcut=0.5_wp*Min(kmaxa_r*celprp(7),kmaxb_r*celprp(8),kmaxc_r*celprp(9))
  rcpcut=rcpcut*1.05_wp*twopi
  rcpct2=rcpcut**2

! initialise temporary stress tensor

  strs = 0.0_wp

! calculate convolution of charge array with gaussian function

  Do l=1,kmaxc

     ll=l-1
     If (l > kmaxc/2) ll=ll-kmaxc
     tmp=twopi*Real(ll,wp)

     rkx1=tmp*rcell(3)
     rky1=tmp*rcell(6)
     rkz1=tmp*rcell(9)

     bb3=Real( bscz(l)*Conjg(bscz(l)),wp )

     Do k=1,kmaxb

        kk=k-1
        If (k > kmaxb/2) kk=kk-kmaxb
        tmp=twopi*Real(kk,wp)

        rkx2=rkx1+tmp*rcell(2)
        rky2=rky1+tmp*rcell(5)
        rkz2=rkz1+tmp*rcell(8)

        bb2=bb3*Real( bscy(k)*Conjg(bscy(k)),wp )

        Do j=1,kmaxa

           jj=j-1
           If (j > kmaxa/2) jj=jj-kmaxa
           tmp=twopi*Real(jj,wp)

           rkx3=rkx2+tmp*rcell(1)
           rky3=rky2+tmp*rcell(4)
           rkz3=rkz2+tmp*rcell(7)

           bb1=bb2*Real( bscx(j)*Conjg(bscx(j)),wp )

           rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

!=================================================================
! For higher order contributions to the stress tensor

           rrkxy=0.0_wp ; rrkxz=0.0_wp
           rrkyx=0.0_wp ; rrkyz=0.0_wp
           rrkzx=0.0_wp ; rrkzy=0.0_wp

           If (rkx3 > 0.0_wp) Then

              rrkyx=rky3/rkx3; rrkzx=rkz3/rkx3

           End If

           If (rky3 > 0.0_wp) Then

              rrkxy=rkx3/rky3; rrkzy=rkz3/rky3

           End If

           If (rkz3 > 0.0_wp) Then

              rrkxz=rkx3/rkz3; rrkyz=rky3/rkz3

           End If

!====================================================================

           If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

! For  monopole contribution to the stress tensor

              vterm=bb1*Exp(ralph*rksq)/rksq*qqq(j,k,l)
              akv=2.0_wp*(1.0_wp/rksq-ralph)*Real( vterm*Conjg(qqq(j,k,l)),wp )

!=======================================================================
! For higher order contributions to the stress tensor

              pterm=2.0_wp*bb1*Exp(ralph*rksq)/rksq*Conjg(qqq(j,k,l))
              sq1=Real(pterm*qt1(j,k,l),wp)
              sq2=Real(pterm*qt2(j,k,l),wp)
              sq3=Real(pterm*qt3(j,k,l),wp)

!=======================================================================

              strs(1)=strs(1)-rkx3*rkx3*akv+sq1
              strs(5)=strs(5)-rky3*rky3*akv+sq2
              strs(9)=strs(9)-rkz3*rkz3*akv+sq3
              strs(2)=strs(2)-rkx3*rky3*akv+rrkxy*sq2
              strs(3)=strs(3)-rkx3*rkz3*akv+rrkxz*sq3
              strs(4)=strs(4)-rkx3*rky3*akv+rrkyx*sq1
              strs(6)=strs(6)-rky3*rkz3*akv+rrkyz*sq3
              strs(7)=strs(7)-rkx3*rkz3*akv+rrkzx*sq1
              strs(8)=strs(8)-rky3*rkz3*akv+rrkzy*sq2

              qqq(j,k,l)=vterm

           Else

              qqq(j,k,l)=(0.0_wp,0.0_wp)

           End If

        End Do
     End Do
  End Do

! scale strs and distribute per node

  strs = strs * scale / Real(mxnode,wp)

! calculate atomic energy

  Call Dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmaxa,kmaxb,kmaxc,1,1.0_wp,buffer,limit)

  Deallocate (buffer, Stat = fail(1))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mforces deallocation failure 0, node: ', idnode
     Call error(0)
  End If

  eng = Real(Sum(qqq*qqc),wp)

! scale eng and distribute per node

  eng = eng * scale / Real(mxnode,wp)

! Second part of the monopole contribution to the stress tensor
! calculate stress tensor (symmetrical, per node)

  strs      = strs      + eng*uni
  stress(1) = stress(1) + strs(1)
  stress(2) = stress(2) + 0.5_wp * (strs(2) + strs(4))
  stress(3) = stress(3) + 0.5_wp * (strs(3) + strs(7))
  stress(4) = stress(4) + 0.5_wp * (strs(2) + strs(4))
  stress(5) = stress(5) + strs(5)
  stress(6) = stress(6) + 0.5_wp * (strs(6) + strs(8))
  stress(7) = stress(7) + 0.5_wp * (strs(3) + strs(7))
  stress(8) = stress(8) + 0.5_wp * (strs(6) + strs(8))
  stress(9) = stress(9) + strs(9)

! distribute energy and virial terms (per node)

  engcpe_rc = eng
  vircpe_rc = -(strs(1)+strs(5)+strs(9))

! infrequent calculations copying

  If (l_cp) Then
     e_rc=engcpe_rc
     v_rc=vircpe_rc
     s_rc=strs
  End If

! calculate atomic forces

  Call spme_mforces(rcell,scale, ixx,iyy,izz,bsddx,bsddy,bsddz,qqq)

  Deallocate (ixx,iyy,izz,it,    Stat = fail(1))
  Deallocate (bdx,bdy,bdz,       Stat = fail(2))
  Deallocate (bsddx,bsddy,bsddz, Stat = fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mforces deallocation failure, node: ', idnode
     Call error(0)
  End If

Contains

  Subroutine spme_mforces(rcell,scale, ixx,iyy,izz,bsddx,bsddy,bsddz,qqq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic forces due to
! multipolar interactions in a periodic system using smooth particle
! mesh ewald method (fourier part)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov february 2016
! amended   - h.a.boateng may 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use configuration, Only : fxx,fyy,fzz

    Implicit None

    Integer,              Intent( In    ) :: ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms)
    Real( Kind = wp ),    Intent( In    ) :: scale,rcell(1:9),                &
                                             bsddx(0:mxspl,1:mxspl,1:mxatms), &
                                             bsddy(0:mxspl,1:mxspl,1:mxatms), &
                                             bsddz(0:mxspl,1:mxspl,1:mxatms)

    Complex( Kind = wp ), Intent( In    ) :: qqq(1:kmaxa,1:kmaxb,1:kmaxc)

    Integer           :: fail,i,j,k,l, jj,kk,ll
    Real( Kind = wp ) :: Dtpbsp,                                                       &
                         ka11,kb22,kc33,ka11sq,kb22sq,kc33sq,                          &
                         ka11cu,kb22cu,kc33cu,kakb,kakc,kbkc,                          &
                         kasqkb,kasqkc,kakbsq,kakcsq,kbsqkc,kbkcsq,kakbkc,             &
                         tmp,gmp,fff(0:3),fix,fiy,fiz,qsum,                            &
                         tix,tiy,tiz,dtp,tq1,tq2,tq3,tq4,tq5,tq6,tq7,tq8,tq9,tq10,     &
                         imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10,           &
                         impx1,impx2,impx3,impx4,impx5,impx6,impx7,impx8,impx9,impx10, &
                         impy1,impy2,impy3,impy4,impy5,impy6,impy7,impy8,impy9,impy10, &
                         impz1,impz2,impz3,impz4,impz5,impz6,impz7,impz8,impz9,impz10

    Real( Kind = wp ) :: imp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)

    Real( Kind = wp ), Dimension( : ), Allocatable :: bdx,bdy,bdz

    fail=0
    Allocate (bdx(0:mxspl),bdy(0:mxspl),bdz(0:mxspl), Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_mforces allocation failure, node: ', idnode
    End If

! Real values of kmax vectors

    ka11 = Real(kmaxa,wp)*rcell(1)
    kb22 = Real(kmaxb,wp)*rcell(5)
    kc33 = Real(kmaxc,wp)*rcell(9)

    ka11sq=ka11*ka11; ka11cu=ka11sq*ka11; kb22sq=kb22*kb22; kb22cu=kb22sq*kb22
    kc33sq=kc33*kc33; kc33cu=kc33sq*kc33; kakb=ka11*kb22; kakc=ka11*kc33
    kbkc=kb22*kc33; kasqkb=ka11sq*kb22; kasqkc=ka11sq*kc33; kakbsq=ka11*kb22sq
    kakcsq=ka11*kc33sq; kbsqkc=kb22sq*kc33; kbkcsq=kb22*kc33sq; kakbkc=ka11*kb22*kc33

    tmp=-2.0_wp*scale

    fff=0.0_wp

    Do i=1,natms

! get the multipoles for site i

       imp=mplgfr(:,i)

       If (mxompl > 0 .and. induce) Then

         imp(2)=imp(2)+indipx(i)
         imp(3)=imp(3)+indipy(i)
         imp(4)=imp(4)+indipz(i)

       End If

       If (Maxval(Abs(imp)) > zero_plus) Then

! scale imp multipoles

          imp=tmp*imp

! get the components for site i infinitesmial rotations

          impx=scale*mprotx(:,i)
          impy=scale*mproty(:,i)
          impz=scale*mprotz(:,i)

          imp1=imp(1); impx1=impx(1); impy1=impy(1); impz1=impz(1)

          If (mxompl >= 1) Then

             imp2=imp(2); imp3=imp(3); imp4=imp(4)

             impx2=impx(2); impx3=impx(3); impx4=impx(4)
             impy2=impy(2); impy3=impy(3); impy4=impy(4)
             impz2=impz(2); impz3=impz(3); impz4=impz(4)

          End If

          If (mxompl == 2) Then

             imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

             impx5=impx(5); impx6=impx(6); impx7=impx(7); impx8=impx(8); impx9=impx(9); impx10=impx(10)
             impy5=impy(5); impy6=impy(6); impy7=impy(7); impy8=impy(8); impy9=impy(9); impy10=impy(10)
             impz5=impz(5); impz6=impz(6); impz7=impz(7); impz8=impz(8); impz9=impz(9); impz10=impz(10)

           End If

! initialise forces & torques

          fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp
          tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp

          Do l=1,mxspl
             ll=izz(i)-l+2

             If (ll < 1) ll=ll+kmaxc

             bdz=bsddz(:,l,i)

             Do k=1,mxspl
                kk=iyy(i)-k+2

                If (kk < 1) kk=kk+kmaxb

                bdy=bsddy(:,k,i)

                Do j=1,mxspl
                   jj=ixx(i)-j+2

                   If (jj < 1) jj=jj+kmaxa

                   bdx=bsddx(:,j,i)

                   qsum=Real(qqq(jj,kk,ll),wp)

                   dtp = imp1*bdx(0)*bdy(0)*bdz(0)

                   tq2 = qsum*ka11*bdx(1)*bdy(0)*bdz(0); tq3 = qsum*kb22*bdx(0)*bdy(1)*bdz(0)
                   tq4 = qsum*kc33*bdx(0)*bdy(0)*bdz(1)

                   fix = fix+imp1*tq2
                   fiy = fiy+imp1*tq3
                   fiz = fiz+imp1*tq4

                   If (mxompl >= 1) Then

! torques

                      tix = tix + impx2*tq2 + impx3*tq3 + impx4*tq4
                      tiy = tiy + impy2*tq2 + impy3*tq3 + impx4*tq4
                      tiz = tiz + impz2*tq2 + impz3*tq3 + impz4*tq4

! forces

                      tq5 = qsum*ka11sq*bdx(2)*bdy(0)*bdz(0); tq6 = qsum*kakb*bdx(1)*bdy(1)*bdz(0)
                      tq7 = qsum*kakc*bdx(1)*bdy(0)*bdz(1);   tq8 = qsum*kb22sq*bdx(0)*bdy(2)*bdz(0)
                      tq9 = qsum*kbkc*bdx(0)*bdy(1)*bdz(1);   tq10= qsum*kc33sq*bdx(0)*bdy(0)*bdz(2)

                      fix = fix+imp2*tq5+imp3*tq6+imp4*tq7
                      fiy = fiy+imp2*tq6+imp3*tq8+imp4*tq9
                      fiz = fiz+imp2*tq7+imp3*tq9+imp4*tq10

                   End If

                   If (mxompl==2) Then

! torques

                      tix = tix + impx5*tq5+impx6*tq6+impx7*tq7+impx8*tq8+impx9*tq9+impx10*tq10
                      tiy = tiy + impy5*tq5+impy6*tq6+impy7*tq7+impy8*tq8+impy9*tq9+impy10*tq10
                      tiz = tiz + impz5*tq5+impz6*tq6+impz7*tq7+impz8*tq8+impz9*tq9+impz10*tq10

! forces
                      fix = fix+qsum*(imp5*ka11cu*bdx(3)*bdy(0)*bdz(0)+imp6*kasqkb*bdx(2)*bdy(1)*bdz(0)+ &
                                      imp7*kasqkc*bdx(2)*bdy(0)*bdz(1)+imp8*kakbsq*bdx(1)*bdy(2)*bdz(0)+ &
                                      imp9*kakbkc*bdx(1)*bdy(1)*bdz(1)+imp10*kakcsq*bdx(1)*bdy(0)*bdz(2))

                      fiy = fiy+qsum*(imp5*kasqkb*bdx(2)*bdy(1)*bdz(0)+imp6*kakbsq*bdx(1)*bdy(2)*bdz(0)+ &
                                      imp7*kakbkc*bdx(1)*bdy(1)*bdz(1)+imp8*kb22cu*bdx(0)*bdy(3)*bdz(0)+ &
                                      imp9*kbsqkc*bdx(0)*bdy(2)*bdz(1)+imp10*kbkcsq*bdx(0)*bdy(1)*bdz(2))

                      fiz = fiz+qsum*(imp5*kasqkc*bdx(2)*bdy(0)*bdz(1)+imp6*kakbkc*bdx(1)*bdy(1)*bdz(1)+ &
                                      imp7*kakcsq*bdx(1)*bdy(0)*bdz(2)+imp8*kbsqkc*bdx(0)*bdy(2)*bdz(1)+ &
                                      imp9*kbkcsq*bdx(0)*bdy(1)*bdz(2)+imp10*kc33cu*bdx(0)*bdy(0)*bdz(3))

                   End If
                End Do
             End Do
          End Do

! accumulate forces

          fff(0)=fff(0)+1.0_wp
          fff(1)=fff(1)+fix
          fff(2)=fff(2)+fiy
          fff(3)=fff(3)+fiz

! load forces

          fxx(i)=fxx(i)+fix
          fyy(i)=fyy(i)+fiy
          fzz(i)=fzz(i)+fiz

! and torque

          mptrqx(i)=mptrqx(i)+0.5_wp*tix
          mptrqy(i)=mptrqy(i)+0.5_wp*tiy
          mptrqz(i)=mptrqz(i)+0.5_wp*tiz

! infrequent calculations copying

          If (l_cp) Then
             fcx(i)=fcx(i)+fix
             fcy(i)=fcy(i)+fiy
             fcz(i)=fcz(i)+fiz
          End If

       End If

    End Do

! remove COM drift arising from SPME approximations

    If (mxnode > 1) Call gsum(fff)
    If (fff(0) > zero_plus) Then
       fff(1:3)=fff(1:3)/fff(0)

       Do i=1,natms
          imp=mplgfr(:,i)
          If (Maxval(Abs(imp)) > zero_plus) Then

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

    Deallocate (bdx,bdy,bdz, Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_mforces dealocation failure, node: ', idnode
    End If

  End Subroutine spme_mforces

End Subroutine ewald_spme_mforces
