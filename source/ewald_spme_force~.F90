Subroutine ewald_spme_forces(alpha,epsq,engcpe_rc,vircpe_rc,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using the smooth particle mesh ewald method
! by Essmann et al. J. Chem. Phys. 103 (1995) 8577
!
! Note: (fourier) reciprocal space terms by global fft sumation;
!       (0) It is crucial to have large enough "limit", at least >=
!           kmaxa*kmaxb*kmaxc size needed for the Reshape operation!!!
!       (1) Dcft3, a FFT function, is needed.  It can be found in the
!           standard ESSL library (add -lessl in LDFLAGS in Makefile).
!       (2) ewald_spme_force~.f90 is not dependent on exchange_grid.f90
!           and parallel_fft.90, which is dependent on gpfa.f90.
!           Therefore, all references in Makefile to these three files
!           are not needed.
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith & i.j.bush february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup
  Use domains, Only : nprx,npry,nprz,idx,idy,idz
  Use configuration,  Only : cell,volm,natms,nlast,chge,xxx,yyy,zzz
  Use ewald_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: alpha,epsq
  Real( Kind = wp ), Intent(   Out ) :: engcpe_rc,vircpe_rc
  Real( Kind = wp ), Intent( InOut ) :: stress(1:9)

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ixb,iyb,izb, ixt,iyt,izt
  Real( Kind = wp ), Save :: kmaxa_r,kmaxb_r,kmaxc_r

  Integer              :: fail(1:4),limit, i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt, inc2,inc3

  Real( Kind = wp )    :: det,rcell(1:9),celprp(1:10),ralph,rvolm,scale,   &
                          rcpcut,rcpct2,strs(1:9),eng,akv,tmp,bb1,bb2,bb3, &
                          rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3

  Complex( Kind = wp ) :: vterm

! uni is the diagonal unit matrix

  Real( Kind = wp )    :: &
     uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

! B-spline coefficients

  Complex( Kind = wp ), Dimension( : ),     Allocatable, Save :: bscx,bscy,bscz
  Complex( Kind = wp ), Dimension( : ),     Allocatable       :: ww1,ww2,ww3

  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: csp,buffer
  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: txx,tyy,tzz
  Integer,              Dimension( : ),     Allocatable       :: ixx,iyy,izz,it
  Real( Kind = wp ),    Dimension( :,: ),   Allocatable       :: bsdx,bsdy,bsdz
  Real( Kind = wp ),    Dimension( :,: ),   Allocatable       :: bspx,bspy,bspz

! FFT workspace arrays

  Real( Kind = wp ),    Dimension( :,:,: ), Allocatable, Save :: qqc
  Complex( Kind = wp ), Dimension( :,:,: ), Allocatable, Save :: qqq


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

! allocate the global charge arrays (NOT deallocated manually)

     Allocate (qqc(1:kmaxa,1:kmaxb,1:kmaxc), Stat = fail(1))
     Allocate (qqq(1:kmaxa,1:kmaxb,1:kmaxc), Stat = fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'qqq,qqc workspace arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! calculate self-interaction correction

     engsic=0.0_wp
     Do i=1,natms
        engsic=engsic+chge(i)**2
     End Do
     If (mxnode > 1) Call gsum(engsic)
     engsic=-r4pie0/epsq * alpha*engsic/sqrpi
  End If

  limit=Max(mxbuff,kmaxa*kmaxb*kmaxc) ! Try 3*kmaxa*kmaxb*kmaxc
  Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),buffer(1:limit),            Stat = fail(1))
  Allocate (ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms),it(1:mxatms),               Stat = fail(2))
  Allocate (bsdx(1:mxspl,1:mxatms),bsdy(1:mxspl,1:mxatms),bsdz(1:mxspl,1:mxatms), Stat = fail(3))
  Allocate (bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms), Stat = fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_forces allocation failure, node: ', idnode
     Call error(0)
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
! bounding is unsafe!!!  However, it doesn't matter when the 3D
! charge density matrix is global as is the case here!

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

  qqc = 0.0_wp

! construct 3D charge array

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
!              bb2=bb3*bspz(l,i)
!
!              Do k=1,mxspl
!                 kk=iyy(i)-k+2
!
!                 If (kk >= iyb .and. kk <= iyt) Then
!                    bb1=bb2*bspy(k,i)
!
!                    Do j=1,mxspl
!                       jj=ixx(i)-j+2
!
!                       If (jj >= ixb .and. jj <= ixt) Then
!                          det=bb1*bspx(j,i)
!
!                          qqc(jj,kk,ll)=qqc(jj,kk,ll)+det
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

        Do ll = llb, llt
           l = izz(i) - ll + 2

           bb2=bb3*bspz(l,i)

           Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              bb1=bb2*bspy(k,i)

              Do jj = jjb, jjt
                 j = ixx(i) - jj + 2

                 det=bb1*bspx(j,i)

                 qqc(jj,kk,ll)=qqc(jj,kk,ll)+det
              End Do
           End Do
        End Do
     End If

  End Do

! Sum up qqc

  If (mxnode > 1) Then
     i=kmaxa*kmaxb*kmaxc
     buffer(1:i) = Reshape( qqc, (/i/) )

     Call gsum(buffer(1:i))

     qqc = Reshape( buffer(1:kmaxa*kmaxb*kmaxc), (/kmaxa,kmaxb,kmaxc/) )
  End If

! load charge array into complex array for FFT

  qqq=Cmplx(qqc , Kind = wp)

! calculate inverse 3D FFT of charge array (in place)

  inc2=kmaxa
  inc3=kmaxa*kmaxb

  Call Dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmaxa,kmaxb,kmaxc,-1,1.0_wp,buffer,limit)

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

           If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

              vterm=bb1*Exp(ralph*rksq)/rksq*qqq(j,k,l)
              akv=2.0_wp*(1.0_wp/rksq-ralph)*Real( vterm*Conjg(qqq(j,k,l)),wp )
              strs(1)=strs(1)-rkx3*rkx3*akv
              strs(5)=strs(5)-rky3*rky3*akv
              strs(9)=strs(9)-rkz3*rkz3*akv
              strs(2)=strs(2)-rkx3*rky3*akv
              strs(3)=strs(3)-rkx3*rkz3*akv
              strs(6)=strs(6)-rky3*rkz3*akv
              qqq(j,k,l)=vterm

           Else

              qqq(j,k,l)=(0.0_wp,0.0_wp)

           End If

        End Do

     End Do

  End Do

! complete strs

  strs(4) = strs(2)
  strs(7) = strs(3)
  strs(8) = strs(6)

! scale strs and distribute per node

  strs = strs * scale / Real(mxnode,wp)

  Call Dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmaxa,kmaxb,kmaxc,1,1.0_wp,buffer,limit)

! calculate atomic energy

  eng = Real(Sum(qqq*qqc),wp)

! scale eng and distribute per node

  eng = eng * scale / Real(mxnode,wp)

! Second part of the monopole contribution to the stress tensor
! calculate stress tensor (symmetrical, per node)

  strs   = strs + eng*uni
  stress = stress + strs

! distribute energy and virial terms (per node)

  engcpe_rc = eng + engsic / Real(mxnode,wp)
  vircpe_rc = -(strs(1)+strs(5)+strs(9))

! infrequent calculations copying

  If (l_cp) Then
     e_rc=engcpe_rc
     v_rc=vircpe_rc
     s_rc=strs
  End If

! calculate atomic forces

  Call spme_forces(rcell,scale, ixx,iyy,izz, bspx,bspy,bspz, bsdx,bsdy,bsdz, qqq)

  Deallocate (ixx,iyy,izz,it, Stat = fail(1))
  Deallocate (bsdx,bsdy,bsdz, Stat = fail(2))
  Deallocate (bspx,bspy,bspz, Stat = fail(3))
  Deallocate (buffer,         Stat = fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

Contains

  Subroutine spme_forces(rcell,scale, ixx,iyy,izz, bspx,bspy,bspz, bsdx,bsdy,bsdz, qqq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic forces in a periodic
! system using smooth particle mesh ewald method (fourier part)
!
! copyright - daresbury laboratory
! author    - w.smith october 1998
! amended   - i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use configuration, Only : fxx,fyy,fzz

    Implicit None

    Integer,              Intent( In    ) :: ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms)
    Real( Kind = wp ),    Intent( In    ) :: scale,rcell(1:9),                &
        bsdx(1:mxspl,1:mxatms),bsdy(1:mxspl,1:mxatms),bsdz(1:mxspl,1:mxatms), &
        bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms)

    Complex( Kind = wp ), Intent( In    ) :: qqq(1:kmaxa,1:kmaxb,1:kmaxc)

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: kmaxa_r,kmaxb_r,kmaxc_r

    Integer           :: i,j,k,l, jj,kk,ll
    Real( Kind = wp ) :: tmp,facx,facy,facz,fff(0:3),fx,fy,fz,fix,fiy,fiz,qsum, &
                         bdxl,bdyl,bdzl,bdxk,bdyk,bdzk,bdxj,bdyj,bdzj

! Real values of kmax vectors

    If (newjob) Then
       newjob = .false.

       kmaxa_r=Real(kmaxa,wp)
       kmaxb_r=Real(kmaxb,wp)
       kmaxc_r=Real(kmaxc,wp)
    End If

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

             If (ll < 1) ll=ll+kmaxc

             bdxl=tmp*facx*bspz(l,i)
             bdyl=tmp*facy*bspz(l,i)
             bdzl=tmp*facz*bsdz(l,i)

             Do k=1,mxspl
                kk=iyy(i)-k+2

                If (kk < 1) kk=kk+kmaxb

                bdxk=bdxl*bspy(k,i)
                bdyk=bdyl*bsdy(k,i)
                bdzk=bdzl*bspy(k,i)

                Do j=1,mxspl
                   jj=ixx(i)-j+2

                   If (jj < 1) jj=jj+kmaxa

                   qsum=Real(qqq(jj,kk,ll),wp)

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

  End Subroutine spme_forces

End Subroutine ewald_spme_forces
