Subroutine ewald_spme_mfield_d(alpha,epsq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field
! due to multipolar interactions in a periodic system using the smooth
! particle mesh ewald method for multipoles
!
! This version allows for extension to arbitrary order
!
! Note: (fourier) reciprocal space terms by global fft sumation;
!       (1) ewald_spme_mfie~d_d.f90 is not dependent on exchange_grid.f90
!           and parallel_fft.90, which is dependent on gpfaf90.
!           Therefore, all references in Makefile to these three files
!           are not needed.
!
! copyright - daresbury laboratory
! author    - h.a.boateng and i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,   Only : idnode
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz,idx,idy,idz
  Use configuration,  Only : cell,natms,nlast,xxx,yyy,zzz
  Use mpoles_module,  Only : mplgfr
  Use ewald_module
  Use parallel_fft

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: alpha,epsq

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ixb,iyb,izb, ixt,iyt,izt
  Real( Kind = wp ), Save :: ixbm1_r,iybm1_r,izbm1_r, &
                             ixtm0_r,iytm0_r,iztm0_r, &
                             kmaxa_r,kmaxb_r,kmaxc_r

  Logical              :: llspl=.true.
  Integer              :: fail(1:4), i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt

  Real( Kind = wp )    :: det,rcell(1:9),celprp(1:10),ralph,                   &
                          rcpcut,rcpct2,tmp,bb1,bb2,bb3,                       &
                          rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3, &
                          dtp,                                                 &
                          imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10,  &
                          bdy0,bdy1,bdy2,bdz0,bdz1,bdz2,                       &
                          ka11,kb22,kc33,ka11sq,kb22sq,kc33sq,kakb,kakc,kbkc

  Real( Kind = wp )    :: imp(1:mximpl)

  Complex( Kind = wp ) :: vterm

! blocking factors for fft

  Integer,           Save :: block_x,block_y,block_z

! B-spline coefficients

  Complex( Kind = wp ), Dimension( : ),     Allocatable, Save :: bscx,bscy,bscz
  Complex( Kind = wp ), Dimension( : ),     Allocatable, Save :: ww1,ww2,ww3

  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: csp
  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: txx,tyy,tzz
  Integer,              Dimension( : ),     Allocatable       :: ixx,iyy,izz,it
  Real( Kind = wp ),    Dimension( :,:,: ), Allocatable       :: bsddx,bsddy,bsddz
  Real( Kind = wp ),    Dimension( :,: ),   Allocatable       :: bspx,bspy,bspz
  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: bdx,bdy,bdz

! context for parallel fft

  Integer,           Save :: context

! indexing arrays for x, y & z as used in parallel fft

  Integer,              Dimension( : ),     Allocatable, Save :: index_x,index_y,index_z

! temporary workspace for parallel fft

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
! allocate the complex exponential arrays (NOT deallocated manually)

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

! deallocate the helper array

     Deallocate (csp,         Stat = fail(1))
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
  End If

  Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),                            Stat = fail(1))
  Allocate (ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms),it(1:mxatms),               Stat = fail(2))
  Allocate (bdx(0:mxspl),bsddx(0:mxspl,1:mxspl,1:mxatms), &
            bdy(0:mxspl),bsddy(0:mxspl,1:mxspl,1:mxatms), &
            bdz(0:mxspl),bsddz(0:mxspl,1:mxspl,1:mxatms),                         Stat = fail(3))
  Allocate (bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms), Stat = fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mfield allocation failure, node: ', idnode
     Call error(0)
  End If

! set working parameters

  ralph=-0.25_wp/alpha**2

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

! If not DD bound in kmax grid space when .not.llvnl = (mxspl1 == mxspl)

     If (mxspl1 == mxspl .and. i <= natms) Then
        If (txx(i) < ixbm1_r .or. txx(i) > ixtm0_r .or. &
            tyy(i) < iybm1_r .or. tyy(i) > iytm0_r .or. &
            tzz(i) < izbm1_r .or. tzz(i) > iztm0_r) llspl=.false.
     End If

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

  Deallocate (txx,tyy,tzz, Stat = fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mfield allocation failure, node: ', idnode
     Call error(0)
  End If

! zero 3D charge array
! DaFT version - only need set local bit to zero

  qqc = 0.0_wp

! construct 3D charge array

  ka11 = kmaxa_r*rcell(1)
  kb22 = kmaxb_r*rcell(5)
  kc33 = kmaxc_r*rcell(9)

  ka11sq = ka11*ka11; kb22sq = kb22*kb22; kc33sq = kc33*kc33
  kakb = ka11*kb22; kakc = ka11*kc33; kbkc = kb22*kc33

! DaFT version - use array that holds only the local data

  Do i=1,nlast

! If a particle is charged and in the MD cell or in its positive halo
! (t(i) >= 0) as the B-splines are negative directionally by propagation

     If (it(i) == 1) Then

! get the multipoles for site i

        imp=mplgfr(:,i)

        imp1=imp(1)

        If (mxompl >= 1) Then
           imp2=imp(2) ; imp3=imp(3) ; imp4=imp(4)
        End If

        If (mxompl == 2) Then
           imp5=imp(5) ; imp6=imp(6) ; imp7=imp(7)
           imp8=imp(8) ; imp9=imp(9) ; imp10=imp(10)
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
           bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

           Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              Do jj = jjb, jjt
                 j = ixx(i) - jj + 2

                 bdx=bsddx(:,j,i)

                 dtp = imp1*bdx(0)*bdy0*bdz0

                 If (mxompl >= 1 ) Then
                    dtp = dtp+imp2*ka11*bdx(1)*bdy0*bdz0
                    dtp = dtp+imp3*kb22*bdx(0)*bdy1*bdz0
                    dtp = dtp+imp4*kc33*bdx(0)*bdy0*bdz1
                 End If

                 If (mxompl == 2) Then
                    dtp = dtp+imp5*ka11sq  *bdx(2)*bdy0*bdz0
                    dtp = dtp+imp8*kb22sq  *bdx(0)*bdy2*bdz0
                    dtp = dtp+imp10*kc33sq *bdx(0)*bdy0*bdz2

                    dtp = dtp+imp6*kakb*bdx(1)*bdy1*bdz0
                    dtp = dtp+imp7*kakc*bdx(1)*bdy0*bdz1
                    dtp = dtp+imp9*kbkc*bdx(0)*bdy1*bdz1
                 End If

                 qqc(jj,kk,ll) = qqc(jj,kk,ll) + dtp
              End Do
           End Do
        End Do

     End If

  End Do

! load charge array into complex array for FFT

  qqq=Cmplx(qqc , Kind = wp)

! calculate inverse 3D FFT of charge array (in place)

  Call dlpfft3(1,kmaxa,kmaxb,kmaxc,ww1,ww2,ww3,qqq)

! set reciprocal space cutoff

  Call dcell(rcell,celprp)

  rcpcut=0.5_wp*Min(kmaxa_r*celprp(7),kmaxb_r*celprp(8),kmaxc_r*celprp(9))
  rcpcut=rcpcut*1.05_wp*twopi
  rcpct2=rcpcut**2

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

              qqq(j,k,l)=vterm

           Else

              qqq(j,k,l)=(0.0_wp,0.0_wp)

           End If

        End Do

     End Do

  End Do

  Call dlpfft3(-1,kmaxa,kmaxb,kmaxc,ww1,ww2,ww3,qqq)

! calculate electrostatic field

  Call spme_mfield(rcell, ixx,iyy,izz,bsddx,bsddy,bsddz, qqq)

  Deallocate (ixx,iyy,izz,it,    Stat = fail(1))
  Deallocate (bsddx,bsddy,bsddz, Stat = fail(2))
  Deallocate (bspx,bspy,bspz,    Stat = fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mfield deallocation failure, node: ', idnode
     Call error(0)
  End If

Contains

  Subroutine spme_mfield(rcell, ixx,iyy,izz, bsddx,bsddy,bsddz, qqq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field due to
! multipolar interactions in periodic system using the smooth particle
! mesh ewald method (fourier part)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov december 2016
! amended   - h.a.boateng december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use mpoles_module, Only : mpfldx,mpfldy,mpfldz

    Implicit None

    Integer,              Intent( In    ) :: ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms)
    Real( Kind = wp ),    Intent( In    ) :: rcell(1:9),                      &
                                             bsddx(0:mxspl,1:mxspl,1:mxatms), &
                                             bsddy(0:mxspl,1:mxspl,1:mxatms), &
                                             bsddz(0:mxspl,1:mxspl,1:mxatms)

    Complex( Kind = wp ), Intent( In    ) :: qqq(1:kmaxa,1:kmaxb,1:kmaxc)

    Integer           :: fail, i,j,k,l, jj,kk,ll

    Real( Kind = wp ) :: ka11,kb22,kc33,fix,fiy,fiz,qsum

    Real( Kind = wp ), Dimension( : ), Allocatable :: bdx,bdy,bdz

    fail=0
    Allocate (bdx(0:mxspl),bdy(0:mxspl),bdz(0:mxspl), Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_mfield allocation failure, node: ', idnode
    End If

! Real values of kmax vectors

    ka11 = Real(kmaxa,wp)*rcell(1)
    kb22 = Real(kmaxb,wp)*rcell(5)
    kc33 = Real(kmaxc,wp)*rcell(9)

    Do i=1,natms

! initialise field

       fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp

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

! field
                fix=fix+qsum*ka11*bdx(1)*bdy(0)*bdz(0)
                fiy=fiy+qsum*kb22*bdx(0)*bdy(1)*bdz(0)
                fiy=fiy+qsum*kc33*bdx(0)*bdy(0)*bdz(1)
             End Do
          End Do
       End Do

! load field

       mpfldx(i)=mpfldx(i)+fix
       mpfldy(i)=mpfldy(i)+fiy
       mpfldz(i)=mpfldz(i)+fiz

    End Do

    Deallocate (bdx,bdy,bdz, Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_mfield dealocation failure, node: ', idnode
    End If

  End Subroutine spme_mfield

End Subroutine ewald_spme_mfield_d
