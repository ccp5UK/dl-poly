Subroutine bspcoe(nospl,kmax1,kmax2,kmax3,csp,bscx,bscy,bscz,ww1,ww2,ww3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to calculate B-spline coefficients for
! Euler exponential splines
!
! copyright - daresbury laboratory
! author    - w.smith july 1998
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module

  Implicit None

  Integer,                 Intent( In    ) :: nospl,kmax1,kmax2,kmax3
  Real( Kind = wp ),       Intent(   Out ) :: csp(1:mxspl)
  Complex( Kind = wp ),    Intent( In    ) :: ww1(1:kmaxa),ww2(1:kmaxb),ww3(1:kmaxc)
  Complex( Kind = wp ),    Intent(   Out ) :: bscx(1:kmaxa),bscy(1:kmaxb),bscz(1:kmaxc)

  Integer              :: i,j,k
  Complex( Kind = wp ) :: ccc

! calculate B-splines at knots

  csp(1)=0.0_wp
  csp(2)=1.0_wp

  Do k=3,nospl
     csp(k)=0.0_wp

     Do j=k,2,-1
        csp(j)=(Real(j-1,wp)*csp(j)+Real(k-j+1,wp)*csp(j-1))/Real(k-1,wp)
     End Do
  End Do

! calculate B-spline coefficients

  Do i=0,kmax1-1
     ccc=(0.0_wp,0.0_wp)

     Do k=0,nospl-2
        ccc=ccc+csp(k+2)*ww1(Mod(i*k,kmax1)+1)
     End Do

     bscx(i+1)=ww1(Mod(i*(nospl-1),kmax1)+1)/ccc
  End Do

  Do i=0,kmax2-1
     ccc=(0.0_wp,0.0_wp)

     Do k=0,nospl-2
        ccc=ccc+csp(k+2)*ww2(Mod(i*k,kmax2)+1)
     End Do

     bscy(i+1)=ww2(Mod(i*(nospl-1),kmax2)+1)/ccc
  End Do

  Do i=0,kmax3-1
     ccc=(0.0_wp,0.0_wp)

     Do k=0,nospl-2
        ccc=ccc+csp(k+2)*ww3(Mod(i*k,kmax3)+1)
     End Do

     bscz(i+1)=ww3(Mod(i*(nospl-1),kmax3)+1)/ccc
  End Do

End Subroutine bspcoe

Subroutine bspgen(natms,nospl,xxx,yyy,zzz,bspx,bspy,bspz,bsdx,bsdy,bsdz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to calculate B-splines for SPME method
!
! copyright - daresbury laboratory
! author    - w.smith july 1998
! amended   - i.t.todorov april 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode
  Use setup_module

  Implicit None

  Integer,                                  Intent( In    ) :: natms,nospl
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: xxx,yyy,zzz

  Real( Kind = wp ), Dimension( 1:mxspl , 1:mxatms ), Intent(   Out ) :: &
                                             bsdx,bsdy,bsdz,bspx,bspy,bspz
  Integer           :: fail,i,j,k
  Real( Kind = wp ) :: aaa,bbb,ccc, rix0,riy0,riz0, jm1_r,k_r,km1_rr

  Real( Kind = wp ), Dimension( : ), Allocatable :: real_no, inv_no

  fail=0
  Allocate (real_no(1:nospl),inv_no(1:nospl), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bspgen allocation failure, node: ', idnode
     Call error(0)
  End If

! construct B-splines

  Do i=1,nospl
     real_no(i) = Real(i,wp)
     inv_no(i)  = 1.0_wp / real_no(i)
  End Do

  Do i=1,natms

! initializing 2nd order B-spline
! for u where (0<u<1) and (1<u<2)

     bspx(1,i)=xxx(i)-Aint(xxx(i),wp)
     bspy(1,i)=yyy(i)-Aint(yyy(i),wp)
     bspz(1,i)=zzz(i)-Aint(zzz(i),wp)

     bspx(2,i)=1.0_wp-bspx(1,i)
     bspy(2,i)=1.0_wp-bspy(1,i)
     bspz(2,i)=1.0_wp-bspz(1,i)

! Now on to calculate order k B-spline values at k
! points where (0<u<k)

     rix0=bspx(1,i)
     riy0=bspy(1,i)
     riz0=bspz(1,i)

     bsdx(1,i)= 1.0_wp
     bsdy(1,i)= 1.0_wp
     bsdz(1,i)= 1.0_wp

     bsdx(2,i)=-1.0_wp
     bsdy(2,i)=-1.0_wp
     bsdz(2,i)=-1.0_wp

     Do k=3,nospl-1 ! Order of B-spline

        bspx(k,i)=0.0_wp
        bspy(k,i)=0.0_wp
        bspz(k,i)=0.0_wp

        k_r   =real_no(k)
        km1_rr=inv_no(k-1)

        Do j=k,2,-1 ! Compute order k B-spline at points {k,k-1,...,1}

           jm1_r=real_no(j-1)

           aaa=rix0+jm1_r
           bbb=riy0+jm1_r
           ccc=riz0+jm1_r

           bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
           bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
           bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

        End Do

        bspx(1,i)=bspx(1,i)*rix0*km1_rr
        bspy(1,i)=bspy(1,i)*riy0*km1_rr
        bspz(1,i)=bspz(1,i)*riz0*km1_rr

     End Do

! Now compute B-splines for order nospl at k points where
! (0<u<nospl)

     k=nospl

     bspx(k,i)=0.0_wp
     bspy(k,i)=0.0_wp
     bspz(k,i)=0.0_wp

     k_r   =real_no(k)
     km1_rr=inv_no(k-1)

     Do j=k,2,-1

! Derivatives of B-splines with order nospl at k-1 points

        bsdx(j,i)=bspx(j,i)-bspx(j-1,i)
        bsdy(j,i)=bspy(j,i)-bspy(j-1,i)
        bsdz(j,i)=bspz(j,i)-bspz(j-1,i)

        jm1_r=real_no(j-1)

        aaa=rix0+jm1_r
        bbb=riy0+jm1_r
        ccc=riz0+jm1_r

        bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
        bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
        bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

     End Do

     bsdx(1,i)=bspx(1,i)
     bsdy(1,i)=bspy(1,i)
     bsdz(1,i)=bspz(1,i)

     bspx(1,i)=bspx(1,i)*rix0*km1_rr
     bspy(1,i)=bspy(1,i)*riy0*km1_rr
     bspz(1,i)=bspz(1,i)*riz0*km1_rr

  End Do

  Deallocate (real_no,inv_no, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bspgen allocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine bspgen

Subroutine bspgen_mpl(natms,nospl,xxx,yyy,zzz,bspx,bspy,bspz,bsddx,bsddy,bsddz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to calculate B-splines for SPME method for
! multipolar interactions
!
! copyright - daresbury laboratory
! author    - h.a.boateng april 2014
! amended   - i.t.todorov april 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
!  Use mpoles_module, Only : ncombk

  Implicit None

  Integer,                                                      Intent( In    ) :: natms,nospl
  Real( Kind = wp ), Dimension( 1:mxatms ),                     Intent( In    ) :: xxx,yyy,zzz

  Real( Kind = wp ), Dimension( 1:mxspl , 1:mxatms ),           Intent(   Out ) :: bspx,bspy,bspz
  Real( Kind = wp ), Dimension( 0:mxspl , 1:mxspl , 1:mxatms ), Intent(   Out ) :: bsddx,bsddy,bsddz

  Integer           :: fail,i,j,k,m,n,p,r,s
  Real( Kind = wp ) :: aaa,bbb,ccc, rix0,riy0,riz0, jm1_r,k_r,km1_rr,ncombk(0:mxspl,0:mxspl) = 0.0_wp
  Real( Kind = wp ) :: tmp,tempx,tempy,tempz,pcombr

  Real( Kind = wp ), Dimension( : ), Allocatable :: real_no, inv_no, pmo_no

  fail=0
  Allocate (real_no(1:nospl),inv_no(1:nospl),pmo_no(1:nospl), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bspgen_mpl allocation failure, node: ', idnode
     Call error(0)
  End If

! initialize derivatives

  bsddx(:,:,:)=0.0_wp
  bsddy(:,:,:)=0.0_wp
  bsddz(:,:,:)=0.0_wp

! construct B-splines

  Do i=1,nospl
     real_no(i) = Real(i,wp)
     inv_no(i)  = 1.0_wp / real_no(i)
     pmo_no(i)  = Real(-1**i,wp)
  End Do

  Do i=1,natms

! initializing 2nd order B-spline
! for u where (0<u<1) and (1<u<2)

     bspx(1,i)=xxx(i)-Aint(xxx(i),wp)
     bspy(1,i)=yyy(i)-Aint(yyy(i),wp)
     bspz(1,i)=zzz(i)-Aint(zzz(i),wp)

     bspx(2,i)=1.0_wp-bspx(1,i)
     bspy(2,i)=1.0_wp-bspy(1,i)
     bspz(2,i)=1.0_wp-bspz(1,i)

! compute the (nospl-2)nd derivatives

     k=2; p=nospl-2

     Do j=1,nospl

        m=Max(0,j-k)
        n=Min(p,j-1)

        tempx=0.0_wp ; tempy=0.0_wp ; tempz=0.0_wp

        Do r=m,n
           s     = j - r
           pcombr= pmo_no(r)*ncombk(p,r)

           tempx = tempx + pcombr*bspx(s,i)
           tempy = tempy + pcombr*bspy(s,i)
           tempz = tempz + pcombr*bspz(s,i)
        End Do

        bsddx(p,j,i)=tempx
        bsddy(p,j,i)=tempy
        bsddz(p,j,i)=tempz

     End Do

! Now on to calculate order k B-spline values at k
! points where (0<u<k)

     rix0=bspx(1,i)
     riy0=bspy(1,i)
     riz0=bspz(1,i)

     Do k=3,nospl-1 ! Order of B-spline

        bspx(k,i)=0.0_wp
        bspy(k,i)=0.0_wp
        bspz(k,i)=0.0_wp

        k_r   =real_no(k)
        km1_rr=inv_no(k-1)

        Do j=k,2,-1 ! Compute order k B-spline at points {k,k-1,...,1}

           jm1_r=real_no(j-1)

           aaa=rix0+jm1_r
           bbb=riy0+jm1_r
           ccc=riz0+jm1_r

           bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
           bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
           bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

        End Do

        bspx(1,i)=bspx(1,i)*rix0*km1_rr
        bspy(1,i)=bspy(1,i)*riy0*km1_rr
        bspz(1,i)=bspz(1,i)*riz0*km1_rr

! compute the (nospl-3)rd to 1st derivatives

        p = nospl-k

        Do j=1,nospl

           m=Max(0,j-k)
           n=Min(p,j-1)

           tempx=0.0_wp ; tempy=0.0_wp ; tempz=0.0_wp

           Do r=m,n
              s     = j - r
              pcombr= pmo_no(r)*ncombk(p,r)

              tempx = tempx + pcombr*bspx(s,i)
              tempy = tempy + pcombr*bspy(s,i)
              tempz = tempz + pcombr*bspz(s,i)
           End Do

           bsddx(p,j,i)=tempx
           bsddy(p,j,i)=tempy
           bsddz(p,j,i)=tempz

        End Do

     End Do

! Now compute B-splines for order nospl at k points where
! (0<u<nospl)

     k=nospl

     bspx(k,i)=0.0_wp
     bspy(k,i)=0.0_wp
     bspz(k,i)=0.0_wp

     k_r   =real_no(k)
     km1_rr=inv_no(k-1)

     Do j=k,2,-1

! B-splines with order nospl at k-1 points

        jm1_r=real_no(j-1)

        aaa=rix0+jm1_r
        bbb=riy0+jm1_r
        ccc=riz0+jm1_r

        bspx(j,i)=(aaa*bspx(j,i)+(k_r-aaa)*bspx(j-1,i))*km1_rr
        bspy(j,i)=(bbb*bspy(j,i)+(k_r-bbb)*bspy(j-1,i))*km1_rr
        bspz(j,i)=(ccc*bspz(j,i)+(k_r-ccc)*bspz(j-1,i))*km1_rr

     End Do

     bspx(1,i)=bspx(1,i)*rix0*km1_rr
     bspy(1,i)=bspy(1,i)*riy0*km1_rr
     bspz(1,i)=bspz(1,i)*riz0*km1_rr

! Now the zeroth derivatives

     bsddx(0,:,i)=bspx(:,i)
     bsddy(0,:,i)=bspy(:,i)
     bsddz(0,:,i)=bspz(:,i)

  End Do

  Deallocate (real_no,inv_no,pmo_no, Stat=fail )
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bspgen_mpl deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine bspgen_mpl

Function Dtpbsp(s1,s2,s3,rcell,bsddx,bsddy,bsddz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Function to compute arbitrary derivatives of the product of three
! b-splines for use with multipolear interactions
!
! copyright - daresbury laboratory
! author    - h.a.boateng april 2014
! amended   - i.t.todorov march 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : imcon
!  Use mpoles_module, Only : ncombk

  Implicit None

  Real( Kind = wp ) :: Dtpbsp

  Integer,                                 Intent( In   ) :: s1,s2,s3
  Real( Kind = wp ),                       Intent( In   ) :: rcell(9)
  Real( Kind = wp ), Dimension( 0:mxspl ), Intent( In   ) :: bsddx,bsddy,bsddz

  Real( Kind = wp ) :: tx,ty,tz,sx,sy,sz,ncombk(0:mxspl,0:mxspl) = 0.0_wp
  Real( Kind = wp ) :: ka11,ka12,ka13,kb21,kb22,kb23,kc31,kc32,kc33
  Integer           :: i,j,k,j1,j2,j3,k1,k2,k3,jj,kk,sk,sk3,sk2,sk1

  Dtpbsp = 0.0_wp

  ka11 = Real(kmaxa,wp)*rcell(1)
  ka12 = Real(kmaxa,wp)*rcell(4)
  ka13 = Real(kmaxa,wp)*rcell(7)
  kb21 = Real(kmaxb,wp)*rcell(2)
  kb22 = Real(kmaxb,wp)*rcell(5)
  kb23 = Real(kmaxb,wp)*rcell(8)
  kc31 = Real(kmaxc,wp)*rcell(3)
  kc32 = Real(kmaxc,wp)*rcell(6)
  kc33 = Real(kmaxc,wp)*rcell(9)

! Typically, the box is orthogonal => only diagonals-ka11,kb22,kc33-are non-zero

  If (imcon /= 3) Then

     Dtpbsp = ka11**s1 * kb22**s2 * kc33**s3 * bsddx(s1) * bsddy(s2) * bsddz(s3)

  Else

     tz = 1.0_wp
     Do k3 = 0, s3
        ty = tz * ncombk(s3,k3); sk3=s3-k3

        Do k2 = 0, s2
           tx = ty * ncombk(s2,k2); sk2=s2-k2

           Do k1 = 0, s1
              kk = k1+k2+k3; sk1=s1-k1; sk=sk1+sk2+sk3

              sz = tx * ncombk(s1,k1)*bsddx(kk)

              Do j3 = 0, sk3
                 sy = sz * ncombk(sk3,j3)*kc33**(sk3-j3)

                 Do j2 = 0, sk2
                    sx = sy * ncombk(sk2,j2)*kc32**(sk2-j2)

                    Do j1 = 0, sk1
                       jj = j1+j2+j3

                       Dtpbsp = Dtpbsp + sx * kc31**(sk1-j1) * ncombk(sk1,j1)*bsddy(jj)*bsddz(sk-jj)

                       sx=sx*kb21
                    End Do

                    sy=sy*kb22
                 End Do

                 sz=sz*kb23
              End Do

              tx=tx*ka11
           End Do

           ty=ty*ka12
        End Do

        tz=tz*ka13
     End Do

  End If

End Function Dtpbsp

Subroutine spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to create complex exponential arrays for b-splines
!
! copyright - daresbury laboratory
! author    - w.smith october 1998
! amended   - i.t.todorov october 2006
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : twopi

  Implicit None

  Integer,              Intent( In    ) :: ndiv1,ndiv2,ndiv3
  Complex( Kind = wp ), Intent(   Out ) :: ww1(1:ndiv1),ww2(1:ndiv2),ww3(1:ndiv3)

  Integer           :: i
  Real( Kind = wp ) :: arg

! initialise complex exponential factors

  ww1(1)=(1.0_wp,0.0_wp)

  Do i=1,ndiv1/2
    arg=(twopi/Real(ndiv1,wp))*Real(i,wp)

    ww1(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
    ww1(ndiv1+1-i)=Conjg(ww1(i+1))
  End Do

  ww2(1)=(1.0_wp,0.0_wp)

  Do i=1,ndiv2/2
    arg=(twopi/Real(ndiv2,wp))*Real(i,wp)

    ww2(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
    ww2(ndiv2+1-i)=Conjg(ww2(i+1))
  End Do

  ww3(1)=(1.0_wp,0.0_wp)

  Do i=1,ndiv3/2
    arg=(twopi/Real(ndiv3,wp))*Real(i,wp)

    ww3(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
    ww3(ndiv3+1-i)=Conjg(ww3(i+1))
  End Do

End Subroutine spl_cexp

Subroutine dlpfft3(isw,ndiv1,ndiv2,ndiv3,ww1,ww2,ww3,aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 3D fast fourier transform routine (dependent upon spl_cexp)
!
! copyright - daresbury laboratory
! author    - w.smith july 1998
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode
  Use setup_module, Only : nrite

  Implicit None

  Integer,                 Intent( In    ) :: isw,ndiv1,ndiv2,ndiv3
  Complex( Kind = wp ),    Intent( InOut ) :: aaa(1:ndiv1,1:ndiv2,1:ndiv3), &
                                              ww1(1:ndiv1),ww2(1:ndiv2),ww3(1:ndiv3)

  Logical, Save :: newjob = .true.
  Integer, Save :: nu1,nu2,nu3

  Integer              :: fail,iii,jjj,kkk,i,j,k,l,jj2,kk1,k12,num
  Real( Kind = wp )    :: tmp
  Complex( Kind = wp ) :: ttt

  Integer, Dimension( : ), Allocatable, Save :: key1,key2,key3

  If (newjob) Then
     newjob = .false.

! check FFT array dimensions

     nu1 = Nint(Log(Real(ndiv1,wp))/Log(2.0_wp))
     nu2 = Nint(Log(Real(ndiv2,wp))/Log(2.0_wp))
     nu3 = Nint(Log(Real(ndiv3,wp))/Log(2.0_wp))

     tmp = Log(Real(ndiv1,wp))*Log(Real(ndiv2,wp))*Log(Real(ndiv3,wp)) / &
           (Log(2.0_wp)**3)

     If (Abs(tmp-Real(nu1*nu2*nu3,wp)) > 1.0e-6_wp) Then
        Write(nrite,'(/,1x,a,20i2)') 'error - FFT array dimension not 2^N ',ndiv1,ndiv2,ndiv3
        Call error(0)
     End If

! allocate reverse bit address arrays

     Allocate (key1(1:ndiv1),key2(1:ndiv2),key3(1:ndiv3), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'dlpfft3 (FFT reverse bit address arrays) allocation failure, node: ', idnode
        Call error(0)
     End If

! set reverse bit address arrays

     Do kkk=1,ndiv1
        iii=0
        jjj=kkk-1

        Do j=1,nu1
           jj2=jjj/2
           iii=2*(iii-jj2)+jjj
           jjj=jj2
        End Do

        key1(kkk)=iii+1
     End Do

     Do kkk=1,ndiv2
        iii=0
        jjj=kkk-1

        Do j=1,nu2
           jj2=jjj/2
           iii=2*(iii-jj2)+jjj
           jjj=jj2
        End Do

        key2(kkk)=iii+1
     End Do

     Do kkk=1,ndiv3
        iii=0
        jjj=kkk-1

        Do j=1,nu3
           jj2=jjj/2
           iii=2*(iii-jj2)+jjj
           jjj=jj2
        End Do

        key3(kkk)=iii+1
     End Do

! initialise complex exponential factors
! Call spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

  End If

! take conjugate of exponentials if required

  If (isw < 0) Then
     Do i=1,ndiv1
        ww1(i)=Conjg(ww1(i))
     End Do

     Do i=1,ndiv2
        ww2(i)=Conjg(ww2(i))
     End Do

     Do i=1,ndiv3
        ww3(i)=Conjg(ww3(i))
     End Do
  End If

! perform fourier transform in X direction

  kkk=0
  num=ndiv1/2
  Do l=1,nu1
     Do While (kkk < ndiv1)
        Do i=1,num
           iii=key1(kkk/num+1)
           kk1=kkk+1
           k12=kk1+num

           Do j=1,ndiv2
              Do k=1,ndiv3
                 ttt=aaa(k12,j,k)*ww1(iii)
                 aaa(k12,j,k)=aaa(kk1,j,k)-ttt
                 aaa(kk1,j,k)=aaa(kk1,j,k)+ttt
              End Do
           End Do

           kkk=kkk+1
        End Do

        kkk=kkk+num
     End Do

     kkk=0
     num=num/2
  End Do

! unscramble the fft using bit address array

  Do kkk=1,ndiv1
     iii=key1(kkk)

     If (iii > kkk) Then
        Do j=1,ndiv2
           Do k=1,ndiv3
              ttt=aaa(kkk,j,k)
              aaa(kkk,j,k)=aaa(iii,j,k)
              aaa(iii,j,k)=ttt
           End Do
        End Do
     End If
  End Do

! perform fourier transform in Y direction

  kkk=0
  num=ndiv2/2
  Do l=1,nu2
     Do While (kkk < ndiv2)
        Do i=1,num
           iii=key2(kkk/num+1)
           kk1=kkk+1
           k12=kk1+num

           Do j=1,ndiv1
              Do k=1,ndiv3
                 ttt=aaa(j,k12,k)*ww2(iii)
                 aaa(j,k12,k)=aaa(j,kk1,k)-ttt
                 aaa(j,kk1,k)=aaa(j,kk1,k)+ttt
              End Do
           End Do

           kkk=kkk+1
        End Do

        kkk=kkk+num
     End Do

     kkk=0
     num=num/2
  End Do

! unscramble the fft using bit address array

  Do kkk=1,ndiv2
     iii=key2(kkk)

     If (iii > kkk) Then
        Do j=1,ndiv1
           Do k=1,ndiv3
              ttt=aaa(j,kkk,k)
              aaa(j,kkk,k)=aaa(j,iii,k)
              aaa(j,iii,k)=ttt
           End Do
        End Do
     End If
  End Do

! perform fourier transform in Z direction

   kkk=0
   num=ndiv3/2
   Do l=1,nu3

      Do While (kkk < ndiv3)
         Do i=1,num
            iii=key3(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num

            Do j=1,ndiv1
               Do k=1,ndiv2
                  ttt=aaa(j,k,k12)*ww3(iii)
                  aaa(j,k,k12)=aaa(j,k,kk1)-ttt
                  aaa(j,k,kk1)=aaa(j,k,kk1)+ttt
               End Do
            End Do

            kkk=kkk+1
         End Do

         kkk=kkk+num
      End Do

      kkk=0
      num=num/2
   End Do

! unscramble the fft using bit address array

   Do kkk=1,ndiv3
      iii=key3(kkk)

      If (iii > kkk) Then
         Do j=1,ndiv1
            Do k=1,ndiv2
               ttt=aaa(j,k,kkk)
               aaa(j,k,kkk)=aaa(j,k,iii)
               aaa(j,k,iii)=ttt
           End Do
        End Do
     End If
  End Do

! restore exponentials to unconjugated values if necessary

  If (isw < 0) Then
     Do i=1,ndiv1
        ww1(i)=Conjg(ww1(i))
     End Do

     Do i=1,ndiv2
        ww2(i)=Conjg(ww2(i))
     End Do

     Do i=1,ndiv3
        ww3(i)=Conjg(ww3(i))
     End Do
  End If

End Subroutine dlpfft3
