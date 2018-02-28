Subroutine ewald_real_mforces_d &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe_rl,vircpe_rl,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using multipoles with the ewald real space
! kernel
!
! copyright - daresbury laboratory
! author    - h.a.boateng february 2014
! amended   - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode,gtime
  Use setup_module
  Use config_module, Only : natms,ltg,list,fxx,fyy,fzz
  Use mpoles_module
  Use ewald_module,  Only : engsic

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_rl,vircpe_rl
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: drewd,rdrewd,alp2,co1,co2,co3,co4,co5,exclcoef, &
                             twzz,twtwz,fozz

  Integer           :: fail,idi,jatm,k,m

  Real( Kind = wp ) :: scl,rsq,rrr,engmpl,fix,fiy,fiz,fx,fy,fz,  &
                       strs1,strs2,strs3,strs5,strs6,strs9,      &
                       ppp,vk0,vk1,vk2,t1,t2,tx0,ty0,tz0,        &
                       siceng,tmpi,tmpj,tix,talp2,               &
                       tiy,tiz,tjx,tjy,tjz,b0,b1,b2,b3,b4,b5,    &
                       alpsqrpi,exparr,ecc,ecd,ecq,edd,          &
                       edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,      &
                       zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,      &
                       tm9,tm10,tm11,tm12,tm13,tm14,tm15,        &
                       dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,     &
                       dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,       &
                       dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,            &
                       dpxzzz,dpxxxx,dpyyyy,dpzzzz,              &
                       dpxxxy,dpxxxz,dpxxyy,dpxxzz,dpxxyz,       &
                       dpxyyy,dpyzzz,dpyyzz,dpxyyz,dpxyzz,       &
                       dpyyyz,dpxxxxx,dpyyyyy,dpzzzzz,dpxxxxy,   &
                       dpxxxxz,dpyyyyz,dpxzzzz,dpyzzzz,dpxyyyy,  &
                       dpxxxyy,dpxxxzz,dpyyyzz,dpyyzzz,dpxxyyy,  &
                       dpxxzzz,dpxyyzz,dpxxyyz,dpxxyzz,dpxyyyz,  &
                       dpxyzzz,dpxxxyz,xyz,                      &
                       thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,      &
                       imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,  &
                       imp9,imp10,jmp1,jmp2,jmp3,jmp4,jmp5,      &
                       jmp6,jmp7,jmp8,jmp9,jmp10,impx1,impx2,    &
                       impx3,impx4,impx5,impx6,impx7,impx8,      &
                       impx9,impx10,jmpx1,jmpx2,jmpx3,jmpx4,     &
                       jmpx5,jmpx6,jmpx7,jmpx8,jmpx9,jmpx10,     &
                       impy1,impy2,impy3,impy4,impy5,impy6,      &
                       impy7,impy8,impy9,impy10,jmpy1,jmpy2,     &
                       jmpy3,jmpy4,jmpy5,jmpy6,jmpy7,jmpy8,      &
                       jmpy9,jmpy10,impz1,impz2,impz3,impz4,     &
                       impz5,impz6,impz7,impz8,impz9,impz10,     &
                       jmpz1,jmpz2,jmpz3,jmpz4,jmpz5,jmpz6,      &
                       jmpz7,jmpz8,jmpz9,jmpz10,                 &
                       tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
  Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
  Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  If (newjob) Then
     newjob = .false.

     fail=0
     Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_real_mforces allocation failure, node: ', idnode
        Call error(0)
     End If

! interpolation interval

     drewd = rcut/Real(mxgele-4,wp)

! reciprocal of interpolation interval

     rdrewd = 1.0_wp/drewd

! generate error function complement tables for ewald sum

     Call erfcgen(rcut,alpha,mxgele,erc,fer)

! coefficients for exponential in recurrence relation

     talp2 = 2.0_wp*alpha*alpha
     alpsqrpi = 1.0_wp/(alpha*sqrpi)

     co1 = talp2*alpsqrpi
     co2 = talp2*co1
     co3 = talp2*co2
     co4 = talp2*co3
     co5 = talp2*co4

     alp2 = alpha*alpha

     exclcoef = r4pie0*alpha /sqrpi/epsq

     twzz=-2.0_wp*alpha**3 *r4pie0/(3.0_wp*sqrpi*epsq)
     twtwz=4.0_wp*alpha**5 *r4pie0/(5.0_wp*sqrpi*epsq)
     fozz=12.0_wp*alpha**5 *r4pie0/(5.0_wp*sqrpi*epsq)
  End If

! initialise potential energy and virial

  engcpe_rl=0.0_wp
  vircpe_rl=0.0_wp
  siceng   =0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity of iatm

  idi=ltg(iatm)

! get the multipoles for site i

  imp=mplgfr(:,iatm)

  If (mxompl > 0 .and. induce) Then

     imp(2)=imp(2)+indipx(iatm)
     imp(3)=imp(3)+indipy(iatm)
     imp(4)=imp(4)+indipz(iatm)

  End If

  tx0 = 0.0_wp ; ty0 = 0.0_wp ; tz0 = 0.0_wp

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) > zero_plus) Then

! get the components for site i infinitesimal rotations

     impx=mprotx(:,iatm)
     impy=mproty(:,iatm)
     impz=mprotz(:,iatm)

     imp1=imp(1); impx1=impx(1); impy1=impy(1); impz1=impz(1)

! self-interaction

     siceng = siceng+exclcoef*imp1*imp1

     If (mxompl >= 1) Then

        imp2=imp(2); imp3=imp(3); imp4=imp(4)

        impx2=impx(2); impx3=impx(3); impx4=impx(4)
        impy2=impy(2); impy3=impy(3); impy4=impy(4)
        impz2=impz(2); impz3=impz(3); impz4=impz(4)

        siceng = siceng - twzz*(imp2*imp2+imp3*imp3+imp4*imp4)
!        tx0    = tx0    - twzz*(impx2*imp2+impx3*imp3+impx4*imp4)
!        ty0    = ty0    - twzz*(impy2*imp2+impy3*imp3+impy4*imp4)
!        tz0    = tz0    - twzz*(impz2*imp2+impz3*imp3+impz4*imp4)

     End If

     If (mxompl == 2) Then

        imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

        impx5=impx(5); impx6=impx(6); impx7=impx(7); impx8=impx(8); impx9=impx(9); impx10=impx(10)
        impy5=impy(5); impy6=impy(6); impy7=impy(7); impy8=impy(8); impy9=impy(9); impy10=impy(10)
        impz5=impz(5); impz6=impz(6); impz7=impz(7); impz8=impz(8); impz9=impz(9); impz10=impz(10)

        siceng = siceng + 2.0_wp*twzz*(imp1*(imp5+imp8+imp10)) + twtwz*(2.0_wp*imp5*(imp8+imp10)+ &
                 2.0_wp*imp8*imp10+imp6*imp6+imp7*imp7+imp9*imp9) + fozz*(imp5*imp5+imp8*imp8+    &
                 imp10*imp10)
!        tx0    = tx0 + twzz*imp1*(impx5+impx8+impx10)+twtwz*(impx5*(imp8+imp10)+impx6*imp6 +     &
!                 impx7*imp7+impx8*(imp5+imp10)+impx9*imp9+impx10*(imp5+imp8))+fozz*(impx5*imp5 + &
!                 impx8*imp8+impx10*imp10)
!        ty0    = ty0 + twzz*imp1*(impy5+impy8+impy10)+twtwz*(impy5*(imp8+imp10)+impy6*imp6 +     &
!                 impy7*imp7+impy8*(imp5+imp10)+impy9*imp9+impy10*(imp5+imp8))+fozz*(impy5*imp5 + &
!                 impy8*imp8+impy10*imp10)
!        tz0    = tz0 + twzz*imp1*(impz5+impz8+impz10)+twtwz*(impz5*(imp8+imp10)+impz6*imp6 +     &
!                 impz7*imp7+impz8*(imp5+imp10)+impz9*imp9+impz10*(imp5+imp8))+fozz*(impz5*imp5 + &
!                 impz8*imp8+impz10*imp10)

     End If

! multipole scaler

     scl=r4pie0/epsq

! rescale multipoles

     imp=imp*scl

     imp1=imp(1)

     If (mxompl >= 1) Then

        imp2=imp(2); imp3=imp(3); imp4=imp(4)

     End If

     If (mxompl == 2) Then

        imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

     End If

! add self-interaction energy

     engcpe_rl = engcpe_rl - siceng
     engsic    = engsic    - siceng

! initialize torques for atom i (temporary)

     tix = .0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

!     tix = tix-tx0
!     tiy = tiy-ty0
!     tiz = tiz-tz0

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

! start of primary loop for forces evaluation

     Do m=1,list(0,iatm)

! atomic index

        jatm=list(m,iatm)

! get the multipoles for site j and the components for its infinitesimal rotations

        jmp=mplgfr(:,jatm)

        If (mxompl > 0 .and. induce) Then

           jmp(2)=jmp(2)+indipx(jatm)
           jmp(3)=jmp(3)+indipy(jatm)
           jmp(4)=jmp(4)+indipz(jatm)

        End If

        jmpx=mprotx(:,jatm)
        jmpy=mproty(:,jatm)
        jmpz=mprotz(:,jatm)

        jmp1=jmp(1); jmpx1=jmpx(1); jmpy1=jmpy(1); jmpz1=jmpz(1)

        If (mxompl >= 1) Then

           jmp2=jmp(2); jmp3=jmp(3); jmp4=jmp(4)

           jmpx2=jmpx(2); jmpx3=jmpx(3); jmpx4=jmpx(4)
           jmpy2=jmpy(2); jmpy3=jmpy(3); jmpy4=jmpy(4)
           jmpz2=jmpz(2); jmpz3=jmpz(3); jmpz4=jmpz(4)

        End If

        If (mxompl == 2) Then

           jmp5=jmp(5); jmp6=jmp(6); jmp7=jmp(7); jmp8=jmp(8); jmp9=jmp(9); jmp10=jmp(10)

           jmpx5=jmpx(5); jmpx6=jmpx(6); jmpx7=jmpx(7); jmpx8=jmpx(8); jmpx9=jmpx(9); jmpx10=jmpx(10)
           jmpy5=jmpy(5); jmpy6=jmpy(6); jmpy7=jmpy(7); jmpy8=jmpy(8); jmpy9=jmpy(9); jmpy10=jmpy(10)
           jmpz5=jmpz(5); jmpz6=jmpz(6); jmpz7=jmpz(7); jmpz8=jmpz(8); jmpz9=jmpz(9); jmpz10=jmpz(10)

        End If

! interatomic distance

        rrr = rrt(m)

! truncation of potential

        If (Maxval(Abs(jmp)) > zero_plus .and. rrr < rcut) Then

! Squared distance

           rsq=rrr**2

           engmpl = 0.0_wp
           ecc = 0.0_wp; ecd = 0.0_wp; ecq = 0.0_wp
           edd = 0.0_wp; edq = 0.0_wp; eqq = 0.0_wp
           fx = 0.0_wp; fy = 0.0_wp; fz = 0.0_wp
           tjx = 0.0_wp; tjy = 0.0_wp; tjz = 0.0_wp

           xx=xxt(m); yy=yyt(m); zz=zzt(m)
           xx2=xx*xx; yy2=yy*yy; zz2=zz*zz

! evaluate the exponential

           exparr = Exp(-alp2*rsq)

! get the value of the kernel using 3pt interpolation

           k   = Int(rrr*rdrewd)
           ppp = rrr*rdrewd - Real(k,wp)

           vk0 = erc(k)
           vk1 = erc(k+1)
           vk2 = erc(k+2)

           t1 = vk0 + (vk1 - vk0)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

! compute recurrence terms

           b0 = (t1 + (t2-t1)*ppp*0.5_wp)
           b1 = (b0        + co1*exparr)/rsq
           b2 = (3.0_wp*b1 + co2*exparr)/rsq
           b3 = (5.0_wp*b2 + co3*exparr)/rsq
           b4 = (7.0_wp*b3 + co4*exparr)/rsq
           b5 = (9.0_wp*b4 + co5*exparr)/rsq

! charge-charge interaction

           ijmp= imp(1)*jmp(1)

           ecc = ijmp*b0

           bijmp = ijmp*b1
           fx  = bijmp*xx
           fy  = bijmp*yy
           fz  = bijmp*zz

! There is no torque for charges

           If (mxompl == 1) Then

! charge-dipole interactions

              tm1 = imp2*jmp1-imp1*jmp2
              tm2 = imp3*jmp1-imp1*jmp3
              tm3 = imp4*jmp1-imp1*jmp4

              tmpi= jmp1*b1
              tmpj= imp1*b1

              dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
              dpxx= xx2*b2-b1
              dpyy= yy2*b2-b1
              dpzz= zz2*b2-b1
              dpxy= xx*yy*b2
              dpxz= xx*zz*b2
              dpyz= yy*zz*b2

              ecd = tm1*dpx + tm2*dpy + tm3*dpz

              fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
              fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
              fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

! dipole-dipole interactions

              dpxxx = xx*(3.0_wp*b2-xx2*b3)
              dpyyy = yy*(3.0_wp*b2-yy2*b3)
              dpzzz = zz*(3.0_wp*b2-zz2*b3)
              dpxxy = yy*(b2-xx2*b3)
              dpxxz = zz*(b2-xx2*b3)
              dpxyy = xx*(b2-yy2*b3)
              dpyyz = zz*(b2-yy2*b3)
              dpxzz = xx*(b2-zz2*b3)
              dpyzz = yy*(b2-zz2*b3)
              dpxyz = -xx*yy*zz*b3

              tm1 = imp2*jmp2
              tm2 = imp2*jmp3+imp3*jmp2
              tm3 = imp2*jmp4+imp4*jmp2
              tm4 = imp3*jmp3
              tm5 = imp3*jmp4+imp4*jmp3
              tm6 = imp4*jmp4

              edd = -(tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz)

              fx = fx + (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
              fy = fy + (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
              fz = fz + (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

! Torques

              tmp2 = jmp1*dpx-jmp2*dpxx-jmp3*dpxy-jmp4*dpxz
              tmp3 = jmp1*dpy-jmp2*dpxy-jmp3*dpyy-jmp4*dpyz
              tmp4 = jmp1*dpz-jmp2*dpxz-jmp3*dpyz-jmp4*dpzz

              tix  = tix + impx2*tmp2 + impx3*tmp3 + impx4*tmp4
              tiy  = tiy + impy2*tmp2 + impy3*tmp3 + impy4*tmp4
              tiz  = tiz + impz2*tmp2 + impz3*tmp3 + impz4*tmp4

              tmp2 = imp1*dpx-imp2*dpxx-imp3*dpxy-imp4*dpxz
              tmp3 = imp1*dpy-imp2*dpxy-imp3*dpyy-imp4*dpyz
              tmp4 = imp1*dpz-imp2*dpxz-imp3*dpyz-imp4*dpzz

              tjx  = tjx + jmpx2*tmp2 + jmpx3*tmp3 + jmpx4*tmp4
              tjy  = tjy + jmpy2*tmp2 + jmpy3*tmp3 + jmpy4*tmp4
              tjz  = tjz + jmpz2*tmp2 + jmpz3*tmp3 + jmpz4*tmp4

           End If

           If (mxompl == 2) Then

! charge-dipole interactions

              tm1 = imp2*jmp1-imp1*jmp2
              tm2 = imp3*jmp1-imp1*jmp3
              tm3 = imp4*jmp1-imp1*jmp4

              tmpi= jmp1*b1
              tmpj= imp1*b1

              dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
              dpxx= xx2*b2-b1
              dpyy= yy2*b2-b1
              dpzz= zz2*b2-b1
              dpxy= xx*yy*b2
              dpxz= xx*zz*b2
              dpyz= yy*zz*b2

              ecd = tm1*dpx + tm2*dpy + tm3*dpz

              fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
              fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
              fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

! dipole-dipole interactions

              dpxxx = xx*(3.0_wp*b2-xx2*b3)
              dpyyy = yy*(3.0_wp*b2-yy2*b3)
              dpzzz = zz*(3.0_wp*b2-zz2*b3)
              dpxxy = yy*(b2-xx2*b3)
              dpxxz = zz*(b2-xx2*b3)
              dpxyy = xx*(b2-yy2*b3)
              dpyyz = zz*(b2-yy2*b3)
              dpxzz = xx*(b2-zz2*b3)
              dpyzz = yy*(b2-zz2*b3)
              dpxyz = -xx*yy*zz*b3

              tm1 = imp2*jmp2
              tm2 = imp2*jmp3+imp3*jmp2
              tm3 = imp2*jmp4+imp4*jmp2
              tm4 = imp3*jmp3
              tm5 = imp3*jmp4+imp4*jmp3
              tm6 = imp4*jmp4

              edd = -(tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz)

! The forces for dipole-dipole interactions will be computed as part of the charge-quadrupole interactions

! charge-quadrupole interactions

              tm1 = imp1*jmp5+imp5*jmp1-tm1; tm2=imp1*jmp6+imp6*jmp1-tm2
              tm3 = imp1*jmp7+imp7*jmp1-tm3; tm4=imp1*jmp8+imp8*jmp1-tm4
              tm5 = imp1*jmp9+imp9*jmp1-tm5; tm6=imp1*jmp10+imp10*jmp1-tm6

              ecq = tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz-edd

              fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
              fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
              fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

! dipole-quadrupole interactions

              thrb2  = 3.0_wp*b2; thrb3 = 3.0_wp*b3

              dpxxxx = thrb2 - xx2*(6.0_wp*b3 - xx2*b4)
              dpyyyy = thrb2 - yy2*(6.0_wp*b3 - yy2*b4)
              dpzzzz = thrb2 - zz2*(6.0_wp*b3 - zz2*b4)
              dpxxxy = yy*xx*(xx2*b4-thrb3)
              dpxxxz = zz*xx*(xx2*b4-thrb3)
              dpxyyy = xx*yy*(yy2*b4-thrb3)
              dpyyyz = yy*zz*(yy2*b4-thrb3)
              dpyzzz = yy*zz*(zz2*b4-thrb3)
              dpxzzz = xx*zz*(zz2*b4-thrb3)
              dpxxyy = b2-(xx2+yy2)*b3+xx2*yy2*b4
              dpxxzz = b2-(xx2+zz2)*b3+xx2*zz2*b4
              dpyyzz = b2-(yy2+zz2)*b3+yy2*zz2*b4
              dpxxyz = yy*zz*(xx2*b4-b3)
              dpxyyz = xx*zz*(yy2*b4-b3)
              dpxyzz = xx*yy*(zz2*b4-b3)

              tm1 = imp2*jmp5-imp5*jmp2
              tm2 = imp2*jmp6-imp6*jmp2+imp3*jmp5-imp5*jmp3
              tm3 = imp2*jmp7-imp7*jmp2+imp4*jmp5-imp5*jmp4
              tm4 = imp2*jmp8-imp8*jmp2+imp3*jmp6-imp6*jmp3
              tm5 = imp2*jmp9-imp9*jmp2+imp3*jmp7-imp7*jmp3+imp4*jmp6-imp6*jmp4
              tm6 = imp2*jmp10-imp10*jmp2+imp4*jmp7-imp7*jmp4
              tm7 = imp3*jmp8-imp8*jmp3
              tm8 = imp3*jmp9-imp9*jmp3+imp4*jmp8-imp8*jmp4
              tm9 = imp3*jmp10-imp10*jmp3+imp4*jmp9-imp9*jmp4
              tm10= imp4*jmp10-imp10*jmp4

              edq = tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz + &
                    tm6*dpxzz+tm7*dpyyy+tm8*dpyyz+tm9*dpyzz+tm10*dpzzz

              fx = fx - (tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz + &
                         tm6*dpxxzz+tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz)

              fy = fy - (tm1*dpxxxy+tm2*dpxxyy+tm3*dpxxyz+tm4*dpxyyy+tm5*dpxyyz + &
                         tm6*dpxyzz+tm7*dpyyyy+tm8*dpyyyz+tm9*dpyyzz+tm10*dpyzzz)

              fz = fz - (tm1*dpxxxz+tm2*dpxxyz+tm3*dpxxzz+tm4*dpxyyz+tm5*dpxyzz + &
                         tm6*dpxzzz+tm7*dpyyyz+tm8*dpyyzz+tm9*dpyzzz+tm10*dpzzzz)

! quadrupole-quadrupole interactions

              fifb3 = 15.0_wp*b3; thrb4 = 3.0_wp*b4; sixb4 = 6.0_wp*b4; tenb4 = 10.0_wp*b4
              xyz   = xx*yy*zz

              dpxxxxx = xx*(xx2*(tenb4-xx2*b5)-fifb3)
              dpyyyyy = yy*(yy2*(tenb4-yy2*b5)-fifb3)
              dpzzzzz = zz*(zz2*(tenb4-zz2*b5)-fifb3)
              dpxxxxy = yy*(xx2*(sixb4-xx2*b5)-thrb3)
              dpxxxxz = zz*(xx2*(sixb4-xx2*b5)-thrb3)
              dpyyyyz = zz*(yy2*(sixb4-yy2*b5)-thrb3)
              dpyzzzz = yy*(zz2*(sixb4-zz2*b5)-thrb3)
              dpxzzzz = xx*(zz2*(sixb4-zz2*b5)-thrb3)
              dpxyyyy = xx*(yy2*(sixb4-yy2*b5)-thrb3)
              dpxxxyy = xx*((xx2+3.0_wp*yy2)*b4-xx2*yy2*b5-thrb3)
              dpxxxzz = xx*((xx2+3.0_wp*zz2)*b4-xx2*zz2*b5-thrb3)
              dpyyyzz = yy*((yy2+3.0_wp*zz2)*b4-yy2*zz2*b5-thrb3)
              dpyyzzz = zz*((3.0_wp*yy2+zz2)*b4-yy2*zz2*b5-thrb3)
              dpxxyyy = yy*((3.0_wp*xx2+yy2)*b4-xx2*yy2*b5-thrb3)
              dpxxzzz = zz*((3.0_wp*xx2+zz2)*b4-xx2*zz2*b5-thrb3)
              dpxyyzz = xx*((yy2+zz2)*b4-yy2*zz2*b5-b3)
              dpxxyyz = zz*((xx2+yy2)*b4-xx2*yy2*b5-b3)
              dpxxyzz = yy*((xx2+zz2)*b4-xx2*zz2*b5-b3)
              dpxyyyz = xyz*(thrb4-yy2*b5)
              dpxyzzz = xyz*(thrb4-zz2*b5)
              dpxxxyz = xyz*(thrb4-xx2*b5)

              tm1 = imp5*jmp5
              tm2 = imp5*jmp6+imp6*jmp5
              tm3 = imp5*jmp7+imp7*jmp5
              tm4 = imp5*jmp8+imp8*jmp5+imp6*jmp6
              tm5 = imp5*jmp9+imp9*jmp5+imp6*jmp7+imp7*jmp6
              tm6 = imp5*jmp10+imp10*jmp5+imp7*jmp7
              tm7 = imp6*jmp8+imp8*jmp6
              tm8 = imp6*jmp9+imp9*jmp6+imp7*jmp8+imp8*jmp7
              tm9 = imp6*jmp10+imp10*jmp6+imp7*jmp9+imp9*jmp7
              tm10= imp7*jmp10+imp10*jmp7
              tm11= imp8*jmp8
              tm12= imp8*jmp9+imp9*jmp8
              tm13= imp8*jmp10+imp10*jmp8+imp9*jmp9
              tm14= imp9*jmp10+imp10*jmp9
              tm15= imp10*jmp10

              eqq = tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz+tm6*dpxxzz +    &
                    tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz+tm11*dpyyyy+tm12*dpyyyz + &
                    tm13*dpyyzz+tm14*dpyzzz+tm15*dpzzzz

              fx = fx - (tm1*dpxxxxx+tm2*dpxxxxy+tm3*dpxxxxz+tm4*dpxxxyy+tm5*dpxxxyz+tm6*dpxxxzz +    &
                         tm7*dpxxyyy+tm8*dpxxyyz+tm9*dpxxyzz+tm10*dpxxzzz+tm11*dpxyyyy+tm12*dpxyyyz + &
                         tm13*dpxyyzz+tm14*dpxyzzz+tm15*dpxzzzz)

              fy = fy - (tm1*dpxxxxy+tm2*dpxxxyy+tm3*dpxxxyz+tm4*dpxxyyy+tm5*dpxxyyz+tm6*dpxxyzz +    &
                         tm7*dpxyyyy+tm8*dpxyyyz+tm9*dpxyyzz+tm10*dpxyzzz+tm11*dpyyyyy+tm12*dpyyyyz + &
                         tm13*dpyyyzz+tm14*dpyyzzz+tm15*dpyzzzz)

              fz = fz - (tm1*dpxxxxz+tm2*dpxxxyz+tm3*dpxxxzz+tm4*dpxxyyz+tm5*dpxxyzz+tm6*dpxxzzz +    &
                         tm7*dpxyyyz+tm8*dpxyyzz+tm9*dpxyzzz+tm10*dpxzzzz+tm11*dpyyyyz+tm12*dpyyyzz + &
                         tm13*dpyyzzz+tm14*dpyzzzz+tm15*dpzzzzz)

! second-order torques

              tmp2 = jmp1*dpx-jmp2*dpxx-jmp3*dpxy-jmp4*dpxz+jmp5*dpxxx+jmp6*dpxxy+jmp7*dpxxz+jmp8*dpxyy + &
                     jmp9*dpxyz+jmp10*dpxzz
              tmp3 = jmp1*dpy-jmp2*dpxy-jmp3*dpyy-jmp4*dpyz+jmp5*dpxxy+jmp6*dpxyy+jmp7*dpxyz+jmp8*dpyyy + &
                     jmp9*dpyyz+jmp10*dpyzz
              tmp4 = jmp1*dpz-jmp2*dpxz-jmp3*dpyz-jmp4*dpzz+jmp5*dpxxz+jmp6*dpxyz+jmp7*dpxzz+jmp8*dpyyz + &
                     jmp9*dpyzz+jmp10*dpzzz

              tmp5 = jmp1*dpxx-jmp2*dpxxx-jmp3*dpxxy-jmp4*dpxxz+jmp5*dpxxxx+jmp6*dpxxxy+jmp7*dpxxxz + &
                     jmp8*dpxxyy+jmp9*dpxxyz+jmp10*dpxxzz
              tmp6 = jmp1*dpxy-jmp2*dpxxy-jmp3*dpxyy-jmp4*dpxyz+jmp5*dpxxxy+jmp6*dpxxyy+jmp7*dpxxyz + &
                     jmp8*dpxyyy+jmp9*dpxyyz+jmp10*dpxyzz
              tmp7 = jmp1*dpxz-jmp2*dpxxz-jmp3*dpxyz-jmp4*dpxzz+jmp5*dpxxxz+jmp6*dpxxyz+jmp7*dpxxzz + &
                     jmp8*dpxyyz+jmp9*dpxyzz+jmp10*dpxzzz
              tmp8 = jmp1*dpyy-jmp2*dpxyy-jmp3*dpyyy-jmp4*dpyyz+jmp5*dpxxyy+jmp6*dpxyyy+jmp7*dpxyyz + &
                     jmp8*dpyyyy+jmp9*dpyyyz+jmp10*dpyyzz
              tmp9 = jmp1*dpyz-jmp2*dpxyz-jmp3*dpyyz-jmp4*dpyzz+jmp5*dpxxyz+jmp6*dpxyyz+jmp7*dpxyzz + &
                     jmp8*dpyyyz+jmp9*dpyyzz+jmp10*dpyzzz
              tmp10= jmp1*dpzz-jmp2*dpxzz-jmp3*dpyzz-jmp4*dpzzz+jmp5*dpxxzz+jmp6*dpxyzz+jmp7*dpxzzz + &
                      jmp8*dpyyzz+jmp9*dpyzzz+jmp10*dpzzzz

              tix = tix + impx2*tmp2+impx3*tmp3+impx4*tmp4+impx5*tmp5+impx6*tmp6+impx7*tmp7+impx8*tmp8 + &
                          impx9*tmp9+impx10*tmp10
              tiy = tiy + impy2*tmp2+impy3*tmp3+impy4*tmp4+impy5*tmp5+impy6*tmp6+impy7*tmp7+impy8*tmp8 + &
                          impy9*tmp9+impy10*tmp10
              tiz = tiz + impz2*tmp2+impz3*tmp3+impz4*tmp4+impz5*tmp5+impz6*tmp6+impz7*tmp7+impz8*tmp8 + &
                          impz9*tmp9+impz10*tmp10

              tmp2 = imp1*dpx-imp2*dpxx-imp3*dpxy-imp4*dpxz+imp5*dpxxx+imp6*dpxxy+imp7*dpxxz+imp8*dpxyy + &
                     imp9*dpxyz+imp10*dpxzz
              tmp3 = imp1*dpy-imp2*dpxy-imp3*dpyy-imp4*dpyz+imp5*dpxxy+imp6*dpxyy+imp7*dpxyz+imp8*dpyyy + &
                     imp9*dpyyz+imp10*dpyzz
              tmp4 = imp1*dpz-imp2*dpxz-imp3*dpyz-imp4*dpzz+imp5*dpxxz+imp6*dpxyz+imp7*dpxzz+imp8*dpyyz + &
                     imp9*dpyzz+imp10*dpzzz

              tmp5 = imp1*dpxx-imp2*dpxxx-imp3*dpxxy-imp4*dpxxz+imp5*dpxxxx+imp6*dpxxxy+imp7*dpxxxz + &
                     imp8*dpxxyy+imp9*dpxxyz+imp10*dpxxzz
              tmp6 = imp1*dpxy-imp2*dpxxy-imp3*dpxyy-imp4*dpxyz+imp5*dpxxxy+imp6*dpxxyy+imp7*dpxxyz + &
                     imp8*dpxyyy+imp9*dpxyyz+imp10*dpxyzz
              tmp7 = imp1*dpxz-imp2*dpxxz-imp3*dpxyz-imp4*dpxzz+imp5*dpxxxz+imp6*dpxxyz+imp7*dpxxzz + &
                     imp8*dpxyyz+imp9*dpxyzz+imp10*dpxzzz
              tmp8 = imp1*dpyy-imp2*dpxyy-imp3*dpyyy-imp4*dpyyz+imp5*dpxxyy+imp6*dpxyyy+imp7*dpxyyz + &
                     imp8*dpyyyy+imp9*dpyyyz+imp10*dpyyzz
              tmp9 = imp1*dpyz-imp2*dpxyz-imp3*dpyyz-imp4*dpyzz+imp5*dpxxyz+imp6*dpxyyz+imp7*dpxyzz + &
                     imp8*dpyyyz+imp9*dpyyzz+imp10*dpyzzz
              tmp10= imp1*dpzz-imp2*dpxzz-imp3*dpyzz-imp4*dpzzz+imp5*dpxxzz+imp6*dpxyzz+imp7*dpxzzz + &
                     imp8*dpyyzz+imp9*dpyzzz+imp10*dpzzzz

              tjx = tjx + jmpx2*tmp2+jmpx3*tmp3+jmpx4*tmp4+jmpx5*tmp5+jmpx6*tmp6+jmpx7*tmp7+jmpx8*tmp8 + &
                          jmpx9*tmp9+jmpx10*tmp10
              tjy = tjy + jmpy2*tmp2+jmpy3*tmp3+jmpy4*tmp4+jmpy5*tmp5+jmpy6*tmp6+jmpy7*tmp7+jmpy8*tmp8 + &
                          jmpy9*tmp9+jmpy10*tmp10
              tjz = tjz + jmpz2*tmp2+jmpz3*tmp3+jmpz4*tmp4+jmpz5*tmp5+jmpz6*tmp6+jmpz7*tmp7+jmpz8*tmp8 + &
                          jmpz9*tmp9+jmpz10*tmp10

           End If

           engmpl = ecc+ecd+ecq+edd+edq+eqq

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

              mptrqx(jatm)=mptrqx(jatm)+tjx
              mptrqy(jatm)=mptrqy(jatm)+tjy
              mptrqz(jatm)=mptrqz(jatm)+tjz

           End If

           If (jatm <= natms .or. idi < ltg(jatm)) Then

! accumulate potential energy

              engcpe_rl = engcpe_rl + engmpl

! calculate virial

              vircpe_rl = vircpe_rl - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

! calculate stress tensor

              strs1 = strs1 + xxt(m)*fx
              strs2 = strs2 + xxt(m)*fy
              strs3 = strs3 + xxt(m)*fz
              strs5 = strs5 + yyt(m)*fy
              strs6 = strs6 + yyt(m)*fz
              strs9 = strs9 + zzt(m)*fz

           End If

        End If

     End Do

! load back forces

     fxx(iatm)=fix
     fyy(iatm)=fiy
     fzz(iatm)=fiz

! and torques due to multipoles

     mptrqx(iatm)=mptrqx(iatm)+scl*tix
     mptrqy(iatm)=mptrqy(iatm)+scl*tiy
     mptrqz(iatm)=mptrqz(iatm)+scl*tiz

! complete stress tensor

     stress(1) = stress(1) + strs1
     stress(2) = stress(2) + strs2
     stress(3) = stress(3) + strs3
     stress(4) = stress(4) + strs2
     stress(5) = stress(5) + strs5
     stress(6) = stress(6) + strs6
     stress(7) = stress(7) + strs3
     stress(8) = stress(8) + strs6
     stress(9) = stress(9) + strs9

  End If

End Subroutine ewald_real_mforces_d
