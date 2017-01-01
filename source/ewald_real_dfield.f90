Subroutine ewald_real_dfield &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,flag)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field
! in a periodic system due to induced dipoles with the ewald real space
! kernel
!
! copyright - daresbury laboratory
! author    - h.a.boateng december 2014
! amended   - i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,gtime
  Use setup_module,  Only : mxgele
  Use config_module, Only : natms,nlast,ltg,list
  Use mpoles_module, Only : indipx,indipy,indipz,mpfldx,mpfldy,mpfldz,&
                            plratm
  Use setup_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm,flag
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: drewd,rdrewd,alp2,co1,co2,exclcoef,awidth

  Integer           :: fail,idi,jatm,k,m,limit
  Real( Kind = wp ) :: scl,rsq,rrr,enempol,fix,fiy,fiz,fx,fy,fz,  &
                       ppp,vk0,vk1,vk2,t1,t2,tensic,tx0,ty0,tz0,  &
                       tmp,tmpi,tmpj,tix,talp2,                   &
                       tiy,tiz,b0,b1,b2,b3,b4,b5,                 &
                       alpsqrpi,exparr,ecc,ecd,ecq,edd,           &
                       edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,       &
                       zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,       &
                       tm9,tm10,tm11,tm12,tm13,tm14,tm15,         &
                       dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,      &
                       dpyz,                                      &
                       thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,       &
                       imp1,imp2,imp3,jmp1,jmp2,jmp3,             &
                       tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,   &
                       p9,tmp10,damp,expdamp,scale3,scale5,iplrz, &
                       rr3,rr5

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

  If (newjob) Then
     newjob = .false.

     fail=0
     Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_real_mfield allocation failure, node: ', idnode
        Call error(0)
     End If

! width of exponential screening for Thole damping

     awidth = 0.39_wp

! interpolation interval

     drewd = rcut/Real(mxgele-4,wp)

! reciprocal ofinterpolation interval

     rdrewd = 1.0_wp/drewd

! generate error function complement tables for ewald sum

     Call erfcgen(rcut,alpha,mxgele,erc,fer)

! coefficients for exponential in recurrence relation

     talp2 = 2.0_wp*alpha*alpha; alpsqrpi = 1.0_wp/(alpha*sqrpi)

     co1 = talp2*alpsqrpi; co2 = talp2*co1

     alp2= alpha*alpha

     exclcoef = r4pie0*alpha/sqrpi/epsq
  End If

! global identity of iatm

  idi=ltg(iatm)

! dipole scaler

  scl=r4pie0/epsq

! get the dipoles for site i

  imp1=indipx(iatm)*scl
  imp2=indipy(iatm)*scl
  imp3=indipz(iatm)*scl
  iplrz=plratm(iatm)

! load field

  fix=mpfldx(iatm)
  fiy=mpfldy(iatm)
  fiz=mpfldz(iatm)

! induced dipole-induced dipole

  If (flag == 0) limit = list(0,iatm)
  If (flag == 1) limit = list(-1,iatm)-list(0,iatm)

! start of primary loop for electrostatic field evaluation

  Do m=1,limit

! interatomic distance

     rrr=rrt(m)
     rsq=rrr**2

! truncation of potential

     If (rrr < rcut) Then
        rr3 = 1.0_wp/(rrr*rsq)
        rr5 = rr3/rsq

! atomic index

        If (flag == 0) Then
           jatm = list(m,iatm)

           scale3 = 1.0_wp
           scale5 = 1.0_wp
        Else
           jatm = list(list(0,iatm)+m,iatm)

           damp    = -awidth*rrr*rsq/Sqrt(iplrz*plratm(jatm))
           expdamp = Exp(damp)
           scale3  = 1.0_wp-expdamp
           scale5  = 1.0_wp-(1.0_wp-damp)*expdamp
        End If

! get the dipoles for site j

        jmp1=indipx(jatm)*scl
        jmp2=indipy(jatm)*scl
        jmp3=indipz(jatm)*scl

! initialise field on atom i

        fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp

        xx  = xxt(m) ; yy  = yyt(m) ; zz  = zzt(m)
        xx2 = xx*xx  ; yy2 = yy*yy  ; zz2 = zz*zz

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

        b0 = (t1 + (t2-t1)*ppp*0.5_wp)
        b1 = (b0        + co1*exparr)/rsq
        b2 = (3.0_wp*b1 + co2*exparr)/rsq

! apply Thole damping

        b1 = b1 - (1.0_wp-scale3)*rr3
        b2 = b2 - 3.0_wp*(1.0_wp-scale5)*rr5

! field due to induce dipole on atom j

        dpxx = xx2*b2-b1
        dpyy = yy2*b2-b1
        dpzz = zz2*b2-b1
        dpxy = xx*yy*b2
        dpxz = xx*zz*b2
        dpyz = yy*zz*b2

! field due to the induced dipole on atom j

        If (Abs(jmp1)+Abs(jmp2)+Abs(jmp3) > zero_plus) Then
           fx = fx + (dpxx*jmp1 + dpxy*jmp2 + dpxz*jmp3)
           fy = fy + (dpxy*jmp1 + dpyy*jmp2 + dpyz*jmp3)
           fz = fz + (dpxz*jmp1 + dpyz*jmp2 + dpzz*jmp3)

           fix = fix + fx
           fiy = fiy + fy
           fiz = fiz + fz
        End If

! electrostatic field at position j due to induce dipole on atom i

        If (flag == 0) Then
           If (jatm <= natms) Then
              If (Abs(imp1)+Abs(imp2)+Abs(imp3) > zero_plus) Then
                 mpfldx(jatm) = mpfldx(jatm) + (dpxx*imp1 + dpxy*imp2 + dpxz*imp3)
                 mpfldy(jatm) = mpfldy(jatm) + (dpxy*imp1 + dpyy*imp2 + dpyz*imp3)
                 mpfldz(jatm) = mpfldz(jatm) + (dpxz*imp1 + dpyz*imp2 + dpzz*imp3)
              End If
           End If
        End If
     End If

  End Do

! load back electrostatic field

  mpfldx(iatm)=fix
  mpfldy(iatm)=fiy
  mpfldz(iatm)=fiz

End Subroutine ewald_real_dfield
