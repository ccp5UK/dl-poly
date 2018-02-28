Subroutine coul_cp_dfield &
           (iatm,rcut,epsq,xxt,yyt,zzt,rrt,flag)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field
! in a periodic system due to induced dipoles for the kernel 1/r
!
! copyright - daresbury laboratory
! author    - h.a.boateng december 2014
! amended   - i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode,gtime
  Use setup_module
  Use config_module, Only : natms,nlast,ltg,list
  Use mpoles_module, Only : indipx,indipy,indipz,mpfldx,mpfldy,mpfldz,&
                            plratm

  Implicit None

  Integer,                                  Intent( In    ) :: iatm,flag
  Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt

  Integer           :: fail,idi,jatm,k,m,limit
  Real( Kind = wp ) :: awidth,scl,                                &
                       rsq,rrr,enempol,fix,fiy,fiz,fx,fy,fz,      &
                       ppp,vk0,vk1,vk2,t1,t2,tensic,tx0,ty0,tz0,  &
                       tmp,tmpi,tmpj,tix,talp2,                   &
                       b0,b1,b2,b3,b4,b5,                         &
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

! width of exponential screening for Thole damping

  awidth = 0.39_wp

! dipole scaler

  scl=r4pie0/epsq

! global identity of iatm

  idi=ltg(iatm)

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
        Else
           jatm = list(list(0,iatm)+m,iatm)
        End If

        damp    = -awidth*rrr*rsq/Sqrt(iplrz*plratm(jatm))
        expdamp = Exp(damp)
        scale3  = 1.0_wp-expdamp
        scale5  = 1.0_wp-(1.0_wp-damp)*expdamp

! get the dipoles for site j

        jmp1=indipx(jatm)*scl
        jmp2=indipy(jatm)*scl
        jmp3=indipz(jatm)*scl

! initialise field on atom i

        fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp

        xx  = xxt(m) ; yy  = yyt(m) ; zz  = zzt(m)
        xx2 = xx*xx  ; yy2 = yy*yy  ; zz2 = zz*zz

! apply Thole damping

        b1 = (1.0_wp-scale3)*rr3
        b2 = 3.0_wp*(1.0_wp-scale5)*rr5

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

End Subroutine coul_cp_dfield
