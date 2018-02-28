Subroutine ewald_real_mfield_d &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field
! in a periodic system due to permanent multipoles with the ewald real space
! kernel
!
! copyright - daresbury laboratory
! author    - i.t.todorov & h.a.boateng december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode,gtime
  Use setup_module
  Use config_module, Only : natms,nlast,ltg,list
  Use mpoles_module, Only : mplgfr,mpfldx,mpfldy,mpfldz

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: drewd,rdrewd,alp2,co1,co2,co3,co4,co5,exclcoef, &
                             twzz,twtwz,fozz,engsic

  Integer           :: fail,idi,jatm,k,m,counter
  Real( Kind = wp ) :: scl,rsq,rrr,enempol,fix,fiy,fiz,fx,fy,fz, &
                       ppp,vk0,vk1,vk2,t1,t2,tensic,tx0,ty0,tz0, &
                       tix,talp2,                                &
                       tiy,tiz,tjx,tjy,tjz,b0,b1,b2,b3,b4,b5,    &
                       alpsqrpi,exparr,ecc,ecd,ecq,edd,          &
                       edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,zz2,  &
                       dpxx,dpyy,dpzz,dpxy,dpxz,                 &
                       dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,       &
                       dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,            &
                       imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,  &
                       imp9,imp10,jmp1,jmp2,jmp3,jmp4,jmp5,      &
                       jmp6,jmp7,jmp8,jmp9,jmp10

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)

  If (newjob) Then
     newjob = .false.

     fail=0
     Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_real_mfield allocation failure, node: ', idnode
        Call error(0)
     End If

! interpolation interval

     drewd = rcut/Real(mxgele-4,wp)

! reciprocal ofinterpolation interval

     rdrewd = 1.0_wp/drewd

! generate error function complement tables for ewald sum

     Call erfcgen(rcut,alpha,mxgele,erc,fer)

! coefficients for exponential in recurrence relation

     talp2    = 2.0_wp*alpha*alpha
     alpsqrpi = 1.0_wp/(alpha*sqrpi)

     co1 = talp2*alpsqrpi
     co2 = talp2*co1
     co3 = talp2*co2
     co4 = talp2*co3
     co5 = talp2*co4

     alp2= alpha*alpha

     exclcoef = r4pie0*alpha /sqrpi/epsq

     twzz=-2.0_wp*alpha**3 *r4pie0/(3.0_wp*sqrpi*epsq)
     twtwz=4.0_wp*alpha**5 *r4pie0/(5.0_wp*sqrpi*epsq)
     fozz=12.0_wp*alpha**5 *r4pie0/(5.0_wp*sqrpi*epsq)
  End If

! global identity of iatm

  idi=ltg(iatm)

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) > zero_plus) Then

! multipole scaler

     scl=r4pie0/epsq

! scale imp multipoles

     imp=mplgfr(:,iatm)*scl

     imp1=imp(1)

     If (mxompl >= 1) Then

        imp2=imp(2); imp3=imp(3); imp4=imp(4)

     End If

     If (mxompl == 2) Then

        imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

     End If

! load field

     fix=mpfldx(iatm)
     fiy=mpfldy(iatm)
     fiz=mpfldz(iatm)

! start of primary loop for electrostatic field evaluation

     Do m=1,list(0,iatm)

! atomic index

        jatm=list(m,iatm)

! get the multipoles for site j

        jmp=mplgfr(:,jatm)

! interatomic distance

        rrr = rrt(m)

! truncation of potential

        If (Maxval(Abs(jmp)) > zero_plus .and. rrr < rcut) Then

! Squared distance

           rsq=rrr**2

! scale jmp multipoles

           jmp=jmp*scl

           jmp1=jmp(1)

           If (mxompl >= 1) Then

              jmp2=jmp(2); jmp3=jmp(3); jmp4=jmp(4)

           End If

           If (mxompl == 2) Then

              jmp5=jmp(5); jmp6=jmp(6); jmp7=jmp(7); jmp8=jmp(8); jmp9=jmp(9); jmp10=jmp(10)

           End If

           fx = 0.0_wp; fy = 0.0_wp; fz = 0.0_wp

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


           bijmp = jmp1*b1
           fx  = bijmp*xx
           fy  = bijmp*yy
           fz  = bijmp*zz

           If (mxompl >= 1) Then

! field due to dipole on atom j

              dpxx= xx2*b2-b1
              dpyy= yy2*b2-b1
              dpzz= zz2*b2-b1
              dpxy= xx*yy*b2
              dpxz= xx*zz*b2
              dpyz= yy*zz*b2

              fx = fx + (dpxx*jmp2 + dpxy*jmp3 + dpxz*jmp4)
              fy = fy + (dpxy*jmp2 + dpyy*jmp3 + dpyz*jmp4)
              fz = fz + (dpxz*jmp2 + dpyz*jmp3 + dpzz*jmp4)

           End If

           If (mxompl == 2) Then

! field due to quadrupole on atom j

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

              fx = fx - (jmp5*dpxxx+jmp6*dpxxy+jmp7*dpxxz+jmp8*dpxyy+jmp9*dpxyz+jmp10*dpxzz)
              fy = fy - (jmp5*dpxxy+jmp6*dpxyy+jmp7*dpxyz+jmp8*dpyyy+jmp9*dpyyz+jmp10*dpyzz)
              fz = fz - (jmp5*dpxxz+jmp6*dpxyz+jmp7*dpxzz+jmp8*dpyyz+jmp9*dpyzz+jmp10*dpzzz)

           End If

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

! electrostatic field at position j due to atom i

           If (jatm <= natms) Then

              tjx = 0.0_wp; tjy=0.0_wp; tjz = 0.0_wp

! electrostatic field due to charge on atom i

              bijmp = imp1*b1
              tjx  = -bijmp*xx
              tjy  = -bijmp*yy
              tjz  = -bijmp*zz

! electrostatic field due to dipole on atom i

              If (mxompl >= 1) Then

                 tjx = tjx + (dpxx*imp2 + dpxy*imp3 + dpxz*imp4)
                 tjy = tjy + (dpxy*imp2 + dpyy*imp3 + dpyz*imp4)
                 tjz = tjz + (dpxz*imp2 + dpyz*imp3 + dpzz*imp4)

              End If

! electrostatic field due to quadrupole on atom i

              If (mxompl == 2) Then

                 tjx = tjx + (imp5*dpxxx+imp6*dpxxy+imp7*dpxxz+imp8*dpxyy+imp9*dpxyz+imp10*dpxzz)
                 tjy = tjy + (imp5*dpxxy+imp6*dpxyy+imp7*dpxyz+imp8*dpyyy+imp9*dpyyz+imp10*dpyzz)
                 tjz = tjz + (imp5*dpxxz+imp6*dpxyz+imp7*dpxzz+imp8*dpyyz+imp9*dpyzz+imp10*dpzzz)

              End If

              mpfldx(jatm) = mpfldx(jatm) + tjx
              mpfldy(jatm) = mpfldy(jatm) + tjy
              mpfldz(jatm) = mpfldz(jatm) + tjz

           End If

        End If

     End Do

! load back electrostatic field

     mpfldx(iatm)=fix
     mpfldy(iatm)=fiy
     mpfldz(iatm)=fiz

  End If

End Subroutine ewald_real_mfield_d
