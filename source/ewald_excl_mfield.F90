Subroutine ewald_excl_mfield &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating electrostatic field
! in a periodic system using multipoles with the ewald real space
! kernel
!
! Note: exclusion correction term
!
! copyright - daresbury laboratory
! author    - h.a.boateng november 2014
! amended   - i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup
  Use configuration, Only : natms,nlast,ltg,ltype,list,fxx,fyy,fzz
  Use mpoles_module, Only : mplmap,mplgfr,mpfldx,mpfldy,mpfldz
  Use ewald_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp
  Real( Kind = wp ), Parameter :: rr3  = 1.0_wp/3.0_wp
  Real( Kind = wp ), Parameter :: r10  = 0.1_wp
  Real( Kind = wp ), Parameter :: r42  = 1.0_wp/42.0_wp
  Real( Kind = wp ), Parameter :: r216 = 1.0_wp/216.0_wp

  Integer           :: limit,idi,ai,jatm,aj,k1,k2,k3,s1,s2,s3,m,n
  Integer           :: k,s11,s21,s31,ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

  Real( Kind = wp ) :: scl,rrr,fix,fiy,fiz,fx,fy,fz,           &
                       strs1,strs2,strs3,strs5,strs6,strs9,    &
                       alpr,alpr2,txyz,erfr,exp1,kx,ky,kz,     &
                       tyz,tt,ttxyz,tmp,tmpi,tmpj,tix,tiy,tiz, &
                       tjx,tjy,tjz,sx,sy,sz,ti,tj,t1,t2,t3


  Real( Kind = wp ) :: d1(-2:mxompl+1,-2:mxompl+1,-2:mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)

! global identity of iatm

  idi=ltg(iatm)

! get the multipoles for site i

  imp=mplgfr(:,iatm)

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) > zero_plus) Then

! multipole scaler

     scl=2.0_wp*alpha*r4pie0/(sqrpi*epsq)

! scale imp multipoles

     imp=imp*scl

! load field

     fix=mpfldx(iatm)
     fiy=mpfldy(iatm)
     fiz=mpfldz(iatm)

! initialize torques for atom i (temporary)

     tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

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

! scale jmp multipoles

           jmp=jmp*scl

! get the value of the kernel using 3pt interpolation

           alpr =rrr*alpha
           alpr2=alpr*alpr

! calculate error function and derivative

           If (alpr < 1.0e-2_wp) Then

! close particles (core-shell units) - small distances limit

              erfr=2.0_wp/sqrpi * &
              (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

           Else

! distant particles - traditional

              exp1=Exp(-(alpha*rrr)**2)
              tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

              erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

           End If

! compute derivatives of kernel

           Call ewald_deriv(-2,mxompl+1,2,erfr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

! calculate field

           fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
           tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

           sz = 1.0_wp
           Do s3=0,mxompl
              s31=s3+1

              sy = sz
              Do s2=0,mxompl-s3
                 s21=s2+1

                 sx = sy
                 Do s1=0,mxompl-s3-s2
                    s11=s1+1

                    n   = s1+s2+s3+1

                    ii  = mplmap(s1,s2,s3)

                    tmp = alpha**n

                    tj  = sx*tmp*jmp(ii)

                    t1  = d1(s11,s2,s3)
                    t2  = d1(s1,s21,s3)
                    t3  = d1(s1,s2,s31)

! field at position i due to atom j

                    If (Abs(jmp(ii)) > zero_plus) Then

                       fx = fx - tj*t1
                       fy = fy - tj*t2
                       fz = fz - tj*t3

                    End If

! field at position j due to atom i

                    If (Abs(imp(ii)) > zero_plus) Then

                       ti  = tmp*imp(ii)

                       tjx = tjx + ti*t1
                       tjy = tjy + ti*t2
                       tjz = tjz + ti*t3

                    End If

                    sx = -sx
                 End Do

                 sy = -sy
              End Do

              sz = -sz
           End Do

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              mpfldx(jatm)=mpfldx(jatm)+tjx
              mpfldy(jatm)=mpfldy(jatm)+tjy
              mpfldz(jatm)=mpfldz(jatm)+tjz

           End If

        End If

     End Do

! load back field

     mpfldx(iatm)=fix
     mpfldy(iatm)=fiy
     mpfldz(iatm)=fiz

  End If

End Subroutine ewald_excl_mfield
