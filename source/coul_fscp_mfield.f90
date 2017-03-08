Subroutine coul_fscp_mfield &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the permanent electrostatic field
! for a system of multipoles assuming a force shifted coulomb potential
! kernel
!
! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
!
! copyright - daresbury laboratory
! author    - h.a.boateng december 2014
! amended   - i.t.todorov december 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : natms,nlast,ltg,ltype,list
  Use mpoles_module, Only : mplmap,mplgfr,mpfldx,mpfldy,mpfldz

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt

  Logical,           Save :: newjob = .true. , damp
  Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                             rdrewd = 0.0_wp , &
                             aa     = 0.0_wp , &
                             bb     = 0.0_wp

  Integer           :: fail,k,n,idi,ai,jatm,aj,s1,s2,s3,m
  Integer           :: s11,s21,s31,ii,jj

  Real( Kind = wp ) :: scl,rrr,fix,fiy,fiz,fx,fy,fz,       &
                       ppp,vk0,vk1,vk2,t1,t2,kx,ky,kz,tyz, &
                       erfcr,tix,tiy,tiz,tjx,tjy,          &
                       tjz,tmp,tmpi,tmpj,ti,tj,t3,sx,sy,sz

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

  Real( Kind = wp ) :: d1(-2:mxompl+1,-2:mxompl+1,-2:mxompl+1)
  Real( Kind = wp ) :: a1(-2:mxompl+1,-2:mxompl+1,-2:mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)

  If (newjob) Then
     newjob = .false.

     If (alpha > zero_plus) Then
        damp = .true.
     Else
        damp = .false.
     End If

     If (damp) Then

! interpolation interval

        drewd = rcut/Real(mxgele-4,wp)

! reciprocal of interpolation interval

        rdrewd = 1.0_wp/drewd

        fail=0
        Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
        If (fail > 0) Then
           Write(nrite,'(/,1x,a,i0)') 'coul_fscp_mfield allocation failure, idnode: ', idnode
           Call error(0)
        End If

! generate error function complement tables for ewald sum

        Call erfcgen(rcut,alpha,mxgele,erc,fer)

! set force and potential shifting parameters (screened terms)

        aa =   fer(mxgele-4)*rcut
        bb = -(erc(mxgele-4)+aa*rcut)

     Else

! set force and potential shifting parameters (screened terms)

        aa =  1.0_wp/rcut**2
        bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)

     End If
  End If

! global identity of iatm

  idi=ltg(iatm)

! get the multipoles for site i

  imp=mplgfr(:,iatm)

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) > zero_plus) Then

! multipole scaler

     scl=r4pie0/epsq

! scale imp multipoles

     imp=imp*scl

! load field

     fix=mpfldx(iatm)
     fiy=mpfldy(iatm)
     fiz=mpfldz(iatm)

! start of primary loop for forces evaluation

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

           If (damp) Then

! get the value of the  real space kernel using 3pt interpolation

              k   = Int(rrr*rdrewd)
              ppp = rrr*rdrewd - Real(k,wp)

              vk0 = erc(k)
              vk1 = erc(k+1)
              vk2 = erc(k+2)

              t1 = vk0 + (vk1 - vk0)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/alpha

! compute derivatives of the ewald real space kernel

              Call ewald_deriv(-2,mxompl+1,1,erfcr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

! scale the derivatives into the right form

              d1 = 2.0_wp*alpha*d1/sqrpi

! compute derivatives of 'r'

              Call coul_deriv(-1,mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

! scale the derivatives of 'r'

              a1 = aa*a1

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

                       tj  = sx*jmp(ii)

                       t1  = tmp*d1(s11,s2,s3)+a1(s11,s2,s3)
                       t2  = tmp*d1(s1,s21,s3)+a1(s1,s21,s3)
                       t3  = tmp*d1(s1,s2,s31)+a1(s1,s2,s31)

! field at position i due to atom j

                       If (Abs(jmp(ii)) > zero_plus) Then

                          fx = fx - tj*t1
                          fy = fy - tj*t2
                          fz = fz - tj*t3

                       End If

! field at position j due to atom i

                       If (Abs(imp(ii)) > zero_plus) Then

                          ti  = imp(ii) ! multiply by sx to account for change in
                                        ! sign for odd derivatives. sx*sx=1.0
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

           Else

! compute derivatives of '1/r'

              Call coul_deriv( 1,mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

! compute derivatives of 'r'

              Call coul_deriv(-1,mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

! scale the derivatives of 'r' and add to d1

              d1 = d1 + a1/rcut**2

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

                       ii = mplmap(s1,s2,s3)
                       tj = sx*jmp(ii)

                       t1 = d1(s11,s2,s3)
                       t2 = d1(s1,s21,s3)
                       t3 = d1(s1,s2,s31)

! field at position i due to atom j

                       If (Abs(jmp(ii)) > zero_plus) Then

                          fx = fx - tj*t1
                          fy = fy - tj*t2
                          fz = fz - tj*t3

                       End If

! field at position j due to atom i

                       If (Abs(imp(ii)) > zero_plus) Then

                          ti  = imp(ii) ! multiply by sx to account for change in
                                        ! sign for odd derivatives. sx*sx=1.0
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

        End If

     End Do

! load back field

     mpfldx(iatm)=fix
     mpfldy(iatm)=fiy
     mpfldz(iatm)=fiz

  End If

End Subroutine coul_fscp_mfield
