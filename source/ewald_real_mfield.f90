Subroutine ewald_real_mfield &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating permanent electrostatic field
! in a periodic system using multipoles with the ewald real space
! kernel
!
! copyright - daresbury laboratory
! author    - h.a.boateng november 2014
! amended   - i.t.todorov september 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,gtime
  Use setup_module
  Use config_module, Only : natms,nlast,ltg,ltype,list
  Use mpoles_module, Only : mplmap,mplgfr,mplfldx,mplfldy,mplfldz

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xdf,ydf,zdf,rrt

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: drewd,rdrewd

  Integer           :: fail,idi,ai,jatm,aj,k,k1,k2,k3,s1,s2,s3,m,n
  Integer           :: s11,s21,s31,ii,jj

  Real( Kind = wp ) :: scl,rrr,fix,fiy,fiz,fx,fy,fz,         &
                       ppp,vk0,vk1,vk2,ti,tj,kx,ky,kz,tyz,   &
                       txyz,erfcr,ttxyz,tmp,tmpi,tmpj,tix,   &
                       tiy,tiz,tjx,tjy,tjz,sx,sy,sz,t1,t2,t3

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

  Real( Kind = wp ) :: d1(0:mxompl+1,0:mxompl+1,0:mxompl+1)
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

! reciprocal of interpolation interval

     rdrewd = 1.0_wp/drewd

! generate error function complement tables for ewald sum

     Call erfcgen(rcut,alpha,mxgele,erc,fer)
  End If

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

     fix=mplfldx(iatm)
     fiy=mplfldy(iatm)
     fiz=mplfldz(iatm)

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

           k   = Int(rrr*rdrewd)
           ppp = rrr*rdrewd - Real(k,wp)

           vk0 = erc(k)
           vk1 = erc(k+1)
           vk2 = erc(k+2)

           t1 = vk0 + (vk1 - vk0)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/alpha

! compute derivatives of kernel

           Call ewald_deriv( 0,mxompl+1,1,erfcr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

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

!  field at position i due to atom j

                    If (Abs(jmp(ii)) > zero_plus) Then

                       fx = fx - tj*t1
                       fy = fy - tj*t2
                       fz = fz - tj*t3

                    End If

!  field at position j due to atom i

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

              mplfldx(jatm)=mplfldx(jatm)+tjx
              mplfldy(jatm)=mplfldy(jatm)+tjy
              mplfldz(jatm)=mplfldz(jatm)+tjz

           End If

        End If

     End Do

! load back field

     mplfldx(iatm)=fix
     mplfldy(iatm)=fiy
     mplfldz(iatm)=fiz

  End If

End Subroutine ewald_real_mfield
