Subroutine coul_dddp_mfield &
           (iatm,rcut,epsq,xxt,yyt,zzt,rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the permanent electrostatic field
! for a system of multipoles with 1/r kernel assuming a distance
! dependend dielectric 'constant'
!
! copyright - daresbury laboratory
! author    - h.a.boateng december 2014
! amended   - i.t.todorov september 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,nlast,ltg,ltype,list
  Use mpoles_module, Only : mplmap,mplgfr,mplfldx,mplfldy,mplfldz

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt

  Integer           :: idi,ai,jatm,aj,s1,s2,s3,m
  Integer           :: s11,s21,s31,ii,jj

  Real( Kind = wp ) :: scl.fix,fiy,fiz,fx,fy,fz,          &
                       kx,ky,kz,tyz,txyz,tix,tiy,tiz,tjx, &
                       tjy,tjz,tmp,t1,t2,t3,tj,ti,sx,sy,sz

  Real( Kind = wp ) :: d1(-2:mxompl+1,-2:mxompl+1,-2:mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)

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

     fix=mplfldx(iatm)
     fiy=mplfldy(iatm)
     fiz=mplfldz(iatm)

! start of primary loop for forces evaluation

     Do m=1,list(0,iatm)

! atomic index

        jatm=list(m,iatm)

! get the multipoles for site j

        jmp=mplgfr(:,jatm)

! truncation of potential - rrt(m) is the interatomic distance

        If (Maxval(Abs(jmp)) > zero_plus .and. rrt(m) < rcut) Then

! scale jmp multipoles

           jmp=jmp*scl

! compute derivatives of kernel

           Call coul_deriv(2,mxompl+1,xxt(m),yyt(m),zzt(m),rrt(m),d1)

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
                                     ! sign for odd derivatives, i.e. sx*sx=1
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

End Subroutine coul_dddp_mfield
