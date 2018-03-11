Subroutine coul_dddp_mforces &
           (iatm,rcut,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using multipoles with 1/r kernel assuming a
! distance dependent dielectric 'constant'
!
! copyright - daresbury laboratory
! author    - h.a.boateng february 2016
! amended   - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module
  Use configuration, Only : natms,ltg,list,fxx,fyy,fzz
  Use mpoles_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Integer           :: idi,jatm,k1,k2,k3,s1,s2,s3,m, &
                       ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

  Real( Kind = wp ) :: scl,engmpl,fix,fiy,fiz,fx,fy,fz,     &
                       strs1,strs2,strs3,strs5,strs6,strs9, &
                       t1,kx,ky,kz,txyz,tix,tiy,tiz,tjx,    &
                       tjy,tjz,tmp,tmpi,tmpj,sx,sy,sz

  Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
  Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
  Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

! initialise potential energy and virial

  engcpe=0.0_wp
  vircpe=0.0_wp

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

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) > zero_plus) Then

! get the components for site i infinitesimal rotations

     impx=mprotx(:,iatm)
     impy=mproty(:,iatm)
     impz=mprotz(:,iatm)

! multipole scaler

     scl=r4pie0/epsq

! scale imp multipoles

     imp=imp*scl

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

! initialize torques for atom i (temporary)

     tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

! start of primary loop for forces evaluation

     Do m=1,list(0,iatm)

! atomic index

        jatm=list(m,iatm)

! get the multipoles for site j

        jmp=mplgfr(:,jatm)

        If (mxompl > 0 .and. induce) Then

           jmp(2)=jmp(2)+indipx(jatm)
           jmp(3)=jmp(3)+indipy(jatm)
           jmp(4)=jmp(4)+indipz(jatm)

        End If

! truncation of potential - rrt(m) is the interatomic distance

        If (Maxval(Abs(jmp)) > zero_plus .and. rrt(m) < rcut) Then

! get the components for site j infinitesimal rotations

           jmpx=mprotx(:,jatm)
           jmpy=mproty(:,jatm)
           jmpz=mprotz(:,jatm)

! compute derivatives of kernel

           Call coul_deriv(2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrt(m),d1)

! calculate forces

           engmpl = 0.0_wp
           fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
           tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

           If (mxompl < 5) Then

              kz = 1.0_wp
              Do k3=0,mxompl

                 ky = kz
                 Do k2=0,mxompl-k3

                    kx = ky
                    Do k1=0,mxompl-k3-k2

                       jj = mplmap(k1,k2,k3)

                       If (Abs(jmp(jj)) > zero_plus) Call explicit_ewald_real_loops &
           (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
           imp,       impx,    impy,    impz,    tix,tiy,tiz, &
           kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
           engmpl,fx,fy,fz)

                       kx = -kx

                    End Do

                    ky = -ky

                 End Do

                 kz = -kz

              End Do

           Else

              kz = 1.0_wp
              Do k3=0,mxompl

                 ky = kz
                 Do k2=0,mxompl-k3

                    kx = ky
                    Do k1=0,mxompl-k3-k2

                       jj = mplmap(k1,k2,k3)

                       If (Abs(jmp(jj)) > zero_plus) Then

                          txyz=kx*jmp(jj)

                          sz = 1.0_wp
                          Do s3=0,mxompl
                             ks3=k3+s3; ks31=ks3+1

                             sy = sz
                             Do s2=0,mxompl-s3
                                ks2=k2+s2; ks21=ks2+1

                                sx = sy
                                Do s1=0,mxompl-s3-s2
                                   ks1=k1+s1; ks11=ks1+1

                                   ii   = mplmap(s1,s2,s3)

                                   tmp  = d1(ks1,ks2,ks3)

                                   tmpi = txyz       * tmp
                                   tmpj = sx*imp(ii) * tmp

                                   t1   = txyz*imp(ii)

! energy

                                   engmpl = engmpl  + t1*tmp

! force

                                   fx     = fx     - t1*d1(ks11,ks2,ks3)
                                   fy     = fy     - t1*d1(ks1,ks21,ks3)
                                   fz     = fz     - t1*d1(ks1,ks2,ks31)

! torque on iatm

                                   tix    = tix    + impx(ii)*tmpi
                                   tiy    = tiy    + impy(ii)*tmpi
                                   tiz    = tiz    + impz(ii)*tmpi

! torque on jatm

                                   tjx    = tjx    + jmpx(jj)*tmpj
                                   tjy    = tjy    + jmpy(jj)*tmpj
                                   tjz    = tjz    + jmpz(jj)*tmpj

                                   sx = -sx
                                End Do

                                sy = -sy
                             End Do

                             sz = -sz
                          End Do

                       End If

                       kx = -kx

                    End Do

                    ky = -ky

                 End Do

                 kz = -kz

              End Do

           End If

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

              engcpe = engcpe + engmpl

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

! virial

     vircpe = -engcpe

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

End Subroutine coul_dddp_mforces
