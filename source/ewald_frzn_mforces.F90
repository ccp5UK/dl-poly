Subroutine ewald_frzn_mforces(rcut,alpha,epsq,engcpe_fr,vircpe_fr,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating corrections to coulombic forces
! in a periodic system arising from multipoles on frozen pairs
!
! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
!       end (and any COM drift removed) but corrections to the stress
!       and the virial are important as they feed into the system
!       pressure response.  Constant volume ensembles (keyens < 20)
!       need this calculation just once! - controlled by lf_fce in
!       ewald_check<-two_body_forces
!
! copyright - daresbury laboratory
! author    - i.t.todorov & h.a.boateng february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module
  Use setup_module
  Use configuration, Only : cell,natms,list,ltg,lfrzn,xxx,yyy,zzz,fxx,fyy,fzz
  Use mpoles_module
  Use ewald_module

  Implicit None

  Real( Kind = wp ),                     Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ),                     Intent(   Out ) :: engcpe_fr,vircpe_fr
  Real( Kind = wp ), Dimension( 1:9 ),   Intent( InOut ) :: stress

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

  Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
  Integer           :: k1,k2,k3,s1,s2,s3,n
  Integer           :: ks1,ks2,ks3,ks11,ks21,ks31,nn,mm

  Real( Kind = wp ) :: scl,det,rcell(1:9),xrr,yrr,zrr,rrr,      &
                       engmpl,erfr,exp1,tt,t1,kx,ky,kz,         &
                       txyz,fx,fy,fz,xss,yss,zss,               &
                       strs1,strs2,strs3,strs5,strs6,strs9,tmp, &
                       alphan,tmpi,tmpj,tix,tiy,tiz,tjx,tjy,tjz,sx,sy,sz

  Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr
  Real( Kind = wp ), Dimension( : ), Allocatable :: xfr,yfr,zfr
  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

  Real( Kind = wp ), Dimension( : ), Allocatable :: mmp(:,:)
  Real( Kind = wp ), Dimension( : ), Allocatable :: mmpx(:,:)
  Real( Kind = wp ), Dimension( : ), Allocatable :: mmpy(:,:)
  Real( Kind = wp ), Dimension( : ), Allocatable :: mmpz(:,:)

  Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
  Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
  Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  If (.not.lf_fce) Then ! All's been done but needs copying
     Do i=1,natms
        fxx(i)=fxx(i)+ffx(i)
        fyy(i)=fyy(i)+ffy(i)
        fzz(i)=fzz(i)+ffz(i)
     End Do

     engcpe_fr=ef_fr
     vircpe_fr=vf_fr
     stress=stress+sf_fr

     If (l_cp) Then
        Do i=1,natms
           fcx(i)=fcx(i)+ffx(i)
           fcy(i)=fcy(i)+ffy(i)
           fcz(i)=fcz(i)+ffz(i)
        End Do

        e_fr=ef_fr
        v_fr=vf_fr
        s_fr=sf_fr
     End If

     Return
  End If

  fail=0
  Allocate (l_ind(1:mxatdm),nz_fr(0:mxnode), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mforces allocation failure, node: ', idnode
     Call error(0)
  End If

  Call invert(cell,rcell,det)

! Initialise contributions

  engcpe_fr=0.0_wp
  vircpe_fr=0.0_wp

  strs1 = 0.0_wp
  strs2 = 0.0_wp
  strs3 = 0.0_wp
  strs5 = 0.0_wp
  strs6 = 0.0_wp
  strs9 = 0.0_wp

! initialize torques for atom i (temporary)

  tix = 0.0_wp; tiy = 0.0_wp; tiz = 0.0_wp

  l_ind=0 ; nz_fr=0
  Do i=1,natms

! If using multipoles then get multipoles in global frame for atom i

     imp=mplgfr(:,i)

     If (lfrzn(i) > 0 .and. Maxval(Abs(imp)) > zero_plus) Then
        nz_fr(idnode+1)=nz_fr(idnode+1)+1
        l_ind(nz_fr(idnode+1))=i
     End If
  End Do
  If (mxnode > 1) Call gsum(nz_fr)
  nz_fr(0) = Sum(nz_fr(0:idnode)) ! Offset

  scl=2.0_wp*alpha*r4pie0/(sqrpi*epsq)
  nzfr = Sum(nz_fr(1:mxnode))     ! Total
  If (nzfr <= 10*mxatms) Then

     Allocate (mmp(1:mximpl,1:nzfr),                                              &
               mmpx(1:mximpl,1:nzfr),mmpy(1:mximpl,1:nzfr),mmpz(1:mximpl,1:nzfr), &
               xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mforces allocation failure 1, node: ', idnode
        Call error(0)
     End If

     mmp=0.0_wp
     mmpx=0.0_wp ; mmpz=0.0_wp ; mmpz=0.0_wp
     xfr=0.0_wp ; yfr=0.0_wp; zfr=0.0_wp
     Do i=1,nz_fr(idnode+1)
        ii=nz_fr(0)+i

! if using multipoles then get them in global frame for atom ii

        mmp(:,ii)=mplgfr(:,l_ind(i))

        If (mxompl > 0 .and. induce) Then

           mmp(2,ii)=mmp(2,ii)+indipx(l_ind(i))
           mmp(3,ii)=mmp(3,ii)+indipy(l_ind(i))
           mmp(4,ii)=mmp(4,ii)+indipz(l_ind(i))

        End If

! get the components for site ii infinitesimal rotations

        mmpx(:,ii)=mprotx(:,l_ind(i))
        mmpy(:,ii)=mproty(:,l_ind(i))
        mmpz(:,ii)=mprotz(:,l_ind(i))

        xfr(ii)=xxx(l_ind(i))
        yfr(ii)=yyy(l_ind(i))
        zfr(ii)=zzz(l_ind(i))
     End Do
     If (mxnode > 1) Then
        Call gsum(mmp)
        Call gsum(mmpx)
        Call gsum(mmpy)
        Call gsum(mmpz)

        Call gsum(xfr)
        Call gsum(yfr)
        Call gsum(zfr)
     End If

     Do i=1,nz_fr(idnode+1)
        ii=nz_fr(0)+i

! get the multipoles for site ii

        imp=mmp(:,ii)*scl

! get the components for site ii infinitesimal rotations

        impx=mmpx(:,ii)
        impy=mmpy(:,ii)
        impz=mmpz(:,ii)

        Do jj=1,nz_fr(0) ! -, on nodes<idnode
           xrr=xfr(ii)-xfr(jj)
           yrr=yfr(ii)-yfr(jj)
           zrr=zfr(ii)-zfr(jj)

           xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
           yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
           zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

           xss=xss-Anint(xss)
           yss=yss-Anint(yss)
           zss=zss-Anint(zss)

           xrr=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
           yrr=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
           zrr=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

! calculate interatomic distance

           rrr=Sqrt(xrr**2+yrr**2+zrr**2)

! get the multipoles for site jj

           jmp=mmp(:,jj)

! get the components for site jj infinitesimal rotations

           jmpx=mmpx(:,jj)
           jmpy=mmpy(:,jj)
           jmpz=mmpz(:,jj)

! calculate error function and derivative

           exp1=Exp(-(alpha*rrr)**2)
           tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

           Call ewald_deriv(-2,2*mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                       nn = mplmap(k1,k2,k3)

                       If (Abs(jmp(nn)) > zero_plus) Call explicit_ewald_real_loops &
           (-2,2*mxompl+1, k1,k2,k3, alpha, d1,               &
           imp,       impx,    impy,    impz,    tix,tiy,tiz, &
           kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
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

                       nn = mplmap(k1,k2,k3)

                       If (Abs(jmp(nn)) > zero_plus) Then

                          txyz=kx*jmp(nn)

                          sz = 1.0_wp
                          Do s3=0,mxompl
                             ks3=k3+s3; ks31=ks3+1

                             sy = sz
                             Do s2=0,mxompl-s3
                             ks2=k2+s2; ks21=ks2+1

                                sx = sy
                                Do s1=0,mxompl-s3-s2
                                   ks1=k1+s1; ks11=ks1+1

                                   n       = ks1+ks2+ks3
                                   alphan  = alpha**n

                                   mm      = mplmap(s1,s2,s3)

                                   tmp     = alphan*d1(ks1,ks2,ks3)

                                   tmpi    = txyz       * tmp
                                   tmpj    = sx*imp(mm) * tmp

                                   t1      = alphan     * txyz*imp(mm)

! energy
                                   engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                                   t1      = t1*alpha

! force

                                   fx      = fx      - t1*d1(ks11,ks2,ks3)
                                   fy      = fy      - t1*d1(ks1,ks21,ks3)
                                   fz      = fz      - t1*d1(ks1,ks2,ks31)

! torque on iatm

                                   tix     = tix     + impx(mm)*tmpi
                                   tiy     = tiy     + impy(mm)*tmpi
                                   tiz     = tiz     + impz(mm)*tmpi

! torque on jatm

                                   tjx     = tjx     + jmpx(nn)*tmpj
                                   tjy     = tjy     + jmpy(nn)*tmpj
                                   tjz     = tjz     + jmpz(nn)*tmpj

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

! calculate forces

           fxx(l_ind(i))=fxx(l_ind(i))-fx
           fyy(l_ind(i))=fyy(l_ind(i))-fy
           fzz(l_ind(i))=fzz(l_ind(i))-fz

! redundant calculations copying

           If (lf_cp) Then
              ffx(l_ind(i))=ffx(l_ind(i))-fx
              ffy(l_ind(i))=ffy(l_ind(i))-fy
              ffz(l_ind(i))=ffz(l_ind(i))-fz
           End If

! infrequent calculations copying

           If (l_cp) Then
              fcx(l_ind(i))=fcx(l_ind(i))-fx
              fcy(l_ind(i))=fcy(l_ind(i))-fy
              fcz(l_ind(i))=fcz(l_ind(i))-fz
           End If
        End Do

        Do j=i+1,nz_fr(idnode+1) ! =, node=idnode (OVERLAP but no SELF)!
           jj=nz_fr(0)+j

           xrr=xfr(ii)-xfr(jj)
           yrr=yfr(ii)-yfr(jj)
           zrr=zfr(ii)-zfr(jj)

           xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
           yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
           zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

           xss=xss-Anint(xss)
           yss=yss-Anint(yss)
           zss=zss-Anint(zss)

           xrr=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
           yrr=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
           zrr=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

! interatomic distance

           rrr=Sqrt(xrr**2+yrr**2+zrr**2)

! get the multipoles for site jj

           jmp=mmp(:,jj)

! get the components for site jj infinitesimal rotations

           jmpx=mmpx(:,jj)
           jmpy=mmpy(:,jj)
           jmpz=mmpz(:,jj)

! calculate error function and derivative

           exp1=Exp(-(alpha*rrr)**2)
           tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

           Call ewald_deriv(-2,2*mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                       nn = mplmap(k1,k2,k3)

                       If (Abs(jmp(nn)) > zero_plus) Call explicit_ewald_real_loops &
           (-2,2*mxompl+1, k1,k2,k3, alpha, d1,               &
           imp,       impx,    impy,    impz,    tix,tiy,tiz, &
           kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
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

                       nn = mplmap(k1,k2,k3)

                       If (Abs(jmp(nn)) > zero_plus) Then

                          txyz=kx*jmp(nn)

                          sz = 1.0_wp
                          Do s3=0,mxompl
                             ks3=k3+s3; ks31=ks3+1

                             sy = sz
                             Do s2=0,mxompl-s3
                             ks2=k2+s2; ks21=ks2+1

                                sx = sy
                                Do s1=0,mxompl-s3-s2
                                   ks1=k1+s1; ks11=ks1+1

                                   n       = ks1+ks2+ks3
                                   alphan  = alpha**n

                                   mm      = mplmap(s1,s2,s3)

                                   tmp     = alphan*d1(ks1,ks2,ks3)

                                   tmpi    = txyz       * tmp
                                   tmpj    = sx*imp(mm) * tmp

                                   t1      = alphan     * txyz*imp(mm)

! energy
                                   engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                                   t1      = t1*alpha

! force

                                   fx      = fx      - t1*d1(ks11,ks2,ks3)
                                   fy      = fy      - t1*d1(ks1,ks21,ks3)
                                   fz      = fz      - t1*d1(ks1,ks2,ks31)

! torque on iatm

                                   tix     = tix     + impx(mm)*tmpi
                                   tiy     = tiy     + impy(mm)*tmpi
                                   tiz     = tiz     + impz(mm)*tmpi

! torque on jatm

                                   tjx     = tjx     + jmpx(nn)*tmpj
                                   tjy     = tjy     + jmpy(nn)*tmpj
                                   tjz     = tjz     + jmpz(nn)*tmpj

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

! calculate forces

           fxx(l_ind(i))=fxx(l_ind(i))-fx
           fyy(l_ind(i))=fyy(l_ind(i))-fy
           fzz(l_ind(i))=fzz(l_ind(i))-fz

           fxx(l_ind(j))=fxx(l_ind(j))+fx
           fyy(l_ind(j))=fyy(l_ind(j))+fy
           fzz(l_ind(j))=fzz(l_ind(j))+fz

           mptrqx(l_ind(j))=mptrqx(l_ind(j))+tjx
           mptrqy(l_ind(j))=mptrqy(l_ind(j))+tjy
           mptrqz(l_ind(j))=mptrqz(l_ind(j))+tjz

! redundant calculations copying

           If (lf_cp) Then
              ffx(l_ind(i))=ffx(l_ind(i))-fx
              ffy(l_ind(i))=ffy(l_ind(i))-fy
              ffz(l_ind(i))=ffz(l_ind(i))-fz

              ffx(l_ind(j))=ffx(l_ind(j))+fx
              ffy(l_ind(j))=ffy(l_ind(j))+fy
              ffz(l_ind(j))=ffz(l_ind(j))+fz
           End If

! infrequent calculations copying

           If (l_cp) Then
              fcx(l_ind(i))=fcx(l_ind(i))-fx
              fcy(l_ind(i))=fcy(l_ind(i))-fy
              fcz(l_ind(i))=fcz(l_ind(i))-fz

              fcx(l_ind(j))=fcx(l_ind(j))+fx
              fcy(l_ind(j))=fcy(l_ind(j))+fy
              fcz(l_ind(j))=fcz(l_ind(j))+fz
           End If

! calculate potential energy and virial

           engcpe_fr = engcpe_fr - engmpl
           vircpe_fr = vircpe_fr - (fx*xrr + fy*yrr + fz*zrr)

! calculate stress tensor

           strs1 = strs1 + xrr*fx
           strs2 = strs2 + xrr*fy
           strs3 = strs3 + xrr*fz
           strs5 = strs5 + yrr*fy
           strs6 = strs6 + yrr*fz
           strs9 = strs9 + zrr*fz
        End Do

        Do jj=nz_fr(0)+nz_fr(idnode+1)+1,nzfr ! +, on nodes>idnode
           xrr=xfr(ii)-xfr(jj)
           yrr=yfr(ii)-yfr(jj)
           zrr=zfr(ii)-zfr(jj)

           xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
           yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
           zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

           xss=xss-Anint(xss)
           yss=yss-Anint(yss)
           zss=zss-Anint(zss)

           xrr=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
           yrr=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
           zrr=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

! interatomic distance

           rrr=Sqrt(xrr**2+yrr**2+zrr**2)

! get the multipoles for site jj

           jmp=mmp(:,jj)

! get the components for site jj infinitesimal rotations

           jmpx=mmpx(:,jj)
           jmpy=mmpy(:,jj)
           jmpz=mmpz(:,jj)

! calculate error function and derivative

           exp1=Exp(-(alpha*rrr)**2)
           tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

           Call ewald_deriv(-2,2*mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                       nn = mplmap(k1,k2,k3)

                       If (Abs(jmp(nn)) > zero_plus) Call explicit_ewald_real_loops &
           (-2,2*mxompl+1, k1,k2,k3, alpha, d1,               &
           imp,       impx,    impy,    impz,    tix,tiy,tiz, &
           kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
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

                       nn = mplmap(k1,k2,k3)

                       If (Abs(jmp(nn)) > zero_plus) Then

                          txyz=kx*jmp(nn)

                          sz = 1.0_wp
                          Do s3=0,mxompl
                             ks3=k3+s3; ks31=ks3+1

                             sy = sz
                             Do s2=0,mxompl-s3
                             ks2=k2+s2; ks21=ks2+1

                                sx = sy
                                Do s1=0,mxompl-s3-s2
                                   ks1=k1+s1; ks11=ks1+1

                                   n       = ks1+ks2+ks3
                                   alphan  = alpha**n

                                   mm      = mplmap(s1,s2,s3)

                                   tmp     = alphan*d1(ks1,ks2,ks3)

                                   tmpi    = txyz       * tmp
                                   tmpj    = sx*imp(mm) * tmp

                                   t1      = alphan     * txyz*imp(mm)

! energy
                                   engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                                   t1      = t1*alpha

! force

                                   fx      = fx      - t1*d1(ks11,ks2,ks3)
                                   fy      = fy      - t1*d1(ks1,ks21,ks3)
                                   fz      = fz      - t1*d1(ks1,ks2,ks31)

! torque on iatm

                                   tix     = tix     + impx(mm)*tmpi
                                   tiy     = tiy     + impy(mm)*tmpi
                                   tiz     = tiz     + impz(mm)*tmpi

! torque on jatm

                                   tjx     = tjx     + jmpx(nn)*tmpj
                                   tjy     = tjy     + jmpy(nn)*tmpj
                                   tjz     = tjz     + jmpz(nn)*tmpj

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

! calculate forces

           fxx(l_ind(i))=fxx(l_ind(i))-fx
           fyy(l_ind(i))=fyy(l_ind(i))-fy
           fzz(l_ind(i))=fzz(l_ind(i))-fz

! redundant calculations copying

           If (lf_cp) Then
              ffx(l_ind(i))=ffx(l_ind(i))-fx
              ffy(l_ind(i))=ffy(l_ind(i))-fy
              ffz(l_ind(i))=ffz(l_ind(i))-fz
           End If

! infrequent calculations copying

           If (l_cp) Then
              fcx(l_ind(i))=fcx(l_ind(i))-fx
              fcy(l_ind(i))=fcy(l_ind(i))-fy
              fcz(l_ind(i))=fcz(l_ind(i))-fz
           End If

! calculate potential energy and virial

           engcpe_fr = engcpe_fr - engmpl
           vircpe_fr = vircpe_fr - (fx*xrr + fy*yrr + fz*zrr)

! calculate stress tensor

           strs1 = strs1 + xrr*fx
           strs2 = strs2 + xrr*fy
           strs3 = strs3 + xrr*fz
           strs5 = strs5 + yrr*fy
           strs6 = strs6 + yrr*fz
           strs9 = strs9 + zrr*fz
        End Do
     End Do

     Deallocate (mmp,mmpx,mmpy,mmpz,xfr,yfr,zfr, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mforces deallocation failure 1, node: ', idnode
        Call error(0)
     End If

  Else

! We resort to approximating N*(N-1)/2 interactions
! with the short-range one from the two body linked cell list

     Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mforces allocation failure 2, node: ', idnode
        Call error(0)
     End If

     Do ii=1,nz_fr(idnode+1)
        i=l_ind(nz_fr(idnode+1))
        idi=ltg(ii)

! get the multipoles for site i

        imp=mplgfr(:,i)*scl

! get the components for site i infinitesimal rotations

        impx=mprotx(:,i)
        impy=mproty(:,i)
        impz=mprotz(:,i)

! Get list limit

        limit=list(-2,i)-list(-1,i)
        If (limit > 0) Then

! calculate interatomic distances

           Do k=1,limit
              j=list(list(-1,i)+k,i)

              xxt(k)=xxx(i)-xxx(j)
              yyt(k)=yyy(i)-yyy(j)
              zzt(k)=zzz(i)-zzz(j)
           End Do

! periodic boundary conditions
!
!           Call images(imcon,cell,limit,xxt,yyt,zzt)

! get distances

           Do k=1,limit
              rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
           End Do

           Do k=1,limit
              j=list(list(-1,i)+k,i)

! get the multipoles for site j

              jmp=mplgfr(:,j)

! get the components for site j infinitesimal rotations

              jmpx=mprotx(:,j)
              jmpy=mproty(:,j)
              jmpz=mprotz(:,j)

! interatomic distance

              rrr=rrt(k)

! truncation of potential

              If (Maxval(Abs(jmp)) > zero_plus .and. rrr < rcut) Then

! calculate error function and derivative

                 exp1=Exp(-(alpha*rrr)**2)
                 tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

                 erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

                 Call ewald_deriv(-2,2*mxompl+1,2,erfr,alpha*xxt(k),alpha*yyt(k),alpha*zzt(k),alpha*rrr,d1)

! calculate forces

                 engmpl = 0.0_wp
                 fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                 tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                 If (mxompl < 5) Then ! fully rewritten by ITT

                    kz = 1.0_wp
                    Do k3=0,mxompl

                       ky = kz
                       Do k2=0,mxompl-k3

                          kx = ky
                          Do k1=0,mxompl-k3-k2

                             nn = mplmap(k1,k2,k3)

                             If (Abs(jmp(nn)) > zero_plus) Call explicit_ewald_real_loops &
           (-2,2*mxompl+1, k1,k2,k3, alpha, d1,               &
           imp,       impx,    impy,    impz,    tix,tiy,tiz, &
           kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
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

                             nn = mplmap(k1,k2,k3)

                             If (Abs(jmp(nn)) > zero_plus) Then

                                txyz=kx*jmp(nn)

                                sz = 1.0_wp
                                Do s3=0,mxompl
                                   ks3=k3+s3; ks31=ks3+1

                                   sy = sz
                                   Do s2=0,mxompl-s3
                                   ks2=k2+s2; ks21=ks2+1

                                      sx = sy
                                      Do s1=0,mxompl-s3-s2
                                         ks1=k1+s1; ks11=ks1+1

                                         n      = ks1+ks2+ks3
                                         alphan = alpha**n

                                         mm     = mplmap(s1,s2,s3)

                                         tmp    = alphan*d1(ks1,ks2,ks3)

                                         tmpi   = txyz       * tmp
                                         tmpj   = sx*imp(mm) * tmp

                                         t1     = alphan     * txyz*imp(mm)

! energy

                                         engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

! force

                                         t1      = t1*alpha

                                         fx      = fx      - t1*d1(ks11,ks2,ks3)
                                         fy      = fy      - t1*d1(ks1,ks21,ks3)
                                         fz      = fz      - t1*d1(ks1,ks2,ks31)

! torque on iatm

                                         tix     = tix     + impx(mm)*tmpi
                                         tiy     = tiy     + impy(mm)*tmpi
                                         tiz     = tiz     + impz(mm)*tmpi

! torque on jatm

                                         tjx     = tjx     + jmpx(nn)*tmpj
                                         tjy     = tjy     + jmpy(nn)*tmpj
                                         tjz     = tjz     + jmpz(nn)*tmpj

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

! calculate forces

                 fxx(i)=fxx(i)-fx
                 fyy(i)=fyy(i)-fy
                 fzz(i)=fzz(i)-fz

! redundant calculations copying

                 If (lf_cp) Then
                    ffx(i)=ffx(i)-fx
                    ffy(i)=ffy(i)-fy
                    ffz(i)=ffz(i)-fz
                 End If

! infrequent calculations copying

                 If (l_cp) Then
                    fcx(i)=fcx(i)-fx
                    fcy(i)=fcy(i)-fy
                    fcz(i)=fcz(i)-fz
                 End If

                 If (j <= natms) Then

                    fxx(j)=fxx(j)+fx
                    fyy(j)=fyy(j)+fy
                    fzz(j)=fzz(j)+fz

                    mptrqx(j)=mptrqx(j)+tjx
                    mptrqy(j)=mptrqy(j)+tjy
                    mptrqz(j)=mptrqz(j)+tjz

! redundant calculations copying

                    If (lf_cp) Then
                       ffx(j)=ffx(j)+fx
                       ffy(j)=ffy(j)+fy
                       ffz(j)=ffz(j)+fz
                    End If

! infrequent calculations copying

                    If (l_cp) Then
                       fcx(j)=fcx(j)+fx
                       fcy(j)=fcy(j)+fy
                       fcz(j)=fcz(j)+fz
                    End If

                 End If

                 If (j <= natms .or. idi < ltg(j)) Then

! calculate potential energy and virial

                    engcpe_fr = engcpe_fr - engmpl
                    vircpe_fr = vircpe_fr - (fx*xxt(k) + fy*yyt(k) + fz*zzt(k))

! calculate stress tensor

                    strs1 = strs1 + xxt(k)*fx
                    strs2 = strs2 + xxt(k)*fy
                    strs3 = strs3 + xxt(k)*fz
                    strs5 = strs5 + yyt(k)*fy
                    strs6 = strs6 + yyt(k)*fz
                    strs9 = strs9 + zzt(k)*fz

                 End If

              End If
           End Do

        End If

!  update torques due to multipoles

        mptrqx(i)=mptrqx(i)+scl*tix
        mptrqy(i)=mptrqy(i)+scl*tiy
        mptrqz(i)=mptrqz(i)+scl*tiz

     End Do

     Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mforces deallocation failure 2, node: ', idnode
        Call error(0)
     End If

  End If

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

! redundant calculations copying

  If (lf_cp) Then
     ef_fr=engcpe_fr
     vf_fr=vircpe_fr

     sf_fr(1) = strs1
     sf_fr(2) = strs2
     sf_fr(3) = strs3
     sf_fr(4) = strs2
     sf_fr(5) = strs5
     sf_fr(6) = strs6
     sf_fr(7) = strs3
     sf_fr(8) = strs6
     sf_fr(9) = strs9
  End If

! infrequent calculations copying

  If (l_cp) Then
     e_fr=engcpe_fr
     v_fr=vircpe_fr

     s_fr(1) = strs1
     s_fr(2) = strs2
     s_fr(3) = strs3
     s_fr(4) = strs2
     s_fr(5) = strs5
     s_fr(6) = strs6
     s_fr(7) = strs3
     s_fr(8) = strs6
     s_fr(9) = strs9
  End If

  Deallocate (l_ind,nz_fr, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mforces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine ewald_frzn_mforces
