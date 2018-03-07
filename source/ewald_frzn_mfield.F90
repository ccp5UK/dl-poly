Subroutine ewald_frzn_mfield(rcut,alpha,epsq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating electrostiac field corrections to
! coulombic forces in a periodic system arising from multipoles on
! frozen pairs
!
! copyright - daresbury laboratory
! author    - h.a.boateng december 2014
! amended   - i.t.todorov december 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module
  Use setup_module
  Use configuration, Only : cell,natms,nlast,list,ltg,lfrzn, &
                            xxx,yyy,zzz,fxx,fyy,fzz
  Use mpoles_module
  Use ewald_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rcut,alpha,epsq

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

  Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
  Integer           :: k1,k2,k3,s1,s2,s3,m,n
  Integer           :: ks1,ks2,ks3,s11,s21,s31,nn,mm

  Real( Kind = wp ) :: scl,det,rcell(1:9),xrr,yrr,zrr,rrr,  &
                       erfr,egamma,exp1,tt,kx,ky,kz,        &
                       tyz,txyz,ttxyz,fx,fy,fz,xss,yss,zss, &
                       strs1,strs2,strs3,strs5,strs6,strs9, &
                       alphan,tjx,tjy,tjz,sx,sy,sz,         &
                       ti,tj,t1,t2,t3

  Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr
  Real( Kind = wp ), Dimension( : ), Allocatable :: xfr,yfr,zfr
  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

  Real( Kind = wp ), Dimension( : ), Allocatable :: mmp(:,:)

  Real( Kind = wp ) :: d1(-2:mxompl+1,-2:mxompl+1,-2:mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)

  fail=0
  Allocate (l_ind(1:mxatdm),nz_fr(0:mxnode), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mfield allocation failure, node: ', idnode
     Call error(0)
  End If

  Call invert(cell,rcell,det)

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

     Allocate (mmp(1:mximpl,1:nzfr),xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mfield allocation failure 1, node: ', idnode
        Call error(0)
     End If

     imp=0.0_wp ; jmp=0.0_wp
     mmp=0.0_wp
     xfr=0.0_wp ; yfr=0.0_wp ; zfr=0.0_wp
     Do i=1,nz_fr(idnode+1)
        ii=nz_fr(0)+i

! if using multipoles then get them in global frame for atom ii

        mmp(:,ii)=mplgfr(:,l_ind(i))

        xfr(ii)=xxx(l_ind(i))
        yfr(ii)=yyy(l_ind(i))
        zfr(ii)=zzz(l_ind(i))
     End Do
     If (mxnode > 1) Then
        Call gsum(mmp)

        Call gsum(xfr)
        Call gsum(yfr)
        Call gsum(zfr)
     End If

     Do i=1,nz_fr(idnode+1)
        ii=nz_fr(0)+i

! get the multipoles for site ii

        imp=mmp(:,ii)*scl

        Do jj=1,nz_fr(0) ! -, on nodes < idnode
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

           jmp=mmp(:,jj)*scl

! calculate error function and derivative

           exp1=Exp(-(alpha*rrr)**2)
           tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

           Call ewald_deriv(-2,mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                    n      = s1+s2+s3+1 !+1 added ??? ITT

                    nn     = mplmap(s1,s2,s3)

                    alphan = alpha**n

                    tj     = sx*alphan*jmp(nn) ! sx added here for change of sign below to work ??? ITT

                    t1  = d1(s11,s2,s3)
                    t2  = d1(s1,s21,s3)
                    t3  = d1(s1,s2,s31)

! field at position ii due to atom jj

                    If (Abs(jmp(nn)) > zero_plus) Then

                       fx = fx - tj*t1
                       fy = fy - tj*t2
                       fz = fz - tj*t3

                    End If

! field at position jj due to atom ii

                    If (Abs(imp(nn)) > zero_plus) Then

                       ti  = sx*alphan*imp(nn) ! multiply by sx to account for change in
                                               ! sign for odd derivatives

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

! calculate field

           mpfldx(l_ind(i))=mpfldx(l_ind(i))-fx
           mpfldy(l_ind(i))=mpfldy(l_ind(i))-fy
           mpfldz(l_ind(i))=mpfldz(l_ind(i))-fz

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

! calculate interatomic distance

           rrr=Sqrt(xrr**2+yrr**2+zrr**2)

! get the multipoles for site jj

           jmp=mmp(:,jj)*scl

! calculate error function and derivative

           exp1=Exp(-(alpha*rrr)**2)
           tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

           Call ewald_deriv(-2,mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                    n      = s1+s2+s3+1 !+1 added ??? ITT

                    nn     = mplmap(s1,s2,s3)

                    alphan = alpha**n

                    tj     = sx*alphan*jmp(nn) ! sx added here for change of sign below to work ??? ITT

                    t1     = d1(s11,s2,s3)
                    t2     = d1(s1,s21,s3)
                    t3     = d1(s1,s2,s31)

! field at position ii due to atom jj

                    If (Abs(jmp(nn)) > zero_plus) Then

                       fx = fx - tj*t1
                       fy = fy - tj*t2
                       fz = fz - tj*t3

                    End If

! field at position jj due to atom ii

                    If (Abs(imp(nn)) > zero_plus) Then

                       ti  = sx*alphan*imp(nn) ! multiply by sx to account for change in
                                               ! sign for odd derivatives

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

! calculate field

           mpfldx(l_ind(i))=mpfldx(l_ind(i))-fx
           mpfldy(l_ind(i))=mpfldy(l_ind(i))-fy
           mpfldz(l_ind(i))=mpfldz(l_ind(i))-fz

           mpfldx(l_ind(j))=mpfldx(l_ind(j))-tjx
           mpfldy(l_ind(j))=mpfldy(l_ind(j))-tjy
           mpfldz(l_ind(j))=mpfldz(l_ind(j))-tjz

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

! calculate interatomic distance

           rrr=Sqrt(xrr**2+yrr**2+zrr**2)

! get the multipoles for site jj

           jmp=mmp(:,jj)*scl

! calculate error function and derivative

           exp1=Exp(-(alpha*rrr)**2)
           tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

           Call ewald_deriv(-2,mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                    n      = s1+s2+s3+1 !+1 added ??? ITT

                    nn     = mplmap(s1,s2,s3)

                    alphan = alpha**n

                    tj     = sx*alphan*jmp(nn) ! sx added here for change of sign below to work ??? ITT

                    t1     = d1(s11,s2,s3)
                    t2     = d1(s1,s21,s3)
                    t3     = d1(s1,s2,s31)

! field at position ii due to atom jj

                    If (Abs(jmp(nn)) > zero_plus) Then

                       fx = fx - tj*t1
                       fy = fy - tj*t2
                       fz = fz - tj*t3

                    End If

! field at position jj due to atom ii

                    If (Abs(imp(nn)) > zero_plus) Then

                       ti  = sx*alphan*imp(nn) ! multiply by sx to account for change in
                                               ! sign for odd derivatives

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

! calculate field

           mpfldx(l_ind(i))=mpfldx(l_ind(i))-fx
           mpfldy(l_ind(i))=mpfldy(l_ind(i))-fy
           mpfldz(l_ind(i))=mpfldz(l_ind(i))-fz

        End Do
     End Do

     Deallocate (mmp,xfr,yfr,zfr, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mfield deallocation failure 1, node: ', idnode
        Call error(0)
     End If

  Else

! We resort to approximating N*(N-1)/2 interactions
! with the short-range one from the two body linked cell list

     Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mfield allocation failure 2, node: ', idnode
        Call error(0)
     End If

     Do ii=1,nz_fr(idnode+1)
        i=l_ind(nz_fr(idnode+1))
        idi=ltg(ii)

! get the multipoles for site i

        imp=mplgfr(:,i)*scl

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

! periodic boundary conditions not needed by LC construction
!
!           Call images(imcon,cell,limit,xxt,yyt,zzt)

! square of distances

           Do k=1,limit
              rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
           End Do

           Do k=1,limit
              j=list(list(-1,i)+k,i)

! get the multipoles for site j

              jmp=mplgfr(:,j)*scl

              rrr=rrt(k)
              If (rrr < rcut) Then

! calculate error function and derivative

                 exp1=Exp(-(alpha*rrr)**2)
                 tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

                 erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(alpha*rrr)

! compute derivatives of kernel

                 Call ewald_deriv(-2,mxompl+1,2,erfr,alpha*xrr,alpha*yrr,alpha*zrr,alpha*rrr,d1)

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

                          n      = s1+s2+s3+1 !+1 added ??? ITT

                          nn     = mplmap(s1,s2,s3)

                          alphan = alpha**n

                          tj     = sx*alphan*jmp(nn) ! sx added here for change of sign below to work ??? ITT

                          t1     = d1(s11,s2,s3)
                          t2     = d1(s1,s21,s3)
                          t3     = d1(s1,s2,s31)

! field at position i due to atom j

                          If (Abs(jmp(nn)) > zero_plus) Then

                             fx = fx - tj*t1
                             fy = fy - tj*t2
                             fz = fz - tj*t3

                          End If

! field at position j due to atom i

                          If (Abs(imp(nn)) > zero_plus) Then

                             ti  = sx*alphan*imp(nn) ! multiply by sx to account for change in
                                                     ! sign for odd derivatives

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

! calculate field

                 mpfldx(l_ind(i))=mpfldx(l_ind(i))-fx
                 mpfldy(l_ind(i))=mpfldy(l_ind(i))-fy
                 mpfldz(l_ind(i))=mpfldz(l_ind(i))-fz

                 If (j <= natms) Then


                    mpfldx(l_ind(i))=mpfldx(l_ind(i))-fx
                    mpfldy(l_ind(i))=mpfldy(l_ind(i))-fy
                    mpfldz(l_ind(i))=mpfldz(l_ind(i))-fz

                 End If

              End If

           End Do

        End If
     End Do

     Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mfield deallocation failure 2, node: ', idnode
        Call error(0)
     End If

  End If

  Deallocate (l_ind,nz_fr, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frzn_mfield deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine ewald_frzn_mfield
