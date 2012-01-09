Subroutine ewald_frozen_forces &
           (imcon,rcut,alpha,epsq,keyens,engcpe_fr,vircpe_fr,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating corrections to coulombic forces
! in a periodic system arising from frozen pairs
!
! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
!       end (and any COM drift removed) but corrections to the stress
!       and the virial are important as they feed into the system
!       pressure response.  No volume changing ensembles (keyens < 20)
!       need this calculation just once!
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use config_module, Only : cell,natms,list,ltg,lfrzn,chge, &
                            xxx,yyy,zzz,fxx,fyy,fzz
  Use ewald_module

  Implicit None

  Integer,                             Intent( In    ) :: imcon,keyens
  Real( Kind = wp ),                   Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ),                   Intent(   Out ) :: engcpe_fr,vircpe_fr
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

  Logical, Save     :: newjob = .true., l_ens_do = .true.

  Logical           :: ll_cp
  Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
  Real( Kind = wp ) :: rcsq,det,rcell(1:9),xrr,yrr,zrr,rrr,rsq, &
                       chgprd,erfr,egamma,exp1,tt,              &
                       fx,fy,fz,xss,yss,zss,                    &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr, gfr
  Real( Kind = wp ), Dimension( : ), Allocatable :: cfr,xfr,yfr,zfr
  Real( Kind = wp ), Dimension( : ), Allocatable :: xdf,ydf,zdf,rsqdf

  fail=0
  Allocate (l_ind(1:mxatdm),nz_fr(0:mxnode), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces allocation failure, node: ', idnode
     Call error(0)
  End If

  If (.not.l_ens_do) Then
     Go To 100
  Else
     If (keyens < 20) Then
        ll_cp=l_cp
        l_cp =.true.
     End If
  End If

  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

     rcsq=rcut**2
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

  l_ind=0 ; nz_fr=0
  Do i=1,natms
     If (lfrzn(i) > 0 .and. Abs(chge(i)) > zero_plus) Then
        nz_fr(idnode+1)=nz_fr(idnode+1)+1
        l_ind(nz_fr(idnode+1))=i
     End If
  End Do
  If (mxnode > 1) Call gsum(nz_fr)
  nz_fr(0) = Sum(nz_fr(0:idnode)) ! Offset

  nzfr = Sum(nz_fr(1:mxnode))     ! Total
  If (nzfr <= 10*mxatms) Then

     Allocate (gfr(1:nzfr),cfr(1:nzfr),xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces allocation failure 1, node: ', idnode
        Call error(0)
     End If

     gfr=0
     cfr=0.0_wp
     xfr=0.0_wp
     yfr=0.0_wp
     zfr=0.0_wp
     Do i=1,nz_fr(idnode+1)
        ii=nz_fr(0)+i

        gfr(ii)=ltg(l_ind(i))
        cfr(ii)=chge(l_ind(i))
        xfr(ii)=xxx(l_ind(i))
        yfr(ii)=yyy(l_ind(i))
        zfr(ii)=zzz(l_ind(i))
     End Do
     If (mxnode > 1) Then
        Call gsum(gfr)
        Call gsum(cfr)
        Call gsum(xfr)
        Call gsum(yfr)
        Call gsum(zfr)
     End If

     Do i=1,nz_fr(idnode+1)
        ii=nz_fr(0)+i
        idi=gfr(ii) ! =ltg(l_ind(i))

        Do jj=1,nz_fr(0)
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

           rsq=xrr**2+yrr**2+zrr**2

           rrr=Sqrt(rsq)
           chgprd=cfr(ii)*cfr(jj)/epsq*r4pie0

! calculate error function and derivative

           exp1 =Exp(-(alpha*rrr)**2)
           tt   =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=chgprd * &
           (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

           egamma=-(erfr-2.0_wp*chgprd*(alpha/sqrpi)*exp1)/rsq

           fx = egamma*xrr
           fy = egamma*yrr
           fz = egamma*zrr

! calculate forces

           fxx(l_ind(i))=fxx(l_ind(i))-fx
           fyy(l_ind(i))=fyy(l_ind(i))-fy
           fzz(l_ind(i))=fzz(l_ind(i))-fz

! infrequent calculations copying

           If (l_cp) Then
              fcx(l_ind(i))=fcx(l_ind(i))-fx
              fcy(l_ind(i))=fcy(l_ind(i))-fy
              fcz(l_ind(i))=fcz(l_ind(i))-fz
           End If

           If (idi < gfr(jj)) Then

! calculate potential energy and virial

              engcpe_fr = engcpe_fr - erfr
              vircpe_fr = vircpe_fr - egamma*rsq

! calculate stress tensor

              strs1 = strs1 + xrr*fx
              strs2 = strs2 + xrr*fy
              strs3 = strs3 + xrr*fz
              strs5 = strs5 + yrr*fy
              strs6 = strs6 + yrr*fz
              strs9 = strs9 + zrr*fz

           End If
        End Do

        Do j=i+1,nz_fr(idnode+1)
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

           rsq=xrr**2+yrr**2+zrr**2

           rrr=Sqrt(rsq)
           chgprd=cfr(ii)*cfr(jj)/epsq*r4pie0

! calculate error function and derivative

           exp1 =Exp(-(alpha*rrr)**2)
           tt   =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=chgprd * &
           (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

           egamma=-(erfr-2.0_wp*chgprd*(alpha/sqrpi)*exp1)/rsq

           fx = egamma*xrr
           fy = egamma*yrr
           fz = egamma*zrr

! calculate forces

           fxx(l_ind(i))=fxx(l_ind(i))-fx
           fyy(l_ind(i))=fyy(l_ind(i))-fy
           fzz(l_ind(i))=fzz(l_ind(i))-fz

! infrequent calculations copying

           If (l_cp) Then
              fcx(l_ind(i))=fcx(l_ind(i))-fx
              fcy(l_ind(i))=fcy(l_ind(i))-fy
              fcz(l_ind(i))=fcz(l_ind(i))-fz
           End If

           If (l_ind(j) <= natms) Then

              fxx(l_ind(j))=fxx(l_ind(j))+fx
              fyy(l_ind(j))=fyy(l_ind(j))+fy
              fzz(l_ind(j))=fzz(l_ind(j))+fz

! infrequent calculations copying

              If (l_cp) Then
                 fcx(l_ind(j))=fcx(l_ind(j))+fx
                 fcy(l_ind(j))=fcy(l_ind(j))+fy
                 fcz(l_ind(j))=fcz(l_ind(j))+fz
              End If

           End If

           If (l_ind(j) <= natms .or. idi < ltg(l_ind(j))) Then

! calculate potential energy and virial

              engcpe_fr = engcpe_fr - erfr
              vircpe_fr = vircpe_fr - egamma*rsq

! calculate stress tensor

              strs1 = strs1 + xrr*fx
              strs2 = strs2 + xrr*fy
              strs3 = strs3 + xrr*fz
              strs5 = strs5 + yrr*fy
              strs6 = strs6 + yrr*fz
              strs9 = strs9 + zrr*fz

           End If
        End Do

        Do jj=nz_fr(0)+nz_fr(idnode+1)+1,nzfr
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

           rsq=xrr**2+yrr**2+zrr**2

           rrr=Sqrt(rsq)
           chgprd=cfr(ii)*cfr(jj)/epsq*r4pie0

! calculate error function and derivative

           exp1 =Exp(-(alpha*rrr)**2)
           tt   =1.0_wp/(1.0_wp+pp*alpha*rrr)

           erfr=chgprd * &
           (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

           egamma=-(erfr-2.0_wp*chgprd*(alpha/sqrpi)*exp1)/rsq

           fx = egamma*xrr
           fy = egamma*yrr
           fz = egamma*zrr

! calculate forces

           fxx(l_ind(i))=fxx(l_ind(i))-fx
           fyy(l_ind(i))=fyy(l_ind(i))-fy
           fzz(l_ind(i))=fzz(l_ind(i))-fz

! infrequent calculations copying

           If (l_cp) Then
              fcx(l_ind(i))=fcx(l_ind(i))-fx
              fcy(l_ind(i))=fcy(l_ind(i))-fy
              fcz(l_ind(i))=fcz(l_ind(i))-fz
           End If

           If (idi < gfr(jj)) Then

! calculate potential energy and virial

              engcpe_fr = engcpe_fr - erfr
              vircpe_fr = vircpe_fr - egamma*rsq

! calculate stress tensor

              strs1 = strs1 + xrr*fx
              strs2 = strs2 + xrr*fy
              strs3 = strs3 + xrr*fz
              strs5 = strs5 + yrr*fy
              strs6 = strs6 + yrr*fz
              strs9 = strs9 + zrr*fz

           End If
        End Do
     End Do

     Deallocate (gfr,cfr,xfr,yfr,zfr, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces deallocation failure 1, node: ', idnode
        Call error(0)
     End If

  Else

! We resort to approximating N*(N-1)/2 interactions
! with the short-range one from the two body linked cell list

     Allocate (xdf(1:mxlist),ydf(1:mxlist),zdf(1:mxlist),rsqdf(1:mxlist), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces allocation failure 2, node: ', idnode
        Call error(0)
     End If

     Do ii=1,nz_fr(idnode+1)
        i=l_ind(nz_fr(idnode+1))

! Get list limit

        limit=list(-2,i)-list(-1,i)
        If (limit > 0) Then

! calculate interatomic distances

           Do k=1,limit
              j=list(list(-1,i)+k,i)

              xdf(k)=xxx(i)-xxx(j)
              ydf(k)=yyy(i)-yyy(j)
              zdf(k)=zzz(i)-zzz(j)
           End Do

! periodic boundary conditions

           Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distances

           Do k=1,limit
              rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
           End Do

           Do k=1,limit
              j=list(list(-1,i)+k,i)

              rsq=rsqdf(k)
              If (rsq < rcsq .and. Abs(chge(j)) > zero_plus) Then
                 rrr=Sqrt(rsq)
                 chgprd=chge(i)*chge(j)/epsq*r4pie0

! calculate error function and derivative

                 exp1 =Exp(-(alpha*rrr)**2)
                 tt   =1.0_wp/(1.0_wp+pp*alpha*rrr)

                 erfr=chgprd * &
                 (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

                 egamma=-(erfr-2.0_wp*chgprd*(alpha/sqrpi)*exp1)/rsq

                 fx = egamma*xrr
                 fy = egamma*yrr
                 fz = egamma*zrr

! calculate forces

                 fxx(i)=fxx(i)-fx
                 fyy(i)=fyy(i)-fy
                 fzz(i)=fzz(i)-fz

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

! infrequent calculations copying

                    If (l_cp) Then
                       fcx(j)=fcx(j)+fx
                       fcy(j)=fcy(j)+fy
                       fcz(j)=fcz(j)+fz
                    End If

                 End If

                 If (j <= natms .or. idi < ltg(j)) Then

! calculate potential energy and virial

                    engcpe_fr = engcpe_fr - erfr
                    vircpe_fr = vircpe_fr - egamma*rsq

! calculate stress tensor

                    strs1 = strs1 + xrr*fx
                    strs2 = strs2 + xrr*fy
                    strs3 = strs3 + xrr*fz
                    strs5 = strs5 + yrr*fy
                    strs6 = strs6 + yrr*fz
                    strs9 = strs9 + zrr*fz

                 End If
              End If
           End Do

        End If
     End Do

     Deallocate (xdf,ydf,zdf,rsqdf, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces deallocation failure 2, node: ', idnode
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

100 Continue
  If (keyens < 20) Then
     If (l_ens_do) Then ! first round
        l_ens_do=.false.

        l_cp=ll_cp
     Else               ! consecutive rounds
        engcpe_fr = e_fr
        vircpe_fr = v_fr

        stress(1) = s_fr(1)
        stress(2) = s_fr(2)
        stress(3) = s_fr(3)
        stress(4) = s_fr(4)
        stress(5) = s_fr(5)
        stress(6) = s_fr(6)
        stress(7) = s_fr(7)
        stress(8) = s_fr(8)
        stress(9) = s_fr(9)
     End If
  End If

  Deallocate (l_ind,nz_fr, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine ewald_frozen_forces
