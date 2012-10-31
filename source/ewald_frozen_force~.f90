Subroutine ewald_frozen_forces &
           (imcon,rcut,alpha,epsq,megfrz,engcpe_fr,vircpe_fr,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating corrections to coulombic forces
! in a periodic system arising from frozen pairs
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use domains_module, Only : map
  Use config_module,  Only : cell,natms,nlast,ltg,lfrzn,chge, &
                            xxx,yyy,zzz,fxx,fyy,fzz
  Use ewald_module

  Implicit None

  Integer,                             Intent( In    ) :: imcon,megfrz
  Real( Kind = wp ),                   Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ),                   Intent(   Out ) :: engcpe_fr,vircpe_fr
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

  Logical,           Save :: newjob = .true. , &
                             unsafe
  Real( Kind = wp ), Save :: rcsq

  Integer           :: fail,i,j,ii,jj,idi,numfrz,allfrz
  Real( Kind = wp ) :: det,rcell(1:9),xrr,yrr,zrr,rrr,rsq, &
                       chgprd,erfr,egamma,exp1,tt,         &
                       fx,fy,fz,xss,yss,zss,               &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Integer, Dimension( : ), Allocatable :: l_index

  fail=0
  Allocate (l_index(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! unsafe to omit repetition check

     unsafe=(Any(map == idnode))

     rcsq=rcut**2
  End If

! list domain+halo non-repeated frozen particles

  numfrz=0
  allfrz=0
  l_index=0
  Do i=1,nlast
     If (lfrzn(i) > 0) Then
        ii=ltg(i)
        If (i > natms) Then
           If (unsafe) Then
              If (Any(ltg(l_index(1:numfrz)) == ii)) Go To 10
           End If
        Else
           allfrz=allfrz+1
        End If

        numfrz=numfrz+1
        l_index(numfrz)=i
10      Continue
     End If
  End Do

  If (mxnode > 1) Call gsum(allfrz)
  If (allfrz /= megfrz) Call error(0)

  Call invert(cell,rcell,det)

! Go through all frozen pairs on this node
! Initialise contributions

  engcpe_fr=0.0_wp
  vircpe_fr=0.0_wp

  strs1 = 0.0_wp
  strs2 = 0.0_wp
  strs3 = 0.0_wp
  strs5 = 0.0_wp
  strs6 = 0.0_wp
  strs9 = 0.0_wp

  Do i=1,numfrz
     ii=l_index(i)

! for all native and non-zero charged frozen particles
! by construction l_index is sorted

     If (ii <= natms .and. Abs(chge(ii)) > zero_plus) Then

        idi=ltg(ii)
        Do j=i+1,numfrz
           jj=l_index(j)

           If (Abs(chge(jj)) > zero_plus) Then

              xrr=xxx(ii)-xxx(jj)
              yrr=yyy(ii)-yyy(jj)
              zrr=zzz(ii)-zzz(jj)

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

              If (rsq < rcsq) Then

                 rrr=Sqrt(rsq)
                 chgprd=chge(ii)*chge(jj)/epsq*r4pie0

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

                 fxx(ii)=fxx(ii)-fx
                 fyy(ii)=fyy(ii)-fy
                 fzz(ii)=fzz(ii)-fz

! infrequent calculations copying

                 If (l_cp) Then
                    fcx(ii)=fcx(ii)-fx
                    fcy(ii)=fcy(ii)-fy
                    fcz(ii)=fcz(ii)-fz
                 End If

                 If (jj <= natms) Then

                    fxx(jj)=fxx(jj)+fx
                    fyy(jj)=fyy(jj)+fy
                    fzz(jj)=fzz(jj)+fz

! infrequent calculations copying

                    If (l_cp) Then
                       fcx(jj)=fcx(jj)+fx
                       fcy(jj)=fcy(jj)+fy
                       fcz(jj)=fcz(jj)+fz
                    End If

                 End If

                 If (jj <= natms .or. idi < ltg(jj)) Then

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

           End If

        End Do

     End If

  End Do

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

  fail=0
  Deallocate (l_index, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_frozen_forces allocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine ewald_frozen_forces
