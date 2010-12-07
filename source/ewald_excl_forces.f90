Subroutine ewald_excl_forces &
           (iatm,alpha,epsq,xdf,ydf,zdf,engcpe_ex,vircpe_ex,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using ewald's method
!
! Note: exclusion correction terms
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov february 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,nlast,lsi,lsa,lfrzn,chge, &
                            lexatm,fxx,fyy,fzz
  Use ewald_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: xdf,ydf,zdf
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

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

  Integer           :: jatm,m,local_index
  Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,alpr,alpr2, &
                       erfr,egamma,exp1,tt,             &
                       fix,fiy,fiz,fx,fy,fz,            &
                       strs1,strs2,strs3,strs5,strs6,strs9

! initialise potential energy and virial

  engcpe_ex=0.0_wp
  vircpe_ex=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! start of primary loop for forces evaluation

  chgea = chge(iatm)

! ignore interaction if the charge is zero

  If (Abs(chgea) > zero_plus) Then

! halve the charge of iatm since exclusions are double counted
! lexatm(:,iatm) points to jatm & lexatm(:,jatm) points to iatm

     chgea = 0.5_wp*chgea*r4pie0/epsq

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

     Do m=1,lexatm(0,iatm)

! atomic index

        jatm=local_index(lexatm(m,iatm),nlast,lsi,lsa)

! particles must be native or natively shared to this node (idnode)

        If (jatm > 0) Then
           chgprd=chge(jatm)

! ignore interaction the charge is zero
! ignore frozen pairs, they must not be dealt with here

           If (Abs(chgprd) > zero_plus .and. lfrzn(iatm)*lfrzn(jatm) == 0) Then

! charge product

              chgprd=chgprd*chgea

! calculate interatomic distance

              rsq=xdf(m)**2+ydf(m)**2+zdf(m)**2

              rrr  =Sqrt(rsq)
              alpr =rrr*alpha
              alpr2=alpr*alpr

! calculate error function and derivative

              If (alpr < 1.0e-2_wp) Then

! close particles (core-shell units) - small distances limit

                 erfr=2.0_wp*chgprd*(alpha/sqrpi) * &
                 (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

                 egamma=-4.0_wp*chgprd*(alpha**3/sqrpi) * &
                 (rr3+alpr2*(-2.0_wp*r10+alpr2*(3.0_wp*r42-4.0_wp*alpr2*r216)))

              Else

! distant particles - traditional

                 exp1 =Exp(-(alpha*rrr)**2)
                 tt   =1.0_wp/(1.0_wp+pp*alpha*rrr)

                 erfr=chgprd * &
                 (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

                 egamma=-(erfr-2.0_wp*chgprd*(alpha/sqrpi)*exp1)/rsq

              End If

! calculate forces

              fx = egamma*xdf(m)
              fy = egamma*ydf(m)
              fz = egamma*zdf(m)

              fix=fix+fx
              fiy=fiy+fy
              fiz=fiz+fz

! infrequent calculations copying

              If (l_cp) Then
                 fcx(iatm)=fcx(iatm)+fx
                 fcy(iatm)=fcy(iatm)+fy
                 fcz(iatm)=fcz(iatm)+fz
              End If

              If (jatm <= natms) Then

                 fxx(jatm)=fxx(jatm)-fx
                 fyy(jatm)=fyy(jatm)-fy
                 fzz(jatm)=fzz(jatm)-fz

! infrequent calculations copying

                 If (l_cp) Then
                    fcx(jatm)=fcx(jatm)-fx
                    fcy(jatm)=fcy(jatm)-fy
                    fcz(jatm)=fcz(jatm)-fz
                 End If

              Else

! account for the half inclusion of force corrections for iatm
! when natms < jatm <= nlast, i.e. iatm and jatm are on different
! nodes but see each other as a node seeing its halo

                 fix=fix+fx
                 fiy=fiy+fy
                 fiz=fiz+fz

! infrequent calculations copying

                 If (l_cp) Then
                    fcx(iatm)=fcx(iatm)+fx
                    fcy(iatm)=fcy(iatm)+fy
                    fcz(iatm)=fcz(iatm)+fz
                 End If

              End If

! calculate potential energy and virial

              engcpe_ex = engcpe_ex - erfr
              vircpe_ex = vircpe_ex - egamma*rsq

! infrequent calculations copying

              If (l_cp) Then
                 e_ex = e_ex - erfr
                 v_ex = v_ex - egamma*rsq
              End If

! calculate stress tensor

              strs1 = strs1 + xdf(m)*fx
              strs2 = strs2 + xdf(m)*fy
              strs3 = strs3 + xdf(m)*fz
              strs5 = strs5 + ydf(m)*fy
              strs6 = strs6 + ydf(m)*fz
              strs9 = strs9 + zdf(m)*fz

           End If
        End If

     End Do

! load back forces

     fxx(iatm)=fix
     fyy(iatm)=fiy
     fzz(iatm)=fiz

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
        s_ex(1) = s_ex(1) + strs1
        s_ex(2) = s_ex(2) + strs2
        s_ex(3) = s_ex(3) + strs3
        s_ex(4) = s_ex(4) + strs2
        s_ex(5) = s_ex(5) + strs5
        s_ex(6) = s_ex(6) + strs6
        s_ex(7) = s_ex(7) + strs3
        s_ex(8) = s_ex(8) + strs6
        s_ex(9) = s_ex(9) + strs9
     End If

  End If

End Subroutine ewald_excl_forces
