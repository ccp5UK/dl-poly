Subroutine metal_forces &
          (iatm,rmet,xxt,yyt,zzt,rrt,engmet,virmet,stress,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating metal energy and force terms
! for EAM and FST interactions using verlet neighbour list
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov january 2016
! contrib   - r.davidchak (eeam) june 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use site_module,   Only : ntpatm
  Use config_module, Only : natms,ltg,ltype,list,fxx,fyy,fzz
  Use metal_module,  Only : ld_met,l2bmet,tabmet,lstmet,ltpmet, &
                            vmet,dmet,dmes,prmmet,rho,rhs,merf,mfer

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engmet,virmet
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Logical,                                  Intent( InOut ) :: safe

  Integer           :: m,idi,ai,ki,jatm,aj,kj, &
                       key,kmn,kmx,k0,keypot,  &
                       k1,k2,nnn,mmm,l,ld
  Real( Kind = wp ) :: fix,fiy,fiz,fx,fy,fz,        &
                       rrr,rsq,rdr,rr1,ppd,eng,     &
                       gk0,gk1,gk2,vk0,vk1,vk2,     &
                       t1,t2,t3,t4,gam1,gam2,       &
                       eps,sig,nnnr,mmmr,           &
                       cc0,cc1,cc2,cc3,cc4,         &
                       aaa,bbb,ccc,ddd,ppp,qqq,     &
                       bet,cut1,cut2,rr0,           &
                       gamma,gamma1,                &
                       gamma2,gamma3,gamm2s,gamm3s, &
                       strs1,strs2,strs3,strs5,strs6,strs9

! initialise potential energy and virial

  engmet=0.0_wp
  virmet=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity and type of atom iatm

  idi=ltg(iatm)
  ai=ltype(iatm)

! load forces

  fix=fxx(iatm)
  fiy=fyy(iatm)
  fiz=fzz(iatm)

! start of primary loop for forces evaluation

  Do m=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(m,iatm)
     aj=ltype(jatm)

     If      (tabmet == 1 .or. tabmet == 3) Then ! EAM & 2BEAM
        ki=ai
        kj=aj
     Else If (tabmet == 2 .or. tabmet == 4) Then ! EEAM & 2BEEAM
        ki=(aj-1)*ntpatm+ai ! aj-ai
        kj=(ai-1)*ntpatm+aj ! ai-aj
     End If

     If (ai > aj) Then
        key=ai*(ai-1)/2 + aj
     Else
        key=aj*(aj-1)/2 + ai
     End If

     k0=lstmet(key)

     If (ld_met) Then
        k1=Max(ai,aj)
        k2=Min(ai,aj)

        kmx=k1*(k1+1)/2
        kmn=k2*(k2+1)/2

        k1=lstmet(kmx)
        k2=lstmet(kmn)
     End If

! interatomic distance

     rrr = rrt(m)

! truncation and validity of metal interaction

     keypot=ltpmet(k0)
     If (keypot >= 0 .and. rrr <= rmet) Then

! Squared distance

        rsq = rrr**2

! Zero energy and force components

        eng   = 0.0_wp
        gamma1= 0.0_wp
        gamma2= 0.0_wp
        gamma3= 0.0_wp
        gamm2s= 0.0_wp
        gamm3s= 0.0_wp
        gamma = 0.0_wp

        If (ld_met) Then ! direct calculation (keypot /= 0)

! Type of analytic potential

           If      (keypot == 1) Then

! finnis-sinclair potentials

              cc0=prmmet(1,k0)
              cc1=prmmet(2,k0)
              cc2=prmmet(3,k0)
              ccc=prmmet(4,k0)
              ddd=prmmet(6,k0)
              bet=prmmet(7,k0)
              cut1=ccc
              cut2=ddd

! calculate pair forces and energies

              If (rrr <= cut1) Then
                 gamma1 = -rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2) * &
                          (rrr-ccc)+(cc1+2.0_wp*cc2*rrr)*(rrr-ccc)**2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng = (cc0+cc1*rrr+cc2*rrr**2)*(rrr-ccc)**2
              End If

! calculate density contributions

              If (rrr <= cut2) &
                 gamma2 = -rrr*(2.0_wp*(rrr-ddd)+3.0_wp*bet*(rrr-ddd)**2/ddd)

              If (ai == aj) Then
                 t1=prmmet(5,k0)**2
                 t2=t1
              Else
                 t1=prmmet(5,k1)**2
                 t2=prmmet(5,k2)**2
              End If

           Else If (keypot == 2) Then

! extended finnis-sinclair potentials

              cc0=prmmet(1,k0)
              cc1=prmmet(2,k0)
              cc2=prmmet(3,k0)
              cc3=prmmet(4,k0)
              cc4=prmmet(5,k0)
              ccc=prmmet(6,k0)
              ddd=prmmet(8,k0)
              bbb=prmmet(9,k0)
              cut1=ccc
              cut2=ddd

! calculate pair forces and energies

              If (rrr <= cut1) Then
                 gamma1 = -rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc) + &
                                (cc1+2.0_wp*cc2*rrr+3.0_wp*cc3*rrr**2+4.0_wp*cc4*rrr**3)*(rrr-ccc)**2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng = (cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc)**2
              End If

! calculate density contributions

              If (rrr <= cut2) &
                 gamma2 = -rrr*(2.0_wp*(rrr-ddd)+4.0_wp*bbb*2*(rrr-ddd)**3)

              If (ai == aj) Then
                 t1=prmmet(7,k0)**2
                 t2=t1
              Else
                 t1=prmmet(7,k1)**2
                 t2=prmmet(7,k2)**2
              End If

           Else If (keypot == 3) Then

! sutton-chen potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              nnn=Nint(prmmet(3,k0)) ; nnnr=Real(nnn,wp)
              mmm=Nint(prmmet(4,k0)) ; mmmr=Real(mmm,wp)

! calculate pair forces and energies

              gamma1=nnnr*eps*(sig/rrr)**nnn
              If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng = gamma1/nnnr

! calculate density contributions

              gamma2=mmmr*(sig/rrr)**mmm

              If (ai == aj) Then
                 t1=(prmmet(1,k0)*prmmet(5,k0))**2
                 t2=t1
              Else
                 t1=(prmmet(1,k1)*prmmet(5,k1))**2
                 t2=(prmmet(1,k2)*prmmet(5,k2))**2
              End If

           Else If (keypot == 4) Then

! gupta potentials

              aaa=prmmet(1,k0)
              rr0=prmmet(2,k0)
              ppp=prmmet(3,k0)
              qqq=prmmet(5,k0)

              cut1=(rrr-rr0)/rr0
              cut2=cut1+1.0_wp

! calculate pair forces and energies

              gamma1=2.0_wp*aaa*Exp(-ppp*cut1)*ppp*cut2
              If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng = gamma1/(ppp*cut2)

! calculate density contributions

              gamma2=2.0_wp*Exp(-2.0_wp*qqq*cut1)*qqq*cut2

              t1=prmmet(4,k0)**2
              t2=t1

           Else If (keypot == 5) Then

! many-body perturbation component only potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              mmm=Nint(prmmet(3,k0)) ; mmmr=Real(mmm,wp)

! no pair forces and energies

!             gamma1=0.0_wp
!              If (jatm <= natms .or. idi < ltg(jatm)) &
!                 eng = 0.0_wp

! calculate density contributions

! interpolation parameters

              rdr = 1.0_wp/merf(4)
              rr1 = rrr - merf(2)
              l   = Min(Nint(rr1*rdr),Nint(merf(1))-1)
              If (l < 5) Then ! catch unsafe value
                 safe=.false.
                 l=6
              End If
              ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

              vk0 = merf(l-1)
              vk1 = merf(l  )
              vk2 = merf(l+1)

              t1 = vk1 + ppp*(vk1 - vk0)
              t2 = vk1 + ppp*(vk2 - vk1)

              gk0 = mfer(l-1)
              gk1 = mfer(l  )
              gk2 = mfer(l+1)

              t3 = gk1 + ppp*(gk1 - gk0)
              t4 = gk1 + ppp*(gk2 - gk1)

              If (ppp < 0.0_wp) Then
                 gam1 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                 gam2 = t3 + 0.5_wp*(t4-t3)*(ppp+1.0_wp)
              Else If (l == 5) Then
                 gam1 = t2
                 gam2 = t4
              Else
                 gam1 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                 gam2 = t4 + 0.5_wp*(t4-t3)*(ppp-1.0_wp)
              End If

              gamma2=(sig/rrr**mmm)*(mmmr*gam1-rrr*gam2)

              If (ai == aj) Then
                 t1=prmmet(1,k0)**2
                 t2=t1
              Else
                 t1=prmmet(1,k1)**2
                 t2=prmmet(1,k2)**2
              End If

           End If

           If (ai > aj) Then
              gamma = (gamma1 - gamma2*(rho(iatm)*t1+rho(jatm)*t2))/rsq
           Else
              gamma = (gamma1 - gamma2*(rho(iatm)*t2+rho(jatm)*t1))/rsq
           End If

           fx = gamma*xxt(m)
           fy = gamma*yyt(m)
           fz = gamma*zzt(m)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

           End If

           If (jatm <= natms .or. idi < ltg(jatm)) Then

! add interaction energy

              engmet = engmet + eng

! add virial

              virmet = virmet - gamma*rsq

! add stress tensor

              strs1 = strs1 + xxt(m)*fx
              strs2 = strs2 + xxt(m)*fy
              strs3 = strs3 + xxt(m)*fz
              strs5 = strs5 + yyt(m)*fy
              strs6 = strs6 + yyt(m)*fz
              strs9 = strs9 + zzt(m)*fz

           End If

        Else ! tabulated calculation

! truncation of potential

           If (Abs(vmet(1,k0,1)) > zero_plus) Then

! interpolation parameters

              If (rrr <= vmet(3,k0,1) .or. & ! Next covers the FST density!
                  (keypot /= 0 .and. rrr <= dmet(3,k0,1))) Then

                 rdr = 1.0_wp/vmet(4,k0,1)
                 rr1 = rrr - vmet(2,k0,1)
                 l   = Min(Nint(rr1*rdr),Nint(vmet(1,k0,1))-1)
                 If (l < 5) Then ! catch unsafe value
                    safe=.false.
                    l=6
                 End If
                 ppp = rr1*rdr - Real(l,wp)

              End If

! calculate pair forces using 3-point interpolation

              If (rrr <= vmet(3,k0,1)) Then

                 gk0 = vmet(l-1,k0,2)
                 gk1 = vmet(l  ,k0,2)
                 gk2 = vmet(l+1,k0,2)

                 t1 = gk1 + ppp*(gk1 - gk0)
                 t2 = gk1 + ppp*(gk2 - gk1)

                 If (ppp < 0.0_wp) Then
                    gamma1 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                 Else If (l == 5) Then
                    gamma1 = t2
                 Else
                    gamma1 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                 End If

! calculate interaction energy using 3-point interpolation

                 If ((jatm <= natms .or. idi < ltg(jatm)) .and. keypot /= 5) Then

                    vk0 = vmet(l-1,k0,1)
                    vk1 = vmet(l  ,k0,1)
                    vk2 = vmet(l+1,k0,1)

                    t1 = vk1 + ppp*(vk1 - vk0)
                    t2 = vk1 + ppp*(vk2 - vk1)

                    If (ppp < 0.0_wp) Then
                       eng = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                    Else If (l == 5) Then
                       eng = t2
                    Else
                       eng = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                    End If

                 End If

              End If

           End If

! calculate embedding forces using 3-point interpolation

           If (keypot == 0) Then ! EAM

! contribution from first metal atom identity

              If (Abs(dmet(1,kj,1)) > zero_plus .and. Nint(dmet(1,ki,1)) > 5) Then
                 If (rrr <= dmet(3,kj,1)) Then

! interpolation parameters

                    rdr = 1.0_wp/dmet(4,kj,1)
                    rr1 = rrr - dmet(2,kj,1)
                    ld  = Min(Nint(rr1*rdr),Nint(dmet(1,kj,1))-1)
                    If (ld < 5) Then ! catch unsafe value: EAM
                       safe=.false.
                       ld=6
                    End If
                    ppd = rr1*rdr - Real(ld,wp)

                    gk0 = dmet(ld-1,kj,2)
                    gk1 = dmet(ld  ,kj,2)
                    gk2 = dmet(ld+1,kj,2)

                    t1 = gk1 + ppd*(gk1 - gk0)
                    t2 = gk1 + ppd*(gk2 - gk1)

                    If (ppd < 0.0_wp) Then
                       gamma2 = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                    Else If (ld == 5) Then
                       gamma2 = t2
                    Else
                       gamma2 = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                    End If

                 End If
              End If

! Now if we have 2B(EAM & EEAM) then do s-band too

              If (l2bmet) Then
                 If (Abs(dmes(1,kj,1)) > zero_plus .and. Nint(dmes(1,ki,1)) > 5) Then
                    If (rrr <= dmes(3,kj,1)) Then

! interpolation parameters

                       rdr = 1.0_wp/dmes(4,kj,1)
                       rr1 = rrr - dmes(2,kj,1)
                       ld  = Min(Nint(rr1*rdr),Nint(dmes(1,kj,1))-1)
                       If (ld < 5) Then ! catch unsafe value: EAM
                          safe=.false.
                          ld=6
                       End If
                       ppd = rr1*rdr - Real(ld,wp)

                       gk0 = dmes(ld-1,kj,2)
                       gk1 = dmes(ld  ,kj,2)
                       gk2 = dmes(ld+1,kj,2)

                       t1 = gk1 + ppd*(gk1 - gk0)
                       t2 = gk1 + ppd*(gk2 - gk1)

                       If (ppd < 0.0_wp) Then
                          gamm2s = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                       Else If (ld == 5) Then
                          gamm2s = t2
                       Else
                          gamm2s = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                       End If

                    End If
                 End If
              End If

! contribution from second metal atom identity

              If (ki == kj) Then

                 gamma3=gamma2
                 If (l2bmet) gamm3s=gamm2s !2B(EAM & EEAM)

              Else

                 If (Abs(dmet(1,ki,1)) > zero_plus .and. Nint(dmet(1,ki,1)) > 5) Then
                    If (rrr <= dmet(3,ki,1)) Then

! interpolation parameters

                       rdr = 1.0_wp/dmet(4,ki,1)
                       rr1 = rrr - dmet(2,ki,1)
                       ld  = Min(Nint(rr1*rdr),Nint(dmet(1,ki,1))-1)
                       If (ld < 5) Then ! catch unsafe value: EAM
                          safe=.false.
                          ld=6
                       End If
                       ppd = rr1*rdr - Real(ld,wp)

                       gk0 = dmet(ld-1,ki,2)
                       gk1 = dmet(ld  ,ki,2)
                       gk2 = dmet(ld+1,ki,2)

                       t1 = gk1 + ppd*(gk1 - gk0)
                       t2 = gk1 + ppd*(gk2 - gk1)

                       If (ppd < 0.0_wp) Then
                          gamma3 = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                       Else If (ld == 5) Then
                          gamma3 = t2
                       Else
                          gamma3 = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                       End If

                    End If
                 End If

                 If (l2bmet) Then !2B(EAM & EEAM)
                    If (Abs(dmes(1,ki,1)) > zero_plus .and. Nint(dmes(1,ki,1)) > 5) Then
                       If (rrr <= dmes(3,ki,1)) Then

! interpolation parameters

                          rdr = 1.0_wp/dmes(4,ki,1)
                          rr1 = rrr - dmes(2,ki,1)
                          ld  = Min(Nint(rr1*rdr),Nint(dmes(1,ki,1))-1)
                          If (ld < 5) Then ! catch unsafe value: EAM
                             safe=.false.
                             ld=6
                          End If
                          ppd = rr1*rdr - Real(ld,wp)

                          gk0 = dmes(ld-1,ki,2)
                          gk1 = dmes(ld  ,ki,2)
                          gk2 = dmes(ld+1,ki,2)

                          t1 = gk1 + ppd*(gk1 - gk0)
                          t2 = gk1 + ppd*(gk2 - gk1)

                          If (ppd < 0.0_wp) Then
                             gamm3s = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                          Else If (ld == 5) Then
                             gamm3s = t2
                          Else
                             gamm3s = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                          End If

                       End If
                    End If
                 End If

              End If

              If (.not.l2bmet) Then
                 gamma=(gamma1 + (gamma2*rho(iatm)+gamma3*rho(jatm)))/rsq
              Else !2B(EAM & EEAM)
                 gamma=(gamma1 + (gamma2*rho(iatm)+gamma3*rho(jatm)) &
                               + (gamm2s*rhs(iatm)+gamm3s*rhs(jatm)))/rsq
              End If

           Else ! FST, interpolation parameters are the same for all force arrays

              If (rrr <= dmet(3,k0,1)) Then ! interpolation parameters covered above

                 gk0 = dmet(l-1,k0,2)
                 gk1 = dmet(l  ,k0,2)
                 gk2 = dmet(l+1,k0,2)

                 t1 = gk1 + ppp*(gk1 - gk0)
                 t2 = gk1 + ppp*(gk2 - gk1)

                 If (ppp < 0.0_wp) Then
                    gamma2 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                 Else If (l == 5) Then
                    gamma2 = t2
                 Else
                    gamma2 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                 End If

              End If

              If (ai > aj) Then
                 gamma = (gamma1 - gamma2*(rho(iatm)*dmet(1,k0,2)+rho(jatm)*dmet(2,k0,2)))/rsq
              Else
                 gamma = (gamma1 - gamma2*(rho(iatm)*dmet(2,k0,2)+rho(jatm)*dmet(1,k0,2)))/rsq
              End If

           End If

           fx = gamma*xxt(m)
           fy = gamma*yyt(m)
           fz = gamma*zzt(m)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

           End If

           If (jatm <= natms .or. idi < ltg(jatm)) Then

! add interaction energy using 3-point interpolation

              engmet = engmet + eng

! add virial

              virmet = virmet - gamma*rsq

! add stress tensor

              strs1 = strs1 + xxt(m)*fx
              strs2 = strs2 + xxt(m)*fy
              strs3 = strs3 + xxt(m)*fz
              strs5 = strs5 + yyt(m)*fy
              strs6 = strs6 + yyt(m)*fz
              strs9 = strs9 + zzt(m)*fz

           End If

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

End Subroutine metal_forces
