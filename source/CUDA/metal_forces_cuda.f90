! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Centre for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

Subroutine metal_forces_helper &
          (iatm,xdf,ydf,zdf,rsqdf,rho,keypot,engmet,virmet,stress,safe)


  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltg,ltype,list,fxx,fyy,fzz
  Use metal_module,  Only : ntpmet,lstmet,ltpmet,vmet,dmet

  Implicit None

  Integer,                                  Intent( In    ) :: iatm,keypot
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: rho
  Real( Kind = wp ),                        Intent(   Out ) :: engmet,virmet
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Logical,                                  Intent( InOut ) :: safe

  Integer           :: m,idi,ai,aj,jatm,key,k0,l,ld
  Real( Kind = wp ) :: fix,fiy,fiz,fx,fy,fz,          &
                       rsq,rdr,rrr,ppp,ppd,           &
                       gk0,gk1,gk2,vk0,vk1,vk2,t1,t2, &
                       gamma,gamma1,gamma2,gamma3,    &
                       strs1,strs2,strs3,strs5,strs6,strs9

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity of atom iatm

  idi=ltg(iatm)

! start of primary loop for forces evaluation

  ai=ltype(iatm)

! load forces

  fix=0.0_wp
  fiy=0.0_wp
  fiz=0.0_wp

!$OMP DO
  Do m=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(m,iatm)
     aj=ltype(jatm)

     If (ai > aj) Then
        key=ai*(ai-1)/2 + aj
     Else
        key=aj*(aj-1)/2 + ai
     End If

     k0=lstmet(key)

! interatomic distance

     rsq = rsqdf(m)

! validity and truncation of potential

     If (ltpmet(k0) >= 0 .and. Abs(vmet(1,k0,1)) > zero_plus .and. rsq <= vmet(3,k0,1)**2) Then

! interpolation parameters

        rdr = 1.0_wp/vmet(4,k0,1)
        rrr = Sqrt(rsq) - vmet(2,k0,1)
        l   = Nint(rrr*rdr)
        ppp = rrr*rdr - Real(l,wp)

! catch unsafe value: EAM

        If (keypot == 0 .and. l < 2) Then
           safe=.false.
           l=2
        End If

! calculate pair forces using 3-point interpolation

        gk0 = vmet(l+3,k0,2)
        gk1 = vmet(l+4,k0,2)
        gk2 = vmet(l+5,k0,2)

        t1 = gk1 + ppp*(gk1 - gk0)
        t2 = gk1 + ppp*(gk2 - gk1)

        If (ppp < 0.0_wp) Then
           gamma1 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
        Else
           gamma1 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
        End If

! calculate embedding forces using 3-point interpolation

        If (keypot == 0) Then ! EAM

! contribution from first metal atom identity

           If (rsq <= dmet(3,aj,1)**2) Then

! interpolation parameters

              rdr = 1.0_wp/dmet(4,aj,1)
              rrr = Sqrt(rsq) - dmet(2,aj,1)
              ld  = Nint(rrr*rdr)
              ppd = rrr*rdr - Real(ld,wp)

! catch unsafe value: EAM

              If (ld < 2) Then
                 safe=.false.
                 ld=2
              End If

              gk0 = dmet(ld+3,aj,2)
              gk1 = dmet(ld+4,aj,2)
              gk2 = dmet(ld+5,aj,2)

              t1 = gk1 + ppd*(gk1 - gk0)
              t2 = gk1 + ppd*(gk2 - gk1)

              If (ppd < 0.0_wp) Then
                 gamma2 = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
              Else
                 gamma2 = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
              End If

           Else

              gamma2=0.0_wp

           End If

! contribution from second metal atom identity

           If (ai == aj) Then

              gamma3=gamma2

           Else

              If (rsq <= dmet(3,ai,1)**2) Then

! interpolation parameters

                 rdr = 1.0_wp/dmet(4,ai,1)
                 rrr = Sqrt(rsq) - dmet(2,ai,1)
                 ld  = Nint(rrr*rdr)
                 ppd = rrr*rdr - Real(ld,wp)

! catch unsafe value: EAM

                 If (ld < 2) Then
                    safe=.false.
                    ld=2
                 End If

                 gk0 = dmet(ld+3,ai,2)
                 gk1 = dmet(ld+4,ai,2)
                 gk2 = dmet(ld+5,ai,2)

                 t1 = gk1 + ppd*(gk1 - gk0)
                 t2 = gk1 + ppd*(gk2 - gk1)

                 If (ppd < 0.0_wp) Then
                    gamma3 = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                 Else
                    gamma3 = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                 End If

              Else

                 gamma3=0.0_wp

              End If

           End If

           gamma=(gamma1+(gamma2*rho(iatm)+gamma3*rho(jatm)))/rsq

        Else ! FST, interpolation parameters are the same for all force arrays

           gk0 = dmet(l+3,k0,2)
           gk1 = dmet(l+4,k0,2)
           gk2 = dmet(l+5,k0,2)

           t1 = gk1 + ppp*(gk1 - gk0)
           t2 = gk1 + ppp*(gk2 - gk1)

           If (ppp < 0.0_wp) Then
              gamma2 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
           Else
              gamma2 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
           End If

           If (ai > aj) Then
              gamma = (gamma1-gamma2*(rho(iatm)*dmet(1,k0,2)+rho(jatm)*dmet(2,k0,2)))/rsq
           Else
              gamma = (gamma1-gamma2*(rho(iatm)*dmet(2,k0,2)+rho(jatm)*dmet(1,k0,2)))/rsq
           End If

        End If

        fx = gamma*xdf(m)
        fy = gamma*ydf(m)
        fz = gamma*zdf(m)

        fix=fix+fx
        fiy=fiy+fy
        fiz=fiz+fz

        If (jatm <= natms) Then

           fxx(jatm)=fxx(jatm)-fx
           fyy(jatm)=fyy(jatm)-fy
           fzz(jatm)=fzz(jatm)-fz

        End If

        If (jatm <= natms .or. idi < ltg(jatm)) Then

! calculate interaction energy using 3-point interpolation

           vk0 = vmet(l+3,k0,1)
           vk1 = vmet(l+4,k0,1)
           vk2 = vmet(l+5,k0,1)

           t1 = vk1 + ppp*(vk1 - vk0)
           t2 = vk1 + ppp*(vk2 - vk1)

           If (ppp < 0.0_wp) Then
              engmet = engmet + t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
           Else
              engmet = engmet + t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
           End If

! calculate virial

           virmet = virmet - gamma*rsq

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
!$OMP END DO NOWAIT

! load back forces

!$OMP ATOMIC
  fxx(iatm)=fxx(iatm) + fix
!$OMP ATOMIC
  fyy(iatm)=fyy(iatm) + fiy
!$OMP ATOMIC
  fzz(iatm)=fzz(iatm) + fiz

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
End Subroutine metal_forces_helper

Subroutine metal_forces &
          (iatm,rmet,xdf,ydf,zdf,rsqdf,rho,engmet,virmet,stress,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating metal energy and force terms
! for EAM and FST interactions using verlet neighbour list
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov june 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltg,ltype,list,fxx,fyy,fzz
  Use metal_module,  Only : ld_met,ntpmet,lstmet,ltpmet,vmet,dmet,prmmet

#ifdef COMPILE_CUDA
  Use dl_poly_cuda_module
#endif

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: rho
  Real( Kind = wp ),                        Intent(   Out ) :: engmet,virmet
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Logical,                                  Intent( InOut ) :: safe

  Logical,           Save :: newjob = .true.
  Integer,           Save :: keypot
  Real( Kind = wp ), Save :: rcsq

  Integer           :: m,idi,ai,aj,jatm,key,kmn,kmx,k0,k1,k2,l,ld
  Real( Kind = wp ) :: fix,fiy,fiz,fx,fy,fz,          &
                       rsq,rdr,rrr,ppd,eng,           &
                       gk0,gk1,gk2,vk0,vk1,vk2,t1,t2, &
                       eps,sig,nnn,mmm,               &
                       cc0,cc1,cc2,cc3,cc4,           &
                       aaa,bbb,ccc,ddd,ppp,qqq,       &
                       bet,cut1,cut2,rr0,             &
                       gamma,gamma1,gamma2,gamma3,    &
                       strs1,strs2,strs3,strs5,strs6,strs9

#ifdef COMPILE_CUDA
  Call start_timing_metal_forces()
#endif

! check on mixing metal types

  If (newjob) Then
     newjob = .false.

     keypot = 0
     Do l=1,ntpmet
        keypot=ltpmet(l)
        If (l > 1) Then
           If (keypot /= ltpmet(l-1)) Call error(92)
        End If
     End Do

! set cutoff condition

     If (ld_met .and. keypot /= 0) rcsq=rmet**2
  End If

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

! global identity of atom iatm

  idi=ltg(iatm)

! start of primary loop for forces evaluation

  ai=ltype(iatm)

! load forces

  fix=fxx(iatm)
  fiy=fyy(iatm)
  fiz=fzz(iatm)

  Do m=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(m,iatm)
     aj=ltype(jatm)

     If (ai > aj) Then
        key=ai*(ai-1)/2 + aj
     Else
        key=aj*(aj-1)/2 + ai
     End If

     k0=lstmet(key)

     If (ld_met .and. keypot /= 0) Then
        k1=Max(ai,aj)
        k2=Min(ai,aj)

        kmx=k1*(k1-1)/2+k2
        kmn=k2*(k2-1)/2+k1

        k1=lstmet(kmx)
        k2=lstmet(kmn)
     End If

! interatomic distance

     rsq = rsqdf(m)

     If (ld_met .and. keypot /= 0) Then ! direct calculation

! validity and truncation of potential

        If (ltpmet(k0) > 0 .and. rsq < rcsq) Then

! interpolation parameters

           rrr = Sqrt(rsq)

           If      (ltpmet(k0) == 1) Then

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

              gamma1=0.0_wp
              eng   =0.0_wp
              If (rrr <= cut1) Then
                 gamma1 = -rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2) * &
                          (rrr-ccc)+(cc1+2.0_wp*cc2*rrr)*(rrr-ccc)**2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng = (cc0+cc1*rrr+cc2*rrr**2)*(rrr-ccc)**2
              End If

! calculate density contributions

              gamma2=0.0_wp
              If (rrr <= cut2) &
                 gamma2 = -rrr*(2.0_wp*(rrr-ddd)+3.0_wp*bet*(rrr-ddd)**2/ddd)

              If (ai == aj) Then
                 t1=prmmet(5,k0)**2
                 t2=t1
              Else
                 t1=prmmet(5,k1)**2
                 t2=prmmet(5,k2)**2
              End If

           Else If (ltpmet(k0) == 2) Then

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

              gamma1=0.0_wp
              eng   =0.0_wp
              If (rrr <= cut1) Then
                 gamma1 = -rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc) + &
                                (cc1+2.0_wp*cc2*rrr+3.0_wp*cc3*rrr**2+4.0_wp*cc4*rrr**3)*(rrr-ccc)**2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng = (cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc)**2
              End If

! calculate density contributions

              gamma2=0.0_wp
              If (rrr <= cut2) &
                 gamma2 = -rrr*(2.0_wp*(rrr-ddd)+4.0_wp*bbb*2*(rrr-ddd)**3)

              If (ai == aj) Then
                 t1=prmmet(7,k0)**2
                 t2=t1
              Else
                 t1=prmmet(7,k1)**2
                 t2=prmmet(7,k2)**2
              End If

           Else If (ltpmet(k0) == 3) Then

! sutton-chen potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              nnn=Nint(prmmet(3,k0))
              mmm=Nint(prmmet(4,k0))

! calculate pair forces and energies

              gamma1=Real(nnn,wp)*eps*(sig/rrr)**nnn
              eng   =0.0_wp
              If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng = gamma1/Real(nnn,wp)

! calculate density contributions

              gamma2=Real(mmm,wp)*(sig/rrr)**mmm

              If (ai == aj) Then
                 t1=(prmmet(1,k0)*prmmet(5,k0))**2
                 t2=t1
              Else
                 t1=(prmmet(1,k1)*prmmet(5,k1))**2
                 t2=(prmmet(1,k2)*prmmet(5,k2))**2
              End If

           Else If (ltpmet(k0) == 4) Then

! gupta potentials

              aaa=prmmet(1,k0)
              rr0=prmmet(2,k0)
              ppp=prmmet(3,k0)
              qqq=prmmet(5,k0)

              cut1=(rrr-rr0)/rr0
              cut2=cut1+1.0_wp

! calculate pair forces and energies

              gamma1=2.0_wp*aaa*Exp(-ppp*cut1)*ppp*cut2
              eng   =0.0_wp
              If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng = gamma1/(ppp*cut2)

! calculate density contributions

              gamma2=2.0_wp*Exp(-2.0_wp*qqq**cut1)*qqq*cut2

              t1=prmmet(4,k0)**2
              t2=t1

           End If

           If (ai > aj) Then
              gamma = (gamma1-gamma2*(rho(iatm)*t1+rho(jatm)*t2))/rsq
           Else
              gamma = (gamma1-gamma2*(rho(iatm)*t2+rho(jatm)*t1))/rsq
           End If

           fx = gamma*xdf(m)
           fy = gamma*ydf(m)
           fz = gamma*zdf(m)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

           End If

           If (jatm <= natms .or. idi < ltg(jatm)) Then

! calculate interaction energy

              engmet = engmet + eng

! calculate virial

              virmet = virmet - gamma*rsq

! calculate stress tensor

              strs1 = strs1 + xdf(m)*fx
              strs2 = strs2 + xdf(m)*fy
              strs3 = strs3 + xdf(m)*fz
              strs5 = strs5 + ydf(m)*fy
              strs6 = strs6 + ydf(m)*fz
              strs9 = strs9 + zdf(m)*fz

           End If

        End If

     Else ! tabulated calculation

! validity and truncation of potential

        If (ltpmet(k0) >= 0 .and. Abs(vmet(1,k0,1)) > zero_plus .and. rsq <= vmet(3,k0,1)**2) Then

! interpolation parameters

           rdr = 1.0_wp/vmet(4,k0,1)
           rrr = Sqrt(rsq) - vmet(2,k0,1)
           l   = Nint(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)

! catch unsafe value: EAM

           If (keypot == 0 .and. l < 2) Then
              safe=.false.
              l=2
           End If

! calculate pair forces using 3-point interpolation

           gk0 = vmet(l+3,k0,2)
           gk1 = vmet(l+4,k0,2)
           gk2 = vmet(l+5,k0,2)

           t1 = gk1 + ppp*(gk1 - gk0)
           t2 = gk1 + ppp*(gk2 - gk1)

           If (ppp < 0.0_wp) Then
              gamma1 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
           Else
              gamma1 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
           End If

! calculate embedding forces using 3-point interpolation

           If (keypot == 0) Then ! EAM

! contribution from first metal atom identity

              If (rsq <= dmet(3,aj,1)**2) Then

! interpolation parameters

                 rdr = 1.0_wp/dmet(4,aj,1)
                 rrr = Sqrt(rsq) - dmet(2,aj,1)
                 ld  = Nint(rrr*rdr)
                 ppd = rrr*rdr - Real(ld,wp)

! catch unsafe value: EAM

                 If (ld < 2) Then
                    safe=.false.
                    ld=2
                 End If

                 gk0 = dmet(ld+3,aj,2)
                 gk1 = dmet(ld+4,aj,2)
                 gk2 = dmet(ld+5,aj,2)

                 t1 = gk1 + ppd*(gk1 - gk0)
                 t2 = gk1 + ppd*(gk2 - gk1)

                 If (ppd < 0.0_wp) Then
                    gamma2 = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                 Else
                    gamma2 = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                 End If

              Else

                 gamma2=0.0_wp

              End If

! contribution from second metal atom identity

              If (ai == aj) Then

                 gamma3=gamma2

              Else

                 If (rsq <= dmet(3,ai,1)**2) Then

! interpolation parameters

                    rdr = 1.0_wp/dmet(4,ai,1)
                    rrr = Sqrt(rsq) - dmet(2,ai,1)
                    ld  = Nint(rrr*rdr)
                    ppd = rrr*rdr - Real(ld,wp)

! catch unsafe value: EAM

                    If (ld < 2) Then
                       safe=.false.
                       ld=2
                    End If

                    gk0 = dmet(ld+3,ai,2)
                    gk1 = dmet(ld+4,ai,2)
                    gk2 = dmet(ld+5,ai,2)

                    t1 = gk1 + ppd*(gk1 - gk0)
                    t2 = gk1 + ppd*(gk2 - gk1)

                    If (ppd < 0.0_wp) Then
                       gamma3 = t1 + 0.5_wp*(t2-t1)*(ppd+1.0_wp)
                    Else
                       gamma3 = t2 + 0.5_wp*(t2-t1)*(ppd-1.0_wp)
                    End If

                 Else

                    gamma3=0.0_wp

                 End If

              End If

              gamma=(gamma1+(gamma2*rho(iatm)+gamma3*rho(jatm)))/rsq

           Else ! FST, interpolation parameters are the same for all force arrays

              gk0 = dmet(l+3,k0,2)
              gk1 = dmet(l+4,k0,2)
              gk2 = dmet(l+5,k0,2)

              t1 = gk1 + ppp*(gk1 - gk0)
              t2 = gk1 + ppp*(gk2 - gk1)

              If (ppp < 0.0_wp) Then
                 gamma2 = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              Else
                 gamma2 = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              End If

              If (ai > aj) Then
                 gamma = (gamma1-gamma2*(rho(iatm)*dmet(1,k0,2)+rho(jatm)*dmet(2,k0,2)))/rsq
              Else
                 gamma = (gamma1-gamma2*(rho(iatm)*dmet(2,k0,2)+rho(jatm)*dmet(1,k0,2)))/rsq
              End If

           End If

           fx = gamma*xdf(m)
           fy = gamma*ydf(m)
           fz = gamma*zdf(m)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

           End If

           If (jatm <= natms .or. idi < ltg(jatm)) Then

! calculate interaction energy using 3-point interpolation

              vk0 = vmet(l+3,k0,1)
              vk1 = vmet(l+4,k0,1)
              vk2 = vmet(l+5,k0,1)

              t1 = vk1 + ppp*(vk1 - vk0)
              t2 = vk1 + ppp*(vk2 - vk1)

              If (ppp < 0.0_wp) Then
                 engmet = engmet + t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              Else
                 engmet = engmet + t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              End If

! calculate virial

              virmet = virmet - gamma*rsq

! calculate stress tensor

              strs1 = strs1 + xdf(m)*fx
              strs2 = strs2 + xdf(m)*fy
              strs3 = strs3 + xdf(m)*fz
              strs5 = strs5 + ydf(m)*fy
              strs6 = strs6 + ydf(m)*fz
              strs9 = strs9 + zdf(m)*fz

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

#ifdef COMPILE_CUDA
  Call stop_timing_metal_forces()
#endif

End Subroutine metal_forces
