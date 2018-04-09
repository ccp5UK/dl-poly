Module metal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global metal interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use setup
  Use site,   Only : ntpatm,unqatm,dens
  Use configuration, Only : natms,ltg,ltype,list,fxx,fyy,fzz,&
                            xxx,yyy,zzz,imcon,volm,nlast,ixyz

  Use comms,  Only : comms_type,gsum,gcheck,gmax,MetLdExp_tag,wp_mpi,gsend,gwait
  Use parse, Only : get_line,get_word,lower_case,word_2_real
  Use domains, Only : map

  Use errors_warnings, Only : error,warning, info
  Use numerics, Only : erfgen_met
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif
  Implicit None

  Logical,                        Save :: ld_met = .false., & ! no direct calculations are opted
                                          ls_met = .false., & ! no embedding over Sqrt(rho) but over rho
                                          l2bmet = .false.    ! no 2B(EAM or EEAM)

  Integer,                        Save :: ntpmet = 0 , & ! number of different metal interactions
                                          tabmet = -1    ! undefined, 0 - no TABEAM, 1 - EAM, 2 - EEAM, 3- 2BEAM, 4 - 2BEEAM


  Integer,           Allocatable, Save :: lstmet(:),ltpmet(:)

  Real( Kind = wp ), Allocatable, Save :: prmmet(:,:)

  Real( Kind = wp ), Allocatable, Save :: elrcm(:),vlrcm(:)

! Possible tabulated calculation arrays

  Real( Kind = wp ), Allocatable, Save :: vmet(:,:,:), dmet(:,:,:),dmes(:,:,:), &
                                          fmet(:,:,:),fmes(:,:,:)

! Atomic density [reused as embedding derivative(s)] helper array(s)

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: rho,rhs

! Many-body perturbation potential error function and derivative arrays

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: merf,mfer


  Public :: allocate_metal_arrays, &
            allocate_metal_table_arrays

Contains

  Subroutine allocate_metal_arrays()

    Integer, Dimension( 1:7 ) :: fail

    If (tabmet == 3 .or. tabmet == 4) l2bmet=.true.

    fail = 0

    Allocate (lstmet(1:mxmet),                  Stat = fail(1))
    Allocate (ltpmet(1:mxmet),                  Stat = fail(2))
    Allocate (prmmet(1:mxpmet,1:mxmet),         Stat = fail(3))
    Allocate (rho(1:Merge(mxatms,0,mxmet > 0)), Stat = fail(4))
    If (l2bmet) & ! the new S-band density
    Allocate (rhs(1:Merge(mxatms,0,mxmet > 0)), Stat = fail(5))
    Allocate (elrcm(0:mxatyp),                  Stat = fail(6))
    Allocate (vlrcm(0:mxatyp),                  Stat = fail(7))

    If (Any(fail > 0)) Call error(1023)

    lstmet = 0
    ltpmet = 0

    prmmet = 0.0_wp

! rho and rhs get initialised in metal_ld_compute!!!

    elrcm  = 0.0_wp
    vlrcm  = 0.0_wp

  End Subroutine allocate_metal_arrays

  Subroutine allocate_metal_table_arrays()

    Integer, Dimension( 1:5 ) :: fail

    fail = 0

    Allocate (vmet(1:mxgmet,1:mxmet, 1:2),    Stat = fail(1))
    Allocate (dmet(1:mxgmet,1:mxmed, 1:2),    Stat = fail(2))
    Allocate (fmet(1:mxgmet,1:mxatyp,1:2),    Stat = fail(3))
    If (tabmet == 3 .or. tabmet == 4) Then ! the new S-band density and embedding
       Allocate (dmes(1:mxgmet,1:mxmds, 1:2), Stat = fail(4))
       Allocate (fmes(1:mxgmet,1:mxatyp,1:2), Stat = fail(5))
    End If

    If (Any(fail > 0)) Call error(1069)

    vmet = 0.0_wp
    dmet = 0.0_wp
    fmet = 0.0_wp
    If (tabmet == 3 .or. tabmet == 4) Then ! the new S-band density and embedding
      dmes = 0.0_wp
      fmes = 0.0_wp
    End If

  End Subroutine allocate_metal_table_arrays

  Subroutine allocate_metal_erf_arrays()

    Integer :: fail

    fail = 0

    Allocate (merf(1:mxgmet),mfer(1:mxgmet), Stat = fail)

    If (fail > 0) Call error(1082)

    merf = 0.0_wp
    mfer = 0.0_wp

  End Subroutine allocate_metal_erf_arrays
  
  Subroutine metal_forces &
          (iatm,rmet,xxt,yyt,zzt,rrt,engmet,virmet,stress,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating metal energy and force terms
! for EAM and FST interactions using verlet neighbour list
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov november 2016
! contrib   - r.davidchak (eeam) june 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
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
                 gamma2 = -rrr*(2.0_wp*(rrr-ddd)+4.0_wp*bbb**2*(rrr-ddd)**3)

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

              If (rrr <= vmet(3,k0,1) .or. & ! Next covers the FST density - merge used to avoid table check beyond bound!
                  (keypot /= 0 .and. rrr <= dmet(3,Merge(k0,1,keypot /= 0),1))) Then

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

Subroutine metal_ld_compute(rmet,elrcm,vlrcm,engden,virden,stress,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local density in metals using
! the verlet neighbour list and sutton-chen potentials
!
! Note: Designed to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov january 2016
! contrib   - r.davidchak (eeam) june 2012
! contrib   - b.palmer (2band) may 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( In    ) :: elrcm,vlrcm
  Real( Kind = wp ),                        Intent(   Out ) :: engden,virden
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Type( comms_type ),                       Intent( InOut ) :: comm

  Logical           :: safe = .true.
  Integer           :: fail,limit,i,j,k,l,k0
  Real( Kind = wp ) :: rhosqr,rdr,rrr,ppp,fk0,fk1,fk2,t1,t2

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

  Character( Len = 256 ) :: message

! check on mixing metal types done in read_field

! initialise energy and virial accumulators

  engden=0.0_wp
  virden=0.0_wp

! initialise density array

  rho=0.0_wp
  If (l2bmet) rhs = 0.0_wp

! All calls below act on rho (rhs)

! calculate local atomic density
! outer loop over atoms

  fail=0
  Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'metal_ld_compute allocation failure'
     Call error(0,message)
  End If

  Do i=1,natms
     limit=list(0,i) ! Get list limit

! calculate interatomic distances

     Do k=1,limit
        j=list(k,i)

        xxt(k)=xxx(i)-xxx(j)
        yyt(k)=yyy(i)-yyy(j)
        zzt(k)=zzz(i)-zzz(j)
     End Do

! periodic boundary conditions not needed by LC construction
!
!     Call images(imcon,cell,limit,xxt,yyt,zzt)

! square of distances

     Do k=1,limit
        rrt(k) = Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
     End Do

! calculate contributions to local density

     If (tabmet > 0) Then         ! EAM contributions
        Call metal_ld_collect_eam(i,rrt,safe)
     Else ! If (tabmet == 0) Then ! FST contributions
        Call metal_ld_collect_fst(i,rmet,rrt,safe)
     End If
  End Do

  Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a,i0)') 'metal_ld_compute allocation failure'
     Call error(0,message)
  End If

! Check safety for densities of EAM and MBPC

  Call gcheck(comm,safe)
  If (.not.safe) Call error(506)

  Do i=1,natms

! calculate density terms to energy and virial

     If (tabmet > 0) Then ! EAM potential

! potential function index

        k0=ltype(i)

! Now start traditional s-band (EAM & EEAM) or d-band for 2B(EAM & EEAM)

! validity of potential

        If (Abs(fmet(1,k0,1)) > zero_plus) Then

! check for unsafe densities (mind start was shifted)

           If (.not.ls_met) Then ! fmet over rho grid
              rhosqr = rho(i)
           Else                  ! fmet over Sqrt(rho) grid
              rhosqr = Sqrt(rho(i))
           End If
           If (rhosqr >= fmet(2,k0,1)+5.0_wp*fmet(4,k0,1)) Then
              If (rhosqr <= fmet(3,k0,1)) Then

! interpolation parameters

                 rdr = 1.0_wp/fmet(4,k0,1)
                 rrr = rhosqr - fmet(2,k0,1)
                 l   = Min(Nint(rrr*rdr),Nint(fmet(1,k0,1))-1)
                 If (l < 5) Then ! catch unsafe value
                    Write(*,*) 'good density range problem: (LTG,RHO) ',ltg(i),rho(i)
                    safe=.false.
                    l=6
                 End If
                 ppp = rrr*rdr - Real(l,wp)

! calculate embedding energy using 3-point interpolation

                 fk0 = fmet(l-1,k0,1)
                 fk1 = fmet(l  ,k0,1)
                 fk2 = fmet(l+1,k0,1)

                 t1 = fk1 + ppp*(fk1 - fk0)
                 t2 = fk1 + ppp*(fk2 - fk1)

                 If (ppp < 0.0_wp) Then
                    engden = engden + t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                 Else If (l == 5) Then
                    engden = engden + t2
                 Else
                    engden = engden + t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                 End If

! calculate derivative of embedding function wrt density
! using 3-point interpolation and STORE/OVERWRITE result in rho array

                 fk0 = fmet(l-1,k0,2)
                 fk1 = fmet(l  ,k0,2)
                 fk2 = fmet(l+1,k0,2)

                 t1 = fk1 + ppp*(fk1 - fk0)
                 t2 = fk1 + ppp*(fk2 - fk1)

                 If (ppp < 0.0_wp) Then
                    If (.not.ls_met) Then ! fmet over rho grid
                       rho(i) = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                    Else                  ! fmet over Sqrt(rho) grid
                       rho(i) = 0.5_wp*(t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp))/rhosqr
                    End If
                 Else If (l == 5) Then
                    If (.not.ls_met) Then ! fmet over rho grid
                       rho(i) = t2
                    Else                  ! fmet over Sqrt(rho) grid
                       rho(i) = 0.5_wp*t2/rhosqr
                    End If
                 Else
                    If (.not.ls_met) Then ! fmet over rho grid
                       rho(i) = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                    Else                  ! fmet over Sqrt(rho) grid
                       rho(i) = 0.5_wp*(t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp))/rhosqr
                    End If
                 End If

              Else ! RLD: assume that fmet(rho(i) > fmet(3,k0,1)) = fmet(rho(i) = fmet(3,k0,1))

                l      = Nint(fmet(1,k0,1))

                engden = engden + fmet(l,k0,1)

                rho(i) = 0.0_wp

              End If
           Else
              Write(*,*) 'bad density range problem: (LTG,RHO) ',ltg(i),rho(i)
              safe=.false.
           End If

        End If

! Atomic density (rho & rhs) are overwritten here in order
! to be reused as embedding derivative(s) helper array(s)
! i.e. hold d_fmet/d_rho for later usage in metal_forces

! Now if we have 2B(EAM & EEAM) then do s-band too

        If (l2bmet) Then

! validity of potential

           If (Abs(fmes(1,k0,1)) > zero_plus) Then

! check for unsafe densities (mind start was shifted)

              If (.not.ls_met) Then ! fmes over rhs grid
                 rhosqr = rhs(i)
              Else                  ! fmes over Sqrt(rhs) grid
                 rhosqr = Sqrt(rhs(i))
              End If
              If (rhosqr >= fmes(2,k0,1)+5.0_wp*fmes(4,k0,1)) Then
                 If (rhosqr <= fmes(3,k0,1)) Then

! interpolation parameters

                    rdr = 1.0_wp/fmes(4,k0,1)
                    rrr = rhosqr - fmes(2,k0,1)
                    l   = Min(Nint(rrr*rdr),Nint(fmes(1,k0,1))-1)
                    If (l < 5) Then ! catch unsafe value
                       Write(*,*) 'good density range problem: (LTG,RHS) ',ltg(i),rhs(i)
                       safe=.false.
                       l=6
                    End If
                    ppp = rrr*rdr - Real(l,wp)

! calculate embedding energy using 3-point interpolation

                    fk0 = fmes(l-1,k0,1)
                    fk1 = fmes(l  ,k0,1)
                    fk2 = fmes(l+1,k0,1)

                    t1 = fk1 + ppp*(fk1 - fk0)
                    t2 = fk1 + ppp*(fk2 - fk1)

                    If (ppp < 0.0_wp) Then
                       engden = engden + t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                    Else If (l == 5) Then
                       engden = engden + t2
                    Else
                       engden = engden + t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                    End If

! calculate derivative of embedding function wrt density
! using 3-point interpolation and STORE/OVERWRITE result in rhs array

                    fk0 = fmes(l-1,k0,2)
                    fk1 = fmes(l  ,k0,2)
                    fk2 = fmes(l+1,k0,2)

                    t1 = fk1 + ppp*(fk1 - fk0)
                    t2 = fk1 + ppp*(fk2 - fk1)

                    If (ppp < 0.0_wp) Then
                       If (.not.ls_met) Then ! fmes over rhs grid
                          rhs(i) = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                       Else                  ! fmes over Sqrt(rhs) grid
                          rhs(i) = 0.5_wp*(t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp))/rhosqr
                       End If
                    Else If (l == 5) Then
                       If (.not.ls_met) Then ! fmes over rhs grid
                          rhs(i) = t2
                       Else                  ! fmes over Sqrt(rhs) grid
                          rhs(i) = 0.5_wp*t2/rhosqr
                       End If
                    Else
                       If (.not.ls_met) Then ! fmes over rhs grid
                          rhs(i) = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                       Else                  ! fmes over Sqrt(rhs) grid
                          rhs(i) = 0.5_wp*(t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp))/rhosqr
                       End If
                    End If

                 Else ! RLD: assume that fmes(rhs(i) > fmes(3,k0,1)) = fmes(rhs(i) = fmes(3,k0,1))

                   l      = Nint(fmes(1,k0,1))

                   engden = engden + fmes(l,k0,1)

                   rhs(i) = 0.0_wp

                 End If
              Else
                 Write(message,'(a,2(i0,1x))') 'bad density range problem: (LTG,RHS) ',ltg(i),rhs(i)
                 Call info(message)
                 safe=.false.
              End If

           End If

        End If

     Else ! If (tabmet == 0) Then FST of metal potentials

        If      (rho(i) > zero_plus) Then

! calculate analytical square root of (density + lrc to it)

           rhosqr = Sqrt(rho(i) + elrcm(ltype(i)))
           engden = engden - rhosqr
           virden = virden + vlrcm(ltype(i))/rhosqr

! store the derivatives of the FST embedding-like function
! (with corrected density) in rho array

           rho(i) = 0.5_wp/rhosqr

        Else If (rho(i) < -zero_plus) Then

! check for unsafe densities (rho was initialised to zero)

           safe=.false.

        End If

     End If

  End Do

! Check safety for densities

  Call gcheck(comm,safe)
  If (.not.safe) Call error(507)

! virial term (averaged per node)

  Call gsum(comm,virden)
  virden=virden/Real(comm%mxnode,wp)

! calculate stress tensor (density contributions are to
! diagonal elements only)

  stress(1)=stress(1)-virden/3.0_wp
  stress(5)=stress(5)-virden/3.0_wp
  stress(9)=stress(9)-virden/3.0_wp

! obtain atomic densities for outer border regions

  Call metal_ld_set_halo(comm)

End Subroutine metal_ld_compute

Subroutine metal_lrc(rmet,elrcm,vlrcm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to evaluate metal long-range corrections to
! pressure and energy in a 3D periodic system
!
! copyright - daresbury laboratory
! author    - w.smith june 1995
! amended   - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent(   Out ) :: elrcm,vlrcm
  Type( comms_type ), Intent( InOut ) :: comm

  Logical, Save     :: newjob = .true.
  Integer           :: i,j,k0,k1,k2,kmet,keypot,nnn,mmm

  Real( Kind = wp ) :: elrc0,elrc1,elrc2,elrcsum,vlrc0,vlrc1,vlrc2, tmp, &
                       eps,sig,nnnr,mmmr,ccc, aaa,rr0,ppp,zet,qqq,eee


! long-range corrections to energy, pressure and density

  elrcm   = 0.0_wp
  vlrcm   = 0.0_wp
  elrcsum = 0.0_wp

  If (imcon /= 0 .and. imcon /= 6) Then
     kmet = 0

     Do i=1,ntpatm
        Do j=1,i

           elrc0=0.0_wp
           elrc1=0.0_wp
           elrc2=0.0_wp

           vlrc0=0.0_wp
           vlrc1=0.0_wp
           vlrc2=0.0_wp

           kmet = kmet + 1
           k0 = lstmet(kmet)

           keypot=ltpmet(k0)
           If      (keypot == 3) Then

! sutton-chen potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              nnn=Nint(prmmet(3,k0)) ; nnnr=Real(nnn,wp)
              mmm=Nint(prmmet(4,k0)) ; mmmr=Real(mmm,wp)
              ccc=prmmet(5,k0)

              elrc0=eps*sig**3*(sig/rmet)**(nnn-3)/(nnnr-3.0_wp)
              vlrc0=nnnr*elrc0

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

              If (i /= j) Then
                 elrc0 = elrc0*2.0_wp
                 vlrc0 = vlrc0*2.0_wp
              End If

              elrcm(0) = elrcm(0) + twopi*volm*dens(i)*dens(j)*elrc0
              vlrcm(0) = vlrcm(0) - twopi*volm*dens(i)*dens(j)*vlrc0

              tmp=sig**3*(sig/rmet)**(mmm-3)/(mmmr-3.0_wp)
              If (i == j) Then
                 elrc1=tmp*(eps*ccc)**2
                 elrcm(i)=elrcm(i)+fourpi*dens(i)*elrc1
                 elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1

                 vlrc1=mmmr*elrc1
                 vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1
              Else
                 k1=lstmet((i*(i+1))/2)
                 k2=lstmet((j*(j+1))/2)

                 elrc1=tmp*(prmmet(1,k1)*prmmet(5,k1))**2
                 elrc2=tmp*(prmmet(1,k2)*prmmet(5,k2))**2
                 elrcm(i)=elrcm(i)+fourpi*dens(j)*elrc1
                 elrcm(j)=elrcm(j)+fourpi*dens(i)*elrc2
                 elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)

                 vlrc1=mmmr*elrc1
                 vlrc2=mmmr*elrc2
                 vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                 vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2
              End If

           Else If (keypot == 4) Then

! gupta potentials

              aaa=prmmet(1,k0)
              rr0=prmmet(2,k0)
              ppp=prmmet(3,k0)
              zet=prmmet(4,k0)
              qqq=prmmet(5,k0)
              eee=Exp(-ppp*(rmet-rr0)/rr0)

              elrc0=2.0_wp*aaa*(rr0/ppp)*(rmet**2+2.0_wp*rmet*(rr0/ppp)+2.0_wp*(rr0/ppp)**2)*eee
              vlrc0=2.0_wp*aaa*rmet**3*eee+3.0_wp*elrc0

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

              If (i /= j) Then
                 elrc0=elrc0*2.0_wp
                 vlrc0=vlrc0*2.0_wp
              End If

              elrcm(0)=elrcm(0)+twopi*volm*dens(i)*dens(j)*elrc0
              vlrcm(0)=vlrcm(0)-twopi*volm*dens(i)*dens(j)*vlrc0

              eee=Exp(-2.0_wp*qqq*(rmet-rr0)/rr0)

              If (i == j) Then
                 elrc1=(rmet**2+2.0_wp*rmet*(0.5_wp*rr0/qqq)+2.0_wp*(0.5_wp*rr0/qqq)**2)*(0.5_wp*rr0/qqq)*eee*zet**2
                 elrcm(i)=elrcm(i)+fourpi*dens(i)*elrc1
                 elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1

                 vlrc1=(rmet**3+3.0_wp*rmet**2*(0.5_wp*rr0/qqq)+6.0_wp*rmet*(0.5_wp*rr0/qqq)**2+(0.5_wp*rr0/qqq)**3)*eee*zet**2
                 vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1
              Else
                 elrc1=(rmet**2+2.0_wp*rmet*(0.5_wp*rr0/qqq)+2.0_wp*(0.5_wp*rr0/qqq)**2)*(0.5_wp*rr0/qqq)*eee*zet**2
                 elrc2=elrc2
                 elrcm(i)=elrcm(i)+fourpi*dens(j)*elrc1
                 elrcm(j)=elrcm(j)+fourpi*dens(i)*elrc2
                 elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)

                 vlrc1=(rmet**3+3.0_wp*rmet**2*(0.5_wp*rr0/qqq)+6.0_wp*rmet*(0.5_wp*rr0/qqq)**2+(0.5_wp*rr0/qqq)**3)*eee*zet**2
                 vlrc2=vlrc1
                 vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                 vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2
              End If

           Else If (keypot == 5) Then

! many-body perturbation component only potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              mmm=Nint(prmmet(3,k0)) ; mmmr=Real(mmm,wp)

! No pairwise contributions for mbpc potentials!!!

!              elrc0=0.0_wp
!              vlrc0=0.0_wp

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

!              If (i /= j) Then
!                 elrc0 = elrc0*2.0_wp
!                 vlrc0 = vlrc0*2.0_wp
!              End If

!              elrcm(0) = elrcm(0) + twopi*volm*dens(i)*dens(j)*elrc0
!              vlrcm(0) = vlrcm(0) - twopi*volm*dens(i)*dens(j)*vlrc0

              tmp=sig/((mmmr-3.0_wp)*rmet**(mmm-3))
              If (i == j) Then
                 elrc1=tmp*eps**2
                 elrcm(i)=elrcm(i)+fourpi*dens(i)*elrc1
                 elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1

                 vlrc1=mmmr*elrc1
                 vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1
              Else
                 k1=lstmet((i*(i+1))/2)
                 k2=lstmet((j*(j+1))/2)

                 elrc1=tmp*prmmet(1,k1)**2
                 elrc2=tmp*prmmet(1,k2)**2
                 elrcm(i)=elrcm(i)+fourpi*dens(j)*elrc1
                 elrcm(j)=elrcm(j)+fourpi*dens(i)*elrc2
                 elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)

                 vlrc1=mmmr*elrc1
                 vlrc2=mmmr*elrc2
                 vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                 vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2
              End If

           End If

        End Do
     End Do
  End If

  If (newjob) Then
     newjob =.false.

     If (comm%idnode == 0) Then
        Write(nrite,"(1p,                                  &
        & 'long-range correction to metal energy    ',e15.6,/,1x, &
        & 'lr correction for metal atom density     ',e15.6,/,1x, &
        & '1st partial lr correction to metal virial',e15.6,/)")  &
           elrcm(0)/engunit,elrcsum/engunit**2,vlrcm(0)/engunit

        Write(nrite,"(1x,'density dependent energy and virial corrections',/)")
        Do i=1,ntpatm
           kmet=lstmet((i*(i+1))/2)
           If (lstmet(kmet) > 0) Write(nrite,"(25x,a8,1p,2e15.6)") &
              unqatm(i),elrcm(i)/engunit,vlrcm(i)/engunit
        End Do
     End If
  End If

End Subroutine metal_lrc

Subroutine metal_table_read(l_top,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABEAM file (for metal EAM & EEAM forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 2006
! amended   - i.t.todorov march 2016
! contrib   - r.davidchak (eeam) june 2012
! contrib   - b.palmer (2band) may 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical, Intent( In    ) :: l_top
  Type( comms_type ), Intent ( InOut ) :: comm

  Logical                :: safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 4   ) :: keyword
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail(1:2),i,j,ipot,numpot,ktype,ngrid, &
                            cp,cd,cds,ce,ces,katom1,katom2,keymet,k0,jtpatm
  Real( Kind = wp )      :: start,finish

  Integer,           Dimension( : ), Allocatable :: cpair, cdens,cdnss, &
                                                    cembed,cembds
  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  Character( Len = 256 ) :: message

  fail=0
  If      (tabmet == 1) Then ! EAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet),                              &
                                              cembed(1:ntpmet),                  Stat=fail(1))
  Else If (tabmet == 2) Then ! EEAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet**2),                           &
                                              cembed(1:ntpmet),                  Stat=fail(1))
  Else If (tabmet == 3) Then ! 2BEAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet),cdnss(1:ntpmet*(ntpmet+1)/2), &
                                              cembed(1:ntpmet),cembds(1:ntpmet), Stat=fail(1))
  Else If (tabmet == 4) Then ! 2BEEAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet**2),cdnss(1:ntpmet**2), &
                                              cembed(1:ntpmet),cembds(1:ntpmet), Stat=fail(1))
  End If
  Allocate (buffer(1:mxgmet),                                                    Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'metal_table_read allocation failure'
     Call error(0,message)
  End If
  cpair=0 ; cp=0
  cdens=0 ; cd=0
  cembed=0 ; ce=0
  If (tabmet == 3 .or. tabmet == 4) Then
    cdnss=0 ; cds=0
    cembds=0 ; ces=0
  End If

  If (comm%idnode == 0) Open(Unit=ntable, File='TABEAM')

! skip header record

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

! read number of potential functions in file

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100
  Call get_word(record,word)
  numpot = Nint(word_2_real(word))

  Do ipot=1,numpot

! read data type, atom labels, number of points, start and end

     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

! identify data type

     Call get_word(record,keyword)
     Call lower_case(keyword)
     If      (keyword == 'pair') Then
          ktype = 1
     Else If (keyword == 'dens' .or. keyword == 'dden') Then
          ktype = 2
     Else If (keyword == 'embe' .or. keyword == 'demb') Then
          ktype = 3
     Else If (keyword == 'sden') Then
          ktype = 4
     Else If (keyword == 'semb') Then
          ktype = 5
     Else
          Call error(151)
     End If

! identify atom types

     Call get_word(record,atom1)
     If (ktype == 1 .or.                                        & ! pair
         (ktype == 2 .and. (tabmet == 2 .or. tabmet == 4)) .or. & ! den for EEAM and dden for 2BEEAM
         (ktype == 4 .and. (tabmet == 3 .or. tabmet == 4))) Then  ! sden for 2B(EAM and EEAM)
        Call get_word(record,atom2)
     Else
        atom2 = atom1
     End If

! data specifiers

     Call get_word(record,word)
     ngrid = Nint(word_2_real(word))
     Call get_word(record,word)
     start  = word_2_real(word)
     Call get_word(record,word)
     finish = word_2_real(word)

! check atom identities

     katom1=0
     katom2=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0) Then
        If (l_top) &
           Write(message,'(a)') '****',atom1,'***',atom2,'**** entry in TABEAM'
        Call error(81,message,.true.)
     End If

! store working parameters

     buffer(1)=Real(ngrid+4,wp) ! as if there are 4 extra elements after finish
     buffer(4)=(finish-start)/Real(ngrid-1,wp)
     buffer(2)=start-5.0_wp*buffer(4)
     buffer(3)=finish

     If (l_top) Then
        Write(message,"(1x,i10,4x,2a8,3x,2a4,2x,i6,1p,3e15.6)") &
        ipot,atom1,atom2,'EAM-',keyword,ngrid,start,finish,buffer(4)
      Call info(message,.true.)
    End IF

! check array dimensions

     If (ngrid+4 > mxgmet) Then
        Call warning(270,Real(ngrid+4,wp),Real(mxgmet,wp),0.0_wp)
        Call error(48)
     End If

     keymet=(Max(katom1,katom2)*(Max(katom1,katom2)-1))/2 + Min(katom1,katom2)
     k0=lstmet(keymet)

! check for undefined potential

     If (k0 == 0) Call error(508)

! read in potential arrays

     Do i=1,(ngrid+3)/4
        j=Min(4,ngrid-(i-1)*4)
        If (comm%idnode == 0) Then
           Read(Unit=ntable, Fmt=*, End=100) buffer(4*i+1:4*i+j)
        Else
           buffer(4*i+1:4*i+j)=0.0_wp
        End If
     End Do

     Call gsum(comm,buffer(5:ngrid+4))

! copy data to internal arrays

     If       (ktype == 1) Then

! pair potential terms

! Set indices

!        k0=lstmet(keymet)

        cp=cp+1
        If (Any(cpair(1:cp-1) == k0)) Then
           Call error(509)
        Else
           cpair(cp)=k0
        End If

        vmet(1,k0,1)=buffer(1)
        vmet(2,k0,1)=buffer(2)
        vmet(3,k0,1)=buffer(3)
        vmet(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             vmet(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             vmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of pair potential function

        Call metal_table_derivatives(k0,buffer,Size(vmet,2),vmet)

! adapt derivatives for use in interpolation

        Do i=5,ngrid+4
           vmet(i,k0,2)=-(Real(i,wp)*buffer(4)+buffer(2))*vmet(i,k0,2)
        End Do

     Else If (ktype == 2) Then

! density function terms
! s-density density function terms for EAM & EEAM
! d-density density function terms for 2B(EAM & EEAM)

! Set indices

        If      (tabmet == 1 .or. tabmet == 3) Then ! EAM
           k0=katom1
        Else If (tabmet == 2 .or. tabmet == 4) Then ! EEAM
           k0=(katom1-1)*ntpatm+katom2
        End If

        cd=cd+1
        If (Any(cdens(1:cd-1) == k0)) Then
           Call error(510)
        Else
           cdens(cd)=k0
        End If

        dmet(1,k0,1)=buffer(1)
        dmet(2,k0,1)=buffer(2)
        dmet(3,k0,1)=buffer(3)
        dmet(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             dmet(i,k0,1)=0.0_wp
           Else
             dmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of density function

        Call metal_table_derivatives(k0,buffer,Size(dmet,2),dmet)

! adapt derivatives for use in interpolation

        dmet(1,k0,2)=0.0_wp
        dmet(2,k0,2)=0.0_wp
        dmet(3,k0,2)=0.0_wp
        dmet(4,k0,2)=0.0_wp

        Do i=5,ngrid+4
           dmet(i,k0,2)=-(Real(i,wp)*buffer(4)+buffer(2))*dmet(i,k0,2)
        End Do

     Else If (ktype == 3) Then

! embedding function terms
! s-density embedding function terms for EAM & EEAM
! d-density embedding function terms for 2B(EAM & EEAM)

! Set indices

        k0=katom1

        ce=ce+1
        If (Any(cembed(1:ce-1) == k0)) Then
           Call error(511)
        Else
           cembed(ce)=k0
        End If

        fmet(1,k0,1)=buffer(1)
        fmet(2,k0,1)=buffer(2)
        fmet(3,k0,1)=buffer(3)
        fmet(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             fmet(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             fmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of embedding function

        Call metal_table_derivatives(k0,buffer,Size(fmet,2),fmet)

     Else If (ktype == 4) Then

! s-density function terms

! The 2BM formalism for alloys allows for a mixed s-band density: rho_{atom1,atom2} /= 0
! (and in general for the EEAM it may be non-symmetric: rho_{atom1,atom2} may be /= rho_{atom2,atom2})
! Some 2BM models rho_{atom1,atom1}=rho_{atom2,atom2}==0 with rho_{atom1,atom2} /= 0
! whereas others choose not to have mixed s-band densities.

! Set indices

        If (tabmet == 3) Then ! 2BMEAM
!           k0=lstmet(keymet)
        Else If (tabmet == 4) Then ! 2BMEEAM
           k0=(katom1-1)*ntpatm+katom2
        End If

        cds=cds+1
        If (Any(cdnss(1:cds-1) == k0)) Then
           Call error(510)
        Else
           cdnss(cds)=k0
        End If

        dmes(1,k0,1)=buffer(1)
        dmes(2,k0,1)=buffer(2)
        dmes(3,k0,1)=buffer(3)
        dmes(4,k0,1)=buffer(4)

        If (Nint(buffer(1)) > 5) Then

           Do i=5,mxgmet
              If (i-4 > ngrid) Then
                 dmes(i,k0,1)=0.0_wp
              Else
                 dmes(i,k0,1)=buffer(i)
              End If
           End Do
! calculate derivative of density function

           Call metal_table_derivatives(k0,buffer,Size(dmes,2),dmes)

! adapt derivatives for use in interpolation

           dmes(1,k0,2)=0.0_wp
           dmes(2,k0,2)=0.0_wp
           dmes(3,k0,2)=0.0_wp
           dmes(4,k0,2)=0.0_wp

           Do i=5,ngrid+4
              dmes(i,k0,2)=-(Real(i,wp)*buffer(4)+buffer(2))*dmes(i,k0,2)
           End Do

        End If

     Else If (ktype == 5) Then

! s-embedding function terms

! Set index

        k0=katom1

        ces=ces+1
        If (Any(cembds(1:ces-1) == k0)) Then
           Call error(511)
        Else
           cembds(ces)=k0
        End If

        fmes(1,k0,1)=buffer(1)
        fmes(2,k0,1)=buffer(2)
        fmes(3,k0,1)=buffer(3)
        fmes(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             fmes(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             fmes(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of embedding function

        Call metal_table_derivatives(k0,buffer,Size(fmes,2),fmes)

     End If

  End Do

  If (comm%idnode == 0) Close(Unit=ntable)
  If (l_top) Then
    Write(message,'(a)') 'potential tables read from TABEAM file'
    Call info(message,.true.)
  End IF

  If      (tabmet == 1 .or. tabmet == 2) Then ! EAM & EEAM
     Deallocate (cpair,cdens,cembed,              Stat=fail(1))
  Else If (tabmet == 3 .or. tabmet == 4) Then ! 2B(EAM & EEAM)
     Deallocate (cpair,cdens,cdnss,cembed,cembds, Stat=fail(1))
  End If
  Deallocate (buffer,                             Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a,i0)') 'metal_table_read deallocation failure'
     Call error(0,message)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine metal_table_read


Subroutine metal_generate(rmet)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating potential energy and force arrays
! for metal potentials
!
! copyright - daresbury laboratory
! author    - w.smith june 2006
! amended   - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  Real( Kind = wp ), Intent( In    ) :: rmet

  Integer           :: i,imet,kmet,keypot,katom1,katom2, &
                       pmet,qmet,nnn,mmm
  Real( Kind = wp ) :: dlrpot,rrr,              &
                       eps,sig,nnnr,mmmr,       &
                       cc0,cc1,cc2,cc3,cc4,     &
                       aaa,bbb,ccc,ddd,ppp,qqq, &
                       bet,cut1,cut2,rr0

! define grid resolution for potential arrays

  dlrpot=rmet/Real(mxgmet-1,wp)

! construct arrays for metal potentials

  kmet=0
  Do katom1=1,ntpatm
     Do katom2=1,katom1
        kmet=kmet+1

! calculate potentials for defined interactions

        imet=lstmet(kmet)
        keypot=ltpmet(imet)
        If (keypot > 0) Then

! store array specification parameters

           vmet(1,imet,1)=Real(mxgmet,wp)
           vmet(2,imet,1)=0.0_wp          ! l_int(min) >= 1
           vmet(3,imet,1)=rmet            ! rmet=rcut
           vmet(4,imet,1)=dlrpot

           Do i=1,4
              vmet(i,imet,2)=vmet(i,imet,1)
              dmet(i,imet,1)=vmet(i,imet,1)
              dmet(i,imet,2)=0.0_wp
           End Do

           If      (keypot == 1) Then

! finnis-sinclair potentials

              cc0=prmmet(1,imet)
              cc1=prmmet(2,imet)
              cc2=prmmet(3,imet)
              ccc=prmmet(4,imet)
              ddd=prmmet(6,imet)
              bet=prmmet(7,imet)
              cut1=ccc+4.0_wp*dlrpot
              cut2=ddd+4.0_wp*dlrpot

              vmet(3,imet,1:2)=cut1
              dmet(3,imet,1)=cut2

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot

                 If (rrr <= cut1) Then
                    vmet(i,imet,1)=(cc0+cc1*rrr+cc2*rrr**2)*(rrr-ccc)**2
                    vmet(i,imet,2)=-rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2) * &
                                    (rrr-ccc)+(cc1+2.0_wp*cc2*rrr)*(rrr-ccc)**2)
                 End If

                 If (rrr <= cut2) Then
                    dmet(i,imet,1)=(rrr-ddd)**2+bet*(rrr-ddd)**3/ddd
                    dmet(i,imet,2)=-rrr*(2.0_wp*(rrr-ddd)+3.0_wp*bet*(rrr-ddd)**2/ddd)
                 End If
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=prmmet(5,imet)**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=prmmet(5,pmet)**2
                 dmet(2,imet,2)=prmmet(5,qmet)**2
              End If

           Else If (keypot == 2) Then

! extended finnis-sinclair potentials

              cc0=prmmet(1,imet)
              cc1=prmmet(2,imet)
              cc2=prmmet(3,imet)
              cc3=prmmet(4,imet)
              cc4=prmmet(5,imet)
              ccc=prmmet(6,imet)
              ddd=prmmet(8,imet)
              bbb=prmmet(9,imet)
              cut1=ccc+4.0_wp*dlrpot
              cut2=ddd+4.0_wp*dlrpot

              vmet(3,imet,1:2)=cut1
              dmet(3,imet,1)=cut2

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot

                 If (rrr <= cut1) Then
                    vmet(i,imet,1)=(cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc)**2
                    vmet(i,imet,2)=-rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc) + &
                                         (cc1+2.0_wp*cc2*rrr+3.0_wp*cc3*rrr**2+4.0_wp*cc4*rrr**3)*(rrr-ccc)**2)
                 End If

                 If (rrr <= cut2) Then
                    dmet(i,imet,1)=(rrr-ddd)**2+bbb**2*(rrr-ddd)**4
                    dmet(i,imet,2)=-rrr*(2.0_wp*(rrr-ddd)+4.0_wp*bbb**2*(rrr-ddd)**3)
                 End If
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=prmmet(7,imet)**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=prmmet(7,pmet)**2
                 dmet(2,imet,2)=prmmet(7,qmet)**2
              End If

           Else If (keypot == 3) Then

! sutton-chen potentials

              eps=prmmet(1,imet)
              sig=prmmet(2,imet)
              nnn=Nint(prmmet(3,imet)) ; nnnr=Real(nnn,wp)
              mmm=Nint(prmmet(4,imet)) ; mmmr=Real(mmm,wp)

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot
                 vmet(i,imet,1)=eps*(sig/rrr)**nnn
                 vmet(i,imet,2)=nnnr*eps*(sig/rrr)**nnn
                 dmet(i,imet,1)=(sig/rrr)**mmm
                 dmet(i,imet,2)=mmmr*(sig/rrr)**mmm
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=(prmmet(1,imet)*prmmet(5,imet))**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=(prmmet(1,pmet)*prmmet(5,pmet))**2
                 dmet(2,imet,2)=(prmmet(1,qmet)*prmmet(5,qmet))**2
              End If

           Else If (keypot == 4) Then

! gupta potentials

              aaa=prmmet(1,imet)
              rr0=prmmet(2,imet)
              ppp=prmmet(3,imet)
              qqq=prmmet(5,imet)

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot

                 cut1=(rrr-rr0)/rr0
                 cut2=cut1+1.0_wp

                 vmet(i,imet,1)=2.0_wp*aaa*Exp(-ppp*cut1)
                 vmet(i,imet,2)=vmet(i,imet,1)*ppp*cut2
                 dmet(i,imet,1)=Exp(-2.0_wp*qqq*cut1)
                 dmet(i,imet,2)=2.0_wp*dmet(i,imet,1)*qqq*cut2
              End Do

              dmet(1,imet,2)=prmmet(4,imet)**2
              dmet(2,imet,2)=dmet(1,imet,2)

           Else If (keypot == 5) Then

! many-body perturbation component only potentials

              eps=prmmet(1,imet)
              sig=prmmet(2,imet)
              mmm=Nint(prmmet(3,imet)) ; mmmr=Real(mmm,wp)

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot
!                 vmet(i,imet,1)=0.0_wp
!                 vmet(i,imet,2)=0.0_wp
                 nnnr=sig/rrr**mmm
                 dmet(i,imet,1)=nnnr*merf(i)
                 dmet(i,imet,2)=mmmr*dmet(i,imet,1)-rrr*nnnr*mfer(i)
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=prmmet(1,imet)**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=prmmet(1,pmet)**2
                 dmet(2,imet,2)=prmmet(1,qmet)**2
              End If

           Else

              Call error(151)

           End If

        End If
     End Do
  End Do

End Subroutine metal_generate

Subroutine metal_generate_erf(rmet)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating erf and fer arrays for
! many-body perturbation component only potentials
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real( Kind = wp ), Intent( In    ) :: rmet

  Integer           :: imet
  Real( Kind = wp ) :: alpha,beta

! Determine alpha and beta for the erf bit of all MBPC potentials

  If (Any(ltpmet(1:ntpmet) == 5)) Then ! all are == 5 == MBPC
     alpha= 0.0_wp
     beta = 0.0_wp
     Do imet=1,ntpmet
        alpha=Max(alpha,Abs(prmmet(6,imet)))
        beta =Max(beta, Abs(prmmet(7,imet)))
     End Do

! If unset then set to defaults

     If (alpha <= zero_plus) alpha=20.0_wp
     If (beta  <= zero_plus) beta =Min(1.5_wp,0.2_wp*rmet)

! Allocate arrays: merf,mfer

     Call allocate_metal_erf_arrays()

! Generate error function and derivative arrays

     Call erfgen_met(rmet,alpha,beta,mxgmet,merf,mfer)

! Translate merf and mfer to the functional form 0.5*{1+erf[alpha(r-beta)]}

     merf(5:mxgmet)=0.5_wp*(1.0_wp+merf(5:mxgmet))
     mfer(5:mxgmet)=0.5_wp*mfer(5:mxgmet)
  End If

End Subroutine metal_generate_erf

Subroutine metal_ld_collect_eam(iatm,rrt,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local atomic density for
! Embedded Atom Model & Extended Embedded Atom Model metal potentials
!
! Note: Designed to be used as part of metal_ld_compute
!
! copyright - daresbury laboratory
! author    - w.smith june 1995
! amended   - i.t.todorov december 2016
! contrib   - r.davidchak (eeam) june 2012
! contrib   - b.palmer (2band) may 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: rrt
  Logical,                                  Intent( InOut ) :: safe

  Integer           :: m,ai,ki,jatm,aj,kj,l,key,k0
  Real( Kind = wp ) :: rrr,rdr,rr1,ppp,vk0,vk1,vk2,t1,t2,density

! global type of itam

  ai=ltype(iatm)

! start of primary loop for density

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

! interatomic distance

     rrr=rrt(m)

! Now start traditional s-band (EAM & EEAM) or d-band for 2B(EAM & EEAM)

! first metal atom density and validity and truncation of potential

     If (Abs(dmet(1,kj,1)) > zero_plus .and. Nint(dmet(1,kj,1)) > 5) Then
        If (rrr <= dmet(3,kj,1)) Then

! interpolation parameters

           rdr = 1.0_wp/dmet(4,kj,1)
           rr1 = rrr - dmet(2,kj,1)
           l   = Min(Nint(rr1*rdr),Nint(dmet(1,kj,1))-1)
           If (l < 5) Then ! catch unsafe value
              safe=.false.
              Write(*,*) 'aaa',l,iatm,jatm,rrr
              l=6
           End If
           ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

           vk0 = dmet(l-1,kj,1)
           vk1 = dmet(l  ,kj,1)
           vk2 = dmet(l+1,kj,1)

           t1 = vk1 + ppp*(vk1 - vk0)
           t2 = vk1 + ppp*(vk2 - vk1)

           If (ppp < 0.0_wp) Then ! density is a positive function!
              density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
           Else If (l == 5) Then
              density = t2
           Else
              density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
           End If

           rho(iatm) = rho(iatm) + density
           If (ki == kj .and. jatm <= natms) rho(jatm) = rho(jatm) + density

        End If
     End If

! second metal atom density and validity and truncation of potential

     If (Abs(dmet(1,ki,1)) > zero_plus .and. Nint(dmet(1,ki,1)) > 5) Then
        If (ki /= kj .and. jatm <= natms) Then
           If (rrr <= dmet(3,ki,1)) Then

! interpolation parameters

              rdr = 1.0_wp/dmet(4,ki,1)
              rr1 = rrr - dmet(2,ki,1)
              l   = Min(Nint(rr1*rdr),Nint(dmet(1,ki,1))-1)
              If (l < 5) Then ! catch unsafe value
                 safe=.false.
                 Write(*,*) 'bbb',l,iatm,jatm,rrr
                 l=6
              End If
              ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

              vk0 = dmet(l-1,ki,1)
              vk1 = dmet(l  ,ki,1)
              vk2 = dmet(l+1,ki,1)

              t1 = vk1 + ppp*(vk1 - vk0)
              t2 = vk1 + ppp*(vk2 - vk1)

              If (ppp < 0.0_wp) Then ! density is a positive function!
                 density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                 If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
              Else If (l == 5) Then
                 density = t2
              Else
                 density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                 If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
              End If

              rho(jatm) = rho(jatm) + density

           End If
        End If
     End If

! Now if we have the 2B(EAM & EEAM) then do the s-band (dmes and rhs are defined)

     If (l2bmet) Then
        If      (tabmet == 3) Then ! 2BEAM

! 2BEAM has symmetric s-densities with respect to atom type
! e.g. rho_(atom1,atom1), rho_(atom1,atom2) = rho_(atom2,atom1), rho_(atom2,atom2)

           key=(Max(ai,aj)*(Max(ai,aj)-1))/2 + Min(ai,aj)
           k0=lstmet(key)

! first metal atom density and validity and truncation of potential

           If (Abs(dmes(1,k0,1)) > zero_plus .and. Nint(dmes(1,k0,1)) > 5) Then
              If (rrr <= dmes(3,k0,1)) Then

! interpolation parameters

                 rdr = 1.0_wp/dmes(4,k0,1)
                 rr1 = rrr - dmes(2,k0,1)
                 l   = Min(Nint(rr1*rdr),Nint(dmes(1,k0,1))-1)
                 If (l < 5) Then ! catch unsafe value
                    safe=.false.
                    Write(*,*) 'ccc',l,iatm,jatm,rrr
                    l=6
                 End If
                 ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

                 vk0 = dmes(l-1,k0,1)
                 vk1 = dmes(l  ,k0,1)
                 vk2 = dmes(l+1,k0,1)

                 t1 = vk1 + ppp*(vk1 - vk0)
                 t2 = vk1 + ppp*(vk2 - vk1)

                 If (ppp < 0.0_wp) Then ! density is a positive function!
                    density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                    If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
                 Else If (l == 5) Then
                    density = t2
                 Else
                    density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                    If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
                 End If

                 rhs(iatm) = rhs(iatm) + density
                 If (jatm <= natms) rhs(jatm) = rhs(jatm) + density

              End If
           End If

        Else If (tabmet == 4) Then ! 2BEEAM

! first metal atom density and validity and truncation of potential

           If (Abs(dmes(1,kj,1)) > zero_plus .and. Nint(dmes(1,kj,1)) > 5) Then
              If (rrr <= dmes(3,kj,1)) Then

! interpolation parameters

                 rdr = 1.0_wp/dmes(4,kj,1)
                 rr1 = rrr - dmes(2,kj,1)
                 l   = Min(Nint(rr1*rdr),Nint(dmes(1,kj,1))-1)
                 If (l < 5) Then ! catch unsafe value
                    safe=.false.
                    Write(*,*) 'ddd',l,iatm,jatm,rrr
                    l=6
                 End If
                 ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

                 vk0 = dmes(l-1,kj,1)
                 vk1 = dmes(l  ,kj,1)
                 vk2 = dmes(l+1,kj,1)

                 t1 = vk1 + ppp*(vk1 - vk0)
                 t2 = vk1 + ppp*(vk2 - vk1)

                 If (ppp < 0.0_wp) Then ! density is a positive function!
                    density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                    If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
                 Else If (l == 5) Then
                    density = t2
                 Else
                    density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                    If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
                 End If

                 rhs(iatm) = rhs(iatm) + density
                 If (ki == kj .and. jatm <= natms) rhs(jatm) = rhs(jatm) + density

              End If
           End If

! second metal atom density and validity and truncation of potential

           If (Abs(dmes(1,ki,1)) > zero_plus .and. Nint(dmes(1,ki,1)) > 5) Then
              If (ki /= kj .and. jatm <= natms) Then
                 If (rrr <= dmes(3,ki,1)) Then

! interpolation parameters

                    rdr = 1.0_wp/dmes(4,ki,1)
                    rr1 = rrr - dmes(2,ki,1)
                    l   = Min(Nint(rr1*rdr),Nint(dmes(1,ki,1))-1)
                    If (l < 5) Then ! catch unsafe value
                       safe=.false.
                       Write(*,*) 'eee',l,iatm,jatm,rrr
                       l=6
                    End If
                    ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

                    vk0 = dmes(l-1,ki,1)
                    vk1 = dmes(l  ,ki,1)
                    vk2 = dmes(l+1,ki,1)

                    t1 = vk1 + ppp*(vk1 - vk0)
                    t2 = vk1 + ppp*(vk2 - vk1)

                    If (ppp < 0.0_wp) Then ! density is a positive function!
                       density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                       If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
                    Else If (l == 5) Then
                       density = t2
                    Else
                       density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                       If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
                    End If

                    rhs(jatm) = rhs(jatm) + density

                 End If
              End If
           End If
        End If
     End If
  End Do

End Subroutine metal_ld_collect_eam

Subroutine metal_ld_collect_fst(iatm,rmet,rrt,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local atomic density for
! Finnis-Sinclair Type metal potentials
!
! Note: Designed to be used as part of metal_ld_compute
!
! copyright - daresbury laboratory
! author    - w.smith june 1995
! amended   - i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: rrt
  Logical,                                  Intent( InOut ) :: safe

  Integer           :: m,ai,aj,jatm,key,kmn,kmx,k0,k1,k2, &
                       keypot,nnn,mmm,l
  Real( Kind = wp ) :: rrr,rdr,rr1,vk0,vk1,vk2, &
                       t1,t2,density,eps,sig,   &
                       cc0,cc1,cc2,cc3,cc4,     &
                       aaa,bbb,ccc,ddd,ppp,qqq, &
                       bet,cut1,cut2,rr0

! global type of itam

  ai=ltype(iatm)

! start of primary loop for density

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

     If (ld_met) Then
        k1=Max(ai,aj)
        k2=Min(ai,aj)

        kmx=k1*(k1+1)/2
        kmn=k2*(k2+1)/2

        k1=lstmet(kmx)
        k2=lstmet(kmn)
     End If

! interatomic distance

     rrr=rrt(m)

! validity and truncation of analytic potential

     keypot=ltpmet(k0)
     If (keypot > 0 .and. rrr <= rmet) Then

! Abs(dmet(1,k0,1)) > zero_plus, as potentials are analytic

        If (ld_met) Then ! direct calculation

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

              density=0.0_wp
              If (rrr <= cut2) density=(rrr-ddd)**2+bet*(rrr-ddd)**3/ddd

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

              density=0.0_wp
              If (rrr <= cut2) density=(rrr-ddd)**2+bbb**2*(rrr-ddd)**4

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
              nnn=Nint(prmmet(3,k0))
              mmm=Nint(prmmet(4,k0))

              density=(sig/rrr)**mmm

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

              density=Exp(-2.0_wp*qqq*(rrr-rr0)/rr0)

              t1=prmmet(4,k0)**2
              t2=t1

           Else If (keypot == 5) Then

! many-body perturbation component only potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              mmm=Nint(prmmet(3,k0))

! interpolation parameters

              rdr = 1.0_wp/merf(4)
              rr1 = rrr - merf(2)
              l   = Min(Nint(rr1*rdr),Nint(merf(1))-1)
              If (l < 5) Then ! catch unsafe value
                 safe=.false.
                 Write(*,*) 'aaa',l,iatm,jatm,rrr
                 l=6
              End If
              ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

              vk0 = merf(l-1)
              vk1 = merf(l  )
              vk2 = merf(l+1)

              t1 = vk1 + ppp*(vk1 - vk0)
              t2 = vk1 + ppp*(vk2 - vk1)

              If (ppp < 0.0_wp) Then
                 density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              Else If (l == 5) Then
                 density = t2
              Else
                 density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              End If
              density=density*sig/rrr**mmm

              If (ai == aj) Then
                 t1=prmmet(1,k0)**2
                 t2=t1
              Else
                 t1=prmmet(1,k1)**2
                 t2=prmmet(1,k2)**2
              End If

           End If

           If (ai > aj) Then
              rho(iatm) = rho(iatm) + density*t1
              If (jatm <= natms) rho(jatm) = rho(jatm) + density*t2
           Else
              rho(iatm) = rho(iatm) + density*t2
              If (jatm <= natms) rho(jatm) = rho(jatm) + density*t1
           End If

        Else ! tabulated calculation

! truncation of potential

           If (rrr <= dmet(3,k0,1)) Then

! interpolation parameters

              rdr = 1.0_wp/dmet(4,k0,1)
              rr1 = rrr - dmet(2,k0,1)
              l   = Min(Nint(rr1*rdr),Nint(dmet(1,k0,1))-1)
              If (l < 5) Then ! catch unsafe value
                 safe=.false.
                 Write(*,*) 'bbb',l,iatm,jatm,rrr
                 l=6
              End If
              ppp = rr1*rdr - Real(l,wp)

! calculate density using 3-point interpolation

              vk0 = dmet(l-1,k0,1)
              vk1 = dmet(l  ,k0,1)
              vk2 = dmet(l+1,k0,1)

              t1 = vk1 + ppp*(vk1 - vk0)
              t2 = vk1 + ppp*(vk2 - vk1)

              If (ppp < 0.0_wp) Then
                 density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              Else If (l == 5) Then
                 density = t2
              Else
                 density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              End If

              If (ai > aj) Then
                 rho(iatm) = rho(iatm) + density*dmet(1,k0,2)
                 If (jatm <= natms) rho(jatm) = rho(jatm) + density*dmet(2,k0,2)
              Else
                 rho(iatm) = rho(iatm) + density*dmet(2,k0,2)
                 If (jatm <= natms) rho(jatm) = rho(jatm) + density*dmet(1,k0,2)
              End If

           End If

        End If

     End If

  End Do

End Subroutine metal_ld_collect_fst

Subroutine metal_table_derivatives(ityp,buffer,v2d,vvv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to calculate numerical derivatives of tabulated
! EAM metal potentials
!
! copyright - daresbury laboratory
! author    - w.smith march 2006
! amended   - i.t.todorov april 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: ityp,v2d
  Real( Kind = wp ), Intent( In    ) :: buffer(1:mxgmet)
  Real( Kind = wp ), Intent( InOut ) :: vvv(1:mxgmet,1:v2d,1:2)

  Integer           :: i,v_end,i_start,i_end
  Real( Kind = wp ) :: delmet,aa0,aa1,aa2,aa3,aa4, &
                       d1y,d2y,d3y,d4y,f0,f1,f2,f3,f4

! interpolation parameters

  vvv(1,ityp,2)=buffer(1)
  vvv(2,ityp,2)=buffer(2)
  vvv(3,ityp,2)=buffer(3)
  vvv(4,ityp,2)=buffer(4)

! construct interpolation table

  delmet=buffer(4)
  v_end=Nint(buffer(1))
  i_start=5    +2
  i_end  =v_end-2
  Do i=i_start,i_end
     aa0=buffer(i)
     If (Abs(aa0) <= zero_plus) Then
        f0=0.0_wp
        f1=0.0_wp
        f2=0.0_wp
        f3=0.0_wp
        f4=0.0_wp
     Else
        f0=buffer(i-2)/aa0
        f1=buffer(i-1)/aa0
        f2=1.0_wp
        f3=buffer(i+1)/aa0
        f4=buffer(i+2)/aa0
     End If

! calculate numerical differences for 5-point interpolation

     d1y=(f1-f0)
     d2y=(f2-f1)-(f1-f0)
     d3y=(f3-f0)+3.0_wp*(f1-f2)
     d4y=(f4-f3)+3.0_wp*(f2-f3)+3.0_wp*(f2-f1)+(f0-f1)

! calculate polynomial coefficients

     aa0=aa0/delmet
     aa4=d4y/24.0_wp
     aa3=(d3y+12.0_wp*aa4)/6.0_wp
     aa2=(d2y+6.0_wp*aa3-14.0_wp*aa4)/2.0_wp
     aa1=d1y+3.0_wp*aa2-7.0_wp*aa3+15.0_wp*aa4

! calculate derivatives

     vvv(i,ityp,2)=aa1*aa0

! derivatives at extremes of range

     If      (i == i_start) Then
        vvv(i_start-2,ityp,2)=(aa1-4.0_wp*aa2+12.0_wp*aa3-32.0_wp*aa4)*aa0
        vvv(i_start-1,ityp,2)=(aa1-2.0_wp*aa2+3.0_wp*aa3-4.0_wp*aa4)*aa0
     Else If (i == i_end  ) Then
        vvv(i_end  +1,ityp,2)=(aa1+2.0_wp*aa2+3.0_wp*aa3+4.0_wp*aa4)*aa0
        vvv(i_end  +2,ityp,2)=(aa1+4.0_wp*aa2+12.0_wp*aa3+32.0_wp*aa4)*aa0
     End If
  End Do

! set derivatives to zero beyond end point of function

  Do i=v_end+3,mxgmet
     vvv(i,ityp,2)=0.0_wp
  End Do

End Subroutine metal_table_derivatives

Subroutine metal_ld_export(mdir,mlast,ixyz0,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export metal density data in domain boundary
! regions for halo formation
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer, Intent( In    ) :: mdir
  Integer, Intent( InOut ) :: mlast,ixyz0(1:mxatms)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical :: safe,lrhs
  Integer :: fail,iadd,limit,iblock,          &
             i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
             jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message
! Number of transported quantities per particle

  If (.not.(tabmet == 3 .or. tabmet == 4)) Then
     lrhs=.false.
     iadd=2
  Else
     lrhs=.true.
     iadd=3
  End If

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'metal_ld_export allocation failure'
     Call error(0,message)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/Merge(2,1,comm%mxnode > 1)

! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
! kxyz - corrected halo reduction factor particles haloing both +&- sides
! jdnode - destination (send to), kdnode - source (receive from)

  kx = 0 ; ky = 0 ; kz = 0
  If      (mdir == -1) Then ! Direction -x
     kx  = 1
     jxyz= 1
     kxyz= 3

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= 3

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= 30

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= 30

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= 300

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= 300

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(47)
  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PARTICLES ON THIS NODE


  Do i=1,mlast

! If the particle is within the remaining 'inverted halo' of this domain

     If (ixyz0(i) > 0) Then

! Get the necessary halo indices

        ix=Mod(ixyz0(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz0(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz0(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

! If safe to proceed

           If ((imove+iadd) <= iblock) Then

! pack particle density and halo indexing

              If (.not.lrhs) Then
                 buffer(imove+1)=rho(i)

! Use the corrected halo reduction factor when the particle is halo to both +&- sides

                 buffer(imove+2)=Real(ixyz0(i)-Merge(jxyz,kxyz,j == jxyz),wp)
              Else
                 buffer(imove+1)=rho(i)
                 buffer(imove+2)=rhs(i)

! Use the corrected halo reduction factor when the particle is halo to both +&- sides

                 buffer(imove+3)=Real(ixyz0(i)-Merge(jxyz,kxyz,j == jxyz),wp)
              End If

           Else

              safe=.false.


           End If

           imove=imove+iadd

        End If

     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  Call gcheck(comm,safe)
  If (.not.safe) Then
     itmp=Merge(2,1,comm%mxnode > 1)*imove
     Call gmax(comm,itmp)
     Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
     Call error(38)
  End If

! exchange information on buffer sizes

  If (comm%mxnode > 1) Then
     Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,MetLdExp_tag,comm%comm,comm%request,comm%ierr)
     Call gsend(comm,imove,jdnode,MetLdExp_tag)
     Call gwait(comm)
  Else
     jmove=imove
  End If

! Check for array bound overflow (can arrays cope with incoming data)

  safe=((mlast+jmove/iadd) <= mxatms)
  Call gcheck(comm,safe)
  If (.not.safe) Then
     itmp=mlast+jmove/iadd
     Call gmax(comm,itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
     Call error(39)
  End If

! exchange buffers between nodes (this is a MUST)

  If (comm%mxnode > 1) Then
     If (jmove > 0) Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,MetLdExp_tag,comm%comm,comm%request,comm%ierr)
     If (imove > 0) Then
       Call gsend(comm,buffer(1:imove),jdnode,MetLdExp_tag)
     End If
     If (jmove > 0) Call gwait(comm)
  End If

! load transferred data

  j=Merge(iblock,0,comm%mxnode > 1)
  Do i=1,jmove/iadd
     mlast=mlast+1

! unpack particle density and remaining halo indexing

     If (.not.lrhs) Then
        rho(mlast) =buffer(j+1)
        ixyz0(mlast)=Nint(buffer(j+2))
     Else
        rho(mlast) =buffer(j+1)
        rhs(mlast) =buffer(j+2)
        ixyz0(mlast)=Nint(buffer(j+3))
     End If

     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'metal_ld_export deallocation failure'
     Call error(0,message)
  End If

End Subroutine metal_ld_export
Subroutine metal_ld_set_halo(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of density data between
! neighbouring domains/nodes
!
! Note: all depends on the ixyz halo array set in set_halo, this assumes
!       that (i) rmet=rcut! as well as (ii) all the error checks in there
!
! copyright - daresbury laboratory
! amended   - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( comms_type ), Intent( InOut ) :: comm

  Logical :: safe
  Integer :: fail,mlast

  Integer, Allocatable :: ixyz0(:)
  Character( Len = 256 ) :: message

  fail = 0
  Allocate (ixyz0(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'metal_ld_set_halo allocation failure'
     Call error(0,message)
  End If
  ixyz0(1:nlast) = ixyz(1:nlast)

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call metal_ld_export(-1,mlast,ixyz0,comm)
  Call metal_ld_export( 1,mlast,ixyz0,comm)

! exchange atom data in -/+ y directions

  Call metal_ld_export(-2,mlast,ixyz0,comm)
  Call metal_ld_export( 2,mlast,ixyz0,comm)

! exchange atom data in -/+ z directions

  Call metal_ld_export(-3,mlast,ixyz0,comm)
  Call metal_ld_export( 3,mlast,ixyz0,comm)

! check atom totals after data transfer

  safe=(mlast == nlast)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(96)

  Deallocate (ixyz0, Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'metal_ld_set_halo deallocation failure'
     Call error(0,message)
  End If

End Subroutine metal_ld_set_halo


End Module metal
