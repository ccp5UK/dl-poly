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
! amended   - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltype,list
  Use metal_module,  Only : ld_met,lstmet,ltpmet,prmmet,dmet,merf,rho

  Implicit None

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

! start of primary loop for density

  ai=ltype(iatm)

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
