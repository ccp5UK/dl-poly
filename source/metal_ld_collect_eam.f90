Subroutine metal_ld_collect_eam(iatm,rsqdf,rho,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local atomic density for
! Embedded Atom Model & Extended Embedded Atom Model metal potentials
!
! Note: Designed to be used as part of metal_ld_compute
!
! copyright - daresbury laboratory
! author    - w.smith june 1995
! amended   - i.t.todorov june 2012
! contrib   - r.davidchak june 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltype,list
  Use metal_module,  Only : tabmet,dmet
  Use site_module,   Only : ntpatm

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: rsqdf
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: rho
  Logical,                                  Intent( InOut ) :: safe

  Integer           :: m,ai,ki,jatm,aj,kj,l
  Real( Kind = wp ) :: rsq,rdr,rrr,ppp,vk0,vk1,vk2,t1,t2,density

! start of primary loop for density

  ai=ltype(iatm)

  Do m=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(m,iatm)
     aj=ltype(jatm)

     If      (tabmet == 1) Then ! EAM
        ki=ai
        kj=aj
     Else If (tabmet == 2) Then ! EEAM
        ki=(aj-1)*ntpatm+ai ! aj-ai
        kj=(ai-1)*ntpatm+aj ! ai-aj
     End If

! interatomic distance

     rsq=rsqdf(m)

! first metal atom density and validity and truncation of potential

     If (Abs(dmet(1,kj,1)) > zero_plus) Then
        If (rsq <= dmet(3,kj,1)**2) Then

! interpolation parameters

           rdr = 1.0_wp/dmet(4,kj,1)
           rrr = Sqrt(rsq) - dmet(2,kj,1)
           l   = Min(Nint(rrr*rdr),Nint(dmet(1,kj,1))-1)
           If (l < 6) Then ! catch unsafe value
              safe=.false.
              l=6
           End If
           ppp = rrr*rdr - Real(l,wp)

! calculate density using 3-point interpolation

           vk0 = dmet(l-1,kj,1)
           vk1 = dmet(l  ,kj,1)
           vk2 = dmet(l+1,kj,1)

           t1 = vk1 + ppp*(vk1 - vk0)
           t2 = vk1 + ppp*(vk2 - vk1)

           If (ppp < 0.0_wp) Then ! density is a positive function!
              density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
           Else
              density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
           End If

           rho(iatm) = rho(iatm) + density
           If (ki == kj .and. jatm <= natms) rho(jatm) = rho(jatm) + density

        End If
     End If

! second metal atom density and validity and truncation of potential

     If (Abs(dmet(1,ki,1)) > zero_plus) Then
        If (ki /= kj .and. jatm <= natms) Then
           If (rsq <= dmet(3,ki,1)**2) Then

! interpolation parameters

              rdr = 1.0_wp/dmet(4,ki,1)
              rrr = Sqrt(rsq) - dmet(2,ki,1)
              l   = Min(Nint(rrr*rdr),Nint(dmet(1,ki,1))-1)
              If (l < 6) Then ! catch unsafe value
                 safe=.false.
                 l=6
              End If
              ppp = rrr*rdr - Real(l,wp)

! calculate density using 3-point interpolation

              vk0 = dmet(l-1,ki,1)
              vk1 = dmet(l  ,ki,1)
              vk2 = dmet(l+1,ki,1)

              t1 = vk1 + ppp*(vk1 - vk0)
              t2 = vk1 + ppp*(vk2 - vk1)

              If (ppp < 0.0_wp) Then ! density is a positive function!
                 density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
                 If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
              Else
                 density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
                 If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
              End If

              rho(jatm) = rho(jatm) + density

           End If
        End If
     End If

  End Do

End Subroutine metal_ld_collect_eam
