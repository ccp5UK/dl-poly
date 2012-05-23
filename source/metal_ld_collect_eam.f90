Subroutine metal_ld_collect_eam(iatm,rsqdf,rho,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local atomic density for
! Embedded Atom Model metal potentials
!
! Note: Designed to be used as part of metal_ld_compute
!
! copyright - daresbury laboratory
! author    - w.smith june 1995
! amended   - i.t.todorov may 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltg,ltype,list
  Use metal_module,  Only : dmet

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: rsqdf
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: rho
  Logical,                                  Intent( InOut ) :: safe

  Integer           :: m,ai,aj,jatm,l
  Real( Kind = wp ) :: rsq,rdr,rrr,ppp,vk0,vk1,vk2,t1,t2,density

! start of primary loop for density

  ai=ltype(iatm)

  Do m=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(m,iatm)
     aj=ltype(jatm)

! interatomic distance

     rsq=rsqdf(m)

! first metal atom density and validity and truncation of potential

     If (Abs(dmet(1,aj,1)) > zero_plus) Then
        If (rsq <= dmet(3,aj,1)**2) Then

! interpolation parameters

           rdr = 1.0_wp/dmet(4,aj,1)
           rrr = Sqrt(rsq) - dmet(2,aj,1)
           l   = Min(Nint(rrr*rdr),Int(dmet(1,aj,1))-1)
           If (l < 2) Then ! catch unsafe value
              safe=.false.
              l=2
           End If
           ppp = rrr*rdr - Real(l,wp)

! calculate density using 3-point interpolation

           vk0 = dmet(l-1,aj,1)
           vk1 = dmet(l  ,aj,1)
           vk2 = dmet(l+1,aj,1)

           t1 = vk1 + ppp*(vk1 - vk0)
           t2 = vk1 + ppp*(vk2 - vk1)

           If (ppp < 0.0_wp) Then
              density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
           Else
              density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
           End If
           If (density < 0.0_wp) density = -density ! for non-smooth descend to zero, or ascend from zero
           rho(iatm) = rho(iatm) + density

           If (rho(iatm) < -zero_plus .or. density < -zero_plus) &
              Write(*,*) 'negative density: (LTG,RHO_SUM,RHO) ',ltg(iatm),rho(iatm),density,l,vk0,vk1,vk2,t1,t2,ppp

           If (ai == aj .and. jatm <= natms) Then
              rho(jatm) = rho(jatm) + density
              If (rho(jatm) < -zero_plus .or. density < -zero_plus) &
                 Write(*,*) 'negative density: (LTG,RHO_SUM,RHO) ',ltg(jatm),rho(jatm),density,l,vk0,vk1,vk2,t1,t2,ppp
           End If

        End If
     End If

! second metal atom density and validity and truncation of potential

     If (Abs(dmet(1,ai,1)) > zero_plus) Then
        If (ai /= aj .and. jatm <= natms .and. rsq <= dmet(3,ai,1)**2) Then

! interpolation parameters

           rdr = 1.0_wp/dmet(4,ai,1)
           rrr = Sqrt(rsq) - dmet(2,ai,1)
           l   = Min(Nint(rrr*rdr),Int(dmet(1,ai,1))-1)
           If (l < 2) Then ! catch unsafe value
              safe=.false.
              l=2
           End If
           ppp = rrr*rdr - Real(l,wp)

! calculate density using 3-point interpolation

           vk0 = dmet(l-1,ai,1)
           vk1 = dmet(l  ,ai,1)
           vk2 = dmet(l+1,ai,1)

           t1 = vk1 + ppp*(vk1 - vk0)
           t2 = vk1 + ppp*(vk2 - vk1)

           If (ppp < 0.0_wp) Then
              density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
           Else
              density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
           End If
           If (density < 0.0_wp) density = -density ! for non-smooth descend to zero, or ascend from zero
           rho(jatm) = rho(jatm) + density

           If (rho(jatm) < -zero_plus .or. density < -zero_plus) &
              Write(*,*) 'negative density: (LTG,RHO_SUM,RHO) ',ltg(jatm),rho(jatm),density,l,vk0,vk1,vk2,t1,t2,ppp

        End If
     End If

  End Do

End Subroutine metal_ld_collect_eam
