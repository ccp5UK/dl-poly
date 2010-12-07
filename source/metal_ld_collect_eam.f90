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
! amended   - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltype,list
  Use metal_module,  Only : dmet

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: rsqdf
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

     If (Abs(dmet(1,aj,1)) > zero_plus .and. rsq <= dmet(3,aj,1)**2) Then

! interpolation parameters

        rdr = 1.0_wp/dmet(4,aj,1)
        rrr = Sqrt(rsq) - dmet(2,aj,1)
        l   = Nint(rrr*rdr)
        ppp = rrr*rdr - Real(l,wp)

! catch unsafe value

        If (l < 2) Then
           safe=.false.
           l=2
        End If

! calculate density using 3-point interpolation

        vk0 = dmet(l+3,aj,1)
        vk1 = dmet(l+4,aj,1)
        vk2 = dmet(l+5,aj,1)

        t1 = vk1 + ppp*(vk1 - vk0)
        t2 = vk1 + ppp*(vk2 - vk1)

        If (ppp < 0.0_wp) Then
           density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
        Else
           density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
        End If

        rho(iatm) = rho(iatm) + density
        If (ai == aj .and. jatm <= natms) rho(jatm) = rho(jatm) + density

     End If

! second metal atom density and validity and truncation of potential

     If (ai /= aj .and. jatm <= natms .and. &
         Abs(dmet(1,ai,1)) > zero_plus .and. rsq <= dmet(3,ai,1)**2) Then

! interpolation parameters

        rdr = 1.0_wp/dmet(4,ai,1)
        rrr = Sqrt(rsq) - dmet(2,ai,1)
        l   = Nint(rrr*rdr)
        ppp = rrr*rdr - Real(l,wp)

! catch unsafe value

        If (l < 2) Then
           safe=.false.
           l=2
        End If

! calculate density using 3-point interpolation

        vk0 = dmet(l+3,ai,1)
        vk1 = dmet(l+4,ai,1)
        vk2 = dmet(l+5,ai,1)

        t1 = vk1 + ppp*(vk1 - vk0)
        t2 = vk1 + ppp*(vk2 - vk1)

        If (ppp < 0.0_wp) Then
           density = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
        Else
           density = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
        End If

        rho(jatm) = rho(jatm) + density

     End If

  End Do

End Subroutine metal_ld_collect_eam
