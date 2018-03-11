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

  Use kinds, only : wp
  Use setup_module
  Use configuration, Only : natms,ltype,list
  Use metal_module,  Only : l2bmet,tabmet,lstmet,dmet,dmes,rho,rhs
  Use site_module,   Only : ntpatm

  Implicit None

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
