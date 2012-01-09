Subroutine metal_ld_compute                &
           (imcon,rmet,keyfce,elrcm,vlrcm, &
           xdf,ydf,zdf,rsqdf,              &
           rho,engden,virden,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local density in metals using
! the verlet neighbour list and sutton-chen potentials
!
! Note: Designed to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov may 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : mxnode,gsum,gcheck
  Use setup_module
  Use config_module, Only : cell,natms,ltype,list,xxx,yyy,zzz
  Use metal_module,  Only : ntpmet,ltpmet,fmet

  Implicit None

  Integer,                                  Intent( In    ) :: imcon,keyfce
  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( In    ) :: elrcm,vlrcm
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( InOut ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent(   Out ) :: rho
  Real( Kind = wp ),                        Intent(   Out ) :: engden,virden
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical, Save     :: newjob = .true.
  Integer, Save     :: keypot

  Logical           :: safe = .true.
  Integer           :: i,limit,j,k,l,k0
  Real( Kind = wp ) :: rhosqr,rdr,rrr,ppp,fk0,fk1,fk2,t1,t2

! check on mixing metal types

  If (newjob) Then
     newjob = .false.

     keypot=0
     Do i=1,ntpmet
        keypot=ltpmet(i)
        If (i > 1) Then
           If (keypot /= ltpmet(i-1)) Call error(92)
        End If
     End Do
  End If

! initialise energy and virial accumulators

  engden=0.0_wp
  virden=0.0_wp

! initialise density array

  rho=0.0_wp

! calculate local atomic density
! outer loop over atoms

  Do i=1,natms

! Get list limit

     limit=list(0,i)

! calculate interatomic distances

     Do k=1,limit
        j=list(k,i)

        xdf(k)=xxx(i)-xxx(j)
        ydf(k)=yyy(i)-yyy(j)
        zdf(k)=zzz(i)-zzz(j)
     End Do

! periodic boundary conditions

     Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distances

     Do k=1,limit
        rsqdf(k) = xdf(k)**2+ydf(k)**2+zdf(k)**2
     End Do

! calculate contributions to local density

     If (keypot == 0) Then ! EAM contributions
        Call metal_ld_collect_eam(i,rsqdf,rho,safe)
     Else                  ! FST contributions
        Call metal_ld_collect_fst(i,rsqdf,rho,rmet)
     End If
  End Do

! Check safety for densities

  If (keypot == 0) Then
     If (mxnode > 1) Call gcheck(safe)
     If (.not.safe) Call error(506)
  End If

  Do i=1,natms

! calculate density terms to energy and virial

     If (keypot == 0) Then ! EAM potential

! potential function index

        k0=ltype(i)

! validity of potential

        If (Abs(fmet(1,k0,1)) > zero_plus) Then

! check for unsafe densities (mind start was shifted)

           If (rho(i) >= fmet(2,k0,1)+fmet(4,k0,1) .and. rho(i) <= fmet(3,k0,1)) Then

! interpolation parameters

              rdr = 1.0_wp/fmet(4,k0,1)
              rrr = rho(i) - fmet(2,k0,1)
              l   = Nint(rrr*rdr)
              ppp = rrr*rdr - Real(l,wp)

! catch unsafe value

              If (l < 1) Then
                 safe=.false.
                 l=2
              End If

! calculate embedding energy using 3-point interpolation

              fk0 = fmet(l+3,k0,1)
              fk1 = fmet(l+4,k0,1)
              fk2 = fmet(l+5,k0,1)

              t1 = fk1 + ppp*(fk1 - fk0)
              t2 = fk1 + ppp*(fk2 - fk1)

              If (ppp < 0.0_wp) Then
                 engden = engden + t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              Else
                 engden = engden + t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              End If

! calculate derivative of embedding function wrt density
! using 3-point interpolation and store result in rho array

              fk0 = fmet(l+3,k0,2)
              fk1 = fmet(l+4,k0,2)
              fk2 = fmet(l+5,k0,2)

              t1 = fk1 + ppp*(fk1 - fk0)
              t2 = fk1 + ppp*(fk2 - fk1)

              If (ppp < 0.0_wp) Then
                 rho(i) = t1 + 0.5_wp*(t2-t1)*(ppp+1.0_wp)
              Else
                 rho(i) = t2 + 0.5_wp*(t2-t1)*(ppp-1.0_wp)
              End If

           Else

              safe=.false.

           End If

        End If

     Else ! FST of metal potentials

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

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(507)

! virial term (averaged per node)

  If (mxnode > 1) Call gsum(virden)
  virden=virden/Real(mxnode,wp)

! calculate stress tensor (density contributions are to
! diagonal elements only)

  stress(1)=stress(1)-virden/3.0_wp
  stress(5)=stress(5)-virden/3.0_wp
  stress(9)=stress(9)-virden/3.0_wp

! obtain atomic densities for outer border regions

  Call metal_ld_set_halo(imcon,rmet,keyfce,rho)

End Subroutine metal_ld_compute
