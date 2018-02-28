Subroutine metal_ld_compute(rmet,elrcm,vlrcm,engden,virden,stress)

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

  Use kinds, only : wp
  Use comms_module,  Only : idnode,mxnode,gsum,gcheck
  Use setup_module
  Use config_module, Only : natms,ltg,ltype,list,xxx,yyy,zzz
  Use metal_module,  Only : ls_met,l2bmet,tabmet,fmet,fmes,rho,rhs

  Implicit None

  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( In    ) :: elrcm,vlrcm
  Real( Kind = wp ),                        Intent(   Out ) :: engden,virden
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical           :: safe = .true.
  Integer           :: fail,limit,i,j,k,l,k0
  Real( Kind = wp ) :: rhosqr,rdr,rrr,ppp,fk0,fk1,fk2,t1,t2

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

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
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_compute allocation failure, node: ', idnode
     Call error(0)
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
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_compute allocation failure, node: ', idnode
     Call error(0)
  End If

! Check safety for densities of EAM and MBPC

  If (mxnode > 1) Call gcheck(safe)
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
                 Write(*,*) 'bad density range problem: (LTG,RHS) ',ltg(i),rhs(i)
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

  Call metal_ld_set_halo()

End Subroutine metal_ld_compute
