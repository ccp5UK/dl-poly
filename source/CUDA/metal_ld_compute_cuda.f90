! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Center for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

Subroutine metal_ld_compute_get_keypot(keypotr)
  Use metal_module,  Only : ntpmet,ltpmet
  Implicit None
  Integer, Intent(Out) :: keypotr
  Integer           :: i
  Logical, Save     :: newjob = .true.
  Integer, Save     :: keypot

  If (newjob) Then
     newjob=.false.

     keypot=0
     Do i=1,ntpmet
        keypot=ltpmet(i)
        If (i > 1) Then
           If (keypot /= ltpmet(i-1)) Call error(92)
        End If
     End Do
  End If
  keypotr = keypot
End Subroutine metal_ld_compute_get_keypot

Subroutine metal_ld_compute         &
           (imcon,rmet,elrcm,vlrcm, &
           xdf,ydf,zdf,rsqdf,       &
           engden,virden,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating local density in metals using
! the verlet neighbour list and sutton-chen potentials
!
! Note: Designed to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov june 2013
! contrib   - r.davidchak (eeam) june 2012
! contrib   - b.palmer (2band) may 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : mxnode,gsum,gcheck
  Use setup_module
  Use config_module, Only : cell,natms,ltg,ltype,list,xxx,yyy,zzz
  Use metal_module,  Only : ls_met,l2bmet,ntpmet,ltpmet, &
                            fmet,fmes,rho,rhs, ld_met, lstmet, vmet, dmet
                            !fmet,fmes,rho,rhs !Modified by BBG: Changed as
                            !above

#ifdef COMPILE_CUDA
  Use dl_poly_cuda_module
#endif

  Implicit None

  Integer,                                  Intent( In    ) :: imcon
  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( In    ) :: elrcm,vlrcm
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( InOut ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ),                        Intent(   Out ) :: engden,virden
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical, Save     :: newjob = .true.
  Integer, Save     :: keypot

  Logical           :: safe = .true.
  Integer           :: limit,i,j,k,l,k0
  Real( Kind = wp ) :: rhosqr,rdr,rrr,ppp,fk0,fk1,fk2,t1,t2

#ifdef COMPILE_CUDA
  Call start_timing_metal_ld_compute()
#endif

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
  If (l2bmet) rhs = 0.0_wp

! All calls below act on rho (rhs)

! calculate local atomic density
! outer loop over atoms

#ifdef COMPILE_CUDA
! The CUDA port implements the features of dl_poly_3 - a subset of those found in dl_poly_4
! In particular, only tabulated calculations are available in metal_ld_compute
! The CUDA acceleration is not called if direct calculation is required
  !If (dl_poly_cuda_offload_metal_ld_compute() .and. ld_met .eqv. .false.) Then !*CHANGE == to .eqv. !Modified by BBG: Line is changed as follows.
  If (dl_poly_cuda_offload_metal_ld_compute() .eqv. .true. .and. ld_met .eqv. .false.) Then !*CHANGE == to .eqv.
     Call metal_ld_compute_cuda_initialise(& 
          0,mxatms,natms,mxgrid,ntpmet,mxmet,mxatdm,mxlist,&
          xxx,yyy,zzz,list,ltype,ltpmet,lstmet,vmet,dmet,cell,rho)
     Call metal_ld_compute_cuda_invoke()
     Call metal_ld_compute_cuda_finalise()
     !Modified by BBG: Above 5 lines were uncommented.
  Else
#endif

  Do i=1,natms
     limit=list(0,i) ! Get list limit

! calculate interatomic distances
!PART 1 - P1
     Do k=1,limit
        j=list(k,i)

        xdf(k)=xxx(i)-xxx(j)
        ydf(k)=yyy(i)-yyy(j)
        zdf(k)=zzz(i)-zzz(j)
     End Do

! periodic boundary conditions
!PART 2 - P2
     Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distances
!PART 3 - P3
     Do k=1,limit
        rsqdf(k) = xdf(k)**2+ydf(k)**2+zdf(k)**2
     End Do

! calculate contributions to local density
!PART 4 - P4
     If (keypot == 0) Then ! EAM contributions
        Call metal_ld_collect_eam(i,rsqdf,safe)
     Else                  ! FST contributions
        Call metal_ld_collect_fst(i,rsqdf,safe,rmet)
     End If
  End Do

#ifdef COMPILE_CUDA
  End If
#endif

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

                   rhs(i) = rhs(i) + 0.0_wp

                 End If
              Else
                 Write(*,*) 'bad density range problem: (LTG,RHS) ',ltg(i),rhs(i)
                 safe=.false.
              End If

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

  Call metal_ld_set_halo()

#ifdef COMPILE_CUDA
  Call stop_timing_metal_ld_compute()
#endif

End Subroutine metal_ld_compute
