Subroutine two_body_forces                        &
           (rcut,rlnk,rvdw,rmet,pdplnc,keyens,    &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           lrdf,nstrdf,leql,nsteql,nstep,         &
           elrc,virlrc,elrcm,vlrcm,               &
           engcpe,vircpe,engsrp,virsrp,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating interatomic forces and rdf
! using the verlet neighbour list
!
! ntpvdw > 0 ------ switch for vdw potentials calculation
! ntpmet > 0 ------ switch for metal local density and potentials
!                   calculations
!
! ELECTROSTATICS KEYS
!
! keyfce = 0 ------ no electrostatics
! keyfce = 2 ------ Ewald sum (ewald_spme,ewald_real,ewald_excl)
! keyfce = 4 ------ distance dependent dielectric potential (coul_dddp)
! keyfce = 6 ------ coulombic 1/r potential (coul_cp)
! keyfce = 8 ------ force-shifted and damped coulombic potential (coul_fscp)
! keyfce =10 ------ reaction field and damped coulombic potential (coul_rfp)
! keyfce =12 ------ direct space Poisson solver (possion_module)
!
! nstfce - the rate at which the k-space contributions of SPME are
!          refreshed.  Once every 1 <= nstfce <= 7 steps.
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
! contrib   - h.a.boateng february 2016
! contrib   - p.s.petkov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use site_module,    Only : ntpatm,unqatm
  Use config_module,  Only : volm,sumchg,natms,list,xxx,yyy,zzz
  Use vnl_module,     Only : l_vnl
  Use ewald_module
  Use mpoles_module,  Only : induce
  Use poisson_module, Only : poisson_forces,poisson_excl_forces,poisson_frzn_forces
  Use vdw_module,     Only : ntpvdw
  Use metal_module,   Only : ntpmet
  Use kim_module
  Use rdf_module,     Only : ncfrdf

  Implicit None

  Logical,                                  Intent( In    ) :: lbook,lrdf,leql
  Integer,                                  Intent( In    ) :: keyens,        &
                                                               keyfce,nstfce, &
                                                               megfrz,nstrdf, &
                                                               nsteql,nstep
  Real( Kind = wp ),                        Intent( In    ) :: rcut,rlnk,rvdw,rmet, &
                                                               pdplnc,alpha,epsq
  Real( Kind = wp ),                        Intent( In    ) :: elrc,virlrc
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( InOut ) :: elrcm,vlrcm
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe, &
                                                               engsrp,virsrp
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress


  Logical,           Save :: new_nz    = .true.
  Real( Kind = wp ), Save :: factor_nz = 0.0_wp

  Logical           :: safe = .true., l_do_rdf
  Integer           :: fail,i,j,k,limit
  Real( Kind = wp ) :: engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl, &
                       engcpe_ex,vircpe_ex,engcpe_fr,vircpe_fr, &
                       engcpe_nz,vircpe_nz,          vircpe_dt, &
                       engden,virden,engmet,virmet,             &
                       engvdw,virvdw,engkim,virkim,             &
                       engacc,viracc,tmp,buffer(0:17)

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

  fail=0
  Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'two_body_forces allocation failure, node: ', idnode
     Call error(0)
  End If

  l_do_rdf = (lrdf .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstrdf) == 0)

! If k-space SPME is evaluated infrequently check whether
! at this timestep to evaluate or "refresh" with old values.
! At restart allocate the "refresh" arrays and force a fresh
! evaluation.  Repeat the same but only for the SPME k-space
! frozen-frozen evaluations in constant volume ensembles only.

  If (keyfce == 2 .or. keyfce == 12) Call ewald_check(keyens,megfrz,nsteql,nstfce,nstep)

! initialise energy and virial accumulators

  engkim = 0.0_wp
  virkim = 0.0_wp


  engden    = 0.0_wp
  virden    = 0.0_wp

  engmet    = 0.0_wp
  virmet    = 0.0_wp

  engvdw    = 0.0_wp
  virvdw    = 0.0_wp

  engsrp    = 0.0_wp
  virsrp    = 0.0_wp


  engcpe_rc = 0.0_wp
  vircpe_rc = 0.0_wp

  engcpe_rl = 0.0_wp
  vircpe_rl = 0.0_wp

  engcpe_ex = 0.0_wp
  vircpe_ex = 0.0_wp

  engcpe_fr = 0.0_wp
  vircpe_fr = 0.0_wp

  engcpe_nz = 0.0_wp
  vircpe_nz = 0.0_wp

  vircpe_dt = 0.0_wp

  engcpe    = 0.0_wp
  vircpe    = 0.0_wp

! Set up non-bonded interaction (verlet) list using link cells

  If ((.not.induce) .and. l_vnl) Call link_cell_pairs(rlnk,pdplnc,lbook,megfrz)

! Calculate all contributions from KIM

  If (kim /= ' ') Then
     Call kim_setup(ntpatm,unqatm,kim)
     Call kim_forces(engkim,virkim,stress)
     Call kim_cleanup()
  End If

  If (ntpmet > 0) Then

! Reset metal long-range corrections (constant pressure/stress only)

     If (keyens >= 20) Call metal_lrc(rmet,elrcm,vlrcm)

! calculate local density in metals

     Call metal_ld_compute(rmet,elrcm,vlrcm,engden,virden,stress)

  End If

! calculate coulombic forces, Ewald sum - fourier contribution

  If (keyfce == 2 .and. l_fce) Then
     If (mximpl > 0) Then
        If (mxompl <= 2) Then
           Call ewald_spme_mforces_d(alpha,epsq,engcpe_rc,vircpe_rc,stress)
        Else
           Call ewald_spme_mforces(alpha,epsq,engcpe_rc,vircpe_rc,stress)
        End If
     Else
        Call ewald_spme_forces(alpha,epsq,engcpe_rc,vircpe_rc,stress)
     End If
  End If

! outer loop over atoms

  Do i=1,natms

! Get list limit

     limit=list(0,i)

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

! distances, thanks to Alin Elena (one too many changes)

     Do k=1,limit
        rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
     End Do

! calculate metal forces and potential

     If (ntpmet > 0) Then
        Call metal_forces(i,rmet,xxt,yyt,zzt,rrt,engacc,viracc,stress,safe)

        engmet=engmet+engacc
        virmet=virmet+viracc
     End If

! calculate short-range force and potential terms

     If (ntpvdw > 0) Then
        Call vdw_forces(i,rvdw,xxt,yyt,zzt,rrt,engacc,viracc,stress)

        engvdw=engvdw+engacc
        virvdw=virvdw+viracc
     End If

!!!!!!!!!!!!!!!!!!!!!!!!!
! COULOMBIC CONTRIBUTIONS
!!!!!!!!!!!!!!!!!!!!!!!!1

     If (mximpl > 0) Then

!!! MULTIPOLAR ATOMIC SITES

        If (keyfce == 2) Then

! calculate coulombic forces, Ewald sum - real space contribution

           If (mxompl <= 2) Then
              Call ewald_real_mforces_d(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)
           Else
              Call ewald_real_mforces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)
           End If

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 4) Then

! distance dependant dielectric potential

           Call coul_dddp_mforces(i,rcut,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 6) Then

! coulombic 1/r potential with no truncation or damping

           Call coul_cp_mforces(i,rcut,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 8) Then

! force-shifted coulomb potentials

           Call coul_fscp_mforces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 10) Then

! reaction field potential

           Call coul_rfp_mforces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        End If

     Else

        If (keyfce == 2) Then

! calculate coulombic forces, Ewald sum - real space contribution

           Call ewald_real_forces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 4) Then

! distance dependant dielectric potential

           Call coul_dddp_forces(i,rcut,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 6) Then

! coulombic 1/r potential with no truncation or damping

           Call coul_cp_forces(i,rcut,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 8) Then

! force-shifted coulomb potentials

           Call coul_fscp_forces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 10) Then

! reaction field potential

           Call coul_rfp_forces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        End If

     End If

! accumulate radial distribution functions

     If (l_do_rdf) Call rdf_collect(i,rcut,rrt)

  End Do

! Poisson solver alternative to Ewald

  If (keyfce == 12) Then
     Call poisson_forces(alpha,epsq,engacc,viracc,stress)

     engcpe_rl=engcpe_rl+engacc
     vircpe_rl=vircpe_rl+viracc
  End If

  If (safe) Then
     tmp=0.0_wp
  Else
     tmp=1.0_wp
  End If

! In the case of excluded interactions
! accumulate further radial distribution functions and/or
! calculate Ewald corrections due to short-range exclusions

  If ( lbook .and. (l_do_rdf .or. (keyfce == 2 .or. keyfce == 12)) ) Then

! outer loop over atoms

     Do i=1,natms

! Get list limit

        limit=list(-1,i)-list(0,i)
        If (limit > 0) Then

! calculate interatomic distances

           Do k=1,limit
              j=list(list(0,i)+k,i)

              xxt(k)=xxx(i)-xxx(j)
              yyt(k)=yyy(i)-yyy(j)
              zzt(k)=zzz(i)-zzz(j)
           End Do

! periodic boundary conditions not needed by LC construction
!
!           Call images(imcon,cell,limit,xxt,yyt,zzt)

! square of distances

           Do k=1,limit
              rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
           End Do

! accumulate radial distribution functions

           If (l_do_rdf) Call rdf_excl_collect(i,rcut,rrt)

           If (keyfce == 2) Then ! Ewald corrections
              If (mximpl > 0) Then
                 If (mxompl <= 2) Then
                    Call ewald_excl_mforces_d(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)
                 Else
                    Call ewald_excl_mforces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)
                 End If
              Else
                 Call ewald_excl_forces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress)
              End If

              engcpe_ex=engcpe_ex+engacc
              vircpe_ex=vircpe_ex+viracc
           End If

        End If

     End Do

  End If

! counter for rdf statistics outside loop structures
! and frozen-frozen rdf completeness

  If (l_do_rdf) Then
     If (megfrz /= 0) Then

! outer loop over atoms

        Do i=1,natms

! Get list limit

           limit=list(-2,i)-list(-1,i)
           If (limit > 0) Then

! calculate interatomic distances

              Do k=1,limit
                 j=list(list(0,i)+k,i)

                 xxt(k)=xxx(i)-xxx(j)
                 yyt(k)=yyy(i)-yyy(j)
                 zzt(k)=zzz(i)-zzz(j)
              End Do

! periodic boundary conditions not needed by LC construction
!
!              Call images(imcon,cell,limit,xxt,yyt,zzt)

! square of distances

              Do k=1,limit
                 rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
              End Do

! accumulate radial distribution functions

              Call rdf_frzn_collect(i,rcut,rrt)
           End If

        End Do

     End If

     ncfrdf = ncfrdf + 1
  End If

  Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'two_body_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

! Further Ewald/Poisson Solver corrections or an infrequent refresh

  If (keyfce == 2 .or. keyfce == 12) Then
     If (l_fce) Then

! frozen pairs corrections to coulombic forces

        If (megfrz /= 0) Then
           If (keyfce == 2) Then ! Ewald
              If (mximpl > 0) Then
                 Call ewald_frzn_mforces(rcut,alpha,epsq,engcpe_fr,vircpe_fr,stress)
              Else
                 Call ewald_frzn_forces(rcut,alpha,epsq,engcpe_fr,vircpe_fr,stress)
              End If
           Else !If (keyfce == 12) Then ! Poisson Solver
              Call poisson_frzn_forces(rcut,epsq,engcpe_fr,vircpe_fr,stress)
           End If
        End If

     Else

! Refresh all Ewald k-space contributions

        Call ewald_refresh(engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stress)

     End If

! non-zero total system charge correction (for the whole system)
! ( Fuchs, Proc. R. Soc., A, 151, (585),1935 )

     If (Abs(sumchg) > 1.0e-6_wp) Then
        If (new_nz) Then
           new_nz = .false.
           factor_nz = -0.5_wp * (pi*r4pie0/epsq) * (sumchg/alpha)**2
        End If

        engcpe_nz=factor_nz/volm
        vircpe_nz=-3.0_wp*engcpe_nz
     End If
  End If

! Find the change of energy produced by the torques on multipoles
! under infinitesimal rotations & convert to Cartesian coordinates

  If (mximpl > 0) Call d_ene_trq_mpoles(vircpe_dt,stress)

! sum up contributions to potentials

  If (mxnode > 1) Then
     buffer( 0) = tmp
     buffer( 1) = engkim
     buffer( 2) = virkim
     buffer( 3) = engden
     buffer( 4) = virden
     buffer( 5) = engmet
     buffer( 6) = virmet
     buffer( 7) = engvdw
     buffer( 8) = virvdw
     buffer( 9) = engcpe_rc
     buffer(10) = vircpe_rc
     buffer(11) = engcpe_rl
     buffer(12) = vircpe_rl
     buffer(13) = engcpe_ex
     buffer(14) = vircpe_ex
     buffer(15) = engcpe_fr
     buffer(16) = vircpe_fr
     buffer(17) = vircpe_dt

     Call gsum(buffer(0:17))

     tmp       = buffer( 0)
     engkim    = buffer( 1)
     virkim    = buffer( 2)
     engden    = buffer( 3)
     virden    = buffer( 4)
     engmet    = buffer( 5)
     virmet    = buffer( 6)
     engvdw    = buffer( 7)
     virvdw    = buffer( 8)
     engcpe_rc = buffer( 9)
     vircpe_rc = buffer(10)
     engcpe_rl = buffer(11)
     vircpe_rl = buffer(12)
     engcpe_ex = buffer(13)
     vircpe_ex = buffer(14)
     engcpe_fr = buffer(15)
     vircpe_fr = buffer(16)
     vircpe_dt = buffer(17)
  End If

  safe=(tmp < 0.5_wp)
  If (.not.safe) Call error(505)

! Self-interaction is constant for the default charges only SPME

  If (mxnode > 1 .and. keyfce == 2) Then ! Sum it up for multipolar SPME
     If (mximpl > 0 .and. mxompl <= 2) Call gsum(engsic)
!     If (idnode == 0) Write(nrite,'(1x,a,1p,e18.10 )') 'Self-interaction term: ',engsic
  End If

! Globalise coulombic contributions: cpe

  engcpe = engcpe_rc + engcpe_rl + engcpe_ex + engcpe_fr + engcpe_nz
  vircpe = vircpe_rc + vircpe_rl + vircpe_ex + vircpe_fr + vircpe_nz + vircpe_dt

! Add non-zero total system charge correction to
! diagonal terms of stress tensor (per node)

  tmp = - vircpe_nz/(3.0_wp*Real(mxnode,wp))
  stress(1) = stress(1) + tmp
  stress(5) = stress(5) + tmp
  stress(9) = stress(9) + tmp

! Globalise short-range, KIM and metal interactions with
! their long-range corrections contributions: srp

  engsrp = engkim + (engden + engmet + elrcm(0)) + (engvdw + elrc)
  virsrp = virkim + (virden + virmet + vlrcm(0)) + (virvdw + virlrc)

! Add long-range corrections to diagonal terms of stress tensor (per node)

  tmp = - (virlrc+vlrcm(0))/(3.0_wp*Real(mxnode,wp))
  stress(1) = stress(1) + tmp
  stress(5) = stress(5) + tmp
  stress(9) = stress(9) + tmp

End Subroutine two_body_forces
