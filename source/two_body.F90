Module two_body
  Use kinds, Only : wp
  Use comms,   Only : comms_type,gsum
  Use setup
  Use site, Only : site_type
  Use configuration,  Only : volm,sumchg,natms,xxx,yyy,zzz
  Use neighbours,     Only : neighbours_type,link_cell_pairs
  Use ewald,           Only : ewald_type
  Use mpole,          Only : keyind
  Use coul_spole,     Only : coul_fscp_forces, coul_rfp_forces, coul_cp_forces, coul_dddp_forces
  Use coul_mpole,    Only : coul_fscp_mforces, coul_rfp_mforces, coul_cp_mforces, &
                             coul_dddp_mforces, coul_chrm_forces, d_ene_trq_mpoles
  Use poisson, Only : poisson_type,poisson_forces,poisson_excl_forces,poisson_frzn_forces
  Use vdw,     Only : vdw_type,vdw_forces
  Use metal,   Only : metal_type,metal_forces,metal_ld_compute,metal_lrc
  Use kim
  Use rdfs,    Only : rdf_type,rdf_collect,rdf_excl_collect,rdf_frzn_collect, &
                      rdf_increase_block_number
  Use errors_warnings, Only : error
  Use ewald_spole, Only : ewald_spme_forces,ewald_real_forces,ewald_frzn_forces, ewald_excl_forces
  Use ewald_mpole, Only : ewald_spme_mforces, ewald_real_mforces,ewald_frzn_mforces,ewald_excl_mforces, &
                         ewald_spme_mforces_d,ewald_real_mforces_d,ewald_excl_mforces_d,ewald_excl_mforces

  Use timer,  Only : timer_type,start_timer,stop_timer
  Use development, Only : development_type
  Use statistics, Only : stats_type
  Use core_shell, Only : core_shell_type
  Implicit None

  Private

  Public :: two_body_forces
Contains

Subroutine two_body_forces                        &
           (pdplnc,ensemble,    &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           leql,nsteql,nstep,         &
           cshell,               &
           stats,ewld,devel,met,pois,neigh,site,vdw,rdf,tmr,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating interatomic forces and rdf%rdf
! using the verlet neighbour list
!
! vdw%n_vdw > 0 ------ switch for vdw potentials calculation
! met%n_potentials > 0 ------ switch for metal local density and potentials
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
! author    - i.t.todorov february 2017
! contrib   - h.a.boateng february 2016
! contrib   - p.s.petkov february 2015
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  Logical,                                  Intent( In    ) :: lbook,leql
  Integer,                                  Intent( In    ) :: ensemble,        &
                                                               keyfce,nstfce, &
                                                               megfrz, &
                                                               nsteql,nstep
  Type( core_shell_type ), Intent( InOut ) :: cshell
  Real( Kind = wp ),                        Intent( In    ) :: pdplnc,alpha,epsq
  Type( stats_type ), Intent( InOut )                       :: stats
  Type( ewald_type ),                       Intent( InOut ) :: ewld
  Type( development_type ),                 Intent( In    ) :: devel
  Type( metal_type ),                       Intent( InOut ) :: met
  Type( poisson_type ),                     Intent( InOut ) :: pois
  Type( neighbours_type ),                  Intent( InOut ) :: neigh
  Type( site_type ),                        Intent( In    ) :: site
  Type( vdw_type ),                         Intent( InOut ) :: vdw
  Type( rdf_type ),                         Intent( InOut ) :: rdf
  Type( timer_type ),                       Intent( InOut ) :: tmr
  Type( comms_type ),                       Intent( InOut ) :: comm


  Logical,           Save :: new_nz    = .true.
  Real( Kind = wp ), Save :: factor_nz = 0.0_wp

  Logical           :: safe = .true., l_do_rdf
  Integer           :: fail,i,j,k,limit
  Real( Kind = wp ) :: engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl, &
                       engcpe_ch,vircpe_ch,engcpe_ex,vircpe_ex, &
                       engcpe_fr,vircpe_fr,engcpe_nz,vircpe_nz, &
                       vircpe_dt,                              &
                       engden,virden,engmet,virmet,             &
                       engvdw,virvdw,engkim,virkim,             &
                       engacc,viracc,tmp,buffer(0:19)

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

  Character( Len = 256 ) :: message
  fail=0
  Allocate (xxt(1:neigh%max_list),yyt(1:neigh%max_list),zzt(1:neigh%max_list),rrt(1:neigh%max_list), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'two_body_forces allocation failure'
     Call error(0,message)
  End If

  l_do_rdf = (rdf%l_collect .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,rdf%freq) == 0)

! If k-space SPME is evaluated infrequently check whether
! at this timestep to evaluate or "refresh" with old values.
! At restart allocate the "refresh" arrays and force a fresh
! evaluation.  Repeat the same but only for the SPME k-space
! frozen-frozen evaluations in constant volume ensembles only.

  If (keyfce == 2 .or. keyfce == 12) Then
    Call ewld%check(ensemble,megfrz,nsteql,nstfce,nstep)
  End If

! initialise energy and virial accumulators

  engkim = 0.0_wp
  virkim = 0.0_wp


  engden    = 0.0_wp
  virden    = 0.0_wp

  engmet    = 0.0_wp
  virmet    = 0.0_wp

  engvdw    = 0.0_wp
  virvdw    = 0.0_wp

  stats%engsrp    = 0.0_wp
  stats%virsrp    = 0.0_wp


  engcpe_rc = 0.0_wp
  vircpe_rc = 0.0_wp

  engcpe_rl = 0.0_wp
  vircpe_rl = 0.0_wp

  engcpe_ch = 0.0_wp
  vircpe_ch = 0.0_wp

  engcpe_ex = 0.0_wp
  vircpe_ex = 0.0_wp

  engcpe_fr = 0.0_wp
  vircpe_fr = 0.0_wp

  engcpe_nz = 0.0_wp
  vircpe_nz = 0.0_wp

  vircpe_dt = 0.0_wp

  stats%engcpe    = 0.0_wp
  stats%vircpe    = 0.0_wp

! Set up non-bonded interaction (verlet) list using link cells
  If (neigh%update) Then
    Call link_cell_pairs(vdw%cutoff,met%rcut,pdplnc,lbook,megfrz,cshell,devel,neigh,tmr,comm)
  End If
! Calculate all contributions from KIM

  If (kimim /= ' ') Then
     Call kim_setup(site%ntype_atom,site%unique_atom,site%site_name,kimim,neigh%max_list,comm)
     Call kim_forces(engkim,virkim,stats%stress,neigh%list,comm)
     Call kim_cleanup(comm)
  End If

  If (met%n_potentials > 0) Then

! Reset metal long-range corrections (constant pressure/stress only)

     If (ensemble >= 20) Call metal_lrc(met,site,comm)

! calculate local density in metals

     Call metal_ld_compute(engden,virden,stats%stress,site%ntype_atom,met,neigh,comm)

  End If

! calculate coulombic forces, Ewald sum - fourier contribution
#ifdef CHRONO
  Call start_timer(tmr%t_longrange)
#endif

  If (keyfce == 2 .and. ewld%l_fce) Then
     If (mximpl > 0) Then
        If (mxompl <= 2) Then
           Call ewald_spme_mforces_d(alpha,epsq,engcpe_rc,vircpe_rc,stats%stress,ewld,comm)
        Else
           Call ewald_spme_mforces(alpha,epsq,engcpe_rc,vircpe_rc,stats%stress,ewld,comm)
        End If
     Else
        Call ewald_spme_forces(alpha,epsq,engcpe_rc,vircpe_rc,stats%stress,ewld,comm)
     End If
  End If
#ifdef CHRONO
  Call stop_timer(tmr%t_longrange)
#endif

#ifdef CHRONO
  Call start_timer(tmr%t_shortrange)
#endif
! outer loop over atoms

  Do i=1,natms

! Get neigh%list limit

     limit=neigh%list(0,i)

! calculate interatomic distances

     Do k=1,limit
        j=neigh%list(k,i)

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

     If (met%n_potentials > 0) Then
        Call metal_forces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,safe,site%ntype_atom,met,neigh)

        engmet=engmet+engacc
        virmet=virmet+viracc
     End If

! calculate short-range force and potential terms

     If (vdw%n_vdw > 0) Then
        Call vdw_forces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,vdw)

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
              Call ewald_real_mforces_d(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,ewld,neigh,comm)
           Else
              Call ewald_real_mforces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,comm)
           End If

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 4) Then

! distance dependant dielectric potential

           Call coul_dddp_mforces(i,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 6) Then

! coulombic 1/r potential with no truncation or damping

           Call coul_cp_mforces(i,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 8) Then

! force-shifted coulomb potentials

           Call coul_fscp_mforces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,comm)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 10) Then

! reaction field potential

           Call coul_rfp_mforces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,comm)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        End If

     Else

        If (keyfce == 2) Then

! calculate coulombic forces, Ewald sum - real space contribution

           Call ewald_real_forces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,comm)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 4) Then

! distance dependant dielectric potential

           Call coul_dddp_forces(i,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 6) Then

! coulombic 1/r potential with no truncation or damping

           Call coul_cp_forces(i,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 8) Then

! force-shifted coulomb potentials

           Call coul_fscp_forces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,comm)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        Else If (keyfce == 10) Then

! reaction field potential

           Call coul_rfp_forces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,comm)

           engcpe_rl=engcpe_rl+engacc
           vircpe_rl=vircpe_rl+viracc

        End If

     End If

! accumulate radial distribution functions

     If (l_do_rdf) Call rdf_collect(i,rrt,neigh,rdf)

  End Do

! Poisson solver alternative to Ewald

  If (keyfce == 12) Then
     Call poisson_forces(alpha,epsq,engacc,viracc,stats%stress,pois,comm)

     engcpe_rl=engcpe_rl+engacc
     vircpe_rl=vircpe_rl+viracc
  End If

! metal potential safety

  If (safe) Then
     tmp=0.0_wp
  Else
     tmp=1.0_wp
  End If

! in the case of bonded interactions 3 possible subcases
! cases for excluded interactions:
! (1) RDF accumulate further the short-range exclusions
! (2) Ewald corrections due to short-range exclusions
! (3) CHARMM core-shell self-induction additions

  If ( lbook .and. (l_do_rdf .or. (keyfce == 2 .or. keyfce == 12) .or. keyind == 1) ) Then
     Do i=1,natms ! outer loop over atoms
        limit=neigh%list(-1,i)-neigh%list(0,i) ! Get neigh%list limit
        If (limit > 0) Then

! calculate interatomic distances

           Do k=1,limit
              j=neigh%list(neigh%list(0,i)+k,i)

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

           If (l_do_rdf) Call rdf_excl_collect(i,rrt,neigh,rdf)

           If (keyfce == 2) Then ! Ewald corrections
              If (mximpl > 0) Then
                 If (mxompl <= 2) Then
                    Call ewald_excl_mforces_d(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)
                 Else
                    Call ewald_excl_mforces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)
                 End If
              Else
                 Call ewald_excl_forces(i,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)
              End If

              engcpe_ex=engcpe_ex+engacc
              vircpe_ex=vircpe_ex+viracc
           End If

! get CHARMM core-shell self-induction contributions

           If (keyind == 1) Then
              If (neigh%list(-3,i)-neigh%list(0,i) > 0) Then
                 Call coul_chrm_forces(i,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh)

                 engcpe_ch=engcpe_ch+engacc
                 vircpe_ch=vircpe_ch+viracc
              End If
           End If

        End If
     End Do
  End If

! counter for rdf%rdf statistics outside loop structures
! and frozen-frozen rdf%rdf completeness

  If (l_do_rdf) Then
     If (megfrz /= 0) Then

! outer loop over atoms

        Do i=1,natms

! Get neigh%list limit

           limit=neigh%list(-2,i)-neigh%list(-1,i)
           If (limit > 0) Then

! calculate interatomic distances

              Do k=1,limit
                 j=neigh%list(neigh%list(-1,i)+k,i)

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

              Call rdf_frzn_collect(i,rrt,neigh,rdf)
           End If

        End Do

     End If

     rdf%n_configs = rdf%n_configs + 1
  End If

  Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'two_body_forces deallocation failure'
     Call error(0,message)
  End If

!Increase rdf%block_number when required
  If (l_do_rdf) Then
    Call rdf_increase_block_number(rdf,nstep)
  End If

! Further Ewald/Poisson Solver corrections or an infrequent refresh

  If (keyfce == 2 .or. keyfce == 12) Then
     If (ewld%l_fce) Then

! frozen pairs corrections to coulombic forces

        If (megfrz /= 0) Then
           If (keyfce == 2) Then ! Ewald
              If (mximpl > 0) Then
                 Call ewald_frzn_mforces(alpha,epsq,engcpe_fr,vircpe_fr,stats%stress,ewld,neigh,comm)
              Else
                 Call ewald_frzn_forces(alpha,epsq,engcpe_fr,vircpe_fr,stats%stress,ewld,neigh,comm)
              End If
           Else !If (keyfce == 12) Then ! Poisson Solver
              Call poisson_frzn_forces(epsq,engcpe_fr,vircpe_fr,stats%stress,ewld,neigh,comm)
           End If
        End If

     Else

! Refresh all Ewald k-space contributions

        Call ewld%refresh(engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stats%stress)

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

  If (mximpl > 0) Call d_ene_trq_mpoles(vircpe_dt,stats%stress)

! sum up contributions to potentials


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
     buffer(13) = engcpe_ch
     buffer(14) = vircpe_ch
     buffer(15) = engcpe_ex
     buffer(16) = vircpe_ex
     buffer(17) = engcpe_fr
     buffer(18) = vircpe_fr
     buffer(19) = vircpe_dt

     Call gsum(comm,buffer(0:17))

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
     engcpe_ch = buffer(13)
     vircpe_ch = buffer(14)
     engcpe_ex = buffer(15)
     vircpe_ex = buffer(16)
     engcpe_fr = buffer(17)
     vircpe_fr = buffer(18)
     vircpe_dt = buffer(19)


  safe=(tmp < 0.5_wp)
  If (.not.safe) Call error(505)

! Self-interaction is constant for the default charges only SPME

  If (keyfce == 2) Then ! Sum it up for multipolar SPME
     If (mximpl > 0 .and. mxompl <= 2) Call gsum(comm,ewld%engsic)
     !Write(message,'(a,1p,e18.10)') 'Self-interaction term: ',engsic
     !Call info(message,.true.)
  End If

! Globalise coulombic contributions: cpe

  stats%engcpe = engcpe_rc + engcpe_rl + engcpe_ch + engcpe_ex + engcpe_fr + engcpe_nz
  stats%vircpe = vircpe_rc + vircpe_rl + vircpe_ch + vircpe_ex + vircpe_fr + vircpe_nz + vircpe_dt

! Add non-zero total system charge correction to
! diagonal terms of stress tensor (per node)

  tmp = - vircpe_nz/(3.0_wp*Real(comm%mxnode,wp))
  stats%stress(1) = stats%stress(1) + tmp
  stats%stress(5) = stats%stress(5) + tmp
  stats%stress(9) = stats%stress(9) + tmp

! Globalise short-range, KIM and metal interactions with
! their long-range corrections contributions: srp

  stats%engsrp = engkim + (engden + engmet + met%elrc(0)) + (engvdw + vdw%elrc)
  stats%virsrp = virkim + (virden + virmet + met%vlrc(0)) + (virvdw + vdw%vlrc)

! Add long-range corrections to diagonal terms of stress tensor (per node)

  tmp = - (vdw%vlrc+met%vlrc(0))/(3.0_wp*Real(comm%mxnode,wp))
  stats%stress(1) = stats%stress(1) + tmp
  stats%stress(5) = stats%stress(5) + tmp
  stats%stress(9) = stats%stress(9) + tmp

#ifdef CHRONO
  Call stop_timer(tmr%t_shortrange)
#endif

End Subroutine two_body_forces
End Module two_body
