Module drivers
  Use angles,               Only: angles_forces,&
                                  angles_type
  Use angular_distribution, Only: adf_calculate,&
                                  adf_type
  Use bonds,                Only: bonds_forces,&
                                  bonds_type
  Use build_book,           Only: build_book_intra
  Use build_excl,           Only: build_excl_intra
  Use comms,                Only: comms_type,&
                                  gmax,&
                                  gsum,&
                                  gtime,&
                                  root_id, gsync
  Use configuration,        Only: check_config,&
                                  configuration_type,&
                                  freeze_atoms,&
                                  coordinate_buffer_type,&
                                  len_atmnam,&
                                  gather_coordinates,&
                                  gather_atomic_names,&
                                  distribute_forces,&
                                  gather_forces
  Use constants,            Only: boltz
  Use constraints,          Only: constraints_quench,&
                                  constraints_type
  Use coord,                Only: checkcoord,&
                                  coord_type,&
                                  init_coord_list
  Use core_shell,           Only: SHELL_ADIABATIC,&
                                  SHELL_RELAXED,&
                                  core_shell_forces,&
                                  core_shell_kinetic,&
                                  core_shell_on_top,&
                                  core_shell_quench,&
                                  core_shell_relax,&
                                  core_shell_type
  Use defects,              Only: defects_type,&
                                  defects_write
  Use deport_data,          Only: mpoles_rotmat_set_halo,&
                                  relocate_particles
  Use development,          Only: development_type
#ifdef DFTBP
  Use dftb_library,         Only: run_dftbplus,&
                                  dftb_geometry_type,&
                                  convert_unit, &
                                  output_dftb_forces
#endif
  Use dihedrals,            Only: dihedrals_forces,&
                                  dihedrals_type
  Use domains,              Only: domains_type
  Use dpd,                  Only: dpd_thermostat
  Use electrostatic,        Only: ELECTROSTATIC_NULL,&
                                  electrostatic_type
  Use errors_warnings,      Only: error,&
                                  info,&
                                  warning
  Use evb,                 Only : evb_type,&
                                  evb_pes,&
                                  read_evb_settings,&
                                  evb_merge_stochastic, &
                                  evb_setzero,&
                                  evb_population,&
                                  evb_check_external,&
                                  evb_check_intermolecular,&
                                  evb_check_intramolecular,&
                                  evb_check_constraints,&
                                  evb_check_configs, &
                                  evb_check_intrinsic, &
                                  evb_check_vdw, &
                                  evb_prevent
  Use ewald,                Only: ewald_type
  Use external_field,       Only: FIELD_NULL,&
                                  external_field_apply,&
                                  external_field_correct,&
                                  external_field_type
  Use filename,             Only: FILE_HISTORF,&
                                  FILE_HISTORY,&
                                  FILE_POPEVB,&
                                  FILE_OUTPUT,&
                                  file_type
  Use flow_control,         Only: RESTART_KEY_CLEAN,&
                                  RESTART_KEY_OLD,&
                                  flow_type,&
                                  DFTB
  Use four_body,            Only: four_body_forces,&
                                  four_body_type
  Use greenkubo,            Only: greenkubo_type,&
                                  vaf_collect,&
                                  vaf_write
  Use halo,                 Only: refresh_halo_positions,&
                                  set_halo_particles
  Use impacts,              Only: impact,&
                                  impact_type
  Use inversions,           Only: inversions_forces,&
                                  inversions_type
  Use io,                   Only: io_type
  Use kim,                  Only: kim_type
  Use kinds,                Only: wi,&
                                  wp
  Use kinetics,             Only: cap_forces,&
                                  getknr,&
                                  kinstresf,&
                                  kinstress,&
                                  kinstrest
  Use langevin,             Only: langevin_forces
  Use metal,                Only: metal_type
  Use minimise,             Only: minimise_relax,&
                                  minimise_type,&
                                  zero_k_optimise
  Use mpole,                Only: mpole_type
  Use msd,                  Only: msd_type,&
                                  msd_write
  Use neighbours,           Only: link_cell_pairs,&
                                  neighbours_type,&
                                  vnl_check
  Use netcdf_wrap,          Only: netcdf_param
  Use npt_berendsen,        Only: npt_b0_vv,&
                                  npt_b1_vv
  Use npt_langevin,         Only: npt_l0_vv,&
                                  npt_l1_vv
  Use npt_mtk,              Only: npt_m0_vv,&
                                  npt_m1_vv
  Use npt_nose_hoover,      Only: npt_h0_vv,&
                                  npt_h1_vv
  Use nst_berendsen,        Only: nst_b0_vv,&
                                  nst_b1_vv
  Use nst_langevin,         Only: nst_l0_vv,&
                                  nst_l1_vv
  Use nst_mtk,              Only: nst_m0_vv,&
                                  nst_m1_vv
  Use nst_nose_hoover,      Only: nst_h0_vv,&
                                  nst_h1_vv
  Use numerics,             Only: seed_type
  Use nve,                  Only: nve_0_vv,&
                                  nve_1_vv
  Use nvt_anderson,         Only: nvt_a0_vv,&
                                  nvt_a1_vv
  Use nvt_berendsen,        Only: nvt_b0_vv,&
                                  nvt_b1_vv
  Use nvt_ekin,             Only: nvt_e0_vv,&
                                  nvt_e1_vv
  Use nvt_gst,              Only: nvt_g0_vv,&
                                  nvt_g1_vv
  Use nvt_langevin,         Only: nvt_l0_vv,&
                                  nvt_l1_vv,&
                                  nvt_l2_vv
  Use nvt_nose_hoover,      Only: nvt_h0_vv,&
                                  nvt_h1_vv
  Use parse,                Only: lower_case,&
                                  strip_blanks
  Use plumed,               Only: plumed_apply,&
                                  plumed_type
  Use pmf,                  Only: pmf_quench,&
                                  pmf_type
  Use poisson,              Only: poisson_type
  Use rdfs,                 Only: rdf_type
  Use rigid_bodies,         Only: rigid_bodies_coms,&
                                  rigid_bodies_quench,&
                                  rigid_bodies_stress,&
                                  rigid_bodies_tags,&
                                  rigid_bodies_type,&
                                  xscale
  Use rsds,                 Only: rsd_type,&
                                  rsd_write
  Use shared_units,         Only: SHARED_UNIT_UPDATE_FORCES,&
                                  update_shared_units
  Use site,                 Only: site_type
  Use statistics,           Only: calculate_heat_flux,&
                                  statistics_collect,&
                                  statistics_connect_frames,&
                                  statistics_connect_set,&
                                  stats_type,&
                                  write_per_part_contribs
  Use stochastic_boundary,  Only: stochastic_boundary_vv
  Use system,               Only: system_revive
  Use temperature,          Only: regauss_temperature,&
                                  scale_temperature,&
                                  set_temperature
  Use tersoff,              Only: tersoff_forces,&
                                  tersoff_type
  Use tethers,              Only: tethers_forces,&
                                  tethers_type
  Use thermostat,           Only: &
                                  DPD_NULL, DPD_SECOND_ORDER, ENS_NPT_BERENDSEN, &
                                  ENS_NPT_BERENDSEN_ANISO, ENS_NPT_LANGEVIN, &
                                  ENS_NPT_LANGEVIN_ANISO, ENS_NPT_MTK, ENS_NPT_MTK_ANISO, &
                                  ENS_NPT_NOSE_HOOVER, ENS_NPT_NOSE_HOOVER_ANISO, ENS_NVE, &
                                  ENS_NVT_ANDERSON, ENS_NVT_BERENDSEN, ENS_NVT_EVANS, &
                                  ENS_NVT_GENTLE, ENS_NVT_LANGEVIN, ENS_NVT_LANGEVIN_INHOMO, &
                                  ENS_NVT_NOSE_HOOVER, VV_FIRST_STAGE, VV_SECOND_STAGE, &
                                  thermostat_type
  Use three_body,           Only: three_body_forces,&
                                  threebody_type
  Use timer,                Only: start_timer,&
                                  stop_timer,&
                                  timer_type
  Use trajectory,           Only: read_history,&
                                  trajectory_type,&
                                  trajectory_write
  Use ttm,                  Only: ttm_type
  Use ttm_track,            Only: ttm_ion_temperature,&
                                  ttm_thermal_diffusion
  Use two_body,             Only: two_body_forces
  Use vdw,                  Only: vdw_type
  Use z_density,            Only: z_density_type

#ifdef EXPERIMENT
  Use meminfo, Only : mem_picture
#endif

  Implicit None
  Private

  Public :: md_vv
  Public :: replay_historf
  Public :: replay_history
  Public :: calculate_dftb_forces

Contains

  Subroutine impact_option(levcfg, nstep, nsteql, rigid, cshell, stats, impa, config, comm)

    Integer(Kind=wi),         Intent(InOut) :: levcfg, nstep, nsteql
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(stats_type),         Intent(InOut) :: stats
    Type(impact_type),        Intent(InOut) :: impa
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256) :: messages(6)

    !!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Apply impact
    ! levcfg == 2 avoids application twice when tmd happens at (re)start for VV

    If (nstep == impa%tmd .and. levcfg == 2) Then
      Write (messages(1), '(a)') ''
      Write (messages(2), '(a)') 'initiating IMPACT:'
      Write (messages(3), '(a,i10)') 'particle (index): ', impa%imd
      Write (messages(4), '(a,i10)') 'timestep (steps): ', impa%tmd
      Write (messages(5), '(a,1p,e12.5)') 'energy   (keV):   ', impa%emd
      Write (messages(6), '(a,1p,3e12.4)') 'v-r(x,y,z):       ', impa%vmx, impa%vmy, impa%vmz
      Call info(messages, 6, .true.)

      If (nstep + 1 <= nsteql) Call warning(380, Real(nsteql, wp), 0.0_wp, 0.0_wp)

      Call impact(rigid, cshell, impa, config, comm)

      ! Correct kinetic stress and energy

      If (rigid%total > 0) Then
        Call kinstresf(stats%strknf, config, comm)
        Call kinstrest(rigid, stats%strknt, comm)

        stats%strkin = stats%strknf + stats%strknt

        stats%engrot = getknr(rigid, comm)
      Else
        Call kinstress(stats%strkin, config, comm)
      End If
      stats%engke = 0.5_wp * (stats%strkin(1) + stats%strkin(5) + stats%strkin(9))
    End If

    !!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
  End Subroutine impact_option

  Subroutine calculate_forces_for_historf(cnfig, flow, io, cshell, cons, pmf, stat, plume, pois, bond, angle, dihedral, &
                              inversion, tether, threebody, neigh, sites, vdws, tersoffs, fourbody, rdf, netcdf, &
                              minim, mpoles, ext_field, rigid, electro, domain, kim_data, msd_data, tmr, files, &
                              green, devel, ewld, met, seed, thermo, crd, comm)

    Type(configuration_type),  Intent(InOut) :: cnfig
    Type(flow_type),           Intent(InOut) :: flow
    Type(io_type),             Intent(InOut) :: io
    Type(core_shell_type),     Intent(InOut) :: cshell
    Type(constraints_type),    Intent(InOut) :: cons
    Type(pmf_type),            Intent(InOut) :: pmf
    Type(stats_type),          Intent(InOut) :: stat
    Type(plumed_type),         Intent(InOut) :: plume
    Type(poisson_type),        Intent(InOut) :: pois
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(site_type),           Intent(InOut) :: sites
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(netcdf_param),        Intent(In   ) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(In   ) :: domain
    Type(kim_type),            Intent(InOut) :: kim_data
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(timer_type),          Intent(InOut) :: tmr
    Type(file_type),           Intent(InOut) :: files(:)
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(development_type),    Intent(InOut) :: devel
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(metal_type),          Intent(InOut) :: met
    Type(seed_type),           Intent(InOut) :: seed
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(coord_type),          Intent(InOut) :: crd
    Type(comms_type),          Intent(InOut) :: comm

    Integer          :: i
    Integer(Kind=wi) :: switch
    Logical          :: ltmp

    !!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! for new simulations when using the relaxed shell model
    ! set shells on top of their cores preventatively

    If ((cshell%megshl > 0 .and. cshell%keyshl == SHELL_RELAXED) .and. &
        (flow%restart_key == RESTART_KEY_CLEAN .and. flow%step == 0 .and. flow%equil_steps > 0)) Then
      Call core_shell_on_top(cshell, cnfig, comm)

      ! Refresh mappings

      Call refresh_mappings(cnfig, flow, cshell, cons, pmf, stat, msd_data, bond, angle, &
                            dihedral, inversion, tether, neigh, sites, mpoles, rigid, &
                            domain, kim_data, ewld, green, &
                            minim, thermo, electro, crd, comm, tmr)
    End If

    100 Continue ! Only used when relaxed is false

    ! Initialise force arrays and possible torques for multipolar electrostatics
    ! and stress tensor (these are all additive in the force subroutines)
    Do i = 1, cnfig%mxatms
      cnfig%parts(i)%fxx = 0.0_wp
      cnfig%parts(i)%fyy = 0.0_wp
      cnfig%parts(i)%fzz = 0.0_wp
    End Do
    If (mpoles%max_mpoles > 0) Then
      mpoles%torque_x = 0.0_wp; mpoles%torque_y = 0.0_wp; mpoles%torque_z = 0.0_wp
    End If

    stat%stress = 0.0_wp

    ! Initialise variables for two_body interaction (including long range, which is not strictly a two_body interaction problem)
    stat%engsrp = 0.0_wp
    stat%virsrp = 0.0_wp
    stat%engcpe = 0.0_wp
    stat%vircpe = 0.0_wp

    ! Set up non-bonded interaction (verlet) list using link cells
    If ((.not. (met%max_metal == 0 .and. electro%key == ELECTROSTATIC_NULL .and. &
                vdws%no_vdw .and. rdf%max_rdf == 0) .or. kim_data%active) .and. neigh%update) Then
      Call link_cell_pairs(vdws%cutoff, met%rcut, flow%book, cnfig%megfrz, cshell, devel, &
                           neigh, mpoles, domain, tmr, cnfig, comm)
    End If

    ! Calculate tersoff forces

    If (tersoffs%n_potential > 0) Call tersoff_forces(tersoffs, stat, neigh, domain, cnfig, comm)

    ! Calculate three-body forces

    If (threebody%ntptbp > 0) Call three_body_forces(stat, threebody, neigh, domain, cnfig, comm)

    ! Calculate four-body forces

    If (fourbody%n_potential > 0) Call four_body_forces(fourbody, stat, neigh, domain, cnfig, comm)

#ifdef CHRONO
    Call start_timer(tmr, 'Bonded Forces')
#endif

    ! Calculate shell model forces

    If (cshell%megshl > 0) Call core_shell_forces(cshell, stat, cnfig, comm)

    ! Calculate tethered atom forces

    If (tether%total > 0) Call tethers_forces(stat, tether, cnfig, comm)

    ! Calculate bond forces

    If (bond%total > 0) Then
      ltmp = (bond%bin_pdf > 0 .and. &
              ((.not. flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
              Mod(flow%step, flow%freq_bond) == 0)

      if (ltmp) then
        switch = 2
      else
        switch = 1
      end if

      Call bonds_forces(switch, stat%engbnd, stat%virbnd, stat%stress, neigh%cutoff, &
                        stat%engcpe, stat%vircpe, bond, mpoles, electro, cnfig, comm)
    End If

    ! Calculate valence angle forces

    If (angle%total > 0) Then
      ltmp = (angle%bin_adf > 0 .and. &
              ((.not. flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
              Mod(flow%step, flow%freq_angle) == 0)

      switch = 1 + Merge(1, 0, ltmp)
      Call angles_forces(switch, stat%engang, stat%virang, stat%stress, angle, cnfig, comm)
    End If

    ! Calculate dihedral forces

    If (dihedral%total > 0) Then
      ltmp = (dihedral%bin_adf > 0 .and. &
              ((.not. flow%equilibration) .or. flow%step >= flow%equil_steps) &
              .and. Mod(flow%step, flow%freq_dihedral) == 0)

      switch = 1 + Merge(1, 0, ltmp)
      Call dihedrals_forces(switch, stat%engdih, stat%virdih, stat%stress, &
                            neigh%cutoff, stat%engcpe, stat%vircpe, stat%engsrp, &
                            stat%virsrp, dihedral, vdws, mpoles, electro, cnfig, comm)
    End If

    ! Calculate inversion forces

    If (inversion%total > 0) Then
      ltmp = (inversion%bin_adf > 0 .and. &
              ((.not. flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
              Mod(flow%step, flow%freq_inversion) == 0)

      switch = 1 + Merge(1, 0, ltmp)
      Call inversions_forces(switch, stat%enginv, stat%virinv, stat%stress, inversion, cnfig, comm)
    End If

#ifdef CHRONO
    Call stop_timer(tmr, 'Bonded Forces')
#endif

    ! Calculate pair-like forces (metal,vdws,electrostatic) and add lrc

    If (.not. (met%max_metal == 0 .and. electro%key == ELECTROSTATIC_NULL .and. &
               vdws%no_vdw .and. rdf%max_rdf == 0) .or. kim_data%active) Then
      Call two_body_forces(thermo%ensemble, flow%book, cnfig%megfrz, &
                           flow%equilibration, flow%equil_steps, flow%step, stat, ewld, met, pois, neigh, sites, vdws, rdf, &
                           mpoles, electro, domain, tmr, kim_data, cnfig, comm)
    End If

    ! Configurational energy
    stat%stpcfg = stat%engcpe + stat%engsrp + stat%engter + stat%engtbp + stat%engfbp + &
                  stat%engshl + stat%engtet + stat%engbnd + stat%engang + stat%engdih + stat%enginv

    ! Apply external field

    If (ext_field%key /= FIELD_NULL) Then
      Call external_field_apply(flow%time, flow%equilibration, flow%equil_steps, flow%step, cshell, stat, rdf, &
                                ext_field, rigid, domain, cnfig, comm)

      ! Add external energy contribution
      stat%stpcfg = stat%stpcfg + stat%engfld
    End If


    ! Apply PLUMED driven dynamics

    If (plume%l_plumed) Then
      Call plumed_apply(cnfig,flow%run_steps,flow%step,stat,plume,comm)
    End If
    ! Apply pseudo thermostat - force cycle (0)

    If (thermo%l_stochastic_boundaries) Then
      Call stochastic_boundary_vv(0, thermo%tstep, flow%step, sites%dof_site, cshell, &
                                  stat, thermo, rigid, domain, cnfig, seed, comm)
    End If

    ! Cap forces in equilibration mode

    If (flow%step <= flow%equil_steps .and. flow%force_cap) Call cap_forces(thermo%temp, cnfig, comm)

    ! Frozen atoms option

    Call freeze_atoms(cnfig)

    ! Minimisation option and Relaxed shell model optimisation

    If (flow%simulation .and. (minim%minimise .or. cshell%keyshl == SHELL_RELAXED)) Then

      If (cshell%keyshl == SHELL_RELAXED) Then
        Call core_shell_relax(flow%strict, rdf%l_collect, &
                              stat%stpcfg, cshell, stat, domain, cnfig, files, comm)
      End If

      If (.not. cshell%relaxed) Go To 200 ! Shells relaxation takes priority over minimisation

      If (minim%minimise .and. flow%step >= 0 .and. flow%step <= flow%run_steps .and. flow%step <= flow%equil_steps) Then
        If (minim%freq == 0 .and. flow%step == 0) Then
          Call minimise_relax(flow%strict .or. cshell%keyshl == SHELL_RELAXED, &
                              rdf%l_collect, thermo%tstep, stat%stpcfg, io, stat, pmf, cons, &
                              netcdf, minim, rigid, domain, cnfig, files, comm)
        Else If (minim%freq > 0 .and. flow%step > 0) Then
          If (Mod(flow%step - flow%equil_steps, minim%freq) == 0) Then
            Call minimise_relax(flow%strict .or. cshell%keyshl == SHELL_RELAXED, &
                                rdf%l_collect, thermo%tstep, stat%stpcfg, io, stat, pmf, cons, &
                                netcdf, minim, rigid, domain, cnfig, files, comm)
          End If
        End If
      End If

      200 Continue

      ! Refresh mappings

      If (.not. (cshell%relaxed .and. minim%relaxed)) Then
        Call refresh_mappings(cnfig, flow, cshell, cons, pmf, stat, msd_data, bond, angle, &
                              dihedral, inversion, tether, neigh, sites, mpoles, rigid, domain, &
                              kim_data, ewld, green, minim, thermo, &
                              electro, crd, comm, tmr)
        Go To 100
      End If
    End If

    ! Get RB COM stress and virial at restart only - also available at w_at_start_vv for cnfig%levcfg==2

    If (flow%newjob) Then
      If (rigid%total > 0) Then
        If (thermo%l_langevin) Then
          Call langevin_forces(flow%step, thermo%temp, thermo%tstep, thermo%chi, &
                               thermo%fxl, thermo%fyl, thermo%fzl, cshell, cnfig, seed)
          If (rigid%share) Then
            Call update_shared_units(cnfig, rigid%list_shared, &
                                     rigid%map_shared, thermo%fxl, thermo%fyl, thermo%fzl, domain, comm)
          End If
          Call rigid_bodies_stress(stat%strcom, cnfig, rigid, comm, thermo%fxl, thermo%fyl, thermo%fzl)
        Else
          Call rigid_bodies_stress(stat%strcom, rigid, cnfig, comm)
        End If
        stat%vircom = -(stat%strcom(1) + stat%strcom(5) + stat%strcom(9))
      End If
    End If

    ! Total virial (excluding constraint, PMF and RB COM virials for npt routines)
    ! Total stress (excluding constraint, PMF, RB COM and kinetic stress for npt routines)
    !
    ! NOTE(1):  virsrp already includes vdws%vlrc and vlrcm(0) and so
    !           does the stress diagonal elements (by minus a third),
    !           engsrp includes vdws%elrc and elrcm(0)
    !
    ! NOTE(2):  virfbp, virinv and virdih are allegedly always zero

    stat%virtot = stat%vircpe + stat%virsrp + stat%virter + stat%virtbp + stat%virfbp + &
                  stat%virshl + stat%virtet + stat%virbnd + stat%virang + stat%virdih + &
                  stat%virinv + stat%virfld

    Call gsum(comm, stat%stress)

    ! If RBs are present update forces on shared ones

    If (rigid%share) Then
      Call update_shared_units(cnfig, rigid%list_shared, &
                               rigid%map_shared, SHARED_UNIT_UPDATE_FORCES, domain, comm)
    End If

    !!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine calculate_forces_for_historf



  Subroutine calculate_forces_evb(evbff, cnfig, flow, io, cshell, cons, pmf, stat, plume, pois, bond, angle, dihedral,&
                                   inversion, tether, threebody, neigh, sites, vdws, tersoffs, fourbody, rdf, netcdf, &
                                   minim, mpoles, ext_field, rigid, electro, domain, kim_data, msd_data, tmr, files,&
                                   green, devel, ewld, met, seed, thermo, crd, comm)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! dl_poly_4 subroutine to calculate energies, forces, stress and virials for a given ionic configuration
   ! This subroutine is an extension of of the former subroutine calculate_forces that allows coupling up to
   ! three FFs are coupled via EVB, and it also works for a single FF
   !
   ! copyright - daresbury laboratory
   ! contribution - i.scivetti January 2020
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),            Intent(InOut) :: evbff
    Type(configuration_type),  Intent(InOut) :: cnfig(:)
    Type(io_type),             Intent(InOut) :: io
    Type(flow_type),           Intent(InOut) :: flow
    Type(core_shell_type),     Intent(InOut) :: cshell(:)
    Type(constraints_type),    Intent(InOut) :: cons(:)
    Type(pmf_type),            Intent(InOut) :: pmf(:)
    Type(stats_type),          Intent(InOut) :: stat(:)
    Type(plumed_type),         Intent(InOut) :: plume(:)
    Type(poisson_type),        Intent(InOut) :: pois(:)
    Type(bonds_type),          Intent(InOut) :: bond(:)
    Type(angles_type),         Intent(InOut) :: angle(:)
    Type(dihedrals_type),      Intent(InOut) :: dihedral(:)
    Type(inversions_type),     Intent(InOut) :: inversion(:)
    Type(tethers_type),        Intent(InOut) :: tether(:)
    Type(threebody_type),      Intent(InOut) :: threebody(:)
    Type(neighbours_type),     Intent(InOut) :: neigh(:)
    Type(site_type),           Intent(InOut) :: sites(:)
    Type(vdw_type),            Intent(InOut) :: vdws(:)
    Type(tersoff_type),        Intent(InOut) :: tersoffs(:)
    Type(four_body_type),      Intent(InOut) :: fourbody(:)
    Type(rdf_type),            Intent(InOut) :: rdf(:)
    Type(netcdf_param),        Intent(In   ) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim(:)
    Type(mpole_type),          Intent(InOut) :: mpoles(:)
    Type(external_field_type), Intent(InOut) :: ext_field(:)
    Type(rigid_bodies_type),   Intent(InOut) :: rigid(:)
    Type(electrostatic_type),  Intent(InOut) :: electro(:)
    Type(domains_type),        Intent(In   ) :: domain(:)
    Type(kim_type),            Intent(InOut) :: kim_data(:)
    Type(msd_type),            Intent(InOut) :: msd_data(:)
    Type(timer_type),          Intent(InOut) :: tmr
    Type(file_type),           Intent(InOut) :: files(:)
    Type(greenkubo_type),      Intent(InOut) :: green(:)
    Type(development_type),    Intent(InOut) :: devel
    Type(ewald_type),          Intent(InOut) :: ewld(:)
    Type(metal_type),          Intent(InOut) :: met(:)
    Type(seed_type),           Intent(InOut) :: seed
    Type(thermostat_type),     Intent(InOut) :: thermo(:)
    Type(coord_type),          Intent(InOut) :: crd(:)
    Type(comms_type)    ,      Intent(InOut) :: comm

    Integer          :: i, ff
    Integer(Kind=wi) :: switch
    Logical          :: ltmp


    ! for new simulations when using the relaxed shell model
    ! set shells on top of their cores preventatively

     Do ff = 1, flow%NUM_FF
       If ( (cshell(ff)%megshl > 0 .and. cshell(ff)%keyshl == SHELL_RELAXED) .and. &
         (flow%restart_key == RESTART_KEY_CLEAN .and. flow%step == 0 .and. flow%equil_steps > 0) ) Then
         Call core_shell_on_top(cshell(ff), cnfig(ff), comm)

      ! Refresh mappings

         Call refresh_mappings(cnfig(ff), flow,cshell(ff), cons(ff), pmf(ff), stat(ff), msd_data(ff), bond(ff), angle(ff), &
          dihedral(ff),inversion(ff),tether(ff),neigh(ff),sites(ff),mpoles(ff),rigid(ff),domain(ff),kim_data(ff), &
          ewld(ff),green(ff),minim(ff),thermo(ff),electro(ff), crd(ff), comm, tmr)
       End If

     End Do


    100  Continue ! Only used when relaxed is false

    ! Initialise force arrays and possible torques for multipolar electrostatics,
    ! stress tensor (these are all additive in the force subroutines) and torques
    Do ff = 1, flow%NUM_FF
      Do i =1, cnfig(ff)%mxatms
        cnfig(ff)%parts(i)%fxx = 0.0_wp
        cnfig(ff)%parts(i)%fyy = 0.0_wp
        cnfig(ff)%parts(i)%fzz = 0.0_wp
      End Do
      stat(ff)%stress = 0.0_wp
      If (mpoles(ff)%max_mpoles > 0) Then
        mpoles(ff)%torque_x=0.0_wp ; mpoles(ff)%torque_y=0.0_wp ; mpoles(ff)%torque_z=0.0_wp
      End If
    End Do

    ! Initialise variables for two_body interaction (including long range, which is not strictly a two_body interaction problem)
    stat%engsrp = 0.0_wp
    stat%virsrp = 0.0_wp
    stat%engcpe = 0.0_wp
    stat%vircpe = 0.0_wp

    Do ff = 1, flow%NUM_FF
    ! Set up non-bonded interaction (verlet) list using link cells
      If ((.not.(met(ff)%max_metal == 0 .and. electro(ff)%key == ELECTROSTATIC_NULL .and. &
          vdws(ff)%no_vdw .and. rdf(ff)%max_rdf == 0) .or. kim_data(ff)%active) .and. neigh(ff)%update) Then
         Call link_cell_pairs(vdws(ff)%cutoff, met(ff)%rcut, flow%book, cnfig(ff)%megfrz, cshell(ff), devel, &
                              neigh(ff), mpoles(ff), domain(ff), tmr, cnfig(ff), comm)
      End If
    End Do

    Do ff = 1, flow%NUM_FF
      ! Calculate tersoff forces

      If (tersoffs(ff)%n_potential > 0) Call tersoff_forces(tersoffs(ff), stat(ff), neigh(ff), domain(ff), cnfig(ff), comm)

      ! Calculate three-body forces

      If (threebody(ff)%ntptbp > 0) Call three_body_forces(stat(ff), threebody(ff), neigh(ff), domain(ff), cnfig(ff), comm)

      ! Calculate four-body forces

      If (fourbody(ff)%n_potential > 0) Call four_body_forces(fourbody(ff), stat(ff), neigh(ff), domain(ff), cnfig(ff), comm)

    End Do

#ifdef CHRONO
    Call start_timer(tmr, 'Bonded Forces')
#endif

    ! Computation of intramolecular interactions
    Do ff = 1, flow%NUM_FF

      ! Calculate shell model forces

      If (cshell(ff)%megshl > 0) Call core_shell_forces(cshell(ff), stat(ff), cnfig(ff), comm)

      ! Calculate tethered atom forces

      If (tether(ff)%total > 0) Call tethers_forces(stat(ff), tether(ff), cnfig(ff), comm)

      ! Calculate bond forces

      If (bond(ff)%total > 0) Then
        ltmp = (bond(ff)%bin_pdf > 0 .and. &
          ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
          Mod(flow%step,flow%freq_bond) == 0)

         switch = 1 + Merge(1, 0, ltmp)
         Call bonds_forces(switch, stat(ff)%engbnd, stat(ff)%virbnd, stat(ff)%stress, neigh(ff)%cutoff, &
                           stat(ff)%engcpe, stat(ff)%vircpe, bond(ff), mpoles(ff), electro(ff), cnfig(ff), comm)
      End If

      ! Calculate valence angle forces

      If(angle(ff)%total > 0) Then
        ltmp = (angle(ff)%bin_adf > 0 .and. &
        ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
        Mod(flow%step,flow%freq_angle) == 0)

        switch = 1 + Merge(1, 0, ltmp)
        Call angles_forces(switch, stat(ff)%engang, stat(ff)%virang, stat(ff)%stress, angle(ff), cnfig(ff), comm)
      End If

      ! Calculate dihedral forces

      If (dihedral(ff)%total > 0) Then
        ltmp = (dihedral(ff)%bin_adf > 0 .and. &
        ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) &
        .and. Mod(flow%step,flow%freq_dihedral) == 0)

        switch = 1 + Merge(1, 0, ltmp)
        Call dihedrals_forces(switch, stat(ff)%engdih, stat(ff)%virdih, stat(ff)%stress, &
                              neigh(ff)%cutoff, stat(ff)%engcpe, stat(ff)%vircpe, stat(ff)%engsrp, &
                              stat(ff)%virsrp, dihedral(ff), vdws(ff), mpoles(ff), electro(ff), cnfig(ff), comm)
      End If

      ! Calculate inversion forces

      If (inversion(ff)%total > 0) Then
        ltmp = (inversion(ff)%bin_adf > 0 .and. &
        ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
         Mod(flow%step,flow%freq_inversion) == 0)

         switch = 1 + Merge(1, 0, ltmp)
         Call inversions_forces(switch, stat(ff)%enginv, stat(ff)%virinv, stat(ff)%stress, inversion(ff), cnfig(ff), comm)
       End If

    End Do

#ifdef CHRONO
    Call stop_timer(tmr,'Bonded Forces')
#endif

    ! Computation of electrostatics, vdW and metal interactions
    Do ff=1,flow%NUM_FF
      ! Calculate pair-like forces (metal,vdws,electrostatic) and add lrc
      If (.not.(met(ff)%max_metal == 0 .and. electro(ff)%key == ELECTROSTATIC_NULL .and. &
        vdws(ff)%no_vdw .and. rdf(ff)%max_rdf == 0) .or. kim_data(ff)%active) Then
        Call two_body_forces(thermo(ff)%ensemble, flow%book, cnfig(ff)%megfrz, &
                             flow%equilibration, flow%equil_steps, flow%step,stat(ff), ewld(ff), met(ff), &
                             pois(ff), neigh(ff), sites(ff), vdws(ff), rdf(ff), mpoles(ff), electro(ff), domain(ff), &
                             tmr, kim_data(ff), cnfig(ff), comm)
      End If
    End Do

    ! Gather info for stress
    Do ff=1,flow%NUM_FF
      Call gsum(comm, stat(ff)%stress)
    End do

    ! Compute configurational energy
    Do ff=1,flow%NUM_FF
      stat(ff)%stpcfg = stat(ff)%engcpe + stat(ff)%engsrp + stat(ff)%engter + stat(ff)%engtbp + stat(ff)%engfbp + &
                        stat(ff)%engshl + stat(ff)%engtet + stat(ff)%engbnd + stat(ff)%engang + stat(ff)%engdih + &
                        stat(ff)%enginv
    End Do

    ! Compute EVB energy, forces and stress tensor only if number of FFs > 1
    If(flow%NUM_FF > 1)Then
      Call evb_pes(evbff,flow,cnfig,stat)
    End If

    ! Apply external field
    ! For standard EVB we prevent such a calculation in evb_check_external
    Do ff=1,flow%NUM_FF
      If (ext_field(1)%key /= FIELD_NULL) Then
        Call external_field_apply(flow%time, flow%equilibration, flow%equil_steps, flow%step, cshell(ff), stat(ff), rdf(ff), &
                                  ext_field(ff), rigid(ff), domain(ff), cnfig(ff), comm)
        stat(ff)%stpcfg = stat(ff)%stpcfg + stat(ff)%engfld
      End If
    End Do

    Do ff=1,flow%NUM_FF
    ! Apply PLUMED driven dynamics
      If (plume(ff)%l_plumed) Then
         Call plumed_apply(cnfig(ff), flow%run_steps, flow%step, stat(ff), plume(ff), comm)
      End If
    !  Apply pseudo thermostat - force cycle (0)
      If (thermo(ff)%l_stochastic_boundaries) Then
        Call stochastic_boundary_vv(0, thermo(ff)%tstep, flow%step, sites(ff)%dof_site, cshell(ff), stat(ff), thermo(ff), &
                                    rigid(ff), domain(ff), cnfig(ff), seed, comm)
      End If
     ! Cap forces in equilibration mode
      If (flow%step <= flow%equil_steps .and. flow%force_cap) Call cap_forces(thermo(ff)%temp, cnfig(ff), comm)
     ! Frozen atoms option
     Call freeze_atoms(cnfig(ff))
    End Do

    ! Minimisation option and Relaxed shell model optimisation
    Do ff = 1, flow%NUM_FF
      If (flow%simulation .and. (minim(ff)%minimise .or. cshell(ff)%keyshl == SHELL_RELAXED)) Then
        If (cshell(ff)%keyshl == SHELL_RELAXED) Then
          Call core_shell_relax(flow%strict,rdf(ff)%l_collect, stat(ff)%stpcfg, cshell(ff), &
                                stat(ff), domain(ff), cnfig(ff), files, comm)
        End If

        If (.not.cshell(ff)%relaxed) Go To 200 ! Shells relaxation takes priority over minimisation

        If (minim(ff)%minimise .and. flow%step >= 0 .and. flow%step <= flow%run_steps .and. flow%step <= flow%equil_steps) Then
          If(minim(ff)%freq == 0 .and. flow%step == 0) Then
            Call minimise_relax(flow%strict .or. cshell(ff)%keyshl == SHELL_RELAXED, &
                                rdf(ff)%l_collect, thermo(ff)%tstep, stat(ff)%stpcfg, io, stat(ff), pmf(ff), cons(ff), &
                                netcdf, minim(ff), rigid(ff), domain(ff), cnfig(ff), files, comm)
          Else If (minim(ff)%freq >  0 .and. flow%step >  0) Then
            If (Mod(flow%step-flow%equil_steps,minim(ff)%freq) == 0) Then
              Call minimise_relax(flow%strict .or. cshell(ff)%keyshl == SHELL_RELAXED, &
                                  rdf(ff)%l_collect, thermo(ff)%tstep, stat(ff)%stpcfg, io, stat(ff), pmf(ff), cons(ff), &
                                  netcdf, minim(ff), rigid(ff), domain(ff), cnfig(ff), files, comm)
            End If
          End If
        End If

      200     Continue

      ! Refresh mappings
         If (.not.(cshell(ff)%relaxed .and. minim(ff)%relaxed)) Then
            Call refresh_mappings(cnfig(ff), flow, cshell(ff), cons(ff), pmf(ff), stat(ff), msd_data(ff), bond(ff), angle(ff), &
                                  dihedral(ff), inversion(ff), tether(ff), neigh(ff), sites(ff), mpoles(ff), rigid(ff), &
                                  domain(ff), kim_data(ff), ewld(ff), green(ff), minim(ff), thermo(ff), electro(ff), crd(ff), &
                                  comm, tmr)
            If(ff == flow%NUM_FF) Go To 100
         End If
      End If
    End Do

    ! Get RB COM stress and virial at restart only - also available at w_at_start_vv for cnfig%levcfg==2
    If (flow%newjob) Then
      Do ff=1,flow%NUM_FF
        If (rigid(ff)%total > 0) Then
          If (thermo(ff)%l_langevin) Then
            Call langevin_forces(flow%step, thermo(ff)%temp, thermo(ff)%tstep, thermo(ff)%chi, &
                                 thermo(ff)%fxl, thermo(ff)%fyl, thermo(ff)%fzl, cshell(ff), cnfig(ff), seed)
            If (rigid(ff)%share) Then
               Call update_shared_units(cnfig(ff), rigid(ff)%list_shared, rigid(ff)%map_shared, thermo(ff)%fxl, &
                                        thermo(ff)%fyl, thermo(ff)%fzl, domain(ff), comm)
            End If
            Call rigid_bodies_stress(stat(ff)%strcom, cnfig(ff), rigid(ff),comm, thermo(ff)%fxl, thermo(ff)%fyl, thermo(ff)%fzl)
          Else
            Call rigid_bodies_stress(stat(ff)%strcom, rigid(ff), cnfig(ff), comm)
          End If
          stat(ff)%vircom = - (stat(ff)%strcom(1) + stat(ff)%strcom(5) + stat(ff)%strcom(9))
        End If
      End Do
    End If

    If(flow%NUM_FF == 1)Then
      ! Total virial (excluding constraint, PMF and RB COM virials for npt routines)
      ! Total stress (excluding constraint, PMF, RB COM and kinetic stress for npt routines)
      !
      ! NOTE(1):  virsrp already includes vdws%vlrc and vlrcm(0) and so
      !           does the stress diagonal elements (by minus a third),
      !           engsrp includes vdws%elrc and elrcm(0)
      !
      ! NOTE(2):  virfbp, virinv and virdih are allegedly always zero
      stat(1)%virtot = stat(1)%vircpe + stat(1)%virsrp + stat(1)%virter + stat(1)%virtbp + stat(1)%virfbp + &
                       stat(1)%virshl + stat(1)%virtet + stat(1)%virbnd + stat(1)%virang + stat(1)%virdih + &
                       stat(1)%virinv + stat(1)%virfld
    EndIf

    ! If coupling terms are non-zero functions, it is not possible wihtin the EVB framework
    ! to have a decomposition of the energy and virial into separate contributions
    ! for each type of interaction (e.g. angles, bonds, dihedrals, etc).
    ! For this reason, we set all these components to zero in the subroutine ebv_setzero.
    !
    ! The user might want to set no coupling terms between the force fields by choosing zero functions
    ! (setting 'const' equal to zero in the SETEVB file).
    ! Only in this case the energy/virial decomposition is possible and evb_setzero is not called
    If(flow%NUM_FF > 1)Then
      If(.Not. evbff%no_coupling)Then
        Call evb_setzero(flow,stat)
      End If
    End If


    Do ff=1,flow%NUM_FF
      ! If RBs are present update forces on shared ones
      If (rigid(ff)%share) Then
        Call update_shared_units(cnfig(ff), rigid(ff)%list_shared, &
                                 rigid(ff)%map_shared, SHARED_UNIT_UPDATE_FORCES, domain(ff), comm)
      End If
    End Do

  End Subroutine calculate_forces_evb


  !> @brief Compute forces using DFTB+ v18.2
  !!
  !! For a given set of atomic coordinates, forces are computed using
  !! the density-functional tight-binding code, DFTB+ v19.2
  !!
  !! @param[inout] comm      Object containing MPI communicator
  !! @param[in]    flow      Object containing MD step/time data
  !! @param[inout] config    Object containing configuration data
  !!                         Returns updated forces in config
  !! @param[in]    devel     Object containing developement data
  !
  Subroutine calculate_dftb_forces(comm, flow, config, devel)
    Type(comms_type),         Intent( InOut ) :: comm
    Type(flow_type),          Intent( In    ) :: flow
    Type(configuration_type), Intent( InOut ) :: config
    Type(development_type),   Intent( In    ) :: devel

#ifdef DFTBP
    !> Gathered coordinates and mpi index arrays
    Type(coordinate_buffer_type) :: gathered
    !> Gathered atom names
    Character(len=len_atmnam), Allocatable :: atmnam(:)
    !> DFTB+ geometry object
    Type(dftb_geometry_type) :: geo
    !> DFTB+ forces
    Real(wp), Allocatable :: forces(:,:)
    !> DFTB+ atomic Mulliken charges
    Real(wp), Allocatable :: atomic_charges(:)
    !>DFTB+ requires all geometry data on every process
    Logical, Parameter :: to_master_only = .false.
    !> Unit conversion factor
    Real(wp) :: unit_factor

    !NOTE: In principle, only need to init/finalise
    !when atoms move to different processes
    !Gather atomic coordinates and names (excluding halo atoms)
    Call gathered%initialise(comm, config%megatm*3)
    Call gather_coordinates(comm, config, to_master_only, gathered)
    If(.not. Allocated(atmnam)) Allocate(atmnam(config%megatm))
    Call gather_atomic_names(comm, config, atmnam)

    !Number of atoms and species types conserved, hence assign once
    If(flow%step == flow%initial_md_step) Then
       Call geo%initialise(config, atmnam)
       Allocate(forces(3,config%megatm))
       Allocate(atomic_charges(config%megatm))
       unit_factor = convert_unit('DFTB+','DLPOLY','force')
    Endif

    Call geo%set_geometry(comm, config, gathered, atmnam)
    If (devel%app_test%dftb_library) Then
       !Call print_DFTB_geometry_data(geo, flow%step)
       Call run_dftbplus(comm, flow, geo, forces, atomic_charges, &
                         run_app_test = devel%app_test%dftb_library)
    Else
       Call run_dftbplus(comm, flow, geo, forces, atomic_charges)
    Endif

    !TODO(Alex) Consider doing this over unit conversion after assignment to config
    !forces(:,:) = forces(:,:) * unit_factor

    !MPI index arrays same for forces as for coordinates
    Call distribute_forces(comm, gathered%mpi, forces, config)
    config%parts(:)%fxx = config%parts(:)%fxx * unit_factor
    config%parts(:)%fyy = config%parts(:)%fyy * unit_factor
    config%parts(:)%fzz = config%parts(:)%fzz * unit_factor

    Call gathered%finalise()

    If(flow%step == flow%run_steps) Then
       Call geo%finalise()
       Deallocate(forces)
       Deallocate(atomic_charges)
    Endif

#endif
  End Subroutine calculate_dftb_forces

  Subroutine refresh_mappings(cnfig, flow, cshell, cons, pmf, stat, msd_data, bond, angle, &
                              dihedral, inversion, tether, neigh, sites, mpoles, rigid, domain, &
                              kim_data, ewld, green, minim, thermo, &
                              electro, crd, comm, tmr)
    Type(configuration_type), Intent(InOut) :: cnfig
    Type(flow_type),          Intent(InOut) :: flow
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(bonds_type),         Intent(InOut) :: bond
    Type(angles_type),        Intent(InOut) :: angle
    Type(dihedrals_type),     Intent(InOut) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(tethers_type),       Intent(InOut) :: tether
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(site_type),          Intent(InOut) :: sites
    Type(mpole_type),         Intent(InOut) :: mpoles
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(kim_type),           Intent(InOut) :: kim_data
    Type(ewald_type),         Intent(InOut) :: ewld
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(minimise_type),      Intent(InOut) :: minim
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(electrostatic_type), Intent(InOut) :: electro
    Type(coord_type),         Intent(InOut) :: crd
    Type(comms_type),         Intent(InOut) :: comm
    Type(timer_type),         Intent(InOut) :: tmr

    !!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
#ifdef CHRONO
    Call start_timer(tmr, 'Refresh Mappings')
#endif

    ! Scale t=0 reference positions

    If (flow%step > 0) Call xscale(cnfig, thermo%tstep, thermo, stat, neigh, rigid, domain, comm)

    ! Check VNL conditioning

    Call vnl_check(flow%strict, cnfig%width, neigh, stat, domain, cnfig, ewld%bspline%num_splines, kim_data, comm)

    If (neigh%update) Then

      ! Relocate atoms to new domains and restore bonding description

      Call relocate_particles(cnfig%dvar, neigh%cutoff_extended, flow%book, &
                              msd_data%l_msd, cnfig%megatm, flow, cshell, cons, pmf, stat,  thermo, green, &
                              bond, angle, dihedral, inversion, tether, neigh, sites, minim, mpoles, &
                              rigid, domain, cnfig, crd, comm)

      ! Exchange atomic data in border regions

      Call set_halo_particles(electro%key, neigh, sites, mpoles, domain, cnfig, ewld, kim_data, comm) ! inducing in here only

      ! Re-tag RBs when called again after the very first flow%time
      ! when it's done in rigid_bodies_setup <- build_book_intra

      If (rigid%on) Then
        Call rigid_bodies_tags(cnfig, rigid, comm)
        Call rigid_bodies_coms(cnfig, rigid%xxx, rigid%yyy, rigid%zzz, rigid)
      End If

    Else

      ! Exchange atomic positions in border regions

      Call refresh_halo_positions(domain, cnfig, kim_data, comm)
    End If

    ! set and halo rotational matrices and their infinitesimal rotations

    If (mpoles%max_mpoles > 0) Then
      Call mpoles_rotmat_set_halo(mpoles, domain, cnfig, comm)
    End If

#ifdef CHRONO
    Call stop_timer(tmr, 'Refresh Mappings')
#endif
    !!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine refresh_mappings

  Subroutine integrate_vv(stage, flow, cnfig, ttm, cshell, cons, pmf, stat, thermo, &
                          sites, vdws, rigid, domain, seed, tmr, neigh, comm)
    Integer,                  Intent(In   ) :: stage
    Type(flow_type),          Intent(InOut) :: flow
    Type(configuration_type), Intent(InOut) :: cnfig
    Type(ttm_type),           Intent(InOut) :: ttm
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(site_type),          Intent(InOut) :: sites
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(seed_type),          Intent(InOut) :: seed
    Type(timer_type),         Intent(InOut) :: tmr
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(comms_type),         Intent(InOut) :: comm

! used for vv stage control

    !!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Sharlow's splittings for VV only (LFV->VV) DPD thermostat - no variable flow%time-stepping!!!
    ! One-off application for first order splitting and symmetric application for second order splitting
    ! Velocity field change + generation of DPD virial & stat%stress due to random and drag forces

    If (thermo%key_dpd /= DPD_NULL .and. stage == VV_FIRST_STAGE) Then
      Call dpd_thermostat(stage, flow%strict, neigh%cutoff, flow%step, thermo%tstep, stat, thermo, &
                          neigh, rigid, domain, cnfig, seed, comm)
    End If

    ! Integrate equations of motion - velocity verlet

    If (.not. rigid%on) Then
      If (thermo%ensemble == ENS_NVE) Then

        ! Microcanonical ensemble

        Call nve_0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%strkin, stat%engke, thermo, &
           cshell, cons, pmf, stat, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_EVANS) Then

        ! Evans thermostat (Gaussian temperature constraints)

        Call nvt_e0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           thermo%chi_t, &
           stat%strkin, stat%engke, thermo, &
           cshell, cons, pmf, stat, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_LANGEVIN) Then

        ! Langevin thermostat (Stochastic Dynamics)

        Call nvt_l0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NVT_ANDERSON) Then

        ! Andersen thermostat (Stochastic Dynamics)

        Call nvt_a0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NVT_BERENDSEN) Then

        ! Berendsen thermostat

        Call nvt_b0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_NOSE_HOOVER) Then

        ! Nose-Hoover thermostat

        Call nvt_h0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_GENTLE) Then

        ! Gentle-Stochastic thermostat

        Call nvt_g0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, cnfig%degfre, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NVT_LANGEVIN_INHOMO) Then

        ! Inhomogeneous (two-temperature)
        ! Langevin thermostat (Stochastic Dynamics)

        Call nvt_l2_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           stat%strkin, stat%engke, &
           ttm, cshell, cons, pmf, stat, thermo, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NPT_LANGEVIN) Then

        ! Langevin thermostat and isotropic barostat

        Call npt_l0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           cnfig%degfre, stat%virtot, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, domain, &
           tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

        ! Berendsen thermostat and isotropic barostat

        Call npt_b0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%virtot, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, domain, &
           tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

        ! Nose-Hoover thermostat and isotropic barostat

        Call npt_h0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, stat%virtot, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_MTK) Then

        ! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

        Call npt_m0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, stat%virtot, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, domain, &
           tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) Then

        ! Langevin thermostat and barostat anisotropic (cell shape varying)

        Call nst_l0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           cnfig%degfre, stat%stress, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, domain, &
           tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

        ! Berendsen thermostat and barostat anisotropic (cell shape varying)

        Call nst_b0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%stress, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, domain, &
           tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

        ! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

        Call nst_h0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, stat%stress, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_MTK_ANISO) Then

        ! MTK thermostat and anisotropic barostat (cell shape varying)

        Call nst_m0_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, stat%stress, &
           stat%consv, &
           stat%strkin, stat%engke, &
           cshell, cons, pmf, stat, thermo, sites, vdws, domain, &
           tmr, cnfig, comm)

      Else

        ! Invalid ensemble option

        Call error(430)

      End If
    Else
      If (thermo%ensemble == ENS_NVE) Then

        ! Microcanonical ensemble

        Call nve_1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, thermo, cshell, cons, pmf, &
           stat, rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_EVANS) Then

        ! Evans thermostat (Gaussian temperature constraints)

        Call nvt_e1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           thermo%chi_t, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, thermo, cshell, cons, pmf, &
           stat, rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_LANGEVIN) Then

        ! Langevin thermostat (Stochastic Dynamics)

        Call nvt_l1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, cshell, cons, pmf, &
           stat, thermo, rigid, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NVT_ANDERSON) Then

        ! Andersen thermostat (Stochastic Dynamics)

        Call nvt_a1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, cshell, cons, pmf, &
           stat, thermo, sites, rigid, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NVT_BERENDSEN) Then

        ! Berendsen thermostat

        Call nvt_b1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, cshell, cons, pmf, &
           stat, thermo, rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_NOSE_HOOVER) Then

        ! Nose-Hoover thermostat

        Call nvt_h1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, cshell, cons, pmf, &
           stat, thermo, rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NVT_GENTLE) Then

        ! Gentle-Stochastic thermostat

        Call nvt_g1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, cnfig%degfre, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, cshell, cons, pmf, &
           stat, thermo, rigid, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NPT_LANGEVIN) Then

        ! Langevin thermostat and isotropic barostat

        Call npt_l1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           cnfig%degfre, cnfig%degrot, stat%virtot, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

        ! Berendsen thermostat and isotropic barostat

        Call npt_b1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%virtot, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

        ! Nose-Hoover thermostat and isotropic barostat

        Call npt_h1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, cnfig%degrot, stat%virtot, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_MTK) Then

        ! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

        Call npt_m1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, cnfig%degrot, stat%virtot, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) Then

        ! Langevin thermostat and barostat anisotropic (cell shape varying)

        Call nst_l1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           flow%step, &
           cnfig%degfre, cnfig%degrot, stat%stress, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%consv, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, seed, comm)

      Else If (thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

        ! Berendsen thermostat and barostat anisotropic (cell shape varying)

        Call nst_b1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           stat%stress, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

        ! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

        Call nst_h1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, cnfig%degrot, stat%stress, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, comm)

      Else If (thermo%ensemble == ENS_NPT_MTK_ANISO) Then

        ! MTK thermostat and anisotropic barostat (cell shape varying)

        Call nst_m1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
           cnfig%degfre, cnfig%degrot, stat%stress, &
           stat%consv, &
           stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
           stat%strcom, stat%vircom, &
           cshell, cons, pmf, stat, thermo, sites, vdws, &
           rigid, domain, tmr, cnfig, comm)

      Else

        ! Invalid ensemble option

        Call error(430)

      End If
    End If

    ! Sharlow's second order splittings for VV only (LFV->VV) DPD thermostat - no variable flow%time-stepping!!!
    ! Symmetric application for second order splitting
    ! Velocity field change + generation of DPD virial & stat%stress due to random and drag forces

    If (thermo%key_dpd == DPD_SECOND_ORDER .and. stage == VV_SECOND_STAGE) Then
      Call dpd_thermostat(stage, flow%strict, neigh%cutoff, flow%step, thermo%tstep, stat, thermo, &
                          neigh, rigid, domain, cnfig, seed, comm)
    End If

    !!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine integrate_vv

  Subroutine kinetic_options(flow, cnfig, cshell, cons, pmf, stat, sites, ext_field, domain, seed, rigid, thermo, comm)

    Type(flow_type),           Intent(InOut) :: flow
    Type(configuration_type),  Intent(InOut) :: cnfig
    Type(core_shell_type),     Intent(InOut) :: cshell
    Type(constraints_type),    Intent(InOut) :: cons
    Type(pmf_type),            Intent(InOut) :: pmf
    Type(stats_type),          Intent(InOut) :: stat
    Type(site_type),           Intent(InOut) :: sites
    Type(external_field_type), Intent(In   ) :: ext_field
    Type(domains_type),        Intent(In   ) :: domain
    Type(seed_type),           Intent(InOut) :: seed
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(comms_type),          Intent(InOut) :: comm

    Logical :: safe

    !!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Apply external field

    If (ext_field%key /= FIELD_NULL) Then
      Call external_field_correct(stat%engfld, ext_field, rigid, cnfig, comm)
    End If

    ! Apply pseudo thermostat - velocity cycle (1)
    If (thermo%l_stochastic_boundaries) Then
      Call stochastic_boundary_vv(1, thermo%tstep, flow%step, sites%dof_site, cshell, &
                                  stat, thermo, rigid, domain, cnfig, seed, comm)
    End If

    ! Apply temperature regaussing

    If (thermo%l_tgaus .and. flow%step <= flow%equil_steps .and. &
        Mod(flow%step - flow%equil_steps, thermo%freq_tgaus) == 0) Then
      thermo%chi_t = 0.0_wp
      thermo%chi_p = 0.0_wp
      thermo%eta = 0.0_wp

      Call regauss_temperature(rigid, domain, cnfig, seed, comm)

      ! quench constraints & PMFs

      If (cons%megcon > 0) Call constraints_quench(cons, stat, domain, cnfig, comm)
      If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake, cons%tolerance, stat, pmf, cnfig, comm)

      ! quench core-shell units in adiabatic model

      If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
        stat%stptmp = 2.0_wp * (stat%engke + stat%engrot) / (boltz * Real(cnfig%degfre, wp))
        Do
          Call scale_temperature(stat%engke + stat%engrot, cnfig%degtra, cnfig%degrot, cnfig%degfre, rigid, cnfig, comm)
          Call core_shell_quench(cnfig, safe, stat%stptmp, cshell, domain, comm)
          If (cons%megcon > 0) Call constraints_quench(cons, stat, domain, cnfig, comm)
          If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake, cons%tolerance, stat, pmf, cnfig, comm)
          If (rigid%total > 0) Call rigid_bodies_quench(rigid, domain, cnfig, comm)
          If (safe) Exit
        End Do
      Else
        Call scale_temperature(stat%engke + stat%engrot, cnfig%degtra, cnfig%degrot, cnfig%degfre, rigid, cnfig, comm)
      End If

      ! Correct kinetic stress and energy

      If (rigid%total > 0) Then
        Call kinstresf(stat%strknf, cnfig, comm)
        Call kinstrest(rigid, stat%strknt, comm)

        stat%strkin = stat%strknf + stat%strknt

        stat%engrot = getknr(rigid, comm)
      Else
        Call kinstress(stat%strkin, cnfig, comm)
      End If
      stat%engke = 0.5_wp * (stat%strkin(1) + stat%strkin(5) + stat%strkin(9))
    End If

    ! Apply temperature scaling

    If (thermo%l_tscale .and. flow%step <= flow%equil_steps .and. Mod(flow%step - flow%equil_steps, thermo%freq_tscale) == 0) Then
      thermo%chi_t = 0.0_wp
      thermo%chi_p = 0.0_wp
      thermo%eta = 0.0_wp

      ! quench constraints & PMFs

      If (cons%megcon > 0) Call constraints_quench(cons, stat, domain, cnfig, comm)
      If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake, cons%tolerance, stat, pmf, cnfig, comm)

      ! quench core-shell units in adiabatic model

      If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
        Do
          Call scale_temperature(thermo%sigma, cnfig%degtra, cnfig%degrot, cnfig%degfre, rigid, cnfig, comm)
          Call core_shell_quench(cnfig, safe, stat%stptmp, cshell, domain, comm)
          If (cons%megcon > 0) Call constraints_quench(cons, stat, domain, cnfig, comm)
          If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake, cons%tolerance, stat, pmf, cnfig, comm)
          If (rigid%total > 0) Call rigid_bodies_quench(rigid, domain, cnfig, comm)
          If (safe) Exit
        End Do
      Else
        Call scale_temperature(thermo%sigma, cnfig%degtra, cnfig%degrot, cnfig%degfre, rigid, cnfig, comm)
      End If

      ! Correct kinetic stress and energy

      If (rigid%total > 0) Then
        Call kinstresf(stat%strknf, cnfig, comm)
        Call kinstrest(rigid, stat%strknt, comm)

        stat%strkin = stat%strknf + stat%strknt

        stat%engrot = getknr(rigid, comm)
      Else
        Call kinstress(stat%strkin, cnfig, comm)
      End If
      stat%engke = 0.5_wp * (stat%strkin(1) + stat%strkin(5) + stat%strkin(9))
    End If

    !!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine kinetic_options

  Subroutine statistics_report(cnfig, ttm, cshell, cons, pmf, stat, msd_data, zdensity, &
                               sites, domain, flow, files, thermo, tmr, green, minim, comm, ff)

    Type(configuration_type), Intent(InOut) :: cnfig
    Type(ttm_type),           Intent(InOut) :: ttm
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(site_type),          Intent(InOut) :: sites
    Type(domains_type),       Intent(In   ) :: domain
    Type(flow_type),          Intent(InOut) :: flow
    Type(file_type),          Intent(InOut) :: files(:)
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(timer_type),         Intent(InOut) :: tmr
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(minimise_type),      Intent(InOut) :: minim
    Type(comms_type),         Intent(InOut) :: comm
    Integer( Kind = wi ), Intent( In     ), Optional :: ff

    Character(Len=256) :: message, messages(5)

!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    Logical                 :: ffpass

    If (present(ff)) then
      If(ff==1)Then
        ffpass=.True.
      Else
        ffpass=.False.
      End If
    Else
      ffpass= .True.
    Endif


    ! Get complete stress tensor

    stat%strtot = stat%strcon + stat%strpmf + stat%stress + stat%strkin + stat%strcom + stat%strdpd

    ! Get core-shell kinetic energy for adiabatic shell model

    If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
      Call core_shell_kinetic(cnfig, stat%shlke, cshell, domain, comm)
    End If

    ! Calculate physical quantities and collect statistics

    Call statistics_collect &
      (cnfig, flow%simulation, flow%equilibration, flow%equil_steps, msd_data%l_msd, &
       flow%restart_key, &
       cnfig%degfre, cnfig%degshl, cnfig%degrot, &
       flow%step, thermo%tstep, flow%time, flow%start_time, &
       cnfig%mxatdm, stat, thermo, zdensity, sites, files, comm, ff)

    ! VV forces evaluation report for 0th or weird restart

    If (cnfig%levcfg == 1) Then
      Call info('forces evaluated at (re)start for VV integration...', .true.)
      Call info('', .true.)
    End If

    ! line-printer output every flow%freq_output steps

    If ((flow%lines == 0 .or. Mod(flow%step,flow%freq_output) == 0) .and. ffpass)Then

      ! Update cpu flow%time

      Call gtime(tmr%elapsed)

      If (flow%new_page()) Then
        Write (messages(1), '(a)') Repeat('-', 130)
        Write (messages(2), '(9x,a4,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
          'step', 'eng_tot', 'temp_tot', 'eng_cfg', 'eng_src', 'eng_cou', 'eng_bnd', 'eng_ang', 'eng_dih', 'eng_tet'
        Write (messages(3), '(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
          'time(ps)', ' eng_pv', 'temp_rot', 'vir_cfg', 'vir_src', 'vir_cou', 'vir_bnd', 'vir_ang', 'vir_con', 'vir_tet'
        Write (messages(4), '(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,7x,a5,8x,a4,7x,a5,5x,a7,7x,a5)') &
          'cpu  (s)', 'volume', 'temp_shl', 'eng_shl', 'vir_shl', 'alpha', 'beta', 'gamma', 'vir_pmf', 'press'
        Write (messages(5), '(a)') Repeat('-', 130)
        Call info(messages, 5, .true.)
      Else If (ttm%l_ttm) Then
        Write (messages(1), '(a)') Repeat('-', 130)
        Call info(messages, 1, .true.)
      End If

      Write (messages(1), '(i13,1p,9e12.4)') flow%step, stat%stpval(1:9)
      Write (messages(2), '(f13.5,1p,9e12.4)') flow%time, stat%stpval(10:18)
      Write (messages(3), '(0p,f13.3,1p,9e12.4)') tmr%elapsed, stat%stpval(19:27)
      Write (messages(4), '(a)') ''
      Call info(messages, 4, .true.)

      Write (messages(1), '(6x,a7,1p,9e12.4)') 'rolling', stat%ravval(1:9)
      Write (messages(2), '(5x,a8,1p,9e12.4)') 'averages', stat%ravval(10:18)
      Write (messages(3), '(13x,9e12.4)') stat%ravval(19:27)
      Write (messages(4), '(a)') Repeat('-', 130)
      Call info(messages, 4, .true.)

      If (flow%step /= 0) Then
        Call flow%line_printed()
      End If

    End If

    ! Reports at end of equilibration period

    If (flow%step == flow%equil_steps .and. ffpass) Then

      If (flow%step > 0) Then
        Call info(Repeat('-', 130), .true.)
        Write (message, '(a,i10)') 'switching off equilibration at step ', flow%step
        Call info(message, .true.)

        If (thermo%l_zero) Then
          thermo%l_zero = .false.
          Write (message, '(a,i10)') 'switching off zero Kelvin optimiser at step ', flow%step
          Call info(message, .true.)
        End If

        If (minim%minimise) Then
          minim%minimise = .false.
          Write (message, '(a,i10)') 'switching off CGM minimiser at step ', flow%step
          Call info(message, .true.)
        End If

        If (thermo%l_tscale) Then
          thermo%l_tscale = .false.
          Write (message, '(a,i10)') 'switching off temperature scaling at step ', flow%step
          Call info(message, .true.)
        End If

        If (thermo%l_tgaus) Then
          thermo%l_tgaus = .false.
          Write (message, '(a,i10)') 'switching off temperature regaussing at step ', flow%step
          Call info(message, .true.)
        End If
      End If

      ! bond & PMF constraint quenching iterative cycles statistics report

      If (cons%megcon > 0) Then
        Call gmax(comm, stat%passcnq(3:5))
        If (stat%passcnq(3) > 0.0_wp) Then
          Write (message, '(2(a,f5.2),4(a,i3))') &
            'constraints quench run statistics per call: average cycles ', &
            stat%passcnq(3), '/', stat%passcnq(3), &
            ' minimum cycles ', Nint(stat%passcnq(4)), '/', Nint(stat%passcnq(4)), &
            ' maximum cycles ', Nint(stat%passcnq(5)), '/', Nint(stat%passcnq(5))
          Call info(message, .true.)
        End If
      End If

      If (pmf%megpmf > 0) Then
        Call gmax(comm, stat%passpmq(3:5))
        If (stat%passpmq(3) > 0.0_wp) Then
          Write (message, '(2(a,f5.2),4(a,i3))') &
            'PMFs quench run statistics per call: average cycles ', &
            stat%passpmq(3), '/', stat%passpmq(3), &
            ' minimum cycles ', Nint(stat%passpmq(4)), '/', Nint(stat%passpmq(4)), &
            ' maximum cycles ', Nint(stat%passpmq(5)), '/', Nint(stat%passpmq(5))
          Call info(message, .true.)
        End If
      End If

      If (flow%step > 0 .or. cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call info(Repeat('-', 130), .true.)
      End If
    End If

    ! Calculate green-kubo properties

    If (green%samp > 0) Then
      Call vaf_collect(cnfig, sites%mxatyp, flow%equilibration, flow%equil_steps, &
                       flow%step, flow%time, green, comm)
    End If

    If (stat%cur%on .and. (Mod(flow%step, flow%freq_output) == 0) .and. (flow%step > flow%equil_steps)) Then

      Call stat%cur%compute(cnfig, flow%time, comm)

    End If

  End Subroutine statistics_report

  Subroutine write_options(cnfig, io, rsdc, cshell, stat, sites, netcdf, domain, traj, files, dfcts, &
                           flow, thermo, msd_data, green, neigh, comm)

    Type(configuration_type), Intent(InOut) :: cnfig
    Type(io_type),            Intent(InOut) :: io
    Type(rsd_type),           Intent(Inout) :: rsdc
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(stats_type),         Intent(InOut) :: stat
    Type(site_type),          Intent(InOut) :: sites
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(domains_type),       Intent(In   ) :: domain
    Type(trajectory_type),    Intent(InOut) :: traj
    Type(file_type),          Intent(InOut) :: files(:)
    Type(defects_type),       Intent(InOut) :: dfcts(:)
    Type(flow_type),          Intent(InOut) :: flow
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(comms_type),         Intent(InOut) :: comm

    !!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

    If (traj%ltraj) Then
      Call trajectory_write(flow%restart_key, flow%step, thermo%tstep, flow%time, io, stat%rsd, netcdf, cnfig, &
                            traj, files, comm)
    End If

    If (dfcts(1)%ldef) Then
      Call defects_write(flow%restart_key, thermo%ensemble, flow%step, thermo%tstep, &
                         flow%time, io, cshell, dfcts(1), &
                         neigh, sites, netcdf, domain, cnfig, files, comm)
      If (dfcts(2)%ldef) Then
        Call defects_write(flow%restart_key, thermo%ensemble, flow%step, thermo%tstep, &
                           flow%time, io, cshell, dfcts(2), &
                           neigh, sites, netcdf, domain, cnfig, files, comm)
      End If
    End If

    If (msd_data%l_msd) Then
      Call msd_write(cnfig, flow%restart_key, cnfig%megatm, flow%step, thermo%tstep, &
                     flow%time, stat%stpval, sites%dof_site, &
                     io, msd_data, files, comm)
    End If

    If (rsdc%lrsd) Then
      Call rsd_write(flow%restart_key, flow%step, thermo%tstep, io, rsdc, flow%time, &
                     cshell, stat%rsd, cnfig, comm)
    End If

    If (green%samp > 0) Then
      Call vaf_write(cnfig, flow%restart_key, flow%step, thermo%tstep, green, sites, comm)
    End If

    !!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine write_options

  Subroutine refresh_output(files, flow, tmr, comm)

    Type(file_type),  Intent(InOut) :: files(:)
    Type(flow_type),  Intent(InOut) :: flow
    Type(timer_type), Intent(InOut) :: tmr
    Type(comms_type), Intent(InOut) :: comm

    Character(Len=10) :: c_out
    Integer           :: i
    Logical           :: l_out

    !!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

    i = 20
    If (flow%step > 0) Then
      If (Mod(flow%step, i * flow%freq_output) == 0 .or. &
          (tmr%elapsed > Real(i * 60, wp) .and. &
           tmr%elapsed - Real(((Int(tmr%elapsed) / (i * 60)) * i * 60), wp) < &
           tmr%elapsed / Real(flow%step, wp))) Then

        If (comm%idnode == 0) Then
          Inquire (File=files(FILE_OUTPUT)%filename, Exist=l_out, Position=c_out)
          Call strip_blanks(c_out)
          Call lower_case(c_out)
          If (l_out .and. c_out(1:6) == 'append') Then
            Call files(FILE_OUTPUT)%close ()
            Open (Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Position='append')
          End If
        End If

      End If
    End If

    !!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine refresh_output

  Subroutine md_vv(cnfig,ttm,io,rsdc,flow,cshell,cons,pmf,stat,thermo,plume, &
                       pois,bond,angle,dihedral,inversion,zdensity,neigh,sites,fourbody,rdf, &
                       netcdf,mpoles,ext_field,rigid,domain,seed,traj,kim_data,files,tmr,&
                       minim,impa,green,ewld,electro,dfcts,&
                       msd_data,tersoffs,tether,threebody,vdws,devel,met, crd, adf, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to drive a MD simulation with EVB. This subroutine
    ! is the EVB version of md_vv adapted to considering one and multiple
    ! force-fields. This subroutine should replace vv_md upon merging EVB changes.
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type),  Intent(InOut) :: cnfig(:)
    Type(ttm_type ),           Intent(InOut) :: ttm(:)
    Type(io_type ),            Intent(InOut) :: io
    Type(rsd_type ),           Intent(InOut) :: rsdc
    Type(flow_type ),          Intent(InOut) :: flow
    Type(constraints_type ),   Intent(InOut) :: cons(:)
    Type(core_shell_type ),    Intent(InOut) :: cshell(:)
    Type(pmf_type ),           Intent(InOut) :: pmf(:)
    Type(stats_type),          Intent(InOut) :: stat(:)
    Type(thermostat_type),     Intent(InOut) :: thermo(:)
    Type(plumed_type),         Intent(InOut) :: plume(:)
    Type(poisson_type),        Intent(InOut) :: pois(:)
    Type(bonds_type ),         Intent(InOut) :: bond(:)
    Type(angles_type ),        Intent(InOut) :: angle(:)
    Type(dihedrals_type ),     Intent(InOut) :: dihedral(:)
    Type(inversions_type ),    Intent(InOut) :: inversion(:)
    Type(z_density_type ),     Intent(InOut) :: zdensity
    Type(neighbours_type ),    Intent(InOut) :: neigh(:)
    Type(site_type ),          Intent(InOut) :: sites(:)
    Type(four_body_type ),     Intent(InOut) :: fourbody(:)
    Type(rdf_type ),           Intent(InOut) :: rdf(:)
    Type(netcdf_param ),       Intent(In   ) :: netcdf
    Type(mpole_type ),         Intent(InOut) :: mpoles(:)
    Type(external_field_type), Intent(InOut) :: ext_field(:)
    Type(rigid_bodies_type ),  Intent(InOut) :: rigid(:)
    Type(domains_type ),       Intent(In   ) :: domain(:)
    Type(seed_type ),          Intent(InOut) :: seed
    Type(trajectory_type ),    Intent(InOut) :: traj
    Type(kim_type ),           Intent(InOut) :: kim_data(:)
    Type(timer_type ),         Intent(InOut) :: tmr
    Type(file_type ),          Intent(InOut) :: files(:)
    Type(defects_type),        Intent(InOut) :: dfcts(:)
    Type(electrostatic_type ), Intent(InOut) :: electro(:)
    Type(ewald_type),          Intent(InOut) :: ewld(:)
    Type(greenkubo_type ),     Intent(InOut) :: green(:)
    Type(impact_type ),        Intent(InOut) :: impa
    Type(minimise_type ),      Intent(InOut) :: minim(:)
    Type(msd_type ),           Intent(InOut) :: msd_data(:)
    Type(tersoff_type),        Intent(InOut) :: tersoffs(:)
    Type(tethers_type),        Intent(InOut) :: tether(:)
    Type(threebody_type ),     Intent(InOut) :: threebody(:)
    Type(vdw_type ),           Intent(InOut) :: vdws(:)
    Type(development_type ),   Intent(InOut) :: devel
    Type(metal_type ),         Intent(InOut) :: met(:)
    Type(coord_type),          Intent(InOut) :: crd(:)
    Type(adf_type),            Intent(InOut) :: adf(:)
    Type(comms_type )    ,     Intent(InOut) :: comm

    Type( evb_type )            :: evbff
    Integer( Kind = wi )        :: ff
    Logical                     :: fregauss

    Integer                     :: heat_flux_unit
    Real(kind=wp), Dimension(3) :: heat_flux


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV_EVB INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Calculate physical quantities at restart
    ! Calculate kinetic tensor and energy

    !!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV_EVB INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    If(flow%NUM_FF>1)Then
      ! Prevent EVB simulation if there are nonsense settings
      Call evb_prevent(plume(1)%l_plumed, kim_data(1)%active, ttm(1)%l_ttm, thermo(1)%key_dpd)
      ! Allocate EVB variables
      Call evbff%init(flow%NUM_FF)
      ! Read EVB settings
      Call read_evb_settings(evbff, flow, sites, files, comm)
      Call info('Start EVB checking for the consistency of:',.True.)
      ! Check consistency of intra-molecular interactions between different force fiels
      ! for atoms that are not part of the EVB site
      Call evb_check_intramolecular(evbff, flow, sites, bond, angle, dihedral, inversion)
      ! Check consistency of inter-molecular interactions
      Call evb_check_intermolecular(evbff, flow, sites, tersoffs, met, threebody, fourbody)
      ! Check external fields
      Call evb_check_external(evbff, flow, ext_field)
      ! Check consistency between config files
      Call evb_check_configs(cnfig, flow, comm)
      ! Check consistency in the constraint specification between different FFs
      Call evb_check_constraints(evbff, cnfig, cons, cshell, tether, sites, flow, rigid, comm)
      ! Check consistency of intrinsic properties for sites
      Call evb_check_intrinsic(evbff,sites,cnfig,flow,comm)
      ! Check consistency of intrinsic properties for sites
      Call evb_check_vdw(evbff,flow, sites, vdws)
      Call info(' ',.True.)
      Call info('EVB checking was successful !',.True.)
      Call info(' ',.True.)
   End If

    ! Calculate kinetic tensor and energy at restart
    Do ff = 1, flow%NUM_FF

      If (rigid(ff)%total > 0) Then
        Call kinstresf(stat(ff)%strknf, cnfig(ff), comm)
        Call kinstrest(rigid(ff), stat(ff)%strknt, comm)
        stat(ff)%strkin = stat(ff)%strknf + stat(ff)%strknt
      Else
        Call kinstress(stat(ff)%strkin, cnfig(ff), comm)
      End If

      stat(ff)%engke = 0.5_wp * (stat(ff)%strkin(1) +stat(ff)%strkin(5) + stat(ff)%strkin(9))

      ! If cnfig%levcfg=2 and RBs are present, update forces on shared ones
      ! and get RB COM stress and virial at restart.  If cnfig%levcfg<2
      ! forces are calculated at (re)start

      If (cnfig(ff)%levcfg == 2) Then
        If (rigid(ff)%total > 0) Then
          If (rigid(ff)%share) Then
            Call update_shared_units(cnfig(ff), rigid(ff)%list_shared, rigid(ff)%map_shared, &
                                     SHARED_UNIT_UPDATE_FORCES,domain(ff),comm)
          End If

          If (thermo(ff)%l_langevin) Then
            Call langevin_forces(flow%step, thermo(ff)%temp, thermo(ff)%tstep, thermo(ff)%chi, &
                                 thermo(ff)%fxl, thermo(ff)%fyl,thermo(ff)%fzl, cshell(ff), cnfig(ff), seed)
            If (rigid(ff)%share) Then
              Call update_shared_units(cnfig(ff), rigid(ff)%list_shared, rigid(ff)%map_shared, &
                                       thermo(ff)%fxl, thermo(ff)%fyl, thermo(ff)%fzl, domain(ff), comm)
            End If
            Call rigid_bodies_stress(stat(ff)%strcom, cnfig(ff), rigid(ff), comm, thermo(ff)%fxl, thermo(ff)%fyl, thermo(ff)%fzl)
          Else
            Call rigid_bodies_stress(stat(ff)%strcom, rigid(ff), cnfig(ff), comm)
          End If
          stat(ff)%vircom = -(stat(ff)%strcom(1) + stat(ff)%strcom(5) + stat(ff)%strcom(9))
        End If
      End If

    End Do

    !!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV_EVB INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


    ! START OF MOLECULAR DYNAMICS CALCULATIONS

    Do While ( (flow%step < flow%run_steps .or. (flow%step == flow%run_steps .and. flow%newjob)) .and. &
        (tmr%job-tmr%elapsed) > tmr%clear_screen )

      ! Apply impact
      if (impa%active) then
        Do ff=1,flow%NUM_FF
          Call impact_option(cnfig(ff)%levcfg, flow%step, &
               flow%equil_steps, rigid(ff), cshell(ff), stat(ff), impa, cnfig(ff), comm)
        End Do
      end if

      ! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart
      ! cnfig%levcfg == 2 avoids application twice when forces are calculated at (re)start

      If (flow%newjob) Then
        If (cnfig(1)%levcfg == 2) Then
          flow%newjob = .false.
          If (flow%restart_key /= RESTART_KEY_OLD) Then
           Call write_options(cnfig(1), io, rsdc, cshell(1), stat(1), sites(1), netcdf, domain(1), traj, files, dfcts,&
                              flow,thermo(1), msd_data(1), green(1), neigh(1), comm)
          End If
          If (flow%step == 0 .and. flow%step == flow%run_steps) Go To 1000
        End If
      End If


      ! DO THAT ONLY IF 0<=flow%step<flow%run_steps AND FORCES ARE PRESENT (cnfig%levcfg=2)

      If (flow%step >= 0 .and. flow%step < flow%run_steps .and. cnfig(1)%levcfg == 2) Then

        ! Increase step counter

        flow%step=flow%step+1

        ! zero Kelvin structure optimisation
        Do ff = 1, flow%NUM_FF
          If (thermo(ff)%l_zero .and. flow%step <= flow%equil_steps .and. &
              Mod(flow%step-flow%equil_steps,thermo(ff)%freq_zero) == 0) Then
              Call zero_k_optimise(stat(ff), rigid(ff), cnfig(ff), comm)
          End If
        End Do

        Do ff = 1, flow%NUM_FF
        ! Switch on electron-phonon coupling only after flow%time offset
        if (ttm(ff)%l_ttm) ttm(ff)%l_epcp = (flow%time >= ttm(ff)%ttmoffset)

        ! Integrate equations of motion - velocity verlet first stage
          Call integrate_vv(VV_FIRST_STAGE, flow,cnfig(ff), ttm(ff), cshell(ff), cons(ff), pmf(ff), stat(ff), &
                            thermo(ff), sites(ff), vdws(ff), rigid(ff), domain(ff), seed,tmr, neigh(ff), comm)

        ! Refresh mappings
          Call refresh_mappings(cnfig(ff), flow, cshell(ff), cons(ff), pmf(ff), stat(ff), msd_data(ff), bond(ff), &
                                angle(ff), dihedral(ff), inversion(ff), tether(ff), neigh(ff), sites(ff), mpoles(ff),&
                                rigid(ff), domain(ff), kim_data(ff), ewld(ff), green(ff), minim(ff), thermo(ff), &
                                electro(ff), crd(ff), comm, tmr)
        EndDo

      End If ! DO THAT ONLY IF 0<=flow%step<flow%run_steps AND FORCES ARE PRESENT (cnfig%levcfg=2)

      Do ff = 1, flow%NUM_FF
      ! If system is to write per-particle data AND write step AND not equilibration
        If (stat(ff)%require_pp .and. Mod(flow%step, stat(ff)%intsta) == 0 .and. flow%step >= flow%equil_steps) Then
#ifndef HALF_HALO
          Call stat(ff)%allocate_per_particle_arrays(cnfig(ff)%natms)
#else /* HALF_HALO */
          Call stat(ff)%allocate_per_particle_arrays(cnfig(ff)%mxatms)
#endif /* HALF_HALO */
        End If
      End Do

      ! Evaluate forces

      If (flow%simulation_method /= DFTB) Then
        Call calculate_forces_evb(evbff, cnfig, flow, io, cshell, cons, pmf, stat, plume, pois, bond, angle, dihedral, &
                                  inversion, tether, threebody,neigh, sites, vdws, tersoffs, fourbody, rdf, netcdf, &
                                  minim, mpoles, ext_field, rigid, electro, domain, kim_data, msd_data, tmr, &
                                  files, green, devel, ewld, met, seed, thermo, crd, comm)
      Else If (flow%simulation_method == DFTB) Then
         Call calculate_dftb_forces(comm, flow, cnfig(1), devel)
         !Output forces for app test
         !TODO(Alex) Remove this in favour of STATIS
!!$#ifdef DFTBP
!!$         If(devel%app_test%dftb_library) Then
!!$            Call output_dftb_forces(comm, flow, cnfig)
!!$         Endif
!!$#endif
      Endif

     ! If system has written per-particle data
      Do ff = 1, flow%NUM_FF
         If (stat(ff)%collect_pp) Then
          heat_flux = calculate_heat_flux(stat(ff), cnfig(ff), comm)

          If (flow%heat_flux .and. comm%idnode == 0) Then
            If(ff==1)Then
              Open (Newunit=heat_flux_unit, File='HEATFLUX', Position='append')
              Write (heat_flux_unit, '(I8.1, 1X, 5(G19.12, 1X))') flow%step, stat(ff)%stptmp, cnfig(ff)%volm, heat_flux
              Close (heat_flux_unit)
            End If
          End If

          If (flow%write_per_particle) Then
            Call write_per_part_contribs(cnfig(ff), comm, stat(ff)%pp_energy, stat(ff)%pp_stress, flow%step)
          End If

          Call stat(ff)%deallocate_per_particle_arrays()
        End If
      End Do

      ! Calculate physical quantities, collect statistics and report at t=0
      If (flow%step == 0) Then
        Do ff = 1, flow%NUM_FF
          If(ff == 1)Then
            Call crd(ff)%init_coordlist(neigh(ff)%max_list, cnfig(ff)%mxatms)
            Call init_coord_list(cnfig(ff), neigh(ff), crd(ff), sites(ff), flow, comm)
            Call checkcoord(cnfig(ff), crd(ff), sites(ff), flow, stat(ff), comm)
            Call adf_calculate(cnfig(ff), sites(ff), flow, crd(ff), adf(ff), comm)
          End If
          Call statistics_report(cnfig(ff), ttm(ff),cshell(ff), cons(ff), pmf(ff), stat(ff), msd_data(ff), zdensity, &
                                 sites(ff), domain(ff), flow,files, thermo(ff), tmr, green(ff), minim(ff), comm, ff)
        EndDo
      End If


      ! DO THAT ONLY IF 0<flow%step<=flow%run_steps AND THIS IS AN OLD JOB (flow%newjob=.false.)

      If (flow%step > 0 .and. flow%step <= flow%run_steps .and. (.not.flow%newjob)) Then

        ! Evolve electronic temperature for two-temperature model
        Do ff = 1, flow%NUM_FF
          If (ttm(ff)%l_ttm) Then
            Call ttm_ion_temperature(ttm(ff), thermo(ff), domain(ff), cnfig(ff), comm)
            Call ttm_thermal_diffusion(thermo(ff)%tstep, flow%time, flow%step, flow%equil_steps, flow%freq_output, &
                                       flow%freq_restart, flow%run_steps,ttm(ff),thermo(ff),domain(ff),comm)
          End If
        End Do


        ! Integrate equations of motion - velocity verlet second stage
        Do ff = 1, flow%NUM_FF
          Call integrate_vv(VV_SECOND_STAGE, flow, cnfig(ff), ttm(ff), cshell(ff), cons(ff), pmf(ff), stat(ff), &
                            thermo(ff), sites(ff), vdws(ff), rigid(ff), domain(ff), seed, tmr, neigh(ff),comm)

        End Do

        ! Apply kinetic options
        fregauss=.False.

        If(thermo(1)%l_tgaus .and. (flow%step <= flow%equil_steps) .and. Mod(flow%step-flow%equil_steps,thermo(1)%freq_tgaus) == 0)&
          fregauss=.True.

        If(thermo(1)%l_stochastic_boundaries .or. fregauss)Then
          Call kinetic_options(flow, cnfig(1), cshell(1), cons(1), pmf(1), stat(1), sites(1), ext_field(1), domain(1), &
                               seed, rigid(1), thermo(1), comm)

          If(flow%NUM_FF>1)Then
             Call evb_merge_stochastic(flow, cnfig, stat, rigid, thermo, cshell, cons, pmf)
          End If
        Else
          Do ff=1,flow%NUM_FF
           Call kinetic_options(flow, cnfig(ff), cshell(ff), cons(ff), pmf(ff), stat(ff), sites(ff), ext_field(ff), domain(ff), &
                                seed, rigid(ff), thermo(ff), comm)
          End Do
        End If

        ! Evaluate and write EVB population
         If(flow%NUM_FF > 1 .and. evbff%population)then
          Call evb_population(evbff, flow, files, comm)
         End If

        ! Update total flow%time of simulation
        flow%time = flow%time + thermo(1)%tstep

        ! Calculate physical quantities, collect statistics and report regularly
        Do ff = 1, flow%NUM_FF
          Call statistics_report(cnfig(ff), ttm(ff), cshell(ff), cons(ff), pmf(ff), stat(ff), msd_data(ff), zdensity, &
                                 sites(ff), domain(ff), flow, files, thermo(ff), tmr, green(ff), minim(ff), comm, ff)
        End Do

        ! Write HISTORY, DEFECTS, MSDTMP & DISPDAT
        Call write_options(cnfig(1), io, rsdc, cshell(1), stat(1), sites(1), netcdf, domain(1), traj, files, dfcts,&
                           flow,thermo(1), msd_data(1), green(1), neigh(1), comm)

        ! Save restart data in event of system crash
        If (Mod(flow%step,flow%freq_restart) == 0 .and. flow%step /= flow%run_steps .and. (.not.devel%l_tor)) Then
        ! Only print the data from FIELD 1
           Call system_revive(neigh(1)%cutoff,flow%step,flow%time,sites(1),io,flow%start_time, &
            stat(1),devel,green(1),thermo(1),bond(1),angle(1),dihedral(1),inversion(1),zdensity,rdf(1), &
            netcdf,cnfig(1),files,comm)
        End If

        Do ff = 1, flow%NUM_FF
          If(ff == 1)Then
            Call init_coord_list(cnfig(ff), neigh(ff), crd(ff), sites(ff), flow, comm)
            Call checkcoord(cnfig(ff), crd(ff), sites(ff), flow, stat(ff), comm)
            Call adf_calculate(cnfig(ff), sites(ff), flow, crd(ff), adf(ff), comm)
          End If
        End Do


      End If ! DO THAT ONLY IF 0<flow%step<=flow%run_steps AND THIS IS AN OLD JOB (flow%newjob=.false.)

      1000 Continue ! Escape forces evaluation at t=0 when flow%step=flow%run_steps=0 and flow%newjob=.false.

      ! Refresh output

      Call refresh_output(files, flow, tmr, comm)

      ! Complete flow%time check

      Call gtime(tmr%elapsed)

      ! Change cnfig%levcfg appropriately
      Do ff = 1, flow%NUM_FF
        If (cnfig(ff)%levcfg == 1) cnfig(ff)%levcfg=2
      End Do

    End Do

      !Close EVB population file if opened
      If (evbff%population_file_open .and. comm%idnode==0) Then
        Call files(FILE_POPEVB)%close()
      End If

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV_EVB INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
#ifdef EXPERIMENT
      call mem_picture(files(FILE_OUTPUT)%unit_no, comm)
#endif


  End Subroutine md_vv

  Subroutine replay_history(cnfig, io, rsdc, flow, cshell, cons, pmf, stat, thermo, msd_data, &
                            met, pois, bond, angle, dihedral, inversion, zdensity, neigh, sites, vdws, rdf, &
                            netcdf, minim, mpoles, ext_field, rigid, electro, domain, seed, traj, kim_data, dfcts, files, &
                            tmr, tether, green, ewld, devel, comm)

    Type(configuration_type),  Intent(InOut) :: cnfig
    Type(io_type),             Intent(InOut) :: io
    Type(rsd_type),            Intent(Inout) :: rsdc
    Type(flow_type),           Intent(InOut) :: flow
    Type(core_shell_type),     Intent(InOut) :: cshell
    Type(constraints_type),    Intent(InOut) :: cons
    Type(pmf_type),            Intent(InOut) :: pmf
    Type(stats_type),          Intent(InOut) :: stat
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(metal_type),          Intent(InOut) :: met
    Type(poisson_type),        Intent(InOut) :: pois
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(site_type),           Intent(InOut) :: sites
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(netcdf_param),        Intent(In   ) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(In   ) :: domain
    Type(seed_type),           Intent(InOut) :: seed
    Type(trajectory_type),     Intent(InOut) :: traj
    Type(kim_type),            Intent(InOut) :: kim_data
    Type(defects_type),        Intent(inOut) :: dfcts(:)
    Type(file_type),           Intent(InOut) :: files(:)
    Type(timer_type),          Intent(InOut) :: tmr
    Type(tethers_type),        Intent(InOut) :: tether
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(development_type),    Intent(InOut) :: devel
    Type(comms_type),          Intent(InOut) :: comm

    Character(Len=10)  :: c_out
    Character(Len=256) :: messages(6)
    Integer            :: exout, i
    Integer(Kind=wi)   :: nstpe, nstph, switch
    Logical            :: l_out
    Real(Kind=wp)      :: tmsh, tsths

! flow%start_time replacement
! flow%step replacements
! exit indicator for reading

    !!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORY INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
    ! Report work

    Call info('*** HISTORY is replayed for recalculation of structural properties ***', .true.)

    ! Stay safe

    If (traj%ltraj) Then
      traj%ltraj = .false.
      Call info('*** warning - aborting printing into HISTORY while reading it ***', .true.)
    End If

    ! Make sure of no equilibration

    flow%equil_steps = 0
    flow%equilibration = .false.

    ! nullify forces
    Do i = 1, cnfig%mxatms
      cnfig%parts(i)%fxx = 0.0_wp
      cnfig%parts(i)%fyy = 0.0_wp
      cnfig%parts(i)%fzz = 0.0_wp
    End Do

    ! nullify all two-body force switches = just do rdf%rdf calculation

    electro%key = ELECTROSTATIC_NULL
    vdws%n_vdw = 0
    met%n_potentials = 0

    ! defect detection for every entry in HISTORY

    dfcts(:)%nsdef = 0
    dfcts(:)%isdef = 1

    ! MSDTMP option for every entry in HISTORY

    msd_data%start = 0
    msd_data%freq = 1

    ! displacement detection for every entry in HISTORY

    rsdc%nsrsd = 0
    rsdc%isrsd = 1

    ! intramolecular PDF analysis for every entry in HISTORF
    ! enforce printing and collection if the calculation exists

    stat%lpana = (cnfig%mxgana > 0)
    flow%freq_bond = 1
    flow%freq_angle = 1
    flow%freq_dihedral = 1
    flow%freq_inversion = 1

    ! rdf%rdf and z-density detection for every entry in HISTORF
    ! enforce printing and collection if the calculation exists

    rdf%l_print = rdf%l_collect; rdf%freq = 1
    zdensity%l_print = zdensity%l_collect; zdensity%frequency = 1

    ! Calculate kinetic tensor and energy at restart as it may not exists later

    If (rigid%total > 0) Then
      Call kinstresf(stat%strknf, cnfig, comm)
      Call kinstrest(rigid, stat%strknt, comm)

      stat%strkin = stat%strknf + stat%strknt
    Else
      Call kinstress(stat%strkin, cnfig, comm)
    End If
    stat%engke = 0.5_wp * (stat%strkin(1) + stat%strkin(5) + stat%strkin(9))

    nstpe = flow%step
    nstph = 0 ! trajectory points counter
    Do
      Call stat%init_connect(cnfig%mxatdm)
      10 Continue
      If (nstph > nstpe) Then
        Call statistics_connect_set(cnfig, neigh%cutoff_extended, cnfig%mxatdm, msd_data%l_msd, stat, domain, comm)
      End If

      ! Make a move - Read a frame

      Call read_history(flow%strict, files(FILE_HISTORY)%filename, cnfig%megatm, cnfig%levcfg, cnfig%dvar, &
                        flow%step, thermo%tstep, flow%time, exout, io, traj, sites, domain, cnfig, files, comm)

      If (traj%restart) Then
        traj%restart = .false.

        flow%start_time = flow%time
        tmsh = 0.0_wp ! flow%start_time substitute
      End If

      If (exout == 0) Then
        nstph = nstph + 1

        If (nstph <= nstpe) Then

          ! Deal with restarts but remember the old cell parameters

          stat%clin = cnfig%cell
          Go To 10

        Else

          ! CHECK MD CONFIGURATION

          Call check_config(cnfig, electro%key, thermo, sites, flow, comm)

          ! First frame positions (for estimates of MSD when cnfig%levcfg==0)

          If (nstph == 1) Then
            Do i = 1, cnfig%natms
              stat%xin(i) = cnfig%parts(i)%xxx
              stat%yin(i) = cnfig%parts(i)%yyy
              stat%zin(i) = cnfig%parts(i)%zzz
            End Do
            stat%clin = cnfig%cell
            !              xin(natms+1: ) = 0.0_wp
            !              yin(natms+1: ) = 0.0_wp
            !              zin(natms+1: ) = 0.0_wp
            Call statistics_connect_set(cnfig, neigh%cutoff_extended, cnfig%mxatdm, msd_data%l_msd, stat, domain, comm)
          End If

          ! get xto/xin/msdtmp arrays sorted

          Call statistics_connect_frames(cnfig, cnfig%megatm, cnfig%mxatdm, msd_data%l_msd, stat, domain, comm)
          Call stat%clean_connect()

          ! SET domain borders and link-cells as default for new jobs
          ! exchange atomic data and positions in border regions

          Call set_halo_particles(electro%key, neigh, sites, mpoles, domain, cnfig, ewld, kim_data, comm)

          ! For any intra-like interaction, construct book keeping arrays and
          ! exclusion arrays for overlapped two-body inter-like interactions

          If (flow%book) Then
            Call build_book_intra(flow%strict, flow%print_topology, flow%simulation, &
                                  flow, cshell, cons, pmf, bond, angle, dihedral, &
                                  inversion, tether, neigh, sites, rigid, domain, cnfig, comm)
            If (flow%exclusions) Then
              Call build_excl_intra(electro%lecx, cshell, cons, bond, angle, dihedral, &
                                    inversion, neigh, rigid, cnfig, comm)
            End If
          End If

          ! Accumulate RDFs if needed (flow%step->nstph)
          ! Make sure RDFs are complete (flow%book=.false. - no exclusion lists)

          If (rdf%l_collect) Then
            ! Initialise variables for two_body interaction (including long range, which is not strictly a two_body interaction problem)
            stat%engsrp = 0.0_wp
            stat%virsrp = 0.0_wp
            stat%engcpe = 0.0_wp
            stat%vircpe = 0.0_wp
            ! Set up non-bonded interaction (verlet) list using link cells
            If (neigh%update) Then
              Call link_cell_pairs(vdws%cutoff, met%rcut, flow%book, cnfig%megfrz, cshell, devel, &
                                   neigh, mpoles, domain, tmr, cnfig, comm)
            End If
            Call two_body_forces(thermo%ensemble, .false., cnfig%megfrz, &
                                 flow%equilibration, flow%equil_steps, nstph, stat, ewld, met, pois, neigh, sites, &
                                 vdws, rdf, mpoles, electro, domain, tmr, kim_data, cnfig, comm)
          End If

          ! Calculate bond forces

          If (bond%total > 0 .and. bond%bin_pdf > 0) Then
            switch = 0
            Call bonds_forces(switch, stat%engbnd, stat%virbnd, stat%stress, &
                              neigh%cutoff, stat%engcpe, stat%vircpe, bond, mpoles, electro, cnfig, comm)
          End If

          ! Calculate valence angle forces

          If (angle%total > 0 .and. angle%bin_adf > 0) Then
            switch = 0
            Call angles_forces(switch, stat%engang, stat%virang, stat%stress, angle, cnfig, comm)
          End If

          ! Calculate dihedral forces

          If (dihedral%total > 0 .and. dihedral%bin_adf > 0) Then
            switch = 0
            Call dihedrals_forces(switch, stat%engdih, stat%virdih, stat%stress, &
                                  neigh%cutoff, stat%engcpe, stat%vircpe, stat%engsrp, stat%virsrp, &
                                  dihedral, vdws, mpoles, electro, cnfig, comm)
          End If

          ! Calculate inversion forces

          If (inversion%total > 0 .and. inversion%bin_adf > 0) Then
            switch = 0
            Call inversions_forces(switch, stat%enginv, stat%virinv, stat%stress, inversion, cnfig, comm)
          End If

          ! Calculate kinetic stress and energy if available

          If (cnfig%levcfg > 0 .and. cnfig%levcfg < 3) Then
            If (rigid%total > 0) Then
              Call rigid_bodies_quench(rigid, domain, cnfig, comm)

              Call kinstresf(stat%strknf, cnfig, comm)
              Call kinstrest(rigid, stat%strknt, comm)

              stat%strkin = stat%strknf + stat%strknt

              stat%engrot = getknr(rigid, comm)
              If (cnfig%levcfg == 2) Then
                Call rigid_bodies_stress(stat%strcom, rigid, cnfig, comm)
                stat%vircom = -(stat%strcom(1) + stat%strcom(5) + stat%strcom(9))
              End If
            Else
              Call kinstress(stat%strkin, cnfig, comm)
            End If
            stat%engke = 0.5_wp * (stat%strkin(1) + stat%strkin(5) + stat%strkin(9))

            ! Apply kinetic options

            Call kinetic_options(flow, cnfig, cshell, cons, pmf, stat, sites, ext_field, domain, seed, rigid, thermo, comm)

            ! Get core-shell kinetic energy for adiabatic shell model

            If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
              Call core_shell_kinetic(cnfig, stat%shlke, cshell, domain, comm)
            End If
          End If

          ! Get complete stress tensor

          stat%strtot = stat%strcon + stat%strpmf + stat%stress + stat%strkin + stat%strcom + stat%strdpd

          ! Calculate physical quantities and collect statistics,
          ! accumulate z-density if needed
          ! (flow%step->nstph,tstep->tsths,flow%start_time->tmsh)

          tsths = Max(thermo%tstep, (flow%time - tmsh) / Real(Merge(nstph - 1, 1, nstph > 2), wp))

          ! Collect VAF if kinetics is available

          Call vaf_collect(cnfig, sites%mxatyp, flow%equilibration, flow%equil_steps, nstph - 1, flow%time, green, comm)

          Call statistics_collect &
            (cnfig, flow%simulation, flow%equilibration, flow%equil_steps, msd_data%l_msd, &
             flow%restart_key, &
             cnfig%degfre, cnfig%degshl, cnfig%degrot, &
             nstph, tsths, flow%time, tmsh, &
             cnfig%mxatdm, stat, thermo, &
             zdensity, sites, files, comm)

          ! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

          If (traj%ltraj) Call trajectory_write(flow%restart_key, flow%step, thermo%tstep, flow%time, &
                                                io, stat%rsd, netcdf, cnfig, traj, files, comm)
          If (dfcts(1)%ldef) Then
            Call defects_write(flow%restart_key, thermo%ensemble, flow%step, thermo%tstep, flow%time, io, cshell, &
                               dfcts(1), neigh, sites, netcdf, domain, cnfig, files, comm)
            If (dfcts(2)%ldef) Then
              Call defects_write(flow%restart_key, thermo%ensemble, flow%step, thermo%tstep, flow%time, &
                                 io, cshell, dfcts(2), neigh, sites, netcdf, domain, cnfig, files, comm)
            End If
          End If
          If (msd_data%l_msd) Then
            Call msd_write(cnfig, flow%restart_key, cnfig%megatm, flow%step, thermo%tstep, flow%time, stat%stpval, &
                           sites%dof_site, io, msd_data, files, comm)
          End If
          If (rsdc%lrsd) Then
            Call rsd_write(flow%restart_key, flow%step, thermo%tstep, io, rsdc, &
                           flow%time, cshell, stat%rsd, cnfig, comm)
          End If
          If (green%samp > 0) Then
            Call vaf_write(cnfig, flow%restart_key, nstph, tsths, green, sites, comm)
          End If

          ! Complete flow%time check

          Call gtime(tmr%elapsed)
          Write (messages(1), '(2(a,i10),a)') &
            'HISTORY step ', flow%step, ' (', nstph, ' entry) processed'
          Write (messages(2), '(a,f12.3,a)') &
            'time elapsed since job start: ', tmr%elapsed, ' sec'
          Call info(messages, 2, .true.)

          ! Save restart data in event of system crash

          If (Mod(nstph, flow%freq_restart) == 0 .and. nstph /= flow%run_steps .and. (.not. devel%l_tor)) Then
            Call system_revive(neigh%cutoff, flow%step, flow%time, sites, io, flow%start_time, &
                               stat, devel, green, thermo, bond, angle, dihedral, inversion, zdensity, &
                               rdf, netcdf, cnfig, files, comm)
          End If

          ! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

          i = 20
          If (Mod(nstph, i) == 0 .or. &
              (tmr%elapsed > Real(i * 60, wp) .and. &
               tmr%elapsed - Real(((Int(tmr%elapsed) / (i * 60)) * i * 60), wp) < &
               tmr%elapsed / Real(nstph, wp))) Then

            If (comm%idnode == 0) Then
              Inquire (File=files(FILE_OUTPUT)%filename, Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                Call files(FILE_OUTPUT)%close ()
                Open (Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Position='append')
              End If
            End If

          End If
        End If

        ! Save last frame positions (for estimates of MSD when cnfig%levcfg==0)

        Do i = 1, cnfig%natms
          stat%xin(i) = cnfig%parts(i)%xxx
          stat%yin(i) = cnfig%parts(i)%yyy
          stat%zin(i) = cnfig%parts(i)%zzz
        End Do
        stat%clin = cnfig%cell
      Else
        Exit
      End If
    End Do

    ! Finish with grace

    If (exout > 0) Then ! normal exit

      ! recover connectivity arrays for REVCON, REVIVE and printing purposes
      ! read_history MUST NOT initialise R,V,F arrays!!!

      cnfig%ltg(1:cnfig%natms) = stat%ltg0(1:cnfig%natms)
      cnfig%lsa(1:cnfig%natms) = stat%lsa0(1:cnfig%natms)
      cnfig%lsi(1:cnfig%natms) = stat%lsi0(1:cnfig%natms)

    Else If (exout < 0) Then ! abnormal exit

      ! If reading HISTORY finished awkwardly
      ! recover positions and generate kinetics

      Do i = 1, cnfig%natms
        cnfig%parts(i)%xxx = stat%xin(i)
        cnfig%parts(i)%yyy = stat%yin(i)
        cnfig%parts(i)%zzz = stat%zin(i)
      End Do
      cnfig%cell = stat%clin

      Call set_temperature(flow%restart_key, flow%step, flow%run_steps, &
                           stat%engrot, sites%dof_site, cshell, stat, cons, pmf, &
                           thermo, minim, rigid, domain, cnfig, seed, comm)

    End If
    Call stat%clean_connect()

    ! Save restart data because of next action (and disallow the same in dl_poly)

    If (.not. devel%l_tor) Then
      Call system_revive(neigh%cutoff, flow%step, flow%time, sites, io, flow%start_time, stat, &
                         devel, green, thermo, bond, angle, dihedral, inversion, zdensity, rdf, netcdf, &
                         cnfig, files, comm)
    End If

    ! step counter is data counter now, so statistics_result is triggered

    flow%step = nstph
    thermo%tstep = tsths

    !!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORY INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine replay_history

  Subroutine replay_historf(cnfig, io, rsdc, flow, cshell, cons, pmf, stat, thermo, plume, &
                            msd_data, bond, angle, dihedral, inversion, zdensity, neigh, sites, vdws, tersoffs, &
                            fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, seed, traj, &
                            kim_data, files, dfcts, tmr, tether, threebody, pois, green, ewld, devel, met, &
                            crd, comm)

    Type(configuration_type),  Intent(InOut) :: cnfig
    Type(io_type),             Intent(InOut) :: io
    Type(rsd_type),            Intent(Inout) :: rsdc
    Type(flow_type),           Intent(InOut) :: flow
    Type(core_shell_type),     Intent(InOut) :: cshell
    Type(constraints_type),    Intent(InOut) :: cons
    Type(pmf_type),            Intent(InOut) :: pmf
    Type(stats_type),          Intent(InOut) :: stat
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(plumed_type),         Intent(InOut) :: plume
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(site_type),           Intent(InOut) :: sites
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(netcdf_param),        Intent(In   ) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(In   ) :: domain
    Type(seed_type),           Intent(InOut) :: seed
    Type(trajectory_type),     Intent(InOut) :: traj
    Type(kim_type),            Intent(InOut) :: kim_data
    Type(file_type),           Intent(InOut) :: files(:)
    Type(defects_type),        Intent(inOut) :: dfcts(:)
    Type(timer_type),          Intent(InOut) :: tmr
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(poisson_type),        Intent(InOut) :: pois
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(development_type),    Intent(InOut) :: devel
    Type(metal_type),          Intent(InOut) :: met
    Type(coord_type),          Intent(InOut) :: crd
    Type(comms_type),          Intent(InOut) :: comm

    Character(Len=10)  :: c_out
    Character(Len=256) :: messages(5)
    Integer            :: exout, i, nstpe, nstph
    Logical            :: l_out
    Real(Kind=wp)      :: tmsh, tsths

! flow%start_time replacement
! flow%step replacements
! exit indicator for reading
!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

    ! Report work

    Write (messages(1), '(a)') ''
    Write (messages(2), '(a)') '*** HISTORF will be replayed in full (with no dynamics)!!!     ***'
    Write (messages(3), '(a)') '*** Large particle displacements between frames within HISTROF ***'
    Write (messages(4), '(a)') '*** w.r.t. CONFIG at start may lead failures in parallel!!!    ***'
    Call info(messages, 4, .true.)

    ! defect detection for every entry in HISTORF

    dfcts(:)%nsdef = 0
    dfcts(:)%isdef = 1

    ! MSDTMP option for every entry in HISTORF

    msd_data%start = 0
    msd_data%freq = 1

    ! displacement detection for every entry in HISTORF

    rsdc%nsrsd = 0
    rsdc%isrsd = 1

    ! intramolecular PDF analysis for every entry in HISTORF
    ! enforce printing and collection if the calculation exists

    stat%lpana = (cnfig%mxgana > 0)
    flow%freq_bond = 1
    flow%freq_angle = 1
    flow%freq_dihedral = 1
    flow%freq_inversion = 1

    ! rdf%rdf and z-density detection for every entry in HISTORF
    ! enforce printing and collection if the calculation exists

    rdf%l_print = rdf%l_collect; rdf%freq = 1
    zdensity%l_print = zdensity%l_collect; zdensity%frequency = 1

    ! Calculate kinetic tensor and energy at restart as it may not exists later

    If (rigid%total > 0) Then
      Call kinstresf(stat%strknf, cnfig, comm)
      Call kinstrest(rigid, stat%strknt, comm)

      stat%strkin = stat%strknf + stat%strknt
    Else
      Call kinstress(stat%strkin, cnfig, comm)
    End If
    stat%engke = 0.5_wp * (stat%strkin(1) + stat%strkin(5) + stat%strkin(9))

    nstpe = flow%step
    nstph = 0 ! HISTORF trajectory points counter
    Do
      Call stat%init_connect(cnfig%mxatdm)
      10 Continue
      If (nstph > nstpe) Call statistics_connect_set(cnfig, neigh%cutoff_extended, cnfig%mxatdm, msd_data%l_msd, stat, domain, comm)

      ! Make a move - Read a frame

      Call read_history(flow%strict, files(FILE_HISTORF)%filename, cnfig%megatm, cnfig%levcfg, cnfig%dvar, &
                        flow%step, thermo%tstep, flow%time, exout, io, traj, sites, domain, cnfig, files, comm)

      If (traj%restart) Then
        traj%restart = .false.

        flow%start_time = flow%time
        tmsh = 0.0_wp ! flow%start_time substitute
      End If

      If (exout == 0) Then
        nstph = nstph + 1

        If (nstph <= nstpe) Then

          ! Deal with restarts but remember the old cell parameters

          stat%clin = cnfig%cell
          Go To 10

        Else

          ! CHECK MD CONFIGURATION

          Call check_config(cnfig, electro%key, thermo, sites, flow, comm)

          ! First frame positions (for estimates of MSD when cnfig%levcfg==0)

          If (nstph == 1) Then
            Do i = 1, cnfig%natms
              stat%xin(i) = cnfig%parts(i)%xxx
              stat%yin(i) = cnfig%parts(i)%yyy
              stat%zin(i) = cnfig%parts(i)%zzz
            End Do
            stat%clin = cnfig%cell
            !              xin(natms+1: ) = 0.0_wp
            !              yin(natms+1: ) = 0.0_wp
            !              zin(natms+1: ) = 0.0_wp
            Call statistics_connect_set(cnfig, neigh%cutoff_extended, cnfig%mxatdm, msd_data%l_msd, stat, domain, comm)
          End If

          ! get xto/xin/msdtmp arrays sorted

          Call statistics_connect_frames(cnfig, cnfig%megatm, cnfig%mxatdm, msd_data%l_msd, stat, domain, comm)
          Call stat%clean_connect()

          ! SET domain borders and link-cells as default for new jobs
          ! exchange atomic data and positions in border regions

          Call set_halo_particles(electro%key, neigh, sites, mpoles, domain, cnfig, ewld, kim_data, comm)

          ! For any intra-like interaction, construct book keeping arrays and
          ! exclusion arrays for overlapped two-body inter-like interactions

          If (flow%book) Then
            Call build_book_intra(flow%strict, flow%print_topology, flow%simulation, &
                                  flow, cshell, cons, pmf, bond, angle, dihedral, &
                                  inversion, tether, neigh, sites, rigid, domain, cnfig, comm)
            If (flow%exclusions) Then
              Call build_excl_intra(electro%lecx, cshell, cons, bond, angle, dihedral, &
                                    inversion, neigh, rigid, cnfig, comm)
            End If
          End If

          ! Evaluate forces, flow%newjob must always be true for vircom evaluation

          Call calculate_forces_for_historf(cnfig, flow, io, cshell, cons, pmf, stat, plume, pois, bond, angle, dihedral, &
                                          inversion, tether, threebody, neigh, sites, vdws, tersoffs, fourbody, rdf, &
                                          netcdf, minim, mpoles, ext_field, rigid, electro, domain, kim_data, msd_data, &
                                          tmr, files, green, devel, ewld, met, seed, thermo, crd, comm)

          ! Evaluate kinetics if available

          If (cnfig%levcfg > 0 .and. cnfig%levcfg < 3) Then
            If (thermo%l_zero .and. &
                flow%step <= flow%equil_steps .and. &
                Mod(flow%step + 1 - flow%equil_steps, thermo%freq_zero) == 0) Then
              Call zero_k_optimise(stat, rigid, cnfig, comm)
            End If

            If (thermo%l_zero .and. flow%step <= flow%equil_steps) Then
              Call zero_k_optimise(stat, rigid, cnfig, comm)
            End If

            ! Calculate kinetic stress and energy if available

            If (rigid%total > 0) Then
              Call rigid_bodies_quench(rigid, domain, cnfig, comm)

              Call kinstresf(stat%strknf, cnfig, comm)
              Call kinstrest(rigid, stat%strknt, comm)

              stat%strkin = stat%strknf + stat%strknt

              stat%engrot = getknr(rigid, comm)
              Call rigid_bodies_stress(stat%strcom, rigid, cnfig, comm)
              stat%vircom = -(stat%strcom(1) + stat%strcom(5) + stat%strcom(9))
            Else
              Call kinstress(stat%strkin, cnfig, comm)
            End If
            stat%engke = 0.5_wp * (stat%strkin(1) + stat%strkin(5) + stat%strkin(9))

            ! Apply kinetic options

            Call kinetic_options(flow, cnfig, cshell, cons, pmf, stat, sites, ext_field, domain, seed, rigid, thermo, comm)

            ! Get core-shell kinetic energy for adiabatic shell model

            If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
              Call core_shell_kinetic(cnfig, stat%shlke, cshell, domain, comm)
            End If
          End If

          ! Get complete stress tensor

          stat%strtot = stat%strcon + stat%strpmf + stat%stress + stat%strkin + stat%strcom + stat%strdpd

          ! Calculate physical quantities and collect statistics,
          ! accumulate z-density if needed
          ! (flow%step->nstph,tstep->tsths,flow%start_time->tmsh)

          tsths = Max(thermo%tstep, (flow%time - tmsh) / Real(Merge(nstph - 1, 1, nstph > 2), wp))

          ! Collect VAF if kinetics is available

          Call vaf_collect(cnfig, sites%mxatyp, flow%equilibration, flow%equil_steps, nstph - 1, flow%time, green, comm)

          Call statistics_collect &
            (cnfig, flow%simulation, flow%equilibration, flow%equil_steps, msd_data%l_msd, &
             flow%restart_key, &
             cnfig%degfre, cnfig%degshl, cnfig%degrot, &
             nstph, tsths, flow%time, tmsh, &
             cnfig%mxatdm, stat, thermo, &
             zdensity, sites, files, comm)

          ! line-printer output
          ! Update cpu flow%time

          Call gtime(tmr%elapsed)
          If (flow%new_page()) Then
            Write (messages(1), '(a)') Repeat('-', 130)
            Write (messages(2), '(9x,a4,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
              'step', 'eng_tot', 'temp_tot', 'eng_cfg', 'eng_src', 'eng_cou', 'eng_bnd', 'eng_ang', 'eng_dih', 'eng_tet'
            Write (messages(3), '(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
              'time(ps)', ' eng_pv', 'temp_rot', 'vir_cfg', 'vir_src', 'vir_cou', 'vir_bnd', 'vir_ang', 'vir_con', 'vir_tet'
            Write (messages(4), '(5x,a8,5x,a6,4x,a8,5x,a7,5x,a7,7x,a5,8x,a4,7x,a5,5x,a7,7x,a5)') &
              'cpu  (s)', 'volume', 'temp_shl', 'eng_shl', 'vir_shl', 'alpha', 'beta', 'gamma', 'vir_pmf', 'press'
            Write (messages(5), '(a)') Repeat('-', 130)
            Call info(messages, 5, .true.)
          End If

          Write (messages(1), '(i13,1p,9e12.4)') flow%step, stat%stpval(1:9)
          Write (messages(2), '(f13.5,1p,9e12.4)') flow%time, stat%stpval(10:18)
          Write (messages(3), '(0p,f13.3,1p,9e12.4)') tmr%elapsed, stat%stpval(19:27)
          Write (messages(4), '(a)') ''
          Call info(messages, 4, .true.)

          Write (messages(1), '(6x,a7,1p,9e12.4)') 'rolling', stat%ravval(1:9)
          Write (messages(2), '(5x,a8,1p,9e12.4)') 'averages', stat%ravval(10:18)
          Write (messages(3), '(13x,9e12.4)') stat%ravval(19:27)
          Write (messages(4), '(a)') Repeat('-', 130)
          Call info(messages, 4, .true.)

          If (nstph /= 0) Then
            Call flow%line_printed()
          End If

          ! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

          If (traj%ltraj) Call trajectory_write(flow%restart_key, flow%step, thermo%tstep, flow%time, &
                                                io, stat%rsd, netcdf, cnfig, traj, files, comm)
          If (dfcts(1)%ldef) Then
            Call defects_write(flow%restart_key, thermo%ensemble, flow%step, thermo%tstep, flow%time, io, cshell, &
                               dfcts(1), neigh, sites, netcdf, domain, cnfig, files, comm)
            If (dfcts(2)%ldef) Then
              Call defects_write(flow%restart_key, thermo%ensemble, flow%step, thermo%tstep, flow%time, io, cshell, &
                                 dfcts(2), neigh, sites, netcdf, domain, cnfig, files, comm)
            End If
          End If
          If (msd_data%l_msd) Then
            Call msd_write(cnfig, flow%restart_key, cnfig%megatm, flow%step, thermo%tstep, flow%time, stat%stpval, &
                           sites%dof_site, io, msd_data, files, comm)
          End If
          If (rsdc%lrsd) Then
            Call rsd_write(flow%restart_key, flow%step, thermo%tstep, io, rsdc, &
                           flow%time, cshell, stat%rsd, cnfig, comm)
          End If
          If (green%samp > 0) Then
            Call vaf_write(cnfig, flow%restart_key, nstph, tsths, green, sites, comm)
          End If

          ! Save restart data in event of system crash

          If (Mod(nstph, flow%freq_restart) == 0 .and. nstph /= flow%run_steps .and. (.not. devel%l_tor)) Then
            Call system_revive(neigh%cutoff, flow%step, flow%time, sites, io, flow%start_time, &
                               stat, devel, green, thermo, bond, angle, dihedral, inversion, zdensity, &
                               rdf, netcdf, cnfig, files, comm)
          End If

          ! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

          i = 20
          If (Mod(nstph, i) == 0 .or. &
              (tmr%elapsed > Real(i * 60, wp) .and. &
               tmr%elapsed - Real(((Int(tmr%elapsed) / (i * 60)) * i * 60), wp) < &
               tmr%elapsed / Real(nstph, wp))) Then

            If (comm%idnode == 0) Then
              Inquire (File=files(FILE_OUTPUT)%filename, Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                Call files(FILE_OUTPUT)%close ()
                Open (Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Position='append')
              End If
            End If

          End If
        End If

        ! Save last frame positions (for estimates of MSD when cnfig%levcfg==0)

        Do i = 1, cnfig%natms
          stat%xin(i) = cnfig%parts(i)%xxx
          stat%yin(i) = cnfig%parts(i)%yyy
          stat%zin(i) = cnfig%parts(i)%zzz
        End Do
        stat%clin = cnfig%cell
      Else
        Exit
      End If
    End Do

    ! Finish with grace

    If (exout > 0) Then ! normal exit

      ! recover connectivity arrays for REVCON, REVIVE and printing purposes
      ! read_history MUST NOT initialise R,V,F arrays!!!

      cnfig%ltg(1:cnfig%natms) = stat%ltg0(1:cnfig%natms)
      cnfig%lsa(1:cnfig%natms) = stat%lsa0(1:cnfig%natms)
      cnfig%lsi(1:cnfig%natms) = stat%lsi0(1:cnfig%natms)

    Else If (exout < 0) Then ! abnormal exit

      ! If reading HISTORY finished awkwardly
      ! recover positions and generate kinetics

      Do i = 1, cnfig%natms
        cnfig%parts(i)%xxx = stat%xin(i)
        cnfig%parts(i)%yyy = stat%yin(i)
        cnfig%parts(i)%zzz = stat%zin(i)
      End Do
      cnfig%cell = stat%clin

      Call set_temperature(flow%restart_key, flow%step, flow%run_steps, &
                           stat%engrot, sites%dof_site, cshell, stat, cons, pmf, &
                           thermo, minim, rigid, domain, cnfig, seed, comm)

    End If
    Call stat%clean_connect()

    ! Save restart data because of next action (and disallow the same in dl_poly)

    If (.not. devel%l_tor) Then
      Call system_revive(neigh%cutoff, flow%step, flow%time, sites, io, flow%start_time, stat, &
                         devel, green, thermo, bond, angle, dihedral, inversion, zdensity, rdf, netcdf, &
                         cnfig, files, comm)
    End If

    ! step counter is data counter now, so statistics_result is triggered

    flow%step = nstph
    thermo%tstep = tsths

    !!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

  End Subroutine replay_historf

End Module drivers
