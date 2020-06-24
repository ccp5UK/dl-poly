Module meta
!> meta-simulation routines
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J. Madge October 2018
!> contrib - a.m.elena march 2019   updated deallocate uniform routine
!> contrib - i. scivetti march 2020 modidfiction for EVB

  Use analysis,                           Only: analysis_result
  Use angles,                             Only: angles_type
  Use angular_distribution,               Only: adf_type
  Use bonds,                              Only: bonds_type
  Use bounds,                             Only: set_bounds,&
                                                set_bounds_new
  Use build_book,                         Only: build_book_intra
  Use build_chrm,                         Only: build_chrm_intra
  Use build_excl,                         Only: build_excl_intra
  Use build_tplg,                         Only: build_tplg_intra
  Use comms,                              Only: comms_type,&
                                                gsum,&
                                                gsync,&
                                                gtime,&
                                                exit_comms,&
                                                root_id
  Use configuration,                      Only: check_config,&
                                                scan_config,&
                                                configuration_type,&
                                                freeze_atoms,&
                                                origin_config,&
                                                scale_config, &
                                                write_config
  Use constants,                          Only: DLP_RELEASE,&
                                                DLP_VERSION
  Use constraints,                        Only: constraints_type
  Use control,                            Only: read_control,&
                                                scan_control_io,&
                                                scan_control_output
  Use coord,                              Only: coord_type
  Use core_shell,                         Only: core_shell_type
  Use defects,                            Only: defects_type
  Use deport_data,                        Only: mpoles_rotmat_set_halo
  Use development,                        Only: build_info,&
                                                development_type,&
                                                scan_development
  Use dihedrals,                          Only: dihedrals_type
  Use domains,                            Only: domains_type
  Use drivers,                            Only: md_vv,&
                                                replay_historf,&
                                                replay_history
  Use electrostatic,                      Only: ELECTROSTATIC_EWALD,&
                                                electrostatic_type
  Use errors_warnings,                    Only: error,&
                                                warning,&
                                                info,&
                                                init_error_system,&
                                                get_print_level,&
                                                set_print_level,&
                                                check_print_level
  Use evb,                                Only: print_evb_banner
  Use ewald,                              Only: ewald_type
  Use external_field,                     Only: external_field_type
  Use ffield,                             Only: read_field,&
                                                scan_field,&
                                                report_topology
  Use filename,                           Only: FILENAME_SIZE,&
                                                FILE_CONTROL,&
                                                FILE_CURRENT,&
                                                FILE_KPOINTS,&
                                                FILE_OUTPUT,&
                                                FILE_REVCON,&
                                                FILE_REVCON_2,&
                                                FILE_REVCON_3,&
                                                FILE_STATS,&
                                                default_filenames,&
                                                file_type
  Use flow_control,                       Only: flow_type
  Use four_body,                          Only: four_body_type
  Use greenkubo,                          Only: greenkubo_type
  Use halo,                               Only: set_halo_particles
  Use impacts,                            Only: impact_type
  Use inversions,                         Only: inversions_type
  Use io,                                 Only: io_type
  Use, Intrinsic :: iso_fortran_env,      Only: error_unit
  Use kim,                                Only: kim_citations,&
                                                kim_setup,&
                                                kim_type
  Use kinds,                              Only: wi,&
                                                wp
  Use kinetics,                           Only: cap_forces
  Use metal,                              Only: metal_type
  Use minimise,                           Only: minimise_type
  Use mpole,                              Only: POLARISATION_CHARMM,&
                                                mpole_type
  Use msd,                                Only: msd_type
  Use neighbours,                         Only: neighbours_type
  Use netcdf_wrap,                        Only: netcdf_param
  Use numerics,                           Only: seed_type
  Use plumed,                             Only: plumed_finalize,&
                                                plumed_init,&
                                                plumed_type
  Use pmf,                                Only: pmf_type
  Use poisson,                            Only: poisson_type
  Use rdfs,                               Only: rdf_type
  Use rigid_bodies,                       Only: rigid_bodies_type
  Use rsds,                               Only: rsd_type
  Use site,                               Only: site_type
  Use statistics,                         Only: statistics_result,&
                                                stats_type
  Use system,                             Only: system_expand,&
                                                system_init,&
                                                system_revive
  Use temperature,                        Only: set_temperature
  Use tersoff,                            Only: tersoff_type
  Use test_configuration,                 Only: run_configuration_tests
  Use tethers,                            Only: tethers_type
  Use thermostat,                         Only: thermostat_type
  Use three_body,                         Only: threebody_type
  Use timer,                              Only: init_timer_system,&
                                                start_timer,&
                                                stop_timer,&
                                                time_elapsed,&
                                                timer_report,&
                                                timer_type
  Use trajectory,                         Only: trajectory_type,&
                                                trajectory_write
  Use ttm,                                Only: allocate_ttm_arrays,&
                                                ttm_system_init,&
                                                ttm_system_revive,&
                                                ttm_table_read,&
                                                ttm_table_scan,&
                                                ttm_type
  Use ttm_track,                          Only: ttm_ion_temperature
  Use ttm_utils,                          Only: peakProfiler,&
                                                peakProfilerElec,&
                                                printElecLatticeStatsToFile,&
                                                printLatticeStatsToFile
  Use vdw,                                Only: vdw_type
  Use z_density,                          Only: z_density_type

  Use new_control, Only : read_new_control,&
       read_io, &
       read_devel, &
       read_units, &
       read_ttm, &
       read_ensemble, &
       read_bond_analysis, &
       read_structure_analysis, &
       read_forcefield,&
       read_run_parameters,&
       read_system_parameters,&
       write_parameters


  ! HACK
  Use control_parameter_module, Only : parameters_hash_table
  Use new_control_old_style, Only : setup_file_io, scan_new_control_output_old, &
       read_new_control_old

  Implicit None
  Private

#if WIN32
  Character(Len=*), Parameter :: null_unit = 'NUL'
#elif UNIX
  Character(Len=*), Parameter :: null_unit = '/dev/null'
#endif

  Public :: molecular_dynamics

Contains

  !>  MD simulation (including EVB)
  Subroutine molecular_dynamics(dlp_world, thermo, ewld, tmr, devel, stats, &
                                    green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, &
                                    tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, &
                                    fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, flow, &
                                    seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
                                    output_filename, crd, adf)

    Type(comms_type),                       Intent(InOut) :: dlp_world(0:)
    Type(thermostat_type), Allocatable,     Intent(InOut) :: thermo(:)
    Type(ewald_type), Allocatable,          Intent(InOut) :: ewld(:)
    Type(timer_type), Allocatable,          Intent(InOut) :: tmr(:)
    Type(development_type), Allocatable,    Intent(InOut) :: devel(:)
    Type(stats_type), Allocatable,          Intent(InOut) :: stats(:)
    Type(greenkubo_type), Allocatable,      Intent(InOut) :: green(:)
    Type(plumed_type), Allocatable,         Intent(InOut) :: plume(:)
    Type(msd_type), Allocatable,            Intent(InOut) :: msd_data(:)
    Type(metal_type), Allocatable,          Intent(InOut) :: met(:)
    Type(poisson_type), Allocatable,        Intent(InOut) :: pois(:)
    Type(impact_type), Allocatable,         Intent(InOut) :: impa(:)
    Type(defects_type), Allocatable,        Intent(InOut) :: dfcts(:, :)
    Type(bonds_type), Allocatable,          Intent(InOut) :: bond(:)
    Type(angles_type), Allocatable,         Intent(InOut) :: angle(:)
    Type(dihedrals_type), Allocatable,      Intent(InOut) :: dihedral(:)
    Type(inversions_type), Allocatable,     Intent(InOut) :: inversion(:)
    Type(tethers_type), Allocatable,        Intent(InOut) :: tether(:)
    Type(threebody_type), Allocatable,      Intent(InOut) :: threebody(:)
    Type(z_density_type), Allocatable,      Intent(InOut) :: zdensity(:)
    Type(constraints_type), Allocatable,    Intent(InOut) :: cons(:)
    Type(neighbours_type), Allocatable,     Intent(InOut) :: neigh(:)
    Type(pmf_type), Allocatable,            Intent(InOut) :: pmfs(:)
    Type(site_type), Allocatable,           Intent(InOut) :: sites(:)
    Type(core_shell_type), Allocatable,     Intent(InOut) :: core_shells(:)
    Type(vdw_type), Allocatable,            Intent(InOut) :: vdws(:)
    Type(tersoff_type), Allocatable,        Intent(InOut) :: tersoffs(:)
    Type(four_body_type), Allocatable,      Intent(InOut) :: fourbody(:)
    Type(rdf_type), Allocatable,            Intent(InOut) :: rdf(:)
    Type(netcdf_param), Allocatable,        Intent(InOut) :: netcdf(:)
    Type(minimise_type), Allocatable,       Intent(InOut) :: minim(:)
    Type(mpole_type), Allocatable,          Intent(InOut) :: mpoles(:)
    Type(external_field_type), Allocatable, Intent(InOut) :: ext_field(:)
    Type(rigid_bodies_type), Allocatable,   Intent(InOut) :: rigid(:)
    Type(electrostatic_type), Allocatable,  Intent(InOut) :: electro(:)
    Type(domains_type), Allocatable,        Intent(InOut) :: domain(:)
    Type(flow_type), Allocatable,           Intent(InOut) :: flow(:)
    Type(seed_type), Allocatable,           Intent(InOut) :: seed(:)
    Type(trajectory_type), Allocatable,     Intent(InOut) :: traj(:)
    Type(kim_type), Allocatable, Target,    Intent(InOut) :: kim_data(:)
    Type(configuration_type), Allocatable,  Intent(InOut) :: config(:)
    Type(io_type), Allocatable,             Intent(InOut) :: ios(:)
    Type(ttm_type), Allocatable,            Intent(InOut) :: ttms(:)
    Type(rsd_type), Allocatable, Target,    Intent(InOut) :: rsdsc(:)
    Type(file_type), Allocatable,           Intent(InOut) :: files(:, :)
    Character(len=1024),       Intent(In   ) :: control_filename
    Character(len=1024),       Intent(In   ) :: output_filename
    Type(coord_type), Allocatable,          Intent(InOut) :: crd(:)
    Type(adf_type), Allocatable,            Intent(InOut) :: adf(:)

    Type(comms_type) :: comm

    ! Allocate type arrays
    Call allocate_types_uniform(flow(1)%NUM_FF, thermo, ewld,tmr, devel, stats, &
                                green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, &
                                tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, &
                                fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro,domain, &
                                seed, traj, kim_data, config, ios, ttms, rsdsc, files, crd, adf)

    comm=dlp_world(0) ! this shall vanish asap w_ are proper things

    Call molecular_dynamics_driver(dlp_world(0:), comm, thermo, ewld, &
                                       tmr(1), devel(1), stats, green, plume, msd_data, met, pois, &
                                       impa(1), dfcts(1,:), bond, angle, dihedral, inversion, tether, &
                                       threebody, zdensity, cons, neigh, pmfs, sites, &
                                       core_shells, vdws, tersoffs, fourbody, rdf, netcdf(1), &
                                       minim, mpoles, ext_field, rigid, electro, domain, flow(1), &
                                       seed(1), traj(1), kim_data, config, ios(1), ttms, rsdsc, files(1,:), &
                                       output_filename, control_filename, crd, adf)

    Call deallocate_types_uniform     (thermo, ewld, tmr, devel, stats, &
                                      green, plume, msd_data, met,pois, impa, dfcts, bond, angle, dihedral, inversion, &
                                      tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, &
                                      fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, &
                                      seed, traj, kim_data, config, ios, ttms, rsdsc, files, crd, adf)

  End Subroutine molecular_dynamics

  Subroutine molecular_dynamics_initialise(dlp_world, comm, thermo, ewld, tmr, devel, &
       stats, green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, &
       inversion, tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, &
       vdws, tersoffs, fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, &
       domain, flow, seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
       output_filename, crd, adf)

    Type(comms_type),          Intent(InOut) :: dlp_world(0:), comm
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(timer_type),          Intent(InOut) :: tmr
    Type(development_type),    Intent(InOut) :: devel
    Type(stats_type),          Intent(InOut) :: stats
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(plumed_type),         Intent(InOut) :: plume
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(metal_type),          Intent(InOut) :: met
    Type(poisson_type),        Intent(InOut) :: pois
    Type(impact_type),         Intent(InOut) :: impa
    Type(defects_type),        Intent(InOut) :: dfcts(2)
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(constraints_type),    Intent(InOut) :: cons
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(pmf_type),            Intent(InOut) :: pmfs
    Type(site_type),           Intent(InOut) :: sites
    Type(core_shell_type),     Intent(InOut) :: core_shells
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(netcdf_param),        Intent(InOut) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(InOut) :: domain
    Type(flow_type),           Intent(InOut) :: flow
    Type(seed_type),           Intent(InOut) :: seed
    Type(trajectory_type),     Intent(InOut) :: traj
    Type(kim_type), Target,    Intent(InOut) :: kim_data
    Type(configuration_type),  Intent(InOut) :: config
    Type(io_type),             Intent(InOut) :: ios
    Type(ttm_type),            Intent(InOut) :: ttms
    Type(rsd_type), Target,    Intent(InOut) :: rsdsc
    Type(file_type),           Intent(InOut) :: files(FILENAME_SIZE)
    Type(coord_type),          Intent(InOut) :: crd
    Type(adf_type),            Intent(InOut) :: adf
    Character(len=1024),       Intent(In   ) :: control_filename
    Character(len=1024),       Intent(In   ) :: output_filename

    Integer            :: megatm, mtangl, mtbond, mtcons, mtdihd, mtinv, mtrgd, &
         mtshl, mtteth
    Integer, Dimension(3) :: link_cell
    Type( parameters_hash_table )            :: params
    Real(Kind=wp) :: xhi, yhi, zhi
    Integer :: i, ifile, ierr
    Logical :: can_parse

    ! Rename control file if argument was passed
    If (Len_Trim(control_filename) > 0 ) Then
       Call files(FILE_CONTROL)%rename(control_filename)
    Else
       Call files(FILE_CONTROL)%rename('CONTROL')
    End If


    ! Temporary error system
    Call init_error_system(error_unit, dlp_world(0))
    call read_new_control(files(FILE_CONTROL), params, comm, can_parse)

    ! Cannot read as new style
    if (.not. can_parse) then
       devel%new_control = .false.

       !! Enable when new becomes standard
       ! call warning('Control file '//trim(files(FILE_CONTROL)%filename)//' is in old style', .true.)
       ! call warning('Please update, as this will be deprecated in future releases', .true.)

       call molecular_dynamics_initialise_old(dlp_world, comm, thermo, ewld, tmr, devel, &
       stats, green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, &
       inversion, tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, &
       vdws, tersoffs, fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, &
       domain, flow, seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
       output_filename, crd, adf)
       return
    end if

    ! Setup io immediately
    call read_io(params, ios, netcdf, files, comm)
    call read_devel(params, devel, tmr, seed)
    call read_units(params)

    if (output_filename /= "") files(FILE_OUTPUT)%filename = output_filename

    do i = 1, FILENAME_SIZE
      select case (files(i)%filename)
      case ("SCREEN")
        files(i)%unit_no = error_unit
      case ("NONE")
        files(i)%filename = null_unit
      end select
    end do

    if (files(FILE_OUTPUT)%unit_no /= error_unit) &
         Open (Newunit=files(FILE_OUTPUT)%unit_no, File=trim(files(FILE_OUTPUT)%filename), Status='replace')
    dlp_world(0)%ou = files(FILE_OUTPUT)%unit_no

    Call init_error_system(files(FILE_OUTPUT)%unit_no, dlp_world(0))
    Call print_banner(dlp_world)
    Call info('', .true.)

    Open(newunit = ifile, file='build.info', STATUS='REPLACE', iostat=ierr)
    if (ierr .ne. 0) Call error(0, 'Error opening build.info')
    Call build_info(ifile)
    if (check_print_level(2)) Call build_info()

#ifdef CHRONO
    ! Start main timer
    Call init_timer_system(tmr, files(FILE_OUTPUT)%unit_no, dlp_world(0))
    Call start_timer(tmr, 'Initialisation')
#endif

    ! scan the FIELD file data

    Call scan_field(megatm, sites, neigh%max_exclude, mtshl, &
                    mtcons, mtrgd, mtteth, mtbond, mtangl, mtdihd, mtinv, &
                    ext_field, core_shells, cons, pmfs, met, bond, angle, dihedral, inversion, tether, threebody, &
                    vdws, tersoffs, fourbody, rdf, mpoles, rigid, kim_data, files, electro, comm)

    ! Get imc_r & set config%dvar

    call params%retrieve('density_variance', config%dvar)
    config%dvar = 1.0_wp + config%dvar

    ! scan CONFIG file data

    Call scan_config(config, megatm, config%dvar, config%levcfg, xhi, yhi, zhi, ios, domain, files, comm)

    ! scan CONTROL file data

    call read_bond_analysis(params, flow, bond, angle, dihedral, inversion, config%mxgana)
    call read_structure_analysis(params, stats, msd_data, rdf, green, zdensity, adf, crd, traj, dfcts, rsdsc)
    call read_forcefield(params, neigh, config, xhi, yhi, zhi, flow, vdws, electro, ewld, mpoles, core_shells, met, &
       kim_Data, bond, threebody, fourbody, tersoffs)
    call read_run_parameters(params, flow, thermo, stats, config%l_ind)
    call read_ttm(params, ttms)
    call read_ensemble(params, thermo, ttms%l_ttm)
    Call read_system_parameters(params, flow, config, thermo, impa, minim, plume, cons, pmfs, ttms%l_ttm)

    ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
    ! (setup and domains)
    Call set_bounds_new(sites, ttms, ios, core_shells, cons, pmfs, stats, green, devel, &
       msd_data, met, bond, angle, dihedral, inversion, tether, threebody, zdensity, &
       neigh, vdws, tersoffs, fourbody, rdf, mpoles, ext_field, &
       rigid, electro, domain, config, ewld, kim_data, files, flow, comm, &
       xhi, yhi, zhi, megatm, mtangl, mtbond, mtcons, mtdihd, mtinv, mtrgd, &
       mtshl, mtteth, link_cell)

    Call write_parameters(ios, netcdf, files, neigh, config, link_cell, flow, stats, thermo, ttms, mpoles, vdws, &
         electro, core_shells, ewld, met, impa, minim, plume, cons, pmfs, bond, angle, dihedral, inversion, &
         msd_data, rdf, green, zdensity, adf, crd, dfcts, traj, rsdsc)

    Call molecular_dynamics_allocate(sites, config, neigh, thermo, vdws, core_shells, cons, pmfs, &
         rigid, tether, bond, angle, dihedral, inversion, mpoles, met, tersoffs, threebody, fourbody, &
         ext_field, rdf, zdensity, stats, green, ttms, domain, ewld, kim_data, comm)

    If (stats%cur%on) Then
      Call config%k%init(files(FILE_KPOINTS)%filename, comm)
      Call stats%cur%init(config%k%n, 200, files(FILE_CURRENT), comm)
    End If

    ! READ SIMULATION FORCE FIELD
    Call read_field(neigh%cutoff, core_shells, pmfs, cons, thermo, met, bond, angle, &
                    dihedral, inversion, tether, threebody, sites, vdws, tersoffs, fourbody, rdf, &
                    mpoles, ext_field, rigid, electro, config, kim_data, files, flow, crd, comm)

    ! If computing rdf errors, we need to initialise the arrays.
    If (rdf%l_errors_jack .or. rdf%l_errors_block) Then
      Call rdf%init_block(flow%run_steps, sites%ntype_atom)
    End If

    ! CHECK MD CONFIGURATION
    Call check_config(config, electro%key, thermo, sites, flow, comm)

#ifdef CHRONO
    Call stop_timer(tmr, 'Initialisation')
#endif

  end Subroutine molecular_dynamics_initialise

  Subroutine molecular_dynamics_initialise_old(dlp_world, comm, thermo, ewld, tmr, devel, &
       stats, green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, &
       inversion, tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, &
       vdws, tersoffs, fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, &
       domain, flow, seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
       output_filename, crd, adf)
    Type(comms_type),          Intent(InOut) :: dlp_world(0:), comm
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(timer_type),          Intent(InOut) :: tmr
    Type(development_type),    Intent(InOut) :: devel
    Type(stats_type),          Intent(InOut) :: stats
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(plumed_type),         Intent(InOut) :: plume
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(metal_type),          Intent(InOut) :: met
    Type(poisson_type),        Intent(InOut) :: pois
    Type(impact_type),         Intent(InOut) :: impa
    Type(defects_type),        Intent(InOut) :: dfcts(2)
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(constraints_type),    Intent(InOut) :: cons
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(pmf_type),            Intent(InOut) :: pmfs
    Type(site_type),           Intent(InOut) :: sites
    Type(core_shell_type),     Intent(InOut) :: core_shells
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(netcdf_param),        Intent(InOut) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(InOut) :: domain
    Type(flow_type),           Intent(InOut) :: flow
    Type(seed_type),           Intent(InOut) :: seed
    Type(trajectory_type),     Intent(InOut) :: traj
    Type(kim_type), Target,    Intent(InOut) :: kim_data
    Type(configuration_type),  Intent(InOut) :: config
    Type(io_type),             Intent(InOut) :: ios
    Type(ttm_type),            Intent(InOut) :: ttms
    Type(rsd_type), Target,    Intent(InOut) :: rsdsc
    Type(file_type),           Intent(InOut) :: files(FILENAME_SIZE)
    Type(coord_type),          Intent(InOut) :: crd
    Type(adf_type),            Intent(InOut) :: adf
    Character(len=1024),       Intent(In   ) :: control_filename
    Character(len=1024),       Intent(In   ) :: output_filename
    Logical            :: lfce

    ! Set default file names
    Call default_filenames(files)
    ! Rename control file if argument was passed
    If (Len_Trim(control_filename) > 0 ) Then
       Call files(FILE_CONTROL)%rename(control_filename)
    End If
    If (Len_Trim(output_filename) > 0 ) Then
       Call files(FILE_OUTPUT)%rename(output_filename)
    End If

    Call scan_development(devel, files, comm)
    ! Open output file, or direct output unit to stderr
    If (.not. devel%l_scr) Then
       Open (Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Status='replace')
    Else
       files(FILE_OUTPUT)%unit_no = error_unit
    End If
    dlp_world(0)%ou = files(FILE_OUTPUT)%unit_no
    Call init_error_system(files(FILE_OUTPUT)%unit_no, dlp_world(0))


    ! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES
    Call scan_control_output(files, comm)

    Call print_banner(dlp_world)

#ifdef DEBUG
    Call build_info()
#endif

    Call scan_control_io(ios, netcdf, files, comm)

#ifdef CHRONO
    ! Start main timer
    Call init_timer_system(tmr, files(FILE_OUTPUT)%unit_no,dlp_world(0))
    Call start_timer(tmr,'Initialisation')
#endif

    ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
    ! (setup and domains)
    Call set_bounds(sites, ttms, ios, core_shells, cons, pmfs, stats, &
      thermo, green, devel, msd_data, met, pois, bond, angle, dihedral, inversion, &
      tether, threebody, zdensity, neigh, vdws, tersoffs, fourbody, rdf, mpoles, ext_field, &
      rigid, electro, domain, config, ewld, kim_data, files, flow, comm)


    Call info('', .true.)
    Call info("*** pre-scanning stage (set_bounds) DONE ***", .true.)
    Call time_elapsed(tmr)

    If(flow%NUM_FF>1)Then
      Call info(' ',.true.)
      Call print_evb_banner(flow%NUM_FF)
    End If

    Do ff=1,flow%NUM_FF
      Call molecular_dynamics_allocate(sites(ff), config(ff), neigh(ff), thermo(ff), vdws(ff), core_shells(ff), cons(ff), pmfs(ff), &
           rigid(ff), tether(ff), bond(ff), angle(ff), dihedral(ff), inversion(ff), mpoles(ff), met(ff), tersoffs(ff), threebody(ff), fourbody(ff), &
           ext_field(ff), rdf(ff), zdensity(ff), stats(ff), green(ff), ttms(ff), domain(ff), ewld(ff), kim_data(ff), comm)
    end Do

    ! READ SIMULATION CONTROL PARAMETERS
    Call read_control(lfce, impa, ttms(1), dfcts, rigid(1), rsdsc(1), core_shells(1), cons(1), pmfs(1), &
         stats(1), thermo(1), green(1), devel, plume(1), msd_data(1), met(1), pois(1), bond(1), angle(1), dihedral(1), &
         inversion(1), zdensity(1), neigh(1), vdws(1), rdf(1), minim(1), mpoles(1), electro(1), ewld(1), &
         seed, traj, files, tmr, config(1), flow, crd(1), adf(1), comm)

    call set_print_level(0)
    Do ff = 2, flow%NUM_FF
      Call read_control(lfce(ff), impa(ff), ttms(ff), dfcts(ff), rigid(ff), rsdsc(ff), core_shells(ff), cons(ff), pmfs(ff), &
           stats(ff), thermo(ff), green(ff), devel(ff), plume(ff), msd_data(ff), met(ff), pois(ff), bond(ff), angle(ff), dihedral(ff), &
           inversion(ff), zdensity(ff), neigh(ff), vdws(ff), rdf(ff), minim(ff), mpoles(ff), electro(ff), ewld(ff), &
           seed(ff), traj(ff), files(ff), tmr(ff), config(ff), flow(ff), crd(ff), adf(ff), comm)
    End Do
    call set_print_level(old_print_level)

    Do ff = 1, flow%NUM_FF
      If (stats(ff)%cur%on) Then
        Call config(ff)%k%init(files(FILE_KPOINTS)%filename, comm)
        Call stats(ff)%cur%init(config(ff)%k%n, 200, files(FILE_CURRENT), comm)
      End If
    End Do

    ! READ SIMULATION FORCE FIELD
    Do ff = 1, flow%NUM_FF
      If(flow%NUM_FF > 1) Then
        write(message,'(i0)') ff
        Call info(" ",.true.)
        Call info("*** DETAILS OF INTERACTIONS FOR FIELD "//trim(message)//" ***",.true.)
      End If
      Call read_field(neigh(ff)%cutoff, core_shells(ff), pmfs(ff), cons(ff), thermo(ff), met(ff), bond(ff), angle(ff), &
        dihedral(ff), inversion(ff), tether(ff), threebody(ff), sites(ff), vdws(ff), tersoffs(ff), fourbody(ff),rdf(ff), &
        mpoles(ff), ext_field(ff), rigid(ff), electro(ff), config(ff), kim_data(ff), files, flow, crd(ff), comm, ff)

      ! If computing rdf errors, we need to initialise the arrays.
      If(rdf(ff)%l_errors_jack .or. rdf(ff)%l_errors_block) then
        Call rdf(ff)%init_block(flow%run_steps, sites(ff)%ntype_atom)
      End If

      ! CHECK MD CONFIGURATION
      Call check_config(config(ff),electro(ff)%key,thermo(ff),sites(ff),flow,comm)
    End Do

#ifdef CHRONO
    Call stop_timer(tmr, 'Initialisation')
#endif

  end Subroutine molecular_dynamics_initialise_old

  Subroutine molecular_dynamics_allocate(sites, config, neigh, thermo, vdws, core_shells, cons, pmfs, &
       & rigid, tether, bond, angle, dihedral, inversion, mpoles, met, tersoffs, threebody, fourbody, &
       & ext_field, rdf, zdensity, stats, green, ttms, domain, ewld, kim_data, comm)
    Type(comms_type),          Intent(InOut) :: comm
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(stats_type),          Intent(InOut) :: stats
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(metal_type),          Intent(InOut) :: met
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(constraints_type),    Intent(InOut) :: cons
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(pmf_type),            Intent(InOut) :: pmfs
    Type(site_type),           Intent(InOut) :: sites
    Type(core_shell_type),     Intent(InOut) :: core_shells
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(domains_type),        Intent(InOut) :: domain
    Type(kim_type), Target,    Intent(InOut) :: kim_data
    Type(configuration_type),  Intent(InOut) :: config
    Type(ttm_type),            Intent(InOut) :: ttms

        ! ALLOCATE SITE & CONFIG
    Call sites%init()
    Call config%init()

    Call neigh%init_list(config%mxatdm)

    ! ALLOCATE DPD ARRAYS
    Call thermo%init_dpd(vdws%max_vdw)

    ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS
    Call core_shells%init(config%mxatdm, sites%mxtmls, config%mxlshp, domain%neighbours)
    Call cons%init(sites%mxtmls, config%mxatdm, config%mxlshp, domain%neighbours)
    Call pmfs%init(sites%mxtmls, config%mxatdm)
    Call rigid%init(config%mxlshp, sites%mxtmls, config%mxatdm, domain%neighbours)
    Call tether%init(sites%mxtmls, config%mxatdm)
    Call bond%init(config%mxatdm, sites%mxtmls)
    Call angle%init(config%mxatdm, sites%mxtmls)
    Call dihedral%init(config%mxatdm, sites%mxtmls)
    Call inversion%init(config%mxatms, sites%mxtmls)
    Call mpoles%init(sites%max_site, neigh%max_exclude, config%mxatdm, ewld%bspline%num_splines, config%mxatms)

    ! ALLOCATE INTER-LIKE INTERACTION ARRAYS
    Call vdws%init()
    Call met%init(config%mxatms, sites%mxatyp)
    Call tersoffs%init(sites%max_site)
    Call threebody%init(sites%max_site)
    Call fourbody%init(sites%max_site)
    Call ext_field%init()

    ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS
    Call rdf%init()
    Call zdensity%init(sites%mxatyp)
    Call stats%init(rigid%max_rigid, config%mxatms, config%mxatdm)
    Call green%init(config%mxatms, sites%mxatyp)

    ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS
    Call allocate_ttm_arrays(ttms, domain, config, comm)
    Call ttm_table_scan(config%mxbuff, ttms, comm)

    ! Setup KIM
    Call kim_setup(kim_data, config%mxatms, config%mxatdm, config%megatm, neigh%max_list, domain%mxbfxp, comm%mxnode)

  end Subroutine molecular_dynamics_allocate

  !> Simple MD driver
  Subroutine molecular_dynamics_driver(dlp_world, comm, thermo, ewld, tmr, devel, &
       stats, green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, &
       inversion, tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, &
       vdws, tersoffs, fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, &
       domain, flow, seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
       output_filename, crd, adf)

    Type(comms_type),          Intent(InOut) :: dlp_world(0:), comm
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(timer_type),          Intent(InOut) :: tmr
    Type(development_type),    Intent(InOut) :: devel
    Type(stats_type),          Intent(InOut) :: stats
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(plumed_type),         Intent(InOut) :: plume
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(metal_type),          Intent(InOut) :: met
    Type(poisson_type),        Intent(InOut) :: pois
    Type(impact_type),         Intent(InOut) :: impa
    Type(defects_type),        Intent(InOut) :: dfcts(2)
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(constraints_type),    Intent(InOut) :: cons
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(pmf_type),            Intent(InOut) :: pmfs
    Type(site_type),           Intent(InOut) :: sites
    Type(core_shell_type),     Intent(InOut) :: core_shells
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(netcdf_param),        Intent(InOut) :: netcdf
    Type(minimise_type),       Intent(InOut) :: minim
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(InOut) :: domain
    Type(flow_type),           Intent(InOut) :: flow
    Type(seed_type),           Intent(InOut) :: seed
    Type(trajectory_type),     Intent(InOut) :: traj
    Type(kim_type), Target,    Intent(InOut) :: kim_data
    Type(configuration_type),  Intent(InOut) :: config
    Type(io_type),             Intent(InOut) :: ios
    Type(ttm_type),            Intent(InOut) :: ttms
    Type(rsd_type), Target,    Intent(InOut) :: rsdsc
    Type(file_type),           Intent(InOut) :: files(FILENAME_SIZE)
    Character(len=1024),       Intent(In   ) :: control_filename
    Character(len=1024),       Intent(In   ) :: output_filename
    Type(coord_type),          Intent(InOut) :: crd
    Type(adf_type),            Intent(InOut) :: adf

    Character(len=256) :: message
    Integer(Kind=wi)   :: vacuum
    Logical            :: lfce

    Call gtime(tmr%elapsed) ! Initialise wall clock time

    call molecular_dynamics_initialise(dlp_world, comm, thermo, ewld, tmr, devel, &
       stats, green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, &
       inversion, tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, &
       vdws, tersoffs, fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, &
       domain, flow, seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
       output_filename, crd, adf)

    Call info('', .true.)
    Call info("*** all reading and connectivity checks DONE ***", .true.)
    Call time_elapsed(tmr)

    If (flow%l_vdw) Then
       If (vdws(1)%l_direct) Then
         Call error(0,"Error l_vdw does not work with vdw direct, remove vdw direct")
       Else
       Call vdws(1)%print(comm)
       Call info("Dumped vdw interaction tables!",.true.)
       Call exit_comms(dlp_world)
       Stop 0
     End If
    End If

    ! devel%l_org: translate CONFIG into CFGORG and exit gracefully
    If (devel%l_org) Then
      Call info('',.true.)
      Call info("*** Translating the MD system along a vector (CONFIG to CFGORG) ***",.true.)

      Do ff=1,flow%NUM_FF
        Call origin_config(config(ff),ios,devel,netcdf,comm)
      End Do

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    ! devel%l_scl: rescale CONFIG to CFGSCL and exit gracefully
    If (devel%l_scl) Then
      Call info('',.true.)
      Call info("*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***",.true.)

      Do ff=1,flow%NUM_FF
        Call scale_config(config(ff),ios,devel,netcdf,comm)
      End Do

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    ! devel%l_his: generate HISTORY and exit gracefully
    If (devel%l_his) Then
      Call info('',.true.)
      Call info("*** Generating a zero timestep HISTORY frame of the MD system ***",.true.)

      Do ff=1,flow%NUM_FF
        Call traj%init(key=0,freq=1,start=0)
        flow%step  = 0                            ! no steps done
        flow%time  = 0.0_wp                       ! time is not relevant
        Call trajectory_write(flow%restart_key,flow%step,thermo(ff)%tstep,flow%time,ios, &
                              stats(ff)%rsd,netcdf,config(ff),traj,files,comm)
      End Do
      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    Do ff=1,flow%NUM_FF
    ! Expand current system if opted for
      If (config(ff)%l_exp) Then
        Call system_expand(flow%strict, neigh(ff)%cutoff, ios, core_shells(ff), &
          cons(ff), bond(ff), angle(ff), dihedral(ff), inversion(ff), sites(ff),  &
          netcdf, rigid(ff), config(ff), files, comm, ff)
      End If

      If(flow%NUM_FF > 1) Then
        write(message,'(i0)') ff
        Call info("*** LONG RANGE INFORMATION FOR FIELD "//trim(message)//" ***",.true.)
      End If
    ! READ REVOLD (thermodynamic and structural data from restart file)
      Call system_init(neigh(ff)%cutoff, flow%restart_key, flow%time, flow%start_time, flow%step, &
        stats(ff), devel, green(ff), thermo(ff), met(ff), bond(ff), angle(ff), dihedral(ff), inversion(ff), &
        zdensity(ff), sites(ff), vdws(ff), rdf(ff), config(ff), files, comm)

    ! SET domain borders and link-config%cells as default for new jobs
    ! exchange atomic data and positions in border regions
      Call set_halo_particles(electro(ff)%key, neigh(ff), sites(ff), mpoles(ff), domain(ff), config(ff), &
        ewld(ff), kim_data(ff), comm)
    End Do

    Call info('',.true.)
    Call info("*** initialisation and haloing DONE ***",.true.)
    Call time_elapsed(tmr)

    ! For any intra-like interaction, construct book keeping arrays and
    ! exclusion arrays for overlapped two-body inter-like interactions
    Do ff = 1, flow%NUM_FF
      If(flow%NUM_FF > 1) Then
        write(message,'(i0)') ff
        Call info("*** TOPOLOGY FOR FIELD "//trim(message)//" ***",.true.)
      End If
      If (flow%book) Then
        Call build_book_intra(flow%strict, flow%print_topology, flow%simulation, &
          flow,core_shells(ff), cons(ff), pmfs(ff), bond(ff), angle(ff), dihedral(ff), inversion(ff), tether(ff), &
          neigh(ff), sites(ff), rigid(ff), domain(ff), config(ff), comm)
        ! Setting newjob_build_book to FALSE if ff == flow%NUM_FF
        If(ff == flow%NUM_FF)Then
          flow%newjob_build_book=.FALSE.
        EndIf
        If (mpoles(ff)%max_mpoles > 0) Then
          Call build_tplg_intra(neigh(ff)%max_exclude, bond(ff), angle(ff), dihedral(ff), inversion(ff), &
            mpoles(ff), config(ff), comm)
        ! multipoles topology for internal coordinate system
          If (mpoles(ff)%key == POLARISATION_CHARMM) Then
            Call build_chrm_intra(neigh(ff)%max_exclude, core_shells(ff), cons(ff), bond(ff), angle(ff), &
              dihedral(ff), inversion(ff), mpoles(ff), rigid(ff), config(ff), comm)
          End If
        ! CHARMM core-shell screened electrostatic induction interactions
        End If
        If (flow%exclusions) Then
          Call build_excl_intra(electro(ff)%lecx, core_shells(ff), cons(ff), bond(ff), angle(ff), dihedral(ff), &
            inversion(ff), neigh(ff), rigid(ff), config(ff), comm)
        End If
      Else
        Call report_topology(config(ff)%megatm, config(ff)%megfrz, config(ff)%atmfre, config(ff)%atmfrz, core_shells(ff), &
                             cons(ff), pmfs(ff), bond(ff), angle(ff), dihedral(ff), inversion(ff), tether(ff), sites(ff), &
                             rigid(ff))

        ! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed
        If (flow%simulation) Then
          Call core_shells(ff)%deallocate_core_shell_tmp_arrays()

          Call cons(ff)%deallocate_constraints_temps()
          Call pmfs(ff)%deallocate_pmf_tmp_arrays()

          Call rigid(ff)%deallocate_temp()

         Call tether(ff)%deallocate_temp()
        End If
      End If
    End Do

    Call info('',.true.)
    Call info("*** bookkeeping DONE ***",.true.)
    Call time_elapsed(tmr)

    Do ff = 1, flow%NUM_FF
      If(flow%NUM_FF > 1) Then
        write(message,'(i0)') ff
        Call info("*** DETAILS OF NEIGHBOUR LIST FOR FIELD "//trim(message)//" ***",.true.)
      End If
      ! set and halo rotational matrices and their infinitesimal rotations
      If (mpoles(ff)%max_mpoles > 0) Then
        Call mpoles_rotmat_set_halo(mpoles(ff),domain(ff),config(ff),comm)
      End If

      ! SET initial system temperature
      Call set_temperature               &
        (flow%restart_key, flow%step, flow%run_steps, &
        stats(ff)%engrot, sites(ff)%dof_site, core_shells(ff), stats(ff), cons(ff), pmfs(ff), thermo(ff), minim(ff), &
        rigid(ff), domain(ff), config(ff), seed, comm)
    End Do

    Call info('',.true.)
    Call info("*** temperature setting DONE ***",.true.)
    Call time_elapsed(tmr)

    ! Read ttm table file and initialise electronic temperature
    ! grid from any available restart file
    Do ff = 1, flow%NUM_FF
      If (ttms(ff)%l_ttm) Then
        Call ttm_table_read(ttms(ff), comm)
        Call ttm_system_init(flow%step, flow%equil_steps, flow%restart_key, 'DUMP_E', flow%time,&
          thermo(ff)%temp, domain(ff), ttms(ff), comm)
      End If

      ! Frozen atoms option
      Call freeze_atoms(config(ff))

      ! Cap forces in equilibration mode
      If (flow%step <= flow%equil_steps .and. flow%force_cap) Call cap_forces(thermo(ff)%temp,config(ff),comm)

      ! PLUMED initialisation or information message
      If (plume(ff)%l_plumed) Call plumed_init(config(ff)%megatm,thermo(ff)%tstep,thermo(ff)%temp,plume(ff),comm)

    End Do

    ! Indicate nodes mapped on vacuum (no particles)
    ! Here the checkis done with natms of CONFIG, later in evb_check_config we check if ntms is the same for all CONFIG files
    vacuum = 0
    If (config(1)%natms == 0) Then
      vacuum = 1
      Call warning('mapped on vacuum (no particles)')
    End If
    Call gsum(comm,vacuum)
    If (vacuum > 0) Then
      Call warning(2, Real(vacuum,wp), Real(comm%mxnode,wp), 0.0_wp)
    End If

    ! start-up time when forces are not recalculated
    s = tmr%elapsed

#ifdef CHRONO
  call start_timer(tmr,'Main Calc')
#endif

    ! Unit testing (in the absence of a unit testing framework)
    If (devel%run_unit_tests) Then

      If (devel%unit_test%configuration) Then
        If(comm%idnode == root_id) Then
          Write(*,*) 'Running unit tests for configuration module'
        Endif
        Call run_configuration_tests()
      End If

      If(comm%idnode == root_id) Write(*,*) 'Unit tests completed'
      Call exit_comms(dlp_world)
      Stop 0
    Endif

    ! Now you can run fast, boy
    If (devel%l_fast) Call gsync(comm, devel%l_fast)

    If (flow%simulation) Then
      Call md_vv(config, ttms, ios, rsdsc(1), flow,core_shells, cons, pmfs, stats, thermo, &
        plume, pois, bond, angle, dihedral, inversion, zdensity(1), neigh, sites, fourbody, rdf, &
        netcdf, mpoles, ext_field, rigid, domain, seed, traj, kim_data, files, tmr, minim, &
        impa, green, ewld, electro, dfcts, msd_data, tersoffs, tether, threebody, vdws, &
        devel, met, crd, adf, comm)
    Else
      If(flow%NUM_FF==1)Then
        If (lfce) Then
          Call replay_historf(config(1), ios,rsdsc(1), flow,core_shells(1), cons(1), pmfs(1), stats(1), &
            thermo(1), plume(1), msd_data(1), bond(1), angle(1), dihedral(1), inversion(1), zdensity(1), neigh(1), &
            sites(1), vdws(1), tersoffs(1), fourbody(1), rdf(1), netcdf, minim(1), mpoles(1), ext_field(1), rigid(1), &
            electro(1), domain(1), seed, traj, kim_data(1), files, dfcts, tmr, tether(1), threebody(1), &
            pois(1), green(1), ewld(1), devel, met(1), crd(1), comm)
        Else
          Call replay_history(config(1), ios, rsdsc(1), flow,core_shells(1), cons(1), pmfs(1), stats(1), &
            thermo(1), msd_data(1), met(1), pois(1), bond(1), angle(1), dihedral(1), inversion(1), zdensity(1), &
            neigh(1), sites(1), vdws(1), rdf(1), netcdf, minim(1), mpoles(1), ext_field(1), rigid(1), electro(1), &
            domain(1), seed, traj, kim_data(1), dfcts, files, tmr, tether(1), green(1), ewld(1), devel, comm)
        EndIf
      Else
        Write(message,'(1x,a)') 'error - replay option is not implemented for the EVB method.'
        Call error(0,message)
      End If
    End If

#ifdef CHRONO
    call stop_timer(tmr,'Main Calc')
    call start_timer(tmr,'Termination')
#endif

    !Close the statis file if we used it.
    If (stats(1)%statis_file_open) Call files(FILE_STATS)%close()

    ! Report termination of the MD simulation
    Write(message,'(3(a,f12.3),a)') 'run terminating... elapsed  cpu time: ', &
      tmr%elapsed , ' sec, job time: ', tmr%job, ' sec, close time: ', tmr%clear_screen, ' sec'
    Call info(message,.true.)

    ! Two-temperature model simulations: calculate final ionic temperatures and
    !print statistics to files (final)
    If (ttms(1)%l_ttm) Then
      Call ttm_ion_temperature (ttms(1), thermo(1), domain(1), config(1), comm)
      Call printElecLatticeStatsToFile('PEAK_E', flow%time, thermo(1)%temp, flow%step, ttms(1)%ttmstats,ttms(1), comm)
      Call peakProfilerElec('LATS_E', flow%step, ttms(1)%ttmtraj, ttms(1), comm)
      Call printLatticeStatsToFile(ttms(1)%tempion, 'PEAK_I', flow%time, flow%step, ttms(1)%ttmstats, ttms(1), comm)
      Call peakProfiler(ttms(1)%tempion, 'LATS_I', flow%step, ttms(1)%ttmtraj, ttms(1), comm)
    End If

    ! Save restart data for real simulations only (final)
    If (flow%simulation .and. (.not.devel%l_tor)) Then
    ! Write REVCON
      Do ff = 1, flow%NUM_FF
        If(ff == 1)Then
          frevc = FILE_REVCON
        Else If(ff == 2)Then
          frevc = FILE_REVCON_2
        Else If(ff == 3)Then
          frevc = FILE_REVCON_3
        End If
        Call write_config(config(ff), files(frevc), 2, flow%step, thermo(ff)%tstep, ios, flow%time, netcdf, comm)
      End Do

     Call system_revive(neigh(1)%cutoff, flow%step, flow%time, sites(1), ios, flow%start_time, stats(1), &
        devel, green(1), thermo(1), bond(1), angle(1), dihedral(1), inversion(1), zdensity(1), rdf(1), netcdf, config(1), &
        files, comm)
      If (ttms(1)%l_ttm)Then
        Call ttm_system_revive ('DUMP_E', flow%step, flow%time, 1, flow%run_steps, ttms(1), comm)
      End If
    End If

    ! Produce summary of simulation
    If (neigh(1)%unconditional_update .and. flow%step > 0) Then
      If (.not.neigh(1)%update) Then ! Include the final skip in skipping statistics
        stats(1)%neighskip(3) = stats(1)%neighskip(2)*stats(1)%neighskip(3)
        stats(1)%neighskip(2) = stats(1)%neighskip(2)+1.0_wp
        stats(1)%neighskip(3) = stats(1)%neighskip(3)/stats(1)%neighskip(2)+stats(1)%neighskip(1)/stats(1)%neighskip(2)
        stats(1)%neighskip(4) = Min(stats(1)%neighskip(1),stats(1)%neighskip(4))
        stats(1)%neighskip(5) = Max(stats(1)%neighskip(1),stats(1)%neighskip(5))
      End If
    End If

    Call statistics_result &
      (config(1), minim(1)%minimise, msd_data(1)%l_msd, &
      flow%run_steps, core_shells(1)%keyshl, cons(1)%megcon, pmfs(1)%megpmf, &
      flow%step,flow%time, flow%start_time, config(1)%mxatdm, neigh(1)%unconditional_update, &
      stats(1), thermo(1), sites(1), comm)

    ! Final anlysis
    Call analysis_result(neigh(1)%cutoff, thermo(1), &
                         bond(1), angle(1), dihedral(1), inversion(1), stats(1), green(1), &
                         zdensity(1), sites(1), rdf(1), config(1), files, comm)

    ! PLUMED finalisation
    If (plume(1)%l_plumed) Call plumed_finalize()

#ifdef CHRONO
    Call stop_timer(tmr,'Termination')
    Call timer_report(tmr, comm)
#endif

    ! Ask for reference in publications

    Call print_citations(electro(1),mpoles(1),ttms(1))

   If (kim_data(1)%active) Then
      Call kim_citations(kim_data(1), comm)
   End If

    ! Get just the one number to compare against

    If (devel%l_eng) Then
      Write(message,'(a,1p,e20.10)') "TOTAL ENERGY: ", stats(1)%stpval(1)
      Call info('',.true.)
      Call info(message,.true.)
    End If

    ! Close output channel

    If (.not.devel%l_scr) Call files(FILE_OUTPUT)%close()

  End Subroutine molecular_dynamics_driver

  !> Allocate all types uniformly, _i.e._ N of every type

  Subroutine allocate_types_uniform(array_size, thermo, ewld, tmr, devel, stats, &
                                    green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, &
                                    tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, &
                                    fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, &
                                    seed, traj, kim_data, config, ios, ttms, rsdsc, files, crd, adf)

    Integer(Kind=wi),                       Intent(In   ) :: array_size
    Type(thermostat_type), Allocatable,     Intent(InOut) :: thermo(:)
    Type(ewald_type), Allocatable,          Intent(InOut) :: ewld(:)
    Type(timer_type), Allocatable,          Intent(InOut) :: tmr(:)
    Type(development_type), Allocatable,    Intent(InOut) :: devel(:)
    Type(stats_type), Allocatable,          Intent(InOut) :: stats(:)
    Type(greenkubo_type), Allocatable,      Intent(InOut) :: green(:)
    Type(plumed_type), Allocatable,         Intent(InOut) :: plume(:)
    Type(msd_type), Allocatable,            Intent(InOut) :: msd_data(:)
    Type(metal_type), Allocatable,          Intent(InOut) :: met(:)
    Type(poisson_type), Allocatable,        Intent(InOut) :: pois(:)
    Type(impact_type), Allocatable,         Intent(InOut) :: impa(:)
    Type(defects_type), Allocatable,        Intent(InOut) :: dfcts(:, :)
    Type(bonds_type), Allocatable,          Intent(InOut) :: bond(:)
    Type(angles_type), Allocatable,         Intent(InOut) :: angle(:)
    Type(dihedrals_type), Allocatable,      Intent(InOut) :: dihedral(:)
    Type(inversions_type), Allocatable,     Intent(InOut) :: inversion(:)
    Type(tethers_type), Allocatable,        Intent(InOut) :: tether(:)
    Type(threebody_type), Allocatable,      Intent(InOut) :: threebody(:)
    Type(z_density_type), Allocatable,      Intent(InOut) :: zdensity(:)
    Type(constraints_type), Allocatable,    Intent(InOut) :: cons(:)
    Type(neighbours_type), Allocatable,     Intent(InOut) :: neigh(:)
    Type(pmf_type), Allocatable,            Intent(InOut) :: pmfs(:)
    Type(site_type), Allocatable,           Intent(InOut) :: sites(:)
    Type(core_shell_type), Allocatable,     Intent(InOut) :: core_shells(:)
    Type(vdw_type), Allocatable,            Intent(InOut) :: vdws(:)
    Type(tersoff_type), Allocatable,        Intent(InOut) :: tersoffs(:)
    Type(four_body_type), Allocatable,      Intent(InOut) :: fourbody(:)
    Type(rdf_type), Allocatable,            Intent(InOut) :: rdf(:)
    Type(netcdf_param), Allocatable,        Intent(InOut) :: netcdf(:)
    Type(minimise_type), Allocatable,       Intent(InOut) :: minim(:)
    Type(mpole_type), Allocatable,          Intent(InOut) :: mpoles(:)
    Type(external_field_type), Allocatable, Intent(InOut) :: ext_field(:)
    Type(rigid_bodies_type), Allocatable,   Intent(InOut) :: rigid(:)
    Type(electrostatic_type), Allocatable,  Intent(InOut) :: electro(:)
    Type(domains_type), Allocatable,        Intent(InOut) :: domain(:)
    Type(seed_type), Allocatable,           Intent(InOut) :: seed(:)
    Type(trajectory_type), Allocatable,     Intent(InOut) :: traj(:)
    Type(kim_type), Allocatable, Target,    Intent(InOut) :: kim_data(:)
    Type(configuration_type), Allocatable,  Intent(InOut) :: config(:)
    Type(io_type), Allocatable,             Intent(InOut) :: ios(:)
    Type(ttm_type), Allocatable,            Intent(InOut) :: ttms(:)
    Type(rsd_type), Allocatable, Target,    Intent(InOut) :: rsdsc(:)
    Type(file_type), Allocatable,           Intent(InOut) :: files(:, :)
    Type(coord_type), Allocatable,          Intent(InOut) :: crd(:)
    Type(adf_type), Allocatable,            Intent(InOut) :: adf(:)

    Allocate (thermo(array_size))
    Allocate (ewld(array_size))
    Allocate (stats(array_size))
    Allocate (green(array_size))
    Allocate (plume(array_size))
    Allocate (msd_data(array_size))
    Allocate (met(array_size))
    Allocate (pois(array_size))
    Allocate (bond(array_size))
    Allocate (angle(array_size))
    Allocate (dihedral(array_size))
    Allocate (inversion(array_size))
    Allocate (tether(array_size))
    Allocate (threebody(array_size))
    Allocate (zdensity(array_size))
    Allocate (cons(array_size))
    Allocate (neigh(array_size))
    Allocate (pmfs(array_size))
    Allocate (sites(array_size))
    Allocate (core_shells(array_size))
    Allocate (vdws(array_size))
    Allocate (tersoffs(array_size))
    Allocate (fourbody(array_size))
    Allocate (rdf(array_size))
    Allocate (minim(array_size))
    Allocate (mpoles(array_size))
    Allocate (ext_field(array_size))
    Allocate (rigid(array_size))
    Allocate (electro(array_size))
    Allocate (domain(array_size))
    Allocate (kim_data(array_size))
    Allocate (config(array_size))
    Allocate (ios(array_size))
    Allocate (ttms(array_size))
    Allocate (rsdsc(array_size))
    Allocate (crd(array_size))
    Allocate (adf(array_size))

    If(array_size == 1)Then
      Allocate(dfcts(array_size,2))
      Allocate(tmr(array_size))
      Allocate(impa(array_size))
      Allocate(devel(array_size))
      Allocate(netcdf(array_size))
      Allocate(seed(array_size))
      Allocate(traj(array_size))
      Allocate(files(array_size,FILENAME_SIZE))
    Else
      Allocate(dfcts(1,2))
      Allocate(tmr(1))
      Allocate(impa(1))
      Allocate(devel(1))
      Allocate(netcdf(1))
      Allocate(seed(1))
      Allocate(traj(1))
      Allocate(files(1,FILENAME_SIZE))
    End If

  End Subroutine allocate_types_uniform

  Subroutine deallocate_types_uniform(thermo, ewld, tmr, devel, stats, &
                                      green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, &
                                      tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, &
                                      fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, &
                                      seed, traj, kim_data, config, ios, ttms, rsdsc, files, crd, adf)

    Type(thermostat_type), Allocatable,     Intent(InOut) :: thermo(:)
    Type(ewald_type), Allocatable,          Intent(InOut) :: ewld(:)
    Type(timer_type), Allocatable,          Intent(InOut) :: tmr(:)
    Type(development_type), Allocatable,    Intent(InOut) :: devel(:)
    Type(stats_type), Allocatable,          Intent(InOut) :: stats(:)
    Type(greenkubo_type), Allocatable,      Intent(InOut) :: green(:)
    Type(plumed_type), Allocatable,         Intent(InOut) :: plume(:)
    Type(msd_type), Allocatable,            Intent(InOut) :: msd_data(:)
    Type(metal_type), Allocatable,          Intent(InOut) :: met(:)
    Type(poisson_type), Allocatable,        Intent(InOut) :: pois(:)
    Type(impact_type), Allocatable,         Intent(InOut) :: impa(:)
    Type(defects_type), Allocatable,        Intent(InOut) :: dfcts(:, :)
    Type(bonds_type), Allocatable,          Intent(InOut) :: bond(:)
    Type(angles_type), Allocatable,         Intent(InOut) :: angle(:)
    Type(dihedrals_type), Allocatable,      Intent(InOut) :: dihedral(:)
    Type(inversions_type), Allocatable,     Intent(InOut) :: inversion(:)
    Type(tethers_type), Allocatable,        Intent(InOut) :: tether(:)
    Type(threebody_type), Allocatable,      Intent(InOut) :: threebody(:)
    Type(z_density_type), Allocatable,      Intent(InOut) :: zdensity(:)
    Type(constraints_type), Allocatable,    Intent(InOut) :: cons(:)
    Type(neighbours_type), Allocatable,     Intent(InOut) :: neigh(:)
    Type(pmf_type), Allocatable,            Intent(InOut) :: pmfs(:)
    Type(site_type), Allocatable,           Intent(InOut) :: sites(:)
    Type(core_shell_type), Allocatable,     Intent(InOut) :: core_shells(:)
    Type(vdw_type), Allocatable,            Intent(InOut) :: vdws(:)
    Type(tersoff_type), Allocatable,        Intent(InOut) :: tersoffs(:)
    Type(four_body_type), Allocatable,      Intent(InOut) :: fourbody(:)
    Type(rdf_type), Allocatable,            Intent(InOut) :: rdf(:)
    Type(netcdf_param), Allocatable,        Intent(InOut) :: netcdf(:)
    Type(minimise_type), Allocatable,       Intent(InOut) :: minim(:)
    Type(mpole_type), Allocatable,          Intent(InOut) :: mpoles(:)
    Type(external_field_type), Allocatable, Intent(InOut) :: ext_field(:)
    Type(rigid_bodies_type), Allocatable,   Intent(InOut) :: rigid(:)
    Type(electrostatic_type), Allocatable,  Intent(InOut) :: electro(:)
    Type(domains_type), Allocatable,        Intent(InOut) :: domain(:)
    Type(seed_type), Allocatable,           Intent(InOut) :: seed(:)
    Type(trajectory_type), Allocatable,     Intent(InOut) :: traj(:)
    Type(kim_type), Allocatable, Target,    Intent(InOut) :: kim_data(:)
    Type(configuration_type), Allocatable,  Intent(InOut) :: config(:)
    Type(io_type), Allocatable,             Intent(InOut) :: ios(:)
    Type(ttm_type), Allocatable,            Intent(InOut) :: ttms(:)
    Type(rsd_type), Allocatable, Target,    Intent(InOut) :: rsdsc(:)
    Type(file_type), Allocatable,           Intent(InOut) :: files(:, :)
    Type(coord_type), Allocatable,          Intent(InOut) :: crd(:)
    Type(adf_type), Allocatable,            Intent(InOut) :: adf(:)

    If (Allocated(angle)) Deallocate (angle)
    If (Allocated(bond)) Deallocate (bond)
    If (Allocated(config)) Deallocate (config)
    If (Allocated(cons)) Deallocate (cons)
    If (Allocated(core_shells)) Deallocate (core_shells)
    If (Allocated(devel)) Deallocate (devel)
    If (Allocated(dfcts)) Deallocate (dfcts)
    If (Allocated(dihedral)) Deallocate (dihedral)
    If (Allocated(domain)) Deallocate (domain)
    If (Allocated(electro)) Deallocate (electro)
    If (Allocated(ewld)) Deallocate (ewld)
    If (Allocated(ext_field)) Deallocate (ext_field)
    If (Allocated(files)) Deallocate (files)
    If (Allocated(fourbody)) Deallocate (fourbody)
    If (Allocated(green)) Deallocate (green)
    If (Allocated(impa)) Deallocate (impa)
    If (Allocated(inversion)) Deallocate (inversion)
    If (Allocated(ios)) Deallocate (ios)
    If (Allocated(kim_data)) Deallocate (kim_data)
    If (Allocated(met)) Deallocate (met)
    If (Allocated(minim)) Deallocate (minim)
    If (Allocated(mpoles)) Deallocate (mpoles)
    If (Allocated(msd_data)) Deallocate (msd_data)
    If (Allocated(neigh)) Deallocate (neigh)
    If (Allocated(netcdf)) Deallocate (netcdf)
    If (Allocated(plume)) Deallocate (plume)
    If (Allocated(pmfs)) Deallocate (pmfs)
    If (Allocated(pois)) Deallocate (pois)
    If (Allocated(rdf)) Deallocate (rdf)
    If (Allocated(rigid)) Deallocate (rigid)
    If (Allocated(rsdsc)) Deallocate (rsdsc)
    If (Allocated(seed)) Deallocate (seed)
    If (Allocated(sites)) Deallocate (sites)
    If (Allocated(stats)) Deallocate (stats)
    If (Allocated(tersoffs)) Deallocate (tersoffs)
    If (Allocated(tether)) Deallocate (tether)
    If (Allocated(thermo)) Deallocate (thermo)
    If (Allocated(threebody)) Deallocate (threebody)
    If (Allocated(tmr)) Deallocate (tmr)
    If (Allocated(traj)) Deallocate (traj)
    If (Allocated(ttms)) Deallocate (ttms)
    If (Allocated(vdws)) Deallocate (vdws)
    If (Allocated(zdensity)) Deallocate (zdensity)
    If (Allocated(crd)) Deallocate (crd)
    If (Allocated(adf)) Deallocate (adf)

  End Subroutine deallocate_types_uniform

  Subroutine print_banner(dlp_world)
    Type(comms_type), Intent(In   ) :: dlp_world(0:)

    Character(Len=*), Parameter :: fmt1 = '(a)', fmt2 = '(a25,a8,a4,a14,a15)', fmt3 = '(a,i10,a)'

    Character(Len=66) :: banner(15)

    Write (banner(1), fmt1) Repeat("*", 66)
    Write (banner(2), fmt1) "*************  stfc/ccp5  program  library  package  ** D ********"
    Write (banner(3), fmt1) "*************  daresbury laboratory general purpose  *** L *******"
    Write (banner(4), fmt1) "**         **  classical molecular dynamics program  **** \ ******"
    Write (banner(5), fmt1) "** DL_POLY **  authors:   i.t.todorov   &   w.smith  ***** P *****"
    Write (banner(6), fmt2) "**         **  version:  ", DLP_VERSION, " /  ", DLP_RELEASE, "  ****** O ****"
    Write (banner(7), fmt3) "*************  execution on  ", dlp_world(0)%mxnode, " process(es)  ******* L ***"
    Write (banner(8), fmt1) "*************  contributors' list:                   ******** Y **"
    Write (banner(9), fmt1) "*************  ------------------------------------  *************"
    Write (banner(10), fmt1) "*************  i.j.bush, h.a.boateng, r.davidchak,   *************"
    Write (banner(11), fmt1) "*************  m.a.seaton, a.v.brukhno, a.m.elena,   *************"
    Write (banner(12), fmt1) "*************  s.l.daraszewicz,g.khara,s.t.murphy    *************"
    Write (banner(13), fmt1) "*************  j.madge,a.b.g.chalk,i.scivetti,       *************"
    Write (banner(14), fmt1) "*************  j.wilkins                             *************"
    Write (banner(15), fmt1) "******************************************************************"
    Call info(banner, 15, .true., level=-1)
  End Subroutine print_banner

  Subroutine print_citations(electro, mpoles, ttms)
    Type(electrostatic_type), Intent(In   ) :: electro
    Type(mpole_type),         Intent(In   ) :: mpoles
    Type(ttm_type),           Intent(In   ) :: ttms

    Character(Len=*), Parameter :: fmt1 = '(a)'

    Character(Len=66) :: banner(14)

    ! Ask for reference in publications

    Call info('', .true.)
    Write (banner(1), fmt1) Repeat("*", 66)
    Write (banner(2), fmt1) '**** Thank you for using the DL_POLY_4 package in your work.  ****'
    Write (banner(3), fmt1) '**** Please, acknowledge our efforts by including the         ****'
    Write (banner(4), fmt1) '**** following references when publishing data obtained using ****'
    Write (banner(5), fmt1) '**** DL_POLY_4:                                               ****'
    Write (banner(6), fmt1) '****   - I.T. Todorov, W. Smith, K. Trachenko & M.T. Dove,    ****'
    Write (banner(7), fmt1) '****     J. Mater. Chem., 16, 1911-1918 (2006),               ****'
    Write (banner(8), fmt1) '****     https://doi.org/10.1039/B517931A                     ****'
    Call info(banner, 8, .true., level=-1)
    If (electro%key == ELECTROSTATIC_EWALD) Then
      Write (banner(1), fmt1) '****   - I.J. Bush, I.T. Todorov & W. Smith,                  ****'
      Write (banner(2), fmt1) '****     Comp. Phys. Commun., 175, 323-329 (2006),            ****'
      Write (banner(3), fmt1) '****     https://doi.org/10.1016/j.cpc.2006.05.001            ****'
      Call info(banner, 3, .true., level=-1)
    End If
    If (mpoles%max_mpoles > 0) Then
      Write (banner(1), fmt1) '****   - H.A. Boateng & I.T. Todorov,                         ****'
      Write (banner(2), fmt1) '****     J. Chem. Phys., 142, 034117 (2015),                  ****'
      Write (banner(3), fmt1) '****     https://doi.org/10.1063/1.4905952                    ****'
      Call info(banner, 3, .true., level=-1)
    End If
    If (ttms%l_ttm) Then
      Write (banner(1), fmt1) '****   - E. Zarkadoula, S.L. Daraszewicz, D.M. Duffy,         ****'
      Write (banner(2), fmt1) '****     M.A. Seaton, I.T. Todorov, K. Nordlund, M.T. Dove &  ****'
      Write (banner(3), fmt1) '****     K. Trachenko                                         ****'
      Write (banner(4), fmt1) '****     J. Phys.: Condens. Matter, 24, 085401 (2014),        ****'
      Write (banner(5), fmt1) '****     https://doi.org/10.1088/0953-8984/26/8/085401        ****'
      Call info(banner, 5, .true., level=-1)
    End If
    Call info(Repeat("*", 66), .true.)
  End Subroutine print_citations


End Module meta
