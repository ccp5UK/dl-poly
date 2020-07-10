Module meta_evb 
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
  Use bounds,                             Only: set_bounds
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
                                                md_vv_evb,&   
                                                replay_historf,&
                                                replay_history
  Use electrostatic,                      Only: ELECTROSTATIC_EWALD,&
                                                electrostatic_type
  Use errors_warnings,                    Only: error,&
                                                info,&
                                                init_error_system,&
                                                warning
  Use ewald,                              Only: ewald_type
  Use external_field,                     Only: external_field_type
  Use ffield,                             Only: read_field,&
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
  Use meta,                               Only: print_citations,& 
                                                print_banner,&
                                                deallocate_types_uniform,&
                                                allocate_types_uniform
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

  Implicit None
  Private

  Public :: evb_molecular_dynamics

Contains

  !>  EVB MD simulation
  Subroutine evb_molecular_dynamics(dlp_world, thermo, ewld, tmr, devel, stats, &
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

    Call evb_molecular_dynamics_driver(dlp_world(0:), comm, thermo, ewld, &
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

  End Subroutine evb_molecular_dynamics
  ! EVB MD driver
  ! This is the equivalent of the subroutine "molecular_dynamics_driver" for the EVB method.
  ! Essentially, each variable (except tmr, impc, defects and ios) is an array of dimension 
  ! equal to the number of force-fields for the EVB method.  
  
  Subroutine evb_molecular_dynamics_driver(dlp_world, comm, thermo, ewld, tmr, devel, &
                                           stats, green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, &
                                           inversion, tether, threebody, zdensity, cons, neigh, pmfs, sites, core_shells, &
                                           vdws, tersoffs, fourbody, rdf, netcdf, minim, mpoles, ext_field, rigid, electro, &
                                           domain, flow,seed, traj, kim_data, config, ios, ttms, rsdsc, files, control_filename, &
                                           output_filename, crd, adf)

    Type(comms_type),            Intent(InOut) :: dlp_world(0:),comm
    Type(thermostat_type),       Intent(InOut) :: thermo(:)
    Type(ewald_type),            Intent(InOut) :: ewld(:)
    Type(timer_type),            Intent(InOut) :: tmr
    Type(development_type),      Intent(InOut) :: devel
    Type(stats_type),            Intent(InOut) :: stats(:)
    Type(greenkubo_type),        Intent(InOut) :: green(:)
    Type(plumed_type),           Intent(InOut) :: plume(:)
    Type(msd_type),              Intent(InOut) :: msd_data(:)
    Type(metal_type),            Intent(InOut) :: met(:)
    Type(poisson_type),          Intent(InOut) :: pois(:)
    Type(impact_type),           Intent(InOut) :: impa
    Type(defects_type),          Intent(InOut) :: dfcts(:)
    Type(bonds_type),            Intent(InOut) :: bond(:)
    Type(angles_type ),          Intent(InOut) :: angle(:)
    Type(dihedrals_type ),       Intent(InOut) :: dihedral(:)
    Type(inversions_type ),      Intent(InOut) :: inversion(:)
    Type(tethers_type ),         Intent(InOut) :: tether(:)
    Type(threebody_type ),       Intent(InOut) :: threebody(:)
    Type(z_density_type ),       Intent(InOut) :: zdensity(:)
    Type(constraints_type ),     Intent(InOut) :: cons(:)
    Type(neighbours_type ),      Intent(InOut) :: neigh(:)
    Type(pmf_type ),             Intent(InOut) :: pmfs(:)
    Type(site_type ),            Intent(InOut) :: sites(:)
    Type(core_shell_type ),      Intent(InOut) :: core_shells(:)
    Type(vdw_type ),             Intent(InOut) :: vdws(:)
    Type(tersoff_type ),         Intent(InOut) :: tersoffs(:)
    Type(four_body_type ),       Intent(InOut) :: fourbody(:)
    Type(rdf_type ),             Intent(InOut) :: rdf(:)
    Type(netcdf_param ),         Intent(InOut) :: netcdf
    Type(minimise_type ),        Intent(InOut) :: minim(:)
    Type(mpole_type ),           Intent(InOut) :: mpoles(:)
    Type(external_field_type ),  Intent(InOut) :: ext_field(:)
    Type(rigid_bodies_type ),    Intent(InOut) :: rigid(:)
    Type(electrostatic_type ),   Intent(InOut) :: electro(:)
    Type(domains_type ),         Intent(InOut) :: domain(:)
    Type(flow_type ),            Intent(InOut) :: flow
    Type(seed_type ),            Intent(InOut) :: seed
    Type(trajectory_type ),      Intent(InOut) :: traj
    Type(kim_type ), Target,     Intent(InOut) :: kim_data(:)
    Type(configuration_type ),   Intent(InOut) :: config(:)
    Type(io_type),               Intent(InOut) :: ios
    Type(ttm_type),              Intent(InOut) :: ttms(:)
    Type(rsd_type ), Target,     Intent(InOut) :: rsdsc(:)
    Type(file_type ),            Intent(InOut) :: files(:)
    Character(Len= 1024),        Intent(In   ) :: output_filename
    Character(Len= 1024),        Intent(In   ) :: control_filename
    Type(coord_type),            Intent(InOut) :: crd(:)
    Type(adf_type),              Intent(InOut) :: adf(:)

    Character(Len=256 ) :: message
    Integer(Kind=wi)    :: vacuum
    Integer(Kind=wi)    :: ff, frevc
    Logical             :: lfce

    Call gtime(tmr%elapsed) ! Initialise wall clock time

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
  
#ifdef CHRONO
    ! Start main timer
    Call init_timer_system(tmr, files(FILE_OUTPUT)%unit_no,dlp_world(0))
    Call start_timer(tmr,'Initialisation')
#endif

    ! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES
    Call scan_control_output(files,comm)

    Call print_banner(dlp_world)

    Call build_info()

    Call scan_control_io(ios,netcdf,files,comm)

    ! Print error message in case there was an error when reading settings evb from CONTROL via read_simtype
    If(flow%evbfail)Then
      Write(message,'(a)') 'error - Either no or wrong specificaton following the keyword "evb" in CONTROL &
                           &for the number of fields to be coupled via. Value MUST an integer larger than 1'      
      call error(0,message)
    End If

    ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
    ! (setup and domains)
    Do ff = 1, flow%NUM_FF
      Call set_bounds ( &
        sites(ff), ttms(ff), ios, core_shells(ff), cons(ff), pmfs(ff), stats(ff), &
        thermo(ff), green(ff), devel, msd_data(ff), met(ff), pois(ff), bond(ff), angle(ff), dihedral(ff), inversion(ff), &
        tether(ff), threebody(ff), zdensity(ff), neigh(ff), vdws(ff), tersoffs(ff), fourbody(ff), rdf(ff), mpoles(ff), & 
        ext_field(ff), rigid(ff), electro(ff), domain(ff), config(ff), ewld(ff), kim_data(ff), files, flow, comm, ff)
    End Do

    Call info('',.true.)
    Call info("*** pre-scanning stage (set_bounds) DONE ***",.true.)
    Call time_elapsed(tmr)

    ! Now the loop is over the number of force fields to be coupled. Variable flow%NUM_FF was assigned in read_simtype 
    Do ff=1,flow%NUM_FF
      ! ALLOCATE SITE & CONFIG
      Call sites(ff)%init()
      Call config(ff)%init()
      Call neigh(ff)%init_list(config(ff)%mxatdm)

      ! ALLOCATE DPD ARRAYS
      Call thermo(ff)%init_dpd(vdws(ff)%max_vdw)
  
      ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS
      Call core_shells(ff)%init(config(ff)%mxatdm, sites(ff)%mxtmls, config(ff)%mxlshp, domain(ff)%neighbours)
      Call cons(ff)%init(sites(ff)%mxtmls, config(ff)%mxatdm, config(ff)%mxlshp, domain(ff)%neighbours)
      Call pmfs(ff)%init(sites(ff)%mxtmls, config(ff)%mxatdm)
      Call rigid(ff)%init(config(ff)%mxlshp, sites(ff)%mxtmls, config(ff)%mxatdm, domain(ff)%neighbours)
      Call tether(ff)%init(sites(ff)%mxtmls, config(ff)%mxatdm)
      Call bond(ff)%init(config(ff)%mxatdm, sites(ff)%mxtmls)
      Call angle(ff)%init(config(ff)%mxatdm, sites(ff)%mxtmls)
      Call dihedral(ff)%init(config(ff)%mxatdm, sites(ff)%mxtmls)
      Call inversion(ff)%init(config(ff)%mxatms, sites(ff)%mxtmls)
      Call mpoles(ff)%init(sites(ff)%max_site, neigh(ff)%max_exclude, config(ff)%mxatdm, &
                          & ewld(ff)%bspline%num_splines, config(ff)%mxatms)
  
      ! ALLOCATE INTER-LIKE INTERACTION ARRAYS
      Call vdws(ff)%init()
      Call met(ff)%init(config(ff)%mxatms, sites(ff)%mxatyp)
      Call tersoffs(ff)%init(sites(ff)%max_site)
      Call threebody(ff)%init(sites(ff)%max_site)
      Call fourbody(ff)%init(sites(ff)%max_site)
  
      ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS
      Call rdf(ff)%init()
      Call zdensity(ff)%init(sites(ff)%mxatyp)
      Call stats(ff)%init(rigid(ff)%max_rigid, config(ff)%mxatms, config(ff)%mxatms)
      Call green(ff)%init(config(ff)%mxatms, sites(ff)%mxatyp)

      ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS
      Call allocate_ttm_arrays(ttms(ff),domain(ff),config(ff),comm)
      Call ttm_table_scan(config(ff)%mxbuff,ttms(ff),comm)
  
      ! Setup KIM
      Call kim_setup(kim_data(ff), config(ff)%mxatms, config(ff)%mxatdm, config(ff)%megatm,& 
        neigh(ff)%max_list, domain(ff)%mxbfxp, comm%mxnode)

      ! External field
      Call ext_field(ff)%init()

    End Do
    
    ! READ SIMULATION CONTROL PARAMETERS
    Do ff = 1, flow%NUM_FF
      Call read_control(lfce, impa, ttms(ff), dfcts, rigid(ff), rsdsc(ff), core_shells(ff), cons(ff), pmfs(ff), &
        stats(ff), thermo(ff), green(ff), devel, plume(ff), msd_data(ff), met(ff), pois(ff), bond(ff), angle(ff), dihedral(ff), &
        inversion(ff), zdensity(ff), neigh(ff), vdws(ff), rdf(ff), minim(ff), mpoles(ff), electro(ff), ewld(ff), &
        seed, traj, files, tmr, config(ff), flow, crd(ff), adf(ff), comm, ff)
    End Do
    
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

    Call info('',.true.)
    Call info("*** all reading and connectivity checks DONE ***",.true.)
    Call time_elapsed(tmr)

#ifdef CHRONO
    Call stop_timer(tmr,'Initialisation')
#endif

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
    Call time_elapsed(tmr)

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
      Call md_vv_evb(config, ttms, ios, rsdsc(1), flow,core_shells, cons, pmfs, stats, thermo, &
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
      bond(1), angle(1), dihedral(1), inversion(1), stats(1), green(1), zdensity(1), sites(1), rdf(1), config(1), comm)

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

  End Subroutine evb_molecular_dynamics_driver

End Module meta_evb
