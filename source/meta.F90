Module meta
!> meta-simulation routines
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J. Madge October 2018
!> contrib - a.m.elena march 2019 updated deallocate uniform routine
  Use, Intrinsic :: iso_fortran_env, Only : error_unit
  Use kinds, Only : wi,wp
  Use comms, Only : comms_type, init_comms, exit_comms, gsync, gtime,gsum
  Use development, Only : development_type,scan_development,build_info
  Use netcdf_wrap, Only : netcdf_param
  Use domains, Only : domains_type
  Use site, Only : site_type
  Use constants, Only : DLP_RELEASE,DLP_VERSION
  Use configuration, Only : configuration_type,check_config, scale_config, origin_config, freeze_atoms
  Use control, Only : read_control,scan_control_output,scan_control_io
  Use neighbours, Only : neighbours_type
  Use core_shell, Only : core_shell_type
  Use pmf, only : pmf_type
  Use rigid_bodies, Only : rigid_bodies_type
  Use minimise, Only : minimise_type
  Use tethers, Only : tethers_type
  Use bonds, Only : bonds_type
  Use angles, Only : angles_type
  Use dihedrals, Only : dihedrals_type
  Use inversions, Only : inversions_type
  Use three_body, Only : threebody_type
  Use mpole, Only : mpole_type,POLARISATION_CHARMM
  Use vdw, Only : vdw_type
  Use metal, Only : metal_type
  Use tersoff, Only : tersoff_type
  Use four_body, Only : four_body_type
  Use kim, Only : kim_type,kim_setup
  Use plumed, Only : plumed_type,plumed_init,plumed_finalize
  Use external_field, Only : external_field_type
  Use rdfs, Only : rdf_type
  Use z_density, Only : z_density_type
  Use statistics, Only : stats_type,statistics_result
  Use greenkubo, Only : greenkubo_type
  Use msd, Only : msd_type
  Use drivers, Only : w_md_vv, w_replay_historf,w_replay_history
  Use errors_warnings, Only : init_error_system,info, warning
  Use ewald, Only : ewald_type
  Use impacts, Only : impact_type
  Use defects, Only : defects_type
  Use halo, Only : set_halo_particles
  Use deport_data, Only : mpoles_rotmat_set_halo
  Use temperature, Only : set_temperature
  Use rsds, Only : rsd_type
  Use trajectory, Only : trajectory_write,trajectory_type
  Use system, Only : system_revive,system_expand,system_init
  Use build_excl, Only : build_excl_intra
  Use build_book, Only : build_book_intra
  Use ffield, Only : read_field,report_topology
  Use bounds, Only : set_bounds
  Use build_tplg, Only : build_tplg_intra
  Use build_chrm, Only : build_chrm_intra
  Use thermostat, Only : thermostat_type
  Use timer, Only  : timer_type, time_elapsed,timer_report, start_timer, stop_timer, init_timer_system
  Use poisson, Only : poisson_type
  Use analysis, Only : analysis_result
  Use constraints, Only : constraints_type
  Use electrostatic, Only : electrostatic_type,ELECTROSTATIC_EWALD
  Use numerics, Only : seed_type
  Use io, Only : io_type
  Use ttm, Only : ttm_type, ttm_system_init,ttm_system_revive,ttm_table_scan,&
    ttm_table_read,allocate_ttm_arrays
  Use ttm_utils, Only : printElecLatticeStatsToFile,printLatticeStatsToFile,&
    peakProfilerElec,peakProfiler
  Use ttm_track, Only : ttm_ion_temperature
  Use filename, Only : file_type,default_filenames,FILE_CONTROL,FILE_OUTPUT, &
    FILE_STATS,FILENAME_SIZE
  Use flow_control, Only : flow_type
  Use kinetics, Only : cap_forces
  Use coord, Only : coord_type
  Implicit None
  Private

  Public :: molecular_dynamics

Contains

  !> A 'simple', single MD simulation
  Subroutine molecular_dynamics(dlp_world,thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain,flow, &
      seed,traj,kim_data,config,ios,ttms,rsdsc,files,control_filename,crd)

    Type(comms_type), Intent(InOut) :: dlp_world(0:)
    Type(thermostat_type), Allocatable, Intent(InOut) :: thermo(:)
    Type(ewald_type), Allocatable, Intent(InOut) :: ewld(:)
    Type(timer_type), Allocatable, Intent(InOut) :: tmr(:)
    Type(development_type), Allocatable, Intent(InOut) :: devel(:)
    Type(stats_type), Allocatable, Intent(InOut) :: stats(:)
    Type(greenkubo_type), Allocatable, Intent(InOut) :: green(:)
    Type(plumed_type), Allocatable, Intent(InOut) :: plume(:)
    Type(msd_type), Allocatable, Intent(InOut) :: msd_data(:)
    Type(metal_type), Allocatable, Intent(InOut) :: met(:)
    Type(poisson_type), Allocatable, Intent(InOut) :: pois(:)
    Type(impact_type), Allocatable, Intent(InOut) :: impa(:)
    Type(defects_type), Allocatable, Intent(InOut) :: dfcts(:,:)
    Type(bonds_type), Allocatable, Intent(InOut) :: bond(:)
    Type( angles_type ), Allocatable, Intent(InOut) :: angle(:)
    Type( dihedrals_type ), Allocatable, Intent(InOut) :: dihedral(:)
    Type( inversions_type ), Allocatable, Intent(InOut) :: inversion(:)
    Type( tethers_type ), Allocatable, Intent(InOut) :: tether(:)
    Type( threebody_type ), Allocatable, Intent(InOut) :: threebody(:)
    Type( z_density_type ), Allocatable, Intent(InOut) :: zdensity(:)
    Type( constraints_type ), Allocatable, Intent(InOut) :: cons(:)
    Type( neighbours_type ), Allocatable, Intent(InOut) :: neigh(:)
    Type( pmf_type ), Allocatable, Intent(InOut) :: pmfs(:)
    Type( site_type ), Allocatable, Intent(InOut) :: sites(:)
    Type( core_shell_type ), Allocatable, Intent(InOut) :: core_shells(:)
    Type( vdw_type ), Allocatable, Intent(InOut) :: vdws(:)
    Type( tersoff_type ), Allocatable, Intent(InOut) :: tersoffs(:)
    Type( four_body_type ), Allocatable, Intent(InOut) :: fourbody(:)
    Type( rdf_type ), Allocatable, Intent(InOut) :: rdf(:)
    Type( netcdf_param ), Allocatable, Intent(InOut) :: netcdf(:)
    Type( minimise_type ), Allocatable, Intent(InOut) :: minim(:)
    Type( mpole_type ), Allocatable, Intent(InOut) :: mpoles(:)
    Type( external_field_type ), Allocatable, Intent(InOut) :: ext_field(:)
    Type( rigid_bodies_type ), Allocatable, Intent(InOut) :: rigid(:)
    Type( electrostatic_type ), Allocatable, Intent(InOut) :: electro(:)
    Type( domains_type ), Allocatable, Intent(InOut) :: domain(:)
    Type( flow_type ), Allocatable, Intent(InOut) :: flow(:)
    Type( seed_type ), Allocatable, Intent(InOut) :: seed(:)
    Type( trajectory_type ), Allocatable, Intent(InOut) :: traj(:)
    Type( kim_type ), Allocatable, Target, Intent(InOut) :: kim_data(:)
    Type( configuration_type ), Allocatable, Intent(InOut) :: config(:)
    Type( io_type), Allocatable, Intent(InOut) :: ios(:)
    Type( ttm_type), Allocatable, Intent(InOut) :: ttms(:)
    Type( rsd_type ), Allocatable, Target, Intent(InOut) :: rsdsc(:)
    Type( file_type ), Allocatable, Intent(InOut) :: files(:,:)
    Type( coord_type), Allocatable, Intent(InOut) :: crd(:)

    Integer( Kind = wi ), Parameter :: TYPE_SIZE = 1
    Type(comms_type) :: comm
    Character( Len = 1024 ) :: control_filename

    ! Allocate type arrays
    Call allocate_types_uniform(TYPE_SIZE,thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
      seed,traj,kim_data,config,ios,ttms,rsdsc,files,crd)

    comm=dlp_world(0) ! this shall vanish asap w_ are proper things

    Call molecular_dynamics_driver(dlp_world(0:),comm,thermo(1),ewld(1), &
      tmr(1),devel(1),stats(1),green(1),plume(1),msd_data(1),met(1),pois(1), &
      impa(1),dfcts(1,:),bond(1),angle(1),dihedral(1),inversion(1),tether(1), &
      threebody(1),zdensity(1),cons(1),neigh(1),pmfs(1),sites(1), &
      core_shells(1),vdws(1),tersoffs(1),fourbody(1),rdf(1),netcdf(1), &
      minim(1),mpoles(1),ext_field(1),rigid(1),electro(1),domain(1),flow(1), &
      seed(1),traj(1),kim_data(1),config(1),ios(1),ttms(1),rsdsc(1),files(1,:), &
      control_filename,crd(1))
    Call deallocate_types_uniform(thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
      seed,traj,kim_data,config,ios,ttms,rsdsc,files)
  End Subroutine molecular_dynamics

  !> Simple MD driver
  Subroutine molecular_dynamics_driver(dlp_world,comm,thermo,ewld,tmr,devel, &
      stats,green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral, &
      inversion,tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells, &
      vdws,tersoffs,fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro, &
      domain,flow,seed,traj,kim_data,config,ios,ttms,rsdsc,files,control_filename,crd)

    Type(comms_type), Intent(InOut) :: dlp_world(0:),comm
    Type(thermostat_type), Intent(InOut) :: thermo
    Type(ewald_type), Intent(InOut) :: ewld
    Type(timer_type), Intent(InOut) :: tmr
    Type(development_type), Intent(InOut) :: devel
    Type(stats_type), Intent(InOut) :: stats
    Type(greenkubo_type), Intent(InOut) :: green
    Type(plumed_type), Intent(InOut) :: plume
    Type(msd_type), Intent(InOut) :: msd_data
    Type(metal_type), Intent(InOut) :: met
    Type(poisson_type), Intent(InOut) :: pois
    Type(impact_type), Intent(InOut) :: impa
    Type(defects_type), Intent(InOut) :: dfcts(2)
    Type(bonds_type), Intent(InOut) :: bond
    Type( angles_type ), Intent(InOut) :: angle
    Type( dihedrals_type ), Intent(InOut) :: dihedral
    Type( inversions_type ), Intent(InOut) :: inversion
    Type( tethers_type ), Intent(InOut) :: tether
    Type( threebody_type ), Intent(InOut) :: threebody
    Type( z_density_type ), Intent(InOut) :: zdensity
    Type( constraints_type ), Intent(InOut) :: cons
    Type( neighbours_type ), Intent(InOut) :: neigh
    Type( pmf_type ), Intent(InOut) :: pmfs
    Type( site_type ), Intent(InOut) :: sites
    Type( core_shell_type ), Intent(InOut) :: core_shells
    Type( vdw_type ), Intent(InOut) :: vdws
    Type( tersoff_type ), Intent(InOut) :: tersoffs
    Type( four_body_type ), Intent(InOut) :: fourbody
    Type( rdf_type ), Intent(InOut) :: rdf
    Type( netcdf_param ), Intent(InOut) :: netcdf
    Type( minimise_type ), Intent(InOut) :: minim
    Type( mpole_type ), Intent(InOut) :: mpoles
    Type( external_field_type ), Intent(InOut) :: ext_field
    Type( rigid_bodies_type ), Intent(InOut) :: rigid
    Type( electrostatic_type ), Intent(InOut) :: electro
    Type( domains_type ), Intent(InOut) :: domain
    Type( flow_type ), Intent(InOut) :: flow
    Type( seed_type ), Intent(InOut) :: seed
    Type( trajectory_type ), Intent(InOut) :: traj
    Type( kim_type ), Target, Intent(InOut) :: kim_data
    Type( configuration_type ), Intent(InOut) :: config
    Type( io_type), Intent(InOut) :: ios
    Type( ttm_type), Intent(InOut) :: ttms
    Type( rsd_type ), Target, Intent(InOut) :: rsdsc
    Type( coord_type), Intent(InOut) :: crd
    Type( file_type ), Intent(InOut) :: files(FILENAME_SIZE)
    character( len = 1024 ), Intent(In) :: control_filename

    character( len = 256 ) :: message

    Integer( Kind = wi ) :: vacuum
    Logical :: lfce

    Call gtime(tmr%elapsed) ! Initialise wall clock time
    
    ! Set default file names
    Call default_filenames(files)
    ! Rename control file if argument was passed
    If (command_argument_count() == 1) Then
      Call files(FILE_CONTROL)%rename(control_filename)
    End If

    Call scan_development(devel,files,comm)
    ! Open output file, or direct output unit to stderr
    If (.not.devel%l_scr) Then
      Open(Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Status='replace')
    Else
      files(FILE_OUTPUT)%unit_no = error_unit
    End If
    dlp_world(0)%ou=files(FILE_OUTPUT)%unit_no
    Call init_error_system(files(FILE_OUTPUT)%unit_no,dlp_world(0))

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

    ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
    ! (setup and domains)
    Call set_bounds ( &
      sites,ttms,ios,core_shells,cons,pmfs,stats, &
      thermo,green,devel,msd_data,met,pois,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,neigh,vdws,tersoffs,fourbody,rdf,mpoles,ext_field, &
      rigid,electro,domain,config,ewld,kim_data,files,flow,comm)

    Call info('',.true.)
    Call info("*** pre-scanning stage (set_bounds) DONE ***",.true.)
    Call time_elapsed(tmr)

    ! ALLOCATE SITE & CONFIG
    Call sites%init(sites%mxtmls,sites%mxatyp)
    Call config%init()
    Call neigh%init_list(config%mxatdm)

    ! ALLOCATE DPD ARRAYS
    Call thermo%init_dpd(vdws%max_vdw)

    ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS
    Call core_shells%init(config%mxatdm,sites%mxtmls,config%mxlshp,domain%neighbours)
    Call cons%init(sites%mxtmls,config%mxatdm,config%mxlshp,domain%neighbours)
    Call pmfs%init(sites%mxtmls,config%mxatdm)
    Call rigid%init(config%mxlshp,sites%mxtmls,config%mxatdm,domain%neighbours)
    Call tether%init(sites%mxtmls,config%mxatdm)
    Call bond%init(config%mxatdm,sites%mxtmls)
    Call angle%init(config%mxatdm,sites%mxtmls)
    Call dihedral%init(config%mxatdm,sites%mxtmls)
    Call inversion%init(config%mxatms,sites%mxtmls)
    Call mpoles%init(sites%max_site,neigh%max_exclude,config%mxatdm,ewld%bspline,config%mxatms)

    ! ALLOCATE INTER-LIKE INTERACTION ARRAYS
    Call vdws%init()
    Call met%init(config%mxatms,sites%mxatyp)
    Call tersoffs%init(sites%max_site)
    Call threebody%init(sites%max_site)
    Call fourbody%init(sites%max_site)
    Call ext_field%init()

    ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS
    Call rdf%init()
    Call zdensity%init(rdf%max_grid,sites%mxatyp)
    Call stats%init(config%mxatms)
    Call green%init(config%mxatms,sites%mxatyp)

    ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS
    Call allocate_ttm_arrays(ttms,domain,config,comm)
    Call ttm_table_scan(config%mxbuff,ttms,comm)

    ! Setup KIM
    Call kim_setup(kim_data,config%mxatms,config%mxatdm,config%megatm,neigh%max_list,domain%mxbfxp,comm%mxnode)

    ! READ SIMULATION CONTROL PARAMETERS
    Call read_control(lfce,impa,ttms,dfcts,rigid,rsdsc,core_shells,cons,pmfs, &
      stats,thermo,green,devel,plume,msd_data,met,pois,bond,angle,dihedral, &
      inversion,zdensity,neigh,vdws,rdf, minim,mpoles,electro,ewld, &
      seed,traj,files,tmr,config,flow,crd,comm)

    ! READ SIMULATION FORCE FIELD
    Call read_field(neigh%cutoff,core_shells,pmfs,cons,thermo,met,bond,angle, &
      dihedral,inversion,tether,threebody,sites,vdws,tersoffs,fourbody,rdf, &
      mpoles,ext_field,rigid,electro,config,kim_data,files,flow,crd,comm)

    ! If computing rdf errors, we need to initialise the arrays.
    If(rdf%l_errors_jack .or. rdf%l_errors_block) then
      Call rdf%init_block(flow%run_steps,sites%ntype_atom)
    End If


    ! CHECK MD CONFIGURATION
    Call check_config(config,electro%key,thermo,sites,flow,comm)

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

      Call origin_config(config,ios,devel,netcdf,comm)

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    ! devel%l_scl: rescale CONFIG to CFGSCL and exit gracefully
    If (devel%l_scl) Then
      Call info('',.true.)
      Call info("*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***",.true.)

      Call scale_config(config,ios,devel,netcdf,comm)

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    ! devel%l_his: generate HISTORY and exit gracefully
    If (devel%l_his) Then
      Call info('',.true.)
      Call info("*** Generating a zero timestep HISTORY frame of the MD system ***",.true.)

      Call traj%init(key=0,freq=1,start=0)
      flow%step  = 0                            ! no steps done
      flow%time  = 0.0_wp                       ! time is not relevant
      Call trajectory_write(flow%restart_key,flow%step,thermo%tstep,flow%time,ios,stats%rsd,netcdf,config,traj,files,comm)

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    ! Expand current system if opted for
    If (config%l_exp) Then
      Call system_expand(flow%strict,neigh%cutoff,ios,core_shells, &
        cons,bond,angle,dihedral,inversion,sites,netcdf,rigid,config,files,comm)
    End If

    ! READ REVOLD (thermodynamic and structural data from restart file)
    Call system_init(neigh%cutoff,flow%restart_key,flow%time,flow%start_time,flow%step, &
      stats,devel,green,thermo,met,bond,angle,dihedral,inversion, &
      zdensity,sites,vdws,rdf,config,files,comm)

    ! SET domain borders and link-config%cells as default for new jobs
    ! exchange atomic data and positions in border regions
    Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,config,ewld,kim_data,comm)

    Call info('',.true.)
    Call info("*** initialisation and haloing DONE ***",.true.)
    Call time_elapsed(tmr)

    ! For any intra-like interaction, construct book keeping arrays and
    ! exclusion arrays for overlapped two-body inter-like interactions
    If (flow%book) Then
      Call build_book_intra(flow%strict,flow%print_topology,flow%simulation, &
        flow,core_shells,cons,pmfs,bond,angle,dihedral,inversion,tether, &
        neigh,sites,rigid,domain,config,comm)
      If (mpoles%max_mpoles > 0) Then
        Call build_tplg_intra(neigh%max_exclude,bond,angle,dihedral,inversion, &
          mpoles,config,comm)
        ! multipoles topology for internal coordinate system
        If (mpoles%key == POLARISATION_CHARMM) Then
          Call build_chrm_intra(neigh%max_exclude,core_shells,cons,bond,angle, &
            dihedral,inversion,mpoles,rigid,config,comm)
        End If
        ! CHARMM core-shell screened electrostatic induction interactions
      End If
      If (flow%exclusions) Then
        Call build_excl_intra(electro%lecx,core_shells,cons,bond,angle,dihedral, &
          inversion,neigh,rigid,config,comm)
      End If
    Else
      Call report_topology(config%megatm,config%megfrz,config%atmfre,config%atmfrz,core_shells,cons, &
        pmfs,bond,angle,dihedral,inversion,tether,sites,rigid)

      ! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed
      If (flow%simulation) Then
        Call core_shells%deallocate_core_shell_tmp_arrays()

        Call cons%deallocate_constraints_temps()
        Call pmfs%deallocate_pmf_tmp_arrays()

        Call rigid%deallocate_temp()

        Call tether%deallocate_temp()
      End If
    End If

    Call info('',.true.)
    Call info("*** bookkeeping DONE ***",.true.)
    Call time_elapsed(tmr)

    ! set and halo rotational matrices and their infinitesimal rotations
    If (mpoles%max_mpoles > 0) Then
      Call mpoles_rotmat_set_halo(mpoles,domain,config,comm)
    End If

    ! SET initial system temperature
    Call set_temperature               &
      (flow%restart_key,flow%step,flow%run_steps, &
      stats%engrot,sites%dof_site,core_shells,stats,cons,pmfs,thermo,minim, &
      rigid,domain,config,seed,comm)

    Call info('',.true.)
    Call info("*** temperature setting DONE ***",.true.)
    Call time_elapsed(tmr)

    ! Read ttm table file and initialise electronic temperature
    ! grid from any available restart file
    If (ttms%l_ttm) Then
      Call ttm_table_read(ttms,comm)
      Call ttm_system_init(flow%step,flow%equil_steps,flow%restart_key,'DUMP_E',flow%time,thermo%temp,domain,ttms,comm)
    End If

    ! Frozen atoms option
    Call freeze_atoms(config)

    ! Cap forces in equilibration mode
    If (flow%step <= flow%equil_steps .and. flow%force_cap) Call cap_forces(thermo%temp,config,comm)

    ! PLUMED initialisation or information message
    If (plume%l_plumed) Call plumed_init(config%megatm,thermo%tstep,thermo%temp,plume,comm)

    ! Print out sample of initial configuration on node zero
    Call print_initial_configuration(config)

    ! Indicate nodes mapped on vacuum (no particles)
    vacuum = 0
    If (config%natms == 0) Then
      vacuum = 1
      Call warning('mapped on vacuum (no particles)')
    End If
    Call gsum(comm,vacuum)
    If (vacuum > 0) Then
      Call warning(2,Real(vacuum,wp),Real(comm%mxnode,wp),0.0_wp)
    End If

    ! start-up time when forces are not recalculated
    Call time_elapsed(tmr)

#ifdef CHRONO
  call start_timer(tmr,'Main Calc')
#endif

    ! Now you can run fast, boy
    If (devel%l_fast) Call gsync(comm,devel%l_fast)

    If (flow%simulation) Then
      Call w_md_vv(config,ttms,ios,rsdsc,flow,core_shells,cons,pmfs,stats,thermo, &
        plume,pois,bond,angle,dihedral,inversion,zdensity,neigh,sites,fourbody,rdf, &
        netcdf,mpoles,ext_field,rigid,domain,seed,traj,kim_data,files,tmr,minim, &
        impa,green,ewld,electro,dfcts,msd_data,tersoffs,tether,threebody,vdws, &
        devel,met,crd,comm)
    Else
      If (lfce) Then
        Call w_replay_historf(config,ios,rsdsc,flow,core_shells,cons,pmfs,stats, &
          thermo,plume,msd_data,bond,angle,dihedral,inversion,zdensity,neigh, &
          sites,vdws,tersoffs,fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid, &
          electro,domain,seed,traj,kim_data,files,dfcts,tmr,tether,threebody, &
          pois,green,ewld,devel,met,crd,comm)
      Else
        Call w_replay_history(config,ios,rsdsc,flow,core_shells,cons,pmfs,stats, &
          thermo,msd_data,met,pois,bond,angle,dihedral,inversion,zdensity,neigh, &
          sites,vdws,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
          seed,traj,kim_data,dfcts,files,tmr,tether,green,ewld,devel,comm)
      End If
    End If


#ifdef CHRONO
    call stop_timer(tmr,'Main Calc')
    call start_timer(tmr,'Termination')
#endif

    !Close the statis file if we used it.
    If (stats%statis_file_open) Call files(FILE_STATS)%close()

    ! Report termination of the MD simulation
    Write(message,'(3(a,f12.3),a)') 'run terminating... elapsed  cpu time: ', &
      tmr%elapsed , ' sec, job time: ', tmr%job, ' sec, close time: ', tmr%clear_screen, ' sec'
    Call info(message,.true.)

    ! Print out sample of final configuration on node zero
    Call print_final_configuration(config)

    ! Two-temperature model simulations: calculate final ionic temperatures and
    !print statistics to files (final)
    If (ttms%l_ttm) Then
      Call ttm_ion_temperature (ttms,thermo,domain,config,comm)
      Call printElecLatticeStatsToFile('PEAK_E', flow%time, thermo%temp, flow%step, ttms%ttmstats,ttms,comm)
      Call peakProfilerElec('LATS_E', flow%step, ttms%ttmtraj,ttms,comm)
      Call printLatticeStatsToFile(ttms%tempion, 'PEAK_I', flow%time, flow%step, ttms%ttmstats,ttms,comm)
      Call peakProfiler(ttms%tempion, 'LATS_I', flow%step, ttms%ttmtraj,ttms,comm)
    End If

    ! Save restart data for real simulations only (final)
    If (flow%simulation .and. (.not.devel%l_tor)) Then
      Call system_revive(neigh%cutoff,flow%step,flow%time,sites,ios,flow%start_time,stats, &
        devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf,config, &
        files,comm)
      If (ttms%l_ttm) Call ttm_system_revive ('DUMP_E',flow%step,flow%time,1,flow%run_steps,ttms,comm)
    End If

    ! Produce summary of simulation
    If (neigh%unconditional_update .and. flow%step > 0) Then
      If (.not.neigh%update) Then ! Include the final skip in skipping statistics
        stats%neighskip(3)=stats%neighskip(2)*stats%neighskip(3)
        stats%neighskip(2)=stats%neighskip(2)+1.0_wp
        stats%neighskip(3)=stats%neighskip(3)/stats%neighskip(2)+stats%neighskip(1)/stats%neighskip(2)
        stats%neighskip(4)=Min(stats%neighskip(1),stats%neighskip(4))
        stats%neighskip(5)=Max(stats%neighskip(1),stats%neighskip(5))
      End If
    End If

    Call statistics_result                                        &
      (config,minim%minimise,msd_data%l_msd, &
      flow%run_steps,core_shells%keyshl,cons%megcon,pmfs%megpmf,              &
      flow%step,flow%time,flow%start_time,config%mxatdm,neigh%unconditional_update,&
      stats,thermo,sites,comm)

    ! Final anlysis
    Call analysis_result(neigh%cutoff,thermo, &
      bond,angle,dihedral,inversion,stats,green,zdensity,sites,rdf,config,comm)

    ! PLUMED finalisation
    If (plume%l_plumed) Call plumed_finalize()

#ifdef CHRONO
    Call stop_timer(tmr,'Termination')
    Call timer_report(tmr, comm)
#endif

    ! Ask for reference in publications

    Call print_citations(electro,mpoles,ttms)

    ! Get just the one number to compare against

    If (devel%l_eng) Then
      Write(message,'(a,1p,e20.10)') "TOTAL ENERGY: ", stats%stpval(1)
      Call info('',.true.)
      Call info(message,.true.)
    End If

    ! Close output channel

    If (.not.devel%l_scr) Call files(FILE_OUTPUT)%close()
  End Subroutine molecular_dynamics_driver

  !> Allocate all types uniformly, _i.e._ N of every type
  Subroutine allocate_types_uniform(array_size,thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
    seed,traj,kim_data,config,ios,ttms,rsdsc,files,crd)

    Integer( Kind = wi ), Intent( In    ) :: array_size
    Type(thermostat_type), Allocatable, Intent(InOut) :: thermo(:)
    Type(ewald_type), Allocatable, Intent(InOut) :: ewld(:)
    Type(timer_type), Allocatable, Intent(InOut) :: tmr(:)
    Type(development_type), Allocatable, Intent(InOut) :: devel(:)
    Type(stats_type), Allocatable, Intent(InOut) :: stats(:)
    Type(greenkubo_type), Allocatable, Intent(InOut) :: green(:)
    Type(plumed_type), Allocatable, Intent(InOut) :: plume(:)
    Type(msd_type), Allocatable, Intent(InOut) :: msd_data(:)
    Type(metal_type), Allocatable, Intent(InOut) :: met(:)
    Type(poisson_type), Allocatable, Intent(InOut) :: pois(:)
    Type(impact_type), Allocatable, Intent(InOut) :: impa(:)
    Type(defects_type), Allocatable, Intent(InOut) :: dfcts(:,:)
    Type(bonds_type), Allocatable, Intent(InOut) :: bond(:)
    Type( angles_type ), Allocatable, Intent(InOut) :: angle(:)
    Type( dihedrals_type ), Allocatable, Intent(InOut) :: dihedral(:)
    Type( inversions_type ), Allocatable, Intent(InOut) :: inversion(:)
    Type( tethers_type ), Allocatable, Intent(InOut) :: tether(:)
    Type( threebody_type ), Allocatable, Intent(InOut) :: threebody(:)
    Type( z_density_type ), Allocatable, Intent(InOut) :: zdensity(:)
    Type( constraints_type ), Allocatable, Intent(InOut) :: cons(:)
    Type( neighbours_type ), Allocatable, Intent(InOut) :: neigh(:)
    Type( pmf_type ), Allocatable, Intent(InOut) :: pmfs(:)
    Type( site_type ), Allocatable, Intent(InOut) :: sites(:)
    Type( core_shell_type ), Allocatable, Intent(InOut) :: core_shells(:)
    Type( vdw_type ), Allocatable, Intent(InOut) :: vdws(:)
    Type( tersoff_type ), Allocatable, Intent(InOut) :: tersoffs(:)
    Type( four_body_type ), Allocatable, Intent(InOut) :: fourbody(:)
    Type( rdf_type ), Allocatable, Intent(InOut) :: rdf(:)
    Type( netcdf_param ), Allocatable, Intent(InOut) :: netcdf(:)
    Type( minimise_type ), Allocatable, Intent(InOut) :: minim(:)
    Type( mpole_type ), Allocatable, Intent(InOut) :: mpoles(:)
    Type( external_field_type ), Allocatable, Intent(InOut) :: ext_field(:)
    Type( rigid_bodies_type ), Allocatable, Intent(InOut) :: rigid(:)
    Type( electrostatic_type ), Allocatable, Intent(InOut) :: electro(:)
    Type( domains_type ), Allocatable, Intent(InOut) :: domain(:)
    Type( seed_type ), Allocatable, Intent(InOut) :: seed(:)
    Type( trajectory_type ), Allocatable, Intent(InOut) :: traj(:)
    Type( kim_type ), Allocatable, Target, Intent(InOut) :: kim_data(:)
    Type( configuration_type ), Allocatable, Intent(InOut) :: config(:)
    Type( io_type), Allocatable, Intent(InOut) :: ios(:)
    Type( ttm_type), Allocatable, Intent(InOut) :: ttms(:)
    Type( rsd_type ), Allocatable, Target, Intent(InOut) :: rsdsc(:)
    Type( file_type ), Allocatable, Intent(InOut) :: files(:,:)
    Type( coord_type), Allocatable, Intent(InOut) :: crd(:)

    Allocate(thermo(array_size))
    Allocate(ewld(array_size))
    Allocate(tmr(array_size))
    Allocate(devel(array_size))
    Allocate(stats(array_size))
    Allocate(green(array_size))
    Allocate(plume(array_size))
    Allocate(msd_data(array_size))
    Allocate(met(array_size))
    Allocate(pois(array_size))
    Allocate(impa(array_size))
    Allocate(dfcts(array_size,2))
    Allocate(bond(array_size))
    Allocate(angle(array_size))
    Allocate(dihedral(array_size))
    Allocate(inversion(array_size))
    Allocate(tether(array_size))
    Allocate(threebody(array_size))
    Allocate(zdensity(array_size))
    Allocate(cons(array_size))
    Allocate(neigh(array_size))
    Allocate(pmfs(array_size))
    Allocate(sites(array_size))
    Allocate(core_shells(array_size))
    Allocate(vdws(array_size))
    Allocate(tersoffs(array_size))
    Allocate(fourbody(array_size))
    Allocate(rdf(array_size))
    Allocate(netcdf(array_size))
    Allocate(minim(array_size))
    Allocate(mpoles(array_size))
    Allocate(ext_field(array_size))
    Allocate(rigid(array_size))
    Allocate(electro(array_size))
    Allocate(domain(array_size))
    !  Allocate(flow(array_size))
    Allocate(seed(array_size))
    Allocate(traj(array_size))
    Allocate(kim_data(array_size))
    Allocate(config(array_size))
    Allocate(ios(array_size))
    Allocate(ttms(array_size))
    Allocate(rsdsc(array_size))
    Allocate(files(array_size,FILENAME_SIZE))
    Allocate(crd(array_size))
  End Subroutine allocate_types_uniform

  Subroutine deallocate_types_uniform(thermo,ewld,tmr,devel,stats, &
    green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
    tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
    fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
    seed,traj,kim_data,config,ios,ttms,rsdsc,files)

    Type( angles_type ), Allocatable, Intent(InOut) :: angle(:)
    Type(bonds_type), Allocatable, Intent(InOut) :: bond(:)
    Type( configuration_type ), Allocatable, Intent(InOut) :: config(:)
    Type( constraints_type ), Allocatable, Intent(InOut) :: cons(:)
    Type( core_shell_type ), Allocatable, Intent(InOut) :: core_shells(:)
    Type(defects_type), Allocatable, Intent(InOut) :: dfcts(:,:)
    Type(development_type), Allocatable, Intent(InOut) :: devel(:)
    Type( dihedrals_type ), Allocatable, Intent(InOut) :: dihedral(:)
    Type( domains_type ), Allocatable, Intent(InOut) :: domain(:)
    Type( electrostatic_type ), Allocatable, Intent(InOut) :: electro(:)
    Type(ewald_type), Allocatable, Intent(InOut) :: ewld(:)
    Type( external_field_type ), Allocatable, Intent(InOut) :: ext_field(:)
    Type( file_type ), Allocatable, Intent(InOut) :: files(:,:)
    Type( four_body_type ), Allocatable, Intent(InOut) :: fourbody(:)
    Type(greenkubo_type), Allocatable, Intent(InOut) :: green(:)
    Type(impact_type), Allocatable, Intent(InOut) :: impa(:)
    Type( inversions_type ), Allocatable, Intent(InOut) :: inversion(:)
    Type( io_type), Allocatable, Intent(InOut) :: ios(:)
    Type( kim_type ), Allocatable, Target, Intent(InOut) :: kim_data(:)
    Type(metal_type), Allocatable, Intent(InOut) :: met(:)
    Type( minimise_type ), Allocatable, Intent(InOut) :: minim(:)
    Type( mpole_type ), Allocatable, Intent(InOut) :: mpoles(:)
    Type(msd_type), Allocatable, Intent(InOut) :: msd_data(:)
    Type( neighbours_type ), Allocatable, Intent(InOut) :: neigh(:)
    Type( netcdf_param ), Allocatable, Intent(InOut) :: netcdf(:)
    Type(plumed_type), Allocatable, Intent(InOut) :: plume(:)
    Type( pmf_type ), Allocatable, Intent(InOut) :: pmfs(:)
    Type(poisson_type), Allocatable, Intent(InOut) :: pois(:)
    Type( rdf_type ), Allocatable, Intent(InOut) :: rdf(:)
    Type( rigid_bodies_type ), Allocatable, Intent(InOut) :: rigid(:)
    Type( rsd_type ), Allocatable, Target, Intent(InOut) :: rsdsc(:)
    Type( seed_type ), Allocatable, Intent(InOut) :: seed(:)
    Type( site_type ), Allocatable, Intent(InOut) :: sites(:)
    Type(stats_type), Allocatable, Intent(InOut) :: stats(:)
    Type( tersoff_type ), Allocatable, Intent(InOut) :: tersoffs(:)
    Type( tethers_type ), Allocatable, Intent(InOut) :: tether(:)
    Type(thermostat_type), Allocatable, Intent(InOut) :: thermo(:)
    Type( threebody_type ), Allocatable, Intent(InOut) :: threebody(:)
    Type(timer_type), Allocatable, Intent(InOut) :: tmr(:)
    Type( trajectory_type ), Allocatable, Intent(InOut) :: traj(:)
    Type( ttm_type), Allocatable, Intent(InOut) :: ttms(:)
    Type( vdw_type ), Allocatable, Intent(InOut) :: vdws(:)
    Type( z_density_type ), Allocatable, Intent(InOut) :: zdensity(:)

    If (Allocated(angle)) Deallocate(angle)
    If (Allocated(bond)) Deallocate(bond)
    If (Allocated(config)) Deallocate(config)
    If (Allocated(cons)) Deallocate(cons)
    If (Allocated(core_shells)) Deallocate(core_shells)
    If (Allocated(devel)) Deallocate(devel)
    If (Allocated(dfcts)) Deallocate(dfcts)
    If (Allocated(dihedral)) Deallocate(dihedral)
    If (Allocated(domain)) Deallocate(domain)
    If (Allocated(electro)) Deallocate(electro)
    If (Allocated(ewld)) Deallocate(ewld)
    If (Allocated(ext_field)) Deallocate(ext_field)
    If (Allocated(files)) Deallocate(files)
    If (Allocated(fourbody)) Deallocate(fourbody)
    If (Allocated(green)) Deallocate(green)
    If (Allocated(impa)) Deallocate(impa)
    If (Allocated(inversion)) Deallocate(inversion)
    If (Allocated(ios)) Deallocate(ios)
    If (Allocated(kim_data)) Deallocate(kim_data)
    If (Allocated(met)) Deallocate(met)
    If (Allocated(minim)) Deallocate(minim)
    If (Allocated(mpoles)) Deallocate(mpoles)
    If (Allocated(msd_data)) Deallocate(msd_data)
    If (Allocated(neigh)) Deallocate(neigh)
    If (Allocated(netcdf)) Deallocate(netcdf)
    If (Allocated(plume)) Deallocate(plume)
    If (Allocated(pmfs)) Deallocate(pmfs)
    If (Allocated(pois)) Deallocate(pois)
    If (Allocated(rdf)) Deallocate(rdf)
    If (Allocated(rigid)) Deallocate(rigid)
    If (Allocated(rsdsc)) Deallocate(rsdsc)
    If (Allocated(seed)) Deallocate(seed)
    If (Allocated(sites)) Deallocate(sites)
    If (Allocated(stats)) Deallocate(stats)
    If (Allocated(tersoffs)) Deallocate(tersoffs)
    If (Allocated(tether)) Deallocate(tether)
    If (Allocated(thermo)) Deallocate(thermo)
    If (Allocated(threebody)) Deallocate(threebody)
    If (Allocated(tmr)) Deallocate(tmr)
    If (Allocated(traj)) Deallocate(traj)
    If (Allocated(ttms)) Deallocate(ttms)
    If (Allocated(vdws)) Deallocate(vdws)
    If (Allocated(zdensity)) Deallocate(zdensity)

  End Subroutine deallocate_types_uniform

  Subroutine print_banner(dlp_world)
    Type(comms_type), Intent(In) :: dlp_world(0:)

    Character(Len=*), Parameter :: fmt1 = '(a)', &
      fmt2 = '(a25,a8,a4,a14,a15)', &
      fmt3 = '(a,i10,a)'
    Character(Len=66)  :: banner(14)

    Write(banner(1),fmt1)  Repeat("*",66)
    Write(banner(2),fmt1)  "*************  stfc/ccp5  program  library  package  ** D ********"
    Write(banner(3),fmt1)  "*************  daresbury laboratory general purpose  *** L *******"
    Write(banner(4),fmt1)  "**         **  classical molecular dynamics program  **** \ ******"
    Write(banner(5),fmt1)  "** DL_POLY **  authors:   i.t.todorov   &   w.smith  ***** P *****"
    Write(banner(6),fmt2)  "**         **  version:  ", DLP_VERSION, " /  ", DLP_RELEASE, "  ****** O ****"
    Write(banner(7),fmt3)  "*************  execution on  ",dlp_world(0)%mxnode," process(es)  ******* L ***"
    Write(banner(8),fmt1)  "*************  contributors' list:                   ******** Y **"
    Write(banner(9),fmt1)  "*************  ------------------------------------  *************"
    Write(banner(10),fmt1) "*************  i.j.bush, h.a.boateng, r.davidchak,   *************"
    Write(banner(11),fmt1) "*************  m.a.seaton, a.v.brukhno, a.m.elena,   *************"
    Write(banner(12),fmt1) "*************  s.l.daraszewicz,g.khara,s.t.murphy    *************"
    Write(banner(13),fmt1) "*************  j.madge,a.b.g.chalk,i.scivetti        *************"
    Write(banner(14),fmt1) "******************************************************************"
    Call info(banner,14,.true.)
  End Subroutine print_banner

  Subroutine print_citations(electro,mpoles,ttms)
    Type(electrostatic_type), Intent(In) :: electro
    Type(mpole_type), Intent(In) :: mpoles
    Type(ttm_type), Intent(In) :: ttms

    Character(Len=*), Parameter :: fmt1 = '(a)'
    Character(Len=66)  :: banner(14)

    Call info('',.true.)
    Write(banner(1),fmt1) Repeat("*",66)
    Write(banner(2),fmt1) "****  Please do cite `J. Mater. Chem.', 16, 1911-1918 (2006)  ****"
    Write(banner(3),fmt1) "****  when publishing research data obtained using DL_POLY_4  ****"
    Write(banner(4),fmt1) Repeat("*",66)
    Call info(banner,4,.true.)

    ! Ask for reference in publications

    Call info('',.true.)
    Write(banner(1),fmt1) Repeat("*",66)
    Write(banner(2),fmt1) '**** Thank you for using the DL_POLY_4 package in your work.  ****'
    Write(banner(3),fmt1) '**** Please, acknowledge our efforts by including the         ****'
    Write(banner(4),fmt1) '**** following references when publishing data obtained using ****'
    Write(banner(5),fmt1) '**** DL_POLY_4:                                               ****'
    Write(banner(6),fmt1) '****   - I.T. Todorov, W. Smith, K. Trachenko & M.T. Dove,    ****'
    Write(banner(7),fmt1) '****     J. Mater. Chem., 16, 1911-1918 (2006),               ****'
    Write(banner(8),fmt1) '****     https://doi.org/10.1039/B517931A                     ****'
    Call info(banner,8,.true.)
    If (electro%key == ELECTROSTATIC_EWALD) Then
      Write(banner(1),fmt1) '****   - I.J. Bush, I.T. Todorov & W. Smith,                  ****'
      Write(banner(2),fmt1) '****     Comp. Phys. Commun., 175, 323-329 (2006),            ****'
      Write(banner(3),fmt1) '****     https://doi.org/10.1016/j.cpc.2006.05.001            ****'
      Call info(banner,3,.true.)
    End If
    If (mpoles%max_mpoles > 0) Then
      Write(banner(1),fmt1) '****   - H.A. Boateng & I.T. Todorov,                         ****'
      Write(banner(2),fmt1) '****     J. Chem. Phys., 142, 034117 (2015),                  ****'
      Write(banner(3),fmt1) '****     https://doi.org/10.1063/1.4905952                    ****'
      Call info(banner,3,.true.)
    End If
    If (ttms%l_ttm) Then
      Write(banner(1),fmt1) '****   - E. Zarkadoula, S.L. Daraszewicz, D.M. Duffy,         ****'
      Write(banner(2),fmt1) '****     M.A. Seaton, I.T. Todorov, K. Nordlund, M.T. Dove &  ****'
      Write(banner(3),fmt1) '****     K. Trachenko                                         ****'
      Write(banner(4),fmt1) '****     J. Phys.: Condens. Matter, 24, 085401 (2014),        ****'
      Call info(banner,4,.true.)
    End If
    Call info(Repeat("*",66),.true.)
  End Subroutine print_citations


  Subroutine print_initial_configuration(config)
    Type(configuration_type), Intent(In) :: config

    ! Print out sample of initial configuration on node zero
    Call info('',.true.)
    Call info('sample of starting configuration on node zero:',.true.)

    Call print_configuration_sample(config)

    Call info('',.true.)
  End Subroutine print_initial_configuration

  Subroutine print_final_configuration(config)
    Type(configuration_type), Intent(In) :: config

    ! Print out sample of initial configuration on node zero
    Call info('',.true.)
    Call info("sample of final configuration on node zero",.true.)

    Call print_configuration_sample(config)

    Call info('',.true.)
  End Subroutine print_final_configuration

  Subroutine print_configuration_sample(config)
    Type(configuration_type), Intent(In) :: config

    Integer(Kind=wi) :: atom, span
    Character(Len=256) :: message

    If (config%levcfg <= 1) Then
      Write(message,'(7x,a1,7x,a4,2(8x,a4),3(7x,a5))') &
        'i', 'x(i)', 'y(i)', 'z(i)', 'vx(i)', 'vy(i)', 'vz(i)'
      Call info(message,.true.)
    End If

    If (config%levcfg == 2) Then
      Write(message,'(7x,a1,7x,a4,2(8x,a4),6(7x,a5))') &
        'i', 'x(i)', 'y(i)', 'z(i)', 'vx(i)', 'vy(i)', 'vz(i)', &
        'fx(i)', 'fy(i)', 'fz(i)'
      Call info(message,.true.)
    End If

    span=(config%natms+19)/20
    If (span > 0) Then
      Do atom=1,config%natms,span
        If (config%levcfg <= 1) Then
          Write(message,'(i8,1p,6e12.4)') config%ltg(atom), &
            config%parts(atom)%xxx,config%parts(atom)%yyy,config%parts(atom)%zzz, &
            config%vxx(atom),config%vyy(atom),config%vzz(atom)
        End If

        If (config%levcfg == 2) Then
          Write(message,"(i8,1p,9e12.4)") config%ltg(atom), &
            config%parts(atom)%xxx,config%parts(atom)%yyy,config%parts(atom)%zzz,&
            config%vxx(atom),config%vyy(atom),config%vzz(atom), &
            config%parts(atom)%fxx,config%parts(atom)%fyy,config%parts(atom)%fzz
        End If
        Call info(message,.true.)
      End Do
    End If
  End Subroutine print_configuration_sample
End Module meta
