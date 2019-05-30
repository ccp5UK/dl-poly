Module evb 
!> evb-simulation routines
!>
!> Copyright - Daresbury Laboratory
!>
!> Author  - i.scivetti  xxx 2019 
!> contrib - a.m.elena   xxx 2019
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
  Use filename, Only : file_type,default_filenames_evb,FILE_CONTROL,FILE_OUTPUT, &
    FILE_STATS,FILENAME_SIZE
  Use flow_control, Only : flow_type
  Use kinetics, Only : cap_forces
  Use meta, Only: print_initial_configuration, print_final_configuration, print_citations, print_banner

  Implicit None
  Private

  Public :: evb_molecular_dynamics

Contains

  !>  EVB MD simulation
  Subroutine evb_molecular_dynamics(dlp_world,thermo,ewld,tmr,devel,stats, &
    green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
    tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
    fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain,flow, &
    seed,traj,kim_data,config,ios,ttms,rsdsc,files,control_filename)

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

    Type(comms_type) :: comm
    Character( Len = 1024 ) :: control_filename

    ! Allocate type arrays
    Call evb_allocate_types_uniform(flow(1)%TYPE_SIZE_FF,thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
      seed,traj,kim_data,config,ios,ttms,rsdsc,files)

    comm=dlp_world(0) ! this shall vanish asap w_ are proper things

    Call evb_molecular_dynamics_driver(dlp_world(0:),comm,thermo,ewld, &
      tmr(1),devel,stats,green,plume,msd_data,met,pois, &
      impa(1),dfcts(1,:),bond,angle,dihedral,inversion,tether, &
      threebody,zdensity,cons,neigh,pmfs,sites, &
      core_shells,vdws,tersoffs,fourbody,rdf,netcdf, &
      minim,mpoles,ext_field,rigid,electro,domain,flow, &
      seed,traj,kim_data,config,ios(1),ttms,rsdsc,files, &
      control_filename)

    Call evb_deallocate_types_uniform(thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
      seed,traj,kim_data,config,ios,ttms,rsdsc,files)

  End Subroutine evb_molecular_dynamics

  !> EVB MD driver
  ! This is the equivalent of the subroutine "molecular_dynamics_driver" for the EVB method.
  ! Essentially, each variable (except tmr, impc, defects and ios) is an array of dimension 
  ! equal to the number of force-fields for the EVB method. As it is, it prints details of the
  ! CONTROL file twice. This should be fixed in a later stage. 
  !
  Subroutine evb_molecular_dynamics_driver(dlp_world,comm,thermo,ewld,tmr,devel, &
    stats,green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral, &
    inversion,tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells, &
    vdws,tersoffs,fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro, &
    domain,flow,seed,traj,kim_data,config,ios,ttms,rsdsc,files,control_filename)

    Type(comms_type), Intent(InOut) :: dlp_world(0:),comm
    Type(thermostat_type), Intent(InOut) :: thermo(:)
    Type(ewald_type), Intent(InOut) :: ewld(:)
    Type(timer_type), Intent(InOut) :: tmr
    Type(development_type), Intent(InOut) :: devel(:)
    Type(stats_type), Intent(InOut) :: stats(:)
    Type(greenkubo_type), Intent(InOut) :: green(:)
    Type(plumed_type), Intent(InOut) :: plume(:)
    Type(msd_type), Intent(InOut) :: msd_data(:)
    Type(metal_type), Intent(InOut) :: met(:)
    Type(poisson_type), Intent(InOut) :: pois(:)
    Type(impact_type), Intent(InOut) :: impa
    Type(defects_type), Intent(InOut) :: dfcts(2)
    Type(bonds_type), Intent(InOut) :: bond(:)
    Type( angles_type ), Intent(InOut) :: angle(:)
    Type( dihedrals_type ), Intent(InOut) :: dihedral(:)
    Type( inversions_type ), Intent(InOut) :: inversion(:)
    Type( tethers_type ), Intent(InOut) :: tether(:)
    Type( threebody_type ), Intent(InOut) :: threebody(:)
    Type( z_density_type ), Intent(InOut) :: zdensity(:)
    Type( constraints_type ), Intent(InOut) :: cons(:)
    Type( neighbours_type ), Intent(InOut) :: neigh(:)
    Type( pmf_type ), Intent(InOut) :: pmfs(:)
    Type( site_type ), Intent(InOut) :: sites(:)
    Type( core_shell_type ), Intent(InOut) :: core_shells(:)
    Type( vdw_type ), Intent(InOut) :: vdws(:)
    Type( tersoff_type ), Intent(InOut) :: tersoffs(:)
    Type( four_body_type ), Intent(InOut) :: fourbody(:)
    Type( rdf_type ), Intent(InOut) :: rdf(:)
    Type( netcdf_param ), Intent(InOut) :: netcdf(:)
    Type( minimise_type ), Intent(InOut) :: minim(:)
    Type( mpole_type ), Intent(InOut) :: mpoles(:)
    Type( external_field_type ), Intent(InOut) :: ext_field(:)
    Type( rigid_bodies_type ), Intent(InOut) :: rigid(:)
    Type( electrostatic_type ), Intent(InOut) :: electro(:)
    Type( domains_type ), Intent(InOut) :: domain(:)
    Type( flow_type ), Intent(InOut) :: flow(:)
    Type( seed_type ), Intent(InOut) :: seed(:)
    Type( trajectory_type ), Intent(InOut) :: traj(:)
    Type( kim_type ), Target, Intent(InOut) :: kim_data(:)
    Type( configuration_type ), Intent(InOut) :: config(:)
    Type( io_type), Intent(InOut) :: ios
    Type( ttm_type), Intent(InOut) :: ttms(:)
    Type( rsd_type ), Target, Intent(InOut) :: rsdsc(:)
    Type( file_type ), Intent(InOut) :: files(:,:)
    character( len = 1024 ), Intent(In) :: control_filename

    character( len = 256 ) :: message

    Integer( Kind = wi ) :: vacuum, ff
    Logical :: lfce


    Call gtime(tmr%elapsed) ! Initialise wall clock time

    ! Set default file names. Here we loop over the numer of fields to be coupled via the EVB method
    Do ff=1,flow(1)%TYPE_SIZE_FF
      Call default_filenames_evb(ff,files(ff,:))
    End Do

    ! Rename control file if argument was passed. No need to consider the loop over other 
    ! fields as there is only one CONTROL file
    If (command_argument_count() == 1) Then
      Do ff=1,flow(1)%TYPE_SIZE_FF      
        Call files(ff,FILE_CONTROL)%rename(control_filename)
      End Do  
    End If

    Do ff=1,flow(1)%TYPE_SIZE_FF
      Call scan_development(devel(ff),files(ff,:),comm)
    End Do

    ! Open output file, or direct output unit to stderr. 
    ! Only one OUTPUT is considered
    If (.not.devel(1)%l_scr) Then
      Open(Newunit=files(1,FILE_OUTPUT)%unit_no, File=files(1,FILE_OUTPUT)%filename, Status='replace')
    Else
      files(1,FILE_OUTPUT)%unit_no = error_unit
    End If

    dlp_world(0)%ou=files(1,FILE_OUTPUT)%unit_no

    Call init_error_system(files(1,FILE_OUTPUT)%unit_no,dlp_world(0))

#ifdef CHRONO
    ! Start main timer
    Call init_timer_system(tmr, files(1,FILE_OUTPUT)%unit_no,dlp_world(0))
    Call start_timer(tmr,'Initialisation')
#endif

    ! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES
    Do ff=1,flow(1)%TYPE_SIZE_FF
      Call scan_control_output(files(ff,:),comm)
    End Do

    Call print_banner(dlp_world)

    Call build_info()

    Do ff=1,flow(1)%TYPE_SIZE_FF
      Call scan_control_io(ios,netcdf(ff),files(ff,:),comm)
    End Do

    ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
    ! (setup and domains)

    Do ff=1,flow(1)%TYPE_SIZE_FF
    Call set_bounds ( &
      sites(ff),ttms(ff),ios,core_shells(ff),cons(ff),pmfs(ff),stats(ff), &
      thermo(ff),green(ff),devel(ff),msd_data(ff),met(ff),pois(ff),bond(ff),angle(ff),dihedral(ff),inversion(ff), &
      tether(ff),threebody(ff),zdensity(ff),neigh(ff),vdws(ff),tersoffs(ff),fourbody(ff),rdf(ff),mpoles(ff), & 
      ext_field(ff),rigid(ff),electro(ff),domain(ff),config(ff),ewld(ff),kim_data(ff),files(ff,:),flow(ff),comm)
    End Do

    Call info('',.true.)
    Call info("*** pre-scanning stage (set_bounds) DONE ***",.true.)
    Call time_elapsed(tmr)

    Do ff=1,flow(1)%TYPE_SIZE_FF
      ! ALLOCATE SITE & CONFIG
      Call sites(ff)%init(sites(ff)%mxtmls,sites(ff)%mxatyp)
      Call config(ff)%init()
      Call neigh(ff)%init_list(config(ff)%mxatdm)
  
      ! ALLOCATE DPD ARRAYS
      Call thermo(ff)%init_dpd(vdws(ff)%max_vdw)
  
      ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS
      Call core_shells(ff)%init(config(ff)%mxatdm,sites(ff)%mxtmls,config(ff)%mxlshp,domain(ff)%neighbours)
      Call cons(ff)%init(sites(ff)%mxtmls,config(ff)%mxatdm,config(ff)%mxlshp,domain(ff)%neighbours)
      Call pmfs(ff)%init(sites(ff)%mxtmls,config(ff)%mxatdm)
      Call rigid(ff)%init(config(ff)%mxlshp,sites(ff)%mxtmls,config(ff)%mxatdm,domain(ff)%neighbours)
      Call tether(ff)%init(sites(ff)%mxtmls,config(ff)%mxatdm)
      Call bond(ff)%init(config(ff)%mxatdm,sites(ff)%mxtmls)
      Call angle(ff)%init(config(ff)%mxatdm,sites(ff)%mxtmls)
      Call dihedral(ff)%init(config(ff)%mxatdm,sites(ff)%mxtmls)
      Call inversion(ff)%init(config(ff)%mxatms,sites(ff)%mxtmls)
      Call mpoles(ff)%init(sites(ff)%max_site,neigh(ff)%max_exclude,config(ff)%mxatdm,ewld(ff)%bspline,config(ff)%mxatms)
  
      ! ALLOCATE INTER-LIKE INTERACTION ARRAYS
      Call vdws(ff)%init()
      Call met(ff)%init(config(ff)%mxatms,sites(ff)%mxatyp)
      Call tersoffs(ff)%init(sites(ff)%max_site)
      Call threebody(ff)%init(sites(ff)%max_site)
      Call fourbody(ff)%init(sites(ff)%max_site)
      Call ext_field(ff)%init()
  
      ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS
      Call rdf(ff)%init()
      Call zdensity(ff)%init(rdf(ff)%max_grid,sites(ff)%mxatyp)
      Call stats(ff)%init(config(ff)%mxatms)
      Call green(ff)%init(config(ff)%mxatms,sites(ff)%mxatyp)
  
      ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS
      Call allocate_ttm_arrays(ttms(ff),domain(ff),config(ff),comm)
      Call ttm_table_scan(config(ff)%mxbuff,ttms(ff),comm)
  
      ! Setup KIM
      Call kim_setup(kim_data(ff),config(ff)%mxatms,config(ff)%mxatdm,config(ff)%megatm,& 
        neigh(ff)%max_list,domain(ff)%mxbfxp,comm%mxnode)
 
      ! READ SIMULATION CONTROL PARAMETERS
      Call read_control(lfce,impa,ttms(ff),dfcts,rigid(ff),rsdsc(ff),core_shells(ff),cons(ff),pmfs(ff), &
        stats(ff),thermo(ff),green(ff),devel(ff),plume(ff),msd_data(ff),met(ff),pois(ff),bond(ff),angle(ff),dihedral(ff), &
        inversion(ff),zdensity(ff),neigh(ff),vdws(ff),rdf(ff),minim(ff),mpoles(ff),electro(ff),ewld(ff), &
        seed(ff),traj(ff),files(ff,:),tmr,config(ff),flow(ff),comm)
  
      ! READ SIMULATION FORCE FIELD
      write(message,'(i0)') ff
      Call info("*** DETAILS OF INTERACTIONS FOR FIELD "//trim(message)//" ***",.true.)
      Call read_field(neigh(ff)%cutoff,core_shells(ff),pmfs(ff),cons(ff),thermo(ff),met(ff),bond(ff),angle(ff), &
        dihedral(ff),inversion(ff),tether(ff),threebody(ff),sites(ff),vdws(ff),tersoffs(ff),fourbody(ff),rdf(ff), &
        mpoles(ff),ext_field(ff),rigid(ff),electro(ff),config(ff),kim_data(ff),files(ff,:),flow(ff),comm)
  
      ! If computing rdf errors, we need to initialise the arrays.
      If(rdf(ff)%l_errors_jack .or. rdf(ff)%l_errors_block) then
        Call rdf(ff)%init_block(flow(ff)%run_steps,sites(ff)%ntype_atom)
      End If
  
      ! CHECK MD CONFIGURATION
      Call check_config(config(ff),electro(ff)%key,thermo(ff),sites(ff),flow(ff),comm)
  
    End Do

    Call info('',.true.)
    Call info("*** all reading and connectivity checks DONE ***",.true.)
    Call time_elapsed(tmr)

#ifdef CHRONO
    Call stop_timer(tmr,'Initialisation')
#endif

    ! devel%l_org: translate CONFIG into CFGORG and exit gracefully
    If (devel(1)%l_org) Then
      Call info('',.true.)
      Call info("*** Translating the MD system along a vector (CONFIG to CFGORG) ***",.true.)
  
      Do ff=1,flow(1)%TYPE_SIZE_FF
        Call origin_config(config(ff),ios,devel(ff),netcdf(ff),comm)
      End Do

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If
  
    ! devel%l_scl: rescale CONFIG to CFGSCL and exit gracefully
    If (devel(1)%l_scl) Then
      Call info('',.true.)
      Call info("*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***",.true.)
  
      Do ff=1,flow(1)%TYPE_SIZE_FF
        Call scale_config(config(1),ios,devel(1),netcdf(1),comm)
      End Do

      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    ! devel%l_his: generate HISTORY and exit gracefully
    If (devel(1)%l_his) Then
      Call info('',.true.)
      Call info("*** Generating a zero timestep HISTORY frame of the MD system ***",.true.)

      Do ff=1,flow(1)%TYPE_SIZE_FF 
        Call traj(ff)%init(key=0,freq=1,start=0)
        flow(ff)%step  = 0                            ! no steps done
        flow(ff)%time  = 0.0_wp                       ! time is not relevant
        Call trajectory_write(flow(ff)%restart_key,flow(ff)%step,thermo(ff)%tstep,flow(ff)%time,ios, &
                              stats(ff)%rsd,netcdf(ff),config(ff),traj(ff),files(ff,:),comm)
      End Do 
      Call info("*** ALL DONE ***",.true.)
      Call time_elapsed(tmr)
    End If

    Do ff=1,flow(1)%TYPE_SIZE_FF
    ! Expand current system if opted for
      If (config(ff)%l_exp) Then
        Call system_expand(flow(ff)%strict,neigh(ff)%cutoff,ios,core_shells(ff), &
          cons(ff),bond(ff),angle(ff),dihedral(ff),inversion(ff),sites(ff),  &
          netcdf(ff),rigid(ff),config(ff),files(ff,:),comm)
      End If

      write(message,'(i0)') ff
      Call info("*** LONG RANGE INFORMATION FOR FIELD "//trim(message)//" ***",.true.)
    ! READ REVOLD (thermodynamic and structural data from restart file)
      Call system_init(neigh(ff)%cutoff,flow(ff)%restart_key,flow(ff)%time,flow(ff)%start_time,flow(ff)%step, &
        stats(ff),devel(ff),green(ff),thermo(ff),met(ff),bond(ff),angle(ff),dihedral(ff),inversion(ff), &
        zdensity(ff),sites(ff),vdws(ff),rdf(ff),config(ff),files(ff,:),comm)

    ! SET domain borders and link-config%cells as default for new jobs
    ! exchange atomic data and positions in border regions
      Call set_halo_particles(electro(ff)%key,neigh(ff),sites(ff),mpoles(ff),domain(ff),config(ff), &
        ewld(ff),kim_data(ff),comm)
    End Do

    Call info('',.true.)
    Call info("*** initialisation and haloing DONE ***",.true.)
    Call time_elapsed(tmr)

    ! For any intra-like interaction, construct book keeping arrays and
    ! exclusion arrays for overlapped two-body inter-like interactions
    Do ff=1,flow(1)%TYPE_SIZE_FF
      write(message,'(i0)') ff
      Call info("*** TOPOLOGY FOR FIELD "//trim(message)//" ***",.true.) 
      If (flow(ff)%book) Then
        Call build_book_intra(flow(ff)%strict,flow(ff)%print_topology,flow(ff)%simulation, & 
        flow(ff),core_shells(ff),cons(ff),pmfs(ff),bond(ff),angle(ff),dihedral(ff),inversion(ff),tether(ff), &
        neigh(ff),sites(ff),rigid(ff),domain(ff),config(ff),comm)
        If (mpoles(ff)%max_mpoles > 0) Then
          Call build_tplg_intra(neigh(ff)%max_exclude,bond(ff),angle(ff),dihedral(ff),inversion(ff), &
            mpoles(ff),config(ff),comm)
        ! multipoles topology for internal coordinate system
          If (mpoles(ff)%key == POLARISATION_CHARMM) Then
            Call build_chrm_intra(neigh(ff)%max_exclude,core_shells(ff),cons(ff),bond(ff),angle(ff), &
              dihedral(ff),inversion(ff),mpoles(ff),rigid(ff),config(ff),comm)
          End If
        ! CHARMM core-shell screened electrostatic induction interactions
        End If
        If (flow(ff)%exclusions) Then
          Call build_excl_intra(electro(ff)%lecx,core_shells(ff),cons(ff),bond(ff),angle(ff),dihedral(ff), &
            inversion(ff),neigh(ff),rigid(ff),config(ff),comm)
        End If

      Else
        Call report_topology(config(ff)%megatm,config(ff)%megfrz,config(ff)%atmfre,config(ff)%atmfrz, &
          core_shells(ff),cons(ff),pmfs(ff),bond(ff),angle(ff),dihedral(ff),inversion(ff),tether(ff),sites(ff),rigid(ff))

        ! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed
        If (flow(ff)%simulation) Then
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

    Do ff=1,flow(1)%TYPE_SIZE_FF
      write(message,'(i0)') ff
      Call info("*** DETAILS OF NEIGHBOUR LIST FOR FIELD "//trim(message)//" ***",.true.) 
      ! set and halo rotational matrices and their infinitesimal rotations
      If (mpoles(ff)%max_mpoles > 0) Then
        Call mpoles_rotmat_set_halo(mpoles(ff),domain(ff),config(ff),comm)
      End If
  
      ! SET initial system temperature
      Call set_temperature               &
        (flow(ff)%restart_key,flow(ff)%step,flow(ff)%run_steps, &
        stats(ff)%engrot,sites(ff)%dof_site,core_shells(ff),stats(ff),cons(ff),pmfs(ff),thermo(ff),minim(ff), &
        rigid(ff),domain(ff),config(ff),seed(ff),comm)
    End Do

    Call info('',.true.)
    Call info("*** temperature setting DONE ***",.true.)
    Call time_elapsed(tmr)

    ! Read ttm table file and initialise electronic temperature
    ! grid from any available restart file
    Do ff=1,flow(1)%TYPE_SIZE_FF
      If (ttms(ff)%l_ttm) Then
        Call ttm_table_read(ttms(ff),comm)
        Call ttm_system_init(flow(ff)%step,flow(ff)%equil_steps,flow(ff)%restart_key,'DUMP_E',flow(ff)%time,& 
          thermo(ff)%temp,domain(ff),ttms(ff),comm)
      End If
  
      ! Frozen atoms option
      Call freeze_atoms(config(ff))
  
      ! Cap forces in equilibration mode
      If (flow(ff)%step <= flow(ff)%equil_steps .and. flow(ff)%force_cap) Call cap_forces(thermo(ff)%temp,config(ff),comm)
  
      ! PLUMED initialisation or information message
      If (plume(ff)%l_plumed) Call plumed_init(config(ff)%megatm,thermo(ff)%tstep,thermo(ff)%temp,plume(ff),comm)
  
    End Do

    ! Print out sample of initial configuration on node zero
    Call print_initial_configuration(config(1))

    ! Indicate nodes mapped on vacuum (no particles)
    vacuum = 0
    If (config(1)%natms == 0) Then
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
    If (devel(1)%l_fast) Call gsync(comm,devel(1)%l_fast)
    ! Note to reviewer: For the time being we only pass one field. Once we are all happy with the structuring
    ! we proceed to change w_md_vv --> w_md_vv_evb, etc
    If (flow(1)%simulation) Then
      Call w_md_vv(config(1),ttms(1),ios,rsdsc(1),flow(1),core_shells(1),cons(1),pmfs(1),stats(1),thermo(1), &
        plume(1),pois(1),bond(1),angle(1),dihedral(1),inversion(1),zdensity(1),neigh(1),sites(1),fourbody(1),rdf(1), &
        netcdf(1),mpoles(1),ext_field(1),rigid(1),domain(1),seed(1),traj(1),kim_data(1),files(1,:),tmr,minim(1), &
        impa,green(1),ewld(1),electro(1),dfcts,msd_data(1),tersoffs(1),tether(1),threebody(1),vdws(1), &
        devel(1),met(1),comm)
    Else
      If (lfce) Then
        Call w_replay_historf(config(1),ios,rsdsc(1),flow(1),core_shells(1),cons(1),pmfs(1),stats(1), &
          thermo(1),plume(1),msd_data(1),bond(1),angle(1),dihedral(1),inversion(1),zdensity(1),neigh(1), &
          sites(1),vdws(1),tersoffs(1),fourbody(1),rdf(1),netcdf(1),minim(1),mpoles(1),ext_field(1),rigid(1), &
          electro(1),domain(1),seed(1),traj(1),kim_data(1),files(1,:),dfcts,tmr,tether(1),threebody(1), &
          pois(1),green(1),ewld(1),devel(1),met(1),comm)
      Else
        Call w_replay_history(config(1),ios,rsdsc(1),flow(1),core_shells(1),cons(1),pmfs(1),stats(1), &
          thermo(1),msd_data(1),met(1),pois(1),bond(1),angle(1),dihedral(1),inversion(1),zdensity(1),neigh(1), &
          sites(1),vdws(1),rdf(1),netcdf(1),minim(1),mpoles(1),ext_field(1),rigid(1),electro(1),domain(1), &
          seed(1),traj(1),kim_data(1),dfcts,files(1,:),tmr,tether(1),green(1),ewld(1),devel(1),comm)
      End If
    End If

#ifdef CHRONO
    call stop_timer(tmr,'Main Calc')
    call start_timer(tmr,'Termination')
#endif

    !Close the statis file if we used it.
    If (stats(1)%statis_file_open) Call files(1,FILE_STATS)%close()

    ! Report termination of the MD simulation
    Write(message,'(3(a,f12.3),a)') 'run terminating... elapsed  cpu time: ', &
      tmr%elapsed , ' sec, job time: ', tmr%job, ' sec, close time: ', tmr%clear_screen, ' sec'
    Call info(message,.true.)

    ! Print out sample of final configuration on node zero
    Call print_final_configuration(config(1))

    ! Two-temperature model simulations: calculate final ionic temperatures and
    !print statistics to files (final)
    If (ttms(1)%l_ttm) Then
      Call ttm_ion_temperature (ttms(1),thermo(1),domain(1),config(1),comm)
      Call printElecLatticeStatsToFile('PEAK_E', flow(1)%time, thermo(1)%temp, flow(1)%step, ttms(1)%ttmstats,ttms(1),comm)
      Call peakProfilerElec('LATS_E', flow(1)%step, ttms(1)%ttmtraj,ttms(1),comm)
      Call printLatticeStatsToFile(ttms(1)%tempion, 'PEAK_I', flow(1)%time, flow(1)%step, ttms(1)%ttmstats,ttms(1),comm)
      Call peakProfiler(ttms(1)%tempion, 'LATS_I', flow(1)%step, ttms(1)%ttmtraj,ttms(1),comm)
    End If

    ! Save restart data for real simulations only (final)
    If (flow(1)%simulation .and. (.not.devel(1)%l_tor)) Then
      Call system_revive(neigh(1)%cutoff,flow(1)%step,flow(1)%time,sites(1),ios,flow(1)%start_time,stats(1), &
        devel(1),green(1),thermo(1),bond(1),angle(1),dihedral(1),inversion(1),zdensity(1),rdf(1),netcdf(1),config(1), &
        files(1,:),comm)
      If (ttms(1)%l_ttm)Then
        Call ttm_system_revive ('DUMP_E',flow(1)%step,flow(1)%time,1,flow(1)%run_steps,ttms(1),comm)
      End If  
    End If

    ! Produce summary of simulation
    If (neigh(1)%unconditional_update .and. flow(1)%step > 0) Then
      If (.not.neigh(1)%update) Then ! Include the final skip in skipping statistics
        stats(1)%neighskip(3)=stats(1)%neighskip(2)*stats(1)%neighskip(3)
        stats(1)%neighskip(2)=stats(1)%neighskip(2)+1.0_wp
        stats(1)%neighskip(3)=stats(1)%neighskip(3)/stats(1)%neighskip(2)+stats(1)%neighskip(1)/stats(1)%neighskip(2)
        stats(1)%neighskip(4)=Min(stats(1)%neighskip(1),stats(1)%neighskip(4))
        stats(1)%neighskip(5)=Max(stats(1)%neighskip(1),stats(1)%neighskip(5))
      End If
    End If

    Call statistics_result                                        &
      (config(1),minim(1)%minimise,msd_data(1)%l_msd, &
      flow(1)%run_steps,core_shells(1)%keyshl,cons(1)%megcon,pmfs(1)%megpmf,              &
      flow(1)%step,flow(1)%time,flow(1)%start_time,config(1)%mxatdm,neigh(1)%unconditional_update,&
      stats(1),thermo(1),sites(1),comm)

    ! Final anlysis
    Call analysis_result(neigh(1)%cutoff,thermo(1), &
      bond(1),angle(1),dihedral(1),inversion(1),stats(1),green(1),zdensity(1),sites(1),rdf(1),config(1),comm)

    ! PLUMED finalisation
    If (plume(1)%l_plumed) Call plumed_finalize()

#ifdef CHRONO
    Call stop_timer(tmr,'Termination')
    Call timer_report(tmr, comm)
#endif

    ! Ask for reference in publications

    Call print_citations(electro(1),mpoles(1),ttms(1))

    ! Get just the one number to compare against

    If (devel(1)%l_eng) Then
      Write(message,'(a,1p,e20.10)') "TOTAL ENERGY: ", stats(1)%stpval(1)
      Call info('',.true.)
      Call info(message,.true.)
    End If

    ! Close output channel

    If (.not.devel(1)%l_scr) Call files(1,FILE_OUTPUT)%close()

  End Subroutine evb_molecular_dynamics_driver

  !> Allocate all types uniformly, _i.e._ N of every type
  Subroutine evb_allocate_types_uniform(array_size,thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion, &
      tether,threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs, &
      fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
      seed,traj,kim_data,config,ios,ttms,rsdsc,files)

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

    Allocate(thermo(array_size))
    Allocate(ewld(array_size))
    Allocate(tmr(1))
    Allocate(devel(array_size))
    Allocate(stats(array_size))
    Allocate(green(array_size))
    Allocate(plume(array_size))
    Allocate(msd_data(array_size))
    Allocate(met(array_size))
    Allocate(pois(array_size))
    Allocate(impa(1))
    Allocate(dfcts(1,2))
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
    Allocate(seed(array_size))
    Allocate(traj(array_size))
    Allocate(kim_data(array_size))
    Allocate(config(array_size))
    Allocate(ios(array_size))
    Allocate(ttms(array_size))
    Allocate(rsdsc(array_size))
    Allocate(files(array_size,FILENAME_SIZE))
  End Subroutine evb_allocate_types_uniform

  Subroutine evb_deallocate_types_uniform(thermo,ewld,tmr,devel,stats, &
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

  End Subroutine evb_deallocate_types_uniform

End Module evb
