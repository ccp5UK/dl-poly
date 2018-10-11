program dl_poly



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 is an stfc/ccp5 program package for the dynamical
  ! simulation of molecular systems.
  !
  ! dl_poly_4 is the property of the stfc daresbury laboratory,
  ! daresbury, warrington wa4 4ad.  no part of the package may
  ! be redistributed to third parties without the consent of
  ! daresbury laboratory.
  !
  ! dl_poly_4 is available free of charge to academic institutions
  ! engaged in non-commercial research only.  potential users not
  ! in this category must consult the ccp5 program librarian at
  ! daresbury to negotiate terms of use.
  !
  ! neither the stfc, daresbury laboratory, ccp5 nor the authors
  ! of this package claim that it is free from errors and do not
  ! accept liability for any loss or damage that may arise from
  ! its use.  it is the users responsibility to verify that the
  ! package dl_poly_4 is fit for the purpose the user intends for
  ! it.
  !
  ! users of this package are recommended to consult the dl_poly_4
  ! user and reference manuals for the full terms and conditions
  ! of its use.
  !
  ! dl_poly_4 is based on dl_poly_3 by i.t.todorov & w.smith.
  !
  ! copyright - daresbury laboratory
  ! authors   - i.t.todorov & w.smith march 2016
  ! contrib   - i.j.bush, h.a.boateng, a.m.elena, a.b.g.chalk
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Use, Intrinsic :: iso_fortran_env, Only : error_unit

  ! SETUP MODULES

  Use kinds, Only : wp,li,wi
  Use comms, Only : comms_type, init_comms, exit_comms, gsync, gtime,gmax,gsum, gbcast 
  Use constants, Only : DLP_RELEASE,DLP_VERSION

  ! PARSE MODULE

  Use parse, Only :   

  ! DEVELOPMENT MODULE

  Use development, Only : development_type,scan_development,build_info

  ! IO & DOMAINS MODULES

  Use netcdf_wrap, Only : netcdf_param
  Use domains, Only : domains_type

  ! SITE & CONFIG MODULES

  Use site, Only : site_type
  Use configuration, Only : configuration_type,check_config, scale_config, origin_config, freeze_atoms,allocate_config_arrays
  Use control, Only : read_control,scan_control_output,scan_control_io

  ! VNL module

  Use neighbours, Only : neighbours_type,vnl_check

  ! DPD module

  Use dpd, Only : dpd_thermostat

  ! INTERACTION MODULES

  Use core_shell, Only : core_shell_type,core_shell_relax,SHELL_RELAXED,SHELL_ADIABATIC,&
    core_shell_kinetic,core_shell_quench,core_shell_on_top,&
    core_shell_forces

  Use pmf, only : pmf_type,pmf_quench

  Use rigid_bodies, Only : rigid_bodies_type,rigid_bodies_quench,rigid_bodies_str_ss, &
    rigid_bodies_str__s,xscale,rigid_bodies_tags, &
    rigid_bodies_coms

  Use tethers, Only : tethers_type, tethers_forces

  Use bonds, Only : bonds_type,allocate_bonds_arrays,bonds_forces
  Use angles, Only : angles_type,allocate_angles_arrays,angles_forces
  Use dihedrals, Only : dihedrals_type,allocate_dihedrals_arrays,dihedrals_forces
  Use inversions, Only : inversions_type,allocate_inversions_arrays,inversions_forces
  Use three_body, Only : threebody_type, allocate_three_body_arrays, three_body_forces

  Use mpole, Only : mpole_type,POLARISATION_CHARMM

  Use vdw, Only : vdw_type
  Use metal, Only : metal_type
  Use tersoff, Only : tersoff_type,tersoff_forces
  Use four_body, Only : four_body_type,four_body_forces

  Use kim, Only : kim_type,kim_setup
  Use plumed, Only : plumed_type,plumed_init,plumed_finalize,plumed_apply

  Use external_field, Only : external_field_type, external_field_apply, &
    external_field_correct, &
    FIELD_NULL, FIELD_ELECTRIC, FIELD_SHEAR_OSCILLATING, &
    FIELD_SHEAR_CONTINUOUS, FIELD_GRAVITATIONAL, &
    FIELD_MAGNETIC, FIELD_SPHERE, FIELD_WALL, &
    FIELD_WALL_PISTON, FIELD_ZRES, FIELD_ZRES_MINUS, &
    FIELD_ZRES_PLUS, FIELD_ELECTRIC_OSCILLATING, &
    FIELD_UMBRELLA

  ! STATISTICS MODULES

  Use rdfs, Only : rdf_type
  Use z_density, Only : z_density_type,allocate_z_density_arrays
  Use statistics, Only : stats_type,allocate_statistics_arrays,&
    statistics_result,statistics_collect,deallocate_statistics_connect, &
    allocate_statistics_connect,statistics_connect_set, &
    statistics_connect_frames
  Use greenkubo, Only : greenkubo_type,allocate_greenkubo_arrays,vaf_compute, &
    vaf_collect,vaf_write

  ! MSD MODULE

  Use msd, Only : msd_type,msd_write

  ! TWO-TEMPERATURE MODEL MODULES

  Use drivers, Only : w_md_vv, w_replay_historf,w_replay_history
  Use errors_warnings, Only : init_error_system,info, warning

  Use minimise, Only : minimise_type,minimise_relax,zero_k_optimise
  Use two_body, Only : two_body_forces
  Use ewald, Only : ewald_type

  ! IMPACT MODULE
  Use impacts, Only : impact_type

  ! DEFECTS MODULE
  Use defects, Only : defects_type

  Use halo, Only : refresh_halo_positions,set_halo_particles
  Use deport_data, Only : mpoles_rotmat_set_halo,relocate_particles
  Use temperature, Only : scale_temperature,regauss_temperature,set_temperature
  Use rsds, Only : rsd_write,rsd_type  
  Use trajectory, Only : trajectory_write,read_history, trajectory_type
  Use system, Only : system_revive,system_expand,system_init
  Use build_excl, Only : build_excl_intra 
  Use build_book, Only : build_book_intra
  Use ffield, Only : read_field,report_topology
  Use bounds, Only : set_bounds
  Use build_tplg, Only : build_tplg_intra
  Use build_chrm, Only : build_chrm_intra
  Use thermostat, Only : thermostat_type, &
    ENS_NVE, ENS_NVT_EVANS, ENS_NVT_LANGEVIN,  &
    ENS_NVT_ANDERSON, ENS_NVT_BERENDSEN, ENS_NVT_NOSE_HOOVER, &
    ENS_NVT_GENTLE, ENS_NVT_LANGEVIN_INHOMO, &
    ENS_NPT_LANGEVIN, ENS_NPT_BERENDSEN, ENS_NPT_NOSE_HOOVER, &
    ENS_NPT_MTK, ENS_NPT_LANGEVIN_ANISO, ENS_NPT_BERENDSEN_ANISO, &
    ENS_NPT_NOSE_HOOVER_ANISO,ENS_NPT_MTK_ANISO
  Use nvt_anderson, Only : nvt_a0_vv, nvt_a1_vv
  Use nvt_berendsen, Only : nvt_b0_vv, nvt_b1_vv, nvt_b0_scl, nvt_b1_scl
  Use nvt_ekin, Only : nvt_e0_vv, nvt_e1_vv, nvt_e0_scl, nvt_e1_scl
  Use nvt_gst, Only : nvt_g0_vv, nvt_g1_vv, nvt_g0_scl, nvt_g1_scl
  Use nvt_langevin, Only : nvt_l0_vv, nvt_l1_vv, nvt_l2_vv
  Use nvt_nose_hoover, Only : nvt_h0_vv, nvt_h1_vv, nvt_h0_scl, nvt_h1_scl
  Use nst_berendsen, Only : nst_b0_vv,nst_b1_vv
  Use nst_langevin, Only : nst_l0_vv,nst_l1_vv
  Use nst_mtk, Only : nst_m0_vv,nst_m1_vv
  Use nst_nose_hoover, Only : nst_h0_vv,nst_h1_vv, nst_h0_scl, nst_h1_scl
  Use npt_berendsen, Only : npt_b0_vv,npt_b1_vv
  Use npt_langevin, Only : npt_l0_vv,npt_l1_vv
  Use npt_mtk, Only : npt_m0_vv,npt_m1_vv
  Use npt_nose_hoover, Only : npt_h0_vv,npt_h1_vv, npt_h0_scl, npt_h1_scl
  Use nve, Only : nve_0_vv, nve_1_vv 
  Use timer, Only  : timer_type, time_elapsed,timer_report,start_timer,stop_timer
  Use poisson, Only : poisson_type
  Use analysis, Only : analysis_result
  Use constraints, Only : constraints_type, constraints_quench
  Use shared_units, Only : update_shared_units, SHARED_UNIT_UPDATE_FORCES
  Use electrostatic, Only : electrostatic_type,ELECTROSTATIC_EWALD,ELECTROSTATIC_NULL
  Use stochastic_boundary, Only : stochastic_boundary_vv
  Use numerics, Only : seed_type
  Use io, Only : io_type
  Use ttm, Only : ttm_type, ttm_system_init,ttm_system_revive,ttm_table_scan,&
    ttm_table_read,allocate_ttm_arrays
  Use ttm_utils, Only : printElecLatticeStatsToFile,printLatticeStatsToFile,&
    peakProfilerElec,peakProfiler
  Use ttm_track, Only : ttm_ion_temperature,ttm_thermal_diffusion
  Use filename, Only : file_type,default_filenames,FILE_CONTROL,FILE_OUTPUT,FILE_STATS
  Use flow_control, Only : flow_type
  Use kinetics, Only : cap_forces
  
  Implicit None

! all your simulation variables
  Type(comms_type), Allocatable :: dlp_world(:),comm
  Type(thermostat_type) :: thermo
  Type(ewald_type) :: ewld
  Type(timer_type) :: tmr
  Type(development_type) :: devel
  Type(stats_type) :: stats
  Type(greenkubo_type) :: green
  Type(plumed_type) :: plume
  Type(msd_type) :: msd_data
  Type(metal_type) :: met
  Type(poisson_type) :: pois
  Type(impact_type) :: impa
  Type(defects_type) :: dfcts(2)
  Type(bonds_type) :: bond
  Type( angles_type ) :: angle
  Type( dihedrals_type ) :: dihedral
  Type( inversions_type ) :: inversion
  Type( tethers_type ) :: tether
  Type( threebody_type ) :: threebody
  Type( z_density_type ) :: zdensity
  Type( constraints_type ) :: cons
  Type( neighbours_type ) :: neigh
  Type( pmf_type ) :: pmfs
  Type( site_type ) :: sites
  Type( core_shell_type ) :: core_shells
  Type( vdw_type ) :: vdws
  Type( tersoff_type ) :: tersoffs
  Type( four_body_type ) :: fourbody
  Type( rdf_type ) :: rdf
  Type( netcdf_param ) :: netcdf
  Type( minimise_type ) :: minim
  Type( mpole_type ) :: mpoles
  Type( external_field_type ) :: ext_field
  Type( rigid_bodies_type ) :: rigid
  Type( electrostatic_type ) :: electro
  Type( domains_type ) :: domain
  Type( flow_type ) :: flow
  Type( seed_type ) :: seed
  Type( trajectory_type ) :: traj
  Type( kim_type ), Target :: kim_data
  Type( configuration_type ) :: config
  Type( io_type) :: ios
  Type( ttm_type) :: ttms
  Type( rsd_type ), Target :: rsdsc
  Type( file_type ), Allocatable :: files(:)
  
  ! Local Variables

  Integer(Kind=wi) :: i,j
  Logical :: lfce
  Character( Len = 256 ) :: message,messages(5)
  Character( Len = 66 )  :: banner(13)
  Character( Len = 1024 ) :: control_filename

  Character( Len = * ), Parameter :: fmt1 = '(a)', &
                                     fmt2 = '(a25,a8,a4,a14,a15)', &
                                     fmt3 = '(a,i10,a)'

  ! SET UP COMMUNICATIONS & CLOCKING

  Allocate(dlp_world(0:0))
  Call init_comms(dlp_world(0))
  !dlp_world(0)%ou=nrite
  !Call init_error_system(nrite,dlp_world(0))
  comm=dlp_world(0) ! this shall vanish asap w_ are proper things
  Call gsync(dlp_world(0))
  Call gtime(tmr%elapsed) ! Initialise wall clock time
  If (dlp_world(0)%idnode == 0) Then
    If (command_argument_count() == 1 ) Then
      Call get_command_argument(1, control_filename)
    End If
  End If

  ! Set default file names
  Call default_filenames(files)
  ! Rename control file if argument was passed
  If (command_argument_count() == 1) Then
    Call files(FILE_CONTROL)%rename(control_filename)
  End If

  Call scan_development(devel,files,comm)

  ! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES

  Call scan_control_output(files,comm)

  If (dlp_world(0)%idnode == 0) Then
    ! Open output file, or direct output unit to stderr
    If (.not.devel%l_scr) Then
      Open(Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Status='replace')
    Else
      files(FILE_OUTPUT)%unit_no = error_unit
    End If
  End If
  Call gbcast(dlp_world(0),files(FILE_OUTPUT)%unit_no,0)
  dlp_world(0)%ou=files(FILE_OUTPUT)%unit_no
  Call init_error_system(files(FILE_OUTPUT)%unit_no,dlp_world(0))

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
  Write(banner(13),fmt1) "******************************************************************"
  Call info(banner,13,.true.)

  Call build_info(devel)

  Call info('',.true.)
  Write(banner(1),fmt1) Repeat("*",66)
  Write(banner(2),fmt1) "****  Please do cite `J. Mater. Chem.', 16, 1911-1918 (2006)  ****"
  Write(banner(3),fmt1) "****  when publishing research data obtained using DL_POLY_4  ****"
  Write(banner(4),fmt1) Repeat("*",66)
  Call info(banner,4,.true.)

  ! TEST I/O

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
  Call time_elapsed(tmr%elapsed)

  ! ALLOCATE SITE & CONFIG

  Call sites%init(sites%mxtmls,sites%mxatyp)
  Call allocate_config_arrays(config)
  Call neigh%init_list(config%mxatdm)

  ! ALLOCATE DPD ARRAYS

  Call thermo%init_dpd(vdws%max_vdw)

  ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS

  Call core_shells%init(config%mxatdm,sites%mxtmls,config%mxlshp,domain%neighbours)

  Call cons%init(sites%mxtmls,config%mxatdm,config%mxlshp,domain%neighbours)
  Call pmfs%init(sites%mxtmls,config%mxatdm)

  Call rigid%init(config%mxlshp,sites%mxtmls,config%mxatdm,domain%neighbours)

  Call tether%init(sites%mxtmls,config%mxatdm)

  Call allocate_bonds_arrays(bond,config%mxatdm,sites%mxtmls)
  Call allocate_angles_arrays(angle,config%mxatdm,sites%mxtmls)
  Call allocate_dihedrals_arrays(dihedral,config%mxatdm,sites%mxtmls)
  Call allocate_inversions_arrays(inversion,config%mxatms,sites%mxtmls)

  Call mpoles%init(sites%max_site,neigh%max_exclude,config%mxatdm,ewld%bspline,config%mxatms)

  ! ALLOCATE INTER-LIKE INTERACTION ARRAYS

  Call vdws%init()
  Call met%init(config%mxatms,sites%mxatyp)
  Call tersoffs%init(sites%max_site)
  Call allocate_three_body_arrays(sites%max_site,threebody)
  Call fourbody%init(sites%max_site)

  Call ext_field%init()

  ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS

  Call rdf%init()
  Call allocate_z_density_arrays(zdensity,rdf%max_grid,sites%mxatyp)
  Call allocate_statistics_arrays(config%mxatms,stats)
  Call allocate_greenkubo_arrays(green,config%mxatms,sites%mxatyp)

  ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS

  Call allocate_ttm_arrays(ttms,domain,config,comm)
  Call ttm_table_scan(config%mxbuff,ttms,comm)

  ! Setup KIM
  Call kim_setup(kim_data,config%mxatms,config%mxatdm,config%megatm,neigh%max_list,domain%mxbfxp,comm%mxnode)

  ! READ SIMULATION CONTROL PARAMETERS

  Call read_control(lfce,impa,ttms,dfcts,rigid,rsdsc,core_shells,cons,pmfs, &
    stats,thermo,green,devel,plume,msd_data,met,pois,bond,angle,dihedral, &
    inversion,zdensity,neigh,vdws,tersoffs,rdf, minim,mpoles,electro,ewld, &
    seed,traj,files,tmr,config,flow,comm)

  ! READ SIMULATION FORCE FIELD

  Call read_field(neigh%cutoff,core_shells,pmfs,cons,thermo,met,bond,angle, &
    dihedral,inversion,tether,threebody,sites,vdws,tersoffs,fourbody,rdf, &
    mpoles,ext_field,rigid,electro,config,kim_data,files,flow,comm)

  ! If computing rdf errors, we need to initialise the arrays.
  If(rdf%l_errors_jack .or. rdf%l_errors_block) then
    Call rdf%init_block(flow%run_steps,sites%ntype_atom)
  End If

  ! CHECK MD CONFIGURATION

  Call check_config(config,electro%key,thermo,sites,flow,comm)

  Call info('',.true.)
  Call info("*** all reading and connectivity checks DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! devel%l_org: translate CONFIG into CFGORG and exit gracefully

  If (devel%l_org) Then
    Call info('',.true.)
    Call info("*** Translating the MD system along a vector (CONFIG to CFGORG) ***",.true.)

    Call origin_config(config,ios,devel,netcdf,comm)

    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! devel%l_scl: rescale CONFIG to CFGSCL and exit gracefully

  If (devel%l_scl) Then
    Call info('',.true.)
    Call info("*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***",.true.)

    Call scale_config(config,ios,devel,netcdf,comm)

    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! devel%l_his: generate HISTORY and exit gracefully

  If (devel%l_his) Then
    Call info('',.true.)
    Call info("*** Generating a zero timestep HISTORY frame of the MD system ***",.true.)

    Call traj%init(key=0,freq=1,start=0)
    flow%step  = 0                            ! no steps done
    flow%time   = 0.0_wp                       ! time is not relevant
    Call trajectory_write(flow%restart_key,flow%step,thermo%tstep,flow%time,ios,stats%rsd,netcdf,config,traj,files,comm)

    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! Expand current system if opted for

  If (config%l_exp) Then
    Call system_expand(flow%strict,neigh%cutoff,ios,core_shells, &
      cons,bond,angle,dihedral,inversion,sites,netcdf,rigid,config,files,comm)
  End If

  ! EXIT gracefully

  If (devel%l_trm) Then
    Call info('',.true.)
    Call info("*** Exiting gracefully ***",.true.)
    Go To 10
  End If

  ! READ REVOLD (thermodynamic and structural data from restart file)

  Call system_init(neigh%cutoff,flow%restart_key,flow%time,flow%start_time,flow%step, &
    core_shells,stats,devel,green,thermo,met,bond,angle,dihedral,inversion, &
    zdensity,sites,vdws,rdf,config,files,comm)

  ! SET domain borders and link-config%cells as default for new jobs
  ! exchange atomic data and positions in border regions

  Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,config,ewld,kim_data,comm)

  Call info('',.true.)
  Call info("*** initialisation and haloing DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! For any intra-like interaction, construct book keeping arrays and
  ! exclusion arrays for overlapped two-body inter-like interactions

  If (flow%book) Then
    Call build_book_intra(flow%strict,flow%print_topology,flow%simulation, &
      flow,core_shells,cons,pmfs,bond,angle,dihedral,inversion,tether, &
      neigh,sites,mpoles,rigid,domain,config,comm)
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
      pmfs,bond,angle,dihedral,inversion,tether,sites,rigid,comm)

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
  Call time_elapsed(tmr%elapsed)

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
  Call time_elapsed(tmr%elapsed)

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

  Call info('',.true.)
  Call info('sample of starting configuration on node zero:',.true.)
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

  j=(config%natms+19)/20
  If (j > 0) Then
    Do i=1,config%natms,j
      If (config%levcfg <= 1) Then
        Write(message,'(i8,1p,6e12.4)') &
          config%ltg(i),config%parts(i)%xxx,config%parts(i)%yyy,config%parts(i)%zzz,config%vxx(i),config%vyy(i),config%vzz(i)
      End If

      If (config%levcfg == 2) Then
        Write(message,"(i8,1p,9e12.4)") &
        config%ltg(i),config%parts(i)%xxx,config%parts(i)%yyy,config%parts(i)%zzz,&
        config%vxx(i),config%vyy(i),config%vzz(i),config%parts(i)%fxx,config%parts(i)%fyy,&
        config%parts(i)%fzz
      End If
      Call info(message,.true.)
    End Do
  End If
  Call info('',.true.)

  ! Indicate nodes mapped on vacuum (no particles)

  j=0
  If (config%natms == 0) Then
    j=1
    Call warning('mapped on vacuum (no particles)')
  End If
  Call gsum(comm,j)
  If (j > 0) Call warning(2,Real(j,wp),Real(comm%mxnode,wp),0.0_wp)

  ! start-up time when forces are not recalculated

  Call time_elapsed(tmr%elapsed)

  ! Now you can run fast, boy

  If (devel%l_fast) Call gsync(comm,devel%l_fast)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  If (flow%simulation) Then
    Call w_md_vv(config,ttms,ios,rsdsc,flow,core_shells,cons,pmfs,stats,thermo, &
      plume,pois,bond,angle,dihedral,inversion,zdensity,neigh,sites,fourbody,rdf, &
      netcdf,mpoles,ext_field,rigid,domain,seed,traj,kim_data,files,tmr,minim, &
      impa,green,ewld,electro,dfcts,msd_data,tersoffs,tether,threebody,vdws, &
      devel,met,comm)
  Else
    If (lfce) Then
      Call w_replay_historf(config,ios,rsdsc,flow,core_shells,cons,pmfs,stats, &
        thermo,plume,msd_data,bond,angle,dihedral,inversion,zdensity,neigh, &
        sites,vdws,tersoffs,fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid, &
        electro,domain,seed,traj,kim_data,files,dfcts,tmr,tether,threebody, &
        pois,green,ewld,devel,met,comm)
    Else
      Call w_replay_history(config,ios,rsdsc,flow,core_shells,cons,pmfs,stats, &
        thermo,msd_data,met,pois,bond,angle,dihedral,inversion,zdensity,neigh, &
        sites,vdws,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain, &
        seed,traj,kim_data,dfcts,files,tmr,tether,green,ewld,devel,threebody,comm)
    End If
  End If

  !Close the statis file if we used it.
  If (stats%statis_file_open) Close(Unit=files(FILE_STATS)%unit_no)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Report termination of the MD simulation

  Write(message,'(3(a,f12.3),a)') 'run terminating... elapsed  cpu time: ', &
    tmr%elapsed , ' sec, job time: ', tmr%job, ' sec, close time: ', tmr%clear_screen, ' sec'
  Call info(message,.true.)

  ! Print out sample of final configuration on node zero

  Call info('',.true.)
  Call info("sample of final configuration on node zero",.true.)
  Write(message,'(7x,a1,7x,a4,2(8x,a4),6(7x,a5))') &
    'i', 'x(i)', 'y(i)', 'z(i)', 'vx(i)', 'vy(i)', 'vz(i)', &
    'fx(i)', 'fy(i)', 'fz(i)'
  Call info(message,.true.)
  j=(config%natms+19)/20
  If (j > 0) Then
    Do i=1,config%natms,j
      If (config%levcfg <= 1) Then
        Write(message,'(i8,1p,6e12.4)') &
          config%ltg(i),config%parts(i)%xxx,config%parts(i)%yyy,config%parts(i)%zzz,config%vxx(i),config%vyy(i),config%vzz(i)
      End If

      If (config%levcfg == 2) Then
        Write(message,"(i8,1p,9e12.4)") &
        config%ltg(i),config%parts(i)%xxx,config%parts(i)%yyy,config%parts(i)%zzz,config%vxx(i),config%vyy(i),&
        config%vzz(i),config%parts(i)%fxx,config%parts(i)%fyy,config%parts(i)%fzz
      End If
      Call info(message,.true.)
    End Do
  End If
  Call info('',.true.)

  ! Two-temperature model simulations: calculate final 
  ! ionic temperatures and print statistics to files
  ! (final)

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
    stats,thermo,green,sites,comm)

  ! Final anlysis
  Call analysis_result(flow%step,neigh%cutoff,thermo, &
    bond,angle,dihedral,inversion,stats,green,zdensity,neigh,sites,rdf,config,comm)

  10 Continue

  ! PLUMED finalisation

  If (plume%l_plumed) Call plumed_finalize()

#ifdef CHRONO
  Call timer_report(tmr,comm)
#endif
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
  Call info(Repeat("*",66),.true.)

  ! Get just the one number to compare against

  If (devel%l_eng) Then
    Write(message,'(a,1p,e20.10)') "TOTAL ENERGY: ", stats%stpval(1)
    Call info('',.true.)
    Call info(message,.true.)
  End If

  ! Close output channel

  If (dlp_world(0)%idnode == 0 .and. (.not.devel%l_scr)) Close(unit=files(FILE_OUTPUT)%unit_no)

  ! Terminate job

  Call gsync(dlp_world(0))
  Call exit_comms(comm)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate(dlp_world)
End Program dl_poly
