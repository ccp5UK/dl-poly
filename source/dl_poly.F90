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
  ! contrib   - i.j.bush, h.a.boateng, m.a.seaton, a.m.elena,
  !             s.l.daraszewicz, g.khara, a.brukhno, a.b.g.chalk
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! SETUP MODULES

  Use comms, Only : comms_type, init_comms, exit_comms, gsync, gtime,gsum

  ! DEVELOPMENT MODULE

  Use development, Only : development_type

  ! IO & DOMAINS MODULES

  Use netcdf_wrap, Only : netcdf_param
  Use domains, Only : domains_type

  ! SITE & CONFIG MODULES

  Use site, Only : site_type
  Use configuration, Only : configuration_type,check_config, scale_config, origin_config, freeze_atoms
  Use control, Only : read_control,scan_control_output,scan_control_io

  ! VNL module

  Use neighbours, Only : neighbours_type

  ! INTERACTION MODULES

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

  Use mpole, Only : mpole_type

  Use vdw, Only : vdw_type
  Use metal, Only : metal_type
  Use tersoff, Only : tersoff_type
  Use four_body, Only : four_body_type

  Use kim, Only : kim_type,kim_setup
  Use plumed, Only : plumed_type,plumed_init,plumed_finalize

  Use external_field, Only : external_field_type

  ! STATISTICS MODULES

  Use rdfs, Only : rdf_type
  Use z_density, Only : z_density_type
  Use statistics, Only : stats_type,statistics_result
  Use greenkubo, Only : greenkubo_type

  ! MSD MODULE

  Use msd, Only : msd_type

  ! TWO-TEMPERATURE MODEL MODULES

  Use drivers, Only : w_md_vv, w_replay_historf,w_replay_history
  Use errors_warnings, Only : init_error_system,info, warning

  Use ewald, Only : ewald_type

  ! IMPACT MODULE
  Use impacts, Only : impact_type

  ! DEFECTS MODULE
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
  Use poisson, Only : poisson_type
  Use analysis, Only : analysis_result
  Use constraints, Only : constraints_type
  Use electrostatic, Only : electrostatic_type
  Use numerics, Only : seed_type
  Use io, Only : io_type
  Use ttm, Only : ttm_type, ttm_system_init,ttm_system_revive,ttm_table_scan,&
    ttm_table_read,allocate_ttm_arrays
  Use ttm_utils, Only : printElecLatticeStatsToFile,printLatticeStatsToFile,&
    peakProfilerElec,peakProfiler
  Use ttm_track, Only : ttm_ion_temperature
  Use filename, Only : file_type,default_filenames
  Use flow_control, Only : flow_type,EVB,FFS,MD
  Use kinetics, Only : cap_forces
  Use timer, Only  : timer_type
  Use meta, Only : molecular_dynamics

  Implicit None

  ! all your simulation variables
  Type(comms_type), Allocatable :: dlp_world(:)
  Type(thermostat_type), Allocatable :: thermo(:)
  Type(ewald_type), Allocatable :: ewld(:)
  Type(timer_type), Allocatable :: tmr(:)
  Type(development_type), Allocatable :: devel(:)
  Type(stats_type), Allocatable :: stats(:)
  Type(greenkubo_type), Allocatable :: green(:)
  Type(plumed_type), Allocatable :: plume(:)
  Type(msd_type), Allocatable :: msd_data(:)
  Type(metal_type), Allocatable :: met(:)
  Type(poisson_type), Allocatable :: pois(:)
  Type(impact_type), Allocatable :: impa(:)
  Type(defects_type), Allocatable :: dfcts(:,:)
  Type(bonds_type), Allocatable :: bond(:)
  Type( angles_type ), Allocatable :: angle(:)
  Type( dihedrals_type ), Allocatable :: dihedral(:)
  Type( inversions_type ), Allocatable :: inversion(:)
  Type( tethers_type ), Allocatable :: tether(:)
  Type( threebody_type ), Allocatable :: threebody(:)
  Type( z_density_type ), Allocatable :: zdensity(:)
  Type( constraints_type ), Allocatable :: cons(:)
  Type( neighbours_type ), Allocatable :: neigh(:)
  Type( pmf_type ), Allocatable :: pmfs(:)
  Type( site_type ), Allocatable :: sites(:)
  Type( core_shell_type ), Allocatable :: core_shells(:)
  Type( vdw_type ), Allocatable :: vdws(:)
  Type( tersoff_type ), Allocatable :: tersoffs(:)
  Type( four_body_type ), Allocatable :: fourbody(:)
  Type( rdf_type ), Allocatable :: rdf(:)
  Type( netcdf_param ), Allocatable :: netcdf(:)
  Type( minimise_type ), Allocatable :: minim(:)
  Type( mpole_type ), Allocatable :: mpoles(:)
  Type( external_field_type ), Allocatable :: ext_field(:)
  Type( rigid_bodies_type ), Allocatable :: rigid(:)
  Type( electrostatic_type ), Allocatable :: electro(:)
  Type( domains_type ), Allocatable :: domain(:)
  Type( flow_type ), Allocatable :: flow(:)
  Type( seed_type ), Allocatable :: seed(:)
  Type( trajectory_type ), Allocatable :: traj(:)
  Type( kim_type ), Allocatable, Target :: kim_data(:)
  Type( configuration_type ), Allocatable :: config(:)
  Type( io_type), Allocatable :: ios(:)
  Type( ttm_type), Allocatable :: ttms(:)
  Type( rsd_type ), Allocatable, Target :: rsdsc(:)
  Type( file_type ), Allocatable :: files(:,:)

  ! Local Variables
  character( len = 1024 ) :: control_filename

  ! SET UP COMMUNICATIONS & CLOCKING

  Allocate(dlp_world(0:0))
  Call init_comms(dlp_world(0))
  !dlp_world(0)%ou=nrite
  !Call init_error_system(nrite,dlp_world(0))
  Call gsync(dlp_world(0))
  If (dlp_world(0)%idnode == 0) Then
    If (command_argument_count() == 1 ) Then
      Call get_command_argument(1, control_filename)
    End If
  End If

  ! temporary stuff this will need to be abstracted 
  Allocate(flow(1))
  flow(1)%simulation_method = MD
  ! Select metasimulation method
  If (flow(1)%simulation_method == MD) Then
    write(0,*) "simulation type: MD" 
    Call molecular_dynamics(dlp_world,thermo,ewld,tmr,devel,stats, &
      green,plume,msd_data,met,pois,impa,dfcts,bond,angle,dihedral,inversion,tether, &
      threebody,zdensity,cons,neigh,pmfs,sites,core_shells,vdws,tersoffs,fourbody, &
      rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain,flow,seed,traj, &
      kim_data,config,ios,ttms,rsdsc,files,control_filename)
  Else If (flow(1)%simulation_method == EVB) Then 
    write(0,*) "simulation type: EVB" 
  Else If (flow(1)%simulation_method == FFS) Then 
    write(0,*) "simulation type: FFS" 
  Else
    write(0,*) "Unknown simulation type" 
  End IF

  ! Terminate job

  Call gsync(dlp_world(0))
  Call exit_comms(dlp_world)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate(flow)
  Deallocate(dlp_world)
End Program dl_poly
