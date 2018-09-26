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



  ! SETUP MODULES

  Use kinds, Only : wp,li,wi
  Use comms, Only : comms_type, init_comms, exit_comms, gsync, gtime 
  Use setup

  ! PARSE MODULE

  Use parse 

  ! DEVELOPMENT MODULE

  Use development, Only : development_type,scan_development,build_info

  ! IO & DOMAINS MODULES

  Use io
  Use netcdf_wrap, Only : netcdf_param
  Use domains, Only : domains_type

  ! SITE & CONFIG MODULES

  Use site, Only : site_type
  Use configuration
  Use kontrol, Only : read_control

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

  ! KINETIC MODULE

  Use kinetics

  ! LANGEVIN MODULE

  Use langevin

  ! TWO-TEMPERATURE MODEL MODULES

  Use ttm
  Use ttm_utils
  
  Use drivers
  Use errors_warnings, Only : init_error_system
  
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
  Use defects, Only : defects_write
  Use trajectory, Only : trajectory_write,read_history
  use system, Only : system_revive,system_expand,system_init
  Use ttm_track, Only : ttm_ion_temperature,ttm_thermal_diffusion
  Use build_excl, Only : build_excl_intra 
  Use build_book, Only : build_book_intra
  Use ffield, Only : read_field,report_topology
  Use kontrol, Only : scan_control_output,scan_control_io, control_type
  Use bounds, Only : set_bounds
  Use build_tplg, Only : build_tplg_intra
  use build_chrm, Only : build_chrm_intra
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
  Use trajectory, Only : trajectory_type
    ! MAIN PROGRAM VARIABLES
  Implicit None

  ! newjob used for trajectory_write &
  !                 defects_write    &
  !                 msd_write        &
  !                 rsd_write        &

  Logical, Save :: newjob = .true.

  ! general flags

  Logical           :: ltmp,l_vv,l_n_e,l_n_v,       &
    l_ind,l_str,l_top,           &
    l_exp,lecx,lfcap,      &
    lvar,leql,lsim,lfce,    &
    lpana, &
    ltraj,             &
    safe,lbook,lexcl,            &
    relaxed_shl = .true.

  Integer           :: i,j,isw,levcfg,nstfce,              &
    nx,ny,nz,                           &
    keyres,nstrun,nsteql,               &
    nstbpo,    &
    mxquat,               &
    nstbnd,nstang,nstdih,nstinv,        &
    ndump,nstep,                 &
    atmfre,atmfrz,megatm,megfrz

  ! Degrees of freedom must be in long integers so we do 2.1x10^9 particles

  Integer(Kind=li)  :: degfre,degshl,degtra,degrot

  ! vdws%elrc,vdws%vlrc - vdw energy and virial are scalars and in vdw

  Real( Kind = wp ) :: tsths,                                     &
    tstep,time,tmst,      &
    dvar,                       &
    rbin,          &
    fmax,                           &
    width,mndis,mxdis,mxstp,     &
    rlx_tol(1:2),                 &
    quattol,                 &
    pdplnc

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
  Type( control_type ) :: flow
  Type( seed_type ) :: seed
  Type( trajectory_type ) :: traj
  Type( kim_type ), Target :: kim_data
  Type( rsd_type ), Target :: rsdsc

  Character( Len = 256 ) :: message,messages(5)
  Character( Len = 66 )  :: banner(13)

  Character( Len = * ), Parameter :: fmt1 = '(a)', &
                                     fmt2 = '(a25,a8,a4,a14,a15)', &
                                     fmt3 = '(a,i10,a)'

  ! SET UP COMMUNICATIONS & CLOCKING

  Allocate(dlp_world(0:0))
  Call init_comms(dlp_world(0))
  dlp_world(0)%ou=nrite
  Call init_error_system(nrite,dlp_world(0))
  comm=dlp_world(0) ! this shall vanish asap w_ are proper things
  Call gsync(dlp_world(0))
  Call gtime(tmr%elapsed) ! Initialise wall clock time
  If (dlp_world(0)%idnode == 0) Then
    If (command_argument_count() == 1 ) Then
      Call get_command_argument(1, control)
    End If
  End If

  Call scan_development(devel,comm)

  ! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES

  Call scan_control_output(comm)

  If (dlp_world(0)%idnode == 0) Then
    If (.not.devel%l_scr) Then
      Open(Unit=nrite, File=Trim(output), Status='replace')
    End If
  End If

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

  Call scan_control_io(netcdf,comm)

  ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
  ! (setup and domains)

  Call set_bounds (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,l_ind, &
    dvar,rbin,nstfce,width,sites%max_site,core_shells,cons,pmfs,stats, &
    thermo,green,devel,msd_data,met,pois,bond,angle,dihedral,inversion, &
    tether,threebody,zdensity,neigh,vdws,tersoffs,fourbody,rdf,mpoles,ext_field, &
    rigid,electro,domain,ewld,kim_data,comm)

  Call info('',.true.)
  Call info("*** pre-scanning stage (set_bounds) DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! ALLOCATE SITE & CONFIG

  Call sites%init(mxtmls,mxatyp)
  Call allocate_config_arrays()
  Call neigh%init_list(mxatdm)

  ! ALLOCATE DPD ARRAYS

  Call thermo%init_dpd(vdws%max_vdw)

  ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS

  Call core_shells%init(mxatdm,mxtmls,mxlshp,domain%neighbours)

  Call cons%init(mxtmls,mxatdm,mxlshp,domain%neighbours)
  Call pmfs%init(mxtmls,mxatdm)

  Call rigid%init(mxtmls,mxatms,domain%neighbours)

  Call tether%init(mxtmls,mxatdm)

  Call allocate_bonds_arrays(bond)
  Call allocate_angles_arrays(angle)
  Call allocate_dihedrals_arrays(dihedral)
  Call allocate_inversions_arrays(inversion)

  Call mpoles%init(sites%max_site,neigh%max_exclude,mxatdm,ewld%bspline,mxatms)

  ! ALLOCATE INTER-LIKE INTERACTION ARRAYS

  Call vdws%init()
  Call met%init(mxatms,mxatyp)
  Call tersoffs%init(sites%max_site)
  Call allocate_three_body_arrays(sites%max_site,threebody)
  Call fourbody%init(sites%max_site)

  Call ext_field%init()

  ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS

  Call rdf%init()
  Call allocate_z_density_arrays(zdensity,rdf%max_grid,mxatyp)
  Call allocate_statistics_arrays(mxatdm,stats)
  Call allocate_greenkubo_arrays(green)

  ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS

  Call allocate_ttm_arrays(domain,comm)
  Call ttm_table_scan(comm)

  ! Setup KIM
  Call kim_setup(kim_data,mxatms,mxatdm,megatm,neigh%max_list,mxbfxp,comm%mxnode)

  ! READ SIMULATION CONTROL PARAMETERS

  Call read_control                                    &
    (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,        &
    rbin,nstfce,width,     &
    l_exp,lecx,lfcap,l_top,          &
    lvar,leql,               &
    lfce,lpana,           &
    ltraj,               &
    nx,ny,nz,impa,                            &
    keyres,                   &
    tstep,mndis,mxdis,mxstp,nstrun,nsteql,      &
    fmax,nstbpo,             &
    rlx_tol,mxquat,quattol,       &
    nstbnd,nstang,nstdih,nstinv,  &
    dfcts,          &
    ndump,pdplnc,rsdsc,core_shells,cons,pmfs,stats,thermo,green,devel,plume,msd_data, &
    met,pois,bond,angle,dihedral,inversion,zdensity,neigh,vdws,tersoffs,rdf, &
    minim,mpoles,electro,ewld,seed,traj,tmr,comm)

  ! READ SIMULATION FORCE FIELD

  Call read_field(l_str,l_top,l_n_v,neigh%cutoff,width,lecx,lbook,lexcl,atmfre, &
    atmfrz,megatm,megfrz,core_shells,pmfs,cons,thermo,met,bond,angle,dihedral, &
    inversion,tether,threebody,sites,vdws,tersoffs,fourbody,rdf,mpoles, &
    ext_field,rigid,electro,kim_data,comm)

  ! If computing rdf errors, we need to initialise the arrays.
  If(rdf%l_errors_jack .or. rdf%l_errors_block) then
    Call rdf%init_block(nstrun,sites%ntype_atom)
  End If

  ! CHECK MD CONFIGURATION

  Call check_config(levcfg,l_str,electro%key,keyres,megatm,thermo,sites,comm)

  Call info('',.true.)
  Call info("*** all reading and connectivity checks DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! devel%l_org: translate CONFIG into CFGORG and exit gracefully

  If (devel%l_org) Then
    Call info('',.true.)
    Call info("*** Translating the MD system along a vector (CONFIG to CFGORG) ***",.true.)

    Call origin_config(megatm,devel,netcdf,comm)

    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! devel%l_scl: rescale CONFIG to CFGSCL and exit gracefully

  If (devel%l_scl) Then
    Call info('',.true.)
    Call info("*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***",.true.)

    Call scale_config(megatm,devel,netcdf,comm)

    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! devel%l_his: generate HISTORY and exit gracefully

  If (devel%l_his) Then
    Call info('',.true.)
    Call info("*** Generating a zero timestep HISTORY frame of the MD system ***",.true.)

    Call traj%init(key=0,freq=1,start=0)
    nstep  = 0                            ! no steps done
    time   = 0.0_wp                       ! time is not relevant
    Call trajectory_write(keyres,megatm,nstep,tstep,time,stats%rsd,netcdf,parts,traj,comm)

    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! Expand current system if opted for

  If (l_exp) Then
    Call system_expand(l_str,neigh%cutoff,nx,ny,nz,megatm,core_shells, &
      cons,bond,angle,dihedral,inversion,sites,netcdf,rigid,parts,comm)
  End If

  ! EXIT gracefully

  If (devel%l_trm) Then
    Call info('',.true.)
    Call info("*** Exiting gracefully ***",.true.)
    Go To 10
  End If

  ! READ REVOLD (thermodynamic and structural data from restart file)

  Call system_init                                                 &
    (levcfg,neigh%cutoff,rbin,keyres,megatm,    &
    time,tmst,nstep,tstep,core_shells,stats,devel, &
    green,thermo,met,bond,angle,dihedral,inversion,zdensity,sites,vdws,rdf,parts,comm)

  ! SET domain borders and link-cells as default for new jobs
  ! exchange atomic data and positions in border regions

  Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,ewld,kim_data,comm)

  Call info('',.true.)
  Call info("*** initialisation and haloing DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! For any intra-like interaction, construct book keeping arrays and
  ! exclusion arrays for overlapped two-body inter-like interactions

  If (lbook) Then
    Call build_book_intra(l_str,l_top,lsim,dvar,megatm,megfrz,atmfre,atmfrz, &
      degrot,degtra,flow,core_shells,cons,pmfs,bond,angle,dihedral,inversion,tether, &
      neigh,sites,mpoles,rigid,domain,parts,comm)
    If (mpoles%max_mpoles > 0) Then
      Call build_tplg_intra(neigh%max_exclude,bond,angle,dihedral,inversion, &
        mpoles,comm)
      ! multipoles topology for internal coordinate system
      If (mpoles%key == POLARISATION_CHARMM) Then
        Call build_chrm_intra(neigh%max_exclude,core_shells,cons,bond,angle, &
          dihedral,inversion,mpoles,rigid,comm)
      End If
       ! CHARMM core-shell screened electrostatic induction interactions
    End If
    If (lexcl) Then
      Call build_excl_intra(lecx,core_shells,cons,bond,angle,dihedral, &
        inversion,neigh,rigid,comm)
    End If
  Else
    Call report_topology(megatm,megfrz,atmfre,atmfrz,core_shells,cons, &
      pmfs,bond,angle,dihedral,inversion,tether,sites,rigid,comm)

    ! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed

    If (lsim) Then
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
    Call mpoles_rotmat_set_halo(mpoles,domain,comm)
  End If

  ! SET initial system temperature

  Call set_temperature               &
    (levcfg,keyres,nstep,nstrun,atmfre,atmfrz,degtra,degrot,degfre,degshl, &
    stats%engrot,sites%dof_site,core_shells,stats,cons,pmfs,thermo,minim, &
    rigid,domain,parts,seed,comm)

  Call info('',.true.)
  Call info("*** temperature setting DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! Read ttm table file and initialise electronic temperature
  ! grid from any available restart file

  If (l_ttm) Then
    Call ttm_table_read(comm)
    Call ttm_system_init(nstep,nsteql,keyres,'DUMP_E',time,thermo%temp,domain,comm)
  End If

  ! Frozen atoms option

  Call freeze_atoms()

  ! Cap forces in equilibration mode

  If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,thermo%temp,parts,comm)

  ! PLUMED initialisation or information message

  If (plume%l_plumed) Call plumed_init(megatm,tstep,thermo%temp,plume,comm)

  ! Print out sample of initial configuration on node zero

  Call info('',.true.)
  Call info('sample of starting configuration on node zero:',.true.)
  If (levcfg <= 1) Then
    Write(message,'(7x,a1,7x,a4,2(8x,a4),3(7x,a5))') &
      'i', 'x(i)', 'y(i)', 'z(i)', 'vx(i)', 'vy(i)', 'vz(i)'
    Call info(message,.true.)
  End If

  If (levcfg == 2) Then
    Write(message,'(7x,a1,7x,a4,2(8x,a4),6(7x,a5))') &
      'i', 'x(i)', 'y(i)', 'z(i)', 'vx(i)', 'vy(i)', 'vz(i)', &
      'fx(i)', 'fy(i)', 'fz(i)'
    Call info(message,.true.)
  End If

  j=(natms+19)/20
  If (j > 0) Then
    Do i=1,natms,j
      If (levcfg <= 1) Then
        Write(message,'(i8,1p,6e12.4)') &
          ltg(i),parts(i)%xxx,parts(i)%yyy,parts(i)%zzz,vxx(i),vyy(i),vzz(i)
      End If

      If (levcfg == 2) Then
        Write(message,"(i8,1p,9e12.4)") &
        ltg(i),parts(i)%xxx,parts(i)%yyy,parts(i)%zzz,vxx(i),vyy(i),vzz(i),parts(i)%fxx,parts(i)%fyy,parts(i)%fzz
      End If
      Call info(message,.true.)
    End Do
  End If
  Call info('',.true.)

  ! Indicate nodes mapped on vacuum (no particles)

  j=0
  If (natms == 0) Then
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


  If (lsim) Then
    Call w_md_vv(mxatdm,rsdsc,flow,core_shells,cons,pmfs,stats,thermo,plume,&
      pois,bond,angle,dihedral,inversion,zdensity,neigh,sites,fourbody,rdf, &
      netcdf,mpoles,ext_field,rigid,domain,seed,traj,kim_data,tmr)
  Else
    If (lfce) Then
      Call w_replay_historf(mxatdm,rsdsc,flow,core_shells,cons,pmfs,stats,thermo,plume,&
        msd_data,bond,angle,dihedral,inversion,zdensity,neigh,sites,vdws,tersoffs, &
        fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain,seed,traj, &
        kim_data,tmr)
    Else
      Call w_replay_history(mxatdm,rsdsc,flow,core_shells,cons,pmfs,stats,thermo,msd_data,&
        met,pois,bond,angle,dihedral,inversion,zdensity,neigh,sites,vdws,rdf, &
        netcdf,minim,mpoles,ext_field,rigid,electro,domain,seed,traj,kim_data)
    End If
  End If

  !Close the statis file if we used it.
  If (stats%statis_file_open) Close(Unit=nstats)

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
  j=(natms+19)/20
  If (j > 0) Then
    Do i=1,natms,j
      If (levcfg <= 1) Then
        Write(message,'(i8,1p,6e12.4)') &
          ltg(i),parts(i)%xxx,parts(i)%yyy,parts(i)%zzz,vxx(i),vyy(i),vzz(i)
      End If

      If (levcfg == 2) Then
        Write(message,"(i8,1p,9e12.4)") &
        ltg(i),parts(i)%xxx,parts(i)%yyy,parts(i)%zzz,vxx(i),vyy(i),vzz(i),parts(i)%fxx,parts(i)%fyy,parts(i)%fzz
      End If
      Call info(message,.true.)
    End Do
  End If
  Call info('',.true.)

  ! Two-temperature model simulations: calculate final 
  ! ionic temperatures and print statistics to files
  ! (final)

  If (l_ttm) Then
    Call ttm_ion_temperature (thermo,domain,comm)
    Call printElecLatticeStatsToFile('PEAK_E', time, thermo%temp, nstep, ttmstats,comm)
    Call peakProfilerElec('LATS_E', nstep, ttmtraj,comm)
    Call printLatticeStatsToFile(tempion, 'PEAK_I', time, nstep, ttmstats,comm)
    Call peakProfiler(tempion, 'LATS_I', nstep, ttmtraj,comm)
  End If

  ! Save restart data for real simulations only (final)

  If (lsim .and. (.not.devel%l_tor)) Then
    Call system_revive &
      (neigh%cutoff,rbin,megatm,nstep,tstep,time,tmst, &
      stats,devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf,comm)
    If (l_ttm) Call ttm_system_revive ('DUMP_E',nstep,time,1,nstrun,comm)
  End If

  ! Produce summary of simulation
  If (neigh%unconditional_update .and. nstep > 0) Then
     If (.not.neigh%update) Then ! Include the final skip in skipping statistics
        stats%neighskip(3)=stats%neighskip(2)*stats%neighskip(3)
        stats%neighskip(2)=stats%neighskip(2)+1.0_wp
        stats%neighskip(3)=stats%neighskip(3)/stats%neighskip(2)+stats%neighskip(1)/stats%neighskip(2)
        stats%neighskip(4)=Min(stats%neighskip(1),stats%neighskip(4))
        stats%neighskip(5)=Max(stats%neighskip(1),stats%neighskip(5))
     End If
  End If

  Call statistics_result                                        &
    (minim%minimise,msd_data%l_msd, &
    nstrun,core_shells%keyshl,cons%megcon,pmfs%megpmf,              &
    nstep,tstep,time,tmst,mxatdm,neigh%unconditional_update,&
    stats,thermo,green,sites,comm)

  ! Final anlysis
  Call analysis_result(lpana, &
                       nstep,tstep,neigh%cutoff,stats%sumval(2),thermo%ensemble, &
                       bond,angle,dihedral,inversion,stats,green,zdensity,neigh,sites,rdf,comm)

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

  If (dlp_world(0)%idnode == 0 .and. (.not.devel%l_scr)) Close(Unit=nrite)

  ! Terminate job

  Call gsync(dlp_world(0))
  Call exit_comms(comm)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate(dlp_world)
Contains

  Subroutine w_calculate_forces(flw,cshell,cons,pmf,stat,plume,pois,bond,angle,dihedral,&
      inversion,tether,threebody,neigh,sites,vdws,tersoffs,fourbody,rdf,netcdf, &
      minim,mpoles,ext_field,rigid,electro,domain,kim_data,tmr)
    Type( control_type ), Intent( InOut ) :: flw
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(plumed_type), Intent(InOut) :: plume
    Type(poisson_type), Intent(InOut) :: pois
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( tethers_type ), Intent( InOut ) :: tether
    Type( threebody_type ), Intent( InOut ) :: threebody
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( tersoff_type ), Intent( InOut )  :: tersoffs
    Type( four_body_type ), Intent( InOut ) :: fourbody
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( netcdf_param ), Intent( In    ) :: netcdf
    Type( minimise_type ), Intent( InOut ) :: minim
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( external_field_type ), Intent( InOut ) :: ext_field
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( kim_type ), Intent( InOut) :: kim_data
    Type( timer_type ), Intent( InOut ) :: tmr
    Include 'w_calculate_forces.F90'
  End Subroutine w_calculate_forces

  Subroutine w_refresh_mappings(flw,cshell,cons,pmf,stat,msd_data,bond,angle, &
    dihedral,inversion,tether,neigh,sites,mpoles,rigid,domain,kim_data)
    Type( control_type ), Intent( InOut ) :: flw
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(msd_type), Intent(InOut) :: msd_data
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( tethers_type ), Intent( InOut ) :: tether
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( site_type ), Intent( InOut ) :: sites
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( kim_type ), Intent( InOut ) :: kim_data
    Include 'w_refresh_mappings.F90'
  End Subroutine w_refresh_mappings

  Subroutine w_integrate_vv(isw,cshell,cons,pmf,stat,thermo,sites,vdws,rigid,domain,seed,tmr)
    Integer, Intent( In    ) :: isw ! used for vv stage control
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( seed_type ), Intent( InOut ) :: seed
    Type( timer_type ), Intent( InOut ) :: tmr
    Include 'w_integrate_vv.F90'
  End Subroutine w_integrate_vv

  Subroutine w_kinetic_options(cshell,cons,pmf,stat,sites,ext_field,domain,seed)
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type( site_type ), Intent( InOut ) :: sites
    Type( external_field_type ), Intent( In    ) :: ext_field
    Type( domains_type ), Intent( In    ) :: domain
    Type( seed_type ), Intent( InOut ) :: seed
    Include 'w_kinetic_options.F90'
  End Subroutine w_kinetic_options

  Subroutine w_statistics_report(mxatdm_,cshell,cons,pmf,stat,msd_data,zdensity, &
      sites,rdf,domain,flw)
    Integer( Kind = wi ), Intent ( In ) :: mxatdm_
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(msd_type), Intent(InOut) :: msd_data
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( site_type ), Intent( InOut ) :: sites
    Type( rdf_type ), Intent( In    ) :: rdf
    Type( domains_type ), Intent( In    ) :: domain
    Type( control_type ), Intent( InOut ) :: flw
    Include 'w_statistics_report.F90'
  End Subroutine w_statistics_report

  Subroutine w_write_options(rsdc,cshell,stat,sites,netcdf,domain,traj)
    Type( rsd_type ), Intent( Inout ) :: rsdc
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type(stats_type), Intent(InOut) :: stat
    Type( site_type ), Intent( InOut ) :: sites
    Type( netcdf_param ), Intent( In    ) :: netcdf
    Type( domains_type ), Intent( In    ) :: domain
    Type( trajectory_type ), Intent( InOut ) :: traj
    Include 'w_write_options.F90'
  End Subroutine w_write_options

  Subroutine w_refresh_output()
    Logical :: l_out
    Character(Len=10) :: c_out
    Include 'w_refresh_output.F90'
  End Subroutine w_refresh_output

  Subroutine w_md_vv(mxatdm_,rsdc,flw,cshell,cons,pmf,stat,thermo,plume,pois,bond,angle, &
    dihedral,inversion,zdensity,neigh,sites,fourbody,rdf,netcdf,mpoles, &
    ext_field,rigid,domain,seed,traj,kim_data,tmr)
    Integer( Kind = wi ), Intent( In ) :: mxatdm_
    Type( rsd_type ), Intent( InOut ) :: rsdc
    Type( control_type ), Intent( InOut ) :: flw
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo
    Type(plumed_type), Intent(InOut) :: plume
    Type(poisson_type), Intent(InOut) :: pois
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( site_type ), Intent( InOut ) :: sites
    Type( four_body_type ), Intent( InOut ) :: fourbody
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( netcdf_param ), Intent( In    ) :: netcdf
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( external_field_type ), Intent( InOut ) :: ext_field
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( seed_type ), Intent( InOut ) :: seed
    Type( trajectory_type ), Intent( InOut ) :: traj
    Type( kim_type ), Intent( InOut ) :: kim_data
    Type( timer_type ), Intent( InOut ) :: tmr
    Include 'w_md_vv.F90'
  End Subroutine w_md_vv

  Subroutine w_replay_history(mxatdm_,rsdc,flw,cshell,cons,pmf,stat,thermo,msd_data, &
      met,pois,bond,angle,dihedral,inversion,zdensity,neigh,sites,vdws,rdf, &
      netcdf,minim,mpoles,ext_field,rigid,electro,domain,seed,traj,kim_data)
    Integer( Kind = wi ), Intent( In  )  :: mxatdm_
    Type( rsd_type ), Intent( Inout ) :: rsdc
    Type( control_type ), Intent( InOut ) :: flw
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo
    Type(msd_type), Intent(InOut) :: msd_data
    Type( metal_type ), Intent( InOut ) :: met
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( netcdf_param ), Intent( In    ) :: netcdf
    Type( minimise_type ), Intent( InOut ) :: minim
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( external_field_type ), Intent( InOut ) :: ext_field
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( seed_type ), Intent( InOut ) :: seed
    Type( trajectory_type ), Intent( InOut ) :: traj
    Type( kim_type ), Intent( InOut ) :: kim_data

    Logical,     Save :: newjb = .true.
    Real( Kind = wp ) :: tmsh        ! tmst replacement
    Integer( Kind = wi )           :: nstpe,nstph ! nstep replacements
    Integer           :: exout       ! exit indicator for reading
    Logical :: l_out
    Character(Len=10) :: c_out

    Include 'w_replay_history.F90'
  End Subroutine w_replay_history

  Subroutine w_replay_historf(mxatdm_,rsdc,flw,cshell,cons,pmf,stat,thermo,plume, &
    msd_data,bond,angle,dihedral,inversion,zdensity,neigh,sites,vdws,tersoffs, &
    fourbody,rdf,netcdf,minim,mpoles,ext_field,rigid,electro,domain,seed,traj, &
    kim_data,tmr)
    Integer( Kind = wi ), Intent( In  )  :: mxatdm_
    Type( rsd_type ), Intent( Inout ) :: rsdc
    Type( control_type ), Intent( InOut ) :: flw
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo
    Type(plumed_type), Intent(InOut) :: plume
    Type(msd_type), Intent(InOut) :: msd_data
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( tersoff_type ), Intent( InOut )  :: tersoffs
    Type( four_body_type ), Intent( InOut ) :: fourbody
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( netcdf_param ), Intent( In    ) :: netcdf
    Type( minimise_type ), Intent( InOut ) :: minim
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( external_field_type ), Intent( InOut ) :: ext_field
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( seed_type ), Intent( InOut ) :: seed
    Type( trajectory_type ), Intent( InOut ) :: traj
    Type( kim_type ), Intent( InOut ) :: kim_data
    Type( timer_type ), Intent( InOut ) :: tmr

    Logical,     Save :: newjb = .true.
    Real( Kind = wp ) :: tmsh        ! tmst replacement
    Integer           :: nstpe,nstph ! nstep replacements
    Integer           :: exout       ! exit indicator for reading
    Logical :: l_out
    Character(Len=10) :: c_out

    Include 'w_replay_historf.F90'
  End Subroutine w_replay_historf


End Program dl_poly
