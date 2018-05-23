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
  Use domains

  ! SITE & CONFIG MODULES

  Use site
  Use configuration
  Use kontrol, Only : read_control

  ! VNL module

  Use vnl

  ! DPD module

  Use dpd

  ! INTERACTION MODULES

  Use core_shell

  Use constraints
  Use pmf

  Use rigid_bodies

  Use tethers

  Use bonds
  Use angles
  Use dihedrals
  Use inversions

  Use mpole

  Use vdw
  Use metal
  Use tersoff
  Use three_body
  Use four_body

  Use kim
  Use plumed

  Use external_field

  ! STATISTICS MODULES

  Use rdfs
  Use z_density
  Use statistics, Only : stats_type,allocate_statistics_arrays,&
    statistics_result,statistics_collect,deallocate_statistics_connect, &
    allocate_statistics_connect,statistics_connect_set, &
    statistics_connect_frames
  Use greenkubo, Only : greenkubo_type,allocate_greenkubo_arrays,vaf_compute, &
                        vaf_collect,vaf_write

  ! MSD MODULE

  Use msd

  ! KINETIC MODULE

  Use kinetics

  ! LANGEVIN MODULE

  Use langevin

  ! TWO-TEMPERATURE MODEL MODULES

  Use ttm
  Use ttm_utils
  
  Use drivers
  Use errors_warnings, Only : init_error_system
  
  Use minimise, Only : passmin, minimise_relax,zero_k_optimise
  Use two_body, Only : two_body_forces
  Use ewald, Only : ewald_type
  Use impacts, Only : impact_type
  Use halo, Only : refresh_halo_positions,set_halo_particles
  Use deport_data, Only : mpoles_rotmat_set_halo,relocate_particles
  Use temperature, Only : scale_temperature,regauss_temperature,set_temperature
  Use rsds, Only : rsd_write
  Use defects, Only : defects_write
  Use trajectory, Only : trajectory_write,read_history
  use system, Only : system_revive,system_expand,system_init
  Use ttm_track, Only : ttm_ion_temperature,ttm_thermal_diffusion
  Use build_excl, Only : build_excl_intra 
  Use build_book, Only : build_book_intra
  Use ffield, Only : read_field,report_topology
  Use kontrol, Only : scan_control_output,scan_control_io
  Use bounds, Only : set_bounds
  Use build_tplg, Only : build_tplg_intra
  use build_chrm, Only : build_chrm_intra

  Use thermostat, Only : thermostat_type
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
    ! MAIN PROGRAM VARIABLES
  Use timer, Only  : timer_type, time_elapsed,timer_report
  Implicit None

  ! newjob used for trajectory_write &
  !                 defects_write    &
  !                 msd_write        &
  !                 rsd_write        &

  Logical, Save :: newjob = .true.

  ! OUTPUT existence

  Logical               :: l_out
  Character( Len = 10 ) :: c_out

  ! lines and page used for printing controls

  Integer       :: lines = 0 , &
    npage = 8

  ! general flags

  Logical           :: ltmp,l_vv,l_n_e,l_n_v,       &
    l_ind,l_str,l_top,           &
    l_exp,lecx,lfcap,      &
    lmin,          &
    lvar,leql,lsim,lfce,    &
    lpana,lrdf,lprdf,lzdn,lpzdn, &
    lvafav,lpvaf,                &
    ltraj,ldef,lrsd,             &
    safe,lbook,lexcl,            &
    relaxed_shl = .true.,        &
    relaxed_min = .true.

  Integer           :: i,j,isw,levcfg,nstfce,              &
    nx,ny,nz,                           &
    keyres,nstrun,nsteql,               &
    keymin,nstmin,                      &
    keyens,nstbpo,    &
    keyfce,mxshak,mxquat,               &
    nstbnd,nstang,nstdih,nstinv,        &
    nstrdf,nstzdn,                      &
    nstmsd,istmsd,nstraj,istraj,keytrj, &
    nsdef,isdef,nsrsd,isrsd,            &
    ndump,nstep,keyshl,                 &
    atmfre,atmfrz,megatm,megfrz,        &
    megshl,megcon,megpmf,megrgd,        &
    megtet,megbnd,megang,megdih,meginv

  ! Degrees of freedom must be in long integers so we do 2.1x10^9 particles

  Integer(Kind=li)  :: degfre,degshl,degtra,degrot

  ! elrc,virlrc - vdw energy and virial are scalars and in vdw
  ! elrcm,vlrcm - metal energy and virial are array-like and in metal_module

  Real( Kind = wp ) :: tsths,                                     &
    tstep,time,tmst,      &
    dvar,rcut,rpad,rlnk,                       &
    rvdw,rmet,rbin,rcter,rctbp,rcfbp,          &
    alpha,epsq,fmax,                           &
    width,mndis,mxdis,mxstp,     &
    rlx_tol(1:2),min_tol(1:2),                 &
    tolnce,quattol,rdef,rrsd,                  &
    pdplnc

  Type(comms_type), Allocatable :: dlp_world(:),comm
  Type(thermostat_type) :: thermo
  Type(ewald_type) :: ewld
  Type(timer_type) :: tmr
  Type(impact_type) :: impa
  Type(development_type) :: devel
  Type(stats_type) :: stats
  Type(greenkubo_type) :: green

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
  If (dlp_world(0)%mxnode > 1) Call gsync(dlp_world(0))
  Call gtime(tmr%elapsed)
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

  Call scan_control_io(comm)

  ! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
  ! (setup and domains)

  Call set_bounds                                     &
    (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,l_ind, &
    dvar,rcut,rpad,rlnk,rvdw,rmet,rbin,nstfce,alpha,width,stats,thermo,green,devel,comm)

  Call gtime(tmr%elapsed)
  Call info('',.true.)
  Call info("*** pre-scanning stage (set_bounds) DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! ALLOCATE SITE & CONFIG

  Call allocate_site_arrays()
  Call allocate_config_arrays()

  ! ALLOCATE DPD ARRAYS

  Call allocate_dpd_arrays(thermo)

  ! ALLOCATE INTRA-LIKE INTERACTION ARRAYS

  Call allocate_core_shell_arrays()

  Call allocate_constraints_arrays()
  Call allocate_pmf_arrays()

  Call allocate_rigid_bodies_arrays()

  Call allocate_tethers_arrays()

  Call allocate_bonds_arrays()
  Call allocate_angles_arrays()
  Call allocate_dihedrals_arrays()
  Call allocate_inversions_arrays()

  Call allocate_mpoles_arrays()

  ! ALLOCATE INTER-LIKE INTERACTION ARRAYS

  Call allocate_vdw_arrays()
  Call allocate_metal_arrays()
  Call allocate_tersoff_arrays()
  Call allocate_three_body_arrays()
  Call allocate_four_body_arrays()

  Call allocate_external_field_arrays()

  ! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS

  Call allocate_rdf_arrays()
  Call allocate_z_density_arrays()
  Call allocate_statistics_arrays(mxatdm,stats)
  Call allocate_greenkubo_arrays(green)

  ! ALLOCATE TWO-TEMPERATURE MODEL ARRAYS

  Call allocate_ttm_arrays(comm)
  Call ttm_table_scan(comm)

  ! READ SIMULATION CONTROL PARAMETERS

  Call read_control                                    &
    (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,        &
    rcut,rpad,rvdw,rbin,nstfce,alpha,width,     &
    l_exp,lecx,lfcap,l_top,lmin,          &
    lvar,leql,               &
    lfce,lpana,lrdf,lprdf,lzdn,lpzdn,           &
    lvafav,lpvaf,ltraj,ldef,lrsd,               &
    nx,ny,nz,impa,                            &
    keyres,                   &
    tstep,mndis,mxdis,mxstp,nstrun,nsteql,      &
    keymin,nstmin,min_tol,                      &
    keyens,&
    fmax,nstbpo,keyfce,epsq,             &
    rlx_tol,mxshak,tolnce,mxquat,quattol,       &
    nstbnd,nstang,nstdih,nstinv,nstrdf,nstzdn,  &
    nstmsd,istmsd,nstraj,istraj,keytrj,         &
    nsdef,isdef,rdef,nsrsd,isrsd,rrsd,          &
    ndump,pdplnc,stats,thermo,green,devel,tmr,comm)

  ! READ SIMULATION FORCE FIELD

  Call read_field                          &
    (l_str,l_top,l_n_v,             &
    rcut,rvdw,rmet,width,epsq, &
    keyens,keyfce,keyshl,           &
    lecx,lbook,lexcl,               &
    rcter,rctbp,rcfbp,              &
    atmfre,atmfrz,megatm,megfrz,    &
    megshl,megcon,megpmf,megrgd,    &
    megtet,megbnd,megang,megdih,meginv,thermo,comm)

  ! If computing rdf errors, we need to initialise the arrays.
  If(l_errors_jack .or. l_errors_block) then
    Call allocate_block_average_array(nstrun)
  End If

  ! If using induced dipoles then read in atomic polarizability

  !  If (induce) Call read_polarity()

  ! CHECK MD CONFIGURATION

  Call check_config(levcfg,l_str,keyens,keyfce,keyres,megatm,thermo,comm)

  Call gtime(tmr%elapsed)
  Call info('',.true.)
  Call info("*** all reading and connectivity checks DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! devel%l_org: translate CONFIG into CFGORG and exit gracefully

  If (devel%l_org) Then
    Call gtime(tmr%elapsed)
    Call info('',.true.)
    Call info("*** Translating the MD system along a vector (CONFIG to CFGORG) ***",.true.)

    Call origin_config(megatm,devel,comm)

    Call gtime(tmr%elapsed)
    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! devel%l_scl: rescale CONFIG to CFGSCL and exit gracefully

  If (devel%l_scl) Then
    Call gtime(tmr%elapsed)
    Call info('',.true.)
    Call info("*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***",.true.)

    Call scale_config(megatm,devel,comm)

    Call gtime(tmr%elapsed)
    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! devel%l_his: generate HISTORY and exit gracefully

  If (devel%l_his) Then
    Call gtime(tmr%elapsed)
    Call info('',.true.)
    Call info("*** Generating a zero timestep HISTORY frame of the MD system ***",.true.)

    ! Nail down necessary parameters

    nstraj = 0 ; istraj = 1 ; keytrj = 0  ! default trajectory
    nstep  = 0                            ! no steps done
    time   = 0.0_wp                       ! time is not relevant
    Call trajectory_write(keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,stats%rsd,comm)

    Call gtime(tmr%elapsed)
    Call info("*** ALL DONE ***",.true.)
    Call time_elapsed(tmr%elapsed)
  End If

  ! Expand current system if opted for

  If (l_exp) Call system_expand(l_str,rcut,nx,ny,nz,megatm,comm)

  ! EXIT gracefully

  If (devel%l_trm) Then
    Call info('',.true.)
    Call info("*** Exiting gracefully ***",.true.)
    Go To 10
  End If

  ! READ REVOLD (thermodynamic and structural data from restart file)

  Call system_init                                                 &
    (levcfg,rcut,rvdw,rbin,rmet,lrdf,lzdn,keyres,megatm,    &
    time,tmst,nstep,tstep,elrc,virlrc,elrcm,vlrcm,stats,devel,green,thermo,comm)

  ! SET domain borders and link-cells as default for new jobs
  ! exchange atomic data and positions in border regions

  Call set_halo_particles(rlnk,keyfce,comm)

  Call gtime(tmr%elapsed)
  Call info('',.true.)
  Call info("*** initialisation and haloing DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! For any intra-like interaction, construct book keeping arrays and
  ! exclusion arrays for overlapped two-body inter-like interactions

  If (lbook) Then
    Call build_book_intra              &
      (l_str,l_top,lsim,dvar,      &
      megatm,megfrz,atmfre,atmfrz, &
      megshl,megcon,megpmf,        &
      megrgd,degrot,degtra,        &
      megtet,megbnd,megang,megdih,meginv,comm)
    If (mximpl > 0) Then
      Call build_tplg_intra(comm) ! multipoles topology for internal coordinate system
      If (keyind == 1) Call build_chrm_intra(comm) ! CHARMM core-shell screened electrostatic induction interactions
    End If
    If (lexcl) Call build_excl_intra(lecx,comm)
  Else
    Call report_topology                &
      (megatm,megfrz,atmfre,atmfrz, &
      megshl,megcon,megpmf,megrgd,  &
      megtet,megbnd,megang,megdih,meginv,comm)

    ! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed

    If (lsim) Then
      Call deallocate_core_shell_arrays()

      Call deallocate_constraints_arrays()
      Call deallocate_pmf_arrays()

      Call deallocate_rigid_bodies_arrays()

      Call deallocate_tethers_arrays()

      Call deallocate_bonds_arrays()
      Call deallocate_angles_arrays()
      Call deallocate_dihedrals_arrays()
      Call deallocate_inversions_arrays()
    End If
  End If

  Call gtime(tmr%elapsed)
  Call info('',.true.)
  Call info("*** bookkeeping DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! set and halo rotational matrices and their infinitesimal rotations

  If (mximpl > 0) Call mpoles_rotmat_set_halo(comm)

  ! SET initial system temperature

  Call set_temperature               &
    (levcfg,keyres,      &
    lmin,nstep,nstrun,nstmin, &
    mxshak,tolnce,keyshl,     &
    atmfre,atmfrz,            &
    megshl,megcon,megpmf,     &
    megrgd,degtra,degrot,     &
    degfre,degshl,stats%engrot,thermo,comm)

  Call gtime(tmr%elapsed)
  Call info('',.true.)
  Call info("*** temperature setting DONE ***",.true.)
  Call time_elapsed(tmr%elapsed)

  ! Read ttm table file and initialise electronic temperature
  ! grid from any available restart file

  If (l_ttm) Then
    Call ttm_table_read(comm)
    Call ttm_system_init(nstep,nsteql,keyres,'DUMP_E',time,thermo%temp,comm)
  End If

  ! Frozen atoms option

  Call freeze_atoms()

  ! Cap forces in equilibration mode

  If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,thermo%temp,comm)

  ! PLUMED initialisation or information message

  If (l_plumed) Call plumed_init(megatm,tstep,thermo%temp,comm)

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
          ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)
      End If

      If (levcfg == 2) Then
        Write(message,"(i8,1p,9e12.4)") &
        ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i),fxx(i),fyy(i),fzz(i)
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

  Call gtime(tmr%elapsed)
  Call time_elapsed(tmr%elapsed)

  ! Now you can run fast, boy

  If (devel%l_fast) Call gsync(comm,devel%l_fast)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  If (lsim) Then
    Call w_md_vv(mxatdm,stats,thermo)
  Else
    If (lfce) Then
      Call w_replay_historf(mxatdm,stats,thermo)
    Else
      Call w_replay_history(mxatdm,stats,thermo)
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
          ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)
      End If

      If (levcfg == 2) Then
        Write(message,"(i8,1p,9e12.4)") &
        ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i),fxx(i),fyy(i),fzz(i)
      End If
      Call info(message,.true.)
    End Do
  End If
  Call info('',.true.)

  ! Two-temperature model simulations: calculate final 
  ! ionic temperatures and print statistics to files
  ! (final)

  If (l_ttm) Then
    Call ttm_ion_temperature (thermo,comm)
    Call printElecLatticeStatsToFile('PEAK_E', time, thermo%temp, nstep, ttmstats,comm)
    Call peakProfilerElec('LATS_E', nstep, ttmtraj,comm)
    Call printLatticeStatsToFile(tempion, 'PEAK_I', time, nstep, ttmstats,comm)
    Call peakProfiler(tempion, 'LATS_I', nstep, ttmtraj,comm)
  End If

  ! Save restart data for real simulations only (final)

  If (lsim .and. (.not.devel%l_tor)) Then
    Call system_revive &
      (rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
      stats,devel,green,thermo,comm)
    If (l_ttm) Call ttm_system_revive ('DUMP_E',nstep,time,1,nstrun,comm)
  End If

  ! Produce summary of simulation

  Call statistics_result                                        &
    (rcut,lmin,lpana,lrdf,lprdf,lzdn,lpzdn,lvafav,lpvaf, &
    nstrun,keyens,keyshl,megcon,megpmf,              &
    nstep,tstep,time,tmst,mxatdm,stats,thermo,green,comm,passmin)

  10 Continue

  ! PLUMED finalisation

  If (l_plumed) Call plumed_finalize()

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
  If (keyfce == 2) Then
    Write(banner(1),fmt1) '****   - I.J. Bush, I.T. Todorov & W. Smith,                  ****'
    Write(banner(2),fmt1) '****     Comp. Phys. Commun., 175, 323-329 (2006),            ****'
    Write(banner(3),fmt1) '****     https://doi.org/10.1016/j.cpc.2006.05.001            ****'
    Call info(banner,3,.true.)
  End If
  If (mximpl > 0) Then
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

  If (dlp_world(0)%mxnode > 1) Call gsync(dlp_world(0))
  Call exit_comms(comm)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate(dlp_world)
Contains

  Subroutine w_calculate_forces(stat)
    Type(stats_type), Intent(InOut) :: stat
    Include 'w_calculate_forces.F90'
  End Subroutine w_calculate_forces

  Subroutine w_refresh_mappings(stat)
    Type(stats_type), Intent(InOut) :: stat

    Include 'w_refresh_mappings.F90'
  End Subroutine w_refresh_mappings

  Subroutine w_integrate_vv(isw,stat,thermo)
    Integer, Intent( In    ) :: isw ! used for vv stage control
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo

    Include 'w_integrate_vv.F90'
  End Subroutine w_integrate_vv

  Subroutine w_kinetic_options(stat)
    Type(stats_type), Intent(InOut) :: stat
    Include 'w_kinetic_options.F90'
  End Subroutine w_kinetic_options

  Subroutine w_statistics_report(mxatdm_,stat)
    Integer( Kind = wi ), Intent ( In ) :: mxatdm_
    Type(stats_type), Intent(InOut) :: stat
    Include 'w_statistics_report.F90'
  End Subroutine w_statistics_report

  Subroutine w_write_options(stat)
    Type(stats_type), Intent(InOut) :: stat
    Include 'w_write_options.F90'
  End Subroutine w_write_options

  Subroutine w_refresh_output()
    Include 'w_refresh_output.F90'
  End Subroutine w_refresh_output

  Subroutine w_md_vv(mxatdm_,stat,thermo)
    Integer( Kind = wi ), Intent ( In ) :: mxatdm_
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo
    Include 'w_md_vv.F90'
  End Subroutine w_md_vv

  Subroutine w_replay_history(mxatdm_,stat,thermo)
    Integer( Kind = wi ), Intent( In  )  :: mxatdm_ 
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo

    Logical,     Save :: newjb = .true.
    Real( Kind = wp ) :: tmsh        ! tmst replacement
    Integer( Kind = wi )           :: nstpe,nstph ! nstep replacements
    Integer           :: exout       ! exit indicator for reading

    Include 'w_replay_history.F90'
  End Subroutine w_replay_history

  Subroutine w_replay_historf(mxatdm_,stat,thermo)
    Integer( Kind = wi ), Intent( In  )  :: mxatdm_ 
    Type(stats_type), Intent(InOut) :: stat
    Type(thermostat_type), Intent(InOut) :: thermo

    Logical,     Save :: newjb = .true.
    Real( Kind = wp ) :: tmsh        ! tmst replacement
    Integer           :: nstpe,nstph ! nstep replacements
    Integer           :: exout       ! exit indicator for reading

    Include 'w_replay_historf.F90'
  End Subroutine w_replay_historf


End Program dl_poly
