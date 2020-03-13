Program dl_poly

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
  ! EVB       - i.scivetti march-october 2019 
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use angles,                             Only: angles_type
  Use angular_distribution,               Only: adf_type
  Use bonds,                              Only: bonds_type
  Use comms,                              Only: comms_type,&
                                                exit_comms,&
                                                gbcast,&
                                                gsync,&
                                                init_comms
  Use configuration,                      Only: configuration_type
  Use constraints,                        Only: constraints_type
  Use coord,                              Only: coord_type
  Use core_shell,                         Only: core_shell_type
  Use defects,                            Only: defects_type
  Use development,                        Only: development_type
  Use dihedrals,                          Only: dihedrals_type
  Use domains,                            Only: domains_type
  Use electrostatic,                      Only: electrostatic_type
  Use errors_warnings,                    Only: init_error_system
  Use ewald,                              Only: ewald_type
  Use external_field,                     Only: external_field_type
  Use filename,                           Only: file_type
  Use flow_control,                       Only: EmpVB,&
                                                FFS,&
                                                MD,&
                                                flow_type,&
                                                read_simtype
  Use four_body,                          Only: four_body_type
  Use greenkubo,                          Only: greenkubo_type
  Use impacts,                            Only: impact_type
  Use inversions,                         Only: inversions_type
  Use io,                                 Only: io_type
  Use, Intrinsic :: iso_fortran_env,      Only: eu => error_unit
  Use kim,                                Only: kim_type
  Use meta,                               Only: molecular_dynamics
  Use meta_evb,                           Only: evb_molecular_dynamics
  Use metal,                              Only: metal_type
  Use minimise,                           Only: minimise_type
  Use mpole,                              Only: mpole_type
  Use msd,                                Only: msd_type
  Use neighbours,                         Only: neighbours_type
  Use netcdf_wrap,                        Only: netcdf_param
  Use numerics,                           Only: seed_type
  Use plumed,                             Only: plumed_type
  Use pmf,                                Only: pmf_type
  Use poisson,                            Only: poisson_type
  Use rdfs,                               Only: rdf_type
  Use rigid_bodies,                       Only: rigid_bodies_type
  Use rsds,                               Only: rsd_type
  Use site,                               Only: site_type
  Use statistics,                         Only: stats_type
  Use tersoff,                            Only: tersoff_type
  Use tethers,                            Only: tethers_type
  Use thermostat,                         Only: thermostat_type
  Use three_body,                         Only: threebody_type
  Use timer,                              Only: timer_type
  Use trajectory,                         Only: trajectory_type
  Use ttm,                                Only: ttm_type
  Use vdw,                                Only: vdw_type
  Use z_density,                          Only: z_density_type

  Implicit None

  ! all your simulation variables
  Type(comms_type), Allocatable          :: dlp_world(:)
  Type(thermostat_type), Allocatable     :: thermo(:)
  Type(ewald_type), Allocatable          :: ewld(:)
  Type(timer_type), Allocatable          :: tmr(:)
  Type(development_type), Allocatable    :: devel(:)
  Type(stats_type), Allocatable          :: stats(:)
  Type(greenkubo_type), Allocatable      :: green(:)
  Type(plumed_type), Allocatable         :: plume(:)
  Type(msd_type), Allocatable            :: msd_data(:)
  Type(metal_type), Allocatable          :: met(:)
  Type(poisson_type), Allocatable        :: pois(:)
  Type(impact_type), Allocatable         :: impa(:)
  Type(defects_type), Allocatable        :: dfcts(:, :)
  Type(bonds_type), Allocatable          :: bond(:)
  Type(angles_type), Allocatable         :: angle(:)
  Type(dihedrals_type), Allocatable      :: dihedral(:)
  Type(inversions_type), Allocatable     :: inversion(:)
  Type(tethers_type), Allocatable        :: tether(:)
  Type(threebody_type), Allocatable      :: threebody(:)
  Type(z_density_type), Allocatable      :: zdensity(:)
  Type(constraints_type), Allocatable    :: cons(:)
  Type(neighbours_type), Allocatable     :: neigh(:)
  Type(pmf_type), Allocatable            :: pmfs(:)
  Type(site_type), Allocatable           :: sites(:)
  Type(core_shell_type), Allocatable     :: core_shells(:)
  Type(vdw_type), Allocatable            :: vdws(:)
  Type(tersoff_type), Allocatable        :: tersoffs(:)
  Type(four_body_type), Allocatable      :: fourbody(:)
  Type(rdf_type), Allocatable            :: rdf(:)
  Type(netcdf_param), Allocatable        :: netcdf(:)
  Type(minimise_type), Allocatable       :: minim(:)
  Type(mpole_type), Allocatable          :: mpoles(:)
  Type(external_field_type), Allocatable :: ext_field(:)
  Type(rigid_bodies_type), Allocatable   :: rigid(:)
  Type(electrostatic_type), Allocatable  :: electro(:)
  Type(domains_type), Allocatable        :: domain(:)
  Type(flow_type), Allocatable           :: flow(:)
  Type(seed_type), Allocatable           :: seed(:)
  Type(trajectory_type), Allocatable     :: traj(:)
  Type(kim_type), Allocatable, Target    :: kim_data(:)
  Type(configuration_type), Allocatable  :: config(:)
  Type(io_type), Allocatable             :: ios(:)
  Type(ttm_type), Allocatable            :: ttms(:)
  Type(rsd_type), Allocatable            :: rsdsc(:)
  Type(file_type), Allocatable           :: files(:, :)
  Type(coord_type), Allocatable          :: crd(:)
  Type(adf_type), Allocatable            :: adf(:)

  ! Local Variables
  Character(len=1024) :: control_filename = '', arg
  Character(len=1024) :: output_filename = ''
  Logical             :: finish
  Integer             :: i

  ! SET UP COMMUNICATIONS & CLOCKING

  Allocate (dlp_world(0:0))
  Call init_comms(dlp_world(0))
  !dlp_world(0)%ou=nrite
  !Call init_error_system(nrite,dlp_world(0))
  Call gsync(dlp_world(0))

  finish = .false.
  If (dlp_world(0)%idnode == 0) Then
    If (command_argument_count() > 0) Then
      i = 0
      Do
        i = i + 1
        Call get_command_argument(i, arg)
        Select Case (Trim (arg))
        Case ('-h')
          Call get_command_argument(0, arg)
          Write (eu, '(a)') "Usage: "//Trim(arg)//" -c CONTROL_FILENAME -o OUTPUT_FILENAME"
          Write (eu, '(a)') "Each of -c or -o options are optional"
          Write (eu, '(a)') "use -h to see this help"
          finish = .true.
          Exit
        Case ('-c')
          i = i + 1
          Call get_command_argument(i, control_filename)
        Case ('-o')
          i = i + 1
          Call get_command_argument(i, output_filename)
        Case default
          Write (eu, *) "No idea what you want, try -h "
          finish = .true.
          Exit
        End Select
        If (i == command_argument_count()) Exit
      End Do
    End If
  End If

  Call gbcast(dlp_world(0), finish, 0)
  If (finish) Then
    Call exit_comms(dlp_world)
    Stop 0
  End If

  ! IS: This has to be abstracted or defined to be of dimension 1 in module flow. 
  Allocate(flow(1)) 

  ! Set the type of calculation to be performed. By default it is the standard DL_POLY 
  ! calculation. Tag evb activates EVB calculation
  Call read_simtype(flow(1),dlp_world(0))

  ! Select metasimulation method
  ! IS: The following two subroutines should be merged into a single one. We separate them
  ! for the time being though.
  If (flow(1)%simulation_method == MD) Then
    Call molecular_dynamics(dlp_world, thermo, ewld, tmr, devel, stats, &
                            green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, tether, &
                            threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, fourbody, &
                            rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, flow, seed, traj, &
                            kim_data, config, ios, ttms, rsdsc, files, output_filename, control_filename, crd, adf)

  Else If (flow(1)%simulation_method == EmpVB) Then
    Write (0, *) "simulation type: EVB"
   Call evb_molecular_dynamics(dlp_world, thermo, ewld, tmr, devel, stats, &
                            green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, tether, &
                            threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, fourbody, &
                            rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, flow, seed, traj, &
                            kim_data, config, ios, ttms, rsdsc, files, output_filename, control_filename, crd, adf)
    Else If (flow(1)%simulation_method == FFS) Then 
      write(0,*) "simulation type: FFS" 
    Else
      Write (0, *) "Unknown simulation type"
  End If

  ! Terminate job

  Call gsync(dlp_world(0))
  Call exit_comms(dlp_world)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate (flow)
  Deallocate (dlp_world)

End Program dl_poly
