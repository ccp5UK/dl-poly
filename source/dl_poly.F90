Program dl_poly

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 is an stfc/ccp5 program package for the dynamical
  ! simulation of molecular systems.
  !
  ! dl_poly_4 is based on dl_poly_3 by i.t.todorov & w.smith.
  !
  ! copyright - daresbury laboratory
  ! license   - LGPL 3.0;  https://www.gnu.org/licenses/lgpl-3.0.en.html
  ! authors   - i.t.todorov & w.smith april 2020
  ! contrib   - i.j.bush, h.a.boateng, m.a.seaton,
  !             a.brukhno, a.m.elena, r.davidchak,
  !             s.l.daraszewicz, g.khara, s.t.murphy,
  !             a.b.g.chalk, i.scivetti
  !
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  ! EVB       - i.scivetti march-october 2019
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use angles,                                 Only: angles_type
  Use angular_distribution,                   Only: adf_type
  Use bonds,                                  Only: bonds_type
  Use comms,                                  Only: comms_type,&
                                                    exit_comms,&
                                                    gbcast,&
                                                    gsync,&
                                                    init_comms
  Use configuration,                          Only: configuration_type
  Use constraints,                            Only: constraints_type
  Use control,                                Only: read_simtype
  Use control_parameter_module,               Only: dump_parameters,&
                                                    parameters_hash_table
  Use coord,                                  Only: coord_type
  Use core_shell,                             Only: core_shell_type
  Use defects,                                Only: defects_type
  Use development,                            Only: development_type
  Use dihedrals,                              Only: dihedrals_type
  Use domains,                                Only: domains_type
  Use electrostatic,                          Only: electrostatic_type
  Use errors_warnings,                        Only: init_error_system
  Use ewald,                                  Only: ewald_type
  Use external_field,                         Only: external_field_type
  Use filename,                               Only: FILENAME_SIZE,&
                                                    FILE_CONTROL,&
                                                    file_type
  Use flow_control,                           Only: EmpVB,&
                                                    FFS,&
                                                    MD_STD,&
                                                    flow_type
  Use four_body,                              Only: four_body_type
  Use greenkubo,                              Only: greenkubo_type
  Use hash,                                   Only: STR_LEN
  Use impacts,                                Only: impact_type
  Use inversions,                             Only: inversions_type
  Use io,                                     Only: io_type
  Use, Intrinsic :: iso_fortran_env,          Only: eu => error_unit,&
                                                    ou => output_unit
  Use kim,                                    Only: kim_type
  Use meta,                                   Only: molecular_dynamics
  Use metal,                                  Only: metal_type
  Use minimise,                               Only: minimise_type
  Use mpole,                                  Only: mpole_type
  Use msd,                                    Only: msd_type
  Use neighbours,                             Only: neighbours_type
  Use new_control,                            Only: initialise_control,&
                                                    read_new_control
  Use numerics,                               Only: seed_type
  Use plumed,                                 Only: plumed_type
  Use pmf,                                    Only: pmf_type
  Use poisson,                                Only: poisson_type
  Use rdfs,                                   Only: rdf_type
  Use rigid_bodies,                           Only: rigid_bodies_type
  Use rsds,                                   Only: rsd_type
  Use site,                                   Only: site_type
  Use statistics,                             Only: stats_type
  Use tersoff,                                Only: tersoff_type
  Use tethers,                                Only: tethers_type
  Use thermostat,                             Only: thermostat_type
  Use three_body,                             Only: threebody_type
  Use timer,                                  Only: timer_type
  Use trajectory,                             Only: trajectory_type
  Use ttm,                                    Only: ttm_type
  Use unit_test,                              Only: testing_type
  Use units,                                  Only: initialise_units
  Use vdw,                                    Only: vdw_type
  Use z_density,                              Only: z_density_type
#ifdef NVIDIA
  Use constants,                              Only: wp, &
                                                    half_minus, &
                                                    half_plus
#endif

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

  Type(parameters_hash_table) :: params
  Type(testing_type) :: tests

  ! Local Variables
  Character(len=1024)           :: control_filename = '', arg
  Character(len=1024)           :: output_filename = ''
  Character(Len=STR_LEN)        :: option
  Character(Len=10)             :: mode
  Logical                       :: finish
  Integer                       :: i, ifile

#ifdef NVIDIA
  half_plus = Nearest(0.5_wp, +1.0_wp)
  half_minus = Nearest(0.5_wp, -1.0_wp)
#endif

  ! SET UP COMMUNICATIONS & CLOCKING

  Allocate (dlp_world(0:0))
  Call init_comms(dlp_world(0))
  !dlp_world(0)%ou=nrite
  !Call init_error_system(nrite,dlp_world(0))
  Call gsync(dlp_world(0))

  ! temporary stuff this will need to be abstracted
  Allocate (flow(1))
  Allocate (devel(1))
  Call initialise_control(params)
  Call initialise_units()

  ! Assume we're running
  flow(1)%simulation = .true.
  ! Assume we're using old format
  finish = .false.
  If (dlp_world(0)%idnode == 0) Then
    If (command_argument_count() > 0) Then
      i = 0
      parse_cmd: Do
        i = i + 1
        Call get_command_argument(i, arg)
        Select Case (arg)
        Case ('-h')
          Call get_command_argument(0, arg)
          Write (eu, '(a)') "Usage: "//Trim(arg)//" -c CONTROL_FILENAME -o OUTPUT_FILENAME"
          Write (eu, '(a)') "Each of -c or -o options are optional"
          Write (eu, '(a)') "use -h to see this help"
          Write (eu, '(a)') "use --help to search help"
          finish = .true.
          Exit
        Case ('--dump')
          i = i + 1
          Call get_command_argument(i, mode)
          If (mode == '') Then
            mode = "default"
          End IF
          Select Case (mode)
          Case ('latexdoc', 'latex', 'python', 'csv', 'test', 'default')
            Continue
          Case Default
            Write (eu, '(a)') 'Bad mode option '//Trim(mode)
            finish = .true.
            Exit
          End Select

          Select Case (mode)
          Case ('latexdoc', 'latex', 'python', 'csv', 'test')
            i = i + 1
            Call get_command_argument(i, arg)
          End Select

          If (Trim(arg) == "SCREEN" .or. arg == '') Then
            ifile = ou
          Else
            Open (newunit=ifile, file=Trim(arg))
          End If

          Call dump_parameters(ifile, params, mode)
          finish = .true.
          Exit
        Case ('--help')
          i = i + 1
          Call get_command_argument(i, arg)
          If (arg == '' ) Then
            mode = 'default'
            Call dump_parameters(ou, params, mode)
          Else
            Call params%help(arg)
          End IF
          finish = .true.
          Exit

        Case ('--keywords', '-k')
          Call params%help()
          finish = .true.
          Exit

        Case ('-c')
          i = i + 1
          Call get_command_argument(i, control_filename)
        Case ('-o')
          i = i + 1
          Call get_command_argument(i, output_filename)
        Case ('-test', '-t')

          Call get_command_argument(i + 1, arg)
          Do While (arg(1:1) /= "-" .and. i < command_argument_count())
            i = i + 1
            Select Case (arg)
            Case ("control")
              tests%control = .true.
            Case ("configuration")
              tests%configuration = .true.
            Case ("units")
              tests%units = .true.
            Case ("vdw")
              tests%vdw = .true.
            Case ("all")
              Call tests%all()
            Case Default
              Write (eu, *) "Invalid test option:", Trim(arg)
              finish = .true.
              Exit parse_cmd
            End Select
            Call get_command_argument(i + 1, arg)
          End Do

          Call init_error_system(ou, dlp_world(0))
          Call tests%run(dlp_world(0))
          finish = .true.

        Case ('--replay', '-r')
          flow(1)%simulation = .false.
        Case default
          Write (eu, *) "No idea what you want, try -h "
          finish = .true.
          Exit parse_cmd
        End Select
        If (i == command_argument_count()) Exit
      End Do parse_cmd
    End If
  End If

  Call gbcast(dlp_world(0), control_filename, 0)
  Call gbcast(dlp_world(0), output_filename, 0)
  Call gbcast(dlp_world(0), finish, 0)
  If (finish) Then
    Call exit_comms(dlp_world)
    Stop 0
  End If

  Allocate (files(1, FILENAME_SIZE))
  ! Rename control file if argument was passed
  If (Len_trim(control_filename) > 0) Then
    Call files(1, FILE_CONTROL)%rename(control_filename)
  Else
    Call files(1, FILE_CONTROL)%rename('CONTROL')
  End If

  ! Temporary error system
  Call init_error_system(eu, dlp_world(0))
  Call read_new_control(files(1, FILE_CONTROL), params, dlp_world(0), devel(1)%new_control)

  If (devel(1)%new_control) Then
    Call params%retrieve('simulation_method', option)
    Select Case (option)
    Case ('md')
      flow(1)%simulation_method = MD_STD
      flow(1)%NUM_FF = 1
    Case ('evb')
      flow(1)%simulation_method = EmpVB
      Call params%retrieve('evb_num_ff', flow(1)%NUM_FF)
    Case ('ffs')
      flow(1)%simulation_method = FFS
    Case Default
      flow(1)%simulation_method = -1
    End Select

  Else ! Cannot read as new style

    ! Set the type of calculation to be performed. By default it is the standard DL_POLY
    ! calculation. Tag evb activates EVB calculation
    Call read_simtype(control_filename, flow(1), dlp_world(0))

  End If

  ! Select metasimulation method
  ! IS: The following two subroutines should be merged into a single one. We separate them
  ! for the time being though.
  Select Case (flow (1)%simulation_method)
  Case (MD_STD, EmpVB)
    Call molecular_dynamics(params, dlp_world, thermo, ewld, tmr, devel, stats, &
                            green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, tether, &
                            threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, fourbody, &
                            rdf, minim, mpoles, ext_field, rigid, electro, domain, flow, seed, traj, &
                            kim_data, config, ios, ttms, rsdsc, files, output_filename, control_filename, crd, adf)
  Case (FFS)
    Write (0, *) "simulation type: FFS"
  Case Default
    Write (0, *) "Unknown simulation type"
  End Select

  ! Terminate job

  Call gsync(dlp_world(0))
  Call exit_comms(dlp_world)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate (flow)
  Deallocate (dlp_world)

End Program dl_poly
