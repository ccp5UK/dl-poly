! YL Oct 2021: adapted from dl_poly.F90 for interfacing DL_POLY 5
!              as a general-purpose driver
Module DLPOLYModule

  Private

  Public :: dl_poly

  Real(Kind=8)         , Public :: dl_poly_energy
  ! fxx/fyy/fzz is not TARGET
  Real(Kind=8), Pointer, Public :: dl_poly_gradients(:,:)
  Integer     , Pointer, Public :: dl_poly_global_indices(:)
  Integer                       :: ff
  Character(len=1024)  , Public :: control_filename = ''
  Character(len=1024)  , Public :: output_filename = ''

  Contains

  Subroutine dl_poly(comm_external)

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
  Use control,                            Only: read_simtype
  Use control_parameter_module,           Only: parameters_hash_table,&
                                                dump_parameters
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
  Use filename,                           Only: file_type,&
                                                FILENAME_SIZE,&
                                                FILE_CONTROL
  Use flow_control,                       Only: EmpVB,&
                                                FFS,&
                                                MD_STD,&
                                                flow_type
  Use four_body,                          Only: four_body_type
  Use greenkubo,                          Only: greenkubo_type
  Use hash,                               Only: STR_LEN
  Use impacts,                            Only: impact_type
  Use inversions,                         Only: inversions_type
  Use io,                                 Only: io_type
  Use, Intrinsic :: iso_fortran_env,      Only: eu => error_unit, &
                                                ou => output_unit
  Use kim,                                Only: kim_type
  ! YL 14/08/2021: needed for interfacing
  Use meta,                               Only: molecular_dynamics
  Use metal,                              Only: metal_type
  Use minimise,                           Only: minimise_type
  Use mpole,                              Only: mpole_type
  Use msd,                                Only: msd_type
  Use neighbours,                         Only: neighbours_type
  Use netcdf_wrap,                        Only: netcdf_param
  Use new_control,                        Only: initialise_control, &
                                                read_new_control
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
  Use unit_test,                          Only: testing_type
  Use vdw,                                Only: vdw_type
  Use z_density,                          Only: z_density_type

  Use units,                              Only: initialise_units

  Implicit None

! External MPI communicator
  Integer       :: comm_external

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

  Type( parameters_hash_table ) :: params
  Type( testing_type ) :: tests

  ! Local Variables
  Character(len=1024)           :: arg
  Character(Len=STR_LEN)        :: option
  Character(Len=10)             :: mode
  Logical                       :: finish
  Integer                       :: i, ifile

  ! SET UP COMMUNICATIONS & CLOCKING

  Allocate (dlp_world(0:0))
  Call init_external_comms(dlp_world(0), comm_external)
  !dlp_world(0)%ou=nrite
  !Call init_error_system(nrite,dlp_world(0))
  Call gsync(dlp_world(0))

  ! temporary stuff this will need to be abstracted
  Allocate (flow(1))
  Allocate(devel(1))
  Call initialise_control(params)
  call initialise_units()

  ! Assume we're running
  flow(1)%simulation = .true.
  ! Assume we're using old format
  finish = .false.

  Call gbcast(dlp_world(0), control_filename, 0)
  Call gbcast(dlp_world(0), output_filename, 0)
  Call gbcast(dlp_world(0), finish, 0)
  If (finish) Then
! YL: do not let it kill the process (`stop`)
    Call exit_external_comms(dlp_world)
!    Stop 0
  End If

  Allocate(files(1,FILENAME_SIZE))
      ! Rename control file if argument was passed
    If (Len_Trim(control_filename) > 0 ) Then
       Call files(1,FILE_CONTROL)%rename(control_filename)
    Else
       Call files(1,FILE_CONTROL)%rename('CONTROL')
    End If


    ! Temporary error system
    Call init_error_system(eu, dlp_world(0))
    call read_new_control(files(1,FILE_CONTROL), params, dlp_world(0), devel(1)%new_control)

    if (devel(1)%new_control) then
      Call params%retrieve('simulation_method', option)
      Select Case (option)
      Case ('md')
        flow(1)%simulation_method=MD_STD
        flow(1)%NUM_FF = 1
      Case ('evb')
        flow(1)%simulation_method=EmpVB
        Call params%retrieve('evb_num_ff', flow(1)%NUM_FF)
      Case ('ffs')
        flow(1)%simulation_method=FFS
      Case Default
        flow(1)%simulation_method=-1
      End Select

    else     ! Cannot read as new style

      ! Set the type of calculation to be performed. By default it is the standard DL_POLY
       ! calculation. Tag evb activates EVB calculation
       Call read_simtype(control_filename, flow(1), dlp_world(0))

    end if



  ! Select metasimulation method
  ! IS: The following two subroutines should be merged into a single one. We separate them
  ! for the time being though.
  Select Case (flow(1)%simulation_method)
  Case (MD_STD, EmpVB)
    Call molecular_dynamics(params, dlp_world, thermo, ewld, tmr, devel, stats, &
                            green, plume, msd_data, met, pois, impa, dfcts, bond, angle, dihedral, inversion, tether, &
                            threebody, zdensity, cons, neigh, pmfs, sites, core_shells, vdws, tersoffs, fourbody, &
                            rdf, netcdf, minim, mpoles, ext_field, rigid, electro, domain, flow, seed, traj, &
                            kim_data, config, ios, ttms, rsdsc, files, output_filename, control_filename, crd, adf)
  Case (FFS)
     write(0,*) "simulation type: FFS"
  Case Default
     Write (0, *) "Unknown simulation type"
  End Select

  ! Terminate job

  Call gsync(dlp_world(0))
  Call exit_external_comms(dlp_world)

  Do ff = 1, flow(1)%NUM_FF
      If (dlp_world(0)%idnode == 0) Then
          dl_poly_energy = stats(ff)%engcpe + stats(ff)%engsrp + stats(ff)%engtbp + stats(ff)%engbnd + stats(ff)%engang + stats(ff)%engdih + stats(ff)%enginv
          ! Print the energy breakdown
          Call print_energies(stats(ff))
      End If
      ! get the gradients
      Allocate(dl_poly_gradients(3,config(ff)%natms))
      If(Associated(dl_poly_gradients)) Then
          dl_poly_gradients(1,:) = config(ff)%parts(1:config(ff)%natms)%fxx
          dl_poly_gradients(2,:) = config(ff)%parts(1:config(ff)%natms)%fyy
          dl_poly_gradients(3,:) = config(ff)%parts(1:config(ff)%natms)%fzz
      End If
      Allocate(dl_poly_global_indices(config(ff)%natms))
      If(Associated(dl_poly_global_indices)) Then
          dl_poly_global_indices(1:config(ff)%natms) = config(ff)%ltg(1:config(ff)%natms)
      End If
  End Do
  ! the deallocation of config and stats are commented out in source/meta.F90
  Deallocate(config, stats)

  ! Create wrappers for the MD cycle in VV, and replay history
  Deallocate (flow)
  Deallocate (dlp_world)

  contains
  Subroutine init_external_comms(comm, comm_external)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for initialising the communication harness
! determining the MPI precision, and node identification and count
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
! contrib   - a.m.elena march 2016
! refactoring:
!           - a.m.elena march-october 2018
!           - j.madge march-october 2018
!           - a.b.g.chalk march-october 2018
!           - i.scivetti march-october 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef WITH_MPI
    Use mpi
#else
    Use mpi_api
#endif
    Use particle, Only: corePart
    Use kinds   , Only: dp, qp, si, sp, wi, wp
    Use comms   , Only: lib_version, mpi_ver, mpi_subver, proc_name, root_id, wp_mpi
    Use, Intrinsic :: iso_fortran_env, Only: error_unit

    Implicit None

    Type(comms_type), Intent(InOut) :: comm
    Integer         , Intent(In)    :: comm_external

    Integer                        :: lname, lversion
    Integer(KIND=MPI_ADDRESS_KIND) :: base, displacements(1:9), extent, lb
    Integer, Dimension(1:9)        :: block_lengths, types
    Type(corePart)                 :: part_array(1:5), part_temp

! YL TODO 14/08/2021: we currently don't support this
!#ifdef WITH_DFTBP
!    Integer, Parameter :: requiredThreading=MPI_THREAD_FUNNELED
!    Integer, Parameter :: errorcode = 0
!    Integer            :: providedThreading
!    Integer            :: idnode
!
!    Call MPI_INIT_THREAD(requiredThreading, providedThreading, comm%ierr)
!    If( providedThreading < requiredThreading )Then
!       Call MPI_COMM_RANK(comm_external, idnode, comm%ierr)
!       if(idnode == root_id)then
!          Write(error_unit,'(1x,a,I1,a,I1)') 'error - MPI library threading support, ',&
!               providedThreading,', is less than that required by DFTB+ v18.2, ',requiredThreading
!          Call MPI_ABORT(comm_external, errorcode, comm%ierr)
!       End if
!    End If
!#else
!!    Call MPI_INIT(comm%ierr)
!#endif

    Call MPI_COMM_DUP(comm_external, comm%comm, comm%ierr)

    ! use iso_fortran_env
    If (wp == sp) Then
      wp_mpi = MPI_REAL
    Else If (wp == dp) Then
      ! MPI_REAL8 is apparently not in the strict MPI2 standard
      ! It is just an optional data type in the FORTRAN Bindings
      wp_mpi = MPI_DOUBLE_PRECISION
    Else If (wp == qp) Then
      wp_mpi = MPI_REAL16
    Else
      Write (0, '(/,1x,a)') 'error - working precision mismatch between FORTRAN90 and MPI implementation'
      Call abort_comms(comm, 1000)
    End If
    Call MPI_COMM_RANK(comm%comm, comm%idnode, comm%ierr)
    Call MPI_COMM_SIZE(comm%comm, comm%mxnode, comm%ierr)

    !Create the transfer type for the corePart type.
    block_lengths(1:9) = 1
    Call MPI_GET_ADDRESS(part_temp%xxx, displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%yyy, displacements(2), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%zzz, displacements(3), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fxx, displacements(4), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fyy, displacements(5), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fzz, displacements(6), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%chge, displacements(7), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%pad1, displacements(8), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%pad2, displacements(9), comm%ierr)
    base = displacements(1)
    displacements(1:9) = displacements(1:9) - base
    types(1:7) = wp_mpi
    types(8:9) = MPI_INTEGER

    Call MPI_TYPE_CREATE_STRUCT(9, block_lengths, displacements, types, comm%part_type, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_type, comm%ierr)

    Call MPI_GET_ADDRESS(part_array(1), displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_array(2), displacements(2), comm%ierr)
    extent = displacements(2) - displacements(1)
    lb = 0
    Call MPI_TYPE_CREATE_RESIZED(comm%part_type, lb, extent, comm%part_array_type, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_array_type, comm%ierr)

    Call MPI_GET_ADDRESS(part_temp%xxx, displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%yyy, displacements(2), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%zzz, displacements(3), comm%ierr)
    base = displacements(1)
    block_lengths(1:3) = 1
    types(1:3) = wp_mpi
    Call MPI_TYPE_CREATE_STRUCT(3, block_lengths, displacements, types, comm%part_type_positions, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_type_positions, comm%ierr)
    Call MPI_GET_ADDRESS(part_array(1), displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_array(2), displacements(2), comm%ierr)
    extent = displacements(2) - displacements(1)
    lb = 0
    Call MPI_TYPE_CREATE_RESIZED(comm%part_type_positions, lb, extent, comm%part_array_type_positions, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_array_type_positions, comm%ierr)

    Call MPI_GET_ADDRESS(part_temp%fxx, displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fyy, displacements(2), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fzz, displacements(3), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%xxx, base, comm%ierr)
    displacements(1:3) = displacements(1:3) - base
    Call MPI_TYPE_CREATE_STRUCT(3, block_lengths, displacements, types, comm%part_type_forces, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_type_forces, comm%ierr)
    Call MPI_GET_ADDRESS(part_array(1), displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_array(2), displacements(2), comm%ierr)
    extent = displacements(2) - displacements(1)
    lb = 0
    Call MPI_TYPE_CREATE_RESIZED(comm%part_type_forces, lb, extent, comm%part_array_type_forces, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_array_type_forces, comm%ierr)
! YL 14/08/2021: I don't think this is needed(?)
!#ifndef OLDMPI
!    Call MPI_GET_PROCESSOR_NAME(proc_name, lname, comm%ierr)
!    Call MPI_GET_VERSION(mpi_ver, mpi_subver, comm%ierr)
!    Call MPI_GET_LIBRARY_VERSION(lib_version, lversion, comm%ierr)
!#endif

  End Subroutine init_external_comms
  Subroutine exit_external_comms(comm)

#ifdef WITH_MPI
    Use mpi
#else
    Use mpi_api , Only: MPI_TYPE_FREE
#endif
    Use kinds   , Only: wi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for exiting the communication harness
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm(:)

    Integer(Kind=wi) :: i

    Do i = 1, Size(comm, dim=1)
      Call MPI_TYPE_FREE(comm(i)%part_array_type, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_type, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_array_type_positions, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_type_positions, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_array_type_forces, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_type_forces, comm(i)%ierr)
! YL 14/08/2021: commented out for interfacing
!      Call MPI_FINALIZE(comm(i)%ierr)
    End Do

  End Subroutine exit_external_comms
  Subroutine print_energies(stats)

    Type(stats_type), intent(in) :: stats

    Write(*,'(/A)') " ----------------------------------------------------------------------"
    Write(*,'(A)') " DL_POLY 5 Interface"
    Write(*,'(A/)') " ----------------------------------------------------------------------"
    
    ! YL NB: internal energy in 10 J/mol!
    Write(*,'(A)') " Breakdown of DL_POLY energies:"
    Write(*,'(31X,A,17X,A)') " a.u.", "kcal/mol"
    ! coulombic
    Write(*,'(A,F23.14,1X,F25.14)') " Coulombic   (engcpe) ", stats%engcpe/627.524955413*10.0/4.184/1000.0, stats%engcpe*10.0/4.184/1000.0
    ! short range
    Write(*,'(A,F23.14,1X,F25.14)') " Short range (engsrp) ", stats%engsrp/627.524955413*10.0/4.184/1000.0, stats%engsrp*10.0/4.184/1000.0
    ! 3-body
    Write(*,'(A,F23.14,1X,F25.14)') " 3-body      (engtbp) ", stats%engtbp/627.524955413*10.0/4.184/1000.0, stats%engtbp*10.0/4.184/1000.0
    ! bond
    Write(*,'(A,F23.14,1X,F25.14)') " Bond        (engbnd) ", stats%engbnd/627.524955413*10.0/4.184/1000.0, stats%engbnd*10.0/4.184/1000.0
    ! angle
    Write(*,'(A,F23.14,1X,F25.14)') " Angle       (engang) ", stats%engang/627.524955413*10.0/4.184/1000.0, stats%engang*10.0/4.184/1000.0
    ! dihedral
    Write(*,'(A,F23.14,1X,F25.14)') " Dihedral    (engdih) ", stats%engdih/627.524955413*10.0/4.184/1000.0, stats%engdih*10.0/4.184/1000.0
    ! inversion
    Write(*,'(A,F23.14,1X,F25.14)') " Inversion   (enginv) ", stats%enginv/627.524955413*10.0/4.184/1000.0, stats%enginv*10.0/4.184/1000.0
    ! engstbnd, enganan, engaat, engub missing in this formula
    Write(*,'(A)') " (Not including engstbnd, enganan, engaat, or engub)"
    Write(*,*)
    Write(*,'(A,F23.14,1X,F25.14)') " Total                ", dl_poly_energy/627.524955413*10.0/4.184/1000.0, dl_poly_energy*10.0/4.184/1000.0
    Write(*,'(A)') " (engcpe + engsrp + engtbp + engbnd + engang + engdih + enginv)"
    Write(*,*)
    Write(*,'(A)') " Other terms:"
    Write(*,'(31X,A,17X,A)') " a.u.", "kcal/mol"
    Write(*,'(A,F23.14,1X,F25.14)') " Kinetics    (engke)  ", stats%engke /627.524955413*10.0/4.184/1000.0, stats%engke *10.0/4.184/1000.0
    Write(*,'(A,F23.14,1X,F25.14)') " Tersoff     (engter) ", stats%engter/627.524955413*10.0/4.184/1000.0, stats%engter*10.0/4.184/1000.0
    Write(*,'(A,F23.14,1X,F25.14)') " 4-body      (engfbp) ", stats%engfbp/627.524955413*10.0/4.184/1000.0, stats%engfbp*10.0/4.184/1000.0
    Write(*,'(A,F23.14,1X,F25.14)') " Core-shell  (engshl) ", stats%engshl/627.524955413*10.0/4.184/1000.0, stats%engshl*10.0/4.184/1000.0
    Write(*,'(A,F23.14,1X,F25.14)') " Tether      (engtet) ", stats%engtet/627.524955413*10.0/4.184/1000.0, stats%engtet*10.0/4.184/1000.0
    Write(*,'(A,F23.14,1X,F25.14)') " Field       (engfld) ", stats%engfld/627.524955413*10.0/4.184/1000.0, stats%engfld*10.0/4.184/1000.0
  End Subroutine

End Subroutine dl_poly

End Module DLPOLYModule
