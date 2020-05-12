Module configuration

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global configuration variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov january 2017
  ! contrib   - i.j.bush march 2010
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use asserts,         Only: assert
  Use comms,           Only: &
                             WriteConf_tag, comm_self, comms_type, gallgather, gallgatherv, &
                             gallreduce, galltoall, galltoallv, gatherv_scatterv_index_arrays, &
                             gbcast, gcheck, ggatherv, gmax, gmin, grecv, gscatter, &
                             gscatter_columns, gscatterv, gsend, gsum, gsync, mode_create, &
                             mode_rdonly, mode_wronly, mpi_distribution_type, offset_kind, op_land
  Use constants,       Only: half_minus,&
                             zero_plus
  Use development,     Only: development_type
  Use domains,         Only: domains_type
  Use electrostatic,   Only: ELECTROSTATIC_EWALD,&
                             ELECTROSTATIC_NULL
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use filename,        Only: FILE_CONFIG,&
                             file_type
  Use flow_control,    Only: RESTART_KEY_CLEAN,&
                             flow_type
  Use io,              Only: &
                             IO_ALLOCATION_ERROR, IO_BASE_COMM_NOT_SET, IO_READ_MASTER, &
                             IO_READ_NETCDF, IO_RESTART, IO_SUBSET_FORCES, IO_SUBSET_POSITIONS, &
                             IO_UNKNOWN_WRITE_LEVEL, IO_UNKNOWN_WRITE_OPTION, &
                             IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
                             IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_NETCDF, &
                             IO_WRITE_UNSORTED_DIRECT, IO_WRITE_UNSORTED_MASTER, &
                             IO_WRITE_UNSORTED_MPIIO, io_close, io_delete, io_finalize, &
                             io_get_parameters, io_get_var, io_init, io_nc_create, io_nc_get_att, &
                             io_nc_get_dim, io_nc_get_var, io_nc_put_var, io_open, io_read_batch, &
                             io_set_parameters, io_type, io_write_batch, io_write_record, &
                             io_write_sorted_file, recsz
  Use kinds,           Only: li,&
                             wi,&
                             wp
  Use kpoints,         Only: kpoints_type
  Use netcdf_wrap,     Only: netcdf_param
  Use numerics,        Only: dcell,&
                             images,&
                             invert,&
                             pbcshift,&
                             shellsort,&
                             shellsort2
  Use parse,           Only: get_line,&
                             get_word,&
                             strip_blanks,&
                             tabs_2_blanks,&
                             word_2_real
  Use particle,        Only: corePart
  Use site,            Only: site_type
  Use thermostat,      Only: CONSTRAINT_NONE,&
                             thermostat_type

  Implicit None

  Private

  Integer, Public, Parameter ::  len_atmnam = 8

  Type, Public :: configuration_type

    Character(Len=72)             :: cfgname = ' ', &
                                     sysname = ' '
    Integer                       :: imcon = -1, &
                                     imc_n = -1, &
                                     natms = 0, &
                                     nlast = 0, &
                                     nfree = 0
    Real(Kind=wp)                 :: cell(1:9) = 0.0_wp, &
                                     volm = 0.0_wp, &
                                     sumchg = 0.0_wp
    Character(Len=len_atmnam), Allocatable       :: atmnam(:)
    Integer, Allocatable          :: lsite(:), ltype(:)
    Integer, Allocatable          :: lfrzn(:), lfree(:)
    Integer, Allocatable          :: lsi(:), lsa(:), ltg(:)
    Integer, Allocatable          :: ixyz(:)
    Integer, Allocatable          :: lstfre(:)
    Real(Kind=wp), Allocatable    :: weight(:) !,chge(:)
    !  Real( Kind = wp ),    Allocatable, Save :: xxx(:),yyy(:),zzz(:)
    Real(Kind=wp), Allocatable    :: vxx(:), vyy(:), vzz(:)
    !  Real( Kind = wp ),    Allocatable, Save :: fxx(:),fyy(:),fzz(:)
    Type(corePart), Allocatable   :: parts(:)
    Logical                       :: newjob_check_config = .true.
    Logical                       :: newjob_totmas = .true.
    Logical                       :: newjob_totmas_r = .true.
    Logical                       :: newjob_meg = .true.
    Real(Kind=wp)                 :: totmas
    Real(Kind=wp)                 :: totmas_r
    Real(Kind=wp)                 :: meg
    Logical, Public               :: l_vom = .true., &
                                     lvom = .true. ! this is confusing and needless complicated
    Integer(Kind=wi), Public      :: mxtana, mxgana, mxbfss, mxbuff
    Integer(Kind=wi), Public      :: mxlshp, mxatms, mxatdm
    ! general flags
    Logical                       :: l_ind, l_exp
    Integer                       :: levcfg, nx, ny, nz, &
                                     atmfre, atmfrz, megatm, megfrz
    ! Degrees of freedom must be in long integers so we do 2.1x10^9 particles
    Integer(Kind=li)              :: degfre, degshl, degtra, degrot
    ! vdws%elrc,vdws%vlrc - vdw energy and virial are scalars and in vdw
    Real(Kind=wp)                 :: dvar, fmax, width
    Type(kpoints_type)            :: k

  Contains

    Private

    Procedure, Public :: chvom
    Procedure, Public :: init => allocate_config_arrays
    Procedure, Public :: init_read => allocate_config_arrays_read
    Final             :: deallocate_config_arrays

  End Type configuration_type

  !> Coordinates of all atoms in the system, stored in packed form,
  !> with mpi index arrays to allow unpacking
  !
  Type, Public :: coordinate_buffer_type
    ! Atomic positions from all processes stored in a packed form
    Real(wp), Public, Allocatable :: coords(:)
    ! MPI indexing required to unpack coordinates stored in packed form
    Type(mpi_distribution_type), Public :: mpi
    ! Atomic names/species
    ! Uses different MPI indexing to coords, hence maybe makes sense
    ! not to include in this type
    !Character(Len=len_atmnam), Public, Allocatable :: atmnam(:)
  Contains
    Procedure, Public :: initialise => initialise_coordinate_buffer_type
    Procedure, Public :: finalise => finalise_coordinate_buffer_type
  End Type coordinate_buffer_type

  Interface reallocate
    Module Procedure reallocate_chr_v
    Module Procedure reallocate_int_v
    Module Procedure reallocate_rwp_v
    Module Procedure reallocate_corePart_v
  End Interface

  Interface getcom
    Module Procedure getcom_parts
    Module Procedure getcom_arrays
  End Interface getcom

  Public :: reallocate
  Public :: check_config
  Public :: read_config, read_config_parallel
  Public :: scan_config
  Public :: scale_config
  Public :: origin_config
  Public :: write_config
  Public :: freeze_atoms
  Public :: getcom, getcom_mol
  Public :: gather_coordinates, gather_atomic_names
  Public :: unpack_gathered_coordinates
  Public :: distribute_forces

Contains

  Pure Subroutine chvom(T, flag)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to change behaviour for COM momentum removal
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov july 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Class(configuration_type), Intent(InOut) :: T
    Logical, Optional,         Intent(In   ) :: flag

    Logical :: lflag

    lflag = T%l_vom
    If (Present(flag)) lflag = flag

    If (lflag) Then
      T%lvom = .true. ! Remove COM momentum
    Else
      T%lvom = .false. ! Don't
    End If

  End Subroutine chvom

  Subroutine reallocate_corePart_v(delta, a, stat)

    Integer,                     Intent(In   ) :: delta
    Type(corePart), Allocatable, Intent(InOut) :: a(:)
    Integer,                     Intent(  Out) :: stat

    Integer                     :: size_crs, size_new, size_old
    Type(corePart), Allocatable :: tmp(:)

    stat = 0
    If (delta == 0) Return

    size_old = Size(a)
    size_new = size_old + delta
    size_crs = Min(size_old, size_new)

    Allocate (tmp(1:size_new), Stat=stat)
    If (stat /= 0) Return

    tmp(1:size_crs) = a(1:size_crs)

    Deallocate (a, Stat=stat)
    Call move_alloc(tmp, a)

  End Subroutine

  Subroutine reallocate_chr_v(delta, a, stat)

    Integer,                       Intent(In   ) :: delta
    Character(Len=*), Allocatable, Intent(InOut) :: a(:)
    Integer,                       Intent(  Out) :: stat

    Character(Len=Len(a)), Allocatable :: tmp(:)
    Integer                            :: size_crs, size_new, size_old

    stat = 0
    If (delta == 0) Return

    size_old = Size(a)
    size_new = size_old + delta
    size_crs = Min(size_old, size_new)

    Allocate (tmp(1:size_new), Stat=stat)
    If (stat /= 0) Return

    tmp(1:size_crs) = a(1:size_crs)

    Deallocate (a, Stat=stat)
    Call move_alloc(tmp, a)

  End Subroutine reallocate_chr_v

  Subroutine reallocate_int_v(delta, a, stat)

    Integer,              Intent(In   ) :: delta
    Integer, Allocatable, Intent(InOut) :: a(:)
    Integer,              Intent(  Out) :: stat

    Integer              :: size_crs, size_new, size_old
    Integer, Allocatable :: tmp(:)

    stat = 0
    If (delta == 0) Return

    size_old = Size(a)
    size_new = size_old + delta
    size_crs = Min(size_old, size_new)

    Allocate (tmp(1:size_new), Stat=stat)
    If (stat /= 0) Return

    tmp(1:size_crs) = a(1:size_crs)

    Deallocate (a, Stat=stat)
    Call move_alloc(tmp, a)

  End Subroutine reallocate_int_v

  Subroutine reallocate_rwp_v(delta, a, stat)

    Integer,                    Intent(In   ) :: delta
    Real(Kind=wp), Allocatable, Intent(InOut) :: a(:)
    Integer,                    Intent(  Out) :: stat

    Integer                    :: size_crs, size_new, size_old
    Real(Kind=wp), Allocatable :: tmp(:)

    stat = 0
    If (delta == 0) Return

    size_old = Size(a)
    size_new = size_old + delta
    size_crs = Min(size_old, size_new)

    Allocate (tmp(1:size_new), Stat=stat)
    If (stat /= 0) Return

    tmp(1:size_crs) = a(1:size_crs)

    Deallocate (a, Stat=stat)
    If (stat /= 0) Return

    Allocate (a(1:size_new), Stat=stat)
    Call move_alloc(tmp, a)

  End Subroutine reallocate_rwp_v

  Subroutine allocate_config_arrays_read(config)

    Class(configuration_type), Intent(InOut) :: config

    Integer :: fail(1:4), i

    fail = 0
    Allocate (config%atmnam(1:config%mxatms), Stat=fail(1))
    Allocate (config%lsi(1:config%mxatms), config%lsa(1:config%mxatms), config%ltg(1:config%mxatms), Stat=fail(2))
    Allocate (config%parts(1:config%mxatms), Stat=fail(3))
    Allocate (config%vxx(1:config%mxatms), config%vyy(1:config%mxatms), config%vzz(1:config%mxatms), Stat=fail(4))
    If (Any(fail > 0)) Call error(1025)

    config%atmnam = ' '
    config%lsi = 0; config%lsa = 0; config%ltg = 0

    config%vxx = 0.0_wp; config%vyy = 0.0_wp; config%vzz = 0.0_wp
    Do i = 1, config%mxatms
      config%parts(i)%xxx = 0.0_wp
      config%parts(i)%yyy = 0.0_wp
      config%parts(i)%zzz = 0.0_wp
      config%parts(i)%fxx = 0.0_wp
      config%parts(i)%fyy = 0.0_wp
      config%parts(i)%fzz = 0.0_wp
      config%parts(i)%chge = 0.0_wp
    End Do

  End Subroutine allocate_config_arrays_read

  Subroutine allocate_config_arrays(config)

    Class(configuration_type), Intent(InOut) :: config

    Integer :: fail(1:5), stat(1:13)

    fail = 0

    Allocate (config%lsite(1:config%mxatms), config%ltype(1:config%mxatms), Stat=fail(1))
    Allocate (config%lfrzn(1:config%mxatms), config%lfree(1:config%mxatms), Stat=fail(2))
    Allocate (config%ixyz(1:config%mxatms), Stat=fail(3))
    Allocate (config%lstfre(1:config%mxatdm), Stat=fail(4))
    Allocate (config%weight(1:config%mxatms), Stat=fail(5))

    If (Any(fail > 0)) Call error(1025)

    config%lsite = 0; config%ltype = 0
    config%lfrzn = 0; config%lfree = 0

    config%ixyz = 0
    config%lstfre = 0

    config%weight = 0.0_wp

    ! Resize the arrays in allocate_config_arrays_read

    stat = 0

    Call reallocate(config%mxatms - Size(config%atmnam), config%atmnam, stat(1))
    Call reallocate(config%mxatms - Size(config%lsi), config%lsi, stat(2))
    Call reallocate(config%mxatms - Size(config%lsa), config%lsa, stat(3))
    Call reallocate(config%mxatms - Size(config%ltg), config%ltg, stat(4))
    Call reallocate(config%mxatms - Size(config%parts), config%parts, stat(5))
    Call reallocate(config%mxatms - Size(config%vxx), config%vxx, stat(6))
    Call reallocate(config%mxatms - Size(config%vyy), config%vyy, stat(7))
    Call reallocate(config%mxatms - Size(config%vzz), config%vzz, stat(8))

    If (Any(stat /= 0)) Call error(1025)

  End Subroutine allocate_config_arrays

  Subroutine deallocate_config_arrays(T)
    Type(configuration_type) :: T

    If (Allocated(T%atmnam)) Deallocate (T%atmnam)
    If (Allocated(T%lsi)) Deallocate (T%lsi)
    If (Allocated(T%lsa)) Deallocate (T%lsa)
    If (Allocated(T%ltg)) Deallocate (T%ltg)
    If (Allocated(T%parts)) Deallocate (T%parts)
    If (Allocated(T%vxx)) Deallocate (T%vxx)
    If (Allocated(T%vyy)) Deallocate (T%vyy)
    If (Allocated(T%vzz)) Deallocate (T%vzz)
    If (Allocated(T%lsite)) Deallocate (T%lsite)
    If (Allocated(T%ltype)) Deallocate (T%ltype)
    If (Allocated(T%lfrzn)) Deallocate (T%lfrzn)
    If (Allocated(T%lfree)) Deallocate (T%lfree)
    If (Allocated(T%ixyz)) Deallocate (T%ixyz)
    If (Allocated(T%lstfre)) Deallocate (T%lstfre)
    If (Allocated(T%weight)) Deallocate (T%weight)

  End Subroutine deallocate_config_arrays

  Subroutine check_config(config, electro_key, thermo, sites, flow, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reporting and checking the configuration data
    ! in CONFIG for: (i) unsuitable simulation options in CONTROL and
    ! (ii) connectivity to FIELD; before connecting the crystallographic
    ! data (positions+) to the topology (sites+), i.e. CONFIG to FIELD
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2015
    ! amended   - i.t.todorov november 2019 (global index printing for errors)
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent(In) :: electro_key
    Type(configuration_type), Intent(InOut) :: config
    Type(thermostat_type), Intent(In) :: thermo
    Type(site_type), Intent(In) :: sites
    Type(flow_type), Intent(In) :: flow
    Type(comms_type), Intent(InOut) :: comm

    Logical                :: safe
    Integer                :: fail, k, l, m, &
                              indatm, totatm, mol_sit, loc_ind
    Real(Kind=wp)      :: rcell(1:9), det

    Character(Len=256) :: message

    Integer, Allocatable :: iwrk(:)

    fail = 0
    If (flow%strict) Then
      Allocate (iwrk(1:config%mxatms), Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'check_config allocation failure'
        Call error(0, message)
      End If
    End If

    If (config%newjob_check_config) Then
      Write (message, "('configuration file name: ',10x,a)") config%cfgname
      Call info(message, .true.)
      Write (message, "('selected image convention',6x,i10)") config%imcon
      Call info(message, .true.)
    End If

    ! Check things for non-periodic systems

    If (config%imcon == 0 .or. config%imcon == 6) Then
      If (electro_key == ELECTROSTATIC_EWALD) Then
        Call warning(220, 0.0_wp, 0.0_wp, 0.0_wp)
      Else If (electro_key /= ELECTROSTATIC_NULL) Then
        Call warning(30, 0.0_wp, 0.0_wp, 0.0_wp)
      End If

      Call warning(260, 0.0_wp, 0.0_wp, 0.0_wp)

      If (thermo%variable_cell) Call error(390)
    End If

    ! Check and adapt image conditions for nst ensembles

    If (thermo%anisotropic_pressure) Then
      If (thermo%iso == CONSTRAINT_NONE) Then
        If (config%imcon == 1 .or. config%imcon == 2) Then
          Call warning(110, Real(config%imcon, wp), 3.0_wp, 0.0_wp)
          config%imcon = 3
        End If
      Else ! thermo%iso > 0
        If (config%imcon == 1) Then
          Call warning(110, Real(config%imcon, wp), 3.0_wp, 0.0_wp)
          config%imcon = 2
        End If
      End If
    End If

    ! Check image condition for pseudo

    If (thermo%l_stochastic_boundaries .and. (config%imcon == 0 .or. config%imcon == 6)) Call error(540)

    Call invert(config%cell, rcell, det)

    ! Specify molecular dynamics simulation cell

    If (config%newjob_check_config) Then
      Write (message, "('simulation cell vectors')")
      Call Info(message, .true.)
      Write (message, "(3f20.10)") config%cell(1:3)
      Call Info(message, .true.)
      Write (message, "(3f20.10)") config%cell(4:6)
      Call Info(message, .true.)
      Write (message, "(3f20.10)") config%cell(7:9)
      Call Info(message, .true.)
      Write (message, "('system volume     ',2x,1p,g22.12)") det
      Call Info(message, .true.)
    End If

    ! Check on validity of CONFIG contents wrt FIELD

    If (flow%restart_key /= RESTART_KEY_CLEAN .and. config%levcfg < 1) Call error(85)

    If (flow%strict) iwrk(1:config%natms) = 0 ! initialise

    ! Safe flag

    safe = .true.

    ! Global atom counter

    totatm = 0

    ! Local atom counter

    indatm = 1

    ! Site counter for molecules

    mol_sit = 0

    ! Check atom names and assign atomic characteristics
    ! Connecting FIELD to CONFIG on local basis

    Do k = 1, sites%ntype_mol

      Do l = 1, sites%num_mols(k)

        Do m = 1, sites%num_site(k)

          ! Increase global atom counter

          totatm = totatm + 1

          ! If a local atom has a global index totatm

          If (config%lsa(indatm) == totatm) Then

            ! Get the local index. mol_sit+m is the global site

            loc_ind = config%lsi(indatm)

            ! Second check: Do particle identities and their order from the topology
            ! (FIELD) match the one found in the crystallographic data (CONFIG)?
            ! Check for unidentified atoms in CONFIG by their existence in FIELD

            If (config%atmnam(loc_ind) /= sites%site_name(mol_sit + m)) Then
              Write (message, "( 'unidentified atom label :',a8,': atom number ',i5)") config%atmnam(loc_ind), totatm
              Call info(message)
              safe = .false.
            End If

            ! Assign global sites, type, weight, charge & frozen status to localised atoms

            config%lsite(loc_ind) = mol_sit + m
            config%ltype(loc_ind) = sites%type_site(mol_sit + m)
            config%weight(loc_ind) = sites%weight_site(mol_sit + m)
            config%parts(loc_ind)%chge = sites%charge_site(mol_sit + m)
            config%lfrzn(loc_ind) = sites%freeze_site(mol_sit + m)
            config%lfree(loc_ind) = sites%free_site(mol_sit + m)

            ! Print global indices for a later check on ordering (mixed indexing)

            If (flow%strict) iwrk(indatm) = totatm ! Populate

            ! Increase local atom counter

            indatm = indatm + 1

          End If

        End Do

      End Do

      ! Increase site counter per molecule

      mol_sit = mol_sit + sites%num_site(k)

    End Do
    indatm = indatm - 1 ! Correct presence number

    ! Global check on successful connection between CONFIG and FIELD

    Call gcheck(comm, safe)
    If (.not. safe) Call error(25)

    ! Check CONFIG indices eligibility

    If (flow%strict) Then
      Call all_inds_present(iwrk, indatm, config%megatm, safe)
      Call gcheck(comm, safe)
      If (.not. safe) Call error(28)

      Deallocate (iwrk, Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'check_config deallocation failure'
        Call error(0, message)
      End If
    End If

    ! For subsequent checks

    If (config%newjob_check_config) config%newjob_check_config = .false.

  Contains

    Subroutine all_inds_present(ind, n_loc, n, all_present)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! dl_poly_4 routine to check if the array IND which is distributed
      ! accross the communicator contains all indices from 1 to N once,
      ! and only once. N_LOC is the number of indices local to this processor
      !
      ! copyright - daresbury laboratory
      ! author    - i.j.bush march 2009
      ! adapted   - i.t.todorov june 2010
      !
      ! refactoring:
      !           - a.m.elena march-october 2018
      !           - j.madge march-october 2018
      !           - a.b.g.chalk march-october 2018
      !           - i.scivetti march-october 2018
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Dimension(:), Intent(In   ) :: ind
    Integer,               Intent(In   ) :: n_loc, n
    Logical,               Intent(  Out) :: all_present

    Integer                            :: fail, i, me, nproc, where_send
    Integer, Allocatable, Dimension(:) :: all_n_loc, displs_recv, displs_send, loc_start, &
                                          local_ind, reorg_ind, to_recv, to_send
    Logical                            :: loc_present

      me = comm%idnode
      nproc = comm%mxnode

      ! Ever the optimist !
      all_present = .true.

      ! If the number on each proc do not add up to N there is no
      ! way things can be correct. will need all the local totals so
      ! gather and sum locally rather than just a global sum as
      ! saves a comm.
      Allocate (all_n_loc(0:nproc - 1), Stat=fail)
      If (fail /= 0) Go To 100

      ! No check on mpi error as all is pointless because
      ! of the way mpi error handlers are usually dealt with
      !
      !    all_n_loc = 0
      !    all_n_loc( comm%idnode ) = n_loc
      !    Call gsum( all_n_loc( 0:nrpocs - 1 ) )
      !
      Call gallgather(comm, n_loc, all_n_loc(:))
      all_present = (Sum(all_n_loc) == n)
      If (.not. all_present) Return

      ! Work out the first index on each proc
      Allocate (loc_start(0:nproc - 1), Stat=fail)
      If (fail /= 0) Go To 100

      loc_start(0) = 1
      Do i = 1, nproc - 1
        loc_start(i) = loc_start(i - 1) + all_n_loc(i - 1)
      End Do

      ! Array to work with local indices
      Allocate (local_ind(1:n_loc), Stat=fail)
      If (fail /= 0) Go To 100

      local_ind = ind(1:n_loc)

      ! Sort the local data
      Call shellsort(n_loc, local_ind)

      ! Work out how much data to send to each processor
      Allocate (to_send(0:nproc - 1), Stat=fail)
      If (fail /= 0) Go To 100

      to_send = 0
      where_send = 0
      Do i = 1, n_loc
        Do While (local_ind(i) > loc_start(where_send) + all_n_loc(where_send) - 1)
          where_send = where_send + 1
        End Do
        to_send(where_send) = to_send(where_send) + 1
      End Do

      ! How much node i sends to me is how much I recv
      Allocate (to_recv(0:nproc - 1), Stat=fail)
      If (fail /= 0) Go To 100

      Call galltoall(comm, to_send(:), 1, to_recv(:))

      ! Work out the displacements in the sending and receiving arrays
      Allocate (displs_send(0:nproc - 1), Stat=fail)
      If (fail /= 0) Go To 100

      Allocate (displs_recv(0:nproc - 1), Stat=fail)
      If (fail /= 0) Go To 100

      displs_send(0) = 0
      Do i = 1, nproc - 1
        displs_send(i) = displs_send(i - 1) + to_send(i - 1)
      End Do

      displs_recv(0) = 0
      Do i = 1, nproc - 1
        displs_recv(i) = displs_recv(i - 1) + to_recv(i - 1)
      End Do

      ! Put the index on the proc that should own it
      Allocate (reorg_ind(1:n_loc), Stat=fail)
      If (fail /= 0) Go To 100

      Call galltoallv(comm, local_ind(:), to_send(:), displs_send(:), &
                      reorg_ind(:), to_recv(:), displs_recv(:))

      ! Sort the reorganized data
      Call shellsort(n_loc, reorg_ind)

      ! Check it starts and ends at the right place for any owned particle
      ! i.e. some domaines my not own any (domain on vacuum)
      If (n_loc > 0) &
        all_present = (reorg_ind(1) == loc_start(me) .and. &
                       reorg_ind(n_loc) == loc_start(me) + n_loc - 1)

      If (all_present) Then
        ! Check it contains all the numbers in between
        Do i = 2, n_loc
          all_present = (reorg_ind(i) - reorg_ind(i - 1) == 1)
          If (.not. all_present) Then
            Exit
          End If
        End Do
      End If

      ! Is everybody happy?
      loc_present = all_present
      Call gallreduce(comm, loc_present, all_present, op_land)

      Deallocate (reorg_ind, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (displs_recv, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (displs_send, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (to_recv, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (to_send, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (local_ind, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (loc_start, Stat=fail)
      If (fail /= 0) Go To 100
      Deallocate (all_n_loc, Stat=fail)
      If (fail /= 0) Go To 100

      Return

      100 Continue
      If (fail > 0) Then
        Write (message, '(a)') 'all_inds_present allocation/deallocation failure'
        Call error(0, message)
      End If

    End Subroutine all_inds_present

  End Subroutine check_config

  Subroutine read_config(config, megatm, levcfg, l_ind, strict, rcut, dvar, xhi, yhi, zhi, dens0, dens, io, domain, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reading CONFIG and getting the average
    ! particle density
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    ! contrib   - a.m.elena february 2017
    ! contrib   - i.t.todorov may 2019 - maximum domain density change
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(io_type), Intent(InOut) :: io
    Integer, Intent(In) :: megatm, levcfg
    Logical, Intent(In) :: l_ind, strict
    Real(Kind=wp), Intent(In) :: rcut, dvar
    Real(Kind=wp), Intent(InOut) :: xhi, yhi, zhi
    Real(Kind=wp), Intent(Out) :: dens0, dens
    Type(configuration_type), Intent(InOut) :: config
    Type(domains_type), Intent(In) :: domain
    Type(file_type), Intent(InOut) :: files(:)
    Type(comms_type), Intent(InOut) :: comm

    Real(Kind=wp) :: cut

    Character(Len=200) :: record
    Character(Len=40) :: word, fname
    Logical                :: safe, l_his, l_xtr, fast
    Integer                :: fail(1:4), i, j, idm, max_fail, min_fail, &
                              icell, ncells, &
                              indatm, nattot, totatm, &
                              ipx, ipy, ipz, nlx, nly, nlz, &
                              ix, iy, iz, jx, jy, jz
    Real(Kind=wp)      :: celprp(1:10), rcell(1:9), celh(1:9), det, &
                          volm, vcell, &
                          sxx, syy, szz, xdc, ydc, zdc, &
                          pda_max, pda_min, pda_ave, &
                          pda_dom_max, pda_dom_min

    ! Some parameters and variables needed by io interfaces

    Integer                       :: fh, io_read
    Integer(Kind=offset_kind) :: top_skip

    Real(Kind=wp), Dimension(:), Allocatable :: pda

    Character(Len=8), Dimension(:), Allocatable :: chbuf
    Integer, Dimension(:), Allocatable :: iwrk
    Real(Kind=wp), Dimension(:), Allocatable :: axx, ayy, azz, &
                                                bxx, byy, bzz, &
                                                cxx, cyy, czz

    Character(Len=256) :: message
    Character(Len=256) :: messages(3)

    safe = .true. ! we start safe
    l_his = .false. ! no HISTORY reading
    l_xtr = .false. ! no CONFIG cell extremes needed

    ! image conditions not compliant with DD and link-cell

    If (config%imcon == 4 .or. config%imcon == 5 .or. config%imcon == 7) Call error(300)

    ! Real space cutoff shortened by 50% but not < 1 Angstrom
    !(or ==rcut_def in scan_control)

    cut = Max(0.5_wp * rcut, 1.0_wp) + 1.0e-6_wp

    ! Get the dimensional properties of the MD cell

    Call dcell(config%cell, celprp)
    volm = celprp(10)

    ! Calculate the number of link-cells per domain in every direction

    nlx = Int(celprp(7) / (cut * domain%nx_real))
    nly = Int(celprp(8) / (cut * domain%ny_real))
    nlz = Int(celprp(9) / (cut * domain%nz_real))

    ncells = nlx * nly * nlz

    ! Check for link cell algorithm violations

    If (ncells == 0) Call error(307)

    ! Amend volume of density cell if cluster, slab or bulk slab
    ! cell dimensional properties overwritten but not needed anyway

    If (config%imcon == 0 .or. config%imcon == 6 .or. config%imc_n == 6) Then
      celh = config%cell

      If (config%imcon == 0) Then
        celh(1) = Max(1.0_wp, xhi)
        celh(5) = Max(1.0_wp, yhi)
        celh(9) = Max(1.0_wp, zhi)
      Else If (config%imcon == 6) Then
        celh(9) = Max(1.0_wp, zhi)
      End If

      Call dcell(celh, celprp)
      volm = celprp(10)
    End If

    vcell = volm / (Real(ncells, wp) * Real(comm%mxnode, wp))

    ! Approximate density and mxatms

    dens = Real(megatm, wp) / volm
    config%mxatms = 10 * Max(1, Nint((dvar**1.7_wp) * dens * vcell * Real((nlx + 3) * (nly + 3) * (nlz + 3), wp)))

    ! Allocate necessary arrays to read CONFIG

    Call config%init_read()

    ! Get type of I/O for reading

    Call io_get_parameters(io, user_method_read=io_read)

    ! Define filename ASCII or netCDF

    If (io_read /= IO_READ_NETCDF) Then
      fname = Trim(files(FILE_CONFIG)%filename)
    Else
      fname = Trim(files(FILE_CONFIG)%filename)//'.nc'
    End If

    ! Define/Detect the FAST reading status

    If (io_read == IO_READ_MASTER) Then

      fast = .false.

    Else If (io_read == IO_READ_NETCDF) Then

      fast = .true.

    Else

      ! Check if the system input file is a new style CONFIG:
      ! (i)  all lines are 72 ASCII characters long with
      !      a UNIX carriage return as end of line;
      ! (ii) LINE2 has the particles total value
      !      after values of levcfg and imcon.
      ! No fall back if users have mangled with further lines

      fast = .true.
      If (comm%idnode == 0) Then

        ! Open CONFIG

        Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)

        ! Read the CONFIG file header (TITLE record)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
          i = i + 1
          safe = .false.
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a1)', Advance='No', IOStat=j, End=10) record(i:i)
          safe = .true.
          If (j < 0) Go To 10
        End Do
        10 Continue
        fast = (fast .and. i == recsz)

        ! Read configuration level and image condition (RECORD2)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
          i = i + 1
          safe = .false.
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a1)', Advance='No', IOStat=j, End=20) record(i:i)
          safe = .true.
          If (j < 0) Go To 20
        End Do
        20 Continue
        fast = (fast .and. i == recsz)

        ! Read particles total value

        Call get_word(record, word); Call get_word(record, word)
        Call get_word(record, word); i = Nint(word_2_real(word, 0.0_wp, strict))
        fast = (fast .and. i == megatm)

      End If
      Call gsync(comm)
      Call gcheck(comm, safe, "enforce")
      Call gcheck(comm, fast, "enforce")

      If (.not. safe) Go To 50

      ! Close CONFIG

      If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()

    End If

    fail = 0

    ! If MASTER read

    If (io_read == IO_READ_MASTER) Then

      Call invert(config%cell, rcell, det)

      ! Open CONFIG and skip the header

      If (comm%idnode == 0) Then
        Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)

        Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*) ! CONFIG file header (TITLE record)
        Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*) ! configuration level and image condition

        If (config%imcon /= 0) Then
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*) ! cell vectors (not defined for imcon=0) but cell
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*) ! is modified in set_bounds for imcon 0 and 6!!!
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*)
        End If
      End If

      Allocate (chbuf(1:config%mxatms), iwrk(1:config%mxatms), Stat=fail(1))
      Allocate (axx(1:config%mxatms), ayy(1:config%mxatms), azz(1:config%mxatms), Stat=fail(2))
      Allocate (bxx(1:config%mxatms), byy(1:config%mxatms), bzz(1:config%mxatms), Stat=fail(3))
      Allocate (cxx(1:config%mxatms), cyy(1:config%mxatms), czz(1:config%mxatms), Stat=fail(4))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'read_config allocation failure'
        Call error(0, message)
      End If

      ! Initialise domain localised atom counter (configuration)
      ! and dispatched atom counter

      config%natms = 0
      indatm = 0
      Do nattot = 1, megatm
        indatm = indatm + 1

        ! Initialise transmission arrays

        chbuf(indatm) = ' '
        iwrk(indatm) = 0

        axx(indatm) = 0.0_wp
        ayy(indatm) = 0.0_wp
        azz(indatm) = 0.0_wp

        If (levcfg > 0) Then
          bxx(indatm) = 0.0_wp
          byy(indatm) = 0.0_wp
          bzz(indatm) = 0.0_wp

          If (levcfg > 1) Then
            cxx(indatm) = 0.0_wp
            cyy(indatm) = 0.0_wp
            czz(indatm) = 0.0_wp
          End If
        End If

        ! Read in transmission arrays

        If (comm%idnode == 0 .and. safe) Then
          record = ' '
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a)', End=30) record
          Call tabs_2_blanks(record); Call strip_blanks(record)
          Call get_word(record, word); chbuf(indatm) = word(1:8)
          If (l_ind) Then
            Call get_word(record, word)
            iwrk(indatm) = Nint(word_2_real(word, 0.0_wp, strict))
            If (iwrk(indatm) /= 0) Then
              iwrk(indatm) = Abs(iwrk(indatm))
            Else
              iwrk(indatm) = nattot
            End If
          Else
            iwrk(indatm) = nattot
          End If

          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) axx(indatm), ayy(indatm), azz(indatm)

          If (levcfg > 0) Then
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) bxx(indatm), byy(indatm), bzz(indatm)
            If (levcfg > 1) Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) cxx(indatm), cyy(indatm), czz(indatm)
          End If
          Go To 40

          30 Continue
          safe = .false. ! catch error

          40 Continue
        End If

        ! Circulate configuration data to all nodes when transmission arrays are filled up

        If (indatm == config%mxatms .or. nattot == megatm) Then

          ! Check if batch was read fine

          Call gcheck(comm, safe)
          If (.not. safe) Go To 50

          ! Ensure all atoms are in prescribed simulation cell (DD bound) and broadcast them
          !
          !           Call pbcshift(imcon,cell,indatm,axx,ayy,azz)

          Call gbcast(comm, chbuf, 0)
          Call gbcast(comm, iwrk, 0)

          Call gbcast(comm, axx, 0)
          Call gbcast(comm, ayy, 0)
          Call gbcast(comm, azz, 0)

          If (levcfg > 0) Then

            Call gbcast(comm, bxx, 0)
            Call gbcast(comm, byy, 0)
            Call gbcast(comm, bzz, 0)

            If (levcfg > 1) Then

              Call gbcast(comm, cxx, 0)
              Call gbcast(comm, cyy, 0)
              Call gbcast(comm, czz, 0)
            End If
          End If

          ! Assign atoms to correct domains

          Do i = 1, indatm
            sxx = rcell(1) * axx(i) + rcell(4) * ayy(i) + rcell(7) * azz(i)
            syy = rcell(2) * axx(i) + rcell(5) * ayy(i) + rcell(8) * azz(i)
            szz = rcell(3) * axx(i) + rcell(6) * ayy(i) + rcell(9) * azz(i)

            ! sxx,syy,szz are in [-0.5,0.5) interval as values as 0.4(9) may pose a problem

            sxx = sxx - Anint(sxx); If (sxx >= half_minus) sxx = -sxx
            syy = syy - Anint(syy); If (syy >= half_minus) syy = -syy
            szz = szz - Anint(szz); If (szz >= half_minus) szz = -szz

            ! fold back coordinates

            axx(i) = config%cell(1) * sxx + config%cell(4) * syy + config%cell(7) * szz
            ayy(i) = config%cell(2) * sxx + config%cell(5) * syy + config%cell(8) * szz
            azz(i) = config%cell(3) * sxx + config%cell(6) * syy + config%cell(9) * szz

            ! assign domain coordinates (call for errors)

            ipx = Int((sxx + 0.5_wp) * domain%nx_real)
            ipy = Int((syy + 0.5_wp) * domain%ny_real)
            ipz = Int((szz + 0.5_wp) * domain%nz_real)

            idm = ipx + domain%nx * (ipy + domain%ny * ipz)
            If (idm < 0 .or. idm > (comm%mxnode - 1)) Then
              Call error(513)
            Else If (idm == comm%idnode) Then
              config%natms = config%natms + 1

              If (config%natms < config%mxatms) Then
                config%atmnam(config%natms) = chbuf(i)
                config%ltg(config%natms) = iwrk(i)

                config%parts(config%natms)%xxx = axx(i)
                config%parts(config%natms)%yyy = ayy(i)
                config%parts(config%natms)%zzz = azz(i)

                If (levcfg > 0) Then
                  config%vxx(config%natms) = bxx(i)
                  config%vyy(config%natms) = byy(i)
                  config%vzz(config%natms) = bzz(i)
                Else
                  config%vxx(config%natms) = 0.0_wp
                  config%vyy(config%natms) = 0.0_wp
                  config%vzz(config%natms) = 0.0_wp
                End If

                If (levcfg > 1) Then
                  config%parts(config%natms)%fxx = cxx(i)
                  config%parts(config%natms)%fyy = cyy(i)
                  config%parts(config%natms)%fzz = czz(i)
                Else
                  config%parts(config%natms)%fxx = 0.0_wp
                  config%parts(config%natms)%fyy = 0.0_wp
                  config%parts(config%natms)%fzz = 0.0_wp
                End If
              Else
                safe = .false.
              End If
            End If
          End Do

          ! Check if all is dispatched fine

          max_fail = config%natms
          min_fail = config%natms
          Call gcheck(comm, safe)
          Call gmax(comm, max_fail)
          Call gmin(comm, min_fail)

          If (.not. safe) Then
            Write (messages(1), '(a,i0)') 'next error due to maximum number of atoms per domain set to : ', config%mxatms
            Write (messages(2), '(2(a,i0))') 'but maximum & minumum numbers of atoms per domain asked for : ', &
              max_fail, ' & ', min_fail
            Write (messages(3), '(a,i0)') 'estimated densvar value for passing this stage safely is : ', &
              Ceiling((dvar * (Real(max_fail, wp) / Real(config%mxatms, wp))**(1.0_wp / 1.7_wp) - 1.0_wp) * 100.0_wp)
            Call info(messages, 3, .true.)
            Call error(45)
          End If

          ! Nullify dispatch counter

          indatm = 0

        End If
      End Do

      ! Close CONFIG

      If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()
      Call gsync(comm)

      Deallocate (chbuf, iwrk, Stat=fail(1))
      Deallocate (axx, ayy, azz, Stat=fail(2))
      Deallocate (bxx, byy, bzz, Stat=fail(3))
      Deallocate (cxx, cyy, czz, Stat=fail(4))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'read_config deallocation failure'
        Call error(0, message)
      End If

      ! If PROPER read

    Else

      ! Open CONFIG

      If (fast) Then
        Call io_set_parameters(io, user_comm=comm%comm)
        Call io_init(io, recsz)
        Call io_open(io, io_read, comm%comm, fname, mode_rdonly, fh)
      Else
        Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)
      End If

      ! top_skip is header size

      If (io_read /= IO_READ_NETCDF) Then
        If (config%imcon == 0) Then
          top_skip = Int(2, offset_kind)
        Else
          top_skip = Int(5, offset_kind)
        End If
      Else
        top_skip = Int(1, offset_kind) ! This is now the frame = 1
      End If

      Call read_config_parallel(config, levcfg, dvar, l_ind, strict, megatm, l_his, l_xtr, &
                                fast, fh, top_skip, xhi, yhi, zhi, io, domain, files, comm)

      ! Close CONFIG

      If (fast) Then
        Call io_close(io, fh)
        Call io_finalize(io)
      Else
        Call files(FILE_CONFIG)%close ()
      End If

    End If

    ! To prevent users from the danger of changing the order of calls
    ! in dl_poly set 'nlast' to the innocent 'natms'

    config%nlast = config%natms

    ! Check if the number of atoms in the system (MD cell) derived by
    ! topology description (FIELD) match the crystallographic (CONFIG) one?

    totatm = config%natms
    Call gsum(comm, totatm)
    If (totatm /= megatm) Call error(58)

    ! Record global atom indices for local sorting (configuration)

    Do i = 1, config%natms
      config%lsi(i) = i
      config%lsa(i) = config%ltg(i)
    End Do
    Call shellsort2(config%natms, config%lsi, config%lsa)

    If (io_read /= IO_READ_MASTER) Then

      ! This section is not strictly necessary.  However, the new read in method
      ! means the atoms are not necessarily in the same order in memory as the
      ! older, slower, method would put them.  This bit makes sure that the order
      ! is so that 'ltg' is strictly monotonically increasing.  This captures the
      ! common case where CONFIG has the 'ltg' values all in order (or not
      ! specified), but there is no easy way for the general case of arbitrary
      ! ordering of the 'ltg' values in CONFIG.  Of course, this makes no
      ! difference to the science and to restarts.  However, for initial runs it
      ! means the initial velocities will not be the same as the old method for
      ! the arbitrary ordering case.

      config%atmnam(1:config%natms) = config%atmnam(config%lsi(1:config%natms))
      config%ltg(1:config%natms) = config%ltg(config%lsi(1:config%natms))

      config%parts(1:config%natms)%xxx = config%parts(config%lsi(1:config%natms))%xxx
      config%parts(1:config%natms)%yyy = config%parts(config%lsi(1:config%natms))%yyy
      config%parts(1:config%natms)%zzz = config%parts(config%lsi(1:config%natms))%zzz

      If (levcfg > 0) Then
        config%vxx(1:config%natms) = config%vxx(config%lsi(1:config%natms))
        config%vyy(1:config%natms) = config%vyy(config%lsi(1:config%natms))
        config%vzz(1:config%natms) = config%vzz(config%lsi(1:config%natms))

        If (levcfg > 1) Then
          config%parts(1:config%natms)%fxx = config%parts(config%lsi(1:config%natms))%fxx
          config%parts(1:config%natms)%fyy = config%parts(config%lsi(1:config%natms))%fyy
          config%parts(1:config%natms)%fzz = config%parts(config%lsi(1:config%natms))%fzz
        End If
      End If
      Do i = 1, config%natms
        config%lsi(i) = i
        config%lsa(i) = config%ltg(i)
      End Do

    End If

    ! READ CONFIG END

    ! PARTICLE DENSITY START
    ! Allocate and initialise particle density array

    Allocate (pda(1:ncells), Stat=fail(1))
    If (fail(1) > 0) Then
      Write (message, '(a)') 'read_config allocation failure'
      Call error(0, message)
    End If
    pda = 0.0_wp

    ! Get the total number of link-cells in MD cell per direction

    xdc = Real(nlx * domain%nx, wp)
    ydc = Real(nly * domain%ny, wp)
    zdc = Real(nlz * domain%nz, wp)

    ! Shifts from global to local link-cell space:
    ! (0,0,0) left-most link-cell on the domain (halo)
    ! (nlx+2*nlp-1,nly+2*nlp-1,nly+2*nlp-1) right-most
    ! link-cell on the domain (halo)

    jx = 1 - nlx * domain%idx
    jy = 1 - nly * domain%idy
    jz = 1 - nlz * domain%idz

    ! Get the inverse cell matrix

    Call invert(config%cell, rcell, celprp(10))

    Do i = 1, config%natms
      sxx = rcell(1) * config%parts(i)%xxx + rcell(4) * config%parts(i)%yyy + rcell(7) * config%parts(i)%zzz
      syy = rcell(2) * config%parts(i)%xxx + rcell(5) * config%parts(i)%yyy + rcell(8) * config%parts(i)%zzz
      szz = rcell(3) * config%parts(i)%xxx + rcell(6) * config%parts(i)%yyy + rcell(9) * config%parts(i)%zzz

      ! Get cell coordinates accordingly

      ix = Int(xdc * (sxx + 0.5_wp)) + jx
      iy = Int(ydc * (syy + 0.5_wp)) + jy
      iz = Int(zdc * (szz + 0.5_wp)) + jz

      ! Put all particles in bounded link-cell space: lower and upper
      ! bounds as 1 <= i_coordinate <= nl_coordinate

      ix = Max(Min(ix, nlx), 1)
      iy = Max(Min(iy, nly), 1)
      iz = Max(Min(iz, nlz), 1)

      ! Hypercube function transformation (counting starts from one
      ! rather than zero /map_domains/

      icell = 1 + (ix - 1) + nlx * ((iy - 1) + nly * (iz - 1))

      pda(icell) = pda(icell) + 1.0_wp
    End Do

    pda_max = 0.0_wp
    pda_min = 100.0_wp
    pda_ave = 0.0_wp
    Do icell = 1, ncells
      pda_max = Max(pda(icell), pda_max)
      pda_min = Min(pda(icell), pda_min)
      pda_ave = pda_ave + pda(icell)
    End Do
    pda_ave = pda_ave / Real(ncells, wp)

    pda_dom_max = pda_ave
    pda_dom_min = pda_ave
    Call gmax(comm, pda_dom_max)
    Call gmin(comm, pda_dom_min)

    Call gmax(comm, pda_max)
    Call gmin(comm, pda_min)

    Call gsum(comm, pda_ave)
    pda_ave = pda_ave / Real(comm%mxnode, wp)

    ! define dens0 & dens

    dens0 = pda_max / vcell ! maximum local density
    dens = pda_dom_max / vcell ! maximum domain density

    Deallocate (pda, Stat=fail(1))
    If (fail(1) > 0) Then
      Write (message, '(a)') 'read_config deallocation failure'
      Call error(0, message)
    End If

    ! PARTICLE DENSITY END

    Return

    ! error exit for CONFIG file read

    50 Continue
    If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()
    Call error(55)

  End Subroutine read_config

  Subroutine read_config_parallel(config, levcfg, dvar, l_ind, strict, megatm, l_his, &
                                  l_xtr, fast, fh, top_skip, xhi, yhi, zhi, io, domain, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reading in the CONFIG data file in parallel
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush & i.t.todorov march 2016
    ! amended   - i.t.todorov february 2018 - non-fast read fix
    ! amended   - i.t.todorov september 2018 - record initialisation
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type),  Intent(InOut) :: config
    Integer,                   Intent(In   ) :: levcfg
    Real(Kind=wp),             Intent(In   ) :: dvar
    Logical,                   Intent(In   ) :: l_ind, strict
    Integer,                   Intent(In   ) :: megatm
    Logical,                   Intent(In   ) :: l_his, l_xtr, fast
    Integer,                   Intent(In   ) :: fh
    Integer(Kind=offset_kind), Intent(In   ) :: top_skip
    Real(Kind=wp),             Intent(  Out) :: xhi, yhi, zhi
    Type(io_type),             Intent(InOut) :: io
    Type(domains_type),        Intent(In   ) :: domain
    Type(file_type),           Intent(InOut) :: files(:)
    Type(comms_type),          Intent(InOut) :: comm

    Character(Len=1), Allocatable, Dimension(:, :) :: rec_buff
    Character(Len=200)                             :: record
    Character(Len=256)                             :: message, messages(3)
    Character(Len=40)                              :: forma, word
    Character(Len=8), Allocatable, Dimension(:)    :: chbuf, chbuf_read, chbuf_scat
    Integer                                        :: ats_per_proc, batsz, Count(1:3), fail(1:8), &
                                                      frame, i, idm, ierr, indatm, io_read, ipx, &
                                                      ipy, ipz, j, k, max_fail, min_fail, &
                                                      my_read_proc_num, n_ats_in_file, n_loc, &
                                                      n_read_procs_use, per_read_proc, &
                                                      recs_per_at, recs_per_proc, recs_to_read, &
                                                      recsz, start(1:3), this_base_proc, &
                                                      this_rec_buff, to_read, which_read_proc, &
                                                      wp_vals_per_at
    Integer(Kind=li)                               :: n_ii, n_jj, n_sk
    Integer(Kind=offset_kind)                      :: n_skip, rec_mpi_io
    Integer, Allocatable, Dimension(:)             :: first_at, iwrk, iwrk_read, iwrk_scat, &
                                                      n_held, orig_first_at, owner_read, where_buff
    Logical                                        :: do_read, safe
    Real(Kind=wp)                                  :: det, rcell(1:9), sxx, syy, szz
    Real(Kind=wp), Allocatable, Dimension(:)       :: axx_read, ayy_read, azz_read, bxx_read, &
                                                      byy_read, bzz_read, cxx_read, cyy_read, &
                                                      czz_read
    Real(Kind=wp), Allocatable, Dimension(:, :)    :: scatter_buffer, scatter_buffer_read

! Some parameters and variables needed by io interfaces
! netCDF

    ! Get reading method, total number of I/O heads and buffer size

    Call io_get_parameters(io, user_method_read=io_read)
    Call io_get_parameters(io, user_n_io_procs_read=n_read_procs_use)
    Call io_get_parameters(io, user_buffer_size_read=batsz)

    fail = 0 ! fail initialisation

    If (levcfg /= 3) Then
      wp_vals_per_at = 3 * (1 + levcfg) ! Scatter buffer sizes
      recs_per_at = 2 + levcfg ! Scatter buffer sizes
    Else
      wp_vals_per_at = 3 ! Scatter buffer sizes
      recs_per_at = 1 ! Scatter buffer sizes
    End If

    ! Note: make 'first_at' and 'orig_first_at' 1 element bigger than strictly
    ! required to make checking at the end of reading much easier and clearer

    Allocate (first_at(0:n_read_procs_use), orig_first_at(0:n_read_procs_use), Stat=fail(1))
    Allocate (chbuf(1:batsz), iwrk(1:batsz), Stat=fail(2))
    Allocate (scatter_buffer(1:wp_vals_per_at, 1:batsz), Stat=fail(3))
    If (Any(fail(1:3) > 0)) Then
      Write (message, '(a)') 'read_config_parallel allocation failure 1'
      Call error(0, message)
    End If

    ! define basic quantities for the parallel ASCII reading

    per_read_proc = comm%mxnode / n_read_procs_use
    do_read = (Mod(comm%idnode, per_read_proc) == 0 .and. comm%idnode < per_read_proc * n_read_procs_use)
    my_read_proc_num = comm%idnode / per_read_proc

    ! Note 'first_at' and 'orig_first_at' have one more element
    ! in the array than strictly required - makes it easier to
    ! check that reading by the last I/O processor has finished

    ats_per_proc = megatm / n_read_procs_use
    Do i = 0, n_read_procs_use
      first_at(i) = i * ats_per_proc + Min(i, megatm - ats_per_proc * n_read_procs_use)
    End Do
    orig_first_at = first_at
    ats_per_proc = Max(1, ats_per_proc) ! Fix it if 0
    recs_per_proc = ats_per_proc * recs_per_at

    ! Catch the case where the first atom belonging to
    ! a read processor does not actually exists - i.e.
    ! I/O procs count > megatm, and limit reading by do_read

    If (my_read_proc_num < n_read_procs_use) &
      do_read = (do_read .and. first_at(my_read_proc_num) < megatm)

    ! Skip to the point of reading

    If (do_read) Then

      n_skip = Int(recs_per_at, offset_kind) * Int(first_at(my_read_proc_num), offset_kind) + &
               top_skip - Int(1, offset_kind)
      If (.not. fast) Then
        n_sk = Int(n_skip, li)
        n_jj = 73 * batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
          Do n_ii = 1_li, n_sk / n_jj
            forma = ' '
            Write (forma, '( "(", i0, "/)" )') n_jj
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, End=100)
          End Do
          n_ii = Mod(n_sk, n_jj) - n_ii + 1_li
          If (n_ii > 0_li) Then
            forma = ' '
            Write (forma, '( "(", i0, "/)" )') n_ii
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, End=100)
          End If
        Else
          forma = ' '
          Write (forma, '( "(", i0, "/)" )') n_sk
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, End=100)
        End If

        recsz = 200
        forma = ' '
        Write (forma, '( "(", i0, "a1)" )') recsz
      Else
        rec_mpi_io = n_skip + Int(1, offset_kind)
        recsz = 73
        If (levcfg == 3) recsz = 35
      End If

      ! Allocate record buffer, reading buffers, scatter buffers and indexing arrays

      If (io_read /= IO_READ_NETCDF) Then
        Allocate (rec_buff(1:recsz, 1:batsz), Stat=fail(1))
      Else
        Allocate (rec_buff(1:Len(chbuf_read), 1:batsz), Stat=fail(1))
      End If
      Allocate (chbuf_read(1:batsz), iwrk_read(1:batsz), Stat=fail(2))
      Allocate (axx_read(1:batsz), ayy_read(1:batsz), azz_read(1:batsz), Stat=fail(3))
      Allocate (bxx_read(1:batsz), byy_read(1:batsz), bzz_read(1:batsz), Stat=fail(4))
      Allocate (cxx_read(1:batsz), cyy_read(1:batsz), czz_read(1:batsz), Stat=fail(5))
      Allocate (scatter_buffer_read(1:wp_vals_per_at, 1:batsz), Stat=fail(6))
      Allocate (chbuf_scat(1:batsz), iwrk_scat(1:batsz), Stat=fail(7))
      Allocate (n_held(0:comm%mxnode - 1), where_buff(0:comm%mxnode - 1), owner_read(1:batsz), Stat=fail(8))
      If (Any(fail(1:8) > 0)) Then
        Write (message, '(a)') 'read_config_parallel allocation failure 2'
        Call error(0, message)
      End If

    Else

      ! It is Illegal to pass unallocated allocatable arrays to routines.
      ! Therefore for arrays that are used by the mpi_scatterv calls
      ! below allocate them to zero size if they are not used on this core

      Allocate (scatter_buffer_read(1:0, 1:0), Stat=fail(1))
      Allocate (chbuf_scat(1:0), iwrk_scat(1:0), Stat=fail(2))
      Allocate (n_held(0:-1), where_buff(0:-1), Stat=fail(3))
      If (Any(fail(1:3) > 0)) Then
        Write (message, '(a)') 'read_config_parallel allocation failure 3'
        Call error(0, message)
      End If

    End If

    ! Initialise extreme box dimensions

    xhi = 0.0_wp
    yhi = 0.0_wp
    zhi = 0.0_wp

    If (.not. l_xtr) Call invert(config%cell, rcell, det)

    ! Initialise domain localised atom counter (configuration),
    ! dispatched atom counter and safe dispatch flag

    config%natms = 0
    indatm = 0
    safe = .true.

    Do k = 1, megatm

      ! Read in transmission arrays

      Readers_only: If (do_read .and. indatm == 0) Then
        to_read = Min(batsz, orig_first_at(my_read_proc_num + 1) - first_at(my_read_proc_num))

        No_netCDF: If (io_read /= IO_READ_NETCDF) Then

          this_rec_buff = 0
          recs_to_read = 0
          Do i = 1, to_read
            If (this_rec_buff == 0) Then
              recs_to_read = Min(Size(rec_buff, Dim=2), (to_read - i + 1) * recs_per_at)
              If (.not. fast) Then
                Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma) rec_buff(:, 1:recs_to_read)
              Else
                Call io_read_batch(io, fh, rec_mpi_io, recs_to_read, rec_buff, ierr)
                rec_mpi_io = rec_mpi_io + Int(recs_to_read, offset_kind)
              End If
            End If

            ! Atom details
            record = ''
            this_rec_buff = this_rec_buff + 1
            Do j = 1, Min(Len(record), Size(rec_buff, Dim=1))
              record(j:j) = rec_buff(j, this_rec_buff)
            End Do
            Call strip_blanks(record)

            Call get_word(record, word); chbuf_read(i) = word(1:8)
            If (l_ind) Then
              Call get_word(record, word)
              iwrk_read(i) = Nint(word_2_real(word, 0.0_wp, strict))
              If (iwrk_read(i) /= 0) Then
                iwrk_read(i) = Abs(iwrk_read(i))
              Else
                iwrk_read(i) = first_at(my_read_proc_num) + i
              End If
            Else
              iwrk_read(i) = first_at(my_read_proc_num) + i
            End If

            If (levcfg == 3) Read (record, Fmt=*, End=100) axx_read(i), ayy_read(i), azz_read(i)

            If (this_rec_buff == recs_to_read) Then
              this_rec_buff = 0
              recs_to_read = Min(Size(rec_buff, Dim=2), (to_read - i + 1) * recs_per_at - 1)
              If (.not. fast) Then
                Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma) rec_buff(:, 1:recs_to_read)
              Else
                Call io_read_batch(io, fh, rec_mpi_io, recs_to_read, rec_buff, ierr)
                rec_mpi_io = rec_mpi_io + Int(recs_to_read, offset_kind)
              End If
            End If

            If (levcfg /= 3) Then

              ! Positions

              this_rec_buff = this_rec_buff + 1
              Do j = 1, Min(Len(record), Size(rec_buff, Dim=1))
                record(j:j) = rec_buff(j, this_rec_buff)
              End Do
              Read (record, Fmt=*, End=100) axx_read(i), ayy_read(i), azz_read(i)
              If (this_rec_buff == recs_to_read) Then
                this_rec_buff = 0
                If (levcfg > 0) Then
                  recs_to_read = Min(Size(rec_buff, Dim=2), (to_read - i + 1) * recs_per_at - 2)
                  If (.not. fast) Then
                    Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma) rec_buff(:, 1:recs_to_read)
                  Else
                    Call io_read_batch(io, fh, rec_mpi_io, recs_to_read, rec_buff, ierr)
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read, offset_kind)
                  End If
                End If
              End If

              ! Velocities

              If (levcfg > 0) Then
                this_rec_buff = this_rec_buff + 1
                Do j = 1, Min(Len(record), Size(rec_buff, Dim=1))
                  record(j:j) = rec_buff(j, this_rec_buff)
                End Do
                Read (record, Fmt=*, End=100) bxx_read(i), byy_read(i), bzz_read(i)
                If (this_rec_buff == recs_to_read) Then
                  this_rec_buff = 0
                  If (levcfg > 1) Then
                    recs_to_read = Min(Size(rec_buff, Dim=2), (to_read - i + 1) * recs_per_at - 3)
                    If (.not. fast) Then
                      Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma) rec_buff(:, 1:recs_to_read)
                    Else
                      Call io_read_batch(io, fh, rec_mpi_io, recs_to_read, rec_buff, ierr)
                      rec_mpi_io = rec_mpi_io + Int(recs_to_read, offset_kind)
                    End If
                  End If
                End If

                ! Forces

                If (levcfg > 1) Then
                  this_rec_buff = this_rec_buff + 1
                  Do j = 1, Min(Len(record), Size(rec_buff, Dim=1))
                    record(j:j) = rec_buff(j, this_rec_buff)
                  End Do
                  Read (record, Fmt=*, End=100) cxx_read(i), cyy_read(i), czz_read(i)
                  If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
                  End If
                End If
              End If

            End If
          End Do

        Else

          If (to_read /= 0) Then
            frame = Int(top_skip, Kind(frame))

            Call io_nc_get_var(io, 'atomnames', fh, rec_buff, [first_at(my_read_proc_num) + 1, frame], [8, to_read, 1])
            Do i = 1, to_read
              Do j = 1, Min(Len(chbuf_read), Size(rec_buff, Dim=1))
                chbuf_read(i) (j:j) = rec_buff(j, i)
              End Do
            End Do
            If (l_ind) Then
              Call io_nc_get_var(io, 'indices', fh, iwrk_read, [first_at(my_read_proc_num) + 1, frame], [to_read, 1])
            End If

            start = [1, first_at(my_read_proc_num) + 1, frame]
            count = [3, to_read, 1]

            Select Case (levcfg)
            Case (0, 3)
              Call io_get_var(io, 'coordinates', fh, start, count, axx_read, ayy_read, azz_read)
            Case (1)
              Call io_get_var(io, 'coordinates', fh, start, count, axx_read, ayy_read, azz_read)
              Call io_get_var(io, 'velocities', fh, start, count, bxx_read, byy_read, bzz_read)
            Case (2)
              Call io_get_var(io, 'coordinates', fh, start, count, axx_read, ayy_read, azz_read)
              Call io_get_var(io, 'velocities', fh, start, count, bxx_read, byy_read, bzz_read)
              Call io_get_var(io, 'forces', fh, start, count, cxx_read, cyy_read, czz_read)
            End Select
          End If

        End If No_netCDF

        If (.not. l_xtr) Then

          ! Ensure all atoms are in prescribed simulation cell (DD bound)
          !
          n_held = 0
          Do i = 1, to_read
            sxx = rcell(1) * axx_read(i) + rcell(4) * ayy_read(i) + rcell(7) * azz_read(i)
            syy = rcell(2) * axx_read(i) + rcell(5) * ayy_read(i) + rcell(8) * azz_read(i)
            szz = rcell(3) * axx_read(i) + rcell(6) * ayy_read(i) + rcell(9) * azz_read(i)

            ! sxx,syy,szz are in [-0.5,0.5) interval as values as 0.4(9) may pose a problem

            sxx = sxx - Anint(sxx); If (sxx >= half_minus) sxx = -sxx
            syy = syy - Anint(syy); If (syy >= half_minus) syy = -syy
            szz = szz - Anint(szz); If (szz >= half_minus) szz = -szz

            ! fold back coordinates

            axx_read(i) = config%cell(1) * sxx + config%cell(4) * syy + config%cell(7) * szz
            ayy_read(i) = config%cell(2) * sxx + config%cell(5) * syy + config%cell(8) * szz
            azz_read(i) = config%cell(3) * sxx + config%cell(6) * syy + config%cell(9) * szz

            ! assign domain coordinates (call for errors)

            ipx = Int((sxx + 0.5_wp) * domain%nx_real)
            ipy = Int((syy + 0.5_wp) * domain%ny_real)
            ipz = Int((szz + 0.5_wp) * domain%nz_real)

            idm = ipx + domain%nx * (ipy + domain%ny * ipz)
            If (idm < 0 .or. idm > (comm%mxnode - 1)) Call error(513)
            owner_read(i) = idm
            n_held(idm) = n_held(idm) + 1
          End Do

          where_buff(0) = 0
          Do i = 1, comm%mxnode - 1
            where_buff(i) = where_buff(i - 1) + n_held(i - 1)
          End Do

          Do i = 1, to_read
            idm = where_buff(owner_read(i))
            idm = idm + 1
            where_buff(owner_read(i)) = idm

            chbuf_scat(idm) = chbuf_read(i)
            iwrk_scat(idm) = iwrk_read(i)

            scatter_buffer_read(1, idm) = axx_read(i)
            scatter_buffer_read(2, idm) = ayy_read(i)
            scatter_buffer_read(3, idm) = azz_read(i)

            If (levcfg /= 3) Then
              If (levcfg > 0) Then
                scatter_buffer_read(4, idm) = bxx_read(i)
                scatter_buffer_read(5, idm) = byy_read(i)
                scatter_buffer_read(6, idm) = bzz_read(i)

                If (levcfg > 1) Then
                  scatter_buffer_read(7, idm) = cxx_read(i)
                  scatter_buffer_read(8, idm) = cyy_read(i)
                  scatter_buffer_read(9, idm) = czz_read(i)
                End If
              End If
            End If
          End Do

          ! If only detecting box dimensions for imcon == 0 or 6 or imc_n == 6

        Else

          ! Get extremes

          xhi = Max(xhi, Maxval(Abs(axx_read(1:to_read))))
          yhi = Max(yhi, Maxval(Abs(ayy_read(1:to_read))))
          zhi = Max(zhi, Maxval(Abs(azz_read(1:to_read))))

        End If
      End If Readers_only

      ! Increase buffer counter and update first_at for
      ! the readers that have something left to read

      indatm = indatm + 1
      If (do_read) Then
        If (first_at(my_read_proc_num) < first_at(my_read_proc_num + 1)) &
          first_at(my_read_proc_num) = first_at(my_read_proc_num) + 1
      End If

      ! Circulate configuration data to all nodes when transmission arrays are filled up
      ! Check against megatm since at low processors counts (i.e. 1) batsz can be > megatm

      Reorganize_buffer: If (indatm == batsz .or. (indatm > 0 .and. k == megatm)) Then

        Extent_2: If (.not. l_xtr) Then

          Do which_read_proc = 0, n_read_procs_use - 1
            If (orig_first_at(which_read_proc) >= megatm) Exit ! for non-reading readers

            this_base_proc = which_read_proc * per_read_proc
            If (comm%idnode == this_base_proc) Then
              where_buff(0) = 0
              Do i = 1, comm%mxnode - 1
                where_buff(i) = where_buff(i - 1) + n_held(i - 1)
              End Do
            End If

            Call gscatter(comm, n_held(:), n_loc, this_base_proc)

            Call gscatterv(comm, chbuf_scat(:), n_held(:), where_buff(:), &
                           chbuf(1:n_loc), this_base_proc)
            Call gscatterv(comm, iwrk_scat(:), n_held(:), where_buff(:), &
                           iwrk(1:n_loc), this_base_proc)
            Call gscatter_columns(comm, scatter_buffer_read(:, :), n_held(:), &
                                  where_buff(:), &
                                  scatter_buffer(1:wp_vals_per_at, 1:n_loc), &
                                  this_base_proc)

            ! Assign atoms to correct domains

            Do i = 1, n_loc
              config%natms = config%natms + 1

              ! Check safety by the upper bound of: atmnam,ltg,xxx,yyy,zzz &
              ! possibly vxx,vyy,vzz & possibly fxx,fyy,fzz as guided by xxx

              If (config%natms <= config%mxatms) Then
                config%atmnam(config%natms) = chbuf(i)
                config%ltg(config%natms) = iwrk(i)

                config%parts(config%natms)%xxx = scatter_buffer(1, i)
                config%parts(config%natms)%yyy = scatter_buffer(2, i)
                config%parts(config%natms)%zzz = scatter_buffer(3, i)

                If (levcfg /= 3) Then
                  If (levcfg > 0) Then
                    config%vxx(config%natms) = scatter_buffer(4, i)
                    config%vyy(config%natms) = scatter_buffer(5, i)
                    config%vzz(config%natms) = scatter_buffer(6, i)
                  Else
                    config%vxx(config%natms) = 0.0_wp
                    config%vyy(config%natms) = 0.0_wp
                    config%vzz(config%natms) = 0.0_wp
                  End If

                  If (levcfg > 1) Then
                    config%parts(config%natms)%fxx = scatter_buffer(7, i)
                    config%parts(config%natms)%fyy = scatter_buffer(8, i)
                    config%parts(config%natms)%fzz = scatter_buffer(9, i)
                  Else
                    config%parts(config%natms)%fxx = 0.0_wp
                    config%parts(config%natms)%fyy = 0.0_wp
                    config%parts(config%natms)%fzz = 0.0_wp
                  End If
                End If
              Else
                safe = .false.
              End If
            End Do
          End Do

          ! Check if all is dispatched fine

          max_fail = config%natms
          min_fail = config%natms
          Call gcheck(comm, safe)
          Call gmax(comm, max_fail)
          Call gmin(comm, min_fail)

          If (.not. safe) Then
            Write (messages(1), '(a,i0)') 'next error due to maximum number of atoms per domain set to : ', config%mxatms
            Write (messages(2), '(2(a,i0))') 'but maximum & minumum numbers of atoms per domain asked for : ', &
              max_fail, ' & ', min_fail
            Write (messages(3), '(a,i0)') 'estimated densvar value for passing this stage safely is : ', &
              Ceiling((dvar * (Real(max_fail, wp) / Real(config%mxatms, wp))**(1.0_wp / 1.7_wp) - 1.0_wp) * 100.0_wp)
            Call info(messages, 3, .true.)
            Call error(45)
          End If

        End If Extent_2

        ! Nullify dispatch counter

        indatm = 0

      End If Reorganize_buffer

    End Do

    ! If only detecting box dimensions for imcon == 0 or 6 or imc_n == 6

    If (l_xtr) Then
      Call gmax(comm, xhi)
      Call gmax(comm, yhi)
      Call gmax(comm, zhi)
    End If

    If (l_his) Then

      ! Skip to the EoFrame of HISTORY when not fast

      If (do_read .and. (.not. fast)) Then
        n_skip = Int(recs_per_at, offset_kind) * Int(megatm - first_at(my_read_proc_num), offset_kind)

        n_sk = Int(n_skip, li)
        n_jj = 73 * batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
          Do n_ii = 1_li, n_sk / n_jj
            forma = ' '
            Write (forma, '( "(", i0, "/)" )') n_jj
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, End=100)
          End Do
          n_ii = Mod(Int(n_skip, li), n_jj)
          If (n_ii > 0_li) Then
            forma = ' '
            Write (forma, '( "(", i0, "/)" )') n_ii
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, End=100)
          End If
        Else
          forma = ' '
          Write (forma, '( "(", i0, "/)" )') n_sk
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, End=100)
        End If
      End If

    Else

      If (do_read) Then

        If (io_read /= IO_READ_NETCDF) Then

          ! The last reader to check for EoFile in CONFIG
          ! and if none is hit to call error to abort

          If (first_at(my_read_proc_num) == megatm) Then
            recs_to_read = 1
            If (.not. fast) Then
              Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=forma, Iostat=ierr) rec_buff(:, 1:recs_to_read)
            Else
              Call io_read_batch(io, fh, rec_mpi_io, 1, rec_buff, ierr)
            End If
            safe = (ierr /= 0)
          End If

        Else

          ! As netCDF files have no real concept of line numbers,
          ! instead check the arrays are the correct size

          Call io_nc_get_dim(io, 'atom', fh, n_ats_in_file)
          safe = n_ats_in_file == megatm

        End If

      End If
      Call gcheck(comm, safe)
      If (.not. safe) Call error(58)

    End If

    If (do_read) Then
      Deallocate (rec_buff, Stat=fail(1))
      Deallocate (chbuf_read, iwrk_read, Stat=fail(2))
      Deallocate (axx_read, ayy_read, azz_read, Stat=fail(3))
      Deallocate (bxx_read, byy_read, bzz_read, Stat=fail(4))
      Deallocate (cxx_read, cyy_read, czz_read, Stat=fail(5))
      Deallocate (owner_read, Stat=fail(6))
      If (Any(fail(1:6) > 0)) Then
        Write (message, '(a)') 'read_config_parallel deallocation failure 2'
        Call error(0, message)
      End If
    End If

    Deallocate (first_at, orig_first_at, Stat=fail(1))
    Deallocate (n_held, where_buff, Stat=fail(2))
    Deallocate (chbuf, chbuf_scat, Stat=fail(3))
    Deallocate (iwrk, iwrk_scat, Stat=fail(4))
    Deallocate (scatter_buffer_read, Stat=fail(5))
    Deallocate (scatter_buffer, Stat=fail(6))
    If (Any(fail(1:6) > 0)) Then
      Write (message, '(a)') 'read_config_parallel deallocation failure 1'
      Call error(0, message)
    End If

    Return

    ! error exit for CONFIG file read

    100 Continue
    Call error(55)
  End Subroutine read_config_parallel

  Subroutine scan_config(config, megatm, dvar, levcfg, xhi, yhi, zhi, &
                         io, domain, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for raw scanning the contents of configuration file
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2014
    ! contrib   - i.j.bush april 2010
    ! contrib   - a.m.elena february 2017
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Integer,                  Intent(In   ) :: megatm
    Real(Kind=wp),            Intent(In   ) :: dvar
    Integer,                  Intent(  Out) :: levcfg
    Real(Kind=wp),            Intent(  Out) :: xhi, yhi, zhi
    Type(io_type),            Intent(InOut) :: io
    Type(domains_type),       Intent(In   ) :: domain
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=200)        :: record
    Character(Len=40)         :: fname, word
    Integer                   :: fh, i, io_read, j, recsz, totatm
    Integer(Kind=offset_kind) :: top_skip
    Logical                   :: fast, l_his, l_ind, l_xtr, safe, strict
    Real(Kind=wp)             :: buffer(1:4), cell_vecs(1:3, 1:3), xxx, yyy, zzz

! Some parameters and variables needed by io interfaces
! default record size

    safe = .true. ! we start safe
    l_ind = .false. ! no indeces needed
    strict = .false. ! not in a strict mode
    l_his = .false. ! not reading history
    l_xtr = .true. ! seeking extreme cell dimensions
    recsz = 73 ! ASCII record size for CONFIG

    ! Get type of I/O for reading

    Call io_get_parameters(io, user_method_read=io_read)

    ! Define filename ASCII or netCDF

    If (io_read /= IO_READ_NETCDF) Then
      fname = Trim(files(FILE_CONFIG)%filename)
    Else
      fname = Trim(files(FILE_CONFIG)%filename)//'nc'
    End If

    ! Check if we have a CONFIG

    If (comm%idnode == 0) Inquire (File=fname, Exist=safe)
    Call gcheck(comm, safe)
    If (.not. safe) Call error(124)

    ! Define/Detect the FAST reading status

    If (io_read == IO_READ_MASTER) Then

      fast = .false.

    Else If (io_read == IO_READ_NETCDF) Then

      fast = .true.

    Else

      ! Check if the system input file is a new style CONFIG:
      ! (i)  all lines are 72 ASCII characters long with
      !      a UNIX carriage return as end of line;
      ! (ii) LINE2 has the particles total value
      !      after values of levcfg and imcon.
      ! No fall back if users have mangled with further lines

      fast = .true.
      If (comm%idnode == 0) Then

        ! Open CONFIG

        Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)

        ! Read the CONFIG file header (TITLE record)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
          i = i + 1
          safe = .false.
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a1)', Advance='No', IOStat=j, End=10) record(i:i)
          safe = .true.
          If (j < 0) Go To 10
        End Do
        10 Continue
        fast = (fast .and. i == recsz)

        ! Read configuration level and image condition (RECORD2)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
          i = i + 1
          safe = .false.
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a1)', Advance='No', IOStat=j, End=20) record(i:i)
          safe = .true.
          If (j < 0) Go To 20
        End Do
        20 Continue
        fast = (fast .and. i == recsz)

        ! Read particles total value

        Call get_word(record, word); Call get_word(record, word)
        Call get_word(record, word); i = Nint(word_2_real(word, 0.0_wp, strict))
        fast = (fast .and. i == megatm)

      End If

      Call gsync(comm)
      Call gcheck(comm, safe, "enforce")
      Call gcheck(comm, fast, "enforce")

      If (.not. safe) Go To 50

      ! Close CONFIG

      If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()

    End If

    !!! SCAN HEADER

    If (io_read /= IO_READ_NETCDF) Then ! ASCII read

      ! Open CONFIG

      If (comm%idnode == 0) Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)

      ! Read TITLE record (file header)

      Call get_line(safe, files(FILE_CONFIG)%unit_no, record, comm)
      If (.not. safe) Go To 50

      Call strip_blanks(record)
      config%cfgname = Trim(record)

      ! Read configuration level and image condition

      Call get_line(safe, files(FILE_CONFIG)%unit_no, record, comm)
      If (.not. safe) Go To 50

      Call get_word(record, word)
      levcfg = Nint(word_2_real(word))

      ! halt execution if configuration level is unsupported

      If (levcfg < 0 .or. levcfg > 2) Call error(517)

      Call get_word(record, word)
      config%imcon = Nint(word_2_real(word))

      ! halt execution if image conventions is unsupported

      If (config%imcon < 0 .or. config%imcon > 7) Call error(514)

      ! specify MD cell (not defined for imcon=0)

      If (config%imcon /= 0) Then
        Call get_line(safe, files(FILE_CONFIG)%unit_no, record, comm)
        If (.not. safe) Go To 50
        Call get_word(record, word)
        config%cell(1) = word_2_real(word)
        Call get_word(record, word)
        config%cell(2) = word_2_real(word)
        Call get_word(record, word)
        config%cell(3) = word_2_real(word)

        Call get_line(safe, files(FILE_CONFIG)%unit_no, record, comm)
        If (.not. safe) Go To 50
        Call get_word(record, word)
        config%cell(4) = word_2_real(word)
        Call get_word(record, word)
        config%cell(5) = word_2_real(word)
        Call get_word(record, word)
        config%cell(6) = word_2_real(word)

        Call get_line(safe, files(FILE_CONFIG)%unit_no, record, comm)
        If (.not. safe) Go To 50
        Call get_word(record, word)
        config%cell(7) = word_2_real(word)
        Call get_word(record, word)
        config%cell(8) = word_2_real(word)
        Call get_word(record, word)
        config%cell(9) = word_2_real(word)
      End If

      ! Close CONFIG

      If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()
      Call gsync(comm)

    Else ! netCDF read

      ! Open CONFIG

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_open(io, io_read, comm%comm, fname, mode_rdonly, fh)

      i = 1 ! For config there is only one frame

      Call io_nc_get_att(io, 'title', fh, config%cfgname)

      Call io_nc_get_var(io, 'datalevel', fh, levcfg, i, 1)
      If (levcfg < 0 .or. levcfg > 2) Call error(517)

      Call io_nc_get_var(io, 'imageconvention', fh, config%imcon, i, 1)
      If (config%imcon < 0 .or. config%imcon > 7) Call error(514)

      Call io_nc_get_var(io, 'cell', fh, cell_vecs, [1, 1, i], [3, 3, 1])
      config%cell = Reshape(cell_vecs, [Size(config%cell)])

      ! Close CONFIG

      Call io_close(io, fh)

    End If

    ! ELABORATE SCAN

    totatm = 0
    xhi = 0.0_wp
    yhi = 0.0_wp
    zhi = 0.0_wp

    If (config%imcon == 0 .or. config%imcon == 6 .or. config%imc_n == 6) Then

      ! If MASTER read

      If (io_read == IO_READ_MASTER) Then

        If (comm%idnode == 0) Then

          ! Open CONFIG

          Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)

          ! Skip the header (we know exists from the basic scan)

          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*)
          Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*)
          If (config%imcon /= 0) Then
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*)
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*)
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*)
          End If

          ! Find the extreme dimensions for the system

          Do
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=40)
            Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30) xxx, yyy, zzz

            If (levcfg > 0) Then
              Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30)
              If (levcfg > 1) Read (Unit=files(FILE_CONFIG)%unit_no, Fmt=*, End=30)
            End If

            totatm = totatm + 1
            xhi = Max(xhi, Abs(xxx))
            yhi = Max(yhi, Abs(yyy))
            zhi = Max(zhi, Abs(zzz))
          End Do

          30 Continue ! catch error
          safe = .false.

        End If

        40 Continue ! catch EoF

        Call gsync(comm)
        Call gcheck(comm, safe, "enforce")

        If (.not. safe) Go To 50

        ! Close CONFIG

        If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()

        buffer(1) = xhi
        buffer(2) = yhi
        buffer(3) = zhi
        buffer(4) = Real(totatm, wp)
        Call gsum(comm, buffer)
        xhi = buffer(1)
        yhi = buffer(2)
        zhi = buffer(3)
        totatm = Nint(buffer(4))
        If (totatm /= megatm) Call error(58)

        ! If PROPER read

      Else

        ! Open CONFIG

        If (fast) Then
          Call io_set_parameters(io, user_comm=comm%comm)
          Call io_init(io, recsz)
          Call io_open(io, io_read, comm%comm, fname, mode_rdonly, fh)
        Else
          Open (Newunit=files(FILE_CONFIG)%unit_no, File=files(FILE_CONFIG)%filename)
        End If

        ! top_skip is header size

        If (io_read /= IO_READ_NETCDF) Then
          If (config%imcon == 0) Then
            top_skip = Int(2, offset_kind)
          Else
            top_skip = Int(5, offset_kind)
          End If
        Else
          top_skip = Int(1, offset_kind) ! This is now the frame = 1
        End If

        Call read_config_parallel(config, levcfg, dvar, l_ind, strict, megatm, l_his, l_xtr, &
                                  fast, fh, top_skip, xhi, yhi, zhi, io, domain, files, comm)

        ! Close CONFIG

        If (fast) Then
          Call io_close(io, fh)
          Call io_finalize(io)
        Else
          Call files(FILE_CONFIG)%close ()
        End If

      End If

    End If

    Return

    ! error exit for CONFIG file read

    50 Continue
    If (comm%idnode == 0) Call files(FILE_CONFIG)%close ()
    Call error(55)

  End Subroutine scan_config

  Subroutine scale_config(config, io, devel, netcdf, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for rescaling the crystallographic information
    ! from CONFIG to new lattice parameters and saving it in CFGSCL
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Type(io_type),            Intent(InOut) :: io
    Type(development_type),   Intent(In   ) :: devel
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(comms_type),         Intent(InOut) :: comm

    Integer         :: i, step
    Real(Kind=wp)   :: det, rcell(1:9), time, tstep, uuu, vvv, www
    Type(file_type) :: cfgscl

    ! Get the inverse cell matrix

    Call invert(config%cell, rcell, det)

    ! Rescale

    Do i = 1, config%natms
      uuu = config%parts(i)%xxx
      vvv = config%parts(i)%yyy
      www = config%parts(i)%zzz

      config%parts(i)%xxx = rcell(1) * uuu + rcell(4) * vvv + rcell(7) * www
      config%parts(i)%yyy = rcell(2) * uuu + rcell(5) * vvv + rcell(8) * www
      config%parts(i)%zzz = rcell(3) * uuu + rcell(6) * vvv + rcell(9) * www

      uuu = config%parts(i)%xxx
      vvv = config%parts(i)%yyy
      www = config%parts(i)%zzz

      config%parts(i)%xxx = devel%cels(1) * uuu + devel%cels(4) * vvv + devel%cels(7) * www
      config%parts(i)%yyy = devel%cels(2) * uuu + devel%cels(5) * vvv + devel%cels(8) * www
      config%parts(i)%zzz = devel%cels(3) * uuu + devel%cels(6) * vvv + devel%cels(9) * www
    End Do

    ! Write REVCON

    Call cfgscl%init('CFGSCL')
    step = 0 ! no steps done
    tstep = 0.0_wp ! no step exists
    time = 0.0_wp ! time is not relevant

    rcell = config%cell; config%cell = devel%cels
    Call write_config(config, cfgscl, devel%lvcfscl, step, tstep, io, time, netcdf, comm)
    config%cell = rcell

  End Subroutine scale_config

  Subroutine write_config(config, cfile, levcfg, step, tstep, io, time, netcdf, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for writing configuration file
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    ! contrib   - i.j.bush
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Type(file_type),          Intent(InOut) :: cfile
    Integer,                  Intent(In   ) :: levcfg, step
    Real(Kind=wp),            Intent(In   ) :: tstep
    Type(io_type),            Intent(InOut) :: io
    Real(Kind=wp),            Intent(In   ) :: time
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(comms_type),         Intent(InOut) :: comm

    Character                                      :: lf
    Character(Len=1), Allocatable, Dimension(:, :) :: chbat
    Character(Len=1024)                            :: fname
    Character(Len=256)                             :: message
    Character(Len=8), Allocatable, Dimension(:)    :: chbuf
    Character(Len=recsz)                           :: record
    Integer                                        :: batsz, fail(1:4), fh, i, ierr, io_write, &
                                                      jatms, jdnode, jj, k
    Integer(Kind=li)                               :: rec, rec1
    Integer(Kind=offset_kind)                      :: rec_mpi_io
    Integer, Allocatable, Dimension(:)             :: iwrk, n_atm
    Logical                                        :: ready
    Real(Kind=wp)                                  :: angles(1:3), cell_vecs(1:3, 1:3), &
                                                      celprp(1:10), lengths(1:3)
    Real(Kind=wp), Allocatable, Dimension(:)       :: axx, ayy, azz, bxx, byy, bzz, cxx, cyy, czz

! record line
! Some parameters and variables needed by io interfaces

    ! Get write method buffer size and line feed character

    Call io_get_parameters(io, user_method_write=io_write)
    Call io_get_parameters(io, user_buffer_size_write=batsz)
    Call io_get_parameters(io, user_line_feed=lf)

    ! Get offsets and define batch

    fail = 0
    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT .or. &
        io_write == IO_WRITE_UNSORTED_MASTER) Then
      Allocate (n_atm(0:comm%mxnode), Stat=fail(1))
      Allocate (chbat(1:recsz, 1:batsz), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'write_config allocation failure 0'
        Call error(0, message)
      End If

      chbat = ' '
      n_atm = 0; n_atm(comm%idnode + 1) = config%natms
      Call gsum(comm, n_atm)
      n_atm(0) = Sum(n_atm(0:comm%idnode))
    End If

    ! Notes:
    ! the MPI-I/O records are numbered from 0 (not 1)
    ! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
    !   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
    !   both declared as: Integer(Kind = offset_kind)

    ! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT) Then

      ! Write header only at start, where just one node is needed
      ! Start of file

      rec_mpi_io = Int(1 - 1, offset_kind)
      jj = 0
      If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm=comm_self)
        Call io_init(io, recsz)
        Call io_delete(io, cfile%filename, comm) ! Sort existence issues
        Call io_open(io, io_write, comm_self, cfile%filename, mode_wronly + mode_create, fh)

        ! Accumulate header

        Write (record, Fmt='(a72,a1)') config%cfgname(1:72), lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        Write (record, Fmt='(4i10,1p,2e16.7,a1)') levcfg, config%imcon, config%megatm, step, tstep, time, lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        ! Accumulate header - optional cell information (if present)

        If (config%imcon > 0) Then
          Do i = 0, 2
            Write (record, Fmt='(3f20.10,a12,a1)') &
              config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3), Repeat(' ', 12), lf
            jj = jj + 1
            Do k = 1, recsz
              chbat(k, jj) = record(k:k)
            End Do
          End Do
        End If

        ! Dump header

        Call io_write_batch(io, fh, rec_mpi_io, jj, chbat)

        Call io_close(io, fh)
        Call io_finalize(io)

      Else

        jj = jj + 2
        If (config%imcon > 0) jj = jj + 3

      End If
      Call gsync(comm)

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_init(io, recsz)
!      Call io_delete(io, cfile%filename, comm)
      Call io_open(io, io_write, comm%comm, cfile%filename, mode_wronly, fh)

      ! Start of file (updated)

      rec_mpi_io = Int(jj, offset_kind) + Int(n_atm(0), offset_kind) * Int(levcfg + 2, offset_kind)
      jj = 0
      Do i = 1, config%natms
        Write (record, Fmt='(a8,i10,a54,a1)') config%atmnam(i), config%ltg(i), Repeat(' ', 54), lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        Write (record, Fmt='(3g20.10,a12,a1)') config%parts(i)%xxx, config%parts(i)%yyy, config%parts(i)%zzz, Repeat(' ', 12), lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        If (levcfg > 0) Then
          Write (record, Fmt='(3g20.10,a12,a1)') config%vxx(i), config%vyy(i), config%vzz(i), Repeat(' ', 12), lf
          jj = jj + 1
          Do k = 1, recsz
            chbat(k, jj) = record(k:k)
          End Do

          If (levcfg > 1) Then
            Write (record, '(3g20.10,a12,a1)') config%parts(i)%fxx, config%parts(i)%fyy, config%parts(i)%fzz, Repeat(' ', 12), lf
            jj = jj + 1
            Do k = 1, recsz
              chbat(k, jj) = record(k:k)
            End Do
          End If
        End If

        ! Dump batch and update start of file

        If (jj + levcfg + 2 >= batsz .or. i == config%natms) Then
          Call io_write_batch(io, fh, rec_mpi_io, jj, chbat)
          rec_mpi_io = rec_mpi_io + Int(jj, offset_kind)
          jj = 0
        End If
      End Do

      Call io_close(io, fh)
      Call io_finalize(io)

      ! UNSORTED Serial Direct Access FORTRAN

    Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

      Allocate (chbuf(1:config%mxatms), iwrk(1:config%mxatms), Stat=fail(1))
      Allocate (axx(1:config%mxatms), ayy(1:config%mxatms), azz(1:config%mxatms), Stat=fail(2))
      Allocate (bxx(1:config%mxatms), byy(1:config%mxatms), bzz(1:config%mxatms), Stat=fail(3))
      Allocate (cxx(1:config%mxatms), cyy(1:config%mxatms), czz(1:config%mxatms), Stat=fail(4))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'write_config allocation failure'
        Call error(0, message)
      End If

      ! node 0 handles I/O
      ! Start of file

      rec = Int(0, li)
      jj = 0
      If (comm%idnode == 0) Then

        ! Write configuration data to new configuration file

        Open (Newunit=cfile%unit_no, File=cfile%filename, Form='formatted', Access='direct', Recl=recsz, Status='replace')

        ! Accumulate header

        Write (record, Fmt='(a72,a1)') config%cfgname(1:72), lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        Write (record, Fmt='(4i10,1p,2e16.7,a1)') levcfg, config%imcon, config%megatm, step, tstep, time, lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        ! Accumulate header - optional cell information (if present)

        If (config%imcon > 0) Then
          Do i = 0, 2
            Write (record, Fmt='(3f20.10,a12,a1)') &
              config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3), Repeat(' ', 12), lf
            jj = jj + 1
            Do k = 1, recsz
              chbat(k, jj) = record(k:k)
            End Do
          End Do
        End If

        ! Dump header and update start of file

        Write (Unit=cfile%unit_no, Fmt='(73a)', Rec=rec + Int(1, li)) (chbat(:, k), k=1, jj)
        rec = Int(jj, li)
        jj = 0

        Do i = 1, config%natms
          iwrk(i) = config%ltg(i)
          chbuf(i) = config%atmnam(i)

          axx(i) = config%parts(i)%xxx
          ayy(i) = config%parts(i)%yyy
          azz(i) = config%parts(i)%zzz

          If (levcfg > 0) Then
            bxx(i) = config%vxx(i)
            byy(i) = config%vyy(i)
            bzz(i) = config%vzz(i)

            If (levcfg > 1) Then
              cxx(i) = config%parts(i)%fxx
              cyy(i) = config%parts(i)%fyy
              czz(i) = config%parts(i)%fzz
            End If
          End If
        End Do

        jatms = config%natms
        ready = .true.
        Do jdnode = 0, comm%mxnode - 1
          If (jdnode > 0) Then
            Call gsend(comm, ready, jdnode, WriteConf_tag)

            Call grecv(comm, jatms, jdnode, WriteConf_tag)
            If (jatms > 0) Then
              Call grecv(comm, chbuf(1:jatms), jdnode, WriteConf_tag)
              Call grecv(comm, iwrk(1:jatms), jdnode, WriteConf_tag)

              Call grecv(comm, axx(1:jatms), jdnode, WriteConf_tag)
              Call grecv(comm, ayy(1:jatms), jdnode, WriteConf_tag)
              Call grecv(comm, azz(1:jatms), jdnode, WriteConf_tag)

              If (levcfg > 0) Then
                Call grecv(comm, bxx(1:jatms), jdnode, WriteConf_tag)
                Call grecv(comm, byy(1:jatms), jdnode, WriteConf_tag)
                Call grecv(comm, bzz(1:jatms), jdnode, WriteConf_tag)

                If (levcfg > 1) Then
                  Call grecv(comm, cxx(1:jatms), jdnode, WriteConf_tag)
                  Call grecv(comm, cyy(1:jatms), jdnode, WriteConf_tag)
                  Call grecv(comm, czz(1:jatms), jdnode, WriteConf_tag)
                End If
              End If
            End If
          End If

          jj = 0
          Do i = 1, jatms
            Write (record, Fmt='(a8,i10,a54,a1)') config%atmnam(i), iwrk(i), Repeat(' ', 54), lf
            jj = jj + 1
            Do k = 1, recsz
              chbat(k, jj) = record(k:k)
            End Do

            Write (record, Fmt='(3g20.10,a12,a1)') axx(i), ayy(i), azz(i), Repeat(' ', 12), lf
            jj = jj + 1
            Do k = 1, recsz
              chbat(k, jj) = record(k:k)
            End Do

            If (levcfg > 0) Then
              Write (record, Fmt='(3g20.10,a12,a1)') bxx(i), byy(i), bzz(i), Repeat(' ', 12), lf
              jj = jj + 1
              Do k = 1, recsz
                chbat(k, jj) = record(k:k)
              End Do

              If (levcfg > 1) Then
                Write (record, Fmt='(3g20.10,a12,a1)') cxx(i), cyy(i), czz(i), Repeat(' ', 12), lf
                jj = jj + 1
                Do k = 1, recsz
                  chbat(k, jj) = record(k:k)
                End Do
              End If
            End If

            ! Dump batch and update start of file

            If (jj + levcfg + 2 >= batsz .or. i == jatms) Then
              Write (Unit=cfile%unit_no, Fmt='(73a)', Rec=rec + Int(1, li)) (chbat(:, k), k=1, jj)
              rec = rec + Int(jj, li)
              jj = 0
            End If
          End Do
        End Do

        Call cfile%close ()

      Else

        Call grecv(comm, ready, 0, WriteConf_tag)

        Call gsend(comm, config%natms, 0, WriteConf_tag)
        If (config%natms > 0) Then
          Call gsend(comm, config%atmnam(1:config%natms), 0, WriteConf_tag)
          Call gsend(comm, config%ltg(1:config%natms), 0, WriteConf_tag)

          Call gsend(comm, config%parts(1:config%natms), 0, WriteConf_tag)
          If (levcfg > 0) Then
            Call gsend(comm, config%vxx(1:config%natms), 0, WriteConf_tag)
            Call gsend(comm, config%vyy(1:config%natms), 0, WriteConf_tag)
            Call gsend(comm, config%vzz(1:config%natms), 0, WriteConf_tag)

          End If
        End If

      End If

      Deallocate (chbuf, iwrk, Stat=fail(1))
      Deallocate (axx, ayy, azz, Stat=fail(2))
      Deallocate (bxx, byy, bzz, Stat=fail(3))
      Deallocate (cxx, cyy, czz, Stat=fail(4))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'write_config deallocation failure'
        Call error(0, message)
      End If

      ! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

    Else If (io_write == IO_WRITE_SORTED_MPIIO .or. &
             io_write == IO_WRITE_SORTED_DIRECT .or. &
             io_write == IO_WRITE_SORTED_NETCDF) Then

      ! name convention

      If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = cfile%filename
      Else
        fname = Trim(cfile%filename)//'.nc'
      End If

      ! Write header only at start, where just one node is needed
      ! Start of file

      rec_mpi_io = Int(1 - 1, offset_kind)
      jj = 0
      If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm=comm_self)
        Call io_init(io, recsz)
        Call io_delete(io, fname, comm) ! Sort existence issues
        If (io_write == IO_WRITE_SORTED_NETCDF) Then
          Call io_nc_create(netcdf, comm_self, fname, config%cfgname, config%megatm)
        End If
        Call io_open(io, io_write, comm_self, fname, mode_wronly + mode_create, fh)

        ! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

          ! Write header

          Write (record, Fmt='(a72,a1)') config%cfgname(1:72), lf
          Call io_write_record(io, fh, Int(jj, offset_kind), record)
          jj = jj + 1

          Write (record, Fmt='(4i10,1p,2e16.7,a1)') levcfg, config%imcon, config%megatm, step, tstep, time, lf
          Call io_write_record(io, fh, Int(jj, offset_kind), record)
          jj = jj + 1

          ! Write optional cell information (if present)

          If (config%imcon > 0) Then
            Do i = 0, 2
              Write (record, '( 3f20.10, a12, a1 )') &
                config%cell(1 + i * 3:3 + i * 3), Repeat(' ', 12), lf
              Call io_write_record(io, fh, Int(jj, offset_kind), record)
              jj = jj + 1
            End Do
          End If

        Else ! netCDF write

          jj = 1 ! For config there is only one frame

          Call io_nc_put_var(io, 'time', fh, time, jj, 1)
          Call io_nc_put_var(io, 'step', fh, step, jj, 1)
          Call io_nc_put_var(io, 'datalevel', fh, levcfg, jj, 1)
          Call io_nc_put_var(io, 'imageconvention', fh, config%imcon, jj, 1)
          Call io_nc_put_var(io, 'timestep', fh, tstep, jj, 1)

          If (config%imcon > 0) Then
            Call dcell(config%cell, celprp) ! get cell properties

            cell_vecs = Reshape(config%cell, [3, 3])

            lengths(1) = celprp(1)
            lengths(2) = celprp(2)
            lengths(3) = celprp(3)

            angles(1) = Acos(celprp(5))
            angles(2) = Acos(celprp(6))
            angles(3) = Acos(celprp(4))
            angles = angles * 180.0_wp / (4.0_wp * Atan(1.0_wp)) ! Convert to degrees

            ! Print

            Call io_nc_put_var(io, 'cell', fh, cell_vecs, [1, 1, jj], [3, 3, 1])
            Call io_nc_put_var(io, 'cell_lengths', fh, lengths, [1, jj], [3, 1])
            Call io_nc_put_var(io, 'cell_angles', fh, angles, [1, jj], [3, 1])
          End If

        End If

        Call io_close(io, fh)
        Call io_finalize(io)

      Else

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
          jj = jj + 2
          If (config%imcon > 0) jj = jj + 3
        Else
          jj = 1
        End If

      End If
      Call gsync(comm)

      ! Write the rest

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_init(io, recsz)
      Call io_open(io, io_write, comm%comm, fname, mode_wronly, fh)

      rec_mpi_io = rec_mpi_io + Int(jj, offset_kind)
      Call io_write_sorted_file(io, fh, levcfg, IO_RESTART, rec_mpi_io, config%natms, &
                                config%ltg, config%atmnam, [0.0_wp], [0.0_wp], config%parts, &
                                config%vxx, config%vyy, config%vzz, IO_SUBSET_POSITIONS + IO_SUBSET_FORCES, ierr)

      If (ierr /= 0) Then
        Select Case (ierr)
        Case (IO_BASE_COMM_NOT_SET)
          Call error(1050)
        Case (IO_ALLOCATION_ERROR)
          Call error(1053)
        Case (IO_UNKNOWN_WRITE_OPTION)
          Call error(1056)
        Case (IO_UNKNOWN_WRITE_LEVEL)
          Call error(1059)
        End Select
      End If

      Call io_close(io, fh)
      Call io_finalize(io)

      ! SORTED Serial Direct Access FORTRAN

    Else If (io_write == IO_WRITE_SORTED_MASTER) Then

      Allocate (chbuf(1:config%mxatms), iwrk(1:config%mxatms), Stat=fail(1))
      Allocate (axx(1:config%mxatms), ayy(1:config%mxatms), azz(1:config%mxatms), Stat=fail(2))
      Allocate (bxx(1:config%mxatms), byy(1:config%mxatms), bzz(1:config%mxatms), Stat=fail(3))
      Allocate (cxx(1:config%mxatms), cyy(1:config%mxatms), czz(1:config%mxatms), Stat=fail(4))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'write_config allocation failure'
        Call error(0, message)
      End If

      ! node 0 handles I/O
      ! Start of file

      rec = Int(0, li)
      If (comm%idnode == 0) Then

        ! Write configuration data to new configuration file

        Open (Newunit=cfile%unit_no, File=cfile%filename, Form='formatted', Access='direct', Recl=recsz, Status='replace')

        ! Write header

        rec = rec + Int(1, li)
        Write (Unit=cfile%unit_no, Fmt='(a72,a1)', Rec=rec) config%cfgname(1:72), lf
        rec = rec + Int(1, li)
        Write (Unit=cfile%unit_no, Fmt='(4i10,1p,2e16.7,a1)', Rec=rec) levcfg, config%imcon, config%megatm, step, tstep, time, lf

        ! Write optional cell information (if present)

        If (config%imcon > 0) Then
          Do i = 0, 2
            rec = rec + Int(1, li)
            Write (Unit=cfile%unit_no, Fmt='(3f20.10,a12,a1)', Rec=rec) &
              config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3), Repeat(' ', 12), lf
          End Do
        End If

        Do i = 1, config%natms
          iwrk(i) = config%ltg(i)
          chbuf(i) = config%atmnam(i)

          axx(i) = config%parts(i)%xxx
          ayy(i) = config%parts(i)%yyy
          azz(i) = config%parts(i)%zzz

          If (levcfg > 0) Then
            bxx(i) = config%vxx(i)
            byy(i) = config%vyy(i)
            bzz(i) = config%vzz(i)

            If (levcfg > 1) Then
              cxx(i) = config%parts(i)%fxx
              cyy(i) = config%parts(i)%fyy
              czz(i) = config%parts(i)%fzz
            End If
          End If
        End Do

        jatms = config%natms
        ready = .true.
        Do jdnode = 0, comm%mxnode - 1
          If (jdnode > 0) Then
            Call gsend(comm, ready, jdnode, WriteConf_tag)

            Call grecv(comm, jatms, jdnode, WriteConf_tag)
            If (jatms > 0) Then
              Call gsend(comm, chbuf(1:jatms), jdnode, WriteConf_tag)
              Call grecv(comm, iwrk(1:jatms), jdnode, WriteConf_tag)

              Call grecv(comm, axx(1:jatms), jdnode, WriteConf_tag)
              Call grecv(comm, ayy(1:jatms), jdnode, WriteConf_tag)
              Call grecv(comm, azz(1:jatms), jdnode, WriteConf_tag)

              If (levcfg > 0) Then
                Call grecv(comm, bxx(1:jatms), jdnode, WriteConf_tag)
                Call grecv(comm, byy(1:jatms), jdnode, WriteConf_tag)
                Call grecv(comm, bzz(1:jatms), jdnode, WriteConf_tag)

                If (levcfg > 1) Then
                  Call grecv(comm, cxx(1:jatms), jdnode, WriteConf_tag)
                  Call grecv(comm, cyy(1:jatms), jdnode, WriteConf_tag)
                  Call grecv(comm, czz(1:jatms), jdnode, WriteConf_tag)
                End If
              End If
            End If
          End If

          Do i = 1, jatms
            rec1 = rec + Int(iwrk(i) - 1, li) * Int(levcfg + 2) + Int(1, li)
            Write (Unit=cfile%unit_no, Fmt='(a8,i10,a54,a1)', Rec=rec1) chbuf(i), iwrk(i), Repeat(' ', 54), lf
            rec1 = rec1 + Int(1, li)
            Write (Unit=cfile%unit_no, Fmt='(3g20.10,a12,a1)', Rec=rec1) axx(i), ayy(i), azz(i), Repeat(' ', 12), lf

            If (levcfg > 0) Then
              rec1 = rec1 + Int(1, li)
              Write (Unit=cfile%unit_no, Fmt='(3g20.10,a12,a1)', Rec=rec1) bxx(i), byy(i), bzz(i), Repeat(' ', 12), lf

              If (levcfg > 1) Then
                rec1 = rec1 + Int(1, li)
                Write (Unit=cfile%unit_no, Fmt='(3g20.10,a12,a1)', Rec=rec1) cxx(i), cyy(i), czz(i), Repeat(' ', 12), lf
              End If
            End If
          End Do
        End Do

        Call cfile%close ()

      Else

        Call grecv(comm, ready, 0, WriteConf_tag)

        Call gsend(comm, config%natms, 0, WriteConf_tag)
        If (config%natms > 0) Then
          Call gsend(comm, config%atmnam(1:config%natms), 0, WriteConf_tag)
          Call gsend(comm, config%ltg(1:config%natms), 0, WriteConf_tag)

          Call gsend(comm, config%parts(1:config%natms), 0, WriteConf_tag)

          If (levcfg > 0) Then
            Call gsend(comm, config%vxx(1:config%natms), 0, WriteConf_tag)
            Call gsend(comm, config%vyy(1:config%natms), 0, WriteConf_tag)
            Call gsend(comm, config%vzz(1:config%natms), 0, WriteConf_tag)

          End If
        End If

      End If

      Deallocate (chbuf, iwrk, Stat=fail(1))
      Deallocate (axx, ayy, azz, Stat=fail(2))
      Deallocate (bxx, byy, bzz, Stat=fail(3))
      Deallocate (cxx, cyy, czz, Stat=fail(4))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'write_config deallocation failure'
        Call error(0, message)
      End If

    End If

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT .or. &
        io_write == IO_WRITE_UNSORTED_MASTER) Then
      Deallocate (n_atm, Stat=fail(1))
      Deallocate (chbat, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'write_config deallocation failure 0'
        Call error(0, message)
      End If
    End If

    Call gsync(comm)

  End Subroutine write_config

  Subroutine getcom_arrays(txx, tyy, tzz, config, com, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to calculate system centre of mass position
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov october 2012
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Dimension(:), Intent(In   ) :: txx, tyy, tzz
    Type(configuration_type),    Intent(InOut) :: config
    Real(Kind=wp), Dimension(:), Intent(  Out) :: com
    Type(comms_type),            Intent(InOut) :: comm

    Integer :: i

    ! total system mass

    If (config%newjob_totmas) Then
      config%newjob_totmas = .false.

      config%totmas = 0.0_wp
      Do i = 1, config%natms
        If (config%lfrzn(i) == 0) config%totmas = config%totmas + config%weight(i)
      End Do

      Call gsum(comm, config%totmas)
    End If

    com = 0.0_wp

    Do i = 1, config%natms
      If (config%lfrzn(i) == 0) Then
        com(1) = com(1) + config%weight(i) * txx(i)
        com(2) = com(2) + config%weight(i) * tyy(i)
        com(3) = com(3) + config%weight(i) * tzz(i)
      End If
    End Do

    Call gsum(comm, com)
    If (config%totmas >= zero_plus) com = com / config%totmas

  End Subroutine getcom_arrays

  Subroutine getcom_parts(config, com, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to calculate system centre of mass position
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov october 2012
    !
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type),      Intent(InOut) :: config
    Real(Kind=wp), Dimension(1:3), Intent(  Out) :: com
    Type(comms_type),              Intent(InOut) :: comm

    Integer :: i

    ! total system mass

    If (config%newjob_totmas) Then
      config%newjob_totmas = .false.

      config%totmas = 0.0_wp
      Do i = 1, config%natms
        If (config%lfrzn(i) == 0) config%totmas = config%totmas + config%weight(i)
      End Do

      Call gsum(comm, config%totmas)
    End If

    com = 0.0_wp

    Do i = 1, config%natms
      If (config%lfrzn(i) == 0) Then
        com(1) = com(1) + config%weight(i) * config%parts(i)%xxx
        com(2) = com(2) + config%weight(i) * config%parts(i)%yyy
        com(3) = com(3) + config%weight(i) * config%parts(i)%zzz
      End If
    End Do

    Call gsum(comm, com)
    If (config%totmas >= zero_plus) com = com / config%totmas

  End Subroutine getcom_parts

  Subroutine getcom_mol(config, istart, ifinish, cmm, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to calculate a molecule's mass and COM
    !
    ! istart  - the global index of the first atom of the molecule
    ! ifinish - the global index of the last atom of the molecule
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Integer,                  Intent(In   ) :: istart, ifinish
    Real(Kind=wp),            Intent(  Out) :: cmm(0:3)
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail, i, j, k
    Real(Kind=wp)              :: mass, r(1:3)
    Real(Kind=wp), Allocatable :: mol(:, :)

    fail = 0
    Allocate (mol(1:(ifinish - istart + 1), 0:3), Stat=fail)
    If (fail > 0) Then
      Write (message, '(1x,a,i0)') 'getcom_mol allocation failure', fail
      Call error(0, message)
    End If

    ! Initialise

    mass = 0.0_wp
    cmm = 0.0_wp

    mol = 0.0_wp
    Do i = 1, config%natms
      j = config%ltg(i)
      If (j >= istart .and. j <= ifinish) Then
        k = j - istart + 1

        mol(k, 0) = config%weight(i)
        mol(k, 1) = config%parts(i)%xxx
        mol(k, 2) = config%parts(i)%yyy
        mol(k, 3) = config%parts(i)%zzz
      End If
    End Do

    Call gsum(comm, mol)

    r(1) = mol(1, 1)
    r(2) = mol(1, 2)
    r(3) = mol(1, 3)

    mol(:, 1) = mol(:, 1) - r(1)
    mol(:, 2) = mol(:, 2) - r(2)
    mol(:, 3) = mol(:, 3) - r(3)

    k = ifinish - istart + 1
    Call images(config%imcon, config%cell, k, mol(:, 1), mol(:, 2), mol(:, 3))

    mol(:, 1) = mol(:, 1) + r(1)
    mol(:, 2) = mol(:, 2) + r(2)
    mol(:, 3) = mol(:, 3) + r(3)

    Do i = 1, k
      mass = mass + mol(i, 0)
      cmm(0) = cmm(0) + mol(i, 0) * Real(1 - config%lfrzn(i), wp)
      cmm(1) = cmm(1) + mol(i, 0) * mol(i, 1)
      cmm(2) = cmm(2) + mol(i, 0) * mol(i, 2)
      cmm(3) = cmm(3) + mol(i, 0) * mol(i, 3)
    End Do

    If (cmm(0) >= zero_plus) cmm(1:3) = cmm(1:3) / mass

    fail = 0
    Deallocate (mol, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'getcom_mol deallocation failure'
      Call error(0, message)
    End If

  End Subroutine getcom_mol

  Subroutine freeze_atoms(config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to quench forces and velocities on 'frozen' atoms
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov july 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config

    Integer :: i

    Do i = 1, config%natms
      If (config%lfrzn(i) /= 0) Then
        config%vxx(i) = 0.0_wp
        config%vyy(i) = 0.0_wp
        config%vzz(i) = 0.0_wp

        config%parts(i)%fxx = 0.0_wp
        config%parts(i)%fyy = 0.0_wp
        config%parts(i)%fzz = 0.0_wp
      End If
    End Do

  End Subroutine freeze_atoms

  Subroutine origin_config(config, io, devel, netcdf, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for translating the origin of the MD box as
    ! defined in CONFIG by the (devel%xorg,devel%yorg,devel%zorg) vector
    ! and saving it in CFGORG
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Type(io_type),            Intent(InOut) :: io
    Type(development_type),   Intent(In   ) :: devel
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(comms_type),         Intent(InOut) :: comm

    Integer         :: i, step
    Real(Kind=wp)   :: time, tstep
    Type(file_type) :: cfgorg

    ! Translate

    Do i = 1, config%natms
      config%parts(i)%xxx = config%parts(i)%xxx + devel%xorg
      config%parts(i)%yyy = config%parts(i)%yyy + devel%yorg
      config%parts(i)%zzz = config%parts(i)%zzz + devel%zorg
    End Do

    ! Restore periodic boundaries

    Call pbcshift(config%imcon, config%cell, config%natms, config%parts)

    ! Write REVCON

    Call cfgorg%init('CFGORG')
    step = 0 ! no steps done
    tstep = 0.0_wp ! no step exists
    time = 0.0_wp ! time is not relevant

    Call write_config(config, cfgorg, devel%lvcforg, step, tstep, io, time, netcdf, comm)

  End Subroutine origin_config

  !> @brief Initialise coordinate_buffer_type
  !!
  !! @param[in]  comm          Object containing MPI communicator
  !! @param[in]  size_coords   Size of 1D packed buffer
  !
  Subroutine initialise_coordinate_buffer_type(this, comm, megatm)
    Class(coordinate_buffer_type), Intent(InOut) :: this
    Type(comms_type),              Intent(In   ) :: comm
    Integer,                       Intent(In   ) :: megatm

    Allocate (this%coords(megatm * 3))
    !Allocate(this%atmnam(megatm))
    Call this%mpi%init(comm)
  End Subroutine initialise_coordinate_buffer_type

  !> @brief Deallocate coordinate_buffer_type
  Subroutine finalise_coordinate_buffer_type(this)
    Class(coordinate_buffer_type), Intent(InOut) :: this

    Deallocate (this%coords)
    !Deallocate(this%atmnam)
    Call this%mpi%finalise()
  End Subroutine finalise_coordinate_buffer_type

  !> @brief Gather atomic coordinates from each process into a packed,
  !!  1D array and broadcast
  !!
  !! @param[inout]  comm             Object containing MPI communicator
  !! @param[in]     config           Object containing configuration data
  !! @param[in]     to_master_only   Logic, broadcast to master only
  !! @param[inout]  gathered         Object containing packed coordinates
  !!                                 and mpi index arrays for unpacking
  !
  Subroutine gather_coordinates(comm, config, to_master_only, gathered)
    Type(comms_type),             Intent(InOut) :: comm
    Type(configuration_type),     Intent(In   ) :: config
    Logical,                      Intent(In   ) :: to_master_only
    Type(coordinate_buffer_type), Intent(InOut) :: gathered

    Character(Len=100)    :: error_message
    Character(Len=5)      :: str1, str2
    Integer               :: i, i2, i3
    Real(wp), Allocatable :: coords(:)

    Call assert(Allocated(gathered%coords), &
                'Packed coordinate array not allocated prior to gathering coordinates')
    Call assert(Size(gathered%coords) == config%megatm * 3, &
                'Packed coordinate array does not contain megatm*3 elements')

    !Create local 1D array to contain x,y & z components of each local set of atoms
    Allocate (coords(config%natms * 3))
    coords = 0._wp

    Do i = 1, config%natms
      i2 = config%natms + i
      i3 = (2 * config%natms) + i
      coords(i) = config%parts(i)%xxx
      coords(i2) = config%parts(i)%yyy
      coords(i3) = config%parts(i)%zzz
    End Do

    !Generate record_size and displacement arrays required by gatherv
    Call gatherv_scatterv_index_arrays(comm, Size(coords), &
                                       gathered%mpi%counts, gathered%mpi%displ)

    !Check total number of atoms in system
    Write (str1, '(I5)') Int(Sum(gathered%mpi%counts) / 3)
    Write (str2, '(I5)') config%megatm
    error_message = 'Disagreement in total number of atoms found from summing all'// &
                    'elements of MPI record size vs megatm: '//str1//','//str2
    Call assert(Int(Sum(gathered%mpi%counts) / 3) == config%megatm, error_message)

    If (to_master_only) Then
      !Gather coords from each process into gathered%coords on master
      Call ggatherv(comm, coords, gathered%mpi%counts, &
                    gathered%mpi%displ, gathered%coords)
    Else
      !Gather coords from each process into gathered%coords on all processes
      Call gallgatherv(comm, coords, gathered%mpi%counts, &
                       gathered%mpi%displ, gathered%coords)
    Endif

  End Subroutine gather_coordinates

  !! @brief Unpack \p gathered%coords into 2D array
  !!
  !! Could be a member function of coordinate_buffer_type
  !!
  !! @param[in]  comm             Object containing MPI communicator
  !! @param[in]  config           Object containing configuration data
  !! @param[in]  gathered         Object containing packed coordinates
  !!                              and mpi index arrays for unpacking
  !! @param[inout] coords         Unpacked atomic coordinates
  !
  Subroutine unpack_gathered_coordinates(comm, config, gathered, coords)
    Type(comms_type),             Intent(In   ) :: comm
    Type(configuration_type),     Intent(In   ) :: config
    Type(coordinate_buffer_type), Intent(In   ) :: gathered
    Real(wp), Allocatable,        Intent(InOut) :: coords(:, :)

    Integer :: cnt, ia, ip, ix, iy, iz, natms_local

    Call assert(Allocated(coords), "coords not allocated")
    Call assert(Size(coords, 1) == 3, "Dim 1 of coords /= 3")
    Call assert(Size(coords, 2) == config%megatm, "Size(coords,2) /= config%megatm")

    cnt = 0
    Do ip = 1, comm%mxnode
      natms_local = Int(gathered%mpi%counts(ip) / 3)
      Do ia = 1, natms_local
        !Composite index cnt == ia+int(MPI_buffer%displ(ip)/3)
        cnt = cnt + 1
        !mpi process offset + xyz offset + local_index
        ix = gathered%mpi%displ(ip) + ia
        iy = gathered%mpi%displ(ip) + natms_local + ia
        iz = gathered%mpi%displ(ip) + natms_local * 2 + ia
        coords(1:3, cnt) = (/gathered%coords(ix), gathered%coords(iy), gathered%coords(iz)/)
      End Do
    End Do

  End Subroutine unpack_gathered_coordinates

  !> @brief Gather atom names from each process into one array and
  !!  broadcast
  !!
  !! Because gathered_atom_name is not packed, one does not require
  !! the corresponding MPI index arrays
  !!
  !! @param[inout]  comm             Object containing MPI communicator
  !! @param[in]     config           Object containing configuration data
  !! @param[in]     to_master_only   Logic, broadcast to master only
  !! @param[inout]  atmnam           Gathered atomic names
  !
  Subroutine gather_atomic_names(comm, config, to_master_only, atmnam)
    Type(comms_type),                       Intent(InOut) :: comm
    Type(configuration_type),               Intent(In   ) :: config
    Logical,                                Intent(In   ) :: to_master_only
    Character(Len=len_atmnam), Allocatable, Intent(InOut) :: atmnam(:)

    Integer, Allocatable :: displ(:), rec_size(:)

    Call assert(Allocated(atmnam), &
                'Gathered atmnam not allocated prior to calling gather_atomic_names')
    Call assert(Size(atmnam) == config%megatm, &
                'Gathered atmnam array does not contain natoms elements')

    !Different to mpi arrays for gathered%coords
    Allocate (rec_size(comm%mxnode), displ(comm%mxnode))
    rec_size = 0
    displ = 0
    Call gatherv_scatterv_index_arrays(comm, config%natms, rec_size, displ)
    atmnam = ''
    Call gallgatherv(comm, config%atmnam, rec_size, displ, atmnam)

  End Subroutine gather_atomic_names

  !> @brief Distribute forces from array of shape(3,megatm) to
  !! config, without using MPI
  !!
  !! Main useage: DFTB+ returns all forces on every process, therefore one isn't
  !! required to use MPI calls to distribute from the DFTB+ output to
  !! DLPOLY config object
  !!
  !! Unit conversion should be handled elsewhere
  !!
  !! @param[in]     comm        Object containing MPI settings and comms
  !! @param[in]     mpi_index   Object containing receive_count and displ
  !!                            index arrays
  !! @param[in]     forces      Forces returned by DFTB+, of size(3,megatm)
  !!
  !! @param[inout]  config      DLPOLY object containing local config data
  !!                            Returned with updated forces
  !
  Subroutine distribute_forces(comm, mpi_index, forces, config)
    Type(comms_type),            Intent(InOut) :: comm
    Type(mpi_distribution_type), Intent(In   ) :: mpi_index
    Real(wp), Allocatable,       Intent(In   ) :: forces(:, :)
    Type(configuration_type),    Intent(InOut) :: config

    Integer :: ia, ia_local, natms, offset, process_id

    !Zero DLPOLY forces from prior MD step
    natms = config%natms
    config%parts(1:natms)%fxx = 0._wp
    config%parts(1:natms)%fyy = 0._wp
    config%parts(1:natms)%fzz = 0._wp

    process_id = comm%idnode
    offset = Int(mpi_index%displ(process_id + 1) / 3)

    Do ia_local = 1, natms
      ia = offset + ia_local
      config%parts(ia_local)%fxx = forces(1, ia)
      config%parts(ia_local)%fyy = forces(2, ia)
      config%parts(ia_local)%fzz = forces(3, ia)
    End Do

  End Subroutine distribute_forces

End Module configuration
