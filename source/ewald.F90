Module ewald

  !!-----------------------------------------------------------------------
  !!
  !! dl_poly_4 module declaring ewald routines and arrays
  !!
  !! copyright - daresbury laboratory
  !! author    - i.t.todorov february 2015
  !! amended   - j.s.wilkins october  2018
  !!
  !!-----------------------------------------------------------------------

  Use bspline,         Only: bspline_type
  Use configuration,   Only: configuration_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc
  Use kinds,           Only: wp
  Use kspace,          Only: kspace_type
  Use spme,            Only: init_spme_data,&
                             spme_component
  Use vdw,             Only: VDW_12_6,&
                             VDW_BORN_HUGGINS_MEYER,&
                             VDW_BUCKINGHAM,&
                             VDW_HYDROGEN_BOND,&
                             VDW_LENNARD_JONES,&
                             VDW_MORSE_12,&
                             VDW_N_M,&
                             VDW_N_M_SHIFT,&
                             vdw_type

  Implicit None

  Public :: ewald_vdw_count, ewald_vdw_init, ewald_vdw_coeffs

  Private

  Type, Public :: ewald_type
    !! Base Ewald type containing data relevant to *all* ewald variants
    Private

    !> Ewald being used in this calculation
    Logical, Public :: active = .false.

    !> Ewald is setup?
    Logical, Public :: initialised = .false.

    !> Ewald is performing polynomial VdW parts
    Logical, Public :: vdw = .false.

    !> VdW sets
    Integer, Allocatable, Dimension(:), Public :: reduced_VdW

    !> Ewald real-space to be calculated direct not via tabulation
    Logical, Public :: direct = .false.

    !> FFT and KSpace info container
    Type(kspace_type), Public :: kspace

    !> Ewald convergence parameter or Coulomb damping parameter (A^-1)
    Real(Kind=wp), Public :: alpha

    ! Merged types here for simplicity
    ! Ewald type containing data relevant to SPME style ewald

    !> Number of potentials to handle
    Integer, Public :: num_pots = 0
    !>
    Logical, Public :: newjob_erf = .true.
    Logical, Public :: newjob_spme = .true.
    Logical, Public :: newjob_spme_gen = .true.

    !> SPME function container
    Type(spme_component), Dimension(:), Allocatable, Public :: spme_data
    !> Bspline container
    Type(bspline_type), Public :: bspline

  Contains

    Final :: deallocate_ewald_type

  End Type ewald_type

Contains

  Subroutine deallocate_ewald_type(ewld)
    Type(ewald_type) :: ewld
    Integer, Dimension(2) :: fail

    fail = 0
    If (Allocated(ewld%reduced_VdW)) Deallocate(ewld%reduced_VdW, stat=fail(1))
    If (Allocated(ewld%spme_data)) Deallocate(ewld%spme_data, stat=fail(2))

    If (Any(fail /= 0)) call error_dealloc('ewald_type', 'deallocate_ewald_type')

  End Subroutine deallocate_ewald_type

  Subroutine ewald_vdw_count(ewld, vdws)
    Type(ewald_type)              :: ewld
    Type(vdw_type), Intent(In   ) :: vdws

    Integer :: fail, ipot, ivdw, keypot

    ! Count number of variables needed (overallocate?)
    ewld%num_pots = 0

    countloop: Do ivdw = 1, vdws%n_vdw
      keypot = vdws%ltp(ivdw)

      ! Need to remove duplicates of the same pot for mapping to linear array
      Do ipot = 1, ewld%num_pots
        If (keypot == ewld%reduced_vdw(ipot)) Cycle countloop
      End Do

      Select Case (keypot)
      Case (VDW_12_6, VDW_LENNARD_JONES, VDW_N_M, VDW_N_M_SHIFT, VDW_HYDROGEN_BOND, VDW_BORN_HUGGINS_MEYER) ! 6-12
        ewld%num_pots = ewld%num_pots + 2
      Case (VDW_MORSE_12, VDW_BUCKINGHAM)
        ewld%num_pots = ewld%num_pots + 1
      Case default
        Call error(0, 'Ewald VdW potential requested but not possible')
      End Select
    End Do countloop

    Allocate (ewld%reduced_vdw(ewld%num_pots), stat=fail)
    If (fail > 0) Call error_alloc('ewld%reduced_vdw', 'two_body_forces')

    ! Reset counter for assignment
    ewld%num_pots = 0

    ! Need to know how many we're doing
    vdwloop: Do ivdw = 1, vdws%n_vdw
      keypot = vdws%ltp(ivdw)

      ! Need to remove duplicates of the same pot for mapping to linear array
      Do ipot = 1, ewld%num_pots
        If (keypot == ewld%reduced_vdw(ipot)) Cycle vdwloop
      End Do

      Select Case (keypot)
      Case (VDW_12_6, VDW_LENNARD_JONES) ! 6-12
        ewld%reduced_vdw(ewld%num_pots + 1:ewld%num_pots + 2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      Case (VDW_N_M, VDW_N_M_SHIFT) ! N-M
        ewld%reduced_vdw(ewld%num_pots + 1:ewld%num_pots + 2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      Case (VDW_BORN_HUGGINS_MEYER) ! 6-8
        ewld%reduced_vdw(ewld%num_pots + 1:ewld%num_pots + 2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      Case (VDW_HYDROGEN_BOND) ! 12-10
        ewld%reduced_vdw(ewld%num_pots + 1:ewld%num_pots + 2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      Case (VDW_MORSE_12) ! 12
        ewld%reduced_vdw(ewld%num_pots + 1:ewld%num_pots + 1) = keypot
        ewld%num_pots = ewld%num_pots + 1
      Case (VDW_BUCKINGHAM) ! 6
        ewld%reduced_vdw(ewld%num_pots + 1:ewld%num_pots + 1) = keypot
        ewld%num_pots = ewld%num_pots + 1
      Case default
        Call error(0, 'Ewald VdW potential requested but not possible')
      End Select
    End Do vdwloop

    If (ewld%num_pots == 0) &
      Call error(0, 'Ewald VdW potential requested, but not possible for given potentials')

  End Subroutine ewald_vdw_count

  Subroutine ewald_vdw_init(ewld)

    Type(ewald_type)              :: ewld

    Integer :: ipot, keypot
    Logical :: skip

    skip = .false.

    Do ipot = 1, ewld%num_pots
      keypot = ewld%reduced_vdw(ipot)

      If (skip) Then
        skip = .false.
        Cycle
      End If

      Select Case (keypot)

      Case (VDW_12_6, VDW_LENNARD_JONES) ! 12-6
        Call init_spme_data(ewld%spme_data(ipot), 12)
        Call init_spme_data(ewld%spme_data(ipot + 1), 6)

        skip = .true.

      Case (VDW_N_M, VDW_N_M_SHIFT) ! N-M
        ! Initialised in ewald_vdw_coeffs

        skip = .true.

      Case (VDW_HYDROGEN_BOND) ! 12-10
        Call init_spme_data(ewld%spme_data(ipot), 12)
        Call init_spme_data(ewld%spme_data(ipot + 1), 10)

        skip = .true.

      Case (VDW_BORN_HUGGINS_MEYER) ! 6-8
        Call error(0, 'Ewald VdW N-M potential requested but not possible (contains exponential)')
        Call init_spme_data(ewld%spme_data(ipot), 6)
        Call init_spme_data(ewld%spme_data(ipot + 1), 8)

        skip = .true.

      Case (VDW_MORSE_12) ! 12
        Call error(0, 'Ewald VdW N-M potential requested but not possible (contains exponential)')
        ewld%spme_data(ipot)%scaling = 1.0_wp
        Call init_spme_data(ewld%spme_data(ipot), 12)

      Case (VDW_BUCKINGHAM) ! -6
        Call error(0, 'Ewald VdW N-M potential requested but not possible (contains exponential)')
        ewld%spme_data(ipot)%scaling = -1.0_wp
        Call init_spme_data(ewld%spme_data(ipot), 6)

      Case default
        Call error(0, 'Ewald VdW potential requested but not possible')
      End Select
    End Do

  End Subroutine ewald_vdw_init

  Subroutine ewald_vdw_coeffs(config, vdws, ewld, vdw_coeffs)

    Type(configuration_type),                    Intent(In   ) :: config
    Type(vdw_type),                              Intent(In   ) :: vdws
    Type(ewald_type),                            Intent(inout) :: ewld
    Real(kind=wp), Allocatable, Dimension(:, :), Intent(  Out) :: vdw_coeffs

    Integer :: fail, ipot, j, k, keypot

    ! Set up potentials

    Allocate (vdw_coeffs(config%mxatms, ewld%num_pots), stat=fail)
    If (fail > 0) Call error_alloc('vdw_coeffs', 'two_body_forces')

    !! JW952
    ! Build coeffs array ( Assume sqrt(i)*sqrt(j) division [see: Darden])
    ! Inefficient, but temporary?
    Do k = 1, config%mxatms
      j = config%ltype(k)
      If (j == 0) Cycle
      j = j * (j - 1) / 2 + j
      j = vdws%list(j)
      keypot = vdws%ltp(j)

      ! Look up pot
      Do ipot = 1, ewld%num_pots
        If (keypot == ewld%reduced_vdw(ipot)) Exit
      End Do
      If (ipot > ewld%num_pots) Call error(0, 'Error in potentials mapping')

      ! Assign params to this atom
      Select Case (keypot)

      Case (VDW_LENNARD_JONES) ! 12-6
        ! 4eps * (A)
        ewld%spme_data(ipot)%scaling = vdws%param(1, j) * 4.0_wp
        ewld%spme_data(ipot + 1)%scaling = -vdws%param(1, j) * 4.0_wp
        vdw_coeffs(k, ipot) = vdws%param(2, j)**6
        vdw_coeffs(k, ipot + 1) = vdws%param(2, j)**3

      Case (VDW_12_6) ! 12-6
        ewld%spme_data(ipot)%scaling = 1.0_wp
        ewld%spme_data(ipot)%scaling = -1.0_wp
        vdw_coeffs(k, ipot) = vdws%param(1, j)
        vdw_coeffs(k, ipot + 1) = vdws%param(2, j)

      Case (VDW_N_M, VDW_N_M_SHIFT) ! N-M
        If (.not. ewld%spme_data(ipot)%initialised) Then
          Call init_spme_data(ewld%spme_data(ipot), Int(vdws%param(2, j)))
          Call init_spme_data(ewld%spme_data(ipot + 1), Int(vdws%param(3, j)))
        End If

        ! Prefac
        ewld%spme_data(ipot)%scaling = vdws%param(1, j) / (vdws%param(2, j) - vdws%param(3, j))

        ewld%spme_data(ipot + 1)%scaling = -vdws%param(2, j) * ewld%spme_data(ipot)%scaling
        ewld%spme_data(ipot)%scaling = vdws%param(3, j) * ewld%spme_data(ipot)%scaling

        vdw_coeffs(k, ipot) = vdws%param(4, j)
        vdw_coeffs(k, ipot + 1) = vdws%param(4, j)

      Case (VDW_BORN_HUGGINS_MEYER) ! 6-8
        vdw_coeffs(k, ipot) = -vdws%param(4, j)
        vdw_coeffs(k, ipot + 1) = -vdws%param(5, j)

      Case (VDW_HYDROGEN_BOND) ! 12-10
        ewld%spme_data(ipot)%scaling = 1.0_wp
        ewld%spme_data(ipot)%scaling = -1.0_wp
        vdw_coeffs(k, ipot) = vdws%param(1, j)
        vdw_coeffs(k, ipot + 1) = vdws%param(2, j)

      Case (VDW_MORSE_12) ! 12
        vdw_coeffs(k, ipot) = vdws%param(1, j)

      Case (VDW_BUCKINGHAM) ! -6
        vdw_coeffs(k, ipot) = vdws%param(1, j)

      Case default
        Call error(0, 'Ewald VdW potential requested but not possible')
      End Select
    End Do

  End Subroutine ewald_vdw_coeffs

End Module ewald
