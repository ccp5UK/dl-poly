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

  Use bspline,         Only : bspline_type
  Use comms,           Only : comms_type
  Use configuration,   Only : configuration_type
  Use domains,         Only : domains_type
  Use errors_warnings, Only : error, error_alloc, error_dealloc
  Use kinds,           Only : wp
  Use kspace,          Only : kspace_type
  Use constants,       Only : twopi
  Use spme,            Only : spme_component
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
    Type( kspace_type ), Public :: kspace

    !> Ewald convergence parameter or Coulomb damping parameter (A^-1)
    Real( Kind = wp ), Public :: alpha

    Logical, Public :: newjob_two_body=.True.
    Logical, Public :: newjob_spme_init=.True.

    ! Merged types here for simplicity
    ! Ewald type containing data relevant to SPME style ewald

    !> Number of potentials to handle
    Integer,                                Public :: num_pots = 0
    !> SPME function container
    Type( spme_component ), Dimension( : ), Allocatable, Public :: spme_data
    !> Bspline container
    Type( bspline_type ),                                Public :: bspline

  End Type ewald_type

contains

  subroutine ewald_vdw_count(ewld, vdws)
    use vdw, only : vdw_type, VDW_NULL, VDW_TAB, VDW_12_6, VDW_LENNARD_JONES, VDW_N_M, &
      VDW_BUCKINGHAM, VDW_BORN_HUGGINS_MEYER, VDW_HYDROGEN_BOND, &
      VDW_N_M_SHIFT, VDW_MORSE, VDW_WCA, VDW_DPD, VDW_AMOEBA, &
      VDW_LENNARD_JONES_COHESIVE, VDW_MORSE_12, VDW_RYDBERG, VDW_ZBL, &
      VDW_ZBL_SWITCH_MORSE, VDW_ZBL_SWITCH_BUCKINGHAM
    type( ewald_type ) :: ewld
    type( vdw_type ), intent( in    ) :: vdws
    integer :: keypot
    integer :: ivdw, ipot

    integer :: fail

    ! Count number of variables needed (overallocate?)
    ewld%num_pots = 0

    countloop:do ivdw=1,vdws%n_vdw
      keypot=vdws%ltp(ivdw)

      ! Need to remove duplicates of the same pot for mapping to linear array
      do ipot = 1, ewld%num_pots
        if (keypot == ewld%reduced_vdw(ipot)) cycle countloop
      end do

      select case (keypot)
      case (VDW_12_6, VDW_LENNARD_JONES, VDW_N_M,VDW_N_M_SHIFT, VDW_HYDROGEN_BOND, VDW_BORN_HUGGINS_MEYER) ! 6-12
        ewld%num_pots = ewld%num_pots + 2
      case (VDW_MORSE_12, VDW_BUCKINGHAM)
        ewld%num_pots = ewld%num_pots + 1
      case default
        call error(0,'Ewald VdW potential requested but not possible')
      end select
    end do countloop


    allocate(ewld%reduced_vdw(ewld%num_pots), stat=fail)
    if (fail>0) call error_alloc('ewld%reduced_vdw','two_body_forces')

    ! Reset counter for assignment
    ewld%num_pots = 0

    ! Need to know how many we're doing
    vdwloop:do ivdw=1,vdws%n_vdw
      keypot=vdws%ltp(ivdw)

      ! Need to remove duplicates of the same pot for mapping to linear array
      do ipot = 1, ewld%num_pots
        if (keypot == ewld%reduced_vdw(ipot)) cycle vdwloop
      end do

      select case (keypot)
      case (VDW_12_6, VDW_LENNARD_JONES) ! 6-12
        ewld%reduced_vdw(ewld%num_pots+1:ewld%num_pots+2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      case (VDW_N_M,VDW_N_M_SHIFT) ! N-M
        ewld%reduced_vdw(ewld%num_pots+1:ewld%num_pots+2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      case (VDW_BORN_HUGGINS_MEYER) ! 6-8
        ewld%reduced_vdw(ewld%num_pots+1:ewld%num_pots+2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      case (VDW_HYDROGEN_BOND) ! 12-10
        ewld%reduced_vdw(ewld%num_pots+1:ewld%num_pots+2) = keypot
        ewld%num_pots = ewld%num_pots + 2
      case (VDW_MORSE_12) ! 12
        ewld%reduced_vdw(ewld%num_pots+1:ewld%num_pots+1) = keypot
        ewld%num_pots = ewld%num_pots + 1
      case (VDW_BUCKINGHAM) ! 6
        ewld%reduced_vdw(ewld%num_pots+1:ewld%num_pots+1) = keypot
        ewld%num_pots = ewld%num_pots + 1
      case default
        call error(0,'Ewald VdW potential requested but not possible')
      end select
    end do vdwloop

    if (ewld%num_pots == 0) &
      call error(0,'Ewald VdW potential requested, but not possible for given potentials')

  end subroutine ewald_vdw_count

  subroutine ewald_vdw_init(ewld, vdws)

    use vdw, only : vdw_type, VDW_NULL, VDW_TAB, VDW_12_6, VDW_LENNARD_JONES, VDW_N_M, &
      VDW_BUCKINGHAM, VDW_BORN_HUGGINS_MEYER, VDW_HYDROGEN_BOND, &
      VDW_N_M_SHIFT, VDW_MORSE, VDW_WCA, VDW_DPD, VDW_AMOEBA, &
      VDW_LENNARD_JONES_COHESIVE, VDW_MORSE_12, VDW_RYDBERG, VDW_ZBL, &
      VDW_ZBL_SWITCH_MORSE, VDW_ZBL_SWITCH_BUCKINGHAM
    use spme, only : init_spme_data

    type ( ewald_type ) :: ewld
    type ( vdw_type ), intent( in    ) :: vdws
    Integer :: keypot
    Integer :: ipot
    Logical :: skip

    skip = .false.

    do ipot = 1, ewld%num_pots
      keypot = ewld%reduced_vdw(ipot)

      if ( skip ) then
        skip = .false.
        cycle
      end if

      select case (keypot)

      case (VDW_12_6, VDW_LENNARD_JONES) ! 12-6
        call init_spme_data( ewld%spme_data(ipot), 12 )
        call init_spme_data( ewld%spme_data(ipot+1), 6 )

        skip = .true.

      case (VDW_N_M,VDW_N_M_SHIFT) ! N-M
        ! Initialised in ewald_vdw_coeffs

        skip = .true.

      case (VDW_HYDROGEN_BOND) ! 12-10
        call init_spme_data( ewld%spme_data(ipot), 12 )
        call init_spme_data( ewld%spme_data(ipot+1), 10 )

        skip = .true.

      case (VDW_BORN_HUGGINS_MEYER) ! 6-8
        call error(0,'Ewald VdW N-M potential requested but not possible (contains exponential)')
        call init_spme_data( ewld%spme_data(ipot), 6 )
        call init_spme_data( ewld%spme_data(ipot+1), 8 )

        skip = .true.

      case (VDW_MORSE_12) ! 12
        call error(0,'Ewald VdW N-M potential requested but not possible (contains exponential)')
        ewld%spme_data(ipot)%scaling = 1.0_wp
        call init_spme_data( ewld%spme_data(ipot), 12 )

      case (VDW_BUCKINGHAM) ! -6
        call error(0,'Ewald VdW N-M potential requested but not possible (contains exponential)')
        ewld%spme_data(ipot)%scaling = -1.0_wp
        call init_spme_data( ewld%spme_data(ipot), 6 )

      case default
        call error(0,'Ewald VdW potential requested but not possible')
      end select
    end do

  end subroutine ewald_vdw_init

  subroutine ewald_vdw_coeffs(config, vdws, ewld, vdw_coeffs )

    use vdw, only : vdw_type, VDW_NULL, VDW_TAB, VDW_12_6, VDW_LENNARD_JONES, VDW_N_M, &
      VDW_BUCKINGHAM, VDW_BORN_HUGGINS_MEYER, VDW_HYDROGEN_BOND, &
      VDW_N_M_SHIFT, VDW_MORSE, VDW_WCA, VDW_DPD, VDW_AMOEBA, &
      VDW_LENNARD_JONES_COHESIVE, VDW_MORSE_12, VDW_RYDBERG, VDW_ZBL, &
      VDW_ZBL_SWITCH_MORSE, VDW_ZBL_SWITCH_BUCKINGHAM
    use spme, only : init_spme_data

    type ( configuration_type ), intent ( in    ) :: config
    type ( ewald_type ),    intent ( inout ) :: ewld
    type ( vdw_type ),           intent ( in    ) :: vdws
    real ( kind = wp ), dimension( :,: ), allocatable, intent (   out ) :: vdw_coeffs

    integer :: keypot
    integer :: j, k, ipot
    integer :: fail

    ! Set up potentials

    allocate(vdw_coeffs(config%mxatms,ewld%num_pots), stat =fail)
    if ( fail > 0 ) call error_alloc('vdw_coeffs','two_body_forces')

    !! JW952
    ! Build coeffs array ( Assume sqrt(i)*sqrt(j) division [see: Darden])
    ! Inefficient, but temporary?
    do k = 1, config%mxatms
      j = config%ltype(k)
      if (j == 0) cycle
      j = j*(j-1)/2 + j
      j = vdws%list(j)
      keypot = vdws%ltp(j)

      ! Look up pot
      do ipot = 1, ewld%num_pots
        if ( keypot == ewld%reduced_vdw(ipot) ) exit
      end do
      if ( ipot > ewld%num_pots ) call error(0,'Error in potentials mapping')

      ! Assign params to this atom
      select case (keypot)

      case (VDW_LENNARD_JONES) ! 12-6
        ! 4eps * (A)
        ewld%spme_data(ipot)%scaling = vdws%param(1,j) * 4.0_wp
        ewld%spme_data(ipot+1)%scaling = -vdws%param(1,j) * 4.0_wp
        vdw_coeffs(k,ipot) = vdws%param(2,j)**6
        vdw_coeffs(k,ipot+1) = vdws%param(2,j)**3

      case (VDW_12_6) ! 12-6
        ewld%spme_data(ipot)%scaling = 1.0_wp
        ewld%spme_data(ipot)%scaling = -1.0_wp
        vdw_coeffs(k,ipot) = vdws%param(1,j)
        vdw_coeffs(k,ipot+1) = vdws%param(2,j)

      case (VDW_N_M,VDW_N_M_SHIFT) ! N-M
        if (.not. ewld%spme_data(ipot)%initialised) then
          call init_spme_data( ewld%spme_data(ipot), int(vdws%param(2,j)) )
          call init_spme_data( ewld%spme_data(ipot+1), int(vdws%param(3,j)) )
        end if

        ! Prefac
        ewld%spme_data(ipot)%scaling = vdws%param(1,j) / (vdws%param(2,j) - vdws%param(3,j))

        ewld%spme_data(ipot+1)%scaling = -vdws%param(2,j) * ewld%spme_data(ipot)%scaling
        ewld%spme_data(ipot)%scaling = vdws%param(3,j) * ewld%spme_data(ipot)%scaling

        vdw_coeffs(k,ipot) = vdws%param(4,j)
        vdw_coeffs(k,ipot+1) = vdws%param(4,j)

      case (VDW_BORN_HUGGINS_MEYER) ! 6-8
        vdw_coeffs(k,ipot) = -vdws%param(4,j)
        vdw_coeffs(k,ipot+1) = -vdws%param(5,j)

      case (VDW_HYDROGEN_BOND) ! 12-10
        ewld%spme_data(ipot)%scaling = 1.0_wp
        ewld%spme_data(ipot)%scaling = -1.0_wp
        vdw_coeffs(k,ipot) = vdws%param(1,j)
        vdw_coeffs(k,ipot+1) = vdws%param(2,j)

      case (VDW_MORSE_12) ! 12
        vdw_coeffs(k,ipot) = vdws%param(1,j)

      case (VDW_BUCKINGHAM) ! -6
        vdw_coeffs(k,ipot) = vdws%param(1,j)

      case default
        call error(0,'Ewald VdW potential requested but not possible')
      end select
    end do

  end subroutine ewald_vdw_coeffs

End Module ewald
