Module z_density

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global z-density variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2014
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: comms_type,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: nzdndt
  Use errors_warnings, Only: error,&
                             info,&
                             error_alloc,&
                             error_dealloc
  Use kinds,           Only: wi,&
    wp,STR_LEN
  Use site,            Only: site_type

  Implicit None

  Private

  !> Type containing z density data
  Type, Public :: z_density_type
    Private

    !> Collection switch
    Logical, Public                    :: l_collect = .false.
    !> Printing switch
    Logical, Public                    :: l_print = .false.
    !> Number of configurations sampled
    Integer(Kind=wi), Public           :: n_samples = 0
    !> Collection frequency in steps
    Integer(Kind=wi), Public           :: frequency = 1
    !> Bin width
    Real(Kind=wp), Public              :: bin_width = 0.0_wp
    !> z density
    Real(Kind=wp), Allocatable, Public :: density(:, :)
    !> Maximum number of Zden grid points
    Integer(Kind=wi), Public           :: max_grid = 0

  Contains
    Private

    Procedure, Public :: init => allocate_z_density_arrays
    Final             :: cleanup
  End Type z_density_type

  Public ::  z_density_collect, z_density_compute

Contains

  Subroutine allocate_z_density_arrays(zdensity, mxatyp)
    Class(z_density_type), Intent(InOut) :: zdensity
    Integer(Kind=wi),      Intent(In   ) :: mxatyp

    Integer :: fail

    fail = 0

    Allocate (zdensity%density(1:zdensity%max_grid, 1:mxatyp), Stat=fail)

    If (fail > 0) Call error_alloc('zdensity%density', 'allocate_z_density_arrays')

    zdensity%density = 0.0_wp

  End Subroutine allocate_z_density_arrays

  Subroutine z_density_collect(zdensity, config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistic for z-density profile
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(configuration_type), Intent(InOut) :: config

    Integer       :: i, k, l
    Real(Kind=wp) :: rdelr, zlen, zleno2

    ! accumulator

    zdensity%n_samples = zdensity%n_samples + 1

    ! length of cell in z direction

    zlen = Abs(config%cell(3)) + Abs(config%cell(6)) + Abs(config%cell(9))

    ! half of z length

    zleno2 = zlen * 0.5_wp

    ! grid interval for density profiles

    rdelr = Real(zdensity%max_grid, wp) / zlen

    ! set up atom iatm type and accumulate statistic

    Do i = 1, config%natms
      k = config%ltype(i)

      l = Min(1 + Int((config%parts(i)%zzz + zleno2) * rdelr), zdensity%max_grid)

      zdensity%density(l, k) = zdensity%density(l, k) + 1.0_wp
    End Do

  End Subroutine z_density_collect

  Subroutine z_density_compute(config, zdensity, sites, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating z-density profile from
    ! accumulated data
    !
    ! copyright - daresbury laboratory
    ! author    - t.forester march 1994
    ! amended   - i.t.todorov march 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(site_type),          Intent(In   ) :: sites
    Type(comms_type),         Intent(InOut) :: comm

    Integer            :: j, k
    Real(Kind=wp)      :: delr, dvolz, factor, rho, rrr, sum, zlen

    ! open Z density file and write headers

    If (comm%idnode == 0) Then
      Open (Unit=nzdndt, File='ZDNDAT', Status='replace')
      Write (nzdndt, '(a)') config%cfgname
      Write (nzdndt, '(2i10)') sites%ntype_atom, zdensity%max_grid
    End If

    ! length of cell in z direction

    zlen = Abs(config%cell(3)) + Abs(config%cell(6)) + Abs(config%cell(9))

    ! grid interval for density profiles

    delr = zlen / Real(zdensity%max_grid, wp)

    ! volume of z-strip

    dvolz = (config%volm / zlen) * delr

    ! normalisation factor

    zdensity%n_samples = Max(zdensity%n_samples, 1)
    factor = 1.0_wp / (Real(zdensity%n_samples, wp) * dvolz)

    ! for every species

    Do k = 1, sites%ntype_atom
      If (comm%idnode == 0) Then
        Write (nzdndt, '(a8)') sites%unique_atom(k)
      End If

      ! global sum of data on all nodes

      Call gsum(comm, zdensity%density(1:zdensity%max_grid, k))

      ! running integration of z-density

      sum = 0.0_wp

      ! loop over distances

      Do j = 1, zdensity%max_grid
        rrr = (Real(j, wp) - 0.5_wp) * delr - zlen * 0.5_wp
        rho = zdensity%density(j, k) * factor
        sum = sum + rho * dvolz

        If (comm%idnode == 0) Then
          Write (nzdndt, "(1p,3e14.6)") rrr, rho, sum
        End If
      End Do
    End Do

    If (comm%idnode == 0) Close (Unit=nzdndt)

  End Subroutine z_density_compute

  Subroutine cleanup(zdensity)
    Type(z_density_type) :: zdensity
    Integer :: fail

    If (Allocated(zdensity%density)) Then
       Deallocate (zdensity%density, stat=fail)
       if (fail /= 0) call error_dealloc('zdensity%density', 'zdensity_cleanup')
    End If
  End Subroutine cleanup
End Module z_density
