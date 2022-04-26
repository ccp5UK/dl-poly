Module units
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle unit conversion in dlpoly
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins april 2020
  !! contrib - a.m.elena april 2022 - add more units schemes
  !!-----------------------------------------------------------------------

  Use kinds,           Only: wp,STR_LEN
  Use constants,       Only: pi, boltz
  Use hash,            Only: hash_table, MAX_KEY
  Use hash,            Only: get_double, get_int, get_complex
  Use errors_warnings, Only: error, error_units, error_alloc, error_dealloc
  Use parse,           Only: lower_case
  Implicit None

  Private

  Integer, Parameter :: ndims = 7
  Character, Dimension(ndims), parameter :: dims = ['M','L','T','t','C','l','a'] ! mass length time temp current luminosity angle

  Integer, Dimension(ndims), Parameter :: &
       dims_charge = [0,0,1,0,1,0,0], &
       dims_energy = [1,2,-2,0,0,0,0], &
       dims_pressure = [1,-1,-2,0,0,0,0], &
       dims_force = [1,1,-2,0,0,0,0], &
       dims_velocity = [0,1,-1,0,0,0,0], &
       dims_power = [1,2,-3,0,0,0,0], &
       dims_surf_ten = [1,0,-2,0,0,0,0], &
       dims_emf = [1,2,-3,0,-1,0,0]

  Type, Private, Extends(hash_table) :: units_hash_table
  Contains
    Generic  , Public  :: get => get_unit
    Procedure, Private :: get_unit
  end type units_hash_table

  Type, Private :: unit_data
    !! Type containing data corresponding to units
    Character(Len=STR_LEN) :: name
    Character(Len=MAX_KEY) :: abbrev
    Real(Kind=wp) :: conversion_to_internal
    Integer, Dimension(ndims) :: dims = [0, 0, 0, 0, 0, 0, 0] ! mass length time temp current luminosity angle
  contains
    Generic, Public :: Operator(*) => unit_mult
    Generic, Public :: Operator(/) => unit_div
    Generic, Public :: Operator(**) => unit_pow
    Procedure, Pass, Private :: unit_mult, unit_div, unit_pow
    Procedure, Nopass :: init => init_unit
  end type unit_data

  Type(unit_data), Private, Parameter :: null_unit = unit_data('', '', 1.0_wp)


  Type, Public :: units_scheme
    !! Basic encapsulation of output units scheme
    Character(Len=STR_LEN), Public :: length, time, mass, charge, energy, &
      temp='K', current='A', luminosity='', angle='rad', &
      pressure, force, velocity, power, surf_ten, emf
    Type(unit_data), Public :: length_unit=null_unit, time_unit=null_unit, &
      mass_unit=null_unit, charge_unit=null_unit, energy_unit=null_unit, &
      luminosity_unit=null_unit, angle_unit=null_unit, pressure_unit=null_unit, &
      force_unit=null_unit, velocity_unit =null_unit, power_unit=null_unit, &
      temp_unit=null_unit, current_unit=null_unit, surf_ten_unit=null_unit,&
      emf_unit=null_unit
  Contains
    Procedure :: test => test_scheme
    Procedure :: recalc => recalculate_scheme_units
    Procedure, Nopass :: init => init_scheme
  End Type units_scheme

  Type(units_scheme), Public, Protected, Save :: internal_units = units_scheme( &
    length = 'internal_l', &
    time = 'internal_t', &
    mass = 'internal_m', &
    charge = 'internal_q', &
    energy = 'internal_e', &
    temp = 'K', &
    current = 'A', &
    luminosity = '', &
    angle = 'rad', &
    pressure = 'internal_p', &
    force = 'internal_f', &
    velocity = 'internal_v', &
    power = 'internal_e/internal_t', &
    surf_ten = 'internal_f/internal_l', &
    emf = 'internal_e/internal_q')

  Type(units_scheme), Public, Protected, Save :: si_units = units_scheme( &
    length = 'm', &
    time = 's', &
    mass = 'kg', &
    charge = 'C', &
    energy = 'J', &
    pressure = 'Pa', &
    force = 'N', &
    velocity = 'm/s', &
    power = 'W', &
    surf_ten = 'N/m', &
    emf = 'V')

  Type(units_scheme), Public, Protected, Save :: atomic_units = units_scheme( &
    length = 'ang', &
    time = 'ps', &
    mass = 'amu', &
    charge = 'q_e', &
    energy = 'e.V', &
    pressure = 'GPa', &
    force = 'e.V/ang', &
    velocity = 'ang/ps', &
    power = 'e.V/ps', &
    surf_ten = 'e.V/ang^2', &
    emf = 'e.V/q_e')

  Type(units_scheme), Public, Protected, Save :: hartree_units = units_scheme( &
    length = 'bohr', &
    time = 'aut', &
    mass = 'm_e', &
    charge = 'q_e', &
    energy = 'Ha', &
    pressure = 'Ha/bohr^3', &
    force = 'Ha/bohr', &
    velocity = 'auv', &
    power = 'Ha/aut', &
    surf_ten = 'Ha/bohr^2', &
    emf = 'Ha/q_e')

  Type(units_scheme), Public, Protected, Save :: kj_units = units_scheme( &
    length = 'ang', &
    time = 'ps', &
    mass = 'amu', &
    charge = 'q_e', &
    energy = 'kJ', &
    pressure = 'GPa', &
    force = 'kJ/ang', &
    velocity = 'ang/ps', &
    power = 'kJ/ps', &
    surf_ten = 'kJ/ang^2', &
    emf = 'kJ/q_e')

  Type(units_scheme), Public, Protected, Save :: kcal_units = units_scheme( &
    length = 'ang', &
    time = 'ps', &
    mass = 'amu', &
    charge = 'q_e', &
    energy = 'kCal', &
    pressure = 'GPa', &
    force = 'kCal/ang', &
    velocity = 'ang/ps', &
    power = 'kCal/ps', &
    surf_ten = 'Cal/ang^2', &
    emf = 'kCal/q_e')

  Type(units_scheme), Public, Protected, Save :: kb_units = units_scheme( &
    length = 'ang', &
    time = 'ps', &
    mass = 'amu', &
    charge = 'q_e', &
    energy = 'e.V', &
    pressure = 'GPa', &
    force = 'e.V/ang', &
    velocity = 'ang/ps', &
    power = 'e.V/ps', &
    surf_ten = 'e.V/ang^2', &
    emf = 'e.V/q_e')

  Type(units_scheme), Public, Protected, Save :: out_units

  Type(units_hash_table), Private, save :: units_table

  Interface convert_units
    module procedure convert_units_str_str
    module procedure convert_units_unit_unit
  end Interface convert_units

  Public :: initialise_units
  Public :: destroy_units
  Public :: convert_units
  Public :: set_timestep

  Public :: init_scheme
  Public :: set_out_units
  Public :: to_out_units


contains

  Subroutine initialise_units()
    !!-----------------------------------------------------------------------
    !!
    !! Initialise units table with all units
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Real(kind=wp), Parameter :: planck_internal = 6.350780668_wp
    Real(kind=wp), Parameter :: electron_charge = 1.0_wp
    Real(kind=wp), Parameter :: coulomb = 6.241509074460763e+18_wp
    Real(kind=wp), Parameter :: avogadro = 6.022140857e23_wp

    Real(kind=wp), Parameter :: metre = 1e10_wp
    Real(kind=wp), Parameter :: angstrom = 1.0_wp
    Real(kind=wp), Parameter :: bohr = 0.52918_wp*angstrom
    Real(kind=wp), Parameter :: inch = 2.54e-2_wp*metre

    Real(kind=wp), Parameter :: joulepmol = 0.1_wp, joule = joulepmol !*avogadro
    Real(kind=wp), Parameter :: calorie = 4.1842_wp*joule

    Real(kind=wp), Parameter :: hartree = 4.359744722e-18_wp*joule

    Real(kind=wp), Parameter :: kilogram = 6.0229552894949e+26_wp
    Real(kind=wp), Parameter :: electron_mass = 9.1093837015e-31_wp*kilogram
    Real(kind=wp), Parameter :: pound = 0.45359237_wp*kilogram

    Real(kind=wp), Parameter :: second = 1e12_wp
    Real(kind=wp), Parameter :: aut = 2.4188843265857e-17_wp*second

    Real(kind=wp), Parameter :: atmosphere = 1.0_wp / 163.882576_wp

    Real(kind=wp), Parameter :: newton = kilogram*metre/second**2
    Real(kind=wp), Parameter :: pascal = newton/metre**2

    Real(kind=wp), Parameter :: gravity = 9.81_wp*metre/second**2

    Real(kind=wp), Parameter :: degree = pi / 180.0_wp
    Real(kind=wp), Parameter :: gradian = pi / 200.0_wp

    call units_table%init(100)

    ! Time

    call units_table%set("internal_t", init_unit(abbrev="internal_t", name="Picosecond", time=1, to_internal=1.0_wp))
    call units_table%set("hr", init_unit(abbrev="hr", name="Hour", time=1, to_internal=3600.0_wp*second))
    call units_table%set("min", init_unit(abbrev="min", name="Minute", time=1, to_internal=60.0_wp*second))
    call units_table%set("s", init_unit(abbrev="s", name="Second", time=1, to_internal=second))
    call units_table%set("aut", init_unit(abbrev="aut", name="Atomic Time Unit", time=1, to_internal=aut))

    ! Length

    call units_table%set("internal_l", init_unit(abbrev="internal_l", name="Angstrom", length=1, to_internal=1.0_wp))
    call units_table%set("ang", init_unit(abbrev="ang", name="Angstrom", length=1, to_internal=angstrom))
    call units_table%set("bohr", init_unit(abbrev="bohr", name="Bohr", length=1, to_internal=bohr))
    call units_table%set("m", init_unit(abbrev="m", name="Metre", length=1, to_internal=metre))
    call units_table%set('in', init_unit(abbrev="in", name="Inch", length=1, to_internal=inch))
    call units_table%set("ft", init_unit(abbrev="ft", name="Foot", length=1, to_internal=12.0_wp*inch))

    ! Mass

    call units_table%set("internal_m", init_unit(abbrev="internal_m", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("da", init_unit(abbrev="Da", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("amu", init_unit(abbrev="amu", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("m_e", init_unit(abbrev="m_e", name="Electron mass", mass=1, to_internal=electron_mass))
    call units_table%set("g", init_unit(abbrev="g", name="Gram", mass=1, to_internal=1.0e-3_wp*kilogram))
    call units_table%set("lb", init_unit(abbrev="lb", name="Pound", mass=1, to_internal=pound))
    call units_table%set("oz", init_unit(abbrev="oz", name="Ounce", mass=1, to_internal=pound/16.0_wp))

    ! Charge

    call units_table%set("internal_q", &
         & init_unit(abbrev="internal_q", name="Elementary charge", current=1, time=1, to_internal=1.0_wp))
    call units_table%set("q_e", &
         & init_unit(abbrev="q_e", name="Elementary charge", current=1, time=1, to_internal=electron_charge))
    call units_table%set("e", &
         & init_unit(abbrev="e", name="Elementary charge", current=1, time=1, to_internal=electron_charge))
    call units_table%set("c", &
         & init_unit(abbrev="C", name="Coulomb", current=1, time=1, to_internal=coulomb))

    ! Current

    call units_table%set("a", &
         & init_unit(abbrev="A", name="Ampere", current=1, to_internal=electron_charge/second))

    ! Energy

    call units_table%set("internal_e", &
         & init_unit(abbrev="internal_e", name="10 J/mol", mass=1, length=2, time=-2, to_internal=1.0_wp))
    call units_table%set("j", init_unit(abbrev="J", name="Joule", mass=1, length=2, time=-2, to_internal=joule))
    call units_table%set("cal", &
         & init_unit(abbrev="Cal", name="Calorie", mass=1, length=2, time=-2, to_internal=calorie))
    call units_table%set("ha", init_unit(abbrev="Ha", name="Hartree", &
         & mass=1, length=2, time=-2, to_internal=hartree))
    call units_table%set("e_h", init_unit(abbrev="Ha", name="Hartree", &
         & mass=1, length=2, time=-2, to_internal=hartree))
    call units_table%set("ry", init_unit(abbrev="Ry", name="Rydberg", &
         & mass=1, length=2, time=-2, to_internal=0.5_wp*hartree))

    ! Temperature

    call units_table%set("k", init_unit(abbrev="K", name="Kelvin", temp=1, to_internal=1.0_wp))

    ! Pressure

    call units_table%set("internal_p", &
         & init_unit(abbrev="internal_p", name="163 atm", mass=1, length=-1, time=-2, to_internal=1.0_wp))
    call units_table%set("atm", init_unit(abbrev="atm", name="Atmosphere", mass=1, length=-1, time=-2, to_internal=atmosphere))
    call units_table%set("pa", init_unit(abbrev="Pa", name="Pascal", mass=1, length=-1, time=-2, to_internal=pascal))

    ! Force

    call units_table%set("internal_f", &
         & init_unit(abbrev="internal_f", name="internal_f", mass=1, length=1, time=-2, to_internal=1.0_wp))
    call units_table%set("n", &
         & init_unit(abbrev="N", name="Newton", mass=1, length=1, time=-2, to_internal=newton))
    call units_table%set("dyn", &
         & init_unit(abbrev="dyn", name="Dyne", mass=1, length=1, time=-2, to_internal=1e-5_wp*newton))

    ! Velocity

    call units_table%set("internal_v", &
         & init_unit(abbrev="ang/ps", name="Angstrom per picosecond", length=1, time=-1, to_internal=1.0_wp))
    call units_table%set("auv", &
         & init_unit(abbrev="auv", name="Atomic Velocity Unit", length=1, time=-1, to_internal=bohr/aut))

    ! Constants

    call units_table%set("grav", &
         & init_unit(abbrev="g_e", name="9.81 m/s^2", length=1, time=-2, to_internal=gravity))
    call units_table%set("k_b", &
         & init_unit(abbrev="k_B", name="Boltzmann constant", length=2, mass=1, time=-2, temp=-1, to_internal=boltz))

    ! Power

    call units_table%set("w", init_unit(abbrev="W", name="Watt", mass=1, length=2, time=-3, to_internal=joule/second))

    ! Voltage

    call units_table%set("v", init_unit(abbrev="V", name="Volt", mass=1, length=2, time=-3, current=-1, &
         & to_internal=joule/coulomb))


    ! Angles

    call units_table%set("rad", init_unit(abbrev="rad", name="Radian", angle=1, to_internal=1.0_wp))
    call units_table%set("deg", init_unit(abbrev="deg", name="Degree", angle=1, to_internal=degree))
    call units_table%set("grad", init_unit(abbrev="grad", name="Gradian", angle=1, to_internal=gradian))

    ! Unitless

    call units_table%set("%", init_unit(abbrev="%", name="%", to_internal=0.01_wp))
    call units_table%set("", init_unit(abbrev="", name="", to_internal=1.0_wp))
    call units_table%set("mol", init_unit(abbrev="mol", name="Mole", to_internal=avogadro))

    ! System Properties

    call units_table%set("temp", init_unit(abbrev="T", name="System Temperature", temp=1, to_internal=1.0_wp))

    call internal_units%recalc()
    call set_out_units(internal_units)

  End Subroutine initialise_units

  Subroutine destroy_units()
    !!-----------------------------------------------------------------------
    !! Deallocate units table
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    call units_table%destroy()
  end Subroutine destroy_units

  Function convert_units_str_str(val, from, to, stat) result(res)
    !!-----------------------------------------------------------------------
    !!
    !! Convert val amount of unit "from" into unit "to"
    !!
    !! On failure --
    !!
    !! if stat present
    !! Set output name to error, stat to true
    !! Else
    !! Error out
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Real(kind=wp), Intent(in) :: val
    Character(Len=*), Intent(In) :: from, to
    Logical, Intent(Out), Optional :: stat
    Real(kind=wp) :: res
    Integer :: i
    Type( unit_data ) :: from_unit, to_unit
    Type( unit_data ) :: output

    if (present(stat)) then
      stat = .true.
    end if

    from_unit = parse_unit_string(from)
    to_unit = parse_unit_string(to)
    output = to_unit / from_unit

    if (any(output%dims /= 0)) then

#ifdef debug
      do i = 1, ndims
        if (from_unit%dims(i) /= 0) write(0, "(A2,i0,'.')", advance='No') dims(i)//"^",from_unit%dims(i)
      end do
      write(0,*)
      do i = 1, ndims
        if (to_unit%dims(i) /= 0) write(0, "(A2,i0,'.')", advance='No') dims(i)//"^",to_unit%dims(i)
      end do
#endif

      if (present(stat)) then
        stat=.false.
        return
      else
        call error_units(from, to)
      end if

    end if

    res = val / output%conversion_to_internal

  end Function convert_units_str_str

  Function convert_units_unit_unit(val, from_unit, to_unit, stat) result(res)
    !!-----------------------------------------------------------------------
    !!
    !! Convert val amount of unit "from" into unit "to"
    !!
    !! On failure --
    !!
    !! if stat present
    !! Set output name to error, stat to true
    !! Else
    !! Error out
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Real(kind=wp), Intent(in) :: val
    Type( unit_data ),Intent(In) :: from_unit, to_unit
    Logical, Intent(Out), Optional :: stat
    Real(kind=wp) :: res
    Integer :: i
    Type( unit_data ) :: output

    if (present(stat)) then
      stat = .true.
    end if

    output = to_unit / from_unit

    if (any(output%dims /= 0)) then

#ifdef debug
      do i = 1, ndims
        if (from_unit%dims(i) /= 0) write(0, "(A2,i0,'.')", advance='No') dims(i)//"^",from_unit%dims(i)
      end do
      write(0,*)
      do i = 1, ndims
        if (to_unit%dims(i) /= 0) write(0, "(A2,i0,'.')", advance='No') dims(i)//"^",to_unit%dims(i)
      end do
#endif

      if (present(stat)) then
        stat=.false.
        return
      else
        call error_units(from_unit%name, to_unit%name)
      end if

    end if

    res = val / output%conversion_to_internal

  end Function convert_units_unit_unit

  Function init_scheme(length, time, mass, charge, energy, &
       temp, current, luminosity, angle, &
       pressure, force, velocity, power, surf_ten, emf) Result(scheme)
    !!-----------------------------------------------------------------------
    !! Initialise a unit scheme, if not otherwise specified it will default to
    !! current out_units
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2022
    !!-----------------------------------------------------------------------

    Type(units_scheme) :: scheme
    Character(Len=*), Intent(In), Optional :: length, time, mass, charge, energy, &
         temp, current, luminosity, angle, &
         pressure, force, velocity, power, surf_ten, emf

    if (present(length)) then
      scheme%length = length
    else
      scheme%length = out_units%length
    end if

    if (present(time)) then
      scheme%time = time
    else
      scheme%time = out_units%time
    end if

    if (present(mass)) then
      scheme%mass = mass
    else
      scheme%mass = out_units%mass
    end if

    if (present(charge)) then
      scheme%charge = charge
    else
      scheme%charge = out_units%charge
    end if

    if (present(energy)) then
      scheme%energy = energy
    else
      scheme%energy = out_units%energy
    end if

    if (present(temp)) then
      scheme%temp = temp
    else
      scheme%temp = out_units%temp
    end if

    if (present(current)) then
      scheme%current = current
    else
      scheme%current = out_units%current
    end if

    if (present(luminosity)) then
      scheme%luminosity = luminosity
    else
      scheme%luminosity = out_units%luminosity
    end if

    if (present(angle)) then
      scheme%angle = angle
    else
      scheme%angle = out_units%angle
    end if

    if (present(pressure)) then
      scheme%pressure = pressure
    else
      scheme%pressure = out_units%pressure
    end if

    if (present(force)) then
      scheme%force = force
    else
      scheme%force = out_units%force
    end if

    if (present(velocity)) then
      scheme%velocity = velocity
    else
      scheme%velocity = out_units%velocity
    end if

    if (present(power)) then
      scheme%power = power
    else
      scheme%power = out_units%power
    end if

    if (present(surf_ten)) then
      scheme%surf_ten = surf_ten
    else
      scheme%surf_ten = out_units%surf_ten
    end if

    if (present(emf)) then
      scheme%emf = emf
    else
      scheme%emf = out_units%emf
    end if

    call scheme%recalc()

  end Function init_scheme

  Subroutine recalculate_scheme_units(scheme)
    !!-----------------------------------------------------------------------
    !!
    !! Test validity of units scheme
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2022
    !!-----------------------------------------------------------------------
    Class(units_scheme), Intent( InOut ) :: scheme

    scheme%length_unit = parse_unit_string(scheme%length)
    scheme%time_unit = parse_unit_string(scheme%time)
    scheme%mass_unit = parse_unit_string(scheme%mass)
    scheme%charge_unit = parse_unit_string(scheme%charge)
    scheme%energy_unit = parse_unit_string(scheme%energy)
    scheme%temp_unit = parse_unit_string(scheme%temp)
    scheme%current_unit = parse_unit_string(scheme%current)
    scheme%luminosity_unit = parse_unit_string(scheme%luminosity)
    scheme%angle_unit = parse_unit_string(scheme%angle)
    scheme%pressure_unit = parse_unit_string(scheme%pressure)
    scheme%force_unit = parse_unit_string(scheme%force)
    scheme%velocity_unit = parse_unit_string(scheme%velocity)
    scheme%power_unit = parse_unit_string(scheme%power)
    scheme%surf_ten_unit = parse_unit_string(scheme%surf_ten)
    scheme%emf_unit = parse_unit_string(scheme%emf)

  end Subroutine recalculate_scheme_units

  Subroutine test_scheme(scheme)
    !!-----------------------------------------------------------------------
    !!
    !! Test validity of units scheme
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2022
    !!-----------------------------------------------------------------------
    Class(units_scheme), Intent( In    ) :: scheme
    Real(kind=wp)          :: test
    Logical                :: stat

    ! Check units validity
    test = convert_units(1.0_wp, scheme%length, internal_units%length, stat)
    If (.not. stat) Call error_units(scheme%length, internal_units%length, 'Invalid output unit for length')

    test = convert_units(1.0_wp, scheme%time, internal_units%time, stat)
    If (.not. stat) Call error_units(scheme%time, internal_units%time, 'Invalid output unit for time')

    test = convert_units(1.0_wp, scheme%mass, internal_units%mass, stat)
    If (.not. stat) Call error_units(scheme%mass, internal_units%mass, 'Invalid output unit for mass')

    test = convert_units(1.0_wp, scheme%charge, internal_units%charge, stat)
    If (.not. stat) Call error_units(scheme%charge, internal_units%charge, 'Invalid output unit for charge')

    test = convert_units(1.0_wp, scheme%energy, internal_units%energy, stat)
    If (.not. stat) Call error_units(scheme%energy, internal_units%energy, 'Invalid output unit for energy')

    test = convert_units(1.0_wp, scheme%pressure, internal_units%pressure, stat)
    If (.not. stat) Call error_units(scheme%pressure, internal_units%pressure, 'Invalid output unit for pressure')

    test = convert_units(1.0_wp, scheme%force, internal_units%force, stat)
    If (.not. stat) Call error_units(scheme%force, internal_units%force, 'Invalid output unit for force')

    test = convert_units(1.0_wp, scheme%velocity, internal_units%velocity, stat)
    If (.not. stat) Call error_units(scheme%velocity, internal_units%velocity, 'Invalid output unit for velocity')

    test = convert_units(1.0_wp, scheme%power, internal_units%power, stat)
    If (.not. stat) Call error_units(scheme%power, internal_units%power, 'Invalid output unit for power')

    test = convert_units(1.0_wp, scheme%surf_ten, internal_units%surf_ten, stat)
    If (.not. stat) Call error_units(scheme%surf_ten, internal_units%surf_ten, 'Invalid output unit for surface tension')

    test = convert_units(1.0_wp, scheme%emf, internal_units%emf, stat)
    If (.not. stat) Call error_units(scheme%emf, internal_units%emf, 'Invalid output unit for emf')

  end Subroutine test_scheme

  Subroutine set_out_units(scheme)
    !!-----------------------------------------------------------------------
    !!
    !! Set an output units scheme
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type(units_scheme), Value :: scheme

    Call scheme%recalc()
    Call scheme%test()

    out_units = scheme

  end Subroutine set_out_units

  Subroutine to_out_units(val, from, res, out_unit, out_unit_full, stat)
    !!-----------------------------------------------------------------------
    !!
    !! Convert var to default output units
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2022
    !!-----------------------------------------------------------------------
    Real(kind=wp), Intent(in) :: val
    Character(Len=*), Intent(In) :: from
    Logical, Intent(Out), Optional :: stat

    Real(kind=wp), Intent(Out) :: res
    Character(Len=STR_LEN), Intent(Out), Optional :: out_unit
    Character(Len=STR_LEN), Intent(Out), Optional :: out_unit_full

    Type( unit_data ) :: from_unit, to_unit

    Type(unit_data), Dimension(ndims) :: out_unit_vals
    Integer :: i, curr_dim

    if (present(stat)) then
      stat = .true.
    end if

    from_unit = parse_unit_string(from)

    if (all(from_unit%dims == dims_energy)) then
      to_unit = out_units%energy_unit

    else if (all(from_unit%dims == dims_charge)) then
      to_unit = out_units%charge_unit

    else if (all(from_unit%dims == dims_pressure)) then
      to_unit = out_units%pressure_unit

    else if (all(from_unit%dims == dims_force)) then
      to_unit = out_units%force_unit

    else if (all(from_unit%dims == dims_velocity)) then
      to_unit = out_units%velocity_unit

    else if (all(from_unit%dims == dims_power)) then
      to_unit = out_units%power_unit

    else if (all(from_unit%dims == dims_surf_ten)) then
      to_unit = out_units%surf_ten_unit

    else if (all(from_unit%dims == dims_emf)) then
      to_unit = out_units%emf_unit

    else

      out_unit_vals = [ &
           out_units%mass_unit, &
           out_units%length_unit, &
           out_units%time_unit, &
           out_units%temp_unit, &
           out_units%current_unit, &
           out_units%luminosity_unit, &
           out_units%angle_unit]

      to_unit = null_unit

      do i = 1, ndims
        to_unit = to_unit * (out_unit_vals(i)**from_unit%dims(i))
      end do
    end if

    res = convert_units(val, from_unit, to_unit, stat)

    if (present(out_unit)) out_unit = to_unit%abbrev
    if (present(out_unit_full)) out_unit_full = to_unit%name

  end Subroutine to_out_units

  Function init_unit(name, abbrev, to_internal, mass, length, time, temp, mol, current, luminosity, angle)
    Type(unit_data) :: init_unit
    Character(Len = *), Intent( In    ) :: name, abbrev
    Real(Kind = wp), Intent( In    ) :: to_internal
    Integer, Optional, Intent( In    ) :: mass, length, time, temp, mol, current, luminosity, angle

    init_unit = unit_data(name=name, abbrev=abbrev, conversion_to_internal=to_internal)
    if (present(mass)) then
      init_unit%dims(1) = mass
    end if

    if (present(length)) then
      init_unit%dims(2) = length
    end if

    if (present(time)) then
      init_unit%dims(3) = time
    end if

    if (present(temp)) then
      init_unit%dims(4) = temp
    end if

    if (present(current)) then
      init_unit%dims(5) = current
    end if

    if (present(luminosity)) then
      init_unit%dims(6) = luminosity
    end if

    if (present(angle)) then
      init_unit%dims(7) = angle
    end if

  end Function init_unit

  Function unit_mult(a, b)
    Type(unit_data) :: unit_mult
    Class(unit_data), Intent( in ) :: a, b

    unit_mult = unit_data( &
         name = trim(a%name)//" "//b%name, &
         abbrev = trim(a%abbrev)//b%abbrev, &
         conversion_to_internal = a%conversion_to_internal*b%conversion_to_internal, &
         dims = a%dims + b%dims &
         )
  end Function unit_mult

  Function unit_div(a, b)
    Type(unit_data) :: unit_div
    Class(unit_data), Intent( in ) :: a, b

    unit_div = unit_data( &
         name = trim(a%name)//" per "//b%name, &
         abbrev = trim(a%abbrev)//"/"//b%abbrev, &
         conversion_to_internal = a%conversion_to_internal/b%conversion_to_internal, &
         dims = a%dims - b%dims &
         )
  end Function unit_div

  Function unit_pow(a, b)
    Type(unit_data) :: unit_pow
    Class(unit_data), Intent( in ) :: a
    Integer, Intent( in ) :: b
    Character(Len=8) :: b_str

    write(b_str, "(i0.1)") b

    select case (b)
    case (0)
      unit_pow = null_unit
    case (1)
      unit_pow = a
    case default
      unit_pow = unit_data( &
           name = trim(a%name)//"^"//trim(b_str), &
           abbrev = trim(a%abbrev)//"^"//trim(b_str), &
           conversion_to_internal = a%conversion_to_internal**b, &
           dims = a%dims*b &
           )
    end select

  end Function unit_pow

  Subroutine handle_decimal_prefix(string, factor)
    !!-----------------------------------------------------------------------
    !!
    !! Handle SI Decimal prefixes in units
    !! (NB removes the prefix from string)
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Character(len=*), Intent(inout) :: string
    Integer :: i
    Type(unit_data), intent(out) :: factor
    Character(Len=STR_LEN) :: tmp
    Character(len=*), Parameter :: prefix_symbol = "YZEPTGMk dcmunpfazy"
    Character(Len=*), Parameter :: number = "1234567890-+"
    Type(unit_data), Dimension(19), Parameter :: prefix = [ &
         & unit_data(name="Yotta", abbrev="Y", conversion_to_internal=1e24_wp), &
         & unit_data(name="Zetta", abbrev="Z", conversion_to_internal=1e21_wp), &
         & unit_data(name="Exa",   abbrev="E", conversion_to_internal=1e18_wp), &
         & unit_data(name="Peta",  abbrev="P", conversion_to_internal=1e15_wp), &
         & unit_data(name="Tera",  abbrev="T", conversion_to_internal=1e12_wp), &
         & unit_data(name="Giga",  abbrev="G", conversion_to_internal=1e9_wp), &
         & unit_data(name="Mega",  abbrev="M", conversion_to_internal=1e6_wp), &
         & unit_data(name="Kilo",  abbrev="k", conversion_to_internal=1e3_wp), &
         & null_unit, &
         & unit_data(name="Deci",  abbrev="d", conversion_to_internal=1e-1_wp), &
         & unit_data(name="Centi", abbrev="c", conversion_to_internal=1e-2_wp), &
         & unit_data(name="Milli", abbrev="m", conversion_to_internal=1e-3_wp), &
         & unit_data(name="Micro", abbrev="u", conversion_to_internal=1e-6_wp), &
         & unit_data(name="Nano",  abbrev="n", conversion_to_internal=1e-9_wp), &
         & unit_data(name="Pico",  abbrev="p", conversion_to_internal=1e-12_wp), &
         & unit_data(name="Femto", abbrev="f", conversion_to_internal=1e-15_wp), &
         & unit_data(name="Atto",  abbrev="a", conversion_to_internal=1e-18_wp), &
         & unit_data(name="Zepto", abbrev="z", conversion_to_internal=1e-21_wp), &
         & unit_data(name="Yocto", abbrev="y", conversion_to_internal=1e-24_wp)]

    factor = null_unit
    tmp = string(2:)
    call lower_case(tmp)
    if (.not. units_table%in(tmp) .and. verify(trim(string(2:)), number) /= 0) then
      i = index(prefix_symbol, string(2:2))
      tmp = string(3:)
      call lower_case(tmp)
      if (i < 1 .or. .not. units_table%in(tmp)) call error(0, "Unit not found "//string(2:))
      factor = prefix(i)
      string = string(1:1) // string(3:)
    end if

  end Subroutine handle_decimal_prefix

  Function parse_unit_string(string) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Convert a unit string into its to_internal representation
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( unit_data ) :: output
    Type( unit_data ) :: tmp_unit
    Character(Len=*), Intent ( in ) :: string
    Character(Len=STR_LEN) :: curr_parse
    Character(Len=STR_LEN), Dimension(:), Allocatable :: parsed
    Character(Len=*), Parameter :: number = "1234567890-+"
    Type(unit_data) :: factor
    Integer :: i
    Integer :: n

    call decompose_unit_string("."//string, parsed)

    output = null_unit

    ! Handle powers first
    do i = size(parsed), 1, -1
      curr_parse = parsed(i)
      if (curr_parse(1:1) == "^") then
        if (verify(trim(curr_parse(2:)), number) /= 0) call error(0, "Non numeric power issued")
        read(curr_parse(2:), "(i8.1)") n
        curr_parse = parsed(i-1)
        call handle_decimal_prefix(curr_parse, factor)
        call lower_case(curr_parse)
        call units_table%get(curr_parse(2:), tmp_unit)
        tmp_unit = (factor*tmp_unit) ** n
        select case (curr_parse(1:1))
        case (".")
          output = output * tmp_unit
        case ("/")
          output = output / tmp_unit
        case default
          call error(0, "Cannot parse unit string"//string)
        end select
        parsed(i) = "."
        parsed(i-1) = "."
      end if
    end do

    do i = 1, size(parsed)
      curr_parse = parsed(i)

      call handle_decimal_prefix(curr_parse, factor)
      call lower_case(curr_parse)
      if (verify(trim(curr_parse(2:)), number) /= 0) then
        call units_table%get(curr_parse(2:), tmp_unit)
      else if (trim(curr_parse(2:)) /= "") then
        tmp_unit = null_unit
        read(curr_parse(2:), *) tmp_unit%conversion_to_internal
      else
        cycle
      end if

      select case (curr_parse(1:1))
      case (".")
        output = output * (factor * tmp_unit)
      case ("/")
        output = output / (factor * tmp_unit)
      case default
        call error(0, "Cannot parse unit string "//string)
      end select
    end do

    output%abbrev = string

  end Function parse_unit_string

  Subroutine decompose_unit_string(string, output)
    !!-----------------------------------------------------------------------
    !!
    !! Split unit string into components separated by product, quotient or exponent
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Character(Len=STR_LEN), Dimension(:), Allocatable :: output
    Character(Len=*), Intent( In ) :: string
    Character(Len=*), Parameter :: alphabet = "-+1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_%'"//'"'
    Character(Len=*), Parameter :: punc = "./^"
    Integer :: nParts
    Integer :: i, j, k
    Integer :: ierr

    ierr = 0
    if (allocated(output)) deallocate(output, stat=ierr)
    if (ierr /= 0) call error_dealloc("output", "decompose_unit_string")

    nParts = count([(verify(string(j:j), punc) == 0, j=1,len(string))])
    allocate(output(nParts), stat=ierr)
    if (ierr /= 0) call error_alloc("output", "decompose_unit_string")

    j = 1
    do i = 1, nParts
      k = verify(string(j+1:), alphabet)
      if (k == 0) then
        k = len(string) + 1
      else
        k = j+k
      end if
      output(i) = string(j:k-1)
      j = k
    end do


  end Subroutine decompose_unit_string

  Subroutine set_timestep(val)
    !!-----------------------------------------------------------------------
    !!
    !! Set timestep unit (must be done post-initialisation)
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Real(kind=wp), Intent( In    ) :: val

    call units_table%set("steps", init_unit(abbrev="steps", name="Timestep", time=1, to_internal=1.0_wp/val))

  end Subroutine set_timestep

  Subroutine get_unit( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get unit type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( units_hash_table ), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    type(unit_data)            , Intent(   Out ) :: val
    type(unit_data), Intent( In    ), Optional :: default
    Class( * ), Pointer :: stuff

    call table%get_cont(key, default, stuff)

    Select Type( stuff )
    Type is ( unit_data )
      val = stuff
    Class Default
      Call error(0, 'Trying to get unit from a not unit')
    End Select
    deallocate(stuff)
    nullify(stuff)

  End Subroutine get_unit

end Module units
