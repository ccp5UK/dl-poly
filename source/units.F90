Module units
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle unit conversion in dlpoly
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins march 2020
  !!-----------------------------------------------------------------------

  Use kinds, only : wp
  Use constants, only : boltz
  Use hashables, only : unit_data, null_unit, init_unit
  Use hash, only: hash_table
  Use errors_warnings, only : error, error_alloc, error_dealloc
  Use parse, only : lower_case
  Implicit None

  Private

  type(hash_table), Private, save :: units_table

  Public :: initialise_units
  Public :: convert_units

contains

  Subroutine initialise_units()
    !!-----------------------------------------------------------------------
    !!
    !! Module to handle unit conversion in dlpoly
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Real(kind=wp), Parameter :: planck_internal = 6.350780668_wp
    Real(kind=wp), Parameter :: electron_charge = 1.6021766208e-19_wp ! C
    Real(kind=wp), Parameter :: avogadro = 6.022140857e23_wp

    Real(kind=wp), Parameter :: metre = 1e10_wp
    Real(kind=wp), Parameter :: bohr = 0.52916_wp

    Real(kind=wp), Parameter :: joulepmol = 0.1_wp
    Real(kind=wp), Parameter :: caloriepmol = 0.41842_wp
    Real(kind=wp), Parameter :: ev = 0.00010364272224984791_wp
    Real(kind=wp), Parameter :: hartree = ev*27.211396641308_wp

    Real(kind=wp), Parameter :: kilogram = 6.022e26_wp

    Real(kind=wp), Parameter :: second = 1e12_wp

    Real(kind=wp), Parameter :: atmosphere = 0.006101929957459297_wp
    Real(kind=wp), Parameter :: pascal = 6.02213665167516e-8_wp

    Real(kind=wp), Parameter :: newton = kilogram*metre/second**2

    Real(kind=wp), Parameter :: gravity = 9.81e-14_wp

    call units_table%init(100)

    ! Time

    call units_table%set("internal_t", init_unit(abbrev="internal_t", name="Picosecond", time=1, to_internal=1.0_wp))
    call units_table%set("hr", init_unit(abbrev="hr", name="Hour", time=1, to_internal=3600.0_wp*second))
    call units_table%set("min", init_unit(abbrev="min", name="Minute", time=1, to_internal=60.0_wp*second))
    call units_table%set("s", init_unit(abbrev="s", name="Second", time=1, to_internal=second))
    ! call units_table%set("ms", init_unit(abbrev="ms", name="Millisecond", time=1, to_internal=1e-3_wp*second))
    ! call units_table%set("us", init_unit(abbrev="us", name="Microsecond", time=1, to_internal=1e-6_wp*second))
    ! call units_table%set("ns", init_unit(abbrev="ns", name="Nanosecond", time=1, to_internal=1e-9_wp*second))
    ! call units_table%set("ps", init_unit(abbrev="ps", name="Picosecond", time=1, to_internal=1e-12_wp*second))
    ! call units_table%set("fs", init_unit(abbrev="fs", name="Femtosecond", time=1, to_internal=1e-15_wp*second))

    ! Length

    call units_table%set("internal_l", init_unit(abbrev="internal_l", name="Angstrom", length=1, to_internal=1.0_wp))
    call units_table%set("ang", init_unit(abbrev="ang", name="Angstrom", length=1, to_internal=1e-10_wp*metre))
    call units_table%set("bohr", init_unit(abbrev="bohr", name="Bohr", length=1, to_internal=bohr))
    call units_table%set("m", init_unit(abbrev="m", name="Metre", length=1, to_internal=metre))
    ! call units_table%set("cm", init_unit(abbrev="cm", name="Centimetre", length=1, to_internal=1e-2_wp*metre))
    ! call units_table%set("nm", init_unit(abbrev="nm", name="Nanometre", length=1, to_internal=1e-9_wp*metre))
    ! call units_table%set("pm", init_unit(abbrev="pm", name="Picometre", length=1, to_internal=1e-12_wp*metre))
    ! call units_table%set("fm", init_unit(abbrev="fm", name="Femtometre", length=1, to_internal=1e-15_wp*metre))
    call units_table%set('"', init_unit(abbrev='"', name="Inch", length=1, to_internal=2.54e-2*metre))

    ! Mass

    call units_table%set("internal_m", init_unit(abbrev="internal_m", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("da", init_unit(abbrev="Da", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("amu", init_unit(abbrev="amu", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    ! call units_table%set("kg", init_unit(abbrev="kg", name="Kilogram", mass=1, to_internal=kilogram))
    call units_table%set("g", init_unit(abbrev="g", name="Gram", mass=1, to_internal=1e-3*kilogram))
    ! call units_table%set("mg", init_unit(abbrev="mg", name="Milligram", mass=1, to_internal=1e-6*kilogram))

    ! Charge

    call units_table%set("internal_q", &
         & init_unit(abbrev="internal_q", name="Elementary charge", current=1, time=1, to_internal=1.0_wp))
    call units_table%set("q_e", &
         & init_unit(abbrev="q_e", name="Elementary charge", current=1, time=1, to_internal=electron_charge))
    call units_table%set("e", &
         & init_unit(abbrev="e", name="Elementary charge", current=1, time=1, to_internal=electron_charge))
    call units_table%set("c", &
         & init_unit(abbrev="C", name="Coulomb", current=1, time=1, to_internal=1.0_wp / electron_charge))

    ! Energy

    call units_table%set("internal_e", &
         & init_unit(abbrev="internal_e", name="10 J/mol", mass=1, length=2, time=-2, to_internal=1.0_wp))
    ! call units_table%set("j/mol", &
    !      & init_unit(abbrev="J/mol", name="Joule per mole", mol=-1, mass=1, length=2, time=-2, to_internal=joulepmol))
    ! call units_table%set("kj/mol", &
    !      & init_unit(abbrev="kJ/mol", name="Kilojoule per mole", &
    !      & mol=-1, mass=1, length=2, time=-2, to_internal=1e3_wp*joulepmol))
    call units_table%set("j", init_unit(abbrev="J", name="Joule", mass=1, length=2, time=-2, to_internal=joulepmol*avogadro))
    ! call units_table%set("kj", &
    !      & init_unit(abbrev="kJ", name="Kilojoule", mass=1, length=2, time=-2, to_internal=1e3_wp*joulepmol*avogadro))
    call units_table%set("cal", &
         & init_unit(abbrev="Cal", name="Calorie", mass=1, length=2, time=-2, to_internal=caloriepmol*avogadro))
    ! call units_table%set("ev", init_unit(abbrev="ev", name="Electron-volt", mass=1, length=2, time=-2, to_internal=ev))
    ! call units_table%set("mev", &
    !      & init_unit(abbrev="meV", name="Millielectron-volt", mass=1, length=2, time=-2, to_internal=1e-3*ev))
    call units_table%set("ha", init_unit(abbrev="Ha", name="Hartree", mass=1, length=2, time=-2, to_internal=hartree))
    ! call units_table%set("mha", &
    !      & init_unit(abbrev="mHa", name="Millihartree", mass=1, length=2, time=-2, to_internal=1e-3*hartree))
    call units_table%set("ry", init_unit(abbrev="Ry", name="Rydberg", mass=1, length=2, time=-2, to_internal=0.5_wp*hartree))
    ! call units_table%set("mry", &
    !      & init_unit(abbrev="mRy", name="Millirydberg", mass=1, length=2, time=-2, to_internal=0.5e-3_wp*hartree))
    call units_table%set("k", init_unit(abbrev="K", name="Kelvin", mass=1, length=2, time=-2, to_internal=boltz))

    ! Pressure

    call units_table%set("internal_p", &
         & init_unit(abbrev="internal_p", name="163 atm", mass=1, length=-1, time=-2, to_internal=1.0_wp))
    call units_table%set("atm", init_unit(abbrev="atm", name="Atmosphere", mass=1, length=-1, time=-2, to_internal=atmosphere))
    ! call units_table%set("katm", &
    !      & init_unit(abbrev="katm", name="Kiloatmosphere", mass=1, length=-1, time=-2, to_internal=1e3_wp*atmosphere))
    ! call units_table%set("matm", &
    !      & init_unit(abbrev="Matm", name="Megaatmosphere", mass=1, length=-1, time=-2, to_internal=1e6_wp*atmosphere))
    ! call units_table%set("gatm", &
    !      & init_unit(abbrev="Gatm", name="Gigaatmosphere", mass=1, length=-1, time=-2, to_internal=1e9_wp*atmosphere))
    call units_table%set("pa", init_unit(abbrev="Pa", name="Pascal", mass=1, length=-1, time=-2, to_internal=pascal))
    ! call units_table%set("kpa", &
    !      & init_unit(abbrev="kPa", name="Kilopascal", mass=1, length=-1, time=-2, to_internal=1e3_wp*pascal))
    ! call units_table%set("mpa", &
    !      & init_unit(abbrev="MPa", name="Megapascal", mass=1, length=-1, time=-2, to_internal=1e6_wp*pascal))
    ! call units_table%set("gpa", &
    !      & init_unit(abbrev="GPa", name="Gigapascal", mass=1, length=-1, time=-2, to_internal=1e9_wp*pascal))
    ! call units_table%set("tpa", &
    !      & init_unit(abbrev="TPa", name="Terapascal", mass=1, length=-1, time=-2, to_internal=1e12_wp*pascal))

    ! Force

    call units_table%set("internal_f", &
         & init_unit(abbrev="internal_f", name="internal_f", mass=1, length=1, time=-2, to_internal=1.0_wp))
    call units_table%set("n", &
         & init_unit(abbrev="N", name="Newton", mass=1, length=1, time=-2, to_internal=newton))
    ! call units_table%set("kn", &
    !      & init_unit(abbrev="kN", name="Kilonewton", mass=1, length=1, time=-2, to_internal=1e3_wp*newton)
    ! call units_table%set("mn", &
    !      & init_unit(abbrev="MN", name="Meganewton", mass=1, length=1, time=-2, to_internal=1e6_wp*newton)
    call units_table%set("dyn", &
         & init_unit(abbrev="dyn", name="Dyne", mass=1, length=1, time=-2, to_internal=1e-5_wp*newton))

    ! Acceleration

    call units_table%set("grav", &
         & init_unit(abbrev="grav", name="gravitation", length=1, time=-2, to_internal=gravity))

    ! Voltage

    call units_table%set("v", &
         & init_unit(abbrev="V", name="Volt", mass=1, length=2, time=-3, current=-1, to_internal=1.0_wp))

    ! Mols

    call units_table%set("mol", init_unit(abbrev="mol", name="mole", mol=1, to_internal=avogadro))

    ! Unitless

    call units_table%set("%", init_unit(abbrev="%", name="%", to_internal=0.01_wp))
    call units_table%set("", init_unit(abbrev="", name="", to_internal=1.0_wp))

  End Subroutine initialise_units

  Function convert_units(val, from, to) result(res)
    Real(kind=wp), Intent(in) :: val
    Real(kind=wp) :: res
    Character(Len=*) :: from, to
    Type( unit_data ) :: from_unit, to_unit
    Type( unit_data ) :: output

    output = output%init("", "", 1.0_wp)
    from_unit = parse_unit_string(from)
    to_unit = parse_unit_string(to)
    output = from_unit / to_unit
    res = val / output%conversion_to_internal

  end Function convert_units

  Subroutine handle_decimal_prefix(string, factor)
    Character(len=*), Intent(inout) :: string
    Integer :: i
    Type(unit_data), intent(out) :: factor
    Character(len=*), Parameter :: prefix_symbol = "YZEPTGMk munpfazy"
    Type(unit_data), Dimension(17), Parameter :: prefix = [(&
         & unit_data(name=prefix_symbol(i:i), abbrev=prefix_symbol(i:i), &
         & conversion_to_internal = 10.0_wp**(i*3)), i = 8, -8, -1)]

    factor = null_unit

    if (.not. units_table%in(string(2:))) then
       i = index(prefix_symbol, string(2:2))
       if (i < 1 .or. .not. units_table%in(string(3:))) call error(0, "Unit not found "//string(2:))
       factor = prefix(i)
       string = string(1:1) // string(3:)
    end if

  end Subroutine handle_decimal_prefix

  Function parse_unit_string(string) result(output)
    Type( unit_data ) :: output
    Type( unit_data ) :: tmp_unit
    Character(Len=*), Intent ( in ) :: string
    Character(Len=256) :: curr_parse
    Character(Len=256), Dimension(:), Allocatable :: parsed
    Character(Len=*), Parameter :: number = "1234567890"
    Type(unit_data) :: factor
    Integer :: i
    Integer :: n

    call decompose_unit_string("."//string, parsed)

    output = null_unit

    ! Handle powers first
    do i = 1, size(parsed)
       if (parsed(i)(1:1) == "^") then
          if (verify(parsed(i)(2:), number) /= 0) call error(0, "Non numeric power issued")
          read(parsed(i)(2:), "(i8.1)") n
          curr_parse = parsed(i-1)
          call handle_decimal_prefix(curr_parse, factor)
          call lower_case(curr_parse)
          call units_table%get(curr_parse(2:), tmp_unit)
          tmp_unit = (factor*tmp_unit) ** n
          select case (parsed(i-1)(1:1))
          case (".")
             output = output * tmp_unit
          case ("/")
             output = output / tmp_unit
          case default
             call error(0, "Cannot parse unit string"//string)
          end select
          parsed(i) = ""
          parsed(i-1) = ""
       end if
    end do

    do i = 1, size(parsed)
       curr_parse = parsed(i)
       call handle_decimal_prefix(curr_parse, factor)
       call lower_case(curr_parse)
       call units_table%get(curr_parse(2:), tmp_unit)

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
    Character(Len=256), Dimension(:), Allocatable :: output
    Character(Len=*), Intent( In ) :: string
    Character(Len=*), Parameter :: alphabet = "1234567890abcdefghijklmnopqrstuvwxyz_"
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




end Module units
