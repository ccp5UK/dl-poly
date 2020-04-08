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
  Use hashables, only : unit_data
  Use hash, only: hash_table
  Use errors_warnings, only : error, error_alloc, error_dealloc
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
    Type(unit_data) :: null_unit
    Real(kind=wp), Parameter :: planck_internal = 6.350780668_wp
    Real(kind=wp), Parameter :: electron_charge = 1.6021766208e-19_wp ! C
    Real(kind=wp), Parameter :: avogadro = 6.022140857e23_wp

    Real(kind=wp), Parameter :: metre = 1e10_wp
    Real(kind=wp), Parameter :: bohr = 0.52916_wp

    Real(kind=wp), Parameter :: joulepmol = 0.1_wp
    Real(kind=wp), Parameter :: caloriepmol = 0.41842_wp
    Real(kind=wp), Parameter :: ev = 0.00010364272224984791_wp
    Real(kind=wp), Parameter :: hartree = ev*27.211396641308_wp

    Real(kind=wp), Parameter :: kilogram = 6.022e26

    Real(kind=wp), Parameter :: second = 1e12_wp

    Real(kind=wp), Parameter :: atmosphere = 0.006101929957459297_wp
    Real(kind=wp), Parameter :: pascal = 6.02213665167516e-8_wp

    call units_table%init(100)

    ! Time

    call units_table%set("internal_t", null_unit%init(abbrev="internal_t", name="Picosecond", time=1, to_internal=1.0_wp))
    call units_table%set("hr", null_unit%init(abbrev="hr", name="Hour", time=1, to_internal=3600.0_wp*second))
    call units_table%set("min", null_unit%init(abbrev="min", name="Minute", time=1, to_internal=60.0_wp*second))
    call units_table%set("s", null_unit%init(abbrev="s", name="Second", time=1, to_internal=second))
    call units_table%set("ms", null_unit%init(abbrev="ms", name="Millisecond", time=1, to_internal=1e-3_wp*second))
    call units_table%set("us", null_unit%init(abbrev="us", name="Microsecond", time=1, to_internal=1e-6_wp*second))
    call units_table%set("ns", null_unit%init(abbrev="ns", name="Nanosecond", time=1, to_internal=1e-9_wp*second))
    call units_table%set("ps", null_unit%init(abbrev="ps", name="Picosecond", time=1, to_internal=1e-12_wp*second))
    call units_table%set("fs", null_unit%init(abbrev="fs", name="Femtosecond", time=1, to_internal=1e-15_wp*second))

    ! Length

    call units_table%set("internal_l", null_unit%init(abbrev="internal_l", name="Angstrom", length=1, to_internal=1.0_wp))
    call units_table%set("ang", null_unit%init(abbrev="ang", name="Angstrom", length=1, to_internal=1e-10_wp*metre))
    call units_table%set("bohr", null_unit%init(abbrev="bohr", name="Bohr", length=1, to_internal=bohr))
    call units_table%set("m", null_unit%init(abbrev="m", name="Metre", length=1, to_internal=metre))
    call units_table%set("cm", null_unit%init(abbrev="cm", name="Centimetre", length=1, to_internal=1e-2_wp*metre))
    call units_table%set("nm", null_unit%init(abbrev="nm", name="Nanometre", length=1, to_internal=1e-9_wp*metre))
    call units_table%set("pm", null_unit%init(abbrev="pm", name="Picometre", length=1, to_internal=1e-12_wp*metre))
    call units_table%set("fm", null_unit%init(abbrev="fm", name="Femtometre", length=1, to_internal=1e-15_wp*metre))
    call units_table%set('"', null_unit%init(abbrev='"', name="Inch", length=1, to_internal=2.54e-2*metre))

    ! Mass

    call units_table%set("internal_m", null_unit%init(abbrev="internal_m", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("da", null_unit%init(abbrev="Da", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("amu", null_unit%init(abbrev="amu", name="Atomic Mass Unit", mass=1, to_internal=1.0_wp))
    call units_table%set("kg", null_unit%init(abbrev="kg", name="Kilogram", mass=1, to_internal=kilogram))
    call units_table%set("g", null_unit%init(abbrev="g", name="Gram", mass=1, to_internal=1e-3*kilogram))
    call units_table%set("mg", null_unit%init(abbrev="mg", name="Milligram", mass=1, to_internal=1e-6*kilogram))

    ! Charge

    call units_table%set("internal_q", &
         & null_unit%init(abbrev="internal_q", name="Elementary charge", current=1, time=1, to_internal=1.0_wp))
    call units_table%set("q_e", &
         & null_unit%init(abbrev="q_e", name="Elementary charge", current=1, time=1, to_internal=electron_charge))
    call units_table%set("e", &
         & null_unit%init(abbrev="e", name="Elementary charge", current=1, time=1, to_internal=electron_charge))
    call units_table%set("c", &
         & null_unit%init(abbrev="C", name="Coulomb", current=1, time=1, to_internal=1.0_wp / electron_charge))

    ! Energy

    call units_table%set("internal_e", &
         & null_unit%init(abbrev="internal_e", name="10 J/mol", mass=1, length=2, time=-2, to_internal=1.0_wp))
    call units_table%set("j/mol", &
         & null_unit%init(abbrev="J/mol", name="Joule per mole", mol=-1, mass=1, length=2, time=-2, to_internal=joulepmol))
    call units_table%set("kj/mol", &
         & null_unit%init(abbrev="kJ/mol", name="Kilojoule per mole", &
         & mol=-1, mass=1, length=2, time=-2, to_internal=1e3_wp*joulepmol))
    call units_table%set("j", null_unit%init(abbrev="J", name="Joule", mass=1, length=2, time=-2, to_internal=joulepmol*avogadro))
    call units_table%set("kj", &
         & null_unit%init(abbrev="kJ", name="Kilojoule", mass=1, length=2, time=-2, to_internal=1e3_wp*joulepmol*avogadro))
    call units_table%set("cal", &
         & null_unit%init(abbrev="Cal", name="Calorie", mass=1, length=2, time=-2, to_internal=caloriepmol*avogadro))
    call units_table%set("ev", null_unit%init(abbrev="ev", name="Electron-volt", mass=1, length=2, time=-2, to_internal=ev))
    call units_table%set("mev", &
         & null_unit%init(abbrev="meV", name="Millielectron-volt", mass=1, length=2, time=-2, to_internal=1e-3*ev))
    call units_table%set("ha", null_unit%init(abbrev="Ha", name="Hartree", mass=1, length=2, time=-2, to_internal=hartree))
    call units_table%set("mha", &
         & null_unit%init(abbrev="mHa", name="Millihartree", mass=1, length=2, time=-2, to_internal=1e-3*hartree))
    call units_table%set("ry", null_unit%init(abbrev="Ry", name="Rydberg", mass=1, length=2, time=-2, to_internal=0.5_wp*hartree))
    call units_table%set("mry", &
         & null_unit%init(abbrev="mRy", name="Millirydberg", mass=1, length=2, time=-2, to_internal=0.5e-3_wp*hartree))
    call units_table%set("k", null_unit%init(abbrev="K", name="Kelvin", mass=1, length=2, time=-2, to_internal=boltz))

    ! Pressure

    call units_table%set("internal_p", &
         & null_unit%init(abbrev="internal_p", name="163 atm", mass=1, length=-1, time=-2, to_internal=1.0_wp))
    call units_table%set("atm", null_unit%init(abbrev="atm", name="Atmosphere", mass=1, length=-1, time=-2, to_internal=atmosphere))
    call units_table%set("katm", &
         & null_unit%init(abbrev="katm", name="Kiloatmosphere", mass=1, length=-1, time=-2, to_internal=1e3_wp*atmosphere))
    call units_table%set("matm", &
         & null_unit%init(abbrev="Matm", name="Megaatmosphere", mass=1, length=-1, time=-2, to_internal=1e6_wp*atmosphere))
    call units_table%set("gatm", &
         & null_unit%init(abbrev="Gatm", name="Gigaatmosphere", mass=1, length=-1, time=-2, to_internal=1e9_wp*atmosphere))
    call units_table%set("pa", null_unit%init(abbrev="Pa", name="Pascal", mass=1, length=-1, time=-2, to_internal=pascal))
    call units_table%set("kpa", &
         & null_unit%init(abbrev="kPa", name="Kilopascal", mass=1, length=-1, time=-2, to_internal=1e3_wp*pascal))
    call units_table%set("mpa", &
         & null_unit%init(abbrev="MPa", name="Megapascal", mass=1, length=-1, time=-2, to_internal=1e6_wp*pascal))
    call units_table%set("gpa", &
         & null_unit%init(abbrev="GPa", name="Gigapascal", mass=1, length=-1, time=-2, to_internal=1e9_wp*pascal))
    call units_table%set("tpa", &
         & null_unit%init(abbrev="TPa", name="Terapascal", mass=1, length=-1, time=-2, to_internal=1e12_wp*pascal))

    ! Mols

    call units_table%set("mol", null_unit%init(abbrev="mol", name="mole", mol=1, to_internal=avogadro))

    ! Unitless

    call units_table%set("", null_unit%init(abbrev="", name="", to_internal=1.0_wp))

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

  Function parse_unit_string(string) result(output)
    Type( unit_data ) :: output
    Type( unit_data ) :: tmp_unit
    Character(Len=*), Intent ( in ) :: string
    Character(Len=256), Dimension(:), Allocatable :: parsed
    Character(Len=*), Parameter :: number = "1234567890"
    Integer :: i
    Integer :: n

    call decompose_unit_string("."//string, parsed)

    ! Handle powers first
    do i = 1, size(parsed)
       if (parsed(i)(1:1) == "^") then
          if (verify(parsed(i)(2:), number) /= 0) call error(0, "Non numeric power issued")
          read(parsed(i)(2:), "(i8.1)") n
          call units_table%get(parsed(i-1)(2:), tmp_unit)
          tmp_unit = tmp_unit ** n
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
       select case (parsed(i)(1:1))
       case (".")
          call units_table%get(parsed(i)(2:), tmp_unit)
          print*, i, tmp_unit
          output = output * tmp_unit
       case ("/")
          call units_table%get(parsed(i)(2:), tmp_unit)
          output = output / tmp_unit
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
