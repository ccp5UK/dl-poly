Module hashables

  Use kinds, only : wp
  Use errors_warnings, only : error
  Implicit None

  Private
  Integer, Parameter :: STR_LEN = 256

  Type, Public :: unit_data
     !! Type containing data corresponding to units
     Integer, Dimension(7) :: dims = [0, 0, 0, 0, 0, 0, 0] ! mass length time temp mol current luminosity
     Character(Len=STR_LEN) :: name
     Character(Len=STR_LEN) :: abbrev
     Real(Kind=wp) :: conversion_to_internal
   contains
     Generic, Public :: Operator(*) => unit_mult
     Generic, Public :: Operator(/) => unit_div
     Generic, Public :: Operator(**) => unit_pow
     Procedure, Pass, Private :: unit_mult, unit_div, unit_pow
     ! Procedure, Pass, Public :: to_internal, from_internal
     Procedure, Nopass :: init => init_unit
  end type unit_data

  Type, Public :: control_parameter
     !! Type containing breakdown of control parameter
     Character(Len=STR_LEN) :: key
     Character(Len=STR_LEN) :: val
     Character(Len=STR_LEN) :: unit_type
     Character(Len=STR_LEN) :: description
     Character(Len=STR_LEN) :: default
  End Type control_parameter

  Interface resolve
     Module Procedure resolve_int
     Module Procedure resolve_double
     Module Procedure resolve_param
     Module Procedure resolve_unit
  End Interface resolve

  Public :: resolve

contains

  type(unit_data) Function unit_mult(a, b)
    Class(unit_data), Intent( in ) :: a, b

    unit_mult = unit_data( &
         name = a%name//" "//b%name, &
         abbrev = a%abbrev//b%abbrev, &
         conversion_to_internal = a%conversion_to_internal*b%conversion_to_internal, &
         dims = a%dims + b%dims &
         )
  end Function unit_mult

  type(unit_data) Function unit_div(a, b)
    Class(unit_data), Intent( in ) :: a, b

    unit_div = unit_data( &
         name = a%name//" per "//b%name, &
         abbrev = a%abbrev//"/"//b%abbrev, &
         conversion_to_internal = a%conversion_to_internal/b%conversion_to_internal, &
         dims = a%dims - b%dims &
         )
  end Function unit_div

  type(unit_data) Function unit_pow(a, b)
    Class(unit_data), Intent( in ) :: a
    Integer, Intent( in ) :: b
    Character(Len=8) :: b_str

    write(b_str, "(i0.1)") b

    unit_pow = unit_data( &
         name = a%name//"^"//trim(b_str), &
         abbrev = a%abbrev//"^"//trim(b_str), &
         conversion_to_internal = a%conversion_to_internal**b, &
         dims = a%dims*b &
    )

  end Function unit_pow

  type(unit_data) Function init_unit(name, abbrev, to_internal, mass, length, time, temp, mol, current, luminosity)
    Character(Len = *), Intent( In    ) :: name, abbrev
    Real(Kind = wp), Intent( In    ) :: to_internal
    Integer, Optional, Intent( In    ) :: mass, length, time, temp, mol, current, luminosity

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

    if (present(mol)) then
       init_unit%dims(5) = mol
    end if

    if (present(current)) then
       init_unit%dims(6) = current
    end if

    if (present(luminosity)) then
       init_unit%dims(7) = luminosity
    end if

  end Function init_unit


  Subroutine resolve_int( a, b )

    Implicit None

    Integer   , Intent(   Out ) :: a
    Class( * ), Intent( In    ) :: b

    Select Type( b )
    Type is ( Integer )
       a = b
    Class Default
       Call error(0, 'Trying to get integer from a not integer')
    End Select

  End Subroutine resolve_int

  Subroutine resolve_double( a, b )

    Implicit None

    Real( wp ), Intent(   Out ) :: a
    Class( * ), Intent( In    ) :: b

    Select Type( b )
    Type is ( Real( wp ) )
       a = b
    Class Default
       Call error(0, 'Trying to get double from a not double')
    End Select

  End Subroutine resolve_double

  Subroutine resolve_param( a, b )

    Implicit None

    Type(control_parameter), Intent(   Out ) :: a
    Class( * ), Intent( In    ) :: b

    Select Type( b )
    Type is ( control_parameter )
       a = b
    Class Default
       Call error(0, 'Trying to get param from a not param')
    End Select

  End Subroutine resolve_param

  Subroutine resolve_unit( a, b )

    Implicit None

    Type(unit_data), Intent(   Out ) :: a
    Class( * ), Intent( In    ) :: b

    Select Type( b )
    Type is ( unit_data )
       a = b
    Class Default
       Call error(0, 'Trying to get unit from a not unit')
    End Select

  End Subroutine resolve_unit

end Module hashables
