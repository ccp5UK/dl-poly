Module control_parameters
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle control parameters in dlpoly
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins april 2020
  !!-----------------------------------------------------------------------

  Use errors_warnings,               Only: error,&
                                           error_units
  Use hash,                          Only: MAX_KEY,&
                                           hash_table
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use kinds,                         Only: STR_LEN,&
                                           wp
  Use parse,                         Only: get_word,&
                                           word_2_real
  Use units,                         Only: convert_units

  Implicit None

  Private

  !> Data types enumeration
  Integer, Parameter, Public :: DATA_NULL = 0, DATA_INT = 1, DATA_FLOAT = 2, DATA_STRING = 3, &
                                DATA_BOOL = 4, DATA_OPTION = 5, DATA_VECTOR3 = 6, DATA_VECTOR6 = 7

  Type, Public, Extends(hash_table) :: parameters_hash_table
Contains
  !> Update get to include params
  Generic, Public  :: get => get_param
  Procedure, Private :: get_param
  !> Set retrieve up to parse stored params
  Generic, Public  :: retrieve => retrieve_option_or_string, retrieve_float, &
       & retrieve_vector_real, retrieve_vector_int, retrieve_int, retrieve_bool
  Procedure, Pass(table), Private :: retrieve_option_or_string, retrieve_float, &
       & retrieve_vector_real, retrieve_vector_int, retrieve_int, retrieve_bool
  !> Check if list of things is set
  Procedure, Public :: is_any_set
  Procedure :: is_set_single, is_all_set, num_set
  Generic, Public :: is_set => is_set_single, is_all_set
  Procedure, Private, Pass :: control_help_single, control_help_all
  Generic, Public :: help => control_help_single, control_help_all
End Type parameters_hash_table

Type, Public :: control_parameter
  !! Type containing breakdown of control parameter
  !> Internal key name
  Character(Len=MAX_KEY) :: key = ""
  !> Long name to be printed on high print_level
  Character(Len=STR_LEN) :: name = ""
  !> Current value -- Initialise to default
  Character(Len=STR_LEN) :: val = " "
  !> User specified units
  Character(Len=MAX_KEY) :: units = ""
  !> Units to be converted to internally
  Character(Len=MAX_KEY) :: internal_units = ""
  !> Information to be printed with help
  Character(Len=STR_LEN) :: description = ""
  !> Control parameter data type (int, float, vector3, vector6, string, bool, option)
  Integer :: data_type = DATA_NULL
  !> Is value set
  Logical :: set = .false.
Contains
  Procedure, Private :: write_control_param
  Generic :: Write (formatted) => write_control_param
End Type control_parameter

Public :: dump_parameters
Public :: control_help_single, control_help_all
Public :: print_set

Contains

Subroutine control_help_single(params, key)
    Class(parameters_hash_table), Intent(In   ) :: params
    Character(Len=*),             Intent(In   ) :: key

    Type(control_parameter) :: param

  Call params%get(key, param)
  Call write_control_param_help(param, error_unit)

End Subroutine control_help_single

Subroutine control_help_all(params)
    Class(parameters_hash_table), Intent(In   ) :: params

    Character(Len=MAX_KEY), Allocatable, Dimension(:) :: keys
    Integer                                           :: i
    Type(control_parameter)                           :: param

  Call params%get_keys(keys)

  Do i = 1, params%used_keys
    Call params%get(keys(i), param)
    Call write_control_param_help(param, error_unit)
  End Do

End Subroutine control_help_all

Subroutine dump_parameters(ifile, params, mode)
  Integer, Intent(In) :: ifile
  Class(parameters_hash_table), Intent(In) :: params
  Character(Len=10), Intent(In) :: mode
  Type(control_parameter) :: param
  Character(Len=MAX_KEY), Dimension(:), Allocatable :: keys
  Character(Len=*), Dimension(0:7), Parameter :: data_name = &
                                                 [Character(Len=7) :: 'NULL', 'INT', 'FLOAT', 'STRING', &
                                                                       'BOOL', 'OPTION', 'VECTOR3', 'VECTOR6']
  Character(Len=*), Dimension(0:7), Parameter :: python_data_name = &
                                                 [Character(Len=10) :: 'None', 'int', 'float', 'str', 'bool', &
                                                                        'str', '(float,)*3', '(float,)*6']
  Integer :: i

  Call params%get_keys(keys)

  Select Case (mode)
  Case ('latexdoc')
    Write (ifile, '(a)') '\documentclass{article}'
    Write (ifile, '(a)') '\usepackage[margin=1cm]{geometry}'
    Write (ifile, '(a)') '\usepackage{longtable}'
    Write (ifile, '(a)') '\begin{document}'
    Write (ifile, '(a)') '\begin{longtable}{l l p{10cm}}'

  Case ('latex')
    Write (ifile, '(a)') '\begin{longtable}{l l p{10cm}}'
  Case ('python')
    Write (ifile, '(a)') 'DLPData.__init__(self, {'
  Case ('csv', 'test')
    Continue
  Case ('default')
    Write (ifile, '(a)') 'Param | type | description | default value | unit '
  Case Default
    Call error(0, 'Bad mode option '//Trim(mode))
  End Select

  Do i = 1, params%used_keys
    Call params%get(keys(i), param)
    Select Case (mode)
    Case ('latex', 'latexdoc')
      ! Escape _
      param%key = escape(param%key)
      param%val = escape(param%val)
      param%description = escape(param%description)

      If (param%val == '') Then
        Write (ifile, '(2(a,1X,"&",1X), a, "\\")') &
          Trim(param%key), Trim(data_name(param%data_type)), &
          Trim(param%description)
      Else If (param%units == '') Then
        Write (ifile, '(2(a,1X,"&",1X), a, 1X, "(default = ", a, ") \\")') &
          Trim(param%key), Trim(data_name(param%data_type)), &
          Trim(param%description), Trim(param%val)
      Else
        Write (ifile, '(2(a,1X,"&",1X), a, 1X, "(default = ", a, 1X, "\verb#", a, "#) \\")') &
          Trim(param%key), Trim(data_name(param%data_type)), &
          Trim(param%description), Trim(param%val), Trim(param%units)
      End If
    Case ('csv')
      Write (ifile, '(5(a,";"))') &
        Trim(param%key), Trim(data_name(param%data_type)), &
        Trim(param%description), Trim(param%val), Trim(param%units)
    Case ('default')
      Write (ifile, '(5(a,"| "))') &
        Trim(param%key), Trim(data_name(param%data_type)), &
        Trim(param%description), Trim(param%val), Trim(param%units)
    Case ('python')
      If (Len_trim(param%units) > 0) Then
        Write (ifile, '("''",a, "'':",1X,"(",a,",",1X,a,"),")') Trim(param%key), Trim(python_data_name(param%data_type)), &
          Trim(python_data_name(DATA_STRING))
      Else
        Write (ifile, '("''",a,"'':",1X,a,",")') Trim(param%key), Trim(python_data_name(param%data_type))
      End If
    Case ('test')
      Select Case (param%data_type)
      Case (DATA_NULL)
        Continue
      Case (DATA_INT)
        param%val = "66666"
      Case (DATA_FLOAT)
        param%val = "6.666"
      Case (DATA_STRING, DATA_OPTION)
        Write (param%val, *) "JUNK"
      Case (DATA_BOOL)
        If (param%val == "on") Then
          param%val = "off"
        Else
          param%val = "on"
        End If
      Case (DATA_VECTOR3)
        param%val = "[ 6.666 6.666 6.666 ]"
      Case (DATA_VECTOR6)
        param%val = "[ 6.666 6.666 6.666 6.666 6.666 6.666 ]"
      End Select
      Write (ifile, '(3(a,1X))') Trim(param%key), Trim(param%val), Trim(param%units)
    End Select

  End Do

  Select Case (mode)
  Case ('latexdoc')
    Write (ifile, '(a)') '\end{longtable}'
    Write (ifile, '(a)') '\end{document}'
  Case ('latex')
    Write (ifile, '(a)') '\end{longtable}'
  Case ('python')
    Write (ifile, '(a)') '})'
  Case ('csv', 'test', 'default')
    Continue
  End Select

  flush(ifile)

Contains

  pure Function escape(string)
    Character(Len=*), Intent(In   ) :: string
    Character(Len=Len(string))      :: escape

    Character :: curr_char
    Integer   :: read_pos, write_pos

    write_pos = 1
    escape = ""
    Do read_pos = 1, Len_trim(string)
      curr_char = string(read_pos:read_pos)
      Select Case (curr_char)
      Case ("_")
        escape(write_pos:write_pos + 1) = "\"//curr_char !"
        write_pos = write_pos + 2
      Case Default
        escape(write_pos:write_pos) = curr_char
        write_pos = write_pos + 1
      End Select

    End Do

  End Function escape

End Subroutine dump_parameters

Subroutine print_set(params)
  Class(parameters_hash_table), Intent(In) :: params
  Type(control_parameter) :: param
  Character(Len=MAX_KEY), Dimension(:), Allocatable :: keys
  Integer :: i

  Call params%get_keys(keys)
  Do i = 1, params%used_keys
    Call params%get(keys(i), param)
    If (param%set) Print *, param
  End Do

End Subroutine print_set

Subroutine write_control_param_help(param, unit)
    Type(control_parameter), Intent(In   ) :: param
    Integer,                 Intent(In   ) :: unit

    Character(Len=*), Dimension(7), Parameter :: Type = [Character(Len=8) ::  "Int", "Real", &
                                                 "String", "Boolean ", "Option", "3-Vector", &
                                                 "6-Vector"]

  Write (unit, '(A,A)') "Key: ", Trim(param%key)
  Write (unit, '(A,A)') "Name: ", Trim(param%name)
  Write (unit, '(A,A,1X,A)') "Default: ", Trim(param%val), Trim(param%units)
  Write (unit, '(A,A)') "Description: ", Trim(param%description)
  Write (unit, '(A,A)') "Type: ", Trim(Type(param%data_type))
  Write (unit, *) ""

End Subroutine write_control_param_help

Subroutine write_control_param(param, unit, iotype, v_list, iostat, iomsg)
  !!-----------------------------------------------------------------------
  !!
  !! Print a friendly representation of the control parameter
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins april 2020
  !!-----------------------------------------------------------------------
    Class(control_parameter), Intent(In   ) :: param
    Integer,                  Intent(In   ) :: unit
    Character(Len=*),         Intent(In   ) :: iotype
    Integer, Dimension(:),    Intent(In   ) :: v_list
    Integer,                  Intent(  Out) :: iostat
    Character(Len=*),         Intent(Inout) :: iomsg

    Integer       :: itmp
    Logical       :: stat
    Real(Kind=wp) :: rtmp, rtmp3(3), rtmp6(6)

  Select Case (param%data_type)
  Case (DATA_FLOAT)
    Read (param%val, *) rtmp
    rtmp = convert_units(rtmp, param%units, param%internal_units, stat)
    If (.not. stat) Call error_units(param%units, param%internal_units, 'Cannot write '//param%key//': bad units')

    Write (unit, '(3(A,1X), "-> ", g15.6e2, 1X, A)') Trim(param%key), Trim(param%val), &
         & Trim(param%units), rtmp, Trim(param%internal_units)
  Case (DATA_INT)
    Read (param%val, *) itmp
    Write (unit, '(3(A,1X), "-> ", i0, 1X, A)') Trim(param%key), Trim(param%val), &
         & Trim(param%units), itmp, Trim(param%internal_units)
  Case (DATA_VECTOR3)
    Read (param%val, *) rtmp3
    Do itmp = 1, 3
      rtmp3(itmp) = convert_units(rtmp3(itmp), param%units, param%internal_units, stat)
      If (.not. stat) Call error_units(param%units, param%internal_units, 'Cannot write '//param%key//': bad units')

    End Do
    Write (unit, '(3(A,1X), "-> [", 3(g15.6e2,1X), "]", 1X, A)') Trim(param%key), Trim(param%val), &
         & Trim(param%units), rtmp3, Trim(param%internal_units)
  Case (DATA_VECTOR6)
    Read (param%val, *) rtmp6
    Do itmp = 1, 6
      rtmp6(itmp) = convert_units(rtmp6(itmp), param%units, param%internal_units, stat)
      If (.not. stat) Call error_units(param%units, param%internal_units, 'Cannot write '//param%key//': bad units')

    End Do

    Write (unit, '(3(A,1X), "-> [", 6(g15.6e2,1X), "]", 1X, A)') Trim(param%key), Trim(param%val), &
         & Trim(param%units), rtmp6, Trim(param%internal_units)
  Case default
    Write (unit, fmt='(3(A,1X))') Trim(param%key), Trim(param%val), Trim(param%units)
  End Select

End Subroutine write_control_param

Subroutine retrieve_option_or_string(table, key, output, required)
    Class(parameters_hash_table)          :: table
    Character(Len=*),       Intent(In   ) :: key
    Character(Len=STR_LEN), Intent(  Out) :: output
    Logical, Optional,      Intent(In   ) :: required

    Type(control_parameter) :: param

  Call table%get(key, param)
  If (Present(required)) Then
    If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
  End If
  output = param%val

End Subroutine retrieve_option_or_string

Subroutine retrieve_float(table, key, output, required)
    Class(parameters_hash_table)     :: table
    Character(Len=*),  Intent(In   ) :: key
    Real(kind=wp),     Intent(  Out) :: output
    Logical, Optional, Intent(In   ) :: required

    Character(Len=STR_LEN)  :: parse, val
    Logical                 :: stat
    Type(control_parameter) :: param

  Call table%get(key, param)
  If (Present(required)) Then
    If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
  End If
  val = param%val
  Call get_word(val, parse)

  output = word_2_real(parse)
  output = convert_units(output, param%units, param%internal_units, stat)
  If (.not. stat) Call error_units(param%units, param%internal_units, "When parsing key: "//Trim(key))

End Subroutine retrieve_float

Subroutine retrieve_vector_real(table, key, output, required)
    Class(parameters_hash_table)               :: table
    Character(Len=*),            Intent(In   ) :: key
    Real(kind=wp), Dimension(:), Intent(  Out) :: output
    Logical, Optional,           Intent(In   ) :: required

    Character(Len=STR_LEN)       :: parse, val
    Integer                      :: i
    Logical                      :: stat
    Real(kind=wp), Dimension(10) :: tmp
    Type(control_parameter)      :: param

  Call table%get(key, param)
  If (Present(required)) Then
    If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
  End If
  val = param%val

  Do i = 1, 10
    Call get_word(val, parse)
    If (parse == "") Exit
    tmp(i) = word_2_real(parse)
    tmp(i) = convert_units(tmp(i), param%units, param%internal_units, stat)
    If (.not. stat) Call error_units(param%units, param%internal_units, "When parsing key: "//Trim(key))

  End Do

  Select Case (param%data_type)
  Case (DATA_VECTOR3)
    If (Size(output) /= 3) Call error(0, "Bad length output vector")
    If (i /= 4) Call error(0, "Bad length input vector")
    output = tmp(1:3)

  Case (DATA_VECTOR6)
    If (Size(output) /= 6) Call error(0, "Bad length output vector")
    Select Case (i)
    Case (7)
      output = tmp(1:6)
    Case (10)
      output = [tmp(1), tmp(5), tmp(9), tmp(2), tmp(3), tmp(4)]
    Case default
      Call error(0, "Bad length input vector")
    End Select

  End Select

End Subroutine retrieve_vector_real

Subroutine retrieve_vector_int(table, key, output, required)
    Class(parameters_hash_table)         :: table
    Character(Len=*),      Intent(In   ) :: key
    Integer, Dimension(:), Intent(  Out) :: output
    Logical, Optional,     Intent(In   ) :: required

    Character(Len=STR_LEN)  :: parse, val
    Integer                 :: i
    Integer, Dimension(9)   :: tmp
    Type(control_parameter) :: param

  Call table%get(key, param)
  If (Present(required)) Then
    If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
  End If
  val = param%val

  Do i = 1, 9
    Call get_word(val, parse)
    If (parse == "") Exit
    tmp(i) = Nint(word_2_real(parse))
  End Do

  Select Case (param%data_type)
  Case (DATA_VECTOR3)
    If (Size(output) /= 3) Call error(0, "Bad length output vector")

    If (i /= 4) Call error(0, "Bad length input vector")
    output = tmp(1:3)

  Case (DATA_VECTOR6)
    If (Size(output) /= 6) Call error(0, "Bad length output vector")
    Select Case (i)
    Case (7)
      output = tmp(1:6)
    Case (10)
      output = [tmp(1), tmp(5), tmp(9), tmp(2), tmp(3), tmp(4)]
    Case default
      Call error(0, "Bad length input vector")
    End Select

  End Select

End Subroutine retrieve_vector_int

Subroutine retrieve_int(table, key, output, required)
    Class(parameters_hash_table)     :: table
    Character(Len=*),  Intent(In   ) :: key
    Integer,           Intent(  Out) :: output
    Logical, Optional, Intent(In   ) :: required

    Real(kind=wp)           :: rtmp
    Type(control_parameter) :: param

  Call table%get(key, param)
  If (Present(required)) Then
    If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
  End If

  Select Case (param%data_type)
  Case (DATA_INT)
    output = Nint(word_2_real(param%val))
  Case (DATA_FLOAT)
    If (param%internal_units /= 'steps') &
      Call error(0, 'Tried to parse physical value to int')
    Call table%retrieve(key, rtmp)
    output = Nint(rtmp)
  End Select

End Subroutine retrieve_int

Subroutine retrieve_bool(table, key, output, required)
    Class(parameters_hash_table)     :: table
    Character(Len=*),  Intent(In   ) :: key
    Logical,           Intent(  Out) :: output
    Logical, Optional, Intent(In   ) :: required

    Type(control_parameter) :: param

  Call table%get(key, param)
  If (Present(required)) Then
    If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
  End If

  Select Case (param%val)
  Case ('on', 'y')
    output = .true.
  Case Default
    output = .false.
  End Select

End Subroutine retrieve_bool

Subroutine get_param(table, key, val, default)
    Class(parameters_hash_table),      Intent(In   ) :: table
    Character(Len=*),                  Intent(In   ) :: key
    Type(control_parameter),           Intent(  Out) :: val
    Type(control_parameter), Optional, Intent(In   ) :: default

    Class(*), Pointer :: stuff

  Call table%get_cont(key, default, stuff)

  Select Type (stuff)
  Type is (control_parameter)
    val = stuff
  Class Default
    Call error(0, 'Trying to get control_param from a not control_param')
  End Select
  Deallocate (stuff)
  Nullify (stuff)

End Subroutine get_param

Function is_set_single(table, key) Result(is_set)
    Class(parameters_hash_table), Intent(In   ) :: table
    Character(Len=*),             Intent(In   ) :: key
    Logical                                     :: is_set

    Type(control_parameter) :: param

  Call table%get_param(key, param)
  is_set = param%set

End Function is_set_single

Function is_all_set(table, key) Result(is_set)
    Class(parameters_hash_table),   Intent(In   ) :: table
    Character(Len=*), Dimension(:), Intent(In   ) :: key
    Logical                                       :: is_set

    Character(Len=STR_LEN)  :: curr_key
    Integer                 :: i
    Type(control_parameter) :: param

  is_set = .true.
  Do i = 1, Size(key)
    curr_key = key(i)
    Call table%get_param(curr_key, param)
    is_set = is_set .and. param%set
  End Do

End Function is_all_set

Function is_any_set(table, key) Result(is_set)
    Class(parameters_hash_table),   Intent(In   ) :: table
    Character(Len=*), Dimension(:), Intent(In   ) :: key
    Logical                                       :: is_set

    Character(Len=STR_LEN)  :: curr_key
    Integer                 :: i
    Type(control_parameter) :: param

  is_set = .false.
  Do i = 1, Size(key)
    curr_key = key(i)
    Call table%get_param(curr_key, param)
    is_set = is_set .or. param%set
  End Do

End Function is_any_set

Function num_set(table, key)
    Class(parameters_hash_table),   Intent(In   ) :: table
    Character(Len=*), Dimension(:), Intent(In   ) :: key
    Integer                                       :: num_set

    Character(Len=STR_LEN)  :: curr_key
    Integer                 :: i
    Type(control_parameter) :: param

  num_set = 0
  Do i = 1, Size(key)
    curr_key = key(i)
    Call table%get_param(curr_key, param)
    If (param%set) num_set = num_set + 1
  End Do

End Function num_set

End Module control_parameters
