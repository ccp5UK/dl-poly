Module control_parameters
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle control parameters in dlpoly
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins april 2020
  !!-----------------------------------------------------------------------

  Use errors_warnings,               Only: error,&
                                           error_units,&
                                           info
  Use hash,                          Only: MAX_KEY,&
                                           hash_table
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use kinds,                         Only: STR_LEN,&
                                           wp
  Use parse,                         Only: get_word,&
                                           word_2_real
  Use units,                         Only: convert_units,&
                                           to_out_units

  Implicit None

  Private

  !> Data types enumeration
  Integer, Parameter, Public :: DATA_NULL = 0, DATA_INT = 1, DATA_FLOAT = 2, DATA_STRING = 3, &
       DATA_BOOL = 4, DATA_OPTION = 5, DATA_VECTOR = 6, STRING_VECTOR = 7


  Type, Public, Extends(hash_table) :: parameters_hash_table
  Contains
    !> Update get to include params
    Generic, Public  :: get => get_param
    Procedure, Private :: get_param
    !> Set retrieve up to parse stored params
    Generic, Public  :: retrieve => retrieve_option_or_string, retrieve_float, &
         & retrieve_vector_real, retrieve_vector_int, retrieve_vector_string, retrieve_int, retrieve_bool
    Procedure, Pass(table), Private :: retrieve_option_or_string, retrieve_float, &
         & retrieve_vector_real, retrieve_vector_int, retrieve_vector_string, retrieve_int, retrieve_bool
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
    !> Control parameter data type (int, float, VECTOR, stringVECTOR, string, bool, option)
    Integer :: data_type = DATA_NULL
    !> length, 1 for scalars, possibly >1 for data_VECTOR and string_VECTOR
    Integer :: length = 0
    !> Is value set
    Logical :: set = .false.
  Contains
    Procedure, Private :: write_control_param
    Generic :: Write (formatted) => write_control_param
  End Type control_parameter

  Character(Len=5), Parameter :: INDENT_STR = '  -- '

  Public :: dump_parameters
  Public :: control_help_single, control_help_all
  Public :: print_set
  Public :: write_param

  Interface write_param
    module procedure write_single_param_logical
    module procedure write_single_param_unit
    module procedure write_single_param_int
    module procedure write_single_param_str
  End Interface write_param

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
         [Character(Len=13) :: 'NULL', 'INT', 'FLOAT', 'STRING', &
         'BOOL', 'OPTION', 'VECTOR', 'STRING_VECTOR']
    Character(Len=*), Dimension(0:7), Parameter :: python_data_name = &
         [Character(Len=10) :: 'None', 'int', 'float', 'str', 'bool', &
         'str', '(float,)*N', '(Char,)*N']
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
        Case (DATA_VECTOR)
          param%val = "[ 6.666 6.666 6.666 ...]"
        Case (STRING_VECTOR)
          param%val = "[ "" "" "" ... ]"
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

    Integer                     :: itmp
    Logical                     :: stat
    Real(Kind=wp)               :: rtmp
    Real(Kind=wp), Allocatable  :: rtmpN(:)

    Select Case (param%data_type)
      Case (DATA_FLOAT)
        Read (param%val, *) rtmp
        rtmp = convert_units(rtmp, param%units, param%internal_units, stat)
        If (.not. stat) Call error_units(param%units, param%internal_units, 'Cannot write '//param%key//': bad units')

        Write (unit, '(3(A,1X), "-> ", g15.6e2, 1X, A)', iostat=iostat, iomsg=iomsg) Trim(param%key), Trim(param%val), &
        & Trim(param%units), rtmp, Trim(param%internal_units)
      Case (DATA_INT)
        Read (param%val, *) itmp
        Write (unit, '(3(A,1X), "-> ", i0, 1X, A)',iostat=iostat, iomsg=iomsg) Trim(param%key), Trim(param%val), &
        & Trim(param%units), itmp, Trim(param%internal_units)

      Case (DATA_VECTOR)
        If (Allocated(rtmpN)) Then 
          Deallocate(rtmpN)
          Allocate(rtmpN(1:param%length))
        End If

        Read (param%val, *) rtmpN
        Do itmp = 1, param%length
          rtmpN(itmp) = convert_units(rtmpN(itmp), param%units, param%internal_units, stat)
        End Do
        Write (unit, '(3(A,1X), "-> [", 3(g15.6e2,1X), "]", 1X, A)', iostat=iostat, iomsg=iomsg) Trim(param%key), Trim(param%val), &
        & Trim(param%units), rtmpN, Trim(param%internal_units)

       Case default
         Write (unit, fmt='(3(A,1X))', iostat=iostat, iomsg=iomsg) Trim(param%key), Trim(param%val), Trim(param%units)
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

  Subroutine retrieve_vector_real(table, key, output, dimension, required)
    Class(parameters_hash_table)               :: table
    Character(Len=*),            Intent(In   ) :: key
    Real(kind=wp), Allocatable,  Intent(  Out) :: output(:)
    Integer, Optional,           Intent(In   ) :: dimension
    Logical, Optional,           Intent(In   ) :: required

    Character(Len=STR_LEN)       :: parse, val
    Integer                      :: i, length
    Logical                      :: stat
    Type(control_parameter)      :: param

    Call table%get(key, param)
    If (Present(required)) Then
      If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
    End If
    val = param%val

    If (Allocated(output)) Then 
      Deallocate(output)
    End If

    ! get length
    length = 0
    Do While (parse /= "")
      Call get_word(val, parse)
      If (parse /= "") Then 
        length = length + 1
      End If
    End Do

    If (Present(dimension)) Then

      If (length == 0) Then
        Call error(0,"Zero length input vector")
      End If 

      If (length /= dimension) Then
        Call error(0,"input vector length not equal to request vector length")
      End If
      
    End If

    val = param%val

    param%length = length

    Allocate(output(1:length))

    Do i = 1, length
      Call get_word(val, parse)
      If (parse == "") Exit
      output(i) = word_2_real(parse)
      output(i) = convert_units(output(i), param%units, param%internal_units, stat)
      If (.not. stat) Call error_units(param%units, param%internal_units, "When parsing key: "//Trim(key))

    End Do

  End Subroutine retrieve_vector_real

  Subroutine retrieve_vector_string(table, key, output, dimension, required)
    Class(parameters_hash_table)                        :: table
    Character(Len=*),                     Intent(In   ) :: key
    Character(Len=STR_LEN), Allocatable,  Intent(  Out) :: output(:)
    Integer, Optional,                    Intent(In   ) :: dimension
    Logical, Optional,                    Intent(In   ) :: required

    Character(Len=STR_LEN)                              :: parse, val
    Integer                                             :: i, length
    Type(control_parameter)                             :: param

    Call table%get(key, param)
    If (Present(required)) Then
      If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
    End If
    val = param%val

    If (Allocated(output)) Then 
      Deallocate(output)
    End If

    ! get length
    length = 0
    Do While (parse /= "")
      Call get_word(val, parse)
      If (parse /= "") Then 
        length = length + 1
      End If
    End Do

    If (Present(dimension)) Then

      If (length == 0) Then
        Call error(0,"Zero length input vector")
      End If 

      If (length /= dimension) Then
        Call error(0,"input vector length not equal to request vector length")
      End If
      
    End If

    val = param%val
    param%length = length

    Allocate(output(1:length))

    Do i = 1, length
      Call get_word(val, parse)
      If (parse == "") Exit
      output(i) = Trim(parse)
    End Do

  End Subroutine retrieve_vector_string

  Subroutine retrieve_vector_int(table, key, output, dimension, required)
    Class(parameters_hash_table)         :: table
    Character(Len=*),      Intent(In   ) :: key
    Integer, Allocatable,  Intent(  Out) :: output(:)
    Integer, Optional,     Intent(In   ) :: dimension
    Logical, Optional,     Intent(In   ) :: required

    Character(Len=STR_LEN)  :: parse, val
    Integer                 :: i, length
    Type(control_parameter) :: param

    Call table%get(key, param)
    If (Present(required)) Then
      If (required .and. .not. param%set) Call error(0, 'Necessary parameter '//Trim(key)//' not set')
    End If
    val = param%val

    If (Allocated(output)) Then 
      Deallocate(output)
    End If

    ! get length
    length = 0
    Do While (parse /= "")
      Call get_word(val, parse)
      If (parse /= "") Then 
        length = length + 1
      End If
    End Do
    
    If (Present(dimension)) Then

      If (length == 0) Then
        Call error(0,"Zero length input vector")
      End If 

      If (length /= dimension) Then
        Call error(0,"input vector length not equal to request vector length")
      End If

    End If

    val = param%val
    param%length = length

    Allocate(output(1:length))

    Do i = 1, length
      Call get_word(val, parse)
      If (parse == "") Exit
      output(i) = Nint(word_2_real(parse))
    End Do

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

  Function get_indent_str(indent) result(indent_str)
    Integer :: indent_str
    Integer, Intent(In), Optional :: indent
    Integer :: indent_in

    indent_in = 0
    If (Present(indent)) Then
      indent_in = indent
    End If

    Select Case (indent_in)
    Case (0)
      indent_str = 0
    Case (1)
      indent_str = 2
    case (2)
      indent_str = 5
    End Select

  End Function get_indent_str

  Subroutine write_single_param_unit(name, val, from, to, indent, level)
    !!-----------------------------------------------------------------------
    !!
    !! Return a string ready to be output in standard form
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins july 2022
    !!-----------------------------------------------------------------------
    Character(Len=*), Intent(In) :: name
    Real(kind=wp), Intent(In) :: val
    Character(Len=*), Intent(In), Optional :: from, to
    Integer, Intent(In), Optional :: indent
    Integer, Intent(In), Optional :: level
    Character(Len=STR_LEN) :: output
    Real(kind=wp) :: res


    Character(Len=STR_LEN) :: out_unit, unit_str

    If (Present(from) .and. Present(to)) Then
      res = convert_units(val, from, to)
      unit_str = ' (' // trim(to) // ')'
    Else If (Present(from)) Then
      Call to_out_units(val, from, res, out_unit)
      unit_str = ' (' // trim(out_unit) // ')'
    Else
      unit_str = ' '
      res = val
    End If

    Write(output, '(1A,1A,1A,": ",G12.5E2)') INDENT_STR(1:get_indent_str(indent)), trim(name), trim(unit_str), res
    Call info(output, .true., level)

  End Subroutine write_single_param_unit

  Subroutine write_single_param_logical(name, val, indent, off_level, on_level)
    !!-----------------------------------------------------------------------
    !!
    !! Return a string ready to be output in standard form
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins july 2022
    !!-----------------------------------------------------------------------
    Character(Len=*), Intent(In) :: name
    Logical, Intent(In) :: val
    Integer, Intent(In), Optional :: indent
    Integer, Intent(In), Optional :: off_level, on_level
    Character(Len=STR_LEN) :: output

    Character(Len=5) :: intent_str



    if (val) then
      Write(output, '(1A,1A,": ",1A)') INDENT_STR(1:get_indent_str(indent)), trim(name), 'ON'
      Call info(output, .true., on_level)
    else

      Write(output, '(1A,1A,": ",1A)') INDENT_STR(1:get_indent_str(indent)), trim(name), 'OFF'
      Call info(output, .true., off_level)
    end if


  End Subroutine write_single_param_logical

  Subroutine write_single_param_int(name, val, indent, level)
    !!-----------------------------------------------------------------------
    !!
    !! Return a string ready to be output in standard form
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins july 2022
    !!-----------------------------------------------------------------------
    Character(Len=*), Intent(In) :: name
    Integer, Intent(In) :: val
    Integer, Intent(In), Optional :: indent
    Integer, Intent(In), Optional :: level
    Character(Len=STR_LEN) :: output

    Write(output, '(1A,1A,": ",i0.1)') INDENT_STR(1:get_indent_str(indent)), trim(name), val
    Call info(output, .true., level)

  End Subroutine write_single_param_int

  Subroutine write_single_param_str(name, val, indent, level)
    !!-----------------------------------------------------------------------
    !!
    !! Return a string ready to be output in standard form
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins july 2022
    !!-----------------------------------------------------------------------
    Character(Len=*), Intent(In) :: name
    Character(Len=*), Intent(In) :: val
    Integer, Intent(In), Optional :: indent
    Integer, Intent(In), Optional :: level
    Character(Len=STR_LEN) :: output

    Write(output, '(1A,1A,": ",1A)') INDENT_STR(1:get_indent_str(indent)), trim(name), trim(val)
    Call info(output, .true., level)

  End Subroutine write_single_param_str


End Module control_parameters
