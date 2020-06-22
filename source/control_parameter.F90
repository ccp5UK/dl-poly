module control_parameter_module

  Use units, only : convert_units
  Use parse, only : word_2_real, get_word, get_line
  Use kinds, only : wp
  Use hash, only : hash_table, MAX_KEY, STR_LEN
  Use errors_warnings, only : error
  Implicit None

  Private

  !> Data types enumeration
  Integer, Parameter, Public :: DATA_NULL=0, DATA_INT=1, DATA_FLOAT=2, DATA_STRING=3, &
       DATA_BOOL=4, DATA_OPTION=5, DATA_VECTOR3=6, DATA_VECTOR6=7

  Type, Public, Extends(hash_table) :: parameters_hash_table
   contains
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
  end type parameters_hash_table

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
   contains
     Procedure, Private :: write_control_param
     Generic :: write(formatted) => write_control_param
  End Type control_parameter

  Public :: dump_parameters
  Public :: control_help_single, control_help_all
  Public :: print_set

contains

  Subroutine control_help_single(params, key)
    Class( parameters_hash_table ), intent( In ) :: params
    Character(Len=*), Intent( In ) :: key
    Type (control_parameter) :: param

    call params%get(key, param)
    call write_control_param_help(param, 0)

  End Subroutine control_help_single

  Subroutine control_help_all(params)
    Class( parameters_hash_table ), intent( In ) :: params
    Type (control_parameter) :: param
    Character(Len=MAX_KEY), Dimension(:), Allocatable :: keys
    Integer :: i

    call params%get_keys(keys)

    do i = 1, params%used_keys
       call params%get(keys(i), param)
       call write_control_param_help(param, 0)
    end do

  End Subroutine control_help_all

  Subroutine dump_parameters(ifile, params, mode)
    Integer, Intent(In) :: ifile
    Class(parameters_hash_table), Intent(In) :: params
    Character(Len=10), Intent(In), Value :: mode
    Type (control_parameter) :: param
    Character(Len=MAX_KEY), Dimension(:), Allocatable :: keys
    Character(Len=*), Dimension(0:7), Parameter :: data_name = &
         [Character(Len=7) :: 'NULL', 'INT', 'FLOAT', 'STRING', 'BOOL', 'OPTION', 'VECTOR3', 'VECTOR6']
    Integer :: i

    call params%get_keys(keys)

    Select Case (mode)
    Case ('latexdoc')
      Write(ifile, '(a)') '\documentclass{article}'
      Write(ifile, '(a)') '\usepackage[margin=1cm]{geometry}'
      Write(ifile, '(a)') '\usepackage{longtable}'
      Write(ifile, '(a)') '\begin{document}'
      Write(ifile, '(a)') '\begin{longtable}{l l p{10cm}}'

    Case ('latex')
      Write(ifile, '(a)') '\begin{longtable}{l l p{10cm}}'
    Case ('csv')
      Continue
    Case Default
      Call error(0, 'Bad mode option '//trim(mode))
    end Select

    do i = 1, params%used_keys
       call params%get(keys(i), param)
       Select Case (mode)
       Case ('latex', 'latexdoc')
         ! Escape _
         param%key = escape(param%key)
         param%val = escape(param%val)
         param%units = escape(param%units)
         param%description = escape(param%description)

         if (param%val == '') then
           Write(ifile, '(2(a,1X,"&",1X), a, "\\")') &
                trim(param%key), trim(data_name(param%data_type)), &
                trim(param%description)
         else if (param%units == '') then
           Write(ifile, '(2(a,1X,"&",1X), a, 1X, "(default = ", a, ") \\")') &
                trim(param%key), trim(data_name(param%data_type)), &
                trim(param%description), trim(param%val)
         else
           Write(ifile, '(2(a,1X,"&",1X), a, 1X, "(default = ", a, 1X, "\verb#", a, "#) \\")') &
                trim(param%key), trim(data_name(param%data_type)), &
                trim(param%description), trim(param%val), trim(param%units)
         end if
       Case ('csv')
         Write(ifile, '(5(a,";"))') &
              trim(param%key), trim(data_name(param%data_type)), &
              trim(param%description), trim(param%val), trim(param%units)
       end Select

    end do

    Select Case (mode)
    Case ('latexdoc')
      Write(ifile, '(a)') '\end{longtable}'
      Write(ifile, '(a)') '\end{document}'
    Case ('latex')
      Write(ifile, '(a)') '\end{longtable}'
    Case ('csv')
      Continue
    end Select

  contains

    pure function escape(string)
      Character(Len=*), Intent(In) :: string
      Character(Len=len(string)) :: escape

      Integer :: read_pos, write_pos
      Character :: curr_char

      write_pos = 1
      escape = ""
      do read_pos = 1, len_trim(string)
        curr_char = string(read_pos:read_pos)
        Select Case(curr_char)
        Case("_")
          escape(write_pos:write_pos+1) = "\"//curr_char !"
          write_pos = write_pos + 2
        Case Default
          escape(write_pos:write_pos) = curr_char
          write_pos = write_pos + 1
        end Select

      end do

    end function escape

  end Subroutine dump_parameters


  Subroutine print_set(params)
    Class( parameters_hash_table ), intent( In ) :: params
    Type (control_parameter) :: param
    Character(Len=MAX_KEY), Dimension(:), Allocatable :: keys
    Integer :: i

    call params%get_keys(keys)
    do i = 1, params%used_keys
       call params%get(keys(i), param)
       if (param%set) print*, param
    end do

  end Subroutine print_set

  Subroutine write_control_param_help(param, unit)
    Type (control_parameter), Intent(In) :: param
    Character(Len=*), Dimension(7), Parameter :: type = [Character(Len=8) :: &
         & "Int", "Real", "String", "Boolean ", "Option", "3-Vector", "6-Vector"]
    Integer, Intent(In) :: unit

    write(unit, '(A,A)') "Key: ", trim(param%key)
    write(unit, '(A,A)') "Name: ", trim(param%name)
    write(unit, '(A,A,1X,A)') "Default: ", trim(param%val), trim(param%units)
    write(unit, '(A,A)') "Description: ", trim(param%description)
    write(unit, '(A)') trim(type(param%data_type))
    write(unit, *) ""

  end Subroutine write_control_param_help

  Subroutine write_control_param(param, unit, iotype, v_list, iostat, iomsg)
    !!-----------------------------------------------------------------------
    !!
    !! Print a friendly representation of the control parameter
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Class (control_parameter), Intent(In) :: param
    Integer, Intent(In) :: unit
    Character (Len=*), Intent(In) :: iotype
    Integer, Intent(In), Dimension(:) :: v_list
    Integer, Intent(Out) :: iostat
    Real(Kind=wp) :: rtmp, rtmp3(3), rtmp6(6)
    Integer :: itmp
    Character (Len=*), Intent(Inout) :: iomsg

    select case (param%data_type)
    case(DATA_FLOAT)
       read(param%val, *) rtmp
       rtmp = convert_units(rtmp, param%units, param%internal_units)
       write(unit, '(3(A,1X), "-> ", g15.6e2, 1X, A)') trim(param%key), trim(param%val), &
            & trim(param%units), rtmp, trim(param%internal_units)
    case(DATA_INT)
       read(param%val, *) itmp
       write(unit, '(3(A,1X), "-> ", i0, 1X, A)') trim(param%key), trim(param%val), &
            & trim(param%units), rtmp, trim(param%internal_units)
    case(DATA_VECTOR3)
       read(param%val, *) rtmp3
       do itmp = 1,3
          rtmp3(itmp) = convert_units(rtmp3(itmp), param%units, param%internal_units)
       end do
       write(unit, '(3(A,1X), "-> [", 3(g15.6e2,1X), "]", 1X, A)') trim(param%key), trim(param%val), &
            & trim(param%units), rtmp, trim(param%internal_units)
    case(DATA_VECTOR6)
       read(param%val, *) rtmp6
       do itmp = 1,6
          rtmp6(itmp) = convert_units(rtmp6(itmp), param%units, param%internal_units)
       end do

       write(unit, '(3(A,1X), "-> [", 6(g15.6e2,1X), "]", 1X, A)') trim(param%key), trim(param%val), &
            & trim(param%units), rtmp, trim(param%internal_units)
    case default
       write(unit, fmt='(3(A,1X))') trim(param%key), trim(param%val), trim(param%units)
    end select

  end Subroutine write_control_param

  Subroutine retrieve_option_or_string(table, key, output, required)
    Class( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Character(Len=STR_LEN), Intent( Out    ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    output = param%val

  end Subroutine retrieve_option_or_string

  Subroutine retrieve_float(table, key, output, required)
    Class( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse, val
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Real(kind=wp), Intent( Out    ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    val = param%val
    call get_word(val, parse)
    output = word_2_real(parse)
    output = convert_units(output, param%units, param%internal_units)

  End Subroutine retrieve_float

  Subroutine retrieve_vector_real(table, key, output, required)
    Class( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse, val
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Real(kind=wp), dimension(9) :: tmp
    Real(kind=wp), dimension(:), Intent( Out    ) :: output

    Integer :: i

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    val = param%val

    do i = 1, 4
       call get_word(val, parse)
       if (parse == "") exit
       tmp(i) = word_2_real(parse)
       tmp(i) = convert_units(tmp(i), param%units, param%internal_units)
    end do

    select case(param%data_type)
    case (DATA_VECTOR3)
       if (size(output) /= 3) call error(0, "Bad length output vector")

       if (i /= 4) call error(0, "Bad length input vector")
       output = tmp(1:3)

    case (DATA_VECTOR6)
       if (size(output) /= 6) call error(0, "Bad length output vector")
       select case (i)
       case(7)
          output = tmp(1:6)
       case(10)
          output = [tmp(1), tmp(5), tmp(9), tmp(2), tmp(3), tmp(4)]
       case default
          call error(0, "Bad length input vector")
       end select

    end select

  End Subroutine retrieve_vector_real

  Subroutine retrieve_vector_int(table, key, output, required)
    Class( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse, val
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Integer, dimension(9) :: tmp
    Integer, dimension(:), Intent( Out    ) :: output

    Integer :: i

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    val = param%val

    do i = 1, 9
       call get_word(val, parse)
       if (parse == "") exit
       tmp(i) = nint(word_2_real(parse))
    end do

    select case(param%data_type)
    case (DATA_VECTOR3)
       if (size(output) /= 3) call error(0, "Bad length output vector")

       if (i /= 4) call error(0, "Bad length input vector")
       output = tmp(1:3)

    case (DATA_VECTOR6)
       if (size(output) /= 6) call error(0, "Bad length output vector")
       select case (i)
       case(7)
          output = tmp(1:6)
       case(10)
          output = [tmp(1), tmp(5), tmp(9), tmp(2), tmp(3), tmp(4)]
       case default
          call error(0, "Bad length input vector")
       end select

    end select

  End Subroutine retrieve_vector_int

  Subroutine retrieve_int(table, key, output, required)
    Class( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type( control_parameter ) :: param
    Real( kind = wp ) :: rtmp
    Logical, Intent( In    ), Optional :: required
    Integer, Intent(   Out ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if

    select case (param%data_type)
    case (DATA_INT)
       output = nint(word_2_real(param%val))
    case (DATA_FLOAT)
       if (param%internal_units /= 'steps') &
            call error(0, 'Tried to parse physical value to int')
       call table%retrieve(key, rtmp)
       output = nint(rtmp)
    end select

  End Subroutine retrieve_int

  Subroutine retrieve_bool(table, key, output, required)
    Class( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Logical, Intent(   Out ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if

    Select Case(param%val)
    Case('on', 'y')
       output = .true.
    Case Default
       output = .false.
    end Select

  End Subroutine retrieve_bool

  Subroutine get_param( table, key, val, default )
    Class( parameters_hash_table ), Intent( In    ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type(control_parameter), Intent(   Out ) :: val
    Type(control_parameter), Intent( In    ), Optional :: default
    Class( * ), Allocatable :: stuff

    stuff = table%get_cont(key, default)

    Select Type( stuff )
    Type is ( control_parameter )
       val = stuff
    Class Default
       Call error(0, 'Trying to get control_param from a not control_param')
    End Select

  End Subroutine get_param

  Function is_set_single(table, key) result(is_set)
    Class( parameters_hash_table ), Intent( In    ) :: table
    Character(Len=*), Intent( In    ) :: key
    Logical :: is_set
    Type( control_parameter ) :: param

    call table%get_param(key, param)
    is_set = param%set

  end Function is_set_single

  Function is_all_set(table, key) result(is_set)
    Class( parameters_hash_table ), Intent( In    ) :: table
    Character(Len=*), dimension(:), Intent( In    ) :: key
    Character(Len=STR_LEN) :: curr_key
    Logical :: is_set
    Type( control_parameter ) :: param
    Integer :: i

    is_set = .true.
    do i = 1, size(key)
       curr_key = key(i)
       call table%get_param(curr_key, param)
       is_set = is_set .and. param%set
    end do

  end Function is_all_set

  Function is_any_set(table, key) result(is_set)
    Class( parameters_hash_table ), Intent( In    ) :: table
    Character(Len=*), dimension(:), Intent( In    ) :: key
    Character(Len=STR_LEN) :: curr_key
    Logical :: is_set
    Type( control_parameter ) :: param
    Integer :: i

    is_set = .false.
    do i = 1, size(key)
       curr_key = key(i)
       call table%get_param(curr_key, param)
       is_set = is_set .or. param%set
    end do

  end Function is_any_set

  Function num_set(table, key)
    Class( parameters_hash_table ), Intent( In    ) :: table
    Character(Len=*), dimension(:), Intent( In    ) :: key
    Character(Len=STR_LEN) :: curr_key
    Integer :: num_set
    Type( control_parameter ) :: param
    Integer :: i

    num_set = 0
    do i = 1, size(key)
       curr_key = key(i)
       call table%get_param(curr_key, param)
       if (param%set) num_set = num_set + 1
    end do

  end Function num_set

end module control_parameter_module
