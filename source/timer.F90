Module timer
  !!------------------------------------------------!
  !!
  !! dl_poly_4 module containing timing routines
  !!
  !! copyright - daresbury laboratory
  !! author    - j.s.wilkins february 2019
  !!
  !!------------------------------------------------!
  Use comms, Only: comms_type,&
                   gmax,&
                   gmin,&
                   grecv,&
                   gsend,&
                   gsum,&
                   gsync,&
                   gtime,&
                   mtime,&
                   timer_tag
  Use, Intrinsic :: iso_fortran_env, Only: eu => error_unit
  Use kinds, Only: wp

  Implicit None

  Private

  Integer, Parameter :: max_depth = 6, max_name = 18

  Type :: node_timer
    !!------------------------------------------------!
    !! Timer
    !!------------------------------------------------!
    Character(Len=max_name) :: name
    Integer           :: id
    Real(Kind=wp) :: max, min, total, last
    Real(Kind=wp) :: start, Stop
    Integer           :: calls
    Logical           :: running = .false.
  End Type node_timer

  Type :: call_stack
    !!------------------------------------------------!
    !! Call stack
    !!------------------------------------------------!
    Character(Len=max_name), Dimension(max_depth) :: name
    Integer :: depth = 0
  End Type call_stack

  Type :: node
    !!------------------------------------------------!
    !! Tree node
    !!------------------------------------------------!
    Type(node_timer) :: time
    Type(timer_tree), Pointer :: tree => null()
    Type(node), Pointer :: child => null()
    Type(node), Pointer :: parent => null()
    Type(node), Pointer :: next_sibling => null()
  End Type node

  Type :: timer_tree
    !!------------------------------------------------!
    !! Tree structure
    !!------------------------------------------------!
    Type(node), Pointer :: head => null()
    Integer :: n_timers = 0
  End Type timer_tree

  Type, Public :: timer_type
    !!------------------------------------------------!
    !! Main timer system
    !!------------------------------------------------!
    Type(timer_tree), Pointer :: tree
    Type(call_stack) :: stack
    Real(Kind=wp) :: elapsed, job, clear_screen
    Logical :: proc_detail = .false.
    Integer :: max_depth = 1
    Integer :: proc_id
    Integer :: out_unit
  End Type timer_type

  Interface timer_write
    Module Procedure timer_write_sing
    Module Procedure timer_write_mul
  End Interface timer_write

  Public :: timer_report
  Public :: timer_last_time
  Public :: start_timer
  Public :: stop_timer
  Public :: init_timer_system
  Public :: dump_call_stack
  Public :: start_timer_path
  Public :: stop_timer_path

  Public :: time_elapsed

Contains

  Subroutine dump_call_stack(stack)
    !!------------------------------------------------!
    !!
    !! Print out a call stack to dump current location
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(call_stack) :: stack

    Integer :: i

    Call timer_write('')
    Call timer_write('Process stack:')
    Do i = 1, stack%depth
      Call timer_write(stack%name(i))
    End Do

  End Subroutine dump_call_stack

  Subroutine push_stack(tmr, name)
    !!------------------------------------------------!
    !!
    !! Add a timer to the call stack
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Character(Len=*), Intent(In   ) :: name

    tmr%stack%depth = tmr%stack%depth + 1
    If (tmr%stack%depth > max_depth) &
      & Call timer_error(tmr, 'Call stack exceeds max depth : recursive call or unended timer?')
    tmr%stack%name(tmr%stack%depth) = name

  End Subroutine push_stack

  Subroutine pop_stack(tmr, name)
    !!------------------------------------------------!
    !!
    !! Remove a timer from the call stack
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Character(Len=*), Intent(In   ) :: name

    If (name /= tmr%stack%name(tmr%stack%depth)) &
      & Call timer_error(tmr, 'Child timer '//name//' ended before parent')
    tmr%stack%name(tmr%stack%depth) = ''
    tmr%stack%depth = tmr%stack%depth - 1

  End Subroutine pop_stack

  Subroutine init_timer_system(tmr, nrite, comm)
    !!------------------------------------------------!
    !!
    !! Initialise a timer system
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(inout) :: tmr
    Integer,          Intent(In   ) :: nrite
    Type(comms_type), Intent(In   ) :: comm

    tmr%proc_id = comm%idnode
    tmr%out_unit = nrite
    Allocate (tmr%tree)
    Allocate (tmr%tree%head)
    tmr%tree%head%tree => tmr%tree
    Call init_timer(tmr%tree%head%time, 'Head')
    Call start_timer(tmr, 'Main')

  End Subroutine init_timer_system

  Function find_timer(tmr, name, stack_in) Result(current)
    !!------------------------------------------------!
    !!
    !! Locate a timer node witin a given timer system
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Character(Len=*)                :: name
    Type(call_stack), Optional      :: stack_in
    Type(node), Pointer             :: current

    Integer          :: depth
    Type(call_stack) :: stack

    If (Present(stack_in)) Then
      stack = stack_in
    Else
      stack = tmr%stack
    End If

    current => tmr%tree%head
    depth = 0

    Do While (depth < stack%depth)
      If (.not. Associated(current%child)) Call timer_error(tmr, 'Call stack does not match call tree (no child)')
      depth = depth + 1
      current => current%child
      Do While (current%time%name /= stack%name(depth))
        If (.not. Associated(current%next_sibling)) &
          & Call timer_error(tmr, 'Call stack does not match call tree (no sibling)')
        current => current%next_sibling
      End Do
    End Do

    If (current%time%name == name) Then
      Continue
    Else If (.not. Associated(current%child)) Then
      Call init_child_node(name, current)
      current => current%child
      Return
    Else
      current => current%child
      Do While (current%time%name /= name)
        If (.not. Associated(current%next_sibling)) Then
          Call init_sibling_node(name, current)
          current => current%next_sibling
          Return
        End If
        current => current%next_sibling
      End Do
    End If

  End Function find_timer

  Subroutine init_child_node(name, parent)
    !!------------------------------------------------!
    !!
    !! Create a timer node which is a child of the parent node
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character(len=*), Intent(In   ) :: name
    Type(node), Target              :: parent

    Allocate (parent%child)
    parent%child%parent => parent
    parent%child%tree => parent%tree
    parent%tree%n_timers = parent%tree%n_timers + 1
    Call init_timer(parent%child%time, name)

  End Subroutine init_child_node

  Subroutine init_sibling_node(name, sibling)
    !!------------------------------------------------!
    !!
    !! Create a timer node which is a child of the parent node
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character(len=*) :: name
    Type(node)       :: sibling

    Type(node), Pointer :: child

    Allocate (sibling%next_sibling)
    child => sibling%next_sibling
    child%parent => sibling%parent

    sibling%next_sibling%tree => sibling%tree
    sibling%tree%n_timers = sibling%tree%n_timers + 1

    Call init_timer(child%time, name)

  End Subroutine init_sibling_node

  Subroutine start_timer(tmr, name, stack)
    !!------------------------------------------------!
    !!
    !! Start a timer running on a given timer system
    !! If stack is supplied bypass standard call stack
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Character(Len=*), Intent(In   ) :: name
    Type(call_stack), Optional      :: stack

    Type(node), Pointer :: current_timer

    current_timer => find_timer(tmr, name, stack)

    If (.not. Present(stack)) Call push_stack(tmr, name)
    Call mtime(current_timer%time%start)
    current_timer%time%running = .true.

  End Subroutine start_timer

  Subroutine stop_timer(tmr, name, stack)
    !!------------------------------------------------!
    !!
    !! Stop a timer running on a given timer system
    !! If stack is supplied bypass standard call stack
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Character(Len=*), Intent(In   ) :: name
    Type(call_stack), Optional      :: stack

    Type(node), Pointer :: current_timer

    current_timer => find_timer(tmr, name, stack)

    If (.not. Present(stack)) Call pop_stack(tmr, name)
    If (.not. current_timer%time%running) &
      & Call timer_error(tmr, 'Timer '//Trim(current_timer%time%name)//' stopped but not running')

    Call mtime(current_timer%time%stop)

    current_timer%time%running = .false.
    current_timer%time%last = current_timer%time%stop - current_timer%time%start
    If (current_timer%time%last > current_timer%time%max) current_timer%time%max = current_timer%time%last
    If (current_timer%time%last < current_timer%time%min) current_timer%time%min = current_timer%time%last
    current_timer%time%total = current_timer%time%total + current_timer%time%last
    current_timer%time%calls = current_timer%time%calls + 1

  End Subroutine stop_timer

  Subroutine start_timer_path(tmr, name_in, start_parents)
    !!------------------------------------------------!
    !!
    !! Start a timer running on a given timer system
    !! ignoring the timer call stack
    !! If start_parents start all non-running timers
    !! which are on the path -- Default TRUE
    !! - Path should be colon separated
    !!
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type),  Intent(InOut) :: tmr
    Character(Len=*),  Intent(In   ) :: name_in
    Logical, Optional, Intent(In   ) :: start_parents

    Character(Len=max_name) :: name
    Integer                 :: depth
    Logical                 :: parents
    Type(call_stack)        :: stack
    Type(node), Pointer     :: current_timer, is_running

    Call timer_split_stack_string(tmr, name_in, stack, name)
    current_timer => find_timer(tmr, name, stack)

    parents = .true.
    If (Present(start_parents)) parents = start_parents

    If (parents) Then
      Do depth = 1, tmr%stack%depth
        tmr%stack%depth = depth - 1
        is_running => find_timer(tmr, stack%name(depth), stack)
        If (.not. is_running%time%running) Call start_timer(tmr, (tmr%stack%name(depth)), stack)
      End Do
      tmr%stack%depth = depth - 1
    End If

    Call start_timer(tmr, name, stack)

  End Subroutine start_timer_path

  Subroutine stop_timer_path(tmr, name_in, stop_parents)
    !!------------------------------------------------!
    !!
    !! Stop a timer running on a given timer system
    !! ignoring the timer call stack
    !! If stop_parents stop all running timers
    !! which are on the path -- Default TRUE
    !! - Path should be colon separated
    !!
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type),  Intent(InOut) :: tmr
    Character(Len=*),  Intent(In   ) :: name_in
    Logical, Optional, Intent(In   ) :: stop_parents

    Character(Len=max_name) :: name
    Integer                 :: depth
    Logical                 :: parents
    Type(call_stack)        :: stack
    Type(node), Pointer     :: is_running

    Call timer_split_stack_string(tmr, name_in, stack, name)

    Call stop_timer(tmr, name, stack)

    parents = .true.
    If (Present(stop_parents)) parents = stop_parents

    If (parents) Then
      Do depth = tmr%stack%depth, 1, -1
        tmr%stack%depth = depth - 1
        is_running => find_timer(tmr, tmr%stack%name(depth), stack)
        If (is_running%time%running) Call stop_timer(tmr, (tmr%stack%name(depth)), stack)
      End Do
    End If

  End Subroutine stop_timer_path

  Subroutine timer_report(tmr, comm)
    !!------------------------------------------------!
    !!
    !! Stop the main timer system and print the full
    !! table of timers to stdout
    !!
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Type(comms_type), Intent(InOut) :: comm

    Character(Len=256), Allocatable, Dimension(:) :: message
    Integer                                       :: ierr, proc
    Type(node), Pointer                           :: current_timer

    Call stop_timer(tmr, 'Main')

    current_timer => tmr%tree%head%child

    Allocate (message(-2:tmr%tree%n_timers + 3), stat=ierr)
    If (ierr > 0) Call timer_error(tmr, 'Error allocating message in timer_print_tree')

    Call timer_print_tree(comm, tmr, current_timer, tmr%max_depth, -1, message)
    Call timer_write(message, tmr)

    If (tmr%proc_detail) Then
      Do proc = 0, comm%mxnode - 1
        If (comm%idnode == proc) Then
          Call timer_print_tree(comm, tmr, current_timer, tmr%max_depth, proc, message)
        End If

        If (proc /= 0) Then
          If (comm%idnode == proc) Call gsend(comm, message, 0, timer_tag)
          If (comm%idnode == 0) Call grecv(comm, message, proc, timer_tag)
        End If

        Call timer_write(message, tmr)
        Call gsync(comm)
      End Do
    End If

  End Subroutine timer_report

  Subroutine timer_print_tree(comm, tmr, init_node, max_depth, proc_id, message)
    !!------------------------------------------------!
    !!
    !! Return a table of the given timer system to
    !! the message variable
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(comms_type),                   Intent(InOut) :: comm
    Type(timer_type),                   Intent(In   ) :: tmr
    Type(node), Target,                 Intent(In   ) :: init_node
    Integer,                            Intent(In   ) :: max_depth, proc_id
    Character(Len=*), Dimension(-2:), Intent(  Out) :: message

    Character(Len=7)    :: depth_symb
    Character(Len=8)    :: proc_string
    Integer             :: depth, itimer, write_node
    Real(Kind=wp)       :: call_av, call_max, call_min, sum_timed, total_av, total_elapsed, &
                           total_max, total_min
    Type(node), Pointer :: current_timer
    Character(Len = 256 ) :: fline,fhead,fcontent

    fline = '(1X, "+", 28("-"), 2("+", 10("-")), 7("+", 12("-")), "+")'
    fhead = '(1X, "|", 12X, "Name", 12X, "| Process  ", "|  Calls   ", "|  Call Min  ", "|  Call Max  ",'//&
          & '"|  Call Ave  ", "|  Tot Min   ", "|   Tot Max  ", "|   Tot Ave  ", "|      %     ", "|")'
    fcontent = '(1X, "|", 1X, A7, 1X, A18, 1X, "|", 1X, A8, 1X, "|", 1X, I8, 1X, "|", 1X, ES10.3, 1X, "|", 1X, ES10.3, 1X,'//&
          & '"|", 1X, ES10.3, 1X, "|", 1X, ES10.3, 1X, "|", 1X, ES10.3, 1X, "|", 1X, ES10.3, 1X, "|", 1X, ES10.3, 1X, "|")'

    message(:) = ''

    sum_timed = 0.0_wp

    If (proc_id < 0) Then
      write_node = 0
      proc_string = "All"
    Else
      write_node = proc_id
      Write (proc_string, '(i8.1)') proc_id
    End If

    current_timer => init_node
    total_elapsed = current_timer%time%total

    depth = 0
    itimer = 0

    ! Write table open and header
    Write (message(-2), Trim(fline))
    Write (message(-1), Trim(fhead))

    Do While (depth > -1)

      If (current_timer%time%running) &
        & Call timer_error(tmr, 'Program terminated while timer '//Trim(current_timer%time%name)//' still running')

      If (depth == 0 .and. current_timer%time%name /= "Main") Then
        depth_symb = Repeat('-', 7)
      Else If (Associated(current_timer%child)) Then
        depth_symb = Repeat(" ", depth)//"|v"
      Else
        depth_symb = Repeat(" ", depth)//"|-"
      End If

      total_min = current_timer%time%total
      total_max = current_timer%time%total
      total_av = current_timer%time%total

      If (proc_id < 0) Then
        Call gmin(comm, total_min)
        Call gmax(comm, total_max)
        Call gsum(comm, total_av)
        total_av = total_av / comm%mxnode
      End If

      If (depth == 1 .and. current_timer%parent%time%name == "Main") sum_timed = sum_timed + total_av

      call_min = current_timer%time%min
      call_max = current_timer%time%max

      If (proc_id < 0) Then
        Call gmin(comm, call_min)
        Call gmax(comm, call_max)
      End If

      call_av = total_av / current_timer%time%calls

      Write (message(itimer), Trim(fcontent)) depth_symb, current_timer%time%name, proc_string, current_timer%time%calls, &
        & call_min, call_max, call_av,  &
        & total_min, total_max, total_av, total_av * 100.0_wp / total_elapsed

      If (Associated(current_timer%child) .and. depth < max_depth) Then
        current_timer => current_timer%child
        depth = depth + 1
      Else If (Associated(current_timer%next_sibling)) Then
        current_timer => current_timer%next_sibling
      Else If (Associated(current_timer%parent)) Then ! Recurse back up
        Do While (Associated(current_timer%parent))
          current_timer => current_timer%parent
          depth = depth - 1
          If (Associated(current_timer%next_sibling)) Then
            current_timer => current_timer%next_sibling
            Exit
          End If
        End Do
      Else
        Exit
      End If
      itimer = itimer + 1

    End Do

    Write (message(itimer), Trim(fcontent)) Repeat('-', 7), "Untimed           ", proc_string, 0, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & total_elapsed - sum_timed, 100.0_wp - sum_timed * 100.0_wp / total_elapsed
    Write (message(itimer + 1), Trim(fline))
    Write (message(itimer + 2), '(a)') ''

    !100 Format(1X, "+", 28("-"), 2("+", 10("-")), 7("+", 16("-")), "+")
    !101 Format(1X, "|", 12X, "Name", 12X, "| Process  ", "|  Calls   ", "|   Call Min     ", "|   Call Max     ",&
    !      & "|   Call Ave     ", "|    Tot Min     ", "|    Tot Max     ", "|    Tot Ave     ", "|       %        ", "|")
    !102 Format(1X, "|", 1X, A7, 1X, A18, 1X, "|", 1X, A8, 1X, "|", 1X, I8, 1X, "|", 1X, ES14.7, 1X, "|", 1X, ES14.7, 1X,&
          !& "|", 1X, ES14.7, 1X, "|", 1X, ES14.7, 1X, "|", 1X, ES14.7, 1X, "|", 1X, ES14.7, 1X, "|", 1X, ES14.7, 1X, "|")

  End Subroutine timer_print_tree

  Subroutine init_timer(current_timer, name)
    !!------------------------------------------------!
    !!
    !! Initialise a node timer
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(node_timer) :: current_timer
    Character(Len=*) :: name

    current_timer%name = name
    current_timer%calls = 0
    current_timer%max = -1.0_wp
    current_timer%min = Huge(1.0_wp)
    current_timer%total = 0.0_wp
    current_timer%last = Huge(1.0_wp)

  End Subroutine init_timer

  Subroutine timer_last_time(tmr, name, screen)
    !!------------------------------------------------!
    !!
    !! Write the length of the previous call to a given
    !! node timer
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr
    Character(Len=*)                :: name
    Logical, Optional               :: screen

    Character(Len=72)   :: message
    Logical             :: to_screen
    Type(node), Pointer :: current_timer

    to_screen = .false.
    If (Present(screen)) to_screen = screen

    current_timer => find_timer(tmr, name)
    If (to_screen) Then
      Write (eu, *) current_timer%time%name, current_timer%time%calls, current_timer%time%last
    Else
      Write (message, '(a,2(1X,i0))') current_timer%time%name, current_timer%time%calls, current_timer%time%last
      Call timer_write(message, tmr)
    End If
  End Subroutine timer_last_time

  Subroutine time_elapsed(tmr)
    !!------------------------------------------------!
    !!
    !! Write the elapsed time since the start of the
    !! program
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(InOut) :: tmr

    Character(Len=256) :: message

    Call gtime(tmr%elapsed)
    Write (message, '(a,f12.3,a)') "time elapsed since job start: ", tmr%elapsed, " sec"
    Call timer_write(message, tmr)

  End Subroutine time_elapsed

  Subroutine timer_write_mul(message, timer_in)
    !!------------------------------------------------!
    !!
    !! Write multiple lines to standard out
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character(Len=*), Dimension(:), Intent(In   ) :: message
    Type(timer_type), Optional,     Intent(In   ) :: timer_in

    Integer :: i

    If (Present(timer_in)) Then

      If (timer_in%proc_id == 0) Then
        Do i = 1, Size(message)
          Write (timer_in%out_unit, '(a)') Trim(message(i))
        End Do
      End If

    Else

      Do i = 1, Size(message)
        Write (eu, '(a)') Trim(message(i))
      End Do

    End If

  End Subroutine timer_write_mul

  Subroutine timer_write_sing(message, timer_in)
    !!------------------------------------------------!
    !!
    !! Write a single line to standard out
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Character(Len=*),           Intent(In   ) :: message
    Type(timer_type), Optional, Intent(In   ) :: timer_in

    If (Present(timer_in)) Then

      If (timer_in%proc_id == 0) Then
        Write (timer_in%out_unit, '(a)') Trim(message)
      End If

    Else

      Write (eu, '(a)') Trim(message)

    End If

  End Subroutine timer_write_sing

  Subroutine timer_error(tmr, message)
    !!------------------------------------------------!
    !!
    !! Report an error with the timer system
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type), Intent(In   ) :: tmr
    Character(len=*), Intent(In   ) :: message

    Call timer_write('')
    ! Call timer_report(dummy_timer, timer_comm)
    Call timer_write(message)
    Call timer_write('')
    Call dump_call_stack(tmr%stack)

    Stop
  End Subroutine timer_error

  Subroutine timer_split_stack_string(tmr, stack_string, newStack, name)
    !!------------------------------------------------!
    !!
    !! Split a colon-separated path string into a stack
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!------------------------------------------------!
    Type(timer_type),        Intent(In   ) :: tmr
    Character(Len=*),        Intent(In   ) :: stack_string
    Type(call_stack),        Intent(  Out) :: newStack
    Character(Len=max_name), Intent(  Out) :: name

    Character(Len=256) :: stack
    Integer            :: cnt, i

    stack = Adjustl(stack_string)
    cnt = 1
    Do i = 1, Len(stack)
      If (stack(i:i) == ":") cnt = cnt + 1
    End Do

    If (cnt > max_depth) Call timer_error(tmr, 'Stack depth greater than max depth in timer_split_stack')

    newStack%depth = 0
    Do While (Index(stack, ':') > 0)
      newStack%depth = newStack%depth + 1
      newStack%name(newStack%depth) = Trim(stack(1:Index(stack, ':') - 1))

      stack(1:Index(stack, ':')) = " "
      stack = Adjustl(stack)

    End Do

    name = Trim(stack)

  End Subroutine timer_split_stack_string

End Module timer
