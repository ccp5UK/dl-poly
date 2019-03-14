Module timer
  !! This module has no header !
  Use kinds, Only : wp
  Use comms, Only : comms_type,gtime,mtime,gmin,gmax,gsum,gsync,gsend,grecv,timer_tag,abort_comms
  Implicit None

  Private

  Integer, Parameter :: max_depth = 6, max_name = 18


  Type :: node_timer
    Character( Len = max_name ) :: name
    Integer           :: id
    Real( Kind = wp ) :: max, min, total, last
    Real( Kind = wp ) :: start, stop
    Integer           :: calls
    Logical           :: running = .false.
  end type node_timer

  Type :: call_stack
    Character ( Len = max_name ), dimension( max_depth ) :: name
    Integer :: depth = 0
  end type call_stack

  Type :: node
    Type ( node_timer ) :: time
    Type ( timer_tree ), pointer :: tree => null()
    Type ( node ), pointer :: child => null()
    Type ( node ), pointer :: parent => null()
    Type ( node ), pointer :: next_sibling => null()
  end type node
  
  Type :: timer_tree
    Type ( node ), pointer :: head => null()
    Integer :: n_timers = 0
  end type timer_tree

  Type, Public :: timer_type
    Type ( timer_tree ), pointer :: tree
    Type ( call_stack ) :: stack
    Real( Kind = wp) :: elapsed,job,clear_screen
    Logical :: proc_detail = .false.
    Integer :: max_depth = 1
    Integer :: proc_id
    Integer :: out_unit
  End Type timer_type
  
  interface timer_write
    module procedure timer_write_sing
    module procedure timer_write_mul
  end interface timer_write

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

  Subroutine dump_call_stack ( stack )
    Type ( call_stack ) :: stack
    Integer :: i

    Call timer_write('')
    Call timer_write('Process stack:')
    Do i = 1, stack%depth
      call timer_write(stack%name(i))
    End Do
    
  End Subroutine dump_call_stack

  Subroutine push_stack ( tmr, name )
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent ( In    ) :: name

    tmr%stack%depth = tmr%stack%depth + 1
    If ( tmr%stack%depth > max_depth ) &
      & Call timer_error(tmr, 'Call stack exceeds max depth : recursive call or unended timer?')
    tmr%stack%name(tmr%stack%depth) = name

  End Subroutine push_stack

  Subroutine pop_stack ( tmr, name )
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    ) :: name

    If ( name /= tmr%stack%name(tmr%stack%depth)) &
      & Call timer_error(tmr, 'Child timer '//name//' ended before parent')
    tmr%stack%name(tmr%stack%depth) = ''
    tmr%stack%depth = tmr%stack%depth - 1

  End Subroutine pop_stack

  Subroutine init_timer_system ( tmr, nrite, comm )
    Type ( timer_type ), intent ( inout ) :: tmr
    Type ( comms_type ), intent ( in    ) :: comm
    Integer, intent ( in ) :: nrite

    tmr%proc_id = comm%idnode
    tmr%out_unit = nrite
    Allocate(tmr%tree)
    Allocate(tmr%tree%head)
    tmr%tree%head%tree => tmr%tree
    Call init_timer ( tmr%tree%head%time, 'Head')
    Call start_timer( tmr, 'Main' )
    
  End Subroutine init_timer_system

  Function find_timer ( tmr, name, stack_in ) result(current)
    Type ( timer_type ), Intent ( InOut ) :: tmr
    Character ( Len = * ) :: name
    Type ( call_stack ), Optional :: stack_in
    Type ( call_stack ) :: stack
    Type ( node ), Pointer :: current
    Integer :: depth

    If (Present(stack_in)) Then
      stack = stack_in
    Else
      stack = tmr%stack
    End If
    
    current => tmr%tree%head
    depth = 0

    Do While ( depth < stack%depth )
      If ( .not. associated(current%child)) call timer_error(tmr, 'Call stack does not match call tree (no child)')
      depth = depth + 1
      current => current%child
      Do While ( current%time%name /= stack%name(depth) )
        If ( .not. associated(current%next_sibling)) &
          & Call timer_error(tmr, 'Call stack does not match call tree (no sibling)')
        current => current%next_sibling
      End Do
    End Do

    If (current%time%name == name ) Then
      Continue
    Else If (.not. associated(current%child)) Then
      Call init_child_node(name, current)
      current => current%child
      Return
    Else
      current => current%child
      do while ( current%time%name /= name )
        If ( .not. Associated(current%next_sibling)) Then
          Call init_sibling_node(name, current)
          current => current%next_sibling
          Return
        End If
        current => current%next_sibling
      End Do
    End If

  End Function find_timer

  Subroutine init_child_node(name, parent)
    Character ( len = * ), Intent( In    ) :: name
    Type ( node ), Target :: parent

    Allocate ( parent%child )
    parent%child%parent => parent
    parent%child%tree => parent%tree
    parent%tree%n_timers = parent%tree%n_timers + 1
    Call init_timer(parent%child%time, name)

  End Subroutine init_child_node

  Subroutine init_sibling_node(name, sibling)
    Character ( len = * ) :: name
    Type ( node ) :: sibling
    Type ( node ), pointer :: child

    Allocate ( sibling%next_sibling )
    child => sibling%next_sibling
    child%parent => sibling%parent

    sibling%next_sibling%tree => sibling%tree
    sibling%tree%n_timers = sibling%tree%n_timers + 1
    
    Call init_timer(child%time, name)

  End Subroutine init_sibling_node

  Subroutine start_timer(tmr, name, stack)
    !! This routine has no header !
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name
    Type ( call_stack ), optional :: stack
    Type ( node ), Pointer :: current_timer

    current_timer => find_timer(tmr, name, stack)

    if (.not. present(stack)) Call push_stack(tmr, name)
    Call mtime(current_timer%time%start)
    current_timer%time%running = .true.

  End Subroutine start_timer

  Subroutine stop_timer(tmr, name, stack)
    !! This routine has no header !
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name
    Type ( call_stack ), optional :: stack
    Type ( node ), Pointer :: current_timer

    current_timer => find_timer(tmr, name, stack)

    if (.not. present(stack)) Call pop_stack(tmr, name)
    If ( .not. current_timer%time%running ) &
      & Call timer_error(tmr, 'Timer '//Trim(current_timer%time%name)//' stopped but not running')

    Call mtime(current_timer%time%stop)

    current_timer%time%running = .false.
    current_timer%time%last = current_timer%time%stop - current_timer%time%start
    If ( current_timer%time%last > current_timer%time%max ) current_timer%time%max = current_timer%time%last
    If ( current_timer%time%last < current_timer%time%min ) current_timer%time%min = current_timer%time%last
    current_timer%time%total = current_timer%time%total + current_timer%time%last
    current_timer%time%calls = current_timer%time%calls + 1

  End Subroutine stop_timer

  Subroutine start_timer_path(tmr, name_in, start_parents)
    !! This routine has no header !
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name_in
    Logical, Intent ( In    ), Optional :: start_parents
    Logical :: parents
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: is_running
    Type ( node ), Pointer :: current_timer
    Integer :: depth

    Call timer_split_stack_string(tmr, name_in, stack, name)
    current_timer => find_timer(tmr, name, stack)
    
    parents = .true.
    if ( present(start_parents) ) parents = start_parents

    if (parents) then
      do depth = 1, tmr%stack%depth
        tmr%stack%depth = depth-1
        is_running => find_timer(tmr, stack%name(depth), stack )
        if (.not. is_running%time%running) call start_timer(tmr, (tmr%stack%name(depth)), stack)
      end do
      tmr%stack%depth = depth-1
    end if
    
    call start_timer(tmr, name, stack)

  End Subroutine start_timer_path

  Subroutine stop_timer_path(tmr, name_in, stop_parents)
    !! This routine has no header !
    Type ( timer_type ), Intent( InOut ) :: tmr
    Character ( Len = * ), Intent( In    )  :: name_in
        Logical, Intent ( In    ), Optional :: stop_parents
    Logical :: parents
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: is_running
    Type ( node ), Pointer :: current_timer
    Integer :: depth

    Call timer_split_stack_string(tmr, name_in, stack, name)

    call stop_timer(tmr, name, stack)

    parents = .true.
    if ( present(stop_parents) ) parents = stop_parents

    if (parents) then
      do depth = tmr%stack%depth, 1, -1
        tmr%stack%depth = depth-1
        is_running => find_timer(tmr, tmr%stack%name(depth), stack )
        if (is_running%time%running) call stop_timer(tmr, (tmr%stack%name(depth)), stack)
      end do
    end if

  End Subroutine stop_timer_path
  
  Subroutine timer_report(tmr,comm)
    !! This routine has no header !
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm
    Type ( node ), Pointer :: current_timer
    Character( Len = 138 ), Dimension(:), Allocatable :: message
    Integer :: proc
    Integer :: ierr
    
    Call stop_timer(tmr, 'Main')

    current_timer => tmr%tree%head%child
    
    Allocate(message(-2:tmr%tree%n_timers+3), stat=ierr)
    If ( ierr > 0 ) call timer_error(tmr, 'Error allocating message in timer_print_tree')

    Call timer_print_tree(comm, tmr, current_timer, tmr%max_depth, -1, message)
    Call timer_write(message, tmr)
    
    If (tmr%proc_detail) Then
      Do proc = 0, comm%mxnode-1
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
    Implicit None
    Type ( comms_type ), Intent ( InOut ) :: comm
    Type ( timer_type ), Intent ( In    ) :: tmr
    Type ( node ), Target, Intent ( In     )  :: init_node
    Character( Len = 138 ), Dimension(-2:), Intent(   Out ) :: message
    Type ( node ), Pointer :: current_timer
    Integer, Intent ( In    ) :: proc_id
    Integer, Intent ( In    ) :: max_depth
    Integer :: write_node
    Real ( Kind = wp ) :: total_min, total_max, total_av
    Real ( Kind = wp ) :: call_min, call_max, call_av
    Real ( Kind = wp ) :: sum_timed, total_elapsed
    Integer :: depth, itimer
    Character( Len = 8   ) :: proc_string
    Character( Len = 7   ) :: depth_symb

    message(:) = ''
    
    sum_timed = 0.0_wp

    If ( proc_id < 0 ) Then
      write_node = 0
      proc_string = "All"
    Else
      write_node = proc_id
      Write(proc_string,'(i8.1)') proc_id
    End If
    
    current_timer => init_node
    total_elapsed = current_timer%time%total
    
    depth = 0
    itimer = 0
    
    ! Write table open and header
    Write(message(-2), 100)
    Write(message(-1), 101)

    Do While (depth > -1)

      If (current_timer%time%running) &
        & Call timer_error(tmr, 'Program terminated while timer '// trim(current_timer%time%name)//' still running')
      
      If ( depth == 0 .and. current_timer%time%name /= "Main") Then
        depth_symb = Repeat('-', 7)
      Else If ( Associated(current_timer%child) ) Then
        depth_symb = Repeat(" ",depth)//"|v"
      Else
        depth_symb = Repeat(" ",depth)//"|-"
      End If

      total_min = current_timer%time%total
      total_max = current_timer%time%total
      total_av  = current_timer%time%total

      If ( proc_id < 0 ) Then
        Call gmin(comm,total_min)
        Call gmax(comm,total_max)
        Call gsum(comm,total_av)
        total_av = total_av / comm%mxnode
      End If
      
      if (depth == 1 .and. current_timer%parent%time%name == "Main") sum_timed = sum_timed + total_av

      call_min  = current_timer%time%min
      call_max  = current_timer%time%max

      If ( proc_id < 0 ) Then
        Call gmin(comm,call_min)
        Call gmax(comm,call_max)
      End If
      
      call_av   = total_av/current_timer%time%calls

      Write(message(itimer), 102 ) depth_symb,current_timer%time%name, proc_string, current_timer%time%calls, &
        & call_min, call_max, call_av,  &
        & total_min, total_max, total_av, total_av*100.0_wp/total_elapsed
      
      If (Associated(current_timer%child) .and. depth < max_depth ) Then
        current_timer => current_timer%child
        depth = depth + 1
      Else If (Associated(current_timer%next_sibling)) Then
        current_timer => current_timer%next_sibling
      Else If (Associated(current_timer%parent)) Then ! Recurse back up
        Do While (Associated(current_timer%parent))
          current_timer => current_timer%parent
          depth = depth - 1
          If (associated(current_timer%next_sibling)) Then
            current_timer => current_timer%next_sibling
            Exit
          End If
        End Do
      Else
        Exit
      End If
      itimer = itimer + 1

    End Do

    Write(message(itimer), 102) Repeat('-',7), "Untimed           ", proc_string, 0, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & 
      & total_elapsed - sum_timed , 100.0_wp - sum_timed*100.0_wp/total_elapsed
    Write(message(itimer+1), 100)
    Write(message(itimer+2), *) ''

100 Format(1X,"+",28("-"),2("+",10("-")),7("+",11("-")),"+")
101 Format(1X,"|",12X,"Name",12X,"| Process  ","|  Calls   ","| Call Min  ","| Call Max  ",&
      & "| Call Ave  ","|  Tot Min  ","|  Tot Max  ","|  Tot Ave  ","|     %     ","|")
102 Format(1X,"|",1X,A7,1X,A18,1X,"|",1X,A8,1X,"|",1X,I8,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,&
      & "|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",2X,F8.4,1X,"|")
    
  End Subroutine timer_print_tree
  
  Subroutine init_timer(current_timer, name)
    !! This routine has no header !
    Type ( node_timer ) :: current_timer
    Character ( Len = * )  :: name

    current_timer%name  = name
    current_timer%calls = 0
    current_timer%max   = -1.0_wp
    current_timer%min   = huge(1.0_wp)
    current_timer%total = 0.0_wp
    current_timer%last  = huge(1.0_wp)

  End Subroutine init_timer

  Subroutine timer_last_time(tmr, name, screen)
    !! This routine has no header !
    Type ( timer_type ), Intent ( InOut ) :: tmr
    Character ( Len = * )  :: name
    Logical, Optional :: screen
    Logical :: to_screen
    Character ( Len = 72 ) :: message
    Type ( node ), pointer :: current_timer

    to_screen = .false.
    If (Present(screen)) to_screen = screen

    current_timer => find_timer(tmr, name)
    If (to_screen) Then
      Write(0,*) current_timer%time%name, current_timer%time%calls, current_timer%time%last
    Else
      Write(message,'(a,2(1X,i0))') current_timer%time%name, current_timer%time%calls, current_timer%time%last
      Call timer_write(message, tmr)
    End If
  End Subroutine timer_last_time

  Subroutine time_elapsed(tmr)
    !! This routine has no header !
    Type ( timer_type ), Intent( InOut ) :: tmr

    Character( Len = 256 ) :: message

    Call gtime(tmr%elapsed)
    Write(message,'(a,f12.3,a)') "time elapsed since job start: ", tmr%elapsed, " sec"
    Call timer_write(message, tmr)
    
  End Subroutine time_elapsed

  Subroutine timer_write_mul(message, timer_in)
    Character( Len = * ), Dimension(:), Intent( In    ) :: message
    Type ( timer_type ), Optional, Intent ( In    ):: timer_in
    Integer :: i

    If (Present(timer_in)) Then

      If ( timer_in%proc_id == 0 ) Then
        Do i = 1, Size(message)
          Write(timer_in%out_unit,*) message(i)
        End Do
      End If
      
    Else

      Do i = 1, Size(message)
        Write(0,*) message(i)
      End Do
      
    End If
    
  End Subroutine timer_write_mul

  Subroutine timer_write_sing(message, timer_in)
    Character( Len = * ), Intent( In    ) :: message
    Type ( timer_type ), Optional, Intent ( In    ):: timer_in
    
    If (Present(timer_in)) Then

      If ( timer_in%proc_id == 0 ) Then
          Write(timer_in%out_unit,*) message
      End If
      
    Else

        Write(0,*) message
      
    End If
    
  End Subroutine timer_write_sing

  Subroutine timer_error(tmr, message)
    Type ( timer_type ), Intent ( In    ):: tmr
    Character ( len = * ), Intent( In    ) :: message
   
    Call timer_write('')
    ! Call timer_report(dummy_timer, timer_comm)
    Call timer_write(message)
    Call timer_write('')
    Call dump_call_stack(tmr%stack)

    Stop
  End Subroutine timer_error

  Subroutine timer_split_stack_string(tmr, stack_string, newStack, name)
    Type ( timer_type ), Intent ( In    ):: tmr
    Character( Len = * ), Intent( In    ) :: stack_string
    Character( Len = 256 ) :: stack
    Type ( call_stack ), Intent(   Out ) :: newStack
    Character( Len = max_name ), Intent(   Out ) :: name
    Integer :: i
    Integer :: cnt
    
    stack = Adjustl(stack_string)
    cnt = 1
    Do i = 1, Len(stack)
      If (stack(i:i) == ":") cnt = cnt + 1
    End Do
   
    If (cnt > max_depth) Call timer_error(tmr, 'Stack depth greater than max depth in timer_split_stack')

    newStack%depth = 0
    Do While (Index(stack,':') > 0)
      newStack%depth = newStack%depth + 1
      newStack%name(newStack%depth) = Trim(stack(1:index(stack,':')-1))

      stack(1:index(stack,':')) = " "
      stack = adjustl(stack)
      
    End Do

    name = Trim(stack)
          
  End Subroutine timer_split_stack_string

End Module timer
