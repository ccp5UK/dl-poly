Module timer
  !! This module has no header !
  Use kinds, Only : wp
  Use comms, Only : comms_type,gtime,mtime,gmin,gmax,gsum,gsync,gsend,grecv,timer_tag,abort_comms
  Implicit None

  Private

  Integer, Parameter :: max_depth = 6, max_name = 18


  Type :: timer_type_new
    Character( Len = max_name ) :: name
    Integer           :: id
    Real( Kind = wp ) :: max, min, total, last
    Real( Kind = wp ) :: start, stop
    Integer           :: calls
    Logical           :: running = .false.
  end type timer_type_new

  Type :: call_stack
    Character ( Len = max_name ), dimension( max_depth ) :: name
    Integer :: depth = 0
  end type call_stack

  Type :: node
    Type ( timer_type_new ) :: time
    Type ( timer_tree ), pointer :: tree => null()
    Type ( node ), pointer :: child => null()
    Type ( node ), pointer :: parent => null()
    Type ( node ), pointer :: next_sibling => null()
  end type node
  
  Type :: timer_tree
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: head => null()
    Integer :: n_timers = 0
  end type timer_tree

  Type, Public :: timer_type
    Type ( timer_tree ) :: tree
    Real( Kind = wp) :: elapsed,job,clear_screen
    Logical :: proc_detail = .false.
    Integer :: max_depth = 1
    Integer :: proc_id
    Integer :: out_unit
  End Type timer_type
  
  Type( call_stack ),                              save :: calls
  Type( timer_tree ), target,                      save :: call_tree

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
    type ( call_stack ), optional :: stack
    integer :: i

    call timer_write('')
    call timer_write('Process stack:')
    if ( .not. present(stack) ) then
      do i = 1, calls%depth
        call timer_write(calls%name(i))
      end do
    else
      do i = 1, stack%depth
        call timer_write(stack%name(i))
      end do
    end if
    
  end Subroutine dump_call_stack

  Subroutine push_stack ( name )
    Character ( Len = * ) :: name

    calls%depth = calls%depth + 1
    if ( calls%depth > max_depth ) &
      & call timer_error('Call stack exceeds max depth : recursive call or unended timer?')
    calls%name(calls%depth) = name

  end Subroutine push_stack

  Subroutine pop_stack ( name )
    Character ( Len = * ) :: name

    If ( name /= calls%name(calls%depth)) &
    & call timer_error('Child timer '//name//' ended before parent')
    calls%name(calls%depth) = ''
    calls%depth = calls%depth - 1

  End Subroutine pop_stack

  subroutine init_timer_system ( timer_in, nrite, comm )
    Type ( timer_type ), intent ( inout ) :: timer_in
    Type ( comms_type ), intent ( in    ) :: comm
    Integer, intent ( in ) :: nrite

    timer_in%proc_id = comm%idnode
    timer_in%out_unit = nrite
    allocate(call_tree%head)
    call_tree%head%tree => call_tree
    call init_timer ( call_tree%head%time, 'Head')
    Call start_timer('Main')
    
  end subroutine init_timer_system

  Function find_timer ( name, stack_in ) result(current)
    Character ( Len = * ) :: name
    Type ( call_stack ), optional :: stack_in
    Type ( call_stack ) :: stack
    Type ( node ), Pointer :: current
    Integer :: depth

    if (present(stack_in)) then
      stack = stack_in
    else
      stack = calls
    end if
    
    current => call_tree%head
    depth = 0

    do while ( depth < stack%depth )
      if ( .not. associated(current%child)) call timer_error('Call stack does not match call tree (no child)')
      depth = depth + 1
      current => current%child
      do while ( current%time%name /= stack%name(depth) )
        if ( .not. associated(current%next_sibling)) &
          & call timer_error('Call stack does not match call tree (no sibling)')
        current => current%next_sibling
      end do
    end do

    if (current%time%name == name ) then
      continue
    else if (.not. associated(current%child)) then
      call init_child_node(name, current)
      current => current%child
      return
    else
      current => current%child
      do while ( current%time%name /= name )
        if ( .not. associated(current%next_sibling)) then
          call init_sibling_node(name, current)
          current => current%next_sibling
          return
        end if
        current => current%next_sibling
      end do
    end if

  end Function find_timer

  subroutine init_child_node(name, parent)
    Character ( len = * ) :: name
    Type ( node ), target :: parent

    allocate ( parent%child )
    parent%child%parent => parent
    parent%child%tree => parent%tree
    parent%tree%n_timers = parent%tree%n_timers + 1
    call init_timer(parent%child%time, name)

  end subroutine init_child_node

  subroutine init_sibling_node(name, sibling)
    Character ( len = * ) :: name
    Type ( node ) :: sibling
    Type ( node ), pointer :: child

    allocate ( sibling%next_sibling )
    child => sibling%next_sibling
    child%parent => sibling%parent

    sibling%next_sibling%tree => sibling%tree
    sibling%tree%n_timers = sibling%tree%n_timers + 1
    
    call init_timer(child%time, name)

  end subroutine init_sibling_node

  Subroutine start_timer(name)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Type ( node ), pointer :: current_timer

    current_timer => find_timer(name)

    call push_stack(name)
    Call mtime(current_timer%time%start)
    current_timer%time%running = .true.

  end Subroutine start_timer

  Subroutine stop_timer(name)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Type ( node ), pointer :: current_timer

    current_timer => find_timer(name)

    call pop_stack(name)
    if ( .not. current_timer%time%running ) call timer_error('Timer '//trim(current_timer%time%name)//' stopped but not running')

    Call mtime(current_timer%time%stop)

    current_timer%time%running = .false.
    current_timer%time%last = current_timer%time%stop - current_timer%time%start
    if ( current_timer%time%last > current_timer%time%max ) current_timer%time%max = current_timer%time%last
    if ( current_timer%time%last < current_timer%time%min ) current_timer%time%min = current_timer%time%last
    current_timer%time%total = current_timer%time%total + current_timer%time%last
    current_timer%time%calls = current_timer%time%calls + 1

  End Subroutine stop_timer

  Subroutine start_timer_path(name_in)
    !! This routine has no header !
    Character ( Len = * )  :: name_in
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: current_timer

    call timer_split_stack_string(name_in, stack, name)
    current_timer => find_timer(name, stack)
    
    Call mtime(current_timer%time%start)
    current_timer%time%running = .true.

  end Subroutine start_timer_path

  Subroutine stop_timer_path(name_in)
    !! This routine has no header !
    Character ( Len = * )  :: name_in
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: current_timer

    call timer_split_stack_string(name_in, stack, name)
    current_timer => find_timer(name, stack)

    if ( .not. current_timer%time%running ) call timer_error('Timer '//trim(current_timer%time%name)//' stopped but not running')

    Call mtime(current_timer%time%stop)

    current_timer%time%running = .false.
    current_timer%time%last = current_timer%time%stop - current_timer%time%start
    if ( current_timer%time%last > current_timer%time%max ) current_timer%time%max = current_timer%time%last
    if ( current_timer%time%last < current_timer%time%min ) current_timer%time%min = current_timer%time%last
    current_timer%time%total = current_timer%time%total + current_timer%time%last
    current_timer%time%calls = current_timer%time%calls + 1

  End Subroutine stop_timer_path
  
  Subroutine timer_report(tmr,comm)
    !! This routine has no header !
    Type( timer_type ), intent( In    ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm
    Type ( node ), pointer :: current_timer
    Character( Len = 138 ), dimension(:), allocatable :: message
    Integer :: proc
    Integer :: ierr
    
    call stop_timer('Main')

    current_timer => call_tree%head%child
    allocate(message(-2:call_tree%n_timers+3), stat=ierr)
    if ( ierr > 0 ) call timer_error('Error allocating message in timer_print_tree')

    call timer_print_tree(comm, tmr, current_timer, tmr%max_depth, -1, message)
    call timer_write(message, tmr)
    
    if (tmr%proc_detail) then
      do proc = 0, comm%mxnode-1
        if (comm%idnode == proc) then
          call timer_print_tree(comm, tmr, current_timer, tmr%max_depth, proc, message)
        end if

        if (proc /= 0) then
          if (comm%idnode == proc) call gsend(comm, message, 0, timer_tag)
          if (comm%idnode == 0) call grecv(comm, message, proc, timer_tag)
        end if
        
        call timer_write(message, tmr)
        call gsync(comm)
      end do
    end if

  End Subroutine timer_report

  Subroutine timer_print_tree(comm, tmr, init_node, max_depth, proc_id, message)
    implicit none
    Type ( comms_type ), intent ( inout ) :: comm
    Type ( timer_type ), intent ( in    ) :: tmr
    Type ( node ), target, intent ( in     )  :: init_node
    Character( Len = 138 ), dimension(-2:), intent(   out ) :: message
    Type ( node ), pointer :: current_timer
    Integer, intent ( in    ) :: proc_id
    Integer, intent ( in    ) :: max_depth
    Integer :: write_node
    Real ( kind = wp ) :: total_min, total_max, total_av
    Real ( kind = wp ) :: call_min, call_max, call_av
    Real ( kind = wp ) :: sum_timed, total_elapsed
    Integer :: depth, itimer
    Character( Len = 8   )                 :: proc_string
    Character( Len = 7   )                 :: depth_symb

    message(:) = ''
    
    sum_timed = 0.0_wp

    if ( proc_id < 0 ) then
      write_node = 0
      proc_string = "All"
    else
      write_node = proc_id
      write(proc_string,'(i8.1)') proc_id
    end if
    
    current_timer => init_node
    total_elapsed = current_timer%time%total
    
    depth = 0
    itimer = 0
    
    ! Write table open and header
    write(message(-2), 100)
    write(message(-1), 101)

    do while (depth > -1)

      if (current_timer%time%running) &
        & Call timer_error('Program terminated while timer '// trim(current_timer%time%name)//' still running')
      
      if ( depth == 0 .and. current_timer%time%name /= "Main") then
        depth_symb = repeat('-', 7)
      else if ( associated(current_timer%child) ) then
        depth_symb = repeat(" ",depth)//"|v"
      else
        depth_symb = repeat(" ",depth)//"|-"
      end if

      total_min = current_timer%time%total
      total_max = current_timer%time%total
      total_av  = current_timer%time%total

      if ( proc_id < 0 ) then
        Call gmin(comm,total_min)
        Call gmax(comm,total_max)
        Call gsum(comm,total_av)
        total_av = total_av / comm%mxnode
      end if
      
      if (depth == 1) sum_timed = sum_timed + total_av

      call_min  = current_timer%time%min
      call_max  = current_timer%time%max

      if ( proc_id < 0 ) then
        Call gmin(comm,call_min)
        Call gmax(comm,call_max)
      end if
      
      call_av   = total_av/current_timer%time%calls

      write(message(itimer), 102 ) depth_symb,current_timer%time%name, proc_string, current_timer%time%calls, &
        & call_min, call_max, call_av,  &
        & total_min, total_max, total_av, total_av*100.0_wp/total_elapsed
      
      if (associated(current_timer%child) .and. depth < max_depth ) then
        current_timer => current_timer%child
        depth = depth + 1
      else if (associated(current_timer%next_sibling)) then
        current_timer => current_timer%next_sibling
      else if (associated(current_timer%parent)) then ! Recurse back up
        do while (associated(current_timer%parent))
          current_timer => current_timer%parent
          depth = depth - 1
          if (associated(current_timer%next_sibling)) then
            current_timer => current_timer%next_sibling
            exit
          end if
        end do
      else
        exit
      end if
      itimer = itimer + 1

    end do

    write(message(itimer), 102) repeat('-',7), "Untimed           ", proc_string, 0, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & 
      & total_elapsed - sum_timed , 100.0_wp - sum_timed*100.0_wp/total_elapsed
    write(message(itimer+1), 100)
    write(message(itimer+2), *) ''

100 format(1X,"+",28("-"),2("+",10("-")),7("+",11("-")),"+")
101 format(1X,"|",12X,"Name",12X,"| Process  ","|  Calls   ","| Call Min  ","| Call Max  ",&
      & "| Call Ave  ","|  Tot Min  ","|  Tot Max  ","|  Tot Ave  ","|     %     ","|")
102 format(1X,"|",1X,A7,1X,A18,1X,"|",1X,A8,1X,"|",1X,I8,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,&
      & "|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",1X,F9.4,1X,"|",2X,F8.4,1X,"|")
    
  end Subroutine timer_print_tree
  
  
  Subroutine init_timer(current_timer, name)
    !! This routine has no header !
    Type ( timer_type_new ) :: current_timer
    Character ( Len = * )  :: name

    current_timer%name  = name
    current_timer%calls = 0
    current_timer%max   = -1.0_wp
    current_timer%min   = huge(1.0_wp)
    current_timer%total = 0.0_wp
    current_timer%last  = huge(1.0_wp)

  end Subroutine init_timer

  Subroutine timer_last_time(tmr, name, screen)
    !! This routine has no header !
    Type ( timer_type ), intent ( in    ) :: tmr
    Character ( Len = * )  :: name
    Logical, Optional :: screen
    Logical :: to_screen
    Character ( Len = 72 ) :: message
    Type ( node ), pointer :: current_timer

    to_screen = .false.
    if (present(screen)) to_screen = screen

    current_timer => find_timer(name)
    if (to_screen) then
      write(0,*) current_timer%time%name, current_timer%time%calls, current_timer%time%last
    else
      write(message,'(a,2(1X,i0))') current_timer%time%name, current_timer%time%calls, current_timer%time%last
      call timer_write(message, tmr)
    end if
  end Subroutine timer_last_time

  Subroutine time_elapsed(tmr)
    !! This routine has no header !
    Type ( timer_type ), Intent( InOut ) :: tmr

    Character( Len = 256 ) :: message

    Call gtime(tmr%elapsed)
    Write(message,'(a,f12.3,a)') "time elapsed since job start: ", tmr%elapsed, " sec"
    Call timer_write(message, tmr)
    
  End Subroutine time_elapsed

  Subroutine timer_write_mul(message, timer_in)
    Character( Len = * ), dimension(:), intent(in) :: message
    Type ( timer_type ), optional, intent ( in    ):: timer_in
    Integer :: i

    if (present(timer_in)) then

      if ( timer_in%proc_id == 0 ) then
        do i = 1, size(message)
          write(timer_in%out_unit,*) message(i)
        end do
      end if
      
    else

      do i = 1, size(message)
        write(0,*) message(i)
      end do
      
    end if
    
  end Subroutine timer_write_mul

  Subroutine timer_write_sing(message, timer_in)
    Character( Len = * ), intent(in) :: message
    Type ( timer_type ), optional, intent ( in    ):: timer_in
    
    if (present(timer_in)) then

      if ( timer_in%proc_id == 0 ) then
          write(timer_in%out_unit,*) message
      end if
      
    else

        write(0,*) message
      
    end if
    
  end Subroutine timer_write_sing

  Subroutine timer_error(message)
    Character ( len = * ) :: message
   
    call timer_write('')
    ! call timer_report(dummy_timer, timer_comm)
    call timer_write(message)
    call timer_write('')
    call dump_call_stack()

    stop
  end Subroutine timer_error

  Subroutine timer_split_stack_string(stack_string, newStack, name)
    Character( Len = * ), intent( in    ) :: stack_string
    Character( Len = 256 ) :: stack
    Type ( call_stack ), intent(   out ) :: newStack
    Character( Len = max_name ), intent(   out ) :: name
    Integer :: i
    Integer :: cnt
    
    stack = adjustl(stack_string)
    cnt = 1
    do i = 1, len(stack)
      if (stack(i:i) == ":") cnt = cnt + 1
    end do
   
    if (cnt > max_depth) call timer_error('Stack depth greater than max depth in timer_split_stack')

    newStack%depth = 0
    do while (index(stack,':') > 0)
      newStack%depth = newStack%depth + 1
      newStack%name(newStack%depth) = trim(stack(1:index(stack,':')-1))

      stack(1:index(stack,':')) = " "
      stack = adjustl(stack)
      
    end do

    name = trim(stack)
          
  end Subroutine timer_split_stack_string

End Module timer
