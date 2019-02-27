Module timer
  !! This module has no header !
  Use kinds, Only : wp
  Use comms, Only : comms_type,gtime,mtime,gmin,gmax,gsum,gsync,gsend,grecv,timer_tag,abort_comms
  Implicit None

  Private

  Integer, Parameter :: max_depth = 6, max_timers = 50, max_name = 18

  Type, Public :: timer_type
    Real( Kind = wp) :: elapsed,job,clear_screen
    Logical :: proc_detail = .false.
    Integer :: max_depth = 1
  End Type timer_type

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
    Integer :: id
    Type ( node ), pointer :: child => null()
    Type ( node ), pointer :: parent => null()
    Type ( node ), pointer :: next_sibling => null()
    Logical :: done = .false.
  end type node

  Type( timer_type ),                              save :: dummy_timer
  Type( call_stack ),                              save :: calls
  Type( node ), target,                            save :: call_tree
  Integer, save :: id = 1

  Type ( comms_type ), save :: timer_comm
  Integer, save :: out_unit

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

  subroutine init_timer_system ( nrite, comm )

    Type ( comms_type ), intent ( in ) :: comm
    Integer, intent ( in ) :: nrite

    timer_comm = comm
    out_unit = nrite
    dummy_timer%max_depth = huge(1)
    dummy_timer%proc_detail = .false.
    call init_timer ( call_tree%time, 'Head')
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
    
    current => call_tree
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
    if ( id > max_timers ) call timer_error('Too many timers active')
    parent%child%id = id
    id = id + 1
    parent%child%parent => parent

    call init_timer(parent%child%time, name)

  end subroutine init_child_node

  subroutine init_sibling_node(name, sibling)
    Character ( len = * ) :: name
    Type ( node ) :: sibling
    Type ( node ), pointer :: child

    allocate ( sibling%next_sibling )
    child => sibling%next_sibling
    if ( id > max_timers ) call timer_error('Too many timers active')
    child%id = id
    id = id + 1
    child%parent => sibling%parent
    call init_timer(child%time, name)

  end subroutine init_sibling_node

  Subroutine start_timer(name)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Type ( node ), pointer :: timer

    timer => find_timer(name)

    call push_stack(name)
    Call mtime(timer%time%start)
    timer%time%running = .true.

  end Subroutine start_timer

  Subroutine stop_timer(name)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Type ( node ), pointer :: timer


    timer => find_timer(name)

    call pop_stack(name)
    if ( .not. timer%time%running ) call timer_error('Timer '//trim(timer%time%name)//' stopped but not running')

    Call mtime(timer%time%stop)

    timer%time%running = .false.
    timer%time%last = timer%time%stop - timer%time%start
    if ( timer%time%last > timer%time%max ) timer%time%max = timer%time%last
    if ( timer%time%last < timer%time%min ) timer%time%min = timer%time%last
    timer%time%total = timer%time%total + timer%time%last
    timer%time%calls = timer%time%calls + 1

  End Subroutine stop_timer

  Subroutine start_timer_path(name_in)
    !! This routine has no header !
    Character ( Len = * )  :: name_in
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: timer
    integer :: i

    call timer_split_stack_string(name_in, stack, name)
    timer => find_timer(name, stack)
    
    Call mtime(timer%time%start)
    timer%time%running = .true.

  end Subroutine start_timer_path

  Subroutine stop_timer_path(name_in)
    !! This routine has no header !
    Character ( Len = * )  :: name_in
    Character ( Len = max_name ) :: name
    Type ( call_stack ) :: stack
    Type ( node ), pointer :: timer

    call timer_split_stack_string(name_in, stack, name)
    timer => find_timer(name, stack)

    if ( .not. timer%time%running ) call timer_error('Timer '//trim(timer%time%name)//' stopped but not running')

    Call mtime(timer%time%stop)

    timer%time%running = .false.
    timer%time%last = timer%time%stop - timer%time%start
    if ( timer%time%last > timer%time%max ) timer%time%max = timer%time%last
    if ( timer%time%last < timer%time%min ) timer%time%min = timer%time%last
    timer%time%total = timer%time%total + timer%time%last
    timer%time%calls = timer%time%calls + 1

  End Subroutine stop_timer_path
  
  Subroutine timer_report(tmr,comm)
    !! This routine has no header !
    Type( timer_type ), intent( In    ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm
    Character( Len = 132 ), dimension(-2:2) :: message
    Real ( kind = wp ) :: call_min, call_max, call_av
    Real ( kind = wp ) :: total_min, total_max, total_av
    Real ( kind = wp ) :: total_elapsed, sum_timed
    Character ( len = 7 ) :: depth_symb
    Type ( node ), pointer :: timer
    Integer :: proc
    Integer :: i, depth

    call stop_timer('Main')
    write(message(-2), 100)
    write(message(-1), 101)
    call timer_write(message(-2:-1))
    
    sum_timed = 0.0_wp

    timer => call_tree%child
    total_elapsed = timer%time%total
    call gmax(comm,total_elapsed)

    i = 0
    depth = 0

    do while (depth > -1)
      if (timer%time%running) Call timer_write('Program terminated while timer '//&
        & trim(timer%time%name)//' still running')

      if ( depth == 0 .and. timer%time%name /= "Main") then
        depth_symb = repeat('-', 7)
      else if ( associated(timer%child) ) then
        depth_symb = repeat(" ",depth)//"|v"
      else
        depth_symb = repeat(" ",depth)//"|-"
      end if

      total_min = timer%time%total
      total_max = timer%time%total
      total_av  = timer%time%total
      Call gmin(comm,total_min)
      Call gmax(comm,total_max)
      Call gsum(comm,total_av)
      total_av = total_av / comm%mxnode
      if (depth == 1) sum_timed = sum_timed + total_av

      call_min  = timer%time%min
      call_max  = timer%time%max
      Call gmin(comm,call_min)
      Call gmax(comm,call_max)
      call_av   = total_av/timer%time%calls

      write(message(0), 102 ) depth_symb,timer%time%name, timer%time%calls, &
        & call_min, call_max, call_av,  &
        & total_min, total_max, total_av, total_av*100.0_wp/total_elapsed
      timer%done = .true.
      call timer_write(message(0))
      
      if (associated(timer%child) .and. depth < tmr%max_depth ) then
        timer => timer%child
        depth = depth + 1
      else if (associated(timer%next_sibling)) then
        timer => timer%next_sibling
      else if (associated(timer%parent)) then ! Recurse back up
        do while (associated(timer%parent))
          timer => timer%parent
          depth = depth - 1
          if (associated(timer%next_sibling)) then
            timer => timer%next_sibling
            exit
          end if
        end do
      else
        exit
      end if

    end do

    write(message(0), 103) "Untimed           ", & 
      & total_elapsed - sum_timed , 100.0_wp - sum_timed*100.0_wp/total_elapsed
    write(message(1), 100)
    call timer_write(message(0:1))

    call timer_write('')

    if (tmr%proc_detail) then
      write(message(-2),200)
      write(message(-1),201)
      call timer_write(message(-2:-1))

      depth = 0 

      do proc = 0, comm%mxnode-1
        if (comm%idnode == proc) then
          write(message(-2), 100)
          write(message(-1), 201)

          call_tree%parent => call_tree
          timer => call_tree
          i = 0
          do while (associated(timer%parent))
            nullify(call_tree%parent)

            if ( associated(timer%child) ) then
              depth_symb = repeat(" ",depth)//"|v"
            else
              depth_symb = repeat(" ",depth)//"|-"
            end if
            total_av  = timer%time%total
            sum_timed = sum_timed + total_av
            call_min  = timer%time%min
            call_max  = timer%time%max
            call_av   = total_av/timer%time%calls
            write(message(i), 202 ) depth_symb,timer%time%name, proc, timer%time%calls, &
              & call_min, call_max, call_av, total_av, total_av*100.0_wp/total_elapsed

            if (associated(timer%child) .and. depth < tmr%max_depth) then
              timer => timer%child
              depth = depth + 1

            else if (associated(timer%next_sibling)) then
              timer => timer%next_sibling

            else
              do while (associated(timer%parent))
                timer => timer%parent
                depth = depth - 1

                if (associated(timer%next_sibling)) then
                  timer => timer%next_sibling
                  exit
                end if
              end do
            end if
            i = i + 1

          End do
          write(message(i), 200)
          if (proc > 0) call gsend(comm, message, 0, timer_tag)
        end if
        if (comm%idnode == 0 .and. proc > 0) then
          call grecv(comm, message, proc, timer_tag)
          call timer_write(message(0:i))
        else if (proc == 0) then
          call timer_write(message(0:i))
        end if

        call gsync(comm)
      end do
    end if

100 format("+",28("-"),4("+",10("-")),3("+",11("-")),"+",9("-"),"+")
101 format("|",12X,"Name",12X,"|   Calls  ","| Call Min ","| Call Max ","| Call Ave ", &
      & "|  Tot Min  ","|  Tot Max  ","|  Tot Ave  ","|    %    ","|")
102 format("|",1X,A7,1X,A18,1X,"|",1X,I8.1,1X,3("|",1X,F8.4,1X),3("|",1X,F9.4,1X),"|",1X,F7.3,1X,"|")
103 format("|",1X,7("-"),1X,A18,1X,4("|",10X),2("|",11X),"|",1X,F9.4,1X,"|",1X,F7.3,1X,"|")

200 format("+",28("-"),5("+",10("-")),"+",11("-"),"+",9("-"),"+")
201 format("|",12X,"Name",12X,"| Process  ","|   Calls  ","| Call Min ","| Call Max ","| Call Ave ", &
      & "|   Total   ","|    %    ","|")
202 format("|",1X,A7,1X,A18,1X,2("|",1X,I8.1,1X),3("|",1X,F8.4,1X),"|",1X,F9.4,1X,"|",1X,F7.3,1X,"|")

  End Subroutine timer_report

  Subroutine init_timer(timer, name)
    !! This routine has no header !
    Type ( timer_type_new ) :: timer
    Character ( Len = * )  :: name

    timer%name  = name
    timer%calls = 0
    timer%max   = -1.0_wp
    timer%min   = huge(1.0_wp)
    timer%total = 0.0_wp
    timer%last  = huge(1.0_wp)

  end Subroutine init_timer

  Subroutine timer_last_time(name, screen)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Logical, Optional :: screen
    Logical :: to_screen
    Character ( Len = 72 ) :: message
    Type ( node ), pointer :: timer

    to_screen = .false.
    if (present(screen)) to_screen = screen

    timer => find_timer(name)
    if (to_screen) then
      write(*,*) timer%time%name, timer%time%calls, timer%time%last
    else
      write(message,'(a,2(1X,i0))') timer%time%name, timer%time%calls, timer%time%last
      call timer_write(message)
    end if
  end Subroutine timer_last_time

  Subroutine time_elapsed(time)
    !! This routine has no header !
    Real( Kind = wp ), Intent( InOut ) :: time

    Character( Len = 256 ) :: message

    Call gtime(time)
    Write(message,'(a,f12.3,a)') "time elapsed since job start: ", time, " sec"
    Call timer_write(message)
  End Subroutine time_elapsed

  Subroutine timer_write_sing(message)
    Character( Len = * ), intent(in) :: message

    if ( timer_comm%idnode == 0 ) &
      write(out_unit,*) message

  end Subroutine timer_write_sing

  Subroutine timer_write_mul(message)
    Character( Len = * ), dimension(:), intent(in) :: message
    Integer :: i

    if ( timer_comm%idnode == 0 ) then
      do i = 1, size(message)
        write(out_unit,*) message(i)
      end do
    end if

  end Subroutine timer_write_mul

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
