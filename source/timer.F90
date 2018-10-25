Module timer
  !! This module has no header !
  Use kinds, Only : wp
  Use comms, Only : comms_type,gtime,mtime,gmin,gmax,gsum,gsync,gsend,grecv,timer_tag
  Use errors_warnings, Only : info, error
  Implicit None

  Private

  Integer, Parameter :: max_depth = 20, max_timers = 20, max_name = 18

  Type, Public :: timer_type
    Real( Kind = wp) :: elapsed,job,clear_screen
    Logical :: proc_detail = .false.
  End Type timer_type

  Type :: timer_type_new
    Character( Len = max_name ) :: name = 'UninitialisedTimer'
    Integer           :: id
    Real( Kind = wp ) :: max, min, total, last
    Real( Kind = wp ) :: start, stop
    Integer           :: calls
    Logical           :: running = .false.
  end type timer_type_new

  Type :: call_stack
    Character ( Len = max_name ), dimension( max_depth ) :: name
    Integer :: depth
  end type call_stack

  Type :: node
    Type ( timer_type_new ) :: time
    Type ( node ), pointer :: child => null()
    Type ( node ), pointer :: parent => null()
    Type ( node ), pointer :: next_sibling => null()
  end type node

  Type( timer_type_new ), dimension(0:max_timers), save :: timers
  Type( call_stack ),                              save :: calls
  Type( node ),                                    save :: call_tree

  Public :: start_timer_new
  Public :: stop_timer_new
  Public :: timer_report_new

  Public :: time_elapsed

Contains

  ! Subroutine push_stack ( name )
  !   Character ( Len = * ) :: name

  !   calls%depth = calls%depth + 1
  !   if ( calls%depth > max_depth ) &
  !     & call error(0,'Call stack exceeds max depth : recursive call or unended timer?')
  !   calls%name(calls%depth) = name

  ! end Subroutine push_stack

  ! Subroutine pop_stack ( name )
  !   Character ( Len = * ) :: name

  !   if ( name /= calls%name(calls%depth)) &
  !     & call error(0,'Child timer ended before parent')
  !   calls%name(calls%depth) = ''
  !   calls%depth = calls%depth - 1

  ! end Subroutine pop_stack

  ! Function find_timer_tree ( name )
  !   Character ( Len = * ) :: name
  !   Type ( node ), Pointer :: current
  !   Integer :: depth

  !   depth = 0
  !   current => call_tree%child

  !   do while ( depth /= calls%depth-1 )

  !     do while ( current%name /= calls(depth)%name)
  !       if ( .not. associated(call_tree%next_sibling)) call error(0,'Call stack does not match call tree')
  !       current => current%next_sibling
  !     end do
  !     depth = depth + 1
  !     current => current%child
  !   end do


  ! end Function find_timer_tree

  Function find_timer(name) result(timer)
    !! This routine has no header !
    Integer:: timer
    Character ( Len = * )  :: name

    do timer = 0, max_timers
      if ( timers(timer)%name == name ) then
        return
      end if
      if ( timers(timer)%name == 'UninitialisedTimer' ) then
        call init_timer_new(timers(timer), name)
        timers(timer)%id = timer
        return
      end if
    end do

    call error(0,'Max Timers Exceeded')

  end Function find_timer

  Subroutine init_timer_new(timer, name)
    !! This routine has no header !
    Type ( timer_type_new ) :: timer
    Character ( Len = * )  :: name

    timer%name  = name
    timer%calls = 0
    timer%max   = -1.0_wp
    timer%min   = huge(1.0_wp)
    timer%total = 0.0_wp
    timer%last  = huge(1.0_wp)

  end Subroutine init_timer_new

  Subroutine start_timer_new(name)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Integer :: timer

    timer = find_timer(name)
    Call mtime(timers(timer)%start)
    timers(timer)%running = .true.

  end Subroutine start_timer_new

  Subroutine stop_timer_new(name)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Integer :: timer

    timer = find_timer(name)
    if ( .not. timers(timer)%running ) call error(0,'Timer '//trim(timers(timer)%name)//' stopped but not running')

    Call mtime(timers(timer)%stop)

    timers(timer)%running = .false.
    timers(timer)%last = timers(timer)%stop - timers(timer)%start
    if ( timers(timer)%last > timers(timer)%max ) timers(timer)%max = timers(timer)%last
    if ( timers(timer)%last < timers(timer)%min ) timers(timer)%min = timers(timer)%last
    timers(timer)%total = timers(timer)%total + timers(timer)%last
    timers(timer)%calls = timers(timer)%calls + 1

  End Subroutine stop_timer_new

  Subroutine timer_last_time(name, screen)
    !! This routine has no header !
    Character ( Len = * )  :: name
    Logical, Optional :: screen
    Logical :: to_screen
    Character ( Len = 72 ) :: message
    Integer :: timer

    if (present(screen)) to_screen = screen

    timer = find_timer(name)
    if (to_screen) then
      write(*,*) timers(timer)%name, timers(timer)%calls, timers(timer)%last
    else
      write(message,'(a,2(1X,i0))') timers(timer)%name, timers(timer)%calls, timers(timer)%last
      call info(message,.true.)
    end if
  end Subroutine timer_last_time

  Subroutine timer_report_new(comm,node_detail)
    !! This routine has no header !
    Type( comms_type ), Intent( InOut ) :: comm
    Logical, optional,  Intent( In    ) :: node_detail
    Character( Len = 132 ), dimension(-2:max_timers + 2) :: message
    Real ( kind = wp ) :: call_min, call_max, call_av
    Real ( kind = wp ) :: total_min, total_max, total_av
    Real ( kind = wp ) :: total_elapsed, sum_timed
    Integer :: node
    Integer :: i

    call gtime(total_elapsed)
    write(message(-2), 100)
    write(message(-1), 101)

    sum_timed = 0.0_wp

    do i = 0, max_timers

      if (timers(i)%name == 'UninitialisedTimer') exit
      if (timers(i)%running) Call error(0,'Program terminated while timer '//trim(timers(i)%name)//' still running')
      total_min = timers(i)%total
      total_max = timers(i)%total
      total_av  = timers(i)%total
      Call gmin(comm,total_min)
      Call gmax(comm,total_max)
      Call gsum(comm,total_av)
      total_av = total_av / comm%mxnode
      if (i > 0) sum_timed = sum_timed + total_av

      call_min  = timers(i)%min
      call_max  = timers(i)%max
      Call gmin(comm,call_min)
      Call gmax(comm,call_max)
      call_av   = total_av/timers(i)%calls
      write(message(i), 102 ) timers(i)%name, timers(i)%calls, &
        & call_min, call_max, call_av,  &
        & total_min, total_max, total_av, total_av*100.0_wp/total_elapsed

    end do

    !! JW952
    ! Untimed does not account for children
    ! write(message(i), 103) "Untimed           ", & 
    !      & total_elapsed - sum_timed , 100.0_wp - sum_timed*100.0_wp/total_elapsed
    write(message(i), 100)
    call info(message,i+3,.true.)

    call info('',.true.)


    if (present(node_detail)) then
      if (node_detail) then
        write(message(-2),200)
        write(message(-1),201)
        call info(message,2,.true.)
        do node = 0, comm%mxnode-1
          if (comm%idnode == node) then
            write(message(-2), 100)
            write(message(-1), 201)
            do i = 0, max_timers

              if (timers(i)%name == 'UninitialisedTimer') exit

              total_av  = timers(i)%total
              sum_timed = sum_timed + total_av
              call_min  = timers(i)%min
              call_max  = timers(i)%max
              call_av   = total_av/timers(i)%calls
              write(message(i), 202 ) timers(i)%name, node, timers(i)%calls, &
                & call_min, call_max, call_av, total_av, total_av*100.0_wp/total_elapsed

            End do
            write(message(i), 200)
            if (node > 0 ) call gsend(comm, message, 0, timer_tag)
          end if
          if (comm%idnode == 0 .and. node > 0) then
            call grecv(comm, message, node, timer_tag)
            call info(message(0:),i+1,.true.)
          else if (node == 0) then
            call info(message(0:),i+1,.true.)
          end if

          call gsync(comm)
        end do
      end if
    end if


100 format("+",20("-"),4("+",10("-")),3("+",11("-")),"+",9("-"),"+")
101 format("|",8X,"Name",8X,"|   Calls  ","| Call Min ","| Call Max ","| Call Ave ", &
      & "|  Tot Min  ","|  Tot Max  ","|  Tot Ave  ","|    %    ","|")
102 format("|",1X,A18,1X,"|",1X,I8.1,1X,3("|",1X,F8.4,1X),3("|",1X,F9.4,1X),"|",1X,F7.3,1X,"|")

200 format("+",20("-"),5("+",10("-")),"+",11("-"),"+",9("-"),"+")
201 format("|",8X,"Name",8X,"| Process  ","|   Calls  ","| Call Min ","| Call Max ","| Call Ave ", &
      & "|   Total   ","|    %    ","|")
202 format("|",1X,A18,1X,2("|",1X,I8.1,1X),3("|",1X,F8.4,1X),"|",1X,F9.4,1X,"|",1X,F7.3,1X,"|")

  End Subroutine timer_report_new

  Subroutine time_elapsed(time)
    !! This routine has no header !
    Real( Kind = wp ), Intent( InOut ) :: time

    Character( Len = 256 ) :: message

    Call gtime(time)
    Write(message,'(a,f12.3,a)') "time elapsed since job start: ", time, " sec"
    Call info(message,.true.)
  End Subroutine time_elapsed

End Module timer
