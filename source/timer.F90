!> Routine timer module, provides procedures for calculating real time spent in
!> various program units
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - A. M. Elena May 2018
Module timer
  Use kinds, Only : wp
  Use comms, Only : comms_type,gtime,mtime,gmin,gmax,gsum
  Use errors_warnings, Only : info
  Implicit None

  Private
  Type, Public :: timer_type
    Real( Kind = wp) :: elapsed,job,clear_screen
    Real( Kind = wp) :: t_shortrange(4) = 0.0_wp
    Real( Kind = wp) :: t_longrange(4) = 0.0_wp
    Real( Kind = wp) :: t_linkcell(4) = 0.0_wp
    Real( Kind = wp) :: t_shake(4) = 0.0_wp
    Real( Kind = wp) :: t_rattle(4) = 0.0_wp
    Real( Kind = wp) :: t_bonded(4) = 0.0_wp
  End Type
  Character( Len = 12), Parameter  :: routine(6) = ['short range',' long range','  link cell','      shake','     rattle',&
    '     bonded']

  Public :: start_timer
  Public :: stop_timer
  Public :: time_elapsed
  Public :: timer_report

Contains

  Subroutine start_timer(t)
    Real( kind = wp ), Intent(   Out ) :: t(4)
    Call mtime(t(1))
  End Subroutine start_timer

  Subroutine stop_timer(t)
    Real( kind = wp ), Intent(   Out ) :: t(4)
    Call mtime(t(2))
    t(3) = t(3) + t(2) - t(1)
    t(4) = t(4) + 1.0_wp
  End Subroutine stop_timer

  Subroutine time_elapsed(time)
    Real( Kind = wp ), Intent( InOut ) :: time

    Character( Len = 256 ) :: message

    Call gtime(time)
    Write(message,'(a,f12.3,a)') "time elapsed since job start: ", time, " sec"
    Call info(message,.true.)
  End Subroutine time_elapsed

  Subroutine timer_report(tmr,comm)
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm
    Character( Len = 72 ) :: message(10)

    Real( Kind = wp ) :: buffer(6),mins(6),maxs(6),avg(6)
    Integer :: i

    Write(message(1),'(5a12,a6,a2)') "| Routine   ","| Calls     ","|   Min [s] ","|   Max [s] ","|    Avg [s]","| [%] ","|"
    buffer=0.0_wp
    Do i = 1 , 3
      buffer(1) = tmr%t_shortrange(3)
      buffer(2) = tmr%t_longrange(3)
      buffer(3) = tmr%t_linkcell(3)
      buffer(4) = tmr%t_shake(3)
      buffer(5) = tmr%t_rattle(3)
      buffer(6) = tmr%t_bonded(3)
      Select Case(i)
       Case(1)
        Call gmin(comm,buffer(1:6))
        mins=buffer
       Case(2)
        Call gmax(comm,buffer(1:6))
        maxs=buffer
       Case(3)
        Call gsum(comm,buffer(1:6))
        avg=buffer
      End Select
    Enddo
    buffer(1) = tmr%t_shortrange(4)
    buffer(2) = tmr%t_longrange(4)
    buffer(3) = tmr%t_linkcell(4)
    buffer(4) = tmr%t_shake(4)
    buffer(5) = tmr%t_rattle(4)
    buffer(6) = tmr%t_bonded(4)
    Do i =1,6
      write(message(i+1),'(1x,a11,i12,es12.5,es12.5,es12.5,f6.2)')routine(i),Int(buffer(i)),mins(i),maxs(i),avg(i)/comm%mxnode,&
        avg(i)/comm%mxnode/tmr%elapsed*100.0_wp
    End Do
    write(message(8),'(1x,a11,36x,es12.5,f6.2)')"Other",tmr%elapsed-sum(avg(1:6)/comm%mxnode),&
      100.0_wp-sum(avg(1:6)/comm%mxnode)/tmr%elapsed*100.0_wp
    Call info(message,8,.true.)

  End Subroutine timer_report

End Module timer
