Module timer
  Use kinds, Only : wp 
  Use comms, Only : comms_type,gtime,mtime,gmin,gmax,gsum
  Use errors_warnings, Only : info
  Implicit None

  Private
  Type, Public :: timer_type
    Real( Kind = wp) :: timelp,timjob,timcls
    Real( Kind = wp) :: t_shortrange(4) = 0.0_wp
    Real( Kind = wp) :: t_longrange(4) = 0.0_wp
    Real( Kind = wp) :: t_linkcell(4) = 0.0_wp
  End Type

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
    Real( Kind = wp ), Intent( In    ) :: time

    Character( Len = 256 ) :: message

    Write(message,'(a,f12.3,a)') "time elapsed since job start: ", time, " sec"
    Call info(message,.true.)
  End Subroutine time_elapsed

  Subroutine timer_report(tmr,comm)
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm
    Character( Len = 72 ) :: message(10)

    Real( Kind = wp ) :: buffer(10),mins(10),maxs(10),avg(10)
    Character( Len = 12) :: routine(3) = ['short range',' long range','  link cell']
    Integer :: i

    Write(message(1),'(5a12,a6,a2)') "| Routine   ","| Calls     ","|   Min [s] ","|   Max [s] ","|    Avg [s]","| [%] ","|"
    buffer=0.0_wp
    Do i = 1 , 3
      buffer(1) = tmr%t_shortrange(3)
      buffer(2) = tmr%t_longrange(3)
      buffer(3) = tmr%t_linkcell(3)
      Select Case(i)
         Case(1)
           Call gmin(comm,buffer(1:3))
           mins=buffer
         Case(2)
           Call gmax(comm,buffer(1:3))
           maxs=buffer
         Case(3)
           Call gsum(comm,buffer(1:3))
           avg=buffer
       End Select
     Enddo
      buffer(1) = tmr%t_shortrange(4)
      buffer(2) = tmr%t_longrange(4)
      buffer(3) = tmr%t_linkcell(4)
     Do i =1,3
        write(message(i+1),'(1x,a11,i12,es12.5,es12.5,es12.5,f6.2)')routine(i),Int(buffer(i)),mins(i),maxs(i),avg(i)/comm%mxnode,&
          avg(i)/comm%mxnode/tmr%timelp*100.0_wp
     End Do 
     write(message(5),'(1x,a11,36x,es12.5,f6.2)')"Other",tmr%timelp-sum(avg(1:3)/comm%mxnode),&
          100.0_wp-sum(avg(1:3)/comm%mxnode)/tmr%timelp*100.0_wp
    Call info(message,5,.true.)

  End Subroutine timer_report

End Module timer
