!!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

     i=20
     If (nstep > 0) Then
        If ( Mod(nstep,i*nstbpo) == 0 .or.                        &
             (tmr%timelp > Real(i*60,wp) .and.                        &
              tmr%timelp-Real( ((Int(tmr%timelp)/(i*60)) * i*60) , wp ) < &
              tmr%timelp/Real( nstep , wp) ) ) Then

           If (comm%idnode == 0) Then
              Inquire(File=Trim(output), Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                 Close(Unit=nrite)
                 Open(Unit=nrite, File=Trim(output), Position='append')
              End If
           End If

        End If
     End If


!!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
