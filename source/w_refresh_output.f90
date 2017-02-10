!!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

     i=20
     If (nstep > 0) Then
        If ( Mod(nstep,i*nstbpo) == 0 .or.                        &
             (timelp > Real(i*60,wp) .and.                        &
              timelp-Real( ((Int(timelp)/(i*60)) * i*60) , wp ) < &
              timelp/Real( nstep , wp) ) ) Then

           If (idnode == 0) Then
              Inquire(File=trim(output), Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                 Close(Unit=nrite)
                 Open(Unit=nrite, File=trim(output), Position='append')
              End If
           End If

        End If
     End If


!!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
