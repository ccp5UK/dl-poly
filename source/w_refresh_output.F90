!!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

     i=20
     If (nstep > 0) Then
        If ( Mod(nstep,i*nstbpo) == 0 .or.                        &
             (tmr%elapsed > Real(i*60,wp) .and.                        &
              tmr%elapsed-Real( ((Int(tmr%elapsed)/(i*60)) * i*60) , wp ) < &
              tmr%elapsed/Real( nstep , wp) ) ) Then

           If (comm%idnode == 0) Then
              Inquire(File=files(FILE_OUTPUT)%filename, Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                 Close(unit=files(FILE_OUTPUT)%unit_no)
                 Open(Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Position='append')
              End If
           End If

        End If
     End If


!!!!!!!!!!!!!!!!!!!!  W_REFRESH_OUTPUT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
