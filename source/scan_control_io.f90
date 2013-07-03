Subroutine scan_control_io()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning the I/O options in the control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
! amended   - i.j.bush october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module,  Only : nread,nrite
  Use config_module, Only : imc_n
  Use parse_module,  Only : get_line,get_word,lower_case,strip_blanks,word_2_real
  Use io_module,     Only : io_set_parameters,        &
                            io_get_parameters,        &
                            io_nc_set_real_precision, &
                            io_nc_compiled,           &
                            IO_READ_MPIIO,            &
                            IO_READ_DIRECT,           &
                            IO_READ_MASTER,           &
                            IO_READ_NETCDF,           &
                            IO_WRITE_UNSORTED_MPIIO,  &
                            IO_WRITE_UNSORTED_DIRECT, &
                            IO_WRITE_UNSORTED_MASTER, &
                            IO_WRITE_SORTED_MPIIO,    &
                            IO_WRITE_SORTED_DIRECT,   &
                            IO_WRITE_SORTED_NETCDF,   &
                            IO_WRITE_SORTED_MASTER

  Implicit None

  Logical                :: carry,safe
  Character( Len = 200 ) :: record,record1
  Character( Len = 40  ) :: word,word1
  Real( Kind = wp )      :: tmp

! Some parameters and variables needed for dealing with I/O options

  Integer            :: io_read,io_write,itmp, err_r
  Logical            :: l_io_r,l_io_w,l_tmp
  Integer, Parameter :: MAX_BATCH_SIZE  = 10000000 !~1GB memory per writer
  Integer, Parameter :: MAX_BUFFER_SIZE =   100000 !~1GB memory per node/domain


! flags

  l_io_r = .false. ! io read  not specified
  l_io_w = .false. ! io write not specified

  safe   = .true.  ! all is safe

! Open the simulation input file

  If (idnode == 0) Inquire(File='CONTROL', Exist=safe)
  If (mxnode > 1) Call gcheck(safe,"enforce")
  If (.not.safe) Then
     Go To 10
  Else
     If (idnode == 0) Open(Unit=nread, File='CONTROL', Status='old')
  End If

! Read TITLE record

  Call get_line(safe,nread,record)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! read slab option
! limiting DD slicing in z direction to 2 for load balancing purposes
! this is really a pre-scan in order to get the MD box dimensions
! when scan_config before the option is read again in scan_control

     If      (word(1:4) == 'slab') Then

        imc_n=6

! io options

     Else If (word(1:2) == 'io' ) Then

        Call get_word( record, word )

        If      (word(1:4) == 'read') Then

           l_io_r  = .true.
           io_read = IO_READ_MPIIO

! get read method

           word1=' ' ; word1=word
           Call get_word( record, word )
           If      ( word( 1:5 ) == 'mpiio'  ) Then
              io_read = IO_READ_MPIIO
              If (idnode == 0) Write(nrite,"(/,1x,'I/O read method: parallel by using MPI-I/O')" )
           Else If ( word( 1:6 ) == 'direct' ) Then
              io_read = IO_READ_DIRECT
              If (idnode == 0) Write(nrite,"(/,1x,'I/O read method: parallel by using direct access')" )
           Else If ( word( 1:6 ) == 'netcdf' ) Then
              io_read = IO_READ_NETCDF
              If (idnode == 0) Write(nrite,"(/,1x,'I/O read method: parallel by using netCDF')" )
           Else If ( word( 1:6 ) == 'master' ) Then
              io_read = IO_READ_MASTER
              If (idnode == 0) Write(nrite,"(/,1x,'I/O read method: serial by using a single master process')" )
           Else
              Call strip_blanks(record)
              If (idnode == 0) Write(nrite,"(/,/,4a)") 'io ',word1(1:Len_Trim(word1)+1),word(1:Len_Trim(word)+1),record
              Call error(3)
           End If

           Call io_set_parameters( user_method_read = io_read )

! get number of readers
! 1 <= readers <= mxnode or be wise, default = 8 or 1 when master
! Note set the default number of readers to 1/4 the number of writers using some (limited)
! empirical evidence from HECToR

           If      (io_read == IO_READ_MPIIO  .or. &
                    io_read == IO_READ_DIRECT .or. &
                    io_read == IO_READ_NETCDF) Then

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real( word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 tmp = Min( Real(mxnode,wp), 2.0_wp*Real(mxnode,wp)**0.5_wp )
                 itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                 Do While ( Mod( mxnode, itmp ) /= 0 )
                    itmp = itmp - 1
                 End Do
                 If (idnode == 0) Write(nrite,"(1x,'I/O readers (assumed) ',9x,i10)") itmp
              Else
                 If (itmp > mxnode) Then
                    tmp = Min( Real(mxnode,wp), 2.0_wp*Real(mxnode,wp)**0.5_wp )
                    itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                    Do While ( Mod( mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    If (idnode == 0) Write(nrite,"(1x,'I/O readers (enforced)',9x,i10)") itmp
                 Else
                    Do While ( Mod( mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    If (idnode == 0) Write(nrite,"(1x,'I/O readers set to    ',9x,i10)") itmp
                 End If
              End If

! the number of readers is now ready to set

              Call io_set_parameters( user_n_io_procs_read = itmp )

! get read batch size
! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
! Note zero or negative values indicate use the default

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real( word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 Call io_get_parameters( user_batch_size_read = itmp )
                 If (idnode == 0) Write(nrite,"(1x,'I/O read batch size (assumed)',2x,i10)") itmp
              Else
                 itmp = Min( itmp, MAX_BATCH_SIZE )
                 Call io_set_parameters( user_batch_size_read = itmp )
                 If (idnode == 0) Write(nrite,"(1x,'I/O read batch size set to   ',2x,i10)") itmp
              End If

           Else If (io_read == IO_READ_MASTER) Then

              If (idnode == 0) Write(nrite,"(1x,'I/O readers (enforced)',9x,i10)") 1

           Else

              tmp = Min( Real(mxnode,wp), 2.0_wp*Real(mxnode,wp)**0.5_wp )
              itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
              If (idnode == 0) Write(nrite,"(1x,'I/O readers (enforced)',9x,i10)") itmp

! the number of readers is now ready to set

              Call io_set_parameters( user_n_io_procs_read = itmp )

           End If

! get read buffer size
! 100 <= buffer <= MAX_BUFFER_SIZE, default 20000
! Note zero or negative values indicate use the default

           Call get_word( record, word )
           itmp = Nint( Abs( word_2_real( word, 0.0_wp ) ) )
           If (itmp == 0) Then
              Call io_get_parameters( user_buffer_size_read = itmp )
              If (idnode == 0) Write(nrite,"(1x,'I/O read buffer size (assumed)',1x,i10)") itmp
           Else
              itmp = Min( Max( itmp,100 ),Min( itmp,MAX_BUFFER_SIZE ))
              Call io_set_parameters( user_buffer_size_read = itmp )
              If (idnode == 0) Write(nrite,"(1x,'I/O read buffer size set to   ',1x,i10)") itmp
           End If

! switch error checking flag for reading

           If (io_read /= IO_READ_MASTER) Then
              Call get_word( record, word )
              l_tmp = ( word( 1:1 ) == 'Y' .or. word( 1:1 ) == 'y' )
              If (.not.l_tmp) Then
                 If (idnode == 0) Write(nrite,"(1x,'I/O parallel read error checking off')")
              Else
                 If (idnode == 0) Write(nrite,"(1x,'I/O parallel read error checking on')")
              End If
              Call io_set_parameters( user_error_check = l_tmp )
           End If

        Else If (word(1:4) == 'writ'  .or. &
                 word(1:5) == 'mpiio' .or. word(1:6) == 'direct' .or. word(1:6) == 'netcdf' .or. word(1:6) == 'master') Then

           l_io_w   = .true.
           io_write = IO_WRITE_SORTED_MPIIO

! for backwards compatibility
! see the second line of the "Else If" above

           If (word(1:4) /= 'writ') Then
              word1=' ' ; word1='write'
           Else
              word1=' ' ; word1=word
              Call get_word( record, word )
           End If

! get write method

           If      ( word( 1:5 ) == 'mpiio'  ) Then
              io_write = IO_WRITE_SORTED_MPIIO
              If (idnode == 0) Write(nrite,"(/,1x,'I/O write method: parallel by using MPI-I/O')" )
           Else If ( word( 1:6 ) == 'direct' ) Then
              io_write = IO_WRITE_SORTED_DIRECT
              If (idnode == 0) Write(nrite,"(/,1x,'I/O write method: parallel by using direct access', &
                                           & /,1x,'*** warning - in parallel this has portability issues !!! ***' )" )
           Else If ( word( 1:6 ) == 'netcdf' ) Then
              io_write = IO_WRITE_SORTED_NETCDF
              Call get_word( record, word ) ! Check if the user wants the "amber-like/32-bit" format
              If ( word( 1:5 ) == 'amber' .or. word( 1:5 ) == '32bit' ) Then
                 ! Use 32-bit quantities in output for real numbers
                 If (idnode == 0) &
                 Write(nrite,"(/,1x,'I/O write method: parallel by using netCDF in the amber-like/32-bit format')" )
                 Call io_nc_set_real_precision( Precision( 1.0 ), Range( 1.0 ), err_r )
              Else
                 ! Use 64-bit quantities in output for real numbers
                 If (idnode == 0) Write(nrite,"(/,1x,'I/O write method: parallel by using netCDF in 64-bit format')" )
                 Call io_nc_set_real_precision( Precision( 1.0d0 ), Range( 1.0d0 ), err_r )
                 record1=' '
                 record1=word(1:Len_Trim(word)+1) // record ! back up
                 record=record1
              End If
           Else If ( word( 1:6 ) == 'master' ) Then
              io_write = IO_WRITE_SORTED_MASTER
              If (idnode == 0) Write(nrite,"(/,1x,'I/O write method: serial by using a single master process')" )
           Else
              Call strip_blanks(record)
              If (idnode == 0) Write(nrite,"(/,/,4a)") 'io ',word1(1:Len_Trim(word1)+1),word(1:Len_Trim(word)+1),record
              Call error(3)
           End If

! get write type

           Call get_word( record, word )
           If      ( word( 1:6 ) == 'unsort' ) Then
              If (idnode == 0) Write(nrite,"(1x,'I/O write type: data sorting off')")
              Select Case( io_write )
              Case( IO_WRITE_SORTED_MPIIO  )
                 io_write = IO_WRITE_UNSORTED_MPIIO
              Case( IO_WRITE_SORTED_DIRECT )
                 io_write = IO_WRITE_UNSORTED_DIRECT
              Case( IO_WRITE_SORTED_MASTER )
                 io_write = IO_WRITE_UNSORTED_MASTER
              End Select
           Else
              If ( word( 1:6 ) == 'sorted' ) Then
                 If (idnode == 0) Write(nrite,"(1x,'I/O write type: data sorting on')")
              Else
                 record1=' '
                 record1=word(1:Len_Trim(word)+1) // record ! back up
                 record=record1
                 If (idnode == 0) Write(nrite,"(1x,'I/O write type: data sorting on (assumed)')")
              End If
           End If

! the write method and type are now ready to set

           Call io_set_parameters( user_method_write = io_write )

! get number of writers
! 1 <= writers <= mxnode or be wise, default = 8 or 1 when master

           If      (io_write == IO_WRITE_SORTED_MPIIO  .or. &
                    io_write == IO_WRITE_SORTED_DIRECT .or. &
                    io_write == IO_WRITE_SORTED_NETCDF ) Then

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real( word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 tmp = Min( Real(mxnode,wp), 8.0_wp*Real(mxnode,wp)**0.5_wp )
                 itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                 Do While ( Mod( mxnode, itmp ) /= 0 )
                    itmp = itmp - 1
                 End Do
                 If (idnode == 0) Write(nrite,"(1x,'I/O writers (assumed) ',9x,i10)") itmp
              Else
                 If (itmp > mxnode) Then
                    tmp = Min( Real(mxnode,wp), 8.0_wp*Real(mxnode,wp)**0.5_wp )
                    itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                    Do While ( Mod( mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    If (idnode == 0) Write(nrite,"(1x,'I/O writers (enforced)',9x,i10)") itmp
                 Else
                    Do While ( Mod( mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    If (idnode == 0) Write(nrite,"(1x,'I/O writers set to    ',9x,i10)") itmp
                 End If
              End If

! the number of writers is now ready to set

              Call io_set_parameters( user_n_io_procs_write = itmp )

! get write batch size
! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
! Note zero or negative values indicate use the default

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real( word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 Call io_get_parameters( user_batch_size_write = itmp )
                 If (idnode == 0) Write(nrite,"(1x,'I/O write batch size (assumed)',1x,i10)") itmp
              Else
                 itmp = Min( itmp, MAX_BATCH_SIZE )
                 Call io_set_parameters( user_batch_size_write = itmp )
                 If (idnode == 0) Write(nrite,"(1x,'I/O write batch size set to   ',1x,i10)") itmp
              End If


           Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. io_write == IO_WRITE_SORTED_MASTER) Then

              If (idnode == 0) Write(nrite,"(1x,'I/O writers (enforced)',9x,i10)") 1

           Else

              tmp = Min( Real(mxnode,wp), 8.0_wp*Real(mxnode,wp)**0.5_wp )
              itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
              If (idnode == 0) Write(nrite,"(1x,'I/O writers (enforced)',9x,i10)") itmp

! the number of writers is now ready to set

              Call io_set_parameters( user_n_io_procs_write = itmp )

           End If

! get write buffer size
! 100 <= buffer <= MAX_BUFFER_SIZE, default 20000
! Note zero or negative values indicate use the default

           Call get_word( record, word )
           itmp = Nint( Abs( word_2_real( word, 0.0_wp ) ) )
           If (itmp == 0) Then
              Call io_get_parameters( user_buffer_size_write = itmp )
              If (idnode == 0) Write(nrite,"(1x,'I/O write buffer size (assumed)',i10)") itmp
           Else
              itmp = Min( Max( itmp,100 ),Min( itmp,MAX_BUFFER_SIZE ))
              Call io_set_parameters( user_buffer_size_write = itmp )
              If (idnode == 0) Write(nrite,"(1x,'I/O write buffer size set to   ',i10)") itmp
           End If

! switch error checking flag for writing

           If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
              Call get_word( record, word )
              l_tmp = ( word( 1:1 ) == 'Y' .or. word( 1:1 ) == 'y' )
              If (.not.l_tmp) Then
                 If (idnode == 0) Write(nrite,"(1x,'I/O parallel write error checking off')")
              Else
                 If (idnode == 0) Write(nrite,"(1x,'I/O parallel write error checking on')")
              End If
              Call io_set_parameters( user_error_check = l_tmp )
           End If

        Else

           Call strip_blanks(record)
           If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
           Call error(3)

        End If

! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     End If

  End Do

  If (idnode == 0) Close(Unit=nread)

!!! IO DEFAULTS

! io option defaults

  If (.not.l_io_r) Then

! read method

     Call io_set_parameters( user_method_read = IO_READ_MPIIO ) ; io_read = IO_READ_MPIIO
     If (idnode == 0) Write(nrite,"(/,1x,'I/O read method: parallel by using MPI-I/O (assumed)')")

! number of readers

     tmp = Min( Real(mxnode,wp), 2.0_wp*Real(mxnode,wp)**0.5_wp )
     itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
     Do While ( Mod( mxnode, itmp ) /= 0 )
        itmp = itmp - 1
     End Do
     Call io_set_parameters( user_n_io_procs_read = itmp )
     If (idnode == 0) Write(nrite,"(1x,'I/O readers (assumed) ',9x,i10)") itmp

! read batch size

     Call io_get_parameters( user_batch_size_read = itmp )
     If (idnode == 0) Write(nrite,"(1x,'I/O read batch size (assumed)',2x,i10)") itmp

! read buffer size

     Call io_get_parameters( user_buffer_size_read = itmp )
     If (idnode == 0) Write(nrite,"(1x,'I/O read buffer size (assumed)',1x,i10)") itmp

! error checking flag for reading

     If (io_read /= IO_READ_MASTER) Then
        Call io_set_parameters( user_error_check = .false. )
        If (idnode == 0) Write(nrite,"(1x,'I/O parallel read error checking off (assumed)')")
     End If

  End If

  If (.not.l_io_w) Then

! write method

     Call io_set_parameters( user_method_write = IO_WRITE_SORTED_MPIIO ) ; io_write = IO_WRITE_SORTED_MPIIO
     If (idnode == 0) Write(nrite,"(/,1x,'I/O write method: parallel by using MPI-I/O (assumed)')")

! write type

     If (idnode == 0) Write(nrite, "(1x,'I/O write type: data sorting on (assumed)')")

! number of writers

     tmp = Min( Real(mxnode,wp), 8.0_wp*Real(mxnode,wp)**0.5_wp )
     itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
     Do While ( Mod( mxnode, itmp ) /= 0 )
        itmp = itmp - 1
     End Do
     Call io_set_parameters( user_n_io_procs_write = itmp )
     If (idnode == 0) Write(nrite,"(1x,'I/O writers (assumed) ',9x,i10)") itmp

! batch size

     Call io_get_parameters( user_batch_size_write = itmp )
     If (idnode == 0) Write(nrite,"(1x,'I/O write batch size (assumed)',1x,i10)") itmp

! write buffer size

     Call io_get_parameters( user_buffer_size_write = itmp )
     If (idnode == 0) Write(nrite,"(1x,'I/O write buffer size (assumed)',i10)") itmp

! error checking flag for writing

     If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
        Call io_set_parameters( user_error_check = .false. )
        If (idnode == 0) Write(nrite,"(1x,'I/O parallel write error checking off (assumed)')")
     End If

  End If

  If (io_write == IO_WRITE_SORTED_NETCDF .or. io_read == IO_READ_NETCDF) Call io_nc_compiled()

  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control_io
