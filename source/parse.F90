Module parse

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module containing tools for parsing textual input
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov may 2004
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use comms,                         Only: comms_type,&
                                           gbcast,&
                                           gcheck,&
                                           gsync
  Use errors_warnings,               Only: error,&
                                           warning
  Use, Intrinsic :: iso_fortran_env, Only: IOSTAT_END
  Use kinds,                         Only: wp

  Implicit None

  Public :: tabs_2_blanks, nls_2_blanks, strip_blanks, get_word, &
            clean_string, lower_case, get_line, word_2_real

Contains

  Subroutine tabs_2_blanks(record)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to convert tabs into blanks in a string
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record

    Integer :: i

    Do i = 1, Len_trim(record)
      If (record(i:i) == Achar(9)) record(i:i) = ' '
    End Do

  End Subroutine tabs_2_blanks

  Subroutine nls_2_blanks_old(record)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to convert new line sequences into blanks in a
    ! string, catching both EoL and CR
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov april 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record

    Integer :: i

    Do i = 1, Len_trim(record)
      If (record(i:i) == Achar(10) .or. &
          record(i:i) == Achar(13)) record(i:i) = ' '
    End Do

  End Subroutine nls_2_blanks_old

  Subroutine nls_2_blanks(record)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to convert new line sequences into blanks in a
    ! string, catching both EoL and CR
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena april 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record

    Integer :: i, j, k

    i = Len_trim(record)
    j = Index(record(1:i), New_line("a"))
    k = Len(New_line("a"))
    Do While (j > 0)
      If (j == i) Then
        record = record(1:j - k)
        j = 0 !Exit
      Else
        record = record(1:j - k)//record(j + 1:i)
        i = i - k
        j = Index(record(1:i), New_line("a")) ! nearly impossible to have a second one
      End If
    End Do

  End Subroutine nls_2_blanks

  Subroutine strip_blanks(record)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to strip blanks from either end of a string
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov july 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record

    record = Adjustl(record)

  End Subroutine strip_blanks

  Subroutine get_word(record, word)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to transfer a word from a string
    !
    ! record loses a word and leading blanks
    ! word fills up with the word as much as it can contain
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record
    Character(Len=*), Intent(  Out) :: word

    Integer :: rec_ind, rec_len, word_ind, word_len
    Logical :: transfer

    ! Strip blanks in record

    Call strip_blanks(record)

    ! Get record and word lengths

    rec_len = Len_trim(record)
    word_len = Len(word)

    ! Initialise counters and word, and keep-transferring boolean

    rec_ind = 0
    word_ind = 0

    word = ' '

    transfer = .true.

    ! Start transferring

    Do While (transfer)

      ! Check for end of record

      If (rec_ind < rec_len) Then

        rec_ind = rec_ind + 1

        ! Check for end of word in record

        If (record(rec_ind:rec_ind) == ' ') transfer = .false.

      Else

        transfer = .false.

      End If

      ! Transfer in word if there is space in word and transfer is true

      If (word_ind < word_len .and. transfer) Then

        word_ind = word_ind + 1

        word(word_ind:word_ind) = record(rec_ind:rec_ind)

        record(rec_ind:rec_ind) = ' '

      Else

        ! Transfer to nothing if there is no space in word and transfer is true

        If (transfer) record(rec_ind:rec_ind) = ' '

      End If

    End Do

  End Subroutine get_word

  Subroutine clean_string(record)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to clean multiple white spacing such as blanks,
    ! tabs and new lines form a string
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov april 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record

    Character(Len=10000) :: record1
    Character(Len=200)   :: word

    ! assign

    record1 = record(1:Len_trim(record))

    ! clean

    Call tabs_2_blanks(record1)
    Call nls_2_blanks(record1)
    Call strip_blanks(record1)

    ! compress for fun

    record = ' '
    Do While (Len_trim(record1) > 0)

      ! read word

      word(1:1) = ' '
      Do While (word(1:1) == ' ')
        Call get_word(record1, word)
      End Do

      ! add cleanly

      If (Len_trim(record) > 0) Then
        record = record(1:Len_trim(record))//' '//word(1:Len_trim(word))
      Else
        record = word(1:Len_trim(word))
      End If
    End Do

  End Subroutine clean_string

  Subroutine lower_case(record)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to lower the character case of a string.
    ! Transportable to non-ASCII machines
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*), Intent(InOut) :: record

    Integer :: i

    Do i = 1, Len(record)
      If (record(i:i) == 'A') Then
        record(i:i) = 'a'
      Else If (record(i:i) == 'B') Then
        record(i:i) = 'b'
      Else If (record(i:i) == 'C') Then
        record(i:i) = 'c'
      Else If (record(i:i) == 'D') Then
        record(i:i) = 'd'
      Else If (record(i:i) == 'E') Then
        record(i:i) = 'e'
      Else If (record(i:i) == 'F') Then
        record(i:i) = 'f'
      Else If (record(i:i) == 'G') Then
        record(i:i) = 'g'
      Else If (record(i:i) == 'H') Then
        record(i:i) = 'h'
      Else If (record(i:i) == 'I') Then
        record(i:i) = 'i'
      Else If (record(i:i) == 'J') Then
        record(i:i) = 'j'
      Else If (record(i:i) == 'K') Then
        record(i:i) = 'k'
      Else If (record(i:i) == 'L') Then
        record(i:i) = 'l'
      Else If (record(i:i) == 'M') Then
        record(i:i) = 'm'
      Else If (record(i:i) == 'N') Then
        record(i:i) = 'n'
      Else If (record(i:i) == 'O') Then
        record(i:i) = 'o'
      Else If (record(i:i) == 'P') Then
        record(i:i) = 'p'
      Else If (record(i:i) == 'Q') Then
        record(i:i) = 'q'
      Else If (record(i:i) == 'R') Then
        record(i:i) = 'r'
      Else If (record(i:i) == 'S') Then
        record(i:i) = 's'
      Else If (record(i:i) == 'T') Then
        record(i:i) = 't'
      Else If (record(i:i) == 'U') Then
        record(i:i) = 'u'
      Else If (record(i:i) == 'V') Then
        record(i:i) = 'v'
      Else If (record(i:i) == 'W') Then
        record(i:i) = 'w'
      Else If (record(i:i) == 'X') Then
        record(i:i) = 'x'
      Else If (record(i:i) == 'Y') Then
        record(i:i) = 'y'
      Else If (record(i:i) == 'Z') Then
        record(i:i) = 'z'
      End If
    End Do

  End Subroutine lower_case

  !  Subroutine get_line(safe,ifile,record,comm)
  !
  !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    !
  !    ! dl_poly_4 subroutine to read a character string on node zero and
  !    ! broadcast it to all other nodes
  !    !
  !    ! copyright - daresbury laboratory
  !    ! author    - i.t.todorov june 2011
  !    ! contrib   - a.m.elena   february 2018
  !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !    Logical,              Intent(   Out ) :: safe
  !    Integer,              Intent( In    ) :: ifile
  !    Character( Len = * ), Intent(   Out ) :: record
  !    Type(comms_type),     Intent( InOut ) :: comm
  !
  !    Integer                              :: i,fail,rec_len,ierr
  !    Integer, Dimension( : ), Allocatable :: line
  !
  !    rec_len = Len(record)
  !
  !    record = ' '
  !    safe = .true.
  !
  !    Call gsync(comm)
  !
  !    If (comm%idnode == 0) Then
  !
  !      Read(Unit=ifile, Fmt='(a)', iostat=ierr) record
  !
  !      If (ierr /= 0 ) Then
  !        safe = .false.
  !      End If
  !    End If
  !
  !    Call gcheck(comm,safe)
  !
  !    If ( .not. safe) Then
  !      Call tabs_2_blanks(record)
  !      Return
  !    End If
  !
  !    If (comm%mxnode > 1) Then
  !      fail = 0
  !
  !      Allocate (line(1:rec_len), Stat = fail)
  !      If (fail > 0) Call error(1011)
  !
  !      line = 0
  !      Do i=1,rec_len
  !        line(i) = Ichar(record(i:i))
  !      End Do
  !
  !      Call gbcast(comm,line,0)
  !
  !      If (comm%idnode > 0) Then
  !        Do i=1,rec_len
  !          record(i:i) = Char(line(i))
  !        End Do
  !      End If
  !
  !      Call tabs_2_blanks(record)
  !      Deallocate (line, Stat = fail)
  !      If (fail > 0) Call error(1012)
  !    End If
  !
  !  End Subroutine get_line

  Subroutine get_line(safe, ifile, record, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to read a character string on node zero and
    ! broadcast it to all other nodes
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2011
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,          Intent(  Out) :: safe
    Integer,          Intent(In   ) :: ifile
    Character(Len=*), Intent(  Out) :: record
    Type(comms_type), Intent(InOut) :: comm

    Integer                            :: fail, i, rec_len
    Integer, Allocatable, Dimension(:) :: line

    rec_len = Len(record)

    fail = 0
    Allocate (line(1:rec_len), Stat=fail)
    If (fail > 0) Call error(1011)

    record = ' '
    safe = .true.

    Call gsync(comm)

    If (comm%idnode == 0) Then

      Read (Unit=ifile, Fmt='(a)', End=100) record

      If (comm%mxnode > 1) Then
        Do i = 1, rec_len
          line(i) = Ichar(record(i:i))
        End Do

        Call gcheck(comm, safe)

        !Call MPI_BCAST(line(1:rec_len), rec_len, MPI_INTEGER, 0, dlp_comm_world, ierr)
        Call gbcast(comm, line, 0)
      End If

      Go To 200

      100 safe = .false.

      Call gcheck(comm, safe)
      If (.not. safe) Go To 200

    Else

      Call gcheck(comm, safe)
      If (.not. safe) Go To 200

      line = 0

      !Call MPI_BCAST(line(1:rec_len), rec_len, MPI_INTEGER, 0, dlp_comm_world, ierr)
      Call gbcast(comm, line, 0)

      Do i = 1, rec_len
        record(i:i) = Char(line(i))
      End Do

    End If

    200 Continue

    Call tabs_2_blanks(record) ; Call strip_blanks(record)

    Deallocate (line, Stat=fail)
    If (fail > 0) Call error(1012)

  End Subroutine get_line

  !> @brief Get number of lines in file
  !!
  !! copyright - daresbury laboratory
  !! author    - a. buccheri 2019
  !!
  !! @param[in]  fname   Name of file
  !! @param[out] nlines  Number of lines in file
  !
  Function number_of_lines(fname) Result(nlines)

    Character(Len=*), Intent(In   ) :: fname
    Integer                         :: nlines

    Character(Len=50) :: error_message
    Integer           :: ios, u

    nlines = 0
    Open (newunit=u, file=Trim(Adjustl(fname)), &
          form='formatted', status='old', iostat=ios)

    If (ios /= 0) Then
      error_message = 'Error openning: '//Trim(Adjustl(fname))
      Call error(0, message=error_message, master_only=.true.)
    End If

    Do
      Read (u, *, iostat=ios)
      If (ios == IOSTAT_END) Exit
      nlines = nlines + 1
    End Do

    Close (u)
  End Function number_of_lines

  Function word_2_real(word, def, report)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 function for extracting real numbers from a character string
    ! with no blanks between the characters of the number.  The optional
    ! argument 'def' suppresses error reporting to return a safe value
    !
    ! (1) Numbers as 2.0e-3/3.d-04 are processable as only one slash is
    !     permitted in the string!
    ! (2) Numbers cannot start or finish with a slash!
    ! (3) A blank string is read as zero!
    ! (4) Numbers must sensible!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len=*),        Intent(In   ) :: word
    Real(Kind=wp), Optional, Intent(In   ) :: def
    Logical, Optional,       Intent(In   ) :: report
    Real(Kind=wp)                          :: word_2_real

    Character(Len=256) :: message
    Character(Len=40)  :: forma
    Integer            :: slash_position, word_end
    Logical            :: l_report
    Real(Kind=wp)      :: denominator

    l_report = .true.
    If (Present(report)) l_report = report

    denominator = 1.0_wp

    word_end = Len_trim(word)
    slash_position = Index(word, '/')

    If (word_end /= 0) Then
      If (slash_position == 1 .or. slash_position == word_end) Go To 30
    Else
      If (Present(def)) Then
        word_2_real = def
      Else
        word_2_real = 0.0_wp
      End If
      Return
    End If

    If (slash_position > 0) Then
      forma = ' '
      Write (forma, '(a,i0,a)') '(f', word_end - slash_position, '.0)'
      Read (word(slash_position + 1:word_end), forma, Err=30) denominator
      word_end = slash_position - 1
    End If

    forma = ' '
    Write (forma, '(a,i0,a)') '(f', word_end, '.0)'
    Read (word(1:word_end), forma, Err=30) word_2_real
    word_2_real = word_2_real / denominator

    Return

    30 Continue
    If (Present(def)) Then
      word_2_real = def
      If (l_report) Then
        Write (message, '(3a,g20.10,a)') &
          "word_2_real defaulted word # ", word(1:word_end), " # to number # ", def, " #"
        Call warning(message, .true.)
      End If
    Else
      word_2_real = 0.0_wp
      Write (message, '(3a)') &
        "word_2_real expected to read a number but found # ", word(1:word_end), " #"
      Call error(1, message, .true.)
    End If

  End Function word_2_real

  Function truncate_real(r)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 function for truncating real numbers to the approximate
    ! precision in decimal digits for the +/-0.___E+/-___ representation,
    ! which is 2*Bit_Size(real)-1
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov october 2005
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(In   ) :: r
    Real(Kind=wp)                :: truncate_real

    Character(Len=100) :: word
    Character(Len=40)  :: forma
    Integer            :: e_position, i, k, word_end

    k = 64 / 4 - 1 ! Bit_Size(0.0_wp)/4 - 1
    Write (forma, '("(0p,e",i0,".",i0,")")') k + 10, k

    word = ' '
    Write (word, forma) r
    Call lower_case(word)
    word_end = Len_trim(word)
    e_position = 0
    e_position = Index(word, 'e')
    Do i = e_position - 3, word_end
      If (i + 3 <= word_end) Then
        word(i:i) = word(i + 3:i + 3)
      Else
        word(i:i) = ' '
      End If
    End Do

    Call strip_blanks(word)
    truncate_real = word_2_real(word)

  End Function truncate_real

End Module parse
