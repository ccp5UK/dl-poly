program test_hash

  Use hash
  Use units
  Use hashables, only : control_parameter, unit_data
  Use kinds, only : wp
  Implicit None

  Type( hash_table ) :: table
  Type( control_parameter ) :: param, param2
  Character(Len=256) :: input, key, val, units
  Type( unit_data ) :: test
  Real( kind = wp ) :: unit_test

  open(unit=50, file="test_new_control")

  call table%init(256)

  do
     key = ''; val = ''; units = ''
     read(50,'(A)', end=99) input
     call get_word(input, key)
     if (trim(key) == 'title') then
        val = adjustl(input)
     else
        call get_word(input, val)
        if (trim(input) /= '') call get_word(input, units)
     end if
     param = control_parameter(key, val, units, 'Hello', '')
     call table%set(key, param)
  end do

99 continue

  call initialise_units()
  unit_test = convert_units(1.0_wp, 'ns', 'ps')
  print*, unit_test

  print*, convert_units(1.0_wp, 'katm', 'MPa')

contains

  Subroutine get_word(record,word)

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

    Character( Len = * ), Intent( InOut ) :: record
    Character( Len = * ), Intent(   Out ) :: word

    Logical :: transfer
    Integer :: rec_len,word_len,rec_ind,word_ind

    ! Strip blanks in record

    record = adjustl(record)

    ! Get record and word lengths

    rec_len  = Len_Trim(record)
    word_len = Len(word)

    ! Initialise counters and word, and keep-transferring boolean

    rec_ind  = 0
    word_ind = 0

    word     = ' '

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

end program test_hash
