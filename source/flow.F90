Module flow
  Use kinds, Only : wi
  Implicit None

  Private

  !> Type containing program flow data
  Type, Public :: flow_type
    Private
    !> Check if is first time we call build_book_intra
    Logical, Public :: newjob_build_book = .true.
    Logical, Public :: oldjob_shared_units = .false.

    ! STDOUT printing control
    !> Number of print events before starting a new 'page'
    Integer(Kind=wi) :: npage = 8
    !> Current number of print events
    Integer(Kind=wi), Public :: lines = 0

    !> Check if first call of md_vv or calculate_forces
    Logical, Public :: newjob = .true.
  Contains
    Procedure, Public :: new_page => flow_type_new_page
    Procedure, Public :: line_printed => flow_type_line_printed
  End Type flow_type

Contains

  Pure Function flow_type_new_page(T) result(new_page)
    Class(flow_type), Intent(In) :: T
    Logical :: new_page

    new_page = Mod(T%lines,T%npage) == 0
  End Function flow_type_new_page

  Subroutine flow_type_line_printed(T)
    Class(flow_type), Intent(InOut) :: T

    T%lines = T%lines + 1
  End Subroutine flow_type_line_printed

End Module flow
