!! Author: Alin M Elena
!! Date: 04-07-2020
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
module table
  use, intrinsic :: iso_fortran_env, only: li => int64,&
                                           wp => real64,&
                                           ou => output_unit

  implicit none

  private

  interface print_row
    module procedure print_row_i
    module procedure print_row_ir
  end interface

  public :: print_header
  public :: print_row
  public :: hline

contains

  pure character(len=250) function hline(w)
    integer, intent(in)    :: w(:)

    character(len=250) :: line
    integer            :: i

    write (line, '(*(a,a),a)') ("+", repeat("-", w(i)), i=1, size(w)), "+"
    hline = line

  end function hline

  subroutine print_header(labels, w, uo)
    character(len=*), intent(in)    :: labels(:)
    integer, intent(in)             :: w(:)
    integer, intent(in), optional    :: uo

    character(len=250) :: fmta, line
    integer            :: i
    integer :: u

    if (present(uo)) then
      u = uo
    else
      u = ou
    endif

    write (fmta, '(a,*(a,i0,a),a)') '(', ("a1,a", w(i), ",", i=1, size(w)), 'a1)'

    line = hline(w)
    write (u, '(a)') trim(line)
    write (u, fmt=trim(fmta)) ("|", trim(labels(i)), i=1, size(labels)), "|"
    write (u, '(a)') trim(line)

  end subroutine print_header

  subroutine print_row_i(c0, d, w,uo)
    character(len=*), intent(in)    :: c0
    integer(kind=li), intent(in)    :: d(:)
    integer, intent(in)             :: w(:)
    integer, intent(in), optional    :: uo

    character(len=250) :: fmta
    !character(len=250) :: line
    integer            :: i

    integer :: u

    if (present(uo)) then
      u = uo
    else
      u = ou
    endif

    write (fmta, '(a,*(a,i0,a),a)') '(', "a1,a", w(1), ",", ("a1,i", w(i), ",", i=2, size(w)), 'a1)'
    write (u, fmt=trim(fmta)) "|", trim(c0), ("|", d(i), i=1, size(d)), "|"
 !   line = hline(w)
 !   write (u, '(a)') trim(line)
  end subroutine print_row_i

  subroutine print_row_ir(c0, d, r, w,uo)
    character(len=*), intent(in)    :: c0
    integer(kind=li), intent(in)    :: d(:)
    real(kind=wp), intent(in)       :: r(:)
    integer, intent(in)             :: w(:)
    integer, intent(in), optional    :: uo

    character(len=250) :: fmta
    !character(len=250) :: line
    integer            :: i
    integer :: u

    if (present(uo)) then
      u = uo
    else
      u = ou
    endif

    write (fmta, '(a,*(a,i0,a),a)') '(', "a1,a", w(1), ",", &
      ("a1,i", w(i + 1), ",", i=1, size(d)), &
      ("a1,f", w(i), ".2,", i=size(d) + 1, size(w)), &
      'a1)'
    write (u, fmt=trim(fmta)) "|", trim(c0), &
      ("|", d(i), i=1, size(d)), &
      ("|", r(i), i=1, size(r)), &
      "|"
!    line = hline(w)
!    write (u, '(a)') trim(line)
  end subroutine print_row_ir
end module table
