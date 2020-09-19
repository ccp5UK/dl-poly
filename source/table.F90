!! Author: Alin M Elena
!! Date: 04-07-2020
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
Module table
  Use, Intrinsic :: iso_fortran_env, Only: li => int64,&
                                           ou => output_unit,&
                                           wp => real64

  Implicit None

  Private

  Interface print_row
    Module Procedure print_row_i
    Module Procedure print_row_ir
  End Interface

  Public :: print_header
  Public :: print_row
  Public :: hline

Contains

  pure Character(len=250) Function hline(w)
    Integer, Intent(In   ) :: w(:)

    Character(len=250) :: line
    Integer            :: i

    Write (line, '(*(a,a),a)') ("+", Repeat("-", w(i)), i=1, Size(w)), "+"
    hline = line

  End Function hline

  Subroutine print_header(labels, w, uo)
    Character(len=*),  Intent(In   ) :: labels(:)
    Integer,           Intent(In   ) :: w(:)
    Integer, Optional, Intent(In   ) :: uo

    Character(len=250) :: fmta, line
    Integer            :: i, u

    If (Present(uo)) Then
      u = uo
    Else
      u = ou
    Endif

    Write (fmta, '(a,*(a,i0,a),a)') '(', ("a1,a", w(i), ",", i=1, Size(w)), 'a1)'

    line = hline(w)
    Write (u, '(a)') Trim(line)
    Write (u, fmt=Trim(fmta)) ("|", Trim(labels(i)), i=1, Size(labels)), "|"
    Write (u, '(a)') Trim(line)

  End Subroutine print_header

  Subroutine print_row_i(c0, d, w, uo)
    Character(len=*),  Intent(In   ) :: c0
    Integer(kind=li),  Intent(In   ) :: d(:)
    Integer,           Intent(In   ) :: w(:)
    Integer, Optional, Intent(In   ) :: uo

    Character(len=250) :: fmta
    Integer            :: i, u

!character(len=250) :: line

    If (Present(uo)) Then
      u = uo
    Else
      u = ou
    Endif

    Write (fmta, '(a,*(a,i0,a),a)') '(', "a1,a", w(1), ",", ("a1,i", w(i), ",", i=2, Size(w)), 'a1)'
    Write (u, fmt=Trim(fmta)) "|", Trim(c0), ("|", d(i), i=1, Size(d)), "|"
    !   line = hline(w)
    !   write (u, '(a)') trim(line)
  End Subroutine print_row_i

  Subroutine print_row_ir(c0, d, r, w, uo)
    Character(len=*),  Intent(In   ) :: c0
    Integer(kind=li),  Intent(In   ) :: d(:)
    Real(kind=wp),     Intent(In   ) :: r(:)
    Integer,           Intent(In   ) :: w(:)
    Integer, Optional, Intent(In   ) :: uo

    Character(len=250) :: fmta
    Integer            :: i, u

!character(len=250) :: line

    If (Present(uo)) Then
      u = uo
    Else
      u = ou
    Endif

    Write (fmta, '(a,*(a,i0,a),a)') '(', "a1,a", w(1), ",", &
      ("a1,i", w(i + 1), ",", i=1, Size(d)), &
      ("a1,f", w(i), ".2,", i=Size(d) + 1, Size(w)), &
      'a1)'
    Write (u, fmt=Trim(fmta)) "|", Trim(c0), &
      ("|", d(i), i=1, Size(d)), &
      ("|", r(i), i=1, Size(r)), &
      "|"
!    line = hline(w)
!    write (u, '(a)') trim(line)
  End Subroutine print_row_ir
End Module table
