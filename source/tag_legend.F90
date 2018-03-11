Subroutine tag_legend(safe,iatm,nt,legend,mxf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to tag the legend arrays recording the
! intra-like descriptions for each atom in a domain
!
! nt should always be supplied positive except for shells!
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use setup_module

  Implicit None

  Logical,                              Intent( InOut) :: safe
  Integer,                              Intent( In   ) :: iatm,nt,mxf
  Integer, Dimension( 0:mxf, 1:mxatdm), Intent( InOut) :: legend

  Logical :: safe_local
  Integer :: last

! Get current length

  last = Abs(legend(0,iatm))

! Get local safety no array overflow

  safe_local = (last < mxf-1)

! Determine global safety

  safe = safe .and. safe_local

  If (safe_local) Then

! Increase length of links and tag: I, iatm, am linked to one more entity
! (particle) in a unit of this legend type with a local domain number 'nt'

     last = last + 1
     If (.not.(nt < 0 .and. last == 1)) Then ! THE NORMAL CASE
        legend(0,iatm) = last
        legend(last,iatm) = nt
     Else ! This should be hit just the once for shell particles only!!!
        legend(0,iatm) = -1
        legend(last,iatm) = -nt
     End If

  Else

! Collect number of offences

     legend(mxf,iatm) = legend(mxf,iatm) + 1

  End If

End Subroutine tag_legend
