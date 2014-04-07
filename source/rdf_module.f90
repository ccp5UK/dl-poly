Module rdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring RDF property variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ncfrdf = 0 , &
                                          ntprdf = 0

  Integer,           Allocatable, Save :: lstrdf(:)

  Real( Kind = wp ), Allocatable, Save :: rdf(:,:)

  Public :: allocate_rdf_arrays

Contains

  Subroutine allocate_rdf_arrays()

    Use setup_module, Only :mxrdf,mxgrdf

    Implicit None

    Integer, Dimension( 1:2 ) :: fail

    fail = 0

    Allocate (lstrdf(1:mxrdf),       Stat = fail(1))
    Allocate (rdf(1:mxgrdf,1:mxrdf), Stat = fail(2))

    If (Any(fail > 0)) Call error(1016)

    lstrdf = 0

    rdf = 0.0_wp

  End Subroutine allocate_rdf_arrays

End Module rdf_module
