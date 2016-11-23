Module rdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring RDF property variables and arrays
! including UPR (umbrella potential restraint) RDF
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ncfrdf = 0 , &
                                          ntprdf = 0

  Integer,                        Save :: ncfupr = 0

  Real( Kind = wp ),              Save :: rupr   = 0.0_wp ! UPR RDF cutoff

  Integer,           Allocatable, Save :: lstrdf(:)

  Real( Kind = wp ), Allocatable, Save :: rdf(:,:),upr(:)

  Public :: allocate_rdf_arrays

Contains

  Subroutine allocate_rdf_arrays()

    Use setup_module, Only :mxrdf,mxgrdf,mxgupr

    Implicit None

    Integer, Dimension( 1:3 ) :: fail

    fail = 0

    Allocate (lstrdf(1:mxrdf),       Stat = fail(1))
    Allocate (rdf(1:mxgrdf,1:mxrdf), Stat = fail(2))
    Allocate (upr(1:mxgupr),         Stat = fail(3))

    If (Any(fail > 0)) Call error(1016)

    lstrdf = 0

    rdf = 0.0_wp ; upr = 0.0_wp

  End Subroutine allocate_rdf_arrays

End Module rdf_module
