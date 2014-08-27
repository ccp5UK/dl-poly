Module pmf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global PMF constraints' variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntpmf  = 0

  Real( Kind = wp ),              Save :: prmpmf = 0.0_wp
  Real( Kind = wp ),              Save :: passpmq(1:5) = (/ & ! QUENCHING per call
                                          0.0_wp         ,  & ! cycles counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average cycles
                                          999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
                                          0.0_wp /)           ! maximum cycles
  Real( Kind = wp ),              Save :: passpmf(1:5,1:2,1:2) = Reshape( (/ & ! dim::1-shake, dim:1:-per-call
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::1-shake, dim:2:-per-tst
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:1:-per-call
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:2:-per-tst
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp /) , (/5,2,2/) )


  Integer,           Allocatable, Save :: numpmf(:),pmffrz(:)
  Integer,           Allocatable, Save :: lstpmf(:,:),listpmf(:,:,:),legpmf(:,:)

  Real( Kind = wp ), Allocatable, Save :: pmfwgt(:,:),pmfwg1(:,:)

  Public :: allocate_pmf_arrays , deallocate_pmf_arrays

Contains

  Subroutine allocate_pmf_arrays()

    Use setup_module, Only : mxtmls,mxtpmf,mxfpmf,mxpmf,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numpmf(1:mxtmls),pmffrz(1:2),                    Stat = fail(1))
    Allocate (lstpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2),          Stat = fail(2))
    Allocate (listpmf(0:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat = fail(3))
    Allocate (legpmf(0:mxfpmf,1:mxatdm),                       Stat = fail(4))
    Allocate (pmfwgt(0:Max(mxtpmf(1),mxtpmf(2)),1:2),          Stat = fail(5))
    Allocate (pmfwg1(0:Max(mxtpmf(1),mxtpmf(2)),1:2),          Stat = fail(6))

    If (Any(fail > 0)) Call error(1036)

    numpmf  = 0
    pmffrz  = 0
    lstpmf  = 0
    listpmf = 0
    legpmf  = 0

    pmfwgt = 0.0_wp
    pmfwg1 = 0.0_wp

  End Subroutine allocate_pmf_arrays

  Subroutine deallocate_pmf_arrays()

    Implicit None
    Integer :: fail

    fail = 0

    Deallocate (numpmf,lstpmf, Stat = fail)

    If (fail > 0) Call error(1037)

  End Subroutine deallocate_pmf_arrays

End module pmf_module
