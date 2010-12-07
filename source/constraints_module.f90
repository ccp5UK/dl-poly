Module constraints_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global constraint bonds variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: lshmv_con = .false.

  Integer,                        Save :: ntcons  = 0 , &
                                          ntcons1 = 0 , &
                                          m_con   = 0


  Integer,           Allocatable, Save :: numcon(:)
  Integer,           Allocatable, Save :: lstcon(:,:),listcon(:,:),legcon(:,:)
  Integer,           Allocatable, Save :: lishp_con(:),lashp_con(:)

  Real( Kind = wp ), Allocatable, Save :: prmcon(:)

  Public :: allocate_constraints_arrays , deallocate_constraints_arrays

Contains

  Subroutine allocate_constraints_arrays()

    Use setup_module, Only : mxtmls,mxtcon,mxcons,mxfcon,mxlshp,mxproc,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numcon(1:mxtmls),                        Stat = fail(1))
    Allocate (lstcon(1:2,1:mxtcon),                    Stat = fail(2))
    Allocate (listcon(0:2,1:mxcons),                   Stat = fail(3))
    Allocate (legcon(0:mxfcon,1:mxatdm),               Stat = fail(4))
    Allocate (lishp_con(1:mxlshp),lashp_con(1:mxproc), Stat = fail(5))
    Allocate (prmcon(1:mxtcon),                        Stat = fail(6))

    If (Any(fail > 0)) Call error(1018)

    numcon  = 0
    lstcon  = 0
    listcon = 0
    legcon  = 0

    lishp_con = 0 ; lashp_con = 0

    prmcon  = 0.0_wp

  End Subroutine allocate_constraints_arrays

  Subroutine deallocate_constraints_arrays()

    Implicit None
    Integer :: fail

    fail = 0

    Deallocate (numcon,lstcon, Stat = fail)

    If (fail > 0) Call error(1032)

  End Subroutine deallocate_constraints_arrays

End module constraints_module
