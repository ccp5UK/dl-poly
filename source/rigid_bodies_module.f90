Module rigid_bodies_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining rigid bodies' (RBs) variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: lshmv_rgd = .false.

  Integer,                        Save :: ntrgd  = 0 , &
                                          ntrgd1 = 0 , &
                                          m_rgd  = 0

! Globalised megrgd & rcut

  Integer,                        Save :: rgdmeg = 0
  Real( Kind = wp ),              Save :: rgdrct = 0.0_wp

  Integer,           Allocatable, Save :: numrgd(:)
  Integer,           Allocatable, Save :: lstrgd(:,:),listrgd(:,:),legrgd(:,:)
  Integer,           Allocatable, Save :: rgdfrz(:,:),rgdind(:,:)
  Integer,           Allocatable, Save :: lishp_rgd(:),lashp_rgd(:)
  Integer,           Allocatable, Save :: indrgd(:,:)

  Real( Kind = wp ), Allocatable, Save :: rgdwgt(:,:),rgdwg1(:,:)
  Real( Kind = wp ), Allocatable, Save :: rgdx(:,:),rgdy(:,:),rgdz(:,:)
  Real( Kind = wp ), Allocatable, Save :: rgdrix(:,:),rgdriy(:,:),rgdriz(:,:)
  Real( Kind = wp ), Allocatable, Save :: rgdaxs(:,:)
  Real( Kind = wp ), Allocatable, Save :: q0(:),q1(:),q2(:),q3(:)
  Real( Kind = wp ), Allocatable, Save :: rgdxxx(:),rgdyyy(:),rgdzzz(:)
  Real( Kind = wp ), Allocatable, Save :: rgdvxx(:),rgdvyy(:),rgdvzz(:)
  Real( Kind = wp ), Allocatable, Save :: rgdoxx(:),rgdoyy(:),rgdozz(:)

  Public :: allocate_rigid_bodies_arrays , deallocate_rigid_bodies_arrays

Contains

  Subroutine allocate_rigid_bodies_arrays()

    Use setup_module, Only : mxtmls,mxtrgd,mxrgd,mxlrgd,mxfrgd,mxlshp,mxproc,mxatdm

    Implicit None

    Integer, Dimension( 1:15 ) :: fail

    fail = 0

    Allocate (numrgd(1:mxtmls),                                                        Stat = fail( 1))
    Allocate (lstrgd(0:mxlrgd,1:mxtrgd),                                               Stat = fail( 2))
    Allocate (listrgd(-1:mxlrgd,1:mxrgd),                                              Stat = fail( 3))
    Allocate (legrgd(0:mxfrgd,1:mxatdm),                                               Stat = fail( 4))
    Allocate (lishp_rgd(1:mxlshp),lashp_rgd(1:mxproc),                                 Stat = fail( 5))
    Allocate (rgdfrz(0:mxlrgd,1:mxtrgd),rgdind(0:mxlrgd,1:mxtrgd),                     Stat = fail( 6))
    Allocate (rgdwgt(0:mxlrgd,1:mxtrgd),rgdwg1(0:mxlrgd,1:mxtrgd),                     Stat = fail( 7))
    Allocate (indrgd(0:mxlrgd,1:mxrgd),                                                Stat = fail( 8))
    Allocate (rgdx(1:mxlrgd,1:mxtrgd),rgdy(1:mxlrgd,1:mxtrgd),rgdz(1:mxlrgd,1:mxtrgd), Stat = fail( 9))
    Allocate (rgdrix(1:2,1:mxtrgd),rgdriy(1:2,1:mxtrgd),rgdriz(1:2,1:mxtrgd),          Stat = fail(10))
    Allocate (rgdaxs(1:9,1:mxtrgd),                                                    Stat = fail(11))
    Allocate (q0(1:mxrgd),q1(1:mxrgd),q2(1:mxrgd),q3(1:mxrgd),                         Stat = fail(12))
    Allocate (rgdxxx(1:mxrgd),rgdyyy(1:mxrgd),rgdzzz(1:mxrgd),                         Stat = fail(13))
    Allocate (rgdvxx(1:mxrgd),rgdvyy(1:mxrgd),rgdvzz(1:mxrgd),                         Stat = fail(14))
    Allocate (rgdoxx(1:mxrgd),rgdoyy(1:mxrgd),rgdozz(1:mxrgd),                         Stat = fail(15))

    If (Any(fail > 0)) Call error(1042)

    numrgd  = 0
    lstrgd  = 0
    listrgd = 0
    legrgd  = 0

    lishp_rgd = 0 ; lashp_rgd = 0

    rgdfrz = 0 ; rgdind = 0 ; indrgd = 0

    rgdwgt = 0.0_wp ; rgdwg1 = 0.0_wp
    rgdx   = 0.0_wp ; rgdy   = 0.0_wp ; rgdz   = 0.0_wp
    rgdrix = 0.0_wp ; rgdriy = 0.0_wp ; rgdriz = 0.0_wp
    rgdaxs = 0.0_wp

    q0 = 0.0_wp ; q1 = 0.0_wp ; q2 = 0.0_wp ; q3 = 0.0_wp

    rgdxxx = 0.0_wp ; rgdyyy = 0.0_wp ; rgdzzz = 0.0_wp
    rgdvxx = 0.0_wp ; rgdvyy = 0.0_wp ; rgdvzz = 0.0_wp
    rgdoxx = 0.0_wp ; rgdoyy = 0.0_wp ; rgdozz = 0.0_wp

  End Subroutine allocate_rigid_bodies_arrays

  Subroutine deallocate_rigid_bodies_arrays()

    Implicit None
    Integer :: fail

    fail = 0

    Deallocate (numrgd,lstrgd, Stat = fail)

    If (fail > 0) Call error(1043)

  End Subroutine deallocate_rigid_bodies_arrays

End module rigid_bodies_module
