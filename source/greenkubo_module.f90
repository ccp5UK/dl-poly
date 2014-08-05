Module greenkubo_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring arrays for green-kubo relations
! calculations to calculate transport properties
!
! copyright - daresbury laboratory
! author    - m.a.seaton june 2014
! contrib   - i.t.todorov july 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer           :: isvaf   = 1 , & ! VAF sampling frequency in steps
                       nsvaf   = 0 , & ! VAF sample size
                       vaftsts =-1 , & ! VAF timestep start
                       vafsamp = 0     ! VAF simultaneously overlapping samples

  Real( Kind = wp ) :: vafcount = 0.0_wp

  Integer,           Allocatable, Save :: vafstep(:)
  Real( Kind = wp ), Allocatable, Save :: vxi(:,:),vyi(:,:),vzi(:,:)
  Real( Kind = wp ), Allocatable, Save :: vafdata(:,:),vaftime(:),vaf(:,:)

!  Real( Kind = wp ), Allocatable, Save :: stxx(:),stxy(:),stxz(:),styy(:),styz(:),stzz(:)
!  Real( Kind = wp ), Allocatable, Save :: gkpot(:),tcond(:),tctime(:)

  Public :: allocate_greenkubo_arrays

Contains

  Subroutine allocate_greenkubo_arrays()

    Use setup_module, Only : mxatms,mxatyp

    Implicit None

    Integer :: i
    Integer, Dimension( 1:7 ) :: fail

    fail = 0

    Allocate (vafstep(1:vafsamp),                                                           Stat = fail(1))
    Allocate (vxi(1:mxatms,1:vafsamp),vyi(1:mxatms,1:vafsamp),vzi(1:mxatms,1:vafsamp),      Stat = fail(2))
    Allocate (vafdata(0:nsvaf,1:vafsamp*(mxatyp+1)),vaftime(0:nsvaf),vaf(0:nsvaf,1:mxatyp), Stat = fail(3))
!    Allocate (gkpot(1:mxatms),                                                              Stat = fail(4))
!    Allocate (stxx(1:mxatms),stxy(1:mxatms),stxz(1:mxatms),                                 Stat = fail(5))
!    Allocate (styy(1:mxatms),styz(1:mxatms),stzz(1:mxatms),                                 Stat = fail(6))
!    Allocate (tcond(0:mxgkstk),tctime(0:mxgkstk),                                           Stat = fail(7))

    If (Any(fail > 0)) Call error(1080)

    vxi = 0.0_wp ; vyi = 0.0_wp ; vzi = 0.0_wp
    vaf = 0.0_wp ; vaftime = 0.0_wp
    vafdata = 0.0_wp

! setup step counters for samples, starting after equilibration

    Do i=1,vafsamp
      vafstep(i)=(1-i)*isvaf-1
    End Do

!    sxx = 0.0_wp ; sxy = 0.0_wp ; sxz = 0.0_wp ; syy = 0.0_wp; syz = 0.0_wp; szz = 0.0_wp
!    tcond = 0.0_wp ; tctime = 0.0_wp

  End Subroutine allocate_greenkubo_arrays

End Module greenkubo_module
