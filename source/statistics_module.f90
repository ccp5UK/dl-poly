Module statistics_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global simulation property variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: numacc = 0 , &
                                          natms0 = 0

  Real( Kind = wp ),              Save :: clin(1:9) = 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: xin(:),yin(:),zin(:)
  Real( Kind = wp ), Allocatable, Save :: xto(:),yto(:),zto(:),rsd(:)

  Real( Kind = wp ), Allocatable, Save :: stpval(:),stpvl0(:),sumval(:),ssqval(:)
  Real( Kind = wp ), Allocatable, Save :: zumval(:),ravval(:),stkval(:,:)

  Integer,           Allocatable, Save :: found(:), found0(:)
  Integer,           Allocatable, Save :: lsi0(:),lsa0(:),ltg0(:)

  Real( Kind = wp ), Allocatable, Save :: xin0(:),yin0(:),zin0(:)
  Real( Kind = wp ), Allocatable, Save :: xto0(:),yto0(:),zto0(:)

  Real( Kind = wp ), Allocatable, Save :: stpval0(:),stpvl00(:),sumval0(:),ssqval0(:)
  Real( Kind = wp ), Allocatable, Save :: zumval0(:),ravval0(:),stkval0(:,:)

  Logical,                        Save :: statis_file_open = .false.

  Public :: allocate_statistics_arrays, allocate_statistics_connect, &
            deallocate_statistics_connect

Contains

  Subroutine allocate_statistics_arrays()

    Use setup_module, Only : mxatdm,mxnstk,mxstak

    Implicit None

    Integer, Dimension( 1:4 ) :: fail

    fail = 0

    Allocate (xin(1:mxatdm),yin(1:mxatdm),zin(1:mxatdm),                           Stat = fail(1))
    Allocate (xto(1:mxatdm),yto(1:mxatdm),zto(1:mxatdm),rsd(1:mxatdm),             Stat = fail(2))
    Allocate (stpval(0:mxnstk),stpvl0(0:mxnstk),sumval(0:mxnstk),ssqval(0:mxnstk), Stat = fail(3))
    Allocate (zumval(0:mxnstk),ravval(0:mxnstk),stkval(1:mxstak,0:mxnstk),         Stat = fail(4))

    If (Any(fail > 0)) Call error(1016)

    xin = 0.0_wp ; yin = 0.0_wp ; zin = 0.0_wp
    xto = 0.0_wp ; yto = 0.0_wp ; zto = 0.0_wp ; rsd = 0.0_wp

    stpval = 0.0_wp ; stpvl0 = 0.0_wp ; sumval = 0.0_wp ; ssqval = 0.0_wp
    zumval = 0.0_wp ; ravval = 0.0_wp ; stkval = 0.0_wp

  End Subroutine allocate_statistics_arrays

  Subroutine allocate_statistics_connect()

    Use setup_module, Only : mxatdm,mxstak

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (found(1:mxatdm),found0(1:mxatdm),                                                Stat = fail(1))
    Allocate (lsi0(1:mxatdm),lsa0(1:mxatdm),ltg0(1:mxatdm),                                    Stat = fail(2))
    Allocate (xin0(1:mxatdm),yin0(1:mxatdm),zin0(1:mxatdm),                                    Stat = fail(3))
    Allocate (xto0(1:mxatdm),yto0(1:mxatdm),zto0(1:mxatdm),                                    Stat = fail(4))
    Allocate (stpval0(1:2*mxatdm),stpvl00(1:2*mxatdm),sumval0(1:2*mxatdm),ssqval0(1:2*mxatdm), Stat = fail(5))
    Allocate (zumval0(1:2*mxatdm),ravval0(1:2*mxatdm),stkval0(1:mxstak,1:2*mxatdm),            Stat = fail(6))

    If (Any(fail > 0)) Call error(1060)

  End Subroutine allocate_statistics_connect

  Subroutine deallocate_statistics_connect()

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Deallocate (found,found0,                    Stat = fail(1))
    Deallocate (lsi0,lsa0,ltg0,                  Stat = fail(2))
    Deallocate (xin0,yin0,zin0,                  Stat = fail(3))
    Deallocate (xto0,yto0,zto0,                  Stat = fail(4))
    Deallocate (stpval0,stpvl00,sumval0,ssqval0, Stat = fail(5))
    Deallocate (zumval0,ravval0,stkval0,         Stat = fail(6))

    If (Any(fail > 0)) Call error(1061)

  End Subroutine deallocate_statistics_connect

End Module statistics_module
