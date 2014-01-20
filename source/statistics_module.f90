Module statistics_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global simulation property variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: numacc = 0 , &
                                          numzdn = 0 , &
                                          numrdf = 0 , &
                                          ntprdf = 0 , &
                                          natms0 = 0

  Real( Kind = wp ),              Save :: pass(1:5) = 0.0_wp


  Integer,           Allocatable, Save :: lstrdf(:)

  Real( Kind = wp ), Allocatable, Save :: xin(:),yin(:),zin(:)
  Real( Kind = wp ), Allocatable, Save :: xto(:),yto(:),zto(:),rsd(:)

  Real( Kind = wp ), Allocatable, Save :: rdf(:,:),zdens(:,:)
  Real( Kind = wp ), Allocatable, Save :: stpval(:),stpvl0(:),sumval(:),ssqval(:)
  Real( Kind = wp ), Allocatable, Save :: zumval(:),ravval(:),stkval(:,:)

  Integer,           Allocatable, Save :: found(:), found0(:)
  Integer,           Allocatable, Save :: lsi0(:),lsa0(:),ltg0(:)

  Real( Kind = wp ), Allocatable, Save :: xin0(:),yin0(:),zin0(:)
  Real( Kind = wp ), Allocatable, Save :: xto0(:),yto0(:),zto0(:)

  Real( Kind = wp ), Allocatable, Save :: stpval0(:),stpvl00(:),sumval0(:),ssqval0(:)
  Real( Kind = wp ), Allocatable, Save :: zumval0(:),ravval0(:),stkval0(:,:)

  Public :: allocate_statistics_arrays, allocate_statistics_connect_arrays

Contains

  Subroutine allocate_statistics_arrays()

    Use setup_module, Only : mxatdm,mxatyp,mxrdf,mxgrdf,mxnstk,mxstak

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (lstrdf(1:mxrdf),                                                     Stat = fail(1))
    Allocate (xin(1:mxatdm),yin(1:mxatdm),zin(1:mxatdm),                           Stat = fail(2))
    Allocate (xto(1:mxatdm),yto(1:mxatdm),zto(1:mxatdm),rsd(1:mxatdm),             Stat = fail(3))
    Allocate (rdf(1:mxgrdf,1:mxrdf),zdens(1:mxgrdf,1:mxatyp),                      Stat = fail(4))
    Allocate (stpval(1:mxnstk),stpvl0(1:mxnstk),sumval(1:mxnstk),ssqval(1:mxnstk), Stat = fail(5))
    Allocate (zumval(1:mxnstk),ravval(1:mxnstk),stkval(1:mxstak,1:mxnstk),         Stat = fail(6))

    If (Any(fail > 0)) Call error(1016)

    lstrdf = 0

    xin = 0.0_wp ; yin = 0.0_wp ; zin = 0.0_wp
    xto = 0.0_wp ; yto = 0.0_wp ; zto = 0.0_wp ; rsd = 0.0_wp

    rdf = 0.0_wp ; zdens = 0.0_wp

    stpval = 0.0_wp ; stpvl0 = 0.0_wp ; sumval = 0.0_wp ; ssqval = 0.0_wp
    zumval = 0.0_wp ; ravval = 0.0_wp ; stkval = 0.0_wp

  End Subroutine allocate_statistics_arrays

  Subroutine allocate_statistics_connect_arrays()

    Use comms_module, Only : mxnode
    Use setup_module, Only : mxatdm,mxnstk,mxstak

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    If (mxnode > 1) Then
       fail = 0

       Allocate (found(1:mxatdm),found0(1:mxatdm),                                                Stat = fail(1))
       Allocate (lsi0(1:mxatdm),lsa0(1:mxatdm),ltg0(1:mxatdm),                                    Stat = fail(2))
       Allocate (xin0(1:mxatdm),yin0(1:mxatdm),zin0(1:mxatdm),                                    Stat = fail(3))
       Allocate (xto0(1:mxatdm),yto0(1:mxatdm),zto0(1:mxatdm),                                    Stat = fail(4))
       Allocate (stpval0(1:2*mxatdm),stpvl00(1:2*mxatdm),sumval0(1:2*mxatdm),ssqval0(1:2*mxatdm), Stat = fail(5))
       Allocate (zumval0(1:2*mxatdm),ravval0(1:2*mxatdm),stkval0(1:mxstak,1:2*mxatdm),            Stat = fail(6))

       If (Any(fail > 0)) Call error(1060)
    End If

  End Subroutine allocate_statistics_connect_arrays

  Subroutine deallocate_statistics_connect_arrays()

    Use comms_module, Only : mxnode
    Use setup_module, Only : mxatdm,mxnstk,mxstak

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    If (mxnode > 1) Then
       fail = 0

       Deallocate (found,found0,                    Stat = fail(1))
       Deallocate (lsi0,lsa0,ltg0,                  Stat = fail(2))
       Deallocate (xin0,yin0,zin0,                  Stat = fail(3))
       Deallocate (xto0,yto0,zto0,                  Stat = fail(4))
       Deallocate (stpval0,stpvl00,sumval0,ssqval0, Stat = fail(5))
       Deallocate (zumval0,ravval0,stkval0,         Stat = fail(6))

       If (Any(fail > 0)) Call error(1061)
    End If

  End Subroutine deallocate_statistics_connect_arrays

End Module statistics_module
