Module core_shell_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global core-shell interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: lshmv_shl = .false.

  Integer,                        Save :: ntshl  = 0 , &
                                          ntshl1 = 0 , &
                                          ntshl2 = 0

  Real( Kind = wp ),              Save :: smax = 0.0_wp
  Real( Kind = wp ),              Save :: passshl(1:5) = (/ &
                                          0.0_wp         ,  & ! cycles counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average cycles
                                          999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
                                          0.0_wp /)           ! maximum cycles


  Integer,           Allocatable, Save :: numshl(:)
  Integer,           Allocatable, Save :: lstshl(:,:),listshl(:,:),legshl(:,:)
  Integer,           Allocatable, Save :: lishp_shl(:),lashp_shl(:)

  Real( Kind = wp ), Allocatable, Save :: prmshl(:,:)

  Public :: allocate_core_shell_arrays , deallocate_core_shell_arrays

Contains

  Subroutine allocate_core_shell_arrays()

    Use setup_module, Only : mxtmls,mxtshl,mxshl,mxfshl,mxlshp,mxproc,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numshl(1:mxtmls),                        Stat = fail(1))
    Allocate (lstshl(1:2,1:mxtshl),                    Stat = fail(2))
    Allocate (listshl(0:2,1:mxshl),                    Stat = fail(3))
    Allocate (legshl(0:mxfshl,1:mxatdm),               Stat = fail(4))
    Allocate (lishp_shl(1:mxlshp),lashp_shl(1:mxproc), Stat = fail(5))
    Allocate (prmshl(1:2,1:mxtshl),                    Stat = fail(6))

    If (Any(fail > 0)) Call error(1005)

    numshl  = 0
    lstshl  = 0
    listshl = 0
    legshl  = 0

    lishp_shl = 0 ; lashp_shl = 0

    prmshl  = 0.0_wp

  End Subroutine allocate_core_shell_arrays

  Subroutine deallocate_core_shell_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numshl,lstshl, Stat = fail)

    If (fail > 0) Call error(1030)

  End Subroutine deallocate_core_shell_arrays

End module core_shell_module
