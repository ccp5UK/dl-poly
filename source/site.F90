Module site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global atomic site variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2008
! contrib   - m.a.seaton june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : Wp

  Implicit None

  Integer,                            Save :: ntpmls = 0 , &
                                              ntpatm = 0 , &
                                              ntpshl = 0


  Character( Len = 40 ), Allocatable, Save :: molnam(:)
  Character( Len = 8 ),  Allocatable, Save :: sitnam(:),unqatm(:),unqshl(:)

  Integer,               Allocatable, Save :: nummols(:),numsit(:),numfrz(:)
  Integer,               Allocatable, Save :: typsit(:),frzsit(:),fresit(:)

  Real( Kind = wp ),     Allocatable, Save :: wgtsit(:),chgsit(:),dofsit(:)
  Real( Kind = wp ),     Allocatable, Save :: numtyp(:),numtypnf(:),dens(:)

  Public :: allocate_site_arrays

Contains

  Subroutine allocate_site_arrays()

    Use setup_module, Only : mxtmls,mxsite,mxatyp

    Implicit None

    Integer, Dimension( 1:7 ) :: fail

    fail = 0

    Allocate (molnam(1:mxtmls),                                    Stat = fail(1))
    Allocate (sitnam(1:mxsite),unqatm(1:mxsite),unqshl(1:mxsite),  Stat = fail(2))
    Allocate (nummols(1:mxtmls),numsit(1:mxtmls),numfrz(1:mxtmls), Stat = fail(3))
    Allocate (typsit(1:mxsite),frzsit(1:mxsite),fresit(1:mxsite),  Stat = fail(4))
    Allocate (wgtsit(1:mxsite),chgsit(1:mxsite),dofsit(1:mxsite),  Stat = fail(5))
    Allocate (numtyp(1:mxatyp),numtypnf(1:mxatyp),dens(1:mxatyp),  Stat = fail(6))

    If (Any(fail > 0)) Call error(1026)

    molnam = ' '
    sitnam = ' ' ; unqatm = ' ' ; unqshl = ' '

    nummols = 0 ; numsit = 0 ; numfrz = 0
    typsit  = 0 ; frzsit = 0 ; fresit = 0

    wgtsit = 0.0_wp ; chgsit = 0.0_wp ; dofsit = 0.0_wp
    numtyp = 0.0_wp ; numtypnf = 0.0_wp ; dens = 0.0_wp

  End Subroutine allocate_site_arrays

End Module site
