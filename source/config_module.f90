Module config_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global configuration variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2010
! contrib   - i.j.bush march 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Character( Len = 72 ), Save :: cfgname = ' ' , &
                                 sysname = ' '

  Integer,               Save :: natms = 0 , &
                                 nlast = 0 , &
                                 nfree = 0 , &
                                 imc_n =-1 , &
                                 gtl_b = 0

  Real( Kind = wp ),     Save :: cell(1:9) = 0.0_wp , &
                                 volm      = 0.0_wp , &
                                 sumchg    = 0.0_wp


  Character( Len = 8 ), Allocatable, Save :: atmnam(:)

  Integer,              Allocatable, Save :: lsite(:),ltype(:)
  Integer,              Allocatable, Save :: lfrzn(:),lfree(:)
  Integer,              Allocatable, Save :: lsi(:),lsa(:),ltg(:)
  Integer,              Allocatable, Save :: lexatm(:,:)
  Integer,              Allocatable, Save :: list(:,:)
  Integer,              Allocatable, Save :: lstfre(:),gtl(:)

  Real( Kind = wp ),    Allocatable, Save :: weight(:),chge(:)
  Real( Kind = wp ),    Allocatable, Save :: xxx(:),yyy(:),zzz(:)
  Real( Kind = wp ),    Allocatable, Save :: vxx(:),vyy(:),vzz(:)
  Real( Kind = wp ),    Allocatable, Save :: fxx(:),fyy(:),fzz(:)

  Public :: reallocate, allocate_config_arrays_read, allocate_config_arrays

  Interface reallocate
     Module Procedure reallocate_chr_v
     Module Procedure reallocate_int_v
     Module Procedure reallocate_rwp_v
  End Interface

  Private :: reallocate_chr_v, reallocate_int_v, reallocate_rwp_v

Contains

  Subroutine reallocate_chr_v( delta, a, stat )

    Implicit None

    Integer,                           Intent( In    ) :: delta
    Character( Len = * ), Allocatable, Intent( InOut ) :: a(:)
    Integer,                           Intent(   Out ) :: stat

    Integer :: size_old, size_new, size_crs

    Character( Len = Len( a ) ), Allocatable :: tmp(:)

    stat = 0
    If ( delta == 0 ) Return

    size_old = Size( a )
    size_new = size_old + delta
    size_crs = Min( size_old, size_new )

    Allocate ( tmp( 1:size_old ), Stat = stat )
    If ( stat /= 0 ) Return

    tmp = a

    Deallocate ( a, Stat = stat )
    If ( stat /= 0 ) Return

    Allocate ( a( 1:size_new ), Stat = stat )
    If ( stat /= 0 ) Return

    If ( size_crs > 0 ) Then
       a( 1:size_crs ) = tmp( 1:size_crs )
       If ( delta > 0 ) a( size_crs+1:size_new ) = ' '
    End If

    Deallocate ( tmp, Stat = stat )

  End Subroutine reallocate_chr_v

  Subroutine reallocate_int_v( delta, a, stat )

    Implicit None

    Integer,              Intent( In    ) :: delta
    Integer, Allocatable, Intent( InOut ) :: a(:)
    Integer,              Intent(   Out ) :: stat

    Integer :: size_old, size_new, size_crs

    Integer, Allocatable :: tmp(:)

    stat = 0
    If ( delta == 0 ) Return

    size_old = Size( a )
    size_new = size_old + delta
    size_crs = Min( size_old, size_new )

    Allocate ( tmp( 1:size_old ), Stat = stat )
    If ( stat /= 0 ) Return

    tmp = a

    Deallocate ( a, Stat = stat )
    If ( stat /= 0 ) Return

    Allocate ( a( 1:size_new ), Stat = stat )
    If ( stat /= 0 ) Return

    If ( size_crs > 0 ) Then
       a( 1:size_crs ) = tmp( 1:size_crs )
       If ( delta > 0 ) a( size_crs+1:size_new ) = 0
    End If

    Deallocate ( tmp, Stat = stat )

  End Subroutine reallocate_int_v

  Subroutine reallocate_rwp_v( delta, a, stat )

    Implicit None

    Integer,                        Intent( In    ) :: delta
    Real( Kind = wp ), Allocatable, Intent( InOut ) :: a(:)
    Integer,                        Intent(   Out ) :: stat

    Integer :: size_old, size_new, size_crs

    Real( Kind = wp ), Allocatable :: tmp(:)

    stat = 0
    If ( delta == 0 ) Return

    size_old = Size( a )
    size_new = size_old + delta
    size_crs = Min( size_old, size_new )

    Allocate ( tmp( 1:size_old ), Stat = stat )
    If ( stat /= 0 ) Return

    tmp = a

    Deallocate ( a, Stat = stat )
    If ( stat /= 0 ) Return

    Allocate ( a( 1:size_new ), Stat = stat )
    If ( stat /= 0 ) Return

    If ( size_crs > 0 ) Then
       a( 1:size_crs ) = tmp( 1:size_crs )
       If ( delta > 0 ) a( size_crs+1:size_new ) = 0.0_wp
    End If

    Deallocate ( tmp, Stat = stat )

  End Subroutine reallocate_rwp_v


  Subroutine allocate_config_arrays_read( isize )

    Implicit None

    Integer, Intent( In    ) :: isize

    Integer :: fail(1:5)

    fail = 0

    Allocate (atmnam(1:isize),                        Stat = fail(1))
    Allocate (lsi(1:isize),lsa(1:isize),ltg(1:isize), Stat = fail(2))
    Allocate (xxx(1:isize),yyy(1:isize),zzz(1:isize), Stat = fail(3))
    Allocate (vxx(1:isize),vyy(1:isize),vzz(1:isize), Stat = fail(4))
    Allocate (fxx(1:isize),fyy(1:isize),fzz(1:isize), Stat = fail(5))

    If (Any(fail > 0)) Call error(1025)

    atmnam = ' '
    lsi = 0 ; lsa = 0 ;ltg = 0

    xxx = 0.0_wp ; yyy = 0.0_wp ; zzz = 0.0_wp
    vxx = 0.0_wp ; vyy = 0.0_wp ; vzz = 0.0_wp
    fxx = 0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp

  End Subroutine allocate_config_arrays_read

  Subroutine allocate_config_arrays()

    Use setup_module, Only : mxatms,mxatdm,mxexcl,mxlist

    Implicit None

    Integer :: fail(1:5),stat(1:13)

    fail = 0

    Allocate (lsite(1:mxatms),ltype(1:mxatms),          Stat = fail(1))
    Allocate (lfrzn(1:mxatms),lfree(1:mxatms),          Stat = fail(2))
    Allocate (lexatm(0:mxexcl,1:mxatdm),                Stat = fail(3))
    Allocate (list(0:mxlist,1:mxatdm),lstfre(1:mxatdm), Stat = fail(4))
    Allocate (weight(1:mxatms),chge(1:mxatms),          Stat = fail(5))

    If (Any(fail > 0)) Call error(1025)

    lsite = 0 ; ltype = 0
    lfrzn = 0 ; lfree = 0

    lexatm = 0
    list = 0 ; lstfre = 0

    weight = 0.0_wp ; chge = 0.0_wp

! Resize the arrays in allocate_config_arrays_read

    stat = 0

    Call reallocate( mxatms - Size( atmnam ), atmnam, stat( 1) )
    Call reallocate( mxatms - Size( lsi    ), lsi,    stat( 2) )
    Call reallocate( mxatms - Size( lsa    ), lsa,    stat( 3) )
    Call reallocate( mxatms - Size( ltg    ), ltg,    stat( 4) )
    Call reallocate( mxatms - Size( xxx    ), xxx,    stat( 5) )
    Call reallocate( mxatms - Size( yyy    ), yyy,    stat( 6) )
    Call reallocate( mxatms - Size( zzz    ), zzz,    stat( 7) )
    Call reallocate( mxatms - Size( vxx    ), vxx,    stat( 8) )
    Call reallocate( mxatms - Size( vyy    ), vyy,    stat( 9) )
    Call reallocate( mxatms - Size( vzz    ), vzz,    stat(10) )
    Call reallocate( mxatms - Size( fxx    ), fxx,    stat(11) )
    Call reallocate( mxatms - Size( fyy    ), fyy,    stat(12) )
    Call reallocate( mxatms - Size( fzz    ), fzz,    stat(13) )

    If ( Any(stat /= 0 )) Call error(1025)

  End Subroutine allocate_config_arrays

End Module config_module
