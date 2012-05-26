Module netcdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module implementing netCDF wrappers for the sorrted parallel
! I/O only in the CONFIG like files: REFERENCE, HISTORY, REVCON & CFGMIN
!
! copyright - daresbury laboratory
! author    - i.j.bush april 2011
! amended   - i.t.todorov april 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Public :: netcdf_create
  Public :: netcdf_open
  Public :: netcdf_close
  Public :: netcdf_set_def
  Public :: netcdf_get_def
  Public :: netcdf_put_var
  Public :: netcdf_get_var
  Public :: netcdf_get_dim
  Public :: netcdf_get_att
  Public :: netcdf_set_real_precision
  Public :: netcdf_get_real_precision
  Public :: netcdf_get_file_real_precision
  Public :: netcdf_compiled

  Type, Public :: netcdf_desc
     Private
     Integer :: ncid
     Integer :: spatial_id, atom_id, frame_id, cell_spatial_id, cell_angular_id, label_id
     Integer :: spatial_var_id, cell_spatial_var_id, cell_angular_var_id
     Integer :: form_id, imcon_id, time_step_id, time_id, step_id, cell_id, cell_lengths_id, cell_angles_id
     Integer :: coords_id, vels_id, forces_id, name_id, index_id, w_id, q_id, rsd_id
     Integer :: dummy_id
  End Type netcdf_desc

  Private

  Interface netcdf_put_var
     Module Procedure netcdf_put_var_rwp_0d
     Module Procedure netcdf_put_var_rwp_1d
     Module Procedure netcdf_put_var_rwp_2d
     Module Procedure netcdf_put_var_int_0d
     Module Procedure netcdf_put_var_int_1d
     Module Procedure netcdf_put_var_chr_1d
     Module Procedure netcdf_put_var_chr_2d
  End Interface

  Interface netcdf_get_var
     Module Procedure netcdf_get_var_rwp_0d
     Module Procedure netcdf_get_var_rwp_1d
     Module Procedure netcdf_get_var_rwp_2d
     Module Procedure netcdf_get_var_int_0d
     Module Procedure netcdf_get_var_int_1d
     Module Procedure netcdf_get_var_chr_1d
     Module Procedure netcdf_get_var_chr_2d
  End Interface

  Interface netcdf_get_att
     Module Procedure netcdf_get_att_int
     Module Procedure netcdf_get_att_chr
  End Interface

Contains

  Subroutine netcdf_create( name, desc, comm, info )

    Implicit None

    Character( Len = * )  , Intent( In    )            :: name
    Type( netcdf_desc )   , Intent(   Out )            :: desc
    Integer               , Intent( In    ) , Optional :: comm
    Integer               , Intent( In    ) , Optional :: info

    desc%ncid                  = 0

    desc%spatial_id            = 0
    desc%atom_id               = 0
    desc%frame_id              = 0
    desc%cell_spatial_id       = 0
    desc%cell_angular_id       = 0
    desc%label_id              = 0

    desc%spatial_var_id        = 0
    desc%cell_spatial_var_id   = 0
    desc%cell_angular_var_id   = 0

    desc%form_id               = 0
    desc%imcon_id              = 0
    desc%time_step_id          = 0
    desc%time_id               = 0
    desc%step_id               = 0
    desc%cell_id               = 0
    desc%cell_lengths_id       = 0
    desc%cell_angles_id        = 0

    desc%coords_id             = 0
    desc%vels_id               = 0
    desc%forces_id             = 0
    desc%name_id               = 0
    desc%index_id              = 0
    desc%w_id                  = 0
    desc%q_id                  = 0
    desc%rsd_id                = 0

    Call netcdf_compiled()

  End Subroutine netcdf_create

  Subroutine netcdf_open( name, desc, comm, info )

    Implicit None

    Character( Len = * )  , Intent( In    )            :: name
    Type( netcdf_desc )   , Intent(   Out )            :: desc
    Integer               , Intent( In    ) , Optional :: comm
    Integer               , Intent( In    ) , Optional :: info

    desc%ncid                  = 0

    desc%spatial_id            = 0
    desc%atom_id               = 0
    desc%frame_id              = 0
    desc%cell_spatial_id       = 0
    desc%cell_angular_id       = 0
    desc%label_id              = 0

    desc%spatial_var_id        = 0
    desc%cell_spatial_var_id   = 0
    desc%cell_angular_var_id   = 0

    desc%form_id               = 0
    desc%imcon_id              = 0
    desc%time_step_id          = 0
    desc%time_id               = 0
    desc%step_id               = 0
    desc%cell_id               = 0
    desc%cell_lengths_id       = 0
    desc%cell_angles_id        = 0

    desc%coords_id             = 0
    desc%vels_id               = 0
    desc%forces_id             = 0
    desc%name_id               = 0
    desc%index_id              = 0
    desc%w_id                  = 0
    desc%q_id                  = 0
    desc%rsd_id                = 0

    Call netcdf_compiled()

  End Subroutine netcdf_open

  Subroutine netcdf_close( desc )

    Implicit None

    Type( netcdf_desc ), Intent( In    ) :: desc

    Call netcdf_compiled()

  End Subroutine netcdf_close

  Subroutine netcdf_set_def( title, n, desc )

    Implicit None

    Character( Len = * )  , Intent( In    ) :: title
    Integer               , Intent( In    ) :: n
    Type( netcdf_desc )   , Intent( InOut ) :: desc

    Call netcdf_compiled()

  End Subroutine netcdf_set_def

  Subroutine netcdf_put_var_rwp_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Real( Kind = wp )     , Intent( In    )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_rwp_0d

  Subroutine netcdf_put_var_rwp_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_rwp_1d

  Subroutine netcdf_put_var_rwp_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( :, : ), Intent( In    )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_rwp_2d

  Subroutine netcdf_put_var_int_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent( In    )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_int_0d

  Subroutine netcdf_put_var_int_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Integer               , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_int_1d

  Subroutine netcdf_put_var_chr_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_chr_1d

  Subroutine netcdf_put_var_chr_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( :, : ), Intent( In    )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

    Call netcdf_compiled()

  End Subroutine netcdf_put_var_chr_2d

  Subroutine netcdf_get_def( desc, title, n )

    Implicit None

    Type( netcdf_desc )   , Intent( InOut )           :: desc
    Character( Len = * )  , Intent(   Out ), Optional :: title
    Integer               , Intent(   Out ), Optional :: n

    If ( Present( title ) ) title = ' '
    If ( Present( n     ) ) n     = 0

    Call netcdf_compiled()

  End Subroutine netcdf_get_def

  Subroutine netcdf_get_dim( what, desc, val )

    Implicit None

    Character( Len = * )  , Intent( In    ) :: what
    Type( netcdf_desc )   , Intent( In    ) :: desc
    Integer               , Intent(   Out ) :: val

    val = 0

    Call netcdf_compiled()

  End Subroutine netcdf_get_dim

  Subroutine netcdf_get_att_int( what, desc, val )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent(   Out )           :: val

    val = 0

    Call netcdf_compiled()

  End Subroutine netcdf_get_att_int

  Subroutine netcdf_get_att_chr( what, desc, val )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Character( Len = * )  , Intent(   Out )           :: val

    val = ' '

    Call netcdf_compiled()

  End Subroutine netcdf_get_att_chr

  Subroutine netcdf_get_var_rwp_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Real( Kind = wp )     , Intent(   Out )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

    val = 0.0_wp

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_rwp_0d

  Subroutine netcdf_get_var_rwp_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

    val = 0.0_wp

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_rwp_1d

  Subroutine netcdf_get_var_rwp_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( :, : ), Intent(   Out )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

    val = 0.0_wp

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_rwp_2d

  Subroutine netcdf_get_var_int_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent(   Out )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

    val = 0

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_int_0d

  Subroutine netcdf_get_var_int_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Integer               , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

    val = 0

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_int_1d

  Subroutine netcdf_get_var_chr_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

    val = ' '

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_chr_1d

  Subroutine netcdf_get_var_chr_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( :, : ), Intent(   Out )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

    val = ' '

    Call netcdf_compiled()

  End Subroutine netcdf_get_var_chr_2d

  Subroutine netcdf_set_real_precision( p, r, error )

    Integer, Intent( In    ) :: p
    Integer, Intent( In    ) :: r
    Integer, Intent(   Out ) :: error

    error = -1

  End Subroutine netcdf_set_real_precision

  Subroutine netcdf_get_real_precision( p, r, error )

    Integer, Intent(   Out ) :: p
    Integer, Intent(   Out ) :: r
    Integer, Intent(   Out ) :: error

    p     = -1
    r     = -1
    error = -1

  End Subroutine netcdf_get_real_precision

  Subroutine netcdf_get_file_real_precision( desc, p, r, error )

    Implicit None

    Type( netcdf_desc )   , Intent( In    ) :: desc
    Integer               , Intent(   Out ) :: p
    Integer               , Intent(   Out ) :: r
    Integer               , Intent(   Out ) :: error

    p = -1
    r = -1
    error = -1

  End Subroutine netcdf_get_file_real_precision

  Subroutine check( status )

    Implicit None

    Integer, Intent( In    ) :: status

    Call netcdf_compiled()

  End Subroutine check

  Subroutine netcdf_compiled()

    Implicit None

    Write(Unit=*, Fmt=*) ' *** NOT netCDF COMPILED BUILD !!! ***'
    Call error(0)
    Stop

  End Subroutine netcdf_compiled

End Module netcdf_module
