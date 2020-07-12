Module netcdf_wrap

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module implementing netCDF wrappers for the sorrted parallel
  ! I/O only in the CONFIG like files: REFERENCE, HISTORY, REVCON & CFGMIN
  !
  ! copyright - daresbury laboratory
  ! author    - i.j.bush april 2011
  ! amended   - i.t.todorov aplril 2011
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use, Intrinsic :: iso_fortran_env, Only : real32,real64
  Use kinds, Only : wp,wi
  Use errors_warnings, Only : error
#ifdef NETCDF
  Use netcdf, Only : NF90_NETCDF4, NF90_CLOBBER, NF90_WRITE,              &
    NF90_GLOBAL, NF90_UNLIMITED, NF90_INDEPENDENT,       &
    NF90_CHAR, NF90_INT, NF90_DOUBLE, NF90_FLOAT,        &
    nf90_create, nf90_create_par,                        &
    nf90_open, nf90_open_par, nf90_close,                &
    nf90_redef, nf90_def_var, nf90_def_dim, nf90_enddef, &
    nf90_put_att, nf90_put_var,                          &
    nf90_get_att, nf90_get_var, nf90_var_par_access,     &
    nf90_inquire_dimension, nf90_inq_dimid,              &
    nf90_inq_varid, nf90_inquire_variable,               &
    nf90_noerr, nf90_strerror
#endif
  Implicit None

  Private

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

  Type, Public :: netcdf_param
    Private
    !> The printing precision for reals. Affects only entities of dimension "atom" size
    Integer( Kind = wi ) :: pp = wp
    !> The netCDF handle corresponding to the real printing precision
#ifdef NETCDF
    Integer( Kind = wi ) :: ncp = NF90_DOUBLE
#endif
  End Type netcdf_param

  Type, Public :: netcdf_desc
    Private
    Integer( Kind = wi ) :: ncid
    Integer( Kind = wi ) :: spatial_id, atom_id, frame_id, cell_spatial_id, cell_angular_id, label_id
    Integer( Kind = wi ) :: spatial_var_id, cell_spatial_var_id, cell_angular_var_id
    Integer( Kind = wi ) :: form_id, imcon_id, time_step_id, time_id, step_id, cell_id, cell_lengths_id, cell_angles_id
    Integer( Kind = wi ) :: coords_id, vels_id, forces_id, name_id, index_id, w_id, q_id, rsd_id
    Integer( Kind = wi ) :: dummy_id
  End Type netcdf_desc


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

    Character( Len = * )  , Intent( In    )            :: name
    Type( netcdf_desc )   , Intent(   Out )            :: desc
    Integer               , Intent( In    ) , Optional :: comm
    Integer               , Intent( In    ) , Optional :: info
#ifdef NETCDF
    If ( .Not. Present( comm ) ) Then
      Call check( nf90_create( name, NF90_NETCDF4 + NF90_CLOBBER, desc%ncid ) )
    Else
      Call check( nf90_create_par( name, NF90_NETCDF4 + NF90_CLOBBER, comm, info, desc%ncid ) )
    End If
#else
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
#endif
  End Subroutine netcdf_create

  Subroutine netcdf_open( name, desc, comm, info )

    Character( Len = * )  , Intent( In    )            :: name
    Type( netcdf_desc )   , Intent(   Out )            :: desc
    Integer               , Intent( In    ) , Optional :: comm
    Integer               , Intent( In    ) , Optional :: info
#ifdef NETCDF
    If ( .Not. Present( comm ) ) Then
      Call check( nf90_open( name, NF90_WRITE, desc%ncid ) )
    Else
      Call check( nf90_open_par( name, NF90_WRITE, comm, info, desc%ncid ) )
    End If
#else
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
#endif
  End Subroutine netcdf_open

  Subroutine netcdf_close( desc )

    Type( netcdf_desc ), Intent( In    ) :: desc

#ifdef NETCDF
    Call check( nf90_close( desc%ncid ) )
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_close

  Subroutine netcdf_set_def( title, n, param, desc )

    Character( Len = * ), Intent( In    ) :: title
    Integer,              Intent( In    ) :: n
    Type( netcdf_param ), Intent( In    ) :: param
    Type( netcdf_desc ),  Intent( InOut ) :: desc

#ifdef NETCDF
    Character( Len = 10 ) :: earth_time
    Character( Len = 10 ) :: earth_date

    Integer :: rubbish

    ! Force into definition mode - i.e. ignore the error that this routine returns if already
    ! in definition mode
    rubbish = nf90_redef( desc%ncid )

    earth_date = Repeat( ' ', Len( earth_date ) )
    earth_time = Repeat( ' ', Len( earth_time ) )

    Call date_and_time( earth_date, earth_time )

    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "title"            , Trim( title )      ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "application"      , "DL_POLY"          ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "program"          , "DL_POLY"          ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "programVersion"   , "4.0"              ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "Conventions"      , "DL_POLY"          ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "ConventionVersion", "0.0"              ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "CreationDate"     , Trim( earth_date ) ) )
    Call check( nf90_put_att( desc%ncid, NF90_GLOBAL, "CreationTime"     , Trim( earth_time ) ) )

    Call check( nf90_def_dim( desc%ncid, "spatial", 3             , desc%spatial_id ) )
    Call check( nf90_def_dim( desc%ncid, "dummy"  , 1             , desc%dummy_id ) )
    Call check( nf90_def_dim( desc%ncid, "atom"   , n             , desc%atom_id    ) )
    Call check( nf90_def_dim( desc%ncid, "frame"  , NF90_UNLIMITED, desc%frame_id   ) )
    Call check( nf90_def_dim( desc%ncid, "label"  , 8             , desc%label_id ) )

    Call check( nf90_def_dim( desc%ncid, "cell_spatial", 3, desc%cell_spatial_id ) )
    Call check( nf90_def_dim( desc%ncid, "cell_angular", 3, desc%cell_angular_id ) )

    Call check( nf90_def_var( desc%ncid, "spatial"     , &
      NF90_CHAR, (/ desc%spatial_id /)               , desc%spatial_var_id      ) )
    Call check( nf90_var_par_access( desc%ncid, desc%spatial_var_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "cell_spatial", &
      NF90_CHAR, (/ desc%cell_spatial_id /)          , desc%cell_spatial_var_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_spatial_var_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "cell_angular", &
      NF90_CHAR, (/ desc%label_id, desc%cell_angular_id /), desc%cell_angular_var_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_angular_var_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "datalevel", NF90_INT, (/ desc%frame_id /), varid = desc%form_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%form_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "imageconvention", NF90_INT, (/ desc%frame_id /), varid = desc%imcon_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%imcon_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "timestep", NF90_DOUBLE, (/ desc%frame_id /), varid = desc%time_step_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%time_step_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%time_step_id, "units", "picosecond" ) )

    Call check( nf90_def_var( desc%ncid, "time", NF90_DOUBLE, (/ desc%frame_id /), desc%time_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%time_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%time_id, "units", "picosecond" ) )

    Call check( nf90_def_var( desc%ncid, "step", NF90_INT, (/ desc%frame_id /), desc%step_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%step_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "cell", NF90_DOUBLE, &
      (/ desc%spatial_id, desc%cell_spatial_id, desc%frame_id /), desc%cell_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%cell_id, "units", "Angstrom" ) )

    Call check( nf90_def_var( desc%ncid, "cell_lengths", NF90_DOUBLE, &
      (/ desc%cell_spatial_id, desc%frame_id /), desc%cell_lengths_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_lengths_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%cell_lengths_id, "units", "Angstrom" ) )

    Call check( nf90_def_var( desc%ncid, "cell_angles", NF90_DOUBLE, &
      (/ desc%cell_angular_id, desc%frame_id /), desc%cell_angles_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_angles_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%cell_angles_id, "units", "degree" ) )

    Call check( nf90_def_var( desc%ncid, "coordinates", param%ncp, &
      (/ desc%spatial_id, desc%atom_id, desc%frame_id /), desc%coords_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%coords_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%coords_id, "units", "Angstrom" ) )

    Call check( nf90_def_var( desc%ncid, "velocities", param%ncp, &
      (/ desc%spatial_id, desc%atom_id, desc%frame_id /), desc%vels_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%vels_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%vels_id, "units", "Angstrom/picosecond" ) )

    Call check( nf90_def_var( desc%ncid, "forces", param%ncp, &
      (/ desc%spatial_id, desc%atom_id, desc%frame_id /), desc%forces_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%forces_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%forces_id, "units", "Dalton*Angstrom/picosecond^2" ) )

    Call check( nf90_def_var( desc%ncid, "atomnames", NF90_CHAR, &
      (/ desc%label_id, desc%atom_id, desc%dummy_id /), desc%name_id  ) )
    Call check( nf90_var_par_access( desc%ncid, desc%name_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "indices"  , NF90_INT , &
      (/ desc%atom_id, desc%dummy_id /), desc%index_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%index_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "masses", param%ncp, &
      (/ desc%atom_id, desc%dummy_id /), desc%w_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%w_id, NF90_INDEPENDENT ) )

    Call check( nf90_put_att( desc%ncid, desc%w_id, "units", "Dalton" ) )

    Call check( nf90_def_var( desc%ncid, "charges", param%ncp, &
      (/ desc%atom_id, desc%dummy_id /), desc%q_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%q_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%q_id, "units", "atomic charge units" ) )

    Call check( nf90_def_var( desc%ncid, "RSD", param%ncp,    &
      (/ desc%atom_id, desc%frame_id /), desc%rsd_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%rsd_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%rsd_id, "units", "Angstrom" ) )

    Call check( nf90_enddef( desc%ncid ) )

    Call check( nf90_put_var( desc%ncid,      desc%spatial_var_id, "xyz" ) )
    Call check( nf90_put_var( desc%ncid, desc%cell_spatial_var_id, "abc" ) )
    Call check( nf90_put_var( desc%ncid, desc%cell_angular_var_id, (/ "alpha", "beta ", "gamma" /) ) )
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_set_def

  Subroutine netcdf_put_var_rwp_0d( what, desc, val, start, count )

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Real( Kind = wp )     , Intent( In    )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'time' )
      id = desc%time_id
    Case( 'timestep' )
      id = desc%time_step_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      Call check( nf90_put_var( desc%ncid, id, val, start = (/ start /) ) )
    Else
      Call check( nf90_put_var( desc%ncid, id, val ) )
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_rwp_0d

  Subroutine netcdf_put_var_rwp_1d( what, desc, val, start, count )

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'cell_lengths' )
      id = desc%cell_lengths_id
    Case( 'cell_angles' )
      id = desc%cell_angles_id
    Case( 'masses' )
      id = desc%w_id
    Case( 'charges' )
      id = desc%q_id
    Case( 'rsd' )
      id = desc%rsd_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val ) )
      End If
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_rwp_1d

  Subroutine netcdf_put_var_rwp_2d( what, desc, val, start, count )

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( :, : ), Intent( In    )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'cell' )
      id = desc%cell_id
    Case( 'coordinates' )
      id = desc%coords_id
    Case( 'velocities' )
      id = desc%vels_id
    Case( 'forces' )
      id = desc%forces_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val ) )
      End If
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_rwp_2d

  Subroutine netcdf_put_var_int_0d( what, desc, val, start, count )

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent( In    )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'step' )
      id = desc%step_id
    Case( 'datalevel' )
      id = desc%form_id
    Case( 'imageconvention' )
      id = desc%imcon_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      Call check( nf90_put_var( desc%ncid, id, val, start = (/ start /) ) )
    Else
      Call check( nf90_put_var( desc%ncid, id, val ) )
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_int_0d

  Subroutine netcdf_put_var_int_1d( what, desc, val, start, count )

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Integer               , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'indices' )
      id = desc%index_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val ) )
      End If
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_int_1d

  Subroutine netcdf_put_var_chr_1d( what, desc, val, start, count )

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'atomnames' )
      id = desc%name_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val ) )
      End If
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_chr_1d

  Subroutine netcdf_put_var_chr_2d( what, desc, val, start, count )

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( :, : ), Intent( In    )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'atomnames' )
      id = desc%name_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, start = (/ 1, start /), count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val, start = (/ 1, start /) ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_put_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_put_var( desc%ncid, id, val ) )
      End If
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_put_var_chr_2d

  Subroutine netcdf_get_def( desc, title, n )

    Type( netcdf_desc )   , Intent( InOut )           :: desc
    Character( Len = * )  , Intent(   Out ), Optional :: title
    Integer               , Intent(   Out ), Optional :: n

#ifdef NETCDF
    Call check( nf90_inq_dimid( desc%ncid, "spatial"     , desc%spatial_id ) )
    Call check( nf90_inq_dimid( desc%ncid, "atom"        , desc%atom_id    ) )
    Call check( nf90_inq_dimid( desc%ncid, "frame"       , desc%frame_id   ) )
    Call check( nf90_inq_dimid( desc%ncid, "label"       , desc%label_id ) )
    Call check( nf90_inq_dimid( desc%ncid, "cell_spatial", desc%cell_spatial_id ) )
    Call check( nf90_inq_dimid( desc%ncid, "cell_angular", desc%cell_angular_id ) )

    If ( Present( n ) ) Then
      Call netcdf_get_dim( 'atom', desc, n )
    End If

    Call check( nf90_inq_varid( desc%ncid, "spatial"        , desc%spatial_var_id      ) )
    Call check( nf90_var_par_access( desc%ncid, desc%spatial_var_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "cell_spatial"   , desc%cell_spatial_var_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_spatial_var_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "cell_angular"   , desc%cell_angular_var_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_angular_var_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "datalevel"      , desc%form_id             ) )
    Call check( nf90_var_par_access( desc%ncid, desc%form_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "imageconvention", desc%imcon_id            ) )
    Call check( nf90_var_par_access( desc%ncid, desc%imcon_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "timestep"       , desc%time_step_id        ) )
    Call check( nf90_var_par_access( desc%ncid, desc%time_step_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "time"           , desc%time_id             ) )
    Call check( nf90_var_par_access( desc%ncid, desc%time_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "step"           , desc%step_id             ) )
    Call check( nf90_var_par_access( desc%ncid, desc%step_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "cell"           , desc%cell_id             ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "cell_lengths"   , desc%cell_lengths_id     ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_lengths_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "cell_angles"    , desc%cell_angles_id      ) )
    Call check( nf90_var_par_access( desc%ncid, desc%cell_angles_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "coordinates"    , desc%coords_id           ) )
    Call check( nf90_var_par_access( desc%ncid, desc%coords_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "velocities"     , desc%vels_id             ) )
    Call check( nf90_var_par_access( desc%ncid, desc%vels_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "forces"         , desc%forces_id           ) )
    Call check( nf90_var_par_access( desc%ncid, desc%forces_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "atomnames"      , desc%name_id             ) )
    Call check( nf90_var_par_access( desc%ncid, desc%name_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "indices"        , desc%index_id            ) )
    Call check( nf90_var_par_access( desc%ncid, desc%index_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "masses"         , desc%w_id                ) )
    Call check( nf90_var_par_access( desc%ncid, desc%w_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "charges"        , desc%q_id                ) )
    Call check( nf90_var_par_access( desc%ncid, desc%q_id, NF90_INDEPENDENT ) )

    Call check( nf90_inq_varid( desc%ncid, "RSD"            , desc%rsd_id              ) )
    Call check( nf90_var_par_access( desc%ncid, desc%rsd_id, NF90_INDEPENDENT ) )

    If ( Present( title ) ) Then
      Call check( nf90_get_att( desc%ncid, NF90_GLOBAL, 'title'          , title ) )
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_def

  Subroutine netcdf_get_dim( what, desc, val )

    Character( Len = * )  , Intent( In    ) :: what
    Type( netcdf_desc )   , Intent( In    ) :: desc
    Integer               , Intent(   Out ) :: val

#ifdef NETCDF
    Integer :: dim_id

    Select Case( what )
    Case( 'spatial' )
      dim_id= desc%spatial_id
    Case( 'atom' )
      dim_id = desc%atom_id
    Case( 'frame' )
      dim_id= desc%frame_id
    Case( 'label' )
      dim_id= desc%label_id
    Case( 'cell_spatial' )
      dim_id= desc%cell_spatial_id
    Case( 'cell_angular' )
      dim_id= desc%cell_angular_id
    End Select

    Call check( nf90_inquire_dimension( desc%ncid, dim_id, len = val ) )
#else
    val=0
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_dim

  Subroutine netcdf_get_att_int( what, desc, val )

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent(   Out )           :: val

#ifdef NETCDF
    Call check( nf90_get_att( desc%ncid, NF90_GLOBAL, what, val ) )
#else
    val=0
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_att_int

  Subroutine netcdf_get_att_chr( what, desc, val )

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Character( Len = * )  , Intent(   Out )           :: val

#ifdef NETCDF
    Call check( nf90_get_att( desc%ncid, NF90_GLOBAL, what, val ) )
#else
    val="*"
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_att_chr

  Subroutine netcdf_get_var_rwp_0d( what, desc, val, start, count )

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Real( Kind = wp )     , Intent(   Out )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'time' )
      id = desc%time_id
    Case( 'timestep' )
      id = desc%time_step_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      Call check( nf90_get_var( desc%ncid, id, val, start = (/ start /) ) )
    Else
      Call check( nf90_get_var( desc%ncid, id, val ) )
    End If
#else
    val=0.0_wp
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_rwp_0d

  Subroutine netcdf_get_var_rwp_1d( what, desc, val, start, count )

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'cell_lengths' )
      id = desc%cell_lengths_id
    Case( 'cell_angles' )
      id = desc%cell_angles_id
    Case( 'masses' )
      id = desc%w_id
    Case( 'charges' )
      id = desc%q_id
    Case( 'rsd' )
      id = desc%rsd_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val ) )
      End If
    End If
#else
    val=0.0_wp
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_rwp_1d

  Subroutine netcdf_get_var_rwp_2d( what, desc, val, start, count )

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( :, : ), Intent(   Out )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'cell' )
      id = desc%cell_id
    Case( 'coordinates' )
      id = desc%coords_id
    Case( 'velocities' )
      id = desc%vels_id
    Case( 'forces' )
      id = desc%forces_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val ) )
      End If
    End If
#else
    val=0.0_wp
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_rwp_2d

  Subroutine netcdf_get_var_int_0d( what, desc, val, start, count )

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent(   Out )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'step' )
      id = desc%step_id
    Case( 'datalevel' )
      id = desc%form_id
    Case( 'imageconvention' )
      id = desc%imcon_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      Call check( nf90_get_var( desc%ncid, id, val, start = (/ start /) ) )
    Else
      Call check( nf90_get_var( desc%ncid, id, val ) )
    End If
#else
    val=0
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_int_0d

  Subroutine netcdf_get_var_int_1d( what, desc, val, start, count )

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Integer               , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'indices' )
      id = desc%index_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val ) )
      End If
    End If
#else
    val=0
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_int_1d

  Subroutine netcdf_get_var_chr_1d( what, desc, val, start, count )

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'atomnames' )
      id = desc%name_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, start = start, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val, start = start ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val ) )
      End If
    End If
#else
    val="*"
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_chr_1d

  Subroutine netcdf_get_var_chr_2d( what, desc, val, start, count )

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( :, : ), Intent(   Out )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

#ifdef NETCDF
    Integer :: id

    Select Case( what )
    Case( 'atomnames' )
      id = desc%name_id
    Case Default
    End Select

    If ( Present( start ) ) Then
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, start = (/ 1, start /), count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val, start = (/ 1, start /) ) )
      End If
    Else
      If ( Present( count ) ) Then
        Call check( nf90_get_var( desc%ncid, id, val, count = count ) )
      Else
        Call check( nf90_get_var( desc%ncid, id, val ) )
      End If
    End If
#else
    val="*"
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_var_chr_2d

  Subroutine netcdf_compiled()

#ifndef NETCDF
    Character( Len = 256 ) :: message
    Write(message,'(a)') ' *** NOT netCDF COMPILED BUILD !!! ***'
    Call error(0,message,.true.)
#endif

  End Subroutine netcdf_compiled

  Subroutine netcdf_set_real_precision( k, param, error )

    Integer, Intent( In    ) :: k
    Type( netcdf_param ), Intent( InOut ) :: param
    Integer, Intent(   Out ) :: error

    error = 0

    param%pp = k

#ifdef NETCDF
    Select Case( k )
    Case( real32 )
      param%ncp = NF90_FLOAT
    Case( real64 )
      param%ncp = NF90_DOUBLE
    Case Default
      error = k
    End Select

#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_set_real_precision

  Subroutine netcdf_get_real_precision( param, p, r, error )

    Type( netcdf_param ), Intent( In    ) :: param
    Integer, Intent(   Out ) :: p
    Integer, Intent(   Out ) :: r
    Integer, Intent(   Out ) :: error

#ifdef NETCDF
    error = 0

    Select Case( param%pp )
    Case( real32 )
      p = Precision( 1.0_real32 )
      r = Range( 1.0_real32 )
    Case( real64 )
      p = Precision( 1.0_real64 )
      r = Range( 1.0_real64 )
    Case Default
#endif
      p = -1
      r = -1
      error = -1
#ifdef NETCDF
    End Select
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_real_precision

  Subroutine netcdf_get_file_real_precision( desc, p, r, error )

    Type( netcdf_desc )   , Intent( In    ) :: desc
    Integer               , Intent(   Out ) :: p
    Integer               , Intent(   Out ) :: r
    Integer               , Intent(   Out ) :: error

#ifdef NETCDF
    Integer :: file_real

    error = 0

    Call check( nf90_inquire_variable( desc%ncid, desc%coords_id, xtype = file_real ) )

    Select Case( file_real )
    Case( NF90_FLOAT )
      p = Precision( 1.0_real32 )
      r = Range( 1.0_real32 )
    Case( NF90_DOUBLE )
      p = Precision( 1.0_real64 )
      r = Range( 1.0_real64 )
    Case Default
#endif
      p = -1
      r = -1
      error = -1
#ifdef NETCDF
    End Select
#else
    Call netcdf_compiled()
#endif
  End Subroutine netcdf_get_file_real_precision


  Subroutine check( status )
    Integer, Intent( In    ) :: status
#ifdef NETCDF
    Character( Len = 256 ) :: message
    If ( status /= nf90_noerr ) Then
      Write( message, '( a )' ) 'NETCDF error: '//Trim( nf90_strerror( status ) )
      Call error(0,message,.true.)
    End If
#else
    Call netcdf_compiled()
#endif
  End Subroutine check

End Module netcdf_wrap
