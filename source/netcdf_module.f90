Module netcdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module implementing netCDF wrappers for the sorrted parallel
! I/O only in the CONFIG like files: REFERENCE, HISTORY, REVCON & CFGMIN
!
! copyright - daresbury laboratory
! author    - i.j.bush october 2010
! amended   - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90, Only : wp
  Use netcdf   , Only : NF90_NETCDF4, NF90_CLOBBER, NF90_WRITE,              &
                        NF90_GLOBAL, NF90_UNLIMITED, NF90_INDEPENDENT,       &
                        NF90_CHAR, NF90_INT, NF90_DOUBLE, NF90_FLOAT,        &
                        nf90_create, nf90_create_par,                        &
                        nf90_open, nf90_open_par, nf90_close,                &
                        nf90_redef, nf90_def_var, nf90_def_dim, nf90_enddef, &
                        nf90_put_att, nf90_put_var,                          &
                        nf90_get_att, nf90_get_var, nf90_var_par_access,     &
                        nf90_inquire_dimension, nf90_inq_dimid,              &
                        nf90_inq_varid, nf90_inquire_variable

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
  End Type netcdf_desc

  Private

  ! The printing precision for reals. Affects only entities of dimension "atom" size
  Integer :: pp = wp
  ! The netcdf handle corresponding to the real printing precision
  Integer :: ncp = NF90_DOUBLE

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

    If ( .Not. Present( comm ) ) Then
       Call check( nf90_create( name, NF90_NETCDF4 + NF90_CLOBBER, desc%ncid ) )
    Else
       Call check( nf90_create_par( name, NF90_NETCDF4 + NF90_CLOBBER, comm, info, desc%ncid ) )
    End If

  End Subroutine netcdf_create

  Subroutine netcdf_open( name, desc, comm, info )

    Implicit None

    Character( Len = * )  , Intent( In    )            :: name
    Type( netcdf_desc )   , Intent(   Out )            :: desc
    Integer               , Intent( In    ) , Optional :: comm
    Integer               , Intent( In    ) , Optional :: info

    If ( .Not. Present( comm ) ) Then
       Call check( nf90_open( name, NF90_WRITE, desc%ncid ) )
    Else
       Call check( nf90_open_par( name, NF90_WRITE, comm, info, desc%ncid ) )
    End If

  End Subroutine netcdf_open

  Subroutine netcdf_close( desc )

    Implicit None

    Type( netcdf_desc ), Intent( In    ) :: desc

    Call check( nf90_close( desc%ncid ) )

  End Subroutine netcdf_close

  Subroutine netcdf_set_def( title, n, desc )

    Implicit None

    Character( Len = * )  , Intent( In    ) :: title
    Integer               , Intent( In    ) :: n
    Type( netcdf_desc )   , Intent( InOut ) :: desc

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

    Call check( nf90_def_var( desc%ncid, "coordinates", ncp, &
         (/ desc%spatial_id, desc%atom_id, desc%frame_id /), desc%coords_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%coords_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%coords_id, "units", "Angstrom" ) )

    Call check( nf90_def_var( desc%ncid, "velocities", ncp, &
         (/ desc%spatial_id, desc%atom_id, desc%frame_id /), desc%vels_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%vels_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%vels_id, "units", "Angstrom/picosecond" ) )

    Call check( nf90_def_var( desc%ncid, "forces", ncp, &
         (/ desc%spatial_id, desc%atom_id, desc%frame_id /), desc%forces_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%forces_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%forces_id, "units", "Dalton*Angstrom/picosecond^2" ) )

    Call check( nf90_def_var( desc%ncid, "atomnames", NF90_CHAR, &
         (/ desc%label_id, desc%atom_id, desc%frame_id /), desc%name_id  ) )
    Call check( nf90_var_par_access( desc%ncid, desc%name_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "indices"  , NF90_INT , &
         (/ desc%atom_id, desc%frame_id /), desc%index_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%index_id, NF90_INDEPENDENT ) )

    Call check( nf90_def_var( desc%ncid, "masses", ncp, &
         (/ desc%atom_id, desc%frame_id /), desc%w_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%w_id, NF90_INDEPENDENT ) )

    Call check( nf90_put_att( desc%ncid, desc%w_id, "units", "Dalton" ) )

    Call check( nf90_def_var( desc%ncid, "charges", ncp, &
         (/ desc%atom_id, desc%frame_id /), desc%q_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%q_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%q_id, "units", "atomic charge units" ) )

    Call check( nf90_def_var( desc%ncid, "RSD", ncp,    &
         (/ desc%atom_id, desc%frame_id /), desc%rsd_id ) )
    Call check( nf90_var_par_access( desc%ncid, desc%rsd_id, NF90_INDEPENDENT ) )
    Call check( nf90_put_att( desc%ncid, desc%rsd_id, "units", "Angstrom" ) )

    Call check( nf90_enddef( desc%ncid ) )

    Call check( nf90_put_var( desc%ncid,      desc%spatial_var_id, "xyz" ) )
    Call check( nf90_put_var( desc%ncid, desc%cell_spatial_var_id, "abc" ) )
    Call check( nf90_put_var( desc%ncid, desc%cell_angular_var_id, (/ "alpha", "beta ", "gamma" /) ) )

  End Subroutine netcdf_set_def

  Subroutine netcdf_put_var_rwp_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Real( Kind = wp )     , Intent( In    )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_rwp_0d

  Subroutine netcdf_put_var_rwp_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_rwp_1d

  Subroutine netcdf_put_var_rwp_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( :, : ), Intent( In    )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_rwp_2d

  Subroutine netcdf_put_var_int_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent( In    )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_int_0d

  Subroutine netcdf_put_var_int_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Integer               , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_int_1d

  Subroutine netcdf_put_var_chr_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( : ), Intent( In    )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_chr_1d

  Subroutine netcdf_put_var_chr_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( :, : ), Intent( In    )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_put_var_chr_2d

  Subroutine netcdf_get_def( desc, title, n )

    Implicit None

    Type( netcdf_desc )   , Intent( InOut )           :: desc
    Character( Len = * )  , Intent(   Out ), Optional :: title
    Integer               , Intent(   Out ), Optional :: n

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

  End Subroutine netcdf_get_def

  Subroutine netcdf_get_dim( what, desc, val )

    Implicit None

    Character( Len = * )  , Intent( In    ) :: what
    Type( netcdf_desc )   , Intent( In    ) :: desc
    Integer               , Intent(   Out ) :: val

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

  End Subroutine netcdf_get_dim

  Subroutine netcdf_get_att_int( what, desc, val )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent(   Out )           :: val

    Call check( nf90_get_att( desc%ncid, NF90_GLOBAL, what, val ) )

  End Subroutine netcdf_get_att_int

  Subroutine netcdf_get_att_chr( what, desc, val )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Character( Len = * )  , Intent(   Out )           :: val

    Call check( nf90_get_att( desc%ncid, NF90_GLOBAL, what, val ) )

  End Subroutine netcdf_get_att_chr

  Subroutine netcdf_get_var_rwp_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Real( Kind = wp )     , Intent(   Out )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_rwp_0d

  Subroutine netcdf_get_var_rwp_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_rwp_1d

  Subroutine netcdf_get_var_rwp_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Real( Kind = wp )     , Dimension( :, : ), Intent(   Out )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_rwp_2d

  Subroutine netcdf_get_var_int_0d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )  , Intent( In    )           :: what
    Type( netcdf_desc )   , Intent( In    )           :: desc
    Integer               , Intent(   Out )           :: val
    Integer               , Intent( In    ), Optional :: start
    Integer               , Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_int_0d

  Subroutine netcdf_get_var_int_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Integer               , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_int_1d

  Subroutine netcdf_get_var_chr_1d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                  , Intent( In    )           :: what
    Type( netcdf_desc )                   , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( : ), Intent(   Out )           :: val
    Integer               , Dimension( : ), Intent( In    ), Optional :: start
    Integer               , Dimension( : ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_chr_1d

  Subroutine netcdf_get_var_chr_2d( what, desc, val, start, count )

    Implicit None

    Character( Len = * )                     , Intent( In    )           :: what
    Type( netcdf_desc )                      , Intent( In    )           :: desc
    Character( Len = * )  , Dimension( :, : ), Intent(   Out )           :: val
    Integer               , Dimension( :    ), Intent( In    ), Optional :: start
    Integer               , Dimension( :    ), Intent( In    ), Optional :: count

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

  End Subroutine netcdf_get_var_chr_2d

  Subroutine netcdf_compiled()

    Implicit None

  End Subroutine netcdf_compiled

  Subroutine netcdf_set_real_precision( p, r, error )

    Integer, Intent( In    ) :: p
    Integer, Intent( In    ) :: r
    Integer, Intent(   Out ) :: error

    error = 0

    pp = Selected_real_kind( p, r )

    !In the following the use of 1.0d0 is deliberate - do NOT change them to use wp !!!!!!!!!
    !The reasons are that:
    ! 1) The working precision is independent of the printing precision
    ! 2) I am not allowed to say things like p = Precision( 1.0_pp ) because the kind must be an
    !    initialization expression, essentially one known at compile time
    ! 3) So I must check against the available real kinds, and the only definitely supported ones
    !    in FORTRAN are Kind( 1.0 ) and Kind( 1.0d0 )
    Select Case( pp )
    Case( Kind( 1.0 ) )
       ncp = NF90_FLOAT
    Case( Kind( 1.0d0 ) )
       ncp = NF90_DOUBLE
    End Select

    If ( pp < 0 ) Then
       error = pp
    End If

  End Subroutine netcdf_set_real_precision

  Subroutine netcdf_get_real_precision( p, r, error )

    Integer, Intent(   Out ) :: p
    Integer, Intent(   Out ) :: r
    Integer, Intent(   Out ) :: error

    error = 0

    !In the following the use of 1.0d0 is deliberate - do NOT change them to use wp !!!!!!!!!
    !The reasons are that:
    ! 1) The working precision is independent of the printing precision
    ! 2) I am not allowed to say things like p = Precision( 1.0_pp ) because the kind must be an
    !    initialization expression, essentially one known at compile time
    ! 3) So I must check against the available real kinds, and the only definitely supported ones
    !    in FORTRAN are Kind( 1.0 ) and Kind( 1.0d0 )
    Select Case( pp )
    Case( Kind( 1.0 ) )
       p = Precision( 1.0 )
       r = Range( 1.0 )
    Case( Kind( 1.0d0 ) )
       p = Precision( 1.0d0 )
       r = Range( 1.0d0 )
    Case Default
       p = -1
       r = -1
       error = -1
    End Select

  End Subroutine netcdf_get_real_precision

  Subroutine netcdf_get_file_real_precision( desc, p, r, error )

    Implicit None

    Type( netcdf_desc )   , Intent( In    ) :: desc
    Integer               , Intent(   Out ) :: p
    Integer               , Intent(   Out ) :: r
    Integer               , Intent(   Out ) :: error

    Integer :: file_real

    error = 0

    Call check( nf90_inquire_variable( desc%ncid, desc%coords_id, xtype = file_real ) )

    Select Case( file_real )
    Case( NF90_FLOAT )
       p = Precision( 1.0 )
       r = Range( 1.0 )
    Case( NF90_DOUBLE )
       p = Precision( 1.0d0 )
       r = Range( 1.0d0 )
    Case Default
       p = -1
       r = -1
       error = -1
    End Select

  End Subroutine netcdf_get_file_real_precision

  Subroutine check( status )

    Use netcdf, Only : nf90_noerr, nf90_strerror

    Implicit None

    Integer, Intent( In    ) :: status

    If ( status /= nf90_noerr ) Then
       Write( Unit=*, Fmt=* ) 'NETCDF error:'
       Write( Unit=*, Fmt=* )
       Write( Unit=*, Fmt='( a )' ) trim( nf90_strerror( status ) )
       Call error(0)
       Stop
    End If

  End Subroutine check

End Module netcdf_module
