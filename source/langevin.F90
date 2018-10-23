Module langevin

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring:
  !           1) Langevin npt, nst and nvt_lfv ensembles switches & arrays
  !           2) gentle stochastic ensembles switch & Gaussian random number
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2014
  ! amended   - i.t.todorov march 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use constants,      Only : boltz
  Use configuration,     Only : configuration_type
  Use core_shell, Only : core_shell_type
  Use ttm,        Only : ttm_type 
  Use ttm_utils,         Only : Gep
  Use numerics, Only : seed_type,box_mueller_saru3
  Use errors_warnings, Only : error
  Use thermostat, Only : thermostat_type
  Implicit None


  Public :: langevin_allocate_arrays, langevin_forces

Contains

  Subroutine langevin_allocate_arrays(thermo,mxatms)
    Type(thermostat_type), Intent(InOut) :: thermo
    Integer(Kind=wi), Intent( In ) :: mxatms

    Integer :: fail

    fail = 0

    Allocate (thermo%fxl(1:mxatms),thermo%fyl(1:mxatms),thermo%fzl(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1041)

    thermo%fxl = 0.0_wp ; thermo%fyl = 0.0_wp ; thermo%fzl = 0.0_wp

  End Subroutine langevin_allocate_arrays

  Subroutine langevin_forces(nstep,temp,tstep,chi,fxr,fyr,fzr,cshell,config,seed,ttm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to generate Langevin random forces consistent
    ! with the target temperature and the Langevin thermostat friction
    !
    ! Note: (1) Random forces do not contribute to the stress and virial
    !           of the system they are accounted by the system pressure.
    !       (2) Random forces do not apply to frozen and massless particles
    !           as well as shells.
    !       (3) Random forces are scaled according to local electronic
    !           temperature (and dynamically determined friction parameter)
    !           if two-temperature model is in use.
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! contrib   - g.khara & m.a.seaton february 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( ttm_type ), Intent( InOut ),Optional :: ttm   
    Integer          , Intent( In    ) :: nstep
    Real( Kind = wp ), Intent( In    ) :: temp,tstep,chi
    Real( Kind = wp ), Intent(   Out ) :: fxr(:),fyr(:),fzr(:)
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( configuration_type ), Intent( InOut ) :: config
    Type(seed_type), Intent(InOut) :: seed
    Integer           :: i,ia,ja,ka,ijk
    Real( Kind = wp ) :: scale,tmp

    If (Present(ttm) .and. ttm%l_ttm .and. nstep>ttm%nstepcpl) Then

      ! Rescale chi for average electronic temperature if using
      ! homogeneous electron-phonon coupling

      Select Case (ttm%gvar)
      Case (0,1)
        ! constant electron-phonon chi parameter and homogeneous
        ! e-p coupling cases
        scale = Sqrt(2.0_wp * chi * boltz / tstep)
        Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0) Then
            Call box_mueller_saru3(seed,config%ltg(i),nstep,fxr(i),fyr(i),fzr(i))
            ia = Floor((config%parts(i)%xxx+ttm%zerocell(1))/ttm%delx) + 1
            ja = Floor((config%parts(i)%yyy+ttm%zerocell(2))/ttm%dely) + 1
            ka = Floor((config%parts(i)%zzz+ttm%zerocell(3))/ttm%delz) + 1
            ijk = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
            tmp = scale*Sqrt(ttm%eltemp(ijk,0,0,0)*config%weight(i))

            fxr(i) = fxr(i)*tmp
            fyr(i) = fyr(i)*tmp
            fzr(i) = fzr(i)*tmp
          Else
            fxr(i) = 0.0_wp
            fyr(i) = 0.0_wp
            fzr(i) = 0.0_wp
          End If
        End Do

      Case (2)
        ! heterogeneous electron-phonon coupling case: calculate individual
        ! chi value for each ionic temperature voxel (ignore input value)
        scale = Sqrt(2.0_wp * boltz / tstep)
        Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0) Then
            Call box_mueller_saru3(seed,config%ltg(i),nstep,fxr(i),fyr(i),fzr(i))
            ia = Floor((config%parts(i)%xxx+ttm%zerocell(1))/ttm%delx) + 1
            ja = Floor((config%parts(i)%yyy+ttm%zerocell(2))/ttm%dely) + 1
            ka = Floor((config%parts(i)%zzz+ttm%zerocell(3))/ttm%delz) + 1
            ijk = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
            tmp = scale*Sqrt(Gep(ttm%eltemp(ijk,0,0,0),ttm)*ttm%eltemp(ijk,0,0,0)*config%weight(i))

            fxr(i) = fxr(i)*tmp
            fyr(i) = fyr(i)*tmp
            fzr(i) = fzr(i)*tmp
          Else
            fxr(i) = 0.0_wp
            fyr(i) = 0.0_wp
            fzr(i) = 0.0_wp
          End If

        End Do

      End Select

    Else

      ! Get scaler to target variance*Sqrt(weight)

      scale = Sqrt(2.0_wp * chi * boltz * temp / tstep)

      ! Make variance = target variance and nullify the rest

      Do i=1,config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0) Then
          Call box_mueller_saru3(seed,config%ltg(i),nstep,fxr(i),fyr(i),fzr(i))

          tmp = scale*Sqrt(config%weight(i))

          fxr(i) = fxr(i)*tmp
          fyr(i) = fyr(i)*tmp
          fzr(i) = fzr(i)*tmp
        Else
          fxr(i) = 0.0_wp
          fyr(i) = 0.0_wp
          fzr(i) = 0.0_wp
        End If
      End Do

    End If

  End Subroutine langevin_forces


End Module langevin
