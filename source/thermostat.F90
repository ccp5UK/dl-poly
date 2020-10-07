!> Thermostat data module
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J. Madge May 2018
!> Amended - i.t.todorov April 2020

Module thermostat

  Use comms,           Only: comms_type,&
                             gmax
  Use errors_warnings, Only: error
  Use kinds,           Only: wi,&
                             wp
  Use particle,        Only: corePart

  Implicit None

  Private

  ! VV stage keys
  !> First stage
  Integer(Kind=wi), Parameter, Public :: VV_FIRST_STAGE = 0
  !> Second stage
  Integer(Kind=wi), Parameter, Public :: VV_SECOND_STAGE = 1

  ! DPD keys
  !> No DPD
  Integer(Kind=wi), Parameter, Public :: DPD_NULL = 0
  !> First order splitting
  Integer(Kind=wi), Parameter, Public :: DPD_FIRST_ORDER = 1
  !> Second order splitting
  Integer(Kind=wi), Parameter, Public :: DPD_SECOND_ORDER = 2

  ! Pseudo thermostat keys
  !> Langevin + direct temperature scaling
  Integer(Kind=wi), Parameter, Public :: PSEUDO_LANGEVIN_DIRECT = 0
  !> Langevin temperature scaling
  Integer(Kind=wi), Parameter, Public :: PSEUDO_LANGEVIN = 1
  !> Gaussian temperature scaling
  Integer(Kind=wi), Parameter, Public :: PSEUDO_GAUSSIAN = 2
  !> Direct temperature scaling
  Integer(Kind=wi), Parameter, Public :: PSEUDO_DIRECT = 3


  ! Thermostat keys
  !> Microcannonical ensemble
  Integer(Kind=wi), Parameter, Public :: ENS_NVE = 0

  !> Cannonical ensemble Evans
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_EVANS = 1
  !> Cannonical ensemble Langevin (stochastic dynamics)
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_LANGEVIN = 10
  !> Cannonical ensemble Anderson
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_ANDERSON = 11
  !> Cannonical ensemble Berendsen
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_BERENDSEN = 12
  !> Cannonical ensemble Nose-Hoover
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_NOSE_HOOVER = 13
  !> Cannonical ensemble gentle stocastic
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_GENTLE = 14
  !> Cannonical ensemble inhomogeneous Langevin (stocastic dynamics)
  Integer(Kind=wi), Parameter, Public :: ENS_NVT_LANGEVIN_INHOMO = 15

  !> Isobaric ensemble isothermal Langevin (stochastic dynamics) (isotropic)
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_LANGEVIN = 20
  !> Isobaric isothermal ensemble Berendsen
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_BERENDSEN = 21
  !> Isobaric isothermal ensemble Nose-Hoover (isotropic) (Melchionna)
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_NOSE_HOOVER = 22
  !> Isobaric isothermal ensemble Martyna-Tuckerman-Klein (isotropic)
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_MTK = 23

  !> Isobaric isothermal ensemble anisotropic Langvein (stochastic dynamics)
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_LANGEVIN_ANISO = 30
  !> Isobaric isothermal ensemble anisotropic Berendsen
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_BERENDSEN_ANISO = 31
  !> Isobaric isothermal ensemble anisotropic Nos-Hoover (Melchionna)
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_NOSE_HOOVER_ANISO = 32
  !> Isobaric isothermal ensemble anisotropic Martyna-Tuckerman-Klein
  Integer(Kind=wi), Parameter, Public :: ENS_NPT_MTK_ANISO = 33

  ! Anisotropic barostat constraint keys
  !> Fully anisotropic, no constraints
  Integer(Kind=wi), Parameter, Public :: CONSTRAINT_NONE = 0
  !> Semi-isotropic, constant normal pressure and surface area
  Integer(Kind=wi), Parameter, Public :: CONSTRAINT_SURFACE_AREA = 1
  !> Semi-isotropic, constant normal pressure and surface tension
  !> (Orthorhombic constraints when thermo%tension = 0)
  Integer(Kind=wi), Parameter, Public :: CONSTRAINT_SURFACE_TENSION = 2
  !> Semi-orthorhombic constrains
  Integer(Kind=wi), Parameter, Public :: CONSTRAINT_SEMI_ORTHORHOMBIC = 3

  !> Type containing thermostat and barostat variables
  Type, Public :: thermostat_type
    Private

    !> Ensemble key
    Integer(Kind=wi), Public           :: ensemble = ENS_NVE
    !> Flag for variable config%cell size e.g. NPT ensembles
    Logical, Public                    :: variable_cell = .false.
    !> Flag for anisotropic pressure
    Logical, Public                    :: anisotropic_pressure = .false.
    !> Simulation temperature
    Real(Kind=wp), Public              :: temp = 0.0_wp
    !> Simulation pressure
    Real(Kind=wp), Public              :: press = 0.0_wp
    !> Simulation stress
    Real(Kind=wp), Public              :: stress(1:9) = 0.0_wp
    !> Average total energy due to equipartition
    Real(Kind=wp), Public              :: sigma = 0.0_wp
    !> Thermostat relaxation time
    Real(Kind=wp), Public              :: tau_t = 0.0_wp
    !> Barostat relxation time
    Real(Kind=wp), Public              :: tau_p = 0.0_wp
    !> Surface tension
    Real(Kind=wp), Public              :: tension = 0.0_wp
    !> Constraint type for anisotropic barostats
    Integer(Kind=wi), Public           :: iso = CONSTRAINT_NONE
    !> Andersen thermostat softness
    Real(Kind=wp), Public              :: soft = 0.0_wp
    !> Langevin switch
    Logical, Public                    :: l_langevin = .false.
    !> Gentle Stochastic dynamics (Langevin) thermostat friction
    Real(Kind=wp), Public              :: gama = 0.0_wp
    !> Stochastic Dynamics (SD Langevin) thermostat friction
    Real(Kind=wp), Public              :: chi = 0.0_wp
    !> Inhomogeneous Stochastic Dynamics (SD Langevin)
    !> thermostat/electron-phonon friction
    Real(Kind=wp), Public              :: chi_ep = 0.5_wp
    !> Inhomogeneous Stochastic Dynamics (SD Langevin)
    !> thermostat/electronic stopping friction
    Real(Kind=wp), Public              :: chi_es = 0.0_wp
    !> Stochastic Dynamics (SD Langevin) barostat friction
    Real(Kind=wp), Public              :: tai = 0.0_wp
    !> Square of cutoff velocity for inhomogeneous Langevin thermostat and ttm
    Real(Kind=wp), Public              :: vel_es2 = 50.0_wp
    !> Instantaneous thermostat friction
    Real(Kind=wp), Public              :: chi_t = 0.0_wp
    !> Instantaneous barostat friction
    Real(Kind=wp), Public              :: chi_p = 0.0_wp
    !> Friction integral for thermostat/barostat
    Real(Kind=wp), Public              :: cint = 0.0_wp
    !> Cell parameter scaling factor for barostats
    Real(Kind=wp), Public              :: eta(1:9) = 0.0_wp
    !> DPD switch
    Integer, Public                    :: key_dpd = DPD_NULL
    !> DPD drag?
    Real(Kind=wp), Allocatable, Public :: gamdpd(:)
    Real(Kind=wp), Allocatable, Public :: sigdpd(:)
    !> Pseudo thermostat switch
    Logical, Public                    :: l_stochastic_boundaries = .false.
    !> Pseudo thermostat type
    Integer, Public                    :: key_pseudo = PSEUDO_LANGEVIN_DIRECT
    !> Pseudo thermostat temperature
    Real(Kind=wp), Public              :: temp_pseudo = 1.0_wp
    !> Pseudo thermostat thickness
    Real(Kind=wp), Public              :: width_pseudo = 2.0_wp
    !> Temperature scaling switch
    Logical, Public                    :: l_tscale = .false.
    !> Temperature scaling frequency
    Integer, Public                    :: freq_tscale = 0
    !> Temperature regaussing switch
    Logical, Public                    :: l_tgaus = .false.
    !> Temperature regaussing frequency
    Integer, Public                    :: freq_tgaus = 0
    !> Zero temperature optimisation switch
    Logical, Public                    :: l_zero = .false.
    !> Zero temperature regaussing frequency
    Integer, Public                    :: freq_zero = 0
    Logical, Public                    :: newjob_0 = .true.
    Logical, Public                    :: newjob_1 = .true.
    Logical, Public                    :: newjob_2 = .true.
    Logical, Public                    :: newjob_sb = .true.
    Logical, Public                    :: newjob_nst_scl_0 = .true.
    Logical, Public                    :: newjob_nst_scl_1 = .true.
    Logical, Public                    :: newjob_npt_scl_0 = .true.
    Logical, Public                    :: newjob_npt_scl_1 = .true.
    Integer, Public                    :: mxiter, mxkit, kit
    Logical, Public                    :: unsafe = .false.
    Real(Kind=wp), Public              :: volm0, elrc0, virlrc0, h_z, cell0(1:9)
    Real(Kind=wp), Public              :: qmass, ceng, pmass, chip0, rf, factor, temp_lang
    Real(Kind=wp), Allocatable, Public :: dens0(:)
    Integer, Public                    :: ntp, stp, rtp
    Real(Kind=wp), Public              :: rcell(1:9), cwx,cwy,cwz, ecwx,ecwy,ecwz, chit_sb = 0.0_wp
    ! q. index arrays and tp. sum arrays
    Integer, Allocatable, Public       :: qn(:), tpn(:)
    Integer, Allocatable, Public       :: qs(:, :), tps(:)
    Integer, Allocatable, Public       :: qr(:), tpr(:)
    Real(Kind=wp), Public              :: fpl(1:9) = 0.0_wp
    Real(Kind=wp), Public, Allocatable :: fxl(:), fyl(:), fzl(:)
    !> variable timestep control
    Logical, Public                    :: lvar = .false.
    Real(Kind=wp), Public              :: tstep = 0.0_wp
    Real(Kind=wp), Public              :: mndis = 0.03_wp, mxdis = 0.10_wp, mxstp = 0.0_wp

  Contains

    Private

    Procedure, Public :: init_dpd => allocate_dpd_arrays
    Final             :: cleanup

  End Type thermostat_type

  Interface adjust_timestep
    Module Procedure adjust_timestep_1
    Module Procedure adjust_timestep_2
  End Interface

  Public :: adjust_timestep

Contains

  Subroutine allocate_dpd_arrays(thermo, max_vdw)

    Class(thermostat_type)          :: thermo
    Integer(Kind=wi), Intent(In   ) :: max_vdw

    Integer :: fail

    If (thermo%key_dpd == DPD_NULL) Return

    fail = 0

    Allocate (thermo%gamdpd(0:max_vdw), thermo%sigdpd(1:max_vdw), stat=fail)

    If (fail > 0) Call error(1081)

    thermo%gamdpd = 0.0_wp
    thermo%sigdpd = 0.0_wp

  End Subroutine allocate_dpd_arrays

  Subroutine cleanup(thermo)

    Type(thermostat_type) :: thermo

    If (Allocated(thermo%gamdpd)) Then
      Deallocate (thermo%gamdpd)
    End If
    If (Allocated(thermo%sigdpd)) Then
      Deallocate (thermo%sigdpd)
    End If
    If (Allocated(thermo%dens0)) Then
      Deallocate (thermo%dens0)
    End If
    If (Allocated(thermo%qn)) Then
      Deallocate (thermo%qn)
    End If
    If (Allocated(thermo%tpn)) Then
      Deallocate (thermo%tpn)
    End If
    If (Allocated(thermo%qs)) Then
      Deallocate (thermo%qs)
    End If
    If (Allocated(thermo%tps)) Then
      Deallocate (thermo%tps)
    End If
    If (Allocated(thermo%qr)) Then
      Deallocate (thermo%qr)
    End If
    If (Allocated(thermo%tpr)) Then
      Deallocate (thermo%tpr)
    End If

  End Subroutine cleanup

  Logical Function adjust_timestep_1(tstep, hstep, rstep, mndis, mxdis, mxstp, natms, parts, &
                                     xxt, yyt, zzt, legshl, message, t, comm)

    ! update maximum distance a particle has travelled
    Real(Kind=wp),    Intent(InOut) :: tstep
    Real(Kind=wp),    Intent(  Out) :: hstep, rstep
    Real(Kind=wp),    Intent(In   ) :: mndis, mxdis, mxstp
    Integer,          Intent(In   ) :: natms
    Type(corePart),   Intent(In   ) :: parts(1:)
    Real(Kind=wp),    Intent(In   ) :: xxt(1:), yyt(1:), zzt(1:)
    Integer,          Intent(In   ) :: legshl(0:, 1:)
    Character(Len=*), Intent(  Out) :: message
    Real(Kind=wp),    Intent(InOut) :: t
    Type(comms_type), Intent(InOut) :: comm

    Integer       :: i
    Logical       :: lv_dn, lv_up
    Real(Kind=wp) :: mxdr

    adjust_timestep_1 = .false.
    lv_up = .false.
    lv_dn = .false.
    mxdr = 0.0_wp
    t = 0.0_wp

    Do i = 1, natms
      If (legshl(0, i) >= 0) &
        mxdr = Max(mxdr, (parts(i)%xxx - xxt(i))**2 + (parts(i)%yyy - yyt(i))**2 + (parts(i)%zzz - zzt(i))**2)
    End Do
    mxdr = Sqrt(mxdr)
    Call gmax(comm, mxdr)

    If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

      ! scale tstep and derivatives

      If (mxdr > mxdis) Then
        lv_up = .true.
        If (lv_dn) Then
          t = Sqrt(4.0_wp / 3.0_wp)
          tstep = 0.75_wp * tstep
          hstep = 0.50_wp * tstep
        Else
          t = Sqrt(2.0_wp)
          tstep = hstep
          hstep = 0.50_wp * tstep
        End If
        Write (message, "('timestep decreased, new timestep is:',3x,1p,e12.4)") tstep
      End If
      If (mxdr < mndis) Then
        lv_dn = .true.
        If (lv_up) Then
          t = Sqrt(2.0_wp / 3.0_wp)
          tstep = 1.50_wp * tstep
          hstep = 0.50_wp * tstep
        Else
          t = Sqrt(0.5_wp)
          hstep = tstep
          tstep = 2.00_wp * tstep
        End If
        If (tstep > mxstp) Then
          t = t * Sqrt(tstep / mxstp)
          tstep = mxstp
          hstep = 0.50_wp * tstep
        End If
        Write (message, "('timestep increased, new timestep is:',3x,1p,e12.4)") tstep
      End If
      rstep = 1.0_wp / tstep

      ! restart vv1

      adjust_timestep_1 = .true.
    End If

  End Function adjust_timestep_1

  Logical Function adjust_timestep_2(tstep, hstep, rstep, qstep, mndis, mxdis, mxstp, natms, parts, &
                                     xxt, yyt, zzt, legshl, message, t, comm)

    ! update maximum distance a particle has travelled
    Real(Kind=wp),    Intent(InOut) :: tstep
    Real(Kind=wp),    Intent(  Out) :: hstep, rstep, qstep
    Real(Kind=wp),    Intent(In   ) :: mndis, mxdis, mxstp
    Integer,          Intent(In   ) :: natms
    Type(corePart),   Intent(In   ) :: parts(1:)
    Real(Kind=wp),    Intent(In   ) :: xxt(1:), yyt(1:), zzt(1:)
    Integer,          Intent(In   ) :: legshl(0:, 1:)
    Character(Len=*), Intent(  Out) :: message
    Real(Kind=wp),    Intent(InOut) :: t
    Type(comms_type), Intent(InOut) :: comm

    Integer       :: i
    Logical       :: lv_dn, lv_up
    Real(Kind=wp) :: mxdr

    adjust_timestep_2 = .false.
    lv_up = .false.
    lv_dn = .false.
    mxdr = 0.0_wp

    Do i = 1, natms
      If (legshl(0, i) >= 0) &
        mxdr = Max(mxdr, (parts(i)%xxx - xxt(i))**2 + (parts(i)%yyy - yyt(i))**2 + (parts(i)%zzz - zzt(i))**2)
    End Do
    mxdr = Sqrt(mxdr)
    Call gmax(comm, mxdr)

    If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

      ! scale tstep and derivatives
      If (mxdr > mxdis) Then
        lv_up = .true.
        If (lv_dn) Then
          t = Sqrt(4.0_wp / 3.0_wp)
          tstep = 0.75_wp * tstep
          hstep = 0.50_wp * tstep
          qstep = 0.50_wp * hstep
        Else
          t = Sqrt(2.0_wp)
          tstep = hstep
          hstep = 0.50_wp * tstep
          qstep = 0.50_wp * hstep
        End If
        Write (message, '(a,1p,e12.4)') &
          'timestep decreased, new timestep is: ', tstep
      End If
      If (mxdr < mndis) Then
        lv_dn = .true.
        If (lv_up) Then
          t = Sqrt(2.0_wp / 3.0_wp)
          tstep = 1.50_wp * tstep
          hstep = 0.50_wp * tstep
          qstep = 0.50_wp * hstep
        Else
          t = Sqrt(0.5_wp)
          qstep = hstep
          hstep = tstep
          tstep = 2.00_wp * tstep
        End If
        If (tstep > mxstp) Then
          t = t * Sqrt(tstep / mxstp)
          tstep = mxstp
          hstep = 0.50_wp * tstep
          qstep = 0.50_wp * hstep
        End If
        Write (message, '(a,1p,e12.4)') &
          'timestep increased, new timestep is: ', tstep
      End If
      rstep = 1.0_wp / tstep

      ! restart vv1

      adjust_timestep_2 = .true.
    End If

  End Function adjust_timestep_2

End Module thermostat
