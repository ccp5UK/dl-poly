Module thermostat
  Use kinds, Only : wp,wi
  Implicit None

  Private

  !> Type containing thermostat and barostat variables
  Type, Public :: thermostat_type
    !> Ensemble key
    Integer( Kind = wi ) :: ensemble
    !> Flag for variable cell size e.g. NPT ensembles
    Logical :: variable_cell = .false.
    !> Flag for anisotropic pressure
    Logical :: anisotropic_pressure = .false.

    !> Simulation temperature
    Real( Kind = wp ) :: temp
    !> Simulation pressure
    Real( Kind = wp ) :: press
    !> Simulation stress
    Real( Kind = wp ) :: stress(1:9)
    !> Average total energy due to equipartition
    Real( Kind = wp ) :: sigma


    !> Thermostat relaxation time
    Real( Kind = wp ) :: tau_t
    !> Barostat relxation time
    Real( Kind = wp ) :: tau_p

    !> Surface tensionsionsion
    Real( Kind = wp ) :: tension

    !> Constraint type for anisotropic barostats
    !>
    !> - 0 fully anisotropic
    !> - 1 semi-isotropic barostat to constant normal pressure & surface area
    !> - 2 semi-isotropic barostat to constant normal pressure & surface tension
    !>     or with orthorhombic constraints (thermo%tension=0.0_wp)
    !> - 3 semi-isotropic barostat with semi-orthorhombic constraints
    Integer( Kind = wi ) :: iso

    !> Andersen thermostat softness
    Real( Kind = wp ) :: soft

    !> Langevin switch
    Logical :: l_langevin

    !> Gentle Stochastic dynamics (Langevin) thermostat friction
    Real( Kind = wp ) :: gama

    !> Stochastic Dynamics (SD Langevin) thermostat friction
    Real( Kind = wp ) :: chi
    !> Inhomogeneous Stochastic Dynamics (SD Langevin)
    !> thermostat/electron-phonon friction
    Real( Kind = wp ) :: chi_ep
    !> Inhomogeneous Stochastic Dynamics (SD Langevin)
    !> thermostat/electronic stopping friction
    Real( Kind = wp ) :: chi_es
    !> Stochastic Dynamics (SD Langevin) barostat friction
    Real( Kind = wp ) :: tai
    !> Square of cutoff velocity for inhomogeneous Langevin thermostat and ttm
    Real( Kind = wp ) :: vel_es2

    !> Instantaneous thermostat friction
    Real( Kind = wp ) :: chi_t
    !> Instantaneous barostat friction
    Real( Kind = wp ) :: chi_p
    !> Friction integral for thermostat/barostat
    Real( Kind = wp ) :: cint
    !> Cell parameter scaling factor for barostats
    Real( Kind = wp ) :: eta(1:9)

    !> DPD switch
    !>
    !> - 0 no DPD
    !> - 1 first order splitting
    !> - 2 second order splitting
    Integer :: key_dpd
    !> DPD drag?
    Real( Kind = wp ), Allocatable :: gamdpd(:)

    !> Pseudo thermostat switch
    Logical :: l_pseudo
    !> Pseudo thermostat type
    !>
    !> - 0 Langevin + direct temperature scaling
    !> - 1 Langevin temperature scaling
    !> - 2 Gaussian temperature scaling
    !> - 3 direct temperature scaling
    Integer :: key_pseudo
    !> Pseudo thermostat temperature
    Real( Kind = wp ) :: temp_pseudo
    !> Pseudo thermostat thickness
    Real( Kind = wp ) :: width_pseudo

    !> Temperature scaling switch
    Logical :: l_tscale
    !> Temperature scaling frequency
    Integer :: freq_tscale

    !> Temperature regaussing switch
    Logical :: l_tgaus
    !> Temperature regaussing frequency
    Integer :: freq_tgaus

    !> Zero temperature optimisation switch
    Logical :: l_zero
    !> Zero temperature regaussing frequency
    Integer :: freq_zero
  Contains
    Private

    Final :: cleanup
  End Type thermostat_type

  ! Thermostat keys
  !> Microcannonical ensemble
  Integer( Kind = wi ), Parameter, Public :: ENS_NVE = 0

  !> Cannonical ensemble Evans
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_EVANS = 1
  !> Cannonical ensemble Langevin (stochastic dynamics)
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_LANGEVIN = 10
  !> Cannonical ensemble Anderson
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_ANDERSON = 11
  !> Cannonical ensemble Berendsen
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_BERENDSEN = 12
  !> Cannonical ensemble Nosé-Hoover
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_NOSE_HOOVER = 13
  !> Cannonical ensemble gentle stocastic
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_GENTLE = 14
  !> Cannonical ensemble inhomogeneous Langevin (stocastic dynamics)
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_LANGEVIN_INHOMO = 15

  !> Isobaric ensemble isothermal Langevin (stochastic dynamics) (isotropic)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_LANGEVIN = 20
  !> Isobaric isothermal ensemble Berendsen
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_BERENDSEN = 21
  !> Isobaric isothermal ensemble Nosé-Hoover (isotropic) (Melchionna)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_NOSE_HOOVER = 22
  !> Isobaric isothermal ensemble Martyna-Tuckerman-Klein (isotropic)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_MTK = 23

  !> Isobaric isothermal ensemble anisotropic Langvein (stochastic dynamics)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_LANGEVIN_ANISO = 30
  !> Isobaric isothermal ensemble anisotropic Berendsen
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_BERENDSEN_ANISO = 31
  !> Isobaric isothermal ensemble anisotropic Nosé-Hoover (Melchionna)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_NOSE_HOOVER_ANISO = 32
  !> Isobaric isothermal ensemble anistropic Martyna-Tuckerman-Klein
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_MTK_ANISO = 33

Contains

  Subroutine cleanup(thermo)
    Type(thermostat_type) :: thermo

    If (Allocated(thermo%gamdpd)) Then
      Deallocate(thermo%gamdpd)
    End If
  End Subroutine cleanup
End Module thermostat
