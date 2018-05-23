Module thermostat
  Use kinds, Only : wp
  Implicit None

  Private

  !> Type containing thermostat and barostat variables
  Type, Public :: thermostat_type
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

  End Type thermostat_type
End Module thermostat
