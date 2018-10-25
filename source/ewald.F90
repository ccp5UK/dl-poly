Module ewald

  !!-----------------------------------------------------------------------
  !!
  !! dl_poly_4 module declaring ewald routines and arrays
  !!
  !! copyright - daresbury laboratory
  !! author    - i.t.todorov february 2015
  !! amended   - j.s.wilkins october  2018
  !!
  !!-----------------------------------------------------------------------

  Use bspline,         Only : bspline_type
  Use comms,           Only : comms_type
  Use configuration,   Only : configuration_type
  Use domains,         Only : domains_type
  Use errors_warnings, Only : error
  Use kinds,           Only : wp
  Use kspace,          Only : kspace_type
  Use constants,       Only : twopi
  Use spme,            Only : spme_component
  Implicit None

  Private

  Type, Public :: ewald_type
    !! Base Ewald type containing data relevant to *all* ewald variants
    Private

    !> Ewald being used in this calculation
    Logical, Public :: active = .false.

    !> Ewald is setup?
    Logical, Public :: initialised = .false.

    !> Ewald is performing polynomial VdW parts
    Logical, Public :: vdw = .false.
    
    !> FFT and KSpace info container
    Type( kspace_type ), Public :: kspace

    !> Frequency of writing per-particle data
    Integer, Public :: pp_write_freq = -1

    !> Ewald convergence parameter or Coulomb damping parameter (A^-1)
    Real( Kind = wp ), Public :: alpha

  End Type ewald_type

  Type, extends (ewald_type), public :: ewald_spme_type
    !! Ewald type containing data relevant to SPME style ewald
    
    !> Number of potentials to handle
    Integer,                                Public :: num_pots = 0
    !> SPME function container
    Type( spme_component ), Dimension( : ), Allocatable, Public :: spme_data
    !> Bspline container
    Type( bspline_type ),                                Public :: bspline
    
  End Type ewald_spme_type  

End Module ewald
