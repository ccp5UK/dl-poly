!> Module for electrostatic data and routines, common between Ewald and 
!> non-Ewald methods
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J.Madge July 2018
Module electrostatic
  Use kinds, Only : wi,wp
  Implicit None

  Private

  ! Electrostatic potential keys
  !> No electrostatics
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_NULL = 0
  !> Ewald Sum
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_EWALD = 1
  !> Distance dependent dielectric potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_DDDP = 2
  !> Coulomb potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_COULOMB = 3
  !> Force-shifted and damped coulomb potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_COULOMB_FORCE_SHIFT = 4
  !> Reaction field and damped coulomb potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_COULOMB_REACTION_FIELD = 5
  !> Direct space Poisson solver
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_POISSON = 6

  !> Type containing electrostatic potential data
  Type, Public :: electrostatic_type
    Private

    !> Electrostatic potential key
    Integer( Kind = wi ), Public :: key = ELECTROSTATIC_NULL

    !> Ewald convergence parameter or Coulomb damping parameter (A^-1)
    Real( Kind = wp ), Public :: alpha
    !> Relative dielectric constant
    Real( Kind = wp ), Public :: eps

    !> Grid points for ewald exclusion potential arrays
    Integer( Kind = wi ), Public :: ewald_exclusion_grid

    Logical,           Public :: newjob_intra = .true. , damp
    Real( Kind = wp ), Public :: aa     = 0.0_wp , &
      bb     = 0.0_wp , &
      rfld0  = 0.0_wp , &
      rfld1  = 0.0_wp , &
      rfld2  = 0.0_wp

    Logical,           Public :: newjob_fscp = .true. , damp_fscp
    Real( Kind = wp ), Public :: drewd_fscp  = 0.0_wp , &
      rdrewd_fscp = 0.0_wp , &
      aa_fscp     = 0.0_wp , &
      bb_fscp     = 0.0_wp
    Real( Kind = wp ), Dimension( : ), Allocatable, Public :: erc_fscp,fer_fscp
    Logical,           Public :: newjob_rfp = .true. , damp_rfp
    Real( Kind = wp ), Public :: b0_rfp     = 0.0_wp , &
      rfld0_rfp  = 0.0_wp , &
      rfld1_rfp  = 0.0_wp , &
      rfld2_rfp  = 0.0_wp , &
      drewd_rfp  = 0.0_wp , &
      rdrewd_rfp = 0.0_wp , &
      aa_rfp     = 0.0_wp , &
      bb_rfp     = 0.0_wp , &
      rcsq_rfp   = 0.0_wp
    Real( Kind = wp ), Dimension( : ), Allocatable, Public :: erc_rfp,fer_rfp
    Logical,           Public :: newjob_m = .true. , damp_m
    Real( Kind = wp ), Public :: aa_m     = 0.0_wp , &
      bb_m     = 0.0_wp , &
      rfld0_m  = 0.0_wp , &
      rfld1_m  = 0.0_wp , &
      rfld2_m  = 0.0_wp
    Logical,           Public :: newjob_mfscp = .true. , damp_mfscp
    Real( Kind = wp ), Public :: drewd_mfscp  = 0.0_wp , &
      rdrewd_mfscp = 0.0_wp , &
      aa_mfscp     = 0.0_wp , &
      bb_mfscp     = 0.0_wp
    Real( Kind = wp ), Dimension( : ), Allocatable, Public :: erc_mfscp,fer_mfscp
    Logical,           Public :: newjob_mrfp = .true. , damp_mrfp
    Real( Kind = wp ), Public :: drewd_mrfp  = 0.0_wp , &
      rdrewd_mrfp = 0.0_wp , &
      aa_mrfp     = 0.0_wp , &
      bb_mrfp     = 0.0_wp , &
      b0_mrfp     = 0.0_wp , &
      rfld0_mrfp  = 0.0_wp , &
      rfld1_mrfp  = 0.0_wp , &
      rfld2_mrfp  = 0.0_wp
    Real( Kind = wp ), Dimension( : ), Allocatable, Public :: erc_mrfp,fer_mrfp

  Contains
    Private
    Final ::cleanup
  End Type electrostatic_type
Contains 
  Subroutine cleanup(T)
    Type(electrostatic_type) :: T

    If (Allocated(T%erc_fscp)) Then 
      Deallocate(T%erc_fscp)
    End If
    If (Allocated(T%fer_fscp)) Then 
      Deallocate(T%fer_fscp)
    End If
    If (Allocated(T%erc_rfp)) Then 
      Deallocate(T%erc_rfp)
    End If
    If (Allocated(T%fer_rfp)) Then 
      Deallocate(T%fer_rfp)
    End If
    If (Allocated(T%erc_mfscp)) Then 
      Deallocate(T%erc_mfscp)
    End If
    If (Allocated(T%fer_mfscp)) Then 
      Deallocate(T%fer_mfscp)
    End If
    If (Allocated(T%erc_mrfp)) Then 
      Deallocate(T%erc_mrfp)
    End If
    If (Allocated(T%fer_mrfp)) Then 
      Deallocate(T%fer_mrfp)
    End If

  End Subroutine cleanup
End Module electrostatic
