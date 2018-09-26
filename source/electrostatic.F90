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
    Logical, Public                :: newjob_spole = .true.
    Real( Kind = wp ), Public :: drewd_spole,rdrewd_spole
    Real( Kind = wp ), Dimension( : ), Allocatable, Public:: erc_spole,fer_spole
    Logical,           Public :: newjob_sspme = .true.
    Integer,           Public :: ixb,iyb,izb, ixt,iyt,izt
    Real( Kind = wp ), Public :: ixbm1_r,iybm1_r,izbm1_r, &
      ixtm0_r,iytm0_r,iztm0_r, &
      kmaxa_r,kmaxb_r,kmaxc_r
    !> blocking factors for splines and fft
    Integer,           Public :: block_x,block_y,block_z
    !> B-spline coefficients
    Complex( Kind = wp ), Dimension( : ),   Allocatable, Public :: bscx,bscy,bscz
    !> indexing arrays for x, y & z as used in parallel fft
    Integer,              Dimension( : ),   Allocatable, Public :: index_x,index_y,index_z
    !> temporary workspace for parallel fft
    Real( Kind = wp ),    Dimension( :,:,: ), Allocatable, Public :: qqc_local
    Complex( Kind = wp ), Dimension( :,:,: ), Allocatable, Public :: qqq_local
    Complex( Kind = wp ), Dimension( :,:,: ), Allocatable, Public :: pfft_work
    !> context for parallel fft
    Integer,           Public :: context
    Logical,           Public :: newjob_fspme= .true.
    Real( Kind = wp ), Public :: kmaxa_lr,kmaxb_lr,kmaxc_lr
!> multipole electrostatics 
    Logical,           Public :: newjob_mpoles = .true.
    Real( Kind = wp ), Public :: drewd_mpoles,rdrewd_mpoles
    Real( Kind = wp ), Dimension( : ), Allocatable, Public :: erc_mpoles,fer_mpoles
!> mpoles

    Logical,           Public :: newjob_mpolesd = .true.
    Real( Kind = wp ), Public :: drewd_mpolesd,rdrewd_mpolesd,alp2,co1,co2,co3,co4,co5,exclcoef, &
      twzz,twtwz,fozz
    Real( Kind = wp ), Dimension( : ), Allocatable, Public :: erc_mpolesd,fer_mpolesd
!> mpoles 2
    Logical,           Public :: newjob_fmpoles = .true.
    Integer,           Public :: ixb_mf,iyb_mf,izb_mf, ixt_mf,iyt_mf,izt_mf
    Real( Kind = wp ), Public :: ixbm1_r_mf,iybm1_r_mf,izbm1_r_mf, &
      ixtm0_r_mf,iytm0_r_mf,iztm0_r_mf, &
      kmaxa_r_mf,kmaxb_r_mf,kmaxc_r_mf, &
      d1_mf(0:8,0:8,0:8)
    ! blocking factors for fft
    Integer,           Public :: block_x_mf,block_y_mf,block_z_mf
    ! B-spline coefficients
    Complex( Kind = wp ), Dimension( : ),       Allocatable, Public :: bscx_mf,bscy_mf,bscz_mf
    ! context for parallel fft
    Integer,           Public :: context_mf
    ! indexing arrays for x, y & z as used in parallel fft
    Integer,              Dimension( : ),       Allocatable, Public :: index_x_mf,index_y_mf,index_z_mf
    ! temporary workspace for parallel fft
    Real( Kind = wp ),    Dimension( :,:,: ),   Allocatable, Public :: qqc_local_mf
    Real( Kind = wp ),    Dimension( :,:,:,: ), Allocatable, Public :: qtc_local_mf
    Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Public :: qqq_local_mf
    Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Public :: qt1_local_mf, &
      qt2_local_mf, &
      qt3_local_mf
    Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Public :: pfft_work_mf

    Logical,           Public :: newjob_mfd = .true.
    Integer,           Public :: ixb_mfd,iyb_mfd,izb_mfd, ixt_mfd,iyt_mfd,izt_mfd
    Real( Kind = wp ), Public :: ixbm1_r_mfd,iybm1_r_mfd,izbm1_r_mfd, &
      ixtm0_r_mfd,iytm0_r_mfd,iztm0_r_mfd, &
      kmaxa_r_mfd,kmaxb_r_mfd,kmaxc_r_mfd
    ! blocking factors for fft
    Integer,           Public :: block_x_mfd,block_y_mfd,block_z_mfd
    ! B-spline coefficients
    Complex( Kind = wp ), Dimension( : ),       Allocatable, Public :: bscx_mfd,bscy_mfd,bscz_mfd
    ! context for parallel fft
    Integer,           Public :: context_mfd
    ! indexing arrays for x, y & z as used in parallel fft
    Integer,              Dimension( : ),       Allocatable, Public :: index_x_mfd,index_y_mfd,index_z_mfd
    ! temporary workspace for parallel fft
    Real( Kind = wp ),    Dimension( :,:,: ),   Allocatable, Public :: qqc_local_mfd
    Real( Kind = wp ),    Dimension( :,:,:,: ), Allocatable, Public :: qtc_local_mfd
    Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Public :: qqq_local_mfd
    Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Public :: qt1_local_mfd, &
      qt2_local_mfd, &
      qt3_local_mfd
    Complex( Kind = wp ), Dimension( :,:,: ),   Allocatable, Public :: pfft_work_mfd
!> exluded mforces
    Logical,           Public :: newjob_emf = .true.
    Real( Kind = wp ), Public :: alp2_emf,co1_emf,co2_emf,co3_emf,co4_emf,co5_emf,co6_emf

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
    If (Allocated(T%erc_spole)) Then 
      Deallocate(T%erc_spole)
    End If
    If (Allocated(T%fer_spole)) Then 
      Deallocate(T%fer_spole)
    End If
    If (Allocated(T%bscx)) Then 
      Deallocate(T%bscx)
    End If
    If (Allocated(T%bscy)) Then 
      Deallocate(T%bscy)
    End If
    If (Allocated(T%bscz)) Then 
      Deallocate(T%bscz)
    End If
    If (Allocated(T%index_x)) Then 
      Deallocate(T%index_x)
    End If
    If (Allocated(T%index_y)) Then 
      Deallocate(T%index_y)
    End If
    If (Allocated(T%index_z)) Then 
      Deallocate(T%index_z)
    End If
    If (Allocated(T%qqq_local)) Then 
      Deallocate(T%qqq_local)
    End If
    If (Allocated(T%qqc_local)) Then 
      Deallocate(T%qqc_local)
    End If
    If (Allocated(T%pfft_work)) Then 
      Deallocate(T%pfft_work)
    End If
    If (Allocated(T%erc_mpoles)) Then 
      Deallocate(T%erc_mpoles)
    End If
    If (Allocated(T%fer_mpoles)) Then 
      Deallocate(T%fer_mpoles)
    End If
    If (Allocated(T%erc_mpolesd)) Then 
      Deallocate(T%erc_mpolesd)
    End If
    If (Allocated(T%fer_mpolesd)) Then 
      Deallocate(T%fer_mpolesd)
    End If

  End Subroutine cleanup
End Module electrostatic
