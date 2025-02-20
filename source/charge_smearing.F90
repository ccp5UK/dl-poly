Module charge_smearing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for evaluating charge smearing corrections for 
  ! SPME based Coulomb interactions 
  !
  ! copyright - daresbury laboratory
  ! author    - b.t.speake July 2024
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only: wp, wi

  Implicit None 

  Private

  Integer(Kind=wi), Parameter, Public :: BETA_ORIGINAL = 0
  Integer(Kind=wi), Parameter, Public :: BETA_OVERLAP = 1
  Integer(Kind=wi), Parameter, Public :: BETA_DISTRIBUTION = 2
  
  Type, Public :: smearing_correction 
    Real(Kind=wp) :: energy, force 
  End Type

  Public :: linear_smearing, linear_smearing_pot, linear_smearing_force 
  Public :: slater_exp_smearing_pot, slater_exp_smearing_force, slater_exp_smearing
  Public :: slater_apprx_smearing_pot, slater_apprx_smearing_force, slater_apprx_smearing

  Contains

  !> Evaluates the linear charge smearing term for both the potential and force 
  !> contributions. The associated correction function for comparison with the 
  !> Coloumbic potential is given as, 
  !>
  !> \begin{equation}
  !> f (r_{ij}) =
  !> \begin{cases}
  !> 1 - \frac{52}{35} \frac{r_{ij}}{R} + \frac{4}{5} \left( \frac{r_{ij}}{R} \right)^3 - \frac{2}{5} \left( \frac{r_{ij}}{R} \right)^5 + \frac{2120}{15603} \left( \frac{r_{ij}}{R} \right)^{6.145} & (r_{ij} < R) \\
  !> \frac{36813504}{11468205} \frac{r_{ij}}{R} \left(1 - \frac{r_{ij}}{2R} \right)^6 & (R \le r_{ij} < 2R) \\
  !> 0 & (R \ge 2R).
  !> \end{cases}
  !> \end{equation}
  Pure Type(smearing_correction) Function linear_smearing(rij, R)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate both the potential and force term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on linear charge smearing. Used when the Ewald evaluation method is 
    ! set to direct. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Real(Kind=wp)             :: pot, force 
    Real(Kind=wp)             :: rijR, rijR2, rijR3, rijR5, rijR6

    rijR = rij / R 

    If (rij < R) Then 
      rijR2 = rijR * rijR             ! ^2
      rijR3 = rijR2 * rijR            ! ^3 
      rijR5 = rijR3 * rijR2           ! ^5 
      rijR6 = rijR**(6.145)           ! ^6.145 

      pot = 1.4857142857142858_wp * rijR ! 52 / 35 = 1.4857142857142858
      pot = pot - (0.8_wp * rijR3) 
      pot = pot + (0.4_wp * rijR5) ! rijR3 * rijR2 ?? 
      pot = pot - (0.13587130679997436_wp * rijR6) ! 2120 / 15603 = 0.13587130679997436 

      force = 1.6_wp * rijR3 
      force = force - (1.6_wp * rijR5) 
      force = force + (0.6990578734858681_wp * rijR6) ! 2597 / 3715 = 0.6990578734858681

    Else  
      rijR5 = (1 - 0.5 * rijR)            ! tmp (1 - rij/2R)
      rijR2 = rijR5 * rijR5               ! ^2
      rijR5 = rijR5  * rijR2 * rijR2      ! ^5 
      rijR6 = rijR2 * rijR2 * rijR2       ! ^6 

      pot = 3.2100493494840734 * rijR * rijR6 

      force = 9.63014804845222_wp * rijR * rijR *rijR5  
    
    End If 
    linear_smearing%energy = pot 
    linear_smearing%force = force 
    Return 
  End Function linear_smearing

  !> Evaluates the linear charge smearing correction term for the Coulomb 
  !> potential. 
  Pure Real(Kind=wp) Function linear_smearing_pot(rij, R) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate the potential term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on linear charge smearing. Used when the Ewald evaluation method is 
    ! set to tabulated. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Real(Kind=wp)             :: pot
    Real(Kind=wp)             :: rijR, rijR2, rijR3, rijR5, rijR6 

    rijR = rij / R 

    If (rij < R) Then 
      rijR2 = rijR * rijR             ! ^2
      rijR3 = rijR2 * rijR            ! ^3 
      rijR5 = rijR3 * rijR2           ! ^5 
      rijR6 = rijR**(6.145)           ! ^6.145 

      pot = 1.4857142857142858_wp * rijR ! 52 / 35 = 1.4857142857142858
      pot = pot - (0.8_wp * rijR3) 
      pot = pot + (0.4_wp * rijR5) ! rijR3 * rijR2 ?? 
      pot = pot - (0.13587130679997436_wp * rijR6) ! 2120 / 15603 = 0.13587130679997436 

    Else  
      rijR5 = (1 - 0.5 * rijR)            ! tmp (1 - rij/2R)
      rijR2 = rijR5 * rijR5               ! ^2
      rijR5 = rijR5  * rijR2 * rijR2      ! ^5 
      rijR6 = rijR2 * rijR2 * rijR2       ! ^6 

      pot = 3.2100493494840734 * rijR * rijR6 
    
    End If 
    linear_smearing_pot = pot 
    Return 
  End Function linear_smearing_pot

  !> Evaluates the linear charge smearing correction for the Coulombic force
  Pure Real(Kind=wp) Function linear_smearing_force(rij, R) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate the force term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on linear charge smearing. Used when the Ewald evaluation method is 
    ! set to tabulated. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Real(Kind=wp)             :: force 
    Real(Kind=wp)             :: rijR, rijR2, rijR3, rijR5, rijR6

    rijR = rij / R 

    If (rij < R) Then 
      rijR2 = rijR * rijR             ! ^2
      rijR3 = rijR2 * rijR            ! ^3 
      rijR5 = rijR3 * rijR2           ! ^5 
      rijR6 = rijR**(6.145)           ! ^6.145 

      force = 1.6_wp * rijR3 
      force = force - (1.6_wp * rijR5) 
      force = force + (0.6990578734858681_wp * rijR6) ! 2597 / 3715 = 0.6990578734858681

    Else  
      rijR5 = (1 - 0.5 * rijR)            ! tmp (1 - rij/2R)
      rijR2 = rijR5 * rijR5               ! ^2
      rijR5 = rijR5  * rijR2 * rijR2      ! ^5 
      rijR6 = rijR2 * rijR2 * rijR2       ! ^6 

      force = 9.63014804845222_wp * rijR * rijR *rijR5  
    
    End If 
    linear_smearing_force = force 
    Return 
  End Function linear_smearing_force

  !> Evaluates the Slater-type charge smearing corrections to both the pairwise 
  !> potential and forces. The related correction function to the standard Coulombic 
  !> potential is given by, 
  !> 
  !> \begin{equation}
  !> f (r_{ij}) = \exp \left(-2 \beta r_{ij} \right) \left(1 + \tfrac{11}{8} \beta r_{ij} + \tfrac{3}{4} \beta^2 r_{ij}^2 + \tfrac{1}{6}\beta^3 r_{ij}^3 \right)
  !> \end{equation}
  Pure Type(smearing_correction) Function slater_exp_smearing(rij, R, beta_type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate both the potential and force term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on slater charge smearing (exact). Used when the Ewald evaluation method
    ! is set to direct. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Integer(Kind=wi), Intent(In) :: beta_type
    Real(Kind=wp) :: br, br2 
    
    Select Case(beta_type)
    Case (BETA_ORIGINAL) 
      br= 1.0_wp / R 
    Case(BETA_OVERLAP)
      br = 5.0_wp / (8.0_wp * R)
    Case(BETA_DISTRIBUTION) 
      br = 1.0_wp / (SQRT(2.0_wp) * R)
    Case Default 
      br = 1.0_wp / R 
    End Select 

    br = br * rij 
    br2 = br * br 

    slater_exp_smearing%energy = 1.0_wp + (1.375_wp*br) + (0.75_wp * br2) + ((1.0_wp/6.0_wp) * br2 * br)
    slater_exp_smearing%energy = slater_exp_smearing%energy * EXP(-2.0_wp * br)

    slater_exp_smearing%force = 1.0_wp + (2.0_wp * br) + (2.0_wp * br2) + ((7.0_wp/6.0_wp) * br2 * br) + &
                                 ((1.0_wp/3.0_wp) * br2 * br2)
    slater_exp_smearing%force = slater_exp_smearing%force * EXP(-2.0_wp * br)
  End Function slater_exp_smearing

  !> Evaluates the Slater-type charge smearing correction to the pairwise Coulombic 
  !> potential. 
  Pure Real(Kind=wp) Function slater_exp_smearing_pot(rij, R, beta_type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate the potential term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on slater charge smearing (exect). Used when the Ewald evaluation method is 
    ! set to tabulated. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Integer(Kind=wi), Intent(In) :: beta_type
    Real(Kind=wp) :: br, br2   


    Select Case(beta_type)
    Case (BETA_ORIGINAL) 
      br= 1.0_wp / R 
    Case(BETA_OVERLAP)
      br = 5.0_wp / (8.0_wp * R)
    Case(BETA_DISTRIBUTION) 
      br = 1.0_wp / (SQRT(2.0_wp) * R)
    Case Default 
      br = 1.0_wp / R 
    End Select 

    br = br * rij 
    br2 = br * br 

    slater_exp_smearing_pot = 1.0_wp + (1.375_wp*br) + (0.75_wp * br2) + ((1.0_wp/6.0_wp) * br2 * br)
    slater_exp_smearing_pot = slater_exp_smearing_pot * EXP(-2.0_wp * br)

  End Function slater_exp_smearing_pot

  !> Evaluates the Slater-type charge smearing correction to the pairwise Coulombic 
  !> force. 
  Pure Real(Kind=wp) Function slater_exp_smearing_force(rij, R, beta_type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate the force term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on slater charge smearing (exact). Used when the Ewald evaluation method is 
    ! set to tabulated. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Integer(Kind=wi), Intent(In) :: beta_type 
    Real(Kind=wp) :: br, br2 

    Select Case(beta_type)
    Case (BETA_ORIGINAL) 
      br= 1.0_wp / R 
    Case(BETA_OVERLAP)
      br = 5.0_wp / (8.0_wp * R)
    Case(BETA_DISTRIBUTION) 
      br = 1.0_wp / (SQRT(2.0_wp) * R)
    Case Default 
      br = 1.0_wp / R 
    End Select 

    br = br * rij 
    br2 = br * br 

    slater_exp_smearing_force = 1.0_wp + (2.0_wp * br) + (2.0_wp * br2) + ((7.0_wp/6.0_wp) * br2 * br) + &
                                 ((1.0_wp/3.0_wp) * br2 * br2)
    slater_exp_smearing_force = slater_exp_smearing_force * EXP(-2.0_wp * br)

  End Function slater_exp_smearing_force

  !> Evaluates the approximate Slater-type charge smearing correction to both the 
  !> pairwise potential and force terms. The related correction function to the 
  !> standard Coulombic potential is given by, 
  !> 
  !> \begin{equation} 
  !> f (r_{ij}) = \exp \left(-2 \beta r_{ij} \right) \left(1 + \beta r_{ij} \right)
  !> \end{equation}
  Pure Type(smearing_correction) Function slater_apprx_smearing(rij, R, beta_type) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate both the potential and force term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on slater charge smearing (approx). Used when the Ewald evaluation method
    ! is set to direct. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Integer(Kind=wi), Intent(In) :: beta_type 
    Real(Kind=wp) :: ebr, br

    Select Case(beta_type)
    Case (BETA_ORIGINAL) 
      br= 1.0_wp / R 
    Case(BETA_OVERLAP)
      br = 5.0_wp / (8.0_wp * R)
    Case(BETA_DISTRIBUTION) 
      br = 1.0_wp / (SQRT(2.0_wp) * R)
    Case Default 
      br = 1.0_wp / R 
    End Select 

    br = br * rij 
    ebr = EXP(-2.0_wp * br)
    
    slater_apprx_smearing%energy = ebr * (1.0_wp + br)
    slater_apprx_smearing%force = ebr * (1.0_wp + 2.0_wp * br * (1.0_wp + br))
  End Function slater_apprx_smearing

  !> Evaluates the approximate Slater-type charge smearing correction to the 
  !> pairwise Coulombic potential. 
  Pure Real(Kind=wp) Function slater_apprx_smearing_pot(rij, R, beta_type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate the potential term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on slater charge smearing (apprx). Used when the Ewald evaluation method is 
    ! set to tabulated. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Integer(Kind=wi), Intent(In) :: beta_type 
    Real(Kind=wp) :: br 

    Select Case(beta_type)
    Case (BETA_ORIGINAL) 
      br= 1.0_wp / R 
    Case(BETA_OVERLAP)
      br = 5.0_wp / (8.0_wp * R)
    Case(BETA_DISTRIBUTION) 
      br = 1.0_wp / (SQRT(2.0_wp) * R)
    Case Default 
      br = 1.0_wp / R 
    End Select 

    br = br * rij 

    slater_apprx_smearing_pot = EXP(-2.0_wp * br) * (1.0_wp + br)

  End Function slater_apprx_smearing_pot

  !> Evaluates the approximate Slater-type charge smearing correction to the 
  !> pairwise Coulombic forces. 
  Pure Real(Kind=wp) Function slater_apprx_smearing_force(rij, R, beta_type) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to evaluate the force term corrections 
    ! to the standard real Ewald component of the Coulomb interactions based 
    ! on slater charge smearing (exact). Used when the Ewald evaluation method is 
    ! set to tabulated. 
    !
    ! copyright - daresbury laboratory
    ! author    - b.t.speake July 2024
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In) :: rij, R 
    Integer(Kind=wi), Intent(In) :: beta_type 
    Real(Kind=wp) :: br 

    Select Case(beta_type)
    Case (BETA_ORIGINAL) 
      br= 1.0_wp / R 
    Case(BETA_OVERLAP)
      br = 5.0_wp / (8.0_wp * R)
    Case(BETA_DISTRIBUTION) 
      br = 1.0_wp / (SQRT(2.0_wp) * R)
    Case Default 
      br = 1.0_wp / R 
    End Select 

    br = br * rij 

    slater_apprx_smearing_force = EXP(-2.0_wp * br) * (1.0_wp + 2.0_wp * br * (1.0_wp + br))

  End Function slater_apprx_smearing_force

End Module 