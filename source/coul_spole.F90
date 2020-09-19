Module coul_spole
  Use kinds,           Only : wp
  Use configuration,   Only : configuration_type
  Use constants,       Only : r4pie0,zero_plus,sqrpi
  Use errors_warnings, Only : error, error_alloc, error_dealloc
  Use numerics,        Only : calc_erfc, calc_erfc_deriv
  Use neighbours,      Only : neighbours_type
  Use electrostatic,   Only : electrostatic_type, &
    ELECTROSTATIC_EWALD,ELECTROSTATIC_DDDP, &
    ELECTROSTATIC_COULOMB,ELECTROSTATIC_COULOMB_FORCE_SHIFT, &
    ELECTROSTATIC_COULOMB_REACTION_FIELD
  Use statistics,      Only : calculate_stress, stats_type
  Implicit None

  Private

  Public :: intra_coul,coul_fscp_forces, coul_rfp_forces, coul_cp_forces, coul_dddp_forces

Contains

  Subroutine intra_coul(rcut,chgprd,rrr,rsq,coul,fcoul,safe,electro)

    !!------------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating bond's or 1-4 dihedral
    !! electrostatics: adjusted by a config%weighting factor
    !!
    !! copyright - daresbury laboratory
    !! amended   - i.t.todorov february 2016
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Real( Kind = wp ), Intent( In    ) :: chgprd,rcut,rrr,rsq
    Real( Kind = wp ), Intent(   Out ) :: coul,fcoul
    Logical,           Intent( InOut ) :: safe
    Type( electrostatic_type ), Intent( InOut ) :: electro

    Logical, save :: newjob = .true.

    Real( Kind = wp ) :: exp1,tt,erc,fer,b0

    Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

    If (newjob) Then
      newjob = .false.

      ! Check for damped force-shifted coulombic and reaction field interactions
      ! and set force and potential shifting parameters dependingly

      If (electro%damp) Then

        erc = calc_erfc(electro%damping*rcut) / rcut
        fer = (erc + electro%damping*calc_erfc_deriv(electro%damping*rcut)) / rcut**2

        electro%force_shift  = fer*rcut
        electro%energy_shift  = -(erc + electro%force_shift*rcut)
      Else If (electro%key == ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then
        electro%force_shift =  1.0_wp/rcut**2
        electro%energy_shift = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)
      End If

      ! set reaction field terms for RFC

      If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then
        b0    = 2.0_wp*(electro%eps - 1.0_wp)/(2.0_wp*electro%eps + 1.0_wp)
        electro%reaction_field(0) = b0/rcut**3
        electro%reaction_field(1) = (1.0_wp + 0.5_wp*b0)/rcut
        electro%reaction_field(2) = 0.5_wp*electro%reaction_field(0)
      End If
    End If

    ! initialise defaults for coulombic energy and force contributions

    coul =0.0_wp
    fcoul=0.0_wp

    ! Electrostatics by ewald sum = direct coulombic

    If (Any([ELECTROSTATIC_EWALD,ELECTROSTATIC_COULOMB] == electro%key)) Then

      coul = chgprd/rrr
      fcoul= coul/rsq

      ! distance dependent dielectric

    Else If (electro%key ==  ELECTROSTATIC_DDDP) Then

      coul = chgprd/rsq
      fcoul= 2.0_wp*coul/rsq

      ! force shifted coulombic and reaction field

    Else If (Any([ELECTROSTATIC_COULOMB_FORCE_SHIFT,ELECTROSTATIC_COULOMB_REACTION_FIELD] == electro%key)) Then

      If (electro%damp) Then ! calculate damping contributions
        exp1= Exp(-(electro%damping*rrr)**2)
        tt  = 1.0_wp/(1.0_wp+pp*electro%damping*rrr)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rrr
        fer = (erc + 2.0_wp*(electro%damping/sqrpi)*exp1)/rsq

        coul = chgprd*(erc + electro%force_shift*rrr + electro%energy_shift)
        fcoul= chgprd*(fer - electro%force_shift/rrr)
      End If

      If (electro%key ==  ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then ! force shifted coulombic
        If (.not.electro%damp) Then ! pure
          coul = chgprd*(1.0_wp/rrr + electro%force_shift*rrr+ electro%energy_shift)
          fcoul= chgprd*(1.0_wp/rsq - electro%force_shift)/rrr
        Else                ! damped
          coul = coul
          fcoul= fcoul
        End If
      Else If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then ! reaction field
        If (.not.electro%damp) Then ! pure
          coul = chgprd*(1.0_wp/rrr + electro%reaction_field(2)*rsq - electro%reaction_field(1))
          fcoul= chgprd*(1.0_wp/rsq/rrr - electro%reaction_field(0))
        Else                ! damped
          coul = coul  + chgprd*(electro%reaction_field(2)*rsq - electro%reaction_field(1))
          fcoul= fcoul + chgprd*(-electro%reaction_field(0))
        End If
      End If

    Else

      safe = .false.

    End If
  End Subroutine intra_coul

  Subroutine coul_fscp_forces(iatm,xxt,yyt,zzt,rrt,engcpe,vircpe,stats,neigh, &
    electro,config)

    !!------------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system assuming a force shifted coulombic potential
    !!
    !! U is proportional to ( 1/r + aa*r  + bb ) such that dU(neigh%cutoff)/dr = 0
    !! therefore aa = 1/(neigh%cutoff)**2 and U(neigh%cutoff) = 0 therefore bb = -2/(neigh%cutoff)
    !!
    !! Note: FS potential can be generalised (R1) by using a damping function
    !! as used for damping the real space coulombic interaction in the
    !! standard Ewald summation.  This generalisation applies when electro%damping > 0.
    !!
    !! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester october 1995
    !! amended   - i.t.todorov november 2014
    !! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                          Intent( In    ) :: iatm
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe,vircpe
    Type( stats_type ),                               Intent( InOut ) :: stats
    Type( electrostatic_type ),                       Intent( InOut ) :: electro
    Type( configuration_type ),                       Intent( InOut ) :: config
    Real( Kind = wp ), Dimension( 9 )                                 :: stress_temp, stress_temp_comp
    Real( Kind = wp ) :: coul

    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,egamma, &
      fix,fiy,fiz,fx,fy,fz

    Logical, save :: newjob = .true.

    If (newjob) Then
      newjob = .false.

      if ( electro%damp ) then
        call electro%erfcgen(neigh%cutoff, electro%damping)

        electro%force_shift =   electro%erfc_deriv%end_sample * neigh%cutoff
        electro%energy_shift = -(electro%erfc%end_sample + electro%force_shift*neigh%cutoff)
      else

        electro%force_shift =  1.0_wp/neigh%cutoff**2
        electro%energy_shift = -2.0_wp/neigh%cutoff ! = -(1.0_wp/neigh%cutoff+aa*neigh%cutoff)

      End If

    End If

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    stress_temp = 0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/electro%eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index and charge

        jatm=neigh%list(m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! Squared distance

          rsq=rrr**2

          ! calculate forces

          If (electro%damp) Then
            egamma = (electro%calc_erfc_deriv(rrr) - electro%force_shift/rrr)*chgprd
          Else
            egamma=chgprd*(1.0_wp/rsq - electro%force_shift)/rrr
          End If

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

            ! calculate potential energy and virial

            If (electro%damp) Then

              ! calculate interaction energy using 3-point interpolation

              coul = (electro%calc_erfc(rrr) + electro%force_shift*rrr + electro%energy_shift)*chgprd
            Else

              coul = chgprd*(1.0_wp/rrr + electro%force_shift*rrr + electro%energy_shift)

            End If

            engcpe = engcpe + coul
            vircpe = vircpe - egamma*rsq

            ! calculate stress tensor

            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stress_temp = stress_temp + stress_temp_comp


#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          If (stats%collect_pp) Then
            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stats%pp_energy(iatm) = stats%pp_energy(iatm) + coul * 0.5_wp
            stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            If (jatm <= config%natms) Then
#endif /* HALF_HALO */
              stats%pp_energy(jatm) = stats%pp_energy(jatm) + coul * 0.5_wp
              stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            End If
#endif /* HALF_HALO */
          End If
          
        End If


      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! complete stress tensor

      stats%stress = stats%stress + stress_temp


    End If

  End Subroutine coul_fscp_forces

  Subroutine coul_rfp_forces(iatm,xxt,yyt,zzt,rrt,engcpe,vircpe,stats,neigh, &
    electro,config)

    !!------------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using reaction field potential (correcting for
    !! the existence of a dipole moment outside rcut)
    !!
    !! Note: RF potential can be generalised (R1) by using a damping function
    !! as used for damping the real space coulombic interaction in the
    !! standard Ewald summation.  This generalisation applies when electro%damping > 0.
    !!
    !! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
    !! R2: M. Neumann, J. Chem. Phys., 82 (12), 5663, (1985)
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester february 1995
    !! amended   - i.t.todorov november 2014
    !! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                          Intent( In    ) :: iatm
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe,vircpe
    Type( stats_type ),                               Intent( InOut ) :: stats
    Type( electrostatic_type ),                       Intent( InOut ) :: electro
    Type( configuration_type ),                       Intent( InOut ) :: config
    Real( Kind = wp ), Dimension( 9 )                                 :: stress_temp, stress_temp_comp
    Real( Kind = wp ) :: coul

    !> Intermediate reaction field variable
    Real( Kind = wp ) :: b0

    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,egamma, &
      fix,fiy,fiz,fx,fy,fz

    Logical, save :: newjob = .true.

    If (newjob) Then
      newjob = .false.

      ! reaction field terms

      b0    = 2.0_wp*(electro%eps - 1.0_wp)/(2.0_wp*electro%eps + 1.0_wp)
      electro%reaction_field(0) = b0/neigh%cutoff**3
      electro%reaction_field(1) = (1.0_wp + 0.5_wp*b0)/neigh%cutoff
      electro%reaction_field(2) = 0.5_wp*electro%reaction_field(0)

      If (electro%damp) Then
        call electro%erfcgen(neigh%cutoff, electro%damping)
        electro%force_shift =   electro%erfc_deriv%end_sample * neigh%cutoff
        electro%energy_shift = -(electro%erfc%end_sample + electro%force_shift*neigh%cutoff)
      End If

    End If

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    stress_temp = 0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/electro%eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index and charge

        jatm=neigh%list(m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! Squared distance

          rsq=rrr**2

          ! calculate forces

          If (electro%damp) Then
            egamma = (electro%calc_erfc_deriv(rrr) - &
              & electro%force_shift/rrr - electro%reaction_field(0))*chgprd
          Else
            egamma=chgprd*(1.0_wp/rsq/rrr - electro%reaction_field(0))
          End If

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

            ! calculate potential energy and virial

            If (electro%damp) Then

              ! calculate interaction energy using 3-point interpolation

              coul = (electro%calc_erfc(rrr) + electro%force_shift*rrr + electro%energy_shift + &
                electro%reaction_field(2)*(rsq-neigh%cutoff_2))*chgprd

            Else

              coul = chgprd*(1.0_wp/rrr + electro%reaction_field(2)*rsq - electro%reaction_field(1))

            End If

            engcpe = engcpe + coul
            vircpe = vircpe - egamma*rsq

            ! calculate stress tensor

            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stress_temp = stress_temp + stress_temp_comp

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          If (stats%collect_pp) Then
            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stats%pp_energy(iatm) = stats%pp_energy(iatm) + coul * 0.5_wp
            stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            If (jatm <= config%natms) Then
#endif /* HALF_HALO */
              stats%pp_energy(jatm) = stats%pp_energy(jatm) + coul * 0.5_wp
              stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            End If
#endif /* HALF_HALO */
          End If

        End If

        
      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! complete stress tensor

      stats%stress = stats%stress + stress_temp
    End If

  End Subroutine coul_rfp_forces

  Subroutine coul_cp_forces(iatm,eps,xxt,yyt,zzt,rrt,engcpe,vircpe,stats,neigh,config)

    !!------------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using 1/r potential with no truncation or damping
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester february 1993
    !! amended   - i.t.todorov november 2014
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                          Intent( In    ) :: iatm
    Real( Kind = wp ),                                Intent( In    ) :: eps
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe,vircpe
    Type( stats_type ),                               Intent( InOut ) :: stats
    Type( configuration_type ),                       Intent( InOut ) :: config
    Real( Kind = wp ), Dimension( 9 )                                 :: stress_temp, stress_temp_comp
    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,coul,fcoul, &
      fix,fiy,fiz,fx,fy,fz

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    stress_temp = 0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index and charge

        jatm=neigh%list(m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! calculate forces

          coul = chgprd/rrr
          fcoul = coul/rrr**2

          fx = fcoul*xxt(m)
          fy = fcoul*yyt(m)
          fz = fcoul*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

            ! calculate potential energy

            engcpe = engcpe + coul

            ! calculate stress tensor

            ! calculate stress tensor

            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stress_temp = stress_temp + stress_temp_comp

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          If (stats%collect_pp) Then
            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stats%pp_energy(iatm) = stats%pp_energy(iatm) + coul * 0.5_wp
            stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            If (jatm <= config%natms) Then
#endif /* HALF_HALO */
              stats%pp_energy(jatm) = stats%pp_energy(jatm) + coul * 0.5_wp
              stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            End If
#endif /* HALF_HALO */
          End If

        End If


      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! virial

      vircpe = -engcpe

      ! complete stress tensor

      stats%stress = stats%stress + stress_temp

    End If
  End Subroutine coul_cp_forces

  Subroutine coul_dddp_forces(iatm,eps,xxt,yyt,zzt,rrt,engcpe,vircpe,stats,neigh,config)

    !!------------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system assuming a distance dependant dielectric
    !! `constant'
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester april 1993
    !! amended   - i.t.todorov november 2014
    !! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                          Intent( In    ) :: iatm
    Real( Kind = wp ),                                Intent( In    ) :: eps
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe,vircpe
    Type( stats_type ),                               Intent( InOut ) :: stats
    Type( configuration_type ),                       Intent( InOut ) :: config
    Real( Kind = wp ), Dimension( 9 )                                 :: stress_temp, stress_temp_comp
    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,rsq,coul,fcoul, &
      fix,fiy,fiz,fx,fy,fz

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    stress_temp = 0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index and charge

        jatm=neigh%list(m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! Squared distance

          rsq=rrr**2

          ! calculate forces

          coul = chgprd/rsq
          fcoul = 2.0_wp*coul/rsq

          fx = fcoul*xxt(m)
          fy = fcoul*yyt(m)
          fz = fcoul*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

            ! calculate potential energy

            engcpe = engcpe + coul

            ! calculate stress tensor

            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stress_temp = stress_temp + stress_temp_comp

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          If (stats%collect_pp) Then
            stress_temp_comp = calculate_stress( [xxt(m), yyt(m), zzt(m)], [fx,fy,fz] )
            stats%pp_energy(iatm) = stats%pp_energy(iatm) + coul * 0.5_wp
            stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            If (jatm <= config%natms) Then
#endif /* HALF_HALO */
              stats%pp_energy(jatm) = stats%pp_energy(jatm) + coul * 0.5_wp
              stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
            End If
#endif /* HALF_HALO */
          End If

        End If


      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! virial

      vircpe = -2.0_wp*engcpe

      ! complete stress tensor

      stats%stress = stats%stress + stress_temp

    End If
  End Subroutine coul_dddp_forces
End Module coul_spole
