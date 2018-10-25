Module coul_spole
  !! This module has no header !
  Use kinds,           Only : wp
  Use configuration,   Only : configuration_type     
  Use constants,       Only : r4pie0,zero_plus,sqrpi
  Use errors_warnings, Only : error, error_alloc, error_dealloc
  Use numerics,        Only : erfcgen, erfc, erfc_deriv, three_p_interp
  Use numerics,        Only : calc_erfc, calc_erfc_deriv
  Use neighbours,      Only : neighbours_type
  Use electrostatic,   Only : electrostatic_type, &
    ELECTROSTATIC_EWALD,ELECTROSTATIC_DDDP, &
    ELECTROSTATIC_COULOMB,ELECTROSTATIC_COULOMB_FORCE_SHIFT, &
    ELECTROSTATIC_COULOMB_REACTION_FIELD

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

        exp1= Exp(-(electro%damping*rcut)**2)
        tt  = 1.0_wp/(1.0_wp+pp*electro%damping*rcut)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rcut
        fer = (erc + 2.0_wp*(electro%damping/sqrpi)*exp1)/rcut**2

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

  Subroutine coul_fscp_forces(iatm,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh, &
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
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( configuration_type ),               Intent( InOut ) :: config


    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,egamma, &
      fix,fiy,fiz,fx,fy,fz,            &
      strs1,strs2,strs3,strs5,strs6,strs9

    Logical, save :: newjob = .true.

    If (newjob) Then
      newjob = .false.

      if ( electro%damp ) then
        call erfcgen(neigh%cutoff, electro%damping, erfc, erfc_deriv)

        electro%force_shift =   erfc_deriv%table(erfc_deriv%nsamples-4)*neigh%cutoff
        electro%energy_shift = -(erfc%table(erfc%nsamples-4)+electro%force_shift*neigh%cutoff)
      else

        electro%force_shift =  1.0_wp/neigh%cutoff**2
        electro%energy_shift = -2.0_wp/neigh%cutoff ! = -(1.0_wp/neigh%cutoff+aa*neigh%cutoff)

      End If

    End If

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

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
            egamma = (three_p_interp(erfc_deriv,rrr) - electro%force_shift/rrr)*chgprd
          Else
            egamma=chgprd*(1.0_wp/rsq - electro%force_shift)/rrr
          End If

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! calculate potential energy and virial

            If (electro%damp) Then

              ! calculate interaction energy using 3-point interpolation

              engcpe = engcpe + (three_p_interp(erfc,rrr) + electro%force_shift*rrr + electro%energy_shift)*chgprd
            Else

              engcpe = engcpe + chgprd*(1.0_wp/rrr + electro%force_shift*rrr + electro%energy_shift)

            End If

            vircpe = vircpe - egamma*rsq

            ! calculate stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! complete stress tensor

      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9

    End If
  End Subroutine coul_fscp_forces

  Subroutine coul_rfp_forces(iatm,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh, &
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
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( configuration_type ),               Intent( InOut ) :: config

    !> Intermediate reaction field variable
    Real( Kind = wp ) :: b0

    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,egamma, &
      fix,fiy,fiz,fx,fy,fz, strs1,strs2,strs3,strs5,strs6,strs9

    Logical, save :: newjob = .true.

    If (newjob) Then
      newjob = .false.

      ! reaction field terms

      b0    = 2.0_wp*(electro%eps - 1.0_wp)/(2.0_wp*electro%eps + 1.0_wp)
      electro%reaction_field(0) = b0/neigh%cutoff**3
      electro%reaction_field(1) = (1.0_wp + 0.5_wp*b0)/neigh%cutoff
      electro%reaction_field(2) = 0.5_wp*electro%reaction_field(0)

      If (electro%damp) Then
        call erfcgen(neigh%cutoff, electro%damping, erfc, erfc_deriv)
        electro%force_shift =   erfc_deriv%table(erfc_deriv%nsamples-4)*neigh%cutoff
        electro%energy_shift = -(erfc%table      (erfc%nsamples      -4)+electro%force_shift*neigh%cutoff)
      End If

    End If

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

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
            egamma = (three_p_interp(erfc_deriv,rrr) - electro%force_shift/rrr - electro%reaction_field(0))*chgprd
          Else
            egamma=chgprd*(1.0_wp/rsq/rrr - electro%reaction_field(0))
          End If

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! calculate potential energy and virial

            If (electro%damp) Then

              ! calculate interaction energy using 3-point interpolation

              engcpe = engcpe + (three_p_interp(erfc,rrr) + electro%force_shift*rrr + electro%energy_shift + &
                electro%reaction_field(2)*(rsq-neigh%cutoff_2))*chgprd

            Else

              engcpe = engcpe + chgprd*(1.0_wp/rrr + electro%reaction_field(2)*rsq - electro%reaction_field(1))

            End If

            vircpe = vircpe - egamma*rsq

            ! calculate stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! complete stress tensor

      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9

    End If
  End Subroutine coul_rfp_forces

  Subroutine coul_cp_forces(iatm,eps,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh,config)

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

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: eps
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( configuration_type ),               Intent( InOut ) :: config

    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,coul,fcoul, &
      fix,fiy,fiz,fx,fy,fz,        &
      strs1,strs2,strs3,strs5,strs6,strs9

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

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

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! calculate potential energy

            engcpe = engcpe + coul

            ! calculate stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

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

      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9

    End If
  End Subroutine coul_cp_forces

  Subroutine coul_dddp_forces(iatm,eps,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh,config)

    !!------------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system assuming a distance dependant dielectric
    !! `constant'
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester april 1993
    !! amended   - i.t.todorov november 2014
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!------------------------------------------------------------------------!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: eps
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( configuration_type ),               Intent( InOut ) :: config

    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,rsq,coul,fcoul, &
      fix,fiy,fiz,fx,fy,fz,            &
      strs1,strs2,strs3,strs5,strs6,strs9

    ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

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

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! calculate potential energy

            engcpe = engcpe + coul

            ! calculate stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

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

      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9

    End If
  End Subroutine coul_dddp_forces
End Module coul_spole
