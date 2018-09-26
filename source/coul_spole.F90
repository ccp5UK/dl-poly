Module coul_spole
  Use kinds,           Only : wp
  Use comms,           Only : comms_type
  Use configuration,   Only : configuration_type     
  Use particle,         Only : corePart
  Use setup,           Only :  r4pie0, zero_plus,  nrite,sqrpi
  Use errors_warnings, Only : error
  Use numerics,        Only : erfcgen
  Use neighbours,      Only : neighbours_type
  Use electrostatic,   Only : electrostatic_type, &
                            ELECTROSTATIC_EWALD,ELECTROSTATIC_DDDP, &
                            ELECTROSTATIC_COULOMB,ELECTROSTATIC_COULOMB_FORCE_SHIFT, &
                            ELECTROSTATIC_COULOMB_REACTION_FIELD,ELECTROSTATIC_POISSON

  Implicit None

  Private

  Public :: intra_coul,coul_fscp_forces, coul_rfp_forces, coul_cp_forces, coul_dddp_forces

  Contains

  Subroutine intra_coul(rcut,chgprd,rrr,rsq,coul,fcoul,safe,electro)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating bond's or 1-4 dihedral
  ! electrostatics: adjusted by a config%weighting factor
  !
  ! copyright - daresbury laboratory
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent( In    ) :: chgprd,rcut,rrr,rsq
    Real( Kind = wp ), Intent(   Out ) :: coul,fcoul
    Logical,           Intent( InOut ) :: safe
    Type( electrostatic_type ), Intent( InOut ) :: electro


    Real( Kind = wp ) :: exp1,tt,erc,fer,b0

    Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

    If (electro%newjob_intra) Then
      electro%newjob_intra = .false.

      ! Check for damped force-shifted coulombic and reaction field interactions
      ! and set force and potential shifting parameters dependingly

      electro%damp=.false.
      If (electro%alpha > zero_plus) Then
        electro%damp=.true.

        exp1= Exp(-(electro%alpha*rcut)**2)
        tt  = 1.0_wp/(1.0_wp+pp*electro%alpha*rcut)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rcut
        fer = (erc + 2.0_wp*(electro%alpha/sqrpi)*exp1)/rcut**2

        electro%aa  = fer*rcut
        electro%bb  = -(erc + electro%aa*rcut)
      Else If (electro%key == ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then
        electro%aa =  1.0_wp/rcut**2
        electro%bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)
      End If

      ! set reaction field terms for RFC

      If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then
        b0    = 2.0_wp*(electro%eps - 1.0_wp)/(2.0_wp*electro%eps + 1.0_wp)
        electro%rfld0 = b0/rcut**3
        electro%rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
        electro%rfld2 = 0.5_wp*electro%rfld0
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
          exp1= Exp(-(electro%alpha*rrr)**2)
          tt  = 1.0_wp/(1.0_wp+pp*electro%alpha*rrr)

          erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rrr
          fer = (erc + 2.0_wp*(electro%alpha/sqrpi)*exp1)/rsq

          coul = chgprd*(erc + electro%aa*rrr + electro%bb)
          fcoul= chgprd*(fer - electro%aa/rrr)
        End If

        If (electro%key ==  ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then ! force shifted coulombic
          If (.not.electro%damp) Then ! pure
            coul = chgprd*(1.0_wp/rrr + electro%aa*rrr+ electro%bb)
            fcoul= chgprd*(1.0_wp/rsq - electro%aa)/rrr
          Else                ! damped
             coul = coul
             fcoul= fcoul
          End If
       Else If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then ! reaction field
          If (.not.electro%damp) Then ! pure
            coul = chgprd*(1.0_wp/rrr + electro%rfld2*rsq - electro%rfld1)
            fcoul= chgprd*(1.0_wp/rsq/rrr - electro%rfld0)
          Else                ! damped
            coul = coul  + chgprd*(electro%rfld2*rsq - electro%rfld1)
            fcoul= fcoul + chgprd*(-electro%rfld0)
          End If
       End If

    Else

       safe = .false.

    End If
  End Subroutine intra_coul

  Subroutine coul_fscp_forces(iatm,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh, &
      electro,config,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system assuming a force shifted coulombic potential
  !
  ! U is proportional to ( 1/r + aa*r  + bb ) such that dU(neigh%cutoff)/dr = 0
  ! therefore aa = 1/(neigh%cutoff)**2 and U(neigh%cutoff) = 0 therefore bb = -2/(neigh%cutoff)
  !
  ! Note: FS potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when electro%alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester october 1995
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( comms_type ),                       Intent( In    ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config


    Integer           :: fail,k,idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,ppp,egamma, &
      fix,fiy,fiz,fx,fy,fz,            &
      vk0,vk1,vk2,gk0,gk1,gk2,t1,t2,   &
      strs1,strs2,strs3,strs5,strs6,strs9

    Character ( Len = 256 ) :: message

    If (electro%newjob_fscp) Then
      electro%newjob_fscp = .false.

      If (electro%alpha > zero_plus) Then
        electro%damp_fscp = .true.
      Else
        electro%damp_fscp = .false.
      End If

      If (electro%damp_fscp) Then

        ! interpolation interval

        electro%drewd_fscp = neigh%cutoff/Real(electro%ewald_exclusion_grid-4,wp)

        ! reciprocal of interpolation interval

        electro%rdrewd_fscp = 1.0_wp/electro%drewd_fscp

        fail=0
        Allocate (electro%erc_fscp(0:electro%ewald_exclusion_grid),electro%fer_fscp(0:electro%ewald_exclusion_grid), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'coul_fscp_forces allocation failure'
          Call error(0,message)
        End If

        ! generate error function complement tables for ewald sum

        Call erfcgen(neigh%cutoff,electro%alpha,electro%ewald_exclusion_grid,electro%erc_fscp,electro%fer_fscp)

        ! set force and potential shifting parameters (screened terms)

        electro%aa_fscp =   electro%fer_fscp(electro%ewald_exclusion_grid-4)*neigh%cutoff
        electro%bb_fscp = -(electro%erc_fscp(electro%ewald_exclusion_grid-4)+electro%aa_fscp*neigh%cutoff)

      Else

        ! set force and potential shifting parameters (screened terms)

        electro%aa_fscp =  1.0_wp/neigh%cutoff**2
        electro%bb_fscp = -2.0_wp/neigh%cutoff ! = -(1.0_wp/neigh%cutoff+aa*neigh%cutoff)

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

             If (electro%damp_fscp) Then
               k   = Int(rrr*electro%rdrewd_fscp)
               ppp = rrr*electro%rdrewd_fscp - Real(k,wp)

               ! calculate forces using 3pt interpolation

               gk0 = electro%fer_fscp(k) ; If (k == 0) gk0 = gk0*rrr
               gk1 = electro%fer_fscp(k+1)
               gk2 = electro%fer_fscp(k+2)

               t1 = gk0 + (gk1 - gk0)*ppp
               t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

               egamma = ((t1 + (t2-t1)*ppp*0.5_wp) - electro%aa_fscp/rrr)*chgprd
             Else
               egamma=chgprd*(1.0_wp/rsq - electro%aa_fscp)/rrr
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

               If (electro%damp_fscp) Then

  ! calculate interaction energy using 3-point interpolation

                 vk0 = electro%erc_fscp(k)
                 vk1 = electro%erc_fscp(k+1)
                 vk2 = electro%erc_fscp(k+2)

                   t1 = vk0 + (vk1 - vk0)*ppp
                   t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                   engcpe = engcpe + ((t1 + (t2-t1)*ppp*0.5_wp) + electro%aa_fscp*rrr + electro%bb_fscp)*chgprd

                 Else

                   engcpe = engcpe + chgprd*(1.0_wp/rrr + electro%aa_fscp*rrr + electro%bb_fscp)

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
      electro,config,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using reaction field potential (correcting for
  ! the existence of a dipole moment outside rcut)
  !
  ! Note: RF potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when electro%alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  ! R2: M. Neumann, J. Chem. Phys., 82 (12), 5663, (1985)
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester february 1995
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( comms_type ),                       Intent( In    ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config



    Integer           :: fail,k,idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,ppp,egamma, &
                         fix,fiy,fiz,fx,fy,fz,            &
                         vk0,vk1,vk2,gk0,gk1,gk2,t1,t2,   &
                         strs1,strs2,strs3,strs5,strs6,strs9

    Character ( Len = 256 )  :: message


    If (electro%newjob_rfp) Then
      electro%newjob_rfp = .false.

      If (electro%alpha > zero_plus) Then
        electro%damp_rfp = .true.
      Else
        electro%damp_rfp = .false.
      End If

      ! reaction field terms

      electro%b0_rfp    = 2.0_wp*(electro%eps - 1.0_wp)/(2.0_wp*electro%eps + 1.0_wp)
      electro%rfld0_rfp = electro%b0_rfp/neigh%cutoff**3
      electro%rfld1_rfp = (1.0_wp + 0.5_wp*electro%b0_rfp)/neigh%cutoff
      electro%rfld2_rfp = 0.5_wp*electro%rfld0_rfp

      If (electro%damp_rfp) Then

        ! interpolation interval

        electro%drewd_rfp = neigh%cutoff/Real(electro%ewald_exclusion_grid-4,wp)

        ! reciprocal of interpolation interval

        electro%rdrewd_rfp = 1.0_wp/electro%drewd_rfp

        fail=0
        Allocate (electro%erc_rfp(0:electro%ewald_exclusion_grid),electro%fer_rfp(0:electro%ewald_exclusion_grid), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'coul_fscp_forces allocation failure'
          Call error(0,message)
        End If

        ! generate error function complement tables for ewald sum

        Call erfcgen(neigh%cutoff,electro%alpha,electro%ewald_exclusion_grid,electro%erc_rfp,electro%fer_rfp)

        ! set force and potential shifting parameters (screened terms)

        electro%aa_rfp =   electro%fer_rfp(electro%ewald_exclusion_grid-4)*neigh%cutoff
        electro%bb_rfp = -(electro%erc_rfp(electro%ewald_exclusion_grid-4)+electro%aa_rfp*neigh%cutoff)

        ! Cutoff squared

        electro%rcsq_rfp = neigh%cutoff**2

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

             If (electro%damp_rfp) Then
               k   = Int(rrr*electro%rdrewd_rfp)
               ppp = rrr*electro%rdrewd_rfp - Real(k,wp)

                ! calculate forces using 3pt interpolation

                gk0 = electro%fer_rfp(k) ; If (k == 0) gk0 = gk0*rrr
                gk1 = electro%fer_rfp(k+1)
                gk2 = electro%fer_rfp(k+2)

                t1 = gk0 + (gk1 - gk0)*ppp
                t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                egamma = ((t1 + (t2-t1)*ppp*0.5_wp) - electro%aa_rfp/rrr - electro%rfld0_rfp)*chgprd
              Else
                egamma=chgprd*(1.0_wp/rsq/rrr - electro%rfld0_rfp)
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

               If (electro%damp_rfp) Then

  ! calculate interaction energy using 3-point interpolation

                  vk0 = electro%erc_rfp(k)
                  vk1 = electro%erc_rfp(k+1)
                  vk2 = electro%erc_rfp(k+2)

                   t1 = vk0 + (vk1 - vk0)*ppp
                   t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                   engcpe = engcpe + ((t1 + (t2-t1)*ppp*0.5_wp) + electro%aa_rfp*rrr + electro%bb_rfp + &
                     electro%rfld2_rfp*(rsq-electro%rcsq_rfp))*chgprd

                 Else

                   engcpe = engcpe + chgprd*(1.0_wp/rrr + electro%rfld2_rfp*rsq - electro%rfld1_rfp)

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using 1/r potential with no truncation or damping
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester february 1993
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system assuming a distance dependant dielectric
  ! `constant'
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester april 1993
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
