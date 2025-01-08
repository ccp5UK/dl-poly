Module dpd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global DPD variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2016
  ! contrib   - m.a.seaton august 2020
  !           - k.a.jonathan september 2024
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: DpdVExp_tag,&
                             comms_type,&
                             gcheck,&
                             girecv,&
                             gmax,&
                             gsend,&
                             gsum,&
                             gwait
  Use configuration,   Only: configuration_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             warning
  Use kinds,           Only: wp,STR_LEN
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: box_mueller_saru2,&
                             images,&
                             seed_type
  Use rigid_bodies,    Only: getrotmat, rigid_bodies_type
  Use shared_units,    Only: SHARED_UNIT_UPDATE_FORCES,&
                             update_shared_units
  Use statistics,      Only: stats_type
  Use thermostat,      Only: DPD_ZEROTH_ORDER,&
                             DPD_FIRST_ORDER,&
                             DPD_SECOND_ORDER,&
                             DPD_MDVV,&
                             DPD_NULL,&
                             VV_FIRST_STAGE,&
                             VV_SECOND_STAGE,&
                             thermostat_type
#ifdef HALF_HALO
  Use numerics,        Only: local_index
#endif /* HALF_HALO */

  Implicit None

  Private

  Public :: dpd_shardlow_integrate, dpd_mdvv_forces

Contains

  Subroutine dpd_shardlow_integrate(stage, l_str, rcut, nstep, tstep, stats, thermo, neigh, rigid, domain, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine applying DPD thermostat in a Shardlow splitting 
    ! manner using the verlet neighbour neigh%list
    !
    ! thermo%key_dpd = DPD_ZEROTH_ORDER for zeroth order splitting
    ! thermo%key_dpd = DPD_FIRST_ORDER for first order splitting
    ! thermo%key_dpd = DPD_SECOND_ORDER for second order splitting
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    ! contrib   - i.t.todorov may 2020 - 'half-halo' VNL
    !           - m.a.seaton august 2020 - preprocessing tags and array sizes
    !           - k.a.jonathan september 2024 - zeroth order splitting,
    !             integration of RBs vels, reduced mass in equations
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: l_str
    Real(Kind=wp),            Intent(In   ) :: rcut
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(In   ) :: tstep
    Type(stats_type),         Intent(InOut) :: stats
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(neighbours_type),    Intent(In   ) :: neigh
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(len=256)                       :: message
    Integer                                  :: ai, aj, fail, i, idi, idj, j, k, key, limit, &
                                                nst_p
    Real(Kind=wp)                            :: dgamma, gamma, gauss, &
                                                hstep, rgamma, rrr, r_sqrt_tstp, scl, scrn, &
                                                tmp, tst_p, rmassij, xdif, ydif, zdif, &
                                                vxdif, vydif, vzdif, rdotv, tstepfrac, &
                                                strsdrag(9), strsrand(9)
    Real(Kind=wp), Allocatable, Dimension(:) :: fdpdx, fdpdy, fdpdz

    If ( (thermo%key_dpd == DPD_NULL) .or. &
         (thermo%key_dpd == DPD_MDVV) .or. &
         (thermo%key_dpd /= DPD_SECOND_ORDER .and. stage == VV_SECOND_STAGE) ) Return

    fail = 0
    Allocate (fdpdx(1:config%mxatms), fdpdy(1:config%mxatms), fdpdz(1:config%mxatms), Stat=fail)

    If (fail > 0) Then
      Write (message, '(a)') 'dpd_shardlow_integrate force array allocation failure'
      Call error(0, message)
    End If

    ! set effective timestep, tst_p, and random number seed, nst_p, wrt splitting order

    If (thermo%key_dpd == DPD_ZEROTH_ORDER) Then

      nst_p = nstep
      ! double tstep accounts for integration over full timestep 
      ! in one pass for zeroth order splitting
      tst_p = 2.0_wp * tstep

    Else If (thermo%key_dpd == DPD_FIRST_ORDER) Then
      
      nst_p = nstep
      tst_p = tstep

    Else ! DPD_SECOND_ORDER

      ! accounts for integration over half timestep twice in second order splitting
      tst_p = 0.5_wp * tstep

      If (stage == VV_FIRST_STAGE) Then

        nst_p = nstep

      Else ! VV_SECOND_STAGE

        ! ensures different random number in second shardlow call for second order splitting
        nst_p = -nstep

      End If

    End If

    ! Set tstep derivatives

    hstep = 0.5_wp * tst_p
    tstepfrac = hstep / tstep
    r_sqrt_tstp = 1.0_wp / Sqrt(tstep)

    ! random force scaled in 2nd order to account for symmetric application over
    ! effectively a timestep of half the size

    If (thermo%key_dpd == DPD_SECOND_ORDER) Then
      r_sqrt_tstp = r_sqrt_tstp * Sqrt(2.0_wp)
    End If

    ! Initialise DPD virial and stress contributions

    If (stage == VV_FIRST_STAGE) Then
      stats%virdpd = 0.0_wp
      stats%strdpdr = 0.0_wp
      stats%strdpdd = 0.0_wp
    End If

    ! Initialise local force and stress arrays

    fdpdx = 0.0_wp
    fdpdy = 0.0_wp
    fdpdz = 0.0_wp

    strsdrag = 0.0_wp
    strsrand = 0.0_wp

    ! Refresh halo velocities

    Call dpd_v_set_halo(domain, config, comm)

    ! First pass of shardlow ignored in zeroth order splitting
    If (thermo%key_dpd /= DPD_ZEROTH_ORDER) Then

      ! FIRST PASS

      ! outer loop over atoms

      Do i = 1, config%natms

        ! primary atom type and global index

        ai = config%ltype(i)
        idi = config%ltg(i)

        ! loop over all valid pairs - ignore massless particles

        limit = Merge(neigh%list(0, i), 0, config%weight(i) > 1.0e-6_wp)

        Do k = 1, limit

          ! Secondary atom index, type, and global index

          j = neigh%list(k,i)
          aj = config%ltype(j)
          idj = config%ltg(j)

          ! Calculate r_ij, |r_ij|

          xdif = config%parts(i)%xxx - config%parts(j)%xxx
          ydif = config%parts(i)%yyy - config%parts(j)%yyy
          zdif = config%parts(i)%zzz - config%parts(j)%zzz

          rrr = Sqrt(xdif**2 + ydif**2 + zdif**2)

          ! Calculate forces for valid pairs

          If (rrr < rcut .and. config%weight(j) > 1.0e-6_wp) Then

            ! Calculate v_ij, and (r_ij . v_ij)

            vxdif = config%vxx(i) - config%vxx(j)
            vydif = config%vyy(i) - config%vyy(j)
            vzdif = config%vzz(i) - config%vzz(j)

            rdotv = xdif*vxdif + ydif*vydif + zdif*vzdif

            ! Get mixing type function - key for interaction strength

            If (ai > aj) Then
              key = ai * (ai - 1) / 2 + aj
            Else
              key = aj * (aj - 1) / 2 + ai
            End If

            ! Get gaussian random number with zero mean (held in gauss var)
            ! Global id check ensure same random number for same pair of particles

            If (idi < idj) Then
              Call box_mueller_saru2(seed, idi, idj, nst_p, gauss, l_str)
            Else
              Call box_mueller_saru2(seed, idj, idi, nst_p, gauss, l_str)
            End If

            ! Screening function, related to the drag and random weight functions:
            ! w_D = scrn**2 * rrr**2
            ! w_R = scrn * rrr

            scrn = (rcut - rrr) / (rrr * rcut)         

            ! Calculate random and drag components

            rgamma = thermo%sigdpd(key) * scrn * gauss * r_sqrt_tstp
            dgamma = - thermo%gamdpd(key) * scrn**2 * rdotv

            ! Total force component gamma_ij such that, summed over all pairs
            ! v_i(t+hstep) = v_i(t) + (hstep / m_i) * (gamma_ij * r_ij)

            gamma = rgamma + dgamma

            ! Update forces

            fdpdx(i) = fdpdx(i) + gamma * xdif
            fdpdy(i) = fdpdy(i) + gamma * ydif
            fdpdz(i) = fdpdz(i) + gamma * zdif
      
#ifndef HALF_HALO
            If (j <= config%natms) Then
#endif /* HALF_HALO */

              fdpdx(j) = fdpdx(j) - gamma * xdif
              fdpdy(j) = fdpdy(j) - gamma * ydif
              fdpdz(j) = fdpdz(j) - gamma * zdif
              
#ifndef HALF_HALO
            End If
#endif /* HALF_HALO */

            !     Assign stress terms (only when second particle is in subdomain
            !     or has larger global particle index than first to avoid double-counting)

#ifndef HALF_HALO
            If (j <= config%natms .or. idi < idj) Then
#endif /* HALF_HALO */

              strsrand(1) = strsrand(1) + rgamma * xdif * xdif * tstepfrac ! random stress_xx
              strsrand(2) = strsrand(2) + rgamma * ydif * xdif * tstepfrac ! random stress_xy
              strsrand(3) = strsrand(3) + rgamma * zdif * xdif * tstepfrac ! random stress_xz
              strsrand(5) = strsrand(5) + rgamma * ydif * ydif * tstepfrac ! random stress_yy
              strsrand(6) = strsrand(6) + rgamma * ydif * zdif * tstepfrac ! random stress_yz
              strsrand(9) = strsrand(9) + rgamma * zdif * zdif * tstepfrac ! random stress_zz
              
              strsdrag(1) = strsdrag(1) + dgamma * xdif * xdif * tstepfrac ! drag stress_xx
              strsdrag(2) = strsdrag(2) + dgamma * ydif * xdif * tstepfrac ! drag stress_xy
              strsdrag(3) = strsdrag(3) + dgamma * zdif * xdif * tstepfrac ! drag stress_xz
              strsdrag(5) = strsdrag(5) + dgamma * ydif * ydif * tstepfrac ! drag stress_yy
              strsdrag(6) = strsdrag(6) + dgamma * ydif * zdif * tstepfrac ! drag stress_yz
              strsdrag(9) = strsdrag(9) + dgamma * zdif * zdif * tstepfrac ! drag stress_zz

#ifndef HALF_HALO
            End If
#endif /* HALF_HALO */

          End If

        End Do

      End Do

#ifdef HALF_HALO
      ! Share the dpd forces collected in the halo with the parent domains

      Call refresh_halo_dpd_forces(domain, config, config%mxatms, fdpdx, fdpdy, fdpdz, comm)

#endif /* HALF_HALO */

      ! Update velocities

      Do i = 1, config%natms
        If (config%lfree(i) == 0) Then
          If (config%weight(i) > 1.0e-6_wp) Then
            tmp = hstep / config%weight(i)
            config%vxx(i) = config%vxx(i) + tmp * fdpdx(i)
            config%vyy(i) = config%vyy(i) + tmp * fdpdy(i)
            config%vzz(i) = config%vzz(i) + tmp * fdpdz(i)
          End If
        End If
      End Do

      ! Share and update dpd forces for any RBs shared across domains
      
      If (rigid%share) Then
        Call update_shared_units(config, rigid%list_shared, &
                                rigid%map_shared, fdpdx, fdpdy, fdpdz, domain, comm)
      End If    

      ! Integrate velocity updates of RBs due to dpd forces

      Call integrate_rbs_dpd(config, rigid, fdpdx, fdpdy, fdpdz, tst_p)

    End If

    ! SECOND PASS

    ! Re-initialise local force arrays

    fdpdx = 0.0_wp
    fdpdy = 0.0_wp
    fdpdz = 0.0_wp

    ! Refresh halo velocities

    Call dpd_v_set_halo(domain, config, comm)

    Do i = 1, config%natms

      ! primary atom type and global index

      ai = config%ltype(i)
      idi = config%ltg(i)

      ! loop over all valid pairs - ignore massless particles

      limit = Merge(neigh%list(0, i), 0, config%weight(i) > 1.0e-6_wp)

      Do k = 1, limit

        ! Secondary atom index, type, and global index

        j = neigh%list(k,i)
        aj = config%ltype(j)
        idj = config%ltg(j)

        ! Calculate r_ij, |r_ij|

        xdif = config%parts(i)%xxx - config%parts(j)%xxx
        ydif = config%parts(i)%yyy - config%parts(j)%yyy
        zdif = config%parts(i)%zzz - config%parts(j)%zzz

        rrr = Sqrt(xdif**2 + ydif**2 + zdif**2)

        ! Calculate forces for valid pairs

        If (rrr < rcut .and. config%weight(j) > 1.0e-6_wp) Then

          ! Calculate v_ij, and (r_ij . v_ij)

          vxdif = config%vxx(i) - config%vxx(j)
          vydif = config%vyy(i) - config%vyy(j)
          vzdif = config%vzz(i) - config%vzz(j)

          rdotv = xdif*vxdif + ydif*vydif + zdif*vzdif

          ! Get mixing type function - key for interaction strength

          If (ai > aj) Then
            key = ai * (ai - 1) / 2 + aj
          Else
            key = aj * (aj - 1) / 2 + ai
          End If

          ! Get gaussian random number with zero mean (held in gauss var)
          ! Global id check ensure same random number for same pair of particles

          If (idi < idj) Then
            Call box_mueller_saru2(seed, idi, idj, nst_p, gauss, l_str)
          Else
            Call box_mueller_saru2(seed, idj, idi, nst_p, gauss, l_str)
          End If

          ! Screening function, related to the drag and random weight functions:
          ! w_D = scrn**2 * rrr**2
          ! w_R = scrn * rrr

          scrn = (rcut - rrr) / (rrr * rcut)  
          
          ! reciprocal reduced mass
          
          rmassij = (config%weight(i) + config%weight(j)) / (config%weight(i) * config%weight(j))

          ! Calculate random and drag components

          tmp = thermo%sigdpd(key) * scrn * gauss * r_sqrt_tstp

          scl = (thermo%gamdpd(key) * scrn**2) / (1.0_wp + hstep * thermo%gamdpd(key) * scrn**2 * rrr**2 * rmassij)

          rgamma = (1.0_wp - scl * 0.25_wp * rmassij * rrr**2 * tst_p) * tmp

          dgamma = - scl * rdotv

          ! Total force component gamma_ij such that, summed over all pairs
          ! v_i(t+hstep) = v_i(t) + (hstep / m_i) * (gamma_ij * r_ij)

          gamma = rgamma + dgamma

          ! Update forces

          fdpdx(i) = fdpdx(i) + gamma * xdif
          fdpdy(i) = fdpdy(i) + gamma * ydif
          fdpdz(i) = fdpdz(i) + gamma * zdif
    
#ifndef HALF_HALO
          If (j <= config%natms) Then
#endif /* HALF_HALO */

            fdpdx(j) = fdpdx(j) - gamma * xdif
            fdpdy(j) = fdpdy(j) - gamma * ydif
            fdpdz(j) = fdpdz(j) - gamma * zdif
          
#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          !     Assign stress terms (only when second particle is in subdomain
          !     or has larger global particle index than first to avoid double-counting)

#ifndef HALF_HALO
          If (j <= config%natms .or. idi < idj) Then
#endif /* HALF_HALO */

            strsrand(1) = strsrand(1) + rgamma * xdif * xdif * tstepfrac ! random stress_xx
            strsrand(2) = strsrand(2) + rgamma * ydif * xdif * tstepfrac ! random stress_xy
            strsrand(3) = strsrand(3) + rgamma * zdif * xdif * tstepfrac ! random stress_xz
            strsrand(5) = strsrand(5) + rgamma * ydif * ydif * tstepfrac ! random stress_yy
            strsrand(6) = strsrand(6) + rgamma * ydif * zdif * tstepfrac ! random stress_yz
            strsrand(9) = strsrand(9) + rgamma * zdif * zdif * tstepfrac ! random stress_zz
            
            strsdrag(1) = strsdrag(1) + dgamma * xdif * xdif * tstepfrac ! drag stress_xx
            strsdrag(2) = strsdrag(2) + dgamma * ydif * xdif * tstepfrac ! drag stress_xy
            strsdrag(3) = strsdrag(3) + dgamma * zdif * xdif * tstepfrac ! drag stress_xz
            strsdrag(5) = strsdrag(5) + dgamma * ydif * ydif * tstepfrac ! drag stress_yy
            strsdrag(6) = strsdrag(6) + dgamma * ydif * zdif * tstepfrac ! drag stress_yz
            strsdrag(9) = strsdrag(9) + dgamma * zdif * zdif * tstepfrac ! drag stress_zz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

        End If

      End Do

    End Do

#ifdef HALF_HALO
    ! Share the dpd forces collected in the halo with the parent domains

    Call refresh_halo_dpd_forces(domain, config, config%mxatms, fdpdx, fdpdy, fdpdz, comm)

#endif /* HALF_HALO */

    ! Update velocities

    Do i = 1, config%natms
      If (config%lfree(i) == 0) Then
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * fdpdx(i)
          config%vyy(i) = config%vyy(i) + tmp * fdpdy(i)
          config%vzz(i) = config%vzz(i) + tmp * fdpdz(i)
        End If
      End If
    End Do

    ! Share and update dpd forces for any RBs shared across domains
      
    If (rigid%share) Then
      Call update_shared_units(config, rigid%list_shared, &
                              rigid%map_shared, fdpdx, fdpdy, fdpdz, domain, comm)
    End If    

    ! Integrate velocity updates of RBs due to dpd forces

    Call integrate_rbs_dpd(config, rigid, fdpdx, fdpdy, fdpdz, tst_p)

    ! Symmetrise and globalise random and drag stresses

    strsrand(4) = strsrand(2) ! random stress_yx
    strsrand(7) = strsrand(3) ! random stress_zx
    strsrand(8) = strsrand(6) ! random stress_zy

    Call gsum(comm, strsrand)

    strsdrag(4) = strsdrag(2) ! drag stress_yx
    strsdrag(7) = strsdrag(3) ! drag stress_zx
    strsdrag(8) = strsdrag(6) ! drag stress_zy

    Call gsum(comm, strsdrag)

    ! Load temporary dpd stress tensors into stats

    stats%strdpdd = stats%strdpdd + strsdrag
    stats%strdpdr = stats%strdpdr + strsrand

    ! Update virial (vir = - Tr(strdpdr) - Tr(strdpdd))

    stats%virdpd = stats%virdpd &
                   - stats%strdpdd(1) - stats%strdpdd(5) - stats%strdpdd(9) &
                   - stats%strdpdr(1) - stats%strdpdr(5) - stats%strdpdr(9)

    Deallocate (fdpdx, fdpdy, fdpdz, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'dpd_shardlow_integrate force array deallocation failure'
      Call error(0, message)
    End If

  End Subroutine dpd_shardlow_integrate

  Subroutine dpd_mdvv_forces(stage, l_str, rcut, nstep, tstep, stats, thermo, neigh, rigid, domain, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine applying DPD thermostat in a traditional
    ! molecular dynamics velocity verlet manner - forces are calculated
    ! using the verlet neighbour neigh%list and added to the conservative
    ! forces (config%parts%fxx,yy,zz) for inclusion in standard VV integration
    !
    ! copyright - daresbury laboratory
    ! author    - k.a.jonathan september 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: l_str
    Real(Kind=wp),            Intent(In   ) :: rcut
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(In   ) :: tstep
    Type(stats_type),         Intent(InOut) :: stats
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(neighbours_type),    Intent(In   ) :: neigh
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(len=256)                       :: message
    Integer                                  :: ai, aj, fail, i, idi, idj, j, k, key, limit
    Real(Kind=wp)                            :: dgamma, gamma, gauss, &
                                                hstep, rgamma, rrr, r_sqrt_tstp, scrn, &
                                                xdif, ydif, zdif, vxdif, vydif, vzdif, rdotv
    Real(Kind=wp), Allocatable, Dimension(:) :: fdpdx, fdpdy, fdpdz

    If (thermo%key_dpd /= DPD_MDVV .or. &
        stage /= VV_SECOND_STAGE) Return

    fail = 0
    Allocate (fdpdx(1:config%mxatdm), fdpdy(1:config%mxatdm), fdpdz(1:config%mxatdm), Stat=fail)

    If (fail > 0) Then
      Write (message, '(a)') 'dpd_mdvv_forces force array allocation failure'
      Call error(0, message)
    End If

    ! Set tstep derivatives

    hstep = 0.5_wp * tstep
    r_sqrt_tstp = 1.0_wp / Sqrt(tstep)

    ! Initialise DPD virial and stress contributions

    stats%virdpd = 0.0_wp
    stats%strdpdr = 0.0_wp
    stats%strdpdd = 0.0_wp

    ! Initialise local force arrays

    fdpdx = 0.0_wp
    fdpdy = 0.0_wp
    fdpdz = 0.0_wp

    ! Refresh halo velocities

    Call dpd_v_set_halo(domain, config, comm)

    ! Calculate forces

    ! outer loop over atoms

    Do i = 1, config%natms

      ! primary atom type and global index

      ai = config%ltype(i)
      idi = config%ltg(i)

      ! loop over all valid pairs - ignore massless particles

      limit = Merge(neigh%list(0, i), 0, config%weight(i) > 1.0e-6_wp)

      Do k = 1, limit

        ! Secondary atom index, type, and global index

        j = neigh%list(k,i)
        aj = config%ltype(j)
        idj = config%ltg(j)

        ! Calculate r_ij, |r_ij|

        xdif = config%parts(i)%xxx - config%parts(j)%xxx
        ydif = config%parts(i)%yyy - config%parts(j)%yyy
        zdif = config%parts(i)%zzz - config%parts(j)%zzz

        rrr = Sqrt(xdif**2 + ydif**2 + zdif**2)

        ! Calculate forces for valid pairs

        If (rrr < rcut .and. config%weight(j) > 1.0e-6_wp) Then

          ! Calculate v_ij, and (r_ij . v_ij)

          vxdif = config%vxx(i) - config%vxx(j)
          vydif = config%vyy(i) - config%vyy(j)
          vzdif = config%vzz(i) - config%vzz(j)

          rdotv = xdif*vxdif + ydif*vydif + zdif*vzdif

          ! Get mixing type function - key for interaction strength

          If (ai > aj) Then
            key = ai * (ai - 1) / 2 + aj
          Else
            key = aj * (aj - 1) / 2 + ai
          End If

          ! Get gaussian random number with zero mean (held in gauss var)
          ! Global id check ensure same random number for same pair of particles

          If (idi < idj) Then
            Call box_mueller_saru2(seed, idi, idj, nstep, gauss, l_str)
          Else
            Call box_mueller_saru2(seed, idj, idi, nstep, gauss, l_str)
          End If

          ! Screening function, related to the drag and random weight functions:
          ! w_D = scrn**2 * rrr**2
          ! w_R = scrn * rrr

          scrn = (rcut - rrr) / (rrr * rcut)         

          ! Calculate random and drag components

          rgamma = thermo%sigdpd(key) * scrn * gauss * r_sqrt_tstp
          dgamma = - thermo%gamdpd(key) * scrn**2 * rdotv

          ! Total force component gamma_ij such that, summed over all pairs
          ! v_i(t+hstep) = v_i(t) + (hstep / m_i) * (gamma_ij * r_ij)

          gamma = rgamma + dgamma

          ! Update forces

          fdpdx(i) = fdpdx(i) + gamma * xdif
          fdpdy(i) = fdpdy(i) + gamma * ydif
          fdpdz(i) = fdpdz(i) + gamma * zdif
    
#ifndef HALF_HALO
          If (j <= config%natms) Then
#endif /* HALF_HALO */

            fdpdx(j) = fdpdx(j) - gamma * xdif
            fdpdy(j) = fdpdy(j) - gamma * ydif
            fdpdz(j) = fdpdz(j) - gamma * zdif
            
#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          !     Assign stress terms (only when second particle is in subdomain
          !     or has larger global particle index than first to avoid double-counting)

#ifndef HALF_HALO
          If (j <= config%natms .or. idi < idj) Then
#endif /* HALF_HALO */

            stats%strdpdr(1) = stats%strdpdr(1) + rgamma * xdif * xdif ! random stress_xx
            stats%strdpdr(2) = stats%strdpdr(2) + rgamma * ydif * xdif ! random stress_xy
            stats%strdpdr(3) = stats%strdpdr(3) + rgamma * zdif * xdif ! random stress_xz
            stats%strdpdr(5) = stats%strdpdr(5) + rgamma * ydif * ydif ! random stress_yy
            stats%strdpdr(6) = stats%strdpdr(6) + rgamma * ydif * zdif ! random stress_yz
            stats%strdpdr(9) = stats%strdpdr(9) + rgamma * zdif * zdif ! random stress_zz
            
            stats%strdpdd(1) = stats%strdpdd(1) + dgamma * xdif * xdif ! drag stress_xx
            stats%strdpdd(2) = stats%strdpdd(2) + dgamma * ydif * xdif ! drag stress_xy
            stats%strdpdd(3) = stats%strdpdd(3) + dgamma * zdif * xdif ! drag stress_xz
            stats%strdpdd(5) = stats%strdpdd(5) + dgamma * ydif * ydif ! drag stress_yy
            stats%strdpdd(6) = stats%strdpdd(6) + dgamma * ydif * zdif ! drag stress_yz
            stats%strdpdd(9) = stats%strdpdd(9) + dgamma * zdif * zdif ! drag stress_zz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

        End If

      End Do

    End Do

#ifdef HALF_HALO
    ! Share the dpd forces collected in the halo with the parent domains

    Call refresh_halo_dpd_forces(domain, config, config%mxatms, fdpdx, fdpdy, fdpdz, comm)

#endif /* HALF_HALO */


    ! Update velocities

    Do i = 1, config%natms
      If (config%weight(i) > 1.0e-6_wp) Then
        config%parts(i)%fxx = config%parts(i)%fxx + fdpdx(i)
        config%parts(i)%fyy = config%parts(i)%fyy + fdpdy(i)
        config%parts(i)%fzz = config%parts(i)%fzz + fdpdz(i)
      End If
    End Do

    ! Update forces on RBs shared across domains

    If (rigid%share) Then
      Call update_shared_units(config, rigid%list_shared, &
                               rigid%map_shared, SHARED_UNIT_UPDATE_FORCES, domain, comm)
    End If

    ! Symmetrise and globalise random and drag stresses

    stats%strdpdr(4) = stats%strdpdr(2) ! random stress_yx
    stats%strdpdr(7) = stats%strdpdr(3) ! random stress_zx
    stats%strdpdr(8) = stats%strdpdr(6) ! random stress_zy

    Call gsum(comm, stats%strdpdr)

    stats%strdpdd(4) = stats%strdpdd(2) ! drag stress_yx
    stats%strdpdd(7) = stats%strdpdd(3) ! drag stress_zx
    stats%strdpdd(8) = stats%strdpdd(6) ! drag stress_zy

    Call gsum(comm, stats%strdpdd)

    ! Update virial (vir = - Tr(strdpdr) - Tr(strdpdd))

    stats%virdpd = stats%virdpd &
                   - stats%strdpdd(1) - stats%strdpdd(5) - stats%strdpdd(9) &
                   - stats%strdpdr(1) - stats%strdpdr(5) - stats%strdpdr(9)

    Deallocate (fdpdx, fdpdy, fdpdz, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'dpd_mdvv_forces force array deallocation failure'
      Call error(0, message)
    End If

  End Subroutine dpd_mdvv_forces

  Subroutine integrate_rbs_dpd(config, rigid, fdpdx, fdpdy, fdpdz, tstep)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, due to velocity contributions of dpd thermostatting forces -
    ! utilises RB integration methods from nve_1_vv, nve.F90
    !
    ! copyright - daresbury laboratory
    ! author    - k.a.jonathan 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(configuration_type), Intent(InOut) :: config
    Real(Kind=wp),            Intent(InOut) :: fdpdx(:), fdpdy(:), fdpdz(:), tstep

    Character(Len=STR_LEN)         :: message
    Integer                    :: fail(1:14), i, i1, i2, irgd, jrgd, krgd, lrgd, rgdtyp
    Real(Kind=wp)              :: fmx, fmy, fmz, hstep, opx, opy, opz, p0, p1, p2, p3, qt0, &
                                  qt1, qt2, qt3, rot(1:9), tmp, tqx, tqy, tqz, trx, try, &
                                  trz, vpx, vpy, vpz, x(1:1), y(1:1), z(1:1)
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), ggx(:), ggy(:), ggz(:), &
                                  q0t(:), q1t(:), q2t(:), q3t(:), rgdoxt(:), rgdoyt(:), &
                                  rgdozt(:), rgdvxt(:), rgdvyt(:), rgdvzt(:), rgdxxt(:), &
                                  rgdyyt(:), rgdzzt(:), vxt(:), vyt(:), vzt(:), xxt(:), yyt(:), &
                                  zzt(:)

    fail = 0
    Allocate (ggx(1:rigid%max_list * rigid%max_rigid), &
      ggy(1:rigid%max_list * rigid%max_rigid), &
      ggz(1:rigid%max_list * rigid%max_rigid), Stat=fail(7))
    Allocate (xxt(1:config%nlast), yyt(1:config%nlast), zzt(1:config%nlast), Stat=fail(8))
    Allocate (vxt(1:config%nlast), vyt(1:config%nlast), vzt(1:config%nlast), Stat=fail(9))
    Allocate (fxt(1:config%nlast), fyt(1:config%nlast), fzt(1:config%nlast), Stat=fail(10))
    Allocate (q0t(1:rigid%max_rigid), &
      q1t(1:rigid%max_rigid), &
      q2t(1:rigid%max_rigid), &
      q3t(1:rigid%max_rigid), Stat=fail(11))
    Allocate (rgdxxt(1:rigid%max_rigid), &
      rgdyyt(1:rigid%max_rigid), &
      rgdzzt(1:rigid%max_rigid), Stat=fail(12))
    Allocate (rgdvxt(1:rigid%max_rigid), &
      rgdvyt(1:rigid%max_rigid), &
      rgdvzt(1:rigid%max_rigid), Stat=fail(13))
    Allocate (rgdoxt(1:rigid%max_rigid), &
      rgdoyt(1:rigid%max_rigid), &
      rgdozt(1:rigid%max_rigid), Stat=fail(14))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'integrate_rbs_dpd allocation failure'
      Call error(0, message)
    End If

    ! Get the RB particles vectors wrt the RB's COM

    krgd = 0
    Do irgd = 1, rigid%n_types
      rgdtyp = rigid%list(0, irgd)

      ! For all good RBs

      lrgd = rigid%list(-1, irgd)
      If (rigid%frozen(0, rgdtyp) < lrgd) Then
        Do jrgd = 1, lrgd
          krgd = krgd + 1

          i = rigid%index_local(jrgd, irgd) ! local index of particle/site

          ! COM distances

          ggx(krgd) = config%parts(i)%xxx - rigid%xxx(irgd)
          ggy(krgd) = config%parts(i)%yyy - rigid%yyy(irgd)
          ggz(krgd) = config%parts(i)%zzz - rigid%zzz(irgd)
        End Do
      End If
    End Do

    ! minimum image convention for bond vectors

    Call images(config%imcon, config%cell, krgd, ggx, ggy, ggz)

    ! timestep derivatives

    hstep = 0.5_wp * tstep

    ! update velocity of RBs

    krgd = 0
    Do irgd = 1, rigid%n_types
      rgdtyp = rigid%list(0, irgd)

      ! For all good RBs

      lrgd = rigid%list(-1, irgd)
      If (rigid%frozen(0, rgdtyp) < lrgd) Then ! Not that it matters

        ! calculate COM force and torque

        fmx = 0.0_wp; fmy = 0.0_wp; fmz = 0.0_wp
        tqx = 0.0_wp; tqy = 0.0_wp; tqz = 0.0_wp
        Do jrgd = 1, lrgd
          krgd = krgd + 1

          i = rigid%index_local(jrgd, irgd) ! local index of particle/site

          ! If the RB has a frozen particle then no net force

          If (rigid%frozen(0, rgdtyp) == 0) Then
            fmx = fmx + fdpdx(i)
            fmy = fmy + fdpdy(i)
            fmz = fmz + fdpdz(i)
          End If

          tqx = tqx + ggy(krgd) * fdpdz(i) - ggz(krgd) * fdpdy(i)
          tqy = tqy + ggz(krgd) * fdpdx(i) - ggx(krgd) * fdpdz(i)
          tqz = tqz + ggx(krgd) * fdpdy(i) - ggy(krgd) * fdpdx(i)
        End Do

        ! If the RB has 2+ frozen particles (ill=1) the net torque
        ! must align along the axis of rotation

        If (rigid%frozen(0, rgdtyp) > 1) Then
          i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
          i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

          x(1) = config%parts(i1)%xxx - config%parts(i2)%xxx
          y(1) = config%parts(i1)%yyy - config%parts(i2)%yyy
          z(1) = config%parts(i1)%zzz - config%parts(i2)%zzz

          Call images(config%imcon, config%cell, 1, x, y, z)

          tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) / (x(1)**2 + y(1)**2 + z(1)**2)
          tqx = x(1) * tmp
          tqy = y(1) * tmp
          tqz = z(1) * tmp
        End If

        ! current rotation matrix

        Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

        ! calculate torque in principal frame

        trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
        try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
        trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

        ! calculate quaternion torques

        qt0 = 2.0_wp * (-rigid%q1(irgd) * trx - rigid%q2(irgd) * try - rigid%q3(irgd) * trz)
        qt1 = 2.0_wp * (rigid%q0(irgd) * trx - rigid%q3(irgd) * try + rigid%q2(irgd) * trz)
        qt2 = 2.0_wp * (rigid%q3(irgd) * trx + rigid%q0(irgd) * try - rigid%q1(irgd) * trz)
        qt3 = 2.0_wp * (-rigid%q2(irgd) * trx + rigid%q1(irgd) * try + rigid%q0(irgd) * trz)

        ! recover quaternion momenta at half time step

        opx = rigid%oxx(irgd) * rigid%rix(1, rgdtyp)
        opy = rigid%oyy(irgd) * rigid%riy(1, rgdtyp)
        opz = rigid%ozz(irgd) * rigid%riz(1, rgdtyp)

        p0 = 2.0_wp * (-rigid%q1(irgd) * opx - rigid%q2(irgd) * opy - rigid%q3(irgd) * opz)
        p1 = 2.0_wp * (rigid%q0(irgd) * opx - rigid%q3(irgd) * opy + rigid%q2(irgd) * opz)
        p2 = 2.0_wp * (rigid%q3(irgd) * opx + rigid%q0(irgd) * opy - rigid%q1(irgd) * opz)
        p3 = 2.0_wp * (-rigid%q2(irgd) * opx + rigid%q1(irgd) * opy + rigid%q0(irgd) * opz)

        ! update quaternion momenta to full step

        p0 = p0 + hstep * qt0
        p1 = p1 + hstep * qt1
        p2 = p2 + hstep * qt2
        p3 = p3 + hstep * qt3

        ! update RB angular & COM velocities to full step

        opx = 0.5_wp * (-rigid%q1(irgd) * p0 + rigid%q0(irgd) * p1 + rigid%q3(irgd) * p2 - rigid%q2(irgd) * p3)
        opy = 0.5_wp * (-rigid%q2(irgd) * p0 - rigid%q3(irgd) * p1 + rigid%q0(irgd) * p2 + rigid%q1(irgd) * p3)
        opz = 0.5_wp * (-rigid%q3(irgd) * p0 + rigid%q2(irgd) * p1 - rigid%q1(irgd) * p2 + rigid%q0(irgd) * p3)

        rigid%oxx(irgd) = opx * rigid%rix(2, rgdtyp)
        rigid%oyy(irgd) = opy * rigid%riy(2, rgdtyp)
        rigid%ozz(irgd) = opz * rigid%riz(2, rgdtyp)

        tmp = hstep / rigid%weight(0, rgdtyp)
        rigid%vxx(irgd) = rigid%vxx(irgd) + tmp * fmx
        rigid%vyy(irgd) = rigid%vyy(irgd) + tmp * fmy
        rigid%vzz(irgd) = rigid%vzz(irgd) + tmp * fmz

        ! update RB members velocities

        Do jrgd = 1, lrgd
          If (rigid%frozen(jrgd, rgdtyp) == 0) Then
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              x(1) = rigid%x(jrgd, rgdtyp)
              y(1) = rigid%y(jrgd, rgdtyp)
              z(1) = rigid%z(jrgd, rgdtyp)

              ! new atomic velocities in body frame

              vpx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
              vpy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
              vpz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

              ! new atomic velocities in lab frame

              config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
              config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
              config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
            End If
          End If
        End Do

      End If
    End Do

    Deallocate (ggx, ggy, ggz, Stat=fail(7))
    Deallocate (xxt, yyt, zzt, Stat=fail(8))
    Deallocate (vxt, vyt, vzt, Stat=fail(9))
    Deallocate (fxt, fyt, fzt, Stat=fail(10))
    Deallocate (q0t, q1t, q2t, q3t, Stat=fail(11))
    Deallocate (rgdxxt, rgdyyt, rgdzzt, Stat=fail(12))
    Deallocate (rgdvxt, rgdvyt, rgdvzt, Stat=fail(13))
    Deallocate (rgdoxt, rgdoyt, rgdozt, Stat=fail(14))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'integrate_rbs_dpd deallocation failure'
      Call error(0, message)
    End If

  End Subroutine integrate_rbs_dpd

  Subroutine dpd_v_export(mdir, mlast, ixyz0, domain, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export metal density data in domain boundary
    ! regions for halo formation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    ! amended   - i.t.todorov may 2020 - simplification for ixyz0
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: mdir
    Integer,                  Intent(InOut) :: mlast, ixyz0(:)
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)                       :: message
    Integer                                  :: fail, i, iadd, iblock, imove, itmp, ix, iy, iz, j, &
                                                jdnode, jmove, jxyz, kdnode, kx, kxyz, ky, kz, &
                                                limit
    Logical                                  :: safe
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer

    ! Number of transported quantities per particle

    iadd = 4

    fail = 0; limit = iadd * domain%mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'dpd_v_export allocation failure'
      Call error(0, message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock = limit / Merge(2, 1, comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! kxyz - corrected halo reduction factor particles haloing both +&- sides
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0; ky = 0; kz = 0
    If (mdir == -1) Then ! Direction -x
      kx = 1
      jxyz = 1
      kxyz = 3

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir == 1) Then ! Direction +x
      kx = 1
      jxyz = 2
      kxyz = 3

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky = 1
      jxyz = 10
      kxyz = 30

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir == 2) Then ! Direction +y
      ky = 1
      jxyz = 20
      kxyz = 30

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz = 1
      jxyz = 100
      kxyz = 300

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir == 3) Then ! Direction +z
      kz = 1
      jxyz = 200
      kxyz = 300

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(152)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove = 0
    jmove = 0

    ! Initialise array overflow flags

    safe = .true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    Do i = 1, mlast

      ! If the particle is within the remaining 'inverted halo' of this domain

      If (ixyz0(i) > 0) Then

        ! Get the necessary halo indices

        ix = Mod(ixyz0(i), 10) ! [0,1,2,3=1+2]
        iy = Mod(ixyz0(i) - ix, 100) ! [0,10,20,30=10+20]
        iz = Mod(ixyz0(i) - (ix + iy), 1000) ! [0,100,200,300=100+200]

        ! Filter the halo index for the selected direction

        j = ix * kx + iy * ky + iz * kz

        ! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j, 3) == 0)) Then

          ! If safe to proceed

          If ((imove + iadd) <= iblock) Then

            ! pack particle velocity and halo indexing

            buffer(imove + 1) = config%vxx(i)
            buffer(imove + 2) = config%vyy(i)
            buffer(imove + 3) = config%vzz(i)

            ! Use the corrected halo reduction factor when the particle is halo to both +&- sides

            buffer(imove + 4) = Real(ixyz0(i) - Merge(jxyz, kxyz, j == jxyz), wp)

          Else

            safe = .false.

          End If
          imove = imove + iadd

        End If

      End If

    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm, safe)
    If (.not. safe) Then
      itmp = Merge(2, 1, comm%mxnode > 1) * imove
      Call gmax(comm, itmp)
      Call warning(150, Real(itmp, wp), Real(limit, wp), 0.0_wp)
      Call error(154)
    End If

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm, jmove, kdnode, DpdVExp_tag)
      Call gsend(comm, imove, jdnode, DpdVExp_tag)
      Call gwait(comm)
    Else
      jmove = imove
    End If

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe = ((mlast + jmove / iadd) <= config%mxatms)
    Call gcheck(comm, safe)
    If (.not. safe) Then
      itmp = mlast + jmove / iadd
      Call gmax(comm, itmp)
      Call warning(160, Real(itmp, wp), Real(config%mxatms, wp), 0.0_wp)
      Call error(156)
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm, buffer(iblock + 1:iblock + jmove), kdnode, DpdVExp_tag)
      End If
      If (imove > 0) Then
        Call gsend(comm, buffer(1:imove), jdnode, DpdVExp_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data

    j = Merge(iblock, 0, comm%mxnode > 1)
    Do i = 1, jmove / iadd
      mlast = mlast + 1

      ! unpack particle velocity and remaining halo indexing

      config%vxx(mlast) = buffer(j + 1)
      config%vyy(mlast) = buffer(j + 2)
      config%vzz(mlast) = buffer(j + 3)
      ixyz0(mlast) = Nint(buffer(j + 4))

      j = j + iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'dpd_v_export deallocation failure'
      Call error(0, message)
    End If

  End Subroutine dpd_v_export

  Subroutine dpd_v_set_halo(domain, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange exchange of velocity data between
    ! neighbouring domains/nodes
    !
    ! copyright - daresbury laboratory
    ! amended   - i.t.todorov november 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    ! amended   - i.t.todorov may 2020 - simplification for ixyz0
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)   :: message
    Integer              :: fail, mlast
    Integer, Allocatable :: ixyz0(:)
    Logical              :: safe

    fail = 0
    Allocate (ixyz0(1:config%mxatms), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'dpd_v_set_halo allocation failure'
      Call error(0, message)
    End If
    ixyz0(1:config%nlast) = config%ixyz(1:config%nlast)

    ! No halo, start with domain only particles

    mlast = config%natms

    ! exchange atom data in -/+ x directions

    Call dpd_v_export(-1, mlast, ixyz0, domain, config, comm)
    Call dpd_v_export(1, mlast, ixyz0, domain, config, comm)

    ! exchange atom data in -/+ y directions

    Call dpd_v_export(-2, mlast, ixyz0, domain, config, comm)
    Call dpd_v_export(2, mlast, ixyz0, domain, config, comm)

    ! exchange atom data in -/+ z directions

    Call dpd_v_export(-3, mlast, ixyz0, domain, config, comm)
    Call dpd_v_export(3, mlast, ixyz0, domain, config, comm)

    ! check atom totals after data transfer

    safe = (mlast == config%nlast)
    Call gcheck(comm, safe)
    If (.not. safe) Call error(96)

    Deallocate (ixyz0, Stat=fail)
  End Subroutine dpd_v_set_halo

#ifdef HALF_HALO
  Subroutine refresh_halo_dpd_forces(domain, config, mxatms, fxx, fyy, fzz, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to refresh the parent domain forces that are calculated
    ! in the halos of neighbouring domains/nodes
    !
    ! Note: all depends on the ixyz halo array set in set_halo
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2020 - helper routine for 'half-halo' VNL
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( domains_type ),       Intent( In    ) :: domain
    Type( configuration_type ), Intent( In    ) :: config
    Type( comms_type ),         Intent( InOut ) :: comm
    Integer,                    Intent( In    ) :: mxatms
    Real( Kind = wp ),          Intent( InOut ) :: fxx(1:mxatms), fyy(1:mxatms), fzz(1:mxatms)

    ! if one defines npdirB and npdirE indices like in export_atomic_forces,
    ! then one should use 'positive' directions earlier than 'negative' ones
    ! i.e the reverse order of the normal coordinate communications in export_atomic_data

    Call export_dpd_forces( 3, domain, config, mxatms, fxx, fyy, fzz, comm) ! x0, y0, z+
    !Call export_atomic_forces(-3, domain, config, mxatms, fxx, fyy, fzz, comm) ! x0, y0, z+ ! one can skip this with 'half-halo' VNL
    Call export_dpd_forces( 2, domain, config, mxatms, fxx, fyy, fzz, comm) ! x0, y+, z0
    Call export_dpd_forces(-2, domain, config, mxatms, fxx, fyy, fzz, comm) ! x0, y+, z0
    Call export_dpd_forces( 1, domain, config, mxatms, fxx, fyy, fzz, comm) ! x-, y0, z0
    Call export_dpd_forces(-1, domain, config, mxatms, fxx, fyy, fzz, comm) ! x+, y0, z0

  End Subroutine refresh_halo_dpd_forces

  Subroutine export_dpd_forces(mdir, domain, config, mxatms, fxx, fyy, fzz, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export dpd forces in domain boundary regions
    ! for halo refresh
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2020 (helper routine for irreducable VNL)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                      Intent( In    ) :: mdir
    Type( domains_type ),         Intent( In    ) :: domain
    Type( configuration_type ),   Intent( In    ) :: config
    Type( comms_type ),           Intent( InOut ) :: comm
    Integer,                      Intent( In    ) :: mxatms
    Real( Kind = wp ),            Intent( InOut ) :: fxx(1:mxatms), fyy(1:mxatms), fzz(1:mxatms)

#ifdef CHECKS
    Character ( Len = 256 )  :: message
#endif
    Logical           :: safe
    Integer           :: fail,iadd,limit,iblock,npdirB,npdirE,mlast,i,j, &
                         jdnode,kdnode,imove,jmove

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    ! Number of transported quantities per particle

    iadd = 4

    fail=0
    limit=iadd*domain%mxbfxp
    Allocate (buffer(1:limit), Stat=fail)
#ifdef CHECKS
    If (fail > 0) Then
      Write(message,'(a)') 'export_dpd_forces allocation failure for buffer'
      Call error(0,message)
    End If
#endif

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock=limit/Merge(2,1,comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with respect to the direction (mdir)
    ! jdnode - destination (send to), kdnode - source (receive from)
    ! in this case call the routine with 'positive' directions earlier than 'negative' ones: mdir = 3(,-3),2,-2,1,-1

    If (mdir == -1) Then ! Direction -x
      jdnode = domain%map(1)
      kdnode = domain%map(2)

      npdirB = config%ixyzM(1)+1
      npdirE = config%ixyzM(2)
    Else If (mdir == 1) Then ! Direction +x
      jdnode = domain%map(2)
      kdnode = domain%map(1)

      npdirB = config%ixyzM(0)+1
      npdirE = config%ixyzM(1)
    Else If (mdir == -2) Then ! Direction -y
      jdnode = domain%map(3)
      kdnode = domain%map(4)

      npdirB = config%ixyzM(3)+1
      npdirE = config%ixyzM(4)
    Else If (mdir == 2) Then ! Direction +y
      jdnode = domain%map(4)
      kdnode = domain%map(3)

      npdirB = config%ixyzM(2)+1
      npdirE = config%ixyzM(3)
    Else If (mdir == -3) Then ! Direction -z
      jdnode = domain%map(5)
      kdnode = domain%map(6)

      npdirB = config%ixyzM(5)+1
      npdirE = config%ixyzM(6)
    Else If (mdir == 3) Then ! Direction +z
      jdnode = domain%map(6)
      kdnode = domain%map(5)

      npdirB = config%ixyzM(4)+1
      npdirE = config%ixyzM(5)
    Else
      Call error(46)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove=0
    jmove=0

    ! Initialise array overflow flags

    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove=0
    jmove=0

    ! Initialise array overflow flags

    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    If (imove+iadd*(npdirE-npdirB) <= iblock) Then
      Do i=npdirB,npdirE
        buffer(imove+1) = Real(config%ltg(i),wp)
        buffer(imove+2) = fxx(i)
        buffer(imove+3) = fyy(i)
        buffer(imove+4) = fzz(i)

        imove=imove+iadd
      End Do
    Else
      safe=.false.
    End If

    ! Check for array bound overflow (have arrays coped with outgoing data)

#ifdef CHECKS
    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=Merge(2,1,comm%mxnode > 1)*imove
      Call gmax(comm,itmp)
      Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
      Call error(54)
    End If
#endif

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm,jmove,kdnode,DpdVExp_tag-1)
      Call gsend(comm,imove,jdnode,DpdVExp_tag-1)
      Call gwait(comm)
    Else
      jmove=imove
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,DpdVExp_tag-1)
      End If
      If (imove > 0 ) Then
        Call gsend(comm,buffer(1:imove),jdnode,DpdVExp_tag-1)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data

    j=Merge(iblock,0,comm%mxnode > 1)

    Do i=1, jmove/iadd
      mlast = local_index(Int(buffer(j+1)),config%nlast,config%lsi,config%lsa)

      fxx(mlast) = fxx(mlast) + buffer(j+2)
      fyy(mlast) = fyy(mlast) + buffer(j+3)
      fzz(mlast) = fzz(mlast) + buffer(j+4)

      j = j + iadd
    End Do

    Deallocate (buffer, Stat=fail)
#ifdef CHECKS
    If (fail > 0) Then
      Write(message,'(a)') 'export_dpd_forces deallocation failure for buffer'
      Call error(0,message)
    End If
#endif

  End Subroutine export_dpd_forces
#endif

End Module dpd
