Module dpd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global DPD variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2016
  ! contrib   - m.a.seaton august 2020
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
  Use kinds,           Only: wp
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: box_mueller_saru2,&
                             seed_type
  Use rigid_bodies,    Only: rigid_bodies_type
  Use shared_units,    Only: SHARED_UNIT_UPDATE_FORCES,&
                             update_shared_units
  Use statistics,      Only: stats_type
  Use thermostat,      Only: DPD_FIRST_ORDER,&
                             DPD_SECOND_ORDER,&
                             VV_FIRST_STAGE,&
                             VV_SECOND_STAGE,&
                             thermostat_type
#ifdef HALF_HALO
  Use numerics,        Only: local_index
#endif /* HALF_HALO */

  Implicit None

  Private

  Public :: dpd_thermostat

Contains

  Subroutine dpd_thermostat(stage, l_str, rcut, nstep, tstep, stats, thermo, neigh, rigid, domain, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine applying DPD thermostat in a Shardlow's VV manner
    ! using the verlet neighbour neigh%list
    !
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
    Integer                                  :: ai, aj, fail(1:2), i, idi, idj, j, k, key, limit, &
                                                nst_p
    Real(Kind=wp)                            :: dgamma, fix, fiy, fiz, fx, fy, fz, gamma, gauss, &
                                                hstep, rgamma, rrr, rstsq, scl, scrn, strs1, &
                                                strs2, strs3, strs5, strs6, strs9, tmp, tst_p
    Real(Kind=wp), Allocatable, Dimension(:) :: fdpdx, fdpdy, fdpdz, rrt, xxt, yyt, zzt

    If (Any(thermo%key_dpd /= [DPD_FIRST_ORDER, DPD_SECOND_ORDER]) .or. &
        (thermo%key_dpd == DPD_FIRST_ORDER .and. stage == VV_SECOND_STAGE)) Return

    fail = 0
    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail(1))
#ifdef HALF_HALO
    Allocate (fdpdx(1:config%mxatms), fdpdy(1:config%mxatms), fdpdz(1:config%mxatms), Stat=fail(2))
#else
    Allocate (fdpdx(1:config%mxatdm), fdpdy(1:config%mxatdm), fdpdz(1:config%mxatdm), Stat=fail(2))
#endif
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'dpd_thermostat allocation failure'
      Call error(0, message)
    End If

    ! set tstep and nstep wrt to order of splitting

    If (thermo%key_dpd == DPD_FIRST_ORDER) Then
      nst_p = nstep
      tst_p = tstep
    Else
      If (stage == VV_FIRST_STAGE) Then
        nst_p = nstep
      Else ! If (stage == VV_SECOND_STAGE) Then
        nst_p = -nstep
      End If
      tst_p = 0.5_wp * tstep
    End If

    ! Set tstep derivatives

    hstep = 0.5_wp * tst_p
    rstsq = 1.0_wp / Sqrt(tst_p)

    ! initialise DPD virial and stress contributions

    If (stage == VV_FIRST_STAGE) Then
      stats%virdpd = 0.0_wp
      stats%strdpd = 0.0_wp
    End If

    ! FIRST PASS

    ! Initialise forces

    fdpdx = 0.0_wp
    fdpdy = 0.0_wp
    fdpdz = 0.0_wp

    ! Refresh halo velocities

    Call dpd_v_set_halo(domain, config, comm)

    ! outer loop over atoms

    Do i = 1, config%natms

      ! Get neigh%list limit

      limit = Merge(neigh%list(0, i), 0, config%weight(i) > 1.0e-6_wp)

      ! calculate interatomic distances

      Do k = 1, limit
        j = neigh%list(k, i)

        xxt(k) = config%parts(i)%xxx - config%parts(j)%xxx
        yyt(k) = config%parts(i)%yyy - config%parts(j)%yyy
        zzt(k) = config%parts(i)%zzz - config%parts(j)%zzz
      End Do

      ! square of distances

      Do k = 1, limit
        rrt(k) = Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
      End Do

      ! initialise stress tensor accumulators

      strs1 = 0.0_wp
      strs2 = 0.0_wp
      strs3 = 0.0_wp
      strs5 = 0.0_wp
      strs6 = 0.0_wp
      strs9 = 0.0_wp

      ! global identity and atomic type of i

      idi = config%ltg(i)
      ai = config%ltype(i)

      ! load forces

      fix = fdpdx(i)
      fiy = fdpdy(i)
      fiz = fdpdz(i)

      ! start of primary loop for forces evaluation

      Do k = 1, limit

        ! secondary atomic index

        j = neigh%list(k, i)

        ! interatomic distance

        rrr = rrt(k)

        ! validity of thermalisation

        If (rrr < rcut .and. config%weight(j) > 1.0e-6_wp) Then

          ! secondary atomic type and global index

          aj = config%ltype(j)
          idj = config%ltg(j)

          ! Get gaussian random number with zero mean

          Call box_mueller_saru2(seed, idi, idj, nst_p, gauss, l_str)

          ! screening function

          scrn = (rcut - rrr) / (rrr * rcut)

          ! Get mixing type function

          If (ai > aj) Then
            key = ai * (ai - 1) / 2 + aj
          Else
            key = aj * (aj - 1) / 2 + ai
          End If

          ! Calculate force component

          rgamma = thermo%sigdpd(key) * scrn * gauss * rstsq

          tmp = thermo%gamdpd(key) * (scrn**2)
          dgamma = -tmp * (xxt(k) * (config%vxx(i) - config%vxx(j)) + &
                           yyt(k) * (config%vyy(i) - config%vyy(j)) + zzt(k) * (config%vzz(i) - config%vzz(j)))

          gamma = rgamma + dgamma

          ! calculate forces

          fx = gamma * xxt(k)
          fy = gamma * yyt(k)
          fz = gamma * zzt(k)

          fix = fix + fx
          fiy = fiy + fy
          fiz = fiz + fz

#ifndef HALF_HALO
          If (j <= config%natms) Then
#endif /* HALF_HALO */

            fdpdx(j) = fdpdx(j) - fx
            fdpdy(j) = fdpdy(j) - fy
            fdpdz(j) = fdpdz(j) - fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (j <= config%natms .or. idi < idj) Then
#endif /* HALF_HALO */

            ! add virial

            stats%virdpd = stats%virdpd - gamma * rrr * rrr

            ! add stress tensor

            strs1 = strs1 + xxt(k) * fx
            strs2 = strs2 + xxt(k) * fy
            strs3 = strs3 + xxt(k) * fz
            strs5 = strs5 + yyt(k) * fy
            strs6 = strs6 + yyt(k) * fz
            strs9 = strs9 + zzt(k) * fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

        End If

      End Do

      ! load back forces

      fdpdx(i) = fix
      fdpdy(i) = fiy
      fdpdz(i) = fiz

      ! complete stress tensor

      stats%strdpd(1) = stats%strdpd(1) + strs1
      stats%strdpd(2) = stats%strdpd(2) + strs2
      stats%strdpd(3) = stats%strdpd(3) + strs3
      stats%strdpd(4) = stats%strdpd(4) + strs2
      stats%strdpd(5) = stats%strdpd(5) + strs5
      stats%strdpd(6) = stats%strdpd(6) + strs6
      stats%strdpd(7) = stats%strdpd(7) + strs3
      stats%strdpd(8) = stats%strdpd(8) + strs6
      stats%strdpd(9) = stats%strdpd(9) + strs9

    End Do

#ifdef HALF_HALO
    ! Share the dpd forces collected in the halo with the parent domains

    Call refresh_halo_dpd_forces(domain, config, config%mxatms, fdpdx, fdpdy, fdpdz, comm)

#endif /* HALF_HALO */
    ! Update velocities or add to conservative forces

    Do i = 1, config%natms
      If (config%lfree(i) == 0) Then
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * fdpdx(i)
          config%vyy(i) = config%vyy(i) + tmp * fdpdy(i)
          config%vzz(i) = config%vzz(i) + tmp * fdpdz(i)
        End If
      Else ! a RB member
        config%parts(i)%fxx = config%parts(i)%fxx + fdpdx(i)
        config%parts(i)%fyy = config%parts(i)%fyy + fdpdy(i)
        config%parts(i)%fzz = config%parts(i)%fzz + fdpdz(i)
      End If
    End Do

    ! SECOND PASS

    ! Refresh halo velocities

    Call dpd_v_set_halo(domain, config, comm)

    ! Initialise forces

    fdpdx = 0.0_wp
    fdpdy = 0.0_wp
    fdpdz = 0.0_wp

    ! outer loop over atoms

    Do i = 1, config%natms

      ! Get neigh%list limit

      limit = Merge(neigh%list(0, i), 0, config%weight(i) > 1.0e-6_wp)

      ! calculate interatomic distances

      Do k = 1, limit
        j = neigh%list(k, i)

        xxt(k) = config%parts(i)%xxx - config%parts(j)%xxx
        yyt(k) = config%parts(i)%yyy - config%parts(j)%yyy
        zzt(k) = config%parts(i)%zzz - config%parts(j)%zzz
      End Do

      ! square of distances

      Do k = 1, limit
        rrt(k) = Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
      End Do

      ! initialise stress tensor accumulators

      strs1 = 0.0_wp
      strs2 = 0.0_wp
      strs3 = 0.0_wp
      strs5 = 0.0_wp
      strs6 = 0.0_wp
      strs9 = 0.0_wp

      ! global identity and atomic type of i

      idi = config%ltg(i)
      ai = config%ltype(i)

      ! load forces

      fix = fdpdx(i)
      fiy = fdpdy(i)
      fiz = fdpdz(i)

      ! start of primary loop for forces evaluation

      Do k = 1, limit

        ! secondary atomic index

        j = neigh%list(k, i)

        ! interatomic distance

        rrr = rrt(k)

        ! validity of thermalisation

        If (rrr < rcut .and. config%weight(j) > 1.0e-6_wp) Then

          ! secondary atomic type and global index

          aj = config%ltype(j)
          idj = config%ltg(j)

          ! Get gaussian random number with zero mean

          Call box_mueller_saru2(seed, idi, idj, nst_p, gauss, l_str)

          ! screening function

          scrn = (rcut - rrr) / (rrr * rcut)

          ! Get mixing type function

          If (ai > aj) Then
            key = ai * (ai - 1) / 2 + aj
          Else
            key = aj * (aj - 1) / 2 + ai
          End If

          ! Calculate force component

          rgamma = thermo%sigdpd(key) * scrn * gauss * rstsq

          tmp = thermo%gamdpd(key) * (scrn**2)
          scl = tmp / (1.0_wp + tmp * tst_p)
          dgamma = -tmp * (xxt(k) * (config%vxx(i) - config%vxx(j)) + &
                           yyt(k) * (config%vyy(i) - config%vyy(j)) + zzt(k) * (config%vzz(i) - config%vzz(j)))

          gamma = rgamma + scl * (dgamma - rgamma)

          ! calculate forces

          fx = gamma * xxt(k)
          fy = gamma * yyt(k)
          fz = gamma * zzt(k)

          fix = fix + fx
          fiy = fiy + fy
          fiz = fiz + fz

#ifndef HALF_HALO
          If (j <= config%natms) Then
#endif /* HALF_HALO */

            fdpdx(j) = fdpdx(j) - fx
            fdpdy(j) = fdpdy(j) - fy
            fdpdz(j) = fdpdz(j) - fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (j <= config%natms .or. idi < idj) Then
#endif /* HALF_HALO */

            ! add virial

            stats%virdpd = stats%virdpd - gamma * rrr * rrr

            ! add stress tensor

            strs1 = strs1 + xxt(k) * fx
            strs2 = strs2 + xxt(k) * fy
            strs3 = strs3 + xxt(k) * fz
            strs5 = strs5 + yyt(k) * fy
            strs6 = strs6 + yyt(k) * fz
            strs9 = strs9 + zzt(k) * fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

        End If

      End Do

      ! load back forces

      fdpdx(i) = fix
      fdpdy(i) = fiy
      fdpdz(i) = fiz

      ! complete stress tensor

      stats%strdpd(1) = stats%strdpd(1) + strs1
      stats%strdpd(2) = stats%strdpd(2) + strs2
      stats%strdpd(3) = stats%strdpd(3) + strs3
      stats%strdpd(4) = stats%strdpd(4) + strs2
      stats%strdpd(5) = stats%strdpd(5) + strs5
      stats%strdpd(6) = stats%strdpd(6) + strs6
      stats%strdpd(7) = stats%strdpd(7) + strs3
      stats%strdpd(8) = stats%strdpd(8) + strs6
      stats%strdpd(9) = stats%strdpd(9) + strs9

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
      Else ! a RB member
        config%parts(i)%fxx = config%parts(i)%fxx + fdpdx(i)
        config%parts(i)%fyy = config%parts(i)%fyy + fdpdy(i)
        config%parts(i)%fzz = config%parts(i)%fzz + fdpdz(i)
      End If
    End Do

    ! Update forces on RBs

    If (rigid%share) Then
      Call update_shared_units(config, rigid%list_shared, &
                               rigid%map_shared, SHARED_UNIT_UPDATE_FORCES, domain, comm)
    End If

    ! globalise stats%virdpd

    Call gsum(comm, stats%virdpd)

    Deallocate (xxt, yyt, zzt, rrt, Stat=fail(1))
    Deallocate (fdpdx, fdpdy, fdpdz, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'dpd_thermostat deallocation failure'
      Call error(0, message)
    End If

  End Subroutine dpd_thermostat

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

    Character(Len=256)                       :: message
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

    Character(Len=256)   :: message
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
