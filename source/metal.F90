Module metal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global metal interaction variables and
  ! arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov december 2014
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: MetLdExp_tag,&
                             comms_type,&
                             gcheck,&
                             girecv,&
                             gmax,&
                             gsend,&
                             gsum,&
                             gwait
  Use configuration,   Only: configuration_type
  Use constants,       Only: engunit,&
                             fourpi,&
                             sqrpi,&
                             twopi,&
                             zero_plus
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use filename,        Only: FILE_TABEAM, &
                             file_type
  Use kinds,           Only: wi,&
                             wp
  Use neighbours,      Only: neighbours_type
  Use parse,           Only: get_line,&
                             get_word,&
                             lower_case,&
                             word_2_real
  Use site,            Only: site_type

  Implicit None

  Private

  !> Type to contain metal interaction variables
  Type, Public :: metal_type
    Private

    !> Metal potential cut off
    Real(Kind=wp), Public                                  :: rcut
    !> Direct calculations switch
    Logical, Public                                        :: l_direct = .false.
    !> Embedding over Sqrt(rho) but over rho switch
    Logical, Public                                        :: l_emb = .false.
    !> 2B(EAM or EEAM) switch
    Logical, Public                                        :: l_2b = .false.
    !> Number of different metal interactions
    Integer(Kind=wi), Public                               :: n_potentials = 0
    !> - 0 = no TABEAM
    !> - 1 = EAM
    !> - 2 = EEAM
    !> - 3 = 2BEAM
    !> - 4 = 2BEEAM
    Integer(Kind=wi), Public                               :: tab = -1
    Integer(Kind=wi), Allocatable, Public                  :: list(:), ltp(:)
    Real(Kind=wp), Allocatable, Public                     :: prm(:, :)
    !> Energy long range correction
    Real(Kind=wp), Allocatable, Public                     :: elrc(:)
    !> Virial long range correction
    Real(Kind=wp), Allocatable, Public                     :: vlrc(:)
    !> Possible tabulated calculation arrays
    Real(Kind=wp), Allocatable, Dimension(:, :, :), Public :: vmet, dmet, dmes, &
                                                              fmet, fmes
    ! Atomic density [reused as embedding derivative(s)] helper array(s)
    Real(Kind=wp), Allocatable, Dimension(:), Public       :: rho, rhs
    !> Maximum number of metal interactions
    Integer(Kind=wi), Public                               :: max_metal
    !> Maximum number of metal interaction parameters
    Integer(Kind=wi), Public                               :: max_param
    !> Maximum number of grid points
    Integer(Kind=wi), Public                               :: maxgrid
    Integer(Kind=wi), Public                               :: max_med
    Integer(Kind=wi), Public                               :: max_mds
    ! Many-body perturbation potential error function and derivative arrays
    Real(Kind=wp), Allocatable, Dimension(:), Public       :: merf, mfer
    Logical                                                :: newjob = .true.

  Contains
    Private

    Procedure, Public :: init => allocate_metal_arrays
    Procedure, Public :: init_table => allocate_metal_table_arrays
    Final             :: cleanup
  End Type metal_Type

  Public :: metal_generate, metal_generate_erf, metal_table_read, &
            metal_lrc, metal_forces, metal_ld_compute, erfgen_met

Contains

  Subroutine allocate_metal_arrays(met, mxatms, mxatyp)
    Class(metal_type)               :: met
    Integer(Kind=wi), Intent(In   ) :: mxatms, mxatyp

    Integer, Dimension(1:7) :: fail

    If (met%tab == 3 .or. met%tab == 4) met%l_2b = .true.

    fail = 0

    Allocate (met%list(1:met%max_metal), stat=fail(1))
    Allocate (met%ltp(1:met%max_metal), stat=fail(2))
    Allocate (met%prm(1:met%max_param, 1:met%max_metal), stat=fail(3))
    Allocate (met%rho(1:Merge(mxatms, 0, met%max_metal > 0)), stat=fail(4))
    ! the new S-band density
    If (met%l_2b) Then
      Allocate (met%rhs(1:Merge(mxatms, 0, met%max_metal > 0)), stat=fail(5))
    End If
    Allocate (met%elrc(0:mxatyp), stat=fail(6))
    Allocate (met%vlrc(0:mxatyp), stat=fail(7))

    If (Any(fail > 0)) Call error(1023)

    met%list = 0
    met%ltp = 0

    met%prm = 0.0_wp

    ! met%rho and met%rhs get initialised in metal_ld_compute!!!

    met%elrc = 0.0_wp
    met%vlrc = 0.0_wp
  End Subroutine allocate_metal_arrays

  Subroutine allocate_metal_table_arrays(met, mxatyp)
    Class(metal_type)               :: met
    Integer(Kind=wi), Intent(In   ) :: mxatyp

    Integer, Dimension(1:5) :: fail

    fail = 0

    Allocate (met%vmet(1:met%maxgrid, 1:met%max_metal, 1:2), stat=fail(1))
    Allocate (met%dmet(1:met%maxgrid, 1:met%max_med, 1:2), stat=fail(2))
    Allocate (met%fmet(1:met%maxgrid, 1:mxatyp, 1:2), stat=fail(3))
    If (met%tab == 3 .or. met%tab == 4) Then ! the new S-band density and embedding
      Allocate (met%dmes(1:met%maxgrid, 1:met%max_mds, 1:2), stat=fail(4))
      Allocate (met%fmes(1:met%maxgrid, 1:mxatyp, 1:2), stat=fail(5))
    End If

    If (Any(fail > 0)) Call error(1069)

    met%vmet = 0.0_wp
    met%dmet = 0.0_wp
    met%fmet = 0.0_wp
    If (met%tab == 3 .or. met%tab == 4) Then ! the new S-band density and embedding
      met%dmes = 0.0_wp
      met%fmes = 0.0_wp
    End If
  End Subroutine allocate_metal_table_arrays

  Subroutine allocate_metal_erf_arrays(met)
    Type(metal_type), Intent(InOut) :: met

    Integer :: fail

    fail = 0

    Allocate (met%merf(1:met%maxgrid), met%mfer(1:met%maxgrid), Stat=fail)

    If (fail > 0) Call error(1082)

    met%merf = 0.0_wp
    met%mfer = 0.0_wp
  End Subroutine allocate_metal_erf_arrays

  Subroutine metal_forces &
    (iatm, xxt, yyt, zzt, rrt, engmet, virmet, stress, safe, ntype_atom, met, neigh, config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating metal energy and force terms
    ! for EAM and FST interactions using verlet neighbour list
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith & i.t.todorov november 2016
    ! contrib   - r.davidchak (eeam) june 2012
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                    Intent(In   ) :: iatm
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(metal_type),                           Intent(InOut) :: met
    Integer(Kind=wi),                           Intent(In   ) :: ntype_atom
    Logical,                                    Intent(InOut) :: safe
    Real(Kind=wp), Dimension(1:9),              Intent(InOut) :: stress
    Real(Kind=wp),                              Intent(  Out) :: virmet, engmet
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: rrt, zzt, yyt, xxt
    Type(configuration_type),                   Intent(InOut) :: config

    Integer       :: ai, aj, idi, jatm, k0, k1, k2, key, keypot, ki, kj, kmn, kmx, l, ld, m, mmm, &
                     nnn
    Real(Kind=wp) :: aaa, bbb, bet, cc0, cc1, cc2, cc3, cc4, ccc, cut1, cut2, ddd, eng, eps, fix, &
                     fiy, fiz, fx, fy, fz, gam1, gam2, gamm2s, gamm3s, gamma, gamma1, gamma2, &
                     gamma3, gk0, gk1, gk2, mmmr, nnnr, ppd, ppp, qqq, rdr, rr0, rr1, rrr, rsq, &
                     sig, strs1, strs2, strs3, strs5, strs6, strs9, t1, t2, t3, t4, vk0, vk1, vk2

    ! initialise potential energy and virial

    engmet = 0.0_wp
    virmet = 0.0_wp

    ! initialise stress tensor accumulators

    strs1 = 0.0_wp
    strs2 = 0.0_wp
    strs3 = 0.0_wp
    strs5 = 0.0_wp
    strs6 = 0.0_wp
    strs9 = 0.0_wp

    ! global identity and type of atom iatm

    idi = config%ltg(iatm)
    ai = config%ltype(iatm)

    ! load forces

    fix = config%parts(iatm)%fxx
    fiy = config%parts(iatm)%fyy
    fiz = config%parts(iatm)%fzz

    ! start of primary loop for forces evaluation

    Do m = 1, neigh%list(0, iatm)

      ! atomic and potential function indices

      jatm = neigh%list(m, iatm)
      aj = config%ltype(jatm)

      If (met%tab == 1 .or. met%tab == 3) Then ! EAM & 2BEAM
        ki = ai
        kj = aj
      Else If (met%tab == 2 .or. met%tab == 4) Then ! EEAM & 2BEEAM
        ki = (aj - 1) * ntype_atom + ai ! aj-ai
        kj = (ai - 1) * ntype_atom + aj ! ai-aj
      End If

      If (ai > aj) Then
        key = ai * (ai - 1) / 2 + aj
      Else
        key = aj * (aj - 1) / 2 + ai
      End If

      k0 = met%list(key)

      If (met%l_direct) Then
        k1 = Max(ai, aj)
        k2 = Min(ai, aj)

        kmx = k1 * (k1 + 1) / 2
        kmn = k2 * (k2 + 1) / 2

        k1 = met%list(kmx)
        k2 = met%list(kmn)
      End If

      ! interatomic distance

      rrr = rrt(m)

      ! truncation and validity of metal interaction

      keypot = met%ltp(k0)
      If (keypot >= 0 .and. rrr <= met%rcut) Then

        ! Squared distance

        rsq = rrr**2

        ! Zero energy and force components

        eng = 0.0_wp
        gamma1 = 0.0_wp
        gamma2 = 0.0_wp
        gamma3 = 0.0_wp
        gamm2s = 0.0_wp
        gamm3s = 0.0_wp
        gamma = 0.0_wp

        If (met%l_direct) Then ! direct calculation (keypot /= 0)

          ! Type of analytic potential

          If (keypot == 1) Then

            ! finnis-sinclair potentials

            cc0 = met%prm(1, k0)
            cc1 = met%prm(2, k0)
            cc2 = met%prm(3, k0)
            ccc = met%prm(4, k0)
            ddd = met%prm(6, k0)
            bet = met%prm(7, k0)
            cut1 = ccc
            cut2 = ddd

            ! calculate pair forces and energies

            If (rrr <= cut1) Then
              gamma1 = -rrr * (2.0_wp * (cc0 + cc1 * rrr + cc2 * rrr**2) * &
                               (rrr - ccc) + (cc1 + 2.0_wp * cc2 * rrr) * (rrr - ccc)**2)

              If (jatm <= config%natms .or. idi < config%ltg(jatm)) &
                eng = (cc0 + cc1 * rrr + cc2 * rrr**2) * (rrr - ccc)**2
            End If

            ! calculate density contributions

            If (rrr <= cut2) &
              gamma2 = -rrr * (2.0_wp * (rrr - ddd) + 3.0_wp * bet * (rrr - ddd)**2 / ddd)

            If (ai == aj) Then
              t1 = met%prm(5, k0)**2
              t2 = t1
            Else
              t1 = met%prm(5, k1)**2
              t2 = met%prm(5, k2)**2
            End If

          Else If (keypot == 2) Then

            ! extended finnis-sinclair potentials

            cc0 = met%prm(1, k0)
            cc1 = met%prm(2, k0)
            cc2 = met%prm(3, k0)
            cc3 = met%prm(4, k0)
            cc4 = met%prm(5, k0)
            ccc = met%prm(6, k0)
            ddd = met%prm(8, k0)
            bbb = met%prm(9, k0)
            cut1 = ccc
            cut2 = ddd

            ! calculate pair forces and energies

            If (rrr <= cut1) Then
              gamma1 = -rrr * (2.0_wp * (cc0 + cc1 * rrr + cc2 * rrr**2 + cc3 * rrr**3 + cc4 * rrr**4) * (rrr - ccc) + &
                               (cc1 + 2.0_wp * cc2 * rrr + 3.0_wp * cc3 * rrr**2 + 4.0_wp * cc4 * rrr**3) * (rrr - ccc)**2)

              If (jatm <= config%natms .or. idi < config%ltg(jatm)) &
                eng = (cc0 + cc1 * rrr + cc2 * rrr**2 + cc3 * rrr**3 + cc4 * rrr**4) * (rrr - ccc)**2
            End If

            ! calculate density contributions

            If (rrr <= cut2) &
              gamma2 = -rrr * (2.0_wp * (rrr - ddd) + 4.0_wp * bbb**2 * (rrr - ddd)**3)

            If (ai == aj) Then
              t1 = met%prm(7, k0)**2
              t2 = t1
            Else
              t1 = met%prm(7, k1)**2
              t2 = met%prm(7, k2)**2
            End If

          Else If (keypot == 3) Then

            ! sutton-chen potentials

            eps = met%prm(1, k0)
            sig = met%prm(2, k0)
            nnn = Nint(met%prm(3, k0)); nnnr = Real(nnn, wp)
            mmm = Nint(met%prm(4, k0)); mmmr = Real(mmm, wp)

            ! calculate pair forces and energies

            gamma1 = nnnr * eps * (sig / rrr)**nnn
            If (jatm <= config%natms .or. idi < config%ltg(jatm)) &
              eng = gamma1 / nnnr

            ! calculate density contributions

            gamma2 = mmmr * (sig / rrr)**mmm

            If (ai == aj) Then
              t1 = (met%prm(1, k0) * met%prm(5, k0))**2
              t2 = t1
            Else
              t1 = (met%prm(1, k1) * met%prm(5, k1))**2
              t2 = (met%prm(1, k2) * met%prm(5, k2))**2
            End If

          Else If (keypot == 4) Then

            ! gupta potentials

            aaa = met%prm(1, k0)
            rr0 = met%prm(2, k0)
            ppp = met%prm(3, k0)
            qqq = met%prm(5, k0)

            cut1 = (rrr - rr0) / rr0
            cut2 = cut1 + 1.0_wp

            ! calculate pair forces and energies

            gamma1 = 2.0_wp * aaa * Exp(-ppp * cut1) * ppp * cut2
            If (jatm <= config%natms .or. idi < config%ltg(jatm)) &
              eng = gamma1 / (ppp * cut2)

            ! calculate density contributions

            gamma2 = 2.0_wp * Exp(-2.0_wp * qqq * cut1) * qqq * cut2

            t1 = met%prm(4, k0)**2
            t2 = t1

          Else If (keypot == 5) Then

            ! many-body perturbation component only potentials

            eps = met%prm(1, k0)
            sig = met%prm(2, k0)
            mmm = Nint(met%prm(3, k0)); mmmr = Real(mmm, wp)

            ! no pair forces and energies

            !             gamma1=0.0_wp
            !              If (jatm <= natms .or. idi < ltg(jatm)) &
            !                 eng = 0.0_wp

            ! calculate density contributions

            ! interpolation parameters

            rdr = 1.0_wp / met%merf(4)
            rr1 = rrr - met%merf(2)
            l = Min(Nint(rr1 * rdr), Nint(met%merf(1)) - 1)
            If (l < 5) Then ! catch unsafe value
              safe = .false.
              l = 6
            End If
            ppp = rr1 * rdr - Real(l, wp)

            ! calculate density using 3-point interpolation

            vk0 = met%merf(l - 1)
            vk1 = met%merf(l)
            vk2 = met%merf(l + 1)

            t1 = vk1 + ppp * (vk1 - vk0)
            t2 = vk1 + ppp * (vk2 - vk1)

            gk0 = met%mfer(l - 1)
            gk1 = met%mfer(l)
            gk2 = met%mfer(l + 1)

            t3 = gk1 + ppp * (gk1 - gk0)
            t4 = gk1 + ppp * (gk2 - gk1)

            If (ppp < 0.0_wp) Then
              gam1 = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
              gam2 = t3 + 0.5_wp * (t4 - t3) * (ppp + 1.0_wp)
            Else If (l == 5) Then
              gam1 = t2
              gam2 = t4
            Else
              gam1 = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
              gam2 = t4 + 0.5_wp * (t4 - t3) * (ppp - 1.0_wp)
            End If

            gamma2 = (sig / rrr**mmm) * (mmmr * gam1 - rrr * gam2)

            If (ai == aj) Then
              t1 = met%prm(1, k0)**2
              t2 = t1
            Else
              t1 = met%prm(1, k1)**2
              t2 = met%prm(1, k2)**2
            End If

          End If

          If (ai > aj) Then
            gamma = (gamma1 - gamma2 * (met%rho(iatm) * t1 + met%rho(jatm) * t2)) / rsq
          Else
            gamma = (gamma1 - gamma2 * (met%rho(iatm) * t2 + met%rho(jatm) * t1)) / rsq
          End If

          fx = gamma * xxt(m)
          fy = gamma * yyt(m)
          fz = gamma * zzt(m)

          fix = fix + fx
          fiy = fiy + fy
          fiz = fiz + fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - fx
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - fy
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! add interaction energy

            engmet = engmet + eng

            ! add virial

            virmet = virmet - gamma * rsq

            ! add stress tensor

            strs1 = strs1 + xxt(m) * fx
            strs2 = strs2 + xxt(m) * fy
            strs3 = strs3 + xxt(m) * fz
            strs5 = strs5 + yyt(m) * fy
            strs6 = strs6 + yyt(m) * fz
            strs9 = strs9 + zzt(m) * fz

          End If

        Else ! tabulated calculation

          ! truncation of potential

          If (Abs(met%vmet(1, k0, 1)) > zero_plus) Then

            ! interpolation parameters

            If (rrr <= met%vmet(3, k0, 1) .or. & ! Next covers the FST density - merge used to avoid table check beyond bound!
                (keypot /= 0 .and. rrr <= met%dmet(3, Merge(k0, 1, keypot /= 0), 1))) Then

              rdr = 1.0_wp / met%vmet(4, k0, 1)
              rr1 = rrr - met%vmet(2, k0, 1)
              l = Min(Nint(rr1 * rdr), Nint(met%vmet(1, k0, 1)) - 1)
              If (l < 5) Then ! catch unsafe value
                safe = .false.
                l = 6
              End If
              ppp = rr1 * rdr - Real(l, wp)

            End If

            ! calculate pair forces using 3-point interpolation

            If (rrr <= met%vmet(3, k0, 1)) Then

              gk0 = met%vmet(l - 1, k0, 2)
              gk1 = met%vmet(l, k0, 2)
              gk2 = met%vmet(l + 1, k0, 2)

              t1 = gk1 + ppp * (gk1 - gk0)
              t2 = gk1 + ppp * (gk2 - gk1)

              If (ppp < 0.0_wp) Then
                gamma1 = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
              Else If (l == 5) Then
                gamma1 = t2
              Else
                gamma1 = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
              End If

              ! calculate interaction energy using 3-point interpolation

              If ((jatm <= config%natms .or. idi < config%ltg(jatm)) .and. keypot /= 5) Then

                vk0 = met%vmet(l - 1, k0, 1)
                vk1 = met%vmet(l, k0, 1)
                vk2 = met%vmet(l + 1, k0, 1)

                t1 = vk1 + ppp * (vk1 - vk0)
                t2 = vk1 + ppp * (vk2 - vk1)

                If (ppp < 0.0_wp) Then
                  eng = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                Else If (l == 5) Then
                  eng = t2
                Else
                  eng = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                End If

              End If

            End If

          End If

          ! calculate embedding forces using 3-point interpolation

          If (keypot == 0) Then ! EAM

            ! contribution from first metal atom identity

            If (Abs(met%dmet(1, kj, 1)) > zero_plus .and. Nint(met%dmet(1, ki, 1)) > 5) Then
              If (rrr <= met%dmet(3, kj, 1)) Then

                ! interpolation parameters

                rdr = 1.0_wp / met%dmet(4, kj, 1)
                rr1 = rrr - met%dmet(2, kj, 1)
                ld = Min(Nint(rr1 * rdr), Nint(met%dmet(1, kj, 1)) - 1)
                If (ld < 5) Then ! catch unsafe value: EAM
                  safe = .false.
                  ld = 6
                End If
                ppd = rr1 * rdr - Real(ld, wp)

                gk0 = met%dmet(ld - 1, kj, 2)
                gk1 = met%dmet(ld, kj, 2)
                gk2 = met%dmet(ld + 1, kj, 2)

                t1 = gk1 + ppd * (gk1 - gk0)
                t2 = gk1 + ppd * (gk2 - gk1)

                If (ppd < 0.0_wp) Then
                  gamma2 = t1 + 0.5_wp * (t2 - t1) * (ppd + 1.0_wp)
                Else If (ld == 5) Then
                  gamma2 = t2
                Else
                  gamma2 = t2 + 0.5_wp * (t2 - t1) * (ppd - 1.0_wp)
                End If

              End If
            End If

            ! Now if we have 2B(EAM & EEAM) then do s-band too

            If (met%l_2b) Then
              If (Abs(met%dmes(1, kj, 1)) > zero_plus .and. Nint(met%dmes(1, ki, 1)) > 5) Then
                If (rrr <= met%dmes(3, kj, 1)) Then

                  ! interpolation parameters

                  rdr = 1.0_wp / met%dmes(4, kj, 1)
                  rr1 = rrr - met%dmes(2, kj, 1)
                  ld = Min(Nint(rr1 * rdr), Nint(met%dmes(1, kj, 1)) - 1)
                  If (ld < 5) Then ! catch unsafe value: EAM
                    safe = .false.
                    ld = 6
                  End If
                  ppd = rr1 * rdr - Real(ld, wp)

                  gk0 = met%dmes(ld - 1, kj, 2)
                  gk1 = met%dmes(ld, kj, 2)
                  gk2 = met%dmes(ld + 1, kj, 2)

                  t1 = gk1 + ppd * (gk1 - gk0)
                  t2 = gk1 + ppd * (gk2 - gk1)

                  If (ppd < 0.0_wp) Then
                    gamm2s = t1 + 0.5_wp * (t2 - t1) * (ppd + 1.0_wp)
                  Else If (ld == 5) Then
                    gamm2s = t2
                  Else
                    gamm2s = t2 + 0.5_wp * (t2 - t1) * (ppd - 1.0_wp)
                  End If

                End If
              End If
            End If

            ! contribution from second metal atom identity

            If (ki == kj) Then

              gamma3 = gamma2
              If (met%l_2b) gamm3s = gamm2s !2B(EAM & EEAM)

            Else

              If (Abs(met%dmet(1, ki, 1)) > zero_plus .and. Nint(met%dmet(1, ki, 1)) > 5) Then
                If (rrr <= met%dmet(3, ki, 1)) Then

                  ! interpolation parameters

                  rdr = 1.0_wp / met%dmet(4, ki, 1)
                  rr1 = rrr - met%dmet(2, ki, 1)
                  ld = Min(Nint(rr1 * rdr), Nint(met%dmet(1, ki, 1)) - 1)
                  If (ld < 5) Then ! catch unsafe value: EAM
                    safe = .false.
                    ld = 6
                  End If
                  ppd = rr1 * rdr - Real(ld, wp)

                  gk0 = met%dmet(ld - 1, ki, 2)
                  gk1 = met%dmet(ld, ki, 2)
                  gk2 = met%dmet(ld + 1, ki, 2)

                  t1 = gk1 + ppd * (gk1 - gk0)
                  t2 = gk1 + ppd * (gk2 - gk1)

                  If (ppd < 0.0_wp) Then
                    gamma3 = t1 + 0.5_wp * (t2 - t1) * (ppd + 1.0_wp)
                  Else If (ld == 5) Then
                    gamma3 = t2
                  Else
                    gamma3 = t2 + 0.5_wp * (t2 - t1) * (ppd - 1.0_wp)
                  End If

                End If
              End If

              If (met%l_2b) Then !2B(EAM & EEAM)
                If (Abs(met%dmes(1, ki, 1)) > zero_plus .and. Nint(met%dmes(1, ki, 1)) > 5) Then
                  If (rrr <= met%dmes(3, ki, 1)) Then

                    ! interpolation parameters

                    rdr = 1.0_wp / met%dmes(4, ki, 1)
                    rr1 = rrr - met%dmes(2, ki, 1)
                    ld = Min(Nint(rr1 * rdr), Nint(met%dmes(1, ki, 1)) - 1)
                    If (ld < 5) Then ! catch unsafe value: EAM
                      safe = .false.
                      ld = 6
                    End If
                    ppd = rr1 * rdr - Real(ld, wp)

                    gk0 = met%dmes(ld - 1, ki, 2)
                    gk1 = met%dmes(ld, ki, 2)
                    gk2 = met%dmes(ld + 1, ki, 2)

                    t1 = gk1 + ppd * (gk1 - gk0)
                    t2 = gk1 + ppd * (gk2 - gk1)

                    If (ppd < 0.0_wp) Then
                      gamm3s = t1 + 0.5_wp * (t2 - t1) * (ppd + 1.0_wp)
                    Else If (ld == 5) Then
                      gamm3s = t2
                    Else
                      gamm3s = t2 + 0.5_wp * (t2 - t1) * (ppd - 1.0_wp)
                    End If

                  End If
                End If
              End If

            End If

            If (.not. met%l_2b) Then
              gamma = (gamma1 + (gamma2 * met%rho(iatm) + gamma3 * met%rho(jatm))) / rsq
            Else !2B(EAM & EEAM)
              gamma = (gamma1 + (gamma2 * met%rho(iatm) + gamma3 * met%rho(jatm)) &
                       + (gamm2s * met%rhs(iatm) + gamm3s * met%rhs(jatm))) / rsq
            End If

          Else ! FST, interpolation parameters are the same for all force arrays

            If (rrr <= met%dmet(3, k0, 1)) Then ! interpolation parameters covered above

              gk0 = met%dmet(l - 1, k0, 2)
              gk1 = met%dmet(l, k0, 2)
              gk2 = met%dmet(l + 1, k0, 2)

              t1 = gk1 + ppp * (gk1 - gk0)
              t2 = gk1 + ppp * (gk2 - gk1)

              If (ppp < 0.0_wp) Then
                gamma2 = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
              Else If (l == 5) Then
                gamma2 = t2
              Else
                gamma2 = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
              End If

            End If

            If (ai > aj) Then
              gamma = (gamma1 - gamma2 * (met%rho(iatm) * met%dmet(1, k0, 2) + met%rho(jatm) * met%dmet(2, k0, 2))) / rsq
            Else
              gamma = (gamma1 - gamma2 * (met%rho(iatm) * met%dmet(2, k0, 2) + met%rho(jatm) * met%dmet(1, k0, 2))) / rsq
            End If

          End If

          fx = gamma * xxt(m)
          fy = gamma * yyt(m)
          fz = gamma * zzt(m)

          fix = fix + fx
          fiy = fiy + fy
          fiz = fiz + fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - fx
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - fy
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! add interaction energy using 3-point interpolation

            engmet = engmet + eng

            ! add virial

            virmet = virmet - gamma * rsq

            ! add stress tensor

            strs1 = strs1 + xxt(m) * fx
            strs2 = strs2 + xxt(m) * fy
            strs3 = strs3 + xxt(m) * fz
            strs5 = strs5 + yyt(m) * fy
            strs6 = strs6 + yyt(m) * fz
            strs9 = strs9 + zzt(m) * fz

          End If

        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = fix
    config%parts(iatm)%fyy = fiy
    config%parts(iatm)%fzz = fiz

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
  End Subroutine metal_forces

  Subroutine metal_ld_compute(engden, virden, stress, ntype_atom, met, neigh, domain, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating local density in metals using
    ! the verlet neighbour list and sutton-chen potentials
    !
    ! Note: Designed to be used as part of two_body_forces
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1998
    ! amended   - i.t.todorov january 2016
    ! contrib   - r.davidchak (eeam) june 2012
    ! contrib   - b.palmer (2band) may 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),                 Intent(  Out) :: engden, virden
    Real(Kind=wp), Dimension(1:9), Intent(InOut) :: stress
    Integer(Kind=wi),              Intent(In   ) :: ntype_atom
    Type(metal_type),              Intent(InOut) :: met
    Type(neighbours_type),         Intent(In   ) :: neigh
    Type(domains_type),            Intent(In   ) :: domain
    Type(configuration_type),      Intent(InOut) :: config
    Type(comms_type),              Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Integer                                  :: fail, i, j, k, k0, l, limit
    Logical                                  :: safe
    Real(Kind=wp)                            :: fk0, fk1, fk2, ppp, rdr, rhosqr, rrr, t1, t2
    Real(Kind=wp), Allocatable, Dimension(:) :: rrt, xxt, yyt, zzt

    ! check on mixing metal types done in read_field

    ! initialise energy and virial accumulators
    safe = .true.

    engden = 0.0_wp
    virden = 0.0_wp

    ! initialise density array

    met%rho = 0.0_wp
    If (met%l_2b) met%rhs = 0.0_wp

    ! All calls below act on met%rho (met%rhs)

    ! calculate local atomic density
    ! outer loop over atoms

    fail = 0
    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'metal_ld_compute allocation failure'
      Call error(0, message)
    End If

    Do i = 1, config%natms
      limit = neigh%list(0, i) ! Get list limit

      ! calculate interatomic distances

      Do k = 1, limit
        j = neigh%list(k, i)

        xxt(k) = config%parts(i)%xxx - config%parts(j)%xxx
        yyt(k) = config%parts(i)%yyy - config%parts(j)%yyy
        zzt(k) = config%parts(i)%zzz - config%parts(j)%zzz
      End Do

      ! periodic boundary conditions not needed by LC construction
      !
      !     Call images(imcon,cell,limit,xxt,yyt,zzt)

      ! square of distances

      Do k = 1, limit
        rrt(k) = Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
      End Do

      ! calculate contributions to local density

      If (met%tab > 0) Then ! EAM contributions
        Call metal_ld_collect_eam(i, rrt, safe, ntype_atom, met, neigh, config)
      Else ! If (met%tab == 0) Then ! FST contributions
        Call metal_ld_collect_fst(i, rrt, safe, met, neigh, config)
      End If
    End Do

    Deallocate (xxt, yyt, zzt, rrt, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a,i0)') 'metal_ld_compute allocation failure'
      Call error(0, message)
    End If

    ! Check safety for densities of EAM and MBPC

    Call gcheck(comm, safe)
    If (.not. safe) Call error(506)

    Do i = 1, config%natms

      ! calculate density terms to energy and virial

      If (met%tab > 0) Then ! EAM potential

        ! potential function index

        k0 = config%ltype(i)

        ! Now start traditional s-band (EAM & EEAM) or d-band for 2B(EAM & EEAM)

        ! validity of potential

        If (Abs(met%fmet(1, k0, 1)) > zero_plus) Then

          ! check for unsafe densities (mind start was shifted)

          If (.not. met%l_emb) Then ! met%fmet over met%rho grid
            rhosqr = met%rho(i)
          Else ! met%fmet over Sqrt(met%rho) grid
            rhosqr = Sqrt(met%rho(i))
          End If
          If (rhosqr >= met%fmet(2, k0, 1) + 5.0_wp * met%fmet(4, k0, 1)) Then
            If (rhosqr <= met%fmet(3, k0, 1)) Then

              ! interpolation parameters

              rdr = 1.0_wp / met%fmet(4, k0, 1)
              rrr = rhosqr - met%fmet(2, k0, 1)
              l = Min(Nint(rrr * rdr), Nint(met%fmet(1, k0, 1)) - 1)
              If (l < 5) Then ! catch unsafe value
                Write (*, *) 'good density range problem: (LTG,RHO) ', config%ltg(i), met%rho(i)
                safe = .false.
                l = 6
              End If
              ppp = rrr * rdr - Real(l, wp)

              ! calculate embedding energy using 3-point interpolation

              fk0 = met%fmet(l - 1, k0, 1)
              fk1 = met%fmet(l, k0, 1)
              fk2 = met%fmet(l + 1, k0, 1)

              t1 = fk1 + ppp * (fk1 - fk0)
              t2 = fk1 + ppp * (fk2 - fk1)

              If (ppp < 0.0_wp) Then
                engden = engden + t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
              Else If (l == 5) Then
                engden = engden + t2
              Else
                engden = engden + t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
              End If

              ! calculate derivative of embedding function wrt density
              ! using 3-point interpolation and STORE/OVERWRITE result in met%rho array

              fk0 = met%fmet(l - 1, k0, 2)
              fk1 = met%fmet(l, k0, 2)
              fk2 = met%fmet(l + 1, k0, 2)

              t1 = fk1 + ppp * (fk1 - fk0)
              t2 = fk1 + ppp * (fk2 - fk1)

              If (ppp < 0.0_wp) Then
                If (.not. met%l_emb) Then ! met%fmet over met%rho grid
                  met%rho(i) = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                Else ! met%fmet over Sqrt(met%rho) grid
                  met%rho(i) = 0.5_wp * (t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)) / rhosqr
                End If
              Else If (l == 5) Then
                If (.not. met%l_emb) Then ! met%fmet over met%rho grid
                  met%rho(i) = t2
                Else ! met%fmet over Sqrt(met%rho) grid
                  met%rho(i) = 0.5_wp * t2 / rhosqr
                End If
              Else
                If (.not. met%l_emb) Then ! met%fmet over met%rho grid
                  met%rho(i) = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                Else ! met%fmet over Sqrt(met%rho) grid
                  met%rho(i) = 0.5_wp * (t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)) / rhosqr
                End If
              End If

            Else ! RLD: assume that met%fmet(met%rho(i) > met%fmet(3,k0,1)) = met%fmet(met%rho(i) = met%fmet(3,k0,1))

              l = Nint(met%fmet(1, k0, 1))

              engden = engden + met%fmet(l, k0, 1)

              met%rho(i) = 0.0_wp

            End If
          Else
            Write (*, *) 'bad density range problem: (LTG,RHO) ', config%ltg(i), met%rho(i)
            safe = .false.
          End If

        End If

        ! Atomic density (met%rho & met%rhs) are overwritten here in order
        ! to be reused as embedding derivative(s) helper array(s)
        ! i.e. hold d_fmet/d_rho for later usage in metal_forces

        ! Now if we have 2B(EAM & EEAM) then do s-band too

        If (met%l_2b) Then

          ! validity of potential

          If (Abs(met%fmes(1, k0, 1)) > zero_plus) Then

            ! check for unsafe densities (mind start was shifted)

            If (.not. met%l_emb) Then ! met%fmes over met%rhs grid
              rhosqr = met%rhs(i)
            Else ! met%fmes over Sqrt(met%rhs) grid
              rhosqr = Sqrt(met%rhs(i))
            End If
            If (rhosqr >= met%fmes(2, k0, 1) + 5.0_wp * met%fmes(4, k0, 1)) Then
              If (rhosqr <= met%fmes(3, k0, 1)) Then

                ! interpolation parameters

                rdr = 1.0_wp / met%fmes(4, k0, 1)
                rrr = rhosqr - met%fmes(2, k0, 1)
                l = Min(Nint(rrr * rdr), Nint(met%fmes(1, k0, 1)) - 1)
                If (l < 5) Then ! catch unsafe value
                  Write (*, *) 'good density range problem: (LTG,RHS) ', config%ltg(i), met%rhs(i)
                  safe = .false.
                  l = 6
                End If
                ppp = rrr * rdr - Real(l, wp)

                ! calculate embedding energy using 3-point interpolation

                fk0 = met%fmes(l - 1, k0, 1)
                fk1 = met%fmes(l, k0, 1)
                fk2 = met%fmes(l + 1, k0, 1)

                t1 = fk1 + ppp * (fk1 - fk0)
                t2 = fk1 + ppp * (fk2 - fk1)

                If (ppp < 0.0_wp) Then
                  engden = engden + t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                Else If (l == 5) Then
                  engden = engden + t2
                Else
                  engden = engden + t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                End If

                ! calculate derivative of embedding function wrt density
                ! using 3-point interpolation and STORE/OVERWRITE result in met%rhs array

                fk0 = met%fmes(l - 1, k0, 2)
                fk1 = met%fmes(l, k0, 2)
                fk2 = met%fmes(l + 1, k0, 2)

                t1 = fk1 + ppp * (fk1 - fk0)
                t2 = fk1 + ppp * (fk2 - fk1)

                If (ppp < 0.0_wp) Then
                  If (.not. met%l_emb) Then ! met%fmes over met%rhs grid
                    met%rhs(i) = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                  Else ! met%fmes over Sqrt(met%rhs) grid
                    met%rhs(i) = 0.5_wp * (t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)) / rhosqr
                  End If
                Else If (l == 5) Then
                  If (.not. met%l_emb) Then ! met%fmes over met%rhs grid
                    met%rhs(i) = t2
                  Else ! met%fmes over Sqrt(met%rhs) grid
                    met%rhs(i) = 0.5_wp * t2 / rhosqr
                  End If
                Else
                  If (.not. met%l_emb) Then ! met%fmes over met%rhs grid
                    met%rhs(i) = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                  Else ! met%fmes over Sqrt(met%rhs) grid
                    met%rhs(i) = 0.5_wp * (t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)) / rhosqr
                  End If
                End If

              Else ! RLD: assume that met%fmes(met%rhs(i) > met%fmes(3,k0,1)) = met%fmes(met%rhs(i) = met%fmes(3,k0,1))

                l = Nint(met%fmes(1, k0, 1))

                engden = engden + met%fmes(l, k0, 1)

                met%rhs(i) = 0.0_wp

              End If
            Else
              Write (message, '(a,2(i0,1x))') 'bad density range problem: (LTG,RHS) ', config%ltg(i), met%rhs(i)
              Call info(message)
              safe = .false.
            End If

          End If

        End If

      Else ! If (met%tab == 0) Then FST of metal potentials

        If (met%rho(i) > zero_plus) Then

          ! calculate analytical square root of (density + lrc to it)

          rhosqr = Sqrt(met%rho(i) + met%elrc(config%ltype(i)))
          engden = engden - rhosqr
          virden = virden + met%vlrc(config%ltype(i)) / rhosqr

          ! store the derivatives of the FST embedding-like function
          ! (with corrected density) in met%rho array

          met%rho(i) = 0.5_wp / rhosqr

        Else If (met%rho(i) < -zero_plus) Then

          ! check for unsafe densities (met%rho was initialised to zero)

          safe = .false.

        End If

      End If

    End Do

    ! Check safety for densities

    Call gcheck(comm, safe)
    If (.not. safe) Call error(507)

    ! virial term (averaged per node)

    Call gsum(comm, virden)
    virden = virden / Real(comm%mxnode, wp)

    ! calculate stress tensor (density contributions are to
    ! diagonal elements only)

    stress(1) = stress(1) - virden / 3.0_wp
    stress(5) = stress(5) - virden / 3.0_wp
    stress(9) = stress(9) - virden / 3.0_wp

    ! obtain atomic densities for outer border regions

    Call metal_ld_set_halo(met, domain, config, comm)
  End Subroutine metal_ld_compute

  Subroutine metal_lrc(met, sites, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to evaluate metal long-range corrections to
    ! pressure and energy in a 3D periodic system
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith june 1995
    ! amended   - i.t.todorov february 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(metal_type),         Intent(InOut) :: met
    Type(site_type),          Intent(In   ) :: sites
    Type(configuration_type), Intent(In   ) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256) :: message, messages(3)
    Integer            :: i, j, k0, k1, k2, keypot, kmet, mmm, nnn
    Real(Kind=wp)      :: aaa, ccc, eee, elrc0, elrc1, elrc2, elrcsum, eps, mmmr, nnnr, ppp, qqq, &
                          rr0, sig, tmp, vlrc0, vlrc1, vlrc2, zet

    ! long-range corrections to energy, pressure and density

    met%elrc = 0.0_wp
    met%vlrc = 0.0_wp
    elrcsum = 0.0_wp

    If (config%imcon /= 0 .and. config%imcon /= 6) Then
      kmet = 0

      Do i = 1, sites%ntype_atom
        Do j = 1, i

          elrc0 = 0.0_wp
          elrc1 = 0.0_wp
          elrc2 = 0.0_wp

          vlrc0 = 0.0_wp
          vlrc1 = 0.0_wp
          vlrc2 = 0.0_wp

          kmet = kmet + 1
          k0 = met%list(kmet)

          keypot = met%ltp(k0)
          If (keypot == 3) Then

            ! sutton-chen potentials

            eps = met%prm(1, k0)
            sig = met%prm(2, k0)
            nnn = Nint(met%prm(3, k0)); nnnr = Real(nnn, wp)
            mmm = Nint(met%prm(4, k0)); mmmr = Real(mmm, wp)
            ccc = met%prm(5, k0)

            elrc0 = eps * sig**3 * (sig / met%rcut)**(nnn - 3) / (nnnr - 3.0_wp)
            vlrc0 = nnnr * elrc0

            ! Self-interaction accounted once, interaction between different species
            ! MUST be accounted twice!!

            If (i /= j) Then
              elrc0 = elrc0 * 2.0_wp
              vlrc0 = vlrc0 * 2.0_wp
            End If

            met%elrc(0) = met%elrc(0) + twopi * config%volm * sites%dens(i) * sites%dens(j) * elrc0
            met%vlrc(0) = met%vlrc(0) - twopi * config%volm * sites%dens(i) * sites%dens(j) * vlrc0

            tmp = sig**3 * (sig / met%rcut)**(mmm - 3) / (mmmr - 3.0_wp)
            If (i == j) Then
              elrc1 = tmp * (eps * ccc)**2
              met%elrc(i) = met%elrc(i) + fourpi * sites%dens(i) * elrc1
              elrcsum = elrcsum + twopi * config%volm * sites%dens(i)**2 * elrc1

              vlrc1 = mmmr * elrc1
              met%vlrc(i) = met%vlrc(i) + twopi * sites%dens(i) * vlrc1
            Else
              k1 = met%list((i * (i + 1)) / 2)
              k2 = met%list((j * (j + 1)) / 2)

              elrc1 = tmp * (met%prm(1, k1) * met%prm(5, k1))**2
              elrc2 = tmp * (met%prm(1, k2) * met%prm(5, k2))**2
              met%elrc(i) = met%elrc(i) + fourpi * sites%dens(j) * elrc1
              met%elrc(j) = met%elrc(j) + fourpi * sites%dens(i) * elrc2
              elrcsum = elrcsum + twopi * config%volm * sites%dens(i) * sites%dens(j) * (elrc1 + elrc2)

              vlrc1 = mmmr * elrc1
              vlrc2 = mmmr * elrc2
              met%vlrc(i) = met%vlrc(i) + twopi * sites%dens(j) * vlrc1
              met%vlrc(j) = met%vlrc(j) + twopi * sites%dens(i) * vlrc2
            End If

          Else If (keypot == 4) Then

            ! gupta potentials

            aaa = met%prm(1, k0)
            rr0 = met%prm(2, k0)
            ppp = met%prm(3, k0)
            zet = met%prm(4, k0)
            qqq = met%prm(5, k0)
            eee = Exp(-ppp * (met%rcut - rr0) / rr0)

            elrc0 = 2.0_wp * aaa * (rr0 / ppp) * (met%rcut**2 + 2.0_wp * met%rcut * (rr0 / ppp) + 2.0_wp * (rr0 / ppp)**2) * eee
            vlrc0 = 2.0_wp * aaa * met%rcut**3 * eee + 3.0_wp * elrc0

            ! Self-interaction accounted once, interaction between different species
            ! MUST be accounted twice!!

            If (i /= j) Then
              elrc0 = elrc0 * 2.0_wp
              vlrc0 = vlrc0 * 2.0_wp
            End If

            met%elrc(0) = met%elrc(0) + twopi * config%volm * sites%dens(i) * sites%dens(j) * elrc0
            met%vlrc(0) = met%vlrc(0) - twopi * config%volm * sites%dens(i) * sites%dens(j) * vlrc0

            eee = Exp(-2.0_wp * qqq * (met%rcut - rr0) / rr0)

            If (i == j) Then
              elrc1 = (met%rcut**2 + 2.0_wp * met%rcut * (0.5_wp * rr0 / qqq) + &
                       2.0_wp * (0.5_wp * rr0 / qqq)**2) * (0.5_wp * rr0 / qqq) * eee * zet**2
              met%elrc(i) = met%elrc(i) + fourpi * sites%dens(i) * elrc1
              elrcsum = elrcsum + twopi * config%volm * sites%dens(i)**2 * elrc1

              vlrc1 = (met%rcut**3 + 3.0_wp * met%rcut**2 * (0.5_wp * rr0 / qqq) + &
                       6.0_wp * met%rcut * (0.5_wp * rr0 / qqq)**2 + (0.5_wp * rr0 / qqq)**3) * eee * zet**2
              met%vlrc(i) = met%vlrc(i) + twopi * sites%dens(i) * vlrc1
            Else
              elrc1 = (met%rcut**2 + 2.0_wp * met%rcut * (0.5_wp * rr0 / qqq) + &
                       2.0_wp * (0.5_wp * rr0 / qqq)**2) * (0.5_wp * rr0 / qqq) * eee * zet**2
              elrc2 = elrc2
              met%elrc(i) = met%elrc(i) + fourpi * sites%dens(j) * elrc1
              met%elrc(j) = met%elrc(j) + fourpi * sites%dens(i) * elrc2
              elrcsum = elrcsum + twopi * config%volm * sites%dens(i) * sites%dens(j) * (elrc1 + elrc2)

              vlrc1 = (met%rcut**3 + 3.0_wp * met%rcut**2 * (0.5_wp * rr0 / qqq) + &
                       6.0_wp * met%rcut * (0.5_wp * rr0 / qqq)**2 + (0.5_wp * rr0 / qqq)**3) * eee * zet**2
              vlrc2 = vlrc1
              met%vlrc(i) = met%vlrc(i) + twopi * sites%dens(j) * vlrc1
              met%vlrc(j) = met%vlrc(j) + twopi * sites%dens(i) * vlrc2
            End If

          Else If (keypot == 5) Then

            ! many-body perturbation component only potentials

            eps = met%prm(1, k0)
            sig = met%prm(2, k0)
            mmm = Nint(met%prm(3, k0)); mmmr = Real(mmm, wp)

            ! No pairwise contributions for mbpc potentials!!!

            !              elrc0=0.0_wp
            !              vlrc0=0.0_wp

            ! Self-interaction accounted once, interaction between different species
            ! MUST be accounted twice!!

            !              If (i /= j) Then
            !                 elrc0 = elrc0*2.0_wp
            !                 vlrc0 = vlrc0*2.0_wp
            !              End If

            !              met%elrc(0) = met%elrc(0) + twopi*volm*sites%dens(i)*sites%dens(j)*elrc0
            !              met%vlrc(0) = met%vlrc(0) - twopi*volm*sites%dens(i)*sites%dens(j)*vlrc0

            tmp = sig / ((mmmr - 3.0_wp) * met%rcut**(mmm - 3))
            If (i == j) Then
              elrc1 = tmp * eps**2
              met%elrc(i) = met%elrc(i) + fourpi * sites%dens(i) * elrc1
              elrcsum = elrcsum + twopi * config%volm * sites%dens(i)**2 * elrc1

              vlrc1 = mmmr * elrc1
              met%vlrc(i) = met%vlrc(i) + twopi * sites%dens(i) * vlrc1
            Else
              k1 = met%list((i * (i + 1)) / 2)
              k2 = met%list((j * (j + 1)) / 2)

              elrc1 = tmp * met%prm(1, k1)**2
              elrc2 = tmp * met%prm(1, k2)**2
              met%elrc(i) = met%elrc(i) + fourpi * sites%dens(j) * elrc1
              met%elrc(j) = met%elrc(j) + fourpi * sites%dens(i) * elrc2
              elrcsum = elrcsum + twopi * config%volm * sites%dens(i) * sites%dens(j) * (elrc1 + elrc2)

              vlrc1 = mmmr * elrc1
              vlrc2 = mmmr * elrc2
              met%vlrc(i) = met%vlrc(i) + twopi * sites%dens(j) * vlrc1
              met%vlrc(j) = met%vlrc(j) + twopi * sites%dens(i) * vlrc2
            End If

          End If

        End Do
      End Do
    End If

    If (met%newjob) Then
      met%newjob = .false.

      Write (messages(1), '(a,1p,e15.6)') &
        'long-range correction to metal energy ', met%elrc(0) / engunit
      Write (messages(2), '(a,1p,e15.6)') &
        'lr correction for metal atom density ', elrcsum / engunit**2
      Write (messages(3), '(a,1p,e15.6)') &
        '1st partial lr correction to metal virial', met%vlrc(0) / engunit
      Call info(messages, 3, .true.)

      Call info('density dependent energy and virial corrections:', .true.)
      If (comm%idnode == 0) Then
        Do i = 1, sites%ntype_atom
          kmet = met%list((i * (i + 1)) / 2)
          If (met%list(kmet) > 0) Then
            Write (message, "(2x,a8,1p,2e15.6)") &
              sites%unique_atom(i), met%elrc(i) / engunit, met%vlrc(i) / engunit
            Call info(message, .true.)
          End If
        End Do
      End If
    End If
  End Subroutine metal_lrc

  Subroutine metal_table_read(l_top, met, sites, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reading potential energy and force arrays
    ! from TABEAM file (for metal EAM & EEAM forces only)
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith march 2006
    ! amended   - i.t.todorov march 2016
    ! contrib   - r.davidchak (eeam) june 2012
    ! contrib   - b.palmer (2band) may 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,          Intent(In   ) :: l_top
    Type(metal_type), Intent(InOut) :: met
    Type(site_type),  Intent(In   ) :: sites
    Type(file_type),  Intent(InOut) :: files(:)
    Type(comms_type), Intent(InOut) :: comm

    Character(Len=200)                       :: record
    Character(Len=256)                       :: message
    Character(Len=4)                         :: keyword
    Character(Len=40)                        :: word
    Character(Len=8)                         :: atom1, atom2
    Integer                                  :: cd, cds, ce, ces, cp, fail(1:2), i, ipot, j, &
                                                jtpatm, k0, katom1, katom2, keymet, ktype, ngrid, &
                                                numpot, ntable
    Integer, Allocatable, Dimension(:)       :: cdens, cdnss, cembds, cembed, cpair
    Logical                                  :: safe
    Real(Kind=wp)                            :: finish, start
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer

    fail = 0
    If (met%tab == 1) Then ! EAM
      Allocate (cpair(1:(met%n_potentials * (met%n_potentials + 1)) / 2), &
                cdens(1:met%n_potentials), &
                cembed(1:met%n_potentials), Stat=fail(1))
    Else If (met%tab == 2) Then ! EEAM
      Allocate (cpair(1:(met%n_potentials * (met%n_potentials + 1)) / 2), &
                cdens(1:met%n_potentials**2), &
                cembed(1:met%n_potentials), Stat=fail(1))
    Else If (met%tab == 3) Then ! 2BEAM
      Allocate (cpair(1:(met%n_potentials * (met%n_potentials + 1)) / 2), &
                cdens(1:met%n_potentials), &
                cdnss(1:met%n_potentials * (met%n_potentials + 1) / 2), &
                cembed(1:met%n_potentials), &
                cembds(1:met%n_potentials), Stat=fail(1))
    Else If (met%tab == 4) Then ! 2BEEAM
      Allocate (cpair(1:(met%n_potentials * (met%n_potentials + 1)) / 2), &
                cdens(1:met%n_potentials**2), &
                cdnss(1:met%n_potentials**2), &
                cembed(1:met%n_potentials), &
                cembds(1:met%n_potentials), Stat=fail(1))
    End If
    Allocate (buffer(1:met%maxgrid), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'metal_table_read allocation failure'
      Call error(0, message)
    End If
    cpair = 0; cp = 0
    cdens = 0; cd = 0
    cembed = 0; ce = 0
    If (met%tab == 3 .or. met%tab == 4) Then
      cdnss = 0; cds = 0
      cembds = 0; ces = 0
    End If

    If (comm%idnode == 0) Then
      Open (Newunit=files(FILE_TABEAM)%unit_no, File=files(FILE_TABEAM)%filename)
      ntable = files(FILE_TABEAM)%unit_no
    End If

    ! skip header record

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Go To 100

    ! read number of potential functions in file

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Go To 100

    Call get_word(record, word)
    numpot = Nint(word_2_real(word))

    Do ipot = 1, numpot

      ! read data type, atom labels, number of points, start and end

      Call get_line(safe, ntable, record, comm)
      If (.not. safe) Go To 100

      ! identify data type

      Call get_word(record, keyword)
      Call lower_case(keyword)
      If (keyword == 'pair') Then
        ktype = 1
      Else If (keyword == 'dens' .or. keyword == 'dden') Then
        ktype = 2
      Else If (keyword == 'embe' .or. keyword == 'demb') Then
        ktype = 3
      Else If (keyword == 'sden') Then
        ktype = 4
      Else If (keyword == 'semb') Then
        ktype = 5
      Else
        Call error(151)
      End If

      ! identify atom types

      Call get_word(record, atom1)
      If (ktype == 1 .or. & ! pair
          (ktype == 2 .and. (met%tab == 2 .or. met%tab == 4)) .or. & ! den for EEAM and dden for 2BEEAM
          (ktype == 4 .and. (met%tab == 3 .or. met%tab == 4))) Then ! sden for 2B(EAM and EEAM)
        Call get_word(record, atom2)
      Else
        atom2 = atom1
      End If

      ! data specifiers

      Call get_word(record, word)
      ngrid = Nint(word_2_real(word))
      Call get_word(record, word)
      start = word_2_real(word)
      Call get_word(record, word)
      finish = word_2_real(word)

      ! check atom identities

      katom1 = 0
      katom2 = 0

      Do jtpatm = 1, sites%ntype_atom
        If (atom1 == sites%unique_atom(jtpatm)) katom1 = jtpatm
        If (atom2 == sites%unique_atom(jtpatm)) katom2 = jtpatm
      End Do

      If (katom1 == 0 .or. katom2 == 0) Then
        If (l_top) &
          Write (message, '(a)') '****', atom1, '***', atom2, '**** entry in TABEAM'
        Call error(81, message, .true.)
      End If

      ! store working parameters

      buffer(1) = Real(ngrid + 4, wp) ! as if there are 4 extra elements after finish
      buffer(4) = (finish - start) / Real(ngrid - 1, wp)
      buffer(2) = start - 5.0_wp * buffer(4)
      buffer(3) = finish

      If (l_top) Then
        Write (message, "(1x,i10,4x,2a8,3x,2a4,2x,i6,1p,3e15.6)") &
          ipot, atom1, atom2, 'EAM-', keyword, ngrid, start, finish, buffer(4)
        Call info(message, .true.)
      End If

      ! check array dimensions

      If (ngrid + 4 > met%maxgrid) Then
        Call warning(270, Real(ngrid + 4, wp), Real(met%maxgrid, wp), 0.0_wp)
        Call error(48)
      End If

      keymet = (Max(katom1, katom2) * (Max(katom1, katom2) - 1)) / 2 + Min(katom1, katom2)
      k0 = met%list(keymet)

      ! check for undefined potential

      If (k0 == 0) Call error(508)

      ! read in potential arrays

      Do i = 1, (ngrid + 3) / 4
        j = Min(4, ngrid - (i - 1) * 4)
        If (comm%idnode == 0) Then
          Read (Unit=ntable, Fmt=*, End=100) buffer(4 * i + 1:4 * i + j)
        Else
          buffer(4 * i + 1:4 * i + j) = 0.0_wp
        End If
      End Do

      Call gsum(comm, buffer(5:ngrid + 4))

      ! copy data to internal arrays

      If (ktype == 1) Then

        ! pair potential terms

        ! Set indices

        !        k0=met%list(keymet)

        cp = cp + 1
        If (Any(cpair(1:cp - 1) == k0)) Then
          Call error(509)
        Else
          cpair(cp) = k0
        End If

        met%vmet(1, k0, 1) = buffer(1)
        met%vmet(2, k0, 1) = buffer(2)
        met%vmet(3, k0, 1) = buffer(3)
        met%vmet(4, k0, 1) = buffer(4)

        Do i = 5, met%maxgrid
          If (i - 4 > ngrid) Then
            met%vmet(i, k0, 1) = 0.0_wp
          Else
            buffer(i) = buffer(i) * engunit
            met%vmet(i, k0, 1) = buffer(i)
          End If
        End Do

        ! calculate derivative of pair potential function

        Call metal_table_derivatives(k0, buffer, Size(met%vmet, 2), met%vmet, met)

        ! adapt derivatives for use in interpolation

        Do i = 5, ngrid + 4
          met%vmet(i, k0, 2) = -(Real(i, wp) * buffer(4) + buffer(2)) * met%vmet(i, k0, 2)
        End Do

      Else If (ktype == 2) Then

        ! density function terms
        ! s-density density function terms for EAM & EEAM
        ! d-density density function terms for 2B(EAM & EEAM)

        ! Set indices

        If (met%tab == 1 .or. met%tab == 3) Then ! EAM
          k0 = katom1
        Else If (met%tab == 2 .or. met%tab == 4) Then ! EEAM
          k0 = (katom1 - 1) * sites%ntype_atom + katom2
        End If

        cd = cd + 1
        If (Any(cdens(1:cd - 1) == k0)) Then
          Call error(510)
        Else
          cdens(cd) = k0
        End If

        met%dmet(1, k0, 1) = buffer(1)
        met%dmet(2, k0, 1) = buffer(2)
        met%dmet(3, k0, 1) = buffer(3)
        met%dmet(4, k0, 1) = buffer(4)

        Do i = 5, met%maxgrid
          If (i - 4 > ngrid) Then
            met%dmet(i, k0, 1) = 0.0_wp
          Else
            met%dmet(i, k0, 1) = buffer(i)
          End If
        End Do

        ! calculate derivative of density function

        Call metal_table_derivatives(k0, buffer, Size(met%dmet, 2), met%dmet, met)

        ! adapt derivatives for use in interpolation

        met%dmet(1, k0, 2) = 0.0_wp
        met%dmet(2, k0, 2) = 0.0_wp
        met%dmet(3, k0, 2) = 0.0_wp
        met%dmet(4, k0, 2) = 0.0_wp

        Do i = 5, ngrid + 4
          met%dmet(i, k0, 2) = -(Real(i, wp) * buffer(4) + buffer(2)) * met%dmet(i, k0, 2)
        End Do

      Else If (ktype == 3) Then

        ! embedding function terms
        ! s-density embedding function terms for EAM & EEAM
        ! d-density embedding function terms for 2B(EAM & EEAM)

        ! Set indices

        k0 = katom1

        ce = ce + 1
        If (Any(cembed(1:ce - 1) == k0)) Then
          Call error(511)
        Else
          cembed(ce) = k0
        End If

        met%fmet(1, k0, 1) = buffer(1)
        met%fmet(2, k0, 1) = buffer(2)
        met%fmet(3, k0, 1) = buffer(3)
        met%fmet(4, k0, 1) = buffer(4)

        Do i = 5, met%maxgrid
          If (i - 4 > ngrid) Then
            met%fmet(i, k0, 1) = 0.0_wp
          Else
            buffer(i) = buffer(i) * engunit
            met%fmet(i, k0, 1) = buffer(i)
          End If
        End Do

        ! calculate derivative of embedding function

        Call metal_table_derivatives(k0, buffer, Size(met%fmet, 2), met%fmet, met)

      Else If (ktype == 4) Then

        ! s-density function terms

        ! The 2BM formalism for alloys allows for a mixed s-band density: rho_{atom1,atom2} /= 0
        ! (and in general for the EEAM it may be non-symmetric: rho_{atom1,atom2} may be /= rho_{atom2,atom2})
        ! Some 2BM models rho_{atom1,atom1}=rho_{atom2,atom2}==0 with rho_{atom1,atom2} /= 0
        ! whereas others choose not to have mixed s-band densities.

        ! Set indices

        If (met%tab == 3) Then ! 2BMEAM
          !           k0=met%list(keymet)
        Else If (met%tab == 4) Then ! 2BMEEAM
          k0 = (katom1 - 1) * sites%ntype_atom + katom2
        End If

        cds = cds + 1
        If (Any(cdnss(1:cds - 1) == k0)) Then
          Call error(510)
        Else
          cdnss(cds) = k0
        End If

        met%dmes(1, k0, 1) = buffer(1)
        met%dmes(2, k0, 1) = buffer(2)
        met%dmes(3, k0, 1) = buffer(3)
        met%dmes(4, k0, 1) = buffer(4)

        If (Nint(buffer(1)) > 5) Then

          Do i = 5, met%maxgrid
            If (i - 4 > ngrid) Then
              met%dmes(i, k0, 1) = 0.0_wp
            Else
              met%dmes(i, k0, 1) = buffer(i)
            End If
          End Do
          ! calculate derivative of density function

          Call metal_table_derivatives(k0, buffer, Size(met%dmes, 2), met%dmes, met)

          ! adapt derivatives for use in interpolation

          met%dmes(1, k0, 2) = 0.0_wp
          met%dmes(2, k0, 2) = 0.0_wp
          met%dmes(3, k0, 2) = 0.0_wp
          met%dmes(4, k0, 2) = 0.0_wp

          Do i = 5, ngrid + 4
            met%dmes(i, k0, 2) = -(Real(i, wp) * buffer(4) + buffer(2)) * met%dmes(i, k0, 2)
          End Do

        End If

      Else If (ktype == 5) Then

        ! s-embedding function terms

        ! Set index

        k0 = katom1

        ces = ces + 1
        If (Any(cembds(1:ces - 1) == k0)) Then
          Call error(511)
        Else
          cembds(ces) = k0
        End If

        met%fmes(1, k0, 1) = buffer(1)
        met%fmes(2, k0, 1) = buffer(2)
        met%fmes(3, k0, 1) = buffer(3)
        met%fmes(4, k0, 1) = buffer(4)

        Do i = 5, met%maxgrid
          If (i - 4 > ngrid) Then
            met%fmes(i, k0, 1) = 0.0_wp
          Else
            buffer(i) = buffer(i) * engunit
            met%fmes(i, k0, 1) = buffer(i)
          End If
        End Do

        ! calculate derivative of embedding function

        Call metal_table_derivatives(k0, buffer, Size(met%fmes, 2), met%fmes, met)

      End If

    End Do

    If (comm%idnode == 0) Call files(FILE_TABEAM)%close ()
    If (l_top) Then
      Write (message, '(a)') 'potential tables read from TABEAM file'
      Call info(message, .true.)
    End If

    If (met%tab == 1 .or. met%tab == 2) Then ! EAM & EEAM
      Deallocate (cpair, cdens, cembed, Stat=fail(1))
    Else If (met%tab == 3 .or. met%tab == 4) Then ! 2B(EAM & EEAM)
      Deallocate (cpair, cdens, cdnss, cembed, cembds, Stat=fail(1))
    End If
    Deallocate (buffer, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a,i0)') 'metal_table_read deallocation failure'
      Call error(0, message)
    End If

    Return

    ! end of file error exit

    100 Continue

    If (comm%idnode == 0) Call files(FILE_TABEAM)%close ()
    Call error(24)

  End Subroutine metal_table_read

  Subroutine metal_generate(ntype_atom, met)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for generating potential energy and force arrays
    ! for metal potentials
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith june 2006
    ! amended   - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer(Kind=wi), Intent(In   ) :: ntype_atom
    Type(metal_type), Intent(InOut) :: met

    Integer       :: i, imet, katom1, katom2, keypot, kmet, mmm, nnn, pmet, qmet
    Real(Kind=wp) :: aaa, bbb, bet, cc0, cc1, cc2, cc3, cc4, ccc, cut1, cut2, ddd, dlrpot, eps, &
                     mmmr, nnnr, ppp, qqq, rr0, rrr, sig

    ! define grid resolution for potential arrays

    dlrpot = met%rcut / Real(met%maxgrid - 1, wp)

    ! construct arrays for metal potentials

    kmet = 0
    Do katom1 = 1, ntype_atom
      Do katom2 = 1, katom1
        kmet = kmet + 1

        ! calculate potentials for defined interactions

        imet = met%list(kmet)
        keypot = met%ltp(imet)
        If (keypot > 0) Then

          ! store array specification parameters

          met%vmet(1, imet, 1) = Real(met%maxgrid, wp)
          met%vmet(2, imet, 1) = 0.0_wp ! l_int(min) >= 1
          met%vmet(3, imet, 1) = met%rcut ! met%rcut=neigh%cutoff
          met%vmet(4, imet, 1) = dlrpot

          Do i = 1, 4
            met%vmet(i, imet, 2) = met%vmet(i, imet, 1)
            met%dmet(i, imet, 1) = met%vmet(i, imet, 1)
            met%dmet(i, imet, 2) = 0.0_wp
          End Do

          If (keypot == 1) Then

            ! finnis-sinclair potentials

            cc0 = met%prm(1, imet)
            cc1 = met%prm(2, imet)
            cc2 = met%prm(3, imet)
            ccc = met%prm(4, imet)
            ddd = met%prm(6, imet)
            bet = met%prm(7, imet)
            cut1 = ccc + 4.0_wp * dlrpot
            cut2 = ddd + 4.0_wp * dlrpot

            met%vmet(3, imet, 1:2) = cut1
            met%dmet(3, imet, 1) = cut2

            Do i = 5, met%maxgrid
              rrr = Real(i, wp) * dlrpot

              If (rrr <= cut1) Then
                met%vmet(i, imet, 1) = (cc0 + cc1 * rrr + cc2 * rrr**2) * (rrr - ccc)**2
                met%vmet(i, imet, 2) = -rrr * (2.0_wp * (cc0 + cc1 * rrr + cc2 * rrr**2) * &
                                               (rrr - ccc) + (cc1 + 2.0_wp * cc2 * rrr) * (rrr - ccc)**2)
              End If

              If (rrr <= cut2) Then
                met%dmet(i, imet, 1) = (rrr - ddd)**2 + bet * (rrr - ddd)**3 / ddd
                met%dmet(i, imet, 2) = -rrr * (2.0_wp * (rrr - ddd) + 3.0_wp * bet * (rrr - ddd)**2 / ddd)
              End If
            End Do

            If (katom1 == katom2) Then
              met%dmet(1, imet, 2) = met%prm(5, imet)**2
              met%dmet(2, imet, 2) = met%dmet(1, imet, 2)
            Else
              pmet = met%list((katom1 * (katom1 + 1)) / 2)
              qmet = met%list((katom2 * (katom2 + 1)) / 2)
              met%dmet(1, imet, 2) = met%prm(5, pmet)**2
              met%dmet(2, imet, 2) = met%prm(5, qmet)**2
            End If

          Else If (keypot == 2) Then

            ! extended finnis-sinclair potentials

            cc0 = met%prm(1, imet)
            cc1 = met%prm(2, imet)
            cc2 = met%prm(3, imet)
            cc3 = met%prm(4, imet)
            cc4 = met%prm(5, imet)
            ccc = met%prm(6, imet)
            ddd = met%prm(8, imet)
            bbb = met%prm(9, imet)
            cut1 = ccc + 4.0_wp * dlrpot
            cut2 = ddd + 4.0_wp * dlrpot

            met%vmet(3, imet, 1:2) = cut1
            met%dmet(3, imet, 1) = cut2

            Do i = 5, met%maxgrid
              rrr = Real(i, wp) * dlrpot

              If (rrr <= cut1) Then
                met%vmet(i, imet, 1) = (cc0 + cc1 * rrr + cc2 * rrr**2 + cc3 * rrr**3 + cc4 * rrr**4) * (rrr - ccc)**2
            met%vmet(i, imet, 2) = -rrr * (2.0_wp * (cc0 + cc1 * rrr + cc2 * rrr**2 + cc3 * rrr**3 + cc4 * rrr**4) * (rrr - ccc) + &
                                        (cc1 + 2.0_wp * cc2 * rrr + 3.0_wp * cc3 * rrr**2 + 4.0_wp * cc4 * rrr**3) * (rrr - ccc)**2)
              End If

              If (rrr <= cut2) Then
                met%dmet(i, imet, 1) = (rrr - ddd)**2 + bbb**2 * (rrr - ddd)**4
                met%dmet(i, imet, 2) = -rrr * (2.0_wp * (rrr - ddd) + 4.0_wp * bbb**2 * (rrr - ddd)**3)
              End If
            End Do

            If (katom1 == katom2) Then
              met%dmet(1, imet, 2) = met%prm(7, imet)**2
              met%dmet(2, imet, 2) = met%dmet(1, imet, 2)
            Else
              pmet = met%list((katom1 * (katom1 + 1)) / 2)
              qmet = met%list((katom2 * (katom2 + 1)) / 2)
              met%dmet(1, imet, 2) = met%prm(7, pmet)**2
              met%dmet(2, imet, 2) = met%prm(7, qmet)**2
            End If

          Else If (keypot == 3) Then

            ! sutton-chen potentials

            eps = met%prm(1, imet)
            sig = met%prm(2, imet)
            nnn = Nint(met%prm(3, imet)); nnnr = Real(nnn, wp)
            mmm = Nint(met%prm(4, imet)); mmmr = Real(mmm, wp)

            Do i = 5, met%maxgrid
              rrr = Real(i, wp) * dlrpot
              met%vmet(i, imet, 1) = eps * (sig / rrr)**nnn
              met%vmet(i, imet, 2) = nnnr * eps * (sig / rrr)**nnn
              met%dmet(i, imet, 1) = (sig / rrr)**mmm
              met%dmet(i, imet, 2) = mmmr * (sig / rrr)**mmm
            End Do

            If (katom1 == katom2) Then
              met%dmet(1, imet, 2) = (met%prm(1, imet) * met%prm(5, imet))**2
              met%dmet(2, imet, 2) = met%dmet(1, imet, 2)
            Else
              pmet = met%list((katom1 * (katom1 + 1)) / 2)
              qmet = met%list((katom2 * (katom2 + 1)) / 2)
              met%dmet(1, imet, 2) = (met%prm(1, pmet) * met%prm(5, pmet))**2
              met%dmet(2, imet, 2) = (met%prm(1, qmet) * met%prm(5, qmet))**2
            End If

          Else If (keypot == 4) Then

            ! gupta potentials

            aaa = met%prm(1, imet)
            rr0 = met%prm(2, imet)
            ppp = met%prm(3, imet)
            qqq = met%prm(5, imet)

            Do i = 5, met%maxgrid
              rrr = Real(i, wp) * dlrpot

              cut1 = (rrr - rr0) / rr0
              cut2 = cut1 + 1.0_wp

              met%vmet(i, imet, 1) = 2.0_wp * aaa * Exp(-ppp * cut1)
              met%vmet(i, imet, 2) = met%vmet(i, imet, 1) * ppp * cut2
              met%dmet(i, imet, 1) = Exp(-2.0_wp * qqq * cut1)
              met%dmet(i, imet, 2) = 2.0_wp * met%dmet(i, imet, 1) * qqq * cut2
            End Do

            met%dmet(1, imet, 2) = met%prm(4, imet)**2
            met%dmet(2, imet, 2) = met%dmet(1, imet, 2)

          Else If (keypot == 5) Then

            ! many-body perturbation component only potentials

            eps = met%prm(1, imet)
            sig = met%prm(2, imet)
            mmm = Nint(met%prm(3, imet)); mmmr = Real(mmm, wp)

            Do i = 5, met%maxgrid
              rrr = Real(i, wp) * dlrpot
              !                 met%vmet(i,imet,1)=0.0_wp
              !                 met%vmet(i,imet,2)=0.0_wp
              nnnr = sig / rrr**mmm
              met%dmet(i, imet, 1) = nnnr * met%merf(i)
              met%dmet(i, imet, 2) = mmmr * met%dmet(i, imet, 1) - rrr * nnnr * met%mfer(i)
            End Do

            If (katom1 == katom2) Then
              met%dmet(1, imet, 2) = met%prm(1, imet)**2
              met%dmet(2, imet, 2) = met%dmet(1, imet, 2)
            Else
              pmet = met%list((katom1 * (katom1 + 1)) / 2)
              qmet = met%list((katom2 * (katom2 + 1)) / 2)
              met%dmet(1, imet, 2) = met%prm(1, pmet)**2
              met%dmet(2, imet, 2) = met%prm(1, qmet)**2
            End If

          Else

            Call error(151)

          End If

        End If
      End Do
    End Do
  End Subroutine metal_generate

  Subroutine metal_generate_erf(met)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for generating erf and fer arrays for
    ! many-body perturbation component only potentials
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(metal_type), Intent(InOut) :: met

    Integer       :: imet
    Real(Kind=wp) :: alpha, beta

    ! Determine alpha and beta for the erf bit of all MBPC potentials

    If (Any(met%ltp(1:met%n_potentials) == 5)) Then ! all are == 5 == MBPC
      alpha = 0.0_wp
      beta = 0.0_wp
      Do imet = 1, met%n_potentials
        alpha = Max(alpha, Abs(met%prm(6, imet)))
        beta = Max(beta, Abs(met%prm(7, imet)))
      End Do

      ! If unset then set to defaults

      If (alpha <= zero_plus) alpha = 20.0_wp
      If (beta <= zero_plus) beta = Min(1.5_wp, 0.2_wp * met%rcut)

      ! Allocate arrays: met%merf,met%mfer

      Call allocate_metal_erf_arrays(met)

      ! Generate error function and derivative arrays

      Call erfgen_met(alpha, beta, met)

      ! Translate met%merf and met%mfer to the functional form 0.5*{1+erf[alpha(r-beta)]}

      met%merf(5:met%maxgrid) = 0.5_wp * (1.0_wp + met%merf(5:met%maxgrid))
      met%mfer(5:met%maxgrid) = 0.5_wp * met%mfer(5:met%maxgrid)
    End If
  End Subroutine metal_generate_erf

  Subroutine metal_ld_collect_eam(iatm, rrt, safe, ntype_atom, met, neigh, config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating local atomic density for
    ! Embedded Atom Model & Extended Embedded Atom Model metal potentials
    !
    ! Note: Designed to be used as part of metal_ld_compute
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith june 1995
    ! amended   - i.t.todorov december 2016
    ! contrib   - r.davidchak (eeam) june 2012
    ! contrib   - b.palmer (2band) may 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                    Intent(In   ) :: iatm
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(metal_type),                           Intent(InOut) :: met
    Integer(Kind=wi),                           Intent(In   ) :: ntype_atom
    Logical,                                    Intent(InOut) :: safe
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: rrt
    Type(configuration_type),                   Intent(InOut) :: config

    Integer       :: ai, aj, jatm, k0, key, ki, kj, l, m
    Real(Kind=wp) :: density, ppp, rdr, rr1, rrr, t1, t2, vk0, vk1, vk2

    ! global type of itam

    ai = config%ltype(iatm)

    ! start of primary loop for density

    Do m = 1, neigh%list(0, iatm)

      ! atomic and potential function indices

      jatm = neigh%list(m, iatm)
      aj = config%ltype(jatm)

      If (met%tab == 1 .or. met%tab == 3) Then ! EAM & 2BEAM
        ki = ai
        kj = aj
      Else If (met%tab == 2 .or. met%tab == 4) Then ! EEAM & 2BEEAM
        ki = (aj - 1) * ntype_atom + ai ! aj-ai
        kj = (ai - 1) * ntype_atom + aj ! ai-aj
      End If

      ! interatomic distance

      rrr = rrt(m)

      ! Now start traditional s-band (EAM & EEAM) or d-band for 2B(EAM & EEAM)

      ! first metal atom density and validity and truncation of potential

      If (Abs(met%dmet(1, kj, 1)) > zero_plus .and. Nint(met%dmet(1, kj, 1)) > 5) Then
        If (rrr <= met%dmet(3, kj, 1)) Then

          ! interpolation parameters

          rdr = 1.0_wp / met%dmet(4, kj, 1)
          rr1 = rrr - met%dmet(2, kj, 1)
          l = Min(Nint(rr1 * rdr), Nint(met%dmet(1, kj, 1)) - 1)
          If (l < 5) Then ! catch unsafe value
            safe = .false.
            Write (*, *) 'aaa', l, iatm, jatm, rrr
            l = 6
          End If
          ppp = rr1 * rdr - Real(l, wp)

          ! calculate density using 3-point interpolation

          vk0 = met%dmet(l - 1, kj, 1)
          vk1 = met%dmet(l, kj, 1)
          vk2 = met%dmet(l + 1, kj, 1)

          t1 = vk1 + ppp * (vk1 - vk0)
          t2 = vk1 + ppp * (vk2 - vk1)

          If (ppp < 0.0_wp) Then ! density is a positive function!
            density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
            If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
          Else If (l == 5) Then
            density = t2
          Else
            density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
            If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
          End If

          met%rho(iatm) = met%rho(iatm) + density
          If (ki == kj .and. jatm <= config%natms) met%rho(jatm) = met%rho(jatm) + density

        End If
      End If

      ! second metal atom density and validity and truncation of potential

      If (Abs(met%dmet(1, ki, 1)) > zero_plus .and. Nint(met%dmet(1, ki, 1)) > 5) Then
        If (ki /= kj .and. jatm <= config%natms) Then
          If (rrr <= met%dmet(3, ki, 1)) Then

            ! interpolation parameters

            rdr = 1.0_wp / met%dmet(4, ki, 1)
            rr1 = rrr - met%dmet(2, ki, 1)
            l = Min(Nint(rr1 * rdr), Nint(met%dmet(1, ki, 1)) - 1)
            If (l < 5) Then ! catch unsafe value
              safe = .false.
              Write (*, *) 'bbb', l, iatm, jatm, rrr
              l = 6
            End If
            ppp = rr1 * rdr - Real(l, wp)

            ! calculate density using 3-point interpolation

            vk0 = met%dmet(l - 1, ki, 1)
            vk1 = met%dmet(l, ki, 1)
            vk2 = met%dmet(l + 1, ki, 1)

            t1 = vk1 + ppp * (vk1 - vk0)
            t2 = vk1 + ppp * (vk2 - vk1)

            If (ppp < 0.0_wp) Then ! density is a positive function!
              density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
              If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
            Else If (l == 5) Then
              density = t2
            Else
              density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
              If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
            End If

            met%rho(jatm) = met%rho(jatm) + density

          End If
        End If
      End If

      ! Now if we have the 2B(EAM & EEAM) then do the s-band (met%dmes and met%rhs are defined)

      If (met%l_2b) Then
        If (met%tab == 3) Then ! 2BEAM

          ! 2BEAM has symmetric s-densities with respect to atom type
          ! e.g. rho_(atom1,atom1), rho_(atom1,atom2) = rho_(atom2,atom1), rho_(atom2,atom2)

          key = (Max(ai, aj) * (Max(ai, aj) - 1)) / 2 + Min(ai, aj)
          k0 = met%list(key)

          ! first metal atom density and validity and truncation of potential

          If (Abs(met%dmes(1, k0, 1)) > zero_plus .and. Nint(met%dmes(1, k0, 1)) > 5) Then
            If (rrr <= met%dmes(3, k0, 1)) Then

              ! interpolation parameters

              rdr = 1.0_wp / met%dmes(4, k0, 1)
              rr1 = rrr - met%dmes(2, k0, 1)
              l = Min(Nint(rr1 * rdr), Nint(met%dmes(1, k0, 1)) - 1)
              If (l < 5) Then ! catch unsafe value
                safe = .false.
                Write (*, *) 'ccc', l, iatm, jatm, rrr
                l = 6
              End If
              ppp = rr1 * rdr - Real(l, wp)

              ! calculate density using 3-point interpolation

              vk0 = met%dmes(l - 1, k0, 1)
              vk1 = met%dmes(l, k0, 1)
              vk2 = met%dmes(l + 1, k0, 1)

              t1 = vk1 + ppp * (vk1 - vk0)
              t2 = vk1 + ppp * (vk2 - vk1)

              If (ppp < 0.0_wp) Then ! density is a positive function!
                density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
              Else If (l == 5) Then
                density = t2
              Else
                density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
              End If

              met%rhs(iatm) = met%rhs(iatm) + density
              If (jatm <= config%natms) met%rhs(jatm) = met%rhs(jatm) + density

            End If
          End If

        Else If (met%tab == 4) Then ! 2BEEAM

          ! first metal atom density and validity and truncation of potential

          If (Abs(met%dmes(1, kj, 1)) > zero_plus .and. Nint(met%dmes(1, kj, 1)) > 5) Then
            If (rrr <= met%dmes(3, kj, 1)) Then

              ! interpolation parameters

              rdr = 1.0_wp / met%dmes(4, kj, 1)
              rr1 = rrr - met%dmes(2, kj, 1)
              l = Min(Nint(rr1 * rdr), Nint(met%dmes(1, kj, 1)) - 1)
              If (l < 5) Then ! catch unsafe value
                safe = .false.
                Write (*, *) 'ddd', l, iatm, jatm, rrr
                l = 6
              End If
              ppp = rr1 * rdr - Real(l, wp)

              ! calculate density using 3-point interpolation

              vk0 = met%dmes(l - 1, kj, 1)
              vk1 = met%dmes(l, kj, 1)
              vk2 = met%dmes(l + 1, kj, 1)

              t1 = vk1 + ppp * (vk1 - vk0)
              t2 = vk1 + ppp * (vk2 - vk1)

              If (ppp < 0.0_wp) Then ! density is a positive function!
                density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
              Else If (l == 5) Then
                density = t2
              Else
                density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
              End If

              met%rhs(iatm) = met%rhs(iatm) + density
              If (ki == kj .and. jatm <= config%natms) met%rhs(jatm) = met%rhs(jatm) + density

            End If
          End If

          ! second metal atom density and validity and truncation of potential

          If (Abs(met%dmes(1, ki, 1)) > zero_plus .and. Nint(met%dmes(1, ki, 1)) > 5) Then
            If (ki /= kj .and. jatm <= config%natms) Then
              If (rrr <= met%dmes(3, ki, 1)) Then

                ! interpolation parameters

                rdr = 1.0_wp / met%dmes(4, ki, 1)
                rr1 = rrr - met%dmes(2, ki, 1)
                l = Min(Nint(rr1 * rdr), Nint(met%dmes(1, ki, 1)) - 1)
                If (l < 5) Then ! catch unsafe value
                  safe = .false.
                  Write (*, *) 'eee', l, iatm, jatm, rrr
                  l = 6
                End If
                ppp = rr1 * rdr - Real(l, wp)

                ! calculate density using 3-point interpolation

                vk0 = met%dmes(l - 1, ki, 1)
                vk1 = met%dmes(l, ki, 1)
                vk2 = met%dmes(l + 1, ki, 1)

                t1 = vk1 + ppp * (vk1 - vk0)
                t2 = vk1 + ppp * (vk2 - vk1)

                If (ppp < 0.0_wp) Then ! density is a positive function!
                  density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
                  If (density < 0.0_wp) density = t1 ! for non-smooth descend to zero, or ascend from zero
                Else If (l == 5) Then
                  density = t2
                Else
                  density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
                  If (density < 0.0_wp) density = t2 ! for non-smooth descend to zero, or ascend from zero
                End If

                met%rhs(jatm) = met%rhs(jatm) + density

              End If
            End If
          End If
        End If
      End If
    End Do
  End Subroutine metal_ld_collect_eam

  Subroutine metal_ld_collect_fst(iatm, rrt, safe, met, neigh, config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating local atomic density for
    ! Finnis-Sinclair Type metal potentials
    !
    ! Note: Designed to be used as part of metal_ld_compute
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith june 1995
    ! amended   - i.t.todorov december 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                    Intent(In   ) :: iatm
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(metal_type),                           Intent(InOut) :: met
    Logical,                                    Intent(InOut) :: safe
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: rrt
    Type(configuration_type),                   Intent(In   ) :: config

    Integer       :: ai, aj, jatm, k0, k1, k2, key, keypot, kmn, kmx, l, m, mmm, nnn
    Real(Kind=wp) :: aaa, bbb, bet, cc0, cc1, cc2, cc3, cc4, ccc, cut1, cut2, ddd, density, eps, &
                     ppp, qqq, rdr, rr0, rr1, rrr, sig, t1, t2, vk0, vk1, vk2

    ! global type of itam

    ai = config%ltype(iatm)

    ! start of primary loop for density

    Do m = 1, neigh%list(0, iatm)

      ! atomic and potential function indices

      jatm = neigh%list(m, iatm)
      aj = config%ltype(jatm)

      If (ai > aj) Then
        key = ai * (ai - 1) / 2 + aj
      Else
        key = aj * (aj - 1) / 2 + ai
      End If

      k0 = met%list(key)

      If (met%l_direct) Then
        k1 = Max(ai, aj)
        k2 = Min(ai, aj)

        kmx = k1 * (k1 + 1) / 2
        kmn = k2 * (k2 + 1) / 2

        k1 = met%list(kmx)
        k2 = met%list(kmn)
      End If

      ! interatomic distance

      rrr = rrt(m)

      ! validity and truncation of analytic potential

      keypot = met%ltp(k0)
      If (keypot > 0 .and. rrr <= met%rcut) Then

        ! Abs(met%dmet(1,k0,1)) > zero_plus, as potentials are analytic

        If (met%l_direct) Then ! direct calculation

          If (keypot == 1) Then

            ! finnis-sinclair potentials

            cc0 = met%prm(1, k0)
            cc1 = met%prm(2, k0)
            cc2 = met%prm(3, k0)
            ccc = met%prm(4, k0)
            ddd = met%prm(6, k0)
            bet = met%prm(7, k0)
            cut1 = ccc
            cut2 = ddd

            density = 0.0_wp
            If (rrr <= cut2) density = (rrr - ddd)**2 + bet * (rrr - ddd)**3 / ddd

            If (ai == aj) Then
              t1 = met%prm(5, k0)**2
              t2 = t1
            Else
              t1 = met%prm(5, k1)**2
              t2 = met%prm(5, k2)**2
            End If

          Else If (keypot == 2) Then

            ! extended finnis-sinclair potentials

            cc0 = met%prm(1, k0)
            cc1 = met%prm(2, k0)
            cc2 = met%prm(3, k0)
            cc3 = met%prm(4, k0)
            cc4 = met%prm(5, k0)
            ccc = met%prm(6, k0)
            ddd = met%prm(8, k0)
            bbb = met%prm(9, k0)
            cut1 = ccc
            cut2 = ddd

            density = 0.0_wp
            If (rrr <= cut2) density = (rrr - ddd)**2 + bbb**2 * (rrr - ddd)**4

            If (ai == aj) Then
              t1 = met%prm(7, k0)**2
              t2 = t1
            Else
              t1 = met%prm(7, k1)**2
              t2 = met%prm(7, k2)**2
            End If

          Else If (keypot == 3) Then

            ! sutton-chen potentials

            eps = met%prm(1, k0)
            sig = met%prm(2, k0)
            nnn = Nint(met%prm(3, k0))
            mmm = Nint(met%prm(4, k0))

            density = (sig / rrr)**mmm

            If (ai == aj) Then
              t1 = (met%prm(1, k0) * met%prm(5, k0))**2
              t2 = t1
            Else
              t1 = (met%prm(1, k1) * met%prm(5, k1))**2
              t2 = (met%prm(1, k2) * met%prm(5, k2))**2
            End If

          Else If (keypot == 4) Then

            ! gupta potentials

            aaa = met%prm(1, k0)
            rr0 = met%prm(2, k0)
            ppp = met%prm(3, k0)
            qqq = met%prm(5, k0)

            density = Exp(-2.0_wp * qqq * (rrr - rr0) / rr0)

            t1 = met%prm(4, k0)**2
            t2 = t1

          Else If (keypot == 5) Then

            ! many-body perturbation component only potentials

            eps = met%prm(1, k0)
            sig = met%prm(2, k0)
            mmm = Nint(met%prm(3, k0))

            ! interpolation parameters

            rdr = 1.0_wp / met%merf(4)
            rr1 = rrr - met%merf(2)
            l = Min(Nint(rr1 * rdr), Nint(met%merf(1)) - 1)
            If (l < 5) Then ! catch unsafe value
              safe = .false.
              Write (*, *) 'aaa', l, iatm, jatm, rrr
              l = 6
            End If
            ppp = rr1 * rdr - Real(l, wp)

            ! calculate density using 3-point interpolation

            vk0 = met%merf(l - 1)
            vk1 = met%merf(l)
            vk2 = met%merf(l + 1)

            t1 = vk1 + ppp * (vk1 - vk0)
            t2 = vk1 + ppp * (vk2 - vk1)

            If (ppp < 0.0_wp) Then
              density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
            Else If (l == 5) Then
              density = t2
            Else
              density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
            End If
            density = density * sig / rrr**mmm

            If (ai == aj) Then
              t1 = met%prm(1, k0)**2
              t2 = t1
            Else
              t1 = met%prm(1, k1)**2
              t2 = met%prm(1, k2)**2
            End If

          End If

          If (ai > aj) Then
            met%rho(iatm) = met%rho(iatm) + density * t1
            If (jatm <= config%natms) met%rho(jatm) = met%rho(jatm) + density * t2
          Else
            met%rho(iatm) = met%rho(iatm) + density * t2
            If (jatm <= config%natms) met%rho(jatm) = met%rho(jatm) + density * t1
          End If

        Else ! tabulated calculation

          ! truncation of potential

          If (rrr <= met%dmet(3, k0, 1)) Then

            ! interpolation parameters

            rdr = 1.0_wp / met%dmet(4, k0, 1)
            rr1 = rrr - met%dmet(2, k0, 1)
            l = Min(Nint(rr1 * rdr), Nint(met%dmet(1, k0, 1)) - 1)
            If (l < 5) Then ! catch unsafe value
              safe = .false.
              Write (*, *) 'bbb', l, iatm, jatm, rrr
              l = 6
            End If
            ppp = rr1 * rdr - Real(l, wp)

            ! calculate density using 3-point interpolation

            vk0 = met%dmet(l - 1, k0, 1)
            vk1 = met%dmet(l, k0, 1)
            vk2 = met%dmet(l + 1, k0, 1)

            t1 = vk1 + ppp * (vk1 - vk0)
            t2 = vk1 + ppp * (vk2 - vk1)

            If (ppp < 0.0_wp) Then
              density = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)
            Else If (l == 5) Then
              density = t2
            Else
              density = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)
            End If

            If (ai > aj) Then
              met%rho(iatm) = met%rho(iatm) + density * met%dmet(1, k0, 2)
              If (jatm <= config%natms) met%rho(jatm) = met%rho(jatm) + density * met%dmet(2, k0, 2)
            Else
              met%rho(iatm) = met%rho(iatm) + density * met%dmet(2, k0, 2)
              If (jatm <= config%natms) met%rho(jatm) = met%rho(jatm) + density * met%dmet(1, k0, 2)
            End If

          End If

        End If

      End If

    End Do
  End Subroutine metal_ld_collect_fst

  Subroutine metal_table_derivatives(ityp, buffer, v2d, vvv, met)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to calculate numerical derivatives of tabulated
    ! EAM metal potentials
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith march 2006
    ! amended   - i.t.todorov april 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent(In   ) :: ityp, v2d
    Type(metal_type), Intent(InOut) :: met
    Real(Kind=wp),    Intent(InOut) :: vvv(1:met%maxgrid, 1:v2d, 1:2)
    Real(Kind=wp),    Intent(In   ) :: buffer(1:met%maxgrid)

    Integer       :: i, i_end, i_start, v_end
    Real(Kind=wp) :: aa0, aa1, aa2, aa3, aa4, d1y, d2y, d3y, d4y, delmet, f0, f1, f2, f3, f4

    ! interpolation parameters

    vvv(1, ityp, 2) = buffer(1)
    vvv(2, ityp, 2) = buffer(2)
    vvv(3, ityp, 2) = buffer(3)
    vvv(4, ityp, 2) = buffer(4)

    ! construct interpolation table

    delmet = buffer(4)
    v_end = Nint(buffer(1))
    i_start = 5 + 2
    i_end = v_end - 2
    Do i = i_start, i_end
      aa0 = buffer(i)
      If (Abs(aa0) <= zero_plus) Then
        f0 = 0.0_wp
        f1 = 0.0_wp
        f2 = 0.0_wp
        f3 = 0.0_wp
        f4 = 0.0_wp
      Else
        f0 = buffer(i - 2) / aa0
        f1 = buffer(i - 1) / aa0
        f2 = 1.0_wp
        f3 = buffer(i + 1) / aa0
        f4 = buffer(i + 2) / aa0
      End If

      ! calculate numerical differences for 5-point interpolation

      d1y = (f1 - f0)
      d2y = (f2 - f1) - (f1 - f0)
      d3y = (f3 - f0) + 3.0_wp * (f1 - f2)
      d4y = (f4 - f3) + 3.0_wp * (f2 - f3) + 3.0_wp * (f2 - f1) + (f0 - f1)

      ! calculate polynomial coefficients

      aa0 = aa0 / delmet
      aa4 = d4y / 24.0_wp
      aa3 = (d3y + 12.0_wp * aa4) / 6.0_wp
      aa2 = (d2y + 6.0_wp * aa3 - 14.0_wp * aa4) / 2.0_wp
      aa1 = d1y + 3.0_wp * aa2 - 7.0_wp * aa3 + 15.0_wp * aa4

      ! calculate derivatives

      vvv(i, ityp, 2) = aa1 * aa0

      ! derivatives at extremes of range

      If (i == i_start) Then
        vvv(i_start - 2, ityp, 2) = (aa1 - 4.0_wp * aa2 + 12.0_wp * aa3 - 32.0_wp * aa4) * aa0
        vvv(i_start - 1, ityp, 2) = (aa1 - 2.0_wp * aa2 + 3.0_wp * aa3 - 4.0_wp * aa4) * aa0
      Else If (i == i_end) Then
        vvv(i_end + 1, ityp, 2) = (aa1 + 2.0_wp * aa2 + 3.0_wp * aa3 + 4.0_wp * aa4) * aa0
        vvv(i_end + 2, ityp, 2) = (aa1 + 4.0_wp * aa2 + 12.0_wp * aa3 + 32.0_wp * aa4) * aa0
      End If
    End Do

    ! set derivatives to zero beyond end point of function

    Do i = v_end + 3, met%maxgrid
      vvv(i, ityp, 2) = 0.0_wp
    End Do
  End Subroutine metal_table_derivatives

  Subroutine metal_ld_export(mdir, mlast, mxatms, ixyz0, met, domain, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export metal density data in domain boundary
    ! regions for halo formation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent(In   ) :: mdir
    Integer,            Intent(InOut) :: mlast
    Integer,            Intent(In   ) :: mxatms
    Integer,            Intent(InOut) :: ixyz0(:)
    Type(metal_type),   Intent(InOut) :: met
    Type(domains_type), Intent(In   ) :: domain
    Type(comms_type),   Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Integer                                  :: fail, i, iadd, iblock, imove, itmp, ix, iy, iz, j, &
                                                jdnode, jmove, jxyz, kdnode, kx, kxyz, ky, kz, &
                                                limit
    Logical                                  :: lrhs, safe
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer

    ! Number of transported quantities per particle

    If (.not. (met%tab == 3 .or. met%tab == 4)) Then
      lrhs = .false.
      iadd = 2
    Else
      lrhs = .true.
      iadd = 3
    End If

    fail = 0; limit = iadd * domain%mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'metal_ld_export allocation failure'
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
      Call error(47)
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

            ! pack particle density and halo indexing

            If (.not. lrhs) Then
              buffer(imove + 1) = met%rho(i)

              ! Use the corrected halo reduction factor when the particle is halo to both +&- sides

              buffer(imove + 2) = Real(ixyz0(i) - Merge(jxyz, kxyz, j == jxyz), wp)
            Else
              buffer(imove + 1) = met%rho(i)
              buffer(imove + 2) = met%rhs(i)

              ! Use the corrected halo reduction factor when the particle is halo to both +&- sides

              buffer(imove + 3) = Real(ixyz0(i) - Merge(jxyz, kxyz, j == jxyz), wp)
            End If

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
      Call error(38)
    End If

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm, jmove, kdnode, MetLdExp_tag)
      Call gsend(comm, imove, jdnode, MetLdExp_tag)
      Call gwait(comm)
    Else
      jmove = imove
    End If

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe = ((mlast + jmove / iadd) <= mxatms)
    Call gcheck(comm, safe)
    If (.not. safe) Then
      itmp = mlast + jmove / iadd
      Call gmax(comm, itmp)
      Call warning(160, Real(itmp, wp), Real(mxatms, wp), 0.0_wp)
      Call error(39)
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm, buffer(iblock + 1:iblock + jmove), kdnode, MetLdExp_tag)
      End If
      If (imove > 0) Then
        Call gsend(comm, buffer(1:imove), jdnode, MetLdExp_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data

    j = Merge(iblock, 0, comm%mxnode > 1)
    Do i = 1, jmove / iadd
      mlast = mlast + 1

      ! unpack particle density and remaining halo indexing

      If (.not. lrhs) Then
        met%rho(mlast) = buffer(j + 1)
        ixyz0(mlast) = Nint(buffer(j + 2))
      Else
        met%rho(mlast) = buffer(j + 1)
        met%rhs(mlast) = buffer(j + 2)
        ixyz0(mlast) = Nint(buffer(j + 3))
      End If

      j = j + iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'metal_ld_export deallocation failure'
      Call error(0, message)
    End If
  End Subroutine metal_ld_export

  Subroutine metal_ld_set_halo(met, domain, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange exchange of density data between
    ! neighbouring domains/nodes
    !
    ! Note: all depends on the ixyz halo array set in set_halo, this assumes
    !       that (i) met%rcut=neigh%cutoff! as well as (ii) all the error checks in there
    !
    ! copyright - daresbury laboratory
    ! amended   - i.t.todorov february 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(metal_type),         Intent(InOut) :: met
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
      Write (message, '(a)') 'metal_ld_set_halo allocation failure'
      Call error(0, message)
    End If
    ixyz0(1:config%nlast) = config%ixyz(1:config%nlast)

    ! No halo, start with domain only particles

    mlast = config%natms

    ! exchange atom data in -/+ x directions

    Call metal_ld_export(-1, mlast, config%mxatms, ixyz0, met, domain, comm)
    Call metal_ld_export(1, mlast, config%mxatms, ixyz0, met, domain, comm)

    ! exchange atom data in -/+ y directions

    Call metal_ld_export(-2, mlast, config%mxatms, ixyz0, met, domain, comm)
    Call metal_ld_export(2, mlast, config%mxatms, ixyz0, met, domain, comm)

    ! exchange atom data in -/+ z directions

    Call metal_ld_export(-3, mlast, config%mxatms, ixyz0, met, domain, comm)
    Call metal_ld_export(3, mlast, config%mxatms, ixyz0, met, domain, comm)

    ! check atom totals after data transfer

    safe = (mlast == config%nlast)
    Call gcheck(comm, safe)
    If (.not. safe) Call error(96)

    Deallocate (ixyz0, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'metal_ld_set_halo deallocation failure'
      Call error(0, message)
    End If
  End Subroutine metal_ld_set_halo

  Subroutine erfgen_met(alpha, beta, met)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine for generating interpolation tables for a magnified
    ! offset error function and its true derivative - for use with MBPC type
    ! of metal-like potentials
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),    Intent(In   ) :: alpha, beta
    Type(metal_type), Intent(InOut) :: met

    Real(Kind=wp), Parameter :: a1 = 0.254829592_wp, a2 = -0.284496736_wp, a3 = 1.421413741_wp, &
                                a4 = -1.453152027_wp, a5 = 1.061405429_wp, pp = 0.3275911_wp

    Integer       :: i, offset
    Real(Kind=wp) :: drmet, exp1, rrr, rsq, tt

    ! look-up tables for mbpc metal interaction

    drmet = met%rcut / Real(met%maxgrid - 1, wp)

    ! store array specification parameters

    met%merf(1) = Real(met%maxgrid, wp)
    met%merf(2) = 0.0_wp ! l_int(min) >= 1
    met%merf(3) = met%rcut ! met%rcut=neigh%cutoff
    met%merf(4) = drmet

    ! offset

    offset = Nint(beta / drmet) + 1

    Do i = 1, 4
      met%mfer(i) = met%merf(i)
    End Do

    Do i = 5, offset - 1
      rrr = -Real(i - offset, wp) * drmet ! make positive
      rsq = rrr * rrr

      tt = 1.0_wp / (1.0_wp + pp * alpha * rrr)
      exp1 = Exp(-(alpha * rrr)**2)

      met%merf(i) = tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1 / rrr - 1.0_wp
      met%mfer(i) = -2.0_wp * (alpha / sqrpi) * exp1
    End Do

    met%merf(offset) = 0.0_wp
    met%mfer(offset) = 2.0_wp * (alpha / sqrpi)

    Do i = offset + 1, met%maxgrid
      rrr = Real(i - offset, wp) * drmet

      tt = 1.0_wp / (1.0_wp + pp * alpha * rrr)
      exp1 = Exp(-(alpha * rrr)**2)

      met%merf(i) = 1.0_wp - tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1 / rrr
      met%mfer(i) = 2.0_wp * (alpha / sqrpi) * exp1
    End Do
  End Subroutine erfgen_met

  Subroutine cleanup(met)
    Type(metal_type) :: met

    If (Allocated(met%list)) Then
      Deallocate (met%list)
    End If
    If (Allocated(met%ltp)) Then
      Deallocate (met%ltp)
    End If
    If (Allocated(met%prm)) Then
      Deallocate (met%prm)
    End If

    If (Allocated(met%elrc)) Then
      Deallocate (met%elrc)
    End If
    If (Allocated(met%vlrc)) Then
      Deallocate (met%vlrc)
    End If

    If (Allocated(met%vmet)) Then
      Deallocate (met%vmet)
    End If
    If (Allocated(met%dmet)) Then
      Deallocate (met%dmet)
    End If
    If (Allocated(met%dmes)) Then
      Deallocate (met%dmes)
    End If
    If (Allocated(met%fmet)) Then
      Deallocate (met%fmet)
    End If
    If (Allocated(met%fmes)) Then
      Deallocate (met%fmes)
    End If

    If (Allocated(met%rho)) Then
      Deallocate (met%rho)
    End If
    If (Allocated(met%rhs)) Then
      Deallocate (met%rhs)
    End If

    If (Allocated(met%merf)) Then
      Deallocate (met%merf)
    End If
    If (Allocated(met%mfer)) Then
      Deallocate (met%mfer)
    End If
  End Subroutine cleanup
End Module metal
