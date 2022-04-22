Module ttm_track

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for utility functions for initial energy deposition
  ! in the electronic subsystem in the two-temperature model (ttm)
  !
  ! copyright - daresbury laboratory
  ! authors   - s.l.daraszewicz & m.a.seaton july 2012
  ! contrib   - g.khara may 2016
  ! contrib   - m.a.seaton february 2017
  ! contrib   - m.a.seaton january 2020
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: Grid3_tag,&
                             Grid4_tag,&
                             comms_type,&
                             gcheck,&
                             gmax,&
                             gmin,&
                             grid1_tag,&
                             grid2_tag,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: boltz,&
                             zero_plus,&
                             eV_to_kB
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use kinds,           Only: wp,STR_LEN
  Use thermostat,      Only: thermostat_type
  Use ttm,             Only: boundaryCond,&
                             boundaryHalo,&
                             depoinit,&
                             eltemp_max,&
                             eltemp_maxKe,&
                             eltemp_mean,&
                             eltemp_min,&
                             eltemp_minKe,&
                             ttm_system_revive,&
                             ttm_type,&
                             TTM_TDEPO_GAUSS, TTM_TDEPO_EXP, TTM_TDEPO_DELTA, TTM_TDEPO_PULSE,&
                             TTM_CE_CONST, TTM_CE_CONST_DYN, TTM_CE_TANH, TTM_CE_TANH_DYN,&
                             TTM_CE_LINEAR, TTM_CE_LINEAR_DYN, TTM_CE_TABULATED, TTM_CE_TABULATED_DYN,&
                             TTM_EPVAR_NULL, TTM_EPVAR_HOMO, TTM_EPVAR_HETERO,&
                             TTM_KE_INFINITE, TTM_KE_CONST, TTM_KE_DRUDE, TTM_KE_TABULATED,&
                             TTM_SDEPO_NULL
  Use ttm_utils,       Only: Ce,&
                             Ke,&
                             alp,&
                             gep,&
                             ked,&
                             peakProfiler,&
                             peakProfilerElec,&
                             printElecLatticeStatsToFile,&
                             printLatticeStatsToFile,&
                             redistribute_te

#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif
  Implicit None

Contains

  Subroutine depoevolve(time, tstep, redtstep, redtstepmx, ttm, comm)

    ! determine how deposition evolves over time
    Real(Kind=wp),    Intent(In   ) :: time, tstep
    Integer,          Intent(In   ) :: redtstep, redtstepmx
    Type(ttm_type),   Intent(InOut) :: ttm
    Type(comms_type), Intent(InOut) :: comm

    Character(Len=14)       :: number
    Character(Len=STR_LEN)      :: messages(2)
    Integer                 :: i, ijk, j, k
    Integer, Dimension(1:3) :: fail = 0
    Logical                 :: deposit
    Real(Kind=wp)           :: adjeng, adjtime, Ce0a, Cemaxa, currenttime, end_Te, energy_diff, &
                               err_tol = 0.01_wp, increase, invbin, lat_I_max, lat_I_min, &
                               lat_I_sum, lat_U_max, lat_U_min, lat_U_sum, newCe, oldCe, sh_Aa, &
                               start_Te

    ! provide atomic density corrections to heat capacities

    Ce0a = ttm%Ce0
    sh_Aa = ttm%sh_A
    Cemaxa = ttm%Cemax
    if (ttm%ttmdyndens) then
      Ce0a = ttm%Ce0 * ttm%cellrho
      sh_Aa = ttm%sh_A * ttm%cellrho
      Cemaxa = ttm%Cemax * ttm%cellrho
    end if

    ! start deposition, reducing size of timestep for thermal diffusion

    currenttime = time - ttm%depostart + tstep / Real(redtstepmx, Kind=wp) * Real(redtstep, Kind=wp)

    Select Case (ttm%tdepoType)
    Case (TTM_TDEPO_GAUSS)
      ! Gaussian temporal deposition
      adjtime = currenttime / ttm%tdepo - ttm%tcdepo
      deposit = (currenttime < 2.0_wp * ttm%tdepo * ttm%tcdepo)
      invbin = tstep / Real(redtstepmx, Kind=wp)
      adjeng = Exp(-0.5_wp * adjtime * adjtime)
    Case (TTM_TDEPO_EXP)
      ! decaying exponential temporal deposition
      adjtime = currenttime / ttm%tdepo
      deposit = (currenttime < ttm%tdepo * ttm%tcdepo)
      invbin = tstep / (ttm%tdepo * Real(redtstepmx, Kind=wp))
      adjeng = Exp(-adjtime)
    Case (TTM_TDEPO_DELTA)
      ! delta temporal deposition (over single diffusion timestep)
      deposit = (currenttime < tstep / (Real(redtstepmx, Kind=wp)))
      invbin = 1.0_wp
      adjeng = 1.0_wp
    Case (TTM_TDEPO_PULSE)
      ! pulse temporal deposition (over ttm%tdepo ps)
      deposit = (currenttime < ttm%tdepo)
      invbin = 1.0_wp
      adjeng = 1.0_wp
    End Select

    ! if (still) depositing energy, add to electronic temperature
    ! grid (active cells only) and ttm%adjust electronic temperatures
    ! accordingly

    If (deposit) Then
      ttm%lat_B(:, :, :) = ttm%lat_U(:, :, :) * ttm%norm * invbin * adjeng
      Select Case (ttm%CeType)
      Case (TTM_CE_CONST, TTM_CE_CONST_DYN)
        ! constant specific heat capacity
        Do k = 1, ttm%ntcell(3)
          Do j = 1, ttm%ntcell(2)
            Do i = 1, ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
              ttm%lat_I(i, j, k) = ttm%lat_I(i, j, k) + ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0)
              energy_diff = ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0) * ttm%rvolume * eV_to_kB
              If (energy_diff > zero_plus) Then
                start_Te = ttm%eltemp(ijk, 0, 0, 0)
                end_Te = start_Te + energy_diff / Ce0a
                ttm%eltemp(ijk, 0, 0, 0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (TTM_CE_TANH, TTM_CE_TANH_DYN)
        ! hyperbolic tangent specific heat capacity
        Do k = 1, ttm%ntcell(3)
          Do j = 1, ttm%ntcell(2)
            Do i = 1, ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
              ttm%lat_I(i, j, k) = ttm%lat_I(i, j, k) + ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0)
              energy_diff = ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0) * ttm%rvolume * eV_to_kB
              If (energy_diff > zero_plus) Then
                start_Te = ttm%eltemp(ijk, 0, 0, 0)
                increase = Cosh(ttm%sh_B * start_Te) * Exp(ttm%sh_B * energy_diff / sh_Aa)
                ! using equivalent function: Acosh(x)=Log(x+Sqrt((x-1.0)*(x+1.0)))
                end_Te = Log(increase + Sqrt((increase - 1.0_wp) * (increase + 1.0_wp))) / ttm%sh_B
                ttm%eltemp(ijk, 0, 0, 0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (TTM_CE_LINEAR, TTM_CE_LINEAR_DYN)
        ! linear specific heat capacity to Fermi temperature
        Do k = 1, ttm%ntcell(3)
          Do j = 1, ttm%ntcell(2)
            Do i = 1, ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
              ttm%lat_I(i, j, k) = ttm%lat_I(i, j, k) + ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0)
              energy_diff = ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0) * ttm%rvolume * eV_to_kB
              If (energy_diff > zero_plus) Then
                start_Te = ttm%eltemp(ijk, 0, 0, 0)
                If (start_Te >= ttm%Tfermi) Then
                  end_Te = start_Te + energy_diff / Cemaxa
                Else
                  end_Te = Sqrt(start_Te * start_Te + 2.0_wp * energy_diff * ttm%Tfermi / Cemaxa)
                  If (end_Te > ttm%Tfermi) end_Te = 0.5_wp * (start_Te * start_Te / ttm%Tfermi + ttm%Tfermi) + energy_diff / Cemaxa
                End If
                ttm%eltemp(ijk, 0, 0, 0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (TTM_CE_TABULATED, TTM_CE_TABULATED_DYN)
        ! tabulated ttm%volumetric heat capacity or more complex
        ! function: find new temperature iteratively by
        ! gradual integration (0.01 kelvin at a time)
        ! and interpolate over last 0.01 kelvin
        Do k = 1, ttm%ntcell(3)
          Do j = 1, ttm%ntcell(2)
            Do i = 1, ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
              ttm%lat_I(i, j, k) = ttm%lat_I(i, j, k) + ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0)
              energy_diff = ttm%lat_B(i, j, k) * ttm%act_ele_cell(ijk, 0, 0, 0) * ttm%rvolume * eV_to_kB
              ! work out change in electronic temperature for energy deposition
              ! (searching first in 0.01 kelvin increments, then interpolate based
              ! on constant heat capacity)
              If (energy_diff > zero_plus) Then
                start_Te = ttm%eltemp(ijk, 0, 0, 0)
                oldCe = Ce(start_Te, ttm)
                Do While (energy_diff > 0.0_wp)
                  newCe = Ce(start_Te + 0.01_wp, ttm)
                  increase = 0.005_wp * (oldCe + newCe)
                  energy_diff = energy_diff - increase
                  start_Te = start_Te + 0.01_wp
                  oldCe = newCe
                End Do
                energy_diff = energy_diff + increase
                start_Te = start_Te - 0.01_wp
                newCe = Ce(start_Te, ttm)
                end_Te = start_Te + 2.0_wp * energy_diff / (oldCe + newCe)
                ttm%eltemp(ijk, 0, 0, 0) = end_Te
              End If
            End Do
          End Do
        End Do

      Case Default
        call error(0, 'Unrecognised TTM heat capacity model')
      End Select

    Else

      ! if at end of deposition, find deposited energy at end of deposition

      lat_I_min = Minval(ttm%lat_I(1:ttm%ntcell(1), 1:ttm%ntcell(2), 1:ttm%ntcell(3)))
      lat_I_max = Maxval(ttm%lat_I(1:ttm%ntcell(1), 1:ttm%ntcell(2), 1:ttm%ntcell(3)))
      lat_I_sum = Sum(ttm%lat_I(1:ttm%ntcell(1), 1:ttm%ntcell(2), 1:ttm%ntcell(3)))
      Call gmin(comm, lat_I_min)
      Call gmax(comm, lat_I_max)
      Call gsum(comm, lat_I_sum)

      ! find energy input into electronic temperature system

      lat_U_min = Minval(ttm%lat_U(1:ttm%ntcell(1), 1:ttm%ntcell(2), 1:ttm%ntcell(3)))
      lat_U_max = Maxval(ttm%lat_U(1:ttm%ntcell(1), 1:ttm%ntcell(2), 1:ttm%ntcell(3)))
      lat_U_sum = Sum(ttm%lat_U(1:ttm%ntcell(1), 1:ttm%ntcell(2), 1:ttm%ntcell(3)))
      Call gmin(comm, lat_U_min)
      Call gmax(comm, lat_U_max)
      Call gsum(comm, lat_U_sum)

      ! check how closely two values match up: if error greater than
      ! tolerance, report discrepancy as warning

      If (comm%idnode == 0) Then
        If (Abs(lat_I_sum - lat_U_sum) > Abs(err_tol * lat_U_sum) .or. &
            Abs(lat_I_max - lat_U_max) > Abs(err_tol * lat_U_max) .or. &
            Abs(lat_I_min - lat_U_min) > Abs(err_tol * lat_U_min)) Then
          Call warning(530, Abs(lat_I_sum - lat_U_sum) / lat_U_sum * 100.0_wp, 0.0_wp, 0.0_wp)
        End If
      End If

      ! report successful completion of energy deposition

      If (currenttime < 1.0_wp) Then
        Write (number, '(f14.3)') currenttime * 1000.0_wp
        Write (messages(1), '(a,es12.5,3a)') &
          'electronic energy deposition of ', lat_I_sum, &
          ' eV completed successfully after ', Trim(Adjustl(number)), ' fs'
      Else
        Write (number, '(f14.6)') currenttime
        Write (messages(1), '(a,es12.5,3a)') &
          'electronic energy deposition of ', lat_I_sum, &
          ' eV completed successfully after ', Trim(Adjustl(number)), ' ps'
      End If
      Write (messages(2), '(a)') Repeat('-', 130)

      Call info(messages, 2, .true.)

      ! switch off tracking and deallocate arrays

      ttm%trackInit = .false.
      ttm%findepo = .true.

      Deallocate (ttm%lat_U, Stat=fail(1))
      Deallocate (ttm%lat_B, Stat=fail(2))
      Deallocate (ttm%lat_I, Stat=fail(3))

      If (Any(fail > 0)) Call error(1090)

    End If

  End Subroutine depoevolve

  Subroutine ttm_ion_temperature(ttm, thermo, domain, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating ion temperatures in simulation
    ! cells for two-temperature model
    !
    ! copyright - daresbury laboratory
    ! author    - s.l.darazewicz & m.a.seaton july 2012
    ! contrib   - g.khara may 2016
    ! contrib   - m.a.seaton september 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(ttm_type),           Intent(InOut) :: ttm
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Integer                                :: fail, i, ia, ii, ijk, ijk1, ijk2, ja, jj, ka, kk, &
                                              natmin
    Integer, Allocatable                   :: buf5(:), ijkatm(:), nat(:)
    Integer, Dimension(8)                  :: req
    Integer, Dimension(MPI_STATUS_SIZE, 8) :: stat
    Real(Kind=wp)                          :: crho, gsadd, tmp, velsq, vx, vy, vz, xxt, yyt, zzt
    Real(Kind=wp), Allocatable             :: buf1(:), buf2(:), buf3(:), buf4(:)

    ! allocate and zero arrays

    Allocate (nat(1:2 * ttm%numcell), ijkatm(1:config%natms), Stat=fail)
    If (fail > 0) Call error(1085)

    nat = 0
    ttm%tempion = 0.0_wp
    ttm%asource = 0.0_wp
    ttm%gsource = 0.0_wp
    ttm%ttmvom = 0.0_wp
    ijkatm = 0

    ! zero variable for dynamic cell density calculations

    crho = 0.0_wp

    ! if heterogeneous, ttm%gsource is an array containing no. of atoms in each cell
    ! otherwise it is an array containing the effective electron-phonon relaxation
    ! strength (known value for constant and homogeneous dynamic calculations)

    Select Case (ttm%gvar)
    Case (TTM_EPVAR_NULL, TTM_EPVAR_HOMO)
      gsadd = thermo%chi_ep
    Case (TTM_EPVAR_HETERO)
      gsadd = 1.0_wp
    End Select

    ! calculate overall momenta of cells for ionic temperature corrections

    If (comm%mxnode > 1) Then
      Do i = 1, config%natms
        xxt = config%parts(i)%xxx
        yyt = config%parts(i)%yyy
        zzt = config%parts(i)%zzz
        ia = Floor(xxt*ttm%grcell(1)+yyt*ttm%grcell(4)+zzt*ttm%grcell(7)+ttm%zerocell(1)) + 1
        ja = Floor(xxt*ttm%grcell(2)+yyt*ttm%grcell(5)+zzt*ttm%grcell(8)+ttm%zerocell(2)) + 1
        ka = Floor(xxt*ttm%grcell(3)+yyt*ttm%grcell(6)+zzt*ttm%grcell(9)+ttm%zerocell(3)) + 1

        ijk = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
        ijkatm(i) = ijk

        tmp = config%weight(i)
        If (config%lfrzn(i) == 0) Then
          ttm%ttmvom(ijk, 1) = ttm%ttmvom(ijk, 1) + tmp * config%vxx(i)
          ttm%ttmvom(ijk, 2) = ttm%ttmvom(ijk, 2) + tmp * config%vyy(i)
          ttm%ttmvom(ijk, 3) = ttm%ttmvom(ijk, 3) + tmp * config%vzz(i)
          ttm%ttmvom(ijk, 4) = ttm%ttmvom(ijk, 4) + tmp
        End If
      End Do

      Allocate (buf1(1:ttm%numcell), buf2(1:ttm%numcell), buf3(1:ttm%numcell), buf4(1:ttm%numcell), Stat=fail)
      If (fail > 0) Call error(1085)
      ! Sum up cell momenta and atomic masses in boundaries for ionic temperature corrections
      ! -z/+z directions
      buf1 = 0.0_wp
      buf2 = 0.0_wp
      buf3 = 0.0_wp
      buf4 = 0.0_wp
      ijk1 = 1
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * ttm%ntcell(3)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgz, domain%map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf4(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgz, domain%map(6), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf4(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      Do i = 1, ttm%numcell
        ttm%ttmvom(i, 1) = ttm%ttmvom(i, 1) + buf1(i)
        ttm%ttmvom(i, 2) = ttm%ttmvom(i, 2) + buf2(i)
        ttm%ttmvom(i, 3) = ttm%ttmvom(i, 3) + buf3(i)
        ttm%ttmvom(i, 4) = ttm%ttmvom(i, 4) + buf4(i)
      End Do
      ! -y/+y directions
      buf1 = 0.0_wp
      buf2 = 0.0_wp
      buf3 = 0.0_wp
      buf4 = 0.0_wp
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (2 * ttm%ntcell(2) + 2)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgy, domain%map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf4(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (2 * ttm%ntcell(2) + 3)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgy, domain%map(4), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf4(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      Do i = 1, ttm%numcell
        ttm%ttmvom(i, 1) = ttm%ttmvom(i, 1) + buf1(i)
        ttm%ttmvom(i, 2) = ttm%ttmvom(i, 2) + buf2(i)
        ttm%ttmvom(i, 3) = ttm%ttmvom(i, 3) + buf3(i)
        ttm%ttmvom(i, 4) = ttm%ttmvom(i, 4) + buf4(i)
      End Do
      ! -x/+x directions
      buf1 = 0.0_wp
      buf2 = 0.0_wp
      buf3 = 0.0_wp
      buf4 = 0.0_wp
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      ijk2 = 1 + ttm%ntcell(1) + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgx, domain%map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf4(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk2 = (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 4)
      ijk1 = 2 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgx, domain%map(2), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf4(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      Do i = 1, ttm%numcell
        ttm%ttmvom(i, 1) = ttm%ttmvom(i, 1) + buf1(i)
        ttm%ttmvom(i, 2) = ttm%ttmvom(i, 2) + buf2(i)
        ttm%ttmvom(i, 3) = ttm%ttmvom(i, 3) + buf3(i)
        ttm%ttmvom(i, 4) = ttm%ttmvom(i, 4) + buf4(i)
      End Do

      Deallocate (buf1, buf2, buf3, buf4, Stat=fail)
      If (fail > 0) Call error(1086)

      ! Send totals within domain cells to boundary halo voxels (needed to calculate correct
      ! velocities for particles in boundary halos, particularly if using VNL)
      ijk1 = 2 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      ijk2 = (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 4)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgx, domain%map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgx, domain%map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk1 = 1 + ttm%ntcell(1) + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgx, domain%map(2), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgx, domain%map(1), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (2 * ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgy, domain%map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgy, domain%map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (2 * ttm%ntcell(2) + 2)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgy, domain%map(4), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgy, domain%map(3), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgz, domain%map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgz, domain%map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * ttm%ntcell(3)
      ijk2 = 1
      Call MPI_ISEND(ttm%ttmvom(ijk1, 1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 2), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 3), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 3), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(ttm%ttmvom(ijk1, 4), 1, ttm%tmpmsgz, domain%map(6), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(ttm%ttmvom(ijk2, 4), 1, ttm%tmpmsgz, domain%map(5), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
    Else
      ! serial version: automatically corrects voxel to be within range for
      ! boundary halo/VNL-drifted particles
      Do i = 1, config%natms
        xxt = config%parts(i)%xxx
        yyt = config%parts(i)%yyy
        zzt = config%parts(i)%zzz
        ia = Modulo(Floor(xxt*ttm%grcell(1)+yyt*ttm%grcell(4)+zzt*ttm%grcell(7)+ttm%zerocell(1)), ttm%ntsys(1)) + 1
        ja = Modulo(Floor(xxt*ttm%grcell(2)+yyt*ttm%grcell(5)+zzt*ttm%grcell(8)+ttm%zerocell(2)), ttm%ntsys(2)) + 1
        ka = Modulo(Floor(xxt*ttm%grcell(3)+yyt*ttm%grcell(6)+zzt*ttm%grcell(9)+ttm%zerocell(3)), ttm%ntsys(3)) + 1
        ijk = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
        ijkatm(i) = ijk
        tmp = config%weight(i)
        If (config%lfrzn(i) == 0) Then
          ttm%ttmvom(ijk, 1) = ttm%ttmvom(ijk, 1) + tmp * config%vxx(i)
          ttm%ttmvom(ijk, 2) = ttm%ttmvom(ijk, 2) + tmp * config%vyy(i)
          ttm%ttmvom(ijk, 3) = ttm%ttmvom(ijk, 3) + tmp * config%vzz(i)
          ttm%ttmvom(ijk, 4) = ttm%ttmvom(ijk, 4) + tmp
        End If
      End Do

    End If

    ! calculate cell velocities

    Do i = 1, ttm%numcell
      If (ttm%ttmvom(i, 4) > zero_plus) Then
        ttm%ttmvom(i, 1:3) = ttm%ttmvom(i, 1:3) / ttm%ttmvom(i, 4)
      Else
        ttm%ttmvom(i, 1:3) = 0.0_wp
      End If
    End Do

    ! calculate ionic temperatures (accounting for cell velocities)
    ! and source terms: electron-phonon (ttm%gsource) and electronic
    ! stopping (ttm%asource)

    Do i = 1, config%natms

      ijk = ijkatm(i)

      vx = config%vxx(i) - ttm%ttmvom(ijk, 1)
      vy = config%vyy(i) - ttm%ttmvom(ijk, 2)
      vz = config%vzz(i) - ttm%ttmvom(ijk, 3)
      velsq = vx * vx + vy * vy + vz * vz
      tmp = config%weight(i)

      ttm%tempion(ijk) = ttm%tempion(ijk) + tmp * velsq

      ttm%gsource(ijk) = ttm%gsource(ijk) + gsadd

      nat(2 * ijk - 1) = nat(2 * ijk - 1) + 1

      If ((velsq > thermo%vel_es2) .and. (thermo%chi_es > zero_plus)) Then
        ttm%asource(ijk) = ttm%asource(ijk) + tmp * velsq
        nat(2 * ijk) = nat(2 * ijk) + 1
      End If

    End Do

    If (comm%mxnode > 1) Then
      Allocate (buf1(1:ttm%numcell), buf2(1:ttm%numcell), buf3(1:ttm%numcell), buf5(1:2 * ttm%numcell), Stat=fail)
      If (fail > 0) Call error(1085)
      ! Sum up boundary values of nat, tempion, asource and gsource within main grid only
      ! -z/+z directions
      buf1 = 0.0_wp
      buf2 = 0.0_wp
      buf3 = 0.0_wp
      buf5 = 0
      ijk1 = 1
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * (ttm%ntcell(3))
      Call MPI_ISEND(ttm%tempion(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%gsource(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%asource(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(nat(2 * ijk1 - 1), 1, ttm%nummsgz, domain%map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf5(2 * ijk2 - 1), 1, ttm%nummsgz, domain%map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      Call MPI_ISEND(ttm%tempion(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%gsource(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%asource(ijk2), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk1), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(nat(2 * ijk2 - 1), 1, ttm%nummsgz, domain%map(6), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf5(2 * ijk1 - 1), 1, ttm%nummsgz, domain%map(5), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ttm%tempion = ttm%tempion + buf1
      ttm%gsource = ttm%gsource + buf2
      ttm%asource = ttm%asource + buf3
      nat = nat + buf5
      ! -y/+y directions
      buf1 = 0.0_wp
      buf2 = 0.0_wp
      buf3 = 0.0_wp
      buf5 = 0
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (2 * ttm%ntcell(2) + 2)
      Call MPI_ISEND(ttm%tempion(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%gsource(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%asource(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(nat(2 * ijk1 - 1), 1, ttm%nummsgy, domain%map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf5(2 * ijk2 - 1), 1, ttm%nummsgy, domain%map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (2 * ttm%ntcell(2) + 3)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%tempion(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%gsource(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%asource(ijk2), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk1), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(nat(2 * ijk2 - 1), 1, ttm%nummsgy, domain%map(4), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf5(2 * ijk1 - 1), 1, ttm%nummsgy, domain%map(3), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ttm%tempion = ttm%tempion + buf1
      ttm%gsource = ttm%gsource + buf2
      ttm%asource = ttm%asource + buf3
      nat = nat + buf5
      ! -x/+x directions
      buf1 = 0.0_wp
      buf2 = 0.0_wp
      buf3 = 0.0_wp
      buf5 = 0
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ijk2 = 1 + ttm%ntcell(1) + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      Call MPI_ISEND(ttm%tempion(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%gsource(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%asource(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(nat(2 * ijk1 - 1), 1, ttm%nummsgx, domain%map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf5(2 * ijk2 - 1), 1, ttm%nummsgx, domain%map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ijk2 = (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 4)
      ijk1 = 2 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 3)
      Call MPI_ISEND(ttm%tempion(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(buf1(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
      Call MPI_ISEND(ttm%gsource(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(buf2(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_ISEND(ttm%asource(ijk2), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(5), comm%ierr)
      Call MPI_IRECV(buf3(ijk1), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(6), comm%ierr)
      Call MPI_ISEND(nat(2 * ijk2 - 1), 1, ttm%nummsgx, domain%map(2), Grid4_tag, comm%comm, req(7), comm%ierr)
      Call MPI_IRECV(buf5(2 * ijk1 - 1), 1, ttm%nummsgx, domain%map(1), Grid4_tag, comm%comm, req(8), comm%ierr)
      Call MPI_WAITALL(8, req, stat, comm%ierr)
      ttm%tempion = ttm%tempion + buf1
      ttm%gsource = ttm%gsource + buf2
      ttm%asource = ttm%asource + buf3
      nat = nat + buf5

      Deallocate (buf1, buf2, buf3, buf5, Stat=fail)
      If (fail > 0) Call error(1086)
    End If

    natmin = Merge(0, ttm%amin - 1, ttm%trackInit)
    ttm%old_ele_cell = ttm%act_ele_cell
    ttm%act_ele_cell = 1.0_wp
    ttm%acell_old = ttm%acell
    ttm%acell = 0

    ! loop through ionic temperature cells in current node
    Do ka = 1, ttm%ntcell(3)
      Do ja = 1, ttm%ntcell(2)
        Do ia = 1, ttm%ntcell(1)
          ijk = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
          ! calculate ionic temperature for all cells with at least
          ! minimum number of particles (1 during deposition, ttm%amin
          ! at all other times), calculating dynamic cell density
          ! (if required) from active cells, removing centre of mass
          ! motion and determining any inactive ionic temperature cells
          If (nat(2 * ijk - 1) > natmin .and. nat(2 * ijk - 1) > 1) Then
            ttm%tempion(ijk) = ttm%tempion(ijk) / (3.0_wp * boltz * Real(nat(2 * ijk - 1), Kind=wp))
            ttm%acell = ttm%acell + 1
            crho = crho + ttm%gsource(ijk)
          Else If (natmin == 0 .and. nat(2 * ijk - 1) == 1) Then
            vx = 0.5_wp * ttm%ttmvom(ijk, 1)
            vy = 0.5_wp * ttm%ttmvom(ijk, 2)
            vz = 0.5_wp * ttm%ttmvom(ijk, 3)
            velsq = ttm%ttmvom(ijk, 4) * (vx * vx + vy * vy + vz * vz)
            ttm%tempion(ijk) = velsq / (3.0_wp * boltz)
            ttm%acell = ttm%acell + 1
            crho = crho + ttm%gsource(ijk)
          Else
            ttm%tempion(ijk) = 0.0_wp
            ttm%act_ele_cell(ijk, 0, 0, 0) = 0.0_wp
          End If
          ! calculate electronic stopping terms (if more than one atom with speed > vel_cs)
          If (nat(2 * ijk) > 0) Then
            ttm%asource(ijk) = ttm%asource(ijk) * thermo%chi_es / boltz
          Else
            ttm%asource(ijk) = 0.0_wp
          End If
          ! calculate electron-phonon coupling terms
          ttm%gsource(ijk) = ttm%gsource(ijk) * 3.0_wp
        End Do
      End Do
    End Do

    ! communicate and work out active ionic temperature cells in boundary halos

    If (comm%mxnode > 1) Then
      Call gsum(comm, ttm%acell)
      Call gsum(comm, crho)
      ijk1 = 2 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ijk2 = 1 + (ttm%ntcell(1) + 1) + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ii = Merge(-1, 0, (domain%idx == domain%nx - 1))
      Call MPI_ISEND(ttm%act_ele_cell(ijk1, 0, 0, 0), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%act_ele_cell(ijk2, ii, 0, 0), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1)) + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ii = Merge(1, 0, (domain%idx == 0))
      Call MPI_ISEND(ttm%act_ele_cell(ijk1, 0, 0, 0), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%act_ele_cell(ijk2, ii, 0, 0), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_WAITALL(4, req, stat, comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2) + 2))
      jj = Merge(-1, 0, (domain%idy == domain%ny - 1))
      Call MPI_ISEND(ttm%act_ele_cell(ijk1, 0, 0, 0), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%act_ele_cell(ijk2, 0, jj, 0), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + (ttm%ntcell(2) + 2))
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      jj = Merge(1, 0, (domain%idy == 0))
      Call MPI_ISEND(ttm%act_ele_cell(ijk1, 0, 0, 0), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%act_ele_cell(ijk2, 0, jj, 0), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_WAITALL(4, req, stat, comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * ((ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1))
      kk = Merge(-1, 0, (domain%idz == domain%nz - 1))
      Call MPI_ISEND(ttm%act_ele_cell(ijk1, 0, 0, 0), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
      Call MPI_IRECV(ttm%act_ele_cell(ijk2, 0, 0, kk), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * ((ttm%ntcell(2) + 2) * ttm%ntcell(3))
      ijk2 = 1
      kk = Merge(1, 0, (domain%idz == 0))
      Call MPI_ISEND(ttm%act_ele_cell(ijk1, 0, 0, 0), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
      Call MPI_IRECV(ttm%act_ele_cell(ijk2, 0, 0, kk), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
      Call MPI_WAITALL(4, req, stat, comm%ierr)
    Else
      Do ka = 1, ttm%ntcell(3)
        Do ja = 1, ttm%ntcell(2)
          ijk1 = 2 + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
          ijk2 = 1 + (ttm%ntcell(1) + 1) + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
          ttm%act_ele_cell(ijk2, -1, 0, 0) = ttm%act_ele_cell(ijk1, 0, 0, 0)
          ijk1 = 1 + ttm%ntcell(1) + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
          ijk2 = 1 + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
          ttm%act_ele_cell(ijk2, 1, 0, 0) = ttm%act_ele_cell(ijk1, 0, 0, 0)
        End Do
      End Do
      Do ka = 1, ttm%ntcell(3)
        Do ia = 0, ttm%ntcell(1) + 1
          ijk1 = 1 + ia + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2) * ka)
          ijk2 = 1 + ia + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2) + 2) * ka)
          ttm%act_ele_cell(ijk2, 0, -1, 0) = ttm%act_ele_cell(ijk1, 0, 0, 0)
          ijk1 = 1 + ia + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + (ttm%ntcell(2) + 2) * ka)
          ijk2 = 1 + ia + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * ka
          ttm%act_ele_cell(ijk2, 0, 1, 0) = ttm%act_ele_cell(ijk1, 0, 0, 0)
        End Do
      End Do
      Do ja = 0, ttm%ntcell(2) + 1
        Do ia = 0, ttm%ntcell(1) + 1
          ijk1 = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2))
          ijk2 = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1))
          ttm%act_ele_cell(ijk2, 0, 0, -1) = ttm%act_ele_cell(ijk1, 0, 0, 0)
          ijk1 = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ttm%ntcell(3))
          ijk2 = 1 + ia + (ttm%ntcell(1) + 2) * ja
          ttm%act_ele_cell(ijk2, 0, 0, 1) = ttm%act_ele_cell(ijk1, 0, 0, 0)
        End Do
      End Do
    End If

    ! cell velocities generally used to correct velocities
    ! used in inhomogeneous Langevin thermostat: zero these
    ! velocities (x- and y-components only) if the user
    ! says otherwise (only for z-component)

    If (.not. ttm%ttmthvel) Then
      ttm%ttmvom = 0.0_wp
    Else If (ttm%ttmthvelz) Then
      ttm%ttmvom(1:ttm%numcell, 1:2) = 0.0_wp
    End If

    ! dynamically calculate cell density for active cells
    ! if requested by user (after deposition stage)

    If (ttm%ttmdyndens .and. ttm%findepo) Then
      Select Case (ttm%gvar)
      Case (TTM_EPVAR_NULL, TTM_EPVAR_HOMO)
        ttm%cellrho = crho / (Real(ttm%acell, Kind=wp) * thermo%chi_ep * ttm%volume)
      Case (TTM_EPVAR_HETERO)
        ttm%cellrho = crho / (Real(ttm%acell, Kind=wp) * ttm%volume)
      End Select
      If (ttm%cellrho > zero_plus) Then
        ttm%rcellrho = 1.0_wp / ttm%cellrho
      Else
        ttm%rcellrho = 0.0_wp
      End If
    End If

    ! optional: send ttm%asource terms back to boundary voxels (only needed in +x, +y, +z directions)

    !  If (mxnode>1) Then
    !    ijk1 = 1 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    !    ijk2 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    !    Call MPI_ISEND (ttm%asource(ijk1), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
    !    Call MPI_IRECV (ttm%asource(ijk2), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
    !    Call MPI_WAITALL (2, req, stat, comm%ierr)
    !    ijk1 = 1 + (ttm%ntcell(1)+2) * (2*ttm%ntcell(2) + 2)
    !    ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
    !    Call MPI_ISEND (ttm%asource(ijk1), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
    !    Call MPI_IRECV (ttm%asource(ijk2), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
    !    Call MPI_WAITALL (2, req, stat, comm%ierr)
    !    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2) * (ttm%ntcell(3))
    !    ijk2 = 1
    !    Call MPI_ISEND (ttm%asource(ijk1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
    !    Call MPI_IRECV (ttm%asource(ijk2), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
    !    Call MPI_WAITALL (2, req, stat, comm%ierr)
    !  End If

    Deallocate (nat, ijkatm, Stat=fail)
    If (fail > 0) Call error(1086)

  End Subroutine ttm_ion_temperature

  Subroutine ttm_thermal_diffusion(tstep, time, nstep, nsteql, nstbpo, ndump, nstrun, &
                                   ttm, thermo, domain, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for iterating thermal diffusion equation for
    ! electrons using numerical integration, after a.rutherford and
    ! d.duffy
    !
    ! copyright - daresbury laboratory
    ! author    - s.l.darazewicz & m.a.seaton july 2012
    ! contrib -   g.khara may 2016
    ! contrib -   g.khara, s.t.murphy & m.a.seaton september 2017
    ! contrib -   m.a.seaton february 2020
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),         Intent(In   ) :: tstep, time
    Integer,               Intent(In   ) :: nstep, nsteql, nstbpo, ndump, nstrun
    Type(ttm_type),        Intent(InOut) :: ttm
    Type(thermostat_type), Intent(In   ) :: thermo
    Type(domains_type),    Intent(In   ) :: domain
    Type(comms_type),      Intent(InOut) :: comm

    Character(Len=STR_LEN)         :: messages(6)
    Integer                    :: fail, i, ii, ijk, j, jj, k, kk, redtstep, redtstepmx, &
                                  ixp, ixm, iyp, iym, izp, izm, ixpiyp, ixpiym, ixmiyp, ixmiym, &
                                  ixpizp, ixpizm, ixmizp, ixmizm, iypizp, iypizm, iymizp, iymizm
    Logical                    :: debug1 = .false., safe
    Real(Kind=wp)              :: actsite, actxm, actxp, actym, actyp, actzm, actzp, &
                                  actxpyp, actxpym, actxmyp, actxmym, actxpzp, actxpzm, &
                                  actxmzp, actxmzm, actypzp, actypzm, actymzp, actymzm, alploc, &
                                  del2av, delx2, dely2, delz2, delxy, delxz, delyz, &
                                  eltempKe, eltempmax, eltempmaxKe, eltempmean, eltempmin, eltempminKe, &
                                  fomAx, fomAy, fomAz, fomAxy, fomAxz, fomAyz, &
                                  fopttstep, maxtstep, mintstep, opttstep, temp, &
                                  d2tdu2, d2tdv2, d2tdw2, d2tdudv, d2tdudw, d2tdvdw
    Real(Kind=wp), Allocatable :: eltemp1(:, :, :, :)

! Debugging flag

    ! Initialise eltemp1 (electronic temperature grid for next timestep) and timestep sizes

 Allocate (eltemp1(1:ttm%numcell, -ttm%eltcell(1):ttm%eltcell(1), -ttm%eltcell(2):ttm%eltcell(2), -ttm%eltcell(3):ttm%eltcell(3)), &
              Stat=fail)
    If (fail > 0) Call error(1087)
    eltemp1 = 0.0_wp
    redtstepmx = 1

    ! Initialise temp
    temp = thermo%temp

    ! deposition stage 1 (initialization):
    ! nstep-nsteql offsets equilibration time

    If ((nstep - nsteql) == 1 .and. (ttm%sdepoType /= TTM_SDEPO_NULL &
         .and. (ttm%dEdX > zero_plus .or. ttm%fluence > zero_plus))) Then
      Call depoinit(time, ttm, comm)
    End If

    ! determine timestep reduction factor (chosen empirically, acts beyond minimum stability condition)

    fopttstep = 0.25_wp

    ! determine maximum/minimum spacings and electronic temperatures

    delx2 = ttm%delu * ttm%delu
    dely2 = ttm%delv * ttm%delv
    delz2 = ttm%delw * ttm%delw
    delxy = ttm%delu * ttm%delv
    delxz = ttm%delu * ttm%delw
    delyz = ttm%delv * ttm%delw
    del2av = (delx2 * dely2 * delz2) / (dely2 * delz2 + delx2 * delz2 + delx2 * dely2)
    Call eltemp_max(eltempmax, ttm, comm)
    Call eltemp_min(eltempmin, ttm, comm)

    ! This section of the code establishes the optimum size of fourier mesh to ensure stability of the electronic
    ! temperature finite difference solver

    Select Case (ttm%KeType)
    Case (TTM_KE_INFINITE)
      ! infinite thermal conductivity
      redtstepmx = 1
      opttstep = tstep
    Case (TTM_KE_CONST)
      ! constant thermal conductivity and non-metal systems
      mintstep = 0.5_wp * del2av / alp(eltempmax, ttm)
      maxtstep = 0.5_wp * del2av / alp(eltempmin, ttm)
      opttstep = fopttstep * Min(mintstep, maxtstep)
      redtstepmx = Max(Ceiling(tstep / opttstep), 2)
    Case (TTM_KE_DRUDE)
      ! Drude-type thermal conductivity
      mintstep = 0.5_wp * del2av * Ce(eltempmax, ttm) / KeD(eltempmax, temp, ttm)
      maxtstep = 0.5_wp * del2av * Ce(eltempmin, ttm) / KeD(eltempmin, temp, ttm)
      opttstep = fopttstep * Min(mintstep, maxtstep)
      redtstepmx = Max(Ceiling(tstep / opttstep), 2)
    Case (TTM_KE_TABULATED)
      Call eltemp_maxKe(temp, eltempmaxKe, ttm, comm)
      Call eltemp_minke(temp, eltempminKe, ttm, comm)
      ! tabulated thermal conductivity
      mintstep = 0.5_wp * del2av * Ce(eltempmax, ttm) / Ke(eltempmaxKe, ttm)
      maxtstep = 0.5_wp * del2av * Ce(eltempmin, ttm) / Ke(eltempminKe, ttm)
      opttstep = fopttstep * Min(mintstep, maxtstep)
      redtstepmx = Max(Ceiling(tstep / opttstep), 2)
    End Select

    ! reduce timestep further for deposition stage

    If (ttm%KeType /= TTM_KE_INFINITE .and. ttm%trackInit) Then
      Select Case (ttm%tdepoType)
      Case (TTM_TDEPO_GAUSS)
        redtstepmx = Max(50, redtstepmx)
      Case (TTM_TDEPO_EXP)
        redtstepmx = Max(10000, redtstepmx)
      Case Default
        redtstepmx = Max(2, redtstepmx)
      End Select
    End If

    fomAx = tstep / (delx2 * Real(redtstepmx, Kind=wp))
    fomAy = tstep / (dely2 * Real(redtstepmx, Kind=wp))
    fomAz = tstep / (delz2 * Real(redtstepmx, Kind=wp))
    fomAxy = 0.25_wp * tstep / (delxy * Real(redtstepmx, Kind=wp))
    fomAxz = 0.25_wp * tstep / (delxz * Real(redtstepmx, Kind=wp))
    fomAyz = 0.25_wp * tstep / (delyz * Real(redtstepmx, Kind=wp))

    ! write information to OUTPUT

    If (Mod(nstep, nstbpo) == 0 .or. nstep == 1) Then
      Write (messages(1), '(a)') 'ttm thermal diffusion timesteps:'
      Write (messages(2), '(4x,a,3x,a,5x,a)') 'optimal/ps', 'actual/ps', 'diff/md'
      Write (messages(3), '(2x,2es12.4,2x,i10)') opttstep, tstep / Real(redtstepmx, Kind=wp), redtstepmx
      If (ttm%ttmdyndens) Then
        Write (messages(4), '(a)') 'active ion temperature cells:'
        Write (messages(5), '(4x,a,2x,a)') 'atom dens.', 'no. of active cells'
        Write (messages(6), '(2x,es12.4,11x,i10)') ttm%cellrho, ttm%acell
        Call info(messages, 6, .true.)
      Else
        Call info(messages, 3, .true.)
      End If
    End If

    ! apply boundary conditions

    Call boundaryHalo(ttm, domain, comm)
    Call boundaryCond(ttm%bcTypeE, temp, ttm, comm)

    ! print statistics to files: electronic and ionic temperatures
    ! (note timestep is subtracted by 1, as these are values at
    !  beginning of MD timestep)

    Call printElecLatticeStatsToFile('PEAK_E', time, temp, nstep - 1, ttm%ttmstats, ttm, comm)
    Call peakProfilerElec('LATS_E', nstep - 1, ttm%ttmtraj, ttm, comm)

    Call printLatticeStatsToFile(ttm%tempion, 'PEAK_I', time, nstep - 1, ttm%ttmstats, ttm, comm)
    Call peakProfiler(ttm%tempion, 'LATS_I', nstep - 1, ttm%ttmtraj, ttm, comm)

    ! debugging option: print electron-phonon and electronic stopping source terms
    !                   (ttm%normally switched off)

    If (debug1) Then
      Call printLatticeStatsToFile(ttm%gsource, 'PEAK_G', time, nstep - 1, ttm%ttmstats, ttm, comm)
      Call peakProfiler(ttm%gsource, 'LATS_G', nstep - 1, ttm%ttmtraj, ttm, comm)
      Call printLatticeStatsToFile(ttm%asource, 'PEAK_A', time, nstep - 1, ttm%ttmstats, ttm, comm)
      Call peakProfiler(ttm%asource, 'LATS_A', nstep - 1, ttm%ttmtraj, ttm, comm)
    End If

    safe = .true.

    ! determine energy redistribution from deactivated ionic temperature voxels for slab geometry

    If (ttm%redistribute) Then
      Call redistribute_te(temp, ttm, domain, comm)
    End If

    ! Adaptive timestep

    Do redtstep = 1, redtstepmx

      ! deposition stage 2 (with boundary conditions)

      If (ttm%trackInit) Then
        Call depoevolve(time, tstep, redtstep, redtstepmx, ttm, comm)
        Call boundaryCond(ttm%bcTypeE, temp, ttm, comm)
      End If

      ! MAIN LOOP
      ! this portion of the code is the main electronic temperature solver, solving
      ! for heat diffusion with two source terms within ionic cells, i.e.
      !
      !   Ce*d(T_e)/dt = d/dx(Ke * d(T_e)/dx) + ttm%volume*(ttm%asource + ttm%gsource*(T_e-T_i))
      !
      ! the partial differential equation is solved using an explicit finite
      ! difference solver: care is needed in choosing timestep to ensure
      ! numerical stability

      Select Case (ttm%KeType)
      Case (TTM_KE_INFINITE)
        ! infinite thermal conductivity case: set all electronic temperatures
        ! to mean value in active cells, to system temperature in inactive cells
        Call eltemp_mean(eltempmean, ttm, comm)
        eltemp1 = eltempmean
        Do ijk = 1, ttm%numcell
          If (ttm%act_ele_cell(ijk, 0, 0, 0) <= zero_plus) eltemp1(ijk, 0, 0, 0) = temp
        End Do

      Case (TTM_KE_CONST)
        ! constant thermal conductivity or non-metal case
        If (ttm%redistribute) Then
          ! system with cell deactivation/energy redistribution
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)

                If (ii > -2 .and. ii < 2 .and. jj > -2 .and. jj < 2 .and. kk > -2 .and. kk < 2) Then
                  ! replace electronic temperatures with values required for energy redistribution
                  Do ijk = 1, ttm%numcell
                    If (ttm%adjust(ijk, ii, jj, kk)) ttm%eltemp(ijk, ii, jj, kk) = ttm%eltemp_adj(ijk, ii, jj, kk)
                  End Do
                  ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                  Do k = 1, ttm%ntcell(3)
                    Do j = 1, ttm%ntcell(2)
                      Do i = 1, ttm%ntcell(1)
                        ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                        ixm = ijk - 1
                        ixp = ijk + 1
                        iym = ijk - (ttm%ntcell(1) + 2)
                        iyp = ijk + (ttm%ntcell(1) + 2)
                        izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        ixpiyp = ijk + ttm%ntcell(1) + 3
                        ixpiym = ijk - ttm%ntcell(1) - 1
                        ixmiyp = ijk + ttm%ntcell(1) + 1
                        ixmiym = ijk - ttm%ntcell(1) - 3
                        ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        actsite = ttm%act_ele_cell(ijk, ii, jj, kk)
                        actxm = actsite * ttm%act_ele_cell(ixm, ii, jj, kk)
                        actxp = actsite * ttm%act_ele_cell(ixp, ii, jj, kk)
                        actym = actsite * ttm%act_ele_cell(iym, ii, jj, kk)
                        actyp = actsite * ttm%act_ele_cell(iyp, ii, jj, kk)
                        actzm = actsite * ttm%act_ele_cell(izm, ii, jj, kk)
                        actzp = actsite * ttm%act_ele_cell(izp, ii, jj, kk)
                        actxpyp = actsite * ttm%act_ele_cell(ixpiyp, ii, jj, kk)
                        actxpym = actsite * ttm%act_ele_cell(ixpiym, ii, jj, kk)
                        actxmyp = actsite * ttm%act_ele_cell(ixmiyp, ii, jj, kk)
                        actxmym = actsite * ttm%act_ele_cell(ixmiym, ii, jj, kk)
                        actxpzp = actsite * ttm%act_ele_cell(ixpizp, ii, jj, kk)
                        actxpzm = actsite * ttm%act_ele_cell(ixpizm, ii, jj, kk)
                        actxmzp = actsite * ttm%act_ele_cell(ixmizp, ii, jj, kk)
                        actxmzm = actsite * ttm%act_ele_cell(ixmizm, ii, jj, kk)
                        actypzp = actsite * ttm%act_ele_cell(iypizp, ii, jj, kk)
                        actypzm = actsite * ttm%act_ele_cell(iypizm, ii, jj, kk)
                        actymzp = actsite * ttm%act_ele_cell(iymizp, ii, jj, kk)
                        actymzm = actsite * ttm%act_ele_cell(iymizm, ii, jj, kk)
                        alploc = alp(ttm%eltemp(ijk, ii, jj, kk), ttm)
                        d2tdu2 = fomAx * alploc * (actxm * (ttm%eltemp(ixm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                                   actxp * (ttm%eltemp(ixp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdv2 = fomAy * alploc * (actym * (ttm%eltemp(iym, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                                   actyp * (ttm%eltemp(iyp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdw2 = fomAz * alploc * (actzm * (ttm%eltemp(izm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                                   actzp * (ttm%eltemp(izp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdudv = fomAxy * alploc * (actxpyp * ttm%eltemp(ixpiyp, ii, jj, kk) - &
                                                     actxpym * ttm%eltemp(ixpiym, ii, jj, kk) - &
                                                     actxmyp * ttm%eltemp(ixmiyp, ii, jj, kk) + &
                                                     actxmym * ttm%eltemp(ixmiym, ii, jj, kk))
                        d2tdudw = fomAxz * alploc * (actxpzp * ttm%eltemp(ixpizp, ii, jj, kk) - &
                                                     actxpzm * ttm%eltemp(ixpizm, ii, jj, kk) - &
                                                     actxmzp * ttm%eltemp(ixmizp, ii, jj, kk) + &
                                                     actxmzm * ttm%eltemp(ixmizm, ii, jj, kk))
                        d2tdvdw = fomAyz * alploc * (actypzp * ttm%eltemp(iypizp, ii, jj, kk) - &
                                                     actypzm * ttm%eltemp(iypizm, ii, jj, kk) - &
                                                     actymzp * ttm%eltemp(iymizp, ii, jj, kk) + &
                                                     actymzm * ttm%eltemp(iymizm, ii, jj, kk))
                        eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                               + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                               + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                      End Do
                    End Do
                  End Do
                Else
                  ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                  Do k = 1, ttm%ntcell(3)
                    Do j = 1, ttm%ntcell(2)
                      Do i = 1, ttm%ntcell(1)
                        ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                        ixm = ijk - 1
                        ixp = ijk + 1
                        iym = ijk - (ttm%ntcell(1) + 2)
                        iyp = ijk + (ttm%ntcell(1) + 2)
                        izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        ixpiyp = ijk + ttm%ntcell(1) + 3
                        ixpiym = ijk - ttm%ntcell(1) - 1
                        ixmiyp = ijk + ttm%ntcell(1) + 1
                        ixmiym = ijk - ttm%ntcell(1) - 3
                        ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        alploc = alp(ttm%eltemp(ijk, ii, jj, kk), ttm)
                        d2tdu2 = fomAx * alploc * (ttm%eltemp(ixm, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk) - &
                                                                        2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                        d2tdv2 = fomAy * alploc * (ttm%eltemp(iym, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk) - &
                                                                        2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                        d2tdw2 = fomAz * alploc * (ttm%eltemp(izm, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk) - &
                                                                        2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                        d2tdudv = fomAxy * alploc * (ttm%eltemp(ixpiyp, ii, jj, kk) - ttm%eltemp(ixpiym, ii, jj, kk) - &
                                                     ttm%eltemp(ixmiyp, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk))
                        d2tdudw = fomAxz * alploc * (ttm%eltemp(ixpizp, ii, jj, kk) - ttm%eltemp(ixpizm, ii, jj, kk) - &
                                                     ttm%eltemp(ixmizp, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk))
                        d2tdvdw = fomAyz * alploc * (ttm%eltemp(iypizp, ii, jj, kk) - ttm%eltemp(iypizm, ii, jj, kk) - &
                                                     ttm%eltemp(iymizp, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk))
                        eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                               + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                               + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                      End Do
                    End Do
                  End Do
                End If

              End Do
            End Do
          End Do
        Else
          ! standard thermal diffusion calculation applies when energy redistribution is not applicable
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)
                Do k = 1, ttm%ntcell(3)
                  Do j = 1, ttm%ntcell(2)
                    Do i = 1, ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                      ixm = ijk - 1
                      ixp = ijk + 1
                      iym = ijk - (ttm%ntcell(1) + 2)
                      iyp = ijk + (ttm%ntcell(1) + 2)
                      izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                      izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                      ixpiyp = ijk + ttm%ntcell(1) + 3
                      ixpiym = ijk - ttm%ntcell(1) - 1
                      ixmiyp = ijk + ttm%ntcell(1) + 1
                      ixmiym = ijk - ttm%ntcell(1) - 3
                      ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                      ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                      ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                      ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                      iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                      iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                      iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                      iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                      alploc = alp(ttm%eltemp(ijk, ii, jj, kk), ttm)
                      d2tdu2 = fomAx * alploc * (ttm%eltemp(ixm, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk) - &
                                                                      2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                      d2tdv2 = fomAy * alploc * (ttm%eltemp(iym, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk) - &
                                                                      2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                      d2tdw2 = fomAz * alploc * (ttm%eltemp(izm, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk) - &
                                                                      2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                      d2tdudv = fomAxy * alploc * (ttm%eltemp(ixpiyp, ii, jj, kk) - ttm%eltemp(ixpiym, ii, jj, kk) - &
                                                   ttm%eltemp(ixmiyp, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk))
                      d2tdudw = fomAxz * alploc * (ttm%eltemp(ixpizp, ii, jj, kk) - ttm%eltemp(ixpizm, ii, jj, kk) - &
                                                   ttm%eltemp(ixmizp, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk))
                      d2tdvdw = fomAyz * alploc * (ttm%eltemp(iypizp, ii, jj, kk) - ttm%eltemp(iypizm, ii, jj, kk) - &
                                                   ttm%eltemp(iymizp, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk))
                      eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                             + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                             + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        End If

      Case (TTM_KE_DRUDE)
        ! Drude-type thermal conductivity case
        If (ttm%redistribute) Then
          ! system with cell deactivation/energy redistribution
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)

                If (ii > -2 .and. ii < 2 .and. jj > -2 .and. jj < 2 .and. kk > -2 .and. kk < 2) Then
                  ! replace electronic temperatures with values required for energy redistribution
                  Do ijk = 1, ttm%numcell
                    If (ttm%adjust(ijk, ii, jj, kk)) ttm%eltemp(ijk, ii, jj, kk) = ttm%eltemp_adj(ijk, ii, jj, kk)
                  End Do
                  ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                  Do k = 1, ttm%ntcell(3)
                    Do j = 1, ttm%ntcell(2)
                      Do i = 1, ttm%ntcell(1)
                        ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                        ixm = ijk - 1
                        ixp = ijk + 1
                        iym = ijk - (ttm%ntcell(1) + 2)
                        iyp = ijk + (ttm%ntcell(1) + 2)
                        izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        ixpiyp = ijk + ttm%ntcell(1) + 3
                        ixpiym = ijk - ttm%ntcell(1) - 1
                        ixmiyp = ijk + ttm%ntcell(1) + 1
                        ixmiym = ijk - ttm%ntcell(1) - 3
                        ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        actsite = ttm%act_ele_cell(ijk, ii, jj, kk) / Ce(ttm%eltemp(ijk, ii, jj, kk), ttm)
                        actxm = actsite * ttm%act_ele_cell(ixm, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixm, ii, jj, kk)), temp, ttm)
                        actxp = actsite * ttm%act_ele_cell(ixp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk)), temp, ttm)
                        actym = actsite * ttm%act_ele_cell(iym, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iym, ii, jj, kk)), temp, ttm)
                        actyp = actsite * ttm%act_ele_cell(iyp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk)), temp, ttm)
                        actzm = actsite * ttm%act_ele_cell(izm, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(izm, ii, jj, kk)), temp, ttm)
                        actzp = actsite * ttm%act_ele_cell(izp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk)), temp, ttm)
                        actxpyp = actsite * ttm%act_ele_cell(ixpiyp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpiyp, ii, jj, kk)), temp, ttm)
                        actxpym = actsite * ttm%act_ele_cell(ixpiym, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpiym, ii, jj, kk)), temp, ttm)
                        actxmyp = actsite * ttm%act_ele_cell(ixmiyp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmiyp, ii, jj, kk)), temp, ttm)
                        actxmym = actsite * ttm%act_ele_cell(ixmiym, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk)), temp, ttm)
                        actxpzp = actsite * ttm%act_ele_cell(ixpizp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpizp, ii, jj, kk)), temp, ttm)
                        actxpzm = actsite * ttm%act_ele_cell(ixpizm, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpizm, ii, jj, kk)), temp, ttm)
                        actxmzp = actsite * ttm%act_ele_cell(ixmizp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmizp, ii, jj, kk)), temp, ttm)
                        actxmzm = actsite * ttm%act_ele_cell(ixmizm, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk)), temp, ttm)
                        actypzp = actsite * ttm%act_ele_cell(iypizp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iypizp, ii, jj, kk)), temp, ttm)
                        actypzm = actsite * ttm%act_ele_cell(iypizm, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iypizm, ii, jj, kk)), temp, ttm)
                        actymzp = actsite * ttm%act_ele_cell(iymizp, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iymizp, ii, jj, kk)), temp, ttm)
                        actymzm = actsite * ttm%act_ele_cell(iymizm, ii, jj, kk) * &
                                          KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk)), temp, ttm)
                        d2tdu2 = fomAx * (actxm * (ttm%eltemp(ixm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                          actxp * (ttm%eltemp(ixp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdv2 = fomAy * (actym * (ttm%eltemp(iym, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                          actyp * (ttm%eltemp(iyp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdw2 = fomAz * (actzm * (ttm%eltemp(izm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                          actzp * (ttm%eltemp(izp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdudv = fomAxy * (actxpyp * ttm%eltemp(ixpiyp, ii, jj, kk) - actxpym * ttm%eltemp(ixpiym, ii, jj, kk) - &
                                            actxmyp * ttm%eltemp(ixmiyp, ii, jj, kk) + actxmym * ttm%eltemp(ixmiym, ii, jj, kk))
                        d2tdudw = fomAxz * (actxpzp * ttm%eltemp(ixpizp, ii, jj, kk) - actxpzm * ttm%eltemp(ixpizm, ii, jj, kk) - &
                                            actxmzp * ttm%eltemp(ixmizp, ii, jj, kk) + actxmzm * ttm%eltemp(ixmizm, ii, jj, kk))
                        d2tdvdw = fomAyz * (actypzp * ttm%eltemp(iypizp, ii, jj, kk) - actypzm * ttm%eltemp(iypizm, ii, jj, kk) - &
                                            actymzp * ttm%eltemp(iymizp, ii, jj, kk) + actymzm * ttm%eltemp(iymizm, ii, jj, kk))
                        eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                               + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                               + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                      End Do
                    End Do
                  End Do
                Else
                  ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                  Do k = 1, ttm%ntcell(3)
                    Do j = 1, ttm%ntcell(2)
                      Do i = 1, ttm%ntcell(1)
                        ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                        ixm = ijk - 1
                        ixp = ijk + 1
                        iym = ijk - (ttm%ntcell(1) + 2)
                        iyp = ijk + (ttm%ntcell(1) + 2)
                        izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        ixpiyp = ijk + ttm%ntcell(1) + 3
                        ixpiym = ijk - ttm%ntcell(1) - 1
                        ixmiyp = ijk + ttm%ntcell(1) + 1
                        ixmiym = ijk - ttm%ntcell(1) - 3
                        ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        actsite = 1.0_wp / Ce(ttm%eltemp(ijk, ii, jj, kk), ttm)
                        actxm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixm, ii, jj, kk)), temp, ttm)
                        actxp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk)), temp, ttm)
                        actym = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iym, ii, jj, kk)), temp, ttm)
                        actyp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk)), temp, ttm)
                        actzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(izm, ii, jj, kk)), temp, ttm)
                        actzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk)), temp, ttm)
                        actxpyp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpiyp, ii, jj, kk)), temp, ttm)
                        actxpym = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpiym, ii, jj, kk)), temp, ttm)
                        actxmyp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmiyp, ii, jj, kk)), temp, ttm)
                        actxmym = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk)), temp, ttm)
                        actxpzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpizp, ii, jj, kk)), temp, ttm)
                        actxpzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpizm, ii, jj, kk)), temp, ttm)
                        actxmzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmizp, ii, jj, kk)), temp, ttm)
                        actxmzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk)), temp, ttm)
                        actypzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iypizp, ii, jj, kk)), temp, ttm)
                        actypzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iypizm, ii, jj, kk)), temp, ttm)
                        actymzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iymizp, ii, jj, kk)), temp, ttm)
                        actymzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk)), temp, ttm)
                        d2tdu2 = fomAx * (actxm * (ttm%eltemp(ixm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                          actxp * (ttm%eltemp(ixp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdv2 = fomAy * (actym * (ttm%eltemp(iym, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                          actyp * (ttm%eltemp(iyp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdw2 = fomAz * (actzm * (ttm%eltemp(izm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                          actzp * (ttm%eltemp(izp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdudv = fomAxy * (actxpyp * ttm%eltemp(ixpiyp, ii, jj, kk) - actxpym * ttm%eltemp(ixpiym, ii, jj, kk) - &
                                            actxmyp * ttm%eltemp(ixmiyp, ii, jj, kk) + actxmym * ttm%eltemp(ixmiym, ii, jj, kk))
                        d2tdudw = fomAxz * (actxpzp * ttm%eltemp(ixpizp, ii, jj, kk) - actxpzm * ttm%eltemp(ixpizm, ii, jj, kk) - &
                                            actxmzp * ttm%eltemp(ixmizp, ii, jj, kk) + actxmzm * ttm%eltemp(ixmizm, ii, jj, kk))
                        d2tdvdw = fomAyz * (actypzp * ttm%eltemp(iypizp, ii, jj, kk) - actypzm * ttm%eltemp(iypizm, ii, jj, kk) - &
                                            actymzp * ttm%eltemp(iymizp, ii, jj, kk) + actymzm * ttm%eltemp(iymizm, ii, jj, kk))
                        eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                               + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                               + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                      End Do
                    End Do
                  End Do
                End If

              End Do
            End Do
          End Do

        Else
          ! standard thermal diffusion calculation applies when energy redistribution is not applicable
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)
                Do k = 1, ttm%ntcell(3)
                  Do j = 1, ttm%ntcell(2)
                    Do i = 1, ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                      ixm = ijk - 1
                      ixp = ijk + 1
                      iym = ijk - (ttm%ntcell(1) + 2)
                      iyp = ijk + (ttm%ntcell(1) + 2)
                      izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                      izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                      ixpiyp = ijk + ttm%ntcell(1) + 3
                      ixpiym = ijk - ttm%ntcell(1) - 1
                      ixmiyp = ijk + ttm%ntcell(1) + 1
                      ixmiym = ijk - ttm%ntcell(1) - 3
                      ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                      ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                      ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                      ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                      iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                      iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                      iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                      iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                      actsite = 1.0_wp / Ce(ttm%eltemp(ijk, ii, jj, kk), ttm)
                      actxm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixm, ii, jj, kk)), temp, ttm)
                      actxp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk)), temp, ttm)
                      actym = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iym, ii, jj, kk)), temp, ttm)
                      actyp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk)), temp, ttm)
                      actzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(izm, ii, jj, kk)), temp, ttm)
                      actzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk)), temp, ttm)
                      actxpyp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpiyp, ii, jj, kk)), temp, ttm)
                      actxpym = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpiym, ii, jj, kk)), temp, ttm)
                      actxmyp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmiyp, ii, jj, kk)), temp, ttm)
                      actxmym = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk)), temp, ttm)
                      actxpzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpizp, ii, jj, kk)), temp, ttm)
                      actxpzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixpizm, ii, jj, kk)), temp, ttm)
                      actxmzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmizp, ii, jj, kk)), temp, ttm)
                      actxmzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk)), temp, ttm)
                      actypzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iypizp, ii, jj, kk)), temp, ttm)
                      actypzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iypizm, ii, jj, kk)), temp, ttm)
                      actymzp = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iymizp, ii, jj, kk)), temp, ttm)
                      actymzm = actsite * KeD(0.5_wp * (ttm%eltemp(ijk, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk)), temp, ttm)
                      d2tdu2 = fomAx * (actxm * (ttm%eltemp(ixm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                        actxp * (ttm%eltemp(ixp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                      d2tdv2 = fomAy * (actym * (ttm%eltemp(iym, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                        actyp * (ttm%eltemp(iyp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                      d2tdw2 = fomAz * (actzm * (ttm%eltemp(izm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                        actzp * (ttm%eltemp(izp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                      d2tdudv = fomAxy * (actxpyp * ttm%eltemp(ixpiyp, ii, jj, kk) - actxpym * ttm%eltemp(ixpiym, ii, jj, kk) - &
                                          actxmyp * ttm%eltemp(ixmiyp, ii, jj, kk) + actxmym * ttm%eltemp(ixmiym, ii, jj, kk))
                      d2tdudw = fomAxz * (actxpzp * ttm%eltemp(ixpizp, ii, jj, kk) - actxpzm * ttm%eltemp(ixpizm, ii, jj, kk) - &
                                          actxmzp * ttm%eltemp(ixmizp, ii, jj, kk) + actxmzm * ttm%eltemp(ixmizm, ii, jj, kk))
                      d2tdvdw = fomAyz * (actypzp * ttm%eltemp(iypizp, ii, jj, kk) - actypzm * ttm%eltemp(iypizm, ii, jj, kk) - &
                                          actymzp * ttm%eltemp(iymizp, ii, jj, kk) + actymzm * ttm%eltemp(iymizm, ii, jj, kk))
                      eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                             + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                             + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        End If

      Case (TTM_KE_TABULATED)
        ! tabulated thermal conductivity: uses local ionic or system temperature to calculate value
        If (ttm%redistribute) Then
          ! system with cell deactivation/energy redistribution
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)

                If (ii > -2 .and. ii < 2 .and. jj > -2 .and. jj < 2 .and. kk > -2 .and. kk < 2) Then
                  ! replace electronic temperatures with values required for energy redistribution
                  Do ijk = 1, ttm%numcell
                    If (ttm%adjust(ijk, ii, jj, kk)) ttm%eltemp(ijk, ii, jj, kk) = ttm%eltemp_adj(ijk, ii, jj, kk)
                  End Do
                  ! calculate thermal diffusion only for active ionic temperature sites (and active neighbours)
                  Do k = 1, ttm%ntcell(3)
                    Do j = 1, ttm%ntcell(2)
                      Do i = 1, ttm%ntcell(1)
                        ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                        ixm = ijk - 1
                        ixp = ijk + 1
                        iym = ijk - (ttm%ntcell(1) + 2)
                        iyp = ijk + (ttm%ntcell(1) + 2)
                        izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        ixpiyp = ijk + ttm%ntcell(1) + 3
                        ixpiym = ijk - ttm%ntcell(1) - 1
                        ixmiyp = ijk + ttm%ntcell(1) + 1
                        ixmiym = ijk - ttm%ntcell(1) - 3
                        ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        actsite = ttm%act_ele_cell(ijk, ii, jj, kk)
                        actxm = actsite * ttm%act_ele_cell(ixm, ii, jj, kk)
                        actxp = actsite * ttm%act_ele_cell(ixp, ii, jj, kk)
                        actym = actsite * ttm%act_ele_cell(iym, ii, jj, kk)
                        actyp = actsite * ttm%act_ele_cell(iyp, ii, jj, kk)
                        actzm = actsite * ttm%act_ele_cell(izm, ii, jj, kk)
                        actzp = actsite * ttm%act_ele_cell(izp, ii, jj, kk)
                        actxpyp = actsite * ttm%act_ele_cell(ixpiyp, ii, jj, kk)
                        actxpym = actsite * ttm%act_ele_cell(ixpiym, ii, jj, kk)
                        actxmyp = actsite * ttm%act_ele_cell(ixmiyp, ii, jj, kk)
                        actxmym = actsite * ttm%act_ele_cell(ixmiym, ii, jj, kk)
                        actxpzp = actsite * ttm%act_ele_cell(ixpizp, ii, jj, kk)
                        actxpzm = actsite * ttm%act_ele_cell(ixpizm, ii, jj, kk)
                        actxmzp = actsite * ttm%act_ele_cell(ixmizp, ii, jj, kk)
                        actxmzm = actsite * ttm%act_ele_cell(ixmizm, ii, jj, kk)
                        actypzp = actsite * ttm%act_ele_cell(iypizp, ii, jj, kk)
                        actypzm = actsite * ttm%act_ele_cell(iypizm, ii, jj, kk)
                        actymzp = actsite * ttm%act_ele_cell(iymizp, ii, jj, kk)
                        actymzm = actsite * ttm%act_ele_cell(iymizm, ii, jj, kk)
                        eltempKe = Merge(ttm%tempion(ijk), temp, (ii == 0 .and. jj == 0 .and. kk == 0))
                        alploc = Ke(eltempKe, ttm) / Ce(ttm%eltemp(ijk, ii, jj, kk), ttm)
                        d2tdu2 = fomAx * alploc * (actxm * (ttm%eltemp(ixm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                                   actxp * (ttm%eltemp(ixp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdv2 = fomAy * alploc * (actym * (ttm%eltemp(iym, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                                   actyp * (ttm%eltemp(iyp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdw2 = fomAz * alploc * (actzm * (ttm%eltemp(izm, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)) + &
                                                   actzp * (ttm%eltemp(izp, ii, jj, kk) - ttm%eltemp(ijk, ii, jj, kk)))
                        d2tdudv = fomAxy * alploc * (actxpyp * ttm%eltemp(ixpiyp, ii, jj, kk) - &
                                                     actxpym * ttm%eltemp(ixpiym, ii, jj, kk) - &
                                                     actxmyp * ttm%eltemp(ixmiyp, ii, jj, kk) + &
                                                     actxmym * ttm%eltemp(ixmiym, ii, jj, kk))
                        d2tdudw = fomAxz * alploc * (actxpzp * ttm%eltemp(ixpizp, ii, jj, kk) - &
                                                     actxpzm * ttm%eltemp(ixpizm, ii, jj, kk) - &
                                                     actxmzp * ttm%eltemp(ixmizp, ii, jj, kk) + &
                                                     actxmzm * ttm%eltemp(ixmizm, ii, jj, kk))
                        d2tdvdw = fomAyz * alploc * (actypzp * ttm%eltemp(iypizp, ii, jj, kk) - &
                                                     actypzm * ttm%eltemp(iypizm, ii, jj, kk) - &
                                                     actymzp * ttm%eltemp(iymizp, ii, jj, kk) + &
                                                     actymzm * ttm%eltemp(iymizm, ii, jj, kk))
                        eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                               + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                               + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                      End Do
                    End Do
                  End Do
                Else
                  ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                  Do k = 1, ttm%ntcell(3)
                    Do j = 1, ttm%ntcell(2)
                      Do i = 1, ttm%ntcell(1)
                        ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                        ixm = ijk - 1
                        ixp = ijk + 1
                        iym = ijk - (ttm%ntcell(1) + 2)
                        iyp = ijk + (ttm%ntcell(1) + 2)
                        izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                        ixpiyp = ijk + ttm%ntcell(1) + 3
                        ixpiym = ijk - ttm%ntcell(1) - 1
                        ixmiyp = ijk + ttm%ntcell(1) + 1
                        ixmiym = ijk - ttm%ntcell(1) - 3
                        ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                        ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                        iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                        iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                        ! note that temperature for thermal conductivity is always system
                        ! temperature for electronic cells away from ionic cells
                        alploc = Ke(temp, ttm) / Ce(ttm%eltemp(ijk, ii, jj, kk), ttm)
                        d2tdu2 = fomAx * alploc * (ttm%eltemp(ixm, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk) - &
                                                                        2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                        d2tdv2 = fomAy * alploc * (ttm%eltemp(iym, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk) - &
                                                                        2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                        d2tdw2 = fomAz * alploc * (ttm%eltemp(izm, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk) - &
                                                                        2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                        d2tdudv = fomAxy * alploc * (ttm%eltemp(ixpiyp, ii, jj, kk) - ttm%eltemp(ixpiym, ii, jj, kk) - &
                                                     ttm%eltemp(ixmiyp, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk))
                        d2tdudw = fomAxz * alploc * (ttm%eltemp(ixpizp, ii, jj, kk) - ttm%eltemp(ixpizm, ii, jj, kk) - &
                                                     ttm%eltemp(ixmizp, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk))
                        d2tdvdw = fomAyz * alploc * (ttm%eltemp(iypizp, ii, jj, kk) - ttm%eltemp(iypizm, ii, jj, kk) - &
                                                     ttm%eltemp(iymizp, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk))
                        eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                               + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                               + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                      End Do
                    End Do
                  End Do
                End If

              End Do
            End Do
          End Do
        Else
          ! standard thermal diffusion calculation applies when energy redistribution is not applicable
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)
                Do k = 1, ttm%ntcell(3)
                  Do j = 1, ttm%ntcell(2)
                    Do i = 1, ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                      ixm = ijk - 1
                      ixp = ijk + 1
                      iym = ijk - (ttm%ntcell(1) + 2)
                      iyp = ijk + (ttm%ntcell(1) + 2)
                      izm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                      izp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
                      ixpiyp = ijk + ttm%ntcell(1) + 3
                      ixpiym = ijk - ttm%ntcell(1) - 1
                      ixmiyp = ijk + ttm%ntcell(1) + 1
                      ixmiym = ijk - ttm%ntcell(1) - 3
                      ixpizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                      ixpizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) + 1
                      ixmizp = ijk + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                      ixmizm = ijk - (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) - 1
                      iypizp = ijk + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                      iypizm = ijk + (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                      iymizp = ijk - (ttm%ntcell(1) + 2) * (1 - (ttm%ntcell(2) + 2))
                      iymizm = ijk - (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
                      eltempKe = Merge(ttm%tempion(ijk), temp, (ii == 0 .and. jj == 0 .and. kk == 0))
                      alploc = Ke(eltempKe, ttm) / Ce(ttm%eltemp(ijk, ii, jj, kk), ttm)
                      d2tdu2 = fomAx * alploc * (ttm%eltemp(ixm, ii, jj, kk) + ttm%eltemp(ixp, ii, jj, kk) - &
                                                                      2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                      d2tdv2 = fomAy * alploc * (ttm%eltemp(iym, ii, jj, kk) + ttm%eltemp(iyp, ii, jj, kk) - &
                                                                      2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                      d2tdw2 = fomAz * alploc * (ttm%eltemp(izm, ii, jj, kk) + ttm%eltemp(izp, ii, jj, kk) - &
                                                                      2.0_wp * ttm%eltemp(ijk, ii, jj, kk))
                      d2tdudv = fomAxy * alploc * (ttm%eltemp(ixpiyp, ii, jj, kk) - ttm%eltemp(ixpiym, ii, jj, kk) - &
                                                   ttm%eltemp(ixmiyp, ii, jj, kk) + ttm%eltemp(ixmiym, ii, jj, kk))
                      d2tdudw = fomAxz * alploc * (ttm%eltemp(ixpizp, ii, jj, kk) - ttm%eltemp(ixpizm, ii, jj, kk) - &
                                                   ttm%eltemp(ixmizp, ii, jj, kk) + ttm%eltemp(ixmizm, ii, jj, kk))
                      d2tdvdw = fomAyz * alploc * (ttm%eltemp(iypizp, ii, jj, kk) - ttm%eltemp(iypizm, ii, jj, kk) - &
                                                   ttm%eltemp(iymizp, ii, jj, kk) + ttm%eltemp(iymizm, ii, jj, kk))
                      eltemp1(ijk, ii, jj, kk) = ttm%eltemp(ijk, ii, jj, kk) + ttm%tdiffw(1) * d2tdu2  + ttm%tdiffw(2) * d2tdv2 &
                                                                             + ttm%tdiffw(3) * d2tdw2  + ttm%tdiffw(4) * d2tdudv &
                                                                             + ttm%tdiffw(5) * d2tdudw + ttm%tdiffw(6) * d2tdvdw
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        End If

      End Select

      ! electron stopping and electron-phonon couplings

      If (ttm%oneway) Then
        If (nstep > ttm%nstepcpl) Then
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              Do i = 1, ttm%ntcell(1)
                ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                If (ttm%act_ele_cell(ijk, 0, 0, 0) > zero_plus) Then
                  alploc = tstep * ttm%rvolume / (Ce(ttm%eltemp(ijk, 0, 0, 0), ttm) * Real(redtstepmx, Kind=wp))
                  ! e-s coupling term
                  eltemp1(ijk, 0, 0, 0) = eltemp1(ijk, 0, 0, 0) + alploc * ttm%asource(ijk)
                  ! e-p coupling term: only use if electronic temperature
                  ! exceeds ionic temperature
                  If (ttm%l_epcp .and. ttm%eltemp(ijk, 0, 0, 0) > ttm%tempion(ijk)) Then
                    Select Case (ttm%gvar)
                    Case (TTM_EPVAR_NULL, TTM_EPVAR_HOMO)
                      eltemp1(ijk, 0, 0, 0) = eltemp1(ijk, 0, 0, 0) - alploc * ttm%gsource(ijk) * &
                                                                      (ttm%eltemp(ijk, 0, 0, 0) - ttm%tempion(ijk))
                    Case (TTM_EPVAR_HETERO)
                      eltemp1(ijk, 0, 0, 0) = eltemp1(ijk, 0, 0, 0) - alploc * ttm%gsource(ijk) * &
                                 Gep(ttm%eltemp(ijk, 0, 0, 0), ttm) * (ttm%eltemp(ijk, 0, 0, 0) - ttm%tempion(ijk))
                    End Select
                  End If
                End If
              End Do
            End Do
          End Do
        End If
      Else
        If (nstep > ttm%nstepcpl) Then
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              Do i = 1, ttm%ntcell(1)
                ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                If (ttm%act_ele_cell(ijk, 0, 0, 0) > zero_plus) Then
                  alploc = tstep * ttm%rvolume / (Ce(ttm%eltemp(ijk, 0, 0, 0), ttm) * Real(redtstepmx, Kind=wp))
                  ! e-s coupling term
                  eltemp1(ijk, 0, 0, 0) = eltemp1(ijk, 0, 0, 0) + alploc * ttm%asource(ijk)
                  ! e-p coupling term
                  If (ttm%l_epcp) Then
                    Select Case (ttm%gvar)
                    Case (TTM_EPVAR_NULL, TTM_EPVAR_HOMO)
                      eltemp1(ijk, 0, 0, 0) = eltemp1(ijk, 0, 0, 0) - alploc * ttm%gsource(ijk) * &
                                                                      (ttm%eltemp(ijk, 0, 0, 0) - ttm%tempion(ijk))
                    Case (TTM_EPVAR_HETERO)
                      eltemp1(ijk, 0, 0, 0) = eltemp1(ijk, 0, 0, 0) - alploc * ttm%gsource(ijk) * &
                                 Gep(ttm%eltemp(ijk, 0, 0, 0), ttm) * (ttm%eltemp(ijk, 0, 0, 0) - ttm%tempion(ijk))
                    End Select
                  End If
                End If
              End Do
            End Do
          End Do
        End If
      End If

      ! update electronic temperatures to ttm%adjusted values

      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do ii = -ttm%eltcell(1), ttm%eltcell(1)
            Do ijk = 1, ttm%numcell
              ttm%eltemp(ijk, ii, jj, kk) = eltemp1(ijk, ii, jj, kk)
            End Do
          End Do
        End Do
      End Do

      ! update boundary halo values and apply boundary conditions

      Call boundaryHalo(ttm, domain, comm)
      Call boundaryCond(ttm%bcTypeE, temp, ttm, comm)

      ! simple stability check for simulation

      If (Any(ttm%eltemp < 0.0_wp)) safe = .false.
      Call gcheck(comm, safe)
      If (.not. safe) Call error(693)

    End Do

    ! Dumping Te file every ndump steps
    Call ttm_system_revive('DUMP_E', nstep, time, ndump, nstrun, ttm, comm)

    Deallocate (eltemp1, Stat=fail)
    If (fail > 0) Call error(1088)
  End Subroutine ttm_thermal_diffusion
End Module ttm_track
