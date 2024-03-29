Module ttm_utils

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for utility functions for the two-temperature model
  ! (ttm)
  !
  ! copyright - daresbury laboratory
  ! authors   - s.l.daraszewicz & m.a.seaton july 2012
  ! contrib   - g.khara may 2016
  ! contrib   - m.a.seaton february 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: Grid1_tag,&
                             Grid2_tag,&
                             comms_type,&
                             gmax,&
                             gmin,&
                             gsum
  Use constants,       Only: zero_plus,&
                             kB_to_eV
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error
  Use kinds,           Only: wp
  Use ttm,             Only: eltemp_max,&
                             eltemp_mean,&
                             eltemp_min,&
                             eltemp_sum,&
                             ttm_type,&
                             TTM_CE_CONST, TTM_CE_TANH, TTM_CE_LINEAR, TTM_CE_TABULATED,&
                             TTM_CE_CONST_DYN, TTM_CE_TANH_DYN, TTM_CE_LINEAR_DYN, TTM_CE_TABULATED_DYN,&
                             TTM_KE_TABULATED, TTM_DE_METAL, TTM_DE_CONST, TTM_DE_RECIP, TTM_DE_TABULATED

#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None

  Private

  Public ::  Ce, Ke, alp, ked, peakProfilerElec, &
            peakProfiler, printElecLatticeStatsToFile, gep, calcchies, &
            printLatticeStatsToFile, redistribute_Te
Contains

  ! indeed this is not a joke...
  Function Kep(T, ttm)

    ! Calcuate interpolated thermal conductivity from tabulated values: given as kB/A^3

    Real(Kind=wp),  Intent(In   ) :: T
    Type(ttm_type), Intent(InOut) :: ttm
    Real(Kind=wp)                 :: Kep

    Call interpolate(ttm%kel, ttm%ketable, T, Kep)

    Return

  End Function Kep

  Function Gep(T, ttm)

    ! Calculate electron-phonon coupling term (friction parameter for inhomogeneous
    ! Langevin thermostat): given as ps^-1

    Real(Kind=wp),  Intent(In   ) :: T
    Type(ttm_type), Intent(InOut) :: ttm
    Real(Kind=wp)                 :: Gep

    Call interpolate(ttm%gel, ttm%gtable, T, Gep)
    Gep = Gep * Merge(ttm%rcellrho, 1.0_wp, ttm%ttmdyndens)

    Return

  End Function Gep

  Function Ce(T, ttm)

    ! Temperature-dependent specific heat capacity (different for metals and insulators)
    ! given as kB/A^3: conversions from kB/atom carried out for cases 0, 1 and 2, allowing
    ! for changes in number of atoms per voxel (ttm%cellrho) if dynamic cell density option
    ! selected

    Real(Kind=wp),  Intent(In   ) :: T
    Type(ttm_type), Intent(InOut) :: ttm
    Real(Kind=wp)                 :: Ce

    Ce = 1.0_wp

    Select Case (ttm%CeType)
    Case (TTM_CE_CONST)
      ! Case 0: constant specific heat capacity (given as kB/A^3)
      Ce = ttm%Ce0
    Case (TTM_CE_TANH)
      ! Case 1: hyperbolic tangent specific heat capacity (given as kB/A^3) -
      !         Ce = ttm%sh_A*Tanh(T*ttm%sh_B*1.0e-4) [kB/atom]
      Ce = ttm%sh_A * Tanh(T * ttm%sh_B)
    Case (TTM_CE_LINEAR)
      ! Case 2: linear specific heat capacity to maximum value at/beyond ttm%Tfermi
      !         (given as kB/A^3)
      Ce = Min(T / ttm%Tfermi, 1.0_wp) * ttm%Cemax
    Case (TTM_CE_CONST_DYN)
      ! Case 4: constant specific heat capacity (converted as kB/A^3)
      Ce = ttm%Ce0 * ttm%cellrho
    Case (TTM_CE_TANH_DYN)
      ! Case 5: hyperbolic tangent specific heat capacity (converted to kB/A^3) -
      !         Ce = ttm%sh_A*Tanh(T*ttm%sh_B*1.0e-4) [kB/atom]
      Ce = ttm%sh_A * ttm%cellrho * Tanh(T * ttm%sh_B)
    Case (TTM_CE_LINEAR_DYN)
      ! Case 6: linear specific heat capacity to maximum value at/beyond ttm%Tfermi
      !         (converted to kB/A^3)
      Ce = Min(T / ttm%Tfermi, 1.0_wp) * ttm%Cemax * ttm%cellrho
    Case (TTM_CE_TABULATED, TTM_CE_TABULATED_DYN)
      ! Case 3: interpolated specific heat capacity from table (given as kB/A^3)
      Call interpolate(ttm%cel, ttm%cetable, T, Ce)
    Case Default
      ! Default case: constant specific heat capacity of 1 kB/atom (convert to kB/A^3)
      Ce = ttm%cellrho
    End Select

  End Function Ce

  Function Ke(T, ttm)

    ! Temperature-dependent lattice thermal conductivity: given as kB/(ps A)

    Real(Kind=wp),  Intent(In   ) :: T
    Type(ttm_type), Intent(InOut) :: ttm
    Real(Kind=wp)                 :: Ke

    Select Case (ttm%KeType)
    Case (TTM_KE_TABULATED)
      ! Case 3: thermal conductivity interpolated from table
      Ke = Kep(T, ttm)
    Case Default
      ! Case 1: constant thermal conductivity
      Ke = ttm%Ka0
    End Select

  End Function Ke

  Function KeD(Te, temp, ttm)

    ! Temperature-dependent Drude-like lattice thermal conductivity: given as kB/(ps A)
    ! (ttm%Ka0 = thermal conductivity at system temperature)

    Real(Kind=wp),  Intent(In   ) :: Te, temp
    Type(ttm_type), Intent(InOut) :: ttm
    Real(Kind=wp)                 :: KeD

    KeD = ttm%Ka0 * Ce(Te, ttm) / Ce(temp, ttm)

  End Function KeD

  Function alp(Te, ttm)

    ! Thermal diffusivity: given as A^2/ps
    ! (different cases for metals and insulators)

    Real(Kind=wp),  Intent(In   ) :: Te
    Type(ttm_type), Intent(InOut) :: ttm
    Real(Kind=wp)                 :: alp

    Select Case (ttm%DeType)
    Case (TTM_DE_METAL)
      ! Case 0: metal system - ratio of conductivity/heat capacity
      alp = Ke(Te, ttm) / Ce(Te, ttm)
    Case (TTM_DE_CONST)
      ! Case 1: non-metal system - constant value (given as A^2/ps)
      alp = ttm%Diff0
    Case (TTM_DE_RECIP)
      ! Case 2: non-metal system - reciprocal function of temperature
      !         up to Fermi temperature, ttm%Diff0 previously scaled with
      !         system temperature (given as A^2/ps)
      alp = ttm%Diff0 / Min(Te, ttm%Tfermi)
    Case (TTM_DE_TABULATED)
      ! Case 3: non-metal system - thermal diffusivity interpolated
      !         from table (given as A^2/ps)
      Call interpolate(ttm%del, ttm%detable, Te, alp)
    Case Default
      alp = 0.0_wp
    End Select

  End Function alp

  Subroutine calcchies(chi, ttm, comm)

    ! Calculate electron-phonon coupling friction term (chi:
    ! for homogeneously coupled system: uses mean electronic temperature
    ! and interpolates from tabulated values given in g.dat file

    Real(Kind=wp),    Intent(  Out) :: chi
    Type(ttm_type),   Intent(InOut) :: ttm
    Type(comms_type), Intent(InOut) :: comm

    Real(Kind=wp) :: eltempav

    Call eltemp_mean(eltempav, ttm, comm)
    Call interpolate(ttm%gel, ttm%gtable, eltempav, chi)

    if (ttm%ttmdyndens) then
      chi = chi * ttm%rcellrho
    end if

  End Subroutine calcchies

  Subroutine peakProfiler(lat, peakfile, nstep, freq, ttm, comm)

    ! prints a slice across y (centred in xz-plane) of a lattice to a file

    Type(ttm_type),                        Intent(InOut) :: ttm
    Integer,                               Intent(In   ) :: freq, nstep
    Character(Len=*),                      Intent(In   ) :: peakfile
    Real(Kind=wp), Dimension(ttm%numcell), Intent(In   ) :: lat
    Type(comms_type),                      Intent(InOut) :: comm

    Integer                                :: i, ijk, iounit = 114, j, k, l
    Real(Kind=wp), Dimension(ttm%ntsys(2)) :: laty

    laty = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep, freq) == 0) Then
        i = ttm%midI(1) - ttm%ntcelloff(1)
        k = ttm%midI(3) - ttm%ntcelloff(3)
        If (i > 0 .and. i <= ttm%ntcell(1) .and. k > 0 .and. k <= ttm%ntcell(3)) Then
          Do j = 1, ttm%ntcell(2)
            l = j + ttm%ntcelloff(2)
            ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + k * (ttm%ntcell(2) + 2))
            laty(l) = lat(ijk)
          End Do
        End If
        Call gsum(comm, laty)
        If (comm%idnode == 0) Then
          If (nstep == 0) Then
            Open (Unit=iounit, File=peakfile, Status='replace')
          Else
            Open (Unit=iounit, File=peakfile, Position='append')
          End If
          Do i = 1, ttm%ntsys(2)
            Write (iounit, Fmt='(9es12.4,1p)', Advance='no') laty(i)
          End Do
          Close (iounit)
        End If
      End If
    End If

  End Subroutine peakProfiler

  Subroutine peakProfilerElec(peakfile, nstep, freq, ttm, comm)

    ! prints a slice across y (centred in xz-plane) of the electronic temperature lattice to a file

    Character(Len=*), Intent(In   ) :: peakfile
    Integer,          Intent(In   ) :: nstep, freq
    Type(ttm_type),   Intent(InOut) :: ttm
    Type(comms_type), Intent(InOut) :: comm

    Integer                                 :: i, ijk, iounit = 114, j, jj, jmax, jmin, k, l
    Real(Kind=wp), Dimension(ttm%eltsys(2)) :: laty

    laty = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep, freq) == 0) Then
        i = ttm%midI(1) - ttm%ntcelloff(1)
        k = ttm%midI(3) - ttm%ntcelloff(3)
        If (i > 0 .and. i <= ttm%ntcell(1) .and. k > 0 .and. k <= ttm%ntcell(3)) Then
          Do jj = -ttm%eltcell(2), ttm%eltcell(2)
            If (ttm%eltcell(2) > 0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3) >= 0) Then
              jmin = ttm%ttmbc(3)
            Else
              jmin = 1
            End If
            If (ttm%eltcell(2) > 0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4) >= 0) Then
              jmax = ttm%ttmbc(4)
            Else
              jmax = ttm%ntcell(2)
            End If
            Do j = jmin, jmax
              l = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
              ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + k * (ttm%ntcell(2) + 2))
              If (l > 0 .and. l <= ttm%eltsys(2)) laty(l) = ttm%eltemp(ijk, 0, jj, 0)
            End Do
          End Do
        End If
        Call gsum(comm, laty)
        If (comm%idnode == 0) Then
          If (nstep == 0) Then
            Open (Unit=iounit, File=peakfile, Status='replace')
          Else
            Open (Unit=iounit, File=peakfile, Position='append')
          End If
          Do i = 1, ttm%eltsys(2)
            Write (iounit, Fmt='(9es12.4,1p)', Advance='no') laty(i)
          End Do
          Close (iounit)
        End If
      End If
    End If

  End Subroutine peakProfilerElec

  Subroutine printLatticeStatsToFile(lat, latfile, time, nstep, freq, ttm, comm)

    ! prints lattice statistics (minimum, maximum, sum) to file

    Type(ttm_type),                        Intent(InOut) :: ttm
    Integer,                               Intent(In   ) :: freq, nstep
    Real(Kind=wp),                         Intent(In   ) :: time
    Character(Len=*),                      Intent(In   ) :: latfile
    Real(Kind=wp), Dimension(ttm%numcell), Intent(In   ) :: lat
    Type(comms_type),                      Intent(InOut) :: comm

    Integer       :: i, ijk, iounit = 113, j, k
    Real(kind=wp) :: lat_max, lat_min, lat_sum, rtotal

    lat_sum = 0.0_wp
    lat_min = 1.0e30_wp
    lat_max = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep, freq) == 0) Then

        ! check over all active ionic temperature config%cells (CIT)
        Do k = 1, ttm%ntcell(3)
          Do j = 1, ttm%ntcell(2)
            Do i = 1, ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
              If (ttm%act_ele_cell(ijk, 0, 0, 0) > zero_plus) Then
                lat_sum = lat_sum + lat(ijk)
                lat_min = Min(lat_min, lat(ijk))
                lat_max = Max(lat_max, lat(ijk))
              End If
            End Do
          End Do
        End Do

        Call gsum(comm, lat_sum)
        Call gmin(comm, lat_min)
        Call gmax(comm, lat_max)

        rtotal = Merge(1.0_wp / Real(ttm%acell, Kind=wp), 0.0_wp, (ttm%acell > 0))

        If (comm%idnode == 0) Then
          If (nstep == 0) Then
            Open (Unit=iounit, File=latfile, Status='replace')
          Else
            Open (Unit=iounit, File=latfile, Position='append')
          End If
          Write (iounit, '(i8,5(1p,es12.4))') nstep, time, lat_min, lat_max, lat_sum * rtotal, lat_sum
          Close (iounit)
        End If

      End If
    End If

  End Subroutine printLatticeStatsToFile

  Subroutine printElecLatticeStatsToFile(latfile, time, temp0, nstep, freq, ttm, comm)

    ! prints electronic temperature lattice statistics (minimum, maximum, sum)
    ! and energy (E = integral of Ce(Te,ttm)*Te between temp0 and Te) to file

    Character(Len=*), Intent(In   ) :: latfile
    Real(Kind=wp),    Intent(In   ) :: time, temp0
    Integer,          Intent(In   ) :: nstep, freq
    Type(ttm_type),   Intent(InOut) :: ttm
    Type(comms_type), Intent(InOut) :: comm

    Integer       :: i, ii, ijk, imax, imin, iounit = 115, j, jj, jmax, jmin, k, kk, kmax, kmin, &
                     lx, ly, lz, n, numint
    Real(Kind=wp) :: Ce0a, Cemaxa, eltmp, lat_max, lat_min, lat_sum, rtotal, sgnplus, sh_Aa, tmp, &
                     totalcell, Ue

    Ce0a = ttm%Ce0 * Merge(ttm%cellrho, 1.0_wp, ttm%ttmdyndens)
    sh_Aa = ttm%sh_A * Merge(ttm%cellrho, 1.0_wp, ttm%ttmdyndens)
    Cemaxa = ttm%Cemax * Merge(ttm%cellrho, 1.0_wp, ttm%ttmdyndens)

    If (freq /= 0) Then
      If (Mod(nstep, freq) == 0) Then

        Call eltemp_sum(lat_sum, ttm, comm)
        Call eltemp_min(lat_min, ttm, comm)
        Call eltemp_max(lat_max, ttm, comm)

        Ue = 0.0_wp
        totalcell = 0.0_wp

        Do kk = -ttm%eltcell(3), ttm%eltcell(3)
          If (ttm%eltcell(3) > 0 .and. kk == -ttm%eltcell(3) .and. ttm%ttmbcmap(5) >= 0) Then
            kmin = ttm%ttmbc(5)
          Else
            kmin = 1
          End If

          If (ttm%eltcell(3) > 0 .and. kk == ttm%eltcell(3) .and. ttm%ttmbcmap(6) >= 0) Then
            kmax = ttm%ttmbc(6)
          Else
            kmax = ttm%ntcell(3)
          End If

          Do jj = -ttm%eltcell(2), ttm%eltcell(2)
            If (ttm%eltcell(2) > 0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3) >= 0) Then
              jmin = ttm%ttmbc(3)
            Else
              jmin = 1
            End If
            If (ttm%eltcell(2) > 0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4) >= 0) Then
              jmax = ttm%ttmbc(4)
            Else
              jmax = ttm%ntcell(2)
            End If

            Do ii = -ttm%eltcell(1), ttm%eltcell(1)
              If (ttm%eltcell(1) > 0 .and. ii == -ttm%eltcell(1) .and. ttm%ttmbcmap(1) >= 0) Then
                imin = ttm%ttmbc(1)
              Else
                imin = 1
              End If
              If (ttm%eltcell(1) > 0 .and. ii == ttm%eltcell(1) .and. ttm%ttmbcmap(2) >= 0) Then
                imax = ttm%ttmbc(2)
              Else
                imax = ttm%ntcell(1)
              End If
              Do k = kmin, kmax
                lz = k + ttm%ntcelloff(3) + (kk + ttm%eltcell(3)) * ttm%ntsys(3) - ttm%zeroE(3)
                Do j = jmin, jmax
                  ly = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
                  Do i = imin, imax
                    lx = i + ttm%ntcelloff(1) + (ii + ttm%eltcell(1)) * ttm%ntsys(1) - ttm%zeroE(1)
                    ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
                    ! integrate electronic heat capacity from temp0 to electronic temperature
                    ! of config%cell if it is within global range and active (within ionic temperature config%cells)
                    If (lx > 0 .and. lx <= ttm%eltsys(1) .and. ly > 0 .and. &
                         ly <= ttm%eltsys(2) .and. lz > 0 .and. &
                         lz <= ttm%eltsys(3)) Then
                      eltmp = ttm%eltemp(ijk, ii, jj, kk)
                      tmp = Merge(1.0_wp, 0.0_wp, &
                           (ttm%act_ele_cell(ijk, 0, 0, 0) > zero_plus .or. (ii /= 0 .or. jj /= 0 .or. kk /= 0)))
                      totalcell = totalcell + tmp

                      Select Case (ttm%CeType)
                      Case (TTM_CE_CONST, TTM_CE_CONST_DYN)
                        ! constant specific heat capacity
                        tmp = tmp * Ce0a * (eltmp - temp0)
                      Case (TTM_CE_TANH, TTM_CE_TANH_DYN)
                        ! hyperbolic tangent specific heat capacity
                        tmp = tmp * sh_Aa * Log(Cosh(ttm%sh_B * eltmp) / Cosh(ttm%sh_B * temp0)) / ttm%sh_B
                      Case (TTM_CE_LINEAR, TTM_CE_LINEAR_DYN)
                        ! linear specific heat capacity to Fermi temperature
                        tmp = tmp * Cemaxa * (0.5_wp * ((Min(ttm%Tfermi, eltmp))**2 - temp0 * temp0) / ttm%Tfermi + &
                             Max(eltmp - ttm%Tfermi, 0.0_wp))
                      Case (TTM_CE_TABULATED, TTM_CE_TABULATED_DYN)
                        ! tabulated ttm%volumetric heat capacity or more complex
                        ! functions: integrate using trapezium rule
                        ! with final interval of <=1 kelvin
                        numint = Floor(tmp * (eltmp - temp0))
                        tmp = 0.0_wp
                        sgnplus = Sign(1.0_wp, Real(numint, Kind=wp))
                        Do n = 1, Abs(numint)
                          tmp = tmp + 0.5_wp * sgnplus * (Ce(temp0 + sgnplus * Real(n - 1, Kind=wp), ttm) + &
                               Ce(temp0 + sgnplus * Real(n, Kind=wp), ttm))
                        End Do
                        tmp = tmp + 0.5_wp * (Ce(temp0 + Real(numint, Kind=wp), ttm) + Ce(eltmp, ttm))
                      Case Default
                        Call error(0, 'Unrecognised TTM heat capacity model')
                      End Select
                      Ue = Ue + tmp * ttm%volume * kB_to_eV
                    End If
                  End Do
                End Do
              End Do

            End Do
          End Do
        End Do

        Call gsum(comm, Ue)
        Call gsum(comm, totalcell)

        rtotal = Merge(1.0_wp / totalcell, 0.0_wp, (totalcell > zero_plus))

        If (comm%idnode == 0) Then
          If (nstep == 0) Then
            Open (Unit=iounit, File=latfile, Status='replace')
          Else
            Open (Unit=iounit, File=latfile, Position='append')
          End If
          Write (iounit, '(i8,6(1p,es12.4))') nstep, time, lat_min, lat_max, lat_sum * rtotal, lat_sum, Ue
          Close (iounit)
        End If

      End If
    End If

  End Subroutine printElecLatticeStatsToFile

  Subroutine interpolate(tabsize, tab, x0, resp)

    ! Interpolation helper subroutine

    Integer,                                  Intent(In   ) :: tabsize
    Real(Kind=wp), Dimension(1:tabsize, 1:2), Intent(In   ) :: tab
    Real(Kind=wp),                            Intent(In   ) :: x0
    Real(Kind=wp),                            Intent(  Out) :: resp

    Integer :: i, j, k

    If (tabsize == 1) Then
      resp = tab(1, 2)
      Return
    End If

    resp = 0.0_wp

    ! if input is within range, find nearest entry in table
    ! below this value by means of bisection - calculate
    ! output by linear interpolation between this entry
    ! and the next

    If ((x0 > tab(1, 1)) .and. (x0 < tab(tabsize, 1))) Then
      i = 1
      j = tabsize
      Do
        k = (i + j) / 2
        If (x0 < tab(k, 1)) Then
          j = k
        Else
          i = k
        End If
        If (i + 1 >= j) Then
          k = j - 1
          Exit
        End If
      End Do
      resp = tab(k, 2) + (x0 - tab(k, 1)) * (tab(k + 1, 2) - tab(k, 2)) / (tab(k + 1, 1) - tab(k, 1))

      ! if input is over maximum value, set output to value at maximum

    Else If (x0 >= tab(tabsize, 1)) Then

      resp = tab(tabsize, 2)

      ! If input is under minimum value, set output to value at minimum

    Else If (x0 <= tab(1, 1)) Then

      resp = tab(1, 2)

      ! If input is completely unusable, terminate with error

    Else

      Call error(683)

    End If

  End Subroutine interpolate

  Subroutine redistribute_Te(temp0, ttm, domain, comm)

    ! Redistribute electronic energy when electronic temperature voxels are closed

    Real(Kind=wp),      Intent(In   ) :: temp0
    Type(ttm_type),     Intent(InOut) :: ttm
    Type(domains_type), Intent(In   ) :: domain
    Type(comms_type),   Intent(In   ) :: comm

    Integer                                :: i, ierr, ii, ijk, ijk1, ijk2, ijkmx, ijkmy, ijkmz, &
                                              ijkpx, ijkpy, ijkpz, j, jj, k, kk, n, numint
    Integer, Dimension(4)                  :: req
    Integer, Dimension(MPI_STATUS_SIZE, 4) :: stat
    Real(Kind=wp)                          :: act_sur_cells, Ce0a, Cemaxa, end_Te, energy_diff, &
                                              energy_per_cell, increase, newCe, oldCe, sgnplus, &
                                              sh_Aa, start_Te, U_e
    Real(Kind=wp), Allocatable             :: buffer(:, :, :, :), energydist(:, :, :, :)

    Ce0a = ttm%Ce0 * Merge(ttm%cellrho, 1.0_wp, ttm%ttmdyndens)
    sh_Aa = ttm%sh_A * Merge(ttm%cellrho, 1.0_wp, ttm%ttmdyndens)
    Cemaxa = ttm%Cemax * Merge(ttm%cellrho, 1.0_wp, ttm%ttmdyndens)

    Allocate (energydist(ttm%numcell, -1:1, -1:1, -1:1))
    energydist = 0.0_wp
    ttm%eltemp_adj = 0.0_wp
    ttm%adjust = .false.

    Do k = 1, ttm%ntcell(3)
      Do j = 1, ttm%ntcell(2)
        Do i = 1, ttm%ntcell(1)
          ijk = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
          ijkpx = 2 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
          ijkmx = i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
          ijkpy = 1 + i + (ttm%ntcell(1) + 2) * (j + 1 + (ttm%ntcell(2) + 2) * k)
          ijkmy = 1 + i + (ttm%ntcell(1) + 2) * (j - 1 + (ttm%ntcell(2) + 2) * k)
          ijkpz = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * (k + 1))
          ijkmz = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * (k - 1))
          If (ttm%act_ele_cell(ijk, 0, 0, 0) <= zero_plus .and. ttm%old_ele_cell(ijk, 0, 0, 0) > zero_plus) Then
            ! Calculate amount of energy left in the config%cell
            ! (Ue = integral of Ce(Te,ttm)d(Te) between temp0 and Te)
            U_e = 0.0_wp
            end_Te = ttm%eltemp(ijk, 0, 0, 0)
            Select Case (ttm%CeType)
            Case (TTM_CE_CONST, TTM_CE_CONST_DYN)
              ! constant specific heat capacity
              U_e = Ce0a * (end_Te - temp0)
            Case (TTM_CE_TANH, TTM_CE_TANH_DYN)
              ! hyperbolic tangent specific heat capacity
              U_e = sh_Aa * Log(Cosh(ttm%sh_B * end_Te) / Cosh(ttm%sh_B * temp0)) / ttm%sh_B
            Case (TTM_CE_LINEAR, TTM_CE_LINEAR_DYN)
              ! linear specific heat capacity to Fermi temperature
              increase = Min(ttm%Tfermi, end_Te)
              U_e = Cemaxa * (0.5_wp * (increase * increase - temp0 * temp0) / ttm%Tfermi + Max(end_Te - ttm%Tfermi, 0.0_wp))
            Case (TTM_CE_TABULATED, TTM_CE_TABULATED_DYN)
              ! tabulated ttm%volumetric heat capacity or more complex
              ! functions: integrate using trapezium rule
              ! with final interval of <=1 kelvin
              numint = Floor(end_Te - temp0)
              sgnplus = Sign(1.0_wp, Real(numint, Kind=wp))
              Do n = 1, Abs(numint)
      U_e = U_e + 0.5_wp * sgnplus * (Ce(temp0 + sgnplus * Real(n - 1, Kind=wp), ttm) + Ce(temp0 + sgnplus * Real(n, Kind=wp), ttm))
              End Do
              U_e = U_e + 0.5_wp * (Ce(temp0 + Real(numint, Kind=wp), ttm) + Ce(end_Te, ttm))
            End Select
            ! Check how many config%cells are connected to the now turned-off config%cell
          act_sur_cells = ttm%act_ele_cell(ijkmx, 0, 0, 0) + ttm%act_ele_cell(ijkpx, 0, 0, 0) + ttm%act_ele_cell(ijkmy, 0, 0, 0) + &
                            ttm%act_ele_cell(ijkpy, 0, 0, 0) + ttm%act_ele_cell(ijkmz, 0, 0, 0) + ttm%act_ele_cell(ijkpz, 0, 0, 0)
            ! Calculate energy redistribution to each config%cell and assign to config%cells
            If (act_sur_cells > zero_plus) energy_per_cell = U_e / act_sur_cells
         If (ttm%act_ele_cell(ijkmx, 0, 0, 0) > zero_plus) energydist(ijkmx, 0, 0, 0) = energydist(ijkmx, 0, 0, 0) + energy_per_cell
         If (ttm%act_ele_cell(ijkpx, 0, 0, 0) > zero_plus) energydist(ijkpx, 0, 0, 0) = energydist(ijkpx, 0, 0, 0) + energy_per_cell
         If (ttm%act_ele_cell(ijkmy, 0, 0, 0) > zero_plus) energydist(ijkmy, 0, 0, 0) = energydist(ijkmy, 0, 0, 0) + energy_per_cell
         If (ttm%act_ele_cell(ijkpy, 0, 0, 0) > zero_plus) energydist(ijkpy, 0, 0, 0) = energydist(ijkpy, 0, 0, 0) + energy_per_cell
         If (ttm%act_ele_cell(ijkmz, 0, 0, 0) > zero_plus) energydist(ijkmz, 0, 0, 0) = energydist(ijkmz, 0, 0, 0) + energy_per_cell
         If (ttm%act_ele_cell(ijkpz, 0, 0, 0) > zero_plus) energydist(ijkpz, 0, 0, 0) = energydist(ijkpz, 0, 0, 0) + energy_per_cell
          End If
        End Do
      End Do
    End Do

    If (comm%mxnode > 1) Then

      Allocate (buffer(ttm%numcell, -1:1, -1:1, -1:1))
      buffer = 0.0_wp

      ! sum up redistributed energies (placing energies for electronic grid
      ! outside ionic grid in boundary halos: removing previous values in halos
      ! before summation)

      ! -x/+x direction
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      ijk2 = 1 + ttm%ntcell(1) + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      If (domain%idx == domain%nx - 1) Then
        ii = -1
      Else
        ii = 0
      End If
      Call MPI_ISEND(energydist(ijk1, 0, 0, 0), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV(buffer(ijk2, ii, 0, 0), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
      ijk1 = (ttm%ntcell(1) + 2) * (2 + (ttm%ntcell(2) + 2))
      ijk2 = 2 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      If (domain%idx == 0) Then
        ii = 1
      Else
        ii = 0
      End If
      Call MPI_ISEND(energydist(ijk1, 0, 0, 0), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_IRECV(buffer(ijk2, ii, 0, 0), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
      Call MPI_WAITALL(4, req, stat, ierr)

      ! -y/+y direction
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + (ttm%ntcell(2) + 2))
      If (domain%idy == domain%ny - 1) Then
        jj = -1
      Else
        jj = 0
      End If
      Call MPI_ISEND(energydist(ijk1, 0, 0, 0), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV(buffer(ijk2, 0, jj, 0), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2) + 2))
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2))
      If (domain%idy == 0) Then
        jj = 1
      Else
        jj = 0
      End If
      Call MPI_ISEND(energydist(ijk1, 0, 0, 0), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_IRECV(buffer(ijk2, 0, jj, 0), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
      Call MPI_WAITALL(4, req, stat, ierr)

      ! -z/+z direction
      ijk1 = 1
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * ttm%ntcell(3)
      If (domain%idz == domain%nz - 1) Then
        kk = -1
      Else
        kk = 0
      End If
      Call MPI_ISEND(energydist(ijk1, 0, 0, 0), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV(buffer(ijk2, 0, 0, kk), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
      ijk1 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1)
      ijk2 = 1 + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2)
      If (domain%idz == 0) Then
        kk = 1
      Else
        kk = 0
      End If
      Call MPI_ISEND(energydist(ijk1, 0, 0, 0), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_IRECV(buffer(ijk2, 0, 0, kk), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
      Call MPI_WAITALL(4, req, stat, ierr)

      energydist = energydist + buffer
      Deallocate (buffer)
    Else
      ! -x/+x direction
      Do k = 1, ttm%ntcell(3)
        Do j = 1, ttm%ntcell(2)
          ijk1 = 1 + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
          ijk2 = 1 + ttm%ntcell(1) + (ttm%ntcell(1) + 2) * (1 + j + (ttm%ntcell(2) + 2) * k)
          energydist(ijk2, -1, 0, 0) = energydist(ijk1, 0, 0, 0)
          ijk1 = (ttm%ntcell(1) + 2) * (1 + j + (ttm%ntcell(2) + 2) * k)
          ijk2 = 2 + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * k)
          energydist(ijk2, 1, 0, 0) = energydist(ijk1, 0, 0, 0)
        End Do
      End Do
      ! -y/+y direction
      Do k = 1, ttm%ntcell(3)
        Do i = 0, ttm%ntcell(1) + 1
          ijk1 = 1 + i + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 2) * k
          ijk2 = 1 + i + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + (ttm%ntcell(2) + 2) * k)
          energydist(ijk2, 0, -1, 0) = energydist(ijk1, 0, 0, 0)
          ijk1 = 1 + i + (ttm%ntcell(1) + 2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2) + 2) * k)
          ijk2 = 1 + i + (ttm%ntcell(1) + 2) * (1 + (ttm%ntcell(2) + 2) * k)
          energydist(ijk2, 0, 1, 0) = energydist(ijk1, 0, 0, 0)
        End Do
      End Do
      ! -z/+z direction
      Do j = 0, ttm%ntcell(2) + 1
        Do i = 0, ttm%ntcell(1) + 1
          ijk1 = 1 + i + (ttm%ntcell(1) + 2) * j
          ijk2 = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * ttm%ntcell(3))
          energydist(ijk2, 0, 0, -1) = energydist(ijk1, 0, 0, 0)
          ijk1 = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2) * (ttm%ntcell(3) + 1))
          ijk2 = 1 + i + (ttm%ntcell(1) + 2) * (j + (ttm%ntcell(2) + 2))
          energydist(ijk2, 0, 1, 1) = energydist(ijk1, 0, 0, 0)
        End Do
      End Do
    End If

    ! work through electronic temperature config%cells within ionic temperature voxels
    ! and in immediate borders: determine electronic temperatures needed to add
    ! redistributed energies and identify relevant config%cells - split according to
    ! available functional form of heat capacity

    Select Case (ttm%CeType)
    Case (TTM_CE_CONST, TTM_CE_CONST_DYN)
      ! constant specific heat capacity
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, ttm%numcell
              If (Abs(energydist(ijk, ii, jj, kk)) > zero_plus) Then
                energy_per_cell = energydist(ijk, ii, jj, kk)
                start_Te = ttm%eltemp(ijk, ii, jj, kk)
                end_Te = start_Te + energy_per_cell / Ce0a
                ttm%eltemp_adj(ijk, ii, jj, kk) = end_Te
                ttm%adjust(ijk, ii, jj, kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    Case (TTM_CE_TANH, TTM_CE_TANH_DYN)
      ! hyperbolic tangent specific heat capacity
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, ttm%numcell
              If (Abs(energydist(ijk, ii, jj, kk)) > zero_plus) Then
                energy_per_cell = energydist(ijk, ii, jj, kk)
                start_Te = ttm%eltemp(ijk, ii, jj, kk)
                increase = Cosh(ttm%sh_B * start_Te) * Exp(ttm%sh_B * energy_per_cell / sh_Aa)
                ! using equivalent function: Acosh(x)=Log(x+Sqrt((x-1.0)*(x+1.0)))
                end_Te = Log(increase + Sqrt((increase - 1.0_wp) * (increase + 1.0_wp))) / ttm%sh_B
                ttm%eltemp_adj(ijk, ii, jj, kk) = end_Te
                ttm%adjust(ijk, ii, jj, kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    Case (TTM_CE_LINEAR, TTM_CE_LINEAR_DYN)
      ! linear specific heat capacity to Fermi temperature
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, ttm%numcell
              If (Abs(energydist(ijk, ii, jj, kk)) > zero_plus) Then
                energy_per_cell = energydist(ijk, ii, jj, kk)
                start_Te = ttm%eltemp(ijk, ii, jj, kk)
                If (energy_per_cell > zero_plus) Then
                  end_Te = Sqrt(start_Te * start_Te + 2.0_wp * energy_per_cell * ttm%Tfermi / Cemaxa)
               If (end_Te > ttm%Tfermi) end_Te = 0.5_wp * (start_Te * start_Te / ttm%Tfermi + ttm%Tfermi) + energy_per_cell / Cemaxa
                Else
                  end_Te = start_Te + energy_per_cell / Cemaxa
                  If (end_Te < ttm%Tfermi) end_Te = Sqrt(ttm%Tfermi * (2.0_wp * (start_Te + energy_per_cell / Cemaxa) - ttm%Tfermi))
                End If
                ttm%eltemp_adj(ijk, ii, jj, kk) = end_Te
                ttm%adjust(ijk, ii, jj, kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    Case (TTM_CE_TABULATED, TTM_CE_TABULATED_DYN)
      ! tabulated ttm%volumetric heat capacity or more complex
      ! function: find new temperature iteratively by
      ! gradual integration (1 kelvin at a time)
      ! and interpolate within last kelvin
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, ttm%numcell
              If (Abs(energydist(ijk, ii, jj, kk)) > zero_plus) Then
                energy_per_cell = energydist(ijk, ii, jj, kk)
                start_Te = ttm%eltemp(ijk, ii, jj, kk)
                sgnplus = Sign(1.0_wp, energy_per_cell)
                energy_diff = sgnplus * energy_per_cell
                oldCe = Ce(start_Te, ttm)
                ! increase/decrease temperature of electronic config%cell by 1 kelvin until
                ! required energy change is reached or exceeded
                Do While (energy_diff >= zero_plus)
                  newCe = Ce(start_Te + sgnplus, ttm)
                  increase = 0.5_wp * (oldCe + newCe)
                  energy_diff = energy_diff - increase
                  start_Te = start_Te + sgnplus
                  oldCe = newCe
                End Do
                ! now interpolate for new temperature by using average specific heat
                ! over surrounding one-kelvin interval
                energy_diff = energy_diff + increase
                start_Te = start_Te - sgnplus
                newCe = Ce(start_Te, ttm)
                end_Te = start_Te + 2.0_wp * sgnplus * energy_diff / (oldCe + newCe)
                ttm%eltemp_adj(ijk, ii, jj, kk) = end_Te
                ttm%adjust(ijk, ii, jj, kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    End Select

    Deallocate (energydist)

  End Subroutine redistribute_Te
End Module ttm_utils
