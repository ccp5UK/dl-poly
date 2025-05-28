Module npt_utils

  !!-----------------------------------------------------------------------
  !!
  !! dl_poly_5 module declaring NPT (flexible cell) utils.
  !!
  !! author    - h.l.devereux may 2025
  !!
  !!-----------------------------------------------------------------------

  
  Use configuration,   Only: configuration_type,&
                             getcom
  Use comms,           Only: comms_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error
  Use kinds,           Only: wp,&
                             STR_LEN
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: pbcshift
  Use rigid_bodies,    Only: rigid_bodies_type,&
                             rigid_bodies_coms
  Use shared_units,    Only: update_shared_units
  Use statistics,      Only: stats_type
  Use thermostat,      Only: thermostat_type,&
                             ENS_NPT_BERENDSEN,&
                             ENS_NPT_BERENDSEN_ANISO,&
                             ENS_NPT_LANGEVIN,&
                             ENS_NPT_LANGEVIN_ANISO,&
                             ENS_NPT_MTK,&
                             ENS_NPT_MTK_ANISO,&
                             ENS_NPT_NOSE_HOOVER,&
                             ENS_NPT_NOSE_HOOVER_ANISO
  Use timer,           Only: timer_type,&
                             start_timer,&
                             stop_timer

  Implicit None

  Public :: xscale

  Contains

  Subroutine xscale(config, tstep, thermo, stats, neigh, rigid, domain, tmr, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to scale initial positions with change in box shape
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Real(Kind=wp),            Intent(In   ) :: tstep
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(stats_type),         Intent(InOut) :: stats
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)         :: message
    Integer                    :: fail, i, irgd, j, jrgd, lrgd
    Real(Kind=wp)              :: a1, a2, a3, a5, a6, a9, b1, b2, b3, b5, b6, b9, com(1:3), scale, &
                                  x, xa, y, ya, z, za
    Real(Kind=wp), Allocatable :: rgdxin(:), rgdyin(:), rgdzin(:)

#ifdef CHRONO
    Call start_timer(tmr, 'xscale')
#endif

    If (.not. thermo%variable_cell) Then 
#ifdef CHRONO
      Call stop_timer(tmr, 'xscale')
#endif
      Return
    End If
    If (.not. rigid%on) Then

      If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

        ! berendsen npt/nst

        If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

          scale = thermo%eta(1)

          Do i = 1, config%natms
            stats%xin(i) = scale * stats%xin(i)
            stats%yin(i) = scale * stats%yin(i)
            stats%zin(i) = scale * stats%zin(i)
          End Do

        Else

          Do i = 1, config%natms
            xa = stats%xin(i) * thermo%eta(1) + stats%yin(i) * thermo%eta(2) + stats%zin(i) * thermo%eta(3)
            ya = stats%xin(i) * thermo%eta(4) + stats%yin(i) * thermo%eta(5) + stats%zin(i) * thermo%eta(6)
            za = stats%xin(i) * thermo%eta(7) + stats%yin(i) * thermo%eta(8) + stats%zin(i) * thermo%eta(9)

            stats%xin(i) = xa
            stats%yin(i) = ya
            stats%zin(i) = za
          End Do

        End If

      Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

        ! hoover npt/nst

        Call getcom(stats%xin, stats%yin, stats%zin, config, com, comm)

        If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

          scale = Exp(tstep * thermo%eta(1))

          Do i = 1, config%natms
            stats%xin(i) = scale * (stats%xin(i) - com(1)) + com(1)
            stats%yin(i) = scale * (stats%yin(i) - com(2)) + com(2)
            stats%zin(i) = scale * (stats%zin(i) - com(3)) + com(3)
          End Do

        Else

          ! second order taylor expansion of Exp(tstep*thermo%eta)

          a1 = tstep * thermo%eta(1)
          a2 = tstep * thermo%eta(2)
          a3 = tstep * thermo%eta(3)
          a5 = tstep * thermo%eta(5)
          a6 = tstep * thermo%eta(6)
          a9 = tstep * thermo%eta(9)

          b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
          b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
          b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
          b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
          b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
          b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

          Do i = 1, config%natms
            xa = stats%xin(i) - com(1)
            ya = stats%yin(i) - com(2)
            za = stats%zin(i) - com(3)

            stats%xin(i) = xa * b1 + ya * b2 + za * b3 + com(1)
            stats%yin(i) = xa * b2 + ya * b5 + za * b6 + com(2)
            stats%zin(i) = xa * b3 + ya * b6 + za * b9 + com(3)
          End Do

        End If

      Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
               thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
               thermo%ensemble == ENS_NPT_MTK .or. &
               thermo%ensemble == ENS_NPT_MTK_ANISO) Then

        ! Langevin and MTK npt/nst

        If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

          scale = Exp(tstep * thermo%eta(1))

          Do i = 1, config%natms
            stats%xin(i) = scale * stats%xin(i)
            stats%yin(i) = scale * stats%yin(i)
            stats%zin(i) = scale * stats%zin(i)
          End Do

        Else

          ! second order taylor expansion of Exp(tstep*thermo%eta)

          a1 = tstep * thermo%eta(1)
          a2 = tstep * thermo%eta(2)
          a3 = tstep * thermo%eta(3)
          a5 = tstep * thermo%eta(5)
          a6 = tstep * thermo%eta(6)
          a9 = tstep * thermo%eta(9)

          b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
          b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
          b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
          b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
          b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
          b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

          Do i = 1, config%natms
            xa = stats%xin(i)
            ya = stats%yin(i)
            za = stats%zin(i)

            stats%xin(i) = xa * b1 + ya * b2 + za * b3
            stats%yin(i) = xa * b2 + ya * b5 + za * b6
            stats%zin(i) = xa * b3 + ya * b6 + za * b9
          End Do

        End If

      End If

      If (.not. neigh%update) Then

        If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

          ! berendsen npt/nst

          If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

            scale = thermo%eta(1)

            Do i = 1, config%natms
              neigh%xbg(i) = scale * neigh%xbg(i)
              neigh%ybg(i) = scale * neigh%ybg(i)
              neigh%zbg(i) = scale * neigh%zbg(i)
            End Do

          Else

            Do i = 1, config%natms
              xa = neigh%xbg(i) * thermo%eta(1) + neigh%ybg(i) * thermo%eta(2) + neigh%zbg(i) * thermo%eta(3)
              ya = neigh%xbg(i) * thermo%eta(4) + neigh%ybg(i) * thermo%eta(5) + neigh%zbg(i) * thermo%eta(6)
              za = neigh%xbg(i) * thermo%eta(7) + neigh%ybg(i) * thermo%eta(8) + neigh%zbg(i) * thermo%eta(9)

              neigh%xbg(i) = xa
              neigh%ybg(i) = ya
              neigh%zbg(i) = za
            End Do

          End If

        Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

          ! hoover npt/nst

          Call getcom(neigh%xbg, neigh%ybg, neigh%zbg, config, com, comm)

          If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

            scale = Exp(tstep * thermo%eta(1))

            Do i = 1, config%natms
              neigh%xbg(i) = scale * (neigh%xbg(i) - com(1)) + com(1)
              neigh%ybg(i) = scale * (neigh%ybg(i) - com(2)) + com(2)
              neigh%zbg(i) = scale * (neigh%zbg(i) - com(3)) + com(3)
            End Do

          Else

            ! second order taylor expansion of Exp(tstep*thermo%eta)

            a1 = tstep * thermo%eta(1)
            a2 = tstep * thermo%eta(2)
            a3 = tstep * thermo%eta(3)
            a5 = tstep * thermo%eta(5)
            a6 = tstep * thermo%eta(6)
            a9 = tstep * thermo%eta(9)

            b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
            b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
            b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
            b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
            b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
            b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

            Do i = 1, config%natms
              xa = neigh%xbg(i) - com(1)
              ya = neigh%ybg(i) - com(2)
              za = neigh%zbg(i) - com(3)

              neigh%xbg(i) = xa * b1 + ya * b2 + za * b3 + com(1)
              neigh%ybg(i) = xa * b2 + ya * b5 + za * b6 + com(2)
              neigh%zbg(i) = xa * b3 + ya * b6 + za * b9 + com(3)
            End Do

          End If

        Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
                 thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
                 thermo%ensemble == ENS_NPT_MTK .or. &
                 thermo%ensemble == ENS_NPT_MTK_ANISO) Then

          ! Langevin and MTK npt/nst

          If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

            scale = Exp(tstep * thermo%eta(1))

            Do i = 1, config%natms
              neigh%xbg(i) = scale * neigh%xbg(i)
              neigh%ybg(i) = scale * neigh%ybg(i)
              neigh%zbg(i) = scale * neigh%zbg(i)
            End Do

          Else

            ! second order taylor expansion of Exp(tstep*thermo%eta)

            a1 = tstep * thermo%eta(1)
            a2 = tstep * thermo%eta(2)
            a3 = tstep * thermo%eta(3)
            a5 = tstep * thermo%eta(5)
            a6 = tstep * thermo%eta(6)
            a9 = tstep * thermo%eta(9)

            b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
            b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
            b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
            b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
            b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
            b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

            Do i = 1, config%natms
              xa = neigh%xbg(i)
              ya = neigh%ybg(i)
              za = neigh%zbg(i)

              neigh%xbg(i) = xa * b1 + ya * b2 + za * b3
              neigh%ybg(i) = xa * b2 + ya * b5 + za * b6
              neigh%zbg(i) = xa * b3 + ya * b6 + za * b9
            End Do

          End If

        End If

      End If

    Else ! RBs exist

      fail = 0
      Allocate (rgdxin(1:rigid%max_rigid), rgdyin(1:rigid%max_rigid), rgdzin(1:rigid%max_rigid), Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'xscale allocation failure'
        Call error(0, message)
      End If

      ! Halo initial RB members positions across onto neighbouring domains
      ! to get initial COMs

      If (rigid%share) Then
        Call update_shared_units(config, rigid%list_shared, &
                                 rigid%map_shared, stats%xin, stats%yin, stats%zin, domain, comm)
      End If
      Call rigid_bodies_coms(config, stats%xin, stats%yin, stats%zin, rgdxin, rgdyin, rgdzin, rigid)

      If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

        ! berendsen npt/nst

        If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

          scale = thermo%eta(1)

          Do j = 1, config%nfree
            i = config%lstfre(j)

            stats%xin(i) = scale * stats%xin(i)
            stats%yin(i) = scale * stats%yin(i)
            stats%zin(i) = scale * stats%zin(i)
          End Do

          Do irgd = 1, rigid%n_types
            x = rgdxin(irgd)
            y = rgdyin(irgd)
            z = rgdzin(irgd)

            rgdxin(irgd) = scale * rgdxin(irgd)
            rgdyin(irgd) = scale * rgdyin(irgd)
            rgdzin(irgd) = scale * rgdzin(irgd)

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd)

              If (i <= config%natms) Then
                stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
              End If
            End Do
          End Do

        Else

          Do j = 1, config%nfree
            i = config%lstfre(j)

            xa = stats%xin(i) * thermo%eta(1) + stats%yin(i) * thermo%eta(2) + stats%zin(i) * thermo%eta(3)
            ya = stats%xin(i) * thermo%eta(4) + stats%yin(i) * thermo%eta(5) + stats%zin(i) * thermo%eta(6)
            za = stats%xin(i) * thermo%eta(7) + stats%yin(i) * thermo%eta(8) + stats%zin(i) * thermo%eta(9)

            stats%xin(i) = xa
            stats%yin(i) = ya
            stats%zin(i) = za
          End Do

          Do irgd = 1, rigid%n_types
            x = rgdxin(irgd)
            y = rgdyin(irgd)
            z = rgdzin(irgd)

            xa = rgdxin(irgd) * thermo%eta(1) + rgdyin(irgd) * thermo%eta(2) + rgdzin(irgd) * thermo%eta(3)
            ya = rgdxin(irgd) * thermo%eta(4) + rgdyin(irgd) * thermo%eta(5) + rgdzin(irgd) * thermo%eta(6)
            za = rgdxin(irgd) * thermo%eta(7) + rgdyin(irgd) * thermo%eta(8) + rgdzin(irgd) * thermo%eta(9)

            rgdxin(irgd) = xa
            rgdyin(irgd) = ya
            rgdzin(irgd) = za

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd)

              If (i <= config%natms) Then
                stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
              End If
            End Do
          End Do

        End If

      Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

        ! hoover npt/nst

        Call getcom(stats%xin, stats%yin, stats%zin, config, com, comm)

        If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

          scale = Exp(tstep * thermo%eta(1))

          Do j = 1, config%nfree
            i = config%lstfre(j)

            stats%xin(i) = scale * (stats%xin(i) - com(1)) + com(1)
            stats%yin(i) = scale * (stats%yin(i) - com(2)) + com(2)
            stats%zin(i) = scale * (stats%zin(i) - com(3)) + com(3)
          End Do

          Do irgd = 1, rigid%n_types
            x = rgdxin(irgd)
            y = rgdyin(irgd)
            z = rgdzin(irgd)

            rgdxin(irgd) = scale * (rgdxin(irgd) - com(1)) + com(1)
            rgdyin(irgd) = scale * (rgdyin(irgd) - com(2)) + com(2)
            rgdzin(irgd) = scale * (rgdzin(irgd) - com(3)) + com(3)

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd)

              If (i <= config%natms) Then
                stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
              End If
            End Do
          End Do

        Else

          ! second order taylor expansion of Exp(tstep*thermo%eta)

          a1 = tstep * thermo%eta(1)
          a2 = tstep * thermo%eta(2)
          a3 = tstep * thermo%eta(3)
          a5 = tstep * thermo%eta(5)
          a6 = tstep * thermo%eta(6)
          a9 = tstep * thermo%eta(9)

          b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
          b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
          b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
          b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
          b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
          b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

          Do j = 1, config%nfree
            i = config%lstfre(j)

            xa = stats%xin(i) - com(1)
            ya = stats%yin(i) - com(2)
            za = stats%zin(i) - com(3)

            stats%xin(i) = xa * b1 + ya * b2 + za * b3 + com(1)
            stats%yin(i) = xa * b2 + ya * b5 + za * b6 + com(2)
            stats%zin(i) = xa * b3 + ya * b6 + za * b9 + com(3)
          End Do

          Do irgd = 1, rigid%n_types
            x = rgdxin(irgd)
            y = rgdyin(irgd)
            z = rgdzin(irgd)

            xa = rgdxin(irgd) - com(1)
            ya = rgdyin(irgd) - com(2)
            za = rgdzin(irgd) - com(3)

            rgdxin(irgd) = xa * b1 + ya * b2 + za * b3 + com(1)
            rgdyin(irgd) = xa * b2 + ya * b5 + za * b6 + com(2)
            rgdzin(irgd) = xa * b3 + ya * b6 + za * b9 + com(3)

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd)

              If (i <= config%natms) Then
                stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
              End If
            End Do
          End Do

        End If

      Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
               thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
               thermo%ensemble == ENS_NPT_MTK .or. &
               thermo%ensemble == ENS_NPT_MTK_ANISO) Then

        ! Langevin and MTK npt/nst

        If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

          scale = Exp(tstep * thermo%eta(1))

          Do j = 1, config%nfree
            i = config%lstfre(j)

            stats%xin(i) = scale * stats%xin(i)
            stats%yin(i) = scale * stats%yin(i)
            stats%zin(i) = scale * stats%zin(i)
          End Do

          Do irgd = 1, rigid%n_types
            x = rgdxin(irgd)
            y = rgdyin(irgd)
            z = rgdzin(irgd)

            rgdxin(irgd) = scale * rgdxin(irgd)
            rgdyin(irgd) = scale * rgdyin(irgd)
            rgdzin(irgd) = scale * rgdzin(irgd)

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd)

              If (i <= config%natms) Then
                stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
              End If
            End Do
          End Do

        Else

          ! second order taylor expansion of Exp(tstep*thermo%eta)

          a1 = tstep * thermo%eta(1)
          a2 = tstep * thermo%eta(2)
          a3 = tstep * thermo%eta(3)
          a5 = tstep * thermo%eta(5)
          a6 = tstep * thermo%eta(6)
          a9 = tstep * thermo%eta(9)

          b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
          b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
          b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
          b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
          b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
          b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

          Do j = 1, config%nfree
            i = config%lstfre(j)

            xa = stats%xin(i)
            ya = stats%yin(i)
            za = stats%zin(i)

            stats%xin(i) = xa * b1 + ya * b2 + za * b3
            stats%yin(i) = xa * b2 + ya * b5 + za * b6
            stats%zin(i) = xa * b3 + ya * b6 + za * b9
          End Do

          Do irgd = 1, rigid%n_types
            x = rgdxin(irgd)
            y = rgdyin(irgd)
            z = rgdzin(irgd)

            xa = rgdxin(irgd)
            ya = rgdyin(irgd)
            za = rgdzin(irgd)

            rgdxin(irgd) = xa * b1 + ya * b2 + za * b3
            rgdyin(irgd) = xa * b2 + ya * b5 + za * b6
            rgdzin(irgd) = xa * b3 + ya * b6 + za * b9

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd)

              If (i <= config%natms) Then
                stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
              End If
            End Do
          End Do

        End If

      End If

      If (.not. neigh%update) Then

        Call rigid_bodies_coms(config, neigh%xbg, neigh%ybg, neigh%zbg, rgdxin, rgdyin, rgdzin, rigid)

        If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

          ! berendsen npt/nst

          If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

            scale = thermo%eta(1)

            Do j = 1, config%nfree
              i = config%lstfre(j)

              neigh%xbg(i) = scale * neigh%xbg(i)
              neigh%ybg(i) = scale * neigh%ybg(i)
              neigh%zbg(i) = scale * neigh%zbg(i)
            End Do

            Do irgd = 1, rigid%n_types
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              rgdxin(irgd) = scale * rgdxin(irgd)
              rgdyin(irgd) = scale * rgdyin(irgd)
              rgdzin(irgd) = scale * rgdzin(irgd)

              lrgd = rigid%list(-1, irgd)
              Do jrgd = 1, lrgd
                i = rigid%index_local(jrgd, irgd)

                If (i <= config%natms) Then
                  neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                  neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                  neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                End If
              End Do
            End Do

          Else

            Do j = 1, config%nfree
              i = config%lstfre(j)

              xa = neigh%xbg(i) * thermo%eta(1) + neigh%ybg(i) * thermo%eta(2) + neigh%zbg(i) * thermo%eta(3)
              ya = neigh%xbg(i) * thermo%eta(4) + neigh%ybg(i) * thermo%eta(5) + neigh%zbg(i) * thermo%eta(6)
              za = neigh%xbg(i) * thermo%eta(7) + neigh%ybg(i) * thermo%eta(8) + neigh%zbg(i) * thermo%eta(9)

              neigh%xbg(i) = xa
              neigh%ybg(i) = ya
              neigh%zbg(i) = za
            End Do

            Do irgd = 1, rigid%n_types
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd) * thermo%eta(1) + rgdyin(irgd) * thermo%eta(2) + rgdzin(irgd) * thermo%eta(3)
              ya = rgdxin(irgd) * thermo%eta(4) + rgdyin(irgd) * thermo%eta(5) + rgdzin(irgd) * thermo%eta(6)
              za = rgdxin(irgd) * thermo%eta(7) + rgdyin(irgd) * thermo%eta(8) + rgdzin(irgd) * thermo%eta(9)

              rgdxin(irgd) = xa
              rgdyin(irgd) = ya
              rgdzin(irgd) = za

              lrgd = rigid%list(-1, irgd)
              Do jrgd = 1, lrgd
                i = rigid%index_local(jrgd, irgd)

                If (i <= config%natms) Then
                  neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                  neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                  neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                End If
              End Do
            End Do

          End If

        Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

          ! hoover npt/nst

          Call getcom(neigh%xbg, neigh%ybg, neigh%zbg, config, com, comm)

          If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

            scale = Exp(tstep * thermo%eta(1))

            Do j = 1, config%nfree
              i = config%lstfre(j)

              neigh%xbg(i) = scale * (neigh%xbg(i) - com(1)) + com(1)
              neigh%ybg(i) = scale * (neigh%ybg(i) - com(2)) + com(2)
              neigh%zbg(i) = scale * (neigh%zbg(i) - com(3)) + com(3)
            End Do

            Do irgd = 1, rigid%n_types
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              rgdxin(irgd) = scale * (rgdxin(irgd) - com(1)) + com(1)
              rgdyin(irgd) = scale * (rgdyin(irgd) - com(2)) + com(2)
              rgdzin(irgd) = scale * (rgdzin(irgd) - com(3)) + com(3)

              lrgd = rigid%list(-1, irgd)
              Do jrgd = 1, lrgd
                i = rigid%index_local(jrgd, irgd)

                If (i <= config%natms) Then
                  neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                  neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                  neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                End If
              End Do
            End Do

          Else

            ! second order taylor expansion of Exp(tstep*thermo%eta)

            a1 = tstep * thermo%eta(1)
            a2 = tstep * thermo%eta(2)
            a3 = tstep * thermo%eta(3)
            a5 = tstep * thermo%eta(5)
            a6 = tstep * thermo%eta(6)
            a9 = tstep * thermo%eta(9)

            b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
            b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
            b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
            b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
            b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
            b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

            Do j = 1, config%nfree
              i = config%lstfre(j)

              xa = neigh%xbg(i) - com(1)
              ya = neigh%ybg(i) - com(2)
              za = neigh%zbg(i) - com(3)

              neigh%xbg(i) = xa * b1 + ya * b2 + za * b3 + com(1)
              neigh%ybg(i) = xa * b2 + ya * b5 + za * b6 + com(2)
              neigh%zbg(i) = xa * b3 + ya * b6 + za * b9 + com(3)
            End Do

            Do irgd = 1, rigid%n_types
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd) - com(1)
              ya = rgdyin(irgd) - com(2)
              za = rgdzin(irgd) - com(3)

              rgdxin(irgd) = xa * b1 + ya * b2 + za * b3 + com(1)
              rgdyin(irgd) = xa * b2 + ya * b5 + za * b6 + com(2)
              rgdzin(irgd) = xa * b3 + ya * b6 + za * b9 + com(3)

              lrgd = rigid%list(-1, irgd)
              Do jrgd = 1, lrgd
                i = rigid%index_local(jrgd, irgd)

                If (i <= config%natms) Then
                  neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                  neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                  neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                End If
              End Do
            End Do

          End If

        Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
                 thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
                 thermo%ensemble == ENS_NPT_MTK .or. &
                 thermo%ensemble == ENS_NPT_MTK_ANISO) Then

          ! Langevin and MTK npt/nst

          If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

            scale = Exp(tstep * thermo%eta(1))

            Do j = 1, config%nfree
              i = config%lstfre(j)

              neigh%xbg(i) = scale * neigh%xbg(i)
              neigh%ybg(i) = scale * neigh%ybg(i)
              neigh%zbg(i) = scale * neigh%zbg(i)
            End Do

            Do irgd = 1, rigid%n_types
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              rgdxin(irgd) = scale * rgdxin(irgd)
              rgdyin(irgd) = scale * rgdyin(irgd)
              rgdzin(irgd) = scale * rgdzin(irgd)

              lrgd = rigid%list(-1, irgd)
              Do jrgd = 1, lrgd
                i = rigid%index_local(jrgd, irgd)

                If (i <= config%natms) Then
                  neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                  neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                  neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                End If
              End Do
            End Do

          Else

            ! second order taylor expansion of Exp(tstep*thermo%eta)

            a1 = tstep * thermo%eta(1)
            a2 = tstep * thermo%eta(2)
            a3 = tstep * thermo%eta(3)
            a5 = tstep * thermo%eta(5)
            a6 = tstep * thermo%eta(6)
            a9 = tstep * thermo%eta(9)

            b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
            b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
            b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
            b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
            b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
            b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

            Do j = 1, config%nfree
              i = config%lstfre(j)

              xa = neigh%xbg(i)
              ya = neigh%ybg(i)
              za = neigh%zbg(i)

              neigh%xbg(i) = xa * b1 + ya * b2 + za * b3
              neigh%ybg(i) = xa * b2 + ya * b5 + za * b6
              neigh%zbg(i) = xa * b3 + ya * b6 + za * b9
            End Do

            Do irgd = 1, rigid%n_types
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd)
              ya = rgdyin(irgd)
              za = rgdzin(irgd)

              rgdxin(irgd) = xa * b1 + ya * b2 + za * b3
              rgdyin(irgd) = xa * b2 + ya * b5 + za * b6
              rgdzin(irgd) = xa * b3 + ya * b6 + za * b9

              lrgd = rigid%list(-1, irgd)
              Do jrgd = 1, lrgd
                i = rigid%index_local(jrgd, irgd)

                If (i <= config%natms) Then
                  neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                  neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                  neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                End If
              End Do
            End Do

          End If

        End If

        ! Halo final RB members positions across onto neighbouring domains

        If (rigid%share) Then
          Call update_shared_units(config, rigid%list_shared, &
                                   rigid%map_shared, neigh%xbg, neigh%ybg, neigh%zbg, domain, comm)
        End If
      End If

      Deallocate (rgdxin, rgdyin, rgdzin, Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'xscale deallocation failure, node'
        Call error(0, message)
      End If

    End If

    Call pbcshift(config%imcon, config%cell, config%natms, stats%xin, stats%yin, stats%zin)
    If (neigh%unconditional_update) Call pbcshift(config%imcon, config%cell, config%natms, neigh%xbg, neigh%ybg, neigh%zbg)

#ifdef CHRONO
    Call stop_timer(tmr, 'xscale')
#endif

  End Subroutine xscale
End Module