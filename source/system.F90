Module system

  Use angles,          Only: angles_type
  Use bonds,           Only: bonds_type
  Use comms,           Only: &
                             Revive_tag, SysExpand_tag, comm_self, comms_type, gbcast, gcheck, &
                             grecv, gsend, gsum, gsync, gtime, mode_create, mode_wronly, &
                             offset_kind
  Use configuration,   Only: configuration_type,&
                             write_config,&
                             IMCON_NOPBC,&
                             IMCON_CUBIC,&
                             IMCON_ORTHORHOMBIC,&
                             IMCON_PARALLELOPIPED,&
                             IMCON_SLAB,&
                             IMCON_TRUNC_OCTO,&
                             IMCON_RHOMBIC_DODEC,&
                             IMCON_HEXAGONAL
  Use constants,       Only: engunit,&
                             nmpldt,&
                             zero_plus
  Use constraints,     Only: constraints_type
  Use core_shell,      Only: core_shell_type
  Use development,     Only: development_type
  Use dihedrals,       Only: dihedrals_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use filename,        Only: FILE_CONFIG,&
                             FILE_FIELD,&
                             FILE_REVCON,&
                             FILE_REVIVE,&
                             FILE_REVOLD,&
                             file_type
  Use flow_control,    Only: RESTART_KEY_NOSCALE,&
                             RESTART_KEY_OLD
  Use greenkubo,       Only: greenkubo_type
  Use inversions,      Only: inversions_type
  Use io,              Only: &
                             IO_ALLOCATION_ERROR, IO_BASE_COMM_NOT_SET, IO_RESTART, &
                             IO_UNKNOWN_WRITE_LEVEL, IO_UNKNOWN_WRITE_OPTION, &
                             IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
                             IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_NETCDF, &
                             IO_WRITE_UNSORTED_DIRECT, IO_WRITE_UNSORTED_MASTER, &
                             IO_WRITE_UNSORTED_MPIIO, io_close, io_delete, io_finalize, &
                             io_get_parameters, io_init, io_nc_create, io_nc_put_var, io_open, &
                             io_set_parameters, io_type, io_write_record, io_write_sorted_file, &
                             recsz
  Use kinds,           Only: li,&
                             wp
  Use metal,           Only: metal_lrc,&
                             metal_type
  Use netcdf_wrap,     Only: netcdf_param
  Use numerics,        Only: dcell,&
                             images
  Use parse,           Only: get_word,&
                             lower_case,&
                             strip_blanks,&
                             tabs_2_blanks,&
                             word_2_real
  Use rdfs,            Only: rdf_type
  Use rigid_bodies,    Only: rigid_bodies_type
  Use site,            Only: site_type
  Use statistics,      Only: stats_type
  Use thermostat,      Only: thermostat_type
  Use vdw,             Only: vdw_lrc,&
                             vdw_type
  Use z_density,       Only: z_density_type

  Implicit None
  Private
  Public :: system_revive
  Public :: system_init
  Public :: system_expand
Contains

  Subroutine system_init(rcut, keyres, time, tmst, nstep, &
                         stats, devel, green, thermo, met, bond, angle, dihedral, inversion, &
                         zdensity, sites, vdws, rdf, config, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reading the REVIVE file data and defining the
    ! initial thermodynamic and structural accumulators
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2016
    ! contrib   - m.a.seaton june 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),            Intent(In   ) :: rcut
    Integer,                  Intent(InOut) :: keyres
    Real(Kind=wp),            Intent(  Out) :: time, tmst
    Integer,                  Intent(  Out) :: nstep
    Type(stats_type),         Intent(InOut) :: stats
    Type(development_type),   Intent(In   ) :: devel
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(metal_type),         Intent(InOut) :: met
    Type(bonds_type),         Intent(InOut) :: bond
    Type(angles_type),        Intent(InOut) :: angle
    Type(dihedrals_type),     Intent(InOut) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(site_type),          Intent(InOut) :: sites
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(configuration_type), Intent(InOut) :: config
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=40) :: forma
    Integer           :: gidx, i, i_tmp, j, k, keyio, l
    Logical           :: l_tmp
    Real(Kind=wp)     :: dmxgusr, dncfang, dncfbnd, dncfdih, dncfinv, dncfrdf, dncfusr, dncfzdn, &
                         dnstep, dnumacc, drusr, dtstep, r_mxnode, xyz(0:6), dvafstep(1:green%samp)

    forma = ' '
    ! Define format for REVOLD reading in ASCII

    If (devel%l_rin) Then
      i = 64 / 4 - 1 ! Bit_Size(0.0_wp)/4 - 1
      j = Max(stats%mxstak * stats%mxnstk, rdf%max_grid * rdf%max_rdf, zdensity%max_grid, rdf%max_grid_usr, &
           & config%mxgana * config%mxtana)

      Write (forma, 10) j / 4 + 1, i + 9, i
      10 Format('(1p,', i0, '(/,4e', i0, '.', i0, 'E3))')
    End If

    ! Initialise read failure flag

    keyio = 0

    50 Continue
    If (keyres /= RESTART_KEY_OLD .or. comm%idnode /= 0) Then

      ! initialise step and time related accumulators

      nstep = 0
      dtstep = 0.0_wp
      time = 0.0_wp
      tmst = 0.0_wp

      ! initialise temperature and pressure coupling parameters
      ! and integral for conserved quantity

      thermo%chi_t = 0.0_wp
      thermo%cint = 0.0_wp
      thermo%chi_p = 0.0_wp
      thermo%eta = 0.0_wp

      ! initialise stats%strcon,stats%stress,stats%virtot and stats%vircon

      stats%virtot = 0.0_wp
      stats%stress = 0.0_wp
      stats%vircon = 0.0_wp
      stats%strcon = 0.0_wp
      stats%virpmf = 0.0_wp
      stats%strpmf = 0.0_wp

      ! initialise accumulator arrays if reading failure occurred

      If (keyio > 0) Then

        stats%numacc = 0
        stats%stpval = 0.0_wp
        stats%stpvl0 = 0.0_wp
        stats%sumval = 0.0_wp
        stats%ssqval = 0.0_wp
        stats%zumval = 0.0_wp
        stats%ravval = 0.0_wp
        stats%stkval = 0.0_wp

        If (rdf%l_collect) Then
          rdf%n_configs = 0
          rdf%rdf = 0.0_wp
        End If

        If (rdf%max_grid_usr > 0) Then
          rdf%cutoff_usr = 0.0_wp
          rdf%n_configs_usr = 0
          rdf%usr = 0
        End If

        If (zdensity%l_collect) Then
          zdensity%n_samples = 0
          zdensity%density = 0.0_wp
        End If

        If (green%samp > 0) Then
          green%vafcount = 0.0_wp
          green%step = 0
          green%vafdata = 0.0_wp
          green%vaf = 0.0_wp
          green%time = 0.0_wp
        End If

        If (bond%bin_pdf > 0) Then
          bond%n_frames = 0
          bond%dst = 0.0_wp
        End If

        If (angle%bin_adf > 0) Then
          angle%n_frames = 0
          angle%dst = 0.0_wp
        End If

        If (dihedral%bin_adf > 0) Then
          dihedral%n_frames = 0
          dihedral%dst = 0.0_wp
        End If

        If (inversion%bin_adf > 0) Then
          inversion%n_frames = 0
          inversion%dst = 0.0_wp
        End If

      End If

    End If

    ! restart simulation and continue

    If (keyres == RESTART_KEY_OLD) Then

      ! If REVOLD doesn't exist then abort (mishmashed REVOLD is handled separately)

      l_tmp = .true.
      If (comm%idnode == 0) Inquire (File=files(FILE_REVOLD)%filename, Exist=l_tmp)
      Call gcheck(comm, l_tmp, "enforce")
      If (.not. l_tmp) Call error(519)

      ! Check REVOLD restart compatibility: rcut,rdf%rbin,config%megatm

      xyz(1:3) = 0.0_wp
      If (comm%idnode == 0) Then
        If (devel%l_rin) Then
          Open (Newunit=files(FILE_REVOLD)%unit_no, file=files(FILE_REVOLD)%filename, form='formatted', IOStat=keyio)
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) xyz(1), xyz(2), xyz(3)
        Else
          Open (Newunit=files(FILE_REVOLD)%unit_no, file=files(FILE_REVOLD)%filename, form='unformatted', IOStat=keyio)
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) xyz(1), xyz(2), xyz(3)
        End If
      End If
      Call gbcast(comm, xyz, 0)
      If (Abs(xyz(1) - rcut) > 1.0e-6_wp .or. Abs(xyz(2) - rdf%rbin) > 1.0e-6_wp .or. &
          Nint(xyz(3)) /= config%megatm) Call error(519)

      ! read the rest of the accumulator data from dump file

      If (comm%idnode == 0) Then
        If (devel%l_rin) Then
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) &
            dnstep, dtstep, time, tmst, dnumacc, thermo%chi_t, thermo%chi_p, thermo%cint
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) thermo%eta
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%stpval
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%stpvl0
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%sumval
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%ssqval
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%zumval
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%ravval
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%stkval
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%strcon
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%strpmf
          Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) stats%stress

          If (rdf%l_collect) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfrdf, rdf%rdf
          End If
          If (rdf%max_grid_usr > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dmxgusr, drusr, dncfusr, rdf%usr
          End If
          If (zdensity%l_collect) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfzdn, zdensity%density
          End If
          If (green%samp > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) green%vafcount
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dvafstep
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) green%vafdata
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) green%vaf
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) green%time
          End If

          If (bond%bin_pdf > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfbnd, bond%dst
          End If
          If (angle%bin_adf > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfang, angle%dst
          End If
          If (dihedral%bin_adf > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfdih, dihedral%dst
          End If
          If (inversion%bin_adf > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfinv, inversion%dst
          End If
        Else
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) &
            dnstep, dtstep, time, tmst, dnumacc, thermo%chi_t, thermo%chi_p, thermo%cint
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) thermo%eta
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%stpval
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%stpvl0
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%sumval
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%ssqval
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%zumval
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%ravval
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%stkval
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%strcon
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%strpmf
          Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) stats%stress

          If (rdf%l_collect) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dncfrdf, rdf%rdf
          If (rdf%max_grid_usr > 0) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dmxgusr, drusr, dncfusr, rdf%usr
          If (zdensity%l_collect) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dncfzdn, zdensity%density
          If (green%samp > 0) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) green%vafcount
            Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dvafstep
            Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) green%vafdata
            Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) green%vaf
            Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) green%time
          End If

          If (bond%bin_pdf > 0) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dncfbnd, bond%dst
          If (angle%bin_adf > 0) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dncfang, angle%dst
          If (dihedral%bin_adf > 0) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dncfdih, dihedral%dst
          If (inversion%bin_adf > 0) Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio, End=100) dncfinv, inversion%dst
        End If

        nstep = Nint(dnstep)

        stats%numacc = Nint(dnumacc)

        If (rdf%l_collect) rdf%n_configs = Nint(dncfrdf)
        If (rdf%max_grid_usr > 0) Then
          rdf%max_grid_usr = Nint(dmxgusr)
          rdf%cutoff_usr = drusr
          rdf%n_configs_usr = Nint(dncfusr)
        End If
        If (zdensity%l_collect) zdensity%n_samples = Nint(dncfzdn)
        If (green%samp > 0) green%step = Nint(dvafstep)

        If (bond%bin_pdf > 0) bond%n_frames = Nint(dncfbnd)
        If (angle%bin_adf > 0) angle%n_frames = Nint(dncfang)
        If (dihedral%bin_adf > 0) dihedral%n_frames = Nint(dncfdih)
        If (inversion%bin_adf > 0) inversion%n_frames = Nint(dncfinv)

        ! calculate stats%virtot = stats%virtot-stats%vircon-stats%virpmf

        stats%vircon = stats%stpval(17) * engunit
        stats%virpmf = stats%stpval(26) * engunit
        stats%virtot = (stats%stpval(12) - stats%stpval(17) - stats%stpval(26)) * engunit
      End If

      100 Continue

      ! If 'restart' is impossible go to 'restart noscale' and reinitialise

      Call gbcast(comm, keyio, 0)
      If (keyio /= 0) Then
        If (comm%idnode == 0) Then
          Call warning(190, 0.0_wp, 0.0_wp, 0.0_wp)
          Call files(FILE_REVOLD)%close ()
        End If
        keyres = RESTART_KEY_NOSCALE
        Go To 50
      End If

      ! broadcast stored variables

      If (comm%mxnode > 1) Then

        Call gbcast(comm, nstep, 0)
        Call gbcast(comm, dtstep, 0)
        Call gbcast(comm, time, 0)
        Call gbcast(comm, tmst, 0)
        Call gbcast(comm, stats%numacc, 0)
        Call gbcast(comm, thermo%chi_t, 0)
        Call gbcast(comm, thermo%chi_p, 0)
        Call gbcast(comm, thermo%cint, 0)
        Call gbcast(comm, thermo%eta, 0)
        Call gbcast(comm, stats%stpval, 0)
        Call gbcast(comm, stats%stpvl0, 0)
        Call gbcast(comm, stats%sumval, 0)
        Call gbcast(comm, stats%ssqval, 0)
        Call gbcast(comm, stats%zumval, 0)
        Call gbcast(comm, stats%ravval, 0)
        Do k = 0, stats%mxnstk
          Call gbcast(comm, stats%stkval(:, k), 0)
        End Do
        Call gbcast(comm, stats%strcon, 0)
        Call gbcast(comm, stats%strpmf, 0)
        Call gbcast(comm, stats%stress, 0)
        Call gbcast(comm, stats%vircon, 0)
        Call gbcast(comm, stats%virpmf, 0)
        Call gbcast(comm, stats%virtot, 0)

        ! Reset timestep

        thermo%tstep = dtstep

        r_mxnode = 1.0_wp / Real(comm%mxnode, wp)

        ! rdf%rdf table - broadcast and normalise

        If (rdf%l_collect) Then
          Call gbcast(comm, rdf%n_configs, 0)
          Do k = 1, rdf%max_rdf
            Call gbcast(comm, rdf%rdf(:, k), 0)
            rdf%rdf(:, k) = rdf%rdf(:, k) * r_mxnode
          End Do
        End If

        ! USR RDF table - broadcast and normalise

        If (rdf%max_grid_usr > 0) Then
          Call gbcast(comm, rdf%max_grid_usr, 0)
          Call gbcast(comm, rdf%cutoff_usr, 0)
          Call gbcast(comm, rdf%n_configs_usr, 0)
          Call gbcast(comm, rdf%usr, 0)
          rdf%usr(:) = rdf%usr(:) * r_mxnode
        End If

        ! z-density table - broadcast and normalise

        If (zdensity%l_collect) Then
          Call gbcast(comm, zdensity%n_samples, 0)
          Do k = 1, sites%mxatyp
            Call gbcast(comm, zdensity%density(:, k), 0)
            zdensity%density(:, k) = zdensity%density(:, k) * r_mxnode
          End Do
        End If

        ! green%vafdata table - broadcast and normalise

        If (green%samp > 0) Then
          Do j = 1, green%samp
            l = (j - 1) * (sites%mxatyp + 1)
            Do k = 1, sites%mxatyp + 1
              Call gbcast(comm, green%vafdata(:, l + k), 0)

              ! avoid normalising timing information

              If (k /= sites%mxatyp + 1) green%vafdata(:, l + k) = green%vafdata(:, l + k) * r_mxnode
            End Do
          End Do
        End If

        ! bonds table - broadcast and normalise

        If (bond%bin_pdf > 0) Then
          Call gbcast(comm, bond%n_frames, 0)
          Do k = 1, bond%ldf(0)
            Call gbcast(comm, bond%dst(:, k), 0)

            bond%dst(:, k) = bond%dst(:, k) * r_mxnode
          End Do
        End If

        ! angles table - broadcast and normalise

        If (angle%bin_adf > 0) Then
          Call gbcast(comm, angle%n_frames, 0)
          Do k = 1, angle%ldf(0)
            Call gbcast(comm, angle%dst(:, k), 0)

            angle%dst(:, k) = angle%dst(:, k) * r_mxnode
          End Do
        End If

        ! dihedrals table - broadcast and normalise

        If (dihedral%bin_adf > 0) Then
          Call gbcast(comm, dihedral%n_frames, 0)
          Do k = 1, dihedral%ldf(0)
            Call gbcast(comm, dihedral%bin_adf, 0)

            dihedral%dst(:, k) = dihedral%dst(:, k) * r_mxnode
          End Do
        End If

        ! inversions table - broadcast and normalise

        If (inversion%bin_adf > 0) Then
          Call gbcast(comm, inversion%n_frames, 0)
          Do k = 1, inversion%ldf(0)
            Call gbcast(comm, inversion%bin_adf, 0)

            inversion%dst(:, k) = inversion%dst(:, k) * r_mxnode
          End Do
        End If
      End If

    End If

    ! initialise initial positions to current positions
    ! and final displacements to zero

    Do i = 1, config%natms
      stats%xin(i) = config%parts(i)%xxx
      stats%yin(i) = config%parts(i)%yyy
      stats%zin(i) = config%parts(i)%zzz

      stats%xto(i) = 0.0_wp
      stats%yto(i) = 0.0_wp
      stats%zto(i) = 0.0_wp
    End Do

    If (keyres == RESTART_KEY_OLD) Then

      ! Error accumulator: keyio is still zero otherwise we cannot get here

      i_tmp = 0

      Do k = 1, config%megatm
        xyz = 0.0_wp

        If (comm%idnode == 0) Then
          If (devel%l_rin) Then
            Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio) &
              xyz(0), xyz(1), xyz(2), xyz(3), xyz(4), xyz(5), xyz(6)
          Else
            Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio) &
              xyz(0), xyz(1), xyz(2), xyz(3), xyz(4), xyz(5), xyz(6)
          End If
        End If
        If (keyio /= 0) i_tmp = 1

        Call gbcast(comm, xyz, 0)
        gidx = Nint(xyz(0))

        ! assign particle initial positions and final displacements
        ! to the corresponding domains

        Do i = 1, config%natms
          If (config%ltg(i) == gidx) Then
            stats%xin(i) = xyz(1)
            stats%yin(i) = xyz(2)
            stats%zin(i) = xyz(3)

            stats%xto(i) = xyz(4)
            stats%yto(i) = xyz(5)
            stats%zto(i) = xyz(6)
          End If
        End Do
      End Do

      ! If 'restart' is impossible go to 'restart noscale' and reinitialise

      Call gbcast(comm, i_tmp, 0)
      If (i_tmp /= 0) Then
        If (comm%idnode == 0) Then
          Call warning(190, 0.0_wp, 0.0_wp, 0.0_wp)
          Call files(FILE_REVOLD)%close ()
        End If
        keyres = RESTART_KEY_NOSCALE
        Go To 50
      End If

      ! Read velocities for VAF calculations if needed

      If (green%samp > 0) Then

        i_tmp = 0

        Do j = 1, green%samp
          Do k = 1, config%megatm
            xyz = 0.0_wp

            If (comm%idnode == 0) Then
              If (devel%l_rin) Then
                Read (Unit=files(FILE_REVOLD)%unit_no, Fmt=forma, Advance='No', IOStat=keyio) &
                  xyz(0), xyz(1), xyz(2), xyz(3)
              Else
                Read (Unit=files(FILE_REVOLD)%unit_no, IOStat=keyio) &
                  xyz(0), xyz(1), xyz(2), xyz(3)
              End If
            End If
            If (keyio /= 0) i_tmp = 1

            Call gbcast(comm, xyz, 0)
            gidx = Nint(xyz(0))

            ! assign particle velocities to the corresponding domains

            Do i = 1, config%natms
              If (config%ltg(i) == gidx) Then
                green%vxi(i, j) = xyz(1)
                green%vyi(i, j) = xyz(2)
                green%vzi(i, j) = xyz(3)
              End If
            End Do
          End Do
        End Do

        ! If 'restart' is impossible go to 'restart noscale' and reinitialise

        Call gbcast(comm, i_tmp, 0)
        If (i_tmp /= 0) Then
          If (comm%idnode == 0) Then
            Call warning(190, 0.0_wp, 0.0_wp, 0.0_wp)
            Call files(FILE_REVOLD)%close ()
          End If
          keyres = RESTART_KEY_NOSCALE
          Go To 50
        End If

      End If

      If (comm%idnode == 0) Call files(FILE_REVOLD)%close ()

    Else

      ! force force and stats%stress recalculation when 'restart' is not on

      config%levcfg = 1

    End If

    ! number densities needed for long-range corrections

    ! evaluate species populations in system (separate totals for non-frozen atoms)

    Do i = 1, config%natms
      k = config%ltype(i)
      sites%num_type(k) = sites%num_type(k) + 1.0_wp
      If (config%lfrzn(i) == 0) sites%num_type_nf(k) = sites%num_type_nf(k) + 1.0_wp
    End Do

    ! global number densities

    Call gsum(comm, sites%num_type(1:sites%ntype_atom))
    Call gsum(comm, sites%num_type_nf(1:sites%ntype_atom))

    ! number densities

    Do i = 1, sites%ntype_atom
      If (sites%num_type(i) > zero_plus) sites%dens(i) = sites%num_type(i) / config%volm
    End Do

    ! Get long-range corrections

    ! vdws%elrc & vdws%vlrc arrays are zeroed in vdws,
    ! no lrc when vdw interactions are force-shifted

    If (vdws%n_vdw > 0 .and. (.not. vdws%l_force_shift)) Call vdw_lrc(sites, vdws, config, comm)

    ! met%elrc & met%vlrc arrays are zeroed in metal_module

    If (met%n_potentials > 0) Call metal_lrc(met, sites, config, comm)

  End Subroutine system_init

  Subroutine system_expand(l_str, rcut, io, cshell, cons, bond, angle, &
                           dihedral, inversion, sites, netcdf, rigid, config, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 utility to expand the MD system by a config%nx*config%ny*config%nz volumetric
    ! replication of its contents along the MD cell lattice vectors,
    ! creating a new matching pair of topology-interaction (FIELD) and
    ! crystallographic (CONFIG) files, preserving FIELD's template intact
    !
    ! supported image conditions: 1,2,3, 6(config%nz==1)
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! contrib   - w.smith, i.j.bush
    ! contrib   - a.m.elena february 2017
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,                  Intent(In   ) :: l_str
    Real(Kind=wp),            Intent(In   ) :: rcut
    Type(io_type),            Intent(InOut) :: io
    Type(core_shell_type),    Intent(In   ) :: cshell
    Type(constraints_type),   Intent(In   ) :: cons
    Type(bonds_type),         Intent(In   ) :: bond
    Type(angles_type),        Intent(In   ) :: angle
    Type(dihedrals_type),     Intent(InOut) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(site_type),          Intent(In   ) :: sites
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(rigid_bodies_type),  Intent(In   ) :: rigid
    Type(configuration_type), Intent(InOut) :: config
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm

    Character                                                    :: lf
    Character(Len=200)                                           :: record, record1
    Character(len=256)                                           :: message, messages(5)
    Character(Len=40)                                            :: fcfg, ffld, fmpl, word
    Real(Kind=wp), Allocatable, Dimension(:)                     :: f1, f2, f3, f4, f5, f6, f7, &
                                                                    f8, f9, x_scaled, xm, &
                                                                    y_scaled, ym, z_scaled, zm
    Real(Kind=wp)                                                :: angles(1:3), c1, c2, c3, c4, &
                                                                    cell_vecs(1:3, 1:3), &
                                                                    celprp(1:10), fx, fy, fz, hwx, &
                                                                    hwy, hwz, lengths(1:3), r, t, &
                                                                    x, x1(1:1), y, y1(1:1), z, &
                                                                    z1(1:1)
    Logical                                                      :: lmpldt, safe, safeg, safel, &
                                                                    safem, safer, safex, safey, &
                                                                    safez
    Integer, Allocatable, Dimension(:, :, :)                     :: i_xyz
    Integer, Allocatable, Dimension(:)                           :: ltg_scaled
    Integer(Kind=offset_kind)                                    :: top_skip
    Integer(Kind=li)                                             :: offset, rec
    Integer                                                      :: at_scaled, fail(1:5), fh, i, &
                                                                    iang, iatm, ibond, icnst, &
                                                                    idih, idm, ierr, iinv, imols, &
                                                                    indatm, indatm1, index, &
                                                                    io_write, irgd, ishls, itmols, &
                                                                    ix, iy, iz, j, jatm, loc_ind, &
                                                                    lrgd, m, mxiter, nall, nangle, &
                                                                    nattot, nbonds, nconst, &
                                                                    ndihed, ninver, nrigid, &
                                                                    nshels, sapmpt, sapmtt, setspc
    Character(Len=recsz)                                         :: record2, record3
    Character(Len=Len(config%atmnam)), Allocatable, Dimension(:) :: atmnam_scaled

! Some parameters and variables needed by io interfaces

    lmpldt = .false.

    fail = 0
    Allocate (f1(1:config%nx), f2(1:config%nx), f3(1:config%nx), Stat=fail(1))
    Allocate (f4(1:config%ny), f5(1:config%ny), f6(1:config%ny), Stat=fail(2))
    Allocate (f7(1:config%nz), f8(1:config%nz), f9(1:config%nz), Stat=fail(3))
    Allocate (i_xyz(1:config%nx, 1:config%ny, 1:config%nz), Stat=fail(4))
    Allocate (xm(1:10 * config%mxatms), ym(1:10 * config%mxatms), zm(1:10 * config%mxatms), Stat=fail(5))

    If (Any(fail > 0)) Then
      Write (message, '(a)') 'system_expand allocation failure'
      Call error(0, message)
    End If

    ! Get write buffer size and line feed character

    Call io_get_parameters(io, user_method_write=io_write)
    Call io_get_parameters(io, user_line_feed=lf)

    ! Print elapsed time and option header

    Call gtime(t)
    Write (message, '(a,f12.3,a)') 'time elapsed since job start: ', t, ' sec'
    Call info(message, .true.)
    Write (messages(1), '(a)') '*** Expanding the MD system by a nx*ny*nz volumetric replication        ***'
    Write (messages(2), '(a)') '*** of its contents along the MD cell lattice vectors, creating         ***'
    Write (messages(3), '(a)') '*** a new matching pair of topology-interaction (FIELD) and             ***'
    Write (messages(4), '(a)') "*** crystallographic (CONFIG) files, preserving FIELD's template intact ***"
    Write (messages(5), '(a,3i5)') '*** Replication dimensions (nx,ny,nz):', config%nx, config%ny, config%nz
    Call info(messages, 5, .true.)

    ! Holt or change execution if imcon is unsupported

    If (config%imcon == IMCON_NOPBC) Call error(570)
    If (config%imcon == IMCON_SLAB .and. config%nz > 1) Then
      config%nz = 1
      Call warning(350, 0.0_wp, 0.0_wp, 0.0_wp)
      Write (message, '(a,3i5)') '*** Replication dimensions (nx,ny,nz):', config%nx, config%ny, config%nz
      Call info(message)
    End If

    ! Create names for the expanded CONFIG and FIELD

    record = ' '; Write (record, '(3(a1,i0))') '_', config%nx, '_', config%ny, '_', config%nz
    fcfg = ' '
    fcfg = Trim(files(FILE_CONFIG)%filename)//record(1:Len_trim(record))
    ffld = ' '
    ffld = Trim(files(FILE_FIELD)%filename)//record(1:Len_trim(record))
    fmpl = ' '
    fmpl = "MPOLES"//record(1:Len_trim(record))

    ! netCDF CONFIG name convention

    If (io_write == IO_WRITE_SORTED_NETCDF) fcfg = fcfg(1:Len_trim(fcfg))//'.nc'

    fx = Real(config%nx, wp)
    fy = Real(config%ny, wp)
    fz = Real(config%nz, wp)

    nall = config%nx * config%ny * config%nz

    ! Define cell vector displacement in z direction

    Do iz = 1, config%nz
      z = Real(2 * iz - config%nz - 1, wp)
      f7(iz) = config%cell(7) * z / 2.0_wp
      f8(iz) = config%cell(8) * z / 2.0_wp
      f9(iz) = config%cell(9) * z / 2.0_wp
    End Do

    ! Define cell vector displacement in y direction

    Do iy = 1, config%ny
      y = Real(2 * iy - config%ny - 1, wp)
      f4(iy) = config%cell(4) * y / 2.0_wp
      f5(iy) = config%cell(5) * y / 2.0_wp
      f6(iy) = config%cell(6) * y / 2.0_wp
    End Do

    ! Define cell vector displacement in x direction

    Do ix = 1, config%nx
      x = Real(2 * ix - config%nx - 1, wp)
      f1(ix) = config%cell(1) * x / 2.0_wp
      f2(ix) = config%cell(2) * x / 2.0_wp
      f3(ix) = config%cell(3) * x / 2.0_wp
    End Do

    ! Define hypercube counter

    Do iz = 1, config%nz
      Do iy = 1, config%ny
        Do ix = 1, config%nx
          i_xyz(ix, iy, iz) = (ix - 1) + config%nx * ((iy - 1) + config%ny * (iz - 1))
        End Do
      End Do
    End Do

    Call dcell(config%cell, celprp) ! get config%cell properties

    ! define half cell widths and bond-length limit

    hwx = celprp(7) / 2.0_wp
    hwy = celprp(8) / 2.0_wp
    hwz = celprp(9) / 2.0_wp
    c1 = Min(rcut / 2.0_wp, 1.75_wp)
    c2 = c1 * 4.0_wp / 3.0_wp
    c3 = c2 * 4.0_wp / 3.0_wp
    c4 = c3 * 4.0_wp / 3.0_wp

    If (comm%idnode == 0) Then

      ! Make sure CONFIG(new) is empty and open it

      If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
          io_write == IO_WRITE_UNSORTED_DIRECT .or. &
          io_write == IO_WRITE_SORTED_MPIIO .or. &
          io_write == IO_WRITE_SORTED_DIRECT .or. &
          io_write == IO_WRITE_SORTED_NETCDF) Then

        Call io_set_parameters(io, user_comm=comm_self)
        Call io_init(io, recsz)
        Call io_delete(io, fcfg(1:Len_trim(fcfg)), comm)
        If (io_write == IO_WRITE_SORTED_NETCDF) Then
          Call io_nc_create(netcdf, comm_self, fcfg(1:Len_trim(fcfg)), config%cfgname, config%megatm * nall)
        End If
        Call io_open(io, io_write, comm_self, fcfg(1:Len_trim(fcfg)), mode_wronly + mode_create, fh)

      Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
               io_write == IO_WRITE_SORTED_MASTER) Then

        Open (Newunit=files(FILE_CONFIG)%unit_no, File=fcfg(1:Len_trim(fcfg)), Status='replace')
        Call files(FILE_CONFIG)%close ()
        Open (Newunit=files(FILE_CONFIG)%unit_no, File=fcfg(1:Len_trim(fcfg)), Form='formatted', Access='direct', Recl=recsz)
      End If

      ! Write configuration file headers

      Write (message, '(2a)') '*** Expanding CONFIG in file ', fcfg(1:Len_trim(fcfg))
      Call info(message, .true.)

      ! check if we expand a cube since not all time we end up with a cube back
      If ((config%imcon == IMCON_CUBIC) .and. &
          ((config%nx /= config%ny) .or. (config%nx /= config%nz) .or. (config%ny /= config%nz))) Then
        config%imcon = IMCON_PARALLELOPIPED
      End If

      If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
          io_write == IO_WRITE_UNSORTED_DIRECT .or. &
          io_write == IO_WRITE_SORTED_MPIIO .or. &
          io_write == IO_WRITE_SORTED_DIRECT) Then

        Write (record2, Fmt='(a72,a1)') config%cfgname(1:72), lf
        Call io_write_record(io, fh, Int(0, offset_kind), record2)

        Write (record2, Fmt='(3i10,a42,a1)') 0, config%imcon, nall * config%megatm, Repeat(' ', 42), lf
        Call io_write_record(io, fh, Int(1, offset_kind), record2)

        Write (record2, Fmt='(3f20.10,a12,a1)') fx * config%cell(1), fx * config%cell(2), fx * config%cell(3), Repeat(' ', 12), lf
        Call io_write_record(io, fh, Int(2, offset_kind), record2)

        Write (record2, Fmt='(3f20.10,a12,a1)') fy * config%cell(4), fy * config%cell(5), fy * config%cell(6), Repeat(' ', 12), lf
        Call io_write_record(io, fh, Int(3, offset_kind), record2)

        Write (record2, Fmt='(3f20.10,a12,a1)') fz * config%cell(7), fz * config%cell(8), fz * config%cell(9), Repeat(' ', 12), lf
        Call io_write_record(io, fh, Int(4, offset_kind), record2)

      Else If (io_write == IO_WRITE_SORTED_NETCDF) Then

        i = 1 ! For config there is only one frame

        Call io_nc_put_var(io, 'time', fh, 0.0_wp, i, 1)
        Call io_nc_put_var(io, 'step', fh, 0, i, 1)
        Call io_nc_put_var(io, 'datalevel', fh, 0, i, 1)
        Call io_nc_put_var(io, 'imageconvention', fh, config%imcon, i, 1)

        cell_vecs(:, 1) = fx * config%cell(1:3)
        cell_vecs(:, 2) = fy * config%cell(4:6)
        cell_vecs(:, 3) = fz * config%cell(7:9)

        lengths(1) = fx * celprp(1)
        lengths(2) = fy * celprp(2)
        lengths(3) = fz * celprp(3)

        angles(1) = Acos(celprp(5))
        angles(2) = Acos(celprp(6))
        angles(3) = Acos(celprp(4))
        angles = angles * 180.0_wp / (4.0_wp * Atan(1.0_wp)) ! Convert to degrees

        Call io_nc_put_var(io, 'cell', fh, cell_vecs, (/1, 1, i/), (/3, 3, 1/))
        Call io_nc_put_var(io, 'cell_lengths', fh, lengths, (/1, i/), (/3, 1/))
        Call io_nc_put_var(io, 'cell_angles', fh, angles, (/1, i/), (/3, 1/))

      Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
               io_write == IO_WRITE_SORTED_MASTER) Then

        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a72,a1)', Rec=Int(1, li)) &
          config%cfgname(1:72), lf
        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(3i10,a42,a1)', Rec=Int(2, li)) &
          0, config%imcon, nall * config%megatm, Repeat(' ', 42), lf
        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(3f20.12,a12,a1)', Rec=Int(3, li)) &
          fx * config%cell(1), fx * config%cell(2), fx * config%cell(3), Repeat(' ', 12), lf
        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(3f20.12,a12,a1)', Rec=Int(4, li)) &
          fy * config%cell(4), fy * config%cell(5), fy * config%cell(6), Repeat(' ', 12), lf
        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(3f20.12,a12,a1)', Rec=Int(5, li)) &
          fz * config%cell(7), fz * config%cell(8), fz * config%cell(9), Repeat(' ', 12), lf

      End If

    End If

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_SORTED_MPIIO .or. &
        io_write == IO_WRITE_SORTED_NETCDF) Then

      If (comm%idnode == 0) Then
        Call io_close(io, fh)
        Call io_finalize(io)
      End If

      Allocate (atmnam_scaled(1:config%natms * nall), ltg_scaled(1:config%natms * nall), Stat=fail(1))
      Allocate (x_scaled(1:config%natms * nall), y_scaled(1:config%natms * nall), z_scaled(1:config%natms * nall), &
                Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'system_expand allocation failure 0 '
        Call error(0, message)
      End If

    End If

    ! Line counter in CONFIG(new):levcfg=0

    offset = Int(5, li)

    ! Global atom counter

    nattot = 0

    ! Local atom counter

    indatm = 1

    ! local counter in scaled system

    at_scaled = 0

    ! running site and topology related intra-indices

    nshels = 0
    nconst = 0
    nrigid = 0
    nbonds = 0
    nangle = 0
    ndihed = 0
    ninver = 0

    Write (message, '(a)') 'Checking topological contiguity of molecules...'
    Call info(message, .true.)

    safeg = .true. ! topology presumed safe

    sapmpt = 0
    Do itmols = 1, sites%ntype_mol
      setspc = sites%num_mols(itmols) * sites%num_site(itmols)

      sapmtt = 0
      Do imols = 1, sites%num_mols(itmols)
        If (sites%num_site(itmols) > 10 * config%mxatms) Call error(0, message)

        ! Grab the coordinates of the atoms constituting this molecule

        indatm1 = indatm
        Do m = 1, sites%num_site(itmols)
          nattot = nattot + 1 ! Increase global atom counter in CONFIG(old)

          If (config%lsa(indatm1) == nattot) Then ! If a local atom has a global index nattot
            loc_ind = config%lsi(indatm1)
            xm(m) = config%parts(loc_ind)%xxx
            ym(m) = config%parts(loc_ind)%yyy
            zm(m) = config%parts(loc_ind)%zzz
            indatm1 = indatm1 + 1 ! Increase local atom counter
          Else
            xm(m) = 0.0_wp
            ym(m) = 0.0_wp
            zm(m) = 0.0_wp
          End If
        End Do
        nattot = nattot - sites%num_site(itmols)

        Call gsum(comm, xm(1:sites%num_site(itmols)))
        Call gsum(comm, ym(1:sites%num_site(itmols)))
        Call gsum(comm, zm(1:sites%num_site(itmols)))

        ! Start unwrapping - not safe at start for each molecule

        indatm1 = nattot - sapmpt - sapmtt
        safe = .false.; mxiter = 0
        Do While ((.not. safe) .and. mxiter < 42) ! meaning of LUEE is the limit
          If (.not. safe) mxiter = mxiter + 1

          If ((mxiter == 42 .and. (.not. safe)) .and. l_str) Then
            Write (message, '(2(a,i10))') 'molecular type #: ', itmols, ' molecule #: ', imols
            Call info(message, .true.)
          End If

          safe = .true.

          safel = .true.
          Do ishls = 1, cshell%numshl(itmols)
            nshels = nshels + 1

            iatm = cshell%lstshl(1, nshels) - indatm1
            jatm = cshell%lstshl(2, nshels) - indatm1

            safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
            safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
            safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
            safer = (safex .and. safey .and. safez)
            If (.not. safer) Then
              x1(1) = xm(jatm) - xm(iatm)
              y1(1) = ym(jatm) - ym(iatm)
              z1(1) = zm(jatm) - zm(iatm)
              Call images(config%imcon, config%cell, 1, x1, y1, z1)
              xm(jatm) = x1(1) + xm(iatm)
              ym(jatm) = y1(1) + ym(iatm)
              zm(jatm) = z1(1) + zm(iatm)
              safer = .true.
            End If
            x = Abs(xm(jatm) - xm(iatm))
            y = Abs(ym(jatm) - ym(iatm))
            z = Abs(zm(jatm) - zm(iatm))
            t = c1
            safex = (x < t)
            safey = (y < t)
            safez = (z < t)
            r = Sqrt(x**2 + y**2 + z**2)
            If (safex .and. safey .and. safez) Then
              safer = (r < t)
            Else
              safer = .false.
            End If
            If ((mxiter == 42 .and. (.not. safer)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (message, '(a,2(f7.2,a))') &
                'possible distance violation: ', r, ' > ', t, ' Angstroms'
              Call info(message, .true.)

              t = c2
              If (r > t) Then
                Write (message, Fmt='(a,f7.2,a)') 'cutoff: ', t, ' Angstroms'
                Call warning(message, .true.)
              End If

              Write (messages(1), '(a,3i10)') &
                'core_shell unit #(local) -> m. type # -> molecule #:', ishls, itmols, imols
              Write (messages(2), '(a,3(1x,l1))') &
                'member :: global index :: x ::      y ::      z', safex, safey, safez
              Write (messages(3), '(a,i10,3f10.1)') &
                'core  ', nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
              Write (messages(4), '(a,i10,3f10.1)') &
                'shell ', nattot + jatm, xm(jatm), ym(jatm), zm(jatm)
              Call info(messages, 4, .true.)
            End If
            safel = (safel .and. safer)
          End Do
          safe = (safe .and. safel)

          safel = .true.
          Do icnst = 1, cons%numcon(itmols)
            nconst = nconst + 1

            iatm = cons%lstcon(1, nconst) - indatm1
            jatm = cons%lstcon(2, nconst) - indatm1

            safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
            safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
            safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
            safer = (safex .and. safey .and. safez)
            If (.not. safer) Then
              x1(1) = xm(jatm) - xm(iatm)
              y1(1) = ym(jatm) - ym(iatm)
              z1(1) = zm(jatm) - zm(iatm)
              Call images(config%imcon, config%cell, 1, x1, y1, z1)
              xm(jatm) = x1(1) + xm(iatm)
              ym(jatm) = y1(1) + ym(iatm)
              zm(jatm) = z1(1) + zm(iatm)
              safer = .true.
            End If
            x = Abs(xm(jatm) - xm(iatm))
            y = Abs(ym(jatm) - ym(iatm))
            z = Abs(zm(jatm) - zm(iatm))
            t = c2
            safex = (x < t)
            safey = (y < t)
            safez = (z < t)
            r = Sqrt(x**2 + y**2 + z**2)
            If (safex .and. safey .and. safez) Then
              safer = (r < t)
            Else
              safer = .false.
            End If
            If ((mxiter == 42 .and. (.not. safer)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (message, '(a,2(f7.2,a))') &
                'possible distance violation: ', r, ' > ', t, ' Angstroms'
              Call info(message, .true.)

              t = c3
              If (r > t) Then
                Write (message, Fmt='(a,f7.2,a)') 'cutoff: ', t, ' Angstroms'
                Call warning(message, .true.)
              End If

              Write (messages(1), '(a,3i10)') &
                'constraint unit #(local) -> m. type # -> molecule #:', icnst, itmols, imols
              Write (messages(2), '(a,3(1x,l1))') &
                'member :: global index :: x ::      y ::      z', safex, safey, safez
              Write (messages(3), '(2i10,3f10.1)') &
                1, nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
              Write (messages(4), '(2i10,3f10.1)') &
                2, nattot + jatm, xm(jatm), ym(jatm), zm(jatm)
              Call info(messages, 4, .true.)
            End If
            safel = (safel .and. safer)
          End Do
          safe = (safe .and. safel)

          safel = .true.
          Do irgd = 1, rigid%num(itmols)
            nrigid = nrigid + 1

            safem = .true.
            lrgd = rigid%lst(0, nrigid)
            Do i = 1, lrgd - 1
              iatm = rigid%lst(i, nrigid) - indatm1
              Do j = i + 1, lrgd
                jatm = rigid%lst(j, nrigid) - indatm1

                safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
                safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
                safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
                safer = (safex .and. safey .and. safez)
                If (.not. safer) Then
                  x1(1) = xm(jatm) - xm(iatm)
                  y1(1) = ym(jatm) - ym(iatm)
                  z1(1) = zm(jatm) - zm(iatm)
                  Call images(config%imcon, config%cell, 1, x1, y1, z1)
                  xm(jatm) = x1(1) + xm(iatm)
                  ym(jatm) = y1(1) + ym(iatm)
                  zm(jatm) = z1(1) + zm(iatm)
                  safer = .true.
                End If
                x = Abs(xm(jatm) - xm(iatm))
                y = Abs(ym(jatm) - ym(iatm))
                z = Abs(zm(jatm) - zm(iatm))
                t = rcut
                safex = (x < t)
                safey = (y < t)
                safez = (z < t)
                r = Sqrt(x**2 + y**2 + z**2)
                If (safex .and. safey .and. safez) Then
                  safer = (r < t)
                Else
                  safer = .false.
                End If
                If ((mxiter == 42 .and. (.not. safer)) .and. l_str) Then
                  Write (message, '(a,2i10,2(f7.2,a))') &
                    'distance violation: ', i, j, r, ' > ', t, ' Angstroms'
                  Call warning(message, .true.)
                End If
                safem = (safem .and. safer)
              End Do
            End Do

            If ((mxiter == 42 .and. (.not. safem)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (messages(1), '(a,3i10)') &
                'rigid body unit #(local) -> m. type # -> molecule #:', irgd, itmols, imols
              Write (messages(2), '(a)') &
                'member :: global index :: x ::      y ::      z'
              Call info(messages, 2, .true.)

              Do i = 1, lrgd
                iatm = rigid%lst(i, nrigid) - indatm1
                Write (message, '(2i10,3f10.1)') i, nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
                Call info(message, .true.)
              End Do
            End If
            safel = (safel .and. safem)
          End Do
          safe = (safe .and. safel)

          safel = .true.
          Do ibond = 1, bond%num(itmols)
            nbonds = nbonds + 1

            iatm = bond%lst(1, nbonds) - indatm1
            jatm = bond%lst(2, nbonds) - indatm1

            safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
            safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
            safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
            safer = (safex .and. safey .and. safez)
            If (.not. safer) Then
              x1(1) = xm(jatm) - xm(iatm)
              y1(1) = ym(jatm) - ym(iatm)
              z1(1) = zm(jatm) - zm(iatm)
              Call images(config%imcon, config%cell, 1, x1, y1, z1)
              xm(jatm) = x1(1) + xm(iatm)
              ym(jatm) = y1(1) + ym(iatm)
              zm(jatm) = z1(1) + zm(iatm)
              safer = .true.
            End If
            x = Abs(xm(jatm) - xm(iatm))
            y = Abs(ym(jatm) - ym(iatm))
            z = Abs(zm(jatm) - zm(iatm))
            If (.not. bond%restrained(nbonds)) Then
              t = c2
            Else
              t = 3.0_wp * c1
            End If
            safex = (x < t)
            safey = (y < t)
            safez = (z < t)
            r = Sqrt(x**2 + y**2 + z**2)
            If (safex .and. safey .and. safez) Then
              safer = (r < t)
            Else
              safer = .false.
            End If
            If ((mxiter == 42 .and. (.not. safer)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (message, '(a,2(f7.2,a))') &
                'possible distance violation: ', r, ' > ', t, ' Angstroms'
              Call info(message, .true.)

              If (.not. bond%restrained(nbonds)) Then
                t = c3
              Else
                t = 3.0_wp * c2
              End If
              If (r > t) Then
                Write (message, Fmt='(a,f7.2,a)') 'cutoff: ', t, ' Angstroms'
                Call warning(message, .true.)
              End If

              Write (messages(1), '(a,3i10)') &
                'bond unit #(local) -> m. type # -> molecule #:', ibond, itmols, imols
              Write (messages(2), '(a,3(1x,l1))') &
                'member :: global index :: x ::      y ::      z', safex, safey, safez
              Write (messages(3), '(2i10,3f10.1)') &
                1, nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
              Write (messages(4), '(2i10,3f10.1)') &
                2, nattot + jatm, xm(jatm), ym(jatm), zm(jatm)
              Call info(messages, 4, .true.)
            End If
            safel = (safel .and. safer)
          End Do
          safe = (safe .and. safel)

          safel = .true.
          Do iang = 1, angle%num(itmols)
            nangle = nangle + 1

            safem = .true.
            Do i = 1, 2
              iatm = angle%lst(i, nangle) - indatm1
              Do j = i + 1, 3
                jatm = angle%lst(j, nangle) - indatm1

                safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
                safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
                safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
                safer = (safex .and. safey .and. safez)
                If (.not. safer) Then
                  x1(1) = xm(jatm) - xm(iatm)
                  y1(1) = ym(jatm) - ym(iatm)
                  z1(1) = zm(jatm) - zm(iatm)
                  Call images(config%imcon, config%cell, 1, x1, y1, z1)
                  xm(jatm) = x1(1) + xm(iatm)
                  ym(jatm) = y1(1) + ym(iatm)
                  zm(jatm) = z1(1) + zm(iatm)
                  safer = .true.
                End If
                x = Abs(xm(jatm) - xm(iatm))
                y = Abs(ym(jatm) - ym(iatm))
                z = Abs(zm(jatm) - zm(iatm))
                If (.not. angle%restrained(nangle)) Then
                  t = c1 * Real(j - i + 1, wp)
                Else
                  t = c3 * Real(j - i + 1, wp)
                End If
                safex = (x < t)
                safey = (y < t)
                safez = (z < t)
                r = Sqrt(x**2 + y**2 + z**2)
                If (safex .and. safey .and. safez) Then
                  safer = (r < t)
                Else
                  safer = .false.
                End If
                If ((mxiter == 42 .and. (.not. safer)) .and. (l_str .and. comm%idnode == 0)) Then
                  Write (message, '(a,2i10,2(f7.2,a))') &
                    'possible distance violation: ', i, j, r, ' > ', t, ' Angstroms'
                  Call info(message, .true.)

                  If (.not. angle%restrained(nangle)) Then
                    t = c2 * Real(j - i + 1, wp)
                  Else
                    t = c4 * Real(j - i + 1, wp)
                  End If
                  If (r > t) Then
                    Write (message, '(a,f7.2,a)') 'cutoff: ', t, ' Angstroms'
                    Call warning(message, .true.)
                  End If
                End If
                safem = (safem .and. safer)
              End Do
            End Do

            If ((mxiter == 42 .and. (.not. safem)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (messages(1), '(a,3i10)') &
                'angle unit #(local) -> m. type # -> molecule #:', iang, itmols, imols
              Write (messages(2), '(a)') &
                'member :: global index :: x ::      y ::      z'
              Call info(messages, 2, .true.)

              Do i = 1, 3
                iatm = angle%lst(i, nangle) - indatm1
                Write (message, '(2i10,3f10.1)') i, nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
                Call info(message, .true.)
              End Do
            End If
            safel = (safel .and. safem)
          End Do
          safe = (safe .and. safel)

          safel = .true.
          Do idih = 1, dihedral%num(itmols)
            ndihed = ndihed + 1

            safem = .true.
            Do i = 1, 3
              iatm = dihedral%lst(i, ndihed) - indatm1
              Do j = i + 1, 4
                jatm = dihedral%lst(j, ndihed) - indatm1

                safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
                safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
                safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
                safer = (safex .and. safey .and. safez)
                If (.not. safer) Then
                  x1(1) = xm(jatm) - xm(iatm)
                  y1(1) = ym(jatm) - ym(iatm)
                  z1(1) = zm(jatm) - zm(iatm)
                  Call images(config%imcon, config%cell, 1, x1, y1, z1)
                  xm(jatm) = x1(1) + xm(iatm)
                  ym(jatm) = y1(1) + ym(iatm)
                  zm(jatm) = z1(1) + zm(iatm)
                  safer = .true.
                End If
                t = (c1 + c2) * Real(j - i + 1, wp) / 2.0_wp
                x = Abs(xm(jatm) - xm(iatm))
                y = Abs(ym(jatm) - ym(iatm))
                z = Abs(zm(jatm) - zm(iatm))
                safex = (x < t)
                safey = (y < t)
                safez = (z < t)
                r = Sqrt(x**2 + y**2 + z**2)
                If (safex .and. safey .and. safez) Then
                  safer = (r < t)
                Else
                  safer = .false.
                End If
                If ((mxiter == 42 .and. (.not. safer)) .and. (l_str .and. comm%idnode == 0)) Then
                  Write (message, '(a,2i10,2(f7.2,a))') &
                    'possible distance violation: ', i, j, r, ' > ', t, ' Angstroms'
                  Call info(message, .true.)

                  t = (c2 + c3) * Real(j - i + 1, wp) / 2.0_wp
                  If (r > t) Then
                    Write (message, Fmt='(a,f7.2,a)') 'cutoff: ', t, ' Angstroms'
                    Call warning(message, .true.)
                  End If
                End If
                safem = (safem .and. safer)
              End Do
            End Do

            If ((mxiter == 42 .and. (.not. safem)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (messages(1), '(a,3i10)') &
                'dihedral unit #(local) -> m. type # -> molecule #:', idih, itmols, imols
              Write (messages(2), '(a)') &
                'member :: global index :: x ::      y ::      z'
              Call info(messages, 2, .true.)

              Do i = 1, 4
                iatm = dihedral%lst(i, ndihed) - indatm1
                Write (message, '(2i10,3f10.1)') i, nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
                Call info(message, .true.)
              End Do
            End If
            safel = (safel .and. safem)
          End Do
          safe = (safe .and. safel)

          safel = .true.
          Do iinv = 1, inversion%num(itmols)
            ninver = ninver + 1

            safem = .true.
            Do i = 1, 3
              iatm = inversion%lst(i, ninver) - indatm1
              Do j = i + 1, 4
                jatm = inversion%lst(j, ninver) - indatm1

                safex = (Abs(xm(jatm) - xm(iatm)) < hwx)
                safey = (Abs(ym(jatm) - ym(iatm)) < hwy)
                safez = (Abs(zm(jatm) - zm(iatm)) < hwz)
                safer = (safex .and. safey .and. safez)
                If (.not. safer) Then
                  x1(1) = xm(jatm) - xm(iatm)
                  y1(1) = ym(jatm) - ym(iatm)
                  z1(1) = zm(jatm) - zm(iatm)
                  Call images(config%imcon, config%cell, 1, x1, y1, z1)
                  xm(jatm) = x1(1) + xm(iatm)
                  ym(jatm) = y1(1) + ym(iatm)
                  zm(jatm) = z1(1) + zm(iatm)
                  safer = .true.
                End If
                x = Abs(xm(jatm) - xm(iatm))
                y = Abs(ym(jatm) - ym(iatm))
                z = Abs(zm(jatm) - zm(iatm))
                t = c2
                safex = (x < t)
                safey = (y < t)
                safez = (z < t)
                r = Sqrt(x**2 + y**2 + z**2)
                If (safex .and. safey .and. safez) Then
                  safer = (r < t)
                Else
                  safer = .false.
                End If
                If ((mxiter == 42 .and. (.not. safer)) .and. (l_str .and. comm%idnode == 0)) Then
                  Write (message, '(a,2i10,2(f7.2,a))') &
                    'possible distance violation: ', i, j, r, ' > ', t, ' Angstroms'
                  Call info(message, .true.)

                  t = c3
                  If (r > t) Then
                    Write (message, Fmt='(a,f7.2,a)') 'cutoff: ', t, ' Angstroms'
                    Call warning(message, .true.)
                  End If
                End If

                safem = (safem .and. safer)
              End Do
            End Do

            If ((mxiter == 42 .and. (.not. safem)) .and. (l_str .and. comm%idnode == 0)) Then
              Write (messages(1), '(a,3i10)') &
                'inversion unit #(local) -> m. type # -> molecule #:', iinv, itmols, imols
              Write (messages(2), '(a)') &
                'member :: global index :: x ::      y ::      z'
              Call info(messages, 2, .true.)

              Do i = 1, 4
                iatm = inversion%lst(i, ninver) - indatm1
                Write (message, '(2i10,3f10.1)') i, nattot + iatm, xm(iatm), ym(iatm), zm(iatm)
                Call info(message, .true.)
              End Do
            End If
            safel = (safel .and. safem)
          End Do
          safe = (safe .and. safel)

          If (((.not. safe) .and. imols <= sites%num_mols(itmols) .and. mxiter < 42) .or. &
              imols < sites%num_mols(itmols)) Then
            nshels = nshels - cshell%numshl(itmols)
            nconst = nconst - cons%numcon(itmols)
            nrigid = nrigid - rigid%num(itmols)
            nbonds = nbonds - bond%num(itmols)
            nangle = nangle - angle%num(itmols)
            ndihed = ndihed - dihedral%num(itmols)
            ninver = ninver - inversion%num(itmols)
          End If
        End Do

        safeg = (safeg .and. safe)

        Do m = 1, sites%num_site(itmols)
          nattot = nattot + 1 ! Increase global atom counter in CONFIG(old)

          If (config%lsa(indatm) == nattot) Then ! If a local atom has a global index nattot

            ! Determine sending node for UN/SORTED MASTER

            If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                io_write == IO_WRITE_SORTED_MASTER) Then
              idm = comm%idnode
              Call gsum(comm, idm)
            End If

            ! Get the local index of the particle

            loc_ind = config%lsi(indatm)

            ! Do particle replication by vector displacements in cyclic (z,y,x) directions

            Do iz = 1, config%nz
              Do iy = 1, config%ny
                Do ix = 1, config%nx

                  x = xm(m) + f1(ix) + f4(iy) + f7(iz)
                  y = ym(m) + f2(ix) + f5(iy) + f8(iz)
                  z = zm(m) + f3(ix) + f6(iy) + f9(iz)

                  ! Write 2 records @ line 'rec+1' and 'rec+2' for particle 'index' in CONFIG(new)

                  index = i_xyz(ix, iy, iz) * setspc + m
                  rec = offset + Int(2, li) * Int(index, li) - Int(2, li)
                  index = index + Int((offset - Int(5, li)) / Int(2, li))

                  If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
                      io_write == IO_WRITE_SORTED_MPIIO .or. &
                      io_write == IO_WRITE_SORTED_NETCDF) Then

                    at_scaled = at_scaled + 1

                    atmnam_scaled(at_scaled) = config%atmnam(loc_ind)

                    ltg_scaled(at_scaled) = index

                    x_scaled(at_scaled) = x
                    y_scaled(at_scaled) = y
                    z_scaled(at_scaled) = z

                  Else

                    Write (record2, Fmt='(a8,i10,a54,a1)') config%atmnam(loc_ind), index, Repeat(' ', 54), lf
                    Write (record3, Fmt='(3g20.12,a12,a1)') x, y, z, Repeat(' ', 12), lf

                    If (io_write == IO_WRITE_UNSORTED_DIRECT .or. &
                        io_write == IO_WRITE_SORTED_DIRECT) Then

                      Call io_write_record(io, fh, Int(rec, offset_kind), record2)
                      rec = rec + Int(1, li)
                      Call io_write_record(io, fh, Int(rec, offset_kind), record3)

                    Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                             io_write == IO_WRITE_SORTED_MASTER) Then

                      If (comm%idnode == 0) Then
                        rec = rec + Int(1, li)
                        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a73)', Rec=rec) record2
                        rec = rec + Int(1, li)
                        Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a73)', Rec=rec) record3
                      Else
                        Call gsend(comm, record2, 0, SysExpand_tag)
                        Call gsend(comm, record3, 0, SysExpand_tag)
                      End If

                    End If

                  End If

                End Do
              End Do
            End Do

            ! Increase local atom counter

            indatm = indatm + 1

          Else

            ! Determine sending node for UN/SORTED MASTER

            If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                io_write == IO_WRITE_SORTED_MASTER) Then
              idm = 0 ! Initialise node number
              Call gsum(comm, idm)

              Do iz = 1, config%nz
                Do iy = 1, config%ny
                  Do ix = 1, config%nx
                    rec = offset + Int(2, li) * (Int(i_xyz(ix, iy, iz), li) * Int(setspc, li) + Int(m, li)) - Int(2, li)

                    If (comm%idnode == 0) Then
                      Call grecv(comm, record2, idm, SysExpand_tag)
                      Call grecv(comm, record3, idm, SysExpand_tag)

                      rec = rec + Int(1, li)
                      Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a73)', Rec=rec) record2
                      rec = rec + Int(1, li)
                      Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a73)', Rec=rec) record3
                    End If
                  End Do
                End Do
              End Do
            End If

          End If
        End Do

        sapmtt = sapmtt + sites%num_site(itmols)
        offset = offset + Int(2, li) * Int(sites%num_site(itmols), li)
      End Do

      sapmpt = sapmpt + sapmtt
      offset = offset + Int(2, li) * Int(nall - 1, li) * Int(setspc, li)
    End Do

    If (.not. safeg) Then
      Call warning('possible topological contiguity failures occurred')
    End If

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_SORTED_MPIIO .or. &
        io_write == IO_WRITE_SORTED_NETCDF) Then

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_init(io, recsz)
      Call io_open(io, io_write, comm%comm, fcfg(1:Len_trim(fcfg)), mode_wronly, fh)

      If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        top_skip = Int(5, offset_kind)
      Else
        top_skip = Int(1, offset_kind) ! netCDF frame
      End If

      Call io_write_sorted_file(io, fh, 0, IO_RESTART, top_skip, at_scaled, &
                                ltg_scaled, atmnam_scaled, &
                                (/0.0_wp/), (/0.0_wp/), (/0.0_wp/), &
                                x_scaled, y_scaled, z_scaled, &
                                (/0.0_wp/), (/0.0_wp/), (/0.0_wp/), &
                                (/0.0_wp/), (/0.0_wp/), (/0.0_wp/), ierr)

      If (ierr /= 0) Then
        Select Case (ierr)
        Case (IO_BASE_COMM_NOT_SET)
          Call error(1050)
        Case (IO_ALLOCATION_ERROR)
          Call error(1053)
        Case (IO_UNKNOWN_WRITE_OPTION)
          Call error(1056)
        Case (IO_UNKNOWN_WRITE_LEVEL)
          Call error(1059)
        End Select
      End If

      Call io_close(io, fh)
      Call io_finalize(io)

      Deallocate (atmnam_scaled, ltg_scaled, Stat=fail(1))
      Deallocate (x_scaled, y_scaled, z_scaled, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'system_expand deallocation failure 0 '
        Call error(0, message)
      End If

    Else If (io_write == IO_WRITE_UNSORTED_DIRECT .or. &
             io_write == IO_WRITE_SORTED_DIRECT) Then

      Call io_close(io, fh)
      Call io_finalize(io)

    Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
             io_write == IO_WRITE_SORTED_MASTER) Then

      Call files(FILE_CONFIG)%close ()

    End If

    Call gtime(t)

    ! Write summary data and proceed with FIELD

    x = 0.5_wp * Min(fx * celprp(7), fy * celprp(8), fz * celprp(9))

    Write (messages(1), '(3a)') '*** ', fcfg(1:Len_trim(fcfg)), ' expansion completed !'
    Write (messages(2), '(a,i10,a)') '*** Size: ', nall * config%megatm, ' particles'
    Write (messages(3), '(a,f10.2,a)') '*** Maximum radius of cutoff: ', x, ' Angstroms'
    Call info(messages, 3, .true.)

    Write (message, '(a,f12.3,a)') 'time elapsed since job start: ', t, ' sec'
    Call info(message, .true.)

    If (comm%idnode == 0) Then
      Open (Newunit=files(FILE_FIELD)%unit_no, File=files(FILE_FIELD)%filename, Status='old')
      Open (Newunit=files(FILE_CONFIG)%unit_no, File=ffld(1:Len_trim(ffld)), Status='replace')
      Write (message, '(2a)') '*** Expanding FIELD in file ', ffld(1:Len_trim(ffld))
      Call info(message, .true.)

      ! omit first line

      record = ' '
      Read (Unit=files(FILE_FIELD)%unit_no, Fmt='(a)', End=10) record
      Call tabs_2_blanks(record); Call strip_blanks(record)
      Write (files(FILE_CONFIG)%unit_no, '(a)') record(1:Len_trim(record))

      ! read and process directives from field file

      Do
        record = ' '
        Read (Unit=files(FILE_FIELD)%unit_no, Fmt='(a)', End=10) record
        Call tabs_2_blanks(record); Call strip_blanks(record)
        record1 = record
        Call get_word(record, word)
        Call lower_case(word)

        If (word(1:5) == 'multi') Then

          ! MPOLES should exist

          lmpldt = .true.

          ! number of molecules of this type

        Else If (word(1:6) == 'nummol') Then

          Call get_word(record, word)
          index = Nint(word_2_real(word))
          Call get_word(record1, word)
          Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a,i10)') word(1:Len_trim(word)), nall * index

          ! close force field file

        Else If (word(1:5) == 'close') Then

          Call gtime(t)

          Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a)') record1(1:Len_trim(record1))
          Call files(FILE_FIELD)%close ()
          Call files(FILE_CONFIG)%close ()
          Write (message, '(3a)') '*** ', ffld(1:Len_trim(ffld)), ' expansion done !'
          Call info(message, .true.)
          Exit

          ! just paste the copy

        Else

          Write (files(FILE_CONFIG)%unit_no, '(a)') record1(1:Len_trim(record1))

        End If
      End Do

      10 Continue

      If (lmpldt) Inquire (File='MPOLES', Exist=lmpldt)
      If (lmpldt) Then
        Open (Unit=nmpldt, File='MPOLES', Status='old')
        Open (Newunit=files(FILE_CONFIG)%unit_no, File=fmpl(1:Len_trim(fmpl)), Status='replace')
        Write (message, '(2a)') '*** Expanding MPOLES in file ', fmpl(1:Len_trim(fmpl))
        Call info(message, .true.)

        ! omit first line

        record = ' '
        Read (Unit=nmpldt, Fmt='(a)', End=10) record
        Call tabs_2_blanks(record); Call strip_blanks(record)
        Write (files(FILE_CONFIG)%unit_no, '(a)') record(1:Len_trim(record))

        ! read and process directives from mpoles file

        Do
          record = ' '
          Read (Unit=nmpldt, Fmt='(a)', End=20) record
          Call tabs_2_blanks(record); Call strip_blanks(record)
          record1 = record
          Call get_word(record, word)
          Call lower_case(word)

          If (word(1:6) == 'nummol') Then

            Call get_word(record, word)
            index = Nint(word_2_real(word))
            Call get_word(record1, word)
            Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a,i10)') word(1:Len_trim(word)), nall * index

            ! close mpoles file

          Else If (word(1:5) == 'close') Then

            Call gtime(t)

            Write (Unit=files(FILE_CONFIG)%unit_no, Fmt='(a)') record1(1:Len_trim(record1))
            Close (Unit=nmpldt)
            Call files(FILE_CONFIG)%close ()
            Write (message, '(3a)') '*** ', fmpl(1:Len_trim(fmpl)), ' expansion done !'
            Call info(message, .true.)
            Exit

            ! just paste the copy

          Else

            Write (files(FILE_CONFIG)%unit_no, '(a)') record1(1:Len_trim(record1))

          End If
        End Do
      End If

      20 Continue

      Write (message, '(a,f12.3,a)') 'time elapsed since job start: ', t, ' sec'
      Call info(message, .true.)
    End If
    Call gsync(comm)

    Deallocate (f1, f2, f3, Stat=fail(1))
    Deallocate (f4, f5, f6, Stat=fail(2))
    Deallocate (f7, f8, f9, Stat=fail(3))
    Deallocate (i_xyz, Stat=fail(4))
    Deallocate (xm, ym, zm, Stat=fail(5))
    If (Any(fail > 0)) Then
      Call info('system_expand dellocation failure ')
      Call error(0, message)
    End If

  End Subroutine system_expand

  Subroutine system_revive(rcut, nstep, time, sites, io, tmst, stats, devel, &
                           green, thermo, bond, angle, dihedral, inversion, zdensity, rdf, netcdf, config, &
                           files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for writing restart files at job termination or
    ! selected intervals in simulation
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith december 1992
    ! amended   - i.t.todorov november 2016
    ! contrib   - m.a.seaton june 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),            Intent(In   ) :: rcut
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(In   ) :: time
    Type(site_type),          Intent(InOut) :: sites
    Type(io_type),            Intent(InOut) :: io
    Real(Kind=wp),            Intent(In   ) :: tmst
    Type(stats_type),         Intent(InOut) :: stats
    Type(development_type),   Intent(In   ) :: devel
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(bonds_type),         Intent(InOut) :: bond
    Type(angles_type),        Intent(InOut) :: angle
    Type(dihedrals_type),     Intent(InOut) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(configuration_type), Intent(InOut) :: config
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Character(Len=40)                        :: forma
    Integer                                  :: fail(1:3), i, j, jatms, jdnode, l, levcfg, nsum
    Integer, Allocatable, Dimension(:)       :: iwrk
    Logical                                  :: ready
    Real(Kind=wp)                            :: r_mxnode
    Real(Kind=wp), Allocatable, Dimension(:) :: axx, ayy, azz, bxx, byy, bzz

    forma = ' '
    fail = 0
    Allocate (iwrk(1:config%mxatms), Stat=fail(1))
    Allocate (axx(1:config%mxatms), ayy(1:config%mxatms), azz(1:config%mxatms), Stat=fail(2))
    Allocate (bxx(1:config%mxatms), byy(1:config%mxatms), bzz(1:config%mxatms), Stat=fail(3))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'system_revive allocation failure '
      Call error(0, message)
    End If

    ! Define format for REVIVE printing in ASCII

    If (devel%l_rout) Then
      i = 64 / 4 - 1 ! Bit_Size(0.0_wp)/4 - 1
      j = Max(stats%mxstak * stats%mxnstk + 1, rdf%max_grid * rdf%max_rdf, rdf%max_grid_usr, config%mxgana * config%mxtana)

      Write (forma, 10) j / 4 + 1, i + 9, i
      10 Format('(1p,', i0, '(/,4e', i0, '.', i0, 'E3))')
    End If

    If (comm%mxnode > 1) Then

      ! globally sum RDF information before saving

      If (rdf%l_collect) Then

        ! maximum rdfs that can be summed in each step

        nsum = config%mxbuff / rdf%max_grid
        If (nsum == 0) Call error(200)

        Do i = 1, rdf%max_rdf, nsum
          Call gsum(comm, rdf%rdf(:, i:Min(i + nsum - 1, rdf%max_rdf)))
        End Do

      End If

      ! globally sum USR RDF information before saving

      If (rdf%max_grid_usr > 0) Call gsum(comm, rdf%usr(1:rdf%max_grid_usr))

      ! globally sum z-density information before saving

      If (zdensity%l_collect) Then

        ! maximum zdensity%density that can be summed in each step

        nsum = config%mxbuff / zdensity%max_grid
        If (nsum == 0) Call error(200)

        Do i = 1, sites%mxatyp, nsum
          Call gsum(comm, zdensity%density(:, i:Min(i + nsum - 1, sites%mxatyp)))
        End Do

      End If

      ! globally sum green%vafdata information before saving

      If (green%samp > 0) Then

        ! maximum green%vafdata that can be summed in each step

        nsum = config%mxbuff / (green%binsize + 1)
        If (nsum == 0) Call error(200)

        Do j = 1, green%samp
          l = (j - 1) * (sites%mxatyp + 1) ! avoid summing up timing information
          Do i = 1, sites%mxatyp, nsum
            Call gsum(comm, green%vafdata(:, l + i:l + Min(i + nsum - 1, sites%mxatyp)))
          End Do
        End Do

      End If

      ! globally sum bonds' distributions information before saving

      If (bond%bin_pdf > 0) Then

        ! maximum bond%dst that can be summed in each step

        nsum = config%mxbuff / (bond%bin_pdf + 1)
        If (nsum == 0) Call error(200)

        Do i = 1, bond%ldf(0), nsum
          Call gsum(comm, bond%dst(:, i:Min(i + nsum - 1, bond%ldf(0))))
        End Do

      End If

      ! globally sum angles' distributions information before saving

      If (angle%bin_adf > 0) Then

        ! maximum angle%dst that can be summed in each step

        nsum = config%mxbuff / (angle%bin_adf + 1)
        If (nsum == 0) Call error(200)

        Do i = 1, angle%ldf(0), nsum
          Call gsum(comm, angle%dst(:, i:Min(i + nsum - 1, angle%ldf(0))))
        End Do

      End If

      ! globally sum dihedrals' distributions information before saving

      If (dihedral%bin_adf > 0) Then

        ! maximum dihedral%dst that can be summed in each step

        nsum = config%mxbuff / (dihedral%bin_adf + 1)
        If (nsum == 0) Call error(200)

        Do i = 1, dihedral%ldf(0), nsum
          Call gsum(comm, dihedral%dst(:, i:Min(i + nsum - 1, dihedral%ldf(0))))
        End Do

      End If

      ! globally sum inversions' distributions information before saving

      If (inversion%bin_adf > 0) Then

        ! maximum inversion%dst that can be summed in each step

        nsum = config%mxbuff / (inversion%bin_adf + 1)
        If (nsum == 0) Call error(200)

        Do i = 1, inversion%ldf(0), nsum
          Call gsum(comm, inversion%dst(:, i:Min(i + nsum - 1, inversion%ldf(0))))
        End Do

      End If

    End If

    ! Write REVCON
    levcfg = 2 ! define level of information in REVCON
    Call write_config(config, files(FILE_REVCON), levcfg, nstep, thermo%tstep, io, time, netcdf, comm)

    ! node 0 handles I/O

    If (comm%idnode == 0) Then

      ! Write accumulator data to dump file

      If (devel%l_rout) Then
        Open (Newunit=files(FILE_REVIVE)%unit_no, File=files(FILE_REVIVE)%filename, Form='formatted', Status='replace')

        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') rcut, rdf%rbin, Real(config%megatm, wp)
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') &
          Real(nstep, wp), thermo%tstep, time, tmst, Real(stats%numacc, wp), thermo%chi_t, thermo%chi_p, thermo%cint
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') thermo%eta
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%stpval
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%stpvl0
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%sumval
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%ssqval
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%zumval
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%ravval
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%stkval
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%strcon
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%strpmf
        Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') stats%stress

        If (rdf%l_collect) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') Real(rdf%n_configs, wp), rdf%rdf
        End If
        If (rdf%max_grid_usr > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') &
            Real(rdf%max_grid_usr), rdf%cutoff_usr, Real(rdf%n_configs_usr, wp), rdf%usr
        End If
        If (zdensity%l_collect) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') &
            Real(zdensity%n_samples, wp), zdensity%density
        End If
        If (green%samp > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') green%vafcount
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') Real(green%step, wp)
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') green%vafdata
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') green%vaf
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') green%time
        End If

        If (bond%bin_pdf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') Real(bond%n_frames, wp), bond%dst
        End If
        If (angle%bin_adf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') Real(angle%n_frames, wp), angle%dst
        End If
        If (dihedral%bin_adf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') Real(dihedral%n_frames, wp), dihedral%dst
        End If
        If (inversion%bin_adf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') Real(inversion%n_frames, wp), inversion%dst
        End If
      Else
        Open (Newunit=files(FILE_REVIVE)%unit_no, File=files(FILE_REVIVE)%filename, Form='unformatted', Status='replace')

        Write (Unit=files(FILE_REVIVE)%unit_no) rcut, rdf%rbin, Real(config%megatm, wp)
        Write (Unit=files(FILE_REVIVE)%unit_no) &
          Real(nstep, wp), thermo%tstep, time, tmst, Real(stats%numacc, wp), thermo%chi_t, thermo%chi_p, thermo%cint
        Write (Unit=files(FILE_REVIVE)%unit_no) thermo%eta
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%stpval
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%stpvl0
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%sumval
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%ssqval
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%zumval
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%ravval
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%stkval
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%strcon
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%strpmf
        Write (Unit=files(FILE_REVIVE)%unit_no) stats%stress

        If (rdf%l_collect) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(rdf%n_configs, wp), rdf%rdf
        End If
        If (rdf%max_grid_usr > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(rdf%max_grid_usr), rdf%cutoff_usr, Real(rdf%n_configs_usr, wp), rdf%usr
        End If
        If (zdensity%l_collect) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(zdensity%n_samples, wp), zdensity%density
        End If
        If (green%samp > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) green%vafcount
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(green%step, wp)
          Write (Unit=files(FILE_REVIVE)%unit_no) green%vafdata
          Write (Unit=files(FILE_REVIVE)%unit_no) green%vaf
          Write (Unit=files(FILE_REVIVE)%unit_no) green%time
        End If

        If (bond%bin_pdf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(bond%n_frames, wp), bond%dst
        End If
        If (angle%bin_adf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(angle%n_frames, wp), angle%dst
        End If
        If (dihedral%bin_adf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(dihedral%n_frames, wp), dihedral%dst
        End If
        If (inversion%bin_adf > 0) Then
          Write (Unit=files(FILE_REVIVE)%unit_no) Real(inversion%n_frames, wp), inversion%dst
        End If
      End If

      ! Write initial position and final displacement data to REVIVE

      jatms = config%natms

      Do i = 1, config%natms
        iwrk(i) = config%ltg(i)

        axx(i) = stats%xin(i)
        ayy(i) = stats%yin(i)
        azz(i) = stats%zin(i)

        bxx(i) = stats%xto(i)
        byy(i) = stats%yto(i)
        bzz(i) = stats%zto(i)
      End Do

      ready = .true.
      Do jdnode = 0, comm%mxnode - 1
        If (jdnode > 0) Then
          Call gsend(comm, ready, jdnode, Revive_tag)

          Call grecv(comm, jatms, jdnode, Revive_tag)

          Call grecv(comm, iwrk(1:jatms), jdnode, Revive_tag)

          Call grecv(comm, axx(1:jatms), jdnode, Revive_tag)
          Call grecv(comm, ayy(1:jatms), jdnode, Revive_tag)
          Call grecv(comm, azz(1:jatms), jdnode, Revive_tag)

          Call grecv(comm, bxx(1:jatms), jdnode, Revive_tag)
          Call grecv(comm, byy(1:jatms), jdnode, Revive_tag)
          Call grecv(comm, bzz(1:jatms), jdnode, Revive_tag)
        End If

        If (devel%l_rout) Then
          Do i = 1, jatms
            Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') &
              Real(iwrk(i), wp), axx(i), ayy(i), azz(i), bxx(i), byy(i), bzz(i)
          End Do
        Else
          Do i = 1, jatms
            Write (Unit=files(FILE_REVIVE)%unit_no) Real(iwrk(i), wp), axx(i), ayy(i), azz(i), bxx(i), byy(i), bzz(i)
          End Do
        End If
      End Do

    Else

      Call grecv(comm, ready, 0, Revive_tag)

      Call gsend(comm, config%natms, 0, Revive_tag)

      Call gsend(comm, config%ltg(1:config%natms), 0, Revive_tag)

      Call gsend(comm, stats%xin(1:config%natms), 0, Revive_tag)
      Call gsend(comm, stats%yin(1:config%natms), 0, Revive_tag)
      Call gsend(comm, stats%zin(1:config%natms), 0, Revive_tag)

      Call gsend(comm, stats%xto(1:config%natms), 0, Revive_tag)
      Call gsend(comm, stats%yto(1:config%natms), 0, Revive_tag)
      Call gsend(comm, stats%zto(1:config%natms), 0, Revive_tag)

    End If

    ! Write initial velocities for VAF calculations if needed

    If (green%samp > 0) Then

      If (comm%idnode == 0) Then

        jatms = config%natms
        Do j = 1, green%samp
          Do i = 1, config%natms
            iwrk(i) = config%ltg(i)
            axx(i) = green%vxi(i, j)
            ayy(i) = green%vyi(i, j)
            azz(i) = green%vzi(i, j)
          End Do

          ready = .true.
          Do jdnode = 0, comm%mxnode - 1
            If (jdnode > 0) Then
              Call gsend(comm, ready, jdnode, Revive_tag)

              Call grecv(comm, jatms, jdnode, Revive_tag)

              Call grecv(comm, iwrk(1:jatms), jdnode, Revive_tag)

              Call grecv(comm, axx(1:jatms), jdnode, Revive_tag)
              Call grecv(comm, ayy(1:jatms), jdnode, Revive_tag)
              Call grecv(comm, azz(1:jatms), jdnode, Revive_tag)
            End If

            If (devel%l_rout) Then
              Do i = 1, jatms
                Write (Unit=files(FILE_REVIVE)%unit_no, Fmt=forma, Advance='No') &
                  Real(iwrk(i), wp), axx(i), ayy(i), azz(i)
              End Do
            Else
              Do i = 1, jatms
                Write (Unit=files(FILE_REVIVE)%unit_no) Real(iwrk(i), wp), axx(i), ayy(i), azz(i)
              End Do
            End If
          End Do
        End Do

      Else

        Do j = 1, green%samp
          Call grecv(comm, ready, 0, Revive_tag)

          Call gsend(comm, config%natms, 0, Revive_tag)

          Call gsend(comm, config%ltg(1:config%natms), 0, Revive_tag)

          Call gsend(comm, green%vxi(1, j), 0, Revive_tag)
          Call gsend(comm, green%vyi(1, j), 0, Revive_tag)
          Call gsend(comm, green%vzi(1, j), 0, Revive_tag)
        End Do

      End If

      Call gsync(comm)

    End If

    Call gsync(comm)
    If (comm%idnode == 0) Call files(FILE_REVIVE)%close ()

    r_mxnode = 1.0_wp / Real(comm%mxnode, wp)

    ! globally divide rdf%rdf data between nodes

    If (rdf%l_collect) rdf%rdf = rdf%rdf * r_mxnode

    ! globally divide USR rdf%rdf data between nodes

    If (rdf%max_grid_usr > 0) rdf%usr = rdf%usr * r_mxnode

    ! globally divide z-density data between nodes

    If (zdensity%l_collect) zdensity%density = zdensity%density * r_mxnode

    ! globally divide bonds' distributions data between nodes

    If (bond%bin_pdf > 0) bond%dst = bond%dst * r_mxnode

    ! globally divide angles' distributions data between nodes

    If (angle%bin_adf > 0) angle%dst = angle%dst * r_mxnode

    ! globally divide dihedrals' distributions data between nodes

    If (dihedral%bin_adf > 0) dihedral%dst = dihedral%dst * r_mxnode

    ! globally divide inversions' distributions data between nodes

    If (inversion%bin_adf > 0) inversion%dst = inversion%dst * r_mxnode

    Deallocate (iwrk, Stat=fail(1))
    Deallocate (axx, ayy, azz, Stat=fail(2))
    Deallocate (bxx, byy, bzz, Stat=fail(3))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'system_revive deallocation failure'
      Call error(0, message)
    End If

  End Subroutine system_revive
End Module system
