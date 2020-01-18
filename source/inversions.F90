Module inversions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global valence invle interaction variables
  ! and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov april 2014
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: comms_type,&
                             gbcast,&
                             gcheck,&
                             gsum,&
                             gsync
  Use configuration,   Only: configuration_type
  Use constants,       Only: boltz,&
                             delth_max,&
                             engunit,&
                             npdfdt,&
                             npdgdt,&
                             ntable,&
                             pi,&
                             zero_plus
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use kinds,           Only: wi,&
                             wp
  Use numerics,        Only: images,&
                             local_index
  Use parse,           Only: get_line,&
                             get_word,&
                             word_2_real
  Use site,            Only: site_type

  Implicit None

  Private

  ! Inversion potential keys
  !> Null potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_NULL = -1
  !> Tabulated potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_TAB = 0
  !> Harmonic potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_HARMONIC = 1
  !> Harmonic cosine potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_HARMONIC_COSINE = 2
  !> Planar potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_PLANAR = 3
  !> Extened planar potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_EXTENDED_PLANAR = 4
  !> Calcite four-body potential
  Integer(Kind=wi), Parameter, Public :: INVERSION_CALCITE = 5

  Type, Public :: inversions_type
    Private

    !> Tabulated potential switch
    Logical, Public :: l_tab = .false.
    !> Core shell switch
    Logical, Public :: l_core_shell = .false.
    !> Number of inversion angle types (potentials)
    Integer(Kind=wi), Public :: n_types = 0
    Integer(Kind=wi), Public :: n_types1 = 0
    !> Number of frames
    Integer(Kind=wi), Public :: n_frames = 0
    !> Total number of inversion angles (all nodes)
    Integer(Kind=wi), Public :: total

    Integer(Kind=wi), Allocatable, Public :: num(:)
    !> Inversion potential keys
    Integer(Kind=wi), Allocatable, Public :: key(:)
    !> Restrained potential flag
    Logical, Allocatable, Public          :: restrained(:)
    !> Atom indices (local)
    Integer(Kind=wi), Allocatable, Public :: lst(:, :)
    !> Atom indices
    Integer(Kind=wi), Allocatable, Public :: list(:, :)
    !> Legend
    Integer(Kind=wi), Allocatable, Public :: legend(:, :)
    !> Angle parameters (force constant, etc.)
    Real(Kind=wp), Allocatable, Public    :: param(:, :)
    ! Possible tabulated calculation arrays
    Integer, Allocatable, Public          :: ltp(:)
    !> Tabulated potential
    Real(Kind=wp), Allocatable, Public    :: tab_potential(:, :)
    !> Tabulated force
    Real(Kind=wp), Allocatable, Public    :: tab_force(:, :)
    ! Possible distribution arrays
    Integer, Allocatable, Public          :: ldf(:), typ(:, :)
    Real(Kind=wp), Allocatable, Public    :: dst(:, :)
    ! Maximums
    !> Maximum number of inversion angle types
    Integer(Kind=wi), Public              :: max_types
    !> Maximum number of inversion angles per node
    Integer(Kind=wi), Public              :: max_angles
    !> Length of legend array
    Integer(Kind=wi), Public              :: max_legend
    !> Maximum number of inversion parameters
    Integer(Kind=wi), Public              :: max_param
    ! Number of bins
    !> Angular distribution function bins
    Integer(Kind=wi), Public              :: bin_adf
    !> Tabulated potential bins
    Integer(Kind=wi), Public              :: bin_tab
  Contains
    Private

    Procedure, Public :: init => allocate_inversions_arrays
    Procedure, Public :: init_pot => allocate_invr_pot_arrays
    Procedure, Public :: init_dst => allocate_invr_dst_arrays
    Final             :: cleanup
  End Type inversions_type

  Public :: inversions_compute, inversions_forces, inversions_table_read

Contains

  Subroutine allocate_inversions_arrays(inversion, mxatdm, mxtmls)
    Class(inversions_type), Intent(InOut) :: inversion
    Integer(kind=wi),       Intent(In   ) :: mxatdm, mxtmls

    Integer, Dimension(9) :: fail

    fail = 0

    Allocate (inversion%num(1:mxtmls), Stat=fail(1))
    Allocate (inversion%key(1:inversion%max_types), Stat=fail(2))
    Allocate (inversion%restrained(1:inversion%max_types), Stat=fail(3))
    Allocate (inversion%lst(1:4, 1:inversion%max_types), Stat=fail(4))
    Allocate (inversion%list(0:4, 1:inversion%max_angles), Stat=fail(5))
    Allocate (inversion%legend(0:inversion%max_legend, 1:mxatdm), Stat=fail(6))
    Allocate (inversion%param(1:inversion%max_param, 1:inversion%max_types), Stat=fail(7))
    If (inversion%l_tab) &
      Allocate (inversion%ltp(0:inversion%max_types), Stat=fail(8))
    If (inversion%bin_adf > 0) &
      Allocate (inversion%ldf(0:inversion%max_types), Stat=fail(9))

    If (Any(fail > 0)) Call error(1021)

    inversion%num = 0
    inversion%key = INVERSION_NULL
    inversion%restrained = .false.
    inversion%lst = 0
    inversion%list = 0
    inversion%legend = 0

    inversion%param = 0.0_wp

    If (inversion%l_tab) &
      inversion%ltp = 0

    If (inversion%bin_adf > 0) &
      inversion%ldf = 0

  End Subroutine allocate_inversions_arrays

  Subroutine allocate_invr_pot_arrays(inversion)
    Class(inversions_type), Intent(InOut) :: inversion

    Integer :: fail(1:2)

    fail = 0

    Allocate (inversion%tab_potential(-1:inversion%bin_tab, 1:inversion%ltp(0)), Stat=fail(1))
    Allocate (inversion%tab_force(-1:inversion%bin_tab, 1:inversion%ltp(0)), Stat=fail(2))

    If (Any(fail > 0)) Call error(1078)

    inversion%tab_potential = 0.0_wp
    inversion%tab_force = 0.0_wp

  End Subroutine allocate_invr_pot_arrays

  Subroutine allocate_invr_dst_arrays(inversion)
    Class(inversions_type), Intent(InOut) :: inversion

    Integer :: fail

    fail = 0

    Allocate (inversion%typ(-1:4, 1:inversion%ldf(0)), inversion%dst(1:inversion%bin_adf, 1:inversion%ldf(0)), Stat=fail)

    If (fail > 0) Call error(1079)

    inversion%typ = 0
    inversion%dst = 0.0_wp

  End Subroutine allocate_invr_dst_arrays

  Subroutine inversions_compute(temp, unique_atom, inversion, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating inversions distribution functions
    ! from accumulated data
    !
    ! copyright - daresbury laboratory
    ! author    - a.v.brukhno & i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),                  Intent(In   ) :: temp
    Character(Len=8), Dimension(:), Intent(In   ) :: unique_atom
    Type(inversions_type),          Intent(InOut) :: inversion
    Type(configuration_type),       Intent(InOut) :: config
    Type(comms_type),               Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail, i, ig, j, kk, ll, ngrid
    Logical                    :: zero
    Real(Kind=wp)              :: coef, delth, dfed, dfed0, dfed1, dfed2, dgr2rad, dgrid, factor, &
                                  factor1, fed, fed0, fed1, fed2, kT2engo, pdfinv, pdfinv1, &
                                  pdfzero, rad2dgr, rdlth, sum, sum1, t1, t2, theta, tmp
    Real(Kind=wp), Allocatable :: dstdinv(:, :), pmf(:), vir(:)

    fail = 0
  Allocate (dstdinv(0:inversion%bin_adf, 1:inversion%ldf(0)), pmf(0:inversion%bin_adf + 2), vir(0:inversion%bin_adf + 2), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'inversions_compute - allocation failure'
      Call error(0, message)
    End If

    ! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

    kT2engo = boltz * temp / engunit

    ! conversion: radians <-> degrees (affects not only angle units but also force units!)

    rad2dgr = 180.0_wp / pi
    dgr2rad = pi / 180.0_wp

    ! grid interval for pdf tables

    delth = pi / Real(inversion%bin_adf, wp)
    rdlth = Real(inversion%bin_adf, wp) / 180.0_wp

    ! resampling grid and grid interval for pmf tables

    ngrid = Max(Nint(360.0_wp / delth_max), inversion%bin_adf, inversion%bin_tab - 4)
    dgrid = pi / Real(ngrid, wp)

    ! loop over all valid PDFs to get valid totals

    kk = 0
    ll = 0
    Do i = 1, inversion%ldf(0)
      If (inversion%typ(0, i) > 0) Then
        kk = kk + 1
        ll = ll + inversion%typ(0, i)
      End If
    End Do

    ! normalisation factor

    factor = 1.0_wp / Real(inversion%n_frames, wp)

    ! the lower bound to nullify the nearly-zero histogram (PDF) values

    pdfzero = 1.0e-5_wp

    Call info('', .true.)
    Call info('INVERSIONS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)', .true.)
    Write (message, '(a,5(1x,i10))') '# bins, cutoff, frames, types: ', &
      inversion%bin_adf, 180, inversion%n_frames, kk, ll
    Call info(message, .true.)

    ! open RDF file and write headers

    If (comm%idnode == 0) Then
      Open (Unit=npdfdt, File='INVDAT', Status='replace')
      Write (npdfdt, '(a)') '# '//config%cfgname
      Write (npdfdt, '(a)') '# INVERSIONS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
      Write (npdfdt, '(a,4(1x,i10))') '# bins, cutoff, frames, types: ', inversion%bin_adf, 180, inversion%n_frames, kk
      Write (npdfdt, '(a)') '#'
      Write (npdfdt, '(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)   @   dTheta_bin = ', delth * rad2dgr
      Write (npdfdt, '(a)') '#'
    End If

    ! loop over all valid PDFs

    j = 0
    Do i = 1, inversion%ldf(0)
      If (inversion%typ(0, i) > 0) Then
        j = j + 1

        Write (message, '(a,4(a8,1x),2(i10,1x))') 'type, index, instances: ', &
          unique_atom(inversion%typ(1, i)), unique_atom(inversion%typ(2, i)), &
          unique_atom(inversion%typ(3, i)), &
          unique_atom(inversion%typ(4, i)), j, inversion%typ(0, i)
        Call info(message, .true.)
        Write (message, '(a,f8.5)') &
          'Theta(degrees)  P_inv(Theta)  Sum_P_inv(Theta)   @   dTheta_bin = ', &
          delth * rad2dgr
        Call info(message, .true.)
        If (comm%idnode == 0) Then
          Write (npdfdt, '(/,a,4(a8,1x),2(i10,1x))') '# type, index, instances: ', &
            unique_atom(inversion%typ(1, i)), unique_atom(inversion%typ(2, i)), &
            unique_atom(inversion%typ(3, i)), &
            unique_atom(inversion%typ(4, i)), j, inversion%typ(0, i)
        End If

        ! global sum of data on all nodes

        Call gsum(comm, inversion%dst(1:inversion%bin_adf, i))

        ! factor in instances (first, pdfinv is normalised to unity)

        factor1 = factor / Real(inversion%typ(0, i), wp)

        ! running integration of pdf

        sum = 0.0_wp

        ! loop over distances

        zero = .true.
        Do ig = 1, inversion%bin_adf
          If (zero .and. ig < (inversion%bin_adf - 3)) zero = (inversion%dst(ig + 2, i) <= 0.0_wp)

          pdfinv = inversion%dst(ig, i) * factor1
          sum = sum + pdfinv

          ! null it if < pdfzero

          If (pdfinv < pdfzero) Then
            pdfinv1 = 0.0_wp
          Else
            pdfinv1 = pdfinv
          End If

          If (sum < pdfzero) Then
            sum1 = 0.0_wp
          Else
            sum1 = sum
          End If

          theta = (Real(ig, wp) - 0.5_wp) * delth

          ! now pdfinv is normalised by the volume element (as to go to unity at infinity in gases and liquids)

          pdfinv = pdfinv * rdlth

          ! print out information

          theta = theta * rad2dgr
          If (.not. zero) Then
            Write (message, "(f11.5,1p,2e14.6)") theta, pdfinv1, sum1
            Call info(message, .true.)
          End If
          If (comm%idnode == 0) Then
            Write (npdfdt, "(f11.5,1p,e14.6)") theta, pdfinv
          End If

          ! We use the non-normalised tail-truncated PDF version,
          ! pdf...1 (not pdf...) in order to exclude the nearly-zero
          ! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
          ! would have poorly-defined noisy "borders/walls"

          dstdinv(ig, i) = pdfinv1 ! PDFs density
        End Do
      Else
        dstdinv(:, i) = 0.0_wp ! PDFs density
      End If
    End Do

    If (comm%idnode == 0) Close (Unit=npdfdt)

    ! open PDF files and write headers

    If (comm%idnode == 0) Then
      Open (Unit=npdgdt, File='INVPMF', Status='replace')
      Write (npdgdt, '(a)') '# '//config%cfgname
      Write (npdgdt, '(a,i10,2f12.5,i10,a,e15.7)') '# ', inversion%bin_adf, &
        delth * Real(inversion%bin_adf, wp) * rad2dgr, delth * rad2dgr, kk, &
        '   conversion factor(kT -> energy units) =', kT2engo

      Open (Unit=npdfdt, File='INVPMF', Status='replace')
      Write (npdfdt, '(a)') '# '//config%cfgname
      Write (npdfdt, '(a,i10,2f12.5,i10,a,e15.7)') '# ', ngrid, &
        dgrid * Real(ngrid, wp) * rad2dgr, dgrid * rad2dgr, kk, &
        '   conversion factor(kT -> energy units) =', kT2engo
    End If

    ! loop over all valid PDFs

    j = 0
    Do i = 1, inversion%ldf(0)
      If (inversion%typ(0, i) > 0) Then
        j = j + 1

        If (comm%idnode == 0) Then
          Write (npdgdt, '(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
            unique_atom(inversion%typ(1, i)), unique_atom(inversion%typ(2, i)), &
            unique_atom(inversion%typ(3, i)), &
            unique_atom(inversion%typ(4, i)), j, inversion%typ(0, i), ' (type, index, instances)'
          Write (npdfdt, '(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
            unique_atom(inversion%typ(1, i)), unique_atom(inversion%typ(2, i)), &
            unique_atom(inversion%typ(3, i)), &
            unique_atom(inversion%typ(4, i)), j, inversion%typ(0, i), ' (type, index, instances)'
        End If

        ! Smoothen and get derivatives

        fed0 = 0.0_wp
        dfed0 = 10.0_wp
        dfed = 10.0_wp

        Do ig = 1, inversion%bin_adf
          tmp = Real(ig, wp) - 0.5_wp
          theta = tmp * delth

          If (dstdinv(ig, i) > zero_plus) Then
            fed = -Log(dstdinv(ig, i)) - fed0
            If (fed0 <= zero_plus) Then
              fed0 = fed
              fed = 0.0_wp
            End If

            If (ig < inversion%bin_adf - 1) Then
              If (dstdinv(ig + 1, i) <= zero_plus .and. dstdinv(ig + 2, i) > zero_plus) &
                dstdinv(ig + 1, i) = 0.5_wp * (dstdinv(ig, i) + dstdinv(ig + 2, i))
            End If
          Else
            fed = 0.0_wp
          End If

          If (ig == 1) Then
            If (dstdinv(ig, i) > zero_plus .and. dstdinv(ig + 1, i) > zero_plus) Then
              dfed = Log(dstdinv(ig + 1, i) / dstdinv(ig, i))
            Else If (dfed > 0.0_wp) Then
              dfed = dfed0
            Else
              dfed = -dfed0
            End If
          Else If (ig == inversion%bin_adf) Then
            If (dstdinv(ig, i) > zero_plus .and. dstdinv(ig - 1, i) > zero_plus) Then
              dfed = Log(dstdinv(ig, i) / dstdinv(ig - 1, i))
            Else If (dfed > 0.0_wp) Then
              dfed = dfed0
            Else
              dfed = -dfed0
            End If
          Else If (dstdinv(ig - 1, i) > zero_plus) Then
            If (dstdinv(ig + 1, i) > zero_plus) Then
              dfed = 0.5_wp * (Log(dstdinv(ig + 1, i) / dstdinv(ig - 1, i)))
            Else
              dfed = 0.5_wp * Log(dstdinv(ig - 1, i))
            End If
          Else If (dstdinv(ig + 1, i) > zero_plus) Then
            dfed = -0.5_wp * Log(dstdinv(ig + 1, i))
          Else If (dfed > 0.0_wp) Then
            dfed = dfed0
          Else
            dfed = -dfed0
          End If

          pmf(ig) = fed
          vir(ig) = dfed

          ! Print

          If (comm%idnode == 0) Write (npdgdt, "(f11.5,1p,2e14.6)") theta * rad2dgr, fed * kT2engo, dfed * kT2engo * dgr2rad / delth
        End Do

        ! Define edges

        pmf(0) = 2.0_wp * pmf(1) - pmf(2)
        vir(0) = 2.0_wp * vir(1) - vir(2)
        pmf(inversion%bin_adf + 1) = 2.0_wp * pmf(inversion%bin_adf) - pmf(inversion%bin_adf - 1)
        vir(inversion%bin_adf + 1) = 2.0_wp * vir(inversion%bin_adf) - vir(inversion%bin_adf - 1)
        pmf(inversion%bin_adf + 2) = 2.0_wp * pmf(inversion%bin_adf + 1) - pmf(inversion%bin_adf)
        vir(inversion%bin_adf + 2) = 2.0_wp * vir(inversion%bin_adf + 1) - vir(inversion%bin_adf)

        ! resample using 3pt interpolation

        Do ig = 1, ngrid
          theta = Real(ig, wp) * dgrid
          ll = Int(theta / delth)

          ! +0.5_wp due to half-a-bin shift in the original data

          coef = theta / delth - Real(ll, wp) + 0.5_wp

          fed0 = pmf(ll)
          fed1 = pmf(ll + 1)
          fed2 = pmf(ll + 2)

          t1 = fed0 + (fed1 - fed0) * coef
          t2 = fed1 + (fed2 - fed1) * (coef - 1.0_wp)

          fed = t1 + (t2 - t1) * coef * 0.5_wp

          dfed0 = vir(ll)
          dfed1 = vir(ll + 1)
          dfed2 = vir(ll + 2)

          t1 = dfed0 + (dfed1 - dfed0) * coef
          t2 = dfed1 + (dfed2 - dfed1) * (coef - 1.0_wp)

          dfed = t1 + (t2 - t1) * coef * 0.5_wp

          If (comm%idnode == 0) &
            Write (npdfdt, "(f11.5,1p,2e14.6)") theta * rad2dgr, fed * kT2engo, dfed * kT2engo * dgr2rad / delth
        End Do
      End If
    End Do

    If (comm%idnode == 0) Then
      Close (Unit=npdgdt)
      Close (Unit=npdfdt)
    End If

    Deallocate (dstdinv, pmf, vir, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'inversions_compute - deallocation failure'
      Call error(0, message)
    End If

  End Subroutine inversions_compute

  Subroutine inversions_forces(isw, enginv, virinv, stress, inversion, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating inversion energy and force terms
    !
    ! isw = 0 - collect statistics
    ! isw = 1 - calculate forces
    ! isw = 2 - do both
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith may 1996
    ! amended   - i.t.todorov february 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                       Intent(In   ) :: isw
    Real(Kind=wp),                 Intent(  Out) :: enginv, virinv
    Real(Kind=wp), Dimension(1:9), Intent(InOut) :: stress
    Type(inversions_type),         Intent(InOut) :: inversion
    Type(configuration_type),      Intent(InOut) :: config
    Type(comms_type),              Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail(1:4), i, ia, ib, ic, id, j, keyi, kk, l
    Integer, Allocatable       :: lstopt(:, :)
    Logical                    :: safe
    Logical, Allocatable       :: lunsafe(:)
    Real(Kind=wp)              :: a, b, buffer(1:2), cos0, cosb, cosc, cosd, fax, fay, faz, fbx, &
                                  fby, fbz, fcx, fcy, fcz, fdx, fdy, fdz, gamb, gamc, gamd, gamma, &
                                  k, m, ppp, pterm, rab2, rac2, rad2, rbc, rcd, rdb, rdelth, rdr, &
                                  rrab, rrac, rrad, rub, rubc, rubd, ruc, rucb, rucd, rud, rudb, &
                                  rudc, rvb, rvbc, rvbd, rvc, rvcb, rvcd, rvd, rvdb, rvdc, strs1, &
                                  strs2, strs3, strs5, strs6, strs9, t1, t2, th0, thb, thc, thd, &
                                  ubn, ubx, uby, ubz, ucn, ucx, ucy, ucz, udn, udx, udy, udz, uu2, &
                                  uun, uuu, uux, uuy, uuz, vbn, vbx, vby, vbz, vcn, vcx, vcy, vcz, &
                                  vdn, vdx, vdy, vdz, vk, vk1, vk2, vterm, wwb, wwc, wwd, xab, &
                                  xac, xad, yab, yac, yad, zab, zac, zad
    Real(Kind=wp), Allocatable :: xdab(:), xdac(:), xdad(:), ydab(:), ydac(:), ydad(:), zdab(:), &
                                  zdac(:), zdad(:)

    fail = 0
    Allocate (lunsafe(1:inversion%max_angles), lstopt(0:4, 1:inversion%max_angles), Stat=fail(1))
    Allocate (xdab(1:inversion%max_angles), ydab(1:inversion%max_angles), zdab(1:inversion%max_angles), Stat=fail(2))
    Allocate (xdac(1:inversion%max_angles), ydac(1:inversion%max_angles), zdac(1:inversion%max_angles), Stat=fail(3))
    Allocate (xdad(1:inversion%max_angles), ydad(1:inversion%max_angles), zdad(1:inversion%max_angles), Stat=fail(4))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'inversions_forces allocation failure'
      Call error(0, message)
    End If

    ! calculate atom separation vectors

    Do i = 1, inversion%n_types
      lunsafe(i) = .false.

      ! indices of atoms involved

      ia = local_index(inversion%list(1, i), config%nlast, config%lsi, config%lsa); lstopt(1, i) = ia
      ib = local_index(inversion%list(2, i), config%nlast, config%lsi, config%lsa); lstopt(2, i) = ib
      ic = local_index(inversion%list(3, i), config%nlast, config%lsi, config%lsa); lstopt(3, i) = ic
      id = local_index(inversion%list(4, i), config%nlast, config%lsi, config%lsa); lstopt(4, i) = id

      lstopt(0, i) = 0
      If (ia > 0 .and. ib > 0 .and. ic > 0 .and. id > 0) Then !Tag
        If (config%lfrzn(ia) * config%lfrzn(ib) * config%lfrzn(ic) * config%lfrzn(id) == 0) Then
          If (ia <= config%natms .or. ib <= config%natms .or. &
              ic <= config%natms .or. id <= config%natms) Then
            lstopt(0, i) = 1
          End If
        End If
      Else ! Detect uncompressed unit
        If (((ia > 0 .and. ia <= config%natms) .or. &
             (ib > 0 .and. ib <= config%natms) .or. &
             (ic > 0 .and. ic <= config%natms) .or. &
             (id > 0 .and. id <= config%natms)) .and. &
            (ia == 0 .or. ib == 0 .or. ic == 0 .or. id == 0)) lunsafe(i) = .true.
      End If

      ! define components of bond vectors

      If (lstopt(0, i) > 0) Then
        xdab(i) = config%parts(ib)%xxx - config%parts(ia)%xxx
        ydab(i) = config%parts(ib)%yyy - config%parts(ia)%yyy
        zdab(i) = config%parts(ib)%zzz - config%parts(ia)%zzz

        ! select potential energy function type

        kk = inversion%list(0, i)
        keyi = inversion%key(kk)

        If (keyi == INVERSION_CALCITE) Then
          xdac(i) = config%parts(ic)%xxx - config%parts(ib)%xxx
          ydac(i) = config%parts(ic)%yyy - config%parts(ib)%yyy
          zdac(i) = config%parts(ic)%zzz - config%parts(ib)%zzz

          xdad(i) = config%parts(id)%xxx - config%parts(ib)%xxx
          ydad(i) = config%parts(id)%yyy - config%parts(ib)%yyy
          zdad(i) = config%parts(id)%zzz - config%parts(ib)%zzz
        Else
          xdac(i) = config%parts(ic)%xxx - config%parts(ia)%xxx
          ydac(i) = config%parts(ic)%yyy - config%parts(ia)%yyy
          zdac(i) = config%parts(ic)%zzz - config%parts(ia)%zzz

          xdad(i) = config%parts(id)%xxx - config%parts(ia)%xxx
          ydad(i) = config%parts(id)%yyy - config%parts(ia)%yyy
          zdad(i) = config%parts(id)%zzz - config%parts(ia)%zzz
        End If
      Else ! (DEBUG)
        xdab(i) = 0.0_wp
        ydab(i) = 0.0_wp
        zdab(i) = 0.0_wp

        xdac(i) = 0.0_wp
        ydac(i) = 0.0_wp
        zdac(i) = 0.0_wp

        xdad(i) = 0.0_wp
        ydad(i) = 0.0_wp
        zdad(i) = 0.0_wp
      End If
    End Do

    ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:inversion%n_types))
    Call gcheck(comm, safe)
    If (.not. safe) Then
      Do j = 0, comm%mxnode - 1
        If (comm%idnode == j) Then
          Do i = 1, inversion%n_types
            If (lunsafe(i)) Then
              Write (message, '(2(a,i10))') 'global unit number', inversion%list(0, i), &
                ' , with a head particle number', inversion%list(1, i)
              Call info(message)
              Call warning('contributes towards next error')
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(134)
    End If

    ! periodic boundary condition

    Call images(config%imcon, config%cell, inversion%n_types, xdab, ydab, zdab)
    Call images(config%imcon, config%cell, inversion%n_types, xdac, ydac, zdac)
    Call images(config%imcon, config%cell, inversion%n_types, xdad, ydad, zdad)

    If (Mod(isw, 3) > 0) Then

      ! Initialise safety flag

      safe = .true.

      ! zero inversion energy accumulator

      enginv = 0.0_wp
      virinv = 0.0_wp

      ! initialise stress tensor accumulators

      strs1 = 0.0_wp
      strs2 = 0.0_wp
      strs3 = 0.0_wp
      strs5 = 0.0_wp
      strs6 = 0.0_wp
      strs9 = 0.0_wp

    End If

    ! Recover bin size and increment counter

    If (Mod(isw, 2) == 0) Then
      rdelth = Real(inversion%bin_adf, wp) / pi
      inversion%n_frames = inversion%n_frames + 1
    End If

    ! loop over all specified inversions

    Do i = 1, inversion%n_types
      If (lstopt(0, i) > 0) Then

        ! indices of atoms involved

        ia = lstopt(1, i)
        ib = lstopt(2, i)
        ic = lstopt(3, i)
        id = lstopt(4, i)

        ! define components of bond vectors

        xab = xdab(i)
        yab = ydab(i)
        zab = zdab(i)
        rab2 = xab * xab + yab * yab + zab * zab
        rrab = 1.0_wp / Sqrt(rab2)

        xac = xdac(i)
        yac = ydac(i)
        zac = zdac(i)
        rac2 = xac * xac + yac * yac + zac * zac
        rrac = 1.0_wp / Sqrt(rac2)

        xad = xdad(i)
        yad = ydad(i)
        zad = zdad(i)
        rad2 = xad * xad + yad * yad + zad * zad
        rrad = 1.0_wp / Sqrt(rad2)

        ! select potential energy function type

        kk = inversion%list(0, i)
        keyi = inversion%key(kk)

        If (keyi == INVERSION_CALCITE) Then

          ! calculate vector normal to plane

          uux = yac * zad - zac * yad
          uuy = zac * xad - xac * zad
          uuz = xac * yad - yac * xad
          uun = 1.0_wp / Sqrt(uux**2 + uuy**2 + uuz**2)
          uux = uun * uux
          uuy = uun * uuy
          uuz = uun * uuz
          uuu = xab * uux + yab * uuy + zab * uuz

        Else

          ! scalar products of bond vectors

          rbc = xab * xac + yab * yac + zab * zac
          rcd = xac * xad + yac * yad + zac * zad
          rdb = xad * xab + yad * yab + zad * zab

          ! calculate bond-angle-plane vectors

          ubx = xac * rrac + xad * rrad
          uby = yac * rrac + yad * rrad
          ubz = zac * rrac + zad * rrad
          ubn = 1.0_wp / Sqrt(ubx**2 + uby**2 + ubz**2)
          ubx = ubn * ubx
          uby = ubn * uby
          ubz = ubn * ubz
          rub = xab * ubx + yab * uby + zab * ubz

          vbx = xac * rrac - xad * rrad
          vby = yac * rrac - yad * rrad
          vbz = zac * rrac - zad * rrad
          vbn = 1.0_wp / Sqrt(vbx**2 + vby**2 + vbz**2)
          vbx = vbn * vbx
          vby = vbn * vby
          vbz = vbn * vbz
          rvb = xab * vbx + yab * vby + zab * vbz
          wwb = Sqrt(rub**2 + rvb**2)

          ucx = xad * rrad + xab * rrab
          ucy = yad * rrad + yab * rrab
          ucz = zad * rrad + zab * rrab
          ucn = 1.0_wp / Sqrt(ucx**2 + ucy**2 + ucz**2)
          ucx = ucn * ucx
          ucy = ucn * ucy
          ucz = ucn * ucz
          ruc = xac * ucx + yac * ucy + zac * ucz

          vcx = xad * rrad - xab * rrab
          vcy = yad * rrad - yab * rrab
          vcz = zad * rrad - zab * rrab
          vcn = 1.0_wp / Sqrt(vcx**2 + vcy**2 + vcz**2)
          vcx = vcn * vcx
          vcy = vcn * vcy
          vcz = vcn * vcz
          rvc = xac * vcx + yac * vcy + zac * vcz
          wwc = Sqrt(ruc**2 + rvc**2)

          udx = xab * rrab + xac * rrac
          udy = yab * rrab + yac * rrac
          udz = zab * rrab + zac * rrac
          udn = 1.0_wp / Sqrt(udx**2 + udy**2 + udz**2)
          udx = udn * udx
          udy = udn * udy
          udz = udn * udz
          rud = xad * udx + yad * udy + zad * udz

          vdx = xab * rrab - xac * rrac
          vdy = yab * rrab - yac * rrac
          vdz = zab * rrab - zac * rrac
          vdn = 1.0_wp / Sqrt(vdx**2 + vdy**2 + vdz**2)
          vdx = vdn * vdx
          vdy = vdn * vdy
          vdz = vdn * vdz
          rvd = xad * vdx + yad * vdy + zad * vdz
          wwd = Sqrt(rud**2 + rvd**2)

          ! calculate inversion angle cosines

          cosb = wwb * rrab; If (Abs(cosb) > 1.0_wp) cosb = Sign(1.0_wp, cosb)
          cosc = wwc * rrac; If (Abs(cosc) > 1.0_wp) cosc = Sign(1.0_wp, cosb)
          cosd = wwd * rrad; If (Abs(cosd) > 1.0_wp) cosd = Sign(1.0_wp, cosb)

          ! accumulate the histogram (distribution)

          If (Mod(isw, 2) == 0 .and. ib <= config%natms) Then
            j = inversion%ldf(kk)

            thb = Acos(cosb)
            l = Min(1 + Int(thb * rdelth), inversion%bin_adf)
            inversion%dst(l, j) = inversion%dst(l, j) + 1.0_wp / 3.0_wp

            thc = Acos(cosc)
            l = Min(1 + Int(thc * rdelth), inversion%bin_adf)
            inversion%dst(l, j) = inversion%dst(l, j) + 1.0_wp / 3.0_wp

            thd = Acos(cosd)
            l = Min(1 + Int(thd * rdelth), inversion%bin_adf)
            inversion%dst(l, j) = inversion%dst(l, j) + 1.0_wp / 3.0_wp
          End If

        End If
        If (isw == 0) Cycle

        ! calculate potential energy and scalar force term

        If (keyi == INVERSION_HARMONIC) Then

          ! harmonic inversion potential

          k = inversion%param(1, kk) / 6.0_wp
          th0 = inversion%param(2, kk)

          thb = Acos(cosb)
          thc = Acos(cosc)
          thd = Acos(cosd)

          pterm = k * ((thb - th0)**2 + (thc - th0)**2 + (thd - th0)**2)
          vterm = 0.0_wp
          gamma = 0.0_wp
          gamb = 0.0_wp
          If (Abs(thb) > 1.0e-10_wp) gamb = 2.0_wp * k * (thb - th0) / Sin(thb)
          gamc = 0.0_wp
          If (Abs(thc) > 1.0e-10_wp) gamc = 2.0_wp * k * (thc - th0) / Sin(thc)
          gamd = 0.0_wp
          If (Abs(thd) > 1.0e-10_wp) gamd = 2.0_wp * k * (thd - th0) / Sin(thd)

        Else If (keyi == INVERSION_HARMONIC_COSINE) Then

          ! harmonic cosine inversion potential

          k = inversion%param(1, kk) / 6.0_wp
          cos0 = inversion%param(2, kk)

          pterm = k * ((cosb - cos0)**2 + (cosc - cos0)**2 + (cosd - cos0)**2)
          vterm = 0.0_wp
          gamma = 0.0_wp
          gamb = -2.0_wp * k * (cosb - cos0)
          gamc = -2.0_wp * k * (cosc - cos0)
          gamd = -2.0_wp * k * (cosd - cos0)

        Else If (keyi == INVERSION_PLANAR) Then

          ! planar inversion potentials

          a = inversion%param(1, kk)

          pterm = a * (1.0_wp - (cosb + cosc + cosd) / 3.0_wp)
          vterm = 0.0_wp
          gamma = 0.0_wp
          gamb = a / 3.0_wp
          gamc = a / 3.0_wp
          gamd = a / 3.0_wp

        Else If (keyi == INVERSION_EXTENDED_PLANAR) Then

          ! extended planar inversion potentials

          k = inversion%param(1, kk) / 6.0_wp
          th0 = inversion%param(2, kk)
          m = inversion%param(3, kk)

          thb = Acos(cosb)
          thc = Acos(cosc)
          thd = Acos(cosd)

          pterm = 3.0_wp * k * (1.0_wp - (Cos(m * thb - th0) + Cos(m * thc - th0) + Cos(m * thd - th0)) / 3.0_wp)
          vterm = 0.0_wp
          gamma = 0.0_wp
          gamb = 0.0_wp
          If (Abs(thb) > 1.0e-10_wp) gamb = k * Sin(m * thb - th0) / Sin(thb)
          gamc = 0.0_wp
          If (Abs(thc) > 1.0e-10_wp) gamc = k * Sin(m * thc - th0) / Sin(thc)
          gamd = 0.0_wp
          If (Abs(thd) > 1.0e-10_wp) gamd = k * Sin(m * thd - th0) / Sin(thd)

        Else If (keyi == INVERSION_CALCITE) Then

          ! planar calcite potential

          a = inversion%param(1, kk)
          b = inversion%param(2, kk)

          uu2 = uuu * uuu
          m = 2.0_wp * a + 4.0_wp * b * uu2

          pterm = uu2 * (a + b * uu2)
          vterm = uu2 * m
          gamma = -uuu * m
          gamb = 0.0_wp
          gamc = 0.0_wp
          gamd = 0.0_wp

        Else If (keyi == INVERSION_TAB) Then

          ! TABINV potential

          pterm = 0.0_wp

          j = inversion%ltp(kk)
          rdr = inversion%tab_force(-1, j) ! 1.0_wp/delpot (in rad^-1)

          thb = Acos(cosb)
          thc = Acos(cosc)
          thd = Acos(cosd)

          l = Int(thb * rdr)
          ppp = thb * rdr - Real(l, wp)

          vk = inversion%tab_potential(l, j)
          vk1 = inversion%tab_potential(l + 1, j)
          vk2 = inversion%tab_potential(l + 2, j)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          pterm = pterm + t1 + (t2 - t1) * ppp * 0.5_wp

          vk = inversion%tab_force(l, j); If (l == 0) vk = vk * thb
          vk1 = inversion%tab_force(l + 1, j)
          vk2 = inversion%tab_force(l + 2, j)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          gamb = -t1 + (t2 - t1) * ppp * 0.5_wp

          l = Int(thc * rdr)
          ppp = thc * rdr - Real(l, wp)

          vk = inversion%tab_potential(l, j)
          vk1 = inversion%tab_potential(l + 1, j)
          vk2 = inversion%tab_potential(l + 2, j)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          pterm = pterm + t1 + (t2 - t1) * ppp * 0.5_wp

          vk = inversion%tab_force(l, j); If (l == 0) vk = vk * thc
          vk1 = inversion%tab_force(l + 1, j)
          vk2 = inversion%tab_force(l + 2, j)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          gamc = -t1 + (t2 - t1) * ppp * 0.5_wp

          l = Int(thd * rdr)
          ppp = thd * rdr - Real(l, wp)

          vk = Merge(inversion%tab_potential(l, j), 0.0_wp, l > 0)
          vk1 = inversion%tab_potential(l + 1, j)
          vk2 = inversion%tab_potential(l + 2, j)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          pterm = pterm + t1 + (t2 - t1) * ppp * 0.5_wp

          vk = inversion%tab_force(l, j); If (l == 0) vk = vk * thd
          vk1 = inversion%tab_force(l + 1, j)
          vk2 = inversion%tab_force(l + 2, j)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          gamd = -t1 + (t2 - t1) * ppp * 0.5_wp

          vterm = 0.0_wp
          gamma = 0.0_wp

        Else

          ! undefined potential

          safe = .false.
          pterm = 0.0_wp
          vterm = 0.0_wp
          gamb = 0.0_wp
          gamc = 0.0_wp
          gamd = 0.0_wp

        End If

        If (inversion%key(kk) == INVERSION_CALCITE) Then

          ! calculate atomic forces

          fax = -gamma * uux
          fay = -gamma * uuy
          faz = -gamma * uuz

          fcx = gamma * uun * ((yad * zab - zad * yab) - uuu * (yad * uuz - zad * uuy))
          fcy = gamma * uun * ((zad * xab - xad * zab) - uuu * (zad * uux - xad * uuz))
          fcz = gamma * uun * ((xad * yab - yad * xab) - uuu * (xad * uuy - yad * uux))

          fdx = gamma * uun * ((yab * zac - zab * yac) - uuu * (zac * uuy - yac * uuz))
          fdy = gamma * uun * ((zab * xac - xab * zac) - uuu * (xac * uuz - zac * uux))
          fdz = gamma * uun * ((xab * yac - yab * xac) - uuu * (yac * uux - xac * uuy))

          fbx = -(fax + fcx + fdx)
          fby = -(fay + fcy + fdy)
          fbz = -(faz + fcz + fdz)

        Else

          ! calculate bond and u,v scalar products

          rubc = xab * ucx + yab * ucy + zab * ucz
          rubd = xab * udx + yab * udy + zab * udz
          rucd = xac * udx + yac * udy + zac * udz
          rucb = xac * ubx + yac * uby + zac * ubz
          rudb = xad * ubx + yad * uby + zad * ubz
          rudc = xad * ucx + yad * ucy + zad * ucz

          rvbc = xab * vcx + yab * vcy + zab * vcz
          rvbd = xab * vdx + yab * vdy + zab * vdz
          rvcd = xac * vdx + yac * vdy + zac * vdz
          rvcb = xac * vbx + yac * vby + zac * vbz
          rvdb = xad * vbx + yad * vby + zad * vbz
          rvdc = xad * vcx + yad * vcy + zad * vcz

          ! calculate atomic forces

          fbx = gamb * (-cosb * xab * rrab**2 + rrab * (rub * ubx + rvb * vbx) / wwb) + &
                (ruc * ucn * rrab * (xac - ruc * ucx - (rbc - ruc * rubc) * xab * rrab**2) - &
                 rvc * vcn * rrab * (xac - rvc * vcx - (rbc - rvc * rvbc) * xab * rrab**2)) * &
                gamc * rrac / wwc + &
                (rud * udn * rrab * (xad - rud * udx - (rdb - rud * rubd) * xab * rrab**2) + &
                 rvd * vdn * rrab * (xad - rvd * vdx - (rdb - rvd * rvbd) * xab * rrab**2)) * &
                gamd * rrad / wwd

          fby = gamb * (-cosb * yab * rrab**2 + rrab * (rub * uby + rvb * vby) / wwb) + &
                (ruc * ucn * rrab * (yac - ruc * ucy - (rbc - ruc * rubc) * yab * rrab**2) - &
                 rvc * vcn * rrab * (yac - rvc * vcy - (rbc - rvc * rvbc) * yab * rrab**2)) * &
                gamc * rrac / wwc + &
                (rud * udn * rrab * (yad - rud * udy - (rdb - rud * rubd) * yab * rrab**2) + &
                 rvd * vdn * rrab * (yad - rvd * vdy - (rdb - rvd * rvbd) * yab * rrab**2)) * &
                gamd * rrad / wwd

          fbz = gamb * (-cosb * zab * rrab**2 + rrab * (rub * ubz + rvb * vbz) / wwb) + &
                (ruc * ucn * rrab * (zac - ruc * ucz - (rbc - ruc * rubc) * zab * rrab**2) - &
                 rvc * vcn * rrab * (zac - rvc * vcz - (rbc - rvc * rvbc) * zab * rrab**2)) * &
                gamc * rrac / wwc + &
                (rud * udn * rrab * (zad - rud * udz - (rdb - rud * rubd) * zab * rrab**2) + &
                 rvd * vdn * rrab * (zad - rvd * vdz - (rdb - rvd * rvbd) * zab * rrab**2)) * &
                gamd * rrad / wwd

          fcx = gamc * (-cosc * xac * rrac**2 + rrac * (ruc * ucx + rvc * vcx) / wwc) + &
                (rud * udn * rrac * (xad - rud * udx - (rcd - rud * rucd) * xac * rrac**2) - &
                 rvd * vdn * rrac * (xad - rvd * vdx - (rcd - rvd * rvcd) * xac * rrac**2)) * &
                gamd * rrad / wwd + &
                (rub * ubn * rrac * (xab - rub * ubx - (rbc - rub * rucb) * xac * rrac**2) + &
                 rvb * vbn * rrac * (xab - rvb * vbx - (rbc - rvb * rvcb) * xac * rrac**2)) * &
                gamb * rrab / wwb

          fcy = gamc * (-cosc * yac * rrac**2 + rrac * (ruc * ucy + rvc * vcy) / wwc) + &
                (rud * udn * rrac * (yad - rud * udy - (rcd - rud * rucd) * yac * rrac**2) - &
                 rvd * vdn * rrac * (yad - rvd * vdy - (rcd - rvd * rvcd) * yac * rrac**2)) * &
                gamd * rrad / wwd + &
                (rub * ubn * rrac * (yab - rub * uby - (rbc - rub * rucb) * yac * rrac**2) + &
                 rvb * vbn * rrac * (yab - rvb * vby - (rbc - rvb * rvcb) * yac * rrac**2)) * &
                gamb * rrab / wwb

          fcz = gamc * (-cosc * zac * rrac**2 + rrac * (ruc * ucz + rvc * vcz) / wwc) + &
                (rud * udn * rrac * (zad - rud * udz - (rcd - rud * rucd) * zac * rrac**2) - &
                 rvd * vdn * rrac * (zad - rvd * vdz - (rcd - rvd * rvcd) * zac * rrac**2)) * &
                gamd * rrad / wwd + &
                (rub * ubn * rrac * (zab - rub * ubz - (rbc - rub * rucb) * zac * rrac**2) + &
                 rvb * vbn * rrac * (zab - rvb * vbz - (rbc - rvb * rvcb) * zac * rrac**2)) * &
                gamb * rrab / wwb

          fdx = gamd * (-cosd * xad * rrad**2 + rrad * (rud * udx + rvd * vdx) / wwd) + &
                (rub * ubn * rrad * (xab - rub * ubx - (rdb - rub * rudb) * xad * rrad**2) - &
                 rvb * vbn * rrad * (xab - rvb * vbx - (rdb - rvb * rvdb) * xad * rrad**2)) * &
                gamb * rrab / wwb + &
                (ruc * ucn * rrad * (xac - ruc * ucx - (rcd - ruc * rudc) * xad * rrad**2) + &
                 rvc * vcn * rrad * (xac - rvc * vcx - (rcd - rvc * rvdc) * xad * rrad**2)) * &
                gamc * rrac / wwc

          fdy = gamd * (-cosd * yad * rrad**2 + rrad * (rud * udy + rvd * vdy) / wwd) + &
                (rub * ubn * rrad * (yab - rub * uby - (rdb - rub * rudb) * yad * rrad**2) - &
                 rvb * vbn * rrad * (yab - rvb * vby - (rdb - rvb * rvdb) * yad * rrad**2)) * &
                gamb * rrab / wwb + &
                (ruc * ucn * rrad * (yac - ruc * ucy - (rcd - ruc * rudc) * yad * rrad**2) + &
                 rvc * vcn * rrad * (yac - rvc * vcy - (rcd - rvc * rvdc) * yad * rrad**2)) * &
                gamc * rrac / wwc

          fdz = gamd * (-cosd * zad * rrad**2 + rrad * (rud * udz + rvd * vdz) / wwd) + &
                (rub * ubn * rrad * (zab - rub * ubz - (rdb - rub * rudb) * zad * rrad**2) - &
                 rvb * vbn * rrad * (zab - rvb * vbz - (rdb - rvb * rvdb) * zad * rrad**2)) * &
                gamb * rrab / wwb + &
                (ruc * ucn * rrad * (zac - ruc * ucz - (rcd - ruc * rudc) * zad * rrad**2) + &
                 rvc * vcn * rrad * (zac - rvc * vcz - (rcd - rvc * rvdc) * zad * rrad**2)) * &
                gamc * rrac / wwc

          fax = -(fbx + fcx + fdx)
          fay = -(fby + fcy + fdy)
          faz = -(fbz + fcz + fdz)

        End If

        If (ia <= config%natms) Then

          ! inversion energy and virial (associated to the head atom)

          enginv = enginv + pterm
          virinv = virinv + vterm

          ! stress tensor calculation for inversion terms

          If (inversion%key(kk) == INVERSION_CALCITE) Then
            strs1 = strs1 + uuu * gamma * uux * uux
            strs2 = strs2 + uuu * gamma * uux * uuy
            strs3 = strs3 + uuu * gamma * uux * uuz
            strs5 = strs5 + uuu * gamma * uuy * uuy
            strs6 = strs6 + uuu * gamma * uuy * uuz
            strs9 = strs9 + uuu * gamma * uuz * uuz
          Else
            strs1 = strs1 + xab * fbx + xac * fcx + xad * fdx
            strs2 = strs2 + yab * fbx + yac * fcx + yad * fdx
            strs3 = strs3 + zab * fbx + zac * fcx + zad * fdx
            strs5 = strs5 + yab * fby + yac * fcy + yad * fdy
            strs6 = strs6 + yab * fbz + yac * fcz + yad * fdz
            strs9 = strs9 + zab * fbz + zac * fcz + zad * fdz
          End If

          config%parts(ia)%fxx = config%parts(ia)%fxx + fax
          config%parts(ia)%fyy = config%parts(ia)%fyy + fay
          config%parts(ia)%fzz = config%parts(ia)%fzz + faz

        End If

        If (ib <= config%natms) Then

          config%parts(ib)%fxx = config%parts(ib)%fxx + fbx
          config%parts(ib)%fyy = config%parts(ib)%fyy + fby
          config%parts(ib)%fzz = config%parts(ib)%fzz + fbz

        End If

        If (ic <= config%natms) Then

          config%parts(ic)%fxx = config%parts(ic)%fxx + fcx
          config%parts(ic)%fyy = config%parts(ic)%fyy + fcy
          config%parts(ic)%fzz = config%parts(ic)%fzz + fcz

        End If

        If (id <= config%natms) Then

          config%parts(id)%fxx = config%parts(id)%fxx + fdx
          config%parts(id)%fyy = config%parts(id)%fyy + fdy
          config%parts(id)%fzz = config%parts(id)%fzz + fdz

        End If

      End If
    End Do

    If (Mod(isw, 3) > 0) Then

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

      ! check for undefined potentials

      Call gcheck(comm, safe)
      If (.not. safe) Call error(449)

      ! global sum of inversion potential and virial

      buffer(1) = enginv
      buffer(2) = virinv
      Call gsum(comm, buffer(1:2))
      enginv = buffer(1)
      virinv = buffer(2)

    End If

    Deallocate (lunsafe, lstopt, Stat=fail(1))
    Deallocate (xdab, ydab, zdab, Stat=fail(2))
    Deallocate (xdac, ydac, zdac, Stat=fail(3))
    Deallocate (xdad, ydad, zdad, Stat=fail(4))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'inversions_forces deallocation failure'
      Call error(0, message)
    End If

  End Subroutine inversions_forces

  Subroutine inversions_table_read(invr_name, inversion, sites, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reading potential energy and force arrays
    ! from TABINV file (for inversion potentials & forces only)
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

    Type(inversions_type), Intent(InOut) :: inversion
    Character(Len=32),     Intent(In   ) :: invr_name(1:inversion%max_types)
    Type(site_type),       Intent(In   ) :: sites
    Type(comms_type),      Intent(InOut) :: comm

    Character(Len=200)         :: record
    Character(Len=256)         :: message, messages(4)
    Character(Len=32)          :: idinvr
    Character(Len=40)          :: word
    Character(Len=8)           :: atom1, atom2, atom3, atom4
    Integer                    :: fail(1:2), i, itinv, jtinv, jtpatm, katom1, katom2, katom3, &
                                  katom4, l, ngrid, rtinv
    Integer, Allocatable       :: read_type(:)
    Logical                    :: remake, safe
    Real(Kind=wp)              :: bufp0, bufv0, delpot, dgr2rad, dlrpot, ppp, rad2dgr, rdr, rrr, &
                                  rrr0, t1, t2, vk, vk1, vk2
    Real(Kind=wp), Allocatable :: bufpot(:), bufvir(:)

    If (comm%idnode == 0) Open (Unit=ntable, File='TABINV')

    ! skip header record

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Go To 100

    ! read mesh resolution not needed for inversion angle dependent
    ! potentials/forces as delpot=180/ngrid running from 0 to 180

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Go To 100

    i = Index(record, '#') ! replace hash as it may occur in
    If (i > 0) record(i:i) = ' ' ! TABINV if it's in .xvg format

    Call get_word(record, word)
    ngrid = Nint(word_2_real(word))

    delpot = 180.0_wp / Real(ngrid, wp)

    dlrpot = 180.0_wp / Real(inversion%bin_tab - 4, wp)

    ! check grid spacing

    safe = .false.
    If (Abs(delpot - dlrpot) < 1.0e-8_wp) Then
      safe = .true.
      delpot = dlrpot
    End If
    If (delpot > delth_max .and. (.not. safe)) Then
      Write (messages(1), '(a,1p,e15.7)') 'expected (maximum) angular increment : ', delth_max
      Write (messages(2), '(a,1p,e15.7)') 'TABINV file actual angular increment : ', delpot
      Write (messages(3), '(a,i10)') 'expected (minimum) number of grid points : ', inversion%bin_tab - 4
      Write (messages(4), '(a,i10)') 'TABINV file actual number of grid points : ', ngrid
      Call info(messages, 4, .true.)
      Call error(22)
    End If
    safe = .true.

    remake = .false.
    If (Abs(1.0_wp - (delpot / dlrpot)) > 1.0e-8_wp) Then
      remake = .true.
      rdr = 1.0_wp / delpot
      Write (message, '(a,i10)') 'TABINV arrays resized for mxgrid = ', inversion%bin_tab - 4
      Call info(message, .true.)
    End If

    ! compare grids dimensions

    If (ngrid < inversion%bin_tab - 4) Then
      Call warning(270, Real(ngrid, wp), Real(inversion%bin_tab - 4, wp), 0.0_wp)
      Call error(48)
    End If

    rad2dgr = 180.0_wp / pi
    dgr2rad = pi / 180.0_wp

    fail = 0
    Allocate (read_type(1:inversion%ltp(0)), Stat=fail(1))
    Allocate (bufpot(0:ngrid), bufvir(0:ngrid), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'error - inversions_table_read allocation failure'
      Call error(0, message)
    End If
    Call inversion%init_pot()

    read_type = 0 ! initialise read_type
    Do rtinv = 1, inversion%ltp(0)
      Call get_line(safe, ntable, record, comm)
      If (.not. safe) Go To 100

      Call get_line(safe, ntable, record, comm)
      If (.not. safe) Go To 100

      i = Index(record, '#') ! replace hash as it may occur in
      If (i > 0) record(i:i) = ' ' ! TABINV if it's in .xvg format

      Call get_word(record, atom1)
      Call get_word(record, atom2)
      Call get_word(record, atom3)
      Call get_word(record, atom4)

      katom1 = 0
      katom2 = 0
      katom3 = 0
      katom4 = 0

      Do jtpatm = 1, sites%ntype_atom
        If (atom1 == sites%unique_atom(jtpatm)) katom1 = jtpatm
        If (atom2 == sites%unique_atom(jtpatm)) katom2 = jtpatm
        If (atom3 == sites%unique_atom(jtpatm)) katom3 = jtpatm
        If (atom4 == sites%unique_atom(jtpatm)) katom4 = jtpatm
      End Do

      If (katom1 == 0 .or. katom2 == 0 .or. katom3 == 0 .or. katom4 == 0) Then
        Call info('****'//atom1//'***'//atom2//'***'//atom3//'***'//atom4// &
                  '**** entry in TABINV')
        Call error(91)
      End If

      ! Construct unique name for the tabulated inversion

      If (Min(katom2, katom3, katom4) == katom2) Then
        If (katom3 <= katom4) Then
          idinvr = atom1//atom2//atom3//atom4
        Else
          idinvr = atom1//atom2//atom4//atom3
        End If
      Else If (Min(katom2, katom3, katom4) == katom3) Then
        If (katom2 <= katom4) Then
          idinvr = atom1//atom3//atom2//atom4
        Else
          idinvr = atom1//atom3//atom4//atom2
        End If
      Else
        If (katom2 <= katom3) Then
          idinvr = atom1//atom4//atom2//atom3
        Else
          idinvr = atom1//atom4//atom3//atom2
        End If
      End If

      ! read potential arrays if potential is defined

      itinv = 0
      Do jtinv = 1, inversion%ltp(0)
        If (invr_name(jtinv) == idinvr) Then
          Do itinv = 1, inversion%max_types
            If (inversion%ltp(itinv) == jtinv) Exit
          End Do
          Exit
        End If
      End Do

      If (itinv == 0) Then ! All(invr_name /= idinvr)
        Call info('****'//atom1//'***'//atom2//'***'//atom3//'***'//atom4// &
                  '**** entry in TABINV')
        Call error(89)
      End If
      If (Any(read_type == jtinv)) Then
        Call info('****'//atom1//'***'//atom2//'***'//atom3//'***'//atom4// &
                  '**** entry in TABINV')
        Call error(172)
      Else
        read_type(jtinv) = jtinv
      End If

      ! read in potential & force arrays

      Do i = 0, 2
        bufpot(i) = 0.0_wp
        bufvir(i) = 0.0_wp
      End Do

      ! read in the zero and/or first & second data elements (potential & virial)

      If (comm%idnode == 0) Then
        rrr = 0.0_wp
        Read (Unit=ntable, Fmt=*, End=100, Err=100) rrr, bufp0, bufv0

        If (rrr > zero_plus) Then ! no zero element data => extrapolate to zero
          If (Abs((rrr - delpot) / delpot) > 1.0e-8_wp) Then
            safe = .false.
            Write (messages(1), '(a,1p,e15.7)') 'TABINV stated  angular increment : ', delpot
            Write (messages(2), '(a,1p,e15.7)') 'TABINV read-in angular increment : ', rrr
            Call info(messages, 2, .true.)
          End If

          bufpot(1) = bufp0
          bufvir(1) = bufv0
          rrr0 = rrr

          Read (Unit=ntable, Fmt=*, End=100, Err=100) rrr, bufp0, bufv0

          If (Abs((rrr - rrr0 - delpot) / delpot) > 1.0e-8_wp) Then
            safe = .false.
            Write (messages(1), '(a,1p,e15.7)') 'TABINV stated  angular increment : ', delpot
            Write (messages(2), '(a,1p,e15.7)') 'TABINV read-in angular increment : ', rrr
            Call info(messages, 2, .true.)
          End If

          bufpot(2) = bufp0
          bufvir(2) = bufv0

          ! linear extrapolation for grid point 0 at distances close to 0

          bufpot(0) = 2.0_wp * bufpot(1) - bufpot(2)
          bufvir(0) = (2.0_wp * bufvir(1) - 0.5_wp * bufvir(2)) / dlrpot
        Else ! zero element data found => read in the first element for checking delr
          bufpot(0) = bufp0
          bufvir(0) = bufv0

          Read (Unit=ntable, Fmt=*, End=100, Err=100) rrr, bufp0, bufv0

          If (Abs((rrr - delpot) / delpot) > 1.0e-8_wp) Then
            safe = .false.
            Write (messages(1), '(a,1p,e15.7)') 'TABINV stated  angular increment : ', delpot
            Write (messages(2), '(a,1p,e15.7)') 'TABINV read-in angular increment : ', rrr
            Call info(messages, 2, .true.)
          End If

          bufpot(1) = bufp0
          bufvir(1) = bufv0

          Read (Unit=ntable, Fmt=*, End=100, Err=100) rrr, bufp0, bufv0

          bufpot(2) = bufp0
          bufvir(2) = bufv0
        End If
      End If

      ! read in potential & force arrays

      Do i = 3, ngrid
        If (comm%idnode == 0) Then
          Read (Unit=ntable, Fmt=*, End=100, Err=100) rrr, bufpot(i), bufvir(i)
        Else
          bufpot(i) = 0.0_wp
          bufvir(i) = 0.0_wp
        End If
      End Do

      Call gbcast(comm, bufpot, 0)
      Call gbcast(comm, bufvir, 0)

      ! reconstruct arrays using 3pt interpolation

      If (remake) Then
        Do i = 1, inversion%bin_tab - 4
          rrr = Real(i, wp) * dlrpot
          l = Int(rrr * rdr)
          ppp = rrr * rdr - Real(l, wp)

          vk = bufpot(l)

          ! linear extrapolation for the grid points just beyond the cutoff

          If (l + 2 > ngrid) Then
            If (l + 1 > ngrid) Then
              vk1 = 2.0_wp * bufpot(l) - bufpot(l - 1)
              vk2 = 2.0_wp * vk1 - bufpot(l)
            Else
              vk1 = bufpot(l + 1)
              vk2 = 2.0_wp * bufpot(l + 1) - bufpot(l)
            End If
          Else
            vk1 = bufpot(l + 1)
            vk2 = bufpot(l + 2)
          End If

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)
          inversion%tab_potential(i, jtinv) = t1 + (t2 - t1) * ppp * 0.5_wp
          inversion%tab_potential(i, jtinv) = inversion%tab_potential(i, jtinv) * engunit ! convert to internal units

          vk = bufvir(l)

          ! linear extrapolation for the grid points just beyond the cutoff

          If (l + 2 > ngrid) Then
            If (l + 1 > ngrid) Then
              vk1 = 2.0_wp * bufvir(l) - bufvir(l - 1)
              vk2 = 2.0_wp * vk1 - bufvir(l)
            Else
              vk1 = bufvir(l + 1)
              vk2 = 2.0_wp * bufvir(l + 1) - bufvir(l)
            End If
          Else
            vk1 = bufvir(l + 1)
            vk2 = bufvir(l + 2)
          End If

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)
          inversion%tab_force(i, jtinv) = t1 + (t2 - t1) * ppp * 0.5_wp
          inversion%tab_force(i, jtinv) = inversion%tab_force(i, jtinv) * engunit * rad2dgr ! convert to internal units
        End Do

        inversion%tab_force(-1, jtinv) = rad2dgr / dlrpot
      Else
        Do i = 1, inversion%bin_tab - 4
          inversion%tab_potential(i, jtinv) = bufpot(i) * engunit ! convert to internal units
          inversion%tab_force(i, jtinv) = bufvir(i) * engunit * rad2dgr ! convert to internal units
        End Do

        ! linear extrapolation for the grid point just beyond the cutoff

        inversion%tab_potential(inversion%bin_tab - 3, jtinv) = &
          2.0_wp * inversion%tab_potential(inversion%bin_tab - 4, jtinv) - &
          inversion%tab_potential(inversion%bin_tab - 5, jtinv)
        inversion%tab_force(inversion%bin_tab - 3, jtinv) = &
          2.0_wp * inversion%tab_force(inversion%bin_tab - 4, jtinv) - &
          inversion%tab_force(inversion%bin_tab - 5, jtinv)

        inversion%tab_force(-1, jtinv) = rad2dgr / delpot
      End If

      ! grid point at 0 and linear extrapolation for the grid point at inversion%bin_tab-2

      inversion%tab_potential(0, jtinv) = bufpot(0)
      inversion%tab_force(0, jtinv) = bufvir(0)

      inversion%tab_potential(inversion%bin_tab - 2, jtinv) = &
        2.0_wp * inversion%tab_potential(inversion%bin_tab - 3, jtinv) - &
        inversion%tab_potential(inversion%bin_tab - 4, jtinv)
      inversion%tab_force(inversion%bin_tab - 2, jtinv) = &
        2.0_wp * inversion%tab_force(inversion%bin_tab - 3, jtinv) - &
        inversion%tab_force(inversion%bin_tab - 4, jtinv)
    End Do

    If (comm%idnode == 0) Then
      Close (Unit=ntable)
    End If
    Call info('potential tables read from TABINV file', .true.)

    ! Break if not safe

    Call gcheck(comm, safe)
    If (.not. safe) Call error(22)

    Deallocate (read_type, Stat=fail(1))
    Deallocate (bufpot, bufvir, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'error - inversions_table_read deallocation failure'
      Call error(0, message)
    End If

    Return

    ! end of file error exit

    100 Continue

    If (comm%idnode == 0) Close (Unit=ntable)
    Call error(24)

  End Subroutine inversions_table_read

  Subroutine cleanup(inversion)
    Type(inversions_type) :: inversion

    If (Allocated(inversion%num)) Then
      Deallocate (inversion%num)
    End If
    If (Allocated(inversion%key)) Then
      Deallocate (inversion%key)
    End If

    If (Allocated(inversion%lst)) Then
      Deallocate (inversion%lst)
    End If
    If (Allocated(inversion%list)) Then
      Deallocate (inversion%list)
    End If
    If (Allocated(inversion%legend)) Then
      Deallocate (inversion%legend)
    End If

    If (Allocated(inversion%param)) Then
      Deallocate (inversion%param)
    End If

    If (Allocated(inversion%ltp)) Then
      Deallocate (inversion%ltp)
    End If
    If (Allocated(inversion%tab_potential)) Then
      Deallocate (inversion%tab_potential)
    End If
    If (Allocated(inversion%tab_force)) Then
      Deallocate (inversion%tab_force)
    End If
  End Subroutine cleanup
End Module inversions
