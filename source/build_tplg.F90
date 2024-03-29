Module build_tplg
  ! SETUP MODULES

  Use angles,          Only: ANGLE_NULL,&
                             angles_type
  Use bonds,           Only: BOND_NULL,&
                             bonds_type
  Use build_excl,      Only: add_exclusion
  Use comms,           Only: comms_type,&
                             gcheck,&
                             gmax
  ! CONFIG MODULE
  Use configuration,   Only: configuration_type
  ! INTERACTION MODULES
  Use dihedrals,       Only: dihedrals_type
  Use errors_warnings, Only: error,&
                             warning
  Use inversions,      Only: inversions_type
  ! MULTIPOLES MODULE
  Use kinds,           Only: wi,&
                             wp
  Use mpole,           Only: mpole_type
  Use numerics,        Only: local_index,&
                             shellsort

  Implicit None

  Private
  Public :: build_tplg_intra

Contains

  Subroutine build_tplg_intra(max_exclude, bond, angle, dihedral, inversion, mpoles, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for constructing the topology related list
    ! of neighbours for the MD system mapped onto this node.  It is presumed
    ! that constraint bonds are put on top of chemical bonds!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer(Kind=wi),         Intent(In   ) :: max_exclude
    Type(bonds_type),         Intent(In   ) :: bond
    Type(angles_type),        Intent(In   ) :: angle
    Type(dihedrals_type),     Intent(In   ) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(mpole_type),         Intent(InOut) :: mpoles
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Integer :: i, ia, ia0, ib, ib0, ibig, ic, ic0, id, id0, j
    Logical :: safe

    ! variables for array bound checking

    ibig = 0
    safe = .true.

    ! include sites on basis of chemical bonds

    Do i = 1, bond%n_types
      If (bond%key(bond%list(0, i)) /= BOND_NULL) Then
        ia = bond%list(1, i)
        ib = bond%list(2, i)

        ia0 = local_index(ia, config%nlast, config%lsi, config%lsa)
        ib0 = local_index(ib, config%nlast, config%lsi, config%lsa)

        If (ia0 > config%natms) ia0 = 0
        If (ib0 > config%natms) ib0 = 0

        ! add atoms to topology list

        If (ia0 > 0) Call add_exclusion(safe, ia0, ib, ibig, mpoles%ltp)
        If (ib0 > 0) Call add_exclusion(safe, ib0, ia, ibig, mpoles%ltp)
      End If
    End Do

    ! include sites on basis of bond angles

    Do i = 1, angle%n_types
      If (angle%key(angle%list(0, i)) /= ANGLE_NULL) Then
        ia = angle%list(1, i)
        ib = angle%list(2, i)
        ic = angle%list(3, i)

        ia0 = local_index(ia, config%nlast, config%lsi, config%lsa)
        ib0 = local_index(ib, config%nlast, config%lsi, config%lsa)
        ic0 = local_index(ic, config%nlast, config%lsi, config%lsa)

        If (ia0 > config%natms) ia0 = 0
        If (ib0 > config%natms) ib0 = 0
        If (ic0 > config%natms) ic0 = 0

        ! add atoms to topology list

        If (ia0 > 0) Then ! ia : ib - ic neighbours
          Call add_exclusion(safe, ia0, ib, ibig, mpoles%ltp)
          Call add_exclusion(safe, ia0, ic, ibig, mpoles%ltp)
        End If

        If (ib0 > 0) Then ! ib : ia - ic neighbours
          Call add_exclusion(safe, ib0, ia, ibig, mpoles%ltp)
          Call add_exclusion(safe, ib0, ic, ibig, mpoles%ltp)
        End If

        If (ic0 > 0) Then ! ic : ia - ib neighbours
          Call add_exclusion(safe, ic0, ia, ibig, mpoles%ltp)
          Call add_exclusion(safe, ic0, ib, ibig, mpoles%ltp)
        End If
      End If
    End Do

    ! include sites on basis of dihedral angles

    Do i = 1, dihedral%n_types
      ia = dihedral%list(1, i)
      ib = dihedral%list(2, i)
      ic = dihedral%list(3, i)
      id = dihedral%list(4, i)

      ia0 = local_index(ia, config%nlast, config%lsi, config%lsa)
      ib0 = local_index(ib, config%nlast, config%lsi, config%lsa)
      ic0 = local_index(ic, config%nlast, config%lsi, config%lsa)
      id0 = local_index(id, config%nlast, config%lsi, config%lsa)

      If (ia0 > config%natms) ia0 = 0
      If (ib0 > config%natms) ib0 = 0
      If (ic0 > config%natms) ic0 = 0
      If (id0 > config%natms) id0 = 0

      ! add atoms to topology list

      If (ia0 > 0) Then ! ia : ib - ic neighbours
        Call add_exclusion(safe, ia0, ib, ibig, mpoles%ltp)
        Call add_exclusion(safe, ia0, ic, ibig, mpoles%ltp)
      End If

      If (ib0 > 0) Then ! ib : ia - ic neighbours
        Call add_exclusion(safe, ib0, ia, ibig, mpoles%ltp)
        Call add_exclusion(safe, ib0, ic, ibig, mpoles%ltp)
      End If

      If (ic0 > 0) Then ! ic : ib - id neighbours
        Call add_exclusion(safe, ic0, ib, ibig, mpoles%ltp)
        Call add_exclusion(safe, ic0, id, ibig, mpoles%ltp)
      End If

      If (id0 > 0) Then ! id : ib - ic neighbours
        Call add_exclusion(safe, id0, ib, ibig, mpoles%ltp)
        Call add_exclusion(safe, id0, ic, ibig, mpoles%ltp)
      End If
    End Do

    ! include sites on basis of inversion potentials

    Do i = 1, inversion%n_types
      ia = inversion%list(1, i)
      ib = inversion%list(2, i)
      ic = inversion%list(3, i)
      id = inversion%list(4, i)

      ia0 = local_index(ia, config%nlast, config%lsi, config%lsa)
      ib0 = local_index(ib, config%nlast, config%lsi, config%lsa)
      ic0 = local_index(ic, config%nlast, config%lsi, config%lsa)
      id0 = local_index(id, config%nlast, config%lsi, config%lsa)

      If (ia0 > config%natms) ia0 = 0
      If (ib0 > config%natms) ib0 = 0
      If (ic0 > config%natms) ic0 = 0
      If (id0 > config%natms) id0 = 0

      ! add atoms to topology list

      If (ia0 > 0) Then ! ia : ib - ic - id neighbours
        Call add_exclusion(safe, ia0, ib, ibig, mpoles%ltp)
        Call add_exclusion(safe, ia0, ic, ibig, mpoles%ltp)
        Call add_exclusion(safe, ia0, id, ibig, mpoles%ltp)
      End If

      If (ib0 > 0) Then ! ib : ia - ic - id neighbours
        Call add_exclusion(safe, ib0, ia, ibig, mpoles%ltp)
        Call add_exclusion(safe, ib0, ic, ibig, mpoles%ltp)
        Call add_exclusion(safe, ib0, id, ibig, mpoles%ltp)
      End If

      If (ic0 > 0) Then ! ic : ia - ib - id neighbours
        Call add_exclusion(safe, ic0, ia, ibig, mpoles%ltp)
        Call add_exclusion(safe, ic0, ib, ibig, mpoles%ltp)
        Call add_exclusion(safe, ic0, id, ibig, mpoles%ltp)
      End If

      If (id0 > 0) Then ! id : ia - ib - ic neighbours
        Call add_exclusion(safe, id0, ia, ibig, mpoles%ltp)
        Call add_exclusion(safe, id0, ib, ibig, mpoles%ltp)
        Call add_exclusion(safe, id0, ic, ibig, mpoles%ltp)
      End If
    End Do

    ! check for exceeded array bounds

    Call gcheck(comm, safe)
    If (.not. safe) Then
      Call gmax(comm, ibig)
      Call warning(250, Real(ibig, wp), Real(max_exclude, wp), 0.0_wp)
      Call error(65)
    End If

    ! sort mpoles%ltp

    Do i = 1, config%natms
      j = mpoles%ltp(0, i)
      If (j > 0) Call shellsort(j, mpoles%ltp(1:j, i))
    End Do
  End Subroutine build_tplg_intra
End Module build_tplg

