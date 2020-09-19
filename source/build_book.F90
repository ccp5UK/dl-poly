Module build_book
  Use angles,          Only: angles_type
  Use bonds,           Only: bonds_type
  Use comms,           Only: comms_type,&
                             gcheck,&
                             gmax
  Use configuration,   Only: configuration_type
  Use constraints,     Only: constraints_type
  Use core_shell,      Only: core_shell_type
  Use dihedrals,       Only: dihedrals_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use ffield,          Only: report_topology
  Use flow_control,    Only: flow_type
  Use inversions,      Only: inversions_type
  Use kinds,           Only: wp
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: local_index
  Use pmf,             Only: pmf_type
  Use rigid_bodies,    Only: rigid_bodies_coms,&
                             rigid_bodies_setup,&
                             rigid_bodies_tags,&
                             rigid_bodies_type,&
                             rigid_bodies_widths
  Use shared_units,    Only: pass_shared_units,&
                             tag_legend
  Use site,            Only: site_type
  Use tethers,         Only: tethers_type

  Implicit None

  Private

  Public :: build_book_intra
  Public :: compress_book_intra

Contains

  Subroutine build_book_intra(l_str, l_top, lsim, flow, &
                              cshell, cons, pmf, bond, angle, dihedral, &
                              inversion, tether, neigh, sites, rigid, domain, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for setting up the domain distributed bookkeeping
    ! of intra-like interactions: core-shell, bond constraints, PMF
    ! constraints, RBs, tethered atoms, chemical bonds, valence angles,
    ! torsion and improper torsion angles, and inversion angles
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,                  Intent(In   ) :: l_str, l_top, lsim
    Type(flow_type),          Intent(InOut) :: flow
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(Inout) :: cons
    Type(pmf_type),           Intent(Inout) :: pmf
    Type(bonds_type),         Intent(InOut) :: bond
    Type(angles_type),        Intent(InOut) :: angle
    Type(dihedrals_type),     Intent(InOut) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(tethers_type),       Intent(InOut) :: tether
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(site_type),          Intent(InOut) :: sites
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)                 :: message, messages(7)
    Integer                            :: fail(1:2), i, iangle, iat0, iatm, ibonds, iconst, &
                                          idihed, iinver, imols, ipmf, irigid, ishels, isite, &
                                          iteths, itmols, itmp(1:9), j, jangle, jat0, jatm, &
                                          jbonds, jconst, jdihed, jinver, jrigid, jshels, jteths, &
                                          jtmp(1:9), kangle, kat0, katm, kbonds, kconst, kdihed, &
                                          kinver, krigid, kshels, kteths, langle, lat0, latm, &
                                          lbonds, lconst, ldihed, linver, lpmf, lrigid, lshels, &
                                          lteths, mat0, matm, mrigid, mshels, nat0, natm, neatm, &
                                          nlapm, nsatm
    Integer, Allocatable, Dimension(:) :: i1pmf, i1pmf0, i2pmf, i2pmf0, irgd, irgd0, iwrk
    Logical                            :: go, safe(1:11)
    Real(Kind=wp)                      :: tmp

    fail = 0
    Allocate (iwrk(1:config%mxatms), Stat=fail(1))
    If (rigid%on) Then
      Allocate (irgd(1:rigid%max_list), &
                irgd0(1:rigid%max_list), Stat=fail(2))
    End If
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'build_book_intra allocation failure'
      Call error(0, message)
    End If

    If (.not. (flow%newjob_build_book .or. lsim)) Then
      Call init_intra(cshell, cons, pmf, bond, angle, dihedral, inversion, tether, neigh, rigid)
    End If

    ! Initialise safety flags

    safe = .true.

    ! Initialise bookkeeping indices

    isite = 0

    ! "i" stands for excess number of intra-like units
    !
    ! "j" stands for running index of intra-like unit
    !
    ! "k" stands for last index of intra-like unit
    !     on last molecule of previous molecule type

    ishels = 0
    jshels = 0
    kshels = 0

    iconst = 0
    jconst = 0
    kconst = 0

    ! No 'jpmf' and 'kpmf' needed since PMF is defined on one and only one
    ! molecular type and therefore 'pmf%ntpmf' is enough and is used locally as 'jpmf'

    ipmf = 0

    irigid = 0
    jrigid = 0
    krigid = 0

    iteths = 0
    jteths = 0
    kteths = 0

    ibonds = 0
    jbonds = 0
    kbonds = 0

    iangle = 0
    jangle = 0
    kangle = 0

    idihed = 0
    jdihed = 0
    kdihed = 0

    iinver = 0
    jinver = 0
    kinver = 0

    ! global atom counter

    nsatm = 0

    ! loop over molecule types in the system

    Do itmols = 1, sites%ntype_mol

      ! loop over molecules of this type

      Do imols = 1, sites%num_mols(itmols)

        ! last atom in the molecule

        neatm = nsatm + sites%num_site(itmols)

        ! number of local atoms on this molecule

        nlapm = 0

        ! From the first till the last atom of this molecule, get the number
        ! of localised atoms on this node and save their local_index in iwrk

        Do iatm = nsatm + 1, neatm
          iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
          If (iat0 > config%natms) iat0 = 0

          If (iat0 > 0) Then
            nlapm = nlapm + 1
            iwrk(nlapm) = iat0
          End If
        End Do

        ! If there are atoms of this molecule on this node, get global indices and
        ! corresponding local_indices of the specific intra-like unit.  If any of
        ! the local indices exists increase the local number of these intra-like
        ! units (local = belongs to this node).  If it's safe to proceed (array
        ! bounds check) record the global number of the intra-like unit and the
        ! global atomic indices of its atoms.  Tag the local numbers of the
        ! intra-like units and count them for each local atom (on domain) that is
        ! a member of them in the intra-like unit's legend array.

        If (nlapm > 0) Then

          ! Construct core-shell list

          Do lshels = 1, cshell%numshl(itmols)
            iatm = cshell%lstshl(1, lshels + kshels) + isite
            jatm = cshell%lstshl(2, lshels + kshels) + isite

            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)

            If (iat0 > config%natms) iat0 = 0
            If (jat0 > config%natms) jat0 = 0

            If (iat0 > 0 .or. jat0 > 0) Then
              jshels = jshels + 1

              If (jshels <= cshell%mxshl) Then
                cshell%listshl(0, jshels) = lshels + kshels
                cshell%listshl(1, jshels) = iatm
                cshell%listshl(2, jshels) = jatm

                If (iat0 > 0) Then
                  Call tag_legend(safe(1), iat0, jshels, cshell%legshl, cshell%mxfshl)
                  If (cshell%legshl(cshell%mxfshl, iat0) > 0) Then
                    Call warning('too many core-shell type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      cshell%mxfshl - 1 + cshell%legshl(cshell%mxfshl, iat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', cshell%mxfshl - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', cshell%lstshl(1, lshels + kshels)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lshels
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (jat0 > 0) Then
                  Call tag_legend(safe(1), jat0, -jshels, cshell%legshl, cshell%mxfshl)
                  If (cshell%legshl(cshell%mxfshl, jat0) > 0) Then
                    Call warning('too many core-shell type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      cshell%mxfshl - 1 + cshell%legshl(cshell%mxfshl, jat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', cshell%mxfshl - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', jatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', cshell%lstshl(2, lshels + kshels)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lshels
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If
              Else
                ishels = ishels + 1
                safe(2) = .false.
              End If
            End If
          End Do

          ! Construct constraint bond list

          Do lconst = 1, cons%numcon(itmols)
            iatm = cons%lstcon(1, lconst + kconst) + isite
            jatm = cons%lstcon(2, lconst + kconst) + isite

            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)

            If (iat0 > config%natms) iat0 = 0
            If (jat0 > config%natms) jat0 = 0

            If (iat0 > 0 .or. jat0 > 0) Then
              jconst = jconst + 1
              If (jconst <= cons%mxcons) Then
                cons%listcon(0, jconst) = lconst + kconst
                cons%listcon(1, jconst) = iatm
                cons%listcon(2, jconst) = jatm

                If (iat0 > 0) Then
                  Call tag_legend(safe(1), iat0, jconst, cons%legcon, cons%mxfcon)
                  If (cons%legcon(cons%mxfcon, iat0) > 0) Then
                    Call warning('too many constraint type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      cons%mxfcon - 1 + cons%legcon(cons%mxfcon, iat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', cons%mxfcon - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', cons%lstcon(1, lconst + kconst)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lconst
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (jat0 > 0) Then
                  Call tag_legend(safe(1), jat0, jconst, cons%legcon, cons%mxfcon)
                  If (cons%legcon(cons%mxfcon, jat0) > 0) Then
                    Call warning('too many constraint type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      cons%mxfcon - 1 + cons%legcon(cons%mxfcon, jat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', cons%mxfcon - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', jatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', cons%lstcon(2, lconst + kconst)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lconst
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If
              Else
                iconst = iconst + 1
                safe(3) = .false.
              End If
            End If
          End Do

          ! Construct PMF bond list
          ! Note: This is executed for only one given molecular
          !       type as only one PMF type per MD system is allowed.

          Do lpmf = 1, pmf%numpmf(itmols) ! pmf%numpmf can only be 1 or 0, so the 'Do' loop is used as an 'If' condition
            Allocate (i1pmf(1:pmf%mxtpmf(1)), i1pmf0(1:pmf%mxtpmf(1)), i2pmf(1:pmf%mxtpmf(2)), i2pmf0(1:pmf%mxtpmf(2)), &
                      Stat=fail(1))
            If (fail(1) > 0) Then
              Write (message, '(a)') 'build_book_intra PMF allocation failure'
              Call error(0, message)
            End If

            i1pmf = 0; i1pmf0 = 0
            Do i = 1, pmf%mxtpmf(1)
              i1pmf(i) = pmf%lstpmf(i, 1) + isite
              i1pmf0(i) = local_index(i1pmf(i), config%nlast, config%lsi, config%lsa)
              If (i1pmf0(i) > config%natms) i1pmf0(i) = 0
            End Do

            i2pmf = 0; i2pmf0 = 0
            Do i = 1, pmf%mxtpmf(2)
              i2pmf(i) = pmf%lstpmf(i, 2) + isite
              i2pmf0(i) = local_index(i2pmf(i), config%nlast, config%lsi, config%lsa)
              If (i2pmf0(i) > config%natms) i2pmf0(i) = 0
            End Do

            If (Any(i1pmf0 > 0) .or. Any(i2pmf0 > 0)) Then
              pmf%ntpmf = pmf%ntpmf + 1
              If (pmf%ntpmf <= pmf%mxpmf) Then

                ! This holds the global PMF index

                pmf%listpmf(0, 1, pmf%ntpmf) = imols

                ! For presence of : PMF unit 1 only - this holds 1
                !                   PMF unit 2 only - this holds 2
                !                   both units 1&2  - this holds 3
                ! It CANNOT and MUST NOT hold ZERO

                pmf%listpmf(0, 2, pmf%ntpmf) = 0

                Do i = 1, pmf%mxtpmf(1)
                  pmf%listpmf(i, 1, pmf%ntpmf) = i1pmf(i)
                  If (i1pmf0(i) > 0) Then
                    Call tag_legend(safe(1), i1pmf0(i), pmf%ntpmf, pmf%legpmf, pmf%mxfpmf)
                    If (pmf%legpmf(pmf%mxfpmf, i1pmf0(i)) > 0) Then
                      Call warning('too many PMF type neighbours')
                      Write (messages(1), '(a,i0)') 'requiring a neigh%list length of: ', &
                        pmf%mxfpmf - 1 + pmf%legpmf(pmf%mxfpmf, i1pmf0(i))
                      Write (messages(2), '(a,i0)') 'but maximum length allowed: ', pmf%mxfpmf - 1
                      Write (messages(3), '(a,i0)') 'for particle (global ID #): ', i1pmf(i)
                      Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', pmf%lstpmf(i, 1)
                      Write (messages(5), '(a,i0)') 'of PMF unit  (1 or 2 only): ', 1
                      Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                      Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                      Call info(messages, 7)
                    End If
                  End If
                End Do
                If (Any(i1pmf0 > 0)) pmf%listpmf(0, 2, pmf%ntpmf) = pmf%listpmf(0, 2, pmf%ntpmf) + 1

                Do i = 1, pmf%mxtpmf(2)
                  pmf%listpmf(i, 2, pmf%ntpmf) = i2pmf(i)
                  If (i2pmf0(i) > 0) Then
                    Call tag_legend(safe(1), i2pmf0(i), pmf%ntpmf, pmf%legpmf, pmf%mxfpmf)
                    If (pmf%legpmf(pmf%mxfpmf, i2pmf0(i)) > 0) Then
                      Call warning('too many PMF type neighbours')
                      Write (messages(1), '(a,i0)') 'requiring a neigh%list length of: ', &
                        pmf%mxfpmf - 1 + pmf%legpmf(pmf%mxfpmf, i2pmf0(i))
                      Write (messages(2), '(a,i0)') 'but maximum length allowed: ', pmf%mxfpmf - 1
                      Write (messages(3), '(a,i0)') 'for particle (global ID #): ', i2pmf(i)
                      Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', pmf%lstpmf(i, 2)
                      Write (messages(5), '(a,i0)') 'of PMF unit  (1 or 2 only): ', 2
                      Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                      Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                      Call info(messages, 7)
                    End If
                  End If
                End Do
                If (Any(i2pmf0 > 0)) pmf%listpmf(0, 2, pmf%ntpmf) = pmf%listpmf(0, 2, pmf%ntpmf) + 2

              Else
                ipmf = ipmf + 1
                safe(4) = .false.
              End If
            End If

            Deallocate (i1pmf, i1pmf0, i2pmf, i2pmf0, Stat=fail(1))
            If (fail(1) > 0) Then
              Write (message, '(a)') 'build_book_intra PMF deallocation failure'
              Call error(0, message)
            End If
          End Do

          ! Construct RBs list

          Do lrigid = 1, rigid%num(itmols)
            mrigid = rigid%lst(0, lrigid + krigid)

            irgd = 0; irgd0 = 0
            Do irigid = 1, mrigid
              irgd(irigid) = rigid%lst(irigid, lrigid + krigid) + isite
              irgd0(irigid) = local_index(irgd(irigid), config%nlast, config%lsi, config%lsa)
              If (irgd0(irigid) > config%natms) irgd0(irigid) = 0
            End Do

            If (Any(irgd0 > 0)) Then
              jrigid = jrigid + 1
              If (jrigid <= rigid%max_rigid) Then
                rigid%list(-1, jrigid) = mrigid
                rigid%list(0, jrigid) = lrigid + krigid
                Do irigid = 1, mrigid
                  rigid%list(irigid, jrigid) = irgd(irigid)
                  If (irgd0(irigid) > 0) Then
                    Call tag_legend(safe(1), irgd0(irigid), jrigid, rigid%legend, rigid%max_frozen)
                    If (rigid%legend(rigid%max_frozen, irgd0(irigid)) > 0) Then
                      Call warning('too many RB type neighbours')
                      Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                        rigid%max_frozen - 1 + rigid%legend(rigid%max_frozen, irgd0(irigid))
                      Write (messages(2), '(a,i0)') 'but maximum length allowed: ', rigid%max_frozen - 1
                      Write (messages(3), '(a,i0)') 'for particle (global ID #): ', irgd(irigid)
                      Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', rigid%lst(irigid, lrigid + krigid)
                      Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lrigid
                      Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                      Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                      Call info(messages, 7)
                    End If
                  End If
                End Do
              Else
                irigid = irigid + 1
                safe(5) = .false.
              End If
            End If
          End Do

          ! Construct tethered atoms interaction list

          Do lteths = 1, tether%numteth(itmols)
            iatm = tether%lsttet(lteths + kteths) + isite
            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            If (iat0 > config%natms) iat0 = 0

            If (iat0 > 0) Then
              jteths = jteths + 1
              If (jteths <= tether%mxteth) Then
                tether%listtet(0, jteths) = lteths + kteths
                tether%listtet(1, jteths) = iatm

                Call tag_legend(safe(1), iat0, jteths, tether%legtet, tether%mxftet)
                If (tether%legtet(tether%mxftet, iat0) > 0) Then
                  Call warning('too many tether type neighbours')
                  Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                    tether%mxftet - 1 + tether%legtet(tether%mxftet, iat0)
                  Write (messages(2), '(a,i0)') 'but maximum length allowed: ', tether%mxftet - 1
                  Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                  Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', tether%lsttet(lteths + kteths)
                  Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lteths
                  Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                  Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                  Call info(messages, 7)
                End If
              Else
                iteths = iteths + 1
                safe(6) = .false.
              End If
            End If
          End Do

          ! Construct chemical bond interaction list

          Do lbonds = 1, bond%num(itmols)
            iatm = bond%lst(1, lbonds + kbonds) + isite
            jatm = bond%lst(2, lbonds + kbonds) + isite

            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)

            If (iat0 > config%natms) iat0 = 0
            If (jat0 > config%natms) jat0 = 0

            If (iat0 > 0 .or. jat0 > 0) Then
              jbonds = jbonds + 1
              If (jbonds <= bond%max_bonds) Then
                bond%list(0, jbonds) = lbonds + kbonds
                bond%list(1, jbonds) = iatm
                bond%list(2, jbonds) = jatm

                If (iat0 > 0) Then
                  Call tag_legend(safe(1), iat0, jbonds, bond%legend, bond%max_legend)
                  If (bond%legend(bond%max_legend, iat0) > 0) Then
                    Call warning('too many bond type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      bond%max_legend - 1 + bond%legend(bond%max_legend, iat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', bond%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', bond%lst(1, lbonds + kbonds)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lbonds
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (jat0 > 0) Then
                  Call tag_legend(safe(1), jat0, jbonds, bond%legend, bond%max_legend)
                  If (bond%legend(bond%max_legend, jat0) > 0) Then
                    Call warning('too many bond type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      bond%max_legend - 1 + bond%legend(bond%max_legend, jat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', bond%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', bond%lst(2, lbonds + kbonds)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', lbonds
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If
              Else
                ibonds = ibonds + 1
                safe(7) = .false.
              End If
            End If

          End Do

          ! Construct valence angle interaction list

          Do langle = 1, angle%num(itmols)
            iatm = angle%lst(1, langle + kangle) + isite
            jatm = angle%lst(2, langle + kangle) + isite
            katm = angle%lst(3, langle + kangle) + isite

            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)
            kat0 = local_index(katm, config%nlast, config%lsi, config%lsa)

            If (iat0 > config%natms) iat0 = 0
            If (jat0 > config%natms) jat0 = 0
            If (kat0 > config%natms) kat0 = 0

            If (iat0 > 0 .or. jat0 > 0 .or. kat0 > 0) Then
              jangle = jangle + 1
              If (jangle <= angle%max_angles) Then
                angle%list(0, jangle) = langle + kangle
                angle%list(1, jangle) = iatm
                angle%list(2, jangle) = jatm
                angle%list(3, jangle) = katm

                If (iat0 > 0) Then
                  Call tag_legend(safe(1), iat0, jangle, angle%legend, angle%max_legend)
                  If (angle%legend(angle%max_legend, iat0) > 0) Then
                    Call warning('too many angle type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      angle%max_legend - 1 + angle%legend(angle%max_legend, iat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', angle%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', angle%lst(1, langle + kangle)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', langle
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (jat0 > 0) Then
                  Call tag_legend(safe(1), jat0, jangle, angle%legend, angle%max_legend)
                  If (angle%legend(angle%max_legend, jat0) > 0) Then
                    Call warning('too many angle type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      angle%max_legend - 1 + angle%legend(angle%max_legend, jat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', angle%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', angle%lst(2, langle + kangle)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', langle
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (kat0 > 0) Then
                  Call tag_legend(safe(1), kat0, jangle, angle%legend, angle%max_legend)
                  If (angle%legend(angle%max_legend, kat0) > 0) Then
                    Call warning('too many angle type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      angle%max_legend - 1 + angle%legend(angle%max_legend, kat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', angle%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', angle%lst(3, langle + kangle)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', langle
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If
              Else
                iangle = iangle + 1
                safe(8) = .false.
              End If
            End If
          End Do

          ! Construct dihedral angle interaction list

          Do ldihed = 1, dihedral%num(itmols)
            iatm = dihedral%lst(1, ldihed + kdihed) + isite
            jatm = dihedral%lst(2, ldihed + kdihed) + isite
            katm = dihedral%lst(3, ldihed + kdihed) + isite
            latm = dihedral%lst(4, ldihed + kdihed) + isite
            If (dihedral%l_core_shell) Then
              matm = dihedral%lst(5, ldihed + kdihed) + isite
              natm = dihedral%lst(6, ldihed + kdihed) + isite
            End If

            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)
            kat0 = local_index(katm, config%nlast, config%lsi, config%lsa)
            lat0 = local_index(latm, config%nlast, config%lsi, config%lsa)
            If (dihedral%l_core_shell) Then
              mat0 = local_index(matm, config%nlast, config%lsi, config%lsa)
              nat0 = local_index(natm, config%nlast, config%lsi, config%lsa)
            Else
              mat0 = 0
              nat0 = 0
            End If

            If (iat0 > config%natms) iat0 = 0
            If (jat0 > config%natms) jat0 = 0
            If (kat0 > config%natms) kat0 = 0
            If (lat0 > config%natms) lat0 = 0
            If (dihedral%l_core_shell) Then
              If (mat0 > config%natms) mat0 = 0
              If (nat0 > config%natms) nat0 = 0
            Else
              mat0 = 0
              nat0 = 0
            End If

            If (iat0 > 0 .or. jat0 > 0 .or. kat0 > 0 .or. lat0 > 0 .or. & ! dihedral%l_core_shell=.false.
                mat0 > 0 .or. nat0 > 0) Then ! lx_dix=.true.
              jdihed = jdihed + 1
              If (jdihed <= dihedral%max_angles) Then
                dihedral%list(0, jdihed) = ldihed + kdihed
                dihedral%list(1, jdihed) = iatm
                dihedral%list(2, jdihed) = jatm
                dihedral%list(3, jdihed) = katm
                dihedral%list(4, jdihed) = latm
                If (dihedral%l_core_shell) Then
                  dihedral%list(5, jdihed) = matm
                  dihedral%list(6, jdihed) = natm
                End If

                If (iat0 > 0) Then
                  Call tag_legend(safe(1), iat0, jdihed, dihedral%legend, dihedral%max_legend)
                  If (dihedral%legend(dihedral%max_legend, iat0) > 0) Then
                    Call warning('too many dihedral type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      dihedral%max_legend - 1 + dihedral%legend(dihedral%max_legend, iat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', dihedral%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', dihedral%lst(1, ldihed + kdihed)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', ldihed
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (jat0 > 0) Then
                  Call tag_legend(safe(1), jat0, jdihed, dihedral%legend, dihedral%max_legend)
                  If (dihedral%legend(dihedral%max_legend, jat0) > 0) Then
                    Call warning('too many dihedral type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      dihedral%max_legend - 1 + dihedral%legend(dihedral%max_legend, jat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', dihedral%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', dihedral%lst(2, ldihed + kdihed)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', ldihed
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (kat0 > 0) Then
                  Call tag_legend(safe(1), kat0, jdihed, dihedral%legend, dihedral%max_legend)
                  If (dihedral%legend(dihedral%max_legend, kat0) > 0) Then
                    Call warning('too many dihedral type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      dihedral%max_legend - 1 + dihedral%legend(dihedral%max_legend, kat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', dihedral%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', dihedral%lst(3, ldihed + kdihed)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', ldihed
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (lat0 > 0) Then
                  Call tag_legend(safe(1), lat0, jdihed, dihedral%legend, dihedral%max_legend)
                  If (dihedral%legend(dihedral%max_legend, lat0) > 0) Then
                    Call warning('too many dihedral type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      dihedral%max_legend - 1 + dihedral%legend(dihedral%max_legend, lat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', dihedral%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', dihedral%lst(4, ldihed + kdihed)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', ldihed
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (mat0 > 0) Then
                  Call tag_legend(safe(1), mat0, jdihed, dihedral%legend, dihedral%max_legend)
                  If (dihedral%legend(dihedral%max_legend, mat0) > 0) Then
                    Call warning('too many dihedral type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      dihedral%max_legend - 1 + dihedral%legend(dihedral%max_legend, mat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', dihedral%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', dihedral%lst(4, ldihed + kdihed)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', ldihed
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (nat0 > 0) Then
                  Call tag_legend(safe(1), nat0, jdihed, dihedral%legend, dihedral%max_legend)
                  If (dihedral%legend(dihedral%max_legend, nat0) > 0) Then
                    Call warning('too many dihedral type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      dihedral%max_legend - 1 + dihedral%legend(dihedral%max_legend, nat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', dihedral%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', dihedral%lst(4, ldihed + kdihed)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', ldihed
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If
              Else
                idihed = idihed + 1
                safe(9) = .false.
              End If
            End If
          End Do

          ! Construct inversion potential interaction list

          Do linver = 1, inversion%num(itmols)
            iatm = inversion%lst(1, linver + kinver) + isite
            jatm = inversion%lst(2, linver + kinver) + isite
            katm = inversion%lst(3, linver + kinver) + isite
            latm = inversion%lst(4, linver + kinver) + isite

            iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
            jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)
            kat0 = local_index(katm, config%nlast, config%lsi, config%lsa)
            lat0 = local_index(latm, config%nlast, config%lsi, config%lsa)

            If (iat0 > config%natms) iat0 = 0
            If (jat0 > config%natms) jat0 = 0
            If (kat0 > config%natms) kat0 = 0
            If (lat0 > config%natms) lat0 = 0

            If (iat0 > 0 .or. jat0 > 0 .or. kat0 > 0 .or. lat0 > 0) Then
              jinver = jinver + 1
              If (jinver <= inversion%max_angles) Then
                inversion%list(0, jinver) = linver + kinver
                inversion%list(1, jinver) = iatm
                inversion%list(2, jinver) = jatm
                inversion%list(3, jinver) = katm
                inversion%list(4, jinver) = latm

                If (iat0 > 0) Then
                  Call tag_legend(safe(1), iat0, jinver, inversion%legend, inversion%max_legend)
                  If (inversion%legend(inversion%max_legend, iat0) > 0) Then
                    Call warning('too many inversion type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      inversion%max_legend - 1 + inversion%legend(inversion%max_legend, iat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', inversion%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', inversion%lst(1, linver + kinver)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', linver
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (jat0 > 0) Then
                  Call tag_legend(safe(1), jat0, jinver, inversion%legend, inversion%max_legend)
                  If (inversion%legend(inversion%max_legend, jat0) > 0) Then
                    Call warning('too many inversion type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      inversion%max_legend - 1 + inversion%legend(inversion%max_legend, jat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', inversion%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', inversion%lst(2, linver + kinver)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', linver
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (kat0 > 0) Then
                  Call tag_legend(safe(1), kat0, jinver, inversion%legend, inversion%max_legend)
                  If (inversion%legend(inversion%max_legend, kat0) > 0) Then
                    Call warning('too many inversion type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      inversion%max_legend - 1 + inversion%legend(inversion%max_legend, kat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', inversion%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', inversion%lst(3, linver + kinver)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', linver
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If

                If (lat0 > 0) Then
                  Call tag_legend(safe(1), lat0, jinver, inversion%legend, inversion%max_legend)
                  If (inversion%legend(inversion%max_legend, lat0) > 0) Then
                    Call warning('too many inversion type neighbours')
                    Write (messages(1), '(a,i0)') 'requiring a list length of: ', &
                      inversion%max_legend - 1 + inversion%legend(inversion%max_legend, lat0)
                    Write (messages(2), '(a,i0)') 'but maximum length allowed: ', inversion%max_legend - 1
                    Write (messages(3), '(a,i0)') 'for particle (global ID #): ', iatm
                    Write (messages(4), '(a,i0)') 'on mol. site (local  ID #): ', inversion%lst(4, linver + kinver)
                    Write (messages(5), '(a,i0)') 'of unit      (local  ID #): ', linver
                    Write (messages(6), '(a,i0)') 'in molecule  (local  ID #): ', imols
                    Write (messages(7), '(a,i0)') 'of type      (       ID #): ', itmols
                    Call info(messages, 7)
                  End If
                End If
              Else
                iinver = iinver + 1
                safe(10) = .false.
              End If
            End If
          End Do

        End If

        isite = isite + sites%num_site(itmols)
        nsatm = neatm

      End Do

      ! Update global unit numbers for all passed molecules so far

      kshels = kshels + cshell%numshl(itmols)

      kconst = kconst + cons%numcon(itmols)
      ! No 'kpmf' needed since PMF is defined on one and only one molecular type

      krigid = krigid + rigid%num(itmols)

      kteths = kteths + tether%numteth(itmols)

      kbonds = kbonds + bond%num(itmols)
      kangle = kangle + angle%num(itmols)
      kdihed = kdihed + dihedral%num(itmols)
      kinver = kinver + inversion%num(itmols)

    End Do

    ! Store array counters for bookkeeping

    cshell%ntshl = jshels

    cons%ntcons = jconst
    ! 'pmf%ntpmf' is updated locally as PMFs are global and one type only

    rigid%n_types = jrigid

    tether%ntteth = jteths

    bond%n_types = jbonds
    angle%n_types = jangle
    dihedral%n_types = jdihed
    inversion%n_types = jinver

    If (cshell%megshl == 0) Then
      cshell%ntshl1 = cshell%ntshl

      cons%ntcons1 = cons%ntcons

      rigid%n_types_book = rigid%n_types

      bond%n_types1 = bond%n_types
      angle%n_types1 = angle%n_types
      dihedral%n_types1 = dihedral%n_types
      inversion%n_types1 = inversion%n_types

      cshell%ntshl2 = cshell%ntshl1

      Go To 400
    End If

    ! Cycle through all constraint, RB, bond, angle, dihedral and inversion
    ! units on this node and record the non-local index particles

    iwrk = 0
    mshels = 0
    Do i = 1, cons%ntcons
      iatm = cons%listcon(1, i)
      jatm = cons%listcon(2, i)

      iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
      jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)

      If (iat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == iatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (jat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == jatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If
    End Do
    Do i = 1, rigid%n_types
      mrigid = rigid%list(-1, i)

      Do lrigid = 1, mrigid
        iatm = rigid%list(lrigid, i)
        iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)

        If (iat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == iatm))) Then
          mshels = mshels + 1
          iwrk(mshels) = iatm
        End If
      End Do
    End Do
    Do i = 1, bond%n_types
      iatm = bond%list(1, i)
      jatm = bond%list(2, i)

      iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
      jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)

      If (iat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == iatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (jat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == jatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If
    End Do
    Do i = 1, angle%n_types
      iatm = angle%list(1, i)
      jatm = angle%list(2, i)
      katm = angle%list(3, i)

      iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
      jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)
      kat0 = local_index(katm, config%nlast, config%lsi, config%lsa)

      If (iat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == iatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (jat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == jatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If

      If (kat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == katm))) Then
        mshels = mshels + 1
        iwrk(mshels) = katm
      End If
    End Do
    Do i = 1, dihedral%n_types
      iatm = dihedral%list(1, i)
      jatm = dihedral%list(2, i)
      katm = dihedral%list(3, i)
      latm = dihedral%list(4, i)
      If (dihedral%l_core_shell) Then
        matm = dihedral%list(5, i)
        natm = dihedral%list(6, i)
      End If

      iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
      jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)
      kat0 = local_index(katm, config%nlast, config%lsi, config%lsa)
      lat0 = local_index(latm, config%nlast, config%lsi, config%lsa)
      If (dihedral%l_core_shell) Then
        mat0 = local_index(matm, config%nlast, config%lsi, config%lsa)
        nat0 = local_index(natm, config%nlast, config%lsi, config%lsa)
      Else
        mat0 = 0
        nat0 = 0
      End If

      If (iat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == iatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (jat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == jatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If

      If (kat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == katm))) Then
        mshels = mshels + 1
        iwrk(mshels) = katm
      End If

      If (lat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == latm))) Then
        mshels = mshels + 1
        iwrk(mshels) = latm
      End If

      If (mat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == matm))) Then
        mshels = mshels + 1
        iwrk(mshels) = matm
      End If

      If (nat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == natm))) Then
        mshels = mshels + 1
        iwrk(mshels) = natm
      End If
    End Do
    Do i = 1, inversion%n_types
      iatm = inversion%list(1, i)
      jatm = inversion%list(2, i)
      katm = inversion%list(3, i)
      latm = inversion%list(4, i)

      iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa)
      jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa)
      kat0 = local_index(katm, config%nlast, config%lsi, config%lsa)
      lat0 = local_index(latm, config%nlast, config%lsi, config%lsa)

      If (iat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == iatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (jat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == jatm))) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If

      If (kat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == katm))) Then
        mshels = mshels + 1
        iwrk(mshels) = katm
      End If
      If (lat0 > config%natms .and. (.not. Any(iwrk(1:mshels) == latm))) Then
        mshels = mshels + 1
        iwrk(mshels) = latm
      End If
    End Do

    If (mshels == 0) Go To 100

    ! Include in local core-shell units' description
    ! non-local units that are connected to partly shared
    ! constraints, RBs, bonds, angles, dihedrals and inversions

    isite = 0 ! initialise bookkeeping indices
    kshels = 0 ! last index of core-shell unit
    nsatm = 0 ! global atom counter
    Do itmols = 1, sites%ntype_mol ! loop over molecule types in the system
      Do imols = 1, sites%num_mols(itmols) ! loop over molecules of this type
        neatm = nsatm + sites%num_site(itmols) ! last atom in the molecule

        ! From the first till the last atom of this molecule, is there a
        ! non-local atom iwrk(1:mshels)

        i = 0 ! There is not
        Do iatm = nsatm + 1, neatm
          If (Any(iwrk(1:mshels) == iatm)) i = i + 1
        End Do

        ! If there is a non-local atom iwrk(1:mshels) on this
        ! molecule on this node, extend cshell%listshl

        If (i > 0) Then

          ! Extend core-shell units interaction list

          Do lshels = 1, cshell%numshl(itmols)
            iatm = cshell%lstshl(1, lshels + kshels) + isite
            jatm = cshell%lstshl(2, lshels + kshels) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm)) Then
              jshels = jshels + 1
              If (jshels <= cshell%mxshl) Then
                cshell%listshl(0, jshels) = lshels + kshels
                cshell%listshl(1, jshels) = iatm
                cshell%listshl(2, jshels) = jatm
              Else
                safe(2) = .false.
              End If
            End If
          End Do

        End If

        isite = isite + sites%num_site(itmols)
        nsatm = neatm
      End Do

      ! Update core-shell units number for all passed molecules so far

      kshels = kshels + cshell%numshl(itmols)
    End Do

    100 Continue

    ! Store first (local+non-local) extended array counter
    ! for bookkeeping and exclusion of core-shell units

    cshell%ntshl1 = jshels

    ! Cycle through all local and partly shared core-shell units on
    ! this node and record the non-local indices of cross-domained
    ! core-shell unit particles

    iwrk = 0
    mshels = 0
    Do i = 1, cshell%ntshl
      iatm = cshell%listshl(1, i)
      jatm = cshell%listshl(2, i)

      iat0 = local_index(iatm, config%nlast, config%lsi, config%lsa) ! This is a core
      jat0 = local_index(jatm, config%nlast, config%lsi, config%lsa) ! This is a shell

      If (iat0 == 0 .and. jat0 == 0) safe(11) = .false.

      ! Core is out

      If (iat0 > config%natms) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      ! Shell is out

      If (jat0 > config%natms) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If
    End Do

    If (mshels == 0) Go To 200

    ! Include in (local constraint, RB, bond, angle, dihedral
    ! and inversion) units' descriptions non-local units that
    ! are connected to partly shared core-shell units

    isite = 0 ! initialise bookkeeping indices
    kconst = 0 ! last index of constraint unit
    krigid = 0 ! last index of RB unit
    kbonds = 0 ! last index of bond unit
    kangle = 0 ! last index of angle unit
    kdihed = 0 ! last index of dihedral unit
    kinver = 0 ! last index of inversion unit
    nsatm = 0 ! global atom counter
    Do itmols = 1, sites%ntype_mol ! loop over molecule types in the system
      Do imols = 1, sites%num_mols(itmols) ! loop over molecules of this type
        neatm = nsatm + sites%num_site(itmols) ! last atom in the molecule

        ! From the first till the last atom of this molecule, is there a
        ! non-local, cross-domained core-shell unit atom

        i = 0 ! There is not
        Do iatm = nsatm + 1, neatm
          If (Any(iwrk(1:mshels) == iatm)) i = i + 1
        End Do

        ! If there is a non-local, cross-domained core-shell unit atom on this
        ! molecule on this node, extend cons%listcon, rigid%list, bond%list and angle%list

        If (i > 0) Then

          ! Extend constraint bond list

          Do lconst = 1, cons%numcon(itmols)
            iatm = cons%lstcon(1, lconst + kconst) + isite
            jatm = cons%lstcon(2, lconst + kconst) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm)) Then
              jconst = jconst + 1
              If (jconst <= cons%mxcons) Then
                cons%listcon(0, jconst) = lconst + kconst
                cons%listcon(1, jconst) = iatm
                cons%listcon(2, jconst) = jatm
              Else
                safe(3) = .false.
              End If
            End If
          End Do

          ! Extend RB list

          Do lrigid = 1, rigid%num(itmols)
            mrigid = rigid%lst(0, lrigid + krigid)

            irgd = 0; irgd0 = 0; go = .false.
            Do irigid = 1, mrigid
              irgd(irigid) = rigid%lst(irigid, lrigid + krigid) + isite
              go = (go .or. Any(iwrk(1:mshels) == irgd(irigid)))
            End Do

            If (go) Then
              jrigid = jrigid + 1
              If (jrigid <= rigid%max_rigid) Then
                rigid%list(-1, jrigid) = rigid%lst(0, lrigid + krigid)
                rigid%list(0, jrigid) = lrigid + krigid
                Do irigid = 1, mrigid
                  rigid%list(irigid, jrigid) = irgd(irigid)
                End Do
              Else
                safe(5) = .false.
              End If
            End If
          End Do

          ! Extend chemical bond interaction list

          Do lbonds = 1, bond%num(itmols)
            iatm = bond%lst(1, lbonds + kbonds) + isite
            jatm = bond%lst(2, lbonds + kbonds) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm)) Then
              jbonds = jbonds + 1
              If (jbonds <= bond%max_bonds) Then
                bond%list(0, jbonds) = lbonds + kbonds
                bond%list(1, jbonds) = iatm
                bond%list(2, jbonds) = jatm
              Else
                safe(7) = .false.
              End If
            End If
          End Do

          ! Extend valence angle interaction list

          Do langle = 1, angle%num(itmols)
            iatm = angle%lst(1, langle + kangle) + isite
            jatm = angle%lst(2, langle + kangle) + isite
            katm = angle%lst(3, langle + kangle) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm) .or. &
                Any(iwrk(1:mshels) == katm)) Then
              jangle = jangle + 1
              If (jangle <= angle%max_angles) Then
                angle%list(0, jangle) = langle + kangle
                angle%list(1, jangle) = iatm
                angle%list(2, jangle) = jatm
                angle%list(3, jangle) = katm
              Else
                safe(8) = .false.
              End If
            End If
          End Do

          ! Extend dihedral angle interaction list

          Do ldihed = 1, dihedral%num(itmols)
            iatm = dihedral%lst(1, ldihed + kdihed) + isite
            jatm = dihedral%lst(2, ldihed + kdihed) + isite
            katm = dihedral%lst(3, ldihed + kdihed) + isite
            latm = dihedral%lst(4, ldihed + kdihed) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm) .or. &
                Any(iwrk(1:mshels) == katm) .or. &
                Any(iwrk(1:mshels) == latm)) Then

              If (dihedral%l_core_shell) Then
                matm = dihedral%lst(5, ldihed + kdihed) + isite
                natm = dihedral%lst(6, ldihed + kdihed) + isite
                If (.not. (Any(iwrk(1:mshels) == jatm) .or. &
                           Any(iwrk(1:mshels) == katm))) Exit
              End If

              jdihed = jdihed + 1
              If (jdihed <= dihedral%max_angles) Then
                dihedral%list(0, jdihed) = ldihed + kdihed
                dihedral%list(1, jdihed) = iatm
                dihedral%list(2, jdihed) = jatm
                dihedral%list(3, jdihed) = katm
                dihedral%list(4, jdihed) = latm
                If (dihedral%l_core_shell) Then
                  dihedral%list(5, jdihed) = matm
                  dihedral%list(6, jdihed) = natm
                End If
              Else
                safe(9) = .false.
              End If
            End If
          End Do

          ! Extend inversion potential interaction list

          Do linver = 1, inversion%num(itmols)
            iatm = inversion%lst(1, linver + kinver) + isite
            jatm = inversion%lst(2, linver + kinver) + isite
            katm = inversion%lst(3, linver + kinver) + isite
            latm = inversion%lst(4, linver + kinver) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm) .or. &
                Any(iwrk(1:mshels) == katm) .or. &
                Any(iwrk(1:mshels) == latm)) Then
              jinver = jinver + 1
              If (jinver <= inversion%max_angles) Then
                inversion%list(0, jinver) = linver + kinver
                inversion%list(1, jinver) = iatm
                inversion%list(2, jinver) = jatm
                inversion%list(3, jinver) = katm
                inversion%list(4, jinver) = latm
              Else
                safe(10) = .false.
              End If
            End If
          End Do

        End If

        isite = isite + sites%num_site(itmols)
        nsatm = neatm
      End Do

      ! Update constraint, RB, bond, angle, dihedral and inversion
      ! units numbers for all passed molecules so far

      kconst = kconst + cons%numcon(itmols)

      krigid = krigid + rigid%num(itmols)

      kbonds = kbonds + bond%num(itmols)
      kangle = kangle + angle%num(itmols)
      kdihed = kdihed + dihedral%num(itmols)
      kinver = kinver + inversion%num(itmols)
    End Do

    200 Continue

    ! Store first extended array counters for bookkeeping and exclusion
    ! of constraint, RB, bond, angle, dihedral and inversion units

    cons%ntcons1 = jconst

    rigid%n_types_book = jrigid

    bond%n_types1 = jbonds
    angle%n_types1 = jangle
    dihedral%n_types1 = jdihed
    inversion%n_types1 = jinver

    ! Cycle through the extended -
    ! constraint, RB, bond, angle, dihedral and inversion units
    ! on this node and record the non-local index particles

    iwrk = 0
    mshels = 0
    Do i = cons%ntcons + 1, cons%ntcons1
      iatm = cons%listcon(1, i)
      jatm = cons%listcon(2, i)

      If (.not. Any(iwrk(1:mshels) == iatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (.not. Any(iwrk(1:mshels) == jatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If
    End Do
    Do i = rigid%n_types + 1, rigid%n_types_book
      mrigid = rigid%list(-1, i)

      Do j = 1, mrigid
        iatm = rigid%list(j, i)
        If (.not. Any(iwrk(1:mshels) == iatm)) Then
          mshels = mshels + 1
          iwrk(mshels) = iatm
        End If
      End Do
    End Do
    Do i = bond%n_types + 1, bond%n_types1
      iatm = bond%list(1, i)
      jatm = bond%list(2, i)

      If (.not. Any(iwrk(1:mshels) == iatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (.not. Any(iwrk(1:mshels) == jatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If
    End Do
    Do i = angle%n_types + 1, angle%n_types1
      iatm = angle%list(1, i)
      jatm = angle%list(2, i)
      katm = angle%list(3, i)

      If (.not. Any(iwrk(1:mshels) == iatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (.not. Any(iwrk(1:mshels) == jatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If

      If (.not. Any(iwrk(1:mshels) == katm)) Then
        mshels = mshels + 1
        iwrk(mshels) = katm
      End If
    End Do
    Do i = dihedral%n_types + 1, dihedral%n_types1
      iatm = dihedral%list(1, i)
      jatm = dihedral%list(2, i)
      katm = dihedral%list(3, i)
      latm = dihedral%list(4, i)
      If (dihedral%l_core_shell) Then
        matm = dihedral%list(5, i)
        natm = dihedral%list(6, i)
      End If

      If (.not. Any(iwrk(1:mshels) == iatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (.not. Any(iwrk(1:mshels) == jatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If

      If (.not. Any(iwrk(1:mshels) == katm)) Then
        mshels = mshels + 1
        iwrk(mshels) = katm
      End If

      If (.not. Any(iwrk(1:mshels) == latm)) Then
        mshels = mshels + 1
        iwrk(mshels) = latm
      End If

      If (dihedral%l_core_shell) Then
        If (.not. Any(iwrk(1:mshels) == matm)) Then
          mshels = mshels + 1
          iwrk(mshels) = matm
        End If

        If (.not. Any(iwrk(1:mshels) == natm)) Then
          mshels = mshels + 1
          iwrk(mshels) = natm
        End If
      End If
    End Do
    Do i = inversion%n_types + 1, inversion%n_types1
      iatm = inversion%list(1, i)
      jatm = inversion%list(2, i)
      katm = inversion%list(3, i)
      latm = inversion%list(4, i)

      If (.not. Any(iwrk(1:mshels) == iatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = iatm
      End If

      If (.not. Any(iwrk(1:mshels) == jatm)) Then
        mshels = mshels + 1
        iwrk(mshels) = jatm
      End If

      If (.not. Any(iwrk(1:mshels) == katm)) Then
        mshels = mshels + 1
        iwrk(mshels) = katm
      End If

      If (.not. Any(iwrk(1:mshels) == latm)) Then
        mshels = mshels + 1
        iwrk(mshels) = latm
      End If
    End Do

    If (mshels == 0) Go To 300

    ! Include in local and non-local - connected to partly shared
    ! constraints, RBs, bonds, angles, dihedrals or inversions -
    ! core-shell unit description of foreign units that are connected to
    ! non-local constraints, RBs, bonds, angles, dihedrals or inversions
    ! by partly shared core-shell units

    isite = 0 ! initialise bookkeeping indices
    kshels = 0 ! last index of core-shell unit
    nsatm = 0 ! global atom counter
    Do itmols = 1, sites%ntype_mol ! loop over molecule types in the system
      Do imols = 1, sites%num_mols(itmols) ! loop over molecules of this type
        neatm = nsatm + sites%num_site(itmols) ! last atom in the molecule

        ! From the first till the last atom of this molecule, is there a
        ! non-local atom iwrk(1:mshels)

        i = 0 ! There is not
        Do iatm = nsatm + 1, neatm
          If (Any(iwrk(1:mshels) == iatm)) i = i + 1
        End Do

        ! If there is a non-local atom iwrk(1:mshels) on this
        ! molecule on this node, extend cshell%listshl

        If (i > 0) Then

          ! Extend core-shell units interaction list

          Do lshels = 1, cshell%numshl(itmols)
            iatm = cshell%lstshl(1, lshels + kshels) + isite
            jatm = cshell%lstshl(2, lshels + kshels) + isite

            If (Any(iwrk(1:mshels) == iatm) .or. &
                Any(iwrk(1:mshels) == jatm)) Then
              jshels = jshels + 1
              If (jshels <= cshell%mxshl) Then
                cshell%listshl(0, jshels) = lshels + kshels
                cshell%listshl(1, jshels) = iatm
                cshell%listshl(2, jshels) = jatm
              Else
                safe(2) = .false.
              End If
            End If
          End Do

        End If

        isite = isite + sites%num_site(itmols)
        nsatm = neatm
      End Do

      ! Update core-shell units number for all passed molecules so far

      kshels = kshels + cshell%numshl(itmols)
    End Do

    300 Continue

    ! Store second (local+non-local+foreign) extended array counter
    ! for bookkeeping and exclusion of core-shell units

    cshell%ntshl2 = jshels

    !  If (dihedral%l_core_shell) dihedral%n_types=dihedral%n_types1 ! extend the dihedrals' set

    400 Continue

    ! error exit for all error conditions (size of work arrays)

    Call gcheck(comm, safe)
    If (.not. safe(1)) Call error(88)

    If (Any(.not. safe)) Then
      itmp(1) = ishels; jtmp(1) = cshell%mxshl
      itmp(2) = iconst; jtmp(2) = cons%mxcons
      itmp(3) = ipmf; jtmp(3) = pmf%mxpmf
      itmp(4) = irigid; jtmp(4) = rigid%max_rigid
      itmp(5) = iteths; jtmp(5) = tether%mxteth
      itmp(6) = ibonds; jtmp(6) = bond%max_bonds
      itmp(7) = iangle; jtmp(7) = angle%max_angles
      itmp(8) = idihed; jtmp(8) = dihedral%max_angles
      itmp(9) = iinver; jtmp(9) = inversion%max_angles

      Call gmax(comm, itmp(1:9))

      tmp = 1.0_wp
      Do i = 1, 9
        tmp = Max(tmp, 1.0_wp + Real(itmp(i), wp) / Real(Max(1, jtmp(i)), wp))
      End Do

      Write (message, '(a,i0)') 'estimated densvar value for passing this stage safely is : ', &
        Nint((config%dvar * tmp - 1.0_wp) * 100.0_wp + 0.5_wp)
      Call warning(message, .true.)
    End If

    If (.not. safe(2)) Call error(59)
    If (.not. safe(3)) Call error(41)
    If (.not. safe(4)) Call error(488)
    If (.not. safe(5)) Call error(640)
    If (.not. safe(6)) Call error(63)
    If (.not. safe(7)) Call error(31)
    If (.not. safe(8)) Call error(51)
    If (.not. safe(9)) Call error(61)
    If (.not. safe(10)) Call error(77)
    If (.not. safe(11)) Call error(64)

    Deallocate (iwrk, Stat=fail(1))
    If (rigid%on) Then
      Deallocate (irgd, irgd0, Stat=fail(2))
    End If
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'build_book_intra deallocation failure'
      Call error(0, message)
    End If

    If (flow%newjob_build_book) Then

      ! Set RB particulars and quaternions

      If (rigid%on) Then
        Call rigid_bodies_setup(l_str, l_top, config%megatm, config%megfrz, config%degtra, config%degrot, &
                                neigh%cutoff, sites, rigid, config, comm)
      End If

      Call report_topology(config%megatm, config%megfrz, config%atmfre, config%atmfrz, &
                           cshell, cons, pmf, bond, angle, dihedral, inversion, tether, sites, rigid)

      ! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed

      If (lsim) Then
        Call cshell%deallocate_core_shell_tmp_arrays()

        Call cons%deallocate_constraints_temps()
        Call pmf%deallocate_pmf_tmp_arrays()

        Call rigid%deallocate_temp()

        Call tether%deallocate_temp()
      End If

    Else

      ! Tag RBs, find their COMs and check their widths to neigh%cutoff (system cutoff)

      Call rigid_bodies_tags(config, rigid, comm)
      Call rigid_bodies_coms(config, rigid%xxx, rigid%yyy, rigid%zzz, rigid)
      Call rigid_bodies_widths(neigh%cutoff, rigid, config, comm)

    End If

    ! Update shared core-shell, constraint and RB units
    ! (pmf data updated by construction)

    If (cshell%megshl > 0 .and. comm%mxnode > 1) Call pass_shared_units &
      (config, cshell%mxshl, Lbound(cshell%listshl, Dim=1), Ubound(cshell%listshl, Dim=1), cshell%ntshl, &
       cshell%listshl, cshell%mxfshl, cshell%legshl, cshell%lshmv_shl, cshell%lishp_shl, cshell%lashp_shl, &
       flow%oldjob_shared_units, domain, comm, &
       rigid%q0, rigid%q1, rigid%q2, rigid%q3, rigid%vxx, rigid%vyy, rigid%vzz, &
       rigid%oxx, rigid%oyy, rigid%ozz)

    If (cons%m_con > 0 .and. comm%mxnode > 1) Call pass_shared_units &
      (config, cons%mxcons, Lbound(cons%listcon, Dim=1), Ubound(cons%listcon, Dim=1), cons%ntcons, cons%listcon, cons%mxfcon, &
       cons%legcon, cons%lshmv_con, &
       cons%lishp_con, cons%lashp_con, flow%oldjob_shared_units, domain, comm, &
       rigid%q0, rigid%q1, rigid%q2, rigid%q3, rigid%vxx, rigid%vyy, rigid%vzz, &
       rigid%oxx, rigid%oyy, rigid%ozz)

    If (rigid%on .and. comm%mxnode > 1) Call pass_shared_units &
      (config, rigid%max_rigid, Lbound(rigid%list, Dim=1), Ubound(rigid%list, Dim=1), rigid%n_types, &
       rigid%list, rigid%max_frozen, rigid%legend, rigid%share, rigid%list_shared, &
       rigid%map_shared, flow%oldjob_shared_units, domain, comm, &
       rigid%q0, rigid%q1, rigid%q2, rigid%q3, rigid%vxx, rigid%vyy, rigid%vzz, &
       rigid%oxx, rigid%oyy, rigid%ozz)

  End Subroutine build_book_intra

  Subroutine compress_book_intra(mx_u, nt_u, b_u, list_u, leg_u, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to prevent bookkeeping arrays from expanding
    ! when execution is on many nodes, mxnode>1, (shells, constraints, PMFs
    ! and RBs are dealt differently pass_shared_units and pmf_units_set)
    !
    ! Note: This routine is to be only called from relocate_particles
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov october 2012
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: mx_u
    Integer,                  Intent(InOut) :: nt_u
    Integer,                  Intent(In   ) :: b_u
    Integer,                  Intent(InOut) :: list_u(0:b_u, 1:mx_u), leg_u(0:, 1:)
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Integer :: i, j, k, l, m
    Logical :: keep_k, keep_nt, ok

    ! is it ok not to do it since it's safe - there's enough buffering space

    ok = .true.
    If (mx_u > 0) Then
      ok = .not. (Real(nt_u, wp) / Real(mx_u, wp) > 0.85_wp)
      Call gcheck(comm, ok)
    End If

    If (.not. ok) Then
      k = 0
      Do While (k < nt_u)
        k = k + 1
        10 Continue

        keep_k = .false.
        Do i = 1, b_u
          keep_k = keep_k .or. (local_index(list_u(i, k), config%natms, config%lsi, config%lsa) /= 0)
        End Do

        If (.not. keep_k) Then
          20 Continue

          ! If the whole unit has moved out of this node - compress list_u and leg_u

          If (k < nt_u) Then

            keep_nt = .false.
            Do i = 1, b_u
              j = local_index(list_u(i, nt_u), config%natms, config%lsi, config%lsa)
              If (j > 0) Then ! For all particles in list_u(1:b_u,nt_u),
                keep_nt = .true. ! [indicate that this unit is being kept]
                m = leg_u(0, j) ! if present on this node, repoint unit
                Do l = 1, m ! 'nt_u' to 'k' in their leg_u array
                  If (leg_u(l, j) == nt_u) leg_u(l, j) = k
                End Do
              End If
            End Do

            If (keep_nt) Then ! Do repointing
              list_u(:, k) = list_u(:, nt_u) ! Copy list content from 'nt_u' to 'k'
              list_u(:, nt_u) = 0 ! Remove list content in 'nt_u'
              nt_u = nt_u - 1 ! Reduce 'nt_u' pointer
            Else
              list_u(:, nt_u) = 0 ! Remove list content in 'nt_u'
              nt_u = nt_u - 1 ! Reduce 'nt_u' pointer

              Go To 20 ! Go back and check again for the new list contents in 'nt_u'
            End If

            Go To 10 ! Go back and check it all again for the new list contents in 'k'

          Else If (k == nt_u) Then

            list_u(:, nt_u) = 0 ! Remove list content in 'k=nt_u'
            nt_u = nt_u - 1 ! Reduce 'nt_u' pointer

          End If
        End If
      End Do
    End If

  End Subroutine compress_book_intra

  Subroutine init_intra(cshell, cons, pmf, bond, angle, dihedral, inversion, tether, neigh, rigid)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for initialising local bookkeeping of intra-like
    ! interactions: core-shell, bond constraints, PMF constraints, RBs,
    ! tethered atoms, chemical bonds, valence angles, torsion and improper
    ! torsion angles, and inversion angles; with exclusions too at the top
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(core_shell_type),   Intent(InOut) :: cshell
    Type(constraints_type),  Intent(InOut) :: cons
    Type(pmf_type),          Intent(InOut) :: pmf
    Type(bonds_type),        Intent(InOut) :: bond
    Type(angles_type),       Intent(InOut) :: angle
    Type(dihedrals_type),    Intent(InOut) :: dihedral
    Type(inversions_type),   Intent(InOut) :: inversion
    Type(tethers_type),      Intent(InOut) :: tether
    Type(neighbours_type),   Intent(InOut) :: neigh
    Type(rigid_bodies_type), Intent(InOut) :: rigid

    ! exclusions locals

    neigh%list_excl = 0

    ! core-shell locals

    cshell%ntshl = 0; cshell%ntshl1 = 0; cshell%ntshl2 = 0
    cshell%listshl = 0
    cshell%legshl = 0

    ! constraints locals

    cons%ntcons = 0; cons%ntcons1 = 0
    cons%listcon = 0
    cons%legcon = 0

    ! PMFs locals

    pmf%ntpmf = 0
    pmf%listpmf = 0
    pmf%legpmf = 0

    ! RBs locals

    rigid%n_types = 0; rigid%n_types_book = 0
    rigid%list = 0
    rigid%legend = 0

    rigid%frozen = 0; rigid%index_global = 0; rigid%index_local = 0

    rigid%weight = 0.0_wp; rigid%weightless = 0.0_wp
    rigid%x = 0.0_wp; rigid%y = 0.0_wp; rigid%z = 0.0_wp
    rigid%rix = 0.0_wp; rigid%riy = 0.0_wp; rigid%riz = 0.0_wp
    rigid%axs = 0.0_wp

    rigid%q0 = 0.0_wp; rigid%q1 = 0.0_wp; rigid%q2 = 0.0_wp; rigid%q3 = 0.0_wp

    rigid%xxx = 0.0_wp; rigid%yyy = 0.0_wp; rigid%zzz = 0.0_wp
    rigid%vxx = 0.0_wp; rigid%vyy = 0.0_wp; rigid%vzz = 0.0_wp
    rigid%oxx = 0.0_wp; rigid%oyy = 0.0_wp; rigid%ozz = 0.0_wp

    ! tethers locals

    tether%ntteth = 0
    tether%listtet = 0
    tether%legtet = 0

    ! bonds locals

    bond%n_types = 0; bond%n_types1 = 0
    bond%list = 0
    bond%legend = 0

    ! angles locals

    angle%n_types = 0; angle%n_types1 = 0
    angle%list = 0
    angle%legend = 0

    ! dihedrals locals

    dihedral%n_types = 0; dihedral%n_types1 = 0
    dihedral%list = 0
    dihedral%legend = 0

    ! inversions locals

    inversion%n_types = 0; inversion%n_types1 = 0
    inversion%list = 0
    inversion%legend = 0

  End Subroutine init_intra
End Module build_book
