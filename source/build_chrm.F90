Module build_chrm
  ! SETUP MODULES

  Use kinds, only : wp,wi
  Use comms,  Only : comms_type,gcheck,gmax,gsum

  Use configuration, Only : configuration_type

  Use core_shell, Only : core_shell_type

  Use rigid_bodies, Only : rigid_bodies_type

  Use bonds, Only : bonds_type
  Use angles, Only : angles_type
  Use dihedrals, Only : dihedrals_type
  Use inversions, Only : inversions_type

  Use mpole, Only : mpole_type,POLARISATION_DEFAULT
  Use numerics, Only : local_index,shellsort
  Use build_excl, Only : add_exclusion
  Use constraints, Only : constraints_type
  Use errors_warnings, Only : error,warning,info
  Implicit None

  Private

  Public :: build_chrm_intra
Contains
  Subroutine build_chrm_intra(max_exclude,cshell,cons,bond,angle,dihedral, &
      inversion,mpoles,rigid,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for constructing the global CHARMMing core-shell
    ! units cross pair interactions list of the part of the MD system mapped
    ! onto this node.  While building the list it is presumed that:
    ! (1) constraint bonds are on top of chemical bonds
    ! (2) all 1-2 and 1-2-3 components of the conventional bonded
    !     interactions are checked out, i.e. 1-2-3 + 2-3-4 in a dihedral are
    !     scanned, and decomposed by primitive pairs, while 1-4 is excluded!
    !     all pairs are accounted for an inversion!
    ! (3) if a bonded interactions is RBed then any possible core-core pairs
    !     are removed and treated as non-contributing
    ! (4) if a bonded interaction is frozen then any possible core-core
    !     pairs are removed and treated as non-contributing
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

    Integer( Kind = wi ), Intent( In    ) :: max_exclude
    Type( constraints_type ), Intent( In    ) :: cons
    Type( bonds_type ), Intent( In    ) :: bond
    Type( angles_type ), Intent( In    ) :: angle
    Type( dihedrals_type ), Intent( In    ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( rigid_bodies_type ), Intent( In    ) :: rigid
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ),  Intent( InOut) :: comm

    Logical :: safe
    Integer :: fail,i,j,k,l,kk,                  &
      ia,ib,ia0,ib0,      &
      ibig,ja,jb,jc,jd,ja0,jb0,jc0,jd0, &
      ka,kb,kc,kd,ka0,kb0,kc0,kd0

    Character( Len = 256 ) :: message

    ! variables for array bound checking

    ibig=0
    safe=.true.

    ! go over the extended list of core-shell units

    Do i=1,cshell%ntshl2
      ia=cshell%listshl(1,i) ! This is the core
      ib=cshell%listshl(2,i) ! This is the shell

      ia0=local_index(ia,config%nlast,config%lsi,config%lsa)
      ib0=local_index(ib,config%nlast,config%lsi,config%lsa)

      If (ia0 > config%natms) ia0=0
      If (ib0 > config%natms) ib0=0

      ! add sites on basis of constraint bonds to core-shell units

      Do kk=1,cons%ntcons1
        ja=cons%listcon(1,kk)
        jb=cons%listcon(2,kk)

        ja0=local_index(ja,config%nlast,config%lsi,config%lsa)
        jb0=local_index(jb,config%nlast,config%lsi,config%lsa)

        If (ja0 > config%natms) ja0=0
        If (jb0 > config%natms) jb0=0

        ! Check for adjacent core or shell so that all possible, even between
        ! adjacent particles in any intra-unit, interactions are excluded.
        ! The checks commented out are not needed due to preventative checks
        ! in read_field by which shells are not allowed to be frozen,
        ! constrained bonded, RBed or tethered.

        ka = 0 ; ka0 = 0
        kb = 0 ; kb0 = 0

        Do j=1,cshell%ntshl2
          If (cshell%listshl(1,j) == ja) Then
            ka=cshell%listshl(2,j)
            ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
            If (ka0 > config%natms) ka0=0
            !           Else If (cshell%listshl(2,j) == ja) Then
            !              ka=cshell%listshl(1,j)
            !              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
            !              If (ka0 > config%natms) ka0=0
          End If

          If (cshell%listshl(1,j) == jb) Then
            kb=cshell%listshl(2,j)
            kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
            If (kb0 > config%natms) kb0=0
            !           Else If (cshell%listshl(2,j) == jb) Then
            !              kb=cshell%listshl(1,j)
            !              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
            !              If (kb0 > config%natms) kb0=0
          End If
        End Do

        If (ka*kb > 0) Then
          If (ia == ja) Then
            If (ib0 > 0) Then
              Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
              If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
            End If
            If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
            If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
              If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
            End If
            If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
            If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
          End If

          !           If (ib == ja) Then
          !              If (ia0 > 0) Then
          !                 Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
          !                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
          !              End If
          !              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
          !              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
          !           End If
          !
          !           If (ib == jb) Then
          !              If (ia0 > 0) Then
          !                 Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
          !                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
          !              End If
          !              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
          !           If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
          !           End If
        End If
      End Do

      ! add sites on basis of bonds to core-shell units

      Do kk=1,bond%n_types1
        If (.not. bond%restrained(bond%list(0,kk))) Then
          ja=bond%list(1,kk)
          jb=bond%list(2,kk)

          ja0=local_index(ja,config%nlast,config%lsi,config%lsa)
          jb0=local_index(jb,config%nlast,config%lsi,config%lsa)

          If (ja0 > config%natms) ja0=0
          If (jb0 > config%natms) jb0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0

          Do j=1,cshell%ntshl2
            If (cshell%listshl(1,j) == ja) Then
              ka=cshell%listshl(2,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            Else If (cshell%listshl(2,j) == ja) Then
              ka=cshell%listshl(1,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            End If

            If (cshell%listshl(1,j) == jb) Then
              kb=cshell%listshl(2,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            Else If (cshell%listshl(2,j) == jb) Then
              kb=cshell%listshl(1,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            End If
          End Do

          If (ka*kb > 0) Then
            If (ia == ja) Then
              If (ib0 > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (ia == jb) Then
              If (ib0 > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (ib == ja) Then
              If (ia0 > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (ib == jb) Then
              If (ia0 > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If
          End If
        End If
      End Do

      ! add sites on basis of valence angles to core-shell units

      Do kk=1,angle%n_types1
        If (.not. angle%restrained(angle%list(0,kk))) Then
          ja=angle%list(1,kk)
          jb=angle%list(2,kk)
          jc=angle%list(3,kk)

          ja0=local_index(ja,config%nlast,config%lsi,config%lsa)
          jb0=local_index(jb,config%nlast,config%lsi,config%lsa)
          jc0=local_index(jc,config%nlast,config%lsi,config%lsa)

          If (ja0 > config%natms) ja0=0
          If (jb0 > config%natms) jb0=0
          If (jc0 > config%natms) jc0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0
          kc = 0 ; kc0 = 0

          Do j=1,cshell%ntshl2
            If (cshell%listshl(1,j) == ja) Then
              ka=cshell%listshl(2,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            Else If (cshell%listshl(2,j) == ja) Then
              ka=cshell%listshl(1,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            End If

            If (cshell%listshl(1,j) == jb) Then
              kb=cshell%listshl(2,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            Else If (cshell%listshl(2,j) == jb) Then
              kb=cshell%listshl(1,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            End If

            If (cshell%listshl(1,j) == jc) Then
              kc=cshell%listshl(2,j)
              kc0=local_index(kc,config%nlast,config%lsi,config%lsa)
              If (kc0 > config%natms) kc0=0
            Else If (cshell%listshl(2,j) == jc) Then
              kc=cshell%listshl(1,j)
              kc0=local_index(kc,config%nlast,config%lsi,config%lsa)
              If (kc0 > config%natms) kc0=0
            End If
          End Do

          If (ia == ja) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jc) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ib == ja) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jb) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jc) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If
          End If
        End If
      End Do

      ! add sites on basis of dihedral angles to core-shell units

      Do kk=1,dihedral%n_types1
        If (.not. dihedral%restrained(dihedral%list(0,kk))) Then
          ja=dihedral%list(1,kk)
          jb=dihedral%list(2,kk)
          jc=dihedral%list(3,kk)
          jd=dihedral%list(4,kk)

          ja0=local_index(ja,config%nlast,config%lsi,config%lsa)
          jb0=local_index(jb,config%nlast,config%lsi,config%lsa)
          jc0=local_index(jc,config%nlast,config%lsi,config%lsa)
          jd0=local_index(jd,config%nlast,config%lsi,config%lsa)

          If (ja0 > config%natms) ja0=0
          If (jb0 > config%natms) jb0=0
          If (jc0 > config%natms) jc0=0
          If (jd0 > config%natms) jd0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0
          kc = 0 ; kc0 = 0
          kd = 0 ; kd0 = 0

          Do j=1,cshell%ntshl2
            If (cshell%listshl(1,j) == ja) Then
              ka=cshell%listshl(2,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            Else If (cshell%listshl(2,j) == ja) Then
              ka=cshell%listshl(1,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            End If

            If (cshell%listshl(1,j) == jb) Then
              kb=cshell%listshl(2,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            Else If (cshell%listshl(2,j) == jb) Then
              kb=cshell%listshl(1,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            End If

            If (cshell%listshl(1,j) == jc) Then
              kc=cshell%listshl(2,j)
              kc0=local_index(kc,config%nlast,config%lsi,config%lsa)
              If (kc0 > config%natms) kc0=0
            Else If (cshell%listshl(2,j) == jc) Then
              kc=cshell%listshl(1,j)
              kc0=local_index(kc,config%nlast,config%lsi,config%lsa)
              If (kc0 > config%natms) kc0=0
            End If

            If (cshell%listshl(1,j) == jd) Then
              kd=cshell%listshl(2,j)
              kd0=local_index(kd,config%nlast,config%lsi,config%lsa)
              If (kd0 > config%natms) kd0=0
            Else If (cshell%listshl(2,j) == jd) Then
              kd=cshell%listshl(1,j)
              kd0=local_index(kd,config%nlast,config%lsi,config%lsa)
              If (kd0 > config%natms) kd0=0
            End If
          End Do

          If (ia == ja) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jc) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jd) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ib == ja) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jb) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jc) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jd) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If
          End If
        End If
      End Do

      ! add sites on basis of inversion angles to core-shell units

      Do kk=1,inversion%n_types1
        If (.not. inversion%restrained(inversion%list(0,kk))) Then
          ja=inversion%list(1,kk)
          jb=inversion%list(2,kk)
          jc=inversion%list(3,kk)
          jd=inversion%list(4,kk)

          ja0=local_index(ja,config%nlast,config%lsi,config%lsa)
          jb0=local_index(jb,config%nlast,config%lsi,config%lsa)
          jc0=local_index(jc,config%nlast,config%lsi,config%lsa)
          jd0=local_index(jd,config%nlast,config%lsi,config%lsa)

          If (ja0 > config%natms) ja0=0
          If (jb0 > config%natms) jb0=0
          If (jc0 > config%natms) jc0=0
          If (jd0 > config%natms) jd0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0
          kc = 0 ; kc0 = 0
          kd = 0 ; kd0 = 0

          Do j=1,cshell%ntshl2
            If (cshell%listshl(1,j) == ja) Then
              ka=cshell%listshl(2,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            Else If (cshell%listshl(2,j) == ja) Then
              ka=cshell%listshl(1,j)
              ka0=local_index(ka,config%nlast,config%lsi,config%lsa)
              If (ka0 > config%natms) ka0=0
            End If

            If (cshell%listshl(1,j) == jb) Then
              kb=cshell%listshl(2,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            Else If (cshell%listshl(2,j) == jb) Then
              kb=cshell%listshl(1,j)
              kb0=local_index(kb,config%nlast,config%lsi,config%lsa)
              If (kb0 > config%natms) kb0=0
            End If

            If (cshell%listshl(1,j) == jc) Then
              kc=cshell%listshl(2,j)
              kc0=local_index(kc,config%nlast,config%lsi,config%lsa)
              If (kc0 > config%natms) kc0=0
            Else If (cshell%listshl(2,j) == jc) Then
              kc=cshell%listshl(1,j)
              kc0=local_index(kc,config%nlast,config%lsi,config%lsa)
              If (kc0 > config%natms) kc0=0
            End If

            If (cshell%listshl(1,j) == jd) Then
              kd=cshell%listshl(2,j)
              kd0=local_index(kd,config%nlast,config%lsi,config%lsa)
              If (kd0 > config%natms) kd0=0
            Else If (cshell%listshl(2,j) == jd) Then
              kd=cshell%listshl(1,j)
              kd0=local_index(kd,config%nlast,config%lsi,config%lsa)
              If (kd0 > config%natms) kd0=0
            End If
          End Do

          If (ia == ja) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jc) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ia == jd) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ib0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,mpoles%charmm)
            End If
          End If

          If (ib == ja) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jb) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jc) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kd,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,mpoles%charmm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,mpoles%charmm)
            End If
          End If

          If (ib == jd) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,ka,ibig,mpoles%charmm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kb,ibig,mpoles%charmm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,mpoles%charmm)
                Call add_exclusion(safe,ia0,kc,ibig,mpoles%charmm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,mpoles%charmm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,mpoles%charmm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,mpoles%charmm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,mpoles%charmm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,mpoles%charmm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,mpoles%charmm)
            End If
          End If
        End If
      End Do
    End Do

    ! check for exceeded array bounds

    Call gcheck(comm,safe)
    If (.not.safe) Then
      Call gmax(comm,ibig)
      Call warning(250,Real(ibig,wp),Real(max_exclude,wp),0.0_wp)
      Call error(65)
    End If

    ! EXCLUDE core sites mapped on the same RB unit

    Do i=1,config%natms                                                 ! on this node only (& below)
      l=mpoles%charmm(0,i)                                             ! end of list tag
      If (l > 0 .and. cshell%legshl(0,i) > 0) Then                     ! this is a qualifying CHARMMing core
        ibig=config%ltg(i)
        Do j=1,rigid%n_types_book                                          ! loop over the extended list of all RB
          k=rigid%list(-1,j)
          If (Any(rigid%list(1:k,j) == ibig)) Then               ! This core resides on a RB
            kk=l                                             ! running index
            Do While (kk > 0)                                ! run down to the beginning of the list
              If (Any(rigid%list(1:k,j) == mpoles%charmm(kk,i))) Then ! any other cores in this RB list
                If (kk < l) mpoles%charmm(kk,i)=mpoles%charmm(l,i)       ! swap with last "valid"
                mpoles%charmm(l,i)=0                              ! invalidate the last entry
                l=l-1                                      ! reduce "the end of list" tag
              End If
              kk=kk-1                                       ! reduce the running index
            End Do
            If (l < mpoles%charmm(0,i)) mpoles%charmm(0,i)=l               ! refresh the end of list tag
            If (l == 0) Exit                                 ! exit to main do loop
          End If
        End Do
      End If
    End Do

    ! EXCLUDE core sites on basis of frozen-frozen interactions

    Do i=1,config%natms                                     ! on this node only (& below)
      l=mpoles%charmm(0,i)                                 ! end of list tag
      If (l > 0 .and. config%lfrzn(i) > 0) Then            ! this is a qualifying CHARMMing core
        kk=l                                       ! running index
        Do While (kk > 0)                          ! run down to the beginning of the list
          k=mpoles%charmm(kk,i)
          j=local_index(k,config%nlast,config%lsi,config%lsa)
          If (config%lfrzn(j) > 0) Then                  ! any other frozen cores
            If (kk < l) mpoles%charmm(kk,i)=mpoles%charmm(l,i) ! swap with last "valid"
            mpoles%charmm(l,i)=0                        ! invalidate the last entry
            l=l-1                                ! reduce "the end of list" tag
          End If
          kk=kk-1                                 ! reduce the running index
        End Do
        If (l < mpoles%charmm(0,i)) mpoles%charmm(0,i)=l         ! refresh the end of list tag
      End If
    End Do

    ! sort mpoles%charmm

    kk=0
    Do i=1,config%natms
      j=mpoles%charmm(0,i)
      If (j > 0) Call shellsort(j,mpoles%charmm(1:j,i))
      kk=kk+j
    End Do
    Call gsum(comm,kk)

    Write(message,'(a,i11)') "Total of CHARMMing core-shell (scaled and dumped) cross-interatcions",kk/2
    Call info('',.true.)
    Call info(message,.true.)

    If (kk == 0) Then
      mpoles%key = POLARISATION_DEFAULT
      Call warning('CHARMM polarisation scheme unapplicable as no pair are detected',.true.)
      Deallocate (mpoles%charmm, Stat=fail)
    End If
  End Subroutine build_chrm_intra
End Module build_chrm
