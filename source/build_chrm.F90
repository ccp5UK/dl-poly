Module build_chrm
  ! SETUP MODULES

  Use kinds, only : wp
  Use comms,  Only : comms_type,gcheck,gmax,gsum
  Use setup_module

  ! CONFIG MODULE

  Use configuration, Only : natms,nlast,lsi,lsa,ltg,lfrzn

  ! INTERACTION MODULES

  Use core_shell

  Use constraints

  Use rigid_bodies_module

  Use bonds
  Use angles
  Use dihedrals_module
  Use inversions_module

  ! MULTIPOLES MODULE

  Use mpoles, Only : keyind,lchatm ! equivalent to lexatm in configuration

  Implicit None

  Private

  Public :: build_chrm_intra
Contains 
  Subroutine build_chrm_intra(comm)

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
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( comms_type ),  Intent( InOut) :: comm

    Logical :: safe
    Integer :: fail,i,j,k,l,kk,                  &
      ia,ib,ic,id,ia0,ib0,ic0,id0,      &
      ibig,ja,jb,jc,jd,ja0,jb0,jc0,jd0, &
      ka,kb,kc,kd,ka0,kb0,kc0,kd0,local_index

    ! variables for array bound checking

    ibig=0
    safe=.true.

    ! go over the extended list of core-shell units

    Do i=1,ntshl2
      ia=listshl(1,i) ! This is the core
      ib=listshl(2,i) ! This is the shell

      ia0=local_index(ia,nlast,lsi,lsa)
      ib0=local_index(ib,nlast,lsi,lsa)

      If (ia0 > natms) ia0=0
      If (ib0 > natms) ib0=0

      ! add sites on basis of constraint bonds to core-shell units

      Do kk=1,ntcons1
        ja=listcon(1,kk)
        jb=listcon(2,kk)

        ja0=local_index(ja,nlast,lsi,lsa)
        jb0=local_index(jb,nlast,lsi,lsa)

        If (ja0 > natms) ja0=0
        If (jb0 > natms) jb0=0

        ! Check for adjacent core or shell so that all possible, even between
        ! adjacent particles in any intra-unit, interactions are excluded.
        ! The checks commented out are not needed due to preventative checks
        ! in read_field by which shells are not allowed to be frozen,
        ! constrained bonded, RBed or tethered.

        ka = 0 ; ka0 = 0
        kb = 0 ; kb0 = 0

        Do j=1,ntshl2
          If (listshl(1,j) == ja) Then
            ka=listshl(2,j)
            ka0=local_index(ka,nlast,lsi,lsa)
            If (ka0 > natms) ka0=0
            !           Else If (listshl(2,j) == ja) Then
            !              ka=listshl(1,j)
            !              ka0=local_index(ka,nlast,lsi,lsa)
            !              If (ka0 > natms) ka0=0
          End If

          If (listshl(1,j) == jb) Then
            kb=listshl(2,j)
            kb0=local_index(kb,nlast,lsi,lsa)
            If (kb0 > natms) kb0=0
            !           Else If (listshl(2,j) == jb) Then
            !              kb=listshl(1,j)
            !              kb0=local_index(kb,nlast,lsi,lsa)
            !              If (kb0 > natms) kb0=0
          End If
        End Do

        If (ka*kb > 0) Then
          If (ia == ja) Then
            If (ib0 > 0) Then
              Call add_exclusion(safe,ib0,jb,ibig,lchatm)
              If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lchatm)
            End If
            If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
            If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              Call add_exclusion(safe,ib0,ja,ibig,lchatm)
              If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lchatm)
            End If
            If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
            If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
          End If

          !           If (ib == ja) Then
          !              If (ia0 > 0) Then
          !                 Call add_exclusion(safe,ia0,jb,ibig,lchatm)
          !                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lchatm)
          !              End If
          !              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
          !              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
          !           End If
          !
          !           If (ib == jb) Then
          !              If (ia0 > 0) Then
          !                 Call add_exclusion(safe,ia0,ja,ibig,lchatm)
          !                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lchatm)
          !              End If
          !              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
          !           If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
          !           End If
        End If
      End Do

      ! add sites on basis of bonds to core-shell units

      Do kk=1,ntbond1
        If (keybnd(listbnd(0,kk)) > 0) Then
          ja=listbnd(1,kk)
          jb=listbnd(2,kk)

          ja0=local_index(ja,nlast,lsi,lsa)
          jb0=local_index(jb,nlast,lsi,lsa)

          If (ja0 > natms) ja0=0
          If (jb0 > natms) jb0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0

          Do j=1,ntshl2
            If (listshl(1,j) == ja) Then
              ka=listshl(2,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            Else If (listshl(2,j) == ja) Then
              ka=listshl(1,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            End If

            If (listshl(1,j) == jb) Then
              kb=listshl(2,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            Else If (listshl(2,j) == jb) Then
              kb=listshl(1,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            End If
          End Do

          If (ka*kb > 0) Then
            If (ia == ja) Then
              If (ib0 > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (ia == jb) Then
              If (ib0 > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (ib == ja) Then
              If (ia0 > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (ib == jb) Then
              If (ia0 > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If
          End If
        End If
      End Do

      ! add sites on basis of valence angles to core-shell units

      Do kk=1,ntangl1
        If (keyang(listang(0,kk)) > 0) Then
          ja=listang(1,kk)
          jb=listang(2,kk)
          jc=listang(3,kk)

          ja0=local_index(ja,nlast,lsi,lsa)
          jb0=local_index(jb,nlast,lsi,lsa)
          jc0=local_index(jc,nlast,lsi,lsa)

          If (ja0 > natms) ja0=0
          If (jb0 > natms) jb0=0
          If (jc0 > natms) jc0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0
          kc = 0 ; kc0 = 0

          Do j=1,ntshl2
            If (listshl(1,j) == ja) Then
              ka=listshl(2,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            Else If (listshl(2,j) == ja) Then
              ka=listshl(1,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            End If

            If (listshl(1,j) == jb) Then
              kb=listshl(2,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            Else If (listshl(2,j) == jb) Then
              kb=listshl(1,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            End If

            If (listshl(1,j) == jc) Then
              kc=listshl(2,j)
              kc0=local_index(kc,nlast,lsi,lsa)
              If (kc0 > natms) kc0=0
            Else If (listshl(2,j) == jc) Then
              kc=listshl(1,j)
              kc0=local_index(kc,nlast,lsi,lsa)
              If (kc0 > natms) kc0=0
            End If
          End Do

          If (ia == ja) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jc) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If
          End If

          If (ib == ja) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jb) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jc) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If
          End If
        End If
      End Do

      ! add sites on basis of dihedral angles to core-shell units

      Do kk=1,ntdihd1
        If (keydih(listdih(0,kk)) > 0) Then
          ja=listdih(1,kk)
          jb=listdih(2,kk)
          jc=listdih(3,kk)
          jd=listdih(4,kk)

          ja0=local_index(ja,nlast,lsi,lsa)
          jb0=local_index(jb,nlast,lsi,lsa)
          jc0=local_index(jc,nlast,lsi,lsa)
          jd0=local_index(jd,nlast,lsi,lsa)

          If (ja0 > natms) ja0=0
          If (jb0 > natms) jb0=0
          If (jc0 > natms) jc0=0
          If (jd0 > natms) jd0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0
          kc = 0 ; kc0 = 0
          kd = 0 ; kd0 = 0

          Do j=1,ntshl2
            If (listshl(1,j) == ja) Then
              ka=listshl(2,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            Else If (listshl(2,j) == ja) Then
              ka=listshl(1,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            End If

            If (listshl(1,j) == jb) Then
              kb=listshl(2,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            Else If (listshl(2,j) == jb) Then
              kb=listshl(1,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            End If

            If (listshl(1,j) == jc) Then
              kc=listshl(2,j)
              kc0=local_index(kc,nlast,lsi,lsa)
              If (kc0 > natms) kc0=0
            Else If (listshl(2,j) == jc) Then
              kc=listshl(1,j)
              kc0=local_index(kc,nlast,lsi,lsa)
              If (kc0 > natms) kc0=0
            End If

            If (listshl(1,j) == jd) Then
              kd=listshl(2,j)
              kd0=local_index(kd,nlast,lsi,lsa)
              If (kd0 > natms) kd0=0
            Else If (listshl(2,j) == jd) Then
              kd=listshl(1,j)
              kd0=local_index(kd,nlast,lsi,lsa)
              If (kd0 > natms) kd0=0
            End If
          End Do

          If (ia == ja) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,lchatm)
                Call add_exclusion(safe,ib0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jc) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,lchatm)
                Call add_exclusion(safe,ib0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jd) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If
          End If

          If (ib == ja) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jb) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,lchatm)
                Call add_exclusion(safe,ia0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jc) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,lchatm)
                Call add_exclusion(safe,ia0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jd) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If
          End If
        End If
      End Do

      ! add sites on basis of inversion angles to core-shell units

      Do kk=1,ntinv1
        If (keyinv(listinv(0,kk)) > 0) Then
          ja=listinv(1,kk)
          jb=listinv(2,kk)
          jc=listinv(3,kk)
          jd=listinv(4,kk)

          ja0=local_index(ja,nlast,lsi,lsa)
          jb0=local_index(jb,nlast,lsi,lsa)
          jc0=local_index(jc,nlast,lsi,lsa)
          jd0=local_index(jd,nlast,lsi,lsa)

          If (ja0 > natms) ja0=0
          If (jb0 > natms) jb0=0
          If (jc0 > natms) jc0=0
          If (jd0 > natms) jd0=0

          ! Check for adjacent core or shell so that all possible, even between
          ! adjacent particles in any intra-unit, interactions are excluded

          ka = 0 ; ka0 = 0
          kb = 0 ; kb0 = 0
          kc = 0 ; kc0 = 0
          kd = 0 ; kd0 = 0

          Do j=1,ntshl2
            If (listshl(1,j) == ja) Then
              ka=listshl(2,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            Else If (listshl(2,j) == ja) Then
              ka=listshl(1,j)
              ka0=local_index(ka,nlast,lsi,lsa)
              If (ka0 > natms) ka0=0
            End If

            If (listshl(1,j) == jb) Then
              kb=listshl(2,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            Else If (listshl(2,j) == jb) Then
              kb=listshl(1,j)
              kb0=local_index(kb,nlast,lsi,lsa)
              If (kb0 > natms) kb0=0
            End If

            If (listshl(1,j) == jc) Then
              kc=listshl(2,j)
              kc0=local_index(kc,nlast,lsi,lsa)
              If (kc0 > natms) kc0=0
            Else If (listshl(2,j) == jc) Then
              kc=listshl(1,j)
              kc0=local_index(kc,nlast,lsi,lsa)
              If (kc0 > natms) kc0=0
            End If

            If (listshl(1,j) == jd) Then
              kd=listshl(2,j)
              kd0=local_index(kd,nlast,lsi,lsa)
              If (kd0 > natms) kd0=0
            Else If (listshl(2,j) == jd) Then
              kd=listshl(1,j)
              kd0=local_index(kd,nlast,lsi,lsa)
              If (kd0 > natms) kd0=0
            End If
          End Do

          If (ia == ja) Then
            If (ib0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,lchatm)
                Call add_exclusion(safe,ib0,kd,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jb) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,lchatm)
                Call add_exclusion(safe,ib0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jc) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ib0,jd,ibig,lchatm)
                Call add_exclusion(safe,ib0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lchatm)
            End If
          End If

          If (ia == jd) Then
            If (ib0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ib0,ja,ibig,lchatm)
                Call add_exclusion(safe,ib0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ib0,jb,ibig,lchatm)
                Call add_exclusion(safe,ib0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ib0,jc,ibig,lchatm)
                Call add_exclusion(safe,ib0,kc,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lchatm)
            End If
          End If

          If (ib == ja) Then
            If (ia0 > 0) Then
              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,lchatm)
                Call add_exclusion(safe,ia0,kd,ibig,lchatm)
              End If
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jb) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,lchatm)
                Call add_exclusion(safe,ia0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jc) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kd > 0) Then
                Call add_exclusion(safe,ia0,jd,ibig,lchatm)
                Call add_exclusion(safe,ia0,kd,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kd > 0) Then
              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lchatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lchatm)
            End If
          End If

          If (ib == jd) Then
            If (ia0 > 0) Then
              If (ka > 0) Then
                Call add_exclusion(safe,ia0,ja,ibig,lchatm)
                Call add_exclusion(safe,ia0,ka,ibig,lchatm)
              End If

              If (kb > 0) Then
                Call add_exclusion(safe,ia0,jb,ibig,lchatm)
                Call add_exclusion(safe,ia0,kb,ibig,lchatm)
              End If

              If (kc > 0) Then
                Call add_exclusion(safe,ia0,jc,ibig,lchatm)
                Call add_exclusion(safe,ia0,kc,ibig,lchatm)
              End If
            End If

            If (ka > 0) Then
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lchatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lchatm)
            End If

            If (kb > 0) Then
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lchatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lchatm)
            End If

            If (kc > 0) Then
              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lchatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lchatm)
            End If
          End If
        End If
      End Do
    End Do

    ! check for exceeded array bounds

    Call gcheck(comm,safe)
    If (.not.safe) Then
      Call gmax(comm,ibig)
      Call warning(250,Real(ibig,wp),Real(mxexcl,wp),0.0_wp)
      Call error(65)
    End If

    ! EXCLUDE core sites mapped on the same RB unit

    Do i=1,natms                                                 ! on this node only (& below)
      l=lchatm(0,i)                                             ! end of list tag
      If (l > 0 .and. legshl(0,i) > 0) Then                     ! this is a qualifying CHARMMing core
        ibig=ltg(i)
        Do j=1,ntrgd1                                          ! loop over the extended list of all RB
          k=listrgd(-1,j)
          If (Any(listrgd(1:k,j) == ibig)) Then               ! This core resides on a RB
            kk=l                                             ! running index
            Do While (kk > 0)                                ! run down to the beginning of the list
              If (Any(listrgd(1:k,j) == lchatm(kk,i))) Then ! any other cores in this RB list
                If (kk < l) lchatm(kk,i)=lchatm(l,i)       ! swap with last "valid"
                lchatm(l,i)=0                              ! invalidate the last entry
                l=l-1                                      ! reduce "the end of list" tag
              End If
              kk=kk-1                                       ! reduce the running index
            End Do
            If (l < lchatm(0,i)) lchatm(0,i)=l               ! refresh the end of list tag
            If (l == 0) Exit                                 ! exit to main do loop
          End If
        End Do
      End If
    End Do

    ! EXCLUDE core sites on basis of frozen-frozen interactions

    Do i=1,natms                                     ! on this node only (& below)
      l=lchatm(0,i)                                 ! end of list tag
      If (l > 0 .and. lfrzn(i) > 0) Then            ! this is a qualifying CHARMMing core
        kk=l                                       ! running index
        Do While (kk > 0)                          ! run down to the beginning of the list
          k=lchatm(kk,i)
          j=local_index(k,nlast,lsi,lsa)
          If (lfrzn(j) > 0) Then                  ! any other frozen cores
            If (kk < l) lchatm(kk,i)=lchatm(l,i) ! swap with last "valid"
            lchatm(l,i)=0                        ! invalidate the last entry
            l=l-1                                ! reduce "the end of list" tag
          End If
          kk=kk-1                                 ! reduce the running index
        End Do
        If (l < lchatm(0,i)) lchatm(0,i)=l         ! refresh the end of list tag
      End If
    End Do

    ! sort lchatm

    kk=0
    Do i=1,natms
      j=lchatm(0,i)
      If (j > 0) Call shellsort(j,lchatm(1:j,i))
      kk=kk+j
    End Do
    Call gsum(comm,kk)

    If (comm%idnode == 0) Write(nrite,'(/,1x,a,i11)') "Total of CHARMMing core-shell (scaled and dumped) cross-interatcions",kk/2

    If (kk == 0) Then
      keyind = 0
      If (comm%idnode == 0) Then 
        Write(nrite,"(1x,a)") "*** warning - CHARMM polarisation scheme unapplicable as no pair are detected !!! ***"
      End If 
      Deallocate (lchatm, Stat=fail)
    End If

  Contains

    Subroutine add_exclusion(safe,ia0,ib,ibig,lchatm)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! dl_poly_4 subroutine to add CHARMMing core-shell atoms to the
      ! CHARMMing core-shell atoms list provided they are not already included
      !
      ! copyright - daresbury laboratory
      ! author    - w.smith march 1999
      ! amended   - i.t.todorov december 2016
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Use setup_module, Only : mxatdm,mxexcl

      Implicit None

      Integer,                                  Intent( In    ) :: ia0,ib
      Logical,                                  Intent( InOut ) :: safe
      Integer,                                  Intent( InOut ) :: ibig
      Integer, Dimension( 0:mxexcl, 1:mxatdm ), Intent( InOut ) :: lchatm

      Logical :: safe_local,l_excluded
      Integer :: last

      ! Get current length

      last = lchatm(0,ia0)

      ! Determine possible exclusion

      l_excluded = Any(lchatm(1:last,ia0) == ib)

      If (.not.l_excluded) Then

        ! Get local safety no array overflow

        safe_local = (last < mxexcl-1)

        ! Determine global safety

        safe = safe .and. safe_local

        If (safe_local) Then

          ! Increase length of the ia0 exclusion list and tag ib in it

          last = last + 1
          lchatm(0,ia0) = last
          lchatm(last,ia0) = ib

          ibig=Max(ibig,last)

        Else

          ! Collect number of offences

          lchatm(mxexcl,ia0) = lchatm(mxexcl,ia0) + 1

          ibig=Max(ibig,last+lchatm(mxexcl,ia0))

        End If

      End If

    End Subroutine add_exclusion

  End Subroutine build_chrm_intra
End Module build_chrm
