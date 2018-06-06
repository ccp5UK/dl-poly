Module build_excl
! SETUP MODULES

  Use kinds, Only : wp
  Use comms,  Only : comms_type,gcheck,gmax
  Use setup

! CONFIG MODULE

  Use configuration, Only : natms,nlast,lsi,lsa,lexatm

! INTERACTION MODULES

  Use core_shell

  Use constraints

  Use rigid_bodies

  Use bonds, Only : bonds_type
  Use angles, Only : angles_type
  Use dihedrals, Only : dihedrals_type
  Use inversions, Only : inversions_type
  Use numerics, Only : local_index,shellsort

  Implicit None

  Private
  Public :: build_excl_intra, add_exclusion

Contains

Subroutine build_excl_intra(lecx,bond,angle,dihedral,inversion,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing the excluded pair interaction
! of the part of the MD system mapped onto this node
!
! keyinteraction < 0 exclusion does not apply
!
! copyright - daresbury laboratory
! author    - w.smith february 1999
! amended   - i.t.todorov july 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Logical,             Intent( In    ) :: lecx
  Type( bonds_type ), Intent( In    ) :: bond
  Type( angles_type ), Intent( In    ) :: angle
  Type( dihedrals_type ), Intent( In    ) :: dihedral
  Type( inversions_type ), Intent( In    ) :: inversion
  Type( comms_type ),  Intent( InOut ) :: comm

  Logical :: safe
  Integer :: fail(1:2),i,j,k,l,kk,             &
             ia,ib,ic,id,ia0,ib0,ic0,id0,      &
             ibig,ja,jb,jc,jd,ja0,jb0,jc0,jd0, &
             ka,kb,kc,kd,ka0,kb0,kc0,kd0

  Integer, Dimension( : , : ), Allocatable :: irgd,irgd0,jrgd,jrgd0
  Character ( Len = 256 )  ::  message

  If (mxrgd > 0) Then
     fail=0
     Allocate (irgd(1:mxlrgd,1:mxrgd),irgd0(1:mxlrgd,1:mxrgd), Stat=fail(1))
     Allocate (jrgd(1:mxlrgd,1:mxrgd),jrgd0(1:mxlrgd,1:mxrgd), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'build_excl_intra allocation failure'
        Call error(0,message)
     End If
  End If

! variables for array bound checking

  ibig=0
  safe=.true.

! exclude sites on basis of bond constraints

  Do i=1,ntcons
     ia=listcon(1,i)
     ib=listcon(2,i)

     ia0=local_index(ia,nlast,lsi,lsa)
     ib0=local_index(ib,nlast,lsi,lsa)

     If (ia0 > natms) ia0=0
     If (ib0 > natms) ib0=0

! add atoms to exclusion list

     If (ia0 > 0) Call add_exclusion(safe,ia0,ib,ibig,lexatm)
     If (ib0 > 0) Call add_exclusion(safe,ib0,ia,ibig,lexatm)
  End Do

! exclude sites on basis of RBs

  Do i=1,ntrgd
     Do j=1,listrgd(-1,i)
        irgd(j,i)=listrgd(j,i)
        irgd0(j,i)=local_index(listrgd(j,i),nlast,lsi,lsa)
        If (irgd0(j,i) > natms) irgd0(j,i)=0
     End Do

! add atoms to exclusion list

     Do j=1,listrgd(-1,i)
        If (irgd0(j,i) > 0) Then
           Do k=1,listrgd(-1,i)
              If (k /= j) Call add_exclusion(safe,irgd0(j,i),irgd(k,i),ibig,lexatm)
           End Do
        End If
     End Do
  End Do

! exclude sites on basis of chemical bonds

  Do i=1,bond%n_types
     If (bond%key(bond%list(0,i)) > 0) Then
        ia=bond%list(1,i)
        ib=bond%list(2,i)

        ia0=local_index(ia,nlast,lsi,lsa)
        ib0=local_index(ib,nlast,lsi,lsa)

        If (ia0 > natms) ia0=0
        If (ib0 > natms) ib0=0

! add atoms to exclusion list

        If (ia0 > 0) Call add_exclusion(safe,ia0,ib,ibig,lexatm)
        If (ib0 > 0) Call add_exclusion(safe,ib0,ia,ibig,lexatm)
     End If
  End Do

! exclude sites on basis of bond angles

  Do i=1,angle%n_types
     If (angle%key(angle%list(0,i)) > 0) Then
        ia=angle%list(1,i)
        ib=angle%list(2,i)
        ic=angle%list(3,i)


        ia0=local_index(ia,nlast,lsi,lsa)
        ib0=local_index(ib,nlast,lsi,lsa)
        ic0=local_index(ic,nlast,lsi,lsa)

        If (ia0 > natms) ia0=0
        If (ib0 > natms) ib0=0
        If (ic0 > natms) ic0=0

! add atoms to exclusion list

        If (ia0 > 0) Then ! ia : ib - ic interactions
           Call add_exclusion(safe,ia0,ib,ibig,lexatm)
           Call add_exclusion(safe,ia0,ic,ibig,lexatm)
        End If

        If (ib0 > 0) Then ! ib : ia - ic interactions
           Call add_exclusion(safe,ib0,ia,ibig,lexatm)
           Call add_exclusion(safe,ib0,ic,ibig,lexatm)
        End If

        If (ic0 > 0) Then ! ic : ia - ib interactions
           Call add_exclusion(safe,ic0,ia,ibig,lexatm)
           Call add_exclusion(safe,ic0,ib,ibig,lexatm)
        End If
     End If
  End Do

! exclude sites on basis of dihedral angles

  Do i=1,dihedral%n_types
     If (dihedral%key(dihedral%list(0,i)) > 0) Then
        ia=dihedral%list(1,i)
        ib=dihedral%list(2,i)
        ic=dihedral%list(3,i)
        id=dihedral%list(4,i)

        ia0=local_index(ia,nlast,lsi,lsa)
        ib0=local_index(ib,nlast,lsi,lsa)
        ic0=local_index(ic,nlast,lsi,lsa)
        id0=local_index(id,nlast,lsi,lsa)

        If (ia0 > natms) ia0=0
        If (ib0 > natms) ib0=0
        If (ic0 > natms) ic0=0
        If (id0 > natms) id0=0

! add atoms to exclusion list

        If (ia0 > 0) Then ! ia : ib - ic - id interactions
           Call add_exclusion(safe,ia0,ib,ibig,lexatm)
           Call add_exclusion(safe,ia0,ic,ibig,lexatm)
           Call add_exclusion(safe,ia0,id,ibig,lexatm)
        End If

        If (ib0 > 0) Then ! ib : ia - ic - id interactions
           Call add_exclusion(safe,ib0,ia,ibig,lexatm)
           Call add_exclusion(safe,ib0,ic,ibig,lexatm)
           Call add_exclusion(safe,ib0,id,ibig,lexatm)
        End If

        If (ic0 > 0) Then ! ic : ia - ib - id interactions
           Call add_exclusion(safe,ic0,ia,ibig,lexatm)
           Call add_exclusion(safe,ic0,ib,ibig,lexatm)
           Call add_exclusion(safe,ic0,id,ibig,lexatm)
        End If

        If (id0 > 0) Then ! id : ia - ib - ic interactions
           Call add_exclusion(safe,id0,ia,ibig,lexatm)
           Call add_exclusion(safe,id0,ib,ibig,lexatm)
           Call add_exclusion(safe,id0,ic,ibig,lexatm)
        End If
     End If
  End Do

! exclude sites on basis of inversion potentials

  Do i=1,inversion%n_types
     If (inversion%key(inversion%list(0,i)) > 0) Then
        ia=inversion%list(1,i)
        ib=inversion%list(2,i)
        ic=inversion%list(3,i)
        id=inversion%list(4,i)

        ia0=local_index(ia,nlast,lsi,lsa)
        ib0=local_index(ib,nlast,lsi,lsa)
        ic0=local_index(ic,nlast,lsi,lsa)
        id0=local_index(id,nlast,lsi,lsa)

        If (ia0 > natms) ia0=0
        If (ib0 > natms) ib0=0
        If (ic0 > natms) ic0=0
        If (id0 > natms) id0=0

! add atoms to exclusion list

        If (ia0 > 0) Then ! ia : ib - ic - id interactions
           Call add_exclusion(safe,ia0,ib,ibig,lexatm)
           Call add_exclusion(safe,ia0,ic,ibig,lexatm)
           Call add_exclusion(safe,ia0,id,ibig,lexatm)
        End If

        If (ib0 > 0) Then ! ib : ia - ic - id interactions
           Call add_exclusion(safe,ib0,ia,ibig,lexatm)
           Call add_exclusion(safe,ib0,ic,ibig,lexatm)
           Call add_exclusion(safe,ib0,id,ibig,lexatm)
        End If

        If (ic0 > 0) Then ! ic : ia - ib - id interactions
           Call add_exclusion(safe,ic0,ia,ibig,lexatm)
           Call add_exclusion(safe,ic0,ib,ibig,lexatm)
           Call add_exclusion(safe,ic0,id,ibig,lexatm)
        End If

        If (id0 > 0) Then ! id : ia - ib - ic interactions
           Call add_exclusion(safe,id0,ia,ibig,lexatm)
           Call add_exclusion(safe,id0,ib,ibig,lexatm)
           Call add_exclusion(safe,id0,ic,ibig,lexatm)
        End If
     End If
  End Do

! exclude sites on basis of core-shell units

  Do i=1,ntshl2
     ia=listshl(1,i) ! This is the core
     ib=listshl(2,i) ! This is the shell

     ia0=local_index(ia,nlast,lsi,lsa)
     ib0=local_index(ib,nlast,lsi,lsa)

     If (ia0 > natms) ia0=0
     If (ib0 > natms) ib0=0

! add atoms to exclusion list

     If (ia0 > 0) Call add_exclusion(safe,ia0,ib,ibig,lexatm)
     If (ib0 > 0) Call add_exclusion(safe,ib0,ia,ibig,lexatm)

! exclude sites on basis of constraint bonds to core-shell units

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

        If (lecx) Then
           Do j=1,ntshl2
              If (listshl(1,j) == ja) Then
                 ka=listshl(2,j)
                 ka0=local_index(ka,nlast,lsi,lsa)
                 If (ka0 > natms) ka0=0
!              Else If (listshl(2,j) == ja) Then
!                 ka=listshl(1,j)
!                 ka0=local_index(ka,nlast,lsi,lsa)
!                 If (ka0 > natms) ka0=0
              End If

              If (listshl(1,j) == jb) Then
                 kb=listshl(2,j)
                 kb0=local_index(kb,nlast,lsi,lsa)
                 If (kb0 > natms) kb0=0
!              Else If (listshl(2,j) == jb) Then
!                 kb=listshl(1,j)
!                 kb0=local_index(kb,nlast,lsi,lsa)
!                 If (kb0 > natms) kb0=0
              End If
           End Do
        End If

        If (ia == ja) Then
           If (ib0 > 0) Then
              Call add_exclusion(safe,ib0,jb,ibig,lexatm)
              If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)
           End If
           If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
           If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)
        End If

        If (ia == jb) Then
           If (ib0 > 0) Then
              Call add_exclusion(safe,ib0,ja,ibig,lexatm)
              If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)
           End If
           If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
           If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)
        End If

!        If (ib == ja) Then
!           If (ia0 > 0) Then
!              Call add_exclusion(safe,ia0,jb,ibig,lexatm)
!              If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)
!           End If
!           If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
!           If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)
!        End If
!
!        If (ib == jb) Then
!           If (ia0 > 0) Then
!              Call add_exclusion(safe,ia0,ja,ibig,lexatm)
!              If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)
!           End If
!           If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
!           If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)
!        End If
     End Do

! exclude sites on basis of RBs to core-shell units

     Do j=1,ntrgd1
        Do k=1,listrgd(-1,j)
           irgd(k,j)=listrgd(k,j)
           irgd0(k,j)=local_index(listrgd(k,j),nlast,lsi,lsa)
           If (irgd0(k,j) > natms) irgd0(k,j)=0
        End Do

! Check for adjacent core or shell so that all possible, even between
! adjacent particles in any intra-unit, interactions are excluded
! The checks commented out are not needed due to preventative checks
! in read_field by which shells are not allowed to be frozen,
! constrained bonded, RBed or tethered.

        jrgd = 0 ; jrgd0 = 0

        If (lecx) Then
           Do kk=1,ntshl2
              Do k=1,listrgd(-1,j)
                 If (irgd(k,j) == listshl(1,kk)) Then
                    jrgd(k,j)=listshl(2,kk)
                    jrgd0(k,j)=local_index(jrgd(k,j),nlast,lsi,lsa)
                    If (jrgd0(k,j) > natms) jrgd0(k,j)=0
!                 Else If (irgd(k,j) == listshl(2,kk)) Then
!                    jrgd(k,j)=listshl(2,kk)
!                    jrgd0(k,j)=local_index(jrgd(k,j),nlast,lsi,lsa)
!                    If (jrgd0(k,j) > natms) jrgd0(k,j)=0
                 End If
              End Do
           End Do
        End If

        Do k=1,listrgd(-1,j)
           If (irgd(k,j) == ia) Then
              Do l=1,listrgd(-1,j)
                 If (l /= k) Then
                    If (ib0 > 0) Then
                       Call add_exclusion(safe,ib0,irgd(l,j),ibig,lexatm)
                       If (jrgd(l,j) > 0) Call add_exclusion(safe,ib0,jrgd(l,j),ibig,lexatm)
                    End If
                    If (irgd0(l,j) > 0) Call add_exclusion(safe,irgd0(l,j),ib,ibig,lexatm)
                    If (jrgd0(l,j) > 0) Call add_exclusion(safe,jrgd0(l,j),ib,ibig,lexatm)
                 End If
              End Do
           End If

!           If (irgd(k,j) == ib) Then
!              Do l=1,listrgd(-1,j)
!                 If (l /= k) Then
!                    If (ia0 > 0) Then
!                       Call add_exclusion(safe,ia0,irgd(l,j),ibig,lexatm)
!                       If (jrgd(l,j) > 0) Call add_exclusion(safe,ia0,jrgd(l,j),ibig,lexatm)
!                    End If
!                    If (irgd0(l,j) > 0) Call add_exclusion(safe,irgd0(l,j),ia,ibig,lexatm)
!                    If (jrgd0(l,j) > 0) Call add_exclusion(safe,jrgd0(l,j),ia,ibig,lexatm)
!                 End If
!              End Do
!           End If
        End Do
     End Do

! exclude sites on basis of bonds to core-shell units

     Do kk=1,bond%n_types1
        If (bond%key(bond%list(0,kk)) > 0) Then
           ja=bond%list(1,kk)
           jb=bond%list(2,kk)

           ja0=local_index(ja,nlast,lsi,lsa)
           jb0=local_index(jb,nlast,lsi,lsa)

           If (ja0 > natms) ja0=0
           If (jb0 > natms) jb0=0

! Check for adjacent core or shell so that all possible, even between
! adjacent particles in any intra-unit, interactions are excluded

           ka = 0 ; ka0 = 0
           kb = 0 ; kb0 = 0

           If (lecx) Then
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
           End If

           If (ia == ja) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)
              End If
              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)
           End If

           If (ia == jb) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)
              End If
              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)
           End If

           If (ib == ja) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)
              End If
              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)
           End If

           If (ib == jb) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)
              End If
              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)
           End If
        End If
     End Do

! exclude sites on basis of valence angles to core-shell units

     Do kk=1,angle%n_types1
        If (angle%key(angle%list(0,kk)) > 0) Then
           ja=angle%list(1,kk)
           jb=angle%list(2,kk)
           jc=angle%list(3,kk)

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

           If (lecx) Then
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
           End If

           If (ia == ja) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)
              End If

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)
           End If

           If (ia == jb) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)
           End If

           If (ia == jc) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)
           End If

           If (ib == ja) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)
              End If

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)
           End If

           If (ib == jb) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)
           End If

           If (ib == jc) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)
           End If
        End If
     End Do

! exclude sites on basis of dihedral angles to core-shell units

     Do kk=1,dihedral%n_types1
        If (dihedral%key(dihedral%list(0,kk)) > 0) Then
           ja=dihedral%list(1,kk)
           jb=dihedral%list(2,kk)
           jc=dihedral%list(3,kk)
           jd=dihedral%list(4,kk)

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

           If (lecx) Then
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
           End If

           If (ia == ja) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ib0,kd,ibig,lexatm)
              End If

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lexatm)
           End If

           If (ia == jb) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ib0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lexatm)
           End If

           If (ia == jc) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ib0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lexatm)
           End If

           If (ia == jd) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)
           End If

           If (ib == ja) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ia0,kd,ibig,lexatm)
              End If

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lexatm)
           End If

           If (ib == jb) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ia0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lexatm)
           End If

           If (ib == jc) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ia0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lexatm)
           End If

           If (ib == jd) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)
           End If
        End If
     End Do

! exclude sites on basis of inversion angles to core-shell units

     Do kk=1,inversion%n_types1
        If (inversion%key(inversion%list(0,kk)) > 0) Then
           ja=inversion%list(1,kk)
           jb=inversion%list(2,kk)
           jc=inversion%list(3,kk)
           jd=inversion%list(4,kk)

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

           If (lecx) Then
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
           End If

           If (ia == ja) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ib0,kd,ibig,lexatm)
              End If

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lexatm)
           End If

           If (ia == jb) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ib0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lexatm)
           End If

           If (ia == jc) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ib0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ib,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ib,ibig,lexatm)
           End If

           If (ia == jd) Then
              If (ib0 > 0) Then
                 Call add_exclusion(safe,ib0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ib0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ib0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ib0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ib0,kc,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ib,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ib,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ib,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ib,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ib,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ib,ibig,lexatm)
           End If

           If (ib == ja) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ia0,kd,ibig,lexatm)
              End If

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lexatm)
           End If

           If (ib == jb) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ia0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lexatm)
           End If

           If (ib == jc) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jd,ibig,lexatm)
                 If (kd > 0) Call add_exclusion(safe,ia0,kd,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jd0 > 0) Call add_exclusion(safe,jd0,ia,ibig,lexatm)
              If (kd0 > 0) Call add_exclusion(safe,kd0,ia,ibig,lexatm)
           End If

           If (ib == jd) Then
              If (ia0 > 0) Then
                 Call add_exclusion(safe,ia0,ja,ibig,lexatm)
                 If (ka > 0) Call add_exclusion(safe,ia0,ka,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jb,ibig,lexatm)
                 If (kb > 0) Call add_exclusion(safe,ia0,kb,ibig,lexatm)

                 Call add_exclusion(safe,ia0,jc,ibig,lexatm)
                 If (kc > 0) Call add_exclusion(safe,ia0,kc,ibig,lexatm)
              End If

              If (ja0 > 0) Call add_exclusion(safe,ja0,ia,ibig,lexatm)
              If (ka0 > 0) Call add_exclusion(safe,ka0,ia,ibig,lexatm)

              If (jb0 > 0) Call add_exclusion(safe,jb0,ia,ibig,lexatm)
              If (kb0 > 0) Call add_exclusion(safe,kb0,ia,ibig,lexatm)

              If (jc0 > 0) Call add_exclusion(safe,jc0,ia,ibig,lexatm)
              If (kc0 > 0) Call add_exclusion(safe,kc0,ia,ibig,lexatm)
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

! sort lexatm

  Do i=1,natms
     j=lexatm(0,i)
     If (j > 0) Call shellsort(j,lexatm(1:j,i))
  End Do

  If (mxrgd > 0) Then
     Deallocate (irgd,irgd0, Stat=fail(1))
     Deallocate (jrgd,jrgd0, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'build_excl_intra deallocation failure'
        Call error(0,message)
     End If
  End If
End Subroutine build_excl_intra

Subroutine add_exclusion(safe,ia0,ib,ibig,lexatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to add excluded atoms to the excluded atom list
! provided they are not already excluded
!
! copyright - daresbury laboratory
! author    - w.smith march 1999
! amended   - i.t.todorov july 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,                                  Intent( In    ) :: ia0,ib
  Logical,                                  Intent( InOut ) :: safe
  Integer,                                  Intent( InOut ) :: ibig
  Integer, Dimension( 0:mxexcl, 1:mxatdm ), Intent( InOut ) :: lexatm

  Logical :: safe_local,l_excluded
  Integer :: last

! Get current length

  last = lexatm(0,ia0)

! Determine possible exclusion

  l_excluded = Any(lexatm(1:last,ia0) == ib)

  If (.not.l_excluded) Then

! Get local safety no array overflow

     safe_local = (last < mxexcl-1)

! Determine global safety

     safe = safe .and. safe_local

     If (safe_local) Then

! Increase length of the ia0 exclusion list and tag ib in it

        last = last + 1
        lexatm(0,ia0) = last
        lexatm(last,ia0) = ib

        ibig=Max(ibig,last)

     Else

! Collect number of offences

        lexatm(mxexcl,ia0) = lexatm(mxexcl,ia0) + 1

        ibig=Max(ibig,last+lexatm(mxexcl,ia0))

     End If

  End If

End Subroutine add_exclusion

End Module build_excl
