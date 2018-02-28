Subroutine build_excl_intra(lecx)

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

! SETUP MODULES

  Use kinds, only : wp
  Use comms_module,  Only : idnode,mxnode,gcheck,gmax
  Use setup_module

! CONFIG MODULE

  Use config_module, Only : natms,nlast,lsi,lsa,lexatm

! INTERACTION MODULES

  Use core_shell_module

  Use constraints_module

  Use rigid_bodies_module

  Use bonds_module
  Use angles_module
  Use dihedrals_module
  Use inversions_module

  Implicit None

  Logical, Intent( In    ) :: lecx

  Logical :: safe
  Integer :: fail(1:2),i,j,k,l,kk,             &
             ia,ib,ic,id,ia0,ib0,ic0,id0,      &
             ibig,ja,jb,jc,jd,ja0,jb0,jc0,jd0, &
             ka,kb,kc,kd,ka0,kb0,kc0,kd0,local_index

  Integer, Dimension( : , : ), Allocatable :: irgd,irgd0,jrgd,jrgd0

  If (mxrgd > 0) Then
     fail=0
     Allocate (irgd(1:mxlrgd,1:mxrgd),irgd0(1:mxlrgd,1:mxrgd), Stat=fail(1))
     Allocate (jrgd(1:mxlrgd,1:mxrgd),jrgd0(1:mxlrgd,1:mxrgd), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'build_excl_intra allocation failure, node: ', idnode
        Call error(0)
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

  Do i=1,ntbond
     If (keybnd(listbnd(0,i)) > 0) Then
        ia=listbnd(1,i)
        ib=listbnd(2,i)

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

  Do i=1,ntangl
     If (keyang(listang(0,i)) > 0) Then
        ia=listang(1,i)
        ib=listang(2,i)
        ic=listang(3,i)


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

  Do i=1,ntdihd
     If (keydih(listdih(0,i)) > 0) Then
        ia=listdih(1,i)
        ib=listdih(2,i)
        ic=listdih(3,i)
        id=listdih(4,i)

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

  Do i=1,ntinv
     If (keyinv(listinv(0,i)) > 0) Then
        ia=listinv(1,i)
        ib=listinv(2,i)
        ic=listinv(3,i)
        id=listinv(4,i)

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

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Call gmax(ibig)
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
        Write(nrite,'(/,1x,a,i0)') 'build_excl_intra deallocation failure, node: ', idnode
        Call error(0)
     End If
  End If

Contains

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

    Use setup_module, Only : mxatdm,mxexcl

    Implicit None

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

End Subroutine build_excl_intra
