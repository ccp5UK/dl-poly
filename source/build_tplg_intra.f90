Subroutine build_tplg_intra()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing the topology related list
! of neighbours for the MD system mapped onto this node
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES

  Use kinds_f90
  Use comms_module,  Only : mxnode,gcheck,gmax
  Use setup_module

! CONFIG MODULE

  Use config_module, Only : natms,nlast,lsi,lsa

! INTERACTION MODULES

  Use bonds_module
  Use angles_module
  Use dihedrals_module
  Use inversions_module

! MULTIPOLES MODULE

  Use mpoles_module, Only : ltpatm ! equivalent to lexatm in config_module

  Implicit None

  Logical :: safe
  Integer :: i,j,ia,ib,ic,id,ia0,ib0,ic0,id0,ibig,local_index

! variables for array bound checking

  ibig=0
  safe=.true.

! exclude sites on basis of chemical bonds

  Do i=1,ntbond
     If (Abs(keybnd(listbnd(0,i))) > 0) Then
        ia=listbnd(1,i)
        ib=listbnd(2,i)

        ia0=local_index(ia,nlast,lsi,lsa)
        ib0=local_index(ib,nlast,lsi,lsa)

        If (ia0 > natms) ia0=0
        If (ib0 > natms) ib0=0

! add atoms to exclusion list

        If (ia0 > 0) Call add_exclusion(safe,ia0,ib,ibig,ltpatm)
        If (ib0 > 0) Call add_exclusion(safe,ib0,ia,ibig,ltpatm)
     End If
  End Do

! exclude sites on basis of bond angles

  Do i=1,ntangl
     If (Abs(keyang(listang(0,i))) > 0) Then
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

        If (ia0 > 0) Then ! ia : ib interactions
           Call add_exclusion(safe,ia0,ib,ibig,ltpatm)
        End If

        If (ib0 > 0) Then ! ib : ia - ic interactions
           Call add_exclusion(safe,ib0,ia,ibig,ltpatm)
           Call add_exclusion(safe,ib0,ic,ibig,ltpatm)
        End If

        If (ic0 > 0) Then ! ic : ib interactions
           Call add_exclusion(safe,ic0,ib,ibig,ltpatm)
        End If
     End If
  End Do

! exclude sites on basis of dihedral angles

  Do i=1,ntdihd
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

     If (ia0 > 0) Then ! ia : ib interactions
        Call add_exclusion(safe,ia0,ib,ibig,ltpatm)
     End If

     If (ib0 > 0) Then ! ib : ia - ic interactions
        Call add_exclusion(safe,ib0,ia,ibig,ltpatm)
        Call add_exclusion(safe,ib0,ic,ibig,ltpatm)
     End If

     If (ic0 > 0) Then ! ic : ib - id interactions
        Call add_exclusion(safe,ic0,ib,ibig,ltpatm)
        Call add_exclusion(safe,ic0,id,ibig,ltpatm)
     End If

     If (id0 > 0) Then ! id : ic interactions
        Call add_exclusion(safe,id0,ic,ibig,ltpatm)
     End If
  End Do

! exclude sites on basis of inversion potentials

  Do i=1,ntinv
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
        Call add_exclusion(safe,ia0,ib,ibig,ltpatm)
        Call add_exclusion(safe,ia0,ic,ibig,ltpatm)
        Call add_exclusion(safe,ia0,id,ibig,ltpatm)
     End If

     If (ib0 > 0) Then ! ib : ia interactions
        Call add_exclusion(safe,ib0,ia,ibig,ltpatm)
     End If

     If (ic0 > 0) Then ! ic : ia interactions
        Call add_exclusion(safe,ic0,ia,ibig,ltpatm)
     End If

     If (id0 > 0) Then ! id : ia interactions
        Call add_exclusion(safe,id0,ia,ibig,ltpatm)
     End If
  End Do

! check for exceeded array bounds

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Call gmax(ibig)
     Call warning(250,Real(ibig,wp),Real(mxexcl,wp),0.0_wp)
     Call error(65)
  End If

! sort ltpatm

  Do i=1,natms
     j=ltpatm(0,i)
     If (j > 0) Call shellsort(j,ltpatm(1:j,i))
  End Do

Contains

  Subroutine add_exclusion(safe,ia0,ib,ibig,ltpatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to add excluded atoms to the excluded atom list
! provided they are not already excluded
!
! copyright - daresbury laboratory
! author    - w.smith march 1999
! amended   - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use setup_module, Only : mxatdm,mxexcl

    Implicit None

    Integer,                                  Intent( In    ) :: ia0,ib
    Logical,                                  Intent( InOut ) :: safe
    Integer,                                  Intent( InOut ) :: ibig
    Integer, Dimension( 0:mxexcl, 1:mxatdm ), Intent( InOut ) :: ltpatm

    Logical :: safe_local,l_excluded
    Integer :: last

! Get current length

    last = ltpatm(0,ia0)

! Determine possible exclusion

    l_excluded = Any(ltpatm(1:last,ia0) == ib)

    If (.not.l_excluded) Then

! Get local safety no array overfloat

       safe_local = (last < mxexcl-1)

! Determine global safety

       safe = safe .and. safe_local

       If (safe_local) Then

! Increase length of the ia0 exclusion list and tag ib in it

          last = last + 1
          ltpatm(0,ia0) = last
          ltpatm(last,ia0) = ib

          ibig=Max(ibig,last)

       Else

! Collect number of offences

          ltpatm(mxexcl,ia0) = ltpatm(mxexcl,ia0) + 1

          ibig=Max(ibig,last+ltpatm(mxexcl,ia0))

       End If

    End If

  End Subroutine add_exclusion

End Subroutine build_tplg_intra
