Module constraints

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for defining global constraint bonds variables and
  ! arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only           : wp
  Use comms,        Only    : comms_type,gsum,gcheck
  Use configuration,   Only : natms,lfrzn,nlast, vxx,vyy,vzz,weight,lsa,lsi, &
    imcon,cell,xxx,yyy,zzz
  Use setup_module, Only    : mxtmls,mxtcon,mxcons,mxfcon,mxlshp,mxproc,mxatdm, &
    mxatms,nrite

  Implicit None

  Logical,                        Save :: lshmv_con = .false.

  Integer,                        Save :: ntcons  = 0 , &
    ntcons1 = 0 , &
    m_con   = 0

  Real( Kind = wp ),              Save :: passcnq(1:5) = (/ & ! QUENCHING per call
    0.0_wp         ,  & ! cycles counter
    0.0_wp         ,  & ! access counter
    0.0_wp         ,  & ! average cycles
    999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
    0.0_wp /)           ! maximum cycles
  Real( Kind = wp ),              Save :: passcon(1:5,1:2,1:2) = Reshape( (/ & ! dim::1-shake, dim:1:-per-call
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::1-shake, dim:2:-per-tst
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:1:-per-call
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:2:-per-tst
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp /) , (/5,2,2/) )


  Integer,           Allocatable, Save :: numcon(:)
  Integer,           Allocatable, Save :: lstcon(:,:),listcon(:,:),legcon(:,:)
  Integer,           Allocatable, Save :: lishp_con(:),lashp_con(:)

  Real( Kind = wp ), Allocatable, Save :: prmcon(:)

  Public :: allocate_constraints_arrays , deallocate_constraints_arrays
  Public :: constraints_pseudo_bonds
  Public :: constraints_quench
  Public :: constraints_tags

Contains

  Subroutine allocate_constraints_arrays()

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numcon(1:mxtmls),                        Stat = fail(1))
    Allocate (lstcon(1:2,1:mxtcon),                    Stat = fail(2))
    Allocate (listcon(0:2,1:mxcons),                   Stat = fail(3))
    Allocate (legcon(0:mxfcon,1:mxatdm),               Stat = fail(4))
    Allocate (lishp_con(1:mxlshp),lashp_con(1:mxproc), Stat = fail(5))
    Allocate (prmcon(1:mxtcon),                        Stat = fail(6))

    If (Any(fail > 0)) Call error(1018)

    numcon  = 0
    lstcon  = 0
    listcon = 0
    legcon  = 0

    lishp_con = 0 ; lashp_con = 0

    prmcon  = 0.0_wp

  End Subroutine allocate_constraints_arrays

  Subroutine deallocate_constraints_arrays()

    Implicit None
    Integer :: fail

    fail = 0

    Deallocate (numcon,lstcon, Stat = fail)

    If (fail > 0) Call error(1032)

  End Subroutine deallocate_constraints_arrays

  Subroutine constraints_pseudo_bonds(lstopt,dxx,dyy,dzz,gxx,gyy,gzz,engcon,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for treating constraint bonds as stiff harmonic
    ! springs for use with the conjugate gradient method (minimise_relax.f90)
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith & i.t.todorov november 2011
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    )  :: lstopt(0:2,1:mxcons)
    Real( Kind = wp ), Intent( In    )  :: dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons)
    Real( Kind = wp ), Intent( InOut )  :: gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms)
    Real( Kind = wp ), Intent(   Out )  :: engcon
    Type( comms_type ), Intent( InOut ) :: comm

    Real( Kind = wp ), Parameter :: rigid=1.0e6_wp

    Integer           :: k,i,j
    Real( Kind = wp ) :: r,r0,ebond,gamma

    engcon=0.0_wp
    Do k=1,ntcons
      If (lstopt(0,k) == 0) Then
        i=lstopt(1,k)
        j=lstopt(2,k)

        ! if a pair is frozen and constraint bonded, it is more frozen
        ! than constrained (users!!!)

        r=Sqrt(dxx(k)**2+dyy(k)**2+dzz(k)**2)
        r0=prmcon(listcon(0,k))

        gamma=rigid*(r-r0)
        ebond=gamma*0.5_wp*(r-r0)
        gamma=gamma/r

        ! Accumulate energy and add forces

        If (i <= natms) Then
          engcon=engcon+ebond

          If (lfrzn(i) == 0) Then
            gxx(i)=gxx(i)-dxx(k)*gamma
            gyy(i)=gyy(i)-dyy(k)*gamma
            gzz(i)=gzz(i)-dzz(k)*gamma
          End If
        End If

        If (j <= natms .and. lfrzn(j) == 0) Then
          gxx(j)=gxx(j)+dxx(k)*gamma
          gyy(j)=gyy(j)+dyy(k)*gamma
          gzz(j)=gzz(j)+dzz(k)*gamma
        End If
      End If
    End Do

    ! global sum of energy

    Call gsum(comm,engcon)

  End Subroutine constraints_pseudo_bonds

  Subroutine constraints_quench(mxshak,tolnce,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for quenching the internal bond energies in the
    ! initial structure of a molecule defined by constraints
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: mxshak
    Real( Kind = wp ), Intent( In    ) :: tolnce
    Type( comms_type), Intent( InOut ) :: comm

    Logical           :: safe
    Integer           :: fail(1:4),i,j,k,icyc
    Real( Kind = wp ) :: dis,amti,amtj,dlj,dli,esig,gamma,gammi,gammj

    Logical,           Allocatable :: lstitr(:)
    Integer,           Allocatable :: lstopt(:,:),listot(:)
    Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)

    fail=0
    Allocate (lstitr(1:mxatms),                          Stat=fail(1))
    Allocate (lstopt(0:2,1:mxcons),listot(1:mxatms),     Stat=fail(2))
    Allocate (dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons), Stat=fail(3))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms), Stat=fail(4))
    If (Any(fail > 0)) Then
      Write(nrite,'(/,1x,a,i0)') 'constraints_quench allocation failure, node: ', comm%idnode
      Call error(0)
    End If

    ! gather velocities of shared atoms

    If (lshmv_con) Call update_shared_units(natms,nlast,lsi,lsa,lishp_con,lashp_con,vxx,vyy,vzz)

    ! construct current constrained bond vectors and listot array (shared
    ! constraint atoms) for iterative (constraints) algorithms

    lstitr(1:natms)=.false. ! initialise lstitr
    Call constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot,comm)

    ! normalise constraint vectors

    Do k=1,ntcons
      If (lstopt(0,k) == 0) Then
        dis=1.0_wp/Sqrt(dxx(k)**2+dyy(k)**2+dzz(k)**2)
        dxx(k)=dxx(k)*dis
        dyy(k)=dyy(k)*dis
        dzz(k)=dzz(k)*dis
      Else ! DEBUG
        !        dxx(k)=0.0_wp
        !        dyy(k)=0.0_wp
        !        dzz(k)=0.0_wp
      End If
    End Do

    ! application of constraint (quench) algorithm
    ! Initialise number of cycles to zero and unsafe passage of the algorithm

    icyc=0
    safe=.false.

    Do While ((.not.safe) .and. icyc < mxshak)
      icyc=icyc+1

      ! initialise velocity correction arrays

      Do i=1,natms
        vxt(i)=0.0_wp
        vyt(i)=0.0_wp
        vzt(i)=0.0_wp
      End Do

      ! calculate velocity constraint corrections

      esig=0.0_wp
      Do k=1,ntcons
        If (lstopt(0,k) == 0) Then
          i=lstopt(1,k)
          j=lstopt(2,k)

          ! if a pair is frozen and constraint bonded, it is more frozen
          ! than constrained (users!!!)

          amti=1.0_wp/weight(i)
          amtj=1.0_wp/weight(j)

          ! no corrections for frozen atoms

          If (lfrzn(i) /= 0) amti=0.0_wp
          If (lfrzn(j) /= 0) amtj=0.0_wp

          ! calculate constraint force parameter - gamma

          gamma = dxx(k)*(vxx(i)-vxx(j)) + dyy(k)*(vyy(i)-vyy(j)) + dzz(k)*(vzz(i)-vzz(j))

          esig=Max(esig,0.5_wp*Abs(gamma))

          gamma = gamma / (amti+amtj)

          ! improve approximate constraint velocity

          If (i <= natms .and. lfrzn(i) == 0) Then
            gammi =-gamma*amti
            vxt(i)=vxt(i)+dxx(k)*gammi
            vyt(i)=vyt(i)+dyy(k)*gammi
            vzt(i)=vzt(i)+dzz(k)*gammi
          End If

          If (j <= natms .and. lfrzn(j) == 0) Then
            gammj = gamma*amtj
            vxt(j)=vxt(j)+dxx(k)*gammj
            vyt(j)=vyt(j)+dyy(k)*gammj
            vzt(j)=vzt(j)+dzz(k)*gammj
          End If
        End If
      End Do

      ! global verification of convergence

      safe=(esig < tolnce)
      Call gcheck(comm,safe,"enforce")

      ! bypass next section and terminate iteration if all tolerances ok

      If (.not.safe) Then

        ! update velocities

        Do k=1,ntcons
          If (lstopt(0,k) == 0) Then
            i=lstopt(1,k)
            j=lstopt(2,k)

            If (i <= natms .and. lfrzn(i) == 0) Then
              dli = 1.0_wp/Real(listot(i),wp)
              vxx(i)=vxx(i)+vxt(i)*dli
              vyy(i)=vyy(i)+vyt(i)*dli
              vzz(i)=vzz(i)+vzt(i)*dli
            End If

            If (j <= natms .and. lfrzn(j) == 0) Then
              dlj = 1.0_wp/Real(listot(j),wp)
              vxx(j)=vxx(j)+vxt(j)*dlj
              vyy(j)=vyy(j)+vyt(j)*dlj
              vzz(j)=vzz(j)+vzt(j)*dlj
            End If
          End If
        End Do

        ! transport velocity updates to other nodes

        If (lshmv_con) Call update_shared_units(natms,nlast,lsi,lsa,lishp_con,lashp_con,vxx,vyy,vzz)

      End If
    End Do

    If (.not.safe) Then ! error exit if quenching fails
      Call error(70)
    Else ! Collect per call passage statistics
      passcnq(1)=Real(icyc-1,wp)
      passcnq(3)=passcnq(2)*passcnq(3)
      passcnq(2)=passcnq(2)+1.0_wp
      passcnq(3)=passcnq(3)/passcnq(2)+passcnq(1)/passcnq(2)
      passcnq(4)=Min(passcnq(1),passcnq(4))
      passcnq(5)=Max(passcnq(1),passcnq(5))
      passcnq(1)=0.0_wp ! Reset
    End If

    Deallocate (lstitr,        Stat=fail(1))
    Deallocate (lstopt,listot, Stat=fail(2))
    Deallocate (dxx,dyy,dzz,   Stat=fail(3))
    Deallocate (vxt,vyt,vzt,   Stat=fail(4))
    If (Any(fail > 0)) Then
      Write(nrite,'(/,1x,a,i0)') 'constraints_quench deallocation failure, node: ', comm%idnode
      Call error(0)
    End If

  End Subroutine constraints_quench

  Subroutine constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for identifying, constructing and indexing
    ! constraints' vectors for iterative (constraints) algorithms
    !
    ! Note: must be used in conjunction with integration algorithms
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( InOut ) :: lstitr(1:mxatms)
    Integer,            Intent(   Out ) :: lstopt(0:2,1:mxcons)
    Real( Kind = wp ),  Intent(   Out ) :: dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons)
    Integer,            Intent(   Out ) :: listot(1:mxatms)
    Type( comms_type ), Intent( InOut ) :: comm

    Logical :: safe
    Integer :: fail,i,j,k,l,local_index

    Logical, Allocatable :: lunsafe(:)

    fail=0
    Allocate (lunsafe(1:mxcons), Stat=fail)
    If (fail > 0) Then
      Write(nrite,'(/,1x,a,i0)') 'constraints_tags allocation failure, node: ', comm%idnode
      Call error(0)
    End If

    ! initialise listot array (shared constraint bond)

    Do k=1,natms
      listot(k)=0
    End Do

    Do k=1,ntcons
      lunsafe(k)=.false.

      ! Opt out bond from action by default

      lstopt(0,k)=1

      ! indices of atoms in bond

      i=local_index(listcon(1,k),nlast,lsi,lsa)
      j=local_index(listcon(2,k),nlast,lsi,lsa)

      ! store indices

      lstopt(1,k)=i
      lstopt(2,k)=j

      ! for all native and natively shared constraints

      If ((i > 0 .and. j > 0) .and. (i <= natms .or. j <= natms)) Then

        ! if a pair is frozen and constraint bonded
        ! it is more frozen than constrained (users!!!)

        If (lfrzn(i)*lfrzn(j) == 0) Then

          ! Select bond for action if product gives zero

          lstopt(0,k)=0

          ! calculate bond vectors

          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)

          ! indicate sharing on local ends of bonds
          ! summed contributions (quench/shake/rattle) for each local
          ! constrained bonded atom must be weighted by listot(atom) -
          ! how many bond constraints an atom is involved in

          If (i <= natms) listot(i)=listot(i)+1
          If (j <= natms) listot(j)=listot(j)+1

        End If

      Else If ((i == 0 .and. j == 0) .or. (i > natms .and. j > natms) .or. &
        (i == 0 .and. j > natms) .or. (j == 0 .and. i > natms)) Then

        ! constraints lying outside or completely in the halo, or partly in the
        ! halo and partly outside it are not considered and zero bond vectors
        ! are assigned (DEBUG)

        !        dxx(k)=0.0_wp
        !        dyy(k)=0.0_wp
        !        dzz(k)=0.0_wp

      Else

        ! weird case - constraints flying apart
        ! Detect uncompressed unit

        lunsafe(k)=.true.

      End If
    End Do

    ! check for incomplete bond data

    safe = .not. Any(lunsafe(1:ntcons))
    Call gcheck(comm,safe)
    If (.not.safe) Then
      Do l=0,comm%mxnode-1
        If (comm%idnode == l) Then
          Do k=1,ntcons
            If (lunsafe(k)) Write(nrite,'(/,1x,a,2(i10,a))')     &
              '*** warning - global unit number', listcon(0,k), &
              ' , with a head particle number', listcon(1,k),   &
              ' contributes towards next error !!! ***'
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(111)
    End If

    ! minimum image convention for bond vectors

    Call images(imcon,cell,ntcons,dxx,dyy,dzz)

    ! update lstitr

    Do k=1,natms
      lstitr(k)=(listot(k) > 0 .and. lfrzn(k) == 0)
    End Do

    Deallocate (lunsafe, Stat=fail)
    If (fail > 0) Then
      Write(nrite,'(/,1x,a,i0)') 'constraints_tags deallocation failure, node: ', comm%idnode
      Call error(0)
    End If

  End Subroutine constraints_tags


End module constraints
