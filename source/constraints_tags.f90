Subroutine constraints_tags(imcon,lstitr,lstopt,dxx,dyy,dzz,listot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for identifying, constructing and indexing
! constraints' vectors for iterative (constraints) algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsync,gcheck
  Use setup_module
  Use config_module,      Only : cell,natms,nlast,lsi,lsa,lfrzn, &
                                 xxx,yyy,zzz
  Use constraints_module, Only : ntcons,listcon


  Implicit None

  Integer,           Intent( In    ) :: imcon
  Logical,           Intent( InOut ) :: lstitr(1:mxatms)
  Integer,           Intent(   Out ) :: lstopt(0:2,1:mxcons)
  Real( Kind = wp ), Intent(   Out ) :: dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons)
  Integer,           Intent(   Out ) :: listot(1:mxatms)

  Logical :: safe
  Integer :: fail,i,j,k,l,local_index

  Logical, Allocatable :: lunsafe(:)

  fail=0
  Allocate (lunsafe(1:mxcons), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_tags allocation failure, node: ', idnode
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
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     Do l=0,mxnode-1
        If (idnode == l) Then
           Do k=1,ntcons
              If (lunsafe(k)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                 '*** warning - global unit number', listcon(0,k), &
                 ' , with a head particle number', listcon(1,k),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
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
     Write(nrite,'(/,1x,a,i0)') 'constraints_tags deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine constraints_tags
