Subroutine rigid_bodies_tags()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for indentifying and indexing RB units
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,  Only : idnode,mxnode,gsync,gcheck
  Use setup_module,  Only : nrite,mxrgd
  Use config_module, Only : natms,nlast,lsi,lsa
  Use rigid_bodies_module

  Implicit None

  Logical :: safe
  Integer :: fail,irgd,jrgd,lrgd,s,i,local_index

  Logical, Allocatable :: lunsafe(:)

  fail=0
  Allocate (lunsafe(1:mxrgd), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_tags allocation failure, node: ', idnode
     Call error(0)
  End If


! Loop over all local to this node RB units and save the indices of the members
! and their presence on the domain in total

  Do irgd=1,ntrgd
     lunsafe(irgd)=.false.
     lrgd=listrgd(-1,irgd)

! Initialise local indices

     indrgd(:,irgd)=0
     Do jrgd=1,lrgd
        s=listrgd(jrgd,irgd)
        i=local_index(s,nlast,lsi,lsa)

        indrgd(jrgd,irgd)=i
        If (i > 0 .and. i <= natms) indrgd(0,irgd)=indrgd(0,irgd)+1
     End Do

! Detect uncompressed unit

     If (Any(indrgd(1:lrgd,irgd) == 0) .and. indrgd(0,irgd) > 0) lunsafe(irgd)=.true.
  End Do

! Check if a RB unit has a diameter > rcut (the system cutoff)
! Check for uncompressed units

  safe = .not. Any(lunsafe(1:ntrgd))
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     Do i=0,mxnode-1
        If (idnode == i) Then
           Do irgd=1,ntrgd
              If (lunsafe(irgd)) Write(nrite,'(/,1x,a,2(i10,a))')        &
                 '*** warning - global unit number', listrgd(0,irgd), &
                 ' , with a head particle number', listrgd(1,irgd),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(642)
  End If

  Deallocate (lunsafe, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_tags deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rigid_bodies_tags
