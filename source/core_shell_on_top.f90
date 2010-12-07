Subroutine core_shell_on_top()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for positioning shells on top of their cores
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,      Only : idnode,mxnode,gsync,gcheck
  Use setup_module,      Only : nrite,mxshl
  Use config_module,     Only : natms,nlast,lsi,lsa,xxx,yyy,zzz
  Use core_shell_module, Only : ntshl,listshl

  Implicit None

  Logical :: safe
  Integer :: fail,i,j,ia,ib,local_index

  Logical, Allocatable :: lunsafe(:)

  fail=0
  Allocate (lunsafe(1:mxshl), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_on_top allocation failure, node: ', idnode
     Call error(0)
  End If


! Coincide shell with thir cores

  Do i=1,ntshl
     lunsafe(i)=.false.

! indices of atoms in a core-shell

     ia=local_index(listshl(1,i),nlast,lsi,lsa) ! This is a core
     ib=local_index(listshl(2,i),nlast,lsi,lsa) ! This is a shell

! For every shell in the domain get the coordinates of its
! coresponding core (which must be in the domain+hello
! area by construction, if not go to a controlled termination)

     If (ib > 0 .and. ib <= natms .and. ia > 0) Then
        xxx(ib)=xxx(ia)
        yyy(ib)=yyy(ia)
        zzz(ib)=zzz(ia)
     End If

! Detect uncompressed unit

     If ( ((ia > 0 .and. ia <= natms) .or.   &
           (ib > 0 .and. ib <= natms)) .and. &
          (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
  End Do

! Check for uncompressed units

  safe = .not. Any(lunsafe(1:ntshl))
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     Do j=0,mxnode-1
        If (idnode == j) Then
           Do i=1,ntshl
              If (lunsafe(i)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                 '*** warning - global unit number', listshl(0,i), &
                 ' , with a head particle number', listshl(1,i),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(100)
  End If

  Deallocate (lunsafe, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_on_top deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine core_shell_on_top
