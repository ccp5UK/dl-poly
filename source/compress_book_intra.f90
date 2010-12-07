Subroutine compress_book_intra(mx_u,nt_u,b_u,list_u,mxf_u,leg_u)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to prevent bookeeping arrays from expanding
! when execution is on many nodes, mxnode>0, (shells, constraints, PMFs
! and RBs are dealt differntly pass_shared_units and pmf_units_set)
!
! Note: This routine is to be only called from relocate_particles
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : gcheck
  Use setup_module
  Use config_module, Only : natms,lsi,lsa

  Implicit None

  Integer, Intent( In    ) :: mxf_u,mx_u,b_u
  Integer, Intent( InOut ) :: nt_u,list_u(0:b_u,1:mx_u),leg_u(0:mxf_u,1:mxatdm)

  Logical :: ok,keep
  Integer :: i,j,k,l,m,local_index

! is it ok not to do it since it's safe - there's enough buffering space

  ok=.true.
  If (mx_u > 0) Then
     ok=.not.(Real(nt_u,wp)/Real(mx_u,wp) > 0.85_wp)
     Call gcheck(ok)
  End If

  If (.not.ok) Then
     k=0
     Do While (k < nt_u)
        k=k+1
10      Continue

        keep=.false.
        Do i=1,b_u
           keep=keep .or. (local_index(list_u(i,k),natms,lsi,lsa) /= 0)
        End Do

        If (.not.keep) Then
           If      (k  < nt_u) Then
              Do i=1,b_u
                 j=local_index(list_u(i,nt_u),natms,lsi,lsa)
                 If (j > 0) Then         ! For all particles in list_u(1:b_u,nt_u),
                    m=leg_u(0,j)         ! if present on this node, repoint unit
                    Do l=1,m             ! 'nt_u' to 'k' in their leg_u array
                       If (leg_u(l,j) == nt_u) leg_u(l,j) = k
                    End Do
                 End If
              End Do

              list_u(:,k)=list_u(:,nt_u) ! Copy list content from 'nt_u' to 'k'
              list_u(:,nt_u)=0           ! Remove list content in 'nt_u'
              nt_u=nt_u-1                ! Reduce 'nt_u' pointer

              Go To 10                   ! Go back and check it all again for the new list content in 'k'
           Else If (k == nt_u) Then
              list_u(:,nt_u)=0           ! Remove list content in 'k=nt_u'
              nt_u=nt_u-1                ! Reduce 'nt_u' pointer
           End If
        End If
     End Do
  End If

End Subroutine compress_book_intra
