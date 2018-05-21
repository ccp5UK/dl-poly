Subroutine pass_shared_units &
           (mx_u,b_l,b_u,nt_u,list_u,mxf_u,leg_u,lshmv,lishp,lashp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for passing information about units between nodes
! and compressing list_u if needed
!
! Note (A): 'Use rigid_bodies_module' is only used when RB compression
!           occurs in order to compress q., p., rgdv.. & rgdo.. arrays.
!           RBs are detected by b_l=-1!
!
! Note (B): Called only in:
!                       (1) build_book_intra when mxnode>1 and meg_u>0
!                           but not fully functional - no compressing
!                       (2) relocate_particles if mxnode>1 and meg_u>0
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module,      Only : map,mop
  Use config_module,       Only : natms,nlast,lsi,lsa
  Use rigid_bodies_module, Only : q0,q1,q2,q3,          &
                                  rgdvxx,rgdvyy,rgdvzz, &
                                  rgdoxx,rgdoyy,rgdozz

  Implicit None

  Integer, Intent( In    ) :: mx_u,b_l,b_u,mxf_u
  Integer, Intent( InOut ) :: nt_u,list_u(b_l:b_u,1:mx_u),leg_u(0:mxf_u,1:mxatdm)

  Logical, Intent(   Out ) :: lshmv
  Integer, Intent(   Out ) :: lishp(1:mxlshp),lashp(1:mxproc)

  Logical, Save :: oldjob = .false.

  Logical :: safe,ok
  Integer :: fail,i,j,k,l,m,n_k,n_nt,k0,l_me,l_out,l_in,jdnode,kdnode,local_index

  Integer, Dimension( : ), Allocatable :: i0,j0,listme,lstout,listin

  fail=0
  Allocate (i0(1:b_u),j0(1:b_u),listme(1:mxatms),lstout(1:mxatms),listin(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pass_shared_units allocation failure, node: ', idnode
     Call error(0)
  End If

! Initialise safety flag

  safe=.true.

! is it ok not to compress the bookkeeping list arrays
! since it's safe - there's enough buffering space

  ok=.true.
  If (oldjob) Then
     If (mx_u > 0) Then
        If (b_l == -1) Then
           ok=.false. ! Ghosting is not allowed for RBs
        Else
           ok=.not.(Real(nt_u,wp)/Real(mx_u,wp) > 0.85_wp)
        End If
        Call gcheck(ok)
     End If
  End If

! initialise lshmv and lishp and lashp arrays

  lshmv=.false. ! sharing flag
!  lishp=0       ! list of shared particles (DEBUG)
!  lashp=0       ! break-down of lishp onto a DD map around idnode (DEBUG)

! zero list arrays and last element pointers

  l_me =0 !; listme=0 ! (DEBUG)
  l_out=0 !; lstout=0 ! (DEBUG)

! For all units

  k=0
  Do While (k < nt_u)
     k=k+1
10   Continue

! define members-per-unit limit

     If (b_l == -1) Then
        n_k=list_u(-1,k) ! This can vary but is >= 2 and <= mxlrgd
     Else
        n_k=b_u          ! This is 2 for core-shell and constraint units
     End If

! get domain only local indices of unit members

     i0=0
     Do i=1,n_k
        i0(i)=local_index(list_u(i,k),nlast,lsi,lsa)
        If (i0(i) > natms) i0(i)=0
     End Do

     If (All(i0(1:n_k) == 0) .and. (.not.ok)) Then
20      Continue

! If the whole unit has moved out of this node - compress list_u and leg_u

        If      (k  < nt_u) Then

! define members-per-unit limit

           If (b_l == -1) Then
              n_nt=list_u(-1,nt_u) ! This can vary but is >= 2 and <= mxlrgd
           Else
              n_nt=b_u             ! This is 2 for core-shell and constraint units
           End If

! get domain only local indices of unit members

           j0=0
           Do i=1,n_nt
              j0(i)=local_index(list_u(i,nt_u),nlast,lsi,lsa)
              If (j0(i) > natms) j0(i)=0
           End Do

           If (Any(j0(1:n_nt) > 0)) Then ! Do repointing

              Do i=1,n_nt
                 j=j0(i)
                 If (j > 0) Then      ! For all particles in list_u(1:n_nt,nt_u),
                    m=Abs(leg_u(0,j)) ! if present on this node, repoint unit [for shell particles legshl(0,:)=-1]
                    Do l=1,m          ! 'nt_u' to 'k' in their leg_u array
                       If (leg_u(l,j) == nt_u) leg_u(l,j) = k
                    End Do
                 End If
              End Do

              list_u(:,k)=list_u(:,nt_u) ! Copy list content from 'nt_u' to 'k'

              If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo..
                 q0(k)=q0(nt_u)
                 q1(k)=q1(nt_u)
                 q2(k)=q2(nt_u)
                 q3(k)=q3(nt_u)

                 rgdvxx(k)=rgdvxx(nt_u)
                 rgdvyy(k)=rgdvyy(nt_u)
                 rgdvzz(k)=rgdvzz(nt_u)

                 rgdoxx(k)=rgdoxx(nt_u)
                 rgdoyy(k)=rgdoyy(nt_u)
                 rgdozz(k)=rgdozz(nt_u)
              End If

              list_u(:,nt_u)=0           ! Remove list content in 'nt_u'

              If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo..
                 q0(nt_u)=0.0_wp
                 q1(nt_u)=0.0_wp
                 q2(nt_u)=0.0_wp
                 q3(nt_u)=0.0_wp

                 rgdvxx(nt_u)=0.0_wp
                 rgdvyy(nt_u)=0.0_wp
                 rgdvzz(nt_u)=0.0_wp

                 rgdoxx(nt_u)=0.0_wp
                 rgdoyy(nt_u)=0.0_wp
                 rgdozz(nt_u)=0.0_wp
              End If

              nt_u=nt_u-1                ! Reduce 'nt_u' pointer

           Else

              list_u(:,nt_u)=0           ! Remove list content in 'nt_u'

              If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo..
                 q0(nt_u)=0.0_wp
                 q1(nt_u)=0.0_wp
                 q2(nt_u)=0.0_wp
                 q3(nt_u)=0.0_wp

                 rgdvxx(nt_u)=0.0_wp
                 rgdvyy(nt_u)=0.0_wp
                 rgdvzz(nt_u)=0.0_wp

                 rgdoxx(nt_u)=0.0_wp
                 rgdoyy(nt_u)=0.0_wp
                 rgdozz(nt_u)=0.0_wp
              End If

              nt_u=nt_u-1                ! Reduce 'nt_u' pointer

              Go To 20 ! Go back and check again for the new list contents in 'nt_u'

           End If

           Go To 10    ! Go back and check it all again for the new list contents in 'k'

        Else If (k == nt_u) Then

           list_u(:,nt_u)=0           ! Remove list content in 'k=nt_u'

           If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo.. @ 'k=nt_u'
              q0(nt_u)=0.0_wp
              q1(nt_u)=0.0_wp
              q2(nt_u)=0.0_wp
              q3(nt_u)=0.0_wp

              rgdvxx(nt_u)=0.0_wp
              rgdvyy(nt_u)=0.0_wp
              rgdvzz(nt_u)=0.0_wp

              rgdoxx(nt_u)=0.0_wp
              rgdoyy(nt_u)=0.0_wp
              rgdozz(nt_u)=0.0_wp
           End If

           nt_u=nt_u-1                ! Reduce 'nt_u' pointer

        End If

     Else If (Any(i0(1:n_k) /= 0) .and. Any(i0(1:n_k) == 0)) Then

! The unit is partly shared between domains

        If (l_me > mxatms .or. l_out > mxatms) Then

! Test and avoid possible array bound overflow

           safe=.false.

        Else

! Construct list of bonded atoms crossing domains

           Do i=1,n_k
              If (i0(i) > 0) Then
                 l_me=l_me+1
                 listme(l_me)=list_u(i,k)
              Else
                 l_out=l_out+1
                 lstout(l_out)=list_u(i,k)
              End If
           End Do

        End If

     End If
  End Do

! check for array overflow

  Call gcheck(safe)
  If (.not.safe) Call error(104)

! check for a mess-up

  If (l_me*l_out == 0 .and. Abs(l_me)+Abs(l_out) /= 0) Call error(118)

! sort and discard multiple entries by resetting pointers

  If (l_me*l_out /= 0) Then
     Call shellsort(l_me ,listme)
     Call shellsort(l_out,lstout)

! The following bit is not needed for for core-shell units
! but is present as this subroutine is for shared use

     j=1
     Do i=2,l_me
        If (listme(j) /= listme(i)) Then
           j=j+1
           listme(j)=listme(i)
        End If
     End Do
     l_me=j

     j=1
     Do i=2,l_out
        If (lstout(j) /= lstout(i)) Then
           j=j+1
           lstout(j)=lstout(i)
        End If
     End Do
     l_out=j
  End If

! synchronise all processors

  Call gsync()

! loop over all processors connected to this one (I'm surrounded by 26)

  m=0
  Do k=1,26

! For all genuinely unique processors (no images or myself)

     If (mop(k) == 0) Then

! inter node communication
! jdnode - destination (send to), kdnode - source (receive from)

        kdnode=map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        jdnode=map(k0)

! transmit length of message

        l_in=0
        Call MPI_IRECV(l_in,1,MPI_INTEGER,kdnode,PassUnit_tag+k,dlp_comm_world,request,ierr)
        Call MPI_SEND(l_out,1,MPI_INTEGER,jdnode,PassUnit_tag+k,dlp_comm_world,ierr)
        Call MPI_WAIT(request,status,ierr)

! transmit atom list of units

        listin=0
        If (l_in  > 0) Call MPI_IRECV(listin,l_in,MPI_INTEGER,kdnode,PassUnit_tag+k,dlp_comm_world,request,ierr)
        If (l_out > 0) Call MPI_SEND(lstout,l_out,MPI_INTEGER,jdnode,PassUnit_tag+k,dlp_comm_world,ierr)
        If (l_in  > 0) Call MPI_WAIT(request,status,ierr)

        Do i=1,l_me
           Do j=1,l_in
              If (listme(i) == listin(j)) Then
                 m=m+1
                 If (m > mxlshp) Then
                    safe=.false.
                 Else
                    lishp(m)=listme(i)
                 End If
              End If
           End Do
        End Do

     End If

     lashp(k)=m

  End Do

! check for array overflow

  Call gcheck(safe)
  If (.not.safe) Call error(103)

! Get sharing flag

  lshmv=(m == 0)
  Call gcheck(lshmv)
  lshmv=.not.lshmv

  If (.not.oldjob) oldjob = .true.

  Deallocate (i0,j0,listme,lstout,listin, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pass_shared_units deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pass_shared_units