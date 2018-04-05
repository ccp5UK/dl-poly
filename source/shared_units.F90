Module shared_units

  Use kinds, Only : wp
  Use comms, Only : comms_type,PassUnit_tag,wp_mpi,UpdShUnit_tag,&
                    gcheck,gsync,gsend
  Use setup
  Use domains,      Only : map,mop
  Use configuration,       Only : natms,nlast,lsi,lsa
  Use numerics, Only : local_index,shellsort
  Use errors_warnings, Only : error
  
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None


  Private
  Public :: update_shared_units
  Public :: update_shared_units_int
  Public :: pass_shared_units
  Public :: tag_legend

  Contains
  
  Subroutine pass_shared_units &
           (mx_u,b_l,b_u,nt_u,list_u,mxf_u,leg_u,lshmv,lishp,lashp,comm,&
           q0,q1,q2,q3,rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for passing information about units between nodes
! and compressing list_u if needed
!
! Note (A): 'Use rigid_bodies' is only used when RB compression
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



  Integer, Intent( In    ) :: mx_u,b_l,b_u,mxf_u
  Integer, Intent( InOut ) :: nt_u,list_u(b_l:b_u,1:mx_u),leg_u(0:mxf_u,1:mxatdm)

  Logical, Intent(   Out ) :: lshmv
  Integer, Intent(   Out ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Type( comms_type ), Intent( InOut ) :: comm
  Real( Kind = wp ), Intent( InOut ),Dimension(*) :: q0,q1,q2,q3,rgdvxx,rgdvyy,&
                                                      rgdvzz,rgdoxx,rgdoyy,rgdozz

  Logical, Save :: oldjob = .false.

  Logical :: safe,ok
  Integer :: fail,i,j,k,l,m,n_k,n_nt,k0,l_me,l_out,l_in,jdnode,kdnode

  Integer, Dimension( : ), Allocatable :: i0,j0,listme,lstout,listin
  Character( Len = 256 ) :: message

  fail=0
  Allocate (i0(1:b_u),j0(1:b_u),listme(1:mxatms),lstout(1:mxatms),listin(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pass_shared_units allocation failure'
     Call error(0,message)
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
        Call gcheck(comm,ok)
     End If
  End If

! initialise lshmv and lishp and lashp arrays

  lshmv=.false. ! sharing flag
!  lishp=0       ! list of shared particles (DEBUG)
!  lashp=0       ! break-down of lishp onto a DD map around comm%idnode (DEBUG)

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

  Call gcheck(comm,safe)
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

  Call gsync(comm)

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
        Call MPI_IRECV(l_in,1,MPI_INTEGER,kdnode,PassUnit_tag+k,comm%comm,comm%request,comm%ierr)
        Call gsend(comm,l_out,jdnode,PassUnit_tag+k)
        Call MPI_WAIT(comm%request,comm%status,comm%ierr)

! transmit atom list of units

        listin=0
        If (l_in  > 0) Call MPI_IRECV(listin,l_in,MPI_INTEGER,kdnode,PassUnit_tag+k,comm%comm,comm%request,comm%ierr)
        If (l_out > 0) Then
          Call gsend(comm,lstout,jdnode,PassUnit_tag+k)
        End If
        If (l_in  > 0) Call MPI_WAIT(comm%request,comm%status,comm%ierr)

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

  Call gcheck(comm,safe)
  If (.not.safe) Call error(103)

! Get sharing flag

  lshmv=(m == 0)
  Call gcheck(comm,lshmv)
  lshmv=.not.lshmv

  If (.not.oldjob) oldjob = .true.

  Deallocate (i0,j0,listme,lstout,listin, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pass_shared_units deallocation failure'
     Call error(0,message)
  End If

End Subroutine pass_shared_units

Subroutine update_shared_units(natms,nlast,lsi,lsa,lishp,lashp,qxx,qyy,qzz,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for passing atomic coordinates, velocities or
! just updates of shared core-shell, constraint and RB units
! between nodes
!
! Note: This subroutine is only called when there is sharing.
!       In case of mxnode=1 there is no sharing by construction.
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: natms,nlast
  Integer,           Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer,           Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Real( Kind = wp ), Intent( InOut ) :: qxx(1:mxatms),qyy(1:mxatms),qzz(1:mxatms)
  Type( comms_type), Intent( InOut ) :: comm

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message
! Number of transported quantities per particle

  iadd=4

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units allocation failure'
     Call error(0,message)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync(comm)

! transfer coordinate update data to all neighbouring nodes (26 surrounding me)

  j0=0
  Do k=1,26

! For all genuinely unique processors (no images or myself)

     If (mop(k) == 0) Then

! Initialise the number of buffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),nlast,lsi,lsa)
              If (m > natms) m=0

! m should always be > 0 (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > 0) Then
                 buffer(i+1)=Real(lishp(j),wp)

                 buffer(i+2)=qxx(m)
                 buffer(i+3)=qyy(m)
                 buffer(i+4)=qzz(m)
              Else
                 safe(2)=.false.
              End If
              i=i+iadd
           Else
              safe(1)=.false.
           End If
        End Do

! inter node communication - right the opposite of that in pass_shared_units
! jdnode - destination (send to), kdnode - source (receive from)

        jdnode=map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=map(k0)

! pass only non-zero length messages

! transmit length of message

        n=0
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,UpdShUnit_tag+k,comm%comm,comm%request,comm%ierr)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call MPI_WAIT(comm%request,comm%status,comm%ierr)

        If (n > 0) Call MPI_IRECV(buffer(i+1),n,wp_mpi,kdnode,UpdShUnit_tag+k,comm%comm,comm%request,comm%ierr)
        If (i > 0) Then
          Call gsend(comm,buffer(1),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call MPI_WAIT(comm%request,comm%status,comm%ierr)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from buffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(Nint(buffer(i+1)),nlast,lsi,lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > natms) Then
                 qxx(m)=buffer(i+2)
                 qyy(m)=buffer(i+3)
                 qzz(m)=buffer(i+4)
              Else
                 safe(2)=.false.
              End If
              i=i+iadd
           Else
              safe(1)=.false.
           End If
        End Do

     End If
     j0=lashp(k)

  End Do

  Call gcheck(comm,safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units deallocation failure'
     Call error(0,message)
  End If

End Subroutine update_shared_units

Subroutine update_shared_units_int(natms,nlast,lsi,lsa,lishp,lashp,iii,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for passing integer atomic property (int) updates
! of shared core-shell, constraint and RB units between nodes
!
! Note: This subroutine is only called when there is sharing.
!       In case of mxnode=1 there is no sharing by construction.
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer, Intent( In    ) :: natms,nlast
  Integer, Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer, Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Integer, Intent( InOut ) :: iii(1:mxatms)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Integer, Dimension( : ), Allocatable :: ibuffer
  Character ( Len = 256 ) :: message
  
! Number of transported quantities per particle

  iadd=2

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (ibuffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units_int allocation failure'
     Call error(0,message)
  End If

! Set ibuffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync(comm)

! transfer coordinate update data to all neighbouring nodes (26 surrounding me)

  j0=0
  Do k=1,26

! For all genuinely unique processors (no images or myself)

     If (mop(k) == 0) Then

! Initialise the number of ibuffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),nlast,lsi,lsa)
              If (m > natms) m=0

! m should always be > 0 (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > 0) Then
                 ibuffer(i+1)=lishp(j)

                 ibuffer(i+2)=iii(m)
              Else
                 safe(2)=.false.
              End If
              i=i+iadd
           Else
              safe(1)=.false.
           End If
        End Do

! inter node communication - right the opposite of that in pass_shared_units
! jdnode - destination (send to), kdnode - source (receive from)

        jdnode=map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=map(k0)

! pass only non-zero length messages

! transmit length of message

        n=0
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,UpdShUnit_tag+k,comm%comm,comm%request,comm%ierr)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call MPI_WAIT(comm%request,comm%status,comm%ierr)

        If (n > 0) Call MPI_IRECV(ibuffer(i+1),n,MPI_INTEGER,kdnode,UpdShUnit_tag+k,comm%comm,comm%request,comm%ierr)
        If (i > 0) Then
          Call gsend(comm,ibuffer(1),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call MPI_WAIT(comm%request,comm%status,comm%ierr)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from ibuffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(ibuffer(i+1),nlast,lsi,lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > natms) Then
                 iii(m)=ibuffer(i+2)
              Else
                 safe(2)=.false.
              End If
              i=i+iadd
           Else
              safe(1)=.false.
           End If
        End Do

     End If
     j0=lashp(k)

  End Do

  Call gcheck(comm,safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  Deallocate (ibuffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units_int deallocation failure'
     Call error(0,message)
  End If

End Subroutine update_shared_units_int

Subroutine update_shared_units_rwp(natms,nlast,lsi,lsa,lishp,lashp,rrr,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for passing real atomic property (rwp) updates
! of shared core-shell, constraint and RB units between nodes
!
! Note: This subroutine is only called when there is sharing.
!       In case of mxnode=1 there is no sharing by construction.
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: natms,nlast
  Integer,           Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer,           Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Real( Kind = wp ), Intent( InOut ) :: rrr(1:mxatms)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character ( Len = 256 )  :: message

! Number of transported quantities per particle

  iadd=2

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units_rwp allocation failure'
     Call error(0,message)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync(comm)

! transfer coordinate update data to all neighbouring nodes (26 surrounding me)

  j0=0
  Do k=1,26

! For all genuinely unique processors (no images or myself)

     If (mop(k) == 0) Then

! Initialise the number of buffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),nlast,lsi,lsa)
              If (m > natms) m=0

! m should always be > 0 (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > 0) Then
                 buffer(i+1)=Real(lishp(j),wp)

                 buffer(i+2)=rrr(m)
              Else
                 safe(2)=.false.
              End If
              i=i+iadd
           Else
              safe(1)=.false.
           End If
        End Do

! inter node communication - right the opposite of that in pass_shared_units
! jdnode - destination (send to), kdnode - source (receive from)

        jdnode=map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=map(k0)

! pass only non-zero length messages

! transmit length of message

        n=0
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,UpdShUnit_tag+k,comm%comm,comm%request,comm%ierr)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call MPI_WAIT(comm%request,comm%status,comm%ierr)

        If (n > 0) Call MPI_IRECV(buffer(i+1),n,wp_mpi,kdnode,UpdShUnit_tag+k,comm%comm,comm%request,comm%ierr)
        If (i > 0) Then
          Call gsend(comm,buffer(1),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call MPI_WAIT(comm%request,comm%status,comm%ierr)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from buffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(Nint(buffer(i+1)),nlast,lsi,lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > natms) Then
                 rrr(m)=buffer(i+2)
              Else
                 safe(2)=.false.
              End If
              i=i+iadd
           Else
              safe(1)=.false.
           End If
        End Do

     End If
     j0=lashp(k)

  End Do

  Call gcheck(comm,safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units_rwp deallocation failure'
     Call error(0,message)
  End If

End Subroutine update_shared_units_rwp

Subroutine tag_legend(safe,iatm,nt,legend,mxf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to tag the legend arrays recording the
! intra-like descriptions for each atom in a domain
!
! nt should always be supplied positive except for shells!
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,                              Intent( InOut) :: safe
  Integer,                              Intent( In   ) :: iatm,nt,mxf
  Integer, Dimension( 0:mxf, 1:mxatdm), Intent( InOut) :: legend

  Logical :: safe_local
  Integer :: last

! Get current length

  last = Abs(legend(0,iatm))

! Get local safety no array overflow

  safe_local = (last < mxf-1)

! Determine global safety

  safe = safe .and. safe_local

  If (safe_local) Then

! Increase length of links and tag: I, iatm, am linked to one more entity
! (particle) in a unit of this legend type with a local domain number 'nt'

     last = last + 1
     If (.not.(nt < 0 .and. last == 1)) Then ! THE NORMAL CASE
        legend(0,iatm) = last
        legend(last,iatm) = nt
     Else ! This should be hit just the once for shell particles only!!!
        legend(0,iatm) = -1
        legend(last,iatm) = -nt
     End If

  Else

! Collect number of offences

     legend(mxf,iatm) = legend(mxf,iatm) + 1

  End If

End Subroutine tag_legend

End Module shared_units
