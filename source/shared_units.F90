Module shared_units

  Use kinds, Only : wp
  Use particle, Only: corePart
  Use comms, Only : comms_type,PassUnit_tag,wp_mpi,UpdShUnit_tag,&
                    gcheck,gsync,gsend,gwait,girecv
  Use domains,      Only : domains_type
  Use configuration,       Only : configuration_type
  Use numerics, Only : local_index,shellsort
  Use errors_warnings, Only : error

  Implicit None

  Integer, Parameter :: SHARED_UNIT_UPDATE_POSITIONS = 1
  Integer, Parameter :: SHARED_UNIT_UPDATE_FORCES    = 2
  Private
  Public :: update_shared_units
  Public :: update_shared_units_int
  Public :: pass_shared_units
  Public :: tag_legend
  Public :: SHARED_UNIT_UPDATE_POSITIONS
  Public :: SHARED_UNIT_UPDATE_FORCES


  Interface update_shared_units
    Module Procedure update_shared_units_arrays
    Module Procedure update_shared_units_parts
  End Interface update_shared_units
  Contains

    Subroutine pass_shared_units(config,mx_u,b_l,b_u,nt_u,list_u,mxf_u,leg_u,lshmv, &
        lishp,lashp,oldjob,domain,comm,q0,q1,q2,q3,vxx,vyy,vzz,oxx,oyy,ozz)

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


  Type( configuration_type ), Intent( InOut ) :: config
  Integer, Intent( In    ) :: mx_u,b_l,b_u,mxf_u
  Integer, Intent( InOut ) :: nt_u,list_u(b_l:b_u,1:mx_u),leg_u(0:mxf_u,1:config%mxatdm)
  Logical, Intent(   Out ) :: lshmv
  Type( domains_type ), Intent( In    ) :: domain
  Integer, Intent(   Out ) :: lishp(:),lashp(1:domain%neighbours)
  Type( comms_type ), Intent( InOut ) :: comm
  Real( Kind = wp ), Intent( InOut ),Dimension(:) :: q0,q1,q2,q3,vxx,vyy,&
    vzz,oxx,oyy,ozz

  Logical, Intent( InOut ) :: oldjob

  Logical :: safe,ok
  Integer :: fail,i,j,k,l,m,n_k,n_nt,k0,l_me,l_out,l_in,jdnode,kdnode

  Integer, Dimension( : ), Allocatable :: i0,j0,listme,lstout,listin
  Character( Len = 256 ) :: message

  fail=0
  Allocate (i0(1:b_u),j0(1:b_u),listme(1:config%mxatms),lstout(1:config%mxatms),listin(1:config%mxatms), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pass_shared_units allocation failure'
     Call error(0,message)
  End If

! Initialise safety flag

  safe=.true.

! is it ok not to compress the bookkeeping neigh%list arrays
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

! zero neigh%list arrays and last element pointers

  l_me =0 !; listme=0 ! (DEBUG)
  l_out=0 !; lstout=0 ! (DEBUG)

! For all units

  k=0
  Do While (k < nt_u)
     k=k+1
10   Continue

! define members-per-unit limit

     If (b_l == -1) Then
        n_k=list_u(-1,k) ! This can vary but is >= 2 and <= rigid%max_list
     Else
        n_k=b_u          ! This is 2 for core-shell and constraint units
     End If

! get domain only local indices of unit members

     i0=0
     Do i=1,n_k
        i0(i)=local_index(list_u(i,k),config%nlast,config%lsi,config%lsa)
        If (i0(i) > config%natms) i0(i)=0
     End Do

     If (All(i0(1:n_k) == 0) .and. (.not.ok)) Then
20      Continue

! If the whole unit has moved out of this node - compress list_u and leg_u

        If      (k  < nt_u) Then

! define members-per-unit limit

           If (b_l == -1) Then
              n_nt=list_u(-1,nt_u) ! This can vary but is >= 2 and <= rigid%max_list
           Else
              n_nt=b_u             ! This is 2 for core-shell and constraint units
           End If

! get domain only local indices of unit members

           j0=0
           Do i=1,n_nt
              j0(i)=local_index(list_u(i,nt_u),config%nlast,config%lsi,config%lsa)
              If (j0(i) > config%natms) j0(i)=0
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

              list_u(:,k)=list_u(:,nt_u) ! Copy neigh%list content from 'nt_u' to 'k'

              If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo..
                 q0(k)=q0(nt_u)
                 q1(k)=q1(nt_u)
                 q2(k)=q2(nt_u)
                 q3(k)=q3(nt_u)

                 vxx(k)=vxx(nt_u)
                 vyy(k)=vyy(nt_u)
                 vzz(k)=vzz(nt_u)

                 oxx(k)=oxx(nt_u)
                 oyy(k)=oyy(nt_u)
                 ozz(k)=ozz(nt_u)
              End If

              list_u(:,nt_u)=0           ! Remove neigh%list content in 'nt_u'

              If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo..
                 q0(nt_u)=0.0_wp
                 q1(nt_u)=0.0_wp
                 q2(nt_u)=0.0_wp
                 q3(nt_u)=0.0_wp

                 vxx(nt_u)=0.0_wp
                 vyy(nt_u)=0.0_wp
                 vzz(nt_u)=0.0_wp

                 oxx(nt_u)=0.0_wp
                 oyy(nt_u)=0.0_wp
                 ozz(nt_u)=0.0_wp
              End If

              nt_u=nt_u-1                ! Reduce 'nt_u' pointer

           Else

              list_u(:,nt_u)=0           ! Remove neigh%list content in 'nt_u'

              If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo..
                 q0(nt_u)=0.0_wp
                 q1(nt_u)=0.0_wp
                 q2(nt_u)=0.0_wp
                 q3(nt_u)=0.0_wp

                 vxx(nt_u)=0.0_wp
                 vyy(nt_u)=0.0_wp
                 vzz(nt_u)=0.0_wp

                 oxx(nt_u)=0.0_wp
                 oyy(nt_u)=0.0_wp
                 ozz(nt_u)=0.0_wp
              End If

              nt_u=nt_u-1                ! Reduce 'nt_u' pointer

              Go To 20 ! Go back and check again for the new neigh%list contents in 'nt_u'

           End If

           Go To 10    ! Go back and check it all again for the new neigh%list contents in 'k'

        Else If (k == nt_u) Then

           list_u(:,nt_u)=0           ! Remove neigh%list content in 'k=nt_u'

           If (b_l == -1) Then        ! RB handling of q., rgdv.. & rgdo.. @ 'k=nt_u'
              q0(nt_u)=0.0_wp
              q1(nt_u)=0.0_wp
              q2(nt_u)=0.0_wp
              q3(nt_u)=0.0_wp

              vxx(nt_u)=0.0_wp
              vyy(nt_u)=0.0_wp
              vzz(nt_u)=0.0_wp

              oxx(nt_u)=0.0_wp
              oyy(nt_u)=0.0_wp
              ozz(nt_u)=0.0_wp
           End If

           nt_u=nt_u-1                ! Reduce 'nt_u' pointer

        End If

     Else If (Any(i0(1:n_k) /= 0) .and. Any(i0(1:n_k) == 0)) Then

! The unit is partly shared between domains

        If (l_me > config%mxatms .or. l_out > config%mxatms) Then

! Test and avoid possible array bound overflow

           safe=.false.

        Else

! Construct neigh%list of bonded atoms crossing domains

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

     If (domain%map_unique(k) == 0) Then

! inter node communication
! jdnode - destination (send to), kdnode - source (receive from)

        kdnode=domain%map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        jdnode=domain%map(k0)

! transmit length of message

        l_in=0
        Call girecv(comm,l_in,kdnode,PassUnit_tag+k)
        Call gsend(comm,l_out,jdnode,PassUnit_tag+k)
        Call gwait(comm)

! transmit atom neigh%list of units

        listin=0
        If (l_in  > 0) Then
          Call girecv(comm,listin(1:l_in),kdnode,PassUnit_tag+k)
        End If
        If (l_out > 0) Then
          Call gsend(comm,lstout(1:l_out),jdnode,PassUnit_tag+k)
        End If
        If (l_in  > 0) Call gwait(comm)

        Do i=1,l_me
           Do j=1,l_in
              If (listme(i) == listin(j)) Then
                 m=m+1
                 If (m > config%mxlshp) Then
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

Subroutine update_shared_units_parts(config,lishp,lashp,subtype,domain,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for passing atomic coordinates, velocities or
! just updates of shared core-shell, constraint and RB units
! between nodes when the values are stored inside the particle struct
!
! Note: This subroutine is only called when there is sharing.
!       In case of mxnode=1 there is no sharing by construction.
!
! copyright - daresbury laboratory
! author    - a.b.g.chalk july 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( configuration_type), Intent( InOut ) :: config
  Type( domains_type ), Intent( In    ) :: domain
  Integer,              Intent( In    ) :: lishp(:),lashp(1:domain%neighbours)
  Integer,              Intent( In    ) :: subtype
  Type( comms_type),    Intent( InOut ) :: comm
  Integer                               :: mpi_type

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message

  ! Number of transported quantities per particle

  iadd=4

  fail=0 ; limit=iadd*domain%mxbfsh ! limit=2*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'update_shared_units allocation failure'
     Call error(0,message)
  End If

  If(subtype == SHARED_UNIT_UPDATE_POSITIONS) Then
    mpi_type = comm%part_array_type_positions
  Else If( subtype == SHARED_UNIT_UPDATE_FORCES) Then
    mpi_type = comm%part_array_type_forces
  Else
    !! TODO Give error code, has to be either positions or forces
    Call error(0)
  End If
! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! transfer coordinate update data to all neighbouring nodes (26 surrounding me)

  j0=0
  Do k=1,26

! For all genuinely unique processors (no images or myself)

     If (domain%map_unique(k) == 0) Then

! Initialise the number of buffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),config%nlast,config%lsi,config%lsa)
              If (m > config%natms) m=0

! m should always be > 0 (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > 0) Then
                 buffer(i+1)=Real(lishp(j),wp)

                 If(subtype == SHARED_UNIT_UPDATE_POSITIONS) Then
                   buffer(i+2)=config%parts(m)%xxx
                   buffer(i+3)=config%parts(m)%yyy
                   buffer(i+4)=config%parts(m)%zzz
                 Else
                   buffer(i+2) = config%parts(m)%fxx
                   buffer(i+3) = config%parts(m)%fyy
                   buffer(i+4) = config%parts(m)%fzz
                 End If
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

        jdnode=domain%map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=domain%map(k0)

! pass only non-zero length messages
  
! transmit length of message

        n=0
        Call girecv(comm,n,kdnode,UpdShUnit_tag+k)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call gwait(comm)

        If (n > 0) Then
          Call girecv(comm,buffer(i+1:i+n),kdnode,UpdShUnit_tag+k)
        End If
        If (i > 0) Then
          Call gsend(comm,buffer(1:i),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call gwait(comm)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from buffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(Nint(buffer(i+1)),config%nlast,config%lsi,config%lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > config%natms) Then
                 If(subtype == SHARED_UNIT_UPDATE_POSITIONS) Then
                   config%parts(m)%xxx=buffer(i+2)
                   config%parts(m)%yyy=buffer(i+3)
                   config%parts(m)%zzz=buffer(i+4)
                 Else
                   config%parts(m)%fxx=buffer(i+2)
                   config%parts(m)%fyy=buffer(i+3)
                   config%parts(m)%fzz=buffer(i+4)
                 End If
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

End Subroutine update_shared_units_parts


Subroutine update_shared_units_arrays(config,lishp,lashp,qxx,qyy,qzz,domain,comm)

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

  Type( configuration_type ), Intent( InOut ) :: config
  Type( domains_type ), Intent( In    ) :: domain
  Integer,           Intent( In    ) :: lishp(:),lashp(1:domain%neighbours)
  !Real( Kind = wp ), Intent( InOut ) :: qxx(1:config%mxatms),qyy(1:config%mxatms),qzz(1:config%mxatms)
  Real( Kind = wp ), Intent( InOut ) :: qxx(:),qyy(:),qzz(:)
  Type( comms_type), Intent( InOut ) :: comm

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message
! Number of transported quantities per particle

  iadd=4

  fail=0 ; limit=iadd*domain%mxbfsh ! limit=2*iblock*iadd
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

     If (domain%map_unique(k) == 0) Then

! Initialise the number of buffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),config%nlast,config%lsi,config%lsa)
              If (m > config%natms) m=0

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

        jdnode=domain%map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=domain%map(k0)

! pass only non-zero length messages

! transmit length of message

        n=0
        Call girecv(comm,n,kdnode,UpdShUnit_tag+k)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call gwait(comm)

        If (n > 0) Then
          Call girecv(comm,buffer(i+1:i+n),kdnode,UpdShUnit_tag+k)
        End If
        If (i > 0) Then
          Call gsend(comm,buffer(1:i),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call gwait(comm)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from buffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(Nint(buffer(i+1)),config%nlast,config%lsi,config%lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > config%natms) Then
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

End Subroutine update_shared_units_arrays

Subroutine update_shared_units_int(config,lishp,lashp,iii,domain,comm)

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

  Type( configuration_type ), Intent( InOut ) :: config
  Type( domains_type ), Intent( In    ) :: domain
  Integer, Intent( In    ) :: lishp(:),lashp(1:domain%neighbours)
  Integer, Intent( InOut ) :: iii(:)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Integer, Dimension( : ), Allocatable :: ibuffer
  Character ( Len = 256 ) :: message
  
! Number of transported quantities per particle

  iadd=2

  fail=0 ; limit=iadd*domain%mxbfsh ! limit=2*iblock*iadd
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

     If (domain%map_unique(k) == 0) Then

! Initialise the number of ibuffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),config%nlast,config%lsi,config%lsa)
              If (m > config%natms) m=0

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

        jdnode=domain%map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=domain%map(k0)

! pass only non-zero length messages

! transmit length of message

        n=0
        Call girecv(comm,n,kdnode,UpdShUnit_tag+k)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call gwait(comm)

        If (n > 0) Then
          Call girecv(comm,ibuffer(i+1:i+n),kdnode,UpdShUnit_tag+k)
        End If
        If (i > 0) Then
          Call gsend(comm,ibuffer(1:i),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call gwait(comm)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from ibuffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(ibuffer(i+1),config%nlast,config%lsi,config%lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > config%natms) Then
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

Subroutine update_shared_units_rwp(config,lishp,lashp,rrr,domain,comm)

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

  Type( configuration_type ), Intent( In    ) :: config
  Type( domains_type ), Intent( In    ) :: domain
  Integer,           Intent( In    ) :: lishp(:),lashp(1:domain%neighbours)
  Real( Kind = wp ), Intent( InOut ) :: rrr(:)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character ( Len = 256 )  :: message

! Number of transported quantities per particle

  iadd=2

  fail=0 ; limit=iadd*domain%mxbfsh ! limit=2*iblock*iadd
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

     If (domain%map_unique(k) == 0) Then

! Initialise the number of buffer elements I am to send

        i=0
        Do j=j0+1,lashp(k)

! If no out of bound so far then carry on

           If (i+iadd <= iblock) Then
              m=local_index(lishp(j),config%nlast,config%lsi,config%lsa)
              If (m > config%natms) m=0

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

        jdnode=domain%map(k)
        If (Mod(k,2) == 1) Then
           k0=k+1
        Else
           k0=k-1
        End If
        kdnode=domain%map(k0)

! pass only non-zero length messages

! transmit length of message

        n=0
        Call girecv(comm,n,kdnode,UpdShUnit_tag+k)
        Call gsend(comm,i,jdnode,UpdShUnit_tag+k)
        Call gwait(comm)

        If (n > 0) Then
          Call girecv(comm,buffer(i+1:i+n),kdnode,UpdShUnit_tag+k)
        End If
        If (i > 0) Then
          Call gsend(comm,buffer(1:i),jdnode,UpdShUnit_tag+k)
        End If
        If (n > 0) Call gwait(comm)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from buffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(Nint(buffer(i+1)),config%nlast,config%lsi,config%lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > config%natms) Then
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
  Integer, Dimension( 0:, 1:), Intent( InOut) :: legend

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
