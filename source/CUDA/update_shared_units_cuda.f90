! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Center for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

#ifdef COMPILE_CUDA
Subroutine update_shared_units_constraints_shake_specific(natms,nlast,lsi,lsa,lishp,lashp,qxx,qyy,qzz)

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module, Only : map,mop
  Use dl_poly_cuda_module

  Implicit None

  Integer,           Intent( In    ) :: natms,nlast
  Integer,           Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer,           Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Real( Kind = wp ), Intent( InOut ) :: qxx(1:mxatms),qyy(1:mxatms),qzz(1:mxatms)

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n,local_index

  Integer :: count
  Logical :: hrcbc
  Integer, Dimension(:), Allocatable :: rilist
  Integer :: rilistcnt

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  Call start_timing_update_shared_units()
! Number of transported quantities per particle

  iadd=4

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units allocation failure, node: ', idnode
     Call error(0)
  End If

! cuda: flag to known if receive atom snooping is desired
  hrcbc = constraints_shake_cuda_has_recv_list_been_constructed()
  If (.not. hrcbc) Then
     Allocate(rilist(0:mxatms-1))
     rilistcnt = 0
  End If
  count = 0

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync()

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
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        Call MPI_SEND(i,1,MPI_INTEGER,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        Call MPI_WAIT(request,status,ierr)

        If (n > 0) Call MPI_IRECV(buffer(i+1),n,wp_mpi,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        If (i > 0) Call MPI_SEND(buffer(1),i,wp_mpi,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        If (n > 0) Call MPI_WAIT(request,status,ierr)

! consolidate transferred data
! I'm to receive data for n/iadd particles (n array elements from buffer(i+1))

        Do j=1,n/iadd
           If (i+iadd <= limit) Then
              m=local_index(Nint(buffer(i+1)),nlast,lsi,lsa)

! m should always be > natms (halo particles have local index > natms)
! (consistency - a particle strictly belongs to only one domain)

              If (m > natms) Then
                 If (.not. hrcbc) Then
                    rilist(rilistcnt) = m
                    rilistcnt = rilistcnt + 1
                 End If
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

  If (.not. hrcbc) Then
     Call constraints_shake_cuda_commit_receive_list(rilist,rilistcnt)
  End If

  Call gcheck(safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  If (.not. hrcbc) Then
     Deallocate(rilist)
  End If

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units deallocation failure, node: ', idnode
     Call error(0)
  End If
  Call stop_timing_update_shared_units()
End Subroutine update_shared_units_constraints_shake_specific
#endif

Subroutine update_shared_units(natms,nlast,lsi,lsa,lishp,lashp,qxx,qyy,qzz)

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

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module, Only : map,mop

#ifdef COMPILE_CUDA
  Use dl_poly_cuda_module
#endif

  Implicit None

  Integer,           Intent( In    ) :: natms,nlast
  Integer,           Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer,           Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Real( Kind = wp ), Intent( InOut ) :: qxx(1:mxatms),qyy(1:mxatms),qzz(1:mxatms)

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n,local_index

  Integer :: count

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

#ifdef COMPILE_CUDA
  Call start_timing_update_shared_units()
#endif

! Number of transported quantities per particle

  iadd=4

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units_int allocation failure, node: ', idnode
     Call error(0)
  End If

  count = 0

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync()

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
                 count = count + 1
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
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        Call MPI_SEND(i,1,MPI_INTEGER,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        Call MPI_WAIT(request,status,ierr)

        If (n > 0) Call MPI_IRECV(buffer(i+1),n,wp_mpi,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        If (i > 0) Call MPI_SEND(buffer(1),i,wp_mpi,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        If (n > 0) Call MPI_WAIT(request,status,ierr)

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

  Call gcheck(safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units deallocation failure, node: ', idnode
     Call error(0)
  End If

#ifdef COMPILE_CUDA
  Call stop_timing_update_shared_units()
#endif

End Subroutine update_shared_units

Subroutine update_shared_units_int(natms,nlast,lsi,lsa,lishp,lashp,iii)

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

  Use comms_module
  Use setup_module
  Use domains_module, Only : map,mop

  Implicit None

  Integer, Intent( In    ) :: natms,nlast
  Integer, Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer, Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Integer, Intent( InOut ) :: iii(1:mxatms)

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n,local_index

  Integer, Dimension( : ), Allocatable :: ibuffer

! Number of transported quantities per particle

  iadd=2

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (ibuffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units allocation failure, node: ', idnode
     Call error(0)
  End If

! Set ibuffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync()

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
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        Call MPI_SEND(i,1,MPI_INTEGER,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        Call MPI_WAIT(request,status,ierr)

        If (n > 0) Call MPI_IRECV(ibuffer(i+1),n,MPI_INTEGER,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        If (i > 0) Call MPI_SEND(ibuffer(1),i,MPI_INTEGER,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        If (n > 0) Call MPI_WAIT(request,status,ierr)

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

  Call gcheck(safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  Deallocate (ibuffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units_int deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine update_shared_units_int

Subroutine update_shared_units_rwp(natms,nlast,lsi,lsa,lishp,lashp,rrr)

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

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module, Only : map,mop

  Implicit None

  Integer,           Intent( In    ) :: natms,nlast
  Integer,           Intent( In    ) :: lsi(1:mxatms),lsa(1:mxatms)
  Integer,           Intent( In    ) :: lishp(1:mxlshp),lashp(1:mxproc)
  Real( Kind = wp ), Intent( InOut ) :: rrr(1:mxatms)

  Logical :: safe(1:2)
  Integer :: fail,iadd,limit,iblock, &
             i,j,k,j0,k0,jdnode,kdnode,m,n,local_index

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

! Number of transported quantities per particle

  iadd=2

  fail=0 ; limit=iadd*mxbfsh ! limit=2*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units_rwp allocation failure, node: ', idnode
     Call error(0)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/2

! Set logical flag for array overflow

  safe=.true.

! synchronise all processors

  Call gsync()

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
        Call MPI_IRECV(n,1,MPI_INTEGER,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        Call MPI_SEND(i,1,MPI_INTEGER,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        Call MPI_WAIT(request,status,ierr)

        If (n > 0) Call MPI_IRECV(buffer(i+1),n,wp_mpi,kdnode,Updshun_tag+k,dlp_comm_world,request,ierr)
        If (i > 0) Call MPI_SEND(buffer(1),i,wp_mpi,jdnode,Updshun_tag+k,dlp_comm_world,ierr)
        If (n > 0) Call MPI_WAIT(request,status,ierr)

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

  Call gcheck(safe)

! check for array overflow

  If (.not.safe(1)) Call error(115)

! check global error condition

  If (.not.safe(2)) Call error(116)

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'update_shared_units_rwp deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine update_shared_units_rwp