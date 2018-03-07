Subroutine pmf_units_set()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting the existence of the PMF units on the
! node where supposedly PMF constraints exist
!
! Note: (1) Deals with listpmf, legpmf and ntpmf
!       (2) Applies only at the end of relocate_particles
!           if megpmf>0 and mxnode>1
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,  Only : idnode,gcheck
  Use setup_module

  Use configuration, Only : natms,lsi,lsa
  Use pmf_module

  Implicit None

  Logical :: safe,ok
  Integer :: fail,ipmf,jpmf,gpmf,local_index

  Integer, Dimension( : ), Allocatable :: i1pmf0,i2pmf0

  fail=0
  Allocate (i1pmf0(1:mxtpmf(1)),i2pmf0(1:mxtpmf(2)), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_units_set allocation failure, node: ', idnode
     Call error(0)
  End If

! Initialise safety flag

  safe=.true.

! is it ok not to compress the bookkeeping list arrays
! since it's safe - there's enough buffering space

  ok=.true.
  If (mxpmf > 0) ok=.not.(Real(ntpmf,wp)/Real(mxpmf,wp) > 0.85_wp)

! Sort out local PMF units presence and legend array

  ipmf=0
  Do While (ipmf < ntpmf)
     ipmf=ipmf+1
10   Continue

! This holds the global PMF index

     gpmf=listpmf(0,1,ipmf)
     If (gpmf > 0) Then

! For presence of : PMF unit 1 only - this holds 1
!                   PMF unit 2 only - this holds 2
!                   both units 1&2  - this holds 3
! It CANNOT and MUST NOT hold ZERO

        listpmf(0,2,ipmf)=0

        Do jpmf=1,mxtpmf(1)
           i1pmf0(jpmf)=local_index(listpmf(jpmf,1,ipmf),natms,lsi,lsa)

           If (i1pmf0(jpmf) > 0) Then

! This identifies which local PMF constraint the particle belongs to

              legpmf(0,i1pmf0(jpmf))=1
              legpmf(1,i1pmf0(jpmf))=ipmf

           End If
        End Do
        If (Any(i1pmf0 > 0)) listpmf(0,2,ipmf)=listpmf(0,2,ipmf)+1

        Do jpmf=1,mxtpmf(2)
           i2pmf0(jpmf)=local_index(listpmf(jpmf,2,ipmf),natms,lsi,lsa)

           If (i2pmf0(jpmf) > 0) Then

! This identifies which local PMF constraint the particle belongs to

              legpmf(0,i2pmf0(jpmf))=1
              legpmf(1,i2pmf0(jpmf))=ipmf

           End If
        End Do
        If (Any(i2pmf0 > 0)) listpmf(0,2,ipmf)=listpmf(0,2,ipmf)+2

! If the PMF has moved to another node, listpmf(0,2,ipmf)=0,
! compress listpmf

        If (listpmf(0,2,ipmf) == 0 .and. (.not.ok)) Then
           If      (ipmf  < ntpmf) Then
              listpmf(:,:,ipmf)=listpmf(:,:,ntpmf) ! Copy list content from 'ipmf' to 'ntpmf'
              listpmf(:,:,ntpmf)=0                 ! Remove list content in 'ntpmf'
              ntpmf=ntpmf-1                        ! Reduce 'ntpmf' pointer

              Go To 10 ! Go back and check it all again for the new list content in 'ipmf'
           Else If (ipmf == ntpmf) Then
              listpmf(:,:,ntpmf)=0                 ! Remove list content in 'ntpmf=ipmf'
              ntpmf=ntpmf-1                        ! Reduce 'ntpmf' pointer
           End If
        End If

     Else

        safe=.false.

     End If
  End Do

! check for inconsistently built local list

  Call gcheck(safe)
  If (.not.safe) Call error(490)

  Deallocate (i1pmf0,i2pmf0, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_units_set deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pmf_units_set
