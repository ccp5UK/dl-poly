Subroutine pmf_vcoms(indpmf,xpmf,ypmf,zpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing PMF units' c.o.m. velocity
! vectors for iterative constraint (pmf_quench,pmf_rattle) algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module
  Use config_module, Only : natms,vxx,vyy,vzz
  Use pmf_module,    Only : ntpmf,listpmf,pmfwg1

  Implicit None

  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent(   Out ) :: xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf)

  Logical                 :: safe

  Integer                 :: fail,gpmf,gpmf1,gpmf2,ipmf,jpmf,iadd, &
                             j,k,l

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  Allocate (buffer(1:(mxtpmf(1)+mxtpmf(2))*(mxpmf+2)), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_vcoms allocation failure, node: ', idnode
     Call error(0)
  End If


! Initialise safety flag

  safe=.true.

! Initialise PMF c.o.m. momentum arrays

  xpmf = 0.0_wp ; ypmf = 0.0_wp ; zpmf = 0.0_wp

! Loop over all global PMF constraints

  gpmf1=1          ! number of passed global PMFs
  buffer=0.0_wp    ! velocities buffer
  iadd=3           ! adding three velocities only per particle
  Do gpmf=1,mxpmf

! Loop over all local to this node PMF constraints matching the global one

     Do ipmf=1,ntpmf
        If (listpmf(0,1,ipmf) == gpmf) Then

! Loop over all PMF units present on my domain (no halo)

           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then

! Loop over their members present on my domain (no halo)

                 Do k=1,mxtpmf(jpmf)
                    j=indpmf(k,jpmf,ipmf)

! Copy particles' velocities to buffer in an orderly manner

                    If (j > 0 .and. j <= natms) Then ! j is a domain particle
                       l=((gpmf-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1)+(k-1))*iadd
                       buffer(l+1)=vxx(j)
                       buffer(l+2)=vyy(j)
                       buffer(l+3)=vzz(j)
                    End If
                 End Do
              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe=(gpmf-gpmf1+2 < (mxpmf+2)/3)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe) .or. gpmf == mxpmf) Then
        If (mxnode > 1) Call gsum(buffer)

        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,ntpmf
              If (listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Get the COM momentum of this unit

                    Do k=1,mxtpmf(jpmf)
                       l=((gpmf2-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1)+(k-1))*iadd

                       xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf) + pmfwg1(k,jpmf)*buffer(l+1)
                       ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf) + pmfwg1(k,jpmf)*buffer(l+2)
                       zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf) + pmfwg1(k,jpmf)*buffer(l+3)
                    End Do

                    xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf)*pmfwg1(0,jpmf)
                    ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf)*pmfwg1(0,jpmf)
                    zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf)*pmfwg1(0,jpmf)

                 End Do
              End If
           End Do
        End Do

        gpmf1=gpmf2
        buffer=0.0_wp
     End If

  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_vcoms deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pmf_vcoms
