Subroutine vaf_collect(lvafav,leql,nsteql,nstep,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistics for velocity
! autocorrelation functions
!
! copyright - daresbury laboratory
! author    - m.a.seaton june 2014
! amended   - i.t.todorov july 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,     Only : idnode,mxnode,gsum
  Use setup_module,     Only : nrite,mxatyp,mxbuff
  Use configuration,    Only : natms,ltype,lfrzn,vxx,vyy,vzz
  Use greenkubo_module, Only : isvaf,nsvaf,vaftsts,vafsamp,vafcount, &
                               vafstep,vxi,vyi,vzi,vafdata,vaftime,vaf

  Implicit None

  Logical,           Intent( In    ) :: lvafav,leql
  Integer,           Intent( In    ) :: nsteql,nstep
  Real( Kind = wp ), Intent( In    ) :: time

  Integer                        :: fail,i,j,k,l,nsum

  Real( Kind = wp ), Allocatable :: vafcoll(:)

  fail=0
  Allocate (vafcoll(1:mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vaf_collect allocation failure, node: ', idnode
     Call error(0)
  End If

! set VAF timestep start

  If (vaftsts < 0) vaftsts=Merge(nsteql,nstep,leql)

! set reference velocities only for first sampling passes

  If (nstep <= vaftsts+(vafsamp-1)*isvaf) Then
     Do j=1,vafsamp
        If (nstep == (vaftsts+(j-1)*isvaf)) Then
           vxi(1:natms,j) = vxx(1:natms)
           vyi(1:natms,j) = vyy(1:natms)
           vzi(1:natms,j) = vzz(1:natms)
        End If
     End Do
  End If

! advance time step, calculate vaf contribution for each sample

  Do j=1,vafsamp
     vafstep(j) = vafstep(j) + 1
     k=vafstep(j)

     If (k >= 0 .and. k <= nsvaf) Then
        vafcoll=0.0_wp
        Do i=1,natms
           If (lfrzn(i) == 0) Then
              l=ltype(i)
              vafcoll(l) = vafcoll(l) + vxx(i)*vxi(i,j)+vyy(i)*vyi(i,j)+vzz(i)*vzi(i,j)
           End If
        End Do

        Do l=1,mxatyp
           vafdata(k,(j-1)*(mxatyp+1)+l) = vafcoll(l)
        End Do

        vafdata(k,j*(mxatyp+1)) = time
     End If

     If (k == nsvaf) Then
        vafcount = vafcount + 1.0_wp

        If (mxnode > 1) Then
           nsum = mxbuff/(nsvaf+1)

           l=(j-1)*(mxatyp+1) ! avoid summing up timing information
           Do i=1,mxatyp,nsum
              Call gsum(vafdata(:,l+i:l+Min(i+nsum-1,mxatyp)))
           End Do
        End If

        If (lvafav) Then ! if time-averaging, add vaf data to sampling array
           Do i=0,nsvaf
              vaftime(i) = vafdata(i,j*(mxatyp+1))
              Do l=1,mxatyp
                 vaf(i,l) = vaf(i,l) + vafdata(i,(j-1)*(mxatyp+1)+l)
              End Do
           End Do
        Else             ! if not time-averaging, move vaf data to sampling array
           vaftime(0:nsvaf) = vafdata(0:nsvaf,j*(mxatyp+1))
           Do l=1,mxatyp
              vaf(0:nsvaf,l) = vafdata(0:nsvaf,(j-1)*(mxatyp+1)+l)
           End Do
        End If
     End If

! reset counter and reference velocities and get vaf
! at end of time span or at start of sample

     If (k == nsvaf .or. (isvaf > nsvaf .and. k == isvaf)) Then
        vafstep(j) = 0

        vxi(1:natms,j) = vxx(1:natms)
        vyi(1:natms,j) = vyy(1:natms)
        vzi(1:natms,j) = vzz(1:natms)

        vafcoll = 0.0_wp
        Do i=1,natms
           If (lfrzn(i) == 0) Then
              l=ltype(i)
              vafcoll(l) = vafcoll(l) + vxi(i,j)*vxi(i,j)+vyi(i,j)*vyi(i,j)+vzi(i,j)*vzi(i,j)
           End If
        End Do

        vafdata(0,j*(mxatyp+1)) = time
        Do l=1,mxatyp
           vafdata(0,(j-1)*(mxatyp+1)+l) = vafcoll(l)
        End Do
     End If
  End Do

  Deallocate (vafcoll, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vaf_collect deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine vaf_collect
