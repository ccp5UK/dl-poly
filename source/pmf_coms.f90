Subroutine pmf_coms(imcon,indpmf,pxx,pyy,pzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing PMF units' COM vectors for
! iterative (pmf_quench,pmf_shake,pmf_rattle) constraint algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck,gsum
  Use setup_module
  Use config_module, Only : cell,natms,xxx,yyy,zzz
  Use pmf_module,    Only : ntpmf,listpmf,pmfwgt

  Implicit None

  Integer,           Intent( In    ) :: imcon
  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)

  Real( Kind = wp ), Intent(   Out ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)

  Logical                 :: safe(1:2)

  Real( Kind = wp )       :: celprp(1:10),width

  Integer                 :: fail(1:3),gpmf,gpmf1,gpmf2,ipmf,jpmf,iadd, &
                             j,k,l,m

  Real( Kind = wp ), Dimension( : ),    Allocatable :: xxt,yyt,zzt
  Real( Kind = wp ), Dimension( :, : ), Allocatable :: xpmf,ypmf,zpmf
  Real( Kind = wp ), Dimension( : ),    Allocatable :: buffer

  fail=0
  Allocate (xxt(1:Max(mxtpmf(1),mxtpmf(2))),yyt(1:Max(mxtpmf(1),mxtpmf(2))),zzt(1:Max(mxtpmf(1),mxtpmf(2))), Stat=fail(1))
  Allocate (xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf),                                           Stat=fail(2))
  Allocate (buffer(1:(mxtpmf(1)+mxtpmf(2))*(mxpmf+2)),                                                       Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_coms allocation failure, node: ', idnode
     Call error(0)
  End If


! Initialise safety flags

  safe=.true.

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)
  width=Min(celprp(7),celprp(8),celprp(9))

! Initialise PMF COMs' and inter-COMs' vector arrays

  xpmf = 0.0_wp ; ypmf = 0.0_wp ; zpmf = 0.0_wp
  pxx  = 0.0_wp ; pyy  = 0.0_wp ; pzz  = 0.0_wp

! Loop over all global PMF constraints

  gpmf1=1          ! number of passed global PMFs
  buffer=0.0_wp    ! coordinates buffer
  iadd=3           ! adding three coordinates only per particle
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

! Copy particles' coordinates to buffer in an orderly manner

                    If (j <= natms) Then
                       l=((gpmf-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1)+(k-1))*iadd
                       buffer(l+1)=xxx(j)
                       buffer(l+2)=yyy(j)
                       buffer(l+3)=zzz(j)
                    End If
                 End Do
              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe(1)=(gpmf-gpmf1+2 < (mxpmf+2)/3)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe(1)) .or. gpmf == mxpmf) Then
        If (mxnode > 1) Call gsum(buffer)

        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,ntpmf
              If (listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Get in the local scope of the unit

                    m=((gpmf2-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1))*iadd
                    Do k=1,mxtpmf(jpmf)
                       l=m+(k-1)*iadd

                       xxt(k)=buffer(l+1)-buffer(m+1)
                       yyt(k)=buffer(l+2)-buffer(m+2)
                       zzt(k)=buffer(l+3)-buffer(m+3)
                    End Do

                    Call images(imcon,cell,mxtpmf(jpmf),xxt,yyt,zzt)

                    Do k=1,mxtpmf(jpmf)
                       safe(2)=safe(2) .and. (Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2) < width/2.0_wp)
                    End Do

! Get the COM of this unit

                    Do k=1,mxtpmf(jpmf)
                       xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf) + pmfwgt(k,jpmf)*xxt(k)
                       ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf) + pmfwgt(k,jpmf)*yyt(k)
                       zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf) + pmfwgt(k,jpmf)*zzt(k)
                    End Do

! Get out of local frame

                    xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf)*pmfwgt(0,jpmf) + buffer(m+1)
                    ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf)*pmfwgt(0,jpmf) + buffer(m+2)
                    zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf)*pmfwgt(0,jpmf) + buffer(m+3)

                 End Do
              End If
           End Do
        End Do

        gpmf1=gpmf2
        buffer=0.0_wp
     End If

! Check if a PMF unit has a diameter > of the minumum of all half-cell width

     If (mxnode > 1) Call gcheck(safe(2))
     If (.not.safe(2)) Call error(492)

   End Do

! Loop over all PMF constraints on this node and calculate PMF constraint vectors

  Do ipmf=1,ntpmf
     pxx(ipmf) = xpmf(2,ipmf) - xpmf(1,ipmf)
     pyy(ipmf) = ypmf(2,ipmf) - ypmf(1,ipmf)
     pzz(ipmf) = zpmf(2,ipmf) - zpmf(1,ipmf)
  End Do

! Minimum image convention for bond vectors

  Call images(imcon,cell,ntpmf,pxx,pyy,pzz)

  Deallocate (xxt,yyt,zzt,    Stat=fail(1))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(2))
  Deallocate (buffer,         Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_coms deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pmf_coms

Subroutine pmf_update(indpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for updating PMF positions in the halo
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module
  Use config_module, Only : natms,xxx,yyy,zzz
  Use pmf_module,    Only : ntpmf,listpmf

  Implicit None

  Integer, Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)

  Logical :: safe
  Integer :: fail,pmfmem,iblock,gpmf,gpmf1,gpmf2,ipmf,jpmf,j,k,l,m

  Logical,           Allocatable :: lstpmf_sh(:,:)
  Real( Kind = wp ), Allocatable :: buffer(:)

  fail=0
  Allocate (lstpmf_sh(1:2,mxpmf),buffer(1:mxbuff), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_update allocation failure, node: ', idnode
     Call error(0)
  End If

  pmfmem=3*(mxtpmf(1)+mxtpmf(2))
  iblock=mxbuff/2

! Loop over all global PMF constraints

  lstpmf_sh=.false. ! no unit is shared initialisation
  gpmf1=1           ! number of passed global PMFs
  buffer=0.0_wp     ! coordinates buffer
  l=0               ! buffer indexing
  Do gpmf=1,mxpmf

! Loop over all local to this node PMF constraints matching the global one

     Do ipmf=1,ntpmf
        If (listpmf(0,1,ipmf) == gpmf) Then

! Loop over all shared PMF units present on my domain (no halo)

           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then

! Loop over their members present on my domain (no halo) and detect shared units

                 Do k=1,mxtpmf(jpmf)
                    j=indpmf(k,jpmf,ipmf)
                    If (j > natms) lstpmf_sh(jpmf,ipmf)=.true.
                 End Do

! Copy shared particles' coordinates to buffer in an orderly manner

                 If (lstpmf_sh(jpmf,ipmf)) Then
                    Do k=1,mxtpmf(jpmf)
                       j=indpmf(k,jpmf,ipmf)

                       If (j <= natms) Then
                          buffer(l+1)=xxx(j)
                          buffer(l+2)=yyy(j)
                          buffer(l+3)=zzz(j)
                       End If
                       l=l+3
                    End Do
                 End If

              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe=((l+(gpmf-gpmf1+1)*pmfmem) < iblock)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe) .or. gpmf == mxpmf) Then
        If (mxnode > 1) Call gsum(buffer)

        m=0
        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,ntpmf
              If (listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Recover shared particles' coordinates from buffer in an orderly manner

                    If (lstpmf_sh(jpmf,ipmf)) Then
                       Do k=1,mxtpmf(jpmf)
                          j=indpmf(k,jpmf,ipmf)

                          If (j > natms) Then
                             xxx(j)=buffer(m+1)
                             yyy(j)=buffer(m+2)
                             zzz(j)=buffer(m+3)
                          End If
                          m=m+3
                       End Do
                    End If

                 End Do
              End If
           End Do
        End Do

        gpmf1=gpmf2
        buffer=0.0_wp
        l=0
     End If

   End Do

  Deallocate (lstpmf_sh,buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_update deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pmf_update
