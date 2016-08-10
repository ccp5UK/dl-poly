Subroutine core_shell_forces(engshl,virshl,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating core-shell model spring energy
! and force terms
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsync,gsum,gcheck
  Use setup_module,      Only : mxshl,nrite
  Use config_module,     Only : imcon,cell,natms,nlast,lsi,lsa,xxx,yyy,zzz,fxx,fyy,fzz
  Use core_shell_module, Only : ntshl,listshl,prmshl

  Implicit None

  Real( Kind = wp ),                   Intent(   Out ) :: engshl,virshl
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Logical           :: safe
  Integer           :: fail(1:2),i,j,ia,ib,kk,local_index
  Real( Kind = wp ) :: rabsq,fx,fy,fz,gamma,omega,r_4_fac, &
                       strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)

  fail=0
  Allocate (lunsafe(1:mxshl),lstopt(0:2,1:mxshl),      Stat=fail(1))
  Allocate (xdab(1:mxshl),ydab(1:mxshl),zdab(1:mxshl), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_forces allocation failure, node: ', idnode
     Call error(0)
  End If

  r_4_fac = 1.0_wp/24.0_wp ! aharmonic shell coefficient = 1/(4!)

! calculate core-shell separation vectors

  Do i=1,ntshl
     lunsafe(i)=.false.

! indices of atoms in a core-shell

     ia=local_index(listshl(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listshl(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0) Then ! Tag
        If (ia <= natms .or. ib <= natms) Then
           lstopt(0,i)=1
        End If
     Else                          ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
              (ib > 0 .and. ib <= natms)) .and. &
             (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
     End If

! components of bond vector

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
     End If
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

! periodic boundary condition

  Call images(imcon,cell,ntshl,xdab,ydab,zdab)

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! zero core-shell energy and virial accumulators

  engshl=0.0_wp
  virshl=0.0_wp

! loop over all specified core-shell units

  Do i=1,ntshl
     If (lstopt(0,i) > 0) Then

! indices of atoms in a core-shell

        ia=lstopt(1,i)
        ib=lstopt(2,i)

! define components of bond vector

        rabsq = xdab(i)**2+ydab(i)**2+zdab(i)**2

! index of potential function parameters

        kk=listshl(0,i)

! calculate scalar constant terms using spring potential function
! and the parameters in array prmshl

        omega=(0.5_wp*prmshl(1,kk)+r_4_fac*prmshl(2,kk)*rabsq)*rabsq
        gamma=prmshl(1,kk)+prmshl(2,kk)*rabsq

! calculate forces

        fx = -gamma*xdab(i)
        fy = -gamma*ydab(i)
        fz = -gamma*zdab(i)

        If (ia <= natms) Then

           fxx(ia)=fxx(ia)+fx
           fyy(ia)=fyy(ia)+fy
           fzz(ia)=fzz(ia)+fz

! calculate core-shell unit energy

           engshl=engshl+omega
           virshl=virshl+gamma*rabsq

! calculate stress tensor

           strs1 = strs1 + xdab(i)*fx
           strs2 = strs2 + xdab(i)*fy
           strs3 = strs3 + xdab(i)*fz
           strs5 = strs5 + ydab(i)*fy
           strs6 = strs6 + ydab(i)*fz
           strs9 = strs9 + zdab(i)*fz

        End If

        If (ib <= natms) Then

           fxx(ib)=fxx(ib)-fx
           fyy(ib)=fyy(ib)-fy
           fzz(ib)=fzz(ib)-fz

        End If

     End If
  End Do

! complete stress tensor

  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9

! sum contributions to potential and virial

  If (mxnode > 1) Then
     buffer(1)=engshl
     buffer(2)=virshl
     Call gsum(buffer(1:2))
     engshl=buffer(1)
     virshl=buffer(2)
  End If

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine core_shell_forces
