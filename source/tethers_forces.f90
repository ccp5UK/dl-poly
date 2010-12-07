Subroutine tethers_forces(imcon,engtet,virtet,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating energy and force terms for
! tethered particles
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsum,gcheck
  Use setup_module,      Only : mxteth,nrite
  Use config_module,     Only : cell,natms,nlast,lsi,lsa,lfrzn, &
                                xxx,yyy,zzz,fxx,fyy,fzz
  Use tethers_module,    Only : ntteth,keytet,listtet,prmtet
  Use statistics_module, Only : xin,yin,zin

  Implicit None

  Integer,                             Intent( In    ) :: imcon
  Real( Kind = wp ),                   Intent(   Out ) :: engtet,virtet
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Logical           :: safe
  Integer           :: fail(1:2),i,ia,kk,local_index
  Real( Kind = wp ) :: rab,rrab,fx,fy,fz,k,k2,k3,k4,rc,gamma,omega, &
                       strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)

  fail=0
  Allocate (lstopt(0:1,1:mxteth),                         Stat=fail(1))
  Allocate (xdab(1:mxteth),ydab(1:mxteth),zdab(1:mxteth), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'tethers_forces allocation failure, node: ', idnode
     Call error(0)
  End If

! Initialise safety flag

  safe=.true.

! calculate tether vector

  Do i=1,ntteth

! indices of tethered atoms

     ia=local_index(listtet(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia

     lstopt(0,i)=0
     If (ia > 0 .and. ia <= natms) Then !Tag
        If (lfrzn(ia) == 0) lstopt(0,i)=1
     End If

! There is no need to check for uncompressed unit since
! a tether is a point, it is either in or out by construction

! components of tether vector

     If (lstopt(0,i) > 0) Then
        xdab(i) = xxx(ia)-xin(ia)
        ydab(i) = yyy(ia)-yin(ia)
        zdab(i) = zzz(ia)-zin(ia)
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
     End If

  End Do

! periodic boundary condition

  Call images(imcon,cell,ntteth,xdab,ydab,zdab)

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! zero tether energy and virial accumulators

  engtet=0.0_wp
  virtet=0.0_wp

! loop over all specified tethered atoms

  Do i=1,ntteth
     If (lstopt(0,i) > 0) Then

! indices of tethered atoms

        ia=lstopt(1,i)

! define components of bond vector

        rab = Sqrt(xdab(i)**2+ydab(i)**2+zdab(i)**2)

! check for possible zero length vector

        If (rab < 1.0e-10_wp) Then
          rrab=0.0_wp
        Else
          rrab = 1.0_wp/rab
        End If

! index of potential function parameters

        kk=listtet(0,i)

! calculate scalar constant terms

        If      (keytet(kk) == 1) Then

! harmonic function

           k=prmtet(1,kk)

           omega=0.5_wp*k*rab**2
           gamma=k

        Else If (keytet(kk) == 2) Then

! restrained harmonic

           k =prmtet(1,kk)
           rc=prmtet(2,kk)

           omega=0.5_wp*k*Min(rab,rc)**2 + k*rc*Sign(Max(rab-rc,0.0_wp),rab)
           gamma=k*Sign(Min(rab,rc),rab)*rrab

        Else If (keytet(kk) == 3) Then

! quartic potential

           k2=prmtet(1,kk)
           k3=prmtet(2,kk)
           k4=prmtet(3,kk)

           omega=rab**2 * (0.5_wp*k2+(k3/3.0_wp)*rab+0.25_wp*k4*rab**2)
           gamma=k2+k3*rab+k4*rab**2

        Else

! undefined potential

           safe=.false.
           omega=0.0_wp
           gamma=0.0_wp

        End If

! calculate forces

        fx = -gamma*xdab(i)
        fy = -gamma*ydab(i)
        fz = -gamma*zdab(i)


        fxx(ia)=fxx(ia)+fx
        fyy(ia)=fyy(ia)+fy
        fzz(ia)=fzz(ia)+fz

! calculate tether energy and virial

        engtet=engtet+omega
        virtet=virtet+gamma*rab*rab

! calculate stress tensor

        strs1 = strs1 + xdab(i)*fx
        strs2 = strs2 + xdab(i)*fy
        strs3 = strs3 + xdab(i)*fz
        strs5 = strs5 + ydab(i)*fy
        strs6 = strs6 + ydab(i)*fz
        strs9 = strs9 + zdab(i)*fz

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

! check for undefined potentials

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(450)

! sum contributions to potential and virial

  If (mxnode > 1) Then
     buffer(1)=engtet
     buffer(2)=virtet
     Call gsum(buffer(1:2))
     engtet=buffer(1)
     virtet=buffer(2)
  End If

  Deallocate (lstopt,         Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'tethers_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine tethers_forces
