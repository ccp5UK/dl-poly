Subroutine bonds_forces(isw,engbnd,virbnd,stress,rcut,keyfce,alpha,epsq,engcpe,vircpe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating chemical bond energy and force
! terms
!
! isw = 0 - collect statistics
! isw = 1 - calculate forces
! isw = 2 - do both
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov march 2016
! contrib   - a.v.brukhno & i.t.todorov april 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsync,gsum,gcheck
  Use setup_module,  Only : nrite,r4pie0,zero_plus,mximpl,mxbond,mxgbnd1
  Use config_module, Only : imcon,cell,natms,nlast,lsi,lsa,lfrzn, &
                            chge,xxx,yyy,zzz,fxx,fyy,fzz
  Use bonds_module,  Only : ntbond,keybnd,listbnd,prmbnd, &
                            ltpbnd,vbnd,gbnd,ncfbnd,rcbnd,ldfbnd,dstbnd

  Implicit None

  Integer,                             Intent( In    ) :: isw
  Real( Kind = wp ),                   Intent(   Out ) :: engbnd,virbnd
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
  Real( Kind = wp ),                   Intent( In    ) :: rcut,alpha,epsq
  Integer,                             Intent( In    ) :: keyfce
  Real( Kind = wp ),                   Intent( InOut ) :: engcpe,vircpe

  Logical           :: safe(1:3)
  Integer           :: fail(1:2),i,j,l,ia,ib,keyb,kk,local_index
  Real( Kind = wp ) :: rab,rab2,fx,fy,fz,gamma,omega,    &
                       term,term1,term2,eps,sig,         &
                       k,k2,k3,k4,r0,dr,dra,dr2,         &
                       e0,rc,a,b,c,rho,delta,chgprd,     &
                       rdelr,rdr,ppp,vk,vk1,vk2,t1,t2,   &
                       viracc,engc12,virc12,buffer(1:4), &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)

  fail=0
  Allocate (lunsafe(1:mxbond),lstopt(0:2,1:mxbond),       Stat=fail(1))
  Allocate (xdab(1:mxbond),ydab(1:mxbond),zdab(1:mxbond), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'bond_forces allocation failure, node: ', idnode
     Call error(0)
  End If


! calculate atom separation vectors

  Do i=1,ntbond
     lunsafe(i)=.false.

! indices of bonded atoms

     ia=local_index(listbnd(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listbnd(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0) Then ! Tag
        If (lfrzn(ia)*lfrzn(ib) == 0) Then
           If (ia <= natms .or. ib <= natms) Then
              lstopt(0,i)=1
           End If
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

  safe(1) = .not. Any(lunsafe(1:ntbond))
  If (mxnode > 1) Call gcheck(safe(1))
  If (.not.safe(1)) Then
     Do j=0,mxnode-1
        If (idnode == j) Then
           Do i=1,ntbond
              If (lunsafe(i)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                 '*** warning - global unit number', listbnd(0,i), &
                 ' , with a head particle number', listbnd(1,i),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(128)
  End If

! periodic boundary condition

  Call images(imcon,cell,ntbond,xdab,ydab,zdab)

! Initialise safety flag

  safe=.true.

  If (Mod(isw,3) > 0) Then

! zero bond energy and virial accumulators

     engbnd=0.0_wp
     virbnd=0.0_wp

! zero scaled 1-2 electrostatic potential accumulators

     engc12=0.0_wp
     virc12=0.0_wp

     viracc=0.0_wp

! initialise stress tensor accumulators

     strs1=0.0_wp
     strs2=0.0_wp
     strs3=0.0_wp
     strs5=0.0_wp
     strs6=0.0_wp
     strs9=0.0_wp

  End If

! Recover bin size and increment counter

  If (Mod(isw,2) == 0) Then
     rdelr  = Real(mxgbnd1,wp)/rcbnd
     ncfbnd = ncfbnd + 1
  End If

! loop over all specified chemical bond potentials

  Do i=1,ntbond
     If (lstopt(0,i) > 0) Then

! indices of bonded atoms

        ia=lstopt(1,i)
        ib=lstopt(2,i)

! define components of bond vector

        rab2 = xdab(i)**2+ydab(i)**2+zdab(i)**2
        rab  = Sqrt(rab2)

! index of potential function parameters

        kk=listbnd(0,i)
        keyb = Abs(keybnd(kk))

! accumulate the histogram (distribution)

        If (Mod(isw,2) == 0 .and. ia <= natms) Then
           j = ldfbnd(kk)
           l = Min(1+Int(rab*rdelr),mxgbnd1)

           dstbnd(l,j) = dstbnd(l,j) + 1.0_wp

           If (rab > rcbnd) safe(3)=.false. ! catch bondbreaking
        End If
        If (isw == 0) Cycle

! calculate scalar constant terms

        If      (keyb == 0) Then

! null interaction

           omega=0.0_wp
           gamma=0.0_wp

        Else If (keyb == 1) Then

! harmonic potential

           k =prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0

           term=k*dr

           omega=term*0.5_wp*dr
           gamma=-term/rab

        Else If (keyb == 2) Then

! Morse potential

           e0=prmbnd(1,kk)
           r0=prmbnd(2,kk)
           k =prmbnd(3,kk)

           term=Exp(-k*(rab-r0))

           omega=e0*term*(term-2.0_wp)
           gamma=-2.0_wp*e0*k*term*(1.0_wp-term)/rab

        Else If (keyb == 3) Then

! 12-6 potential

           a=prmbnd(1,kk)
           b=prmbnd(2,kk)

           term=rab**(-6)

           omega=term*(a*term-b)
           gamma=6.0_wp*term*(2.0_wp*a*term-b)/rab2

        Else If (keyb == 4) Then

! Lennard-Jones potential

           eps=prmbnd(1,kk)
           sig=prmbnd(2,kk)

           term=(sig/rab)**6

           omega=4.0_wp*eps*term*(term-1.0_wp)
           gamma=24.0_wp*eps*term*(2.0_wp*term-1.0_wp)/rab2

        Else If (keyb == 5) Then

! restrained harmonic

           k =prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0
           dra=Abs(dr)
           rc=prmbnd(3,kk)

           omega=k*(0.5_wp*Min(dra,rc)**2 + rc*Max(dra-rc,0.0_wp))
           gamma=-k*Sign(Min(dra,rc),dr)/rab

        Else If (keyb == 6) Then

! quartic potential

           k2=prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0
           k3=prmbnd(3,kk)
           k4=prmbnd(4,kk)

           dr2=dr**2

           omega=dr2 * (0.5_wp*k2+k3*dr/3.0_wp+0.25_wp*k4*dr2)
           gamma=-dr*(k2+k3*dr+k4*dr2)/rab

        Else If (keyb == 7) Then

! Buckingham exp-6 potential

           a  =prmbnd(1,kk)
           rho=prmbnd(2,kk)
           c  =prmbnd(3,kk)

           term1=a*Exp(-rab/rho)
           term2=-c/rab**6

           omega=term1+term2
           gamma=(term1/rho+6.0_wp*term2/rab)/rab

        Else If (keyb == 8) Then

           omega=0.0_wp
           gamma=0.0_wp

! scaled charge product times dielectric constants

           chgprd=prmbnd(1,kk)*chge(ia)*chge(ib)*r4pie0/epsq
           If ((Abs(chgprd) > zero_plus .or. mximpl > 0) .and. keyfce > 0) Then
              If (mximpl > 0) Then
                 Call intra_mcoul(keyfce,rcut,alpha,epsq,ia,ib,chgprd, &
                      rab,xdab(i),ydab(i),zdab(i),omega,viracc,fx,fy,fz,safe(1))
              Else
                 Call intra_coul(keyfce,rcut,alpha,epsq,chgprd,rab,rab2,omega,gamma,safe(1))

                 fx = gamma*xdab(i)
                 fy = gamma*xdab(i)
                 fz = gamma*xdab(i)

                 viracc = -gamma*rab2
              End If

! correct electrostatic energy and virial

              If (ia <= natms) Then
                 engc12 = engc12 + omega
                 virc12 = virc12 + viracc
              End If
           End If

        Else If (keyb == 9) Then

! extended FENE (Finite Extensive Non-linear Elastic) potential

           k    =prmbnd(1,kk)
           r0   =prmbnd(2,kk)
           delta=prmbnd(3,kk)
           dr   =rab-delta

           term =1.0_wp-(dr/r0)**2
           If (term > 0.0_wp) Then
              omega=-0.5_wp*k*r0**2 * Log(term)
              gamma=-k*dr/term/rab
           Else
              safe(2)=.false.
              omega=0.0_wp
              gamma=0.0_wp
           End If

        Else If (keyb == 10) Then

! MM3-bond-stretch potential

           k =prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0

           e0=2.55_wp*dr

           term=k*dr*(1.0_wp-e0*(1.0_wp-(7.0_wp/12.0_wp)*e0))

           omega=term*dr
           gamma=(4.0_wp*term+k*dr*(-2.0_wp+e0))/rab

        Else If (keyb == 20) Then

! TABBND potential

           j = ltpbnd(kk)
           If (rab <= vbnd(-1,j)) Then ! rab <= cutpot
              rdr = gbnd(-1,j) ! 1.0_wp/delpot

              l   = Int(rab*rdr)
              ppp = rab*rdr - Real(l,wp)

              vk  = vbnd(l,j)
              vk1 = vbnd(l+1,j)
              vk2 = vbnd(l+2,j)

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              omega = t1 + (t2-t1)*ppp*0.5_wp

              vk  = gbnd(l,j) ; If (l == 0) vk = vk*rab
              vk1 = gbnd(l+1,j)
              vk2 = gbnd(l+2,j)

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rab
           Else ! bond breaking
              safe(3)=.false.
              omega=0.0_wp
              gamma=0.0_wp
           End If

        Else

! undefined potential

           safe(1)=.false.
           omega=0.0_wp
           gamma=0.0_wp

        End If

! calculate forces and virial additions

        If (keyb /= 8) Then
           fx = gamma*xdab(i)
           fy = gamma*ydab(i)
           fz = gamma*zdab(i)

           viracc=-gamma*rab2
        End If

! add forces

        If (ia <= natms) Then

           fxx(ia)=fxx(ia)+fx
           fyy(ia)=fyy(ia)+fy
           fzz(ia)=fzz(ia)+fz

! calculate bond energy and virial

           engbnd=engbnd+omega
           virbnd=virbnd+viracc

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

  If (Mod(isw,3) > 0) Then

! sum contributions to potential and virial

     If (mxnode > 1) Then
        buffer(1) = engbnd
        buffer(2) = virbnd
        buffer(3) = engc12
        buffer(4) = virc12
        Call gsum(buffer(1:4))
        engbnd = buffer(1)
        virbnd = buffer(2)
        engc12 = buffer(3)
        virc12 = buffer(4)
     End If

     engbnd = engbnd - engc12
     virbnd = virbnd - virc12

     engcpe = engcpe + engc12
     vircpe = vircpe + virc12

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

  End If

! check for undefined potentials + extras

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe(1)) Call error(444)
  If (.not.safe(2)) Call error(655)
  If (.not.safe(3)) Call error(660)

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'bond_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine bonds_forces
