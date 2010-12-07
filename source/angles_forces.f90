Subroutine angles_forces(imcon,engang,virang,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bond angle energy and force terms
!
! copyright - daresbury laboratory
! author    - w.smith may 1992
! amended   - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsync,gsum,gcheck
  Use setup_module,  Only : mxangl,nrite,pi
  Use config_module, Only : cell,natms,nlast,lsi,lsa,lfrzn, &
                            xxx,yyy,zzz,fxx,fyy,fzz
  Use angles_module, Only : ntangl,keyang,listang,prmang

  Implicit None

  Integer,                             Intent( In    ) :: imcon
  Real( Kind = wp ),                   Intent(   Out ) :: engang,virang
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Logical           :: safe
  Integer           :: fail(1:3),i,j,ia,ib,ic,keya,kk,local_index
  Real( Kind = wp ) :: xab,yab,zab,rab,rrab, xbc,ybc,zbc,rbc,rrbc, &
                       theta,cost,sint,rsint,                      &
                       fxa,fxc,fya, fyc,fza,fzc,                   &
                       k,k2,k3,k4,theta0,dtheta,dthpi,dth0pi,dth,  &
                       rho,rho1,rho2,switch,a,b,c,delta,m,dr1,dr2, &
                       pterm,vterm,gamma,gamsa,gamsc,              &
                       strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Real( Kind = wp ), Allocatable :: xdbc(:),ydbc(:),zdbc(:)

  fail=0
  Allocate (lunsafe(1:mxangl),lstopt(0:3,1:mxangl),       Stat=fail(1))
  Allocate (xdab(1:mxangl),ydab(1:mxangl),zdab(1:mxangl), Stat=fail(2))
  Allocate (xdbc(1:mxangl),ydbc(1:mxangl),zdbc(1:mxangl), Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'angles_forces allocation failure, node: ', idnode
     Call error(0)
  End If


! calculate atom separation vectors

  Do i=1,ntangl
     lunsafe(i)=.false.

! indices of angle bonded atoms

     ia=local_index(listang(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listang(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib
     ic=local_index(listang(3,i),nlast,lsi,lsa) ; lstopt(3,i)=ic

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0 .and. ic > 0) Then ! Tag
        If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic) == 0) Then
           If (ia <= natms .or. ib <= natms .or. ic <= natms) Then
              lstopt(0,i)=1
            End If
        End If
     Else                                       ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
              (ib > 0 .and. ib <= natms) .or.   &
              (ic > 0 .and. ic <= natms)) .and. &
             (ia == 0 .or. ib == 0 .or. ic == 0) ) lunsafe(i)=.true.
     End If

! define components of bond vectors

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)

        xdbc(i)=xxx(ic)-xxx(ib)
        ydbc(i)=yyy(ic)-yyy(ib)
        zdbc(i)=zzz(ic)-zzz(ib)
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
!
!        xdbc(i)=0.0_wp
!        ydbc(i)=0.0_wp
!        zdbc(i)=0.0_wp
     End If
  End Do

! Check for uncompressed units

  safe = .not. Any(lunsafe(1:ntangl))
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     Do j=0,mxnode-1
        If (idnode == j) Then
           Do i=1,ntangl
              If (lunsafe(i)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                 '*** warning - global unit number', listang(0,i), &
                 ' , with a head particle number', listang(1,i),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(130)
  End If

! periodic boundary condition

  Call images(imcon,cell,ntangl,xdab,ydab,zdab)
  Call images(imcon,cell,ntangl,xdbc,ydbc,zdbc)

! Initialise safety flag

  safe=.true.

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! zero angle energy accumulator

  engang=0.0_wp
  virang=0.0_wp

! loop over all specified angle potentials

  Do i=1,ntangl
     If (lstopt(0,i) > 0) Then

! indices of bonded atoms

        ia=lstopt(1,i)
        ib=lstopt(2,i)
        ic=lstopt(3,i)

! define components of first bond vector

        rab = Sqrt(xdab(i)**2+ydab(i)**2+zdab(i)**2)
        rrab = 1.0_wp/rab

        xab=xdab(i)*rrab
        yab=ydab(i)*rrab
        zab=zdab(i)*rrab

! define components of second bond vector

        rbc = Sqrt(xdbc(i)**2+ydbc(i)**2+zdbc(i)**2)
        rrbc = 1.0_wp/rbc

        xbc=xdbc(i)*rrbc
        ybc=ydbc(i)*rrbc
        zbc=zdbc(i)*rrbc

! index of potential function parameters

        kk=listang(0,i)

! determine bond angle and calculate potential energy

        cost=(xab*xbc+yab*ybc+zab*zbc)
        If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)
        theta=Acos(cost)
        sint=Max(1.0e-10_wp,Sqrt(1.0_wp-cost**2))
        rsint=1.0_wp/sint

        keya = Abs(keyang(kk))

        If      (keya == 1) Then

! harmonic potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0

           pterm=k*dtheta

           gamma=pterm*rsint
           pterm=pterm*0.5_wp*dtheta
           vterm=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp

        Else If (keya == 2) Then

! quartic potential

           k2    =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           k3    =prmang(3,kk)
           k4    =prmang(4,kk)

           gamma=dtheta*(k2+k3*dtheta+k4*dtheta**2)*rsint
           pterm=0.5_wp*k2*dtheta**2+(k3/3.0_wp)*dtheta**3+0.25*k4*dtheta**4
           vterm=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp

        Else If (keya == 3) Then

! truncated Harmonic potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho   =prmang(3,kk)
           switch=-(rab**8+rbc**8)/rho**8

           pterm=k*dtheta*Exp(switch)

           gamma=pterm*rsint
           pterm=pterm*0.5_wp*dtheta
           vterm=pterm*8.0_wp*switch
           gamsa=pterm*8.0_wp*rab**7/rho**8
           gamsc=pterm*8.0_wp*rbc**7/rho**8

        Else If (keya == 4) Then

! screened Harmonic potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho1  =prmang(3,kk)
           rho2  =prmang(4,kk)
           switch=-(rab/rho1+rbc/rho2)

           pterm=k*dtheta*Exp(switch)

           gamma=pterm*rsint
           pterm=pterm*0.5_wp*dtheta
           vterm=pterm*switch
           gamsa=pterm/rho1
           gamsc=pterm/rho2

        Else If (keya == 5) Then

! screened vessal potential (type 1)

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dth0pi=theta0-pi
           dthpi =theta -pi
           dth   =dth0pi**2-dthpi**2
           rho1  =prmang(3,kk)
           rho2  =prmang(4,kk)
           switch=-(rab/rho1+rbc/rho2)

           pterm=(k*dth/(2.0_wp*dth0pi**2)) * Exp(switch)

           gamma=pterm*dthpi*rsint
           pterm=pterm*0.25_wp*dth
           vterm=pterm*switch
           gamsa=pterm/rho1
           gamsc=pterm/rho2

        Else If (keya == 6) Then

! truncated Vessal potential (type 2)

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           dth0pi=theta0-pi
           dthpi =theta -pi
           a     =prmang(3,kk)
           rho   =prmang(4,kk)
           switch=-(rab**8+rbc**8)/rho**8

           pterm=k*dtheta*Exp(switch)

           gamma=pterm * (theta**(a-1.0_wp) * (dthpi+dth0pi) *                               &
                 ((a+4.0_wp)*theta**2 - 2.0_wp*pi*(a+2.0_wp)*theta - a*theta0*(dth0pi-pi)) + &
                 a*pi**(a-1.0_wp) * dth0pi**3) * rsint

           pterm=pterm * dtheta * (theta**a * (dthpi+dth0pi)**2 + 0.5_wp*a * pi**(a-1.0_wp) * dth0pi**3)
           vterm=pterm*8.0_wp*switch
           gamsa=pterm*8.0_wp*rab**7/rho**8
           gamsc=pterm*8.0_wp*rbc**7/rho**8

        Else If (keya == 7) Then

! harmonic cosine potential (note cancellation of sint in gamma)

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=Cos(theta)-Cos(theta0)

           pterm=k*dtheta

           gamma=-pterm
           pterm=pterm*0.5_wp*dtheta
           vterm=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp

        Else If (keya == 8) Then

! ordinary cosine potential

           k    =prmang(1,kk)
           delta=prmang(2,kk)
           m    =prmang(3,kk)
           a    =m*theta-delta

           gamma=k*m*Sin(a)*rsint
           pterm=k*(1.0_wp+Cos(a))
           vterm=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp

        Else If (keya == 9) Then

! MM3 stretch-bend potential

           a     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho1  =prmang(3,kk)
           dr1   =rab-rho1
           rho2  =prmang(4,kk)
           dr2   =rbc-rho2

           pterm=a*dr1*dr2

           gamma=pterm*rsint
           pterm=pterm*dtheta
           gamsa=-pterm/dr1
           gamsc=-pterm/dr2
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 10) Then

! compass stretch-stretch potential

           a     =prmang(1,kk)
           rho1  =prmang(2,kk)
           dr1   =rab-rho1
           rho2  =prmang(3,kk)
           dr2   =rbc-rho2

           pterm=a*dr1*dr2

           gamma=0.0_wp
           gamsa=-a*dr2
           gamsc=-a*dr1
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 11) Then

! compass stretch-bend potential

           a     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho1  =prmang(3,kk)
           dr1   =rab-rho1

           pterm=a*dr1

           gamma=pterm*rsint
           pterm=pterm*dtheta
           gamsa=-a*dtheta
           gamsc=0.0_wp
           vterm=-gamsa*rab

        Else If (keya == 12) Then

! combined compass angle potential with 3 coupling terms

           a     =prmang(1,kk)
           b     =prmang(2,kk)
           c     =prmang(3,kk)
           theta0=prmang(4,kk)
           dtheta=theta-theta0
           rho1  =prmang(5,kk)
           dr1   =rab-rho1
           rho2  =prmang(6,kk)
           dr2   =rbc-rho2

           m=b*dr1+c*dr2

           pterm=a*dr1*dr2 + dtheta*m

           gamma=k*rsint
           gamsa=-a*dr2-b*dtheta
           gamsc=-a*dr1-c*dtheta
           vterm=-(gamsa*rab+gamsc*rbc)

        Else

! undefined potential

           safe=.false.
           gamma=0.0_wp
           pterm=0.0_wp
           vterm=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp

        End If

! calculate atomic forces

        fxa = gamma*(xbc-xab*cost)*rrab+gamsa*xab
        fya = gamma*(ybc-yab*cost)*rrab+gamsa*yab
        fza = gamma*(zbc-zab*cost)*rrab+gamsa*zab

        fxc = gamma*(xab-xbc*cost)*rrbc+gamsc*xbc
        fyc = gamma*(yab-ybc*cost)*rrbc+gamsc*ybc
        fzc = gamma*(zab-zbc*cost)*rrbc+gamsc*zbc

        If (ia <= natms) Then

           fxx(ia)=fxx(ia)+fxa
           fyy(ia)=fyy(ia)+fya
           fzz(ia)=fzz(ia)+fza

        End If

        If (ib <= natms) Then

! energy and virial (associated to the head atom)

           engang=engang+pterm
           virang=virang+vterm

! calculate stress tensor (associated to the head atom)

           strs1 = strs1 + rab*xab*fxa + rbc*xbc*fxc
           strs2 = strs2 + rab*xab*fya + rbc*xbc*fyc
           strs3 = strs3 + rab*xab*fza + rbc*xbc*fzc
           strs5 = strs5 + rab*yab*fya + rbc*ybc*fyc
           strs6 = strs6 + rab*yab*fza + rbc*ybc*fzc
           strs9 = strs9 + rab*zab*fza + rbc*zbc*fzc

           fxx(ib)=fxx(ib)-fxa-fxc
           fyy(ib)=fyy(ib)-fya-fyc
           fzz(ib)=fzz(ib)-fza-fzc

        End If

        If (ic <= natms) Then

           fxx(ic)=fxx(ic)+fxc
           fyy(ic)=fyy(ic)+fyc
           fzz(ic)=fzz(ic)+fzc

        End If

     End If
  End Do

! check for undefined potentials

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(440)

! global sum of angular potential and virial

  If (mxnode > 1) Then
     buffer(1)=engang
     buffer(2)=virang
     Call gsum(buffer(1:2))
     engang=buffer(1)
     virang=buffer(2)
  End If

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

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  Deallocate (xdbc,ydbc,zdbc, Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'angles_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine angles_forces
