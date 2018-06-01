Module tethers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global tether interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,      Only : comms_type,gsum,gcheck
  Use configuration,     Only : imcon,cell,natms,nlast,lsi,lsa,lfrzn, &
                                xxx,yyy,zzz,fxx,fyy,fzz
  Use statistics, Only : stats_type
  Use setup, Only : mxtmls,mxtteth,mxteth,mxftet,mxpteth,mxatdm,nrite
  Use errors_warnings, Only : error, warning
  Use numerics, Only : images,local_index


  Implicit None
  Private

  Integer, Public,                   Save :: ntteth = 0


  Integer,Public,       Allocatable, Save :: numteth(:),keytet(:)
  Integer,Public,       Allocatable, Save :: lsttet(:),listtet(:,:),legtet(:,:)

  Real( Kind = wp ),Public, Allocatable, Save :: prmtet(:,:)

  Public :: allocate_tethers_arrays , deallocate_tethers_arrays
  Public :: tethers_forces

Contains

  Subroutine allocate_tethers_arrays()

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numteth(1:mxtmls),           Stat = fail(1))
    Allocate (keytet(1:mxtteth),           Stat = fail(2))
    Allocate (lsttet(1:mxtteth),           Stat = fail(3))
    Allocate (listtet(0:1,1:mxteth),       Stat = fail(4))
    Allocate (legtet(0:mxftet,1:mxatdm),   Stat = fail(5))
    Allocate (prmtet(1:mxpteth,1:mxtteth), Stat = fail(6))

    If (Any(fail > 0)) Call error(1017)

    numteth = 0
    keytet  = 0
    lsttet  = 0
    listtet = 0
    legtet  = 0

    prmtet  = 0.0_wp

  End Subroutine allocate_tethers_arrays

  Subroutine deallocate_tethers_arrays()

    Integer:: fail

    fail = 0

    Deallocate (numteth,lsttet, Stat = fail)

    If (fail > 0) Call error(1031)

  End Subroutine deallocate_tethers_arrays

  Subroutine tethers_forces(stats,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating energy and force terms for
! tethered particles
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Type( stats_type ), Intent( InOut ) :: stats
  Type( comms_type ), Intent( InOut ) :: comm

  Logical           :: safe
  Integer           :: fail(1:2),i,ia,kk
  Real( Kind = wp ) :: rab,rrab,fx,fy,fz,k,k2,k3,k4,rc,gamma,omega, &
                       strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Character( Len = 256 ) :: message
  fail=0
  Allocate (lstopt(0:1,1:mxteth),                         Stat=fail(1))
  Allocate (xdab(1:mxteth),ydab(1:mxteth),zdab(1:mxteth), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'tethers_forces allocation failure'
     Call error(0,message)
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
        xdab(i) = xxx(ia)-stats%xin(ia)
        ydab(i) = yyy(ia)-stats%yin(ia)
        zdab(i) = zzz(ia)-stats%zin(ia)
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

  stats%engtet=0.0_wp
  stats%virtet=0.0_wp

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

        stats%engtet=stats%engtet+omega
        stats%virtet=stats%virtet+gamma*rab*rab

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

  stats%stress(1) = stats%stress(1) + strs1
  stats%stress(2) = stats%stress(2) + strs2
  stats%stress(3) = stats%stress(3) + strs3
  stats%stress(4) = stats%stress(4) + strs2
  stats%stress(5) = stats%stress(5) + strs5
  stats%stress(6) = stats%stress(6) + strs6
  stats%stress(7) = stats%stress(7) + strs3
  stats%stress(8) = stats%stress(8) + strs6
  stats%stress(9) = stats%stress(9) + strs9

! check for undefined potentials

  Call gcheck(comm,safe)
  If (.not.safe) Call error(450)

! sum contributions to potential and virial

     buffer(1)=stats%engtet
     buffer(2)=stats%virtet
     Call gsum(comm,buffer(1:2))
     stats%engtet=buffer(1)
     stats%virtet=buffer(2)

  Deallocate (lstopt,         Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'tethers_forces deallocation failure'
     Call error(0,message)
  End If

End Subroutine tethers_forces

  
End Module tethers