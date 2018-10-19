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

  Use kinds, Only : wp, wi
  Use comms,      Only : comms_type,gsum,gcheck
  Use configuration,     Only : configuration_type
  Use particle,   Only : corePart
  Use statistics, Only : stats_type
  Use errors_warnings, Only : error, warning
  Use numerics, Only : images,local_index
  Implicit None

  Private

  !> Type containing tethers data
  Type, Public :: tethers_type
    Private
    !> Number of types of tether
    Integer( Kind =  wi ), Public :: ntteth = 0
    Integer( Kind =  wi ), Public :: mxtteth,mxteth,mxftet,mxpteth
    Integer( Kind =  wi ), Allocatable, Public :: numteth(:),keytet(:)
    Integer( Kind =  wi ), Allocatable, Public :: lsttet(:),listtet(:,:),legtet(:,:)

    !> Tether parameters
    Real( Kind = wp ), Allocatable, Public :: prmtet(:,:)

    !> Total number of tehters
    Integer( Kind = wi ), Public :: total
  Contains
    Private

    Procedure, Public :: init => allocate_tethers_arrays
    Procedure, Public :: deallocate_temp => deallocate_tethers_arrays
    Final :: cleanup
  End Type

  Public :: tethers_forces

Contains

  Subroutine allocate_tethers_arrays(tether,mxtmls,mxatdm)
  Class(tethers_type) :: tether
    Integer( Kind = wi ), Intent( In    ) :: mxtmls,mxatdm

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (tether%numteth(1:mxtmls),           Stat = fail(1))
    Allocate (tether%keytet(1:tether%mxtteth),   Stat = fail(2))
    Allocate (tether%lsttet(1:tether%mxtteth),   Stat = fail(3))
    Allocate (tether%listtet(0:1,1:tether%mxteth),       Stat = fail(4))
    Allocate (tether%legtet(0:tether%mxftet,1:mxatdm),   Stat = fail(5))
    Allocate (tether%prmtet(1:tether%mxpteth,1:tether%mxtteth), Stat = fail(6))

    If (Any(fail > 0)) Call error(1017)

    tether%numteth = 0
    tether%keytet  = 0
    tether%lsttet  = 0
    tether%listtet = 0
    tether%legtet  = 0

    tether%prmtet  = 0.0_wp
  End Subroutine allocate_tethers_arrays

  Subroutine deallocate_tethers_arrays(tether)
  Class(tethers_type) :: tether

    Integer :: fail(1:2)

    fail = 0

    If (Allocated(tether%numteth)) Then
      Deallocate (tether%numteth, Stat = fail(1))
    End If
    If (Allocated(tether%lsttet)) Then
      Deallocate (tether%lsttet, Stat = fail(2))
    End If

    If (Any(fail > 0)) Call error(1031)
  End Subroutine deallocate_tethers_arrays

  Subroutine cleanup(tether)
    Type( tethers_type ) :: tether

    If (Allocated(tether%numteth)) Then
      Deallocate (tether%numteth)
    End If

    If (Allocated(tether%keytet)) Then
      Deallocate (tether%keytet)
    End If

    If (Allocated(tether%lsttet)) Then
      Deallocate (tether%lsttet)
    End If
    If (Allocated(tether%listtet)) Then
      Deallocate (tether%listtet)
    End If
    If (Allocated(tether%legtet)) Then
      Deallocate (tether%legtet)
    End If
  End Subroutine cleanup

  Subroutine tethers_forces(stats,tether,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating energy and force terms for
    ! tethered particles
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type( stats_type )   , Intent( InOut ) :: stats
    Type( tethers_type ) , Intent( InOut ) :: tether
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type )   , Intent( InOut ) :: comm

    Logical           :: safe
    Integer           :: fail(1:2),i,ia,kk
    Real( Kind = wp ) :: rab,rrab,fx,fy,fz,k,k2,k3,k4,rc,gamma,omega, &
      strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

    Integer,           Allocatable :: lstopt(:,:)
    Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
    Character( Len = 256 ) :: message
    fail=0

    Allocate (lstopt(0:1,1:tether%mxteth),                         Stat=fail(1))
    Allocate (xdab(1:tether%mxteth),ydab(1:tether%mxteth),zdab(1:tether%mxteth), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'tethers_forces allocation failure'
      Call error(0,message)
    End If

    ! Initialise safety flag

    safe=.true.

    ! calculate tether vector

    Do i=1,tether%ntteth

      ! indices of tethered atoms

      ia=local_index(tether%listtet(1,i),config%nlast,config%lsi,config%lsa) ; lstopt(1,i)=ia

      lstopt(0,i)=0
      If (ia > 0 .and. ia <= config%natms) Then !Tag
        If (config%lfrzn(ia) == 0) lstopt(0,i)=1
      End If

      ! There is no need to check for uncompressed unit since
      ! a tether is a point, it is either in or out by construction

      ! components of tether vector

      If (lstopt(0,i) > 0) Then
        xdab(i) = config%parts(ia)%xxx-stats%xin(ia)
        ydab(i) = config%parts(ia)%yyy-stats%yin(ia)
        zdab(i) = config%parts(ia)%zzz-stats%zin(ia)
      Else ! (DEBUG)
        xdab(i)=0.0_wp
        ydab(i)=0.0_wp
        zdab(i)=0.0_wp
      End If

    End Do

    ! periodic boundary condition

    Call images(config%imcon,config%cell,tether%ntteth,xdab,ydab,zdab)

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

    Do i=1,tether%ntteth
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

        kk=tether%listtet(0,i)

        ! calculate scalar constant terms

        If(tether%keytet(kk) == 1) Then

          ! harmonic function

          k=tether%prmtet(1,kk)

          omega=0.5_wp*k*rab**2
          gamma=k

        Else If (tether%keytet(kk) == 2) Then

          ! restrained harmonic

          k =tether%prmtet(1,kk)
          rc=tether%prmtet(2,kk)

          omega=0.5_wp*k*Min(rab,rc)**2 + k*rc*Sign(Max(rab-rc,0.0_wp),rab)
          gamma=k*Sign(Min(rab,rc),rab)*rrab

        Else If (tether%keytet(kk) == 3) Then

          ! quartic potential

          k2=tether%prmtet(1,kk)
          k3=tether%prmtet(2,kk)
          k4=tether%prmtet(3,kk)

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


        config%parts(ia)%fxx=config%parts(ia)%fxx+fx
        config%parts(ia)%fyy=config%parts(ia)%fyy+fy
        config%parts(ia)%fzz=config%parts(ia)%fzz+fz

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
