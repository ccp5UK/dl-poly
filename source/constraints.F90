Module constraints

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for defining global constraint bonds variables and
  ! arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2014
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp,wi
  Use comms,           Only : comms_type,gsum,gcheck,gsync
  Use configuration,   Only : configuration_type
  Use pmf,             Only : pmf_shake_vv,pmf_rattle,pmf_type
  Use constants,       Only : zero_plus
  Use errors_warnings, Only : error,warning,info
  Use shared_units,    Only : update_shared_units,SHARED_UNIT_UPDATE_POSITIONS
  Use numerics,        Only : images,local_index
  Use statistics,      Only : stats_type
  Use timer,           Only : start_timer, stop_timer, timer_type
  Use domains,         Only : domains_type
  Use thermostat,      Only : thermostat_type
  Implicit None

  Private

  Type,      Public :: constraints_type
    Private

    Logical,           Public              :: lshmv_con = .false.
    Integer,           Public              :: ntcons  = 0
    Integer,           Public              :: ntcons1 = 0
    Integer,           Public              :: m_con   = 0
    Integer,           Public              :: megcon
    Integer,           Public              :: mxtcon,mxcons,mxfcon
    Integer,           Public              :: max_iter_shake
    Real( Kind = wp ), Public              :: tolerance
    Integer,           Allocatable, Public :: numcon(:)
    Integer,           Allocatable, Public :: lstcon(:,:),listcon(:,:),legcon(:,:)
    Integer,           Allocatable, Public :: lishp_con(:),lashp_con(:)
    Integer,           Allocatable         :: lstopt(:,:),listot(:)
    Real( Kind = wp ), Allocatable         :: dxx(:),dyy(:),dzz(:)
    Real( Kind = wp ), Allocatable, Public :: prmcon(:)

  Contains
    Private

    Procedure, Public :: init => allocate_constraints_arrays
    Procedure, Public :: deallocate_constraints_temps
    Procedure, Public :: allocate_work
    Procedure, Public :: deallocate_work
    Final             :: deallocate_constraints_arrays
 
  End Type

  Public :: constraints_pseudo_bonds
  Public :: constraints_quench
  Public :: constraints_tags
  Public :: constraints_shake_vv
  Public :: constraints_rattle
  Public :: apply_shake
  Public :: apply_rattle

Contains

  Subroutine allocate_work(T,n)
  Class(constraints_type) :: T
    Integer, Intent( In ) :: n
    Integer               :: fail(2)
    Character(Len=100)    :: message

    If (T%megcon > 0) Then
      Allocate (T%lstopt(0:2,1:T%mxcons),T%listot(1:n),          Stat=fail( 1)) 
      Allocate (T%dxx(1:T%mxcons),T%dyy(1:T%mxcons),T%dzz(1:T%mxcons),      Stat=fail( 2))

    End If
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'failed to allocate work arrays for constraints'
      Call error(0,message)
    End If
  End Subroutine allocate_work

  Subroutine deallocate_work(T)
  Class(constraints_type) :: T

    Integer            :: fail(2)
    Character(Len=100) :: message

    If (T%megcon > 0) Then
      Deallocate (T%lstopt,T%listot,  Stat=fail( 1))
      Deallocate (T%dxx,T%dyy,T%dzz,    Stat=fail( 2))
    End If
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'failed to deallocate work arrays for constraints'
      Call error(0,message)
    End If
  End Subroutine deallocate_work

  Subroutine allocate_constraints_arrays(T,mxtmls,mxatdm,mxlshp,neighbours)
  Class(constraints_type) :: T
    Integer(kind=wi), Intent( In ) :: mxtmls,mxatdm,neighbours,mxlshp

    Integer :: fail(7)

    fail = 0

    Allocate (T%numcon(1:mxtmls),                        Stat = fail(1))
    Allocate (T%lstcon(1:2,1:T%mxtcon),                    Stat = fail(2))
    Allocate (T%listcon(0:2,1:T%mxcons),                   Stat = fail(3))
    Allocate (T%legcon(0:T%mxfcon,1:mxatdm),               Stat = fail(4))
    Allocate (T%lishp_con(1:mxlshp), Stat = fail(5))
    Allocate (T%lashp_con(1:neighbours), Stat = fail(6))
    Allocate (T%prmcon(1:T%mxtcon),                        Stat = fail(7))

    If (Any(fail > 0)) Call error(1018)

    T%numcon  = 0
    T%lstcon  = 0
    T%listcon = 0
    T%legcon  = 0

    T%lishp_con = 0 ; T%lashp_con = 0

    T%prmcon  = 0.0_wp

  End Subroutine allocate_constraints_arrays

  Subroutine deallocate_constraints_temps(T)
  Class(constraints_type) :: T
    Integer :: fail(2)
    If (Allocated(T%numcon)) Deallocate (T%numcon, Stat = fail(1))
    If (Allocated(T%lstcon)) Deallocate (T%lstcon, Stat = fail(2))
    If (Any(fail > 0)) Call error(1032)

  End Subroutine deallocate_constraints_temps
  Subroutine deallocate_constraints_arrays(T)
    Type(constraints_type) :: T

    Integer :: fail(7)

    fail = 0
    If (Allocated(T%numcon)) Deallocate (T%numcon, Stat = fail(1))
    If (Allocated(T%lstcon)) Deallocate (T%lstcon, Stat = fail(2))
    If (Allocated(T%numcon)) Deallocate (T%numcon, Stat = fail(1))
    If (Allocated(T%lstcon)) Deallocate (T%lstcon, Stat = fail(2))
    If (Allocated(T%listcon)) Deallocate (T%listcon, Stat = fail(3))
    If (Allocated(T%legcon)) Deallocate (T%legcon, Stat = fail(4))
    If (Allocated(T%lishp_con)) Deallocate (T%lishp_con, Stat = fail(5))
    If (Allocated(T%lashp_con)) Deallocate (T%lashp_con, Stat = fail(6))
    If (Allocated(T%prmcon)) Deallocate (T%prmcon, Stat = fail(7))

    If (Any(fail > 0)) Call error(1032)

  End Subroutine deallocate_constraints_arrays

  Subroutine constraints_pseudo_bonds(gxx,gyy,gzz,stat,cons,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for treating constraint bonds as stiff harmonic
    ! springs for use with the conjugate gradient method (minimise_relax.f90)
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith & i.t.todorov november 2011
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    Real( Kind = wp ),          Intent( InOut ) :: gxx(:),gyy(:),gzz(:)
    Type( stats_type),          Intent( InOut ) :: stat
    Type( constraints_type),    Intent( InOut ) :: cons
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ),         Intent( InOut ) :: comm

    Real( Kind = wp ), Parameter :: rigid=1.0e6_wp

    Integer           :: k,i,j
    Real( Kind = wp ) :: r,r0,ebond,gamma

    stat%engcon=0.0_wp
    Do k=1,cons%ntcons
      If (cons%lstopt(0,k) == 0) Then
        i=cons%lstopt(1,k)
        j=cons%lstopt(2,k)

        ! if a pair is frozen and constraint bonded, it is more frozen
        ! than constrained (users!!!)

        r=Sqrt(cons%dxx(k)**2+cons%dyy(k)**2+cons%dzz(k)**2)
        r0=cons%prmcon(cons%listcon(0,k))

        gamma=rigid*(r-r0)
        ebond=gamma*0.5_wp*(r-r0)
        gamma=gamma/r

        ! Accumulate energy and add forces

        If (i <= config%natms) Then
          stat%engcon=stat%engcon+ebond

          If (config%lfrzn(i) == 0) Then
            gxx(i)=gxx(i)-cons%dxx(k)*gamma
            gyy(i)=gyy(i)-cons%dyy(k)*gamma
            gzz(i)=gzz(i)-cons%dzz(k)*gamma
          End If
        End If

        If (j <= config%natms .and. config%lfrzn(j) == 0) Then
          gxx(j)=gxx(j)+cons%dxx(k)*gamma
          gyy(j)=gyy(j)+cons%dyy(k)*gamma
          gzz(j)=gzz(j)+cons%dzz(k)*gamma
        End If
      End If
    End Do

    ! global sum of energy

    Call gsum(comm,stat%engcon)

  End Subroutine constraints_pseudo_bonds

  Subroutine constraints_quench(cons,stat,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for quenching the internal bond energies in the
    ! initial structure of a molecule defined by constraints
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( constraints_type ),   Intent( InOut ) :: cons
    Type( stats_type),          Intent( InOut ) :: stat
    Type( domains_type ),       Intent( In    ) :: domain
    Type( comms_type),          Intent( InOut ) :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Logical                :: safe
    Integer                :: fail(1:2),i,j,k,icyc
    Real( Kind = wp )      :: dis,amti,amtj,dlj,dli,esig,gamma,gammi,gammj
    Character( Len = 256 ) :: message

    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)

    fail=0
    Allocate (lstitr(1:config%mxatms),                          Stat=fail(1))
    Allocate (vxt(1:config%mxatms),vyt(1:config%mxatms),vzt(1:config%mxatms), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'constraints_quench allocation failure'
      Call error(0,message)
    End If
    Call cons%allocate_work(config%mxatms)

    ! gather velocities of shared atoms

    If (cons%lshmv_con) Then
      Call update_shared_units(config,cons%lishp_con, &
        cons%lashp_con,config%vxx,config%vyy,config%vzz,domain,comm)
    End If

    ! construct current constrained bond vectors and listot array (shared
    ! constraint atoms) for iterative (constraints) algorithms

    lstitr(1:config%natms)=.false. ! initialise lstitr
    Call constraints_tags(lstitr,cons,config,comm)

    ! normalise constraint vectors

    Do k=1,cons%ntcons
      If (cons%lstopt(0,k) == 0) Then
        dis=1.0_wp/Sqrt(cons%dxx(k)**2+cons%dyy(k)**2+cons%dzz(k)**2)
        cons%dxx(k)=cons%dxx(k)*dis
        cons%dyy(k)=cons%dyy(k)*dis
        cons%dzz(k)=cons%dzz(k)*dis
      Else ! DEBUG
        cons%dxx(k)=0.0_wp
        cons%dyy(k)=0.0_wp
        cons%dzz(k)=0.0_wp
      End If
    End Do

    ! application of constraint (quench) algorithm
    ! Initialise number of cycles to zero and unsafe passage of the algorithm

    icyc=0
    safe=.false.

    Do While ((.not.safe) .and. icyc < cons%max_iter_shake)
      icyc=icyc+1

      ! initialise velocity correction arrays

      Do i=1,config%natms
        vxt(i)=0.0_wp
        vyt(i)=0.0_wp
        vzt(i)=0.0_wp
      End Do

      ! calculate velocity constraint corrections

      esig=0.0_wp
      Do k=1,cons%ntcons
        If (cons%lstopt(0,k) == 0) Then
          i=cons%lstopt(1,k)
          j=cons%lstopt(2,k)

          ! if a pair is frozen and constraint bonded, it is more frozen
          ! than constrained (users!!!)

          amti=1.0_wp/config%weight(i)
          amtj=1.0_wp/config%weight(j)

          ! no corrections for frozen atoms

          If (config%lfrzn(i) /= 0) amti=0.0_wp
          If (config%lfrzn(j) /= 0) amtj=0.0_wp

          ! calculate constraint force parameter - gamma

          gamma = cons%dxx(k)*(config%vxx(i)-config%vxx(j)) + cons%dyy(k)*&
            (config%vyy(i)-config%vyy(j)) + cons%dzz(k)*(config%vzz(i)-config%vzz(j))

          esig=Max(esig,0.5_wp*Abs(gamma))

          gamma = gamma / (amti+amtj)

          ! improve approximate constraint velocity

          If (i <= config%natms .and. config%lfrzn(i) == 0) Then
            gammi =-gamma*amti
            vxt(i)=vxt(i)+cons%dxx(k)*gammi
            vyt(i)=vyt(i)+cons%dyy(k)*gammi
            vzt(i)=vzt(i)+cons%dzz(k)*gammi
          End If

          If (j <= config%natms .and. config%lfrzn(j) == 0) Then
            gammj = gamma*amtj
            vxt(j)=vxt(j)+cons%dxx(k)*gammj
            vyt(j)=vyt(j)+cons%dyy(k)*gammj
            vzt(j)=vzt(j)+cons%dzz(k)*gammj
          End If
        End If
      End Do

      ! global verification of convergence

      safe=(esig < cons%tolerance)
      Call gcheck(comm,safe,"enforce")

      ! bypass next section and terminate iteration if all tolerances ok

      If (.not.safe) Then

        ! update velocities

        Do k=1,cons%ntcons
          If (cons%lstopt(0,k) == 0) Then
            i=cons%lstopt(1,k)
            j=cons%lstopt(2,k)

            If (i <= config%natms .and. config%lfrzn(i) == 0) Then
              dli = 1.0_wp/Real(cons%listot(i),wp)
              config%vxx(i)=config%vxx(i)+vxt(i)*dli
              config%vyy(i)=config%vyy(i)+vyt(i)*dli
              config%vzz(i)=config%vzz(i)+vzt(i)*dli
            End If

            If (j <= config%natms .and. config%lfrzn(j) == 0) Then
              dlj = 1.0_wp/Real(cons%listot(j),wp)
              config%vxx(j)=config%vxx(j)+vxt(j)*dlj
              config%vyy(j)=config%vyy(j)+vyt(j)*dlj
              config%vzz(j)=config%vzz(j)+vzt(j)*dlj
            End If
          End If
        End Do

        ! transport velocity updates to other nodes
        If (cons%lshmv_con) Then
          Call update_shared_units(config,cons%lishp_con, &
            cons%lashp_con,config%vxx,config%vyy,config%vzz,domain,comm)
        End If
      End If
    End Do

    If (.not.safe) Then ! error exit if quenching fails
      Call error(70)
    Else ! Collect per call passage statistics
      stat%passcnq(1)=Real(icyc-1,wp)
      stat%passcnq(3)=stat%passcnq(2)*stat%passcnq(3)
      stat%passcnq(2)=stat%passcnq(2)+1.0_wp
      stat%passcnq(3)=stat%passcnq(3)/stat%passcnq(2)+stat%passcnq(1)/stat%passcnq(2)
      stat%passcnq(4)=Min(stat%passcnq(1),stat%passcnq(4))
      stat%passcnq(5)=Max(stat%passcnq(1),stat%passcnq(5))
      stat%passcnq(1)=0.0_wp ! Reset
    End If

    Deallocate (lstitr,        Stat=fail(1))
    Deallocate (vxt,vyt,vzt,   Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'constraints_quench deallocation failure'
      Call error(0,message)
    End If
    Call cons%deallocate_work()
  End Subroutine constraints_quench

  Subroutine constraints_tags(lstitr,cons,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for identifying, constructing and indexing
    ! constraints' vectors for iterative (constraints) algorithms
    !
    ! Note: must be used in conjunction with integration algorithms
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,                    Intent( InOut ) :: lstitr(:)
    Type( constraints_type ),   Intent( InOut ) :: cons
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ),         Intent( InOut ) :: comm

    Logical                :: safe
    Integer                :: fail,i,j,k,l
    Character( Len = 256 ) :: message

    Logical, Allocatable :: lunsafe(:)


    fail=0
    Allocate (lunsafe(1:cons%mxcons), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'constraints_tags allocation failure'
      Call error(0,message)
    End If

    ! initialise cons%listot array (shared constraint bond)

    Do k=1,config%natms
      cons%listot(k)=0
    End Do

    Do k=1,cons%ntcons
      lunsafe(k)=.false.

      ! Opt out bond from action by default

      cons%lstopt(0,k)=1

      ! indices of atoms in bond

      i=local_index(cons%listcon(1,k),config%nlast,config%lsi,config%lsa)
      j=local_index(cons%listcon(2,k),config%nlast,config%lsi,config%lsa)

      ! store indices

      cons%lstopt(1,k)=i
      cons%lstopt(2,k)=j

      ! for all native and natively shared constraints

      If ((i > 0 .and. j > 0) .and. (i <= config%natms .or. j <= config%natms)) Then

        ! if a pair is frozen and constraint bonded
        ! it is more frozen than constrained (users!!!)

        If (config%lfrzn(i)*config%lfrzn(j) == 0) Then

          ! Select bond for action if product gives zero

          cons%lstopt(0,k)=0

          ! calculate bond vectors

          cons%dxx(k)=config%parts(i)%xxx-config%parts(j)%xxx
          cons%dyy(k)=config%parts(i)%yyy-config%parts(j)%yyy
          cons%dzz(k)=config%parts(i)%zzz-config%parts(j)%zzz

          ! indicate sharing on local ends of bonds
          ! summed contributions (quench/shake/rattle) for each local
          ! constrained bonded atom must be config%weighted by cons%listot(atom) -
          ! how many bond constraints an atom is involved in

          If (i <= config%natms) cons%listot(i)=cons%listot(i)+1
          If (j <= config%natms) cons%listot(j)=cons%listot(j)+1

        End If

      Else If ((i == 0 .and. j == 0) .or. (i > config%natms .and. j > config%natms) .or. &
        (i == 0 .and. j > config%natms) .or. (j == 0 .and. i > config%natms)) Then

        ! constraints lying outside or completely in the halo, or partly in the
        ! halo and partly outside it are not considered and zero bond vectors
        ! are assigned (DEBUG)

        cons%dxx(k)=0.0_wp
        cons%dyy(k)=0.0_wp
        cons%dzz(k)=0.0_wp

      Else

        ! weird case - constraints flying apart
        ! Detect uncompressed unit

        lunsafe(k)=.true.

      End If
    End Do

    ! check for incomplete bond data

    safe = .not. Any(lunsafe(1:cons%ntcons))
    Call gcheck(comm,safe)
    If (.not.safe) Then
      Do l=0,comm%mxnode-1
        If (comm%idnode == l) Then
          Do k=1,cons%ntcons
            If (lunsafe(k)) Then
              Write(message,'(2(a,i10))')     &
                'global unit number', cons%listcon(0,k), &
                ' , with a head particle number', cons%listcon(1,k)
              Call info(message)
              Call warning('contributes towards next error')
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(111)
    End If

    ! minimum image convention for bond vectors

    Call images(config%imcon,config%cell,cons%ntcons,cons%dxx,cons%dyy,cons%dzz)

    ! update lstitr

    Do k=1,config%natms
      lstitr(k)=(cons%listot(k) > 0 .and. config%lfrzn(k) == 0)
    End Do

    Deallocate (lunsafe, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'constraints_tags deallocation failure'
      Call error(0,message)
    End If

  End Subroutine constraints_tags

  Subroutine constraints_rattle(tstep,lfst,lcol,config,stat,cons,domain,tmr,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for applying constraint corrections to the
    ! velocities of constrained atoms
    !
    ! Note: must be used in conjunction with integration algorithms
    !       VV applicable ONLY
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),          Intent( In    ) :: tstep
    Logical,                    Intent( In    ) :: lfst,lcol
    Type( configuration_type ), Intent( InOut ) :: config
    Type( stats_type ),         Intent( InOut ) :: stat
    Type( constraints_type ),   Intent( InOut ) :: cons
    Type( domains_type ),       Intent( In    ) :: domain
    Type( timer_type ),         Intent( InOut ) :: tmr
    Type( comms_type ),         Intent( InOut ) :: comm

    Logical                :: safe
    Integer                :: fail,i,j,k,icyc
    Real( Kind = wp )      :: dis,amti,amtj,dli,dlj,esig,gamma,gammi,gammj
    Character( Len = 256 ) :: message

    Real( Kind = wp ), Dimension( : ), Allocatable :: vxt,vyt,vzt


#ifdef CHRONO
    Call start_timer(tmr, 'Rattle Constraints')
#endif
    fail=0
    Allocate (vxt(1:config%mxatms),vyt(1:config%mxatms),vzt(1:config%mxatms), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'constraints_rattle allocation failure'
      Call error(0,message)
    End If

    ! normalise constraint vectors on first pass outside

    If (lfst) Then
      Do k=1,cons%ntcons
        If (cons%lstopt(0,k) == 0) Then
          dis=1.0_wp/Sqrt(cons%dxx(k)**2+cons%dyy(k)**2+cons%dzz(k)**2)
          cons%dxx(k)=cons%dxx(k)*dis
          cons%dyy(k)=cons%dyy(k)*dis
          cons%dzz(k)=cons%dzz(k)*dis
        Else ! DEBUG
          cons%dxx(k)=0.0_wp
          cons%dyy(k)=0.0_wp
          cons%dzz(k)=0.0_wp
        End If
      End Do
    End If

    ! application of constraint (rattle) algorithm
    ! Initialise number of cycles to zero and unsafe passage of the algorithm

    safe=.false.
    icyc=0
    Do While ((.not.safe) .and. icyc < cons%max_iter_shake)
      icyc=icyc+1

      ! update velocities globally: transport velocity updates of shared atoms to other nodes

      If (cons%lshmv_con) Then
        Call update_shared_units(config,cons%lishp_con, &
          cons%lashp_con,config%vxx,config%vyy,config%vzz,domain,comm)
      End If

      ! initialise velocity correction arrays

      Do i=1,config%natms
        vxt(i)=0.0_wp
        vyt(i)=0.0_wp
        vzt(i)=0.0_wp
      End Do

      ! calculate velocity constraint corrections

      esig=0.0_wp
      Do k=1,cons%ntcons
        If (cons%lstopt(0,k) == 0) Then
          i=cons%lstopt(1,k)
          j=cons%lstopt(2,k)

          amti=tstep/config%weight(i)
          amtj=tstep/config%weight(j)

          ! no corrections for frozen atoms

          If (config%lfrzn(i) /= 0) amti=0.0_wp
          If (config%lfrzn(j) /= 0) amtj=0.0_wp

          ! calculate constraint force parameter - gamma

          gamma = cons%dxx(k)*(config%vxx(i)-config%vxx(j)) + cons%dyy(k)&
            *(config%vyy(i)-config%vyy(j)) + cons%dzz(k)*(config%vzz(i)-config%vzz(j))

          esig=Max(esig,0.5_wp*tstep*Abs(gamma))

          gamma = gamma / (amti+amtj)

          ! improve approximate constraint velocity and force

          If (i <= config%natms .and. config%lfrzn(i) == 0) Then
            gammi =-gamma*amti
            vxt(i)=vxt(i)+cons%dxx(k)*gammi
            vyt(i)=vyt(i)+cons%dyy(k)*gammi
            vzt(i)=vzt(i)+cons%dzz(k)*gammi
          End If

          If (j <= config%natms .and. config%lfrzn(j) == 0) Then
            gammj = gamma*amtj
            vxt(j)=vxt(j)+cons%dxx(k)*gammj
            vyt(j)=vyt(j)+cons%dyy(k)*gammj
            vzt(j)=vzt(j)+cons%dzz(k)*gammj
          End If
        End If
      End Do

      ! global verification of convergence

      safe=(esig < cons%tolerance)
      Call gcheck(comm,safe,"enforce")

      ! bypass next section and terminate iteration if all tolerances ok

      If (.not.safe) Then

        ! update velocities locally

        Do k=1,cons%ntcons
          If (cons%lstopt(0,k) == 0) Then
            i=cons%lstopt(1,k)
            j=cons%lstopt(2,k)

            If (i <= config%natms .and. config%lfrzn(i) == 0) Then
              dli = 1.0_wp/Real(cons%listot(i),wp)
              config%vxx(i)=config%vxx(i)+vxt(i)*dli
              config%vyy(i)=config%vyy(i)+vyt(i)*dli
              config%vzz(i)=config%vzz(i)+vzt(i)*dli
            End If

            If (j <= config%natms .and. config%lfrzn(j) == 0) Then
              dlj = 1.0_wp/Real(cons%listot(j),wp)
              config%vxx(j)=config%vxx(j)+vxt(j)*dlj
              config%vyy(j)=config%vyy(j)+vyt(j)*dlj
              config%vzz(j)=config%vzz(j)+vzt(j)*dlj
            End If
          End If
        End Do

      End If
    End Do

    If (.not.safe) Then ! error exit for non-convergence
      Call error(515)
    Else ! Collect per call and per step passage statistics
      stat%passcon(1,1,2)=Real(icyc-1,wp)
      stat%passcon(3,1,2)=stat%passcon(2,1,2)*stat%passcon(3,1,2)
      stat%passcon(2,1,2)=stat%passcon(2,1,2)+1.0_wp
      stat%passcon(3,1,2)=stat%passcon(3,1,2)/stat%passcon(2,1,2)+stat%passcon(1,1,2)/stat%passcon(2,1,2)
      stat%passcon(4,1,2)=Min(stat%passcon(1,1,2),stat%passcon(4,1,2))
      stat%passcon(5,1,2)=Max(stat%passcon(1,1,2),stat%passcon(5,1,2))

      stat%passcon(1,2,2)=stat%passcon(1,2,2)+stat%passcon(1,1,2)
      If (lcol) Then ! Collect
        stat%passcon(3,2,2)=stat%passcon(2,2,2)*stat%passcon(3,2,2)
        stat%passcon(2,2,2)=stat%passcon(2,2,2)+1.0_wp
        stat%passcon(3,2,2)=stat%passcon(3,2,2)/stat%passcon(2,2,2)+stat%passcon(1,2,2)/stat%passcon(2,2,2)
        stat%passcon(4,2,2)=Min(stat%passcon(1,2,2),stat%passcon(4,2,2))
        stat%passcon(5,2,2)=Max(stat%passcon(1,2,2),stat%passcon(5,2,2))
        stat%passcon(1,2,2)=0.0_wp ! Reset
      End If
      stat%passcon(1,1,2)=0.0_wp ! Reset
    End If

    Deallocate (vxt,vyt,vzt, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'constraints_rattle deallocation failure'
      Call error(0,message)
    End If
#ifdef CHRONO
    Call stop_timer(tmr, 'Rattle Constraints')
#endif

  End Subroutine constraints_rattle


  Subroutine constraints_shake_vv(tstep,config,str,vir,stat,cons,domain,tmr,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for applying bond constraint corrections after
    ! unconstrained motion
    !
    ! Note: must be used in conjunction with integration algorithms
    !       VV compliant
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Real( Kind = wp ),          Intent( In    ) :: tstep
    Type( configuration_type ), Intent( InOut ) :: config
    Real( Kind = wp ),          Intent( InOut ) :: vir, str(:)
    Type( constraints_type),    Intent( InOut ) :: cons
    Type( stats_type),          Intent( InOut ) :: stat
    Type( domains_type ),       Intent( In    ) :: domain
    Type( timer_type ),         Intent( InOut ) :: tmr
    TYpe( comms_type ),         Intent( InOut ) :: comm

    Logical           :: safe
    Integer           :: fail(1:2),i,j,k,icyc
    Real( Kind = wp ) :: amti,amtj,dli,dlj,gamma,gammi,gammj,tstep2

    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt
    Real( Kind = wp ), Dimension( : ), Allocatable :: dxt,dyt,dzt,dt2,esig

    Character(Len=256) :: message
#ifdef CHRONO
    Call start_timer(tmr, 'SHAKE Constraints')
#endif
    fail=0
    Allocate (xxt(1:config%mxatms),yyt(1:config%mxatms),zzt(1:config%mxatms),                              Stat=fail(1))
    Allocate (dxt(1:cons%mxcons),dyt(1:cons%mxcons),dzt(1:cons%mxcons),dt2(1:cons%mxcons),esig(1:cons%mxcons), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'constraints_shake allocation failure'
      Call error(0,message)
    End If


    ! Initialise constraint virial and stress

    vir=0.0_wp
    str=0.0_wp

    ! squared timestep

    tstep2 = tstep*tstep

    ! application of constraint (shake) algorithm
    ! start bond vectors are cons%dxx = xxx(i) - xxx(j) etc.
    ! Initialise number of cycles to zero and unsafe passage of the algorithm

    safe=.false.
    icyc=0
    Do While ((.not.safe) .and. icyc < cons%max_iter_shake)
      icyc=icyc+1

      ! update positions globally: transport position updates of shared atoms to other nodes

      If (cons%lshmv_con) Then
        Call update_shared_units(config,cons%lishp_con, &
          cons%lashp_con,SHARED_UNIT_UPDATE_POSITIONS,domain,comm)
      End If

      ! calculate temporary bond vector

      Do k=1,cons%ntcons
        If (cons%lstopt(0,k) == 0) Then
          i=cons%lstopt(1,k)
          j=cons%lstopt(2,k)

          dxt(k)=config%parts(i)%xxx-config%parts(j)%xxx
          dyt(k)=config%parts(i)%yyy-config%parts(j)%yyy
          dzt(k)=config%parts(i)%zzz-config%parts(j)%zzz
        Else ! DEBUG
          dxt(k)=0.0_wp
          dyt(k)=0.0_wp
          dzt(k)=0.0_wp
        End If
      End Do

      ! periodic boundary condition

      Call images(config%imcon,config%cell,cons%ntcons,dxt,dyt,dzt)

      ! calculate maximum error in bondlength and
      ! do a global verification of convergence

      safe=.true.
      Do k=1,cons%ntcons
        If (cons%lstopt(0,k) == 0) Then
          dt2(k) =dxt(k)**2+dyt(k)**2+dzt(k)**2 - cons%prmcon(cons%listcon(0,k))**2
          esig(k)=0.5_wp*Abs(dt2(k))
          safe=(safe .and. (esig(k) < cons%tolerance*cons%prmcon(cons%listcon(0,k))))
        Else
          dt2(k) =0.0_wp
          esig(k)=0.0_wp
        End If
      End Do
      Call gcheck(comm,safe,"enforce")

      ! bypass next section and terminate iteration if all tolerances ok

      If (.not.safe) Then

        ! initialise position correction arrays

        Do i=1,config%natms
          xxt(i)=0.0_wp
          yyt(i)=0.0_wp
          zzt(i)=0.0_wp
        End Do

        ! calculate constraint forces

        Do k=1,cons%ntcons
          If (cons%lstopt(0,k) == 0) Then
            i=cons%lstopt(1,k)
            j=cons%lstopt(2,k)

            amti=tstep2/config%weight(i)
            amtj=tstep2/config%weight(j)

            ! no corrections for frozen atoms

            If (config%lfrzn(i) /= 0) amti=0.0_wp
            If (config%lfrzn(j) /= 0) amtj=0.0_wp

            ! calculate constraint force parameter

            gamma = dt2(k) / ((amti+amtj)*(cons%dxx(k)*dxt(k)+cons%dyy(k)*dyt(k)+cons%dzz(k)*dzt(k)))

            If (i <= config%natms) Then

              ! accumulate bond stress

              str(1) =str(1) - gamma*cons%dxx(k)*cons%dxx(k)
              str(2) =str(2) - gamma*cons%dxx(k)*cons%dyy(k)
              str(3) =str(3) - gamma*cons%dxx(k)*cons%dzz(k)
              str(5) =str(5) - gamma*cons%dyy(k)*cons%dyy(k)
              str(6) =str(6) - gamma*cons%dyy(k)*cons%dzz(k)
              str(9) =str(9) - gamma*cons%dzz(k)*cons%dzz(k)

              ! calculate atomic position constraint corrections

              If (config%lfrzn(i) == 0) Then
                gammi =-0.5_wp*gamma*amti
                xxt(i)=xxt(i)+cons%dxx(k)*gammi
                yyt(i)=yyt(i)+cons%dyy(k)*gammi
                zzt(i)=zzt(i)+cons%dzz(k)*gammi
              End If

            End If

            If (j <= config%natms .and. config%lfrzn(j) == 0) Then
              gammj = 0.5_wp*gamma*amtj
              xxt(j)=xxt(j)+cons%dxx(k)*gammj
              yyt(j)=yyt(j)+cons%dyy(k)*gammj
              zzt(j)=zzt(j)+cons%dzz(k)*gammj
            End If
          End If
        End Do

        ! update positions locally

        Do k=1,cons%ntcons
          If (cons%lstopt(0,k) == 0) Then
            i=cons%lstopt(1,k)
            j=cons%lstopt(2,k)

            ! apply position corrections if non-frozen

            If (i <= config%natms .and. config%lfrzn(i) == 0) Then
              dli = 1.0_wp/Real(cons%listot(i),wp)
              config%parts(i)%xxx=config%parts(i)%xxx+xxt(i)*dli
              config%parts(i)%yyy=config%parts(i)%yyy+yyt(i)*dli
              config%parts(i)%zzz=config%parts(i)%zzz+zzt(i)*dli
            End If

            If (j <= config%natms .and. config%lfrzn(j) == 0) Then
              dlj = 1.0_wp/Real(cons%listot(j),wp)
              config%parts(j)%xxx=config%parts(j)%xxx+xxt(j)*dlj
              config%parts(j)%yyy=config%parts(j)%yyy+yyt(j)*dlj
              config%parts(j)%zzz=config%parts(j)%zzz+zzt(j)*dlj
            End If
          End If
        End Do

      End If
    End Do

    If (.not.safe) Then ! error exit for non-convergence
      Do i=0,comm%mxnode-1
        If (comm%idnode == i) Then
          Do k=1,cons%ntcons
            If (esig(k) >= cons%tolerance*cons%prmcon(cons%listcon(0,k))) Then
              Write(message,'(3(a,i10))') &
                'global constraint number', cons%listcon(0,k),  &
                ' , with particle numbers:', cons%listcon(1,k), &
                ' &', cons%listcon(2,k)
              Call info(message)
              Write(message,'(a,f8.2,a,1p,e12.4)') &
                'converges to a length of ', &
                Sqrt(dt2(k)+cons%prmcon(cons%listcon(0,k))**2), &
                ' Angstroms with a factor', esig(k)/cons%prmcon(cons%listcon(0,k))
              Call info(message)
              Write(message,'(a)') 'contributes towards next error'
              Call warning(message)
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(105)
    Else ! Collect per call and per step passage statistics
      stat%passcon(1,1,1)=Real(icyc-1,wp)
      stat%passcon(3,1,1)=stat%passcon(2,1,1)*stat%passcon(3,1,1)
      stat%passcon(2,1,1)=stat%passcon(2,1,1)+1.0_wp
      stat%passcon(3,1,1)=stat%passcon(3,1,1)/stat%passcon(2,1,1)+stat%passcon(1,1,1)/stat%passcon(2,1,1)
      stat%passcon(4,1,1)=Min(stat%passcon(1,1,1),stat%passcon(4,1,1))
      stat%passcon(5,1,1)=Max(stat%passcon(1,1,1),stat%passcon(5,1,1))

      stat%passcon(1,2,1)=stat%passcon(1,2,1)+stat%passcon(1,1,1)
      stat%passcon(1,1,1)=0.0_wp ! Reset
    End If

    ! global sum of stress tensor

    Call gsum(comm,str)

    ! complete stress tensor (symmetrise)

    str(4) = str(2)
    str(7) = str(3)
    str(8) = str(6)

    ! total constraint virial

    vir=-(str(1)+str(5)+str(9))

    Deallocate (xxt,yyt,zzt,          Stat=fail(1))
    Deallocate (dxt,dyt,dzt,dt2,esig, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'constraints_shake deallocation failure'
      Call error(0,message)
    End If

#ifdef CHRONO
    Call stop_timer(tmr, 'SHAKE Constraints')
#endif
  End Subroutine constraints_shake_vv

  Subroutine apply_rattle(tstep,kit,pmf,cons,stat,domain,tmr,config,comm)

    Integer,                    Intent( In    ) :: kit
    Type( configuration_type ), Intent( InOut ) :: config
    Real( Kind = wp ),          Intent( In    ) :: tstep
    Type( stats_type),          Intent( InOut ) :: stat
    Type( pmf_type ),           Intent( InOut ) :: pmf
    Type( constraints_type),    Intent( InOut ) :: cons
    Type( domains_type ),       Intent( In    ) :: domain
    Type( timer_type ),         Intent( InOut ) :: tmr
    Type( comms_type ),         Intent( InOut ) :: comm

    Logical :: lfst, lcol
    Integer :: i

    Do i=1,kit
      lfst = (i == 1)
      lcol = (i == kit)

      If (cons%megcon > 0) Then
        Call constraints_rattle(tstep,lfst,lcol,config,stat,cons,domain,tmr,comm)
      End IF

      If (pmf%megpmf > 0) Then
        Call pmf_rattle(cons%max_iter_shake,cons%tolerance,tstep,lfst,lcol, &
          config,stat,pmf,comm)
      End If
    End Do
  End Subroutine apply_rattle

  Subroutine apply_shake(tstep,oxt,oyt,ozt,lstitr,stat,pmf,cons, &
      domain,tmr,config,thermo,comm)
 
    Logical,                   Intent( In    ) :: lstitr(:)
    Real( Kind = wp ),         Intent( InOut ) :: oxt(:),oyt(:),ozt(:)
    Real( Kind = wp ),         Intent( In    ) :: tstep
    Type( stats_type),         Intent( InOut ) :: stat
    Type( pmf_type ),          Intent( InOut ) :: pmf
    Type( constraints_type),   Intent( InOut ) :: cons
    Type( domains_type ),      Intent( In    ) :: domain
    Type( timer_type ),        Intent( InOut ) :: tmr
    Type( thermostat_type ),   Intent( InOut ) :: thermo
    Type( comms_type ),        Intent( InOut ) :: comm
    Type( configuration_type), Intent( InOut ) :: config
    ! constraint virial and stress tensor

    Logical           :: safe
    Integer(kind=wi)  :: i,j
    Real( Kind = wp ) :: xt,yt,zt,vir,str(1:9),hstep,rstep,tmp
    
    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep

    ! SHAKE procedures

    safe=.false.
    thermo%kit =0

    ! store integrated positions

    Do j=1,config%nfree
      i=config%lstfre(j)

      If (lstitr(i)) Then
        oxt(i)=config%parts(i)%xxx
        oyt(i)=config%parts(i)%yyy
        ozt(i)=config%parts(i)%zzz
      End If
    End Do

    Do While ((.not.safe) .and. thermo%kit <= thermo%mxkit)
      thermo%kit=thermo%kit+1

      If (cons%megcon > 0) Then

        ! apply constraint correction: stat%vircon,stat%strcon - constraint virial,stress

        Call constraints_shake_vv(tstep,config,str,vir,stat,cons,domain,tmr,comm)

        ! constraint virial and stress tensor

        stat%vircon=stat%vircon+vir
        stat%strcon=stat%strcon+str

        safe=.true.
      End If

      If (pmf%megpmf > 0) Then

        ! apply PMF correction: virpmf,strpmf - PMF constraint virial,stress

        Call pmf_shake_vv  &
          (cons%max_iter_shake,cons%tolerance,tstep, &
          config,str,vir,stat,pmf,comm)

        ! PMF virial and stress tensor

        stat%virpmf=stat%virpmf+vir
        stat%strpmf=stat%strpmf+str

        safe=(Abs(vir) <= zero_plus)
      End If
    End Do

    If (.not.safe) Call error(478)

    ! Collect per step passage statistics for bond and pmf constraints

    If (cons%megcon > 0) Then
      stat%passcon(3,2,1)=stat%passcon(2,2,1)*stat%passcon(3,2,1)
      stat%passcon(2,2,1)=stat%passcon(2,2,1)+1.0_wp
      stat%passcon(3,2,1)=stat%passcon(3,2,1)/stat%passcon(2,2,1)+stat%passcon(1,2,1)/stat%passcon(2,2,1)
      stat%passcon(4,2,1)=Min(stat%passcon(1,2,1),stat%passcon(4,2,1))
      stat%passcon(5,2,1)=Max(stat%passcon(1,2,1),stat%passcon(5,2,1))
      stat%passcon(1,2,1)=0.0_wp ! Reset
    End If

    If (pmf%megpmf > 0) Then
      stat%passpmf(3,2,1)=stat%passpmf(2,2,1)*stat%passpmf(3,2,1)
      stat%passpmf(2,2,1)=stat%passpmf(2,2,1)+1.0_wp
      stat%passpmf(3,2,1)=stat%passpmf(3,2,1)/stat%passpmf(2,2,1)+stat%passpmf(1,2,1)/stat%passpmf(2,2,1)
      stat%passpmf(4,2,1)=Min(stat%passpmf(1,2,1),stat%passpmf(4,2,1))
      stat%passpmf(5,2,1)=Max(stat%passpmf(1,2,1),stat%passpmf(5,2,1))
      stat%passpmf(1,2,1)=0.0_wp ! Reset
    End If

    ! calculate velocity and force correction

    Do j=1,config%nfree
      i=config%lstfre(j)

      If (lstitr(i)) Then
        xt=(config%parts(i)%xxx-oxt(i))*rstep
        yt=(config%parts(i)%yyy-oyt(i))*rstep
        zt=(config%parts(i)%zzz-ozt(i))*rstep

        config%vxx(i)=config%vxx(i)+xt
        config%vyy(i)=config%vyy(i)+yt
        config%vzz(i)=config%vzz(i)+zt

        tmp=config%weight(i)/hstep
        xt=xt*tmp
        yt=yt*tmp
        zt=zt*tmp

        config%parts(i)%fxx=config%parts(i)%fxx+xt
        config%parts(i)%fyy=config%parts(i)%fyy+yt
        config%parts(i)%fzz=config%parts(i)%fzz+zt
      End If
    End Do
  End Subroutine apply_shake
End module constraints
