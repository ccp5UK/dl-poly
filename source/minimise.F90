Module minimise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring minimisation routine arrays
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

  Use comms,           Only: comms_type,&
                             gmax,&
                             gsum
  Use configuration,   Only: configuration_type,&
                             getcom,&
                             write_config,&
                             IMCON_NOPBC,&
                             IMCON_CUBIC,&
                             IMCON_ORTHORHOMBIC,&
                             IMCON_PARALLELOPIPED,&
                             IMCON_SLAB,&
                             IMCON_TRUNC_OCTO,&
                             IMCON_RHOMBIC_DODEC,&
                             IMCON_HEXAGONAL
  Use constants,       Only: engunit,&
                             zero_plus
  Use constraints,     Only: constraints_pseudo_bonds,&
                             constraints_tags,&
                             constraints_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use filename,        Only: FILE_OUTPUT,&
                             file_type
  Use io,              Only: io_type
  Use kinds,           Only: wi,&
                             wp
  Use kinetics,        Only: getkin,&
                             getknf,&
                             getknr,&
                             getknt,&
                             getvom,&
                             kinstresf,&
                             kinstress,&
                             kinstrest
  Use netcdf_wrap,     Only: netcdf_param
  Use numerics,        Only: images,&
                             invert
  Use parse,           Only: lower_case,&
                             strip_blanks
  Use pmf,             Only: pmf_pseudo_bonds,&
                             pmf_tags,&
                             pmf_type
  Use rigid_bodies,    Only: getrotmat,&
                             q_setup,&
                             rigid_bodies_move,&
                             rigid_bodies_split_torque,&
                             rigid_bodies_type
  Use shared_units,    Only: SHARED_UNIT_UPDATE_POSITIONS,&
                             update_shared_units
  Use statistics,      Only: stats_type

  Implicit None

  Private

  ! Minimisation keys
  !> no minimisation
  Integer(Kind=wi), Parameter, Public :: MIN_NULL = -1
  !> Force minimisation
  Integer(Kind=wi), Parameter, Public :: MIN_FORCE = 1
  !> Energy minimisation
  Integer(Kind=wi), Parameter, Public :: MIN_ENERGY = 2
  !> Distance minimisation
  Integer(Kind=wi), Parameter, Public :: MIN_DISTANCE = 3

  Integer(Kind=wi), Parameter, Public :: NO_OPTIMISATION = 0
  Integer(Kind=wi), Parameter, Public :: LINE_SEARCH = 1
  Integer(Kind=wi), Parameter, Public :: LINEAR_EXTRAPOLATION = 2

  Type, Public :: minimise_type
    Private

    !> Minimisation switch
    Logical, Public                    :: minimise
    !> Minimisation key
    Integer(Kind=wi), Public           :: key
    !> Relaxed indicator?
    Logical, Public                    :: relaxed = .true.
    !> Transport switch
    Logical, Public                    :: transport = .false.
    !> Minimisation fequency (steps)
    Integer(Kind=wi), Public           :: freq
    !> Minimisation tolerance
    Real(Kind=wp), Public              :: tolerance
    !> Conjugate gradients method step length
    Real(Kind=wp), Public              :: step_length
    !> Coordinate arrays?
    Real(Kind=wp), Allocatable, Public :: oxx(:), oyy(:), ozz(:)
    Logical                            :: newjob = .true., l_rdf, l_mov
    Character(Len=8)                   :: word
    Integer                            :: keyopt
    Real(Kind=wp)                      :: min_pass
    Real(Kind=wp)                      :: total, grad_tol, eng_tol, dist_tol, step, &
                         eng_0, eng_min, eng, eng0, eng1, eng2, &
                         grad, grad0, grad1, grad2, onorm, sgn, stride, gamma
  Contains
    Private

    Procedure, Public :: init => allocate_minimise_arrays
    Final             :: deallocate_minimise_arrays
  End Type minimise_type

  Public :: minimise_relax, zero_k_optimise

Contains

  Subroutine allocate_minimise_arrays(T, mxatms)
    Class(minimise_type)            :: T
    Integer(Kind=wi), Intent(In   ) :: mxatms

    Integer :: fail

    fail = 0

    Allocate (T%oxx(1:mxatms), T%oyy(1:mxatms), T%ozz(1:mxatms), Stat=fail)

    If (fail > 0) Call error(1038)

    T%oxx = 0.0_wp
    T%oyy = 0.0_wp
    T%ozz = 0.0_wp
  End Subroutine allocate_minimise_arrays

  Subroutine deallocate_minimise_arrays(T)
    Type(minimise_type) :: T

    If (Allocated(T%oxx)) Then
      Deallocate (T%oxx)
    End If
    If (Allocated(T%oyy)) Then
      Deallocate (T%oyy)
    End If
    If (Allocated(T%ozz)) Then
      Deallocate (T%ozz)
    End If
  End Subroutine deallocate_minimise_arrays

  Subroutine minimise_relax(l_str, rdf_collect, tstep, stpcfg, io, stats, &
                            pmf, cons, netcdf, minim, rigid, domain, config, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for minimising the structure by using conjugate
    ! gradient method (CGM).
    !
    ! Note: minimisation type and criterion:
    !       minim%key = MIN_FORCE : absolute force
    !       minim%key = MIN_ENERGY : relative energy
    !       minim%key = MIN_DISTANCE : absolute displacement
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & w.smith february 2014
    ! contrib   - a.m.elena february 2017
    ! contrib   - i.t.todorov february 2017
    ! contrib   - i.scivetti april 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,                  Intent(In   ) :: l_str
    Logical,                  Intent(InOut) :: rdf_collect
    Real(Kind=wp),            Intent(In   ) :: tstep, stpcfg
    Type(io_type),            Intent(InOut) :: io
    Type(stats_type),         Intent(InOut) :: stats
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(constraints_type),   Intent(InOut) :: cons
    Type(netcdf_param),       Intent(In   ) :: netcdf
    Type(minimise_type),      Intent(InOut) :: minim
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm

    Integer, Parameter :: mxpass = 1000

    Character(Len=10)          :: c_out
    Character(Len=256)         :: message
    Integer                    :: fail(1:8), i, j, levcfg
    Logical                    :: l_out
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp), Allocatable :: gxx(:), gyy(:), gzz(:), txx(:), tyy(:), tzz(:), uxx(:), uyy(:), &
                                  uzz(:)
    Type(file_type)            :: minfile

! OUTPUT existence
! Optimisation iteration and convergence limits
! Constraints and PMFs arrays

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms), Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
    End If
    If (rigid%total > 0) Then
      Allocate (txx(1:config%mxatms), tyy(1:config%mxatms), tzz(1:config%mxatms), Stat=fail(6))
      Allocate (uxx(1:config%mxatms), uyy(1:config%mxatms), uzz(1:config%mxatms), Stat=fail(7))
    End If
    Allocate (gxx(1:config%mxatms), gyy(1:config%mxatms), gzz(1:config%mxatms), Stat=fail(8))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'minimise_relax allocation failure'
      Call error(0, message)
    End If

    If (minim%newjob) Then
      minim%newjob = .false.

      ! At start no optimisation has been attempted yet

      minim%keyopt = NO_OPTIMISATION

      ! At start the minimum energy is defined as zero

      minim%eng_min = Huge(0.0_wp)

      ! Passage accumulators are initialised in minimise_module
      ! stats%passmin(1) - cycles counter
      ! stats%passmin(2) - access counter
      ! stats%passmin(3) - average cycles
      ! stats%passmin(4) - minimum cycles
      ! stats%passmin(5) - maximum cycles

      ! minim%total number of active particles (excluding frozen sites and massless shells)

      minim%total = 0.0_wp
      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. (config%weight(i) > 1.0e-6_wp .or. config%lfree(i) == 1)) &
          minim%total = minim%total + 1.0_wp
      End Do
      Call gsum(comm, minim%total)
    End If

    ! Step length for relaxation

    If (minim%step_length > zero_plus) Then

      ! Optionally specified

      minim%step = minim%step_length

    Else

      ! default if unspecified

      minim%step = tstep**2

      ! enlarged depending on functionality if defaulted

      If (cons%megcon == 0 .and. pmf%megpmf == 0) Then
        If (rigid%total == 0) Then
          minim%step = 10.0_wp * minim%step
        Else
          minim%step = 5.0_wp * minim%step
        End If
      End If

    End If
    If (minim%keyopt == NO_OPTIMISATION) Then

      ! Initial configuration energy

      minim%eng_0 = stpcfg

      ! Allocate working arrays
      Call minim%init(config%mxatms)

      ! No minimisation is yet attempted

      minim%relaxed = .false.

      ! No RB move is yet attempted

      minim%l_mov = .false.

      ! Minimum needed for a pass for this minimisation cycle

      minim%min_pass = Huge(1.0_wp)

      ! Avoid rdf calculation redundancy

      minim%l_rdf = rdf_collect
      If (rdf_collect) rdf_collect = .false.

      ! Determine optimisation

      If (minim%key == MIN_FORCE) Then
        minim%word = 'force   '
      Else If (minim%key == MIN_ENERGY) Then
        minim%word = 'energy  '
      Else If (minim%key == MIN_DISTANCE) Then
        minim%word = 'distance'
      End If

      ! Print header

      If (l_str) Then
        Write (message, '(3(1x,a),6x,a,10x,a,10x,a,11x,a,5x,a,1p,e11.4,3x,a,e11.4)') &
          'Minimising', minim%word, 'pass', 'eng_tot', 'minim%grad_tol', 'minim%eng_tol', 'minim%dist_tol', 'tol=', &
          minim%tolerance, 'minim%step=', minim%step
        Call info(message, .true.)
        Write (message, "(1x,130('-'))")
        Call info(message, .true.)
      End If
    End If

    ! Load original forces

    Do i = 1, config%natms
      gxx(i) = config%parts(i)%fxx
      gyy(i) = config%parts(i)%fyy
      gzz(i) = config%parts(i)%fzz
    End Do

    ! Minimised energy is current configuration energy

    minim%eng = stpcfg

    ! Calculate pseudo forces and energy for constraint bonds and PMFs

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      If (cons%megcon > 0) Then
        Call constraints_tags(lstitr, cons, config, comm)
        Call constraints_pseudo_bonds(gxx, gyy, gzz, stats, cons, config, comm)
        minim%eng = minim%eng + stats%engcon
      End If

      If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr, pmf, config, comm)
        Call pmf_pseudo_bonds(gxx, gyy, gzz, stats, pmf, config, comm)
        minim%eng = minim%eng + stats%engpmf
      End If
    End If

    ! Average forces over all members of a RB and split torques accordingly

    If (rigid%total > 0) Then
      If (rigid%share) Then
        Call update_shared_units(config, rigid%list_shared, &
                                 rigid%map_shared, gxx, gyy, gzz, domain, comm)
      End If
      Call rigid_bodies_split_torque(gxx, gyy, gzz, txx, tyy, tzz, uxx, uyy, uzz, rigid, config)
    End If

    ! Initialise/get minim%eng_tol & verify relaxed condition

    minim%eng_tol = 0.0_wp
    If (minim%keyopt /= NO_OPTIMISATION) Then
      minim%eng_tol = Abs(1.0_wp - minim%eng2 / minim%eng)
      If (minim%key == MIN_ENERGY) minim%relaxed = (minim%eng_tol < minim%tolerance)
    End If

    ! Current gradient (modulus of the minim%total force)
    ! massless shells and frozen particles have zero forces!

    minim%grad = 0.0_wp
    Do i = 1, config%natms
      minim%grad = minim%grad + gxx(i)**2 + gyy(i)**2 + gzz(i)**2
    End Do
    Call gsum(comm, minim%grad)
    minim%grad = Sqrt(minim%grad)

    ! Get minim%grad_tol & verify relaxed condition

    minim%grad_tol = minim%grad / minim%total
    If (minim%key == MIN_FORCE) minim%relaxed = (minim%grad_tol < minim%tolerance)

    ! Initialise minim%dist_tol

    minim%dist_tol = 0.0_wp

    ! CHECK FOR CONVERGENCE

    If (.not. minim%relaxed) Then

      ! Increment main passage counter

      stats%passmin(1) = stats%passmin(1) + 1.0_wp

      ! minim%min_pass = Min(minim%min_pass,._tol)

      If (minim%key == MIN_FORCE) Then
        minim%min_pass = Min(minim%min_pass, minim%grad_tol)
      Else If (minim%key == MIN_ENERGY) Then
        If (minim%keyopt /= NO_OPTIMISATION) minim%min_pass = Min(minim%min_pass, minim%eng_tol)
      Else If (minim%key == MIN_DISTANCE) Then
        minim%min_pass = Min(minim%min_pass, minim%dist_tol)
      End If

      ! If in mxpass iterations we are not there, give up but
      ! allow for thermo%tension-fold boost in iteration cycle length
      ! for the very first MD minim%step

      If (Nint(stats%passmin(2)) == 0) Then
        If (Nint(stats%passmin(1)) >= 10 * mxpass) Then
          Call warning(330, minim%tolerance, minim%min_pass, 0.0_wp)
          Call error(474)
        End If
      Else
        If (Nint(stats%passmin(1)) >= mxpass) Then
          Call warning(330, minim%tolerance, minim%min_pass, 0.0_wp)
          Call error(474)
        End If
      End If

    Else

      Go To 100

    End If

    If (minim%keyopt == NO_OPTIMISATION) Then

      ! Original configuration energy

      minim%eng0 = minim%eng
      minim%eng1 = minim%eng
      minim%eng2 = minim%eng

      ! Original gradient (modulus of the minim%total force)

      minim%onorm = minim%grad
      minim%grad0 = minim%grad
      minim%grad2 = minim%grad

      ! Set original search direction

      Do i = 1, config%natms
        minim%oxx(i) = gxx(i)
        minim%oyy(i) = gyy(i)
        minim%ozz(i) = gzz(i)
      End Do

      minim%keyopt = LINE_SEARCH
      minim%sgn = 1.0_wp
      minim%stride = minim%sgn * minim%step

    Else If (minim%keyopt == LINE_SEARCH) Then

      ! Line search along chosen direction

      minim%eng1 = minim%eng0
      minim%eng2 = minim%eng

      minim%grad1 = minim%grad2
      minim%grad2 = 0.0_wp
      Do i = 1, config%natms
        minim%grad2 = minim%grad2 + minim%oxx(i) * gxx(i) + minim%oyy(i) * gyy(i) + minim%ozz(i) * gzz(i)
      End Do
      Call gsum(comm, minim%grad2)
      minim%grad2 = minim%sgn * minim%grad2 / minim%onorm

      ! Linear extrapolation to minimum

      If (minim%grad2 < 0.0_wp) Then ! BACK UP FROM THIS DIRECTION
        minim%keyopt = LINEAR_EXTRAPOLATION
        minim%stride = minim%sgn * minim%step * minim%grad2 / (minim%grad1 - minim%grad2)
      Else ! CARRY ON IN THIS DIRECTION
        minim%stride = minim%sgn * minim%step
      End If

    Else If (minim%keyopt == LINEAR_EXTRAPOLATION) Then

      ! Construct conjugate search vector

      minim%eng1 = minim%eng2
      minim%eng2 = minim%eng

      minim%gamma = (minim%grad / minim%grad0)**2
      minim%grad0 = minim%grad
      minim%grad2 = 0.0_wp
      minim%onorm = 0.0_wp
      Do i = 1, config%natms
        minim%oxx(i) = gxx(i) + minim%gamma * minim%oxx(i)
        minim%oyy(i) = gyy(i) + minim%gamma * minim%oyy(i)
        minim%ozz(i) = gzz(i) + minim%gamma * minim%ozz(i)

        minim%onorm = minim%onorm + minim%oxx(i)**2 + minim%oyy(i)**2 + minim%ozz(i)**2
        minim%grad2 = minim%grad2 + minim%oxx(i) * gxx(i) + minim%oyy(i) * gyy(i) + minim%ozz(i) * gzz(i)
      End Do
      Call gsum(comm, minim%onorm)
      minim%onorm = Sqrt(minim%onorm)
      Call gsum(comm, minim%grad2)
      minim%grad2 = minim%grad2 / minim%onorm
      minim%sgn = Sign(1.0_wp, minim%grad2)
      minim%grad2 = minim%sgn * minim%grad2

      minim%keyopt = LINE_SEARCH
      minim%stride = minim%sgn * minim%step

    End If

    ! Move particles to their new positions accordingly

    If (rigid%total > 0) Then

      ! active free particles

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%parts(i)%xxx = config%parts(i)%xxx + minim%stride * minim%oxx(i)
          config%parts(i)%yyy = config%parts(i)%yyy + minim%stride * minim%oyy(i)
          config%parts(i)%zzz = config%parts(i)%zzz + minim%stride * minim%ozz(i)
          minim%dist_tol = Max(minim%dist_tol, minim%oxx(i)**2 + minim%oyy(i)**2 + minim%ozz(i)**2)
        End If
      End Do
      minim%dist_tol = Sqrt(minim%dist_tol) * Abs(minim%stride)

      ! RB particles

      Call rigid_bodies_move(minim%stride, minim%oxx, minim%oyy, minim%ozz, &
                             txx, tyy, tzz, uxx, uyy, uzz, minim%dist_tol, rigid, config)
      minim%l_mov = .true.

    Else

      ! active particles

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%parts(i)%xxx = config%parts(i)%xxx + minim%stride * minim%oxx(i)
          config%parts(i)%yyy = config%parts(i)%yyy + minim%stride * minim%oyy(i)
          config%parts(i)%zzz = config%parts(i)%zzz + minim%stride * minim%ozz(i)
          minim%dist_tol = Max(minim%dist_tol, minim%oxx(i)**2 + minim%oyy(i)**2 + minim%ozz(i)**2)
        End If
      End Do
      minim%dist_tol = Sqrt(minim%dist_tol) * Abs(minim%stride)

    End If
    Call gmax(comm, minim%dist_tol)

    If (minim%key == MIN_DISTANCE) minim%relaxed = (minim%dist_tol < minim%tolerance)

    ! Fit headers in and Close and Open OUTPUT at every 25th print-out

    i = Nint(stats%passmin(1))
    If (l_str) Then
      Write (message, '(1x,i23,1p,4e18.8)') i - 1, minim%eng / engunit, minim%grad_tol, minim%eng_tol, minim%dist_tol
      Call info(message, .true.)
      If (Mod(i, 25) == 0) Then
        Write (message, "(1x,130('-'))")
        Call info(message, .true.)
        Write (message, '(3(1x,a),6x,a,10x,a,10x,a,11x,a,5x,a,1p,e11.4,3x,a,e11.4)') &
          'Minimising', minim%word, 'pass', 'eng_tot', 'minim%grad_tol', 'minim%eng_tol', 'minim%dist_tol', 'tol=', &
          minim%tolerance, 'minim%step=', minim%step
        Call info(message, .true.)
        Write (message, "(1x,130('-'))")
        Call info(message, .true.)

        If (comm%idnode == 0) Then
          Inquire (File=files(FILE_OUTPUT)%filename, Exist=l_out, Position=c_out)
          Call strip_blanks(c_out)
          Call lower_case(c_out)
          If (l_out .and. c_out(1:6) == 'append') Then
            Call files(FILE_OUTPUT)%close ()
            Open (Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Position='append')
          End If
        End If
      End If
    End If

    100 Continue

    minim%transport = (.not. minim%relaxed) ! Transportation flag
    If (minim%relaxed) Then

      ! Final/Only printout

      i = Nint(stats%passmin(1))
      If (.not. l_str) Then
        Write (message, '(3(1x,a),5x,a,10x,a,10x,a,11x,a,5x,a,1p,e11.4,3x,a,e11.4)') &
       'Minimised', minim%word, 'passes', 'eng_tot', 'minim%grad_tol', 'minim%eng_tol', 'minim%dist_tol', 'tol=', minim%tolerance, &
          'minim%step=', minim%step
        Call info(message, .true.)
        Write (message, "(1x,130('-'))")
        Call info(message, .true.)
      End If
      Write (message, '(1x,i23,1p,4e18.8)') i, minim%eng / engunit, minim%grad_tol, minim%eng_tol, minim%dist_tol
      Call info(message, .true.)
      Write (message, "(1x,130('-'))")
      Call info(message, .true.)

      ! Collect passage statistics

      stats%passmin(3) = stats%passmin(2) * stats%passmin(3)
      stats%passmin(2) = stats%passmin(2) + 1.0_wp
      stats%passmin(3) = stats%passmin(3) / stats%passmin(2) + stats%passmin(1) / stats%passmin(2)
      stats%passmin(4) = Min(stats%passmin(1), stats%passmin(4))
      stats%passmin(5) = Max(stats%passmin(1), stats%passmin(5))

      ! Rewind minim%keyopt and main passage counter

      minim%keyopt = NO_OPTIMISATION
      stats%passmin(1) = 0.0_wp

      ! Resume rdf%rdf calculations

      If (minim%l_rdf) rdf_collect = minim%l_rdf

      ! Deallocate working arrays

      Call deallocate_minimise_arrays(minim)

      ! Dump the lowest energy configuration

      If (minim%eng < minim%eng_min) Then
        minim%eng_min = minim%eng

        Call minfile%init('CFGMIN')
        levcfg = 0 ! define level of information in file

        Call write_config(config, minfile, levcfg, i - 1, minim%eng_min / engunit, &
                          io, minim%eng_0 / engunit, netcdf, comm)
      End If

      ! setup new quaternions

      If (minim%l_mov) Then
        If (rigid%share) Then
          Call update_shared_units(config, rigid%list_shared, &
                                   rigid%map_shared, SHARED_UNIT_UPDATE_POSITIONS, domain, comm)
        End If
        Call q_setup(rigid, config, comm)
      End If

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Deallocate (lstitr, Stat=fail(1))
      Call cons%deallocate_work()
      Call pmf%Deallocate_work()
    End If
    If (rigid%total > 0) Then
      Deallocate (txx, tyy, tzz, Stat=fail(6))
      Deallocate (uxx, uyy, uzz, Stat=fail(7))
    End If
    Deallocate (gxx, gyy, gzz, Stat=fail(8))
    If (Any(fail > 0)) Then
      Write (message, '(a,i0)') 'minimise_relax deallocation failure'
      Call error(0, message)
    End If

  End Subroutine minimise_relax

  Subroutine zero_k_optimise(stats, rigid, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine for zero Kelvin temperature optimization
    !
    ! in preparation for integration of equations of motion:
    !
    ! the free particle velocity, V, is set in the direction of the force,
    ! F, when V.F > 0 :: V=F*[(V.F)/(F.F)] and set to zero if V.F < 0 :: V=0.
    ! the same rational is extended to RB dynamics where, additionally to
    ! applying the same strategy to the RB COM velocity change upon the COM
    ! force, alos the angular velocity of the RB, W, is set in the direction
    ! of the torque, T, for when W.T > 0 :: W=T*[(W.T)/(T.T)] and set to zero
    ! if W.T < 0 :: W=0.
    !
    ! care must be taken for:
    !     - to remove any COM motion generation
    !     - to zero angular momentum about centre of mass for non-periodic
    !       system (imcon=0)
    !     - to ensure preservation of thermostat's instantaneous kinetic
    !       energy components and remove spurious temperature fluctuations
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(stats_type),         Intent(InOut) :: stats
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Integer                                  :: fail, i, i1, i2, irgd, j, jrgd, krgd, lrgd, rgdtyp
    Real(Kind=wp)                            :: amx, amy, amz, com(1:3), e_f, e_r, e_t, engkf, &
                                                engkt, fmx, fmy, fmz, fsq, rot(1:9), rotinv(1:9), &
                                                scale, tmp, tmp1, tqx, tqy, tqz, trx, try, trz, &
                                                vdotf, vom(1:3), vpx, vpy, vpz, wxx, wyy, wzz, &
                                                x(1:1), y(1:1), z(1:1)
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer, ggx, ggy, ggz

    ! preserve magnitudes of old instantaneous energies in order to scale back to them

    e_t = 0.5_wp * (stats%strknt(1) + stats%strknt(5) + stats%strknt(9)) ! RBs translational
    e_r = stats%engrot ! RBs rotational
    e_f = stats%engke - e_t ! FPs (free particles)

    ! initialise RB energy components

    engkf = 0.0_wp
    engkt = 0.0_wp

    If (rigid%total > 0) Then
      fail = 0
      Allocate (ggx(1:rigid%max_list * rigid%max_rigid), &
                ggy(1:rigid%max_list * rigid%max_rigid), &
                ggz(1:rigid%max_list * rigid%max_rigid), Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'zero_k_optimise allocation failure'
        Call error(0, message)
      End If

      ! Get the RB particles vectors wrt the RB's COM

      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        ! For all good RBs

        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! COM distances

            ggx(krgd) = config%parts(i)%xxx - rigid%xxx(irgd)
            ggy(krgd) = config%parts(i)%yyy - rigid%yyy(irgd)
            ggz(krgd) = config%parts(i)%zzz - rigid%zzz(irgd)
          End Do
        End If
      End Do

      ! minimum image convention for bond vectors

      Call images(config%imcon, config%cell, krgd, ggx, ggy, ggz)

      ! Free particles

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then

          ! take component of velocity in direction of force

          vdotf = config%vxx(i) * config%parts(i)%fxx + config%vyy(i) * config%parts(i)%fyy + &
                  config%vzz(i) * config%parts(i)%fzz

          If (vdotf < 0.0_wp) Then
            config%vxx(i) = 0.0_wp
            config%vyy(i) = 0.0_wp
            config%vzz(i) = 0.0_wp
          Else
            fsq = config%parts(i)%fxx**2 + config%parts(i)%fyy**2 + config%parts(i)%fzz**2
            scale = vdotf / Max(1.0e-10_wp, fsq)

            config%vxx(i) = config%parts(i)%fxx * scale
            config%vyy(i) = config%parts(i)%fyy * scale
            config%vzz(i) = config%parts(i)%fzz * scale
          End If

        End If
      End Do

      ! RBs

      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then

          ! current rotation matrix

          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

          ! calculate COM force and torque

          fmx = 0.0_wp; fmy = 0.0_wp; fmz = 0.0_wp
          tqx = 0.0_wp; tqy = 0.0_wp; tqz = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force

            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + config%parts(i)%fxx
              fmy = fmy + config%parts(i)%fyy
              fmz = fmz + config%parts(i)%fzz
            End If

            tqx = tqx + ggy(krgd) * config%parts(i)%fzz - ggz(krgd) * config%parts(i)%fyy
            tqy = tqy + ggz(krgd) * config%parts(i)%fxx - ggx(krgd) * config%parts(i)%fzz
            tqz = tqz + ggx(krgd) * config%parts(i)%fyy - ggy(krgd) * config%parts(i)%fxx
          End Do

          ! If the RB has 2+ frozen particles (ill=1) the net torque
          ! must align along the axis of rotation

          If (rigid%frozen(0, rgdtyp) > 1) Then
            i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
            i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

            x(1) = config%parts(i1)%xxx - config%parts(i2)%xxx
            y(1) = config%parts(i1)%yyy - config%parts(i2)%yyy
            z(1) = config%parts(i1)%zzz - config%parts(i2)%zzz

            Call images(config%imcon, config%cell, 1, x, y, z)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) / (x(1)**2 + y(1)**2 + z(1)**2)
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp
          End If

          ! calculate torque in principal frame

          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          If (rigid%frozen(0, rgdtyp) == 0) Then

            ! take component of velocity in direction of force

            vdotf = rigid%vxx(irgd) * fmx + rigid%vyy(irgd) * fmy + rigid%vzz(irgd) * fmz

            If (vdotf < 0.0_wp) Then
              rigid%vxx(irgd) = 0.0_wp
              rigid%vyy(irgd) = 0.0_wp
              rigid%vzz(irgd) = 0.0_wp
            Else
              fsq = fmx**2 + fmy**2 + fmz**2
              scale = vdotf / Max(1.0e-10_wp, fsq)

              rigid%vxx(irgd) = fmx * scale
              rigid%vyy(irgd) = fmy * scale
              rigid%vzz(irgd) = fmz * scale
            End If

          End If

          ! take component of the angular velocity in direction of
          ! the angular acceleration (torque./RI.)

          trx = trx * rigid%rix(2, rgdtyp)
          try = try * rigid%riy(2, rgdtyp)
          trz = trz * rigid%riz(2, rgdtyp)

          vdotf = rigid%oxx(irgd) * trx + rigid%oyy(irgd) * try + rigid%ozz(irgd) * trz

          If (vdotf < 0.0_wp) Then
            rigid%oxx(irgd) = 0.0_wp
            rigid%oyy(irgd) = 0.0_wp
            rigid%ozz(irgd) = 0.0_wp
          Else
            fsq = trx**2 + try**2 + trz**2
            scale = vdotf / Max(1.0e-10_wp, fsq)

            rigid%oxx(irgd) = trx * scale
            rigid%oyy(irgd) = try * scale
            rigid%ozz(irgd) = trz * scale
          End If

          ! update RB members velocities

          Do jrgd = 1, lrgd
            If (rigid%frozen(jrgd, rgdtyp) == 0) Then
              i = rigid%index_local(jrgd, irgd) ! local index of particle/site

              If (i <= config%natms) Then
                x(1) = rigid%x(jrgd, rgdtyp)
                y(1) = rigid%y(jrgd, rgdtyp)
                z(1) = rigid%z(jrgd, rgdtyp)

                ! new atomic velocities in body frame

                vpx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
                vpy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
                vpz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

                ! new atomic velocities in lab frame

                config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
                config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
                config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
              End If
            End If
          End Do
        End If
      End Do

      ! Subtract COM velocity

      Call getvom(vom, rigid, config, comm)

      ! remove centre of mass motion

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        If (rigid%frozen(0, rgdtyp) == 0) Then
          rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
          rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
          rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

          lrgd = rigid%list(-1, irgd)
          Do jrgd = 1, lrgd
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              config%vxx(i) = config%vxx(i) - vom(1)
              config%vyy(i) = config%vyy(i) - vom(2)
              config%vzz(i) = config%vzz(i) - vom(3)
            End If
          End Do
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstresf(stats%strknf, config, comm)
      Call kinstrest(rigid, stats%strknt, comm)

      stats%strkin = stats%strknf + stats%strknt
      stats%engke = 0.5_wp * (stats%strkin(1) + stats%strkin(5) + stats%strkin(9))

      ! update rotational energy

      stats%engrot = getknr(rigid, comm)

      Deallocate (ggx, ggy, ggz, Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'zero_k_optimise deallocation failure'
        Call error(0, message)
      End If
    Else
      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then

          ! take component of velocity in direction of force

          vdotf = config%vxx(i) * config%parts(i)%fxx + config%vyy(i) * config%parts(i)%fyy + &
                  config%vzz(i) * config%parts(i)%fzz

          If (vdotf < 0.0_wp) Then
            config%vxx(i) = 0.0_wp
            config%vyy(i) = 0.0_wp
            config%vzz(i) = 0.0_wp
          Else
            fsq = config%parts(i)%fxx**2 + config%parts(i)%fyy**2 + config%parts(i)%fzz**2
            scale = vdotf / Max(1.0e-10_wp, fsq)

            config%vxx(i) = config%parts(i)%fxx * scale
            config%vyy(i) = config%parts(i)%fyy * scale
            config%vzz(i) = config%parts(i)%fzz * scale
          End If

        End If
      End Do

      ! Subtract COM velocity

      Call getvom(vom, config, comm)

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      ! Update kinetic stress and energy

      Call kinstress(stats%strkin, config, comm)
      stats%engke = 0.5_wp * (stats%strkin(1) + stats%strkin(5) + stats%strkin(9))
    End If

    ! zero angular momentum about centre of mass - non-periodic system

    If (config%imcon == IMCON_NOPBC) Then
      fail = 0
      Allocate (buffer(1:12), Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'zero_k_optimise allocation failure'
        Call error(0, message)
      End If

      ! initialise RB energy components

      engkf = 0.0_wp
      engkt = 0.0_wp

      ! calculate centre of mass position

      Call getcom(config, com, comm)

      If (rigid%total > 0) Then

        ! move to centre of mass origin

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx - com(1)
            config%parts(i)%yyy = config%parts(i)%yyy - com(2)
            config%parts(i)%zzz = config%parts(i)%zzz - com(3)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            rigid%xxx(irgd) = rigid%xxx(irgd) - com(1)
            rigid%yyy(irgd) = rigid%yyy(irgd) - com(2)
            rigid%zzz(irgd) = rigid%zzz(irgd) - com(3)
          End If
        End Do

        ! angular momentum accumulators

        amx = 0.0_wp
        amy = 0.0_wp
        amz = 0.0_wp

        ! rotational inertia accumulators

        rot = 0.0_wp

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            amx = amx + config%weight(i) * (config%parts(i)%yyy * config%vzz(i) - config%parts(i)%zzz * config%vyy(i))
            amy = amy + config%weight(i) * (config%parts(i)%zzz * config%vxx(i) - config%parts(i)%xxx * config%vzz(i))
            amz = amz + config%weight(i) * (config%parts(i)%xxx * config%vyy(i) - config%parts(i)%yyy * config%vxx(i))

            tmp = config%parts(i)%xxx**2 + config%parts(i)%yyy**2 + config%parts(i)%zzz**2
            rot(1) = rot(1) + config%weight(i) * (config%parts(i)%xxx * config%parts(i)%xxx - tmp)
            rot(2) = rot(2) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%yyy
            rot(3) = rot(3) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%zzz
            rot(5) = rot(5) + config%weight(i) * (config%parts(i)%yyy * config%parts(i)%yyy - tmp)
            rot(6) = rot(6) + config%weight(i) * config%parts(i)%yyy * config%parts(i)%zzz
            rot(9) = rot(9) + config%weight(i) * (config%parts(i)%zzz * config%parts(i)%zzz - tmp)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            lrgd = rigid%list(-1, irgd)

            tmp1 = rigid%weight(0, rgdtyp) * Real(rigid%index_local(0, irgd), wp) / Real(lrgd, wp)

            amx = amx + tmp1 * (rigid%yyy(irgd) * rigid%vzz(irgd) - rigid%zzz(irgd) * rigid%vyy(irgd))
            amy = amy + tmp1 * (rigid%zzz(irgd) * rigid%vxx(irgd) - rigid%xxx(irgd) * rigid%vzz(irgd))
            amz = amz + tmp1 * (rigid%xxx(irgd) * rigid%vyy(irgd) - rigid%yyy(irgd) * rigid%vxx(irgd))

            tmp = rigid%xxx(irgd)**2 + rigid%yyy(irgd)**2 + rigid%zzz(irgd)**2

            rot(1) = rot(1) + tmp1 * (rigid%xxx(irgd) * rigid%xxx(irgd) - tmp)
            rot(2) = rot(2) + tmp1 * rigid%xxx(irgd) * rigid%yyy(irgd)
            rot(3) = rot(3) + tmp1 * rigid%xxx(irgd) * rigid%zzz(irgd)
            rot(5) = rot(5) + tmp1 * (rigid%yyy(irgd) * rigid%yyy(irgd) - tmp)
            rot(6) = rot(6) + tmp1 * rigid%yyy(irgd) * rigid%zzz(irgd)
            rot(9) = rot(9) + tmp1 * (rigid%zzz(irgd) * rigid%zzz(irgd) - tmp)
          End If
        End Do

        ! complete rotational inertia matrix

        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)

        ! global sum of rotation

        buffer(1) = amx
        buffer(2) = amy
        buffer(3) = amz
        Do i = 1, 9
          buffer(i + 3) = rot(i)
        End Do

        Call gsum(comm, buffer)

        amx = buffer(1)
        amy = buffer(2)
        amz = buffer(3)
        Do i = 1, 9
          rot(i) = buffer(i + 3)
        End Do

        ! invert rotational inertia matrix

        Call invert(rot, rotinv, tmp)

        ! correction to angular velocity

        wxx = rotinv(1) * amx + rotinv(2) * amy + rotinv(3) * amz
        wyy = rotinv(4) * amx + rotinv(5) * amy + rotinv(6) * amz
        wzz = rotinv(7) * amx + rotinv(8) * amy + rotinv(9) * amz

        ! correction to linear velocity

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) + (wyy * config%parts(i)%zzz - wzz * config%parts(i)%yyy)
            config%vyy(i) = config%vyy(i) + (wzz * config%parts(i)%xxx - wxx * config%parts(i)%zzz)
            config%vzz(i) = config%vzz(i) + (wxx * config%parts(i)%yyy - wyy * config%parts(i)%xxx)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            x(1) = (wyy * rigid%zzz(irgd) - wzz * rigid%yyy(irgd))
            y(1) = (wzz * rigid%xxx(irgd) - wxx * rigid%zzz(irgd))
            z(1) = (wxx * rigid%yyy(irgd) - wyy * rigid%xxx(irgd))

            rigid%vxx(irgd) = rigid%vxx(irgd) + x(1)
            rigid%vyy(irgd) = rigid%vyy(irgd) + y(1)
            rigid%vzz(irgd) = rigid%vzz(irgd) + z(1)

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd) ! local index of particle/site

              If (i <= config%natms) Then
                config%vxx(i) = config%vxx(i) + x(1)
                config%vyy(i) = config%vyy(i) + y(1)
                config%vzz(i) = config%vzz(i) + z(1)
              End If
            End Do
          End If
        End Do

        ! get kinetic energy

        engkf = getknf(config%vxx, config%vyy, config%vzz, config, comm)
        engkt = getknt(rigid, comm)
        stats%engke = engkf + engkt

        ! reset positions to original reference frame

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx + com(1)
            config%parts(i)%yyy = config%parts(i)%yyy + com(2)
            config%parts(i)%zzz = config%parts(i)%zzz + com(3)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            rigid%xxx(irgd) = rigid%xxx(irgd) + com(1)
            rigid%yyy(irgd) = rigid%yyy(irgd) + com(2)
            rigid%zzz(irgd) = rigid%zzz(irgd) + com(3)
          End If
        End Do

      Else

        ! move to centre of mass origin

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx - com(1)
            config%parts(i)%yyy = config%parts(i)%yyy - com(2)
            config%parts(i)%zzz = config%parts(i)%zzz - com(3)
          End If
        End Do

        ! angular momentum accumulators

        amx = 0.0_wp
        amy = 0.0_wp
        amz = 0.0_wp

        ! rotational inertia accumulators

        rot = 0.0_wp

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            amx = amx + config%weight(i) * (config%parts(i)%yyy * config%vzz(i) - &
                                            config%parts(i)%zzz * config%vyy(i))
            amy = amy + config%weight(i) * (config%parts(i)%zzz * config%vxx(i) - &
                                            config%parts(i)%xxx * config%vzz(i))
            amz = amz + config%weight(i) * (config%parts(i)%xxx * config%vyy(i) - &
                                            config%parts(i)%yyy * config%vxx(i))

            tmp = config%parts(i)%xxx**2 + config%parts(i)%yyy**2 + config%parts(i)%zzz**2
            rot(1) = rot(1) + config%weight(i) * (config%parts(i)%xxx * config%parts(i)%xxx - tmp)
            rot(2) = rot(2) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%yyy
            rot(3) = rot(3) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%zzz
            rot(5) = rot(5) + config%weight(i) * (config%parts(i)%yyy * config%parts(i)%yyy - tmp)
            rot(6) = rot(6) + config%weight(i) * config%parts(i)%yyy * config%parts(i)%zzz
            rot(9) = rot(9) + config%weight(i) * (config%parts(i)%zzz * config%parts(i)%zzz - tmp)
          End If
        End Do

        ! complete rotational inertia matrix

        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)

        ! global sum of rotation

        buffer(1) = amx
        buffer(2) = amy
        buffer(3) = amz
        Do i = 1, 9
          buffer(i + 3) = rot(i)
        End Do

        Call gsum(comm, buffer)

        amx = buffer(1)
        amy = buffer(2)
        amz = buffer(3)
        Do i = 1, 9
          rot(i) = buffer(i + 3)
        End Do

        ! invert rotational inertia matrix

        Call invert(rot, rotinv, tmp)

        ! correction to angular velocity

        wxx = rotinv(1) * amx + rotinv(2) * amy + rotinv(3) * amz
        wyy = rotinv(4) * amx + rotinv(5) * amy + rotinv(6) * amz
        wzz = rotinv(7) * amx + rotinv(8) * amy + rotinv(9) * amz

        ! correction to linear velocity

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) + (wyy * config%parts(i)%zzz - wzz * config%parts(i)%yyy)
            config%vyy(i) = config%vyy(i) + (wzz * config%parts(i)%xxx - wxx * config%parts(i)%zzz)
            config%vzz(i) = config%vzz(i) + (wxx * config%parts(i)%yyy - wyy * config%parts(i)%xxx)
          End If
        End Do

        ! get kinetic energy

        stats%engke = getkin(config, config%vxx, config%vyy, config%vzz, comm)

        ! reset positions to original reference frame

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx + com(1)
            config%parts(i)%yyy = config%parts(i)%yyy + com(2)
            config%parts(i)%zzz = config%parts(i)%zzz + com(3)
          End If
        End Do

      End If

      Deallocate (buffer, Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'zero_k_optimise deallocation failure'
        Call error(0, message)
      End If
    End If

    ! to remove spurious temperature fluctuations ensure preservation
    ! of thermostat's instantaneous translational and rotational energies
    ! equipartitioning of DoFs may be lost as transfers of energy between
    ! translational and rotational DoFs may happen

    ! apply temperature components scaling

    engkf = stats%engke - engkt
    If (engkf + engkt + stats%engrot > 1.0e-6_wp .and. e_f + e_t + e_r > 1.0e-6_wp) Then
      If (rigid%total > 0) Then
        tmp = Sqrt((e_f + e_t + e_r) / (engkf + engkt + stats%engrot))
        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) * tmp
            config%vyy(i) = config%vyy(i) * tmp
            config%vzz(i) = config%vzz(i) * tmp
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then

            ! new angular velocity

            rigid%oxx(irgd) = rigid%oxx(irgd) * tmp
            rigid%oyy(irgd) = rigid%oyy(irgd) * tmp
            rigid%ozz(irgd) = rigid%ozz(irgd) * tmp

            ! new translational velocity

            If (rigid%frozen(0, rgdtyp) == 0) Then
              rigid%vxx(irgd) = rigid%vxx(irgd) * tmp
              rigid%vyy(irgd) = rigid%vyy(irgd) * tmp
              rigid%vzz(irgd) = rigid%vzz(irgd) * tmp
            End If

            ! new rotational matrix

            Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

            Do jrgd = 1, lrgd
              If (rigid%frozen(jrgd, rgdtyp) == 0) Then ! Apply restrictions
                i = rigid%index_local(jrgd, irgd) ! local index of particle/site

                If (i <= config%natms) Then
                  x(1) = rigid%x(jrgd, rgdtyp)
                  y(1) = rigid%y(jrgd, rgdtyp)
                  z(1) = rigid%z(jrgd, rgdtyp)

                  ! site velocity in body frame

                  wxx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
                  wyy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
                  wzz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

                  ! new atomic velocities in lab frame

                  config%vxx(i) = rot(1) * wxx + rot(2) * wyy + rot(3) * wzz + rigid%vxx(irgd)
                  config%vyy(i) = rot(4) * wxx + rot(5) * wyy + rot(6) * wzz + rigid%vyy(irgd)
                  config%vzz(i) = rot(7) * wxx + rot(8) * wyy + rot(9) * wzz + rigid%vzz(irgd)
                End If
              End If
            End Do
          End If
        End Do

        ! update kinetic energy and stress

        Call kinstresf(stats%strknf, config, comm)
        Call kinstrest(rigid, stats%strknt, comm)

        stats%strkin = stats%strknf + stats%strknt
        stats%engke = 0.5_wp * (stats%strkin(1) + stats%strkin(5) + stats%strkin(9))

        ! update rotational energy

        stats%engrot = getknr(rigid, comm)
      Else
        tmp = Sqrt(e_f / stats%engke)
        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) * tmp
            config%vyy(i) = config%vyy(i) * tmp
            config%vzz(i) = config%vzz(i) * tmp
          End If
        End Do

        ! Update kinetic stress and energy

        Call kinstress(stats%strkin, config, comm)
        stats%engke = 0.5_wp * (stats%strkin(1) + stats%strkin(5) + stats%strkin(9))
      End If
    Else ! zero them and let's see if we can crash
      Do i = 1, config%natms
        config%vxx(i) = 0.0_wp; config%vyy(i) = 0.0_wp; config%vzz(i) = 0.0_wp
      End Do

      If (rigid%total > 0) Then
        Do irgd = 1, rigid%n_types
          rigid%vxx(irgd) = 0.0_wp; rigid%vyy(irgd) = 0.0_wp; rigid%vzz(irgd) = 0.0_wp
          rigid%oxx(irgd) = 0.0_wp; rigid%oyy(irgd) = 0.0_wp; rigid%ozz(irgd) = 0.0_wp
        End Do
      End If
    End If

  End Subroutine zero_k_optimise
End Module minimise
