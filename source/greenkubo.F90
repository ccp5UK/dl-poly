Module greenkubo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring arrays for green-kubo relations
  ! calculations to calculate transport properties
  !
  ! copyright - daresbury laboratory
  ! author    - m.a.seaton june 2014
  ! contrib   - i.t.todorov july 2014
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: comms_type,&
                             gcheck,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: nvafdt,&
                             zero_plus
  Use errors_warnings, Only: error,&
                             info
  Use flow_control,    Only: RESTART_KEY_OLD
  Use kinds,           Only: wi,&
                             wp
  Use site,            Only: site_type

  Implicit None

  Private

  Type, Public :: greenkubo_type
    Private

    !> VAF printing switch
    Logical, Public          :: l_print
    !> VAF time averaging switch
    Logical, Public          :: l_average
    !> VAF sampling frequency in steps
    Integer(Kind=wi), Public :: freq = 1
    !> VAF sample size
    Integer(Kind=wi), Public :: binsize = 0
    !> VAF timestep start
    Integer(Kind=wi), Public :: t_start = -1
    !> VAF simultaneously overlapping samples
    Integer(Kind=wi), Public :: samp = 0

    Real(Kind=wp), Public :: vafcount = 0.0_wp
    Integer(Kind=wi), Allocatable, Public :: step(:)
    Real(Kind=wp), Allocatable, Public :: vxi(:, :), vyi(:, :), vzi(:, :)
    Real(Kind=wp), Allocatable, Public :: vafdata(:, :), time(:), vaf(:, :)
    Logical :: newjob
    !    Real( Kind = wp ), Allocatable :: stxx(:),stxy(:),stxz(:),styy(:),styz(:),stzz(:)
    !    Real( Kind = wp ), Allocatable :: gkpot(:),tcond(:),tctime(:)
  Contains

    Private
    Procedure, Public :: init => allocate_greenkubo_arrays
    Final             :: cleanup
  End Type greenkubo_type

  Public :: vaf_compute, vaf_collect, vaf_write

Contains

  Subroutine allocate_greenkubo_arrays(green, mxatms, mxatyp)
    Class(greenkubo_Type), Intent(InOut) :: green
    Integer,               Intent(In   ) :: mxatms, mxatyp

    Integer                 :: i
    Integer, Dimension(1:7) :: fail

    fail = 0

    Allocate (green%step(1:green%samp), Stat=fail(1))
    Allocate (green%vxi(1:mxatms, 1:green%samp), &
              green%vyi(1:mxatms, 1:green%samp), &
              green%vzi(1:mxatms, 1:green%samp), Stat=fail(2))
    Allocate (green%vafdata(0:green%binsize, 1:green%samp * (mxatyp + 1)), &
              green%time(0:green%binsize), &
              green%vaf(0:green%binsize, 1:mxatyp), Stat=fail(3))
    !    Allocate (green%gkpot(1:mxatms),                                                              Stat = fail(4))
    !    Allocate (green%stxx(1:mxatms),green%stxy(1:mxatms),green%stxz(1:mxatms),                                 Stat = fail(5))
    !    Allocate (green%styy(1:mxatms),green%styz(1:mxatms),green%stzz(1:mxatms),                                 Stat = fail(6))
    !    Allocate (green%tcond(0:mxgkstk),green%tctime(0:mxgkstk),                                           Stat = fail(7))

    If (Any(fail > 0)) Call error(1080)

    green%vxi = 0.0_wp; green%vyi = 0.0_wp; green%vzi = 0.0_wp
    green%vaf = 0.0_wp; green%time = 0.0_wp
    green%vafdata = 0.0_wp

    ! setup green%step counters for samples, starting after equilibration

    Do i = 1, green%samp
      green%step(i) = (1 - i) * green%freq - 1
    End Do

    !    sxx = 0.0_wp ; sxy = 0.0_wp ; sxz = 0.0_wp ; syy = 0.0_wp; syz = 0.0_wp; szz = 0.0_wp
    !    green%tcond = 0.0_wp ; green%tctime = 0.0_wp

  End Subroutine allocate_greenkubo_arrays

  Subroutine vaf_collect(config, mxatyp, leql, nsteql, nstep, time, green, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistics for velocity
    ! autocorrelation functions
    !
    ! copyright - daresbury laboratory
    ! author    - m.a.seaton june 2014
    ! amended   - i.t.todorov july 2014
    ! amended   - m.a.seaton march 2020
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Integer(Kind=wi),         Intent(In   ) :: mxatyp
    Logical,                  Intent(In   ) :: leql
    Integer,                  Intent(In   ) :: nsteql, nstep
    Real(Kind=wp),            Intent(In   ) :: time
    Type(greenkubo_Type),     Intent(InOut) :: green
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail, i, j, k, l, nsum
    Real(Kind=wp), Allocatable :: vafcoll(:)

    fail = 0
    Allocate (vafcoll(1:mxatyp), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vaf_collect allocation failure'
      Call error(0, message)
    End If

    ! set VAF timestep start

    If (green%t_start < 0) green%t_start = Merge(nsteql, nstep, leql)

    ! advance time green%step, calculate green%vaf contribution for each sample

    Do j = 1, green%samp
      green%step(j) = green%step(j) + 1
      k = green%step(j)

      If (k > 0 .and. k <= green%binsize) Then
        vafcoll = 0.0_wp
        Do i = 1, config%natms
          If (config%lfrzn(i) == 0) Then
            l = config%ltype(i)
       vafcoll(l) = vafcoll(l) + config%vxx(i) * green%vxi(i, j) + config%vyy(i) * green%vyi(i, j) + config%vzz(i) * green%vzi(i, j)
          End If
        End Do

        Do l = 1, mxatyp
          green%vafdata(k, (j - 1) * (mxatyp + 1) + l) = vafcoll(l)
        End Do

        green%vafdata(k, j * (mxatyp + 1)) = time
      End If

      ! get velocity autocorrelation function at end of time span for sample
      ! and reset green%step for next sample

      If (k == green%binsize) Then
        green%vafcount = green%vafcount + 1.0_wp

        If (comm%mxnode > 1) Then
          nsum = config%mxbuff / (green%binsize + 1)

          l = (j - 1) * (mxatyp + 1) ! avoid summing up timing information
          Do i = 1, mxatyp, nsum
            Call gsum(comm, green%vafdata(:, l + i:l + Min(i + nsum - 1, mxatyp)))
          End Do
        End If

        If (green%l_average) Then ! if time-averaging, add green%vaf data to sampling array
          Do i = 0, green%binsize
            green%time(i) = green%vafdata(i, j * (mxatyp + 1))
            Do l = 1, mxatyp
              green%vaf(i, l) = green%vaf(i, l) + green%vafdata(i, (j - 1) * (mxatyp + 1) + l)
            End Do
          End Do
        Else ! if not time-averaging, move green%vaf data to sampling array
          green%time(0:green%binsize) = green%vafdata(0:green%binsize, j * (mxatyp + 1))
          Do l = 1, mxatyp
            green%vaf(0:green%binsize, l) = green%vafdata(0:green%binsize, (j - 1) * (mxatyp + 1) + l)
          End Do
        End If
        green%step(j) = Min (0, green%binsize-green%freq)
      End If

      ! reset reference velocities and calculate first VAF value
      ! at start of sample

      If (green%step(j)==0) Then
        green%vxi(1:config%natms, j) = config%vxx(1:config%natms)
        green%vyi(1:config%natms, j) = config%vyy(1:config%natms)
        green%vzi(1:config%natms, j) = config%vzz(1:config%natms)

        vafcoll = 0.0_wp
        Do i = 1, config%natms
          If (config%lfrzn(i) == 0) Then
            l = config%ltype(i)
 vafcoll(l) = vafcoll(l) + green%vxi(i, j) * green%vxi(i, j) + green%vyi(i, j) * green%vyi(i, j) + green%vzi(i, j) * green%vzi(i, j)
          End If
        End Do

        green%vafdata(0, j * (mxatyp + 1)) = time
        Do l = 1, mxatyp
          green%vafdata(0, (j - 1) * (mxatyp + 1) + l) = vafcoll(l)
        End Do
      End If
    End Do

    Deallocate (vafcoll, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vaf_collect deallocation failure'
      Call error(0, message)
    End If

  End Subroutine vaf_collect

  Subroutine vaf_compute(tstep, num_type_nf, mxatyp, green)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating velocity autocorrelation
    ! functions from accumulated data
    !
    ! copyright - daresbury laboratory
    ! author    - m.a.seaton june 2014
    ! amended   - i.t.todorov july 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),               Intent(In   ) :: tstep
    Real(Kind=wp), Dimension(:), Intent(In   ) :: num_type_nf
    Integer(Kind=wi),            Intent(In   ) :: mxatyp
    Type(greenkubo_type),        Intent(In   ) :: green

    Character(Len=256) :: message
    Integer            :: i
    Real(Kind=wp)      :: factor, gvaf, numt, ovaf, time0, timei

    Call info('velocity autocorrelation functions', .true.)
    If (green%l_average) Then
      Write (message, '(a,i8,a)') 'calculated using ', Nint(green%vafcount), ' samples'
    Else
      Write (message, '(a,i8,a,f10.4,a)') 'calculated using sample ', &
        Nint(green%vafcount), ' starting at ', green%time(0), ' ps'
    End If
    Call info(message, .true.)

    time0 = green%time(0)
    numt = Sum(num_type_nf(1:mxatyp))
    factor = 1.0_wp / Sum(green%vaf(0, 1:mxatyp))
    ovaf = Sum(green%vaf(0, 1:mxatyp)) / Real(numt, Kind=wp)
    If (green%l_average) ovaf = ovaf / green%vafcount

    Write (message, '(a,1p,e16.8)') 'absolute value at origin (3kT/m) = ', ovaf
    Call info(message, .true.)

    ! loop over time steps

    Do i = 0, green%binsize

      ! determine time

      If (green%l_average) Then
        timei = tstep * Real(i, Kind=wp)
      Else
        timei = green%time(i) - time0
      End If

      ! null it if number of particles is zero

      If (numt > zero_plus) Then
        gvaf = Sum(green%vaf(i, 1:mxatyp)) * factor
      Else ! THIS CAN NEVER HAPPEN AS TOTAL NON-FROZEN PARTICLES IS ALWAYS > 1
        gvaf = 0.0_wp
      End If

      ! print out information

      Write (message, '(f10.4,1p,e14.6)') timei, gvaf
      Call info(message, .true.)

    End Do

  End Subroutine vaf_compute

  Subroutine vaf_write(config, keyres, nstep, tstep, green, sites, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for writing VAFDAT file at selected intervals
    ! in simulation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & m.a.seaton february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Integer,                  Intent(In   ) :: keyres, nstep
    Real(Kind=wp),            Intent(In   ) :: tstep
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(site_type),          Intent(In   ) :: sites
    Type(comms_type),         Intent(InOut) :: comm

    Integer       :: i, j
    Logical       :: lexist
    Real(Kind=wp) :: factor, gvaf, numt, ovaf, time0, timei

    If (green%t_start < 0) Return
    If (.not. (nstep >= green%t_start + green%binsize .and. Mod(nstep - green%t_start - green%binsize, green%freq) == 0)) Return

    If (green%newjob) Then
      green%newjob = .false.

      ! If the keyres = RESTART_KEY_OLD, is VAFDAT old (does it exist) and
      ! how many frames and records are in there

      Do i = 1, sites%mxatyp
        lexist = .true.
        If (keyres == RESTART_KEY_OLD) Then
          If (comm%idnode == 0) Inquire (File='VAFDAT_'//sites%unique_atom(i), Exist=lexist)
          Call gcheck(comm, lexist)
        Else
          lexist = .false.
        End If

        ! Generate file if non-existent or replace it if outputting time-averaged VAF

        If ((.not. lexist) .or. green%l_average) Then
          If (comm%idnode == 0) Then
            Open (Unit=nvafdt, File='VAFDAT_'//sites%unique_atom(i), Status='replace')
            Write (nvafdt, '(a)') config%cfgname
            Close (Unit=nvafdt)
          End If
        End If
      End Do
    End If

    time0 = green%time(0)

    ! loop over species types

    Do j = 1, sites%mxatyp
      numt = sites%num_type_nf(j)
      factor = 1.0_wp
      If (Abs(green%vaf(0, j)) > 0.0e-6_wp) factor = 1.0_wp / green%vaf(0, j)
      ovaf = green%vaf(0, j) / Real(numt, Kind=wp)
      If (green%l_average) ovaf = ovaf / green%vafcount

      ! replace file (if outputting time-averaged VAF); prepare for data set

      If (comm%idnode == 0) Then
        If (green%l_average) Then
          Open (Unit=nvafdt, File='VAFDAT_'//sites%unique_atom(j), Status='replace')
          Write (nvafdt, '(a)') config%cfgname
          Close (Unit=nvafdt)
        End If
        Open (Unit=nvafdt, File='VAFDAT_'//sites%unique_atom(j), Position='append')
        Write (nvafdt, '(a8,i10,1p,e16.8,f20.6)') sites%unique_atom(j), green%binsize, ovaf, time0
      End If

      ! Then loop over time steps

      Do i = 0, green%binsize

        ! determine time

        If (green%l_average) Then
          timei = tstep * Real(i, Kind=wp)
        Else
          timei = green%time(i) - time0
        End If

        ! null it if number of particles is zero

        If (numt > zero_plus) Then
          gvaf = green%vaf(i, j) * factor
        Else
          gvaf = 0.0_wp
        End If

        ! print out information

        If (comm%idnode == 0) Write (nvafdt, "(1p,2e14.6)") timei, gvaf

      End Do

      If (comm%idnode == 0) Then
        If (.not. green%l_average) Write (nvafdt, '(2/)')
        Close (Unit=nvafdt)
      End If
    End Do
  End Subroutine vaf_write

  !> Deallocate greenkubo arrays
  Subroutine cleanup(green)
    Type(greenkubo_type) :: green

    If (Allocated(green%step)) Then
      Deallocate (green%step)
    End If

    If (Allocated(green%vxi)) Then
      Deallocate (green%vxi)
    End If
    If (Allocated(green%vyi)) Then
      Deallocate (green%vyi)
    End If
    If (Allocated(green%vzi)) Then
      Deallocate (green%vzi)
    End If

    If (Allocated(green%vafdata)) Then
      Deallocate (green%vafdata)
    End If
    If (Allocated(green%time)) Then
      Deallocate (green%time)
    End If
    If (Allocated(green%vaf)) Then
      Deallocate (green%vaf)
    End If
  End Subroutine cleanup
End Module greenkubo
