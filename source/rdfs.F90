Module rdfs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring RDF property variables and arrays
  ! including USR (umbrella sampling restraint) RDF
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov november 2016
  ! contrib   - a.b.g.chalk january 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use site, Only : site_type
  Use configuration, Only : configuration_type
  Use comms,  Only : comms_type,gsum
  Use constants,  Only : fourpi,boltz,delr_max,nrdfdt,npdfdt,npdgdt, &
    engunit,zero_plus
  Use errors_warnings, Only : error,info
  Use neighbours, Only : neighbours_type

  Implicit None

  Private


  !> Type containing RDF and USR RDF data
  Type, Public :: rdf_type
    Private

    !> RDF collection flag
    Logical, Public :: l_collect
    !> RDF recording flag

    Logical, Public :: l_print
    !> Total number of blocks?
    Integer( Kind = wi ), Public :: num_blocks = 25

    !> RDF collection frequency (in steps)
    Integer( Kind = wi ), Public :: freq

    !> Number of configurations used in RDF calculation
    Integer( Kind = wi ), Public :: n_configs = 0
    !> Number of rdf look up pairs
    Integer( Kind = wi ), Public :: n_pairs = 0

    !> Number of configurations used in USR RDF calculation
    Integer( Kind = wi ), Public :: n_configs_usr = 0
    !> USR RDF cutoff
    Real( Kind = wp ), Public :: cutoff_usr   = 0.0_wp

    Integer( Kind = wi ), Allocatable, Public :: list(:)

    !> RDF array
    Real( Kind = wp ), Allocatable, Public :: rdf(:,:)
    !> USR RDF array
    Real( Kind = wp ), Allocatable, Public :: usr(:)

    !> Block averages
    Real( Kind = wp ), Allocatable :: block_averages(:,:,:,:)
    !> Size of a block
    Integer( Kind = wi ) :: block_size
    !> Current block number
    Integer( Kind = wi ) :: block_number = 1

    !> Temporary rdf array
    Real( Kind = wp ), Allocatable :: tmp_rdf(:,:,:)
    Logical :: tmp_rdf_sync = .false.

    Logical, Public :: l_errors_block = .false.
    Logical, Public :: l_errors_jack = .false.

    !> Maximum number of RDF pairs
    Integer( Kind = wi ), Public :: max_rdf
    !> Maximum number of RDF grid points
    Integer( Kind = wi ), Public :: max_grid
    !> Maximum number of USR RDF grid points
    Integer( Kind = wi ), Public :: max_grid_usr
    Real( Kind = wp ), Public :: rbin

  Contains
    Private

    Procedure, Public :: init => allocate_rdf_arrays
    Procedure, Public :: init_block => allocate_block_average_array
    Final :: cleanup
  End Type rdf_type

  Public :: rdf_compute, calculate_errors, calculate_errors_jackknife, &
    usr_compute, usr_collect, rdf_collect, rdf_excl_collect, &
    rdf_frzn_collect, rdf_increase_block_number

Contains

  Subroutine allocate_rdf_arrays(T)
  Class( rdf_type ) :: T
    Integer, Dimension(1:3) :: fail

    fail = 0

    Allocate (T%list(1:T%max_rdf), stat=fail(1))
    Allocate (T%rdf(1:T%max_grid,1:T%max_rdf), stat=fail(2))
    Allocate (T%usr(1:T%max_grid_usr), stat=fail(3))

    If (Any(fail > 0)) Call error(1016)

    T%list = 0

    T%rdf = 0.0_wp
    T%usr = 0.0_wp
  End Subroutine allocate_rdf_arrays

  Subroutine allocate_block_average_array(T,nstrun,ntype_atom)
  Class( rdf_type) :: T
    Integer( Kind = wi ), Intent( In    ) :: nstrun
    Integer( Kind = wi ), Intent( In    ) :: ntype_atom

    Integer :: temp1, temp2
    Integer, Dimension(1:2) :: fail

    T%block_size = nstrun/(T%num_blocks-1)
    if(T%block_size < 2) then
      T%block_size = 2
    endif

    temp1 = T%max_rdf + 16-Mod(T%max_rdf,16)
    temp2 = T%max_grid + 16-Mod(T%max_grid,16)

    Allocate(T%block_averages(1:ntype_atom,1:ntype_atom,1:T%max_grid,1:T%num_blocks+1), Stat = fail(1))
    Allocate(T%tmp_rdf( 1:temp2,1:temp1, 1:T%num_blocks+1 ), Stat = fail(2))

    If (Any(fail > 0)) Call error(1016)
    T%block_averages = 0.0_wp
    T%tmp_rdf = 0.0_wp
  End Subroutine allocate_block_average_array

  Subroutine rdf_collect(iatm,rrt,neigh,config,rdf)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistic for radial
    ! distribution functions
    !
    ! Note: to be used as part of two_body_forces
    !
    ! copyright - daresbury laboratory
    ! author    - t.forester march 1994
    ! amended   - i.t.todorov november 2014
    ! contrib   - a.b.g.chalk january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type( neighbours_type), Intent( In    ) :: neigh
    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: rrt
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( configuration_type ), Intent( InOut ) :: config

    Integer                 :: idi,jatm,ai,aj,keyrdf,kk,ll,m
    Real( Kind = wp )       :: rdelr,rrr

    ! set cutoff condition for pair forces and grid interval for rdf%rdf tables
    rdelr= Real(rdf%max_grid,wp)/neigh%cutoff

    ! global identity and type of iatm
    idi=config%ltg(iatm)
    ai=config%ltype(iatm)

    ! start of primary loop for rdf%rdf accumulation
    Do m=1,neigh%list(0,iatm)
      ! atomic and type indices
      jatm=neigh%list(m,iatm)
      aj=config%ltype(jatm)

      If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
        ! rdf%rdf function indices
        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=rdf%list(keyrdf)
        ! only for valid interactions specified for a look up
        If (kk > 0 .and. kk <= rdf%n_pairs) Then
          ! apply truncation of potential
          rrr=rrt(m)
          If (rrr < neigh%cutoff) Then
            ll=Min(1+Int(rrr*rdelr),rdf%max_grid)
            ! accumulate correlation
            rdf%rdf(ll,kk) = rdf%rdf(ll,kk) + 1.0_wp
            If(rdf%l_errors_block .or. rdf%l_errors_jack) Then
              rdf%tmp_rdf(ll,kk,rdf%block_number) = rdf%tmp_rdf(ll,kk,rdf%block_number) + 1.0_wp
            End If
          End If
        End If
      End If

    End Do

  End Subroutine rdf_collect

  Subroutine rdf_compute(lpana,rcut,temp,sites,rdf,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating radial distribution functions
    ! from accumulated data
    !
    ! copyright - daresbury laboratory
    ! author    - t.forester & i.t.todorov march 2016
    ! contrib   - a.v.brukhno january 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Logical          , Intent( In    ) :: lpana
    Real( Kind = wp ), Intent( In    ) :: rcut,temp
    Type( site_type ), Intent( In    ) :: sites
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type(comms_type), Intent( InOut )  :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Logical           :: zero
    Integer           :: fail,ngrid,i,ia,ib,kk,ig,ll
    Real( Kind = wp ) :: kT2engo,delr,rdlr,dgrid,pdfzero,      &
      factor1,rrr,dvol,gofr,gofr1,sum,sum1, &
      fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

    Real( Kind = wp ), Allocatable :: dstdrdf(:,:)
    Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

    Character ( Len = 256 )  :: message,messages(2)

    If (lpana) Then
      fail = 0
      Allocate (dstdrdf(0:rdf%max_grid,1:rdf%n_pairs),pmf(0:rdf%max_grid+2),vir(0:rdf%max_grid+2), Stat = fail)
      If (fail > 0) Then
        Write(message,'(a)') 'rdf_compute - allocation failure'
        Call error(0,message)
      End If
    End If

    ! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

    kT2engo = boltz*temp/engunit

    ! grid interval for rdf%rdf tables

    delr = rcut/Real(rdf%max_grid,wp)
    rdlr = 1.0_wp/delr

    ! resampling grid and grid interval for rdf%rdf tables

    ngrid = Max(Nint(rcut/delr_max),rdf%max_grid)
    dgrid = rcut/Real(ngrid,wp)

    Write(messages(1),'(a)') 'radial distribution functions:'
    Write(messages(2),'(2x,a,i8,a)') 'calculated using ', rdf%n_configs, ' configurations'
    Call info(messages,2,.true.)

    ! open RDF file and Write headers

    If (comm%idnode == 0) Then
      Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
      Write(nrdfdt,'(a)') config%cfgname
      Write(nrdfdt,'(2i10)') rdf%n_pairs,rdf%max_grid
    End If

    ! the lower bound to nullify the nearly-zero histogram (PDF) values

    pdfzero = 1.0e-6_wp

    ! for all possible unique type-to-type pairs

    Do ia=1,sites%ntype_atom
      Do ib=ia,sites%ntype_atom

        ! number of the interaction by its rdf%rdf key
        If (ia == ib .and. sites%num_type(ia) < 2) Cycle

        kk=rdf%list(ib*(ib-1)/2+ia)

        ! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= rdf%n_pairs) Then
          Write(messages(1),'(2x,a,2(1x,a8))') 'g(r): ',sites%unique_atom(ia),sites%unique_atom(ib)
          Write(messages(2),'(8x,a1,6x,a4,9x,a4)') 'r','g(r)','n(r)'
          Call info(messages,2,.true.)
          If (comm%idnode == 0) Then
            Write(nrdfdt,'(2a8)') sites%unique_atom(ia),sites%unique_atom(ib)
          End If

          ! global sum of data on all nodes

          Call gsum(comm,rdf%rdf(1:rdf%max_grid,kk))

          ! normalisation factor

          factor1=config%volm*sites%dens(ia)*sites%dens(ib)*Real(rdf%n_configs,wp)
          If (ia == ib) factor1=factor1*0.5_wp*(1.0_wp-1.0_wp/sites%num_type(ia))

          ! running integration of rdf%rdf

          sum=0.0_wp

          ! loop over distances

          zero=.true.
          Do i=1,rdf%max_grid
            If (zero .and. i < (rdf%max_grid-3)) zero=(rdf%rdf(i+2,kk) <= 0.0_wp)

            gofr= rdf%rdf(i,kk)/factor1
            sum = sum + gofr*sites%dens(ib)

            rrr = (Real(i,wp)-0.5_wp)*delr
            dvol= fourpi*delr*(rrr**2+delr**2/12.0_wp)
            gofr= gofr/dvol

            ! zero it if < pdfzero

            If (gofr < pdfzero) Then
              gofr1 = 0.0_wp
            Else
              gofr1 = gofr
            End If

            If (sum < pdfzero) Then
              sum1 = 0.0_wp
            Else
              sum1 = sum
            End If

            ! print out information

            If (.not.zero) Then
              Write(message,'(f10.4,1p,2e14.6)') rrr,gofr1,sum1
              Call info(message,.true.)
            End If
            If (comm%idnode == 0) Then
              Write(nrdfdt,"(1p,2e14.6)") rrr,gofr
            End If

            ! We use the non-normalised tail-truncated RDF version,
            ! rdf%rdf...1 (not pdf...) in order to exclude the nearly-zero
            ! rdf%rdf... noise in PMF, otherwise the PMF = -ln(PDF)
            ! would have poorly-defined noisy "borders/walls"

            If (lpana) dstdrdf(i,kk) = gofr1 ! RDFs density
          End Do
        Else
          If (lpana) dstdrdf(:,kk) = 0.0_wp ! RDFs density
        End If

      End Do
    End Do

    If (comm%idnode == 0) Close(Unit=nrdfdt)

    ! Only when PDF analysis is requested

    If (lpana) Then

      ! open PDF files and write headers

      If (comm%idnode == 0) Then
        Open(Unit=npdgdt, File='VDWPMF', Status='replace')
        Write(npdgdt,'(a)') '# '//config%cfgname
        Write(npdgdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',delr*Real(rdf%max_grid,wp),rdf%max_grid,delr,rdf%n_pairs, &
          '   conversion factor(kT -> energy units) =',kT2engo

        Open(Unit=npdfdt, File='VDWTAB', Status='replace')
        Write(npdfdt,'(a)') '# '//config%cfgname
        Write(npdfdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',dgrid*Real(ngrid,wp),ngrid,dgrid,rdf%n_pairs, &
          '   conversion factor(kT -> energy units) =',kT2engo
      End If

      ! loop over all valid RDFs

      Do ia=1,sites%ntype_atom
        Do ib=ia,sites%ntype_atom

          If (ia == ib .and. sites%num_type(ia) < 2) Cycle
          ! number of the interaction by its rdf%rdf key

          kk=rdf%list(ib*(ib-1)/2+ia)

          ! only for valid interactions specified for a look up

          If (kk > 0 .and. kk <= rdf%n_pairs) Then
            If (comm%idnode == 0) Then
              Write(npdgdt,'(/,a2,2a8)') '# ',sites%unique_atom(ia),sites%unique_atom(ib)
              Write(npdfdt,'(/,a2,2a8)') '# ',sites%unique_atom(ia),sites%unique_atom(ib)
            End If

            ! Smoothen and get derivatives

            ! RDFs -> 1 at long distances, so we do not shift the PMFs
            ! but instead put a cap over those, -Log(pdfzero) - hence,
            ! the upper bound for the PMF due to the zero RDF values

            fed0  = -Log(pdfzero)
            dfed0 = 10.0_wp
            dfed  = 10.0_wp

            Do ig=1,rdf%max_grid
              tmp = Real(ig,wp)-0.5_wp
              rrr = tmp*delr

              If (dstdrdf(ig,kk) > zero_plus) Then
                fed = -Log(dstdrdf(ig,kk)+pdfzero) !-fed0
                If (fed0 <= zero_plus) Then
                  fed0 = fed
                  fed  = fed0
                  !                       fed = 0.0_wp
                End If

                If (ig < rdf%max_grid-1) Then
                  If (dstdrdf(ig+1,kk) <= zero_plus .and. dstdrdf(ig+2,kk) > zero_plus) &
                    dstdrdf(ig+1,kk) = 0.5_wp*(dstdrdf(ig,kk)+dstdrdf(ig+2,kk))
                End If
              Else
                fed = fed0
                !                    fed = 0.0_wp
              End If

              If      (ig == 1) Then
                If      (dstdrdf(ig,kk) > zero_plus .and. dstdrdf(ig+1,kk) > zero_plus) Then
                  dfed = Log(dstdrdf(ig+1,kk)/dstdrdf(ig,kk))
                Else If (dfed > 0.0_wp) Then
                  dfed = dfed0
                Else
                  dfed =-dfed0
                End If
              Else If (ig == rdf%max_grid) Then
                If      (dstdrdf(ig,kk) > zero_plus .and. dstdrdf(ig-1,kk) > zero_plus) Then
                  dfed = Log(dstdrdf(ig,kk)/dstdrdf(ig-1,kk))
                Else If (dfed > 0.0_wp) Then
                  dfed = dfed0
                Else
                  dfed =-dfed0
                End If
              Else If (dstdrdf(ig-1,kk) > zero_plus) Then
                If (dstdrdf(ig+1,kk) > zero_plus) Then
                  dfed = 0.5_wp*(Log(dstdrdf(ig+1,kk)/dstdrdf(ig-1,kk)))
                Else
                  dfed = 0.5_wp*Log(dstdrdf(ig-1,kk))
                End If
              Else If (dstdrdf(ig+1,kk) > zero_plus) Then
                dfed =-0.5_wp*Log(dstdrdf(ig+1,kk))
              Else If (dfed > 0.0_wp) Then
                dfed = dfed0
              Else
                dfed =-dfed0
              End If

              pmf(ig) = fed
              vir(ig) = dfed

              ! Print

              If (comm%idnode == 0) &
                Write(npdgdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*tmp
            End Do

            ! Define edges

            pmf(0)        = 2.0_wp*pmf(1)       -pmf(2)
            vir(0)        = 2.0_wp*vir(1)       -vir(2)
            pmf(rdf%max_grid+1) = 2.0_wp*pmf(rdf%max_grid)  -pmf(rdf%max_grid-1)
            vir(rdf%max_grid+1) = 2.0_wp*vir(rdf%max_grid)  -vir(rdf%max_grid-1)
            pmf(rdf%max_grid+2) = 2.0_wp*pmf(rdf%max_grid+1)-pmf(rdf%max_grid)
            vir(rdf%max_grid+2) = 2.0_wp*vir(rdf%max_grid+1)-vir(rdf%max_grid)

            ! resample using 3pt interpolation

            Do ig=1,ngrid
              rrr = Real(ig,wp)*dgrid
              ll = Int(rrr*rdlr)

              ! +0.5_wp due to half-a-bin shift in the original (bin-centered) grid

              coef = rrr*rdlr-Real(ll,wp)+0.5_wp

              fed0 = pmf(ll)
              fed1 = pmf(ll+1)
              fed2 = pmf(ll+2)

              t1 = fed0 + (fed1 - fed0)*coef
              t2 = fed1 + (fed2 - fed1)*(coef - 1.0_wp)

              fed = t1 + (t2-t1)*coef*0.5_wp

              dfed0 = vir(ll)
              dfed1 = vir(ll+1)
              dfed2 = vir(ll+2)

              t1 = dfed0 + (dfed1 - dfed0)*coef
              t2 = dfed1 + (dfed2 - dfed1)*(coef - 1.0_wp)

              dfed = t1 + (t2-t1)*coef*0.5_wp

              If (comm%idnode == 0) &
                Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*rrr*rdlr
            End Do
          End If
        End Do
      End Do

      If (comm%idnode == 0) Then
        Close(Unit=npdgdt)
        Close(Unit=npdfdt)
      End If

      Deallocate (dstdrdf,pmf,vir, Stat = fail)
      If (fail > 0) Then
        Write(message,'(a)') 'rdf_compute - deallocation failure'
        Call error(0,message)
      End If
    End If

  End Subroutine rdf_compute

  Subroutine calculate_block(temp,rcut,neigh,sites,config,rdf)
    Real( Kind = wp ), Intent(in)            :: temp, rcut
    Type( neighbours_type), Intent( In    ) :: neigh
    Type( site_type ), Intent( In    ) :: sites
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( configuration_type ), Intent( InOut ) :: config

    Real( Kind = wp ), Dimension( 1:neigh%max_list ) :: rrt, xxt, yyt, zzt
    Real( Kind = wp )                        :: kT2engo, delr, rdlr, dgrid, pdfzero, factor1, rrr,dvol,gofr,gofr1

    Integer :: i, loopend, j, ia, ib, ngrid, kk, k, limit, jj
    Logical :: zero

    kT2engo = boltz*temp/engunit
    ! grid interval for rdf%rdf tables
    delr = rcut/Real(rdf%max_grid,wp)
    rdlr = 1.0_wp/delr
    ! resampling grid and grid interval for rdf%rdf tables
    ngrid = Max(Nint(rcut/delr_max),rdf%max_grid)
    dgrid = rcut/Real(ngrid,wp)
    pdfzero = 1.0e-9_wp
    Do ia=1,sites%ntype_atom
      Do ib=ia,sites%ntype_atom

        If (ia == ib .and. sites%num_type(ia) < 2) Cycle
        ! number of the interaction by its rdf%rdf key
        kk=rdf%list(ib*(ib-1)/2+ia)
        ! only for valid interactions specified for a look up
        ! global sum of data on all nodes
        ! normalisation factor
        factor1=config%volm*sites%dens(ia)*sites%dens(ib)*Real(rdf%n_configs,wp)
        If (ia == ib) factor1=factor1*0.5_wp*(1.0_wp-1.0_wp/sites%num_type(ia))
        ! loop over distances
        zero=.true.
        Do i=1,rdf%max_grid
          If (zero .and. i < (rdf%max_grid-3)) zero=(rdf%tmp_rdf(i+2,kk, rdf%block_number) <= 0.0_wp)
          gofr= rdf%tmp_rdf(i,kk, rdf%block_number)/factor1
          rrr = (Real(i,wp)-0.5_wp)*delr
          dvol= fourpi*delr*(rrr**2+delr**2/12.0_wp)
          gofr= gofr/dvol
          ! zero it if < pdfzero
          If (gofr < pdfzero) Then
            gofr1 = 0.0_wp
          Else
            gofr1 = gofr
          End If
          ! Store information to compute block average
          rdf%block_averages(ia, ib,i, rdf%block_number) = rdf%block_averages(ia,ib,i, rdf%block_number) + gofr1  
        End Do
      End Do
    End Do

  End Subroutine calculate_block

  Subroutine calculate_errors(temp, rcut, num_steps, neigh, sites, rdf, config, comm)

    Real( Kind = wp ), Intent( In )                      :: temp, rcut
    Type( neighbours_type ), Intent( In    ) :: neigh
    Type( site_type ), Intent( In    ) :: sites
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type(comms_type), Intent( InOut )                    :: comm
    Type( configuration_type ), Intent( InOut )          :: config

    Real( Kind = wp )                                    :: test1, delr
    Real( kind = wp ), Dimension( :, : , :), Allocatable :: averages, errors
    Real( kind = wp)                                     :: i_nr_blocks, s

    Integer, Intent(In) :: num_steps
    Integer             :: nr_blocks, i, j,k ,l, ierr(2), a, b, ia, ib, kk

    Character ( Len = 256 )  :: messages(2)

    test1 = 0.0_wp
    rdf%block_number = 1
    ierr = 0

    If(comm%mxnode > 1 .and. (.not. rdf%tmp_rdf_sync)) Then
      Do i=1, rdf%num_blocks+1
        Call gsum(comm,rdf%tmp_rdf(:,:,i))
      End Do
      rdf%tmp_rdf_sync = .True.
    End If

    Allocate(averages(sites%ntype_atom,sites%ntype_atom, rdf%max_grid), stat = ierr(1))
    Allocate(errors(sites%ntype_atom,sites%ntype_atom, rdf%max_grid), stat = ierr(2))
    If(Any(ierr>0)) Then
      Call error(1084)
    End If
    averages = 0.0_wp
    errors = 0.0_wp

    !Compute the rdf for each of the blocks
    Do nr_blocks=1, rdf%num_blocks+1
      Call calculate_block(temp, rcut,neigh,sites,config,rdf)
    End Do
    rdf%block_number = nr_blocks

    !Compute the errors.
    i_nr_blocks = 1.0_wp / Real(nr_blocks, wp)
    Do k=1, nr_blocks
      Do l=1, rdf%max_grid
        Do j=1, sites%ntype_atom
          Do i=1, sites%ntype_atom
            averages(i,j,l) = averages(i,j,l) + rdf%block_averages(i,j,l,k) * i_nr_blocks
          End Do
        End Do
      End Do
    End Do

    i_nr_blocks = 1.0_wp / Real(nr_blocks *(nr_blocks-1), wp)
    Do i=1, nr_blocks
      Do k=1, sites%ntype_atom
        Do j=1, sites%ntype_atom
          Do l=1, rdf%max_grid
            errors(j,k,l) = errors(j,k,l) + ( (rdf%block_averages(j,k,l,i) - averages(j,k,l))**2 * i_nr_blocks )
          End Do
        End Do
      End Do
    End Do

    Do l=1, rdf%max_grid
      Do j=1, sites%ntype_atom
        Do i = 1, sites%ntype_atom
          averages(i,j,l) = averages(i,j,l) * Real(nr_blocks,wp)
        End Do
      End Do
    End Do

    !output errors
    If (comm%idnode == 0) Then
      Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
      Write(nrdfdt,'(a)') config%cfgname
      Write(nrdfdt,'(2i10)') rdf%n_pairs,rdf%max_grid

      delr = rcut/Real(rdf%max_grid,wp)
      Do j =1, sites%ntype_atom
        Do k = j, sites%ntype_atom
          If (j == k .and. sites%num_type(j) < 2) Cycle
          kk=rdf%list(k*(k-1)/2+j)
          If (kk > 0 .and. kk <= rdf%n_pairs) Then
            Write(messages(1),'(2x,a,2(1x,a8))') 'g(r): ',sites%unique_atom(j),sites%unique_atom(k)
            Write(messages(2),'(8x,a1,6x,a4,9x,a4)') 'r','g(r)','n(r)'
            Call info(messages,2,.true.)
            Write(nrdfdt,'(2a8)') sites%unique_atom(j),sites%unique_atom(k)
            Do i=1,rdf%max_grid
              Write(nrdfdt,"(1p,2e14.6,2e14.6)") ((Real(i,wp)-0.5_wp)*delr),averages(j,k,i),errors(j,k,i)
            End Do
          End If
        End Do
      End Do
      Close(Unit=nrdfdt)
    End If
    Deallocate(averages,stat=ierr(1))
    Deallocate(errors,stat=ierr(2))

    If (Any(ierr>0)) Then
      Call error(1084)
    End If
  End Subroutine calculate_errors

  Subroutine calculate_errors_jackknife(temp,rcut,num_steps,neigh,sites,rdf,config,comm)

    Real( Kind = wp ), Intent(In)                        :: temp, rcut
    Type( neighbours_type ), Intent( In    ) :: neigh
    Type( site_type ), Intent( In    ) :: sites
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type(comms_type), Intent( InOut )                    :: comm
    Type( configuration_type ), Intent( InOut )          :: config

    Real( Kind = wp )                                    :: test1
    Real( Kind = wp ), Dimension( :, : , :), Allocatable :: averages, errors
    Real(Kind = wp)                                      :: i_nr_blocks, delr, s

    Integer, Intent(in) :: num_steps
    Integer             :: nr_blocks, i, j,k ,l, ierr(2), a, b, kk

    Character( Len = 256 ) :: messages(2)

    test1 = 0.0_wp
    rdf%block_number = 1
    ierr = 0 
    If(comm%mxnode > 1 .and. .not. rdf%tmp_rdf_sync) Then
      Do i=1, rdf%num_blocks+1
        Call gsum(comm,rdf%tmp_rdf(:,:,i))
      End Do
      rdf%tmp_rdf_sync = .True.
    End If

    Allocate(averages(sites%ntype_atom,sites%ntype_atom, rdf%max_grid), stat = ierr(1))
    Allocate(errors(sites%ntype_atom,sites%ntype_atom, rdf%max_grid), stat = ierr(2))
    If( Any(ierr > 0)) Then
      Call error(1084)
    End If
    averages = 0.0_wp
    errors = 0.0_wp
    rdf%block_averages =0.0_wp

    !Compute the rdf%rdf for each of the blocks
    Do nr_blocks=1,rdf%num_blocks+1
      Call calculate_block(temp, rcut,neigh,sites,config,rdf)
    End Do
    rdf%block_number = nr_blocks
    i_nr_blocks = 1.0_wp / Real(nr_blocks, wp)

    Do k=1, nr_blocks
      Do l=1, rdf%max_grid
        Do j=1, sites%ntype_atom
          Do i=1, sites%ntype_atom
            averages(i,j,l) = averages(i,j,l) + rdf%block_averages(i,j,l,k) 
          End Do
        End Do
      End Do
    End Do


    i_nr_blocks = 1.0_wp / Real(nr_blocks-1, wp)
    !Create jackknife bins
    Do k=1, nr_blocks
      Do l=1, rdf%max_grid
        Do j=1, sites%ntype_atom
          Do i=1, sites%ntype_atom
            rdf%block_averages(i,j,l,k) = (averages(i,j,l) - rdf%block_averages(i,j,l,k)) * i_nr_blocks
          End Do
        End Do
      End Do
    End Do

    !Average
    i_nr_blocks = 1.0_wp / Real(nr_blocks,wp)
    Do l=1, rdf%max_grid
      Do j=1, sites%ntype_atom
        Do i=1, sites%ntype_atom
          averages(i,j,l) = averages(i,j,l) * i_nr_blocks
        End Do
      End Do
    End Do

    !Errors
    !Compute the errors
    i_nr_blocks = Real((nr_blocks-1), wp) / Real(nr_blocks, wp)
    Do i=1, nr_blocks
      Do k=1, sites%ntype_atom
        Do j=1, sites%ntype_atom
          Do l=1, rdf%max_grid
            errors(j,k,l) = errors(j,k,l) + ( (rdf%block_averages(j,k,l,i) - averages(j,k,l))**2 * i_nr_blocks )
          End Do
        End Do
      End Do
    End Do

    Do l=1, rdf%max_grid
      Do j=1, sites%ntype_atom
        Do i = 1, sites%ntype_atom
          averages(i,j,l) = averages(i,j,l)*Real(nr_blocks,wp)
        End Do
      End Do
    End Do

    !output errors
    If (comm%idnode == 0) Then
      Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
      Write(nrdfdt,'(a)') config%cfgname
      Write(nrdfdt,'(2i10)') rdf%n_pairs,rdf%max_grid

      delr = rcut/Real(rdf%max_grid,wp)
      Do j =1, sites%ntype_atom
        Do k = j, sites%ntype_atom
          If (j == k .and. sites%num_type(j) < 2) Cycle
          kk=rdf%list(k*(k-1)/2+j)
          If (kk > 0 .and. kk <= rdf%n_pairs) Then
            Write(messages(1),'(2x,a,2(1x,a8))') 'g(r): ',sites%unique_atom(j),sites%unique_atom(k)
            Write(messages(2),'(8x,a1,6x,a4,9x,a4)') 'r','g(r)','n(r)'
            Call info(messages,2,.true.)
            Write(nrdfdt,'(2a8)') sites%unique_atom(j),sites%unique_atom(k)
            Do i=1,rdf%max_grid
              Write(nrdfdt,"(1p,2e14.6,2e14.6)") ((Real(i,wp)-0.5_wp)*delr),averages(j,k,i),errors(j,k,i)
            End Do
          End If
        End Do
      End Do
      Close(Unit=nrdfdt)
    End If
    Deallocate(averages,stat=ierr(1))
    Deallocate(errors,stat=ierr(2))

    If (Any(ierr>0)) Then
      Call error(1084)
    End If
  End Subroutine calculate_errors_jackknife

  Subroutine rdf_excl_collect(iatm,rrt,neigh,config,rdf)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistic for radial
    ! distribution functions of excluded pairs
    !
    ! Note: to be used as part of two_body_forces
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2014
    ! contrib   - a.b.g.chalk january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( neighbours_type), Intent( In    ) :: neigh
    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: rrt
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( configuration_type ), Intent( InOut ) :: config

    Integer                 :: limit,idi,jatm,ai,aj,keyrdf,kk,ll,m
    Real( Kind = wp )       :: rdelr,rrr

    ! set cutoff condition for pair forces and grid interval for rdf%rdf tables
    rdelr= Real(rdf%max_grid,wp)/neigh%cutoff

    ! global identity and type of iatm
    idi=config%ltg(iatm)
    ai=config%ltype(iatm)

    ! Get neigh%list limit
    limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

    ! start of primary loop for rdf%rdf accumulation
    Do m=1,limit

      ! atomic and type indices
      jatm=neigh%list(neigh%list(0,iatm)+m,iatm)
      aj=config%ltype(jatm)
      If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
        ! rdf%rdf function indices
        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=rdf%list(keyrdf)
        ! only for valid interactions specified for a look up
        If (kk > 0 .and. kk <= rdf%n_pairs) Then
          ! apply truncation of potential
          rrr=rrt(m)
          If (rrr < neigh%cutoff) Then
            ll=Min(1+Int(rrr*rdelr),rdf%max_grid)
            ! accumulate correlation
            rdf%rdf(ll,kk) = rdf%rdf(ll,kk) + 1.0_wp
            If(rdf%l_errors_block .or. rdf%l_errors_jack) Then
              rdf%tmp_rdf(ll,kk,rdf%block_number) = rdf%tmp_rdf(ll,kk,rdf%block_number) + 1.0_wp
            End If
          End If
        End If
      End If

    End Do

  End Subroutine rdf_excl_collect

  Subroutine rdf_frzn_collect(iatm,rrt,neigh,config,rdf)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistic for radial
    ! distribution functions of frozen pairs
    !
    ! Note: to be used as part of two_body_forces
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2014
    ! contrib   - a.b.g.chalk january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( neighbours_type), Intent( In    ) :: neigh
    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: rrt
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( configuration_type ), Intent( InOut ) :: config

    Integer                 :: limit,idi,jatm,ai,aj,keyrdf,kk,ll,m
    Real( Kind = wp )       :: rdelr,rrr

    ! set cutoff condition for pair forces and grid interval for rdf%rdf tables
    rdelr= Real(rdf%max_grid,wp)/neigh%cutoff

    ! global identity and type of iatm
    idi=config%ltg(iatm)
    ai=config%ltype(iatm)

    ! Get neigh%list limit
    limit=neigh%list(-2,iatm)-neigh%list(-1,iatm)

    ! start of primary loop for rdf%rdf accumulation
    Do m=1,limit
      ! atomic and type indices
      jatm=neigh%list(neigh%list(-1,iatm)+m,iatm)
      aj=config%ltype(jatm)

      If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
        ! rdf%rdf function indices
        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=rdf%list(keyrdf)
        ! only for valid interactions specified for a look up
        If (kk > 0 .and. kk <= rdf%n_pairs) Then
          ! apply truncation of potential
          rrr=rrt(m)
          If (rrr < neigh%cutoff) Then
            ll=Min(1+Int(rrr*rdelr),rdf%max_grid)
            ! accumulate correlation
            rdf%rdf(ll,kk) = rdf%rdf(ll,kk) + 1.0_wp
            If(rdf%l_errors_block .or. rdf%l_errors_jack) Then
              rdf%tmp_rdf(ll,kk,rdf%block_number) = rdf%tmp_rdf(ll,kk,rdf%block_number) + 1.0_wp
            End If
          End If
        End If
      End If

    End Do

  End Subroutine rdf_frzn_collect

  Subroutine usr_collect(rrt,rdf)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistic for USR RDFs
    !
    ! Note: to be used in external_field_apply
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & a.brukhno november 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent( In    ) :: rrt
    Type( rdf_type ), Intent( InOut ) :: rdf

    Integer           :: ll
    Real( Kind = wp ) :: rdelr

    ! set cutoff condition for pair forces and grid interval for rdf%rdf tables

    rdelr= Real(rdf%max_grid_usr,wp)/rdf%cutoff_usr

    If (rrt < rdf%cutoff_usr) Then ! apply truncation of potential
      ll=Min(1+Int(rrt*rdelr),rdf%max_grid_usr)
      rdf%usr(ll) = rdf%usr(ll) + 1.0_wp ! accumulate correlation
      rdf%n_configs_usr = rdf%n_configs_usr + 1        ! Increment sample
    End If

  End Subroutine usr_collect

  Subroutine usr_compute(rdf,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating radial distribution function
    ! from accumulated data for umbrella sampled two COMs separation (ushr)
    !
    ! to be used in exernal_field_apply & statistics_result
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & a.brukhno november 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( comms_type ), Intent( InOut ) :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Integer           :: i
    Real( Kind = wp ) :: delr,rdlr,factor1,rrr,dvol,gofr,sum0,sum1

    ! grid interval for rdf%rdf tables

    delr = rdf%cutoff_usr/Real(rdf%max_grid_usr,wp)
    rdlr = 1.0_wp/delr

    ! open RDF file and Write headers

    If (comm%idnode == 0) Then
      Open(Unit=nrdfdt, File='USRDAT', Status='replace')
      Write(nrdfdt,'(2a)') '# '//config%cfgname
      Write(nrdfdt,'(a)')  "# RDF for the two fragments' COMs (umbrella sampling)"
      Write(nrdfdt,'(a,i10,f12.6,i10,e15.6,/)') '# bins, cutoff, frames, volume: ', &
        rdf%max_grid_usr,rdf%cutoff_usr,rdf%n_configs_usr,config%volm
      Write(nrdfdt,'(a)') '#'
    End If

    ! global sum of data on all nodes

    Call gsum(comm,rdf%usr(1:rdf%max_grid_usr))

    ! get normalisation factor

    factor1 = Sum(rdf%usr(1:rdf%max_grid_usr))

    ! running integration of rdf%rdf

    sum0 = 0.0_wp
    sum1 = 0.0_wp

    ! loop over distances

    Do i=1,rdf%max_grid_usr
      gofr = rdf%usr(i)/factor1
      sum0 = sum0 + gofr

      rrr  = (Real(i,wp)-0.5_wp)*delr
      dvol = fourpi*delr*(rrr**2+delr**2/12.0_wp)
      gofr = gofr*config%volm/dvol
      sum1 = sum1 + gofr

      ! print out information

      If (comm%idnode == 0) Write(nrdfdt,"(1p,4e15.6)") rrr,gofr
      !     If (idnode == 0) Write(nrdfdt,"(1p,4e15.6)") rrr,gofr,sum0,sum1
    End Do

    If (comm%idnode == 0) Close(Unit=nrdfdt)

    ! distribute rdf%usr between nodes

    rdf%usr(:) = rdf%usr(:) / Real(comm%mxnode,wp)

  End Subroutine usr_compute

  !> Increase block number when required
  Subroutine rdf_increase_block_number(rdf,nstep)
    Type( rdf_type ), Intent( InOut ) :: rdf
    Integer( Kind = wi ), Intent( In    ) :: nstep

    If ((rdf%l_errors_block .or. rdf%l_errors_jack) .and. &
      mod(nstep, rdf%block_size) == 0) Then
      rdf%block_number = rdf%block_number + 1
    End If
  End Subroutine rdf_increase_block_number

  Subroutine cleanup(T)
    Type( rdf_type ) :: T

    If (Allocated(T%list)) Then
      Deallocate(T%list)
    End If
    If (Allocated(T%rdf)) Then
      Deallocate(T%rdf)
    End If
    If (Allocated(T%usr)) Then
      Deallocate(T%usr)
    End If
  End Subroutine cleanup
End Module rdfs
