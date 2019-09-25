Module two_body
  Use kinds, Only : wp
  Use comms,   Only : comms_type,gsum
  Use constants, Only : pi, r4pie0
  Use site, Only : site_type
  Use configuration,  Only : configuration_type
  Use neighbours,     Only : neighbours_type,link_cell_pairs
  Use spme,  Only : init_spme_data, spme_self_interaction
  Use ewald,           Only : ewald_type, ewald_type, ewald_vdw_init, ewald_vdw_count, ewald_vdw_coeffs
  Use mpole,          Only : mpole_type,POLARISATION_CHARMM
  Use coul_spole,     Only : coul_fscp_forces, coul_rfp_forces, coul_cp_forces, coul_dddp_forces
  Use coul_mpole,    Only : coul_fscp_mforces, coul_rfp_mforces, coul_cp_mforces, &
    coul_dddp_mforces, coul_chrm_forces, d_ene_trq_mpoles
  Use poisson, Only : poisson_type,poisson_forces,poisson_excl_forces,poisson_frzn_forces
  Use vdw,     Only : vdw_type,vdw_forces, &
    VDW_NULL, VDW_TAB, VDW_12_6, VDW_LENNARD_JONES, VDW_N_M, &
    VDW_BUCKINGHAM, VDW_BORN_HUGGINS_MEYER, VDW_HYDROGEN_BOND, &
    VDW_N_M_SHIFT, VDW_MORSE, VDW_WCA, VDW_DPD, VDW_AMOEBA, &
    VDW_LENNARD_JONES_COHESIVE, VDW_MORSE_12, VDW_RYDBERG, VDW_ZBL, &
    VDW_ZBL_SWITCH_MORSE, VDW_ZBL_SWITCH_BUCKINGHAM
  Use numerics, Only : calc_erfc, calc_erfc_deriv
  Use metal,   Only : metal_type,metal_forces,metal_ld_compute,metal_lrc
  Use kim,     Only : kim_type,kim_energy_and_forces
  Use rdfs,    Only : rdf_type,rdf_collect,rdf_excl_collect,rdf_frzn_collect, &
    rdf_increase_block_number
  Use errors_warnings, Only : error,error_alloc,error_dealloc
  Use ewald_spole, Only : ewald_spme_forces, ewald_real_forces_gen, ewald_frzn_forces, ewald_excl_forces, &
    & ewald_real_forces_coul
  ! Use ewald_mpole, Only : ewald_spme_mforces, ewald_real_mforces,ewald_frzn_mforces,ewald_excl_mforces, &
  !                        ewald_spme_mforces_d,ewald_real_mforces_d,ewald_excl_mforces_d,ewald_excl_mforces

  Use timer,      Only : timer_type, start_timer, stop_timer
  Use development, Only : development_type
  Use statistics, Only : stats_type
  Use core_shell, Only : core_shell_type
  Use electrostatic, Only : electrostatic_type, &
    ELECTROSTATIC_EWALD,ELECTROSTATIC_DDDP, &
    ELECTROSTATIC_COULOMB,ELECTROSTATIC_COULOMB_FORCE_SHIFT, &
    ELECTROSTATIC_COULOMB_REACTION_FIELD,ELECTROSTATIC_POISSON
  Use domains, Only : domains_type
  Implicit None

  Private

  Public :: two_body_forces
Contains

  Subroutine two_body_forces(ensemble,    &
    lbook,megfrz, &
    leql,nsteql,nstep,         &
    cshell,stats,ewld,devel,met,pois,neigh,sites,vdws,rdf,mpoles,electro, &
    domain,tmr,kim_data,config,comm)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating interatomic forces and rdf%rdf
    !! using the verlet neighbour list
    !!
    !! vdws%n_vdw > 0 ------ switch for vdw potentials calculation
    !! met%n_potentials > 0 ------ switch for metal local density and potentials
    !!                   calculations
    !!
    !! nstfce - the rate at which the k-space contributions of SPME are
    !!          refreshed.  Once every 1 <= nstfce <= 7 steps.
    !!
    !! copyright - daresbury laboratory
    !! author    - i.t.todorov february 2017
    !! contrib   - h.a.boateng february 2016
    !! contrib   - p.s.petkov february 2015
    !! contrib   - a.b.g.chalk january 2017
    !! contrig   - j.s. wilkins october 2018
    !! refactoring:
    !!           - a.m.elena march-october 2018
    !!           - j.madge march-october 2018
    !!           - a.b.g.chalk march-october 2018
    !!           - i.scivetti march-october 2018
    !!
    !!----------------------------------------------------------------------!



    Logical,                                  Intent( In    ) :: lbook,leql
    Integer,                                  Intent( In    ) :: ensemble,        &
      megfrz, &
      nsteql,nstep
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( stats_type ), Intent( InOut )                       :: stats
    Type( ewald_type ),                       Intent( InOut ) :: ewld
    Type( development_type ),                 Intent( In    ) :: devel
    Type( metal_type ),                       Intent( InOut ) :: met
    Type( poisson_type ),                     Intent( InOut ) :: pois
    Type( neighbours_type ),                  Intent( InOut ) :: neigh
    Type( site_type ),                        Intent( In    ) :: sites
    Type( vdw_type ),                         Intent( InOut ) :: vdws
    Type( rdf_type ),                         Intent( InOut ) :: rdf
    Type( mpole_type ),                       Intent( InOut ) :: mpoles
    Type( timer_type ),                       Intent( InOut ) :: tmr
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( kim_type ), Intent( InOut ) :: kim_data
    Type( configuration_type ),               Intent( InOut ) :: config
    Type( comms_type ),                       Intent( InOut ) :: comm


    Real( Kind = wp ) :: factor_nz

    Logical           :: safe, l_do_rdf
    Integer           :: fail,i,j,k,limit
    Real( Kind = wp ) :: engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl, &
      engcpe_ch,vircpe_ch,engcpe_ex,vircpe_ex, &
      engcpe_fr,vircpe_fr,engcpe_nz,vircpe_nz, &
      vircpe_dt,                              &
      engden,virden,engmet,virmet,             &
      engvdw,virvdw,engkim,virkim,             &
      engacc,viracc,tmp,buffer(0:23)

    Real( kind = wp ) :: engvdw_rc, engvdw_rl
    Real( kind = wp ) :: virvdw_rc, virvdw_rl

    ! Array to remap charge/disp to 2D array
    Real( kind = wp ), dimension( : ), allocatable :: coul_coeffs

    Real( kind = wp ), dimension( :,: ), allocatable :: vdw_coeffs

    Integer, dimension(:), allocatable, save :: reduced_VdW
    Integer :: ipot, ivdw, keypot
    Logical :: skip

    Logical, save :: newjob = .true.

    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

    call start_timer('Two-Body Init')

    safe = .True.
    fail=0

    Allocate (xxt(1:neigh%max_list),yyt(1:neigh%max_list),zzt(1:neigh%max_list),rrt(1:neigh%max_list), Stat=fail)
    If (fail > 0) call error_alloc('distance arrays','two_body_forces')

    l_do_rdf = (rdf%l_collect .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,rdf%freq) == 0)

    ! If k-space SPME is evaluated infrequently check whether
    ! at this timestep to evaluate or "refresh" with old values.
    ! At restart allocate the "refresh" arrays and force a fresh
    ! evaluation.  Repeat the same but only for the SPME k-space
    ! frozen-frozen evaluations in constant volume ensembles only.

    ! Coulomb
    if (ewld%vdw .and. .not. ewld%active) call error(0,'Ewald VdW requested but ewald not enabled')
    if (ewld%active) then
      allocate(coul_coeffs(config%mxatms), stat=fail)
      if (fail>0) call error_alloc('coul_coeffs','two_body_forces')
      coul_coeffs = config%parts(:)%chge

      call electro%erfc%init(ewld%alpha*neigh%cutoff, calc_erfc)
      call electro%erfc_deriv%init(ewld%alpha*neigh%cutoff, calc_erfc_deriv)

      write(0, *) electro%erfc%calc(0.1_wp), calc_erfc(0.1_wp)
      write(0, *) electro%erfc%calc(3.1_wp), calc_erfc(3.1_wp)
      write(0, *) electro%erfc%calc(6.1_wp), calc_erfc(6.1_wp)

      if ( newjob ) then

          ! Assume no VDW
          ewld%num_pots = 0

          if ( ewld%vdw ) call ewald_vdw_count(ewld, vdws, reduced_VdW)

          allocate ( ewld%spme_data(0:ewld%num_pots), stat=fail )
          if ( fail > 0 ) call error_alloc('ewld%spme_data','two_body_forces')

          if ( electro%key == ELECTROSTATIC_EWALD ) then
            ewld%spme_data(0)%scaling = r4pie0/electro%eps
            call init_spme_data( ewld%spme_data(0), 1 )

            if ( electro%multipolar ) then
              call spme_self_interaction(ewld%alpha, config%natms, coul_coeffs, comm, ewld%spme_data(0), electro%mpoles)
            else
              call spme_self_interaction(ewld%alpha, config%natms, coul_coeffs, comm, ewld%spme_data(0))
            end if

          else
            allocate ( ewld%spme_data(1:ewld%num_pots), stat=fail )
            if ( fail > 0 ) call error_alloc('ewld%spme_data','two_body_forces')
          end if

          if ( ewld%vdw ) then
            call ewald_vdw_init(ewld, vdws, reduced_VdW)
            call ewald_vdw_coeffs(config, vdws, ewld, reduced_vdw, vdw_coeffs)

            do ipot = 1, ewld%num_pots
              call spme_self_interaction(ewld%alpha, config%natms, vdw_coeffs(:,ipot), comm, ewld%spme_data(ipot))
            end do

            ! Disable long-range corrections (calculating long range explicitly)
            vdws%elrc = 0.0_wp
            vdws%vlrc = 0.0_wp

          end if

        newjob = .false.

      end if


    end if

    ! initialise energy and virial accumulators

    engkim = 0.0_wp
    virkim = 0.0_wp

    engden    = 0.0_wp
    virden    = 0.0_wp

    engmet    = 0.0_wp
    virmet    = 0.0_wp

    engvdw    = 0.0_wp
    virvdw    = 0.0_wp
    engvdw_rc = 0.0_wp
    virvdw_rc = 0.0_wp
    engvdw_rl = 0.0_wp
    virvdw_rl = 0.0_wp

    stats%engsrp    = 0.0_wp
    stats%virsrp    = 0.0_wp


    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    engcpe_rl = 0.0_wp
    vircpe_rl = 0.0_wp

    engcpe_ch = 0.0_wp
    vircpe_ch = 0.0_wp

    engcpe_ex = 0.0_wp
    vircpe_ex = 0.0_wp

    engcpe_fr = 0.0_wp
    vircpe_fr = 0.0_wp

    engcpe_nz = 0.0_wp
    vircpe_nz = 0.0_wp

    vircpe_dt = 0.0_wp

    stats%engcpe    = 0.0_wp
    stats%vircpe    = 0.0_wp

#ifdef CHRONO
    Call stop_timer('Two-Body Init')
#endif

    ! Set up non-bonded interaction (verlet) list using link cells
    If (neigh%update) Then
      Call link_cell_pairs(vdws%cutoff,met%rcut,lbook,megfrz,cshell,devel, &
        neigh,mpoles,domain,tmr,config,comm)
    End If

    ! Calculate all contributions from KIM
    If (kim_data%active) Then
#ifdef CHRONO
    Call start_timer('KIM')
#endif
      Call kim_energy_and_forces(kim_data,config%natms,config%nlast,config%parts, &
        neigh%list,domain%map,config%lsite,config%lsi,config%lsa,config%ltg, &
        sites%site_name,engkim,virkim,stats%stress,comm)
#ifdef CHRONO
      Call stop_timer('KIM')
#endif
    End If

    If (met%n_potentials > 0) Then

      ! Reset metal long-range corrections (constant pressure/stress only)

      If (ensemble >= 20) Call metal_lrc(met,sites,config,comm)

      ! calculate local density in metals

      Call metal_ld_compute(engden,virden,stats%stress,sites%ntype_atom,met,neigh, &
        domain,config,comm)
    End If

    ! calculate coulombic forces, Ewald sum - fourier contribution
#ifdef CHRONO
    Call start_timer('Long Range')
#endif

    If (electro%key == ELECTROSTATIC_EWALD) Then

        call ewald_spme_forces(ewld,ewld%spme_data(0),electro,domain,config,comm, &
          & coul_coeffs,nstep,stats,engcpe_rc,vircpe_rc)

    End If

#ifdef CHRONO
    Call stop_timer('Long Range')
#endif

    if (ewld%vdw) then
      call start_timer('SPME Order-n')

      do ipot = 1, ewld%num_pots

        call start_timer('Long Range')
        call ewald_spme_forces(ewld,ewld%spme_data(ipot),electro,domain,config,comm, &
          & vdw_coeffs(:,ipot),nstep,stats,engacc,viracc)
        call stop_timer('Long Range')

        engvdw_rc=engvdw_rc+engacc
        virvdw_rc=virvdw_rc+viracc


      end do

      call stop_timer('SPME Order-n')
    end if

#ifdef CHRONO
    Call start_timer('Short Range')
#endif
    ! outer loop over atoms

    Do i=1,config%natms

      ! Get neigh%list limit

      limit=neigh%list(0,i)

      ! calculate interatomic distances

      Do k=1,limit
        j=neigh%list(k,i)

        xxt(k)=config%parts(i)%xxx-config%parts(j)%xxx
        yyt(k)=config%parts(i)%yyy-config%parts(j)%yyy
        zzt(k)=config%parts(i)%zzz-config%parts(j)%zzz
      End Do

      ! periodic boundary conditions not needed by LC construction
      !
      !     Call images(imcon,cell,limit,xxt,yyt,zzt)

      ! distances, thanks to Alin Elena (one too many changes)

      Do k=1,limit
        rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
      End Do

      ! calculate metal forces and potential

      If (met%n_potentials > 0) Then
        Call metal_forces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,safe,sites%ntype_atom,met,neigh,config)

        engmet=engmet+engacc
        virmet=virmet+viracc
      End If

      ! calculate short-range force and potential terms

      If (vdws%n_vdw > 0) Then

        if ( .not. ewld%vdw ) then

          Call vdw_forces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats,neigh,vdws,config)
          engvdw=engvdw+engacc
          virvdw=virvdw+viracc

        else
          do ipot = 1, ewld%num_pots

            Call ewald_real_forces_gen(ewld%alpha,ewld%spme_data(ipot),neigh,config,stats, &
              & vdw_coeffs(:,ipot),i,xxt,yyt,zzt,rrt,engacc,viracc)

            engvdw_rl=engvdw_rl+engacc
            virvdw_rl=virvdw_rl+viracc

          end do
        end if

      End If

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! COULOMBIC CONTRIBUTIONS
      !!!!!!!!!!!!!!!!!!!!!!!!!

      If (mpoles%max_mpoles > 0) Then

!!! MULTIPOLAR ATOMIC SITES

        If (electro%key == ELECTROSTATIC_EWALD) Then

          ! calculate coulombic forces, Ewald sum - real space contribution

          Call ewald_real_forces_coul(electro, ewld%alpha,ewld%spme_data(0),neigh,config,stats, &
            & i,xxt,yyt,zzt,rrt,engacc,viracc)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_DDDP) Then

          ! distance dependant dielectric potential

          Call coul_dddp_mforces(i,electro%eps,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,mpoles,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_COULOMB) Then

          ! coulombic 1/r potential with no truncation or damping

          Call coul_cp_mforces(i,electro%eps,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,mpoles,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then

          ! force-shifted coulomb potentials

          Call coul_fscp_mforces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,mpoles,electro,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then

          ! reaction field potential

          Call coul_rfp_mforces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,mpoles,electro,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        End If

      Else

        If (electro%key == ELECTROSTATIC_EWALD) Then

          ! calculate coulombic forces, Ewald sum - real space contribution

          Call ewald_real_forces_coul(electro, ewld%alpha,ewld%spme_data(0),neigh,config,stats, &
            & i,xxt,yyt,zzt,rrt,engacc,viracc)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_DDDP) Then

          ! distance dependant dielectric potential

          Call coul_dddp_forces(i,electro%eps,xxt,yyt,zzt,rrt,engacc,viracc,stats,neigh,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_COULOMB) Then

          ! coulombic 1/r potential with no truncation or damping

          Call coul_cp_forces(i,electro%eps,xxt,yyt,zzt,rrt,engacc,viracc,stats,neigh,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then

          ! force-shifted coulomb potentials

          Call coul_fscp_forces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats,neigh,electro,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        Else If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then

          ! reaction field potential

          Call coul_rfp_forces(i,xxt,yyt,zzt,rrt,engacc,viracc,stats,neigh,electro,config)

          engcpe_rl=engcpe_rl+engacc
          vircpe_rl=vircpe_rl+viracc

        End If

      End If

      ! accumulate radial distribution functions

      If (l_do_rdf) Call rdf_collect(i,rrt,neigh,config,rdf)

    End Do

    ! Poisson solver alternative to Ewald

    If (electro%key == ELECTROSTATIC_POISSON) Then
      Call poisson_forces(engacc,viracc,stats%stress,pois,electro,domain,config,comm)
      engcpe_rl=engcpe_rl+engacc
      vircpe_rl=vircpe_rl+viracc
    End If

    ! metal potential safety

    If (safe) Then
      tmp=0.0_wp
    Else
      tmp=1.0_wp
    End If

    ! in the case of bonded interactions 3 possible subcases
    ! cases for excluded interactions:
    ! (1) RDF accumulate further the short-range exclusions
    ! (2) Ewald corrections due to short-range exclusions
    ! (3) CHARMM core-shell self-induction additions

    If ( lbook .and. &
      (l_do_rdf .or. (Any([ELECTROSTATIC_EWALD,ELECTROSTATIC_POISSON] == electro%key)) &
      .or. mpoles%key == POLARISATION_CHARMM) ) Then
      Do i=1,config%natms ! outer loop over atoms
        limit=neigh%list(-1,i)-neigh%list(0,i) ! Get neigh%list limit
        If (limit > 0) Then

          ! calculate interatomic distances

          Do k=1,limit
            j=neigh%list(neigh%list(0,i)+k,i)

            xxt(k)=config%parts(i)%xxx-config%parts(j)%xxx
            yyt(k)=config%parts(i)%yyy-config%parts(j)%yyy
            zzt(k)=config%parts(i)%zzz-config%parts(j)%zzz
          End Do

          ! periodic boundary conditions not needed by LC construction
          !
          !           Call images(imcon,cell,limit,xxt,yyt,zzt)

          ! square of distances

          Do k=1,limit
            rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
          End Do

          ! accumulate radial distribution functions

          If (l_do_rdf) Call rdf_excl_collect(i,rrt,neigh,config,rdf)

          If (electro%key == ELECTROSTATIC_EWALD) Then ! Ewald corrections

            Call ewald_excl_forces(ewld,ewld%spme_data(0),neigh,electro,config,coul_coeffs, &
              i,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress)

            engcpe_ex=engcpe_ex+engacc
            vircpe_ex=vircpe_ex+viracc
          End If

          ! get CHARMM core-shell self-induction contributions

          If (mpoles%key == POLARISATION_CHARMM) Then
            If (neigh%list(-3,i)-neigh%list(0,i) > 0) Then
              Call coul_chrm_forces(i,electro%eps,xxt,yyt,zzt,rrt,engacc,viracc,stats%stress,neigh,mpoles,config)

              engcpe_ch=engcpe_ch+engacc
              vircpe_ch=vircpe_ch+viracc
            End If
          End If

        End If
      End Do
    End If

    ! counter for rdf%rdf statistics outside loop structures
    ! and frozen-frozen rdf%rdf completeness

    If (l_do_rdf) Then
      If (megfrz /= 0) Then

        ! outer loop over atoms

        Do i=1,config%natms

          ! Get neigh%list limit

          limit=neigh%list(-2,i)-neigh%list(-1,i)
          If (limit > 0) Then

            ! calculate interatomic distances

            Do k=1,limit
              j=neigh%list(neigh%list(-1,i)+k,i)

              xxt(k)=config%parts(i)%xxx-config%parts(j)%xxx
              yyt(k)=config%parts(i)%yyy-config%parts(j)%yyy
              zzt(k)=config%parts(i)%zzz-config%parts(j)%zzz
            End Do

            ! periodic boundary conditions not needed by LC construction
            !
            !              Call images(imcon,cell,limit,xxt,yyt,zzt)

            ! square of distances

            Do k=1,limit
              rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
            End Do

            ! accumulate radial distribution functions

            Call rdf_frzn_collect(i,rrt,neigh,config,rdf)
          End If

        End Do

      End If

      rdf%n_configs = rdf%n_configs + 1
    End If

#ifdef CHRONO
    Call stop_timer('Short Range')
    Call start_timer('Two-Body Final')
#endif

    Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
    If (fail > 0) call error_dealloc('distance arrays','two_body_forces')

    !Increase rdf%block_number when required
    If (l_do_rdf) Then
      Call rdf_increase_block_number(rdf,nstep)
    End If

    ! Further Ewald/Poisson Solver corrections or an infrequent refresh

    If (any([ELECTROSTATIC_EWALD,ELECTROSTATIC_POISSON] == electro%key)) Then
      If (ELECTROSTATIC_EWALD == electro%key) Then
        deallocate(coul_coeffs,stat=fail)
        if (fail>0) call error_dealloc('coul_coeffs','two_body_forces')
        if ( ewld%vdw ) then
          deallocate(vdw_coeffs,stat=fail)
          if (fail>0) call error_dealloc('vdw_coeffs','two_body_forces')
        end if
      end If

      ! non-zero total system charge correction (for the whole system)
      ! ( Fuchs, Proc. R. Soc., A, 151, (585),1935 )
      If (Abs(config%sumchg) > 1.0e-6_wp) Then
        factor_nz = -0.5_wp * (pi*r4pie0/electro%eps) * (config%sumchg/electro%damping)**2

        engcpe_nz=factor_nz/config%volm
        vircpe_nz=-3.0_wp*engcpe_nz
      End If
    End If

    ! Find the change of energy produced by the torques on multipoles
    ! under infinitesimal rotations & convert to Cartesian coordinates

    If (mpoles%max_mpoles > 0) Then
      Call d_ene_trq_mpoles(vircpe_dt,stats%stress,mpoles,config)
    End If

    ! sum up contributions to domain,potentials


    buffer( 0) = tmp
    buffer( 1) = engkim
    buffer( 2) = virkim
    buffer( 3) = engden
    buffer( 4) = virden
    buffer( 5) = engmet
    buffer( 6) = virmet
    buffer( 7) = engvdw
    buffer( 8) = virvdw
    buffer( 9) = engcpe_rc
    buffer(10) = vircpe_rc
    buffer(11) = engcpe_rl
    buffer(12) = vircpe_rl
    buffer(13) = engcpe_ch
    buffer(14) = vircpe_ch
    buffer(15) = engcpe_ex
    buffer(16) = vircpe_ex
    buffer(17) = engcpe_fr
    buffer(18) = vircpe_fr
    buffer(19) = vircpe_dt
    buffer(20) = engvdw_rl
    buffer(21) = engvdw_rc
    buffer(22) = virvdw_rl
    buffer(23) = virvdw_rc

    Call gsum(comm,buffer(0:23))

    tmp       = buffer( 0)
    engkim    = buffer( 1)
    virkim    = buffer( 2)
    engden    = buffer( 3)
    virden    = buffer( 4)
    engmet    = buffer( 5)
    virmet    = buffer( 6)
    engvdw    = buffer( 7)
    virvdw    = buffer( 8)
    engcpe_rc = buffer( 9)
    vircpe_rc = buffer(10)
    engcpe_rl = buffer(11)
    vircpe_rl = buffer(12)
    engcpe_ch = buffer(13)
    vircpe_ch = buffer(14)
    engcpe_ex = buffer(15)
    vircpe_ex = buffer(16)
    engcpe_fr = buffer(17)
    vircpe_fr = buffer(18)
    vircpe_dt = buffer(19)
    engvdw_rl = buffer(20)
    engvdw_rc = buffer(21)
    virvdw_rl = buffer(22)
    virvdw_rc = buffer(23)


    safe=(tmp < 0.5_wp)
    If (.not.safe) Call error(505)

    ! Self-interaction is constant for the default charges only SPME

    ! If (electro%key == ELECTROSTATIC_EWALD) Then ! Sum it up for multipolar SPME
    !    If (mpoles%max_mpoles > 0 .and. mpoles%max_order <= 2) Call gsum(comm,ewld%spme_data(0)%self_interaction)
    !Write(message,'(a,1p,e18.10)') 'Self-interaction term: ',engsic
    !Call info(message,.true.)
    ! End If

    ! Globalise coulombic contributions: cpe

    stats%engcpe = engcpe_rc + engcpe_rl + engcpe_ch + engcpe_ex + engcpe_fr + engcpe_nz
    stats%vircpe = vircpe_rc + vircpe_rl + vircpe_ch + vircpe_ex + vircpe_fr + vircpe_nz + vircpe_dt

    ! Add non-zero total system charge correction to
    ! diagonal terms of stress tensor (per node)

    tmp = - vircpe_nz/(3.0_wp*Real(comm%mxnode,wp))
    stats%stress(1) = stats%stress(1) + tmp
    stats%stress(5) = stats%stress(5) + tmp
    stats%stress(9) = stats%stress(9) + tmp

    ! Globalise short-range, KIM and metal interactions with
    ! their long-range corrections contributions: srp

    stats%engsrp = engkim + (engden + engmet + met%elrc(0)) + (engvdw + vdws%elrc)
    stats%virsrp = virkim + (virden + virmet + met%vlrc(0)) + (virvdw + vdws%vlrc)

    ! Add long-range corrections to diagonal terms of stress tensor (per node)

    tmp = - (vdws%vlrc+met%vlrc(0))/(3.0_wp*Real(comm%mxnode,wp))
    stats%stress(1) = stats%stress(1) + tmp
    stats%stress(5) = stats%stress(5) + tmp
    stats%stress(9) = stats%stress(9) + tmp

#ifdef CHRONO
    Call stop_timer('Two-Body Final')
#endif

  End Subroutine two_body_forces
End Module two_body
