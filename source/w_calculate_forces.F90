!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

! for new simulations when using the relaxed shell model
! set shells on top of their cores preventatively

     If ( (cshell%megshl > 0 .and. cshell%keyshl == SHELL_RELAXED) .and. &
          (flow%restart_key == 0 .and. flow%step == 0 .and. flow%equil_steps > 0) ) Then
        Call core_shell_on_top(cshell,config,comm)

! Refresh mappings

        Call w_refresh_mappings(flow,cshell,cons,pmf,stat,msd_data,bond,angle, &
          dihedral,inversion,tether,neigh,sites,mpoles,rigid,domain,kim_data)
     End If

100  Continue ! Only used when relaxed is false

! Initialise force arrays and possible torques for multipolar electrostatics
! and stress tensor (these are all additive in the force subroutines)
     Do i=1,cnfig%mxatms
       config%parts(i)%fxx=0.0_wp
       config%parts(i)%fyy=0.0_wp
       config%parts(i)%fzz=0.0_wp
     End Do
     If (mpoles%max_mpoles > 0) Then
        mpoles%torque_x=0.0_wp ; mpoles%torque_y=0.0_wp ; mpoles%torque_z=0.0_wp
     End If
     stat%stress = 0.0_wp

! Calculate pair-like forces (metal,vdws,electrostatic) and add lrc

     If (.not.(met%max_metal == 0 .and. electro%key == ELECTROSTATIC_NULL .and. &
       l_n_v .and. rdf%max_rdf == 0) .or. kim_data%active) Then
       Call two_body_forces(thermo%ensemble,flow%book,megfrz, &
         flow%equilibration,flow%equil_steps,flow%step,cshell,stat,ewld,devel,met,pois,neigh,sites,vdws,rdf, &
         mpoles,electro,domain,tmr,kim_data,config,comm)
     End If

! Calculate tersoff forces

     If (tersoffs%n_potential > 0) Call tersoff_forces(tersoffs,stat,neigh,domain,config,comm)

! Calculate three-body forces

     If (threebody%ntptbp > 0) Call three_body_forces(stat,threebody,neigh,domain,config,comm)

! Calculate four-body forces

     If (fourbody%n_potential > 0) Call four_body_forces(fourbody,stat,neigh,domain,config,comm)
     call start_timer(tmr%t_bonded)
! Calculate shell model forces

     If (cshell%megshl > 0) Call core_shell_forces(cshell,stat,config,comm)

! Calculate tethered atom forces

     If (tether%total > 0) Call tethers_forces(stat,tether,config,comm)

! Calculate bond forces

     If (bond%total > 0) Then
        ltmp = (bond%bin_pdf > 0 .and. &
          ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
          Mod(flow%step,flow%freq_bond) == 0)

        flow%isw = 1 + Merge(1,0,ltmp)
        Call bonds_forces(flow%isw,stat%engbnd,stat%virbnd,stat%stress,neigh%cutoff, &
          stat%engcpe,stat%vircpe,bond,mpoles,electro,config,comm)
     End If

! Calculate valence angle forces

     If (angle%total > 0) Then
        ltmp = (angle%bin_adf > 0 .and. &
          ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
          Mod(flow%step,flow%freq_angle) == 0)

        flow%isw = 1 + Merge(1,0,ltmp)
        Call angles_forces(flow%isw,stat%engang,stat%virang,stat%stress,angle,config,comm)
     End If

! Calculate dihedral forces

     If (dihedral%total > 0) Then
        ltmp = (dihedral%bin_adf > 0 .and. &
          ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) &
          .and. Mod(flow%step,flow%freq_dihedral) == 0)

        flow%isw = 1 + Merge(1,0,ltmp)
        Call dihedrals_forces(flow%isw,stat%engdih,stat%virdih,stat%stress, &
           neigh%cutoff,stat%engcpe,stat%vircpe,stat%engsrp, &
           stat%virsrp,dihedral,vdws,mpoles,electro,config,comm)
     End If

! Calculate inversion forces

     If (inversion%total > 0) Then
        ltmp = (inversion%bin_adf > 0 .and. &
          ((.not.flow%equilibration) .or. flow%step >= flow%equil_steps) .and. &
          Mod(flow%step,flow%freq_inversion) == 0)

        flow%isw = 1 + Merge(1,0,ltmp)
        Call inversions_forces(flow%isw,stat%enginv,stat%virinv,stat%stress,inversion,config,comm)
     End If
     call stop_timer(tmr%t_bonded)

! Apply external field

     If (ext_field%key /= FIELD_NULL) Then
       Call external_field_apply(flow%time,flow%equilibration,flow%equil_steps,flow%step,cshell,stat,rdf, &
         ext_field,rigid,domain,config,comm)
     End If

! Apply PLUMED driven dynamics

     If (plume%l_plumed) Then
        stat%stpcfg =stat%engcpe + stat%engsrp + stat%engter + stat%engtbp + stat%engfbp + &
                 stat%engshl + stat%engtet + stat%engfld +                   &
                 stat%engbnd + stat%engang + stat%engdih + stat%enginv
        Call plumed_apply(config,flow%run_steps,flow%step,stat,plume,comm)
     End If
! Apply pseudo thermostat - force cycle (0)

     If (thermo%l_stochastic_boundaries) Then
       Call stochastic_boundary_vv(0,thermo%tstep,flow%step,sites%dof_site,cshell,stat,thermo,rigid,domain,config,seed,comm)
     End If

! Cap forces in equilibration mode

     If (flow%step <= flow%equil_steps .and. flow%force_cap) Call cap_forces(fmax,thermo%temp,config,comm)

! Frozen atoms option

     Call freeze_atoms(config)

! Minimisation option and Relaxed shell model optimisation

     If (flow%simulation .and. (minim%minimise .or. cshell%keyshl == SHELL_RELAXED)) Then
        stat%stpcfg = stat%engcpe + stat%engsrp + stat%engter + stat%engtbp + stat%engfbp + &
                 stat%engshl + stat%engtet + stat%engfld +                   &
                 stat%engbnd + stat%engang + stat%engdih + stat%enginv

        If (cshell%keyshl == SHELL_RELAXED) Then
          Call core_shell_relax(flow%strict,rdf%l_collect, &
            stat%stpcfg,cshell,stat,domain,config,files,comm)
        End If

        If (.not.cshell%relaxed) Go To 200 ! Shells relaxation takes priority over minimisation

        If (minim%minimise .and. flow%step >= 0 .and. flow%step <= flow%run_steps .and. flow%step <= flow%equil_steps) Then
          If      (minim%freq == 0 .and. flow%step == 0) Then
            Call minimise_relax(flow%strict .or. cshell%keyshl == SHELL_RELAXED, &
              rdf%l_collect,megatm,thermo%tstep,stat%stpcfg,io,stat,pmf,cons, &
              netcdf,minim,rigid,domain,config,files,comm)
          Else If (minim%freq >  0 .and. flow%step >  0) Then
            If (Mod(flow%step-flow%equil_steps,minim%freq) == 0) Then
              Call minimise_relax(flow%strict .or. cshell%keyshl == SHELL_RELAXED, &
                rdf%l_collect,megatm,thermo%tstep,stat%stpcfg,io,stat,pmf,cons, &
                netcdf,minim,rigid,domain,config,files,comm)
            End If
          End If
        End If

200     Continue

! Refresh mappings

        If (.not.(cshell%relaxed .and. minim%relaxed)) Then
          Call w_refresh_mappings(flow,cshell,cons,pmf,stat,msd_data,bond,angle, &
             dihedral,inversion,tether,neigh,sites,mpoles,rigid,domain,kim_data)
           Go To 100
        End If
     End If

! Get RB COM stress and virial at restart only - also available at w_at_start_vv for levcfg==2

     If (flow%newjob) Then
        If (rigid%total > 0) Then
           If (thermo%l_langevin) Then
              Call langevin_forces(flow%step,thermo%temp,thermo%tstep,thermo%chi, &
                thermo%fxl,thermo%fyl,thermo%fzl,cshell,config,seed)
              If (rigid%share) Then
                Call update_shared_units(config,rigid%list_shared, &
                  rigid%map_shared,thermo%fxl,thermo%fyl,thermo%fzl,domain,comm)
              End If
              Call rigid_bodies_str__s(stat%strcom,config,rigid,comm,thermo%fxl,thermo%fyl,thermo%fzl)
           Else
              Call rigid_bodies_str_ss(stat%strcom,rigid,config,comm)
           End If
           stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
        End If
     End If

! Total virial (excluding constraint, PMF and RB COM virials for npt routines)
! Total stress (excluding constraint, PMF, RB COM and kinetic stress for npt routines)
!
! NOTE(1):  virsrp already includes vdws%vlrc and vlrcm(0) and so
!           does the stress diagonal elements (by minus a third),
!           engsrp includes vdws%elrc and elrcm(0)
!
! NOTE(2):  virfbp, virinv and virdih are allegedly always zero

     stat%virtot = stat%vircpe + stat%virsrp + stat%virter + stat%virtbp + stat%virfbp + &
              stat%virshl + stat%virtet + stat%virbnd + stat%virang + stat%virdih + &
              stat%virinv + stat%virfld

     Call gsum(comm,stat%stress)

! If RBs are present update forces on shared ones

     If (rigid%share) Then
       Call update_shared_units(config,rigid%list_shared, &
         rigid%map_shared,SHARED_UNIT_UPDATE_FORCES,domain,comm)
     End If


!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
