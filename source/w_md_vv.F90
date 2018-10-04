!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate physical quantities at restart
! Calculate kinetic tensor and energy

  !!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate kinetic tensor and energy at restart

  If (rigid%total > 0) Then
     Call kinstresf(config%vxx,config%vyy,config%vzz,stat%strknf,config,comm)
     Call kinstrest(rigid,stat%strknt,comm)

     stat%strkin=stat%strknf+stat%strknt
  Else
     Call kinstress(config%vxx,config%vyy,config%vzz,stat%strkin,config,comm)
  End If
  stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

! If levcfg=2 and RBs are present, update forces on shared ones
! and get RB COM stress and virial at restart.  If levcfg<2
! forces are calculated at (re)start

  If (levcfg == 2) Then
    If (rigid%total > 0) Then
      If (rigid%share) Then
        Call update_shared_units(config,rigid%list_shared,rigid%map_shared,SHARED_UNIT_UPDATE_FORCES,domain,comm)
      End If

      If (thermo%l_langevin) Then
        Call langevin_forces(flow%step,thermo%temp,thermo%tstep,thermo%chi,thermo%fxl,thermo%fyl,thermo%fzl,cshell,config,seed)
        If (rigid%share) Then
          Call update_shared_units(config,rigid%list_shared,rigid%map_shared,&
              thermo%fxl,thermo%fyl,thermo%fzl,domain,comm)
        End If
        Call rigid_bodies_str__s(stat%strcom,config,rigid,comm,thermo%fxl,thermo%fyl,thermo%fzl)
      Else
        Call rigid_bodies_str_ss(stat%strcom,rigid,config,comm)
      End If

      stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
    End If
  End If


!!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! START OF MOLECULAR DYNAMICS CALCULATIONS

  Do While ( (flow%step < flow%run_steps .or. (flow%step == flow%run_steps .and. flow%newjob)) .and. &
             (tmr%job-tmr%elapsed) > tmr%clear_screen )

! Apply impact

     Call w_impact_option(levcfg,flow%step,flow%equil_steps,rigid,cshell,stat,impa,config,comm)

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart
! levcfg == 2 avoids application twice when forces are calculated at (re)start

     If (flow%newjob) Then
       If (levcfg == 2) Then
          flow%newjob = .false.

           If (flow%restart_key /= 1) Then
             Call w_write_options(io,rsdc,cshell,stat,sites,netcdf,domain,traj,files)
           End If

           If (flow%step == 0 .and. flow%step == flow%run_steps) Go To 1000
        End If
     End If

! DO THAT ONLY IF 0<=flow%step<flow%run_steps AND FORCES ARE PRESENT (levcfg=2)

     If (flow%step >= 0 .and. flow%step < flow%run_steps .and. levcfg == 2) Then

! Increase step counter

        flow%step=flow%step+1

! zero Kelvin structure optimisation

        If (thermo%l_zero .and. flow%step <= flow%equil_steps .and. Mod(flow%step-flow%equil_steps,thermo%freq_zero) == 0) Then
          Call zero_k_optimise(stat,rigid,config,comm)
        End If

! Switch on electron-phonon coupling only after flow%time offset

        ttm%l_epcp = (flow%time >= ttm%ttmoffset)

! Integrate equations of motion - velocity verlet first stage

        Call w_integrate_vv(0,ttm,cshell,cons,pmf,stat,thermo,sites,vdws,rigid,domain,seed,tmr)

! Refresh mappings

        Call w_refresh_mappings(flow,cshell,cons,pmf,stat,msd_data,bond,angle, &
          dihedral,inversion,tether,neigh,sites,mpoles,rigid,domain,kim_data)

      End If ! DO THAT ONLY IF 0<=flow%step<flow%run_steps AND FORCES ARE PRESENT (levcfg=2)

! Evaluate forces

      Call w_calculate_forces(cnfig,flow,io,cshell,cons,pmf,stat,plume,pois,bond,angle,dihedral,&
       inversion,tether,threebody,neigh,sites,vdws,tersoffs,fourbody,rdf,netcdf, &
       minim,mpoles,ext_field,rigid,electro,domain,kim_data,tmr)

! Calculate physical quantities, collect statistics and report at t=0

    If (flow%step == 0) Then
      Call w_statistics_report(cnfig%mxatdm,cshell,cons,pmf,stat,msd_data,zdensity, &
        sites,rdf,domain,flow,files)
    End If

! DO THAT ONLY IF 0<flow%step<=flow%run_steps AND THIS IS AN OLD JOB (flow%newjob=.false.)

     If (flow%step > 0 .and. flow%step <= flow%run_steps .and. (.not.flow%newjob)) Then

! Evolve electronic temperature for two-temperature model

       If (ttm%l_ttm) Then
         Call ttm_ion_temperature(ttm,thermo,domain,config,comm)
          Call ttm_thermal_diffusion(thermo%tstep,flow%time,flow%step,flow%equil_steps,flow%freq_output,flow%freq_restart, &
            flow%run_steps,ttm,thermo,domain,comm)
        End If

! Integrate equations of motion - velocity verlet second stage

        Call w_integrate_vv(1,ttm,cshell,cons,pmf,stat,thermo,sites,vdws,rigid,domain,seed,tmr)

! Apply kinetic options

        Call w_kinetic_options(cshell,cons,pmf,stat,sites,ext_field,domain,seed)

! Update total flow%time of simulation

        flow%time = flow%time + thermo%tstep

! Calculate physical quantities, collect statistics and report regularly

        Call w_statistics_report(cnfig%mxatdm,cshell,cons,pmf,stat,msd_data,zdensity, &
          sites,rdf,domain,flow,files)

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        Call w_write_options(io,rsdc,cshell,stat,sites,netcdf,domain,traj,files)

! Save restart data in event of system crash

        If (Mod(flow%step,flow%freq_restart) == 0 .and. flow%step /= flow%run_steps .and. (.not.devel%l_tor)) Then
          Call system_revive(neigh%cutoff,megatm,flow%step,flow%time,sites,io,flow%start_time, &
            stat,devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf, &
            netcdf,config,files,comm)
        End If

     End If ! DO THAT ONLY IF 0<flow%step<=flow%run_steps AND THIS IS AN OLD JOB (flow%newjob=.false.)

1000 Continue ! Escape forces evaluation at t=0 when flow%step=flow%run_steps=0 and flow%newjob=.false.

! Refresh output

     Call w_refresh_output()

! Complete flow%time check

     Call gtime(tmr%elapsed)

! Change levcfg appropriately

     If (levcfg == 1) levcfg=2

  End Do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
