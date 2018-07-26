!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate physical quantities at restart
! Calculate kinetic tensor and energy

  !!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate kinetic tensor and energy at restart

  If (rigid%total > 0) Then
     Call kinstresf(vxx,vyy,vzz,stat%strknf,comm)
     Call kinstrest(rigid,stat%strknt,comm)

     stat%strkin=stat%strknf+stat%strknt
  Else
     Call kinstress(vxx,vyy,vzz,stat%strkin,comm)
  End If
  stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

! If levcfg=2 and RBs are present, update forces on shared ones
! and get RB COM stress and virial at restart.  If levcfg<2
! forces are calculated at (re)start

  If (levcfg == 2) Then
    If (rigid%total > 0) Then
      If (rigid%share) Then
        Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,parts,SHARED_UNIT_UPDATE_FORCES,domain,comm)
      End If

      If (thermo%l_langevin) Then
        Call langevin_forces(nstep,thermo%temp,tstep,thermo%chi,fxl,fyl,fzl,cshell,parts,seed)
        If (rigid%share) Then
          Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,fxl,fyl,fzl,domain,comm)
        End If
        Call rigid_bodies_str__s(stat%strcom,parts(:),rigid,comm,fxl,fyl,fzl)
      Else
        Call rigid_bodies_str_ss(stat%strcom,rigid,parts,comm)
      End If

      stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
    End If
  End If


!!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! START OF MOLECULAR DYNAMICS CALCULATIONS

  Do While ( (nstep < nstrun .or. (nstep == nstrun .and. newjob)) .and. &
             (tmr%job-tmr%elapsed) > tmr%clear_screen )

! Apply impact

     Call w_impact_option(levcfg,nstep,nsteql,rigid,cshell,stat,impa,comm)

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart
! levcfg == 2 avoids application twice when forces are calculated at (re)start

     If (newjob) Then
        If (levcfg == 2) Then
           newjob = .false.

           If (keyres /= 1) Then
             Call w_write_options(cshell,stat,sites,netcdf,domain,traj)
           End If

           If (nstep == 0 .and. nstep == nstrun) Go To 1000
        End If
     End If

! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

     If (nstep >= 0 .and. nstep < nstrun .and. levcfg == 2) Then

! Increase step counter

        nstep=nstep+1

! zero Kelvin structure optimisation

        If (thermo%l_zero .and. nstep <= nsteql .and. Mod(nstep-nsteql,thermo%freq_zero) == 0) Then
          Call zero_k_optimise(stat,rigid,parts,comm)
        End If

! Switch on electron-phonon coupling only after time offset

        l_epcp = (time >= ttmoffset)

! Integrate equations of motion - velocity verlet first stage

        Call w_integrate_vv(0,cshell,cons,pmf,stat,thermo,sites,vdws,rigid,domain,seed,tmr)

! Refresh mappings

        Call w_refresh_mappings(flw,cshell,cons,pmf,stat,msd_data,bond,angle, &
          dihedral,inversion,tether,neigh,sites,mpoles,rigid,domain,kim_data)

      End If ! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

! Evaluate forces

      Call w_calculate_forces(flw,cshell,cons,pmf,stat,plume,pois,bond,angle,dihedral,&
       inversion,tether,threebody,neigh,sites,vdws,tersoffs,fourbody,rdf,netcdf, &
       minim,mpoles,ext_field,rigid,electro,domain,kim_data,tmr)

! Calculate physical quantities, collect statistics and report at t=0

    If (nstep == 0) Then
      Call w_statistics_report(mxatdm_,cshell,cons,pmf,stat,msd_data,zdensity,sites,rdf,domain)
    End If

! DO THAT ONLY IF 0<nstep<=nstrun AND THIS IS AN OLD JOB (newjob=.false.)

     If (nstep > 0 .and. nstep <= nstrun .and. (.not.newjob)) Then

! Evolve electronic temperature for two-temperature model

        If (l_ttm) Then
          Call ttm_ion_temperature(thermo,domain,comm)
          Call ttm_thermal_diffusion(tstep,time,nstep,nsteql,nstbpo,ndump, &
            nstrun,lines,npage,thermo,domain,comm)
        End If

! Integrate equations of motion - velocity verlet second stage

        Call w_integrate_vv(1,cshell,cons,pmf,stat,thermo,sites,vdws,rigid,domain,seed,tmr)

! Apply kinetic options

        Call w_kinetic_options(cshell,cons,pmf,stat,sites,ext_field,domain,seed)

! Update total time of simulation

        time = time + tstep

! Calculate physical quantities, collect statistics and report regularly

        Call w_statistics_report(mxatdm_,cshell,cons,pmf,stat,msd_data,zdensity, &
          sites,rdf,domain)

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        Call w_write_options(cshell,stat,sites,netcdf,domain,traj)

! Save restart data in event of system crash

        If (Mod(nstep,ndump) == 0 .and. nstep /= nstrun .and. (.not.devel%l_tor)) &
           Call system_revive                                 &
           (neigh%cutoff,rbin,megatm,nstep,tstep,time,tmst, &
           stat,devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf,comm)

     End If ! DO THAT ONLY IF 0<nstep<=nstrun AND THIS IS AN OLD JOB (newjob=.false.)

1000 Continue ! Escape forces evaluation at t=0 when nstep=nstrun=0 and newjob=.false.

! Refresh output

     Call w_refresh_output()

! Complete time check

     Call gtime(tmr%elapsed)

! Change levcfg appropriately

     If (levcfg == 1) levcfg=2

  End Do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
