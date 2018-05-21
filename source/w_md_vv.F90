!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate physical quantities at restart
! Calculate kinetic tensor and energy

  !!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate kinetic tensor and energy at restart

  If (megrgd > 0) Then
     Call kinstresf(vxx,vyy,vzz,strknf,comm)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

     strkin=strknf+strknt
  Else
     Call kinstress(vxx,vyy,vzz,strkin,comm)
  End If
  engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! If levcfg=2 and RBs are present, update forces on shared ones
! and get RB COM stress and virial at restart.  If levcfg<2
! forces are calculated at (re)start

  If (levcfg == 2) Then
     If (megrgd > 0) Then
        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz,comm)

        If (thermo%l_langevin) Then
           Call langevin_forces(nstep,thermo%temp,tstep,thermo%chi,fxl,fyl,fzl)
           If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxl,fyl,fzl,comm)
           Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl,comm)
        Else
           Call rigid_bodies_str_ss(strcom,comm)
        End If

        vircom=-(strcom(1)+strcom(5)+strcom(9))
     End If
  End If


!!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! START OF MOLECULAR DYNAMICS CALCULATIONS

  Do While ( (nstep < nstrun .or. (nstep == nstrun .and. newjob)) .and. &
             (tmr%job-tmr%elapsed) > tmr%clear_screen )

! Apply impact

     Call w_impact_option(levcfg,nstep,nsteql,engke,engrot,megrgd,strkin,strknf,strknt,impa,comm)

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart
! levcfg == 2 avoids application twice when forces are calculated at (re)start

     If (newjob) Then
        If (levcfg == 2) Then
           newjob = .false.

           If (keyres /= 1) Call w_write_options()

           If (nstep == 0 .and. nstep == nstrun) Go To 1000
        End If
     End If

! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

     If (nstep >= 0 .and. nstep < nstrun .and. levcfg == 2) Then

! Increase step counter

        nstep=nstep+1

! zero Kelvin structure optimisation

        If (thermo%l_zero .and. nstep <= nsteql .and. Mod(nstep-nsteql,thermo%freq_zero) == 0) &
           Call zero_k_optimise(strkin,strknf,strknt,engke,engrot,comm)

! Switch on electron-phonon coupling only after time offset

        l_epcp = (time >= ttmoffset)

! Integrate equations of motion - velocity verlet first stage

        Call w_integrate_vv(0)

! Refresh mappings

        Call w_refresh_mappings()

     End If ! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

! Evaluate forces

     Call w_calculate_forces()

! Calculate physical quantities, collect statistics and report at t=0

     If (nstep == 0) Call w_statistics_report()

! DO THAT ONLY IF 0<nstep<=nstrun AND THIS IS AN OLD JOB (newjob=.false.)

     If (nstep > 0 .and. nstep <= nstrun .and. (.not.newjob)) Then

! Evolve electronic temperature for two-temperature model

        If (l_ttm) Then
          Call ttm_ion_temperature(vel_es2,thermo,comm)
          Call ttm_thermal_diffusion(tstep,time,nstep,nsteql,nstbpo,ndump,nstrun,lines,npage,thermo,comm)
        End If

! Integrate equations of motion - velocity verlet second stage

        Call w_integrate_vv(1)

! Apply kinetic options

        Call w_kinetic_options()

! Update total time of simulation

        time = time + tstep

! Calculate physical quantities, collect statistics and report regularly

        Call w_statistics_report()

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        Call w_write_options()

! Save restart data in event of system crash

        If (Mod(nstep,ndump) == 0 .and. nstep /= nstrun .and. (.not.devel%l_tor)) &
           Call system_revive                                 &
           (rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress,devel,comm)

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
