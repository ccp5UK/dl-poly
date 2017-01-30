!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_LFV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Fix levcfg=2 as forces are always calculated at restart in LFV
! Calculate kinetic stress tensor and energy at t=0 for reporting

  Call w_at_start_lfv()

! START OF MOLECULAR DYNAMICS CALCULATIONS

  Do While ( (nstep < nstrun .or. (nstep == nstrun .and. newjob)) .and. &
             (timjob-timelp) > timcls )

! Apply impact

     Call w_impact_option()

! Evaluate forces

     Call w_calculate_forces()

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart

     If (newjob) Then
        newjob = .false.

        If (keyres /= 1) Call w_write_options()
     End If

! Calculate physical quantities, collect statistics and report at t=0

     If (nstep == 0) Call w_statistics_report()

! DO THAT ONLY IF 0<=nstep<nstrun

     If (nstep >= 0 .and. nstep < nstrun) Then

! Increase step counter

        nstep=nstep+1

! zero Kelvin structure optimisation

        If (lzero .and. nstep <= nsteql .and. Mod(nstep-nsteql,nstzero) == 0) &
           Call zero_k_optimise(strkin,strknf,strknt,engke,engrot)

! Integrate equations of motion - leap-frog verlet

        Call w_integrate_lfv()

! Refresh mappings

        Call w_refresh_mappings()

! Apply kinetic options

        Call w_kinetic_options()

! Update total time of simulation

        time = time + tstep

! Calculate physical quantities, collect statistics and report regularly

        Call w_statistics_report()

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        Call w_write_options()

! Save restart data in event of system crash

        If (Mod(nstep,ndump) == 0 .and. nstep /= nstrun .and. (.not.l_tor)) &
           Call system_revive                                 &
           (rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

     End If ! DO THAT ONLY IF 0<=nstep<nstrun

! Refresh output

     Call w_refresh_output()

! Complete time check

     Call gtime(timelp)

  End Do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_LFV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
