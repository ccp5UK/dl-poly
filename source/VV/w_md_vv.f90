!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate physical quantities at restart
! Calculate kinetic tensor and energy

  Call w_at_start_vv()

! START OF MOLECULAR DYNAMICS CALCULATIONS

  Do While ( (nstep < nstrun .or. (nstep == nstrun .and. newjob)) .and. &
             (timjob-timelp) > timcls )

! Apply impact

     Call w_impact_option()

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart
! levcfg == 2 avoids application twice when forces are calculated at (re)start

     If (newjob) Then
        If (levcfg == 2) Then
           newjob = .false.

           If (keyres /= 1) Call w_write_options()
        End If
     End If

! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

     If (nstep >= 0 .and. nstep < nstrun .and. levcfg == 2) Then

! Increase step counter

        nstep=nstep+1

! zero Kelvin structure optimisation

        If (lzero .and. nstep <= nsteql) Call zero_k_optimise(strkin,strknf,strknt,engke,engrot)

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

        If (Mod(nstep,ndump) == 0 .and. nstep /= nstrun .and. (.not.l_tor)) &
           Call system_revive                                       &
           (imcon,rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

     End If ! DO THAT ONLY IF 0<nstep<=nstrun AND THIS IS AN OLD JOB (newjob=.false.)

! Refresh output

     Call w_refresh_output()

! Complete time check

     Call gtime(timelp)

! Change levcfg after restart if forces and stress are
! (re)calculated and print start-up time

     If (levcfg == 1) Then
        levcfg=2
        If (idnode == 0 .and. lines /= 0) Then
           Write(nrite,'(/,1x,a)') "full force evaluation at restart..."
           Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
        End If
     End If

  End Do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  W_MD_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
