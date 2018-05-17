!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Sharlow's splittings for VV only (LFV->VV) DPD thermostat - no variable time-stepping!!!
! One-off application for first order splitting and symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stress due to random and drag forces

        If (thermo%key_dpd > 0 .and. thermo%key_dpd*isw == 0) Call dpd_thermostat(isw,l_str,rcut,nstep,tstep,comm)

! Integrate equations of motion - velocity verlet
! first (isw == 0) or second (isw == 1) stage

        If (m_rgd == 0) Then
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_0_vv                   &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           chit,                              &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,thermo%chi,                    &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,keyshl,thermo%tau_t,thermo%soft,       &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,                   &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,degfre,                 &
           sigma,thermo%tau_t,thermo%gama,chit,cint,         &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 15) Then

! Inhomogeneous (two-temperature)
! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l2_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,thermo%chi_ep,thermo%chi_es,vel_es2,  &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,comm)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%chi,                         &
           thermo%press,thermo%tai,nstep,chip,eta,          &
           degfre,virtot,                     &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,                   &
           thermo%press,thermo%tau_p,chip,eta,               &
           virtot,                            &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%tau_p,chip,eta,               &
           degfre,virtot,                     &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%tau_p,chip,eta,               &
           degfre,virtot,                     &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%chi,                         &
           thermo%press,thermo%stress,thermo%tai,nstep,chip,eta,   &
           degfre,thermo%iso,thermo%tension,stress,             &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,                   &
           thermo%press,thermo%stress,thermo%tau_p,chip,eta,        &
           thermo%iso,thermo%tension,stress,                    &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%stress,thermo%tau_p,chip,eta,        &
           degfre,thermo%iso,thermo%tension,stress,             &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%stress,thermo%tau_p,chip,eta,        &
           degfre,thermo%iso,thermo%tension,stress,             &
           consv,                             &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc,comm)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        Else
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_1_vv                   &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           chit,                              &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,thermo%chi,                    &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,keyshl,thermo%tau_t,thermo%soft,       &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,                   &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,thermo%temp,degfre,                 &
           sigma,thermo%tau_t,thermo%gama,chit,cint,         &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,comm)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%chi,                         &
           thermo%press,thermo%tai,nstep,chip,eta,          &
           degfre,degrot,virtot,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,                   &
           thermo%press,thermo%tau_p,chip,eta,               &
           virtot,                            &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%tau_p,chip,eta,               &
           degfre,degrot,virtot,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%tau_p,chip,eta,               &
           degfre,degrot,virtot,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%chi,                         &
           thermo%press,thermo%stress,thermo%tai,nstep,chip,eta,   &
           degfre,degrot,thermo%iso,thermo%tension,stress,      &
           strkin,strknf,strknt,engke,engrot, &
           consv,                             &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,                   &
           thermo%press,thermo%stress,thermo%tau_p,chip,eta,        &
           thermo%iso,thermo%tension,stress,                    &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%stress,thermo%tau_p,chip,eta,        &
           degfre,degrot,thermo%iso,thermo%tension,stress,      &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,thermo%tau_t,chit,cint,              &
           thermo%press,thermo%stress,thermo%tau_p,chip,eta,        &
           degfre,degrot,thermo%iso,thermo%tension,stress,      &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom,                     &
           elrc,virlrc,comm)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        End If

! Sharlow's second order splittings for VV only (LFV->VV) DPD thermostat - no variable time-stepping!!!
! Symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stress due to random and drag forces

        If (thermo%key_dpd > 0 .and. thermo%key_dpd*isw == 2) Call dpd_thermostat(isw,l_str,rcut,nstep,tstep,comm)


!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
