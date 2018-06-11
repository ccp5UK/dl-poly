!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Sharlow's splittings for VV only (LFV->VV) DPD thermostat - no variable time-stepping!!!
! One-off application for first order splitting and symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stat%stress due to random and drag forces

        If (thermo%key_dpd > 0 .and. thermo%key_dpd*isw == 0) Then
          Call dpd_thermostat(isw,l_str,neigh%cutoff,nstep,tstep,stat,thermo,comm)
        End If

! Integrate equations of motion - velocity verlet
! first (isw == 0) or second (isw == 1) stage

        If (m_rgd == 0) Then
           If      (thermo%ensemble == ENS_NVE) Then

! Microcanonical ensemble

              Call nve_0_vv                   &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_EVANS) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           thermo%chi_t,                              &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_LANGEVIN) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,                    &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_ANDERSON) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,keyshl,       &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_BERENDSEN) Then

! Berendsen thermostat

              Call nvt_b0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_NOSE_HOOVER) Then

! Nose-Hoover thermostat

              Call nvt_h0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_GENTLE) Then

! Gentle-Stochastic thermostat

              Call nvt_g0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,degfre,                 &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_LANGEVIN_INHOMO) Then

! Inhomogeneous (two-temperature)
! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l2_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,  &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN) Then

! Langevin thermostat and isotropic barostat

              Call npt_l0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,          &
           degfre,stat%virtot,                     &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%virtot,                            &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,stat%virtot,                     &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,stat%virtot,                     &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,   &
           degfre,stat%stress,             &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%stress,                    &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,stat%stress,             &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,stat%stress,             &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           megpmf,stat%strpmf,stat%virpmf,              &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        Else
           If      (thermo%ensemble == ENS_NVE) Then

! Microcanonical ensemble

              Call nve_1_vv                   &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_EVANS) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           thermo%chi_t,                              &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_LANGEVIN) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,                    &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_ANDERSON) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,keyshl,       &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_BERENDSEN) Then

! Berendsen thermostat

              Call nvt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_NOSE_HOOVER) Then

! Nose-Hoover thermostat

              Call nvt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NVT_GENTLE) Then

! Gentle-Stochastic thermostat

              Call nvt_g1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,degfre,                 &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,          &
           degfre,degrot,stat%virtot,              &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%virtot,                            &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,degrot,stat%virtot,              &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,degrot,stat%virtot,              &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,   &
           degfre,degrot,stat%stress,      &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%consv,                             &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           stat%stress,                    &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,degrot,stat%stress,      &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,degrot,stat%stress,      &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           megpmf,stat%strpmf,stat%virpmf,              &
           stat%strcom,stat%vircom,                     &
           elrc,virlrc,cons,stat,thermo,tmr,comm)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        End If

! Sharlow's second order splittings for VV only (LFV->VV) DPD thermostat - no variable time-stepping!!!
! Symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stat%stress due to random and drag forces

        If (thermo%key_dpd > 0 .and. thermo%key_dpd*isw == 2) Then
          Call dpd_thermostat(isw,l_str,neigh%cutoff,nstep,tstep,stat,thermo,comm)
        End If


!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
