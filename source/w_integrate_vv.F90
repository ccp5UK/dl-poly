!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Sharlow's splittings for VV only (LFV->VV) DPD thermostat - no variable flow%time-stepping!!!
! One-off application for first order splitting and symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stat%stress due to random and drag forces

      If (thermo%key_dpd > 0 .and. thermo%key_dpd*stage == 0) Then
        Call dpd_thermostat(stage,flow%strict,neigh%cutoff,flow%step,thermo%tstep,stat,thermo, &
          neigh,rigid,domain,cnfig,seed,comm)
      End If

! Integrate equations of motion - velocity verlet
! first (stage == 0) or second (stage == 1) stage

      If (.not. rigid%on) Then
         If (thermo%ensemble == ENS_NVE) Then

! Microcanonical ensemble

            Call nve_0_vv                   &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         stat%strkin,stat%engke,thermo,               &
         cshell,cons,pmf,stat,domain,tmr,cnfig,comm)

         Else If (thermo%ensemble == ENS_NVT_EVANS) Then

! Evans thermostat (Gaussian temperature constraints)

            Call nvt_e0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         thermo%chi_t,                              &
         stat%strkin,stat%engke,thermo,               &
         cshell,cons,pmf,stat,domain,tmr,cnfig,comm)

         Else If (thermo%ensemble == ENS_NVT_LANGEVIN) Then

! Langevin thermostat (Stochastic Dynamics)

            Call nvt_l0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         flow%step,                    &
         stat%strkin,stat%engke,                      &
         cshell,cons,pmf,stat,thermo,domain,tmr,cnfig,seed,comm)

         Else If (thermo%ensemble == ENS_NVT_ANDERSON) Then

! Andersen thermostat (Stochastic Dynamics)

            Call nvt_a0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         flow%step,       &
         stat%strkin,stat%engke,                      &
         cshell,cons,pmf,stat,thermo,sites,domain,tmr,cnfig,seed,comm)

         Else If (thermo%ensemble == ENS_NVT_BERENDSEN) Then

! Berendsen thermostat

            Call nvt_b0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         stat%strkin,stat%engke,                      &
         cshell,cons,pmf,stat,thermo,domain,tmr,cnfig,comm)

         Else If (thermo%ensemble == ENS_NVT_NOSE_HOOVER) Then

! Nose-Hoover thermostat

            Call nvt_h0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         stat%consv,                             &
         stat%strkin,stat%engke,                      &
         cshell,cons,pmf,stat,thermo,domain,tmr,cnfig,comm)

         Else If (thermo%ensemble == ENS_NVT_GENTLE) Then

! Gentle-Stochastic thermostat

            Call nvt_g0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         flow%step,cnfig%degfre,                 &
         stat%consv,                             &
         stat%strkin,stat%engke,                      &
         cshell,cons,pmf,stat,thermo,domain,tmr,cnfig,seed,comm)

         Else If (thermo%ensemble == ENS_NVT_LANGEVIN_INHOMO) Then

! Inhomogeneous (two-temperature)
! Langevin thermostat (Stochastic Dynamics)

            Call nvt_l2_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         flow%step,  &
              stat%strkin,stat%engke,                      &
              ttm,cshell,cons,pmf,stat,thermo,domain,tmr,cnfig,seed,comm)

         Else If (thermo%ensemble == ENS_NPT_LANGEVIN) Then

! Langevin thermostat and isotropic barostat

            Call npt_l0_vv                  &
         (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
         flow%step,          &
         cnfig%degfre,stat%virtot,                     &
         stat%consv,                             &
         stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,domain,&
           tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%virtot,                            &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,domain, &
           tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,stat%virtot,                     &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,rigid,&
           domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,stat%virtot,                     &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,domain,&
           tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           flow%step,   &
           cnfig%degfre,stat%stress,             &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,domain,&
           tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%stress,                    &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,domain,&
           tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,stat%stress,             &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,rigid,&
           domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,stat%stress,             &
           stat%consv,                             &
           stat%strkin,stat%engke,                      &
           cshell,cons,pmf,stat,thermo,sites,vdws,domain,&
           tmr,cnfig,comm)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        Else
           If      (thermo%ensemble == ENS_NVE) Then

! Microcanonical ensemble

              Call nve_1_vv                   &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,thermo,cshell,cons,pmf,&
           stat,rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NVT_EVANS) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           thermo%chi_t,                              &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,thermo,cshell,cons,pmf,&
           stat,rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NVT_LANGEVIN) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           flow%step,                    &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,cshell,cons,pmf,&
           stat,thermo,rigid,domain,tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NVT_ANDERSON) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           flow%step,       &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,cshell,cons,pmf,&
           stat,thermo,sites,rigid,domain,tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NVT_BERENDSEN) Then

! Berendsen thermostat

              Call nvt_b1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,cshell,cons,pmf,&
           stat,thermo,rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NVT_NOSE_HOOVER) Then

! Nose-Hoover thermostat

              Call nvt_h1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,cshell,cons,pmf,&
           stat,thermo,rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NVT_GENTLE) Then

! Gentle-Stochastic thermostat

              Call nvt_g1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           flow%step,cnfig%degfre,                 &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,cshell,cons,pmf,&
           stat,thermo,rigid,domain,tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           flow%step,          &
           cnfig%degfre,cnfig%degrot,stat%virtot,              &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%virtot,                            &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,cnfig%degrot,stat%virtot,              &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,cnfig%degrot,stat%virtot,              &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           flow%step,   &
           cnfig%degfre,cnfig%degrot,stat%stress,      &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%consv,                             &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,seed,comm)

           Else If (thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           stat%stress,                    &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,cnfig%degrot,stat%stress,      &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,comm)

           Else If (thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_vv                  &
           (stage,thermo%lvar,thermo%mndis,thermo%mxdis,thermo%mxstp,thermo%tstep, &
           cnfig%degfre,cnfig%degrot,stat%stress,      &
           stat%consv,                             &
           stat%strkin,stat%strknf,stat%strknt,stat%engke,stat%engrot, &
           stat%strcom,stat%vircom,                     &
           cshell,cons,pmf,stat,thermo,sites,vdws,&
           rigid,domain,tmr,cnfig,comm)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        End If

! Sharlow's second order splittings for VV only (LFV->VV) DPD thermostat - no variable flow%time-stepping!!!
! Symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stat%stress due to random and drag forces

        If (thermo%key_dpd > 0 .and. thermo%key_dpd*stage == 2) Then
          Call dpd_thermostat(stage,flow%strict,neigh%cutoff,flow%step,thermo%tstep,stat,thermo, &
            neigh,rigid,domain,cnfig,seed,comm)
        End If


!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
