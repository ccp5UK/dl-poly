!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Apply external field

        If (ext_field%key /= FIELD_NULL) Then
          Call external_field_correct(stat%engfld,ext_field,rigid,config,comm)
        End If

! Apply pseudo thermostat - velocity cycle (1)
        If (thermo%l_stochastic_boundaries) Then
          Call stochastic_boundary_vv(1,thermo%tstep,nstep,sites%dof_site,cshell,stat,thermo,rigid,domain,config,seed,comm)
        End If

! Apply temperature regaussing

        If (thermo%l_tgaus .and. nstep <= nsteql .and. Mod(nstep-nsteql,thermo%freq_tgaus) == 0) Then
           thermo%chi_t = 0.0_wp
           thermo%chi_p = 0.0_wp
           thermo%eta  = 0.0_wp

           Call regauss_temperature(rigid,domain,config,seed,comm)

! quench constraints & PMFs

           If (cons%megcon > 0) Call constraints_quench(cons,stat,domain,config,comm)
           If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,config,comm)

! quench core-shell units in adiabatic model

           If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
              stat%stptmp = 2.0_wp*(stat%engke+stat%engrot) / (boltz*Real(degfre,wp))
              Do
                 Call scale_temperature(stat%engke+stat%engrot,degtra,degrot,degfre,rigid,config,comm)
                 Call core_shell_quench(config,safe,stat%stptmp,cshell,domain,comm)
                 If (cons%megcon > 0) Call constraints_quench(cons,stat,domain,config,comm)
                 If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,config,comm)
                 If (rigid%total > 0) Call rigid_bodies_quench(rigid,domain,config,comm)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(stat%engke+stat%engrot,degtra,degrot,degfre,rigid,config,comm)
           End If

! Correct kinetic stress and energy

           If (rigid%total > 0) Then
              Call kinstresf(config%vxx,config%vyy,config%vzz,stat%strknf,config,comm)
              Call kinstrest(rigid,stat%strknt,comm)

              stat%strkin=stat%strknf+stat%strknt

              stat%engrot=getknr(rigid,comm)
           Else
              Call kinstress(config%vxx,config%vyy,config%vzz,stat%strkin,config,comm)
           End If
           stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))
        End If

! Apply temperature scaling

        If (thermo%l_tscale .and. nstep <= nsteql .and. Mod(nstep-nsteql,thermo%freq_tscale) == 0) Then
           thermo%chi_t = 0.0_wp
           thermo%chi_p = 0.0_wp
           thermo%eta  = 0.0_wp

! quench constraints & PMFs

           If (cons%megcon > 0) Call constraints_quench(cons,stat,domain,config,comm)
           If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,config,comm)

! quench core-shell units in adiabatic model

           If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
              Do
                 Call scale_temperature(thermo%sigma,degtra,degrot,degfre,rigid,config,comm)
                 Call core_shell_quench(config,safe,stat%stptmp,cshell,domain,comm)
                 If (cons%megcon > 0) Call constraints_quench(cons,stat,domain,config,comm)
                 If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,config,comm)
                 If (rigid%total > 0) Call rigid_bodies_quench(rigid,domain,config,comm)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(thermo%sigma,degtra,degrot,degfre,rigid,config,comm)
           End If

! Correct kinetic stress and energy

           If (rigid%total > 0) Then
              Call kinstresf(config%vxx,config%vyy,config%vzz,stat%strknf,config,comm)
              Call kinstrest(rigid,stat%strknt,comm)

              stat%strkin=stat%strknf+stat%strknt

              stat%engrot=getknr(rigid,comm)
           Else
              Call kinstress(config%vxx,config%vyy,config%vzz,stat%strkin,config,comm)
           End If
           stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))
        End If


!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
