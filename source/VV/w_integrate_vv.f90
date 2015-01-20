!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Sharlow's splittings for VV only (LFV->VV) DPD thermostat - no variable time-stepping!!!
! One-off application for first order splitting and symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stress due to random and drag forces

        If (keydpd > 0 .and. keydpd*isw == 0) Call dpd_thermostat(isw,l_str,imcon,rcut,nstep,tstep)

! Integrate equations of motion - velocity verlet
! first (isw == 0) or second (isw == 1) stage

        If (m_rgd == 0) Then
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_0_vv                                &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           chit)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l0_vv                                        &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,chi,strkin,engke, &
           nstep,imcon,mxshak,tolnce,megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,         &
           keyshl,taut,soft,strkin,engke,                  &
           nstep,imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit,cint,consv)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g0_vv                                    &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,strkin,engke, &
           nstep,imcon,mxshak,tolnce,megcon,strcon,vircon,      &
           megpmf,strpmf,virpmf,                                &
           degfre,sigma,taut,gama,chit,cint,consv)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           nstep,imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,chi,consv,                         &
           press,tai,chip,eta,virtot,                      &
           elrc,virlrc)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit,                                &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,taut,chit,cint,consv,              &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,taut,chit,cint,consv,              &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           nstep,imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,chi,consv,                     &
           press,strext,ten,tai,chip,eta,stress,           &
           elrc,virlrc)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,sigma,taut,chit,                            &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,taut,chit,cint,consv,          &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,taut,chit,cint,consv,          &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

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
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           chit,                              &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           chi,                                    &
           strkin,strknf,strknt,engke,engrot,      &
           nstep,imcon,mxshak,tolnce,              &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           keyshl,taut,soft,                       &
           strkin,strknf,strknt,engke,engrot,      &
           nstep,imcon,mxshak,tolnce,              &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,consv,        &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           degfre,sigma,taut,gama,chit,cint,consv, &
           strkin,strknf,strknt,engke,engrot,      &
           nstep,imcon,mxshak,tolnce,              &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,chi,consv,            &
           degrot,press,tai,chip,eta,         &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           nstep,imcon,mxshak,tolnce,         &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           press,taup,chip,eta,               &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,        &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,        &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           iso,degfre,sigma,chi,consv,        &
           degrot,press,tai,chip,eta,         &
           stress,strext,ten,elrc,virlrc,     &
           strkin,strknf,strknt,engke,engrot, &
           nstep,imcon,mxshak,tolnce,         &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           iso,sigma,taut,chit,               &
           press,taup,chip,eta,               &
           stress,strext,ten,elrc,virlrc,     &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_vv                      &
           (isw,lvar,mndis,mxdis,mxstp,tstep,     &
           iso,degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,            &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,                   &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_vv                      &
           (isw,lvar,mndis,mxdis,mxstp,tstep,     &
           iso,degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,            &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,                   &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        End If

! Sharlow's second order splittings for VV only (LFV->VV) DPD thermostat - no variable time-stepping!!!
! Symmetric application for second order splitting
! Velocity field change + generation of DPD virial & stress due to random and drag forces

        If (keydpd > 0 .and. keydpd*isw == 2) Call dpd_thermostat(isw,l_str,imcon,rcut,nstep,tstep)


!!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
