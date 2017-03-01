!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_LFV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Integrate equations of motion - leap-frog verlet

        If (m_rgd == 0) Then
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_0_lfv              &
           (lvar,mndis,mxdis,mxstp,tstep, &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           chit,                          &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,chi,                &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,keyshl,taut,soft,   &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,               &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,          &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 14) Then

! Nose-Hoover thermostat

              Call nvt_g0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,degfre,             &
           sigma,taut,gama,chit,cint,     &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

! Inhomogeneous (two-temperature) 
! Langevin thermostat (Stochastic Dynamics)
		
           Else If (keyens == 15) Then

              Call nvt_l2_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,chi_ep,chi_es,      &
           vel_es2,strkin,engke,          &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,chi,                     &
           press,tai,nstep,chip,eta,      &
           degfre,virtot,                 &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,               &
           press,taup,chip,eta,           &
           virtot,                        &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,          &
           press,taup,chip,eta,           &
           degfre,virtot,                 &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,          &
           press,taup,chip,eta,           &
           degfre,virtot,                 &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_lfv               &
           (lvar,mndis,mxdis,mxstp,tstep,   &
           sigma,chi,                       &
           press,strext,tai,nstep,chip,eta, &
           degfre,iso,ten,stress,           &
           consv,                           &
           strkin,engke,                    &
           mxshak,tolnce,                   &
           megcon,strcon,vircon,            &
           megpmf,strpmf,virpmf,            &
           elrc,virlrc)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,               &
           press,strext,taup,chip,eta,    &
           iso,ten,stress,                &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,          &
           press,strext,taup,chip,eta,    &
           degfre,iso,ten,stress,         &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_lfv             &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,          &
           press,strext,taup,chip,eta,    &
           degfre,iso,ten,stress,         &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

           Else

! Invalid ensemble option

           Call error(430)

           End If
        Else
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_1_lfv                  &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           chit,                              &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           nstep,temp,chi,                    &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           nstep,temp,keyshl,taut,soft,       &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,                   &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,cint,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens == 14) Then

! Nose-Hoover thermostat

              Call nvt_g1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           nstep,temp,degfre,                 &
           sigma,taut,gama,chit,cint,         &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,chi,                         &
           press,tai,nstep,chip,eta,          &
           degfre,degrot,virtot,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,                   &
           press,taup,chip,eta,               &
           virtot,                            &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,cint,              &
           press,taup,chip,eta,               &
           degfre,degrot,virtot,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,cint,              &
           press,taup,chip,eta,               &
           degfre,degrot,virtot,              &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,chi,                         &
           press,strext,tai,nstep,chip,eta,   &
           degfre,degrot,iso,ten,stress,      &
           strkin,strknf,strknt,engke,engrot, &
           consv,                             &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,                   &
           press,strext,taup,chip,eta,        &
           iso,ten,stress,                    &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,cint,              &
           press,strext,taup,chip,eta,        &
           degfre,degrot,iso,ten,stress,      &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_lfv                 &
           (lvar,mndis,mxdis,mxstp,tstep,     &
           sigma,taut,chit,cint,              &
           press,strext,taup,chip,eta,        &
           degfre,degrot,iso,ten,stress,      &
           consv,                             &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           mxquat,quattol,strcom,vircom,      &
           elrc,virlrc)

           Else

! Invalid ensemble option

           Call error(430)

           End If
        End If


!!!!!!!!!!!!!!!!!!!!!  W_INTEGRATE_LFV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
