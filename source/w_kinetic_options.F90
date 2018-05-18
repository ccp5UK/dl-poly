!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Apply external field

        If (keyfld > 0) Call external_field_correct(engfld,comm)

! Apply pseudo thermostat - velocity cycle (1)

        If (thermo%l_pseudo) Then
              Call pseudo_vv                            &
           (1,keyshl,keyens,tstep, &
           nstep,strkin,strknf,strknt,engke,engrot,thermo,comm)
        End If

! Apply temperature regaussing

        If (thermo%l_tgaus .and. nstep <= nsteql .and. Mod(nstep-nsteql,thermo%freq_tgaus) == 0) Then
           chit = 0.0_wp
           chip = 0.0_wp
           eta  = 0.0_wp

           Call regauss_temperature(megrgd,comm)

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(mxshak,tolnce,comm)
           If (megpmf > 0) Call pmf_quench(mxshak,tolnce,comm)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              stptmp = 2.0_wp*(engke+engrot) / (boltz*Real(degfre,wp))
              Do
                 Call scale_temperature(engke+engrot,degtra,degrot,degfre,comm)
                 Call core_shell_quench(safe,stptmp,comm)
                 If (megcon > 0) Call constraints_quench(mxshak,tolnce,comm)
                 If (megpmf > 0) Call pmf_quench(mxshak,tolnce,comm)
                 If (megrgd > 0) Call rigid_bodies_quench(comm)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(engke+engrot,degtra,degrot,degfre,comm)
           End If

! Correct kinetic stress and energy

           If (megrgd > 0) Then
              Call kinstresf(vxx,vyy,vzz,strknf,comm)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)
           Else
              Call kinstress(vxx,vyy,vzz,strkin,comm)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If

! Apply temperature scaling

        If (thermo%l_tscale .and. nstep <= nsteql .and. Mod(nstep-nsteql,thermo%freq_tscale) == 0) Then
           chit = 0.0_wp
           chip = 0.0_wp
           eta  = 0.0_wp

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(mxshak,tolnce,comm)
           If (megpmf > 0) Call pmf_quench(mxshak,tolnce,comm)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              Do
                 Call scale_temperature(sigma,degtra,degrot,degfre,comm)
                 Call core_shell_quench(safe,stptmp,comm)
                 If (megcon > 0) Call constraints_quench(mxshak,tolnce,comm)
                 If (megpmf > 0) Call pmf_quench(mxshak,tolnce,comm)
                 If (megrgd > 0) Call rigid_bodies_quench(comm)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(sigma,degtra,degrot,degfre,comm)
           End If

! Correct kinetic stress and energy

           If (megrgd > 0) Then
              Call kinstresf(vxx,vyy,vzz,strknf,comm)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)
           Else
              Call kinstress(vxx,vyy,vzz,strkin,comm)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If


!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
