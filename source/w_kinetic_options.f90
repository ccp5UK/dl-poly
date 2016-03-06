!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Apply external field

        If (keyfld > 0) Call external_field_correct()

! Apply pseudo thermostat - velocity cycle (1)

        If (lpse) Then
           If (l_vv) Then
              Call pseudo_vv                            &
           (1,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           nstep,strkin,strknf,strknt,engke,engrot)
           Else
              Call pseudo_lfv                           &
           (1,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           nstep,strkin,strknf,strknt,engke,engrot)
           End If
        End If

! Apply temperature regaussing

        If (ltgaus .and. nstep <= nsteql .and. Mod(nstep-nsteql,nstgaus) == 0) Then
           chit = 0.0_wp
           chip = 0.0_wp
           eta  = 0.0_wp

           Call regauss_temperature(megrgd)

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(mxshak,tolnce)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              stptmp = 2.0_wp*(engke+engrot) / (boltz*Real(degfre,wp))
              Do
                 Call scale_temperature(engke+engrot,degtra,degrot,degfre)
                 Call core_shell_quench(safe,stptmp)
                 If (megcon > 0) Call constraints_quench(mxshak,tolnce)
                 If (megpmf > 0) Call pmf_quench(mxshak,tolnce)
                 If (megrgd > 0) Call rigid_bodies_quench()
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(engke+engrot,degtra,degrot,degfre)
           End If

! Correct kinetic stress and energy

           If (megrgd > 0) Then
              Call kinstresf(vxx,vyy,vzz,strknf)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz)
           Else
              Call kinstress(vxx,vyy,vzz,strkin)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If

! Apply temperature scaling

        If (ltscal .and. nstep <= nsteql .and. Mod(nstep-nsteql,nstscal) == 0) Then
           chit = 0.0_wp
           chip = 0.0_wp
           eta  = 0.0_wp

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(mxshak,tolnce)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              Do
                 Call scale_temperature(sigma,degtra,degrot,degfre)
                 Call core_shell_quench(safe,stptmp)
                 If (megcon > 0) Call constraints_quench(mxshak,tolnce)
                 If (megpmf > 0) Call pmf_quench(mxshak,tolnce)
                 If (megrgd > 0) Call rigid_bodies_quench()
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(sigma,degtra,degrot,degfre)
           End If

! Correct kinetic stress and energy

           If (megrgd > 0) Then
              Call kinstresf(vxx,vyy,vzz,strknf)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz)
           Else
              Call kinstress(vxx,vyy,vzz,strkin)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If


!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
