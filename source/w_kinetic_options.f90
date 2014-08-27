!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Apply external field

        If (keyfld > 0) Call external_field_correct(imcon)

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

           Call regauss_temperature(imcon,megrgd)

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              stptmp = 2.0_wp*(engke+engrot) / (boltz*Real(degfre,wp))
              Do
                 Call scale_temperature(imcon,engke+engrot,degtra,degrot,degfre)
                 Call core_shell_quench(safe,stptmp)
                 If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
                 If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)
                 If (megrgd > 0) Call rigid_bodies_quench(imcon)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(imcon,engke+engrot,degtra,degrot,degfre)
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

           If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              Do
                 Call scale_temperature(imcon,sigma,degtra,degrot,degfre)
                 Call core_shell_quench(safe,stptmp)
                 If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
                 If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)
                 If (megrgd > 0) Call rigid_bodies_quench(imcon)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(imcon,sigma,degtra,degrot,degfre)
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

! bond & PMF onstraint quenching iterative cycles statistics report

        If (nstep == nsteql) Then
           If (megcon > 0) Then
              Call gmax(passcnq(3:5))
              If (passcnq(3) > 0.0_wp .and. idnode == 0) Write(nrite,"(//,                            &
                 & ' constraints quench run statistics per call: average cycles ', f5.2, ' / ', f5.2, &
                 & ' minimum cycles ', i3, ' / ', i3, ' maximum cycles ', i3, ' / ', i3)")            &
                 passcnq(3),passcnq(3),Nint(passcnq(4)),Nint(passcnq(4)),Nint(passcnq(5)),Nint(passcnq(5))
           End If

           If (megpmf > 0) Then
              Call gmax(passpmq(3:5))
              If (passpmq(3) > 0.0_wp .and. idnode == 0) Write(nrite,"(//,                     &
                 & ' PMFs quench run statistics per call: average cycles ', f5.2, ' / ', f5.2, &
                 & ' minimum cycles ', i3, ' / ', i3, ' maximum cycles ', i3, ' / ', i3)")     &
                 passpmq(3),passpmq(3),Nint(passpmq(4)),Nint(passpmq(4)),Nint(passpmq(5)),Nint(passpmq(5))
           End If
        End If


!!!!!!!!!!!!!!!!!!!  W_KINETIC_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
