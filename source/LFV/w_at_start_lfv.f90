!!!!!!!!!!!!!!!!!!!!!!  w_AT_START_LFV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate kinetic stress tensor and energy at t=0 for reporting

  If (nstep == 0) Then
     If (megrgd > 0) Then
        Call kinstresf(vxx,vyy,vzz,strknf)
        Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

        strkin=strknf+strknt
     Else
        Call kinstress(vxx,vyy,vzz,strkin)
     End If
     engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
  End If

! Fix levcfg=2 as forces are always calculated at (re)start in LFV

  If (levcfg == 1) levcfg=2


!!!!!!!!!!!!!!!!!!!!!!  w_AT_START_LFV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
