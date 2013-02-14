!!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Apply impact
! levcfg == 2 avoids application twice when tmd happens at (re)start for VV

     If (nstep == tmd .and. levcfg == 2) Then
        If (idnode == 0) Write(nrite,"(/,          &
           & /,1x,'initiating IMPACT:',            &
           & /,1x,46('-'),                         &
           & /,1x,'particle (index)',15x,i10,      &
           & /,1x,'timestep (steps)',15x,i10,      &
           & /,1x,'energy   (keV)  ',18x,1p,e12.4, &
           & /,1x,'v-r(x,y,z)',1p,3e12.4,          &
           & /,1x,46('-'))") imd,tmd,emd,vmx,vmy,vmz

        If (nstep+1 <= nsteql) Call warning(380,Real(nsteql,wp),0.0_wp,0.0_wp)

        Call impact(imd,emd,vmx,vmy,vmz,megrgd)

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


!!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
