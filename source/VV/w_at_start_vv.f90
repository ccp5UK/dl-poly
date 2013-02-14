!!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Calculate kinetic tensor and energy at restart

  If (megrgd > 0) Then
     Call kinstresf(vxx,vyy,vzz,strknf)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

     strkin=strknf+strknt
  Else
     Call kinstress(vxx,vyy,vzz,strkin)
  End If
  engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! If levcfg=2 and RBs are present, update forces on shared ones
! and get RB COM stress and virial at restart.  If levcfg<2
! forces are calculated at (re)start

  If (levcfg == 1) Then
     If (megrgd > 0) Then
        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz)

        If (l_lan .and. (.not.l_lan_s)) Then
           Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl)
        Else
           Call rigid_bodies_str_ss(strcom)
        End If

        vircom=-(strcom(1)+strcom(5)+strcom(9))
     End If
  End If


!!!!!!!!!!!!!!!!!!!!!!!  W_AT_START_VV INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
