!!!!!!!!!!!!!!!!!!!!!  REPLAY HISTORY INCLUSION MODULE  !!!!!!!!!!!!!!!!


  If (idnode == 0) &
     Write(nrite,"(/,/,1x,'*** HISTORY is replayed for recalculation of structural properties ***')")

! MSD option cancelled

  If (l_msd) Then
     l_msd=.false.
     If (idnode == 0) Write(nrite,"(/,1x,'*** MSDTMP file option canceled ***')")
  End If

! Make sure RDFs are complete (no exclusion lists)

  lbook=.false.

! nullify forces

  fxx = 0.0_wp
  fyy = 0.0_wp
  fzz = 0.0_wp

! nullify all two-body force switches = just do rdf calculation

  keyfce = 0
  ntpvdw = 0
  ntpmet = 0

! defect detection for every entry in HISTORY

  nsdef = 0
  isdef = 1

! rdf and z-density detection for every entry in HISTORY

  nsteql = 0
  leql   = .false.
  nstrdf = 1
  nstzdn = 1

  isw = 0 ! trajectory points counter
  Do While ( .true. )
     Call read_history(l_str,megatm,levcfg,imcon,nstep,tstep,time)

     If (newjob) Then
        newjob = .false.
        tmst=time
        tmsth=0.0_wp ! tmst substitute
     End If

     If (nstep >= 0) Then

! Exchange atomic data in border regions

        Call set_halo_particles(imcon,rcut,keyfce,lbook)

        isw=isw+1

! Accumulate RDFs if needed (nstep->isw)

        If (lrdf) Call two_body_forces            &
           (imcon,rcut,rvdw,rmet,keyens,          &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           lrdf,nstrdf,leql,nsteql,isw,           &
           elrc,virlrc,elrcm,vlrcm,               &
           engcpe,vircpe,engsrp,virsrp,stress)

! Calculate kinetic stress and energy if available

        If (levcfg > 0 .and. levcfg < 3) Then
           If (megrgd > 0) Then
              Call rigid_bodies_quench(imcon)

              Call kinstresf(vxx,vyy,vzz,strknf)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz)
           Else
              Call kinstress(vxx,vyy,vzz,strkin)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        Else
           If (megrgd > 0) Then
              rgdvxx(1:ntrgd) = 0.0_wp
              rgdvyy(1:ntrgd) = 0.0_wp
              rgdvzz(1:ntrgd) = 0.0_wp

              rgdoxx(1:ntrgd) = 0.0_wp
              rgdoyy(1:ntrgd) = 0.0_wp
              rgdozz(1:ntrgd) = 0.0_wp

              strknt = 0.0_wp
              engrot = 0.0_wp
           End If
           vxx(1:natms) = 0.0_wp
           vyy(1:natms) = 0.0_wp
           vzz(1:natms) = 0.0_wp

           strknf = 0.0_wp
           strkin = 0.0_wp
           engke = 0.0_wp
        End If

! Calculate physical quantities and collect statistics,
! accumulate z-density if needed (nstep->isw,tmst->tmsth)

        Call statistics_collect              &
           (leql,nsteql,lzdn,nstzdn,         &
           keyres,keyens,intsta,imcon,       &
           degfre,degshl,degrot,             &
           isw,tstep,time,tmsth,             &
           engcpe,vircpe,engsrp,virsrp,      &
           engter,virter,                    &
           engtbp,virtbp,engfbp,virfbp,      &
           engshl,virshl,shlke,              &
           vircon,virpmf,                    &
           engtet,virtet,engfld,virfld,      &
           engbnd,virbnd,engang,virang,      &
           engdih,virdih,enginv,virinv,      &
           engke,engrot,consv,vircom,strtot, &
           stpeng,stpvir,stpcfg,stpeth,      &
           stptmp,stpprs,stpvol)

! Write defects data

        If (ldef) Call defects_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)

! Complete time check

        Call gtime(timelp)
        If (idnode == 0) Then
           Write(nrite,"(1x,'HISTORY step',i10,' (',i10,' entry) processed')") nstep,isw
           Write(nrite,'(1x,"time elapsed since job start: ", f12.3, " sec")') timelp
        End If

! Close and Open OUTPUT at about 'i'th print-out or 'i' minunte intervals

        i=20
        If ( Mod(isw,i) == 0 .or.                                 &
             (timelp > Real(i*60,wp) .and.                        &
              timelp-Real( ((Int(timelp)/(i*60)) * i*60) , wp ) < &
              timelp/Real( isw , wp) ) ) Then

           If (idnode == 0) Then
              Inquire(File='OUTPUT', Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                 Close(Unit=nrite)
                 Open(Unit=nrite, File='OUTPUT', Position='append')
              End If
           End If

        End If
     Else
        Go To 10
     End If
  End Do

10 Continue

! If reading HISTORY finished awkwardly

  If (nstep == -2) Then
     Do i=1,natms
        xxx(i)=xin(i)
        yyy(i)=yin(i)
        zzz(i)=zin(i)
     End Do

     Call set_temperature             &
           (levcfg,imcon,temp,keyres, &
           lmin,nstep,nstrun,nstmin,  &
           mxshak,tolnce,keyshl,      &
           atmfre,atmfrz,             &
           megshl,megcon,megpmf,      &
           megrgd,degtra,degrot,      &
           degfre,degshl,sigma,engrot)
  End If

! step counter is data counter now, so statistics_result is triggered

  nstep=isw


!!!!!!!!!!!!!!!!!!!!!  REPLAY HISTORY INCLUSION MODULE  !!!!!!!!!!!!!!!!
