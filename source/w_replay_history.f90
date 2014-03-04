!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORY INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Report work

  If (idnode == 0) &
     Write(nrite,"(/,/,1x,'*** HISTORY is replayed for recalculation of structural properties ***')")

! Stay safe

  If (ltraj) Then
     ltraj = .false.

     If (idnode == 0) &
     Write(nrite,"(/,/,1x,'*** warning - aborting printing into HISTORY while reading it ***')")
  End If

! Make sure of no equilibration

  nsteql = 0
  leql   = .false.

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

! MSDTMP option for every entry in HISTORY

  nstmsd = 0
  istmsd = 1

! displacement detection for every entry in HISTORY

  nsrsd = 0
  isrsd = 1

! rdf and z-density detection for every entry in HISTORY
! impose printing if calculation exists

  lprdf=lrdf ; nstrdf = 1
  lpzdn=lzdn ; nstzdn = 1

! Calculate kinetic tensor and energy at restart as it may not exists later

  If (megrgd > 0) Then
     Call kinstresf(vxx,vyy,vzz,strknf)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

     strkin=strknf+strknt
  Else
     Call kinstress(vxx,vyy,vzz,strkin)
  End If
  engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  nstpe = nstep
  nstph = 0 ! trajectory points counter
  Do
     Call allocate_statistics_connect()
10   Continue
     If (nstep >= nstpe) Call statistics_connect_set(imcon,rlnk)

! Read a frame

     Call read_history(l_str,"HISTORY",megatm,levcfg,imcon,dvar,nstep,tstep,time,exout)

     If (newjb) Then
        newjb = .false.

        tmst=time
        tmsh=0.0_wp ! tmst substitute
     End If

     If (nstep >= 0 .and. exout >= 0) Then
        nstph=nstph+1

        If (nstep < nstpe) Go To 10 ! Deal with restarts

! CHECK MD CONFIGURATION

        Call check_config &
           (levcfg,imcon,l_str,lpse,keyens,iso,keyfce,keyres,megatm)

! get xto/xin/msdtmp arrays sorted

        Call statistics_connect_frames()
        Call deallocate_statistics_connect()

! SET domain borders and link-cells as default for new jobs
! exchange atomic data and positions in border regions

        Call set_halo_particles(imcon,rlnk,keyfce)

! For any intra-like interaction, construct book keeping arrays and
! exclusion arrays for overlapped two-body inter-like interactions

        If (lbook) Then
           Call build_book_intra              &
           (lsim,megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,              &
           megrgd,degrot,degtra,              &
           megtet,megbnd,megang,megdih,meginv)
           If (lexcl) Call build_excl_intra(lecx)
        End If

! Accumulate RDFs if needed (nstep->nstph)
! Make sure RDFs are complete (lbook=.false. - no exclusion lists)

        If (lrdf) Call two_body_forces            &
           (imcon,rcut,rlnk,rvdw,rmet,keyens,     &
           alpha,epsq,keyfce,nstfce,.false.,megfrz, &
           lrdf,nstrdf,leql,nsteql,nstph,         &
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
              If (levcfg == 2) Then
                 Call rigid_bodies_str_ss(strcom)
                 vircom=-(strcom(1)+strcom(5)+strcom(9))
              End If
           Else
              Call kinstress(vxx,vyy,vzz,strkin)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If

! Calculate physical quantities and collect statistics,
! accumulate z-density if needed
! (nstep->nstph,tstep->tsths,tmst->tmsh)

        tsths=Max(tstep ,(time-tmsh) / Real(Merge( nstph-1, 1, nstph > 2), wp))

        Call statistics_collect            &
           (leql,nsteql,lzdn,nstzdn,       &
           keyres,keyens,iso,intsta,imcon, &
           degfre,degshl,degrot,           &
           nstph,tsths,time,tmsh,          &
           engcpe,vircpe,engsrp,virsrp,    &
           engter,virter,                  &
           engtbp,virtbp,engfbp,virfbp,    &
           engshl,virshl,shlke,            &
           vircon,virpmf,                  &
           engtet,virtet,engfld,virfld,    &
           engbnd,virbnd,engang,virang,    &
           engdih,virdih,enginv,virinv,    &
           engke,engrot,consv,vircom,      &
           strtot,press,strext,            &
           stpeng,stpvir,stpcfg,stpeth,    &
           stptmp,stpprs,stpvol)

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        Call w_write_options()

! Complete time check

        Call gtime(timelp)
        If (idnode == 0) Then
           Write(nrite,"(1x,'HISTORY step',i10,' (',i10,' entry) processed')") nstep,nstph
           Write(nrite,'(1x,"time elapsed since job start: ", f12.3, " sec")') timelp
        End If

! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

        i=20
        If ( Mod(nstph,i) == 0 .or.                               &
             (timelp > Real(i*60,wp) .and.                        &
              timelp-Real( ((Int(timelp)/(i*60)) * i*60) , wp ) < &
              timelp/Real( nstph , wp) ) ) Then

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

        If (exout /= 0) Exit
     Else
        Exit
     End If
  End Do

! If reading HISTORY finished awkwardly

  If (exout < 0) Then
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

  nstep=nstph


!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORY INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
