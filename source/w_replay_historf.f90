!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Report work

  If (idnode == 0) Write(nrite,"(/,3(/,1x,a))") &
     '*** HISTORF is replayed in full as statics (no dynamics)!!! ***', &
     '*** Big frame differences as particle move distance within HISTROF ***', &
     '*** and w.r.t. CONFIG at start may lead failures in parallel !!! ***'

! defect detection for every entry in HISTORF

  nsdef = 0
  isdef = 1

! MSDTMP option for every entry in HISTORF

  nstmsd = 0
  istmsd = 1

! displacement detection for every entry in HISTORF

  nsrsd = 0
  isrsd = 1

! intramolecular PDF analysis for every entry in HISTORF
! enforce printing and collection if the calculation exists

  lpana=(mxgana > 0)
  nstbnd = 1
  nstang = 1
  nstdih = 1
  nstinv = 1

! rdf and z-density detection for every entry in HISTORF
! enforce printing and collection if the calculation exists

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
  nstph = 0 ! HISTORF trajectory points counter
  Do
     Call allocate_statistics_connect()
10   Continue
     If (nstep >= nstpe) Call statistics_connect_set(imcon,rlnk)

! Make a move - Read a frame

     Call read_history(l_str,"HISTORF",megatm,levcfg,imcon,dvar,nstep,tstep,time,exout)

     If (newjb) Then
        newjb = .false.

        tmst=time
        tmsh=0.0_wp ! tmst substitute
     End If

     If (nstep >= 0 .and. exout >= 0) Then
        nstph=nstph+1

        If (nstep < nstpe) Go To 10 ! Deal with restarts

! First frame positions (for estimates of MSD when levcfg==0)

        If (nstph == 1) Then
           Do i=1,natms
              xin(i)=xxx(i)
              yin(i)=yyy(i)
              zin(i)=zzz(i)
           End Do
        End If

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
           Call build_book_intra        &
           (lsim,dvar,                  &
           megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,        &
           megrgd,degrot,degtra,        &
           megtet,megbnd,megang,megdih,meginv)
           If (lexcl) Call build_excl_intra(lecx)
        End If

! Evaluate forces, newjob must always be true for vircom evaluation

        Call w_calculate_forces()

! Evaluate kinetics if available

        If (levcfg > 0 .and. levcfg < 3) Then
           If (lzero .and. nstep <= nsteql) Call zero_k_optimise(strkin,strknf,strknt,engke,engrot)

! Calculate kinetic stress and energy if available

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

! Apply kinetic options

           Call w_kinetic_options()
        End If

! Get complete stress tensor

        strtot = strcon + strpmf + stress + strkin + strcom + strdpd

! Get core-shell kinetic energy for adiabatic shell model

        If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke)

! Calculate physical quantities and collect statistics,
! accumulate z-density if needed
! (nstep->nstph,tstep->tsths,tmst->tmsh)

        tsths=Max(tstep ,(time-tmsh) / Real(Merge( nstph-1, 1, nstph > 2), wp))

! Collect VAF if kinetics is available

        If (levcfg > 0 .and. levcfg < 3) Call vaf_collect(lvafav,leql,nsteql,nstph-1,time)

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

! line-printer output
! Update cpu time

        Call gtime(timelp)

        If (idnode == 0) Then
           If (Mod(lines,npage) == 0) Write(nrite,"(1x,130('-'),/,/,    &
              & 10x,'step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',     &
              & 5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',    &
              & 5x,'eng_dih',5x,'eng_tet',/,6x,'time(ps)',5x,' eng_pv', &
              & 4x,'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',   &
              & 5x,'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,  &
              & 6x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',   &
              & 5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',           &
              & 5x,'vir_pmf',7x,'press',/,/,1x,130('-'))")

           Write(nrite,"(1x,i13,1p,9e12.4,/,0p,f14.5,1p,9e12.4,         &
                & /,1x,0p,f13.3,1p,9e12.4)") nstep,(stpval(i),i=1,9),   &
                time,(stpval(i),i=10,18),timelp,(stpval(i),i=19,27)

           Write(nrite,"(/,7x,'rolling',1p,9e12.4,/,6x,'averages',      &
                & 1p,9e12.4,/,14x,9e12.4)") (ravval(i),i=1,27)

           Write(nrite,"(1x,130('-'))")
        End If

        If (nstph /= 0) lines=lines+1

! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types if kinetics exists

        If (ltraj) Call trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)
        If (ldef) Call defects_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)
        If (l_msd) Call msd_write &
           (keyres,nstmsd,istmsd,megatm,nstep,tstep,time)
        If (lrsd) Call rsd_write &
           (imcon,keyres,nsrsd,isrsd,rrsd,nstep,tstep,time)
        If (levcfg > 0 .and. levcfg < 3) Then
           If (vafsamp > 0) Call vaf_write & ! (nstep->nstph,tstep->tsths,tmst->tmsh)
           (lvafav,keyres,nstph,tsths)
        End If

! Save restart data in event of system crash

        If (Mod(nstph,ndump) == 0 .and. nstph /= nstrun .and. (.not.l_tor)) &
           Call system_revive                                       &
           (imcon,rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

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

! Save last frame positions (for estimates of MSD when levcfg==0)

     Do i=1,natms
        xin(i)=xxx(i)
        yin(i)=yyy(i)
        zin(i)=zzz(i)
     End Do
  End Do

! If reading HISTORY finished awkwardly
! recover positions and generate kinetics

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
  tstep=tsths


!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
