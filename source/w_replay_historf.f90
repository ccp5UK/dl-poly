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

! rdf and z-density detection for every entry in HISTORF
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
  nstph = 0 ! HISTORF trajectory points counter
  Do
10 Continue

! Make a move

     Call read_history(l_str,"HISTORF",megatm,levcfg,imcon,nstep,tstep,time)

     If (newjb) Then
        newjb = .false.

        tmst=time
        tmsh=0.0_wp ! tmst substitute
     End If

     If (nstep >= 0) Then
        nstph=nstph+1

        If (nstep <= nstpe) Go To 10 ! Deal with restarts

! Refresh mappings

        Call w_refresh_mappings()

! Evaluate forces, newjob must always be true for vircom evaluation

        Call w_calculate_forces()

! Evaluate kinetics

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

        strtot = strcon + strpmf + stress + strkin + strcom

! Get core-shell kinetic energy for adiabatic shell model

        If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke)

! Calculate physical quantities and collect statistics,
! accumulate z-density if needed (nstep->nstph,tmst->tmsh)

        Call statistics_collect            &
           (leql,nsteql,lzdn,nstzdn,       &
           keyres,keyens,iso,intsta,imcon, &
           degfre,degshl,degrot,           &
           nstph,tstep,time,tmsh,          &
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
              & 5x,'step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',      &
              & 5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',    &
              & 5x,'eng_dih',5x,'eng_tet',/,1x,'time(ps)',5x,' eng_pv', &
              & 4x,'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',   &
              & 5x,'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,  &
              & 1x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',   &
              & 5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',           &
              & 5x,'vir_pmf',7x,'press',/,/,1x,130('-'))")

           Write(nrite,"(1x,i8,1p,9e12.4,/,0p,f9.5,1p,9e12.4,           &
                & /,1x,0p,f8.1,1p,9e12.4)") nstep,(stpval(i),i=1,9),    &
                time,(stpval(i),i=10,18),timelp,(stpval(i),i=19,27)

           Write(nrite,"(/,2x,'rolling',1p,9e12.4,/,1x,'averages',      &
                & 1p,9e12.4,/,9x,9e12.4)") (ravval(i),i=1,27)

           Write(nrite,"(1x,130('-'))")
        End If

        If (nstph /= 0) lines=lines+1

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        Call w_write_options()

! Save restart data in event of system crash

        If (Mod(nstph,ndump) == 0 .and. nstph /= nstrun .and. (.not.l_tor)) &
           Call system_revive                                       &
           (imcon,rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

! Close and Open OUTPUT at about 'i'th print-out or 'i' minunte intervals

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
     Else
        Exit
     End If
  End Do

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

  nstep=nstph


!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!