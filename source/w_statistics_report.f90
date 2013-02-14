!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Get complete stress tensor

        strtot = strcon + strpmf + stress + strkin + strcom

! Get core-shell kinetic energy for adiabatic shell model

        If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke)

! Calculate physical quantities and collect statistics

        Call statistics_collect              &
           (leql,nsteql,lzdn,nstzdn,         &
           keyres,keyens,iso,intsta,imcon,   &
           degfre,degshl,degrot,             &
           nstep,tstep,time,tmst,            &
           engcpe,vircpe,engsrp,virsrp,      &
           engter,virter,                    &
           engtbp,virtbp,engfbp,virfbp,      &
           engshl,virshl,shlke,              &
           vircon,virpmf,                    &
           engtet,virtet,engfld,virfld,      &
           engbnd,virbnd,engang,virang,      &
           engdih,virdih,enginv,virinv,      &
           engke,engrot,consv,vircom,        &
           strtot,press,strext,              &
           stpeng,stpvir,stpcfg,stpeth,      &
           stptmp,stpprs,stpvol)

! line-printer output every nstbpo steps

        If (lines == 0 .or. Mod(nstep,nstbpo) == 0) Then

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

           If (nstep == 0) Then ! Print start-up time
              Write(nrite,'(/,/,/,1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
           Else ! Exclude t=0 reporting
              lines=lines+1
           End If
        End If

! Report end of equilibration period

        If (ltscal .and. nstep == nsteql) Then
           ltscal=.false.
           If (idnode == 0) Write(nrite,"(/,/,1x,                         &
              & 'switching off temperature scaling at step ',i0,/,/,/,1x, &
              & 130('-'))") nstep
        End If


!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
