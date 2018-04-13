!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Get complete stress tensor

strtot = strcon + strpmf + stress + strkin + strcom + strdpd

! Get core-shell kinetic energy for adiabatic shell model

If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke,comm)

! Calculate physical quantities and collect statistics

Call statistics_collect           &
  (lsim,leql,nsteql,lzdn,nstzdn, &
  keyres,keyens,iso,intsta,      &
  degfre,degshl,degrot,          &
  nstep,tstep,time,tmst,         &
  engcpe,vircpe,engsrp,virsrp,   &
  engter,virter,                 &
  engtbp,virtbp,engfbp,virfbp,   &
  engshl,virshl,shlke,           &
  vircon,virpmf,                 &
  engtet,virtet,engfld,virfld,   &
  engbnd,virbnd,engang,virang,   &
  engdih,virdih,enginv,virinv,   &
  engke,engrot,consv,vircom,     &
  strtot,press,strext,           &
  stpeng,stpvir,stpcfg,stpeth,   &
  stptmp,stpprs,stpvol,comm,virdpd)

! VV forces evaluation report for 0th or weird restart

If (l_vv .and. levcfg == 1 .and. comm%idnode == 0) &
  Write(nrite,'(1x,a,/)') "forces evaluated at (re)start for VV integration..."

! line-printer output every nstbpo steps

If (lines == 0 .or. Mod(nstep,nstbpo) == 0) Then

  ! Update cpu time

  Call gtime(timelp)

  If (Mod(lines,npage) == 0) Then 
    Write(messages(1),"(1x,130('-'))")
    Write(messages(2),"(10x,'step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',     &
      & 5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',    &
      & 5x,'eng_dih',5x,'eng_tet')")
    Write(messages(3),"(6x,'time(ps)',5x,' eng_pv', &
      & 4x,'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',   &
      & 5x,'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet')")
    Write(messages(4), "(6x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',   &
      & 5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',           &
      & 5x,'vir_pmf',7x,'press')")
    Write(messages(5),"(1x,130('-'))")
    Call info(messages,5,.true.)
  End IF

  Write(messages(1),'(1x,i13,1p,9e12.4)')nstep,stpval(1:9)
  Write(messages(2),'(f14.5,1p,9e12.4)')time,stpval(10:18)
  Write(messages(3),'(1x,0p,f13.3,1p,9e12.4)') timelp,stpval(19:27)
  Write(messages(4),'(a)')''
  Call info(messages,4,.true.)

  Write(messages(1),"(7x,'rolling',1p,9e12.4)") ravval(1:9)
  Write(messages(2),"(6x,'averages',1p,9e12.4)") ravval(10:18)
  Write(messages(3),"(14x,9e12.4)") ravval(19:27)
  Write(messages(4),"(1x,130('-'))")
  Call info(messages,4,.true.)

  If (nstep /= 0) lines=lines+1

End If

! Reports at end of equilibration period

If (nstep == nsteql) Then

  If (nstep > 0) Then
    If (lzero) Then
      lzero=.false.
      If (comm%idnode == 0) Write(nrite,"(/,1x,a,i0)") &
        'switching off zero Kelvin optimiser at step ',nstep
    End If

    If (lmin) Then
      lmin=.false.
      If (comm%idnode == 0) Write(nrite,"(/,1x,a,i0)") &
        'switching off CGM minimiser at step ',nstep
    End If

    If (ltscal) Then
      ltscal=.false.
      If (comm%idnode == 0) Write(nrite,"(/,1x,a,i0)") &
        'switching off temperature scaling at step ',nstep
    End If

    If (ltgaus) Then
      ltgaus=.false.
      If (comm%idnode == 0) Write(nrite,"(/,1x,a,i0)") &
        'switching off temperature regaussing at step ',nstep
    End If
  End If

  ! bond & PMF constraint quenching iterative cycles statistics report

  If (megcon > 0) Then
    Call gmax(comm,passcnq(3:5))
    If (passcnq(3) > 0.0_wp .and. comm%idnode == 0) Write(nrite,"(/,                             &
      & ' constraints quench run statistics per call: average cycles ', f5.2, ' / ', f5.2, &
      & ' minimum cycles ', i3, ' / ', i3, ' maximum cycles ', i3, ' / ', i3)")          &
      passcnq(3),passcnq(3),Nint(passcnq(4)),Nint(passcnq(4)),Nint(passcnq(5)),Nint(passcnq(5))
  End If

  If (megpmf > 0) Then
    Call gmax(comm,passpmq(3:5))
    If (passpmq(3) > 0.0_wp .and. comm%idnode == 0) Write(nrite,"(/,                      &
      & ' PMFs quench run statistics per call: average cycles ', f5.2, ' / ', f5.2, &
      & ' minimum cycles ', i3, ' / ', i3, ' maximum cycles ', i3, ' / ', i3)")   &
      passpmq(3),passpmq(3),Nint(passpmq(4)),Nint(passpmq(4)),Nint(passpmq(5)),Nint(passpmq(5))
  End If

  If ((nstep > 0 .or. megcon > 0 .or. megpmf > 0) .and. comm%idnode == 0) &
    Write(nrite,"(/,1x,130('-'))")
End If

! Calculate green-kubo properties

If (vafsamp > 0) Call vaf_collect(lvafav,leql,nsteql,nstep,time,comm)


!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
