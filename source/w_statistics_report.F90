!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Get complete stress tensor

strtot = strcon + strpmf + stress + strkin + strcom + strdpd

! Get core-shell kinetic energy for adiabatic shell model

If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke,comm)

! Calculate physical quantities and collect statistics

Call statistics_collect           &
  (lsim,leql,nsteql,lzdn,nstzdn, &
  keyres,keyens,      &
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
  strtot,           &
  stpeng,stpvir,stpcfg,stpeth,   &
  stptmp,stpprs,stpvol,stat,thermo,comm,virdpd)

! VV forces evaluation report for 0th or weird restart

If (l_vv .and. levcfg == 1) Then
  Call info('forces evaluated at (re)start for VV integration...',.true.)
  Call info('',.true.)
End If

! line-printer output every nstbpo steps

If (lines == 0 .or. Mod(nstep,nstbpo) == 0) Then

  ! Update cpu time

  Call gtime(tmr%elapsed)

  If (Mod(lines,npage) == 0) Then
    Write(messages(1),'(a)') Repeat('-',130)
    Write(messages(2),'(9x,a4,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
     'step','eng_tot','temp_tot','eng_cfg','eng_src','eng_cou','eng_bnd','eng_ang','eng_dih','eng_tet'
    Write(messages(3),'(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
     'time(ps)',' eng_pv','temp_rot','vir_cfg','vir_src','vir_cou','vir_bnd','vir_ang','vir_con','vir_tet'
    Write(messages(4), '(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,7x,a5,8x,a4,7x,a5,5x,a7,7x,a5)') &
      'cpu  (s)','volume','temp_shl','eng_shl','vir_shl','alpha','beta','gamma','vir_pmf','press'
    Write(messages(5),'(a)') Repeat('-',130)
    Call info(messages,5,.true.)
  End If

  Write(messages(1),'(i13,1p,9e12.4)')nstep,stat%stpval(1:9)
  Write(messages(2),'(f13.5,1p,9e12.4)')time,stat%stpval(10:18)
  Write(messages(3),'(0p,f13.3,1p,9e12.4)') tmr%elapsed,stat%stpval(19:27)
  Write(messages(4),'(a)')''
  Call info(messages,4,.true.)

  Write(messages(1),'(6x,a7,1p,9e12.4)') 'rolling',stat%ravval(1:9)
  Write(messages(2),'(5x,a8,1p,9e12.4)') 'averages',stat%ravval(10:18)
  Write(messages(3),'(13x,9e12.4)') stat%ravval(19:27)
  Write(messages(4),'(a)') Repeat('-',130)
  Call info(messages,4,.true.)

  If (nstep /= 0) lines=lines+1

End If

! Reports at end of equilibration period

If (nstep == nsteql) Then

  If (nstep > 0) Then
    Call info(repeat('-',130),.true.)
    If (thermo%l_zero) Then
      thermo%l_zero=.false.
      Write(message,'(a,i10)') 'switching off zero Kelvin optimiser at step ',nstep
      Call info(message,.true.)
    End If

    If (lmin) Then
      lmin=.false.
      Write(message,'(a,i10)') 'switching off CGM minimiser at step ',nstep
      Call info(message,.true.)
    End If

    If (thermo%l_tscale) Then
      thermo%l_tscale=.false.
      Write(message,'(a,i10)') 'switching off temperature scaling at step ',nstep
      Call info(message,.true.)
    End If

    If (thermo%l_tgaus) Then
      thermo%l_tgaus=.false.
      Write(message,'(a,i10)') 'switching off temperature regaussing at step ',nstep
      Call info(message,.true.)
    End If
  End If

  ! bond & PMF constraint quenching iterative cycles statistics report

  If (megcon > 0) Then
    Call gmax(comm,passcnq(3:5))
    If (passcnq(3) > 0.0_wp) Then
      Write(message,'(2(a,f5.2),4(a,i3))') &
        'constraints quench run statistics per call: average cycles ', &
        passcnq(3),'/',passcnq(3), &
        ' minimum cycles ',Nint(passcnq(4)),'/',Nint(passcnq(4)), &
        ' maximum cycles ',Nint(passcnq(5)),'/',Nint(passcnq(5))
      Call info(message,.true.)
    End If
  End If

  If (megpmf > 0) Then
    Call gmax(comm,passpmq(3:5))
    If (passpmq(3) > 0.0_wp) Then
      Write(message,'(2(a,f5.2),4(a,i3))') &
        'PMFs quench run statistics per call: average cycles ', &
        passpmq(3),'/',passpmq(3), &
        ' minimum cycles ',Nint(passpmq(4)),'/',Nint(passpmq(4)), &
        ' maximum cycles ',Nint(passpmq(5)),'/',Nint(passpmq(5))
      Call info(message,.true.)
    End If
  End If

  If (nstep > 0 .or. megcon > 0 .or. megpmf > 0) Then
    Call info(repeat('-',130),.true.)
  End If
End If

! Calculate green-kubo properties

If (vafsamp > 0) Call vaf_collect(lvafav,leql,nsteql,nstep,time,comm)


!!!!!!!!!!!!!!!!!  W_STATISTICS_REPORT INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
