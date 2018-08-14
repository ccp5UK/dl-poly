!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

If (ltraj) Then
  Call trajectory_write(keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time, &
    stat%rsd,netcdf,parts,comm)
End If

If(dfcts(1)%ldef) Then
  Call defects_write(keyres,thermo%ensemble,nstep,tstep,time,cshell,dfcts(1), &
    neigh,sites,netcdf,domain,parts,comm)
  If (dfcts(2)%ldef) Then
    Call defects_write(keyres,thermo%ensemble,nstep,tstep,time,cshell,dfcts(2), &
      neigh,sites,netcdf,domain,parts,comm)
  End If
End If

If (msd_data%l_msd) Then
  Call msd_write(keyres,megatm,nstep,tstep,time,stat%stpval,sites%dof_site, &
    msd_data,comm)
End If

If (lrsd) Then
  Call rsd_write(keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,cshell,stat%rsd,parts,comm)
End If

If (green%samp > 0) Then
  Call vaf_write(keyres,nstep,tstep,green,sites,comm)
End If

!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
