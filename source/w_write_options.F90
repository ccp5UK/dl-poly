!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

If (traj%ltraj) Then
  Call trajectory_write(flow%restart_key,megatm,flow%step,thermo%tstep,flow%time,io,stat%rsd,netcdf,config, &
    traj,files,comm)
End If

If(dfcts(1)%ldef) Then
  Call defects_write(flow%restart_key,thermo%ensemble,flow%step,thermo%tstep,flow%time,io,cshell,dfcts(1), &
    neigh,sites,netcdf,domain,config,files,comm)
  If (dfcts(2)%ldef) Then
    Call defects_write(flow%restart_key,thermo%ensemble,flow%step,thermo%tstep,flow%time,io,cshell,dfcts(2), &
      neigh,sites,netcdf,domain,config,files,comm)
  End If
End If

If (msd_data%l_msd) Then
  Call msd_write(config,flow%restart_key,megatm,flow%step,thermo%tstep,flow%time,stat%stpval,sites%dof_site, &
    io,msd_data,files,comm)
End If

If (rsdc%lrsd) Then
  Call rsd_write(flow%restart_key,flow%step,thermo%tstep,io,rsdc,flow%time,cshell,stat%rsd,config,comm)
End If

If (green%samp > 0) Then
  Call vaf_write(config,flow%restart_key,flow%step,thermo%tstep,green,sites,comm)
End If

!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
