!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

        If (ltraj) Call trajectory_write &
           (keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,stat%rsd,comm)
        If(dfcts(1)%ldef)Then
           Call defects_write &
          (keyres,thermo%ensemble,nstep,tstep,time,cshell,dfcts(1),neigh,site,comm)
           If (dfcts(2)%ldef)Then
             Call defects_write &
             (keyres,thermo%ensemble,nstep,tstep,time,cshell,dfcts(2),neigh,site,comm)
           End If
        End If   
        If (msd_data%l_msd) Call msd_write &
          (keyres,megatm,nstep,tstep,time,stat%stpval,site%dof_site,msd_data,comm)
        If (lrsd) Call rsd_write &
           (keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,cshell,stat%rsd,comm)
        If (green%samp > 0) Call vaf_write &
           (keyres,nstep,tstep,green,site,comm)


!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
