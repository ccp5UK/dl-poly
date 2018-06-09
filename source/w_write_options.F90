!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

        If (ltraj) Call trajectory_write &
           (keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,stat%rsd,comm)
        If(dfcts(1)%ldef)Then
           Call defects_write &
          (neigh%cutoff,keyres,thermo%ensemble,nstep,tstep,time,dfcts(1),comm)
           If (dfcts(2)%ldef)Then
             Call defects_write &
             (neigh%cutoff,keyres,thermo%ensemble,nstep,tstep,time,dfcts(2),comm)
           End If
        End If   
        If (msd_data%l_msd) Call msd_write &
           (keyres,megatm,nstep,tstep,time,stat%stpval,msd_data,comm)
        If (lrsd) Call rsd_write &
           (keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,stat%rsd,comm)
        If (green%samp > 0) Call vaf_write &
           (keyres,nstep,tstep,green,comm)


!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
