!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

        If (ltraj) Call trajectory_write &
           (keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,stat%rsd,comm)
        If (ldef) Call defects_write &
           (rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time,comm)
        If (msd_data%l_msd) Call msd_write &
           (keyres,megatm,nstep,tstep,time,stat%stpval,msd_data,comm)
        If (lrsd) Call rsd_write &
           (keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,stat%rsd,comm)
        If (green%samp > 0) Call vaf_write &
           (lvafav,keyres,nstep,tstep,green,comm)


!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
