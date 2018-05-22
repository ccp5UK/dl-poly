!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

        If (ltraj) Call trajectory_write &
           (keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,comm)
        If (ldef) Call defects_write &
           (rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time,comm)
        If (l_msd) Call msd_write &
           (keyres,nstmsd,istmsd,megatm,nstep,tstep,time,stpval,comm)
        If (lrsd) Call rsd_write &
           (keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,comm)
        If (green%samp > 0) Call vaf_write &
           (lvafav,keyres,nstep,tstep,green,comm)


!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
