!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

        If (ltraj) Call trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)
        If (ldef) Call defects_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)
        If (l_msd) Call msd_write &
           (keyres,nstmsd,istmsd,megatm,nstep,tstep,time)
        If (lrsd) Call rsd_write &
           (imcon,keyres,nsrsd,isrsd,rrsd,nstep,tstep,time)
        If (vafsamp > 0) Call vaf_write &
           (lvafav,keyres,leql,nsteql,nstep,tstep)


!!!!!!!!!!!!!!!!!!!!!  W_WRITE_OPTIONS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
