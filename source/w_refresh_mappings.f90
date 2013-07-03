!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Scale t=0 reference positions

        Call xscale(imcon,m_rgd,keyens,tstep,eta)

! Relocate atoms to new domains and restore bonding description

        Call relocate_particles      &
           (imcon,rcut,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

        Call set_halo_particles(imcon,rcut,keyfce)


!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
