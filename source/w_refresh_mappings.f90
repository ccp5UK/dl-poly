!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Relocate atoms to new domains and restore bonding description

        Call relocate_particles      &
           (imcon,rcut,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

        Call set_halo_particles(imcon,rcut,keyfce)

! Re-tag RBs when called again after the very first time
! when it's done in rigid_bodies_setup <- build_book_intra

        If (m_rgd > 0) Then
           Call rigid_bodies_tags()
           Call rigid_bodies_coms(imcon,xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz)
        End If


!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
