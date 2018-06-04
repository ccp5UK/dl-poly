!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Scale t=0 reference positions

        If (nstep > 0) Call xscale(m_rgd,tstep,thermo,stat,comm)

! Check VNL conditioning

        Call vnl_check(l_str,rcut,rpad,rlnk,width,comm)

        If (l_vnl) Then

! Relocate atoms to new domains and restore bonding description

           Call relocate_particles  &
           (dvar,rlnk,lbook,msd_data%l_msd,megatm, &
           megshl,m_con,megpmf,     &
           m_rgd,megtet,            &
           megang,megdih,meginv,stat,ewld,thermo,green,bond,comm)

! Exchange atomic data in border regions

           Call set_halo_particles(rlnk,keyfce,comm) ! inducing in here only

! Re-tag RBs when called again after the very first time
! when it's done in rigid_bodies_setup <- build_book_intra

           If (m_rgd > 0) Then
              Call rigid_bodies_tags(comm)
              Call rigid_bodies_coms(xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz,comm)
           End If

        Else

! Exchange atomic positions in border regions

           Call refresh_halo_positions(comm)

        End If

! set and halo rotational matrices and their infinitesimal rotations

        If (mximpl > 0) Call mpoles_rotmat_set_halo(comm)

!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
