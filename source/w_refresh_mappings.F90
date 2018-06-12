!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Scale t=0 reference positions

        If (nstep > 0) Call xscale(m_rgd,tstep,thermo,stat,neigh,comm)

! Check VNL conditioning

        Call vnl_check(l_str,width,neigh,comm)

        If (neigh%update) Then

! Relocate atoms to new domains and restore bonding description

           Call relocate_particles  &
           (dvar,neigh%cutoff_extended,lbook,msd_data%l_msd,megatm, &
           megshl,     &
           m_rgd,megtet,            &
           cons, pmf,& 
           stat,ewld,thermo,green,bond,angle,dihedral,inversion,tether,comm)

! Exchange atomic data in border regions

           Call set_halo_particles(keyfce,neigh,comm) ! inducing in here only

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
