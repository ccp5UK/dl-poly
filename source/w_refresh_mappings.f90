!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Scale t=0 reference positions

        If (nstep > 0) Call xscale(m_rgd,keyens,tstep,eta)

! Check VNL conditioning

        Call vnl_check(l_str,rcut,rpad,rlnk,width)

        If (l_vnl) Then

! Relocate atoms to new domains and restore bonding description

           Call relocate_particles  &
           (dvar,rlnk,lbook,megatm, &
           megshl,m_con,megpmf,     &
           m_rgd,megtet,            &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

           Call set_halo_particles(rlnk,keyfce) ! inducing in here only

! Re-tag RBs when called again after the very first time
! when it's done in rigid_bodies_setup <- build_book_intra

           If (m_rgd > 0) Then
              Call rigid_bodies_tags()
              Call rigid_bodies_coms(xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz)
           End If

        Else

! Exchange atomic positions in border regions

           Call refresh_halo_positions()

        End If

! set and halo rotational matrices and their infinitesimal rotations

        If (mximpl > 0) Call mpoles_rotmat_set_halo()

!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
