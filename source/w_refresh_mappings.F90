!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Scale t=0 reference positions

If (nstep > 0) Call xscale(tstep,thermo,stat,neigh,rigid,domain,comm)

! Check VNL conditioning

Call vnl_check(l_str,width,neigh,stat,domain,parts,comm)

If (neigh%update) Then

  ! Relocate atoms to new domains and restore bonding description

  Call relocate_particles(dvar,neigh%cutoff_extended,lbook, &
    msd_data%l_msd,megatm,cshell,cons,pmf,stat,ewld,thermo,green, &
    bond,angle,dihedral,inversion,tether,neigh,sites,minim,mpoles, &
    rigid,domain,comm)

  ! Exchange atomic data in border regions

  Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,comm) ! inducing in here only

  ! Re-tag RBs when called again after the very first time
  ! when it's done in rigid_bodies_setup <- build_book_intra

  If (rigid%on) Then
    Call rigid_bodies_tags(rigid,comm)
    Call rigid_bodies_coms(parts,rigid%xxx,rigid%yyy,rigid%zzz,rigid,comm)
  End If

Else

  ! Exchange atomic positions in border regions

  Call refresh_halo_positions(domain,comm)
End If

! set and halo rotational matrices and their infinitesimal rotations

If (mpoles%max_mpoles > 0) Then
  Call mpoles_rotmat_set_halo(mpoles,domain,comm)
End If

!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
