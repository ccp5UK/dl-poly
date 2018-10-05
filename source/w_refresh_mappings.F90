!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Scale t=0 reference positions

If (flow%step > 0) Call xscale(cnfig,thermo%tstep,thermo,stat,neigh,rigid,domain,comm)

! Check VNL conditioning

Call vnl_check(flow%strict,cnfig%width,neigh,stat,domain,cnfig,ewld%bspline,kim_data,comm)

If (neigh%update) Then

  ! Relocate atoms to new domains and restore bonding description

  Call relocate_particles(cnfig%dvar,neigh%cutoff_extended,flow%book, &
    msd_data%l_msd,cnfig%megatm,flow,cshell,cons,pmf,stat,ewld,thermo,green, &
    bond,angle,dihedral,inversion,tether,neigh,sites,minim,mpoles, &
    rigid,domain,cnfig,comm)

  ! Exchange atomic data in border regions

  Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,cnfig,ewld,kim_data,comm) ! inducing in here only

  ! Re-tag RBs when called again after the very first flow%time
  ! when it's done in rigid_bodies_setup <- build_book_intra

  If (rigid%on) Then
    Call rigid_bodies_tags(cnfig,rigid,comm)
    Call rigid_bodies_coms(cnfig,rigid%xxx,rigid%yyy,rigid%zzz,rigid,comm)
  End If

Else

  ! Exchange atomic positions in border regions

  Call refresh_halo_positions(domain,cnfig,kim_data,comm)
End If

! set and halo rotational matrices and their infinitesimal rotations

If (mpoles%max_mpoles > 0) Then
  Call mpoles_rotmat_set_halo(mpoles,domain,cnfig,comm)
End If

!!!!!!!!!!!!!!!!!!  W_REFRESH_MAPPINGS INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
