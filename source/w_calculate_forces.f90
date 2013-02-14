!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! for new simulations when using the relaxed shell model
! set shells on top of their cores preventatively

     If (megshl > 0 .and. keyshl == 2 .and. nstep == 0 .and. nsteql > 0) Then
        Call core_shell_on_top()

! Relocate atoms to new domains and restore bonding description

        Call relocate_particles      &
           (imcon,rcut,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

        Call set_halo_particles(imcon,rcut,keyfce,lbook)
     End If

100  Continue ! Only used when relaxed is false

! Initialise force arrays and stress tensor (these are all additive
! in the force subroutines)

     fxx = 0.0_wp
     fyy = 0.0_wp
     fzz = 0.0_wp

     stress = 0.0_wp

! Calculate pair-like forces (metal,vdw,electrostatic) and add lrc

     If (.not.(mxmet == 0 .and. keyfce == 0 .and. l_n_v .and. l_n_r)) &
        Call two_body_forces                      &
           (imcon,rcut,rvdw,rmet,keyens,          &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           lrdf,nstrdf,leql,nsteql,nstep,         &
           elrc,virlrc,elrcm,vlrcm,               &
           engcpe,vircpe,engsrp,virsrp,stress)

! Calculate tersoff forces

     If (ntpter > 0) Call tersoff_forces(imcon,rcter,engter,virter,stress)

! Calculate three-body forces

     If (ntptbp > 0) Call three_body_forces(imcon,rctbp,engtbp,virtbp,stress)

! Calculate four-body forces

     If (ntpfbp > 0) Call four_body_forces(imcon,rcfbp,engfbp,virfbp,stress)

! Calculate shell model forces

     If (megshl > 0) Call core_shell_forces(imcon,engshl,virshl,stress)

! Calculate tethered atom forces

     If (megtet > 0) Call tethers_forces(imcon,engtet,virtet,stress)

! Calculate bond forces

     If (megbnd > 0) Call bonds_forces(imcon,engbnd,virbnd,stress, &
                       rcut,keyfce,alpha,epsq,engcpe,vircpe)

! Calculate valence angle forces

     If (megang > 0) Call angles_forces(imcon,engang,virang,stress)

! Calculate dihedral forces

     If (megdih > 0) Call dihedrals_forces(imcon,engdih,virdih,stress, &
           rcut,rvdw,keyfce,alpha,epsq,engcpe,vircpe,engsrp,virsrp)

! Calculate inversion forces

     If (meginv > 0) Call inversions_forces(imcon,enginv,virinv,stress)

! Apply external field

     If (keyfld > 0) Call external_field_apply(imcon,keyshl,tstep,engfld,virfld)

! Apply pseudo thermostat - force cycle (0)

     If (lpse) Then
        If (l_vv) Then
           Call pseudo_vv                               &
           (0,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           strkin,strknf,strknt,engke,engrot)
        Else
           Call pseudo_lfv                              &
           (0,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           strkin,strknf,strknt,engke,engrot)
        End If
     End If

! Cap forces in equilibration mode

     If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,temp)

! Frozen atoms option

     Call freeze_atoms()

! Minimisation option and Relaxed shell model optimisation

     If (lmin .or. keyshl == 2) Then
        stpcfg = engcpe + engsrp + engter + engtbp + engfbp + &
                 engshl + engtet + engfld +                   &
                 engbnd + engang + engdih + enginv

        If (keyshl == 2) Call core_shell_relax(l_str,relaxed_shl,lrdf,rlx_tol,megshl,stpcfg)

        If (.not.relaxed_shl) Go To 200 ! Shells relaxation takes priority over minimisation

        If (lmin .and. nstep >= 0 .and. nstep <= nstrun .and. nstep <= nsteql) Then
           If      (nstmin == 0 .and. nstep == 0) Then
              Call minimise_relax &
           (l_str,relaxed_min,lrdf,imcon,megatm,megcon,megpmf,megrgd, &
           keymin,min_tol,tstep,stpcfg)
           Else If (nstmin >  0 .and. nstep >  0) Then
              If (Mod(nstep-nsteql,nstmin) == 0) Call minimise_relax &
           (l_str,relaxed_min,lrdf,imcon,megatm,megcon,megpmf,megrgd, &
           keymin,min_tol,tstep,stpcfg)
           End If
        End If

200     Continue

        If (.not.(relaxed_shl .and. relaxed_min)) Then

! Relocate atoms to new domains and restore bonding description

           Call relocate_particles   &
           (imcon,rcut,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

           Call set_halo_particles(imcon,rcut,keyfce,lbook)

           Go To 100

        End If
     End If

! Get RB COM stress and virial at restart only

     If (newjob) Then
        If (megrgd > 0) Then
           If (l_lan .and. (.not.l_lan_s)) Then ! Only meaningful for l_vv
              Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl)
           Else                                 ! Covers both VV and LFV
              Call rigid_bodies_str_ss(strcom)
           End If
           vircom=-(strcom(1)+strcom(5)+strcom(9))
        End If
     End If

! Total virial (excluding constraint, PMF and RB COM virials for npt routines)
! Total stress (excluding constraint, PMF, RB COM and kinetic stress for npt routines)
!
! NOTE(1):  virsrp already includes virlrc and vlrcm(0) and so
!           does the stress diagonal elements (by minus a third),
!           engsrp includes elrc and elrcm(0)
!
! NOTE(2):  virfbp, virinv and virdih are allegedly always zero

     virtot = vircpe + virsrp + virter + virtbp + virfbp + &
              virshl + virtet + virbnd + virang + virdih + virinv + virfld

     If (mxnode > 1) Call gsum(stress)

! If RBs are present update forces on shared ones

     If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz)


!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!