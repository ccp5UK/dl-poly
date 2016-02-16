!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! for new simulations when using the relaxed shell model
! set shells on top of their cores preventatively

     If (megshl > 0 .and. keyshl == 2 .and. nstep == 0 .and. nsteql > 0) Then

        Call core_shell_on_top()

! Check VNL conditioning

        Call vnl_check(l_str,imcon,m_rgd,rcut,rpad,rlnk,width)

        If (l_vnl) Then

! Relocate atoms to new domains and restore bonding description

           Call relocate_particles   &
           (imcon,rlnk,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

           Call set_halo_particles(imcon,rlnk,keyfce)

! Re-tag RBs when called again after the very first time
! when it's done in rigid_bodies_setup <- build_book_intra

           If (m_rgd > 0) Then
              Call rigid_bodies_tags()
              Call rigid_bodies_coms(imcon,xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz)
           End If

        Else

! Exchange atomic positions in border regions

           Call refresh_halo_positions()

        End If

     End If

100  Continue ! Only used when relaxed is false

! Initialise force arrays and stress tensor (these are all additive
! in the force subroutines)

     fxx = 0.0_wp
     fyy = 0.0_wp
     fzz = 0.0_wp

     stress = 0.0_wp

! Calculate pair-like forces (metal,vdw,electrostatic) and add lrc

     If (.not.(mxmet == 0 .and. keyfce == 0 .and. l_n_v .and. mxrdf == 0 .and. kim == ' ')) &
        Call two_body_forces                      &
           (imcon,rcut,rlnk,rvdw,rmet,keyens,     &
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

     If (megbnd > 0) Then
        ltmp = (mxgbnd1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstbnd) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call bonds_forces(isw,imcon,engbnd,virbnd,stress, &
                       rcut,keyfce,alpha,epsq,engcpe,vircpe)
     End If

! Calculate valence angle forces

     If (megang > 0) Then
        ltmp = (mxgang1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstang) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call angles_forces(isw,imcon,engang,virang,stress)
     End If

! Calculate dihedral forces

     If (megdih > 0) Then
        ltmp = (mxgdih1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstdih) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call dihedrals_forces(isw,imcon,engdih,virdih,stress, &
           rcut,rvdw,keyfce,alpha,epsq,engcpe,vircpe,engsrp,virsrp)
     End If

! Calculate inversion forces

     If (meginv > 0) Then
        ltmp = (mxginv1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstinv) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call inversions_forces(isw,imcon,enginv,virinv,stress)
     End If

! Apply external field

     If (keyfld > 0) Call external_field_apply(imcon,keyshl,time,engfld,virfld)

     If (l_plumed) Call plumed_apply(xxx,yyy,zzz,nstrun,nstep,stpcfg,stress)

! Apply pseudo thermostat - force cycle (0)

     If (lpse) Then
        If (l_vv) Then
           Call pseudo_vv                               &
           (0,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           nstep,strkin,strknf,strknt,engke,engrot)
        Else
           Call pseudo_lfv                              &
           (0,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           nstep,strkin,strknf,strknt,engke,engrot)
        End If
     End If

! Cap forces in equilibration mode

     If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,temp)

! Frozen atoms option

     Call freeze_atoms()

! Minimisation option and Relaxed shell model optimisation

     If (lsim .and. (lmin .or. keyshl == 2)) Then
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

! Refresh mappings

        If (.not.(relaxed_shl .and. relaxed_min)) Then
           Call w_refresh_mappings()

           Go To 100
        End If
     End If

! Get RB COM stress and virial at restart only - also available at w_at_start_vv for levcfg==2

     If (newjob) Then
        If (megrgd > 0) Then
           If (l_lan) Then
              Call langevin_forces(nstep,temp,tstep,chi,fxl,fyl,fzl)
              If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxl,fyl,fzl)
              Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl)
           Else
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
