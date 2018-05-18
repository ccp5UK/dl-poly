!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

! for new simulations when using the relaxed shell model
! set shells on top of their cores preventatively

     If ( (megshl > 0 .and. keyshl == 2) .and. &
          (keyres == 0 .and. nstep == 0 .and. nsteql > 0) ) Then
        Call core_shell_on_top(comm)

! Refresh mappings

        Call w_refresh_mappings()
     End If

100  Continue ! Only used when relaxed is false

! Initialise force arrays and possible torques for multipolar electrostatics
! and stress tensor (these are all additive in the force subroutines)

     fxx = 0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp
     If (mximpl > 0) Then
        mptrqx=0.0_wp ; mptrqy=0.0_wp ; mptrqz=0.0_wp
     End If
     stress = 0.0_wp

! Calculate pair-like forces (metal,vdw,electrostatic) and add lrc

     If (.not.(mxmet == 0 .and. keyfce == 0 .and. l_n_v .and. mxrdf == 0 .and. kimim == ' ')) &
        Call two_body_forces                      &
           (rcut,rlnk,rvdw,rmet,pdplnc,keyens,    &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           lrdf,nstrdf,leql,nsteql,nstep,         &
           elrc,virlrc,elrcm,vlrcm,               &
           engcpe,vircpe,engsrp,virsrp,stress,ewld,comm)

! Calculate tersoff forces

     If (ntpter > 0) Call tersoff_forces(rcter,engter,virter,stress,comm)

! Calculate three-body forces

     If (ntptbp > 0) Call three_body_forces(rctbp,engtbp,virtbp,stress,comm)

! Calculate four-body forces

     If (ntpfbp > 0) Call four_body_forces(rcfbp,engfbp,virfbp,stress,comm)

! Calculate shell model forces

     If (megshl > 0) Call core_shell_forces(engshl,virshl,stress,comm)

! Calculate tethered atom forces

     If (megtet > 0) Call tethers_forces(engtet,virtet,stress,comm)

! Calculate bond forces

     If (megbnd > 0) Then
        ltmp = (mxgbnd1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstbnd) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call bonds_forces(isw,engbnd,virbnd,stress,rcut,keyfce,alpha,epsq,engcpe,vircpe,comm)
     End If

! Calculate valence angle forces

     If (megang > 0) Then
        ltmp = (mxgang1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstang) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call angles_forces(isw,engang,virang,stress,comm)
     End If

! Calculate dihedral forces

     If (megdih > 0) Then
        ltmp = (mxgdih1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstdih) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call dihedrals_forces(isw,engdih,virdih,stress, &
           rcut,rvdw,keyfce,alpha,epsq,engcpe,vircpe,engsrp,virsrp,comm)
     End If

! Calculate inversion forces

     If (meginv > 0) Then
        ltmp = (mxginv1 > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstinv) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call inversions_forces(isw,enginv,virinv,stress,comm)
     End If

! Apply external field

     If (keyfld > 0) Call external_field_apply(keyshl,time,leql,nsteql,nstep,engfld,virfld,comm)

! Apply PLUMED driven dynamics

     If (l_plumed) Then
        stpcfg = engcpe + engsrp + engter + engtbp + engfbp + &
                 engshl + engtet + engfld +                   &
                 engbnd + engang + engdih + enginv
        Call plumed_apply(xxx,yyy,zzz,nstrun,nstep,stpcfg,stress,comm)
     End If
! Apply pseudo thermostat - force cycle (0)

     If (thermo%l_pseudo) Then
           Call pseudo_vv                               &
           (0,keyshl,keyens,tstep, &
           nstep,strkin,strknf,strknt,engke,engrot,thermo,comm)
     End If

! Cap forces in equilibration mode

     If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,thermo%temp,comm)

! Frozen atoms option

     Call freeze_atoms()

! Minimisation option and Relaxed shell model optimisation

     If (lsim .and. (lmin .or. keyshl == 2)) Then
        stpcfg = engcpe + engsrp + engter + engtbp + engfbp + &
                 engshl + engtet + engfld +                   &
                 engbnd + engang + engdih + enginv

        If (keyshl == 2) Call core_shell_relax(l_str,relaxed_shl,lrdf,rlx_tol,megshl,stpcfg,comm)

        If (.not.relaxed_shl) Go To 200 ! Shells relaxation takes priority over minimisation

        If (lmin .and. nstep >= 0 .and. nstep <= nstrun .and. nstep <= nsteql) Then
           If      (nstmin == 0 .and. nstep == 0) Then
              Call minimise_relax &
           (l_str .or. keyshl == 2,relaxed_min,lrdf,megatm,megcon,megpmf,megrgd, &
           keymin,min_tol,tstep,stpcfg,comm)
           Else If (nstmin >  0 .and. nstep >  0) Then
              If (Mod(nstep-nsteql,nstmin) == 0) Call minimise_relax &
           (l_str .or. keyshl == 2,relaxed_min,lrdf,megatm,megcon,megpmf,megrgd, &
           keymin,min_tol,tstep,stpcfg,comm)
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
           If (thermo%l_langevin) Then
              Call langevin_forces(nstep,thermo%temp,tstep,thermo%chi,fxl,fyl,fzl)
              If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxl,fyl,fzl,comm)
              Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl,comm)
           Else
              Call rigid_bodies_str_ss(strcom,comm)
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

     Call gsum(comm,stress)

! If RBs are present update forces on shared ones

     If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz,comm)


!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
