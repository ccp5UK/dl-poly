!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

! for new simulations when using the relaxed shell model
! set shells on top of their cores preventatively

     If ( (cshell%megshl > 0 .and. cshell%keyshl == SHELL_RELAXED) .and. &
          (keyres == 0 .and. nstep == 0 .and. nsteql > 0) ) Then
        Call core_shell_on_top(cshell,comm)

! Refresh mappings

        Call w_refresh_mappings(cshell,cons,pmf,stat,msd_data,bond,angle, &
          dihedral,inversion,tether,neigh,sites,mpole,rigid)
     End If

100  Continue ! Only used when relaxed is false

! Initialise force arrays and possible torques for multipolar electrostatics
! and stress tensor (these are all additive in the force subroutines)

     fxx = 0.0_wp ; fyy = 0.0_wp ; fzz = 0.0_wp
     If (mpole%max_mpoles > 0) Then
        mpole%torque_x=0.0_wp ; mpole%torque_y=0.0_wp ; mpole%torque_z=0.0_wp
     End If
     stat%stress = 0.0_wp

! Calculate pair-like forces (metal,vdws,electrostatic) and add lrc

     If (.not.(met%max_metal == 0 .and. electro%key == ELECTROSTATIC_NULL .and. &
       l_n_v .and. rdf%max_rdf == 0 .and. kimim == ' ')) Then
       Call two_body_forces(pdplnc,thermo%ensemble,nstfce,lbook,megfrz, &
         leql,nsteql,nstep,cshell,stat,ewld,devel,met,pois,neigh,sites,vdws,rdf, &
         mpole,electro,tmr,comm)
     End If

! Calculate tersoff forces

     If (tersoffs%n_potential > 0) Call tersoff_forces(tersoffs,stat,neigh,comm)

! Calculate three-body forces

     If (threebody%ntptbp > 0) Call three_body_forces(stat,threebody,neigh,comm)

! Calculate four-body forces

     If (fourbody%n_potential > 0) Call four_body_forces(fourbody,stat,neigh,comm)
     call start_timer(tmr%t_bonded)
! Calculate shell model forces

     If (cshell%megshl > 0) Call core_shell_forces(cshell,stat,comm)

! Calculate tethered atom forces

     If (tether%total > 0) Call tethers_forces(stat,tether,comm)

! Calculate bond forces

     If (bond%total > 0) Then
        ltmp = (bond%bin_pdf > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstbnd) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call bonds_forces(isw,stat%engbnd,stat%virbnd,stat%stress,neigh%cutoff, &
          stat%engcpe,stat%vircpe,bond,mpole,electro,comm)
     End If

! Calculate valence angle forces

     If (angle%total > 0) Then
        ltmp = (angle%bin_adf > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstang) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call angles_forces(isw,stat%engang,stat%virang,stat%stress,angle,comm)
     End If

! Calculate dihedral forces

     If (dihedral%total > 0) Then
        ltmp = (dihedral%bin_adf > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstdih) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call dihedrals_forces(isw,stat%engdih,stat%virdih,stat%stress, &
           neigh%cutoff,stat%engcpe,stat%vircpe,stat%engsrp, &
           stat%virsrp,dihedral,vdws,mpole,electro,comm)
     End If

! Calculate inversion forces

     If (inversion%total > 0) Then
        ltmp = (inversion%bin_adf > 0 .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstinv) == 0)

        isw = 1 + Merge(1,0,ltmp)
        Call inversions_forces(isw,stat%enginv,stat%virinv,stat%stress,inversion,comm)
     End If
     call stop_timer(tmr%t_bonded)

! Apply external field

     If (ext_field%key /= FIELD_NULL) Then
       Call external_field_apply(time,leql,nsteql,nstep,cshell,stat,rdf, &
         ext_field,rigid,comm)
     End If

! Apply PLUMED driven dynamics

     If (plume%l_plumed) Then
        stat%stpcfg =stat%engcpe + stat%engsrp + stat%engter + stat%engtbp + stat%engfbp + &
                 stat%engshl + stat%engtet + stat%engfld +                   &
                 stat%engbnd + stat%engang + stat%engdih + stat%enginv
        Call plumed_apply(xxx,yyy,zzz,nstrun,nstep,stat,plume,comm)
     End If
! Apply pseudo thermostat - force cycle (0)

     If (thermo%l_pseudo) Then
       Call pseudo_vv(0,tstep,nstep,sites%dof_site,cshell,stat,thermo,rigid,comm)
     End If

! Cap forces in equilibration mode

     If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,thermo%temp,comm)

! Frozen atoms option

     Call freeze_atoms()

! Minimisation option and Relaxed shell model optimisation

     If (lsim .and. (minimise%minimise .or. cshell%keyshl == SHELL_RELAXED)) Then
        stat%stpcfg = stat%engcpe + stat%engsrp + stat%engter + stat%engtbp + stat%engfbp + &
                 stat%engshl + stat%engtet + stat%engfld +                   &
                 stat%engbnd + stat%engang + stat%engdih + stat%enginv

        If (cshell%keyshl == SHELL_RELAXED) Then
          Call core_shell_relax(l_str,relaxed_shl,rdf%l_collect,rlx_tol,stat%stpcfg,cshell,stat,comm)
        End If

        If (.not.relaxed_shl) Go To 200 ! Shells relaxation takes priority over minimisation

        If (minimise%minimise .and. nstep >= 0 .and. nstep <= nstrun .and. nstep <= nsteql) Then
          If      (minimise%freq == 0 .and. nstep == 0) Then
            Call minimise_relax(l_str .or. cshell%keyshl == SHELL_RELAXED, &
              rdf%l_collect,megatm,pmf%megpmf,tstep,stat%stpcfg,stat,pmf,cons, &
              netcdf,minimise,rigid,comm)
          Else If (minimise%freq >  0 .and. nstep >  0) Then
            If (Mod(nstep-nsteql,minimise%freq) == 0) Then
              Call minimise_relax(l_str .or. cshell%keyshl == SHELL_RELAXED, &
                rdf%l_collect,megatm,pmf%megpmf,tstep,stat%stpcfg,stat,pmf,cons, &
                netcdf,minimise,rigid,comm)
            End If
          End If
        End If

200     Continue

! Refresh mappings

        If (.not.(relaxed_shl .and. minimise%relaxed)) Then
           Call w_refresh_mappings(cshell,cons,pmf,stat,msd_data,bond,angle, &
             dihedral,inversion,tether,neigh,sites,mpole,rigid)
           Go To 100
        End If
     End If

! Get RB COM stress and virial at restart only - also available at w_at_start_vv for levcfg==2

     If (newjob) Then
        If (rigid%total > 0) Then
           If (thermo%l_langevin) Then
              Call langevin_forces(nstep,thermo%temp,tstep,thermo%chi,fxl,fyl,fzl,cshell)
              If (rigid%share) Then
                Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,fxl,fyl,fzl,comm)
              End If
              Call rigid_bodies_str__s(stat%strcom,fxx+fxl,fyy+fyl,fzz+fzl,rigid,comm)
           Else
              Call rigid_bodies_str_ss(stat%strcom,rigid,comm)
           End If
           stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
        End If
     End If

! Total virial (excluding constraint, PMF and RB COM virials for npt routines)
! Total stress (excluding constraint, PMF, RB COM and kinetic stress for npt routines)
!
! NOTE(1):  virsrp already includes vdws%vlrc and vlrcm(0) and so
!           does the stress diagonal elements (by minus a third),
!           engsrp includes vdws%elrc and elrcm(0)
!
! NOTE(2):  virfbp, virinv and virdih are allegedly always zero

     stat%virtot = stat%vircpe + stat%virsrp + stat%virter + stat%virtbp + stat%virfbp + &
              stat%virshl + stat%virtet + stat%virbnd + stat%virang + stat%virdih + &
              stat%virinv + stat%virfld

     Call gsum(comm,stat%stress)

! If RBs are present update forces on shared ones

     If (rigid%share) Then
       Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,fxx,fyy,fzz,comm)
     End If


!!!!!!!!!!!!!!!!!!  W_CALCULATE_FORCES INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
