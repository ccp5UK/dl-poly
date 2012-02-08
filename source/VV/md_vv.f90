!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MD_VV INCLUSION MODULE  !!!!!!!!!!!!!!!!


! Calculate physical quantities at restart
! Calculate kinetic tensor and energy

  If (megrgd > 0) Then
     Call kinstresf(vxx,vyy,vzz,strknf)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

     strkin=strknf+strknt
  Else
     Call kinstress(vxx,vyy,vzz,strkin)
  End If
  engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START OF MOLECULAR DYNAMICS CALCULATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Do While ( (nstep < nstrun .or. (nstep == nstrun .and. newjob)) .and. &
             (timjob-timelp) > timcls )

! Apply impact

     If (nstep == tmd) Then
        If (idnode == 0) Write(nrite,"(/,          &
           & /,1x,'initiating IMPACT:',            &
           & /,1x,46('-'),                         &
           & /,1x,'particle (index)',15x,i10,      &
           & /,1x,'timestep (steps)',15x,i10,      &
           & /,1x,'energy   (keV)  ',18x,1p,e12.4, &
           & /,1x,'v-r(x,y,z)',1p,3e12.4,          &
           & /,1x,46('-'))") imd,tmd,emd,vmx,vmy,vmz

        If (nstep+1 <= nsteql) Call warning(380,Real(nsteql,wp),0.0_wp,0.0_wp)

        Call impact(imd,emd,vmx,vmy,vmz,megrgd)

! Correct kinetic stress and energy

        If (megrgd > 0) Then
           Call kinstresf(vxx,vyy,vzz,strknf)
           Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

           strkin=strknf+strknt

           engrot=getknr(rgdoxx,rgdoyy,rgdozz)
        Else
           Call kinstress(vxx,vyy,vzz,strkin)
        End If
        engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
     End If

! VV - repeat after FORCES (endwhile)!!!

     If (newjob .and. levcfg == 2) Then
        newjob = .false.

! If RBs are present update forces on shared ones

        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz)

! Get RB COM stress and virial at restart

        If (megrgd > 0) Then
           If (l_lan .and. (.not.l_lan_s)) Then
              Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl)
           Else
              Call rigid_bodies_str_ss(strcom)
           End If
           vircom=-(strcom(1)+strcom(5)+strcom(9))
        End If

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart

        If (ltraj .and. keyres /= 1) Call trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)
        If (ldef .and. keyres /= 1) Call defects_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)
        If (l_msd .and. keyres /= 1) Call msd_write &
           (keyres,nstmsd,istmsd,megatm,nstep,tstep,time)
        If (lrsd .and. keyres /= 1) Call rsd_write &
           (imcon,keyres,nsrsd,isrsd,rrsd,nstep,tstep,time)
     End If

! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

     If (nstep >= 0 .and. nstep < nstrun .and. levcfg == 2) Then

! Increase step counter

        nstep=nstep+1

! zero Kelvin structure optimisation

        If (lzero .and. nstep <= nsteql) Call zero_k_optimise(strkin,strknf,strknt,engke,engrot)

! Integrate equations of motion - velocity verlet first stage

        isw = 0

        If (m_rgd == 0) Then
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_0_vv                                &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           chit)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l0_vv                                        &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,chi,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,                &
           megpmf,strpmf,virpmf)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a0_vv                         &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,   &
           keyshl,taut,soft,strkin,engke,            &
           imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit,cint,consv)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g0_vv                                    &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,            &
           megpmf,strpmf,virpmf,                                &
           sigma,taut,gama,chit,cint,consv)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,chi,consv,                         &
           press,tai,chip,eta,virtot,                      &
           elrc,virlrc)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit,                                &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,taut,chit,cint,consv,              &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,taut,chit,cint,consv,              &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,chi,consv,                     &
           press,strext,ten,tai,chip,eta,stress,           &
           elrc,virlrc)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,sigma,taut,chit,                            &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,taut,chit,cint,consv,          &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,taut,chit,cint,consv,          &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        Else
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_1_vv                   &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           chit,                              &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           chi,                                    &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           keyshl,taut,soft,                       &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,consv,        &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           sigma,taut,gama,chit,cint,consv,        &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,chi,consv,            &
           degrot,press,tai,chip,eta,         &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           press,taup,chip,eta,               &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,        &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,        &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           iso,degfre,sigma,chi,consv,        &
           degrot,press,tai,chip,eta,         &
           stress,strext,ten,elrc,virlrc,     &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           iso,sigma,taut,chit,               &
           press,taup,chip,eta,               &
           stress,strext,ten,elrc,virlrc,     &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_vv                      &
           (isw,lvar,mndis,mxdis,mxstp,tstep,     &
           iso,degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,            &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,                   &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_vv                      &
           (isw,lvar,mndis,mxdis,mxstp,tstep,     &
           iso,degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,            &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,                   &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        End If

! Scale t=0 reference positions

        Call xscale(imcon,m_rgd,keyens,tstep,eta)

! Relocate atoms to new domains and restore bonding description

        Call relocate_particles      &
           (imcon,rcut,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

! Exchange atomic data in border regions

        Call set_halo_particles(imcon,rcut,keyfce,lbook)

     End If ! DO THAT ONLY IF 0<=nstep<nstrun AND FORCES ARE PRESENT (levcfg=2)

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

     If (lpse) Call pseudo_vv                           &
           (0,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           strkin,strknf,strknt,engke,engrot)

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

! Calculate physical quantities and print at the very start

     If (nstep == 0) Then

! Calculate total stress tensor

        strtot = strcon + strpmf + stress + strkin + strcom

! Get core-shell kinetic energy for adiabatic shell model

        If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke)

        Call statistics_collect              &
           (leql,nsteql,lzdn,nstzdn,         &
           keyres,keyens,iso,intsta,imcon,   &
           degfre,degshl,degrot,             &
           nstep,tstep,time,tmst,            &
           engcpe,vircpe,engsrp,virsrp,      &
           engter,virter,                    &
           engtbp,virtbp,engfbp,virfbp,      &
           engshl,virshl,shlke,              &
           vircon,virpmf,                    &
           engtet,virtet,engfld,virfld,      &
           engbnd,virbnd,engang,virang,      &
           engdih,virdih,enginv,virinv,      &
           engke,engrot,consv,vircom,        &
           strtot,press,strext,              &
           stpeng,stpvir,stpcfg,stpeth,      &
           stptmp,stpprs,stpvol)

        Call gtime(timelp)

        If (idnode == 0) Then
           Write(nrite,"(/,/,1x,130('-'),/,/,                              &
                 & 5x,'step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',      &
                 & 5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',    &
                 & 5x,'eng_dih',5x,'eng_tet',/,1x,'time(ps)',5x,' eng_pv', &
                 & 4x,'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',   &
                 & 5x,'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,  &
                 & 1x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',   &
                 & 5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',           &
                 & 5x,'vir_pmf',7x,'press',/,/,1x,130('-'))")

           Write(nrite,"(1x,i8,1p,9e12.4,/,0p,f9.5,1p,9e12.4,              &
                   & /,1x,0p,f8.1,1p,9e12.4)") nstep,(stpval(i),i=1,9),    &
                   time,(stpval(i),i=10,18),timelp,(stpval(i),i=19,27)

           Write(nrite,"(1x,130('-'))")
        End If

     End If

! DO THAT ONLY IF 0<nstep<=nstrun AND THIS IS AN OLD JOB (newjob=.false.)

     If (nstep > 0 .and. nstep <= nstrun .and. (.not.newjob)) Then

! Integrate equations of motion - velocity verlet second pass

        isw = 1

        If (m_rgd == 0) Then
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_0_vv                                &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           chit)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l0_vv                                        &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,chi,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,                &
           megpmf,strpmf,virpmf)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a0_vv                         &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,   &
           keyshl,taut,soft,strkin,engke,            &
           imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit,cint,consv)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g0_vv                                    &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,            &
           megpmf,strpmf,virpmf,                                &
           sigma,taut,gama,chit,cint,consv)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,chi,consv,                         &
           press,tai,chip,eta,virtot,                      &
           elrc,virlrc)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           sigma,taut,chit,                                &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,taut,chit,cint,consv,              &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           degfre,sigma,taut,chit,cint,consv,              &
           press,taup,chip,eta,virtot,                     &
           elrc,virlrc)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,chi,consv,                     &
           press,strext,ten,tai,chip,eta,stress,           &
           elrc,virlrc)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,sigma,taut,chit,                            &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,taut,chit,cint,consv,          &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m0_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,       &
           megpmf,strpmf,virpmf,                           &
           iso,degfre,sigma,taut,chit,cint,consv,          &
           press,strext,ten,taup,chip,eta,stress,          &
           elrc,virlrc)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        Else
           If      (keyens ==  0) Then

! Microcanonical ensemble

              Call nve_1_vv                   &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens ==  1) Then

! Evans thermostat (Gaussian temperature constraints)

              Call nvt_e1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           chit,                              &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 10) Then

! Langevin thermostat (Stochastic Dynamics)

              Call nvt_l1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           chi,                                    &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 11) Then

! Andersen thermostat (Stochastic Dynamics)

              Call nvt_a1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           keyshl,taut,soft,                       &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 12) Then

! Berendsen thermostat

              Call nvt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 13) Then

! Nose-Hoover thermostat

              Call nvt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,consv,        &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 14) Then

! Gentle-Stochastic thermostat

              Call nvt_g1_vv                       &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           sigma,taut,gama,chit,cint,consv,        &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

           Else If (keyens == 20) Then

! Langevin thermostat and isotropic barostat

              Call npt_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,chi,consv,            &
           degrot,press,tai,chip,eta,         &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 21) Then

! Berendsen thermostat and isotropic barostat

              Call npt_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           press,taup,chip,eta,               &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 22) Then

! Nose-Hoover thermostat and isotropic barostat

              Call npt_h1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,        &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 23) Then

! Martyna-Tuckerman-Klein (MTK) thermostat and isotropic barostat

              Call npt_m1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,        &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 30) Then

! Langevin thermostat and barostat anisotropic (cell shape varying)

              Call nst_l1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           iso,degfre,sigma,chi,consv,        &
           degrot,press,tai,chip,eta,         &
           stress,strext,ten,elrc,virlrc,     &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 31) Then

! Berendsen thermostat and barostat anisotropic (cell shape varying)

              Call nst_b1_vv                  &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           iso,sigma,taut,chit,               &
           press,taup,chip,eta,               &
           stress,strext,ten,elrc,virlrc,     &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

           Else If (keyens == 32) Then

! Nose-Hoover thermostat and anisotropic barostat (cell shape varying)

              Call nst_h1_vv                      &
           (isw,lvar,mndis,mxdis,mxstp,tstep,     &
           iso,degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,            &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,                   &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

           Else If (keyens == 33) Then

! MTK thermostat and anisotropic barostat (cell shape varying)

              Call nst_m1_vv                      &
           (isw,lvar,mndis,mxdis,mxstp,tstep,     &
           iso,degfre,sigma,taut,chit,cint,consv, &
           degrot,press,taup,chip,eta,            &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,                   &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

           Else

! Invalid ensemble option

              Call error(430)

           End If
        End If

! Apply external field

        If (keyfld > 0) Call external_field_correct(imcon)

! Apply pseudo thermostat - velocity cycle (1)

        If (lpse) Call pseudo_vv                        &
           (1,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           strkin,strknf,strknt,engke,engrot)

! Apply temperature regaussing

        If (ltgaus .and. nstep <= nsteql .and. Mod(nstep-nsteql,nstgaus) == 0) Then
           chit = 0.0_wp
           chip = 0.0_wp
           eta  = 0.0_wp

           Call regauss_temperature(imcon,megrgd)

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              stptmp = 2.0_wp*(engke+engrot) / (boltz*Real(degfre,wp))
              Do
                 Call scale_temperature(imcon,engke+engrot,degtra,degrot,degfre)
                 Call core_shell_quench(safe,stptmp)
                 If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
                 If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)
                 If (megrgd > 0) Call rigid_bodies_quench(imcon)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(imcon,engke+engrot,degtra,degrot,degfre)
           End If

! Correct kinetic stress and energy

           If (megrgd > 0) Then
              Call kinstresf(vxx,vyy,vzz,strknf)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz)
           Else
              Call kinstress(vxx,vyy,vzz,strkin)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If

! Apply temperature scaling

        If (ltscal .and. nstep <= nsteql .and. Mod(nstep-nsteql,nstscal) == 0) Then
           chit = 0.0_wp
           chip = 0.0_wp
           eta  = 0.0_wp

! quench constraints & PMFs

           If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)

! quench core-shell units in adiabatic model

           If (megshl > 0 .and. keyshl == 1) Then
              Do
                 Call scale_temperature(imcon,sigma,degtra,degrot,degfre)
                 Call core_shell_quench(safe,stptmp)
                 If (megcon > 0) Call constraints_quench(imcon,mxshak,tolnce)
                 If (megpmf > 0) Call pmf_quench(imcon,mxshak,tolnce)
                 If (megrgd > 0) Call rigid_bodies_quench(imcon)
                 If (safe) Exit
              End Do
           Else
              Call scale_temperature(imcon,sigma,degtra,degrot,degfre)
           End If

! Correct kinetic stress and energy

           If (megrgd > 0) Then
              Call kinstresf(vxx,vyy,vzz,strknf)
              Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

              strkin=strknf+strknt

              engrot=getknr(rgdoxx,rgdoyy,rgdozz)
           Else
              Call kinstress(vxx,vyy,vzz,strkin)
           End If
           engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))
        End If

! Get complete stress tensor

        strtot = strcon + strpmf + stress + strkin + strcom

! Get core-shell kinetic energy for adiabatic shell model

        If (megshl > 0 .and. keyshl == 1) Call core_shell_kinetic(shlke)

! Update total time of simulation

        time = time + tstep

! Calculate physical quantities and collect statistics

        Call statistics_collect              &
           (leql,nsteql,lzdn,nstzdn,         &
           keyres,keyens,iso,intsta,imcon,   &
           degfre,degshl,degrot,             &
           nstep,tstep,time,tmst,            &
           engcpe,vircpe,engsrp,virsrp,      &
           engter,virter,                    &
           engtbp,virtbp,engfbp,virfbp,      &
           engshl,virshl,shlke,              &
           vircon,virpmf,                    &
           engtet,virtet,engfld,virfld,      &
           engbnd,virbnd,engang,virang,      &
           engdih,virdih,enginv,virinv,      &
           engke,engrot,consv,vircom,        &
           strtot,press,strext,              &
           stpeng,stpvir,stpcfg,stpeth,      &
           stptmp,stpprs,stpvol)

! line-printer output every nstbpo steps

        If (lines == 0 .or. Mod(nstep,nstbpo) == 0) Then

! Update cpu time

           Call gtime(timelp)

           If (idnode == 0) Then
              If (Mod(lines,npage) == 0) Write(nrite,"(1x,130('-'),/,/,    &
                 & 5x,'step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',      &
                 & 5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',    &
                 & 5x,'eng_dih',5x,'eng_tet',/,1x,'time(ps)',5x,' eng_pv', &
                 & 4x,'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',   &
                 & 5x,'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,  &
                 & 1x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',   &
                 & 5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',           &
                 & 5x,'vir_pmf',7x,'press',/,/,1x,130('-'))")

              Write(nrite,"(1x,i8,1p,9e12.4,/,0p,f9.5,1p,9e12.4,           &
                   & /,1x,0p,f8.1,1p,9e12.4)") nstep,(stpval(i),i=1,9),    &
                   time,(stpval(i),i=10,18),timelp,(stpval(i),i=19,27)

              Write(nrite,"(/,2x,'rolling',1p,9e12.4,/,1x,'averages',      &
                   & 1p,9e12.4,/,9x,9e12.4)") (ravval(i),i=1,27)

              Write(nrite,"(1x,130('-'))")
           End If

           lines=lines+1
        End If

! Report end of equilibration period

        If (ltscal .and. nstep == nsteql) Then
           ltscal=.false.
           If (idnode == 0) Write(nrite,"(/,/,1x,                         &
              & 'switching off temperature scaling at step ',i6,/,/,/,1x, &
              & 130('-'))") nstep
        End If

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT

        If (ltraj) Call trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)
        If (ldef) Call defects_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)
        If (l_msd) Call msd_write &
           (keyres,nstmsd,istmsd,megatm,nstep,tstep,time)
        If (lrsd) Call rsd_write &
           (imcon,keyres,nsrsd,isrsd,rrsd,nstep,tstep,time)

! Save restart data in event of system crash

        If (Mod(nstep,ndump) == 0 .and. nstep /= nstrun .and. (.not.l_tor)) &
           Call system_revive                                       &
           (imcon,rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

     End If ! DO THAT ONLY IF 0<nstep<=nstrun AND THIS IS AN OLD JOB (newjob=.false.)

     If (newjob) Then
        newjob = .false.

! Get RB COM stress and virial at restart

        If (megrgd > 0) Then
           If (l_lan .and. (.not.l_lan_s)) Then
              Call rigid_bodies_str__s(strcom,fxx+fxl,fyy+fyl,fzz+fzl)
           Else
              Call rigid_bodies_str_ss(strcom)
           End If
           vircom=-(strcom(1)+strcom(5)+strcom(9))
        End If

! Write HISTORY, DEFECTS, MSDTMP & DISPDAT if needed immediately after restart

        If (ltraj .and. keyres /= 1) Call trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)
        If (ldef .and. keyres /= 1) Call defects_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)
        If (l_msd .and. keyres /= 1) Call msd_write &
           (keyres,nstmsd,istmsd,megatm,nstep,tstep,time)
        If (lrsd .and. keyres /= 1) Call rsd_write &
           (imcon,keyres,nsrsd,isrsd,rrsd,nstep,tstep,time)
     End If

! Complete time check

     Call gtime(timelp)

! Change levcfg after restart if forces and stress are
! (re)calculated and print start-up time

     If (levcfg == 1) Then
        levcfg=2
        If (idnode == 0) Write(nrite,'(/,/,/,1x, &
           & "time elapsed since job start: ", f12.3, " sec",/)') timelp
     End If

! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

     i=20
     If (nstep > 0) Then
        If ( Mod(nstep,i*nstbpo) == 0 .or.                        &
             (timelp > Real(i*60,wp) .and.                        &
              timelp-Real( ((Int(timelp)/(i*60)) * i*60) , wp ) < &
              timelp/Real( nstep , wp) ) ) Then

           If (idnode == 0) Then
              Inquire(File='OUTPUT', Exist=l_out, Position=c_out)
              Call strip_blanks(c_out)
              Call lower_case(c_out)
              If (l_out .and. c_out(1:6) == 'append') Then
                 Close(Unit=nrite)
                 Open(Unit=nrite, File='OUTPUT', Position='append')
              End If
           End If

        End If
     End If

  End Do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END OF MOLECULAR DYNAMICS CALCULATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MD_VV INCLUSION MODULE  !!!!!!!!!!!!!!!!
