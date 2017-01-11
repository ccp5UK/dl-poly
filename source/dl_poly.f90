Program dl_poly



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 is an stfc/ccp5 program package for the dynamical
! simulation of molecular systems.
!
! dl_poly_4 is the property of the stfc daresbury laboratory,
! daresbury, warrington wa4 4ad.  no part of the package may
! be redistributed to third parties without the consent of
! daresbury laboratory.
!
! dl_poly_4 is available free of charge to academic institutions
! engaged in non-commercial research only.  potential users not
! in this category must consult the ccp5 program librarian at
! daresbury to negotiate terms of use.
!
! neither the stfc, daresbury laboratory, ccp5 nor the authors
! of this package claim that it is free from errors and do not
! accept liability for any loss or damage that may arise from
! its use.  it is the users responsibility to verify that the
! package dl_poly_4 is fit for the purpose the user intends for
! it.
!
! users of this package are recommended to consult the dl_poly_4
! user and reference manuals for the full terms and conditions
! of its use.
!
! dl_poly_4 is based on dl_poly_3 by i.t.todorov & w.smith.
!
! copyright - daresbury laboratory
! authors   - i.t.todorov & w.smith march 2016
! contrib   - i.j.bush, h.a.boateng, a.m.elena, a.b.g.chalk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! SETUP MODULES

  Use kinds_f90
  Use comms_module
  Use setup_module

! IO MODULE

  Use io_module

! SITE & CONFIG MODULES

  Use site_module
  Use config_module

! VNL module

  Use vnl_module

! DPD module

  Use dpd_module

! INTERACTION MODULES

  Use core_shell_module

  Use constraints_module
  Use pmf_module

  Use rigid_bodies_module

  Use tethers_module

  Use bonds_module
  Use angles_module
  Use dihedrals_module
  Use inversions_module

  Use mpoles_module

  Use vdw_module
  Use metal_module
  Use tersoff_module
  Use three_body_module
  Use four_body_module

  Use kim_module
  Use plumed_module

  Use external_field_module

! STATISTICS MODULES

  Use rdf_module
  Use z_density_module
  Use statistics_module
  Use greenkubo_module

! PARSE MODULE

  Use parse_module

! MSD MODULE

  Use msd_module

! KINETIC MODULE

  Use kinetic_module

! LANGEVIN MODULE

  Use langevin_module

! DEVELOPMENT MODULE

  Use development_module

! MAIN PROGRAM VARIABLES

  Implicit None

! newjob used for trajectory_write &
!                 defects_write    &
!                 msd_write        &
!                 rsd_write        &

  Logical, Save :: newjob = .true.

! OUTPUT existence

  Logical               :: l_out
  Character( Len = 10 ) :: c_out

! lines and page used for printing controls

  Integer       :: lines = 0 , &
                   npage = 8

! general flags

  Logical           :: ltmp,l_vv,l_n_e,l_n_v,       &
                       l_ind,l_str,l_top,           &
                       l_exp,lecx,lfcap,lzero,      &
                       lmin,ltgaus,ltscal,          &
                       lvar,leql,lpse,lsim,lfce,    &
                       lpana,lrdf,lprdf,lzdn,lpzdn, &
                       lvafav,lpvaf,                &
                       ltraj,ldef,lrsd,             &
                       safe,lbook,lexcl,            &
                       relaxed_shl = .true.,        &
                       relaxed_min = .true.

  Integer           :: i,j,isw,levcfg,nstfce,              &
                       nx,ny,nz,imd,tmd,                   &
                       keyres,nstrun,nsteql,               &
                       keymin,nstmin,nstgaus,nstscal,      &
                       keyens,iso,intsta,keypse,nstbpo,    &
                       keyfce,mxshak,mxquat,               &
                       nstbnd,nstang,nstdih,nstinv,        &
                       nstrdf,nstzdn,                      &
                       nstmsd,istmsd,nstraj,istraj,keytrj, &
                       nsdef,isdef,nsrsd,isrsd,            &
                       ndump,nstep,keyshl,                 &
                       atmfre,atmfrz,megatm,megfrz,        &
                       megshl,megcon,megpmf,megrgd,        &
                       megtet,megbnd,megang,megdih,meginv

! Degrees of freedom must be in long integers so we do 2.1x10^9 particles

  Integer(Kind=ip)  :: degfre,degshl,degtra,degrot

! elrc,virlrc - vdw energy and virial are scalars and in vdw_module
! elrcm,vlrcm - metal energy and virial are array-like and in metal_module

  Real( Kind = wp ) :: tsths,                                     &
                       timelp,timjob,timcls,tstep,time,tmst,      &
                       dvar,rcut,rpad,rlnk,                       &
                       rvdw,rmet,rbin,rcter,rctbp,rcfbp,          &
                       alpha,epsq,fmax,                           &
                       width,mndis,mxdis,mxstp,wthpse,tmppse,     &
                       rlx_tol,min_tol,tolnce,quattol,rdef,rrsd,  &
                       pdplnc,emd,vmx,vmy,vmz,temp,sigma,         &
                       press,strext(1:9),ten,                     &
                       taut,chi,soft,gama,taup,tai,               &
                       chit,eta(1:9),chip,cint,consv,             &
                       strtot(1:9),virtot,                        &
                       strkin(1:9),engke,strknf(1:9),strknt(1:9), &
                       engrot,strcom(1:9),vircom,                 &
                       engcpe,vircpe,engsrp,virsrp,               &
                       engter,virter,engtbp,virtbp,engfbp,virfbp, &
                       engshl,shlke,virshl,                       &
                       strcon(1:9),vircon,strpmf(1:9),virpmf,     &
                       stress(1:9),engtet,virtet,                 &
                       engbnd,virbnd,engang,virang,               &
                       engdih,virdih,enginv,virinv,               &
                       engfld,virfld,                             &
                       stptmp,stpprs,stpvol,stpcfg,stpeng,stpeth,stpvir

! SET UP COMMUNICATIONS & CLOCKING

  Call init_comms()
  If (mxnode > 1) Call gsync()
  Call gtime(timelp)

  Call scan_development()

! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES

  If (idnode == 0) Then
     If (.not.l_scr) Open(Unit=nrite, File='OUTPUT', Status='replace')

     Write(nrite,'(5(1x,a,/),(1x,a25,a8,a4,a14,a15/),1x,a,i10,a,/,5(1x,a,/))')  &
          "******************************************************************", &
          "*************  stfc/ccp5  program  library  package  ** D ********", &
          "*************  daresbury laboratory general purpose  *** L *******", &
          "**         **  classical molecular dynamics program  **** \ ******", &
          "** DL_POLY **  authors:   i.t.todorov   &   w.smith  ***** P *****", &
          "**         **  version:  ", DLP_VERSION,                 " /  "    , &
                                       DLP_RELEASE,          "  ****** O ****", &
          "*************  execution on  ",mxnode," process(es)  ******* L ***", &
          "*************  contributors' list:                   ******** Y **", &
          "*************  ------------------------------------  *************", &
          "*************  i.j.bush, h.a.boateng, r.davidchak,   *************", &
          "*************  m.a.seaton, a.v.brukhno, a.m.elena    *************", &
          "******************************************************************"

     Call build_info()

     Write(nrite,'(7(1x,a,/))') &
          "******************************************************************", &
          "****  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ****", &
          "****  Please do cite `J. Mater. Chem.', 16, 1911-1918 (2006)  ****", &
          "****  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ****", &
          "****  when publishing research data obtained using DL_POLY_4  ****", &
          "****  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ****", &
          "******************************************************************"
  End If

! TEST I/O

  Call scan_control_io()

! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
! (setup_module and domains_module)

  Call set_bounds                                     &
           (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,l_ind, &
           dvar,rcut,rpad,rlnk,rvdw,rmet,rbin,nstfce,alpha,width)

  Call gtime(timelp)
  If (idnode == 0) Then
     Write(nrite,'(/,1x,a)') "*** pre-scanning stage (set_bounds) DONE ***"
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
  End If

! ALLOCATE SITE & CONFIG

  Call allocate_site_arrays()
  Call allocate_config_arrays()

! ALLOCATE DPD ARRAYS

  Call allocate_dpd_arrays()

! ALLOCATE INTRA-LIKE INTERACTION ARRAYS

  Call allocate_core_shell_arrays()

  Call allocate_constraints_arrays()
  Call allocate_pmf_arrays()

  Call allocate_rigid_bodies_arrays()

  Call allocate_tethers_arrays()

  Call allocate_bonds_arrays()
  Call allocate_angles_arrays()
  Call allocate_dihedrals_arrays()
  Call allocate_inversions_arrays()

  Call allocate_mpoles_arrays()

! ALLOCATE INTER-LIKE INTERACTION ARRAYS

  Call allocate_vdw_arrays()
  Call allocate_metal_arrays()
  Call allocate_tersoff_arrays()
  Call allocate_three_body_arrays()
  Call allocate_four_body_arrays()

  Call allocate_external_field_arrays()

! ALLOCATE RDF, Z-DENSITY, STATISTICS & GREEN-KUBO ARRAYS

  Call allocate_rdf_arrays()
  Call allocate_z_density_arrays()
  Call allocate_statistics_arrays()
  Call allocate_greenkubo_arrays()

! READ SIMULATION CONTROL PARAMETERS

  Call read_control                                    &
           (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,        &
           rcut,rpad,rvdw,rbin,nstfce,alpha,width,     &
           l_exp,lecx,lfcap,l_top,lzero,lmin,          &
           ltgaus,ltscal,lvar,leql,lpse,               &
           lfce,lpana,lrdf,lprdf,lzdn,lpzdn,           &
           lvafav,lpvaf,ltraj,ldef,lrsd,               &
           nx,ny,nz,imd,tmd,emd,vmx,vmy,vmz,           &
           temp,press,strext,keyres,                   &
           tstep,mndis,mxdis,mxstp,nstrun,nsteql,      &
           keymin,nstmin,min_tol,nstgaus,nstscal,      &
           keyens,iso,taut,chi,soft,gama,taup,tai,ten, &
           keypse,wthpse,tmppse,                       &
           fmax,nstbpo,intsta,keyfce,epsq,             &
           rlx_tol,mxshak,tolnce,mxquat,quattol,       &
           nstbnd,nstang,nstdih,nstinv,nstrdf,nstzdn,  &
           nstmsd,istmsd,nstraj,istraj,keytrj,         &
           nsdef,isdef,rdef,nsrsd,isrsd,rrsd,          &
           ndump,pdplnc,timjob,timcls)

! READ SIMULATION FORCE FIELD

  Call read_field                       &
           (l_str,l_top,l_n_v,          &
           rcut,rvdw,rmet,width,temp,   &
           keyens,keyfce,keyshl,        &
           lecx,lbook,lexcl,            &
           rcter,rctbp,rcfbp,           &
           atmfre,atmfrz,megatm,megfrz, &
           megshl,megcon,megpmf,megrgd, &
           megtet,megbnd,megang,megdih,meginv)


If(l_jack .or. l_block) then
  Call allocate_block_average_array(nstrun)
End If

! If using induced dipoles then read in atomic polarizability

!  If (induce) Call read_polarity()

! CHECK MD CONFIGURATION

  Call check_config(levcfg,l_str,lpse,keyens,iso,keyfce,keyres,megatm)

  Call gtime(timelp)
  If (idnode == 0) Then
     Write(nrite,'(/,1x,a)') "*** all reading and connectivity checks DONE ***"
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
  End If

! l_org: translate CONFIG into CFGORG and exit gracefully

  If (l_org) Then
     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(/,1x,a)') "*** Translating the MD system along a vector (CONFIG to CFGORG) ***"
        Write(nrite,'(1x,a)') "*** ... ***"
     End If

     Call origin_config(megatm)

     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ALL DONE ***"
        Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
     End If
  End If

! l_scl: rescale CONFIG to CFGSCL and exit gracefully

  If (l_scl) Then
     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(/,1x,a)') "*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***"
        Write(nrite,'(1x,a)') "*** ... ***"
     End If

     Call scale_config(megatm)

     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ALL DONE ***"
        Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
     End If
  End If

! l_his: generate HISTORY and exit gracefully

  If (l_his) Then
     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(/,1x,a)') "*** Generating a zero timestep HISTORY frame of the MD system ***"
        Write(nrite,'(1x,a)') "*** ... ***"
     End If

! Nail down necessary parameters

     nstraj = 0 ; istraj = 1 ; keytrj = 0  ! default trajectory
     nstep  = 0                            ! no steps done
     time   = 0.0_wp                       ! time is not relevant
     Call trajectory_write(keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)

     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ALL DONE ***"
        Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
     End If
  End If

! Expand current system if opted for

  If (l_exp) Call system_expand(l_str,rcut,nx,ny,nz,megatm)

! EXIT gracefully

  If (l_trm) Then
     If (idnode == 0) Write(nrite,'(/,1x,a)') "*** Exiting gracefully ***"
     Go To 10
  End If

! READ REVOLD (thermodynamic and structural data from restart file)

  Call system_init                                                 &
           (levcfg,rcut,rvdw,rbin,rmet,lrdf,lzdn,keyres,megatm,    &
           time,tmst,nstep,tstep,chit,cint,chip,eta,virtot,stress, &
           vircon,strcon,virpmf,strpmf,elrc,virlrc,elrcm,vlrcm)

! SET domain borders and link-cells as default for new jobs
! exchange atomic data and positions in border regions

  Call set_halo_particles(rlnk,keyfce)

  Call gtime(timelp)
  If (idnode == 0) Then
     Write(nrite,'(/,1x,a)') "*** initialisation and haloing DONE ***"
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
  End If

! For any intra-like interaction, construct book keeping arrays and
! exclusion arrays for overlapped two-body inter-like interactions

  If (lbook) Then
     Call build_book_intra              &
           (l_str,l_top,lsim,dvar,      &
           megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,        &
           megrgd,degrot,degtra,        &
           megtet,megbnd,megang,megdih,meginv)
     If (mximpl > 0) Call build_tplg_intra()
     If (lexcl) Call build_excl_intra(lecx)
  Else
     Call report_topology                &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,megrgd,  &
           megtet,megbnd,megang,megdih,meginv)

! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS if no longer needed

     If (lsim) Then
        Call deallocate_core_shell_arrays()

        Call deallocate_constraints_arrays()
        Call deallocate_pmf_arrays()

        Call deallocate_rigid_bodies_arrays()

        Call deallocate_tethers_arrays()

        Call deallocate_bonds_arrays()
        Call deallocate_angles_arrays()
        Call deallocate_dihedrals_arrays()
        Call deallocate_inversions_arrays()
     End If
  End If

  Call gtime(timelp)
  If (idnode == 0) Then
     Write(nrite,'(/,1x,a)') "*** bookkeeping DONE ***"
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
  End If

! set and halo rotational matrices and their infinitesimal rotations

  If (mximpl > 0) Call mpoles_rotmat_set_halo()

! Make first check on VNL conditioning

  Call vnl_check(l_str,m_rgd,rcut,rpad,rlnk,width)

! SET initial system temperature

  Call set_temperature               &
           (levcfg,temp,keyres,      &
           lmin,nstep,nstrun,nstmin, &
           mxshak,tolnce,keyshl,     &
           atmfre,atmfrz,            &
           megshl,megcon,megpmf,     &
           megrgd,degtra,degrot,     &
           degfre,degshl,sigma,engrot)

  Call gtime(timelp)
  If (idnode == 0) Then
     Write(nrite,'(/,1x,a)') "*** temperature setting DONE ***"
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
  End If

! Frozen atoms option

  Call freeze_atoms()

! Cap forces in equilibration mode

  If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,temp)

! PLUMED initialisation or information message

  If (l_plumed) Call plumed_init(megatm,tstep,temp)

! Print out sample of initial configuration on node zero

  If (idnode == 0) Then
     Write(nrite,"(/,/,1x,'sample of starting configuration on node zero',/)")

     If (levcfg <= 1) Write(nrite,"(8x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)', &
        & 7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',/,/)")

     If (levcfg == 2) Write(nrite,"(8x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)', &
        & 7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',                               &
        & 7x,'fx(i)',7x,'fy(i)',7x,'fz(i)',/,/)")

     j=(natms+19)/20
     If (j > 0) Then
        Do i=1,natms,j
           If (levcfg <= 1) Write(nrite,"(1x,i8,1p,3e12.4,3e12.4,3e12.4)") &
              ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)

           If (levcfg == 2) Write(nrite,"(1x,i8,1p,3e12.4,3e12.4,3e12.4)") &
              ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i),fxx(i),fyy(i),fzz(i)
        End Do
     End If
  End If

! Indicate nodes mapped on vacuum (no particles)

  j=0
  If (natms == 0) Then
     j=1
     Write(nrite,'(/,1x,a,i0,a,/)') '*** warning - node ', idnode, ' mapped on vacuum (no particles) !!! ***'
  End If
  If (mxnode > 1) Call gsum(j)
  If (j > 0) Call warning(2,Real(j,wp),Real(mxnode,wp),0.0_wp)

! Initialise kinetic stress and energy contributions,
! energy(or stress) and virial accumulators for rigid bodies,
! electrostatics, short-range potentials, tersoff potentials
! three- and four-body potentials, core-shell and tether
! units, bonds, angles, dihedrals, inversions and field terms

  strkin = 0.0_wp ; engke  = 0.0_wp
  strknf = 0.0_wp ; strknt = 0.0_wp


  strcom = 0.0_wp
  vircom = 0.0_wp


  engcpe = 0.0_wp
  vircpe = 0.0_wp

  engsrp = 0.0_wp
  virsrp = 0.0_wp

  engter = 0.0_wp
  virter = 0.0_wp

  engtbp = 0.0_wp
  virtbp = 0.0_wp

  engfbp = 0.0_wp
  virfbp = 0.0_wp


  shlke  = 0.0_wp
  engshl = 0.0_wp
  virshl = 0.0_wp

  engtet = 0.0_wp
  virtet = 0.0_wp


  engbnd = 0.0_wp
  virbnd = 0.0_wp

  engang = 0.0_wp
  virang = 0.0_wp

  engdih = 0.0_wp
  virdih = 0.0_wp

  enginv = 0.0_wp
  virinv = 0.0_wp


  engfld = 0.0_wp
  virfld = 0.0_wp

! Initialise conserved quantity (other than K + U)

  consv = 0.0_wp


! start-up time when forces are not recalculated

  Call gtime(timelp)
  If (idnode == 0) &
     Write(nrite,'(/,/,/,1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp

! Now you can run fast, boy

  If (l_fast) Call gsync(l_fast)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  If (lsim) Then
     If (l_vv) Then
        Call w_md_vv()
     Else
        Call w_md_lfv()
     End If
  Else
     If (lfce) Then
        Call w_replay_historf()
     Else
        Call w_replay_history()
     End If
  End If


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Report termination of the MD simulation

  If (idnode == 0) Write(nrite,"(/,/,1x,'run terminating...  ',  &
     & 'elapsed cpu time: ', f12.3, ' sec, job time: ', f12.3,   &
     & ' sec, close time: ', f12.3, ' sec',/)") timelp,timjob,timcls

! Print out sample of final configuration on node zero

  If (idnode == 0) Then
     Write(nrite,"(/,/,1x,'sample of final configuration on node zero',/)")

     Write(nrite,"(8x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)', &
          & 7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',            &
          & 7x,'fx(i)',7x,'fy(i)',7x,'fz(i)',/,/)")

     j=(natms+19)/20
     If (j > 0) Then
        Do i=1,natms,j
           If (levcfg <= 1) Write(nrite,"(1x,i8,1p,3e12.4,3e12.4,3e12.4)") &
              ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)

           If (levcfg == 2) Write(nrite,"(1x,i8,1p,3e12.4,3e12.4,3e12.4)") &
              ltg(i),xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i),fxx(i),fyy(i),fzz(i)
        End Do
     End If
  End If

! Save restart data (final)

  If (.not.l_tor) Call system_revive                          &
           (rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

! Produce summary of simulation

  If (.not.lsim) tstep=tsths ! tstep for 'replay history'

  Call statistics_result                                        &
           (rcut,lmin,lpana,lrdf,lprdf,lzdn,lpzdn,lvafav,lpvaf, &
           nstrun,keyens,keyshl,megcon,megpmf,iso,              &
           press,strext,nstep,tstep,time,tmst)

10 Continue

! PLUMED finalisation

  If (l_plumed) Call plumed_finalize()

! Ask for reference in publications

  If (idnode == 0) Then
     Write(nrite,'(/,/,6(1x,a,/),1x,a)') &
  "*************************************************************************************************************************", &
  "**************                                                                                             **************", &
  "**************  Thank you for using the DL_POLY_4 package in your work.  Please, acknowledge our efforts   **************", &
  "**************                                                                                             **************", &
  "**************  by including the following references when publishing data obtained using DL_POLY_4:       **************", &
  "**************                                                                                             **************", &
  "**************  I.T. Todorov, W. Smith, K. Trachenko & M.T. Dove, `J. Mater. Chem.', 16, 1911-1918 (2006)  **************"

     If (keyfce == 2) Write(nrite,'(1x,a)') &
  "**************  I.J. Bush, I.T. Todorov & W. Smith, `Comp. Phys. Commun.', 175, 323-329 (2006)             **************"

     If (mximpl > 0)  Write(nrite,'(1x,a)') &
  "**************  H.A. Boateng & I.T. Todorov, `J. Chem. Phys.', 142, 034117 (2015)                          **************"

     Write(nrite,'(2(1x,a,/))') &
  "**************                                                                                             **************", &
  "*************************************************************************************************************************"
  End If

! Get just the one number to compare against

  If (idnode == 0 .and. l_eng) Write(nrite,"(/,1x,a,1p,e20.10)") "TOTAL ENERGY: ", stpval(1)

! Close output channel

  If (idnode == 0 .and. (.not.l_scr)) Close(Unit=nrite)

! Terminate job

  If (mxnode > 1) Call gsync()
  Call exit_comms()

! Create wrappers for the MD cycle in VV, LFV and replay history

Contains

  Subroutine w_impact_option()
    Include 'w_impact_option.f90'
  End Subroutine w_impact_option

  Subroutine w_calculate_forces()
    Include 'w_calculate_forces.f90'
  End Subroutine w_calculate_forces

  Subroutine w_refresh_mappings()
    Include 'w_refresh_mappings.f90'
  End Subroutine w_refresh_mappings

  Subroutine w_at_start_vv()
    Include 'w_at_start_vv.f90'
  End Subroutine w_at_start_vv

  Subroutine w_integrate_vv(isw)
    Integer, Intent( In    ) :: isw ! used for vv stage control

    Include 'w_integrate_vv.f90'
  End Subroutine w_integrate_vv

  Subroutine w_at_start_lfv()
    Include 'w_at_start_lfv.f90'
  End Subroutine w_at_start_lfv

  Subroutine w_integrate_lfv()
    Include 'w_integrate_lfv.f90'
  End Subroutine w_integrate_lfv

  Subroutine w_kinetic_options()
    Include 'w_kinetic_options.f90'
  End Subroutine w_kinetic_options

  Subroutine w_statistics_report()
    Include 'w_statistics_report.f90'
  End Subroutine w_statistics_report

  Subroutine w_write_options()
    Include 'w_write_options.f90'
  End Subroutine w_write_options

  Subroutine w_refresh_output()
    Include 'w_refresh_output.f90'
  End Subroutine w_refresh_output

  Subroutine w_md_vv()
    Include 'w_md_vv.f90'
  End Subroutine w_md_vv

  Subroutine w_md_lfv()
    Include 'w_md_lfv.f90'
  End Subroutine w_md_lfv

  Subroutine w_replay_history()
    Logical,     Save :: newjb = .true.
    Real( Kind = wp ) :: tmsh        ! tmst replacement
    Integer           :: nstpe,nstph ! nstep replacements
    Integer           :: exout       ! exit indicator for reading

    Include 'w_replay_history.f90'
  End Subroutine w_replay_history

  Subroutine w_replay_historf()
    Logical,     Save :: newjb = .true.
    Real( Kind = wp ) :: tmsh        ! tmst replacement
    Integer           :: nstpe,nstph ! nstep replacements
    Integer           :: exout       ! exit indicator for reading

    Include 'w_replay_historf.f90'
  End Subroutine w_replay_historf

End Program dl_poly
