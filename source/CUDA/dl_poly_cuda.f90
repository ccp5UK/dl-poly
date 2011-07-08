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
! authors   - i.t.todorov & w.smith 2010
! contributors: i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Centre for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.



! SETUP MODULES

  Use kinds_f90
  Use comms_module
  Use setup_module

! IO MODULE

  Use io_module

! SITE & CONFIG MODULES

  Use site_module
  Use config_module

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

  Use vdw_module
  Use metal_module
  Use tersoff_module
  Use three_body_module
  Use four_body_module

  Use external_field_module

! STATISTICS MODULE

  Use statistics_module

! PARSE MODULE

  Use parse_module

! MSD MODULE

  Use msd_module

! KINETIC MODULE

  Use kinetic_module

! DEVELOPMENT MODULE

  Use development_module

#ifdef COMPILE_CUDA
! C INTERFACES
  Use dl_poly_cuda_module
#endif

! MAIN PROGRAM VARIABLES

  Implicit None

! newjob used for trajectory_write & defects_write

  Logical, Save :: newjob = .true.

! OUTPUT existence

  Logical               :: l_out
  Character( Len = 10 ) :: c_out

! lines and page used for printing controls

  Integer       :: lines = 0 , &
                   npage = 8

! general flags

  Logical           :: l_vv,l_n_e,l_n_r,l_n_v, &
                       l_ind,l_str,l_top,l_exp,&
                       lecx,lfcap,lzero,       &
                       lmin,ltgaus,ltscal,     &
                       lvar,leql,lpse,lsim,    &
                       lrdf,lprdf,lzdn,lpzdn,  &
                       ltraj,ldef,lrsd,        &
                       safe,lbook,lexcl,       &
                       relaxed_shl = .true.,   &
                       relaxed_min = .true.

! 'isw' is used for vv stage control

  Integer           :: i,j,isw,levcfg,imcon,nstfce,        &
                       nx,ny,nz,imd,tmd,                   &
                       keyres,nstrun,nsteql,               &
                       keymin,nstmin,nstgaus,nstscal,      &
                       keyens,iso,intsta,keypse,nstbpo,    &
                       keyfce,mxshak,mxquat,nstrdf,nstzdn, &
                       nstmsd,istmsd,nstraj,istraj,keytrj, &
                       nsdef,isdef,nsrsd,isrsd,            &
                       ndump,nstep,keyshl,                 &
                       atmfre,atmfrz,megatm,megfrz,        &
                       megshl,megcon,megpmf,megrgd,        &
                       megtet,megbnd,megang,megdih,meginv

! Degrees of freedom must be in long integers so we do 2.1x10^9 particles

  Integer(Kind=ip)  :: degfre,degshl,degtra,degrot

! elrcm,vlrcm - metal energy and virial are array-like and in metal_module

  Real( Kind = wp ) :: timelp,timjob,timcls,tstep,time,tmst,tmsth,     &
                       alpha,epsq,fmax,                                &
                       rcut,rvdw,rmet,rbin,rcter,rctbp,rcfbp,          &
                       width,mndis,mxdis,wthpse,tmppse,                &
                       rlx_tol,min_tol,tolnce,quattol,rdef,rrsd,       &
                       emd,vmx,vmy,vmz,temp,sigma,                     &
                       press,strext(1:9),ten,                          &
                       taut,soft,taup,chi,tai,                         &
                       chit,eta(1:9),chip,cint,consv,                  &
                       strtot(1:9),virtot,elrc,virlrc,                 &
                       strkin(1:9),engke,strknf(1:9),strknt(1:9),      &
                       engrot,strcom(1:9),vircom,                      &
                       engcpe,vircpe,engsrp,virsrp,                    &
                       engter,virter,engtbp,virtbp,engfbp,virfbp,      &
                       engshl,shlke,virshl,                            &
                       strcon(1:9),vircon,strpmf(1:9),virpmf,          &
                       stress(1:9),engtet,virtet,                      &
                       engbnd,virbnd,engang,virang,                    &
                       engdih,virdih,enginv,virinv,                    &
                       engfld,virfld,                                  &
                       stptmp,stpprs,stpvol,stpcfg,stpeng,stpeth,stpvir

! SET UP COMMUNICATIONS & CLOCKING

  Call init_comms()
  If (mxnode > 1) Call gsync()
  Call gtime(timelp)

  Call scan_development()

#ifdef COMPILE_CUDA
! CK: Probably this is a good place to initialise the cuda-related stuff
  Call dl_poly_cuda_initialise1(.true., .false.)
#endif

! OPEN MAIN OUTPUT CHANNEL & PRINT HEADER AND MACHINE RESOURCES

  If (idnode == 0) Then
     If (.not.l_scr) Open(Unit=nrite, File='OUTPUT', Status='replace')

     Write(nrite,'(7(1x,a,/),1x,a,i12,a,/,(1x,a,/))')                           &
          "******************************************************************", &
          "*************  stfc/ccp5  program  library  package  ** D ********", &
          "*************  daresbury laboratory general purpose  *** L *******", &
          "**         **  classical molecular dynamics program  **** \ ******", &
          "** DL_POLY **  authors:   i.t.todorov   &   w.smith  ***** P *****", &
          "**         **  contributors:  i.j.bush               ****** O ****", &
          "*************  version:  4.01.1   /   november 2010  ******* L ***", &
          "*************  Execution on ", mxnode, "    node(s)  ******** Y **", &
          "******************************************************************"
  End If

#ifdef COMPILE_CUDA
  Call start_timing_total()
#endif

! TEST I/O

  Call scan_control_io()

! DETERMINE ARRAYS' BOUNDS LIMITS & DOMAIN DECOMPOSITIONING
! (setup_module and domains_module)

  Call set_bounds                                            &
           (levcfg,imcon,l_vv,l_str,l_n_e,l_n_r,l_n_v,l_ind, &
           rcut,rvdw,rmet,rbin,nstfce,alpha,width)

! ALLOCATE SITE & CONFIG ARRAYS

  Call allocate_site_arrays()
  Call allocate_config_arrays()

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

! ALLOCATE INTER-LIKE INTERACTION ARRAYS

  Call allocate_vdw_arrays()
  Call allocate_metal_arrays()
  Call allocate_tersoff_arrays()
  Call allocate_three_body_arrays()
  Call allocate_four_body_arrays()

  Call allocate_external_field_arrays()

! ALLOCATE STATISTICS ARRAYS

  Call allocate_statistics_arrays()

! READ SIMULATION CONTROL PARAMETERS

  Call read_control                               &
           (levcfg,l_vv,l_str,l_n_e,l_n_r,l_n_v,  &
           rcut,rvdw,rbin,nstfce,alpha,width,     &
           l_exp,lecx,lfcap,l_top,lzero,lmin,     &
           ltgaus,ltscal,lvar,leql,lpse,          &
           lsim,lrdf,lprdf,lzdn,lpzdn,            &
           ltraj,ldef,lrsd,                       &
           nx,ny,nz,imd,tmd,emd,vmx,vmy,vmz,      &
           temp,press,strext,keyres,              &
           tstep,mndis,mxdis,nstrun,nsteql,       &
           keymin,nstmin,min_tol,nstgaus,nstscal, &
           keyens,iso,taut,soft,taup,chi,tai,ten, &
           keypse,wthpse,tmppse,                  &
           fmax,nstbpo,intsta,keyfce,epsq,        &
           rlx_tol,mxshak,tolnce,mxquat,quattol,  &
           nstrdf,nstzdn,                         &
           nstmsd,istmsd,nstraj,istraj,keytrj,    &
           nsdef,isdef,rdef,nsrsd,isrsd,rrsd,     &
           ndump,timjob,timcls)

! READ SIMULATION FORCE FIELD

  Call read_field                          &
           (imcon,l_n_v,l_str,l_top,       &
           rcut,rvdw,rmet,width,keyfce,    &
           lbook,lexcl,keyshl,             &
           rcter,rctbp,rcfbp,              &
           atmfre,atmfrz,megatm,megfrz,    &
           megshl,megcon,megpmf,megrgd,    &
           megtet,megbnd,megang,megdih,meginv)

! CHECK MD CONFIGURATION

  Call check_config(levcfg,imcon,l_str,lpse,keyens,keyfce,keyres,megatm)

! l_scl: rescale CONFIG to CFGSCL and exit gracefully

  If (l_scl) Then
     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(/,/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
        Write(nrite,'(1x,a)') "*** Rescaling the MD system lattice (CONFIG to CFGSCL) ***"
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ... ***"
     End If

     Call scale_config(imcon,megatm)

     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ALL DONE ***"
        Write(nrite,'(1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
     End If
  End If

! l_his: generate HISTORY and exit gracefully

  If (l_his) Then
     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(/,/,1x, "time elapsed since job start: ", f12.3, " sec")') timelp
        Write(nrite,'(1x,a)') "*** Generating a zero timestep HISTORY frame of the MD system ***"
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ... ***"
     End If

! Nail down necessary parameters

     nstraj = 0 ; istraj = 1 ; keytrj = 0  ! default trajectory
     nstep  = 0                            ! no steps done
     time   = 0.0_wp                       ! time is not relevant
     Call trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)

     Call gtime(timelp)
     If (idnode == 0) Then
        Write(nrite,'(1x,a)') "*** ... ***"
        Write(nrite,'(1x,a)') "*** ALL DONE ***"
        Write(nrite,'(1x, "time elapsed since job start: ", f12.3, " sec",/)') timelp
     End If
  End If

! Expand current system if opted for

  If (l_exp) Call system_expand(imcon,nx,ny,nz,megatm)

! EXIT gracefully

  If (l_trm) Then
     Write(nrite,'(1x,a)') "*** Exiting gracefully ***"
     Go To 10
  End If

! READ REVOLD (thermodynamic and structural data from restart file)

  Call system_init                                           &
           (levcfg,imcon,rcut,rvdw,rbin,rmet,                &
           lrdf,lzdn,keyres,megatm,                          &
           time,tmst,nstep,chit,cint,chip,eta,virtot,stress, &
           vircon,strcon,virpmf,strpmf,elrc,virlrc,elrcm,vlrcm)

! SET domain borders and link-cells as default for new jobs
! exchange atomic data and positions in border regions

  Call set_halo_particles(imcon,rcut,keyfce,lbook)

! For any intra-like interaction, construct book keeping arrays and
! exclusion arrays for overlapped two-body inter-like interactions

  If (lbook) Then
     Call build_book_intra               &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,         &
           megrgd,degrot,degtra,         &
           megtet,megbnd,megang,megdih,meginv)
     If (lexcl) Call build_excl_intra(lecx)
  Else
     Call report_topology                &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,megrgd,  &
           megtet,megbnd,megang,megdih,meginv)

! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS since no longer needed

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

! SET initial system temperature

  Call set_temperature                &
           (levcfg,imcon,temp,keyres, &
           lmin,nstep,nstrun,nstmin,  &
           mxshak,tolnce,keyshl,      &
           atmfre,atmfrz,             &
           megshl,megcon,megpmf,      &
           megrgd,degtra,degrot,      &
           degfre,degshl,sigma,engrot)

! Frozen atoms option

  Call freeze_atoms()

! Cap forces in equilibration mode

  If (nstep <= nsteql .and. lfcap) Call cap_forces(fmax,temp)

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
     Call warning(1,Real(idnode,wp),0.0_wp,0.0_wp)
  End If
  If (mxnode > 1) Call gsum(j)
  If (j > 0) Call warning(2,Real(j,wp),Real(mxnode,wp),0.0_wp)

! Initialize kinetic stress and energy contributions,
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
  If (idnode == 0) Write(nrite,'(/,/,/,1x, &
     & "time elapsed since job start: ", f12.3, " sec",/)') timelp

#ifdef COMPILE_CUDA
! cuda: do the second stage of the cuda initialisations now that the
! dl poly data structures have been populated.
  Call dl_poly_cuda_initialise2(mxatms, mxcons)
! Check for unimplemented functionality in CUDA port
  Call dl_poly_cuda_check_offload_conditions(keyfce, imcon)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  If (lsim) Then
     If (l_vv) Then
        Call md_vv()
     Else
        Call md_lfv()
     End If
  Else
     Call replay_history()
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

  Call system_revive                                                &
           (imcon,rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

! Produce summary of simulation

#ifdef COMPILE_CUDA
! CK: This is probably a good place to finalise the cuda stuff:
  Call dl_poly_cuda_finalise()
#endif

  Call statistics_result                &
           (rcut,lrdf,lprdf,lzdn,lpzdn, &
           nstrun,keyens,keyshl,iso,    &
           press,strext,nstep,tstep,time,tmst)

10 Continue

! Ask for reference in publications

  If (idnode == 0) Write(nrite,'(/,/,9(1x,a,/))') &
     "*************************************************************************************************************************", &
     "**************                                                                                             **************", &
     "**************  Thank you for using the DL_POLY_4 package in your work.  Please, acknowledge our efforts   **************", &
     "**************                                                                                             **************", &
     "**************  by including the following reference when publishing data obtained using DL_POLY_4:        **************", &
     "**************                                                                                             **************", &
     "**************  I.T. Todorov, W. Smith, K. Trachenko & M.T. Dove, `J. Mater. Chem.', 16, 1911-1918 (2006)  **************", &
     "**************                                                                                             **************", &
     "*************************************************************************************************************************"

  If (idnode == 0 .and. l_eng) Write(nrite,"(/,1x,a,1p,e20.10)") "TOTAL ENERGY: ", stpval(1)

#ifdef COMPILE_CUDA
  Call dump_timings()
#endif

! Close output channel

  If (idnode == 0 .and. (.not.l_scr)) Close(Unit=nrite)

! Terminate job

  If (mxnode > 1) Call gsync()
  Call exit_comms()

! Create interfaces to md_step in either Verlet flavour

Contains

  Subroutine md_vv()
    Include 'md_vv.f90'
  End Subroutine md_vv

  Subroutine md_lfv()
    Include 'md_lfv.f90'
  End Subroutine md_lfv

  Subroutine replay_history()
    Include 'replay_history.f90'
  End Subroutine replay_history

End Program dl_poly
