!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORY INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Report work

  Call info('*** HISTORY is replayed for recalculation of structural properties ***',.true.)

! Stay safe

  If (ltraj) Then
    ltraj = .false.
    Call info('*** warning - aborting printing into HISTORY while reading it ***',.true.)
  End If

! Make sure of no equilibration

  nsteql = 0
  leql   = .false.

! nullify forces
  Do i=1,mxatms
    parts(i)%fxx=0.0_wp
    parts(i)%fyy=0.0_wp
    parts(i)%fzz=0.0_wp
  End Do

! nullify all two-body force switches = just do rdf%rdf calculation

  electro%key = ELECTROSTATIC_NULL
  vdws%n_vdw = 0
  met%n_potentials = 0

! defect detection for every entry in HISTORY

  dfcts(:)%nsdef = 0
  dfcts(:)%isdef = 1

! MSDTMP option for every entry in HISTORY

  msd_data%start = 0
  msd_data%freq = 1

! displacement detection for every entry in HISTORY

  nsrsd = 0
  isrsd = 1

! intramolecular PDF analysis for every entry in HISTORF
! enforce printing and collection if the calculation exists

  lpana=(mxgana > 0)
  nstbnd = 1
  nstang = 1
  nstdih = 1
  nstinv = 1

! rdf%rdf and z-density detection for every entry in HISTORF
! enforce printing and collection if the calculation exists

  rdf%l_print=rdf%l_collect ; rdf%freq = 1
  zdensity%l_print=zdensity%l_collect ; zdensity%frequency = 1

! Calculate kinetic tensor and energy at restart as it may not exists later

  If (rigid%total > 0) Then
     Call kinstresf(vxx,vyy,vzz,stat%strknf,comm)
     Call kinstrest(rigid,stat%strknt,comm)

     stat%strkin=stat%strknf+stat%strknt
  Else
     Call kinstress(vxx,vyy,vzz,stat%strkin,comm)
  End If
  stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

  nstpe = nstep
  nstph = 0 ! trajectory points counter
  Do
     Call allocate_statistics_connect(mxatdm,stat)
10   Continue
     If (nstph > nstpe) Then
       Call statistics_connect_set(neigh%cutoff_extended,mxatdm_,msd_data%l_msd,stat,domain,parts,comm)
     End If

! Make a move - Read a frame

     Call read_history(l_str,Trim(history),megatm,levcfg,dvar,nstep,tstep,time,exout,traj,sites,domain,parts,comm)

     If (newjb) Then
        newjb = .false.

        tmst=time
        tmsh=0.0_wp ! tmst substitute
     End If

     If (exout == 0) Then
        nstph=nstph+1

        If (nstph <= nstpe) Then

! Deal with restarts but remember the old cell parameters

           stat%clin=cell
           Go To 10

        Else

! CHECK MD CONFIGURATION

           Call check_config &
           (levcfg,l_str,electro%key,keyres,megatm,thermo,sites,comm)

! First frame positions (for estimates of MSD when levcfg==0)

           If (nstph == 1) Then
              Do i=1,natms
                 stat%xin(i)=parts(i)%xxx
                 stat%yin(i)=parts(i)%yyy
                 stat%zin(i)=parts(i)%zzz
              End Do
              stat%clin=cell
!              xin(natms+1: ) = 0.0_wp
!              yin(natms+1: ) = 0.0_wp
!              zin(natms+1: ) = 0.0_wp
              Call statistics_connect_set(neigh%cutoff_extended,mxatdm,msd_data%l_msd,stat,domain,parts,comm)
           End If

! get xto/xin/msdtmp arrays sorted

           Call statistics_connect_frames(megatm,mxatdm_,msd_data%l_msd,stat,domain,comm)
           Call deallocate_statistics_connect(stat)

! SET domain borders and link-cells as default for new jobs
! exchange atomic data and positions in border regions

           Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,ewld,kim_data,comm)

! For any intra-like interaction, construct book keeping arrays and
! exclusion arrays for overlapped two-body inter-like interactions

           If (lbook) Then
             Call build_book_intra(l_str,l_top,lsim,dvar,megatm,megfrz,atmfre, &
               atmfrz,degrot,degtra,flw,cshell,cons,pmf,bond,angle,dihedral, &
               inversion,tether,neigh,sites,mpoles,rigid,domain,parts,comm)
             If (lexcl) Then
               Call build_excl_intra(lecx,cshell,cons,bond,angle,dihedral, &
                 inversion,neigh,rigid,comm)
             End If
           End If

! Accumulate RDFs if needed (nstep->nstph)
! Make sure RDFs are complete (lbook=.false. - no exclusion lists)

           If (rdf%l_collect) Then
             Call two_body_forces(pdplnc,thermo%ensemble,nstfce,.false.,megfrz, &
               leql,nsteql,nstph,cshell,stat,ewld,devel,met,pois,neigh,sites, &
               vdws,rdf,mpoles,electro,domain,parts,kim_data,tmr,comm)
           End If

! Calculate bond forces

           If (bond%total > 0 .and. bond%bin_pdf > 0) Then
              isw = 0
              Call bonds_forces(isw,stat%engbnd,stat%virbnd,stat%stress, &
              neigh%cutoff,stat%engcpe,stat%vircpe,bond,mpoles,electro,parts,comm)
           End If

! Calculate valence angle forces

           If (angle%total > 0 .and. angle%bin_adf > 0) Then
              isw = 0
              Call angles_forces(isw,stat%engang,stat%virang,stat%stress,angle,parts,comm)
           End If

! Calculate dihedral forces

           If (dihedral%total > 0 .and. dihedral%bin_adf > 0) Then
              isw = 0
              Call dihedrals_forces(isw,stat%engdih,stat%virdih,stat%stress, &
                neigh%cutoff,stat%engcpe,stat%vircpe,stat%engsrp,stat%virsrp, &
                dihedral,vdws,mpoles,electro,parts,comm)
           End If

! Calculate inversion forces

           If (inversion%total > 0 .and. inversion%bin_adf > 0) Then
              isw = 0
              Call inversions_forces(isw,stat%enginv,stat%virinv,stat%stress,inversion,parts,comm)
           End If

! Calculate kinetic stress and energy if available

           If (levcfg > 0 .and. levcfg < 3) Then
              If (rigid%total > 0) Then
                Call rigid_bodies_quench(rigid,domain,parts,comm)

                 Call kinstresf(vxx,vyy,vzz,stat%strknf,comm)
                 Call kinstrest(rigid,stat%strknt,comm)

                 stat%strkin=stat%strknf+stat%strknt

                 stat%engrot=getknr(rigid,comm)
                 If (levcfg == 2) Then
                    Call rigid_bodies_str_ss(stat%strcom,rigid,parts,comm)
                    stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
                 End If
              Else
                 Call kinstress(vxx,vyy,vzz,stat%strkin,comm)
              End If
              stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

! Apply kinetic options

              Call w_kinetic_options(cshell,cons,pmf,stat,sites,ext_field,domain,seed)

! Get core-shell kinetic energy for adiabatic shell model

              If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
                Call core_shell_kinetic(stat%shlke,cshell,domain,comm)
              End If
           End If

! Get complete stress tensor

           stat%strtot = stat%strcon + stat%strpmf + stat%stress + stat%strkin + stat%strcom + stat%strdpd

! Calculate physical quantities and collect statistics,
! accumulate z-density if needed
! (nstep->nstph,tstep->tsths,tmst->tmsh)

           tsths=Max(tstep ,(time-tmsh) / Real(Merge( nstph-1, 1, nstph > 2), wp))

! Collect VAF if kinetics is available

           Call vaf_collect(leql,nsteql,nstph-1,time,green,comm)

           Call statistics_collect        &
           (lsim,leql,nsteql,msd_data%l_msd, &
           keyres,      &
           degfre,degshl,degrot,          &
           nstph,tsths,time,tmsh,         &
           mxatdm_,rdf%max_grid,stat,thermo,&
           zdensity,sites,parts,comm)

! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

           If (ltraj) Call trajectory_write(keyres,megatm,nstep,tstep,time, &
             stat%rsd,netcdf,parts,traj,comm)
           If (dfcts(1)%ldef)Then
             Call defects_write(keyres,thermo%ensemble,nstep,tstep,time,cshell, &
               dfcts(1),neigh,sites,netcdf,domain,parts,comm)
             If (dfcts(2)%ldef)Then
               Call defects_write(keyres,thermo%ensemble,nstep,tstep,time, &
                 cshell,dfcts(2),neigh,sites,netcdf,domain,parts,comm)
             End If
           End If
           If (msd_data%l_msd) Call msd_write &
             (keyres,megatm,nstep,tstep,time,stat%stpval,sites%dof_site,msd_data,comm)
           If (lrsd) Call rsd_write &
           (keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,cshell,stat%rsd,parts,comm)
           If (green%samp > 0) Call vaf_write & ! (nstep->nstph,tstep->tsths,tmst->tmsh)
           (keyres,nstph,tsths,green,sites,comm)

! Complete time check

           Call gtime(tmr%elapsed)
           Write(messages(1),'(2(a,i10),a)') &
             'HISTORY step ',nstep,' (',nstph,' entry) processed'
           Write(messages(2),'(a,f12.3,a)') &
             'time elapsed since job start: ',tmr%elapsed,' sec'
           Call info(messages,2,.true.)

! Save restart data in event of system crash

           If (Mod(nstph,ndump) == 0 .and. nstph /= nstrun .and. (.not.devel%l_tor)) &
              Call system_revive                              &
           (neigh%cutoff,rbin,megatm,nstep,tstep,time,tmst, &
           stat,devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf,comm)

! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

           i=20
           If ( Mod(nstph,i) == 0 .or.                               &
                (tmr%elapsed > Real(i*60,wp) .and.                        &
                 tmr%elapsed-Real( ((Int(tmr%elapsed)/(i*60)) * i*60) , wp ) < &
                 tmr%elapsed/Real( nstph , wp) ) ) Then

              If (comm%idnode == 0) Then
                 Inquire(File=Trim(output), Exist=l_out, Position=c_out)
                 Call strip_blanks(c_out)
                 Call lower_case(c_out)
                 If (l_out .and. c_out(1:6) == 'append') Then
                    Close(Unit=nrite)
                    Open(Unit=nrite, File=Trim(output), Position='append')
                 End If
              End If

           End If
        End If

! Save last frame positions (for estimates of MSD when levcfg==0)

        Do i=1,natms
           stat%xin(i)=parts(i)%xxx
           stat%yin(i)=parts(i)%yyy
           stat%zin(i)=parts(i)%zzz
        End Do
        stat%clin=cell
     Else
        Exit
     End If
  End Do

! Finish with grace

  If      (exout > 0) Then ! normal exit

! recover connectivity arrays for REVCON, REVIVE and printing purposes
! read_history MUST NOT initialise R,V,F arrays!!!

     ltg(1:natms) = stat%ltg0(1:natms)
     lsa(1:natms) = stat%lsa0(1:natms)
     lsi(1:natms) = stat%lsi0(1:natms)

  Else If (exout < 0) Then ! abnormal exit

! If reading HISTORY finished awkwardly
! recover positions and generate kinetics

     Do i=1,natms
        parts(i)%xxx=stat%xin(i)
        parts(i)%yyy=stat%yin(i)
        parts(i)%zzz=stat%zin(i)
     End Do
     cell=stat%clin

     Call set_temperature(levcfg,keyres,nstep,nstrun,atmfre,atmfrz,degtra, &
       degrot,degfre,degshl,stat%engrot,sites%dof_site,cshell,stat,cons,pmf, &
       thermo,minim,rigid,domain,parts,seed,comm)

  End If
  Call deallocate_statistics_connect(stat)

! Save restart data because of next action (and disallow the same in dl_poly)

  If (.not. devel%l_tor) Call system_revive                         &
           (neigh%cutoff,rbin,megatm,nstep,tstep,time,tmst, &
           stat,devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf,comm)

! step counter is data counter now, so statistics_result is triggered

  nstep=nstph
  tstep=tsths


!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORY INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
