!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!


! Report work

  Write(messages(1),'(a)') ''
  Write(messages(2),'(a)') '*** HISTORF will be replayed in full (with no dynamics)!!!     ***'
  Write(messages(3),'(a)') '*** Large particle displacements between frames within HISTROF ***'
  Write(messages(4),'(a)') '*** w.r.t. CONFIG at start may lead failures in parallel!!!    ***'
  Call info(messages,3,.true.)

! defect detection for every entry in HISTORF

  dfcts(:)%nsdef = 0
  dfcts(:)%isdef = 1

! MSDTMP option for every entry in HISTORF

  msd_data%start = 0
  msd_data%freq = 1

! displacement detection for every entry in HISTORF

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
  nstph = 0 ! HISTORF trajectory points counter
  Do
     Call allocate_statistics_connect(mxatdm,stat)
10   Continue
     If (nstph > nstpe) Call statistics_connect_set(neigh%cutoff_extended,mxatdm_,msd_data%l_msd,stat,comm)

! Make a move - Read a frame

     Call read_history(l_str,Trim(historf),megatm,levcfg,dvar,nstep,tstep,time,exout,sites,comm)

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

           Call check_config(levcfg,l_str,electro%key,keyres,megatm,thermo,sites,comm)

! First frame positions (for estimates of MSD when levcfg==0)

           If (nstph == 1) Then
              Do i=1,natms
                 stat%xin(i)=xxx(i)
                 stat%yin(i)=yyy(i)
                 stat%zin(i)=zzz(i)
              End Do
              stat%clin=cell
!              xin(natms+1: ) = 0.0_wp
!              yin(natms+1: ) = 0.0_wp
!              zin(natms+1: ) = 0.0_wp
              Call statistics_connect_set(neigh%cutoff_extended,mxatdm_,msd_data%l_msd,stat,comm)
           End If

! get xto/xin/msdtmp arrays sorted

           Call statistics_connect_frames(megatm,mxatdm_,msd_data%l_msd,stat,comm)
           Call deallocate_statistics_connect(stat)

! SET domain borders and link-cells as default for new jobs
! exchange atomic data and positions in border regions

           Call set_halo_particles(electro%key,neigh,sites,mpole,comm)

! For any intra-like interaction, construct book keeping arrays and
! exclusion arrays for overlapped two-body inter-like interactions

           If (lbook) Then
             Call build_book_intra(l_str,l_top,lsim,dvar,megatm,megfrz,atmfre, &
               atmfrz,degrot,degtra,cshell,cons,pmf,bond,angle,dihedral, &
               inversion,tether,neigh,sites,mpole,rigid,comm)
              If (lexcl) Then
                Call build_excl_intra(lecx,cshell,cons,bond,angle,dihedral, &
                  inversion,neigh,rigid,comm)
              End If
           End If

! Evaluate forces, newjob must always be true for vircom evaluation

           Call w_calculate_forces(cshell,cons,pmf,stat,plume,pois,bond,angle,dihedral, &
             inversion,tether,threebody,neigh,sites,vdws,tersoff,fourbody,rdf, &
             netcdf,minimise,mpole,ext_field,rigid,electro,tmr)

! Evaluate kinetics if available

           If (levcfg > 0 .and. levcfg < 3) Then
              If (thermo%l_zero .and. nstep <= nsteql .and. Mod(nstep+1-nsteql,thermo%freq_zero) == 0) Then
                Call zero_k_optimise(stat,rigid,comm)
              End If

              If (thermo%l_zero .and. nstep <= nsteql) Then
                Call zero_k_optimise(stat,rigid,comm)
              End If

! Calculate kinetic stress and energy if available

              If (rigid%total > 0) Then
                 Call rigid_bodies_quench(rigid,comm)

                 Call kinstresf(vxx,vyy,vzz,stat%strknf,comm)
                 Call kinstrest(rigid,stat%strknt,comm)

                 stat%strkin=stat%strknf+stat%strknt

                 stat%engrot=getknr(rigid,comm)
                 Call rigid_bodies_str_ss(stat%strcom,rigid,comm)
                 stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
              Else
                 Call kinstress(vxx,vyy,vzz,stat%strkin,comm)
              End If
              stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

! Apply kinetic options

              Call w_kinetic_options(cshell,cons,pmf,stat,sites,ext_field)

! Get core-shell kinetic energy for adiabatic shell model

              If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Call core_shell_kinetic(stat%shlke,cshell,comm)
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
           mxatdm_,rdf%max_grid,stat,thermo,zdensity,sites,comm)

! line-printer output
! Update cpu time

           Call gtime(tmr%elapsed)
           If (Mod(lines,npage) == 0) Then
             Write(messages(1),'(a)') Repeat('-',130)
             Write(messages(2),'(9x,a4,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
              'step','eng_tot','temp_tot','eng_cfg','eng_src','eng_cou','eng_bnd','eng_ang','eng_dih','eng_tet'
             Write(messages(3),'(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
              'time(ps)',' eng_pv','temp_rot','vir_cfg','vir_src','vir_cou','vir_bnd','vir_ang','vir_con','vir_tet'
             Write(messages(4), '(5x,a8,5x,a6,4x,a8,5x,a7,5x,a7,7x,a5,8x,a4,7x,a5,5x,a7,7x,a5)') &
               'cpu  (s)','volume','temp_shl','eng_shl','vir_shl','alpha','beta','gamma','vir_pmf','press'
             Write(messages(5),'(a)') Repeat('-',130)
             Call info(messages,5,.true.)
           End If

           Write(messages(1),'(i13,1p,9e12.4)')nstep,stat%stpval(1:9)
           Write(messages(2),'(f13.5,1p,9e12.4)')time,stat%stpval(10:18)
           Write(messages(3),'(0p,f13.3,1p,9e12.4)') tmr%elapsed,stat%stpval(19:27)
           Write(messages(4),'(a)')''
           Call info(messages,4,.true.)

           Write(messages(1),'(6x,a7,1p,9e12.4)') 'rolling',stat%ravval(1:9)
           Write(messages(2),'(5x,a8,1p,9e12.4)') 'averages',stat%ravval(10:18)
           Write(messages(3),'(13x,9e12.4)') stat%ravval(19:27)
           Write(messages(4),'(a)') Repeat('-',130)
           Call info(messages,4,.true.)

           If (nstph /= 0) lines=lines+1

! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

           If (ltraj) Call trajectory_write &
           (keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time,stat%rsd,netcdf,comm)
           If (dfcts(1)%ldef)Then
             Call defects_write &
             (keyres,thermo%ensemble,nstep,tstep,time,cshell,dfcts(1),neigh,sites,netcdf,comm)
             If (dfcts(2)%ldef)Then
               Call defects_write &
               (keyres,thermo%ensemble,nstep,tstep,time,cshell,dfcts(2),neigh,sites,netcdf,comm)
             End If
           End If
           If (msd_data%l_msd) Call msd_write &
             (keyres,megatm,nstep,tstep,time,stat%stpval,sites%dof_site,msd_data,comm)
           If (lrsd) Call rsd_write &
           (keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,cshell,stat%rsd,comm)
           If (green%samp > 0) Call vaf_write & ! (nstep->nstph,tstep->tsths,tmst->tmsh)
           (keyres,nstph,tsths,green,sites,comm)

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
           stat%xin(i)=xxx(i)
           stat%yin(i)=yyy(i)
           stat%zin(i)=zzz(i)
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
        xxx(i)=stat%xin(i)
        yyy(i)=stat%yin(i)
        zzz(i)=stat%zin(i)
     End Do
     cell=stat%clin

     Call set_temperature(levcfg,keyres,nstep,nstrun,atmfre,atmfrz,degtra, &
       degrot,degfre,degshl,stat%engrot,sites%dof_site,cshell,stat,cons,pmf, &
       thermo,minimise,rigid,comm)

  End If
  Call deallocate_statistics_connect(stat)

! Save restart data because of next action (and disallow the same in dl_poly)

  If (.not. devel%l_tor) Call system_revive                         &
           (neigh%cutoff,rbin,megatm,nstep,tstep,time,tmst, &
           stat,devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf,comm)

! step counter is data counter now, so statistics_result is triggered

  nstep=nstph
  tstep=tsths


!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
