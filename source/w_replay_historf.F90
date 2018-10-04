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

  rsdc%nsrsd = 0
  rsdc%isrsd = 1

! intramolecular PDF analysis for every entry in HISTORF
! enforce printing and collection if the calculation exists

  stats%lpana=(config%mxgana > 0)
  flow%freq_bond = 1
  flow%freq_angle = 1
  flow%freq_dihedral = 1
  flow%freq_inversion = 1

! rdf%rdf and z-density detection for every entry in HISTORF
! enforce printing and collection if the calculation exists

  rdf%l_print=rdf%l_collect ; rdf%freq = 1
  zdensity%l_print=zdensity%l_collect ; zdensity%frequency = 1

! Calculate kinetic tensor and energy at restart as it may not exists later

  If (rigid%total > 0) Then
     Call kinstresf(config%vxx,config%vyy,config%vzz,stat%strknf,config,comm)
     Call kinstrest(rigid,stat%strknt,comm)

     stat%strkin=stat%strknf+stat%strknt
  Else
     Call kinstress(config%vxx,config%vyy,config%vzz,stat%strkin,config,comm)
  End If
  stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

  nstpe = flow%step
  nstph = 0 ! HISTORF trajectory points counter
  Do
     Call allocate_statistics_connect(cnfig%mxatdm,stat)
10   Continue
     If (nstph > nstpe) Call statistics_connect_set(config,neigh%cutoff_extended,cnfig%mxatdm,msd_data%l_msd,stat,domain,comm)

! Make a move - Read a frame

     Call read_history(flow%strict,files(FILE_HISTORF)%filename,megatm,levcfg,dvar, &
       flow%step,thermo%tstep,flow%time,exout,io,traj,sites,domain,config,files,comm)

     If (traj%restart) Then
        traj%restart = .false.

        flow%start_time=flow%time
        tmsh=0.0_wp ! flow%start_time substitute
     End If

     If (exout == 0) Then
        nstph=nstph+1

        If (nstph <= nstpe) Then

! Deal with restarts but remember the old cell parameters

           stat%clin=config%cell
           Go To 10

        Else

! CHECK MD CONFIGURATION

           Call check_config(config,levcfg,electro%key,megatm,thermo,sites,flow,comm)

! First frame positions (for estimates of MSD when levcfg==0)

           If (nstph == 1) Then
              Do i=1,config%natms
                 stat%xin(i)=config%parts(i)%xxx
                 stat%yin(i)=config%parts(i)%yyy
                 stat%zin(i)=config%parts(i)%zzz
              End Do
              stat%clin=config%cell
!              xin(natms+1: ) = 0.0_wp
!              yin(natms+1: ) = 0.0_wp
!              zin(natms+1: ) = 0.0_wp
              Call statistics_connect_set(config,neigh%cutoff_extended,cnfig%mxatdm,msd_data%l_msd,stat,domain,comm)
           End If

! get xto/xin/msdtmp arrays sorted

           Call statistics_connect_frames(config,megatm,cnfig%mxatdm,msd_data%l_msd,stat,domain,comm)
           Call deallocate_statistics_connect(stat)

! SET domain borders and link-cells as default for new jobs
! exchange atomic data and positions in border regions

           Call set_halo_particles(electro%key,neigh,sites,mpoles,domain,config,ewld,kim_data,comm)

! For any intra-like interaction, construct book keeping arrays and
! exclusion arrays for overlapped two-body inter-like interactions

           If (flow%book) Then
             Call build_book_intra(flow%strict,flow%print_topology,flow%simulation,dvar,megatm,megfrz,atmfre, &
               atmfrz,degrot,degtra,flow,cshell,cons,pmf,bond,angle,dihedral, &
               inversion,tether,neigh,sites,mpoles,rigid,domain,config,comm)
              If (flow%exclusions) Then
                Call build_excl_intra(electro%lecx,cshell,cons,bond,angle,dihedral, &
                  inversion,neigh,rigid,config,comm)
              End If
           End If

! Evaluate forces, flow%newjob must always be true for vircom evaluation

           Call w_calculate_forces(cnfig,flow,io,cshell,cons,pmf,stat,plume,pois,bond,angle,dihedral, &
             inversion,tether,threebody,neigh,sites,vdws,tersoffs,fourbody,rdf, &
             netcdf,minim,mpoles,ext_field,rigid,electro,domain,kim_data,tmr)

! Evaluate kinetics if available

           If (levcfg > 0 .and. levcfg < 3) Then
              If (thermo%l_zero .and. &
                flow%step <= flow%equil_steps .and. &
                Mod(flow%step+1-flow%equil_steps,thermo%freq_zero) == 0) Then
                Call zero_k_optimise(stat,rigid,config,comm)
              End If

              If (thermo%l_zero .and. flow%step <= flow%equil_steps) Then
                Call zero_k_optimise(stat,rigid,config,comm)
              End If

! Calculate kinetic stress and energy if available

              If (rigid%total > 0) Then
                Call rigid_bodies_quench(rigid,domain,config,comm)

                 Call kinstresf(config%vxx,config%vyy,config%vzz,stat%strknf,config,comm)
                 Call kinstrest(rigid,stat%strknt,comm)

                 stat%strkin=stat%strknf+stat%strknt

                 stat%engrot=getknr(rigid,comm)
                 Call rigid_bodies_str_ss(stat%strcom,rigid,config,comm)
                 stat%vircom=-(stat%strcom(1)+stat%strcom(5)+stat%strcom(9))
              Else
                 Call kinstress(config%vxx,config%vyy,config%vzz,stat%strkin,config,comm)
              End If
              stat%engke = 0.5_wp*(stat%strkin(1)+stat%strkin(5)+stat%strkin(9))

! Apply kinetic options

              Call w_kinetic_options(cshell,cons,pmf,stat,sites,ext_field,domain,seed)

! Get core-shell kinetic energy for adiabatic shell model

              If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC) Then
                Call core_shell_kinetic(config,stat%shlke,cshell,domain,comm)
              End If
           End If

! Get complete stress tensor

           stat%strtot = stat%strcon + stat%strpmf + stat%stress + stat%strkin + stat%strcom + stat%strdpd

! Calculate physical quantities and collect statistics,
! accumulate z-density if needed
! (flow%step->nstph,tstep->tsths,flow%start_time->tmsh)

           tsths=Max(thermo%tstep ,(flow%time-tmsh) / Real(Merge( nstph-1, 1, nstph > 2), wp))

! Collect VAF if kinetics is available

           Call vaf_collect(config,sites%mxatyp,flow%equilibration,flow%equil_steps,nstph-1,flow%time,green,comm)

           Call statistics_collect        &
           (config,flow%simulation,flow%equilibration,flow%equil_steps,msd_data%l_msd, &
           flow%restart_key,      &
           degfre,degshl,degrot,          &
           nstph,tsths,flow%time,tmsh,         &
           cnfig%mxatdm,rdf%max_grid,stat,thermo,&
           zdensity,sites,files,comm)

! line-printer output
! Update cpu flow%time

           Call gtime(tmr%elapsed)
           If (flow%new_page()) Then
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

           Write(messages(1),'(i13,1p,9e12.4)')flow%step,stat%stpval(1:9)
           Write(messages(2),'(f13.5,1p,9e12.4)')flow%time,stat%stpval(10:18)
           Write(messages(3),'(0p,f13.3,1p,9e12.4)') tmr%elapsed,stat%stpval(19:27)
           Write(messages(4),'(a)')''
           Call info(messages,4,.true.)

           Write(messages(1),'(6x,a7,1p,9e12.4)') 'rolling',stat%ravval(1:9)
           Write(messages(2),'(5x,a8,1p,9e12.4)') 'averages',stat%ravval(10:18)
           Write(messages(3),'(13x,9e12.4)') stat%ravval(19:27)
           Write(messages(4),'(a)') Repeat('-',130)
           Call info(messages,4,.true.)

           If (nstph /= 0) Then
             Call flow%line_printed()
           End If

! Write HISTORY, DEFECTS, MSDTMP, DISPDAT & VAFDAT_atom-types

           If (traj%ltraj) Call trajectory_write(flow%restart_key,megatm,flow%step,thermo%tstep,flow%time, &
             io,stat%rsd,netcdf,config,traj,files,comm)
           If (dfcts(1)%ldef)Then
             Call defects_write(flow%restart_key,thermo%ensemble,flow%step,thermo%tstep,flow%time,io,cshell, &
               dfcts(1),neigh,sites,netcdf,domain,config,files,comm)
             If (dfcts(2)%ldef)Then
               Call defects_write(flow%restart_key,thermo%ensemble,flow%step,thermo%tstep,flow%time,io,cshell, &
                 dfcts(2),neigh,sites,netcdf,domain,config,files,comm)
             End If
           End If
           If (msd_data%l_msd) Then
             Call msd_write(config,flow%restart_key,megatm,flow%step,thermo%tstep,flow%time,stat%stpval, &
               sites%dof_site,io,msd_data,files,comm)
           End If
           If (rsdc%lrsd) Call rsd_write(flow%restart_key,flow%step,thermo%tstep,io,rsdc,flow%time,cshell,stat%rsd,config,comm)
           If (green%samp > 0) Call vaf_write & ! (flow%step->nstph,tstep->tsths,flow%start_time->tmsh)
           (config,flow%restart_key,nstph,tsths,green,sites,comm)

! Save restart data in event of system crash

            If (Mod(nstph,flow%freq_restart) == 0 .and. nstph /= flow%run_steps .and. (.not.devel%l_tor)) Then
              Call system_revive(neigh%cutoff,megatm,flow%step,flow%time,sites,io,flow%start_time, &
                stat,devel,green,thermo,bond,angle,dihedral,inversion,zdensity, &
                rdf,netcdf,config,files,comm)
            End IF

! Close and Open OUTPUT at about 'i'th print-out or 'i' minute intervals

           i=20
           If ( Mod(nstph,i) == 0 .or.                               &
                (tmr%elapsed > Real(i*60,wp) .and.                        &
                 tmr%elapsed-Real( ((Int(tmr%elapsed)/(i*60)) * i*60) , wp ) < &
                 tmr%elapsed/Real( nstph , wp) ) ) Then

              If (comm%idnode == 0) Then
                 Inquire(File=files(FILE_OUTPUT)%filename, Exist=l_out, Position=c_out)
                 Call strip_blanks(c_out)
                 Call lower_case(c_out)
                 If (l_out .and. c_out(1:6) == 'append') Then
                    Close(unit=files(FILE_OUTPUT)%unit_no)
                    Open(Newunit=files(FILE_OUTPUT)%unit_no, File=files(FILE_OUTPUT)%filename, Position='append')
                 End If
              End If

           End If
        End If

! Save last frame positions (for estimates of MSD when levcfg==0)

        Do i=1,config%natms
           stat%xin(i)=config%parts(i)%xxx
           stat%yin(i)=config%parts(i)%yyy
           stat%zin(i)=config%parts(i)%zzz
        End Do
        stat%clin=config%cell
     Else
        Exit
     End If
  End Do

! Finish with grace

  If      (exout > 0) Then ! normal exit

! recover connectivity arrays for REVCON, REVIVE and printing purposes
! read_history MUST NOT initialise R,V,F arrays!!!

     config%ltg(1:config%natms) = stat%ltg0(1:config%natms)
     config%lsa(1:config%natms) = stat%lsa0(1:config%natms)
     config%lsi(1:config%natms) = stat%lsi0(1:config%natms)

  Else If (exout < 0) Then ! abnormal exit

! If reading HISTORY finished awkwardly
! recover positions and generate kinetics

     Do i=1,config%natms
        config%parts(i)%xxx=stat%xin(i)
        config%parts(i)%yyy=stat%yin(i)
        config%parts(i)%zzz=stat%zin(i)
     End Do
     config%cell=stat%clin

     Call set_temperature(levcfg,flow%restart_key,flow%step,flow%run_steps,atmfre,atmfrz,degtra, &
       degrot,degfre,degshl,stat%engrot,sites%dof_site,cshell,stat,cons,pmf, &
       thermo,minim,rigid,domain,config,seed,comm)

  End If
  Call deallocate_statistics_connect(stat)

! Save restart data because of next action (and disallow the same in dl_poly)

  If (.not. devel%l_tor) Then
    Call system_revive(neigh%cutoff,megatm,flow%step,flow%time,sites,io,flow%start_time,stat, &
      devel,green,thermo,bond,angle,dihedral,inversion,zdensity,rdf,netcdf, &
      config,files,comm)
  End If

! step counter is data counter now, so statistics_result is triggered

  flow%step=nstph
  thermo%tstep=tsths


!!!!!!!!!!!!!!!!!!!!  W_REPLAY HISTORF INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
