Module drivers
  use kinds, Only : wp,wi
  Use comms, Only : comms_type,gsum
  use kinetics, Only : kinstresf, kinstrest, kinstress,getknr,getvom, getknr
  Use configuration, Only : configuration_type
  Use particle, Only : corePart
  Use rigid_bodies, Only : rigid_bodies_type,getrotmat
  Use setup, Only : boltz,mxatms,zero_plus
  Use angles, Only : angles_type
  Use dihedrals, Only : dihedrals_type
  Use inversions, Only : inversions_type
  Use core_shell,  Only : core_shell_type,SHELL_ADIABATIC 

  Use impacts, Only : impact_type, impact
  Use errors_warnings, Only : error,warning,info
  Use shared_units, Only : update_shared_units,update_shared_units_int
  Use numerics, Only : local_index,images,dcell,invert,box_mueller_saru3
  Use thermostat, Only : thermostat_type
  Use statistics, Only : stats_type
  Use domains, Only : domains_type
  Implicit None
  Private
  Public :: w_impact_option
Contains

  Subroutine w_impact_option(levcfg,nstep,nsteql,rigid,cshell,stats,impa,config,comm)

    Integer( Kind = wi ),   Intent( InOut ) :: levcfg,nstep,nsteql
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type(stats_type)   ,   Intent( InOut ) :: stats
    Type(core_shell_type)   ,   Intent( InOut ) :: cshell
    Type(impact_type)   ,   Intent( InOut ) :: impa
    Type( configuration_type ), Intent( InOut ) :: config
    Type(comms_type)    ,   Intent( InOut ) :: comm

    Character( Len = 256 ) :: messages(6)

!!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

! Apply impact
! levcfg == 2 avoids application twice when tmd happens at (re)start for VV

     If (nstep == impa%tmd .and. levcfg == 2) Then
       Write(messages(1),'(a)') ''
       Write(messages(2),'(a)') 'initiating IMPACT:'
       Write(messages(3),'(a,i10)') 'particle (index): ', impa%imd
       Write(messages(4),'(a,i10)') 'timestep (steps): ', impa%tmd
       Write(messages(5),'(a,1p,e12.5)') 'energy   (keV):   ', impa%emd
       Write(messages(6),'(a,1p,3e12.4)') 'v-r(x,y,z):       ', impa%vmx, impa%vmy, impa%vmz
       Call info(messages,6,.true.)

        If (nstep+1 <= nsteql) Call warning(380,Real(nsteql,wp),0.0_wp,0.0_wp)

        Call impact(rigid,cshell,impa,config,comm)

! Correct kinetic stress and energy

        If (rigid%total > 0) Then
           Call kinstresf(config%vxx,config%vyy,config%vzz,stats%strknf,config,comm)
           Call kinstrest(rigid,stats%strknt,comm)

           stats%strkin=stats%strknf+stats%strknt

           stats%engrot=getknr(rigid,comm)
        Else
           Call kinstress(config%vxx,config%vyy,config%vzz,stats%strkin,config,comm)
        End If
        stats%engke = 0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))
     End If

!!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
  End Subroutine w_impact_option

!  Subroutine w_refresh_mappings()
!    Include 'w_refresh_mappings.F90'
!  End Subroutine w_refresh_mappings
!
!  Subroutine w_at_start_vv()
!    Include 'w_at_start_vv.F90'
!  End Subroutine w_at_start_vv
!
!  Subroutine w_integrate_vv(isw)
!    Integer, Intent( In    ) :: isw ! used for vv stage control
!
!    Include 'w_integrate_vv.F90'
!  End Subroutine w_integrate_vv
!
!  Subroutine w_kinetic_options()
!    Include 'w_kinetic_options.F90'
!  End Subroutine w_kinetic_options
!
!  Subroutine w_statistics_report()
!    Include 'w_statistics_report.F90'
!  End Subroutine w_statistics_report
!
!  Subroutine w_write_options()
!    Include 'w_write_options.F90'
!  End Subroutine w_write_options
!
!  Subroutine w_refresh_output()
!    Include 'w_refresh_output.F90'
!  End Subroutine w_refresh_output
!
!  Subroutine w_md_vv()
!    Include 'w_md_vv.F90'
!  End Subroutine w_md_vv
!
!  Subroutine w_replay_history()
!    Logical,     Save :: newjb = .true.
!    Real( Kind = wp ) :: tmsh        ! tmst replacement
!    Integer           :: nstpe,nstph ! nstep replacements
!    Integer           :: exout       ! exit indicator for reading
!
!    Include 'w_replay_history.F90'
!  End Subroutine w_replay_history
!
!  Subroutine w_replay_historf()
!    Logical,     Save :: newjb = .true.
!    Real( Kind = wp ) :: tmsh        ! tmst replacement
!    Integer           :: nstpe,nstph ! nstep replacements
!    Integer           :: exout       ! exit indicator for reading
!
!    Include 'w_replay_historf.F90'
!  End Subroutine w_replay_historf
End Module drivers
