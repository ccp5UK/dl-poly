Module nvt_dpd_shardlow  
    Use comms,           Only: comms_type
    Use configuration,   Only: configuration_type
    Use constraints,     Only: constraints_type
    Use core_shell,      Only: core_shell_type
    Use domains,         Only: domains_type
    Use dpd,             Only: dpd_shardlow_integrate
    Use flow_control,    Only: flow_type
    Use neighbours,      Only: neighbours_type
    Use numerics,        Only: seed_type
    Use nve,             Only: nve_0_vv, nve_1_vv
    Use pmf,             Only: pmf_type
    Use rigid_bodies,    Only: rigid_bodies_type
    Use statistics,      Only: stats_type
    Use thermostat,      Only: DPD_SECOND_ORDER,&
                               VV_FIRST_STAGE,&
                               VV_SECOND_STAGE,&
                               thermostat_type
    Use timer,           Only: timer_type
  
    Implicit None
  
    Private
  
    Public :: nvt_0_dpd_shardlow, nvt_1_dpd_shardlow
  
  Contains
  
    Subroutine nvt_0_dpd_shardlow(stage, flow, neigh, thermo, stat, rigid, domain, cnfig, seed, comm, cshell, cons, pmf, tmr)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! dl_poly_4 subroutine for integrating newtonian equations of motion in
      ! molecular dynamics - velocity verlet with dissipative particle
      ! dynamics thermostat - thermostatting dpd forces integrated using the
      ! Shardlow splitting technique
      !
      ! copyright - daresbury laboratory
      ! author    - k.a.jonathan september 2024
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Integer,                  Intent(In   ) :: stage
      Type(stats_type),         Intent(InOut) :: stat
      Type(thermostat_type),    Intent(InOut) :: thermo
      Type(neighbours_type),    Intent(In   ) :: neigh
      Type(rigid_bodies_type),  Intent(InOut) :: rigid
      Type(domains_type),       Intent(In   ) :: domain
      Type(configuration_type), Intent(InOut) :: cnfig
      Type(seed_type),          Intent(InOut) :: seed
      Type(comms_type),         Intent(InOut) :: comm
      Type(core_shell_type),    Intent(InOut) :: cshell
      Type(constraints_type),   Intent(InOut) :: cons
      Type(pmf_type),           Intent(InOut) :: pmf
      Type(timer_type),         Intent(InOut) :: tmr
      Type(flow_type),          Intent(InOut) :: flow
  
      If (stage == VV_FIRST_STAGE) Then

        ! one-off application for zeroth and first order splitting, initial application for second order
        ! velocity field change + generation of DPD virial & stat%stress due to random and drag forces
        Call dpd_shardlow_integrate(stage, flow%strict, flow%step, thermo%tstep, stat, thermo, &
                                    neigh, rigid, domain, cnfig, seed, comm)

        ! integrate equations of motion - velocity verlet stage 1
        Call nve_0_vv &
            (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
             stat%strkin, stat%engke, thermo, &
             cshell, cons, pmf, stat, domain, tmr, cnfig, comm)

      Else ! VV_SECOND_STAGE

        ! integrate equations of motion - velocity verlet stage 2
        Call nve_0_vv &
            (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
             stat%strkin, stat%engke, thermo, &
             cshell, cons, pmf, stat, domain, tmr, cnfig, comm)

        ! symmetric application for second order splitting only
        ! velocity field change + generation of DPD virial & stat%stress due to random and drag forces
        If (thermo%key_dpd == DPD_SECOND_ORDER) Then
            Call dpd_shardlow_integrate(stage, flow%strict, flow%step, thermo%tstep, stat, thermo, &
                                        neigh, rigid, domain, cnfig, seed, comm)
        End If

      End If
  
    End Subroutine nvt_0_dpd_shardlow
  
    Subroutine nvt_1_dpd_shardlow(stage, flow, neigh, thermo, stat, rigid, domain, cnfig, seed, comm, cshell, cons, pmf, tmr)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
      ! RBs, equations of motion in molecular dynamics
      ! - velocity verlet with dissipative particle dynamics thermostat - 
      ! thermostatting dpd forces integrated using the Shardlow splitting 
      ! technique
      !
      ! copyright - daresbury laboratory
      ! author    - k.a.jonathan september 2024
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Integer,                  Intent(In   ) :: stage
      Type(stats_type),         Intent(InOut) :: stat
      Type(thermostat_type),    Intent(InOut) :: thermo
      Type(neighbours_type),    Intent(In   ) :: neigh
      Type(rigid_bodies_type),  Intent(InOut) :: rigid
      Type(domains_type),       Intent(In   ) :: domain
      Type(configuration_type), Intent(InOut) :: cnfig
      Type(seed_type),          Intent(InOut) :: seed
      Type(comms_type),         Intent(InOut) :: comm
      Type(core_shell_type),    Intent(InOut) :: cshell
      Type(constraints_type),   Intent(InOut) :: cons
      Type(pmf_type),           Intent(InOut) :: pmf
      Type(timer_type),         Intent(InOut) :: tmr
      Type(flow_type),          Intent(InOut) :: flow
  
      If (stage == VV_FIRST_STAGE) Then

        ! one-off application for zeroth and first order splitting, initial application for second order
        ! velocity field change + generation of DPD virial & stat%stress due to random and drag forces
        Call dpd_shardlow_integrate(stage, flow%strict, flow%step, thermo%tstep, stat, thermo, &
                                    neigh, rigid, domain, cnfig, seed, comm)

        ! integrate equations of motion - velocity verlet stage 1
        Call nve_1_vv &
            (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
             stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
             stat%strcom, stat%vircom, thermo, cshell, cons, pmf, &
             stat, rigid, domain, tmr, cnfig, comm)

      Else ! VV_SECOND_STAGE

        ! integrate equations of motion - velocity verlet stage 2
        Call nve_1_vv &
          (stage, thermo%lvar, thermo%mndis, thermo%mxdis, thermo%mxstp, thermo%tstep, &
          stat%strkin, stat%strknf, stat%strknt, stat%engke, stat%engrot, &
          stat%strcom, stat%vircom, thermo, cshell, cons, pmf, &
          stat, rigid, domain, tmr, cnfig, comm)

        ! symmetric application for second order splitting only
        ! velocity field change + generation of DPD virial & stat%stress due to random and drag forces
        If (thermo%key_dpd == DPD_SECOND_ORDER) Then
          Call dpd_shardlow_integrate(stage, flow%strict, flow%step, thermo%tstep, stat, thermo, &
                                      neigh, rigid, domain, cnfig, seed, comm)
        End If

      End If
  
    End Subroutine nvt_1_dpd_shardlow

End Module nvt_dpd_shardlow