Module thermostat
  Use kinds, Only : wp,wi
  Use comms, Only : gmax,comms_type
  Use errors_warnings, Only : error
  Use particle, Only: corePart
  Implicit None

  Private

  !> Type containing thermostat and barostat variables
  Type, Public :: thermostat_type
    Private

    !> Ensemble key
    Integer( Kind = wi ), Public :: ensemble
    !> Flag for variable config%cell size e.g. NPT ensembles
    Logical, Public :: variable_cell = .false.
    !> Flag for anisotropic pressure
    Logical, Public :: anisotropic_pressure = .false.

    !> Simulation temperature
    Real( Kind = wp ), Public :: temp
    !> Simulation pressure
    Real( Kind = wp ), Public :: press
    !> Simulation stress
    Real( Kind = wp ), Public :: stress(1:9)
    !> Average total energy due to equipartition
    Real( Kind = wp ), Public :: sigma


    !> Thermostat relaxation time
    Real( Kind = wp ), Public :: tau_t
    !> Barostat relxation time
    Real( Kind = wp ), Public :: tau_p

    !> Surface tension
    Real( Kind = wp ), Public :: tension

    !> Constraint type for anisotropic barostats
    Integer( Kind = wi ), Public :: iso

    !> Andersen thermostat softness
    Real( Kind = wp ), Public :: soft

    !> Langevin switch
    Logical, Public :: l_langevin = .false.

    !> Gentle Stochastic dynamics (Langevin) thermostat friction
    Real( Kind = wp ), Public :: gama

    !> Stochastic Dynamics (SD Langevin) thermostat friction
    Real( Kind = wp ), Public :: chi
    !> Inhomogeneous Stochastic Dynamics (SD Langevin)
    !> thermostat/electron-phonon friction
    Real( Kind = wp ), Public :: chi_ep
    !> Inhomogeneous Stochastic Dynamics (SD Langevin)
    !> thermostat/electronic stopping friction
    Real( Kind = wp ), Public :: chi_es
    !> Stochastic Dynamics (SD Langevin) barostat friction
    Real( Kind = wp ), Public :: tai
    !> Square of cutoff velocity for inhomogeneous Langevin thermostat and ttm
    Real( Kind = wp ), Public :: vel_es2

    !> Instantaneous thermostat friction
    Real( Kind = wp ), Public :: chi_t
    !> Instantaneous barostat friction
    Real( Kind = wp ), Public :: chi_p
    !> Friction integral for thermostat/barostat
    Real( Kind = wp ), Public :: cint
    !> Cell parameter scaling factor for barostats
    Real( Kind = wp ), Public :: eta(1:9)

    !> DPD switch
    !>
    !> - 0 no DPD
    !> - 1 first order splitting
    !> - 2 second order splitting
    Integer, Public :: key_dpd = 0
    !> DPD drag?
    Real( Kind = wp ), Allocatable, Public :: gamdpd(:)
    Real( Kind = wp ), Allocatable, Public :: sigdpd(:)

    !> Pseudo thermostat switch
    Logical, Public :: l_stochastic_boundaries
    !> Pseudo thermostat type
    !>
    !> - 0 Langevin + direct temperature scaling
    !> - 1 Langevin temperature scaling
    !> - 2 Gaussian temperature scaling
    !> - 3 direct temperature scaling
    Integer, Public :: key_pseudo
    !> Pseudo thermostat temperature
    Real( Kind = wp ), Public :: temp_pseudo
    !> Pseudo thermostat thickness
    Real( Kind = wp ), Public :: width_pseudo

    !> Temperature scaling switch
    Logical, Public :: l_tscale
    !> Temperature scaling frequency
    Integer, Public :: freq_tscale

    !> Temperature regaussing switch
    Logical, Public :: l_tgaus
    !> Temperature regaussing frequency
    Integer, Public :: freq_tgaus

    !> Zero temperature optimisation switch
    Logical, Public :: l_zero
    !> Zero temperature regaussing frequency
    Integer, Public :: freq_zero

    Logical, Public :: newjob = .true.
    Logical, Public :: newjob_nst_scl = .true.
    Logical, Public :: newjob_npt_scl = .true.
    Integer, Public :: mxiter,mxkit,kit
    Logical, Public :: unsafe = .false.
    Real( Kind = wp ), Public :: volm0,elrc0,virlrc0,h_z,cell0(1:9)
    Real( Kind = wp ), Public :: qmass,ceng,pmass,chip0,rf,factor,temp_lang
    Real( Kind = wp ), Allocatable, Public :: dens0(:)
    Integer,   Public :: ntp,stp,rtp
    Real( Kind = wp ), Public :: rcell(1:9),sx,sy,sz,chit_sb = 0.0_wp
! q. index arrays and tp. sum arrays

    Integer, Allocatable, Public   :: qn(:),tpn(:)
    Integer, Allocatable, Public   :: qs(:,:),tps(:)
    Integer, Allocatable, Public   :: qr(:),tpr(:)
    
    Real( Kind = wp ), Public              :: fpl(1:9) = 0.0_wp
    Real( Kind = wp ), Public, Allocatable :: fxl(:),fyl(:),fzl(:)

    !> variable timestep control
    Logical, Public :: lvar
    Real( Kind = wp ), Public :: tstep
    Real( Kind = wp ), Public :: mndis,mxdis,mxstp

  Contains
    Private

    Procedure, Public :: init_dpd => allocate_dpd_arrays
    Final :: cleanup
  End Type thermostat_type

  ! Thermostat keys
  !> Microcannonical ensemble
  Integer( Kind = wi ), Parameter, Public :: ENS_NVE = 0

  !> Cannonical ensemble Evans
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_EVANS = 1
  !> Cannonical ensemble Langevin (stochastic dynamics)
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_LANGEVIN = 10
  !> Cannonical ensemble Anderson
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_ANDERSON = 11
  !> Cannonical ensemble Berendsen
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_BERENDSEN = 12
  !> Cannonical ensemble Nosé-Hoover
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_NOSE_HOOVER = 13
  !> Cannonical ensemble gentle stocastic
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_GENTLE = 14
  !> Cannonical ensemble inhomogeneous Langevin (stocastic dynamics)
  Integer( Kind = wi ), Parameter, Public :: ENS_NVT_LANGEVIN_INHOMO = 15

  !> Isobaric ensemble isothermal Langevin (stochastic dynamics) (isotropic)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_LANGEVIN = 20
  !> Isobaric isothermal ensemble Berendsen
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_BERENDSEN = 21
  !> Isobaric isothermal ensemble Nosé-Hoover (isotropic) (Melchionna)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_NOSE_HOOVER = 22
  !> Isobaric isothermal ensemble Martyna-Tuckerman-Klein (isotropic)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_MTK = 23

  !> Isobaric isothermal ensemble anisotropic Langvein (stochastic dynamics)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_LANGEVIN_ANISO = 30
  !> Isobaric isothermal ensemble anisotropic Berendsen
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_BERENDSEN_ANISO = 31
  !> Isobaric isothermal ensemble anisotropic Nosé-Hoover (Melchionna)
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_NOSE_HOOVER_ANISO = 32
  !> Isobaric isothermal ensemble anistropic Martyna-Tuckerman-Klein
  Integer( Kind = wi ), Parameter, Public :: ENS_NPT_MTK_ANISO = 33

  ! Anisotropic barostat constraint keys
  !> Fully anisotropic, no constraints
  Integer( Kind = wi ), Parameter, Public :: CONSTRAINT_NONE = 0
  !> Semi-isotropic, constant normal pressure and surface area
  Integer( Kind = wi ), Parameter, Public :: CONSTRAINT_SURFACE_AREA = 1
  !> Semi-isotropic, constant normal pressure and surface tension
  !> (Orthorhombic constraints when thermo%tension = 0)
  Integer( Kind = wi ), Parameter, Public :: CONSTRAINT_SURFACE_TENSION = 2
  !> Semi-orthorhombic constrains
  Integer( Kind = wi ), Parameter, Public :: CONSTRAINT_SEMI_ORTHORHOMBIC = 3

  Interface adjust_timestep
    Module Procedure adjust_timestep_1
    Module Procedure adjust_timestep_2
  End Interface

  Public :: adjust_timestep

Contains

  Subroutine allocate_dpd_arrays(thermo,max_vdw)
    Class( thermostat_type ) :: thermo
    Integer( Kind = wi ), Intent( In    ) :: max_vdw

    Integer :: fail

    If (thermo%key_dpd == 0) Return

    fail = 0

    Allocate (thermo%gamdpd(0:max_vdw),thermo%sigdpd(1:max_vdw), stat=fail)

    If (fail > 0) Call error(1081)

    thermo%gamdpd = 0.0_wp
    thermo%sigdpd = 0.0_wp
  End Subroutine allocate_dpd_arrays

  Subroutine cleanup(thermo)
    Type(thermostat_type) :: thermo

    If (Allocated(thermo%gamdpd)) Then
      Deallocate(thermo%gamdpd)
    End If
    If (Allocated(thermo%sigdpd)) Then
      Deallocate(thermo%sigdpd)
    End If
    If (Allocated(thermo%dens0)) Then
      Deallocate(thermo%dens0)
    End If
    If (Allocated(thermo%qn)) Then
      Deallocate(thermo%qn)
    End If
    If (Allocated(thermo%tpn)) Then
      Deallocate(thermo%tpn)
    End If
    If (Allocated(thermo%qs)) Then
      Deallocate(thermo%qs)
    End If
    If (Allocated(thermo%tps)) Then
      Deallocate(thermo%tps)
    End If
    If (Allocated(thermo%qr)) Then
      Deallocate(thermo%qr)
    End If
    If (Allocated(thermo%tpr)) Then
      Deallocate(thermo%tpr)
    End If
  End Subroutine cleanup

  Logical Function adjust_timestep_1(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,parts,&
    xxt,yyt,zzt,legshl,message,t,comm)
        ! update maximum distance a particle has travelled
        Real( Kind = wp ), Intent( InOut ) :: tstep
        Real( Kind = wp ), Intent(   Out ) :: hstep,rstep
        Type( corePart ),  Intent( In    ) :: parts(1:)
        Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp,&
                    xxt(1:),yyt(1:),zzt(1:)
        Integer,           Intent( In    ) :: legshl(0:,1:),natms          
        Character( Len = *), Intent(   Out ) :: message
        Real( Kind =  wp ), Intent( InOut) :: t
        Type( comms_type), Intent( InOut ) :: comm

        Real( Kind = wp ) :: mxdr
        Logical :: lv_up 
        Logical :: lv_dn
        Integer :: i



        adjust_timestep_1 = .False.
        lv_up = .False.
        lv_dn = .False.
        mxdr = 0.0_wp
        t=0.0_wp
        
        Do i=1,natms
          If (legshl(0,i) >= 0) &
            mxdr=Max(mxdr,(parts(i)%xxx-xxt(i))**2 + (parts(i)%yyy-yyt(i))**2 + (parts(i)%zzz-zzt(i))**2)
        End Do
        mxdr=Sqrt(mxdr)
        Call gmax(comm,mxdr)

        If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

          ! scale tstep and derivatives

          If (mxdr > mxdis) Then
            lv_up = .true.
            If (lv_dn) Then
              t = Sqrt(4.0_wp/3.0_wp)
              tstep = 0.75_wp*tstep
              hstep = 0.50_wp*tstep
            Else
              t = Sqrt(2.0_wp)
              tstep = hstep
              hstep = 0.50_wp*tstep
            End If
            Write(message,"( &
              & 'timestep decreased, new timestep is:',3x,1p,e12.4)") tstep
          End If
          If (mxdr < mndis) Then
            lv_dn = .true.
            If (lv_up) Then
              t = Sqrt(2.0_wp/3.0_wp)
              tstep = 1.50_wp*tstep
              hstep = 0.50_wp*tstep
            Else
              t = Sqrt(0.5_wp)
              hstep = tstep
              tstep = 2.00_wp*tstep
            End If
            If (tstep > mxstp) Then
              t = t*Sqrt(tstep/mxstp)
              tstep = mxstp
              hstep = 0.50_wp*tstep
            End If
            Write(message,"( &
              & 'timestep increased, new timestep is:',3x,1p,e12.4)") tstep
          End If
          rstep = 1.0_wp/tstep

          ! restart vv1

          adjust_timestep_1 = .True.
        End If
    End Function adjust_timestep_1
  Logical Function adjust_timestep_2(tstep,hstep,rstep,qstep,mndis,mxdis,mxstp,natms,parts,&
    xxt,yyt,zzt,legshl,message,t,comm)
        ! update maximum distance a particle has travelled
        Real( Kind = wp ), Intent( InOut ) :: tstep
        Real( Kind = wp ), Intent(   Out ) :: hstep,rstep,qstep
        Type( corePart ),  Intent( In    ) :: parts(1:)
        Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp,&
                    xxt(1:),yyt(1:),zzt(1:)
        Integer,           Intent( In    ) :: legshl(0:,1:),natms          
        Character( Len = *), Intent(   Out ) :: message
        Real( Kind =  wp ), Intent( InOut) :: t
        Type( comms_type), Intent( InOut ) :: comm

        Real( Kind = wp ) :: mxdr
        Logical :: lv_up 
        Logical :: lv_dn
        Integer :: i



        adjust_timestep_2 = .False.
        lv_up = .False.
        lv_dn = .False.
        mxdr = 0.0_wp
        
        Do i=1,natms
          If (legshl(0,i) >= 0) &
            mxdr=Max(mxdr,(parts(i)%xxx-xxt(i))**2 + (parts(i)%yyy-yyt(i))**2 + (parts(i)%zzz-zzt(i))**2)
        End Do
        mxdr=Sqrt(mxdr)
        Call gmax(comm,mxdr)

        If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

          ! scale tstep and derivatives
             If (mxdr > mxdis) Then
                lv_up = .true.
                If (lv_dn) Then
                   t = Sqrt(4.0_wp/3.0_wp)
                   tstep = 0.75_wp*tstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                Else
                   t = Sqrt(2.0_wp)
                   tstep = hstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                End If
                Write(message,'(a,1p,e12.4)') &
                  'timestep decreased, new timestep is: ', tstep
             End If
             If (mxdr < mndis) Then
                lv_dn = .true.
                If (lv_up) Then
                   t = Sqrt(2.0_wp/3.0_wp)
                   tstep = 1.50_wp*tstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                Else
                   t = Sqrt(0.5_wp)
                   qstep = hstep
                   hstep = tstep
                   tstep = 2.00_wp*tstep
                End If
                If (tstep > mxstp) Then
                   t = t*Sqrt(tstep/mxstp)
                   tstep = mxstp
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                End If
                Write(message,'(a,1p,e12.4)') &
                  'timestep increased, new timestep is: ', tstep
             End If
             rstep = 1.0_wp/tstep

          ! restart vv1

          adjust_timestep_2 = .True.
        End If
    End Function adjust_timestep_2
End Module thermostat
