Module kontrol
  Use kinds, only : wp,real32,real64
  Use comms,      Only : comms_type,gcheck
  Use timer,      Only : timer_type
  Use configuration,     Only : sysname
  Use mpole,     Only : mpole_type,POLARISATION_DEFAULT,POLARISATION_CHARMM
  Use langevin,   Only : langevin_allocate_arrays
  Use bonds,      Only : bonds_type
  Use angles,     Only : angles_type
  Use dihedrals,  Only : dihedrals_type
  Use inversions, Only : inversions_type
  Use vdw, Only : vdw_type,MIX_NULL,MIX_LORENTZ_BERTHELOT,MIX_FENDER_HASLEY, &
                  MIX_HALGREN,MIX_HOGERVORST,MIX_WALDMAN_HAGLER,MIX_TANG_TOENNIES, &
                  MIX_FUNCTIONAL
  Use tersoff, Only : tersoff_type
  Use metal,      Only : metal_type
  Use poisson,    Only : poisson_type
  Use msd,        Only : msd_type
  Use kinetics,  Only : l_vom
  Use plumed,   Only : plumed_type
  Use setup,       Only : nread,control,pi,zero_plus,seed, &
                          output,field,config,statis, &
                          history,historf,revive,revcon,revold
  Use parse,       Only : get_line,get_word,lower_case,word_2_real
  
  Use kim,         Only : kimim,rkim
  Use greenkubo,   Only : greenkubo_type
  Use rdfs,        Only : rdf_type
  Use development, Only : development_type
  Use ttm
  Use impacts,     Only : impact_type
  Use defects,     Only : defects_type
  
  Use io,     Only : io_set_parameters,        &
                            io_get_parameters,        &
                            io_nc_set_real_precision, &
                            io_nc_compiled,           &
                            IO_READ_MPIIO,            &
                            IO_READ_DIRECT,           &
                            IO_READ_MASTER,           &
                            IO_READ_NETCDF,           &
                            IO_WRITE_UNSORTED_MPIIO,  &
                            IO_WRITE_UNSORTED_DIRECT, &
                            IO_WRITE_UNSORTED_MASTER, &
                            IO_WRITE_SORTED_MPIIO,    &
                            IO_WRITE_SORTED_DIRECT,   &
                            IO_WRITE_SORTED_NETCDF,   &
                            IO_WRITE_SORTED_MASTER
  Use netcdf_wrap, Only : netcdf_param
  Use numerics, Only : dcell, invert
  Use thermostat, Only : thermostat_type, &
                         ENS_NVE, ENS_NVT_EVANS, ENS_NVT_LANGEVIN,  &
                         ENS_NVT_ANDERSON, ENS_NVT_BERENDSEN, ENS_NVT_NOSE_HOOVER, &
                         ENS_NVT_GENTLE, ENS_NVT_LANGEVIN_INHOMO, &
                         ENS_NPT_LANGEVIN, ENS_NPT_BERENDSEN, ENS_NPT_NOSE_HOOVER, &
                         ENS_NPT_MTK, ENS_NPT_LANGEVIN_ANISO, ENS_NPT_BERENDSEN_ANISO, &
                         ENS_NPT_NOSE_HOOVER_ANISO,ENS_NPT_MTK_ANISO
  Use statistics, Only : stats_type
  USe z_density, Only : z_density_type
  Use constraints, Only : constraints_type
  Use pmf, Only : pmf_type
  Use neighbours, Only : neighbours_type
  Use core_shell, Only : core_shell_type
  Use minimise, Only : minimise_type

  Implicit None
  Private
  Public :: read_control
  Public :: scan_control_output
  Public :: scan_control_io
  Public :: scan_control
  Public :: scan_control_pre
  
  Contains


Subroutine read_control                                &
           (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,        &
           rbin,nstfce,alpha,width,     &
           l_exp,lecx,lfcap,l_top,          &
           lvar,leql,               &
           lfce,lpana,           &
           ltraj,lrsd,               &
           nx,ny,nz,impa,                            &
           keyres,                   &
           tstep,mndis,mxdis,mxstp,nstrun,nsteql,      &
           fmax,nstbpo,keyfce,epsq,             &
           rlx_tol,mxquat,quattol,       &
           nstbnd,nstang,nstdih,nstinv,  &
           nstraj,istraj,keytrj,         &
           dfcts,nsrsd,isrsd,rrsd,          &
           ndump,pdplnc,cshell,cons,pmf,stats,thermo,green,devel,plume,msd_data,met, &
           pois,bond,angle,dihedral,inversion,zdensity,neigh,vdw,tersoff,rdf,minimise,mpole,tmr,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the simulation control parameters
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2017
! contrib   - i.j.bush february 2014
! contrib   - a.v.brukhno march 2014
! contrib   - m.a.seaton june 2014
! contrib   - h.a.boateng february 2015
! contrib   - p.s.petkov february 2015
! contrib   - a.m.elena september 2015
! contrib   - a.m.elena february 2017
! contrib   - g.khara & m.a.seaton march 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Logical,                Intent( In    ) :: l_str,lsim,l_vv,l_n_e,l_n_v
  Integer,                Intent( In    ) :: levcfg
  Integer,                Intent( InOut ) :: nstfce
  Real( Kind = wp ),      Intent( In    ) :: rbin,width
  Real( Kind = wp ),      Intent( InOut ) :: alpha

  Logical,                Intent(   Out ) :: l_exp,lecx,            &
                                             lfcap,l_top,           &
                                             lvar,leql,lfce,   &
                                             lpana,                 &
                                             ltraj,lrsd


  Integer,                Intent(   Out ) :: nx,ny,nz,             &
                                             keyres,nstrun,        &
                                             nsteql,       &
                                             nstbpo,        &
                                             keyfce,        &
                                             mxquat,        &
                                             nstbnd,nstang,        &
                                             nstdih,nstinv,        &
                                             nstraj,istraj,keytrj, &
                                             nsrsd,isrsd,          &
                                             ndump

  Real( Kind = wp ),      Intent(   Out ) :: tstep,mndis,mxdis,mxstp,    &
                                              quattol,&
                                             fmax,epsq,rlx_tol(1:2),     &
                                             rrsd,pdplnc
  Type( pmf_type ), Intent (   InOut )   :: pmf
  Type( core_shell_type ), Intent (   In  )   :: cshell
  Type( constraints_type ), Intent (   InOut )   :: cons
  Type( stats_type ), Intent (   InOut )   :: stats
  Type( impact_type ),     Intent(   Out ) :: impa
  Type ( thermostat_type), Intent( InOut ) :: thermo
  Type( development_type ), Intent( InOut ) :: devel
  Type( greenkubo_type ), Intent( InOut ) :: green
  Type( plumed_type ), Intent( InOut ) :: plume
  Type( msd_type ), Intent( InOut ) :: msd_data
  Type( metal_type ), Intent( InOut ) :: met
  Type( poisson_type ), Intent( InOut ) :: pois
  Type( bonds_type ), Intent( InOut ) :: bond
  Type( angles_type ), Intent( In    ) :: angle
  Type( dihedrals_type ), Intent( In    ) :: dihedral
  Type( inversions_type ), Intent( InOut ) :: inversion
  Type( z_density_type ), Intent( InOut ) :: zdensity
  Type( neighbours_type ), Intent( In    ) :: neigh
  Type( vdw_type ), Intent( InOut ) :: vdw
  Type( tersoff_type ), Intent( In    )  :: tersoff
  Type( rdf_type ), Intent( InOut ) :: rdf
  Type( minimise_type ), Intent( InOut ) :: minimise
  Type( mpole_type ), Intent( InOut ) :: mpole
  Type( timer_type ),      Intent( InOut ) :: tmr
  Type( defects_type ),    Intent( InOut ) :: dfcts(:)
  Type( comms_type ),     Intent( InOut )  :: comm


  Logical                                 :: limp,lvv,lens,lforc,     &
                                             ltemp,l_0,lpres,lstrext, &
                                             lstep,lplumed,safe,      &
                                             l_timjob,l_timcls

  Character( Len = 200 )                  :: record
  Character( Len = 40  )                  :: word,word1,word2,word3,akey

  Integer                                 :: i,j,k,itmp,nstana,grdana,grdbnd,grdang, &
                                             grddih,grdinv,nstall

  Real( Kind = wp )                       :: rcell(1:9),rcut1,rpad1,rvdw1,tmp,eps0,tol,rcb_d,prmps(1:4)

  Character( Len = 256 ) :: message,messages(7)
  Character( Len = 80 )  :: banner(9)


! initialise system control variables and their logical switches

! default expansion option

  l_exp = .false.
  nx    = 1
  ny    = 1
  nz    = 1

! defaults for direct evaluation, force-shifting of VDW interactions
! and type of mixing for undefined cross interaction of certain type
!
! vdw%l_direct = .false. ! (initialised in vdw)
! vdw%l_force_shift = .false. ! (initialised in vdw)
  vdw%mixing = MIX_NULL       ! (initialised in vdw)
!
! defaults for direct evaluation of metal interactions
!
! met%l_direct = .false. ! (initialised in metal_module)

! default impact option: option applied, particle index,
! timestep of impact, energy of impact, (3) direction of impact

  limp = .false.
  impa%imd  = 0
  impa%tmd  = -1
  impa%emd  = 0.0_wp
  impa%vmx  = 0.0_wp
  impa%vmy  = 0.0_wp
  impa%vmz  = 0.0_wp

! temperature & pressure (stress) switches and default values

  ltemp   = .false.
  lpres   = .false.
  lstrext = .false.
  thermo%temp    = 0.0_wp
  thermo%press   = 0.0_wp
  thermo%stress  = 0.0_wp

! default restart key (general)

  keyres = 0

! timestep switch and default value

  lstep = .false.
  tstep = 0.0_wp

! variable timestep switches and default minimum and maximum
! distances for variable timestep

  lvar  = .false.
  mndis = 0.03_wp
  mxdis = 0.10_wp
  mxstp = 0.0_wp

! total number of steps to run

  nstrun = 0

! number of steps for equilibration and default for exclusion of
! statistics collection during those steps

  nsteql = 0
  leql   = .true.

! switch for pseudo thermostat (not applied), type of scaling
! (default 0 where 0 - Langevin+direct, 1 - Langevin, 2 - gauss, 3 - direct )
! and minimum width of the thermostatted boundaries in Angs
! minimum temperature of the thermostat

  thermo%l_pseudo   = .false.
  thermo%key_pseudo = 0
  thermo%width_pseudo = 2.0_wp
  thermo%temp_pseudo = 1.0_wp

! default switch for conjugate gradient minimisation during equilibration

  minimise%minimise   = .false.
  minimise%key = -1
  minimise%freq = 0
  minimise%tolerance = 0.0_wp
  minimise%step_length = -1.0_wp

! default switch for regaussing temperature and default number of
! steps when to be applied

  thermo%l_tgaus  = .false.
  thermo%freq_tgaus = 0

! default switch for temperature scaling and default number of
! steps when to be applied

  thermo%l_tscale  = .false.
  thermo%freq_tscale = 0

! default switch for zero temperature optimisation and default number of
! steps when to be applied

  thermo%l_zero   = .false.
  thermo%freq_zero = 0
  l_0     = .false. ! T/=10K

! default integration type (VV), ensemble switch (not defined) and key

  lvv    = .true.
  lens   = .false.
  thermo%ensemble = ENS_NVE

! default thermostat and barostat friction time constants

  thermo%tau_t   = 0.0_wp ! thermostat relaxation time
  thermo%chi    = 0.0_wp ! Stochastic Dynamics (SD Langevin) thermostat friction
  thermo%chi_ep = 0.5_wp ! Inhomogeneous Stochastic Dynamics (SD Langevin) 
                  ! thermostat/electron-phonon friction
  thermo%chi_es = 0.0_wp ! Inhomogeneous Stochastic Dynamics (SD Langevin)
                  ! thermostat/electronic stopping friction
  thermo%soft   = 0.0_wp ! Softness for Andersen thermostat
  thermo%gama   = 0.0_wp ! Stochastic (Langevin) friction on a thermostat
  thermo%tau_p   = 0.0_wp ! barostat relaxation time
  thermo%tai    = 0.0_wp ! Stochastic Dynamics (SD Langevin) barostat friction
  thermo%iso    = 0      ! no semi-isotropic feature
  thermo%tension    = 0.0_wp ! surface tension

! default value for inhomogeneous Langevin thermostat/
! two-temperature model source term cutoff velocity

  thermo%vel_es2 = 50.0_wp

! default value for accounting extended coulombic exclusion
! (not accounted)

  lecx = .false.

! default switch for force capping and cap value

  lfcap = .false.
  fmax  = 1000.0_wp

! default switch for printing topology

  l_top = .true.

! default switch for removing COM momentum for ensembles
!
!  l_vom = .true. ! initialised in kinetics

! defaults for: force key = no electrostatics,
! specified force field = not yet, relative dielectric constant = 1

  keyfce = 0
  lforc  = .false.
  epsq   = 1.0_wp

! default maximum number of iterations and maximum tolerance
! for constraint algorithms

  cons%max_iter_shake = 250
  cons%tolerance = 1.0e-6_wp

! default maximum number of iterations and maximum tolerance
! for LFV quaternion integration algorithms

  mxquat =100
  quattol=1.0e-8_wp

! Default relaxed shell model tolerance and optional CGM step

  rlx_tol(1:2) = (/ 1.0_wp , -1.0_wp /)

! default switch for two-temperature model (TTM) calculations:
! already determined its use in scan_control but repeating
! here to output message

  l_ttm = .false.

! default values for (i) electronic specific heat capacities,
! (ii) thermal conductivity, (iii) thermal diffusivity,
! (iv) atomic density (converting specific heats to volumetric values)

  Ce0     = 1.0_wp
  sh_A    = 0.0_wp
  sh_B    = 0.0_wp
  Cemax   = 0.0_wp
  Tfermi  = 0.0_wp

  Ka0     = 0.0_wp

  Diff0   = 0.0_wp

  cellrho = 0.0_wp

! default initial stopping power to be deposited in
! electronic system (standard cascade) and laser 
! fluence and penetration depth

  dEdX    = 0.0_wp
  fluence = 0.0_wp
  pdepth  = 0.0_wp

! default values for (i) spatial, (ii) temporal
! energy deposition

  sdepoType = 0
  sig       = 1.0_wp
  sigmax    = 5

  tdepoType = 1
  tdepo     = 1.0e-3_wp
  tcdepo    = 5

! default boundary conditions for electronic temperature

  bcTypeE = 3

! default minimum number of atoms required per voxel cell
! to calculate ionic temperatures and one-way electron-phonon
! coupling in thermal diffusion and thermostat

  amin     = 1
  oneway   = .false.

! default values for time step frequencies to output
! (i) statistical data and (ii) electronic/ionic temperature
! grid values

  ttmstats = 0
  ttmtraj  = 0

! default value for time to start electron-phonon coupling

  ttmoffset = 0.0_wp

! default switch for dynamic cell density calculations

  ttmdyndens = .false.

! proceed normal simulation

  lfce = .false. ! don't recalculate forces based on history positions

! PLUMED read detector

  lplumed=.true.

! default switches for calculation and printing of intramolecular analysis

  nstana = 0 ; grdana = 0
  nstbnd = 0 ; grdbnd = 0 ; rcb_d  = 0.0_wp
  nstang = 0 ; grdang = 0
  nstdih = 0 ; grddih = 0
  nstinv = 0 ; grdinv = 0
  lpana  = .false.

! default switch for calculation of rdfs, default number of steps
! when to be collected and default switch for printing them

  rdf%l_collect   = .false.
  rdf%freq = 1
  rdf%l_print  = .false.

! default switch for calculation of z-density profile, default number of steps
! when to be collected and default switch for printing it

  zdensity%l_collect   = .false.
  zdensity%frequency = 1
  zdensity%l_print  = .false.

! default switches for calculation of velocity autocorrelation functions:
! time-averaging and printing

  green%l_average = .true.
  green%l_print  = .false.

! default for data printing interval

  nstbpo = 100

! default for statistics file interval

  stats%intsta = 100

! default switch for MSD outputing and defaults for
! (i) step to start at, (ii) every step after to be collected

  msd_data%start = 0
  msd_data%freq = 1

! default switch for trajectory outputting and defaults for
! (i) step to start at, (ii) every step after to be collected,
! (iii) level of information to output

  ltraj  = .false.
  nstraj = 0
  istraj = 1
  keytrj = 0

! default switch for defects outputting and defaults for
! (i) step to start at, (ii) every step after to be collected,
! (iii) default value for accepting an atom has defected

  dfcts(:)%ldef   =.false.
  dfcts(:)%nsdef  = 0
  dfcts(:)%isdef  = 1
  dfcts(:)%rdef   = Min(0.75_wp,neigh%cutoff/3.0_wp)
  dfcts(:)%newjob = .True.

! default switch for displacements outputting and defaults for
! (i) step to start at, (ii) every step after to be collected,
! (iii) cutoff value for accepting an atom has moved

  lrsd   =.false.
  nsrsd  = 0
  isrsd  = 1
  rrsd   = 0.15_wp

! default value for data dumping interval

  ndump = 1000

! default value for the particle density per link cell limit
! below which subcelling (decreasing link-cell dimensions) stops

  pdplnc = 50.0_wp

! default times for job execution and output dump

  tmr%job = 0.0_wp ; l_timjob=.false.
  tmr%clear_screen = 0.0_wp ; l_timcls=.false.

! major cutoff, padding and vdw cutoff defaults

  rcut1 = 0.0_wp
  rpad1 = 0.0_wp
  rvdw1 = 0.0_wp

! open the simulation control file

  If (comm%idnode == 0) Open(Unit=nread, File = Trim(control), Status = 'old')

! read simulation control name

  Call get_line(safe,nread,sysname,comm)
  If (.not.safe) Go To 1000
  Call strip_blanks(sysname)

  If (.not.safe) Go To 1000

  Write(banner(1),'(a)') ''
  Write(banner(2),'(a)') Repeat('*',80)
  Write(banner(3),'(a4,a72,a4)') '*** ', 'title:'//Repeat(' ',66), ' ***'
  Write(banner(4),'(a4,a72,a4)') '*** ', sysname, ' ***'
  Write(banner(5),'(a)') Repeat('*',80)
  Write(banner(6),'(a)') ''
  Call info(banner,6,.true.)

  Call info('simulation control parameters',.true.)

! read and process directives from CONTROL file

  Do

     Call get_line(safe,nread,record,comm)
     If (.not.safe) Go To 1000
     Call lower_case(record)
     Call get_word(record,word)

! record is commented out

     If (word(1:1) == '#' .or. word(1:1) == ' ') Then

! read DEVELOPMENT options

     Else If (word(1:5) == 'l_scr') Then
!        l_scr = .true. ! done in scan_development
        Call info('%%% OUTPUT redirected to the default output (screen) !!! %%%',.true.)
     Else If (word(1:6) == 'l_fast') Then
!        l_fast = .true. ! done in scan_development
        Call info('%%% speed up by avoiding global safety checks !!! %%%',.true.)
     Else If (word(1:5) == 'devel%l_eng') Then
        devel%l_eng = .true.
        Call info('%%% OUTPUT contains an extra last line with E_tot !!! %%%',.true.)
     Else If (word(1:6) == 'devel%l_rout') Then
        devel%l_rout = .true.
        Call info('%%% REVIVE writing in ASCII opted !!! %%%',.true.)
     Else If (word(1:5) == 'devel%l_rin') Then
        devel%l_rin = .true.
        Call info('%%% REVOLD reading in ASCII opted !!! %%%',.true.)
     Else If (word(1:5) == 'devel%l_org') Then
        devel%l_org = .true.
        devel%l_trm  = .true.

        Call info('%%% translate CONFIG along a vector into CFGORG after reading input & terminate !!! %%%',.true.)
        Call info('%%% vector and config level read as follows: %%%',.true.)

        Call get_word(record,word)
        devel%xorg = word_2_real(word)
        Call get_word(record,word)
        devel%yorg = word_2_real(word)
        Call get_word(record,word)
        devel%zorg = word_2_real(word)

        Call get_word(record,word)
        devel%lvcforg = Min( Int(Abs(word_2_real(word,0.0_wp))) , levcfg)

        Write(messages(1),'(a)') '%%%'
        Write(messages(2),'(a,3f10.3,a)') '%%% vector(x,y,x) ', devel%xorg, devel%yorg, devel%zorg, ' %%%'
        Write(messages(3),'(a,i0,a)') '%%% CFGORG level ', devel%lvcforg, ' %%%'
        Call info(messages,3,.true.)

     Else If (word(1:5) == 'devel%l_scl') Then
        Call info('%%% rescale CONFIG to CFGSCL, after reading input & terminate !!! %%%',.true.)
        Call info('%%% config level and new cell vectors to rescale to (read in a CONFIG-like manner): %%%',.true.)

        Call get_word(record,word)
        devel%lvcfscl = Min( Int(Abs(word_2_real(word,0.0_wp))) , levcfg)

        itmp=0
        Do i=1,3
           Call get_line(safe,nread,record,comm)
           Do j=1,3
              Call get_word(record,word)
              itmp=itmp+1
              devel%cels(itmp)=word_2_real(word)
           End Do
        End Do

        Call invert(devel%cels,rcell,tmp)

        Write(messages(1),'(a)') '%%%'
        Write(messages(2),'(1x,a,i0,a)')        '%%% CFGSCL level ', devel%lvcfscl, ' %%%'
        Write(messages(3),'(1x,a,3f20.10,a)')   '%%% ', devel%cels(1:3), ' %%%'
        Write(messages(4),'(1x,a,3f20.10,a)')   '%%% ', devel%cels(4:6), ' %%%'
        Write(messages(5),'(1x,a,3f20.10,a)')   '%%% ', devel%cels(7:9), ' %%%'
        Write(messages(6),'(1x,a)')             '%%%'
        Write(messages(7),'(1x,a,1p,g22.12,a)') '%%% CFGSCL volume ', tmp, '%%%'
        Call info(messages,7,.true.)

        If (tmp > zero_plus) Then
           devel%l_scl = .true.
           devel%l_trm  = .true.
        Else
           Call info('%%% OPTION ABORTED DUE TO ZERO VOLUME !!! %%%',.true.)
           devel%l_trm  = .true.
        End If
     Else If (word(1:5) == 'devel%l_his') Then
        devel%l_his = .true.
        devel%l_trm = .true.
        Call info('%%% generate HISTORY after reading input & terminate !!! %%%',.true.)
     Else If (word(1:5) == 'l_tim') Then
!        l_tim = .true.  ! done in scan_development
        Call info('%%% generate detailed timing !!! %%%',.true.)
     Else If (word(1:5) == 'devel%l_tor') Then
        devel%l_tor = .true.
        Call info('%%% Turn off production of REVCON & REVIVE !!! %%%',.true.)
     Else If (word(1:5) == 'devel%l_trm') Then
        devel%l_trm = .true.
        Call info('%%% Terminate gracefully before initialisation !!! %%%',.true.)
     Else If (word(1:5) == 'devel%l_dis') Then
        devel%l_dis = .true.
        devel%r_dis = Min( devel%r_dis , word_2_real(word,0.1_wp) )
        Call info('%%% Turn on the check on minimum separation distance between VNL pairs at re/start !!! %%%',.true.)
        Write(message,'(a,1p,e12.4)') '%%% separation criterion (Angstroms) %%%', devel%r_dis
        Call info(message,.true.)

! read VDW options

     Else If (word(1:3) == 'vdw') Then
        Call get_word(record,word1)

        If      (word1(1:6) == 'direct') Then

! direct evaluation option

           vdw%l_direct = .true.
           Call info('vdw direct option on',.true.)

        Else If (word1(1:6) == 'mixing') Then

! mixing type keywords

           Call info('vdw cross terms mixing opted (for undefined mixed potentials)',.true.)
           Call info('mixing is limited to potentials of the same type only',.true.)
           Call info('mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)',.true.)

           Call get_word(record,word2)

           If      (word2(1:4) == 'lore') Then
              vdw%mixing = MIX_LORENTZ_BERTHELOT
              Call info('type of mixing selected - Lorentz–Berthelot :: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2',.true.)
           Else If (word2(1:4) == 'fend') Then
              vdw%mixing = MIX_FENDER_HASLEY
              Call info('type of mixing selected - Fender-Halsey :: e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2',.true.)
           Else If (word2(1:4) == 'hoge') Then
              vdw%mixing = MIX_HOGERVORST
              Call info('type of mixing selected - Hogervorst (good hope) :: ' &
                //'e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)',.true.)
           Else If (word2(1:4) == 'halg') Then
              vdw%mixing = MIX_HALGREN
              Call info('type of mixing selected - Halgren HHG :: ' &
                //'e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)',.true.)
           Else If (word2(1:4) == 'wald') Then
              vdw%mixing = MIX_WALDMAN_HAGLER
              Call info('type of mixing selected - Waldman–Hagler :: ' &
                //'e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ;s_ij=[(s_i^6+s_j^6)/2]^(1/6)',.true.)
           Else If (word2(1:4) == 'tang') Then
              vdw%mixing = MIX_TANG_TOENNIES
              Call info('type of mixing selected - Tang-Toennies :: ' &
                //' e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13',.true.)
              Call info(Repeat(' ',43)//'s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)',.true.)
           Else If (word2(1:4) == 'func') Then
              vdw%mixing = MIX_FUNCTIONAL
              Call info('type of mixing selected - Functional :: ' &
                //'e_ij=3 * (e_i*e_j)^(1/2) * (s_i*s_j)^3 / SUM_L=0^2{[(s_i^3+s_j^3)^2 / (4*(s_i*s_j)^L)]^(6/(6-2L))}',.true.)
              Call info(Repeat(' ',40)//'s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}',.true.)
           Else
              Call strip_blanks(record)
              Write(message,'(4a)') word(1:Len_Trim(word)+1), &
                word1(1:Len_Trim(word1)+1),word2(1:Len_Trim(word2)+1),record
              Call info(message,.true.)
              Call error(3)
           End If

        Else If (word1(1:5) == 'shift') Then
! force-shifting option
           vdw%l_force_shift = .true.
           Call info('vdw force-shifting option on',.true.)
        Else
           Call strip_blanks(record)
           Write(message,'(3a)') word(1:Len_Trim(word)+1), &
             word1(1:Len_Trim(word1)+1),record
           Call info(message,.true.)
           Call error(3)
        End If

     Else If (word(1:5) == 'metal') Then

        Call get_word(record,word)
        If      (word(1:6) == 'direct') Then
! read metal direct evaluation option
           Call info('metal direct option on',.true.)
           If (met%tab > 0) Then
              Call warning(480,0.0_wp,0.0_wp,0.0_wp)
           Else
              met%l_direct = .true.
           End If
        Else If (word(1:7) == 'sqrtrho') Then
! read metal sqrtrho interpolation option for EAM embeding function in TABEAM
           Call info('metal sqrtrho option on',.true.)
           If (met%tab > 0) Then
              met%l_emb = .true.
           Else
              Call warning(490,0.0_wp,0.0_wp,0.0_wp)
           End If
        End If

! read slab option (dealt with in scan_control<-set_bounds,
! affecting map_domains<-set_bounds)

     Else If (word(1:4) == 'slab') Then

        Call info('slab option on',.true.)

! io options (dealt with in scan_control<-set_bounds)

     Else If (word(1:2) == 'io' ) Then

! read expansion option

     Else If (word(1:5) == 'nfold') Then

        l_exp = .true.
        Call get_word(record,word)
        nx = Max(1,Nint(Abs(word_2_real(word))))
        Call get_word(record,word)
        ny = Max(1,Nint(Abs(word_2_real(word))))
        Call get_word(record,word)
        nz = Max(1,Nint(Abs(word_2_real(word))))
        Write(message,'(a,9x,3i5)') 'system expansion opted',nx,ny,nz
        Call info(message,.true.)

! read impact option

     Else If (word(1:6) == 'impact') Then

        Call get_word(record,word)
        impa%imd = Max(1,Nint(Abs(word_2_real(word))))
        Call get_word(record,word)
        impa%tmd = Nint(Abs(word_2_real(word)))

        Call get_word(record,word)
        impa%emd = Abs(word_2_real(word))
        Call get_word(record,word)
        impa%vmx = word_2_real(word)
        Call get_word(record,word)
        impa%vmy = word_2_real(word)
        Call get_word(record,word)
        impa%vmz = word_2_real(word)

        If (Sqrt(impa%vmx**2+impa%vmy**2+impa%vmz**2) <= zero_plus) Then
           impa%vmx = 1.0_wp
           impa%vmy = 1.0_wp
           impa%vmz = 1.0_wp
        End If

        Write(messages(1),'(a)') ''
        Write(messages(2),'(a,i10)') 'particle (index)',impa%imd
        Write(messages(3),'(a,i10)') 'timestep (steps)',impa%tmd
        Write(messages(4),'(a,1p,e12.4)') 'energy   (keV)  ',impa%emd
        Write(messages(5),'(a,1p,3e12.4)') 'v-r(x,y,z)      ',impa%vmx,impa%vmy,impa%vmz
        Call info(messages,5,.true.)

        If (limp) Call error(600)
        limp = .true.

! read seeding option

     Else If (word(1:4) == 'seed') Then

        lseed=.true.

        Call get_word(record,word)
        seed(1)=Nint(Abs(word_2_real(word)))
        Call get_word(record,word)
        seed(2)=Nint(Abs(word_2_real(word)))
        Call get_word(record,word)
        seed(3)=Nint(Abs(word_2_real(word)))

        Write(message,'(a,3i5)') 'radomisation seeds supplied: ',seed
        Call info(message,.true.)

! read temperature

     Else If (word(1:4) == 'temp') Then

        ltemp = .true.
        Call get_word(record,word)
        thermo%temp = Abs(word_2_real(word))
        Write(message,'(a,1p,e12.4)') 'simulation temperature (K)  ',thermo%temp
        Call info(message,.true.)

! read zero temperature optimisation

     Else If (word(1:4) == 'zero') Then

        thermo%l_zero = .true.

! Check defaults

        Call get_word(record,word)
        l_0 = (word(1:4) == 'fire')
        thermo%freq_zero = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        If (word(1:5) == 'every') Call get_word(record,word)
        thermo%freq_zero = Max(thermo%freq_zero,Abs(Nint(word_2_real(word,0.0_wp))))

        Call info('zero K optimisation on (during equilibration)',.true.)
        Write(message,'(a,i10)') 'temperature regaussing interval',thermo%freq_zero

        If (l_0) Then
           If (comm%idnode == 0) &
           Call info('fire option on - actual temperature will reset to 10 Kelvin if no target tempreature is specified',.true.)
        Else
           ltemp  = .true.
           thermo%temp = 10.0_wp
           Call info('fire option off - actual temperature reset to 10 Kelvin',.true.)
        End If

! read pressure

     Else If (word(1:4) == 'pres') Then

        Call get_word(record,word)

! read stress (6 components only: xx,yy,zz,xy,xz,yz - it's forced symmetric)

        If (word(1:6) == 'tensor') Then
           lstrext=.true.

           Call get_word(record,word)
           thermo%stress(1) = word_2_real(word)
           Call get_word(record,word)
           thermo%stress(5) = word_2_real(word)
           Call get_word(record,word)
           thermo%stress(9) = word_2_real(word)
           Call get_word(record,word)
           thermo%stress(2) = word_2_real(word)
           thermo%stress(4) = thermo%stress(2)
           Call get_word(record,word)
           thermo%stress(3) = word_2_real(word)
           thermo%stress(7) = thermo%stress(3)
           Call get_word(record,word)
           thermo%stress(6) = word_2_real(word)
           thermo%stress(8) = thermo%stress(6)

           Call info('simulation pressure tensor (katms):',.true.)
           Write(messages(1),'(3f20.10)') thermo%stress(1:3)
           Write(messages(2),'(3f20.10)') thermo%stress(4:6)
           Write(messages(3),'(3f20.10)') thermo%stress(7:9)
           Call info(messages,3,.true.)

! convert from katms to internal units of pressure

           thermo%stress = thermo%stress/prsunt
        Else
           lpres=.true.

           thermo%press = word_2_real(word)

           Write(message,'(a,1p,e12.4)') 'simulation pressure (katms) ',thermo%press

! convert from katms to internal units of pressure

           thermo%press = thermo%press/prsunt
        End If

! read restart

     Else If (word(1:7) == 'restart') Then

        Call get_word(record,word)

        If (word(1:7) == 'noscale' .or. word(1:7) == 'unscale') Then
           keyres = 3
           Call info('unscaled restart requested (starting a new simulation)',.true.)
        Else If (word(1:5) == 'scale') Then
           keyres = 2
           Call info('scaled restart requested (starting a new simulation)',.true.)
        Else
           keyres = 1
           Call info('restart requested (continuing an old simulation)',.true.)
        End If

! read timestep options

     Else If (word(1:8) == 'timestep' .or. word(1:8) == 'variable') Then

        lstep = .true.
        Call get_word(record,word1)

        If      (word(1:8) == 'timestep' .and. word1(1:8) /= 'variable') Then
           tstep = word_2_real(word1)
        Else If ( (word(1:8) == 'timestep' .and. word1(1:8) == 'variable') .or. &
                  (word(1:8) == 'variable' .and. word1(1:8) == 'timestep') ) Then
           lvar = .true.
           Call get_word(record,word)
           tstep = word_2_real(word)
        Else
           Call strip_blanks(record)
           Write(message,'(3a)') word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call info(message,.true.)
           Call error(3)
        End If

! read minimum and maximum distance for variable timestep

     Else If (word(1:6) == 'mindis' .or. &
              word(1:6) == 'maxdis' .or. &
              word(1:6) == 'mxstep') Then

        If (word(1:6) == 'mindis') Then
           Call get_word(record,word)
           mndis=Abs(word_2_real(word))
        End If
        If (word(1:6) == 'maxdis') Then
           Call get_word(record,word)
           mxdis=Abs(word_2_real(word))
        End If
        If (word(1:6) == 'mxstep') Then
           Call get_word(record,word)
           mxstp=Abs(word_2_real(word))
        End If

! read number of timesteps

     Else If (word(1:5) == 'steps') Then
        Call get_word(record,word)
        nstrun = Nint(word_2_real(word))
        Write(message,'(a,i10)') 'selected number of timesteps ',nstrun
        Call info(message,.true.)
! read number of equilibration timesteps
     Else If (word(1:5) == 'equil') Then
        Call get_word(record,word)
        If (word(1:5) == 'steps') Call get_word(record,word)
        nsteql = Abs(Nint(word_2_real(word)))
        Write(message,'(a,i10)') 'equilibration period (steps) ', nsteql
        Call info(message,.true.)
! read collection option
     Else If (word(1:7) == 'collect') Then
        leql = .false.
        Call info('equilibration included in overall averages',.true.)
! read pseudo thermostat option
     Else If (word(1:6) == 'pseudo') Then

        Call get_word(record,word)
        If      (word(1:4) == 'lang'  ) Then
           thermo%key_pseudo = 1
           Call get_word(record,word)
        Else If (word(1:5) == 'gauss') Then
           thermo%key_pseudo = 2
           Call get_word(record,word)
        Else If (word(1:6) == 'direct') Then
           thermo%key_pseudo = 3
           Call get_word(record,word)
        End If

! thermo%width_pseudo = 2 Angs by default

        tmp = Abs(word_2_real(word))
        If (width/4.0_wp > thermo%width_pseudo) Then
           thermo%l_pseudo = .true.
           If (comm%idnode == 0) Then
              Call info('pseudo thermostat attached to MD cell boundary',.true.)
              If      (thermo%key_pseudo == 0) Then
                Call info('thermostat control: Langevin + direct temperature scaling',.true.)
              Else If (thermo%key_pseudo == 1) Then
                Call info('thermostat control: Langevin temperature scaling',.true.)
              Else If (thermo%key_pseudo == 2) Then
                Call info('thermostat control: gaussian temperature scaling',.true.)
              Else If (thermo%key_pseudo == 3) Then
                Call info('thermostat control: direct temperature scaling',.true.)
              End If
              Write(message,'(a,1p,e12.4)') 'thermostat thickness (Angs) ',tmp
           End If

           If (width/4.0_wp > tmp .and. tmp >= thermo%width_pseudo) Then
              thermo%width_pseudo = tmp
           Else
              Call info('thermostat thickness insufficient - reset to 2 Angs',.true.)
           End If
        Else
           Call warning(280,thermo%width_pseudo,width,0.0_wp)
           Call error(530)
        End If

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (tmp <= zero_plus) Then
           thermo%temp_pseudo = thermo%temp
        Else
           thermo%temp_pseudo = Max(thermo%temp_pseudo,tmp)
        End If
        Write(message,'(a,1p,e12.4)') 'thermostat temperature (K) ',thermo%temp_pseudo
        Call info(message,.true.)

! read minimiser option

     Else If (word(1:5) == 'minim' .or. word(1:5) == 'optim') Then

        minimise%minimise=.true.
        word2=' ' ; word2=word
        Call get_word(record,word)

        If      (word(1:4) == 'forc') Then
           minimise%key=0
           word1='force   '
        Else If (word(1:4) == 'ener') Then
           minimise%key=1
           word1='energy  '
        Else If (word(1:4) == 'dist') Then
           minimise%key=2
           word1='distance'
        Else
           Call strip_blanks(record)
           Write(message,'(4a)') word2(1:Len_Trim(word2)+1),' ',word(1:Len_Trim(word)+1),record
           Call info(message,.true.)
           Call error(590)
        End If

        If (word2(1:5) == 'minim') Then
           Call get_word(record,word)
           minimise%freq = Abs(Nint(word_2_real(word,0.0_wp)))
        End If

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))

        itmp=0
        If      (minimise%key == 0) Then
           If (tmp < 1.0_wp .or. tmp > 1000.0_wp) Then
              minimise%tolerance=50.0_wp
              itmp=1
           Else
              minimise%tolerance=tmp
           End If
        Else If (minimise%key == 1) Then
           If (tmp < zero_plus .or. tmp > 0.01_wp) Then
              minimise%tolerance=0.005_wp
              itmp=1
           Else
             minimise%tolerance = tmp
             minimise%step_length = tmp
           End If
        Else If (minimise%key == 2) Then
           If (tmp < 1.0e-6_wp .or. tmp > 0.1_wp) Then
              minimise%tolerance=0.005_wp
              itmp=1
           Else
              minimise%tolerance=tmp
           End If
        End If

        If (itmp == 1) Call warning(360,tmp,minimise%tolerance,0.0_wp)

        Call get_word(record,word3)
        minimise%step_length = word_2_real(word3,-1.0_wp)

        If (word2(1:5) == 'minim') Then
           Write(messages(1),'(a)') 'minimisation option on (during equilibration)'
           Write(messages(2),'(a,a8)') 'minimisation criterion        ',word1(1:8)
           Write(messages(3),'(a,i10)') 'minimisation frequency (steps)',minimise%freq
           Write(messages(4),'(a,1p,e12.4)') 'minimisation tolerance        ',minimise%tolerance
           Call info(messages,4,.true.)
           If (minimise%step_length > zero_plus) Then
             Write(message,'(a,1p,e12.4)') 'minimisation CGM step         ',minimise%step_length
             Call info(message,.true.)
           End If
        Else
           Write(messages(1),'(a)') 'optimisation at start'
           Write(messages(2),'(a,a8)') 'optimisation criterion        ',word(1:8)
           Write(messages(4),'(a,1p,e12.4)') 'optimisation tolerance        ',minimise%tolerance
           Call info(messages,3,.true.)
           If (minimise%step_length > zero_plus) Then
             Write(message,'(a,1p,e12.4)')'optimisation CGM step         ',minimise%step_length
           End If
        End If

! read regauss option

     Else If (word(1:6) == 'regaus') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'thermo%temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'thermo%temp') Call get_word(record,word)
        thermo%freq_tgaus = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        thermo%l_tgaus =.true.
        Write(message,'(a,i10)') 'temperature regaussing interval ', thermo%freq_tgaus
        Call info(message,.true.)

! read temperature scaling option

     Else If (word(1:5) == 'scale') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'thermo%temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'thermo%temp') Call get_word(record,word)
        thermo%freq_tscale = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        thermo%l_tscale =.true.
        Call info('temperature scaling on (during equilibration)',.true.)
        Write(message,'(a,i10)') 'temperature scaling interval ',thermo%freq_tscale
        Call info(message,.true.)

! read polarisation option

     Else If (word(1:5) == 'polar') Then

        Call get_word(record,word)
        If (word(1:6) == 'scheme' .or. word(1:4) == 'type') Call get_word(record,word)
        If (word(1:6) == 'scheme' .or. word(1:4) == 'type') Call get_word(record,word)
        If (word(1:6) == 'charmm') Then
           If (word(1:6) == 'thole') Then
              Call get_word(record,word)
              If (word(1:4) == 'dump' .or. word(1:6) == 'factor') Call get_word(record,word)
              If (word(1:4) == 'dump' .or. word(1:6) == 'factor') Call get_word(record,word)
              mpole%thole = Abs(word_2_real(word,0.0_wp))
           End If

           Write(message,'(a,f5.2)') &
             'CHARMM polarisation scheme selected with optional atomic thole dumping of ', &
             mpole%thole
           Call info(message,.true.)
           If (mpole%max_mpoles == 0 ) Then
             Call warning('scheme deselected due to switched off electrostatics',.true.)
           End If
           If (cshell%mxshl == 0) Then
             Call warning('scheme disabled due to lack of core-shell defined interatcions',.true.)
           End If

           If (mpole%max_mpoles == 0 .or. cshell%mxshl == 0) Then
!              mpole%key=POLARISATION_DEFAULT ! done in scan_control
           Else
              lecx = .true. ! enable extended coulombic exclusion
              Call info('Extended Coulombic eXclusion activated for CHARMM polarisation',.true.)
           End If
        End If

! read integration flavour

     Else If (word(1:8) == 'integrat') Then

        Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:8) == 'leapfrog') lvv=.false.

! thermo%key_dpd detected in scan_control

        If (thermo%key_dpd > 0 .and. (lvv .neqv. l_vv)) Then
          Call warning('Leapfrog Verlet selected integration defaulted to Velocity Verlet for DPD thermostats',.true.)
        End If

! read ensemble

     Else If (word(1:8) == 'ensemble') Then

        If (l_vv) Then
          Call info('Integration : Velocity Verlet',.true.)
        Else
          Call info('Integration : Leapfrog Verlet',.true.)
        End If

        Call get_word(record,word)

        If      (word(1:3) == 'nve' .or. word(1:3) == 'pmf') Then

           thermo%ensemble = ENS_NVE

           Call info('Ensemble : NVE (Microcanonical)',.true.)

           If (lens) Call error(414)
           lens=.true.

        Else If (word(1:3) == 'nvt') Then

           Call get_word(record,word)

           If      (word(1:5) == 'evans') Then

              thermo%ensemble = ENS_NVT_EVANS

              Call info('Ensemble : NVT Evans (Isokinetic)',.true.)
              Call info('Gaussian temperature constraints in use',.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:4) == 'lang') Then

              thermo%ensemble = ENS_NVT_LANGEVIN
              If (.not.l_vv) thermo%l_langevin = .true.

              Call get_word(record,word)
              thermo%chi = Abs(word_2_real(word))

              Call info('Ensemble : NVT Langevin (Stochastic Dynamics)',.true.)
              Write(message,'(a,1p,e12.4)') 'thermostat friction ', thermo%chi
              Call info(message,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:5) == 'ander') Then

              thermo%ensemble = ENS_NVT_ANDERSON

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%soft = Abs(word_2_real(word))
              If (thermo%soft > 1.0_wp) thermo%soft=1.0_wp/thermo%soft

              Write(messages(1),'(a)') 'Ensemble : NVT Andersen'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'softness (dimensionless)',thermo%soft
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              thermo%ensemble = ENS_NVT_BERENDSEN

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))

              Call info('Ensemble : NVT Berendsen',.true.)
              Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Call info(message,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              thermo%ensemble = ENS_NVT_NOSE_HOOVER

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))

              Call info('Ensemble : NVT Nose-Hoover',.true.)
              Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Call info(message,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'gst') Then

              thermo%ensemble = ENS_NVT_GENTLE

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))

              Call get_word(record,word)
              thermo%gama = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NVT gentle stochastic thermostat'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'friction on thermostat  (ps^-1) ',thermo%gama
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:4) == 'ttm' .or. word(1:6) == 'inhomo') Then

              thermo%ensemble = ENS_NVT_LANGEVIN_INHOMO
              If (.not.l_vv) thermo%l_langevin = .true.

              Call get_word(record,word)
              thermo%chi_ep  = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%chi_es  = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%vel_es2 = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)'
              Write(messages(2),'(a,1p,e12.4)') 'e-phonon friction (ps^-1) ',thermo%chi_ep
              Write(messages(3),'(a,1p,e12.4)') 'e-stopping friction (ps^-1) ',thermo%chi_es
              Write(messages(4),'(a,1p,e12.4)') 'e-stopping velocity (A ps^-1) ',thermo%vel_es2
              Call info(messages,4,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'dpd') Then

              Call info('Ensemble : NVT dpd (Dissipative Particle Dynamics)',.true.)

! thermo%key_dpd determined in scan_control

              If      (thermo%key_dpd == 1) Then
                 thermo%ensemble = ENS_NVE ! equivalence to doing NVE with some extra fiddling before VV(0)
                 Call info("Ensemble type : Shardlow's first order splitting (S1)",.true.)
              Else If (thermo%key_dpd == 2) Then
                 thermo%ensemble = ENS_NVE ! equivalence to doing NVE with some extra fiddling before VV(0) and after VV(1)
                 Call info("Ensemble type : Shardlow's second order splitting (S2)",.true.)
              Else
                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call error(436)
              End If

              Call get_word(record,word)
              thermo%gamdpd(0) = Abs(word_2_real(word,0.0_wp))

              If (thermo%gamdpd(0) > zero_plus) Then
                 Write(message,'(a,1p,e12.4)') 'drag coefficient (Dalton/ps) ', thermo%gamdpd(0)
                 Call info(message,.true.)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else

              Call strip_blanks(record)
              Write(message,'(2a)') word(1:Len_Trim(word)+1),record
              Call info(message,.true.)
              Call error(436)

           End If

        Else If (word(1:3) == 'npt') Then

           thermo%variable_cell = .true.

           Call get_word(record,word)

           If (word(1:4) == 'lang') Then

              thermo%ensemble = ENS_NPT_LANGEVIN
              thermo%l_langevin = .true.

              Call get_word(record,word)
              thermo%chi = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tai = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Langevin (Stochastic Dynamics)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat friction (ps^-1)',thermo%chi
              Write(messages(3),'(a,1p,e12.4)') 'barostat friction (ps^-1)',thermo%tai
              Call info(messages,3,.true.)

!                 thermo%tau_t=1/(2.0_wp*pi*thermo%chi)
!                 thermo%tau_p=1/(2.0_wp*pi*thermo%tai)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              thermo%ensemble = ENS_NPT_BERENDSEN

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tau_p = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Berendsen'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              thermo%ensemble = ENS_NPT_NOSE_HOOVER

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tau_p = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Nose-Hoover (Melchionna)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'mtk') Then

              thermo%ensemble = ENS_NPT_MTK

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tau_p = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Martyna-Tuckerman-Klein'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else

              Call strip_blanks(record)
              Write(message,'(2a)') word(1:Len_Trim(word)+1),record
              Call info(message,.true.)
              Call error(436)

           End If

        Else If (word(1:3) == 'nst') Then

           thermo%variable_cell = .true.
           thermo%anisotropic_pressure = .true.

           Call get_word(record,word)

           If (word(1:4) == 'lang') Then

              thermo%ensemble = ENS_NPT_LANGEVIN_ANISO
              thermo%l_langevin = .true.

              Call get_word(record,word)
              thermo%chi = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tai = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Langevin (Stochastic Dynamics)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat friction (ps^-1)',thermo%chi
              Write(messages(3),'(a,1p,e12.4)') 'barostat friction (ps^-1)',thermo%tai
              Call info(messages,3,.true.)

!                 thermo%tau_t=thermo%chi
!                 thermo%tau_p=2.0_wp*pi/thermo%tai

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 thermo%iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 thermo%iso=2
                 Call get_word(record,word)
                 thermo%tension = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', thermo%tension
                 Call info(messages,3,.true.)
                 thermo%tension=thermo%tension/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    thermo%iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (thermo%iso >= 1 .and. thermo%iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis',.true.)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              thermo%ensemble = ENS_NPT_BERENDSEN_ANISO

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tau_p = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Berendsen'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
              Call info(messages,3,.true.)

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 thermo%iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 thermo%iso=2
                 Call get_word(record,word)
                 thermo%tension = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', thermo%tension
                 Call info(messages,3,.true.)
                 thermo%tension=thermo%tension/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    thermo%iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (thermo%iso >= 1 .and. thermo%iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis')
              End If

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              thermo%ensemble = ENS_NPT_NOSE_HOOVER_ANISO

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tau_p = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Nose-Hoover (Melchionna)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
              Call info(messages,3,.true.)

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 thermo%iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 thermo%iso=2
                 Call get_word(record,word)
                 thermo%tension = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', thermo%tension
                 Call info(messages,3,.true.)
                 thermo%tension=thermo%tension/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    thermo%iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (thermo%iso >= 1 .and. thermo%iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis',.true.)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'mtk') Then

              thermo%ensemble = ENS_NPT_MTK_ANISO

              Call get_word(record,word)
              thermo%tau_t = Abs(word_2_real(word))
              Call get_word(record,word)
              thermo%tau_p = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Martyna-Tuckerman-Klein'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
              Call info(messages,3,.true.)

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 thermo%iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 thermo%iso=2
                 Call get_word(record,word)
                 thermo%tension = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', thermo%tension
                 Call info(messages,3,.true.)
                 thermo%tension=thermo%tension/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    thermo%iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    thermo%iso=3
                    Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
                 Else
                    Call strip_blanks(record)
                    Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                    Call info(message,.true.)
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (thermo%iso >= 1 .and. thermo%iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis',.true.)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else

              Call strip_blanks(record)
              Write(message,'(2a)') word(1:Len_Trim(word)+1),record
              Call info(message,.true.)
              Call error(436)

           End If

        Else

           Call strip_blanks(record)
           Write(message,'(2a)') word(1:Len_Trim(word)+1),record
           Call info(message,.true.)
           Call error(436)

        End If

! For Langevin ensembles that require arrays

        If (thermo%l_langevin) Call langevin_allocate_arrays()

! read density variation option

     Else If (word(1:7) == 'densvar') Then

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))

        Write(message,'(a,1p,e12.4)') 'density variation allowance (%) ', tmp
        Call info(message,.true.)

! read real space cutoff

     Else If (word(1:3) == 'cut' .or. word(1:4) == 'rcut') Then


        Call get_word(record,word)
        rcut1 = Abs(word_2_real(word))

        Write(message,'(a,1p,e12.4)') 'real space cutoff (Angs) ', rcut1
        Call info(message,.true.)

! read real space cutoff padding

     Else If (word(1:3) == 'pad' .or. word(1:4) == 'rpad') Then

        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        rpad1 = Abs(word_2_real(word))
        Write(message,'(a,1p,e12.4)') 'cutoff padding (Angs) ', rpad1
        Call info(message,.true.)

! read vdw cutoff (short-range potentials)

     Else If (word(1:4) == 'rvdw') Then

        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        rvdw1 = Abs(word_2_real(word))

        Write(message,'(a,1p,e12.4)') 'vdw cutoff (Angs) ', rvdw1
        Call info(message,.true.)

! read Ewald sum parameters

     Else If (word(1:5) == 'ewald' .or. word(1:4) == 'spme') Then

        Call get_word(record,word)

        If (word(1:5) == 'evalu') Then

! This is sorted in set_bounds -> scan_control

        Else

           keyfce = 2

           Call info('Electrostatics : Smooth Particle Mesh Ewald',.true.)

           If (word(1:9) == 'precision') Then
              Call get_word(record,word)
              tmp = Abs(word_2_real(word))
              Write(message,'(a,1p,e12.4)') 'Ewald sum precision ',tmp
              Call info(message,.true.)
           End If

! This is sorted in set_bounds -> scan_control

           Write(messages(1),'(a,1p,e12.4)') 'Ewald convergence parameter (A^-1) ',alpha
           Write(messages(2),'(a,3i5)') 'Ewald kmax1 kmax2 kmax3   (x2) ',kmaxa1,kmaxb1,kmaxc1
           If (kmaxa /= kmaxa1 .or. kmaxb /= kmaxb1 .or. kmaxc /= kmaxc1) Then
             Write(messages(3),'(a,3i5)') 'DaFT adjusted kmax values (x2) ',kmaxa,kmaxb,kmaxc
             Write(messages(4),'(a,1p,i5)') 'B-spline interpolation order ',mxspl
             Call info(messages,4,.true.)
           Else
             Write(messages(3),'(a,1p,i5)') 'B-spline interpolation order ',mxspl
             Call info(messages,3,.true.)
           End If

! Print infrequent k-space SPME evaluation

           If      (nstfce == 0) Then
              Call warning(370,Real(nstfce,wp),1.0_wp,0.0_wp)
              nstfce=1
           Else If (nstfce > 10) Then
              Call warning(370,Real(nstfce,wp),4.0_wp,0.0_wp)
              nstfce=4
           End If
           If (nstfce >= 1) Then
              Write(message,'(a,1p,i5)') 'k-space evaluation interval (steps)',nstfce
              Call info(message,.true.)
           End If

           If (lforc) Call error(416)
           lforc=.true.

        End If

! read distance dependent dielectric option

     Else If (word(1:6) == 'distan') Then

        keyfce = 4
        Call info('Electrostatics : Distance Dependent Dielectric',.true.)

        If (lforc) Call error(416)
        lforc=.true.

! read coulombic potential option

     Else If (word(1:4) == 'coul') Then

        keyfce = 6
        Call info('Electrostatics : Coulombic Potential',.true.)

        If (lforc) Call error(416)
        lforc=.true.

! read force-shifted coulombic potential option

     Else If (word(1:5) == 'shift') Then

        keyfce = 8
        Call info('Electrostatics : Force-Shifted Coulombic Potential',.true.)

        Call get_word(record,word)

        If      (word(1:4) == 'damp') Then
           Call get_word(record,word)
           alpha = Abs(word_2_real(word))
           Write(message,'(a,1p,e12.4)') 'damping parameter (A^-1) ', alpha
           Call info(message,.true.)
        Else If (word(1:9) == 'precision') Then
           Call get_word(record,word)
           eps0 = Abs(word_2_real(word))
           Write(message,'(a,1p,e12.4)') 'precision parameter ', eps0
           Call info(message,.true.)
           eps0 = Max(Min(eps0,0.5_wp),1.0e-20_wp)
           tol = Sqrt(Abs(Log(eps0*neigh%cutoff)))
           alpha = Sqrt(Abs(Log(eps0*neigh%cutoff*tol)))/neigh%cutoff
           Write(message,'(a,1p,e12.4)') 'damping parameter (A^-1) derived ', alpha
           Call info(message,.true.)
        End If
        If (alpha > zero_plus) Then
           Call info('Fennell damping applied',.true.)
           If (neigh%cutoff < 12.0_wp) Call warning(7,neigh%cutoff,12.0_wp,0.0_wp)
        End If

        If (lforc) Call error(416)
        lforc=.true.

! read reaction field option

     Else If (word(1:8) == 'reaction') Then

        keyfce = 10
        Call info('Electrostatics : Reaction Field',.true.)

        If (word(1:5) == 'field') Call get_word(record,word)
        Call get_word(record,word)

        If      (word(1:4) == 'damp') Then
           Call get_word(record,word)
           alpha = Abs(word_2_real(word))
           Write(message,'(a,1p,e12.4)') 'damping parameter (A^-1) ', alpha
           Call info(message,.true.)
        Else If (word(1:9) == 'precision') Then
           Call get_word(record,word)
           eps0 = Abs(word_2_real(word))
           Write(message,'(a,1p,e12.4)') 'precision parameter ', eps0
           Call info(message,.true.)
           eps0 = Max(Min(eps0,0.5_wp),1.0e-20_wp)
           tol = Sqrt(Abs(Log(eps0*neigh%cutoff)))
           alpha = Sqrt(Abs(Log(eps0*neigh%cutoff*tol)))/neigh%cutoff
           Write(message,'(a,1p,e12.4)') 'damping parameter (A^-1) derived ', alpha
           Call info(message,.true.)
        End If
        If (alpha > zero_plus) Then
           Call info('Fennell damping applied',.true.)
           If (neigh%cutoff < 12.0_wp) Call warning(7,neigh%cutoff,12.0_wp,0.0_wp)
        End If

        If (lforc) Call error(416)
        lforc=.true.

     Else If (word(1:5) == 'poiss' .or. word(1:5) == 'psolv' ) Then

        keyfce = 12
        Call info('Electrostatics : Poisson equation solver',.true.)

        prmps=0.0_wp
        Do i=1,4
           Call get_word(record,word)

           If (word(1:5) == 'delta') Then   ! spacing
              Call get_word(record,word)
              prmps(1)=Abs(word_2_real(word))
           End If

           If (word(1:3) == 'eps') Then     ! tolerance
              Call get_word(record,word)
              prmps(2)=Abs(word_2_real(word))
           End If

           If (word(1:6) == 'maxits') Then  ! max number of iteration
              Call get_word(record,word)
              prmps(3)=Abs(word_2_real(word))
           End If

           If (word(1:7) == 'jmaxits') Then ! max number Jacobian iterations
              Call get_word(record,word)
              prmps(4)=Abs(word_2_real(word))
           End If
        End Do

        Write(messages(1),'(a,1p,e12.4)') 'gridspacing parameter (A) ',prmps(1)
        Write(messages(2),'(a,1p,e12.4)') 'convergance epsilon ',prmps(2)
        Write(messages(3),'(a,1p,i5)') 'max # of Psolver iterations ',Nint(prmps(3))
        Write(messages(4),'(a,1p,i5)') 'max # of Jacobi  iterations ',Nint(prmps(4))
        Call info(messages,4,.true.)

        If ( Abs(prmps(1)-1.0_wp/alpha) > 1.0e-6_wp .or. Abs(prmps(2)-pois%eps) > 1.0e-6_wp .or. &
             Nint(prmps(3)) == 0 .or. Nint(prmps(4)) == 0 ) Then
           Call warning('parameters reset to safe defaults occurred',.true.)
           Write(messages(1),'(a,1p,e12.4)') 'gridspacing parameter (A) ',1.0_wp/alpha
           Write(messages(2),'(a,1p,e12.4)') 'convergance epsilon ',pois%eps
           Write(messages(3),'(a,1p,i5)') 'max # of Psolver iterations ',pois%mxitcg
           Write(messages(4),'(a,1p,i5)') 'max # of Jacobi  iterations ',pois%mxitjb
           Call info(messages,4,.true.)
        End If

        If (lforc) Call error(416)
        lforc=.true.

! read relative dielectric constant

     Else If (word(1:3) == 'eps') Then

        Call get_word(record,word)
        If (word(1:8) == 'constant') Call get_word(record,word)
        epsq = word_2_real(word)
        Write(message,'(a,1p,e12.4)') 'relative dielectric constant ',epsq
        Call info(message,.true.)

! read option for accounting for extended coulombic exclusion

     Else If (word(1:5) == 'exclu') Then

        lecx = .true.
        Call info('Extended Coulombic eXclusion opted for',.true.)

! read force capping option

     Else If (word(1:3) == 'cap') Then

        lfcap = .true.

        Call get_word(record,word)
        If (word(1:5) == 'force') Call get_word(record,word)

        tmp = Abs(word_2_real(word))
        If (tmp > zero_plus) fmax=tmp
        Write(messages(1),'(a)') 'force capping on (during equilibration)'
        Write(messages(2),'(a,1p,e12.4)') 'force capping limit (kT/Angs)',fmax
        Call info(messages,2,.true.)

! read 'no vdw', 'no elec', 'no ind' and 'no vafav' options

     Else If (word(1:2) == 'no') Then

        Call get_word(record,word1)

        If      (word1(1:3) == 'vdw' ) Then

        Else If (word1(1:4) == 'elec') Then

        Else If (word1(1:3) == 'ind' ) Then

        Else If (word1(1:3) == 'str' ) Then

           Call info('no strict option on',.true.)
           Write(banner(1),'(a)') '*** It skips printing inessential information in OUTPUT such as many       ***'
           Write(banner(2),'(a)') '*** warnings, FIELD digested information and full iteration cycles         ***'
           Write(banner(3),'(a)') '*** information from CGM based routines!  However, it also assumes some,   ***'
           Write(banner(4),'(a)') '*** deemed safe, defaults for some specified as well as unspecified by the ***'
           Write(banner(5),'(a)') '*** user options, that may or may not be needed for the simulation to run! ***'
           Write(banner(6),'(a)') '*** The defaults are deemed to deliver safer passage as well as optimal    ***'
           Write(banner(7),'(a)') '*** performance without sacrificing on accuracy!  While it may, by chance, ***'
           Write(banner(8),'(a)') '*** help to pass previously failing runs it may as well lead to a run      ***'
           Write(banner(9),'(a)') '*** failure without warnings!  Beware, avoid usage if uncertain!           ***'
           Call info(banner,9,.true.)

        Else If (word1(1:3) == 'top' ) Then

           Call info('no topology option on (avoids printing extended FIELD topology in OUTPUT)',.true.)

           l_top = .false.

        Else If (word1(1:5) == 'vafav') Then

           green%l_average = .false.

        Else If (word1(1:3) == 'vom' ) Then ! "no vom" should be used with TTM

           If (.not. l_ttm) Then
             Call info('"no vom" option auto-switched on - COM momentum removal will be abandoned',.true.)
             Call warning('this may lead to a build up of the COM momentum ' &
               //'and a manifestation of the "flying ice-cube" effect',.true.)
           End If

           l_vom    = .false.

        Else If (word1(1:4) == 'link') Then ! NON-TRANSFERABLE OPTION FROM DL_POLY_2

           Call warning(38,0.0_wp,0.0_wp,0.0_wp)

        Else

           Call strip_blanks(record)
           Write(message,'(3a)') word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call info(message,.true.)
           Call error(3)

        End If

! read tolerance for relaxed shell model

     Else If (word(1:6) == 'rlxtol') Then

        Call get_word(record,word)
        rlx_tol(1) = Max(1.0_wp,Abs(word_2_real(word)))
        Write(message,'(a,1p,e12.4)') 'relaxed shell model CGM tolerance ',rlx_tol(1)
        Call info(message,.true.)

        Call get_word(record,word1)
        rlx_tol(2) = word_2_real(word1,-1.0_wp)
        If (rlx_tol(2) > zero_plus) Then
           Write(message,'(a,1p,e12.4)') 'relaxed shell model CGM step ',rlx_tol(2)
           Call info(message,.true.)
        End If


! read maximum number of iterations in constraint algorithms

     Else If (word(1:6) == 'mxshak') Then

        Call get_word(record,word)
        cons%max_iter_shake = Abs(Nint(word_2_real(word)))

! read tolerance for constraint algorithms

     Else If (word(1:5) == 'shake') Then

        Call get_word(record,word)
        If (word(1:9) == 'tolerance') Call get_word(record,word)
        cons%tolerance = Abs(word_2_real(word))

! read maximum number of iterations in LFV quaternion integration algorithms

     Else If (word(1:6) == 'mxquat') Then

        Call get_word(record,word)
        mxquat = Abs(Nint(word_2_real(word)))

! read tolerance in LFV quaternion integration algorithms

     Else If (word(1:6) == 'quater') Then

        Call get_word(record,word)
        If (word(1:9) == 'tolerance') Call get_word(record,word)
        quattol = Abs(word_2_real(word))

! read two-temperature model (ttm) specific flags

     Else If (word(1:3) == 'ttm') Then

        ! detecting l_ttm purely to print message on first occasion

        If (.not. l_ttm) Then
          l_ttm = .true.
          Call info('Two Temperature Model (TTM) opted for',.true.)
        End If

        Call get_word(record,word1)

        If (word1(1:4) == 'ncit') Then

        ! number of coarse-grained ion temperature cells (CIT):
        ! already determined in scan_control

           Write(messages(1),'(a,3(1x,i8))') 'ionic temperature grid size (x,y,z):',ntsys(1:3)
           Write(messages(2),'(a,3(1x,f8.4))') 'temperature grid size (x,y,z):',delx,dely,delz
           Write(messages(3),'(a,f10.4)') 'average number of atoms/cell: ',sysrho*volume
           Call info(messages,3,.true.)

        Else If (word1(1:4) == 'ncet') Then

        ! number of coarse-grained electronic temperature cells (CET):
        ! already determined in scan_control

           Write(message,'(a,3(1x,i8))') 'electronic temperature grid size (x,y,z):', &
             eltsys(1:3)
           Call info(message,.true.)

        Else If (word1(1:5) == 'metal') Then

        ! sets properties of electronic subsystem as a metal:
        ! already determined in scan_control

          Call info('electronic subsystem represents metal: thermal conductivity required',.true.)

        Else If (word1(1:8) == 'nonmetal') Then

        ! sets properties of electronic subsystem as a non-metal

          Call info('electronic subsystem represents non-metal: thermal diffusivity required',.true.)

        Else If (word1(1:7) == 'ceconst') Then

        ! electronic specific heat capacity given as constant value

          Call get_word(record,word)
          Ce0 = word_2_real(word)
          Call info('electronic specific heat capacity set to constant value',.true.)
          Write(message,'(a,1p,e12.4)') 'electronic s.h.c. (kB/atom) ', Ce0
          Call info(message,.true.)

        Else If (word1(1:6) == 'cetanh') Then

        ! electronic specific heat capacity given as tanh function

          Call get_word(record,word)
          sh_A = word_2_real(word)
          Call get_word(record,word)
          sh_B = word_2_real(word)
          Write(messages(1),'(a)') 'electronic specific heat capacity set to hyperbolic tangent function'
          Write(messages(2),'(a,1p,e12.4)') 'constant term A (kB/atom) ',sh_A
          Write(messages(3),'(a,1p,e12.4)') 'emperature term B (K^-1) ',sh_B
          Call info(messages,3,.true.)

        Else If (word1(1:5) == 'celin') Then

        ! electronic specific heat capacity given as linear function
        ! up to Fermi temperature, constant afterwards

          Call get_word(record,word)
          Cemax = word_2_real(word)
          Call get_word(record,word)
          Tfermi = word_2_real(word)
          Write(messages(1),'(a)') 'electronic specific heat capacity set to linear function up to Fermi temperature'
          Write(messages(2),'(a,1p,e12.4)') 'max. electronic s.h.c. (kB/atom) ',Cemax
          Write(messages(3),'(a,1p,e12.4)') 'Fermi temperature (K)',Tfermi
          Call info(messages,3,.true.)

        Else If (word1(1:5) == 'cetab') Then

        ! electronic volumetric heat capacity given in tabulated form

          Call info('electronic volumetric heat capacity given as tabulated function of temperature',.true.)

        Else If (word1(1:5) == 'keinf') Then

        ! infinite electronic thermal conductivity

          Call info('electronic thermal conductivity set to infinity',.true.)

        Else If (word1(1:7) == "keconst") Then

        ! electronic thermal conductivity given as constant value

          Call get_word(record,word)
          Ka0 = word_2_real(word)
          Write(messages(1),'(a)') 'electronic thermal conductivity set to constant value'
          Write(messages(2),'(a,1p,e12.4)') 'electronic t.c. (W m^-1 K^-1) ',Ka0
          Call info(messages,2,.true.)

        Else If (word1(1:7) == 'kedrude') Then

        ! electronic thermal conductivity given as drude model (propertional to
        ! electronic temperature, giving t.c. at system temperature)

          Call get_word(record,word)
          Ka0 = word_2_real(word)
          Write(messages(1),'(a)') 'electronic thermal conductivity set to drude model'
          Write(messages(2),'(a,1p,e12.4)') 't.c. at system thermo%temp. (W m^-1 K^-1) ',Ka0
          Call info(messages,2,.true.)

        Else If (word1(1:5) == 'ketab') Then

        ! electronic thermal conductivity given in tabulated form

          Write(messages(1),'(a)') 'electronic thermal conductivity given as tabulated function of temperature:'
          Write(messages(2),'(a)') 'uses ionic or system temperature to calculate cell conductivity value'
          Write(messages(3),'(a)') 'for thermal diffusion equation'
          Call info(messages,3,.true.)

        Else If (word1(1:4) == 'diff' .or. word1(1:7)=='deconst') Then

        ! electronic thermal diffusivity given as constant value
        ! (for non-metal systems)

          Call get_word(record,word)
          Diff0 = word_2_real(word)
          Write(messages(1),'(a)') 'electronic thermal diffusivity set to constant value'
          Write(messages(2),'(a,1p,e12.4)') 'electronic t.d. (m^2 s^-1) ',Diff0
          Call info(messages,2,.true.)

        Else If (word1(1:7) == 'derecip') Then

        ! electronic thermal diffusivity given as reciprocal function
        ! of temperature (up to Fermi temperature), constant afterwards

          Call get_word(record,word)
          Diff0 = word_2_real(word)
          Call get_word(record,word)
          Tfermi = word_2_real(word)
          Write(messages(1),'(a)') 'electronic thermal diffusivity set to reciprocal function up to Fermi temperature'
          Write(messages(2),'(a,1p,e12.4)') 'datum electronic t.d. (m^2 s^-1) ',Diff0
          Write(messages(3),'(a,1p,e12.4)') 'Fermi temperature (K) ',Tfermi
          Call info(messages,3,.true.)

        Else If (word1(1:4) == 'detab') Then

        ! electronic thermal diffusivity given in tabulated form

          Call info('electronic thermal diffusivity given as tabulated function of temperature',.true.)

        Else If (word1(1:8) == 'atomdens') Then

        ! user-specified atomic density, used to convert specific
        ! heat capacities to volumetric values

          Call get_word(record,word)
          cellrho = word_2_real(word)
          Write(message,'(a,f10.4)') 'user-specified atomic density (A^-3) ',cellrho
          Call info(message,.true.)

        Else If (word1(1:7) == 'dyndens') Then

        ! dynamic calculation of atom density in active cells during
        ! TTM calculations, used to convert specific heat capacities
        ! to volumetric values

          ttmdyndens = .true.
          Call info('dynamic calculations of average atomic density in active ionic cells',.true.)

        Else If (word1(1:4) == 'amin') Then

        ! minimum number of atoms needed per ionic temperature cell
        ! to give definable ionic temperature (default = 1): smaller
        ! number deactivates ionic and electronic temperature cells
        ! (by default, electronic energies are not redistributed)

          Call get_word(record,word)
          amin = Abs(Nint(word_2_real(word)))
          Write(message,'(a,1p,i8)') 'min. atom no. for ionic cells ',amin
          Call info(message,.true.)

        Else If (word1(1:6) == 'redist') Then

        ! redistribution of electronic energy from deactivated cells 
        ! to active neighbours

          If (redistribute) Then
            Write(messages(1),'(a)') 'redistributing energy from deactivated electronic cells into active neighbours'
            Write(messages(2),'(a)') '(requires at least one electronic temperature cell beyond ionic cells)'
            Call info(messages,2,.true.)
          End If

        Else If (word1(1:4) == 'dedx') Then

        ! electronic stopping power of projectile entering electronic system

          Call get_word(record,word)
          dEdX = word_2_real(word)
          Write(message,'(a,1p,e12.4)') 'elec. stopping power (eV/nm) ',dEdX
          Call info(message,.true.)

        Else If (word1(1:6) == 'sgauss' .or. word1(1:5) == 'thermo%sigma') Then

        ! gaussian spatial distribution for initial energy deposition into
        ! electronic system

          sdepoType = 1
          Call get_word(record,word)
          sig = word_2_real(word)
          Call get_word(record,word)
          sigmax = word_2_real(word)
          Write(messages(1),'(a)') 'initial gaussian spatial energy deposition in electronic system'
          Write(messages(2),'(a,1p,e12.4)') 'thermo%sigma of distribution (nm) ',sig
          Write(messages(3),'(a,1p,e12.4)') 'distribution cutoff (nm) ',sigmax*sig
          Call info(messages,3,.true.)

        Else If (word1(1:5) == 'sflat') Then

        ! homogeneous spatial distribution for initial energy deposition into
        ! electronic system

          sdepoType = 2
          Call info('initial homogeneous (flat) spatial energy deposition in electronic system',.true.)

        Else If (word1(1:5) == 'laser') Then

        ! homogeneous spatial distribution for initial energy deposition into
        ! electronic system due to laser: setting absorbed fluence and
        ! penetration depth

          sdepoType = 2
          Call get_word(record,word)
          fluence = word_2_real(word)
          Call get_word(record,word)
          pdepth = word_2_real(word)
          Call get_word(record,word)
          If (word(1:4) == 'zdep') sdepoType = 3
          Select Case (sdepoType)
          Case (2)
            Write(messages(1),'(a)') 'initial homogeneous (flat) spatial ' &
              //'energy deposition in electronic system due to laser'
            Write(messages(2),'(a,1p,e12.4)') &
              'absorbed fluence (mJ cm^-2) ',fluence
            Write(messages(3),'(a,1p,e12.4)') 'penetration depth (nm) ',pdepth
          Case (3)
            Write(messages(1),'(a)') 'initial xy-homogeneous, z-exponential ' &
              //'decaying spatial energy deposition in electronic system due to laser'
            Write(messages(2),'(a,1p,e12.4)') &
              'absorbed fluence at surface (mJ cm^-2) ',fluence
            Write(messages(3),'(a,1p,e12.4)') 'penetration depth (nm) ',pdepth
          End Select
          Call info(messages,3,.true.)

        Else If (word1(1:5) == 'gauss') Then

        ! gaussian temporal distribution for energy deposition into
        ! electronic system

          tdepoType = 1
          Call get_word(record,word)
          tdepo = word_2_real(word)
          Call get_word(record,word)
          tcdepo = word_2_real(word)
          Write(messages(1),'(a)') 'gaussian temporal energy deposition in electronic system'
          Write(messages(2),'(a,1p,e12.4)') 'thermo%sigma of distribution (ps) ',tdepo
          Write(messages(3),'(a,1p,e12.4)') 'distribution cutoff (ps) ',2.0_wp*tcdepo*tdepo
          Call info(messages,3,.true.)

        Else If (word1(1:5) == 'nexp') Then

        ! decaying exponential temporal distribution for energy deposition into
        ! electronic system

          tdepoType = 2
          Call get_word(record,word)
          tdepo = word_2_real(word)
          Call get_word(record,word)
          tcdepo = word_2_real(word)
          Write(messages(1),'(a)') 'decaying exponential temporal energy deposition in electronic system'
          Write(messages(2),'(a,1p,e12.4)') 'tau of distribution (ps) ',tdepo
          Write(messages(3),'(a,1p,e12.4)') 'distribution cutoff (ps) ',tcdepo*tdepo
          Call info(messages,3,.true.)

        Else If (word1(1:5) == 'delta') Then

        ! dirac delta temporal distribution for energy deposition into
        ! electronic system

          tdepoType = 3
          Call info('dirac delta temporal energy deposition in electronic system',.true.)

        Else If (word1(1:5) == 'pulse') Then

        ! square pulse temporal distribution for energy deposition into
        ! electronic system (defaults to dirac delta if pulse duration
        ! set to zero)

          tdepoType = 4
          Call get_word(record,word)
          tdepo = word_2_real(word)
          If (tdepo<=zero_plus) Then
            tdepoType = 3
            Call info('square pulse temporal energy deposition in electronic' &
              //'system of zero duration: being treated as dirac delta')
          Else
            Write(messages(1),'(a)') 'square pulse temporal energy deposition in electronic system'
            Write(messages(2),'(a,1p,e12.4)') 'pulse duration (ps) ',tdepo
            Call info(messages,2,.true.)
          End If

        Else If (word1(1:4) == 'varg') Then

        ! variable electron-phonon coupling constant (thermo%chi_ep) based on
        ! tabular electronic stopping terms (in g.dat file): option to
        ! apply value homogeneously across system (based on average 
        ! electronic temperature) or heterogeneously (using local 
        ! electronic temperature for each voxel)

          Select Case (gvar)
          Case (1)
            Write(messages(1),'(a)') 'variable electron-phonon coupling values to be applied homogeneously'
          Case (2)
            Write(messages(1),'(a)') 'variable electron-phonon coupling values to be applied heterogeneously'
          End Select
          Write(messages(2),'(a)') '(overrides value given for ensemble, required tabulated stopping'
          Write(messages(3),'(a)') 'terms in g.dat file)'
          Call info(messages,3,.true.)

        Else If (word1(1:3) == 'bcs') Then

        ! electronic temperature boundary conditions

          Call get_word(record,word)

          If (word(1:8) == 'periodic') Then
            bcTypeE = 1
            Call info('electronic temperature boundary conditions set as periodic',.true.)
          Else If (word(1:6) == 'dirich') Then
            bcTypeE = 2
            Write(messages(1),'(a)') 'electronic temperature boundary conditions set as dirichlet:'
            Write(messages(2),'(a)') 'setting boundaries to system temperature'
            Call info(messages,2,.true.)
          Else If (word(1:7) == 'neumann') Then
            bcTypeE = 3
            Write(messages(1),'(a)') 'electronic temperature boundary conditions set as neumann:'
            Write(messages(2),'(a)') 'zero energy flux at boundaries'
            Call info(messages,2,.true.)
          Else If (word(1:8) == 'xydirich') Then
            bcTypeE = 4
            Write(messages(1),'(a)') 'electronic temperature boundary conditions set as dirichlet (xy), neumann (z):'
            Write(messages(2),'(a)') 'system temperature at x and y boundaries'
            Write(messages(3),'(a)') 'zero energy flux at z boundaries'
            Call info(messages,3,.true.)
          Else If (word(1:5) == 'robin') Then
            bcTypeE = 5
            Call get_word(record,word)
            fluxout = word_2_real(word)
            Write(messages(1),'(a)') 'electronic temperature boundary conditions set as robin:'
            Write(messages(2),'(a,1p,e11.4)') 'temperature leakage at boundaries of ',fluxout
            Call info(messages,2,.true.)
          Else If (word(1:7) == 'xyrobin') Then
            bcTypeE = 6
            Call get_word(record,word)
            fluxout = word_2_real(word)
            Write(messages(1),'(a)') 'electronic temperature boundary conditions set as robin (xy), neumann (z):'
            Write(messages(2),'(a,1p,e11.4)') 'temperature leakage at x and y boundaries of ',fluxout
            Write(messages(3),'(a)') 'zero energy flux at z boundaries'
            Call info(messages,3,.true.)
          End If

        Else If (word1(1:6) == 'offset') Then

        ! time offset in coupling electronic and ionic systems

          Call get_word(record,word)
          ttmoffset = word_2_real(word)
          Write(message,'(a,1p,e12.4)') 'electron-ion coupling offset (ps) ',ttmoffset
          Call info(message,.true.)

        Else If (word1(1:6) == 'oneway') Then

        ! one-way electron-phonon coupling in thermostat and thermal
        ! diffusion: only apply when electronic temperature exceeds
        ! ionic temperature

          oneway = .true.
          Call info('one-way electron-phonon coupling option switched on',.true.)

        Else If (word1(1:5) == 'stats') Then

        ! ttm statistics (minimum and maximum ionic/electronic temperatures,
        ! electronic energy) file option and output frequency

          Call get_word(record,word)
          ttmstats = Abs(Nint(word_2_real(word)))
          Write(messages(1),'(a)') 'ttm statistics file option on'
          Write(messages(2),'(a,i10)') 'ttm statistics file interval ',ttmstats
          Call info(messages,2,.true.)

        Else If (word1(1:4) == 'traj') Then

        ! ttm trajectory (one-dimensional ionic and electronic 
        ! temperature profile) file option and output frequency

          Call get_word(record,word)
          ttmtraj = Abs(Nint(word_2_real(word)))
          Write(messages(1),'(a)') 'ttm trajectory (temperature profile) file option on'
          Write(messages(2),'(a,i10)') 'ttm trajectory file interval',ttmtraj
          Call info(messages,2,.true.)

        End If

! read replay history option

     Else If (word(1:6) == 'replay') Then

!        lsim = .false. ! done in scan_control
        Call get_word(record,word)
        If (word(1:4) == 'hist') Call get_word(record,word)
        If (word(1:5) == 'force') lfce=.true.

! read binsize option

     Else If (word(1:7) == 'binsize') Then

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (Abs(rbin-tmp) > 1.0e-6_wp) Call warning(340,tmp,neigh%cutoff/4.0_wp,rbin)

! read analysis (intramolecular distributions calculation) option

     Else If (word(1:3) == 'ana') Then

        Call get_word(record,word1)
        akey = word1(1:3)

        If (akey /= 'all' .and. akey /= 'bon' .and. akey /= 'ang' .and. &
            akey /= 'dih' .and. akey /= 'inv') Then
           Call strip_blanks(record)
           Write(message,'(3a)') word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call info(message,.true.)
           Call error(3)
        End If

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        i=Abs(Nint(word_2_real(word,0.0_wp))) ! frequency

        Call get_word(record,word)
        If (word(1:5) == 'nbins' .or. word(1:5) == 'ngrid' .or. word(1:4) == 'grid') Then
           Call get_word(record,word)
           j=Abs(Nint(word_2_real(word))) ! grid size
        Else
           j=0
        End If

        If      (akey == 'all') Then
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           tmp=Abs(word_2_real(word)) ! bond length

           nstana=Max(nstana,i)
           grdana=Max(grdana,j)
           rcb_d =Max(rcb_d,tmp)
        Else If (akey == 'bon') Then
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           tmp=Abs(word_2_real(word)) ! bond length

           nstbnd=Max(nstbnd,i)
           grdbnd=j
           rcb_d =Max(rcb_d,tmp)
        Else If (akey == 'ang') Then
           nstang=Max(nstang,i)
           grdang=j
        Else If (akey == 'dih') Then
           nstdih=Max(nstdih,i)
           grddih=j
        Else If (akey == 'inv') Then
           nstinv=Max(nstinv,i)
           grdinv=j
        End If
        grdana=Max(grdana,grdbnd,grdang,grddih,grdinv)

! read rdf%rdf calculation option

     Else If (word(1:3) == 'rdf%rdf') Then

        rdf%l_collect = .true.

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        rdf%freq = Abs(Nint(word_2_real(word,1.0_wp)))

! read z-density profile option

     Else If (word(1:4) == 'zden') Then

        zdensity%l_collect = .true.

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        zdensity%frequency = Abs(Nint(word_2_real(word,1.0_wp)))

! read vaf calculation option dealt with in scan_control<-set_bounds

     Else If (word(1:3) == 'vaf') Then

! read thermal conductivity calculation option
!
!     Else If (word(1:5) == 'therm') Then
!
!        ltcond = .true.
!
!        Call get_word(record,word)
!        If (word(1:4) == 'cond' .or. word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:4) == 'over') &
!           Call get_word(record,word)
!        If (word(1:4) == 'cond' .or. word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:4) == 'over') &
!           Call get_word(record,word)
!        If (word(1:4) == 'cond' .or. word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:4) == 'over') &
!           Call get_word(record,word)
!        nsttcond = Max(Abs(Nint(word_2_real(word))),1)

! read print options

     Else If (word(1:5) == 'print') Then

        Call get_word(record,word)

        If      (word(1:3) == 'ana' ) Then
           lpana = .true.
        Else If (word(1:3) == 'rdf%rdf' ) Then
           rdf%l_print = .true.
        Else If (word(1:4) == 'zden') Then
           zdensity%l_print = .true.
        Else If (word(1:3) == 'vaf') Then
           green%l_print = .true.
        Else
           If (word(1:5) == 'every') Call get_word(record,word)
           nstbpo = Abs(Nint(word_2_real(word,1.0_wp)))
           Write(message,'(a,i10)') 'data printing interval (steps) ',nstbpo
           Call info(message,.true.)
        End If

! read stack option (reading done in set_bounds -> scan_control)

     Else If (word(1:5) == 'stack') Then

        Write(message,'(a,i10)') 'data stacking interval (steps) ',stats%mxstak
        Call info(message,.true.)

! read statistics printing option

     Else If (word(1:4) == 'stat') Then

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        stats%intsta = Nint(word_2_real(word))
        Write(message,'(a,i10)') 'statistics file interval ',stats%intsta
        Call info(message,.true.)

! read MSDTMP printing option

     Else If (word(1:6) == 'msdtmp') Then

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        msd_data%start = Max(msd_data%start,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        msd_data%freq = Max(msd_data%freq,itmp)

        Write(messages(1),'(a)') 'MSDTMP file option on'
        Write(messages(2),'(2x,a,i10)') 'MSDTMP file start ',msd_data%start
        Write(messages(3),'(2x,a,i10)') 'MSDTMP file interval ',msd_data%freq
        Call info(messages,3,.true.)

! read trajectory printing option

     Else If (word(1:4) == 'traj') Then

        ltraj = .true.

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        nstraj = Max(nstraj,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        istraj = Max(istraj,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        keytrj = Max(keytrj,itmp)

        Write(messages(1),'(a)') 'trajectory file option on'
        Write(messages(2),'(2x,a,i10)') 'trajectory file start ',nstraj
        Write(messages(3),'(2x,a,i10)') 'trajectory file interval ',istraj
        Write(messages(4),'(2x,a,i10)') 'trajectory file info key ',keytrj
        Call info(messages,4,.true.)

        If (keytrj > 3) Call error(517)
        If (keytrj == 3) Then
          Call warning('trajectory file info key == 3 generates HISTORY in an unindexed and consize manner',.true.)
        End If

! read defects trajectory printing option

     Else If (word(1:4) == 'defe') Then

        dfcts(1)%ldef = .true.

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        dfcts(1)%nsdef = Max(dfcts(1)%nsdef,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        dfcts(1)%isdef = Max(dfcts(1)%isdef,itmp)

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (tmp >= Min(0.3_wp,neigh%cutoff/3.0_wp) .and. tmp <= Min(3.5_wp,neigh%cutoff/2.0_wp)) Then
           dfcts(1)%rdef = tmp ! 3.43 Angs is the Cs VDW radius - largest possible
        Else
           Call warning(310,tmp,dfcts(1)%rdef,0.0_wp)
        End If


        Write(messages(1),'(a)') 'defects file option on'
        Write(messages(2),'(2x,a,i10)') 'defects file start ',dfcts(1)%nsdef
        Write(messages(3),'(2x,a,i10)') 'defects file interval ',dfcts(1)%isdef
        Write(messages(4),'(2x,a,1p,e12.4)') 'defects distance condition (Angs) ',dfcts(1)%rdef
        Call info(messages,4,.true.)
        dfcts(1)%reffile = 'REFERENCE'
        dfcts(1)%deffile = 'DEFECTS'

! REFERENCE1 forcing
       Call get_word(record,word)
       If (word(1:5) == 'extra') Then
         dfcts(2)%ldef=.true.
         Call info('defects1 file option on',.true.)
         dfcts(2)%nsdef =dfcts(1)%nsdef
         dfcts(2)%isdef =dfcts(1)%isdef
         dfcts(2)%rdef  =dfcts(1)%rdef
         dfcts(2)%newjob = .True.
         ! Name REFERENCE and DEFECTS files
         dfcts(2)%reffile = 'REFERENCE1'
         dfcts(2)%deffile = 'DEFECTS1'
       End If


! read displacements trajectory printing option

     Else If (word(1:4) == 'disp') Then

        lrsd = .true.

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        nsrsd = Max(nsrsd,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        isrsd = Max(isrsd,itmp)

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (tmp > rrsd) Then
           rrsd = tmp
        Else
           Call warning(470,tmp,rrsd,0.0_wp)
        End If

        Write(messages(1),'(a)') 'displacements file option on'
        Write(messages(2),'(2x,a,i10)') 'DISPDAT file start ',nsrsd
        Write(messages(3),'(2x,a,i10)') 'DISPDAT file interval ',isrsd
        Write(messages(4),'(2x,a,1p,e12.4)') 'DISPDAT distance condition (Angs) ',rrsd
        Call info(messages,4,.true.)

! read DL_POLY_2/Classic delr Verlet shell strip cutoff option (compatibility)
! as DL_POLY_4 real space cutoff padding option

     Else If (word(1:4) == 'delr') Then

        Call warning(35,0.0_wp,0.0_wp,0.0_wp)
        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        rpad1 = 0.25_wp * Abs(word_2_real(word))
        Write(message,'(a,1p,e12.4)') 'cutoff padding (Angs) ',rpad1
        Call info(message,.true.)

! read DL_POLY_2/Classic multiple timestep option (compatibility)
! as DL_POLY_4 infrequent k-space SPME evaluation option

     Else If (word(1:4) == 'mult') Then

        Call warning(36,0.0_wp,0.0_wp,0.0_wp)

!!!! OTHER NON-TRANSFERABLE OPTIONS FROM DL_POLY_2/Classic !!!!
! read primary cutoff option for multiple timestepping

     Else If (word(1:4) == 'prim') Then

        Call warning(34,0.0_wp,0.0_wp,0.0_wp)

! read all pairs option

     Else If (word(1:3) == 'all') Then

        Call warning(37,0.0_wp,0.0_wp,0.0_wp)

!!!! OTHER NON-TRANSFERABLE OPTIONS FROM DL_POLY_2/Classic !!!!

! read data dumping interval

     Else If (word(1:4) == 'dump') Then

        Call get_word(record,word)
        If (word(1:4) == 'data' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:4) == 'data' .or. word(1:5) == 'every') Call get_word(record,word)
        ndump = Max(Abs(Nint(word_2_real(word))),1)

! default for particle density per link cell below
! which decreasing link-cell size (subcelling) stops

     Else If (word(1:7) == 'subcell') Then

        Call get_word(record,word)
        If (word(1:4) == 'dens' .or. word(1:6) == 'thresh') Call get_word(record,word)
        If (word(1:4) == 'dens' .or. word(1:6) == 'thresh') Call get_word(record,word)
        pdplnc = Max(Abs(word_2_real(word)),1.0_wp) ! disallow any less than 1

! read machine time for simulation run (in seconds)

     Else If (word(1:3) == 'job') Then

        Call get_word(record,word1)
        If (word1(1:4) == 'time') Then

           l_timjob=.true.

           Call get_word(record,word)
           tmr%job = word_2_real(word)

        Else

           Call strip_blanks(record)
           Write(message,'(3a)') word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call info(message,.true.)
           Call error(3)

        End If

! read close-down time allowance

     Else If (word(1:5) == 'close') Then

        Call get_word(record,word1)
        If (word1(1:4) == 'time') Then

           l_timcls=.true.

           Call get_word(record,word)
           tmr%clear_screen = word_2_real(word)

        Else

           Call strip_blanks(record)
           Write(message,'(3a)') word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call info(message,.true.)
           Call error(3)

        End If

! close control file

     Else If (word(1:6) == 'finish') Then

        Go To 2000

     Else If (word(1:6) == 'plumed') Then

        If (lplumed) Then
           plume%l_plumed=.true.

           Call get_word(record,word)
           If (word(1:3) == 'off') Then
              plume%l_plumed=.false.
              lplumed=.false.
           End If

           If (word(1:5) == 'input') Then
              Call get_word(record,word)
              plume%input=Trim(word)
           End If

           If (word(1:3) == 'log') Then
              Call get_word(record,word)
              plume%logfile=Trim(word)
           End If

           If (word(1:9) == 'precision') Then
              Call get_word(record,word)
              plume%prec=Abs(Nint(word_2_real(word,1.0_wp)))
           End If

           If (word(1:7) == 'restart') Then
              plume%restart=1

              Call get_word(record,word)

              If ((word(1:3) == 'yes') .or. (word(1:1) == 'y')) Then
                 plume%restart=1
              End If

              If ((word(1:2) == 'no') .or. (word(1:1) == 'n')) Then
                 plume%restart=0
              End If
           End If
        End If

     Else

        Call strip_blanks(record)
        Write(message,'(2a)') word(1:Len_Trim(word)+1),record
        Call info(message,.true.)
        Call error(3)

     End If

  End Do

! no finish record in CONTROL file

  If (comm%idnode == 0) Close(Unit=nread)
  Call error(17)

! unexpected end of file

1000 Continue

  If (comm%idnode == 0) Close(Unit=nread)
  Call error(53)

! safe termination of reading CONTROL

2000 Continue
  If (comm%idnode == 0) Close(Unit=nread)

!!! FIXES !!!
! fix on step-dependent options

  If (minimise%freq  == 0) minimise%freq  = nsteql+1
  If (thermo%freq_zero == 0) thermo%freq_zero = nsteql+1
  If (thermo%freq_tgaus == 0) thermo%freq_tgaus = nsteql+1
  If (thermo%freq_tscale == 0) thermo%freq_tscale = nsteql+1

!!! REPORTS !!!
! report restart

  If (keyres == 0) Then
     Call info('clean start requested',.true.)
  Else If (levcfg == 0) Then
     Call warning(200,0.0_wp,0.0_wp,0.0_wp)
     keyres=0
  End If

! report default ensemble if none is specified:
! inhomogeneous Langevin if two-temperature model
! is in use, NVE if not

  If (.not.lens) Then
    Call warning(130,0.0_wp,0.0_wp,0.0_wp)
    If (l_vv) Then
      Call info('Integration : Velocity Verlet',.true.)
    Else
      Call info('Integration : Leapfrog Verlet',.true.)
    End If
    If (l_ttm) Then
      Write(messages(1),'(a)') 'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)'
      Write(messages(2),'(a,1p,e12.4)') 'e-phonon friction (ps^-1) ',thermo%chi_ep
      Write(messages(3),'(a,1p,e12.4)') 'e-stopping friction (ps^-1) ',thermo%chi_es
      Write(messages(4),'(a,1p,e12.4)') 'e-stopping velocity (A ps^-1)',thermo%vel_es2
      Call info(messages,4,.true.)
    Else
      Call info('Ensemble : NVE (Microcanonical)',.true.)
    End If
    If (l_ttm) thermo%ensemble = ENS_NVT_LANGEVIN_INHOMO
    lens=.true.
  End If

! report replacement of specified ensemble with inhomogeneous
! Langevin if two-temperature model is in use, replacing
! default electron-phonon friction value with thermo%chi from
! standard Langevin thermostat (if supplied), and use of
! thermal velocities only for thermostat

  If (l_ttm .and. thermo%ensemble/=ENS_NVT_LANGEVIN_INHOMO) Then
    Call warning(130,0.0_wp,0.0_wp,0.0_wp)
    If (thermo%ensemble==ENS_NVT_LANGEVIN .or. &
        thermo%ensemble==ENS_NPT_LANGEVIN .or. &
        thermo%ensemble==ENS_NPT_LANGEVIN_ANISO .and. thermo%chi>zero_plus) Then
      thermo%chi_ep = thermo%chi
    End If

    Write(messages(1),'(a)') 'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)'
    Write(messages(2),'(a,1p,e12.4)') 'e-phonon friction (ps^-1) ',thermo%chi_ep
    Write(messages(3),'(a,1p,e12.4)') 'e-stopping friction (ps^-1) ',thermo%chi_es
    Write(messages(4),'(a,1p,e12.4)') 'e-stopping velocity (A ps^-1)',thermo%vel_es2
    Call info(messages,4,.true.)

    If (ttmthvel) Then
      Call info('applying to thermal velocities in all directions',.true.)
    Else If (ttmthvelz) Then
      Call info('applying to total velocities in x and y directions, thermal velocities in z direction',.true.)
    Else
      Call info('applying to total velocities in all directions',.true.)
    End If
    thermo%ensemble = ENS_NVT_LANGEVIN_INHOMO

  End If

! report iteration length and tolerance condition for constraints and PMF algorithms

  If ((cons%mxcons > 0 .or. pmf%mxpmf > 0) .and. comm%idnode == 0) Then
     Write(messages(1),'(a,i10)') 'iterations for shake/rattle ',cons%max_iter_shake
     Write(messages(2),'(a,1p,e12.4)') 'tolerance for shake/rattle (Angs) ',cons%tolerance
     Call info(messages,2,.true.)
  End If

! report electrostatics

  If (l_n_e) Then
     keyfce=0
     Call info('Electrostatics switched off!!!',.true.)
  Else If (keyfce == 0) Then
     Call info('Electrostatics : None Assumed',.true.)
  End If

! report for extended coulombic exclusion if needed

  If (keyfce /= 0) Then
     If (lecx) Then
        Call info('Extended Coulombic eXclusion : YES',.true.)
     Else
        Call info('Extended Coulombic eXclusion : NO',.true.)
     End If
  End If

! report if neigh%cutoff is reset (measures taken in scan_config -
! neigh%cutoff is the maximum cutoff needed in the system)

  If (Abs(neigh%cutoff-rcut1) > 1.0e-6_wp) Then
    Write(message,'(a,1p,e12.4)') 'real space cutoff reset to (Angs) ',neigh%cutoff
    Call info(message,.true.)
  End If

! report if neigh%padding is reset (measures taken in scan_config & set_bounds -
! neigh%padding is the cutoff padding needed the conditional VNL update)

  If (Abs(neigh%padding-rpad1) > 1.0e-6_wp) Then
    Write(message,'(a,1p,e12.4)') 'cutoff padding reset to (Angs) ',neigh%padding
    Call info(message,.true.)
  End If

! report vdw

  If (l_n_v) Then
    Call info('vdw potential terms switched off',.true.)
  End If

! report if vdw%cutoff is reset (measures taken in scan_config)

  If ((.not.l_n_v) .and. Abs(vdw%cutoff-rvdw1) > 1.0e-6_wp) Then
    Write(message,'(a,1p,e12.4)') 'vdw cutoff reset to (Angs) ',vdw%cutoff
    Call info(message,.true.)
  End If

! report timestep

  If (lvar) Then

     If (thermo%key_dpd > 0) Then
        lvar=.false.
        Call warning('variable timestep unavalable in DPD themostats',.true.)
        Write(message,'(a,1p,e12.4)') 'fixed simulation timestep (ps) ',tstep
        Call info(message,.true.)
     Else
        If (mxdis >= 2.5_wp*mndis .and. mndis > 0.0_wp) Then
          Write(messages(1),'(a,1p,e12.4)') 'variable simulation timestep (ps) ',tstep
          Write(messages(2),'(a)') 'controls for variable timestep:'
          Write(messages(3),'(2x,a,1p,e12.4)') 'minimum distance Dmin (Angs) ',mndis
          Write(messages(4),'(2x,a,1p,e12.4)') 'maximum distance Dmax (Angs) ',mxdis
          Call info(messages,4,.true.)

           If (mxstp > zero_plus) Then
              Write(message,'(a,1p,e12.4)') 'timestep ceiling mxstp (ps) ',mxstp
              Call info(message,.true.)
              tstep=Min(tstep,mxstp)
           Else
              mxstp=Huge(1.0_wp)
           End If
        Else
           Call warning(140,mndis,mxdis,0.0_wp)
           Call error(518)
        End If
     End If

   Else If (lstep) Then
     Write(message,'(a,1p,e12.4)') 'fixed simulation timestep (ps) ',tstep
     Call info(message,.true.)
   End If

! report no vom option: its use recommended with ttm

  If (.not.l_vom .and. .not.l_ttm) Then
    Call info('no vom option on - COM momentum removal will be abandoned',.true.)
    Call warning('this may lead to a build up of the COM momentum and a' &
      //'manifestation of the "flying ice-cube" effect',.true.)
  Else If (l_vom .and. l_ttm) Then
    Call info('no vom option off - COM momentum removal will be used',.true.)
    Call warning('this may lead to incorrect dynamic behaviour for' &
      //'two-temperature model: COM momentum removal recommended',.true.)
  End If

! report intramolecular analysis options

  If (lpana .or. mxgana > 0) Then
     If (mxgana == 0) Then
        Call info('no intramolecular distribution collection requested',.true.)
     Else
        If (bond%bin_pdf > 0 .and. angle%bin_adf > 0 .and. &
            dihedral%bin_adf > 0 .and. inversion%bin_adf > 0) Then
           Call info('full intramolecular distribution collection requested (all=bnd/ang/dih/inv):',.true.)
        Else
           Call info('intramolecular distribution collection requested for:',.true.)
        End If

        i=Max(1,nstana,nstbnd,nstang,nstdih,nstinv)
        nstall=Min(i , Merge(nstana , i , nstana > 0), &
                       Merge(nstbnd , i , nstbnd > 0), &
                       Merge(nstang , i , nstang > 0), &
                       Merge(nstdih , i , nstdih > 0), &
                       Merge(nstinv , i , nstinv > 0))

        If (bond%bin_pdf > 0) Then
           If (nstbnd == 0 .or. (nstbnd > nstana .and. nstana > 0)) Then
              nstbnd = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdbnd /= bond%bin_pdf)
           k=Merge(1, 0, Abs(bond%rcut-rcb_d) > 1.0e-3_wp)
           Write(message,'(2(a,i10),a,f7.2,a)') &
             'bonds      - collection every ',nstbnd,' step(s); ngrid = ', &
             bond%bin_pdf,' points; cutoff = ',bond%rcut, ' Angs'
           Call info(message,.true.)
           If (i+j+k > 1) Then
             Write(message,'(3(a,i10))') &
               'bonds      - reset values at  ', i,'                  ', j, &
               '                 ', k
             Call info(message,.true.)
           End If
        End If

        If (angle%bin_adf > 0) Then
           If (nstang == 0 .or. (nstang > nstana .and. nstana > 0)) Then
              nstang = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdang /= angle%bin_adf)
           Write(message,'(2(a,i10),a)') &
              'angles     - collection every ',nstang,' step(s); ngrid = ',angle%bin_adf,' points'
           Call info(message,.true.)
           If (i+j > 1) Then
             Write(message,'(2(a,i10))') &
               'angles     - reset values at  ',     i,'                  ',     j
             Call info(message,.true.)
           End If
        End If

        If (dihedral%bin_adf > 0) Then
           If (nstdih == 0 .or. (nstdih > nstana .and. nstana > 0)) Then
              nstdih = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grddih /= dihedral%bin_adf)
           Write(message,'(2(a,i10),a)') &
             'dihedrals  - collection every ',nstdih,' step(s); ngrid = ',dihedral%bin_adf,' points'
           Call info(message,.true.)
           If (i+j > 1) Then
             Write(message,'(2(a,i10))') &
               'dihedrals  - reset values at  ',     i,'                  ',     j
             Call info(message,.true.)
           End If
        End If

        If (inversion%bin_adf > 0) Then
           If (nstinv == 0 .or. (nstinv > nstana .and. nstana > 0)) Then
              nstinv = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdinv /= inversion%bin_adf)
           Write(message,'(2(a,i10),a)') &
             'inversions - collection every ',nstinv,' step(s); ngrid = ',inversion%bin_adf,' points'
           Call info(message,.true.)
           If (i+j > 1) Then
             Write(message,'(2(a,i10))') &
               'inversions - reset values at  ',     i,'                  ',     j
             Call info(message,.true.)
           End If
        End If
     End If

     If (lpana) Then
        Call info('probability distribution analysis printing requested',.true.)
     Else
        Call info('no probability distribution analysis printing requested',.true.)
     End If
  End If

! For safety make them /= 0

  nstbnd=Max(1,nstbnd)
  nstang=Max(1,nstang)
  nstdih=Max(1,nstdih)
  nstinv=Max(1,nstinv)

! report rdf

  If (rdf%l_collect .or. rdf%l_print) Then
     If (rdf%l_collect) Then
        Write(messages(1),'(a)') 'rdf collection requested:'
        Write(messages(2),'(2x,a,i10)') 'rdf collection interval ',rdf%freq
        Write(messages(3),'(2x,a,1p,e12.4)') 'rdf binsize (Angstroms) ',rbin
        Call info(messages,3,.true.)
     Else
        Call info('no rdf%rdf collection requested',.true.)
     End If

     If (rdf%l_print) Then
        Call info('rdf%rdf printing requested',.true.)
     Else
        If (lpana) Then
           Call info('rdf%rdf printing triggered due to a PDA printing request',.true.)
           rdf%l_print=lpana
        Else
           Call info('no rdf%rdf printing requested',.true.)
        End If
     End If

     If (rdf%max_rdf == 0) Then
        Call info('no rdf%rdf pairs specified in FIELD',.true.)
     Else
        Call info('rdf%rdf pairs specified in FIELD',.true.)
     End If

     If ((.not.rdf%l_collect) .or. rdf%max_rdf == 0) Then
        Call info('rdf%rdf routines not to be activated',.true.)
        rdf%l_collect=.false.
        rdf%l_print=.false.
     End If
  End If

! report zden

  If (zdensity%l_collect .or. zdensity%l_print) Then
     If (zdensity%l_collect) Then
        Write(messages(1),'(a)') 'z-density profiles requested:'
        Write(messages(2),'(2x,a,i10)') 'z-density collection interval ',zdensity%frequency
        Write(messages(3),'(2x,a,1p,e12.4)') 'z-density binsize (Angstroms) ',rbin
        Call info(messages,3,.true.)
     Else
        Call info('no z-density profiles requested',.true.)
     End If

     If (zdensity%l_print) Then
        Call info('z-density printing requested',.true.)
     Else
        Call info('no z-density printing requested',.true.)
     End If

     If (.not.zdensity%l_collect) Then
        Call info('z-density routines not to be activated',.true.)
        zdensity%l_print=.false.
     End If
  End If

! report vaf

  If (green%samp > 0 .or. green%l_print) Then
     If (green%samp > 0) Then
        Write(messages(1),'(a)') 'vaf profiles requested:'
        Write(messages(2),'(2x,a,i10)') 'vaf collection interval ',green%freq
        Write(messages(3),'(2x,a,i10)') 'vaf binsize  ',green%binsize
        Call info(messages,3,.true.)
     Else
        Call info('no vaf collection requested',.true.)
     End If

     If (green%l_print) Then
        Call info('vaf printing requested',.true.)
     Else
        Call info('no vaf printing requested',.true.)
     End If

     If (green%l_average) Then
        Call info('time-averaged vaf profile',.true.)
     Else
        Call info('instantaneous vaf profiles',.true.)
     End If
  End If

! report thermal conductivity

!  If (ltcond) Then
!    If (comm%idnode == 0) Then
!      Write(messages(1),'(a)') 'thermal conductivities requested'
!      Write(messages(2),'(a,i10)') 'heat current collection binsize',nsttcond
!      Call info(messages,2,.true.)
!    End If
!  End If

! report data dumping interval, subcelling threshold density and job times

  Write(messages(1),'(a,i10)') 'data dumping interval (steps) ',ndump
  Write(messages(2),'(a,1p,e12.4)') 'subcelling threshold density ',pdplnc
  Write(messages(3),'(a,1p,e12.4)') 'allocated job run time (s) ',tmr%job
  Write(messages(4),'(a,1p,e12.4)') 'allocated job close time (s) ',tmr%clear_screen
  Call info(messages,4,.true.)

! report replay history

  If (.not.lsim) Then
     If (lfce) Then
        Write(messages(1),'(a)') '*** HISTORF will be replayed with full force recalculation ***'
        Write(messages(2),'(a)') '*** There is no actual dynamics/integration!!!             ***'
        Call info(messages,2,.true.)
     Else
        Write(messages(1),'(a)') '*** HISTORY will be replayed (no actual simulation) ***'
        Write(messages(2),'(a)') '*** with structural properties will be recalculated ***'
        Call info(messages,2,.true.)
! abort if there's no structural property to recalculate
        If (.not.(rdf%l_collect .or. zdensity%l_collect .or. dfcts(1)%ldef .or. msd_data%l_msd .or. lrsd .or. (mxgana > 0))) Then
          Call error(580)
        End If
     End If

     If (keyres /= 0) Then
        keyres=0 ! Force clean restart
        Call info('clean start enforced',.true.)
     End If
  End If

!!! RESORT TO DEFAULTS IF NEED BE !!!

  If      (nstrun == 0) Then !!! DRY RUN
     ltemp = .true. ! zero is ok
     lpres = .true. ! zero is ok

     lstep = .true. ! zero is not ok
     If (tstep <= zero_plus) Then
        tstep = 0.001_wp
        Write(message,'(a,1p,e12.4)') 'default simulation timestep (ps) ',tstep
        Call info(message,.true.)
     End If
  Else If (.not.l_str) Then !!! NO STRICT
     If (.not.ltemp) Then ! Simulation temperature
        ltemp=.true.
        thermo%temp=300.0_wp
        Write(message,'(a,1p,e12.4)') 'default simulation temperature (K) ',thermo%temp
        Call info(message,.true.)
     End If

     If (.not.lpres) Then ! Simulation pressure
        lpres=.true.
        thermo%press=0.0_wp
        Write(message,'(a,1p,e12.4)') 'default simulation pressure (katms) ',thermo%press*prsunt
        Call info(message,.true.)
     End If

     If (.not.lstep) Then ! Simulation timestep
        lstep = .true.
        tstep = 0.001_wp
        Write(message,'(a,1p,e12.4)') 'default simulation timestep (ps) ',tstep
        Call info(message,.true.)
     End If

     If ((.not.l_timjob) .and. (.not.l_timcls)) Then ! Job times
        tmr%job=1.0e8_wp
        tmr%clear_screen=1.0e7_wp
        Write(messages(1),'(a,1p,e12.4)') 'allocated job run time (s) ',tmr%job
        Write(messages(1),'(a,1p,e12.4)') 'allocated job close time (s) ',tmr%clear_screen
        Call info(messages,2,.true.)
     Else If ((.not.l_timjob) .and. l_timcls) Then
        tmr%job=100.0_wp*tmr%clear_screen
        Write(message,'(a,1p,e12.4)') 'allocated job run time (s) ',tmr%job
        Call info(message,.true.)
     Else If (l_timjob .and. (.not.l_timcls)) Then
        tmr%clear_screen=0.01_wp*tmr%job
        Write(message,'(a,1p,e12.4)') 'allocated job close time (s) ',tmr%clear_screen
        Call info(message,.true.)
     End If

  End If

  If (l_0 .and. (.not. ltemp)) Then ! zero K over zero fire
     thermo%temp = 10.0_wp
     Write(message,'(a,1p,e12.4)') 'default simulation temperature (K) ',thermo%temp
     Call info(message,.true.)
  End If

  thermo%vel_es2 = thermo%vel_es2 * thermo%vel_es2 ! square of cutoff velocity for inhomogeneous Langevin thermostat and ttm
  amin = Max (amin, 1)        ! minimum number of atoms for ttm ionic temperature cell

!!! ERROR CHECKS !!!
! Temperature

  If ((.not.ltemp) .or. (nstrun > 0 .and. thermo%temp < 1.0_wp)) Call error(380)

! Timestep

  If (.not.lstep) Call error(381)

! check settings in Langevin ensembles

  If ((thermo%ensemble == ENS_NVT_LANGEVIN .or. &
       thermo%ensemble == ENS_NPT_LANGEVIN .or. &
       thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) .and. &
      thermo%chi <= zero_plus) Then
    Call error(462)
  End IF

  If (thermo%ensemble == ENS_NVT_LANGEVIN_INHOMO ) Then
    If (gvar==0 .and. thermo%chi_ep <= zero_plus) Call error(462)
    If (Abs(thermo%chi_es) <= zero_plus) Then
      Call info('assuming no electronic stopping in inhomogeneous Langevin thermostat',.true.)
    End If
  End If

  If ((thermo%ensemble == ENS_NPT_LANGEVIN .or. &
       thermo%ensemble == ENS_NPT_LANGEVIN_ANISO) .and. &
      thermo%tai <= zero_plus) Call error(463)

! check settings in ensembles with thermo%tau_t

  If (Any([ENS_NVT_ANDERSON,ENS_NVT_BERENDSEN,ENS_NVT_NOSE_HOOVER, &
           ENS_NPT_BERENDSEN,ENS_NPT_NOSE_HOOVER,ENS_NPT_MTK, &
           ENS_NPT_BERENDSEN_ANISO,ENS_NPT_NOSE_HOOVER_ANISO,ENS_NPT_MTK_ANISO] &
          == thermo%ensemble) .and. thermo%tau_t <= 0.0_wp) Then
    Call error(464)
  End If

! check settings in ensembles with thermo%press

  If (thermo%variable_cell) Then
     If (.not. thermo%anisotropic_pressure) Then
        If (.not.lpres) Then
           If (lstrext) Then
              thermo%press=(thermo%stress(1)+thermo%stress(5)+thermo%stress(9))/3.0_wp
              thermo%stress=0.0_wp

              Write(messages(1),'(a)') 'tensorial system pressure specified for an npt ensemble simulation'
              Write(messages(2),'(a)') 'scalar pressure derived from pressure tensor as p = Trace[P]/3'
              Write(messages(3),'(a)') 'tensorial pressure to be zeroed (discarded)'
              Write(messages(4),'(a,1p,e12.4)') 'simulation pressure (katms) ',thermo%press*prsunt
              Call info(messages,4,.true.)
           Else
              Call error(387)
           End If
        Else
           If (lstrext) Then
              thermo%stress=0.0_wp

              Write(messages(1),'(a)') 'both tensorial and scalar system pressure specified for an npt ensemble simulation'
              Write(messages(2),'(a)') 'tensorial pressure directive is ignored'
              Write(messages(3),'(a)') 'tensorial pressure to be zeroed (discarded)'
              Call info(messages,3,.true.)
           End If
        End If
     Else If (thermo%anisotropic_pressure) Then
        If (.not.lstrext) Then
           If (.not.lpres) Call error(387)
        Else
           If (lpres) Then
              Write(messages(1),'(a)') 'both tensorial and scalar system pressure specified for an nst ensemble simulation'
              Write(messages(2),'(a)') 'scalar pressure directive is ignored'
              Call info(messages,2,.true.)

! Define new scalar pressure and zero trace pressure tensor

              thermo%press=(thermo%stress(1)+thermo%stress(5)+thermo%stress(9))/3.0_wp
              thermo%stress(1)=thermo%stress(1)-thermo%press
              thermo%stress(5)=thermo%stress(5)-thermo%press
              thermo%stress(9)=thermo%stress(9)-thermo%press
           End If
        End If
     End If
     If (thermo%ensemble /= ENS_NPT_LANGEVIN .and. &
         thermo%ensemble /= ENS_NPT_LANGEVIN_ANISO .and. &
         thermo%tau_p <= 0.0_wp) Then
       Call error(466)
     End If
  End If

! Two-temperature model: calculate atomic density (if not
! already specified and electron-phonon friction
! conversion factor (to calculate thermo%chi_ep from G_ep values)

  If (cellrho<=zero_plus) cellrho = sysrho
  If (cellrho>zero_plus) Then
    rcellrho = 1.0_wp/cellrho
  Else
    rcellrho = 0.0_wp
  End If

  epc_to_chi = 1.0e-12_wp*Jm3K_to_kBA3/3.0_wp
  If (.not. ttmdyndens) epc_to_chi = epc_to_chi*rcellrho

! Check sufficient parameters are specified for TTM electronic specific
! heats, thermal conductivity/diffusivity, energy loss and laser deposition

  If (ttmdyndens) CeType = CeType + 4
  Select Case (CeType)
  Case (0)
  ! constant electronic specific heat: will convert from kB/atom to kB/A^3
  ! by multiplication of atomic density
    Ce0 = Ce0*cellrho
  Case (1)
  ! hyperbolic tangent electronic specific heat: multiplier will be converted
  ! from kB/atom to kB/A^3, temperature term (K^-1) is now scaled by 10^-4
    If (Abs(sh_A) <= zero_plus .or. Abs(sh_B) <= zero_plus) Call error(671)
    sh_A = sh_A*cellrho
    sh_B = sh_B*1.0e-4_wp
  Case (2)
  ! linear electronic specific heat to Fermi temperature: maximum
  ! value will be converted from kB/atom to kB/A^3
    If (Abs(Tfermi) <= zero_plus .or. Abs(Cemax) <= zero_plus) Call error(671)
    Cemax = Cemax*cellrho
  Case (4)
  ! constant electronic specific heat: will convert from kB/atom to kB/A^3
  ! by multiplication of atomic density
  Case (5)
  ! hyperbolic tangent electronic specific heat: multiplier will be converted
  ! from kB/atom to kB/A^3, temperature term (K^-1) is now scaled by 10^-4
    If (Abs(sh_A) <= zero_plus .or. Abs(sh_B) <= zero_plus) Call error(671)
    sh_B = sh_B*1.0e-4_wp
  Case (6)
  ! linear electronic specific heat to Fermi temperature: maximum
  ! value will be converted from kB/atom to kB/A^3
    If (Abs(Tfermi) <= zero_plus .or. Abs(Cemax) <= zero_plus) Call error(671)
  End Select

  Select Case (KeType)
  ! constant and Drude thermal conductivity: convert from W m^-1 K^-1
  ! to kB ps^-1 A^-1
  Case (1,2)
    If (isMetal .and. Abs(Ka0) <= zero_plus) Call error(672)
    Ka0 = Ka0*JKms_to_kBAps
  End Select

  Select Case (DeType)
  Case (1)
  ! constant thermal diffusivity: convert from m^2 s^-1 to A^2 ps^-1
    If (.not. isMetal .and. Abs(Diff0) <= zero_plus) Call error(673)
    Diff0 = Diff0*1.0e8_wp
  Case (2)
  ! reciprocal thermal diffusivity: convert from m^2 s^-1 to A^2 ps^-1
  ! and Diff0 scaled with system temperature
    If (.not. isMetal .and. Abs(Diff0) <= zero_plus .or. Abs(Tfermi) <= zero_plus) Call error(673)
    Diff0 = Diff0*thermo%temp*1.0e8_wp
  End Select

  ! spatial deposition (gaussian) standard deviation: convert from nm to A
  sig = sig*10.0_wp

  ! penetration depth: convert from nm to A
  If ((sdepoType == 2 .and. Abs(dEdX) <= zero_plus .and. Abs(pdepth-1.0_wp)<=zero_plus) .or. &
      (sdepoType == 3 .and. Abs(pdepth-1.0_wp) <= zero_plus)) &
  Call warning(510,0.0_wp,0.0_wp,0.0_wp)
  pdepth = 10.0_wp*pdepth


  ! electronic stopping power: convert from eV/nm to eV/A
  ! fluence: convert from mJ cm^-2 to eV A^-2
  If (sdepoType>0 .and. sdepoType<3 .and. (Abs(dEdx) <= zero_plus .or. Abs(fluence)<zero_plus)) &
  Call warning(515,0.0_wp,0.0_wp,0.0_wp)

  fluence = fluence*mJcm2_to_eVA2
  dEdX = 0.1_wp*dEdX

End Subroutine read_control

Subroutine scan_control                                    &
           (rcter, &
           mxrgd,imcon,imc_n,cell,xhi,yhi,zhi,             &
           mxgana,         &
           l_str,lsim,l_vv,l_n_e,l_n_r,lzdn,l_n_v,l_ind,   &
           rbin,                          &
           nstfce,mxspl,alpha,kmaxa1,kmaxb1,kmaxc1,cshell,stats,  &
           thermo,green,devel,msd_data,met,pois,bond,angle, &
           dihedral,inversion,zdensity,neigh,vdw,tersoff,rdf,mpole,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of the control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2017
! contrib   - i.j.bush february 2014
! contrib   - a.v.brukhno & i.t.todorov april 2014 (intramolecular TPs & PDFs)
! contrib   - m.a.seaton june 2014 (VAF)
! contrib   - p.s.petkov february 2015
! contrib   - a.m.elena february 2017
! contrib   - m.a.seaton march 2017 (TTM)
! contrib   - a.b.g.chalk march 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( InOut ) :: l_n_e
  Logical,           Intent(   Out ) :: l_str,lsim,l_vv,l_n_r,lzdn,l_n_v,l_ind
  Integer,           Intent( In    ) :: mxrgd,imcon
  Integer,           Intent( InOut ) :: imc_n
  Integer,           Intent(   Out ) :: mxgana, &
                                        nstfce,mxspl,kmaxa1,kmaxb1,kmaxc1
  Real( Kind = wp ), Intent( In    ) :: xhi,yhi,zhi,rcter
  Real( Kind = wp ), Intent( InOut ) :: cell(1:9)
  Real( Kind = wp ), Intent(   Out ) :: rbin,alpha
  Type( core_shell_type ), Intent (   In  )   :: cshell
  Type( stats_type ), Intent( InOut ) :: stats
  Type( thermostat_type ), Intent( InOut ) :: thermo
  Type( development_type ), Intent( InOut ) :: devel
  Type( greenkubo_type ), Intent( InOut ) :: green
  Type( msd_type ), Intent( InOut ) :: msd_data
  Type( metal_type ), Intent( InOut ) :: met
  Type( poisson_type ), Intent( InOut ) :: pois
  Type( bonds_type ), Intent( InOut ) :: bond
  Type( angles_type ), Intent( InOut ) :: angle
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( inversions_type ), Intent( InOut ) :: inversion
  Type( z_density_type ), Intent( InOut ) :: zdensity
  Type( neighbours_type ), Intent( InOut ) :: neigh
  Type( vdw_type ), Intent( InOut ) :: vdw
  Type( tersoff_type ), Intent( In    )  :: tersoff
  Type( rdf_type ), Intent( InOut ) :: rdf
  Type( mpole_type ), Intent( InOut ) :: mpole
  Type( comms_type ), Intent( InOut ) :: comm

  Logical                :: carry,safe,la_ana,la_bnd,la_ang,la_dih,la_inv, &
                            lrcut,lrpad,lrvdw,lrmet,lelec,lvdw,lmet,l_n_m,lter,l_exp
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,word1,word2,akey
  Integer                :: i,itmp,nstrun,mxspl2
  Real( Kind = wp )      :: celprp(1:10),cut,eps0,fac,tol,tol1

  Integer,           Parameter :: mxspl_def = 8,        & ! default spline for SPME (4 & 6 possible)
                                  mxspl_min = 3           ! minimum spline order, needed for derivatives of forces
  Real( Kind = wp ), Parameter :: rcut_def  = 1.0_wp  , & ! minimum real space cutoff
                                  rbin_def  = 0.05_wp , & ! default bin size (RDF/USR & z-density)
                                  rcbnd_def = 2.5_wp      ! minimum bond length for bond analysis

! default reading indices options

  l_ind=.true.

! strict flag

  l_str = .true.

! replay history option (real dynamics = no replay)

  lsim = .true. ! don't replay history

! slab option default

  imc_n = imcon

! integration flavour - velocity verlet assumed

  l_vv = .true.

! default switches for intramolecular analysis grids

  la_ana = .false. ; mxgana  = 0
  la_bnd = .false. ; bond%bin_pdf = 0
  la_ang = .false. ; angle%bin_adf = 0
  la_dih = .false. ; dihedral%bin_adf = 0
  la_inv = .false. ; inversion%bin_adf = 0

! electrostatics and no electrostatics, rdf and no rdf, vdw and no vdw,
! metal and no metal, tersoff and no tersoff interactions,
! cutoff and padding, and binsize defaults

  lelec = .false.
! l_n_e is now first determined in scan_field l_n_e = (.false.)

  rdf%l_collect  = (rdf%max_rdf > 0)
  l_n_r = .not.rdf%l_collect

  lvdw  = (vdw%max_vdw > 0)
  l_n_v = .false.
  lrvdw = .false. ! Even though it vdw%cutoff may have been read from TABLE

  lmet  = (met%max_metal > 0)
  l_n_m = .not.lmet
  lrmet = (met%rcut > 1.0e-6_wp)

  lter  = (tersoff%max_ter > 0)

  lrcut = .false.
  neigh%cutoff  = 0.0_wp

  lrpad = .false.
  neigh%padding  = 0.0_wp

  rbin  = rbin_def

! Frequency of the SPME k-space evaluation

  nstfce = -1 ! None defined

! Ewald/Poisson Solver sum parameters defaults

  mxspl = 0
  alpha = 0.0_wp
  kmaxa1 = 0
  kmaxb1 = 0
  kmaxc1 = 0

! default number of steps and expansion option

  nstrun = 0
  l_exp = .false.

! default stack size

  stats%mxstak = 1

! default switch for two-temperature model (ttm)

  l_ttm = .false.

! default switch for redistribution of energy from
! deactivated electronic cells for ttm

  redistribute = .false.

! default switches for removing centre-of-mass motion
! when applying inhomogeneous Langevin thermostat with ttm

  ttmthvel  = .true.
  ttmthvelz = .false.

! default values for ttm ionic and electronic voxel grid sizes

  ntsys(3)  = 10
  eltsys(1) = 50
  eltsys(2) = 50
  eltsys(3) = 50

! default ttm heat capacity, conductivity and diffusivity types

  CeType = 0
  KeType = 0
  DeType = 0

! default ttm electron-phonon coupling type

  gvar = 0

! Set safe flag

  safe=.true.

! Open the simulation input file

  If (comm%idnode == 0) Inquire(File=Trim(control), Exist=safe)
  Call gcheck(comm,safe,"enforce")
  If (.not.safe) Then
     Go To 10
  Else
     If (comm%idnode == 0) Open(Unit=nread, File=Trim(control), Status='old')
  End If

! First Pass.  Get cutoff distances, stacksize and density variation.

  Call get_line(safe,nread,record,comm)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record,comm)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! record is commented out

     If      (word(1:1) == '#' .or. word(1:1) == ' ') Then

! read slab option (limiting DD slicing in z direction to 2)

     Else If (word(1:4) == 'slab') Then

        If (imcon /= 0 .and. imcon /= 6) imc_n=6

! read real space cut off

     Else If (word(1:3) == 'cut' .or. word(1:4) == 'rcut') Then

        lrcut = .true.
        Call get_word(record,word)
        neigh%cutoff = Abs(word_2_real(word))
        lrcut = (neigh%cutoff > zero_plus) ! if zero or nothing is entered

! read real space cut off

     Else If (word(1:3) == 'pad' .or. word(1:4) == 'rpad') Then

        lrpad = .true.
        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        neigh%padding = Max(neigh%padding,Abs(word_2_real(word)))
        lrpad = (neigh%padding > zero_plus) ! if zero or nothing is entered

! read vdw cutoff

     Else If (word(1:4) == 'rvdw') Then

        lrvdw=.true.
        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        If (vdw%cutoff > 1.0e-6_wp) Then
           vdw%cutoff = Min(vdw%cutoff,word_2_real(word))
        Else
           vdw%cutoff = Abs(word_2_real(word))
        End If
        lrvdw = (vdw%cutoff > zero_plus) ! if zero or nothing is entered

! read binsize option

     Else If (word(1:7) == 'binsize') Then

        Call get_word(record,word)
        rbin = Abs(word_2_real(word))
        zdensity%bin_width = rbin

! read dpd ensembles option

     Else If (word(1:8) == 'ensemble') Then

        Call get_word(record,word)
        If (word(1:3) == 'nvt') Then
           Call get_word(record,word)
           If (word(1:3) == 'dpd') Then
              If      (word(1:5) == 'dpds1') Then
                 thermo%key_dpd = 1
              Else If (word(1:5) == 'dpds2') Then
                 thermo%key_dpd = 2
              End If
           End If
        End If

! read replay history option

     Else If (word(1:6) == 'replay') Then

        lsim = .false.

! read number of timesteps

     Else If (word(1:5) == 'steps') Then

        Call get_word(record,word)
        nstrun = Nint(Abs(word_2_real(word)))

! read expansion option

     Else If (word(1:5) == 'nfold') Then

        l_exp = .true.

! read stack size

     Else If (word(1:5) == 'stack') Then

        Call get_word(record,word)
        If (word(1:4) == 'size') Call get_word(record,word)

        stats%mxstak = Nint(Abs(word_2_real(word)))

! read MSD option

     Else If (word(1:6) == 'msdtmp') Then

        msd_data%l_msd = .true.

! read VAF option and sample frequency and binsize - defaults in greenkubo_module

     Else If (word(1:3) == 'vaf') Then

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        green%freq = Abs(Nint(word_2_real(word,0.0_wp)))
        If (green%freq == 0) green%freq=50

        Call get_word(record,word)
        If (word(1:3) == 'bin' .or. word(1:5) == 'size') Call get_word(record,word)
        If (word(1:3) == 'bin' .or. word(1:5) == 'size') Call get_word(record,word)
        green%binsize = Abs(Nint(word_2_real(word,0.0_wp)))

        If (green%binsize == 0) green%binsize=Merge(2*green%freq,100,green%freq >= 100)
        green%samp = Ceiling(Real(green%binsize,wp)/Real(green%freq,wp))

! read DL_POLY_2/Classic delr Verlet shell strip cutoff option (compatibility)
! as DL_POLY_4 real space cutoff padding option

     Else If (word(1:4) == 'delr') Then

        lrpad = .true.
        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        neigh%padding = Max(neigh%padding,0.25_wp*Abs(word_2_real(word)))
        lrpad = (neigh%padding > zero_plus) ! if zero or nothing is entered

! read DL_POLY_2/Classic multiple timestep option (compatibility)
! as DL_POLY_4 infrequent k-space SPME evaluation option

     Else If (word(1:4) == 'mult') Then

        Call get_word(record,word)
        If (word(1:5) == 'times' .or. word(1:4) == 'step') Call get_word(record,word)
        nstfce=Max(nstfce,Nint(Abs(word_2_real(word))))

! read electrostatics

     Else If (word(1:5) == 'ewald' .or. word(1:4) == 'spme') Then

        Call get_word(record,word)

        If (word(1:5) == 'evalu') Then

! infrequent k-space SPME evaluation

           Call get_word(record,word)
           If (word(1:5) == 'every') Call get_word(record,word)
           nstfce=Max(nstfce,Nint(Abs(word_2_real(word))))

        Else

           lelec = .true.

        End If

     Else If (word(1:5) == 'poiss' .or. word(1:5) == 'psolv') Then

        lelec = .true.

     Else If (word(1:6) == 'distan') Then

        lelec = .true.

     Else If (word(1:4) == 'coul') Then

        lelec = .true.

     Else If (word(1:5) == 'shift') Then

        lelec = .true.

     Else If (word(1:8) == 'reaction') Then

        lelec = .true.

! read "no vdw", "no elec" and "no str" options

     Else If (word(1:5) == 'polar') Then

        Call get_word(record,word)
        If (word(1:6) == 'scheme' .or. word(1:4) == 'type') Call get_word(record,word)
        If (word(1:6) == 'scheme' .or. word(1:4) == 'type') Call get_word(record,word)
        If (word(1:6) == 'charmm' .and. cshell%mxshl > 0) mpole%key=POLARISATION_CHARMM

     Else If (word(1:2) == 'no') Then

        Call get_word(record,word)

        If (word(1:3) == 'vdw') Then

           l_n_v = .true.

        Else If (word(1:4) == 'elec') Then

           l_n_e = .true.

        Else If (word(1:3) == 'ind' ) Then

           l_ind=.false.
           Call info('no index (reading in CONFIG) option on',.true.)

        Else If (word(1:3) == 'str' ) Then

           l_str=.false.

        End If

! read integration flavour

     Else If (word(1:8) == 'integrat') Then

        Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:8) == 'leapfrog') l_vv=.false.

! read analysis (intramolecular distribution calculation) option

     Else If (word(1:3) == 'ana') Then

        la_ana = .true.

        Call get_word(record,word)
        akey = word(1:3)

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)

        Call get_word(record,word)
        If (word(1:5) == 'nbins' .or. word(1:5) == 'ngrid' .or. word(1:4) == 'grid') Call get_word(record,word)

        If      (akey == 'all') Then
           la_bnd = .true.
           la_ang = .true.
           la_dih = .true.
           la_inv = .true.

           mxgana = Abs(Nint(word_2_real(word)))
           bond%bin_pdf = Max(bond%bin_pdf,mxgana)
           angle%bin_adf = Max(angle%bin_adf,mxgana)
           dihedral%bin_adf = Max(dihedral%bin_adf,mxgana)
           inversion%bin_adf = Max(inversion%bin_adf,mxgana)

           Call get_word(record,word) ! AB: for "rbnd"/"rmax"/"max"/figure
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           bond%rcut=Max(bond%rcut,word_2_real(word,0.0_wp))
        Else If (akey == 'bon') Then
           la_bnd = .true.

           bond%bin_pdf = Max(bond%bin_pdf,Abs(Nint(word_2_real(word))))

           Call get_word(record,word) ! AB: for "rbnd"/"rmax"/"max"/figure
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           bond%rcut=Max(bond%rcut,word_2_real(word,0.0_wp))
        Else If (akey == 'ang') Then
           la_ang = .true.

           angle%bin_adf = Max(angle%bin_adf,Abs(Nint(word_2_real(word))))
        Else If (akey == 'dih') Then
           la_dih = .true.

           dihedral%bin_adf = Max(dihedral%bin_adf,Abs(Nint(word_2_real(word))))
        Else If (akey == 'inv') Then
           la_inv = .true.

           inversion%bin_adf = Max(inversion%bin_adf,Abs(Nint(word_2_real(word))))
        End If

! read rdf calculation option

     Else If (word(1:3) == 'rdf') Then

        rdf%l_collect = .true.
        l_n_r = .not.rdf%l_collect

! read z-density profile option

     Else If (word(1:4) == 'zden') Then

        lzdn = .true.

! read two-temperature model options

     Else If (word(1:3) == 'ttm') Then

        l_ttm = .true.

        Call get_word(record,word)

        If (word(1:4) == 'ncit') Then

        ! number of coarse-grained ion temperature cells (CIT)
        ! in z-direction: geometry of system determines 
        ! CITs in x- and y-directions

          Call get_word(record,word)
          ntsys(3) = Abs(Nint(word_2_real(word)))

        Else If (word(1:4) == 'ncet') Then

        ! number of coarse-grained electronic temperature cells
        ! (CET) in x-, y- and z-directions

          Call get_word(record,word)
          eltsys(1) = Abs(Nint(word_2_real(word)))
          Call get_word(record,word)
          eltsys(2) = Abs(Nint(word_2_real(word)))
          Call get_word(record,word)
          eltsys(3) = Abs(Nint(word_2_real(word)))

        Else If (word1(1:5) == 'metal') Then

        ! sets properties of electronic subsystem as a metal

          isMetal = .true.

        Else If (word1(1:8) == 'nonmetal') Then

        ! sets properties of electronic subsystem as a non-metal

          isMetal = .false.

        Else If (word(1:7) == 'ceconst') Then

        ! electronic specific heat capacity given as constant value

          CeType = 0

        Else If (word(1:6) == 'cetanh') Then

        ! electronic specific heat capacity given as tanh function

          CeType = 1

        Else If (word(1:5) == 'celin') Then

        ! electronic specific heat capacity given as linear function
        ! up to Fermi temperature, constant afterwards

          CeType = 2

        Else If (word(1:5) == 'cetab') Then

        ! electronic volumetric heat capacity given in tabulated form

          CeType = 3

        Else If (word(1:5) == 'keinf') Then

        ! infinite electronic thermal conductivity

          DeType = 0
          KeType = 0
          isMetal = .true.

        Else If (word(1:7) == "keconst") Then

        ! electronic thermal conductivity given as constant value

          DeType = 0
          KeType = 1
          isMetal = .true.

        Else If (word(1:7) == 'kedrude') Then

        ! electronic thermal conductivity given as drude model (propertional to
        ! electronic temperature, giving t.c. at system temperature)

          DeType = 0
          KeType = 2
          isMetal = .true.

        Else If (word(1:5) == 'ketab') Then

        ! electronic thermal conductivity given in tabulated form

          DeType = 0
          KeType = 3
          isMetal = .true.

        Else If (word(1:4) == 'diff' .or. word(1:7)=='deconst') Then

        ! electronic thermal diffusivity given as constant value
        ! (for non-metal systems)

          KeType = 1
          DeType = 1
          isMetal = .false.

        Else If (word(1:7) == 'derecip') Then

        ! electronic thermal diffusivity given as reciprocal function
        ! of temperature (up to Fermi temperature), constant afterwards

          KeType = 1
          DeType = 2
          isMetal = .false.

        Else If (word(1:4) == 'detab') Then

        ! electronic thermal diffusivity given in tabulated form

          KeType = 1
          DeType = 3
          isMetal = .false.

        Else If (word(1:4) == 'varg') Then

        ! variable electron-phonon coupling constant (thermo%chi_ep) based on
        ! tabular electronic stopping terms (in g.dat file): option to
        ! apply value homogeneously across system (based on average 
        ! electronic temperature) or heterogeneously (using local 
        ! electronic temperature for each voxel)

          Call get_word(record,word)
          If (word(1:4) == 'homo') Then
            gvar = 1
          Else If (word(1:6) == 'hetero') Then
            gvar = 2
          End If

        Else If (word(1:7) == 'nothvel') Then

        ! option to switch off centre-of-mass motion correction to
        ! velocities used in inhomogeneous Langevin thermostat

          ttmthvel = .false.
          ttmthvelz = .false.

        Else If (word(1:6) == 'thvelz') Then

        ! option to only apply centre-of-mass motion corrections to
        ! z components of velocities, when used in inhomogeneous
        ! Langevin thermostat

          ttmthvelz = .true.
          ttmthvel = .true.

        Else If (word(1:6) == 'redist') Then

        ! redistribute electronic energy from any deactivated cells
        ! to active neighbouring cells: note that this requires overlap
        ! of electronic and ionic temperature grids in at least
        ! x- and y-directions

          redistribute = .true.

        End If

! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     Else If (word(1:6) == 'errors') Then
       Call lower_case(record)
       Call get_word(record,word)
       If(word(1:4) == 'jack') Then
          rdf%l_errors_jack = .TRUE.
       Else
          rdf%l_errors_block = .TRUE.
       End If
     End If

  End Do

! in case of intramolecular distribution analysis

  If (la_ana) Then
     If (mxgana > 0) Then
        bond%bin_pdf = Max(bond%bin_pdf,mxgana)
        angle%bin_adf = Max(angle%bin_adf,mxgana)
        dihedral%bin_adf = Max(dihedral%bin_adf,mxgana)
        inversion%bin_adf = Max(inversion%bin_adf,mxgana)
     End If

! switch indicators for set_bounds

     If (la_bnd) Then
        If (bond%bin_pdf == 0) bond%bin_pdf = -1
        bond%rcut=Max(bond%rcut,rcbnd_def)
     End If
     If (la_ang .and. angle%bin_adf == 0) angle%bin_adf = -1
     If (la_dih .and. dihedral%bin_adf == 0) dihedral%bin_adf = -1
     If (la_inv .and. inversion%bin_adf == 0) inversion%bin_adf = -1

! mxgana by construction equals the largest possible grid
! or 1 (positive) as an indicator for analysis

     mxgana=Max(1,bond%bin_pdf,angle%bin_adf,dihedral%bin_adf,inversion%bin_adf)
  End If

! Sort electrostatics

  If (lelec) Then
     If (l_n_e) lelec = .not.l_n_e
  Else
     l_n_e = .true.
  End If

! reinitialise multipolar electrostatics indicators

  If (l_n_e) Then
     mpole%max_mpoles = 0
     mpole%max_order = 0
     mpole%key = POLARISATION_DEFAULT
  End If

! Sort vdw

  If (lvdw) Then
     If (.not.lrvdw) Then
        lrvdw = (vdw%cutoff > 1.0e-6_wp)
        vdw%cutoff = Min(vdw%cutoff,Max(neigh%cutoff,rcut_def))
     End If

     If (l_n_v) lvdw = .not.l_n_v
  Else
     l_n_v = .true.
  End If

! Sort neigh%cutoff as the maximum of all valid cutoffs

  neigh%cutoff=Max(neigh%cutoff,vdw%cutoff,met%rcut,rkim,2.0_wp*Max(rcter,bond%rcut)+1.0e-6_wp)

  If (comm%idnode == 0) Rewind(nread)

! Second Pass.  Sort out cutoffs, cell parameters and Ewald/Poisson Solver precision.

  Call get_line(safe,nread,record,comm)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record,comm)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! record is commented out

     If      (word(1:1) == '#' .or. word(1:1) == ' ') Then

     Else If (lelec .and. ((word(1:5) == 'ewald' .or. word(1:4) == 'spme') .or. &
                           (word(1:5) == 'poiss' .or. word(1:5) == 'psolv'))) Then

! Double the kmax size if specified "ewald sum" instead of "spme sum"

        itmp=0
        If       (word(1:5) == 'ewald') Then
           itmp=2
        Else If  (word(1:4) == 'spme' ) Then
           itmp=1
        End If

        Call get_word(record,word)

        If      (word(1:5) == 'evalu')     Then

        Else

! neigh%cutoff MUST be >= rcut_def

           If (neigh%cutoff < rcut_def) neigh%cutoff=rcut_def

! define cut

           cut=neigh%cutoff+1.0e-6_wp

! fix cell vectors for image conditions with discontinuities

           If (imcon == 0) Then

              cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
              cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

              cell(2) = 0.0_wp
              cell(3) = 0.0_wp
              cell(4) = 0.0_wp
              cell(6) = 0.0_wp
              cell(7) = 0.0_wp
              cell(8) = 0.0_wp

           Else If (imcon == 6) Then

              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

           End If

           If (itmp > 0) Then ! Ewald or SPME

              If (word(1:9) == 'precision') Then

                 Call dcell(cell,celprp)

                 Call get_word(record,word)
                 eps0 = Abs(word_2_real(word))
                 eps0 = Max(Min(eps0,0.5_wp),1.0e-20_wp)

                 Call get_word(record,word)
                 mxspl = Abs(Nint(word_2_real(word)))

                 tol = Sqrt(Abs(Log(eps0*neigh%cutoff)))
                 alpha = Sqrt(Abs(Log(eps0*neigh%cutoff*tol)))/neigh%cutoff
                 tol1 = Sqrt(-Log(eps0*neigh%cutoff*(2.0_wp*tol*alpha)**2))

                 fac = 1.0_wp
                 If (imcon == 4 .or. imcon == 5 .or. imcon == 7) fac = 2.0_wp**(1.0_wp/3.0_wp)

                 kmaxa1 = 2*Nint(0.25_wp + fac*celprp(7)*alpha*tol1/pi)
                 kmaxb1 = 2*Nint(0.25_wp + fac*celprp(8)*alpha*tol1/pi)
                 kmaxc1 = 2*Nint(0.25_wp + fac*celprp(9)*alpha*tol1/pi)

! neigh%cutoff is needed directly for the SPME and it MUST exist

                 If (.not.lrcut) Call error(433)

              Else

                 If (word(1:3) == 'sum') Call get_word(record,word)
                 alpha = Abs(word_2_real(word))

                 Call get_word(record,word)
                 kmaxa1 = itmp*Nint(Abs(word_2_real(word)))

                 Call get_word(record,word)
                 kmaxb1 = itmp*Nint(Abs(word_2_real(word)))

                 Call get_word(record,word)
                 kmaxc1 = itmp*Nint(Abs(word_2_real(word)))

                 Call get_word(record,word)
                 mxspl = Nint(Abs(word_2_real(word)))

! Sanity check for ill defined ewald sum parameters 1/8*2*2*2 == 1

                 tol=alpha*Real(kmaxa1,wp)*Real(kmaxa1,wp)*Real(kmaxa1,wp)
                 If (Nint(tol) < 1) Call error(9)

! neigh%cutoff is not needed directly for the SPME but it's needed
! for the link-cell division of the domains
! let's not fail here if no cutoff is specified

              End If

! Get default spline order or one driven by multipolar sums if none is specified
! Only even order splines are allowed so pick the even=odd+1 if resulting in odd!!!

              If (mxspl == 0) Then
                 mxspl  = mxspl_def+mpole%max_order
                 mxspl2 = mxspl
              Else
                 mxspl  = Max(mxspl,mxspl_min)
                 mxspl2 = mxspl+mpole%max_order
                 mxspl2 = 2*Ceiling(0.5_wp*Real(mxspl2,wp))
              End If
              mxspl=2*Ceiling(0.5_wp*Real(mxspl,wp))
              mxspl=Max(mxspl,mxspl2)

           Else !If (itmp == 0) Then ! Poisson Solver

              Do i=1,4
                 If (word(1:5) == 'delta') Then   ! spacing
                    Call get_word(record,word)
                    alpha=1.0_wp/Abs(word_2_real(word))
                 End If

                 If (word(1:3) == 'eps') Then     ! tolerance
                    Call get_word(record,word)
                    pois%eps=Abs(word_2_real(word))
                 End If

                 If (word(1:6) == 'maxits') Then  ! max number of iteration
                    Call get_word(record,word)
                    pois%mxitcg=Nint(Abs(word_2_real(word)))
                 End If

                 If (word(1:7) == 'jmaxits') Then ! max number Jacobian iterations
                    Call get_word(record,word)
                    pois%mxitjb=Nint(Abs(word_2_real(word)))
                 End If

                 Call get_word(record,word)
              End Do

! Check for undefined and ill defined parameters

!             0.1 Angs <= delta=1/alpha <= Min(3 Angs,neigh%cutoff/3) - 3 grid points within a link-cell
              If (alpha > 10.0_wp)                       alpha = 10.0_wp
              If (alpha < 1.0_wp/Min(3.0_wp,cut/3.0_wp)) alpha = 1.0_wp/Min(3.0_wp,cut/3.0_wp)

              If (pois%mxitcg == 0) pois%mxitcg = 1000 ! default
              If (pois%mxitjb == 0) pois%mxitjb = 1000 ! default

! Derive grid spacing represented as a k-vector

              Call dcell(cell,celprp)
              kmaxa1 = Nint(celprp(7)*alpha)
              kmaxb1 = Nint(celprp(8)*alpha)
              kmaxc1 = Nint(celprp(9)*alpha)

! Define stencil halo size (7) as a spline order (7/2==3)

              mxspl = 3

           End If

        End If

     Else If (word(1:6) == 'finish') Then

! Sort vdw%cutoff

        If ((.not.lrvdw) .and. lvdw) Then
           If (lrcut) Then
              lrvdw=.true.
              vdw%cutoff=neigh%cutoff
           Else
              Call error(402)
           End If
        End If

! Sort met%rcut

        If ((.not.lrmet) .and. lmet) Then
           If (lrcut .or. lrvdw) Then
              lrmet=.true.
              met%rcut=Max(neigh%cutoff,vdw%cutoff)
           Else
              Call error(382)
           End If
        End If

! Sort neigh%cutoff by a reset sequence
! neigh%cutoff may be >= rcut_def but lrcut may still be .false.
! mxspl = 0 is an indicator for no SPME or Poisson Solver electrostatics in CONTROL

        If (mxspl /= 0) Then ! SPME or Poisson Solver

! (1) to Max(neigh%cutoff,Max(cell_width*mxspl/kmax),mxspl*delta) satisfying SPME b-splines
! propagation width or the Poisson Solver extra halo relation to cutoff
! delta=1/alpha is the grid spacing and mxspl is the grid length needed for the
! 3 haloed stencil of differentiation

           If (.not.lrcut) Then
              lrcut=.true.

              Call dcell(cell,celprp)
              neigh%cutoff=Max( neigh%cutoff, Merge(Real(mxspl,wp)/alpha,                          &
                                    Max(celprp(7)*Real(mxspl,wp)/Real(kmaxa1,wp),  &
                                        celprp(8)*Real(mxspl,wp)/Real(kmaxb1,wp),  &
                                        celprp(9)*Real(mxspl,wp)/Real(kmaxc1,wp)), &
                                    itmp == 0) )
           End If

! Reset vdw%cutoff, met%rcut and neigh%cutoff when only tersoff potentials are opted for

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r) Then
              vdw%cutoff=0.0_wp
              met%rcut=0.0_wp
              If (.not.l_str) Then
                 If (mxrgd == 0) Then ! compensate for Max(Size(RBs))>vdw%cutoff
                    neigh%cutoff=2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp
                 Else
                    neigh%cutoff=Max(neigh%cutoff,2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp)
                 End If
              End If
           End If

        Else

! no SPME electrostatics is specified but neigh%cutoff is still needed for
! domain decompositioning and link-celling
! It is needed for the rest of the types of electrostatics

           If ((.not.lrcut) .and. lelec) Call error(382)

! So there is neigh%cutoff and some kind of electrostatics(-: or neither

! Reset neigh%cutoff to something sensible if sensible is an option

           If ( ((.not.lrcut) .or. (.not.l_str)) .and. &
                (lrvdw .or. lrmet .or. lter .or. kimim /= ' ') ) Then
              lrcut=.true.
              If (mxrgd == 0) Then ! compensate for Max(Size(RBs))>vdw%cutoff
                 neigh%cutoff=Max(vdw%cutoff,met%rcut,rkim,2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp)
              Else
                 neigh%cutoff=Max(neigh%cutoff,vdw%cutoff,met%rcut,rkim,2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp)
              End If
           End If

! Reset vdw%cutoff and met%rcut when only tersoff potentials are opted for and
! possibly reset neigh%cutoff to 2.0_wp*rcter+1.0e-6_wp (leaving room for failure)

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r .and. kimim == ' ') Then
              vdw%cutoff=0.0_wp
              met%rcut=0.0_wp
              If (.not.l_str) Then
                 lrcut=.true.
                 If (mxrgd == 0) Then ! compensate for Max(Size(RBs))>vdw%cutoff
                    neigh%cutoff=2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp
                 Else
                    neigh%cutoff=Max(neigh%cutoff,2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp)
                 End If
              End If
           End If

! neigh%cutoff must exist

           If (.not.lrcut) Call error(382)

! define cut

           cut=neigh%cutoff+1.0e-6_wp

! fix cell vectors for image conditions with discontinuities

           If (imcon == 0) Then

              cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
              cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

              cell(2) = 0.0_wp
              cell(3) = 0.0_wp
              cell(4) = 0.0_wp
              cell(6) = 0.0_wp
              cell(7) = 0.0_wp
              cell(8) = 0.0_wp

           Else If (imcon == 6) Then

              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

           End If

        End If

! Sort met%rcut=neigh%cutoff if metal interactions are in play, even if
! they are defined by EAM since met%rcut can be /= neigh%cutoff in such
! instances, this can break the NLAST check in metal_ld_set_halo

        If (lmet) met%rcut = neigh%cutoff

! Sort vdw%cutoff=neigh%cutoff if VDW interactions are in play

        If (lvdw .and. vdw%cutoff > neigh%cutoff) vdw%cutoff = neigh%cutoff

! Sort rbin as now neigh%cutoff is already pinned down

        If (rbin < 1.0e-05_wp .or. rbin > neigh%cutoff/4.0_wp) rbin = Min(rbin_def,neigh%cutoff/4.0_wp)

        carry=.false.

     End If

  End Do

  If (comm%idnode == 0) Close(Unit=nread)

! Enforce VV for DPD thermostat

  If (thermo%key_dpd > 0) l_vv = .true.

! When not having dynamics or prepared to terminate
! expanding and not running the small system prepare to exit gracefully

  devel%l_trm = (l_exp .and. nstrun == 0)
  If (((.not.lsim) .or. devel%l_trm) .and. lrpad) neigh%padding=0.0_wp

  rdf%l_errors_block = rdf%l_errors_block .and. rdf%l_collect
  rdf%l_errors_jack = rdf%l_errors_jack .and. rdf%l_collect
  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control

Subroutine scan_control_pre(imc_n,dvar,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning the imc_n & dvar options in the
! control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( InOut ) :: imc_n
  Real( Kind = wp ), Intent(   Out ) :: dvar
  Type( comms_type ), Intent( InOut ) :: comm

  Logical                :: carry,safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word

  safe   = .true.  ! all is safe

! density variation parameter default

  dvar = 1.0_wp

  If (comm%idnode == 0) Inquire(File=Trim(control), Exist=safe)
  Call gcheck(comm,safe,"enforce")
  If (.not.safe) Then
     Go To 10
  Else
     If (comm%idnode == 0) Open(Unit=nread, File=Trim(control), Status='old')
  End If

! Read TITLE record

  Call get_line(safe,nread,record,comm)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record,comm)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! read density variation option
! this is really a pre-scan in order to get the MD box dimensions
! from scan_config giving it failure estimates for when reading

     If      (word(1:7) == 'densvar') Then

        Call get_word(record,word)
        dvar = Abs(word_2_real(word))
        dvar = 1.0_wp + Abs(dvar)/100.0_wp

! read slab option
! limiting DD slicing in z direction to 2 for load balancing purposes
! this is really a pre-scan in order to get the MD box dimensions
! from scan_config before the option is read again in scan_control

     Else If (word(1:4) == 'slab') Then

        imc_n=6

! io options

! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     End If

  End Do

  If (comm%idnode == 0) Close(Unit=nread)

  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control_pre

Subroutine scan_control_io(netcdf,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning the I/O options in the control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
! amended   - i.j.bush october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( netcdf_param ), Intent( InOut ) :: netcdf
  Type( comms_type ), Intent( InOut ) :: comm

  Logical                :: carry,safe
  Character( Len = 200 ) :: record,record1
  Character( Len = 40  ) :: word,word1
  Real( Kind = wp )      :: tmp

! Some parameters and variables needed for dealing with I/O options

  Integer            :: io_read,io_write,itmp, err_r
  Logical            :: l_io_r,l_io_w,l_tmp
  Integer, Parameter :: MAX_BATCH_SIZE  = 10000000 !~1GB memory per writer
  Integer, Parameter :: MAX_BUFFER_SIZE =   100000 !~1GB memory per node/domain

  Character( Len = 256 ) :: message


! flags

  l_io_r = .false. ! io read  not specified
  l_io_w = .false. ! io write not specified

  safe   = .true.  ! all is safe

! Open the simulation input file

  If (comm%idnode == 0) Inquire(File=Trim(control), Exist=safe)
  Call gcheck(comm,safe,"enforce")
  If (.not.safe) Then
     Go To 10
  Else
     If (comm%idnode == 0) Open(Unit=nread, File=Trim(control), Status='old')
  End If

! Read TITLE record

  Call get_line(safe,nread,record,comm)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record,comm)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! io options

     If      (word(1:2) == 'io' ) Then

        Call get_word( record, word )

        If      (word(1:4) == 'read') Then

           l_io_r  = .true.
           io_read = IO_READ_MPIIO

! get read method

           Call info('',.true.)
           word1=' ' ; word1=word
           Call get_word( record, word )
           If      ( word( 1:5 ) == 'mpiio'  ) Then
              io_read = IO_READ_MPIIO
              Call info('I/O read method: parallel by using MPI-I/O',.true.)
           Else If ( word( 1:6 ) == 'direct' ) Then
              io_read = IO_READ_DIRECT
              Call info('I/O read method: parallel by using direct access',.true.)
           Else If ( word( 1:6 ) == 'netcdf' ) Then
              io_read = IO_READ_NETCDF
              Call info('I/O read method: parallel by using netCDF',.true.)
           Else If ( word( 1:6 ) == 'master' ) Then
              io_read = IO_READ_MASTER
              Call info('I/O read method: serial by using a single master process',.true.)
           Else
              Call strip_blanks(record)
              Write(message,'(4a)') 'io ',word1(1:Len_Trim(word1)+1),word(1:Len_Trim(word)+1),record
              Call info(message,.true.)
              Call error(3)
           End If

           Call io_set_parameters( user_method_read = io_read )

! get number of readers
! 1 <= readers <= mxnode or be wise, default = 8 or 1 when master
! Note set the default number of readers to 1/4 the number of writers using some (limited)
! empirical evidence from HECToR

           If      (io_read == IO_READ_MPIIO  .or. &
                    io_read == IO_READ_DIRECT .or. &
                    io_read == IO_READ_NETCDF) Then

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real(word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
                 itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                 Do While ( Mod( comm%mxnode, itmp ) /= 0 )
                    itmp = itmp - 1
                 End Do
                 Write(message,'(a,i10)') 'I/O readers (assumed) ',itmp
                 Call info(message,.true.)
              Else
                 If (itmp > comm%mxnode) Then
                    tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
                    itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                    Do While ( Mod( comm%mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    Write(message,'(a,i10)') 'I/O readers (enforced) ',itmp
                    Call info(message,.true.)
                 Else
                    Do While ( Mod( comm%mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    Write(message,'(a,i10)') 'I/O readers set to ',itmp
                    Call info(message,.true.)
                 End If
              End If

! the number of readers is now ready to set

              Call io_set_parameters( user_n_io_procs_read = itmp )

! get read batch size
! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
! Note zero or negative values indicate use the default

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real(word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 Call io_get_parameters( user_batch_size_read = itmp )
                 Write(message,'(a,i10)') 'I/O read batch size (assumed) ', itmp
                 Call info(message,.true.)
              Else
                 itmp = Min( itmp, MAX_BATCH_SIZE )
                 Call io_set_parameters( user_batch_size_read = itmp )
                 Write(message,'(a,i10)') 'I/O read batch size set to ', itmp
                 Call info(message,.true.)
              End If
           Else If (io_read == IO_READ_MASTER) Then
              Write(message,'(a,i10)') 'I/O readers (enforced) ', 1
              Call info(message,.true.)
           Else
              tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
              itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
              Write(message,'(a,i10)') 'I/O readers (enforced) ', itmp
              Call info(message,.true.)
! the number of readers is now ready to set
              Call io_set_parameters( user_n_io_procs_read = itmp )
           End If

! get read buffer size
! 100 <= buffer <= MAX_BUFFER_SIZE, default 20000
! Note zero or negative values indicate use the default

           Call get_word( record, word )
           itmp = Nint( Abs( word_2_real(word, 0.0_wp ) ) )
           If (itmp == 0) Then
              Call io_get_parameters( user_buffer_size_read = itmp )
              Write(message,'(a,i10)') 'I/O read buffer size (assumed) ',itmp
              Call info(message,.true.)
           Else
              itmp = Min( Max( itmp,100 ),Min( itmp,MAX_BUFFER_SIZE ))
              Call io_set_parameters( user_buffer_size_read = itmp )
              Write(message,'(a,i10)') 'I/O read buffer size set to ',itmp
              Call info(message,.true.)
           End If

! switch error checking flag for reading

           If (io_read /= IO_READ_MASTER) Then
              Call get_word( record, word )
              l_tmp = ( word( 1:1 ) == 'Y' .or. word( 1:1 ) == 'y' )
              If (.not.l_tmp) Then
                 Call info('I/O parallel read error checking off',.true.)
              Else
                 Call info('I/O parallel read error checking on',.true.)
              End If
              Call io_set_parameters( user_error_check = l_tmp )
           End If

        Else If (word(1:4) == 'writ'  .or. &
                 word(1:5) == 'mpiio' .or. word(1:6) == 'direct' .or. word(1:6) == 'netcdf' .or. word(1:6) == 'master') Then

           l_io_w   = .true.
           io_write = IO_WRITE_SORTED_MPIIO

! for backwards compatibility
! see the second line of the "Else If" above

           If (word(1:4) /= 'writ') Then
              word1=' ' ; word1='write'
           Else
              word1=' ' ; word1=word
              Call get_word( record, word )
           End If

! get write method

           Call info('',.true.)
           If      ( word( 1:5 ) == 'mpiio'  ) Then
              io_write = IO_WRITE_SORTED_MPIIO
              Call info('I/O write method: parallel by using MPI-I/O',.true.)
           Else If ( word( 1:6 ) == 'direct' ) Then
              io_write = IO_WRITE_SORTED_DIRECT
              Call info('I/O write method: parallel by using direct access',.true.)
              Call warning('in parallel this I/O write method has portability issues',.true.)
           Else If ( word( 1:6 ) == 'netcdf' ) Then
              io_write = IO_WRITE_SORTED_NETCDF
              Call get_word( record, word ) ! Check if the user wants the "amber-like/32-bit" format
              If ( word( 1:5 ) == 'amber' .or. word( 1:5 ) == '32bit' ) Then
                 ! Use 32-bit quantities in output for real numbers
                 If (comm%idnode == 0) &
                 Call info('I/O write method: parallel by using netCDF in the amber-like/32-bit format',.true.)
                 Call io_nc_set_real_precision( real32, netcdf, err_r )
              Else
                 ! Use 64-bit quantities in output for real numbers
                 Call info('I/O write method: parallel by using netCDF in 64-bit format',.true.)
                 Call io_nc_set_real_precision( real64, netcdf, err_r )
                 record1=' '
                 record1=word(1:Len_Trim(word)+1) // record ! back up
                 record=record1
              End If
           Else If ( word( 1:6 ) == 'master' ) Then
              io_write = IO_WRITE_SORTED_MASTER
              Call info('I/O write method: serial by using a single master process',.true.)
           Else
              Call strip_blanks(record)
              Write(message,'(4a)') 'io ',word1(1:Len_Trim(word1)+1),word(1:Len_Trim(word)+1),record
              Call info(message,.true.)
              Call error(3)
           End If

! get write type

           Call get_word( record, word )
           If      ( word( 1:6 ) == 'unsort' ) Then
              Call info('I/O write type: data sorting off',.true.)
              Select Case( io_write )
              Case( IO_WRITE_SORTED_MPIIO  )
                 io_write = IO_WRITE_UNSORTED_MPIIO
              Case( IO_WRITE_SORTED_DIRECT )
                 io_write = IO_WRITE_UNSORTED_DIRECT
              Case( IO_WRITE_SORTED_MASTER )
                 io_write = IO_WRITE_UNSORTED_MASTER
              End Select
           Else
              If ( word( 1:6 ) == 'sorted' ) Then
                 Call info('I/O write type: data sorting on',.true.)
              Else
                 record1=' '
                 record1=word(1:Len_Trim(word)+1) // record ! back up
                 record=record1
                 Call info('I/O write type: data sorting on (assumed)',.true.)
              End If
           End If

! the write method and type are now ready to set

           Call io_set_parameters( user_method_write = io_write )

! get number of writers
! 1 <= writers <= mxnode or be wise, default = 8 or 1 when master

           If      (io_write == IO_WRITE_SORTED_MPIIO  .or. &
                    io_write == IO_WRITE_SORTED_DIRECT .or. &
                    io_write == IO_WRITE_SORTED_NETCDF ) Then

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real(word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
                 itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                 Do While ( Mod( comm%mxnode, itmp ) /= 0 )
                    itmp = itmp - 1
                 End Do
                 Write(message,'(a,i10)') 'I/O writers (assumed) ',itmp
                 Call info(message,.true.)
              Else
                 If (itmp > comm%mxnode) Then
                    tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
                    itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
                    Do While ( Mod( comm%mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    Write(message,'(a,i10)') 'I/O writers (enforced) ',itmp
                    Call info(message,.true.)
                 Else
                    Do While ( Mod( comm%mxnode, itmp ) /= 0 )
                       itmp = itmp - 1
                    End Do
                    Write(message,'(a,i10)') 'I/O writers set to ',itmp
                    Call info(message,.true.)
                 End If
              End If

! the number of writers is now ready to set

              Call io_set_parameters( user_n_io_procs_write = itmp )

! get write batch size
! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
! Note zero or negative values indicate use the default

              Call get_word( record, word )
              itmp = Nint( Abs( word_2_real(word, 0.0_wp ) ) )
              If (itmp == 0) Then
                 Call io_get_parameters( user_batch_size_write = itmp )
                 Write(message,'(a,i10)') 'I/O write batch size (assumed) ',itmp
                 Call info(message,.true.)
              Else
                 itmp = Min( itmp, MAX_BATCH_SIZE )
                 Call io_set_parameters( user_batch_size_write = itmp )
                 Write(message,'(a,i10)') 'I/O write batch size set to ',itmp
                 Call info(message,.true.)
              End If
           Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. io_write == IO_WRITE_SORTED_MASTER) Then
              Write(message,'(a,i10)') 'I/O writers (enforced) ',1
              Call info(message,.true.)
           Else
              tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
              itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
              Write(message,'(a,i10)') 'I/O writers (enforced) ',itmp
              Call info(message,.true.)
! the number of writers is now ready to set
              Call io_set_parameters( user_n_io_procs_write = itmp )
           End If

! get write buffer size
! 100 <= buffer <= MAX_BUFFER_SIZE, default 20000
! Note zero or negative values indicate use the default

           Call get_word( record, word )
           itmp = Nint( Abs( word_2_real(word, 0.0_wp ) ) )
           If (itmp == 0) Then
              Call io_get_parameters( user_buffer_size_write = itmp )
              Write(message,'(a,i10)') 'I/O write buffer size (assumed) ',itmp
              Call info(message,.true.)
           Else
              itmp = Min( Max( itmp,100 ),Min( itmp,MAX_BUFFER_SIZE ))
              Call io_set_parameters( user_buffer_size_write = itmp )
              Write(message,'(a,i10)') 'I/O write buffer size set to ',itmp
              Call info(message,.true.)
           End If

! switch error checking flag for writing

           If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
              Call get_word( record, word )
              l_tmp = ( word( 1:1 ) == 'Y' .or. word( 1:1 ) == 'y' )
              If (.not.l_tmp) Then
                 Call info('I/O parallel write error checking off',.true.)
              Else
                 Call info('I/O parallel write error checking on',.true.)
              End If
              Call io_set_parameters( user_error_check = l_tmp )
           End If

        Else If ((word(1:6) == 'output') .or. (word(1:6) == 'config') .or. &
          (word(1:5) == 'field') .or. (word(1:7) == 'statis') .or. (word(1:7) == 'history') &
          .or. (word(1:7) == 'historf') .or. (word(1:6) == 'revive') .or. &
          (word(1:6) == 'revcon') .or. (word(1:6) == 'revold')) Then

          If (word(1:6) == 'output') Then
            Call info('OUTPUT file is '//Trim(output),.true.)
          Else If (word(1:6) == 'config') Then
            Call info('CONFIG file is '//Trim(config),.true.)
          Else If (word(1:5) == 'field') Then
            Call info('FIELD file is '//Trim(field),.true.)
          Else If (word(1:6) == 'statis') Then
            Call info('STATIS file is '//Trim(statis),.true.)
          Else If (word(1:7) == 'history') Then
            Call info('HISTORY file is '//Trim(history),.true.)
          Else If (word(1:7) == 'historf') Then
            Call info('HISTORF file is '//Trim(historf),.true.)
          Else If (word(1:6) == 'revive') Then
            Call info('REVIVE file is '//Trim(revive),.true.)
          Else If (word(1:6) == 'revcon') Then
            Call info('REVCON file is '//Trim(revcon),.true.)
          Else If (word(1:6) == 'revold') Then
            Call info('REVOLD file is '//Trim(revold),.true.)
          End If
! close control file

        Else
           Call strip_blanks(record)
           Write(message,'(2a)') word(1:Len_Trim(word)+1),record
           Call info(message,.true.)
           Call error(3)
        End If

! read finish

     Else If (word(1:6) == 'finish') Then
        carry=.false.
     End If
  End Do

  If (comm%idnode == 0) Close(Unit=nread)

!!! IO DEFAULTS

! io option defaults

  If (.not.l_io_r) Then

! read method

     Call io_set_parameters( user_method_read = IO_READ_MPIIO ) ; io_read = IO_READ_MPIIO
     Call info('',.true.)
     Call info('I/O read method: parallel by using MPI-I/O (assumed)',.true.)

! number of readers

     tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
     itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
     Do While ( Mod( comm%mxnode, itmp ) /= 0 )
        itmp = itmp - 1
     End Do
     Call io_set_parameters( user_n_io_procs_read = itmp )
     Write(message,'(a,i10)') 'I/O readers (assumed) ',itmp
     Call info(message,.true.)

! read batch size

     Call io_get_parameters( user_batch_size_read = itmp )
     Write(message,'(a,i10)') 'I/O read batch size (assumed) ',itmp
     Call info(message,.true.)

! read buffer size

     Call io_get_parameters( user_buffer_size_read = itmp )
     Write(message,'(a,i10)') 'I/O read buffer size (assumed) ',itmp
     Call info(message,.true.)

! error checking flag for reading

     If (io_read /= IO_READ_MASTER) Then
        Call io_set_parameters( user_error_check = .false. )
        Call info('I/O parallel read error checking off (assumed)',.true.)
     End If

  End If

  If (.not.l_io_w) Then

! write method

     Call io_set_parameters( user_method_write = IO_WRITE_SORTED_MPIIO ) ; io_write = IO_WRITE_SORTED_MPIIO
     call info('',.true.)
     Call info('I/O write method: parallel by using MPI-I/O (assumed)',.true.)

! write type

     Call info('I/O write type: data sorting on (assumed)',.true.)

! number of writers

     tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
     itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
     Do While ( Mod( comm%mxnode, itmp ) /= 0 )
        itmp = itmp - 1
     End Do
     Call io_set_parameters( user_n_io_procs_write = itmp )
     Write(message,'(a,i10)') 'I/O writers (assumed) ',itmp
     Call info(message,.true.)

! batch size

     Call io_get_parameters( user_batch_size_write = itmp )
     Write(message,'(a,i10)') 'I/O write batch size (assumed) ',itmp
     Call info(message,.true.)

! write buffer size

     Call io_get_parameters( user_buffer_size_write = itmp )
     Write(message,'(a,i10)') 'I/O write buffer size (assumed) ',itmp
     Call info(message,.true.)

! error checking flag for writing

     If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
        Call io_set_parameters( user_error_check = .false. )
        Call info('I/O parallel write error checking off (assumed)',.true.)
     End If

  End If

  If (io_write == IO_WRITE_SORTED_NETCDF .or. io_read == IO_READ_NETCDF) Call io_nc_compiled()
  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control_io

Subroutine scan_control_output(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning the I/O filenames in the control file
!
! copyright - daresbury laboratory
! author    - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( comms_type ), Intent( InOut ) :: comm
  Logical                :: carry,safe
  Character( Len = 200 ) :: record,rec_case_sensitive
  Character( Len = 40  ) :: word,word1,wordo

  safe   = .true.  ! all is safe

! Open the simulation input file

  If (comm%idnode == 0) Inquire(File=Trim(control), Exist=safe)
  Call gcheck(comm,safe,"enforce")
  If (.not.safe) Then
    Open(Unit=nrite, File=Trim(output), Status='replace')
    Call error(126)
  Else
     If (comm%idnode == 0) Open(Unit=nread, File=Trim(control), Status='old')
  End If

! Read TITLE record

  Call get_line(safe,nread,record,comm)
  If (.not.safe) Then
    Open(Unit=nrite, File=Trim(output), Status='replace')
    Call error(17)
  End If

  carry = .true.
  Do While (carry)
     Call get_line(safe,nread,record,comm)
     If (.not.safe) Then
       Call error(17)
     End If
     rec_case_sensitive = record
     Call lower_case(record)

     Call get_word( record, word )
     Call get_word( rec_case_sensitive, wordo )
     If      (word(1:2) == 'io' ) Then

        Call get_word( record, word1 )
        Call get_word( rec_case_sensitive, wordo )
        If (word1(1:6) == 'output') Then
          Call get_word( rec_case_sensitive, output )

        Else If (word1(1:6) == 'config') Then
          Call get_word( rec_case_sensitive, config )

        Else If (word1(1:5) == 'field') Then
          Call get_word( rec_case_sensitive, field )

        Else If (word1(1:6) == 'statis') Then
          Call get_word( rec_case_sensitive, statis )

        Else If (word1(1:7) == 'history') Then
          Call get_word( rec_case_sensitive, history )

        Else If (word1(1:7) == 'historf') Then
          Call get_word( rec_case_sensitive, historf )

        Else If (word1(1:6) == 'revive') Then
          Call get_word( rec_case_sensitive, revive )

        Else If (word1(1:6) == 'revcon') Then
          Call get_word( rec_case_sensitive, revcon )

        Else If (word1(1:6) == 'revold') Then
          Call get_word( rec_case_sensitive, revold )
        End If
! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     End If

  End Do

  If (comm%idnode == 0) Close(Unit=nread)

End Subroutine scan_control_output

End Module kontrol
