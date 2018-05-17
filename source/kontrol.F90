Module kontrol
  Use kinds, only : wp
  Use comms,      Only : comms_type,gcheck
  Use timer,      Only : timer_type
  Use configuration,     Only : sysname
  Use mpole,     Only : thole
  Use dpd,        Only : keydpd,gamdpd
  Use langevin,   Only : l_lan,l_gst,langevin_allocate_arrays
  Use bonds,      Only : rcbnd
  Use vdw,        Only : ld_vdw,ls_vdw,mxtvdw
  Use metal,      Only : ld_met,ls_met,tabmet
  Use poisson,    Only : eps,mxitcg,mxitjb
  Use msd,        Only : l_msd
  Use defects,   Only : l_dfx
  Use kinetics,  Only : l_vom
  Use plumed,   Only : l_plumed, plumed_input, plumed_log, &
                              plumed_precision, plumed_restart
  Use setup,       Only : nread,control,pi,zero_plus,seed, &
                                            output,field,config,statis, &
                                  history,historf,revive,revcon,revold
  Use parse,       Only : get_line,get_word,lower_case,word_2_real
  
  Use kim,         Only : kimim,rkim
  Use greenkubo,   Only : isvaf,nsvaf,vafsamp
  Use rdfs,         Only : l_errors_jack, l_errors_block
  Use development, Only : l_trm,l_eng, l_rout,l_dis,r_dis,l_tor,&
                          l_his,l_scl,l_rin,l_org,xorg,yorg,zorg,&
                          lvcforg,lvcfscl,cels
  Use ttm
  
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
  
  Use numerics, Only : dcell, invert
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
           rcut,rpad,rvdw,rbin,nstfce,alpha,width,     &
           l_exp,lecx,lfcap,l_top,lzero,lmin,          &
           ltgaus,ltscal,lvar,leql,lpse,               &
           lfce,lpana,lrdf,lprdf,lzdn,lpzdn,           &
           lvafav,lpvaf,ltraj,ldef,lrsd,               &
           nx,ny,nz,imd,tmd,emd,vmx,vmy,vmz,           &
           temp,press,strext,keyres,                   &
           tstep,mndis,mxdis,mxstp,nstrun,nsteql,      &
           keymin,nstmin,min_tol,                      &
           nstzero,nstgaus,nstscal,                    &
           keyens,iso,taut,chi,chi_ep,chi_es,soft,gama,&
           taup,tai,ten,vel_es2,keypse,wthpse,tmppse,  &
           fmax,nstbpo,intsta,keyfce,epsq,             &
           rlx_tol,mxshak,tolnce,mxquat,quattol,       &
           nstbnd,nstang,nstdih,nstinv,nstrdf,nstzdn,  &
           nstmsd,istmsd,nstraj,istraj,keytrj,         &
           nsdef,isdef,rdef,nsrsd,isrsd,rrsd,          &
           ndump,pdplnc,tmr,comm)

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
  Real( Kind = wp ),      Intent( In    ) :: rcut,rpad,rvdw,rbin,width
  Real( Kind = wp ),      Intent( InOut ) :: alpha

  Logical,                Intent(   Out ) :: l_exp,lecx,            &
                                             lfcap,l_top,           &
                                             lmin,lzero,            &
                                             ltgaus,ltscal,         &
                                             lvar,leql,lpse,lfce,   &
                                             lpana,                 &
                                             lrdf,lprdf,lzdn,lpzdn, &
                                             lvafav,lpvaf,          &
                                             ltraj,ldef,lrsd


  Integer,                Intent(   Out ) :: nx,ny,nz,imd,tmd,     &
                                             keyres,nstrun,        &
                                             nsteql,nstzero,       &
                                             keymin,nstmin,        &
                                             nstgaus,nstscal,      &
                                             keyens,iso,           &
                                             keypse,nstbpo,        &
                                             intsta,keyfce,        &
                                             mxshak,mxquat,        &
                                             nstbnd,nstang,        &
                                             nstdih,nstinv,        &
                                             nstrdf,nstzdn,        &
                                             nstmsd,istmsd,        &
                                             nstraj,istraj,keytrj, &
                                             nsdef,isdef,          &
                                             nsrsd,isrsd,          &
                                             ndump

  Real( Kind = wp ),      Intent(   Out ) :: emd,vmx,vmy,vmz,            &
                                             temp,press,strext(1:9),     &
                                             tstep,mndis,mxdis,mxstp,    &
                                             taut,chi,chi_ep,chi_es,soft,&
                                             gama,taup,tai,ten,vel_es2,  &
                                             wthpse,tmppse,min_tol(1:2), &
                                             fmax,epsq,rlx_tol(1:2),     &
                                             tolnce,quattol,             &
                                             rdef,rrsd,pdplnc           
  Type( timer_type ),     Intent( InOut)   :: tmr
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
! ld_vdw = .false. ! (initialised in vdw)
! ls_vdw = .false. ! (initialised in vdw)
  mxtvdw = 0       ! (initialised in vdw)
!
! defaults for direct evaluation of metal interactions
!
! ld_met = .false. ! (initialised in metal_module)

! default impact option: option applied, particle index,
! timestep of impact, energy of impact, (3) direction of impact

  limp = .false.
  imd  = 0
  tmd  = -1
  emd  = 0.0_wp
  vmx  = 0.0_wp
  vmy  = 0.0_wp
  vmz  = 0.0_wp

! temperature & pressure (stress) switches and default values

  ltemp   = .false.
  lpres   = .false.
  lstrext = .false.
  temp    = 0.0_wp
  press   = 0.0_wp
  strext  = 0.0_wp

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

  lpse   = .false.
  keypse = 0
  wthpse = 2.0_wp
  tmppse = 1.0_wp

! default switch for conjugate gradient minimisation during equilibration

  lmin   = .false.
  keymin = -1
  nstmin = 0
  min_tol(1:2) = (/ 0.0_wp , -1.0_wp /) ! tolerance, optional CGM step

! default switch for regaussing temperature and default number of
! steps when to be applied

  ltgaus  = .false.
  nstgaus = 0

! default switch for temperature scaling and default number of
! steps when to be applied

  ltscal  = .false.
  nstscal = 0

! default switch for zero temperature optimisation and default number of
! steps when to be applied

  lzero   = .false.
  nstzero = 0
  l_0     = .false. ! T/=10K

! default integration type (VV), ensemble switch (not defined) and key

  lvv    = .true.
  lens   = .false.
  keyens = 0

! default thermostat and barostat friction time constants

  taut   = 0.0_wp ! thermostat relaxation time
  chi    = 0.0_wp ! Stochastic Dynamics (SD Langevin) thermostat friction
  chi_ep = 0.5_wp ! Inhomogeneous Stochastic Dynamics (SD Langevin) 
                  ! thermostat/electron-phonon friction
  chi_es = 0.0_wp ! Inhomogeneous Stochastic Dynamics (SD Langevin)
                  ! thermostat/electronic stopping friction
  soft   = 0.0_wp ! Softness for Andersen thermostat
  gama   = 0.0_wp ! Stochastic (Langevin) friction on a thermostat
  taup   = 0.0_wp ! barostat relaxation time
  tai    = 0.0_wp ! Stochastic Dynamics (SD Langevin) barostat friction
  iso    = 0      ! no semi-isotropic feature
  ten    = 0.0_wp ! surface tension

! default value for inhomogeneous Langevin thermostat/
! two-temperature model source term cutoff velocity

  vel_es2 = 50.0_wp

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

  mxshak = 250
  tolnce = 1.0e-6_wp

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

  lrdf   = .false.
  nstrdf = 1
  lprdf  = .false.

! default switch for calculation of z-density profile, default number of steps
! when to be collected and default switch for printing it

  lzdn   = .false.
  nstzdn = 1
  lpzdn  = .false.

! default switches for calculation of velocity autocorrelation functions:
! time-averaging and printing

  lvafav = .true.
  lpvaf  = .false.

! default for data printing interval

  nstbpo = 100

! default for statistics file interval

  intsta = 100

! default switch for MSD outputing and defaults for
! (i) step to start at, (ii) every step after to be collected

  nstmsd = 0
  istmsd = 1

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

  ldef   =.false.
  nsdef  = 0
  isdef  = 1
  rdef   = Min(0.75_wp,rcut/3.0_wp)

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
     Else If (word(1:5) == 'l_eng') Then
        l_eng = .true.
        Call info('%%% OUTPUT contains an extra last line with E_tot !!! %%%',.true.)
     Else If (word(1:6) == 'l_rout') Then
        l_rout = .true.
        Call info('%%% REVIVE writing in ASCII opted !!! %%%',.true.)
     Else If (word(1:5) == 'l_rin') Then
        l_rin = .true.
        Call info('%%% REVOLD reading in ASCII opted !!! %%%',.true.)
     Else If (word(1:5) == 'l_org') Then
        l_org = .true.
        l_trm  = .true.

        Call info('%%% translate CONFIG along a vector into CFGORG after reading input & terminate !!! %%%',.true.)
        Call info('%%% vector and config level read as follows: %%%',.true.)

        Call get_word(record,word)
        xorg = word_2_real(word)
        Call get_word(record,word)
        yorg = word_2_real(word)
        Call get_word(record,word)
        zorg = word_2_real(word)

        Call get_word(record,word)
        lvcforg = Min( Int(Abs(word_2_real(word,0.0_wp))) , levcfg)

        Write(messages(1),'(a)') '%%%'
        Write(messages(2),'(a,3f10.3,a)') '%%% vector(x,y,x) ', xorg, yorg, zorg, ' %%%'
        Write(messages(3),'(a,i0,a)') '%%% CFGORG level ', lvcforg, ' %%%'
        Call info(messages,3,.true.)

     Else If (word(1:5) == 'l_scl') Then
        Call info('%%% rescale CONFIG to CFGSCL, after reading input & terminate !!! %%%',.true.)
        Call info('%%% config level and new cell vectors to rescale to (read in a CONFIG-like manner): %%%',.true.)

        Call get_word(record,word)
        lvcfscl = Min( Int(Abs(word_2_real(word,0.0_wp))) , levcfg)

        itmp=0
        Do i=1,3
           Call get_line(safe,nread,record,comm)
           Do j=1,3
              Call get_word(record,word)
              itmp=itmp+1
              cels(itmp)=word_2_real(word)
           End Do
        End Do

        Call invert(cels,rcell,tmp)

        Write(messages(1),'(a)') '%%%'
        Write(messages(2),'(1x,a,i0,a)')        '%%% CFGSCL level ', lvcfscl, ' %%%'
        Write(messages(3),'(1x,a,3f20.10,a)')   '%%% ', cels(1:3), ' %%%'
        Write(messages(4),'(1x,a,3f20.10,a)')   '%%% ', cels(4:6), ' %%%'
        Write(messages(5),'(1x,a,3f20.10,a)')   '%%% ', cels(7:9), ' %%%'
        Write(messages(6),'(1x,a)')             '%%%'
        Write(messages(7),'(1x,a,1p,g22.12,a)') '%%% CFGSCL volume ', tmp, '%%%'
        Call info(messages,7,.true.)

        If (tmp > zero_plus) Then
           l_scl = .true.
           l_trm  = .true.
        Else
           Call info('%%% OPTION ABORTED DUE TO ZERO VOLUME !!! %%%',.true.)
           l_trm  = .true.
        End If
     Else If (word(1:5) == 'l_his') Then
        l_his = .true.
        l_trm = .true.
        Call info('%%% generate HISTORY after reading input & terminate !!! %%%',.true.)
     Else If (word(1:5) == 'l_tim') Then
!        l_tim = .true.  ! done in scan_development
        Call info('%%% generate detailed timing !!! %%%',.true.)
     Else If (word(1:5) == 'l_tor') Then
        l_tor = .true.
        Call info('%%% Turn off production of REVCON & REVIVE !!! %%%',.true.)
     Else If (word(1:5) == 'l_trm') Then
        l_trm = .true.
        Call info('%%% Terminate gracefully before initialisation !!! %%%',.true.)
     Else If (word(1:5) == 'l_dis') Then
        l_dis = .true.
        r_dis = Min( r_dis , word_2_real(word,0.1_wp) )
        Call info('%%% Turn on the check on minimum separation distance between VNL pairs at re/start !!! %%%',.true.)
        Write(message,'(a,1p,e12.4)') '%%% separation criterion (Angstroms) %%%', r_dis
        Call info(message,.true.)

! read VDW options

     Else If (word(1:3) == 'vdw') Then
        Call get_word(record,word1)

        If      (word1(1:6) == 'direct') Then

! direct evaluation option

           ld_vdw = .true.
           Call info('vdw direct option on',.true.)

        Else If (word1(1:6) == 'mixing') Then

! mixing type keywords

           Call info('vdw cross terms mixing opted (for undefined mixed potentials)',.true.)
           Call info('mixing is limited to potentials of the same type only',.true.)
           Call info('mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)',.true.)

           Call get_word(record,word2)

           If      (word2(1:4) == 'lore') Then
              mxtvdw = 1
              Call info('type of mixing selected - Lorentz–Berthelot :: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2',.true.)
           Else If (word2(1:4) == 'fend') Then
              mxtvdw = 2
              Call info('type of mixing selected - Fender-Halsey :: e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2',.true.)
           Else If (word2(1:4) == 'hoge') Then
              mxtvdw = 3
              Call info('type of mixing selected - Hogervorst (good hope) :: ' &
                //'e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)',.true.)
           Else If (word2(1:4) == 'halg') Then
              mxtvdw = 4
              Call info('type of mixing selected - Halgren HHG :: ' &
                //'e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)',.true.)
           Else If (word2(1:4) == 'wald') Then
              mxtvdw = 5
              Call info('type of mixing selected - Waldman–Hagler :: ' &
                //'e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ;s_ij=[(s_i^6+s_j^6)/2]^(1/6)',.true.)
           Else If (word2(1:4) == 'tang') Then
              mxtvdw = 6
              Call info('type of mixing selected - Tang-Toennies :: ' &
                //' e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13',.true.)
              Call info(Repeat(' ',43)//'s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)',.true.)
           Else If (word2(1:4) == 'func') Then
              mxtvdw = 7
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
           ls_vdw = .true.
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
           If (tabmet > 0) Then
              Call warning(480,0.0_wp,0.0_wp,0.0_wp)
           Else
              ld_met = .true.
           End If
        Else If (word(1:7) == 'sqrtrho') Then
! read metal sqrtrho interpolation option for EAM embeding function in TABEAM
           Call info('metal sqrtrho option on',.true.)
           If (tabmet > 0) Then
              ls_met = .true.
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
        imd = Max(1,Nint(Abs(word_2_real(word))))
        Call get_word(record,word)
        tmd = Nint(Abs(word_2_real(word)))

        Call get_word(record,word)
        emd = Abs(word_2_real(word))
        Call get_word(record,word)
        vmx = word_2_real(word)
        Call get_word(record,word)
        vmy = word_2_real(word)
        Call get_word(record,word)
        vmz = word_2_real(word)

        If (Sqrt(vmx**2+vmy**2+vmz**2) <= zero_plus) Then
           vmx = 1.0_wp
           vmy = 1.0_wp
           vmz = 1.0_wp
        End If

        Write(messages(1),'(a)') ''
        Write(messages(2),'(a,i10)') 'particle (index)',imd
        Write(messages(3),'(a,i10)') 'timestep (steps)',tmd
        Write(messages(4),'(a,1p,e12.4)') 'energy   (keV)  ',emd
        Write(messages(5),'(a,1p,3e12.4)') 'v-r(x,y,z)      ',vmx,vmy,vmz
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
        temp = Abs(word_2_real(word))
        Write(message,'(a,1p,e12.4)') 'simulation temperature (K)  ',temp
        Call info(message,.true.)

! read zero temperature optimisation

     Else If (word(1:4) == 'zero') Then

        lzero = .true.

! Check defaults

        Call get_word(record,word)
        l_0 = (word(1:4) == 'fire')
        nstzero = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        If (word(1:5) == 'every') Call get_word(record,word)
        nstzero = Max(nstzero,Abs(Nint(word_2_real(word,0.0_wp))))

        Call info('zero K optimisation on (during equilibration)',.true.)
        Write(message,'(a,i10)') 'temperature regaussing interval',nstzero

        If (l_0) Then
           If (comm%idnode == 0) &
           Call info('fire option on - actual temperature will reset to 10 Kelvin if no target tempreature is specified',.true.)
        Else
           ltemp  = .true.
           temp = 10.0_wp
           Call info('fire option off - actual temperature reset to 10 Kelvin',.true.)
        End If

! read pressure

     Else If (word(1:4) == 'pres') Then

        Call get_word(record,word)

! read stress (6 components only: xx,yy,zz,xy,xz,yz - it's forced symmetric)

        If (word(1:6) == 'tensor') Then
           lstrext=.true.

           Call get_word(record,word)
           strext(1) = word_2_real(word)
           Call get_word(record,word)
           strext(5) = word_2_real(word)
           Call get_word(record,word)
           strext(9) = word_2_real(word)
           Call get_word(record,word)
           strext(2) = word_2_real(word)
           strext(4) = strext(2)
           Call get_word(record,word)
           strext(3) = word_2_real(word)
           strext(7) = strext(3)
           Call get_word(record,word)
           strext(6) = word_2_real(word)
           strext(8) = strext(6)

           Call info('simulation pressure tensor (katms):',.true.)
           Write(messages(1),'(3f20.10)') strext(1:3)
           Write(messages(2),'(3f20.10)') strext(4:6)
           Write(messages(3),'(3f20.10)') strext(7:9)
           Call info(messages,3,.true.)

! convert from katms to internal units of pressure

           strext = strext/prsunt
        Else
           lpres=.true.

           press = word_2_real(word)

           Write(message,'(a,1p,e12.4)') 'simulation pressure (katms) ',press

! convert from katms to internal units of pressure

           press = press/prsunt
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
           keypse = 1
           Call get_word(record,word)
        Else If (word(1:5) == 'gauss') Then
           keypse = 2
           Call get_word(record,word)
        Else If (word(1:6) == 'direct') Then
           keypse = 3
           Call get_word(record,word)
        End If

! wthpse = 2 Angs by default

        tmp = Abs(word_2_real(word))
        If (width/4.0_wp > wthpse) Then
           lpse = .true.
           If (comm%idnode == 0) Then
              Call info('pseudo thermostat attached to MD cell boundary',.true.)
              If      (keypse == 0) Then
                Call info('thermostat control: Langevin + direct temperature scaling',.true.)
              Else If (keypse == 1) Then
                Call info('thermostat control: Langevin temperature scaling',.true.)
              Else If (keypse == 2) Then
                Call info('thermostat control: gaussian temperature scaling',.true.)
              Else If (keypse == 3) Then
                Call info('thermostat control: direct temperature scaling',.true.)
              End If
              Write(message,'(a,1p,e12.4)') 'thermostat thickness (Angs) ',tmp
           End If

           If (width/4.0_wp > tmp .and. tmp >= wthpse) Then
              wthpse = tmp
           Else
              Call info('thermostat thickness insufficient - reset to 2 Angs',.true.)
           End If
        Else
           Call warning(280,wthpse,width,0.0_wp)
           Call error(530)
        End If

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (tmp <= zero_plus) Then
           tmppse = temp
        Else
           tmppse = Max(tmppse,tmp)
        End If
        Write(message,'(a,1p,e12.4)') 'thermostat temperature (K) ',tmppse
        Call info(message,.true.)

! read minimiser option

     Else If (word(1:5) == 'minim' .or. word(1:5) == 'optim') Then

        lmin=.true.
        word2=' ' ; word2=word
        Call get_word(record,word)

        If      (word(1:4) == 'forc') Then
           keymin=0
           word1='force   '
        Else If (word(1:4) == 'ener') Then
           keymin=1
           word1='energy  '
        Else If (word(1:4) == 'dist') Then
           keymin=2
           word1='distance'
        Else
           Call strip_blanks(record)
           Write(message,'(4a)') word2(1:Len_Trim(word2)+1),' ',word(1:Len_Trim(word)+1),record
           Call info(message,.true.)
           Call error(590)
        End If

        If (word2(1:5) == 'minim') Then
           Call get_word(record,word)
           nstmin = Abs(Nint(word_2_real(word,0.0_wp)))
        End If

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))

        itmp=0
        If      (keymin == 0) Then
           If (tmp < 1.0_wp .or. tmp > 1000.0_wp) Then
              min_tol(1)=50.0_wp
              itmp=1
           Else
              min_tol(1)=tmp
           End If
        Else If (keymin == 1) Then
           If (tmp < zero_plus .or. tmp > 0.01_wp) Then
              min_tol(1)=0.005_wp
              itmp=1
           Else
              min_tol=tmp
           End If
        Else If (keymin == 2) Then
           If (tmp < 1.0e-6_wp .or. tmp > 0.1_wp) Then
              min_tol(1)=0.005_wp
              itmp=1
           Else
              min_tol(1)=tmp
           End If
        End If

        If (itmp == 1) Call warning(360,tmp,min_tol(1),0.0_wp)

        Call get_word(record,word3)
        min_tol(2) = word_2_real(word3,-1.0_wp)

        If (word2(1:5) == 'minim') Then
           Write(messages(1),'(a)') 'minimisation option on (during equilibration)'
           Write(messages(2),'(a,a8)') 'minimisation criterion        ',word1(1:8)
           Write(messages(3),'(a,i10)') 'minimisation frequency (steps)',nstmin
           Write(messages(4),'(a,1p,e12.4)') 'minimisation tolerance        ',min_tol(1)
           Call info(messages,4,.true.)
           If (min_tol(2) > zero_plus) Then
             Write(message,'(a,1p,e12.4)') 'minimisation CGM step         ',min_tol(2)
             Call info(message,.true.)
           End If
        Else
           Write(messages(1),'(a)') 'optimisation at start'
           Write(messages(2),'(a,a8)') 'optimisation criterion        ',word(1:8)
           Write(messages(4),'(a,1p,e12.4)') 'optimisation tolerance        ',min_tol(1)
           Call info(messages,3,.true.)
           If (min_tol(2) > zero_plus) Then
             Write(message,'(a,1p,e12.4)')'optimisation CGM step         ',min_tol(2)
           End If
        End If

! read regauss option

     Else If (word(1:6) == 'regaus') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        nstgaus = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        ltgaus =.true.
        Write(message,'(a,i10)') 'temperature regaussing interval ', nstgaus
        Call info(message,.true.)

! read temperature scaling option

     Else If (word(1:5) == 'scale') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        nstscal = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        ltscal =.true.
        Call info('temperature scaling on (during equilibration)',.true.)
        Write(message,'(a,i10)') 'temperature scaling interval ',nstscal
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
              thole = Abs(word_2_real(word,0.0_wp))
           End If

           Write(message,'(a,f5.2)') &
             'CHARMM polarisation scheme selected with optional atomic thole dumping of ', &
             thole
           Call info(message,.true.)
           If (mximpl == 0 ) Then
             Call warning('scheme deselected due to switched off electrostatics',.true.)
           End If
           If (mxshl == 0) Then
             Call warning('scheme disabled due to lack of core-shell defined interatcions',.true.)
           End If

           If (mximpl == 0 .or. mxshl == 0) Then
!              keyind=0 ! done in scan_control
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

! keydpd detected in scan_control

        If (keydpd > 0 .and. (lvv .neqv. l_vv)) Then
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

           keyens = 0

           Call info('Ensemble : NVE (Microcanonical)',.true.)

           If (lens) Call error(414)
           lens=.true.

        Else If (word(1:3) == 'nvt') Then

           Call get_word(record,word)

           If      (word(1:5) == 'evans') Then

              keyens = 1

              Call info('Ensemble : NVT Evans (Isokinetic)',.true.)
              Call info('Gaussian temperature constraints in use',.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:4) == 'lang') Then

              keyens = 10
              If (.not.l_vv) l_lan = .true.

              Call get_word(record,word)
              chi = Abs(word_2_real(word))

              Call info('Ensemble : NVT Langevin (Stochastic Dynamics)',.true.)
              Write(message,'(a,1p,e12.4)') 'thermostat friction ', chi
              Call info(message,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:5) == 'ander') Then

              keyens = 11

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              soft = Abs(word_2_real(word))
              If (soft > 1.0_wp) soft=1.0_wp/soft

              Write(messages(1),'(a)') 'Ensemble : NVT Andersen'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'softness (dimensionless)',soft
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              keyens = 12

              Call get_word(record,word)
              taut = Abs(word_2_real(word))

              Call info('Ensemble : NVT Berendsen',.true.)
              Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Call info(message,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              keyens = 13

              Call get_word(record,word)
              taut = Abs(word_2_real(word))

              Call info('Ensemble : NVT Nose-Hoover',.true.)
              Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Call info(message,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'gst') Then

              keyens = 14
              l_gst = .true.

              Call get_word(record,word)
              taut = Abs(word_2_real(word))

              Call get_word(record,word)
              gama = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NVT gentle stochastic thermostat'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'friction on thermostat  (ps^-1) ',gama
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:4) == 'ttm' .or. word(1:6) == 'inhomo') Then

              keyens = 15
              If (.not.l_vv) l_lan = .true.

              Call get_word(record,word)
              chi_ep  = Abs(word_2_real(word))
              Call get_word(record,word)
              chi_es  = Abs(word_2_real(word))
              Call get_word(record,word)
              vel_es2 = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)'
              Write(messages(2),'(a,1p,e12.4)') 'e-phonon friction (ps^-1) ',chi_ep
              Write(messages(3),'(a,1p,e12.4)') 'e-stopping friction (ps^-1) ',chi_es
              Write(messages(4),'(a,1p,e12.4)') 'e-stopping velocity (A ps^-1) ',vel_es2
              Call info(messages,4,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'dpd') Then

              Call info('Ensemble : NVT dpd (Dissipative Particle Dynamics)',.true.)

! keydpd determined in scan_control

              If      (keydpd == 1) Then
                 keyens = 0 ! equivalence to doing NVE with some extra fiddling before VV(0)
                 Call info("Ensemble type : Shardlow's first order splitting (S1)",.true.)
              Else If (keydpd == 2) Then
                 keyens = 0 ! equivalence to doing NVE with some extra fiddling before VV(0) and after VV(1)
                 Call info("Ensemble type : Shardlow's second order splitting (S2)",.true.)
              Else
                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call error(436)
              End If

              Call get_word(record,word)
              gamdpd(0) = Abs(word_2_real(word,0.0_wp))

              If (gamdpd(0) > zero_plus) Then
                 Write(message,'(a,1p,e12.4)') 'drag coefficient (Dalton/ps) ', gamdpd(0)
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

           Call get_word(record,word)

           If (word(1:4) == 'lang') Then

              keyens = 20
              l_lan = .true.

              Call get_word(record,word)
              chi = Abs(word_2_real(word))
              Call get_word(record,word)
              tai = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Langevin (Stochastic Dynamics)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat friction (ps^-1)',chi
              Write(messages(3),'(a,1p,e12.4)') 'barostat friction (ps^-1)',tai
              Call info(messages,3,.true.)

!                 taut=1/(2.0_wp*pi*chi)
!                 taup=1/(2.0_wp*pi*tai)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              keyens = 21

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Berendsen'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',taup
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              keyens = 22

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Nose-Hoover (Melchionna)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',taup
              Call info(messages,3,.true.)

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'mtk') Then

              keyens = 23

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT isotropic Martyna-Tuckerman-Klein'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',taup
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

           Call get_word(record,word)

           If (word(1:4) == 'lang') Then

              keyens = 30
              l_lan = .true.

              Call get_word(record,word)
              chi = Abs(word_2_real(word))
              Call get_word(record,word)
              tai = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Langevin (Stochastic Dynamics)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat friction (ps^-1)',chi
              Write(messages(3),'(a,1p,e12.4)') 'barostat friction (ps^-1)',tai
              Call info(messages,3,.true.)

!                 taut=chi
!                 taup=2.0_wp*pi/tai

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', ten
                 Call info(messages,3,.true.)
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
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
                    iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    iso=3
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
              If (iso >= 1 .and. iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis',.true.)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              keyens = 31

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Berendsen'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',taup
              Call info(messages,3,.true.)

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', ten
                 Call info(messages,3,.true.)
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
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
                    iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    iso=3
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
              If (iso >= 1 .and. iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis')
              End If

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              keyens = 32

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Nose-Hoover (Melchionna)'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',taup
              Call info(messages,3,.true.)

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', ten
                 Call info(messages,3,.true.)
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
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
                    iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    iso=3
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
              If (iso >= 1 .and. iso <= 2) Then
                 Call warning('semi-isotropic ensembles are only correct for infinite' &
                   //'interfaces placed perpendicularly to the z axis',.true.)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'mtk') Then

              keyens = 33

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Martyna-Tuckerman-Klein'
              Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',taut
              Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',taup
              Call info(messages,3,.true.)

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
                 Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
                 Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
                 Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', ten
                 Call info(messages,3,.true.)
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
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
                    iso=2
                    Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
                 Else If (word(1:4) == 'semi') Then
                    iso=3
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
              If (iso >= 1 .and. iso <= 2) Then
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

        If (l_lan) Call langevin_allocate_arrays()

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
           tol = Sqrt(Abs(Log(eps0*rcut)))
           alpha = Sqrt(Abs(Log(eps0*rcut*tol)))/rcut
           Write(message,'(a,1p,e12.4)') 'damping parameter (A^-1) derived ', alpha
           Call info(message,.true.)
        End If
        If (alpha > zero_plus) Then
           Call info('Fennell damping applied',.true.)
           If (rcut < 12.0_wp) Call warning(7,rcut,12.0_wp,0.0_wp)
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
           tol = Sqrt(Abs(Log(eps0*rcut)))
           alpha = Sqrt(Abs(Log(eps0*rcut*tol)))/rcut
           Write(message,'(a,1p,e12.4)') 'damping parameter (A^-1) derived ', alpha
           Call info(message,.true.)
        End If
        If (alpha > zero_plus) Then
           Call info('Fennell damping applied',.true.)
           If (rcut < 12.0_wp) Call warning(7,rcut,12.0_wp,0.0_wp)
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

        If ( Abs(prmps(1)-1.0_wp/alpha) > 1.0e-6_wp .or. Abs(prmps(2)-eps) > 1.0e-6_wp .or. &
             Nint(prmps(3)) == 0 .or. Nint(prmps(4)) == 0 ) Then
           Call warning('parameters reset to safe defaults occurred',.true.)
           Write(messages(1),'(a,1p,e12.4)') 'gridspacing parameter (A) ',1.0_wp/alpha
           Write(messages(2),'(a,1p,e12.4)') 'convergance epsilon ',eps
           Write(messages(3),'(a,1p,i5)') 'max # of Psolver iterations ',mxitcg
           Write(messages(4),'(a,1p,i5)') 'max # of Jacobi  iterations ',mxitjb
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

!     Else If (word(1:6) == 'induce') Then
!
!        Write(messages(1),'(a)') 'Employing induced dipoles'
!        Write(messages(2),'(a)') 'Induced dipole conjugate gradient :'
!        Write(messages(3),'(a,i0)') 'max number of steps = ',politer
!        Write(messages(4),'(a,e7.3') 'convergence criterion = ',convcrit
!        Call info(messages,4,.true.)
!
!     Else If (word(1:4) == 'gear') Then
!
!       Write(message,'(a,i0,a)') "Using gear predictor with ",numcof," points"
!       Call info(message,.true.)
!
!     Else If (word(1:4) == 'aspc') Then
!
!       Write(message,'(a,i0,a)') Using always stable predictor corrector with ",numcof," points"
!       Call info(message,.true.)
!
!     Else If (word(1:5) == 'lstsq') Then
!
!       Write(message,'(a,i0,a)' "Using least squares predictor with ",numcof," points"
!       Call info(message,.true.)

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

           lvafav = .false.

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
        mxshak = Abs(Nint(word_2_real(word)))

! read tolerance for constraint algorithms

     Else If (word(1:5) == 'shake') Then

        Call get_word(record,word)
        If (word(1:9) == 'tolerance') Call get_word(record,word)
        tolnce = Abs(word_2_real(word))

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
          Write(messages(2),'(a,1p,e12.4)') 't.c. at system temp. (W m^-1 K^-1) ',Ka0
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

        Else If (word1(1:6) == 'sgauss' .or. word1(1:5) == 'sigma') Then

        ! gaussian spatial distribution for initial energy deposition into
        ! electronic system

          sdepoType = 1
          Call get_word(record,word)
          sig = word_2_real(word)
          Call get_word(record,word)
          sigmax = word_2_real(word)
          Write(messages(1),'(a)') 'initial gaussian spatial energy deposition in electronic system'
          Write(messages(2),'(a,1p,e12.4)') 'sigma of distribution (nm) ',sig
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
          Write(messages(2),'(a,1p,e12.4)') 'sigma of distribution (ps) ',tdepo
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

        ! variable electron-phonon coupling constant (chi_ep) based on
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
        If (Abs(rbin-tmp) > 1.0e-6_wp) Call warning(340,tmp,rcut/4.0_wp,rbin)

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

! read rdf calculation option

     Else If (word(1:3) == 'rdf') Then

        lrdf = .true.

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        nstrdf = Abs(Nint(word_2_real(word,1.0_wp)))

! read z-density profile option

     Else If (word(1:4) == 'zden') Then

        lzdn = .true.

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        nstzdn = Abs(Nint(word_2_real(word,1.0_wp)))

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
        Else If (word(1:3) == 'rdf' ) Then
           lprdf = .true.
        Else If (word(1:4) == 'zden') Then
           lpzdn = .true.
        Else If (word(1:3) == 'vaf') Then
           lpvaf = .true.
        Else
           If (word(1:5) == 'every') Call get_word(record,word)
           nstbpo = Abs(Nint(word_2_real(word,1.0_wp)))
           Write(message,'(a,i10)') 'data printing interval (steps) ',nstbpo
           Call info(message,.true.)
        End If

! read stack option (reading done in set_bounds -> scan_control)

     Else If (word(1:5) == 'stack') Then

        Write(message,'(a,i10)') 'data stacking interval (steps) ',mxstak
        Call info(message,.true.)

! read statistics printing option

     Else If (word(1:4) == 'stat') Then

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        intsta = Nint(word_2_real(word))
        Write(message,'(a,i10)') 'statistics file interval ',intsta
        Call info(message,.true.)

! read MSDTMP printing option

     Else If (word(1:6) == 'msdtmp') Then

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        nstmsd = Max(nstmsd,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        istmsd = Max(istmsd,itmp)

        Write(messages(1),'(a)') 'MSDTMP file option on'
        Write(messages(2),'(2x,a,i10)') 'MSDTMP file start ',nstmsd
        Write(messages(3),'(2x,a,i10)') 'MSDTMP file interval ',istmsd
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

        ldef = .true.

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        nsdef = Max(nsdef,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        isdef = Max(isdef,itmp)

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (tmp >= Min(0.3_wp,rcut/3.0_wp) .and. tmp <= Min(3.5_wp,rcut/2.0_wp)) Then
           rdef = tmp ! 3.43 Angs is the Cs VDW radius - largest possible
        Else
           Call warning(310,tmp,rdef,0.0_wp)
        End If

        Write(messages(1),'(a)') 'defects file option on'
        Write(messages(2),'(2x,a,i10)') 'defects file start ',nsdef
        Write(messages(3),'(2x,a,i10)') 'defects file interval ',isdef
        Write(messages(4),'(2x,a,1p,e12.4)') 'defects distance condition (Angs) ',rdef
        Call info(messages,4,.true.)

! REFERENCE1 forcing

        Call get_word(record,word)
        If (word(1:5) == 'extra') Then
           l_dfx=.true.
           Call info('defects1 file option on',.true.)
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
           l_plumed=.true.

           Call get_word(record,word)
           If (word(1:3) == 'off') Then
              l_plumed=.false.
              lplumed=.false.
           End If

           If (word(1:5) == 'input') Then
              Call get_word(record,word)
              plumed_input=Trim(word)
           End If

           If (word(1:3) == 'log') Then
              Call get_word(record,word)
              plumed_log=Trim(word)
           End If

           If (word(1:9) == 'precision') Then
              Call get_word(record,word)
              plumed_precision=Abs(Nint(word_2_real(word,1.0_wp)))
           End If

           If (word(1:7) == 'restart') Then
              plumed_restart=1

              Call get_word(record,word)

              If ((word(1:3) == 'yes') .or. (word(1:1) == 'y')) Then
                 plumed_restart=1
              End If

              If ((word(1:2) == 'no') .or. (word(1:1) == 'n')) Then
                 plumed_restart=0
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

  If (nstmin  == 0) nstmin  = nsteql+1
  If (nstzero == 0) nstzero = nsteql+1
  If (nstgaus == 0) nstgaus = nsteql+1
  If (nstscal == 0) nstscal = nsteql+1

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
      Write(messages(2),'(a,1p,e12.4)') 'e-phonon friction (ps^-1) ',chi_ep
      Write(messages(3),'(a,1p,e12.4)') 'e-stopping friction (ps^-1) ',chi_es
      Write(messages(4),'(a,1p,e12.4)') 'e-stopping velocity (A ps^-1)',vel_es2
      Call info(messages,4,.true.)
    Else
      Call info('Ensemble : NVE (Microcanonical)',.true.)
    End If
    If (l_ttm) keyens = 15
    lens=.true.
  End If

! report replacement of specified ensemble with inhomogeneous
! Langevin if two-temperature model is in use, replacing
! default electron-phonon friction value with chi from
! standard Langevin thermostat (if supplied), and use of
! thermal velocities only for thermostat

  If (l_ttm .and. keyens/=15) Then
    Call warning(130,0.0_wp,0.0_wp,0.0_wp)
    If (keyens==10 .or. keyens==20 .or. keyens==30 .and. chi>zero_plus) Then
      chi_ep = chi
    End If

    Write(messages(1),'(a)') 'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)'
    Write(messages(2),'(a,1p,e12.4)') 'e-phonon friction (ps^-1) ',chi_ep
    Write(messages(3),'(a,1p,e12.4)') 'e-stopping friction (ps^-1) ',chi_es
    Write(messages(4),'(a,1p,e12.4)') 'e-stopping velocity (A ps^-1)',vel_es2
    Call info(messages,4,.true.)

    If (ttmthvel) Then
      Call info('applying to thermal velocities in all directions',.true.)
    Else If (ttmthvelz) Then
      Call info('applying to total velocities in x and y directions, thermal velocities in z direction',.true.)
    Else
      Call info('applying to total velocities in all directions',.true.)
    End If
    keyens = 15

  End If

! report iteration length and tolerance condition for constraints and PMF algorithms

  If ((mxcons > 0 .or. mxpmf > 0) .and. comm%idnode == 0) Then
     Write(messages(1),'(a,i10)') 'iterations for shake/rattle ',mxshak
     Write(messages(2),'(a,1p,e12.4)') 'tolerance for shake/rattle (Angs) ',tolnce
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

! report if rcut is reset (measures taken in scan_config -
! rcut is the maximum cutoff needed in the system)

  If (Abs(rcut-rcut1) > 1.0e-6_wp) Then
    Write(message,'(a,1p,e12.4)') 'real space cutoff reset to (Angs) ',rcut
    Call info(message,.true.)
  End If

! report if rpad is reset (measures taken in scan_config & set_bounds -
! rpad is the cutoff padding needed the conditional VNL update)

  If (Abs(rpad-rpad1) > 1.0e-6_wp) Then
    Write(message,'(a,1p,e12.4)') 'cutoff padding reset to (Angs) ',rpad
    Call info(message,.true.)
  End If

! report vdw

  If (l_n_v) Then
    Call info('vdw potential terms switched off',.true.)
  End If

! report if rvdw is reset (measures taken in scan_config)

  If ((.not.l_n_v) .and. Abs(rvdw-rvdw1) > 1.0e-6_wp) Then
    Write(message,'(a,1p,e12.4)') 'vdw cutoff reset to (Angs) ',rvdw
    Call info(message,.true.)
  End If

! report timestep

  If (lvar) Then

     If (keydpd > 0) Then
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
        If (mxgbnd1 > 0 .and. mxgang1 > 0 .and. &
            mxgdih1 > 0 .and. mxginv1 > 0) Then
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

        If (mxgbnd1 > 0) Then
           If (nstbnd == 0 .or. (nstbnd > nstana .and. nstana > 0)) Then
              nstbnd = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdbnd /= mxgbnd1)
           k=Merge(1, 0, Abs(rcbnd-rcb_d) > 1.0e-3_wp)
           Write(message,'(2(a,i10),a,f7.2,a)') &
             'bonds      - collection every ',nstbnd,' step(s); ngrid = ', &
             mxgbnd1,' points; cutoff = ',rcbnd, ' Angs'
           Call info(message,.true.)
           If (i+j+k > 1) Then
             Write(message,'(3(a,i10))') &
               'bonds      - reset values at  ', i,'                  ', j, &
               '                 ', k
             Call info(message,.true.)
           End If
        End If

        If (mxgang1 > 0) Then
           If (nstang == 0 .or. (nstang > nstana .and. nstana > 0)) Then
              nstang = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdang /= mxgang1)
           Write(message,'(2(a,i10),a)') &
              'angles     - collection every ',nstang,' step(s); ngrid = ',mxgang1,' points'
           Call info(message,.true.)
           If (i+j > 1) Then
             Write(message,'(2(a,i10))') &
               'angles     - reset values at  ',     i,'                  ',     j
             Call info(message,.true.)
           End If
        End If

        If (mxgdih1 > 0) Then
           If (nstdih == 0 .or. (nstdih > nstana .and. nstana > 0)) Then
              nstdih = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grddih /= mxgdih1)
           Write(message,'(2(a,i10),a)') &
             'dihedrals  - collection every ',nstdih,' step(s); ngrid = ',mxgdih1,' points'
           Call info(message,.true.)
           If (i+j > 1) Then
             Write(message,'(2(a,i10))') &
               'dihedrals  - reset values at  ',     i,'                  ',     j
             Call info(message,.true.)
           End If
        End If

        If (mxginv1 > 0) Then
           If (nstinv == 0 .or. (nstinv > nstana .and. nstana > 0)) Then
              nstinv = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdinv /= mxginv1)
           Write(message,'(2(a,i10),a)') &
             'inversions - collection every ',nstinv,' step(s); ngrid = ',mxginv1,' points'
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

  If (lrdf .or. lprdf) Then
     If (lrdf) Then
        Write(messages(1),'(a)') 'rdf collection requested:'
        Write(messages(2),'(2x,a,i10)') 'rdf collection interval ',nstrdf
        Write(messages(3),'(2x,a,1p,e12.4)') 'rdf binsize (Angstroms) ',rbin
        Call info(messages,3,.true.)
     Else
        Call info('no rdf collection requested',.true.)
     End If

     If (lprdf) Then
        Call info('rdf printing requested',.true.)
     Else
        If (lpana) Then
           Call info('rdf printing triggered due to a PDA printing request',.true.)
           lprdf=lpana
        Else
           Call info('no rdf printing requested',.true.)
        End If
     End If

     If (mxrdf == 0) Then
        Call info('no rdf pairs specified in FIELD',.true.)
     Else
        Call info('rdf pairs specified in FIELD',.true.)
     End If

     If ((.not.lrdf) .or. mxrdf == 0) Then
        Call info('rdf routines not to be activated',.true.)
        lrdf=.false.
        lprdf=.false.
     End If
  End If

! report zden

  If (lzdn .or. lpzdn) Then
     If (lzdn) Then
        Write(messages(1),'(a)') 'z-density profiles requested:'
        Write(messages(2),'(2x,a,i10)') 'z-density collection interval ',nstzdn
        Write(messages(3),'(2x,a,1p,e12.4)') 'z-density binsize (Angstroms) ',rbin
        Call info(messages,3,.true.)
     Else
        Call info('no z-density profiles requested',.true.)
     End If

     If (lpzdn) Then
        Call info('z-density printing requested',.true.)
     Else
        Call info('no z-density printing requested',.true.)
     End If

     If (.not.lzdn) Then
        Call info('z-density routines not to be activated',.true.)
        lpzdn=.false.
     End If
  End If

! report vaf

  If (vafsamp > 0 .or. lpvaf) Then
     If (vafsamp > 0) Then
        Write(messages(1),'(a)') 'vaf profiles requested:'
        Write(messages(2),'(2x,a,i10)') 'vaf collection interval ',isvaf
        Write(messages(3),'(2x,a,i10)') 'vaf binsize  ',nsvaf
        Call info(messages,3,.true.)
     Else
        Call info('no vaf collection requested',.true.)
     End If

     If (lpvaf) Then
        Call info('vaf printing requested',.true.)
     Else
        Call info('no vaf printing requested',.true.)
     End If

     If (lvafav) Then
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
        If (.not.(lrdf .or. lzdn .or. ldef .or. l_msd .or. lrsd .or. (mxgana > 0))) Call error(580)
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
        temp=300.0_wp
        Write(message,'(a,1p,e12.4)') 'default simulation temperature (K) ',temp
        Call info(message,.true.)
     End If

     If (.not.lpres) Then ! Simulation pressure
        lpres=.true.
        press=0.0_wp
        Write(message,'(a,1p,e12.4)') 'default simulation pressure (katms) ',press*prsunt
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
     temp = 10.0_wp
     Write(message,'(a,1p,e12.4)') 'default simulation temperature (K) ',temp
     Call info(message,.true.)
  End If

  vel_es2 = vel_es2 * vel_es2 ! square of cutoff velocity for inhomogeneous Langevin thermostat and ttm
  amin = Max (amin, 1)        ! minimum number of atoms for ttm ionic temperature cell

!!! ERROR CHECKS !!!
! Temperature

  If ((.not.ltemp) .or. (nstrun > 0 .and. temp < 1.0_wp)) Call error(380)

! Timestep

  If (.not.lstep) Call error(381)

! check settings in Langevin ensembles

  If ((keyens == 10 .or. keyens == 20 .or. keyens == 30) .and. &
      chi <= zero_plus) Call error(462)

  If (keyens == 15 ) Then
    If (gvar==0 .and. chi_ep <= zero_plus) Call error(462)
    If (Abs(chi_es) <= zero_plus) Then
      Call info('assuming no electronic stopping in inhomogeneous Langevin thermostat',.true.)
    End If
  End If

  If ((keyens == 20 .or. keyens == 30) .and. &
      tai <= zero_plus) Call error(463)

! check settings in ensembles with taut

  If (((keyens >= 11 .and. keyens <= 13) .or. &
       (keyens >= 21 .and. keyens <= 23) .or. &
       (keyens >= 31 .and. keyens <= 33)) .and. taut <= 0.0_wp) Call error(464)

! check settings in ensembles with press

  If ((keyens >= 20 .and. keyens <= 23) .or. &
      (keyens >= 30 .and. keyens <= 33)) Then
     If      (keyens >= 20 .and. keyens <= 23) Then
        If (.not.lpres) Then
           If (lstrext) Then
              press=(strext(1)+strext(5)+strext(9))/3.0_wp
              strext=0.0_wp

              Write(messages(1),'(a)') 'tensorial system pressure specified for an npt ensemble simulation'
              Write(messages(2),'(a)') 'scalar pressure derived from pressure tensor as p = Trace[P]/3'
              Write(messages(3),'(a)') 'tensorial pressure to be zeroed (discarded)'
              Write(messages(4),'(a,1p,e12.4)') 'simulation pressure (katms) ',press*prsunt
              Call info(messages,4,.true.)
           Else
              Call error(387)
           End If
        Else
           If (lstrext) Then
              strext=0.0_wp

              Write(messages(1),'(a)') 'both tensorial and scalar system pressure specified for an npt ensemble simulation'
              Write(messages(2),'(a)') 'tensorial pressure directive is ignored'
              Write(messages(3),'(a)') 'tensorial pressure to be zeroed (discarded)'
              Call info(messages,3,.true.)
           End If
        End If
     Else If (keyens >= 30 .and. keyens <= 33) Then
        If (.not.lstrext) Then
           If (.not.lpres) Call error(387)
        Else
           If (lpres) Then
              Write(messages(1),'(a)') 'both tensorial and scalar system pressure specified for an nst ensemble simulation'
              Write(messages(2),'(a)') 'scalar pressure directive is ignored'
              Call info(messages,2,.true.)

! Define new scalar pressure and zero trace pressure tensor

              press=(strext(1)+strext(5)+strext(9))/3.0_wp
              strext(1)=strext(1)-press
              strext(5)=strext(5)-press
              strext(9)=strext(9)-press
           End If
        End If
     End If
     If (keyens /= 20 .and. keyens /= 30 .and. taup <= 0.0_wp) Call error(466)
  End If

! Two-temperature model: calculate atomic density (if not
! already specified and electron-phonon friction
! conversion factor (to calculate chi_ep from G_ep values)

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
    Diff0 = Diff0*temp*1.0e8_wp
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
           (rcbnd,mxrdf,mxvdw,rvdw,mxmet,rmet,mxter,rcter, &
           mxrgd,imcon,imc_n,cell,xhi,yhi,zhi,             &
           mxgana,mxgbnd1,mxgang1,mxgdih1,mxginv1,         &
           l_str,lsim,l_vv,l_n_e,l_n_r,lzdn,l_n_v,l_ind,   &
           rcut,rpad,rbin,mxstak,                          &
           mxshl,mxompl,mximpl,keyind,                     &
           nstfce,mxspl,alpha,kmaxa1,kmaxb1,kmaxc1,comm)

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
  Integer,           Intent( In    ) :: mxrdf,mxvdw,mxmet,mxter,mxrgd,imcon,mxshl
  Integer,           Intent( InOut ) :: imc_n,mxompl,mximpl,keyind
  Integer,           Intent(   Out ) :: mxgana,mxgbnd1,mxgang1,mxgdih1,mxginv1, &
                                        mxstak,nstfce,mxspl,kmaxa1,kmaxb1,kmaxc1
  Real( Kind = wp ), Intent( In    ) :: xhi,yhi,zhi,rcter
  Real( Kind = wp ), Intent( InOut ) :: rvdw,rmet,rcbnd,cell(1:9)
  Real( Kind = wp ), Intent(   Out ) :: rcut,rpad,rbin,alpha
  Type( comms_type ), Intent( InOut ) :: comm

  Logical                :: carry,safe,la_ana,la_bnd,la_ang,la_dih,la_inv, &
                            lrcut,lrpad,lrvdw,lrmet,lelec,lrdf,lvdw,lmet,l_n_m,lter,l_exp
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
  la_bnd = .false. ; mxgbnd1 = 0
  la_ang = .false. ; mxgang1 = 0
  la_dih = .false. ; mxgdih1 = 0
  la_inv = .false. ; mxginv1 = 0

! electrostatics and no electrostatics, rdf and no rdf, vdw and no vdw,
! metal and no metal, tersoff and no tersoff interactions,
! cutoff and padding, and binsize defaults

  lelec = .false.
! l_n_e is now first determined in scan_field l_n_e = (.false.)

  lrdf  = (mxrdf > 0)
  l_n_r = .not.lrdf

  lvdw  = (mxvdw > 0)
  l_n_v = .false.
  lrvdw = .false. ! Even though it rvdw may have been read from TABLE

  lmet  = (mxmet > 0)
  l_n_m = .not.lmet
  lrmet = (rmet > 1.0e-6_wp)

  lter  = (mxter > 0)

  lrcut = .false.
  rcut  = 0.0_wp

  lrpad = .false.
  rpad  = 0.0_wp

  rbin  = rbin_def

! Frequency of the SPME k-space evaluation

  nstfce = -1 ! None defined

! Ewald/Poisson Solver sum parameters defaults

  mxspl = 0
  alpha = 0.0_wp
  kmaxa1 = 0
  kmaxb1 = 0
  kmaxc1 = 0

!  induce = .false.

! default number of steps and expansion option

  nstrun = 0
  l_exp = .false.

! default stack size

  mxstak = 1

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
        rcut = Abs(word_2_real(word))
        lrcut = (rcut > zero_plus) ! if zero or nothing is entered

! read real space cut off

     Else If (word(1:3) == 'pad' .or. word(1:4) == 'rpad') Then

        lrpad = .true.
        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        rpad = Max(rpad,Abs(word_2_real(word)))
        lrpad = (rpad > zero_plus) ! if zero or nothing is entered

! read vdw cutoff

     Else If (word(1:4) == 'rvdw') Then

        lrvdw=.true.
        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        If (rvdw > 1.0e-6_wp) Then
           rvdw = Min(rvdw,word_2_real(word))
        Else
           rvdw = Abs(word_2_real(word))
        End If
        lrvdw = (rvdw > zero_plus) ! if zero or nothing is entered

! read binsize option

     Else If (word(1:7) == 'binsize') Then

        Call get_word(record,word)
        rbin = Abs(word_2_real(word))

! read dpd ensembles option

     Else If (word(1:8) == 'ensemble') Then

        Call get_word(record,word)
        If (word(1:3) == 'nvt') Then
           Call get_word(record,word)
           If (word(1:3) == 'dpd') Then
              If      (word(1:5) == 'dpds1') Then
                 keydpd = 1
              Else If (word(1:5) == 'dpds2') Then
                 keydpd = 2
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

        mxstak = Nint(Abs(word_2_real(word)))

! read MSD option

     Else If (word(1:6) == 'msdtmp') Then

        l_msd = .true.

! read VAF option and sample frequency and binsize - defaults in greenkubo_module

     Else If (word(1:3) == 'vaf') Then

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        isvaf = Abs(Nint(word_2_real(word,0.0_wp)))
        If (isvaf == 0) isvaf=50

        Call get_word(record,word)
        If (word(1:3) == 'bin' .or. word(1:5) == 'size') Call get_word(record,word)
        If (word(1:3) == 'bin' .or. word(1:5) == 'size') Call get_word(record,word)
        nsvaf = Abs(Nint(word_2_real(word,0.0_wp)))

        If (nsvaf == 0) nsvaf=Merge(2*isvaf,100,isvaf >= 100)
        vafsamp = Ceiling(Real(nsvaf,wp)/Real(isvaf,wp))

! read DL_POLY_2/Classic delr Verlet shell strip cutoff option (compatibility)
! as DL_POLY_4 real space cutoff padding option

     Else If (word(1:4) == 'delr') Then

        lrpad = .true.
        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        rpad = Max(rpad,0.25_wp*Abs(word_2_real(word)))
        lrpad = (rpad > zero_plus) ! if zero or nothing is entered

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

!     Else If (word(1:6) == 'induce') Then
!
!        induce=.true.
!        Call get_word(record,word1)
!        politer = Min(500,Nint(Abs(word_2_real(word1))))
!
!        Call get_word(record,word2)
!        convcrit = Abs(word_2_real(word2)) !Min(0.001_wp,(Abs(word_2_real(word2))))

!     Else If (word(1:4) == 'gear') Then
!
!         gear=.true.
!         Call get_word(record,word)
!         numcof = Max(0,1+Nint(Abs(word_2_real(word))))
!         numcof = Min(numcof,7)
!
!     Else If (word(1:4) == 'aspc') Then
!
!         aspc=.true.
!         Call get_word(record,word)
!         numcof = Max(0,1+Nint(Abs(word_2_real(word))))
!         numcof = Min(numcof,7)
!
!     Else If (word(1:5) == 'lstsq') Then
!
!         lstsq=.true.
!         Call get_word(record,word)
!         numcof = Max(0,1+Nint(Abs(word_2_real(word))))

! read "no vdw", "no elec" and "no str" options

     Else If (word(1:5) == 'polar') Then

        Call get_word(record,word)
        If (word(1:6) == 'scheme' .or. word(1:4) == 'type') Call get_word(record,word)
        If (word(1:6) == 'scheme' .or. word(1:4) == 'type') Call get_word(record,word)
        If (word(1:6) == 'charmm' .and. mxshl > 0) keyind=1

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
           mxgbnd1 = Max(mxgbnd1,mxgana)
           mxgang1 = Max(mxgang1,mxgana)
           mxgdih1 = Max(mxgdih1,mxgana)
           mxginv1 = Max(mxginv1,mxgana)

           Call get_word(record,word) ! AB: for "rbnd"/"rmax"/"max"/figure
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           rcbnd=Max(rcbnd,word_2_real(word,0.0_wp))
        Else If (akey == 'bon') Then
           la_bnd = .true.

           mxgbnd1 = Max(mxgbnd1,Abs(Nint(word_2_real(word))))

           Call get_word(record,word) ! AB: for "rbnd"/"rmax"/"max"/figure
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           rcbnd=Max(rcbnd,word_2_real(word,0.0_wp))
        Else If (akey == 'ang') Then
           la_ang = .true.

           mxgang1 = Max(mxgang1,Abs(Nint(word_2_real(word))))
        Else If (akey == 'dih') Then
           la_dih = .true.

           mxgdih1 = Max(mxgdih1,Abs(Nint(word_2_real(word))))
        Else If (akey == 'inv') Then
           la_inv = .true.

           mxginv1 = Max(mxginv1,Abs(Nint(word_2_real(word))))
        End If

! read rdf calculation option

     Else If (word(1:3) == 'rdf') Then

        lrdf = .true.
        l_n_r = .not.lrdf

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

        ! variable electron-phonon coupling constant (chi_ep) based on
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
          l_errors_jack = .TRUE.
       Else
          l_errors_block = .TRUE.
       End If
     End If

  End Do

! in case of intramolecular distribution analysis

  If (la_ana) Then
     If (mxgana > 0) Then
        mxgbnd1 = Max(mxgbnd1,mxgana)
        mxgang1 = Max(mxgang1,mxgana)
        mxgdih1 = Max(mxgdih1,mxgana)
        mxginv1 = Max(mxginv1,mxgana)
     End If

! switch indicators for set_bounds

     If (la_bnd) Then
        If (mxgbnd1 == 0) mxgbnd1 = -1
        rcbnd=Max(rcbnd,rcbnd_def)
     End If
     If (la_ang .and. mxgang1 == 0) mxgang1 = -1
     If (la_dih .and. mxgdih1 == 0) mxgdih1 = -1
     If (la_inv .and. mxginv1 == 0) mxginv1 = -1

! mxgana by construction equals the largest possible grid
! or 1 (positive) as an indicator for analysis

     mxgana=Max(1,mxgbnd1,mxgang1,mxgdih1,mxginv1)
  End If

! Sort electrostatics

  If (lelec) Then
     If (l_n_e) lelec = .not.l_n_e
  Else
     l_n_e = .true.
  End If

! reinitialise multipolar electrostatics indicators

  If (l_n_e) Then
     mximpl = 0
     mxompl = 0
     keyind = 0
  End If

! Sort vdw

  If (lvdw) Then
     If (.not.lrvdw) Then
        lrvdw = (rvdw > 1.0e-6_wp)
        rvdw = Min(rvdw,Max(rcut,rcut_def))
     End If

     If (l_n_v) lvdw = .not.l_n_v
  Else
     l_n_v = .true.
  End If

! Sort rcut as the maximum of all valid cutoffs

  rcut=Max(rcut,rvdw,rmet,rkim,2.0_wp*Max(rcter,rcbnd)+1.0e-6_wp)

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

! rcut MUST be >= rcut_def

           If (rcut < rcut_def) rcut=rcut_def

! define cut

           cut=rcut+1.0e-6_wp

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

                 tol = Sqrt(Abs(Log(eps0*rcut)))
                 alpha = Sqrt(Abs(Log(eps0*rcut*tol)))/rcut
                 tol1 = Sqrt(-Log(eps0*rcut*(2.0_wp*tol*alpha)**2))

                 fac = 1.0_wp
                 If (imcon == 4 .or. imcon == 5 .or. imcon == 7) fac = 2.0_wp**(1.0_wp/3.0_wp)

                 kmaxa1 = 2*Nint(0.25_wp + fac*celprp(7)*alpha*tol1/pi)
                 kmaxb1 = 2*Nint(0.25_wp + fac*celprp(8)*alpha*tol1/pi)
                 kmaxc1 = 2*Nint(0.25_wp + fac*celprp(9)*alpha*tol1/pi)

! rcut is needed directly for the SPME and it MUST exist

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

! rcut is not needed directly for the SPME but it's needed
! for the link-cell division of the domains
! let's not fail here if no cutoff is specified

              End If

! Get default spline order or one driven by multipolar sums if none is specified
! Only even order splines are allowed so pick the even=odd+1 if resulting in odd!!!

              If (mxspl == 0) Then
                 mxspl  = mxspl_def+mxompl
                 mxspl2 = mxspl
              Else
                 mxspl  = Max(mxspl,mxspl_min)
                 mxspl2 = mxspl+mxompl
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
                    eps=Abs(word_2_real(word))
                 End If

                 If (word(1:6) == 'maxits') Then  ! max number of iteration
                    Call get_word(record,word)
                    mxitcg=Nint(Abs(word_2_real(word)))
                 End If

                 If (word(1:7) == 'jmaxits') Then ! max number Jacobian iterations
                    Call get_word(record,word)
                    mxitjb=Nint(Abs(word_2_real(word)))
                 End If

                 Call get_word(record,word)
              End Do

! Check for undefined and ill defined parameters

!             0.1 Angs <= delta=1/alpha <= Min(3 Angs,rcut/3) - 3 grid points within a link-cell
              If (alpha > 10.0_wp)                       alpha = 10.0_wp
              If (alpha < 1.0_wp/Min(3.0_wp,cut/3.0_wp)) alpha = 1.0_wp/Min(3.0_wp,cut/3.0_wp)

              If (mxitcg == 0) mxitcg = 1000 ! default
              If (mxitjb == 0) mxitjb = 1000 ! default

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

! Sort rvdw

        If ((.not.lrvdw) .and. lvdw) Then
           If (lrcut) Then
              lrvdw=.true.
              rvdw=rcut
           Else
              Call error(402)
           End If
        End If

! Sort rmet

        If ((.not.lrmet) .and. lmet) Then
           If (lrcut .or. lrvdw) Then
              lrmet=.true.
              rmet=Max(rcut,rvdw)
           Else
              Call error(382)
           End If
        End If

! Sort rcut by a reset sequence
! rcut may be >= rcut_def but lrcut may still be .false.
! mxspl = 0 is an indicator for no SPME or Poisson Solver electrostatics in CONTROL

        If (mxspl /= 0) Then ! SPME or Poisson Solver

! (1) to Max(rcut,Max(cell_width*mxspl/kmax),mxspl*delta) satisfying SPME b-splines
! propagation width or the Poisson Solver extra halo relation to cutoff
! delta=1/alpha is the grid spacing and mxspl is the grid length needed for the
! 3 haloed stencil of differentiation

           If (.not.lrcut) Then
              lrcut=.true.

              Call dcell(cell,celprp)
              rcut=Max( rcut, Merge(Real(mxspl,wp)/alpha,                          &
                                    Max(celprp(7)*Real(mxspl,wp)/Real(kmaxa1,wp),  &
                                        celprp(8)*Real(mxspl,wp)/Real(kmaxb1,wp),  &
                                        celprp(9)*Real(mxspl,wp)/Real(kmaxc1,wp)), &
                                    itmp == 0) )
           End If

! Reset rvdw, rmet and rcut when only tersoff potentials are opted for

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r) Then
              rvdw=0.0_wp
              rmet=0.0_wp
              If (.not.l_str) Then
                 If (mxrgd == 0) Then ! compensate for Max(Size(RBs))>rvdw
                    rcut=2.0_wp*Max(rcbnd,rcter)+1.0e-6_wp
                 Else
                    rcut=Max(rcut,2.0_wp*Max(rcbnd,rcter)+1.0e-6_wp)
                 End If
              End If
           End If

        Else

! no SPME electrostatics is specified but rcut is still needed for
! domain decompositioning and link-celling
! It is needed for the rest of the types of electrostatics

           If ((.not.lrcut) .and. lelec) Call error(382)

! So there is rcut and some kind of electrostatics(-: or neither

! Reset rcut to something sensible if sensible is an option

           If ( ((.not.lrcut) .or. (.not.l_str)) .and. &
                (lrvdw .or. lrmet .or. lter .or. kimim /= ' ') ) Then
              lrcut=.true.
              If (mxrgd == 0) Then ! compensate for Max(Size(RBs))>rvdw
                 rcut=Max(rvdw,rmet,rkim,2.0_wp*Max(rcbnd,rcter)+1.0e-6_wp)
              Else
                 rcut=Max(rcut,rvdw,rmet,rkim,2.0_wp*Max(rcbnd,rcter)+1.0e-6_wp)
              End If
           End If

! Reset rvdw and rmet when only tersoff potentials are opted for and
! possibly reset rcut to 2.0_wp*rcter+1.0e-6_wp (leaving room for failure)

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r .and. kimim == ' ') Then
              rvdw=0.0_wp
              rmet=0.0_wp
              If (.not.l_str) Then
                 lrcut=.true.
                 If (mxrgd == 0) Then ! compensate for Max(Size(RBs))>rvdw
                    rcut=2.0_wp*Max(rcbnd,rcter)+1.0e-6_wp
                 Else
                    rcut=Max(rcut,2.0_wp*Max(rcbnd,rcter)+1.0e-6_wp)
                 End If
              End If
           End If

! rcut must exist

           If (.not.lrcut) Call error(382)

! define cut

           cut=rcut+1.0e-6_wp

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

! Sort rmet=rcut if metal interactions are in play, even if
! they are defined by EAM since rmet can be /= rcut in such
! instances, this can break the NLAST check in metal_ld_set_halo

        If (lmet) rmet = rcut

! Sort rvdw=rcut if VDW interactions are in play

        If (lvdw .and. rvdw > rcut) rvdw = rcut

! Sort rbin as now rcut is already pinned down

        If (rbin < 1.0e-05_wp .or. rbin > rcut/4.0_wp) rbin = Min(rbin_def,rcut/4.0_wp)

        carry=.false.

     End If

  End Do

  If (comm%idnode == 0) Close(Unit=nread)

! Enforce VV for DPD thermostat

  If (keydpd > 0) l_vv = .true.

! When not having dynamics or prepared to terminate
! expanding and not running the small system prepare to exit gracefully

  l_trm = (l_exp .and. nstrun == 0)
  If (((.not.lsim) .or. l_trm) .and. lrpad) rpad=0.0_wp

  l_errors_block = l_errors_block .and. lrdf
  l_errors_jack = l_errors_jack .and. lrdf
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

Subroutine scan_control_io(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning the I/O options in the control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
! amended   - i.j.bush october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical                :: carry,safe
  Character( Len = 200 ) :: record,record1
  Character( Len = 40  ) :: word,word1
  Real( Kind = wp )      :: tmp
  Type( comms_type ), Intent( InOut ) :: comm

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
                 Call io_nc_set_real_precision( Precision( 1.0 ), Range( 1.0 ), err_r )
              Else
                 ! Use 64-bit quantities in output for real numbers
                 Call info('I/O write method: parallel by using netCDF in 64-bit format',.true.)
                 Call io_nc_set_real_precision( Precision( 1.0d0 ), Range( 1.0d0 ), err_r )
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
