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
           ndump,pdplnc,timjob,timcls)

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

  Use kinds_f90
  Use comms_module,      Only : idnode
  Use setup_module
  Use config_module,     Only : sysname
  Use mpoles_module,     Only : thole
  Use dpd_module,        Only : keydpd,gamdpd
  Use langevin_module,   Only : l_lan,l_gst,langevin_allocate_arrays
  Use parse_module
  Use bonds_module,      Only : rcbnd
  Use vdw_module,        Only : ld_vdw,ls_vdw,mxtvdw
  Use metal_module,      Only : ld_met,ls_met,tabmet
  Use poisson_module,    Only : eps,mxitcg,mxitjb
  Use msd_module,        Only : l_msd
  Use defects1_module,   Only : l_dfx
  Use rdf_module,        Only : l_block,l_jack,num_blocks
  Use greenkubo_module
  Use ttm_module

  Use kinetic_module,  Only : l_vom
  Use plumed_module,   Only : l_plumed, plumed_input, plumed_log, &
                              plumed_precision, plumed_restart
  Use development_module

  Implicit None

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
                                             rdef,rrsd,pdplnc,           &
                                             timjob,timcls


  Logical                                 :: limp,lvv,lens,lforc,     &
                                             ltemp,l_0,lpres,lstrext, &
                                             lstep,lplumed,safe,      &
                                             l_timjob,l_timcls

  Character( Len = 200 )                  :: record
  Character( Len = 40  )                  :: word,word1,word2,word3,akey

  Integer                                 :: i,j,k,itmp,nstana,grdana,grdbnd,grdang, &
                                             grddih,grdinv,nstall

  Real( Kind = wp )                       :: rcell(1:9),rcut1,rpad1,rvdw1,tmp,eps0,tol,rcb_d,prmps(1:4)


! initialise system control variables and their logical switches

! default expansion option

  l_exp = .false.
  nx    = 1
  ny    = 1
  nz    = 1

! defaults for direct evaluation, force-shifting of VDW interactions
! and type of mixing for undefined cross interaction of certain type
!
! ld_vdw = .false. ! (initialised in vdw_module)
! ls_vdw = .false. ! (initialised in vdw_module)
  mxtvdw = 0       ! (initialised in vdw_module)
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
!  l_vom = .true. ! initialised in kinetic_module

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

  timjob = 0.0_wp ; l_timjob=.false.
  timcls = 0.0_wp ; l_timcls=.false.

! major cutoff, padding and vdw cutoff defaults

  rcut1 = 0.0_wp
  rpad1 = 0.0_wp
  rvdw1 = 0.0_wp

! open the simulation control file

  If (idnode == 0) Open(Unit=nread, File = Trim(control), Status = 'old')

! read simulation control name

  Call get_line(safe,nread,sysname)
  If (.not.safe) Go To 1000
  Call strip_blanks(sysname)

  If (.not.safe) Go To 1000

  If (idnode == 0) Write(nrite,"(/,3(1x,130('*'),/),1x,     &
     & 24('*'),5x,a72,5x,24('*'),/,3(1x,130('*'),/),/,/,1x, &
     & 'SIMULATION CONTROL PARAMETERS')") sysname

! read and process directives from CONTROL file

  Do

     Call get_line(safe,nread,record)
     If (.not.safe) Go To 1000
     Call lower_case(record)
     Call get_word(record,word)

! record is commented out

     If (word(1:1) == '#' .or. word(1:1) == ' ') Then

! read DEVELOPMENT options

     Else If (word(1:5) == 'l_scr') Then

!        l_scr = .true. ! done in scan_development
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% OUTPUT redirected to the default output (screen) !!! %%%"

     Else If (word(1:6) == 'l_fast') Then

!        l_fast = .true. ! done in scan_development
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% speed up by avoiding global safety checks !!! %%%"

     Else If (word(1:5) == 'l_eng') Then

        l_eng = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% OUTPUT contains an extra last line with E_tot !!! %%%"

     Else If (word(1:6) == 'l_rout') Then

        l_rout = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% REVIVE writing in ASCII opted !!! %%%"

     Else If (word(1:5) == 'l_rin') Then

        l_rin = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% REVOLD reading in ASCII opted !!! %%%"

     Else If (word(1:5) == 'l_org') Then

        l_org = .true.
        l_trm  = .true.

        If (idnode == 0) Write(nrite,"(2(/,1x,a))")                                                   &
           "%%% translate CONFIG along a vector into CFGORG after reading input & terminate !!! %%%", &
           "%%% vector and config level read as follows: %%%"

        Call get_word(record,word)
        xorg = word_2_real(word)
        Call get_word(record,word)
        yorg = word_2_real(word)
        Call get_word(record,word)
        zorg = word_2_real(word)

        Call get_word(record,word)
        lvcforg = Min( Int(Abs(word_2_real(word,0.0_wp))) , levcfg)

        If (idnode == 0) Then
           Write(nrite,"(1x,a)")          '%%% '
           Write(nrite,"(1x,a,3f10.3,a)") '%%% vector(x,y,x) ', xorg, yorg, zorg, ' %%%'
           Write(nrite,"(1x,a,i0,a)")     '%%% CFGORG level ', lvcforg, ' %%%'
        End If

     Else If (word(1:5) == 'l_scl') Then

        If (idnode == 0) Write(nrite,"(2(/,1x,a))")                                 &
           "%%% rescale CONFIG to CFGSCL, after reading input & terminate !!! %%%", &
           "%%% config level and new cell vectors to rescale to (read in a CONFIG-like manner): %%%"

        Call get_word(record,word)
        lvcfscl = Min( Int(Abs(word_2_real(word,0.0_wp))) , levcfg)

        itmp=0
        Do i=1,3
           Call get_line(safe,nread,record)
           Do j=1,3
              Call get_word(record,word)
              itmp=itmp+1
              cels(itmp)=word_2_real(word)
           End Do
        End Do

        Call invert(cels,rcell,tmp)

        If (idnode == 0) Then
           Write(nrite,"(1x,a)")             '%%% '
           Write(nrite,"(1x,a,i0,a)")        '%%% CFGSCL level ', lvcfscl, ' %%%'
           Write(nrite,"(1x,a,3f20.10,a)")   '%%% ', cels(1:3), ' %%%'
           Write(nrite,"(1x,a,3f20.10,a)")   '%%% ', cels(4:6), ' %%%'
           Write(nrite,"(1x,a,3f20.10,a)")   '%%% ', cels(7:9), ' %%%'
           Write(nrite,"(1x,a)")             '%%% '
           Write(nrite,"(1x,a,1p,g22.12,a)") '%%% CFGSCL volume ', tmp, '%%%'
        End If

        If (tmp > zero_plus) Then
           l_scl = .true.
           l_trm  = .true.
        Else
           If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% OPTION ABORTED DUE TO ZERO VOLUME !!! %%%"
           l_trm  = .true.
        End If

     Else If (word(1:5) == 'l_his') Then

        l_his = .true.
        l_trm = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% generate HISTORY after reading input & terminate !!! %%%"

     Else If (word(1:5) == 'l_tim') Then

!        l_tim = .true.  ! done in scan_development
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% generate detailed timing !!! %%%"

     Else If (word(1:5) == 'l_tor') Then

        l_tor = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% Turn off production of REVCON & REVIVE !!! %%%"

     Else If (word(1:5) == 'l_trm') Then

        l_trm = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% Terminate gracefully before initialisation !!! %%%"

     Else If (word(1:5) == 'l_dis') Then

        l_dis = .true.
        r_dis = Min( r_dis , word_2_real(word,0.1_wp) )
        If (idnode == 0) Write(nrite,"(2(/,1x,a),1p,e12.4)")                            &
           "%%% Turn on the check on minimum separation distance between VNL pairs at re/start !!! %%%", &
           "%%% separation criterion (Angstroms) %%% ", r_dis

! read VDW options

     Else If (word(1:3) == 'vdw') Then

        Call get_word(record,word1)

        If      (word1(1:6) == 'direct') Then

! direct evaluation option

           ld_vdw = .true.
           If (idnode == 0) Write(nrite,"(/,1x,a)") "vdw direct option on"

        Else If (word1(1:6) == 'mixing') Then

! mixing type keywords

           If (idnode == 0) Write(nrite,"(3(/,1x,a))") "vdw cross terms mixing opted (for undefined mixed potentials)", &
                                                       "mixing is limited to potentials of the same type only",         &
                                                       "mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)"

           Call get_word(record,word2)

           If      (word2(1:4) == 'lore') Then

              mxtvdw = 1
              If (idnode == 0) Write(nrite,"(1x,a)") &
  "type of mixing selected - Lorentz–Berthelot :: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2"

           Else If (word2(1:4) == 'fend') Then

              mxtvdw = 2
              If (idnode == 0) Write(nrite,"(1x,a)") &
"type of mixing selected - Fender-Halsey :: e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2"

           Else If (word2(1:4) == 'hoge') Then

              mxtvdw = 3
              If (idnode == 0) Write(nrite,"(1x,a)") &
"type of mixing selected - Hogervorst (good hope) :: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)"

           Else If (word2(1:4) == 'halg') Then

              mxtvdw = 4
              If (idnode == 0) Write(nrite,"(1x,a)") &
"type of mixing selected - Halgren HHG :: e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)"

           Else If (word2(1:4) == 'wald') Then

              mxtvdw = 5
              If (idnode == 0) Write(nrite,"(1x,a)") &
"type of mixing selected - Waldman–Hagler :: e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ; s_ij=[(s_i^6+s_j^6)/2]^(1/6)"

           Else If (word2(1:4) == 'tang') Then

              mxtvdw = 6
              If (idnode == 0) Write(nrite,"(1x,a,/,1x,a)") &
"type of mixing selected - Tang-Toennies :: e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13 ;", &
"                                           s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)"

           Else If (word2(1:4) == 'func') Then

              mxtvdw = 7
              If (idnode == 0) Write(nrite,"(1x,a,a,/,1x,a)")                &
"type of mixing selected - Functional :: e_ij=3 * (e_i*e_j)^(1/2) * ",       &
"(s_i*s_j)^3 / SUM_L=0^2{[(s_i^3+s_j^3)^2 / (4*(s_i*s_j)^L)]^(6/(6-2L))} ;", &
"                                        s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}"

           Else

              Call strip_blanks(record)
              If (idnode == 0) Write(nrite,"(/,/,4a)") &
                 word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),word2(1:Len_Trim(word2)+1),record
              Call error(3)

           End If

        Else If (word1(1:5) == 'shift') Then

! force-shifting option

           ls_vdw = .true.
           If (idnode == 0) Write(nrite,"(/,1x,a)") "vdw force-shifting option on"

        Else

           Call strip_blanks(record)
           If (idnode == 0) Write(nrite,"(/,/,3a)") word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call error(3)

        End If

     Else If (word(1:5) == 'metal') Then

        Call get_word(record,word)
        If      (word(1:6) == 'direct') Then

! read metal direct evaluation option

           If (idnode == 0) Write(nrite,"(/,1x,a)") "metal direct option on"
           If (tabmet > 0) Then
              Call warning(480,0.0_wp,0.0_wp,0.0_wp)
           Else
              ld_met = .true.
           End If

        Else If (word(1:7) == 'sqrtrho') Then

! read metal sqrtrho interpolation option for EAM embeding function in TABEAM

           If (idnode == 0) Write(nrite,"(/,1x,a)") "metal sqrtrho option on"
           If (tabmet > 0) Then
              ls_met = .true.
           Else
              Call warning(490,0.0_wp,0.0_wp,0.0_wp)
           End If

        End If

! read slab option (dealt with in scan_control<-set_bounds,
! affecting map_domains<-set_bounds)

     Else If (word(1:4) == 'slab') Then

        If (idnode == 0) Write(nrite,"(/,1x,a)") "slab option on"

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
        If (idnode == 0) Write(nrite,"(/,1x,'system expansion opted',9x,3i5)") nx,ny,nz

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

        If (idnode == 0) Write(nrite,"(/,1x,'impact option on', &
           & /,1x,'particle (index)',15x,i10,                   &
           & /,1x,'timestep (steps)',15x,i10,                   &
           & /,1x,'energy   (keV)  ',18x,1p,e12.4,              &
           & /,1x,'v-r(x,y,z)',1p,3e12.4)") imd,tmd,emd,vmx,vmy,vmz

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

        If (idnode == 0) &
           Write(nrite,"(/,1x,'radomisation seeds supplied',/,1x,'(seed1,seed2,seed3)  ',10x,3i5)") seed

! read temperature

     Else If (word(1:4) == 'temp') Then

        ltemp = .true.
        Call get_word(record,word)
        temp = Abs(word_2_real(word))
        If (idnode == 0) Write(nrite,"(/,1x,'simulation temperature (K)  ',6x,1p,e12.4)") temp

! read zero temperature optimisation

     Else If (word(1:4) == 'zero') Then

        lzero = .true.

! Check defaults

        Call get_word(record,word)
        l_0 = (word(1:4) == 'fire')
        nstzero = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        If (word(1:5) == 'every') Call get_word(record,word)
        nstzero = Max(nstzero,Abs(Nint(word_2_real(word,0.0_wp))))

        If (idnode == 0) Write(nrite,"(/,1x,'zero K optimisation on (during equilibration)', &
           & /,1x,'temperature regaussing interval',i10)") nstzero

        If (l_0) Then
           If (idnode == 0) &
  Write(nrite,"(1x,a)") 'fire option on - actual temperature will reset to 10 Kelvin if no target tempreature is specified'
        Else
           ltemp  = .true.
           temp = 10.0_wp
           If (idnode == 0) &
  Write(nrite,"(1x,a)") 'fire option off - actual temperature reset to 10 Kelvin'
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


           If (idnode == 0) Then
              Write(nrite,"(/,1x,'simulation pressure tensor (katms)'/)")
              Write(nrite,"(3f20.10)") strext
           End If

! convert from katms to internal units of pressure

           strext = strext/prsunt

        Else

           lpres=.true.

           press = word_2_real(word)

           If (idnode == 0) Write(nrite,"(/,1x,'simulation pressure (katms)  ',5x,1p,e12.4)") press

! convert from katms to internal units of pressure

           press = press/prsunt

        End If

! read restart

     Else If (word(1:7) == 'restart') Then

        Call get_word(record,word)

        If (word(1:7) == 'noscale' .or. word(1:7) == 'unscale') Then

           keyres = 3
           If (idnode == 0) Write(nrite,"(/,1x,'unscaled restart requested (starting a new simulation)')")

        Else If (word(1:5) == 'scale') Then

           keyres = 2
           If (idnode == 0) Write(nrite,"(/,1x,'scaled restart requested (starting a new simulation)')")

        Else

           keyres = 1
           If (idnode == 0) Write(nrite,"(/,1x,'restart requested (continuing an old simulation)')")

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
           If (idnode == 0) Write(nrite,"(/,/,3a)") word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
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
        If (idnode == 0) Write(nrite,"(/,1x,'selected number of timesteps',3x,i10)") nstrun

! read number of equilibration timesteps

     Else If (word(1:5) == 'equil') Then

        Call get_word(record,word)
        If (word(1:5) == 'steps') Call get_word(record,word)
        nsteql = Abs(Nint(word_2_real(word)))
        If (idnode == 0) Write(nrite,"(/,1x,'equilibration period (steps)',3x,i10)") nsteql

! read collection option

     Else If (word(1:7) == 'collect') Then

        leql = .false.
        If (idnode == 0) Write(nrite,"(/,1x,'equilibration included in overall averages')")

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
           If (idnode == 0) Then
              Write(nrite,"(/,1x,'pseudo thermostat attached to MD cell boundary')")
              If      (keypse == 0) Then
                      Write(nrite,'(1x,a)') "thermostat control: Langevin + direct temperature scaling"
              Else If (keypse == 1) Then
                      Write(nrite,'(1x,a)') "thermostat control: Langevin temperature scaling"
              Else If (keypse == 2) Then
                      Write(nrite,'(1x,a)') "thermostat control: gaussian temperature scaling"
              Else If (keypse == 3) Then
                      Write(nrite,'(1x,a)') "thermostat control: direct temperature scaling"
              End If
              Write(nrite,"(1x,'thermostat thickness (Angs) ',6x,1p,e12.4)") tmp
           End If

           If (width/4.0_wp > tmp .and. tmp >= wthpse) Then
              wthpse = tmp
           Else
              If (idnode == 0) Write(nrite,"(1x,'thermostat thickness insufficient - reset to 2 Angs')")
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
        If (idnode == 0) Write(nrite,"(1x,'thermostat temperature (K)',8x,1p,e12.4)") tmppse

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
           If (idnode == 0) Write(nrite,"(/,/,1x,4a)") word2(1:Len_Trim(word2)+1),' ',word(1:Len_Trim(word)+1),record
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
           If (idnode == 0) Write(nrite,                                &
              & "(/,1x,'minimisation option on (during equilibration)', &
              &   /,1x,'minimisation criterion        ',1x,a8,          &
              &   /,1x,'minimisation frequency (steps)',1x,i10,         &
              &   /,1x,'minimisation tolerance        ',4x,1p,e12.4)")  &
              word1(1:8),nstmin,min_tol(1)
           If (min_tol(2) > zero_plus .and. idnode == 0) Write(nrite,   &
              & "(  1x,'minimisation CGM step         ',4x,1p,e12.4)") min_tol(2)
        Else
           If (idnode == 0) Write(nrite,                               &
              & "(/,1x,'optimisation at start',                        &
              &   /,1x,'optimisation criterion        ',1x,a8,         &
              &   /,1x,'optimisation tolerance        ',4x,1p,e12.4)") &
              word1(1:8),min_tol(1)
           If (min_tol(2) > zero_plus .and. idnode == 0) Write(nrite,  &
              & "(  1x,'optimisation CGM step         ',4x,1p,e12.4)") min_tol(2)
        End If



! read regauss option

     Else If (word(1:6) == 'regaus') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        nstgaus = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        ltgaus =.true.
        If (idnode == 0) Write(nrite,"(/,1x,'regauss temperature on (during equilibration)', &
           & /,1x,'temperature regaussing interval',i10)") nstgaus

! read temperature scaling option

     Else If (word(1:5) == 'scale') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        nstscal = Max(1,Abs(Nint(word_2_real(word,0.0_wp))))

        ltscal =.true.
        If (idnode == 0) Write(nrite,"(/,1x,'temperature scaling on (during equilibration)', &
           & /,1x,'temperature scaling interval',3x,i10)") nstscal

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
           If (idnode == 0) Then
  Write(nrite,"(/,1x,'CHARMM polarisation scheme selected with optional atomic thole dumping of ',f5.2)") thole
              If (mximpl == 0) &
  Write(nrite,"(1x,a)") "*** warning - scheme deselected due to switched off electrostatics !!! ***"
              If (mxshl == 0) &
  Write(nrite,"(1x,a)") "*** warning - scheme disabled due to lack of core-shell defined interatcions !!! ***"
           End If

           If (mximpl == 0 .or. mxshl == 0) Then
!              keyind=0 ! done in scan_control
           Else
              lecx = .true. ! enable extended coulombic exclusion
              If (idnode == 0) &
  Write(nrite,"(1x,'Extended Coulombic eXclusion activated for CHARMM polarisation')")
           End If
        End If

! read integration flavour

     Else If (word(1:8) == 'integrat') Then

        Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:8) == 'leapfrog') lvv=.false.

! keydpd detected in scan_control

        If (keydpd > 0 .and. (lvv .neqv. l_vv) .and. idnode == 0) Write(nrite,"(/,1x,a)") &
           "*** warning - Leapfrog Verlet selected integration defaulted to Velocity Verlet for DPD thermostats !!! ***"

! read ensemble

     Else If (word(1:8) == 'ensemble') Then

        If (idnode == 0) Then
           If (l_vv) Then
              Write(nrite,"(/,1x,'Integration : Velocity Verlet')")
           Else
              Write(nrite,"(/,1x,'Integration : Leapfrog Verlet')")
           End If
        End If

        Call get_word(record,word)

        If      (word(1:3) == 'nve' .or. word(1:3) == 'pmf') Then

           keyens = 0

           If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVE (Microcanonical)')")

           If (lens) Call error(414)
           lens=.true.

        Else If (word(1:3) == 'nvt') Then

           Call get_word(record,word)

           If      (word(1:5) == 'evans') Then

              keyens = 1

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT Evans (Isokinetic)', &
                 & /,1x,'Gaussian temperature constraints in use')")

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:4) == 'lang') Then

              keyens = 10
              If (.not.l_vv) l_lan = .true.

              Call get_word(record,word)
              chi = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT Langevin (Stochastic Dynamics)', &
                 & /,1x,'thermostat friction     (ps^-1)',3x,1p,e12.4)") chi

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:5) == 'ander') Then

              keyens = 11

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              soft = Abs(word_2_real(word))
              If (soft > 1.0_wp) soft=1.0_wp/soft

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT Andersen', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,       &
                 & /,1x,'softness        (dimensionless)',3x,1p,e12.4)") taut,soft

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              keyens = 12

              Call get_word(record,word)
              taut = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT Berendsen', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4)") taut

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              keyens = 13

              Call get_word(record,word)
              taut = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT Nose-Hoover', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4)") taut

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'gst') Then

              keyens = 14
              l_gst = .true.

              Call get_word(record,word)
              taut = Abs(word_2_real(word))

              Call get_word(record,word)
              gama = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT gentle stochastic thermostat', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                         &
                 & /,1x,'friction on thermostat  (ps^-1)',3x,1p,e12.4)") taut,gama

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

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)', &
                 & /,1x,'e-phonon friction       (ps^-1)',3x,1p,e12.4, &
                 & /,1x,'e-stopping friction     (ps^-1)',3x,1p,e12.4, &
                 & /,1x,'e-stopping velocity   (A ps^-1)',3x,1p,e12.4)") chi_ep,chi_es,vel_es2

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'dpd') Then

              If (idnode == 0) Write(nrite,"(1x,a)") "Ensemble : NVT dpd (Dissipative Particle Dynamics)"

! keydpd determined in scan_control

              If      (keydpd == 1) Then
                 keyens = 0 ! equivalence to doing NVE with some extra fiddling before VV(0)
                 If (idnode == 0) Write(nrite,"(1x,a)") "Ensemble type : Shardlow's first order splitting (S1)"
              Else If (keydpd == 2) Then
                 keyens = 0 ! equivalence to doing NVE with some extra fiddling before VV(0) and after VV(1)
                 If (idnode == 0) Write(nrite,"(1x,a)") "Ensemble type : Shardlow's second order splitting (S2)"
              Else
                 Call strip_blanks(record)
                 If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                 Call error(436)
              End If

              Call get_word(record,word)
              gamdpd(0) = Abs(word_2_real(word,0.0_wp))

              If (gamdpd(0) > zero_plus) Then
                 If (idnode == 0) Write(nrite,"(1x,'drag coefficient (Dalton/ps)',6x,1p,e12.4)") gamdpd(0)
              End If

              If (lens) Call error(414)
              lens=.true.

           Else

              Call strip_blanks(record)
              If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
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

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT isotropic Langevin (Stochastic Dynamics)', &
                 & /,1x,'thermostat friction     (ps^-1)',3x,1p,e12.4,                                     &
                 & /,1x,'barostat friction       (ps^-1)',3x,1p,e12.4)") chi,tai

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

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT isotropic Berendsen', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                &
                 & /,1x,'barostat relaxation time   (ps)',3x,1p,e12.4)") taut,taup

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              keyens = 22

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT isotropic Nose-Hoover (Melchionna)', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                               &
                 & /,1x,'barostat relaxation time   (ps)',3x,1p,e12.4)") taut,taup

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'mtk') Then

              keyens = 23

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT isotropic Martyna-Tuckerman-Klein', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                              &
                 & /,1x,'barostat relaxation time   (ps)',3x,1p,e12.4)") taut,taup

              If (lens) Call error(414)
              lens=.true.

           Else

              Call strip_blanks(record)
              If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
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

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT anisotropic Langevin (Stochastic Dynamics)', &
                 & /,1x,'thermostat friction     (ps^-1)',3x,1p,e12.4,                                       &
                 & /,1x,'barostat friction       (ps^-1)',3x,1p,e12.4)") chi,tai

!                 taut=chi
!                 taup=2.0_wp*pi/tai

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 If (idnode == 0) Write(nrite,"(2(/,1x,a))")                     &
                    'semi-isotropic barostat : constant normal pressure (Pn) &', &
                    '       (N-Pn-A-T)       : constant surface area (A)'
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 If (idnode == 0) Write(nrite,"(3(/,1x,a),1p,e11.4)")             &
                    'semi-isotropic barostat : constant normal pressure (Pn) &',  &
                    '     (N-Pn-gamma-T)     : constant surface tension (gamma)', &
                    'sumulation surface tension (dyn/cm)', ten
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    iso=2
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : orthorhombic MD cell constraints'
                 Else If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (iso >= 1 .and. iso <= 2 .and. idnode == 0) Write(nrite,'(2(/,1x,a))') &
                 '*** warning - semi-isotropic ensembles are only correct for ***',     &
                 '*** infinite interfaces placed perpendicularly to the z axis !!! ***'

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'ber') Then

              keyens = 31

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT anisotropic Berendsen', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                  &
                 & /,1x,'barostat relaxation time   (ps)',3x,1p,e12.4)") taut,taup

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 If (idnode == 0) Write(nrite,"(2(/,1x,a))")                     &
                    'semi-isotropic barostat : constant normal pressure (Pn) &', &
                    '       (N-Pn-A-T)       : constant surface area (A)'
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 If (idnode == 0) Write(nrite,"(3(/,1x,a),1p,e11.4)")             &
                    'semi-isotropic barostat : constant normal pressure (Pn) &',  &
                    '     (N-Pn-gamma-T)     : constant surface tension (gamma)', &
                    'sumulation surface tension (dyn/cm)', ten
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    iso=2
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : orthorhombic MD cell constraints'
                 Else If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (iso >= 1 .and. iso <= 2 .and. idnode == 0) Write(nrite,'(2(/,1x,a))') &
                 '*** warning - semi-isotropic ensembles are only correct for ***',     &
                 '*** infinite interfaces placed perpendicularly to the z axis !!! ***'

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:6) == 'hoover') Then

              keyens = 32

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT anisotropic Nose-Hoover (Melchionna)', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                                 &
                 & /,1x,'barostat relaxation time   (ps)',3x,1p,e12.4)") taut,taup

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 If (idnode == 0) Write(nrite,"(2(/,1x,a))")                     &
                    'semi-isotropic barostat : constant normal pressure (Pn) &', &
                    '       (N-Pn-A-T)       : constant surface area (A)'
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 If (idnode == 0) Write(nrite,"(3(/,1x,a),1p,e11.4)")             &
                    'semi-isotropic barostat : constant normal pressure (Pn) &',  &
                    '     (N-Pn-gamma-T)     : constant surface tension (gamma)', &
                    'sumulation surface tension (dyn/cm)', ten
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    iso=2
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : orthorhombic MD cell constraints'
                 Else If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (iso >= 1 .and. iso <= 2 .and. idnode == 0) Write(nrite,'(2(/,1x,a))') &
                 '*** warning - semi-isotropic ensembles are only correct for ***',     &
                 '*** infinite interfaces placed perpendicularly to the z axis !!! ***'

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:3) == 'mtk') Then

              keyens = 33

              Call get_word(record,word)
              taut = Abs(word_2_real(word))
              Call get_word(record,word)
              taup = Abs(word_2_real(word))

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NPT anisotropic Martyna-Tuckerman-Klein', &
                 & /,1x,'thermostat relaxation time (ps)',3x,1p,e12.4,                                &
                 & /,1x,'barostat relaxation time   (ps)',3x,1p,e12.4)") taut,taup

              Call get_word(record,word)
              If      (word(1:4) == 'area') Then
                 iso=1
                 If (idnode == 0) Write(nrite,"(2(/,1x,a))")                     &
                    'semi-isotropic barostat : constant normal pressure (Pn) &', &
                    '       (N-Pn-A-T)       : constant surface area (A)'
              Else If (word(1:4) == 'tens') Then
                 iso=2
                 Call get_word(record,word)
                 ten = Abs(word_2_real(word))
                 If (idnode == 0) Write(nrite,"(3(/,1x,a),1p,e11.4)")             &
                    'semi-isotropic barostat : constant normal pressure (Pn) &',  &
                    '     (N-Pn-gamma-T)     : constant surface tension (gamma)', &
                    'sumulation surface tension (dyn/cm)', ten
                 ten=ten/tenunt

                 Call get_word(record,word)
                 If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else If (Len_Trim(word) > 0) Then
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (word(1:4) == 'orth') Then
                 Call get_word(record,word)
                 If (Len_Trim(word) == 0) Then
                    iso=2
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : orthorhombic MD cell constraints'
                 Else If (word(1:4) == 'semi') Then
                    iso=3
                    If (idnode == 0) Write(nrite,"(1x,a)") &
                       'semi-isotropic barostat : semi-orthorhombic MD cell constraints'
                 Else
                    Call strip_blanks(record)
                    If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                    Call warning(460,0.0_wp,0.0_wp,0.0_wp)
                 End If
              Else If (Len_Trim(word) > 0 ) Then
                 Call strip_blanks(record)
                 If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                 Call warning(460,0.0_wp,0.0_wp,0.0_wp)
              End If
              If (iso >= 1 .and. iso <= 2 .and. idnode == 0) Write(nrite,'(2(/,1x,a))') &
                 '*** warning - semi-isotropic ensembles are only correct for ***',     &
                 '*** infinite interfaces placed perpendicularly to the z axis !!! ***'

              If (lens) Call error(414)
              lens=.true.

           Else

              Call strip_blanks(record)
              If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
              Call error(436)

           End If

        Else

           Call strip_blanks(record)
           If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
           Call error(436)

        End If

! For Langevin ensembles that require arrays

        If (l_lan) Call langevin_allocate_arrays()

! read density variation option

     Else If (word(1:7) == 'densvar') Then

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,1x,'density variation allowance (%)',3x,1p,e12.4)") tmp

! read real space cutoff

     Else If (word(1:3) == 'cut' .or. word(1:4) == 'rcut') Then


        Call get_word(record,word)
        rcut1 = Abs(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,1x,'real space cutoff (Angs)    ',6x,1p,e12.4)") rcut1

! read real space cutoff padding

     Else If (word(1:3) == 'pad' .or. word(1:4) == 'rpad') Then

        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        rpad1 = Abs(word_2_real(word))
        If (idnode == 0) Write(nrite,"(/,1x,'cutoff padding (Angs)       ',6x,1p,e12.4)") rpad1

! read vdw cutoff (short-range potentials)

     Else If (word(1:4) == 'rvdw') Then

        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        rvdw1 = Abs(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,1x,'vdw cutoff (Angs) ',16x,1p,e12.4)") rvdw1

! read Ewald sum parameters

     Else If (word(1:5) == 'ewald' .or. word(1:4) == 'spme') Then

        Call get_word(record,word)

        If (word(1:5) == 'evalu') Then

! This is sorted in set_bounds -> scan_control

        Else

           keyfce = 2

           If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : Smooth Particle Mesh Ewald')")

           If (word(1:9) == 'precision') Then
              Call get_word(record,word)
              tmp = Abs(word_2_real(word))
              If (idnode == 0) Write(nrite,"(1x,'Ewald sum precision         ',6x,1p,e12.4)") tmp
           End If

! This is sorted in set_bounds -> scan_control

           If (idnode == 0) Then
              Write(nrite,"(1x,'Ewald convergence parameter (A^-1)',1p,e12.4)") alpha
              Write(nrite,"(1x,'Ewald kmax1 kmax2 kmax3   (x2)',1x,3i5)") kmaxa1,kmaxb1,kmaxc1
              If (kmaxa /= kmaxa1 .or. kmaxb /= kmaxb1 .or. kmaxc /= kmaxc1) &
              Write(nrite,"(1x,'DaFT adjusted kmax values (x2)',1x,3i5)") kmaxa,kmaxb,kmaxc
              Write(nrite,"(1x,'B-spline interpolation order',8x,1p,i5)") mxspl
           End If

! Print infrequent k-space SPME evaluation

           If      (nstfce == 0) Then
              Call warning(370,Real(nstfce,wp),1.0_wp,0.0_wp)
              nstfce=1
           Else If (nstfce > 10) Then
              Call warning(370,Real(nstfce,wp),4.0_wp,0.0_wp)
              nstfce=4
           End If
           If (nstfce >= 1 .and. idnode == 0) &
              Write(nrite,"(1x,'k-space evaluation interval (steps)',1x,1p,i5)") nstfce

           If (lforc) Call error(416)
           lforc=.true.

        End If

! read distance dependent dielectric option

     Else If (word(1:6) == 'distan') Then

        keyfce = 4
        If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : Distance Dependent Dielectric')")

        If (lforc) Call error(416)
        lforc=.true.

! read coulombic potential option

     Else If (word(1:4) == 'coul') Then

        keyfce = 6
        If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : Coulombic Potential')")

        If (lforc) Call error(416)
        lforc=.true.

! read force-shifted coulombic potential option

     Else If (word(1:5) == 'shift') Then

        keyfce = 8
        If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : Force-Shifted Coulombic Potential')")

        Call get_word(record,word)

        If      (word(1:4) == 'damp') Then
           Call get_word(record,word)
           alpha = Abs(word_2_real(word))
           If (idnode == 0) Write(nrite,"(1x,'damping parameter (A^-1)',10x,1p,e12.4)") alpha
        Else If (word(1:9) == 'precision') Then
           Call get_word(record,word)
           eps0 = Abs(word_2_real(word))
           If (idnode == 0) Write(nrite,"(1x,'precision parameter     ',10x,1p,e12.4)") eps0
           eps0 = Max(Min(eps0,0.5_wp),1.0e-20_wp)
           tol = Sqrt(Abs(Log(eps0*rcut)))
           alpha = Sqrt(Abs(Log(eps0*rcut*tol)))/rcut
           If (idnode == 0) Write(nrite,"(1x,'damping parameter (A^-1) derived',2x,1p,e12.4)") alpha
        End If
        If (alpha > zero_plus) Then
           If (idnode == 0) Write(nrite,"(1x,'Fennell damping applied')")
           If (rcut < 12.0_wp) Call warning(7,rcut,12.0_wp,0.0_wp)
        End If

        If (lforc) Call error(416)
        lforc=.true.

! read reaction field option

     Else If (word(1:8) == 'reaction') Then

        keyfce = 10
        If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : Reaction Field')")

        If (word(1:5) == 'field') Call get_word(record,word)
        Call get_word(record,word)

        If      (word(1:4) == 'damp') Then
           Call get_word(record,word)
           alpha = Abs(word_2_real(word))
           If (idnode == 0) Write(nrite,"(1x,'damping parameter (A^-1)',10x,1p,e12.4)") alpha
        Else If (word(1:9) == 'precision') Then
           Call get_word(record,word)
           eps0 = Abs(word_2_real(word))
           If (idnode == 0) Write(nrite,"(1x,'precision parameter     ',10x,1p,e12.4)") eps0
           eps0 = Max(Min(eps0,0.5_wp),1.0e-20_wp)
           tol = Sqrt(Abs(Log(eps0*rcut)))
           alpha = Sqrt(Abs(Log(eps0*rcut*tol)))/rcut
           If (idnode == 0) Write(nrite,"(1x,'damping parameter (A^-1) derived',2x,1p,e12.4)") alpha
        End If
        If (alpha > zero_plus) Then
           If (idnode == 0) Write(nrite,"(1x,'Fennell damping applied')")
           If (rcut < 12.0_wp) Call warning(7,rcut,12.0_wp,0.0_wp)
        End If

        If (lforc) Call error(416)
        lforc=.true.

     Else If (word(1:5) == 'poiss' .or. word(1:5) == 'psolv' ) Then

        keyfce = 12
        If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : Poisson equation solver')")

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

        If (idnode == 0) Then
           Write(nrite,"(1x,'gridspacing parameter (A)',9x,1p,e12.4)") prmps(1)
           Write(nrite,"(1x,'convergance epsilon      ',9x,1p,e12.4)") prmps(2)
           Write(nrite,"(1x,'max # of Psolver iterations',9x,1p,i5)") Nint(prmps(3))
           Write(nrite,"(1x,'max # of Jacobi  iterations',9x,1p,i5)") Nint(prmps(4))

           If ( Abs(prmps(1)-1.0_wp/alpha) > 1.0e-6_wp .or. Abs(prmps(2)-eps) > 1.0e-6_wp .or. &
                Nint(prmps(3)) == 0 .or. Nint(prmps(4)) == 0 ) Then
              Write(nrite,"(/,1x,a)") "*** warning - parameters reset to safe defaults occurred !!! ***"
              Write(nrite,"(1x,'gridspacing parameter (A)',9x,1p,e12.4)") 1.0_wp/alpha
              Write(nrite,"(1x,'convergance epsilon      ',9x,1p,e12.4)") eps
              Write(nrite,"(1x,'max # of Psolver iterations',9x,1p,i5)") mxitcg
              Write(nrite,"(1x,'max # of Jacobi  iterations',9x,1p,i5)") mxitjb
           End If
        End If

        If (lforc) Call error(416)
        lforc=.true.

! read relative dielectric constant

     Else If (word(1:3) == 'eps') Then

        Call get_word(record,word)
        If (word(1:8) == 'constant') Call get_word(record,word)
        epsq = word_2_real(word)
        If (idnode == 0) Write(nrite,"(/,1x,'relative dielectric constant',6x,1p,e12.4)") epsq

!     Else If (word(1:6) == 'induce') Then
!
!        If (idnode == 0) Write(nrite,"(/,1x,a)") "Employing induced dipoles"
!        If (idnode == 0) Write(nrite,"(/,1x,a)") "Induced dipole conjugate gradient :"
!        If (idnode == 0) Write(nrite,"(/,1x,a,i0)") "        max number of steps = ", politer
!        If (idnode == 0) Write(nrite,"(/,1x,a,e7.3)") "        convergence criterion = ", convcrit
!
!     Else If (word(1:4) == 'gear') Then
!
!       If (idnode == 0) Write(nrite,"(/,1x,a,io,a)") "Using gear predictor with ",numcof," points"
!
!     Else If (word(1:4) == 'aspc') Then
!
!       If (idnode == 0) Write(nrite,"(/,1x,a,io,a)") "Using always stable predictor corrector with ",numcof," points"
!
!     Else If (word(1:5) == 'lstsq') Then
!
!       If (idnode == 0) Write(nrite,"(/,1x,a,io,a)") "Using least squares predictor with ",numcof," points"

! read option for accounting for extended coulombic exclusion

     Else If (word(1:5) == 'exclu') Then

        lecx = .true.
        If (idnode == 0) Write(nrite,"(/,1x,'Extended Coulombic eXclusion opted for')")

! read force capping option

     Else If (word(1:3) == 'cap') Then

        lfcap = .true.

        Call get_word(record,word)
        If (word(1:5) == 'force') Call get_word(record,word)

        tmp = Abs(word_2_real(word))
        If (tmp > zero_plus) fmax=tmp
        If (idnode == 0) Write(nrite,"(/,1x,'force capping on (during equilibration)', &
           & /,1x,'force capping limit (kT/Angs)',5x,1p,e12.4)") fmax

! read 'no vdw', 'no elec', 'no ind' and 'no vafav' options

     Else If (word(1:2) == 'no') Then

        Call get_word(record,word1)

        If      (word1(1:3) == 'vdw' ) Then

        Else If (word1(1:4) == 'elec') Then

        Else If (word1(1:3) == 'ind' ) Then

        Else If (word1(1:3) == 'str' ) Then

           If (idnode == 0) Write(nrite,"(8(/,1x,a))") "no strict option on" ,                                             &
              "*** Warning *** It skips printing inessential information in OUTPUT such as many warnings, FIELD",          &
              "*** wArning *** digested information and full iteration cycles information from CGM based routines!",       &
              "*** waRning *** However, it also assumes some, deemed safe, defaults for some specified as well as",        &
              "*** warNing *** unspecified by the user options, that may or may not be needed for the simulation to run!", &
              "*** warnIng *** The defaults are deemed to deliver safer passage as well as optimal performance without",   &
              "*** warniNg *** sacrificing on accuracy!  While it may, by chance, help to pass previously failing runs",   &
              "*** warninG *** it may as well lead to a run failure without warnings!  Beware, avoid usage if uncertain!"

        Else If (word1(1:3) == 'top' ) Then

           If (idnode == 0) Write(nrite,"(/,1x,a)") "no topology option on (avoids printing extended FIELD topology in OUTPUT)"

           l_top = .false.

        Else If (word1(1:5) == 'vafav') Then

           lvafav = .false.

        Else If (word1(1:3) == 'vom' ) Then ! "no vom" should be used with TTM

           If (idnode == 0 .and. .not.l_ttm) Write(nrite,"(3(/,1x,a))")                    &
              '"no vom" option auto-switched on - COM momentum removal will be abandoned', &
              '*** warning - this may lead to a build up of the COM momentum and ***',     &
              '***           a manifestation of the "flying ice-cube" effect !!! ***'

           l_vom    = .false.

        Else If (word1(1:4) == 'link') Then ! NON-TRANSFERABLE OPTION FROM DL_POLY_2

           Call warning(38,0.0_wp,0.0_wp,0.0_wp)

        Else

           Call strip_blanks(record)
           If (idnode == 0) Write(nrite,"(/,/,3a)") word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call error(3)

        End If

! read tolerance for relaxed shell model

     Else If (word(1:6) == 'rlxtol') Then

        Call get_word(record,word)
        rlx_tol(1) = Max(1.0_wp,Abs(word_2_real(word)))
        If (idnode == 0) Write(nrite,"(/,1x,'relaxed shell model CGM tolerance',1x,1p,e12.4)") rlx_tol(1)

        Call get_word(record,word1)
        rlx_tol(2) = word_2_real(word1,-1.0_wp)
        If (rlx_tol(2) > zero_plus .and. idnode == 0) &
                         Write(nrite,"(/,1x,'relaxed shell model CGM step     ',1x,1p,e12.4)") rlx_tol(2)

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
          If (idnode == 0) Write(nrite,"(/,1x,'Two Temperature Model (TTM) opted for')")
        End If

        Call get_word(record,word1)

        If (word1(1:4) == 'ncit') Then

        ! number of coarse-grained ion temperature cells (CIT):
        ! already determined in scan_control

           If (idnode == 0) Write(nrite,"(/,1x,'ionic temperature grid size      (x,y,z): ',3(2x,i8),&
                                         &/,1x,'temperature cell size (A)        (x,y,z): ',3(2x,f8.4),&
                                         &/,1x,'average number of atoms/cell              ',f10.4)") &
                                           ntsys(1),ntsys(2),ntsys(3),delx,dely,delz,sysrho*volume

        Else If (word1(1:4) == 'ncet') Then

        ! number of coarse-grained electronic temperature cells (CET):
        ! already determined in scan_control

           If (idnode == 0 ) Write(nrite,"(/,1x,'electronic temperature grid size (x,y,z): ',3(2x,i8))") &
                                           eltsys(1),eltsys(2),eltsys(3)

        Else If (word1(1:5) == 'metal') Then

        ! sets properties of electronic subsystem as a metal:
        ! already determined in scan_control

          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic subsystem represents metal: thermal conductivity required')")
          End If

        Else If (word1(1:8) == 'nonmetal') Then

        ! sets properties of electronic subsystem as a non-metal

          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic subsystem represents non-metal: thermal diffusivity required')")
          End If

        Else If (word1(1:7) == 'ceconst') Then

        ! electronic specific heat capacity given as constant value

          Call get_word(record,word)
          Ce0 = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic specific heat capacity set to constant value')")
            Write(nrite,"(1x,'electronic s.h.c. (kB/atom)',7x,1p,e12.4)") Ce0
          End If

        Else If (word1(1:6) == 'cetanh') Then

        ! electronic specific heat capacity given as tanh function

          Call get_word(record,word)
          sh_A = word_2_real(word)
          Call get_word(record,word)
          sh_B = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic specific heat capacity set to hyperbolic tangent function')")
            Write(nrite,"(1x,'constant term A    (kB/atom)',6x,1p,e12.4)") sh_A
            Write(nrite,"(1x,'temperature term B    (K^-1)',6x,1p,e12.4)") sh_B
          End If

        Else If (word1(1:5) == 'celin') Then

        ! electronic specific heat capacity given as linear function
        ! up to Fermi temperature, constant afterwards

          Call get_word(record,word)
          Cemax = word_2_real(word)
          Call get_word(record,word)
          Tfermi = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic specific heat capacity set to linear function up to Fermi temperature')")
            Write(nrite,"(1x,'max. electronic s.h.c. (kB/atom)',2x,1p,e12.4)") Cemax
            Write(nrite,"(1x,'Fermi temperature            (K)',2x,1p,e12.4)") Tfermi
          End If

        Else If (word1(1:5) == 'cetab') Then

        ! electronic volumetric heat capacity given in tabulated form

          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic volumetric heat capacity given as tabulated function of temperature')")
          End If

        Else If (word1(1:5) == 'keinf') Then

        ! infinite electronic thermal conductivity

          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal conductivity set to infinity')")
          End If

        Else If (word1(1:7) == "keconst") Then

        ! electronic thermal conductivity given as constant value

          Call get_word(record,word)
          Ka0 = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal conductivity set to constant value')")
            Write(nrite,"(1x,'electronic t.c. (W m^-1 K^-1)',5x,1p,e12.4)") Ka0
          End If

        Else If (word1(1:7) == 'kedrude') Then

        ! electronic thermal conductivity given as drude model (propertional to
        ! electronic temperature, giving t.c. at system temperature)

          Call get_word(record,word)
          Ka0 = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal conductivity set to drude model')")
            Write(nrite,"(1x,'t.c. at system temp. (W m^-1 K^-1)',1p,e12.4)") Ka0
          End If

        Else If (word1(1:5) == 'ketab') Then

        ! electronic thermal conductivity given in tabulated form

          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal conductivity given as tabulated function of temperature:',&
                        & /,1x,'uses ionic or system temperature to calculate cell conductivity value',&
                        & /,1x,'for thermal diffusion equation')")
          End If

        Else If (word1(1:4) == 'diff' .or. word1(1:7)=='deconst') Then

        ! electronic thermal diffusivity given as constant value
        ! (for non-metal systems)

          Call get_word(record,word)
          Diff0 = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal diffusivity set to constant value')")
            Write(nrite,"(1x,'electronic t.d. (m^2 s^-1)',8x,1p,e12.4)") Diff0
          End If

        Else If (word1(1:7) == 'derecip') Then

        ! electronic thermal diffusivity given as reciprocal function
        ! of temperature (up to Fermi temperature), constant afterwards

          Call get_word(record,word)
          Diff0 = word_2_real(word)
          Call get_word(record,word)
          Tfermi = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal diffusivity set to reciprocal function up to Fermi temperature')")
            Write(nrite,"(1x,'datum electronic t.d. (m^2 s^-1)',2x,1p,e12.4)") Diff0
            Write(nrite,"(1x,'Fermi temperature            (K)',2x,1p,e12.4)") Tfermi
          End If

        Else If (word1(1:4) == 'detab') Then

        ! electronic thermal diffusivity given in tabulated form

          If (idnode == 0) Then
            Write(nrite,"(/,1x,'electronic thermal diffusivity given as tabulated function of temperature')")
          End If

        Else If (word1(1:8) == 'atomdens') Then

        ! user-specified atomic density, used to convert specific
        ! heat capacities to volumetric values

          Call get_word(record,word)
          cellrho = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'user-specified atomic density (A^-3)',6x,f10.4)") cellrho
          End If

        Else If (word1(1:7) == 'dyndens') Then

        ! dynamic calculation of atom density in active cells during
        ! TTM calculations, used to convert specific heat capacities
        ! to volumetric values

          ttmdyndens = .true.
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'dynamic calculations of average atomic density in active ionic cells')")
          End If

        Else If (word1(1:4) == 'amin') Then

        ! minimum number of atoms needed per ionic temperature cell
        ! to give definable ionic temperature (default = 1): smaller
        ! number deactivates ionic and electronic temperature cells
        ! (by default, electronic energies are not redistributed)

          Call get_word(record,word)
          amin = Abs(Nint(word_2_real(word)))
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'min. atom no. for ionic cells',5x,1p,i8)") amin
          End If

        Else If (word1(1:6) == 'redist') Then

        ! redistribution of electronic energy from deactivated cells 
        ! to active neighbours

          If (redistribute .and. idnode == 0) Then
            Write(nrite,"(/,1x,'redistributing energy from deactivated electronic cells into active neighbours',&
                        & /,1x,'(requires at least one electronic temperature cell beyond ionic cells)')")
          End If

        Else If (word1(1:4) == 'dedx') Then

        ! electronic stopping power of projectile entering electronic system

          Call get_word(record,word)
          dEdX = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'elec. stopping power (eV/nm)',6x,1p,e12.4)") dEdX
          End If

        Else If (word1(1:6) == 'sgauss' .or. word1(1:5) == 'sigma') Then

        ! gaussian spatial distribution for initial energy deposition into
        ! electronic system

          sdepoType = 1
          Call get_word(record,word)
          sig = word_2_real(word)
          Call get_word(record,word)
          sigmax = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'initial gaussian spatial energy deposition in electronic system')")
            Write(nrite,"(1x,'sigma of distribution (nm)',8x,1p,e12.4)") sig
            Write(nrite,"(1x,'distribution cutoff   (nm)',8x,1p,e12.4)") sigmax*sig
          End If

        Else If (word1(1:5) == 'sflat') Then

        ! homogeneous spatial distribution for initial energy deposition into
        ! electronic system

          sdepoType = 2
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'initial homogeneous (flat) spatial energy deposition in electronic system')")
          End If

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
          If (idnode == 0) Then
            Select Case (sdepoType)
            Case (2)
              Write(nrite,"(/,1x,'initial homogeneous (flat) spatial energy deposition in electronic system due to laser')")
              Write(nrite,"(1x,'absorbed fluence (mJ cm^-2)',7x,1p,e12.4)") fluence
              Write(nrite,"(1x,'penetration depth      (nm)',7x,1p,e12.4)") pdepth
            Case (3)
              Write(nrite,"(/,1x,'initial xy-homogeneous, z-exponential decaying spatial energy deposition in',/&
                            &,1x,'electronic system due to laser')")
              Write(nrite,"(1x,'absorbed fluence at surface (mJ cm^-2)',7x,1p,e12.4)") fluence
              Write(nrite,"(1x,'penetration depth                 (nm)',7x,1p,e12.4)") pdepth
            End Select
          End If

        Else If (word1(1:5) == 'gauss') Then

        ! gaussian temporal distribution for energy deposition into
        ! electronic system

          tdepoType = 1
          Call get_word(record,word)
          tdepo = word_2_real(word)
          Call get_word(record,word)
          tcdepo = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'gaussian temporal energy deposition in electronic system')")
            Write(nrite,"(1x,'sigma of distribution (ps)',8x,1p,e12.4)") tdepo
            Write(nrite,"(1x,'distribution cutoff   (ps)',8x,1p,e12.4)") 2.0_wp*tcdepo*tdepo
          End If

        Else If (word1(1:5) == 'nexp') Then

        ! decaying exponential temporal distribution for energy deposition into
        ! electronic system

          tdepoType = 2
          Call get_word(record,word)
          tdepo = word_2_real(word)
          Call get_word(record,word)
          tcdepo = word_2_real(word)
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'decaying exponential temporal energy deposition in electronic system')")
            Write(nrite,"(1x,'tau of distribution (ps)',10x,1p,e12.4)") tdepo
            Write(nrite,"(1x,'distribution cutoff (ps)',10x,1p,e12.4)") tcdepo*tdepo
          End If

        Else If (word1(1:5) == 'delta') Then

        ! dirac delta temporal distribution for energy deposition into
        ! electronic system

          tdepoType = 3
          If (idnode == 0) Then
            Write(nrite,"(/,1x,'dirac delta temporal energy deposition in electronic system')")
          End If

        Else If (word1(1:5) == 'pulse') Then

        ! square pulse temporal distribution for energy deposition into
        ! electronic system (defaults to dirac delta if pulse duration
        ! set to zero)

          tdepoType = 4
          Call get_word(record,word)
          tdepo = word_2_real(word)
          If (tdepo<=zero_plus) Then
            tdepoType = 3
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'square pulse temporal energy deposition in electronic system',/&
                             &1x,'of zero duration: being treated as dirac delta')")
            End If
          Else
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'square pulse temporal energy deposition in electronic system')")
              Write(nrite,"(1x,'pulse duration (ps)',15x,1p,e12.4)") tdepo
            End If
          End If

        Else If (word1(1:4) == 'varg') Then

        ! variable electron-phonon coupling constant (chi_ep) based on
        ! tabular electronic stopping terms (in g.dat file): option to
        ! apply value homogeneously across system (based on average 
        ! electronic temperature) or heterogeneously (using local 
        ! electronic temperature for each voxel)

          Select Case (gvar)
          Case (1)
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'variable electron-phonon coupling values to be applied homogeneously',&
                          & /,1x,'(overrides value given for ensemble, required tabulated stopping',&
                           & /,2x,'terms in g.dat file)')")
            End If
          Case (2)
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'variable electron-phonon coupling values to be applied heterogeneously',&
                          & /,1x,'(overrides value given for ensemble, required tabulated stopping',&
                           & /,2x,'terms in g.dat file)')")
            End If
          End Select

        Else If (word1(1:3) == 'bcs') Then

        ! electronic temperature boundary conditions

          Call get_word(record,word)

          If (word(1:8) == 'periodic') Then
            bcTypeE = 1
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'electronic temperature boundary conditions set as periodic')")
            End If
          Else If (word(1:6) == 'dirich') Then
            bcTypeE = 2
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'electronic temperature boundary conditions set as dirichlet:',&
                          & /,1x,'setting boundaries to system temperature')")
            End If
          Else If (word(1:7) == 'neumann') Then
            bcTypeE = 3
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'electronic temperature boundary conditions set as neumann:',&
                          & /,1x,'zero energy flux at boundaries')")
            End If
          Else If (word(1:8) == 'xydirich') Then
            bcTypeE = 4
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'electronic temperature boundary conditions set as dirichlet (xy), neumann (z):',&
                          & /,1x,'system temperature at x and y boundaries',&
                          & /,1x,'zero energy flux at z boundaries')")
            End If
          Else If (word(1:5) == 'robin') Then
            bcTypeE = 5
            Call get_word(record,word)
            fluxout = word_2_real(word)
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'electronic temperature boundary conditions set as robin:',&
                          & /,1x,'temperature leakage at boundaries of ',1p,e11.4)") fluxout
            End If
          Else If (word(1:7) == 'xyrobin') Then
            bcTypeE = 6
            Call get_word(record,word)
            fluxout = word_2_real(word)
            If (idnode == 0) Then
              Write(nrite,"(/,1x,'electronic temperature boundary conditions set as robin (xy), neumann (z):',&
                          & /,1x,'temperature leakage at x and y boundaries of ',1p,e11.4,&
                          & /,1x,'zero energy flux at z boundaries')") fluxout
            End If
          End If

        Else If (word1(1:6) == 'offset') Then

        ! time offset in coupling electronic and ionic systems

          Call get_word(record,word)
          ttmoffset = word_2_real(word)
          If (idnode == 0) Write(nrite,"(/,1x,'electron-ion coupling offset (ps)',1x,1p,e12.4)") ttmoffset

        Else If (word1(1:6) == 'oneway') Then

        ! one-way electron-phonon coupling in thermostat and thermal
        ! diffusion: only apply when electronic temperature exceeds
        ! ionic temperature

          oneway = .true.
          If (idnode == 0) Write(nrite,"(/,1x,'one-way electron-phonon coupling option switched on')")

        Else If (word1(1:5) == 'stats') Then

        ! ttm statistics (minimum and maximum ionic/electronic temperatures,
        ! electronic energy) file option and output frequency

          Call get_word(record,word)
          ttmstats = Abs(Nint(word_2_real(word)))
          If (idnode == 0) Write(nrite,"(/,1x,'ttm statistics file option on', &
             & /,1x,'ttm statistics file interval',3x,i10)") ttmstats

        Else If (word1(1:4) == 'traj') Then

        ! ttm trajectory (one-dimensional ionic and electronic 
        ! temperature profile) file option and output frequency

          Call get_word(record,word)
          ttmtraj = Abs(Nint(word_2_real(word)))
          If (idnode == 0) Write(nrite,"(/,1x,'ttm trajectory (temperature profile) file option on', &
             & /,1x,'ttm trajectory file interval',3x,i10)") ttmtraj

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
           If (idnode == 0) Write(nrite,"(/,/,3a)") word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
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
        If (word(1:6) == 'errors') Then
! read if we're doing rdf error analysis
          Call get_word(record,word)
          If(word(1:4) == 'jack') Then
            l_jack = .TRUE.
            Call get_word(record,word)
            itmp = Nint(word_2_real(word, 1.0_wp))
            If(itmp > 1) num_blocks = itmp
          Else
            l_block = .TRUE.
            Call get_word(record,word)
            itmp = Nint(word_2_real(word, 1.0_wp))
            If(itmp > 1) num_blocks = itmp
          End If
        End If

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
           If (idnode == 0) Write(nrite,"(/,1x,'data printing interval (steps)',1x,i10)") nstbpo
        End If

! read stack option (reading done in set_bounds -> scan_control)

     Else If (word(1:5) == 'stack') Then

        If (idnode == 0) Write(nrite,"(/,1x,'data stacking interval (steps)',1x,i10)") mxstak

! read statistics printing option

     Else If (word(1:4) == 'stat') Then

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        intsta = Nint(word_2_real(word))
        If (idnode == 0) Write(nrite,"(/,1x,'statistics file interval    ',3x,i10)") intsta

! read MSDTMP printing option

     Else If (word(1:6) == 'msdtmp') Then

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        nstmsd = Max(nstmsd,itmp)

        Call get_word(record,word)
        itmp = Abs(Nint(word_2_real(word)))
        istmsd = Max(istmsd,itmp)

        If (idnode == 0) Write(nrite,"(/,1x,'MSDTMP file option on', &
           & /,1x,'MSDTMP file start    ',10x,i10,                   &
           & /,1x,'MSDTMP file interval ',10x,i10)") nstmsd,istmsd

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

        If (idnode == 0) Write(nrite,"(/,1x,'trajectory file option on', &
           & /,1x,'trajectory file start       ',3x,i10,                 &
           & /,1x,'trajectory file interval    ',3x,i10,                 &
           & /,1x,'trajectory file info key    ',3x,i10)") nstraj,istraj,keytrj

        If (keytrj > 3) Call error(517)
        If (keytrj == 3 .and. idnode == 0) Write(nrite,'(2(/,1x,a))')   &
           '%%% warning - trajectory file info key == 3 generates %%%', &
           '%%% HISTORY in an unindexed and consize manner !!! %%%'

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

        If (idnode == 0) Write(nrite, "(/,1x,'defects file option on    ', &
           & /,1x,'defects file start        ',5x,i10,                     &
           & /,1x,'defects file interval     ',5x,i10,                     &
           & /,1x,'defects distance condition (Angs) ',1p,e12.4)") nsdef,isdef,rdef

! REFERENCE1 forcing

        Call get_word(record,word)
        If (word(1:5) == 'extra') Then
           l_dfx=.true.
           If (idnode == 0) Write(nrite, "(/,1x,'%%% defects1 file option on %%%')")
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

        If (idnode == 0) Write(nrite, "(/,1x,'displacements file option on', &
           & /,1x,'DISPDAT file start        ',5x,i10,                       &
           & /,1x,'DISPDAT file interval     ',5x,i10,                       &
           & /,1x,'DISPDAT distance condition (Angs)' ,1p,e12.4)") nsrsd,isrsd,rrsd

! read DL_POLY_2/Classic delr Verlet shell strip cutoff option (compatibility)
! as DL_POLY_4 real space cutoff padding option

     Else If (word(1:4) == 'delr') Then

        Call warning(35,0.0_wp,0.0_wp,0.0_wp)
        Call get_word(record,word) ; If (word(1:5) == 'width') Call get_word(record,word)
        rpad1 = 0.25_wp * Abs(word_2_real(word))
        If (idnode == 0) Write(nrite,"(1x,'cutoff padding (Angs)       ',6x,1p,e12.4)") rpad1

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
           timjob = word_2_real(word)

        Else

           Call strip_blanks(record)
           If (idnode == 0) Write(nrite,"(/,/,3a)") word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
           Call error(3)

        End If

! read close-down time allowance

     Else If (word(1:5) == 'close') Then

        Call get_word(record,word1)
        If (word1(1:4) == 'time') Then

           l_timcls=.true.

           Call get_word(record,word)
           timcls = word_2_real(word)

        Else

           Call strip_blanks(record)
           If (idnode == 0) Write(nrite,"(/,/,3a)") word(1:Len_Trim(word)+1),word1(1:Len_Trim(word1)+1),record
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
        If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
        Call error(3)

     End If

  End Do

! no finish record in CONTROL file

  If (idnode == 0) Close(Unit=nread)
  Call error(17)

! unexpected end of file

1000 Continue

  If (idnode == 0) Close(Unit=nread)
  Call error(53)

! safe termination of reading CONTROL

2000 Continue
  If (idnode == 0) Close(Unit=nread)

!!! FIXES !!!
! fix on step-dependent options

  If (nstmin  == 0) nstmin  = nsteql+1
  If (nstzero == 0) nstzero = nsteql+1
  If (nstgaus == 0) nstgaus = nsteql+1
  If (nstscal == 0) nstscal = nsteql+1

!!! REPORTS !!!
! report restart

  If (keyres == 0) Then
     If (idnode == 0) Write(nrite,"(/,1x,'clean start requested')")
  Else If (levcfg == 0) Then
     Call warning(200,0.0_wp,0.0_wp,0.0_wp)
     keyres=0
  End If

! report default ensemble if none is specified:
! inhomogeneous Langevin if two-temperature model
! is in use, NVE if not

  If (.not.lens) Then
     Call warning(130,0.0_wp,0.0_wp,0.0_wp)
     If (idnode == 0) Then
        If (l_vv) Then
           Write(nrite,"(/,1x,'Integration : Velocity Verlet')")
        Else
           Write(nrite,"(/,1x,'Integration : Leapfrog Verlet')")
        End If
        If (l_ttm) Then
          Write(nrite,"(1x,'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)', &
                    & /,1x,'e-phonon friction       (ps^-1)',3x,1p,e12.4, &
                    & /,1x,'e-stopping friction     (ps^-1)',3x,1p,e12.4, &
                    & /,1x,'e-stopping velocity   (A ps^-1)',3x,1p,e12.4)") chi_ep,chi_es,vel_es2
        Else
          Write(nrite,"(1x,'Ensemble : NVE (Microcanonical)')")
        End If
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
     If (keyens==10 .or. keyens==20 .or. keyens==30 .and. chi>zero_plus) chi_ep = chi
     If (idnode == 0) Then
       Write(nrite,"(1x,'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)', &
                 & /,1x,'e-phonon friction       (ps^-1)',3x,1p,e12.4, &
                 & /,1x,'e-stopping friction     (ps^-1)',3x,1p,e12.4, &
                 & /,1x,'e-stopping velocity   (A ps^-1)',3x,1p,e12.4)") chi_ep,chi_es,vel_es2
       If (ttmthvel) Then
         Write(nrite,"(/,1x,'applying to thermal velocities in all directions')")
       Else If (ttmthvelz) Then
         Write(nrite,"(/,1x,'applying to total velocities in x and y directions,', &
                     & /,1x,'thermal velocities in z direction')")
       Else
         Write(nrite,"(/,1x,'applying to total velocities in all directions')")
       End If
     End If
     keyens = 15

  End If

! report iteration length and tolerance condition for constraints and PMF algorithms

  If ((mxcons > 0 .or. mxpmf > 0) .and. idnode == 0) Then
     Write(nrite,"(/,1x,'iterations for shake/rattle ',3x,i10)") mxshak
     Write(nrite,"(1x,'tolerance for shake/rattle (Angs) ',1p,e12.4)") tolnce
  End If

! report electrostatics

  If (l_n_e) Then
     keyfce=0
     If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics switched off!!!')")
  Else If (keyfce == 0) Then
     If (idnode == 0) Write(nrite,"(/,1x,'Electrostatics : None Assumed')")
  End If

! report for extended coulombic exclusion if needed

  If (keyfce /= 0) Then
     If (lecx) Then
        If (idnode == 0) Write(nrite,"(/,1x,'Extended Coulombic eXclusion : YES')")
     Else
        If (idnode == 0) Write(nrite,"(/,1x,'Extended Coulombic eXclusion : NO')")
     End If
  End If

! report if rcut is reset (measures taken in scan_config -
! rcut is the maximum cutoff needed in the system)

  If (Abs(rcut-rcut1) > 1.0e-6_wp .and. idnode == 0) &
     Write(nrite,"(/,1x,'real space cutoff reset to (Angs) ',1p,e12.4)") rcut

! report if rpad is reset (measures taken in scan_config & set_bounds -
! rpad is the cutoff padding needed the conditional VNL update)

  If (Abs(rpad-rpad1) > 1.0e-6_wp .and. idnode == 0) &
     Write(nrite,"(/,1x,'cutoff padding reset to (Angs)  ',2x,1p,e12.4)") rpad

! report vdw

  If (l_n_v .and. idnode == 0) Write(nrite,"(/,1x,'vdw potential terms switched off')")

! report if rvdw is reset (measures taken in scan_config)

  If ((.not.l_n_v) .and. Abs(rvdw-rvdw1) > 1.0e-6_wp .and. idnode == 0) &
     Write(nrite,"(/,1x,'vdw cutoff reset to (Angs)',8x,1p,e12.4)") rvdw

! report timestep

  If (lvar) Then

     If (keydpd > 0) Then
        lvar=.false.
        If (idnode == 0) Then
           Write(nrite,'(/,1x,a)') "*** warning - variable timestep unavalable in DPD themostats, defaulting to:"
           Write(nrite,"(/,1x,'fixed simulation timestep (ps)   ',1x,1p,e12.4)") tstep
        End If
     Else
        If (mxdis >= 2.5_wp*mndis .and. mndis > 0.0_wp) Then
           If (idnode == 0) Then
              Write(nrite,"(/,1x,'variable simulation timestep (ps)',1x,1p,e12.4)") tstep
              Write(nrite,"(/,1x,a,2(/,1x,a,5x,1p,e12.4))") &
              "controls for variable timestep",             &
              "minimum distance Dmin (Angs) ",mndis,        &
              "maximum distance Dmax (Angs) ",mxdis
           End If

           If (mxstp > zero_plus) Then
              If (idnode == 0) Write(nrite,"(1x,a,7x,1p,e12.4)") &
              "timestep ceiling mxstp (ps)",mxstp
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

     If (idnode == 0) &
        Write(nrite,"(/,1x,'fixed simulation timestep (ps)   ',1x,1p,e12.4)") tstep

  End If

! report no vom option: its use recommended with ttm

  If (.not.l_vom .and. .not.l_ttm) Then
     If (idnode == 0) Write(nrite,"(3(/,1x,a))")                                 &
        'no vom option on - COM momentum removal will be abandoned',             &
        '*** warning - this may lead to a build up of the COM momentum and ***', &
        '***           a manifestation of the "flying ice-cube" effect !!! ***'
  Else If (l_vom .and. l_ttm) Then
     If (idnode == 0) Write(nrite,"(3(/,1x,a))")                                         &
        'no vom option off - COM momentum removal will be used',                         &
        '*** warning - this may lead to incorrect dynamic behaviour for            ***', &
        '***           two-temperature model: COM momentum removal recommended !!! ***'
  End If

! report intramolecular analysis options

  If (lpana .or. mxgana > 0) Then
     If (mxgana == 0) Then
        If (idnode == 0) Write(nrite,"(/,1x,a)") 'no intramolecular distribution collection requested'
     Else
        If (mxgbnd1 > 0 .and. mxgang1 > 0 .and. &
            mxgdih1 > 0 .and. mxginv1 > 0) Then
           If (idnode == 0) &
              Write(nrite,"(/,1x,a)") 'full intramolecular distribution collection requested (all=bnd/ang/dih/inv):'
        Else
           If (idnode == 0) Write(nrite,"(/,1x,a)") 'intramolecular distribution collection requested for:'
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
           If (idnode == 0) Write(nrite,"(1x,2(a,i10),a,f7.2,a)") &
              'bonds      - collection every ',nstbnd,' step(s); ngrid = ',mxgbnd1,' points; cutoff = ',rcbnd, ' Angs'
           If (i+j+k > 1 .and. idnode == 0) Write(nrite,"(1x,3(a,i10))") &
              'bonds      - reset values at  ',     i,'                  ',     j,'                 ',    k
        End If

        If (mxgang1 > 0) Then
           If (nstang == 0 .or. (nstang > nstana .and. nstana > 0)) Then
              nstang = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdang /= mxgang1)
           If (idnode == 0) Write(nrite,"(1x,2(a,i10),a)") &
              'angles     - collection every ',nstang,' step(s); ngrid = ',mxgang1,' points'
           If (i+j > 1 .and. idnode == 0) Write(nrite,"(1x,2(a,i10))") &
              'angles     - reset values at  ',     i,'                  ',     j
        End If

        If (mxgdih1 > 0) Then
           If (nstdih == 0 .or. (nstdih > nstana .and. nstana > 0)) Then
              nstdih = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grddih /= mxgdih1)
           If (idnode == 0) Write(nrite,"(1x,2(a,i10),a)") &
              'dihedrals  - collection every ',nstdih,' step(s); ngrid = ',mxgdih1,' points'
           If (i+j > 1 .and. idnode == 0) Write(nrite,"(1x,2(a,i10))") &
              'dihedrals  - reset values at  ',     i,'                  ',     j
        End If

        If (mxginv1 > 0) Then
           If (nstinv == 0 .or. (nstinv > nstana .and. nstana > 0)) Then
              nstinv = Merge(nstana , nstall , nstana > 0)
              i = 1
           Else
              i = 0
           End If
           j=Merge(1, 0, grdinv /= mxginv1)
           If (idnode == 0) Write(nrite,"(1x,2(a,i10),a)") &
              'inversions - collection every ',nstinv,' step(s); ngrid = ',mxginv1,' points'
           If (i+j > 1 .and. idnode == 0) Write(nrite,"(1x,2(a,i10))") &
              'inversions - reset values at  ',     i,'                  ',     j
        End If
     End If

     If (lpana) Then
        If (idnode == 0) Write(nrite,"(/,1x,a)") 'probability distribution analysis printing requested'
     Else
        If (idnode == 0) Write(nrite,"(/,1x,a)") 'no probability distribution analysis printing requested'
     End If
  End If

! For safety make them /= 0

  nstbnd=Max(1,nstbnd)
  nstang=Max(1,nstang)
  nstdih=Max(1,nstdih)
  nstinv=Max(1,nstinv)

! report rdf
  l_errors_block = l_errors_block .and. lrdf
  l_errors_jack = l_errors_jack .and. lrdf

  If (lrdf .or. lprdf) Then
     If (lrdf) Then
        If (idnode == 0) Then
           Write(nrite,"(/,1x,'rdf collection requested')")
           Write(nrite,"(  1x,'rdf collection interval',8x,i10)") nstrdf
           Write(nrite,"(  1x,'rdf binsize (Angstroms)',11x,1p,e12.4)") rbin
        End If
     Else
        If (idnode == 0) Write(nrite,"(/,1x,'no rdf collection requested')")
     End If

     If (lprdf) Then
        If (idnode == 0) Write(nrite,"(1x,'rdf printing requested')")
     Else
        If (lpana) Then
           If (idnode == 0) Write(nrite,"(1x,'rdf printing triggered due to a PDA printing request')")
           lprdf=lpana
        Else
           If (idnode == 0) Write(nrite,"(1x,'no rdf printing requested')")
        End If
     End If

     If (mxrdf == 0) Then
        If (idnode == 0) Write(nrite,"(1x,'no rdf pairs specified in FIELD')")
     Else
        If (idnode == 0) Write(nrite,"(1x,'rdf pairs specified in FIELD')")
     End If

     If ((.not.lrdf) .or. mxrdf == 0) Then
        If (idnode == 0) Write(nrite,"(1x,'rdf routines not to be activated')")
        lrdf=.false.
        lprdf=.false.
     End If
  End If

! report zden

  If (lzdn .or. lpzdn) Then
     If (lzdn) Then
        If (idnode == 0) Then
           Write(nrite,"(/,1x,'z-density profiles requested')")
           Write(nrite,"(  1x,'z-density collection interval',2x,i10)") nstzdn
           Write(nrite,"(  1x,'z-density binsize (Angstroms)',5x,1p,e12.4)") rbin
        End If
     Else
        If (idnode == 0) Write(nrite,"(/,1x,'no z-density profiles requested')")
     End If

     If (lpzdn) Then
        If (idnode == 0) Write(nrite,"(1x,'z-density printing requested')")
     Else
        If (idnode == 0) Write(nrite,"(1x,'no z-density printing requested')")
     End If

     If (.not.lzdn) Then
        If (idnode == 0) Write(nrite,"(  1x,'z-density routines not to be activated')")
        lpzdn=.false.
     End If
  End If

! report vaf

  If (vafsamp > 0 .or. lpvaf) Then
     If (vafsamp > 0) Then
        If (idnode == 0) Then
           Write(nrite,"(/,1x,'vaf profiles requested')")
           Write(nrite,"(  1x,'vaf collection frequency',7x,i10)") isvaf
           Write(nrite,"(  1x,'vaf collection binsize  ',7x,i10)") nsvaf
        End If
     Else
        If (idnode == 0) Write(nrite,"(/,1x,'no vaf collection requested')")
     End If

     If (lpvaf) Then
        If (idnode == 0) Write(nrite,"(1x,'vaf printing requested')")
     Else
        If (idnode == 0) Write(nrite,"(1x,'no vaf printing requested')")
     End If

     If (lvafav) Then
        If (idnode == 0) Write(nrite,"(1x,'time-averaged vaf profile')")
     Else
        If (idnode == 0) Write(nrite,"(1x,'instantaneous vaf profiles')")
     End If
  End If

! report thermal conductivity

!  If (ltcond) Then
!    If (idnode == 0) Then
!      Write(nrite,"(/,1x,'thermal conductivities requested')")
!      Write(nrite,"(  1x,'heat current collection binsize',8x,i10)") nsttcond
!    End If
!  End If

! report data dumping interval, subcelling threshold density and job times

  If (idnode == 0) Then
     Write(nrite,"(/,1x,'data dumping interval (steps)',2x,i10)") ndump
     Write(nrite,"(/,1x,'subcelling threshold density',6x,1p,e12.4)") pdplnc
     Write(nrite,"(/,1x,'allocated job run time   (s)',6x,1p,e12.4)") timjob
     Write(nrite,"(/,1x,'allocated job close time (s)',6x,1p,e12.4)") timcls
  End If

! report replay history

  If (.not.lsim) Then
     If (lfce) Then
        If (idnode == 0) Then
           Write(nrite,"(/,1x,'*** HISTORF will be replayed with full force recalculation ***', &
                       & /,1x,'*** There is no actual dynamics/integration!!! ***')")
        End If
     Else
        If (idnode == 0) Then
           Write(nrite,"(/,1x,'*** HISTORY will be replayed (no actual simulation) ***', &
                       & /,1x,'*** with structural properties will be recalculated ***')")
        End If

! abort if there's no structural property to recalculate

        If (.not.(lrdf .or. lzdn .or. ldef .or. l_msd .or. lrsd .or. (mxgana > 0))) Call error(580)
     End If

     If (keyres /= 0) Then
        keyres=0 ! Force clean restart
        If (idnode == 0) Write(nrite,"(/,1x,'clean start enforced')")
     End If
  End If

!!! RESORT TO DEFAULTS IF NEED BE !!!

  If      (nstrun == 0) Then !!! DRY RUN

     ltemp = .true. ! zero is ok
     lpres = .true. ! zero is ok

     lstep = .true. ! zero is not ok
     If (tstep <= zero_plus) Then
        tstep = 0.001_wp
        If (idnode == 0) &
           Write(nrite,"(/,1x,'default simulation timestep (ps) ',1x,1p,e12.4)") tstep
     End If

  Else If (.not.l_str) Then !!! NO STRICT

     If (.not.ltemp) Then ! Simulation temperature
        ltemp=.true.
        temp=300.0_wp
        If (idnode == 0) &
           Write(nrite,"(/,1x,'default simulation temperature (K)',1p,e12.4)") temp
     End If

     If (.not.lpres) Then ! Simulation pressure
        lpres=.true.
        press=0.0_wp
        If (idnode == 0) &
           Write(nrite,"(/,1x,'default simulation pressure (katms)',1p,e11.4)") press*prsunt
     End If

     If (.not.lstep) Then ! Simulation timestep
        lstep = .true.
        tstep = 0.001_wp
        If (idnode == 0) &
           Write(nrite,"(/,1x,'default simulation timestep (ps) ',1x,1p,e12.4)") tstep
     End If

     If ((.not.l_timjob) .and. (.not.l_timcls)) Then ! Job times
        timjob=1.0e8_wp
        timcls=1.0e7_wp

        If (idnode == 0) Then
           Write(nrite,"(/,1x,'allocated job run time   (s)',6x,1p,e12.4)") timjob
           Write(nrite,"(/,1x,'allocated job close time (s)',6x,1p,e12.4)") timcls
        End If
     Else If ((.not.l_timjob) .and. l_timcls) Then
        timjob=100.0_wp*timcls
        If (idnode == 0) &
           Write(nrite,"(/,1x,'allocated job run time   (s)',6x,1p,e12.4)") timjob
     Else If (l_timjob .and. (.not.l_timcls)) Then
        timcls=0.01_wp*timjob
        If (idnode == 0) &
           Write(nrite,"(/,1x,'allocated job close time (s)',6x,1p,e12.4)") timcls
     End If

  End If

  If (l_0 .and. (.not. ltemp)) Then ! zero K over zero fire
     temp = 10.0_wp
     If (idnode == 0) &
        Write(nrite,"(/,1x,'default simulation temperature (K)',1p,e12.4)") temp
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
    If (Abs(chi_es) <= zero_plus .and. idnode == 0) &
      Write(nrite,"(1x,'assuming no electronic stopping in inhomogeneous Langevin thermostat')")
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

              If (idnode == 0) Then
                 Write(nrite,"(/,1x,'tensorial system pressure specified for an npt ensemble simulation')")
                 Write(nrite,"(  1x,'scalar pressure derived from pressure tensor as p = Trace[P]/3')")
                 Write(nrite,"(  1x,'tensorial pressure to be zeroed (discarded)')")
                 Write(nrite,"(/,1x,'simulation pressure (katms)  ',5x,1p,e12.4)") press*prsunt
              End If
           Else
              Call error(387)
           End If
        Else
           If (lstrext) Then
              strext=0.0_wp

              If (idnode == 0) Then
                 Write(nrite,"(/,1x,'both tensorial and scalar system pressure specified for an npt ensemble simulation')")
                 Write(nrite,"(  1x,'tensorial pressure directive is ignored')")
                 Write(nrite,"(  1x,'tensorial pressure to be zeroed (discarded)')")
              End If
           End If
        End If
     Else If (keyens >= 30 .and. keyens <= 33) Then
        If (.not.lstrext) Then
           If (.not.lpres) Call error(387)
        Else
           If (lpres) Then
              If (idnode == 0) Then
                 Write(nrite,"(/,1x,'both tensorial and scalar system pressure specified for an nst ensemble simulation')")
                 Write(nrite,"(  1x,'scalar pressure directive is ignored')")
              End If

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
    If (Abs(sh_A) <= zero_plus .or. Abs(sh_B) <= zero_plus) Call error(681)
    sh_A = sh_A*cellrho
    sh_B = sh_B*1.0e-4_wp
  Case (2)
  ! linear electronic specific heat to Fermi temperature: maximum
  ! value will be converted from kB/atom to kB/A^3
    If (Abs(Tfermi) <= zero_plus .or. Abs(Cemax) <= zero_plus) Call error(681)
    Cemax = Cemax*cellrho
  Case (4)
  ! constant electronic specific heat: will convert from kB/atom to kB/A^3
  ! by multiplication of atomic density
  Case (5)
  ! hyperbolic tangent electronic specific heat: multiplier will be converted
  ! from kB/atom to kB/A^3, temperature term (K^-1) is now scaled by 10^-4
    If (Abs(sh_A) <= zero_plus .or. Abs(sh_B) <= zero_plus) Call error(681)
    sh_B = sh_B*1.0e-4_wp
  Case (6)
  ! linear electronic specific heat to Fermi temperature: maximum
  ! value will be converted from kB/atom to kB/A^3
    If (Abs(Tfermi) <= zero_plus .or. Abs(Cemax) <= zero_plus) Call error(681)
  End Select

  Select Case (KeType)
  ! constant and Drude thermal conductivity: convert from W m^-1 K^-1
  ! to kB ps^-1 A^-1
  Case (1,2)
    If (isMetal .and. Abs(Ka0) <= zero_plus) Call error(682)
    Ka0 = Ka0*JKms_to_kBAps
  End Select

  Select Case (DeType)
  Case (1)
  ! constant thermal diffusivity: convert from m^2 s^-1 to A^2 ps^-1
    If (.not. isMetal .and. Abs(Diff0) <= zero_plus) Call error(683)
    Diff0 = Diff0*1.0e8_wp
  Case (2)
  ! reciprocal thermal diffusivity: convert from m^2 s^-1 to A^2 ps^-1
  ! and Diff0 scaled with system temperature
    If (.not. isMetal .and. Abs(Diff0) <= zero_plus .or. Abs(Tfermi) <= zero_plus) Call error(683)
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
