Subroutine read_control                                &
           (levcfg,l_vv,l_str,l_n_e,l_n_r,l_n_v,       &
           rcut,rvdw,rbin,nstfce,alpha,width,          &
           l_exp,lecx,lfcap,l_top,lzero,lmin,          &
           ltgaus,ltscal,lvar,leql,lpse,               &
           lsim,lrdf,lprdf,lzdn,lpzdn,                 &
           ltraj,ldef,lrsd,                            &
           nx,ny,nz,imd,tmd,emd,vmx,vmy,vmz,           &
           temp,press,strext,keyres,                   &
           tstep,mndis,mxdis,mxstp,nstrun,nsteql,      &
           keymin,nstmin,min_tol,nstgaus,nstscal,      &
           keyens,iso,taut,chi,soft,gama,taup,tai,ten, &
           keypse,wthpse,tmppse,                       &
           fmax,nstbpo,intsta,keyfce,epsq,             &
           rlx_tol,mxshak,tolnce,mxquat,quattol,       &
           nstrdf,nstzdn,                              &
           nstmsd,istmsd,nstraj,istraj,keytrj,         &
           nsdef,isdef,rdef,nsrsd,isrsd,rrsd,          &
           ndump,timjob,timcls)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the simulation control parameters
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,    Only : idnode
  Use setup_module
  Use config_module,   Only : sysname
  Use langevin_module, Only : l_lan,l_gst,langevin_allocate_arrays
  Use parse_module
  Use vdw_module,      Only : ld_vdw,ls_vdw
  Use metal_module,    Only : ld_met,ls_met,tabmet
  Use defects1_module, Only : l_dfx

  Use development_module

  Implicit None

  Logical,                Intent( In    ) :: l_vv,l_str,l_n_e,l_n_r,l_n_v
  Integer,                Intent( In    ) :: levcfg
  Integer,                Intent( InOut ) :: nstfce
  Real( Kind = wp ),      Intent( In    ) :: rcut,rvdw,rbin,width
  Real( Kind = wp ),      Intent( InOut ) :: alpha

  Logical,                Intent(   Out ) :: l_exp,lecx,             &
                                             lfcap,l_top,lzero,lmin, &
                                             ltgaus,ltscal,          &
                                             lvar,leql,lpse,lsim,    &
                                             lrdf,lprdf,lzdn,lpzdn,  &
                                             ltraj,ldef,lrsd


  Integer,                Intent(   Out ) :: nx,ny,nz,imd,tmd,       &
                                             keyres,nstrun,nsteql,   &
                                             keymin,nstmin,          &
                                             nstgaus,nstscal,        &
                                             keyens,iso,             &
                                             keypse,nstbpo,          &
                                             intsta,keyfce,          &
                                             mxshak,mxquat,          &
                                             nstrdf,nstzdn,          &
                                             nstmsd,istmsd,          &
                                             nstraj,istraj,keytrj,   &
                                             nsdef,isdef,            &
                                             nsrsd,isrsd,            &
                                             ndump

  Real( Kind = wp ),      Intent(   Out ) :: emd,vmx,vmy,vmz,         &
                                             temp,press,strext(1:9),  &
                                             tstep,mndis,mxdis,mxstp, &
                                             taut,chi,soft,gama,      &
                                             taup,tai,ten,            &
                                             wthpse,tmppse,min_tol,   &
                                             fmax,epsq,rlx_tol,       &
                                             tolnce,quattol,          &
                                             rdef,rrsd,timjob,timcls


  Logical                                 :: limp,lens,lforc,        &
                                             lpres,lstrext,          &
                                             lstep,ltemp,safe,       &
                                             l_timjob,l_timcls

  Character( Len = 200 )                  :: record
  Character( Len = 40  )                  :: word,word1,word2

  Integer                                 :: i,j,itmp

  Real( Kind = wp )                       :: rcell(1:9),rcut1,rvdw1,tmp


! initialise system control variables and their logical switches

! default expansion option

  l_exp = .false.
  nx    = 1
  ny    = 1
  nz    = 1

! defaults for direct evaluation and force-shifting of VDW interactions
!
! ld_vdw = .false. ! (initialised in vdw_module)
! ls_vdw = .false. ! (initialised in vdw_module)
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
! (default 0 where 0 - Langevin+direct, 1 - Langevin and 2 - direct) and
! minimum width of the thermostatted boundaries in Angs
! minimum temperature of the thermostat

  lpse   = .false.
  keypse = 0
  wthpse = 2.0_wp
  tmppse = 1.0_wp

! default switch for conjugate gradient minimisation during equilibration

  lmin   = .false.
  keymin = -1
  nstmin = 0
  min_tol= 0.0_wp

! default switch for regaussing temperature and default number of
! steps when to be applied

  ltgaus  = .false.
  nstgaus = 0

! default switch for temperature scaling and default number of
! steps when to be applied

  ltscal  = .false.
  nstscal = 0

! default switch for zero temperature optimisation

  lzero = .false.

! default ensemble switch (not defined) and key

  lens   = .false.
  keyens = 0

! default thermostat and barostat friction time constants

  taut = 0.0_wp ! thermostat relaxation time
  chi  = 0.0_wp ! Stochastic Dynamics (SD Langevin) thermostat friction
  soft = 0.0_wp ! Softness for Andersen thermostat
  gama = 0.0_wp ! Stochastic (Langevin) friction on a thermostat
  taup = 0.0_wp ! barostat relaxation time
  tai  = 0.0_wp ! Stochastic Dynamics (SD Langevin) barostat friction
  iso  = 0      ! no semi-isotropic feature
  ten  = 0.0_wp ! surface tension

! default value for accounting extended coulombic exclusion
! (not accounted)

  lecx = .false.

! default switch for force capping and cap value

  lfcap = .false.
  fmax  = 1000.0_wp

! default switch for printing topology

  l_top = .true.

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

! Default relaxed shell model tolerance

  rlx_tol = 1.0_wp

! proceed normal simulation = don't replay history

  lsim = .true.

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

! default times for job execution and output dump

  timjob = 0.0_wp ; l_timjob=.false.
  timcls = 0.0_wp ; l_timcls=.false.

! cutoff and vdw cutoff defaults

  rcut1 = 0.0_wp
  rvdw1 = 0.0_wp

! open the simulation control file

  If (idnode == 0) Open(Unit=nread, File = 'CONTROL', Status = 'old')

! read simulation control name

  Call get_line(safe,nread,sysname)
  If (.not.safe) Go To 1000
  Call strip_blanks(sysname)

  If (.not.safe) Go To 1000

  If (idnode == 0) Write(nrite,"(/,3(1x,130('*'),/),1x,     &
     & 24('*'),5x,a72,5x,24('*'),/,3(1x,130('*'),/),/,/,1x, &
     & 'SIMULATION CONTROL PARAMETERS')") sysname

! read and process directives from CONTROL file

  Do While (.true.)

     Call get_line(safe,nread,record)
     If (.not.safe) Go To 1000
     Call lower_case(record)
     Call get_word(record,word)

! record is commented out

     If (word(1:1) == '#' .or. word(1:1) == ' ') Then

! read DEVELOPMENT options

     Else If (word(1:5) == 'l_scr') Then

        l_scr = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% OUTPUT redirected to the default output (screen) !!! %%%"

     Else If (word(1:5) == 'l_eng') Then

        l_eng = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% OUTPUT contains an extra last line with E_tot !!! %%%"

     Else If (word(1:6) == 'l_rout') Then

        l_rout = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% REVIVE writing in ASCII opted !!! %%%"

     Else If (word(1:5) == 'l_rin') Then

        l_rin = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% REVOLD reading in ASCII opted !!! %%%"

     Else If (word(1:5) == 'l_his') Then

        l_his = .true.
        l_trm = .true.
        If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% generate HISTORY after reading input & terminate !!! %%%"

     Else If (word(1:6) == 'l_scl') Then

        If (idnode == 0) Write(nrite,"(2(/,1x,a))")                                 &
           "%%% rescale CONFIG to CFGSCL, after reading input & terminate !!! %%%", &
           "%%% new cell vectors to rescale to: (read in a CONFIG-like manner) %%%"

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
           Write(nrite,"(1x,a4)") '%%% '
           Write(nrite,"(1x,a4,3f20.10)") '%%% ',cels(1:3)
           Write(nrite,"(1x,a4,3f20.10)") '%%% ',cels(4:6)
           Write(nrite,"(1x,a4,3f20.10)") '%%% ',cels(7:9)
           Write(nrite,"(1x,a4)") '%%% '
           Write(nrite,"(1x,'%%% CFGSCL volume ',2x,1p,g22.12)") tmp
        End If

        If (tmp > zero_plus) Then
           l_scl = .true.
           l_trm  = .true.
        Else
           If (idnode == 0) Write(nrite,"(/,1x,a)") "%%% OPTION ABORTED DUE TO ZERO VOLUME !!! %%%"
           l_trm  = .true.
        End If

     Else If (word(1:5) == 'l_tim') Then

        l_tim = .true.
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
           "%%% separation criterion (Angsroms) %%% ", r_dis

! read VDW options

     Else If (word(1:3) == 'vdw') Then

        Call get_word(record,word1)

        If      (word1(1:6) == 'direct') Then

! direct evaluation option

           ld_vdw = .true.
           If (idnode == 0) Write(nrite,"(/,1x,a)") "vdw direct option on"

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

        If (idnode == 0) Write(nrite,"(/,1x,'radomisation seeds supplied', &
           & /,1x,'particle (index)',15x,i10,                              &
           & /,1x,'(seed1,seed2)',15x,'(',i5,',',i5,')')") seed

! read temperature

     Else If (word(1:4) == 'temp') Then

        ltemp = .true.
        Call get_word(record,word)
        If (.not.lzero) Then
           temp = Abs(word_2_real(word))
           If (idnode == 0) Write(nrite,"(/,1x,'simulation temperature (K)  ',6x,1p,e12.4)") temp
        End If

! read zero temperature optimisation

     Else If (word(1:4) == 'zero') Then

        lzero  = .true.
        ltemp  = .true.
        ltscal = .false.
        temp   = 10.0_wp

        If (idnode == 0) Write(nrite,"(/,1x,'zero K optimisation requested', &
           & /,1x,'actual temperature reset to 10 Kelvin')")

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
        Else If (word(1:6) == 'direct') Then
           keypse = 2
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
                      Write(nrite,'(1x,a)') "thermostat control: direct temperature scaling"
              End If
              Write(nrite,"(1x,'thermostat thickness (Ang)',8x,1p,e12.4)") tmp
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
           nstmin = Abs(Nint(word_2_real(word)))
        End If

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))

        itmp=0
        If      (keymin == 0) Then
           If (tmp < 1.0_wp .or. tmp > 1000.0_wp) Then
              min_tol=50.0_wp
              itmp=1
           Else
              min_tol=tmp
           End If
        Else If (keymin == 1) Then
           If (tmp < zero_plus .or. tmp > 0.01_wp) Then
              min_tol=0.005_wp
              itmp=1
           Else
              min_tol=tmp
           End If
        Else If (keymin == 2) Then
           If (tmp < 1.0e-6_wp .or. tmp > 0.1_wp) Then
              min_tol=0.005_wp
              itmp=1
           Else
              min_tol=tmp
           End If
        End If

        If (itmp == 1) Call warning(360,tmp,min_tol,0.0_wp)

        If (word2(1:5) == 'minim') Then
           If (idnode == 0) Write(nrite,                                &
              & "(/,1x,'minimisation option on (during equilibration)', &
              &   /,1x,'minimisation criterion        ',1x,a8,          &
              &   /,1x,'minimisation frequency (steps)',1x,i10,         &
              &   /,1x,'minimisation tolerance        ',4x,1p,e12.4)")  &
              word1(1:8),nstmin,min_tol
        Else
           If (idnode == 0) Write(nrite,                               &
              & "(/,1x,'optimisation at start',                        &
              &   /,1x,'optimisation criterion        ',1x,a8,         &
              &   /,1x,'optimisation tolerance        ',4x,1p,e12.4)") &
              word1(1:8),min_tol
        End If

! read regauss option

     Else If (word(1:6) == 'regaus') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        nstgaus = Abs(Nint(word_2_real(word)))

        ltgaus =.true.
        If (idnode == 0) Write(nrite,"(/,1x,'regauss temperature on (during equilibration)', &
           & /,1x,'temperature regaussing interval',i10)") nstgaus

! read temperature scaling option

     Else If (word(1:5) == 'scale') Then

        Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        If (word(1:5) == 'every' .or. word(1:4) == 'temp') Call get_word(record,word)
        nstscal = Abs(Nint(word_2_real(word)))

        ltscal =.true.
        If (idnode == 0) Write(nrite,"(/,1x,'temperature scaling on (during equilibration)', &
           & /,1x,'temperature scaling interval',3x,i10)") nstscal

! read integration flavour

     Else If (word(1:8) == 'integrat') Then

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

           If (word(1:5) == 'evans') Then

              keyens = 1

              If (idnode == 0) Write(nrite,"(1x,'Ensemble : NVT Evans (Isokinetic)', &
                 & /,1x,'Gaussian temperature constraints in use')")

              If (lens) Call error(414)
              lens=.true.

           Else If (word(1:4) == 'lang') Then

              keyens = 10

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

              Call langevin_allocate_arrays()

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

              Call langevin_allocate_arrays()

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

! read density variation option

     Else If (word(1:7) == 'densvar') Then

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,1x,'density variation allowance (%)',3x,1p,e12.4)") tmp


! read real space cutoff

     Else If (word(1:3) == 'cut' .or. word(1:4) == 'rcut') Then


        Call get_word(record,word)
        rcut1 = Abs(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,1x,'real space cutoff (Ang)     ',6x,1p,e12.4)") rcut1

! read vdw cutoff (short-range potentials)

     Else If (word(1:4) == 'rvdw') Then

        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        rvdw1 = Abs(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,1x,'vdw cutoff (Ang)',18x,1p,e12.4)") rvdw1

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

        If (word(1:4) == 'damp') Then
           Call get_word(record,word)
           alpha = Abs(word_2_real(word))

           If (idnode == 0) Write(nrite,"(1x,'Damping parameter (A^-1)',10x,1p,e12.4)") alpha
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

        If (word(1:4) == 'damp') Then
           Call get_word(record,word)
           alpha = Abs(word_2_real(word))

           If (idnode == 0) Write(nrite,"(1x,'Damping parameter (A^-1)',10x,1p,e12.4)") alpha
        End If
        If (alpha > zero_plus) Then
           If (idnode == 0) Write(nrite,"(1x,'Fennell damping applied')")
           If (rcut < 12.0_wp) Call warning(7,rcut,12.0_wp,0.0_wp)
        End If

        If (lforc) Call error(416)
        lforc=.true.

! read relative dielectric constant

     Else If (word(1:3) == 'eps') Then

        Call get_word(record,word)
        If (word(1:8) == 'constant') Call get_word(record,word)
        epsq = word_2_real(word)
        If (idnode == 0) Write(nrite,"(/,1x,'relative dielectric constant',6x,1p,e12.4)") epsq

! read option for accounting for extended coulombic exclusion

     Else If (word(1:5) == 'exclu') Then

        lecx = .true.

! read force capping option

     Else If (word(1:3) == 'cap') Then

        lfcap = .true.

        Call get_word(record,word)
        If (word(1:5) == 'force') Call get_word(record,word)

        tmp = Abs(word_2_real(word))
        If (tmp > zero_plus) fmax=tmp
        If (idnode == 0) Write(nrite,"(/,1x,'force capping on (during equilibration)', &
           & /,1x,'force capping limit (kT/Ang)',6x,1p,e12.4)") fmax

! read 'no vdw', 'no elec' and 'no ind' options

     Else If (word(1:2) == 'no') Then

        Call get_word(record,word1)

        If      (word1(1:3) == 'vdw' ) Then

        Else If (word1(1:4) == 'elec') Then

        Else If (word1(1:3) == 'ind' ) Then

        Else If (word1(1:3) == 'str' ) Then

           If (idnode == 0) Write(nrite,"(/,1x,a)") "no strict option on (avoids printing unessential warnings in OUTPUT)"

        Else If (word1(1:3) == 'top' ) Then

           If (idnode == 0) Write(nrite,"(/,1x,a)") "no topology option on (avoids printing extended FIELD topology in OUTPUT)"

           l_top = .false.

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
        rlx_tol = Max(1.0_wp,Abs(word_2_real(word)))
        If (idnode == 0) Write(nrite,"(/,1x,'tolerance for relaxed shell model',1x,1p,e12.4)") rlx_tol

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

! read replay history option

     Else If (word(1:6) == 'replay') Then

        lsim = .false.

! read binsize option

     Else If (word(1:7) == 'binsize') Then

        Call get_word(record,word)
        tmp = Abs(word_2_real(word))
        If (Abs(rbin-tmp) > 1.0e-6_wp) Call warning(340,tmp,rcut/4.0_wp,rbin)

! read rdf calculation option

     Else If (word(1:3) == 'rdf') Then

        lrdf = .true.

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        nstrdf = Max(Abs(Nint(word_2_real(word))),1)

! read z-density profile option

     Else If (word(1:4) == 'zden') Then

        lzdn = .true.

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        nstzdn = Max(Abs(Nint(word_2_real(word))),1)

! read print options

     Else If (word(1:5) == 'print') Then

        Call get_word(record,word)

        If      (word(1:3) == 'rdf') Then
           lprdf = .true.
        Else If (word(1:4) == 'zden') Then
           lpzdn = .true.
        Else
           If (word(1:5) == 'every') Call get_word(record,word)
           nstbpo = Nint(word_2_real(word))
           nstbpo = Max(nstbpo,1)
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
        If (keytrj == 3) Write(nrite,'(2(/,1x,a))')                     &
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
           & /,1x,'defects distance condition (Ang)',2x,1p,e12.4)") nsdef,isdef,rdef

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
           & /,1x,'DISPDAT distance condition (Ang)',2x,1p,e12.4)") nsrsd,isrsd,rrsd

! read DL_POLY_2 multiple timestep option (compatibility)
! as DL_POLY_4 infrequent k-space SPME evaluation option

     Else If (word(1:4) == 'mult') Then

        Call warning(36,0.0_wp,0.0_wp,0.0_wp)

!!!! OTHER NON-TRANSFERABLE OPTIONS FROM DL_POLY_2 !!!!
! read primary cutoff option for multiple timestepping

     Else If (word(1:4) == 'prim') Then

        Call warning(34,0.0_wp,0.0_wp,0.0_wp)

! read delr Verlet shell strip cutoff

     Else If (word(1:4) == 'delr') Then

        Call warning(35,0.0_wp,0.0_wp,0.0_wp)

! read all pairs option

     Else If (word(1:3) == 'all') Then

        Call warning(37,0.0_wp,0.0_wp,0.0_wp)

!!!! OTHER NON-TRANSFERABLE OPTIONS FROM DL_POLY_2 !!!!

! read data dumping interval

     Else If (word(1:4) == 'dump') Then

        Call get_word(record,word)
        If (word(1:4) == 'data' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:4) == 'data' .or. word(1:5) == 'every') Call get_word(record,word)
        ndump = Max(Abs(Nint(word_2_real(word))),1)

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

  If (nstscal == 0) nstscal = nsteql+1
  If (nstgaus == 0) nstgaus = nsteql+1

!!! REPORTS !!!
! report restart

  If (keyres == 0) Then
     If (idnode == 0) Write(nrite,"(/,1x,'clean start requested')")
  Else If (levcfg == 0) Then
     Call warning(200,0.0_wp,0.0_wp,0.0_wp)
     keyres=0
  End If

! report default ensemble if none is specified

  If (.not.lens) Then
     Call warning(130,0.0_wp,0.0_wp,0.0_wp)
     If (idnode == 0) Then
        If (l_vv) Then
           Write(nrite,"(/,1x,'Integration : Velocity Verlet')")
        Else
           Write(nrite,"(/,1x,'Integration : Leapfrog Verlet')")
        End If
        Write(nrite,"(1x,'Ensemble : NVE (Microcanonical)')")
     End If
     lens=.true.
  End If

! report iteration length and tolerance condition for constraints and PMF algorithms

  If ((mxcons > 0 .or. mxpmf > 0) .and. idnode == 0) Then
     Write(nrite,"(/,1x,'iterations for shake/rattle ',3x,i10)") mxshak
     Write(nrite,"(1x,'tolerance for shake/rattle (Ang)',2x,1p,e12.4)") tolnce
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

! report if rcut reset (measures taken in scan_config -
! rcut is the maximum cutoff needed in the system)

  If (Abs(rcut-rcut1) > 1.0e-6_wp .and. idnode == 0) &
     Write(nrite,"(/,1x,'real space cutoff reset to (Ang)',2x,1p,e12.4)") rcut

! report vdw

  If (l_n_v .and. idnode == 0) Write(nrite,"(/,1x,'vdw potential terms switched off')")

! report if rvdw is reset (measures taken in scan_config)

  If ((.not.l_n_v) .and. Abs(rvdw-rvdw1) > 1.0e-6_wp .and. idnode == 0) &
     Write(nrite,"(/,1x,'vdw cutoff reset to (Ang) ',8x,1p,e12.4)") rvdw

! report timestep

  If (lvar) Then

     If (mxdis >= 2.5_wp*mndis .and. mndis > 0.0_wp) Then
        If (idnode == 0) Then
           Write(nrite,"(/,1x,'variable simulation timestep (ps)',1x,1p,e12.4)") tstep

           Write(nrite,"(/,1x,a,2(/,1x,a,7x,1p,e12.4))") &
           "controls for variable timestep",             &
           "minimum distance Dmin (Ang)",mndis,          &
           "maximum distance Dmax (Ang)",mxdis
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

  Else If (lstep) Then

     If (idnode == 0) &
        Write(nrite,"(/,1x,'fixed simulation timestep (ps)   ',1x,1p,e12.4)") tstep

  End If

! report rdf

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
        If (idnode == 0) Write(nrite,"(1x,'no rdf printing requested')")
     End If

     If (l_n_r) Then
        If (idnode == 0) Write(nrite,"(1x,'no rdf pairs specified in FIELD')")
     Else
        If (idnode == 0) Write(nrite,"(1x,'rdf pairs specified in FIELD')")
     End If

     If ((.not.lrdf) .or. l_n_r) Then
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
           Write(nrite,"(  1x,'z-density binsize (Angsroms)',6x,1p,e12.4)") rbin
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

! report data dumping interval and job times

  If (idnode == 0) Then
     Write(nrite,"(/,1x,'data dumping interval (steps)',2x,i10)") ndump
     Write(nrite,"(/,1x,'allocated job run time   (s)',6x,1p,e12.4)") timjob
     Write(nrite,"(/,1x,'allocated job close time (s)',6x,1p,e12.4)") timcls
  End If

! report replay history

  If (.not.lsim) Then
     lprdf=lrdf; lpzdn=lzdn ! Force printing
     If (idnode == 0) Then
        Write(nrite,"(/,1x,'*** HISTORY will be replayed (no actual simulation) ***', &
                    & /,1x,'*** and structural properties will be recalculated  ***')")
     End If

! abort if there's no structural property to recalculate

     If (.not.(lrdf .or. lzdn .or. ldef .or. lrsd)) Call error(580)

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

!!! ERROR CHECKS !!!
! Temperature

  If ((.not.ltemp) .or. (nstrun > 0 .and. temp < 1.0_wp)) Call error(380)

! Timestep

  If (.not.lstep) Call error(381)

! check settings in Langevin ensembles

  If ((keyens == 10 .or. keyens == 20 .or. keyens == 30) .and. &
      chi <= zero_plus) Call error(462)

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

End Subroutine read_control
