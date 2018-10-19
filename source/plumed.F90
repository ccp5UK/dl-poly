Module plumed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global PLUMED variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - a.m.elena september 2015
  ! contrib   - i.t.todorov march 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use comms,  Only : comms_type
  Use constants,  Only : boltz,  DLP_VERSION
  Use configuration, Only : configuration_type
  Use errors_warnings, Only : error,warning,info
  Use statistics, Only : stats_type
  Use particle, Only : corePart

  Implicit None

  Private

  !> Type to store PLUMED variables
  Type, Public :: plumed_type
    Private

    !> PLUMED switch
    Logical,                Public :: l_plumed = .false.
    Character( Len = 125 ), Public :: input    = "PLUMED"
    Character( Len = 125 ), Public :: logfile  = "OUTPUT.PLUMED"

    !> DL_POLY precision
    Integer, Public :: prec = wp
    !> default no
    Integer( Kind = wi ),   Public :: restart     = 0

    !> PLUMED energy
    Real( Kind = wp ) :: eng
    !> PLUMED virial
    Real( Kind = wp ) :: virial(1:9)

    Integer( Kind = wi ) :: version    = 0, &
      stop       = 0, &
      has_plumed = 0
  End Type plumed_type

  ! PLUMED parameters
  !> DLPOLY_Internal(10J/mol) /kJ/mol
  Real( Kind = wp ), Parameter, Public :: plumed_energyUnits = 0.01_wp
  !> Angstrtom/nanometer
  Real( Kind = wp ), Parameter, Public :: plumed_lengthUnits = 0.1_wp
  !> picosecond
  Real( Kind = wp ), Parameter, Public :: plumed_timeUnits   = 1.0_wp

  Character( Len =   1 ), Parameter  :: sn=Char(0)

  Private :: plumed_print_about
  Private :: plumed_message
  Public  :: plumed_init
  Public  :: plumed_finalize
  Public  :: plumed_apply

Contains

  Subroutine plumed_init(megatm,tstep,temp,plume,comm)

    Integer,           Intent( In    ) :: megatm
    Real( Kind = wp ), Intent( In    ) :: tstep,temp
    Type(plumed_type), Intent( InOut ) :: plume
    Type(comms_type),  Intent( InOut ) :: comm

#ifdef PLUMED
    Call plumed_f_installed(plume%has_plumed)

    If (plume%has_plumed > 0) Then
      Call plumed_f_gcreate()
      Call plumed_f_gcmd("getApiVersion"//sn,plume%version)
      Call plumed_f_gcmd("setRealPrecision"//sn,plume%prec)
      Call plumed_f_gcmd("setMDEnergyUnits"//sn,plumed_energyUnits)
      Call plumed_f_gcmd("setMDLengthUnits"//sn,plumed_lengthUnits)
      Call plumed_f_gcmd("setMDTimeUnits"//sn,plumed_timeUnits)
      Call plumed_f_gcmd("setMPIFComm"//sn,comm%comm)
      ! Ideally would change file names here into names that can be controlled by user
      ! from control
      Call plumed_f_gcmd("setPlumedDat"//sn,Trim(plume%input)//sn)
      Call plumed_f_gcmd("setLogFile"//sn,Trim(plume%logfile)//sn)
      Call plumed_f_gcmd("setNatoms"//sn,megatm)
      ! The name should be updated when there are new releases of dlpoly
      Call plumed_f_gcmd("setMDEngine"//sn,"DL_POLY "//DLP_VERSION//sn)
      Call plumed_f_gcmd("setTimestep"//sn,tstep)
      Call plumed_f_gcmd("setKbT"//sn,temp*boltz)
      Call plumed_f_gcmd("init"//sn,0)
    Else
      Call error(0,'internal PLUMED library failure',.true.)
    End If
#endif

    Call plumed_print_about(plume,comm)

  End Subroutine plumed_init

  Subroutine plumed_print_about(plume,comm)
    Type(plumed_type), Intent( In    ) :: plume
    Type(comms_type),  Intent( InOut ) :: comm
#ifdef PLUMED
    Character( Len = 256 ) :: message,messages(9),banner(15)

    Write(banner(1),'(a)')  ""
    Write(banner(2),'(a)')  "***_____________________________________ "
    Write(banner(3),'(a)')  "***|     .--                            |"
    Write(banner(4),'(a)')  "***|  _/ o)  \                          |"
    Write(banner(5),'(a)')  "***| `""'.'=- |                          |"
    Write(banner(6),'(a)')  "***|     )   (                          |"
    Write(banner(7),'(a)')  "***|    /     `=._                      |"
    Write(banner(8),'(a)')  "***|   |    .-   .`=._                  |"
    Write(banner(9),'(a)')  "***|   |   ;   .' ' _.:===.             |"
    Write(banner(10),'(a)') "***|   ',   '-====""``""""~^`'=.,          |"
    Write(banner(11),'(a)') "***|     '.,_____,,..=-""'""""','=,       |"
    Write(banner(12),'(a)') "***|      /` /`               ', '= .   |"
    Write(banner(13),'(a)') "***|     /,=='-,                 ''= )  |"
    Write(banner(14),'(a)') "***|   ;==`-,                           |"
    Write(banner(15),'(a)') "***--------------------------------------"
    Call info(banner,15,.true.)

    Write(messages(1),'(a)')        "*** Activating PLUMED Extension. ***"
    Write(messages(2),'(a)')        "*** Using PLUMED input file: "//Trim(plume%input)
    Write(messages(3),'(a)')        "*** Using PLUMED log file: "//Trim(plume%logfile)
    Write(messages(4),'(a,i0)')     "*** Using PLUMED API version: ",plume%version
    Write(messages(5),'(a,i0)')     "*** Using PLUMED Real precision: ", plume%prec
    Write(messages(6),'(a,es15.6)') "*** Using PLUMED energy conversion factor: ", plumed_energyUnits
    Write(messages(7),'(a,es15.6)') "*** Using PLUMED length conversion factor: ", plumed_lengthUnits
    Write(messages(8),'(a,es15.6)') "*** Using PLUMED time conversion factor: ", plumed_timeUnits
    Write(messages(9),'(a,i0)')     "*** Using PLUMED restart (0: no, 1: yes): ", plume%restart
    Call info(messages,9,.true.)
#else
    Call plumed_message()
#endif

  End Subroutine plumed_print_about

  Subroutine plumed_apply(config,nstrun,nstep,stats,plume,comm)

    Integer,           Intent( In    ) :: nstep
    Integer,           Intent(   Out ) :: nstrun

    Type(configuration_type),    Intent( InOut ) :: config
    Type(stats_type),            Intent( InOut ) :: stats
    Type(plumed_type),           Intent( InOut ) :: plume
    Type(comms_type),            Intent( InOut ) :: comm

#ifdef PLUMED
    Character( Len = 256 ) :: message
    Real( Kind = wp ),Dimension(:), Allocatable :: tx,ty,tz,tfx,tfy,tfz,tchge
    Integer :: fail(1:2), i
    Allocate(tx(1:config%mxatms),ty(1:config%mxatms),tz(1:config%mxatms),tchge(1:config%mxatms),stat=fail(1))
    Allocate(tfx(1:config%mxatms),tfy(1:config%mxatms),tfz(1:config%mxatms),stat=fail(2))
    !    If(Any(fail)) Call error(0)
    Do i = 1,config%mxatms
      tx(i) = config%parts(i)%xxx
      ty(i) = config%parts(i)%yyy
      tz(i) = config%parts(i)%zzz
      tfx(i) = config%parts(i)%fxx
      tfy(i) = config%parts(i)%fyy
      tfz(i) = config%parts(i)%fzz
      tchge(i) = config%parts(i)%chge
    End Do

    Call plumed_f_gcmd("setAtomsNlocal"//sn,config%natms)
    Call plumed_f_gcmd("setAtomsFGatindex"//sn,config%ltg)
    Call plumed_f_gcmd("setStep"//sn,nstep)
    Call plumed_f_gcmd("setMasses"//sn,config%weight)
    Call plumed_f_gcmd("setCharges"//sn,tchge)
    Call plumed_f_gcmd("setPositionsX"//sn,tx)
    Call plumed_f_gcmd("setPositionsY"//sn,ty)
    Call plumed_f_gcmd("setPositionsZ"//sn,tz)
    Call plumed_f_gcmd("setBox"//sn,config%cell)
    plume%eng = stats%stpcfg / real(comm%mxnode)
    Call plumed_f_gcmd("setEnergy"//sn,plume%eng)
    Call plumed_f_gcmd("setForcesX"//sn,tfx)
    Call plumed_f_gcmd("setForcesY"//sn,tfy)
    Call plumed_f_gcmd("setForcesZ"//sn,tfz)
    plume%virial = -stats%stress
    Call plumed_f_gcmd("setVirial"//sn,plume%virial)
    stats%stress = -plume%virial
    Call plumed_f_gcmd("setStopFlag"//sn,plume%stop)
    Call plumed_f_gcmd("calc"//sn )

    If (plume%stop /= 0) Then
      Write(message,'(a,i0)') 'DL_POLY was stopped cleanly by PLUMED at step: ',nstep
      Call warning(message,.true.)
      nstrun=nstep
    End If
    Do i=1,config%mxatms
      config%parts(i)%xxx=tx(i)
      config%parts(i)%yyy=ty(i)
      config%parts(i)%zzz=tz(i)
      config%parts(i)%fxx=tfx(i)
      config%parts(i)%fyy=tfy(i)
      config%parts(i)%fzz=tfz(i)
      config%parts(i)%chge=tchge(i)
    End Do
    Deallocate(tx,ty,tz,tfx,tfy,tfz,tchge)
#else
    nstrun=nstep

    Call plumed_message()
#endif

  End Subroutine plumed_apply

  Subroutine plumed_finalize()

#ifdef PLUMED
    Call plumed_f_gfinalize()
#endif

  End Subroutine plumed_finalize

  Subroutine plumed_message()
#ifndef PLUMED
    Call error(0,'PLUMED directive found in CONTROL but PLUMED not available',.true.)
#endif
  End Subroutine plumed_message

End Module plumed
