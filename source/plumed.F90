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

  Use kinds, Only : wp
  Use comms,  Only : comms_type
  Use setup,  Only : boltz, mxatms, DLP_VERSION
  Use configuration, Only : cell,natms,weight,ltg,chge,fxx,fyy,fzz
  Use errors_warnings, Only : error,warning,info
  Use statistics, Only : stats_type

  Implicit None

  Private

  Logical,                Public :: l_plumed           = .false.
  Character( Len = 125 ), Public :: plumed_input       = "PLUMED"
  Character( Len = 125 ), Public :: plumed_log         = "OUTPUT.PLUMED"

  Real( Kind = wp ),      Public :: plumed_energyUnits = 0.01_wp ! DLPOLY_Internal(10J/mol) /kJ/mol
  Real( Kind = wp ),      Public :: plumed_lengthUnits = 0.1_wp  ! Angstrtom/nanometer
  Real( Kind = wp ),      Public :: plumed_timeUnits   = 1.0_wp  ! picosecond
  Integer,                Public :: plumed_precision   = wp      ! DL_POLY precision
  Integer,                Public :: plumed_restart     = 0       ! default no

  Character( Len =   1 ), Parameter  :: sn=Char(0)

  Integer           :: plumed_version = 0, &
                       plumed_stop    = 0, &
                       has_plumed     = 0

  Real( Kind = wp ) :: plumed_eng,plumed_virial(1:9)

  Private :: plumed_print_about
  Private :: plumed_message
  Public  :: plumed_init
  Public  :: plumed_finalize
  Public  :: plumed_apply

Contains

  Subroutine plumed_init(megatm,tstep,temp,comm)

    Integer,           Intent( In    ) :: megatm
    Real( Kind = wp ), Intent( In    ) :: tstep,temp
    Type(comms_type),  Intent( InOut ) :: comm

#ifdef PLUMED
    Call plumed_f_installed(has_plumed)

    If (has_plumed > 0) Then
       Call plumed_f_gcreate()
       Call plumed_f_gcmd("getApiVersion"//sn,plumed_version)
       Call plumed_f_gcmd("setRealPrecision"//sn,plumed_precision)
       Call plumed_f_gcmd("setMDEnergyUnits"//sn,plumed_energyUnits)
       Call plumed_f_gcmd("setMDLengthUnits"//sn,plumed_lengthUnits)
       Call plumed_f_gcmd("setMDTimeUnits"//sn,plumed_timeUnits)
       Call plumed_f_gcmd("setMPIFComm"//sn,comm%comm)
! Ideally would change file names here into names that can be controlled by user
! from control
       Call plumed_f_gcmd("setPlumedDat"//sn,Trim(plumed_input)//sn)
       Call plumed_f_gcmd("setLogFile"//sn,Trim(plumed_log)//sn)
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

    Call plumed_print_about(comm)

  End Subroutine plumed_init

  Subroutine plumed_print_about(comm)

    Type(comms_type), Intent( InOut ) :: comm
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
    Write(messages(2),'(a)')        "*** Using PLUMED input file: "//Trim(plumed_input)
    Write(messages(3),'(a)')        "*** Using PLUMED log file: "//Trim(plumed_log)
    Write(messages(4),'(a,i0)')     "*** Using PLUMED API version: ",plumed_version
    Write(messages(5),'(a,i0)')     "*** Using PLUMED Real precision: ", plumed_precision
    Write(messages(6),'(a,es15.6)') "*** Using PLUMED energy conversion factor: ", plumed_energyUnits
    Write(messages(7),'(a,es15.6)') "*** Using PLUMED length conversion factor: ", plumed_lengthUnits
    Write(messages(8),'(a,es15.6)') "*** Using PLUMED time conversion factor: ", plumed_timeUnits
    Write(messages(9),'(a,i0)')     "*** Using PLUMED restart (0: no, 1: yes): ", plumed_restart
    Call info(messages,9,.true.)
#else
    Call plumed_message()
#endif

  End Subroutine plumed_print_about

  Subroutine plumed_apply(xxx,yyy,zzz,nstrun,nstep,stats,comm)

    Integer,           Intent( In    ) :: nstep
    Integer,           Intent(   Out ) :: nstrun

    Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
    Type(stats_type), Intent( InOut )  :: stats
    Type(comms_type), Intent( InOut )  :: comm

#ifdef PLUMED
    Character( Len = 256 ) :: message

    Call plumed_f_gcmd("setAtomsNlocal"//sn,natms)
    Call plumed_f_gcmd("setAtomsFGatindex"//sn,ltg)
    Call plumed_f_gcmd("setStep"//sn,nstep)
    Call plumed_f_gcmd("setMasses"//sn,weight)
    Call plumed_f_gcmd("setCharges"//sn,chge)
    Call plumed_f_gcmd("setPositionsX"//sn,xxx)
    Call plumed_f_gcmd("setPositionsY"//sn,yyy)
    Call plumed_f_gcmd("setPositionsZ"//sn,zzz)
    Call plumed_f_gcmd("setBox"//sn,cell)
    plumed_eng = stats%stpcfg / real(comm%mxnode)
    Call plumed_f_gcmd("setEnergy"//sn,plumed_eng)
    Call plumed_f_gcmd("setForcesX"//sn,fxx)
    Call plumed_f_gcmd("setForcesY"//sn,fyy)
    Call plumed_f_gcmd("setForcesZ"//sn,fzz)
    plumed_virial = -stats%stress
    Call plumed_f_gcmd("setVirial"//sn,plumed_virial)
    stats%stress = -plumed_virial
    Call plumed_f_gcmd("setStopFlag"//sn,plumed_stop)
    Call plumed_f_gcmd("calc"//sn )

    If (plumed_stop /= 0) Then
       Write(message,'(a,i0)') 'DL_POLY was stopped cleanly by PLUMED at step: ',nstep
       Call warning(message,.true.)
       nstrun=nstep
    End If
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
