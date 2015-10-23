Module plumed_module
  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,dlp_comm_world
  Use setup_module,  Only : nrite, boltz, mxatms, dl_poly_version
  Use config_module, Only : cell,natms,weight,ltg,chge,fxx,fyy,fzz
  Implicit None

  private

  Logical, public                :: l_plumed = .false.
  Character( len = 125 ), public :: plumed_input = "PLUMED"
  Character( len = 125 ), public :: plumed_log = "OUTPUT.PLUMED"

  Real ( Kind = wp ), public     :: plumed_energyUnits = 0.01_wp ! DLPOLY_Internal(10J/mol) /kJ/mol
  Real ( Kind = wp ), public     :: plumed_lengthUnits = 0.1_wp ! Angstrtom/nanometer
  Real ( Kind = wp ), public     :: plumed_timeUnits = 1.0_wp ! picosecond
  Integer ( Kind = ip ), public  :: plumed_precision = wp ! DL_POLY preicision
  Integer ( Kind = ip ), public  :: plumed_restart = 0 ! default no

  Real ( Kind = wp )  :: plumed_eng 
  Integer ( Kind=ip ) :: plumed_version=0, plumed_stop=0,has_plumed=0
  Real ( Kind = wp )  :: plumed_virial(1:9)
  character(len=1),parameter  :: sn=char(0)

  public :: plumed_print_about
  public :: plumed_init
  public :: plumed_finalize
  public :: plumed_apply
Contains 
  
  Subroutine plumed_init(megatm,tstep,temp)

    Integer, Intent(In)           :: megatm
    Real(kind = wp ), Intent(in)  :: tstep, temp
#ifdef PLUMED
    Call plumed_f_installed(has_plumed)

    if (has_plumed > 0) then
      Call plumed_f_gcreate()
      Call plumed_f_gcmd("getApiVersion"//sn,plumed_version)
      Call plumed_f_gcmd("setRealPrecision"//sn,plumed_precision)
      Call plumed_f_gcmd("setMDEnergyUnits"//sn,plumed_energyUnits)
      Call plumed_f_gcmd("setMDLengthUnits"//sn,plumed_lengthUnits)
      Call plumed_f_gcmd("setMDTimeUnits"//sn,plumed_timeUnits)
      Call plumed_f_gcmd("setMPIFComm"//sn,dlp_comm_world)
! Ideally would change file names here into names that can be controlled by user
! from control
      Call plumed_f_gcmd("setPlumedDat"//sn,trim(plumed_input)//sn)
      Call plumed_f_gcmd("setLogFile"//sn,trim(plumed_log)//sn)
      Call plumed_f_gcmd("setNatoms"//sn,megatm)
! The name should be updated when there are new releases of dlpoly
      Call plumed_f_gcmd("setMDEngine"//sn,"DL_POLY "//dl_poly_version//sn)
      Call plumed_f_gcmd("setTimestep"//sn,tstep)
      call plumed_f_gcmd("setKbT"//sn,temp*boltz)  
      Call plumed_f_gcmd("init"//sn,0)
    Else
       Call error(1083)
    End If
#endif
  End Subroutine plumed_init

  Subroutine plumed_print_about()

#ifdef PLUMED
   If (idnode == 0) Then
     Write(nrite,'(a)')""
     Write(nrite,'(14(a42,/))')&
       "***_____________________________________ ",&
       "***|     .--                            |",&
       "***|  _/ o)  \                          |",&
       "***| `""'.'=- |                          |",&
       "***|     )   (                          |",&
       "***|    /     `=._                      |",&
       "***|   |    .-   .`=._                  |",&
       "***|   |   ;   .' ' _.:===.             |",&
       "***|   ',   '-====""``""""~^`'=.,          |",&
       "***|     '.,_____,,,..=-""'""""','=,       |",&
       "***|      /` /`               ', '= .   |",&
       "***|     /,=='-,                 ''= )  |",&
       "***|   ;==`-,                           |",&
       "***--------------------------------------"
     Write(nrite,'(a)') "***Activating PLUMED Extension.***"
     Write(nrite,'(a)') "***Using PLUMED input file: "//trim(plumed_input)
     Write(nrite,'(a)') "***Using PLUMED log file: "//trim(plumed_log)
     Write(nrite,'(a,i0)') "***Using PLUMED API version: ",plumed_version
     Write(nrite,'(a,i0)') "***Using PLUMED Real precision: ",plumed_precision
     Write(nrite,'(a,es15.6)') "***Using PLUMED energy conversion factor: ",plumed_energyUnits
     Write(nrite,'(a,es15.6)') "***Using PLUMED length conversion factor: ",plumed_lengthUnits
     Write(nrite,'(a,es15.6)') "***Using PLUMED time conversion factor: ",plumed_timeUnits
     Write(nrite,'(a,i0)') "***Using PLUMED restart (0: no, 1: yes): ",plumed_restart
   End If
#endif
    
  End Subroutine plumed_print_about

  Subroutine plumed_apply(xxx,yyy,zzz,nstrun,nstep,stpcfg,stress)
     Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
     Real( Kind = wp ), Intent( InOut ) :: stress(1:9),stpcfg
     Integer, Intent( In    )            :: nstep
     Integer, Intent(    Out)            :: nstrun

#ifdef PLUMED
     call plumed_f_gcmd("setAtomsNlocal"//sn,natms)
     call plumed_f_gcmd("setAtomsFGatindex"//sn,ltg)
     call plumed_f_gcmd("setStep"//sn,nstep)
     call plumed_f_gcmd("setMasses"//sn,weight)
     call plumed_f_gcmd("setCharges"//sn,chge)
     call plumed_f_gcmd("setPositionsX"//sn,xxx)
     call plumed_f_gcmd("setPositionsY"//sn,yyy)
     call plumed_f_gcmd("setPositionsZ"//sn,zzz)
     call plumed_f_gcmd("setBox"//sn,cell)
     plumed_eng = stpcfg / real(mxnode)
     call plumed_f_gcmd("setEnergy"//sn,plumed_eng)
     call plumed_f_gcmd("setForcesX"//sn,fxx)
     call plumed_f_gcmd("setForcesY"//sn,fyy)
     call plumed_f_gcmd("setForcesZ"//sn,fzz)
     plumed_virial = -stress
     call plumed_f_gcmd("setVirial"//sn,plumed_virial)
     stress = -plumed_virial
     call plumed_f_gcmd("setStopFlag"//sn,plumed_stop)
     call plumed_f_gcmd("calc"//sn )

     If( plumed_stop /= 0 ) then
       If(idnode == 0) write(nrite,'(a,i0,a)')"*** warning DL_POLY was stopped cleanly by PLUMED at step: ",nstep," *** "
       nstrun=nstep
     End If
#endif
  End Subroutine plumed_apply

  Subroutine plumed_finalize()
#ifdef PLUMED
    If ( idnode == 0 ) Call plumed_f_gfinalize()
#endif
  End Subroutine plumed_finalize

End Module plumed_module
