!This file is manually from plumed_module.F90 by 
!gfortran -E plumed_module.F90 > plumed_modul~.F90
! some lines may start with #, replace it by !
! 1 "plumed_module.F90"
! 1 "<built-in>"
! 1 "<command-line>"
! 1 "plumed_module.F90"
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
! 60 "plumed_module.F90"
  End Subroutine plumed_init

  Subroutine plumed_print_about()

! 93 "plumed_module.F90"
    
  End Subroutine plumed_print_about

  Subroutine plumed_apply(xxx,yyy,zzz,nstrun,nstep,stpcfg,stress)
     Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
     Real( Kind = wp ), Intent( InOut ) :: stress(1:9),stpcfg
     Integer, Intent( In    )            :: nstep
     Integer, Intent(    Out)            :: nstrun

! 128 "plumed_module.F90"
  End Subroutine plumed_apply

  Subroutine plumed_finalize()



  End Subroutine plumed_finalize

End Module plumed_module
