! This file is generated manually from plumed_module.F90 by
! gfortran -E -P plumed_module.F90 > plumed_modul~.F90


Module plumed_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global PLUMED variables and arrays
!
! copyright - daresbury laboratory
! author    - a.m.elena september 2015
! contrib   - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode,mxnode,dlp_comm_world
  Use setup_module,  Only : nrite, boltz, mxatms, DLP_VERSION
  Use config_module, Only : cell,natms,weight,ltg,chge,fxx,fyy,fzz

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

  Subroutine plumed_init(megatm,tstep,temp)

    Integer,           Intent( In    ) :: megatm
    Real( Kind = wp ), Intent( In    ) :: tstep,temp


    Call plumed_print_about()

  End Subroutine plumed_init

  Subroutine plumed_print_about()

    Call plumed_message()

  End Subroutine plumed_print_about

  Subroutine plumed_apply(xxx,yyy,zzz,nstrun,nstep,stpcfg,stress)

    Integer,           Intent( In    ) :: nstep
    Integer,           Intent(   Out ) :: nstrun

    Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Real( Kind = wp ), Intent( In    ) :: stpcfg

    nstrun=nstep

    Call plumed_message()

  End Subroutine plumed_apply

  Subroutine plumed_finalize()


  End Subroutine plumed_finalize

  Subroutine plumed_message()

    Use comms_module, Only : idnode
    Use setup_module, Only : nrite

    Implicit None

    If (idnode == 0) Write(nrite,'(1x,a)') "*** warning - PLUMED directive found in CONTROL but PLUMED not available !!! ***"
    Call error(0)
  End Subroutine plumed_message

End Module plumed_module
