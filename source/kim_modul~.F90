! This file is generated manually from kim_module.F90 by
! gfortran -E -P kim_module.F90 > kim_modul~.F90








Module kim_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global KIM interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - r.s.elliott march 2015
! contrib   - h.boateng & i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use, Intrinsic :: iso_c_binding
  Use kinds_f90

  Implicit None

  Character( Len = 200 ), Save :: kim  = ' '      ! KIM IM type for dl_poly
  Real( Kind = wp ),      Save :: rkim = 0.0_wp   ! KIM cutoff for dl_poly
  Integer,                Save :: idhalo(0:2,1:6) ! KIM halo indicator

  Public :: kim_cutoff,  &
            kim_setup,   &
            kim_forces,  &
            kim_cleanup


Contains

  Subroutine kim_cutoff(num_types,model_types,model_name,cutoff)

!-------------------------------------------------------------------------------
!
! kim_cutoff
!
! This function extracts the cutoff distance for the specified KIM model_name
!
!-------------------------------------------------------------------------------

    Implicit None

    Integer( Kind = c_int ), Intent( In    ) :: num_types
    Character( Len = * ),    Intent( In    ) :: model_types(1:num_types)
    Character( Len = * ),    Intent( In    ) :: model_name
    Real( Kind = wp ),       Intent(   Out ) :: cutoff

    cutoff = 0.0_wp
    Call kim_message()
  End Subroutine  kim_cutoff

  Subroutine kim_setup(num_types,model_types,model_name)

!-------------------------------------------------------------------------------
!
! kim_setup
!
! This function creates the interface via the KIM API to KIM model_name
!
!-------------------------------------------------------------------------------

    Implicit None

    Character(Len = *),      Intent( In    ) :: model_name
    Integer( Kind = c_int ), Intent( In    ) :: num_types
    Character( Len = * ),    Intent( In    ) :: model_types(1:num_types)

    Call kim_message()
  End Subroutine kim_setup

  Subroutine kim_cleanup()

!-------------------------------------------------------------------------------
!
! kim_cleanup
!
! This function releases resources and cleans up the KIM API interface
!
!-------------------------------------------------------------------------------

    Implicit None

    Call kim_message()
  End Subroutine kim_cleanup

  Subroutine kim_forces(engkim,virkim,stress)

!-------------------------------------------------------------------------------
!
! kim_forces
!
! This function updates the data needed by the KIM Model and then asks for
! the energy, forces, and virial from the KIM Model, distributes force
! contributions appropriately to all the processors and updates virial
! and pressure values.
!
!-------------------------------------------------------------------------------

    Implicit None

    Real( Kind = wp ), Intent( InOut ) :: engkim
    Real( Kind = wp ), Intent( InOut ) :: virkim
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)

    engkim = 0.0_wp
    virkim = 0.0_wp
    stress = 0.0_wp

    Call kim_message()
  End Subroutine kim_forces


  Subroutine kim_message()

    Use comms_module, Only : idnode
    Use setup_module, Only : nrite

    Implicit None

    If (idnode == 0) Write(nrite,'(1x,a)') "No openKIM crosscompiled"
    Call error(0)
  End Subroutine kim_message

End Module kim_module
