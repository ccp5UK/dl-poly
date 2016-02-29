Module kim_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module implementing the kim dummy

! copyright - daresbury laboratory
! amended   - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Character( Len = 200 ), Save :: kim  = ' '      ! KIM IM type for dl_poly
  Real( Kind = wp ),      Save :: rkim = 0.0_wp   ! KIM cutoff for dl_poly
  Integer,                Save :: idhalo(0:2,1:6) ! KIM halo indicator

  Public :: kim_cutoff,kim_setup,kim_forces,kim_cleanup

Contains

  Subroutine kim_cutoff(num_types,model_types,model_name,cutoff)

    Implicit None

    Integer,              Intent( In    ) :: num_types
    Character( Len = * ), Intent( In    ) :: model_types(1:num_types)
    Character( Len = * ), Intent( In    ) :: model_name
    Real(Kind = wp),      Intent(   Out ) :: cutoff

    cutoff = 0.0_wp

    Call kim_message()

  End Subroutine kim_cutoff

  Subroutine kim_setup(num_types,model_types,model_name)

    Implicit None

    Integer,              Intent( In    ) :: num_types
    Character( Len = * ), Intent( In    ) :: model_types(1:num_types)
    Character( Len = * ), Intent( In    ) :: model_name

    Call kim_message()

  End Subroutine kim_setup

  Subroutine kim_forces(engkim,virkim,stress)

    Implicit None

    Real(Kind = wp), Intent(   Out ) :: engkim,virkim,stress(1:9)

    engkim = 0.0_wp
    virkim = 0.0_wp
    stress = 0.0_wp

    Call kim_message()

  End Subroutine kim_forces

  Subroutine kim_cleanup()

    Implicit None

    Call kim_message()

  End Subroutine kim_cleanup

 Subroutine kim_message()

    Use comms_module, Only : idnode
    Use setup_module, Only : nrite

    Implicit None

    If (idnode == 0) Write(nrite,'(1x,a)') "No openKIM crosscompiled"
    Call error(0)

  End Subroutine kim_message

End Module kim_module
