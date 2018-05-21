Module langevin_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring:
!           1) Langevin npt, nst and nvt_lfv ensembles switches & arrays
!           2) gentle stochastic ensembles switch & Gaussian random number
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
! amended   - i.t.todorov march 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: l_lan    = .false., &
                                          l_gst    = .false.

  Real( Kind = wp ),              Save :: fpl(1:9) = 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: fxl(:),fyl(:),fzl(:)

  Public :: langevin_allocate_arrays

Contains

  Subroutine langevin_allocate_arrays

    Use setup_module, Only : mxatms

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (fxl(1:mxatms),fyl(1:mxatms),fzl(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1041)

    fxl = 0.0_wp ; fyl = 0.0_wp ; fzl = 0.0_wp

  End Subroutine langevin_allocate_arrays

End Module langevin_module