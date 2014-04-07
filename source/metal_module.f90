Module metal_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global metal interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: ld_met = .false., & ! no direct calculations are opted
                                          ls_met = .false., & ! no embedding over Sqrt(rho) but over rho
                                          l2bmet = .false.    ! no 2B(EAM or EEAM)

  Integer,                        Save :: ntpmet = 0 , & ! number of different metal interactions
                                          tabmet = -1    ! undefined, 0 - no TABEAM, 1 - EAM, 2 - EEAM, 3- 2BEAM, 4 - 2BEEAM


  Integer,           Allocatable, Save :: lstmet(:),ltpmet(:)

  Real( Kind = wp ), Allocatable, Save :: prmmet(:,:)

  Real( Kind = wp ), Allocatable, Save :: elrcm(:),vlrcm(:)

! Possible tabulated calculation arrays

  Real( Kind = wp ), Allocatable, Save :: vmet(:,:,:), dmet(:,:,:),dmes(:,:,:), &
                                          fmet(:,:,:),fmes(:,:,:)

! Atomic density [reused as embedding derivative(s)] helper array(s)

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: rho,rhs

  Public :: allocate_metal_arrays, &
            allocate_metal_table_arrays

Contains

  Subroutine allocate_metal_arrays()

    Use setup_module, Only : mxmet,mxpmet,mxatms,mxatyp

    Integer, Dimension( 1:7 ) :: fail

    If (tabmet == 3 .or. tabmet == 4) l2bmet=.true.

    fail = 0

    Allocate (lstmet(1:mxmet),                  Stat = fail(1))
    Allocate (ltpmet(1:mxmet),                  Stat = fail(2))
    Allocate (prmmet(1:mxpmet,1:mxmet),         Stat = fail(3))
    Allocate (rho(1:Merge(mxatms,0,mxmet > 0)), Stat = fail(4))
    If (l2bmet) & ! the new S-band density
    Allocate (rhs(1:Merge(mxatms,0,mxmet > 0)), Stat = fail(5))
    Allocate (elrcm(0:mxatyp),                  Stat = fail(6))
    Allocate (vlrcm(0:mxatyp),                  Stat = fail(7))

    If (Any(fail > 0)) Call error(1023)

    lstmet = 0
    ltpmet = 0

    prmmet = 0.0_wp

! rho and rhs get initialised in metal_ld_compute!!!

    elrcm  = 0.0_wp
    vlrcm  = 0.0_wp

  End Subroutine allocate_metal_arrays

  Subroutine allocate_metal_table_arrays()

    Use setup_module, Only : mxmet,mxmed,mxmds,mxatyp,mxgmet

    Implicit None

    Integer, Dimension( 1:5 ) :: fail

    fail = 0

    Allocate (vmet(1:mxgmet,1:mxmet, 1:2),    Stat = fail(1))
    Allocate (dmet(1:mxgmet,1:mxmed, 1:2),    Stat = fail(2))
    Allocate (fmet(1:mxgmet,1:mxatyp,1:2),    Stat = fail(3))
    If (tabmet == 3 .or. tabmet == 4) Then ! the new S-band density and embedding
       Allocate (dmes(1:mxgmet,1:mxmds, 1:2), Stat = fail(4))
       Allocate (fmes(1:mxgmet,1:mxatyp,1:2), Stat = fail(5))
    End If

    If (Any(fail > 0)) Call error(1069)

    vmet = 0.0_wp
    dmet = 0.0_wp
    fmet = 0.0_wp
    If (tabmet == 3 .or. tabmet == 4) Then ! the new S-band density and embedding
      dmes = 0.0_wp
      fmes = 0.0_wp
    End If

  End Subroutine allocate_metal_table_arrays

End Module metal_module
