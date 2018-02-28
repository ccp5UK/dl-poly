Subroutine metal_generate_erf(rmet)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating erf and fer arrays for
! many-body perturbation component only potentials
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module, Only : mxgmet,zero_plus
  Use metal_module, Only : ntpmet,ltpmet,prmmet,merf,mfer, &
                           allocate_metal_erf_arrays

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rmet

  Integer           :: imet
  Real( Kind = wp ) :: alpha,beta

! Determine alpha and beta for the erf bit of all MBPC potentials

  If (Any(ltpmet(1:ntpmet) == 5)) Then ! all are == 5 == MBPC
     alpha= 0.0_wp
     beta = 0.0_wp
     Do imet=1,ntpmet
        alpha=Max(alpha,Abs(prmmet(6,imet)))
        beta =Max(beta, Abs(prmmet(7,imet)))
     End Do

! If unset then set to defaults

     If (alpha <= zero_plus) alpha=20.0_wp
     If (beta  <= zero_plus) beta =Min(1.5_wp,0.2_wp*rmet)

! Allocate arrays: merf,mfer

     Call allocate_metal_erf_arrays()

! Generate error function and derivative arrays

     Call erfgen_met(rmet,alpha,beta,mxgmet,merf,mfer)

! Translate merf and mfer to the functional form 0.5*{1+erf[alpha(r-beta)]}

     merf(5:mxgmet)=0.5_wp*(1.0_wp+merf(5:mxgmet))
     mfer(5:mxgmet)=0.5_wp*mfer(5:mxgmet)
  End If

End Subroutine metal_generate_erf
