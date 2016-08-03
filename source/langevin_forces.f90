Subroutine langevin_forces(nstep,temp,tstep,chi,fxr,fyr,fzr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to generate Langevin random forces consistent
! with the target temperature and the Langevin thermostat friction
!
! Note: (1) Random forces do not contribute to the stress and virial
!           of the system they are accounted by the system pressure.
!       (2) Random forces do not apply to frozen and massless particles
!           as well as shells.
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,      Only : boltz,mxatms
  Use site_module,       Only : ntpshl,unqshl
  Use config_module,     Only : natms,ltg,lfrzn,atmnam,weight

  Implicit None

  Integer          , Intent( In    ) :: nstep
  Real( Kind = wp ), Intent( In    ) :: temp,tstep,chi
  Real( Kind = wp ), Intent(   Out ) :: fxr(1:mxatms),fyr(1:mxatms),fzr(1:mxatms)

  Integer           :: i
  Real( Kind = wp ) :: scale,tmp

! Get scaler to target variance*Sqrt(weight)

  scale = Sqrt(2.0_wp * chi * boltz * temp / tstep)

! Make variance = target variance and nullify the rest

  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. (.not.Any(unqshl(1:ntpshl) == atmnam(i)))) Then
        Call box_mueller_saru3(ltg(i),nstep,fxr(i),fyr(i),fzr(i))

        tmp = scale*Sqrt(weight(i))

        fxr(i) = fxr(i)*tmp
        fyr(i) = fyr(i)*tmp
        fzr(i) = fzr(i)*tmp
     Else
        fxr(i) = 0.0_wp
        fyr(i) = 0.0_wp
        fzr(i) = 0.0_wp
     End If
  End Do

End Subroutine langevin_forces
