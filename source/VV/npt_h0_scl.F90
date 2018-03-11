Subroutine npt_h0_scl &
           (sw,tstep,degfre,pmass,chit,volm,press,vircon,virtot, &
           vxx,vyy,vzz,chip,engke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NPT barostat
!
! sw=1 coupling to NVT thermostat for MTK ensembles,
!                                   factor=3/degfre
!
! sw=0 coupling to NVT thermostat for Nose-Hoover ensembles,
!                                   factor=0
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, li
  Use setup_module
  Use configuration, Only : natms

  Implicit None

  Integer,           Intent( In    ) :: sw
  Integer(Kind=li),  Intent( In    ) :: degfre

  Real( Kind = wp ), Intent( In    ) :: tstep,pmass,chit,volm,press,vircon,virtot
  Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: chip,engke


  Logical,           Save :: newjob = .true.

  Integer                 :: i

  Real( Kind = wp ), Save :: factor
  Real( Kind = wp )       :: scale

  Real( Kind = wp ) :: hstep,qstep

  If (newjob) Then
     newjob = .false.

     factor = 0.0_wp
     If (sw == 1) factor = 3.0_wp/Real(degfre,wp)
  End If

! timestep derivatives

  hstep=0.5_wp*tstep
  qstep=0.5_wp*hstep

! thermostat chip to 1/4*tstep

  chip = chip*Exp(-qstep*chit)

! barostat chip to 1/2*tstep

  chip = chip + hstep*( (2.0_wp*(1.0_wp+factor)*engke-vircon-virtot) - 3.0_wp*press*volm ) / pmass

! thermostat chip to 2/4*tstep

  chip = chip*Exp(-qstep*chit)

! barostat the velocities to full 1*tstep

  scale=Exp(-tstep*(1.0_wp+factor)*chip)
  Do i=1,natms
     vxx(i)=scale*vxx(i)
     vyy(i)=scale*vyy(i)
     vzz(i)=scale*vzz(i)
  End Do

! barostat the energy consequently

  engke=engke*scale**2

! thermostat chip to 3/4*tstep

  chip = chip*Exp(-qstep*chit)

! barostat chip to full (1/2 + 1/2)*tstep

  chip = chip + hstep*( (2.0_wp*(1.0_wp+factor)*engke-vircon-virtot) - 3.0_wp*press*volm ) / pmass

! thermostat chip to full (4/4)*tstep

  chip = chip*Exp(-qstep*chit)

End Subroutine npt_h0_scl
