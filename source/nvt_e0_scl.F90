Subroutine nvt_e0_scl(isw,tstep,fxx,fyy,fzz,vxx,vyy,vzz,chit,engke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NVT Evans thermostat
! when no singled rigid bodies are present
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,  Only : comms_type, gsu,
  Use setup,  Only : mxatms
  Use configuration, Only : natms,weight

  Implicit None

  Integer,                                  Intent( In    ) :: isw
  Real( Kind = wp ),                        Intent( In    ) :: tstep
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: fxx,fyy,fzz
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ),                        Intent(   Out ) :: chit,engke

  Integer           :: i
  Real( Kind = wp ) :: buffer(1:2),vdotf,scale

  engke = 0.0_wp
  vdotf = 0.0_wp
  Do i=1,natms
     engke = engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
     vdotf = vdotf+vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
  End Do

  If (mxnode > 1) Then
     buffer(1) = engke
     buffer(2) = vdotf
     Call gsum(buffer)
     engke = buffer(1)
     vdotf = buffer(2)
  End If

! velocity friction and temperature scaling coefficient
! for Evans thermostat at tstep

  chit = vdotf/engke

! get corrected energy

  engke = 0.5_wp*engke

  If (isw == 0) Return

! thermostat velocities

  scale=Exp(-chit*tstep)
  Do i=1,natms
     vxx(i)=vxx(i)*scale
     vyy(i)=vyy(i)*scale
     vzz(i)=vzz(i)*scale
  End Do

! thermostat kinetic energy

  engke = engke*scale**2

End Subroutine nvt_e0_scl
