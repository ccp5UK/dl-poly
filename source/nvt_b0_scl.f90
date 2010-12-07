Subroutine nvt_b0_scl(isw,tstep,sigma,taut,vxx,vyy,vzz,chit,strkin,engke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NVT Berendsen thermostat
! when no singled rigid bodies are present
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,   Only : mxatms
  Use config_module,  Only : natms
  Use kinetic_module, Only : kinstress

  Implicit None

  Integer,                                  Intent( In    ) :: isw
  Real( Kind = wp ),                        Intent( In    ) :: tstep,sigma,taut
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ), Dimension( 1:9 ),      Intent(   Out ) :: strkin
  Real( Kind = wp ),                        Intent(   Out ) :: chit,engke

  Integer           :: i
  Real( Kind = wp ) :: tmp

! get kinetic energy and stress

  Call kinstress(vxx,vyy,vzz,strkin)
  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! temperature scaling coefficient - taut is the decay constant

  chit=Sqrt(1.0_wp+tstep/taut*(sigma/engke-1.0_wp))

  If (isw == 0) Return

! thermostat velocities

  Do i=1,natms
     vxx(i)=vxx(i)*chit
     vyy(i)=vyy(i)*chit
     vzz(i)=vzz(i)*chit
  End Do

! thermostat kinetic stress and energy

  tmp=chit**2

  strkin=strkin*tmp
  engke=engke*tmp

End Subroutine nvt_b0_scl
