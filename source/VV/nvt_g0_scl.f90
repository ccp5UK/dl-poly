Subroutine nvt_g0_scl &
           (tstep,ceng,qmass,temp,gama,r_0,pmass,chip, &
           vxx,vyy,vzz,chit,cint,engke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NVT Nose-Hoover thermostat
! with a Langevin process
!
! Note: coupling to NPT barostat included as factor=pmass*chip^2
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module,  Only : natms
  Use kinetic_module, Only : getkin

  Implicit None

  Real( Kind = wp ),                        Intent( In    ) :: tstep,ceng,qmass, &
                                                               temp,gama,r_0,    &
                                                               pmass,chip
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ),                        Intent( InOut ) :: chit,cint
  Real( Kind = wp ),                        Intent(   Out ) :: engke

  Integer           :: i
  Real( Kind = wp ) :: hstep,qstep,factor,scale,fex


! timestep derivative and factor

  hstep  = 0.5_wp*tstep
  qstep  = 0.5_wp*hstep
  factor = pmass*chip**2

! update chi(=cint) to 1/4*tstep

  cint=cint + qstep*chit

! calculate kinetic energy

  engke=getkin(vxx,vyy,vzz)

  fex=Exp(-gama*hstep)

! update chit to 1/2*tstep

  chit=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
       hstep*(2.0_wp*engke + factor - ceng)/qmass

! update chi(=cint) to 3/4*tstep

  cint=cint + hstep*chit

! thermostat the velocities to 1*tstep

  scale=Exp(-tstep*chit)
  Do i=1,natms
     vxx(i)=scale*vxx(i)
     vyy(i)=scale*vyy(i)
     vzz(i)=scale*vzz(i)
  End Do

! thermostat the energy consequently

  engke=engke*scale**2

! update chit to full (2/2)*tstep

  chit=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
       hstep*(2.0_wp*engke + factor - ceng)/qmass

! update chi(=cint) to 4/4*tstep

  cint=cint + qstep*chit

End Subroutine nvt_g0_scl
