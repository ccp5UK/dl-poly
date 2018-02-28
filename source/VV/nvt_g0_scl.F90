Subroutine nvt_g0_scl &
           (tstep,degfre,isw,nstep,ceng,qmass,temp,gama,pmass,chip, &
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
! amended   - i.t.todorov march 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, li
  Use setup_module
  Use config_module,  Only : natms
  Use kinetic_module, Only : getkin

  Implicit None

  Real( Kind = wp ),                        Intent( In    ) :: tstep,ceng,qmass, &
                                                               temp,gama,pmass,chip
  Integer(Kind=li),                         Intent( In    ) :: degfre
  Integer,                                  Intent( In    ) :: isw,nstep
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ),                        Intent( InOut ) :: chit,cint
  Real( Kind = wp ),                        Intent(   Out ) :: engke

  Integer           :: i
  Real( Kind = wp ) :: hstep,qstep,factor,scale,fex,r_0


! timestep derivative and factor

  hstep  = 0.5_wp*tstep
  qstep  = 0.5_wp*hstep
  factor = pmass*chip**2

! update chi(=cint) to 1/4*tstep

  cint=cint + qstep*chit

! calculate kinetic energy

  engke=getkin(vxx,vyy,vzz)

  fex=Exp(-gama*hstep)

! generate a Gaussian random number for use in the
! Langevin process on the thermostat friction

  Call box_mueller_saru2(Int(degfre/3_li),nstep-1,2*isw+1,r_0,.true.)

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

! generate a Gaussian random number for use in the
! Langevin process on the thermostat friction

  Call box_mueller_saru2(Int(degfre/3_li),nstep-1,2*isw+2,r_0,.true.)

! update chit to full (2/2)*tstep

  chit=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
       hstep*(2.0_wp*engke + factor - ceng)/qmass

! update chi(=cint) to 4/4*tstep

  cint=cint + qstep*chit

End Subroutine nvt_g0_scl
