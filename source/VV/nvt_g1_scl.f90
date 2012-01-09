Subroutine nvt_g1_scl &
           (tstep,ceng,qmass,temp,gama,r_0,pmass,chip, &
           vxx,vyy,vzz,                                &
           rgdvxx,rgdvyy,rgdvzz,                       &
           rgdoxx,rgdoyy,rgdozz,                       &
           chit,cint,engke,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NVT Nose-Hoover thermostat
! with a Langevin process when singled RBs are present
!
! Note: coupling to NPT barostat included as factor=pmass*chip^2
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module,       Only : nfree,lstfre
  Use rigid_bodies_module, Only : ntrgd
  Use kinetic_module,      Only : getknf,getknt,getknr

  Implicit None

  Real( Kind = wp ),                        Intent( In    ) :: tstep,ceng,qmass, &
                                                               temp,gama,r_0,    &
                                                               pmass,chip
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ), Dimension( 1:mxrgd ),  Intent( InOut ) :: rgdvxx,rgdvyy,rgdvzz
  Real( Kind = wp ), Dimension( 1:mxrgd ),  Intent( InOut ) :: rgdoxx,rgdoyy,rgdozz
  Real( Kind = wp ),                        Intent( InOut ) :: chit,cint
  Real( Kind = wp ),                        Intent(   Out ) :: engke,engrot

  Integer           :: i,j,irgd
  Real( Kind = wp ) :: engkf,engkt,hstep,qstep,factor,scale,fex


! timestep derivative and factor

  hstep  = 0.5_wp*tstep
  qstep  = 0.5_wp*hstep
  factor = pmass*chip**2

! update chi(=cint) to 1/4*tstep

  cint=cint + qstep*chit

! calculate kinetic energy contributions and rotational energy

  engkf=getknf(vxx,vyy,vzz)
  engkt=getknt(rgdvxx,rgdvyy,rgdvzz)

  engke=engkf+engkt

  engrot=getknr(rgdoxx,rgdoyy,rgdozz)

  fex=Exp(-gama*hstep)

! update chit to 1/2*tstep

  chit=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
       hstep*(2.0_wp*(engke+engrot) + factor - ceng)/qmass

! update chi(=cint) to 3/4*tstep

  cint=cint + hstep*chit

! thermostat the velocities to 1*tstep

  scale=Exp(-tstep*chit)
  Do j=1,nfree
     i=lstfre(j)

     vxx(i)=scale*vxx(i)
     vyy(i)=scale*vyy(i)
     vzz(i)=scale*vzz(i)
  End Do

  Do irgd=1,ntrgd
     rgdvxx(irgd)=scale*rgdvxx(irgd)
     rgdvyy(irgd)=scale*rgdvyy(irgd)
     rgdvzz(irgd)=scale*rgdvzz(irgd)

     rgdoxx(irgd)=scale*rgdoxx(irgd)
     rgdoyy(irgd)=scale*rgdoyy(irgd)
     rgdozz(irgd)=scale*rgdozz(irgd)
  End Do

! thermostat the energy consequently

  engke=engke*scale**2
  engrot=engrot*scale**2

! update chit to full (2/2)*tstep

  chit=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
       hstep*(2.0_wp*(engke+engrot) + factor - ceng)/qmass

! update chi(=cint) to 4/4*tstep

  cint=cint + qstep*chit

End Subroutine nvt_g1_scl
