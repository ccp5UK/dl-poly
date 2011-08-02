Subroutine npt_h1_scl &
           (sw,tstep,degfre,degrot,pmass,chit,volm,press,vircon,virtot,vircom, &
           vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,chip,engke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NPT barostat
! when singled RBs are present
!
! sw=1 coupling to NVT thermostat for MTK ensembles,
!                                   factor=3/(degfre-degrot)
!
! sw=0 coupling to NVT thermostat for Nose-Hoover ensembles,
!                                   factor=0
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module,       Only : nfree,lstfre
  Use rigid_bodies_module, Only : ntrgd

  Implicit None

  Integer,           Intent( In    ) :: sw
  Integer(Kind=ip),  Intent( In    ) :: degfre,degrot

  Real( Kind = wp ), Intent( In    ) :: tstep,pmass,chit,volm,press,vircon,virtot,vircom
  Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: rgdvxx(1:mxrgd),rgdvyy(1:mxrgd),rgdvzz(1:mxrgd)
  Real( Kind = wp ), Intent( InOut ) :: chip,engke


  Logical,           Save :: newjob = .true.

  Integer                 :: i,j,irgd

  Real( Kind = wp ), Save :: factor
  Real( Kind = wp )       :: scale

  Real( Kind = wp ) :: hstep,qstep


  If (newjob) Then
     newjob = .false.

     factor = 0.0_wp
     If (sw == 1) factor = 3.0_wp/Real(degfre-degrot,wp)
  End If

! timestep derivatives

  hstep=0.5_wp*tstep
  qstep=0.5_wp*hstep

! thermostat chip to 1/4*tstep

  chip = chip*Exp(-qstep*chit)

! barostat chip to 1/2*tstep

  chip = chip + hstep*( (2.0_wp*(1.0_wp+factor)*engke-vircon-virtot-vircom) - 3.0_wp*press*volm ) / pmass

! thermostat chip to 2/4*tstep

  chip = chip*Exp(-qstep*chit)

! barostat the velocities to full 1*tstep

  scale=Exp(-tstep*(1.0_wp+factor)*chip)
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
  End Do

! barostat the energy consequently

  engke=engke*scale**2

! thermostat chip to 3/4*tstep

  chip = chip*Exp(-qstep*chit)

! barostat chip to full (1/2 + 1/2)*tstep

  chip = chip + hstep*( (2.0_wp*(1.0_wp+factor)*engke-vircon-virtot-vircom) - 3.0_wp*press*volm ) / pmass

! thermostat chip to full (4/4)*tstep

  chip = chip*Exp(-qstep*chit)

End Subroutine npt_h1_scl
