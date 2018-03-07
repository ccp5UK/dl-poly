Subroutine nvt_e1_scl &
           (isw,tstep,fxx,fyy,fzz,vxx,vyy,vzz,        &
           rgdfxx,rgdfyy,rgdfzz,rgdtxx,rgdtyy,rgdtzz, &
           rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz, &
           chit,engke,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NVT Evans thermostat
! when singled rigid bodies are present
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : mxnode,gsum
  Use setup_module,        Only : mxatms,mxrgd
  Use configuration,       Only : nfree,lstfre,weight
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,rgdwgt, &
                                  listrgd,indrgd,      &
                                  rgdrix,rgdriy,rgdriz

  Implicit None

  Integer,                                  Intent( In    ) :: isw
  Real( Kind = wp ),                        Intent( In    ) :: tstep
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: fxx,fyy,fzz
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ), Dimension( 1:mxrgd  ), Intent( In    ) :: rgdfxx,rgdfyy,rgdfzz, &
                                                               rgdtxx,rgdtyy,rgdtzz
  Real( Kind = wp ), Dimension( 1:mxrgd  ), Intent( InOut ) :: rgdvxx,rgdvyy,rgdvzz, &
                                                               rgdoxx,rgdoyy,rgdozz
  Real( Kind = wp ),                        Intent(   Out ) :: chit,engke,engrot

  Real( Kind = wp ) :: buffer(1:4),vdotf,odott,scale,tmp
  Integer           :: i,j,irgd,lrgd,rgdtyp

  engke = 0.0_wp
  vdotf = 0.0_wp
  engrot= 0.0_wp
  odott = 0.0_wp

! Free particles

  Do j=1,nfree
     i=lstfre(j)

     engke = engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
     vdotf = vdotf+vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
  End Do

! RBs

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        tmp=Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

        If (rgdfrz(0,rgdtyp) == 0) Then
           engke = engke+tmp*rgdwgt(0,rgdtyp)*(rgdvxx(irgd)**2+rgdvyy(irgd)**2+rgdvzz(irgd)**2)
           vdotf = vdotf+tmp*(rgdvxx(irgd)*rgdfxx(irgd)+rgdvyy(irgd)*rgdfyy(irgd)+rgdvzz(irgd)*rgdfzz(irgd))
        End If

        engrot= engrot+tmp*(rgdrix(1,rgdtyp)*rgdoxx(irgd)**2+rgdriy(1,rgdtyp)*rgdoyy(irgd)**2+rgdriz(1,rgdtyp)*rgdozz(irgd)**2)
        odott = odott+tmp*(rgdoxx(irgd)*rgdtxx(irgd)+rgdoyy(irgd)*rgdtyy(irgd)+rgdozz(irgd)*rgdtzz(irgd))
     End If
  End Do

  If (mxnode > 1) Then
     buffer(1) = engke
     buffer(2) = vdotf
     buffer(3) = engrot
     buffer(4) = odott
     Call gsum(buffer)
     engke = buffer(1)
     vdotf = buffer(2)
     engrot= buffer(3)
     odott = buffer(4)
  End If

! velocity friction and temperature scaling coefficient
! for Evans thermostat at tstep

  chit = (vdotf+odott)/(engke+engrot)

! get corrected energies

  engke = 0.5_wp*engke
  engrot = 0.5_wp*engrot

  If (isw == 0) Return

! thermostat velocities

  scale=Exp(-chit*tstep)

  Do j=1,nfree
     i=lstfre(j)

     vxx(i)=vxx(i)*scale
     vyy(i)=vyy(i)*scale
     vzz(i)=vzz(i)*scale
  End Do

  Do irgd=1,ntrgd
     rgdvxx(irgd)=rgdvxx(irgd)*scale
     rgdvyy(irgd)=rgdvyy(irgd)*scale
     rgdvzz(irgd)=rgdvzz(irgd)*scale

     rgdoxx(irgd)=rgdoxx(irgd)*scale
     rgdoyy(irgd)=rgdoyy(irgd)*scale
     rgdozz(irgd)=rgdozz(irgd)*scale
  End Do

! thermostat kinetic and rotational energy

  tmp=scale**2
  engke = engke*tmp
  engrot = engrot*tmp

End Subroutine nvt_e1_scl
