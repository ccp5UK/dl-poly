Subroutine nvt_b1_scl &
           (isw,tstep,sigma,taut,vxx,vyy,vzz,         &
           rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz, &
           chit,strkin,strknf,strknt,engke,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NVT Berendsen thermostat
! when singled rigid bodies are present
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module,        Only : mxatms,mxrgd
  Use config_module,       Only : nfree,lstfre
  Use kinetic_module,      Only : kinstresf,kinstrest,getknr
  Use rigid_bodies_module, Only : ntrgd

  Implicit None

  Integer,                                  Intent( In    ) :: isw
  Real( Kind = wp ),                        Intent( In    ) :: tstep,sigma,taut
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
  Real( Kind = wp ), Dimension( 1:mxrgd  ), Intent( InOut ) :: rgdvxx,rgdvyy,rgdvzz, &
                                                               rgdoxx,rgdoyy,rgdozz
  Real( Kind = wp ), Dimension( 1:9 ),      Intent(   Out ) :: strkin,strknf,strknt
  Real( Kind = wp ),                        Intent(   Out ) :: chit,engke,engrot

  Integer           :: i,j,irgd
  Real( Kind = wp ) :: tmp

! update kinetic energy and stress

  Call kinstresf(vxx,vyy,vzz,strknf)
  Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

  strkin=strknf+strknt
  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! update rotational energy

  engrot=getknr(rgdoxx,rgdoyy,rgdozz)

! temperature scaling coefficient - taut is the decay constant

  chit=Sqrt(1.0_wp+tstep/taut*(sigma/(engke+engrot)-1.0_wp))

  If (isw == 0) Return

! thermostat velocities

  Do j=1,nfree
     i=lstfre(j)

     vxx(i)=vxx(i)*chit
     vyy(i)=vyy(i)*chit
     vzz(i)=vzz(i)*chit
  End Do

  Do irgd=1,ntrgd
     rgdvxx(irgd)=rgdvxx(irgd)*chit
     rgdvyy(irgd)=rgdvyy(irgd)*chit
     rgdvzz(irgd)=rgdvzz(irgd)*chit

     rgdoxx(irgd)=rgdoxx(irgd)*chit
     rgdoyy(irgd)=rgdoyy(irgd)*chit
     rgdozz(irgd)=rgdozz(irgd)*chit
  End Do

! thermostat kinetic stress and energy

  tmp=chit**2

  strknf=strknf*tmp
  strknt=strknt*tmp
  strkin=strknf+strknt

  engke=engke*tmp
  engrot=engrot*tmp

End Subroutine nvt_b1_scl
