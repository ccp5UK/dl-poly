Subroutine rigid_bodies_split_torque(imcon,gxx,gyy,gzz,txx,tyy,tzz,uxx,uyy,uzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for resolving RBs' torques into equivalent atomic
! forces suitable for conjugate gradient minimisation
!
! copyright - daresbury laboratory
! author    - w.smith may 2006
! adapted   - i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode
  Use setup_module,       Only : nrite,mxatms,mxrgd,mxlrgd
  Use config_module,      Only : cell,xxx,yyy,zzz
  Use rigid_bodies_module

  Implicit None

  Integer,           Intent( In    ) :: imcon
  Real( Kind = wp ), Intent( InOut ) :: gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms), &
                                        uxx(1:mxatms),uyy(1:mxatms),uzz(1:mxatms)

  Integer           :: fail,i,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ) :: fmx,fmy,fmz,x(1:1),y(1:1),z(1:1),tqx,tqy,tqz,torque,tmp,qqq

  Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)

  fail = 0
  Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_split_torque allocation failure, node: ', idnode
     Call error(0)
  End If

! Get a RB particles vectors wrt the RB's COM

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

! COM distances

           ggx(krgd)=xxx(i)-rgdxxx(irgd)
           ggy(krgd)=yyy(i)-rgdyyy(irgd)
           ggz(krgd)=zzz(i)-rgdzzz(irgd)
        End Do
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,ggx,ggy,ggz)

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then

! calculate net force and torque on the RB - assumption must hold here

        fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
        tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp

        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

! If the RB has a frozen particle then no net force

           If (rgdfrz(0,rgdtyp) == 0) Then
              fmx=fmx+gxx(i)
              fmy=fmy+gyy(i)
              fmz=fmz+gzz(i)
           End If

           tqx=tqx+ggy(krgd)*gzz(i)-ggz(krgd)*gyy(i)
           tqy=tqy+ggz(krgd)*gxx(i)-ggx(krgd)*gzz(i)
           tqz=tqz+ggx(krgd)*gyy(i)-ggy(krgd)*gxx(i)
        End Do

! offload new forces on RB members

        tmp=1.0_wp/Real(lrgd,wp)
        fmx=fmx*tmp
        fmy=fmy*tmp
        fmz=fmz*tmp

        Do jrgd=1,lrgd
           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (rgdfrz(0,rgdtyp) == 0) Then
              gxx(i)=fmx
              gyy(i)=fmy
              gzz(i)=fmz
           End If

        End Do

! If the RB has 2+ frozen particles (ill=1) the net torque
! must align along the axis of rotation

        If (rgdfrz(0,rgdtyp) > 1) Then
           i1=indrgd(rgdind(1,rgdtyp),irgd)
           i2=indrgd(rgdind(2,rgdtyp),irgd)

           x(1)=xxx(i1)-xxx(i2)
           y(1)=yyy(i1)-yyy(i2)
           z(1)=zzz(i1)-zzz(i2)

           Call images(imcon,cell,1,x,y,z)

           tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)/(x(1)**2+y(1)**2+z(1)**2)
           tqx=x(1)*tmp
           tqy=y(1)*tmp
           tqz=z(1)*tmp
        End If

! Get magnitude of torque

        torque=Sqrt(tqx**2+tqy**2+tqz**2)

! Calculate unit vectors for new site forces

        krgd=krgd-lrgd
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site
           If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
              txx(i)=ggy(krgd)*tqz-tqy*ggz(krgd)
              tyy(i)=ggz(krgd)*tqx-tqz*ggx(krgd)
              tzz(i)=ggx(krgd)*tqy-tqx*ggy(krgd)

              tmp=Sqrt(txx(i)**2+tyy(i)**2+tzz(i)**2)
              If (tmp > 1.0e-10_wp) Then
                 tmp=1.0_wp/tmp
                 txx(i)=txx(i)*tmp
                 tyy(i)=tyy(i)*tmp
                 tzz(i)=tzz(i)*tmp
              Else
                 txx(i)=0.0_wp
                 tyy(i)=0.0_wp
                 tzz(i)=0.0_wp
              End If
           Else
              txx(i)=0.0_wp
              tyy(i)=0.0_wp
              tzz(i)=0.0_wp
           End If
        End Do

! construct unit vectors for new radial positions of
! RB's sites with respect to the axis of rotation

        tmp=1.0_wp
        If (torque > 1.0e-10_wp) tmp=1.0_wp/torque

        Do jrgd=1,lrgd
           If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
              i=indrgd(jrgd,irgd) ! local index of particle/site

              uxx(i)=(tyy(i)*tqz-tqy*tzz(i))*tmp
              uyy(i)=(tzz(i)*tqx-tqz*txx(i))*tmp
              uzz(i)=(txx(i)*tqy-tqx*tyy(i))*tmp
           Else
              uxx(i)=0.0_wp
              uyy(i)=0.0_wp
              uzz(i)=0.0_wp
           End If
        End Do

! scale unit vectors to working lengths

        qqq=0.0_wp

        krgd=krgd-lrgd
        Do jrgd=1,lrgd
           krgd=krgd+1

           If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
              i=indrgd(jrgd,irgd) ! local index of particle/site

              tmp=ggx(krgd)*uxx(i)+ggy(krgd)*uyy(i)+ggz(krgd)*uzz(i)

              txx(i)=txx(i)*tmp
              tyy(i)=tyy(i)*tmp
              tzz(i)=tzz(i)*tmp

              uxx(i)=uxx(i)*tmp
              uyy(i)=uyy(i)*tmp
              uzz(i)=uzz(i)*tmp

              qqq=qqq+tmp**2
           End If
        End Do

! scale new site forces so that all RB members rotate
! to the same angle in their plains of rotation

        tmp=0.0_wp
        If (qqq > 1.0e-10_wp) tmp=torque/qqq

        Do jrgd=1,lrgd
           If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
              i=indrgd(jrgd,irgd) ! local index of particle/site

              txx(i)=txx(i)*tmp
              tyy(i)=tyy(i)*tmp
              tzz(i)=tzz(i)*tmp
           End If
        End Do

     Else

! Initialise unit vectors for new site forces and radial location
! in the plane of rotation wrt the axis of rotation (torque vector)

        Do jrgd=1,lrgd
           i=indrgd(jrgd,irgd) ! local index of particle/site

           txx(i)=0.0_wp ; tyy(i)=0.0_wp ; tzz(i)=0.0_wp
           uxx(i)=0.0_wp ; uyy(i)=0.0_wp ; uzz(i)=0.0_wp
        End Do

     End If
  End Do

  Deallocate (ggx,ggy,ggz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_split_torque deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rigid_bodies_split_torque
