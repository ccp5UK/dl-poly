Subroutine rigid_bodies_quench()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to convert atomic velocities to RB COM and
! angular velocity
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : imcon,cell,natms,nlast,lsi,lsa,xxx,yyy,zzz,vxx,vyy,vzz
  Use rigid_bodies_module

  Implicit None

  Integer           :: fail,i,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ) :: weight,rot(1:9),wxx,wyy,wzz,x(1:1),y(1:1),z(1:1),tmp

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_quench allocation failure, node: ', idnode
     Call error(0)
  End If

! Halo velocity field across onto neighbouring domains

  If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz)

! translate atomic velocities to COM velocity & angular velocity
! frozen velocities and massless sites weights are assumed zero!!!

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     rgdvxx(irgd)=0.0_wp ; rgdvyy(irgd)=0.0_wp ; rgdvzz(irgd)=0.0_wp

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters
        Do jrgd=1,lrgd
           krgd=krgd+1

! local index and mass of particle/site

           i=indrgd(jrgd,irgd)
           weight=rgdwgt(jrgd,rgdtyp)

! COM distances

           gxx(krgd)=xxx(i)-rgdxxx(irgd)
           gyy(krgd)=yyy(i)-rgdyyy(irgd)
           gzz(krgd)=zzz(i)-rgdzzz(irgd)

! If the RB has a frozen particle then no net COM momentum

           If (rgdfrz(0,rgdtyp) == 0) Then
              rgdvxx(irgd)=rgdvxx(irgd)+weight*vxx(i)
              rgdvyy(irgd)=rgdvyy(irgd)+weight*vyy(i)
              rgdvzz(irgd)=rgdvzz(irgd)+weight*vzz(i)
           End If
        End Do

! COM velocity
! If the RB has a frozen particle then no net COM momentum

        If (rgdfrz(0,rgdtyp) == 0) Then
           tmp=1.0_wp/rgdwgt(0,rgdtyp)
           rgdvxx(irgd)=rgdvxx(irgd)*tmp
           rgdvyy(irgd)=rgdvyy(irgd)*tmp
           rgdvzz(irgd)=rgdvzz(irgd)*tmp
        End If
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Get RBs' angular momenta

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     rgdoxx(irgd)=0.0_wp ; rgdoyy(irgd)=0.0_wp; rgdozz(irgd)=0.0_wp

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

! new rotational matrix

        Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! angular momentum accumulators

        wxx=0.0_wp
        wyy=0.0_wp
        wzz=0.0_wp

        Do jrgd=1,lrgd
           krgd=krgd+1

! local index and mass of particle/site - assumption must hold here

           i=indrgd(jrgd,irgd)
           weight=rgdwgt(jrgd,rgdtyp)

           wxx=wxx+weight*(gyy(krgd)*vzz(i)-gzz(krgd)*vyy(i))
           wyy=wyy+weight*(gzz(krgd)*vxx(i)-gxx(krgd)*vzz(i))
           wzz=wzz+weight*(gxx(krgd)*vyy(i)-gyy(krgd)*vxx(i))
        End Do

! If the RB has 2+ frozen particles (ill=1) the net angular momentum
! must align along the axis of rotation keeping its magnitude

        If (rgdfrz(0,rgdtyp) > 1) Then
           i1=indrgd(rgdind(1,rgdtyp),irgd)
           i2=indrgd(rgdind(2,rgdtyp),irgd)

           x(1)=xxx(i1)-xxx(i2)
           y(1)=yyy(i1)-yyy(i2)
           z(1)=zzz(i1)-zzz(i2)

           Call images(imcon,cell,1,x,y,z)

           tmp=Sign( Sqrt((wxx**2+wyy**2+wzz**2)/(x(1)**2+y(1)**2+z(1)**2)) , &
                     Nearest((x(1)*wxx+y(1)*wyy+z(1)*wzz),-1.0_wp) )

           wxx=x(1)*tmp
           wyy=y(1)*tmp
           wzz=z(1)*tmp
        End If

! angular velocity in body fixed frame

        rgdoxx(irgd)=(rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rgdrix(2,rgdtyp)
        rgdoyy(irgd)=(rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rgdriy(2,rgdtyp)
        rgdozz(irgd)=(rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rgdriz(2,rgdtyp)

        Do jrgd=1,lrgd
           If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
              i=indrgd(jrgd,irgd) ! local index of particle/site

              If (i <= natms) Then

! site velocity in body frame

                 x(1)=rgdx(jrgd,rgdtyp)
                 y(1)=rgdy(jrgd,rgdtyp)
                 z(1)=rgdz(jrgd,rgdtyp)

                 wxx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                 wyy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                 wzz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

! new atomic velocities in lab frame

                 vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+rgdvxx(irgd)
                 vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+rgdvyy(irgd)
                 vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+rgdvzz(irgd)

              End If
           End If
        End Do
     End If
  End Do

  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_quench deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rigid_bodies_quench

Subroutine rigid_bodies_q_ench(qr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to convert atomic velocities to RB COM and
! angular velocity
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : imcon,cell,natms,nlast,lsi,lsa,xxx,yyy,zzz,vxx,vyy,vzz
  Use rigid_bodies_module

  Implicit None

  Integer, Intent( In    ) :: qr(1:mxrgd)

  Integer           :: fail,i,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ) :: weight,rot(1:9),wxx,wyy,wzz,x(1:1),y(1:1),z(1:1),tmp

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_q_ench allocation failure, node: ', idnode
     Call error(0)
  End If

! Halo velocity field across onto neighbouring domains

  If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz)

! translate atomic velocities to COM velocity & angular velocity
! frozen velocities and massless sites weights are assumed zero!!!

  krgd=0
  Do irgd=1,ntrgd
     If (qr(irgd) == 1) Then
        rgdtyp=listrgd(0,irgd)

        rgdvxx(irgd)=0.0_wp ; rgdvyy(irgd)=0.0_wp ; rgdvzz(irgd)=0.0_wp

! For all good RBs

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters
           Do jrgd=1,lrgd
              krgd=krgd+1

! local index and mass of particle/site

              i=indrgd(jrgd,irgd)
              weight=rgdwgt(jrgd,rgdtyp)

! COM distances

              gxx(krgd)=xxx(i)-rgdxxx(irgd)
              gyy(krgd)=yyy(i)-rgdyyy(irgd)
              gzz(krgd)=zzz(i)-rgdzzz(irgd)

! If the RB has a frozen particle then no net COM momentum

              If (rgdfrz(0,rgdtyp) == 0) Then
                 rgdvxx(irgd)=rgdvxx(irgd)+weight*vxx(i)
                 rgdvyy(irgd)=rgdvyy(irgd)+weight*vyy(i)
                 rgdvzz(irgd)=rgdvzz(irgd)+weight*vzz(i)
              End If
           End Do

! COM velocity
! If the RB has a frozen particle then no net COM momentum

           If (rgdfrz(0,rgdtyp) == 0) Then
              tmp=1.0_wp/rgdwgt(0,rgdtyp)
              rgdvxx(irgd)=rgdvxx(irgd)*tmp
              rgdvyy(irgd)=rgdvyy(irgd)*tmp
              rgdvzz(irgd)=rgdvzz(irgd)*tmp
           End If
        End If
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Get RBs' angular momenta

  krgd=0
  Do irgd=1,ntrgd
     If (qr(irgd) == 1) Then
        rgdtyp=listrgd(0,irgd)

        rgdoxx(irgd)=0.0_wp ; rgdoyy(irgd)=0.0_wp; rgdozz(irgd)=0.0_wp

! For all good RBs

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

! new rotational matrix

           Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! angular momentum accumulators

           wxx=0.0_wp
           wyy=0.0_wp
           wzz=0.0_wp

           Do jrgd=1,lrgd
              krgd=krgd+1

! local index and mass of particle/site - assumption must hold here

              i=indrgd(jrgd,irgd)
              weight=rgdwgt(jrgd,rgdtyp)

              wxx=wxx+weight*(gyy(krgd)*vzz(i)-gzz(krgd)*vyy(i))
              wyy=wyy+weight*(gzz(krgd)*vxx(i)-gxx(krgd)*vzz(i))
              wzz=wzz+weight*(gxx(krgd)*vyy(i)-gyy(krgd)*vxx(i))
           End Do

! If the RB has 2+ frozen particles (ill=1) the net angular momentum
! must align along the axis of rotation keeping its magnitude

           If (rgdfrz(0,rgdtyp) > 1) Then
              i1=indrgd(rgdind(1,rgdtyp),irgd)
              i2=indrgd(rgdind(2,rgdtyp),irgd)

              x(1)=xxx(i1)-xxx(i2)
              y(1)=yyy(i1)-yyy(i2)
              z(1)=zzz(i1)-zzz(i2)

              Call images(imcon,cell,1,x,y,z)

              tmp=Sign( Sqrt((wxx**2+wyy**2+wzz**2)/(x(1)**2+y(1)**2+z(1)**2)) , &
                        Nearest((x(1)*wxx+y(1)*wyy+z(1)*wzz),-1.0_wp) )

              wxx=x(1)*tmp
              wyy=y(1)*tmp
              wzz=z(1)*tmp
           End If

! angular velocity in body fixed frame

           rgdoxx(irgd)=(rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rgdrix(2,rgdtyp)
           rgdoyy(irgd)=(rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rgdriy(2,rgdtyp)
           rgdozz(irgd)=(rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rgdriz(2,rgdtyp)

           Do jrgd=1,lrgd
              If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                 i=indrgd(jrgd,irgd) ! local index of particle/site

                 If (i <= natms) Then

! site velocity in body frame

                    x(1)=rgdx(jrgd,rgdtyp)
                    y(1)=rgdy(jrgd,rgdtyp)
                    z(1)=rgdz(jrgd,rgdtyp)

                    wxx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                    wyy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                    wzz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

! new atomic velocities in lab frame

                    vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+rgdvxx(irgd)
                    vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+rgdvyy(irgd)
                    vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+rgdvzz(irgd)

                 End If
              End If
           End Do
        End If
     End If
  End Do

  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_q_ench deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rigid_bodies_q_ench
