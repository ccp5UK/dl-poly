Subroutine scale_temperature(sigma,degtra,degrot,degfre)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to scale the instantaneous system temperature
! to the target temperature
!
! Note: zeroes angular momentum in non-periodic systems, frozen
! particles are not considered
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use config_module,  Only : imcon,natms,nfree,lfrzn,lstfre, &
                             weight,xxx,yyy,zzz,vxx,vyy,vzz
  Use rigid_bodies_module
  Use kinetic_module, Only : getcom,getvom,getkin,getknf,getknt,getknr

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: sigma
  Integer(Kind=ip),  Intent( In    ) :: degtra,degrot,degfre


  Integer           :: fail,i,j,irgd,jrgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: engke,engkf,engkt,engrot,         &
                       amx,amy,amz,wxx,wyy,wzz,tmp,tmp1, &
                       com(1:3),vom(1:3),rot(1:9),rotinv(1:9),x,y,z

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer



! recover megrgd

  megrgd=rgdmeg

! remove centre of mass motion

  If (megrgd > 0) Then
     Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz)

     Do j=1,nfree
        i=lstfre(j)

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do

     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

        If (rgdfrz(0,rgdtyp) == 0) Then
           rgdvxx(irgd) = rgdvxx(irgd) - vom(1)
           rgdvyy(irgd) = rgdvyy(irgd) - vom(2)
           rgdvzz(irgd) = rgdvzz(irgd) - vom(3)

           lrgd=listrgd(-1,irgd)
           Do jrgd=1,lrgd
              i=indrgd(jrgd,irgd) ! local index of particle/site

              If (i <= natms) Then
                 vxx(i) = vxx(i) - vom(1)
                 vyy(i) = vyy(i) - vom(2)
                 vzz(i) = vzz(i) - vom(3)
              End If
           End Do
        End If
     End Do
  Else
     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do
  End If

! zero angular momentum about centre of mass - non-periodic system

  If (imcon == 0) Then
     fail=0
     Allocate (buffer(1:12), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'scale_temperature allocation failure, node: ', idnode
        Call error(0)
     End If

! calculate centre of mass position

     Call getcom(xxx,yyy,zzz,com)

     If (megrgd > 0) Then

! move to centre of mass origin

        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              xxx(i) = xxx(i) - com(1)
              yyy(i) = yyy(i) - com(2)
              zzz(i) = zzz(i) - com(3)
           End If
        End Do

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           If (rgdfrz(0,rgdtyp) == 0) Then
              rgdxxx(irgd) = rgdxxx(irgd) - com(1)
              rgdyyy(irgd) = rgdyyy(irgd) - com(2)
              rgdzzz(irgd) = rgdzzz(irgd) - com(3)
           End If
        End Do

! angular momentum accumulators

        amx = 0.0_wp
        amy = 0.0_wp
        amz = 0.0_wp

! rotational inertia accumulators

        rot = 0.0_wp

        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              amx = amx + weight(i)*(yyy(i)*vzz(i) - zzz(i)*vyy(i))
              amy = amy + weight(i)*(zzz(i)*vxx(i) - xxx(i)*vzz(i))
              amz = amz + weight(i)*(xxx(i)*vyy(i) - yyy(i)*vxx(i))

              tmp = xxx(i)**2 + yyy(i)**2 + zzz(i)**2
              rot(1) = rot(1) + weight(i)*(xxx(i)*xxx(i) - tmp)
              rot(2) = rot(2) + weight(i)* xxx(i)*yyy(i)
              rot(3) = rot(3) + weight(i)* xxx(i)*zzz(i)
              rot(5) = rot(5) + weight(i)*(yyy(i)*yyy(i) - tmp)
              rot(6) = rot(6) + weight(i)* yyy(i)*zzz(i)
              rot(9) = rot(9) + weight(i)*(zzz(i)*zzz(i) - tmp)
           End If
        End Do

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           If (rgdfrz(0,rgdtyp) == 0) Then
              lrgd=listrgd(-1,irgd)

              tmp1=rgdwgt(0,rgdtyp)*Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

              amx = amx + tmp1*(rgdyyy(irgd)*rgdvzz(irgd) - rgdzzz(irgd)*rgdvyy(irgd))
              amy = amy + tmp1*(rgdzzz(irgd)*rgdvxx(irgd) - rgdxxx(irgd)*rgdvzz(irgd))
              amz = amz + tmp1*(rgdxxx(irgd)*rgdvyy(irgd) - rgdyyy(irgd)*rgdvxx(irgd))

              tmp = rgdxxx(irgd)**2 + rgdyyy(irgd)**2 + rgdzzz(irgd)**2

              rot(1) = rot(1) + tmp1*(rgdxxx(irgd)*rgdxxx(irgd) - tmp)
              rot(2) = rot(2) + tmp1* rgdxxx(irgd)*rgdyyy(irgd)
              rot(3) = rot(3) + tmp1* rgdxxx(irgd)*rgdzzz(irgd)
              rot(5) = rot(5) + tmp1*(rgdyyy(irgd)*rgdyyy(irgd) - tmp)
              rot(6) = rot(6) + tmp1* rgdyyy(irgd)*rgdzzz(irgd)
              rot(9) = rot(9) + tmp1*(rgdzzz(irgd)*rgdzzz(irgd) - tmp)
           End If
        End Do

! complete rotational inertia matrix

        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)

! global sum of rotation

        If (mxnode > 1) Then
           buffer(1) = amx
           buffer(2) = amy
           buffer(3) = amz
           Do i=1,9
              buffer(i+3) = rot(i)
           End Do

           Call gsum(buffer)

           amx =  buffer(1)
           amy =  buffer(2)
           amz =  buffer(3)
           Do i=1,9
              rot(i) = buffer(i+3)
           End Do
        End If

! invert rotational inertia matrix

        Call invert(rot,rotinv,tmp)

! correction to angular velocity

        wxx = rotinv(1)*amx + rotinv(2)*amy + rotinv(3)*amz
        wyy = rotinv(4)*amx + rotinv(5)*amy + rotinv(6)*amz
        wzz = rotinv(7)*amx + rotinv(8)*amy + rotinv(9)*amz

! correction to linear velocity

        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              vxx(i) = vxx(i) + (wyy*zzz(i) - wzz*yyy(i))
              vyy(i) = vyy(i) + (wzz*xxx(i) - wxx*zzz(i))
              vzz(i) = vzz(i) + (wxx*yyy(i) - wyy*xxx(i))
           End If
        End Do

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           If (rgdfrz(0,rgdtyp) == 0) Then
              x=(wyy*rgdzzz(irgd) - wzz*rgdyyy(irgd))
              y=(wzz*rgdxxx(irgd) - wxx*rgdzzz(irgd))
              z=(wxx*rgdyyy(irgd) - wyy*rgdxxx(irgd))

              rgdvxx(irgd) = rgdvxx(irgd) + x
              rgdvyy(irgd) = rgdvyy(irgd) + y
              rgdvzz(irgd) = rgdvzz(irgd) + z

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! local index of particle/site

                 If (i <= natms) Then
                    vxx(i) = vxx(i) + x
                    vyy(i) = vyy(i) + y
                    vzz(i) = vzz(i) + z
                 End If
              End Do
           End If
        End Do

! reset positions to original reference frame

        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              xxx(i) = xxx(i) + com(1)
              yyy(i) = yyy(i) + com(2)
              zzz(i) = zzz(i) + com(3)
           End If
        End Do

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           If (rgdfrz(0,rgdtyp) == 0) Then
              rgdxxx(irgd) = rgdxxx(irgd) + com(1)
              rgdyyy(irgd) = rgdyyy(irgd) + com(2)
              rgdzzz(irgd) = rgdzzz(irgd) + com(3)
           End If
        End Do

     Else

! move to centre of mass origin

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              xxx(i) = xxx(i) - com(1)
              yyy(i) = yyy(i) - com(2)
              zzz(i) = zzz(i) - com(3)
           End If
        End Do

! angular momentum accumulators

        amx = 0.0_wp
        amy = 0.0_wp
        amz = 0.0_wp

! rotational inertia accumulators

        rot = 0.0_wp

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              amx = amx + weight(i)*(yyy(i)*vzz(i) - zzz(i)*vyy(i))
              amy = amy + weight(i)*(zzz(i)*vxx(i) - xxx(i)*vzz(i))
              amz = amz + weight(i)*(xxx(i)*vyy(i) - yyy(i)*vxx(i))

              tmp = xxx(i)**2 + yyy(i)**2 + zzz(i)**2
              rot(1) = rot(1) + weight(i)*(xxx(i)*xxx(i) - tmp)
              rot(2) = rot(2) + weight(i)* xxx(i)*yyy(i)
              rot(3) = rot(3) + weight(i)* xxx(i)*zzz(i)
              rot(5) = rot(5) + weight(i)*(yyy(i)*yyy(i) - tmp)
              rot(6) = rot(6) + weight(i)* yyy(i)*zzz(i)
              rot(9) = rot(9) + weight(i)*(zzz(i)*zzz(i) - tmp)
           End If
        End Do

! complete rotational inertia matrix

        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)

! global sum of rotation

        If (mxnode > 1) Then
           buffer(1) = amx
           buffer(2) = amy
           buffer(3) = amz
           Do i=1,9
              buffer(i+3) = rot(i)
           End Do

           Call gsum(buffer)

           amx =  buffer(1)
           amy =  buffer(2)
           amz =  buffer(3)
           Do i=1,9
              rot(i) = buffer(i+3)
           End Do
        End If

! invert rotational inertia matrix

        Call invert(rot,rotinv,tmp)

! correction to angular velocity

        wxx = rotinv(1)*amx + rotinv(2)*amy + rotinv(3)*amz
        wyy = rotinv(4)*amx + rotinv(5)*amy + rotinv(6)*amz
        wzz = rotinv(7)*amx + rotinv(8)*amy + rotinv(9)*amz

! correction to linear velocity

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              vxx(i) = vxx(i) + (wyy*zzz(i) - wzz*yyy(i))
              vyy(i) = vyy(i) + (wzz*xxx(i) - wxx*zzz(i))
              vzz(i) = vzz(i) + (wxx*yyy(i) - wyy*xxx(i))
           End If
        End Do

! reset positions to original reference frame

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              xxx(i) = xxx(i) + com(1)
              yyy(i) = yyy(i) + com(2)
              zzz(i) = zzz(i) + com(3)
           End If
        End Do

     End If

     Deallocate (buffer, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'scale_temperature deallocation failure, node: ', idnode
        Call error(0)
     End If
  End If

! ensure equipartitioning - all degrees of freedom are equal
! calculate energy: free particles and RB translational and rotational

  engrot=0.0_wp
  If (megrgd > 0) Then
     engkf=getknf(vxx,vyy,vzz)
     engkt=getknt(rgdvxx,rgdvyy,rgdvzz)

     engrot=getknr(rgdoxx,rgdoyy,rgdozz)

! temporary replacement for small engrot

     tmp1=Max(engrot,1.0e-6_wp)

! Scale rotational energy to translational energy
! according to their respective DoF

     If (degtra > Int(0,ip)) Then ! engkt > 0 (degrot > 0 and tmp1(engrot) > 0)
        tmp=Sqrt((engkt*Real(degrot,wp))/(tmp1*Real(degtra,wp)))

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              rgdoxx(irgd)=rgdoxx(irgd)*tmp
              rgdoyy(irgd)=rgdoyy(irgd)*tmp
              rgdozz(irgd)=rgdozz(irgd)*tmp
           End If
        End Do
        engrot=engrot*tmp**2
        tmp1=Max(engrot,1.0e-6_wp)
     End If

! Scale the energy per DoF of the RBs to that of a DoF of the free particles

     If (degfre-degtra-degrot > Int(0,ip)) Then ! engkf > 0
        tmp=Sqrt((engkf*Real(degrot,wp))/(tmp1*Real(degfre-degtra-degrot,wp)))

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              rgdoxx(irgd)=rgdoxx(irgd)*tmp
              rgdoyy(irgd)=rgdoyy(irgd)*tmp
              rgdozz(irgd)=rgdozz(irgd)*tmp

              If (rgdfrz(0,rgdtyp) == 0) Then
                 rgdvxx(irgd)=rgdvxx(irgd)*tmp
                 rgdvyy(irgd)=rgdvyy(irgd)*tmp
                 rgdvzz(irgd)=rgdvzz(irgd)*tmp
              End If
           End If
        End Do
        engrot=engrot*tmp**2
        If (degtra > Int(0,ip)) engkt=engkt*tmp**2
     End If

     engke=engkf+engkt
  Else
     engke=getkin(vxx,vyy,vzz)
  End If

! apply temperature scaling

  If (engke+engrot > 1.0e-6_wp .and. sigma > zero_plus) Then
     tmp=Sqrt(sigma/(engke+engrot))

     If (megrgd > 0) Then
        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              vxx(i)=vxx(i)*tmp
              vyy(i)=vyy(i)*tmp
              vzz(i)=vzz(i)*tmp
           End If
        End Do

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then

! new angular velocity

              rgdoxx(irgd)=rgdoxx(irgd)*tmp
              rgdoyy(irgd)=rgdoyy(irgd)*tmp
              rgdozz(irgd)=rgdozz(irgd)*tmp

! new translational velocity

              If (rgdfrz(0,rgdtyp) == 0) Then
                 rgdvxx(irgd)=rgdvxx(irgd)*tmp
                 rgdvyy(irgd)=rgdvyy(irgd)*tmp
                 rgdvzz(irgd)=rgdvzz(irgd)*tmp
              End If

! new rotational matrix

              Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

              Do jrgd=1,lrgd
                 If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                    i=indrgd(jrgd,irgd) ! local index of particle/site

                    If (i <= natms) Then
                       x=rgdx(jrgd,rgdtyp)
                       y=rgdy(jrgd,rgdtyp)
                       z=rgdz(jrgd,rgdtyp)

! site velocity in body frame

                       wxx=rgdoyy(irgd)*z-rgdozz(irgd)*y
                       wyy=rgdozz(irgd)*x-rgdoxx(irgd)*z
                       wzz=rgdoxx(irgd)*y-rgdoyy(irgd)*x

! new atomic velocities in lab frame

                       vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+rgdvxx(irgd)
                       vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+rgdvyy(irgd)
                       vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+rgdvzz(irgd)
                    End If
                 End If
              End Do
           End If
        End Do
     Else
        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              vxx(i)=vxx(i)*tmp
              vyy(i)=vyy(i)*tmp
              vzz(i)=vzz(i)*tmp
           End If
        End Do
     End If
  Else ! sigma must be zero
     Do i=1,natms
        vxx(i)=0.0_wp ; vyy(i)=0.0_wp ; vzz(i)=0.0_wp
     End Do

     If (megrgd > 0) Then
        Do irgd=1,ntrgd
           rgdvxx(irgd)=0.0_wp ; rgdvyy(irgd)=0.0_wp ; rgdvzz(irgd)=0.0_wp
           rgdoxx(irgd)=0.0_wp ; rgdoyy(irgd)=0.0_wp ; rgdozz(irgd)=0.0_wp
        End Do
     End If
  End If

End Subroutine scale_temperature
