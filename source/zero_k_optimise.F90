Subroutine zero_k_optimise(strkin,strknf,strknt,engke,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for zero Kelvin temperature optimization
!
! in preparation for integration of equations of motion:
!
! the free particle velocity, V, is set in the direction of the force,
! F, when V.F > 0 :: V=F*[(V.F)/(F.F)] and set to zero if V.F < 0 :: V=0.
! the same rational is extended to RB dynamics where, additionally to
! applying the same strategy to the RB COM velocity change upon the COM
! force, alos the angular velocity of the RB, W, is set in the direction
! of the torque, T, for when W.T > 0 :: W=T*[(W.T)/(T.T)] and set to zero
! if W.T < 0 :: W=0.
!
! care must be taken for:
!     - to remove any COM motion generation
!     - to zero angular momentum about centre of mass for non-periodic
!       system (imcon=0)
!     - to ensure preservation of thermostat's instantaneous kinetic
!       energy components and remove spurious temperature fluctuations
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module,   Only : mxrgd,mxlrgd,nrite
  Use configuration,  Only : imcon,cell,natms,nfree,lfrzn,lstfre,weight, &
                             xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use kinetic_module, Only : getcom,getvom,getkin,getknf,getknt,getknr, &
                             kinstress,kinstresf,kinstrest

  Implicit None

  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                        strknf(1:9),strknt(1:9),engrot

  Integer           :: fail,i,j,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: e_f,e_t,e_r,engkf,engkt,scale,vdotf,fsq, &
                       amx,amy,amz,wxx,wyy,wzz,tmp,tmp1,        &
                       com(1:3),vom(1:3),rot(1:9),rotinv(1:9),  &
                       x(1:1),y(1:1),z(1:1),                    &
                       fmx,fmy,fmz,tqx,tqy,tqz,trx,try,trz,     &
                       vpx,vpy,vpz

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer,ggx,ggy,ggz


! preserve magnitudes of old instantaneous energies in order to scale back to them

  e_t=0.5_wp*(strknt(1)+strknt(5)+strknt(9)) ! RBs translational
  e_r=engrot                                 ! RBs rotational
  e_f=engke-e_t                              ! FPs (free particles)

! initialise RB energy components

  engkf=0.0_wp
  engkt=0.0_wp

! recover megrgd

  megrgd=rgdmeg

  If (megrgd > 0) Then
     fail=0
     Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'zero_k_optimise allocation failure, node: ', idnode
        Call error(0)
     End If

! Get the RB particles vectors wrt the RB's COM

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

! Free particles

     Do j=1,nfree
        i=lstfre(j)

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then

! take component of velocity in direction of force

           vdotf = vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)

           If (vdotf < 0.0_wp) Then
              vxx(i) = 0.0_wp
              vyy(i) = 0.0_wp
              vzz(i) = 0.0_wp
           Else
              fsq = fxx(i)**2+fyy(i)**2+fzz(i)**2
              scale = vdotf/Max(1.0e-10_wp,fsq)

              vxx(i) = fxx(i)*scale
              vyy(i) = fyy(i)*scale
              vzz(i) = fzz(i)*scale
           End If

        End If
     End Do

! RBs

     krgd=0
     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then

! current rotation matrix

           Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! calculate COM force and torque

           fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
           tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
           Do jrgd=1,lrgd
              krgd=krgd+1

              i=indrgd(jrgd,irgd) ! local index of particle/site

! If the RB has a frozen particle then no net force

              If (rgdfrz(0,rgdtyp) == 0) Then
                 fmx=fmx+fxx(i)
                 fmy=fmy+fyy(i)
                 fmz=fmz+fzz(i)
              End If

              tqx=tqx+ggy(krgd)*fzz(i)-ggz(krgd)*fyy(i)
              tqy=tqy+ggz(krgd)*fxx(i)-ggx(krgd)*fzz(i)
              tqz=tqz+ggx(krgd)*fyy(i)-ggy(krgd)*fxx(i)
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

! calculate torque in principal frame

           trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
           try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
           trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

           If (rgdfrz(0,rgdtyp) == 0) Then

! take component of velocity in direction of force

              vdotf = rgdvxx(irgd)*fmx+rgdvyy(irgd)*fmy+rgdvzz(irgd)*fmz

              If (vdotf < 0.0_wp) Then
                 rgdvxx(irgd) = 0.0_wp
                 rgdvyy(irgd) = 0.0_wp
                 rgdvzz(irgd) = 0.0_wp
              Else
                 fsq = fmx**2+fmy**2+fmz**2
                 scale = vdotf/Max(1.0e-10_wp,fsq)

                 rgdvxx(irgd) = fmx*scale
                 rgdvyy(irgd) = fmy*scale
                 rgdvzz(irgd) = fmz*scale
              End If

           End If

! take component of the angular velocity in direction of
! the angular acceleration (torque./RI.)

           trx=trx*rgdrix(2,rgdtyp)
           try=try*rgdriy(2,rgdtyp)
           trz=trz*rgdriz(2,rgdtyp)

           vdotf = rgdoxx(irgd)*trx+rgdoyy(irgd)*try+rgdozz(irgd)*trz

           If (vdotf < 0.0_wp) Then
              rgdoxx(irgd) = 0.0_wp
              rgdoyy(irgd) = 0.0_wp
              rgdozz(irgd) = 0.0_wp
           Else
              fsq = trx**2+try**2+trz**2
              scale = vdotf/Max(1.0e-10_wp,fsq)

              rgdoxx(irgd) = trx*scale
              rgdoyy(irgd) = try*scale
              rgdozz(irgd) = trz*scale
           End If

! update RB members velocities

           Do jrgd=1,lrgd
              If (rgdfrz(jrgd,rgdtyp) == 0) Then
                 i=indrgd(jrgd,irgd) ! local index of particle/site

                 If (i <= natms) Then
                    x(1)=rgdx(jrgd,rgdtyp)
                    y(1)=rgdy(jrgd,rgdtyp)
                    z(1)=rgdz(jrgd,rgdtyp)

! new atomic velocities in body frame

                    vpx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                    vpy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                    vpz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

! new atomic velocities in lab frame

                    vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rgdvxx(irgd)
                    vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rgdvyy(irgd)
                    vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rgdvzz(irgd)
                 End If
              End If
           End Do
        End If
     End Do

! Subtract COM velocity

     Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz)

! remove centre of mass motion

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

! update kinetic energy and stress

     Call kinstresf(vxx,vyy,vzz,strknf)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

     strkin=strknf+strknt
     engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! update rotational energy

     engrot=getknr(rgdoxx,rgdoyy,rgdozz)

     Deallocate (ggx,ggy,ggz, Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'zero_k_optimise deallocation failure, node: ', idnode
        Call error(0)
     End If
  Else
     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then

! take component of velocity in direction of force

           vdotf = vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)

           If (vdotf < 0.0_wp) Then
              vxx(i) = 0.0_wp
              vyy(i) = 0.0_wp
              vzz(i) = 0.0_wp
           Else
              fsq = fxx(i)**2+fyy(i)**2+fzz(i)**2
              scale = vdotf/Max(1.0e-10_wp,fsq)

              vxx(i) = fxx(i)*scale
              vyy(i) = fyy(i)*scale
              vzz(i) = fzz(i)*scale
           End If

        End If
     End Do

! Subtract COM velocity

     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do

! Update kinetic stress and energy

     Call kinstress(vxx,vyy,vzz,strkin)
     engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))
  End If

! zero angular momentum about centre of mass - non-periodic system

  If (imcon == 0) Then
     fail=0
     Allocate (buffer(1:12), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'zero_k_optimise allocation failure, node: ', idnode
        Call error(0)
     End If

! initialise RB energy components

     engkf=0.0_wp
     engkt=0.0_wp

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
              x(1)=(wyy*rgdzzz(irgd) - wzz*rgdyyy(irgd))
              y(1)=(wzz*rgdxxx(irgd) - wxx*rgdzzz(irgd))
              z(1)=(wxx*rgdyyy(irgd) - wyy*rgdxxx(irgd))

              rgdvxx(irgd) = rgdvxx(irgd) + x(1)
              rgdvyy(irgd) = rgdvyy(irgd) + y(1)
              rgdvzz(irgd) = rgdvzz(irgd) + z(1)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! local index of particle/site

                 If (i <= natms) Then
                    vxx(i) = vxx(i) + x(1)
                    vyy(i) = vyy(i) + y(1)
                    vzz(i) = vzz(i) + z(1)
                 End If
              End Do
           End If
        End Do

! get kinetic energy

        engkf=getknf(vxx,vyy,vzz)
        engkt=getknt(rgdvxx,rgdvyy,rgdvzz)
        engke=engkf+engkt

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

! get kinetic energy

        engke=getkin(vxx,vyy,vzz)

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
        Write(nrite,'(/,1x,a,i0)') 'zero_k_optimise deallocation failure, node: ', idnode
        Call error(0)
     End If
  End If

! to remove spurious temperature fluctuations ensure preservation
! of thermostat's instantaneous translational and rotational energies
! equipartitioning of DoFs may be lost as transfers of energy between
! translational and rotational DoFs may happen

! apply temperature components scaling

  engkf=engke-engkt
  If (engkf+engkt+engrot > 1.0e-6_wp .and. e_f+e_t+e_r > 1.0e-6_wp) Then
     If (megrgd > 0) Then
        tmp=Sqrt((e_f+e_t+e_r)/(engkf+engkt+engrot))
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
                       x(1)=rgdx(jrgd,rgdtyp)
                       y(1)=rgdy(jrgd,rgdtyp)
                       z(1)=rgdz(jrgd,rgdtyp)

! site velocity in body frame

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

! update kinetic energy and stress

        Call kinstresf(vxx,vyy,vzz,strknf)
        Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

        strkin=strknf+strknt
        engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! update rotational energy

        engrot=getknr(rgdoxx,rgdoyy,rgdozz)
     Else
        tmp=Sqrt(e_f/engke)
        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              vxx(i)=vxx(i)*tmp
              vyy(i)=vyy(i)*tmp
              vzz(i)=vzz(i)*tmp
           End If
        End Do

! Update kinetic stress and energy

        Call kinstress(vxx,vyy,vzz,strkin)
        engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))
     End If
  Else ! zero them and let's see if we can crash
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

End Subroutine zero_k_optimise
