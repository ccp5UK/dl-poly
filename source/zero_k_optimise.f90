Subroutine zero_k_optimise(strkin,strknf,strknt,engke,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for zero Kelvin temperature optimization
!
! if V.F < 0 then velocity is set to zero else velocity is set in the
! direction of the force as F*[(V.F)/(F.F)] in preparation for
! integration of equations of motion
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode
  Use setup_module,        Only : mxrgd,mxlrgd,nrite
  Use config_module,       Only : cell,natms,nfree,lfrzn,lstfre,weight, &
                                  xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use kinetic_module,      Only : getvom,getknr,kinstress,kinstresf,kinstrest

  Implicit None

  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                        strknf(1:9),strknt(1:9),engrot

  Integer           :: fail,i,j,i1,i2,      &
                       irgd,jrgd,krgd,lrgd, &
                       rgdtyp,imcon,megrgd
  Real( Kind = wp ) :: scale,vdotf,fsq,vom(1:3),tmp,        &
                       x(1:1),y(1:1),z(1:1),rot(1:9),       &
                       fmx,fmy,fmz,tqx,tqy,tqz,trx,try,trz, &
                       vpx,vpy,vpz

  Real( Kind = wp ), Dimension( : ), Allocatable :: ggx,ggy,ggz


! recover megrgd

  megrgd=rgdmeg

  If (megrgd > 0) Then

     fail=0
     Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'zero_k_optimise allocation failure, node: ', idnode
        Call error(0)
     End If

! recover imcon

     imcon=rgdimc

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

End Subroutine zero_k_optimise
