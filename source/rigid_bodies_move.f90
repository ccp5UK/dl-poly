Subroutine rigid_bodies_move(stride,oxx,oyy,ozz,txx,tyy,tzz,uxx,uyy,uzz,dist_tol)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for updating positions of atoms in RBs
! during a conjugate gradient minimisation
!
! copyright - daresbury laboratory
! author    - w.smith may 2006
! adapted   - i.t.todorov october 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,       Only : mxatms
  Use config_module,      Only : natms,xxx,yyy,zzz
  Use rigid_bodies_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: stride,                                    &
                                        oxx(1:mxatms),oyy(1:mxatms),ozz(1:mxatms), &
                                        txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms), &
                                        uxx(1:mxatms),uyy(1:mxatms),uzz(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: dist_tol

  Integer           :: i,irgd,jrgd,lrgd,rgdtyp
  Real( Kind = wp ) :: ttt,uuu,the,dtol,coz,zin,x,y,z

  dtol=0.0_wp
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
              i=indrgd(jrgd,irgd) ! local index of particle/site

              If (i <= natms) Then

! translational motion

                 If (rgdfrz(0,rgdtyp) == 0) Then
                    x=stride*oxx(i)
                    y=stride*oyy(i)
                    z=stride*ozz(i)
                 Else
                    x=0.0_wp
                    y=0.0_wp
                    z=0.0_wp
                 End If

! calculate this RB's particle distance to the axis of rotation

                 uuu=Sqrt(uxx(i)**2+uyy(i)**2+uzz(i)**2)

! add rotational motion for big enough distances

                 If (uuu > 1.0e-10_wp) Then

! force magnitude

                    ttt=Sqrt(txx(i)**2+tyy(i)**2+tzz(i)**2)

! angle of rotation

                    the=(ttt/uuu)*stride

                    coz=Cos(the)-1.0_wp
                    zin=(Sin(the)/the)*stride

                    x=x+coz*uxx(i)+zin*txx(i)
                    y=y+coz*uyy(i)+zin*tyy(i)
                    z=z+coz*uzz(i)+zin*tzz(i)
                 End If

! Add motion

                 xxx(i)=xxx(i)+x
                 yyy(i)=yyy(i)+y
                 zzz(i)=zzz(i)+z

                 dtol=Max(dtol,x**2+y**2+z**2)
              End If

           End If
        End Do
     End If
  End Do

  dist_tol=Max(dist_tol,Sqrt(dtol))

End Subroutine rigid_bodies_move
