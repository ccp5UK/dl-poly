Subroutine xscale(m_rgd,keyens,tstep,eta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to scale initial positions with change in box shape
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : idnode
  Use setup_module
  Use config_module,       Only : imcon,cell,natms,nlast,nfree,lstfre,lsi,lsa
  Use vnl_module,          Only : llvnl,l_vnl,xbg,ybg,zbg
  Use statistics_module,   Only : xin,yin,zin
  Use rigid_bodies_module, Only : ntrgd,listrgd,indrgd,lshmv_rgd,lishp_rgd,lashp_rgd
  Use kinetic_module,      Only : getcom

  Implicit None

  Integer,           Intent( In    ) :: m_rgd,keyens
  Real( Kind = wp ), Intent( In    ) :: tstep,eta(1:9)

  Integer           :: fail,i,j,irgd,jrgd,lrgd
  Real( Kind = wp ) :: a1,a2,a3,a5,a6,a9,b1,b2,b3,b5,b6,b9,scale, &
                       xa,ya,za,x,y,z,com(1:3)

  Real( Kind = wp ), Allocatable :: rgdxin(:),rgdyin(:),rgdzin(:)

  If (keyens < 20) Return

  If (m_rgd == 0) Then

     If (keyens == 21 .or. keyens == 31) Then

! berendsen npt/nst

        If (keyens == 21) Then

           scale = eta(1)

           Do i=1,natms
              xin(i) = scale*xin(i)
              yin(i) = scale*yin(i)
              zin(i) = scale*zin(i)
           End Do

        Else

           Do i=1,natms
              xa = xin(i)*eta(1)+yin(i)*eta(2)+zin(i)*eta(3)
              ya = xin(i)*eta(4)+yin(i)*eta(5)+zin(i)*eta(6)
              za = xin(i)*eta(7)+yin(i)*eta(8)+zin(i)*eta(9)

              xin(i) = xa
              yin(i) = ya
              zin(i) = za
           End Do

        End If

     Else If (keyens == 22 .or. keyens == 32) Then

! hoover npt/nst

        Call getcom(xin,yin,zin,com)

        If (keyens == 22) Then

           scale = Exp(tstep*eta(1))

           Do i=1,natms
              xin(i) = scale*(xin(i)-com(1))+com(1)
              yin(i) = scale*(yin(i)-com(2))+com(2)
              zin(i) = scale*(zin(i)-com(3))+com(3)
           End Do

        Else

! second order taylor expansion of Exp(tstep*eta)

           a1 = tstep*eta(1)
           a2 = tstep*eta(2)
           a3 = tstep*eta(3)
           a5 = tstep*eta(5)
           a6 = tstep*eta(6)
           a9 = tstep*eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do i=1,natms
              xa = xin(i)-com(1)
              ya = yin(i)-com(2)
              za = zin(i)-com(3)

              xin(i) = xa*b1 + ya*b2 + za*b3 + com(1)
              yin(i) = xa*b2 + ya*b5 + za*b6 + com(2)
              zin(i) = xa*b3 + ya*b6 + za*b9 + com(3)
           End Do

        End If

     Else If (keyens == 20 .or. keyens == 30 .or. &
              keyens == 23 .or. keyens == 33) Then

! Langevin and MTK npt/nst

        If (keyens == 20 .or. keyens == 23) Then

           scale = Exp(tstep*eta(1))

           Do i=1,natms
              xin(i) = scale*xin(i)
              yin(i) = scale*yin(i)
              zin(i) = scale*zin(i)
           End Do

        Else

! second order taylor expansion of Exp(tstep*eta)

           a1 = tstep*eta(1)
           a2 = tstep*eta(2)
           a3 = tstep*eta(3)
           a5 = tstep*eta(5)
           a6 = tstep*eta(6)
           a9 = tstep*eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do i=1,natms
              xa = xin(i)
              ya = yin(i)
              za = zin(i)

              xin(i) = xa*b1 + ya*b2 + za*b3
              yin(i) = xa*b2 + ya*b5 + za*b6
              zin(i) = xa*b3 + ya*b6 + za*b9
           End Do

        End If

     End If

     If (.not.l_vnl) Then

        If (keyens == 21 .or. keyens == 31) Then

! berendsen npt/nst

           If (keyens == 21) Then

              scale = eta(1)

              Do i=1,natms
                 xbg(i) = scale*xbg(i)
                 ybg(i) = scale*ybg(i)
                 zbg(i) = scale*zbg(i)
              End Do

           Else

              Do i=1,natms
                 xa = xbg(i)*eta(1)+ybg(i)*eta(2)+zbg(i)*eta(3)
                 ya = xbg(i)*eta(4)+ybg(i)*eta(5)+zbg(i)*eta(6)
                 za = xbg(i)*eta(7)+ybg(i)*eta(8)+zbg(i)*eta(9)

                 xbg(i) = xa
                 ybg(i) = ya
                 zbg(i) = za
              End Do

           End If

        Else If (keyens == 22 .or. keyens == 32) Then

! hoover npt/nst

           Call getcom(xbg,ybg,zbg,com)

           If (keyens == 22) Then

              scale = Exp(tstep*eta(1))

              Do i=1,natms
                 xbg(i) = scale*(xbg(i)-com(1))+com(1)
                 ybg(i) = scale*(ybg(i)-com(2))+com(2)
                 zbg(i) = scale*(zbg(i)-com(3))+com(3)
              End Do

           Else

! second order taylor expansion of Exp(tstep*eta)

              a1 = tstep*eta(1)
              a2 = tstep*eta(2)
              a3 = tstep*eta(3)
              a5 = tstep*eta(5)
              a6 = tstep*eta(6)
              a9 = tstep*eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do i=1,natms
                 xa = xbg(i)-com(1)
                 ya = ybg(i)-com(2)
                 za = zbg(i)-com(3)

                 xbg(i) = xa*b1 + ya*b2 + za*b3 + com(1)
                 ybg(i) = xa*b2 + ya*b5 + za*b6 + com(2)
                 zbg(i) = xa*b3 + ya*b6 + za*b9 + com(3)
              End Do

           End If

        Else If (keyens == 20 .or. keyens == 30 .or. &
                 keyens == 23 .or. keyens == 33) Then

! Langevin and MTK npt/nst

           If (keyens == 20 .or. keyens == 23) Then

              scale = Exp(tstep*eta(1))

              Do i=1,natms
                 xbg(i) = scale*xbg(i)
                 ybg(i) = scale*ybg(i)
                 zbg(i) = scale*zbg(i)
              End Do

           Else

! second order taylor expansion of Exp(tstep*eta)

              a1 = tstep*eta(1)
              a2 = tstep*eta(2)
              a3 = tstep*eta(3)
              a5 = tstep*eta(5)
              a6 = tstep*eta(6)
              a9 = tstep*eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do i=1,natms
                 xa = xbg(i)
                 ya = ybg(i)
                 za = zbg(i)

                 xbg(i) = xa*b1 + ya*b2 + za*b3
                 ybg(i) = xa*b2 + ya*b5 + za*b6
                 zbg(i) = xa*b3 + ya*b6 + za*b9
              End Do

           End If

        End If

     End If

  Else ! RBs exist

     fail = 0
     Allocate (rgdxin(1:mxrgd),rgdyin(1:mxrgd),rgdzin(1:mxrgd), Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'xscale allocation failure, node: ', idnode
        Call error(0)
     End If

! Halo initial RB members positions across onto neighbouring domains
! to get initial COMs

     If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,xin,yin,zin)
     Call rigid_bodies_coms(xin,yin,zin,rgdxin,rgdyin,rgdzin)

     If (keyens == 21 .or. keyens == 31) Then

! berendsen npt/nst

        If (keyens == 21) Then

           scale = eta(1)

           Do j=1,nfree
              i=lstfre(j)

              xin(i) = scale*xin(i)
              yin(i) = scale*yin(i)
              zin(i) = scale*zin(i)
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              rgdxin(irgd) = scale*rgdxin(irgd)
              rgdyin(irgd) = scale*rgdyin(irgd)
              rgdzin(irgd) = scale*rgdzin(irgd)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    xin(i) = xin(i) - x + rgdxin(irgd)
                    yin(i) = yin(i) - y + rgdyin(irgd)
                    zin(i) = zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        Else

           Do j=1,nfree
              i=lstfre(j)

              xa = xin(i)*eta(1)+yin(i)*eta(2)+zin(i)*eta(3)
              ya = xin(i)*eta(4)+yin(i)*eta(5)+zin(i)*eta(6)
              za = xin(i)*eta(7)+yin(i)*eta(8)+zin(i)*eta(9)

              xin(i) = xa
              yin(i) = ya
              zin(i) = za
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd)*eta(1)+rgdyin(irgd)*eta(2)+rgdzin(irgd)*eta(3)
              ya = rgdxin(irgd)*eta(4)+rgdyin(irgd)*eta(5)+rgdzin(irgd)*eta(6)
              za = rgdxin(irgd)*eta(7)+rgdyin(irgd)*eta(8)+rgdzin(irgd)*eta(9)

              rgdxin(irgd) = xa
              rgdyin(irgd) = ya
              rgdzin(irgd) = za

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    xin(i) = xin(i) - x + rgdxin(irgd)
                    yin(i) = yin(i) - y + rgdyin(irgd)
                    zin(i) = zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        End If

     Else If (keyens == 22 .or. keyens == 32) Then

! hoover npt/nst

        Call getcom(xin,yin,zin,com)

        If (keyens == 22) Then

           scale = Exp(tstep*eta(1))

           Do j=1,nfree
              i=lstfre(j)

              xin(i) = scale*(xin(i)-com(1))+com(1)
              yin(i) = scale*(yin(i)-com(2))+com(2)
              zin(i) = scale*(zin(i)-com(3))+com(3)
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              rgdxin(irgd) = scale*(rgdxin(irgd)-com(1))+com(1)
              rgdyin(irgd) = scale*(rgdyin(irgd)-com(2))+com(2)
              rgdzin(irgd) = scale*(rgdzin(irgd)-com(3))+com(3)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    xin(i) = xin(i) - x + rgdxin(irgd)
                    yin(i) = yin(i) - y + rgdyin(irgd)
                    zin(i) = zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        Else

! second order taylor expansion of Exp(tstep*eta)

           a1 = tstep*eta(1)
           a2 = tstep*eta(2)
           a3 = tstep*eta(3)
           a5 = tstep*eta(5)
           a6 = tstep*eta(6)
           a9 = tstep*eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do j=1,nfree
              i=lstfre(j)

              xa = xin(i)-com(1)
              ya = yin(i)-com(2)
              za = zin(i)-com(3)

              xin(i) = xa*b1 + ya*b2 + za*b3 + com(1)
              yin(i) = xa*b2 + ya*b5 + za*b6 + com(2)
              zin(i) = xa*b3 + ya*b6 + za*b9 + com(3)
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd)-com(1)
              ya = rgdyin(irgd)-com(2)
              za = rgdzin(irgd)-com(3)

              rgdxin(irgd) = xa*b1 + ya*b2 + za*b3 + com(1)
              rgdyin(irgd) = xa*b2 + ya*b5 + za*b6 + com(2)
              rgdzin(irgd) = xa*b3 + ya*b6 + za*b9 + com(3)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    xin(i) = xin(i) - x + rgdxin(irgd)
                    yin(i) = yin(i) - y + rgdyin(irgd)
                    zin(i) = zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        End If

     Else If (keyens == 20 .or. keyens == 30 .or. &
              keyens == 23 .or. keyens == 33) Then

! Langevin and MTK npt/nst

        If (keyens == 20 .or. keyens == 23) Then

           scale = Exp(tstep*eta(1))

           Do j=1,nfree
              i=lstfre(j)

              xin(i) = scale*xin(i)
              yin(i) = scale*yin(i)
              zin(i) = scale*zin(i)
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              rgdxin(irgd) = scale*rgdxin(irgd)
              rgdyin(irgd) = scale*rgdyin(irgd)
              rgdzin(irgd) = scale*rgdzin(irgd)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    xin(i) = xin(i) - x + rgdxin(irgd)
                    yin(i) = yin(i) - y + rgdyin(irgd)
                    zin(i) = zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        Else

! second order taylor expansion of Exp(tstep*eta)

           a1 = tstep*eta(1)
           a2 = tstep*eta(2)
           a3 = tstep*eta(3)
           a5 = tstep*eta(5)
           a6 = tstep*eta(6)
           a9 = tstep*eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do j=1,nfree
              i=lstfre(j)

              xa = xin(i)
              ya = yin(i)
              za = zin(i)

              xin(i) = xa*b1 + ya*b2 + za*b3
              yin(i) = xa*b2 + ya*b5 + za*b6
              zin(i) = xa*b3 + ya*b6 + za*b9
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd)
              ya = rgdyin(irgd)
              za = rgdzin(irgd)

              rgdxin(irgd) = xa*b1 + ya*b2 + za*b3
              rgdyin(irgd) = xa*b2 + ya*b5 + za*b6
              rgdzin(irgd) = xa*b3 + ya*b6 + za*b9

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    xin(i) = xin(i) - x + rgdxin(irgd)
                    yin(i) = yin(i) - y + rgdyin(irgd)
                    zin(i) = zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        End If

     End If

     If (.not.l_vnl) Then

        Call rigid_bodies_coms(xbg,ybg,zbg,rgdxin,rgdyin,rgdzin)

        If (keyens == 21 .or. keyens == 31) Then

! berendsen npt/nst

           If (keyens == 21) Then

              scale = eta(1)

              Do j=1,nfree
                 i=lstfre(j)

                 xbg(i) = scale*xbg(i)
                 ybg(i) = scale*ybg(i)
                 zbg(i) = scale*zbg(i)
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 rgdxin(irgd) = scale*rgdxin(irgd)
                 rgdyin(irgd) = scale*rgdyin(irgd)
                 rgdzin(irgd) = scale*rgdzin(irgd)

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       xbg(i) = xbg(i) - x + rgdxin(irgd)
                       ybg(i) = ybg(i) - y + rgdyin(irgd)
                       zbg(i) = zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           Else

              Do j=1,nfree
                 i=lstfre(j)

                 xa = xbg(i)*eta(1)+ybg(i)*eta(2)+zbg(i)*eta(3)
                 ya = xbg(i)*eta(4)+ybg(i)*eta(5)+zbg(i)*eta(6)
                 za = xbg(i)*eta(7)+ybg(i)*eta(8)+zbg(i)*eta(9)

                 xbg(i) = xa
                 ybg(i) = ya
                 zbg(i) = za
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 xa = rgdxin(irgd)*eta(1)+rgdyin(irgd)*eta(2)+rgdzin(irgd)*eta(3)
                 ya = rgdxin(irgd)*eta(4)+rgdyin(irgd)*eta(5)+rgdzin(irgd)*eta(6)
                 za = rgdxin(irgd)*eta(7)+rgdyin(irgd)*eta(8)+rgdzin(irgd)*eta(9)

                 rgdxin(irgd) = xa
                 rgdyin(irgd) = ya
                 rgdzin(irgd) = za

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       xbg(i) = xbg(i) - x + rgdxin(irgd)
                       ybg(i) = ybg(i) - y + rgdyin(irgd)
                       zbg(i) = zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           End If

        Else If (keyens == 22 .or. keyens == 32) Then

! hoover npt/nst

           Call getcom(xbg,ybg,zbg,com)

           If (keyens == 22) Then

              scale = Exp(tstep*eta(1))

              Do j=1,nfree
                 i=lstfre(j)

                 xbg(i) = scale*(xbg(i)-com(1))+com(1)
                 ybg(i) = scale*(ybg(i)-com(2))+com(2)
                 zbg(i) = scale*(zbg(i)-com(3))+com(3)
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 rgdxin(irgd) = scale*(rgdxin(irgd)-com(1))+com(1)
                 rgdyin(irgd) = scale*(rgdyin(irgd)-com(2))+com(2)
                 rgdzin(irgd) = scale*(rgdzin(irgd)-com(3))+com(3)

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       xbg(i) = xbg(i) - x + rgdxin(irgd)
                       ybg(i) = ybg(i) - y + rgdyin(irgd)
                       zbg(i) = zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           Else

! second order taylor expansion of Exp(tstep*eta)

              a1 = tstep*eta(1)
              a2 = tstep*eta(2)
              a3 = tstep*eta(3)
              a5 = tstep*eta(5)
              a6 = tstep*eta(6)
              a9 = tstep*eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do j=1,nfree
                 i=lstfre(j)

                 xa = xbg(i)-com(1)
                 ya = ybg(i)-com(2)
                 za = zbg(i)-com(3)

                 xbg(i) = xa*b1 + ya*b2 + za*b3 + com(1)
                 ybg(i) = xa*b2 + ya*b5 + za*b6 + com(2)
                 zbg(i) = xa*b3 + ya*b6 + za*b9 + com(3)
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 xa = rgdxin(irgd)-com(1)
                 ya = rgdyin(irgd)-com(2)
                 za = rgdzin(irgd)-com(3)

                 rgdxin(irgd) = xa*b1 + ya*b2 + za*b3 + com(1)
                 rgdyin(irgd) = xa*b2 + ya*b5 + za*b6 + com(2)
                 rgdzin(irgd) = xa*b3 + ya*b6 + za*b9 + com(3)

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       xbg(i) = xbg(i) - x + rgdxin(irgd)
                       ybg(i) = ybg(i) - y + rgdyin(irgd)
                       zbg(i) = zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           End If

        Else If (keyens == 20 .or. keyens == 30 .or. &
                 keyens == 23 .or. keyens == 33) Then

! Langevin and MTK npt/nst

           If (keyens == 20 .or. keyens == 23) Then

              scale = Exp(tstep*eta(1))

              Do j=1,nfree
                 i=lstfre(j)

                 xbg(i) = scale*xbg(i)
                 ybg(i) = scale*ybg(i)
                 zbg(i) = scale*zbg(i)
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 rgdxin(irgd) = scale*rgdxin(irgd)
                 rgdyin(irgd) = scale*rgdyin(irgd)
                 rgdzin(irgd) = scale*rgdzin(irgd)

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       xbg(i) = xbg(i) - x + rgdxin(irgd)
                       ybg(i) = ybg(i) - y + rgdyin(irgd)
                       zbg(i) = zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           Else

! second order taylor expansion of Exp(tstep*eta)

              a1 = tstep*eta(1)
              a2 = tstep*eta(2)
              a3 = tstep*eta(3)
              a5 = tstep*eta(5)
              a6 = tstep*eta(6)
              a9 = tstep*eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do j=1,nfree
                 i=lstfre(j)

                 xa = xbg(i)
                 ya = ybg(i)
                 za = zbg(i)

                 xbg(i) = xa*b1 + ya*b2 + za*b3
                 ybg(i) = xa*b2 + ya*b5 + za*b6
                 zbg(i) = xa*b3 + ya*b6 + za*b9
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 xa = rgdxin(irgd)
                 ya = rgdyin(irgd)
                 za = rgdzin(irgd)

                 rgdxin(irgd) = xa*b1 + ya*b2 + za*b3
                 rgdyin(irgd) = xa*b2 + ya*b5 + za*b6
                 rgdzin(irgd) = xa*b3 + ya*b6 + za*b9

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       xbg(i) = xbg(i) - x + rgdxin(irgd)
                       ybg(i) = ybg(i) - y + rgdyin(irgd)
                       zbg(i) = zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           End If

        End If

! Halo final RB members positions across onto neighbouring domains

        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,xbg,ybg,zbg)

     End If

     Deallocate (rgdxin,rgdyin,rgdzin, Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'xscale deallocation failure, node: ', idnode
        Call error(0)
     End If

  End If

  Call pbcshift(imcon,cell,natms,xin,yin,zin)
  If (llvnl) Call pbcshift(imcon,cell,natms,xbg,ybg,zbg)

End Subroutine xscale
