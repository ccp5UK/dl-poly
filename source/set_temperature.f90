Subroutine set_temperature           &
           (levcfg,temp,keyres,      &
           lmin,nstep,nstrun,nstmin, &
           mxshak,tolnce,keyshl,     &
           atmfre,atmfrz,            &
           megshl,megcon,megpmf,     &
           megrgd,degtra,degrot,     &
           degfre,degshl,sigma,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting the initial system temperature
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsum
  Use setup_module
  Use site_module,        Only : dofsit
  Use config_module,      Only : imcon,natms,nlast,nfree,lsite,  &
                                 lsi,lsa,ltg,lfrzn,lfree,lstfre, &
                                 atmnam,weight,vxx,vyy,vzz
  Use dpd_module,         Only : keydpd
  Use rigid_bodies_module
  Use core_shell_module,  Only : ntshl,listshl,legshl,lshmv_shl,lishp_shl,lashp_shl
  Use kinetic_module,     Only : l_vom,getknr,chvom,getvom

  Implicit None

  Logical,           Intent( In    ) :: lmin
  Integer,           Intent( In    ) :: keyres,nstep,  &
                                        nstrun,nstmin, &
                                        mxshak,keyshl, &
                                        atmfre,atmfrz, &
                                        megshl,        &
                                        megcon,megpmf, &
                                        megrgd
  Real( Kind = wp ), Intent( In    ) :: temp,tolnce

  Integer,           Intent( InOut ) :: levcfg
  Integer(Kind=ip),  Intent( InOut ) :: degtra,degrot

  Integer(Kind=ip),  Intent(   Out ) :: degfre,degshl
  Real( Kind = wp ), Intent(   Out ) :: sigma,engrot

  Logical           :: no_min_0,safe
  Integer           :: fail(1:2),i,j,k,ntp,   &
                       stp,i1,i2,local_index, &
                       irgd,jrgd,lrgd,rgdtyp
  Integer(Kind=ip)  :: com,con,frz,meg,non
  Real( Kind = wp ) :: tmp,vom(1:3),                  &
                       x(1:1),y(1:1),z(1:1),rot(1:9), &
                       vpx,vpy,vpz


! q. index arrays and tp. sum arrays

  Integer,           Allocatable :: qn(:),tpn(:)
  Integer,           Allocatable :: qs(:,:),tps(:)

! initialise rotational and translational DoF if no RB are present
! or re-initialise if all are frozen (does no harm)

  If (megrgd == 0) Then
     degtra=Int(0,ip)
     degrot=Int(0,ip)
  End If

! Degrees Of Freedom (DoF) - all free particles

  meg=Int(3,ip)*Int(atmfre,ip)

! lost to frozen free atoms

  frz=Int(3,ip)*Int(atmfrz,ip)

! 3 lost for fixing COM translation

  If (l_vom .and. keydpd == 0) Then
     com=Int(3,ip)
  Else
     com=Int(0,ip)
  End If

! 3 lost for fixing angular momentum about origin
! (non-periodic systems only)

  non=Int(0,ip)
  If (imcon == 0) non=Int(3,ip)

! lost to shells

  degshl=Int(3,ip)*Int(megshl,ip)

! lost to constrained atoms and PMF constraints

  con=Int(megcon)+Int(megpmf)

! TOTAL DoF

  degfre = meg - com - non - frz - degshl - con + degrot + degtra

! Report DoF

  If (idnode == 0) Then
     Write(nrite,"(/,/,1x,'degrees of freedom break-down list')")
     Write(nrite,"(    1x,'----------------------------------')")
     Write(nrite,"(    1x,'free particles        ',i12)") meg
     Write(nrite,"(    1x,'centre of mass        ',i12)") -com
     Write(nrite,"(    1x,'non-periodicity       ',i12)") -non
     Write(nrite,"(    1x,'frozen free particles ',i12)") -frz
     Write(nrite,"(    1x,'shell-pseudo          ',i12)") -degshl
     Write(nrite,"(    1x,'constrained           ',i12)") -con
     Write(nrite,"(    1x,'RB translational      ',i12)") degtra
     Write(nrite,"(    1x,'RB rotational         ',i12)") degrot
     Write(nrite,"(    1x,'----------------------------------')")
     Write(nrite,"(    1x,'total (real)          ',i12)") degfre
  End If

! Check DoF distribution

  tmp=0.0_wp
  Do i=1,natms
     If (dofsit(lsite(i)) > zero_plus) & ! Omit shells' negative DoFs
        tmp=tmp+dofsit(lsite(i))
  End Do
  If (mxnode > 1) Call gsum(tmp)
  If (Nint(tmp,ip)-non-com /= degfre) Call error(360)

! warn for much restrain on the system [ degfre <= (com+non+frz+con) ]
! and catch an over-restrained system

  If (degfre <= Int(com+non+frz+con,ip)) &
     Call warning(210,Real(degfre,wp),Real((com+non+frz+con),wp),0.0_wp)
  If (degfre < Int(1,ip)) Call error(350)

! desired kinetic energy

  sigma=0.5_wp*Real(degfre,wp)*boltz*temp

  If (keyres == 0) Then

     Allocate (qn(1:mxatms),tpn(0:mxnode-1),    Stat=fail(1))
     Allocate (qs(0:2,1:mxshl),tps(0:mxnode-1), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'set_temperature allocation failure, node: ', idnode
        Call error(0)
     End If

! tpn(idnode) number of particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen particles

     qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
     qs(0:2,1:ntshl) = 0 ! unqualified core-shell unit with a local shell

     j = 0
     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. legshl(1,i) >= 0) Then
           j = j + 1
           qn(i) = 1
        End If
     End Do
     tpn(idnode) = j
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tpn(i) = 0
        End Do
        Call gsum(tpn)
     End If
     ntp = Sum(tpn)

! Save core-shell internal energy vectors due to possible hits
! tps(idnode) number of thermostatted core-shell units on this node (idnode)
! stp - grand total of core-shell units to thermostat

     j = 0
     If (keyshl == 1) Then ! just for the adiabatic shell model
        If (lshmv_shl) Then ! refresh the q array for shared core-shell units
           qn(natms+1:nlast) = 0
           Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,qn)
           Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)
        End If

        If (ntshl > 0) Then
           Do k=1,ntshl
              i1=local_index(listshl(1,k),nlast,lsi,lsa)
              i2=local_index(listshl(2,k),nlast,lsi,lsa)

              If (qn(i1) == 1 .and. i2 > 0 .and. i2 <= natms) Then
                 j = j + 1

                 qs(0,k)=1
                 qs(1,k)=i1
                 qs(2,k)=i2
              End If
           End Do
        End If
     End If
     tps(idnode) = j
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tps(i) = 0
        End Do
        Call gsum(tps)
     End If
     stp = Sum(tps)

     If (megrgd > 0) Then

        k = 0
        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! local index of particle/site
                 If (i <= natms) tpn(idnode) = tpn(idnode) - 1 ! Less free particles are hit
              End Do

              i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
              i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum
              If (rgdfrz(0,rgdtyp) == 0) Then
                 If (i1 <= natms) k = k + 1
              End If
              If (i2 <= natms) k = k + 1
           End If
        End Do
! tpn(idnode) number of thermostatted free particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen free particles to thermostat
        If (mxnode > 1) Then
           Do i=0,mxnode-1
              If (i /= idnode) tpn(i) = 0
           End Do
           Call gsum(tpn)
        End If
        ntp = Sum(tpn)

! generate starting velocities

        Do i=1,natms

! frozen and massless atoms are either motionless
! (with no actual DoF - relaxed shells, frozen sites)
! or despite their DoF, their motion is defined by
! the particles with masses in a RB.  The rest are
! to have gaussian distribution of their Ekin
! (adiabatic shells, constraints, PMFs and RBs are
! to be sorted out later by quenching)

           If (qn(i) == 1 .and. lfree(i) == 0) Then
              Call box_mueller_saru3(ltg(i),0,vxx(i),vyy(i),vzz(i))

! Get scaler to target variance/Sqrt(weight)

              tmp = 1.0_wp/Sqrt(weight(i))
              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
              i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

              If (rgdfrz(0,rgdtyp) == 0 .and. i1 <= natms) Then
                 Call box_mueller_saru3(ltg(i1),0,vxx(i1),vyy(i1),vzz(i1))

! Get scaler to target variance/Sqrt(weight)

                 tmp = 1.0_wp/Sqrt(rgdwgt(0,rgdtyp))
                 vxx(i1) = vxx(i1)*tmp
                 vyy(i1) = vyy(i1)*tmp
                 vzz(i1) = vzz(i1)*tmp
              End If

              If (i2 <= natms) Then
                 Call box_mueller_saru3(ltg(i2),0,vxx(i2),vyy(i2),vzz(i2))

! Get scaler to target variance/Sqrt(weight) -
! 3 different reciprocal moments of inertia

                 vxx(i2) = vxx(i2)*Sqrt(rgdrix(2,rgdtyp))
                 vyy(i2) = vyy(i2)*Sqrt(rgdriy(2,rgdtyp))
                 vzz(i2) = vzz(i2)*Sqrt(rgdriz(2,rgdtyp))
              End If
           End If
        End Do

! Update shared RBs' velocities

        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz)

! calculate new RBs' COM and angular velocities

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
              i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

              If (rgdfrz(0,rgdtyp) == 0) Then
                 rgdvxx(irgd) = vxx(i1)
                 rgdvyy(irgd) = vyy(i1)
                 rgdvzz(irgd) = vzz(i1)
              End If

              rgdoxx(irgd) = vxx(i2)
              rgdoyy(irgd) = vyy(i2)
              rgdozz(irgd) = vzz(i2)
              If (i2 <= natms) Then
                 If (lfrzn(i2) > 0 .or. weight(i) < 1.0e-6_wp) Then
                    vxx(i2) = 0.0_wp
                    vyy(i2) = 0.0_wp
                    vzz(i2) = 0.0_wp
                 End If
              End If

! get new rotation matrix

              Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! update RB members new velocities

              lrgd=listrgd(-1,irgd)
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

     Else ! no RBs present in the system

! generate starting velocities

        Do i=1,natms

! frozen and massless atoms are either motionless
! (with no actual DoF - relaxed shells, frozen sites)
! The rest are to have gaussian distribution of their Ekin
! (adiabatic shells, constraints, PMFs are
! to be sorted out later by quenching)

           If (qn(i) == 1) Then
              Call box_mueller_saru3(ltg(i),0,vxx(i),vyy(i),vzz(i))

! Get scaler to target variance/Sqrt(weight)

              tmp = 1.0_wp/Sqrt(weight(i))
              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

     End If

     If (stp > 0) Then
        If (lshmv_shl) Then ! refresh the q array for shared core-shell units
           qn(natms+1:nlast) = 0
           Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,qn)
           Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)
        End If

        If (tps(idnode) > 0) Then
           Do k=1,ntshl
              If (qs(0,k) == 1) Then
                 i1=qs(1,k)
                 i2=qs(2,k)

                 vxx(i2)=vxx(i1)
                 vyy(i2)=vyy(i1)
                 vzz(i2)=vzz(i1)
              End If
           End Do
        End If
     End If

     Deallocate (qn,tpn, Stat=fail(1))
     Deallocate (qs,tps, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'set_temperature deallocation failure, node: ', idnode
        Call error(0)
     End If

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

  Else

! quench RBs

     If (megrgd > 0) Call rigid_bodies_quench()

  End If

! levcfg must be equalised to ONE (velocities now exist)

  levcfg=1

! scale velocities for keyres=0 and keyres=2

  If (keyres == 0 .or. keyres == 2) Then

! Detect pure molecular statics only == CGM minimisation at zero timestep
! and no rescaling

     no_min_0 = .not.(lmin .and. nstmin == 0 .and. nstep == 0 .and. nstrun == 0 .and. keyres == 2)

! quench constraints & PMFs

     If (no_min_0) Then
        If (megcon > 0) Call constraints_quench(mxshak,tolnce)
        If (megpmf > 0) Call pmf_quench(mxshak,tolnce)
     End If

! quench core-shell units in adiabatic model

     If (megshl > 0 .and. keyshl == 1 .and. no_min_0) Then
        Do
           Call scale_temperature(sigma,degtra,degrot,degfre)
           Call core_shell_quench(safe,temp)
           If (megcon > 0) Call constraints_quench(mxshak,tolnce)
           If (megpmf > 0) Call pmf_quench(mxshak,tolnce)
           If (megrgd > 0) Call rigid_bodies_quench()
           If (safe) Exit
        End Do
     Else
        Call scale_temperature(sigma,degtra,degrot,degfre)
     End If

  End If

  If (l_vom) Then

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
  Else                 ! make getvom always return 0, which is not good for
     Call chvom(l_vom) ! standard MD as the flying ice-cub effect may happen
  End If               ! and/or T(MD) is loses its microscopic meaning!

! Initialise engrot and if RBs exist calculate it

  engrot=0.0_wp
  If (megrgd > 0) engrot=getknr(rgdoxx,rgdoyy,rgdozz)

End Subroutine set_temperature
