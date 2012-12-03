Subroutine pseudo_lfv                                     &
           (isw,keyshl,keyens,keypse,wthpse,tmppse,tstep, &
           strkin,strknf,strknt,engke,engrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to scale the velocities of particles from
! boundary thermostat layer to the target temperature
!
! leapfrog verlet version
!
! Note: (1) This algorithm breaks true ensembles!!!
! Additionally, for Langevin temperature control (keypse=1):
!       (2) Pseudo-randomness of forces is lost - depends on DD!!!
!       (3) Random forces do not contribute to the stress and virial
!           of the system (but are picked up in the pressure).
!       (4) Random forces do not apply to frozen and massless particles
!           as well as shells.
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsum
  Use setup_module,       Only : boltz,nrite,mxatms,mxshl, &
                                 mxtrgd,mxrgd,mxlrgd,zero_plus
  Use site_module,        Only : dofsit,ntpshl,unqshl
  Use config_module
  Use rigid_bodies_module
  Use core_shell_module,  Only : ntshl,listshl,lshmv_shl,lishp_shl,lashp_shl
  Use kinetic_module,     Only : getvom,getknr,kinstress,kinstresf,kinstrest

  Implicit None

  Integer,           Intent( In    ) :: isw,keyshl,keyens,keypse
  Real( Kind = wp ), Intent( In    ) :: tstep,wthpse,tmppse
  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                        strknf(1:9),strknt(1:9),engrot

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ntp,stp,rtp,megrgd,imcon
  Integer                 :: fail(1:3),matms,local_index, &
                             i,j,k,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp )       :: celprp(1:10),ssx,ssy,ssz,      &
                             vom(1:3),scale,tmp,mxdr,       &
                             x(1:1),y(1:1),z(1:1),rot(1:9), &
                             fmx,fmy,fmz,tqx,tqy,tqz,       &
                             trx,try,trz,vpx,vpy,vpz,       &
                             tkin,vdotf,trot,odott,buffer(1:4)
  Real( Kind = wp ), Save :: rcell(1:9),sx,sy,sz,chit = 0.0_wp

! q. index arrays and tp. sum arrays

  Integer, Allocatable, Save     :: qn(:),tpn(:)
  Integer, Allocatable, Save     :: qs(:,:),tps(:)
  Integer, Allocatable, Save     :: qr(:),tpr(:)

! Gaussian random forces arrays

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)

  Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
  Real( Kind = wp ), Allocatable :: rgdfxx(:),rgdfyy(:),rgdfzz(:)
  Real( Kind = wp ), Allocatable :: rgdtxx(:),rgdtyy(:),rgdtzz(:)

  fail = 0

! set matms

  matms=nlast
  If (mxnode == 1) matms=natms

  If (isw == 0) Then

! Random forces cycle - thermostat-system decoupling
! Recalculate the following for a newjob

     If (newjob) Then
        newjob = .false.

! recover megrgd

        megrgd=rgdmeg

! recover imcon

        imcon=rgdimc

        Allocate (qn(1:mxatms),tpn(0:mxnode-1),    Stat=fail(1))
        Allocate (qs(0:2,1:mxshl),tps(0:mxnode-1), Stat=fail(2))
        Allocate (qr(1:mxrgd),tpr(0:mxnode-1),     Stat=fail(3))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'pseudo (q. and tp.) allocation failure, node: ', idnode
           Call error(0)
        End If

        Call invert(cell,rcell,celprp(10))
        Call dcell(cell,celprp)

! The origin of the Coordinate System is in the middle of the MD box (remember)
! Get the real coordinates of the close end edge of the thermostat as if the
! origin of the Coordinate System was in the left-most corner of the MD box

        ssx=wthpse*celprp(1)/celprp(7)
        ssy=wthpse*celprp(2)/celprp(8)
        ssz=wthpse*celprp(3)/celprp(9)

! 1. Get the boundary thermostat thicknesses in fractional coordinates
! 2. sx,sy,sz are intervalled as [-0.5,+0.5) {as by construction are (0,+0.25]}
! 3. Outline the edge beyond which a particle belongs to the thermostat
!    0.5*thicknesses[MD box] - thickness[boundary thermostat]

        sx=rcell(1)*ssx+rcell(4)*ssy+rcell(7)*ssz ; sx=sx-Anint(sx) ; sx=0.5_wp-sx
        sy=rcell(2)*ssx+rcell(5)*ssy+rcell(8)*ssz ; sy=sy-Anint(sy) ; sy=0.5_wp-sy
        sz=rcell(3)*ssx+rcell(6)*ssy+rcell(9)*ssz ; sz=sz-Anint(sz) ; sz=0.5_wp-sz

! qualify non-shell, non-frozen free particles (n) for a random kick
! qualify RBs (r) by them and derive related shells (s)

! tpn(idnode) number of thermostatted particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen particles to thermostat

        qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
        qs(0:2,1:ntshl) = 0 ! unqualified core-shell unit with a local shell
        qr(1:ntrgd)     = 0 ! unqualified RB

        j = 0
        Do i=1,natms

! For all particles on this domain get how far they are
! from the edges of the MD box box

           ssx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i) ; ssx=Abs(ssx-Anint(ssx))
           ssy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i) ; ssy=Abs(ssy-Anint(ssy))
           ssz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i) ; ssz=Abs(ssz-Anint(ssz))

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. (.not.Any(unqshl(1:ntpshl) == atmnam(i))) .and. &
               (ssx >= sx .or. ssy >= sy .or. ssz >= sz)) Then
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
     Else ! rebuild tpn
        j = 0
        Do i=1,natms
           If (qn(i) == 1) j = j + 1
        End Do
        tpn(idnode) = j
        If (mxnode > 1) Then
           Do i=0,mxnode-1
              If (i /= idnode) tpn(i) = 0
           End Do
           Call gsum(tpn)
        End If
        ntp = Sum(tpn)
     End If

     If (ntp == 0) Return
     If (chit < 1.0e-6_wp) Return ! Avoid thermostat overheating

! Allocate random force array of length j

     j=tpn(idnode)
     Allocate (xxt(1:j),yyt(1:j),zzt(1:j), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'pseudo (forces) allocation failure, node: ', idnode
        Call error(0)
     End If

! Here we become node-dependent (using of uni and gauss) - i.e.
! pseudo-randomness depends on the DD mapping which depends on
! the number of nodes and system size
!
! Get gaussian distribution (unit variance)

     Call gauss(j,xxt,yyt,zzt)

! Get scaler to target variance*Sqrt(weight)

     scale = Sqrt(2.0_wp * chit * boltz * tmppse / tstep)

     vom = 0.0_wp
     j = 0
     Do i=1,natms
        If (qn(i) == 1) Then
           j = j + 1

           tmp = scale*Sqrt(weight(i))

           xxt(j) = xxt(j)*tmp
           yyt(j) = yyt(j)*tmp
           zzt(j) = zzt(j)*tmp

           vom(1) = vom(1) + xxt(j)
           vom(2) = vom(2) + yyt(j)
           vom(3) = vom(3) + zzt(j)
        End If
     End Do
     If (mxnode > 1) Call gsum(vom)
     vom = vom / Real(ntp,wp)

! Add random force and remove thermostat COM force

     j = 0
     Do i=1,natms
        If (qn(i) == 1) Then
           j = j + 1

           fxx(i) = fxx(i) + xxt(j) - vom(1)
           fyy(i) = fyy(i) + yyt(j) - vom(2)
           fzz(i) = fzz(i) + zzt(j) - vom(3)
        End If
     End Do

! Deallocate gaussian force array

     Deallocate (xxt,yyt,zzt, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'pseudo (forces) deallocation failure, node: ', idnode
        Call error(0)
     End If

  Else

! Recalculate the following for NPT and NST integration

     If (keyens >= 20) Then
        Call invert(cell,rcell,celprp(10))
        Call dcell(cell,celprp)

! The origin of the Coordinate System is in the middle of the MD box (remember)
! Get the real coordinates of the close end edge of the thermostat as if the
! origin of the Coordinate System was in the left-most corner of the MD box

        ssx=wthpse*celprp(1)/celprp(7)
        ssy=wthpse*celprp(2)/celprp(8)
        ssz=wthpse*celprp(3)/celprp(9)

! 1. Get the boundary thermostat thicknesses in fractional coordinates
! 2. sx,sy,sz are intervalled as [-0.5,+0.5) {as by construction are (0,+0.25]}
! 3. Outline the edge beyond which a particle belongs to the thermostat
!    0.5*thiknesses[MD box] - thickness[boundary thermostat]

        sx=rcell(1)*ssx+rcell(4)*ssy+rcell(7)*ssz ; sx=sx-Anint(sx) ; sx=0.5_wp-sx
        sy=rcell(2)*ssx+rcell(5)*ssy+rcell(8)*ssz ; sy=sy-Anint(sy) ; sy=0.5_wp-sy
        sz=rcell(3)*ssx+rcell(6)*ssy+rcell(9)*ssz ; sz=sz-Anint(sz) ; sz=0.5_wp-sz
     End If

! qualify non-shell, non-frozen free particles (n) for a random kick
! qualify RBs (r) by them and derive related shells (s)

! tpn(idnode) number of thermostatted particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen particles to thermostat

     qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
     qs(0:2,1:ntshl) = 0 ! unqualified core-shell unit with a local shell
     qr(1:ntrgd)     = 0 ! unqualified RB

     j = 0
     Do i=1,natms

! For all particles on this domain get how far they are
! from the edges of the MD box box

        ssx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i) ; ssx=Abs(ssx-Anint(ssx))
        ssy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i) ; ssy=Abs(ssy-Anint(ssy))
        ssz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i) ; ssz=Abs(ssz-Anint(ssz))

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. (.not.Any(unqshl(1:ntpshl) == atmnam(i))) .and. &
            (ssx >= sx .or. ssy >= sy .or. ssz >= sz)) Then
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

     If (ntp == 0) Return

! tps(idnode) number of thermostatted core-shell units on this node (idnode)
! stp - grand total of core-shell units to thermostat

     j = 0
     If (keyshl == 1) Then
        If (lshmv_shl) Then ! refresh the q array for shared core-shell units
           qn(natms+1:nlast) = 0
           Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,qn)
        End If

        If (ntshl > 0) Then
           Do k=1,ntshl
              i1=local_index(listshl(1,k),matms,lsi,lsa)
              i2=local_index(listshl(2,k),matms,lsi,lsa)

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

     j = 0 ! no qualified good RB (one qualified RB is enough to trigger all)
     Do i=1,matms
        If (qn(i) == 1) Then
           If (lfree(i) == 1) j = j + 1
        End If
     End Do
     Call gsum(j)

! tpr(idnode) number of thermostatted RB units on this node (idnode)
! rtp - grand total of RB units to thermostat
! (can be larger than megrgd due to sharing)

     k = 0
     If (j > 0) Then
        If (lshmv_rgd) Then
           qn(natms+1:nlast) = 0 ! refresh the q array for shared RB units
           Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,qn)
        End If

        j = 0
        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! local index of particle/site
                 If (qn(i) == 1) Then
                    If (qr(irgd) == 0) Then ! An overall hit is registered
                       qr(irgd) = 1
                       j = j + 1
                    End If

                    If (i <= natms) tpn(idnode) = tpn(idnode) - 1 ! Less free particles are hit
                 End If
              End Do

              If (qr(irgd) == 1) Then ! accounting for a random kick on the RB
                 i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
                 i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum
                 If (rgdfrz(0,rgdtyp) == 0) Then
                    If (i1 <= natms) k = k + 1
                 End If
                 If (i2 <= natms) k = k + 1
              End If
           End If
        End Do
     End If
! tpn(idnode) number of thermostatted free particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen free particles to thermostat
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tpn(i) = 0
        End Do
        Call gsum(tpn)
     End If
     ntp = Sum(tpn)
     tpr(idnode) = j
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tpr(i) = 0
        End Do
        Call gsum(tpr)
     End If
     rtp = Sum(tpr)

! Velocity scaling cycle - thermostating.  k = local, ntp = global
! number of particles within thermostat layers

     If (keypse == 0 .or. keypse == 1)  Then ! Apply LANGEVIN temperature scaling

! Allocate random velocities array of length k

        Allocate (xxt(1:tpn(idnode)+k),yyt(1:tpn(idnode)+k),zzt(1:tpn(idnode)+k), Stat=fail(1))
        If (fail(1) > 0) Then
           Write(nrite,'(/,1x,a,i0)') 'pseudo (velocities) allocation failure, node: ', idnode
           Call error(0)
        End If

! Here we become node-dependent (using of uni and gauss) - i.e.
! pseudo-randomness depends on the DD mapping which depends on
! the number of nodes and system size
!
! Get gaussian distribution (unit variance)

        Call gauss(tpn(idnode)+k,xxt,yyt,zzt)

        tkin = 0.0_wp
        mxdr = 0.0_wp
        j = 0
        Do i=1,natms
           If (qn(i) == 1 .and. lfree(i) == 0) Then
              If (dofsit(lsite(i)) > zero_plus) mxdr = mxdr + dofsit(lsite(i))

              j = j + 1

! Get scaler to target variance/Sqrt(weight)

              tmp = 1.0_wp/Sqrt(weight(i))
              xxt(j) = xxt(j)*tmp
              yyt(j) = yyt(j)*tmp
              zzt(j) = zzt(j)*tmp

              tkin = tkin + weight(i)*(xxt(j)**2+yyt(j)**2+zzt(j)**2)
           End If
        End Do

        If (rtp > 0) Then
           Do irgd=1,ntrgd
              If (qr(irgd) == 1) Then
                 rgdtyp=listrgd(0,irgd)

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd) ! particle index
                    If (i <= natms) Then
                       If (dofsit(lsite(i)) > zero_plus) mxdr = mxdr + dofsit(lsite(i))
                    End If
                 End Do

                 i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
                 i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

                 If (rgdfrz(0,rgdtyp) == 0 .and. i1 <= natms) Then
                    j = j + 1

                    tmp = 1.0_wp/Sqrt(rgdwgt(0,rgdtyp))
                    vxx(i1) = xxt(j)*tmp
                    vyy(i1) = yyt(j)*tmp
                    vzz(i1) = zzt(j)*tmp

                    tkin = tkin + rgdwgt(0,rgdtyp)*(vxx(i1)**2+vyy(i1)**2+vzz(i1)**2)
                 End If

                 If (i2 <= natms) Then
                    j = j + 1

                    vxx(i2) = xxt(j)*Sqrt(rgdrix(2,rgdtyp))
                    vyy(i2) = yyt(j)*Sqrt(rgdriy(2,rgdtyp))
                    vzz(i2) = zzt(j)*Sqrt(rgdriz(2,rgdtyp))

                    tkin = tkin + (rgdrix(1,rgdtyp)*vxx(i2)**2+rgdriy(1,rgdtyp)*vyy(i2)**2+rgdriz(1,rgdtyp)*vzz(i2)**2)
                 End If
              End If
           End Do
        End If

        If (mxnode > 1) Call gsum(tkin)
        If (tkin <= zero_plus) tkin = 1.0_wp
        If (mxnode > 1) Call gsum(mxdr)

! Scale to target temperature and apply thermostat

        scale = Sqrt(mxdr * boltz * tmppse / tkin)

! Scale velocity within the thermostat layer to the gaussian velocities
! scaled with the variance for the target temperature

        j = 0
        Do i=1,natms
           If (qn(i) == 1 .and. lfree(i) == 0) Then
              j = j + 1

              tmp = scale * Sqrt( (xxt(j)**2+yyt(j)**2+zzt(j)**2) / &
                                  (vxx(i)**2+vyy(i)**2+vzz(i)**2) )

              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

! Deallocate gaussian random velocities array

        Deallocate (xxt,yyt,zzt, Stat=fail(1))
        If (fail(1) > 0) Then
           Write(nrite,'(/,1x,a,i0)') 'pseudo (velocities) deallocation failure, node: ', idnode
           Call error(0)
        End If

        If (rtp > 0) Then

! Update shared RBs' velocities

           If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz)

! calculate new RBs' COM and angular velocities

           Do irgd=1,ntrgd
              If (qr(irgd) == 1) Then
                 rgdtyp=listrgd(0,irgd)

                 i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
                 i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

                 If (rgdfrz(0,rgdtyp) == 0) Then
                    tmp = scale * Sqrt( (vxx(i1)**2+vyy(i1)**2+vzz(i1)**2) / &
                                        (rgdvxx(irgd)**2+rgdvyy(irgd)**2+rgdvzz(irgd)**2) )
                    rgdvxx(irgd) = rgdvxx(irgd)*tmp
                    rgdvyy(irgd) = rgdvyy(irgd)*tmp
                    rgdvzz(irgd) = rgdvzz(irgd)*tmp
                 End If

                 tmp = scale * Sqrt( (rgdrix(1,rgdtyp)*vxx(i2)**2+      &
                                      rgdriy(1,rgdtyp)*vyy(i2)**2+      &
                                      rgdriz(1,rgdtyp)*vzz(i2)**2) /    &
                                     (rgdrix(1,rgdtyp)*rgdoxx(irgd)**2+ &
                                      rgdriy(1,rgdtyp)*rgdoyy(irgd)**2+ &
                                      rgdriz(1,rgdtyp)*rgdozz(irgd)**2) )
                 rgdoxx(irgd) = rgdoxx(irgd)*tmp
                 rgdoyy(irgd) = rgdoyy(irgd)*tmp
                 rgdozz(irgd) = rgdozz(irgd)*tmp
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

        End If

! Thermalise the shells on hit cores

        If (stp > 0) Then
           If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)

           If (tps(idnode) > 0) Then
              j = 0
              Do k=1,ntshl
                 If (qs(0,k) == 1) Then
                    j = j + 1

                    i1=qs(1,k)
                    i2=qs(2,k)

                    vxx(i2)=vxx(i1)
                    vyy(i2)=vyy(i1)
                    vzz(i2)=vzz(i1)
                 End If
              End Do
           End If
        End If

        If (megrgd > 0) Then

! remove system centre of mass velocity (random momentum walk)

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

! Get friction for the Langevin thermostat as if gaussian thermal constraints were applied

           tkin  = 0.0_wp
           vdotf = 0.0_wp
           trot  = 0.0_wp
           odott = 0.0_wp

! Free particles

           Do j=1,nfree
              i=lstfre(j)

              If (qn(i) == 1) Then
                 vxx(i) = vxx(i)
                 vyy(i) = vyy(i)
                 vzz(i) = vzz(i)

                 tkin  = tkin  + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
                 vdotf = vdotf + vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
              End If
           End Do

! Shells

           If (stp > 0) Then
              If (tps(idnode) > 0) Then
                 Do k=1,ntshl
                    If (qs(0,k) == 1) Then
                       i2=qs(2,k)

                       tkin  = tkin  + weight(i2)*(vxx(i2)**2+vyy(i2)**2+vzz(i2)**2)
                       vdotf = vdotf + vxx(i2)*fxx(i2)+vyy(i2)*fyy(i2)+vzz(i2)*fzz(i2)
                    End If
                 End Do
              End If
           End If

! RBs

           If (rtp > 0) Then
              Allocate (ggx(1:mxlrgd*Max(mxrgd,mxtrgd)),ggy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
                        ggz(1:mxlrgd*Max(mxrgd,mxtrgd)),                 Stat=fail(1))
              Allocate (rgdfxx(1:mxrgd),rgdfyy(1:mxrgd),rgdfzz(1:mxrgd), Stat=fail(2))
              Allocate (rgdtxx(1:mxrgd),rgdtyy(1:mxrgd),rgdtzz(1:mxrgd), Stat=fail(3))
              If (Any(fail > 0)) Then
                 Write(nrite,'(/,1x,a,i0)') 'pseudo (RB) allocation failure, node: ', idnode
                 Call error(0)
              End If

! Get the RB particles vectors wrt the RB's COM

              krgd=0
              Do irgd=1,ntrgd
                 If (qr(irgd) == 1) Then
                    rgdtyp=listrgd(0,irgd)
                    lrgd=listrgd(-1,irgd)

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

! Get RB force and torque

              krgd=0
              Do irgd=1,ntrgd
                 If (qr(irgd) == 1) Then
                    rgdtyp=listrgd(0,irgd)
                    lrgd=listrgd(-1,irgd)

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

! current rotation matrix

                    Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! calculate torque in principal frame

                    trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
                    try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
                    trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

! store COM force and torque

                    rgdfxx(irgd)=fmx
                    rgdfyy(irgd)=fmy
                    rgdfzz(irgd)=fmz

                    rgdtxx(irgd)=trx
                    rgdtyy(irgd)=try
                    rgdtzz(irgd)=trz
                 Else
                    rgdfxx(irgd)=0.0_wp
                    rgdfyy(irgd)=0.0_wp
                    rgdfzz(irgd)=0.0_wp

                    rgdtxx(irgd)=0.0_wp
                    rgdtyy(irgd)=0.0_wp
                    rgdtzz(irgd)=0.0_wp
                 End If
              End Do

! Energetics

              Do irgd=1,ntrgd
                 If (qr(irgd) == 1) Then
                    rgdtyp=listrgd(0,irgd)
                    lrgd=listrgd(-1,irgd)

                    tmp=Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

                    If (rgdfrz(0,rgdtyp) == 0) Then
                       tkin  = tkin  + tmp * &
                       rgdwgt(0,rgdtyp)*(rgdvxx(irgd)**2+rgdvyy(irgd)**2+rgdvzz(irgd)**2)
                       vdotf = vdotf + tmp * &
                       (rgdvxx(irgd)*rgdfxx(irgd)+rgdvyy(irgd)*rgdfyy(irgd)+rgdvzz(irgd)*rgdfzz(irgd))
                    End If

                    trot  = trot  + tmp * &
                    (rgdrix(1,rgdtyp)*rgdoxx(irgd)**2+rgdriy(1,rgdtyp)*rgdoyy(irgd)**2+rgdriz(1,rgdtyp)*rgdozz(irgd)**2)
                    odott = odott + tmp * &
                    (rgdoxx(irgd)*rgdtxx(irgd)+rgdoyy(irgd)*rgdtyy(irgd)+rgdozz(irgd)*rgdtzz(irgd))
                 End If
              End Do

              Deallocate (ggx,ggy,ggz,          Stat=fail(1))
              Deallocate (rgdfxx,rgdfyy,rgdfzz, Stat=fail(2))
              Deallocate (rgdtxx,rgdtyy,rgdtzz, Stat=fail(3))
              If (Any(fail > 0)) Then
                 Write(nrite,'(/,1x,a,i0)') 'pseudo (RB) deallocation failure, node: ', idnode
                 Call error(0)
              End If
           End If

           If (mxnode > 1) Then
              buffer(1) = tkin
              buffer(2) = vdotf
              buffer(3) = trot
              buffer(4) = odott
              Call gsum(buffer(1:4))
              tkin  = buffer(1)
              vdotf = buffer(2)
              trot  = buffer(3)
              odott = buffer(4)
           End If
           tmp=tkin+trot
           If (tmp <= zero_plus) tmp = 1.0_wp
           chit = (vdotf+odott)/tmp

        Else

! remove system centre of mass velocity (random momentum walk)

           Call getvom(vom,vxx,vyy,vzz)

           Do i=1,natms
              If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                 vxx(i) = vxx(i) - vom(1)
                 vyy(i) = vyy(i) - vom(2)
                 vzz(i) = vzz(i) - vom(3)
              End If
           End Do

! Get friction for the Langevin thermostat as if gaussian thermal constraints were applied

           tkin  = 0.0_wp
           vdotf = 0.0_wp

           Do i=1,natms
              If (qn(i) == 1) Then
                 vxx(i) = vxx(i)
                 vyy(i) = vyy(i)
                 vzz(i) = vzz(i)

                 tkin   = tkin  + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
                 vdotf  = vdotf + vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
              End If
           End Do

! Shells

           If (stp > 0) Then
              If (tps(idnode) > 0) Then
                 Do k=1,ntshl
                    If (qs(0,k) == 1) Then
                       i2=qs(2,k)

                       tkin  = tkin  + weight(i2)*(vxx(i2)**2+vyy(i2)**2+vzz(i2)**2)
                       vdotf = vdotf + vxx(i2)*fxx(i2)+vyy(i2)*fyy(i2)+vzz(i2)*fzz(i2)
                    End If
                 End Do
              End If
           End If

           If (mxnode > 1) Then
              buffer(1) = tkin
              buffer(2) = vdotf
              Call gsum(buffer(1:2))
              tkin  = buffer(1)
              vdotf = buffer(2)
           End If
           If (tkin <= zero_plus) tkin = 1.0_wp
           chit = vdotf/tkin

        End If

     End If

     If (keypse == 0 .or. keypse == 2) Then ! Apply DIRECT temperature scaling

! Targeted energy

        scale = boltz * tmppse

! Scale velocity within the thermostat layer

        Do i=1,natms
           If (qn(i) == 1 .and. lfree(i) == 0) Then

! Get particle kinetic energy and produce a scaler to target temperature

              tmp = Sqrt(scale * dofsit(lsite(i)) / (weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)))

              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

        If (rtp > 0) Then

! calculate new RBs' COM and angular velocities

           Do irgd=1,ntrgd
              If (qr(irgd) == 1) Then
                 rgdtyp=listrgd(0,irgd)

                 i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
                 i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

                 mxdr = 3.0_wp
                 If      (rgdfrz(0,rgdtyp) == 0) Then
                    tmp = Sqrt( scale * mxdr / &
                    (rgdwgt(0,irgd)*(rgdvxx(irgd)**2+rgdvyy(irgd)**2+rgdvzz(irgd)**2)) )
                    rgdvxx(irgd) = rgdvxx(irgd)*tmp
                    rgdvyy(irgd) = rgdvyy(irgd)*tmp
                    rgdvzz(irgd) = rgdvzz(irgd)*tmp
                 Else If (rgdfrz(0,rgdtyp) >  1) Then
                    mxdr = 1.0_wp
                 End If

                 tmp = Sqrt( scale * mxdr / &
                 (rgdrix(1,rgdtyp)*rgdoxx(irgd)**2+ &
                  rgdriy(1,rgdtyp)*rgdoyy(irgd)**2+ &
                  rgdriz(1,rgdtyp)*rgdozz(irgd)**2) )
                 rgdoxx(irgd) = rgdoxx(irgd)*tmp
                 rgdoyy(irgd) = rgdoyy(irgd)*tmp
                 rgdozz(irgd) = rgdozz(irgd)*tmp

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

        End If

! Thermalise the shells on hit cores

        If (stp > 0) Then
           If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)

           If (tps(idnode) > 0) Then
              j = 0
              Do k=1,ntshl
                 If (qs(0,k) == 1) Then
                    j = j + 1

                    i1=qs(1,k)
                    i2=qs(2,k)

                    vxx(i2)=vxx(i1)
                    vyy(i2)=vyy(i1)
                    vzz(i2)=vzz(i1)
                 End If
              End Do
           End If
        End If

! remove system centre of mass velocity (random momentum walk)

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

     End If

! Update total kinetic stress and energy

     If (megrgd > 0) Then
        Call kinstresf(vxx,vyy,vzz,strknf)
        Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

        strkin=strknf+strknt

! update rotational energy

        engrot=getknr(rgdoxx,rgdoyy,rgdozz)
     Else
        Call kinstress(vxx,vyy,vzz,strkin)
     End If
     engke = 0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  End If

End Subroutine pseudo_lfv
