Module stochastic_boundary
  Use kinds, Only : wp,wi
  Use comms, Only : comms_type,gsum
  use kinetics, Only : kinstresf, kinstrest, kinstress,getknr,getvom, getknr
  Use configuration, Only : vxx,vyy,vzz,cell,natms,&
    weight,lfrzn,lstfre,lsite,lfree,ltg,nlast,nfree,&
    lsi,lsa,imcon
  Use errors_warnings, Only : error,warning,info
  Use shared_units, Only : update_shared_units,update_shared_units_int
  Use numerics, Only : local_index,images,dcell,invert,box_mueller_saru3
  Use thermostat, Only : thermostat_type
  Use statistics, Only : stats_type
  Use domains, Only : domains_type
  Use core_shell,  Only : core_shell_type,SHELL_ADIABATIC 
  Use particle, Only : corePart
  Use rigid_bodies, Only : rigid_bodies_type,getrotmat
  Use setup, Only : mxatms,zero_plus,boltz

  Implicit None
  Private
  
  Public :: stochastic_boundary_vv
  
  Contains
  
  Subroutine stochastic_boundary_vv(isw,tstep,nstep,dof_site,cshell,stats,thermo,rigid,domain,parts,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to scale the velocities of particles from
! boundary thermostat layer to the target temperature
!
! velocity verlet version
!
! Note: (1) This algorithm breaks true ensembles!!!
! Additionally, for Langevin temperature control (thermo%key_pseudo=1):
!       (2) Random forces do not contribute to the stress and virial
!           of the thermo%system (but are picked up in the pressure).
!       (3) Random forces do not apply to frozen and massless particles
!           as well as to shells.
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: isw,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep
  Real( Kind = wp ), Dimension(:), Intent( In    ) :: dof_site
  Type( stats_type), Intent( InOut ) :: stats
  Type( core_shell_type), Intent( InOut ) :: cshell
  Type( thermostat_type ), Intent( InOut    ) :: thermo
  Type( rigid_bodies_type ), Intent( InOut ) :: rigid
  Type( domains_type ), Intent( In    ) :: domain
  Type( comms_type), Intent( InOut ) :: comm
  Type (corePart ), Intent( InOut ) :: parts(:)

  Integer                 :: fail(1:3),matms, &
                             i,j,k,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp )       :: celprp(1:10),ssx,ssy,ssz,      &
                             vom(1:3),scale,tmp,mxdr,       &
                             x(1:1),y(1:1),z(1:1),rot(1:9), &
                             fmx,fmy,fmz,tqx,tqy,tqz,       &
                             trx,try,trz,vpx,vpy,vpz,       &
                             tkin,vdotf,trot,odott,buffer(1:4)


! Gaussian random forces arrays

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)

  Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
  Real( Kind = wp ), Allocatable :: rgdfxx(:),rgdfyy(:),rgdfzz(:)
  Real( Kind = wp ), Allocatable :: rgdtxx(:),rgdtyy(:),rgdtzz(:)

  Character ( Len = 256 ) :: message

  fail = 0

! set matms

  matms=nlast
  If (comm%mxnode == 1) matms=natms

  If (isw == 0) Then

! Random forces cycle - thermostat-thermo%system decoupling
! Recalculate the following for NPT and NST integration

    If (thermo%variable_cell .or. thermo%newjob) Then
       If (thermo%newjob) Then
         thermo%newjob = .false.

           Allocate (thermo%qn(1:mxatms),thermo%tpn(0:comm%mxnode-1),    Stat=fail(1))
           Allocate (thermo%qs(0:2,1:cshell%mxshl),thermo%tps(0:comm%mxnode-1), Stat=fail(2))
           Allocate (thermo%qr(1:rigid%max_rigid),thermo%tpr(0:comm%mxnode-1),     Stat=fail(3))
           If (Any(fail > 0)) Then
              Write(message,'(a)') 'pseudo (q. and tp.) allocation failure'
              Call error(0,message)
           End If
        End If

        Call invert(cell,thermo%rcell,celprp(10))
        Call dcell(cell,celprp)

! The origin of the Coordinate System is in the middle of the MD box (remember)
! Get the real coordinates of the close end edge of the thermostat as if the
! origin of the Coordinate System was in the left-most corner of the MD box

        ssx=thermo%width_pseudo*celprp(1)/celprp(7)
        ssy=thermo%width_pseudo*celprp(2)/celprp(8)
        ssz=thermo%width_pseudo*celprp(3)/celprp(9)

! 1. Get the boundary thermostat thicknesses in fractional coordinates
! 2. thermo%sx,thermo%sy,thermo%sz are intervalled as [-0.5,+0.5) {as by construction are (0,+0.25]}
! 3. Outline the edge beyond which a particle belongs to the thermostat
!    0.5*thicknesses[MD box] - thickness[boundary thermostat]

        thermo%sx=thermo%rcell(1)*ssx+thermo%rcell(4)*ssy+thermo%rcell(7)*ssz ; thermo%sx=thermo%sx-Anint(thermo%sx)  
        thermo%sx=0.5_wp-thermo%sx
        thermo%sy=thermo%rcell(2)*ssx+thermo%rcell(5)*ssy+thermo%rcell(8)*ssz ; thermo%sy=thermo%sy-Anint(thermo%sy) 
        thermo%sy=0.5_wp-thermo%sy
        thermo%sz=thermo%rcell(3)*ssx+thermo%rcell(6)*ssy+thermo%rcell(9)*ssz ; thermo%sz=thermo%sz-Anint(thermo%sz) 
        thermo%sz=0.5_wp-thermo%sz
     End If

! qualify non-shell, non-frozen free particles (n) for a random kick
! qualify RBs (r) by them and derive related shells (s)

! thermo%tpn(comm%idnode) number of thermostatted particles on this node (comm%idnode)
! thermo%ntp - grand total of non-shell, non-frozen particles to thermostat

     thermo%qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
     thermo%qs(0:2,1:cshell%ntshl) = 0 ! unqualified core-shell unit with a local shell
     thermo%qr(1:rigid%n_types)     = 0 ! unqualified RB

     j = 0
     Do i=1,natms

! For all particles on this domain get how far they are
! from the origin of the MD box

        ssx=thermo%rcell(1)*parts(i)%xxx+thermo%rcell(4)*parts(i)%yyy+thermo%rcell(7)*parts(i)%zzz ; ssx=Abs(ssx-Anint(ssx))
        ssy=thermo%rcell(2)*parts(i)%xxx+thermo%rcell(5)*parts(i)%yyy+thermo%rcell(8)*parts(i)%zzz ; ssy=Abs(ssy-Anint(ssy))
        ssz=thermo%rcell(3)*parts(i)%xxx+thermo%rcell(6)*parts(i)%yyy+thermo%rcell(9)*parts(i)%zzz ; ssz=Abs(ssz-Anint(ssz))

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0 .and. &
            (ssx >= thermo%sx .or. ssy >= thermo%sy .or. ssz >= thermo%sz)) Then
           j = j + 1
           thermo%qn(i) = 1
        End If
     End Do
     thermo%tpn(comm%idnode) = j
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) thermo%tpn(i) = 0
        End Do
        Call gsum(comm,thermo%tpn)
     End If
     thermo%ntp = Sum(thermo%tpn)

     If (thermo%ntp == 0) Return
     If (thermo%chit_sb < 1.0e-6_wp .or. thermo%key_pseudo > 1) Return ! Avoid thermostat overheating

! Allocate random force array of length j

     j=thermo%tpn(comm%idnode)
     Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'pseudo (forces) allocation failure'
        Call error(0,message)
     End If

! Get gaussian distribution
! Get scaler to target variance*Sqrt(weight)

     scale = Sqrt(2.0_wp * thermo%chit_sb * boltz * thermo%temp_pseudo / tstep)

     vom = 0.0_wp
     Do i=1,natms
        If (thermo%qn(i) == 1) Then

! Get gaussian distribution (unit variance)

           Call box_mueller_saru3(ltg(i),nstep,xxt(i),yyt(i),zzt(i))

! Get scaler to target variance*Sqrt(weight)

           tmp = scale*Sqrt(weight(i))

           xxt(i) = xxt(i)*tmp
           yyt(i) = yyt(i)*tmp
           zzt(i) = zzt(i)*tmp

           vom(1) = vom(1) + xxt(i)
           vom(2) = vom(2) + yyt(i)
           vom(3) = vom(3) + zzt(i)

        End If
     End Do
     Call gsum(comm,vom)
     vom = vom / Real(thermo%ntp,wp)

! Add random force and remove thermostat COM force

     Do i=1,natms
        If (thermo%qn(i) == 1) Then
           parts(i)%fxx = parts(i)%fxx + xxt(i) - vom(1)
           parts(i)%fyy = parts(i)%fyy + yyt(i) - vom(2)
           parts(i)%fzz = parts(i)%fzz + zzt(i) - vom(3)
        End If
     End Do

! Deallocate gaussian force array

     Deallocate (xxt,yyt,zzt, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'pseudo (forces) deallocation failure'
        Call error(0,message)
     End If

  Else

     If (thermo%ntp == 0) Return

     j = 0 ! no qualified good RB (one qualified RB is enough to trigger all)
     Do i=1,matms
        If (thermo%qn(i) == 1) Then
           If (lfree(i) == 1) j = j + 1
        End If
     End Do
     Call gsum(comm,j)

! thermo%tpr(comm%idnode) number of thermostatted RB units on this node (comm%idnode)
! thermo%rtp - grand total of RB units to thermostat
! (can be larger than rigid%total due to sharing)

     k = 0
     If (j > 0) Then
        If (rigid%share) Then
           thermo%qn(natms+1:nlast) = 0 ! refresh the q array for shared RB units
           Call update_shared_units_int(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,thermo%qn,domain,comm)
        End If

        j = 0
        Do irgd=1,rigid%n_types
           rgdtyp=rigid%list(0,irgd)

! For all good RBs

           lrgd=rigid%list(-1,irgd)
           If (rigid%frozen(0,rgdtyp) < lrgd) Then
              Do jrgd=1,lrgd
                 i=rigid%index_local(jrgd,irgd) ! local index of particle/site
                 If (thermo%qn(i) == 1) Then
                    If (thermo%qr(irgd) == 0) Then ! An overall hit is registered
                       thermo%qr(irgd) = 1
                       j = j + 1
                    End If

                    If (i <= natms) thermo%tpn(comm%idnode) = thermo%tpn(comm%idnode) - 1 ! Less free particles are hit
                 End If
              End Do

              If (thermo%qr(irgd) == 1) Then ! accounting for a random kick on the RB
                 i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                 i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum
                 If (rigid%frozen(0,rgdtyp) == 0) Then
                    If (i1 <= natms) k = k + 1
                 End If
                 If (i2 <= natms) k = k + 1
              End If
           End If
        End Do
     End If
! thermo%tpn(comm%idnode) number of thermostatted free particles on this node (comm%idnode)
! thermo%ntp - grand total of non-shell, non-frozen free particles to thermostat
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) thermo%tpn(i) = 0
        End Do
        Call gsum(comm,thermo%tpn)
     End If
     thermo%ntp = Sum(thermo%tpn)
     thermo%tpr(comm%idnode) = j
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) thermo%tpr(i) = 0
        End Do
        Call gsum(comm,thermo%tpr)
     End If
     thermo%rtp = Sum(thermo%tpr)

! thermo%tps(comm%idnode) number of thermostatted core-shell units on this node (comm%idnode)
! thermo%stp - grand total of core-shell units to thermostat

     j = 0
     If (cshell%keyshl == SHELL_ADIABATIC) Then
        If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
           thermo%qn(natms+1:nlast) = 0
           Call update_shared_units_int(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,thermo%qn,domain,comm)
        End If

        If (cshell%ntshl > 0) Then
           Do k=1,cshell%ntshl
              i1=local_index(cshell%listshl(1,k),matms,lsi,lsa)
              i2=local_index(cshell%listshl(2,k),matms,lsi,lsa)

              If (thermo%qn(i1) == 1 .and. i2 > 0 .and. i2 <= natms) Then
                 j = j + 1

                 thermo%qs(0,k)=1
                 thermo%qs(1,k)=i1
                 thermo%qs(2,k)=i2
              End If
           End Do
        End If
     End If
     thermo%tps(comm%idnode) = j
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) thermo%tps(i) = 0
        End Do
        Call gsum(comm,thermo%tps)
     End If
     thermo%stp = Sum(thermo%tps)

! Velocity scaling cycle - thermostatting.  k = local, thermo%ntp = global
! number of particles within thermostat layers

     If (thermo%key_pseudo < 3)  Then ! Apply LANGEVIN temperature scaling

! Allocate random velocities array of length k

        Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail(1))
        If (fail(1) > 0) Then
           Write(message,'(a,i0)') 'pseudo (velocities) allocation failure'
           Call error(0,message)
        End If

! Here we become node-dependent (using of uni and gauss) - i.e.
! pseudo-randomness depends on the DD mapping which depends on
! the number of nodes and thermo%system size
!
! Get gaussian distribution

        tkin = 0.0_wp
        mxdr = 0.0_wp
        Do i=1,natms
           If (thermo%qn(i) == 1 .and. lfree(i) == 0) Then
              If (dof_site(lsite(i)) > zero_plus) mxdr = mxdr + dof_site(lsite(i))

! Get gaussian distribution (unit variance)

              Call box_mueller_saru3(ltg(i),nstep,xxt(i),yyt(i),zzt(i))

! Get scaler to target variance/Sqrt(weight)

              tmp = 1.0_wp/Sqrt(weight(i))

              xxt(i) = xxt(i)*tmp
              yyt(i) = yyt(i)*tmp
              zzt(i) = zzt(i)*tmp

              tkin = tkin + weight(i)*(xxt(i)**2+yyt(i)**2+zzt(i)**2)
           End If
        End Do

        If (thermo%rtp > 0) Then
           Do irgd=1,rigid%n_types
              If (thermo%qr(irgd) == 1) Then
                 rgdtyp=rigid%list(0,irgd)

                 lrgd=rigid%list(-1,irgd)
                 Do jrgd=1,lrgd
                    i=rigid%index_local(jrgd,irgd) ! particle index
                    If (i <= natms) Then
                       If (dof_site(lsite(i)) > zero_plus) mxdr = mxdr + dof_site(lsite(i))
                    End If
                 End Do

                 i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                 i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                 If (rigid%frozen(0,rgdtyp) == 0 .and. i1 <= natms) Then

! Get gaussian distribution (unit variance)

                    Call box_mueller_saru3(ltg(i1),nstep,xxt(i1),yyt(i1),zzt(i1))

! Get scaler to target variance/Sqrt(weight)

                    tmp = 1.0_wp/Sqrt(rigid%weight(0,rgdtyp))

                    xxt(i1) = xxt(i1)*tmp
                    yyt(i1) = yyt(i1)*tmp
                    zzt(i1) = zzt(i1)*tmp

                    tkin = tkin + rigid%weight(0,rgdtyp)*(xxt(i1)**2+yyt(i1)**2+zzt(i1)**2)

                 End If

                 If (i2 <= natms) Then

! Get gaussian distribution (unit variance)

                    Call box_mueller_saru3(ltg(i2),nstep,xxt(i2),yyt(i2),zzt(i2))

! Get scaler to target variance/Sqrt(weight) -
! 3 different reciprocal moments of inertia

                    xxt(i2) = xxt(i2)*Sqrt(rigid%rix(2,rgdtyp))
                    yyt(i2) = yyt(i2)*Sqrt(rigid%riy(2,rgdtyp))
                    zzt(i2) = zzt(i2)*Sqrt(rigid%riz(2,rgdtyp))

                    tkin = tkin + (rigid%rix(1,rgdtyp)*xxt(i2)**2+ &
                      rigid%riy(1,rgdtyp)*yyt(i2)**2+ &
                      rigid%riz(1,rgdtyp)*zzt(i2)**2)

                 End If
              End If
           End Do
        End If

        Call gsum(comm,tkin)
        If (tkin <= zero_plus) tkin = 1.0_wp
        Call gsum(comm,mxdr)

! Scale to target temperature and apply thermostat

        scale = Sqrt(mxdr * boltz * thermo%temp_pseudo / tkin)

! Scale velocity within the thermostat layer to the gaussian velocities
! scaled with the variance for the target temperature

        Do i=1,natms
           If (thermo%qn(i) == 1 .and. lfree(i) == 0) Then
              tmp = scale * Sqrt( (xxt(i)**2+yyt(i)**2+zzt(i)**2) / &
                                  (vxx(i)**2+vyy(i)**2+vzz(i)**2) )

              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

        If (thermo%rtp > 0) Then

! Update shared RBs' velocities

           If (rigid%share) Then
             Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,xxt,yyt,zzt,domain,comm)
           End If

! calculate new RBs' COM and angular velocities

           Do irgd=1,rigid%n_types
              If (thermo%qr(irgd) == 1) Then
                 rgdtyp=rigid%list(0,irgd)

                 i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                 i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                 If (rigid%frozen(0,rgdtyp) == 0) Then
                    tmp = scale * Sqrt( (xxt(i1)**2+yyt(i1)**2+zzt(i1)**2) / &
                                        (rigid%vxx(irgd)**2+rigid%vyy(irgd)**2+rigid%vzz(irgd)**2) )
                    rigid%vxx(irgd) = rigid%vxx(irgd)*tmp
                    rigid%vyy(irgd) = rigid%vyy(irgd)*tmp
                    rigid%vzz(irgd) = rigid%vzz(irgd)*tmp
                 End If

                 tmp = scale * Sqrt( (rigid%rix(1,rgdtyp)*xxt(i2)**2+      &
                                      rigid%riy(1,rgdtyp)*yyt(i2)**2+      &
                                      rigid%riz(1,rgdtyp)*zzt(i2)**2) /    &
                                     (rigid%rix(1,rgdtyp)*rigid%oxx(irgd)**2+ &
                                      rigid%riy(1,rgdtyp)*rigid%oyy(irgd)**2+ &
                                      rigid%riz(1,rgdtyp)*rigid%ozz(irgd)**2) )
                 rigid%oxx(irgd) = rigid%oxx(irgd)*tmp
                 rigid%oyy(irgd) = rigid%oyy(irgd)*tmp
                 rigid%ozz(irgd) = rigid%ozz(irgd)*tmp

! get new rotation matrix

                 Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

! update RB members new velocities

                 lrgd=rigid%list(-1,irgd)
                 Do jrgd=1,lrgd
                    If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                       i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                       If (i <= natms) Then
                          x(1)=rigid%x(jrgd,rgdtyp)
                          y(1)=rigid%y(jrgd,rgdtyp)
                          z(1)=rigid%z(jrgd,rgdtyp)

! new atomic velocities in body frame

                          vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                          vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                          vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

! new atomic velocities in lab frame

                          vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                          vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                          vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                       End If
                    End If
                 End Do
              End If
           End Do

        End If

! Deallocate gaussian random velocities array

        Deallocate (xxt,yyt,zzt, Stat=fail(1))
        If (fail(1) > 0) Then
           Write(message,'(a)') 'pseudo (velocities) deallocation failure'
           Call error(0,message)
        End If

! Thermalise the shells on hit cores

        If (thermo%stp > 0) Then
           If (cshell%lshmv_shl) Then
             Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,vxx,vyy,vzz,domain,comm)
           End If

           If (thermo%tps(comm%idnode) > 0) Then
              j = 0
              Do k=1,cshell%ntshl
                 If (thermo%qs(0,k) == 1) Then
                    j = j + 1

                    i1=thermo%qs(1,k)
                    i2=thermo%qs(2,k)

                    vxx(i2)=vxx(i1)
                    vyy(i2)=vyy(i1)
                    vzz(i2)=vzz(i1)
                 End If
              End Do
           End If
        End If

        If (rigid%total > 0) Then

! remove thermo%system centre of mass velocity (random momentum walk)

           Call getvom(vom,vxx,vyy,vzz,rigid,comm)

           Do j=1,nfree
              i=lstfre(j)

              If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                 vxx(i) = vxx(i) - vom(1)
                 vyy(i) = vyy(i) - vom(2)
                 vzz(i) = vzz(i) - vom(3)
              End If
           End Do

           Do irgd=1,rigid%n_types
              rgdtyp=rigid%list(0,irgd)

              If (rigid%frozen(0,rgdtyp) == 0) Then
                 rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
                 rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
                 rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

                 lrgd=rigid%list(-1,irgd)
                 Do jrgd=1,lrgd
                    i=rigid%index_local(jrgd,irgd) ! local index of particle/site

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

              If (thermo%qn(i) == 1) Then
                 tkin  = tkin  + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
                 vdotf = vdotf + vxx(i)*parts(i)%fxx+vyy(i)*parts(i)%fyy+vzz(i)*parts(i)%fzz
              End If
           End Do

! RBs

           If (thermo%rtp > 0) Then
              Allocate (ggx(1:rigid%max_list*rigid%max_rigid), &
                ggy(1:rigid%max_list*rigid%max_rigid), &
                ggz(1:rigid%max_list*rigid%max_rigid), Stat=fail(1))
              Allocate (rgdfxx(1:rigid%max_rigid), &
                rgdfyy(1:rigid%max_rigid), &
                rgdfzz(1:rigid%max_rigid), Stat=fail(2))
              Allocate (rgdtxx(1:rigid%max_rigid), &
                rgdtyy(1:rigid%max_rigid), &
                rgdtzz(1:rigid%max_rigid), Stat=fail(3))
              If (Any(fail > 0)) Then
                 Write(message,'(a)') 'pseudo (RB) allocation failure'
                 Call error(0,message)
              End If

! Get the RB particles vectors wrt the RB's COM

              krgd=0
              Do irgd=1,rigid%n_types
                 If (thermo%qr(irgd) == 1) Then
                    rgdtyp=rigid%list(0,irgd)
                    lrgd=rigid%list(-1,irgd)

                    Do jrgd=1,lrgd
                       krgd=krgd+1

                       i=rigid%index_local(jrgd,irgd) ! local index of particle/site

! COM distances

                       ggx(krgd)=parts(i)%xxx-rigid%xxx(irgd)
                       ggy(krgd)=parts(i)%yyy-rigid%yyy(irgd)
                       ggz(krgd)=parts(i)%zzz-rigid%zzz(irgd)
                    End Do
                 End If
              End Do

! minimum image convention for bond vectors

              Call images(imcon,cell,krgd,ggx,ggy,ggz)

! Get RB force and torque

              krgd=0
              Do irgd=1,rigid%n_types
                 If (thermo%qr(irgd) == 1) Then
                    rgdtyp=rigid%list(0,irgd)
                    lrgd=rigid%list(-1,irgd)

! calculate COM force and torque

                    fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
                    tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
                    Do jrgd=1,lrgd
                       krgd=krgd+1

                       i=rigid%index_local(jrgd,irgd) ! local index of particle/site

! If the RB has a frozen particle then no net force

                       If (rigid%frozen(0,rgdtyp) == 0) Then
                          fmx=fmx+parts(i)%fxx
                          fmy=fmy+parts(i)%fyy
                          fmz=fmz+parts(i)%fzz
                       End If

                       tqx=tqx+ggy(krgd)*parts(i)%fzz-ggz(krgd)*parts(i)%fyy
                       tqy=tqy+ggz(krgd)*parts(i)%fxx-ggx(krgd)*parts(i)%fzz
                       tqz=tqz+ggx(krgd)*parts(i)%fyy-ggy(krgd)*parts(i)%fxx
                    End Do

! If the RB has 2+ frozen particles (ill=1) the net torque
! must align along the axis of rotation

                    If (rigid%frozen(0,rgdtyp) > 1) Then
                       i1=rigid%index_local(rigid%index_global(1,rgdtyp),irgd)
                       i2=rigid%index_local(rigid%index_global(2,rgdtyp),irgd)

                       x(1)=parts(i1)%xxx-parts(i2)%xxx
                       y(1)=parts(i1)%yyy-parts(i2)%yyy
                       z(1)=parts(i1)%zzz-parts(i2)%zzz

                       Call images(imcon,cell,1,x,y,z)

                       tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)/(x(1)**2+y(1)**2+z(1)**2)
                       tqx=x(1)*tmp
                       tqy=y(1)*tmp
                       tqz=z(1)*tmp
                    End If

! current rotation matrix

                    Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

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

              Do irgd=1,rigid%n_types
                 If (thermo%qr(irgd) == 1) Then
                    rgdtyp=rigid%list(0,irgd)
                    lrgd=rigid%list(-1,irgd)

                    tmp=Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)

                    If (rigid%frozen(0,rgdtyp) == 0) Then
                       tkin  = tkin  + tmp * &
                       rigid%weight(0,rgdtyp)*(rigid%vxx(irgd)**2+ &
                       rigid%vyy(irgd)**2+rigid%vzz(irgd)**2)
                       vdotf = vdotf + tmp * &
                       (rigid%vxx(irgd)*rgdfxx(irgd)+ &
                       rigid%vyy(irgd)*rgdfyy(irgd)+ &
                       rigid%vzz(irgd)*rgdfzz(irgd))
                    End If

                    trot  = trot  + tmp * &
                    (rigid%rix(1,rgdtyp)*rigid%oxx(irgd)**2+ &
                    rigid%riy(1,rgdtyp)*rigid%oyy(irgd)**2+ &
                    rigid%riz(1,rgdtyp)*rigid%ozz(irgd)**2)
                    odott = odott + tmp * &
                    (rigid%oxx(irgd)*rgdtxx(irgd)+ &
                    rigid%oyy(irgd)*rgdtyy(irgd)+ &
                    rigid%ozz(irgd)*rgdtzz(irgd))
                 End If
              End Do

              Deallocate (ggx,ggy,ggz,          Stat=fail(1))
              Deallocate (rgdfxx,rgdfyy,rgdfzz, Stat=fail(2))
              Deallocate (rgdtxx,rgdtyy,rgdtzz, Stat=fail(3))
              If (Any(fail > 0)) Then
                 Write(message,'(a)') 'pseudo (RB) deallocation failure'
                 Call error(0,message)
              End If
           End If

! Shells

           If (thermo%stp > 0) Then
              If (thermo%tps(comm%idnode) > 0) Then
                 Do k=1,cshell%ntshl
                    If (thermo%qs(0,k) == 1) Then
                       i2=thermo%qs(2,k)

                       tkin  = tkin  + weight(i2)*(vxx(i2)**2+vyy(i2)**2+vzz(i2)**2)
                       vdotf = vdotf + vxx(i2)*parts(i2)%fxx+vyy(i2)*parts(i2)%fyy+vzz(i2)*parts(i2)%fzz
                    End If
                 End Do
              End If
           End If

              buffer(1) = tkin
              buffer(2) = vdotf
              buffer(3) = trot
              buffer(4) = odott
              Call gsum(comm,buffer(1:4))
              tkin  = buffer(1)
              vdotf = buffer(2)
              trot  = buffer(3)
              odott = buffer(4)

           tmp=tkin+trot
           If (tmp <= zero_plus) tmp = 1.0_wp
           thermo%chit_sb = (vdotf+odott)/tmp

        Else

! remove thermo%system centre of mass velocity (random momentum walk)

           Call getvom(vom,vxx,vyy,vzz,comm)

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
              If (thermo%qn(i) == 1) Then
                 tkin   = tkin  + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
                 vdotf  = vdotf + vxx(i)*parts(i)%fxx+vyy(i)*parts(i)%fyy+vzz(i)*parts(i)%fzz
              End If
           End Do

! Shells

           If (thermo%stp > 0) Then
              If (thermo%tps(comm%idnode) > 0) Then
                 Do k=1,cshell%ntshl
                    If (thermo%qs(0,k) == 1) Then
                       i2=thermo%qs(2,k)

                       tkin  = tkin  + weight(i2)*(vxx(i2)**2+vyy(i2)**2+vzz(i2)**2)
                       vdotf = vdotf + vxx(i2)*parts(i2)%fxx+vyy(i2)*parts(i2)%fyy+vzz(i2)*parts(i2)%fzz
                    End If
                 End Do
              End If
           End If

              buffer(1) = tkin
              buffer(2) = vdotf
              Call gsum(comm,buffer(1:2))
              tkin  = buffer(1)
              vdotf = buffer(2)

           If (tkin <= zero_plus) tkin = 1.0_wp
           thermo%chit_sb = vdotf/tkin

        End If

     End If

     If (thermo%key_pseudo == 0 .or. thermo%key_pseudo == 3) Then ! Apply DIRECT temperature scaling

! Targeted energy

        scale = boltz * thermo%temp_pseudo

! Scale velocity within the thermostat layer

        Do i=1,natms
           If (thermo%qn(i) == 1 .and. lfree(i) == 0) Then

! Get particle kinetic energy and produce a scaler to target temperature

              tmp = Sqrt(scale * dof_site(lsite(i)) / (weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)))

              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

        If (thermo%rtp > 0) Then

! calculate new RBs' COM and angular velocities

           Do irgd=1,rigid%n_types
              If (thermo%qr(irgd) == 1) Then
                 rgdtyp=rigid%list(0,irgd)

                 i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                 i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                 mxdr = 3.0_wp
                 If      (rigid%frozen(0,rgdtyp) == 0) Then
                    tmp = Sqrt( scale * mxdr / &
                    (rigid%weight(0,irgd)*(rigid%vxx(irgd)**2+rigid%vyy(irgd)**2+rigid%vzz(irgd)**2)) )
                    rigid%vxx(irgd) = rigid%vxx(irgd)*tmp
                    rigid%vyy(irgd) = rigid%vyy(irgd)*tmp
                    rigid%vzz(irgd) = rigid%vzz(irgd)*tmp
                 Else If (rigid%frozen(0,rgdtyp) == 1) Then
                    mxdr = 2.0_wp
                 Else If (rigid%frozen(0,rgdtyp) > 1) Then
                    mxdr = 1.0_wp
                 End If

                 tmp = Sqrt( scale * mxdr / &
                 (rigid%rix(1,rgdtyp)*rigid%oxx(irgd)**2+ &
                  rigid%riy(1,rgdtyp)*rigid%oyy(irgd)**2+ &
                  rigid%riz(1,rgdtyp)*rigid%ozz(irgd)**2) )
                 rigid%oxx(irgd) = rigid%oxx(irgd)*tmp
                 rigid%oyy(irgd) = rigid%oyy(irgd)*tmp
                 rigid%ozz(irgd) = rigid%ozz(irgd)*tmp

! get new rotation matrix

                 Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

! update RB members new velocities

                 lrgd=rigid%list(-1,irgd)
                 Do jrgd=1,lrgd
                    If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                       i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                       If (i <= natms) Then
                          x(1)=rigid%x(jrgd,rgdtyp)
                          y(1)=rigid%y(jrgd,rgdtyp)
                          z(1)=rigid%z(jrgd,rgdtyp)

! new atomic velocities in body frame

                          vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                          vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                          vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

! new atomic velocities in lab frame

                          vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                          vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                          vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                       End If
                    End If
                 End Do
              End If
           End Do

        End If

! Thermalise the shells on hit cores

        If (thermo%stp > 0) Then
           If (cshell%lshmv_shl) Then
             Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,vxx,vyy,vzz,domain,comm)
           End If

           If (thermo%tps(comm%idnode) > 0) Then
              Do k=1,cshell%ntshl
                 If (thermo%qs(0,k) == 1) Then
                    i1=thermo%qs(1,k)
                    i2=thermo%qs(2,k)

                    vxx(i2)=vxx(i1)
                    vyy(i2)=vyy(i1)
                    vzz(i2)=vzz(i1)
                 End If
              End Do
           End If
        End If

! remove thermo%system centre of mass velocity (random momentum walk)

        If (rigid%total > 0) Then

           Call getvom(vom,vxx,vyy,vzz,rigid,comm)

           Do j=1,nfree
              i=lstfre(j)

              If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                 vxx(i) = vxx(i) - vom(1)
                 vyy(i) = vyy(i) - vom(2)
                 vzz(i) = vzz(i) - vom(3)
              End If
           End Do

           Do irgd=1,rigid%n_types
              rgdtyp=rigid%list(0,irgd)

              If (rigid%frozen(0,rgdtyp) == 0) Then
                 rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
                 rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
                 rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

                 lrgd=rigid%list(-1,irgd)
                 Do jrgd=1,lrgd
                    i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                    If (i <= natms) Then
                       vxx(i) = vxx(i) - vom(1)
                       vyy(i) = vyy(i) - vom(2)
                       vzz(i) = vzz(i) - vom(3)
                    End If
                 End Do
              End If
           End Do

        Else

           Call getvom(vom,vxx,vyy,vzz,comm)

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

     If (rigid%total > 0) Then
        Call kinstresf(vxx,vyy,vzz,stats%strknf,comm)
        Call kinstrest(rigid,stats%strknt,comm)

        stats%strkin=stats%strknf+stats%strknt

! update rotational energy

        stats%engrot=getknr(rigid,comm)
     Else
        Call kinstress(vxx,vyy,vzz,stats%strkin,comm)
     End If
     stats%engke = 0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))

  End If

End Subroutine stochastic_boundary_vv
   

End Module stochastic_boundary
