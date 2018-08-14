Module drivers
  use kinds, Only : wp,wi
  Use comms, Only : comms_type,gsum
  use kinetics, Only : kinstresf, kinstrest, kinstress,getknr,getvom, getknr
  Use configuration, Only : vxx,vyy,vzz,cell,natms,&
                            weight,lfrzn,lstfre,lsite,lfree,ltg,nlast,nfree,&
                            lsi,lsa,imcon
  Use particle, Only : corePart
  Use rigid_bodies, Only : rigid_bodies_type,getrotmat
  Use setup, Only : boltz,mxatms,zero_plus
  Use angles, Only : angles_type
  Use dihedrals, Only : dihedrals_type
  Use inversions, Only : inversions_type
  Use core_shell,  Only : core_shell_type,SHELL_ADIABATIC 

  Use impacts, Only : impact_type, impact
  Use errors_warnings, Only : error,warning,info
  Use shared_units, Only : update_shared_units,update_shared_units_int
  Use numerics, Only : local_index,images,dcell,invert,box_mueller_saru3
  Use thermostat, Only : thermostat_type
  Use statistics, Only : stats_type
  Use domains, Only : domains_type
  Implicit None
  Private
  Public :: w_impact_option
  Public :: pseudo_vv

Contains

  Subroutine w_impact_option(levcfg,nstep,nsteql,rigid,cshell,stats,impa,comm)

    Integer( Kind = wi ),   Intent( InOut ) :: levcfg,nstep,nsteql
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type(stats_type)   ,   Intent( InOut ) :: stats
    Type(core_shell_type)   ,   Intent( InOut ) :: cshell
    Type(impact_type)   ,   Intent( InOut ) :: impa
    Type(comms_type)    ,   Intent( InOut ) :: comm

    Character( Len = 256 ) :: messages(6)

!!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!

! Apply impact
! levcfg == 2 avoids application twice when tmd happens at (re)start for VV

     If (nstep == impa%tmd .and. levcfg == 2) Then
       Write(messages(1),'(a)') ''
       Write(messages(2),'(a)') 'initiating IMPACT:'
       Write(messages(3),'(a,i10)') 'particle (index): ', impa%imd
       Write(messages(4),'(a,i10)') 'timestep (steps): ', impa%tmd
       Write(messages(5),'(a,1p,e12.5)') 'energy   (keV):   ', impa%emd
       Write(messages(6),'(a,1p,3e12.4)') 'v-r(x,y,z):       ', impa%vmx, impa%vmy, impa%vmz
       Call info(messages,6,.true.)

        If (nstep+1 <= nsteql) Call warning(380,Real(nsteql,wp),0.0_wp,0.0_wp)

        Call impact(rigid,cshell,impa,comm)

! Correct kinetic stress and energy

        If (rigid%total > 0) Then
           Call kinstresf(vxx,vyy,vzz,stats%strknf,comm)
           Call kinstrest(rigid,stats%strknt,comm)

           stats%strkin=stats%strknf+stats%strknt

           stats%engrot=getknr(rigid,comm)
        Else
           Call kinstress(vxx,vyy,vzz,stats%strkin,comm)
        End If
        stats%engke = 0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))
     End If

!!!!!!!!!!!!!!!!!!!!!  W_IMPACT_OPTION INCLUSION  !!!!!!!!!!!!!!!!!!!!!!
  End Subroutine w_impact_option

!  Subroutine w_refresh_mappings()
!    Include 'w_refresh_mappings.F90'
!  End Subroutine w_refresh_mappings
!
!  Subroutine w_at_start_vv()
!    Include 'w_at_start_vv.F90'
!  End Subroutine w_at_start_vv
!
!  Subroutine w_integrate_vv(isw)
!    Integer, Intent( In    ) :: isw ! used for vv stage control
!
!    Include 'w_integrate_vv.F90'
!  End Subroutine w_integrate_vv
!
!  Subroutine w_kinetic_options()
!    Include 'w_kinetic_options.F90'
!  End Subroutine w_kinetic_options
!
!  Subroutine w_statistics_report()
!    Include 'w_statistics_report.F90'
!  End Subroutine w_statistics_report
!
!  Subroutine w_write_options()
!    Include 'w_write_options.F90'
!  End Subroutine w_write_options
!
!  Subroutine w_refresh_output()
!    Include 'w_refresh_output.F90'
!  End Subroutine w_refresh_output
!
!  Subroutine w_md_vv()
!    Include 'w_md_vv.F90'
!  End Subroutine w_md_vv
!
!  Subroutine w_replay_history()
!    Logical,     Save :: newjb = .true.
!    Real( Kind = wp ) :: tmsh        ! tmst replacement
!    Integer           :: nstpe,nstph ! nstep replacements
!    Integer           :: exout       ! exit indicator for reading
!
!    Include 'w_replay_history.F90'
!  End Subroutine w_replay_history
!
!  Subroutine w_replay_historf()
!    Logical,     Save :: newjb = .true.
!    Real( Kind = wp ) :: tmsh        ! tmst replacement
!    Integer           :: nstpe,nstph ! nstep replacements
!    Integer           :: exout       ! exit indicator for reading
!
!    Include 'w_replay_historf.F90'
!  End Subroutine w_replay_historf

Subroutine pseudo_vv(isw,tstep,nstep,dof_site,cshell,stats,thermo,rigid,domain,parts,comm)

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
!           of the system (but are picked up in the pressure).
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
  Type( thermostat_type ), Intent( In    ) :: thermo
  Type( rigid_bodies_type ), Intent( InOut ) :: rigid
  Type( domains_type ), Intent( In    ) :: domain
  Type( comms_type), Intent( InOut ) :: comm
  Type (corePart ), Intent( InOut ) :: parts(:)

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ntp,stp,rtp
  Integer                 :: fail(1:3),matms, &
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

  Character ( Len = 256 ) :: message

  fail = 0

! set matms

  matms=nlast
  If (comm%mxnode == 1) matms=natms

  If (isw == 0) Then

! Random forces cycle - thermostat-system decoupling
! Recalculate the following for NPT and NST integration

     If (thermo%variable_cell .or. newjob) Then
        If (newjob) Then
           newjob = .false.

           Allocate (qn(1:mxatms),tpn(0:comm%mxnode-1),    Stat=fail(1))
           Allocate (qs(0:2,1:cshell%mxshl),tps(0:comm%mxnode-1), Stat=fail(2))
           Allocate (qr(1:rigid%max_rigid),tpr(0:comm%mxnode-1),     Stat=fail(3))
           If (Any(fail > 0)) Then
              Write(message,'(a)') 'pseudo (q. and tp.) allocation failure'
              Call error(0,message)
           End If
        End If

        Call invert(cell,rcell,celprp(10))
        Call dcell(cell,celprp)

! The origin of the Coordinate System is in the middle of the MD box (remember)
! Get the real coordinates of the close end edge of the thermostat as if the
! origin of the Coordinate System was in the left-most corner of the MD box

        ssx=thermo%width_pseudo*celprp(1)/celprp(7)
        ssy=thermo%width_pseudo*celprp(2)/celprp(8)
        ssz=thermo%width_pseudo*celprp(3)/celprp(9)

! 1. Get the boundary thermostat thicknesses in fractional coordinates
! 2. sx,sy,sz are intervalled as [-0.5,+0.5) {as by construction are (0,+0.25]}
! 3. Outline the edge beyond which a particle belongs to the thermostat
!    0.5*thicknesses[MD box] - thickness[boundary thermostat]

        sx=rcell(1)*ssx+rcell(4)*ssy+rcell(7)*ssz ; sx=sx-Anint(sx) ; sx=0.5_wp-sx
        sy=rcell(2)*ssx+rcell(5)*ssy+rcell(8)*ssz ; sy=sy-Anint(sy) ; sy=0.5_wp-sy
        sz=rcell(3)*ssx+rcell(6)*ssy+rcell(9)*ssz ; sz=sz-Anint(sz) ; sz=0.5_wp-sz
     End If

! qualify non-shell, non-frozen free particles (n) for a random kick
! qualify RBs (r) by them and derive related shells (s)

! tpn(comm%idnode) number of thermostatted particles on this node (comm%idnode)
! ntp - grand total of non-shell, non-frozen particles to thermostat

     qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
     qs(0:2,1:cshell%ntshl) = 0 ! unqualified core-shell unit with a local shell
     qr(1:rigid%n_types)     = 0 ! unqualified RB

     j = 0
     Do i=1,natms

! For all particles on this domain get how far they are
! from the origin of the MD box

        ssx=rcell(1)*parts(i)%xxx+rcell(4)*parts(i)%yyy+rcell(7)*parts(i)%zzz ; ssx=Abs(ssx-Anint(ssx))
        ssy=rcell(2)*parts(i)%xxx+rcell(5)*parts(i)%yyy+rcell(8)*parts(i)%zzz ; ssy=Abs(ssy-Anint(ssy))
        ssz=rcell(3)*parts(i)%xxx+rcell(6)*parts(i)%yyy+rcell(9)*parts(i)%zzz ; ssz=Abs(ssz-Anint(ssz))

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0 .and. &
            (ssx >= sx .or. ssy >= sy .or. ssz >= sz)) Then
           j = j + 1
           qn(i) = 1
        End If
     End Do
     tpn(comm%idnode) = j
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) tpn(i) = 0
        End Do
        Call gsum(comm,tpn)
     End If
     ntp = Sum(tpn)

     If (ntp == 0) Return
     If (chit < 1.0e-6_wp .or. thermo%key_pseudo > 1) Return ! Avoid thermostat overheating

! Allocate random force array of length j

     j=tpn(comm%idnode)
     Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(message,'(a)') 'pseudo (forces) allocation failure'
        Call error(0,message)
     End If

! Get gaussian distribution
! Get scaler to target variance*Sqrt(weight)

     scale = Sqrt(2.0_wp * chit * boltz * thermo%temp_pseudo / tstep)

     vom = 0.0_wp
     Do i=1,natms
        If (qn(i) == 1) Then

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
     vom = vom / Real(ntp,wp)

! Add random force and remove thermostat COM force

     Do i=1,natms
        If (qn(i) == 1) Then
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

     If (ntp == 0) Return

     j = 0 ! no qualified good RB (one qualified RB is enough to trigger all)
     Do i=1,matms
        If (qn(i) == 1) Then
           If (lfree(i) == 1) j = j + 1
        End If
     End Do
     Call gsum(comm,j)

! tpr(comm%idnode) number of thermostatted RB units on this node (comm%idnode)
! rtp - grand total of RB units to thermostat
! (can be larger than rigid%total due to sharing)

     k = 0
     If (j > 0) Then
        If (rigid%share) Then
           qn(natms+1:nlast) = 0 ! refresh the q array for shared RB units
           Call update_shared_units_int(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,qn,domain,comm)
        End If

        j = 0
        Do irgd=1,rigid%n_types
           rgdtyp=rigid%list(0,irgd)

! For all good RBs

           lrgd=rigid%list(-1,irgd)
           If (rigid%frozen(0,rgdtyp) < lrgd) Then
              Do jrgd=1,lrgd
                 i=rigid%index_local(jrgd,irgd) ! local index of particle/site
                 If (qn(i) == 1) Then
                    If (qr(irgd) == 0) Then ! An overall hit is registered
                       qr(irgd) = 1
                       j = j + 1
                    End If

                    If (i <= natms) tpn(comm%idnode) = tpn(comm%idnode) - 1 ! Less free particles are hit
                 End If
              End Do

              If (qr(irgd) == 1) Then ! accounting for a random kick on the RB
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
! tpn(comm%idnode) number of thermostatted free particles on this node (comm%idnode)
! ntp - grand total of non-shell, non-frozen free particles to thermostat
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) tpn(i) = 0
        End Do
        Call gsum(comm,tpn)
     End If
     ntp = Sum(tpn)
     tpr(comm%idnode) = j
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) tpr(i) = 0
        End Do
        Call gsum(comm,tpr)
     End If
     rtp = Sum(tpr)

! tps(comm%idnode) number of thermostatted core-shell units on this node (comm%idnode)
! stp - grand total of core-shell units to thermostat

     j = 0
     If (cshell%keyshl == SHELL_ADIABATIC) Then
        If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
           qn(natms+1:nlast) = 0
           Call update_shared_units_int(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,qn,domain,comm)
        End If

        If (cshell%ntshl > 0) Then
           Do k=1,cshell%ntshl
              i1=local_index(cshell%listshl(1,k),matms,lsi,lsa)
              i2=local_index(cshell%listshl(2,k),matms,lsi,lsa)

              If (qn(i1) == 1 .and. i2 > 0 .and. i2 <= natms) Then
                 j = j + 1

                 qs(0,k)=1
                 qs(1,k)=i1
                 qs(2,k)=i2
              End If
           End Do
        End If
     End If
     tps(comm%idnode) = j
     If (comm%mxnode > 1) Then
        Do i=0,comm%mxnode-1
           If (i /= comm%idnode) tps(i) = 0
        End Do
        Call gsum(comm,tps)
     End If
     stp = Sum(tps)

! Velocity scaling cycle - thermostatting.  k = local, ntp = global
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
! the number of nodes and system size
!
! Get gaussian distribution

        tkin = 0.0_wp
        mxdr = 0.0_wp
        Do i=1,natms
           If (qn(i) == 1 .and. lfree(i) == 0) Then
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

        If (rtp > 0) Then
           Do irgd=1,rigid%n_types
              If (qr(irgd) == 1) Then
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
           If (qn(i) == 1 .and. lfree(i) == 0) Then
              tmp = scale * Sqrt( (xxt(i)**2+yyt(i)**2+zzt(i)**2) / &
                                  (vxx(i)**2+vyy(i)**2+vzz(i)**2) )

              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

        If (rtp > 0) Then

! Update shared RBs' velocities

           If (rigid%share) Then
             Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared,rigid%map_shared,xxt,yyt,zzt,domain,comm)
           End If

! calculate new RBs' COM and angular velocities

           Do irgd=1,rigid%n_types
              If (qr(irgd) == 1) Then
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

        If (stp > 0) Then
           If (cshell%lshmv_shl) Then
             Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,vxx,vyy,vzz,domain,comm)
           End If

           If (tps(comm%idnode) > 0) Then
              j = 0
              Do k=1,cshell%ntshl
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

        If (rigid%total > 0) Then

! remove system centre of mass velocity (random momentum walk)

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

              If (qn(i) == 1) Then
                 tkin  = tkin  + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
                 vdotf = vdotf + vxx(i)*parts(i)%fxx+vyy(i)*parts(i)%fyy+vzz(i)*parts(i)%fzz
              End If
           End Do

! RBs

           If (rtp > 0) Then
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
                 If (qr(irgd) == 1) Then
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
                 If (qr(irgd) == 1) Then
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
                 If (qr(irgd) == 1) Then
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

           If (stp > 0) Then
              If (tps(comm%idnode) > 0) Then
                 Do k=1,cshell%ntshl
                    If (qs(0,k) == 1) Then
                       i2=qs(2,k)

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
           chit = (vdotf+odott)/tmp

        Else

! remove system centre of mass velocity (random momentum walk)

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
              If (qn(i) == 1) Then
                 tkin   = tkin  + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
                 vdotf  = vdotf + vxx(i)*parts(i)%fxx+vyy(i)*parts(i)%fyy+vzz(i)*parts(i)%fzz
              End If
           End Do

! Shells

           If (stp > 0) Then
              If (tps(comm%idnode) > 0) Then
                 Do k=1,cshell%ntshl
                    If (qs(0,k) == 1) Then
                       i2=qs(2,k)

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
           chit = vdotf/tkin

        End If

     End If

     If (thermo%key_pseudo == 0 .or. thermo%key_pseudo == 3) Then ! Apply DIRECT temperature scaling

! Targeted energy

        scale = boltz * thermo%temp_pseudo

! Scale velocity within the thermostat layer

        Do i=1,natms
           If (qn(i) == 1 .and. lfree(i) == 0) Then

! Get particle kinetic energy and produce a scaler to target temperature

              tmp = Sqrt(scale * dof_site(lsite(i)) / (weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)))

              vxx(i) = vxx(i)*tmp
              vyy(i) = vyy(i)*tmp
              vzz(i) = vzz(i)*tmp
           End If
        End Do

        If (rtp > 0) Then

! calculate new RBs' COM and angular velocities

           Do irgd=1,rigid%n_types
              If (qr(irgd) == 1) Then
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

        If (stp > 0) Then
           If (cshell%lshmv_shl) Then
             Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,vxx,vyy,vzz,domain,comm)
           End If

           If (tps(comm%idnode) > 0) Then
              Do k=1,cshell%ntshl
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

! remove system centre of mass velocity (random momentum walk)

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

End Subroutine pseudo_vv
End Module drivers
