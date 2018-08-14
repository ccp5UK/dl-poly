Module temperature
  Use kinds,           Only : wp, li
  Use comms,           Only : comms_type,gsum
  Use setup,           Only : nrite,boltz,mxatms,zero_plus
  Use site, Only : site_type
  Use configuration,   Only : imcon,natms,nlast,nfree,lsite,  &
                              lsi,lsa,ltg,lfrzn,lfree,lstfre, &
                              weight,vxx,vyy,vzz
  Use particle,        Only : corePart
  Use rigid_bodies,    Only : rigid_bodies_type,getrotmat,rigid_bodies_quench
  Use constraints,     Only : constraints_type, constraints_quench
  Use pmf,             Only : pmf_quench, pmf_type
  Use core_shell,      Only : core_shell_type,core_shell_quench,SHELL_ADIABATIC
  Use kinetics,        Only : l_vom,chvom,getcom,getvom,getkin,getknf,getknt,getknr
  Use numerics,        Only : invert,uni,local_index,box_mueller_saru3
  use shared_units,    Only : update_shared_units,update_shared_units_int
  Use errors_warnings, Only : error,warning,info
  Use thermostat,      Only : thermostat_type
  Use statistics, Only : stats_type
  Use minimise, Only : minimise_type
  Use domains, Only : domains_type
  Implicit None

  Private

  Public :: set_temperature, regauss_temperature, scale_temperature

Contains

  Subroutine set_temperature(levcfg,keyres,nstep,nstrun,atmfre,atmfrz,degtra, &
      degrot,degfre,degshl,engrot,dof_site,cshell,stat,cons,pmf,thermo,minim, &
      rigid,domain,parts,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for setting the initial system temperature
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov january 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: nstep,nstrun, &
                                           atmfre,atmfrz
    Integer,            Intent( InOut ) :: keyres,levcfg
    Integer(Kind=li),   Intent( InOut ) :: degtra,degrot
    Integer(Kind=li),   Intent(   Out ) :: degfre,degshl
    Real( Kind = wp ),  Intent(   Out ) :: engrot
    Real( Kind = wp ), Dimension(:), Intent( In    ) :: dof_site
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( stats_type ), Intent( InOut ) :: stat
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( minimise_type ), Intent( In    ) :: minim
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( corePart ),   Intent( InOut ) :: parts(:)
    Type( comms_type ), Intent( InOut ) :: comm

    Logical           :: no_min_0,safe
    Integer           :: fail(1:2),i,j,k,ntp,   &
                         stp,i1,i2, &
                         irgd,jrgd,lrgd,rgdtyp
    Integer(Kind=li)  :: com,con,frz,meg,non
    Real( Kind = wp ) :: tmp,vom(1:3),engk,engf,engt,engr, &
                         x(1:1),y(1:1),z(1:1),rot(1:9),    &
                         vpx,vpy,vpz


  ! q. index arrays and tp. sum arrays

    Integer,           Allocatable :: qn(:),tpn(:)
    Integer,           Allocatable :: qs(:,:),tps(:)

    Character ( Len = 256 )  ::  message,messages(10)

  ! initialise rotational and translational DoF if no RB are present
  ! or re-initialise if all are frozen (does no harm)

    If (rigid%total == 0) Then
       degtra=Int(0,li)
       degrot=Int(0,li)
    End If

  ! Degrees Of Freedom (DoF) - all free particles

    meg=Int(3,li)*Int(atmfre,li)

  ! lost to frozen free atoms

    frz=Int(3,li)*Int(atmfrz,li)

  ! 3 lost for fixing COM translation

    If (l_vom .and. thermo%key_dpd == 0) Then
       com=Int(3,li)
    Else
       com=Int(0,li)
    End If

  ! 3 lost for fixing angular momentum about origin
  ! (non-periodic systems only)

    non=Int(0,li)
    If (imcon == 0) non=Int(3,li)

  ! lost to shells

    degshl=Int(3,li)*Int(cshell%megshl,li)

  ! lost to constrained atoms and PMF constraints

    con=Int(cons%megcon)+Int(pmf%megpmf)

  ! TOTAL DoF

    degfre = meg - com - non - frz - degshl - con + degrot + degtra

  ! Report DoF

    Write(messages(1),'(a)') 'degrees of freedom break-down neigh%list:'
    Write(messages(2),'(2x,a,i12)') 'free particles        ',meg
    Write(messages(3),'(2x,a,i12)') 'centre of mass        ',-com
    Write(messages(4),'(2x,a,i12)') 'non-periodicity       ',-non
    Write(messages(5),'(2x,a,i12)') 'frozen free particles ',-frz
    Write(messages(6),'(2x,a,i12)') 'shell-pseudo          ',-degshl
    Write(messages(7),'(2x,a,i12)') 'constrained           ',-con
    Write(messages(8),'(2x,a,i12)') 'RB translational      ',degtra
    Write(messages(9),'(2x,a,i12)') 'RB rotational         ',degrot
    Write(messages(10),'(2x,a,i12)') 'total (real)          ',degfre
    Call info(messages,10,.true.)

  ! Check DoF distribution

    tmp=0.0_wp
    Do i=1,natms
       If (dof_site(lsite(i)) > zero_plus) & ! Omit shells' negative DoFs
          tmp=tmp+dof_site(lsite(i))
    End Do
    Call gsum(comm,tmp)
    If (Nint(tmp,li)-non-com /= degfre) Call error(360)

  ! warn for much restrain on the system [ degfre <= (com+non+frz+con) ]
  ! and catch an over-restrained system

    If (degfre <= Int(com+non+frz+con,li)) &
       Call warning(210,Real(degfre,wp),Real((com+non+frz+con),wp),0.0_wp)
    If (degfre < Int(1,li)) Call error(350)

  ! desired kinetic energy

    thermo%sigma=0.5_wp*Real(degfre,wp)*boltz*thermo%temp

  ! avoid user defined 0K field to break up anything

    If (rigid%total > 0) Then
       engf=getknf(vxx,vyy,vzz,comm)/Real(Max(1_li,degfre),wp)
       engt=getknt(rigid,comm)/Real(Max(1_li,degtra),wp)

       engr=getknr(rigid,comm)/Real(Max(1_li,degrot),wp)
       engk=engf+engt
    Else
       engk=getkin(vxx,vyy,vzz,comm)/Real(Max(1_li,degfre),wp)
    End If
    If (thermo%sigma > 1.0e-6_wp .and. engk < 1.0e-6_wp .and. (keyres /= 0 .and. nstrun /= 0)) Then
       Call warning('0K velocity field detected in CONFIG with a restart at non 0K temperature in CONTROL',.true.)
       Call info('*** clean start enforced ***',.true.)

       keyres = 0
    End If

    If (keyres == 0) Then

       Allocate (qn(1:mxatms),tpn(0:comm%mxnode-1),    Stat=fail(1))
       Allocate (qs(0:2,1:cshell%mxshl),tps(0:comm%mxnode-1), Stat=fail(2))
       If (Any(fail > 0)) Then
          Write(message,'(a)') 'set_temperature allocation failure'
          Call error(0,message)
       End If

  ! tpn(idnode) number of particles on this node (idnode)
  ! ntp - grand total of non-shell, non-frozen particles

       qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
       qs(0:2,1:cshell%ntshl) = 0 ! unqualified core-shell unit with a local shell

       j = 0
       Do i=1,natms
          If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0) Then
             j = j + 1
             qn(i) = 1
          End If
       End Do
       tpn(comm%idnode) = j
       Do i=0,comm%mxnode-1
          If (i /= comm%idnode) tpn(i) = 0
       End Do
       Call gsum(comm,tpn)
       ntp = Sum(tpn)

  ! Save core-shell internal energy vectors due to possible hits
  ! tps(idnode) number of thermostatted core-shell units on this node (idnode)
  ! stp - grand total of core-shell units to thermostat

       j = 0
       If (cshell%keyshl == SHELL_ADIABATIC) Then ! just for the adiabatic shell model
          If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
             qn(natms+1:nlast) = 0
             Call update_shared_units_int(natms,nlast,lsi,lsa,cshell%lishp_shl, &
               cshell%lashp_shl,qn,domain,comm)
             Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl, &
               cshell%lashp_shl,vxx,vyy,vzz,domain,comm)
          End If

          If (cshell%ntshl > 0) Then
             Do k=1,cshell%ntshl
                i1=local_index(cshell%listshl(1,k),nlast,lsi,lsa)
                i2=local_index(cshell%listshl(2,k),nlast,lsi,lsa)

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
       Do i=0,comm%mxnode-1
          If (i /= comm%idnode) tps(i) = 0
       End Do
       Call gsum(comm,tps)
       stp = Sum(tps)

       If (rigid%total > 0) Then

          k = 0
          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then
                Do jrgd=1,lrgd
                   i=rigid%index_local(jrgd,irgd) ! local index of particle/site
                   If (i <= natms) tpn(comm%idnode) = tpn(comm%idnode) - 1 ! Less free particles are hit
                End Do

                i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum
                If (rigid%frozen(0,rgdtyp) == 0) Then
                   If (i1 <= natms) k = k + 1
                End If
                If (i2 <= natms) k = k + 1
             End If
          End Do
  ! tpn(idnode) number of thermostatted free particles on this node (idnode)
  ! ntp - grand total of non-shell, non-frozen free particles to thermostat
          Do i=0,comm%mxnode-1
             If (i /= comm%idnode) tpn(i) = 0
          End Do
          Call gsum(comm,tpn)
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

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then
                i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                If (rigid%frozen(0,rgdtyp) == 0 .and. i1 <= natms) Then
                   Call box_mueller_saru3(ltg(i1),0,vxx(i1),vyy(i1),vzz(i1))

  ! Get scaler to target variance/Sqrt(weight)

                   tmp = 1.0_wp/Sqrt(rigid%weight(0,rgdtyp))
                   vxx(i1) = vxx(i1)*tmp
                   vyy(i1) = vyy(i1)*tmp
                   vzz(i1) = vzz(i1)*tmp
                End If

                If (i2 <= natms) Then
                   Call box_mueller_saru3(ltg(i2),0,vxx(i2),vyy(i2),vzz(i2))

  ! Get scaler to target variance/Sqrt(weight) -
  ! 3 different reciprocal moments of inertia

                   vxx(i2) = vxx(i2)*Sqrt(rigid%rix(2,rgdtyp))
                   vyy(i2) = vyy(i2)*Sqrt(rigid%riy(2,rgdtyp))
                   vzz(i2) = vzz(i2)*Sqrt(rigid%riz(2,rgdtyp))
                End If
             End If
          End Do

  ! Update shared RBs' velocities

          If (rigid%share) Then
            Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared, &
              rigid%map_shared,vxx,vyy,vzz,domain,comm)
          End If

  ! calculate new RBs' COM and angular velocities

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then
                i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                If (rigid%frozen(0,rgdtyp) == 0) Then
                   rigid%vxx(irgd) = vxx(i1)
                   rigid%vyy(irgd) = vyy(i1)
                   rigid%vzz(irgd) = vzz(i1)
                End If

                rigid%oxx(irgd) = vxx(i2)
                rigid%oyy(irgd) = vyy(i2)
                rigid%ozz(irgd) = vzz(i2)
                If (i2 <= natms) Then
                   If (lfrzn(i2) > 0 .or. weight(i) < 1.0e-6_wp) Then
                      vxx(i2) = 0.0_wp
                      vyy(i2) = 0.0_wp
                      vzz(i2) = 0.0_wp
                   End If
                End If

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
          If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
             qn(natms+1:nlast) = 0
             Call update_shared_units_int(natms,nlast,lsi,lsa,cshell%lishp_shl, &
               cshell%lashp_shl,qn,domain,comm)
             Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl, &
               cshell%lashp_shl,vxx,vyy,vzz,domain,comm)
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

       Deallocate (qn,tpn, Stat=fail(1))
       Deallocate (qs,tps, Stat=fail(2))
       If (Any(fail > 0)) Then
          Write(message,'(a)') 'set_temperature deallocation failure'
          Call error(0,message)
       End If

  ! remove centre of mass motion

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

    Else

  ! quench RBs

       If (rigid%total > 0) Then
         Call rigid_bodies_quench(rigid,domain,parts,comm)
       End If

    End If

  ! levcfg must be equalised to ONE (velocities now exist)

    levcfg=1

  ! scale velocities for keyres=0 and keyres=2

    If (keyres == 0 .or. keyres == 2) Then

  ! Detect pure molecular statics only == CGM minimisation at zero timestep
  ! and no rescaling

       no_min_0 = .not.(minim%minimise .and. minim%freq == 0 .and. nstep == 0 .and. nstrun == 0 .and. keyres == 2)

  ! quench constraints & PMFs

       If (no_min_0) Then
          If (cons%megcon > 0) Then
            Call constraints_quench(cons,stat,domain,parts,comm)
          End If
          If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,parts,comm)
       End If

  ! quench core-shell units in adiabatic model

       If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC .and. no_min_0) Then
          Do
             Call scale_temperature(thermo%sigma,degtra,degrot,degfre,rigid,parts,comm)
             Call core_shell_quench(safe,thermo%temp,cshell,domain,comm)
             If (cons%megcon > 0) Then
               Call constraints_quench(cons,stat,domain,parts,comm)
             End If
             If (pmf%megpmf > 0) Then
               Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,parts,comm)
             End If
             If (rigid%total > 0) Then
               Call rigid_bodies_quench(rigid,domain,parts,comm)
             End If
             If (safe) Exit
          End Do
       Else
          Call scale_temperature(thermo%sigma,degtra,degrot,degfre,rigid,parts,comm)
       End If

    End If

    If (l_vom) Then

  ! remove centre of mass motion

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
    Else                 ! make getvom always return 0, which is not good for
       Call chvom(l_vom) ! standard MD as the flying ice-cub effect may happen
    End If               ! and/or T(MD) is loses its microscopic meaning!

  ! Initialise engrot and if RBs exist calculate it

    engrot=0.0_wp
    If (rigid%total > 0) engrot=getknr(rigid,comm)

  End Subroutine set_temperature

  Subroutine regauss_temperature(rigid,domain,parts,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to regauss the instantaneous system temperature
  ! by random pairwise swaps of the energy scaled momenta of dynamically
  ! active particles (no massless shells or massless RB members)
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( corePart ),   Intent( InOut ) :: parts(:)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: fail,i,j,k,l,is,irgd,jrgd,lrgd,rgdtyp
    Real( Kind = wp ) :: vom(1:3),tmp

    Integer, Allocatable :: ind(:),pair(:,:)

    Character ( Len = 256 )  :: message

    fail=0
    Allocate (ind(1:natms),pair(1:2,natms/2), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'regauss_temperature allocation failure'
       Call error(0,message)
    End If

  ! Create and index array containing the indices of the
  ! dynamically active particles and zeros for the inactive

    Do i=1,natms
       If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
          ind(i)=i
       Else
          ind(i)=0
       End If
    End Do

  ! Compress the index array

    i=1
    j=natms
    Do While (i < j)
       Do While (ind(j) == 0 .and. j > i)
          j=j-1
       End Do

       If (i < j) Then
          If (ind(i) == 0) Then
             ind(i)=ind(j)
             ind(j)=0
             j=j-1
          Else
             i=i+1
          End If
       End If
    End Do
    k=j

  ! Create a non-overlapping array of random pairs
  ! by exhausting the index array

    Do i=1,k/2
       Do l=1,2
          is=1+Int(Real(j,wp)*uni(comm))
          pair(l,i)=ind(is)
          ind(is)=ind(j)
          ind(j)=0
          j=j-1
       End Do
    End Do

  ! Swap particles energies in the pair array

    Do i=1,k/2
       j=pair(1,i)
       l=pair(2,i)

       vom(1)=vxx(j)
       vom(2)=vyy(j)
       vom(3)=vzz(j)

       tmp=Sqrt(weight(l)/weight(j))
       vxx(j)=vxx(l)*tmp
       vyy(j)=vyy(l)*tmp
       vzz(j)=vzz(l)*tmp

       tmp=1.0_wp/tmp
       vxx(l)=vom(1)*tmp
       vyy(l)=vom(2)*tmp
       vzz(l)=vom(3)*tmp
    End Do

    If (rigid%total > 0) Then

  ! quench RBs

      Call rigid_bodies_quench(rigid,domain,parts,comm)

  ! remove centre of mass motion

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

  ! remove centre of mass motion

       Call getvom(vom,vxx,vyy,vzz,comm)

       Do i=1,natms
          If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
             vxx(i) = vxx(i) - vom(1)
             vyy(i) = vyy(i) - vom(2)
             vzz(i) = vzz(i) - vom(3)
          End If
       End Do

    End If

    Deallocate (ind,pair, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'regauss_temperature deallocation failure'
       Call error(0,message)
    End If

  End Subroutine regauss_temperature

  Subroutine scale_temperature(sigma,degtra,degrot,degfre,rigid,parts,comm)

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

    Real( Kind = wp ),   Intent( In    ) :: sigma
    Integer( Kind=li ),  Intent( In    ) :: degtra,degrot,degfre
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( corePart ),    Intent( InOut ) :: parts(:)
    Type ( comms_type ), Intent( InOut ) :: comm


    Integer           :: fail,i,j,irgd,jrgd,lrgd,rgdtyp
    Real( Kind = wp ) :: engke,engkf,engkt,engrot,         &
                         amx,amy,amz,wxx,wyy,wzz,tmp,tmp1, &
                         com(1:3),vom(1:3),rot(1:9),rotinv(1:9),x,y,z

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    Character ( Len = 256 )  ::  message

  ! remove centre of mass motion

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

  ! zero angular momentum about centre of mass - non-periodic system

    If (imcon == 0) Then
       fail=0
       Allocate (buffer(1:12), Stat=fail)
       If (fail > 0) Then
          Write(message,'(a)') 'scale_temperature allocation failure'
          Call error(0,message)
       End If

  ! calculate centre of mass position

       Call getcom(parts,com,comm)

       If (rigid%total > 0) Then

  ! move to centre of mass origin

          Do j=1,nfree
             i=lstfre(j)

             If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                parts(i)%xxx = parts(i)%xxx - com(1)
                parts(i)%yyy = parts(i)%yyy - com(2)
                parts(i)%zzz = parts(i)%zzz - com(3)
             End If
          End Do

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             If (rigid%frozen(0,rgdtyp) == 0) Then
                rigid%xxx(irgd) = rigid%xxx(irgd) - com(1)
                rigid%yyy(irgd) = rigid%yyy(irgd) - com(2)
                rigid%zzz(irgd) = rigid%zzz(irgd) - com(3)
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
                amx = amx + weight(i)*(parts(i)%yyy*vzz(i) - parts(i)%zzz*vyy(i))
                amy = amy + weight(i)*(parts(i)%zzz*vxx(i) - parts(i)%xxx*vzz(i))
                amz = amz + weight(i)*(parts(i)%xxx*vyy(i) - parts(i)%yyy*vxx(i))

                tmp = parts(i)%xxx**2 + parts(i)%yyy**2 + parts(i)%zzz**2
                rot(1) = rot(1) + weight(i)*(parts(i)%xxx*parts(i)%xxx - tmp)
                rot(2) = rot(2) + weight(i)* parts(i)%xxx*parts(i)%yyy
                rot(3) = rot(3) + weight(i)* parts(i)%xxx*parts(i)%zzz
                rot(5) = rot(5) + weight(i)*(parts(i)%yyy*parts(i)%yyy - tmp)
                rot(6) = rot(6) + weight(i)* parts(i)%yyy*parts(i)%zzz
                rot(9) = rot(9) + weight(i)*(parts(i)%zzz*parts(i)%zzz - tmp)
             End If
          End Do

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             If (rigid%frozen(0,rgdtyp) == 0) Then
                lrgd=rigid%list(-1,irgd)

                tmp1=rigid%weight(0,rgdtyp)*Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)

                amx = amx + tmp1*(rigid%yyy(irgd)*rigid%vzz(irgd) - rigid%zzz(irgd)*rigid%vyy(irgd))
                amy = amy + tmp1*(rigid%zzz(irgd)*rigid%vxx(irgd) - rigid%xxx(irgd)*rigid%vzz(irgd))
                amz = amz + tmp1*(rigid%xxx(irgd)*rigid%vyy(irgd) - rigid%yyy(irgd)*rigid%vxx(irgd))

                tmp = rigid%xxx(irgd)**2 + rigid%yyy(irgd)**2 + rigid%zzz(irgd)**2

                rot(1) = rot(1) + tmp1*(rigid%xxx(irgd)*rigid%xxx(irgd) - tmp)
                rot(2) = rot(2) + tmp1* rigid%xxx(irgd)*rigid%yyy(irgd)
                rot(3) = rot(3) + tmp1* rigid%xxx(irgd)*rigid%zzz(irgd)
                rot(5) = rot(5) + tmp1*(rigid%yyy(irgd)*rigid%yyy(irgd) - tmp)
                rot(6) = rot(6) + tmp1* rigid%yyy(irgd)*rigid%zzz(irgd)
                rot(9) = rot(9) + tmp1*(rigid%zzz(irgd)*rigid%zzz(irgd) - tmp)
             End If
          End Do

  ! complete rotational inertia matrix

          rot(4) = rot(2)
          rot(7) = rot(3)
          rot(8) = rot(6)

  ! global sum of rotation

          buffer(1) = amx
          buffer(2) = amy
          buffer(3) = amz
          Do i=1,9
             buffer(i+3) = rot(i)
          End Do

          Call gsum(comm,buffer)

          amx =  buffer(1)
          amy =  buffer(2)
          amz =  buffer(3)
          Do i=1,9
             rot(i) = buffer(i+3)
          End Do

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
                vxx(i) = vxx(i) + (wyy*parts(i)%zzz - wzz*parts(i)%yyy)
                vyy(i) = vyy(i) + (wzz*parts(i)%xxx - wxx*parts(i)%zzz)
                vzz(i) = vzz(i) + (wxx*parts(i)%yyy - wyy*parts(i)%xxx)
             End If
          End Do

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             If (rigid%frozen(0,rgdtyp) == 0) Then
                x=(wyy*rigid%zzz(irgd) - wzz*rigid%yyy(irgd))
                y=(wzz*rigid%xxx(irgd) - wxx*rigid%zzz(irgd))
                z=(wxx*rigid%yyy(irgd) - wyy*rigid%xxx(irgd))

                rigid%vxx(irgd) = rigid%vxx(irgd) + x
                rigid%vyy(irgd) = rigid%vyy(irgd) + y
                rigid%vzz(irgd) = rigid%vzz(irgd) + z

                lrgd=rigid%list(-1,irgd)
                Do jrgd=1,lrgd
                   i=rigid%index_local(jrgd,irgd) ! local index of particle/site

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
                parts(i)%xxx = parts(i)%xxx + com(1)
                parts(i)%yyy = parts(i)%yyy + com(2)
                parts(i)%zzz = parts(i)%zzz + com(3)
             End If
          End Do

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             If (rigid%frozen(0,rgdtyp) == 0) Then
                rigid%xxx(irgd) = rigid%xxx(irgd) + com(1)
                rigid%yyy(irgd) = rigid%yyy(irgd) + com(2)
                rigid%zzz(irgd) = rigid%zzz(irgd) + com(3)
             End If
          End Do

       Else

  ! move to centre of mass origin

          Do i=1,natms
             If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                parts(i)%xxx = parts(i)%xxx - com(1)
                parts(i)%yyy = parts(i)%yyy - com(2)
                parts(i)%zzz = parts(i)%zzz - com(3)
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
                amx = amx + weight(i)*(parts(i)%yyy*vzz(i) - parts(i)%zzz*vyy(i))
                amy = amy + weight(i)*(parts(i)%zzz*vxx(i) - parts(i)%xxx*vzz(i))
                amz = amz + weight(i)*(parts(i)%xxx*vyy(i) - parts(i)%yyy*vxx(i))

                tmp = parts(i)%xxx**2 + parts(i)%yyy**2 + parts(i)%zzz**2
                rot(1) = rot(1) + weight(i)*(parts(i)%xxx*parts(i)%xxx - tmp)
                rot(2) = rot(2) + weight(i)* parts(i)%xxx*parts(i)%yyy
                rot(3) = rot(3) + weight(i)* parts(i)%xxx*parts(i)%zzz
                rot(5) = rot(5) + weight(i)*(parts(i)%yyy*parts(i)%yyy - tmp)
                rot(6) = rot(6) + weight(i)* parts(i)%yyy*parts(i)%zzz
                rot(9) = rot(9) + weight(i)*(parts(i)%zzz*parts(i)%zzz - tmp)
             End If
          End Do

  ! complete rotational inertia matrix

          rot(4) = rot(2)
          rot(7) = rot(3)
          rot(8) = rot(6)

  ! global sum of rotation

          buffer(1) = amx
          buffer(2) = amy
          buffer(3) = amz
          Do i=1,9
             buffer(i+3) = rot(i)
          End Do

          Call gsum(comm,buffer)

          amx =  buffer(1)
          amy =  buffer(2)
          amz =  buffer(3)
          Do i=1,9
             rot(i) = buffer(i+3)
          End Do

  ! invert rotational inertia matrix

          Call invert(rot,rotinv,tmp)

  ! correction to angular velocity

          wxx = rotinv(1)*amx + rotinv(2)*amy + rotinv(3)*amz
          wyy = rotinv(4)*amx + rotinv(5)*amy + rotinv(6)*amz
          wzz = rotinv(7)*amx + rotinv(8)*amy + rotinv(9)*amz

  ! correction to linear velocity

          Do i=1,natms
             If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                vxx(i) = vxx(i) + (wyy*parts(i)%zzz - wzz*parts(i)%yyy)
                vyy(i) = vyy(i) + (wzz*parts(i)%xxx - wxx*parts(i)%zzz)
                vzz(i) = vzz(i) + (wxx*parts(i)%yyy - wyy*parts(i)%xxx)
             End If
          End Do

  ! reset positions to original reference frame

          Do i=1,natms
             If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                parts(i)%xxx = parts(i)%xxx + com(1)
                parts(i)%yyy = parts(i)%yyy + com(2)
                parts(i)%zzz = parts(i)%zzz + com(3)
             End If
          End Do

       End If

       Deallocate (buffer, Stat=fail)
       If (fail > 0) Then
          Write(message,'(a)') 'scale_temperature deallocation failure'
          Call error(0,message)
       End If
    End If

  ! ensure equipartitioning - all degrees of freedom are equal
  ! calculate energy: free particles and RB translational and rotational

    engrot=0.0_wp
    If (rigid%total > 0) Then
       engkf=getknf(vxx,vyy,vzz,comm)
       engkt=getknt(rigid,comm)

       engrot=getknr(rigid,comm)

  ! temporary replacement for small engrot

       tmp1=Max(engrot,1.0e-6_wp)

  ! Scale rotational energy to translational energy
  ! according to their respective DoF

       If (degtra > Int(0,li)) Then ! engkt > 0 (degrot > 0 and tmp1(engrot) > 0)
          tmp=Sqrt((engkt*Real(degrot,wp))/(tmp1*Real(degtra,wp)))

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then
                rigid%oxx(irgd)=rigid%oxx(irgd)*tmp
                rigid%oyy(irgd)=rigid%oyy(irgd)*tmp
                rigid%ozz(irgd)=rigid%ozz(irgd)*tmp
             End If
          End Do
          engrot=engrot*tmp**2
          tmp1=Max(engrot,1.0e-6_wp)
       End If

  ! Scale the energy per DoF of the RBs to that of a DoF of the free particles

       If (degfre-degtra-degrot > Int(0,li)) Then ! engkf > 0
          tmp=Sqrt((engkf*Real(degrot,wp))/(tmp1*Real(degfre-degtra-degrot,wp)))

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then
                rigid%oxx(irgd)=rigid%oxx(irgd)*tmp
                rigid%oyy(irgd)=rigid%oyy(irgd)*tmp
                rigid%ozz(irgd)=rigid%ozz(irgd)*tmp

                If (rigid%frozen(0,rgdtyp) == 0) Then
                   rigid%vxx(irgd)=rigid%vxx(irgd)*tmp
                   rigid%vyy(irgd)=rigid%vyy(irgd)*tmp
                   rigid%vzz(irgd)=rigid%vzz(irgd)*tmp
                End If
             End If
          End Do
          engrot=engrot*tmp**2
          If (degtra > Int(0,li)) engkt=engkt*tmp**2
       End If

       engke=engkf+engkt
    Else
       engke=getkin(vxx,vyy,vzz,comm)
    End If

  ! apply temperature scaling

    If (engke+engrot > 1.0e-6_wp .and. sigma > zero_plus) Then
       tmp=Sqrt(sigma/(engke+engrot))

       If (rigid%total > 0) Then
          Do j=1,nfree
             i=lstfre(j)

             If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
                vxx(i)=vxx(i)*tmp
                vyy(i)=vyy(i)*tmp
                vzz(i)=vzz(i)*tmp
             End If
          End Do

          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then

  ! new angular velocity

                rigid%oxx(irgd)=rigid%oxx(irgd)*tmp
                rigid%oyy(irgd)=rigid%oyy(irgd)*tmp
                rigid%ozz(irgd)=rigid%ozz(irgd)*tmp

  ! new translational velocity

                If (rigid%frozen(0,rgdtyp) == 0) Then
                   rigid%vxx(irgd)=rigid%vxx(irgd)*tmp
                   rigid%vyy(irgd)=rigid%vyy(irgd)*tmp
                   rigid%vzz(irgd)=rigid%vzz(irgd)*tmp
                End If

  ! new rotational matrix

                Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

                Do jrgd=1,lrgd
                   If (rigid%frozen(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                      i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                      If (i <= natms) Then
                         x=rigid%x(jrgd,rgdtyp)
                         y=rigid%y(jrgd,rgdtyp)
                         z=rigid%z(jrgd,rgdtyp)

  ! site velocity in body frame

                         wxx=rigid%oyy(irgd)*z-rigid%ozz(irgd)*y
                         wyy=rigid%ozz(irgd)*x-rigid%oxx(irgd)*z
                         wzz=rigid%oxx(irgd)*y-rigid%oyy(irgd)*x

  ! new atomic velocities in lab frame

                         vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+rigid%vxx(irgd)
                         vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+rigid%vyy(irgd)
                         vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+rigid%vzz(irgd)
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

       If (rigid%total > 0) Then
          Do irgd=1,rigid%n_types
             rigid%vxx(irgd)=0.0_wp ; rigid%vyy(irgd)=0.0_wp ; rigid%vzz(irgd)=0.0_wp
             rigid%oxx(irgd)=0.0_wp ; rigid%oyy(irgd)=0.0_wp ; rigid%ozz(irgd)=0.0_wp
          End Do
       End If
    End If
  End Subroutine scale_temperature
End Module temperature
