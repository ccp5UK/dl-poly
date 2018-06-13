Module temperature
  Use kinds,           Only : wp, li
  Use comms,           Only : comms_type,gsum
  Use setup,           Only : nrite,boltz,mxatms,mxshl,zero_plus
  Use site, Only : site_type
  Use configuration,   Only : imcon,natms,nlast,nfree,lsite,  &
                              lsi,lsa,ltg,lfrzn,lfree,lstfre, &
                              weight,vxx,vyy,vzz,xxx,yyy,zzz
  Use rigid_bodies,    Only : rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz, &
                              rgdxxx,rgdyyy,rgdzzz,rgdx,rgdy,rgdz, &
                              rgdrix,rgdriy,rgdriz,rgdwgt,q0,q1,q2,q3, &
                              ntrgd,rgdmeg,lashp_rgd,lishp_rgd,lshmv_rgd, &
                              rgdfrz,listrgd,indrgd,getrotmat,rigid_bodies_quench
  Use constraints,     Only : constraints_type, constraints_quench
  Use pmf,             Only : pmf_quench, pmf_type
  Use core_shell,      Only : ntshl,listshl,legshl,lshmv_shl,lishp_shl, &
                              lashp_shl,core_shell_quench
  Use kinetics,        Only : l_vom,chvom,getcom,getvom,getkin,getknf,getknt,getknr
  Use numerics,        Only : invert,uni,local_index,box_mueller_saru3
  use shared_units,    Only : update_shared_units,update_shared_units_int
  Use errors_warnings, Only : error,warning,info
  Use thermostat,      Only : thermostat_type
  Use statistics, Only : stats_type

  Implicit None

  Private

  Public :: set_temperature, regauss_temperature, scale_temperature

Contains

  Subroutine set_temperature           &
             (levcfg,keyres,      &
             lmin,nstep,nstrun,nstmin, &
             keyshl,     &
             atmfre,atmfrz,            &
             megshl,     &
             megrgd,degtra,degrot,     &
             degfre,degshl,engrot,dof_site,stat,cons,pmf,thermo,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for setting the initial system temperature
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov january 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( In    ) :: lmin
    Integer,            Intent( In    ) :: nstep,nstrun,nstmin, &
                                           keyshl,       &
                                           atmfre,atmfrz,       &
                                           megshl,              &
                                           megrgd

    Integer,            Intent( InOut ) :: keyres,levcfg
    Integer(Kind=li),   Intent( InOut ) :: degtra,degrot

    Integer(Kind=li),   Intent(   Out ) :: degfre,degshl
    Real( Kind = wp ),  Intent(   Out ) :: engrot
    Real( Kind = wp ), Dimension(:), Intent( In    ) :: dof_site
    Type( stats_type ), Intent( InOut ) :: stat
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( constraints_type ), Intent( InOut ) :: cons
    Type( thermostat_type ), Intent( InOut ) :: thermo
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

    If (megrgd == 0) Then
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

    degshl=Int(3,li)*Int(megshl,li)

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

    If (megrgd > 0) Then
       engf=getknf(vxx,vyy,vzz,comm)/Real(Max(1_li,degfre),wp)
       engt=getknt(rgdvxx,rgdvyy,rgdvzz,comm)/Real(Max(1_li,degtra),wp)

       engr=getknr(rgdoxx,rgdoyy,rgdozz,comm)/Real(Max(1_li,degrot),wp)
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
       Allocate (qs(0:2,1:mxshl),tps(0:comm%mxnode-1), Stat=fail(2))
       If (Any(fail > 0)) Then
          Write(message,'(a)') 'set_temperature allocation failure'
          Call error(0,message)
       End If

  ! tpn(idnode) number of particles on this node (idnode)
  ! ntp - grand total of non-shell, non-frozen particles

       qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
       qs(0:2,1:ntshl) = 0 ! unqualified core-shell unit with a local shell

       j = 0
       Do i=1,natms
          If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. legshl(0,i) >= 0) Then
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
       If (keyshl == 1) Then ! just for the adiabatic shell model
          If (lshmv_shl) Then ! refresh the q array for shared core-shell units
             qn(natms+1:nlast) = 0
             Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,qn,comm)
             Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz,comm)
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
       tps(comm%idnode) = j
       Do i=0,comm%mxnode-1
          If (i /= comm%idnode) tps(i) = 0
       End Do
       Call gsum(comm,tps)
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
                   If (i <= natms) tpn(comm%idnode) = tpn(comm%idnode) - 1 ! Less free particles are hit
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

          If (lshmv_rgd) Then
            Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz,comm)
          End If

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
             Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,qn,comm)
             Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz,comm)
          End If

          If (tps(comm%idnode) > 0) Then
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
          Write(message,'(a)') 'set_temperature deallocation failure'
          Call error(0,message)
       End If

  ! remove centre of mass motion

       If (megrgd > 0) Then
          Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,comm)

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

       If (megrgd > 0) Then
         Call rigid_bodies_quench(comm)
       End If

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
          If (cons%megcon > 0) Then 
            Call constraints_quench(cons,stat,comm)
          End If
          If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,comm)
       End If

  ! quench core-shell units in adiabatic model

       If (megshl > 0 .and. keyshl == 1 .and. no_min_0) Then
          Do
             Call scale_temperature(thermo%sigma,degtra,degrot,degfre,comm)
             Call core_shell_quench(safe,thermo%temp,comm)
             If (cons%megcon > 0) Then
               Call constraints_quench(cons,stat,comm)
             End If
             If (pmf%megpmf > 0) Then
               Call pmf_quench(cons%max_iter_shake,cons%tolerance,stat,pmf,comm)
             End If
             If (megrgd > 0) Then
               Call rigid_bodies_quench(comm)
             End If
             If (safe) Exit
          End Do
       Else
          Call scale_temperature(thermo%sigma,degtra,degrot,degfre,comm)
       End If

    End If

    If (l_vom) Then

  ! remove centre of mass motion

       If (megrgd > 0) Then
          Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,comm)

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
    If (megrgd > 0) engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)

  End Subroutine set_temperature

  Subroutine regauss_temperature(megrgd,comm)

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

    Integer,            Intent( In    ) :: megrgd
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

    If (megrgd > 0) Then

  ! quench RBs

       Call rigid_bodies_quench(comm)

  ! remove centre of mass motion

       Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,comm)

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

  Subroutine scale_temperature(sigma,degtra,degrot,degfre,comm)

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
    Type ( comms_type ), Intent( InOut ) :: comm


    Integer           :: fail,i,j,irgd,jrgd,lrgd,rgdtyp,megrgd
    Real( Kind = wp ) :: engke,engkf,engkt,engrot,         &
                         amx,amy,amz,wxx,wyy,wzz,tmp,tmp1, &
                         com(1:3),vom(1:3),rot(1:9),rotinv(1:9),x,y,z

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    Character ( Len = 256 )  ::  message

  ! recover megrgd

    megrgd=rgdmeg

  ! remove centre of mass motion

    If (megrgd > 0) Then
       Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,comm)

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

       Call getcom(xxx,yyy,zzz,com,comm)

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
          Write(message,'(a)') 'scale_temperature deallocation failure'
          Call error(0,message)
       End If
    End If

  ! ensure equipartitioning - all degrees of freedom are equal
  ! calculate energy: free particles and RB translational and rotational

    engrot=0.0_wp
    If (megrgd > 0) Then
       engkf=getknf(vxx,vyy,vzz,comm)
       engkt=getknt(rgdvxx,rgdvyy,rgdvzz,comm)

       engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)

  ! temporary replacement for small engrot

       tmp1=Max(engrot,1.0e-6_wp)

  ! Scale rotational energy to translational energy
  ! according to their respective DoF

       If (degtra > Int(0,li)) Then ! engkt > 0 (degrot > 0 and tmp1(engrot) > 0)
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

       If (degfre-degtra-degrot > Int(0,li)) Then ! engkf > 0
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
          If (degtra > Int(0,li)) engkt=engkt*tmp**2
       End If

       engke=engkf+engkt
    Else
       engke=getkin(vxx,vyy,vzz,comm)
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
End Module temperature
