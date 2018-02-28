Subroutine dpd_thermostat(isw,l_str,rcut,nstep,tstep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine applying DPD thermostat in a Shardlow's VV manner
! using the verlet neighbour list
!
! isw=isw(VV) : by stages 0 for VV1 and 1 for VV2
! keydpd = 1 for first order splitting
! keydpd = 2 for second order splitting
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : idnode,mxnode,gsum
  Use setup_module,        Only : nrite,mxlist,mxatdm
  Use config_module,       Only : natms,nlast,lsi,lsa,ltg,ltype,lfree, &
                                  list,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use dpd_module
  Use rigid_bodies_module, Only : lshmv_rgd,lishp_rgd,lashp_rgd

  Implicit None

  Logical,           Intent( In    ) :: l_str
  Integer,           Intent( In    ) :: isw,nstep
  Real( Kind = wp ), Intent( In    ) :: rcut,tstep


  Integer           :: fail(1:2),nst_p,i,j,k,limit,idi,idj,ai,aj,key
  Real( Kind = wp ) :: tst_p,rstsq,hstep,fix,fiy,fiz, &
                       rrr,scrn,gauss,tmp,scl,        &
                       rgamma,dgamma,gamma,fx,fy,fz,  &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt
  Real( Kind = wp ), Dimension( : ), Allocatable :: fdpdx,fdpdy,fdpdz

  If (keydpd /= 1 .or. keydpd /= 2 .or. keydpd*isw == 1) Return

  fail=0
  Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat = fail(1))
  Allocate (fdpdx(1:mxatdm),fdpdy(1:mxatdm),fdpdz(1:mxatdm),         Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'dpd_thermostat allocation failure, node: ', idnode
     Call error(0)
  End If

! set tstep and nstep wrt to order of splitting

  If (keydpd == 1) Then
     nst_p = nstep
     tst_p = tstep
  Else
     If (isw == 0) Then
        nst_p = nstep
     Else ! If (isw == 1) Then
        nst_p = -nstep
     End If
     tst_p = 0.5_wp*tstep
  End If

! Set tstep derivatives

  hstep = 0.5_wp*tst_p
  rstsq = 1.0_wp/Sqrt(tst_p)

! initialise DPD virial and stress contributions

  If (isw == 0) Then
     virdpd = 0.0_wp
     strdpd = 0.0_wp
  End If

! FIRST PASS

! Initialise forces

  fdpdx = 0.0_wp
  fdpdy = 0.0_wp
  fdpdz = 0.0_wp

! Refresh halo velocities

  Call dpd_v_set_halo()

! outer loop over atoms

  Do i=1,natms

! Get list limit

     limit=Merge(list(0,i),0,weight(i) > 1.0e-6_wp)

! calculate interatomic distances

     Do k=1,limit
        j=list(k,i)

        xxt(k)=xxx(i)-xxx(j)
        yyt(k)=yyy(i)-yyy(j)
        zzt(k)=zzz(i)-zzz(j)
     End Do

! square of distances

     Do k=1,limit
        rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
     End Do

! initialise stress tensor accumulators

     strs1=0.0_wp
     strs2=0.0_wp
     strs3=0.0_wp
     strs5=0.0_wp
     strs6=0.0_wp
     strs9=0.0_wp

! global identity and atomic type of i

     idi=ltg(i)
     ai=ltype(i)

! load forces

     fix=fdpdx(i)
     fiy=fdpdy(i)
     fiz=fdpdz(i)

! start of primary loop for forces evaluation

     Do k=1,limit

! secondary atomic index

        j=list(k,i)

! interatomic distance

        rrr = rrt(k)

! validity of thermalisation

        If (rrr < rcut .and. weight(j) > 1.0e-6_wp) Then

! secondary atomic type and global index

           aj=ltype(j)
           idj=ltg(j)

! Get gaussian random number with zero mean

           Call box_mueller_saru2(idi,idj,nst_p,gauss,l_str)

! screening function

           scrn = (rcut-rrr)/(rrr*rcut)

! Get mixing type function

           If (ai > aj) Then
              key=ai*(ai-1)/2 + aj
           Else
              key=aj*(aj-1)/2 + ai
           End If

! Calculate force component

           rgamma =  sigdpd(key) * scrn      * gauss * rstsq

           tmp    =  gamdpd(key) * (scrn**2)
           dgamma = -tmp * ( xxt(k)*(vxx(i)-vxx(j)) + yyt(k)*(vyy(i)-vyy(j)) + zzt(k)*(vzz(i)-vzz(j)) )

           gamma=rgamma+dgamma

! calculate forces

           fx = gamma*xxt(k)
           fy = gamma*yyt(k)
           fz = gamma*zzt(k)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (j <= natms) Then

              fdpdx(j)=fdpdx(j)-fx
              fdpdy(j)=fdpdy(j)-fy
              fdpdz(j)=fdpdz(j)-fz

           End If

           If (j <= natms .or. idi < idj) Then

! add virial

              virdpd = virdpd - gamma*rrr*rrr

! add stress tensor

              strs1 = strs1 + xxt(k)*fx
              strs2 = strs2 + xxt(k)*fy
              strs3 = strs3 + xxt(k)*fz
              strs5 = strs5 + yyt(k)*fy
              strs6 = strs6 + yyt(k)*fz
              strs9 = strs9 + zzt(k)*fz

           End If

        End If

     End Do

! load back forces

     fdpdx(i)=fix
     fdpdy(i)=fiy
     fdpdz(i)=fiz

! complete stress tensor

     strdpd(1) = strdpd(1) + strs1
     strdpd(2) = strdpd(2) + strs2
     strdpd(3) = strdpd(3) + strs3
     strdpd(4) = strdpd(4) + strs2
     strdpd(5) = strdpd(5) + strs5
     strdpd(6) = strdpd(6) + strs6
     strdpd(7) = strdpd(7) + strs3
     strdpd(8) = strdpd(8) + strs6
     strdpd(9) = strdpd(9) + strs9

  End Do

! Update velocities or add to conservative forces

  Do i=1,natms
     If (lfree(i) == 0) Then
        If (weight(i) > 1.0e-6_wp) Then
           tmp=hstep/weight(i)
           vxx(i)=vxx(i)+tmp*fdpdx(i)
           vyy(i)=vyy(i)+tmp*fdpdy(i)
           vzz(i)=vzz(i)+tmp*fdpdz(i)
        End If
     Else ! a RB member
        fxx(i)=fxx(i)+fdpdx(i)
        fyy(i)=fyy(i)+fdpdy(i)
        fzz(i)=fzz(i)+fdpdz(i)
     End If
  End Do

! SECOND PASS

! Refresh halo velocities

  Call dpd_v_set_halo()

! Initialise forces

  fdpdx = 0.0_wp
  fdpdy = 0.0_wp
  fdpdz = 0.0_wp

! outer loop over atoms

  Do i=1,natms

! Get list limit

     limit=Merge(list(0,i),0,weight(i) > 1.0e-6_wp)

! calculate interatomic distances

     Do k=1,limit
        j=list(k,i)

        xxt(k)=xxx(i)-xxx(j)
        yyt(k)=yyy(i)-yyy(j)
        zzt(k)=zzz(i)-zzz(j)
     End Do

! square of distances

     Do k=1,limit
        rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
     End Do

! initialise stress tensor accumulators

     strs1=0.0_wp
     strs2=0.0_wp
     strs3=0.0_wp
     strs5=0.0_wp
     strs6=0.0_wp
     strs9=0.0_wp

! global identity and atomic type of i

     idi=ltg(i)
     ai=ltype(i)

! load forces

     fix=fdpdx(i)
     fiy=fdpdy(i)
     fiz=fdpdz(i)

! start of primary loop for forces evaluation

     Do k=1,limit

! secondary atomic index

        j=list(k,i)

! interatomic distance

        rrr = rrt(k)

! validity of thermalisation

        If (rrr < rcut .and. weight(j) > 1.0e-6_wp) Then

! secondary atomic type and global index

           aj=ltype(j)
           idj=ltg(j)

! Get gaussian random number with zero mean

           Call box_mueller_saru2(idi,idj,nst_p,gauss,l_str)

! screening function

           scrn = (rcut-rrr)/(rrr*rcut)

! Get mixing type function

           If (ai > aj) Then
              key=ai*(ai-1)/2 + aj
           Else
              key=aj*(aj-1)/2 + ai
           End If

! Calculate force component

           rgamma =  sigdpd(key) * scrn      * gauss * rstsq

           tmp    =  gamdpd(key) * (scrn**2)
           scl    =  tmp / (1.0_wp+tmp*tst_p)
           dgamma = -tmp * ( xxt(k)*(vxx(i)-vxx(j)) + yyt(k)*(vyy(i)-vyy(j)) + zzt(k)*(vzz(i)-vzz(j)) )

           gamma=rgamma + scl*(dgamma-rgamma)

! calculate forces

           fx = gamma*xxt(k)
           fy = gamma*yyt(k)
           fz = gamma*zzt(k)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (j <= natms) Then

              fdpdx(j)=fdpdx(j)-fx
              fdpdy(j)=fdpdy(j)-fy
              fdpdz(j)=fdpdz(j)-fz

           End If

           If (j <= natms .or. idi < idj) Then

! add virial

              virdpd = virdpd - gamma*rrr*rrr

! add stress tensor

              strs1 = strs1 + xxt(k)*fx
              strs2 = strs2 + xxt(k)*fy
              strs3 = strs3 + xxt(k)*fz
              strs5 = strs5 + yyt(k)*fy
              strs6 = strs6 + yyt(k)*fz
              strs9 = strs9 + zzt(k)*fz

           End If

        End If

     End Do

! load back forces

     fdpdx(i)=fix
     fdpdy(i)=fiy
     fdpdz(i)=fiz

! complete stress tensor

     strdpd(1) = strdpd(1) + strs1
     strdpd(2) = strdpd(2) + strs2
     strdpd(3) = strdpd(3) + strs3
     strdpd(4) = strdpd(4) + strs2
     strdpd(5) = strdpd(5) + strs5
     strdpd(6) = strdpd(6) + strs6
     strdpd(7) = strdpd(7) + strs3
     strdpd(8) = strdpd(8) + strs6
     strdpd(9) = strdpd(9) + strs9

  End Do

! Update velocities

  Do i=1,natms
     If (lfree(i) == 0) Then
        If (weight(i) > 1.0e-6_wp) Then
           tmp=hstep/weight(i)
           vxx(i)=vxx(i)+tmp*fdpdx(i)
           vyy(i)=vyy(i)+tmp*fdpdy(i)
           vzz(i)=vzz(i)+tmp*fdpdz(i)
        End If
     Else ! a RB member
        fxx(i)=fxx(i)+fdpdx(i)
        fyy(i)=fyy(i)+fdpdy(i)
        fzz(i)=fzz(i)+fdpdz(i)
     End If
  End Do

! Update forces on RBs

  If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz)

! globalise virdpd

  If (mxnode > 1) Call gsum(virdpd)

  Deallocate (xxt,yyt,zzt,rrt,   Stat = fail(1))
  Deallocate (fdpdx,fdpdy,fdpdz, Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'dpd_thermostat deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine dpd_thermostat
