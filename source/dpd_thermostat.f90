Subroutine dpd_thermostat(l_str,imcon,rcut,nstep,tstep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine applying DPD thermostat in a Shardlow's VV manner
! using the verlet neighbour list
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode,mxnode,gsum
  Use setup_module,        Only : nrite,mxlist,mxatdm
  Use config_module,       Only : cell,natms,nlast,lsi,lsa,ltg,ltype,lfree, &
                                  list,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use dpd_module
  Use rigid_bodies_module, Only : lshmv_rgd,lishp_rgd,lashp_rgd

  Implicit None

  Logical,           Intent( In    ) :: l_str
  Integer,           Intent( In    ) :: imcon,nstep
  Real( Kind = wp ), Intent( In    ) :: rcut,tstep


  Integer           :: fail(1:2),i,j,k,limit,idi,idj,ai,aj,key
  Real( Kind = wp ) :: rcsq,trsq,hstep,fix,fiy,fiz,  &
                       rsq,rrr,scrn,gauss,tmp,scl,   &
                       rgamma,dgamma,gamma,fx,fy,fz, &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Real( Kind = wp ), Dimension( : ), Allocatable :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ), Dimension( : ), Allocatable :: fdpdx,fdpdy,fdpdz

  fail=0
  Allocate (xdf(1:mxlist),ydf(1:mxlist),zdf(1:mxlist),rsqdf(1:mxlist), Stat = fail(1))
  Allocate (fdpdx(1:mxatdm),fdpdy(1:mxatdm),fdpdz(1:mxatdm),           Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'dpd_thermostat allocation failure, node: ', idnode
     Call error(0)
  End If

! set cutoff condition snd tstep factors

  rcsq = rcut**2
  hstep= 0.5_wp*tstep
  trsq = 1.0_wp/Sqrt(tstep)

! initialise DPD virial and stress contributions

  virdpd = 0.0_wp
  strdpd = 0.0_wp

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

        xdf(k)=xxx(i)-xxx(j)
        ydf(k)=yyy(i)-yyy(j)
        zdf(k)=zzz(i)-zzz(j)
     End Do

! periodic boundary conditions

     Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distances

     Do k=1,limit
        rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
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

        rsq = rsqdf(k)

! validity of thermalisation

        If (rsq < rcsq .and. weight(j) > 1.0e-6_wp) Then

! secondary atomic type and global index

           aj=ltype(j)
           idj=ltg(j)

! Get gaussian random number with zero mean

           Call box_mueller_saru2(i,j,nstep,gauss,l_str)

! Get separation distance and screening function

           rrr = Sqrt(rsq)
           scrn = (rcut-rrr)/(rrr*rcut)

! Get mixing type function

           If (ai > aj) Then
              key=ai*(ai-1)/2 + aj
           Else
              key=aj*(aj-1)/2 + ai
           End If

! Calculate force component

           rgamma =  sigdpd(key) * scrn      * gauss * trsq

           tmp    =  gamdpd(key) * (scrn**2)
           dgamma = -tmp * ( xdf(k)*(vxx(i)-vxx(j)) + ydf(k)*(vyy(i)-vyy(j)) + zdf(k)*(vzz(i)-vzz(j)) )

           gamma=rgamma+dgamma

! calculate forces

           fx = gamma*xdf(k)
           fy = gamma*ydf(k)
           fz = gamma*zdf(k)

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

              virdpd = virdpd - gamma*rsq

! add stress tensor

              strs1 = strs1 + xdf(k)*fx
              strs2 = strs2 + xdf(k)*fy
              strs3 = strs3 + xdf(k)*fz
              strs5 = strs5 + ydf(k)*fy
              strs6 = strs6 + ydf(k)*fz
              strs9 = strs9 + zdf(k)*fz

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

        xdf(k)=xxx(i)-xxx(j)
        ydf(k)=yyy(i)-yyy(j)
        zdf(k)=zzz(i)-zzz(j)
     End Do

! periodic boundary conditions

     Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distances

     Do k=1,limit
        rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
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

        rsq = rsqdf(k)

! validity of thermalisation

        If (rsq < rcsq .and. weight(j) > 1.0e-6_wp) Then

! secondary atomic type and global index

           aj=ltype(j)
           idj=ltg(j)

! Get gaussian random number with zero mean

           Call box_mueller_saru2(i,j,nstep,gauss,l_str)

! Get separation distance and screening function

           rrr = Sqrt(rsq)
           scrn = (rcut-rrr)/(rrr*rcut)

! Get mixing type function

           If (ai > aj) Then
              key=ai*(ai-1)/2 + aj
           Else
              key=aj*(aj-1)/2 + ai
           End If

! Calculate force component

           rgamma =  sigdpd(key) * scrn      * gauss * trsq

           tmp    =  gamdpd(key) * (scrn**2)
           scl    =  tmp / (1.0_wp+tmp*tstep)
           dgamma = -tmp * ( xdf(k)*(vxx(i)-vxx(j)) + ydf(k)*(vyy(i)-vyy(j)) + zdf(k)*(vzz(i)-vzz(j)) )

           gamma=rgamma + scl*(dgamma-rgamma)

! calculate forces

           fx = gamma*xdf(k)
           fy = gamma*ydf(k)
           fz = gamma*zdf(k)

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

              virdpd = virdpd - gamma*rsq

! add stress tensor

              strs1 = strs1 + xdf(k)*fx
              strs2 = strs2 + xdf(k)*fy
              strs3 = strs3 + xdf(k)*fz
              strs5 = strs5 + ydf(k)*fy
              strs6 = strs6 + ydf(k)*fz
              strs9 = strs9 + zdf(k)*fz

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

  Deallocate (xdf,ydf,zdf,rsqdf, Stat = fail(1))
  Deallocate (fdpdx,fdpdy,fdpdz, Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'dpd_thermostat deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine dpd_thermostat
