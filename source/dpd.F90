Module dpd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global DPD variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  
  Use comms,        Only : comms_type,gsum,gcheck,gmax,DpdVExp_tag,wp_mpi, &
                           gsend,gwait,girecv
  Use setup,        Only : nrite,mxlist,mxatdm,mxatms,mxbfxp,mxvdw
  Use configuration,       Only : natms,nlast,lsi,lsa,ltg,ltype,lfree, &
                                  list,weight,xxx,yyy,zzz,vxx,vyy,vzz, &
                                  fxx,fyy,fzz, ixyz
  Use rigid_bodies, Only : lshmv_rgd,lishp_rgd,lashp_rgd
  Use domains

  Use shared_units,    Only : update_shared_units
  Use errors_warnings, Only : error, warning
  Use numerics,        Only : box_mueller_saru2

  Implicit None

  Integer,           Save :: keydpd = 0 ! no DPD

  Real( Kind = wp ), Save :: virdpd      = 0.0_wp , &
                             strdpd(1:9) = 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: gamdpd(:),sigdpd(:)

  Public :: allocate_dpd_arrays

Contains

  Subroutine allocate_dpd_arrays()

    Integer :: fail

    If (keydpd == 0) Return

    fail = 0

    Allocate (gamdpd(0:mxvdw),sigdpd(1:mxvdw), Stat = fail)

    If (fail > 0) Call error(1081)

    gamdpd = 0.0_wp ; sigdpd = 0.0_wp

  End Subroutine allocate_dpd_arrays

  Subroutine dpd_thermostat(isw,l_str,rcut,nstep,tstep,comm)

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

    Logical,           Intent( In    ) :: l_str
    Integer,           Intent( In    ) :: isw,nstep
    Real( Kind = wp ), Intent( In    ) :: rcut,tstep
    Type( comms_type ), Intent( InOut ) :: comm


    Integer           :: fail(1:2),nst_p,i,j,k,limit,idi,idj,ai,aj,key
    Real( Kind = wp ) :: tst_p,rstsq,hstep,fix,fiy,fiz, &
      rrr,scrn,gauss,tmp,scl,        &
      rgamma,dgamma,gamma,fx,fy,fz,  &
      strs1,strs2,strs3,strs5,strs6,strs9

    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt
    Real( Kind = wp ), Dimension( : ), Allocatable :: fdpdx,fdpdy,fdpdz
    Character ( len = 256 ) :: message


    If (keydpd /= 1 .or. keydpd /= 2 .or. keydpd*isw == 1) Return

    fail=0
    Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat = fail(1))
    Allocate (fdpdx(1:mxatdm),fdpdy(1:mxatdm),fdpdz(1:mxatdm),         Stat = fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'dpd_thermostat allocation failure'
      Call error(0,message)
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

    Call dpd_v_set_halo(comm)

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

    Call dpd_v_set_halo(comm)

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

    If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxx,fyy,fzz,comm)

    ! globalise virdpd

    Call gsum(comm,virdpd)

    Deallocate (xxt,yyt,zzt,rrt,   Stat = fail(1))
    Deallocate (fdpdx,fdpdy,fdpdz, Stat = fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'dpd_thermostat deallocation failure'
      Call error(0,message)
    End If

  End Subroutine dpd_thermostat

  Subroutine dpd_v_export(mdir,mlast,ixyz0,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export metal density data in domain boundary
    ! regions for halo formation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Integer,           Intent( In    ) :: mdir
    Integer,           Intent( InOut ) :: mlast,ixyz0(1:mxatms)
    Type( comms_type ), Intent( InOut ) :: comm

    Logical           :: safe
    Integer           :: fail,iadd,limit,iblock,          &
      i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
      jdnode,kdnode,imove,jmove,itmp

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    Character ( Len = 256 )   ::  message
    ! Number of transported quantities per particle

    iadd=4

    fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'dpd_v_export allocation failure'
      Call error(0,message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock=limit/Merge(2,1,comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! kxyz - corrected halo reduction factor particles haloing both +&- sides
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0 ; ky = 0 ; kz = 0
    If      (mdir == -1) Then ! Direction -x
      kx  = 1
      jxyz= 1
      kxyz= 3

      jdnode = map(1)
      kdnode = map(2)
    Else If (mdir ==  1) Then ! Direction +x
      kx  = 1
      jxyz= 2
      kxyz= 3

      jdnode = map(2)
      kdnode = map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky  = 1
      jxyz= 10
      kxyz= 30

      jdnode = map(3)
      kdnode = map(4)
    Else If (mdir ==  2) Then ! Direction +y
      ky  = 1
      jxyz= 20
      kxyz= 30

      jdnode = map(4)
      kdnode = map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz  = 1
      jxyz= 100
      kxyz= 300

      jdnode = map(5)
      kdnode = map(6)
    Else If (mdir ==  3) Then ! Direction +z
      kz  = 1
      jxyz= 200
      kxyz= 300

      jdnode = map(6)
      kdnode = map(5)
    Else
      Call error(152)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove=0
    jmove=0

    ! Initialise array overflow flags

    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE


    Do i=1,mlast

      ! If the particle is within the remaining 'inverted halo' of this domain

      If (ixyz0(i) > 0) Then

        ! Get the necessary halo indices

        ix=Mod(ixyz0(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz0(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz0(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

        ! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

        ! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

          ! If safe to proceed

          If ((imove+iadd) <= iblock) Then

            ! pack particle velocity and halo indexing

            buffer(imove+1)=vxx(i)
            buffer(imove+2)=vyy(i)
            buffer(imove+3)=vzz(i)

            ! Use the corrected halo reduction factor when the particle is halo to both +&- sides

            buffer(imove+4)=Real(ixyz0(i)-Merge(jxyz,kxyz,j == jxyz),wp)

          Else

            safe=.false.

          End If
          imove=imove+iadd

        End If

      End If

    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=Merge(2,1,comm%mxnode > 1)*imove
      Call gmax(comm,itmp)
      Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
      Call error(154)
    End If

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm,jmove,kdnode,DpdVExp_tag)
      Call gsend(comm,imove,jdnode,DpdVExp_tag)
      Call gwait(comm)
    Else
      jmove=imove
    End If

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe=((mlast+jmove/iadd) <= mxatms)
    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=mlast+jmove/iadd
      Call gmax(comm,itmp)
      Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
      Call error(156)
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,DpdVExp_tag)
      End If
      If (imove > 0) Then
        Call gsend(comm,buffer(1:imove),jdnode,DpdVExp_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data

    j=Merge(iblock,0,comm%mxnode > 1)
    Do i=1,jmove/iadd
      mlast=mlast+1

      ! unpack particle velocity and remaining halo indexing

      vxx(mlast)=buffer(j+1)
      vyy(mlast)=buffer(j+2)
      vzz(mlast)=buffer(j+3)
      ixyz0(mlast)=Nint(buffer(j+4))

      j=j+iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'dpd_v_export deallocation failure'
      Call error(0,message)
    End If

  End Subroutine dpd_v_export

  Subroutine dpd_v_set_halo(comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange exchange of velocity data between
    ! neighbouring domains/nodes
    !
    ! copyright - daresbury laboratory
    ! amended   - i.t.todorov november 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( comms_type ), Intent( InOut ) :: comm
    Logical :: safe
    Integer :: fail,mlast

    Integer, Allocatable :: ixyz0(:)
    Character ( Len = 256 )  ::  message

    fail = 0
    Allocate (ixyz0(1:mxatms), Stat = fail)
    If (fail > 0) Then
      Write(message,'(a)') 'dpd_v_set_halo allocation failure'
      Call error(0,message)
    End If
    ixyz0(1:nlast) = ixyz(1:nlast)

    ! No halo, start with domain only particles

    mlast=natms

    ! exchange atom data in -/+ x directions

    Call dpd_v_export(-1,mlast,ixyz0,comm)
    Call dpd_v_export( 1,mlast,ixyz0,comm)

    ! exchange atom data in -/+ y directions

    Call dpd_v_export(-2,mlast,ixyz0,comm)
    Call dpd_v_export( 2,mlast,ixyz0,comm)

    ! exchange atom data in -/+ z directions

    Call dpd_v_export(-3,mlast,ixyz0,comm)
    Call dpd_v_export( 3,mlast,ixyz0,comm)

    ! check atom totals after data transfer

    safe=(mlast == nlast)
    Call gcheck(comm,safe)
    If (.not.safe) Call error(96)

    Deallocate (ixyz0, Stat = fail)
    If (fail > 0) Then
      Write(message,'(a)') 'dpd_v_set_halo deallocation failure'
      Call error(0,message)
    End If

  End Subroutine dpd_v_set_halo
End Module dpd
