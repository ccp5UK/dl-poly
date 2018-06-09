Module rigid_bodies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining rigid bodies' (RBs) variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp,wi,li
  Use comms,           Only : comms_type,gsum,gmin,gmax,gsync,gcheck
  Use setup,           Only : mxtmls,mxtrgd,mxrgd,mxlrgd,mxfrgd,mxlshp,mxproc,mxatdm, &
                              mxatms,zero_plus
  Use site
  Use configuration,   Only : imcon,cell,natms,nlast,lsi,lsa,xxx,yyy,zzz,vxx,vyy,vzz, &
                              ltg,lsite,lfrzn,fxx,fyy,fzz,nfree,lstfre,getcom
  Use neighbours,      Only : neighbours_type
  Use statistics,      Only : stats_type
  Use numerics,        Only : images, jacobi, invert
  Use shared_units,    Only : update_shared_units
  Use numerics,        Only : images,local_index,pbcshift
  Use errors_warnings, Only : info, error, warning
  Use thermostat,      Only : thermostat_type, &
                              ENS_NPT_BERENDSEN, ENS_NPT_BERENDSEN_ANISO, &
                              ENS_NPT_LANGEVIN, ENS_NPT_LANGEVIN_ANISO, &
                              ENS_NPT_NOSE_HOOVER, ENS_NPT_NOSE_HOOVER_ANISO, &
                              ENS_NPT_MTK, ENS_NPT_MTK_ANISO

  Implicit None

  Logical,                        Save :: lshmv_rgd = .false.

  Integer,                        Save :: ntrgd  = 0 , &
                                          ntrgd1 = 0 , &
                                          m_rgd  = 0

! Globalised megrgd & rcut

  Integer,                        Save :: rgdmeg = 0
  Real( Kind = wp ),              Save :: rgdrct = 0.0_wp

  Integer,           Allocatable, Save :: numrgd(:)
  Integer,           Allocatable, Save :: lstrgd(:,:),listrgd(:,:),legrgd(:,:)
  Integer,           Allocatable, Save :: rgdfrz(:,:),rgdind(:,:)
  Integer,           Allocatable, Save :: lishp_rgd(:),lashp_rgd(:)
  Integer,           Allocatable, Save :: indrgd(:,:)

  Real( Kind = wp ), Allocatable, Save :: rgdwgt(:,:),rgdwg1(:,:)
  Real( Kind = wp ), Allocatable, Save :: rgdx(:,:),rgdy(:,:),rgdz(:,:)
  Real( Kind = wp ), Allocatable, Save :: rgdrix(:,:),rgdriy(:,:),rgdriz(:,:)
  Real( Kind = wp ), Allocatable, Save :: rgdaxs(:,:)
  Real( Kind = wp ), Allocatable, Save :: q0(:),q1(:),q2(:),q3(:)
  Real( Kind = wp ), Allocatable, Save :: rgdxxx(:),rgdyyy(:),rgdzzz(:)
  Real( Kind = wp ), Allocatable, Save :: rgdvxx(:),rgdvyy(:),rgdvzz(:)
  Real( Kind = wp ), Allocatable, Save :: rgdoxx(:),rgdoyy(:),rgdozz(:)

  Public :: allocate_rigid_bodies_arrays , deallocate_rigid_bodies_arrays, &
            rigid_bodies_stress, rigid_bodies_stre_s, rigid_bodies_str_ss, &
            rigid_bodies_str__s, getrotmat, q_setup


Contains

  Subroutine allocate_rigid_bodies_arrays()

    Integer, Dimension( 1:15 ) :: fail

    fail = 0

    Allocate (numrgd(1:mxtmls),                                                        Stat = fail( 1))
    Allocate (lstrgd(0:mxlrgd,1:mxtrgd),                                               Stat = fail( 2))
    Allocate (listrgd(-1:mxlrgd,1:mxrgd),                                              Stat = fail( 3))
    Allocate (legrgd(0:mxfrgd,1:mxatdm),                                               Stat = fail( 4))
    Allocate (lishp_rgd(1:mxlshp),lashp_rgd(1:mxproc),                                 Stat = fail( 5))
    Allocate (rgdfrz(0:mxlrgd,1:mxtrgd),rgdind(0:mxlrgd,1:mxtrgd),                     Stat = fail( 6))
    Allocate (rgdwgt(0:mxlrgd,1:mxtrgd),rgdwg1(0:mxlrgd,1:mxtrgd),                     Stat = fail( 7))
    Allocate (indrgd(0:mxlrgd,1:mxrgd),                                                Stat = fail( 8))
    Allocate (rgdx(1:mxlrgd,1:mxtrgd),rgdy(1:mxlrgd,1:mxtrgd),rgdz(1:mxlrgd,1:mxtrgd), Stat = fail( 9))
    Allocate (rgdrix(1:2,1:mxtrgd),rgdriy(1:2,1:mxtrgd),rgdriz(1:2,1:mxtrgd),          Stat = fail(10))
    Allocate (rgdaxs(1:9,1:mxtrgd),                                                    Stat = fail(11))
    Allocate (q0(1:mxrgd),q1(1:mxrgd),q2(1:mxrgd),q3(1:mxrgd),                         Stat = fail(12))
    Allocate (rgdxxx(1:mxrgd),rgdyyy(1:mxrgd),rgdzzz(1:mxrgd),                         Stat = fail(13))
    Allocate (rgdvxx(1:mxrgd),rgdvyy(1:mxrgd),rgdvzz(1:mxrgd),                         Stat = fail(14))
    Allocate (rgdoxx(1:mxrgd),rgdoyy(1:mxrgd),rgdozz(1:mxrgd),                         Stat = fail(15))

    If (Any(fail > 0)) Call error(1042)

    numrgd  = 0
    lstrgd  = 0
    listrgd = 0
    legrgd  = 0

    lishp_rgd = 0 ; lashp_rgd = 0

    rgdfrz = 0 ; rgdind = 0 ; indrgd = 0

    rgdwgt = 0.0_wp ; rgdwg1 = 0.0_wp
    rgdx   = 0.0_wp ; rgdy   = 0.0_wp ; rgdz   = 0.0_wp
    rgdrix = 0.0_wp ; rgdriy = 0.0_wp ; rgdriz = 0.0_wp
    rgdaxs = 0.0_wp

    q0 = 0.0_wp ; q1 = 0.0_wp ; q2 = 0.0_wp ; q3 = 0.0_wp

    rgdxxx = 0.0_wp ; rgdyyy = 0.0_wp ; rgdzzz = 0.0_wp
    rgdvxx = 0.0_wp ; rgdvyy = 0.0_wp ; rgdvzz = 0.0_wp
    rgdoxx = 0.0_wp ; rgdoyy = 0.0_wp ; rgdozz = 0.0_wp

  End Subroutine allocate_rigid_bodies_arrays

  Subroutine deallocate_rigid_bodies_arrays()

    Integer :: fail

    fail = 0

    Deallocate (numrgd,lstrgd, Stat = fail)

    If (fail > 0) Call error(1043)

  End Subroutine deallocate_rigid_bodies_arrays

  Subroutine rigid_bodies_coms(xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for constructing RBs coms
  !
  ! Note: it assumes that all RBs' members are present and fresh
  ! (even those in the halo of the domain)
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent( In    ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
    Real( Kind = wp ),  Intent(   Out ) :: rgdxxx(1:mxrgd),rgdyyy(1:mxrgd),rgdzzz(1:mxrgd)
    Type( comms_type ), Intent( In    ) :: comm

    Integer           :: fail,irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ) :: tmp

    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Character ( Len = 256 ) :: message

    fail = 0
    Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_coms allocation failure'
       Call error(0,message)
    End If

  ! Loop over all local RB units and get in the local scope of the unit

    krgd=0
    Do irgd=1,ntrgd
       lrgd=listrgd(-1,irgd)
       Do jrgd=1,lrgd
          krgd=krgd+1

          gxx(krgd) = xxx(indrgd(jrgd,irgd)) - xxx(indrgd(1,irgd))
          gyy(krgd) = yyy(indrgd(jrgd,irgd)) - yyy(indrgd(1,irgd))
          gzz(krgd) = zzz(indrgd(jrgd,irgd)) - zzz(indrgd(1,irgd))
       End Do
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Get the COM vector

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       rgdxxx(irgd)=0.0_wp
       rgdyyy(irgd)=0.0_wp
       rgdzzz(irgd)=0.0_wp

       lrgd=listrgd(-1,irgd)
       Do jrgd=1,lrgd
          krgd=krgd+1

          rgdxxx(irgd) = rgdxxx(irgd) + rgdwg1(jrgd,rgdtyp)*gxx(krgd)
          rgdyyy(irgd) = rgdyyy(irgd) + rgdwg1(jrgd,rgdtyp)*gyy(krgd)
          rgdzzz(irgd) = rgdzzz(irgd) + rgdwg1(jrgd,rgdtyp)*gzz(krgd)
       End Do
    End Do

    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       tmp=1.0_wp/rgdwg1(0,rgdtyp)

       rgdxxx(irgd) = rgdxxx(irgd)*tmp + xxx(indrgd(1,irgd))
       rgdyyy(irgd) = rgdyyy(irgd)*tmp + yyy(indrgd(1,irgd))
       rgdzzz(irgd) = rgdzzz(irgd)*tmp + zzz(indrgd(1,irgd))
    End Do

    Deallocate (gxx,gyy,gzz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_coms deallocation failure'
       Call error(0,message)
    End If

  End Subroutine rigid_bodies_coms

  Subroutine rigid_bodies_move(stride,oxx,oyy,ozz,txx,tyy,tzz,uxx,uyy,uzz,dist_tol)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine for updating positions of atoms in RBs
  ! during a conjugate gradient minimisation
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith may 2006
  ! adapted   - i.t.todorov october 2012
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent( In    ) :: stride,                                    &
                                          oxx(1:mxatms),oyy(1:mxatms),ozz(1:mxatms), &
                                          txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms), &
                                          uxx(1:mxatms),uyy(1:mxatms),uzz(1:mxatms)
    Real( Kind = wp ), Intent( InOut ) :: dist_tol

    Integer           :: i,irgd,jrgd,lrgd,rgdtyp
    Real( Kind = wp ) :: ttt,uuu,the,dtol,coz,zin,x,y,z

    dtol=0.0_wp
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

  ! For all good RBs

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                i=indrgd(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then

  ! translational motion

                   If (rgdfrz(0,rgdtyp) == 0) Then
                      x=stride*oxx(i)
                      y=stride*oyy(i)
                      z=stride*ozz(i)
                   Else
                      x=0.0_wp
                      y=0.0_wp
                      z=0.0_wp
                   End If

  ! calculate this RB's particle distance to the axis of rotation

                   uuu=Sqrt(uxx(i)**2+uyy(i)**2+uzz(i)**2)

  ! add rotational motion for big enough distances

                   If (uuu > 1.0e-10_wp) Then

  ! force magnitude

                      ttt=Sqrt(txx(i)**2+tyy(i)**2+tzz(i)**2)

  ! angle of rotation

                      the=(ttt/uuu)*stride

                      coz=Cos(the)-1.0_wp
                      zin=(Sin(the)/the)*stride

                      x=x+coz*uxx(i)+zin*txx(i)
                      y=y+coz*uyy(i)+zin*tyy(i)
                      z=z+coz*uzz(i)+zin*tzz(i)
                   End If

  ! Add motion

                   xxx(i)=xxx(i)+x
                   yyy(i)=yyy(i)+y
                   zzz(i)=zzz(i)+z

                   dtol=Max(dtol,x**2+y**2+z**2)

                End If
             End If
          End Do
       End If
    End Do

    dist_tol=Max(dist_tol,Sqrt(dtol))
  End Subroutine rigid_bodies_move

  Subroutine rigid_bodies_quench(comm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to convert atomic velocities to RB COM and
  ! angular velocity
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: fail,i,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ) :: weight,rot(1:9),wxx,wyy,wzz,x(1:1),y(1:1),z(1:1),tmp

    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Character ( Len = 256 )        :: message


    fail = 0
    Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_quench allocation failure'
       Call error(0,message)
    End If

  ! Halo velocity field across onto neighbouring domains

    If (lshmv_rgd) Then
       Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz,comm)
    End If

  ! translate atomic velocities to COM velocity & angular velocity
  ! frozen velocities and massless sites weights are assumed zero!!!

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       rgdvxx(irgd)=0.0_wp ; rgdvyy(irgd)=0.0_wp ; rgdvzz(irgd)=0.0_wp

  ! For all good RBs

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters
          Do jrgd=1,lrgd
             krgd=krgd+1

  ! local index and mass of particle/site

             i=indrgd(jrgd,irgd)
             weight=rgdwgt(jrgd,rgdtyp)

  ! COM distances

             gxx(krgd)=xxx(i)-rgdxxx(irgd)
             gyy(krgd)=yyy(i)-rgdyyy(irgd)
             gzz(krgd)=zzz(i)-rgdzzz(irgd)

  ! If the RB has a frozen particle then no net COM momentum

             If (rgdfrz(0,rgdtyp) == 0) Then
                rgdvxx(irgd)=rgdvxx(irgd)+weight*vxx(i)
                rgdvyy(irgd)=rgdvyy(irgd)+weight*vyy(i)
                rgdvzz(irgd)=rgdvzz(irgd)+weight*vzz(i)
             End If
          End Do

  ! COM velocity
  ! If the RB has a frozen particle then no net COM momentum

          If (rgdfrz(0,rgdtyp) == 0) Then
             tmp=1.0_wp/rgdwgt(0,rgdtyp)
             rgdvxx(irgd)=rgdvxx(irgd)*tmp
             rgdvyy(irgd)=rgdvyy(irgd)*tmp
             rgdvzz(irgd)=rgdvzz(irgd)*tmp
          End If
       End If
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Get RBs' angular momenta

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       rgdoxx(irgd)=0.0_wp ; rgdoyy(irgd)=0.0_wp; rgdozz(irgd)=0.0_wp

  ! For all good RBs

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

  ! new rotational matrix

          Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

  ! angular momentum accumulators

          wxx=0.0_wp
          wyy=0.0_wp
          wzz=0.0_wp

          Do jrgd=1,lrgd
             krgd=krgd+1

  ! local index and mass of particle/site - assumption must hold here

             i=indrgd(jrgd,irgd)
             weight=rgdwgt(jrgd,rgdtyp)

             wxx=wxx+weight*(gyy(krgd)*vzz(i)-gzz(krgd)*vyy(i))
             wyy=wyy+weight*(gzz(krgd)*vxx(i)-gxx(krgd)*vzz(i))
             wzz=wzz+weight*(gxx(krgd)*vyy(i)-gyy(krgd)*vxx(i))
          End Do

  ! If the RB has 2+ frozen particles (ill=1) the net angular momentum
  ! must align along the axis of rotation keeping its magnitude

          If (rgdfrz(0,rgdtyp) > 1) Then
             i1=indrgd(rgdind(1,rgdtyp),irgd)
             i2=indrgd(rgdind(2,rgdtyp),irgd)

             x(1)=xxx(i1)-xxx(i2)
             y(1)=yyy(i1)-yyy(i2)
             z(1)=zzz(i1)-zzz(i2)

             Call images(imcon,cell,1,x,y,z)

             tmp=Sign( Sqrt((wxx**2+wyy**2+wzz**2)/(x(1)**2+y(1)**2+z(1)**2)) , &
                       Nearest((x(1)*wxx+y(1)*wyy+z(1)*wzz),-1.0_wp) )

             wxx=x(1)*tmp
             wyy=y(1)*tmp
             wzz=z(1)*tmp
          End If

  ! angular velocity in body fixed frame

          rgdoxx(irgd)=(rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rgdrix(2,rgdtyp)
          rgdoyy(irgd)=(rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rgdriy(2,rgdtyp)
          rgdozz(irgd)=(rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rgdriz(2,rgdtyp)

          Do jrgd=1,lrgd
             If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                i=indrgd(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then

  ! site velocity in body frame

                   x(1)=rgdx(jrgd,rgdtyp)
                   y(1)=rgdy(jrgd,rgdtyp)
                   z(1)=rgdz(jrgd,rgdtyp)

                   wxx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                   wyy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                   wzz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

  ! new atomic velocities in lab frame

                   vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+rgdvxx(irgd)
                   vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+rgdvyy(irgd)
                   vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+rgdvzz(irgd)

                End If
             End If
          End Do
       End If
    End Do

    Deallocate (gxx,gyy,gzz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_quench deallocation failure'
       Call error(0,message)
    End If
  End Subroutine rigid_bodies_quench

  Subroutine rigid_bodies_q_ench(qr,comm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to convert atomic velocities to RB COM and
  ! angular velocity
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: qr(1:mxrgd)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: fail,i,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ) :: weight,rot(1:9),wxx,wyy,wzz,x(1:1),y(1:1),z(1:1),tmp

    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Character ( Len = 256 )        :: message


    fail = 0
    Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_q_ench allocation failure'
       Call error(0,message)
    End If

  ! Halo velocity field across onto neighbouring domains

    If (lshmv_rgd) Then
     Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,vxx,vyy,vzz,comm)
    End If

  ! translate atomic velocities to COM velocity & angular velocity
  ! frozen velocities and massless sites weights are assumed zero!!!

    krgd=0
    Do irgd=1,ntrgd
       If (qr(irgd) == 1) Then
          rgdtyp=listrgd(0,irgd)

          rgdvxx(irgd)=0.0_wp ; rgdvyy(irgd)=0.0_wp ; rgdvzz(irgd)=0.0_wp

  ! For all good RBs

          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters
             Do jrgd=1,lrgd
                krgd=krgd+1

  ! local index and mass of particle/site

                i=indrgd(jrgd,irgd)
                weight=rgdwgt(jrgd,rgdtyp)

  ! COM distances

                gxx(krgd)=xxx(i)-rgdxxx(irgd)
                gyy(krgd)=yyy(i)-rgdyyy(irgd)
                gzz(krgd)=zzz(i)-rgdzzz(irgd)

  ! If the RB has a frozen particle then no net COM momentum

                If (rgdfrz(0,rgdtyp) == 0) Then
                   rgdvxx(irgd)=rgdvxx(irgd)+weight*vxx(i)
                   rgdvyy(irgd)=rgdvyy(irgd)+weight*vyy(i)
                   rgdvzz(irgd)=rgdvzz(irgd)+weight*vzz(i)
                End If
             End Do

  ! COM velocity
  ! If the RB has a frozen particle then no net COM momentum

             If (rgdfrz(0,rgdtyp) == 0) Then
                tmp=1.0_wp/rgdwgt(0,rgdtyp)
                rgdvxx(irgd)=rgdvxx(irgd)*tmp
                rgdvyy(irgd)=rgdvyy(irgd)*tmp
                rgdvzz(irgd)=rgdvzz(irgd)*tmp
             End If
          End If
       End If
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Get RBs' angular momenta

    krgd=0
    Do irgd=1,ntrgd
       If (qr(irgd) == 1) Then
          rgdtyp=listrgd(0,irgd)

          rgdoxx(irgd)=0.0_wp ; rgdoyy(irgd)=0.0_wp; rgdozz(irgd)=0.0_wp

  ! For all good RBs

          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

  ! new rotational matrix

             Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

  ! angular momentum accumulators

             wxx=0.0_wp
             wyy=0.0_wp
             wzz=0.0_wp

             Do jrgd=1,lrgd
                krgd=krgd+1

  ! local index and mass of particle/site - assumption must hold here

                i=indrgd(jrgd,irgd)
                weight=rgdwgt(jrgd,rgdtyp)

                wxx=wxx+weight*(gyy(krgd)*vzz(i)-gzz(krgd)*vyy(i))
                wyy=wyy+weight*(gzz(krgd)*vxx(i)-gxx(krgd)*vzz(i))
                wzz=wzz+weight*(gxx(krgd)*vyy(i)-gyy(krgd)*vxx(i))
             End Do

  ! If the RB has 2+ frozen particles (ill=1) the net angular momentum
  ! must align along the axis of rotation keeping its magnitude

             If (rgdfrz(0,rgdtyp) > 1) Then
                i1=indrgd(rgdind(1,rgdtyp),irgd)
                i2=indrgd(rgdind(2,rgdtyp),irgd)

                x(1)=xxx(i1)-xxx(i2)
                y(1)=yyy(i1)-yyy(i2)
                z(1)=zzz(i1)-zzz(i2)

                Call images(imcon,cell,1,x,y,z)

                tmp=Sign( Sqrt((wxx**2+wyy**2+wzz**2)/(x(1)**2+y(1)**2+z(1)**2)) , &
                          Nearest((x(1)*wxx+y(1)*wyy+z(1)*wzz),-1.0_wp) )

                wxx=x(1)*tmp
                wyy=y(1)*tmp
                wzz=z(1)*tmp
             End If

  ! angular velocity in body fixed frame

             rgdoxx(irgd)=(rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rgdrix(2,rgdtyp)
             rgdoyy(irgd)=(rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rgdriy(2,rgdtyp)
             rgdozz(irgd)=(rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rgdriz(2,rgdtyp)

             Do jrgd=1,lrgd
                If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                   i=indrgd(jrgd,irgd) ! local index of particle/site

                   If (i <= natms) Then

  ! site velocity in body frame

                      x(1)=rgdx(jrgd,rgdtyp)
                      y(1)=rgdy(jrgd,rgdtyp)
                      z(1)=rgdz(jrgd,rgdtyp)

                      wxx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                      wyy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                      wzz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

  ! new atomic velocities in lab frame

                      vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+rgdvxx(irgd)
                      vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+rgdvyy(irgd)
                      vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+rgdvzz(irgd)

                   End If
                End If
             End Do
          End If
       End If
    End Do

    Deallocate (gxx,gyy,gzz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_q_ench deallocation failure'
       Call error(0,message)
    End If
  End Subroutine rigid_bodies_q_ench

  Subroutine rigid_bodies_setup(l_str,l_top,megatm,megfrz,megrgd,degtra,degrot,comm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for constructing RBs' rotational inertia tesnors
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( In    ) :: l_str,l_top
    Integer,            Intent( In    ) :: megatm
    Integer,            Intent( InOut ) :: megfrz,megrgd
    Integer(Kind=li),   Intent( InOut ) :: degtra,degrot
    Type( comms_type ), Intent( InOut ) :: comm

    Logical           :: l_print,safe,pass1,pass2
    Integer           :: fail(1:2),irgd,jrgd,krgd,lrgd,rgdtyp, &
                         i,ill,i1,i2,i3, nsite,itmols,nrigid,frzrgd, &
                         ifrz, rotrgd,trargd,iatm1,isite1,ntmp
    Real( Kind = wp ) :: rcut,tmp,weight,                            &
                         rotinr(1:3,1:3),rot1(1:3,1:3),rot(1:9),     &
                         rotall,rotxyz,aa(1:9),bb(1:9),det,rsq,dettest

    Integer,           Allocatable :: allrgd(:),fstrgd(:),lstsit(:)
    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Real( Kind = wp ), Allocatable :: buffer(:)
    Character ( Len = 256 )        :: message,messages(2)

    fail = 0 ; ntmp = mxlrgd*Max(mxrgd,mxtrgd)
    Allocate (allrgd(1:mxtrgd),fstrgd(1:mxrgd),      Stat = fail(1))
    Allocate (gxx(1:ntmp),gyy(1:ntmp), gzz(1:ntmp),  Stat = fail(2))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'rigid_bodies_setup allocation failure'
       Call error(0,message)
    End If


  ! Recover/localise rcut

    rcut=rgdrct

  ! Initialise safety flag

    safe=.true.

  ! decide on detailed printing

    l_print = (l_str .and. l_top)

  ! Tag RBs, find their COMs and check their widths to rcut (system cutoff)

    Call rigid_bodies_tags(comm)
    Call rigid_bodies_coms(xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz,comm)
    Call rigid_bodies_widths(rcut,comm)

  ! Find as many as possible different groups of RB units on this domain
  ! and qualify a representative by the oldest copy of the very first one

    allrgd=0        ! Initialise presence counter (un-encountered yet)
    fstrgd=megatm+1 ! Initialise order of presence (outside particle range)
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       If (allrgd(rgdtyp) == 0) allrgd(rgdtyp)=1

       If (allrgd(rgdtyp) == 1) Then
          i1=indrgd(1,irgd) ! local index of first member
          iatm1=ltg(i1)     ! global index of first member
          If (iatm1 < fstrgd(rgdtyp)) fstrgd(rgdtyp) = iatm1
       End If
    End Do

  ! Loop over all local, only domain present representatives of unique RB types
  ! and get principal axis systems of these RB unit types

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       If (allrgd(rgdtyp) == 1 .and. fstrgd(rgdtyp) == ltg(indrgd(1,irgd))) Then

  ! Get in the local scope of the unit if not fully frozen

          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) < lrgd) Then ! If not fully frozen
             rgdind(0,rgdtyp)=0             ! Not fully frozen (yet)

             Do jrgd=1,lrgd
                krgd=krgd+1

                i=indrgd(jrgd,irgd) ! local index of particle/site

                gxx(krgd)=xxx(i)-rgdxxx(irgd)
                gyy(krgd)=yyy(i)-rgdyyy(irgd)
                gzz(krgd)=zzz(i)-rgdzzz(irgd)
             End Do
          Else                              ! Fully frozen (as if frozen point particle)
             rgdind(0,rgdtyp)=5
             rgdind(1,rgdtyp)=1
             rgdind(2,rgdtyp)=1
             rgdind(3,rgdtyp)=1
          End If
       End If
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Get RB members internal coordinates for these unique RB unit types
  ! that are not fully frozen

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       If (allrgd(rgdtyp) == 1 .and. fstrgd(rgdtyp) == ltg(indrgd(1,irgd))) Then
          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) < lrgd) Then ! If not fully frozen
             rotinr=0.0_wp

             Do jrgd=1,lrgd
                krgd=krgd+1

  ! A slight shortcut here as we care about the diagonaliser only
  ! M = D(diagonal(d,d,d)) +/- R(rest) = d*1 +/- R.  The diagonaliser
  ! C does not care about d*1 & the +/- as C^-1*M*C == d*1 +/- C^-1*R*C

                weight=Real(1-rgdfrz(jrgd,rgdtyp),wp)*rgdwgt(jrgd,rgdtyp)

                rotinr(1,1)=rotinr(1,1)+weight*gxx(krgd)**2
                rotinr(2,1)=rotinr(2,1)+weight*gxx(krgd)*gyy(krgd)
                rotinr(3,1)=rotinr(3,1)+weight*gxx(krgd)*gzz(krgd)
                rotinr(2,2)=rotinr(2,2)+weight*gyy(krgd)**2
                rotinr(3,2)=rotinr(3,2)+weight*gyy(krgd)*gzz(krgd)
                rotinr(3,3)=rotinr(3,3)+weight*gzz(krgd)**2
             End Do
             rotinr(1,2)=rotinr(2,1)
             rotinr(1,3)=rotinr(3,1)
             rotinr(2,3)=rotinr(3,2)

  ! Diagonalise to get eigen values and vectors

             Call jacobi(3,rotinr,rot1)

             rot(1)=rot1(1,1)
             rot(2)=rot1(1,2)
             rot(3)=rot1(1,3)
             rot(4)=rot1(2,1)
             rot(5)=rot1(2,2)
             rot(6)=rot1(2,3)
             rot(7)=rot1(3,1)
             rot(8)=rot1(3,2)
             rot(9)=rot1(3,3)

             krgd=krgd-lrgd
             Do jrgd=1,lrgd
                krgd=krgd+1

                rgdx(jrgd,rgdtyp)=rot(1)*gxx(krgd)+rot(4)*gyy(krgd)+rot(7)*gzz(krgd)
                rgdy(jrgd,rgdtyp)=rot(2)*gxx(krgd)+rot(5)*gyy(krgd)+rot(8)*gzz(krgd)
                rgdz(jrgd,rgdtyp)=rot(3)*gxx(krgd)+rot(6)*gyy(krgd)+rot(9)*gzz(krgd)
             End Do
          End If
       End If
    End Do

  ! OUT OF DOMAIN SCOPE & DIVE IN GLOBAL SCOPE

  ! Globalise internal coordinates and rgdind for all RB unit types
  ! by broadcasting from the lowest rank processor that holds a copy
  ! of the original representative via global summation.

    Allocate (buffer(1:mxtrgd*(4+3*mxlrgd)), Stat = fail(1))
    If (fail(1) > 0) Then
       Write(message,'(a)') 'rigid_bodies_setup allocation failure 1'
       Call error(0,message)
    End If

    krgd=0
    Do irgd=1,mxtrgd
       lrgd=lstrgd(0,irgd)

       iatm1=fstrgd(irgd)
       Call gmin(comm,iatm1)

       If (allrgd(irgd) == 1 .and. fstrgd(irgd) == iatm1) Then
          ntmp=comm%idnode
       Else
          ntmp=comm%mxnode
       End If
       Call gmin(comm,ntmp)

       If (comm%idnode == ntmp) Then
          buffer(krgd+1)=Real(rgdind(0,irgd),wp)
          buffer(krgd+2)=Real(rgdind(1,irgd),wp)
          buffer(krgd+3)=Real(rgdind(2,irgd),wp)
          buffer(krgd+4)=Real(rgdind(3,irgd),wp)
          krgd=krgd+4

          Do jrgd=1,lrgd
             buffer(krgd+1)=rgdx(jrgd,irgd)
             buffer(krgd+2)=rgdy(jrgd,irgd)
             buffer(krgd+3)=rgdz(jrgd,irgd)
             krgd=krgd+3
          End Do
       Else
          buffer(krgd+1)=0.0_wp
          buffer(krgd+2)=0.0_wp
          buffer(krgd+3)=0.0_wp
          buffer(krgd+4)=0.0_wp
          krgd=krgd+4

          Do jrgd=1,lrgd
             buffer(krgd+1)=0.0_wp
             buffer(krgd+2)=0.0_wp
             buffer(krgd+3)=0.0_wp
             krgd=krgd+3
          End Do
       End If
    End Do

    Call gsum(comm,buffer(1:krgd))

    krgd=0
    Do irgd=1,mxtrgd
       rgdind(0,irgd)=Nint(buffer(krgd+1))
       rgdind(1,irgd)=Nint(buffer(krgd+2))
       rgdind(2,irgd)=Nint(buffer(krgd+3))
       rgdind(3,irgd)=Nint(buffer(krgd+4))
       krgd=krgd+4

       lrgd=lstrgd(0,irgd)
       Do jrgd=1,lrgd
          rgdx(jrgd,irgd)=buffer(krgd+1)
          rgdy(jrgd,irgd)=buffer(krgd+2)
          rgdz(jrgd,irgd)=buffer(krgd+3)
          krgd=krgd+3
       End Do
    End Do

    Deallocate (buffer, Stat = fail(1))
    If (fail(1) > 0) Then
       Write(message,'(a)') 'rigid_bodies_setup deallocation failure 1'
       Call error(0,message)
    End If

    Deallocate (allrgd,fstrgd, Stat = fail(1))
    Deallocate (gxx,gyy,gzz,   Stat = fail(2))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'rigid_bodies_setup deallocation failure'
       Call error(0,message)
    End If

    Do irgd=1,mxtrgd
       lrgd=lstrgd(0,irgd)
       If (rgdfrz(0,irgd) < lrgd) Then ! If not fully frozen
          Do jrgd=1,lrgd

  ! Impose rounding

             If (Abs(rgdx(jrgd,irgd)) < 1.0e-8_wp) rgdx(jrgd,irgd)=0.0_wp
             If (Abs(rgdy(jrgd,irgd)) < 1.0e-8_wp) rgdy(jrgd,irgd)=0.0_wp
             If (Abs(rgdz(jrgd,irgd)) < 1.0e-8_wp) rgdz(jrgd,irgd)=0.0_wp

  ! rotational inertia tensor of group type

             weight=Real(1-rgdfrz(jrgd,irgd),wp)*rgdwgt(jrgd,irgd)

             rgdrix(1,irgd) = rgdrix(1,irgd) + weight*(rgdy(jrgd,irgd)**2+rgdz(jrgd,irgd)**2)
             rgdriy(1,irgd) = rgdriy(1,irgd) + weight*(rgdz(jrgd,irgd)**2+rgdx(jrgd,irgd)**2)
             rgdriz(1,irgd) = rgdriz(1,irgd) + weight*(rgdx(jrgd,irgd)**2+rgdy(jrgd,irgd)**2)

          End Do

  ! set axis system such that: Ixx >= Iyy >= Izz

          rotxyz=Max(rgdrix(1,irgd),rgdriy(1,irgd),rgdriz(1,irgd))

          If (rotxyz >= rgdrix(1,irgd)) Then
             If (rotxyz <= rgdriy(1,irgd)) Then
                Do jrgd=1,lrgd
                   tmp=rgdx(jrgd,irgd)
                   rgdx(jrgd,irgd)=rgdy(jrgd,irgd)
                   rgdy(jrgd,irgd)=-tmp
                End Do
                rgdriy(1,irgd)=rgdrix(1,irgd)
                rgdrix(1,irgd)=rotxyz
             Else If (rotxyz <= rgdriz(1,irgd)) Then
                Do jrgd=1,lrgd
                   tmp=rgdx(jrgd,irgd)
                   rgdx(jrgd,irgd)=rgdz(jrgd,irgd)
                   rgdz(jrgd,irgd)=-tmp
                End Do
                rgdriz(1,irgd)=rgdrix(1,irgd)
                rgdrix(1,irgd)=rotxyz
             End If
          End If

          If (rgdriz(1,irgd) > rgdriy(1,irgd)) Then
             Do jrgd=1,lrgd
                tmp=rgdy(jrgd,irgd)
                rgdy(jrgd,irgd)=rgdz(jrgd,irgd)
                rgdz(jrgd,irgd)=-tmp
             End Do
             tmp=rgdriz(1,irgd)
             rgdriz(1,irgd)=rgdriy(1,irgd)
             rgdriy(1,irgd)=tmp
          End If

          rotall=rgdrix(1,irgd)+rgdriy(1,irgd)+rgdriz(1,irgd)
          If (rotall <= 1.0e-5_wp) rotall=1.0_wp

  ! test for type of unit (point/linear/bulk RB == ill=2/1/0)
  ! and get reciprocal of RI in RB unit internal frame of axis

          ill=0
          If (rgdrix(1,irgd)/rotall < 1.0e-5_wp) Then
             ill=ill+1
          Else
             rgdrix(2,irgd)=1.0_wp/rgdrix(1,irgd)
          End If
          If (rgdriy(1,irgd)/rotall < 1.0e-5_wp) Then
             ill=ill+1
          Else
             rgdriy(2,irgd)=1.0_wp/rgdriy(1,irgd)
          End If
          If (rgdriz(1,irgd)/rotall < 1.0e-5_wp) Then
             ill=ill+1
          Else
             rgdriz(2,irgd)=1.0_wp/rgdriz(1,irgd)
          End If

          rgdind(0,irgd)=ill

          If (ill > 1) Then

  ! point molecules and one particle RBs are not allowed by default!
  ! also, partly frozen RBs with only massless unfrozen particles are
  ! caught in read_field!!!

             safe=.false.
             Exit

          Else If (ill == 1) Then

             If      (rgdfrz(0,irgd) == 0) Then

  ! linear unfrozen molecule

                rgdind(1,irgd)=1
                rgdind(2,irgd)=2

             Else If (rgdfrz(0,irgd) >  1) Then

  ! RB with 2+ frozen sites in line (not possible for 1 frozen site only)

                i=0
                Do jrgd=1,lrgd
                   If (rgdfrz(jrgd,irgd) == 1) Then
                      i=i+1
                      rgdind(i,irgd)=jrgd
                      If (i == 3) Exit
                   End If
                End Do

             End If

             If (rgdind(3,irgd) == 0) rgdind(3,irgd)=rgdind(1,irgd)

             i1=rgdind(1,irgd)
             i2=rgdind(2,irgd)

             aa(1)=rgdx(i1,irgd)-rgdx(i2,irgd)
             aa(4)=rgdy(i1,irgd)-rgdy(i2,irgd)
             aa(7)=rgdz(i1,irgd)-rgdz(i2,irgd)
             rsq=Sqrt(aa(1)**2+aa(4)**2+aa(7)**2)

             If      (Abs(aa(7)/rsq) > 0.5_wp) Then
                rsq=Sqrt(aa(4)**2+aa(7)**2)
                aa(2)= 0.0_wp
                aa(5)= aa(7)/rsq
                aa(8)=-aa(4)/rsq
             Else If (Abs(aa(4)/rsq) > 0.5_wp) Then
                rsq=Sqrt(aa(4)**2+aa(1)**2)
                aa(2)=-aa(4)/rsq
                aa(5)= aa(1)/rsq
                aa(8)= 0.0_wp
             Else If (Abs(aa(1)/rsq) > 0.5_wp) Then
                rsq=Sqrt(aa(1)**2+aa(7)**2)
                aa(2)=-aa(7)/rsq
                aa(5)= 0.0_wp
                aa(8)= aa(1)/rsq
             End If

             aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
             aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
             aa(9)=aa(1)*aa(5)-aa(4)*aa(2)

             Call invert(aa,bb,det)

  ! Check aa validity

             If (rgdfrz(0,irgd) /= 1 .and. Abs(det) < 1.0e-5_wp) Then
                safe=.false.
                Exit
             End If

  ! Store tensor

             Do i=1,9
               rgdaxs(i,irgd)=bb(i)
             End Do

          Else If (ill == 0) Then

  ! (1) non-linear molecule or (2) RB with 3+ frozen sites
  ! as at least 3 not in line (as if all were)

             If (rgdfrz(0,irgd) > 1) Then
                rgdind(0,irgd)=4
                rgdind(1,irgd)=1
                rgdind(2,irgd)=1
                rgdind(3,irgd)=1
             Else
                i1=1
                i2=1
                i3=1

                pass1=.true.
                dettest=1.0e-1_wp

                Do While (pass1 .and. i2 < lrgd-1)

                   i2=i2+1
                   i3=i2
                   pass2=.true.

                   Do While (pass2 .and. i3 < lrgd)

                      i3=i3+1

                      aa(1)=rgdx(i1,irgd)-rgdx(i2,irgd)
                      aa(4)=rgdy(i1,irgd)-rgdy(i2,irgd)
                      aa(7)=rgdz(i1,irgd)-rgdz(i2,irgd)

                      aa(2)=rgdx(i1,irgd)-rgdx(i3,irgd)
                      aa(5)=rgdy(i1,irgd)-rgdy(i3,irgd)
                      aa(8)=rgdz(i1,irgd)-rgdz(i3,irgd)

                      aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
                      aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
                      aa(9)=aa(1)*aa(5)-aa(4)*aa(2)

  ! invert matrix

                      Call invert(aa,bb,det)

  ! check on size of determinant - to see if the 3 sites are
  ! too close to being linear for safety.

                      pass2=(Abs(det) < dettest)

                   End Do

                   pass1=(Abs(det) < dettest)

                End Do

  ! Check aa validity

                If (Abs(det) < dettest) Then
                   safe=.false.
                   Exit
                End If

  ! store indices used

                rgdind(1,irgd)=i1
                rgdind(2,irgd)=i2
                rgdind(3,irgd)=i3

  ! store coefficients

                Do i=1,9
                  rgdaxs(i,irgd)=bb(i)
                End Do
             End If

          End If
       End If
    End Do

  ! Report the naughty RB type, number of unit in molecule type and abort

    If (.not.safe) Then
       nrigid=0
    S: Do itmols=1,ntpmls
          Do i=1,numrgd(itmols)
             nrigid=nrigid+1
             If (nrigid == irgd) Exit S
          End Do
       End Do S

       Call warning(307,Real(irgd,wp),Real(itmols,wp),Real(ill,wp))
       Call error(650)
    End If

  ! Sort out degrees of freedom (redundancy & corrections)
  ! correct frzsit,dofsit,rgdwgt,rgdwg1 if needed

    Allocate (lstsit(0:mxlrgd*mxtrgd), Stat = fail(1))
    If (fail(1) > 0) Then
       Write(message,'(a)') 'rigid_bodies_setup allocation failure 2'
       Call error(0,message)
    End If
    lstsit=0

  ! initialise rotational and translational DoF

    degtra=Int(0,li)
    degrot=Int(0,li)

    nsite =0
    nrigid=0
    Do itmols=1,ntpmls
       ifrz=0

       frzrgd=0
       trargd=0
       rotrgd=0

       Do irgd=1,numrgd(itmols)
          nrigid=nrigid+1

          lrgd=lstrgd(0,nrigid)
          ill=rgdind(0,nrigid)

          If (ill == 1) Then

             If (rgdfrz(0,nrigid) == 0) Then

  ! linear molecules lose one axis of rotation but com moves about

                rotrgd=rotrgd+2
                trargd=trargd+3

             Else

  ! not fully frozen RB with 2+ frozen sites in line have
  ! no COM momentum and rotate only around 1 axis

                rotrgd=rotrgd+1

             End If

          Else If (ill == 0) Then

             If (rgdfrz(0,nrigid) == 0) Then

  ! proper unfrozen RB with 3 rot DoFs (rot axis)
  ! and 3 tra DoF (COM moves about)

                trargd=trargd+3
                rotrgd=rotrgd+3

              Else

  ! not fully frozen RB with 1 frozen site
  ! no COM momentum with 3 rot DoFs (rot axis)

                rotrgd=rotrgd+3

              End If

          Else If (ill == 4) Then

  ! As if fully frozen RBs - must get fully frozen then

             rgdind(0,nrigid)=5
             rgdind(1,rgdtyp)=1
             rgdind(2,rgdtyp)=1
             rgdind(3,rgdtyp)=1

             rgdfrz(0,nrigid)=lrgd
             rgdwgt(0,nrigid)=0.0_wp
             rgdwg1(0,nrigid)=Real(lrgd,wp)

             rgdx(:,nrigid)=0.0_wp
             rgdy(:,nrigid)=0.0_wp
             rgdz(:,nrigid)=0.0_wp

             rgdrix(:,nrigid)=0.0_wp
             rgdriy(:,nrigid)=0.0_wp
             rgdriz(:,nrigid)=0.0_wp

             Call warning(305,Real(irgd,wp),Real(itmols,wp),0.0_wp)

             frzrgd=frzrgd+1

             Do jrgd=1,lrgd
                iatm1=lstrgd(jrgd,nrigid)
                isite1=nsite+iatm1

                If (frzsit(isite1) == 0) Then
                   ifrz=ifrz+1

                   frzsit(isite1)=1
                   dofsit(isite1)=0.0_wp

                   lstsit(0)=lstsit(0)+1
                   lstsit(lstsit(0))=isite1
                End If
             End Do

          End If
       End Do

       megfrz=megfrz+ifrz*nummols(itmols)

       megrgd=megrgd-frzrgd*nummols(itmols)
       degtra=degtra+Int(nummols(itmols),li)*Int(trargd,li)
       degrot=degrot+Int(nummols(itmols),li)*Int(rotrgd,li)

       nsite=nsite+numsit(itmols)
    End Do

  ! In case of any refreezing changes refresh the local list of frozen atoms

    If (lstsit(0) > 0) Then
       Do i=1,nlast
          lfrzn(i)=frzsit(lsite(i))
       End Do
    End If

    Deallocate (lstsit, Stat = fail(1))
    If (fail(1) > 0) Then
       Write(message,'(a)') 'rigid_bodies_setup deallocation failure 2'
       Call error(0,message)
    End If

  ! summarise results

    If (comm%idnode==0 .and. l_print) Then
       Call info('summary of rigid body set up',.true.)

       nrigid=0
       Do itmols=1,ntpmls
          Write(message,'(2x,a,i6)') 'in molecule',itmols
          Call info(message,.true.)

          If (numrgd(itmols) == 0) Then
            Write(message,'(2x,a)') 'no rigid bodies specified'
            Call info(message,.true.)
          End If

          Do i=1,numrgd(itmols)
             If (Mod(nrigid,5) == 0) Then
               Write(messages(1),'(2x,a)') &
                 'type :: members :: frozen status :: unfrozen mass :: ' &
                 //'translational DoF :: rotational DoF'
               Write(messages(2),'(4x,a)') &
                 'rotational inertia:        x                   y                   z'
               Call info(messages,2,.true.)
             End If
             nrigid=nrigid+1

             lrgd=lstrgd(0,nrigid)
             ifrz=rgdfrz(0,nrigid)
             ill =rgdind(0,nrigid)

             If      (ifrz < lrgd) Then
                If      (ill == 1) Then ! Linear RB
                   If      (ifrz == 0) Then
                      trargd=3
                      rotrgd=2
                   Else If (ifrz >  1) Then
                      trargd=0
                      rotrgd=1
                   End If
                Else If (ill == 0) Then ! Proper RB
                   If      (ifrz == 0) Then
                      trargd=3
                      rotrgd=3
                   Else If (ifrz == 1) Then
                      trargd=0
                      rotrgd=3
                   End If
                End If
             Else If (ifrz == lrgd) Then
                trargd=0
                rotrgd=0
             End If

             Write(messages(1),'(i5,2x,i6,9x,i6,9x,f13.6,8x,i6,14x,i6)') &
               nrigid,lrgd,ifrz,rgdwgt(0,nrigid),trargd,rotrgd
             Write(messages(2),'(18x,3f20.10)') &
               rgdrix(1,nrigid),rgdriy(1,nrigid),rgdriz(1,nrigid)
             Call info(messages,2,.true.)
             If (lrgd > ifrz) Then
               Write(message,'(6x,a)') &
                 'member  ::   coordinates:        x                   y                   z'
               Call info(message,.true.)
                Do jrgd=1,lrgd
                  Write(message,'(3x,i6,17x,3f20.10)') &
                    jrgd,rgdx(jrgd,nrigid),rgdy(jrgd,nrigid),rgdz(jrgd,nrigid)
                  Call info(message,.true.)
                End Do
             End If
          End Do
       End Do
    End If
    If (l_print) Then
      Call info('',.true.)
    End If

  ! equalise sites DoF due to participating in a good RB (not fully frozen)

    nsite =0
    nrigid=0
    Do itmols=1,ntpmls
       Do irgd=1,numrgd(itmols)
          nrigid=nrigid+1

          lrgd=lstrgd(0,nrigid)
          ill=rgdind(0,nrigid)

  ! 6 = 3(rot) + 3(tra) DoF per RB

          If (ill == 1) Then

             ntmp=0
             Do jrgd=1,lrgd
                iatm1=lstrgd(jrgd,nrigid)
                isite1=nsite+iatm1

                If (dofsit(isite1) > zero_plus) ntmp=ntmp+1
             End Do

             If (rgdfrz(0,nrigid) == 0) Then

  ! linear molecule losing one axis of rotation - losing 1(rot) DoF

                Do jrgd=1,lrgd
                   iatm1=lstrgd(jrgd,nrigid)
                   isite1=nsite+iatm1

                   If (dofsit(isite1) > zero_plus) dofsit(isite1)=5.0_wp/Real(ntmp,wp)
                End Do

             Else

  ! RB with 2+ frozen sites in line with members restricted to
  ! a circular line around one axis - losing 2(rot) & 3(tra) DoF

                Do jrgd=1,lrgd
                   iatm1=lstrgd(jrgd,nrigid)
                   isite1=nsite+iatm1

                   If (dofsit(isite1) > zero_plus) dofsit(isite1)=1.0_wp/Real(ntmp,wp)
                End Do

             End If

          Else If (ill == 0) Then

             ntmp=0
             Do jrgd=1,lrgd
                iatm1=lstrgd(jrgd,nrigid)
                isite1=nsite+iatm1

                If (dofsit(isite1) > zero_plus) ntmp=ntmp+1
             End Do

             If (rgdfrz(0,nrigid) == 0) Then

  ! Proper RB

                Do jrgd=1,lrgd
                   iatm1=lstrgd(jrgd,nrigid)
                   isite1=nsite+iatm1

                   If (dofsit(isite1) > zero_plus) dofsit(isite1)=6.0_wp/Real(ntmp,wp)
                End Do

             Else If (rgdfrz(0,nrigid) == 1) Then

  ! RB with 1 frozen site with members restricted
  ! to a spherical surface - losing 3(tra) DoF

                Do jrgd=1,lrgd
                   iatm1=lstrgd(jrgd,nrigid)
                   isite1=nsite+iatm1

                   If (dofsit(isite1) > zero_plus) dofsit(isite1)=3.0_wp/Real(ntmp,wp)
                End Do

             Else

                Do jrgd=1,lrgd
                   iatm1=lstrgd(jrgd,nrigid)
                   isite1=nsite+iatm1

                   If (dofsit(isite1) > zero_plus) safe=.false.
                End Do

             End If

          End If
       End Do

  ! Report the naughty RB type, number of unit in molecule type and abort

       If (.not.safe) Then
          Call warning(307,Real(irgd,wp),Real(itmols,wp),Real(ill,wp))
          Call error(644)
       End If

       nsite=nsite+numsit(itmols)
    End Do

  ! OUT OF GLOBAL SCOPE & BACK IN DOMAIN SCOPE

  ! set-up quaternions

    Call q_setup(comm)

  End Subroutine rigid_bodies_setup

  Subroutine rigid_bodies_split_torque(gxx,gyy,gzz,txx,tyy,tzz,uxx,uyy,uzz,comm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for resolving RBs' torques into equivalent atomic
  ! forces suitable for conjugate gradient minimisation
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith may 2006
  ! adapted   - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent( InOut ) :: gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms)
    Real( Kind = wp ),  Intent(   Out ) :: txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms), &
                                           uxx(1:mxatms),uyy(1:mxatms),uzz(1:mxatms)
    Type( comms_type ), Intent( In    ) :: comm

    Integer           :: fail,i,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ) :: fmx,fmy,fmz,x(1:1),y(1:1),z(1:1),tqx,tqy,tqz,torque,tmp,qqq

    Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
    Character ( Len = 256 )        :: message

    fail = 0
    Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_split_torque allocation failure'
       Call error(0,message)
    End If

  ! Get a RB particles vectors wrt the RB's COM

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

  ! For all good RBs

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
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

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

  ! For all good RBs

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then

  ! calculate net force and torque on the RB - assumption must hold here

          fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
          tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp

          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

             If (rgdfrz(0,rgdtyp) == 0) Then
                fmx=fmx+gxx(i)
                fmy=fmy+gyy(i)
                fmz=fmz+gzz(i)
             End If

             tqx=tqx+ggy(krgd)*gzz(i)-ggz(krgd)*gyy(i)
             tqy=tqy+ggz(krgd)*gxx(i)-ggx(krgd)*gzz(i)
             tqz=tqz+ggx(krgd)*gyy(i)-ggy(krgd)*gxx(i)
          End Do

  ! offload new forces on RB members

          If (rgdfrz(0,rgdtyp) == 0) Then
             tmp=1.0_wp/Real(lrgd,wp)

             fmx=fmx*tmp
             fmy=fmy*tmp
             fmz=fmz*tmp
          Else
             fmx=0.0_wp
             fmy=0.0_wp
             fmz=0.0_wp
          End If

          Do jrgd=1,lrgd
             i=indrgd(jrgd,irgd) ! local index of particle/site

             gxx(i)=fmx
             gyy(i)=fmy
             gzz(i)=fmz
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

  ! Get magnitude of torque

          torque=Sqrt(tqx**2+tqy**2+tqz**2)

  ! Calculate unit vectors for new site forces

          krgd=krgd-lrgd
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site
             If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                txx(i)=ggy(krgd)*tqz-tqy*ggz(krgd)
                tyy(i)=ggz(krgd)*tqx-tqz*ggx(krgd)
                tzz(i)=ggx(krgd)*tqy-tqx*ggy(krgd)

                tmp=Sqrt(txx(i)**2+tyy(i)**2+tzz(i)**2)
                If (tmp > 1.0e-10_wp) Then
                   tmp=1.0_wp/tmp
                   txx(i)=txx(i)*tmp
                   tyy(i)=tyy(i)*tmp
                   tzz(i)=tzz(i)*tmp
                Else
                   txx(i)=0.0_wp
                   tyy(i)=0.0_wp
                   tzz(i)=0.0_wp
                End If
             Else
                txx(i)=0.0_wp
                tyy(i)=0.0_wp
                tzz(i)=0.0_wp
             End If
          End Do

  ! construct unit vectors for new radial positions of
  ! RB's sites with respect to the axis of rotation

          tmp=1.0_wp
          If (torque > 1.0e-10_wp) tmp=1.0_wp/torque

          Do jrgd=1,lrgd
             If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                i=indrgd(jrgd,irgd) ! local index of particle/site

                uxx(i)=(tyy(i)*tqz-tqy*tzz(i))*tmp
                uyy(i)=(tzz(i)*tqx-tqz*txx(i))*tmp
                uzz(i)=(txx(i)*tqy-tqx*tyy(i))*tmp
             Else
                uxx(i)=0.0_wp
                uyy(i)=0.0_wp
                uzz(i)=0.0_wp
             End If
          End Do

  ! scale unit vectors to working lengths

          qqq=0.0_wp

          krgd=krgd-lrgd
          Do jrgd=1,lrgd
             krgd=krgd+1

             If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                i=indrgd(jrgd,irgd) ! local index of particle/site

                tmp=ggx(krgd)*uxx(i)+ggy(krgd)*uyy(i)+ggz(krgd)*uzz(i)

                txx(i)=txx(i)*tmp
                tyy(i)=tyy(i)*tmp
                tzz(i)=tzz(i)*tmp

                uxx(i)=uxx(i)*tmp
                uyy(i)=uyy(i)*tmp
                uzz(i)=uzz(i)*tmp

                qqq=qqq+tmp**2
             End If
          End Do

  ! scale new site forces so that all RB members rotate
  ! to the same angle in their plains of rotation

          tmp=0.0_wp
          If (qqq > 1.0e-10_wp) tmp=torque/qqq

          Do jrgd=1,lrgd
             If (rgdfrz(jrgd,rgdtyp) == 0) Then ! Apply restrictions
                i=indrgd(jrgd,irgd) ! local index of particle/site

                txx(i)=txx(i)*tmp
                tyy(i)=tyy(i)*tmp
                tzz(i)=tzz(i)*tmp
             End If
          End Do

       Else

  ! Initialise unit vectors for new site forces and radial location
  ! in the plane of rotation wrt the axis of rotation (torque vector)

          Do jrgd=1,lrgd
             i=indrgd(jrgd,irgd) ! local index of particle/site

             txx(i)=0.0_wp ; tyy(i)=0.0_wp ; tzz(i)=0.0_wp
             uxx(i)=0.0_wp ; uyy(i)=0.0_wp ; uzz(i)=0.0_wp
             gxx(i)=0.0_wp ; gyy(i)=0.0_wp ; gzz(i)=0.0_wp
          End Do

       End If
    End Do

    Deallocate (ggx,ggy,ggz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_split_torque deallocation failure'
       Call error(0,message)
    End If
  End Subroutine rigid_bodies_split_torque

  Subroutine rigid_bodies_stress(strcom,ggx,ggy,ggz,comm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to calculate RB contributions to the atomic stress
  ! tensor
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov july 2013
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)
    Real( Kind = wp ), Intent( In    ) :: ggx(1:mxlrgd*mxrgd), &
                                          ggy(1:mxlrgd*mxrgd), &
                                          ggz(1:mxlrgd*mxrgd)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: i,irgd,jrgd,krgd,lrgd,rgdtyp

  ! Initialise stress

    strcom=0.0_wp

  ! convert atomic virial to molecular
  ! note convention: virial(atom-atom) = -sum(Ri.Fi)
  ! : virial(com-com) = -sum(Rcom.Fcom) so
  ! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

             If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
                strcom(1)=strcom(1)-ggx(krgd)*fxx(i)
                strcom(2)=strcom(2)-ggx(krgd)*fyy(i)
                strcom(3)=strcom(3)-ggx(krgd)*fzz(i)
                strcom(4)=strcom(4)-ggy(krgd)*fxx(i)
                strcom(5)=strcom(5)-ggy(krgd)*fyy(i)
                strcom(6)=strcom(6)-ggy(krgd)*fzz(i)
                strcom(7)=strcom(7)-ggz(krgd)*fxx(i)
                strcom(8)=strcom(8)-ggz(krgd)*fyy(i)
                strcom(9)=strcom(9)-ggz(krgd)*fzz(i)
             End If
          End Do
       End If
    End Do

    Call gsum(comm,strcom)

  ! Symmetrise

     strcom(2)=0.5_wp*(strcom(2)+strcom(4))
     strcom(4)=strcom(2)
     strcom(3)=0.5_wp*(strcom(3)+strcom(7))
     strcom(7)=strcom(3)
     strcom(6)=0.5_wp*(strcom(6)+strcom(8))
     strcom(8)=strcom(6)

  End subroutine rigid_bodies_stress

  Subroutine rigid_bodies_stre_s(strcom,ggx,ggy,ggz,fxx,fyy,fzz,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to calculate RB contributions to the atomic stress
  ! tensor
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov july 2013
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent(   Out ) :: strcom(1:9)
    Real( Kind = wp ),  Intent( In    ) :: ggx(1:mxlrgd*mxrgd), &
                                           ggy(1:mxlrgd*mxrgd), &
                                           ggz(1:mxlrgd*mxrgd)
    Real( Kind = wp ),  Intent( In    ) :: fxx(1:mxatms),fyy(1:mxatms),fzz(1:mxatms)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: i,irgd,jrgd,krgd,lrgd,rgdtyp

  ! Initialise stress

    strcom=0.0_wp

  ! convert atomic virial to molecular
  ! note convention: virial(atom-atom) = -sum(Ri.Fi)
  ! : virial(com-com) = -sum(Rcom.Fcom) so
  ! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

             If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
                strcom(1)=strcom(1)-ggx(krgd)*fxx(i)
                strcom(2)=strcom(2)-ggx(krgd)*fyy(i)
                strcom(3)=strcom(3)-ggx(krgd)*fzz(i)
                strcom(4)=strcom(4)-ggy(krgd)*fxx(i)
                strcom(5)=strcom(5)-ggy(krgd)*fyy(i)
                strcom(6)=strcom(6)-ggy(krgd)*fzz(i)
                strcom(7)=strcom(7)-ggz(krgd)*fxx(i)
                strcom(8)=strcom(8)-ggz(krgd)*fyy(i)
                strcom(9)=strcom(9)-ggz(krgd)*fzz(i)
             End If
          End Do
       End If
    End Do

    Call gsum(comm,strcom)

  ! Symmetrise

     strcom(2)=0.5_wp*(strcom(2)+strcom(4))
     strcom(4)=strcom(2)
     strcom(3)=0.5_wp*(strcom(3)+strcom(7))
     strcom(7)=strcom(3)
     strcom(6)=0.5_wp*(strcom(6)+strcom(8))
     strcom(8)=strcom(6)

  End Subroutine rigid_bodies_stre_s

  Subroutine rigid_bodies_str_ss(strcom,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to calculate RB contributions to the atomic stress
  ! tensor
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent(   Out ) :: strcom(1:9)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: fail,i,irgd,jrgd,krgd,lrgd,rgdtyp

    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Character ( Len = 256 )        :: message

    fail = 0
    Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_stress allocation failure'
       Call error(0,message)
    End If

  ! Loop over all local RB units and get in the local scope of the unit

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

  ! COM distances

             gxx(krgd)=xxx(i)-rgdxxx(irgd)
             gyy(krgd)=yyy(i)-rgdyyy(irgd)
             gzz(krgd)=zzz(i)-rgdzzz(irgd)
          End Do
       End If
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Initialise stress

    strcom=0.0_wp

  ! convert atomic virial to molecular
  ! note convention: virial(atom-atom) = -sum(Ri.Fi)
  ! : virial(com-com) = -sum(Rcom.Fcom) so
  ! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

             If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
                strcom(1)=strcom(1)-gxx(krgd)*fxx(i)
                strcom(2)=strcom(2)-gxx(krgd)*fyy(i)
                strcom(3)=strcom(3)-gxx(krgd)*fzz(i)
                strcom(4)=strcom(4)-gyy(krgd)*fxx(i)
                strcom(5)=strcom(5)-gyy(krgd)*fyy(i)
                strcom(6)=strcom(6)-gyy(krgd)*fzz(i)
                strcom(7)=strcom(7)-gzz(krgd)*fxx(i)
                strcom(8)=strcom(8)-gzz(krgd)*fyy(i)
                strcom(9)=strcom(9)-gzz(krgd)*fzz(i)
             End If
          End Do
       End If
    End Do

    Call gsum(comm,strcom)

  ! Symmetrise

     strcom(2)=0.5_wp*(strcom(2)+strcom(4))
     strcom(4)=strcom(2)
     strcom(3)=0.5_wp*(strcom(3)+strcom(7))
     strcom(7)=strcom(3)
     strcom(6)=0.5_wp*(strcom(6)+strcom(8))
     strcom(8)=strcom(6)

    fail = 0
    Deallocate (gxx,gyy,gzz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_stress deallocation failure'
       Call error(0,message)
    End If

  End subroutine rigid_bodies_str_ss

  Subroutine rigid_bodies_str__s(strcom,fxx,fyy,fzz,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to calculate RB contributions to the atomic stress
  ! tensor
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent(   Out ) :: strcom(1:9)
    Real( Kind = wp ),  Intent( In    ) :: fxx(1:mxatms),fyy(1:mxatms),fzz(1:mxatms)
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: fail,i,irgd,jrgd,krgd,lrgd,rgdtyp

    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Character ( Len = 256 )  :: message
   

    fail = 0
    Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_stress allocation failure'
       Call error(0,message)
    End If

  ! Loop over all local RB units and get in the local scope of the unit

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

  ! COM distances

             gxx(krgd)=xxx(i)-rgdxxx(irgd)
             gyy(krgd)=yyy(i)-rgdyyy(irgd)
             gzz(krgd)=zzz(i)-rgdzzz(irgd)
          End Do
       End If
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Initialise stress

    strcom=0.0_wp

  ! convert atomic virial to molecular
  ! note convention: virial(atom-atom) = -sum(Ri.Fi)
  ! : virial(com-com) = -sum(Rcom.Fcom) so
  ! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

    krgd=0
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=indrgd(jrgd,irgd) ! local index of particle/site

             If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
                strcom(1)=strcom(1)-gxx(krgd)*fxx(i)
                strcom(2)=strcom(2)-gxx(krgd)*fyy(i)
                strcom(3)=strcom(3)-gxx(krgd)*fzz(i)
                strcom(4)=strcom(4)-gyy(krgd)*fxx(i)
                strcom(5)=strcom(5)-gyy(krgd)*fyy(i)
                strcom(6)=strcom(6)-gyy(krgd)*fzz(i)
                strcom(7)=strcom(7)-gzz(krgd)*fxx(i)
                strcom(8)=strcom(8)-gzz(krgd)*fyy(i)
                strcom(9)=strcom(9)-gzz(krgd)*fzz(i)
             End If
          End Do
       End If
    End Do

    Call gsum(comm,strcom)

  ! Symmetrise

     strcom(2)=0.5_wp*(strcom(2)+strcom(4))
     strcom(4)=strcom(2)
     strcom(3)=0.5_wp*(strcom(3)+strcom(7))
     strcom(7)=strcom(3)
     strcom(6)=0.5_wp*(strcom(6)+strcom(8))
     strcom(8)=strcom(6)

    fail = 0
    Deallocate (gxx,gyy,gzz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_stress deallocation failure'
       Call error(0,message)
    End If
  End subroutine rigid_bodies_str__s

  Subroutine rigid_bodies_tags(comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for identifying and indexing RB units
  !
  ! Note: must be used in conjunction with integration algorithms
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( comms_type ), Intent( InOut ) :: comm

    Logical :: safe
    Integer :: fail,irgd,jrgd,lrgd,s,i

    Logical, Allocatable :: lunsafe(:)
    Character ( Len = 256 ) :: message


    fail=0
    Allocate (lunsafe(1:mxrgd), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_tags allocation failure'
       Call error(0,message)
    End If

  ! Loop over all local to this node RB units and save the indices of the members
  ! and their presence on the domain in total

    Do irgd=1,ntrgd
       lunsafe(irgd)=.false.
       lrgd=listrgd(-1,irgd)

  ! Initialise local indices

       indrgd(:,irgd)=0
       Do jrgd=1,lrgd
          s=listrgd(jrgd,irgd)
          i=local_index(s,nlast,lsi,lsa)

          indrgd(jrgd,irgd)=i
          If (i > 0 .and. i <= natms) indrgd(0,irgd)=indrgd(0,irgd)+1
       End Do

  ! Detect uncompressed unit

       If (Any(indrgd(1:lrgd,irgd) == 0) .and. indrgd(0,irgd) > 0) lunsafe(irgd)=.true.
    End Do

  ! Check if a RB unit has a diameter > rcut (the system cutoff)
  ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:ntrgd))
    Call gcheck(comm,safe)
    If (.not.safe) Then
       Do i=0,comm%mxnode-1
          If (comm%idnode == i) Then
             Do irgd=1,ntrgd
                If (lunsafe(irgd)) Then
                  Write(message,'(a,2(i10,a))') &
                    'global unit number', listrgd(0,irgd), &
                    ' , with a head particle number', listrgd(1,irgd), &
                    ' contributes towards next error'
                  Call warning(message)
                End If
             End Do
          End If
          Call gsync(comm)
       End Do
       Call error(642)
    End If

    Deallocate (lunsafe, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_tags deallocation failure, node'
       Call error(0,message)
    End If
  End Subroutine rigid_bodies_tags

  Subroutine rigid_bodies_widths(rcut,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for checking RBs' widths compliency to < rcut
  ! (the system cutoff)
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent( In    ) :: rcut
    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: fail,irgd,jrgd,krgd,lrgd,mrgd,nrgd,rgdtyp
    Real( Kind = wp ) :: d,width

    Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
    Character ( Len = 256 ) ::  message

    fail = 0
    Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_widths allocation failure'
       Call error(0,message)
    End If

  ! Loop over all local RB units and get in the local scope of the unit

    krgd=0
    Do irgd=1,ntrgd
       lrgd=listrgd(-1,irgd)

       Do jrgd=1,lrgd
          krgd=krgd+1

          gxx(krgd) = xxx(indrgd(jrgd,irgd)) - xxx(indrgd(1,irgd))
          gyy(krgd) = yyy(indrgd(jrgd,irgd)) - yyy(indrgd(1,irgd))
          gzz(krgd) = zzz(indrgd(jrgd,irgd)) - zzz(indrgd(1,irgd))
       End Do
    End Do

  ! minimum image convention for bond vectors

    Call images(imcon,cell,krgd,gxx,gyy,gzz)

  ! Get the COM vector and check of diameter safety of units

    krgd=0
    width=0.0_wp
    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       Do jrgd=1,lrgd
          krgd=krgd+1

          Do mrgd=jrgd+1,lrgd
             nrgd=krgd+mrgd-jrgd

             d=Sqrt((gxx(krgd)-gxx(nrgd))**2+(gyy(krgd)-gyy(nrgd))**2+(gzz(krgd)-gzz(nrgd))**2)
             width=Max(width,d)

             If (d > rcut) Then
               Write(message,'(a,3i5,2f8.3)') &
                 'RB type, members(1,2), width > cutoff', &
                 rgdtyp,jrgd,mrgd,width,rcut
               Call info(message)
             End If
          End Do
       End Do
    End Do

  ! Check if a RB unit has a diameter > the cutoff

    Call gmax(comm,width)
    If (width > rcut) Then
       Call warning(8,width,rcut,0.0_wp)
       Call error(642)
    End If

    Deallocate (gxx,gyy,gzz, Stat = fail)
    If (fail > 0) Then
       Write(message,'(a)') 'rigid_bodies_widths deallocation failure'
       Call error(0,message)
    End If
  End Subroutine rigid_bodies_widths

  Subroutine xscale(m_rgd,tstep,thermo,stats,neigh,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to scale initial positions with change in box shape
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: m_rgd
  Real( Kind = wp ), Intent( In    ) :: tstep
  Type( thermostat_type), Intent( InOut ) :: thermo
  Type( stats_type), Intent( InOut ) :: stats
  Type( neighbours_type ), Intent( InOut ) :: neigh
  Type( comms_type), Intent( InOut ) :: comm

  Integer           :: fail,i,j,irgd,jrgd,lrgd
  Real( Kind = wp ) :: a1,a2,a3,a5,a6,a9,b1,b2,b3,b5,b6,b9,scale, &
                       xa,ya,za,x,y,z,com(1:3)

  Real( Kind = wp ), Allocatable :: rgdxin(:),rgdyin(:),rgdzin(:)
  Character ( Len = 256 )        :: message

  If (.not. thermo%variable_cell) Return

  If (m_rgd == 0) Then

     If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! berendsen npt/nst

        If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

           scale = thermo%eta(1)

           Do i=1,natms
              stats%xin(i) = scale*stats%xin(i)
              stats%yin(i) = scale*stats%yin(i)
              stats%zin(i) = scale*stats%zin(i)
           End Do

        Else

           Do i=1,natms
              xa = stats%xin(i)*thermo%eta(1)+stats%yin(i)*thermo%eta(2)+stats%zin(i)*thermo%eta(3)
              ya = stats%xin(i)*thermo%eta(4)+stats%yin(i)*thermo%eta(5)+stats%zin(i)*thermo%eta(6)
              za = stats%xin(i)*thermo%eta(7)+stats%yin(i)*thermo%eta(8)+stats%zin(i)*thermo%eta(9)

              stats%xin(i) = xa
              stats%yin(i) = ya
              stats%zin(i) = za
           End Do

        End If

     Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! hoover npt/nst

        Call getcom(stats%xin,stats%yin,stats%zin,com,comm)

        If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

           scale = Exp(tstep*thermo%eta(1))

           Do i=1,natms
              stats%xin(i) = scale*(stats%xin(i)-com(1))+com(1)
              stats%yin(i) = scale*(stats%yin(i)-com(2))+com(2)
              stats%zin(i) = scale*(stats%zin(i)-com(3))+com(3)
           End Do

        Else

! second order taylor expansion of Exp(tstep*thermo%eta)

           a1 = tstep*thermo%eta(1)
           a2 = tstep*thermo%eta(2)
           a3 = tstep*thermo%eta(3)
           a5 = tstep*thermo%eta(5)
           a6 = tstep*thermo%eta(6)
           a9 = tstep*thermo%eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do i=1,natms
              xa = stats%xin(i)-com(1)
              ya = stats%yin(i)-com(2)
              za = stats%zin(i)-com(3)

              stats%xin(i) = xa*b1 + ya*b2 + za*b3 + com(1)
              stats%yin(i) = xa*b2 + ya*b5 + za*b6 + com(2)
              stats%zin(i) = xa*b3 + ya*b6 + za*b9 + com(3)
           End Do

        End If

     Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
              thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
              thermo%ensemble == ENS_NPT_MTK .or. &
              thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! Langevin and MTK npt/nst

        If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

           scale = Exp(tstep*thermo%eta(1))

           Do i=1,natms
              stats%xin(i) = scale*stats%xin(i)
              stats%yin(i) = scale*stats%yin(i)
              stats%zin(i) = scale*stats%zin(i)
           End Do

        Else

! second order taylor expansion of Exp(tstep*thermo%eta)

           a1 = tstep*thermo%eta(1)
           a2 = tstep*thermo%eta(2)
           a3 = tstep*thermo%eta(3)
           a5 = tstep*thermo%eta(5)
           a6 = tstep*thermo%eta(6)
           a9 = tstep*thermo%eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do i=1,natms
              xa = stats%xin(i)
              ya = stats%yin(i)
              za = stats%zin(i)

              stats%xin(i) = xa*b1 + ya*b2 + za*b3
              stats%yin(i) = xa*b2 + ya*b5 + za*b6
              stats%zin(i) = xa*b3 + ya*b6 + za*b9
           End Do

        End If

     End If

     If (.not.neigh%update) Then

        If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! berendsen npt/nst

           If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

              scale = thermo%eta(1)

              Do i=1,natms
                 neigh%xbg(i) = scale*neigh%xbg(i)
                 neigh%ybg(i) = scale*neigh%ybg(i)
                 neigh%zbg(i) = scale*neigh%zbg(i)
              End Do

           Else

              Do i=1,natms
                 xa = neigh%xbg(i)*thermo%eta(1)+neigh%ybg(i)*thermo%eta(2)+neigh%zbg(i)*thermo%eta(3)
                 ya = neigh%xbg(i)*thermo%eta(4)+neigh%ybg(i)*thermo%eta(5)+neigh%zbg(i)*thermo%eta(6)
                 za = neigh%xbg(i)*thermo%eta(7)+neigh%ybg(i)*thermo%eta(8)+neigh%zbg(i)*thermo%eta(9)

                 neigh%xbg(i) = xa
                 neigh%ybg(i) = ya
                 neigh%zbg(i) = za
              End Do

           End If

        Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! hoover npt/nst

           Call getcom(neigh%xbg,neigh%ybg,neigh%zbg,com,comm)

           If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

              scale = Exp(tstep*thermo%eta(1))

              Do i=1,natms
                 neigh%xbg(i) = scale*(neigh%xbg(i)-com(1))+com(1)
                 neigh%ybg(i) = scale*(neigh%ybg(i)-com(2))+com(2)
                 neigh%zbg(i) = scale*(neigh%zbg(i)-com(3))+com(3)
              End Do

           Else

! second order taylor expansion of Exp(tstep*thermo%eta)

              a1 = tstep*thermo%eta(1)
              a2 = tstep*thermo%eta(2)
              a3 = tstep*thermo%eta(3)
              a5 = tstep*thermo%eta(5)
              a6 = tstep*thermo%eta(6)
              a9 = tstep*thermo%eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do i=1,natms
                 xa = neigh%xbg(i)-com(1)
                 ya = neigh%ybg(i)-com(2)
                 za = neigh%zbg(i)-com(3)

                 neigh%xbg(i) = xa*b1 + ya*b2 + za*b3 + com(1)
                 neigh%ybg(i) = xa*b2 + ya*b5 + za*b6 + com(2)
                 neigh%zbg(i) = xa*b3 + ya*b6 + za*b9 + com(3)
              End Do

           End If

        Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
                 thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
                 thermo%ensemble == ENS_NPT_MTK .or. &
                 thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! Langevin and MTK npt/nst

           If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

              scale = Exp(tstep*thermo%eta(1))

              Do i=1,natms
                 neigh%xbg(i) = scale*neigh%xbg(i)
                 neigh%ybg(i) = scale*neigh%ybg(i)
                 neigh%zbg(i) = scale*neigh%zbg(i)
              End Do

           Else

! second order taylor expansion of Exp(tstep*thermo%eta)

              a1 = tstep*thermo%eta(1)
              a2 = tstep*thermo%eta(2)
              a3 = tstep*thermo%eta(3)
              a5 = tstep*thermo%eta(5)
              a6 = tstep*thermo%eta(6)
              a9 = tstep*thermo%eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do i=1,natms
                 xa = neigh%xbg(i)
                 ya = neigh%ybg(i)
                 za = neigh%zbg(i)

                 neigh%xbg(i) = xa*b1 + ya*b2 + za*b3
                 neigh%ybg(i) = xa*b2 + ya*b5 + za*b6
                 neigh%zbg(i) = xa*b3 + ya*b6 + za*b9
              End Do

           End If

        End If

     End If

  Else ! RBs exist

     fail = 0
     Allocate (rgdxin(1:mxrgd),rgdyin(1:mxrgd),rgdzin(1:mxrgd), Stat = fail)
     If (fail > 0) Then
        Write(message,'(a)') 'xscale allocation failure'
        Call error(0,message)
     End If

! Halo initial RB members positions across onto neighbouring domains
! to get initial COMs

     If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,stats%xin,stats%yin,stats%zin,comm)
     Call rigid_bodies_coms(stats%xin,stats%yin,stats%zin,rgdxin,rgdyin,rgdzin,comm)

     If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! berendsen npt/nst

        If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

           scale = thermo%eta(1)

           Do j=1,nfree
              i=lstfre(j)

              stats%xin(i) = scale*stats%xin(i)
              stats%yin(i) = scale*stats%yin(i)
              stats%zin(i) = scale*stats%zin(i)
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
                    stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                    stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                    stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        Else

           Do j=1,nfree
              i=lstfre(j)

              xa = stats%xin(i)*thermo%eta(1)+stats%yin(i)*thermo%eta(2)+stats%zin(i)*thermo%eta(3)
              ya = stats%xin(i)*thermo%eta(4)+stats%yin(i)*thermo%eta(5)+stats%zin(i)*thermo%eta(6)
              za = stats%xin(i)*thermo%eta(7)+stats%yin(i)*thermo%eta(8)+stats%zin(i)*thermo%eta(9)

              stats%xin(i) = xa
              stats%yin(i) = ya
              stats%zin(i) = za
           End Do

           Do irgd=1,ntrgd
              x = rgdxin(irgd)
              y = rgdyin(irgd)
              z = rgdzin(irgd)

              xa = rgdxin(irgd)*thermo%eta(1)+rgdyin(irgd)*thermo%eta(2)+rgdzin(irgd)*thermo%eta(3)
              ya = rgdxin(irgd)*thermo%eta(4)+rgdyin(irgd)*thermo%eta(5)+rgdzin(irgd)*thermo%eta(6)
              za = rgdxin(irgd)*thermo%eta(7)+rgdyin(irgd)*thermo%eta(8)+rgdzin(irgd)*thermo%eta(9)

              rgdxin(irgd) = xa
              rgdyin(irgd) = ya
              rgdzin(irgd) = za

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd)

                 If (i <= natms) Then
                    stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                    stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                    stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        End If

     Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! hoover npt/nst

        Call getcom(stats%xin,stats%yin,stats%zin,com,comm)

        If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

           scale = Exp(tstep*thermo%eta(1))

           Do j=1,nfree
              i=lstfre(j)

              stats%xin(i) = scale*(stats%xin(i)-com(1))+com(1)
              stats%yin(i) = scale*(stats%yin(i)-com(2))+com(2)
              stats%zin(i) = scale*(stats%zin(i)-com(3))+com(3)
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
                    stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                    stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                    stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        Else

! second order taylor expansion of Exp(tstep*thermo%eta)

           a1 = tstep*thermo%eta(1)
           a2 = tstep*thermo%eta(2)
           a3 = tstep*thermo%eta(3)
           a5 = tstep*thermo%eta(5)
           a6 = tstep*thermo%eta(6)
           a9 = tstep*thermo%eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do j=1,nfree
              i=lstfre(j)

              xa = stats%xin(i)-com(1)
              ya = stats%yin(i)-com(2)
              za = stats%zin(i)-com(3)

              stats%xin(i) = xa*b1 + ya*b2 + za*b3 + com(1)
              stats%yin(i) = xa*b2 + ya*b5 + za*b6 + com(2)
              stats%zin(i) = xa*b3 + ya*b6 + za*b9 + com(3)
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
                    stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                    stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                    stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        End If

     Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
              thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
              thermo%ensemble == ENS_NPT_MTK .or. &
              thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! Langevin and MTK npt/nst

        If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

           scale = Exp(tstep*thermo%eta(1))

           Do j=1,nfree
              i=lstfre(j)

              stats%xin(i) = scale*stats%xin(i)
              stats%yin(i) = scale*stats%yin(i)
              stats%zin(i) = scale*stats%zin(i)
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
                    stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                    stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                    stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        Else

! second order taylor expansion of Exp(tstep*thermo%eta)

           a1 = tstep*thermo%eta(1)
           a2 = tstep*thermo%eta(2)
           a3 = tstep*thermo%eta(3)
           a5 = tstep*thermo%eta(5)
           a6 = tstep*thermo%eta(6)
           a9 = tstep*thermo%eta(9)

           b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
           b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
           b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
           b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
           b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
           b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

           Do j=1,nfree
              i=lstfre(j)

              xa = stats%xin(i)
              ya = stats%yin(i)
              za = stats%zin(i)

              stats%xin(i) = xa*b1 + ya*b2 + za*b3
              stats%yin(i) = xa*b2 + ya*b5 + za*b6
              stats%zin(i) = xa*b3 + ya*b6 + za*b9
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
                    stats%xin(i) = stats%xin(i) - x + rgdxin(irgd)
                    stats%yin(i) = stats%yin(i) - y + rgdyin(irgd)
                    stats%zin(i) = stats%zin(i) - z + rgdzin(irgd)
                 End If
              End Do
           End Do

        End If

     End If

     If (.not.neigh%update) Then

        Call rigid_bodies_coms(neigh%xbg,neigh%ybg,neigh%zbg,rgdxin,rgdyin,rgdzin,comm)

        If (thermo%ensemble == ENS_NPT_BERENDSEN .or. thermo%ensemble == ENS_NPT_BERENDSEN_ANISO) Then

! berendsen npt/nst

           If (thermo%ensemble == ENS_NPT_BERENDSEN) Then

              scale = thermo%eta(1)

              Do j=1,nfree
                 i=lstfre(j)

                 neigh%xbg(i) = scale*neigh%xbg(i)
                 neigh%ybg(i) = scale*neigh%ybg(i)
                 neigh%zbg(i) = scale*neigh%zbg(i)
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
                       neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                       neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                       neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           Else

              Do j=1,nfree
                 i=lstfre(j)

                 xa = neigh%xbg(i)*thermo%eta(1)+neigh%ybg(i)*thermo%eta(2)+neigh%zbg(i)*thermo%eta(3)
                 ya = neigh%xbg(i)*thermo%eta(4)+neigh%ybg(i)*thermo%eta(5)+neigh%zbg(i)*thermo%eta(6)
                 za = neigh%xbg(i)*thermo%eta(7)+neigh%ybg(i)*thermo%eta(8)+neigh%zbg(i)*thermo%eta(9)

                 neigh%xbg(i) = xa
                 neigh%ybg(i) = ya
                 neigh%zbg(i) = za
              End Do

              Do irgd=1,ntrgd
                 x = rgdxin(irgd)
                 y = rgdyin(irgd)
                 z = rgdzin(irgd)

                 xa = rgdxin(irgd)*thermo%eta(1)+rgdyin(irgd)*thermo%eta(2)+rgdzin(irgd)*thermo%eta(3)
                 ya = rgdxin(irgd)*thermo%eta(4)+rgdyin(irgd)*thermo%eta(5)+rgdzin(irgd)*thermo%eta(6)
                 za = rgdxin(irgd)*thermo%eta(7)+rgdyin(irgd)*thermo%eta(8)+rgdzin(irgd)*thermo%eta(9)

                 rgdxin(irgd) = xa
                 rgdyin(irgd) = ya
                 rgdzin(irgd) = za

                 lrgd=listrgd(-1,irgd)
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd)

                    If (i <= natms) Then
                       neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                       neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                       neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           End If

        Else If (thermo%ensemble == ENS_NPT_NOSE_HOOVER .or. thermo%ensemble == ENS_NPT_NOSE_HOOVER_ANISO) Then

! hoover npt/nst

           Call getcom(neigh%xbg,neigh%ybg,neigh%zbg,com,comm)

           If (thermo%ensemble == ENS_NPT_NOSE_HOOVER) Then

              scale = Exp(tstep*thermo%eta(1))

              Do j=1,nfree
                 i=lstfre(j)

                 neigh%xbg(i) = scale*(neigh%xbg(i)-com(1))+com(1)
                 neigh%ybg(i) = scale*(neigh%ybg(i)-com(2))+com(2)
                 neigh%zbg(i) = scale*(neigh%zbg(i)-com(3))+com(3)
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
                       neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                       neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                       neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           Else

! second order taylor expansion of Exp(tstep*thermo%eta)

              a1 = tstep*thermo%eta(1)
              a2 = tstep*thermo%eta(2)
              a3 = tstep*thermo%eta(3)
              a5 = tstep*thermo%eta(5)
              a6 = tstep*thermo%eta(6)
              a9 = tstep*thermo%eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do j=1,nfree
                 i=lstfre(j)

                 xa = neigh%xbg(i)-com(1)
                 ya = neigh%ybg(i)-com(2)
                 za = neigh%zbg(i)-com(3)

                 neigh%xbg(i) = xa*b1 + ya*b2 + za*b3 + com(1)
                 neigh%ybg(i) = xa*b2 + ya*b5 + za*b6 + com(2)
                 neigh%zbg(i) = xa*b3 + ya*b6 + za*b9 + com(3)
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
                       neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                       neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                       neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           End If

        Else If (thermo%ensemble == ENS_NPT_LANGEVIN .or. &
                 thermo%ensemble == ENS_NPT_LANGEVIN_ANISO .or. &
                 thermo%ensemble == ENS_NPT_MTK .or. &
                 thermo%ensemble == ENS_NPT_MTK_ANISO) Then

! Langevin and MTK npt/nst

           If (thermo%ensemble == ENS_NPT_LANGEVIN .or. thermo%ensemble == ENS_NPT_MTK) Then

              scale = Exp(tstep*thermo%eta(1))

              Do j=1,nfree
                 i=lstfre(j)

                 neigh%xbg(i) = scale*neigh%xbg(i)
                 neigh%ybg(i) = scale*neigh%ybg(i)
                 neigh%zbg(i) = scale*neigh%zbg(i)
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
                       neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                       neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                       neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           Else

! second order taylor expansion of Exp(tstep*thermo%eta)

              a1 = tstep*thermo%eta(1)
              a2 = tstep*thermo%eta(2)
              a3 = tstep*thermo%eta(3)
              a5 = tstep*thermo%eta(5)
              a6 = tstep*thermo%eta(6)
              a9 = tstep*thermo%eta(9)

              b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
              b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
              b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
              b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
              b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
              b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

              Do j=1,nfree
                 i=lstfre(j)

                 xa = neigh%xbg(i)
                 ya = neigh%ybg(i)
                 za = neigh%zbg(i)

                 neigh%xbg(i) = xa*b1 + ya*b2 + za*b3
                 neigh%ybg(i) = xa*b2 + ya*b5 + za*b6
                 neigh%zbg(i) = xa*b3 + ya*b6 + za*b9
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
                       neigh%xbg(i) = neigh%xbg(i) - x + rgdxin(irgd)
                       neigh%ybg(i) = neigh%ybg(i) - y + rgdyin(irgd)
                       neigh%zbg(i) = neigh%zbg(i) - z + rgdzin(irgd)
                    End If
                 End Do
              End Do

           End If

        End If

! Halo final RB members positions across onto neighbouring domains

        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,neigh%xbg,neigh%ybg,neigh%zbg,comm)

     End If

     Deallocate (rgdxin,rgdyin,rgdzin, Stat = fail)
     If (fail > 0) Then
        Write(message,'(a)') 'xscale deallocation failure, node'
        Call error(0,message)
     End If

  End If

  Call pbcshift(imcon,cell,natms,stats%xin,stats%yin,stats%zin)
  If (neigh%unconditional_update) Call pbcshift(imcon,cell,natms,neigh%xbg,neigh%ybg,neigh%zbg)

End Subroutine xscale


!!!!!!!!!!!!!!!!!!!!! THIS IS QUATERNIONS_CONTAINER !!!!!!!!!!!!!!!!!!!!
!
! Subroutine q_setup - sets quaternions for RB dynamics
!
! Subroutine getrotmat - constructs rotation matrix from quaternions
!
! Subroutine no_squish - implements the symplectic no_squish quaternion
!                        algorithm
!
! Subroutine q_update - update the quaternions in LFV scope
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine q_setup(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting up RBs' quaternions
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( comms_type ), Intent( InOut ) :: comm
  Integer           :: fail,irgd,jrgd,krgd,lrgd,rgdtyp, &
                       ill,i1,i2,i3,itmp
  Real( Kind = wp ) :: rot(1:9),aa(1:9),rsq,tol, &
                       aq,bq,cq,dq,eq,fq,gq,hq,rnorm, x,y,z

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  Character( Len = 256 ) :: message

  fail = 0
  Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'q_setup allocation failure'
     Call error(0,message)
  End If

! quaternions for all RB on this domain

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        i1=indrgd(rgdind(1,rgdtyp),irgd)
        i2=indrgd(rgdind(2,rgdtyp),irgd)
        i3=indrgd(rgdind(3,rgdtyp),irgd)

! group basis vectors

        aa(1)=xxx(i1)-xxx(i2)
        aa(4)=yyy(i1)-yyy(i2)
        aa(7)=zzz(i1)-zzz(i2)

! minimum image convention for bond vectors

        Call images(imcon,cell,1,aa(1),aa(4),aa(7))

        ill=rgdind(0,rgdtyp)
        If (ill == 0) Then
           aa(2)=xxx(i1)-xxx(i3)
           aa(5)=yyy(i1)-yyy(i3)
           aa(8)=zzz(i1)-zzz(i3)
        Else
           rsq=Sqrt(aa(1)**2+aa(4)**2+aa(7)**2)
           If      (Abs(aa(7)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(4)**2+aa(7)**2)
              aa(2)= 0.0_wp
              aa(5)= aa(7)/rsq
              aa(8)=-aa(4)/rsq
           Else If (Abs(aa(4)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(4)**2+aa(1)**2)
              aa(2)=-aa(4)/rsq
              aa(5)= aa(1)/rsq
              aa(8)= 0.0_wp
           Else If (Abs(aa(1)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(1)**2+aa(7)**2)
              aa(2)=-aa(7)/rsq
              aa(5)= 0.0_wp
              aa(8)= aa(1)/rsq
           End If
        End If

! minimum image convention for bond vectors

        Call images(imcon,cell,1,aa(2),aa(5),aa(8))

        aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
        aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
        aa(9)=aa(1)*aa(5)-aa(4)*aa(2)

! group rotational matrix

        rot(1)=rgdaxs(1,rgdtyp)*aa(1)+rgdaxs(4,rgdtyp)*aa(2)+rgdaxs(7,rgdtyp)*aa(3)
        rot(2)=rgdaxs(2,rgdtyp)*aa(1)+rgdaxs(5,rgdtyp)*aa(2)+rgdaxs(8,rgdtyp)*aa(3)
        rot(3)=rgdaxs(3,rgdtyp)*aa(1)+rgdaxs(6,rgdtyp)*aa(2)+rgdaxs(9,rgdtyp)*aa(3)
        rot(4)=rgdaxs(1,rgdtyp)*aa(4)+rgdaxs(4,rgdtyp)*aa(5)+rgdaxs(7,rgdtyp)*aa(6)
        rot(5)=rgdaxs(2,rgdtyp)*aa(4)+rgdaxs(5,rgdtyp)*aa(5)+rgdaxs(8,rgdtyp)*aa(6)
        rot(6)=rgdaxs(3,rgdtyp)*aa(4)+rgdaxs(6,rgdtyp)*aa(5)+rgdaxs(9,rgdtyp)*aa(6)
        rot(7)=rgdaxs(1,rgdtyp)*aa(7)+rgdaxs(4,rgdtyp)*aa(8)+rgdaxs(7,rgdtyp)*aa(9)
        rot(8)=rgdaxs(2,rgdtyp)*aa(7)+rgdaxs(5,rgdtyp)*aa(8)+rgdaxs(8,rgdtyp)*aa(9)
        rot(9)=rgdaxs(3,rgdtyp)*aa(7)+rgdaxs(6,rgdtyp)*aa(8)+rgdaxs(9,rgdtyp)*aa(9)

! determine quaternions from rotational matrix

        aq=rot(1)+rot(5)
        bq=rot(2)-rot(4)
        cq=rot(6)-rot(8)
        dq=rot(2)+rot(4)
        eq=rot(3)+rot(7)
        fq=rot(6)+rot(8)
        gq=rot(3)-rot(7)
        hq=rot(1)-rot(5)

        q0(irgd)=0.5_wp*Sqrt(aq+Sqrt(aq*aq+bq*bq))

        If (q0(irgd) > 1.0e-4_wp) Then
           q1(irgd)=-0.25_wp*cq/q0(irgd)
           q2(irgd)= 0.25_wp*gq/q0(irgd)
           q3(irgd)=-0.25_wp*bq/q0(irgd)
        Else
           q1(irgd)=0.5_wp*Sqrt(hq+Sqrt(hq*hq+dq*dq))

           If (q1(irgd) > 1.0e-4_wp) Then
              q2(irgd)=0.25_wp*dq/q1(irgd)
              q3(irgd)=0.25_wp*eq/q1(irgd)
           Else
              q2(irgd)=0.5_wp*Sqrt(-hq+Sqrt(hq*hq+dq*dq))

              If (q2(irgd) > 1.0e-4_wp) Then
                 q3(irgd)=0.25_wp*fq/q2(irgd)
              Else
                 q3(irgd)=1.0_wp
              End If
           End If
        End If

! normalise quaternions

        rnorm=1.0_wp/Sqrt(q0(irgd)**2+q1(irgd)**2+q2(irgd)**2+q3(irgd)**2)
        q0(irgd)=rnorm*q0(irgd)
        q1(irgd)=rnorm*q1(irgd)
        q2(irgd)=rnorm*q2(irgd)
        q3(irgd)=rnorm*q3(irgd)
     Else
        q0(irgd)=0.0_wp
        q1(irgd)=0.0_wp
        q2(irgd)=0.0_wp
        q3(irgd)=1.0_wp
     End If
  End Do

! Check of quaternion set up with atomic positions

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

! new rotational matrix

        Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

        Do jrgd=1,lrgd
           krgd=krgd+1

           x=rot(1)*rgdx(jrgd,rgdtyp)+rot(2)*rgdy(jrgd,rgdtyp)+rot(3)*rgdz(jrgd,rgdtyp)+rgdxxx(irgd)
           y=rot(4)*rgdx(jrgd,rgdtyp)+rot(5)*rgdy(jrgd,rgdtyp)+rot(6)*rgdz(jrgd,rgdtyp)+rgdyyy(irgd)
           z=rot(7)*rgdx(jrgd,rgdtyp)+rot(8)*rgdy(jrgd,rgdtyp)+rot(9)*rgdz(jrgd,rgdtyp)+rgdzzz(irgd)

           gxx(krgd)=xxx(indrgd(jrgd,irgd))-x
           gyy(krgd)=yyy(indrgd(jrgd,irgd))-y
           gzz(krgd)=zzz(indrgd(jrgd,irgd))-z
        End Do

     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! test quaternion setup

  tol=1.0e-2_wp
  rsq=0.0_wp
  krgd=0
  ill=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then ! Keep consistent with above

        itmp=0
        Do jrgd=1,lrgd
           krgd=krgd+1

           rnorm=gxx(krgd)**2+gyy(krgd)**2+gzz(krgd)**2

           If (rnorm > tol) Then
              itmp=1
              Write(message,'(a,2i7,0p,2f7.3)')                                                            &
                   '*** warning - q_setup failure for RB local_id member, with norm > tolerance : ', &
                   irgd,jrgd,rnorm,tol
              Call info(message)
           End If

           rsq=Max(rsq,rnorm)
        End Do
        If (itmp == 1) ill=ill+1

     End If
  End Do
  Call gmax(comm,rsq)
  If (rsq > tol) Then
     Call gsum(comm,ill)
     Write(message,'(a,i7,2f7.3)') &
        '*** warning - q_setup failure for RBs total, with Max(norm) > tolerance : ', ill,rsq,tol
     Call error(648,message,.true.)
  End If

  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'q_setup deallocation failure'
     Call error(0,message)
  End If

End Subroutine q_setup

Subroutine getrotmat(q0,q1,q2,q3,rot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to construct rotation matrix from quaternions using
! x convention for euler angles
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real( Kind = wp ), Intent( In    ) :: q0,q1,q2,q3
  Real( Kind = wp ), Intent(   Out ) :: rot(1:9)

  rot(1)=q0**2+q1**2-q2**2-q3**2
  rot(2)=2.0_wp*(q1*q2-q0*q3)
  rot(3)=2.0_wp*(q1*q3+q0*q2)
  rot(4)=2.0_wp*(q1*q2+q0*q3)
  rot(5)=q0**2-q1**2+q2**2-q3**2
  rot(6)=2.0_wp*(q2*q3-q0*q1)
  rot(7)=2.0_wp*(q1*q3-q0*q2)
  rot(8)=2.0_wp*(q2*q3+q0*q1)
  rot(9)=q0**2-q1**2-q2**2+q3**2

End Subroutine getrotmat

Subroutine no_squish                    &
           (tstep,rgdrix,rgdriy,rgdriz, &
           q0,q1,q2,q3,p0,p1,p2,p3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine routine to implement the symplectic no_squish
! quaternion algorithm of miller et al j.chem.phys 116 (2002) 8649
!
! copyright - daresbury laboratory
! author    - m.leslie january 2004
! adapted   - i.t.todorov september 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real( Kind = wp ), Intent( In    ) :: tstep,rgdrix,rgdriy,rgdriz
  Real( Kind = wp ), Intent( InOut ) :: q0,q1,q2,q3,p0,p1,p2,p3

  Integer, Parameter :: mrot = 10 ! Rotational timesteps
  Integer            :: m
  Real( Kind = wp )  :: rotstep,zetax,zetay,zetaz,cs,sn, &
                        qn1(0:3),pq2(0:3), &
                        qn2(0:3),pq3(0:3), &
                        qn3(0:3),pq4(0:3), rnorm

! rotation: iterate over mrot rotational time steps

  rotstep=tstep/Real(mrot,wp) ! rotational time step
  Do m=1,mrot
     zetaz=0.125_wp*rgdriz*rotstep* &
           ( -p0*q3+p1*q2-          &
              p2*q1+p3*q0 )
     cs=Cos(zetaz)
     sn=Sin(zetaz)

     qn1(0)=cs*q0-sn*q3
     qn1(1)=cs*q1+sn*q2
     qn1(2)=cs*q2-sn*q1
     qn1(3)=cs*q3+sn*q0

     pq2(0)=cs*p0-sn*p3
     pq2(1)=cs*p1+sn*p2
     pq2(2)=cs*p2-sn*p1
     pq2(3)=cs*p3+sn*p0

     zetay=0.125_wp*rgdriy*rotstep*        &
           ( -pq2(0)*qn1(2)-pq2(1)*qn1(3)+ &
              pq2(2)*qn1(0)+pq2(3)*qn1(1) )
     cs=Cos(zetay)
     sn=Sin(zetay)

     qn2(0)=cs*qn1(0)-sn*qn1(2)
     qn2(1)=cs*qn1(1)-sn*qn1(3)
     qn2(2)=cs*qn1(2)+sn*qn1(0)
     qn2(3)=cs*qn1(3)+sn*qn1(1)

     pq3(0)=cs*pq2(0)-sn*pq2(2)
     pq3(1)=cs*pq2(1)-sn*pq2(3)
     pq3(2)=cs*pq2(2)+sn*pq2(0)
     pq3(3)=cs*pq2(3)+sn*pq2(1)

     zetax=0.250_wp*rgdrix*rotstep*        &
           ( -pq3(0)*qn2(1)+pq3(1)*qn2(0)+ &
              pq3(2)*qn2(3)-pq3(3)*qn2(2) )
     cs=Cos(zetax)
     sn=Sin(zetax)

     qn3(0)=cs*qn2(0)-sn*qn2(1)
     qn3(1)=cs*qn2(1)+sn*qn2(0)
     qn3(2)=cs*qn2(2)+sn*qn2(3)
     qn3(3)=cs*qn2(3)-sn*qn2(2)

     pq4(0)=cs*pq3(0)-sn*pq3(1)
     pq4(1)=cs*pq3(1)+sn*pq3(0)
     pq4(2)=cs*pq3(2)+sn*pq3(3)
     pq4(3)=cs*pq3(3)-sn*pq3(2)

     zetay=0.125_wp*rgdriy*rotstep*        &
           ( -pq4(0)*qn3(2)-pq4(1)*qn3(3)+ &
              pq4(2)*qn3(0)+pq4(3)*qn3(1) )
     cs=Cos(zetay)
     sn=Sin(zetay)

     qn2(0)=cs*qn3(0)-sn*qn3(2)
     qn2(1)=cs*qn3(1)-sn*qn3(3)
     qn2(2)=cs*qn3(2)+sn*qn3(0)
     qn2(3)=cs*qn3(3)+sn*qn3(1)

     pq3(0)=cs*pq4(0)-sn*pq4(2)
     pq3(1)=cs*pq4(1)-sn*pq4(3)
     pq3(2)=cs*pq4(2)+sn*pq4(0)
     pq3(3)=cs*pq4(3)+sn*pq4(1)

     zetaz=0.125_wp*rgdriz*rotstep*        &
           ( -pq3(0)*qn2(3)+pq3(1)*qn2(2)- &
              pq3(2)*qn2(1)+pq3(3)*qn2(0) )
     cs=Cos(zetaz)
     sn=Sin(zetaz)

     q0=cs*qn2(0)-sn*qn2(3)
     q1=cs*qn2(1)+sn*qn2(2)
     q2=cs*qn2(2)-sn*qn2(1)
     q3=cs*qn2(3)+sn*qn2(0)

     p0=cs*pq3(0)-sn*pq3(3)
     p1=cs*pq3(1)+sn*pq3(2)
     p2=cs*pq3(2)-sn*pq3(1)
     p3=cs*pq3(3)+sn*pq3(0)
  End Do

! stay normalised (long term drift)

  rnorm=1.0_wp/Sqrt(q0**2+q1**2+q2**2+q3**2)

  q0=q0*rnorm
  q1=q1*rnorm
  q2=q2*rnorm
  q3=q3*rnorm

End Subroutine no_squish

subroutine q_update                                       &
           (tstep,oxp,oyp,ozp,oxq,oyq,ozq,mxquat,quattol, &
           q0,q1,q2,q3,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to update the quaternions in LFV scope
!
! copyright - daresbury laboratory
! author    - t.forester october 1993
! adapted   - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( InOut ) :: safe
  Integer,           Intent( In    ) :: mxquat
  Real( Kind = wp ), Intent( In    ) :: tstep,oxp,oyp,ozp,oxq,oyq,ozq,quattol
  Real( Kind = wp ), Intent( InOut ) :: q0,q1,q2,q3

  Integer           :: itq
  Real( Kind = wp ) :: eps,rnorm,qn0,qn1,qn2,qn3, &
                       qn0a,qn1a,qn2a,qn3a,qn0b,qn1b,qn2b,qn3b

! first estimate of new quaternions (laboratory frame)

  qn0=q0+(-q1*oxp - q2*oyp - q3*ozp)*tstep*0.5_wp
  qn1=q1+( q0*oxp - q3*oyp + q2*ozp)*tstep*0.5_wp
  qn2=q2+( q3*oxp + q0*oyp - q1*ozp)*tstep*0.5_wp
  qn3=q3+(-q2*oxp + q1*oyp + q0*ozp)*tstep*0.5_wp

! first iteration of new quaternions (lab fixed)

  qn0b=0.0_wp ; qn1b=0.0_wp ; qn2b=0.0_wp ; qn3b=0.0_wp

  itq=0 ; eps=2.0_wp*quattol
  Do While (itq < mxquat .and. eps > quattol)
     qn0a=0.5_wp*(-q1 *oxp - q2 *oyp - q3 *ozp) + &
          0.5_wp*(-qn1*oxq - qn2*oyq - qn3*ozq)
     qn1a=0.5_wp*( q0 *oxp - q3 *oyp + q2 *ozp) + &
          0.5_wp*( qn0*oxq - qn3*oyq + qn2*ozq)
     qn2a=0.5_wp*( q3 *oxp + q0 *oyp - q1 *ozp) + &
          0.5_wp*( qn3*oxq + qn0*oyq - qn1*ozq)
     qn3a=0.5_wp*(-q2 *oxp + q1 *oyp + q0 *ozp) + &
          0.5_wp*(-qn2*oxq + qn1*oyq + qn0*ozq)

     qn0=q0+0.5_wp*qn0a*tstep
     qn1=q1+0.5_wp*qn1a*tstep
     qn2=q2+0.5_wp*qn2a*tstep
     qn3=q3+0.5_wp*qn3a*tstep

     rnorm=1.0_wp/Sqrt(qn0**2+qn1**2+qn2**2+qn3**2)

     qn0=qn0*rnorm
     qn1=qn1*rnorm
     qn2=qn2*rnorm
     qn3=qn3*rnorm

! convergence test

     eps=Sqrt(((qn0a-qn0b)**2+(qn1a-qn1b)**2+(qn2a-qn2b)**2+(qn3a-qn3b)**2)*tstep**2)

     qn0b=qn0a
     qn1b=qn1a
     qn2b=qn2a
     qn3b=qn3a

     itq=itq+1
  End Do

! store new quaternions

  q0=qn0
  q1=qn1
  q2=qn2
  q3=qn3

! Check safety

  If (itq == mxquat) safe=.false.

End Subroutine q_update

End module rigid_bodies


