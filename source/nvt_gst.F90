Module nvt_gst
  Use kinds, Only : wp, li
  Use comms,       Only : comms_type,gmax
  Use setup,           Only : mxpmf,mxtpmf,zero_plus,boltz
  Use configuration,      Only : imcon,cell,natms,nlast,nfree, &
                                 lstfre,weight,                &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use domains,     Only : map
  Use kinetics,     Only : kinstress,kinstresf,kinstrest,getkin,getknf,getknt,getknr
  Use core_shell,  Only : legshl
  Use constraints, Only : constraints_tags,apply_shake,&
                          constraints_rattle, constraints_type
  Use pmf,         Only : passpmf,pmf_tags,pmf_rattle
  Use rigid_bodies,    Only : lashp_rgd,lishp_rgd,lshmv_rgd,mxatms,mxlrgd, &
                              ntrgd,rgdx,rgdy,rgdz,rgdxxx,rgdyyy,rgdzzz, &
                              rgdoxx,rgdoyy,rgdozz,rgdvxx,rgdvyy,rgdvzz, &
                              q0,q1,q2,q3,indrgd,listrgd,rgdfrz,rgdwgt, &
                              rgdrix,rgdriy,rgdriz,mxrgd,rgdind,getrotmat, &
                              no_squish,rigid_bodies_stress
  Use numerics, Only : images,box_mueller_saru2
  Use errors_warnings, Only : error,info
  Use thermostat, Only : thermostat_type
  Use statistics, Only : stats_type
  Use timer, Only :  timer_type
  Implicit None

  Private

  Public :: nvt_g0_vv, nvt_g1_vv, nvt_g0_scl, nvt_g1_scl

Contains

  Subroutine nvt_g0_vv                          &
             (isw,lvar,mndis,mxdis,mxstp,tstep, &
             nstep,degfre,                 &
             consv,                             &
             strkin,engke,                      &
             megpmf,strpmf,virpmf,cons,stat,thermo,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian equations of motion in
  ! molecular dynamics - velocity verlet with gentle ergodic thermostat -
  ! a Nose-Hoover thermostat chained to a Langevin thermostat
  !
  ! reference1: B. Leimkuhler, E. Noorizadeh, F. Theil
  !             J. Stat. Phys. (2009) 135: 261-277
  !
  ! reference2: A. Samoletov, M.A.J. Chaplain, C.P. Dettmann
  !             J. Stat. Phys. (2007) 128, 1321-1336

  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  ! amended   - i.t.todorov march 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: isw

    Logical,           Intent( In    ) :: lvar
    Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ), Intent( InOut ) :: tstep

    Integer,           Intent( In    ) :: nstep
    Integer(Kind=li),  Intent( In    ) :: degfre

    Real( Kind = wp ), Intent(   Out ) :: consv

    Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke

    Integer,           Intent( In    ) :: megpmf
    Real( Kind = wp ), Intent( InOut ) ::  &
                                          strpmf(1:9),virpmf
    Type( stats_type), Intent( InOut ) :: stat
    Type( constraints_type), Intent( InOut ) :: cons
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true.
    Logical                 :: safe,lcol,lfst,lv_up,lv_dn
    Integer,           Save :: mxkit,kit
    Integer                 :: fail(1:9),i
    Real( Kind = wp ), Save :: qmass,ceng
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: chitdr,cintdr
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

    Integer,           Allocatable :: lstopt(:,:),listot(:)
    Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
    Integer,           Allocatable :: indpmf(:,:,:)
    Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)
    Character( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
       If (cons%megcon > 0) Then
          Allocate (lstopt(0:2,1:cons%mxcons),listot(1:mxatms),          Stat=fail(2))
          Allocate (dxx(1:cons%mxcons),dyy(1:cons%mxcons),dzz(1:cons%mxcons),      Stat=fail(3))
       End If
       If (megpmf > 0) Then
          Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat=fail(4))
          Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),         Stat=fail(5))
       End If
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
    End If
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail(7))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail(8))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_g0 allocation failure'
       Call error(0,message)
    End If


    If (newjob) Then
       newjob = .false.

  ! inertia parameter for Nose-Hoover thermostat

       qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
       ceng  = 2.0_wp*thermo%sigma

  ! set number of constraint+pmf shake iterations

       If (cons%megcon > 0 .or.  megpmf > 0) mxkit=1
       If (cons%megcon > 0 .and. megpmf > 0) mxkit=cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Then
         Call constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot,cons,comm)
       End If

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (megpmf > 0) Then
         Call pmf_tags(lstitr,indpmf,pxx,pyy,pzz,comm)
       End If
    End If

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep
    lv_up = .false.
    lv_dn = .false.

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! store initial values

       Do i=1,natms
          xxt(i) = xxx(i)
          yyt(i) = yyy(i)
          zzt(i) = zzz(i)

          vxt(i) = vxx(i)
          vyt(i) = vyy(i)
          vzt(i) = vzz(i)

          fxt(i) = fxx(i)
          fyt(i) = fyy(i)
          fzt(i) = fzz(i)
       End Do

  ! store initial values of integration variables

       chitdr=thermo%chi_t
       cintdr=thermo%cint

  100  Continue

  ! constraint virial and stress tensor

       If (cons%megcon > 0) Then
          stat%vircon=0.0_wp
          stat%strcon=0.0_wp
       End If

  ! PMF virial and stress tensor

       If (megpmf > 0) Then
          virpmf=0.0_wp
          strpmf=0.0_wp
       End If

  ! integrate and apply nvt_g0_scl thermostat - 1/2 step

       Call nvt_g0_scl &
             (hstep,degfre,isw,nstep,ceng,qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,engke,thermo,comm)

  ! update velocity and position

       Do i=1,natms
          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*fxt(i)
             vyy(i)=vyy(i)+tmp*fyt(i)
             vzz(i)=vzz(i)+tmp*fzt(i)

             xxx(i)=xxt(i)+tstep*vxx(i)
             yyy(i)=yyt(i)+tstep*vyy(i)
             zzz(i)=zzt(i)+tstep*vzz(i)
          End If
       End Do

  ! SHAKE procedures

       If (cons%megcon > 0 .or. megpmf > 0) Then
        Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,dxx,dyy,dzz,pxx,pyy,pzz,&
          listot,lstitr,lstopt,indpmf,&
          megpmf,virpmf,strpmf,stat,cons,tmr,comm)
       End If

  ! check timestep for variable timestep

       If (lvar) Then

  ! update maximum distance a particle has travelled

          mxdr = 0.0_wp
          Do i=1,natms
             If (legshl(0,i) >= 0) &
                mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
          End Do
          mxdr=Sqrt(mxdr)
          Call gmax(comm,mxdr)

          If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

  ! scale tstep and derivatives

             If (mxdr > mxdis) Then
                lv_up = .true.
                If (lv_dn) Then
                   tstep = 0.75_wp*tstep
                   hstep = 0.50_wp*tstep
                Else
                   tstep = hstep
                   hstep = 0.50_wp*tstep
                End If
                Write(message,"( &
                  & 'timestep decreased, new timestep is:',3x,1p,e12.4,/)") tstep
                Call info(message,.true.)
             End If
             If (mxdr < mndis) Then
                lv_dn = .true.
                If (lv_up) Then
                   tstep = 1.50_wp*tstep
                   hstep = 0.50_wp*tstep
                Else
                   hstep = tstep
                   tstep = 2.00_wp*tstep
                End If
                If (tstep > mxstp) Then
                   tstep = mxstp
                   hstep = 0.50_wp*tstep
                End If
                Write(message,"( &
                  & 'timestep increased, new timestep is:',3x,1p,e12.4,/)") tstep
                Call info(message,.true.)
             End If
             rstep = 1.0_wp/tstep

  ! restore initial conditions

             thermo%chi_t=chitdr
             thermo%cint=cintdr

             Do i=1,natms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! second stage of velocity verlet algorithm

    Else

  ! update velocity

       Do i=1,natms
          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*fxx(i)
             vyy(i)=vyy(i)+tmp*fyy(i)
             vzz(i)=vzz(i)+tmp*fzz(i)
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. megpmf > 0) Then
          Do i=1,kit
             lfst = (i == 1)
             lcol = (i == kit)

             If (cons%megcon > 0) Then
               Call constraints_rattle &
                 (tstep,lfst,lcol, &
                 lstopt,dxx,dyy,dzz,listot,      &
                 vxx,vyy,vzz,stat,cons,tmr,comm)
             End If

             If (megpmf > 0) Then
               Call pmf_rattle &
                 (cons%max_iter_shake,cons%tolerance,tstep,lfst,lcol, &
                 indpmf,pxx,pyy,pzz,             &
                 vxx,vyy,vzz,comm)
             End If
          End Do
       End If

  ! integrate and apply nvt_g0_scl thermostat - 1/2 step

       Call nvt_g0_scl &
             (hstep,degfre,isw,nstep,ceng,qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,engke,thermo,comm)

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*qmass*thermo%chi_t**2 + ceng*thermo%cint

  ! kinetic contribution to stress tensor

       Call kinstress(vxx,vyy,vzz,strkin,comm)

    End If

    If (cons%megcon > 0 .or. megpmf > 0) Then
       Deallocate (lstitr,           Stat=fail(1))
       If (cons%megcon > 0) Then
          Deallocate (lstopt,listot, Stat=fail(2))
          Deallocate (dxx,dyy,dzz,   Stat=fail(3))
       End If
       If (megpmf > 0) Then
          Deallocate (indpmf,        Stat=fail(4))
          Deallocate (pxx,pyy,pzz,   Stat=fail(5))
       End If
       Deallocate (oxt,oyt,ozt,      Stat=fail(6))
    End If
    Deallocate (xxt,yyt,zzt,         Stat=fail(7))
    Deallocate (vxt,vyt,vzt,         Stat=fail(8))
    Deallocate (fxt,fyt,fzt,         Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_g0 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nvt_g0_vv

  Subroutine nvt_g1_vv                          &
             (isw,lvar,mndis,mxdis,mxstp,tstep, &
             nstep,degfre,                 &
             consv,                             &
             strkin,strknf,strknt,engke,engrot, &
             megpmf,strpmf,virpmf,              &
             strcom,vircom,cons,stat,thermo,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
  ! RBs, equations of motion in molecular dynamics
  ! - velocity verlet with gentle ergodic thermostat -
  ! a Nose-Hoover thermostat chained to a Langevin thermostat -
  !
  ! reference1: B. Leimkuhler, E. Noorizadeh, F. Theil
  !             J. Stat. Phys. (2009) 135: 261-277
  !
  ! reference2: A. Samoletov, M.A.J. Chaplain, C.P. Dettmann
  !             J. Stat. Phys. (2007) 128, 1321-1336

  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  ! amended   - i.t.todorov march 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: isw

    Logical,           Intent( In    ) :: lvar
    Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ), Intent( InOut ) :: tstep

    Integer,           Intent( In    ) :: nstep
    Integer(Kind=li),  Intent( In    ) :: degfre

    Real( Kind = wp ), Intent(   Out ) :: consv

    Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                          strknf(1:9),strknt(1:9),engrot

    Integer,           Intent( In    ) :: megpmf
    Real( Kind = wp ), Intent( InOut ) ::  &
                                          strpmf(1:9),virpmf

    Real( Kind = wp ), Intent( InOut ) :: strcom(1:9),vircom
    Type( stats_type), Intent( InOut ) :: stat
    Type( constraints_type), Intent( InOut ) :: cons
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true. , &
                               unsafe = .false.
    Logical                 :: safe,lcol,lfst,lv_up,lv_dn
    Integer,           Save :: mxkit,kit
    Integer                 :: fail(1:14),matms,i,j,i1,i2, &
                               irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ), Save :: qmass,ceng
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: chitdr,cintdr
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp
    Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                               opx,opy,opz,fmx,fmy,fmz,       &
                               tqx,tqy,tqz,trx,try,trz,       &
                               qt0,qt1,qt2,qt3,p0,p1,p2,p3,   &
                               vpx,vpy,vpz


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

    Integer,           Allocatable :: lstopt(:,:),listot(:)
    Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
    Integer,           Allocatable :: indpmf(:,:,:)
    Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

    Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
    Real( Kind = wp ), Allocatable :: q0t(:),q1t(:),q2t(:),q3t(:)
    Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
    Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
    Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)
    Character( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       If (cons%megcon > 0) Then
          Allocate (lstopt(0:2,1:cons%mxcons),listot(1:mxatms),          Stat=fail( 2))
          Allocate (dxx(1:cons%mxcons),dyy(1:cons%mxcons),dzz(1:cons%mxcons),      Stat=fail( 3))
       End If
       If (megpmf > 0) Then
          Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat=fail( 4))
          Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),         Stat=fail( 5))
       End If
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail( 6))
    End If
    Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), &
                                                                    Stat=fail( 7))
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
    Allocate (q0t(1:mxrgd),q1t(1:mxrgd),q2t(1:mxrgd),q3t(1:mxrgd),  Stat=fail(11))
    Allocate (rgdxxt(1:mxrgd),rgdyyt(1:mxrgd),rgdzzt(1:mxrgd),      Stat=fail(12))
    Allocate (rgdvxt(1:mxrgd),rgdvyt(1:mxrgd),rgdvzt(1:mxrgd),      Stat=fail(13))
    Allocate (rgdoxt(1:mxrgd),rgdoyt(1:mxrgd),rgdozt(1:mxrgd),      Stat=fail(14))
    If (Any(fail > 0)) Then
       Write(message,'(a,i0)') 'nvt_g1 allocation failure'
       Call error(0,message)
    End If


    If (newjob) Then
       newjob = .false.

  ! inertia parameter for Nose-Hoover thermostat

       qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
       ceng  = 2.0_wp*thermo%sigma

  ! set number of constraint+pmf shake iterations

       If (cons%megcon > 0 .or.  megpmf > 0) mxkit=1
       If (cons%megcon > 0 .and. megpmf > 0) mxkit=cons%max_iter_shake

  ! unsafe positioning due to possibly locally shared RBs

       unsafe=(Any(map == comm%idnode))
    End If

  ! set matms

    matms=nlast
    If (comm%mxnode == 1) matms=natms

    If (cons%megcon > 0 .or. megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Then
         Call constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot,cons,comm)
       End If

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (megpmf > 0) Then
         Call pmf_tags(lstitr,indpmf,pxx,pyy,pzz,comm)
       End If
    End If

  ! Get the RB particles vectors wrt the RB's COM

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

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep
    lv_up = .false.
    lv_dn = .false.

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! store initial values

       Do i=1,matms
          xxt(i) = xxx(i)
          yyt(i) = yyy(i)
          zzt(i) = zzz(i)

          vxt(i) = vxx(i)
          vyt(i) = vyy(i)
          vzt(i) = vzz(i)

          fxt(i) = fxx(i)
          fyt(i) = fyy(i)
          fzt(i) = fzz(i)
       End Do

       Do irgd=1,ntrgd
          q0t(irgd)=q0(irgd)
          q1t(irgd)=q1(irgd)
          q2t(irgd)=q2(irgd)
          q3t(irgd)=q3(irgd)

          rgdxxt(irgd) = rgdxxx(irgd)
          rgdyyt(irgd) = rgdyyy(irgd)
          rgdzzt(irgd) = rgdzzz(irgd)

          rgdvxt(irgd) = rgdvxx(irgd)
          rgdvyt(irgd) = rgdvyy(irgd)
          rgdvzt(irgd) = rgdvzz(irgd)

          rgdoxt(irgd) = rgdoxx(irgd)
          rgdoyt(irgd) = rgdoyy(irgd)
          rgdozt(irgd) = rgdozz(irgd)
       End Do

  ! store initial values of integration variables

       chitdr=thermo%chi_t
       cintdr=thermo%cint

  100  Continue

  ! constraint virial and stress tensor

       If (cons%megcon > 0) Then
          stat%vircon=0.0_wp
          stat%strcon=0.0_wp
       End If

  ! PMF virial and stress tensor

       If (megpmf > 0) Then
          virpmf=0.0_wp
          strpmf=0.0_wp
       End If

  ! integrate and apply nvt_g1_scl thermostat - 1/2 step

       Call nvt_g1_scl &
             (hstep,degfre,isw,nstep,ceng,qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,                                                &
             rgdvxx,rgdvyy,rgdvzz,                                       &
             rgdoxx,rgdoyy,rgdozz,                                       &
             engke,engrot,thermo,comm)

  ! update velocity and position of FPs

       Do j=1,nfree
          i=lstfre(j)

          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*fxt(i)
             vyy(i)=vyy(i)+tmp*fyt(i)
             vzz(i)=vzz(i)+tmp*fzt(i)

             xxx(i)=xxt(i)+tstep*vxx(i)
             yyy(i)=yyt(i)+tstep*vyy(i)
             zzz(i)=zzt(i)+tstep*vzz(i)
          End If
       End Do

  ! SHAKE procedures

       If (cons%megcon > 0 .or. megpmf > 0) Then
        Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,dxx,dyy,dzz,pxx,pyy,pzz,&
          listot,lstitr,lstopt,indpmf,&
          megpmf,virpmf,strpmf,stat,cons,tmr,comm)
       End If

  ! update velocity and position of RBs

       krgd=0
       Do irgd=1,ntrgd
          rgdtyp=listrgd(0,irgd)

  ! For all good RBs

          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) < lrgd) Then

  ! calculate COM force and torque

             fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
             tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=indrgd(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rgdfrz(0,rgdtyp) == 0) Then
                   fmx=fmx+fxt(i)
                   fmy=fmy+fyt(i)
                   fmz=fmz+fzt(i)
                End If

                tqx=tqx+ggy(krgd)*fzt(i)-ggz(krgd)*fyt(i)
                tqy=tqy+ggz(krgd)*fxt(i)-ggx(krgd)*fzt(i)
                tqz=tqz+ggx(krgd)*fyt(i)-ggy(krgd)*fxt(i)
             End Do

  ! If the RB has 2+ frozen particles (ill=1) the net torque
  ! must align along the axis of rotation

             If (rgdfrz(0,rgdtyp) > 1) Then
                i1=indrgd(rgdind(1,rgdtyp),irgd)
                i2=indrgd(rgdind(2,rgdtyp),irgd)

                x(1)=xxt(i1)-xxt(i2)
                y(1)=yyt(i1)-yyt(i2)
                z(1)=zzt(i1)-zzt(i2)

                Call images(imcon,cell,1,x,y,z)

                tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)/(x(1)**2+y(1)**2+z(1)**2)
                tqx=x(1)*tmp
                tqy=y(1)*tmp
                tqz=z(1)*tmp
             End If

  ! current rotation matrix

             Call getrotmat(q0t(irgd),q1t(irgd),q2t(irgd),q3t(irgd),rot)

  ! calculate torque in principal frame

             trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
             try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
             trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

  ! calculate quaternion torques

             qt0=2.0_wp*(-q1t(irgd)*trx-q2t(irgd)*try-q3t(irgd)*trz)
             qt1=2.0_wp*( q0t(irgd)*trx-q3t(irgd)*try+q2t(irgd)*trz)
             qt2=2.0_wp*( q3t(irgd)*trx+q0t(irgd)*try-q1t(irgd)*trz)
             qt3=2.0_wp*(-q2t(irgd)*trx+q1t(irgd)*try+q0t(irgd)*trz)

  ! recover quaternion momenta

             opx=rgdoxx(irgd)*rgdrix(1,rgdtyp)
             opy=rgdoyy(irgd)*rgdriy(1,rgdtyp)
             opz=rgdozz(irgd)*rgdriz(1,rgdtyp)

             p0=2.0_wp*(-q1t(irgd)*opx-q2t(irgd)*opy-q3t(irgd)*opz)
             p1=2.0_wp*( q0t(irgd)*opx-q3t(irgd)*opy+q2t(irgd)*opz)
             p2=2.0_wp*( q3t(irgd)*opx+q0t(irgd)*opy-q1t(irgd)*opz)
             p3=2.0_wp*(-q2t(irgd)*opx+q1t(irgd)*opy+q0t(irgd)*opz)

  ! update quaternion momenta to half step

             p0=p0+hstep*qt0
             p1=p1+hstep*qt1
             p2=p2+hstep*qt2
             p3=p3+hstep*qt3

  ! rotate RB quaternions - update q to full timestep & amend p
  ! and get new rotation matrix

             Call no_squish                                             &
             (tstep,rgdrix(2,rgdtyp),rgdriy(2,rgdtyp),rgdriz(2,rgdtyp), &
             q0(irgd),q1(irgd),q2(irgd),q3(irgd),p0,p1,p2,p3)
             Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

  ! update RB angular & COM velocities to half step

             opx=0.5_wp*(-q1(irgd)*p0+q0(irgd)*p1+q3(irgd)*p2-q2(irgd)*p3)
             opy=0.5_wp*(-q2(irgd)*p0-q3(irgd)*p1+q0(irgd)*p2+q1(irgd)*p3)
             opz=0.5_wp*(-q3(irgd)*p0+q2(irgd)*p1-q1(irgd)*p2+q0(irgd)*p3)

             rgdoxx(irgd)=opx*rgdrix(2,rgdtyp)
             rgdoyy(irgd)=opy*rgdriy(2,rgdtyp)
             rgdozz(irgd)=opz*rgdriz(2,rgdtyp)

             tmp=hstep/rgdwgt(0,rgdtyp)
             rgdvxx(irgd)=rgdvxx(irgd)+tmp*fmx
             rgdvyy(irgd)=rgdvyy(irgd)+tmp*fmy
             rgdvzz(irgd)=rgdvzz(irgd)+tmp*fmz

  ! update RB COM to full step

             rgdxxx(irgd)=rgdxxt(irgd)+tstep*rgdvxx(irgd)
             rgdyyy(irgd)=rgdyyt(irgd)+tstep*rgdvyy(irgd)
             rgdzzz(irgd)=rgdzzt(irgd)+tstep*rgdvzz(irgd)

  ! update RB members positions and halfstep velocities

             Do jrgd=1,lrgd
                i=indrgd(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then
                   If (rgdfrz(jrgd,rgdtyp) == 0) Then
                      x(1)=rgdx(jrgd,rgdtyp)
                      y(1)=rgdy(jrgd,rgdtyp)
                      z(1)=rgdz(jrgd,rgdtyp)

  ! new atomic positions

                      xxx(i)=rot(1)*x(1)+rot(2)*y(1)+rot(3)*z(1) + rgdxxx(irgd)
                      yyy(i)=rot(4)*x(1)+rot(5)*y(1)+rot(6)*z(1) + rgdyyy(irgd)
                      zzz(i)=rot(7)*x(1)+rot(8)*y(1)+rot(9)*z(1) + rgdzzz(irgd)

  ! new atomic velocities in body frame

                      vpx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                      vpy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                      vpz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

  ! DD bound positions

                      If (unsafe) Then
                         x(1)=xxx(i)-xxt(i)
                         y(1)=yyy(i)-yyt(i)
                         z(1)=zzz(i)-zzt(i)
                         Call images(imcon,cell,1,x,y,z)
                         xxx(i)=x(1)+xxt(i)
                         yyy(i)=y(1)+yyt(i)
                         zzz(i)=z(1)+zzt(i)
                      End If

  ! new atomic velocities in lab frame

                      vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rgdvxx(irgd)
                      vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rgdvyy(irgd)
                      vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rgdvzz(irgd)
                   End If
                End If
             End Do

          End If
       End Do

  ! check timestep for variable timestep

       If (lvar) Then

  ! update maximum distance a particle has travelled

          mxdr = 0.0_wp
          Do i=1,natms
             If (legshl(0,i) >= 0) &
                mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
          End Do
          mxdr=Sqrt(mxdr)
          Call gmax(comm,mxdr)

          If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

  ! scale tstep and derivatives

             If (mxdr > mxdis) Then
                lv_up = .true.
                If (lv_dn) Then
                   tstep = 0.75_wp*tstep
                   hstep = 0.50_wp*tstep
                Else
                   tstep = hstep
                   hstep = 0.50_wp*tstep
                End If
                Write(message,"( &
                  & 'timestep decreased, new timestep is:',3x,1p,e12.4,/)") tstep
                Call info(message,.true.)
             End If
             If (mxdr < mndis) Then
                lv_dn = .true.
                If (lv_up) Then
                   tstep = 1.50_wp*tstep
                   hstep = 0.50_wp*tstep
                Else
                   hstep = tstep
                   tstep = 2.00_wp*tstep
                End If
                If (tstep > mxstp) Then
                   tstep = mxstp
                   hstep = 0.50_wp*tstep
                End If
                Write(message,"( &
                  & 'timestep increased, new timestep is:',3x,1p,e12.4,/)") tstep
                Call info(message,.true.)
             End If
             rstep = 1.0_wp/tstep

  ! restore initial conditions

             thermo%chi_t=chitdr
             thermo%cint=cintdr

             Do i=1,matms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do

             Do irgd=1,ntrgd
                q0(irgd)=q0t(irgd)
                q1(irgd)=q1t(irgd)
                q2(irgd)=q2t(irgd)
                q3(irgd)=q3t(irgd)

                rgdvxx(irgd) = rgdvxt(irgd)
                rgdvyy(irgd) = rgdvyt(irgd)
                rgdvzz(irgd) = rgdvzt(irgd)

                rgdoxx(irgd) = rgdoxt(irgd)
                rgdoyy(irgd) = rgdoyt(irgd)
                rgdozz(irgd) = rgdozt(irgd)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! second stage of velocity verlet algorithm

    Else

  ! update velocity of FPs

       Do j=1,nfree
          i=lstfre(j)

          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*fxx(i)
             vyy(i)=vyy(i)+tmp*fyy(i)
             vzz(i)=vzz(i)+tmp*fzz(i)
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. megpmf > 0) Then
          Do i=1,kit
             lfst = (i == 1)
             lcol = (i == kit)

             If (cons%megcon > 0) Then
               Call constraints_rattle &
                 (tstep,lfst,lcol, &
                 lstopt,dxx,dyy,dzz,listot,      &
                 vxx,vyy,vzz,stat,cons,tmr,comm)
             End If

             If (megpmf > 0) Then
               Call pmf_rattle &
                 (cons%max_iter_shake,cons%tolerance,tstep,lfst,lcol, &
                 indpmf,pxx,pyy,pzz,             &
                 vxx,vyy,vzz,comm)
             End If
          End Do
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stress(strcom,ggx,ggy,ggz,comm)
       vircom=-(strcom(1)+strcom(5)+strcom(9))

  ! update velocity of RBs

       krgd=0
       Do irgd=1,ntrgd
          rgdtyp=listrgd(0,irgd)

  ! For all good RBs

          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

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

  ! calculate quaternion torques

             qt0=2.0_wp*(-q1(irgd)*trx-q2(irgd)*try-q3(irgd)*trz)
             qt1=2.0_wp*( q0(irgd)*trx-q3(irgd)*try+q2(irgd)*trz)
             qt2=2.0_wp*( q3(irgd)*trx+q0(irgd)*try-q1(irgd)*trz)
             qt3=2.0_wp*(-q2(irgd)*trx+q1(irgd)*try+q0(irgd)*trz)

  ! recover quaternion momenta at half time step

             opx=rgdoxx(irgd)*rgdrix(1,rgdtyp)
             opy=rgdoyy(irgd)*rgdriy(1,rgdtyp)
             opz=rgdozz(irgd)*rgdriz(1,rgdtyp)

             p0=2.0_wp*(-q1(irgd)*opx-q2(irgd)*opy-q3(irgd)*opz)
             p1=2.0_wp*( q0(irgd)*opx-q3(irgd)*opy+q2(irgd)*opz)
             p2=2.0_wp*( q3(irgd)*opx+q0(irgd)*opy-q1(irgd)*opz)
             p3=2.0_wp*(-q2(irgd)*opx+q1(irgd)*opy+q0(irgd)*opz)

  ! update quaternion momenta to full step

             p0=p0+hstep*qt0
             p1=p1+hstep*qt1
             p2=p2+hstep*qt2
             p3=p3+hstep*qt3

  ! update RB angular & COM velocities to full step

             opx=0.5_wp*(-q1(irgd)*p0+q0(irgd)*p1+q3(irgd)*p2-q2(irgd)*p3)
             opy=0.5_wp*(-q2(irgd)*p0-q3(irgd)*p1+q0(irgd)*p2+q1(irgd)*p3)
             opz=0.5_wp*(-q3(irgd)*p0+q2(irgd)*p1-q1(irgd)*p2+q0(irgd)*p3)

             rgdoxx(irgd)=opx*rgdrix(2,rgdtyp)
             rgdoyy(irgd)=opy*rgdriy(2,rgdtyp)
             rgdozz(irgd)=opz*rgdriz(2,rgdtyp)

             tmp=hstep/rgdwgt(0,rgdtyp)
             rgdvxx(irgd)=rgdvxx(irgd)+tmp*fmx
             rgdvyy(irgd)=rgdvyy(irgd)+tmp*fmy
             rgdvzz(irgd)=rgdvzz(irgd)+tmp*fmz

  ! update RB members velocities

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

  ! integrate and apply nvt_g1_scl thermostat - 1/2 step

       Call nvt_g1_scl &
             (hstep,degfre,isw,nstep,ceng,qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,                                                &
             rgdvxx,rgdvyy,rgdvzz,                                       &
             rgdoxx,rgdoyy,rgdozz,                                       &
             engke,engrot,thermo,comm)

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*qmass*thermo%chi_t**2 + ceng*thermo%cint

  ! update kinetic energy and stress

       Call kinstresf(vxx,vyy,vzz,strknf,comm)
       Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

       strkin=strknf+strknt

    End If

    If (cons%megcon > 0 .or. megpmf > 0) Then
       Deallocate (lstitr,            Stat=fail( 1))
       If (cons%megcon > 0) Then
          Deallocate (lstopt,listot,  Stat=fail( 2))
          Deallocate (dxx,dyy,dzz,    Stat=fail( 3))
       End If
       If (megpmf > 0) Then
          Deallocate (indpmf,         Stat=fail( 4))
          Deallocate (pxx,pyy,pzz,    Stat=fail( 5))
       End If
       Deallocate (oxt,oyt,ozt,       Stat=fail( 6))
    End If
    Deallocate (ggx,ggy,ggz,          Stat=fail( 7))
    Deallocate (xxt,yyt,zzt,          Stat=fail( 8))
    Deallocate (vxt,vyt,vzt,          Stat=fail( 9))
    Deallocate (fxt,fyt,fzt,          Stat=fail(10))
    Deallocate (q0t,q1t,q2t,q3t,      Stat=fail(11))
    Deallocate (rgdxxt,rgdyyt,rgdzzt, Stat=fail(12))
    Deallocate (rgdvxt,rgdvyt,rgdvzt, Stat=fail(13))
    Deallocate (rgdoxt,rgdoyt,rgdozt, Stat=fail(14))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_g1 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nvt_g1_vv

  Subroutine nvt_g0_scl &
             (tstep,degfre,isw,nstep,ceng,qmass,pmass,chip, &
             vxx,vyy,vzz,engke,thermo,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to integrate and apply NVT Nose-Hoover thermostat
  ! with a Langevin process
  !
  ! Note: coupling to NPT barostat included as factor=pmass*chip^2
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov january 2012
  ! amended   - i.t.todorov march 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),                        Intent( In    ) :: tstep,ceng,qmass, &
                                                                 pmass,chip
    Integer(Kind=li),                         Intent( In    ) :: degfre
    Integer,                                  Intent( In    ) :: isw,nstep
    Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
    Real( Kind = wp ),                        Intent(   Out ) :: engke
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: i
    Real( Kind = wp ) :: hstep,qstep,factor,scale,fex,r_0


  ! timestep derivative and factor

    hstep  = 0.5_wp*tstep
    qstep  = 0.5_wp*hstep
    factor = pmass*chip**2

  ! update thermo%chi(=thermo%cint) to 1/4*tstep

    thermo%cint=thermo%cint + qstep*thermo%chi_t

  ! calculate kinetic energy

    engke=getkin(vxx,vyy,vzz,comm)

    fex=Exp(-thermo%gama*hstep)

  ! generate a Gaussian random number for use in the
  ! Langevin process on the thermostat friction

    Call box_mueller_saru2(Int(degfre/3_li),nstep-1,2*isw+1,r_0,.true.)

  ! update thermo%chi_t to 1/2*tstep

    thermo%chi_t=fex*thermo%chi_t + Sqrt((1.0_wp-fex**2) * boltz*thermo%temp/qmass)*r_0 + &
         hstep*(2.0_wp*engke + factor - ceng)/qmass

  ! update thermo%chi(=thermo%cint) to 3/4*tstep

    thermo%cint=thermo%cint + hstep*thermo%chi_t

  ! thermostat the velocities to 1*tstep

    scale=Exp(-tstep*thermo%chi_t)
    Do i=1,natms
       vxx(i)=scale*vxx(i)
       vyy(i)=scale*vyy(i)
       vzz(i)=scale*vzz(i)
    End Do

  ! thermostat the energy consequently

    engke=engke*scale**2

  ! generate a Gaussian random number for use in the
  ! Langevin process on the thermostat friction

    Call box_mueller_saru2(Int(degfre/3_li),nstep-1,2*isw+2,r_0,.true.)

  ! update thermo%chi_t to full (2/2)*tstep

    thermo%chi_t=fex*thermo%chi_t + Sqrt((1.0_wp-fex**2) * boltz*thermo%temp/qmass)*r_0 + &
         hstep*(2.0_wp*engke + factor - ceng)/qmass

  ! update thermo%chi(=thermo%cint) to 4/4*tstep

    thermo%cint=thermo%cint + qstep*thermo%chi_t

  End Subroutine nvt_g0_scl

  Subroutine nvt_g1_scl &
             (tstep,degfre,isw,nstep,ceng,qmass,pmass,chip, &
             vxx,vyy,vzz,                                             &
             rgdvxx,rgdvyy,rgdvzz,                                    &
             rgdoxx,rgdoyy,rgdozz,                                    &
             engke,engrot,thermo,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to integrate and apply NVT Nose-Hoover thermostat
  ! with a Langevin process when singled RBs are present
  !
  ! Note: coupling to NPT barostat included as factor=pmass*chip^2
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov january 2012
  ! amended   - i.t.todorov march 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),                        Intent( In    ) :: tstep,ceng,qmass, &
                                                                 pmass,chip
    Integer(Kind=li),                         Intent( In    ) :: degfre
    Integer,                                  Intent( In    ) :: isw,nstep
    Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: vxx,vyy,vzz
    Real( Kind = wp ), Dimension( 1:mxrgd ),  Intent( InOut ) :: rgdvxx,rgdvyy,rgdvzz
    Real( Kind = wp ), Dimension( 1:mxrgd ),  Intent( InOut ) :: rgdoxx,rgdoyy,rgdozz
    Real( Kind = wp ),                        Intent(   Out ) :: engke,engrot
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( comms_type ),                       Intent( InOut ) :: comm

    Integer           :: i,j,irgd
    Real( Kind = wp ) :: engkf,engkt,hstep,qstep,factor,scale,fex,r_0


  ! timestep derivative and factor

    hstep  = 0.5_wp*tstep
    qstep  = 0.5_wp*hstep
    factor = pmass*chip**2

  ! update thermo%chi(=thermo%cint) to 1/4*tstep

    thermo%cint=thermo%cint + qstep*thermo%chi_t

  ! calculate kinetic energy contributions and rotational energy

    engkf=getknf(vxx,vyy,vzz,comm)
    engkt=getknt(rgdvxx,rgdvyy,rgdvzz,comm)

    engke=engkf+engkt

    engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)

    fex=Exp(-thermo%gama*hstep)

  ! generate a Gaussian random number for use in the
  ! Langevin process on the thermostat friction

    Call box_mueller_saru2(Int(degfre/3_li),nstep-1,2*isw+1,r_0,.true.)

  ! update thermo%chi_t to 1/2*tstep

    thermo%chi_t=fex*thermo%chi_t + Sqrt((1.0_wp-fex**2) * boltz*thermo%temp/qmass)*r_0 + &
         hstep*(2.0_wp*(engke+engrot) + factor - ceng)/qmass

  ! update thermo%chi(=thermo%cint) to 3/4*tstep

    thermo%cint=thermo%cint + hstep*thermo%chi_t

  ! thermostat the velocities to 1*tstep

    scale=Exp(-tstep*thermo%chi_t)
    Do j=1,nfree
       i=lstfre(j)

       vxx(i)=scale*vxx(i)
       vyy(i)=scale*vyy(i)
       vzz(i)=scale*vzz(i)
    End Do

    Do irgd=1,ntrgd
       rgdvxx(irgd)=scale*rgdvxx(irgd)
       rgdvyy(irgd)=scale*rgdvyy(irgd)
       rgdvzz(irgd)=scale*rgdvzz(irgd)

       rgdoxx(irgd)=scale*rgdoxx(irgd)
       rgdoyy(irgd)=scale*rgdoyy(irgd)
       rgdozz(irgd)=scale*rgdozz(irgd)
    End Do

  ! thermostat the energy consequently

    engke=engke*scale**2
    engrot=engrot*scale**2

  ! generate a Gaussian random number for use in the
  ! Langevin process on the thermostat friction

    Call box_mueller_saru2(Int(degfre/3_li),nstep-1,2*isw+2,r_0,.true.)

  ! update thermo%chi_t to full (2/2)*tstep

    thermo%chi_t=fex*thermo%chi_t + Sqrt((1.0_wp-fex**2) * boltz*thermo%temp/qmass)*r_0 + &
         hstep*(2.0_wp*(engke+engrot) + factor - ceng)/qmass

  ! update thermo%chi(=thermo%cint) to 4/4*tstep

    thermo%cint=thermo%cint + qstep*thermo%chi_t

  End Subroutine nvt_g1_scl
End Module nvt_gst
