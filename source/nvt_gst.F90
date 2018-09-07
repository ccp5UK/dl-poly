Module nvt_gst
  Use kinds, Only : wp, li
  Use comms,       Only : comms_type,gmax
  Use setup,           Only : zero_plus,boltz,mxatms
  Use configuration,      Only : imcon,cell,natms,nlast,nfree, &
                                 lstfre,weight,                &
                                 vxx,vyy,vzz
  Use particle,    Only : corePart
  Use domains,     Only : domains_type
  Use kinetics,     Only : kinstress,kinstresf,kinstrest,getkin,getknf,getknt,getknr
  Use constraints, Only : constraints_tags,apply_shake,&
                          apply_rattle, constraints_type
  Use pmf,         Only : pmf_tags,pmf_type
  Use rigid_bodies,    Only : rigid_bodies_type,getrotmat,no_squish,rigid_bodies_stress
  Use numerics, Only : images,box_mueller_saru2
  Use errors_warnings, Only : error,info
  Use core_shell, Only : core_shell_type
  Use statistics, Only : stats_type
  Use timer, Only :  timer_type
  Use thermostat, Only : adjust_timestep,thermostat_type
  Use core_shell, Only : core_shell_type
  Implicit None

  Private

  Public :: nvt_g0_vv, nvt_g1_vv, nvt_g0_scl, nvt_g1_scl

Contains

  Subroutine nvt_g0_vv(isw,lvar,mndis,mxdis,mxstp,tstep,nstep,degfre,consv, &
      strkin,engke,cshell,cons,pmf,stat,thermo,domain,tmr,parts,comm)

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
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( corePart ),   Intent( InOut ) :: parts(:)
    Type( comms_type ), Intent( InOut ) :: comm


    Logical                 :: safe,lcol,lfst
    Integer                 :: fail(1:9),i
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: chitdr,cintdr
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)
    Character( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
       Call cons%allocate_work(mxatms)
Call pmf%allocate_work()
Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
    End If
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail(7))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail(8))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_g0 allocation failure'
       Call error(0,message)
    End If


    If (thermo%newjob) Then
       thermo%newjob = .false.

  ! inertia parameter for Nose-Hoover thermostat

       thermo%qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
       thermo%ceng  = 2.0_wp*thermo%sigma

  ! set number of constraint+pmf shake iterations

       If (cons%megcon > 0 .or.  pmf%megpmf > 0) thermo%mxkit=1
       If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit=cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Then
         Call constraints_tags(lstitr,cons,parts,comm)
       End If

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0) Then
         Call pmf_tags(lstitr,pmf,parts,comm)
       End If
    End If

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! store initial values

       Do i=1,natms
          xxt(i) = parts(i)%xxx
          yyt(i) = parts(i)%yyy
          zzt(i) = parts(i)%zzz

          vxt(i) = vxx(i)
          vyt(i) = vyy(i)
          vzt(i) = vzz(i)

          fxt(i) = parts(i)%fxx
          fyt(i) = parts(i)%fyy
          fzt(i) = parts(i)%fzz
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

       If (pmf%megpmf > 0) Then
          stat%virpmf=0.0_wp
          stat%strpmf=0.0_wp
       End If

  ! integrate and apply nvt_g0_scl thermostat - 1/2 step

       Call nvt_g0_scl &
             (hstep,degfre,isw,nstep,thermo%ceng,thermo%qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,engke,thermo,comm)

  ! update velocity and position

       Do i=1,natms
          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*fxt(i)
             vyy(i)=vyy(i)+tmp*fyt(i)
             vzz(i)=vzz(i)+tmp*fzt(i)

             parts(i)%xxx=xxt(i)+tstep*vxx(i)
             parts(i)%yyy=yyt(i)+tstep*vyy(i)
             parts(i)%zzz=zzt(i)+tstep*vzz(i)
          End If
       End Do

  ! SHAKE procedures

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep,thermo%mxkit,thermo%kit,oxt,oyt,ozt,&
          lstitr,stat,pmf,cons,domain,tmr,parts,comm)
       End If

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,parts,&
 xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
            Call info(message,.true.)

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
             vxx(i)=vxx(i)+tmp*parts(i)%fxx
             vyy(i)=vyy(i)+tmp*parts(i)%fyy
             vzz(i)=vzz(i)+tmp*parts(i)%fzz
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,thermo%kit,pmf,cons,stat,domain,tmr,comm)
       End If

  ! integrate and apply nvt_g0_scl thermostat - 1/2 step

       Call nvt_g0_scl &
             (hstep,degfre,isw,nstep,thermo%ceng,thermo%qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,engke,thermo,comm)

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*thermo%qmass*thermo%chi_t**2 + thermo%ceng*thermo%cint

  ! kinetic contribution to stress tensor

       Call kinstress(vxx,vyy,vzz,strkin,comm)

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,           Stat=fail(1))
       Call cons%deallocate_work()
Call pmf%deallocate_work()
Deallocate (oxt,oyt,ozt,       Stat=fail( 6))
    End If
    Deallocate (xxt,yyt,zzt,         Stat=fail(7))
    Deallocate (vxt,vyt,vzt,         Stat=fail(8))
    Deallocate (fxt,fyt,fzt,         Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_g0 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nvt_g0_vv

  Subroutine nvt_g1_vv(isw,lvar,mndis,mxdis,mxstp,tstep,nstep,degfre,consv, &
      strkin,strknf,strknt,engke,engrot,strcom,vircom,cshell,cons,pmf,stat, &
      thermo,rigid,domain,tmr,parts,comm)

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
    Real( Kind = wp ), Intent( InOut ) :: strcom(1:9),vircom
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( corePart ),   Intent( InOut ) :: parts(:)
    Type( comms_type ), Intent( InOut ) :: comm


    Logical                 :: safe,lcol,lfst
    Integer                 :: fail(1:14),matms,i,j,i1,i2, &
                               irgd,jrgd,krgd,lrgd,rgdtyp
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
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       Call cons%allocate_work(mxatms)
Call pmf%allocate_work()
Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
    End If
    Allocate (ggx(1:rigid%max_list*rigid%max_rigid), &
      ggy(1:rigid%max_list*rigid%max_rigid), &
      ggz(1:rigid%max_list*rigid%max_rigid), Stat=fail( 7))
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
    Allocate (q0t(1:rigid%max_rigid), &
      q1t(1:rigid%max_rigid), &
      q2t(1:rigid%max_rigid), &
      q3t(1:rigid%max_rigid),  Stat=fail(11))
    Allocate (rgdxxt(1:rigid%max_rigid), &
      rgdyyt(1:rigid%max_rigid), &
      rgdzzt(1:rigid%max_rigid),      Stat=fail(12))
    Allocate (rgdvxt(1:rigid%max_rigid), &
      rgdvyt(1:rigid%max_rigid), &
      rgdvzt(1:rigid%max_rigid),      Stat=fail(13))
    Allocate (rgdoxt(1:rigid%max_rigid), &
      rgdoyt(1:rigid%max_rigid), &
      rgdozt(1:rigid%max_rigid),      Stat=fail(14))
    If (Any(fail > 0)) Then
       Write(message,'(a,i0)') 'nvt_g1 allocation failure'
       Call error(0,message)
    End If


    If (thermo%newjob) Then
       thermo%newjob = .false.

  ! inertia parameter for Nose-Hoover thermostat

       thermo%qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
       thermo%ceng  = 2.0_wp*thermo%sigma

  ! set number of constraint+pmf shake iterations

       If (cons%megcon > 0 .or.  pmf%megpmf > 0) thermo%mxkit=1
       If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit=cons%max_iter_shake

  ! thermo%unsafe positioning due to possibly locally shared RBs

       thermo%unsafe=(Any(domain%map == comm%idnode))
    End If

  ! set matms

    matms=nlast
    If (comm%mxnode == 1) matms=natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Then
         Call constraints_tags(lstitr,cons,parts,comm)
       End If

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0) Then
         Call pmf_tags(lstitr,pmf,parts,comm)
       End If
    End If

  ! Get the RB particles vectors wrt the RB's COM

    krgd=0
    Do irgd=1,rigid%n_types
       rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

       lrgd=rigid%list(-1,irgd)
       If (rigid%frozen(0,rgdtyp) < lrgd) Then
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

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! store initial values

       Do i=1,matms
          xxt(i) = parts(i)%xxx
          yyt(i) = parts(i)%yyy
          zzt(i) = parts(i)%zzz

          vxt(i) = vxx(i)
          vyt(i) = vyy(i)
          vzt(i) = vzz(i)

          fxt(i) = parts(i)%fxx
          fyt(i) = parts(i)%fyy
          fzt(i) = parts(i)%fzz
       End Do

       Do irgd=1,rigid%n_types
          q0t(irgd)=rigid%q0(irgd)
          q1t(irgd)=rigid%q1(irgd)
          q2t(irgd)=rigid%q2(irgd)
          q3t(irgd)=rigid%q3(irgd)

          rgdxxt(irgd) = rigid%xxx(irgd)
          rgdyyt(irgd) = rigid%yyy(irgd)
          rgdzzt(irgd) = rigid%zzz(irgd)

          rgdvxt(irgd) = rigid%vxx(irgd)
          rgdvyt(irgd) = rigid%vyy(irgd)
          rgdvzt(irgd) = rigid%vzz(irgd)

          rgdoxt(irgd) = rigid%oxx(irgd)
          rgdoyt(irgd) = rigid%oyy(irgd)
          rgdozt(irgd) = rigid%ozz(irgd)
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

       If (pmf%megpmf > 0) Then
          stat%virpmf=0.0_wp
          stat%strpmf=0.0_wp
       End If

  ! integrate and apply nvt_g1_scl thermostat - 1/2 step

       Call nvt_g1_scl &
             (hstep,degfre,isw,nstep,thermo%ceng,thermo%qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,                                                &
             engke,engrot,thermo,rigid,comm)

  ! update velocity and position of FPs

       Do j=1,nfree
          i=lstfre(j)

          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*fxt(i)
             vyy(i)=vyy(i)+tmp*fyt(i)
             vzz(i)=vzz(i)+tmp*fzt(i)

             parts(i)%xxx=xxt(i)+tstep*vxx(i)
             parts(i)%yyy=yyt(i)+tstep*vyy(i)
             parts(i)%zzz=zzt(i)+tstep*vzz(i)
          End If
       End Do

  ! SHAKE procedures

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep,thermo%mxkit,thermo%kit,oxt,oyt,ozt,&
          lstitr,stat,pmf,cons,domain,tmr,parts,comm)
       End If

  ! update velocity and position of RBs

       krgd=0
       Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) < lrgd) Then

  ! calculate COM force and torque

             fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
             tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rigid%frozen(0,rgdtyp) == 0) Then
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

             If (rigid%frozen(0,rgdtyp) > 1) Then
                i1=rigid%index_local(rigid%index_global(1,rgdtyp),irgd)
                i2=rigid%index_local(rigid%index_global(2,rgdtyp),irgd)

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

             opx=rigid%oxx(irgd)*rigid%rix(1,rgdtyp)
             opy=rigid%oyy(irgd)*rigid%riy(1,rgdtyp)
             opz=rigid%ozz(irgd)*rigid%riz(1,rgdtyp)

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
             (tstep,rigid%rix(2,rgdtyp),rigid%riy(2,rgdtyp),rigid%riz(2,rgdtyp), &
             rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),p0,p1,p2,p3)
             Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

  ! update RB angular & COM velocities to half step

             opx=0.5_wp*(-rigid%q1(irgd)*p0+rigid%q0(irgd)*p1+rigid%q3(irgd)*p2-rigid%q2(irgd)*p3)
             opy=0.5_wp*(-rigid%q2(irgd)*p0-rigid%q3(irgd)*p1+rigid%q0(irgd)*p2+rigid%q1(irgd)*p3)
             opz=0.5_wp*(-rigid%q3(irgd)*p0+rigid%q2(irgd)*p1-rigid%q1(irgd)*p2+rigid%q0(irgd)*p3)

             rigid%oxx(irgd)=opx*rigid%rix(2,rgdtyp)
             rigid%oyy(irgd)=opy*rigid%riy(2,rgdtyp)
             rigid%ozz(irgd)=opz*rigid%riz(2,rgdtyp)

             tmp=hstep/rigid%weight(0,rgdtyp)
             rigid%vxx(irgd)=rigid%vxx(irgd)+tmp*fmx
             rigid%vyy(irgd)=rigid%vyy(irgd)+tmp*fmy
             rigid%vzz(irgd)=rigid%vzz(irgd)+tmp*fmz

  ! update RB COM to full step

             rigid%xxx(irgd)=rgdxxt(irgd)+tstep*rigid%vxx(irgd)
             rigid%yyy(irgd)=rgdyyt(irgd)+tstep*rigid%vyy(irgd)
             rigid%zzz(irgd)=rgdzzt(irgd)+tstep*rigid%vzz(irgd)

  ! update RB members positions and halfstep velocities

             Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then
                   If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                      x(1)=rigid%x(jrgd,rgdtyp)
                      y(1)=rigid%y(jrgd,rgdtyp)
                      z(1)=rigid%z(jrgd,rgdtyp)

  ! new atomic positions

                      parts(i)%xxx=rot(1)*x(1)+rot(2)*y(1)+rot(3)*z(1) + rigid%xxx(irgd)
                      parts(i)%yyy=rot(4)*x(1)+rot(5)*y(1)+rot(6)*z(1) + rigid%yyy(irgd)
                      parts(i)%zzz=rot(7)*x(1)+rot(8)*y(1)+rot(9)*z(1) + rigid%zzz(irgd)

  ! new atomic velocities in body frame

                      vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                      vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                      vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

  ! DD bound positions

                      If (thermo%unsafe) Then
                         x(1)=parts(i)%xxx-xxt(i)
                         y(1)=parts(i)%yyy-yyt(i)
                         z(1)=parts(i)%zzz-zzt(i)
                         Call images(imcon,cell,1,x,y,z)
                         parts(i)%xxx=x(1)+xxt(i)
                         parts(i)%yyy=y(1)+yyt(i)
                         parts(i)%zzz=z(1)+zzt(i)
                      End If

  ! new atomic velocities in lab frame

                      vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                      vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                      vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                   End If
                End If
             End Do

          End If
       End Do

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,parts,&
 xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
            Call info(message,.true.)


  ! restore initial conditions

             thermo%chi_t=chitdr
             thermo%cint=cintdr

             Do i=1,matms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do

             Do irgd=1,rigid%n_types
                rigid%q0(irgd)=q0t(irgd)
                rigid%q1(irgd)=q1t(irgd)
                rigid%q2(irgd)=q2t(irgd)
                rigid%q3(irgd)=q3t(irgd)

                rigid%vxx(irgd) = rgdvxt(irgd)
                rigid%vyy(irgd) = rgdvyt(irgd)
                rigid%vzz(irgd) = rgdvzt(irgd)

                rigid%oxx(irgd) = rgdoxt(irgd)
                rigid%oyy(irgd) = rgdoyt(irgd)
                rigid%ozz(irgd) = rgdozt(irgd)
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
             vxx(i)=vxx(i)+tmp*parts(i)%fxx
             vyy(i)=vyy(i)+tmp*parts(i)%fyy
             vzz(i)=vzz(i)+tmp*parts(i)%fzz
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,thermo%kit,pmf,cons,stat,domain,tmr,comm)
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stress(strcom,ggx,ggy,ggz,rigid,parts,comm)
       vircom=-(strcom(1)+strcom(5)+strcom(9))

  ! update velocity of RBs

       krgd=0
       Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) < lrgd) Then ! Not that it matters

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

  ! calculate quaternion torques

             qt0=2.0_wp*(-rigid%q1(irgd)*trx-rigid%q2(irgd)*try-rigid%q3(irgd)*trz)
             qt1=2.0_wp*( rigid%q0(irgd)*trx-rigid%q3(irgd)*try+rigid%q2(irgd)*trz)
             qt2=2.0_wp*( rigid%q3(irgd)*trx+rigid%q0(irgd)*try-rigid%q1(irgd)*trz)
             qt3=2.0_wp*(-rigid%q2(irgd)*trx+rigid%q1(irgd)*try+rigid%q0(irgd)*trz)

  ! recover quaternion momenta at half time step

             opx=rigid%oxx(irgd)*rigid%rix(1,rgdtyp)
             opy=rigid%oyy(irgd)*rigid%riy(1,rgdtyp)
             opz=rigid%ozz(irgd)*rigid%riz(1,rgdtyp)

             p0=2.0_wp*(-rigid%q1(irgd)*opx-rigid%q2(irgd)*opy-rigid%q3(irgd)*opz)
             p1=2.0_wp*( rigid%q0(irgd)*opx-rigid%q3(irgd)*opy+rigid%q2(irgd)*opz)
             p2=2.0_wp*( rigid%q3(irgd)*opx+rigid%q0(irgd)*opy-rigid%q1(irgd)*opz)
             p3=2.0_wp*(-rigid%q2(irgd)*opx+rigid%q1(irgd)*opy+rigid%q0(irgd)*opz)

  ! update quaternion momenta to full step

             p0=p0+hstep*qt0
             p1=p1+hstep*qt1
             p2=p2+hstep*qt2
             p3=p3+hstep*qt3

  ! update RB angular & COM velocities to full step

             opx=0.5_wp*(-rigid%q1(irgd)*p0+rigid%q0(irgd)*p1+rigid%q3(irgd)*p2-rigid%q2(irgd)*p3)
             opy=0.5_wp*(-rigid%q2(irgd)*p0-rigid%q3(irgd)*p1+rigid%q0(irgd)*p2+rigid%q1(irgd)*p3)
             opz=0.5_wp*(-rigid%q3(irgd)*p0+rigid%q2(irgd)*p1-rigid%q1(irgd)*p2+rigid%q0(irgd)*p3)

             rigid%oxx(irgd)=opx*rigid%rix(2,rgdtyp)
             rigid%oyy(irgd)=opy*rigid%riy(2,rgdtyp)
             rigid%ozz(irgd)=opz*rigid%riz(2,rgdtyp)

             tmp=hstep/rigid%weight(0,rgdtyp)
             rigid%vxx(irgd)=rigid%vxx(irgd)+tmp*fmx
             rigid%vyy(irgd)=rigid%vyy(irgd)+tmp*fmy
             rigid%vzz(irgd)=rigid%vzz(irgd)+tmp*fmz

  ! update RB members velocities

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

  ! integrate and apply nvt_g1_scl thermostat - 1/2 step

       Call nvt_g1_scl &
             (hstep,degfre,isw,nstep,thermo%ceng,thermo%qmass,0.0_wp,0.0_wp, &
             vxx,vyy,vzz,                                                &
             engke,engrot,thermo,rigid,comm)

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*thermo%qmass*thermo%chi_t**2 + thermo%ceng*thermo%cint

  ! update kinetic energy and stress

       Call kinstresf(vxx,vyy,vzz,strknf,comm)
       Call kinstrest(rigid,strknt,comm)

       strkin=strknf+strknt

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,            Stat=fail( 1))
       Call cons%deallocate_work()
Call pmf%deallocate_work()
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
             engke,engrot,thermo,rigid,comm)

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
    Real( Kind = wp ),                        Intent(   Out ) :: engke,engrot
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
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
    engkt=getknt(rigid,comm)

    engke=engkf+engkt

    engrot=getknr(rigid,comm)

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

    Do irgd=1,rigid%n_types
       rigid%vxx(irgd)=scale*rigid%vxx(irgd)
       rigid%vyy(irgd)=scale*rigid%vyy(irgd)
       rigid%vzz(irgd)=scale*rigid%vzz(irgd)

       rigid%oxx(irgd)=scale*rigid%oxx(irgd)
       rigid%oyy(irgd)=scale*rigid%oyy(irgd)
       rigid%ozz(irgd)=scale*rigid%ozz(irgd)
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
