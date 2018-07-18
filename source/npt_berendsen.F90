Module npt_berendsen
  Use kinds,           Only : wp
  Use comms,           Only : comms_type,gmax
  Use configuration,   Only : imcon,cell,volm,natms,nlast,nfree, &
                              lfrzn,lstfre,weight,               &
                              xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use domains,         Only : domains_type
  Use setup
  Use site, Only : site_type
  Use kinetics,        Only : getvom
  Use constraints,     Only : constraints_type,apply_rattle,&
            constraints_tags,apply_shake
  Use statistics, Only : stats_type
  Use pmf,             Only : pmf_tags,pmf_type
  Use rigid_bodies, Only : rigid_bodies_type,getrotmat,no_squish,rigid_bodies_stress
  Use nvt_berendsen,   Only : nvt_b0_scl,nvt_b1_scl
  Use errors_warnings, Only : error,info
  Use thermostat, Only : thermostat_type
  Use core_shell, Only : core_shell_type
  Use timer, Only: timer_type
  Use thermostat, Only : adjust_timestep
  Use vdw, Only : vdw_type
  Use numerics, Only : images
  Implicit None

  Private

  Public :: npt_b0_vv, npt_b1_vv

Contains

  Subroutine npt_b0_vv(isw,lvar,mndis,mxdis,mxstp,tstep, &
             virtot,                            &
             strkin,engke,                      &
             cshell,cons,pmf,stat,thermo,sites,vdws,domain,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian equations of motion in
  ! molecular dynamics - velocity verlet with Berendsen thermostat and
  ! isotropic pressure control (not symplectic)
  !
  ! isothermal compressibility (beta) set to that of liquid water
  ! = 0.007372 dl_poly units
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: isw
    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep
    Real( Kind = wp ),  Intent( In    ) :: virtot
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( pmf_type), Intent( InOut ) :: pmf
    Type( constraints_type), Intent( InOut ) :: cons
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true.
    Logical                 :: safe,lcol,lfst
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:9),iter,i
    Real( Kind = wp ), Save :: volm0,elrc0,virlrc0
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: czero(1:9)
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                               scale,vom(1:3),pr

    Real( Kind = wp ), Parameter :: beta = 7.3728e-3_wp


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len = 256 ) :: message

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
       Write(message,'(a)') 'npt_b0 allocation failure'
       Call error(0,message)
    End If


    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       volm0   = volm
       elrc0   = vdws%elrc
       virlrc0 = vdws%vlrc

       Allocate (dens0(1:mxatyp), Stat=fail(1))
       If (fail(1) > 0) Then
          Write(message,'(a)') 'dens0 allocation failure'
          Call error(0,message)
       End If
       Do i=1,sites%ntype_atom
          dens0(i) = sites%dens(i)
       End Do

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+12
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0)Then
         Call constraints_tags(lstitr,cons,comm)
       End If


  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0)Then
         Call pmf_tags(lstitr,pmf,comm)
       End If
    End If

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep

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

  ! store temporary cell parameters

       czero=cell

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

  ! iterate forces, stat%vircon and thermo%chi_p

       Do iter=1,mxiter

  ! Berendsen barostat and thermostat are not coupled
  ! calculate system pressure: iterate stat%vircon and stat%virpmf

          pr = (2.0_wp*engke-virtot-stat%vircon-stat%virpmf) / (3.0_wp*volm)

  ! pressure control variable

          thermo%chi_p = 1.0_wp + beta*tstep*(pr-thermo%press)/thermo%tau_p
          scale = thermo%chi_p**(1.0_wp/3.0_wp)

  ! update velocity and position

          Do i=1,natms
             If (weight(i) > 1.0e-6_wp) Then
                tmp=hstep/weight(i)
                vxx(i)=vxt(i)+tmp*fxx(i)
                vyy(i)=vyt(i)+tmp*fyy(i)
                vzz(i)=vzt(i)+tmp*fzz(i)

                xxx(i)=scale*xxt(i)+tstep*vxx(i)
                yyy(i)=scale*yyt(i)+tstep*vyy(i)
                zzz(i)=scale*zzt(i)+tstep*vzz(i)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
! update cell parameters: isotropic

             cell=scale*czero
             Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
               lstitr,stat,pmf,cons,domain,tmr,comm)
          End If

       End Do

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,xxx,yyy,zzz,&
 xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
            Call info(message,.true.)


  ! restore initial conditions

             Do i=1,natms
                fxx(i) = fxt(i)
                fyy(i) = fyt(i)
                fzz(i) = fzt(i)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! update cell parameters: isotropic

       If (cons%megcon == 0 .and. pmf%megpmf == 0) cell=scale*czero

  ! update volume: thermo%chi_p=scale^3

       volm=thermo%chi_p*volm

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       vdws%elrc=elrc0*tmp
       vdws%vlrc=virlrc0*tmp
       Do i=1,sites%ntype_atom
          sites%dens(i)=dens0(i)*tmp
       End Do

  ! construct a 'mock' scaling tensor for xscale

       Do i=2,8
          thermo%eta(i)=0.0_wp
       End Do
       thermo%eta(1)=scale
       thermo%eta(5)=scale
       thermo%eta(9)=scale

  ! second pass of velocity verlet algorithm

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

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,kit,pmf,cons,stat,domain,tmr,comm)
       End If

  ! integrate and apply nvt_b0_scl thermostat - full step

       Call nvt_b0_scl(1,tstep,vxx,vyy,vzz,strkin,engke,thermo,comm)

  ! remove system centre of mass velocity

       Call getvom(vom,vxx,vyy,vzz,comm)

       Do i=1,natms
          If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
             vxx(i)=vxx(i)-vom(1)
             vyy(i)=vyy(i)-vom(2)
             vzz(i)=vzz(i)-vom(3)
          End If
       End Do

  ! update kinetic energy and stress

       Call nvt_b0_scl(0,tstep,vxx,vyy,vzz,strkin,engke,thermo,comm)

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,           Stat=fail(1))
       Call cons%deallocate_work()
       Call pmf%deallocate_work()
       Deallocate (oxt,oyt,ozt,      Stat=fail(6))
    End If
    Deallocate (xxt,yyt,zzt,         Stat=fail(7))
    Deallocate (vxt,vyt,vzt,         Stat=fail(8))
    Deallocate (fxt,fyt,fzt,         Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'npt_b0 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine npt_b0_vv

  Subroutine npt_b1_vv(isw,lvar,mndis,mxdis,mxstp,tstep, &
             virtot,                            &
             strkin,strknf,strknt,engke,engrot, &
             strcom,vircom,                     &
             cshell,cons,pmf,stat,thermo,sites,vdws,rigid,domain,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
  ! RBs, equations of motion in molecular dynamics
  ! - velocity verlet with Berendsen thermostat and isotropic pressure
  ! control (not symplectic)
  !
  ! isothermal compressibility (beta) set to that of liquid water
  ! = 0.007372 dl_poly units
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: isw
    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep
    Real( Kind = wp ),  Intent( In    ) :: virtot
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke, &
                                           strknf(1:9),strknt(1:9),engrot
    Real( Kind = wp ),  Intent( InOut ) :: strcom(1:9),vircom
    Type( pmf_type), Intent( InOut ) :: pmf
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true. , &
                               unsafe = .false.
    Logical                 :: safe,lcol,lfst
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:15),matms,iter,i,j,i1,i2, &
                               irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ), Save :: volm0,elrc0,virlrc0
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: czero(1:9)
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                               scale,vom(1:3),pr
    Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                               opx,opy,opz,fmx,fmy,fmz,       &
                               tqx,tqy,tqz,trx,try,trz,       &
                               qt0,qt1,qt2,qt3,               &
                               vpx,vpy,vpz

    Real( Kind = wp ), Parameter :: beta = 7.3728e-3_wp


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

    Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
    Real( Kind = wp ), Allocatable :: q0t(:),q1t(:),q2t(:),q3t(:)
    Real( Kind = wp ), Allocatable :: p0(:),p1(:),p2(:),p3(:)
    Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
    Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
    Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       Call cons%allocate_work(mxatms)
       Call pmf%allocate_work()
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail( 6))
    End If
    Allocate (ggx(1:rigid%max_list*rigid%max_rigid), &
      ggy(1:rigid%max_list*rigid%max_rigid), &
      ggz(1:rigid%max_list*rigid%max_rigid), &
                                                                    Stat=fail( 7))
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
    Allocate (q0t(1:rigid%max_rigid), &
      q1t(1:rigid%max_rigid), &
      q2t(1:rigid%max_rigid), &
      q3t(1:rigid%max_rigid),  Stat=fail(11))
    Allocate (p0(1:rigid%max_rigid), &
      p1(1:rigid%max_rigid), &
      p2(1:rigid%max_rigid), &
      p3(1:rigid%max_rigid),      Stat=fail(12))
    Allocate (rgdxxt(1:rigid%max_rigid), &
      rgdyyt(1:rigid%max_rigid), &
      rgdzzt(1:rigid%max_rigid),      Stat=fail(13))
    Allocate (rgdvxt(1:rigid%max_rigid), &
      rgdvyt(1:rigid%max_rigid), &
      rgdvzt(1:rigid%max_rigid),      Stat=fail(14))
    Allocate (rgdoxt(1:rigid%max_rigid), &
      rgdoyt(1:rigid%max_rigid), &
      rgdozt(1:rigid%max_rigid),      Stat=fail(15))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'npt_b1 allocation failure'
       Call error(0,message)
    End If


    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       volm0   = volm
       elrc0   = vdws%elrc
       virlrc0 = vdws%vlrc

       Allocate (dens0(1:mxatyp), Stat=fail(1))
       If (fail(1) > 0) Then
          Write(message,'(a)') 'dens0 allocation failure'
          Call error(0,message)
       End If
       Do i=1,sites%ntype_atom
          dens0(i) = sites%dens(i)
       End Do

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+12
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake

  ! unsafe positioning due to possibly locally shared RBs

       unsafe=(Any(domain%map == comm%idnode))
    End If

  ! set matms

    matms=nlast
    If (comm%mxnode == 1) matms=natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0)Then
         Call constraints_tags(lstitr,cons,comm)
       End If
  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0)Then
         Call pmf_tags(lstitr,pmf,comm)
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

             ggx(krgd)=xxx(i)-rigid%xxx(irgd)
             ggy(krgd)=yyy(i)-rigid%yyy(irgd)
             ggz(krgd)=zzz(i)-rigid%zzz(irgd)
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

  ! store temporary cell parameters

       czero=cell

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

                Call images(imcon,czero,1,x,y,z)

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

  ! recover quaternion momenta at start of time step

             opx=rgdoxt(irgd)*rigid%rix(1,rgdtyp)
             opy=rgdoyt(irgd)*rigid%riy(1,rgdtyp)
             opz=rgdozt(irgd)*rigid%riz(1,rgdtyp)

             p0(irgd)=2.0_wp*(-q1t(irgd)*opx-q2t(irgd)*opy-q3t(irgd)*opz)
             p1(irgd)=2.0_wp*( q0t(irgd)*opx-q3t(irgd)*opy+q2t(irgd)*opz)
             p2(irgd)=2.0_wp*( q3t(irgd)*opx+q0t(irgd)*opy-q1t(irgd)*opz)
             p3(irgd)=2.0_wp*(-q2t(irgd)*opx+q1t(irgd)*opy+q0t(irgd)*opz)

  ! update quaternion momenta to half step

             p0(irgd)=p0(irgd)+hstep*qt0
             p1(irgd)=p1(irgd)+hstep*qt1
             p2(irgd)=p2(irgd)+hstep*qt2
             p3(irgd)=p3(irgd)+hstep*qt3

  ! update RB COM velocities to half step

             tmp=hstep/rigid%weight(0,rgdtyp)
             rigid%vxx(irgd)=rgdvxt(irgd)+tmp*fmx
             rigid%vyy(irgd)=rgdvyt(irgd)+tmp*fmy
             rigid%vzz(irgd)=rgdvzt(irgd)+tmp*fmz
          End If
       End Do

  ! iterate forces, stat%vircon and thermo%chi_p

       Do iter=1,mxiter

  ! Berendsen barostat and thermostat are not coupled
  ! calculate system pressure: iterate stat%vircon and stat%virpmf

          pr = (2.0_wp*engke-virtot-stat%vircon-stat%virpmf-vircom) / (3.0_wp*volm)

  ! pressure control variable

          thermo%chi_p = 1.0_wp + beta*tstep*(pr-thermo%press)/thermo%tau_p
          scale = thermo%chi_p**(1.0_wp/3.0_wp)

  ! update cell parameters: isotropic

          cell=scale*czero

  ! update velocity and position of FPs

          Do j=1,nfree
             i=lstfre(j)

             If (weight(i) > 1.0e-6_wp) Then
                tmp=hstep/weight(i)
                vxx(i)=vxt(i)+tmp*fxx(i)
                vyy(i)=vyt(i)+tmp*fyy(i)
                vzz(i)=vzt(i)+tmp*fzz(i)

                xxx(i)=scale*xxt(i)+tstep*vxx(i)
                yyy(i)=scale*yyt(i)+tstep*vyy(i)
                zzz(i)=scale*zzt(i)+tstep*vzz(i)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
            Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
             lstitr,stat,pmf,cons,domain,tmr,comm)
          End If

       End Do

  ! update velocity and position of RBs

       Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) < lrgd) Then

  ! rotate RB quaternions - update q to full timestep & amend p
  ! and get new rotation matrix

             Call no_squish                                             &
             (tstep,rigid%rix(2,rgdtyp),rigid%riy(2,rgdtyp),rigid%riz(2,rgdtyp), &
             rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd), &
             p0(irgd),p1(irgd),p2(irgd),p3(irgd))
             Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

  ! update RB angular velocities to half step

             opx=0.5_wp*(-rigid%q1(irgd)*p0(irgd)+ &
               rigid%q0(irgd)*p1(irgd)+ &
               rigid%q3(irgd)*p2(irgd)- &
               rigid%q2(irgd)*p3(irgd))
             opy=0.5_wp*(-rigid%q2(irgd)*p0(irgd)- &
               rigid%q3(irgd)*p1(irgd)+ &
               rigid%q0(irgd)*p2(irgd)+ &
               rigid%q1(irgd)*p3(irgd))
             opz=0.5_wp*(-rigid%q3(irgd)*p0(irgd)+ &
               rigid%q2(irgd)*p1(irgd)- &
               rigid%q1(irgd)*p2(irgd)+ &
               rigid%q0(irgd)*p3(irgd))

             rigid%oxx(irgd)=opx*rigid%rix(2,rgdtyp)
             rigid%oyy(irgd)=opy*rigid%riy(2,rgdtyp)
             rigid%ozz(irgd)=opz*rigid%riz(2,rgdtyp)

  ! update RB COM to full step

             rigid%xxx(irgd)=scale*rgdxxt(irgd)+tstep*rigid%vxx(irgd)
             rigid%yyy(irgd)=scale*rgdyyt(irgd)+tstep*rigid%vyy(irgd)
             rigid%zzz(irgd)=scale*rgdzzt(irgd)+tstep*rigid%vzz(irgd)

  ! update RB members positions and halfstep velocities

             Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then
                   If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                      x(1)=rigid%x(jrgd,rgdtyp)
                      y(1)=rigid%y(jrgd,rgdtyp)
                      z(1)=rigid%z(jrgd,rgdtyp)

  ! new atomic positions

                      xxx(i)=rot(1)*x(1)+rot(2)*y(1)+rot(3)*z(1) + rigid%xxx(irgd)
                      yyy(i)=rot(4)*x(1)+rot(5)*y(1)+rot(6)*z(1) + rigid%yyy(irgd)
                      zzz(i)=rot(7)*x(1)+rot(8)*y(1)+rot(9)*z(1) + rigid%zzz(irgd)

  ! new atomic velocities in body frame

                      vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                      vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                      vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

  ! DD bound positions

                      If (unsafe) Then
                         vxx(i)=scale*xxt(i)
                         vyy(i)=scale*yyt(i)
                         vzz(i)=scale*zzt(i)

                         x(1)=xxx(i)-vxx(i)
                         y(1)=yyy(i)-vyy(i)
                         z(1)=zzz(i)-vzz(i)
                         Call images(imcon,cell,1,x,y,z)
                         xxx(i)=x(1)+vxx(i)
                         yyy(i)=y(1)+vyy(i)
                         zzz(i)=z(1)+vzz(i)
                      End If

  ! new atomic velocities in lab frame

                      vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                      vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                      vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                   Else
                      x(1)=rigid%xxx(irgd)-rgdxxt(irgd)
                      y(1)=rigid%yyy(irgd)-rgdyyt(irgd)
                      z(1)=rigid%zzz(irgd)-rgdzzt(irgd)
                      If (unsafe) Call images(imcon,cell,1,x,y,z) ! DD bound positions
                      xxx(i)=xxt(i)+x(1)
                      yyy(i)=yyt(i)+y(1)
                      zzz(i)=zzt(i)+z(1)
                   End If
                End If
             End Do

          Else

  ! update RB COM to full step

             rigid%xxx(irgd)=scale*rgdxxt(irgd)
             rigid%yyy(irgd)=scale*rgdyyt(irgd)
             rigid%zzz(irgd)=scale*rgdzzt(irgd)

             Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then
                   x(1)=rigid%xxx(irgd)-rgdxxt(irgd)
                   y(1)=rigid%yyy(irgd)-rgdyyt(irgd)
                   z(1)=rigid%zzz(irgd)-rgdzzt(irgd)
                   If (unsafe) Call images(imcon,cell,1,x,y,z) ! DD bound positions
                   xxx(i)=xxt(i)+x(1)
                   yyy(i)=yyt(i)+y(1)
                   zzz(i)=zzt(i)+z(1)
                End If
             End Do

          End If
       End Do

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,xxx,yyy,zzz,&
 xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
            Call info(message,.true.)


  ! restore initial conditions

             Do i=1,matms
                fxx(i) = fxt(i)
                fyy(i) = fyt(i)
                fzz(i) = fzt(i)
             End Do

             Do irgd=1,rigid%n_types
                rigid%q0(irgd)=q0t(irgd)
                rigid%q1(irgd)=q1t(irgd)
                rigid%q2(irgd)=q2t(irgd)
                rigid%q3(irgd)=q3t(irgd)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! update volume: thermo%chi_p=scale^3

       volm=thermo%chi_p*volm

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       vdws%elrc=elrc0*tmp
       vdws%vlrc=virlrc0*tmp
       Do i=1,sites%ntype_atom
          sites%dens(i)=dens0(i)*tmp
       End Do

  ! construct a 'mock' scaling tensor for xscale

       Do i=2,8
          thermo%eta(i)=0.0_wp
       End Do
       thermo%eta(1)=scale
       thermo%eta(5)=scale
       thermo%eta(9)=scale

  ! second pass of velocity verlet algorithm

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

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,kit,pmf,cons,stat,domain,tmr,comm)
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stress(strcom,ggx,ggy,ggz,rigid,comm)
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

             If (rigid%frozen(0,rgdtyp) > 1) Then
                i1=rigid%index_local(rigid%index_global(1,rgdtyp),irgd)
                i2=rigid%index_local(rigid%index_global(2,rgdtyp),irgd)

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

             p0(irgd)=2.0_wp*(-rigid%q1(irgd)*opx-rigid%q2(irgd)*opy-rigid%q3(irgd)*opz)
             p1(irgd)=2.0_wp*( rigid%q0(irgd)*opx-rigid%q3(irgd)*opy+rigid%q2(irgd)*opz)
             p2(irgd)=2.0_wp*( rigid%q3(irgd)*opx+rigid%q0(irgd)*opy-rigid%q1(irgd)*opz)
             p3(irgd)=2.0_wp*(-rigid%q2(irgd)*opx+rigid%q1(irgd)*opy+rigid%q0(irgd)*opz)

  ! update quaternion momenta to full step

             p0(irgd)=p0(irgd)+hstep*qt0
             p1(irgd)=p1(irgd)+hstep*qt1
             p2(irgd)=p2(irgd)+hstep*qt2
             p3(irgd)=p3(irgd)+hstep*qt3

  ! update RB angular & COM velocities to full step

             opx=0.5_wp*(-rigid%q1(irgd)*p0(irgd)+ &
               rigid%q0(irgd)*p1(irgd)+ &
               rigid%q3(irgd)*p2(irgd)- &
               rigid%q2(irgd)*p3(irgd))
             opy=0.5_wp*(-rigid%q2(irgd)*p0(irgd)- &
               rigid%q3(irgd)*p1(irgd)+ &
               rigid%q0(irgd)*p2(irgd)+ &
               rigid%q1(irgd)*p3(irgd))
             opz=0.5_wp*(-rigid%q3(irgd)*p0(irgd)+ &
               rigid%q2(irgd)*p1(irgd)- &
               rigid%q1(irgd)*p2(irgd)+ &
               rigid%q0(irgd)*p3(irgd))

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

  ! integrate and apply nvt_b1_scl thermostat - full step

       Call nvt_b1_scl &
             (1,tstep,vxx,vyy,vzz,           &
             strkin,strknf,strknt,engke,engrot,thermo,rigid,comm)

  ! remove system centre of mass velocity

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

  ! update kinetic energy and stress

       Call nvt_b1_scl &
             (0,tstep,vxx,vyy,vzz,           &
             strkin,strknf,strknt,engke,engrot,thermo,rigid,comm)

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
    Deallocate (p0,p1,p2,p3,          Stat=fail(12))
    Deallocate (rgdxxt,rgdyyt,rgdzzt, Stat=fail(13))
    Deallocate (rgdvxt,rgdvyt,rgdvzt, Stat=fail(14))
    Deallocate (rgdoxt,rgdoyt,rgdozt, Stat=fail(15))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'npt_b1 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine npt_b1_vv
End Module npt_berendsen
