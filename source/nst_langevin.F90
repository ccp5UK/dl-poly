Module nst_langevin
  Use kinds,           Only : wp, li
  Use comms,           Only : comms_type,gmax
  Use setup
  Use site, Only : site_type
  Use configuration,   Only : imcon,cell,volm,natms,nlast,nfree,  &
                              lsi,lsa,lfrzn,lstfre,weight,        &
                              xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use domains,         Only : map
  Use langevin,        Only : fxl,fyl,fzl,fpl
  Use kinetics,        Only : getvom,getkin,getknt,getknr,getknf, &
                              kinstress,kinstresf,kinstrest
  Use constraints,     Only : apply_rattle,apply_shake,&
                              constraints_tags, constraints_type
  Use pmf,             Only : pmf_tags,pmf_type 
  Use rigid_bodies     
  Use errors_warnings, Only : error,info
  Use shared_units,    Only : update_shared_units
  Use numerics,        Only : dcell, mat_mul,box_mueller_saru6    
  Use langevin,        Only : langevin_forces
  Use nst_nose_hoover, Only : nst_h0_scl,nst_h1_scl
  Use thermostat, Only : thermostat_type
Use core_shell, Only : core_shell_type
  Use statistics, Only : stats_type
  Use timer, Only : timer_type
Use thermostat, Only : adjust_timestep
Use core_shell, Only : core_shell_type
  Implicit None

  Private

  Public :: nst_l0_vv, nst_l1_vv

Contains

  Subroutine nst_l0_vv                          &
             (isw,lvar,mndis,mxdis,mxstp,tstep, &
             nstep,   &
             degfre,stress,             &
             consv,                             &
             strkin,engke,                      &
             elrc,virlrc,cshell,cons,pmf,stat,thermo,site,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian equations of motion in
  ! molecular dynamics - velocity verlet with Langevin thermostat and
  ! barostat (anisotropic pressure control) (symplectic)
  !
  ! Parrinello-Rahman type: changing cell shape
  !
  ! thermo%iso=0 fully anisotropic barostat
  ! thermo%iso=1 semi-isotropic barostat to constant normal pressure & surface area
  ! thermo%iso=2 semi-isotropic barostat to constant normal pressure & surface tension
  !                               or with orthorhombic constraints (thermo%tension=0.0_wp)
  ! thermo%iso=3 semi-isotropic barostat with semi-orthorhombic constraints
  !
  ! reference1: D. Quigley and M.I.J. Probert
  !             J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
  ! reference2: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: isw

    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep

    Integer,            Intent( In    ) :: nstep

    Integer(Kind=li),   Intent( In    ) :: degfre
    Real( Kind = wp ),  Intent( In    ) :: stress(1:9)

    Real( Kind = wp ),  Intent(   Out ) :: consv

    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke


    Real( Kind = wp ),  Intent( InOut ) :: elrc,virlrc
    Type( stats_type), Intent( InOut ) :: stat
Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: site
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true.
    Logical                 :: safe,lcol,lfst
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:9),iter,i
    Real( Kind = wp ), Save :: volm0,elrc0,virlrc0,h_z
    Real( Kind = wp ), Save :: temp,pmass
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: eta0(1:9),engke0
    Real( Kind = wp )       :: cell0(1:9),vzero,celprp(1:10)
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),str1(1:9),mxdr,tmp, &
                               scale,vom(1:3),aaa(1:9),bbb(1:9)

  ! uni1 is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
    uni1(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len = 256 )        :: message

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
       Write(message,'(a)') 'nst_l0 allocation failure'
       Call error(0,message)
    End If


  ! timestep derivatives

    hstep = 0.5_wp*tstep
    qstep = 0.5_wp*hstep
    rstep = 1.0_wp/tstep

    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       volm0   = volm
       elrc0   = elrc
       virlrc0 = virlrc

       Allocate (dens0(1:mxatyp), Stat=fail(1))
       If (fail(1) > 0) Then
          Write(message,'(a)') 'dens0 allocation failure'
          Call error(0,message)
       End If
       Do i=1,site%ntype_atom
          dens0(i) = site%dens(i)
       End Do

  ! Sort thermo%eta for thermo%iso>=1
  ! Initialise and get h_z for thermo%iso>1

       h_z=0
       If      (thermo%iso == 1) Then
          thermo%eta(1:8) = 0.0_wp
       Else If (thermo%iso >  1) Then
          thermo%eta(2:4) = 0.0_wp
          thermo%eta(6:8) = 0.0_wp

          Call dcell(cell,celprp)
          h_z=celprp(9)
       End If

  ! inertia parameter for barostat

       temp  = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
       pmass = ((2.0_wp*thermo%sigma + 3.0_wp*boltz*temp)/3.0_wp) / (2.0_wp*pi*thermo%tai)**2

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+3
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake

  ! Generate Langevin forces for particles and
  ! Langevin tensor force for barostat piston

       fpl=0.0_wp
       Call box_mueller_saru6(Int(degfre/3_li),nstep-1,fpl(1),fpl(2),fpl(3),fpl(4),fpl(5),fpl(6))
       tmp=Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)
       fpl(1:6)=fpl(1:6)*tmp
       fpl(9)=fpl(4)                                 ! Distribute independent
       fpl(4)=fpl(2) ; fpl(7)=fpl(3) ; fpl(8)=fpl(6) ! Symmetrise
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

  ! store current integration variables

       cell0=cell
       vzero=volm
       eta0 =thermo%eta
       engke0=engke

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

       Do iter=1,mxiter

  ! integrate and apply Langevin thermostat - 1/4 step

          scale=Exp(-qstep*thermo%chi)
          Do i=1,natms
             vxx(i)=scale*vxx(i)
             vyy(i)=scale*vyy(i)
             vzz(i)=scale*vzz(i)
          End Do
          engke=engke*scale**2

  ! constraint+pmf virial and stress

          vir=stat%vircon+stat%virpmf
          str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h0_scl barostat - 1/2 step
  ! augment str to include the random tensorial force on the barostat

          str1=str+fpl
          Call nst_h0_scl &
             (1,hstep,degfre,pmass,thermo%tai,volm, &
             h_z,str1,stress,       &
             vxx,vyy,vzz,strkin,engke,thermo,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

          scale=Exp(-qstep*thermo%chi)
          Do i=1,natms
             vxx(i)=scale*vxx(i)
             vyy(i)=scale*vyy(i)
             vzz(i)=scale*vzz(i)
          End Do
          engke=engke*scale**2

  ! update velocities

          Do i=1,natms
             If (weight(i) > 1.0e-6_wp) Then
                tmp=hstep/weight(i)
                vxx(i)=vxx(i)+tmp*(fxx(i)+fxl(i))
                vyy(i)=vyy(i)+tmp*(fyy(i)+fyl(i))
                vzz(i)=vzz(i)+tmp*(fzz(i)+fzl(i))
             End If
          End Do

  ! scale cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

          aaa=tstep*thermo%eta
          Call mat_mul(aaa,aaa,bbb)
          aaa=uni1+aaa+0.5_wp*bbb
          Call mat_mul(aaa,cell0,cell)

  ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

          thermo%chi_p = thermo%eta(1)+thermo%eta(5)+thermo%eta(9)
          volm = volm*Exp(tstep*thermo%chi_p)
          thermo%chi_p = thermo%chi_p / 3.0_wp

  ! update positions: second order taylor expansion of Exp(tstep*thermo%eta)

          Do i=1,natms
             If (weight(i) > 1.0e-6_wp) Then
                xxx(i)=tstep*vxx(i)+xxt(i)*aaa(1)+yyt(i)*aaa(2)+zzt(i)*aaa(3)
                yyy(i)=tstep*vyy(i)+xxt(i)*aaa(2)+yyt(i)*aaa(5)+zzt(i)*aaa(6)
                zzz(i)=tstep*vzz(i)+xxt(i)*aaa(3)+yyt(i)*aaa(6)+zzt(i)*aaa(9)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
           Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
          lstitr,&
          stat,pmf,cons,tmr,comm)
          End If

  ! restore original integration parameters as well as
  ! velocities if iter < mxiter
  ! in the next iteration stat%strcon and stat%strpmf are freshly new

          If (iter < mxiter) Then
             volm=vzero
             thermo%eta =eta0
             engke=engke0

             Do i=1,natms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do
          End If
       End Do

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,xxx,yyy,zzz,&
 xxt,yyt,zzt,cshell%legshl,message,tmp,comm)) Then 
            Call info(message,.true.)

  ! scale Langevin random forces

             Do i=1,natms
                fxl(i) = fxl(i)*tmp
                fyl(i) = fyl(i)*tmp
                fzl(i) = fzl(i)*tmp
             End Do
             fpl = fpl*tmp

  ! restore initial conditions

             volm=vzero
             thermo%eta =eta0
             engke=engke0

             Do i=1,natms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)

                fxx(i) = fxt(i)
                fyy(i) = fyt(i)
                fzz(i) = fzt(i)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       elrc=elrc0*tmp
       virlrc=virlrc0*tmp
       Do i=1,site%ntype_atom
          site%dens(i)=dens0(i)*tmp
       End Do

  ! get h_z for thermo%iso>1

       If (thermo%iso > 1) Then
          Call dcell(cell,celprp)
          h_z=celprp(9)
       End If

  ! second stage of velocity verlet algorithm

    Else

  ! Generate Langevin forces for particles and
  ! Langevin tensor force for barostat piston

       Call langevin_forces(nstep,temp,tstep,thermo%chi,fxl,fyl,fzl,cshell)

       fpl=0.0_wp
       Call box_mueller_saru6(Int(degfre/3_li),nstep,fpl(1),fpl(2),fpl(3),fpl(4),fpl(5),fpl(6))
       tmp=Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)
       fpl(1:6)=fpl(1:6)*tmp
       fpl(9)=fpl(4)                                 ! Distribute independent
       fpl(4)=fpl(2) ; fpl(7)=fpl(3) ; fpl(8)=fpl(6) ! Symmetrise

  ! update velocity

       Do i=1,natms
          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*(fxx(i)+fxl(i))
             vyy(i)=vyy(i)+tmp*(fyy(i)+fyl(i))
             vzz(i)=vzz(i)+tmp*(fzz(i)+fzl(i))
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_rattle(tstep,kit,&
                          pmf,cons,stat,tmr,comm)
       End If

  ! update kinetic energy

       engke=getkin(vxx,vyy,vzz,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
       Do i=1,natms
          vxx(i)=scale*vxx(i)
          vyy(i)=scale*vyy(i)
          vzz(i)=scale*vzz(i)
       End Do
       engke=engke*scale**2

  ! constraint+pmf virial and stress

       vir=stat%vircon+stat%virpmf
       str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h0_scl barostat - 1/2 step
  ! augment str to include the random tensorial force on the barostat

       str1=str+fpl
       Call nst_h0_scl &
             (1,hstep,degfre,pmass,thermo%tai,volm, &
             h_z,str1,stress,       &
             vxx,vyy,vzz,strkin,engke,thermo,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
       Do i=1,natms
          vxx(i)=scale*vxx(i)
          vyy(i)=scale*vyy(i)
          vzz(i)=scale*vzz(i)
       End Do

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

       tmp = ( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
         thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*pmass*tmp + thermo%press*volm

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

       Call kinstress(vxx,vyy,vzz,strkin,comm)
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

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
       Write(message,'(a)') 'nst_l0 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nst_l0_vv

  Subroutine nst_l1_vv                          &
             (isw,lvar,mndis,mxdis,mxstp,tstep, &
             nstep,   &
             degfre,degrot,stress,      &
             strkin,strknf,strknt,engke,engrot, &
             consv,                             &
             strcom,vircom,                     &
             elrc,virlrc,cshell,cons,pmf,stat,thermo,site,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
  ! RBs, equations of motion in molecular dynamics
  ! - velocity verlet with Langevin thermostat and
  ! barostat (anisotropic pressure control) (symplectic)
  !
  ! Parrinello-Rahman type: changing cell shape
  !
  ! thermo%iso=0 fully anisotropic barostat
  ! thermo%iso=1 semi-isotropic barostat to constant normal pressure & surface area
  ! thermo%iso=2 semi-isotropic barostat to constant normal pressure & surface tension
  !                               or with orthorhombic constraints (thermo%tension=0.0_wp)
  ! thermo%iso=3 semi-isotropic barostat with semi-orthorhombic constraints
  !
  ! reference1: D. Quigley and M.I.J. Probert
  !             J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
  ! reference2: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: isw

    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep

    Integer,            Intent( In    ) :: nstep

    Integer(Kind=li),   Intent( In    ) :: degfre,degrot
    Real( Kind = wp ),  Intent( In    ) :: stress(1:9)

    Real( Kind = wp ),  Intent(   Out ) :: consv

    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke, &
                                           strknf(1:9),strknt(1:9),engrot


    Real( Kind = wp ),  Intent( InOut ) :: strcom(1:9),vircom

    Real( Kind = wp ),  Intent( InOut ) :: elrc,virlrc
    Type( stats_type), Intent( InOut ) :: stat
Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: site
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true. , &
                               unsafe = .false.
    Logical                 :: safe,lcol,lfst
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:14),matms,iter,i,j,i1,i2, &
                               irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ), Save :: volm0,elrc0,virlrc0,h_z
    Real( Kind = wp ), Save :: pmass,temp
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: eta0(1:9),engke0,engrot0,engknf,engknt
    Real( Kind = wp )       :: cell0(1:9),vzero,celprp(1:10)
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),str1(1:9),mxdr,tmp, &
                               scale,vom(1:3),aaa(1:9),bbb(1:9)
    Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                               opx,opy,opz,fmx,fmy,fmz,       &
                               tqx,tqy,tqz,trx,try,trz,       &
                               fmxl,fmyl,fmzl,                &
                               tqxl,tqyl,tqzl,trxl,tryl,trzl, &
                               qt0l,qt1l,qt2l,qt3l,           &
                               qt0,qt1,qt2,qt3,p0,p1,p2,p3,   &
                               vpx,vpy,vpz

  ! uni1 is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
    uni1(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


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

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len = 256 )        :: message


    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       Call cons%allocate_work(mxatms)
Call pmf%allocate_work()
Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
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
       Write(message,'(a)') 'nst_l1 allocation failure'
       Call error(0,message)
    End If


  ! timestep derivatives

    hstep = 0.5_wp*tstep
    qstep = 0.5_wp*hstep
    rstep = 1.0_wp/tstep

    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       volm0   = volm
       elrc0   = elrc
       virlrc0 = virlrc

       Allocate (dens0(1:mxatyp), Stat=fail(1))
       If (fail(1) > 0) Then
          Write(message,'(a)') 'dens0 allocation failure'
          Call error(0,message)
       End If
       Do i=1,site%ntype_atom
          dens0(i) = site%dens(i)
       End Do

  ! Sort thermo%eta for thermo%iso>=1
  ! Initialise and get h_z for thermo%iso>1

       h_z=0
       If      (thermo%iso == 1) Then
          thermo%eta(1:8) = 0.0_wp
       Else If (thermo%iso >  1) Then
          thermo%eta(2:4) = 0.0_wp
          thermo%eta(6:8) = 0.0_wp

          Call dcell(cell,celprp)
          h_z=celprp(9)
       End If

  ! inertia parameter for barostat

       temp  = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
       pmass = ((Real(degfre-degrot,wp) + 3.0_wp)/3.0_wp)*boltz*temp / (2.0_wp*pi*thermo%tai)**2

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+3
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake

  ! unsafe positioning due to possibly locally shared RBs

       unsafe=(Any(map == comm%idnode))

  ! Generate Langevin forces for particles and
  ! Langevin tensor force for barostat piston

       fpl=0.0_wp
       Call box_mueller_saru6(Int(degfre/3_li),nstep-1,fpl(1),fpl(2),fpl(3),fpl(4),fpl(5),fpl(6))
       tmp=Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)
       fpl(1:6)=fpl(1:6)*tmp
       fpl(9)=fpl(4)                                 ! Distribute independent
       fpl(4)=fpl(2) ; fpl(7)=fpl(3) ; fpl(8)=fpl(6) ! Symmetrise
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

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! Globalise Langevin random forces for shared RBs

       If (lshmv_rgd)Then
        Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxl,fyl,fzl,comm)
       EndIf

  ! Get strcom & vircom when starting afresh now done in w_calculate_forces
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

  ! store current integration variables

       cell0=cell
       vzero=volm
       eta0 =thermo%eta
       engke0=engke
       engrot0=engrot

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

       Do iter=1,mxiter

  ! integrate and apply Langevin thermostat - 1/4 step

          scale=Exp(-qstep*thermo%chi)
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
          engke=engke*scale**2
          engrot=engrot*scale**2

  ! constraint+pmf virial and stress

          vir=stat%vircon+stat%virpmf
          str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h1_scl barostat - 1/2 step
  ! augment str to include the random tensorial force on the barostat

          str1=str+fpl
          Call nst_h1_scl &
             (1,hstep,degfre,degrot,pmass,thermo%tai,volm,  &
             h_z,str1,stress,strcom,        &
             vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,strkin,strknf,strknt,engke,thermo,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

          scale=Exp(-qstep*thermo%chi)
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
          engke=engke*scale**2
          engrot=engrot*scale**2

  ! update velocity of FPs

          Do j=1,nfree
             i=lstfre(j)

             If (weight(i) > 1.0e-6_wp) Then
                tmp=hstep/weight(i)
                vxx(i)=vxx(i)+tmp*(fxx(i)+fxl(i))
                vyy(i)=vyy(i)+tmp*(fyy(i)+fyl(i))
                vzz(i)=vzz(i)+tmp*(fzz(i)+fzl(i))
             End If
          End Do

  ! scale cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

          aaa=tstep*thermo%eta
          Call mat_mul(aaa,aaa,bbb)
          aaa=uni1+aaa+0.5_wp*bbb
          Call mat_mul(aaa,cell0,cell)

  ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

          thermo%chi_p = thermo%eta(1)+thermo%eta(5)+thermo%eta(9)
          volm = volm*Exp(tstep*thermo%chi_p)
          thermo%chi_p = thermo%chi_p / 3.0_wp

  ! update position of FPs: second order taylor expansion of Exp(tstep*thermo%eta)

          Do j=1,nfree
             i=lstfre(j)

             If (weight(i) > 1.0e-6_wp) Then
                xxx(i)=tstep*vxx(i)+xxt(i)*aaa(1)+yyt(i)*aaa(2)+zzt(i)*aaa(3)
                yyy(i)=tstep*vyy(i)+xxt(i)*aaa(2)+yyt(i)*aaa(5)+zzt(i)*aaa(6)
                zzz(i)=tstep*vzz(i)+xxt(i)*aaa(3)+yyt(i)*aaa(6)+zzt(i)*aaa(9)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
           Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
          lstitr,&
          stat,pmf,cons,tmr,comm)
          End If

  ! restore original integration parameters as well as
  ! velocities if iter < mxiter
  ! in the next iteration stat%strcon and stat%strpmf are freshly new

          If (iter < mxiter) Then
             volm=vzero
             thermo%eta =eta0
             engke=engke0
             engrot=engrot0

             Do j=1,nfree
                i=lstfre(j)

                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do

             Do irgd=1,ntrgd
                rgdvxx(irgd) = rgdvxt(irgd)
                rgdvyy(irgd) = rgdvyt(irgd)
                rgdvzz(irgd) = rgdvzt(irgd)

                rgdoxx(irgd) = rgdoxt(irgd)
                rgdoyy(irgd) = rgdoyt(irgd)
                rgdozz(irgd) = rgdozt(irgd)
             End Do
          End If
       End Do

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

             fmxl=0.0_wp ; fmyl=0.0_wp ; fmzl=0.0_wp
             tqxl=0.0_wp ; tqyl=0.0_wp ; tqzl=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=indrgd(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rgdfrz(0,rgdtyp) == 0) Then
                   fmx=fmx+fxt(i)
                   fmy=fmy+fyt(i)
                   fmz=fmz+fzt(i)

                   fmxl=fmxl+fxl(i)
                   fmyl=fmyl+fyl(i)
                   fmzl=fmzl+fzl(i)
                End If

                tqx=tqx+ggy(krgd)*fzt(i)-ggz(krgd)*fyt(i)
                tqy=tqy+ggz(krgd)*fxt(i)-ggx(krgd)*fzt(i)
                tqz=tqz+ggx(krgd)*fyt(i)-ggy(krgd)*fxt(i)

                tqxl=tqxl+ggy(krgd)*fzl(i)-ggz(krgd)*fyl(i)
                tqyl=tqyl+ggz(krgd)*fxl(i)-ggx(krgd)*fzl(i)
                tqzl=tqzl+ggx(krgd)*fyl(i)-ggy(krgd)*fxl(i)
             End Do

  ! If the RB has 2+ frozen particles (ill=1) the net torque
  ! must align along the axis of rotation

             If (rgdfrz(0,rgdtyp) > 1) Then
                i1=indrgd(rgdind(1,rgdtyp),irgd)
                i2=indrgd(rgdind(2,rgdtyp),irgd)

                x(1)=xxt(i1)-xxt(i2)
                y(1)=yyt(i1)-yyt(i2)
                z(1)=zzt(i1)-zzt(i2)

                Call images(imcon,cell0,1,x,y,z)

                mxdr=1.0_wp/(x(1)**2+y(1)**2+z(1)**2)

                tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)*mxdr
                tqx=x(1)*tmp
                tqy=y(1)*tmp
                tqz=z(1)*tmp

                tmp=(x(1)*tqxl+y(1)*tqyl+z(1)*tqzl)*mxdr
                tqxl=x(1)*tmp
                tqyl=y(1)*tmp
                tqzl=z(1)*tmp
             End If

  ! current rotation matrix

             Call getrotmat(q0t(irgd),q1t(irgd),q2t(irgd),q3t(irgd),rot)

  ! calculate torque in principal frame

             trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
             try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
             trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

             trxl=tqxl*rot(1)+tqyl*rot(4)+tqzl*rot(7)
             tryl=tqxl*rot(2)+tqyl*rot(5)+tqzl*rot(8)
             trzl=tqxl*rot(3)+tqyl*rot(6)+tqzl*rot(9)

  ! calculate quaternion torques

             qt0=2.0_wp*(-q1t(irgd)*trx-q2t(irgd)*try-q3t(irgd)*trz)
             qt1=2.0_wp*( q0t(irgd)*trx-q3t(irgd)*try+q2t(irgd)*trz)
             qt2=2.0_wp*( q3t(irgd)*trx+q0t(irgd)*try-q1t(irgd)*trz)
             qt3=2.0_wp*(-q2t(irgd)*trx+q1t(irgd)*try+q0t(irgd)*trz)

             qt0l=2.0_wp*(-q1t(irgd)*trxl-q2t(irgd)*tryl-q3t(irgd)*trzl)
             qt1l=2.0_wp*( q0t(irgd)*trxl-q3t(irgd)*tryl+q2t(irgd)*trzl)
             qt2l=2.0_wp*( q3t(irgd)*trxl+q0t(irgd)*tryl-q1t(irgd)*trzl)
             qt3l=2.0_wp*(-q2t(irgd)*trxl+q1t(irgd)*tryl+q0t(irgd)*trzl)

  ! recover quaternion momenta

             opx=rgdoxx(irgd)*rgdrix(1,rgdtyp)
             opy=rgdoyy(irgd)*rgdriy(1,rgdtyp)
             opz=rgdozz(irgd)*rgdriz(1,rgdtyp)

             p0=2.0_wp*(-q1t(irgd)*opx-q2t(irgd)*opy-q3t(irgd)*opz)
             p1=2.0_wp*( q0t(irgd)*opx-q3t(irgd)*opy+q2t(irgd)*opz)
             p2=2.0_wp*( q3t(irgd)*opx+q0t(irgd)*opy-q1t(irgd)*opz)
             p3=2.0_wp*(-q2t(irgd)*opx+q1t(irgd)*opy+q0t(irgd)*opz)

  ! update quaternion momenta to half step

             p0=p0+hstep*(qt0+qt0l)
             p1=p1+hstep*(qt1+qt1l)
             p2=p2+hstep*(qt2+qt2l)
             p3=p3+hstep*(qt3+qt3l)

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
             rgdvxx(irgd)=rgdvxx(irgd)+tmp*(fmx+fmxl)
             rgdvyy(irgd)=rgdvyy(irgd)+tmp*(fmy+fmyl)
             rgdvzz(irgd)=rgdvzz(irgd)+tmp*(fmz+fmzl)

  ! update RB COM to full step

             rgdxxx(irgd) = tstep*rgdvxx(irgd) + rgdxxt(irgd)*aaa(1) + rgdyyt(irgd)*aaa(2) + rgdzzt(irgd)*aaa(3)
             rgdyyy(irgd) = tstep*rgdvyy(irgd) + rgdxxt(irgd)*aaa(2) + rgdyyt(irgd)*aaa(5) + rgdzzt(irgd)*aaa(6)
             rgdzzz(irgd) = tstep*rgdvzz(irgd) + rgdxxt(irgd)*aaa(3) + rgdyyt(irgd)*aaa(6) + rgdzzt(irgd)*aaa(9)

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
                         vxx(i)=xxt(i)*aaa(1)+yyt(i)*aaa(2)+zzt(i)*aaa(3)
                         vyy(i)=xxt(i)*aaa(2)+yyt(i)*aaa(5)+zzt(i)*aaa(6)
                         vzz(i)=xxt(i)*aaa(3)+yyt(i)*aaa(6)+zzt(i)*aaa(9)

                         x(1)=xxx(i)-vxx(i)
                         y(1)=yyy(i)-vyy(i)
                         z(1)=zzz(i)-vzz(i)
                         Call images(imcon,cell,1,x,y,z)
                         xxx(i)=x(1)+vxx(i)
                         yyy(i)=y(1)+vyy(i)
                         zzz(i)=z(1)+vzz(i)
                      End If

  ! new atomic velocities in lab frame

                      vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rgdvxx(irgd)
                      vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rgdvyy(irgd)
                      vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rgdvzz(irgd)
                   Else
                      x(1)=rgdxxx(irgd)-rgdxxt(irgd)
                      y(1)=rgdyyy(irgd)-rgdyyt(irgd)
                      z(1)=rgdzzz(irgd)-rgdzzt(irgd)
                      If (unsafe) Call images(imcon,cell,1,x,y,z) ! DD bound positions
                      xxx(i)=xxt(i)+x(1)
                      yyy(i)=yyt(i)+y(1)
                      zzz(i)=zzt(i)+z(1)
                   End If
                End If
             End Do

          Else

  ! update RB COM to full step

             rgdxxx(irgd) = rgdxxt(irgd)*aaa(1) + rgdyyt(irgd)*aaa(2) + rgdzzt(irgd)*aaa(3)
             rgdyyy(irgd) = rgdxxt(irgd)*aaa(2) + rgdyyt(irgd)*aaa(5) + rgdzzt(irgd)*aaa(6)
             rgdzzz(irgd) = rgdxxt(irgd)*aaa(3) + rgdyyt(irgd)*aaa(6) + rgdzzt(irgd)*aaa(9)

             Do jrgd=1,lrgd
                i=indrgd(jrgd,irgd) ! local index of particle/site

                If (i <= natms) Then
                   x(1)=rgdxxx(irgd)-rgdxxt(irgd)
                   y(1)=rgdyyy(irgd)-rgdyyt(irgd)
                   z(1)=rgdzzz(irgd)-rgdzzt(irgd)
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
 xxt,yyt,zzt,cshell%legshl,message,tmp,comm)) Then 
            Call info(message,.true.)

  ! scale Langevin random forces

             Do i=1,matms
                fxl(i) = fxl(i)*tmp
                fyl(i) = fyl(i)*tmp
                fzl(i) = fzl(i)*tmp
             End Do
             fpl = fpl*tmp

  ! restore initial conditions

             volm=vzero
             thermo%eta =eta0
             engke=engke0
             engrot=engrot0

             Do i=1,matms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)

                fxx(i) = fxt(i)
                fyy(i) = fyt(i)
                fzz(i) = fzt(i)
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

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       elrc=elrc0*tmp
       virlrc=virlrc0*tmp
       Do i=1,site%ntype_atom
          site%dens(i)=dens0(i)*tmp
       End Do

  ! get h_z for thermo%iso>1

       If (thermo%iso > 1) Then
          Call dcell(cell,celprp)
          h_z=celprp(9)
       End If

  ! second stage of velocity verlet algorithm

    Else

  ! Generate Langevin forces for particles and
  ! Langevin tensor force for barostat piston

       Call langevin_forces(nstep,temp,tstep,thermo%chi,fxl,fyl,fzl,cshell)
       If (lshmv_rgd)Then
         Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,fxl,fyl,fzl,comm)
       EndIf

       fpl=0.0_wp
       Call box_mueller_saru6(Int(degfre/3_li),nstep,fpl(1),fpl(2),fpl(3),fpl(4),fpl(5),fpl(6))
       tmp=Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)
       fpl(1:6)=fpl(1:6)*tmp
       fpl(9)=fpl(4)                                 ! Distribute independent
       fpl(4)=fpl(2) ; fpl(7)=fpl(3) ; fpl(8)=fpl(6) ! Symmetrise

  ! update velocity of FP

       Do j=1,nfree
          i=lstfre(j)

          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*(fxx(i)+fxl(i))
             vyy(i)=vyy(i)+tmp*(fyy(i)+fyl(i))
             vzz(i)=vzz(i)+tmp*(fzz(i)+fzl(i))
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_rattle(tstep,kit,&
                          pmf,cons,stat,tmr,comm)
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stre_s(strcom,ggx,ggy,ggz,fxx+fxl,fyy+fyl,fzz+fzl,comm)
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

             fmxl=0.0_wp ; fmyl=0.0_wp ; fmzl=0.0_wp
             tqxl=0.0_wp ; tqyl=0.0_wp ; tqzl=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=indrgd(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rgdfrz(0,rgdtyp) == 0) Then
                   fmx=fmx+fxx(i)
                   fmy=fmy+fyy(i)
                   fmz=fmz+fzz(i)

                   fmxl=fmxl+fxl(i)
                   fmyl=fmyl+fyl(i)
                   fmzl=fmzl+fzl(i)
                End If

                tqx=tqx+ggy(krgd)*fzz(i)-ggz(krgd)*fyy(i)
                tqy=tqy+ggz(krgd)*fxx(i)-ggx(krgd)*fzz(i)
                tqz=tqz+ggx(krgd)*fyy(i)-ggy(krgd)*fxx(i)

                tqxl=tqxl+ggy(krgd)*fzl(i)-ggz(krgd)*fyl(i)
                tqyl=tqyl+ggz(krgd)*fxl(i)-ggx(krgd)*fzl(i)
                tqzl=tqzl+ggx(krgd)*fyl(i)-ggy(krgd)*fxl(i)
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

                mxdr=1.0_wp/(x(1)**2+y(1)**2+z(1)**2)

                tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)*mxdr
                tqx=x(1)*tmp
                tqy=y(1)*tmp
                tqz=z(1)*tmp

                tmp=(x(1)*tqxl+y(1)*tqyl+z(1)*tqzl)*mxdr
                tqxl=x(1)*tmp
                tqyl=y(1)*tmp
                tqzl=z(1)*tmp
             End If

  ! current rotation matrix

             Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

  ! calculate torque in principal frame

             trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
             try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
             trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

             trxl=tqxl*rot(1)+tqyl*rot(4)+tqzl*rot(7)
             tryl=tqxl*rot(2)+tqyl*rot(5)+tqzl*rot(8)
             trzl=tqxl*rot(3)+tqyl*rot(6)+tqzl*rot(9)

  ! calculate quaternion torques

             qt0=2.0_wp*(-q1(irgd)*trx-q2(irgd)*try-q3(irgd)*trz)
             qt1=2.0_wp*( q0(irgd)*trx-q3(irgd)*try+q2(irgd)*trz)
             qt2=2.0_wp*( q3(irgd)*trx+q0(irgd)*try-q1(irgd)*trz)
             qt3=2.0_wp*(-q2(irgd)*trx+q1(irgd)*try+q0(irgd)*trz)

             qt0l=2.0_wp*(-q1(irgd)*trxl-q2(irgd)*tryl-q3(irgd)*trzl)
             qt1l=2.0_wp*( q0(irgd)*trxl-q3(irgd)*tryl+q2(irgd)*trzl)
             qt2l=2.0_wp*( q3(irgd)*trxl+q0(irgd)*tryl-q1(irgd)*trzl)
             qt3l=2.0_wp*(-q2(irgd)*trxl+q1(irgd)*tryl+q0(irgd)*trzl)

  ! recover quaternion momenta at half time step

             opx=rgdoxx(irgd)*rgdrix(1,rgdtyp)
             opy=rgdoyy(irgd)*rgdriy(1,rgdtyp)
             opz=rgdozz(irgd)*rgdriz(1,rgdtyp)

             p0=2.0_wp*(-q1(irgd)*opx-q2(irgd)*opy-q3(irgd)*opz)
             p1=2.0_wp*( q0(irgd)*opx-q3(irgd)*opy+q2(irgd)*opz)
             p2=2.0_wp*( q3(irgd)*opx+q0(irgd)*opy-q1(irgd)*opz)
             p3=2.0_wp*(-q2(irgd)*opx+q1(irgd)*opy+q0(irgd)*opz)

  ! update quaternion momenta to full step

             p0=p0+hstep*(qt0+qt0l)
             p1=p1+hstep*(qt1+qt1l)
             p2=p2+hstep*(qt2+qt2l)
             p3=p3+hstep*(qt3+qt3l)

  ! update RB angular & COM velocities to full step

             opx=0.5_wp*(-q1(irgd)*p0+q0(irgd)*p1+q3(irgd)*p2-q2(irgd)*p3)
             opy=0.5_wp*(-q2(irgd)*p0-q3(irgd)*p1+q0(irgd)*p2+q1(irgd)*p3)
             opz=0.5_wp*(-q3(irgd)*p0+q2(irgd)*p1-q1(irgd)*p2+q0(irgd)*p3)

             rgdoxx(irgd)=opx*rgdrix(2,rgdtyp)
             rgdoyy(irgd)=opy*rgdriy(2,rgdtyp)
             rgdozz(irgd)=opz*rgdriz(2,rgdtyp)

             tmp=hstep/rgdwgt(0,rgdtyp)
             rgdvxx(irgd)=rgdvxx(irgd)+tmp*(fmx+fmxl)
             rgdvyy(irgd)=rgdvyy(irgd)+tmp*(fmy+fmyl)
             rgdvzz(irgd)=rgdvzz(irgd)+tmp*(fmz+fmzl)

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

  ! update kinetic energy

       engknf=getknf(vxx,vyy,vzz,comm)
       engknt=getknt(rgdvxx,rgdvyy,rgdvzz,comm)

       engke=engknf+engknt

  ! update rotational energy

       engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
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
       engke=engke*scale**2
       engrot=engrot*scale**2

  ! constraint+pmf virial and stress

       vir=stat%vircon+stat%virpmf
       str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h1_scl barostat - 1/2 step
  ! augment str to include the random tensorial force on the barostat

       str1=str+fpl
       Call nst_h1_scl &
             (1,hstep,degfre,degrot,pmass,thermo%tai,volm,  &
             h_z,str1,stress,strcom,        &
             vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,strkin,strknf,strknt,engke,thermo,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
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
       engke=engke*scale**2
       engrot=engrot*scale**2

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

       tmp = ( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
         thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*pmass*tmp + thermo%press*volm

  ! remove system centre of mass velocity

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

  ! update kinetic energy and stress

       Call kinstresf(vxx,vyy,vzz,strknf,comm)
       Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

       strkin=strknf+strknt
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

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
       Write(message,'(a)') 'nst_l1 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nst_l1_vv
End Module nst_langevin
