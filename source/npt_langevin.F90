Module npt_langevin
  Use kinds,         Only : wp, li
  Use comms,         Only : comms_type,gmax
  Use domains,       Only : domains_type
  Use site, Only : site_type
  Use setup
  Use configuration, Only : imcon,cell,volm,natms,nlast,nfree,  &
                            lsi,lsa,lfrzn,lstfre,weight,        &
                            vxx,vyy,vzz
  Use particle,     Only : corePart
  Use rigid_bodies, Only : rigid_bodies_type,getrotmat,no_squish,rigid_bodies_stre_s
  Use langevin,      Only : fxl,fyl,fzl,fpl
  Use kinetics,      Only : getvom,getknf,getknt,getknr,getkin, &
                            kinstress,kinstresf,kinstrest
  Use shared_units,    Only : update_shared_units
  Use errors_warnings, Only : error,info
  Use numerics,        Only : box_mueller_saru1,images
  Use constraints,     Only : constraints_type,constraints_tags, apply_shake, apply_rattle
  Use pmf,             Only : pmf_tags,pmf_type
  Use npt_nose_hoover, Only : npt_h0_scl,npt_h0_scl,npt_h1_scl
  Use langevin,        Only : langevin_forces
  Use thermostat, Only : thermostat_type
  Use core_shell, Only : core_shell_type
  Use statistics, Only : stats_type
  Use timer, Only : timer_type
  Use thermostat, Only : adjust_timestep
  Use vdw, Only : vdw_type
  Implicit None

  Private

  Public :: npt_l0_vv, npt_l1_vv

Contains

  Subroutine npt_l0_vv(isw,lvar,mndis,mxdis,mxstp,tstep, &
             nstep,          &
             degfre,virtot,                     &
             consv,                             &
             strkin,engke,                      &
             cshell,cons,pmf,stat,thermo,sites, &
             vdws,domain,tmr,parts,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian equations of motion in
  ! molecular dynamics - velocity verlet with Langevin thermostat and
  ! barostat (isotropic pressure control) (symplectic)
  !
  ! isotropic cell fluctuations
  !
  ! reference: D. Quigley and M.I.J. Probert
  !            J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
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
    Real( Kind = wp ),  Intent( In    ) :: virtot
    Real( Kind = wp ),  Intent(   Out ) :: consv
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( corePart ),   Intent( InOut ) :: parts(:)
    Type( comms_type ), Intent( InOut ) :: comm

    Logical,           Save :: newjob = .true.
    Logical                 :: safe,lcol,lfst
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:9),iter,i
    Real( Kind = wp ), Save :: cell0(1:9),volm0,elrc0,virlrc0
    Real( Kind = wp ), Save :: temp,pmass
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: chip0,engke0
    Real( Kind = wp )       :: vzero
    Real( Kind = wp )       :: xt,yt,zt,vir,vir1,str(1:9),mxdr,tmp, &
                               scale,vom(1:3)


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
       Write(message,'(a)') 'npt_l  Use statistics, Only : stats_type0 allocation failure'
       Call error(0,message)
    End If


  ! timestep derivatives

    hstep = 0.5_wp*tstep
    qstep = 0.5_wp*hstep
    rstep = 1.0_wp/tstep

    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       cell0   = cell
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

  ! inertia parameter for barostat

       temp  = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
       pmass = (2.0_wp*thermo%sigma + 3.0_wp*boltz*temp) / (2.0_wp*pi*thermo%tai)**2

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+3
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake

  ! Langevin forces for particles are now generated in w_calculate_forces
  ! Generate Langevin pseudo-tensor force for barostat piston

       fpl=0.0_wp
       Call box_mueller_saru1(Int(degfre/3_li),nstep-1,tmp)
       tmp=tmp*Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)/3.0_wp
       fpl(1)=tmp
       fpl(5)=tmp
       fpl(9)=tmp
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Call constraints_tags(lstitr,cons,parts,comm)

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0) Call pmf_tags(lstitr,pmf,parts,comm)
    End If

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

  ! store current integration variables

       vzero=volm
       chip0=thermo%chi_p
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

  ! integrate and apply npt_h0_scl barostat - 1/2 step
  ! augment vir to include the random force on the barostat

          vir1=vir-3.0_wp*fpl(1)
          Call npt_h0_scl &
             (1,hstep,degfre,pmass,thermo%tai,volm,vir1,virtot, &
             vxx,vyy,vzz,engke,stat,thermo)

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
                vxx(i)=vxx(i)+tmp*(parts(i)%fxx+fxl(i))
                vyy(i)=vyy(i)+tmp*(parts(i)%fyy+fyl(i))
                vzz(i)=vzz(i)+tmp*(parts(i)%fzz+fzl(i))
             End If
          End Do

  ! update volume

          volm=volm*Exp(3.0_wp*tstep*thermo%chi_p)

  ! scale cell vectors - isotropic

          scale=(volm/volm0)**(1.0_wp/3.0_wp)
          cell=cell0*scale

  ! update positions

          scale=Exp(tstep*thermo%chi_p)
          Do i=1,natms
             If (weight(i) > 1.0e-6_wp) Then
                parts(i)%xxx=scale*xxt(i)+tstep*vxx(i)
                parts(i)%yyy=scale*yyt(i)+tstep*vyy(i)
                parts(i)%zzz=scale*zzt(i)+tstep*vzz(i)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
            Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
             lstitr,stat,pmf,cons,domain,tmr,parts,comm)
          End If

  ! restore original integration parameters as well as
  ! velocities if iter < mxiter
  ! in the next iteration stat%vircon and stat%virpmf are freshly new

          If (iter < mxiter) Then
             volm=vzero
             thermo%chi_p=chip0
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
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,parts,&
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
             thermo%chi_p=chip0
             engke=engke0

             Do i=1,natms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)

                parts(i)%fxx = fxt(i)
                parts(i)%fyy = fyt(i)
                parts(i)%fzz = fzt(i)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       vdws%elrc=elrc0*tmp
       vdws%vlrc=virlrc0*tmp
       Do i=1,sites%ntype_atom
          sites%dens(i)=dens0(i)*tmp
       End Do

  ! second stage of velocity verlet algorithm

    Else

  ! Generate Langevin forces for particles and
  ! Langevin pseudo-tensor force for barostat piston

       Call langevin_forces(nstep,temp,tstep,thermo%chi,fxl,fyl,fzl,cshell,parts)

       fpl=0.0_wp
       Call box_mueller_saru1(Int(degfre/3_li),nstep,tmp)
       tmp=tmp*Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)/3.0_wp
       fpl(1)=tmp
       fpl(5)=tmp
       fpl(9)=tmp

  ! update velocity

       Do i=1,natms
          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*(parts(i)%fxx+fxl(i))
             vyy(i)=vyy(i)+tmp*(parts(i)%fyy+fyl(i))
             vzz(i)=vzz(i)+tmp*(parts(i)%fzz+fzl(i))
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,kit,pmf,cons,stat,domain,tmr,comm)
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

  ! integrate and apply npt_h0_scl barostat - 1/2 step
  ! augment vir to include the random force on the barostat

       vir1=vir-3.0_wp*fpl(1)
       Call npt_h0_scl &
             (1,hstep,degfre,pmass,thermo%tai,volm,vir1,virtot, &
             vxx,vyy,vzz,engke,stat,thermo)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
       Do i=1,natms
          vxx(i)=scale*vxx(i)
          vyy(i)=scale*vyy(i)
          vzz(i)=scale*vzz(i)
       End Do
       engke=engke*scale**2

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*pmass*thermo%chi_p**2 + thermo%press*volm

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

  ! construct a 'mock' scaling tensor for xscale

    Do i=2,8
       thermo%eta(i)=0.0_wp
    End Do
    thermo%eta(1)=thermo%chi_p
    thermo%eta(5)=thermo%chi_p
    thermo%eta(9)=thermo%chi_p

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
       Write(message,'(a)') 'npt_l0 deallocation failure, node'
       Call error(0,message)
    End If
  End Subroutine npt_l0_vv

  Subroutine npt_l1_vv(isw,lvar,mndis,mxdis,mxstp,tstep, &
             nstep,          &
             degfre,degrot,virtot,              &
             consv,                             &
             strkin,strknf,strknt,engke,engrot, &
             strcom,vircom,                     &
             cshell,cons,pmf,stat,thermo,sites, &
             vdws,rigid,domain,tmr,parts,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
  ! RBs, equations of motion in molecular dynamics
  ! - velocity verlet with Langevin thermostat and
  ! barostat (isotropic pressure control) (symplectic)
  !
  ! isotropic cell fluctuations
  !
  ! reference: D. Quigley and M.I.J. Probert
  !            J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
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
    Real( Kind = wp ),  Intent( In    ) :: virtot
    Real( Kind = wp ),  Intent(   Out ) :: consv
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke, &
                                           strknf(1:9),strknt(1:9),engrot
    Real( Kind = wp ),  Intent( InOut ) :: strcom(1:9),vircom
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( corePart ),   Intent( InOut ) :: parts(:)
    Type( comms_type ), Intent( InOut ) :: comm


    Logical,           Save :: newjob = .true. , &
                               unsafe = .false.
    Logical                 :: safe,lcol,lfst
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:14),matms,iter,i,j,i1,i2, &
                               irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ), Save :: cell0(1:9),volm0,elrc0,virlrc0
    Real( Kind = wp ), Save :: temp,pmass
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: chip0,engke0,engrot0,engknf,engknt
    Real( Kind = wp )       :: czero(1:9),vzero
    Real( Kind = wp )       :: xt,yt,zt,vir,vir1,str(1:9),mxdr,tmp, &
                               scale,vom(1:3)
    Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                               opx,opy,opz,fmx,fmy,fmz,       &
                               tqx,tqy,tqz,trx,try,trz,       &
                               fmxl,fmyl,fmzl,                &
                               tqxl,tqyl,tqzl,trxl,tryl,trzl, &
                               qt0l,qt1l,qt2l,qt3l,           &
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

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len = 256 )        :: message
 

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
       Write(message,'(a)') 'npt_l1 allocation failure, node'
       Call error(0,message)
    End If


  ! timestep derivatives

    hstep = 0.5_wp*tstep
    qstep = 0.5_wp*hstep
    rstep = 1.0_wp/tstep

    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       cell0   = cell
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

  ! inertia parameter for barostat

       temp  = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
       pmass = (Real(degfre-degrot,wp) + 3.0_wp)*boltz*temp / (2.0_wp*pi*thermo%tai)**2

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+3
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake

  ! unsafe positioning due to possibly locally shared RBs

       unsafe=(Any(domain%map == comm%idnode))

  ! Langevin forces for particles are now generated in w_calculate_forces
  ! Generate Langevin pseudo-tensor force for barostat piston

       fpl=0.0_wp
       Call box_mueller_saru1(Int(degfre/3_li),nstep-1,tmp)
       tmp=tmp*Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)/3.0_wp
       fpl(1)=tmp
       fpl(5)=tmp
       fpl(9)=tmp
    End If

  ! set matms

    matms=nlast
    If (comm%mxnode == 1) matms=natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Call constraints_tags(lstitr,cons,parts,comm)

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0) Call pmf_tags(lstitr,pmf,parts,comm)
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

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! Globalise Langevin random forces for shared RBs

       If (rigid%share)Then
         Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared, &
           rigid%map_shared,fxl,fyl,fzl,domain,comm)
       EndIf

  ! Get strcom & vircom when starting afresh now done in w_calculate_forces
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

  ! store current integration variables

       czero=cell
       vzero=volm
       chip0=thermo%chi_p
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

          Do irgd=1,rigid%n_types
             rigid%vxx(irgd)=scale*rigid%vxx(irgd)
             rigid%vyy(irgd)=scale*rigid%vyy(irgd)
             rigid%vzz(irgd)=scale*rigid%vzz(irgd)

             rigid%oxx(irgd)=scale*rigid%oxx(irgd)
             rigid%oyy(irgd)=scale*rigid%oyy(irgd)
             rigid%ozz(irgd)=scale*rigid%ozz(irgd)
          End Do
          engke=engke*scale**2
          engrot=engrot*scale**2

  ! constraint+pmf virial and stress

          vir=stat%vircon+stat%virpmf
          str=stat%strcon+stat%strpmf

  ! integrate and apply npt_h1_scl barostat - 1/2 step
  ! augment vir to include the random force on the barostat

          vir1=vir-3.0_wp*fpl(1)
          Call npt_h1_scl &
             (1,hstep,degfre,degrot,pmass,thermo%tai,volm,vir1,virtot,vircom, &
             vxx,vyy,vzz,engke,stat,rigid,thermo)

  ! integrate and apply Langevin thermostat - 1/4 step

          scale=Exp(-qstep*thermo%chi)
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
          engke=engke*scale**2
          engrot=engrot*scale**2

  ! update velocity of FPs

          Do j=1,nfree
             i=lstfre(j)

             If (weight(i) > 1.0e-6_wp) Then
                tmp=hstep/weight(i)
                vxx(i)=vxx(i)+tmp*(parts(i)%fxx+fxl(i))
                vyy(i)=vyy(i)+tmp*(parts(i)%fyy+fyl(i))
                vzz(i)=vzz(i)+tmp*(parts(i)%fzz+fzl(i))
             End If
          End Do

  ! update volume

          volm=volm*Exp(3.0_wp*tstep*thermo%chi_p)

  ! scale cell vectors - isotropic

          scale=(volm/volm0)**(1.0_wp/3.0_wp)
          cell=cell0*scale

  ! update position of FPs

          scale=Exp(tstep*thermo%chi_p)
          Do j=1,nfree
             i=lstfre(j)

             If (weight(i) > 1.0e-6_wp) Then
                parts(i)%xxx=scale*xxt(i)+tstep*vxx(i)
                parts(i)%yyy=scale*yyt(i)+tstep*vyy(i)
                parts(i)%zzz=scale*zzt(i)+tstep*vzz(i)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
            Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
             lstitr,stat,pmf,cons,domain,tmr,parts,comm)
          End If

  ! restore original integration parameters as well as
  ! velocities if iter < mxiter
  ! in the next iteration stat%vircon and stat%virpmf are freshly new

          If (iter < mxiter) Then
             volm=vzero
             thermo%chi_p=chip0
             engke=engke0
             engrot=engrot0

             Do j=1,nfree
                i=lstfre(j)

                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do

             Do irgd=1,rigid%n_types
                rigid%vxx(irgd) = rgdvxt(irgd)
                rigid%vyy(irgd) = rgdvyt(irgd)
                rigid%vzz(irgd) = rgdvzt(irgd)

                rigid%oxx(irgd) = rgdoxt(irgd)
                rigid%oyy(irgd) = rgdoyt(irgd)
                rigid%ozz(irgd) = rgdozt(irgd)
             End Do
          End If
       End Do

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

             fmxl=0.0_wp ; fmyl=0.0_wp ; fmzl=0.0_wp
             tqxl=0.0_wp ; tqyl=0.0_wp ; tqzl=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rigid%frozen(0,rgdtyp) == 0) Then
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

             If (rigid%frozen(0,rgdtyp) > 1) Then
                i1=rigid%index_local(rigid%index_global(1,rgdtyp),irgd)
                i2=rigid%index_local(rigid%index_global(2,rgdtyp),irgd)

                x(1)=xxt(i1)-xxt(i2)
                y(1)=yyt(i1)-yyt(i2)
                z(1)=zzt(i1)-zzt(i2)

                Call images(imcon,czero,1,x,y,z)

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

             opx=rigid%oxx(irgd)*rigid%rix(1,rgdtyp)
             opy=rigid%oyy(irgd)*rigid%riy(1,rgdtyp)
             opz=rigid%ozz(irgd)*rigid%riz(1,rgdtyp)

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
             rigid%vxx(irgd)=rigid%vxx(irgd)+tmp*(fmx+fmxl)
             rigid%vyy(irgd)=rigid%vyy(irgd)+tmp*(fmy+fmyl)
             rigid%vzz(irgd)=rigid%vzz(irgd)+tmp*(fmz+fmzl)

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

                      parts(i)%xxx=rot(1)*x(1)+rot(2)*y(1)+rot(3)*z(1) + rigid%xxx(irgd)
                      parts(i)%yyy=rot(4)*x(1)+rot(5)*y(1)+rot(6)*z(1) + rigid%yyy(irgd)
                      parts(i)%zzz=rot(7)*x(1)+rot(8)*y(1)+rot(9)*z(1) + rigid%zzz(irgd)

  ! new atomic velocities in body frame

                      vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                      vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                      vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

  ! DD bound positions

                      If (unsafe) Then
                         vxx(i)=scale*xxt(i)
                         vyy(i)=scale*yyt(i)
                         vzz(i)=scale*zzt(i)

                         x(1)=parts(i)%xxx-vxx(i)
                         y(1)=parts(i)%yyy-vyy(i)
                         z(1)=parts(i)%zzz-vzz(i)
                         Call images(imcon,cell,1,x,y,z)
                         parts(i)%xxx=x(1)+vxx(i)
                         parts(i)%yyy=y(1)+vyy(i)
                         parts(i)%zzz=z(1)+vzz(i)
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
                      parts(i)%xxx=xxt(i)+x(1)
                      parts(i)%yyy=yyt(i)+y(1)
                      parts(i)%zzz=zzt(i)+z(1)
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
                   parts(i)%xxx=xxt(i)+x(1)
                   parts(i)%yyy=yyt(i)+y(1)
                   parts(i)%zzz=zzt(i)+z(1)
                End If
             End Do

          End If
       End Do

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,natms,parts,&
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
             thermo%chi_p=chip0
             engke=engke0
             engrot=engrot0

             Do i=1,matms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)

                parts(i)%fxx = fxt(i)
                parts(i)%fyy = fyt(i)
                parts(i)%fzz = fzt(i)
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

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       vdws%elrc=elrc0*tmp
       vdws%vlrc=virlrc0*tmp
       Do i=1,sites%ntype_atom
          sites%dens(i)=dens0(i)*tmp
       End Do

  ! second stage of velocity verlet algorithm

    Else

  ! Generate Langevin forces for particles and
  ! Langevin pseudo-tensor force for barostat piston

       Call langevin_forces(nstep,temp,tstep,thermo%chi,fxl,fyl,fzl,cshell,parts)
       If (rigid%share)Then
         Call update_shared_units(natms,nlast,lsi,lsa,rigid%list_shared, &
           rigid%map_shared,fxl,fyl,fzl,domain,comm)
       EndIf

       fpl=0.0_wp
       Call box_mueller_saru1(Int(degfre/3_li),nstep,tmp)
       tmp=tmp*Sqrt(2.0_wp*thermo%tai*boltz*temp*pmass*rstep)/3.0_wp
       fpl(1)=tmp
       fpl(5)=tmp
       fpl(9)=tmp

  ! update velocity of FPs

       Do j=1,nfree
          i=lstfre(j)

          If (weight(i) > 1.0e-6_wp) Then
             tmp=hstep/weight(i)
             vxx(i)=vxx(i)+tmp*(parts(i)%fxx+fxl(i))
             vyy(i)=vyy(i)+tmp*(parts(i)%fyy+fyl(i))
             vzz(i)=vzz(i)+tmp*(parts(i)%fzz+fzl(i))
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,kit,pmf,cons,stat,domain,tmr,comm)
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stre_s(strcom,ggx,ggy,ggz,parts,rigid,comm,fxl,fyl,fzl)
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

             fmxl=0.0_wp ; fmyl=0.0_wp ; fmzl=0.0_wp
             tqxl=0.0_wp ; tqyl=0.0_wp ; tqzl=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rigid%frozen(0,rgdtyp) == 0) Then
                   fmx=fmx+parts(i)%fxx
                   fmy=fmy+parts(i)%fyy
                   fmz=fmz+parts(i)%fzz

                   fmxl=fmxl+fxl(i)
                   fmyl=fmyl+fyl(i)
                   fmzl=fmzl+fzl(i)
                End If

                tqx=tqx+ggy(krgd)*parts(i)%fzz-ggz(krgd)*parts(i)%fyy
                tqy=tqy+ggz(krgd)*parts(i)%fxx-ggx(krgd)*parts(i)%fzz
                tqz=tqz+ggx(krgd)*parts(i)%fyy-ggy(krgd)*parts(i)%fxx

                tqxl=tqxl+ggy(krgd)*fzl(i)-ggz(krgd)*fyl(i)
                tqyl=tqyl+ggz(krgd)*fxl(i)-ggx(krgd)*fzl(i)
                tqzl=tqzl+ggx(krgd)*fyl(i)-ggy(krgd)*fxl(i)
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

             Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

  ! calculate torque in principal frame

             trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
             try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
             trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

             trxl=tqxl*rot(1)+tqyl*rot(4)+tqzl*rot(7)
             tryl=tqxl*rot(2)+tqyl*rot(5)+tqzl*rot(8)
             trzl=tqxl*rot(3)+tqyl*rot(6)+tqzl*rot(9)

  ! calculate quaternion torques

             qt0=2.0_wp*(-rigid%q1(irgd)*trx-rigid%q2(irgd)*try-rigid%q3(irgd)*trz)
             qt1=2.0_wp*( rigid%q0(irgd)*trx-rigid%q3(irgd)*try+rigid%q2(irgd)*trz)
             qt2=2.0_wp*( rigid%q3(irgd)*trx+rigid%q0(irgd)*try-rigid%q1(irgd)*trz)
             qt3=2.0_wp*(-rigid%q2(irgd)*trx+rigid%q1(irgd)*try+rigid%q0(irgd)*trz)

             qt0l=2.0_wp*(-rigid%q1(irgd)*trxl-rigid%q2(irgd)*tryl-rigid%q3(irgd)*trzl)
             qt1l=2.0_wp*( rigid%q0(irgd)*trxl-rigid%q3(irgd)*tryl+rigid%q2(irgd)*trzl)
             qt2l=2.0_wp*( rigid%q3(irgd)*trxl+rigid%q0(irgd)*tryl-rigid%q1(irgd)*trzl)
             qt3l=2.0_wp*(-rigid%q2(irgd)*trxl+rigid%q1(irgd)*tryl+rigid%q0(irgd)*trzl)

  ! recover quaternion momenta at half time step

             opx=rigid%oxx(irgd)*rigid%rix(1,rgdtyp)
             opy=rigid%oyy(irgd)*rigid%riy(1,rgdtyp)
             opz=rigid%ozz(irgd)*rigid%riz(1,rgdtyp)

             p0=2.0_wp*(-rigid%q1(irgd)*opx-rigid%q2(irgd)*opy-rigid%q3(irgd)*opz)
             p1=2.0_wp*( rigid%q0(irgd)*opx-rigid%q3(irgd)*opy+rigid%q2(irgd)*opz)
             p2=2.0_wp*( rigid%q3(irgd)*opx+rigid%q0(irgd)*opy-rigid%q1(irgd)*opz)
             p3=2.0_wp*(-rigid%q2(irgd)*opx+rigid%q1(irgd)*opy+rigid%q0(irgd)*opz)

  ! update quaternion momenta to full step

             p0=p0+hstep*(qt0+qt0l)
             p1=p1+hstep*(qt1+qt1l)
             p2=p2+hstep*(qt2+qt2l)
             p3=p3+hstep*(qt3+qt3l)

  ! update RB angular & COM velocities to full step

             opx=0.5_wp*(-rigid%q1(irgd)*p0+rigid%q0(irgd)*p1+rigid%q3(irgd)*p2-rigid%q2(irgd)*p3)
             opy=0.5_wp*(-rigid%q2(irgd)*p0-rigid%q3(irgd)*p1+rigid%q0(irgd)*p2+rigid%q1(irgd)*p3)
             opz=0.5_wp*(-rigid%q3(irgd)*p0+rigid%q2(irgd)*p1-rigid%q1(irgd)*p2+rigid%q0(irgd)*p3)

             rigid%oxx(irgd)=opx*rigid%rix(2,rgdtyp)
             rigid%oyy(irgd)=opy*rigid%riy(2,rgdtyp)
             rigid%ozz(irgd)=opz*rigid%riz(2,rgdtyp)

             tmp=hstep/rigid%weight(0,rgdtyp)
             rigid%vxx(irgd)=rigid%vxx(irgd)+tmp*(fmx+fmxl)
             rigid%vyy(irgd)=rigid%vyy(irgd)+tmp*(fmy+fmyl)
             rigid%vzz(irgd)=rigid%vzz(irgd)+tmp*(fmz+fmzl)

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

  ! update kinetic energy

       engknf=getknf(vxx,vyy,vzz,comm)
       engknt=getknt(rigid,comm)

       engke=engknf+engknt

  ! update rotational energy

       engrot=getknr(rigid,comm)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
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
       engke=engke*scale**2
       engrot=engrot*scale**2

  ! constraint+pmf virial and stress

       vir=stat%vircon+stat%virpmf
       str=stat%strcon+stat%strpmf

  ! integrate and apply npt_h1_scl barostat - 1/2 step
  ! augment vir to include the random force on the barostat

       vir1=vir-3.0_wp*fpl(1)
       Call npt_h1_scl &
             (1,hstep,degfre,degrot,pmass,thermo%tai,volm,vir1,virtot,vircom, &
             vxx,vyy,vzz,engke,stat,rigid,thermo)

  ! integrate and apply Langevin thermostat - 1/4 step

       scale=Exp(-qstep*thermo%chi)
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
       engke=engke*scale**2
       engrot=engrot*scale**2

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*pmass*thermo%chi_p**2 + thermo%press*volm

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

       Call kinstresf(vxx,vyy,vzz,strknf,comm)
       Call kinstrest(rigid,strknt,comm)

       strkin=strknf+strknt
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

    End If

  ! construct a 'mock' scaling tensor for xscale

    Do i=2,8
       thermo%eta(i)=0.0_wp
    End Do
    thermo%eta(1)=thermo%chi_p
    thermo%eta(5)=thermo%chi_p
    thermo%eta(9)=thermo%chi_p

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
       Write(message,'(a)') 'npt_l1 deallocation failure'
       Call error(0,message)
    End If
  End Subroutine npt_l1_vv
End Module npt_langevin
