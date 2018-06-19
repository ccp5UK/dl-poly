Module minimise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring minimisation routine arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use comms,           Only : comms_type,gsum,gmax
  Use setup,           Only : engunit,nrite,output, &
                              mxatms,zero_plus, &
                              mxrgd,mxlrgd
  Use configuration,   Only : natms,nlast,nfree,          &
                              lsi,lsa,lfrzn,lfree,lstfre, &
                              weight,xxx,yyy,zzz,fxx,fyy,fzz, &
                              imcon,cell,vxx,vyy,vzz, &
                              write_config,getcom
  Use rigid_bodies,    Only : lshmv_rgd,lishp_rgd,lashp_rgd,rgdxxx,rgdyyy,rgdzzz,&
                              listrgd,rgdoxx,rgdoyy,rgdozz,rgdvxx,rgdvyy,rgdvzz,&
                              rgdx,rgdy,rgdz,rgdfrz,rgdrix,indrgd,rgdwgt,rgdmeg,&
                              rgdriz,rgdriy,ntrgd,q0,q1,q2,q3,rgdind,&
                              q_setup,getrotmat
  Use parse,           Only : strip_blanks,lower_case

  Use kinetics,        Only : getvom,getkin,getknf,getknt,getknr, &
                              kinstress,kinstresf,kinstrest
  Use numerics,        Only : images,invert
  Use pmf,             Only : pmf_tags,pmf_pseudo_bonds,pmf_type
  Use shared_units,    Only : update_shared_units
  Use rigid_bodies,    Only : rigid_bodies_split_torque,rigid_bodies_move
  Use errors_warnings, Only : error,warning,info
  Use statistics, Only : stats_type
  Use constraints, Only : constraints_type,constraints_tags,constraints_pseudo_bonds
  Use netcdf_wrap, Only : netcdf_param

  Implicit None

  Logical,                        Save :: l_x = .false.

  Real( Kind = wp ),              Save :: passmin(1:5) = (/ &
                                          0.0_wp         ,  & ! cycles counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average cycles
                                          999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
                                          0.0_wp /)           ! maximum cycles


  Real( Kind = wp ), Allocatable, Save :: oxx(:),oyy(:),ozz(:)

  Public :: allocate_minimise_arrays,deallocate_minimise_arrays, minimise_relax


Contains

  Subroutine allocate_minimise_arrays()

    Integer :: fail

    fail = 0

    Allocate (oxx(1:mxatms),oyy(1:mxatms),ozz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1038)

    oxx = 0.0_wp ; oyy = 0.0_wp ; ozz = 0.0_wp

  End Subroutine allocate_minimise_arrays

  Subroutine deallocate_minimise_arrays()

    Integer :: fail

    fail = 0

    Deallocate (oxx,oyy,ozz, Stat = fail)

    If (fail > 0) Call error(1039)

  End Subroutine deallocate_minimise_arrays

  Subroutine minimise_relax &
           (l_str,relaxed,rdf_collect,megatm,megpmf,megrgd, &
           keymin,min_tol,tstep,stpcfg,stat,pmf,cons,netcdf,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for minimising the structure by using conjugate
! gradient method (CGM).
!
! Note: minimisation type and criterion:
!       keymin=0 : absolute force
!       keymin=1 : relative energy
!       keymin=2 : absolute displacement
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith february 2014
! contrib   - a.m.elena february 2017
! contrib   - i.t.todorov february 2017
! contrib   - i.scivetti april 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: l_str
  Logical,           Intent( InOut ) :: relaxed,rdf_collect
  Integer,           Intent( In    ) :: megatm, &
                                        megpmf,megrgd,keymin
  Real( Kind = wp ), Intent( In    ) :: min_tol(1:2),tstep,stpcfg
  Type( stats_type ), Intent(InOut) :: stat
  Type( pmf_type ), Intent( InOut ) :: pmf
  Type( constraints_type ), Intent(InOut) :: cons
  Type( netcdf_param ), Intent( In    ) :: netcdf
  Type( comms_type ), Intent( inOut ) :: comm

  Logical,              Save :: newjob = .true. , l_rdf, l_mov
  Character( Len = 8 ), Save :: word
  Character( Len = 6 )       :: name
  Integer,              Save :: keyopt
  Integer                    :: fail(1:8),i,j,levcfg
  Real( Kind = wp ),    Save :: total,grad_tol,eng_tol,dist_tol,step,           &
                                eng_0,eng_min,eng,eng0,eng1,eng2, &
                                grad,grad0,grad1,grad2,onorm,sgn,stride,gamma

! OUTPUT existence

  Logical               :: l_out
  Character( Len = 10 ) :: c_out

! Optimisation iteration and convergence limits

  Integer, Parameter      :: mxpass = 1000
  Real( Kind = wp ), Save :: min_pass

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

! Constraints and PMFs arrays

  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: txx(:),tyy(:),tzz(:)
  Real( Kind = wp ), Allocatable :: uxx(:),uyy(:),uzz(:)
  Character( Len = 256 ) :: message


  fail=0
  If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
     Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
     call cons%allocate_work(mxatms)
     Call pmf%allocate_work()
  End If
  If (megrgd > 0) Then
     Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),         Stat=fail(6))
     Allocate (uxx(1:mxatms),uyy(1:mxatms),uzz(1:mxatms),         Stat=fail(7))
  End If
  Allocate (gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms),            Stat=fail(8))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'minimise_relax allocation failure'
     Call error(0,message)
  End If

  If (newjob) Then
     newjob = .false.

! At start no optimisation has been attempted yet

     keyopt = 0

! At start the minimum energy is defined as zero

     eng_min = Huge(0.0_wp)

! Passage accumulators are initialised in minimise_module
! passmin(1) - cycles counter
! passmin(2) - access counter
! passmin(3) - average cycles
! passmin(4) - minimum cycles
! passmin(5) - maximum cycles

! total number of active particles (excluding frozen sites and massless shells)

     total=0.0_wp
     Do i=1,natms
        If (lfrzn(i) == 0 .and. (weight(i) > 1.0e-6_wp .or. lfree(i) == 1)) &
           total=total+1.0_wp
     End Do
     Call gsum(comm,total)
  End If

! Step length for relaxation

  If (min_tol(2) > zero_plus) Then

! Optionally specified

     step=min_tol(2)

  Else

! default if unspecified

     step=tstep**2

! enlarged depending on functionality if defaulted

     If (cons%megcon == 0 .and. megpmf == 0) Then
        If (megrgd == 0) Then
           step=10.0_wp*step
        Else
           step=5.0_wp*step
        End If
     End If

  End If

  If (keyopt == 0) Then

! Initial configuration energy

     eng_0=stpcfg

! Allocate working arrays

     Call allocate_minimise_arrays()

! No minimisation is yet attempted

     relaxed=.false.

! No RB move is yet attempted

     l_mov=.false.

! Minimum needed for a pass for this minimisation cycle

     min_pass = Huge(1.0_wp)

! Avoid rdf calculation redundancy

     l_rdf=rdf_collect
     If (rdf_collect) rdf_collect=.false.

! Determine optimisation

     If      (keymin == 0) Then
        word='force   '
     Else If (keymin == 1) Then
        word='energy  '
     Else If (keymin == 2) Then
        word='distance'
     End If

! Print header

    If (l_str) Then
      Write(message,'(3(1x,a),6x,a,10x,a,10x,a,11x,a,5x,a,1p,e11.4,3x,a,e11.4)') &
        'Minimising',word,'pass','eng_tot','grad_tol','eng_tol','dist_tol','tol=', min_tol(1),'step=',step
      Call info(message,.true.)
      Write(message,"(1x,130('-'))")
      Call info(message,.true.)
    End If
  End If

! Load original forces

  Do i=1,natms
     gxx(i)=fxx(i)
     gyy(i)=fyy(i)
     gzz(i)=fzz(i)
  End Do

! Minimised energy is current configuration energy

  eng=stpcfg

! Calculate pseudo forces and energy for constraint bonds and PMFs

  If (cons%megcon > 0 .or. megpmf > 0) Then
     lstitr(1:natms)=.false. ! initialise lstitr

     If (cons%megcon > 0) Then
        Call constraints_tags(lstitr,cons,comm)
        Call constraints_pseudo_bonds(gxx,gyy,gzz,stat,cons,comm)
        eng=eng+stat%engcon
     End If

     If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr,pmf,comm)
        Call pmf_pseudo_bonds(gxx,gyy,gzz,stat,pmf,comm)
        eng=eng+stat%engpmf
     End If
  End If

! Average forces over all members of a RB and split torques accordingly

  If (megrgd > 0) Then
     If (lshmv_rgd) Then
       Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,gxx,gyy,gzz,comm)
     End If
     Call rigid_bodies_split_torque(gxx,gyy,gzz,txx,tyy,tzz,uxx,uyy,uzz,comm)
  End If

! Initialise/get eng_tol & verify relaxed condition

  eng_tol=0.0_wp
  If (keyopt > 0) Then
     eng_tol=Abs(1.0_wp-eng2/eng)
     If (keymin == 1) relaxed=(eng_tol < min_tol(1))
  End If

! Current gradient (modulus of the total force)
! massless shells and frozen particles have zero forces!

  grad=0.0_wp
  Do i=1,natms
     grad=grad+gxx(i)**2+gyy(i)**2+gzz(i)**2
  End Do
  Call gsum(comm,grad)
  grad=Sqrt(grad)

! Get grad_tol & verify relaxed condition

  grad_tol=grad/total
  If (keymin == 0) relaxed=(grad_tol < min_tol(1))

! Initialise dist_tol

  dist_tol=0.0_wp

! CHECK FOR CONVERGENCE

  If (.not.relaxed) Then

! Increment main passage counter

     passmin(1)=passmin(1)+1.0_wp

! min_pass = Min(min_pass,._tol)

     If      (keymin == 0) Then
        min_pass = Min(min_pass,grad_tol)
     Else If (keymin == 1) Then
        If (keyopt > 0) min_pass = Min(min_pass,eng_tol)
     Else If (keymin == 2) Then
        min_pass = Min(min_pass,dist_tol)
     End If

! If in mxpass iterations we are not there, give up but
! allow for thermo%tension-fold boost in iteration cycle length
! for the very first MD step

     If (Nint(passmin(2)) == 0) Then
        If (Nint(passmin(1)) >= 10*mxpass) Then
           Call warning(330,min_tol(1),min_pass,0.0_wp)
           Call error(474)
        End If
     Else
        If (Nint(passmin(1)) >= mxpass) Then
           Call warning(330,min_tol(1),min_pass,0.0_wp)
           Call error(474)
        End If
     End If

  Else

     Go To 100

  End If

  If      (keyopt == 0) Then

! Original configuration energy

     eng0=eng
     eng1=eng
     eng2=eng

! Original gradient (modulus of the total force)

     onorm=grad
     grad0=grad
     grad2=grad

! Set original search direction

     Do i=1,natms
        oxx(i)=gxx(i)
        oyy(i)=gyy(i)
        ozz(i)=gzz(i)
     End Do

     keyopt=1
     sgn=1.0_wp
     stride=sgn*step

  Else If (keyopt == 1) Then

! Line search along chosen direction

     eng1=eng0
     eng2=eng

     grad1=grad2
     grad2=0.0_wp
     Do i=1,natms
        grad2=grad2+oxx(i)*gxx(i)+oyy(i)*gyy(i)+ozz(i)*gzz(i)
     End Do
     Call gsum(comm,grad2)
     grad2=sgn*grad2/onorm

! Linear extrapolation to minimum

     If (grad2 < 0.0_wp) Then ! BACK UP FROM THIS DIRECTION
        keyopt=2
        stride=sgn*step*grad2/(grad1-grad2)
     Else                     ! CARRY ON IN THIS DIRECTION
        stride=sgn*step
     End If

  Else If (keyopt == 2) Then

! Construct conjugate search vector

     eng1=eng2
     eng2=eng

     gamma=(grad/grad0)**2
     grad0=grad
     grad2=0.0_wp
     onorm=0.0_wp
     Do i=1,natms
        oxx(i)=gxx(i)+gamma*oxx(i)
        oyy(i)=gyy(i)+gamma*oyy(i)
        ozz(i)=gzz(i)+gamma*ozz(i)

        onorm=onorm+oxx(i)**2+oyy(i)**2+ozz(i)**2
        grad2=grad2+oxx(i)*gxx(i)+oyy(i)*gyy(i)+ozz(i)*gzz(i)
     End Do
     Call gsum(comm,onorm)
     onorm=Sqrt(onorm)
     Call gsum(comm,grad2)
     grad2=grad2/onorm
     sgn=Sign(1.0_wp,grad2)
     grad2=sgn*grad2

     keyopt=1
     stride=sgn*step

  End If

! Move particles to their new positions accordingly

  If (megrgd > 0) Then

! active free particles

     Do j=1,nfree
        i=lstfre(j)

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           xxx(i)=xxx(i)+stride*oxx(i)
           yyy(i)=yyy(i)+stride*oyy(i)
           zzz(i)=zzz(i)+stride*ozz(i)
           dist_tol=Max(dist_tol,oxx(i)**2+oyy(i)**2+ozz(i)**2)
        End If
     End Do
     dist_tol=Sqrt(dist_tol)*Abs(stride)

! RB particles

     Call rigid_bodies_move(stride,oxx,oyy,ozz,txx,tyy,tzz,uxx,uyy,uzz,dist_tol)
     l_mov=.true.

  Else

! active particles

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           xxx(i)=xxx(i)+stride*oxx(i)
           yyy(i)=yyy(i)+stride*oyy(i)
           zzz(i)=zzz(i)+stride*ozz(i)
           dist_tol=Max(dist_tol,oxx(i)**2+oyy(i)**2+ozz(i)**2)
        End If
     End Do
     dist_tol=Sqrt(dist_tol)*Abs(stride)

  End If
  Call gmax(comm,dist_tol)

  If (keymin == 2) relaxed=(dist_tol < min_tol(1))

! Fit headers in and Close and Open OUTPUT at every 25th print-out

  i=Nint(passmin(1))
  If (l_str) Then
    Write(message,'(1x,i23,1p,4e18.8)') i-1,eng/engunit,grad_tol,eng_tol,dist_tol
    Call info(message,.true.)
    If (Mod(i,25) == 0) Then
      Write(message,"(1x,130('-'))")
      Call info(message,.true.)
      Write(message,'(3(1x,a),6x,a,10x,a,10x,a,11x,a,5x,a,1p,e11.4,3x,a,e11.4)') &
        'Minimising',word,'pass','eng_tot','grad_tol','eng_tol','dist_tol','tol=', min_tol(1),'step=',step
      Call info(message,.true.)
      Write(message,"(1x,130('-'))")
      Call info(message,.true.)

      If (comm%idnode == 0) Then
        Inquire(File=Trim(output), Exist=l_out, Position=c_out)
        Call strip_blanks(c_out)
        Call lower_case(c_out)
        If (l_out .and. c_out(1:6) == 'append') Then
          Close(Unit=nrite)
          Open(Unit=nrite, File=Trim(output), Position='append')
        End If
      End If
    End If
  End If

100 Continue

  l_x=(.not.relaxed) ! Transportation flag
  If (relaxed) Then

! Final/Only printout

     i=Nint(passmin(1))
     If (.not.l_str) Then
       Write(message,'(3(1x,a),5x,a,10x,a,10x,a,11x,a,5x,a,1p,e11.4,3x,a,e11.4)') &
         'Minimised',word,'passes','eng_tot','grad_tol','eng_tol','dist_tol','tol=', min_tol(1),'step=',step
       Call info(message,.true.)
       Write(message,"(1x,130('-'))")
       Call info(message,.true.)
     End If
     Write(message,'(1x,i23,1p,4e18.8)') i,eng/engunit,grad_tol,eng_tol,dist_tol
     Call info(message,.true.)
     Write(message,"(1x,130('-'))")
     Call info(message,.true.)

! Collect passage statistics

     passmin(3)=passmin(2)*passmin(3)
     passmin(2)=passmin(2)+1.0_wp
     passmin(3)=passmin(3)/passmin(2)+passmin(1)/passmin(2)
     passmin(4)=Min(passmin(1),passmin(4))
     passmin(5)=Max(passmin(1),passmin(5))

! Rewind keyopt and main passage counter

     keyopt =0
     passmin(1)=0.0_wp

! Resume rdf%rdf calculations

     If (l_rdf) rdf_collect=l_rdf

! Deallocate working arrays

     Call deallocate_minimise_arrays()

! Dump the lowest energy configuration

     If (eng < eng_min) Then
        eng_min=eng

        name = 'CFGMIN' ! file name
        levcfg = 0      ! define level of information in file

        Call write_config(name,levcfg,megatm,i-1,eng_min/engunit,eng_0/engunit,netcdf,comm)
     End If

! setup new quaternions

     If (l_mov) Then
        If (lshmv_rgd) Then
          Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,xxx,yyy,zzz,comm)
        End If
        Call q_setup(comm)
     End If

  End If

  If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
     Deallocate (lstitr,           Stat=fail(1))
     call cons%deallocate_work()
     Call pmf%Deallocate_work()
  End If
  If (megrgd > 0) Then
     Deallocate (txx,tyy,tzz,      Stat=fail(6))
     Deallocate (uxx,uyy,uzz,      Stat=fail(7))
  End If
  Deallocate (gxx,gyy,gzz,         Stat=fail(8))
  If (Any(fail > 0)) Then
     Write(message,'(a,i0)') 'minimise_relax deallocation failure'
     Call error(0,message)
  End If

End Subroutine minimise_relax

Subroutine zero_k_optimise(stats,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for zero Kelvin temperature optimization
!
! in preparation for integration of equations of motion:
!
! the free particle velocity, V, is set in the direction of the force,
! F, when V.F > 0 :: V=F*[(V.F)/(F.F)] and set to zero if V.F < 0 :: V=0.
! the same rational is extended to RB dynamics where, additionally to
! applying the same strategy to the RB COM velocity change upon the COM
! force, alos the angular velocity of the RB, W, is set in the direction
! of the torque, T, for when W.T > 0 :: W=T*[(W.T)/(T.T)] and set to zero
! if W.T < 0 :: W=0.
!
! care must be taken for:
!     - to remove any COM motion generation
!     - to zero angular momentum about centre of mass for non-periodic
!       system (imcon=0)
!     - to ensure preservation of thermostat's instantaneous kinetic
!       energy components and remove spurious temperature fluctuations
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( stats_type ), Intent( InOut ) :: stats
  Type( comms_type ), Intent( InOut ) :: comm

  Integer           :: fail,i,j,i1,i2,irgd,jrgd,krgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: e_f,e_t,e_r,engkf,engkt,scale,vdotf,fsq, &
                       amx,amy,amz,wxx,wyy,wzz,tmp,tmp1,        &
                       com(1:3),vom(1:3),rot(1:9),rotinv(1:9),  &
                       x(1:1),y(1:1),z(1:1),                    &
                       fmx,fmy,fmz,tqx,tqy,tqz,trx,try,trz,     &
                       vpx,vpy,vpz

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer,ggx,ggy,ggz
  Character( Len = 256 ) :: message


! preserve magnitudes of old instantaneous energies in order to scale back to them

  e_t=0.5_wp*(stats%strknt(1)+stats%strknt(5)+stats%strknt(9)) ! RBs translational
  e_r=stats%engrot                                 ! RBs rotational
  e_f=stats%engke-e_t                              ! FPs (free particles)

! initialise RB energy components

  engkf=0.0_wp
  engkt=0.0_wp

! recover megrgd

  megrgd=rgdmeg

  If (megrgd > 0) Then
     fail=0
     Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), Stat=fail)
     If (fail > 0) Then
        Write(message,'(a)') 'zero_k_optimise allocation failure'
        Call error(0,message)
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

! Free particles

     Do j=1,nfree
        i=lstfre(j)

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then

! take component of velocity in direction of force

           vdotf = vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)

           If (vdotf < 0.0_wp) Then
              vxx(i) = 0.0_wp
              vyy(i) = 0.0_wp
              vzz(i) = 0.0_wp
           Else
              fsq = fxx(i)**2+fyy(i)**2+fzz(i)**2
              scale = vdotf/Max(1.0e-10_wp,fsq)

              vxx(i) = fxx(i)*scale
              vyy(i) = fyy(i)*scale
              vzz(i) = fzz(i)*scale
           End If

        End If
     End Do

! RBs

     krgd=0
     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then

! current rotation matrix

           Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

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

! calculate torque in principal frame

           trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
           try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
           trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

           If (rgdfrz(0,rgdtyp) == 0) Then

! take component of velocity in direction of force

              vdotf = rgdvxx(irgd)*fmx+rgdvyy(irgd)*fmy+rgdvzz(irgd)*fmz

              If (vdotf < 0.0_wp) Then
                 rgdvxx(irgd) = 0.0_wp
                 rgdvyy(irgd) = 0.0_wp
                 rgdvzz(irgd) = 0.0_wp
              Else
                 fsq = fmx**2+fmy**2+fmz**2
                 scale = vdotf/Max(1.0e-10_wp,fsq)

                 rgdvxx(irgd) = fmx*scale
                 rgdvyy(irgd) = fmy*scale
                 rgdvzz(irgd) = fmz*scale
              End If

           End If

! take component of the angular velocity in direction of
! the angular acceleration (torque./RI.)

           trx=trx*rgdrix(2,rgdtyp)
           try=try*rgdriy(2,rgdtyp)
           trz=trz*rgdriz(2,rgdtyp)

           vdotf = rgdoxx(irgd)*trx+rgdoyy(irgd)*try+rgdozz(irgd)*trz

           If (vdotf < 0.0_wp) Then
              rgdoxx(irgd) = 0.0_wp
              rgdoyy(irgd) = 0.0_wp
              rgdozz(irgd) = 0.0_wp
           Else
              fsq = trx**2+try**2+trz**2
              scale = vdotf/Max(1.0e-10_wp,fsq)

              rgdoxx(irgd) = trx*scale
              rgdoyy(irgd) = try*scale
              rgdozz(irgd) = trz*scale
           End If

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

! Subtract COM velocity

     Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,comm)

! remove centre of mass motion

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

     Call kinstresf(vxx,vyy,vzz,stats%strknf,comm)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,stats%strknt,comm)

     stats%strkin=stats%strknf+stats%strknt
     stats%engke=0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))

! update rotational energy

     stats%engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)

     Deallocate (ggx,ggy,ggz, Stat=fail)
     If (fail > 0) Then
        Write(message,'(a)') 'zero_k_optimise deallocation failure'
        Call error(0,message)
     End If
  Else
     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then

! take component of velocity in direction of force

           vdotf = vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)

           If (vdotf < 0.0_wp) Then
              vxx(i) = 0.0_wp
              vyy(i) = 0.0_wp
              vzz(i) = 0.0_wp
           Else
              fsq = fxx(i)**2+fyy(i)**2+fzz(i)**2
              scale = vdotf/Max(1.0e-10_wp,fsq)

              vxx(i) = fxx(i)*scale
              vyy(i) = fyy(i)*scale
              vzz(i) = fzz(i)*scale
           End If

        End If
     End Do

! Subtract COM velocity

     Call getvom(vom,vxx,vyy,vzz,comm)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do

! Update kinetic stress and energy

     Call kinstress(vxx,vyy,vzz,stats%strkin,comm)
     stats%engke=0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))
  End If

! zero angular momentum about centre of mass - non-periodic system

  If (imcon == 0) Then
     fail=0
     Allocate (buffer(1:12), Stat=fail)
     If (fail > 0) Then
        Write(message,'(a)') 'zero_k_optimise allocation failure'
        Call error(0,message)
     End If

! initialise RB energy components

     engkf=0.0_wp
     engkt=0.0_wp

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
              x(1)=(wyy*rgdzzz(irgd) - wzz*rgdyyy(irgd))
              y(1)=(wzz*rgdxxx(irgd) - wxx*rgdzzz(irgd))
              z(1)=(wxx*rgdyyy(irgd) - wyy*rgdxxx(irgd))

              rgdvxx(irgd) = rgdvxx(irgd) + x(1)
              rgdvyy(irgd) = rgdvyy(irgd) + y(1)
              rgdvzz(irgd) = rgdvzz(irgd) + z(1)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! local index of particle/site

                 If (i <= natms) Then
                    vxx(i) = vxx(i) + x(1)
                    vyy(i) = vyy(i) + y(1)
                    vzz(i) = vzz(i) + z(1)
                 End If
              End Do
           End If
        End Do

! get kinetic energy

        engkf=getknf(vxx,vyy,vzz,comm)
        engkt=getknt(rgdvxx,rgdvyy,rgdvzz,comm)
        stats%engke=engkf+engkt

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

! get kinetic energy

        stats%engke=getkin(vxx,vyy,vzz,comm)

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
        Write(message,'(a)') 'zero_k_optimise deallocation failure'
        Call error(0,message)
     End If
  End If

! to remove spurious temperature fluctuations ensure preservation
! of thermostat's instantaneous translational and rotational energies
! equipartitioning of DoFs may be lost as transfers of energy between
! translational and rotational DoFs may happen

! apply temperature components scaling

  engkf=stats%engke-engkt
  If (engkf+engkt+stats%engrot > 1.0e-6_wp .and. e_f+e_t+e_r > 1.0e-6_wp) Then
     If (megrgd > 0) Then
        tmp=Sqrt((e_f+e_t+e_r)/(engkf+engkt+stats%engrot))
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
                       x(1)=rgdx(jrgd,rgdtyp)
                       y(1)=rgdy(jrgd,rgdtyp)
                       z(1)=rgdz(jrgd,rgdtyp)

! site velocity in body frame

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

! update kinetic energy and stress

        Call kinstresf(vxx,vyy,vzz,stats%strknf,comm)
        Call kinstrest(rgdvxx,rgdvyy,rgdvzz,stats%strknt,comm)

        stats%strkin=stats%strknf+stats%strknt
        stats%engke=0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))

! update rotational energy

        stats%engrot=getknr(rgdoxx,rgdoyy,rgdozz,comm)
     Else
        tmp=Sqrt(e_f/stats%engke)
        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              vxx(i)=vxx(i)*tmp
              vyy(i)=vyy(i)*tmp
              vzz(i)=vzz(i)*tmp
           End If
        End Do

! Update kinetic stress and energy

        Call kinstress(vxx,vyy,vzz,stats%strkin,comm)
        stats%engke=0.5_wp*(stats%strkin(1)+stats%strkin(5)+stats%strkin(9))
     End If
  Else ! zero them and let's see if we can crash
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

End Subroutine zero_k_optimise


End Module minimise
