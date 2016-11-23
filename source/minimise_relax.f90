Subroutine minimise_relax &
           (l_str,relaxed,lrdf,megatm,megcon,megpmf,megrgd, &
           keymin,min_tol,tstep,stpcfg)

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
! contrib   - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode,mxnode,gsum,gmax
  Use setup_module,        Only : nrite,mxatms,mxcons,mxtpmf, & 
                                  mxpmf,engunit,output,zero_plus
  Use config_module,       Only : natms,nlast,nfree,          &
                                  lsi,lsa,lfrzn,lfree,lstfre, &
                                  weight,xxx,yyy,zzz,fxx,fyy,fzz
  Use rigid_bodies_module, Only : lshmv_rgd,lishp_rgd,lashp_rgd
  Use parse_module,        Only : strip_blanks,lower_case
  Use minimise_module

  Implicit None

  Logical,           Intent( In    ) :: l_str
  Logical,           Intent( InOut ) :: relaxed,lrdf
  Integer,           Intent( In    ) :: megatm,megcon, &
                                        megpmf,megrgd,keymin
  Real( Kind = wp ), Intent( In    ) :: min_tol(1:2),tstep,stpcfg

  Logical,              Save :: newjob = .true. , l_rdf, l_mov
  Character( Len = 8 ), Save :: word
  Character( Len = 6 )       :: name
  Integer,              Save :: keyopt
  Integer                    :: fail(1:8),i,j,levcfg
  Real( Kind = wp ),    Save :: total,grad_tol,eng_tol,dist_tol,step,           &
                                eng_0,eng_min,engcon,engpmf,eng,eng0,eng1,eng2, &
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
  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)
  Real( Kind = wp ), Allocatable :: txx(:),tyy(:),tzz(:)
  Real( Kind = wp ), Allocatable :: uxx(:),uyy(:),uzz(:)


  fail=0
  If (megcon > 0 .or. megpmf > 0) Then
     Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
     If (megcon > 0) Then
        Allocate (lstopt(0:2,1:mxcons),listot(1:mxatms),          Stat=fail(2))
        Allocate (dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons),      Stat=fail(3))
     End If
     If (megpmf > 0) Then
        Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat=fail(4))
        Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),         Stat=fail(5))
     End If
  End If
  If (megrgd > 0) Then
     Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),         Stat=fail(6))
     Allocate (uxx(1:mxatms),uyy(1:mxatms),uzz(1:mxatms),         Stat=fail(7))
  End If
  Allocate (gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms),            Stat=fail(8))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'minimise_relax allocation failure, node: ', idnode
     Call error(0)
  End If

  If (newjob) Then
     newjob = .false.

! At start no optimisation has been attempted yet

     keyopt = 0

! At start the minimum energy is defined as zero

     eng_min = 0.0_wp

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
     If (mxnode > 1) Call gsum(total)
  End If

! Step length for relaxation

  If (min_tol(2) > zero_plus) Then

! Optionally specified

     step=min_tol(2)

  Else

! default if unspecified

     step=tstep**2

! enlarged depending on functionality if defaulted

     If (megcon == 0 .and. megpmf == 0) Then
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

     l_rdf=lrdf
     If (lrdf) lrdf=.false.

! Determine optimisation

     If      (keymin == 0) Then
        word='force   '
     Else If (keymin == 1) Then
        word='energy  '
     Else If (keymin == 2) Then
        word='distance'
     End If

! Print header

     If (l_str .and. idnode == 0) Then
        Write(nrite, Fmt=*)
        Write(nrite,'(3(1x,a),6x,a,10x,a,10x,a,11x,a,9x,a,1p,e12.4)') &
             'Minimising',word,'pass','eng_tot','grad_tol','eng_tol','dist_tol','tol=', min_tol(1)
        Write(nrite,"(1x,130('-'))")
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

  If (megcon > 0 .or. megpmf > 0) Then
     lstitr(1:natms)=.false. ! initialise lstitr

     If (megcon > 0) Then
        Call constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot)
        Call constraints_pseudo_bonds(lstopt,dxx,dyy,dzz,gxx,gyy,gzz,engcon)
        eng=eng+engcon
     End If

     If (megpmf > 0) Then
        Call pmf_tags(lstitr,indpmf,pxx,pyy,pzz)
        Call pmf_pseudo_bonds(indpmf,pxx,pyy,pzz,gxx,gyy,gzz,engpmf)
        eng=eng+engpmf
     End If
  End If

! Average forces over all members of a RB and split torques accordingly

  If (megrgd > 0) Then
     If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,gxx,gyy,gzz)
     Call rigid_bodies_split_torque(gxx,gyy,gzz,txx,tyy,tzz,uxx,uyy,uzz)
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
  If (mxnode > 1) Call gsum(grad)
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
! allow for ten-fold boost in iteration cycle length
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
     If (mxnode > 1) Call gsum(grad2)
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
     If (mxnode > 1) Call gsum(onorm)
     onorm=Sqrt(onorm)
     If (mxnode > 1) Call gsum(grad2)
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
  If (mxnode > 1) Call gmax(dist_tol)

  If (keymin == 2) relaxed=(dist_tol < min_tol(1))

! Fit headers in and Close and Open OUTPUT at every 25th print-out

  i=Nint(passmin(1))
  If (l_str .and. idnode == 0) Then
     Write(nrite,'(1x,i23,1p,4e18.8)') i-1,eng/engunit,grad_tol,eng_tol,dist_tol
     If (Mod(i,25) == 0) Then
        Write(nrite,"(1x,130('-'))")
        Write(nrite,'(3(1x,a),6x,a,10x,a,10x,a,11x,a,9x,a,1p,e12.4)') &
             'Minimising',word,'pass','eng_tot','grad_tol','eng_tol','dist_tol','tol=', min_tol(1)
        Write(nrite,"(1x,130('-'))")

        If (idnode == 0) Then
           Inquire(File=trim(output), Exist=l_out, Position=c_out)
           Call strip_blanks(c_out)
           Call lower_case(c_out)
           If (l_out .and. c_out(1:6) == 'append') Then
              Close(Unit=nrite)
              Open(Unit=nrite, File=trim(output), Position='append')
           End If
        End If
     End If
  End If

100 Continue

  l_x=(.not.relaxed) ! Transportation flag
  If (relaxed) Then

! Final/Only printout

     i=Nint(passmin(1))
     If (idnode == 0) Then
        If (.not.l_str) Then
           Write(nrite, Fmt=*)
           Write(nrite,'(3(1x,a),5x,a,10x,a,10x,a,11x,a,9x,a,1p,e12.4)') &
                'Minimised',word,'passes','eng_tot','grad_tol','eng_tol','dist_tol','tol=', min_tol(1)
           Write(nrite,"(1x,130('-'))")
        End If
        Write(nrite,'(1x,i23,1p,4e18.8)') i,eng/engunit,grad_tol,eng_tol,dist_tol
        Write(nrite, Fmt=*)
        Write(nrite,"(1x,130('-'))")
     End If

! Collect passage statistics

     passmin(3)=passmin(2)*passmin(3)
     passmin(2)=passmin(2)+1.0_wp
     passmin(3)=passmin(3)/passmin(2)+passmin(1)/passmin(2)
     passmin(4)=Min(passmin(1),passmin(4))
     passmin(5)=Max(passmin(1),passmin(5))

! Rewind keyopt and main passage counter

     keyopt =0
     passmin(1)=0.0_wp

! Resume rdf calculations

     If (l_rdf) lrdf=l_rdf

! Deallocate working arrays

     Call deallocate_minimise_arrays()

! Dump the lowest energy configuration

     If (eng < eng_min) Then
        eng_min=eng

        name = 'CFGMIN' ! file name
        levcfg = 0      ! define level of information in file

        Call write_config(name,levcfg,megatm,i-1,eng_min/engunit,eng_0/engunit)
     End If

! setup new quaternions

     If (l_mov) Then
        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,xxx,yyy,zzz)
        Call q_setup()
     End If

  End If

  If (megcon > 0 .or. megpmf > 0) Then
     Deallocate (lstitr,           Stat=fail(1))
     If (megcon > 0) Then
        Deallocate (lstopt,listot, Stat=fail(2))
        Deallocate (dxx,dyy,dzz,   Stat=fail(3))
     End If
     If (megpmf > 0) Then
        Deallocate (indpmf,        Stat=fail(4))
        Deallocate (pxx,pyy,pzz,   Stat=fail(5))
     End If
  End If
  If (megrgd > 0) Then
     Deallocate (txx,tyy,tzz,      Stat=fail(6))
     Deallocate (uxx,uyy,uzz,      Stat=fail(7))
  End If
  Deallocate (gxx,gyy,gzz,         Stat=fail(8))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'minimise_relax deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine minimise_relax
