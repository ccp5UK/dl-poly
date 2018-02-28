Subroutine core_shell_relax(l_str,relaxed,lrdf,rlx_tol,megshl,stpcfg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for relaxing shells to zero force using conjugate
! gradient method
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith august 2014
! contrib   - a.m.elena february 2017
! contrib   - i.t.todorov february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,   Only : idnode,mxnode,gsum,gmax
  Use setup_module,   Only : engunit,nrite,output, &
                             mxatms,mxatdm,mxshl,zero_plus
  Use config_module,  Only : imcon,cell,natms,nlast,lsi,lsa, &
                             xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use parse_module,   Only : strip_blanks,lower_case
  Use kinetic_module, Only : freeze_atoms
  Use core_shell_module

  Implicit None

  Logical,           Intent( In    ) :: l_str
  Logical,           Intent( InOut ) :: relaxed,lrdf
  Integer,           Intent( In    ) :: megshl
  Real( Kind = wp ), Intent( In    ) :: rlx_tol(1:2),stpcfg

  Logical,           Save :: newjob = .true. , l_rdf
  Integer,           Save :: keyopt
  Integer                 :: fail(1:2),i,ia,ib,jshl,local_index
  Real( Kind = wp ), Save :: grad_tol,eng_tol,dist_tol(1:2),   &
                             step,eng,eng0,eng1,eng2,          &
                             grad,grad0,grad1,grad2,onorm,sgn, &
                             stride,gamma,x(1),y(1),z(1),fff(0:3)

! OUTPUT existence

  Logical               :: l_out
  Character( Len = 10 ) :: c_out

! Optimisation iteration and convergence limits

  Integer,      Parameter :: mxpass = 100

  Real( Kind = wp ), Save :: grad_pass

  Integer,           Allocatable       :: lstopt(:,:),lst_sh(:)
  Real( Kind = wp ), Allocatable       :: fxt(:),fyt(:),fzt(:)
  Real( Kind = wp ), Allocatable, Save :: oxt(:),oyt(:),ozt(:)

  fail=0
  Allocate (lstopt(1:2,1:mxshl),lst_sh(1:mxatms),      Stat=fail(1))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_relax allocation failure, node: ', idnode
     Call error(0)
  End If

  If (newjob) Then
     newjob = .false.

! At start no optimisation has been attempted yet

     keyopt = 0

! Passage accumulators are initialised in core_shell_module
! passshl(1) - cycles counter
! passshl(2) - access counter
! passshl(3) - average cycles
! passshl(4) - minimum cycles
! passshl(5) - maximum cycles

     Allocate (oxt(1:mxshl),oyt(1:mxshl),ozt(1:mxshl), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'core_shell_relax allocation failure SAVE, node: ', idnode
        Call error(0)
     End If
  End If

! Step length for relaxation

  If (rlx_tol(2) > zero_plus) Then

! Optionally specified

     step=rlx_tol(2)

  Else

! default if unspecified

     step=0.5_wp/smax

  End If

  If (keyopt == 0) Then

! No relaxation is yet attempted

     relaxed=.false.

! Minimum needed for a pass for this minimisation cycle

     grad_pass = Huge(1.0_wp)

! Avoid rdf calculation redundancy

     l_rdf=lrdf
     If (lrdf) lrdf=.false.

! Print header

     If (l_str .and. idnode == 0) Then
        Write(nrite, Fmt=*)
        Write(nrite,'(1x,a,3x,a,6x,a,11x,a,8x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
  'Relaxing shells to cores:','pass','eng_tot','grad_tol','dis_tol','dcs_max','tol=',rlx_tol(1),'step=',step
        Write(nrite,"(1x,130('-'))")
     End If
  End If

! gather new forces on shared shells

  If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,fxx,fyy,fzz)

! Load shell forces on cores (cores don't move during the shell relaxation)

  fxt(1:mxshl)=0.0_wp ; fyt(1:mxshl)=0.0_wp ; fzt(1:mxshl)=0.0_wp
  jshl=0
  Do i=1,ntshl
     ia=local_index(listshl(1,i),nlast,lsi,lsa)
     ib=local_index(listshl(2,i),nlast,lsi,lsa)
     If (ia > 0 .and. ia <= natms) Then ! THERE IS AN ib>0 FOR SURE
        jshl=jshl+1
        fxt(jshl)=fxx(ib)
        fyt(jshl)=fyy(ib)
        fzt(jshl)=fzz(ib)
     End If
     lstopt(1,i)=ia
     lstopt(2,i)=ib
  End Do

! Current configuration energy

  eng=stpcfg

! Initialise/get eng_tol & verify relaxed condition

  eng_tol=0.0_wp
  If (keyopt > 0) Then
     eng_tol=Abs(1.0_wp-eng2/eng)
  End If
! Current gradient (modulus of the total force on shells)

  grad=0.0_wp
  Do i=1,jshl
     grad=grad+fxt(i)**2+fyt(i)**2+fzt(i)**2
  End Do
  If (mxnode > 1) Call gsum(grad)
  grad=Sqrt(grad)

! Get grad_tol & verify relaxed condition

  grad_tol=grad/Real(megshl,wp)
  relaxed=(grad_tol < rlx_tol(1))

! Initialise dist_tol

  dist_tol=0.0_wp

! CHECK FOR CONVERGENCE

  If (.not.relaxed) Then

! Increment main passage counter

     passshl(1)=passshl(1)+1.0_wp

! Minimum for passing

     grad_pass = Min(grad_pass,grad_tol)

! If in mxpass iterations we are not there, give up but
! allow for ten-fold boost in iteration cycle length
! for the very first MD step

     If (Nint(passshl(2)) == 0) Then
        If (Nint(passshl(1)) >= 10*mxpass) Then
           Call warning(330,rlx_tol(1),grad_pass,0.0_wp)
           Call error(474)
        End If
     Else
        If (Nint(passshl(1)) >= mxpass) Then
           Call warning(330,rlx_tol(1),grad_pass,0.0_wp)
           Call error(474)
        End If
     End If

  Else

     Go To 100

  End If

  If      (keyopt == 0) Then

! Original configuration energy

     eng0=eng
     eng2=eng

! Original gradient (modulus of the total force on shells)

     onorm=grad
     grad0=grad
     grad2=grad

! Set original search direction

     oxt=0.0_wp ; oyt=0.0_wp ; ozt=0.0_wp
     Do i=1,jshl
        oxt(i)=fxt(i)
        oyt(i)=fyt(i)
        ozt(i)=fzt(i)
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
     Do i=1,jshl
        grad2=grad2+oxt(i)*fxt(i)+oyt(i)*fyt(i)+ozt(i)*fzt(i)
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
     Do i=1,jshl
        oxt(i)=fxt(i)+gamma*oxt(i)
        oyt(i)=fyt(i)+gamma*oyt(i)
        ozt(i)=fzt(i)+gamma*ozt(i)

        onorm=onorm+oxt(i)**2+oyt(i)**2+ozt(i)**2
        grad2=grad2+oxt(i)*fxt(i)+oyt(i)*fyt(i)+ozt(i)*fzt(i)
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

! Load original shell forces on their cores in DD representation

  fxt(1:mxatdm)=0.0_wp ; fyt(1:mxatdm)=0.0_wp ; fzt(1:mxatdm)=0.0_wp
  jshl=0
  Do i=1,ntshl
     ia=lstopt(1,i)
     If (ia > 0 .and. ia <= natms) Then
        jshl=jshl+1
        fxt(ia)=oxt(jshl)
        fyt(ia)=oyt(jshl)
        fzt(ia)=ozt(jshl)
     End If
  End Do

! Exchange original shell forces on shared cores across domains

  If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,fxt,fyt,fzt)

! Move shells accordingly to their new positions

  Do i=1,ntshl
     ia=lstopt(1,i)
     ib=lstopt(2,i)
     If (ia > 0 .and. (ib > 0 .and. ib <= natms)) Then
        xxx(ib)=xxx(ib)+stride*fxt(ia)
        yyy(ib)=yyy(ib)+stride*fyt(ia)
        zzz(ib)=zzz(ib)+stride*fzt(ia)
        dist_tol(1)=Max(dist_tol(1),fxt(ia)**2+fyt(ia)**2+fzt(ia)**2) ! - shell move
        x(1)=xxx(ib)-xxx(ia) ; y(1)=yyy(ib)-yyy(ia) ; z(1)=zzz(ib)-zzz(ia)
        Call images(imcon,cell,1,x,y,z)
        dist_tol(2)=Max(dist_tol(2),x(1)**2+y(1)**2+z(1)**2) ! - core-shell separation
     End If
  End Do
  dist_tol=Sqrt(dist_tol)
  dist_tol(1)=dist_tol(1)*Abs(stride)
  If (mxnode > 1) Call gmax(dist_tol)

! Fit headers in and Close and Open OUTPUT at every 25th print-out

  i=Nint(passshl(1))
  If (l_str .and. idnode == 0) Then
     Write(nrite,'(1x,i31,1x,1p,2e18.8,4x,f7.4,4x,f7.4,12x,e18.8)') i-1,stpcfg/engunit,grad_tol,dist_tol(1),dist_tol(2),eng_tol
     If (Mod(i,25) == 0) Then
        Write(nrite,"(1x,130('-'))")
        Write(nrite,'(1x,a,3x,a,6x,a,11x,a,9x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
  'Relaxing shells to cores:','pass','eng_tot','grad_tol','ds_tol','dcs_max','tol=',rlx_tol(1),'step=',step
        Write(nrite,"(1x,130('-'))")

        If (idnode == 0) Then
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

  If (relaxed) Then

! Final printout

     i=Nint(passshl(1))
     If (idnode == 0) Then
        If (.not.l_str) Then
           Write(nrite, Fmt=*)
           Write(nrite,'(1x,a,4x,a,6x,a,11x,a,8x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
  'Relaxed shells to cores:','pass','eng_tot','grad_tol','dis_tol','dcs_max','tol=',rlx_tol(1),'step=',step
           Write(nrite,"(1x,130('-'))")
        End If
        Write(nrite,'(1x,i31,1x,1p,2e18.8,4x,f7.4,4x,f7.4,12x,e18.8)') &
             i-1,stpcfg/engunit,grad_tol,dist_tol(1),dist_tol(2),eng_tol
        Write(nrite, Fmt=*)
        Write(nrite,"(1x,130('-'))")
     End If

! Collect passage statistics

     passshl(3)=passshl(2)*passshl(3)
     passshl(2)=passshl(2)+1.0_wp
     passshl(3)=passshl(3)/passshl(2)+passshl(1)/passshl(2)
     passshl(4)=Min(passshl(1),passshl(4))
     passshl(5)=Max(passshl(1),passshl(5))

! Rewind keyopt and main passage counter

     keyopt =0
     passshl(1)=0.0_wp

! Resume rdf calculations

     If (l_rdf) lrdf=l_rdf

! Zero shells' velocities and forces and redistribute
! the residual force to the rest of the system to prevent
! COM force generation

     lst_sh(1:natms)=0
     fff(0)=Real(natms,wp)
     fff(1:3)=0.0_wp
     Do i=1,ntshl
        ib=lstopt(2,i)
        If (ib > 0) Then
           lst_sh(ib)=1
           fff(0)=fff(0)-1.0_wp
           fff(1)=fff(1)+fxx(ib) ; fxx(ib)=0.0_wp ; vxx(ib)=0.0_wp
           fff(2)=fff(2)+fyy(ib) ; fyy(ib)=0.0_wp ; vyy(ib)=0.0_wp
           fff(3)=fff(3)+fzz(ib) ; fzz(ib)=0.0_wp ; vzz(ib)=0.0_wp
        End If
     End Do
     If (mxnode > 1) Call gsum(fff)
     fff(1:3)=fff(1:3)/fff(0)
     Do i=1,natms
        If (lst_sh(i) == 0) Then
           fxx(i)=fxx(i)+fff(1)
           fyy(i)=fyy(i)+fff(2)
           fzz(i)=fzz(i)+fff(3)
        End If
     End Do

! Frozen atoms option

     Call freeze_atoms()
  End If

  Deallocate (lstopt,lst_sh, Stat=fail(1))
  Deallocate (fxt,fyt,fzt,   Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_relax deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine core_shell_relax
