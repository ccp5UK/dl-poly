Subroutine core_shell_relax(l_str,relaxed,lrdf,rlx_tol,megshl,stpcfg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for relaxing shells to zero force using conjugate
! gradient method
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith march 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsum
  Use setup_module,      Only : nrite,mxatms,mxshl,engunit
  Use config_module,     Only : natms,nlast,lsi,lsa, &
                                xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use core_shell_module
  Use statistics_module, Only : pass
  Use kinetic_module,    Only : freeze_atoms

  Implicit None

  Logical,           Intent( In    ) :: l_str
  Logical,           Intent( InOut ) :: relaxed,lrdf
  Integer,           Intent( In    ) :: megshl
  Real( Kind = wp ), Intent( In    ) :: rlx_tol,stpcfg

  Logical,           Save :: newjob = .true. , l_rdf
  Integer,           Save :: keyopt
  Integer                 :: fail(1:2),i,ia,ib,jshl,local_index
  Real( Kind = wp ), Save :: grad_tol,step,eng,eng0,eng1,eng2, &
                             grad,grad0,grad1,grad2,onorm,sgn, &
                             stride,gamma,fff(0:3)

! Optimisation iteration and convergence limits

  Integer,      Parameter :: mxpass = 100

  Real( Kind = wp ), Save :: grad_pass

  Integer,           Allocatable       :: lstopt(:,:),lst_sh(:)
  Real( Kind = wp ), Allocatable       :: fxt(:),fyt(:),fzt(:)
  Real( Kind = wp ), Allocatable, Save :: oxt(:),oyt(:),ozt(:)

  fail=0
  Allocate (lstopt(1:2,1:mxshl),lst_sh(1:mxatms),       Stat=fail(1))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'core_shell_relax allocation failure, node: ', idnode
     Call error(0)
  End If

  If (newjob) Then
     newjob = .false.

! At start no optimisation has been attempted yet

     keyopt = 0

! Passage accumulators are zeroed in statistics_module
! pass(1) - cycles counter
! pass(2) - access counter
! pass(3) - average cycles
! pass(4) - minimum cycles
! pass(5) - maximum cycles

     Allocate (oxt(1:mxshl),oyt(1:mxshl),ozt(1:mxshl), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'core_shell_relax allocation failure SAVE, node: ', idnode
        Call error(0)
     End If
  End If

! step length for relaxation

  step=0.5_wp/smax

  If (keyopt == 0) Then

! No relaxation is yet attempted

     relaxed=.false.

! Minimum grad_pass = 10*rlx_tol

     grad_pass = 10.0_wp*rlx_tol

! Avoid rdf calculation redundancy

     l_rdf=lrdf
     If (lrdf) lrdf=.false.

! Print header

     If (l_str .and. idnode == 0) Then
        Write(nrite, Fmt=*)
        Write(nrite,'(1x,a,6x,a,10x,a,11x,a,9x,a,1p,e12.4)') &
             'Relaxing shells to cores:','pass','eng_tot','grad_tol','tol=', rlx_tol
        Write(nrite,"(1x,130('-'))")
     End If
  End If

! gather new forces on shared shells

  If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,fxx,fyy,fzz)

! Load shell forces on cores (cores don't move during the shell relaxation)

  fxt=0.0_wp ; fyt=0.0_wp ; fzt=0.0_wp
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

! Current gradient (modulus of the total force on shells)

  grad=0.0_wp
  Do i=1,jshl
     grad=grad+fxt(i)**2+fyt(i)**2+fzt(i)**2
  End Do
  If (mxnode > 1) Call gsum(grad)
  grad=Sqrt(grad)

! Get grad_tol & verify relaxed condition

  grad_tol=grad/Real(megshl,wp)
  relaxed=(grad_tol < rlx_tol)

! CHECK FOR CONVERGENCE

  If (.not.relaxed) Then

! Increment main passage counter

     pass(1)=pass(1)+1.0_wp

! Minimum grad_pass = 10*rlx_tol

     grad_pass = Min(grad_pass,grad_tol)

! If in mxpass iterations we are not there, give up but
! allow for ten-fold boost in iteration cycle length
! for the very first MD step

     If (Nint(pass(2)) == 0) Then
        If (Nint(pass(1)) >= 10*mxpass) Then
           Call warning(330,rlx_tol,grad_pass,0.0_wp)
           Call error(474)
        End If
     Else
        If (Nint(pass(1)) >= mxpass) Then
           Call warning(330,rlx_tol,grad_pass,0.0_wp)
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

  fxt=0.0_wp ; fyt=0.0_wp ; fzt=0.0_wp
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
     End If
  End Do

  i=Nint(pass(1))
  If (l_str .and. idnode == 0) Then
     Write(nrite,'(1x,i34,4x,1p,2e18.8)') i-1,stpcfg/engunit,grad_tol
     If (Mod(i,25) == 0) Then
        Write(nrite,"(1x,130('-'))")
        Write(nrite,'(1x,a,6x,a,10x,a,11x,a,9x,a,1p,e12.4)') &
             'Relaxing shells to cores:','pass','eng_tot','grad_tol','tol=', rlx_tol
        Write(nrite,"(1x,130('-'))")
     End If
  End If

100 Continue

  If (relaxed) Then

! Final printout

     i=Nint(pass(1))
     If (l_str .and. idnode == 0) Then
        Write(nrite,'(1x,i34,4x,1p,2e18.8)') i,stpcfg/engunit,grad_tol
        Write(nrite, Fmt=*)
        Write(nrite,"(1x,130('-'))")
     End If

! Collect passage statistics

     pass(3)=pass(2)*pass(3)
     pass(2)=pass(2)+1.0_wp
     pass(3)=pass(3)/pass(2)+pass(1)/pass(2)
     pass(4)=Min(pass(1),pass(4))
     pass(5)=Max(pass(1),pass(5))

! Rewind keyopt and main passage counter

     keyopt =0
     pass(1)=0.0_wp

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
