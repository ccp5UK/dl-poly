Subroutine vnl_check(l_str,rcut,rpad,rlnk,width)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine implementing the VNL conditional update checks
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2013
! contrib   - i.j.bush february 2014
! contrib   - i.t.todorov january 2017
! contrib   - i.t.todorov july 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gmax
  Use setup_module,   Only : nrite,mxspl,mxatdm,half_minus
  Use domains_module, Only : r_nprx,r_npry,r_nprz
  Use config_module,  Only : imcon,cell,natms,xxx,yyy,zzz
  Use vnl_module

  Implicit None

  Logical,           Intent ( In    ) :: l_str
  Real( Kind = wp ), Intent ( In    ) :: rcut
  Real( Kind = wp ), Intent ( InOut ) :: rpad,rlnk,width

  Logical, Save :: newstart=.true.

  Logical           :: safe
  Integer           :: fail,ilx,ily,ilz,i,ii,j
  Real( Kind = wp ) :: tol,cut,test,celprp(1:10)

  Real( Kind = wp ), Allocatable :: x(:),y(:),z(:),r(:)

  If (.not.llvnl) Return

! Checks

  fail = 0
  Allocate (x(1:mxatdm),y(1:mxatdm),z(1:mxatdm),r(1:mxatdm), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vnl_check allocation failure, node: ', idnode
     Call error(0)
  End If

  x(1:natms) = xxx(1:natms) - xbg(1:natms)
  y(1:natms) = yyy(1:natms) - ybg(1:natms)
  z(1:natms) = zzz(1:natms) - zbg(1:natms)

  Call images(imcon,cell,natms,x,y,z)

  If (natms > 0) r(1:natms) = Sqrt(x(1:natms)**2 + y(1:natms)**2 + z(1:natms)**2)

! search for violations: local domain then global

  tol=Merge( Maxval( r(1:natms) ) , 0.0_wp , natms > 0)

  Deallocate (x,y,z,r, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vnl_check deallocation failure, node: ', idnode
     Call error(0)
  End If

  If (mxnode > 1) Call gmax(tol)

! decide on update if unsafe

  l_vnl = (tol >= half_minus*rpad)

! get the dimensional properties of the MD cell

  Call dcell(cell,celprp)
  width=Min(celprp(7),celprp(8),celprp(9))

! define cut as in link_cell_pairs

  cut=rlnk+1.0e-6_wp

! calculate link cell dimensions per node

  ilx=Int(r_nprx*celprp(7)/cut)
  ily=Int(r_npry*celprp(8)/cut)
  ilz=Int(r_nprz*celprp(9)/cut)

  tol=Min(0.05_wp,0.005_wp*rcut)                                        ! tolerance
  test = 0.02_wp * Merge( 1.0_wp, 2.0_wp, mxspl > 0)                    ! 2% (w/ SPME) or 4% (w/o SPME)
  cut=Min(r_nprx*celprp(7),r_npry*celprp(8),r_nprz*celprp(9))-1.0e-6_wp ! domain size

  If (ilx*ily*ilz == 0) Then
     If (cut < rcut) Then
        If (idnode == 0) Write(nrite,*) '*** warning - rcut <= Min(domain width) < rlnk = rcut + rpad !!! ***'
        Call error(307)
     Else ! rpad is defined & in 'no strict' mode
        If (cut < rlnk) Then
           If (l_str) Then
              If (idnode == 0) Write(nrite,*) '*** warning - rcut <= Min(domain width) < rlnk = rcut + rpad !!! ***'
              Call error(307)
           Else
              If (cut >= rcut) Then ! Re-set rpad with some slack
                 rpad = Min( 0.95_wp * (cut - rcut) , test * rcut)
                 rpad = Real( Int( 100.0_wp * rpad , wp ) ) / 100.0_wp
                 If (rpad < tol) rpad = 0.0_wp ! Don't bother
                 rlnk = rcut + rpad
                 l_vnl=.true.
              End If
           End If
        End If
     End If
  Else ! push the limits when up for update in a 'no strict' regime
     If (l_vnl .and. (.not.l_str)) Then ! Try to re-set rpad with some slack
        If (Int(Real(Min(ilx,ily,ilz),wp)/(1.0_wp+test)) >= 2) Then
           cut = test * rcut
        Else
           If (mxnode > 1) Then
              cut = Min( 0.95_wp * ( Min ( r_nprx * celprp(7) / Real(ilx,wp) , &
                                           r_npry * celprp(8) / Real(ily,wp) , &
                                           r_nprz * celprp(9) / Real(ilz,wp) ) &
                                     - rcut - 1.0e-6_wp ) , test * rcut )
           Else ! catch & handle exception
              cut = 0.95_wp * (0.5_wp*width - rcut - 1.0e-6_wp)
           End If
        End If
        cut = Real( Int( 100.0_wp * cut ) , wp ) / 100.0_wp
        If ((.not.(cut < tol)) .and. cut-rpad > 0.005_wp) Then ! Do bother
  If (idnode == 0) Write(nrite,'(/,1x,2(a,f5.2),a,/)') 'cutoff padding reset from ', rpad, ' Angs to ', cut, ' Angs'
           rpad = cut
           rlnk = rcut + rpad
        End If
     End If
  End If

  If (l_vnl) Then ! Deal with skipping statistics
     skipvnl(3)=skipvnl(2)*skipvnl(3)
     skipvnl(2)=skipvnl(2)+1.0_wp
     skipvnl(3)=skipvnl(3)/skipvnl(2)+skipvnl(1)/skipvnl(2)
     If (.not.newstart) Then ! avoid first compulsory force evaluation
        skipvnl(4)=Min(skipvnl(1),skipvnl(4))
     Else
        newstart=.false.
     End If
     skipvnl(5)=Max(skipvnl(1),skipvnl(5))

     skipvnl(1) = 0.0_wp              ! Reset here, checkpoit set by vnl_set_check in set_halo_particles
  Else            ! Enjoy telephoning
     skipvnl(1) = skipvnl(1) + 1.0_wp ! Increment, telephony done for xxx,yyy,zzz in set_halo_positions
  End If

End Subroutine vnl_check
