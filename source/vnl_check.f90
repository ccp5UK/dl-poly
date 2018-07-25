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
  Use comms_module,   Only : idnode,mxnode,gcheck
  Use setup_module,   Only : nrite,mxspl,mxatms
  Use domains_module, Only : r_nprx,r_npry,r_nprz
  Use config_module,  Only : imcon,cell,natms,nlast,list, &
                             xxx,yyy,zzz
  Use vnl_module

  Implicit None

  Logical,           Intent ( In    ) :: l_str
  Real( Kind = wp ), Intent ( In    ) :: rcut
  Real( Kind = wp ), Intent ( InOut ) :: rpad,rlnk,width

  Logical, Save :: newstart=.true.

  Logical           :: safe
  Integer           :: fail,ilx,ily,ilz,i,ii,j
  Real( Kind = wp ) :: tol,cut,test,xx(1),yy(1),zz(1),celprp(1:10)

  Real( Kind = wp ), Allocatable :: x(:),y(:),z(:),r(:)

  If (.not.llvnl) Return

! Checks

  fail = 0
  Allocate (x(1:mxatms),y(1:mxatms),z(1:mxatms),r(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vnl_check allocation failure, node: ', idnode
     Call error(0)
  End If

  x(1:nlast) = xxx(1:nlast) - xbg(1:nlast)
  y(1:nlast) = yyy(1:nlast) - ybg(1:nlast)
  z(1:nlast) = zzz(1:nlast) - zbg(1:nlast)

  Call images(imcon,cell,nlast,x,y,z)

  If (nlast > 0) r(1:nlast) = Sqrt(x(1:nlast)**2 + y(1:nlast)**2 + z(1:nlast)**2)

! Start search for violations

  tol=0.5_wp*rpad
  safe=.true.
Q:Do i=1,natms
     If (r(i) > tol) Then
        If (r(i) > rpad) Then ! limit full particle loss
           safe=.false.
           Exit Q
        End If

! Limit pair visibility loss

        Do ii=1,list(-2,i)
           j=list(ii,i) ! may be in the halo!!!
           xx(1)=x(i)-x(j)
           yy(1)=y(i)-y(j)
           zz(1)=z(i)-z(j)
           cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
           If (cut > rpad) Then
              xx(1)=xbg(i)-xbg(j)
              yy(1)=ybg(i)-ybg(j)
              zz(1)=zbg(i)-zbg(j)
              cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
              If (cut < rcut) Then
                 xx(1)=xxx(i)-xxx(j)
                 yy(1)=yyy(i)-yyy(j)
                 zz(1)=zzz(i)-zzz(j)
                 cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
                 If (cut > rlnk) Then
                    safe=.false. ! pair distance violation
                    Exit Q
                 End If
              End If
           End If
        End Do
     End If
  End Do Q

  Deallocate (x,y,z,r, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vnl_check deallocation failure, node: ', idnode
     Call error(0)
  End If

  If (mxnode > 1) Call gcheck(safe,"enforce")
  l_vnl=.not.safe

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
