Module vnl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring VNL conditional update variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
! contrib   - i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,   Only : comms_type,gcheck
  Use setup,   Only : nrite,mxspl,mxatms
  Use domains, Only : r_nprx,r_npry,r_nprz
  Use configuration,  Only : imcon,cell,natms,nlast,list, &
                             xxx,yyy,zzz
  Use errors_warnings, Only : error
  Use numerics, Only : dcell,images
  Implicit None

  Logical,                        Save :: l_vnl = .true., & ! Do update
                                          llvnl = .false.   ! Unconditional VNL

  Real( Kind = wp ),              Save :: skipvnl(1:5) = (/ &
                                          0.0_wp         ,  & ! skips counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average skips
                                          999999999.0_wp ,  & ! minimum skips : ~Huge(1)
                                          0.0_wp /)           ! maximum skips


  Real( Kind = wp ), Allocatable, Save :: xbg(:),ybg(:),zbg(:)

  Public :: vnl_set_check

Contains

Subroutine vnl_check(l_str,rcut,rpad,rlnk,width,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine implementing the VNL conditional update checks
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
! contrib   - i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent ( In    ) :: l_str
  Real( Kind = wp ), Intent ( In    ) :: rcut
  Real( Kind = wp ), Intent ( InOut ) :: rpad,rlnk,width
  Type( comms_type ), Intent ( InOut ) :: comm

  Logical, Save :: newstart=.true.

  Logical           :: safe
  Integer           :: fail,ilx,ily,ilz,i,ii,j
  Real( Kind = wp ) :: cut,test,tol,celprp(1:10),xx(1),yy(1),zz(1)

  Real( Kind = wp ), Allocatable :: x(:),y(:),z(:),r(:)

  Character( Len = 256 ) :: message
  If (.not.llvnl) Return

! Checks

  fail = 0
  Allocate (x(1:mxatms),y(1:mxatms),z(1:mxatms),r(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(message,'(/,1x,a)') 'vnl_check allocation failure'
     Call error(0,message)
  End If

  x(1:nlast) = xxx(1:nlast) - xbg(1:nlast)
  y(1:nlast) = yyy(1:nlast) - ybg(1:nlast)
  z(1:nlast) = zzz(1:nlast) - zbg(1:nlast)

  Call images(imcon,cell,nlast,x,y,z)

  If (nlast > 0) r(1:nlast) = Sqrt(x(1:nlast)**2 + y(1:nlast)**2 + z(1:nlast)**2)

  safe=.true.
Q:Do i=1,natms
     If (r(i) > 0.5_wp*rpad) Then
        Do ii=1,list(-2,i)
           j=list(ii,i) ! may be in the halo!!!
           If (r(j) > 0.5_wp*rpad) Then
              xx(1)=x(i)-x(j)
              yy(1)=y(i)-y(j)
              zz(1)=z(i)-z(j)
              cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
              If (cut > rpad) Then
                 xx(1)=xxx(i)-xxx(j)
                 yy(1)=yyy(i)-yyy(j)
                 zz(1)=zzz(i)-zzz(j)
                 cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
                 If (cut > rlnk) Then
                    safe=.false. ! strong violation
                    Exit Q
                 End If
              End If
           End If
        End Do
     End If
  End Do Q

  Deallocate (x,y,z,r, Stat = fail)
  If (fail > 0) Then
     Write(message,'(/,1x,a)') 'vnl_check deallocation failure'
     Call error(0,message)
  End If

  Call gcheck(comm,safe,"enforce")
  l_vnl=.not.safe

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp) ! cut=rlnk+1.0e-6_wp in link_cell_pairs
  width=Min(celprp(7),celprp(8),celprp(9))

! define cut

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
        Write(message,'(a)') '*** warning - rcut <= Min(domain width) < rlnk = rcut + rpad !!! ***'
        Call error(307,message,.true.)
     Else ! rpad is defined & in 'no strict' mode
        If (cut < rlnk) Then
           If (l_str) Then
              Write(message,'(a)') '*** warning - rcut <= Min(domain width) < rlnk = rcut + rpad !!! ***'
              Call error(307,message,.true.)
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
           If (comm%mxnode > 1) Then
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
  If (comm%idnode == 0) Write(nrite,'(/,1x,2(a,f5.2),a,/)') 'cutoff padding reset from ', rpad, ' Angs to ', cut, ' Angs'
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

Subroutine vnl_set_check(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for (re)setting the conditional VNL checkpoint -
! xbg,ybg,zbg at the end of set_halo_particles
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type ( comms_type ), Intent( InOut ) :: comm
  Logical, Save :: newjob=.true.

  Integer :: fail
  Character( Len = 256 ):: message
  If (.not.llvnl) Return

  If (newjob) Then ! Init set
     newjob = .false.

     fail = 0
     Allocate (xbg(1:mxatms),ybg(1:mxatms),zbg(1:mxatms), Stat = fail)
     If (fail > 0) Then
        Write(message,'(/,1x,a)') 'vnl_set_check allocation failure'
        Call error(0,message)
     End If

! CVNL state and skippage accumulators are initialised in vnl_module
!
!    l_vnl = .true.
!    skipvnl(1) - cycles counter
!    skipvnl(2) - access counter
!    skipvnl(3) - average cycles
!    skipvnl(4) - minimum cycles
!    skipvnl(5) - maximum cycles
  End If

! set tracking point

  xbg(1:nlast)=xxx(1:nlast)
  ybg(1:nlast)=yyy(1:nlast)
  zbg(1:nlast)=zzz(1:nlast)

End Subroutine vnl_set_check


End Module vnl
