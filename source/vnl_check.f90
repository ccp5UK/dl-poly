Subroutine vnl_check(l_str,m_rgd,rcut,rpad,rlnk,width)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine implementing the VNL conditional update checks
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2016
! contrib   - i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gmax
  Use setup_module,   Only : nrite,mxspl,mxatms,mxatdm
  Use domains_module, Only : r_nprx,r_npry,r_nprz
  Use config_module,  Only : imcon,cell,natms,nlast,xxx,yyy,zzz
  Use vnl_module

  Implicit None

  Logical,           Intent ( In    ) :: l_str
  Integer,           Intent ( In    ) :: m_rgd
  Real( Kind = wp ), Intent ( In    ) :: rcut
  Real( Kind = wp ), Intent ( InOut ) :: rpad,rlnk,width

  Logical,     Save :: newjob=.true.,newstart=.true.
  Integer,     Save :: mxsize,alist
  Integer           :: fail,ilx,ily,ilz
  Real( Kind = wp ) :: max_disp,cut,test,tol,celprp(1:10)

  Real( Kind = wp ), Allocatable :: x(:),y(:),z(:)

  If (.not.llvnl) Return

  If (newjob) Then ! Init set

     newjob = .false.

! In principle the upper bound of the tracking arrays should be 'mxatdm'.
! However, when rigid bodies are present in the system (m_rgd > 0) calls
! to update_shared_units and rigid_bodies_coms are made in xscale and
! the upper bound should then be increased to 'mxatms'.  Similarly, the
! the assignment range changes from 'natms' to 'nlast'

     If (m_rgd == 0) Then
        mxsize=mxatdm
        alist =natms
     Else
        mxsize=mxatms
        alist =nlast
     End If

     fail = 0
     Allocate (xbg(1:mxsize),ybg(1:mxsize),zbg(1:mxsize), Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'vnl_check_0 allocation failure, node: ', idnode
        Call error(0)
     End If

     xbg(1:alist)=xxx(1:alist)
     ybg(1:alist)=yyy(1:alist)
     zbg(1:alist)=zzz(1:alist)

! CVNL state and skippage accumulators are initialised in vnl_module
!
! l_vnl = .true.
! skipvnl(1) - cycles counter
! skipvnl(2) - access counter
! skipvnl(3) - average cycles
! skipvnl(4) - minimum cycles
! skipvnl(5) - maximum cycles

  Else ! Checks

     fail = 0
     Allocate (x(1:mxatdm),y(1:mxatdm),z(1:mxatdm), Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'vnl_check allocation failure, node: ', idnode
        Call error(0)
     End If

     x(1:natms) = xxx(1:natms) - xbg(1:natms)
     y(1:natms) = yyy(1:natms) - ybg(1:natms)
     z(1:natms) = zzz(1:natms) - zbg(1:natms)

     Call images(imcon,cell,natms,x,y,z)

     max_disp = Merge( Sqrt( Maxval( x(1:natms)**2 + y(1:natms)**2 + z(1:natms)**2 ) ) , 0.0_wp , natms > 0 )

     Deallocate (x,y,z, Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'vnl_check deallocation failure, node: ', idnode
        Call error(0)
     End If

     If (mxnode > 1) Call gmax(max_disp)

     l_vnl = ( max_disp >= 0.5_wp*rpad )

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
                 cut = Min( 0.85_wp * ( Min ( r_nprx * celprp(7) / Real(ilx,wp) , &
                                              r_npry * celprp(8) / Real(ily,wp) , &
                                              r_nprz * celprp(9) / Real(ilz,wp) ) &
                                        - rcut - 1.0e-6_wp ) , test * rcut )
              Else ! catch & handle exception
                 cut = 0.85_wp * (0.5_wp*width - rcut - 1.0e-6_wp)
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

     If (l_vnl) Then

! Deal with skipping statistics

        skipvnl(3)=skipvnl(2)*skipvnl(3)
        skipvnl(2)=skipvnl(2)+1.0_wp
        skipvnl(3)=skipvnl(3)/skipvnl(2)+skipvnl(1)/skipvnl(2)
        If (.not.newstart) Then ! avoid first compulsory force evaluation
           skipvnl(4)=Min(skipvnl(1),skipvnl(4))
        Else
           newstart=.false.
        End If
        skipvnl(5)=Max(skipvnl(1),skipvnl(5))

! Reset

        skipvnl(1) = 0.0_wp

        xbg(1:alist)=xxx(1:alist)
        ybg(1:alist)=yyy(1:alist)
        zbg(1:alist)=zzz(1:alist)

     Else

! Enjoy telephoning

        skipvnl(1) = skipvnl(1) + 1.0_wp

     End If

  End If

End Subroutine vnl_check
