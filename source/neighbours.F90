!> Module containing neighbour list routines and variables
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J.Madge June 2018
Module neighbours
  Use kinds, Only : wp
  Use comms,   Only : comms_type,gcheck
  Use setup,   Only : nrite,mxspl,mxatms
  Use domains, Only : r_nprx,r_npry,r_nprz
  Use configuration,  Only : imcon,cell,natms,nlast,list, &
                             xxx,yyy,zzz
  Use errors_warnings, Only : error,info
  Use numerics, Only : dcell,images
  Implicit None

  Private

  !> Type containing neighbours data
  Type, Public :: neighbours_type
    Private

    ! Verlet neighbour list data
    !> Update cells flag
    Logical, Public :: update = .true.
    !> Unconditional update flag
    Logical, Public :: unconditional_update = .false.

    !> Skips, elements are as follows
    !>
    !> - 1 skips counter
    !> - 2 access counter
    !> - 3 average skips
    !> - 4 minimum skips ~Huge(1)
    !> - 5 maximum skips
    Real( Kind = wp ), Public :: skip(1:5) = [0.0_wp,0.0_wp,0.0_wp, &
                                              999999999.0_wp,0.0_wp]

    !> Tracking points for Verlet neighbour list
    Real( Kind = wp ), Allocatable, Public :: xbg(:),ybg(:),zbg(:)

    !> Largest vdw cutoff, defines Verlet neighbour list radius
    Real ( Kind = wp ), Public :: cutoff
    !> Padding around cutoff
    Real ( Kind = wp ), Public :: padding
    !> Actual Verlet neighbour list cutoff (cutoff+padding)
    Real ( Kind = wp ), Public :: cutoff_extended
  Contains
    Private

    Final :: cleanup
  End Type neighbours_type

  Public :: vnl_check,vnl_set_check

Contains

  !> Perform Verlet neighbour list update check
  !>
  !> Copyright - Daresbury Laboratory
  !>
  !> Author    - I.T.Todorov january 2017
  !>
  !> Contrib   - I.J.Bush february 2014
  Subroutine vnl_check(l_str,width,neigh,comm)
    Logical,           Intent ( In    ) :: l_str
    Real( Kind = wp ), Intent ( InOut ) :: width
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( comms_type ), Intent ( InOut ) :: comm

    Logical, Save :: newstart=.true.

    Logical           :: safe
    Integer           :: fail,ilx,ily,ilz,i,ii,j
    Real( Kind = wp ) :: cut,test,tol,celprp(1:10),xx(1),yy(1),zz(1)

    Real( Kind = wp ), Allocatable :: x(:),y(:),z(:),r(:)

    Character( Len = 256 ) :: message
    If (.not.neigh%unconditional_update) Return

    ! Checks
    fail = 0
    Allocate (x(1:mxatms),y(1:mxatms),z(1:mxatms),r(1:mxatms), Stat = fail)
    If (fail > 0) Then
      Write(message,'(a)') 'vnl_check allocation failure'
      Call error(0,message)
    End If

    x(1:nlast) = xxx(1:nlast) - neigh%xbg(1:nlast)
    y(1:nlast) = yyy(1:nlast) - neigh%ybg(1:nlast)
    z(1:nlast) = zzz(1:nlast) - neigh%zbg(1:nlast)

    Call images(imcon,cell,nlast,x,y,z)

    If (nlast > 0) r(1:nlast) = Sqrt(x(1:nlast)**2 + y(1:nlast)**2 + z(1:nlast)**2)

    safe=.true.
    Q:Do i=1,natms
      If (r(i) > 0.5_wp*neigh%padding) Then
        Do ii=1,list(-2,i)
          j=list(ii,i) ! may be in the halo!!!
          If (r(j) > 0.5_wp*neigh%padding) Then
            xx(1)=x(i)-x(j)
            yy(1)=y(i)-y(j)
            zz(1)=z(i)-z(j)
            cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
            If (cut > neigh%padding) Then
              xx(1)=xxx(i)-xxx(j)
              yy(1)=yyy(i)-yyy(j)
              zz(1)=zzz(i)-zzz(j)
              cut=Sqrt((xx(1))**2+yy(1)**2+zz(1)**2)
              If (cut > neigh%cutoff_extended) Then
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
      Write(message,'(a)') 'vnl_check deallocation failure'
      Call error(0,message)
    End If

    Call gcheck(comm,safe,"enforce")
    neigh%update=.not.safe

    ! Get the dimensional properties of the MD cell

    Call dcell(cell,celprp) ! cut=neigh%cutoff_extended+1.0e-6_wp in link_cell_pairs
    width=Min(celprp(7),celprp(8),celprp(9))

    ! define cut

    cut=neigh%cutoff_extended+1.0e-6_wp

    ! calculate link cell dimensions per node

    ilx=Int(r_nprx*celprp(7)/cut)
    ily=Int(r_npry*celprp(8)/cut)
    ilz=Int(r_nprz*celprp(9)/cut)

    tol=Min(0.05_wp,0.005_wp*neigh%cutoff)                                        ! tolerance
    test = 0.02_wp * Merge( 1.0_wp, 2.0_wp, mxspl > 0)                    ! 2% (w/ SPME) or 4% (w/o SPME)
    cut=Min(r_nprx*celprp(7),r_npry*celprp(8),r_nprz*celprp(9))-1.0e-6_wp ! domain size

    If (ilx*ily*ilz == 0) Then
      If (cut < neigh%cutoff) Then
        Write(message,'(a)') 'neigh%cutoff <= Min(domain width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
        Call error(307,message,.true.)
      Else ! neigh%padding is defined & in 'no strict' mode
        If (cut < neigh%cutoff_extended) Then
          If (l_str) Then
            Write(message,'(a)') 'neigh%cutoff <= Min(domain width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
            Call error(307,message,.true.)
          Else
            If (cut >= neigh%cutoff) Then ! Re-set neigh%padding with some slack
              neigh%padding = Min( 0.95_wp * (cut - neigh%cutoff) , test * neigh%cutoff)
              neigh%padding = Real( Int( 100.0_wp * neigh%padding , wp ) ) / 100.0_wp
              If (neigh%padding < tol) neigh%padding = 0.0_wp ! Don't bother
              neigh%cutoff_extended = neigh%cutoff + neigh%padding
              neigh%update=.true.
            End If
          End If
        End If
      End If
    Else ! push the limits when up for update in a 'no strict' regime
      If (neigh%update .and. (.not.l_str)) Then ! Try to re-set neigh%padding with some slack
        If (Int(Real(Min(ilx,ily,ilz),wp)/(1.0_wp+test)) >= 2) Then
          cut = test * neigh%cutoff
        Else
          If (comm%mxnode > 1) Then
            cut = Min( 0.95_wp * ( Min ( r_nprx * celprp(7) / Real(ilx,wp) , &
              r_npry * celprp(8) / Real(ily,wp) , &
              r_nprz * celprp(9) / Real(ilz,wp) ) &
              - neigh%cutoff - 1.0e-6_wp ) , test * neigh%cutoff )
          Else ! catch & handle exception
            cut = 0.95_wp * (0.5_wp*width - neigh%cutoff - 1.0e-6_wp)
          End If
        End If
        cut = Real( Int( 100.0_wp * cut ) , wp ) / 100.0_wp
        If ((.not.(cut < tol)) .and. cut-neigh%padding > 0.005_wp) Then ! Do bother
          Write(message,'(2(a,f5.2),a)') 'cutoff padding reset from ', neigh%padding, &
            ' Angs to ', cut, ' Angs'
          Call info(message,.true.)
          neigh%padding = cut
          neigh%cutoff_extended = neigh%cutoff + neigh%padding
        End If
      End If
    End If

    If (neigh%update) Then ! Deal with skipping statistics
      neigh%skip(3)=neigh%skip(2)*neigh%skip(3)
      neigh%skip(2)=neigh%skip(2)+1.0_wp
      neigh%skip(3)=neigh%skip(3)/neigh%skip(2)+neigh%skip(1)/neigh%skip(2)
      If (.not.newstart) Then ! avoid first compulsory force evaluation
        neigh%skip(4)=Min(neigh%skip(1),neigh%skip(4))
      Else
        newstart=.false.
      End If
      neigh%skip(5)=Max(neigh%skip(1),neigh%skip(5))

      neigh%skip(1) = 0.0_wp              ! Reset here, checkpoit set by vnl_set_check in set_halo_particles
    Else            ! Enjoy telephoning
      neigh%skip(1) = neigh%skip(1) + 1.0_wp ! Increment, telephony done for xxx,yyy,zzz in set_halo_positions
    End If
  End Subroutine vnl_check

  !> (Re)set the conditional Verlet neighbour list checkpoint -
  !> neigh%xbg,neigh%ybg,neigh%zbg at the end of set_halo_particles
  !>
  !> Copyright - Daresbury Laboratory
  !>
  !> Author    - I.T.Todorov january 2017
  Subroutine vnl_set_check(neigh,comm)
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type ( comms_type ), Intent( InOut ) :: comm
    Logical, Save :: newjob=.true.

    Integer :: fail
    Character( Len = 256 ):: message
    If (.not.neigh%unconditional_update) Return

    If (newjob) Then ! Init set
      newjob = .false.

      fail = 0
      Allocate (neigh%xbg(1:mxatms),neigh%ybg(1:mxatms),neigh%zbg(1:mxatms), Stat = fail)
      If (fail > 0) Then
        Write(message,'(a)') 'vnl_set_check allocation failure'
        Call error(0,message)
      End If

      ! CVNL state and skippage accumulators are initialised in vnl_module
      !
      !    neigh%update = .true.
      !    neigh%skip(1) - cycles counter
      !    neigh%skip(2) - access counter
      !    neigh%skip(3) - average cycles
      !    neigh%skip(4) - minimum cycles
      !    neigh%skip(5) - maximum cycles
    End If

    ! set tracking point
    neigh%xbg(1:nlast)=xxx(1:nlast)
    neigh%ybg(1:nlast)=yyy(1:nlast)
    neigh%zbg(1:nlast)=zzz(1:nlast)
  End Subroutine vnl_set_check

  !> Dealloate neighbours type
  Subroutine cleanup(neigh)
    Type( neighbours_type ) :: neigh

    If (Allocated(neigh%xbg)) Then
      Deallocate(neigh%xbg)
    End If
    If (Allocated(neigh%ybg)) Then
      Deallocate(neigh%ybg)
    End If
    If (Allocated(neigh%zbg)) Then
      Deallocate(neigh%zbg)
    End If
  End Subroutine cleanup
End Module neighbours
