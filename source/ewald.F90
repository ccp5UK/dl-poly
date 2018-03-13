Module ewald

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring ewald routines arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp

  Implicit None

  Private

  Type, Public :: ewald_type
    Private
    Logical, Public :: l_fce  = .true. , l_cp  = .false. , &
                       lf_fce = .true. , lf_cp = .false.

    Real( Kind = wp ), Public :: engsic = 0.0_wp ,                             &
                          e_rc = 0.0_wp , v_rc = 0.0_wp , s_rc(1:9) = 0.0_wp , &
                          e_fr = 0.0_wp , v_fr = 0.0_wp , s_fr(1:9) = 0.0_wp , &
                          ef_fr= 0.0_wp , vf_fr= 0.0_wp , sf_fr(1:9)= 0.0_wp

    Real( Kind = wp ), Public, Allocatable :: fcx(:),fcy(:),fcz(:)
    Real( Kind = wp ), Public, Allocatable :: ffx(:),ffy(:),ffz(:)

    Logical :: newjob_kall = .true., &
               newjob_kfrz = .true.

    Contains
      Private
      Procedure, Public :: check => ewald_check
      Procedure, Public :: refresh => ewald_refresh
      Final :: ewald_deallocate
  End Type ewald_type

Contains

  Subroutine ewald_allocate_kall_arrays(T)
    Use setup_module, Only : mxatms

    Class(ewald_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%fcx(1:mxatms),T%fcy(1:mxatms),T%fcz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1040)

    T%fcx = 0.0_wp
    T%fcy = 0.0_wp
    T%fcz = 0.0_wp
  End Subroutine ewald_allocate_kall_arrays

  Subroutine ewald_allocate_kfrz_arrays(T)
    Use setup_module, Only : mxatms

    Class(ewald_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%ffx(1:mxatms),T%ffy(1:mxatms),T%ffz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1070)

    T%ffx = 0.0_wp
    T%ffy = 0.0_wp
    T%ffz = 0.0_wp
  End Subroutine ewald_allocate_kfrz_arrays

  Subroutine ewald_check(T,keyens,megfrz,nsteql,nstfce,nstep)
    Class(ewald_type) :: T

    Integer, Intent( In    ) :: keyens,megfrz, &
                                nsteql,nstfce,nstep

! Full frozen-frozen evaluation is TRUE by default
! otherwise a "refresh" is applied.

    If (keyens < 20 .and. megfrz > 0) Then
       T%lf_fce = .false.

       If (T%newjob_kfrz) Then
          T%newjob_kfrz = .false.

! At restart allocate the "refresh" k-space frozen-frozen SPME arrays

          Call ewald_allocate_kfrz_arrays(T)

! Allow copying into these arrays

          T%lf_cp=.true. ! == (keyens < 20 .and. megfrz > 0)

! Force the full frozen force evaluation

          T%lf_fce=.true.
       End If

! Reinitialise

       If (T%lf_fce) Then
         T%ffx = 0.0_wp
         T%ffy = 0.0_wp
         T%ffz = 0.0_wp

         T%ef_fr = 0.0_wp
         T%vf_fr = 0.0_wp
         T%sf_fr = 0.0_wp
       End If
    End If

! Full evaluation is TRUE by default and at any 'nstfce' timestep
! otherwise a "refresh" is applied.

    If (nstfce > 1) Then
       T%l_fce = ( (nstep <  nsteql .and. Mod(nstep,nstfce) == 0) .or. &
                 (nstep >= nsteql .and. Mod(nstep-nsteql,nstfce) == 0) )

       If (T%newjob_kall) Then
          T%newjob_kall = .false.

! At restart allocate the "refresh" k-space all SPME arrays

          Call ewald_allocate_kall_arrays(T)

! Allow copying into these arrays

          T%l_cp=.true. ! == (nstfce > 1)

! Force the full force evaluation

          T%l_fce=.true.
       End If

! Reinitialise

       If (T%l_fce) Then
         T%fcx = 0.0_wp
         T%fcy = 0.0_wp
         T%fcz = 0.0_wp

         T%e_rc = 0.0_wp
         T%v_rc = 0.0_wp
         T%s_rc = 0.0_wp

         T%e_fr = 0.0_wp
         T%v_fr = 0.0_wp
         T%s_fr = 0.0_wp
       End If
    End If

  End Subroutine ewald_check

  Subroutine ewald_refresh(T,engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stress)
    Use configuration, Only : natms,fxx,fyy,fzz

    Class( ewald_type )                :: T
    Real( Kind = wp ), Intent( InOut ) :: engcpe_rc,vircpe_rc, &
                                          engcpe_fr,vircpe_fr, &
                                          stress(1:9)

    Integer :: i

    Do i=1,natms
       fxx(i)=fxx(i)+T%fcx(i)
       fyy(i)=fyy(i)+T%fcy(i)
       fzz(i)=fzz(i)+T%fcz(i)
    End Do

    engcpe_rc=T%e_rc
    vircpe_rc=T%v_rc
    stress=stress+T%s_rc

    engcpe_fr=T%e_fr
    vircpe_fr=T%v_fr
    stress=stress+T%s_fr
  End Subroutine ewald_refresh

  Subroutine ewald_deallocate(T)
  !> Deallocate ewald_type arrays
    Type(ewald_type) :: T

    If (Allocated(T%fcx)) Then
      Deallocate(T%fcx)
    End If
    If (Allocated(T%fcy)) Then
      Deallocate(T%fcy)
    End If
    If (Allocated(T%fcz)) Then
      Deallocate(T%fcz)
    End If

    If (Allocated(T%ffx)) Then
      Deallocate(T%ffx)
    End If
    If (Allocated(T%ffy)) Then
      Deallocate(T%ffy)
    End If
    If (Allocated(T%ffz)) Then
      Deallocate(T%ffz)
    End If
  End Subroutine ewald_deallocate

End Module ewald
