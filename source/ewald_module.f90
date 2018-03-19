Module ewald_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring ewald routines arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: l_fce  = .true. , l_cp  = .false. , &
                                          lf_fce = .true. , lf_cp = .false.

  Real( Kind = wp ),              Save :: engsic = 0.0_wp ,                                    &
                                          e_rc = 0.0_wp , v_rc = 0.0_wp , s_rc(1:9) = 0.0_wp , &
                                          e_fr = 0.0_wp , v_fr = 0.0_wp , s_fr(1:9) = 0.0_wp , &
                                          ef_fr= 0.0_wp , vf_fr= 0.0_wp , sf_fr(1:9)= 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: fcx(:),fcy(:),fcz(:)
  Real( Kind = wp ), Allocatable, Save :: ffx(:),ffy(:),ffz(:)

  Public :: ewald_refresh

Contains

  Subroutine ewald_allocate_kall_arrays

    Use setup_module, Only : mxatms

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (fcx(1:mxatms),fcy(1:mxatms),fcz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1040)

    fcx = 0.0_wp ; fcy = 0.0_wp ; fcz = 0.0_wp

  End Subroutine ewald_allocate_kall_arrays

  Subroutine ewald_allocate_kfrz_arrays

    Use setup_module, Only : mxatms

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (ffx(1:mxatms),ffy(1:mxatms),ffz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1070)

    ffx = 0.0_wp ; ffy = 0.0_wp ; ffz = 0.0_wp

  End Subroutine ewald_allocate_kfrz_arrays

  Subroutine ewald_check(keyens,megfrz,nsteql,nstfce,nstep)

    Implicit None

    Integer, Intent( In    ) :: keyens,megfrz, &
                                nsteql,nstfce,nstep

    Logical, Save :: newjob_kall = .true., &
                     newjob_kfrz = .true.

! Full frozen-frozen evaluation is TRUE by default
! otherwise a "refresh" is applied.

    If (keyens < 20 .and. megfrz > 0) Then
       lf_fce = .false.

       If (newjob_kfrz) Then
          newjob_kfrz = .false.

! At restart allocate the "refresh" k-space frozen-frozen SPME arrays

          Call ewald_allocate_kfrz_arrays()

! Allow copying into these arrays

          lf_cp=.true. ! == (keyens < 20 .and. megfrz > 0)

! Force the full frozen force evaluation

          lf_fce=.true.
       End If

! Reinitialise

       If (lf_fce) Then
          ffx = 0.0_wp ; ffy = 0.0_wp ; ffz = 0.0_wp

          ef_fr = 0.0_wp ; vf_fr = 0.0_wp ; sf_fr = 0.0_wp
       End If
    End If

! Full evaluation is TRUE by default and at any 'nstfce' timestep
! otherwise a "refresh" is applied.

    If (nstfce > 1) Then
       l_fce = ( (nstep <  nsteql .and. Mod(nstep,nstfce) == 0) .or. &
                 (nstep >= nsteql .and. Mod(nstep-nsteql,nstfce) == 0) )

       If (newjob_kall) Then
          newjob_kall = .false.

! At restart allocate the "refresh" k-space all SPME arrays

          Call ewald_allocate_kall_arrays()

! Allow copying into these arrays

          l_cp=.true. ! == (nstfce > 1)

! Force the full force evaluation

          l_fce=.true.
       End If

! Reinitialise

       If (l_fce) Then
          fcx = 0.0_wp ; fcy = 0.0_wp ; fcz = 0.0_wp

          e_rc = 0.0_wp ; v_rc = 0.0_wp ; s_rc = 0.0_wp
          e_fr = 0.0_wp ; v_fr = 0.0_wp ; s_fr = 0.0_wp
       End If
    End If

  End Subroutine ewald_check

  Subroutine ewald_refresh(engcpe_rc,vircpe_rc,engcpe_fr,vircpe_fr,stress)

    Use config_module, Only : natms,fxx,fyy,fzz

    Implicit None

    Real( Kind = wp ), Intent( InOut ) :: engcpe_rc,vircpe_rc, &
                                          engcpe_fr,vircpe_fr, &
                                          stress(1:9)

    Integer :: i

    Do i=1,natms
       fxx(i)=fxx(i)+fcx(i)
       fyy(i)=fyy(i)+fcy(i)
       fzz(i)=fzz(i)+fcz(i)
    End Do

    engcpe_rc=e_rc ; vircpe_rc=v_rc ; stress=stress+s_rc
    engcpe_fr=e_fr ; vircpe_fr=v_fr ; stress=stress+s_fr

  End Subroutine ewald_refresh

End Module ewald_module
