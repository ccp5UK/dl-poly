Module ewald_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring ewald routines arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: l_fce = .true. , &
                                          l_cp  = .false.

  Real( Kind = wp ),              Save :: e_rc = 0.0_wp , v_rc = 0.0_wp , s_rc(1:9) = 0.0_wp , &
                                          e_ex = 0.0_wp , v_ex = 0.0_wp , s_ex(1:9) = 0.0_wp , &
                                          e_fr = 0.0_wp , v_fr = 0.0_wp , s_fr(1:9) = 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: fcx(:),fcy(:),fcz(:)

  Public :: ewald_allocate_arrays,ewald_refresh

Contains

  Subroutine ewald_allocate_arrays

    Use setup_module, Only : mxatms

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (fcx(1:mxatms),fcy(1:mxatms),fcz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1040)

    fcx = 0.0_wp ; fcy = 0.0_wp ; fcz = 0.0_wp

  End Subroutine ewald_allocate_arrays

  Subroutine ewald_check(nstep,nsteql,nstfce)

    Implicit None

    Integer, Intent( In    ) :: nstep,nsteql,nstfce

    Logical, Save :: newjob = .true.

! Full evaluation is TRUE by default and at any 'nstfce' timestep
! otherwise a "refresh" is applied.

    If (nstfce <= 1) Return

    l_fce = ( (nstep <  nsteql .and. Mod(nstep,nstfce) == 0) .or. &
              (nstep >= nsteql .and. Mod(nstep-nsteql,nstfce) == 0) )

    If (newjob) Then
       newjob = .false.

! At restart allocate the "refresh" k-space SPME arrays

       Call ewald_allocate_arrays()

! Allow copying into these arrrays

       l_cp=.true. ! == (nstfce > 1)

! Force the full force evaluation

       l_fce=.true.
    End If

! Reinitialise

    If (l_cp) Then
       fcx = 0.0_wp ; fcy = 0.0_wp ; fcx = 0.0_wp

       e_rc = 0.0_wp ; v_rc = 0.0_wp ; s_rc = 0.0_wp
       e_ex = 0.0_wp ; v_ex = 0.0_wp ; s_ex = 0.0_wp
       e_fr = 0.0_wp ; v_fr = 0.0_wp ; s_fr = 0.0_wp
    End If

  End Subroutine ewald_check

  Subroutine ewald_refresh(engcpe_rc,vircpe_rc,engcpe_ex,vircpe_ex,engcpe_fr,vircpe_fr,stress)

    Use config_module, Only : natms,fxx,fyy,fzz

    Implicit None

    Real( Kind = wp ), Intent( InOut ) :: engcpe_rc,vircpe_rc, &
                                          engcpe_ex,vircpe_ex, &
                                          engcpe_fr,vircpe_fr, &
                                          stress(1:9)

    Integer :: i

    Do i=1,natms
       fxx(i)=fxx(i)+fcx(i)
       fyy(i)=fyy(i)+fcy(i)
       fzz(i)=fzz(i)+fcz(i)
    End Do

    engcpe_rc=e_rc ; vircpe_rc=v_rc ; stress=stress+s_rc
    engcpe_ex=e_ex ; vircpe_ex=v_ex ; stress=stress+s_ex
    engcpe_fr=e_fr ; vircpe_fr=v_fr ; stress=stress+s_fr

  End Subroutine ewald_refresh

End Module ewald_module
