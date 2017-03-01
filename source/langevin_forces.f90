Subroutine langevin_forces(nstep,temp,tstep,chi,fxr,fyr,fzr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to generate Langevin random forces consistent
! with the target temperature and the Langevin thermostat friction
!
! Note: (1) Random forces do not contribute to the stress and virial
!           of the system they are accounted by the system pressure.
!       (2) Random forces do not apply to frozen and massless particles
!           as well as shells.
!       (3) Random forces are scaled according to local electronic
!           temperature (and dynamically determined friction parameter)
!           if two-temperature model is in use.
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
! contrib   - g.khara & m.a.seaton february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,      Only : boltz,mxatms
  Use config_module,     Only : natms,ltg,lfrzn,weight,xxx,yyy,zzz
  Use core_shell_module, Only : legshl
  Use ttm_module,        Only : eltemp,zerocell,ntcell,delx,dely,delz,gvar,l_ttm,nstepcpl
  Use ttm_utils,         Only : Gep

  Implicit None

  Integer          , Intent( In    ) :: nstep
  Real( Kind = wp ), Intent( In    ) :: temp,tstep,chi
  Real( Kind = wp ), Intent(   Out ) :: fxr(1:mxatms),fyr(1:mxatms),fzr(1:mxatms)

  Integer           :: i,ia,ja,ka,ijk
  Real( Kind = wp ) :: scale,tmp

  If (l_ttm .and. nstep>nstepcpl) Then

! Rescale chi for average electronic temperature if using
! homogeneous electron-phonon coupling

    Select Case (gvar)
    Case (0,1)
    ! constant electron-phonon chi parameter and homogeneous
    ! e-p coupling cases
      scale = Sqrt(2.0_wp * chi * boltz / tstep)
      Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. legshl(0,i) >= 0) Then
          Call box_mueller_saru3(ltg(i),nstep,fxr(i),fyr(i),fzr(i))
          ia = Floor((xxx(i)+zerocell(1))/delx) + 1
          ja = Floor((yyy(i)+zerocell(2))/dely) + 1
          ka = Floor((zzz(i)+zerocell(3))/delz) + 1
          ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
          tmp = scale*Sqrt(eltemp(ijk,0,0,0)*weight(i))

          fxr(i) = fxr(i)*tmp
          fyr(i) = fyr(i)*tmp
          fzr(i) = fzr(i)*tmp
       Else
          fxr(i) = 0.0_wp
          fyr(i) = 0.0_wp
          fzr(i) = 0.0_wp
        End If
      End Do

    Case (2)
    ! heterogeneous electron-phonon coupling case: calculate individual
    ! chi value for each ionic temperature voxel (ignore input value)
      scale = Sqrt(2.0_wp * boltz / tstep)
      Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. legshl(0,i) >= 0) Then
          Call box_mueller_saru3(ltg(i),nstep,fxr(i),fyr(i),fzr(i))
          ia = Floor((xxx(i)+zerocell(1))/delx) + 1
          ja = Floor((yyy(i)+zerocell(2))/dely) + 1
          ka = Floor((zzz(i)+zerocell(3))/delz) + 1
          ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
          tmp = scale*Sqrt(Gep(eltemp(ijk,0,0,0))*eltemp(ijk,0,0,0)*weight(i))

          fxr(i) = fxr(i)*tmp
          fyr(i) = fyr(i)*tmp
          fzr(i) = fzr(i)*tmp
       Else
          fxr(i) = 0.0_wp
          fyr(i) = 0.0_wp
          fzr(i) = 0.0_wp
        End If

      End Do

    End Select

  Else

! Get scaler to target variance*Sqrt(weight)

    scale = Sqrt(2.0_wp * chi * boltz * temp / tstep)

! Make variance = target variance and nullify the rest

    Do i=1,natms
       If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. legshl(0,i) >= 0) Then
          Call box_mueller_saru3(ltg(i),nstep,fxr(i),fyr(i),fzr(i))

          tmp = scale*Sqrt(weight(i))

          fxr(i) = fxr(i)*tmp
          fyr(i) = fyr(i)*tmp
          fzr(i) = fzr(i)*tmp
       Else
          fxr(i) = 0.0_wp
          fyr(i) = 0.0_wp
          fzr(i) = 0.0_wp
       End If
    End Do

  End If

End Subroutine langevin_forces
