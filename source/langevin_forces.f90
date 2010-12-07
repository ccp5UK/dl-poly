Subroutine langevin_forces(temp,tstep,chi,fxr,fyr,fzr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to generate Langevin random forces consistent
! with the target temperature and the Langevin thermostat friction
!
! Note: (1) This algorithm breaks reversibility due to the random
!           generastion of forces.
!       (2) Random forces do not contribute to the stress and virial
!           of the system they are accounted by the system pressure.
!       (3) Random forces do not apply to frozen and massless particles
!           as well as shells.
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsum
  Use setup_module,      Only : nrite,boltz,mxatms,zero_plus
  Use config_module,     Only : natms,lfrzn,ltg,weight
  Use core_shell_module, Only : ntshl,listshl

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: temp,tstep,chi
  Real( Kind = wp ), Intent(   Out ) :: fxr(1:mxatms),fyr(1:mxatms),fzr(1:mxatms)

  Integer                 :: fail,i
  Real( Kind = wp )       :: vom(0:3),scale,tmp

  Integer, Allocatable :: q(:)


  fail = 0
  Allocate (q(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'langevin_forces allocation failure, node: ', idnode
     Call error(0)
  End If

! Here we become node-dependent - i.e. pseudo-randomness depends
! on the DD mapping which depends on the number of nodes and system size
!
! Get gaussian distribution (unit variance)

  Call gauss(natms,fxr,fyr,fzr)

! Get scaler to target variance*Sqrt(weight)

  scale = Sqrt(2.0_wp * chi * boltz * temp / tstep)

! Make variance = target variance and nulify the rest and
! calculate force COM correction

  vom = 0.0_wp
  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. All(listshl(2,1:ntshl) /= ltg(i))) Then
        q(i) = 1

        tmp = scale*Sqrt(weight(i))

        fxr(i) = fxr(i)*tmp
        fyr(i) = fyr(i)*tmp
        fzr(i) = fzr(i)*tmp

        vom(0) = vom(0) + 1.0_wp
        vom(1) = vom(1) + fxr(i)
        vom(2) = vom(2) + fyr(i)
        vom(3) = vom(3) + fzr(i)
     Else
        q(i) = 0

        fxr(i) = 0.0_wp
        fyr(i) = 0.0_wp
        fzr(i) = 0.0_wp
     End If
  End Do
  If (mxnode > 1) Call gsum(vom)
  If (vom(0) > zero_plus) vom = vom / vom(0)

! Remove COM force

  Do i=1,natms
     If (q(i) == 1) Then
        fxr(i) = fxr(i) - vom(1)
        fyr(i) = fyr(i) - vom(2)
        fzr(i) = fzr(i) - vom(3)
     End If
  End Do

  Deallocate (q, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'langevin_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine langevin_forces
