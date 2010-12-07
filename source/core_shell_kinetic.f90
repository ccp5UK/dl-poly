Subroutine core_shell_kinetic(shlke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the internal kinetic energy of
! core-shell units in the shell polarisation model
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : mxnode,gsum
  Use config_module,     Only : natms,nlast,lsi,lsa,weight,vxx,vyy,vzz
  Use core_shell_module, Only : ntshl,listshl,lshmv_shl,lishp_shl,lashp_shl

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: shlke

  Integer           :: i,j,k,local_index
  Real( Kind = wp ) :: rmu,rvx,rvy,rvz

! gather velocities of shared particles

  If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)

! initialise energy

  shlke=0.0_wp

! loop over all specified core-shell pairs

  Do k=1,ntshl

! indices of atoms involved

     i=local_index(listshl(1,k),nlast,lsi,lsa)
     j=local_index(listshl(2,k),nlast,lsi,lsa)

! for all native and natively shared core-shell units

     If ((i > 0 .and. j > 0) .and. (i <= natms .or. j <= natms)) Then

! calculate reduced mass

        rmu=(weight(i)*weight(j))/(weight(i)+weight(j))

! frozen particles' velocities are zero
! calculate constraint vector normal

        rvx=vxx(j)-vxx(i)
        rvy=vyy(j)-vyy(i)
        rvz=vzz(j)-vzz(i)

! calculate core-shell internal kinetic energy

        If (i <= natms .and. j <= natms) Then

! for native core-shell units - full contribution

           shlke=shlke+0.5_wp*rmu*(rvx*rvx+rvy*rvy+rvz*rvz)

        Else

! If ( (i <= natms .and. j > natms) .or. (j <= natms .and. i > natms) ) Then
! for shared core-shell units - halved contribution

           shlke=shlke+0.25_wp*rmu*(rvx*rvx+rvy*rvy+rvz*rvz)

        End If

     End If

  End Do

! global sum of core-shell internal kinetic energy

  If (mxnode > 1) Call gsum(shlke)

End Subroutine core_shell_kinetic
