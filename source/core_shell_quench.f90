Subroutine core_shell_quench(safe,temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for quenching the internal bond energies of ions
! defined by shell model
!
! copyright - daresbury laboratory
! author    - w.smith august 1999
! amended   - i.t.todorov november 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : mxnode,gcheck
  Use setup_module
  Use config_module,     Only : natms,nlast,lsi,lsa,lfrzn,weight,vxx,vyy,vzz
  Use core_shell_module, Only : ntshl,listshl,lshmv_shl,lishp_shl,lashp_shl

  Implicit None

  Logical,           Intent(   Out ) :: safe
  Real( Kind = wp ), Intent( In    ) :: temp

  Logical,          :: safek
  Integer           :: ia,ib,k,local_index
  Real( Kind = wp ) :: dvx,dvy,dvz,pke,rmu,scl,tke,tmx,tmy,tmz

! Initialise safe flag

  safe=.true

! gather velocities of shared particles

  If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)

! permitted core-shell internal kinetic energy

  pke=boltz*temp*1.0e-4_wp

! amend core and shell velocities preserving total momentum

  Do k=1,ntshl
     ia=local_index(listshl(1,k),nlast,lsi,lsa)
     ib=local_index(listshl(2,k),nlast,lsi,lsa)

     If ((ia > 0 .and. ib > 0) .and. (ia <= natms .or. ib <= natms)) Then
        rmu=(weight(ia)*weight(ib))/(weight(ia)+weight(ib))

! frozen atoms have zero velocity (no worries)

        dvx=vxx(ib)-vxx(ia)
        dvy=vyy(ib)-vyy(ia)
        dvz=vzz(ib)-vzz(ia)

        tke=rmu*(dvx*dvx+dvy*dvy+dvz*dvz)

        safek = .not.( tke > pke .and. (tke-pke)/pke > 1.0e-3_wp )
        If (.not.safek) Then
           scl=Sqrt(pke/tke)

           tmx=weight(ia)*vxx(ia)+weight(ib)*vxx(ib)
           tmy=weight(ia)*vyy(ia)+weight(ib)*vyy(ib)
           tmz=weight(ia)*vzz(ia)+weight(ib)*vzz(ib)

! no corrections for frozen cores

           If (lfrzn(ia) == 0) Then
              vxx(ia)=tmx/(weight(ia)+weight(ib))-scl*rmu*dvx/weight(ia)
              vyy(ia)=tmy/(weight(ia)+weight(ib))-scl*rmu*dvy/weight(ia)
              vzz(ia)=tmz/(weight(ia)+weight(ib))-scl*rmu*dvz/weight(ia)

              vxx(ib)=tmx/(weight(ia)+weight(ib))+scl*rmu*dvx/weight(ib)
              vyy(ib)=tmy/(weight(ia)+weight(ib))+scl*rmu*dvy/weight(ib)
              vzz(ib)=tmz/(weight(ia)+weight(ib))+scl*rmu*dvz/weight(ib)
           Else
              vxx(ib)=tmx/(weight(ia)+weight(ib))+scl*rmu*dvx/weight(ib)
              vyy(ib)=tmy/(weight(ia)+weight(ib))+scl*rmu*dvy/weight(ib)
              vzz(ib)=tmz/(weight(ia)+weight(ib))+scl*rmu*dvz/weight(ib)
           End If

           safe=.false.
        End If
     End If
  End Do

  If (mxnode > 1) Call gcheck(safe)

End Subroutine core_shell_quench
