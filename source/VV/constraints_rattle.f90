Subroutine constraints_rattle         &
           (mxshak,tolnce,tstep,lcol, &
           lstopt,dxx,dyy,dzz,listot, &
           vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying constraint corrections to the
! velocities of constrained atoms
!
! Note: must be used in conjunction with integration algorithms
!       VV applicable ONLY
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module
  Use config_module, Only : natms,nlast,lsi,lsa,lfrzn,weight
  Use constraints_module

  Implicit None

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce,tstep
  Logical,           Intent( In    ) :: lcol
  Integer,           Intent( In    ) :: lstopt(0:2,1:mxcons)
  Real( Kind = wp ), Intent( In    ) :: dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons)
  Integer,           Intent( In    ) :: listot(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)

  Logical           :: safe
  Integer           :: fail(1:2),i,j,k,icyc
  Real( Kind = wp ) :: amti,amtj,dli,dlj,esig,gamma,gammi,gammj

  Real( Kind = wp ), Dimension( : ), Allocatable :: vxt,vyt,vzt
  Real( Kind = wp ), Dimension( : ), Allocatable :: dt

  fail=0
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms), Stat=fail(1))
  Allocate (dt(1:mxcons),                              Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_rattle allocation failure, node: ', idnode
     Call error(0)
  End If


! initialise dt

  dt(1:ntcons)=0.0_wp

! application of constraint (rattle) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  safe=.false.
  icyc=0
  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! update velocities globally: transport velocity updates of shared atoms to other nodes

     If (lshmv_con) Call update_shared_units(natms,nlast,lsi,lsa,lishp_con,lashp_con,vxx,vyy,vzz)

! initialise velocity correction arrays

     Do i=1,natms
        vxt(i)=0.0_wp
        vyt(i)=0.0_wp
        vzt(i)=0.0_wp
     End Do

! calculate velocity constraint corrections

     esig=0.0_wp
     Do k=1,ntcons
        If (lstopt(0,k) == 0) Then
           i=lstopt(1,k)
           j=lstopt(2,k)

           amti=tstep/weight(i)
           amtj=tstep/weight(j)

! no corrections for frozen atoms

           If (lfrzn(i) /= 0) amti=0.0_wp
           If (lfrzn(j) /= 0) amtj=0.0_wp

! calculate constraint force parameter - gamma

           If (icyc == 1) dt(k) = Sqrt( dxx(k)**2 + dyy(k)**2 + dzz(k)**2 )

           gamma = ( dxx(k)*(vxx(i)-vxx(j)) + dyy(k)*(vyy(i)-vyy(j)) + dzz(k)*(vzz(i)-vzz(j)) ) / dt(k)

           esig=Max(esig,0.5_wp*tstep*Abs(gamma)*dt(k)/prmcon(listcon(0,k)))

           gamma = gamma / ( (amti+amtj) * dt(k) )

! improve approximate constraint velocity and force

           If (i <= natms .and. lfrzn(i) == 0) Then
              gammi =-gamma*amti
              vxt(i)=vxt(i)+dxx(k)*gammi
              vyt(i)=vyt(i)+dyy(k)*gammi
              vzt(i)=vzt(i)+dzz(k)*gammi
           End If

           If (j <= natms .and. lfrzn(j) == 0) Then
              gammj = gamma*amtj
              vxt(j)=vxt(j)+dxx(k)*gammj
              vyt(j)=vyt(j)+dyy(k)*gammj
              vzt(j)=vzt(j)+dzz(k)*gammj
           End If
        End If
     End Do

! global verification of convergence

     safe=(esig < tolnce)
     If (mxnode > 1) Call gcheck(safe,"enforce")

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then

! update velocities locally

        Do k=1,ntcons
           If (lstopt(0,k) == 0) Then
              i=lstopt(1,k)
              j=lstopt(2,k)

              If (i <= natms .and. lfrzn(i) == 0) Then
                 dli = 1.0_wp/Real(listot(i),wp)
                 vxx(i)=vxx(i)+vxt(i)*dli
                 vyy(i)=vyy(i)+vyt(i)*dli
                 vzz(i)=vzz(i)+vzt(i)*dli
              End If

              If (j <= natms .and. lfrzn(j) == 0) Then
                 dlj = 1.0_wp/Real(listot(j),wp)
                 vxx(j)=vxx(j)+vxt(j)*dlj
                 vyy(j)=vyy(j)+vyt(j)*dlj
                 vzz(j)=vzz(j)+vzt(j)*dlj
              End If
           End If
        End Do

     End If
  End Do

  If (.not.safe) Then ! error exit for non-convergence
     Call error(515)
  Else ! Collect per call and per step passage statistics
     passcon(1,1,2)=icyc-1
     passcon(3,1,2)=passcon(2,1,2)*passcon(3,1,2)
     passcon(2,1,2)=passcon(2,1,2)+1
     passcon(3,1,2)=passcon(3,1,2)/passcon(2,1,2)+passcon(1,1,2)/passcon(2,1,2)
     passcon(4,1,2)=Min(passcon(1,1,2),passcon(4,1,2))
     passcon(5,1,2)=Max(passcon(1,1,2),passcon(5,1,2))

     passcon(1,2,2)=passcon(1,2,2)+passcon(1,1,2)
     If (lcol) Then ! Collect
        passcon(3,2,2)=passcon(2,2,2)*passcon(3,2,2)
        passcon(2,2,2)=passcon(2,2,2)+1
        passcon(3,2,2)=passcon(3,2,2)/passcon(2,2,2)+passcon(1,2,2)/passcon(2,2,2)
        passcon(4,2,2)=Min(passcon(1,2,2),passcon(4,2,2))
        passcon(5,2,2)=Max(passcon(1,2,2),passcon(5,2,2))
        passcon(1,2,2)=0.0_wp ! Reset
     End If
     passcon(1,1,2)=0.0_wp ! Reset
  End If

  Deallocate (vxt,vyt,vzt, Stat=fail(1))
  Deallocate (dt,          Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_rattle deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine constraints_rattle
