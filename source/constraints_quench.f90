Subroutine constraints_quench(imcon,mxshak,tolnce)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for quenching the internal bond energies in the
! initial structure of a molecule defined by constraints
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module
  Use config_module
  Use constraints_module

  Implicit None

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce

  Logical           :: safe
  Integer           :: fail(1:4),i,j,k,icyc
  Real( Kind = wp ) :: amti,amtj,dlj,dli,esig,gamma,gammi,gammj

  Logical,           Allocatable :: lstitr(:)
  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:),dt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)

  fail=0
  Allocate (lstitr(1:mxatms),                                       Stat=fail(1))
  Allocate (lstopt(1:2,1:mxcons),listot(1:mxatms),                  Stat=fail(2))
  Allocate (dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons),dt(1:mxcons), Stat=fail(3))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),              Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_quench allocation failure, node: ', idnode
     Call error(0)
  End If

! gather velocities of shared atoms

  If (lshmv_con) Call update_shared_units(natms,nlast,lsi,lsa,lishp_con,lashp_con,vxx,vyy,vzz)

! construct current constrained bond vectors and listot array (shared
! constraint atoms) for iterative (constraints) algorithms

  lstitr(1:natms)=.false. ! initialise lstitr
  Call constraints_tags(imcon,lstitr,lstopt,dxx,dyy,dzz,listot)

! application of constraint (quench) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  icyc=0
  safe=.false.

  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! initialise velocity correction arrays

     Do i=1,natms
        vxt(i)=0.0_wp
        vyt(i)=0.0_wp
        vzt(i)=0.0_wp
     End Do

! calculate velocity constraint corrections

     esig=0.0_wp
     Do k=1,ntcons
        i=lstopt(1,k)
        j=lstopt(2,k)

! for all constrained particles, native and shared

        If ( (i > 0 .and. j > 0) .and. (i <= natms .or. j <= natms) &
             .and. lfrzn(i)*lfrzn(j) == 0 ) Then

! if a pair is frozen and constraint bonded, it is more frozen
! than constrained (users!!!)

           amti=1.0_wp/weight(i)
           amtj=1.0_wp/weight(j)

! no corrections for frozen atoms

           If (lfrzn(i) /= 0) amti=0.0_wp
           If (lfrzn(j) /= 0) amtj=0.0_wp

! no corrections for frozen atoms

           If (lfrzn(i) /= 0) amti=0.0_wp
           If (lfrzn(j) /= 0) amtj=0.0_wp

           If (icyc == 1) Then
              dt(k) = Sqrt( dxx(k)**2 + dyy(k)**2 + dzz(k)**2 )
           End If

! calculate constraint force parameter - gamma

           gamma = (dxx(k)*(vxx(i)-vxx(j)) + dyy(k)*(vyy(i)-vyy(j)) + dzz(k)*(vzz(i)-vzz(j)))/dt(k)

           esig=Max(esig,0.5_wp*tolnce*Abs(gamma)*dt(k)/prmcon(listcon(0,k)))

           gamma = gamma / ( (amti+amtj) * dt(k) )

! improve approximate constraint velocity

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
     If (mxnode > 1) Call gcheck(safe)

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then

! update velocities

        Do k=1,ntcons
           i=lstopt(1,k)
           j=lstopt(2,k)

! for all constrained particles, native and shared

           If ( (i > 0 .and. j > 0) .and. (i <= natms .or. j <= natms) &
                .and. lfrzn(i)*lfrzn(j) == 0 ) Then

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

! transport velocity updates to other nodes

        If (lshmv_con) Call update_shared_units(natms,nlast,lsi,lsa,lishp_con,lashp_con,vxx,vyy,vzz)

     End If
  End Do

! error exit if quenching fails

  If (.not.safe) Call error(70)

  Deallocate (lstitr,         Stat=fail(1))
  Deallocate (lstopt,listot,  Stat=fail(2))
  Deallocate (dxx,dyy,dzz,dt, Stat=fail(3))
  Deallocate (vxt,vyt,vzt,    Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_quench deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine constraints_quench
