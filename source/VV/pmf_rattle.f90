Subroutine pmf_rattle                 &
           (mxshak,tolnce,tstep,lcol, &
           indpmf,pxx,pyy,pzz,        &
           vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying PMF constraint corrections after
! possible constrained motion
!
! Note: must be used in conjunction with integration algorithms
!       VV compliant
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module
  Use config_module, Only : natms,lfrzn
  Use pmf_module

  Implicit None

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce,tstep
  Logical,           Intent( In    ) :: lcol
  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)

  Logical,           Save :: newjob = .true.

  Logical                 :: safe
  Integer                 :: fail(1:2),ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: amt(1:2),dis,esig,gamma,gamm(1:2)

  Real( Kind = wp ), Dimension( : )   , Allocatable :: vxt,vyt,vzt
  Real( Kind = wp ), Dimension( :, : ), Allocatable :: xpmf,ypmf,zpmf

  fail=0
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),             Stat=fail(1))
  Allocate (xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_rattle allocation failure, node: ', idnode
     Call error(0)
  End If

! Get PMF units' reciprocal masses

  If (newjob) Then
     newjob = .false.

     Do jpmf=1,2
        If (pmffrz(jpmf) == mxtpmf(jpmf)) Then
           amt(jpmf)=0.0_wp
        Else
           amt(jpmf)=pmfwg1(0,jpmf)
        End If
     End Do
  End If

! normalise PMF constraint vectors

  Do ipmf=1,ntpmf
     dis=1.0_wp/Sqrt(pxx(ipmf)**2+pyy(ipmf)**2+pzz(ipmf)**2)
     pxx(ipmf)=pxx(ipmf)*dis
     pyy(ipmf)=pyy(ipmf)*dis
     pzz(ipmf)=pzz(ipmf)*dis
  End Do

! application of PMF constraint (rattle) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  safe=.false.
  icyc=0
  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! initialise velocity correction arrays

     Do l=1,natms
        vxt(l)=0.0_wp
        vyt(l)=0.0_wp
        vzt(l)=0.0_wp
     End Do

! calculate temporary COM velocity of each unit

     Call pmf_vcoms(indpmf,xpmf,ypmf,zpmf)

! calculate PMF velocity corrections

     esig=0.0_wp
     Do ipmf=1,ntpmf

! calculate constraint force parameter - gamma

        gamma = pxx(ipmf)*(xpmf(1,ipmf)-xpmf(2,ipmf)) + &
                pyy(ipmf)*(ypmf(1,ipmf)-ypmf(2,ipmf)) + &
                pzz(ipmf)*(zpmf(1,ipmf)-zpmf(2,ipmf))

        esig=Max(esig,0.5_wp*tstep*Abs(gamma))

        gamma = gamma / (amt(1)+amt(2))

        Do jpmf=1,2

! If this unit is present on my domain

           If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
              gamm(jpmf) = Real(1-2*Mod(jpmf,2),wp)*gamma*amt(jpmf)
              Do k=1,mxtpmf(jpmf)
                 l=indpmf(k,jpmf,ipmf)

! improve approximate PMF particles velocity and force (if non-frozen)

                 If (l > 0 .and. l <= natms) Then ! l is a domain particle
                    If (lfrzn(l) == 0) Then
                       vxt(l)=vxt(l)+pxx(ipmf)*gamm(jpmf)
                       vyt(l)=vyt(l)+pyy(ipmf)*gamm(jpmf)
                       vzt(l)=vzt(l)+pzz(ipmf)*gamm(jpmf)
                    End If
                 End If
              End Do
           End If

        End Do

     End Do

! global verification of convergence

     safe=(esig < tolnce)
     If (mxnode > 1) Call gcheck(safe,"enforce")

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then
        Do ipmf=1,ntpmf
           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
                 Do k=1,mxtpmf(jpmf)
                    l=indpmf(k,jpmf,ipmf)
                    If (l > 0 .and. l <= natms) Then ! l is a domain particle
                       If (lfrzn(l) == 0) Then
                          vxx(l)=vxx(l)+vxt(l)
                          vyy(l)=vyy(l)+vyt(l)
                          vzz(l)=vzz(l)+vzt(l)
                       End If
                    End If
                 End Do
              End If
           End Do
        End Do
     End If
  End Do

  If (.not.safe) Then ! error exit for non-convergence
     Call error(499)
  Else ! Collect per call and per step passage statistics
     passpmf(1,1,2)=icyc-1
     passpmf(3,1,2)=passpmf(2,1,2)*passpmf(3,1,2)
     passpmf(2,1,2)=passpmf(2,1,2)+1
     passpmf(3,1,2)=passpmf(3,1,2)/passpmf(2,1,2)+passpmf(1,1,2)/passpmf(2,1,2)
     passpmf(4,1,2)=Min(passpmf(1,1,2),passpmf(4,1,2))
     passpmf(5,1,2)=Max(passpmf(1,1,2),passpmf(5,1,2))

     passpmf(1,2,2)=passpmf(1,2,2)+passpmf(1,1,2)
     If (lcol) Then ! Collect
        passpmf(3,2,2)=passpmf(2,2,2)*passpmf(3,2,2)
        passpmf(2,2,2)=passpmf(2,2,2)+1
        passpmf(3,2,2)=passpmf(3,2,2)/passpmf(2,2,2)+passpmf(1,2,2)/passpmf(2,2,2)
        passpmf(4,2,2)=Min(passpmf(1,2,2),passpmf(4,2,2))
        passpmf(5,2,2)=Max(passpmf(1,2,2),passpmf(5,2,2))
        passpmf(1,2,2)=0.0_wp ! Reset
     End If
     passpmf(1,1,2)=0.0_wp ! Reset
  End If

  Deallocate (vxt,vyt,vzt,    Stat=fail(1))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_rattle deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pmf_rattle
