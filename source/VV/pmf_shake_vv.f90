Subroutine pmf_shake_vv                &
           (imcon,mxshak,tolnce,tstep, &
           indpmf,pxx,pyy,pzz,         &
           xxx,yyy,zzz,strpmf,virpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying PMF constraint corrections after
! possible constrained motion
!
! Note: must be used in conjunction with integration algorithms
!       VV compliant
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsync,gcheck,gsum
  Use setup_module
  Use config_module, Only : natms,lfrzn
  Use pmf_module

  Implicit None

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tstep,tolnce
  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent( In    ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: strpmf(1:9),virpmf

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: rmass_pmf_unit(1:2),dis,dis2

  Logical                 :: safe
  Integer                 :: fail,ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: amt(1:2),gamma,gamm(1:2),tstep2,tmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: pxt,pyt,pzt,pt2,esig

  fail=0
  Allocate (pxt(1:mxpmf),pyt(1:mxpmf),pzt(1:mxpmf),pt2(1:mxpmf),esig(1:mxpmf), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_shake allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! Get reciprocal PMF units' masses

     Do jpmf=1,2
        If (pmffrz(jpmf) == mxtpmf(jpmf)) Then
           rmass_pmf_unit(jpmf)=0.0_wp
        Else
           rmass_pmf_unit(jpmf)=pmfwg1(0,jpmf)
        End If
     End Do

! set PMF constraint parameters

     dis=prmpmf
     dis2=dis**2
  End If

! squared timestep and reciprocal masses

  tstep2 = tstep*tstep
  amt = tstep2*rmass_pmf_unit

! Initialise constraint virial and stress

  virpmf=0.0_wp
  strpmf=0.0_wp

! application of PMF constraint (shake) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  safe=.false.
  icyc=0
  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! calculate temporary PMF units' COM vectors

     Call pmf_coms(imcon,indpmf,pxt,pyt,pzt)

! calculate maximum error in bondlength

     Do ipmf=1,ntpmf
        pt2(ipmf) =pxt(ipmf)**2+pyt(ipmf)**2+pzt(ipmf)**2 - dis2
        esig(ipmf)=0.5_wp*Abs(pt2(ipmf))/dis
     End Do

! global verification of convergence

     safe=(Maxval(esig(1:ntpmf)) < tolnce)
     If (mxnode > 1) Call gcheck(safe)

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then

! calculate PMF constraint forces

        Do ipmf=1,ntpmf

! calculate PMF constraint force parameter

           gamma = -pt2(ipmf) / &
                   ((amt(1)+amt(2))*(pxx(ipmf)*pxt(ipmf)+pyy(ipmf)*pyt(ipmf)+pzz(ipmf)*pzt(ipmf)))
           tmp   = gamma / Real(mxtpmf(1)+mxtpmf(2),wp)

           Do jpmf=1,2

! If this unit is present on my domain

              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
                 gamm(jpmf) = 0.5_wp*Real(1-2*Mod(jpmf,2),wp)*gamma*amt(jpmf)

                 Do k=1,mxtpmf(jpmf)
                    l=indpmf(k,jpmf,ipmf)

! for non-frozen domain particles accumulate PMF constraint stress
! and atomic position corrections

                    If (l <= natms .and. lfrzn(l) == 0) Then
                       strpmf(1) = strpmf(1) - tmp*pxx(ipmf)*pxx(ipmf)
                       strpmf(2) = strpmf(2) - tmp*pxx(ipmf)*pyy(ipmf)
                       strpmf(3) = strpmf(3) - tmp*pxx(ipmf)*pzz(ipmf)
                       strpmf(5) = strpmf(5) - tmp*pyy(ipmf)*pyy(ipmf)
                       strpmf(6) = strpmf(6) - tmp*pyy(ipmf)*pzz(ipmf)
                       strpmf(9) = strpmf(9) - tmp*pzz(ipmf)*pzz(ipmf)

                       xxx(l)=xxx(l)+pxx(ipmf)*gamm(jpmf)
                       yyy(l)=yyy(l)+pyy(ipmf)*gamm(jpmf)
                       zzz(l)=zzz(l)+pzz(ipmf)*gamm(jpmf)
                    End If
                 End Do
              End If

           End Do

        End Do

     End If
  End Do

! error exit for non-convergence

  If (.not.safe) Then
     Do k=0,mxnode-1
        If (idnode == k) Then
           Do ipmf=1,ntpmf
              If (esig(ipmf) >= tolnce .and. idnode == 0)                   &
                 Write(nrite,'(/,1x,3(a,i10),a,/,a,f8.2,a,1p,e12.4,a)')     &
                   '*** warning - global PMF constraint number', ipmf,      &
                   ' , with head particle numbers, U1:', listpmf(1,1,ipmf), &
                   ' & U2:', listpmf(1,2,ipmf), ' ,',                       &
                   ' , converges to a length of', Sqrt(pt2(ipmf)+dis2),     &
                   ' Angstoms with factor', esig(ipmf), ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(498)
  End If

! global sum of stress tensor

  If (mxnode > 1) Call gsum(strpmf)

! complete stress tensor (symmetrise)

  strpmf(4) = strpmf(2)
  strpmf(7) = strpmf(3)
  strpmf(8) = strpmf(6)

! total PMF constraint virial

  virpmf=-(strpmf(1)+strpmf(5)+strpmf(9))

  Deallocate (pxt,pyt,pzt,pt2,esig, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'pmf_shake deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine pmf_shake_vv
