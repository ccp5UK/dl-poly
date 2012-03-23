Subroutine constraints_shake_lfv        &
           (imcon,mxshak,tolnce,tstep, &
           lstopt,dxx,dyy,dzz,listot,  &
           xxx,yyy,zzz,strcon,vircon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying bond constraint corrections after
! unconstrained motion
!
! Note: must be used in conjunction with integration algorithms
!       LFV compliant
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Centre for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY 4 is. You should have received these sources from the
! STFC Daresbury Laboratory.


  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsync,gcheck,gsum
  Use setup_module
  Use config_module, Only : cell,natms,nlast,lsi,lsa,lfrzn,weight
  Use constraints_module
! cuda: needed to pre-calculate what positions we will be sending.
  Use domains_module, Only : mop

#ifdef COMPILE_CUDA
  Use dl_poly_cuda_module
  Use iso_c_binding
#endif

  Implicit None

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce,tstep
  Integer,           Intent( In    ) :: lstopt(0:2,1:mxcons)
  Real( Kind = wp ), Intent( In    ) :: dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons)
  Integer,           Intent( In    ) :: listot(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: strcon(1:9)
  Real( Kind = wp ), Intent(   Out ) :: vircon

  Logical           :: safe
  Integer           :: fail(1:3),i,j,k,icyc

  Real( Kind = wp ) :: amti,amtj,dli,dlj,gamma,gammi,gammj,tstep2
  Real( Kind = wp ) :: esigmax

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt
  Real( Kind = wp ), Dimension( : ), Allocatable :: dxt,dyt,dzt,dt2,esig

#ifdef COMPILE_CUDA
! CUDA: Using the same cuda sources for acceleration of both lfv and vv routines
! define a switch here to tell the port what version we're in (lfv or vv)
  Logical( c_bool ) :: is_vv
  is_vv = .false.
#endif

  fail=0
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),                              Stat=fail(1))
  Allocate (dxt(1:mxcons),dyt(1:mxcons),dzt(1:mxcons),dt2(1:mxcons),esig(1:mxcons), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_shake allocation failure, node: ', idnode
     Call error(0)
  End If

#ifdef COMPILE_CUDA
  Call start_timing_constraints_shake()
#endif

! Initialise constraint virial and stress

  vircon=0.0_wp
  strcon=0.0_wp

! squared timestep

  tstep2 = tstep*tstep

! application of constraint (shake) algorithm
! start bond vectors are dxx = xxx(i) - xxx(j) etc.
! Initialise number of cycles to zero and unsafe passage of the algorithm

  icyc=0
  safe=.false.


#ifdef COMPILE_CUDA
  If (dl_poly_cuda_offload_constraints_shake() .and. &
       dl_poly_cuda_is_cuda_capable()          .and. &
       dl_poly_cuda_constraints_shake_ntcons_enough_to_offload(ntcons)) Then

     Call constraints_shake_cuda_initialise(&
          ntcons, mxcons, mxatms, natms,imcon,&
          lsi,lsa,lishp_con,lashp_con,mop,mxbuff,nlast,&
          lstopt, lfrzn, listcon, listot,&
          prmcon, weight,&
          dxx, dyy, dzz, dxt, dyt, dzt, dt2, tstep2, tolnce, cell,&
          xxx, yyy, zzz, xxt, yyt, zzt, strcon, is_vv)

     Do While ((.not.safe) .and. icyc < mxshak)
        icyc=icyc+1

        If (lshmv_con) Then
           Call update_shared_units_constraints_shake_specific(&
                natms,nlast,lsi,lsa,lishp_con,lashp_con,xxx,yyy,zzz)
           Call constraints_shake_cuda_gather_hs_scatter_dv()
        End If
        Call constraints_shake_cuda_invoke_th(esigmax)

        safe=(esigmax < tolnce)
        If (mxnode > 1) Call gcheck(safe)

        If (.not.safe) Then
           Call constraints_shake_cuda_invoke_bh()
           Call constraints_shake_cuda_invoke_correct_positions()
           Call constraints_shake_cuda_gather_dv_scatter_hs()
        End If
     End Do

     Call constraints_shake_cuda_invoke_strcon_extraction()
     Call constraints_shake_cuda_copy_local_positions_to_host()
     Call constraints_shake_cuda_finalise()


! error exit for non-convergence

     If (.not.safe) Call error(105)

  Else
#endif

   Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! update positions globally: transport position updates of shared atoms to other nodes

     If (lshmv_con) Then
        Call update_shared_units(&
             natms,nlast,lsi,lsa,lishp_con,lashp_con,xxx,yyy,zzz)
     End If

! calculate temporary bond vector

     Do k=1,ntcons
        If (lstopt(0,k) == 0) Then
           i=lstopt(1,k)
           j=lstopt(2,k)

           dxt(k)=xxx(i)-xxx(j)
           dyt(k)=yyy(i)-yyy(j)
           dzt(k)=zzz(i)-zzz(j)
        Else ! DEBUG
!           dxt(k)=0.0_wp
!           dyt(k)=0.0_wp
!           dzt(k)=0.0_wp
        End If
     End Do

! periodic boundary condition

     Call images(imcon,cell,ntcons,dxt,dyt,dzt)

! calculate maximum error in bondlength

     Do k=1,ntcons
        If (lstopt(0,k) == 0) Then
           dt2(k) =dxt(k)**2+dyt(k)**2+dzt(k)**2 - prmcon(listcon(0,k))**2
           esig(k)=0.5_wp*Abs(dt2(k))/prmcon(listcon(0,k))
        Else
           dt2(k) =0.0_wp
           esig(k)=0.0_wp
        End If
     End Do

! global verification of convergence

     safe=(Maxval(esig(1:ntcons)) < tolnce)
     If (mxnode > 1) Call gcheck(safe)

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then

! initialise position correction arrays

        Do i=1,natms
           xxt(i)=0.0_wp
           yyt(i)=0.0_wp
           zzt(i)=0.0_wp
        End Do

! calculate constraint forces

        Do k=1,ntcons
           If (lstopt(0,k) == 0) Then
               i=lstopt(1,k)
               j=lstopt(2,k)

               amti=tstep2/weight(i)
               amtj=tstep2/weight(j)

! no corrections for frozen atoms

               If (lfrzn(i) /= 0) amti=0.0_wp
               If (lfrzn(j) /= 0) amtj=0.0_wp

! calculate constraint force parameter

              gamma = 0.5_wp * dt2(k) / ((amti+amtj)*(dxx(k)*dxt(k)+dyy(k)*dyt(k)+dzz(k)*dzt(k)))

              If (i <= natms) Then

! accumulate bond stress

                 strcon(1) = strcon(1) - gamma*dxx(k)*dxx(k)
                 strcon(2) = strcon(2) - gamma*dxx(k)*dyy(k)
                 strcon(3) = strcon(3) - gamma*dxx(k)*dzz(k)
                 strcon(5) = strcon(5) - gamma*dyy(k)*dyy(k)
                 strcon(6) = strcon(6) - gamma*dyy(k)*dzz(k)
                 strcon(9) = strcon(9) - gamma*dzz(k)*dzz(k)

! calculate atomic position constraint corrections

                 If (lfrzn(i) == 0) Then
                    gammi =-gamma*amti
                    xxt(i)=xxt(i)+dxx(k)*gammi
                    yyt(i)=yyt(i)+dyy(k)*gammi
                    zzt(i)=zzt(i)+dzz(k)*gammi
                 End If

              End If

              If (j <= natms .and. lfrzn(j) == 0) Then
                 gammj = gamma*amtj
                 xxt(j)=xxt(j)+dxx(k)*gammj
                 yyt(j)=yyt(j)+dyy(k)*gammj
                 zzt(j)=zzt(j)+dzz(k)*gammj
              End If
           End If
        End Do

! update positions locally

        Do k=1,ntcons
           If (lstopt(0,k) == 0) Then
              i=lstopt(1,k)
              j=lstopt(2,k)

! apply position corrections if non-frozen

              If (i <= natms .and. lfrzn(i) == 0) Then
                 dli = 1.0_wp/Real(listot(i),wp)
                 xxx(i)=xxx(i)+xxt(i)*dli
                 yyy(i)=yyy(i)+yyt(i)*dli
                 zzz(i)=zzz(i)+zzt(i)*dli
              End If

              If (j <= natms .and. lfrzn(j) == 0) Then
                 dlj = 1.0_wp/Real(listot(j),wp)
                 xxx(j)=xxx(j)+xxt(j)*dlj
                 yyy(j)=yyy(j)+yyt(j)*dlj
                 zzz(j)=zzz(j)+zzt(j)*dlj
              End If
           End If
        End Do
     End If

  End Do

! error exit for non-convergence

  If (.not.safe) Then
     Do i=0,mxnode-1
        If (idnode == i) Then
           Do k=1,ntcons
              If (esig(k) >= tolnce)                                          &
                 Write(nrite,'(/,1x,3(a,i10),a,/,a,f8.2,a,1p,e12.4,a)')       &
                      '*** warning - global constraint number', listcon(0,k), &
                      ' , with particle numbers:', listcon(1,k),              &
                      ' &', listcon(2,k), ' ,', ' converges to a length of ', &
                      Sqrt(dt2(k)+prmcon(listcon(0,k))**2),                   &
                      ' Angstroms with a factor', esig(k),                    &
                      ' ,contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(105)
  End If

#ifdef COMPILE_CUDA
  End If
#endif

! global sum of stress tensor

  If (mxnode > 1) Call gsum(strcon)

! complete stress tensor (symmetrise)

  strcon(4) = strcon(2)
  strcon(7) = strcon(3)
  strcon(8) = strcon(6)

! total constraint virial

  vircon=-(strcon(1)+strcon(5)+strcon(9))

  Deallocate (xxt,yyt,zzt,          Stat=fail(1))
  Deallocate (dxt,dyt,dzt,dt2,esig, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'constraints_shake deallocation failure, node:', idnode
     Call error(0)
  End If

#ifdef COMPILE_CUDA
  Call stop_timing_constraints_shake()
#endif

End Subroutine constraints_shake_lfv
