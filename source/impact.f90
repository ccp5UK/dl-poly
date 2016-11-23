Subroutine impact(imd,emd,vmx,vmy,vmz,megrgd)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting impact on a particle
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,        Only : eu_ev
  Use comms_module,        Only : mxnode,gcheck
  Use config_module,       Only : natms,nlast,nfree,          &
                                  lfrzn,lfree,lstfre,lsi,lsa, &
                                  weight,vxx,vyy,vzz
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd, &
                                  rgdvxx,rgdvyy,rgdvzz
  Use core_shell_module,   Only : ntshl,listshl
  Use kinetic_module,      Only : getvom,l_vom,chvom

  Implicit None

  Integer,           Intent( In    ) :: imd,megrgd
  Real( Kind = wp ), Intent( In    ) :: emd,vmx,vmy,vmz

  Logical           :: safe = .true.

  Integer           :: local_index,i,j,irgd,jrgd,lrgd,rgdtyp
  Real( Kind = wp ) :: tmp,vom(1:3)


! Apply impact to the selected non-frozen, non-massless and non-shell particle

  i=local_index(imd,nlast,lsi,lsa)
  If (i > 0 .and. i <= natms) Then
     If (lfrzn(i) == 0 .and. lfree(i) == 0 .and. All(listshl(2,1:ntshl) /= imd)) Then
        tmp=Sqrt(2000.0_wp*emd*eu_ev/weight(i)/(vmx**2+vmy**2+vmz**2)) !emd is in keV=1000*eu_ev
        vxx(i)=tmp*vmx
        vyy(i)=tmp*vmy
        vzz(i)=tmp*vmz
     Else
        safe=.false.
     End If
  End If

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(610)

  Call chvom(.true.) ! Enable COM momentum removal

! remove centre of mass motion

  If (megrgd > 0) Then
     Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz)

     Do j=1,nfree
        i=lstfre(j)

        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do

     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

        If (rgdfrz(0,rgdtyp) == 0) Then
           rgdvxx(irgd) = rgdvxx(irgd) - vom(1)
           rgdvyy(irgd) = rgdvyy(irgd) - vom(2)
           rgdvzz(irgd) = rgdvzz(irgd) - vom(3)

           lrgd=listrgd(-1,irgd)
           Do jrgd=1,lrgd
              i=indrgd(jrgd,irgd) ! local index of particle/site

              If (i <= natms) Then
                 vxx(i) = vxx(i) - vom(1)
                 vyy(i) = vyy(i) - vom(2)
                 vzz(i) = vzz(i) - vom(3)
              End If
           End Do
        End If
     End Do
  Else
     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do
  End If

  Call chvom(l_vom) ! default to specification

End Subroutine impact
