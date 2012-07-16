Subroutine external_field_apply(imcon,keyshl,tstep,engfld,virfld)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for application of an external field
!
! Note: Only one field at a time is allowed
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module,  Only : pi,nrite,mxshl,mxatms
  Use config_module, Only : cell,natms,nfree,nlast,lsi,lsa,ltg, &
                            lfrzn,lstfre,weight,chge,           &
                            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use core_shell_module
  Use external_field_module

  Implicit None

  Integer,           Intent( In    ) :: imcon,keyshl
  Real( Kind = wp ), Intent( In    ) :: tstep ! for oscilating fields
  Real( Kind = wp ), Intent(   Out ) :: engfld,virfld

  Integer           :: i,j,ia,ib,fail(1:2),local_index, &
                       irgd,jrgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: gamma,rrr,rz,zdif,vxt,vyt,vzt,tmp, &
                       x(1:1),y(1:1),z(1:1)

  Integer,           Allocatable :: lstopt(:,:),list(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

! Recover megrgd

  megrgd=rgdmeg

! energy and virial accumulators

  engfld=0.0_wp
  virfld=0.0_wp

  If (keyfld == 1) Then

! electric field: prmfld(1-3) are field components

     Do i=1,natms
        If (lfrzn(i) == 0) Then
           fxx(i)=fxx(i) + chge(i)*prmfld(1)
           fyy(i)=fyy(i) + chge(i)*prmfld(2)
           fzz(i)=fzz(i) + chge(i)*prmfld(3)
        End If
     End Do

  Else If (keyfld == 2) Then

! oscillating shear: orthorhombic box:  Fx=a*Cos(b.2.pi.z/L)

     If (imcon /= 1 .and. imcon /= 2) Go To 10

     rz=2.0_wp*pi/cell(9)

     Do i=1,natms
        If (lfrzn(i) == 0) fxx(i)=fxx(i) + prmfld(1)*Cos(prmfld(2)*zzz(i)*rz)
     End Do

  Else If (keyfld == 3) Then

! continuous shear of walls : 2D periodic box (imcon=6)

     If (imcon /= 6) Go To 10

! shear rate=prmfld(1) angstrom per ps for non-frozen
! and non-weightless atoms at Abs(z) > prmfld(2)

     If (megrgd > 0) Then

! FPs
        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. Abs(zzz(i)) > prmfld(2)) &
           vxx(i)=0.5_wp*Sign(prmfld(1),zzz(i))
        End Do

! RBs

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) == 0) Then
              x=rgdxxx(irgd) ; y=rgdyyy(irgd) ; z=rgdzzz(irgd)
              Call images(imcon,cell,1,x,y,z)

              rz=z(1)
              If (Abs(rz) > prmfld(2)) Then
                 tmp=0.5_wp*Sign(prmfld(1),rz)
                 vxt=tmp-rgdvxx(irgd)

                 rgdvxx(irgd)=tmp
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd) ! local index of particle/site

                    If (i <= natms) vxx(i)=vxx(i)+vxt
                 End Do
              End If
           End If
        End Do

     Else

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. Abs(zzz(i)) > prmfld(2)) &
           vxx(i)=0.5_wp*Sign(prmfld(1),zzz(i))
        End Do

     End If

  Else If (keyfld == 4) Then

! gravitational field: field components given by prmfld(1-3)

     If (keyshl == 1) Then

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + prmfld(1)*weight(i)
              fyy(i)=fyy(i) + prmfld(2)*weight(i)
              fzz(i)=fzz(i) + prmfld(3)*weight(i)
           End If
        End Do

     Else

        fail=0
        Allocate (lstopt(1:2,1:mxshl),                       Stat=fail(1))
        Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms), Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply allocation failure, node: ', idnode
           Call error(0)
        End If

        Do i=1,ntshl
           lstopt(1,i)=local_index(listshl(1,i),nlast,lsi,lsa)
           lstopt(2,i)=local_index(listshl(2,i),nlast,lsi,lsa)
        End Do

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              oxt(i)=prmfld(1)*weight(i)
              oyt(i)=prmfld(2)*weight(i)
              ozt(i)=prmfld(3)*weight(i)
           Else ! for the sake of massless sites of RBs
              oxt(i)=0.0_wp
              oyt(i)=0.0_wp
              ozt(i)=0.0_wp
           End If
        End Do

        If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,oxt,oyt,ozt)

! Transfer cores' forces to shells

        Do i=1,ntshl
           ia=lstopt(1,i)
           ib=lstopt(2,i)
           If (ia > 0 .and. (ib > 0 .and. ib <= natms)) Then
              oxt(ib)=oxt(ia)
              oyt(ib)=oyt(ia)
              ozt(ib)=ozt(ia)
           End If
        End Do

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + oxt(i)
              fyy(i)=fyy(i) + oyt(i)
              fzz(i)=fzz(i) + ozt(i)
           End If
        End Do

        Deallocate (lstopt,      Stat=fail(1))
        Deallocate (oxt,oyt,ozt, Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply deallocation failure, node: ', idnode
           Call error(0)
        End If

     End If

  Else If (keyfld == 5) Then

! magnetic field: field components given by prmfld(1-3)

     If (keyshl == 1) Then

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + (vyy(i)*prmfld(3)-vzz(i)*prmfld(2))*chge(i)
              fyy(i)=fyy(i) + (vzz(i)*prmfld(1)-vxx(i)*prmfld(3))*chge(i)
              fzz(i)=fzz(i) + (vxx(i)*prmfld(2)-vyy(i)*prmfld(1))*chge(i)
           End If
        End Do

     Else

        fail=0
        Allocate (lstopt(1:2,1:mxshl),                       Stat=fail(1))
        Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms), Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply allocation failure, node: ', idnode
           Call error(0)
        End If

        Do i=1,ntshl
           lstopt(1,i)=local_index(listshl(1,i),nlast,lsi,lsa)
           lstopt(2,i)=local_index(listshl(2,i),nlast,lsi,lsa)
        End Do

! cores' velocities

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              oxt(i)=vxx(i)
              oyt(i)=vyy(i)
              ozt(i)=vzz(i)
           Else
              oxt(i)=0.0_wp
              oyt(i)=0.0_wp
              ozt(i)=0.0_wp
           End If
        End Do

        If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,oxt,oyt,ozt)

! Transfer cores' velocities to shells

        Do i=1,ntshl
           ia=lstopt(1,i)
           ib=lstopt(2,i)
           If (ia > 0 .and. (ib > 0 .and. ib <= natms)) Then
              oxt(ib)=oxt(ia)
              oyt(ib)=oyt(ia)
              ozt(ib)=ozt(ia)
           End If
        End Do

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + (oyt(i)*prmfld(3)-ozt(i)*prmfld(2))*chge(i)
              fyy(i)=fyy(i) + (ozt(i)*prmfld(1)-oxt(i)*prmfld(3))*chge(i)
              fzz(i)=fzz(i) + (oxt(i)*prmfld(2)-oyt(i)*prmfld(1))*chge(i)
           End If
        End Do

        Deallocate (lstopt,      Stat=fail(1))
        Deallocate (oxt,oyt,ozt, Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply deallocation failure, node: ', idnode
           Call error(0)
        End If

     End If

  Else If (keyfld == 6) Then

! containing sphere : r^(-n) potential

     Do i=1,natms
        If (lfrzn(i) == 0) Then
           rrr=Sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
           If (rrr > prmfld(4)) Then
              rrr=prmfld(2)-rrr
              If (rrr < 0.0_wp) rrr=0.1_wp

              gamma=prmfld(1)*rrr**(-prmfld(3))
              engfld=engfld + gamma
              gamma=-prmfld(3)*gamma/(rrr*rrr)

              fxx(i)=fxx(i) + gamma*xxx(i)
              fyy(i)=fyy(i) + gamma*yyy(i)
              fzz(i)=fzz(i) + gamma*zzz(i)
           End If
        End If
     End Do

     virfld=-9.0_wp*engfld

  Else If (keyfld == 7) Then

! repulsive wall (harmonic) starting at z0

     Do i=1,natms
        If (lfrzn(i) == 0 .and. prmfld(3)*zzz(i) > prmfld(3)*prmfld(2)) Then
           zdif=zzz(i)-prmfld(2)
           gamma=-prmfld(1)*zdif

           fzz(i)=fzz(i) + gamma
           engfld=engfld - gamma*zdif/2.0_wp
        End If
     End Do

!  Else If (keyfld == 8) Then
!
!! extension to oscilating electric field: prmfld(1-3) are field components
!! prmfld(4) is the oscillating frequency defined in ps^-1!
!! 'time' must put through Calling Sequence and md_*.f90 amended accordingly!
!
!     tmp=Sin(time*prmfld(4)*2.0_wp*pi)
!     Do i=1,natms
!        If (lfrzn(i) == 0) Then
!           fxx(i)=fxx(i) + chge(i)*prmfld(1)*tmp
!           fyy(i)=fyy(i) + chge(i)*prmfld(2)*tmp
!           fzz(i)=fzz(i) + chge(i)*prmfld(3)*tmp
!        End If
!     End Do
!
!   Else If (keyfld == 9) Then
!
!! extension to exzn- external field
!! prmfld(1) is the first atom of the water molecules to be restrained
!! prmfld(2) is the last atom of the water molecules to be restrained
!! prmfld(3) < prmfld(4) are the two limits, and prmfld(5) is the force
!! constant.  This will keep water away from the membrane region, and allow
!! the DMPC to slowly fill in the gaps before I let water molecules loose.
!!
!! P.-L. Chau  2007/07/25
!
!      ia = Nint(prmfld(1))
!      ib = Nint(prmfld(2))
!      Do i=1,natms
!         If ((ltg(i) >= ia .and. ltg(i) <= ib) .and. &
!             (zzz(i) > prmfld(3) .and. zzz(i) < prmfld(4)) .and. &
!             lfrzn(i) == 0) Then
!            tmp = prmfld(4) + prmfld(3)
!            If (zzz(i) <  tmp) zdif = zzz(i) - prmfld(3)
!            If (zzz(i) >= tmp) zdif = prmfld(4) - zzz(i)
!            gamma=-prmfld(5)*zdif
!            fzz(i)=fzz(i) + gamma
!            engfld=engfld - 0.5_wp*gamma*zdif
!         End If
!      End Do
!
!   Else If (keyfld == 10) Then
!
!! extension to exzn+ external field
!! prmfld(1) is the first atom of the water molecules to be restrained
!! prmfld(2) is the last atom of the water molecules to be restrained
!! prmfld(3) < prmfld(4) are the two limits, and prmfld(5) is the force
!! constant.  This will keep water within this region.
!!
!! P.-L. Chau  2007/10/10
!
!      ia = Nint(prmfld(1))
!      ib = Nint(prmfld(2))
!      Do i=1,natms
!         If ((ltg(i) >= ia .and. ltg(i) <= ib) .and. &
!             (zzz(i) < prmfld(3) .and. zzz(i) > prmfld(4)) .and. &
!             lfrzn(i) == 0) Then
!            If (zzz(i) < prmfld(3)) zdif = zzz(i) - prmfld(3)
!            If (zzz(i) > prmfld(4)) zdif = prmfld(4) - zzz(i)
!            gamma=-prmfld(5)*zdif
!            fzz(i)=fzz(i) + gamma
!            engfld=engfld - 0.5_wp*gamma*zdif
!         End If
!      End Do
!
!   Else If (keyfld == 11) Then
!
!! extension to zres external field: restrain molecule z-position
!! it assumes that the molecule diameter is less than the cutoff
!! prmfld(1) is number of first atom of restrained molecule
!! prmfld(2) is number of last atom of restrained molecule
!! prmfld(3) is the restraining constant and prmfld(4) < prmfld(5)
!! are the two limits in z-direction - z-min and z-max
!!
!! P.-L. Chau  2008/10/27
!
!      Allocate (list(1:mxatms),                            Stat=fail(1))
!      Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms), Stat=fail(2))
!      If (Any(fail > 0)) Then
!         Write(nrite,'(/,1x,a,i0)') 'external_field_apply allocation failure, node: ', idnode
!         Call error(0)
!      End If
!
!      ia = Nint(prmfld(1))
!      ib = Nint(prmfld(2))
!
!! Get in the local scope of the molecule
!
!      j=0         ! total number of atoms in molecule
!      tmp=0.0_wp  ! reciprocal of the molecule's unfrozen weight
!      Do i=1,nlast
!         If (ltg(i) >= ia .and. ltg(i) <= ib) Then
!            If (lfrzn(i) == 0) Then
!               j = j + 1
!
!               list(j) = i
!
!               tmp = tmp + weight(i)
!
!               oxt(j) = xxx(i) - xxx(list(1))
!               oyt(j) = yyy(i) - yyy(list(1))
!               ozt(j) = zzz(i) - zzz(list(1))
!            End If
!         End If
!      End Do
!      tmp=1.0_wp/tmp
!
!! minimum image convention for bond vectors
!
!      Call images(imcon,cell,j,oxt,oyt,ozt)
!
!! Get COM vector
!
!      vxt = 0.0_wp ; vyt = 0.0_wp ; vzt = 0.0_wp
!      Do i=1,j
!         vxt = vxt + weight(list(i)) * oxt(i)
!         vyt = vyt + weight(list(i)) * oyt(i)
!         vzt = vzt + weight(list(i)) * ozt(i)
!      End Do
!      vxt = vxt*tmp + xxx(list(1))
!      vyt = vyt*tmp + yyy(list(1))
!      vzt = vzt*tmp + zzz(list(1))
!
!      Do i=1,j
!         If (vzt < prmfld(4) .and. vzt > prmfld(5) .and. &
!             list(i) <= natms) Then
!            If (vzt < prmfld(4)) zdif = vzt - prmfld(4)
!            If (vzt > prmfld(5)) zdif = prmfld(5) - vzt
!            gamma=-prmfld(3)*zdif*weight(list(i))*tmp
!            fzz(list(i))=fzz(list(i)) + gamma
!            engfld=engfld - 0.5_wp*gamma*zdif
!         End If
!      End Do
!
!      Deallocate (list,        Stat=fail(1))
!      Deallocate (oxt,oyt,ozt, Stat=fail(2))
!      If (Any(fail > 0)) Then
!         Write(nrite,'(/,1x,a,i0)') 'external_field_apply deallocation failure, node: ', idnode
!         Call error(0)
!      End If

  Else

! unidentified field potential error exit

     Call error(454)

  End If

10 Continue

End Subroutine external_field_apply
