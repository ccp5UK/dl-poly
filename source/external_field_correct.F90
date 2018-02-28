Subroutine external_field_correct(engfld)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for correcting an external field application
!
! Note: Only one field at a time is allowed
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
! amnded    - i.t.todorov september 2015 : gsum engfld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : mxnode,gsum
  Use config_module, Only : imcon,cell,natms,nfree,ltg, &
                            lfrzn,lstfre,weight,  &
                            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use external_field_module

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: engfld

  Integer           :: i,j,ia,ib, irgd,jrgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: rz,vxt,tmp,rtmp(1:2), &
                       x(1:1),y(1:1),z(1:1)

! Recover megrgd

  megrgd=rgdmeg

  If (keyfld == 3) Then

! continuous shear of walls : 2D periodic box (imcon=6)

     If (imcon /= 6) Return

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

  Else If (keyfld == 8) Then

     engfld = 0.0_wp

! xpist - piston wall pushing down along the X=bxc direction
! prmfld(1) is the first atom of the layer of molecules (membrane) to be pushed
! prmfld(2) is the last atom of the layer of molecules (membrane) to be pushed
! prmfld(3) is the pressure applied to the layer of molecules (membrane) in the
! +X=bxc direction - i.e. left to right.  The layer plane is defined as _|_ bxc

     If (imcon /= 1 .and. imcon /= 2) Return

     ia = Nint(prmfld(1))
     ib = Nint(prmfld(2))

     rtmp=0.0_wp  ! average velocity and force per atom in x direction of the piston
     Do i=1,natms ! preserve momentum and velocity in the direction of the push
        If (ltg(i) >= ia .and. ltg(i) <= ib) Then
           rtmp(1)=rtmp(1)+weight(i)*vxx(i) ; rtmp(2)=rtmp(2)+fxx(i)
           vyy(i) = 0.0_wp                  ; fyy(i) = 0.0_wp
           vzz(i) = 0.0_wp                  ; fzz(i) = 0.0_wp
        End If
     End Do
     If (mxnode > 1) Call gsum(rtmp) ! net velocity and force to ensure solid wall behaviour

     rtmp(1)=rtmp(1)/mass             ! averaged velocity per particle
     rtmp(2)=(rtmp(2)+prmfld(3))/mass ! averaged acceleration of the slab

     Do i=1,natms
        If (ltg(i) >= ia .and. ltg(i) <= ib) Then
           engfld=weight(i)*(rtmp(1)-vxx(i))**2 ! must change E_kin to reflect solidity
           vxx(i)=rtmp(1)
           fxx(i)=rtmp(2)*weight(i) ! force per particle
        End If
     End Do

     engfld=0.5_wp*engfld

     If (mxnode > 1) Call gsum(engfld)

  End If

End Subroutine external_field_correct
