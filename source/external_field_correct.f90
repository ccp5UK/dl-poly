Subroutine external_field_correct(imcon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for correcting an external field application
!
! Note: Only one field at a time is allowed
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use config_module, Only : cell,natms,nfree,ltg, &
                            lfrzn,lstfre,weight,  &
                            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use external_field_module

  Implicit None

  Integer,           Intent( In    ) :: imcon

  Integer           :: i,j,ia,ib, irgd,jrgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: rz,vxt,tmp, &
                       x(1:1),y(1:1),z(1:1),celprp(10)

! Recover megrgd

  megrgd=rgdmeg

  If (keyfld == 3) Then

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

  Else If (keyfld == 8) Then

! zpist - piston wall pushing down along the Z=axb direction
! prmfld(1) is the first atom of the layer of molecules (membrane) to be pushed
! prmfld(2) is the last atom of the layer of molecules (membrane) to be pushed
! prmfld(3) is the pressure applied to the layer of molecules (membrane) in the
! -Z=-axb direction - i.e. top to bottom.  The layer plane is defined as _|_ axb

     ia = Nint(prmfld(1))
     ib = Nint(prmfld(2))

     tmp=0.0_wp ! average force per atom of the piston
     Do i=1,natms
        If ((ltg(i) >= ia .and. ltg(i) <= ib)) Then
           vxx(i) = 0.0_wp ; fxx(i) = 0.0_wp
           vyy(i) = 0.0_wp ; fyy(i) = 0.0_wp
           vzz(i) = 0.0_wp ; tmp=tmp+fzz(i)
        End If
     End Do
     If (mxnode > 1) Call gsum(tmp)
     tmp=tmp/Real(ib-ia+1) ! solid wall behaviour is ensured

     Call dcell(cell,celprp)
     tmp=tmp-prmfld(3)*(celprp(10)/celprp(9))

     Do i=1,natms
        If ((ltg(i) >= ia .and. ltg(i) <= ib)) fzz(i)=tmp
     End Do

  End If

10 Continue

End Subroutine external_field_correct
