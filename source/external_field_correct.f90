Subroutine external_field_correct(imcon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for correcting an external field application
!
! Note: Only one field at a time is allowed
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use config_module, Only : cell,natms,nfree,    &
                            lfrzn,lstfre,weight, &
                            zzz,vxx
  Use rigid_bodies_module
  Use external_field_module

  Implicit None

  Integer,           Intent( In    ) :: imcon

  Integer           :: i,j, irgd,jrgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: rz,vxt,tmp, &
                       x(1:1),y(1:1),z(1:1)

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

  End If

10 Continue

End Subroutine external_field_correct
