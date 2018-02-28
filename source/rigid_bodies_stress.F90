Subroutine rigid_bodies_stress(strcom,ggx,ggy,ggz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate RB contributions to the atomic stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : mxnode,gsum
  Use setup_module,        Only : mxrgd,mxlrgd
  Use config_module,       Only : natms,lfrzn,fxx,fyy,fzz
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)
  Real( Kind = wp ), Intent( In    ) :: ggx(1:mxlrgd*mxrgd), &
                                        ggy(1:mxlrgd*mxrgd), &
                                        ggz(1:mxlrgd*mxrgd)

  Integer :: i,irgd,jrgd,krgd,lrgd,rgdtyp

! Initialise stress

  strcom=0.0_wp

! convert atomic virial to molecular
! note convention: virial(atom-atom) = -sum(Ri.Fi)
! : virial(com-com) = -sum(Rcom.Fcom) so
! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
              strcom(1)=strcom(1)-ggx(krgd)*fxx(i)
              strcom(2)=strcom(2)-ggx(krgd)*fyy(i)
              strcom(3)=strcom(3)-ggx(krgd)*fzz(i)
              strcom(4)=strcom(4)-ggy(krgd)*fxx(i)
              strcom(5)=strcom(5)-ggy(krgd)*fyy(i)
              strcom(6)=strcom(6)-ggy(krgd)*fzz(i)
              strcom(7)=strcom(7)-ggz(krgd)*fxx(i)
              strcom(8)=strcom(8)-ggz(krgd)*fyy(i)
              strcom(9)=strcom(9)-ggz(krgd)*fzz(i)
           End If
        End Do
     End If
  End Do

  If (mxnode > 1) Call gsum(strcom)

! Symmetrise

   strcom(2)=0.5_wp*(strcom(2)+strcom(4))
   strcom(4)=strcom(2)
   strcom(3)=0.5_wp*(strcom(3)+strcom(7))
   strcom(7)=strcom(3)
   strcom(6)=0.5_wp*(strcom(6)+strcom(8))
   strcom(8)=strcom(6)

End subroutine rigid_bodies_stress

Subroutine rigid_bodies_stre_s(strcom,ggx,ggy,ggz,fxx,fyy,fzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate RB contributions to the atomic stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : mxnode,gsum
  Use setup_module,        Only : mxatms,mxrgd,mxlrgd
  Use config_module,       Only : natms,lfrzn
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)
  Real( Kind = wp ), Intent( In    ) :: ggx(1:mxlrgd*mxrgd), &
                                        ggy(1:mxlrgd*mxrgd), &
                                        ggz(1:mxlrgd*mxrgd)
  Real( Kind = wp ), Intent( In    ) :: fxx(1:mxatms),fyy(1:mxatms),fzz(1:mxatms)

  Integer :: i,irgd,jrgd,krgd,lrgd,rgdtyp

! Initialise stress

  strcom=0.0_wp

! convert atomic virial to molecular
! note convention: virial(atom-atom) = -sum(Ri.Fi)
! : virial(com-com) = -sum(Rcom.Fcom) so
! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
              strcom(1)=strcom(1)-ggx(krgd)*fxx(i)
              strcom(2)=strcom(2)-ggx(krgd)*fyy(i)
              strcom(3)=strcom(3)-ggx(krgd)*fzz(i)
              strcom(4)=strcom(4)-ggy(krgd)*fxx(i)
              strcom(5)=strcom(5)-ggy(krgd)*fyy(i)
              strcom(6)=strcom(6)-ggy(krgd)*fzz(i)
              strcom(7)=strcom(7)-ggz(krgd)*fxx(i)
              strcom(8)=strcom(8)-ggz(krgd)*fyy(i)
              strcom(9)=strcom(9)-ggz(krgd)*fzz(i)
           End If
        End Do
     End If
  End Do

  If (mxnode > 1) Call gsum(strcom)

! Symmetrise

   strcom(2)=0.5_wp*(strcom(2)+strcom(4))
   strcom(4)=strcom(2)
   strcom(3)=0.5_wp*(strcom(3)+strcom(7))
   strcom(7)=strcom(3)
   strcom(6)=0.5_wp*(strcom(6)+strcom(8))
   strcom(8)=strcom(6)

End subroutine rigid_bodies_stre_s

Subroutine rigid_bodies_str_ss(strcom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate RB contributions to the atomic stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : idnode,mxnode,gsum
  Use setup_module,        Only : nrite,mxrgd,mxlrgd
  Use config_module,       Only : imcon,cell,natms,lfrzn,xxx,yyy,zzz,fxx,fyy,fzz
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd, &
                                  rgdxxx,rgdyyy,rgdzzz

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)

  Integer :: fail,i,irgd,jrgd,krgd,lrgd,rgdtyp

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_stress allocation failure, node: ', idnode
     Call error(0)
  End If

! Loop over all local RB units and get in the local scope of the unit

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

! COM distances

           gxx(krgd)=xxx(i)-rgdxxx(irgd)
           gyy(krgd)=yyy(i)-rgdyyy(irgd)
           gzz(krgd)=zzz(i)-rgdzzz(irgd)
        End Do
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Initialise stress

  strcom=0.0_wp

! convert atomic virial to molecular
! note convention: virial(atom-atom) = -sum(Ri.Fi)
! : virial(com-com) = -sum(Rcom.Fcom) so
! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
              strcom(1)=strcom(1)-gxx(krgd)*fxx(i)
              strcom(2)=strcom(2)-gxx(krgd)*fyy(i)
              strcom(3)=strcom(3)-gxx(krgd)*fzz(i)
              strcom(4)=strcom(4)-gyy(krgd)*fxx(i)
              strcom(5)=strcom(5)-gyy(krgd)*fyy(i)
              strcom(6)=strcom(6)-gyy(krgd)*fzz(i)
              strcom(7)=strcom(7)-gzz(krgd)*fxx(i)
              strcom(8)=strcom(8)-gzz(krgd)*fyy(i)
              strcom(9)=strcom(9)-gzz(krgd)*fzz(i)
           End If
        End Do
     End If
  End Do

  If (mxnode > 1) Call gsum(strcom)

! Symmetrise

   strcom(2)=0.5_wp*(strcom(2)+strcom(4))
   strcom(4)=strcom(2)
   strcom(3)=0.5_wp*(strcom(3)+strcom(7))
   strcom(7)=strcom(3)
   strcom(6)=0.5_wp*(strcom(6)+strcom(8))
   strcom(8)=strcom(6)

  fail = 0
  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_stress deallocation failure, node: ', idnode
     Call error(0)
  End If

End subroutine rigid_bodies_str_ss

Subroutine rigid_bodies_str__s(strcom,fxx,fyy,fzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate RB contributions to the atomic stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : idnode,mxnode,gsum
  Use setup_module,        Only : nrite,mxatms,mxrgd,mxlrgd
  Use config_module,       Only : imcon,cell,natms,lfrzn,xxx,yyy,zzz
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd, &
                                  rgdxxx,rgdyyy,rgdzzz

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)
  Real( Kind = wp ), Intent( In    ) :: fxx(1:mxatms),fyy(1:mxatms),fzz(1:mxatms)

  Integer :: fail,i,irgd,jrgd,krgd,lrgd,rgdtyp

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_stress allocation failure, node: ', idnode
     Call error(0)
  End If

! Loop over all local RB units and get in the local scope of the unit

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

! COM distances

           gxx(krgd)=xxx(i)-rgdxxx(irgd)
           gyy(krgd)=yyy(i)-rgdyyy(irgd)
           gzz(krgd)=zzz(i)-rgdzzz(irgd)
        End Do
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Initialise stress

  strcom=0.0_wp

! convert atomic virial to molecular
! note convention: virial(atom-atom) = -sum(Ri.Fi)
! : virial(com-com) = -sum(Rcom.Fcom) so
! virial(com-com) = virial(atom-atom)+sum((Ri-Rcom).Fi)

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i > 0 .and. i <= natms .and. lfrzn(i) == 0) Then
              strcom(1)=strcom(1)-gxx(krgd)*fxx(i)
              strcom(2)=strcom(2)-gxx(krgd)*fyy(i)
              strcom(3)=strcom(3)-gxx(krgd)*fzz(i)
              strcom(4)=strcom(4)-gyy(krgd)*fxx(i)
              strcom(5)=strcom(5)-gyy(krgd)*fyy(i)
              strcom(6)=strcom(6)-gyy(krgd)*fzz(i)
              strcom(7)=strcom(7)-gzz(krgd)*fxx(i)
              strcom(8)=strcom(8)-gzz(krgd)*fyy(i)
              strcom(9)=strcom(9)-gzz(krgd)*fzz(i)
           End If
        End Do
     End If
  End Do

  If (mxnode > 1) Call gsum(strcom)

! Symmetrise

   strcom(2)=0.5_wp*(strcom(2)+strcom(4))
   strcom(4)=strcom(2)
   strcom(3)=0.5_wp*(strcom(3)+strcom(7))
   strcom(7)=strcom(3)
   strcom(6)=0.5_wp*(strcom(6)+strcom(8))
   strcom(8)=strcom(6)

  fail = 0
  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_stress deallocation failure, node: ', idnode
     Call error(0)
  End If

End subroutine rigid_bodies_str__s
