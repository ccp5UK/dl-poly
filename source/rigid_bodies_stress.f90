Subroutine rigid_bodies_stress(strcom,ggx,ggy,ggz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate RB contributions to the atomic stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : mxnode,gsum
  Use setup_module,        Only : mxtrgd,mxrgd,mxlrgd
  Use config_module,       Only : natms,fxx,fyy,fzz
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)
  Real( Kind = wp ), Intent( In    ) :: ggx(1:mxlrgd*Max(mxrgd,mxtrgd)), &
                                        ggy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
                                        ggz(1:mxlrgd*Max(mxrgd,mxtrgd))

  Integer :: i,irgd,jrgd,krgd,lrgd,rgdtyp

! Initilise stress

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
           If (indrgd(jrgd,irgd) > 0 .and. indrgd(jrgd,irgd) <= natms .and. &
               rgdfrz(jrgd,rgdtyp) == 0) Then
              i=indrgd(jrgd,irgd) ! local index of particle/site

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
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : mxnode,gsum
  Use setup_module,        Only : mxatms,mxtrgd,mxrgd,mxlrgd
  Use config_module,       Only : natms
  Use rigid_bodies_module, Only : ntrgd,rgdfrz,listrgd,indrgd

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)
  Real( Kind = wp ), Intent( In    ) :: ggx(1:mxlrgd*Max(mxrgd,mxtrgd)), &
                                        ggy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
                                        ggz(1:mxlrgd*Max(mxrgd,mxtrgd))
  Real( Kind = wp ), Intent( In    ) :: fxx(1:mxatms),fyy(1:mxatms),fzz(1:mxatms)

  Integer :: i,irgd,jrgd,krgd,lrgd,rgdtyp

! Initilise stress

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
           If (indrgd(jrgd,irgd) > 0 .and. indrgd(jrgd,irgd) <= natms .and. &
               rgdfrz(jrgd,rgdtyp) == 0) Then
              i=indrgd(jrgd,irgd) ! local index of particle/site

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
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode,mxnode,gsum
  Use setup_module,        Only : nrite,mxtrgd,mxrgd,mxlrgd
  Use config_module,       Only : cell,natms,xxx,yyy,zzz,fxx,fyy,fzz
  Use rigid_bodies_module, Only : ntrgd,rgdimc,rgdfrz,listrgd,indrgd

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: strcom(1:9)

  Integer :: fail,i,irgd,jrgd,krgd,lrgd,rgdtyp,imcon

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*Max(mxrgd,mxtrgd)),gyy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
            gzz(1:mxlrgd*Max(mxrgd,mxtrgd)), Stat = fail)
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

           gxx(krgd) = xxx(indrgd(jrgd,irgd)) - xxx(indrgd(1,irgd))
           gyy(krgd) = yyy(indrgd(jrgd,irgd)) - yyy(indrgd(1,irgd))
           gzz(krgd) = zzz(indrgd(jrgd,irgd)) - zzz(indrgd(1,irgd))
        End Do
     End If
  End Do

! recover imcon

  imcon=rgdimc

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Initilise stress

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
           If (indrgd(jrgd,irgd) > 0 .and. indrgd(jrgd,irgd) <= natms .and. &
               rgdfrz(jrgd,rgdtyp) == 0) Then
              i=indrgd(jrgd,irgd) ! local index of particle/site

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
