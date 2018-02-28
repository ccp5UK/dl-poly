Subroutine rigid_bodies_coms(xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing RBs coms
!
! Note: it assumes that all RBs' members are present and fresh
! (even those in the halo of the domain)
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,        Only : idnode
  Use setup_module,        Only : nrite,mxatms,mxrgd,mxlrgd
  Use config_module,       Only : imcon,cell
  Use rigid_bodies_module, Only : ntrgd,listrgd,indrgd,rgdwg1

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: rgdxxx(1:mxrgd),rgdyyy(1:mxrgd),rgdzzz(1:mxrgd)

  Integer           :: fail,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ) :: tmp

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*mxrgd),gyy(1:mxlrgd*mxrgd),gzz(1:mxlrgd*mxrgd), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_coms allocation failure, node: ', idnode
     Call error(0)
  End If

! Loop over all local RB units and get in the local scope of the unit

  krgd=0
  Do irgd=1,ntrgd
     lrgd=listrgd(-1,irgd)
     Do jrgd=1,lrgd
        krgd=krgd+1

        gxx(krgd) = xxx(indrgd(jrgd,irgd)) - xxx(indrgd(1,irgd))
        gyy(krgd) = yyy(indrgd(jrgd,irgd)) - yyy(indrgd(1,irgd))
        gzz(krgd) = zzz(indrgd(jrgd,irgd)) - zzz(indrgd(1,irgd))
     End Do
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Get the COM vector

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     rgdxxx(irgd)=0.0_wp
     rgdyyy(irgd)=0.0_wp
     rgdzzz(irgd)=0.0_wp

     lrgd=listrgd(-1,irgd)
     Do jrgd=1,lrgd
        krgd=krgd+1

        rgdxxx(irgd) = rgdxxx(irgd) + rgdwg1(jrgd,rgdtyp)*gxx(krgd)
        rgdyyy(irgd) = rgdyyy(irgd) + rgdwg1(jrgd,rgdtyp)*gyy(krgd)
        rgdzzz(irgd) = rgdzzz(irgd) + rgdwg1(jrgd,rgdtyp)*gzz(krgd)
     End Do
  End Do

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     tmp=1.0_wp/rgdwg1(0,rgdtyp)

     rgdxxx(irgd) = rgdxxx(irgd)*tmp + xxx(indrgd(1,irgd))
     rgdyyy(irgd) = rgdyyy(irgd)*tmp + yyy(indrgd(1,irgd))
     rgdzzz(irgd) = rgdzzz(irgd)*tmp + zzz(indrgd(1,irgd))
  End Do

  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_coms deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rigid_bodies_coms
