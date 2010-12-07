Subroutine rigid_bodies_widths(imcon,rcut)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for checking RBs' widths compliency to < rcut
! (the system cutoff)
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode,mxnode,gmax
  Use setup_module
  Use config_module,       Only : cell,xxx,yyy,zzz
  Use rigid_bodies_module, Only : ntrgd,listrgd,indrgd

  Implicit None

  Integer,           Intent( In    ) :: imcon
  Real( Kind = wp ), Intent( In    ) :: rcut

  Integer           :: fail,irgd,jrgd,krgd,lrgd,mrgd,nrgd,rgdtyp
  Real( Kind = wp ) :: d,width

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*Max(mxrgd,mxtrgd)),gyy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
            gzz(1:mxlrgd*Max(mxrgd,mxtrgd)), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_widths allocation failure, node: ', idnode
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

! Get the COM vector and check of diameter safety of units

  krgd=0
  width=0.0_wp
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     lrgd=listrgd(-1,irgd)
     Do jrgd=1,lrgd
        krgd=krgd+1

        Do mrgd=jrgd+1,lrgd
           nrgd=krgd+mrgd-jrgd

           d=Sqrt((gxx(krgd)-gxx(nrgd))**2+(gyy(krgd)-gyy(nrgd))**2+(gzz(krgd)-gzz(nrgd))**2)
           width=Max(width,d)

           If (d > rcut) Write(nrite, Fmt='(1x,a,4i5,2f8.3)')  &
              'idnode, RB type, members(1,2), width > cutoff', &
              idnode,rgdtyp,jrgd,mrgd,width,rcut
        End Do
     End Do
  End Do

! Check if a RB unit has a diameter > the cutoff

  If (mxnode > 1) Call gmax(width)
  If (width > rcut) Then
     Call warning(8,width,rcut,0.0_wp)
     Call error(642)
  End If

  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_widths deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rigid_bodies_widths
