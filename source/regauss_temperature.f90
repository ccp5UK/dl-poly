Subroutine regauss_temperature(imcon,megrgd)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to regauss the instantenious system temperature
! by random pairwise swaps of the energy scaled momenta of dynamically
! active particles (no massless shells or massless RB members)
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode
  Use setup_module,       Only : nrite
  Use config_module,      Only : natms,nfree,lfrzn,lstfre, &
                                 weight,vxx,vyy,vzz
  Use rigid_bodies_module
  Use kinetic_module,     Only : getvom

  Implicit None

  Integer, Intent( In    ) :: imcon,megrgd

  Integer           :: fail,i,j,k,l,is,irgd,jrgd,lrgd,rgdtyp
  Real( Kind = wp ) :: uni,vom(1:3),tmp

  Integer, Allocatable :: ind(:),pair(:,:)

  fail=0
  Allocate (ind(1:natms),pair(1:2,natms/2), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'regauss_temperature allocation failure, node: ', idnode
     Call error(0)
  End If

! Create and index array containing the indices of the
! dynamically active particles and zeros for the inactive

  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        ind(i)=i
     Else
        ind(i)=0
     End If
  End Do

! Compress the index array

  i=1
  j=natms
  Do While (i < j)
     Do While (ind(j) == 0 .and. j > i)
        j=j-1
     End Do

     If (i < j) Then
        If (ind(i) == 0) Then
           ind(i)=ind(j)
           ind(j)=0
           j=j-1
        Else
           i=i+1
        End If
     End If
  End Do
  k=j

! Create a non-overlapping array of random pairs
! by exhausting the index array

  Do i=1,k/2
     Do l=1,2
        is=1+Int(Real(j,wp)*uni())
        pair(l,i)=ind(is)
        ind(is)=ind(j)
        ind(j)=0
        j=j-1
     End Do
  End Do

! Swap particles energies in the pair array

  Do i=1,k/2
     j=pair(1,i)
     l=pair(2,i)

     vom(1)=vxx(j)
     vom(2)=vyy(j)
     vom(3)=vzz(j)

     tmp=Sqrt(weight(l)/weight(j))
     vxx(j)=vxx(l)*tmp
     vyy(j)=vyy(l)*tmp
     vzz(j)=vzz(l)*tmp

     tmp=1.0_wp/tmp
     vxx(l)=vom(1)*tmp
     vyy(l)=vom(2)*tmp
     vzz(l)=vom(3)*tmp
  End Do

  If (megrgd > 0) Then

! quench RBs

     Call rigid_bodies_quench(imcon)

! remove centre of mass motion

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

! remove centre of mass motion

     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i) = vxx(i) - vom(1)
           vyy(i) = vyy(i) - vom(2)
           vzz(i) = vzz(i) - vom(3)
        End If
     End Do

  End If

  Deallocate (ind,pair, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'regauss_temperature deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine regauss_temperature
