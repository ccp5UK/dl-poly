Subroutine bonds_compute()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bonds distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use site_module
  Use config_module
  Use bonds_module

  Implicit None

  Logical            :: zero
  Integer            :: fail,i,j,ig,kk,ll
  Real( Kind = wp )  :: kT2engo,delr,rdlr,factor,factor1, &
                        rrr,dvol,gofr,sum,gofr1,sum1,fed,dfed,tmp

  Real( Kind = wp ), Allocatable :: dstdbnd(:,:),pmf(:),vir(:)

  fail = 0
  Allocate (dstdbnd(0:mxgbnd,1:ldfbnd(0)), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bonds_compute - allocation failure, node: ', idnode
     Call error(0)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc) - kT2engu = boltz*temp

  kT2engo = 1.0_wp/engunit

! grid interval for pdf tables

  delr = rcbnd/Real(mxgbnd,wp)
  rdlr = 1.0_wp/delr

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        kk=kk+1
        ll=ll+typbnd(0,i)
     End If
  End Do

! normalisation factor

  factor = Real(kk,wp)/Real(ll*numbnd,wp)

  If (idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'BOND PDFs'
     Write(nrite,'(/,1x,a,2(i10,1x),f5.2,1x,i10)') 'types, grid, cutoff, frames: ',kk,mxgbnd,rcbnd,numbnd
  End If

! open RDF file and write headers

  If (idnode == 0) Then
     Open(Unit=npdfdt, File='BNDDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,2(i10,1x),f5.2,1x,i10)') '# types, grid, cutoff, frames: ',kk,mxgbnd,rcbnd,numbnd
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# r(Angs)  P_bond(r)  P_bond(r)/r^2   @   dr_bin = ',delr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(nrite,'(/,a,2(a8,1x),2(i10,1x))') 'Pb(r) : ',unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)
           Write(nrite,'(/,8x,a,6x,a,9x,a,/)') 'r','P(r)','Sum_P(r)'

           Write(npdfdt,'(/,a,2(a8,1x),2(i10,1x))') '# id, type, totals: ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)
        End If

! global sum of data on all nodes

        If (mxnode > 1) Call gsum(dstbnd(1:mxgbnd,i))

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxgbnd
           If (zero .and. ig < (mxgbnd-3)) zero=(dstbnd(ig+2,i) <= 0.0_wp)

! factor in degeneracy

           factor1=factor/Real(typbnd(0,i),wp)

           gofr= dstbnd(ig,i)*factor1
           sum = sum + gofr

           rrr = (Real(ig,wp)-0.5_wp)*delr
           dvol= 4.0_wp*pi*delr*(rrr*rrr+delr*delr/12.0_wp)
           gofr= gofr*rdlr

! null it if < 1.0e-6_wp

           If (gofr < 1.0e-6_wp) Then
              gofr1 = 0.0_wp
           Else
              gofr1 = gofr
           End If

           If (sum < 1.0e-6_wp) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

! print out information

           If (idnode == 0) Then
              If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") rrr,gofr1,sum1
              Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,gofr,gofr/dvol
           End If

           dstdbnd(ig,i) = gofr/dvol ! PDFs density
        End Do
     Else
        dstdbnd(:,i) = 0 ! PDFs density
     End If
  End Do

  If (idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=ntable, File='PMFBND', Status='replace')
     Write(ntable,'(a)') '# '//cfgname
     Write(ntable,'(a,f11.5,2i10,a,e15.7)') &
          '# ',delr,mxgbnd,kk,' conversion factor: kT -> energy units =',kT2engo

     Open(Unit=npdfdt, File='PDFBND', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,2(i10,1x),f5.2,1x,i10)') '# types, grid, cutoff, frames: ',kk,mxgbnd,rcbnd,numbnd
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# r(Angs)  P_bond(r)  P_bond(r)/r^2   @   dr_bin = ',delr
     Write(npdfdt,'(a)') '#'
  End If

!  Allocate (pmf(0:mxgbnd+2),vir(0:mxgbnd+2), Stat = fail)
!  If (fail > 0) Then
!     Write(nrite,'(/,1x,a,i0)') 'bonds_compute - allocation failure 2, node: ', idnode
!     Call error(0)
!  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Write(ntable,'(/,a,2(a8,1x),2(i10,1x))') &
           '# id, type, totals: ',unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)

! Smoothen and get derivatives

        Do ig=1,mxgbnd
           tmp = Real(ig,wp)-0.5_wp
           rrr = tmp*delr
           dvol= 4.0_wp*pi*delr*(rrr*rrr+delr*delr/12.0_wp)

           If (dstdbnd(ig,i) > zero_plus) Then
              fed = -Log(dstdbnd(ig,i))
              If (ig < mxgbnd-1) Then
                 If (dstdbnd(ig+1,i) <= zero_plus .and. dstdbnd(ig+1,i) > zero_plus) &
                    dstdbnd(ig+1,i) = 0.5_wp*(dstdbnd(ig,i)+dstdbnd(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If (ig == 1) Then
              If (fed < 0.0_wp .and. dstdbnd(ig+1,i) > zero_plus) Then
                 dfed = (Log(dstdbnd(ig+1,i))-fed)
              Else
                 dfed = 0.0_wp
              End If
           Else If (ig == mxgbnd) Then
              If (fed < 0.0_wp .and. dstdbnd(ig-1,i) > zero_plus) Then
                 dfed = (Log(dstdbnd(ig-1,i))-fed)
              Else
                 dfed = 0.0_wp
              End If
           Else If (dstdbnd(ig-1,i) > zero_plus) Then
              If (dstdbnd(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstdbnd(ig+1,i))-Log(dstdbnd(ig-1,i)))
              Else
                 dfed =-0.5_wp*Log(dstdbnd(ig-1,i))
              End If
           Else If (dstdbnd(ig+1,i) > zero_plus) Then
              dfed = 0.5_wp*Log(dstdbnd(ig+1,i))
           Else
              dfed = 0.0_wp
           End If

! print

           Write(ntable,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*tmp
           Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,dstdbnd(ig,i)*dvol,dstdbnd(ig,i)
!           pmf(ig) = fed  ! no engunit associated
!           vir(ig) = dfed ! no engunit & distance/bin associated
        End Do

! Retouch if needed to resample
! resampling can be done via 3pt interpolation
!
!        pmf(0)        = 2.0_wp*pmf(1)-pmf(2)
!        vir(0)        = 2.0_wp*vir(1)-vir(2)
!        pmf(mxgbnd+1) = 2.0_wp*pmf(mxgbnd)-pmf(mxgbnd-1)
!        vir(mxgbnd+1) = 2.0_wp*vir(mxgbnd)-vir(mxgbnd-1)
!        pmf(mxgbnd+2) = 2.0_wp*pmf(mxgbnd+1)-pmf(mxgbnd)
!        vir(mxgbnd+2) = 2.0_wp*vir(mxgbnd+1)-vir(mxgbnd)
     End If
  End Do

!  Deallocate (pmf,vir, Stat = fail)
!  If (fail > 0) Then
!     Write(nrite,'(/,1x,a,i0)') 'bonds_compute - deallocation failure 2, node: ', idnode
!     Call error(0)
!  End If

  If (idnode == 0) Then
     Close(Unit=ntable)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstdbnd, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bonds_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine bonds_compute
