Subroutine bonds_compute(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bonds distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : fourpi,nrite,npdfdt,npdgdt,mxgbnd,mxgbnd1,engunit,boltz,zero_plus
  Use site_module,   Only : unqatm
  Use config_module, Only : cfgname
  Use bonds_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: temp

  Logical           :: zero
  Integer           :: fail,i,j,ig,kk,ll,ngrid
  Real( Kind = wp ) :: kT2engo,delr,rdlr,dgrid,factor,factor1, &
                       rrr,dvol,pdfbnd,sum,pdfbnd1,sum1,fed0,fed,dfed,dfed0,tmp, &
                       fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstdbnd(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  fail = 0
  Allocate (dstdbnd(0:mxgbnd1,1:ldfbnd(0)), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bonds_compute - allocation failure, node: ', idnode
     Call error(0)
  End If

  ngrid = Max(1000,mxgbnd-4)

  fail = 0
  Allocate (pmf(0:ngrid), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,2(a,i0))') 'bonds_compute.f90 - allocation failure for pmf(0:',ngrid,'), node: ', idnode
     Call error(0)
  End If

  fail = 0
  Allocate (vir(0:ngrid), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,2(a,i0))') 'bonds_compute.f90 - allocation failure for vir(0:',ngrid,'), node: ', idnode
     Call error(0)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! grid interval for pdf tables

  delr = rcbnd/Real(mxgbnd1,wp)
  rdlr = 1.0_wp/delr

! resampling grid interval for pmf tables

  dgrid = rcbnd/Real(ngrid,wp)

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

  factor = 1.0_wp/Real(ncfbnd,wp)

  If (idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'BONDS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)'
     Write(nrite,'(/,1x,a,i10,1x,f8.3,3(1x,i10))') &
           'bins, cutoff, frames, types: ',mxgbnd1,rcbnd,ncfbnd,kk,ll
  End If

! open RDF file and write headers

  If (idnode == 0) Then
     Open(Unit=npdfdt, File='BNDDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a)') '# BONDS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dr_bin'
     Write(npdfdt,'(a,i10,1x,f8.3,2(1x,i10))') '# bins, cutoff, frames, types: ',mxgbnd1,rcbnd,ncfbnd,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# r(Angstroms)  PDF_norm(r)  PDF_norm(r)/dVol(r)   @   dr_bin = ',delr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(nrite,'(/,1x,a,2(a8,1x),2(i10,1x))') 'type, index, instances: ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)
           Write(nrite,'(/,1x,a,f8.5,/)') 'r(Angstroms)  P_bond(r)  Sum_P_bond(r)   @   dr_bin = ',delr

           Write(npdfdt,'(/,a,2(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)
        End If

! global sum of data on all nodes

        If (mxnode > 1) Call gsum(dstbnd(1:mxgbnd1,i))

! factor in instances (first, pdfbnd is normalised to unity)

        factor1=factor/Real(typbnd(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxgbnd1
           If (zero .and. ig < (mxgbnd1-3)) zero=(dstbnd(ig+2,i) <= 0.0_wp)

           pdfbnd= dstbnd(ig,i)*factor1
           sum = sum + pdfbnd

! null it if < 1.0e-5_wp

           If (pdfbnd < 1.0e-5_wp) Then
              pdfbnd1 = 0.0_wp
           Else
              pdfbnd1 = pdfbnd
           End If

           If (sum < 1.0e-5_wp) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           rrr = (Real(ig,wp)-0.5_wp)*delr
           dvol= fourpi*delr*(rrr*rrr+delr*delr/12.0_wp)

! now pdfbnd is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfbnd= pdfbnd*rdlr

! print out information

           If (idnode == 0) Then
              If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") rrr,pdfbnd1,sum1
              Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,pdfbnd,pdfbnd/rdlr/dvol
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

           dstdbnd(ig,i) = pdfbnd1/dvol ! PDFs density
        End Do
     Else
        dstdbnd(:,i) = 0 ! PDFs density
     End If
  End Do

  If (idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=npdgdt, File='BNDPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',delr*mxgbnd1,mxgbnd1,delr,kk, & 
          '   conversion factor(kT -> energy units) =',kT2engo

     Open(Unit=npdfdt, File='BNDTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',dgrid*ngrid,ngrid,dgrid,kk, & 
          '   conversion factor(kT -> energy units) =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        j=j+1

        If (idnode == 0)  Then
           Write(npdgdt,'(/,a,2(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i), &
                ' (type, index, instances)'
           Write(npdfdt,'(/,a,2(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i), &
                ' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxgbnd1
           tmp = Real(ig,wp)-0.5_wp
           rrr = tmp*delr

           If (dstdbnd(ig,i) > zero_plus) Then
              fed = -Log(dstdbnd(ig,i))-fed0
              If (fed0 <= zero_plus ) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < mxgbnd1-1) Then
                 If (dstdbnd(ig+1,i) <= zero_plus .and. dstdbnd(ig+2,i) > zero_plus) &
                    dstdbnd(ig+1,i) = 0.5_wp*(dstdbnd(ig,i)+dstdbnd(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstdbnd(ig,i) > zero_plus .and. dstdbnd(ig+1,i) > zero_plus) Then
                 dfed = Log(dstdbnd(ig+1,i)/dstdbnd(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == mxgbnd1) Then
              If      (dstdbnd(ig,i) > zero_plus .and. dstdbnd(ig-1,i) > zero_plus) Then
                 dfed = Log(dstdbnd(ig,i)/dstdbnd(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstdbnd(ig-1,i) > zero_plus) Then
              If (dstdbnd(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstdbnd(ig+1,i)/dstdbnd(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstdbnd(ig-1,i))
              End If
           Else If (dstdbnd(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstdbnd(ig+1,i))
           Else If (dfed > 0.0_wp) Then
              dfed = dfed0
           Else
              dfed =-dfed0
           End If

           pmf(ig) = fed
           vir(ig) = dfed

! Print
           If (idnode == 0) &
              Write(npdgdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*tmp
        End Do

        pmf(0)    = 2.0_wp*pmf(1)-pmf(2)
        vir(0)    = 2.0_wp*vir(1)-vir(2)
        pmf(ig)   = 2.0_wp*fed-pmf(ig-2)
        vir(ig)   = 2.0_wp*dfed-vir(ig-2)
        pmf(ig+1) = 2.0_wp*pmf(ig)-pmf(ig-1)
        vir(ig+1) = 2.0_wp*vir(ig)-vir(ig-1)
        
!       resample using 3pt interpolation
        
        Do ig=1,ngrid

           rrr = Real(ig,wp)*dgrid
           ll = Int(rrr/delr)

           if( ll > mxgbnd1 ) go to 113

! +0.5_wp due to half-a-bin shift in the original data
           coef = rrr/delr-Real(ll,wp)+0.5_wp

           fed0 = pmf(ll)
           fed1 = pmf(ll+1)
           fed2 = pmf(ll+2)

           t1 = fed0 + (fed1 - fed0)*coef
           t2 = fed1 + (fed2 - fed1)*(coef - 1.0_wp)

           fed = t1 + (t2-t1)*coef*0.5_wp

           dfed0 = vir(ll)
           dfed1 = vir(ll+1)
           dfed2 = vir(ll+2)

           t1 = dfed0 + (dfed1 - dfed0)*coef
           t2 = dfed1 + (dfed2 - dfed1)*(coef - 1.0_wp)

           dfed = t1 + (t2-t1)*coef*0.5_wp

           Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*rrr/delr
        End Do

113     if( ig < ngrid+1 ) write(nrite,*) &
             'bonds_compute():: PMF cut at r(max) = ',rrr,' ',ig,' ',ll,' ',i

     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstdbnd,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'bonds_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine bonds_compute
