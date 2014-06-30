Subroutine rdf_compute(rcut,temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating radial distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - t.forester march 1994
! amended   - i.t.todorov june 2014
! amended   - a.v.brukhno june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : fourpi,boltz,nrite,nrdfdt,npdfdt,npdgdt, &
                            mxgrdf,engunit,zero_plus
  Use site_module,   Only : ntpatm,unqatm,dens
  Use config_module, Only : cfgname,volm
  Use rdf_module


  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rcut,temp

  Logical           :: zero
  Integer           :: fail,i,ig,ia,ib,kk,ll,ngrid
  Real( Kind = wp ) :: kT2engo,delr,dvol,dgrid,factor,                &
                       gofr,gofr1,rrr,sum,sum1,                       &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2, &
                       coef,t1,t2,pdfzero,zeroln

  Real( Kind = wp ), Allocatable :: dstdrdf(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  fail = 0
  Allocate (dstdrdf(0:mxgrdf,1:ntprdf),pmf(0:mxgrdf+2),vir(0:mxgrdf+2), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rdf_compute - allocation failure, node: ', idnode
     Call error(0)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! open RDF file and Write headers

  If (idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf
  End If

! grid interval for rdf tables

  delr = rcut/Real(mxgrdf,wp)

! resampling grid and grid interval for pmf tables

  ngrid = Max(1000,mxgrdf-4)
  dgrid = rcut/Real(ngrid,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-6_wp

  If (idnode == 0) Write(nrite,"(/,/,12X,'RADIAL DISTRIBUTION FUNCTIONS',/,/, &
     & ' calculated using ',i10,' configurations')") ncfrdf

! for all possible unique type-to-type pairs

  ll=0
  Do ia=1,ntpatm
     Do ib=ia,ntpatm

! number of the interaction by its rdf key

        kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then
           ll=ll+1

           If (idnode == 0) Then
              Write(nrite,"(/,' g(r)  :',2(1x,a8),/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(ia),unqatm(ib)
              Write(nrdfdt,'(2a8)') unqatm(ia),unqatm(ib)
           End If

! global sum of data on all nodes

           If (mxnode > 1) Call gsum(rdf(1:mxgrdf,kk))

! normalisation factor

           factor=volm*dens(ia)*dens(ib)*Real(ncfrdf,wp)
           If (ia == ib) factor=factor*0.5_wp ! *(1.0_wp-1.0_wp/numtyp(ia))

! running integration of rdf

           sum=0.0_wp

! loop over distances

           zero=.true.
           Do i=1,mxgrdf
              If (zero .and. i < (mxgrdf-3)) zero=(rdf(i+2,kk) <= 0.0_wp)

              gofr= rdf(i,kk)/factor
              sum = sum + gofr*dens(ib)

              rrr = (Real(i,wp)-0.5_wp)*delr
              dvol= fourpi*delr*(rrr**2+delr**2/12.0_wp)
              gofr= gofr/dvol

! zero it if < pdfzero

              If (gofr < pdfzero) Then
                 gofr1 = 0.0_wp
              Else
                 gofr1 = gofr
              End If

              If (sum < pdfzero) Then
                 sum1 = 0.0_wp
              Else
                 sum1 = sum
              End If

! print out information

              If (idnode == 0) Then
                 If (.not.zero) Write(nrite,"(f10.4,1p,2e14.6)") rrr,gofr1,sum1
                 Write(nrdfdt,"(1p,2e14.6)") rrr,gofr
              End If

! We use the non-normalised tail-truncated PDF version,
! gofr1 (not gofr) in order to exclude the nearly-zero
! gofr noise in PMF, otherwise the PMF = -ln(pdfzero)
! would have poorly-defined noisy "borders/walls"

              dstdrdf(i,kk) = gofr1
           End Do
        Else
           dstdrdf(:,kk) = 0
        End If

     End Do
  End Do

  If (idnode == 0) Close(Unit=nrdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=npdgdt, File='RDFPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',delr*mxgrdf,mxgrdf,delr,ll, &
          '   conversion factor(kT -> energy units) =',kT2engo

     Open(Unit=npdfdt, File='RDFTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',dgrid*ngrid,ngrid,dgrid,ll, &
          '   conversion factor(kT -> energy units) =',kT2engo
  End If

! loop over all valid PDFs

  Do ia=1,ntpatm
     Do ib=ia,ntpatm

! number of the interaction by its rdf key

        kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then

           If (idnode == 0) Then
              Write(nrite,"(/,' g(r)  :',2(1x,a8),/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(ia),unqatm(ib)
              Write(npdgdt,'(/,a2,2a8)') '# ',unqatm(ia),unqatm(ib)
              Write(npdfdt,'(/,a2,2a8)') '# ',unqatm(ia),unqatm(ib)
           End If

! Smoothen and get derivatives

! RDFs -> 1 at long distances, so we do not shift the PMFs
! but instead put a cap over those, -Log(pdfzero) - hence,
! the upper bound for the PMF due to the zero RDF values

           fed0  = -Log(pdfzero)
           dfed0 = 10.0_wp
           dfed  = 10.0_wp

           Do ig=1,mxgrdf
              tmp = Real(ig,wp)-0.5_wp
              rrr = tmp*delr

              If (dstdrdf(ig,kk) > zero_plus) Then
                 fed = -Log(dstdrdf(ig,kk)+pdfzero)
                 If (fed0 <= zero_plus ) Then
                    fed0 = fed
                    fed  = fed0
                 End If

                 If (ig < mxgrdf-1) Then
                    If (dstdrdf(ig+1,kk) <= zero_plus .and. dstdrdf(ig+2,kk) > zero_plus) &
                         dstdrdf(ig+1,kk) = 0.5_wp*(dstdrdf(ig,kk)+dstdrdf(ig+2,kk))
                 End If
              Else
                 fed = fed0
              End If

              If      (ig == 1) Then
                 If      (dstdrdf(ig,kk) > zero_plus .and. dstdrdf(ig+1,kk) > zero_plus) Then
                    dfed = Log(dstdrdf(ig+1,kk)/dstdrdf(ig,kk))
                 Else If (dfed > 0.0_wp) Then
                    dfed = dfed0
                 Else
                    dfed =-dfed0
                 End If
              Else If (ig == mxgrdf) Then
                 If      (dstdrdf(ig,kk) > zero_plus .and. dstdrdf(ig-1,kk) > zero_plus) Then
                    dfed = Log(dstdrdf(ig,kk)/dstdrdf(ig-1,kk))
                 Else If (dfed > 0.0_wp) Then
                    dfed = dfed0
                 Else
                    dfed =-dfed0
                 End If
              Else If (dstdrdf(ig-1,kk) > zero_plus) Then
                 If (dstdrdf(ig+1,kk) > zero_plus) Then
                    dfed = 0.5_wp*(Log(dstdrdf(ig+1,kk)/dstdrdf(ig-1,kk)))
                 Else
                    dfed = 0.5_wp*Log(dstdrdf(ig-1,kk))
                 End If
              Else If (dstdrdf(ig+1,kk) > zero_plus) Then
                 dfed =-0.5_wp*Log(dstdrdf(ig+1,kk))
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

! Define edges

           pmf(0)         = 2.0_wp*pmf(1)        -pmf(2)
           vir(0)         = 2.0_wp*vir(1)        -vir(2)
           pmf(mxgrdf+1) = 2.0_wp*pmf(mxgrdf)  -pmf(mxgrdf-1)
           vir(mxgrdf+1) = 2.0_wp*vir(mxgrdf)  -vir(mxgrdf-1)
           pmf(mxgrdf+2) = 2.0_wp*pmf(mxgrdf+1)-pmf(mxgrdf)
           vir(mxgrdf+2) = 2.0_wp*vir(mxgrdf+1)-vir(mxgrdf)

! resample using 3pt interpolation

           Do ig=1,ngrid
              rrr = Real(ig,wp)*dgrid
              ll = Int(rrr/delr)

! +0.5_wp due to half-a-bin shift in the original (bin-centered) grid

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
        End If
     End Do
  End Do

  If (idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstdrdf,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rdf_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rdf_compute
