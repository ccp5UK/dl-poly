Subroutine rdf_compute(lpana,rcut,temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating radial distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - t.forester & i.t.todorov march 2016
! contrib   - a.v.brukhno january 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : fourpi,boltz,delr_max,nrite,nrdfdt,npdfdt,npdgdt, &
                            mxgrdf,engunit,zero_plus
  Use site_module,   Only : ntpatm,unqatm,numtyp,dens
  Use config_module, Only : cfgname,volm
  Use rdf_module

  Implicit None

  Logical          , Intent( In    ) :: lpana
  Real( Kind = wp ), Intent( In    ) :: rcut,temp

  Logical           :: zero
  Integer           :: fail,ngrid,i,ia,ib,kk,ig,ll
  Real( Kind = wp ) :: kT2engo,delr,rdlr,dgrid,pdfzero,      &
                       factor1,rrr,dvol,gofr,gofr1,sum,sum1, &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstdrdf(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  If (lpana) Then
     fail = 0
     Allocate (dstdrdf(0:mxgrdf,1:ntprdf),pmf(0:mxgrdf+2),vir(0:mxgrdf+2), Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'rdf_compute - allocation failure, node: ', idnode
        Call error(0)
     End If
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! grid interval for rdf tables

  delr = rcut/Real(mxgrdf,wp)
  rdlr = 1.0_wp/delr

! resampling grid and grid interval for rdf tables

  ngrid = Max(Nint(rcut/delr_max),mxgrdf)
  dgrid = rcut/Real(ngrid,wp)

  If (idnode == 0) Write(nrite,"(/,/,12x,'RADIAL DISTRIBUTION FUNCTIONS',/,/, &
     & ' calculated using ',i8,' configurations')") ncfrdf

! open RDF file and Write headers

  If (idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf
  End If

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-6_wp

! for all possible unique type-to-type pairs

  Do ia=1,ntpatm
     Do ib=ia,ntpatm

! number of the interaction by its rdf key

        kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then
           If (idnode == 0) Then
              Write(nrite,"(/,' g(r)  :',2(1x,a8),/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(ia),unqatm(ib)
              Write(nrdfdt,'(2a8)') unqatm(ia),unqatm(ib)
           End If

! global sum of data on all nodes

           If (mxnode > 1) Call gsum(rdf(1:mxgrdf,kk))

! normalisation factor

           factor1=volm*dens(ia)*dens(ib)*Real(ncfrdf,wp)
           If (ia == ib) factor1=factor1*0.5_wp*(1.0_wp-1.0_wp/numtyp(ia))

! running integration of rdf

           sum=0.0_wp

! loop over distances

           zero=.true.
           Do i=1,mxgrdf
              If (zero .and. i < (mxgrdf-3)) zero=(rdf(i+2,kk) <= 0.0_wp)

              gofr= rdf(i,kk)/factor1
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

! We use the non-normalised tail-truncated RDF version,
! rdf...1 (not pdf...) in order to exclude the nearly-zero
! rdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

                If (lpana) dstdrdf(i,kk) = gofr1 ! RDFs density
           End Do
        Else
           If (lpana) dstdrdf(:,kk) = 0.0_wp ! RDFs density
        End If

     End Do
  End Do

  If (idnode == 0) Close(Unit=nrdfdt)

! Only when PDF analysis is requested

  If (lpana) Then

! open PDF files and write headers

     If (idnode == 0) Then
        Open(Unit=npdgdt, File='VDWPMF', Status='replace')
        Write(npdgdt,'(a)') '# '//cfgname
        Write(npdgdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',delr*Real(mxgrdf,wp),mxgrdf,delr,ntprdf, &
             '   conversion factor(kT -> energy units) =',kT2engo

        Open(Unit=npdfdt, File='VDWTAB', Status='replace')
        Write(npdfdt,'(a)') '# '//cfgname
        Write(npdfdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',dgrid*Real(ngrid,wp),ngrid,dgrid,ntprdf, &
          '   conversion factor(kT -> energy units) =',kT2engo
     End If

! loop over all valid RDFs

     Do ia=1,ntpatm
        Do ib=ia,ntpatm

! number of the interaction by its rdf key

           kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

           If (kk > 0 .and. kk <= ntprdf) Then
              If (idnode == 0) Then
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
                    fed = -Log(dstdrdf(ig,kk)+pdfzero) !-fed0
                    If (fed0 <= zero_plus) Then
                       fed0 = fed
                       fed  = fed0
!                       fed = 0.0_wp
                    End If

                    If (ig < mxgrdf-1) Then
                       If (dstdrdf(ig+1,kk) <= zero_plus .and. dstdrdf(ig+2,kk) > zero_plus) &
                          dstdrdf(ig+1,kk) = 0.5_wp*(dstdrdf(ig,kk)+dstdrdf(ig+2,kk))
                    End If
                 Else
                    fed = fed0
!                    fed = 0.0_wp
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

              pmf(0)        = 2.0_wp*pmf(1)       -pmf(2)
              vir(0)        = 2.0_wp*vir(1)       -vir(2)
              pmf(mxgrdf+1) = 2.0_wp*pmf(mxgrdf)  -pmf(mxgrdf-1)
              vir(mxgrdf+1) = 2.0_wp*vir(mxgrdf)  -vir(mxgrdf-1)
              pmf(mxgrdf+2) = 2.0_wp*pmf(mxgrdf+1)-pmf(mxgrdf)
              vir(mxgrdf+2) = 2.0_wp*vir(mxgrdf+1)-vir(mxgrdf)

! resample using 3pt interpolation

              Do ig=1,ngrid
                 rrr = Real(ig,wp)*dgrid
                 ll = Int(rrr*rdlr)

! +0.5_wp due to half-a-bin shift in the original (bin-centered) grid

                 coef = rrr*rdlr-Real(ll,wp)+0.5_wp

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

                 If (idnode == 0) &
                    Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*rrr*rdlr
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
  End If

End Subroutine rdf_compute

Subroutine usr_compute()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating radial distribution function
! from accumulated data for umbrella sampled two COMs separation (ushr)
!
! to be used in exernal_field_apply & statistics_result
!
! copyright - daresbury laboratory
! author    - i.t.todorov & a.brukhno november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : fourpi,nrdfdt,mxgusr
  Use config_module, Only : cfgname,volm
  Use rdf_module,    Only : ncfusr,rusr,usr

  Implicit None

  Integer           :: i
  Real( Kind = wp ) :: delr,rdlr,factor1,rrr,dvol,gofr,sum0,sum1

! grid interval for rdf tables

  delr = rusr/Real(mxgusr,wp)
  rdlr = 1.0_wp/delr

! open RDF file and Write headers

  If (idnode == 0) Then
     Open(Unit=nrdfdt, File='USRDAT', Status='replace')
     Write(nrdfdt,'(2a)') '# '//cfgname
     Write(nrdfdt,'(a)')  "# RDF for the two fragments' COMs (umbrella sampling)"
     Write(nrdfdt,'(a,i10,f12.6,i10,e15.6,/)') '# bins, cutoff, frames, volume: ',mxgusr,rusr,ncfusr,volm
     Write(nrdfdt,'(a)') '#'
  End If

! global sum of data on all nodes

  If (mxnode > 1) Call gsum(usr(1:mxgusr))

! get normalisation factor

  factor1 = Sum(usr(1:mxgusr))

! running integration of rdf

  sum0 = 0.0_wp
  sum1 = 0.0_wp

! loop over distances

  Do i=1,mxgusr
     gofr = usr(i)/factor1
     sum0 = sum0 + gofr

     rrr  = (Real(i,wp)-0.5_wp)*delr
     dvol = fourpi*delr*(rrr**2+delr**2/12.0_wp)
     gofr = gofr*volm/dvol
     sum1 = sum1 + gofr

! print out information

     If (idnode == 0) Write(nrdfdt,"(1p,4e15.6)") rrr,gofr
!     If (idnode == 0) Write(nrdfdt,"(1p,4e15.6)") rrr,gofr,sum0,sum1
  End Do

  If (idnode == 0) Close(Unit=nrdfdt)

! distribute usr between nodes

  If (mxnode > 1) usr(:) = usr(:) / Real(mxnode,wp)

End Subroutine usr_compute
