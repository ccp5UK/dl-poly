Subroutine inversions_compute(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating inversions distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : pi,boltz,delth_max,nrite,npdfdt,npdgdt, &
                            mxginv,mxginv1,engunit,zero_plus
  Use site_module,   Only : unqatm
  Use config_module, Only : cfgname
  Use inversions_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: temp

  Logical           :: zero
  Integer           :: fail,ngrid,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,rad2dgr,dgr2rad,delth,rdlth,dgrid,factor,pdfzero, &
                       factor1,theta,pdfinv,sum,pdfinv1,sum1,                    &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstdinv(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  fail = 0
  Allocate (dstdinv(0:mxginv1,1:ldfinv(0)),pmf(0:mxginv1+2),vir(0:mxginv1+2), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'inversions_compute - allocation failure, node: ', idnode
     Call error(0)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! conversion: radians <-> degrees (affects not only angle units but also force units!)

  rad2dgr = 180.0_wp/pi
  dgr2rad = pi/180.0_wp

! grid interval for pdf tables

  delth = pi/Real(mxginv1,wp)
  rdlth = Real(mxginv1,wp)/180.0_wp

! resampling grid and grid interval for pmf tables

  ngrid = Max(Nint(360.0_wp/delth_max),mxginv1,mxginv-4)
  dgrid = pi/Real(ngrid,wp)

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,ldfinv(0)
     If (typinv(0,i) > 0) Then
        kk=kk+1
        ll=ll+typinv(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(ncfinv,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-5_wp

  If (idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'INVERSIONS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)'
     Write(nrite,'(/,1x,a,5(1x,i10))') &
            '# bins, cutoff, frames, types: ',mxginv1,180,ncfinv,kk,ll
  End If

! open RDF file and write headers

  If (idnode == 0) Then
     Open(Unit=npdfdt, File='INVDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a)') '# INVERSIONS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
     Write(npdfdt,'(a,4(1x,i10))') '# bins, cutoff, frames, types: ',mxginv1,180,ncfinv,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfinv(0)
     If (typinv(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(nrite,'(/,1x,a,4(a8,1x),2(i10,1x))') 'type, index, instances: ', &
                unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)
           Write(nrite,'(/,1x,a,f8.5)') 'Theta(degrees)  P_inv(Theta)  Sum_P_inv(Theta)   @   dTheta_bin = ',delth*rad2dgr

           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)
        End If

! global sum of data on all nodes

        If (mxnode > 1) Call gsum(dstinv(1:mxginv1,i))

! factor in instances (first, pdfinv is normalised to unity)

        factor1=factor/Real(typinv(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxginv1
           If (zero .and. ig < (mxginv1-3)) zero=(dstinv(ig+2,i) <= 0.0_wp)

           pdfinv = dstinv(ig,i)*factor1
           sum = sum + pdfinv

! null it if < pdfzero

           If (pdfinv < pdfzero) Then
              pdfinv1 = 0.0_wp
           Else
              pdfinv1 = pdfinv
           End If

           If (sum < pdfzero) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           theta = (Real(ig,wp)-0.5_wp)*delth

! now pdfinv is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfinv = pdfinv*rdlth

! print out information

           theta  = theta*rad2dgr
           If (idnode == 0) Then
              If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") theta,pdfinv1,sum1
              Write(npdfdt,"(f11.5,1p,e14.6)") theta,pdfinv
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

           dstdinv(ig,i) = pdfinv1 ! PDFs density
        End Do
     Else
        dstdinv(:,i) = 0.0_wp ! PDFs density
     End If
  End Do

  If (idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=npdgdt, File='INVPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',mxginv1,delth*mxginv1*rad2dgr,delth*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo

     Open(Unit=npdfdt, File='INVPMF', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',ngrid,dgrid*ngrid*rad2dgr,dgrid*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfinv(0)
     If (typinv(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(npdgdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)), &
                unqatm(typinv(4,i)),j,typinv(0,i),' (type, index, instances)'
           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)), &
                unqatm(typinv(4,i)),j,typinv(0,i),' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxginv1
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth

           If (dstdinv(ig,i) > zero_plus) Then
              fed = -Log(dstdinv(ig,i))-fed0
              If (fed0 <= zero_plus) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < mxginv1-1) Then
                 If (dstdinv(ig+1,i) <= zero_plus .and. dstdinv(ig+2,i) > zero_plus) &
                    dstdinv(ig+1,i) = 0.5_wp*(dstdinv(ig,i)+dstdinv(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstdinv(ig,i) > zero_plus .and. dstdinv(ig+1,i) > zero_plus) Then
                 dfed = Log(dstdinv(ig+1,i)/dstdinv(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == mxginv1) Then
              If      (dstdinv(ig,i) > zero_plus .and. dstdinv(ig-1,i) > zero_plus) Then
                 dfed = Log(dstdinv(ig,i)/dstdinv(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstdinv(ig-1,i) > zero_plus) Then
              If (dstdinv(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstdinv(ig+1,i)/dstdinv(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstdinv(ig-1,i))
              End If
           Else If (dstdinv(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstdinv(ig+1,i))
           Else If (dfed > 0.0_wp) Then
              dfed = dfed0
           Else
              dfed =-dfed0
           End If

           pmf(ig) = fed
           vir(ig) = dfed

! Print

           If (idnode == 0) Write(npdgdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do

! Define edges

        pmf(0)         = 2.0_wp*pmf(1)        -pmf(2)
        vir(0)         = 2.0_wp*vir(1)        -vir(2)
        pmf(mxginv1+1) = 2.0_wp*pmf(mxginv1)  -pmf(mxginv1-1)
        vir(mxginv1+1) = 2.0_wp*vir(mxginv1)  -vir(mxginv1-1)
        pmf(mxginv1+2) = 2.0_wp*pmf(mxginv1+1)-pmf(mxginv1)
        vir(mxginv1+2) = 2.0_wp*vir(mxginv1+1)-vir(mxginv1)

! resample using 3pt interpolation

        Do ig=1,ngrid
           theta = Real(ig,wp)*dgrid
           ll = Int(theta/delth)

! +0.5_wp due to half-a-bin shift in the original data

           coef = theta/delth-Real(ll,wp)+0.5_wp

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
              Write(npdfdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do
     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstdinv,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'inversions_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine inversions_compute
