Subroutine inversions_compute(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating inversions distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : pi,nrite,npdfdt,npdgdt,mxginv1,engunit,boltz,zero_plus
  Use site_module,   Only : unqatm
  Use config_module, Only : cfgname
  Use inversions_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: temp

  Logical           :: zero
  Integer           :: fail,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,delth,rdlth,factor,factor1,rad2dgr,dgr2rad, &
                       theta,pdfinv,sum,pdfinv1,sum1,fed0,fed,dfed,dfed0,tmp

  Real( Kind = wp ), Allocatable :: dstdinv(:,:)

  fail = 0
  Allocate (dstdinv(0:mxginv1,1:ldfinv(0)), Stat = fail)
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
  rdlth = Real(mxginv1,wp)/360.0_wp

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

  If (idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'inversions PDFs'
     Write(nrite,'(/,1x,a,4(1x,i10))') '# types, bins, cutoff, frames:',kk,mxginv1,180,ncfinv
  End If

! open RDF file and write headers

  If (idnode == 0) Then
     Open(Unit=npdfdt, File='invDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,4(1x,i10))') '# types, bins, cutoff, frames:',kk,mxginv1,180,ncfinv
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  Pn_inv(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfinv(0)
     If (typinv(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(nrite,'(/,1x,a,4(a8,1x),2(i10,1x))') 'id, type, totals: ', &
                unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)
           Write(nrite,'(/,1x,a,f8.5)') 'Theta(degrees)  P_inv(Theta)  Sum_P_inv(Theta)   @   dTheta_bin = ',delth*rad2dgr

           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x))') '# id, type, totals: ', &
                unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)
        End If

! global sum of data on all nodes

        If (mxnode > 1) Call gsum(dstinv(1:mxginv1,i))

! factor in degeneracy (first, pdfinv is normalised to unity)

        factor1=factor/Real(typinv(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxginv1
           If (zero .and. ig < (mxginv1-3)) zero=(dstinv(ig+2,i) <= 0.0_wp)

           pdfinv = dstinv(ig,i)*factor1
           sum = sum + pdfinv

! null it if < 1.0e-5_wp

           If (pdfinv < 1.0e-5_wp) Then
              pdfinv1 = 0.0_wp
           Else
              pdfinv1 = pdfinv
           End If

           If (sum < 1.0e-5_wp) Then
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
        dstdinv(:,i) = 0 ! PDFs density
     End If
  End Do

  If (idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=npdgdt, File='PMFinv', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,f11.5,2i10,a,e15.7)') '# ',delth,mxginv1,kk,' conversion factor: kT -> energy units =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfinv(0)
     If (typinv(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Write(npdgdt,'(/,a,4(a8,1x),2(i10,1x))') '# id, type, degeneracy: ', &
           unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxginv1
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth

           If (dstdinv(ig,i) > zero_plus) Then
              fed = -Log(dstdinv(ig,i))-fed0
              If (fed0 <= zero_plus ) Then
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

! Print

           If (idnode == 0) Write(npdgdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do
     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=npdgdt)
  End If

  Deallocate (dstdinv, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'inversions_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine inversions_compute
