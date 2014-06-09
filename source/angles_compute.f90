Subroutine angles_compute(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating angles distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : pi,nrite,npdfdt,npdgdt,mxgang1,engunit,boltz,zero_plus
  Use site_module,   Only : unqatm
  Use config_module, Only : cfgname
  Use angles_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: temp

  Logical           :: zero
  Integer           :: fail,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,delth,rdlth,factor,factor1,rad2dgr,dgr2rad, &
                       theta,sinth,rsint,pdfang,sum,pdfang1,sum1,fed0,fed,dfed,dfed0,tmp

  Real( Kind = wp ), Allocatable :: dstdang(:,:)

  fail = 0
  Allocate (dstdang(0:mxgang1,1:ldfang(0)), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'angles_compute - allocation failure, node: ', idnode
     Call error(0)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! conversion: radians <-> degrees (affects not only angle units but also force units!)

  rad2dgr = 180.0_wp/pi
  dgr2rad = pi/180.0_wp

! grid interval for pdf tables

  delth = pi/Real(mxgang1,wp)
  rdlth = Real(mxgang1,wp)/180.0_wp

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,ldfang(0)
     If (typang(0,i) > 0) Then
        kk=kk+1
        ll=ll+typang(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(ncfang,wp)

  If (idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'ANGLES PDFs'
     Write(nrite,'(/,1x,a,4(1x,i10))') '# types, bins, cutoff, frames:',kk,mxgang1,180,ncfang
  End If

! open RDF file and write headers

  If (idnode == 0) Then
     Open(Unit=npdfdt, File='ANGDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,4(1x,i10))') '# types, bins, cutoff, frames:',kk,mxgang1,180,ncfang
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  Pn_ang(Theta)  Pn_ang(Theta)/Sin(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfang(0)
     If (typang(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(nrite,'(/,1x,a,3(a8,1x),2(i10,1x))')  'id, type, totals: ', &
                unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i)
           Write(nrite,'(/,1x,a,f8.5,/)') 'Theta(degrees)  P_ang(Theta)  Sum_P_ang(Theta)   @   dTheta_bin = ',delth*rad2dgr

           Write(npdfdt,'(/,a,3(a8,1x),2(i10,1x))') '# id, type, totals: ', &
                unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i)
        End If

! global sum of data on all nodes

        If (mxnode > 1) Call gsum(dstang(1:mxgang1,i))

! factor in degeneracy (first, pdfang is normalised to unity)

        factor1=factor/Real(typang(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxgang1
           If (zero .and. ig < (mxgang1-3)) zero=(dstang(ig+2,i) <= 0.0_wp)

           pdfang = dstang(ig,i)*factor1
           sum = sum + pdfang

! null it if < 1.0e-4_wp

           If (pdfang < 1.0e-4_wp) Then
              pdfang1 = 0.0_wp
           Else
              pdfang1 = pdfang
           End If

           If (sum < 1.0e-4_wp) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           theta = (Real(ig,wp)-0.5_wp)*delth

! Jacobian rsint = 1.0_wp/sinth must only be used for the PMF calculations below,
! which require renormalising: dstdang(ig,i) = pdfang1*rsint

           sinth = Max(1.0e-10_wp,Sin(theta))
           rsint = 1.0_wp/sinth

! now pdfang is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfang = pdfang*rdlth

! print out information

           theta  = theta*rad2dgr
           If (idnode == 0) Then
              If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") theta,pdfang1,sum1
              Write(npdfdt,"(f11.5,1p,2e14.6)") theta,pdfang,pdfang*rsint
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF/sinth)
! would have poorly-defined noisy "borders/walls"

           dstdang(ig,i) = pdfang1*rsint ! PDFs density
        End Do
     Else
        dstdang(:,i) = 0 ! PDFs density
     End If
  End Do

  If (idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=npdgdt, File='PMFANG', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,f11.5,2i10,a,e15.7)') '# ',delth,mxgang1,kk,' conversion factor: kT -> energy units =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfang(0)
     If (typang(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Write(npdgdt,'(/,a,3(a8,1x),2(i10,1x))') '# id, type, degeneracy: ', &
           unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i)

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxgang1
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth

           If (dstdang(ig,i) > zero_plus) Then
              fed = -Log(dstdang(ig,i))-fed0
              If (fed0 <= zero_plus ) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < mxgang1-1) Then
                 If (dstdang(ig+1,i) <= zero_plus .and. dstdang(ig+2,i) > zero_plus) &
                    dstdang(ig+1,i) = 0.5_wp*(dstdang(ig,i)+dstdang(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstdang(ig,i) > zero_plus .and. dstdang(ig+1,i) > zero_plus) Then
                 dfed = Log(dstdang(ig+1,i)/dstdang(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == mxgang1) Then
              If      (dstdang(ig,i) > zero_plus .and. dstdang(ig-1,i) > zero_plus) Then
                 dfed = Log(dstdang(ig,i)/dstdang(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstdang(ig-1,i) > zero_plus) Then
              If (dstdang(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstdang(ig+1,i)/dstdang(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstdang(ig-1,i))
              End If
           Else If (dstdang(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstdang(ig+1,i))
           Else If (dfed > 0.0_wp) Then
              dfed = dfed0
           Else
              dfed =-dfed0
           End If

! Print
           If (idnode == 0) &
              Write(npdgdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do
     End If
  End Do

  If (idnode == 0) Close(Unit=npdgdt)

  Deallocate (dstdang, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'angles_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine angles_compute
