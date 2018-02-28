Subroutine dihedrals_compute(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating dihedrals distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : pi,twopi,boltz,delth_max,nrite,npdfdt,npdgdt, &
                            mxgdih,mxgdih1,engunit,zero_plus
  Use site_module,   Only : unqatm
  Use config_module, Only : cfgname
  Use dihedrals_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: temp

  Logical           :: zero
  Integer           :: fail,ngrid,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,rad2dgr,dgr2rad,delth,rdlth,dgrid,factor,pdfzero, &
                       factor1,theta,pdfdih,sum,pdfdih1,sum1,                    &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstddih(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  fail = 0
  Allocate (dstddih(0:mxgdih1,1:ldfdih(0)),pmf(0:mxgdih1+2),vir(0:mxgdih1+2), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'dihedrals_compute - allocation failure, node: ', idnode
     Call error(0)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! conversion: radians <-> degrees (affects not only dihedral units but also force units!)

  rad2dgr = 180.0_wp/pi
  dgr2rad = pi/180.0_wp

! grid interval for pdf/pmf tables

  delth = twopi/Real(mxgdih1,wp)
  rdlth = Real(mxgdih1,wp)/360.0_wp

! resampling grid and grid interval for pmf tables

  ngrid = Max(Nint(360.0_wp/delth_max),mxgdih1,mxgdih-4)
  dgrid = twopi/Real(ngrid,wp)

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,ldfdih(0)
     If (typdih(0,i) > 0) Then
        kk=kk+1
        ll=ll+typdih(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(ncfdih,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-5_wp

  If (idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'DIHEDRALS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)'
     Write(nrite,'(/,1x,a,i10,1x,a,2(i0,a),3(1x,i10))') &
           '# bins, range, frames, types: ',mxgdih1,'[',-180,',',180,']',ncfdih,kk,ll
  End If

! open RDF file and write headers

  If (idnode == 0) Then
     Open(Unit=npdfdt, File='DIHDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a)') '# DIHEDRALS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
     Write(npdfdt,'(a,4(1x,i10))') '# bins, cutoff, frames, types: ',mxgdih1,360,ncfdih,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfdih(0)
     If (typdih(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(nrite,'(/,1x,a,4(a8,1x),2(i10,1x))') 'type, index, instances: ', &
                unqatm(typdih(1,i)),unqatm(typdih(2,i)),unqatm(typdih(3,i)),unqatm(typdih(4,i)),j,typdih(0,i)
           Write(nrite,'(/,1x,a,f8.5)') 'Theta(degrees)  P_dih(Theta)  Sum_P_dih(Theta)   @   dTheta_bin = ',delth*rad2dgr

           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                unqatm(typdih(1,i)),unqatm(typdih(2,i)),unqatm(typdih(3,i)),unqatm(typdih(4,i)),j,typdih(0,i)
        End If

! global sum of data on all nodes

        If (mxnode > 1) Call gsum(dstdih(1:mxgdih1,i))

! factor in instances (first, pdfdih is normalised to unity)

        factor1=factor/Real(typdih(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxgdih1
           If (zero .and. ig < (mxgdih1-3)) zero=(dstdih(ig+2,i) <= 0.0_wp)

           pdfdih = dstdih(ig,i)*factor1
           sum = sum + pdfdih

! null it if < pdfzero

           If (pdfdih < pdfzero) Then
              pdfdih1 = 0.0_wp
           Else
              pdfdih1 = pdfdih
           End If

           If (sum < pdfzero) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           theta = (Real(ig,wp)-0.5_wp)*delth-pi

! now pdfdih is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfdih = pdfdih*rdlth

! print out information

           theta  = theta*rad2dgr
           If (idnode == 0) Then
              If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") theta,pdfdih1,sum1
              Write(npdfdt,"(f11.5,1p,e14.6)") theta,pdfdih
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

           dstddih(ig,i) = pdfdih1 ! PDFs density
        End Do
     Else
        dstddih(:,i) = 0.0_wp ! PDFs density
     End If
  End Do

  If (idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (idnode == 0) Then
     Open(Unit=npdgdt, File='DIHPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',mxgdih1,delth*Real(mxgdih1,wp)*rad2dgr,delth*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo

     Open(Unit=npdfdt, File='DIHTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',ngrid,dgrid*Real(ngrid,wp)*rad2dgr,dgrid*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfdih(0)
     If (typdih(0,i) > 0) Then
        j=j+1

        If (idnode == 0) Then
           Write(npdgdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typdih(1,i)),unqatm(typdih(2,i)),unqatm(typdih(3,i)), &
                unqatm(typdih(4,i)),j,typdih(0,i),' (type, index, instances)'
           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typdih(1,i)),unqatm(typdih(2,i)),unqatm(typdih(3,i)), &
                unqatm(typdih(4,i)),j,typdih(0,i),' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxgdih1
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth-pi

           If (dstddih(ig,i) > zero_plus) Then
              fed = -Log(dstddih(ig,i))-fed0
              If (fed0 <= zero_plus) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < mxgdih1-1) Then
                 If (dstddih(ig+1,i) <= zero_plus .and. dstddih(ig+2,i) > zero_plus) &
                    dstddih(ig+1,i) = 0.5_wp*(dstddih(ig,i)+dstddih(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstddih(ig,i) > zero_plus .and. dstddih(ig+1,i) > zero_plus) Then
                 dfed = Log(dstddih(ig+1,i)/dstddih(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == mxgdih1) Then
              If      (dstddih(ig,i) > zero_plus .and. dstddih(ig-1,i) > zero_plus) Then
                 dfed = Log(dstddih(ig,i)/dstddih(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstddih(ig-1,i) > zero_plus) Then
              If (dstddih(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstddih(ig+1,i)/dstddih(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstddih(ig-1,i))
              End If
           Else If (dstddih(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstddih(ig+1,i))
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

! Cyclic grid

        pmf(0)         = 0.5_wp*(pmf(1)-pmf(mxgdih1))
        vir(0)         = 0.5_wp*(vir(1)-vir(mxgdih1))
        pmf(mxgdih1+1) = pmf(0)
        vir(mxgdih1+1) = vir(0)
        pmf(mxgdih1+2) = pmf(1)
        vir(mxgdih1+2) = vir(1)

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
              Write(npdfdt,"(f11.5,1p,2e14.6)") (theta-pi)*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do
     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstddih,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'dihedrals_compute - deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine dihedrals_compute
