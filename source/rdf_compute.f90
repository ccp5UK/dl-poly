Subroutine rdf_compute(rcut)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating radial distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - t.forester march 1994
! amended   - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : mxgrdf,nrite,nrdfdt,pi
  Use site_module,   Only : ntpatm,unqatm,dens
  Use config_module, Only : cfgname,volm
  Use rdf_module


  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rcut

  Logical           :: zero
  Integer           :: i,ia,ib,kk
  Real( Kind = wp ) :: delr,dvol,factor,gofr,gofr1,rrr,sum,sum1

  If (idnode == 0) Write(nrite,"(/,/,12X,'RADIAL DISTRIBUTION FUNCTIONS',/,/, &
     & 'calculated using ',i10,' configurations')") ncfrdf

! open RDF file and Write headers

  If (idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf
  End If

! grid interval for rdf tables

  delr=rcut/Real(mxgrdf,wp)

! for all possible unique type-to-type pairs

  Do ia=1,ntpatm
     Do ib=ia,ntpatm

! number of the interaction by its rdf key

        kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then
           If (idnode == 0) Then
              Write(nrite,"(/,'g(r)  :',2a8,/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(ia),unqatm(ib)
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
              dvol= 4.0_wp*pi*delr*(rrr**2+delr**2/12.0_wp)
              gofr= gofr/dvol

! zero it if < 1.0e-6_wp

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
                 If (.not.zero) Write(nrite,"(f10.4,1p,2e14.6)") rrr,gofr1,sum1
                 Write(nrdfdt,"(1p,2e14.6)") rrr,gofr
              End If
           End Do
        End If

     End Do
  End Do

  If (idnode == 0) Close(Unit=nrdfdt)

End Subroutine rdf_compute
