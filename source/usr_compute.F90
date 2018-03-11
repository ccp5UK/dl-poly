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

  Use kinds, only : wp
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : fourpi,nrdfdt,mxgusr
  Use configuration, Only : cfgname,volm
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
