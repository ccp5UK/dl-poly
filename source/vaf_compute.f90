Subroutine vaf_compute(lvafav,tstep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating velocity autocorrelation
! functions from accumulated data
!
! copyright - daresbury laboratory
! author    - m.a.seaton june 2014
! amended   - i.t.todorov july 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,     Only : idnode
  Use setup_module,     Only : nrite,mxatyp,zero_plus
  Use site_module,      Only : numtypnf
  Use greenkubo_module, Only : nsvaf,vafcount,vafstep,vaftime,vaf

  Implicit None

  Logical,           Intent( In    ) :: lvafav
  Real( Kind = wp ), Intent( In    ) :: tstep

  Integer           :: i
  Real( Kind = wp ) :: factor,gvaf,ovaf,time0,timei,numt

  If (idnode == 0) Then
     If (lvafav) Then
        Write(nrite,"(/,/,12x,'VELOCITY AUTOCORRELATION FUNCTIONS',/,/, &
             & 'calculated using ',i8,' samples')") Nint(vafcount)
     Else
        Write(nrite,"(/,/,12x,'VELOCITY AUTOCORRELATION FUNCTIONS',/,/, &
             & 'calculated using sample ',i8,' starting at ',f10.4,' ps')") Nint(vafcount), vaftime(0)
     End If
  End If

  time0 = vaftime(0)
  numt = Sum(numtypnf(1:mxatyp))
  factor = 1.0_wp/Sum(vaf(0,1:mxatyp))
  ovaf = Sum(vaf(0,1:mxatyp))/Real(numt,Kind=wp)
  If (lvafav) ovaf = ovaf/vafcount

  If (idnode == 0) Write(nrite,"(12x,'absolute value at origin (3kT/m) = ',1p,e16.8)") ovaf

! loop over time steps

  Do i=0,nsvaf

! determine time

     If (lvafav) Then
        timei = tstep*Real(i,Kind=wp)
     Else
        timei = vaftime(i)-time0
     End If

! null it if number of particles is zero

     If (numt > zero_plus) Then
        gvaf = Sum(vaf(i,1:mxatyp))*factor
     Else ! THIS CAN NEVER HAPPEN AS TOTAL NON-FROZEN PARTICLES IS ALWAYS > 1
        gvaf = 0.0_wp
     End If

! print out information

     If (idnode == 0) Write(nrite,"(12x,f10.4,1p,e14.6)") timei,gvaf

  End Do

End Subroutine vaf_compute
