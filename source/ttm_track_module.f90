Module ttm_track_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for utility functions for initial energy deposition
! in the electronic subsystem in the two-temperature model (ttm)
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton july 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use ttm_module
  Use ttm_utils
	
  Implicit None
	
  Real( Kind = wp ), Allocatable, Dimension (:,:,:) :: lat_U,lat_B,lat_I
  Real( Kind = wp ) :: norm

  Logical :: trackInit = .false.

Contains

  Subroutine depoinit(time)

! determine initial energy deposition to electronic system,
! both temporally and spatially

    Implicit None

    Real ( Kind = wp ), Intent( In ) :: time
    Integer, Dimension( 1:3 ) :: fail
    Character ( Len = 14 ) :: number

    fail = 0

    Allocate (lat_U (0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1), Stat = fail(1))
    Allocate (lat_B (0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1), Stat = fail(2))
    Allocate (lat_I (0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1), Stat = fail(3))

    If (Any(fail>0)) Call error(1089)
		
    lat_U(:,:,:) = 0.0_wp ! spatial deposition (eV)
    lat_B(:,:,:) = 0.0_wp ! temporal deposition of lat_U (eV)
    lat_I(:,:,:) = 0.0_wp ! sum of temporal deposition of lat_B (eV)

! spatial distribution of track

    Select Case (sdepoType)
    Case (1)
    ! Gaussian spatial deposition
      Call gaussianTrack(lat_U)
    Case (2)
    ! Constant (flat) spatial deposition
      Call uniformDist(lat_U)
    Case (3)
    ! xy-flat, z-exp spatial deposition
      Call uniformDistZexp(lat_U)
    End Select

    trackInit = .true.                           ! switch on flag indicating track initialisation is in progress
    If (depostart<=zero_plus) depostart = time   ! time (ps) when deposition starts, i.e. current time

! temporal deposition of track: calculate time normalisation factor

    Select Case (tdepoType)
!   type=1: gauss(t)
    Case (1)
    ! Gaussian temporal deposition
      norm = 1.0_wp/(sqrpi*rt2*tdepo)
      depoend = depostart+2.0_wp*tcdepo*tdepo
    Case (2)
    ! decaying exponential temporal deposition
      norm = 1.0_wp/(1.0_wp-Exp(-tcdepo))
      depoend = depostart+2.0_wp*tcdepo*tdepo
    Case (3)
    ! delta temporal deposition
      norm = 1.0_wp
      depoend = depostart
    Case (4)
    ! pulse temporal deposition
      norm = 1.0_wp/tdepo
      depoend = depostart+tdepo
    End Select

    ! report start of energy deposition

    If (idnode == 0) Then
      Write(number, '(f14.5)') depostart
      Write(nrite,"(/,6x,a,a,a,/)") &
        'electronic energy deposition starting at time = ',Trim(Adjustl(number)),' ps'
      Write(nrite,"(1x,130('-'))")
    End If

  End Subroutine depoinit

  Subroutine depoevolve(time,tstep,redtstep,redtstepmx)

! determine how deposition evolves over time

    Implicit None
    Real( Kind = wp ), Intent ( In ) :: tstep,time
    Integer, Intent ( In ) :: redtstepmx,redtstep
		
    Real( Kind = wp ) :: lat_I_max, lat_I_sum, lat_I_min
    Real( Kind = wp ) :: lat_U_min, lat_U_max, lat_U_sum, invbin
    Real( Kind = wp ) :: currenttime,adjtime,adjeng,err_tol = 0.01_wp
    Real( Kind = wp ) :: energy_diff,oldCe,newCe,start_Te,end_Te,increase
    Real( Kind = wp ) :: Ce0a, sh_Aa, Cemaxa
    Integer :: i,j,k,ijk
    Integer, Dimension( 1:3 ) :: fail = 0
    Character ( Len = 14 ) :: number
    Logical :: deposit

    ! provide atomic density corrections to heat capacities

    Ce0a   = Ce0  *Merge(cellrho,1.0_wp,ttmdyndens)
    sh_Aa  = sh_A *Merge(cellrho,1.0_wp,ttmdyndens)
    Cemaxa = Cemax*Merge(cellrho,1.0_wp,ttmdyndens)

    ! start deposition, reducing size of timestep for thermal diffusion

    currenttime = time-depostart+tstep/Real(redtstepmx,Kind=wp)*Real(redtstep,Kind=wp)

    Select Case (tdepoType)
    Case (1)
    ! Gaussian temporal deposition
      adjtime = currenttime/tdepo-tcdepo
      deposit = (currenttime<2.0_wp*tdepo*tcdepo)
      invbin = tstep/Real(redtstepmx,Kind=wp)
      adjeng = Exp(-0.5_wp*adjtime*adjtime)
    Case (2)
    ! decaying exponential temporal deposition
      adjtime = currenttime/tdepo
      deposit = (currenttime<tdepo*tcdepo)
      invbin = tstep/(tdepo*Real(redtstepmx,Kind=wp))
      adjeng = Exp(-adjtime)
    Case (3)
    ! delta temporal deposition (over single diffusion timestep)
      deposit = (currenttime<tstep/(Real(redtstepmx,Kind=wp)))
      invbin = 1.0_wp
      adjeng = 1.0_wp
    Case (4)
    ! pulse temporal deposition (over tdepo ps)
      deposit = (currenttime<tdepo)
      invbin = 1.0_wp
      adjeng = 1.0_wp
    End Select

    ! if (still) depositing energy, add to electronic temperature
    ! grid (active cells only) and adjust electronic temperatures 
    ! accordingly

    If (deposit) Then
      lat_B(:,:,:) = lat_U(:,:,:)*norm*invbin*adjeng
      Select Case (CeType)
      Case (0,4)
      ! constant specific heat capacity
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2)*k)
              lat_I(i,j,k) = lat_I(i,j,k)+lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)
              energy_diff = lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)*rvolume*eV_to_kB
              If (energy_diff>zero_plus) Then
                start_Te = eltemp(ijk,0,0,0)
                end_Te = start_Te + energy_diff/Ce0a
                eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (1,5)
      ! hyperbolic tangent specific heat capacity
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2)*k)
              lat_I(i,j,k) = lat_I(i,j,k)+lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)
              energy_diff = lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)*rvolume*eV_to_kB
              If (energy_diff>zero_plus) Then
                start_Te = eltemp(ijk,0,0,0)
                increase = Cosh(sh_B*start_Te)*Exp(sh_B*energy_diff/sh_Aa)
                ! using equivalent function: Acosh(x)=Log(x+Sqrt((x-1.0)*(x+1.0)))
                end_Te = Log(increase+Sqrt((increase-1.0_wp)*(increase+1.0_wp)))/sh_B
                eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (2,6)
      ! linear specific heat capacity to Fermi temperature
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2)*k)
              lat_I(i,j,k) = lat_I(i,j,k)+lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)
              energy_diff = lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)*rvolume*eV_to_kB
              If (energy_diff>zero_plus) Then
                start_Te = eltemp(ijk,0,0,0)
                If (start_Te>=Tfermi) Then
                  end_Te = start_Te + energy_diff/Cemaxa
                Else
                  end_Te = Sqrt(start_Te*start_Te+2.0_wp*energy_diff*Tfermi/Cemaxa)
                  If (end_Te>Tfermi) end_Te = 0.5_wp*(start_Te*start_Te/Tfermi+Tfermi)+energy_diff/Cemaxa
                End If
                eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case Default
      ! tabulated volumetric heat capacity or more complex
      ! function: find new temperature iteratively by 
      ! gradual integration (0.01 kelvin at a time)
      ! and interpolate over last 0.01 kelvin
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2)*k)
              lat_I(i,j,k) = lat_I(i,j,k)+lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)
              energy_diff = lat_B(i,j,k)*act_ele_cell(ijk,0,0,0)*rvolume*eV_to_kB
              ! work out change in electronic temperature for energy deposition
              ! (searching first in 0.01 kelvin increments, then interpolate based
              ! on constant heat capacity)
              If (energy_diff > zero_plus) Then
                start_Te = eltemp(ijk,0,0,0)
                oldCe = Ce(start_Te)
                Do While (energy_diff > 0.0_wp)
                  newCe = Ce(start_Te+0.01_wp)
                  increase = 0.005_wp*(oldCe+newCe)
                  energy_diff = energy_diff - increase
                  start_Te = start_Te + 0.01_wp
                  oldCe = newCe
                End Do
                energy_diff = energy_diff + increase
                start_Te = start_Te - 0.01_wp
                newCe = Ce(start_Te)
                end_Te = start_Te + 2.0_wp * energy_diff / (oldCe+newCe)
                eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      End Select

    Else

    ! if at end of deposition, find deposited energy at end of deposition

      lat_I_min = Minval(lat_I(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      lat_I_max = Maxval(lat_I(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      lat_I_sum = Sum(lat_I(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      If (mxnode>1) Then
        Call gmin(lat_I_min)
        Call gmax(lat_I_max)
        Call gsum(lat_I_sum)
      End If

    ! find energy input into electronic temperature system

      lat_U_min = Minval(lat_U(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      lat_U_max = Maxval(lat_U(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      lat_U_sum = Sum(lat_U(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      If (mxnode>1) Then
        Call gmin(lat_U_min)
        Call gmax(lat_U_max)
        Call gsum(lat_U_sum)
      End If

    ! check how closely two values match up: if error greater than
    ! tolerance, report discrepancy as warning

      If (idnode == 0) Then
        If (Abs(lat_I_sum-lat_U_sum) > Abs(err_tol*lat_U_sum) .or. &
            Abs(lat_I_max-lat_U_max) > Abs(err_tol*lat_U_max) .or. &
            Abs(lat_I_min-lat_U_min) > Abs(err_tol*lat_U_min)) Then
          Call warning(530,Abs(lat_I_sum-lat_U_sum)/lat_U_sum*100.0_wp,0.0_wp,0.0_wp)
        End If
      End If

    ! report successful completion of energy deposition

      If (idnode == 0) Then
        If (currenttime<1.0_wp) Then
          Write(number, '(f14.3)') currenttime*1000.0_wp
          Write(nrite,"(/,6x,a,es12.5,a,a,a,/)") &
            'electronic energy deposition of ',lat_I_sum,' eV completed successfully after ',Trim(Adjustl(number)),' fs'
        Else
          Write(number, '(f14.6)') currenttime
          Write(nrite,"(/,6x,a,es12.5,a,a,a,/)") &
            'electronic energy deposition of ',lat_I_sum,' eV completed successfully after ',Trim(Adjustl(number)),' ps'
        End If
        Write(nrite,"(1x,130('-'))")
      End If

    ! switch off tracking and deallocate arrays

      trackInit = .false.
      findepo = .true.

      Deallocate(lat_U, Stat = fail(1))
      Deallocate(lat_B, Stat = fail(2))
      Deallocate(lat_I, Stat = fail(3))

      If (Any(fail > 0)) Call error(1090)

    End If
	
  End Subroutine depoevolve
	

  Subroutine uniformDist(lat_in)

! implement constant (homogeneous) spatial deposition
	
    Implicit None
    Real( Kind = wp ), Intent ( Inout ), Dimension(0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1) :: lat_in
    Real( Kind = wp ) :: dEdV

    ! express deposition energy per unit volume (eV/A^3):
    ! note penetration depth will be non-zero if laser is
    ! in use, otherwise use dE/dX value

    If (pdepth>zero_plus) Then
      dEdV = fluence/pdepth
    Else
      dEdV = dEdX/(Real(ntsys(1),Kind=wp)*Real(ntsys(2),Kind=wp)*delx*dely)
    End If
		
    ! homogeneous excitation: each temperature cell receives
    ! the same energy

    lat_in(1:ntcell(1),1:ntcell(2),1:ntcell(3))=dEdV*volume

  End Subroutine uniformDist

  Subroutine uniformDistZexp(lat_in)

! implement constant (homogeneous) spatial deposition
! in x and y-directions, exponential decay of fluence
! in z-direction (only with laser)
	
    Implicit None
    Real( Kind = wp ), Intent ( Inout ), Dimension(0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1) :: lat_in
    Real( Kind = wp ) :: dEdVmax, dEdV, zz, rpdepth
    Integer :: k

    ! express maximum deposition energy per unit volume (eV/A^3)

    If (pdepth>zero_plus) Then
      rpdepth = 1.0_wp/pdepth
    Else
      rpdepth = 0.0_wp
    End If
    dEdVmax = fluence*rpdepth

    ! loop through z-cells: calculate stopping power per
    ! cell based on z-position (maximum at z=0, grid centre)
    ! and assign across all x and y points in plane

    Do k = 1, ntcell(3)
      zz = Abs(Real(k+ntcelloff(3)-midI(3),Kind=wp))
      dEdV = dEdVmax*Exp(-zz*delz*rpdepth)
      lat_in(1:ntcell(1),1:ntcell(2),k) = dEdV*volume
    End Do

  End Subroutine uniformDistZexp

  Subroutine gaussianTrack(lat_in)

! implement gaussian spatial deposition

    Implicit None
    Real ( Kind = wp ), Intent ( Inout ), Dimension(0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1) :: lat_in
    Real ( Kind = wp ) :: normdEdX,realdEdx,sigmamx,sigmamy,sig2x,sig2y,sigcellx,sigcelly
    Real ( Kind = wp ) :: ii,jj,ii2,jj2,iip2,jjp2,iim2,jjm2
    Integer :: i,j,sgmx,sgmy
    Logical :: cutwarn=.false.

    lat_in(:,:,:) = 0.0_wp

    ! converting stopping power to a value per cell (in z-direction)

    normdEdX = dEdX*delz

    ! find extents of gaussian in x and y directions

    sigcellx = sig/delx
    sigcelly = sig/dely
    sig2x = 2.0_wp*sigcellx*sigcellx
    sig2y = 2.0_wp*sigcelly*sigcelly
    sigmamx = sigmax*sigcellx
    sigmamy = sigmax*sigcelly
    
    ! if cutoff larger than ionic temperature grid,
    ! warn of deposition errors

    If (sigmamx > Ceiling(ntsys(1)/2.0_wp)) Then
      sigmamx = Ceiling(ntsys(1)/2.0_wp)
      cutwarn = .true.
    End If
    If (sigmamy > Ceiling(ntsys(2)/2.0_wp)) Then
      sigmamy = Ceiling(ntsys(2)/2.0_wp)
      cutwarn = .true.
    End If

    If (idnode == 0 .and. cutwarn) Then
      Call warning(535,0.0_wp,0.0_wp,0.0_wp)
    End If

    sgmx = Nint(sigmamx)
    sgmy = Nint(sigmamy)

    ! apply five-point linear stencil for gaussian track:
    ! stencil modified to (hopefully!) deposit correct overall energy

    Do j=1,ntcell(2)
      jj = Real(j+ntcelloff(2)-midI(2),Kind=wp)
      jj2 = -jj*jj/sig2y
      jjp2 = -(jj+0.5_wp)*(jj+0.5_wp)/sig2y
      jjm2 = -(jj-0.5_wp)*(jj-0.5_wp)/sig2y
      Do i=1,ntcell(1)
        ii = Real(i+ntcelloff(1)-midI(1),Kind=wp)
        ii2 = -ii*ii/sig2x
        iip2 = -(ii+0.5_wp)*(ii+0.5_wp)/sig2x
        iim2 = -(ii-0.5_wp)*(ii-0.5_wp)/sig2x
        If (Abs(ii)<=sgmx .and. Abs(jj)<=sgmy) Then
          lat_in(i,j,1:ntcell(3)) = 0.2_wp*normdEdX/(2.0_wp*pi*sigcellx*sigcelly)*&
                                    (Exp(ii2+jj2)+Exp(iim2+jj2)+Exp(iim2+jj2)+Exp(ii2+jjp2)+Exp(ii2+jjm2))
        End If
      End Do
    End Do

    ! calculate deposited energy for comparison with specified value
    ! (note that stopping power is in z-direction)

    realdEdx = Sum(lat_in(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
    If (mxnode>1) Call gsum(realdEdx)
    realdEdx = realdEdx/(Real(ntsys(3),Kind=wp)*delz)

    ! check if lattice sum equals the expected value

    If (idnode == 0 .and. Abs((realdEdx-dEdX)/dEdX) > 0.01_wp) Then
      Call warning(540,Abs(realdEdx-dEdX)/dEdX*100_wp,0.0_wp,0.0_wp)
    End If
    
  End Subroutine gaussianTrack

End Module ttm_track_module
