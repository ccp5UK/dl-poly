Module ttm_track

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

  Use kinds, Only : wp
  Use setup
  Use ttm
  Use ttm_utils
  Use comms, Only : comms_type,Grid4_tag,Grid3_tag
  Use configuration
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None


Contains

  Subroutine depoevolve(time,tstep,redtstep,redtstepmx,comm)

! determine how deposition evolves over time

    Real( Kind = wp ), Intent ( In ) :: tstep,time
    Integer, Intent ( In ) :: redtstepmx,redtstep
    Type( comms_type), Intent( InOut ) :: comm

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
      If (comm%mxnode>1) Then
        Call gmin(comm,lat_I_min)
        Call gmax(comm,lat_I_max)
        Call gsum(comm,lat_I_sum)
      End If

    ! find energy input into electronic temperature system

      lat_U_min = Minval(lat_U(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      lat_U_max = Maxval(lat_U(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      lat_U_sum = Sum(lat_U(1:ntcell(1),1:ntcell(2),1:ntcell(3)))
      If (comm%mxnode>1) Then
        Call gmin(comm,lat_U_min)
        Call gmax(comm,lat_U_max)
        Call gsum(comm,lat_U_sum)
      End If

    ! check how closely two values match up: if error greater than
    ! tolerance, report discrepancy as warning

      If (comm%idnode == 0) Then
        If (Abs(lat_I_sum-lat_U_sum) > Abs(err_tol*lat_U_sum) .or. &
            Abs(lat_I_max-lat_U_max) > Abs(err_tol*lat_U_max) .or. &
            Abs(lat_I_min-lat_U_min) > Abs(err_tol*lat_U_min)) Then
          Call warning(530,Abs(lat_I_sum-lat_U_sum)/lat_U_sum*100.0_wp,0.0_wp,0.0_wp)
        End If
      End If

    ! report successful completion of energy deposition

      If (comm%idnode == 0) Then
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

  Subroutine gaussianTrack(lat_in, comm)

! implement gaussian spatial deposition

    Real ( Kind = wp ), Intent ( Inout ), Dimension(0:ntcell(1)+1,0:ntcell(2)+1,0:ntcell(3)+1) :: lat_in
    Type( comms_type), Intent( InOut) :: comm
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

    If (comm%idnode == 0 .and. cutwarn) Then
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
    If (comm%mxnode>1) Call gsum(comm,realdEdx)
    realdEdx = realdEdx/(Real(ntsys(3),Kind=wp)*delz)

    ! check if lattice sum equals the expected value

    If (comm%idnode == 0 .and. Abs((realdEdx-dEdX)/dEdX) > 0.01_wp) Then
      Call warning(540,Abs(realdEdx-dEdX)/dEdX*100_wp,0.0_wp,0.0_wp)
    End If
    
  End Subroutine gaussianTrack

  Subroutine ttm_ion_temperature(chi_ep,chi_es,vel_es2,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating ion temperatures in simulation
! cells for two-temperature model
!
! copyright - daresbury laboratory
! author    - s.l.darazewicz & m.a.seaton july 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real ( Kind = wp ), Intent ( In ) :: chi_ep,chi_es,vel_es2
  Type ( comms_type ), Intent( InOut) :: comm
  Integer :: ia,ja,ka,ijk,ijk1,ijk2,i,ii,jj,kk
  Real ( Kind = wp ) :: velsq,tmp,gsadd,vx,vy,vz,crho
  Integer :: fail, natmin
  Real ( Kind = wp ), Allocatable :: buf1(:),buf2(:),buf3(:),buf4(:)
  Integer, Allocatable :: nat(:),buf5(:),ijkatm(:)
  Integer, Dimension(8) :: req
  Integer, Dimension(MPI_STATUS_SIZE,8) :: stat

  ! allocate and zero arrays

  Allocate (nat (1:2*numcell), ijkatm (1:natms), Stat=fail)
  If (fail > 0) Call error(1085)

  nat = 0
  tempion = 0.0_wp
  asource = 0.0_wp
  gsource = 0.0_wp
  ttmvom = 0.0_wp
  ijkatm = 0

! zero variable for dynamic cell density calculations

  crho = 0.0_wp

! if heterogeneous, gsource is an array containing no. of atoms in each cell
! otherwise it is an array containing the effective electron-phonon relaxation 
! strength (known value for constant and homogeneous dynamic calculations)

  Select Case (gvar)
  Case (0,1)
    gsadd = chi_ep
  Case (2)
    gsadd = 1.0_wp
  End Select

! calculate overall momenta of cells for ionic temperature corrections

  Do i=1,natms

    ia = Floor((xxx(i)+zerocell(1))/delx) + 1
    ja = Floor((yyy(i)+zerocell(2))/dely) + 1
    ka = Floor((zzz(i)+zerocell(3))/delz) + 1

    ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
    ijkatm (i) = ijk

    tmp = weight(i)
    If (lfrzn(i) == 0) Then
      ttmvom(ijk,1) = ttmvom(ijk,1) + tmp*vxx(i)
      ttmvom(ijk,2) = ttmvom(ijk,2) + tmp*vyy(i)
      ttmvom(ijk,3) = ttmvom(ijk,3) + tmp*vzz(i)
      ttmvom(ijk,4) = ttmvom(ijk,4) + tmp
    End If

  End Do

  If (comm%mxnode>1) Then
    Allocate (buf1(1:numcell), buf2(1:numcell), buf3(1:numcell), buf4(1:numcell), Stat=fail)
    If (fail>0) Call error(1085)
    ! Sum up cell momenta and atomic masses in boundaries for ionic temperature corrections
    ! -z/+z directions
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf4 = 0.0_wp
    ijk1 = 1
    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2) * ntcell(3)
    Call MPI_ISEND (ttmvom(ijk1,1), 1, tmpmsgz, map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgz, map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,2), 1, tmpmsgz, map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgz, map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,3), 1, tmpmsgz, map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgz, map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,4), 1, tmpmsgz, map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, tmpmsgz, map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,1), 1, tmpmsgz, map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, tmpmsgz, map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,2), 1, tmpmsgz, map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, tmpmsgz, map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,3), 1, tmpmsgz, map(6), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, tmpmsgz, map(5), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,4), 1, tmpmsgz, map(6), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, tmpmsgz, map(5), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Do i=1,numcell
      ttmvom (i,1) = ttmvom (i,1) + buf1(i)
      ttmvom (i,2) = ttmvom (i,2) + buf2(i)
      ttmvom (i,3) = ttmvom (i,3) + buf3(i)
      ttmvom (i,4) = ttmvom (i,4) + buf4(i)
    End Do
    ! -y/+y directions
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf4 = 0.0_wp
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
    ijk2 = 1 + (ntcell(1)+2) * (2*ntcell(2)+2)
    Call MPI_ISEND (ttmvom(ijk1,1), 1, tmpmsgy, map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgy, map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,2), 1, tmpmsgy, map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgy, map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,3), 1, tmpmsgy, map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgy, map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,4), 1, tmpmsgy, map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, tmpmsgy, map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,1), 1, tmpmsgy, map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, tmpmsgy, map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,2), 1, tmpmsgy, map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, tmpmsgy, map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,3), 1, tmpmsgy, map(4), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, tmpmsgy, map(3), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,4), 1, tmpmsgy, map(4), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, tmpmsgy, map(3), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Do i=1,numcell
      ttmvom(i,1) = ttmvom(i,1) + buf1(i)
      ttmvom(i,2) = ttmvom(i,2) + buf2(i)
      ttmvom(i,3) = ttmvom(i,3) + buf3(i)
      ttmvom(i,4) = ttmvom(i,4) + buf4(i)
    End Do
    ! -x/+x directions
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf4 = 0.0_wp
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+3)
    ijk2 = 1 + ntcell(1) + (ntcell(1)+2) * (ntcell(2)+3)
    Call MPI_ISEND (ttmvom(ijk1,1), 1, tmpmsgx, map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgx, map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,2), 1, tmpmsgx, map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgx, map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,3), 1, tmpmsgx, map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgx, map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk1,4), 1, tmpmsgx, map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, tmpmsgx, map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,1), 1, tmpmsgx, map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, tmpmsgx, map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,2), 1, tmpmsgx, map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, tmpmsgx, map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,3), 1, tmpmsgx, map(2), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, tmpmsgx, map(1), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttmvom(ijk2,4), 1, tmpmsgx, map(2), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, tmpmsgx, map(1), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Do i=1,numcell
      ttmvom(i,1) = ttmvom(i,1) + buf1(i)
      ttmvom(i,2) = ttmvom(i,2) + buf2(i)
      ttmvom(i,3) = ttmvom(i,3) + buf3(i)
      ttmvom(i,4) = ttmvom(i,4) + buf4(i)
    End Do

    Deallocate (buf1, buf2, buf3, buf4, Stat=fail)
    If (fail>0) Call error(1086)
  End If

! calculate cell velocities

  Do i=1,numcell
    If (ttmvom(i,4)>zero_plus) Then
      ttmvom(i,1:3)=ttmvom(i,1:3)/ttmvom(i,4)
    Else
      ttmvom(i,1:3)=0.0_wp
    End If
  End Do

! calculate ionic temperatures (accounting for cell velocities)
! and source terms: electron-phonon (gsource) and electronic 
! stopping (asource)

  Do i=1,natms

    ijk = ijkatm(i)

    vx=vxx(i)-ttmvom(ijk,1)
    vy=vyy(i)-ttmvom(ijk,2)
    vz=vzz(i)-ttmvom(ijk,3)
    velsq = vx*vx+vy*vy+vz*vz
    tmp = weight(i)

    tempion(ijk) = tempion(ijk) + tmp*velsq

    gsource(ijk) = gsource(ijk) + gsadd

    nat(2*ijk-1) = nat(2*ijk-1) + 1

    If ((velsq > vel_es2) .and. (chi_es > zero_plus)) Then
      asource(ijk) = asource(ijk) + tmp*velsq
      nat(2*ijk) = nat(2*ijk) + 1
    End If

  End Do

  If (comm%mxnode>1) Then
    Allocate (buf1(1:numcell), buf2(1:numcell), buf3(1:numcell), buf5 (1:2*numcell), Stat=fail)
    If (fail>0) Call error(1085)
    ! Sum up boundary values of nat, tempion, asource and gsource
    ! -z direction (+z direction not needed)
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf5 = 0
    ijk1 = 1
    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2) * (ntcell(3))
    Call MPI_ISEND (tempion(ijk1) , 1, tmpmsgz, map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgz, map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (gsource(ijk1) , 1, tmpmsgz, map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgz, map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (asource(ijk1) , 1, tmpmsgz, map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgz, map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, nummsgz, map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, nummsgz, map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    tempion = tempion + buf1
    gsource = gsource + buf2
    asource = asource + buf3
    nat = nat + buf5
    ! -y direction  (+y direction not needed)
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf5 = 0
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
    ijk2 = 1 + (ntcell(1)+2) * (2*ntcell(2) + 2)
    Call MPI_ISEND (tempion(ijk1) , 1, tmpmsgy, map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgy, map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (gsource(ijk1) , 1, tmpmsgy, map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgy, map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (asource(ijk1) , 1, tmpmsgy, map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgy, map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, nummsgy, map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, nummsgy, map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    tempion = tempion + buf1
    gsource = gsource + buf2
    asource = asource + buf3
    nat = nat + buf5
    ! -x direction  (+x direction not needed)
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf5 = 0
    ijk1 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + ntcell(1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    Call MPI_ISEND (tempion(ijk1) , 1, tmpmsgx, map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgx, map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (gsource(ijk1) , 1, tmpmsgx, map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgx, map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (asource(ijk1) , 1, tmpmsgx, map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgx, map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, nummsgx, map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, nummsgx, map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    tempion = tempion + buf1
    gsource = gsource + buf2
    asource = asource + buf3
    nat = nat + buf5

    Deallocate (buf1, buf2, buf3, buf5, Stat=fail)
    If (fail>0) Call error(1086)
  End If

  natmin = Merge (0, amin-1, trackInit)
  old_ele_cell = act_ele_cell
  act_ele_cell = 1.0_wp
  acell_old = acell
  acell = 0

! loop through ionic temperature cells in current node
  Do ka = 1, ntcell(3)
    Do ja = 1, ntcell(2)
      Do ia = 1, ntcell(1)
        ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
        ! calculate ionic temperature for all cells with at least
        ! minimum number of particles (1 during deposition, amin 
        ! at all other times), calculating dynamic cell density
        ! (if required) from active cells, removing centre of mass
        ! motion and determining any inactive ionic temperature cells
        If (nat(2*ijk-1)>natmin .and. nat(2*ijk-1)>1) Then
          tempion(ijk) = tempion(ijk)/(3.0_wp*boltz*Real(nat(2*ijk-1),Kind=wp))
          acell = acell + 1
          crho = crho + gsource(ijk)
        Else If (natmin==0 .and. nat(2*ijk-1)==1) Then
          vx = 0.5_wp*ttmvom(ijk,1)
          vy = 0.5_wp*ttmvom(ijk,2)
          vz = 0.5_wp*ttmvom(ijk,3)
          velsq = ttmvom(ijk,4)*(vx*vx+vy*vy+vz*vz)
          tempion(ijk) = velsq/(3.0_wp*boltz)
          acell = acell + 1
          crho = crho + gsource(ijk)
        Else
          tempion(ijk) = 0.0_wp
          act_ele_cell(ijk,0,0,0) = 0.0_wp
        End If
        ! calculate electronic stopping terms (if more than one atom with speed > vel_cs)
        If (nat(2*ijk)>0) Then
          asource(ijk) = asource(ijk)*chi_es/boltz
        Else
          asource(ijk) = 0.0_wp
        End If
        ! calculate electron-phonon coupling terms
        gsource(ijk) = gsource(ijk)*3.0_wp
      End Do
    End Do
  End Do

! communicate and work out active ionic temperature cells in boundary halos

  If (comm%mxnode>1) Then
    Call gsum (acell)
    Call gsum (crho)
    ijk1 = 2 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ii = Merge (-1,0,(idx==nprx-1))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgx, map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,ii,0,0), 1, tmpmsgx, map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
    ijk1 = 1 + (ntcell(1)) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ii = Merge (1,0,(idx==0))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgx, map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,ii,0,0), 1, tmpmsgx, map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_WAITALL (4, req, stat, comm%ierr)
    ijk1 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2) + 1 + (ntcell(2)+2))
    jj = Merge (-1,0,(idy==npry-1))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgy, map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,jj,0), 1, tmpmsgy, map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
    jj = Merge (1,0,(idy==0))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgy, map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,jj,0), 1, tmpmsgy, map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_WAITALL (4, req, stat, comm%ierr)
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
    ijk2 = 1 + (ntcell(1)+2) * ((ntcell(2)+2) * (ntcell(3) + 1))
    kk = Merge (-1,0,(idz==nprz-1))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgz, map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,0,kk), 1, tmpmsgz, map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
    ijk1 = 1 + (ntcell(1)+2) * ((ntcell(2)+2) * ntcell(3))
    ijk2 = 1
    kk = Merge (1,0,(idz==0))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgz, map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,0,kk), 1, tmpmsgz, map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_WAITALL (4, req, stat, comm%ierr)
  Else
    Do ka = 1, ntcell(3)
      Do ja = 1, ntcell(2)
        ijk1 = 2 + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
        ijk2 = 1 + (ntcell(1)+1) + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
        act_ele_cell(ijk2,-1,0,0) = act_ele_cell(ijk1,0,0,0)
        ijk1 = 1 + ntcell(1) + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
        ijk2 = 1 + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
        act_ele_cell(ijk2,1,0,0) = act_ele_cell(ijk1,0,0,0)
      End Do
    End Do
    Do ka = 1, ntcell(3)
      Do ia = 0, ntcell(1)+1
        ijk1 = 1 + ia + (ntcell(1)+2) * (1 + (ntcell(2)+2) * ka)
        ijk2 = 1 + ia + (ntcell(1)+2) * (ntcell(2)+1 + (ntcell(2)+2) * ka)
        act_ele_cell(ijk2,0,-1,0) = act_ele_cell(ijk1,0,0,0)
        ijk1 = 1 + ia + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2) * ka)
        ijk2 = 1 + ia + (ntcell(1)+2) * (ntcell(2)+2) * ka
        act_ele_cell(ijk2,0,1,0) = act_ele_cell(ijk1,0,0,0)
      End Do
    End Do
    Do ja = 0, ntcell(2)+1
      Do ia = 0, ntcell(1)+1
        ijk1 = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2))
        ijk2 = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * (ntcell(3)+1))
        act_ele_cell(ijk2,0,0,-1) = act_ele_cell(ijk1,0,0,0)
        ijk1 = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ntcell(3))
        ijk2 = 1 + ia + (ntcell(1)+2) * ja
        act_ele_cell(ijk2,0,0,1) = act_ele_cell(ijk1,0,0,0)
      End Do
    End Do
  End If

! cell velocities generally used to correct velocities
! used in inhomogeneous Langevin thermostat: zero these
! velocities (x- and y-components only) if the user
! says otherwise (only for z-component)

  If (.not. ttmthvel) Then
    ttmvom = 0.0_wp
  Else If (ttmthvelz) Then
    ttmvom(1:numcell,1:2) = 0.0_wp
  End If

! dynamically calculate cell density for active cells
! if requested by user (after deposition stage)

  If (ttmdyndens .and. findepo) Then
    Select Case (gvar)
    Case (0,1)
      cellrho = crho / (Real(acell,Kind=wp)*chi_ep*volume)
    Case (2)
      cellrho = crho / (Real(acell,Kind=wp)*volume)
    End Select
    If (cellrho>zero_plus) Then
      rcellrho = 1.0_wp/cellrho
    Else
      rcellrho = 0.0_wp
    End If
  End If

! optional: send asource terms back to boundary voxels (only needed in +x, +y, +z directions)

!  If (mxnode>1) Then
!    ijk1 = 1 + ntcell(1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
!    ijk2 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
!    Call MPI_ISEND (asource(ijk1), 1, tmpmsgx, map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
!    Call MPI_IRECV (asource(ijk2), 1, tmpmsgx, map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
!    Call MPI_WAITALL (2, req, stat, comm%ierr)
!    ijk1 = 1 + (ntcell(1)+2) * (2*ntcell(2) + 2)
!    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
!    Call MPI_ISEND (asource(ijk1), 1, tmpmsgy, map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
!    Call MPI_IRECV (asource(ijk2), 1, tmpmsgy, map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
!    Call MPI_WAITALL (2, req, stat, comm%ierr)
!    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2) * (ntcell(3))
!    ijk2 = 1
!    Call MPI_ISEND (asource(ijk1), 1, tmpmsgz, map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
!    Call MPI_IRECV (asource(ijk2), 1, tmpmsgz, map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
!    Call MPI_WAITALL (2, req, stat, comm%ierr)
!  End If

  Deallocate (nat, ijkatm, Stat=fail)
  If (fail > 0) Call error(1086)

End Subroutine ttm_ion_temperature

Subroutine ttm_thermal_diffusion (tstep,time,nstep,nsteql,temp,nstbpo,ndump,nstrun,lines,npage,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for iterating thermal diffusion equation for
! electrons using numerical integration, after a.rutherford and
! d.duffy
!
! copyright - daresbury laboratory
! author    - s.l.darazewicz & m.a.seaton july 2012
! contrib -   g.khara may 2016
! contrib -   g.khara, s.t.murphy & m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer, Intent( In ) :: ndump,nstbpo,nsteql,nstep,nstrun,lines,npage
  Real ( Kind = wp ), Intent( In ) :: temp,tstep,time
  Type( comms_type), Intent( InOut ) :: comm

  Real ( Kind = wp ), Allocatable :: eltemp1(:,:,:,:)
  Real ( Kind = wp ) :: fomAx,fomAy,fomAz,mintstep,maxtstep,opttstep,delx2,dely2,delz2
  Real ( Kind = wp ) :: fopttstep,del2av,eltempmax,eltempmin,eltempmean,eltempKe,eltempmaxKe,eltempminKe
  Real ( Kind = wp ) :: actsite, actxm, actxp, actym, actyp, actzm, actzp, alploc
  Integer :: i,j,k,ii,jj,kk,ijk
  Logical :: safe
  Integer :: fail,redtstepmx,redtstep

! Debugging flag

  Logical :: debug1=.false.

! Initialise eltemp1 (electronic temperature grid for next timestep) and timestep sizes

  Allocate (eltemp1(1:numcell,-eltcell(1):eltcell(1),-eltcell(2):eltcell(2),-eltcell(3):eltcell(3)), Stat = fail)
  If (fail>0) Call error(1087)
  eltemp1 = 0.0_wp
  redtstepmx = 1

! deposition stage 1 (initialization):
! nstep-nsteql offsets equilibration time

  If ((nstep-nsteql)==1 .and. (sdepoType>0 .and. (dEdX>zero_plus .or. fluence>zero_plus))) Call depoinit (time,comm)

! determine timestep reduction factor (chosen empirically, acts beyond minimum stability condition)

  fopttstep = 0.25_wp

! determine maximum/minimum spacings and electronic temperatures

  delx2 = delx*delx
  dely2 = dely*dely
  delz2 = delz*delz
  del2av = (delx2*dely2*delz2)/(dely2*delz2+delx2*delz2+delx2*dely2)
  Call eltemp_max (eltempmax,comm)
  Call eltemp_min (eltempmin,comm)

! This section of the code establishes the optimum size of fourier mesh to ensure stability of the electronic
! temperature finite difference solver

  Select Case (KeType)
  Case (0)
! infinite thermal conductivity
    redtstepmx = 1
    opttstep = tstep
  Case (1)
! constant thermal conductivity and non-metal systems
    mintstep = 0.5_wp*del2av/alp(eltempmax)
    maxtstep = 0.5_wp*del2av/alp(eltempmin)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  Case (2)
! Drude-type thermal conductivity
    mintstep = 0.5_wp*del2av*Ce(eltempmax)/KeD(eltempmax,temp)
    maxtstep = 0.5_wp*del2av*Ce(eltempmin)/KeD(eltempmin,temp)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  Case (3)
  Call eltemp_maxKe (temp, eltempmaxKe,comm)
  Call eltemp_minKe (temp, eltempminKe,comm)
! tabulated thermal conductivity
    mintstep = 0.5_wp*del2av*Ce(eltempmax)/Ke(eltempmaxKe)
    maxtstep = 0.5_wp*del2av*Ce(eltempmin)/Ke(eltempminKe)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  End Select

! reduce timestep further for deposition stage

  If (KeType>0 .and. trackInit) Then
    Select Case (tdepoType)
    Case (1)
      redtstepmx = Max(50, redtstepmx)
    Case (2)
      redtstepmx = Max(10000, redtstepmx)
    Case Default
      redtstepmx = Max(2, redtstepmx)
    End Select
  End If

  fomAx = tstep/(delx2*Real(redtstepmx,Kind=wp))
  fomAy = tstep/(dely2*Real(redtstepmx,Kind=wp))
  fomAz = tstep/(delz2*Real(redtstepmx,Kind=wp))

! write information to OUTPUT

  If (mod(nstep,nstbpo) == 0 .or. nstep == 1) Then
    If (comm%idnode == 0) Then
      Write(nrite,'(6x,"ttm thermal diffusion timesteps:",2x,"optimal/ps",3x,"actual/ps",5x,"diff/md")')
      Write(nrite,'(38x,es12.4,es12.4,2x,i10)') opttstep, tstep/Real(redtstepmx,Kind=wp), redtstepmx
      If (ttmdyndens) Then
        Write(nrite,'(6x,"active ion temperature cells:",5x,"atom dens.",5x,"no. of active cells")')
        Write(nrite,'(38x,es12.4,14x,i10)') cellrho,acell
      End If
      If(nstep>1 .and. Mod(lines,npage)/=0) Write(nrite,"(1x,130('-'))")
    End If
  End If

! apply boundary conditions

  Call boundaryHalo (comm)
  Call boundaryCond (bcTypeE, temp,comm)

! print statistics to files: electronic and ionic temperatures
! (note timestep is subtracted by 1, as these are values at
!  beginning of MD timestep)

  Call printElecLatticeStatsToFile('PEAK_E', time, temp, nstep-1, ttmstats,comm)
  Call peakProfilerElec('LATS_E', nstep-1, ttmtraj,comm)

  Call printLatticeStatsToFile(tempion, 'PEAK_I', time, nstep-1, ttmstats,comm)
  Call peakProfiler(tempion, 'LATS_I', nstep-1, ttmtraj,comm)

! debugging option: print electron-phonon and electronic stopping source terms
!                   (normally switched off)

  If (debug1) Then
    Call printLatticeStatsToFile(gsource, 'PEAK_G', time, nstep-1, ttmstats,comm)
    Call peakProfiler(gsource, 'LATS_G', nstep-1, ttmtraj,comm)
    Call printLatticeStatsToFile(asource, 'PEAK_A', time, nstep-1, ttmstats,comm)
    Call peakProfiler(asource, 'LATS_A', nstep-1, ttmtraj,comm)
  End If

  safe=.true.

! determine energy redistribution from deactivated ionic temperature voxels for slab geometry

  If (redistribute) Call redistribute_Te (temp,comm)

! Adaptive timestep

  Do redtstep = 1, redtstepmx

! deposition stage 2 (with boundary conditions)

    If (trackInit) Then
      Call depoevolve(time, tstep, redtstep, redtstepmx,comm)
      Call boundaryCond(bcTypeE, temp,comm)
    End If

! MAIN LOOP
! this portion of the code is the main electronic temperature solver, solving
! for heat diffusion with two source terms within ionic cells, i.e.
!
!   Ce*d(T_e)/dt = d/dx(Ke * d(T_e)/dx) + volume*(asource + gsource*(T_e-T_i))
!
! the partial differential equation is solved using an explicit finite
! difference solver: care is needed in choosing timestep to ensure
! numerical stability

    Select Case (KeType)
    Case (0)
! infinite thermal conductivity case: set all electronic temperatures
! to mean value in active cells, to system temperature in inactive cells
      Call eltemp_mean(eltempmean,comm)
      eltemp1 = eltempmean
      Do ijk=1,numcell
        If(act_ele_cell(ijk,0,0,0)<=zero_plus) eltemp1(ijk,0,0,0) = temp
      End Do

    Case (1)
! constant thermal conductivity or non-metal case 
      If (redistribute) Then
      ! system with cell deactivation/energy redistribution
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1,numcell
                  If (adjust (ijk,ii,jj,kk)) eltemp(ijk,ii,jj,kk) = eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      actsite = act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*act_ele_cell(ijk-1,ii,jj,kk)
                      actxp = actsite*act_ele_cell(ijk+1,ii,jj,kk)
                      actym = actsite*act_ele_cell(ijk-(ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*act_ele_cell(ijk+(ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*act_ele_cell(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*act_ele_cell(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      alploc = alp(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      alploc = alp(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do
      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)
              Do k=1,ntcell(3)
                Do j=1,ntcell(2)
                  Do i=1,ntcell(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    alploc = alp(eltemp(ijk,ii,jj,kk))
                    eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                      fomAx*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAx*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    Case (2)
! Drude-type thermal conductivity case
      If (redistribute) Then
      ! system with cell deactivation/energy redistribution
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1, numcell
                  If (adjust (ijk,ii,jj,kk)) eltemp(ijk,ii,jj,kk) = eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      actsite = act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*act_ele_cell (ijk-1,ii,jj,kk)
                      actxp = actsite*act_ele_cell (ijk+1,ii,jj,kk)
                      actym = actsite*act_ele_cell (ijk-(ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*act_ele_cell (ijk+(ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*act_ele_cell (ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*act_ele_cell (ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-1,ii,jj,kk)),temp)*&
                        (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+1,ii,jj,kk)),temp)*&
                        (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-1,ii,jj,kk)),temp)*&
                        (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk)) +&
                        fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+1,ii,jj,kk)),temp)*&
                        (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do

      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)
              Do k=1,ntcell(3)
                Do j=1,ntcell(2)
                  Do i=1,ntcell(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                      fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-1,ii,jj,kk)),temp)*&
                      (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+1,ii,jj,kk)),temp)*&
                      (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    Case (3)
! tabulated thermal conductivity: uses local ionic or system temperature to calculate value
      If (redistribute) Then
      ! system with cell deactivation/energy redistribution
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1, numcell
                  If (adjust (ijk,ii,jj,kk)) eltemp(ijk,ii,jj,kk) = eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temperature sites (and active neighbours)
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      actsite = act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*act_ele_cell(ijk-1,ii,jj,kk)
                      actxp = actsite*act_ele_cell(ijk+1,ii,jj,kk)
                      actym = actsite*act_ele_cell(ijk-(ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*act_ele_cell(ijk+(ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*act_ele_cell(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*act_ele_cell(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      eltempKe = Merge(tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                      alploc = Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      ! note that temperature for thermal conductivity is always system
                      ! temperature for electronic cells away from ionic cells
                      alploc = Ke(temp)/Ce(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do
      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)
              Do k=1,ntcell(3)
                Do j=1,ntcell(2)
                  Do i=1,ntcell(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    eltempKe = Merge(tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                    eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                      fomAx*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAx*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    End Select

! electron stopping and electron-phonon couplings

    If (oneway) Then
      If (nstep > nstepcpl) Then
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
              If (act_ele_cell(ijk,0,0,0)>zero_plus) Then
                ! e-s coupling term
                eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)+tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*asource(ijk)
                ! e-p coupling term: only use if electronic temperature 
                ! exceeds ionic temperature
                If (l_epcp .and. eltemp(ijk,0,0,0)>tempion(ijk)) Then
                  Select Case (gvar)
                  Case (0,1)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))
                  Case (2)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))*&
                                                                                   Gep(eltemp(ijk,0,0,0))
                  End Select
                End If
              End If
            End Do
          End Do
        End Do
      End If
    Else
      If (nstep > nstepcpl) Then
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
              If (act_ele_cell(ijk,0,0,0)>zero_plus) Then
                ! e-s coupling term
                eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)+tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*asource(ijk)
                ! e-p coupling term
                If (l_epcp) Then
                  Select Case (gvar)
                  Case (0,1)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))
                  Case (2)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))*&
                                                                                   Gep(eltemp(ijk,0,0,0))
                  End Select
                End If
              End If
            End Do
          End Do
        End Do
      End If
    End If

! update electronic temperatures to adjusted values

    Do kk=-eltcell(3),eltcell(3)
      Do jj=-eltcell(2),eltcell(2)
        Do ii=-eltcell(1),eltcell(1)
          Do ijk=1,numcell
            eltemp(ijk,ii,jj,kk) = eltemp1(ijk,ii,jj,kk)
          End Do
        End Do
      End Do
    End Do

! update boundary halo values and apply boundary conditions

    Call boundaryHalo (comm)
    Call boundaryCond (bcTypeE, temp,comm)

! simple stability check for simulation

    If (Any(eltemp < 0.0_wp)) safe = .false.
    Call gcheck(safe)
    If (.not. safe) Call error (683)

  End Do

! Dumping Te file every ndump steps
  Call ttm_system_revive ('DUMP_E',nstep,time,ndump,nstrun,comm)

  Deallocate (eltemp1, Stat = fail)
  If (fail>0) Call error(1088)

End Subroutine ttm_thermal_diffusion

  
End Module ttm_track
