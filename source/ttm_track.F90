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
  Use ttm, Only : ttm_type,&
    eltemp_min, eltemp_max,eltemp_maxKe,eltemp_minKe,eltemp_mean,ttm_system_revive,&
    boundaryHalo,boundaryCond,depoinit
  Use ttm_utils, Only : Ce,Ke,alp,ked,peakProfilerElec,&
    peakProfiler,printElecLatticeStatsToFile,gep,printLatticeStatsToFile,redistribute_te

  Use comms, Only : comms_type,Grid4_tag,Grid3_tag,gsum, &
    grid1_tag,grid2_tag,gmin,gmax,gcheck
  Use configuration, Only : configuration_type
  Use errors_warnings, Only : error,warning,info
  Use thermostat, Only : thermostat_type
  Use domains, Only : domains_type
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None


Contains

  Subroutine depoevolve(time,tstep,redtstep,redtstepmx,ttm,comm)

! determine how deposition evolves over time
    Type( ttm_type ), Intent ( InOut ) :: ttm
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

    Character ( Len = 256 ) :: message

    ! provide atomic density corrections to heat capacities

    Ce0a   = ttm%Ce0  *Merge(ttm%cellrho,1.0_wp,ttm%ttmdyndens)
    sh_Aa  = ttm%sh_A *Merge(ttm%cellrho,1.0_wp,ttm%ttmdyndens)
    Cemaxa = ttm%Cemax*Merge(ttm%cellrho,1.0_wp,ttm%ttmdyndens)

    ! start deposition, reducing size of timestep for thermal diffusion

    currenttime = time-ttm%depostart+tstep/Real(redtstepmx,Kind=wp)*Real(redtstep,Kind=wp)

    Select Case (ttm%tdepoType)
    Case (1)
    ! Gaussian temporal deposition
      adjtime = currenttime/ttm%tdepo-ttm%tcdepo
      deposit = (currenttime<2.0_wp*ttm%tdepo*ttm%tcdepo)
      invbin = tstep/Real(redtstepmx,Kind=wp)
      adjeng = Exp(-0.5_wp*adjtime*adjtime)
    Case (2)
    ! decaying exponential temporal deposition
      adjtime = currenttime/ttm%tdepo
      deposit = (currenttime<ttm%tdepo*ttm%tcdepo)
      invbin = tstep/(ttm%tdepo*Real(redtstepmx,Kind=wp))
      adjeng = Exp(-adjtime)
    Case (3)
    ! delta temporal deposition (over single diffusion timestep)
      deposit = (currenttime<tstep/(Real(redtstepmx,Kind=wp)))
      invbin = 1.0_wp
      adjeng = 1.0_wp
    Case (4)
    ! pulse temporal deposition (over ttm%tdepo ps)
      deposit = (currenttime<ttm%tdepo)
      invbin = 1.0_wp
      adjeng = 1.0_wp
    End Select

    ! if (still) depositing energy, add to electronic temperature
    ! grid (active config%cells only) and ttm%adjust electronic temperatures 
    ! accordingly

    If (deposit) Then
      ttm%lat_B(:,:,:) = ttm%lat_U(:,:,:)*ttm%norm*invbin*adjeng
      Select Case (ttm%CeType)
      Case (0,4)
      ! constant specific heat capacity
        Do k=1,ttm%ntcell(3)
          Do j=1,ttm%ntcell(2)
            Do i=1,ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2)*k)
              ttm%lat_I(i,j,k) = ttm%lat_I(i,j,k)+ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)
              energy_diff = ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)*ttm%rvolume*ttm%eV_to_kB
              If (energy_diff>zero_plus) Then
                start_Te = ttm%eltemp(ijk,0,0,0)
                end_Te = start_Te + energy_diff/Ce0a
                ttm%eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (1,5)
      ! hyperbolic tangent specific heat capacity
        Do k=1,ttm%ntcell(3)
          Do j=1,ttm%ntcell(2)
            Do i=1,ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2)*k)
              ttm%lat_I(i,j,k) = ttm%lat_I(i,j,k)+ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)
              energy_diff = ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)*ttm%rvolume*ttm%eV_to_kB
              If (energy_diff>zero_plus) Then
                start_Te = ttm%eltemp(ijk,0,0,0)
                increase = Cosh(ttm%sh_B*start_Te)*Exp(ttm%sh_B*energy_diff/sh_Aa)
                ! using equivalent function: Acosh(x)=Log(x+Sqrt((x-1.0)*(x+1.0)))
                end_Te = Log(increase+Sqrt((increase-1.0_wp)*(increase+1.0_wp)))/ttm%sh_B
                ttm%eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case (2,6)
      ! linear specific heat capacity to Fermi temperature
        Do k=1,ttm%ntcell(3)
          Do j=1,ttm%ntcell(2)
            Do i=1,ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2)*k)
              ttm%lat_I(i,j,k) = ttm%lat_I(i,j,k)+ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)
              energy_diff = ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)*ttm%rvolume*ttm%eV_to_kB
              If (energy_diff>zero_plus) Then
                start_Te = ttm%eltemp(ijk,0,0,0)
                If (start_Te>=ttm%Tfermi) Then
                  end_Te = start_Te + energy_diff/Cemaxa
                Else
                  end_Te = Sqrt(start_Te*start_Te+2.0_wp*energy_diff*ttm%Tfermi/Cemaxa)
                  If (end_Te>ttm%Tfermi) end_Te = 0.5_wp*(start_Te*start_Te/ttm%Tfermi+ttm%Tfermi)+energy_diff/Cemaxa
                End If
                ttm%eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      Case Default
      ! tabulated ttm%volumetric heat capacity or more complex
      ! function: find new temperature iteratively by 
      ! gradual integration (0.01 kelvin at a time)
      ! and interpolate over last 0.01 kelvin
      Do k=1,ttm%ntcell(3)
          Do j=1,ttm%ntcell(2)
            Do i=1,ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2)*k)
              ttm%lat_I(i,j,k) = ttm%lat_I(i,j,k)+ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)
              energy_diff = ttm%lat_B(i,j,k)*ttm%act_ele_cell(ijk,0,0,0)*ttm%rvolume*ttm%eV_to_kB
              ! work out change in electronic temperature for energy deposition
              ! (searching first in 0.01 kelvin increments, then interpolate based
              ! on constant heat capacity)
              If (energy_diff > zero_plus) Then
                start_Te = ttm%eltemp(ijk,0,0,0)
                oldCe = Ce(start_Te,ttm)
                Do While (energy_diff > 0.0_wp)
                  newCe = Ce(start_Te+0.01_wp,ttm)
                  increase = 0.005_wp*(oldCe+newCe)
                  energy_diff = energy_diff - increase
                  start_Te = start_Te + 0.01_wp
                  oldCe = newCe
                End Do
                energy_diff = energy_diff + increase
                start_Te = start_Te - 0.01_wp
                newCe = Ce(start_Te,ttm)
                end_Te = start_Te + 2.0_wp * energy_diff / (oldCe+newCe)
                ttm%eltemp(ijk,0,0,0) = end_Te
              End If
            End Do
          End Do
        End Do
      End Select

    Else

    ! if at end of deposition, find deposited energy at end of deposition

      lat_I_min = Minval(ttm%lat_I(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
      lat_I_max = Maxval(ttm%lat_I(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
      lat_I_sum = Sum(ttm%lat_I(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
      Call gmin(comm,lat_I_min)
      Call gmax(comm,lat_I_max)
      Call gsum(comm,lat_I_sum)

    ! find energy input into electronic temperature system

      lat_U_min = Minval(ttm%lat_U(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
      lat_U_max = Maxval(ttm%lat_U(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
      lat_U_sum = Sum(ttm%lat_U(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
      Call gmin(comm,lat_U_min)
      Call gmax(comm,lat_U_max)
      Call gsum(comm,lat_U_sum)

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

      If (currenttime<1.0_wp) Then
        Write(number, '(f14.3)') currenttime*1000.0_wp
        Write(message,'(a,es12.5,3a)') &
          'electronic energy deposition of ',lat_I_sum, &
          ' eV completed successfully after ',Trim(Adjustl(number)),' fs'
      Else
        Write(number, '(f14.6)') currenttime
        Write(message,'(a,es12.5,3a)') &
          'electronic energy deposition of ',lat_I_sum, &
          ' eV completed successfully after ',Trim(Adjustl(number)),' ps'
      End If
      Call info(message,.true.)

    ! switch off tracking and deallocate arrays

      ttm%trackInit = .false.
      ttm%findepo = .true.

      Deallocate(ttm%lat_U, Stat = fail(1))
      Deallocate(ttm%lat_B, Stat = fail(2))
      Deallocate(ttm%lat_I, Stat = fail(3))

      If (Any(fail > 0)) Call error(1090)

    End If

  End Subroutine depoevolve

  Subroutine ttm_ion_temperature(ttm,thermo,domain,config,comm)

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
    Type( ttm_type ), Intent( InOut ) :: ttm
    Type( thermostat_type ), Intent( In    ) :: thermo
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut) :: comm
    Type( configuration_type ), Intent( InOut ) :: config

    Integer :: ia,ja,ka,ijk,ijk1,ijk2,i,ii,jj,kk
    Real ( Kind = wp ) :: velsq,tmp,gsadd,vx,vy,vz,crho
  Integer :: fail, natmin
  Real ( Kind = wp ), Allocatable :: buf1(:),buf2(:),buf3(:),buf4(:)
  Integer, Allocatable :: nat(:),buf5(:),ijkatm(:)
  Integer, Dimension(8) :: req
  Integer, Dimension(MPI_STATUS_SIZE,8) :: stat

  ! allocate and zero arrays

  Allocate (nat (1:2*ttm%numcell), ijkatm (1:config%natms), Stat=fail)
  If (fail > 0) Call error(1085)

  nat = 0
  ttm%tempion = 0.0_wp
  ttm%asource = 0.0_wp
  ttm%gsource = 0.0_wp
  ttm%ttmvom = 0.0_wp
  ijkatm = 0

! zero variable for dynamic cell density calculations

  crho = 0.0_wp

! if heterogeneous, ttm%gsource is an array containing no. of atoms in each cell
! otherwise it is an array containing the effective electron-phonon relaxation 
! strength (known value for constant and homogeneous dynamic calculations)

  Select Case (ttm%gvar)
  Case (0,1)
    gsadd = thermo%chi_ep
  Case (2)
    gsadd = 1.0_wp
  End Select

! calculate overall momenta of cells for ionic temperature corrections

  Do i=1,config%natms

    ia = Floor((config%parts(i)%xxx+ttm%zerocell(1))/ttm%delx) + 1
    ja = Floor((config%parts(i)%yyy+ttm%zerocell(2))/ttm%dely) + 1
    ka = Floor((config%parts(i)%zzz+ttm%zerocell(3))/ttm%delz) + 1

    ijk = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
    ijkatm (i) = ijk

    tmp = config%weight(i)
    If (config%lfrzn(i) == 0) Then
      ttm%ttmvom(ijk,1) = ttm%ttmvom(ijk,1) + tmp*config%vxx(i)
      ttm%ttmvom(ijk,2) = ttm%ttmvom(ijk,2) + tmp*config%vyy(i)
      ttm%ttmvom(ijk,3) = ttm%ttmvom(ijk,3) + tmp*config%vzz(i)
      ttm%ttmvom(ijk,4) = ttm%ttmvom(ijk,4) + tmp
    End If

  End Do

  If (comm%mxnode>1) Then
    Allocate (buf1(1:ttm%numcell), buf2(1:ttm%numcell), buf3(1:ttm%numcell), buf4(1:ttm%numcell), Stat=fail)
    If (fail>0) Call error(1085)
    ! Sum up config%cell momenta and atomic masses in boundaries for ionic temperature corrections
    ! -z/+z directions
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf4 = 0.0_wp
    ijk1 = 1
    ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2) * ttm%ntcell(3)
    Call MPI_ISEND (ttm%ttmvom(ijk1,1), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,2), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,3), 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,4), 1, ttm%tmpmsgz, domain%map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,2), 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,3), 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,4), 1, ttm%tmpmsgz, domain%map(6), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, ttm%tmpmsgz, domain%map(5), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Do i=1,ttm%numcell
      ttm%ttmvom (i,1) = ttm%ttmvom (i,1) + buf1(i)
      ttm%ttmvom (i,2) = ttm%ttmvom (i,2) + buf2(i)
      ttm%ttmvom (i,3) = ttm%ttmvom (i,3) + buf3(i)
      ttm%ttmvom (i,4) = ttm%ttmvom (i,4) + buf4(i)
    End Do
    ! -y/+y directions
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf4 = 0.0_wp
    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
    ijk2 = 1 + (ttm%ntcell(1)+2) * (2*ttm%ntcell(2)+2)
    Call MPI_ISEND (ttm%ttmvom(ijk1,1), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,2), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,3), 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,4), 1, ttm%tmpmsgy, domain%map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,1), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,2), 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,3), 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,4), 1, ttm%tmpmsgy, domain%map(4), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, ttm%tmpmsgy, domain%map(3), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Do i=1,ttm%numcell
      ttm%ttmvom(i,1) = ttm%ttmvom(i,1) + buf1(i)
      ttm%ttmvom(i,2) = ttm%ttmvom(i,2) + buf2(i)
      ttm%ttmvom(i,3) = ttm%ttmvom(i,3) + buf3(i)
      ttm%ttmvom(i,4) = ttm%ttmvom(i,4) + buf4(i)
    End Do
    ! -x/+x directions
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf4 = 0.0_wp
    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+3)
    ijk2 = 1 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+3)
    Call MPI_ISEND (ttm%ttmvom(ijk1,1), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,2), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,3), 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk1,4), 1, ttm%tmpmsgx, domain%map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,1), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,2), 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,3), 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (ttm%ttmvom(ijk2,4), 1, ttm%tmpmsgx, domain%map(2), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, ttm%tmpmsgx, domain%map(1), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    Do i=1,ttm%numcell
      ttm%ttmvom(i,1) = ttm%ttmvom(i,1) + buf1(i)
      ttm%ttmvom(i,2) = ttm%ttmvom(i,2) + buf2(i)
      ttm%ttmvom(i,3) = ttm%ttmvom(i,3) + buf3(i)
      ttm%ttmvom(i,4) = ttm%ttmvom(i,4) + buf4(i)
    End Do

    Deallocate (buf1, buf2, buf3, buf4, Stat=fail)
    If (fail>0) Call error(1086)
  End If

! calculate cell velocities

  Do i=1,ttm%numcell
    If (ttm%ttmvom(i,4)>zero_plus) Then
      ttm%ttmvom(i,1:3)=ttm%ttmvom(i,1:3)/ttm%ttmvom(i,4)
    Else
      ttm%ttmvom(i,1:3)=0.0_wp
    End If
  End Do

! calculate ionic temperatures (accounting for cell velocities)
! and source terms: electron-phonon (ttm%gsource) and electronic 
! stopping (ttm%asource)

  Do i=1,config%natms

    ijk = ijkatm(i)

    vx=config%vxx(i)-ttm%ttmvom(ijk,1)
    vy=config%vyy(i)-ttm%ttmvom(ijk,2)
    vz=config%vzz(i)-ttm%ttmvom(ijk,3)
    velsq = vx*vx+vy*vy+vz*vz
    tmp = config%weight(i)

    ttm%tempion(ijk) = ttm%tempion(ijk) + tmp*velsq

    ttm%gsource(ijk) = ttm%gsource(ijk) + gsadd

    nat(2*ijk-1) = nat(2*ijk-1) + 1

    If ((velsq > thermo%vel_es2) .and. (thermo%chi_es > zero_plus)) Then
      ttm%asource(ijk) = ttm%asource(ijk) + tmp*velsq
      nat(2*ijk) = nat(2*ijk) + 1
    End If

  End Do

  If (comm%mxnode>1) Then
    Allocate (buf1(1:ttm%numcell), buf2(1:ttm%numcell), buf3(1:ttm%numcell), buf5 (1:2*ttm%numcell), Stat=fail)
    If (fail>0) Call error(1085)
    ! Sum up boundary values of nat, ttm%tempion, ttm%asource and ttm%gsource
    ! -z direction (+z direction not needed)
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf5 = 0
    ijk1 = 1
    ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2) * (ttm%ntcell(3))
    Call MPI_ISEND (ttm%tempion(ijk1) , 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%gsource(ijk1) , 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%asource(ijk1) , 1, ttm%tmpmsgz, domain%map(5), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, ttm%tmpmsgz, domain%map(6), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, ttm%nummsgz, domain%map(5), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, ttm%nummsgz, domain%map(6), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    ttm%tempion = ttm%tempion + buf1
    ttm%gsource = ttm%gsource + buf2
    ttm%asource = ttm%asource + buf3
    nat = nat + buf5
    ! -y direction  (+y direction not needed)
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf5 = 0
    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
    ijk2 = 1 + (ttm%ntcell(1)+2) * (2*ttm%ntcell(2) + 2)
    Call MPI_ISEND (ttm%tempion(ijk1) , 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%gsource(ijk1) , 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%asource(ijk1) , 1, ttm%tmpmsgy, domain%map(3), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, ttm%tmpmsgy, domain%map(4), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, ttm%nummsgy, domain%map(3), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, ttm%nummsgy, domain%map(4), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    ttm%tempion = ttm%tempion + buf1
    ttm%gsource = ttm%gsource + buf2
    ttm%asource = ttm%asource + buf3
    nat = nat + buf5
    ! -x direction  (+x direction not needed)
    buf1 = 0.0_wp
    buf2 = 0.0_wp
    buf3 = 0.0_wp
    buf5 = 0
    ijk1 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    ijk2 = 1 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    Call MPI_ISEND (ttm%tempion(ijk1) , 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
    Call MPI_ISEND (ttm%gsource(ijk1) , 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_ISEND (ttm%asource(ijk1) , 1, ttm%tmpmsgx, domain%map(1), Grid3_tag, comm%comm, req(5), comm%ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, ttm%tmpmsgx, domain%map(2), Grid3_tag, comm%comm, req(6), comm%ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, ttm%nummsgx, domain%map(1), Grid4_tag, comm%comm, req(7), comm%ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, ttm%nummsgx, domain%map(2), Grid4_tag, comm%comm, req(8), comm%ierr)
    Call MPI_WAITALL (8, req, stat, comm%ierr)
    ttm%tempion = ttm%tempion + buf1
    ttm%gsource = ttm%gsource + buf2
    ttm%asource = ttm%asource + buf3
    nat = nat + buf5

    Deallocate (buf1, buf2, buf3, buf5, Stat=fail)
    If (fail>0) Call error(1086)
  End If

  natmin = Merge (0, ttm%amin-1, ttm%trackInit)
  ttm%old_ele_cell = ttm%act_ele_cell
  ttm%act_ele_cell = 1.0_wp
  ttm%acell_old = ttm%acell
  ttm%acell = 0

! loop through ionic temperature cells in current node
  Do ka = 1, ttm%ntcell(3)
    Do ja = 1, ttm%ntcell(2)
      Do ia = 1, ttm%ntcell(1)
        ijk = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
        ! calculate ionic temperature for all config%cells with at least
        ! minimum number of particles (1 during deposition, ttm%amin 
        ! at all other times), calculating dynamic config%cell density
        ! (if required) from active config%cells, removing centre of mass
        ! motion and determining any inactive ionic temperature config%cells
        If (nat(2*ijk-1)>natmin .and. nat(2*ijk-1)>1) Then
          ttm%tempion(ijk) = ttm%tempion(ijk)/(3.0_wp*boltz*Real(nat(2*ijk-1),Kind=wp))
          ttm%acell = ttm%acell + 1
          crho = crho + ttm%gsource(ijk)
        Else If (natmin==0 .and. nat(2*ijk-1)==1) Then
          vx = 0.5_wp*ttm%ttmvom(ijk,1)
          vy = 0.5_wp*ttm%ttmvom(ijk,2)
          vz = 0.5_wp*ttm%ttmvom(ijk,3)
          velsq = ttm%ttmvom(ijk,4)*(vx*vx+vy*vy+vz*vz)
          ttm%tempion(ijk) = velsq/(3.0_wp*boltz)
          ttm%acell = ttm%acell + 1
          crho = crho + ttm%gsource(ijk)
        Else
          ttm%tempion(ijk) = 0.0_wp
          ttm%act_ele_cell(ijk,0,0,0) = 0.0_wp
        End If
        ! calculate electronic stopping terms (if more than one atom with speed > vel_cs)
        If (nat(2*ijk)>0) Then
          ttm%asource(ijk) = ttm%asource(ijk)*thermo%chi_es/boltz
        Else
          ttm%asource(ijk) = 0.0_wp
        End If
        ! calculate electron-phonon coupling terms
        ttm%gsource(ijk) = ttm%gsource(ijk)*3.0_wp
      End Do
    End Do
  End Do

! communicate and work out active ionic temperature cells in boundary halos

  If (comm%mxnode>1) Then
    Call gsum (comm,ttm%acell)
    Call gsum (comm,crho)
    ijk1 = 2 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    ijk2 = 1 + (ttm%ntcell(1)+1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    ii = Merge (-1,0,(domain%idx==domain%nx-1))
    Call MPI_ISEND (ttm%act_ele_cell(ijk1,0,0,0),  1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (ttm%act_ele_cell(ijk2,ii,0,0), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(2), comm%ierr)
    ijk1 = 1 + (ttm%ntcell(1)) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    ijk2 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    ii = Merge (1,0,(domain%idx==0))
    Call MPI_ISEND (ttm%act_ele_cell(ijk1,0,0,0),  1, ttm%tmpmsgx, domain%map(2), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (ttm%act_ele_cell(ijk2,ii,0,0), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_WAITALL (4, req, stat, comm%ierr)
    ijk1 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
    ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2)+2))
    jj = Merge (-1,0,(domain%idy==domain%ny-1))
    Call MPI_ISEND (ttm%act_ele_cell(ijk1,0,0,0),  1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (ttm%act_ele_cell(ijk2,0,jj,0), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(2), comm%ierr)
    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + (ttm%ntcell(2)+2))
    ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
    jj = Merge (1,0,(domain%idy==0))
    Call MPI_ISEND (ttm%act_ele_cell(ijk1,0,0,0),  1, ttm%tmpmsgy, domain%map(4), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (ttm%act_ele_cell(ijk2,0,jj,0), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_WAITALL (4, req, stat, comm%ierr)
    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
    ijk2 = 1 + (ttm%ntcell(1)+2) * ((ttm%ntcell(2)+2) * (ttm%ntcell(3) + 1))
    kk = Merge (-1,0,(domain%idz==domain%nz-1))
    Call MPI_ISEND (ttm%act_ele_cell(ijk1,0,0,0),  1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(1), comm%ierr)
    Call MPI_IRECV (ttm%act_ele_cell(ijk2,0,0,kk), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(2), comm%ierr)
    ijk1 = 1 + (ttm%ntcell(1)+2) * ((ttm%ntcell(2)+2) * ttm%ntcell(3))
    ijk2 = 1
    kk = Merge (1,0,(domain%idz==0))
    Call MPI_ISEND (ttm%act_ele_cell(ijk1,0,0,0),  1, ttm%tmpmsgz, domain%map(6), Grid2_tag, comm%comm, req(3), comm%ierr)
    Call MPI_IRECV (ttm%act_ele_cell(ijk2,0,0,kk), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, comm%comm, req(4), comm%ierr)
    Call MPI_WAITALL (4, req, stat, comm%ierr)
  Else
    Do ka = 1, ttm%ntcell(3)
      Do ja = 1, ttm%ntcell(2)
        ijk1 = 2 + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
        ijk2 = 1 + (ttm%ntcell(1)+1) + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
        ttm%act_ele_cell(ijk2,-1,0,0) = ttm%act_ele_cell(ijk1,0,0,0)
        ijk1 = 1 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
        ijk2 = 1 + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ka)
        ttm%act_ele_cell(ijk2,1,0,0) = ttm%act_ele_cell(ijk1,0,0,0)
      End Do
    End Do
    Do ka = 1, ttm%ntcell(3)
      Do ia = 0, ttm%ntcell(1)+1
        ijk1 = 1 + ia + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2) * ka)
        ijk2 = 1 + ia + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+1 + (ttm%ntcell(2)+2) * ka)
        ttm%act_ele_cell(ijk2,0,-1,0) = ttm%act_ele_cell(ijk1,0,0,0)
        ijk1 = 1 + ia + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + (ttm%ntcell(2)+2) * ka)
        ijk2 = 1 + ia + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2) * ka
        ttm%act_ele_cell(ijk2,0,1,0) = ttm%act_ele_cell(ijk1,0,0,0)
      End Do
    End Do
    Do ja = 0, ttm%ntcell(2)+1
      Do ia = 0, ttm%ntcell(1)+1
        ijk1 = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2))
        ijk2 = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * (ttm%ntcell(3)+1))
        ttm%act_ele_cell(ijk2,0,0,-1) = ttm%act_ele_cell(ijk1,0,0,0)
        ijk1 = 1 + ia + (ttm%ntcell(1)+2) * (ja + (ttm%ntcell(2)+2) * ttm%ntcell(3))
        ijk2 = 1 + ia + (ttm%ntcell(1)+2) * ja
        ttm%act_ele_cell(ijk2,0,0,1) = ttm%act_ele_cell(ijk1,0,0,0)
      End Do
    End Do
  End If

! cell velocities generally used to correct velocities
! used in inhomogeneous Langevin thermostat: zero these
! velocities (x- and y-components only) if the user
! says otherwise (only for z-component)

  If (.not. ttm%ttmthvel) Then
    ttm%ttmvom = 0.0_wp
  Else If (ttm%ttmthvelz) Then
    ttm%ttmvom(1:ttm%numcell,1:2) = 0.0_wp
  End If

! dynamically calculate cell density for active cells
! if requested by user (after deposition stage)

  If (ttm%ttmdyndens .and. ttm%findepo) Then
    Select Case (ttm%gvar)
    Case (0,1)
      ttm%cellrho = crho / (Real(ttm%acell,Kind=wp)*thermo%chi_ep*ttm%volume)
    Case (2)
      ttm%cellrho = crho / (Real(ttm%acell,Kind=wp)*ttm%volume)
    End Select
    If (ttm%cellrho>zero_plus) Then
      ttm%rcellrho = 1.0_wp/ttm%cellrho
    Else
      ttm%rcellrho = 0.0_wp
    End If
  End If

! optional: send ttm%asource terms back to boundary voxels (only needed in +x, +y, +z directions)

!  If (mxnode>1) Then
!    ijk1 = 1 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
!    ijk2 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
!    Call MPI_ISEND (ttm%asource(ijk1), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, comm%comm, req(1), comm%ierr)
!    Call MPI_IRECV (ttm%asource(ijk2), 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, comm%comm, req(2), comm%ierr)
!    Call MPI_WAITALL (2, req, stat, comm%ierr)
!    ijk1 = 1 + (ttm%ntcell(1)+2) * (2*ttm%ntcell(2) + 2)
!    ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
!    Call MPI_ISEND (ttm%asource(ijk1), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, comm%comm, req(1), comm%ierr)
!    Call MPI_IRECV (ttm%asource(ijk2), 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, comm%comm, req(2), comm%ierr)
!    Call MPI_WAITALL (2, req, stat, comm%ierr)
!    ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2) * (ttm%ntcell(3))
!    ijk2 = 1
!    Call MPI_ISEND (ttm%asource(ijk1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, comm%comm, req(1), comm%ierr)
!    Call MPI_IRECV (ttm%asource(ijk2), 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, comm%comm, req(2), comm%ierr)
!    Call MPI_WAITALL (2, req, stat, comm%ierr)
!  End If

  Deallocate (nat, ijkatm, Stat=fail)
  If (fail > 0) Call error(1086)

End Subroutine ttm_ion_temperature

Subroutine ttm_thermal_diffusion(tstep,time,nstep,nsteql,nstbpo,ndump,nstrun, &
  ttm,thermo,domain,comm)

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

  Type(ttm_type), Intent( InOut ) :: ttm
  Integer, Intent( In ) :: ndump,nstbpo,nsteql,nstep,nstrun
  Real ( Kind = wp ), Intent( In ) :: tstep,time
  Type( thermostat_type ), Intent( In    ) :: thermo
  Type( domains_type ), Intent( In    ) :: domain
  Type( comms_type), Intent( InOut ) :: comm

  Real ( Kind = wp ), Allocatable :: eltemp1(:,:,:,:)
  Real ( Kind = wp ) :: fomAx,fomAy,fomAz,mintstep,maxtstep,opttstep,delx2,dely2,delz2
  Real ( Kind = wp ) :: fopttstep,del2av,eltempmax,eltempmin,eltempmean,eltempKe,eltempmaxKe,eltempminKe
  Real ( Kind = wp ) :: actsite, actxm, actxp, actym, actyp, actzm, actzp, alploc
  Real ( Kind = wp ) :: temp
  Integer :: i,j,k,ii,jj,kk,ijk
  Logical :: safe
  Integer :: fail,redtstepmx,redtstep

! Debugging flag

  Logical :: debug1=.false.

  Character( Len = 256 ) :: messages(6)

! Initialise eltemp1 (electronic temperature grid for next timestep) and timestep sizes

  Allocate (eltemp1(1:ttm%numcell,-ttm%eltcell(1):ttm%eltcell(1),-ttm%eltcell(2):ttm%eltcell(2),-ttm%eltcell(3):ttm%eltcell(3)), &
    Stat = fail)
  If (fail>0) Call error(1087)
  eltemp1 = 0.0_wp
  redtstepmx = 1

! Initialise temp
  temp = thermo%temp

! deposition stage 1 (initialization):
! nstep-nsteql offsets equilibration time

  If ((nstep-nsteql)==1 .and. (ttm%sdepoType>0 .and. (ttm%dEdX>zero_plus .or. ttm%fluence>zero_plus))) Call depoinit (time,ttm,comm)

! determine timestep reduction factor (chosen empirically, acts beyond minimum stability condition)

  fopttstep = 0.25_wp

! determine maximum/minimum spacings and electronic temperatures

  delx2 = ttm%delx*ttm%delx
  dely2 = ttm%dely*ttm%dely
  delz2 = ttm%delz*ttm%delz
  del2av = (delx2*dely2*delz2)/(dely2*delz2+delx2*delz2+delx2*dely2)
  Call eltemp_max (eltempmax,ttm,comm)
  Call eltemp_min (eltempmin,ttm,comm)

! This section of the code establishes the optimum size of fourier mesh to ensure stability of the electronic
! temperature finite difference solver

  Select Case (ttm%KeType)
  Case (0)
! infinite thermal conductivity
    redtstepmx = 1
    opttstep = tstep
   Case (1)
! constant thermal conductivity and non-metal systems
    mintstep = 0.5_wp*del2av/alp(eltempmax,ttm)
    maxtstep = 0.5_wp*del2av/alp(eltempmin,ttm)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
   Case (2)
! Drude-type thermal conductivity
    mintstep = 0.5_wp*del2av*Ce(eltempmax,ttm)/KeD(eltempmax,temp,ttm)
    maxtstep = 0.5_wp*del2av*Ce(eltempmin,ttm)/KeD(eltempmin,temp,ttm)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
   Case (3)
    Call eltemp_maxKe (temp, eltempmaxKe,ttm,comm)
    Call eltemp_minke (temp, eltempminKe,ttm,comm)
! tabulated thermal conductivity
  mintstep = 0.5_wp*del2av*Ce(eltempmax,ttm)/Ke(eltempmaxKe,ttm)
  maxtstep = 0.5_wp*del2av*Ce(eltempmin,ttm)/Ke(eltempminKe,ttm)
  opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  End Select

! reduce timestep further for deposition stage

  If (ttm%KeType>0 .and. ttm%trackInit) Then
    Select Case (ttm%tdepoType)
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
    Write(messages(1),'(a)') 'ttm thermal diffusion timesteps:'
    Write(messages(2),'(4x,a,3x,a,5x,a)') 'optimal/ps','actual/ps','diff/md'
    Write(messages(3),'(2x,2es12.4,2x,i10)') opttstep, tstep/Real(redtstepmx,Kind=wp),redtstepmx
    If (ttm%ttmdyndens) Then
      Write(messages(4),'(a)') 'active ion temperature config%cells:'
      Write(messages(5),'(4x,a,2x,a)') 'atom dens.','no. of active config%cells'
      Write(messages(6),'(2x,es12.4,11x,i10)') ttm%cellrho,ttm%acell
      Call info(messages,6,.true.)
    Else
      Call info(messages,3,.true.)
    End If
  End If

! apply boundary conditions

  Call boundaryHalo (ttm,domain,comm)
  Call boundaryCond (ttm%bcTypeE,temp,ttm,comm)

! print statistics to files: electronic and ionic temperatures
! (note timestep is subtracted by 1, as these are values at
!  beginning of MD timestep)

  Call printElecLatticeStatsToFile('PEAK_E', time, temp, nstep-1, ttm%ttmstats,ttm,comm)
  Call peakProfilerElec('LATS_E', nstep-1, ttm%ttmtraj,ttm,comm)

  Call printLatticeStatsToFile(ttm%tempion, 'PEAK_I', time, nstep-1, ttm%ttmstats,ttm,comm)
  Call peakProfiler(ttm%tempion, 'LATS_I', nstep-1, ttm%ttmtraj,ttm,comm)

! debugging option: print electron-phonon and electronic stopping source terms
!                   (ttm%normally switched off)

  If (debug1) Then
    Call printLatticeStatsToFile(ttm%gsource, 'PEAK_G', time, nstep-1, ttm%ttmstats,ttm,comm)
    Call peakProfiler(ttm%gsource, 'LATS_G', nstep-1, ttm%ttmtraj,ttm,comm)
    Call printLatticeStatsToFile(ttm%asource, 'PEAK_A', time, nstep-1, ttm%ttmstats,ttm,comm)
    Call peakProfiler(ttm%asource, 'LATS_A', nstep-1, ttm%ttmtraj,ttm,comm)
  End If

  safe=.true.

! determine energy redistribution from deactivated ionic temperature voxels for slab geometry

  If (ttm%redistribute) Then
    Call redistribute_te (temp,ttm,domain,comm)
  End If

! Adaptive timestep

  Do redtstep = 1, redtstepmx

! deposition stage 2 (with boundary conditions)

    If (ttm%trackInit) Then
      Call depoevolve(time, tstep, redtstep, redtstepmx,ttm,comm)
      Call boundaryCond(ttm%bcTypeE, temp,ttm,comm)
    End If

! MAIN LOOP
! this portion of the code is the main electronic temperature solver, solving
! for heat diffusion with two source terms within ionic cells, i.e.
!
!   Ce*d(T_e)/dt = d/dx(Ke * d(T_e)/dx) + ttm%volume*(ttm%asource + ttm%gsource*(T_e-T_i))
!
! the partial differential equation is solved using an explicit finite
! difference solver: care is needed in choosing timestep to ensure
! numerical stability

    Select Case (ttm%KeType)
    Case (0)
! infinite thermal conductivity case: set all electronic temperatures
! to mean value in active cells, to system temperature in inactive cells
      Call eltemp_mean(eltempmean,ttm,comm)
      eltemp1 = eltempmean
      Do ijk=1,ttm%numcell
        If(ttm%act_ele_cell(ijk,0,0,0)<=zero_plus) eltemp1(ijk,0,0,0) = temp
      End Do

    Case (1)
! constant thermal conductivity or non-metal case 
      If (ttm%redistribute) Then
      ! system with config%cell deactivation/energy redistribution
        Do kk=-ttm%eltcell(3),ttm%eltcell(3)
          Do jj=-ttm%eltcell(2),ttm%eltcell(2)
            Do ii=-ttm%eltcell(1),ttm%eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1,ttm%numcell
                  If (ttm%adjust (ijk,ii,jj,kk)) ttm%eltemp(ijk,ii,jj,kk) = ttm%eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                Do k=1,ttm%ntcell(3)
                  Do j=1,ttm%ntcell(2)
                    Do i=1,ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      actsite = ttm%act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*ttm%act_ele_cell(ijk-1,ii,jj,kk)
                      actxp = actsite*ttm%act_ele_cell(ijk+1,ii,jj,kk)
                      actym = actsite*ttm%act_ele_cell(ijk-(ttm%ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*ttm%act_ele_cell(ijk+(ttm%ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*ttm%act_ele_cell(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*ttm%act_ele_cell(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)
                      alploc = alp(ttm%eltemp(ijk,ii,jj,kk),ttm)
                      eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*alploc*(ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*alploc*(ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic config%cells away from ionic config%cells
                Do k=1,ttm%ntcell(3)
                  Do j=1,ttm%ntcell(2)
                    Do i=1,ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      alploc = alp(ttm%eltemp(ijk,ii,jj,kk),ttm)
                      eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                        fomAx*alploc*(ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAx*alploc*(ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do
      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-ttm%eltcell(3),ttm%eltcell(3)
          Do jj=-ttm%eltcell(2),ttm%eltcell(2)
            Do ii=-ttm%eltcell(1),ttm%eltcell(1)
              Do k=1,ttm%ntcell(3)
                Do j=1,ttm%ntcell(2)
                  Do i=1,ttm%ntcell(1)
                    ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                    alploc = alp(ttm%eltemp(ijk,ii,jj,kk),ttm)
                    eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                      fomAx*alploc*(ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAx*alploc*(ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAy*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAy*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAz*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAz*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    Case (2)
! Drude-type thermal conductivity case
      If (ttm%redistribute) Then
      ! system with config%cell deactivation/energy redistribution
        Do kk=-ttm%eltcell(3),ttm%eltcell(3)
          Do jj=-ttm%eltcell(2),ttm%eltcell(2)
            Do ii=-ttm%eltcell(1),ttm%eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1, ttm%numcell
                  If (ttm%adjust (ijk,ii,jj,kk)) ttm%eltemp(ijk,ii,jj,kk) = ttm%eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                Do k=1,ttm%ntcell(3)
                  Do j=1,ttm%ntcell(2)
                    Do i=1,ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      actsite = ttm%act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*ttm%act_ele_cell (ijk-1,ii,jj,kk)
                      actxp = actsite*ttm%act_ele_cell (ijk+1,ii,jj,kk)
                      actym = actsite*ttm%act_ele_cell (ijk-(ttm%ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*ttm%act_ele_cell (ijk+(ttm%ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*ttm%act_ele_cell (ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*ttm%act_ele_cell (ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)
                      eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-1,ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAx*actxp*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+1,ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAy*actym*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAy*actyp*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAz*actzm*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-(ttm%ntcell(1)+2)*&
                        (ttm%ntcell(2)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/&
                        Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAz*actzp*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+(ttm%ntcell(1)+2)*&
                        (ttm%ntcell(2)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/&
                        Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic config%cells away from ionic config%cells
                Do k=1,ttm%ntcell(3)
                  Do j=1,ttm%ntcell(2)
                    Do i=1,ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                        fomAx*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-1,ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm) +&
                        fomAx*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+1,ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAy*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAy*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAz*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-(ttm%ntcell(1)+2)*&
                        (ttm%ntcell(2)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/&
                        Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                        fomAz*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+(ttm%ntcell(1)+2)*&
                        (ttm%ntcell(2)+2),ii,jj,kk)),temp,ttm)*&
                        (ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/&
                        Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do

      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-ttm%eltcell(3),ttm%eltcell(3)
          Do jj=-ttm%eltcell(2),ttm%eltcell(2)
            Do ii=-ttm%eltcell(1),ttm%eltcell(1)
              Do k=1,ttm%ntcell(3)
                Do j=1,ttm%ntcell(2)
                  Do i=1,ttm%ntcell(1)
                    ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                    eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                      fomAx*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-1,ii,jj,kk)),temp,ttm)*&
                      (ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                      fomAx*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+1,ii,jj,kk)),temp,ttm)*&
                      (ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                      fomAy*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)),temp,ttm)*&
                      (ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                      fomAy*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)),temp,ttm)*&
                      (ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+&
                      fomAz*KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk-(ttm%ntcell(1)+2)*&
                      (ttm%ntcell(2)+2),ii,jj,kk)),temp,ttm)*&
                      (ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/&
                      Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)+fomAz*&
                      KeD(0.5_wp*(ttm%eltemp(ijk,ii,jj,kk)+ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)),temp,ttm)*&
                      (ttm%eltemp(ijk+(ttm%ntcell(1)+2)*&
                      (ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    Case (3)
! tabulated thermal conductivity: uses local ionic or system temperature to calculate value
      If (ttm%redistribute) Then
      ! system with config%cell deactivation/energy redistribution
        Do kk=-ttm%eltcell(3),ttm%eltcell(3)
          Do jj=-ttm%eltcell(2),ttm%eltcell(2)
            Do ii=-ttm%eltcell(1),ttm%eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1, ttm%numcell
                  If (ttm%adjust (ijk,ii,jj,kk)) ttm%eltemp(ijk,ii,jj,kk) = ttm%eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temperature sites (and active neighbours)
                Do k=1,ttm%ntcell(3)
                  Do j=1,ttm%ntcell(2)
                    Do i=1,ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      actsite = ttm%act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*ttm%act_ele_cell(ijk-1,ii,jj,kk)
                      actxp = actsite*ttm%act_ele_cell(ijk+1,ii,jj,kk)
                      actym = actsite*ttm%act_ele_cell(ijk-(ttm%ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*ttm%act_ele_cell(ijk+(ttm%ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*ttm%act_ele_cell(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*ttm%act_ele_cell(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)
                      eltempKe = Merge(ttm%tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                      alploc = Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)
                      eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*alploc*(ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*alploc*(ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic config%cells away from ionic config%cells
                Do k=1,ttm%ntcell(3)
                  Do j=1,ttm%ntcell(2)
                    Do i=1,ttm%ntcell(1)
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      ! note that temperature for thermal conductivity is always system
                      ! temperature for electronic config%cells away from ionic config%cells
                      alploc = Ke(temp,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)
                      eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                        fomAx*alploc*(ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAx*alploc*(ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do
      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-ttm%eltcell(3),ttm%eltcell(3)
          Do jj=-ttm%eltcell(2),ttm%eltcell(2)
            Do ii=-ttm%eltcell(1),ttm%eltcell(1)
              Do k=1,ttm%ntcell(3)
                Do j=1,ttm%ntcell(2)
                  Do i=1,ttm%ntcell(1)
                    ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                    eltempKe = Merge(ttm%tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                    eltemp1(ijk,ii,jj,kk) = ttm%eltemp(ijk,ii,jj,kk)+&
                      fomAx*Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)*&
                      (ttm%eltemp(ijk-1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAx*Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)*&
                      (ttm%eltemp(ijk+1,ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAy*Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)*&
                      (ttm%eltemp(ijk-(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAy*Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)*&
                      (ttm%eltemp(ijk+(ttm%ntcell(1)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAz*Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)*&
                      (ttm%eltemp(ijk-(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))+&
                      fomAz*Ke(eltempKe,ttm)/Ce(ttm%eltemp(ijk,ii,jj,kk),ttm)*&
                      (ttm%eltemp(ijk+(ttm%ntcell(1)+2)*(ttm%ntcell(2)+2),ii,jj,kk)-ttm%eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    End Select

! electron stopping and electron-phonon couplings

    If (ttm%oneway) Then
      If (nstep > ttm%nstepcpl) Then
        Do k=1,ttm%ntcell(3)
          Do j=1,ttm%ntcell(2)
            Do i=1,ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
              If (ttm%act_ele_cell(ijk,0,0,0)>zero_plus) Then
                ! e-s coupling term
                eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)+tstep*ttm%rvolume/(Ce(ttm%eltemp(ijk,0,0,0),ttm)*&
                  Real(redtstepmx,Kind=wp))*ttm%asource(ijk)
                ! e-p coupling term: only use if electronic temperature 
                ! exceeds ionic temperature
                If (ttm%l_epcp .and. ttm%eltemp(ijk,0,0,0)>ttm%tempion(ijk)) Then
                  Select Case (ttm%gvar)
                   Case (0,1)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                      tstep*ttm%rvolume/(Ce(ttm%eltemp(ijk,0,0,0),ttm)*Real(redtstepmx,Kind=wp))*&
                      ttm%gsource(ijk)*(ttm%eltemp(ijk,0,0,0)-ttm%tempion(ijk))
                   Case (2)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                      tstep*ttm%rvolume/(Ce(ttm%eltemp(ijk,0,0,0),ttm)*Real(redtstepmx,Kind=wp))*&
                      ttm%gsource(ijk)*(ttm%eltemp(ijk,0,0,0)-ttm%tempion(ijk))*&
                      Gep(ttm%eltemp(ijk,0,0,0),ttm)
                  End Select
                End If
              End If
            End Do
          End Do
        End Do
      End If
    Else
      If (nstep > ttm%nstepcpl) Then
        Do k=1,ttm%ntcell(3)
          Do j=1,ttm%ntcell(2)
            Do i=1,ttm%ntcell(1)
              ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
              If (ttm%act_ele_cell(ijk,0,0,0)>zero_plus) Then
                ! e-s coupling term
                eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)+tstep*ttm%rvolume/(Ce(ttm%eltemp(ijk,0,0,0),ttm)*&
                  Real(redtstepmx,Kind=wp))*ttm%asource(ijk)
                ! e-p coupling term
                If (ttm%l_epcp) Then
                  Select Case (ttm%gvar)
                   Case (0,1)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                      tstep*ttm%rvolume/(Ce(ttm%eltemp(ijk,0,0,0),ttm)*Real(redtstepmx,Kind=wp))*ttm%gsource(ijk)*&
                      (ttm%eltemp(ijk,0,0,0)-ttm%tempion(ijk))
                   Case (2)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                      tstep*ttm%rvolume/(Ce(ttm%eltemp(ijk,0,0,0),ttm)*Real(redtstepmx,Kind=wp))*ttm%gsource(ijk)*&
                      (ttm%eltemp(ijk,0,0,0)-ttm%tempion(ijk))*&
                      Gep(ttm%eltemp(ijk,0,0,0),ttm)
                  End Select
                End If
              End If
            End Do
          End Do
        End Do
      End If
    End If

! update electronic temperatures to ttm%adjusted values

    Do kk=-ttm%eltcell(3),ttm%eltcell(3)
      Do jj=-ttm%eltcell(2),ttm%eltcell(2)
        Do ii=-ttm%eltcell(1),ttm%eltcell(1)
          Do ijk=1,ttm%numcell
            ttm%eltemp(ijk,ii,jj,kk) = eltemp1(ijk,ii,jj,kk)
          End Do
        End Do
      End Do
    End Do

! update boundary halo values and apply boundary conditions

    Call boundaryHalo(ttm,domain,comm)
    Call boundaryCond (ttm%bcTypeE, temp,ttm,comm)

! simple stability check for simulation

    If (Any(ttm%eltemp < 0.0_wp)) safe = .false.
    Call gcheck(comm,safe)
    If (.not. safe) Call error (683)

  End Do

! Dumping Te file every ndump steps
  Call ttm_system_revive ('DUMP_E',nstep,time,ndump,nstrun,ttm,comm)

  Deallocate (eltemp1, Stat = fail)
  If (fail>0) Call error(1088)
End Subroutine ttm_thermal_diffusion
End Module ttm_track
