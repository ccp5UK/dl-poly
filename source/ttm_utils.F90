Module ttm_utils

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for utility functions for the two-temperature model
! (ttm)
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton july 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use setup
  Use ttm
  Use comms, Only : comms_type,Grid1_tag,Grid2_tag,gsum,gmax,gmin
  Use domains, Only : idx,idy,idz,nprx,npry,nprz
  Use errors_warnings, Only : error
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif
  Implicit None

Contains

  Function Kep(T)

  ! Calcuate interpolated thermal conductivity from tabulated values: given as kB/A^3

    Implicit None

    Real (Kind = wp), Intent(In) :: T
    Real (Kind = wp)             :: Kep

    Call interpolate(kel, ketable, T, Kep)

    return

  End Function Kep

  Function Gep(T)

  ! Calculate electron-phonon coupling term (friction parameter for inhomogeneous
  ! Langevin thermostat): given as ps^-1

    Real (Kind = wp), Intent(In) :: T
    Real (Kind = wp)             :: Gep

    Call interpolate (gel, gtable, T, Gep)
    Gep = Gep * Merge(rcellrho,1.0_wp,ttmdyndens)

    return

  End Function Gep

  Function Ce(T)

! Temperature-dependent specific heat capacity (different for metals and insulators)
! given as kB/A^3: conversions from kB/atom carried out for cases 0, 1 and 2, allowing
! for changes in number of atoms per voxel (cellrho) if dynamic cell density option
! selected

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: T
    Real ( Kind = wp )                :: Ce

    Ce = 1.0_wp

    Select Case (CeType)
    Case (0)
      ! Case 0: constant specific heat capacity (given as kB/A^3)
      Ce = Ce0
    Case (1)
      ! Case 1: hyperbolic tangent specific heat capacity (given as kB/A^3) -
      !         Ce = sh_A*Tanh(T*sh_B*1.0e-4) [kB/atom]
      Ce = sh_A*Tanh(T*sh_B)
    Case (2)
      ! Case 2: linear specific heat capacity to maximum value at/beyond Tfermi
      !         (given as kB/A^3)
      Ce = Min(T/Tfermi,1.0_wp)*Cemax
    Case (4)
      ! Case 4: constant specific heat capacity (converted as kB/A^3)
      Ce = Ce0*cellrho
    Case (5)
      ! Case 5: hyperbolic tangent specific heat capacity (converted to kB/A^3) -
      !         Ce = sh_A*Tanh(T*sh_B*1.0e-4) [kB/atom]
     Ce = sh_A*cellrho*Tanh(T*sh_B)
    Case (6)
      ! Case 6: linear specific heat capacity to maximum value at/beyond Tfermi
      !         (converted to kB/A^3)
      Ce = Min(T/Tfermi,1.0_wp)*Cemax*cellrho
    Case (3,7)
      ! Case 3: interpolated specific heat capacity from table (given as kB/A^3)
      Call interpolate(cel, cetable, T, Ce)
    Case Default
      ! Default case: constant specific heat capacity of 1 kB/atom (convert to kB/A^3)
      Ce = cellrho
    End Select

  End Function Ce

  Function Ke(T)

! Temperature-dependent lattice thermal conductivity: given as kB/(ps A)

    Implicit None
    
    Real ( Kind = wp ) :: Ke
    Real ( Kind = wp ), Intent ( In ) :: T    

    Select Case(KeType)
    Case (3)
    ! Case 3: thermal conductivity interpolated from table
      Ke = Kep(T)
    Case Default
    ! Case 1: constant thermal conductivity
      Ke = Ka0
    End Select

  End Function Ke
  
  
  Function KeD(Te, temp)

! Temperature-dependent Drude-like lattice thermal conductivity: given as kB/(ps A)
! (Ka0 = thermal conductivity at system temperature)

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: Te, temp
    Real ( Kind = wp )                :: KeD

    KeD = Ka0 * Ce(Te) / Ce(temp)
    
  End Function KeD
  
  Function alp(Te)

! Thermal diffusivity: given as A^2/ps
! (different cases for metals and insulators)

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: Te
    Real ( Kind = wp )                :: alp

    Select Case (DeType)
    Case (0)
    ! Case 0: metal system - ratio of conductivity/heat capacity
      alp = Ke(Te) / Ce(Te)
    Case (1)
    ! Case 1: non-metal system - constant value (given as A^2/ps)
      alp = Diff0
    Case (2)
    ! Case 2: non-metal system - reciprocal function of temperature
    !         up to Fermi temperature, Diff0 previously scaled with
    !         system temperature (given as A^2/ps)
      alp = Diff0/Min(Te,Tfermi)
    Case (3)
    ! Case 3: non-metal system - thermal diffusivity interpolated 
    !         from table (given as A^2/ps)
      Call interpolate(del, detable, Te, alp)
    Case Default
      alp = 0.0_wp
    End Select

  End Function alp

  Subroutine calcchies(chi_ep,comm)

! Calculate electron-phonon coupling friction term (chi_ep)
! for homogeneously coupled system: uses mean electronic temperature
! and interpolates from tabulated values given in g.dat file

    Real ( Kind = wp ), Intent ( Inout ) :: chi_ep
    Type (comms_type), Intent ( In )     :: comm
    Real ( Kind = wp )                   :: eltempav = 0.0_wp
    Real ( Kind = wp )                   :: epc = 0.0_wp

    Call eltemp_mean (eltempav, comm)

    Call interpolate (gel, gtable, eltempav, epc)

    chi_ep = epc*Merge(rcellrho, 1.0_wp, ttmdyndens)

  End Subroutine calcchies

  Subroutine boundaryHalo (comm)

! fills halo regions of electronic temperature lattice from neighbouring sections
! (periodic boundary conditions)

    Type(comms_type), Intent(In) :: comm
    Integer :: i,ii,iii1,iii2,j,jj,jjj1,jjj2,k,kk,kkk1,kkk2,ijk1,ijk2
    Integer, Dimension(4) :: req
    Integer, Allocatable :: stats(:,:)

    Integer :: ierr

    Allocate (stats(1:MPI_STATUS_SIZE,1:4))
    If (comm%mxnode>1) Then

      Do kk = -eltcell(3), eltcell(3)
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            If (idx==nprx-1) Then
              iii1 = Mod(ii+3*eltcell(1), (eltcell(1)*2+1)) - eltcell(1)
            Else
              iii1 = ii
            End If
            If (idx==0) Then
              iii2 = Mod(ii+eltcell(1)+1, (eltcell(1)*2+1)) - eltcell(1)
            Else
              iii2 = ii
            End If
            ijk1 = 2 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
            ijk2 = 1 + (ntcell(1)+1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
            Call MPI_ISEND (eltemp(ijk1,ii,jj,kk)  , 1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
            Call MPI_IRECV (eltemp(ijk2,iii1,jj,kk), 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
            ijk1 = 1 + (ntcell(1)) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
            ijk2 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
            Call MPI_ISEND (eltemp(ijk1,ii,jj,kk)  , 1, tmpmsgx, map(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
            Call MPI_IRECV (eltemp(ijk2,iii2,jj,kk), 1, tmpmsgx, map(1), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
            Call MPI_WAITALL (4, req, stats, ierr)
          End Do
        End Do
      End Do

      Do kk = -eltcell(3), eltcell(3)
        Do ii = -eltcell(1), eltcell(1)
          Do jj = -eltcell(2), eltcell(2)
            If (idy==npry-1) Then
              jjj1 = Mod(jj+3*eltcell(2), (eltcell(2)*2+1)) - eltcell(2)
            Else
              jjj1 = jj
            End If
            If (idy==0) Then
              jjj2 = Mod(jj+eltcell(2)+1, (eltcell(2)*2+1)) - eltcell(2)
            Else
              jjj2 = jj
            End If
            ijk1 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
            ijk2 = 1 + (ntcell(1)+2) * (ntcell(2) + 1 + (ntcell(2)+2))
            Call MPI_ISEND (eltemp(ijk1,ii,jj,kk)  , 1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
            Call MPI_IRECV (eltemp(ijk2,ii,jjj1,kk), 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
            ijk1 = 1 + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2))
            ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
            Call MPI_ISEND (eltemp(ijk1,ii,jj,kk)  , 1, tmpmsgy, map(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
            Call MPI_IRECV (eltemp(ijk2,ii,jjj2,kk), 1, tmpmsgy, map(3), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
            Call MPI_WAITALL (4, req, stats, ierr)
          End Do
        End Do
      End Do

      Do jj = -eltcell(2), eltcell(2)
        Do ii = -eltcell(1), eltcell(1)
          Do kk = -eltcell(3), eltcell(3)
            If (idz==nprz-1) Then
              kkk1 = Mod(kk+3*eltcell(3), (eltcell(3)*2+1)) - eltcell(3)
            Else
              kkk1 = kk
            End If
            If (idz==0) Then
              kkk2 = Mod(kk+eltcell(3)+1, (eltcell(3)*2+1)) - eltcell(3)
            Else
              kkk2 = kk
            End If
            ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
            ijk2 = 1 + (ntcell(1)+2) * ((ntcell(2)+2) * (ntcell(3) + 1))
            Call MPI_ISEND (eltemp(ijk1,ii,jj,kk)  , 1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
            Call MPI_IRECV (eltemp(ijk2,ii,jj,kkk1), 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
            ijk1 = 1 + (ntcell(1)+2) * ((ntcell(2)+2) * ntcell(3))
            ijk2 = 1
            Call MPI_ISEND (eltemp(ijk1,ii,jj,kk)  , 1, tmpmsgz, map(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
            Call MPI_IRECV (eltemp(ijk2,ii,jj,kkk2), 1, tmpmsgz, map(5), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
            Call MPI_WAITALL (4, req, stats, ierr)
          End Do
        End Do
      End Do

    Else

      Do kk = -eltcell(3), eltcell(3)
        kkk1 = Mod(kk+3*eltcell(3), (eltcell(3)*2+1)) - eltcell(3)
        kkk2 = Mod(kk+eltcell(3)+1, (eltcell(3)*2+1)) - eltcell(3)
        Do jj = -eltcell(2), eltcell(2)
          jjj1 = Mod(jj+3*eltcell(2), (eltcell(2)*2+1)) - eltcell(2)
          jjj2 = Mod(jj+eltcell(2)+1, (eltcell(2)*2+1)) - eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            iii1 = Mod(ii+3*eltcell(1), (eltcell(1)*2+1)) - eltcell(1)
            iii2 = Mod(ii+eltcell(1)+1, (eltcell(1)*2+1)) - eltcell(1)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 2 + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                ijk2 = 2 + ntcell(1) + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                eltemp(ijk2,iii1,jj,kk) = eltemp(ijk1,ii,jj,kk)
                ijk1 = 1 + ntcell(1) + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                ijk2 = 1 + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                eltemp(ijk2,iii2,jj,kk) = eltemp(ijk1,ii,jj,kk)
              End Do
            End Do            
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (1 + (ntcell(2)+2) * k)
                ijk2 = 1 + i + (ntcell(1)+2) * (ntcell(2) + 1 + (ntcell(2)+2) * k)
                eltemp(ijk2,ii,jjj1,kk) = eltemp(ijk1,ii,jj,kk)
                ijk1 = 1 + i + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2) * k)
                ijk2 = 1 + i + (ntcell(1)+2) * (ntcell(2)+2) * k
                eltemp(ijk2,ii,jjj2,kk) = eltemp(ijk1,ii,jj,kk)
              End Do
            End Do
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2))
                ijk2 = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * (ntcell(3)+1))
                eltemp(ijk2,ii,jj,kkk1) = eltemp(ijk1,ii,jj,kk)
                ijk1 = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * ntcell(3))
                ijk2 = 1 + i + (ntcell(1)+2) * j
                eltemp(ijk2,ii,jj,kkk2) = eltemp(ijk1,ii,jj,kk)
              End Do
            End Do
          End Do
        End Do
      End Do

    End If

    Deallocate (stats)

  End Subroutine boundaryHalo

  Subroutine boundaryCond (key, temp,comm)

! appends halo regions of entire electronic temperature lattice with appropriate boundary conditions


    Real ( Kind = wp ), Intent ( In ) :: temp
    Integer,            Intent ( In ) :: key
    Type(comms_type),   Intent ( In ) :: comm

    Integer :: i,ii,j,jj,k,kk,ijk1,ijk2,ierr
    Integer, Dimension(4) :: req
    Integer, Dimension(MPI_STATUS_SIZE,4) :: stat


    Select Case (key)
  ! Periodic boundary conditions
    Case (1)
      If (ttmbcmap(1)>=0 .or. ttmbcmap(2)>=0) Then
        If (comm%mxnode>1) Then
          Do kk = -eltcell(3), eltcell(3)
            Do jj = -eltcell(2), eltcell(2)
              If (ttmbcmap(1)>=0) Then
                ijk1 = 1 + ttmbc(1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
                ijk2 = 1 + (ttmbc(1) - 1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
                ii = -eltcell(1)
                Call MPI_ISEND (eltemp(ijk1,ii,jj,kk) , 1, tmpmsgx, ttmbcmap(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
                Call MPI_IRECV (eltemp(ijk2,-ii,jj,kk), 1, tmpmsgx, ttmbcmap(1), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
              End If
              If (ttmbcmap(2)>=0) Then
                ijk1 = 1 + ttmbc(2) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
                ijk2 = 1 + (ttmbc(2) + 1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
                ii = eltcell(1)
                Call MPI_ISEND (eltemp(ijk1,ii,jj,kk) , 1, tmpmsgx, ttmbcmap(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
                Call MPI_IRECV (eltemp(ijk2,-ii,jj,kk), 1, tmpmsgx, ttmbcmap(2), Grid1_tag, MPI_COMM_WORLD, req(4), ierr)
              End If
              Call MPI_WAITALL (4, req, stat, ierr)
            End Do
          End Do
        Else
          Do kk = -eltcell(3), eltcell(3)
            Do jj = -eltcell(2), eltcell(2)
              Do k = 1, ntcell(3)
                Do j = 1, ntcell(2)
                  ijk1 = 1 + ttmbc(1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                  ijk2 = 1 + (ttmbc(2) + 1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                  eltemp(ijk2,eltcell(1),jj,kk) = eltemp(ijk1,-eltcell(1),jj,kk)
                  ijk1 = 1 + ttmbc(2) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                  ijk2 = 1 + (ttmbc(1) - 1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                  eltemp(ijk2,-eltcell(1),jj,kk) = eltemp(ijk1,eltcell(1),jj,kk)
                End Do
              End Do
            End Do
          End Do
        End If
      End If

      If (ttmbcmap(3)>=0 .or. ttmbcmap(4)>=0) Then
        If (comm%mxnode>1) Then
          Do kk = -eltcell(3), eltcell(3)
            Do ii = -eltcell(1), eltcell(1)
              If (ttmbcmap(3)>=0) Then
                ijk1 = 1 + (ntcell(1)+2) * (ttmbc(3) + (ntcell(2)+2))
                ijk2 = 1 + (ntcell(1)+2) * (ttmbc(3) - 1 + (ntcell(2)+2))
                jj = -eltcell(2)
                Call MPI_ISEND (eltemp(ijk1,ii,jj,kk) , 1, tmpmsgy, ttmbcmap(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
                Call MPI_IRECV (eltemp(ijk2,ii,-jj,kk), 1, tmpmsgy, ttmbcmap(3), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
              End If
              If (ttmbcmap(4)>=0) Then
                ijk1 = 1 + (ntcell(1)+2) * (ttmbc(4) + (ntcell(2)+2))
                ijk2 = 1 + (ntcell(1)+2) * (ttmbc(4) + 1 + (ntcell(2)+2))
                jj = eltcell(2)
                Call MPI_ISEND (eltemp(ijk1,ii,jj,kk) , 1, tmpmsgy, ttmbcmap(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
                Call MPI_IRECV (eltemp(ijk2,ii,-jj,kk), 1, tmpmsgy, ttmbcmap(4), Grid1_tag, MPI_COMM_WORLD, req(4), ierr)
              End If
              Call MPI_WAITALL (4, req, stat, ierr)
            End Do
          End Do
        Else
          Do kk = -eltcell(3), eltcell(3)
            Do ii = -eltcell(1), eltcell(1)
              Do k = 1, ntcell(3)
                Do i = 0, ntcell(1)+1
                  ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(3) + k * (ntcell(2)+2))
                  ijk2 = 1 + i + (ntcell(1)+2) * (ttmbc(4) + 1 + k * (ntcell(2)+2))
                  eltemp(ijk2,ii,eltcell(2),kk) = eltemp(ijk1,ii,-eltcell(2),kk)
                  ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(4) + k * (ntcell(2)+2))
                  ijk2 = 1 + i + (ntcell(1)+2) * (ttmbc(3) - 1 + k * (ntcell(2)+2))
                  eltemp(ijk2,ii,-eltcell(2),kk) = eltemp(ijk1,ii,eltcell(2),kk)
                End Do
              End Do
            End Do
          End Do
        End If
      End If

      If (ttmbcmap(5)>=0 .or. ttmbcmap(6)>=0) Then
        If (comm%mxnode>1) Then
          Do jj = -eltcell(2), eltcell(2)
            Do ii = -eltcell(1), eltcell(1)
              If (ttmbcmap(5)>=0) Then
                ijk1 = 1 + (ntcell(1)+2) * (ttmbc(5) * (ntcell(2)+2))
                ijk2 = 1 + (ntcell(1)+2) * ((ttmbc(5) - 1) * (ntcell(2)+2))
                kk = -eltcell(3)
                Call MPI_ISEND (eltemp(ijk1,ii,jj,kk) , 1, tmpmsgz, ttmbcmap(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
                Call MPI_IRECV (eltemp(ijk2,ii,jj,-kk), 1, tmpmsgz, ttmbcmap(5), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
              End If
              If (ttmbcmap(6)>=0) Then
                ijk1 = 1 + (ntcell(1)+2) * (ttmbc(6) * (ntcell(2)+2))
                ijk2 = 1 + (ntcell(1)+2) * ((ttmbc(6) + 1) * (ntcell(2)+2))
                kk = eltcell(3)
                Call MPI_ISEND (eltemp(ijk1,ii,jj,kk) , 1, tmpmsgz, ttmbcmap(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
                Call MPI_IRECV (eltemp(ijk2,ii,jj,-kk), 1, tmpmsgz, ttmbcmap(6), Grid1_tag, MPI_COMM_WORLD, req(4), ierr)
              End If
              Call MPI_WAITALL (4, req, stat, ierr)
            End Do
          End Do
        Else
          Do jj = -eltcell(2), eltcell(2)
            Do ii = -eltcell(1), eltcell(1)
              Do j = 0, ntcell(2)+1
                Do i = 0, ntcell(1)+1
                  ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(5) * (ntcell(2)+2))
                  ijk2 = 1 + i + (ntcell(1)+2) * (j + (ttmbc(6) + 1) * (ntcell(2)+2))
                  eltemp(ijk2,ii,jj,eltcell(3)) = eltemp(ijk1,ii,jj,-eltcell(3))
                  ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(6) * (ntcell(2)+2))
                  ijk2 = 1 + i + (ntcell(1)+2) * (j + (ttmbc(5) - 1) * (ntcell(2)+2))
                  eltemp(ijk2,ii,jj,-eltcell(3)) = eltemp(ijk1,ii,jj,eltcell(3))
                End Do
              End Do
            End Do
          End Do
        End If

      End If

  ! Infinite sink/source (Dirichlet) boundary conditions
    Case (2)
      If (ttmbcmap(1)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk2 = 1 + (ttmbc(1)-1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                eltemp(ijk2,-eltcell(1),jj,kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(2)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk2 = 1 + (ttmbc(2)+1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                eltemp(ijk2,eltcell(1),jj,kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(3)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk2 = 1 + i + (ntcell(1)+2) * ((ttmbc(3)-1) + k * (ntcell(2)+2))
                eltemp(ijk2,ii,-eltcell(2),kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(4)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk2 = 1 + i + (ntcell(1)+2) * ((ttmbc(4)+1) + k * (ntcell(2)+2))
                eltemp(ijk2,ii,eltcell(2),kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(5)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk2 = 1 + i + (ntcell(1)+2) * (j + (ttmbc(5)-1) * (ntcell(2)+2))
                eltemp(ijk2,ii,jj,-eltcell(3)) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(6)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk2 = 1 + i + (ntcell(1)+2) * (j + (ttmbc(6)+1) * (ntcell(2)+2))
                eltemp(ijk2,ii,jj,eltcell(3)) = temp
              End Do
            End Do
          End Do
        End Do
      End If

  ! 'Confined' (von Neumann) boundary conditions
    Case (3)
      If (ttmbcmap(1)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 1 + ttmbc(1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                ijk2 = ijk1 - 1
                eltemp(ijk2,-eltcell(1),jj,kk) = eltemp(ijk1,-eltcell(1),jj,kk)
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(2)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 1 + ttmbc(2) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                ijk2 = ijk1 + 1
                eltemp(ijk2,eltcell(1),jj,kk) = eltemp(ijk1,eltcell(1),jj,kk)
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(3)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(3) + k * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)
                eltemp(ijk2,ii,-eltcell(2),kk) = eltemp(ijk1,ii,-eltcell(2),kk)
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(4)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(4) + k * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)
                eltemp(ijk2,ii,eltcell(2),kk) = eltemp(ijk1,ii,eltcell(2),kk)
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(5)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(5) * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,-eltcell(3)) = eltemp(ijk1,ii,jj,-eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(6)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(6) * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,eltcell(3)) = eltemp(ijk1,ii,jj,eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

  ! Mixed case: Infinite sink/source (Dirichlet) boundaries in x/y-directions
  !             'Confined' (von Neumann) boundary in z-direction
    Case (4)
      If (ttmbcmap(1)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk2 = 1 + (ttmbc(1)-1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                eltemp(ijk2,-eltcell(1),jj,kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(2)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk2 = 1 + (ttmbc(2)+1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                eltemp(ijk2,eltcell(1),jj,kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(3)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk2 = 1 + i + (ntcell(1)+2) * ((ttmbc(3)-1) + k * (ntcell(2)+2))
                eltemp(ijk2,ii,-eltcell(2),kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(4)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk2 = 1 + i + (ntcell(1)+2) * ((ttmbc(4)+1) + k * (ntcell(2)+2))
                eltemp(ijk2,ii,eltcell(2),kk) = temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(5)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(5) * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,-eltcell(3)) = eltemp(ijk1,ii,jj,-eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(6)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(6) * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,eltcell(3)) = eltemp(ijk1,ii,jj,eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

  ! Robin boundary conditions
    Case (5)
      If (ttmbcmap(1)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 1 + ttmbc(1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                ijk2 = ijk1 - 1
                eltemp(ijk2,-eltcell(1),jj,kk) = fluxout*(eltemp(ijk1,-eltcell(1),jj,kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(2)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 1 + ttmbc(2) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                ijk2 = ijk1 + 1
                eltemp(ijk2,eltcell(1),jj,kk) = fluxout*(eltemp(ijk1,eltcell(1),jj,kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(3)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(3) + k * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)
                eltemp(ijk2,ii,-eltcell(2),kk) = fluxout*(eltemp(ijk1,ii,-eltcell(2),kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(4)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(4) + k * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)
                eltemp(ijk2,ii,eltcell(2),kk) = fluxout*(eltemp(ijk1,ii,eltcell(2),kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(5)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(5) * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,-eltcell(3)) = fluxout*(eltemp(ijk1,ii,jj,-eltcell(3))-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(6)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(6) * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,eltcell(3)) = fluxout*(eltemp(ijk1,ii,jj,eltcell(3))-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

  ! Mixed case: Robin boundaries in x/y-directions
  !             'Confined' (von Neumann) boundary in z-direction
    Case (6)
      If (ttmbcmap(1)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 1 + ttmbc(1) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                ijk2 = ijk1 - 1
                eltemp(ijk2,-eltcell(1),jj,kk) = fluxout*(eltemp(ijk1,-eltcell(1),jj,kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(2)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do jj = -eltcell(2), eltcell(2)
            Do k = 1, ntcell(3)
              Do j = 1, ntcell(2)
                ijk1 = 1 + ttmbc(2) + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
                ijk2 = ijk1 + 1
                eltemp(ijk2,eltcell(1),jj,kk) = fluxout*(eltemp(ijk1,eltcell(1),jj,kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(3)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(3) + k * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)
                eltemp(ijk2,ii,-eltcell(2),kk) = fluxout*(eltemp(ijk1,ii,-eltcell(2),kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(4)>=0) Then
        Do kk = -eltcell(3), eltcell(3)
          Do ii = -eltcell(1), eltcell(1)
            Do k = 1, ntcell(3)
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (ttmbc(4) + k * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)
                eltemp(ijk2,ii,eltcell(2),kk) = fluxout*(eltemp(ijk1,ii,eltcell(2),kk)-temp) + temp
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(5)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(5) * (ntcell(2)+2))
                ijk2 = ijk1 - (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,-eltcell(3)) = eltemp(ijk1,ii,jj,-eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

      If (ttmbcmap(6)>=0) Then
        Do jj = -eltcell(2), eltcell(2)
          Do ii = -eltcell(1), eltcell(1)
            Do j = 0, ntcell(2)+1
              Do i = 0, ntcell(1)+1
                ijk1 = 1 + i + (ntcell(1)+2) * (j + ttmbc(6) * (ntcell(2)+2))
                ijk2 = ijk1 + (ntcell(1)+2)*(ntcell(2)+2)
                eltemp(ijk2,ii,jj,eltcell(3)) = eltemp(ijk1,ii,jj,eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

    End Select

  End Subroutine boundaryCond

  Subroutine peakProfiler(lat, peakfile, nstep, freq,comm)

! prints a slice across y (centred in xz-plane) of a lattice to a file

    Real ( Kind = wp ), Dimension (numcell), Intent( In ) :: lat
    Real ( Kind = wp ), Dimension (ntsys(2)) :: laty
    Integer, Intent( In ) :: nstep,freq
    Character ( Len = * ), Intent( In ) :: peakfile
    Type(comms_type), Intent ( InOut ) :: comm
    Integer :: iounit = 114
    Integer :: i,j,k,l,ijk

    laty = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0) Then
        i = midI(1) - ntcelloff(1)
        k = midI(3) - ntcelloff(3)
        If (i>0 .and. i<=ntcell(1) .and. k>0 .and. k<=ntcell(3)) Then
          Do j=1,ntcell(2)
            l = j + ntcelloff(2)
            ijk = 1 + i + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
            laty(l) = lat(ijk)
          End Do
        End If
        Call gsum (comm,laty)
        If (comm%idnode==0) Then
          If (nstep==0) Then
            Open (Unit=iounit, File=peakfile, Status='replace')
          Else
            Open (Unit=iounit, File=peakfile, Position='append')
          End If
          Do i=1,ntsys(2)
            Write(iounit,Fmt='(9es12.4,1p)', Advance='no') laty(i)
          End Do
          Close (iounit)
        End If
      End If
    End If

  End Subroutine peakProfiler
 
  Subroutine peakProfilerElec(peakfile, nstep, freq, comm)

! prints a slice across y (centred in xz-plane) of the electronic temperature lattice to a file

    Real ( Kind = wp ), Dimension (eltsys(2)) :: laty
    Integer, Intent( In ) :: nstep,freq
    Character ( Len = * ), Intent( In ) :: peakfile
    Type(comms_type), Intent( InOut ) :: comm
    Integer :: iounit = 114
    Integer :: i,j,jj,jmin,jmax,k,l,ijk

    laty = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0) Then
        i = midI(1) - ntcelloff(1)
        k = midI(3) - ntcelloff(3)
        If (i>0 .and. i<=ntcell(1) .and. k>0 .and. k<=ntcell(3)) Then
          Do jj=-eltcell(2),eltcell(2)
            If (eltcell(2)>0 .and. jj==-eltcell(2) .and. ttmbcmap(3)>=0) Then
              jmin = ttmbc(3)
            Else
              jmin = 1
            End If
            If (eltcell(2)>0 .and. jj==eltcell(2) .and. ttmbcmap(4)>=0) Then
              jmax = ttmbc(4)
            Else
              jmax = ntcell(2)
            End If
            Do j=jmin,jmax
              l = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
              ijk = 1 + i + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
              If (l>0 .and. l<=eltsys(2)) laty(l) = eltemp(ijk,0,jj,0)
            End Do
          End Do
        End If
        Call gsum (comm,laty)
        If (comm%idnode==0) Then
          If (nstep==0) Then
            Open (Unit=iounit, File=peakfile, Status='replace')
          Else
            Open (Unit=iounit, File=peakfile, Position='append')
          End If
          Do i=1,eltsys(2)
            Write(iounit,Fmt='(9es12.4,1p)', Advance='no') laty(i)
          End Do
          Close (iounit)
        End If
      End If
    End If

  End Subroutine peakProfilerElec

  Subroutine printLatticeStatsToFile(lat, latfile, time, nstep, freq, comm)

! prints lattice statistics (minimum, maximum, sum) to file

   
    Real ( Kind = wp ), Dimension(numcell), Intent ( In ) :: lat
    Integer, Intent( In ) :: nstep,freq
    Real ( Kind = wp ), Intent ( In ) :: time
    Character ( Len = * ), Intent ( In ) :: latfile
    Type(comms_type), Intent ( InOut ) :: comm

    Real (kind = wp) :: lat_sum,lat_min,lat_max,rtotal
    Integer :: iounit = 113
    Integer :: i,j,k,ijk

    lat_sum=0.0_wp
    lat_min=1.0e30_wp
    lat_max=0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0) Then

        ! check over all active ionic temperature cells (CIT)
        Do k = 1,ntcell(3)
          Do j = 1,ntcell(2)
            Do i = 1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
              If (act_ele_cell(ijk,0,0,0)>zero_plus) Then
                lat_sum = lat_sum + lat(ijk)
                lat_min = Min(lat_min, lat(ijk))
                lat_max = Max(lat_max, lat(ijk))
              End If
            End Do
          End Do
        End Do

        Call gsum(comm,lat_sum)
        Call gmin(comm,lat_min)
        Call gmax(comm,lat_max)

        rtotal = Merge(1.0_wp/Real(acell,Kind=wp),0.0_wp,(acell>0))

        If (comm%idnode == 0) Then
          If (nstep==0) Then
            Open(Unit=iounit, File=latfile, Status='replace')
          Else
            Open(Unit=iounit, File=latfile, Position='append')
          End If
          Write(iounit,'(i8,5(1p,es12.4))') nstep, time, lat_min, lat_max, lat_sum*rtotal, lat_sum
          Close(iounit)
        End If

      End If
    End If
      
  End Subroutine printLatticeStatsToFile
  
  Subroutine printElecLatticeStatsToFile(latfile, time, temp0, nstep, freq, comm)

! prints electronic temperature lattice statistics (minimum, maximum, sum) 
! and energy (E = integral of Ce(Te)*Te between temp0 and Te) to file

    Integer, Intent( In ) :: nstep,freq
    Real ( Kind = wp ), Intent ( In ) :: time,temp0
    Character ( Len = * ), Intent ( In ) :: latfile
    Type( comms_type), Intent ( InOut ) :: comm

    Real ( Kind = wp ) :: lat_sum,lat_min,lat_max,Ue,tmp,sgnplus,eltmp,totalcell,rtotal
    Real ( Kind = wp ) :: Ce0a,sh_Aa,Cemaxa
    Integer :: iounit = 115
    Integer :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax
    Integer :: ijk,numint,n,lx,ly,lz

    Ce0a   = Ce0  *Merge(cellrho,1.0_wp,ttmdyndens)
    sh_Aa  = sh_A *Merge(cellrho,1.0_wp,ttmdyndens)
    Cemaxa = Cemax*Merge(cellrho,1.0_wp,ttmdyndens)

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0) Then

        Call eltemp_sum(lat_sum,comm)
        Call eltemp_min(lat_min,comm)
        Call eltemp_max(lat_max,comm)

        Ue = 0.0_wp
        totalcell = 0.0_wp

        Do kk = -eltcell(3), eltcell(3)
          If (eltcell(3)>0 .and. kk == -eltcell(3) .and. ttmbcmap(5)>=0) Then
            kmin = ttmbc(5)
          Else
            kmin = 1
          End If

          If (eltcell(3)>0 .and. kk == eltcell(3) .and. ttmbcmap(6)>=0) Then
            kmax = ttmbc(6)
          Else
            kmax = ntcell(3)
          End If

          Do jj = -eltcell(2), eltcell(2)
            If (eltcell(2)>0 .and. jj == -eltcell(2) .and. ttmbcmap(3)>=0) Then
              jmin = ttmbc(3)
            Else
              jmin = 1
            End If
            If (eltcell(2)>0 .and. jj == eltcell(2) .and. ttmbcmap(4)>=0) Then
              jmax = ttmbc(4)
            Else
              jmax = ntcell(2)
            End If

            Do ii = -eltcell(1), eltcell(1)
              If (eltcell(1)>0 .and. ii == -eltcell(1) .and. ttmbcmap(1)>=0) Then
                imin = ttmbc(1)
              Else
                imin = 1
              End If
              If (eltcell(1)>0 .and. ii == eltcell(1) .and. ttmbcmap(2)>=0) Then
                imax = ttmbc(2)
              Else
                imax = ntcell(1)
              End If
              Do k = kmin, kmax
                lz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
                Do j = jmin, jmax
                  ly = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
                  Do i = imin, imax
                    lx = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    ! integrate electronic heat capacity from temp0 to electronic temperature
                    ! of cell if it is within global range and active (within ionic temperature cells)
                    If (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3)) Then
                      eltmp = eltemp(ijk,ii,jj,kk)
                      tmp = Merge(1.0_wp,0.0_wp,(act_ele_cell(ijk,0,0,0)>zero_plus .or. (ii/=0 .or. jj/=0 .or. kk/=0)))
                      totalcell = totalcell + tmp
                      Select Case (CeType)
                      Case (0,4)
                      ! constant specific heat capacity
                        tmp = tmp*Ce0a*(eltmp-temp0)
                      Case (1,5)
                      ! hyperbolic tangent specific heat capacity
                        tmp = tmp*sh_Aa*Log(Cosh(sh_B*eltmp)/Cosh(sh_B*temp0))/sh_B
                      Case (2,6)
                      ! linear specific heat capacity to Fermi temperature
                        tmp = tmp*Cemaxa*(0.5_wp*((Min(Tfermi,eltmp))**2-temp0*temp0)/Tfermi+Max(eltmp-Tfermi,0.0_wp))
                      Case Default
                      ! tabulated volumetric heat capacity or more complex 
                      ! functions: integrate using trapezium rule
                      ! with final interval of <=1 kelvin
                        numint = Floor(tmp*(eltmp-temp0))
                        tmp = 0.0_wp
                        sgnplus = Sign(1.0_wp,Real(numint,Kind=wp))
                        Do n=1,Abs(numint)
                          tmp = tmp+0.5_wp*sgnplus*(Ce(temp0+sgnplus*Real(n-1,Kind=wp))+Ce(temp0+sgnplus*Real(n,Kind=wp)))
                        End Do
                        tmp = tmp+0.5_wp*(Ce(temp0+Real(numint,Kind=wp))+Ce(eltmp))
                      End Select
                      Ue = Ue+tmp*volume*kB_to_eV
                    End If
                  End Do
                End Do
              End Do

            End Do
          End Do
        End Do

        Call gsum(comm,Ue)
        Call gsum(comm,totalcell)

        rtotal = Merge(1.0_wp/totalcell,0.0_wp,(totalcell>zero_plus))

        If (comm%idnode == 0) Then
          If (nstep==0) Then
            Open(Unit=iounit, File=latfile, Status='replace')
          Else
            Open(Unit=iounit, File=latfile, Position='append')
          End If
          Write(iounit,'(i8,6(1p,es12.4))') nstep, time, lat_min, lat_max, lat_sum*rtotal, lat_sum, Ue
          Close(iounit)
        End If

      End If
    End If
      
  End Subroutine printElecLatticeStatsToFile
  
  
  Subroutine interpolate(tabsize, tab, x0, resp)

! Interpolation helper subroutine

    Integer,                                      Intent( In )  :: tabsize
    Real( Kind = wp ), Dimension (1:tabsize,1:2), Intent( In )  :: tab
    Real( Kind = wp ),                            Intent( In )  :: x0
    Real( Kind = wp ),                            Intent( Out ) :: resp

    Integer :: i,j,k, ierr

    If (tabsize == 1) Then
      resp = tab(1,2)
      Return
    End If

    resp = 0.0_wp

! if input is within range, find nearest entry in table
! below this value by means of bisection - calculate
! output by linear interpolation between this entry
! and the next

    If ((x0>tab(1,1)) .and. (x0<tab(tabsize,1))) Then
      i=1
      j=tabsize
      Do
        k=(i+j)/2
        If (x0 < tab(k,1)) Then
          j=k
        Else
          i=k
        End If
        If (i+1 >= j) Then
          k = j-1
          Exit
        End If
      End Do
      resp = tab(k,2)+(x0-tab(k,1))*(tab(k+1,2)-tab(k,2))/(tab(k+1,1)-tab(k,1))

! if input is over maximum value, set output to value at maximum

    Else If (x0>=tab(tabsize,1)) Then

      resp = tab(tabsize,2)

! If input is under minimum value, set output to value at minimum

    Else If (x0<=tab(1,1)) Then

      resp = tab(1,2)

! If input is completely unusable, terminate with error

    Else

      Call error(683)

    End If

  End Subroutine interpolate

  Subroutine redistribute_Te (temp0,comm)

! Redistribute electronic energy when electronic temperature voxels are closed

    Real( Kind = wp ), Intent ( In ) :: temp0
    Type( comms_type), Intent ( In ) :: comm
    Integer :: i, j, k, ii, jj, kk, ijk, ijk1, ijk2, ijkpx, ijkmx, ijkpy, ijkmy, &
               ijkpz, ijkmz, n, numint, ierr
    Real( Kind = wp ) :: energy_per_cell, U_e, start_Te, end_Te, increase, oldCe, newCe
    Real( Kind = wp ) :: energy_diff, act_sur_cells, sgnplus, Ce0a, sh_Aa, Cemaxa

    Integer, Dimension(4) :: req
    Integer, Dimension(MPI_STATUS_SIZE,4) :: stat
    Real( Kind = wp ), Allocatable :: energydist (:,:,:,:), buffer (:,:,:,:)

    Ce0a   = Ce0  *Merge(cellrho,1.0_wp,ttmdyndens)
    sh_Aa  = sh_Aa*Merge(cellrho,1.0_wp,ttmdyndens)
    Cemaxa = Cemax*Merge(cellrho,1.0_wp,ttmdyndens)

    Allocate (energydist (numcell,-1:1,-1:1,-1:1))
    energydist = 0.0_wp
    eltemp_adj = 0.0_wp
    adjust = .false.

    Do k = 1, ntcell(3)
      Do j = 1, ntcell(2)
        Do i = 1, ntcell(1)
          ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          ijkpx = 2 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          ijkmx = i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          ijkpy = 1 + i + (ntcell(1)+2) * (j + 1 + (ntcell(2)+2) * k)
          ijkmy = 1 + i + (ntcell(1)+2) * (j - 1 + (ntcell(2)+2) * k)
          ijkpz = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * (k + 1))
          ijkmz = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * (k - 1))
          If (act_ele_cell(ijk,0,0,0)<=zero_plus .and. old_ele_cell(ijk,0,0,0)>zero_plus) Then
        ! Calculate amount of energy left in the cell 
        ! (Ue = integral of Ce(Te)d(Te) between temp0 and Te)
            U_e = 0.0_wp
            end_Te = eltemp(ijk,0,0,0)
            Select Case (CeType)
            Case (0,4)
            ! constant specific heat capacity
              U_e = Ce0a*(end_Te-temp0)
            Case (1,5)
            ! hyperbolic tangent specific heat capacity
              U_e = sh_Aa*Log(Cosh(sh_B*end_Te)/Cosh(sh_B*temp0))/sh_B
            Case (2,6)
            ! linear specific heat capacity to Fermi temperature
              increase = Min(Tfermi,end_Te)
              U_e = Cemaxa*(0.5_wp*(increase*increase-temp0*temp0)/Tfermi+Max(end_Te-Tfermi,0.0_wp))
            Case Default
            ! tabulated volumetric heat capacity or more complex 
            ! functions: integrate using trapezium rule
            ! with final interval of <=1 kelvin
              numint = Floor(end_Te - temp0)
              sgnplus = Sign (1.0_wp, Real(numint, Kind=wp))
              Do n = 1,Abs(numint)
                U_e = U_e + 0.5_wp * sgnplus * (Ce(temp0+sgnplus*Real(n-1, Kind=wp))+Ce(temp0+sgnplus*Real(n, Kind=wp)))
              End Do
              U_e = U_e + 0.5_wp * (Ce(temp0+Real(numint, Kind=wp))+Ce(end_Te))
            End Select
            ! Check how many cells are connected to the now turned-off cell
            act_sur_cells = act_ele_cell(ijkmx,0,0,0) + act_ele_cell(ijkpx,0,0,0) + act_ele_cell(ijkmy,0,0,0) + &
                            act_ele_cell(ijkpy,0,0,0) + act_ele_cell(ijkmz,0,0,0) + act_ele_cell(ijkpz,0,0,0)
            ! Calculate energy redistribution to each cell and assign to cells
            If (act_sur_cells>zero_plus) energy_per_cell = U_e / act_sur_cells
            If (act_ele_cell(ijkmx,0,0,0)>zero_plus) energydist (ijkmx,0,0,0) = energydist (ijkmx,0,0,0) + energy_per_cell
            If (act_ele_cell(ijkpx,0,0,0)>zero_plus) energydist (ijkpx,0,0,0) = energydist (ijkpx,0,0,0) + energy_per_cell
            If (act_ele_cell(ijkmy,0,0,0)>zero_plus) energydist (ijkmy,0,0,0) = energydist (ijkmy,0,0,0) + energy_per_cell
            If (act_ele_cell(ijkpy,0,0,0)>zero_plus) energydist (ijkpy,0,0,0) = energydist (ijkpy,0,0,0) + energy_per_cell
            If (act_ele_cell(ijkmz,0,0,0)>zero_plus) energydist (ijkmz,0,0,0) = energydist (ijkmz,0,0,0) + energy_per_cell
            If (act_ele_cell(ijkpz,0,0,0)>zero_plus) energydist (ijkpz,0,0,0) = energydist (ijkpz,0,0,0) + energy_per_cell
          End If
        End Do
      End Do
    End Do

    If(comm%mxnode>1) Then

      Allocate (buffer (numcell,-1:1,-1:1,-1:1))
      buffer = 0.0_wp

      ! sum up redistributed energies (placing energies for electronic grid
      ! outside ionic grid in boundary halos: removing previous values in halos 
      ! before summation)

      ! -x/+x direction
      ijk1 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
      ijk2 = 1 + ntcell(1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
      If (idx==nprx-1) Then
        ii = -1
      Else
        ii = 0
      End If
      Call MPI_ISEND (energydist(ijk1,0,0,0), 1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV (buffer(ijk2,ii,0,0)   , 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
      ijk1 = (ntcell(1)+2) * (2 + (ntcell(2)+2))
      ijk2 = 2 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
      If (idx==0) Then
        ii = 1
      Else
        ii = 0
      End If
      Call MPI_ISEND (energydist(ijk1,0,0,0), 1, tmpmsgx, map(2), Grid2_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV (buffer(ijk2,ii,0,0)   , 1, tmpmsgx, map(1), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
      Call MPI_WAITALL (4, req, stat, ierr)

      ! -y/+y direction
      ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
      ijk2 = 1 + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2))
      If (idy==npry-1) Then
        jj = -1
      Else
        jj = 0
      End If
      Call MPI_ISEND (energydist(ijk1,0,0,0), 1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV (buffer(ijk2,0,jj,0)   , 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
      ijk1 = 1 + (ntcell(1)+2) * (ntcell(2) + 1 + (ntcell(2)+2))
      ijk2 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
      If (idy==0) Then
        jj = 1
      Else
        jj = 0
      End If
      Call MPI_ISEND (energydist(ijk1,0,0,0), 1, tmpmsgy, map(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_IRECV (buffer(ijk2,0,jj,0)   , 1, tmpmsgy, map(3), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
      Call MPI_WAITALL (4, req, stat, ierr)

      ! -z/+z direction
      ijk1 = 1
      ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2) * ntcell(3)
      If (idz==nprz-1) Then
        kk = -1
      Else
        kk = 0
      End If
      Call MPI_ISEND (energydist(ijk1,0,0,0), 1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_IRECV (buffer(ijk2,0,0,kk)   , 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
      ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2) * (ntcell(3)+1)
      ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
      If (idz==0) Then
        kk = 1
      Else
        kk = 0
      End If
      Call MPI_ISEND (energydist(ijk1,0,0,0), 1, tmpmsgz, map(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_IRECV (buffer(ijk2,0,0,kk)   , 1, tmpmsgz, map(5), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
      Call MPI_WAITALL (4, req, stat, ierr)

      energydist = energydist + buffer
      Deallocate (buffer)
    Else
    ! -x/+x direction
      Do k = 1, ntcell(3)
        Do j = 1, ntcell(2)
          ijk1 = 1 + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          ijk2 = 1 + ntcell(1) + (ntcell(1)+2) * (1 + j + (ntcell(2)+2) * k)
          energydist (ijk2,-1,0,0) = energydist (ijk1,0,0,0)
          ijk1 = (ntcell(1)+2) * (1 + j + (ntcell(2)+2) * k)
          ijk2 = 2 + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          energydist (ijk2,1,0,0) = energydist (ijk1,0,0,0)
        End Do
      End Do
    ! -y/+y direction
      Do k = 1, ntcell(3)
        Do i = 0, ntcell(1)+1
          ijk1 = 1 + i + (ntcell(1)+2) * (ntcell(2)+2) * k
          ijk2 = 1 + i + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2) * k)
          energydist (ijk2,0,-1,0) = energydist (ijk1,0,0,0)
          ijk1 = 1 + i + (ntcell(1)+2) * (ntcell(2)+1 + (ntcell(2)+2) * k)
          ijk2 = 1 + i + (ntcell(1)+2) * (1 + (ntcell(2)+2) * k)
          energydist (ijk2,0,1,0) = energydist (ijk1,0,0,0)
        End Do
      End Do
    ! -z/+z direction
      Do j = 0, ntcell(2)+1
        Do i = 0, ntcell(1)+1
          ijk1 = 1 + i + (ntcell(1)+2) * j
          ijk2 = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * ntcell(3))
          energydist (ijk2,0,0,-1) = energydist (ijk1,0,0,0)
          ijk1 = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * (ntcell(3)+1))
          ijk2 = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2))
          energydist (ijk2,0,1,1) = energydist (ijk1,0,0,0)
        End Do
      End Do
    End If

    ! work through electronic temperature cells within ionic temperature voxels
    ! and in immediate borders: determine electronic temperatures needed to add
    ! redistributed energies and identify relevant cells - split according to
    ! available functional form of heat capacity

    Select Case (CeType)
    Case (0,4)
    ! constant specific heat capacity
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, numcell
              If (Abs(energydist(ijk,ii,jj,kk))>zero_plus) Then
                energy_per_cell = energydist(ijk,ii,jj,kk)
                start_Te = eltemp(ijk,ii,jj,kk)
                end_Te = start_Te+energy_per_cell/Ce0a
                eltemp_adj(ijk,ii,jj,kk) = end_Te
                adjust(ijk,ii,jj,kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    Case (1,5)
    ! hyperbolic tangent specific heat capacity
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, numcell
              If (Abs(energydist(ijk,ii,jj,kk))>zero_plus) Then
                energy_per_cell = energydist(ijk,ii,jj,kk)
                start_Te = eltemp(ijk,ii,jj,kk)
                increase = Cosh(sh_B*start_Te)*Exp(sh_B*energy_per_cell/sh_Aa)
                ! using equivalent function: Acosh(x)=Log(x+Sqrt((x-1.0)*(x+1.0)))
                end_Te = Log(increase+Sqrt((increase-1.0_wp)*(increase+1.0_wp)))/sh_B
                eltemp_adj(ijk,ii,jj,kk) = end_Te
                adjust(ijk,ii,jj,kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    Case (2,6)
    ! linear specific heat capacity to Fermi temperature
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, numcell
              If (Abs(energydist(ijk,ii,jj,kk))>zero_plus) Then
                energy_per_cell = energydist(ijk,ii,jj,kk)
                start_Te = eltemp(ijk,ii,jj,kk)
                If (energy_per_cell>zero_plus) Then
                  end_Te = Sqrt(start_Te*start_Te+2.0_wp*energy_per_cell*Tfermi/Cemaxa)
                  If (end_Te>Tfermi) end_Te = 0.5_wp*(start_Te*start_Te/Tfermi+Tfermi)+energy_per_cell/Cemaxa
                Else
                  end_Te = start_Te + energy_per_cell/Cemaxa
                  If (end_Te<Tfermi) end_Te = Sqrt(Tfermi*(2.0_wp*(start_Te+energy_per_cell/Cemaxa)-Tfermi))
                End If
                eltemp_adj(ijk,ii,jj,kk) = end_Te
                adjust(ijk,ii,jj,kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    Case Default
    ! tabulated volumetric heat capacity or more complex
    ! function: find new temperature iteratively by
    ! gradual integration (1 kelvin at a time)
    ! and interpolate within last kelvin
      Do kk = -1, 1
        Do jj = -1, 1
          Do ii = -1, 1
            Do ijk = 1, numcell
              If (Abs(energydist(ijk,ii,jj,kk))>zero_plus) Then
                energy_per_cell = energydist(ijk,ii,jj,kk)
                start_Te = eltemp(ijk,ii,jj,kk)
                sgnplus = Sign(1.0_wp,energy_per_cell)
                energy_diff = sgnplus*energy_per_cell
                oldCe = Ce(start_Te)
                ! increase/decrease temperature of electronic cell by 1 kelvin until
                ! required energy change is reached or exceeded
                Do While(energy_diff>=zero_plus)
                  newCe = Ce(start_Te+sgnplus)
                  increase = 0.5_wp*(oldCe+newCe)
                  energy_diff = energy_diff-increase
                  start_Te = start_Te+sgnplus
                  oldCe = newCe
                End Do
                ! now interpolate for new temperature by using average specific heat
                ! over surrounding one-kelvin interval
                energy_diff = energy_diff+increase
                start_Te = start_Te-sgnplus
                newCe = Ce(start_Te)
                end_Te = start_Te+2.0_wp*sgnplus*energy_diff/(oldCe+newCe)
                eltemp_adj(ijk,ii,jj,kk) = end_Te
                adjust(ijk,ii,jj,kk) = .true.
              End If
            End Do
          End Do
        End Do
      End Do
    End Select

    Deallocate (energydist)

  End Subroutine redistribute_Te

End Module ttm_utils
