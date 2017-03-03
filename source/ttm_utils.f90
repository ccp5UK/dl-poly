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

  Use setup_module
  Use ttm_module
  Use comms_module
  Use domains_module, Only : idx,idy,idz,nprx,npry,nprz
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

! [GK] Define function to calculate e-p coupling for given electronic temperature
  Function Gep(T)

  ! Calculate electron-phonon coupling term (friction parameter for inhomogeneous
  ! Langevin thermostat): given as ps^-1
    Implicit None

    Real (Kind = wp), Intent(In) :: T
    Real (Kind = wp)             :: Gep

    Call interpolate(gel, gtable, T, Gep)

    return

  End Function Gep

  Function Ce(T)

! Temperature-dependent specific heat capacity (different for metals and insulators)
! given as kB/A^3: all conversions from inputs carried out during setup of ttm

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: T
    Real ( Kind = wp )                :: Ce

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
    Case (3)
      ! Case 3: interpolated specific heat capacity from table (given as kB/A^3)
      Call interpolate(cel, cetable, T, Ce)
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
! (different cases for metals and insulators: constant value for latter)

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: Te
    Real ( Kind = wp )                :: alp

    If (isMetal) Then
      alp = Ke(Te) / Ce(Te)
    Else
      alp = Diff0
    End If

  End Function alp

  Subroutine calcchies(chi_ep)

! Calculate electron-phonon coupling friction term (chi_ep)
! for homogeneously coupled system: uses mean electronic temperature
! and interpolates from tabulated values given in g.dat file

    Implicit None

    Real ( Kind = wp ), Intent ( Inout ) :: chi_ep
    Real ( Kind = wp )                   :: eltempav = 0.0_wp
    Real ( Kind = wp )                   :: epc = 0.0_wp

    Call eltemp_mean (eltempav)

    Call interpolate (gel, gtable, eltempav, epc)

    chi_ep = epc
				
  End Subroutine calcchies

  Subroutine boundaryHalo ()

! fills halo regions of electronic temperature lattice from neighbouring sections
! (periodic boundary conditions)

    Implicit None
    Integer :: i,ii,iii1,iii2,j,jj,jjj1,jjj2,k,kk,kkk1,kkk2,ijk1,ijk2
    Integer, Dimension(4) :: req
    Integer, Allocatable :: stats(:,:)

    Allocate (stats(1:MPI_STATUS_SIZE,1:4))
    If (mxnode>1) Then

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

  Subroutine boundaryCond (key, temp)

! appends halo regions of entire electronic temperature lattice with appropriate boundary conditions

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: temp
    Integer,            Intent ( In ) :: key

    Integer :: i,ii,j,jj,k,kk,ijk1,ijk2
    Integer, Dimension(4) :: req
    Integer, Dimension(MPI_STATUS_SIZE,4) :: stat


    Select Case (key)
  ! Periodic boundary conditions
    Case (1)
      If (ttmbcmap(1)>=0 .or. ttmbcmap(2)>=0) Then
        If (mxnode>1) Then
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
        If (mxnode>1) Then
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
        If (mxnode>1) Then
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

  Subroutine peakProfiler(lat, peakfile, nstep, freq)

! prints a slice across y (centred in xz-plane) of a lattice to a file

    Implicit None

    Real ( Kind = wp ), Dimension (numcell), Intent( In ) :: lat
    Real ( Kind = wp ), Dimension (ntsys(2)) :: laty
    Integer, Intent( In ) :: nstep,freq
    Character ( Len = * ), Intent( In ) :: peakfile
    Integer :: iounit = 114
    Integer :: i,j,k,l,ijk

    laty = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0 .or. nstep==1) Then
        i = midI(1) - ntcelloff(1)
        k = midI(3) - ntcelloff(3)
        If (i>0 .and. i<=ntcell(1) .and. k>0 .and. k<=ntcell(3)) Then
          Do j=1,ntcell(2)
            l = j + ntcelloff(2)
            ijk = 1 + i + (ntcell(1)+2) * (j + k * (ntcell(2)+2))
            laty(l) = lat(ijk)
          End Do
        End If
        If (mxnode>1) Call gsum (laty)
        If (idnode==0) Then
          If (nstep==1) Then
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
 
  Subroutine peakProfilerElec(peakfile, nstep, freq)

! prints a slice across y (centred in xz-plane) of the electronic temperature lattice to a file

    Implicit None

    Real ( Kind = wp ), Dimension (eltsys(2)) :: laty
    Integer, Intent( In ) :: nstep,freq
    Character ( Len = * ), Intent( In ) :: peakfile
    Integer :: iounit = 114
    Integer :: i,j,jj,jmin,jmax,k,l,ijk

    laty = 0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0 .or. nstep==1) Then
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
        If (mxnode>1) Call gsum (laty)
        If (idnode==0) Then
          If (nstep==1) Then
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

  Subroutine printLatticeStatsToFile(lat, latfile, time, nstep, freq)

! prints lattice statistics (minimum, maximum, sum) to file

    Implicit None

    Real ( Kind = wp ), Dimension(numcell), Intent ( In ) :: lat
    Integer, Intent( In ) :: nstep,freq
    Real ( Kind = wp ), Intent ( In ) :: time
    Character ( Len = * ), Intent ( In ) :: latfile

    Real (kind = wp) :: lat_sum,lat_min,lat_max,rtotal
    Integer :: iounit = 113
    Integer :: i,j,k,ijk

    lat_sum=0.0_wp
    lat_min=1.0e30_wp
    lat_max=0.0_wp

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0 .or. nstep==1) Then

        Do k = 1,ntcell(3)
          Do j = 1,ntcell(2)
            Do i = 1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
              lat_sum = lat_sum + lat(ijk)
              lat_min = Min(lat_min, lat(ijk))
              lat_max = Max(lat_max, lat(ijk))
            End Do
          End Do
        End Do

        If (mxnode>1) Then
          Call gsum(lat_sum)
          Call gmin(lat_min)
          Call gmax(lat_max)
        End If

        rtotal = 1.0_wp/Real(acell,Kind=wp)

        If (idnode == 0) Then
          If (nstep==1) Then
            Open(Unit=iounit, File=latfile, Status='replace')
          Else
            Open(Unit=iounit, File=latfile, Position='append')
          End If
          Write(iounit,'(i8,1p,es12.4,1p,es12.4,1p,es12.4,1p,es12.4,1p,es12.4)') &
                        nstep, time, lat_min, lat_max, lat_sum*rtotal, lat_sum
          Close(iounit)
        End If

      End If
    End If
      
  End Subroutine printLatticeStatsToFile
  
  Subroutine printElecLatticeStatsToFile(latfile, time, temp0, nstep, freq)

! prints electronic temperature lattice statistics (minimum, maximum, sum) 
! and energy (E = integral of Ce(Te)*Te between temp0 and Te) to file

    Implicit None

    Integer, Intent( In ) :: nstep,freq
    Real ( Kind = wp ), Intent ( In ) :: time,temp0
    Character ( Len = * ), Intent ( In ) :: latfile

    Real ( Kind = wp ) :: lat_sum,lat_min,lat_max,Ue,tmp,sgnplus
    Integer :: iounit = 115
    Integer :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax
    Integer :: ijk,numint,n,lx,ly,lz

    If (freq /= 0) Then
      If (Mod(nstep,freq)==0 .or. nstep==1) Then

        Call eltemp_sum(lat_sum)
        Call eltemp_min(lat_min)
        Call eltemp_max(lat_max)

        Ue = 0.0_wp

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
                      numint = Merge(Floor(eltemp(ijk,ii,jj,kk)-temp0),0,&
                                     (act_ele_cell(ijk,0,0,0)>zero_plus .or. (ii/=0 .or. jj/=0 .or. kk/=0)))
                      tmp = 0.0_wp
                      sgnplus = Sign (1.0_wp, Real(numint, Kind=wp))
                      Do n = 1,Abs(numint)
                        tmp = tmp + 0.5_wp * sgnplus * (Ce(temp0+sgnplus*Real(n-1, Kind=wp))+Ce(temp0+sgnplus*Real(n, Kind=wp)))
                      End Do
                      tmp = tmp + 0.5_wp * (Ce(temp0+Real(numint, Kind=wp))+Ce(eltemp(ijk,ii,jj,kk)))
                      Ue = Ue + tmp*volume*kB_to_eV
                    End If
                  End Do
                End Do
              End Do

            End Do
          End Do
        End Do

        If (mxnode>1) Call gsum(Ue)

        If (idnode == 0) Then
          If (nstep==1) Then
            Open(Unit=iounit, File=latfile, Status='replace')
          Else
            Open(Unit=iounit, File=latfile, Position='append')
          End If
          Write(iounit,'(i8,1p,es12.4,1p,es12.4,1p,es12.4,1p,es12.4,1p,es12.4)') nstep, time, lat_min, lat_max, lat_sum, Ue
          Close(iounit)
        End If

      End If
    End If
      
  End Subroutine printElecLatticeStatsToFile
  
  Subroutine eltemp_sum (eltempsum)

    Implicit None

    Real ( Kind = wp ), Intent ( Out ) :: eltempsum
    Real ( Kind = wp )                 :: tmp
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempsum = 0.0_wp

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
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
            Do j = jmin, jmax
              ly = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
              Do i = imin, imax
                lx = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                lrange = (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3))
                ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                tmp = eltemp(ijk,ii,jj,kk) * Merge (act_ele_cell (ijk,0,0,0), 1.0_wp, lcentre)
                If (lrange) eltempav = eltempav + tmp
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    If (mxnode>1) Call gsum (eltempsum)

  End Subroutine eltemp_sum

  Subroutine eltemp_mean (eltempav)

    Implicit None

    Real ( Kind = wp ), Intent ( Out ) :: eltempav
    Real ( Kind = wp )                 :: tmp
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz,acl
    Logical                            :: lrange,lcentre

    eltempav = 0.0_wp

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
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
            Do j = jmin, jmax
              ly = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
              Do i = imin, imax
                lx = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                lrange = (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3))
                ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                tmp = eltemp(ijk,ii,jj,kk) * Merge (act_ele_cell (ijk,0,0,0), 1.0_wp, lcentre)
                If (lrange) eltempav = eltempav + tmp
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    If (mxnode>1) Call gsum (eltempav)

    acl = eltsys(1)*eltsys(2)*eltsys(3) + acell - ntsys(1)*ntsys(2)*ntsys(3)
    eltempav = eltempav / Real (acl, Kind = wp)

  End Subroutine eltemp_mean

  Subroutine eltemp_maxKe (eltempmax)

    Implicit None

    Real ( Kind = wp ), Intent ( Out ) :: eltempmax
    Real ( Kind = wp )                 :: eltempKe
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempmax = 0.0_wp

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
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
            Do j = jmin, jmax
              ly = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
              Do i = imin, imax
                lx = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                lrange = (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3))
                If (lcentre) lrange = (lrange .and. (act_ele_cell(ijk,0,0,0)>zero_plus))
                eltempKe = (tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                If (lrange) eltempmax = Max (eltempmax, eltempKe)
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    If (mxnode>1) Call gmax (eltempmax)

  End Subroutine eltemp_maxKe
	
  Subroutine eltemp_max (eltempmax)

    Implicit None

    Real ( Kind = wp ), Intent ( Out ) :: eltempmax
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempmax = 0.0_wp

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
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
            Do j = jmin, jmax
              ly = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
              Do i = imin, imax
                lx = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                lrange = (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3))
                If (lcentre) lrange = (lrange .and. (act_ele_cell(ijk,0,0,0)>zero_plus))
                If (lrange) eltempmax = Max (eltempmax, eltemp(ijk,ii,jj,kk))
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    If (mxnode>1) Call gmax (eltempmax)

  End Subroutine eltemp_max

  Subroutine eltemp_minKe (eltempmin)

    Implicit None

    Real ( Kind = wp ), Intent ( Out ) :: eltempmin
    Real ( Kind = wp )                 :: eltempKe
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange

    eltempmin = 1.0e30_wp

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
                lrange = (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3))
                ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                eltempKe = (tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                If (lrange) eltempmin = Min (eltempmin, eltempKe)
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    If (mxnode>1) Call gmin (eltempmin)

  End Subroutine eltemp_minKe
	
	
  Subroutine eltemp_min (eltempmin)

    Implicit None

    Real ( Kind = wp ), Intent ( Out ) :: eltempmin
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange

    eltempmin = 1.0e30_wp

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
                lrange = (lx>0 .and. lx<=eltsys(1) .and. ly>0 .and. ly<=eltsys(2) .and. lz>0 .and. lz<=eltsys(3))
                If (ii==0 .and. jj==0 .and. kk==0) lrange
                If (lrange) eltempmin = Min (eltempmin, eltemp(ijk,ii,jj,kk))
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    If (mxnode>1) Call gmin (eltempmin)

  End Subroutine eltemp_min
	
	
  Subroutine interpolate(tabsize, tab, x0, resp)

! Interpolation helper subroutine

    Implicit None

    Integer,                                      Intent( In )  :: tabsize
    Real( Kind = wp ), Dimension (1:tabsize,1:2), Intent( In )  :: tab
    Real( Kind = wp ),                            Intent( In )  :: x0
    Real( Kind = wp ),                            Intent( Out ) :: resp

    Integer :: i,j,k

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

  Subroutine redistribute_Te (temp0)

! redistribute electronic energy when electronic temperature voxels are closed

    Implicit None

    Real( Kind = wp ), Intent ( In ) :: temp0
    Integer :: i, j, k, ii, jj, kk, ijk, ijk1, ijk2, ijkpx, ijkmx, ijkpy, ijkmy, ijkpz, ijkmz, n, numint
    Real( Kind = wp ) :: energy_per_cell, U_e, start_Te, end_Te, increase, oldCe, newCe
    Real( Kind = wp ) :: energy_diff, tmp2, act_sur_cells, sgnplus

    Integer, Dimension(4) :: req
    Integer, Dimension(MPI_STATUS_SIZE,4) :: stat
    Real( Kind = wp ), Allocatable :: energydist (:,:,:,:), buffer (:,:,:,:)

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
		! Calculate amount of energy left in the cell (Ue = integral of Ce(Te)d(Te) between 300K and Te)
            tmp2 = 0.0_wp
            numint = Floor(eltemp(ijk,0,0,0) - temp0)
            sgnplus = Sign (1.0_wp, Real(numint, Kind=wp))
            Do n = 1,Abs(numint)
              tmp2 = tmp2 + 0.5_wp * sgnplus * (Ce(temp0+sgnplus*Real(n-1, Kind=wp))+Ce(temp0+sgnplus*Real(n, Kind=wp)))
            End Do
            tmp2 = tmp2 + 0.5_wp * (Ce(temp0+Real(numint, Kind=wp))+Ce(eltemp(ijk,ii,jj,kk)))
            U_e = tmp2*volume*kB_to_eV
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

    If(mxnode>1) Then

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
    ! redistributed energies and identify relevant cells

    Do kk = -1, 1
      Do jj = -1, 1
        Do ii = -1, 1
          Do ijk = 1, numcell
            If (Abs (energydist (ijk,ii,jj,kk))>zero_plus) Then
              energy_per_cell = energydist (ijk,ii,jj,kk)
              sgnplus = Sign (1.0_wp, energy_per_cell)
              start_Te = eltemp(ijk,ii,jj,kk)
              energy_diff = sgnplus*energy_per_cell
              oldCe = Ce(start_Te)
              ! increase/decrease temperature of electronic cell by 1 kelvin until
              ! required energy change is reached or exceeded
              Do While (energy_diff>=zero_plus)
                newCe = Ce(start_Te+sgnplus)
                increase = 0.5_wp * (oldCe+newCe)
                energy_diff = energy_diff - increase*volume*kB_to_eV
                start_Te = start_Te + sgnplus
                oldCe = newCe
              End Do
              ! now interpolate for new temperature by using average specific heat
              ! over surrounding one-kelvin interval
              energy_diff = energy_diff + increase*volume*kB_to_eV
              start_Te = start_Te - sgnplus
              newCe = Ce(start_Te)
              end_Te = start_Te + 2.0_wp * sgnplus * energy_diff / (oldCe+newCe)
              eltemp_adj(ijk,ii,jj,kk) = end_Te
              adjust(ijk,ii,jj,kk) = .true.
            End If
          End Do
        End Do
      End Do
    End Do

    Deallocate (energydist)

  End Subroutine redistribute_Te


End Module ttm_utils
