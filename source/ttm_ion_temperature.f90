Subroutine ttm_ion_temperature(chi_ep,chi_es,vel_es2)

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

  Use setup_module
  Use comms_module
  Use config_module
  Use ttm_module
  Use ttm_track_module

  Implicit None

  Real ( Kind = wp ), Intent ( In ) :: chi_ep,chi_es,vel_es2
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

  If (mxnode>1) Then
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
    Call MPI_ISEND (ttmvom(ijk1,1), 1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (ttmvom(ijk1,2), 1, tmpmsgz, map(5), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgz, map(6), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (ttmvom(ijk1,3), 1, tmpmsgz, map(5), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgz, map(6), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (ttmvom(ijk1,4), 1, tmpmsgz, map(5), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, tmpmsgz, map(6), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
    Call MPI_ISEND (ttmvom(ijk2,1), 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (ttmvom(ijk2,2), 1, tmpmsgz, map(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, tmpmsgz, map(5), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (ttmvom(ijk2,3), 1, tmpmsgz, map(6), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, tmpmsgz, map(5), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (ttmvom(ijk2,4), 1, tmpmsgz, map(6), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, tmpmsgz, map(5), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
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
    Call MPI_ISEND (ttmvom(ijk1,1), 1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (ttmvom(ijk1,2), 1, tmpmsgy, map(3), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgy, map(4), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (ttmvom(ijk1,3), 1, tmpmsgy, map(3), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgy, map(4), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (ttmvom(ijk1,4), 1, tmpmsgy, map(3), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, tmpmsgy, map(4), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
    Call MPI_ISEND (ttmvom(ijk2,1), 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (ttmvom(ijk2,2), 1, tmpmsgy, map(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, tmpmsgy, map(3), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (ttmvom(ijk2,3), 1, tmpmsgy, map(4), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, tmpmsgy, map(3), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (ttmvom(ijk2,4), 1, tmpmsgy, map(4), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, tmpmsgy, map(3), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
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
    Call MPI_ISEND (ttmvom(ijk1,1), 1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (ttmvom(ijk1,2), 1, tmpmsgx, map(1), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgx, map(2), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (ttmvom(ijk1,3), 1, tmpmsgx, map(1), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgx, map(2), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (ttmvom(ijk1,4), 1, tmpmsgx, map(1), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf4(ijk2)    , 1, tmpmsgx, map(2), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
    Call MPI_ISEND (ttmvom(ijk2,1), 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk1)    , 1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (ttmvom(ijk2,2), 1, tmpmsgx, map(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk1)    , 1, tmpmsgx, map(1), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (ttmvom(ijk2,3), 1, tmpmsgx, map(2), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk1)    , 1, tmpmsgx, map(1), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (ttmvom(ijk2,4), 1, tmpmsgx, map(2), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf4(ijk1)    , 1, tmpmsgx, map(1), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
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

  If (mxnode>1) Then
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
    Call MPI_ISEND (tempion(ijk1) , 1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (gsource(ijk1) , 1, tmpmsgz, map(5), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgz, map(6), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (asource(ijk1) , 1, tmpmsgz, map(5), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgz, map(6), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, nummsgz, map(5), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, nummsgz, map(6), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
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
    Call MPI_ISEND (tempion(ijk1) , 1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (gsource(ijk1) , 1, tmpmsgy, map(3), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgy, map(4), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (asource(ijk1) , 1, tmpmsgy, map(3), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgy, map(4), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, nummsgy, map(3), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, nummsgy, map(4), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
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
    Call MPI_ISEND (tempion(ijk1) , 1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (buf1(ijk2)    , 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    Call MPI_ISEND (gsource(ijk1) , 1, tmpmsgx, map(1), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (buf2(ijk2)    , 1, tmpmsgx, map(2), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_ISEND (asource(ijk1) , 1, tmpmsgx, map(1), Grid3_tag, MPI_COMM_WORLD, req(5), ierr)
    Call MPI_IRECV (buf3(ijk2)    , 1, tmpmsgx, map(2), Grid3_tag, MPI_COMM_WORLD, req(6), ierr)
    Call MPI_ISEND (nat(2*ijk1-1) , 1, nummsgx, map(1), Grid4_tag, MPI_COMM_WORLD, req(7), ierr)
    Call MPI_IRECV (buf5(2*ijk2-1), 1, nummsgx, map(2), Grid4_tag, MPI_COMM_WORLD, req(8), ierr)
    Call MPI_WAITALL (8, req, stat, ierr)
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

  If (mxnode>1) Then
    Call gsum (acell)
    Call gsum (crho)
    ijk1 = 2 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+1) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ii = Merge (-1,0,(idx==nprx-1))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,ii,0,0), 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    ijk1 = 1 + (ntcell(1)) + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ii = Merge (1,0,(idx==0))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgx, map(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,ii,0,0), 1, tmpmsgx, map(1), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_WAITALL (4, req, stat, ierr)
    ijk1 = 1 + (ntcell(1)+2) * (1 + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2) + 1 + (ntcell(2)+2))
    jj = Merge (-1,0,(idy==npry-1))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,jj,0), 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2) + (ntcell(2)+2))
    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
    jj = Merge (1,0,(idy==0))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgy, map(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,jj,0), 1, tmpmsgy, map(3), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_WAITALL (4, req, stat, ierr)
    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
    ijk2 = 1 + (ntcell(1)+2) * ((ntcell(2)+2) * (ntcell(3) + 1))
    kk = Merge (-1,0,(idz==nprz-1))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,0,kk), 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
    ijk1 = 1 + (ntcell(1)+2) * ((ntcell(2)+2) * ntcell(3))
    ijk2 = 1
    kk = Merge (1,0,(idz==0))
    Call MPI_ISEND (act_ele_cell(ijk1,0,0,0),  1, tmpmsgz, map(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
    Call MPI_IRECV (act_ele_cell(ijk2,0,0,kk), 1, tmpmsgz, map(5), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
    Call MPI_WAITALL (4, req, stat, ierr)
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
!    Call MPI_ISEND (asource(ijk1), 1, tmpmsgx, map(2), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
!    Call MPI_IRECV (asource(ijk2), 1, tmpmsgx, map(1), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
!    Call MPI_WAITALL (2, req, stat, ierr)
!    ijk1 = 1 + (ntcell(1)+2) * (2*ntcell(2) + 2)
!    ijk2 = 1 + (ntcell(1)+2) * (ntcell(2)+2)
!    Call MPI_ISEND (asource(ijk1), 1, tmpmsgy, map(4), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
!    Call MPI_IRECV (asource(ijk2), 1, tmpmsgy, map(3), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
!    Call MPI_WAITALL (2, req, stat, ierr)
!    ijk1 = 1 + (ntcell(1)+2) * (ntcell(2)+2) * (ntcell(3))
!    ijk2 = 1
!    Call MPI_ISEND (asource(ijk1), 1, tmpmsgz, map(6), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
!    Call MPI_IRECV (asource(ijk2), 1, tmpmsgz, map(5), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
!    Call MPI_WAITALL (2, req, stat, ierr)
!  End If

  Deallocate (nat, ijkatm, Stat=fail)
  If (fail > 0) Call error(1086)

End Subroutine ttm_ion_temperature

