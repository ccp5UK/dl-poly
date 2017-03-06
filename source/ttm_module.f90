Module ttm_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining arrays and initial parameters for 
! two-temperature model (ttm)
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton may 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use comms_module
  Use config_module, only : cell
  Use domains_module

  Implicit None

  Real( Kind = wp ), Allocatable :: eltemp(:,:,:,:),eltemp_adj(:,:,:,:)
  Real( Kind = wp ), Allocatable :: asource(:),tempion(:),gsource(:)
  Real( Kind = wp ), Allocatable :: act_ele_cell(:,:,:,:),old_ele_cell(:,:,:,:)
  Logical          , Allocatable :: adjust(:,:,:,:)

  Real( Kind = wp ), Allocatable :: cetable(:,:),gtable(:,:),ketable(:,:)

  Integer :: ntsys(3),eltsys(3)
  Integer :: ntcell(3),eltcell(3)
  Integer :: ntcelloff(3),midI(3),midE(3),zeroE(3)

  Integer :: tmpmsgx,tmpmsgy,tmpmsgz
  Integer :: nummsgx,nummsgy,nummsgz

  Real ( Kind = wp ) :: delx,dely,delz,volume
  Real ( Kind = wp ) :: zerocell(3)
  Integer :: numcell
  Integer :: ttmbc(6),ttmbcmap(6)

  Logical :: l_ttm,isMetal,l_epcp,deactivation
  Integer :: CeType,bcTypeE,ttmstats,ttmtraj,tdepoType,sdepoType,KeType,gvar
  Real ( Kind = wp ) :: fluxout,ttmoffset
  Real ( Kind = wp ) :: sh_A,sh_B,Ka0,Ce0
  Real ( Kind = wp ) :: Cemax,Tfermi,Diff0
  Real ( Kind = wp ) :: dEdX,sig,sigmax,tdepo,tcdepo

! DEBUG (TODO)
  Real ( Kind = wp ) :: epstart
  Integer :: keyres0,nstepcpl = 0

  Integer :: cel,gel,kel
  Integer :: acell,acell_old,amin

  Real ( Kind = wp ), Save :: Jm3K_to_kBA3,JKms_to_kBAps,kB_to_eV,mJcm2_to_eVA2
  Real ( Kind = wp ), Save :: cellrho,epc_to_chi
  Real ( Kind = wp ) :: fluence,pdepth
  Real ( Kind = wp ) :: epthreshold = 1.1_wp

  Public :: allocate_ttm_arrays , deallocate_ttm_arrays

Contains

  Subroutine allocate_ttm_arrays(temp,megatm)

    Implicit None

    Real ( Kind = wp ), Intent ( In ) :: temp
    Integer, Intent ( In ) :: megatm
    Real ( Kind = wp ) :: start, finish
    Integer, Dimension ( 1:7 ) :: fail
    Integer :: i,numbc,numbcmap
    Integer :: basicslice,oneslicex,oneslicey,oneslicez
    Integer :: bbasicslice,boneslicex,boneslicey,boneslicez
    Integer ( Kind = MPI_ADDRESS_KIND ) :: dbleth,inleth,lb1,lb2
    Integer ( Kind = MPI_ADDRESS_KIND ) :: xlth,ylth,bxlth,bylth

    fail = 0

! Setup constants based on fundamental values (found in
! setup_module.f90)

    JKms_to_kBAps = 10.0_wp/(boltz*tenunt)    ! convert W m^-1 K^-1 to kB A^-1 ps^-1
    Jm3K_to_kBA3  = 1.0e-7_wp/(boltz*tenunt)  ! convert J m^-3 K^-1 to kB A^-3
    kB_to_eV      = boltz/eu_ev               ! convert kB to eV
    mJcm2_to_eVA2 = 1.0e4_wp/(eu_ev*tenunt)   ! convert mJ cm^-2 to eV A^-2

! Calculate atomic density and conversion factor to calculate
! electron-phonon friction factor (chi_ep)

    cellrho    = Real(megatm,Kind=wp)/(cell(1)*cell(5)*cell(9))
    epc_to_chi = 1.0e-12_wp*Jm3K_to_kBa3/(3.0_wp*cellrho)

! Convert inputs for ion temperature grid (number in z-direction)
! into numbers for x- and y-directions, grid spacings and cell volume

    delz = cell(9)/Real(ntsys(3),wp)
    ntsys(1) = Nint(cell(1)/delz)
    ntsys(2) = Nint(cell(5)/delz) 
    delx = cell(1)/Real(ntsys(1),wp)
    dely = cell(5)/Real(ntsys(2),wp)
    volume = delx*dely*delz

! Check number of electronic temperature cells is greater than/
! equal to number of ionic temperature cells

    If (Any(eltsys<ntsys)) Call error(670)

! If slab option selected, check for sufficient electronic temperature
! cells to redistribute energy when ionic tmeperature cells are switched off:
! if not available, switch off this option

    If (deactivation .and. (eltsys(1)<ntsys(1)+2 .or. eltsys(2)<ntsys(2)+2 .or. eltsys(3)<ntsys(3)+2)) Then
      Call warning(500,0.0_wp,0.0_wp,0.0_wp)
      deactivation = .false.
    End If

! Print sizes of ionic and electronic temperature cells and grids,
! and average number of particles per ionic temperature cell

    If (idnode == 0) Then
      Write(nrite,'(/,1x,a,3(2x,f8.4))') "temperature cell size (A)        (x,y,z): ",delx,dely,delz
      Write(nrite,'(1x,a,3(2x,i8))')     "ionic temperature grid size      (x,y,z): ",ntsys(1),ntsys(2),ntsys(3)
      Write(nrite,'(1x,a,3(2x,i8))')     "electronic temperature grid size (x,y,z): ",eltsys(1),eltsys(2),eltsys(3)
      Write(nrite,'(1x,a,f10.4)')        "average no. of atoms per cell           : ",cellrho*volume
    End If

! Check sufficient parameters are specified for electronic specific
! heats, thermal conductivity, diffusivity, energy loss and laser deposition
! and rescale to calculate efficiently (in functions Ce, Ke, KeD etc.)

    Select Case (CeType)
    Case (0)
    ! constant electronic specific heat: converted from kB/atom to kB/A^3
    ! by multiplication of atomic density
      Ce0 = Ce0*cellrho
    Case (1)
    ! hyperbolic tangent electronic specific heat: multiplier converted 
    ! from kB/atom to kB/A^3, temperature term (K^-1) scaled by 10^-4
      If (Abs(sh_A) <= zero_plus .or. Abs(sh_B) <= zero_plus) Call error(671)
      sh_A = sh_A*cellrho
      sh_B = sh_B*1.0e-4_wp
    Case (2)
    ! linear electronic specific heat to Fermi temperature: maximum
    ! value converted from kB/atom to kB/A^3
      If (Abs(Tfermi) <= zero_plus .or. Abs(Cemax) <= zero_plus) Call error(671)
      Cemax = Cemax*cellrho
    End Select

    Select Case (KeType)
    ! constant and Drude thermal conductivity: converted from W m^-1 K^-1
    ! to kB ps^-1 A^-1
    Case (1,2)
      If (isMetal .and. Abs(Ka0) <= zero_plus) Call error(672)
      Ka0 = Ka0*JKms_to_kBAps
    End Select

    ! thermal diffusivity: converted from m^2 s^-1 to A^2 ps^-1
    If (.not. isMetal .and. Abs(Diff0) <= zero_plus) Call error(673)
    Diff0 = Diff0*1.0e8_wp

    ! spatial deposition (gaussian) standard deviation: convert from nm to A
    sig = sig*10.0_wp

    ! penetration depth: convert from nm to A
    If (sdepoType == 2 .and. (Abs(dEdX) <= zero_plus .or. Abs(pdepth-1.0_wp) <= zero_plus)) &
    Call warning(510,0.0_wp,0.0_wp,0.0_wp)
    pdepth = 10.0_wp*pdepth

    ! fluence: convert from mJ cm^-2 to eV A^-2
    fluence = fluence*mJcm2_to_eVA2

    ! electronic stopping power: convert from eV/nm to eV/A
    If (Abs(dEdx) <= zero_plus) Call warning(515,0.0_wp,0.0_wp,0.0_wp)
    dEdX = 0.1_wp*dEdX

! Determine number of ion temperature cells for domain and
! offsets for ion temperature determination

    start = cell(1)*Real(idx,wp)*r_nprx
    finish = cell(1)*Real(idx+1,wp)*r_nprx
    ntcell(1) = Ceiling(finish/delx) - Ceiling(start/delx)
    ntcelloff(1) = Ceiling(start/delx)
    zerocell(1) = 0.5_wp*cell(1) - delx*Real(Ceiling(start/delx),wp)

    start = cell(5)*Real(idy,wp)*r_npry
    finish = cell(5)*Real(idy+1,wp)*r_npry
    ntcell(2) = Ceiling(finish/dely) - Ceiling(start/dely)
    ntcelloff(2) = Ceiling(start/dely)
    zerocell(2) = 0.5_wp*cell(5) - dely*Real(Ceiling(start/dely),wp)

    start = cell(9)*Real(idz,wp)*r_nprz
    finish = cell(9)*Real(idz+1,wp)*r_nprz
    ntcell(3) = Ceiling(finish/delz) - Ceiling(start/delz)
    ntcelloff(3) = Ceiling(start/delz)
    zerocell(3) = 0.5_wp*cell(9) - delz*Real(Ceiling(start/delz),wp)

    numcell = (ntcell(1)+2)*(ntcell(2)+2)*(ntcell(3)+2)

! Determine mid-values for ion and electronic temperature grid

    midI(:) = INT((ntsys(:)+1)/2)
    midE(:) = INT((eltsys(:)+1)/2)

! Determine number of multiple ion temperature grids for electronic
! temperature grids

    eltcell(1) = Ceiling(Real(eltsys(1)-ntsys(1), Kind=wp)/Real(2*ntsys(1), Kind=wp))
    eltcell(2) = Ceiling(Real(eltsys(2)-ntsys(2), Kind=wp)/Real(2*ntsys(2), Kind=wp))
    eltcell(3) = Ceiling(Real(eltsys(3)-ntsys(3), Kind=wp)/Real(2*ntsys(3), Kind=wp))


! Determine positions of boundaries for electronic temperature grids

!   -x boundary
    numbc = -(eltsys(1)-ntsys(1))/2
    numbc = MOD(numbc+ntsys(1)*(eltcell(1)+1),ntsys(1)) + 1
    zeroE(1) = numbc - 1
    numbcmap = (eltsys(1)-ntsys(1))/2
    numbcmap = MOD(numbcmap+ntsys(1)*(eltcell(1)+1)-1,ntsys(1)) + 1
    If (numbc>ntcelloff(1) .and. numbc<=(ntcelloff(1)+ntcell(1))) Then
      ttmbc(1) = numbc - ntcelloff(1)
      Do i = 0, nprx-1
        start = cell(1)*Real(i,wp)*r_nprx
        finish = cell(1)*Real(i+1,wp)*r_nprx
        If (numbcmap>Ceiling(start/delx) .and. numbcmap<=Ceiling(finish/delx)) ttmbcmap(1) = idcube(i,idy,idz)
      End Do
    Else
      ttmbc(1) = 0
      ttmbcmap(1) = -1
    End If

!   +x boundary
    numbc = (eltsys(1)-ntsys(1))/2
    numbc = MOD(numbc+ntsys(1)*(eltcell(1)+1)-1,ntsys(1)) + 1
    numbcmap = -(eltsys(1)-ntsys(1))/2
    numbcmap = MOD(numbcmap+ntsys(1)*(eltcell(1)+1),ntsys(1)) + 1
    If (numbc>ntcelloff(1) .and. numbc<=(ntcelloff(1)+ntcell(1))) Then
      ttmbc(2) = numbc - ntcelloff(1)
      Do i = 0, nprx-1
        start = cell(1)*Real(i,wp)*r_nprx
        finish = cell(1)*Real(i+1,wp)*r_nprx
        If (numbcmap>Ceiling(start/delx) .and. numbcmap<=Ceiling(finish/delx)) ttmbcmap(2) = idcube(i,idy,idz)
      End Do
    Else
      ttmbc(2) = 0
      ttmbcmap(2) = -1
    End If
    
!   -y boundary
    numbc = -(eltsys(2)-ntsys(2))/2
    numbc = MOD(numbc+ntsys(2)*(eltcell(2)+1),ntsys(2)) + 1
    zeroE(2) = numbc - 1
    numbcmap = (eltsys(2)-ntsys(2))/2
    numbcmap = MOD(numbcmap+ntsys(2)*(eltcell(2)+1)-1,ntsys(2)) + 1
    If (numbc>ntcelloff(2) .and. numbc<=(ntcelloff(2)+ntcell(2))) Then
      ttmbc(3) = numbc - ntcelloff(2)
      Do i = 0, npry-1
        start = cell(5)*Real(i,wp)*r_npry
        finish = cell(5)*Real(i+1,wp)*r_npry
        If (numbcmap>Ceiling(start/dely) .and. numbcmap<=Ceiling(finish/dely)) ttmbcmap(3) = idcube(idx,i,idz)
      End Do
    Else
      ttmbc(3) = 0
      ttmbcmap(3) = -1
    End If

!   +y boundary
    numbc = (eltsys(2)-ntsys(2))/2
    numbc = MOD(numbc+ntsys(2)*(eltcell(2)+1)-1,ntsys(2)) + 1
    numbcmap = -(eltsys(2)-ntsys(2))/2
    numbcmap = MOD(numbcmap+ntsys(2)*(eltcell(2)+1),ntsys(2)) + 1
    If (numbc>ntcelloff(2) .and. numbc<=(ntcelloff(2)+ntcell(2))) Then
      ttmbc(4) = numbc - ntcelloff(2)
      Do i = 0, npry-1
        start = cell(5)*Real(i,wp)*r_npry
        finish = cell(5)*Real(i+1,wp)*r_npry
        If (numbcmap>Ceiling(start/dely) .and. numbcmap<=Ceiling(finish/dely)) ttmbcmap(4) = idcube(idx,i,idz)
      End Do
    Else
      ttmbc(4) = 0
      ttmbcmap(4) = -1
    End If

!   -z boundary
    numbc = -(eltsys(3)-ntsys(3))/2
    numbc = MOD(numbc+ntsys(3)*(eltcell(3)+1),ntsys(3)) + 1
    zeroE(3) = numbc - 1
    numbcmap = (eltsys(3)-ntsys(3))/2
    numbcmap = MOD(numbcmap+ntsys(3)*(eltcell(3)+1)-1,ntsys(3)) + 1
    If (numbc>ntcelloff(3) .and. numbc<=(ntcelloff(3)+ntcell(3))) Then
      ttmbc(5) = numbc - ntcelloff(3)
      Do i = 0, nprz-1
        start = cell(9)*Real(i,wp)*r_nprz
        finish = cell(9)*Real(i+1,wp)*r_nprz
        If (numbcmap>Ceiling(start/delz) .and. numbcmap<=Ceiling(finish/delz)) ttmbcmap(5) = idcube(idx,idy,i)
      End Do
    Else
      ttmbc(5) = 0
      ttmbcmap(5) = -1
    End If

!   +z boundary
    numbc = (eltsys(3)-ntsys(3))/2
    numbc = MOD(numbc+ntsys(3)*(eltcell(3)+1)-1,ntsys(3)) + 1
    numbcmap = -(eltsys(3)-ntsys(3))/2
    numbcmap = MOD(numbcmap+ntsys(3)*(eltcell(3)+1),ntsys(3)) + 1
    If (numbc>ntcelloff(3) .and. numbc<=(ntcelloff(3)+ntcell(3))) Then
      ttmbc(6) = numbc - ntcelloff(3)
      Do i = 0, nprz-1
        start = cell(9)*Real(i,wp)*r_nprz
        finish = cell(9)*Real(i+1,wp)*r_nprz
        If (numbcmap>Ceiling(start/delz) .and. numbcmap<=Ceiling(finish/delz)) ttmbcmap(6) = idcube(idx,idy,i)
      End Do
    Else
      ttmbc(6) = 0
      ttmbcmap(6) = -1
    End If

! Derived MPI datatypes for communication of temperatures (MPI 2.x+)

    If (mxnode>1) Then
      Call MPI_TYPE_GET_EXTENT (wp_mpi, lb1, dbleth, ierr)
      Call MPI_TYPE_GET_EXTENT (MPI_INTEGER, lb2, inleth, ierr)
      Call MPI_TYPE_VECTOR (1, 1, 1, wp_mpi, basicslice, ierr)
      Call MPI_TYPE_VECTOR (1, 2, 2, MPI_INTEGER, bbasicslice, ierr)
      xlth = (ntcell(1) + 2) * dbleth
      ylth = (ntcell(2) + 2) * xlth
      Call MPI_TYPE_CREATE_HVECTOR (ntcell(2)    , 1, xlth, basicslice, oneslicex, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (ntcell(3)    , 1, ylth, oneslicex, tmpmsgx, ierr)
      Call MPI_TYPE_COMMIT (tmpmsgx, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (1, ntcell(1) + 2, xlth, basicslice, oneslicey, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (ntcell(3)    , 1, ylth, oneslicey, tmpmsgy, ierr)
      Call MPI_TYPE_COMMIT (tmpmsgy, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (1, ntcell(1) + 2, xlth, basicslice, oneslicez, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (1, ntcell(2) + 2, ylth, oneslicez, tmpmsgz, ierr)
      Call MPI_TYPE_COMMIT (tmpmsgz, ierr)
      bxlth = 2 * (ntcell(1) + 2) * inleth
      bylth = (ntcell(2) + 2) * bxlth
      Call MPI_TYPE_CREATE_HVECTOR (ntcell(2)    , 1, bxlth, bbasicslice, boneslicex, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (ntcell(3)    , 1, bylth, boneslicex, nummsgx, ierr)
      Call MPI_TYPE_COMMIT (nummsgx, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (1, ntcell(1) + 2, bxlth, bbasicslice, boneslicey, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (ntcell(3)    , 1, bylth, boneslicey, nummsgy, ierr)
      Call MPI_TYPE_COMMIT (nummsgy, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (1, ntcell(1) + 2, bxlth, bbasicslice, boneslicez, ierr)
      Call MPI_TYPE_CREATE_HVECTOR (1, ntcell(2) + 2, bylth, boneslicez, nummsgz, ierr)
      Call MPI_TYPE_COMMIT (nummsgz, ierr)
    Else
      tmpmsgx = 0; tmpmsgy = 0; tmpmsgz = 0
      nummsgx = 0; nummsgy = 0; nummsgz = 0
    End If

! Array allocation and initialization

    Allocate (eltemp(1:numcell,-eltcell(1):eltcell(1),-eltcell(2):eltcell(2),-eltcell(3):eltcell(3))    , Stat = fail(1))
    Allocate (asource(1:numcell)                                                                        , Stat = fail(2))
    Allocate (tempion(1:numcell)                                                                        , Stat = fail(3))
    Allocate (gsource(1:numcell)                                                                        , Stat = fail(4))
    Allocate (eltemp_adj(1:numcell,-eltcell(1):eltcell(1),-eltcell(2):eltcell(2),-eltcell(3):eltcell(3)), Stat = fail(5))
    Allocate (act_ele_cell(1:numcell,-1:1,-1:1,-1:1), old_ele_cell(1:numcell,-1:1,-1:1,-1:1)            , Stat = fail(6))
    Allocate (adjust(1:numcell,-1:1,-1:1,-1:1)                                                          , Stat = fail(7))

    If (Any(fail > 0)) Call error(1083)

    eltemp(:,:,:,:) = temp
    eltemp_adj(:,:,:,:) = 0.0_wp
    gsource(:) = 0.0_wp
    asource(:) = 0.0_wp
    tempion(:) = 0.0_wp
    act_ele_cell(:,:,:,:) = 1.0_wp
    old_ele_cell(:,:,:,:) = 1.0_wp
    acell = ntsys(1)*ntsys(2)*ntsys(3)
    acell_old = acell
    adjust = .false.
    cel = 0
    gel = 0
    kel = 0

    keyres0 = 1

  End Subroutine allocate_ttm_arrays

  Subroutine deallocate_ttm_arrays()

    Implicit None

    Integer, Dimension ( 1:4 ) :: fail

    fail = 0

    Deallocate (eltemp,eltemp_adj,asource,tempion,gsource,act_ele_cell,old_ele_cell,adjust, Stat = fail(1))
    If (kel>0) Deallocate(ketable,                                                          Stat = fail(2))
    If (cel>0) Deallocate(cetable,                                                          Stat = fail(3))
    If (gel>0) Deallocate(gtable,                                                           Stat = fail(4))

    If (Any(fail > 0)) Call error(1084)

  End Subroutine deallocate_ttm_arrays

  Function idcube(i,j,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 hypercube mapping function
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2006
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Integer                  :: idcube

    idcube = i + nprx * ( j + npry * k )

  End Function idcube

End Module ttm_module
