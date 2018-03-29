Module ttm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining arrays and initial parameters for 
! two-temperature model (ttm)
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton may 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use setup
  Use configuration,   Only : cell
  Use domains
  Use comms,           Only : wp_mpi,comms_type,gsum,gmin,gmax,gcheck,gsync, &
                              grid1_tag,grid2_tag
  Use parse,           Only : tabs_2_blanks, get_line, get_word, &
                              strip_blanks, word_2_real
  Use errors_warnings, Only : error,warning,info
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None

  Real( Kind = wp ), Allocatable :: eltemp(:,:,:,:),eltemp_adj(:,:,:,:)
  Real( Kind = wp ), Allocatable :: asource(:),tempion(:),ttmvom(:,:),gsource(:)
  Real( Kind = wp ), Allocatable :: act_ele_cell(:,:,:,:),old_ele_cell(:,:,:,:)
  Logical          , Allocatable :: adjust(:,:,:,:)

  Real( Kind = wp ), Allocatable :: cetable(:,:),gtable(:,:),detable(:,:),ketable(:,:)

  Integer :: ntsys(3),eltsys(3)
  Integer :: ntcell(3),eltcell(3)
  Integer :: ntcelloff(3),midI(3),midE(3),zeroE(3)

  Integer :: tmpmsgx,tmpmsgy,tmpmsgz
  Integer :: nummsgx,nummsgy,nummsgz

  Real ( Kind = wp ) :: delx,dely,delz,volume,rvolume
  Real ( Kind = wp ) :: zerocell(3)
  Integer :: numcell
  Integer :: ttmbc(6),ttmbcmap(6)

  Logical :: l_ttm,isMetal,l_epcp,redistribute,ttmthvel,ttmthvelz,oneway,ttmdyndens,findepo
  Integer :: CeType,KeType,DeType,gvar,bcTypeE,ttmstats,ttmtraj,tdepoType,sdepoType
  Real ( Kind = wp ) :: fluxout,ttmoffset,depostart,depoend
  Real ( Kind = wp ) :: sh_A,sh_B,Ka0,Ce0
  Real ( Kind = wp ) :: Cemax,Tfermi,Diff0
  Real ( Kind = wp ) :: dEdX,sig,sigmax,tdepo,tcdepo

! DEBUG (TODO)
  Real ( Kind = wp ) :: epstart
  Integer :: keyres0,nstepcpl = 0

  Integer :: cel,gel,del,kel
  Integer :: acell,acell_old,amin

  Real ( Kind = wp ), Save :: Jm3K_to_kBA3,JKms_to_kBAps,kB_to_eV,eV_to_kB,mJcm2_to_eVA2,epc_to_chi
  Real ( Kind = wp ) :: cellrho,rcellrho,sysrho
  Real ( Kind = wp ) :: fluence,pdepth
  Real ( Kind = wp ) :: epthreshold = 1.1_wp

  Real( Kind = wp ), Allocatable, Dimension (:,:,:) :: lat_U,lat_B,lat_I
  Real( Kind = wp ) :: norm

  Logical :: trackInit = .false.
  
  Public :: allocate_ttm_arrays , deallocate_ttm_arrays

Contains

  Subroutine allocate_ttm_arrays(comm)


    Type(comms_type), Intent(In) :: comm
    Real ( Kind = wp ) :: start, finish
    Integer, Dimension ( 1:7 ) :: fail
    Integer :: i,numbc,numbcmap
    Integer :: basicslice,oneslicex,oneslicey,oneslicez
    Integer :: bbasicslice,boneslicex,boneslicey,boneslicez
    Integer ( Kind = MPI_ADDRESS_KIND ) :: dbleth,inleth,lb1,lb2
    Integer ( Kind = MPI_ADDRESS_KIND ) :: xlth,ylth,bxlth,bylth
    Integer                             :: ierr

    fail = 0

    depostart = 0.0_wp
    depoend = 0.0_wp

! Setup constants based on fundamental values (found in
! setup.f90)

    JKms_to_kBAps = 10.0_wp/(boltz*tenunt)    ! convert W m^-1 K^-1 to kB A^-1 ps^-1
    Jm3K_to_kBA3  = 1.0e-7_wp/(boltz*tenunt)  ! convert J m^-3 K^-1 to kB A^-3
    kB_to_eV      = boltz/eu_ev               ! convert kB to eV
    eV_to_kB      = eu_ev/boltz               ! convert eV to kB
    mJcm2_to_eVA2 = 1.0e4_wp/(eu_ev*tenunt)   ! convert mJ cm^-2 to eV A^-2

    If (l_ttm) Then

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

      If (comm%mxnode>1) Then
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
      Allocate (asource(1:numcell),tempion(1:numcell),gsource(1:numcell)                                  , Stat = fail(2))
      Allocate (ttmvom(1:numcell,1:4)                                                                     , Stat = fail(3))
      Allocate (eltemp_adj(1:numcell,-eltcell(1):eltcell(1),-eltcell(2):eltcell(2),-eltcell(3):eltcell(3)), Stat = fail(4))
      Allocate (act_ele_cell(1:numcell,-1:1,-1:1,-1:1), old_ele_cell(1:numcell,-1:1,-1:1,-1:1)            , Stat = fail(5))
      Allocate (adjust(1:numcell,-1:1,-1:1,-1:1)                                                          , Stat = fail(6))

      If (Any(fail > 0)) Call error(1083)

      eltemp(:,:,:,:)       = 0.0_wp
      eltemp_adj(:,:,:,:)   = 0.0_wp
      gsource(:)            = 0.0_wp
      asource(:)            = 0.0_wp
      tempion(:)            = 0.0_wp
      ttmvom(:,:)           = 0.0_wp
      act_ele_cell(:,:,:,:) = 1.0_wp
      old_ele_cell(:,:,:,:) = 1.0_wp
      acell                 = ntsys(1)*ntsys(2)*ntsys(3)
      acell_old             = acell
      adjust                = .false.
      findepo               = .false.
      cel                   = 0
      kel                   = 0
      del                   = 0
      gel                   = 0

    End If

    keyres0 = 1

  End Subroutine allocate_ttm_arrays

  Subroutine deallocate_ttm_arrays()

    Integer, Dimension ( 1:5 ) :: fail

    fail = 0

    Deallocate (eltemp,eltemp_adj,asource,tempion,gsource,ttmvom,act_ele_cell,old_ele_cell,adjust, Stat = fail(1))
    If (kel>0) Deallocate(ketable,                                                                 Stat = fail(2))
    If (cel>0) Deallocate(cetable,                                                                 Stat = fail(3))
    If (del>0) Deallocate(detable,                                                                 Stat = fail(4))
    If (gel>0) Deallocate(gtable,                                                                  Stat = fail(5))

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

    Integer, Intent( In    ) :: i,j,k

    Integer                  :: idcube

    idcube = i + nprx * ( j + npry * k )

  End Function idcube
  
  Subroutine eltemp_sum (eltempsum,comm)

! Find sum of electronic temperatures over all active CET voxels

    Real ( Kind = wp ), Intent (   Out ) :: eltempsum
    Type ( comms_type), Intent ( InOut ) :: comm
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
                If (lrange) eltempsum = eltempsum + tmp
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    Call gsum (comm,eltempsum)

  End Subroutine eltemp_sum

  Subroutine eltemp_mean (eltempav,comm)

! Find mean electronic temperature over all active CET voxels

    Real ( Kind = wp ), Intent ( Out ) :: eltempav
    Real ( Kind = wp )                 :: tmp,acl
    Type( comms_type ), Intent ( InOut )    :: comm
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempav = 0.0_wp
    acl = 0.0_wp

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
                tmp = Merge (act_ele_cell (ijk,0,0,0), 1.0_wp, lcentre)
                If (lrange) Then
                  eltempav = eltempav + eltemp(ijk,ii,jj,kk) * tmp
                  acl = acl + tmp
                End If
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    Call gsum (comm, eltempav)
    Call gsum (comm, acl)

    If (acl>zero_plus) eltempav = eltempav/acl

  End Subroutine eltemp_mean

  Subroutine eltemp_maxKe (temp, eltempmax, comm)

! Find maximum temperature for calculating tabulated
! thermal conductivities (ionic or system) over all 
! active CET voxels (note that system temperature 
! applies over all CET voxels that do not overlap
! CIT voxels)

    Real ( Kind = wp ), Intent ( In )  :: temp
    Real ( Kind = wp ), Intent ( Out ) :: eltempmax
    Type( comms_type ), Intent ( InOut )  :: comm
    Real ( Kind = wp )                 :: eltempKe
    Integer                            :: i,j,k,ijk

    eltempmax = 0.0_wp

    Do k = 1, ntcell(3)
      Do j = 1, ntcell(2)
        Do i = 1, ntcell(1)
          ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          eltempKe = tempion(ijk)
          If (act_ele_cell(ijk,0,0,0)>zero_plus) eltempmax = Max (eltempmax, eltempKe)
        End Do
      End Do
    End Do
    eltempmax = Max (eltempmax, temp)

    Call gmax (comm,eltempmax)

  End Subroutine eltemp_maxKe

  Subroutine eltemp_max (eltempmax,comm)

! Find maximum electronic temperature over all
! active CET cells

    Real ( Kind = wp ), Intent ( Out ) :: eltempmax
    Type( comms_type ), Intent ( InOut )  :: comm
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

    Call gmax (comm,eltempmax)

  End Subroutine eltemp_max

  Subroutine eltemp_minKe (temp, eltempmin, comm)

! Find minimum temperature for calculating tabulated
! thermal conductivities (ionic or system) over all 
! active CET voxels (note that system temperature
! applies over all CET voxels that do not overlap
! CIT voxels)

    Type( comms_type ), Intent ( InOut )  :: comm
    Real ( Kind = wp ), Intent ( In )  :: temp
    Real ( Kind = wp ), Intent ( Out ) :: eltempmin
    Real ( Kind = wp )                 :: eltempKe
    Integer                            :: i,j,k,ijk

    eltempmin = 1.0e30_wp

    Do k = 1, ntcell(3)
      Do j = 1, ntcell(2)
        Do i = 1, ntcell(1)
          ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
          eltempKe = tempion (ijk)
          If (act_ele_cell(ijk,0,0,0)>zero_plus) eltempmin = Min (eltempmin, eltempKe)
        End Do
      End Do
    End Do

    eltempmin = Min (eltempmin, temp)

    Call gmin (comm,eltempmin)

  End Subroutine eltemp_minKe

  Subroutine eltemp_min (eltempmin,comm)

! Find minimum electronic temperature over all
! active CET cells

    Real( Kind = wp ), Intent ( Out ) :: eltempmin
    Type( comms_type ), Intent ( InOut ) :: comm
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
                If (ii==0 .and. jj==0 .and. kk==0) lrange = (lrange .and. (act_ele_cell(ijk,0,0,0)>zero_plus))
                If (lrange) eltempmin = Min (eltempmin, eltemp(ijk,ii,jj,kk))
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    Call gmin (comm,eltempmin)

  End Subroutine eltemp_min
  
    Subroutine depoinit(time,comm)

! determine initial energy deposition to electronic system,
! both temporally and spatially

    Real ( Kind = wp ), Intent( In ) :: time
    Type( comms_type), Intent( InOut ) :: comm
    Integer, Dimension( 1:3 ) :: fail
    Character ( Len = 14 ) :: number
    Character( Len = 256 ) :: message

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
      Call gaussianTrack(lat_U,comm)
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

    Write(number, '(f14.5)') depostart
    Write(message,"(6x,a,a,a)") &
      'electronic energy deposition starting at time = ',Trim(Adjustl(number)),' ps'
    Call info(message,.true.)
    Write(message,"(1x,130('-'))")
    Call info(message,.true.)

  End Subroutine depoinit
  
  Subroutine ttm_system_init(nstep,nsteql,keyres,dumpfile,time,temp,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing electronic temperature restart files
! at job termination or selected intervals in simulation
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,             Intent ( In ) :: keyres,nstep,nsteql
  Real ( Kind = wp ),  Intent ( In ) :: temp,time
  Character (Len = *), Intent ( In ) :: dumpfile
  Type( comms_type ), Intent( InOut ) :: comm

  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Logical                :: safe,l_tmp = .true.
  Integer                :: nxx,nyy,nzz,nstp,i,ix,iy,iz,ii,jj,kk,ijk,ipos(3)
  Real ( Kind = wp )     :: eltmp,tme,lat_sum,lat_max,lat_min
  Integer                :: iounit = 225
  Character( Len = 256 ) :: message

! check existence of readable restart file (DUMP_E)

  If (comm%idnode == 0) Inquire(File=dumpfile, Exist=l_tmp)
  Call gcheck(comm,l_tmp)
  If ((.not. l_tmp) .and. keyres==keyres0) Call error(684)

! if restarting simulation, read restart file

  If (l_tmp .and. keyres==keyres0) Then

    If (comm%idnode==0) Open (Unit=iounit, File=dumpfile)
    Call get_line(safe,iounit,record,comm); If (.not.safe) Goto 100
    Call get_word(record,word) ; nxx=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; nyy=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; nzz=Nint(word_2_real(word,0.0_wp))
    Call get_line(safe,iounit,record,comm); If (.not.safe) Goto 100
    Call get_word(record,word) ; nstp=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; tme=word_2_real(word,0.0_wp)
    Call get_word(record,word) ; depostart=word_2_real(word,0.0_wp)
    Call get_word(record,word) ; depoend=word_2_real(word,0.0_wp)
    ! check size of electronic temperature grid matches with size given in CONTROL file
    If (nxx/=eltsys(1) .or. nyy/=eltsys(2) .or. nzz/=eltsys(3)) Call error(685)
    ! check restart file is at same timestep as restart
    ! (can proceed if not, but need to warn user)
    If (nstp/=nstep .or. Abs(tme-time)>zero_plus) Call warning(520,0.0_wp,0.0_wp,0.0_wp)
    ! read in each line, find appropriate grid cell and assign
    ! electronic temperature if processor has that cell
    Do i=1,eltsys(1)*eltsys(2)*eltsys(3)
      Call get_line(safe,iounit,record,comm); If (.not.safe) Goto 100
      Call get_word(record,word) ; ipos(1)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; ipos(2)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; ipos(3)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; eltmp=word_2_real(word,0.0_wp)
      ix = ipos(1) + midI(1) - 1
      iy = ipos(2) + midI(2) - 1
      iz = ipos(3) + midI(3) - 1
      ii = Floor(Real(ix,Kind=wp)/Real(ntsys(1),Kind=wp))
      jj = Floor(Real(iy,Kind=wp)/Real(ntsys(2),Kind=wp))
      kk = Floor(Real(iz,Kind=wp)/Real(ntsys(3),Kind=wp))
      ix = Mod(ix+ntsys(1)*eltcell(1),ntsys(1)) + 1 - ntcelloff(1)
      iy = Mod(iy+ntsys(2)*eltcell(2),ntsys(2)) + 1 - ntcelloff(2)
      iz = Mod(iz+ntsys(3)*eltcell(3),ntsys(3)) + 1 - ntcelloff(3)
      If (ix>0 .and. ix<=ntcell(1) .and. iy>0 .and. iy<=ntcell(2) .and. iz>0 .and. iz<=ntcell(3)) Then
        ijk = 1 + ix + (ntcell(1)+2) * (iy + (ntcell(2)+2) * iz)
        eltemp(ijk,ii,jj,kk) = eltmp
      End If
    End Do
    ! fill boundary halo values and deal with required boundary conditions
    Call boundaryHalo (comm)
    Call boundaryCond (bcTypeE, temp, comm)
    ! check whether or not energy deposition has happened yet
    If(nstep>nsteql .and. time>=depostart .and. time<depoend) Then
      Call depoinit(time,comm)
    Else If (time>=depoend) Then
      findepo = .true.
    End If

    ! report successful reading and minimum, maximum and sums of
    ! electronic temperatures
    Call eltemp_sum (lat_sum,comm)
    Call eltemp_max (lat_max,comm)
    Call eltemp_min (lat_min,comm)
    Write(message,'(a)') 'electronic temperatures read from DUMP_E file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature (K) = ",ES11.4,&
               &/,1x,"maximum temperature (K) = ",ES11.4,&
               &/,1x,"sum of temperatures (K) = ",ES11.4)') &
               lat_min, lat_max, lat_sum
    Call info(message,.true.)
    If (comm%idnode==0) Then
      Close (iounit)
    End If

  Else

! if not restarting simulation, set electronic temperature grid
! to system temperature

    eltemp = temp

  End If

  Return

! Abnormal exit from electronic temperature dump file read

100 Continue

  Write(message,"(a)") dumpfile, ' data mishmash detected'
  Call error(686,message,.true.)
  Return

End Subroutine ttm_system_init

Subroutine ttm_system_revive    &
           (dumpfile,nstep,time,freq,nstrun,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing electronic temperature restart files
! at job termination or selected intervals in simulation
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton september 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Character (Len = *), Intent ( In ) :: dumpfile
  Integer, Intent ( In ) :: nstep,freq,nstrun
  Real(Kind=wp), Intent ( In ) :: time
  Type(comms_type), Intent(InOut) :: comm
  Integer :: iounit = 117
  Integer :: id,ii,jj,kk,imin,jmin,kmin,imax,jmax,kmax,i,j,k,ijk,ix,iy,iz
  Logical :: lrange

  If (freq /=0) Then
    If (Mod(nstep,freq)==0 .or. nstep==nstrun) Then

      If (comm%idnode==0) Then
        Open(Unit=iounit, File=dumpfile, Status='replace')
        Write(iounit,'(3i8)') eltsys(1),eltsys(2),eltsys(3)
        Write(iounit,'(i12,3(2x,es24.15))') nstep,time,depostart,depoend
        Close(iounit)
      End If
      Call gsync(comm)

      Do id=0,comm%mxnode-1
        If (comm%idnode==id) Then
          Open(Unit=iounit, File=dumpfile, Status='old', Position='append')
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
                  iz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
                  Do j = jmin, jmax
                    iy = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
                    Do i = imin, imax
                      ix = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                      lrange = (ix>0 .and. ix<=eltsys(1) .and. iy>0 .and. iy<=eltsys(2) .and. iz>0 .and. iz<=eltsys(3))
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      If (lrange) Write(iounit,'(3i8,2x,es24.15)') ix-midE(1),iy-midE(2),iz-midE(3),eltemp(ijk,ii,jj,kk)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          Close(iounit)
        End If
        Call gsync(comm)
      End Do

    End If 
  End If

End Subroutine ttm_system_revive


Subroutine ttm_table_read(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading specific heat capacity and coupling
! constant table files
!
! copyright - daresbury laboratory
! author    - m.a.seaton may 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Logical                :: safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: i
  Real( Kind = wp )      :: vk1,vk2
  Type( comms_type ), Intent( InOut ) :: comm
  Character( Len = 256 ) :: message

! read thermal conductivity data

  If (KeType == 3) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')

    i = 0
    Do While(i<kel)

      Call get_line(safe,ntable,record,comm)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          ketable(i,1) = vk1
          ketable(i,2) = vk2
        End If
      End If

    End Do

    If (comm%idnode==0) Then
      Close(Unit=ntable)
    End If
    Write(message,'(a)') 'thermal conductivity table read from Ke.dat file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature            (K) = ",ES12.4,&
               &/,1x,"maximum temperature            (K) = ",ES12.4,&
               &/,1x,"minimum t.c. value   (W m^-1 K^-1) = ",ES12.4,&
               &/,1x,"maximum t.c. value   (W m^-1 K^-1) = ",ES12.4)') &
               Minval(ketable(:,1)),Maxval(ketable(:,1)),Minval(ketable(:,2)),Maxval(ketable(:,2))
    Call info(message,.true.)

! convert thermal conductivity values from W m^-1 K^-1 to kB A^-1 ps^-1

    ketable(1:kel,2) = ketable(1:kel,2)*JKms_to_kBAps

  End If

! read volumetric heat capacity data

  If (CeType == 3 .or. CeType == 7) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')

    i = 0
    Do While(i<cel)

      Call get_line(safe,ntable,record,comm)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          cetable(i,1) = vk1
          cetable(i,2) = vk2
        End If
      End If

    End Do

    If (comm%idnode==0) Then
      Close(Unit=ntable)
    End If
    Write(message,'(a)') 'electronic volumetric heat capacity table read from Ce.dat file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature            (K) = ",ES12.4)')Minval(cetable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"maximum temperature            (K) = ",ES12.4)')Maxval(cetable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"minimum v.h.c. value (J m^-3 K^-1) = ",ES12.4)')Minval(cetable(:,2))
    Call info(message,.true.)
    Write(message,'(1x,"maximum v.h.c. value (J m^-3 K^-1) = ",ES12.4)') Maxval(cetable(:,2))
    Call info(message,.true.)

! convert volumetric heat capacity values from J m^-3 K^-1 to kB A^-3

    cetable(1:cel,2) = cetable(1:cel,2)*Jm3K_to_kBA3

  End If

! read thermal diffusivity data

  If (DeType == 3) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')

    i = 0
    Do While(i<del)

      Call get_line(safe,ntable,record,comm)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          detable(i,1) = vk1
          detable(i,2) = vk2
        End If
      End If

    End Do

    If (comm%idnode==0) Then
      Close(Unit=ntable)
    End If
    Write(message,'(a)') 'thermal diffusivity table read from De.dat file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature            (K) = ",ES12.4,&
               &/,1x,"maximum temperature            (K) = ",ES12.4,&
               &/,1x,"minimum diffusivity value  (m^2/s) = ",ES12.4,&
               &/,1x,"maximum diffusivity value  (m^2/s) = ",ES12.4)') &
               Minval(detable(:,1)),Maxval(detable(:,1)),Minval(detable(:,2)),Maxval(detable(:,2))
    Call info(message,.true.)

! convert thermal diffusivity values from m^2 s^-1 to A^2 ps^-1

    detable(1:del,2) = detable(1:del,2)*1e8_wp

  End If

! read coupling constant data

  If (gvar>0) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')

    i = 0
    Do While(i<gel)

      Call get_line(safe,ntable,record,comm)
      If (.not.safe) Then
        Go To 100
      Else
        Call get_word(record,word)
        vk1 = word_2_real(word)
        Call get_word(record,word)
        vk2 = word_2_real(word)
        If (vk1>=zero_plus) Then
          i=i+1
          gtable(i,1) = vk1
          gtable(i,2) = vk2
        End If
      End If

    End Do

    If (comm%idnode==0) Then
      Close(Unit=ntable)
    End If
    Write(message,'(a)') 'electron-phonon coupling table read from g.dat file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature            (K) = ",ES12.4)') Minval(gtable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"maximum temperature            (K) = ",ES12.4)') Maxval(gtable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"minimum e-p value    (W m^-3 K^-1) = ",ES12.4)')Minval(gtable(:,2))
    Call info(message,.true.)
    Write(message,'(1x,"maximum e-p value    (W m^-3 K^-1) = ",ES12.4)') Maxval(gtable(:,2))
    Call info(message,.true.)

! convert electron-phonon coupling values from W m^-3 K^-1 to ps^-1

    gtable(1:gel,2) = gtable(1:gel,2)*epc_to_chi

  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Then
    Close(Unit=ntable)
  End If
  Call error(682,master_only=.true.)

End Subroutine ttm_table_read

Subroutine ttm_table_scan(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for scanning specific heat capacity,
! thermal conductivity and electron-phonon coupling
! constant table files to determine numbers of data points
!
! copyright - daresbury laboratory
! author    - m.a.seaton may 2012
! contrib   - g.khara    may 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Type(comms_type), Intent(InOut) :: comm
  Logical                :: safe,lexist
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: fail
  Real( Kind = wp )      :: vk1,vk2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  Character( Len = 256 ) :: message

  If (l_ttm) Then

    fail=0
    Allocate (buffer(1:mxbuff), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'ttm_table_scan allocation failure'
       Call error(0,message)
    End If

! check existence of thermal conductivity table file

    If (KeType == 3) Then

      Inquire (File='Ke.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 100
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')
      End If

! determine number of lines of data to read

      kel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 5
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) kel=kel+1
        End If

      End Do
5  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (kel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(675)
      Else
        Allocate (ketable(1:kel,2), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        ketable(:,:) = 0.0_wp
      End If

    End If

! check existence of specific heat capacity table file

    If (CeType == 3) Then

      Inquire (File='Ce.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 200
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')
      End If

! determine number of lines of data to read

      cel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 10
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) cel=cel+1
        End If

      End Do
10  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (cel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(677)
      Else
        Allocate (cetable(1:cel,2), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        cetable(:,:) = 0.0_wp
      End If

    End If

! check existence of thermal diffusivity table file

    If (DeType == 3) Then

      Inquire (File='De.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 300
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')
      End If

! determine number of lines of data to read

      del = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 15
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) del=del+1
        End If

      End Do
15  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (del>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(679)
      Else
        Allocate (detable(1:del,2), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        detable(:,:) = 0.0_wp
      End If

    End If

! check existence of coupling constant table file

    If (gvar>0) Then

      Inquire (File='g.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 400
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')
      End If

! determine number of lines of data to read

      gel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 20
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) gel=gel+1
        End If

      End Do
20  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (gel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(681)
      Else
        Allocate (gtable(1:gel,2), Stat=fail) ! [GK] array length corrected
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        gtable(:,:) = 0.0_wp
      End If

    End If

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'ttm_table_scan deallocation failure'
       Call error(0,message)
    End If

  End If

  Return

! end of Ke.dat file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(674)

! end of Ce.dat file error exit

200 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(676)

! end of g.dat file error exit

300 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(678)

400 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(680)

End Subroutine ttm_table_scan

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

Subroutine uniformDist(lat_in)

! implement constant (homogeneous) spatial deposition

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
  Call gsum(comm,realdEdx)
  realdEdx = realdEdx/(Real(ntsys(3),Kind=wp)*delz)

  ! check if lattice sum equals the expected value

  If (comm%idnode == 0 .and. Abs((realdEdx-dEdX)/dEdX) > 0.01_wp) Then
    Call warning(540,Abs(realdEdx-dEdX)/dEdX*100_wp,0.0_wp,0.0_wp)
  End If
  
End Subroutine gaussianTrack
End Module ttm
