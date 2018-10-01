Module ttm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining arrays and initial parameters for 
! two-temperature model(ttm)
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton may 2012
! contrib   - g.khara may 2016
! contrib   - m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use setup
  Use configuration,   Only : configuration_type
  Use domains,         Only : domains_type,idcube
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
  Private

  Type, Public :: ttm_type

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

    Real ( Kind = wp ) :: Jm3K_to_kBA3,JKms_to_kBAps,kB_to_eV,eV_to_kB,mJcm2_to_eVA2,epc_to_chi
    Real ( Kind = wp ) :: cellrho,rcellrho,sysrho
    Real ( Kind = wp ) :: fluence,pdepth
    Real ( Kind = wp ) :: epthreshold = 1.1_wp

    Real( Kind = wp ), Allocatable, Dimension (:,:,:) :: lat_U,lat_B,lat_I
    Real( Kind = wp ) :: norm

    Logical :: trackInit = .false.
  End Type ttm_type

    Public :: allocate_ttm_arrays , deallocate_ttm_arrays
    Public :: eltemp_min, eltemp_max,&
      eltemp_maxKe,eltemp_minKe,eltemp_mean,eltemp_sum
    Public :: ttm_system_init,ttm_system_revive,ttm_table_read,&
      ttm_table_scan,boundaryHalo,boundaryCond,depoinit

  Contains

    Subroutine allocate_ttm_arrays(ttm,domain,config,comm)
    Type ( ttm_type ), Intent ( InOut ) :: ttm
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ), Intent( InOut ) :: config
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

    ttm%depostart = 0.0_wp
    ttm%depoend = 0.0_wp

! Setup constants based on fundamental values (found in
! setup.f90)

    ttm%JKms_to_kBAps = 10.0_wp/(boltz*tenunt)    ! convert W m^-1 K^-1 to kB A^-1 ps^-1
    ttm%Jm3K_to_kBA3  = 1.0e-7_wp/(boltz*tenunt)  ! convert J m^-3 K^-1 to kB A^-3
    ttm%kB_to_eV      = boltz/eu_ev               ! convert kB to eV
    ttm%eV_to_kB      = eu_ev/boltz               ! convert eV to kB
    ttm%mJcm2_to_eVA2 = 1.0e4_wp/(eu_ev*tenunt)   ! convert mJ cm^-2 to eV A^-2

    If (ttm%l_ttm) Then

! Determine number of ion temperature ttm%cells for domain and
! offsets for ion temperature determination

      start = config%cell(1)*Real(domain%idx,wp)*domain%nx_recip
      finish = config%cell(1)*Real(domain%idx+1,wp)*domain%nx_recip
      ttm%ntcell(1) = Ceiling(finish/ttm%delx) - Ceiling(start/ttm%delx)
      ttm%ntcelloff(1) = Ceiling(start/ttm%delx)
      ttm%zerocell(1) = 0.5_wp*config%cell(1) - ttm%delx*Real(Ceiling(start/ttm%delx),wp)

      start = config%cell(5)*Real(domain%idy,wp)*domain%ny_recip
      finish = config%cell(5)*Real(domain%idy+1,wp)*domain%ny_recip
      ttm%ntcell(2) = Ceiling(finish/ttm%dely) - Ceiling(start/ttm%dely)
      ttm%ntcelloff(2) = Ceiling(start/ttm%dely)
      ttm%zerocell(2) = 0.5_wp*config%cell(5) - ttm%dely*Real(Ceiling(start/ttm%dely),wp)

      start = config%cell(9)*Real(domain%idz,wp)*domain%nz_recip
      finish = config%cell(9)*Real(domain%idz+1,wp)*domain%nz_recip
      ttm%ntcell(3) = Ceiling(finish/ttm%delz) - Ceiling(start/ttm%delz)
      ttm%ntcelloff(3) = Ceiling(start/ttm%delz)
      ttm%zerocell(3) = 0.5_wp*config%cell(9) - ttm%delz*Real(Ceiling(start/ttm%delz),wp)

      ttm%numcell = (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)*(ttm%ntcell(3)+2)

! Determine mid-values for ion and electronic temperature grid

      ttm%midI(:) = INT((ttm%ntsys(:)+1)/2)
      ttm%midE(:) = INT((ttm%eltsys(:)+1)/2)

! Determine number of multiple ion temperature grids for electronic
! temperature grids

      ttm%eltcell(1) = Ceiling(Real(ttm%eltsys(1)-ttm%ntsys(1), Kind=wp)/Real(2*ttm%ntsys(1), Kind=wp))
      ttm%eltcell(2) = Ceiling(Real(ttm%eltsys(2)-ttm%ntsys(2), Kind=wp)/Real(2*ttm%ntsys(2), Kind=wp))
      ttm%eltcell(3) = Ceiling(Real(ttm%eltsys(3)-ttm%ntsys(3), Kind=wp)/Real(2*ttm%ntsys(3), Kind=wp))

! Determine positions of boundaries for electronic temperature grids

!   -x boundary
      numbc = -(ttm%eltsys(1)-ttm%ntsys(1))/2
      numbc = MOD(numbc+ttm%ntsys(1)*(ttm%eltcell(1)+1),ttm%ntsys(1)) + 1
      ttm%zeroE(1) = numbc - 1
      numbcmap = (ttm%eltsys(1)-ttm%ntsys(1))/2
      numbcmap = MOD(numbcmap+ttm%ntsys(1)*(ttm%eltcell(1)+1)-1,ttm%ntsys(1)) + 1
      If (numbc>ttm%ntcelloff(1) .and. numbc<=(ttm%ntcelloff(1)+ttm%ntcell(1))) Then
        ttm%ttmbc(1) = numbc - ttm%ntcelloff(1)
        Do i = 0, domain%nx-1
          start = config%cell(1)*Real(i,wp)*domain%nx_recip
          finish = config%cell(1)*Real(i+1,wp)*domain%nx_recip
          If (numbcmap>Ceiling(start/ttm%delx) .and. numbcmap<=Ceiling(finish/ttm%delx)) Then
            ttm%ttmbcmap(1) = idcube(i,domain%idy,domain%idz,domain)
          End If
        End Do
      Else
        ttm%ttmbc(1) = 0
        ttm%ttmbcmap(1) = -1
      End If

!   +x boundary
      numbc = (ttm%eltsys(1)-ttm%ntsys(1))/2
      numbc = MOD(numbc+ttm%ntsys(1)*(ttm%eltcell(1)+1)-1,ttm%ntsys(1)) + 1
      numbcmap = -(ttm%eltsys(1)-ttm%ntsys(1))/2
      numbcmap = MOD(numbcmap+ttm%ntsys(1)*(ttm%eltcell(1)+1),ttm%ntsys(1)) + 1
      If (numbc>ttm%ntcelloff(1) .and. numbc<=(ttm%ntcelloff(1)+ttm%ntcell(1))) Then
        ttm%ttmbc(2) = numbc - ttm%ntcelloff(1)
        Do i = 0, domain%nx-1
          start = config%cell(1)*Real(i,wp)*domain%nx_recip
          finish = config%cell(1)*Real(i+1,wp)*domain%nx_recip
          If (numbcmap>Ceiling(start/ttm%delx) .and. numbcmap<=Ceiling(finish/ttm%delx)) Then
            ttm%ttmbcmap(2) = idcube(i,domain%idy,domain%idz,domain)
          End If
        End Do
      Else
        ttm%ttmbc(2) = 0
        ttm%ttmbcmap(2) = -1
      End If
    
!   -y boundary
      numbc = -(ttm%eltsys(2)-ttm%ntsys(2))/2
      numbc = MOD(numbc+ttm%ntsys(2)*(ttm%eltcell(2)+1),ttm%ntsys(2)) + 1
      ttm%zeroE(2) = numbc - 1
      numbcmap = (ttm%eltsys(2)-ttm%ntsys(2))/2
      numbcmap = MOD(numbcmap+ttm%ntsys(2)*(ttm%eltcell(2)+1)-1,ttm%ntsys(2)) + 1
      If (numbc>ttm%ntcelloff(2) .and. numbc<=(ttm%ntcelloff(2)+ttm%ntcell(2))) Then
        ttm%ttmbc(3) = numbc - ttm%ntcelloff(2)
        Do i = 0, domain%ny-1
          start = config%cell(5)*Real(i,wp)*domain%ny_recip
          finish = config%cell(5)*Real(i+1,wp)*domain%ny_recip
          If (numbcmap>Ceiling(start/ttm%dely) .and. numbcmap<=Ceiling(finish/ttm%dely)) Then
            ttm%ttmbcmap(3) = idcube(domain%idx,i,domain%idz,domain)
          End If
        End Do
      Else
        ttm%ttmbc(3) = 0
        ttm%ttmbcmap(3) = -1
      End If

!   +y boundary
      numbc = (ttm%eltsys(2)-ttm%ntsys(2))/2
      numbc = MOD(numbc+ttm%ntsys(2)*(ttm%eltcell(2)+1)-1,ttm%ntsys(2)) + 1
      numbcmap = -(ttm%eltsys(2)-ttm%ntsys(2))/2
      numbcmap = MOD(numbcmap+ttm%ntsys(2)*(ttm%eltcell(2)+1),ttm%ntsys(2)) + 1
      If (numbc>ttm%ntcelloff(2) .and. numbc<=(ttm%ntcelloff(2)+ttm%ntcell(2))) Then
        ttm%ttmbc(4) = numbc - ttm%ntcelloff(2)
        Do i = 0, domain%ny-1
          start = config%cell(5)*Real(i,wp)*domain%ny_recip
          finish = config%cell(5)*Real(i+1,wp)*domain%ny_recip
          If (numbcmap>Ceiling(start/ttm%dely) .and. numbcmap<=Ceiling(finish/ttm%dely)) Then
            ttm%ttmbcmap(4) = idcube(domain%idx,i,domain%idz,domain)
          End If
        End Do
      Else
        ttm%ttmbc(4) = 0
        ttm%ttmbcmap(4) = -1
      End If

!   -z boundary
      numbc = -(ttm%eltsys(3)-ttm%ntsys(3))/2
      numbc = MOD(numbc+ttm%ntsys(3)*(ttm%eltcell(3)+1),ttm%ntsys(3)) + 1
      ttm%zeroE(3) = numbc - 1
      numbcmap = (ttm%eltsys(3)-ttm%ntsys(3))/2
      numbcmap = MOD(numbcmap+ttm%ntsys(3)*(ttm%eltcell(3)+1)-1,ttm%ntsys(3)) + 1
      If (numbc>ttm%ntcelloff(3) .and. numbc<=(ttm%ntcelloff(3)+ttm%ntcell(3))) Then
        ttm%ttmbc(5) = numbc - ttm%ntcelloff(3)
        Do i = 0, domain%nz-1
          start = config%cell(9)*Real(i,wp)*domain%nz_recip
          finish = config%cell(9)*Real(i+1,wp)*domain%nz_recip
          If (numbcmap>Ceiling(start/ttm%delz) .and. numbcmap<=Ceiling(finish/ttm%delz)) Then
            ttm%ttmbcmap(5) = idcube(domain%idx,domain%idy,i,domain)
          End If
        End Do
      Else
        ttm%ttmbc(5) = 0
        ttm%ttmbcmap(5) = -1
      End If

!   +z boundary
      numbc = (ttm%eltsys(3)-ttm%ntsys(3))/2
      numbc = MOD(numbc+ttm%ntsys(3)*(ttm%eltcell(3)+1)-1,ttm%ntsys(3)) + 1
      numbcmap = -(ttm%eltsys(3)-ttm%ntsys(3))/2
      numbcmap = MOD(numbcmap+ttm%ntsys(3)*(ttm%eltcell(3)+1),ttm%ntsys(3)) + 1
      If (numbc>ttm%ntcelloff(3) .and. numbc<=(ttm%ntcelloff(3)+ttm%ntcell(3))) Then
        ttm%ttmbc(6) = numbc - ttm%ntcelloff(3)
        Do i = 0, domain%nz-1
          start = config%cell(9)*Real(i,wp)*domain%nz_recip
          finish = config%cell(9)*Real(i+1,wp)*domain%nz_recip
          If (numbcmap>Ceiling(start/ttm%delz) .and. numbcmap<=Ceiling(finish/ttm%delz)) Then
            ttm%ttmbcmap(6) = idcube(domain%idx,domain%idy,i,domain)
          End If
        End Do
      Else
        ttm%ttmbc(6) = 0
        ttm%ttmbcmap(6) = -1
      End If

! Derived MPI datatypes for communication of temperatures (MPI 2.x+)

      If (comm%mxnode>1) Then
        Call MPI_TYPE_GET_EXTENT (wp_mpi, lb1, dbleth, ierr)
        Call MPI_TYPE_GET_EXTENT (MPI_INTEGER, lb2, inleth, ierr)
        Call MPI_TYPE_VECTOR (1, 1, 1, wp_mpi, basicslice, ierr)
        Call MPI_TYPE_VECTOR (1, 2, 2, MPI_INTEGER, bbasicslice, ierr)
        xlth = (ttm%ntcell(1) + 2) * dbleth
        ylth = (ttm%ntcell(2) + 2) * xlth
        Call MPI_TYPE_CREATE_HVECTOR (ttm%ntcell(2)    , 1, xlth, basicslice, oneslicex, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (ttm%ntcell(3)    , 1, ylth, oneslicex, ttm%tmpmsgx, ierr)
        Call MPI_TYPE_COMMIT (ttm%tmpmsgx, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (1, ttm%ntcell(1) + 2, xlth, basicslice, oneslicey, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (ttm%ntcell(3)    , 1, ylth, oneslicey, ttm%tmpmsgy, ierr)
        Call MPI_TYPE_COMMIT (ttm%tmpmsgy, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (1, ttm%ntcell(1) + 2, xlth, basicslice, oneslicez, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (1, ttm%ntcell(2) + 2, ylth, oneslicez, ttm%tmpmsgz, ierr)
        Call MPI_TYPE_COMMIT (ttm%tmpmsgz, ierr)
        bxlth = 2 * (ttm%ntcell(1) + 2) * inleth
        bylth = (ttm%ntcell(2) + 2) * bxlth
        Call MPI_TYPE_CREATE_HVECTOR (ttm%ntcell(2)    , 1, bxlth, bbasicslice, boneslicex, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (ttm%ntcell(3)    , 1, bylth, boneslicex, ttm%nummsgx, ierr)
        Call MPI_TYPE_COMMIT (ttm%nummsgx, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (1, ttm%ntcell(1) + 2, bxlth, bbasicslice, boneslicey, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (ttm%ntcell(3)    , 1, bylth, boneslicey, ttm%nummsgy, ierr)
        Call MPI_TYPE_COMMIT (ttm%nummsgy, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (1, ttm%ntcell(1) + 2, bxlth, bbasicslice, boneslicez, ierr)
        Call MPI_TYPE_CREATE_HVECTOR (1, ttm%ntcell(2) + 2, bylth, boneslicez, ttm%nummsgz, ierr)
        Call MPI_TYPE_COMMIT (ttm%nummsgz, ierr)
      Else
        ttm%tmpmsgx = 0; ttm%tmpmsgy = 0; ttm%tmpmsgz = 0
        ttm%nummsgx = 0; ttm%nummsgy = 0; ttm%nummsgz = 0
      End If

! Array allocation and initialization

      Allocate (ttm%eltemp(1:ttm%numcell,-ttm%eltcell(1):ttm%eltcell(1),-ttm%eltcell(2):ttm%eltcell(2),&
        -ttm%eltcell(3):ttm%eltcell(3))    , Stat = fail(1))
      Allocate (ttm%asource(1:ttm%numcell),ttm%tempion(1:ttm%numcell),ttm%gsource(1:ttm%numcell) &
        , Stat = fail(2))
      Allocate (ttm%ttmvom(1:ttm%numcell,1:4)                                                                     , Stat = fail(3))
      Allocate (ttm%eltemp_adj(1:ttm%numcell,-ttm%eltcell(1):ttm%eltcell(1),&
        -ttm%eltcell(2):ttm%eltcell(2),-ttm%eltcell(3):ttm%eltcell(3)), Stat = fail(4))
      Allocate (ttm%act_ele_cell(1:ttm%numcell,-1:1,-1:1,-1:1), &
      ttm%old_ele_cell(1:ttm%numcell,-1:1,-1:1,-1:1)            , Stat = fail(5))
      Allocate (ttm%adjust(1:ttm%numcell,-1:1,-1:1,-1:1)                                                          , Stat = fail(6))

      If (Any(fail > 0)) Call error(1083)

      ttm%eltemp(:,:,:,:)       = 0.0_wp
      ttm%eltemp_adj(:,:,:,:)   = 0.0_wp
      ttm%gsource(:)            = 0.0_wp
      ttm%asource(:)            = 0.0_wp
      ttm%tempion(:)            = 0.0_wp
      ttm%ttmvom(:,:)           = 0.0_wp
      ttm%act_ele_cell(:,:,:,:) = 1.0_wp
      ttm%old_ele_cell(:,:,:,:) = 1.0_wp
      ttm%acell                 = ttm%ntsys(1)*ttm%ntsys(2)*ttm%ntsys(3)
      ttm%acell_old             = ttm%acell
      ttm%adjust                = .false.
      ttm%findepo               = .false.
      ttm%cel                   = 0
      ttm%kel                   = 0
      ttm%del                  = 0
      ttm%gel                   = 0

    End If

    ttm%keyres0 = 1

  End Subroutine allocate_ttm_arrays

  Subroutine deallocate_ttm_arrays(ttm)

    Type ( ttm_type ), Intent ( InOut ) :: ttm
    Integer, Dimension ( 1:5 ) :: fail

    fail = 0

    Deallocate (ttm%eltemp,ttm%eltemp_adj,ttm%asource,ttm%tempion,ttm%gsource,&
      ttm%ttmvom,ttm%act_ele_cell,ttm%old_ele_cell,ttm%adjust, Stat = fail(1))
    If (ttm%kel>0) Deallocate(ttm%ketable,                                                                 Stat = fail(2))
    If (ttm%cel>0) Deallocate(ttm%cetable,                                                                 Stat = fail(3))
    If (ttm%del>0) Deallocate(ttm%detable,                                                                 Stat = fail(4))
    If (ttm%gel>0) Deallocate(ttm%gtable,                                                                  Stat = fail(5))

    If (Any(fail > 0)) Call error(1084)

  End Subroutine deallocate_ttm_arrays

  Subroutine eltemp_sum (eltempsum,ttm,comm)

! Find sum of electronic temperatures over all active CET voxels

    Type ( ttm_type ), Intent ( InOut ) :: ttm
    Real ( Kind = wp ), Intent (   Out ) :: eltempsum
    Type ( comms_type), Intent ( InOut ) :: comm
    Real ( Kind = wp )                 :: tmp
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempsum = 0.0_wp

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      If (ttm%eltcell(3)>0 .and. kk == -ttm%eltcell(3) .and. ttm%ttmbcmap(5)>=0) Then
        kmin = ttm%ttmbc(5)
      Else
        kmin = 1
      End If
      If (ttm%eltcell(3)>0 .and. kk == ttm%eltcell(3) .and. ttm%ttmbcmap(6)>=0) Then
        kmax = ttm%ttmbc(6)
      Else
        kmax = ttm%ntcell(3)
      End If

      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        If (ttm%eltcell(2)>0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3)>=0) Then
          jmin = ttm%ttmbc(3)
        Else
          jmin = 1
        End If
        If (ttm%eltcell(2)>0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4)>=0) Then
          jmax = ttm%ttmbc(4)
        Else
          jmax = ttm%ntcell(2)
        End If

        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          If (ttm%eltcell(1)>0 .and. ii == -ttm%eltcell(1) .and. ttm%ttmbcmap(1)>=0) Then
            imin = ttm%ttmbc(1)
          Else
            imin = 1
          End If
          If (ttm%eltcell(1)>0 .and. ii == ttm%eltcell(1) .and. ttm%ttmbcmap(2)>=0) Then
            imax = ttm%ttmbc(2)
          Else
            imax = ttm%ntcell(1)
          End If
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ttm%ntcelloff(3) + (kk + ttm%eltcell(3)) * ttm%ntsys(3) - ttm%zeroE(3)
            Do j = jmin, jmax
              ly = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
              Do i = imin, imax
                lx = i + ttm%ntcelloff(1) + (ii + ttm%eltcell(1)) * ttm%ntsys(1) - ttm%zeroE(1)
                lrange = (lx>0 .and. lx<=ttm%eltsys(1) .and. ly>0 .and. ly<=ttm%eltsys(2) .and. lz>0 .and. lz<=ttm%eltsys(3))
                ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                tmp = ttm%eltemp(ijk,ii,jj,kk) * Merge (ttm%act_ele_cell (ijk,0,0,0), 1.0_wp, lcentre)
                If (lrange) eltempsum = eltempsum + tmp
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    Call gsum (comm,eltempsum)

  End Subroutine eltemp_sum

  Subroutine eltemp_mean (eltempav,ttm,comm)

! Find mean electronic temperature over all active CET voxels

    Type ( ttm_type ), Intent ( InOut ) :: ttm
    Real ( Kind = wp ), Intent ( Out ) :: eltempav
    Real ( Kind = wp )                 :: tmp,acl
    Type( comms_type ), Intent ( InOut )    :: comm
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempav = 0.0_wp
    acl = 0.0_wp

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      If (ttm%eltcell(3)>0 .and. kk == -ttm%eltcell(3) .and. ttm%ttmbcmap(5)>=0) Then
        kmin = ttm%ttmbc(5)
      Else
        kmin = 1
      End If
      If (ttm%eltcell(3)>0 .and. kk == ttm%eltcell(3) .and. ttm%ttmbcmap(6)>=0) Then
        kmax = ttm%ttmbc(6)
      Else
        kmax = ttm%ntcell(3)
      End If

      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        If (ttm%eltcell(2)>0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3)>=0) Then
          jmin = ttm%ttmbc(3)
        Else
          jmin = 1
        End If
        If (ttm%eltcell(2)>0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4)>=0) Then
          jmax = ttm%ttmbc(4)
        Else
          jmax = ttm%ntcell(2)
        End If

        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          If (ttm%eltcell(1)>0 .and. ii == -ttm%eltcell(1) .and. ttm%ttmbcmap(1)>=0) Then
            imin = ttm%ttmbc(1)
          Else
            imin = 1
          End If
          If (ttm%eltcell(1)>0 .and. ii == ttm%eltcell(1) .and. ttm%ttmbcmap(2)>=0) Then
            imax = ttm%ttmbc(2)
          Else
            imax = ttm%ntcell(1)
          End If
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ttm%ntcelloff(3) + (kk + ttm%eltcell(3)) * ttm%ntsys(3) - ttm%zeroE(3)
            Do j = jmin, jmax
              ly = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
              Do i = imin, imax
                lx = i + ttm%ntcelloff(1) + (ii + ttm%eltcell(1)) * ttm%ntsys(1) - ttm%zeroE(1)
                lrange = (lx>0 .and. lx<=ttm%eltsys(1) .and. ly>0 .and. ly<=ttm%eltsys(2) .and. lz>0 .and. lz<=ttm%eltsys(3))
                ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                tmp = Merge (ttm%act_ele_cell (ijk,0,0,0), 1.0_wp, lcentre)
                If (lrange) Then
                  eltempav = eltempav + ttm%eltemp(ijk,ii,jj,kk) * tmp
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

  Subroutine eltemp_maxKe (temp, eltempmax, ttm,comm)

! Find maximum temperature for calculating tabulated
! thermal conductivities (ionic or system) over all 
! active CET voxels (note that system temperature 
! applies over all CET voxels that do not overlap
! CIT voxels)
   
    Type ( ttm_type ), Intent ( InOut ) :: ttm
    Real ( Kind = wp ), Intent ( In )  :: temp
    Real ( Kind = wp ), Intent ( Out ) :: eltempmax
    Type( comms_type ), Intent ( InOut )  :: comm
    Real ( Kind = wp )                 :: eltempKe
    Integer                            :: i,j,k,ijk

    eltempmax = 0.0_wp

    Do k = 1, ttm%ntcell(3)
      Do j = 1, ttm%ntcell(2)
        Do i = 1, ttm%ntcell(1)
          ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
          eltempKe = ttm%tempion(ijk)
          If (ttm%act_ele_cell(ijk,0,0,0)>zero_plus) eltempmax = Max (eltempmax, eltempKe)
        End Do
      End Do
    End Do
    eltempmax = Max (eltempmax, temp)

    Call gmax (comm,eltempmax)

  End Subroutine eltemp_maxKe

  Subroutine eltemp_max (eltempmax,ttm,comm)

! Find maximum electronic temperature over all
! active CET ttm%cells
    Type ( ttm_type ), Intent ( InOut ) :: ttm
    Real ( Kind = wp ), Intent ( Out ) :: eltempmax
    Type( comms_type ), Intent ( InOut )  :: comm
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange,lcentre

    eltempmax = 0.0_wp

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      If (ttm%eltcell(3)>0 .and. kk == -ttm%eltcell(3) .and. ttm%ttmbcmap(5)>=0) Then
        kmin = ttm%ttmbc(5)
      Else
        kmin = 1
      End If
      If (ttm%eltcell(3)>0 .and. kk == ttm%eltcell(3) .and. ttm%ttmbcmap(6)>=0) Then
        kmax = ttm%ttmbc(6)
      Else
        kmax = ttm%ntcell(3)
      End If

      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        If (ttm%eltcell(2)>0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3)>=0) Then
          jmin = ttm%ttmbc(3)
        Else
          jmin = 1
        End If
        If (ttm%eltcell(2)>0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4)>=0) Then
          jmax = ttm%ttmbc(4)
        Else
          jmax = ttm%ntcell(2)
        End If

        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          If (ttm%eltcell(1)>0 .and. ii == -ttm%eltcell(1) .and. ttm%ttmbcmap(1)>=0) Then
            imin = ttm%ttmbc(1)
          Else
            imin = 1
          End If
          If (ttm%eltcell(1)>0 .and. ii == ttm%eltcell(1) .and. ttm%ttmbcmap(2)>=0) Then
            imax = ttm%ttmbc(2)
          Else
            imax = ttm%ntcell(1)
          End If
          lcentre = (ii==0 .and. jj==0 .and. kk==0)
          Do k = kmin, kmax
            lz = k + ttm%ntcelloff(3) + (kk + ttm%eltcell(3)) * ttm%ntsys(3) - ttm%zeroE(3)
            Do j = jmin, jmax
              ly = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
              Do i = imin, imax
                lx = i + ttm%ntcelloff(1) + (ii + ttm%eltcell(1)) * ttm%ntsys(1) - ttm%zeroE(1)
                ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                lrange = (lx>0 .and. lx<=ttm%eltsys(1) .and. ly>0 .and. ly<=ttm%eltsys(2) .and. lz>0 .and. lz<=ttm%eltsys(3))
                If (lcentre) lrange = (lrange .and. (ttm%act_ele_cell(ijk,0,0,0)>zero_plus))
                If (lrange) eltempmax = Max (eltempmax, ttm%eltemp(ijk,ii,jj,kk))
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    Call gmax (comm,eltempmax)

  End Subroutine eltemp_max

  Subroutine eltemp_minKe (temp, eltempmin, ttm,comm)

! Find minimum temperature for calculating tabulated
! thermal conductivities (ionic or system) over all 
! active CET voxels (note that system temperature
! applies over all CET voxels that do not overlap
! CIT voxels)
    Type( ttm_type ), Intent( InOut ) :: ttm
    Type( comms_type ), Intent ( InOut )  :: comm
    Real ( Kind = wp ), Intent ( In )  :: temp
    Real ( Kind = wp ), Intent ( Out ) :: eltempmin
    Real ( Kind = wp )                 :: eltempKe
    Integer                            :: i,j,k,ijk

    eltempmin = 1.0e30_wp

    Do k = 1, ttm%ntcell(3)
      Do j = 1, ttm%ntcell(2)
        Do i = 1, ttm%ntcell(1)
          ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
          eltempKe = ttm%tempion (ijk)
          If (ttm%act_ele_cell(ijk,0,0,0)>zero_plus) eltempmin = Min (eltempmin, eltempKe)
        End Do
      End Do
    End Do

    eltempmin = Min (eltempmin, temp)

    Call gmin (comm,eltempmin)

  End Subroutine eltemp_minKe

  Subroutine eltemp_min (eltempmin,ttm,comm)

! Find minimum electronic temperature over all
! active CET ttm%cells

    Type( ttm_type ), Intent( InOut ) :: ttm
    Real( Kind = wp ), Intent ( Out ) :: eltempmin
    Type( comms_type ), Intent ( InOut ) :: comm
    Integer                            :: i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ijk,lx,ly,lz
    Logical                            :: lrange

    eltempmin = 1.0e30_wp

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      If (ttm%eltcell(3)>0 .and. kk == -ttm%eltcell(3) .and. ttm%ttmbcmap(5)>=0) Then
        kmin = ttm%ttmbc(5)
      Else
        kmin = 1
      End If
      If (ttm%eltcell(3)>0 .and. kk == ttm%eltcell(3) .and. ttm%ttmbcmap(6)>=0) Then
        kmax = ttm%ttmbc(6)
      Else
        kmax = ttm%ntcell(3)
      End If

      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        If (ttm%eltcell(2)>0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3)>=0) Then
          jmin = ttm%ttmbc(3)
        Else
          jmin = 1
        End If
        If (ttm%eltcell(2)>0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4)>=0) Then
          jmax = ttm%ttmbc(4)
        Else
          jmax = ttm%ntcell(2)
        End If

        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          If (ttm%eltcell(1)>0 .and. ii == -ttm%eltcell(1) .and. ttm%ttmbcmap(1)>=0) Then
            imin = ttm%ttmbc(1)
          Else
            imin = 1
          End If
          If (ttm%eltcell(1)>0 .and. ii == ttm%eltcell(1) .and. ttm%ttmbcmap(2)>=0) Then
            imax = ttm%ttmbc(2)
          Else
            imax = ttm%ntcell(1)
          End If
          Do k = kmin, kmax
            lz = k + ttm%ntcelloff(3) + (kk + ttm%eltcell(3)) * ttm%ntsys(3) - ttm%zeroE(3)
            Do j = jmin, jmax
              ly = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
              Do i = imin, imax
                lx = i + ttm%ntcelloff(1) + (ii + ttm%eltcell(1)) * ttm%ntsys(1) - ttm%zeroE(1)
                ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                lrange = (lx>0 .and. lx<=ttm%eltsys(1) .and. ly>0 .and. ly<=ttm%eltsys(2) .and. lz>0 .and. lz<=ttm%eltsys(3))
                If (ii==0 .and. jj==0 .and. kk==0) lrange = (lrange .and. (ttm%act_ele_cell(ijk,0,0,0)>zero_plus))
                If (lrange) eltempmin = Min (eltempmin, ttm%eltemp(ijk,ii,jj,kk))
              End Do
            End Do
          End Do

        End Do
      End Do
    End Do

    Call gmin (comm,eltempmin)

  End Subroutine eltemp_min
  
  Subroutine depoinit(time,ttm,comm)

! determine initial energy deposition to electronic system,
! both temporally and spatially
    Type( ttm_type ), Intent( inOut ) :: ttm
    Real ( Kind = wp ), Intent( In ) :: time
    Type( comms_type), Intent( InOut ) :: comm
    Integer, Dimension( 1:3 ) :: fail
    Character ( Len = 14 ) :: number
    Character( Len = 256 ) :: message

    fail = 0

    Allocate (ttm%lat_U (0:ttm%ntcell(1)+1,0:ttm%ntcell(2)+1,0:ttm%ntcell(3)+1), Stat = fail(1))
    Allocate (ttm%lat_B (0:ttm%ntcell(1)+1,0:ttm%ntcell(2)+1,0:ttm%ntcell(3)+1), Stat = fail(2))
    Allocate (ttm%lat_I (0:ttm%ntcell(1)+1,0:ttm%ntcell(2)+1,0:ttm%ntcell(3)+1), Stat = fail(3))

    If (Any(fail>0)) Call error(1089)

    ttm%lat_U(:,:,:) = 0.0_wp ! spatial deposition (eV)
    ttm%lat_B(:,:,:) = 0.0_wp ! temporal deposition of ttm%lat_U (eV)
    ttm%lat_I(:,:,:) = 0.0_wp ! sum of temporal deposition of ttm%lat_B (eV)

! spatial distribution of track

    Select Case (ttm%sdepoType)
    Case (1)
    ! Gaussian spatial deposition
      Call gaussianTrack(ttm%lat_U,ttm,comm)
     Case (2)
      ! Constant (flat) spatial deposition
      Call uniformDist(ttm%lat_U,ttm)
     Case (3)
      ! xy-flat, z-exp spatial deposition
      Call uniformDistZexp(ttm%lat_U,ttm)
    End Select

    ttm%trackInit = .true.                           ! switch on flag indicating track initialisation is in progress
    If (ttm%depostart<=zero_plus) ttm%depostart = time   ! time (ps) when deposition starts, i.e. current time

! temporal deposition of track: calculate time ttm%normalisation factor

    Select Case (ttm%tdepoType)
!   type=1: gauss(t)
    Case (1)
    ! Gaussian temporal deposition
      ttm%norm = 1.0_wp/(sqrpi*rt2*ttm%tdepo)
      ttm%depoend = ttm%depostart+2.0_wp*ttm%tcdepo*ttm%tdepo
    Case (2)
    ! decaying exponential temporal deposition
      ttm%norm = 1.0_wp/(1.0_wp-Exp(-ttm%tcdepo))
      ttm%depoend = ttm%depostart+2.0_wp*ttm%tcdepo*ttm%tdepo
    Case (3)
    ! delta temporal deposition
      ttm%norm = 1.0_wp
      ttm%depoend = ttm%depostart
    Case (4)
    ! pulse temporal deposition
      ttm%norm = 1.0_wp/ttm%tdepo
      ttm%depoend = ttm%depostart+ttm%tdepo
    End Select

    ! report start of energy deposition

    Write(number, '(f14.5)') ttm%depostart
    Write(message,"(6x,a,a,a)") &
      'electronic energy deposition starting at time = ',Trim(Adjustl(number)),' ps'
    Call info(message,.true.)
    Write(message,"(1x,130('-'))")
    Call info(message,.true.)

  End Subroutine depoinit

  Subroutine ttm_system_init(nstep,nsteql,keyres,dumpfile,time,temp,domain,ttm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing electronic temperature restart files
! at job termination or selected intervals in simulation
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( ttm_type ), Intent( InOut )   :: ttm 
    Integer,             Intent ( In ) :: keyres,nstep,nsteql
  Real ( Kind = wp ),  Intent ( In ) :: temp,time
  Character (Len = *), Intent ( In ) :: dumpfile
  Type( domains_type ), Intent( In    ) :: domain
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
  If ((.not. l_tmp) .and. keyres==ttm%keyres0) Call error(684)

! if restarting simulation, read restart file

  If (l_tmp .and. keyres==ttm%keyres0) Then

    If (comm%idnode==0) Open (Unit=iounit, File=dumpfile)
    Call get_line(safe,iounit,record,comm); If (.not.safe) Goto 100
    Call get_word(record,word) ; nxx=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; nyy=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; nzz=Nint(word_2_real(word,0.0_wp))
    Call get_line(safe,iounit,record,comm); If (.not.safe) Goto 100
    Call get_word(record,word) ; nstp=Nint(word_2_real(word,0.0_wp))
    Call get_word(record,word) ; tme=word_2_real(word,0.0_wp)
    Call get_word(record,word) ; ttm%depostart=word_2_real(word,0.0_wp)
    Call get_word(record,word) ; ttm%depoend=word_2_real(word,0.0_wp)
    ! check size of electronic temperature grid matches with size given in CONTROL file
    If (nxx/=ttm%eltsys(1) .or. nyy/=ttm%eltsys(2) .or. nzz/=ttm%eltsys(3)) Call error(685)
    ! check restart file is at same timestep as restart
    ! (can proceed if not, but need to warn user)
    If (nstp/=nstep .or. Abs(tme-time)>zero_plus) Call warning(520,0.0_wp,0.0_wp,0.0_wp)
    ! read in each line, find appropriate grid config%cell and assign
    ! electronic temperature if processor has that config%cell
    Do i=1,ttm%eltsys(1)*ttm%eltsys(2)*ttm%eltsys(3)
      Call get_line(safe,iounit,record,comm); If (.not.safe) Goto 100
      Call get_word(record,word) ; ipos(1)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; ipos(2)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; ipos(3)=Nint(word_2_real(word,0.0_wp))
      Call get_word(record,word) ; eltmp=word_2_real(word,0.0_wp)
      ix = ipos(1) + ttm%midI(1) - 1
      iy = ipos(2) + ttm%midI(2) - 1
      iz = ipos(3) + ttm%midI(3) - 1
      ii = Floor(Real(ix,Kind=wp)/Real(ttm%ntsys(1),Kind=wp))
      jj = Floor(Real(iy,Kind=wp)/Real(ttm%ntsys(2),Kind=wp))
      kk = Floor(Real(iz,Kind=wp)/Real(ttm%ntsys(3),Kind=wp))
      ix = Mod(ix+ttm%ntsys(1)*ttm%eltcell(1),ttm%ntsys(1)) + 1 - ttm%ntcelloff(1)
      iy = Mod(iy+ttm%ntsys(2)*ttm%eltcell(2),ttm%ntsys(2)) + 1 - ttm%ntcelloff(2)
      iz = Mod(iz+ttm%ntsys(3)*ttm%eltcell(3),ttm%ntsys(3)) + 1 - ttm%ntcelloff(3)
      If (ix>0 .and. ix<=ttm%ntcell(1) .and. iy>0 .and. iy<=ttm%ntcell(2) .and. iz>0 .and. iz<=ttm%ntcell(3)) Then
        ijk = 1 + ix + (ttm%ntcell(1)+2) * (iy + (ttm%ntcell(2)+2) * iz)
        ttm%eltemp(ijk,ii,jj,kk) = eltmp
      End If
    End Do
    ! fill boundary halo values and deal with required boundary conditions
    Call boundaryHalo(ttm,domain,comm)
    Call boundaryCond (ttm%bcTypeE, temp, ttm,comm)
    ! check whether or not energy deposition has happened yet
    If(nstep>nsteql .and. time>=ttm%depostart .and. time<ttm%depoend) Then
      Call depoinit(time,ttm,comm)
    Else If (time>=ttm%depoend) Then
      ttm%findepo = .true.
    End If

    ! report successful reading and minimum, maximum and sums of
    ! electronic temperatures
    Call eltemp_sum (lat_sum,ttm,comm)
    Call eltemp_max (lat_max,ttm,comm)
    Call eltemp_min (lat_min,ttm,comm)
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

    ttm%eltemp = temp

  End If

  Return

! Abttm%normal exit from electronic temperature dump file read

100 Continue

  Write(message,"(a)") dumpfile, ' data mishmash detected'
  Call error(686,message,.true.)
  Return

End Subroutine ttm_system_init

Subroutine ttm_system_revive    &
  (dumpfile,nstep,time,freq,nstrun,ttm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing electronic temperature restart files
! at job termination or selected intervals in simulation
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton september 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( ttm_type ), Intent( InOut )   :: ttm 
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
        Write(iounit,'(3i8)') ttm%eltsys(1),ttm%eltsys(2),ttm%eltsys(3)
        Write(iounit,'(i12,3(2x,es24.15))') nstep,time,ttm%depostart,ttm%depoend
        Close(iounit)
      End If
      Call gsync(comm)

      Do id=0,comm%mxnode-1
        If (comm%idnode==id) Then
          Open(Unit=iounit, File=dumpfile, Status='old', Position='append')
          Do kk = -ttm%eltcell(3), ttm%eltcell(3)
            If (ttm%eltcell(3)>0 .and. kk == -ttm%eltcell(3) .and. ttm%ttmbcmap(5)>=0) Then
              kmin = ttm%ttmbc(5)
            Else
              kmin = 1
            End If
            If (ttm%eltcell(3)>0 .and. kk == ttm%eltcell(3) .and. ttm%ttmbcmap(6)>=0) Then
              kmax = ttm%ttmbc(6)
            Else
              kmax = ttm%ntcell(3)
            End If
            Do jj = -ttm%eltcell(2), ttm%eltcell(2)
              If (ttm%eltcell(2)>0 .and. jj == -ttm%eltcell(2) .and. ttm%ttmbcmap(3)>=0) Then
                jmin = ttm%ttmbc(3)
              Else
                jmin = 1
              End If
              If (ttm%eltcell(2)>0 .and. jj == ttm%eltcell(2) .and. ttm%ttmbcmap(4)>=0) Then
                jmax = ttm%ttmbc(4)
              Else
                jmax = ttm%ntcell(2)
              End If
              Do ii = -ttm%eltcell(1), ttm%eltcell(1)
                If (ttm%eltcell(1)>0 .and. ii == -ttm%eltcell(1) .and. ttm%ttmbcmap(1)>=0) Then
                  imin = ttm%ttmbc(1)
                Else
                  imin = 1
                End If
                If (ttm%eltcell(1)>0 .and. ii == ttm%eltcell(1) .and. ttm%ttmbcmap(2)>=0) Then
                  imax = ttm%ttmbc(2)
                Else
                  imax = ttm%ntcell(1)
                End If
                Do k = kmin, kmax
                  iz = k + ttm%ntcelloff(3) + (kk + ttm%eltcell(3)) * ttm%ntsys(3) - ttm%zeroE(3)
                  Do j = jmin, jmax
                    iy = j + ttm%ntcelloff(2) + (jj + ttm%eltcell(2)) * ttm%ntsys(2) - ttm%zeroE(2)
                    Do i = imin, imax
                      ix = i + ttm%ntcelloff(1) + (ii + ttm%eltcell(1)) * ttm%ntsys(1) - ttm%zeroE(1)
                      lrange = (ix>0 .and. ix<=ttm%eltsys(1) .and. iy>0 .and. iy<=ttm%eltsys(2) .and. iz>0 .and. iz<=ttm%eltsys(3))
                      ijk = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
                      If (lrange) Write(iounit,'(3i8,2x,es24.15)') ix-ttm%midE(1),iy-ttm%midE(2),iz-ttm%midE(3),&
                        ttm%eltemp(ijk,ii,jj,kk)
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


Subroutine ttm_table_read(ttm,comm)

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
  Type( ttm_type ), Intent( InOut )   :: ttm 
  Logical                :: safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: i
  Real( Kind = wp )      :: vk1,vk2
  Type( comms_type ), Intent( InOut ) :: comm
  Character( Len = 256 ) :: message

! read thermal conductivity data

  If (ttm%KeType == 3) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')

    i = 0
    Do While(i<ttm%kel)

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
          ttm%ketable(i,1) = vk1
          ttm%ketable(i,2) = vk2
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
               Minval(ttm%ketable(:,1)),Maxval(ttm%ketable(:,1)),Minval(ttm%ketable(:,2)),Maxval(ttm%ketable(:,2))
    Call info(message,.true.)

! convert thermal conductivity values from W m^-1 K^-1 to kB A^-1 ps^-1

    ttm%ketable(1:ttm%kel,2) = ttm%ketable(1:ttm%kel,2)*ttm%JKms_to_kBAps

  End If

! read volumetric heat capacity data

  If (ttm%CeType == 3 .or. ttm%CeType == 7) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')

    i = 0
    Do While(i<ttm%cel)

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
          ttm%cetable(i,1) = vk1
          ttm%cetable(i,2) = vk2
        End If
      End If

    End Do

    If (comm%idnode==0) Then
      Close(Unit=ntable)
    End If
    Write(message,'(a)') 'electronic volumetric heat capacity table read from Ce.dat file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature            (K) = ",ES12.4)')Minval(ttm%cetable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"maximum temperature            (K) = ",ES12.4)')Maxval(ttm%cetable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"minimum v.h.c. value (J m^-3 K^-1) = ",ES12.4)')Minval(ttm%cetable(:,2))
    Call info(message,.true.)
    Write(message,'(1x,"maximum v.h.c. value (J m^-3 K^-1) = ",ES12.4)') Maxval(ttm%cetable(:,2))
    Call info(message,.true.)

! convert volumetric heat capacity values from J m^-3 K^-1 to kB A^-3

    ttm%cetable(1:ttm%cel,2) = ttm%cetable(1:ttm%cel,2)*ttm%Jm3K_to_kBA3

  End If

! read thermal diffusivity data

  If (ttm%DeType == 3) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')

    i = 0
    Do While(i<ttm%del)

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
          ttm%detable(i,1) = vk1
          ttm%detable(i,2) = vk2
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
               Minval(ttm%detable(:,1)),Maxval(ttm%detable(:,1)),Minval(ttm%detable(:,2)),Maxval(ttm%detable(:,2))
    Call info(message,.true.)

! convert thermal diffusivity values from m^2 s^-1 to A^2 ps^-1

    ttm%detable(1:ttm%del,2) = ttm%detable(1:ttm%del,2)*1e8_wp

  End If

! read coupling constant data

  If (ttm%gvar>0) Then

    If (comm%idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')

    i = 0
    Do While(i<ttm%gel)

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
          ttm%gtable(i,1) = vk1
          ttm%gtable(i,2) = vk2
        End If
      End If

    End Do

    If (comm%idnode==0) Then
      Close(Unit=ntable)
    End If
    Write(message,'(a)') 'electron-phonon coupling table read from g.dat file for two-temperature model'
    Call info(message,.true.)
    Write(message,'(1x,"minimum temperature            (K) = ",ES12.4)') Minval(ttm%gtable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"maximum temperature            (K) = ",ES12.4)') Maxval(ttm%gtable(:,1))
    Call info(message,.true.)
    Write(message,'(1x,"minimum e-p value    (W m^-3 K^-1) = ",ES12.4)')Minval(ttm%gtable(:,2))
    Call info(message,.true.)
    Write(message,'(1x,"maximum e-p value    (W m^-3 K^-1) = ",ES12.4)') Maxval(ttm%gtable(:,2))
    Call info(message,.true.)

! convert electron-phonon coupling values from W m^-3 K^-1 to ps^-1

    ttm%gtable(1:ttm%gel,2) = ttm%gtable(1:ttm%gel,2)*ttm%epc_to_chi

  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Then
    Close(Unit=ntable)
  End If
  Call error(682,master_only=.true.)

End Subroutine ttm_table_read

Subroutine ttm_table_scan(ttm,comm)

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


  Type( ttm_type ), Intent( InOut )   :: ttm 
  Type(comms_type), Intent(InOut) :: comm
  Logical                :: safe,lexist
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: fail
  Real( Kind = wp )      :: vk1,vk2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  Character( Len = 256 ) :: message

  If (ttm%l_ttm) Then

    fail=0
    Allocate (buffer(1:mxbuff), Stat=fail)
    If (fail > 0) Then
       Write(message,'(a)') 'ttm_table_scan allocation failure'
       Call error(0,message)
    End If

! check existence of thermal conductivity table file

    If (ttm%KeType == 3) Then

      Inquire (File='Ke.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 100
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='Ke.dat', Status='old')
      End If

! determine number of lines of data to read

      ttm%kel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 5
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) ttm%kel=ttm%kel+1
        End If

      End Do
5  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (ttm%kel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(675)
      Else
        Allocate (ttm%ketable(1:ttm%kel,2), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        ttm%ketable(:,:) = 0.0_wp
      End If

    End If

! check existence of specific heat capacity table file

    If (ttm%CeType == 3) Then

      Inquire (File='Ce.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 200
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='Ce.dat', Status='old')
      End If

! determine number of lines of data to read

      ttm%cel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 10
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) ttm%cel=ttm%cel+1
        End If

      End Do
10  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (ttm%cel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(677)
      Else
        Allocate (ttm%cetable(1:ttm%cel,2), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        ttm%cetable(:,:) = 0.0_wp
      End If

    End If

! check existence of thermal diffusivity table file

    If (ttm%DeType == 3) Then

      Inquire (File='De.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 300
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='De.dat', Status='old')
      End If

! determine number of lines of data to read

      ttm%del= 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 15
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) ttm%del=ttm%del+1
        End If

      End Do
15  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (ttm%del>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(679)
      Else
        Allocate (ttm%detable(1:ttm%del,2), Stat=fail)
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        ttm%detable(:,:) = 0.0_wp
      End If

    End If

! check existence of coupling constant table file

    If (ttm%gvar>0) Then

      Inquire (File='g.dat', Exist=lexist)
      Call gcheck(comm,lexist)

      If (.not.lexist) Then
        Go To 400
      Else
        If (comm%idnode == 0) Open(Unit=ntable, File='g.dat', Status='old')
      End If

! determine number of lines of data to read

      ttm%gel = 0
      Do While(.true.)

        Call get_line(safe,ntable,record,comm)
        If (.not.safe) Then
          Go To 20
        Else
          Call get_word(record,word)
          vk1 = word_2_real(word)
          Call get_word(record,word)
          vk2 = word_2_real(word)
          If (vk1>=zero_plus) ttm%gel=ttm%gel+1
        End If

      End Do
20  Continue

      If (comm%idnode == 0) Close(Unit=ntable)

! check number of data lines and allocate array

      safe = (ttm%gel>0)
      Call gcheck(comm,safe)
      If (.not. safe) Then
        Call error(681)
      Else
        Allocate (ttm%gtable(1:ttm%gel,2), Stat=fail) ! [GK] array length corrected
        If (fail > 0) Then
          Write(message,'(a)') 'ttm_table_scan allocation failure'
          Call error(0,message)
        End If
        ttm%gtable(:,:) = 0.0_wp
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

! fills halo regions of electronic temperature lattice from neighbouring sections
! (periodic boundary conditions)
Subroutine boundaryHalo(ttm,domain,comm)
  Type( ttm_type ), Intent( InOut )   :: ttm 
  Type( domains_type ), Intent( In    ) :: domain
  Type(comms_type), Intent(In) :: comm

  Integer :: i,ii,iii1,iii2,j,jj,jjj1,jjj2,k,kk,kkk1,kkk2,ijk1,ijk2
  Integer, Dimension(4) :: req
  Integer, Allocatable :: stats(:,:)

  Integer :: ierr

  Allocate (stats(1:MPI_STATUS_SIZE,1:4))
  If (comm%mxnode>1) Then

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          If (domain%idx==domain%nx-1) Then
            iii1 = Mod(ii+3*ttm%eltcell(1), (ttm%eltcell(1)*2+1)) - ttm%eltcell(1)
          Else
            iii1 = ii
          End If
          If (domain%idx==0) Then
            iii2 = Mod(ii+ttm%eltcell(1)+1, (ttm%eltcell(1)*2+1)) - ttm%eltcell(1)
          Else
            iii2 = ii
          End If
          ijk1 = 2 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
          ijk2 = 1 + (ttm%ntcell(1)+1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
          Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk)  , 1, ttm%tmpmsgx, domain%map(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
          Call MPI_IRECV (ttm%eltemp(ijk2,iii1,jj,kk), 1, ttm%tmpmsgx, domain%map(2), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
          ijk1 = 1 + (ttm%ntcell(1)) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
          ijk2 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
          Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk)  , 1, ttm%tmpmsgx, domain%map(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
          Call MPI_IRECV (ttm%eltemp(ijk2,iii2,jj,kk), 1, ttm%tmpmsgx, domain%map(1), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
          Call MPI_WAITALL (4, req, stats, ierr)
        End Do
      End Do
    End Do

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      Do ii = -ttm%eltcell(1), ttm%eltcell(1)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          If (domain%idy==domain%ny-1) Then
            jjj1 = Mod(jj+3*ttm%eltcell(2), (ttm%eltcell(2)*2+1)) - ttm%eltcell(2)
          Else
            jjj1 = jj
          End If
          If (domain%idy==0) Then
            jjj2 = Mod(jj+ttm%eltcell(2)+1, (ttm%eltcell(2)*2+1)) - ttm%eltcell(2)
          Else
            jjj2 = jj
          End If
          ijk1 = 1 + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
          ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2)+2))
          Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk)  , 1, ttm%tmpmsgy, domain%map(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
          Call MPI_IRECV (ttm%eltemp(ijk2,ii,jjj1,kk), 1, ttm%tmpmsgy, domain%map(4), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
          ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + (ttm%ntcell(2)+2))
          ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
          Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk)  , 1, ttm%tmpmsgy, domain%map(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
          Call MPI_IRECV (ttm%eltemp(ijk2,ii,jjj2,kk), 1, ttm%tmpmsgy, domain%map(3), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
          Call MPI_WAITALL (4, req, stats, ierr)
        End Do
      End Do
    End Do

    Do jj = -ttm%eltcell(2), ttm%eltcell(2)
      Do ii = -ttm%eltcell(1), ttm%eltcell(1)
        Do kk = -ttm%eltcell(3), ttm%eltcell(3)
          If (domain%idz==domain%nz-1) Then
            kkk1 = Mod(kk+3*ttm%eltcell(3), (ttm%eltcell(3)*2+1)) - ttm%eltcell(3)
          Else
            kkk1 = kk
          End If
          If (domain%idz==0) Then
            kkk2 = Mod(kk+ttm%eltcell(3)+1, (ttm%eltcell(3)*2+1)) - ttm%eltcell(3)
          Else
            kkk2 = kk
          End If
          ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2)
          ijk2 = 1 + (ttm%ntcell(1)+2) * ((ttm%ntcell(2)+2) * (ttm%ntcell(3) + 1))
          Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk)  , 1, ttm%tmpmsgz, domain%map(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
          Call MPI_IRECV (ttm%eltemp(ijk2,ii,jj,kkk1), 1, ttm%tmpmsgz, domain%map(6), Grid1_tag, MPI_COMM_WORLD, req(2), ierr)
          ijk1 = 1 + (ttm%ntcell(1)+2) * ((ttm%ntcell(2)+2) * ttm%ntcell(3))
          ijk2 = 1
          Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk)  , 1, ttm%tmpmsgz, domain%map(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
          Call MPI_IRECV (ttm%eltemp(ijk2,ii,jj,kkk2), 1, ttm%tmpmsgz, domain%map(5), Grid2_tag, MPI_COMM_WORLD, req(4), ierr)
          Call MPI_WAITALL (4, req, stats, ierr)
        End Do
      End Do
    End Do

  Else

    Do kk = -ttm%eltcell(3), ttm%eltcell(3)
      kkk1 = Mod(kk+3*ttm%eltcell(3), (ttm%eltcell(3)*2+1)) - ttm%eltcell(3)
      kkk2 = Mod(kk+ttm%eltcell(3)+1, (ttm%eltcell(3)*2+1)) - ttm%eltcell(3)
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        jjj1 = Mod(jj+3*ttm%eltcell(2), (ttm%eltcell(2)*2+1)) - ttm%eltcell(2)
        jjj2 = Mod(jj+ttm%eltcell(2)+1, (ttm%eltcell(2)*2+1)) - ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          iii1 = Mod(ii+3*ttm%eltcell(1), (ttm%eltcell(1)*2+1)) - ttm%eltcell(1)
          iii2 = Mod(ii+ttm%eltcell(1)+1, (ttm%eltcell(1)*2+1)) - ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 2 + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
              ijk2 = 2 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
              ttm%eltemp(ijk2,iii1,jj,kk) = ttm%eltemp(ijk1,ii,jj,kk)
              ijk1 = 1 + ttm%ntcell(1) + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
              ijk2 = 1 + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * k)
              ttm%eltemp(ijk2,iii2,jj,kk) = ttm%eltemp(ijk1,ii,jj,kk)
            End Do
          End Do            
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2) * k)
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + 1 + (ttm%ntcell(2)+2) * k)
              ttm%eltemp(ijk2,ii,jjj1,kk) = ttm%eltemp(ijk1,ii,jj,kk)
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ntcell(2) + (ttm%ntcell(2)+2) * k)
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ntcell(2)+2) * k
              ttm%eltemp(ijk2,ii,jjj2,kk) = ttm%eltemp(ijk1,ii,jj,kk)
            End Do
          End Do
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2))
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * (ttm%ntcell(3)+1))
              ttm%eltemp(ijk2,ii,jj,kkk1) = ttm%eltemp(ijk1,ii,jj,kk)
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ntcell(2)+2) * ttm%ntcell(3))
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * j
              ttm%eltemp(ijk2,ii,jj,kkk2) = ttm%eltemp(ijk1,ii,jj,kk)
            End Do
          End Do
        End Do
      End Do
    End Do

  End If

  Deallocate (stats)

End Subroutine boundaryHalo

Subroutine boundaryCond (key, temp,ttm,comm)

! appends halo regions of entire electronic temperature lattice with appropriate boundary conditions


  Type( ttm_type ), Intent( InOut )   :: ttm 
  Real ( Kind = wp ), Intent ( In ) :: temp
  Integer,            Intent ( In ) :: key
  Type(comms_type),   Intent ( In ) :: comm

  Integer :: i,ii,j,jj,k,kk,ijk1,ijk2,ierr
  Integer, Dimension(4) :: req
  Integer, Dimension(MPI_STATUS_SIZE,4) :: stat


  Select Case (key)
! Periodic boundary conditions
  Case (1)
    If (ttm%ttmbcmap(1)>=0 .or. ttm%ttmbcmap(2)>=0) Then
      If (comm%mxnode>1) Then
        Do kk = -ttm%eltcell(3), ttm%eltcell(3)
          Do jj = -ttm%eltcell(2), ttm%eltcell(2)
            If (ttm%ttmbcmap(1)>=0) Then
              ijk1 = 1 + ttm%ttmbc(1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
              ijk2 = 1 + (ttm%ttmbc(1) - 1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
              ii = -ttm%eltcell(1)
              Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk) , 1, ttm%tmpmsgx, ttm%ttmbcmap(1), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
              Call MPI_IRECV (ttm%eltemp(ijk2,-ii,jj,kk), 1, ttm%tmpmsgx, ttm%ttmbcmap(1), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
            End If
            If (ttm%ttmbcmap(2)>=0) Then
              ijk1 = 1 + ttm%ttmbc(2) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
              ijk2 = 1 + (ttm%ttmbc(2) + 1) + (ttm%ntcell(1)+2) * (1 + (ttm%ntcell(2)+2))
              ii = ttm%eltcell(1)
              Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk) , 1, ttm%tmpmsgx, ttm%ttmbcmap(2), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
              Call MPI_IRECV (ttm%eltemp(ijk2,-ii,jj,kk), 1, ttm%tmpmsgx, ttm%ttmbcmap(2), Grid1_tag, MPI_COMM_WORLD, req(4), ierr)
            End If
            Call MPI_WAITALL (4, req, stat, ierr)
          End Do
        End Do
      Else
        Do kk = -ttm%eltcell(3), ttm%eltcell(3)
          Do jj = -ttm%eltcell(2), ttm%eltcell(2)
            Do k = 1, ttm%ntcell(3)
              Do j = 1, ttm%ntcell(2)
                ijk1 = 1 + ttm%ttmbc(1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
                ijk2 = 1 + (ttm%ttmbc(2) + 1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
                ttm%eltemp(ijk2,ttm%eltcell(1),jj,kk) = ttm%eltemp(ijk1,-ttm%eltcell(1),jj,kk)
                ijk1 = 1 + ttm%ttmbc(2) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
                ijk2 = 1 + (ttm%ttmbc(1) - 1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
                ttm%eltemp(ijk2,-ttm%eltcell(1),jj,kk) = ttm%eltemp(ijk1,ttm%eltcell(1),jj,kk)
              End Do
            End Do
          End Do
        End Do
      End If
    End If

    If (ttm%ttmbcmap(3)>=0 .or. ttm%ttmbcmap(4)>=0) Then
      If (comm%mxnode>1) Then
        Do kk = -ttm%eltcell(3), ttm%eltcell(3)
          Do ii = -ttm%eltcell(1), ttm%eltcell(1)
            If (ttm%ttmbcmap(3)>=0) Then
              ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) + (ttm%ntcell(2)+2))
              ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) - 1 + (ttm%ntcell(2)+2))
              jj = -ttm%eltcell(2)
              Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk) , 1, ttm%tmpmsgy, ttm%ttmbcmap(3), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
              Call MPI_IRECV (ttm%eltemp(ijk2,ii,-jj,kk), 1, ttm%tmpmsgy, ttm%ttmbcmap(3), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
            End If
            If (ttm%ttmbcmap(4)>=0) Then
              ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + (ttm%ntcell(2)+2))
              ijk2 = 1 + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + 1 + (ttm%ntcell(2)+2))
              jj = ttm%eltcell(2)
              Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk) , 1, ttm%tmpmsgy, ttm%ttmbcmap(4), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
              Call MPI_IRECV (ttm%eltemp(ijk2,ii,-jj,kk), 1, ttm%tmpmsgy, ttm%ttmbcmap(4), Grid1_tag, MPI_COMM_WORLD, req(4), ierr)
            End If
            Call MPI_WAITALL (4, req, stat, ierr)
          End Do
        End Do
      Else
        Do kk = -ttm%eltcell(3), ttm%eltcell(3)
          Do ii = -ttm%eltcell(1), ttm%eltcell(1)
            Do k = 1, ttm%ntcell(3)
              Do i = 0, ttm%ntcell(1)+1
                ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) + k * (ttm%ntcell(2)+2))
                ijk2 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + 1 + k * (ttm%ntcell(2)+2))
                ttm%eltemp(ijk2,ii,ttm%eltcell(2),kk) = ttm%eltemp(ijk1,ii,-ttm%eltcell(2),kk)
                ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + k * (ttm%ntcell(2)+2))
                ijk2 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) - 1 + k * (ttm%ntcell(2)+2))
                ttm%eltemp(ijk2,ii,-ttm%eltcell(2),kk) = ttm%eltemp(ijk1,ii,ttm%eltcell(2),kk)
              End Do
            End Do
          End Do
        End Do
      End If
    End If

    If (ttm%ttmbcmap(5)>=0 .or. ttm%ttmbcmap(6)>=0) Then
      If (comm%mxnode>1) Then
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do ii = -ttm%eltcell(1), ttm%eltcell(1)
            If (ttm%ttmbcmap(5)>=0) Then
              ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ttmbc(5) * (ttm%ntcell(2)+2))
              ijk2 = 1 + (ttm%ntcell(1)+2) * ((ttm%ttmbc(5) - 1) * (ttm%ntcell(2)+2))
              kk = -ttm%eltcell(3)
              Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk) , 1, ttm%tmpmsgz, ttm%ttmbcmap(5), Grid1_tag, MPI_COMM_WORLD, req(1), ierr)
              Call MPI_IRECV (ttm%eltemp(ijk2,ii,jj,-kk), 1, ttm%tmpmsgz, ttm%ttmbcmap(5), Grid2_tag, MPI_COMM_WORLD, req(2), ierr)
            End If
            If (ttm%ttmbcmap(6)>=0) Then
              ijk1 = 1 + (ttm%ntcell(1)+2) * (ttm%ttmbc(6) * (ttm%ntcell(2)+2))
              ijk2 = 1 + (ttm%ntcell(1)+2) * ((ttm%ttmbc(6) + 1) * (ttm%ntcell(2)+2))
              kk = ttm%eltcell(3)
              Call MPI_ISEND (ttm%eltemp(ijk1,ii,jj,kk) , 1, ttm%tmpmsgz, ttm%ttmbcmap(6), Grid2_tag, MPI_COMM_WORLD, req(3), ierr)
              Call MPI_IRECV (ttm%eltemp(ijk2,ii,jj,-kk), 1, ttm%tmpmsgz, ttm%ttmbcmap(6), Grid1_tag, MPI_COMM_WORLD, req(4), ierr)
            End If
            Call MPI_WAITALL (4, req, stat, ierr)
          End Do
        End Do
      Else
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do ii = -ttm%eltcell(1), ttm%eltcell(1)
            Do j = 0, ttm%ntcell(2)+1
              Do i = 0, ttm%ntcell(1)+1
                ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(5) * (ttm%ntcell(2)+2))
                ijk2 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ttmbc(6) + 1) * (ttm%ntcell(2)+2))
                ttm%eltemp(ijk2,ii,jj,ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,-ttm%eltcell(3))
                ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(6) * (ttm%ntcell(2)+2))
                ijk2 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ttmbc(5) - 1) * (ttm%ntcell(2)+2))
                ttm%eltemp(ijk2,ii,jj,-ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,ttm%eltcell(3))
              End Do
            End Do
          End Do
        End Do
      End If

    End If

! Infinite sink/source (Dirichlet) boundary conditions
  Case (2)
    If (ttm%ttmbcmap(1)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk2 = 1 + (ttm%ttmbc(1)-1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,-ttm%eltcell(1),jj,kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(2)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk2 = 1 + (ttm%ttmbc(2)+1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ttm%eltcell(1),jj,kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(3)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * ((ttm%ttmbc(3)-1) + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ii,-ttm%eltcell(2),kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(4)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * ((ttm%ttmbc(4)+1) + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ii,ttm%eltcell(2),kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(5)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ttmbc(5)-1) * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ii,jj,-ttm%eltcell(3)) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(6)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * (j + (ttm%ttmbc(6)+1) * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ii,jj,ttm%eltcell(3)) = temp
            End Do
          End Do
        End Do
      End Do
    End If

! 'Confined' (von Neumann) boundary conditions
  Case (3)
    If (ttm%ttmbcmap(1)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 1 + ttm%ttmbc(1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - 1
              ttm%eltemp(ijk2,-ttm%eltcell(1),jj,kk) = ttm%eltemp(ijk1,-ttm%eltcell(1),jj,kk)
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(2)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 1 + ttm%ttmbc(2) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + 1
              ttm%eltemp(ijk2,ttm%eltcell(1),jj,kk) = ttm%eltemp(ijk1,ttm%eltcell(1),jj,kk)
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(3)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)
              ttm%eltemp(ijk2,ii,-ttm%eltcell(2),kk) = ttm%eltemp(ijk1,ii,-ttm%eltcell(2),kk)
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(4)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)
              ttm%eltemp(ijk2,ii,ttm%eltcell(2),kk) = ttm%eltemp(ijk1,ii,ttm%eltcell(2),kk)
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(5)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(5) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,-ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,-ttm%eltcell(3))
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(6)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(6) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,ttm%eltcell(3))
            End Do
          End Do
        End Do
      End Do
    End If

! Mixed case: Infinite sink/source (Dirichlet) boundaries in x/y-directions
!             'Confined' (von Neumann) boundary in z-direction
  Case (4)
    If (ttm%ttmbcmap(1)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk2 = 1 + (ttm%ttmbc(1)-1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,-ttm%eltcell(1),jj,kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(2)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk2 = 1 + (ttm%ttmbc(2)+1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ttm%eltcell(1),jj,kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(3)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * ((ttm%ttmbc(3)-1) + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ii,-ttm%eltcell(2),kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(4)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk2 = 1 + i + (ttm%ntcell(1)+2) * ((ttm%ttmbc(4)+1) + k * (ttm%ntcell(2)+2))
              ttm%eltemp(ijk2,ii,ttm%eltcell(2),kk) = temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(5)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(5) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,-ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,-ttm%eltcell(3))
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(6)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(6) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,ttm%eltcell(3))
            End Do
          End Do
        End Do
      End Do
    End If

! Robin boundary conditions
  Case (5)
    If (ttm%ttmbcmap(1)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 1 + ttm%ttmbc(1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - 1
              ttm%eltemp(ijk2,-ttm%eltcell(1),jj,kk) = ttm%fluxout*(ttm%eltemp(ijk1,-ttm%eltcell(1),jj,kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(2)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 1 + ttm%ttmbc(2) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + 1
              ttm%eltemp(ijk2,ttm%eltcell(1),jj,kk) = ttm%fluxout*(ttm%eltemp(ijk1,ttm%eltcell(1),jj,kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(3)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)
              ttm%eltemp(ijk2,ii,-ttm%eltcell(2),kk) = ttm%fluxout*(ttm%eltemp(ijk1,ii,-ttm%eltcell(2),kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(4)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)
              ttm%eltemp(ijk2,ii,ttm%eltcell(2),kk) = ttm%fluxout*(ttm%eltemp(ijk1,ii,ttm%eltcell(2),kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(5)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(5) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,-ttm%eltcell(3)) = ttm%fluxout*(ttm%eltemp(ijk1,ii,jj,-ttm%eltcell(3))-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(6)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(6) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,ttm%eltcell(3)) = ttm%fluxout*(ttm%eltemp(ijk1,ii,jj,ttm%eltcell(3))-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

! Mixed case: Robin boundaries in x/y-directions
!             'Confined' (von Neumann) boundary in z-direction
  Case (6)
    If (ttm%ttmbcmap(1)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 1 + ttm%ttmbc(1) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - 1
              ttm%eltemp(ijk2,-ttm%eltcell(1),jj,kk) = ttm%fluxout*(ttm%eltemp(ijk1,-ttm%eltcell(1),jj,kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(2)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do jj = -ttm%eltcell(2), ttm%eltcell(2)
          Do k = 1, ttm%ntcell(3)
            Do j = 1, ttm%ntcell(2)
              ijk1 = 1 + ttm%ttmbc(2) + (ttm%ntcell(1)+2) * (j + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + 1
              ttm%eltemp(ijk2,ttm%eltcell(1),jj,kk) = ttm%fluxout*(ttm%eltemp(ijk1,ttm%eltcell(1),jj,kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(3)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(3) + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)
              ttm%eltemp(ijk2,ii,-ttm%eltcell(2),kk) = ttm%fluxout*(ttm%eltemp(ijk1,ii,-ttm%eltcell(2),kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(4)>=0) Then
      Do kk = -ttm%eltcell(3), ttm%eltcell(3)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do k = 1, ttm%ntcell(3)
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (ttm%ttmbc(4) + k * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)
              ttm%eltemp(ijk2,ii,ttm%eltcell(2),kk) = ttm%fluxout*(ttm%eltemp(ijk1,ii,ttm%eltcell(2),kk)-temp) + temp
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(5)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(5) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 - (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,-ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,-ttm%eltcell(3))
            End Do
          End Do
        End Do
      End Do
    End If

    If (ttm%ttmbcmap(6)>=0) Then
      Do jj = -ttm%eltcell(2), ttm%eltcell(2)
        Do ii = -ttm%eltcell(1), ttm%eltcell(1)
          Do j = 0, ttm%ntcell(2)+1
            Do i = 0, ttm%ntcell(1)+1
              ijk1 = 1 + i + (ttm%ntcell(1)+2) * (j + ttm%ttmbc(6) * (ttm%ntcell(2)+2))
              ijk2 = ijk1 + (ttm%ntcell(1)+2)*(ttm%ntcell(2)+2)
              ttm%eltemp(ijk2,ii,jj,ttm%eltcell(3)) = ttm%eltemp(ijk1,ii,jj,ttm%eltcell(3))
            End Do
          End Do
        End Do
      End Do
    End If

  End Select

End Subroutine boundaryCond

Subroutine uniformDist(lat_in,ttm)

! implement constant (homogeneous) spatial deposition

  Type( ttm_type ), Intent( InOut )   :: ttm 
  Real( Kind = wp ), Intent ( Inout ), Dimension(0:ttm%ntcell(1)+1,0:ttm%ntcell(2)+1,0:ttm%ntcell(3)+1) :: lat_in
  Real( Kind = wp ) :: dEdV

  ! express deposition energy per unit volume (eV/A^3):
  ! note penetration depth will be non-zero if laser is
  ! in use, otherwise use dE/dX value

  If (ttm%pdepth>zero_plus) Then
    dEdV = ttm%fluence/ttm%pdepth
  Else
    dEdV = ttm%dEdX/(Real(ttm%ntsys(1),Kind=wp)*Real(ttm%ntsys(2),Kind=wp)*ttm%delx*ttm%dely)
  End If

  ! homogeneous excitation: each temperature config%cell receives
  ! the same energy

  lat_in(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3))=dEdV*ttm%volume

End Subroutine uniformDist

Subroutine uniformDistZexp(lat_in,ttm)

! implement constant (homogeneous) spatial deposition
! in x and y-directions, exponential decay of fluence
! in z-direction (only with laser)

  Type( ttm_type ), Intent( InOut )   :: ttm 
  Real( Kind = wp ), Intent ( Inout ), Dimension(0:ttm%ntcell(1)+1,0:ttm%ntcell(2)+1,0:ttm%ntcell(3)+1) :: lat_in
  Real( Kind = wp ) :: dEdVmax, dEdV, zz, rpdepth
  Integer :: k

  ! express maximum deposition energy per unit volume (eV/A^3)

  If (ttm%pdepth>zero_plus) Then
    rpdepth = 1.0_wp/ttm%pdepth
  Else
    rpdepth = 0.0_wp
  End If
  dEdVmax = ttm%fluence*rpdepth

  ! loop through z-config%cells: calculate stopping power per
  ! config%cell based on z-position (maximum at z=0, grid centre)
  ! and assign across all x and y points in plane

  Do k = 1, ttm%ntcell(3)
    zz = Abs(Real(k+ttm%ntcelloff(3)-ttm%midI(3),Kind=wp))
    dEdV = dEdVmax*Exp(-zz*ttm%delz*rpdepth)
    lat_in(1:ttm%ntcell(1),1:ttm%ntcell(2),k) = dEdV*ttm%volume
  End Do

End Subroutine uniformDistZexp

Subroutine gaussianTrack(lat_in, ttm,comm)

! implement gaussian spatial deposition

  Type( ttm_type ), Intent( InOut )   :: ttm 
  Real ( Kind = wp ), Intent ( Inout ), Dimension(0:ttm%ntcell(1)+1,0:ttm%ntcell(2)+1,0:ttm%ntcell(3)+1) :: lat_in
  Type( comms_type), Intent( InOut) :: comm
  Real ( Kind = wp ) :: normdEdX,realdEdx,sigmamx,sigmamy,sig2x,sig2y,sigcellx,sigcelly
  Real ( Kind = wp ) :: ii,jj,ii2,jj2,iip2,jjp2,iim2,jjm2
  Integer :: i,j,sgmx,sgmy
  Logical :: cutwarn=.false.

  lat_in(:,:,:) = 0.0_wp

  ! converting stopping power to a value per config%cell (in z-direction)

  normdEdX = ttm%dEdX*ttm%delz

  ! find extents of gaussian in x and y directions

  sigcellx = ttm%sig/ttm%delx
  sigcelly = ttm%sig/ttm%dely
  sig2x = 2.0_wp*sigcellx*sigcellx
  sig2y = 2.0_wp*sigcelly*sigcelly
  sigmamx = ttm%sigmax*sigcellx
  sigmamy = ttm%sigmax*sigcelly
  
  ! if cutoff larger than ionic temperature grid,
  ! warn of deposition errors

  If (sigmamx > Ceiling(ttm%ntsys(1)/2.0_wp)) Then
    sigmamx = Ceiling(ttm%ntsys(1)/2.0_wp)
    cutwarn = .true.
  End If
  If (sigmamy > Ceiling(ttm%ntsys(2)/2.0_wp)) Then
    sigmamy = Ceiling(ttm%ntsys(2)/2.0_wp)
    cutwarn = .true.
  End If

  If (comm%idnode == 0 .and. cutwarn) Then
    Call warning(535,0.0_wp,0.0_wp,0.0_wp)
  End If

  sgmx = Nint(sigmamx)
  sgmy = Nint(sigmamy)

  ! apply five-point linear stencil for gaussian track:
  ! stencil modified to (hopefully!) deposit correct overall energy

  Do j=1,ttm%ntcell(2)
    jj = Real(j+ttm%ntcelloff(2)-ttm%midI(2),Kind=wp)
    jj2 = -jj*jj/sig2y
    jjp2 = -(jj+0.5_wp)*(jj+0.5_wp)/sig2y
    jjm2 = -(jj-0.5_wp)*(jj-0.5_wp)/sig2y
    Do i=1,ttm%ntcell(1)
      ii = Real(i+ttm%ntcelloff(1)-ttm%midI(1),Kind=wp)
      ii2 = -ii*ii/sig2x
      iip2 = -(ii+0.5_wp)*(ii+0.5_wp)/sig2x
      iim2 = -(ii-0.5_wp)*(ii-0.5_wp)/sig2x
      If (Abs(ii)<=sgmx .and. Abs(jj)<=sgmy) Then
        lat_in(i,j,1:ttm%ntcell(3)) = 0.2_wp*normdEdX/(2.0_wp*pi*sigcellx*sigcelly)*&
                                  (Exp(ii2+jj2)+Exp(iim2+jj2)+Exp(iim2+jj2)+Exp(ii2+jjp2)+Exp(ii2+jjm2))
      End If
    End Do
  End Do

  ! calculate deposited energy for comparison with specified value
  ! (note that stopping power is in z-direction)

  realdEdx = Sum(lat_in(1:ttm%ntcell(1),1:ttm%ntcell(2),1:ttm%ntcell(3)))
  Call gsum(comm,realdEdx)
  realdEdx = realdEdx/(Real(ttm%ntsys(3),Kind=wp)*ttm%delz)

  ! check if lattice sum equals the expected value

  If (comm%idnode == 0 .and. Abs((realdEdx-ttm%dEdX)/ttm%dEdX) > 0.01_wp) Then
    Call warning(540,Abs(realdEdx-ttm%dEdX)/ttm%dEdX*100_wp,0.0_wp,0.0_wp)
  End If
End Subroutine gaussianTrack
End Module ttm
