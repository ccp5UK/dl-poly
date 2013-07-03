Subroutine read_config &
           (megatm,levcfg,imcon,l_ind,l_str,rcut,dvar,xhi,yhi,zhi,dens0,dens)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading CONFIG and getting the average
! particle density
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,   Only : nconf,nrite,half_minus
  Use config_module,  Only : imc_n,cell,allocate_config_arrays_read, &
                             natms,nlast,atmnam,lsi,lsa,ltg, &
                             xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use domains_module
  Use parse_module,   Only : tabs_2_blanks, &
                             strip_blanks, get_word, word_2_real
  Use io_module,      Only : io_set_parameters,io_get_parameters,     &
                             io_init, io_open, io_close, io_finalize, &
                             IO_READ_MASTER, IO_READ_NETCDF

  Implicit None

  Integer,           Intent( In    ) :: megatm,levcfg,imcon
  Logical,           Intent( In    ) :: l_ind,l_str
  Real( Kind = wp ), Intent( In    ) :: rcut,dvar
  Real( Kind = wp ), Intent( InOut ) :: xhi,yhi,zhi
  Real( Kind = wp ), Intent(   Out ) :: dens0,dens

  Real( Kind = wp ) :: cut

  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,fname
  Logical                :: safe  = .true.  , &
                            l_his = .false. , &
                            l_xtr = .false. , &
                            fast
  Integer                :: fail(1:4),i,j,idm,icell,ncells,   &
                            mxatms,indatm,nattot,totatm,      &
                            ipx,ipy,ipz,nlx,nly,nlz,          &
                            ix,iy,iz,jx,jy,jz,read_buffer_size
  Real( Kind = wp )      :: celprp(1:10),rcell(1:9),celh(1:9),det, &
                            volm,vcell,                            &
                            sxx,syy,szz,xdc,ydc,zdc,               &
                            pda_max,pda_min,pda_ave,               &
                            pda_dom_max,pda_dom_min

! Some parameters and variables needed by io_module interfaces

  Integer                           :: recsz = 73 ! default record size
  Integer                           :: fh, io_read
  Integer( Kind = MPI_OFFSET_KIND ) :: top_skip

  Real( Kind = wp ),    Dimension( : ), Allocatable :: pda

  Character( Len = 8 ), Dimension( : ), Allocatable :: chbuf
  Integer,              Dimension( : ), Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ), Allocatable :: axx,ayy,azz, &
                                                       bxx,byy,bzz, &
                                                       cxx,cyy,czz


! image conditions not compliant with DD and link-cell

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Real space cutoff shortened by 50% but not < 1 Angstrom
!(=rcut_def in scan_control)

  cut=Max(0.5_wp*rcut,1.0_wp)+1.0e-6_wp

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)
  volm = celprp(10)

! Calculate the number of link-cells per domain in every direction

  nlx=Int(celprp(7)/(cut*nprx_r))
  nly=Int(celprp(8)/(cut*npry_r))
  nlz=Int(celprp(9)/(cut*nprz_r))

  ncells=nlx*nly*nlz

! Check for link cell algorithm violations

  If (ncells == 0) Call error(307)

! Amend volume of density cell if cluster, slab or bulk slab
! cell dimensional properties overwritten but not needed anyway

  If (imcon == 0 .or. imcon == 6 .or. imc_n == 6) Then
     celh=cell

     If (imcon == 0) Then
        celh(1) = Max(1.0_wp,xhi)
        celh(5) = Max(1.0_wp,yhi)
        celh(9) = Max(1.0_wp,zhi)
     Else If (imcon == 6) Then
        celh(9) = Max(1.0_wp,zhi)
     End If

     Call dcell(celh,celprp)
     volm = celprp(10)
  End If

  vcell = volm / (Real(ncells,wp) * Real(mxnode,wp))

! Approximate density and mxatms

  dens = Real(megatm,wp) / volm
  mxatms = Max(1 , Nint( (dvar**1.7_wp) * dens*vcell * Real((nlx+3)*(nly+3)*(nlz+3),wp)))

! Allocate necessary arrays to read CONFIG

  Call allocate_config_arrays_read(mxatms)

! Get type of I/O for reading

  Call io_get_parameters( user_method_read = io_read )

! Define filename ASCII or netCDF

  If (io_read /= IO_READ_NETCDF) Then
     fname='CONFIG'
  Else
     fname='CONFIG.nc'
  End If

! Define/Detect the FAST reading status

  If      (io_read == IO_READ_MASTER) Then

     fast = .false.

  Else If (io_read == IO_READ_NETCDF) Then

     fast = .true.

  Else

! Check if the system input file is a new style CONFIG:
! (i)  all lines are 72 ASCII characters long with
!      a UNIX carriage return as end of line;
! (ii) LINE2 has the particles total value
!      after values of levcfg and imcon.
! No fall back if users have mangled with further lines

     fast = .true.
     If (idnode == 0) Then

! Open CONFIG

        Open(Unit=nconf, File=fname)

! Read the CONFIG file header (TITLE record)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
           i = i + 1
           safe = .false.
           Read(Unit=nconf, Fmt='(a1)', Advance='No', IOStat=j, End=10) record(i:i)
           safe = .true.
           If (j < 0) Go To 10
        End Do
10      Continue
        fast = (fast .and. i == recsz)

! Read configuration level and image condition (RECORD2)

        i = 0 ! position counter
        j = 0 ! IOStat return
        record = ' '
        Do
           i = i + 1
           safe = .false.
           Read(Unit=nconf, Fmt='(a1)', Advance='No', IOStat=j, End=20) record(i:i)
           safe = .true.
           If (j < 0) Go To 20
        End Do
20      Continue
        fast = (fast .and. i == recsz)

! Read particles total value

        Call get_word(record,word) ; Call get_word(record,word)
        Call get_word(record,word) ; i=Nint(word_2_real(word,0.0_wp,l_str))
        fast = (fast .and. i == megatm)

     End If
     If (mxnode > 1) Then
        Call gsync()
        Call gcheck(safe,"enforce")
        Call gcheck(fast,"enforce")
     End If
     If (.not.safe) Go To 50

! Close CONFIG

     If (idnode == 0) Close(Unit=nconf)

  End If

  fail = 0

! If MASTER read

  If (io_read == IO_READ_MASTER) Then

     Call invert(cell,rcell,det)

! Open CONFIG and skip the header

     If (idnode == 0) Then
        Open(Unit=nconf, File=fname)

        Read(Unit=nconf, Fmt=*)    ! CONFIG file header (TITLE record)
        Read(Unit=nconf, Fmt=*)    ! configuration level and image condition

        If (imcon /= 0) Then
           Read(Unit=nconf, Fmt=*) ! cell vectors (not defined for imcon=0) but cell
           Read(Unit=nconf, Fmt=*) ! is modified in set_bounds for imcon 0 and 6!!!
           Read(Unit=nconf, Fmt=*)
        End If
     End If

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_config allocation failure, node: ', idnode
        Call error(0)
     End If

! Initialise domain localised atom counter (config_module)
! and dispatched atom counter

     natms =0
     indatm=0
     Do nattot=1,megatm
        indatm=indatm+1

! Initialise transmission arrays

        chbuf(indatm)=' '
        iwrk(indatm)=0

        axx(indatm)=0.0_wp
        ayy(indatm)=0.0_wp
        azz(indatm)=0.0_wp

        If (levcfg > 0) Then
           bxx(indatm)=0.0_wp
           byy(indatm)=0.0_wp
           bzz(indatm)=0.0_wp

           If (levcfg > 1) Then
              cxx(indatm)=0.0_wp
              cyy(indatm)=0.0_wp
              czz(indatm)=0.0_wp
           End If
        End If

! Read in transmission arrays

        If (idnode == 0 .and. safe) Then
           record=' '
           Read(Unit=nconf, Fmt='(a)', End=30) record
           Call tabs_2_blanks(record) ; Call strip_blanks(record)
           Call get_word(record,word) ; chbuf(indatm)=word(1:8)
           If (l_ind) Then
              Call get_word(record,word)
              iwrk(indatm)=Nint(word_2_real(word,0.0_wp,l_str))
              If (iwrk(indatm) /= 0) Then
                 iwrk(indatm)=Abs(iwrk(indatm))
              Else
                 iwrk(indatm)=nattot
              End If
           Else
              iwrk(indatm)=nattot
           End If

           Read(Unit=nconf, Fmt=*, End=30) axx(indatm),ayy(indatm),azz(indatm)

           If (levcfg > 0) Then
              Read(Unit=nconf, Fmt=*, End=30) bxx(indatm),byy(indatm),bzz(indatm)
              If (levcfg > 1) Read(Unit=nconf, Fmt=*, End=30) cxx(indatm),cyy(indatm),czz(indatm)
           End If
           Go To 40

30         Continue
           safe=.false. ! catch error

40         Continue
        End If

! Circulate configuration data to all nodes when transmission arrays are filled up

        If (indatm == mxatms .or. nattot == megatm) Then

! Check if batch was read fine

           If (mxnode > 1) Call gcheck(safe)
           If (.not.safe) Go To 50

! Ensure all atoms are in prescribed simulation cell (DD bound) and broadcast them

           Call pbcshift(imcon,cell,indatm,axx,ayy,azz)

           If (mxnode > 1) Then
              Call MPI_BCAST(chbuf,indatm*8,MPI_CHARACTER,0,dlp_comm_world,ierr)
              Call MPI_BCAST(iwrk,indatm,MPI_INTEGER,0,dlp_comm_world,ierr)

              Call MPI_BCAST(axx,indatm,wp_mpi,0,dlp_comm_world,ierr)
              Call MPI_BCAST(ayy,indatm,wp_mpi,0,dlp_comm_world,ierr)
              Call MPI_BCAST(azz,indatm,wp_mpi,0,dlp_comm_world,ierr)

              If (levcfg > 0) Then
                 Call MPI_BCAST(bxx,indatm,wp_mpi,0,dlp_comm_world,ierr)
                 Call MPI_BCAST(byy,indatm,wp_mpi,0,dlp_comm_world,ierr)
                 Call MPI_BCAST(bzz,indatm,wp_mpi,0,dlp_comm_world,ierr)

                 If (levcfg > 1) Then
                    Call MPI_BCAST(cxx,indatm,wp_mpi,0,dlp_comm_world,ierr)
                    Call MPI_BCAST(cyy,indatm,wp_mpi,0,dlp_comm_world,ierr)
                    Call MPI_BCAST(czz,indatm,wp_mpi,0,dlp_comm_world,ierr)
                 End If
              End If
           End If

! Assign atoms to correct domains

           Do i=1,indatm
              sxx=rcell(1)*axx(i)+rcell(4)*ayy(i)+rcell(7)*azz(i)
              syy=rcell(2)*axx(i)+rcell(5)*ayy(i)+rcell(8)*azz(i)
              szz=rcell(3)*axx(i)+rcell(6)*ayy(i)+rcell(9)*azz(i)

! sxx,syy,szz are in [-0.5,0.5) interval as values as 0.4(9) may pose a problem

              sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
              syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
              szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

! fold back coordinates

              axx(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
              ayy(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
              azz(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz

! assign domain coordinates (call for errors)

              ipx=Int((sxx+0.5_wp)*nprx_r)
              ipy=Int((syy+0.5_wp)*npry_r)
              ipz=Int((szz+0.5_wp)*nprz_r)

! check for errors

              If (ipx == nprx .or. ipy == npry .or. ipz == nprz) Call error(513)
! assign domain

              idm=ipx+nprx*(ipy+npry*ipz)
              If (idm == idnode) Then
                 natms=natms+1

                 If (natms < mxatms) Then
                    atmnam(natms)=chbuf(i)
                    ltg(natms)=iwrk(i)

                    xxx(natms)=axx(i)
                    yyy(natms)=ayy(i)
                    zzz(natms)=azz(i)

                    If (levcfg > 0) Then
                       vxx(natms)=bxx(i)
                       vyy(natms)=byy(i)
                       vzz(natms)=bzz(i)
                    Else
                       vxx(natms)=0.0_wp
                       vyy(natms)=0.0_wp
                       vzz(natms)=0.0_wp
                    End If

                    If (levcfg > 1) Then
                       fxx(natms)=cxx(i)
                       fyy(natms)=cyy(i)
                       fzz(natms)=czz(i)
                    Else
                       fxx(natms)=0.0_wp
                       fyy(natms)=0.0_wp
                       fzz(natms)=0.0_wp
                    End If
                 Else
                    safe=.false.
                 End If
              End If
           End Do

! Check if all is dispatched fine

           If (mxnode > 1) Call gcheck(safe)
           If (.not.safe) Call error(45)

! Nullify dispatch counter

           indatm=0

        End If
     End Do

! Close CONFIG

     If (idnode == 0) Close(Unit=nconf)
     If (mxnode > 1) Call gsync()

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_config deallocation failure, node: ', idnode
        Call error(0)
     End If

! If PROPER read

  Else

! Open CONFIG

     If (fast) Then
        Call io_set_parameters( user_comm = dlp_comm_world )
        Call io_init( recsz )
        Call io_open( io_read, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )
     Else
        Open(Unit=nconf, File=fname)
     End If

! top_skip is header size

     If (io_read /= IO_READ_NETCDF) Then
        If (imcon == 0) Then
           top_skip = Int(2,MPI_OFFSET_KIND)
        Else
           top_skip = Int(5,MPI_OFFSET_KIND)
        End If
     Else
        top_skip = Int(1,MPI_OFFSET_KIND) ! This is now the frame = 1
     End If

     read_buffer_size = mxatms
     Call read_config_parallel                                     &
           (levcfg, imcon, l_ind, l_str, megatm, read_buffer_size, &
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi)

! Close CONFIG

     If (fast) Then
        Call io_close( fh )
        Call io_finalize
     Else
        Close(Unit=nconf)
     End If

  End If

! To prevent users from the danger of changing the order of calls
! in dl_poly set 'nlast' to the innocent 'natms'

  nlast=natms

! First check: Does the number of atoms in the system (MD cell) derived by
! topology description (FIELD) match the crystallographic (CONFIG) one?
! Check number of atoms in system (CONFIG = FIELD)

  totatm=natms
  If (mxnode > 1) Call gsum(totatm)
  If (totatm /= megatm) Call error(58)

! Record global atom indices for local sorting (config_module)

  Do i=1,natms
     lsi(i)=i
     lsa(i)=ltg(i)
  End Do
  Call shellsort2(natms,lsi,lsa)

  If (io_read /= IO_READ_MASTER) Then

! This section is not strictly necessary.  However, the new read in method
! means the atoms are not necessarily in the same order in memory as the
! older, slower, method would put them.  This bit makes sure that the order
! is so that 'ltg' is strictly monotonically increasing.  This captures the
! common case where CONFIG has the 'ltg' values all in order (or not
! specified), but there is no easy way for the general case of arbitrary
! ordering of the 'ltg' values in CONFIG.  Of course, this makes no
! difference to the science and to restarts.  However, for initial runs it
! means the initial velocities will not be the same as the old method for
! the arbitrary ordering case.

     atmnam( 1:natms ) = atmnam( lsi( 1:natms ) )
     ltg( 1:natms ) = ltg( lsi( 1:natms ) )

     xxx( 1:natms ) = xxx( lsi( 1:natms ) )
     yyy( 1:natms ) = yyy( lsi( 1:natms ) )
     zzz( 1:natms ) = zzz( lsi( 1:natms ) )

     If (levcfg > 0) Then
        vxx( 1:natms ) = vxx( lsi( 1:natms ) )
        vyy( 1:natms ) = vyy( lsi( 1:natms ) )
        vzz( 1:natms ) = vzz( lsi( 1:natms ) )

        If (levcfg > 1) Then
           fxx( 1:natms ) = fxx( lsi( 1:natms ) )
           fyy( 1:natms ) = fyy( lsi( 1:natms ) )
           fzz( 1:natms ) = fzz( lsi( 1:natms ) )
        End If
     End If
     Do i=1,natms
        lsi(i)=i
        lsa(i)=ltg(i)
     End Do

  End If

! READ CONFIG END

! PARTICLE DENSITY START
! Allocate and initialise particle density array

  fail=0
  Allocate (pda(1:ncells), Stat=fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'read_config allocation failure, node: ', idnode
     Call error(0)
  End If
  pda=0.0_wp

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! Shifts from global to local link-cell space:
! (0,0,0) left-most link-cell on the domain (halo)
! (nlx+2*nlp-1,nly+2*nlp-1,nly+2*nlp-1) right-most
! link-cell on the domain (halo)

  jx=1-nlx*idx
  jy=1-nly*idy
  jz=1-nlz*idz

! Get the inverse cell matrix

  Call invert(cell,rcell,celprp(10))

  Do i=1,natms
     sxx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     syy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     szz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! Get cell coordinates accordingly

     ix = Int(xdc*(sxx+0.5_wp)) + jx
     iy = Int(ydc*(syy+0.5_wp)) + jy
     iz = Int(zdc*(szz+0.5_wp)) + jz

! Put all particles in bounded link-cell space: lower and upper
! bounds as 1 <= i_coordinate <= nl_coordinate

     ix = Max( Min( ix , nlx) , 1)
     iy = Max( Min( iy , nly) , 1)
     iz = Max( Min( iz , nlz) , 1)

! Hypercube function transformation (counting starts from one
! rather than zero /map_domains/

     icell=1+(ix-1)+nlx*((iy-1)+nly*(iz-1))

     pda(icell)=pda(icell)+1.0_wp
  End Do

  pda_max=  0.0_wp
  pda_min=100.0_wp
  pda_ave=  0.0_wp
  Do icell=1,ncells
     pda_max=Max(pda(icell),pda_max)
     pda_min=Min(pda(icell),pda_min)
     pda_ave=pda_ave+pda(icell)
  End Do
  pda_ave=pda_ave/Real(ncells,wp)

  pda_dom_max=pda_ave
  pda_dom_min=pda_ave
  If (mxnode > 1) Then
     Call gmax(pda_dom_max)
     Call gmin(pda_dom_min)

     Call gmax(pda_max)
     Call gmin(pda_min)

     Call gsum(pda_ave)
     pda_ave=pda_ave/Real(mxnode,wp)
  End If

! Approximation for maximum global density by
! the inter-domain imbalance of domain density

  dens0=pda_max/vcell
  If (mxnode > 1) Then
     If (Nint(pda_dom_min) == 0) Then
        dens = dens0 ! domain(s) matched on vacuum (take no risk)
     Else If (1.15_wp*pda_dom_min > pda_dom_max) Then
        dens = pda_ave/vcell
     Else If (1.25_wp*pda_dom_min > pda_dom_max) Then
        dens = pda_dom_max/vcell
     Else
        dens = dens0 ! too big an imbalance (take no risk)
     End If
  Else
     dens = 1.25_wp*pda_ave/vcell ! allow 25% imbalance
  End If

  Deallocate (pda, Stat=fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'read_config deallocation failure, node: ', idnode
     Call error(0)
  End If

! PARTICLE DENSITY END

  Return

! error exit for CONFIG file read

50 Continue
  If (idnode == 0) Close(Unit=nconf)
  Call error(55)

End Subroutine read_config
