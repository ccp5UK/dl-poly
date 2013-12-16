Subroutine read_config_parallel                  &
           (levcfg, imcon, l_ind, l_str, megatm, &
            l_his, l_xtr, fast, fh, top_skip, xhi, yhi, zhi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the CONFIG data file in parallel
!
! copyright - daresbury laboratory
! author    - i.j.bush & i.t.todorov december 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,   Only : nrite, nconf, half_minus
  Use domains_module, Only : nprx,npry,nprz,nprx_r,npry_r,nprz_r
  Use config_module,  Only : cell,natms,atmnam,ltg, &
                             xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use parse_module,   Only : strip_blanks, get_word, word_2_real
  Use io_module,      Only : io_get_parameters, io_read_batch, &
                             io_nc_get_dim, io_nc_get_var, IO_READ_NETCDF

  Implicit None

  Logical,                           Intent( In    ) :: l_ind,l_str,l_his,fast,l_xtr
  Integer,                           Intent( In    ) :: levcfg,imcon,megatm,fh
  Integer( Kind = MPI_OFFSET_KIND ), Intent( In    ) :: top_skip
  Real( Kind = wp ),                 Intent(   Out ) :: xhi,yhi,zhi

  Logical                :: safe,do_read
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,forma
  Integer                :: fail(1:8),i,j,k,                &
                            idm,ipx,ipy,ipz,indatm,         &
                            n_read_procs_use,per_read_proc, &
                            my_read_proc_num,ats_per_proc,  &
                            recs_per_at,recs_per_proc,      &
                            wp_vals_per_at,n_loc,           &
                            to_read,which_read_proc,this_base_proc
  Integer( Kind = ip )   :: n_sk,n_ii,n_jj
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz

! Some parameters and variables needed by io_module interfaces

  Integer                           :: io_read
  Integer                           :: recsz, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io, n_skip
  Integer                           :: this_rec_buff, recs_to_read
  Integer                           :: n_ats_in_file

! netCDF

  Integer :: frame, start(1:3), count(1:3)

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx_read,ayy_read,azz_read, &
                                                          bxx_read,byy_read,bzz_read, &
                                                          cxx_read,cyy_read,czz_read

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf_read,chbuf_scat
  Integer,              Dimension( : ),    Allocatable :: iwrk_read,iwrk_scat

  Integer,              Dimension( : ),    Allocatable :: first_at,orig_first_at
  Integer,              Dimension( : ),    Allocatable :: n_held,where_buff
  Integer,              Dimension( : ),    Allocatable :: owner_read

  Real( Kind = wp ),    Dimension( :, : ), Allocatable :: scatter_buffer_read
  Real( Kind = wp ),    Dimension( :, : ), Allocatable :: scatter_buffer

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: rec_buff


! Get reading method, total number of I/O heads and buffer size

  Call io_get_parameters( user_method_read      = io_read          )
  Call io_get_parameters( user_n_io_procs_read  = n_read_procs_use )
  Call io_get_parameters( user_buffer_size_read = batsz            )

  fail = 0 ! fail initialisation

  If (levcfg /= 3) Then
     wp_vals_per_at = 3 * (1+levcfg) ! Scatter buffer sizes
     recs_per_at    = 2 + levcfg     ! Scatter buffer sizes
  Else
     wp_vals_per_at = 3 ! Scatter buffer sizes
     recs_per_at    = 1 ! Scatter buffer sizes
  End If

! Note: make 'first_at' and 'orig_first_at' 1 element bigger than strictly
! required to make checking at the end of reading much easier and clearer

  Allocate(first_at(0:n_read_procs_use),orig_first_at(0:n_read_procs_use), Stat=fail(1))
  Allocate(chbuf(1:batsz),iwrk(1:batsz),                                   Stat=fail(2))
  Allocate(scatter_buffer(1:wp_vals_per_at,1:batsz),                       Stat=fail(3))
  If (Any(fail(1:3) > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'read_config_parallel allocation failure 1, node: ', idnode
     Call error(0)
  End If

! define basic quantities for the parallel ASCII reading

  per_read_proc = mxnode / n_read_procs_use
  do_read = (Mod( idnode, per_read_proc ) == 0 .and. idnode < per_read_proc * n_read_procs_use)
  my_read_proc_num = idnode / per_read_proc

! Note 'first_at' and 'orig_first_at' have one more element
! in the array than strictly required - makes it easier to
! check that reading by the last I/O processor has finished

  ats_per_proc = megatm / n_read_procs_use
  Do i=0,n_read_procs_use
     first_at(i) = i*ats_per_proc + Min(i,megatm-ats_per_proc*n_read_procs_use)
  End Do
  orig_first_at = first_at
  ats_per_proc = Max(1,ats_per_proc) ! Fix it if 0
  recs_per_proc = ats_per_proc * recs_per_at

! Catch the case where the first atom belonging to
! a read processor does not actually exists - i.e.
! I/O procs count > megatm, and limit reading by do_read

  If (my_read_proc_num < n_read_procs_use) &
     do_read = (do_read .and. first_at(my_read_proc_num) < megatm)

! Skip to the point of reading

  If (do_read) Then

     n_skip = Int(recs_per_at,MPI_OFFSET_KIND) * Int(first_at(my_read_proc_num),MPI_OFFSET_KIND) + &
              top_skip-Int(1,MPI_OFFSET_KIND)

     If (.not.fast) Then
        n_sk=Int(n_skip,ip)
        n_jj=73*batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
           Do n_ii=1_ip,n_sk/n_jj
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_jj
              Read(Unit=nconf, Fmt=forma, End=100)
           End Do
           n_ii=Mod(Int(n_skip,ip),n_jj)
           If (n_ii > 0_ip) Then
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_ii
              Read(Unit=nconf, Fmt=forma, End=100)
           End If
        Else
           forma=' '
           Write(forma,'( "(", i0, "/)" )') n_sk
           Read(Unit=nconf, Fmt=forma, End=100)
        End If

        recsz=200
        forma=' '
        Write(forma,'( "(", i0, "a1)" )') recsz
     Else
        rec_mpi_io = n_skip + Int(1,MPI_OFFSET_KIND)
        recsz=73
        If (levcfg == 3) recsz = 35
     End If

! Allocate record buffer, reading buffers, scatter buffers and indexing arrays

     If (io_read /= IO_READ_NETCDF) Then
        Allocate(rec_buff(1:recsz,1:batsz),                                  Stat=fail(1))
     Else
        Allocate(rec_buff(1:Len( chbuf_read ),1:batsz),                      Stat=fail(1))
     End If
     Allocate(chbuf_read(1:batsz),iwrk_read(1:batsz),                        Stat=fail(2))
     Allocate(axx_read(1:batsz),ayy_read(1:batsz),azz_read(1:batsz),         Stat=fail(3))
     Allocate(bxx_read(1:batsz),byy_read(1:batsz),bzz_read(1:batsz),         Stat=fail(4))
     Allocate(cxx_read(1:batsz),cyy_read(1:batsz),czz_read(1:batsz),         Stat=fail(5))
     Allocate(scatter_buffer_read(1:wp_vals_per_at,1:batsz),                 Stat=fail(6))
     Allocate(chbuf_scat(1:batsz),iwrk_scat(1:batsz),                        Stat=fail(7))
     Allocate(n_held(0:mxnode-1),where_buff(0:mxnode-1),owner_read(1:batsz), Stat=fail(8))
     If (Any(fail(1:8) > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_config_parallel allocation failure 2, node: ', idnode
        Call error(0)
     End If

  Else

! It is Illegal to pass unallocated allocatable arrays to routines.
! Therefore for arrays that are used by the mpi_scatterv calls
! below allocate them to zero size if they are not used on this core

     Allocate(scatter_buffer_read(1:0,1:0),   Stat=fail(1))
     Allocate(chbuf_scat(1:0),iwrk_scat(1:0), Stat=fail(2))
     Allocate(n_held(0:-1),where_buff(0:-1),  Stat=fail(3))
     If (Any(fail(1:3) > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_config_parallel allocation failure 3, node: ', idnode
        Call error(0)
     End If

  End If

! Initialise extreme box dimensions

  xhi = 0.0_wp
  yhi = 0.0_wp
  zhi = 0.0_wp

  If (.not.l_xtr) Call invert(cell,rcell,det)

! Initialise domain localised atom counter (config_module),
! dispatched atom counter and safe dispatch flag

  natms =0
  indatm=0
  safe  =.true.

  Do k=1,megatm

! Read in transmission arrays

     Readers_only: If (do_read .and. indatm == 0) Then
        to_read = Min(batsz,orig_first_at(my_read_proc_num+1)-first_at(my_read_proc_num))

        No_netCDF: If (io_read /= IO_READ_NETCDF) Then

           this_rec_buff = 0
           recs_to_read  = 0
           Do i=1,to_read
              If (this_rec_buff == 0) Then
                 recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at )
                 If (.not.fast) Then
                    Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                 End If
              End If

! Atom details

              this_rec_buff = this_rec_buff + 1
              Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                 record( j:j ) = rec_buff( j, this_rec_buff )
              End Do
              Call strip_blanks(record)
              Call get_word(record,word) ; chbuf_read(i)=word(1:8)
              If (l_ind) Then
                 Call get_word(record,word)
                 iwrk_read(i)=Nint(word_2_real(word,0.0_wp,l_str))
                 If (iwrk_read(i) /= 0) Then
                    iwrk_read(i)=Abs(iwrk_read(i))
                 Else
                    iwrk_read(i)=first_at(my_read_proc_num)+i
                 End If
              Else
                 iwrk_read(i)=first_at(my_read_proc_num)+i
              End If

              If (levcfg == 3) Read(record, Fmt=*, End=100) axx_read(i),ayy_read(i),azz_read(i)

              If (this_rec_buff == recs_to_read) Then
                 this_rec_buff = 0
                 recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 1 )
                 If (.not.fast) Then
                    Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                 End If
              End If

              If (levcfg /= 3) Then

! Positions

                 this_rec_buff = this_rec_buff + 1
                 Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                    record( j:j ) = rec_buff( j, this_rec_buff )
                 End Do
                 Read(record, Fmt=*, End=100) axx_read(i),ayy_read(i),azz_read(i)
                 If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
                    If (levcfg > 0) Then
                       recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 2 )
                       If (.not.fast) Then
                          Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                       Else
                          Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                          rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                       End If
                    End If
                 End If

! Velocities

                 If (levcfg > 0) Then
                    this_rec_buff = this_rec_buff + 1
                    Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                       record( j:j ) = rec_buff( j, this_rec_buff )
                    End Do
                    Read(record, Fmt=*, End=100) bxx_read(i),byy_read(i),bzz_read(i)
                    If (this_rec_buff == recs_to_read) Then
                       this_rec_buff = 0
                       If (levcfg > 1) Then
                          recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 3 )
                          If (.not.fast) Then
                             Read(Unit=nconf, Fmt=forma) rec_buff( :, 1:recs_to_read )
                          Else
                             Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                             rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                          End If
                       End If
                    End If

! Forces

                    If (levcfg > 1) Then
                       this_rec_buff = this_rec_buff + 1
                       Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                          record( j:j ) = rec_buff( j, this_rec_buff )
                       End Do
                       Read(record, Fmt=*, End=100) cxx_read(i),cyy_read(i),czz_read(i)
                       If (this_rec_buff == recs_to_read) Then
                          this_rec_buff = 0
                       End If
                    End If
                 End If

              End If
           End Do

        Else

           If (to_read /= 0) Then
              frame = Int(top_skip,Kind(frame))

              Call io_nc_get_var( 'atomnames', fh, rec_buff, (/ first_at( my_read_proc_num ) + 1, frame /), (/ 8, to_read, 1 /) )
              Do i = 1, to_read
                 Do j = 1, Min( Len( chbuf_read ), Size( rec_buff, Dim = 1 ) )
                    chbuf_read( i )( j:j ) = rec_buff( j, i )
                 End Do
              End Do
              If (l_ind) Then
                 Call io_nc_get_var( 'indices', fh, iwrk_read , (/ first_at( my_read_proc_num ) + 1, frame /), (/ to_read, 1 /) )
              End If

              start = (/ 1, first_at( my_read_proc_num ) + 1, frame /)
              count = (/ 3, to_read, 1 /)

              Select Case( levcfg )
              Case( 0, 3 )
                 Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
              Case( 1 )
                 Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
                 Call get_var( 'velocities' , fh, start, count, bxx_read, byy_read, bzz_read )
              Case( 2 )
                 Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
                 Call get_var( 'velocities' , fh, start, count, bxx_read, byy_read, bzz_read )
                 Call get_var( 'forces'     , fh, start, count, cxx_read, cyy_read, czz_read )
              End Select
           End If

        End If No_netCDF

        If (.not.l_xtr) Then

! Ensure all atoms are in prescribed simulation cell (DD bound)
!
!           Call pbcshift(imcon,cell,to_read,axx_read,ayy_read,azz_read)

           n_held=0
           Do i=1,to_read
              sxx=rcell(1)*axx_read(i)+rcell(4)*ayy_read(i)+rcell(7)*azz_read(i)
              syy=rcell(2)*axx_read(i)+rcell(5)*ayy_read(i)+rcell(8)*azz_read(i)
              szz=rcell(3)*axx_read(i)+rcell(6)*ayy_read(i)+rcell(9)*azz_read(i)

! sxx,syy,szz are in [-0.5,0.5) interval as values as 0.4(9) may pose a problem

              sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
              syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
              szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

! fold back coordinates

              axx_read(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
              ayy_read(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
              azz_read(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz

! assign domain coordinates (call for errors)

              ipx=Int((sxx+0.5_wp)*nprx_r)
              ipy=Int((syy+0.5_wp)*npry_r)
              ipz=Int((szz+0.5_wp)*nprz_r)

              idm=ipx+nprx*(ipy+npry*ipz)
              If (idm < 0 .or. idm > (mxnode-1)) Call error(513)
              owner_read(i) = idm
              n_held(idm) = n_held(idm)+1
           End Do

           where_buff(0)=0
           Do i=1,mxnode-1
              where_buff(i) = where_buff(i-1) + n_held(i-1)
           End Do

           Do i=1,to_read
              idm = where_buff(owner_read(i))
              idm = idm+1
              where_buff(owner_read(i)) = idm

              chbuf_scat(idm) = chbuf_read(i)
              iwrk_scat(idm)  = iwrk_read(i)

              scatter_buffer_read(1,idm) = axx_read(i)
              scatter_buffer_read(2,idm) = ayy_read(i)
              scatter_buffer_read(3,idm) = azz_read(i)

              If (levcfg /= 3) Then
                 If (levcfg > 0) Then
                    scatter_buffer_read(4,idm) = bxx_read(i)
                    scatter_buffer_read(5,idm) = byy_read(i)
                    scatter_buffer_read(6,idm) = bzz_read(i)

                    If (levcfg > 1) Then
                       scatter_buffer_read(7,idm) = cxx_read(i)
                       scatter_buffer_read(8,idm) = cyy_read(i)
                       scatter_buffer_read(9,idm) = czz_read(i)
                    End If
                 End If
              End If
           End Do

! If only detecting box dimensions for imcon == 0 or 6

        Else

! Get extremes

           xhi = Max( xhi, Maxval( Abs( axx_read( 1:to_read ) ) ) )
           yhi = Max( yhi, Maxval( Abs( ayy_read( 1:to_read ) ) ) )
           zhi = Max( zhi, Maxval( Abs( azz_read( 1:to_read ) ) ) )

        End If
     End If Readers_only

! Increase buffer counter and update first_at for
! the readers that have something left to read

     indatm = indatm+1
     If (do_read) Then
        If (first_at(my_read_proc_num) < first_at(my_read_proc_num+1)) &
             first_at(my_read_proc_num) = first_at(my_read_proc_num)+1
     End If

! Circulate configuration data to all nodes when transmission arrays are filled up
! Check against megatm since at low processors counts (i.e. 1) batsz can be > megatm

     Reorganize_buffer: If (indatm == batsz .or. (indatm > 0 .and. k == megatm)) Then

        Extent_2: If (.not.l_xtr) Then

           Do which_read_proc = 0 , n_read_procs_use-1
              If (orig_first_at(which_read_proc) >= megatm) Exit ! for non-reading readers

              this_base_proc = which_read_proc * per_read_proc
              If (idnode == this_base_proc) Then
                 where_buff(0) = 0
                 Do i=1,mxnode-1
                    where_buff(i) = where_buff(i-1) + n_held(i-1)
                 End Do
              End If

              Call MPI_SCATTER( n_held, 1, MPI_INTEGER, n_loc, 1, MPI_INTEGER, this_base_proc, &
                   dlp_comm_world, ierr )

              Call MPI_SCATTERV( chbuf_scat, 8 * n_held, 8 * where_buff, MPI_CHARACTER, &
                                 chbuf     , 8 * n_loc ,                 MPI_CHARACTER, &
                                 this_base_proc, dlp_comm_world, ierr )

              Call MPI_SCATTERV( iwrk_scat ,     n_held,     where_buff, MPI_INTEGER, &
                                 iwrk      ,     n_loc ,                 MPI_INTEGER, &
                                 this_base_proc, dlp_comm_world, ierr )

              Call MPI_SCATTERV( scatter_buffer_read, wp_vals_per_at * n_held, wp_vals_per_at * where_buff, wp_mpi, &
                                 scatter_buffer     , wp_vals_per_at * n_loc ,                              wp_mpi, &
                                 this_base_proc, dlp_comm_world, ierr )

! Assign atoms to correct domains

Dispatch:     Do i=1,n_loc
                 natms=natms+1

! Check safety by the upper bound of: atmnam,ltg,xxx,yyy,zzz &
! possibly vxx,vyy,vzz & possibly fxx,fyy,fzz as guided by xxx

                 If (natms > Size( xxx )) Then
                    safe=.false.
                    Exit Dispatch
                 End If

                 atmnam(natms)=chbuf(i)
                 ltg(natms)=iwrk(i)

                 xxx(natms)=scatter_buffer(1,i)
                 yyy(natms)=scatter_buffer(2,i)
                 zzz(natms)=scatter_buffer(3,i)

                 If (levcfg /=3 ) Then
                    If (levcfg > 0) Then
                       vxx(natms)=scatter_buffer(4,i)
                       vyy(natms)=scatter_buffer(5,i)
                       vzz(natms)=scatter_buffer(6,i)
                    Else
                       vxx(natms)=0.0_wp
                       vyy(natms)=0.0_wp
                       vzz(natms)=0.0_wp
                    End If

                    If (levcfg > 1) Then
                       fxx(natms)=scatter_buffer(7,i)
                       fyy(natms)=scatter_buffer(8,i)
                       fzz(natms)=scatter_buffer(9,i)
                    Else
                       fxx(natms)=0.0_wp
                       fyy(natms)=0.0_wp
                       fzz(natms)=0.0_wp
                    End If
                 End If
              End Do Dispatch
           End Do

! Check if all is dispatched fine

           If (mxnode > 1) Call gcheck(safe)
           If (.not.safe) Call error(45)

        End If Extent_2

! Nullify dispatch counter

        indatm=0

     End If Reorganize_buffer

  End Do

! If only detecting box dimensions for imcon == 0 or 6

  If (l_xtr) Then
     Call gmax(xhi)
     Call gmax(yhi)
     Call gmax(zhi)
  End If

  If (l_his) Then

! Skip to the EoFrame of HISTORY when not fast

     If (do_read .and. (.not.fast)) Then
        n_skip = Int(recs_per_at,MPI_OFFSET_KIND) * Int(megatm-first_at(my_read_proc_num),MPI_OFFSET_KIND)

        n_sk=Int(n_skip,ip)
        n_jj=73*batsz ! Assuming average max line length of 73
        If (n_sk > n_jj) Then
           Do n_ii=1_ip,n_sk/n_jj
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_jj
              Read(Unit=nconf, Fmt=forma, End=100)
           End Do
           n_ii=Mod(Int(n_skip,ip),n_jj)
           If (n_ii > 0_ip) Then
              forma=' '
              Write(forma,'( "(", i0, "/)" )') n_ii
              Read(Unit=nconf, Fmt=forma, End=100)
           End If
        Else
           forma=' '
           Write(forma,'( "(", i0, "/)" )') n_sk
           Read(Unit=nconf, Fmt=forma, End=100)
        End If
     End If

  Else

     If (do_read) Then

        If ( io_read /= IO_READ_NETCDF) Then

! The last reader to check for EoFile in CONFIG
! and if none is hit to call error to abort

           If (first_at(my_read_proc_num) == megatm) Then
              recs_to_read = 1
              If (.not.fast) Then
                 Read(Unit=nconf, Fmt=forma, Iostat=ierr) rec_buff( :, 1:recs_to_read )
              Else
                 Call io_read_batch( fh, rec_mpi_io, 1, rec_buff, ierr )
              End If
              safe = (ierr /= 0)
           End If

        Else

! As netCDF files have no real concept of line numbers,
! instead check the arrays are the correct size

           Call io_nc_get_dim( 'atom', fh, n_ats_in_file )
           safe = n_ats_in_file == megatm

        End If

     End If
     If (mxnode > 1) Call gcheck(safe)
     If (.not.safe) Call error(58)

  End If

  If (do_read) Then
     Deallocate(rec_buff,                   Stat=fail(1))
     Deallocate(chbuf_read,iwrk_read,       Stat=fail(2))
     Deallocate(axx_read,ayy_read,azz_read, Stat=fail(3))
     Deallocate(bxx_read,byy_read,bzz_read, Stat=fail(4))
     Deallocate(cxx_read,cyy_read,czz_read, Stat=fail(5))
     Deallocate(owner_read,                 Stat=fail(6))
     If (Any(fail(1:6) > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'read_config_parallel deallocation failure 2, node: ', idnode
        Call error(0)
     End If
  End If

  Deallocate (first_at,orig_first_at, Stat=fail(1))
  Deallocate (n_held,where_buff,      Stat=fail(2))
  Deallocate (chbuf,chbuf_scat,       Stat=fail(3))
  Deallocate (iwrk,iwrk_scat,         Stat=fail(4))
  Deallocate (scatter_buffer_read,    Stat=fail(5))
  Deallocate (scatter_buffer,         Stat=fail(6))
  If (Any(fail(1:6) > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'read_config_parallel deallocation failure 1, node: ', idnode
     Call error(0)
  End If

  Return

! error exit for CONFIG file read

100 Continue
  Call error(55)

Contains

  Subroutine get_var( what, fh, start, count, x, y, z )

    Implicit None

    Character( Len = * )             , Intent( In    ) :: what
    Integer                          , Intent( In    ) :: fh
    Integer   ,        Dimension( : ), Intent( In    ) :: start
    Integer   ,        Dimension( : ), Intent( In    ) :: count
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: x
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: y
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: z

    Real( Kind = wp ), Dimension( :, : ), Allocatable :: buff

    Integer :: to_read
    Integer :: fail
    Integer :: i

    to_read = count( 2 )

    Allocate (buff( 1:3, 1:to_read ), Stat=fail)
    If (fail /= 0) Then
       Write( nrite, '(/,1x,a,i0)') 'read_config_parallel allocation failure 4, node: ', idnode
       Call error( 0 )
    End If

    Call io_nc_get_var( what, fh, buff, start, count )

    Do i = 1, to_read
       x( i ) = buff( 1, i )
       y( i ) = buff( 2, i )
       z( i ) = buff( 3, i )
    End Do

    Deallocate (buff, Stat=fail)
    If (fail /= 0) Then
       Write( nrite, '(/,1x,a,i0)') 'read_config_parallel allocation failure 4, node: ', idnode
       Call error( 0 )
    End If

  End Subroutine get_var

End Subroutine read_config_parallel
