Subroutine defects_reference_read_parallel              &
           (lvcfgr, imconr, celr, l_ind, l_str, megref, &
            fast, fh, top_skip, nrefs, namr, indr, xr, yr, zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the REFERENCE data file
! in parallel
!
! copyright - daresbury laboratory
! author    - i.j.bush & i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,   Only : nrite, nrefdt, half_minus,mxatms
  Use domains_module, Only : nprx,npry,nprz,nprx_r,npry_r,nprz_r
  Use parse_module,   Only : strip_blanks, get_word, word_2_real
  Use io_module,      Only : io_get_parameters, io_read_batch, &
                             io_nc_get_var, IO_READ_NETCDF

  Implicit None

  Logical,                           Intent( In    ) :: l_ind,l_str,fast
  Integer,                           Intent( In    ) :: lvcfgr,imconr,megref,fh
  Integer( Kind = MPI_OFFSET_KIND ), Intent( In    ) :: top_skip
  Real( Kind = wp ),                 Intent( In    ) :: celr(1:9)
  Character( Len = 8 ),              Intent(   Out ) :: namr(1:mxatms)
  Integer,                           Intent(   Out ) :: indr(1:mxatms),nrefs
  Real( Kind = wp ),                 Intent(   Out ) :: xr(1:mxatms),yr(1:mxatms),zr(1:mxatms)

  Logical                :: safe,do_read
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,forma
  Integer                :: fail(1:6),i,j,k,                &
                            idm,ipx,ipy,ipz,indatm,         &
                            n_read_procs_use,per_read_proc, &
                            my_read_proc_num,ats_per_proc,  &
                            recs_per_at,recs_per_proc,      &
                            wp_vals_per_at,n_loc,           &
                            to_read,which_read_proc,this_base_proc
  Real( Kind = wp )      :: rcell(1:9),det,sxx,syy,szz

! Some parameters and variables needed by io_module interfaces

  Integer                           :: io_read
  Integer                           :: recsz, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io, n_skip
  Integer                           :: this_rec_buff, recs_to_read

! netCDF

  Integer :: frame, start(1:3), count(1:3)

  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx_read,ayy_read,azz_read

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

  wp_vals_per_at = 3              ! Scatter buffer sizes
  recs_per_at    = 2 + lvcfgr     ! Scatter buffer sizes

! Note: make 'first_at' and 'orig_first_at' 1 element bigger than strictly
! required to make checking at the end of reading much easier and clearer

  Allocate(first_at(0:n_read_procs_use),orig_first_at(0:n_read_procs_use), Stat=fail(1))
  Allocate(chbuf(1:batsz),iwrk(1:batsz),                                   Stat=fail(2))
  Allocate(scatter_buffer(1:wp_vals_per_at,1:batsz),                       Stat=fail(3))
  If (Any(fail(1:3) > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_read_parallel allocation failure 1, node: ', idnode
     Call error(0)
  End If

! define basic quantities for the parallel ASCII reading

  per_read_proc = mxnode / n_read_procs_use
  do_read = (Mod( idnode, per_read_proc ) == 0 .and. idnode < per_read_proc * n_read_procs_use)
  my_read_proc_num = idnode / per_read_proc

! Note 'first_at' and 'orig_first_at' have one more element
! in the array than strictly required - makes it easier to
! check that reading by the last I/O processor has finished

  ats_per_proc = megref / n_read_procs_use
  Do i=0,n_read_procs_use
     first_at(i) = i*ats_per_proc + Min(i,megref-ats_per_proc*n_read_procs_use)
  End Do
  orig_first_at = first_at
  ats_per_proc = Max(1,ats_per_proc) ! Fix it if 0
  recs_per_proc = ats_per_proc * recs_per_at

! Catch the case where the first atom belonging to
! a read processor does not actually exists - i.e.
! I/O procs count > megref, and limit reading by do_read

  If (my_read_proc_num < n_read_procs_use) &
     do_read = (do_read .and. first_at(my_read_proc_num) < megref)

! Skip to the point of reading

  If (do_read) Then

     n_skip = Int(recs_per_at,MPI_OFFSET_KIND) * Int(first_at(my_read_proc_num),MPI_OFFSET_KIND) + &
              top_skip-Int(1,MPI_OFFSET_KIND)

     If (.not.fast) Then
        forma=' '
        Write(forma,'( "(", i0, "/)" )') Int(n_skip,ip)
        Read(Unit=nrefdt, Fmt=forma, End=100)
        recsz=200
        forma=' '
        Write(forma,'( "(", i0, "a1)" )') recsz
     Else
        rec_mpi_io = n_skip + Int(1,MPI_OFFSET_KIND)
        recsz=73
     End If

! Allocate record buffer, reading buffers, scatter buffers and indexing arrays

     If (io_read /= IO_READ_NETCDF) Then
        Allocate(rec_buff(1:recsz,1:batsz),                                  Stat=fail(1))
     Else
        Allocate(rec_buff(1:Len( chbuf_read ),1:batsz),                      Stat=fail(1))
     End If
     Allocate(chbuf_read(1:batsz),iwrk_read(1:batsz),                        Stat=fail(2))
     Allocate(axx_read(1:batsz),ayy_read(1:batsz),azz_read(1:batsz),         Stat=fail(3))
     Allocate(scatter_buffer_read(1:wp_vals_per_at,1:batsz),                 Stat=fail(4))
     Allocate(chbuf_scat(1:batsz),iwrk_scat(1:batsz),                        Stat=fail(5))
     Allocate(n_held(0:mxnode-1),where_buff(0:mxnode-1),owner_read(1:batsz), Stat=fail(6))
     If (Any(fail(1:6) > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_read_parallel allocation failure 2, node: ', idnode
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
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_read_parallel allocation failure 3, node: ', idnode
        Call error(0)
     End If

  End If

  Call invert(celr,rcell,det)

! Initialise domain localised atom counter (defects._module),
! dispatched atom counter and safe dispatch flag

  nrefs =0
  indatm=0
  safe  =.true.

  Do k=1,megref

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
                    Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
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

              If (this_rec_buff == recs_to_read) Then
                 this_rec_buff = 0
                 recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 1 )
                 If (.not.fast) Then
                    Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                 Else
                    Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                    rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                 End If
              End If

! Positions

              this_rec_buff = this_rec_buff + 1
              Do j = 1, Min( Len( record ), Size( rec_buff, Dim = 1 ) )
                 record( j:j ) = rec_buff( j, this_rec_buff )
              End Do
              Read(record, Fmt=*, End=100) axx_read(i),ayy_read(i),azz_read(i)
              If (this_rec_buff == recs_to_read) Then
                 this_rec_buff = 0
                 If (lvcfgr > 0) Then
                    recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 2 )
                    If (.not.fast) Then
                       Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                    Else
                       Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                       rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                    End If
                 End If
              End If

! Velocities

              If (lvcfgr > 0) Then
                 this_rec_buff = this_rec_buff + 1
                 If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
                    If (lvcfgr > 1) Then
                       recs_to_read = Min( Size( rec_buff, Dim = 2 ), ( to_read - i + 1 ) * recs_per_at - 3 )
                       If (.not.fast) Then
                          Read(Unit=nrefdt, Fmt=forma) rec_buff( :, 1:recs_to_read )
                       Else
                          Call io_read_batch( fh, rec_mpi_io, recs_to_read, rec_buff, ierr )
                          rec_mpi_io = rec_mpi_io + Int(recs_to_read,MPI_OFFSET_KIND)
                       End If
                    End If
                 End If
              End If

! Forces

              If (lvcfgr > 1) Then
                 this_rec_buff = this_rec_buff + 1

                 If (this_rec_buff == recs_to_read) Then
                    this_rec_buff = 0
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

              Call get_var( 'coordinates', fh, start, count, axx_read, ayy_read, azz_read )
           End If

        End If No_netCDF

! Assign atoms positions in fractional coordinates to the correct domains
! (DD bounding)

        Call pbcshift(imconr,celr,to_read,axx_read,ayy_read,azz_read)
        n_held=0
        Do i=1,to_read
           sxx=rcell(1)*axx_read(i)+rcell(4)*ayy_read(i)+rcell(7)*azz_read(i)
           syy=rcell(2)*axx_read(i)+rcell(5)*ayy_read(i)+rcell(8)*azz_read(i)
           szz=rcell(3)*axx_read(i)+rcell(6)*ayy_read(i)+rcell(9)*azz_read(i)

           sxx=sxx-Anint(sxx) ; If (sxx >= half_minus) sxx=-sxx
           syy=syy-Anint(syy) ; If (syy >= half_minus) syy=-syy
           szz=szz-Anint(szz) ; If (szz >= half_minus) szz=-szz

           axx_read(i)=sxx
           ayy_read(i)=syy
           azz_read(i)=szz

! sxx,syy,szz are in [-0.5,0.5) interval and values as 0.4(9) may pose a problem

           ipx=Int((sxx+0.5_wp)*nprx_r)
           ipy=Int((syy+0.5_wp)*npry_r)
           ipz=Int((szz+0.5_wp)*nprz_r)

! check for errors

           If (ipx == nprx .or. ipy == npry .or. ipz == nprz) Call error(555)

! assign domain

           idm=ipx+nprx*(ipy+npry*ipz)
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
        End Do

     End If Readers_only

! Increase buffer counter and update first_at for
! the readers that have something left to read

     indatm = indatm+1
     If (do_read) Then
        If (first_at(my_read_proc_num) < first_at(my_read_proc_num+1)) &
             first_at(my_read_proc_num) = first_at(my_read_proc_num)+1
     End If

! Circulate configuration data to all nodes when transmission arrays are filled up
! Check against megref since at low processors counts (i.e. 1) batsz can be > megref

     Reorganize_buffer: If (indatm == batsz .or. (indatm > 0 .and. k == megref)) Then

        Do which_read_proc = 0 , n_read_procs_use-1
           If (orig_first_at(which_read_proc) >= megref) Exit ! for non-reading readers

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

Dispatch:  Do i=1,n_loc
              nrefs=nrefs+1

! Check safety by the upper bound of: namr,indr,xr,yr,zr &

              If (nrefs > Size( xr )) Then
                 safe=.false.
                 Exit Dispatch
              End If

              namr(nrefs)=chbuf(i)
              indr(nrefs)=iwrk(i)

              xr(nrefs)=scatter_buffer(1,i)
              yr(nrefs)=scatter_buffer(2,i)
              zr(nrefs)=scatter_buffer(3,i)
           End Do Dispatch
        End Do

! Check if all is dispatched fine

        If (mxnode > 1) Call gcheck(safe)
        If (.not.safe) Call error(556)

! Nullify dispatch counter

        indatm=0

     End If Reorganize_buffer

  End Do

! The last reader to check for EoFile in REFERENCE
! and if none is hit to call error to abort

  If (do_read) Then
     If (first_at(my_read_proc_num) == megref) Then
        recs_to_read = 1
        If (.not.fast) Then
           Read(Unit=nrefdt, Fmt=forma, Iostat=ierr) rec_buff( :, 1:recs_to_read )
        Else
           Call io_read_batch( fh, rec_mpi_io, 1, rec_buff, ierr )
        End If
        safe = (ierr /= 0)
     End If
  End If
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(58)

  If (do_read) Then
     Deallocate(rec_buff,                   Stat=fail(1))
     Deallocate(chbuf_read,iwrk_read,       Stat=fail(2))
     Deallocate(axx_read,ayy_read,azz_read, Stat=fail(3))
     Deallocate(owner_read,                 Stat=fail(4))
     If (Any(fail(1:4) > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_read_parallel deallocation failure 2, node: ', idnode
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
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_read_parallel deallocation failure 1, node: ', idnode
     Call error(0)
  End If

  Return

! error exit for REFERNCE file read

100 Continue
  Call error(554)

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
       Write( nrite, '(/,1x,a,i0)') 'defects_reference_read_parallel allocation failure 4, node: ', idnode
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
       Write( nrite, '(/,1x,a,i0)') 'defects_reference_read_parallel allocation failure 4, node: ', idnode
       Call error( 0 )
    End If

  End Subroutine get_var

End Subroutine defects_reference_read_parallel
