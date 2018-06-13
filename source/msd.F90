Module msd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring msd routines variables
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use kinds, Only : wp,li,wi
  Use comms, Only : comms_type, gcheck,MsdWrite_tag,gsum,wp_mpi,gsync,gbcast, &
                    gsend,grecv,offset_kind,comm_self,mode_wronly
  Use setup
  Use site, Only : site_type
  Use configuration,     Only : cfgname,natms,atmnam,lsite,ltg, &
                                weight

  Use parse,      Only : tabs_2_blanks, get_word, word_2_real
  Use io,         Only : io_set_parameters,             &
                                io_get_parameters,             &
                                io_init,                       &
                                io_open, io_write_record,      &
                                io_write_batch,                &
                                io_write_sorted_file,          &
                                io_close, io_finalize,         &
                                IO_MSDTMP,                     &
                                IO_BASE_COMM_NOT_SET,          &
                                IO_ALLOCATION_ERROR,           &
                                IO_UNKNOWN_WRITE_OPTION,       &
                                IO_UNKNOWN_WRITE_LEVEL,        &
                                IO_WRITE_UNSORTED_MPIIO,       &
                                IO_WRITE_UNSORTED_DIRECT,      &
                                IO_WRITE_UNSORTED_MASTER,      &
                                IO_WRITE_SORTED_MPIIO,         &
                                IO_WRITE_SORTED_DIRECT,        &
                                IO_WRITE_SORTED_NETCDF,        &
                                IO_WRITE_SORTED_MASTER
  Use errors_warnings, Only : error
  Implicit None

  Private

  Type, Public :: msd_type
    Private

    !> MSD recording switch
    Logical, Public :: l_msd = .false.

    !> Step to begin recording MSD
    Integer( Kind = wi ), Public :: start
    !> Frequency to record MSD (steps)
    Integer( Kind = wi ), Public :: freq
  End Type msd_type


  Public :: msd_write
  Contains

    Subroutine msd_write(keyres,megatm,nstep,tstep,time,stpval,dof_site,msd_data,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing MSDTMP file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: keyres, &
                                        megatm,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep,time
  Real( Kind = wp ), Intent( InOut ) :: stpval(:)
  Real( Kind = wp ), Dimension(:), Intent( In    ) :: dof_site
  Type( msd_type ), Intent( Inout ) :: msd_data
  Type( comms_type), Intent( InOut ) :: comm

  Integer, Parameter :: recsz = 53 ! default record size

  Logical,               Save :: newjob = .true. , &
                                fast   = .true.
  Character( Len = 40 ), Save :: fname
  Integer(Kind=li),      Save :: rec = 0_li , &
                                frm = 0_li

  Logical                :: lexist,safe,ready
  Character( Len = 40 )  :: word
  Integer                :: fail(1:2),i,jj,k,jdnode,jatms
  Real( Kind = wp )      :: buffer(1:2),tmp

! Some parameters and variables needed by io interfaces

  Integer                       :: fh, io_write, batsz
  Integer( Kind = offset_kind ) :: rec_mpi_io
  Character( Len = recsz )      :: record
  Character                     :: lf

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: ddd,eee
  Integer :: ierr
  Character( Len = 256 ) :: message

  If (.not.(nstep >= msd_data%start .and. Mod(nstep-msd_data%start,msd_data%freq) == 0)) Return

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! netCDF not implemented for MSDTMP.  Switch to DEFAULT temporarily.

  If (io_write == IO_WRITE_SORTED_NETCDF) io_write = IO_WRITE_SORTED_MPIIO

  If (newjob) Then
    newjob = .false.

! name convention

    fname = 'MSDTMP'

! If keyres=1, is MSDTMP old (does it exist) and
! how many frames and records are in there

    lexist=.true.
    If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File=fname, Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
    Else
        lexist=.false.
    End If

! Generate file is non-existent

10   Continue
    If (.not.lexist) Then

        If (comm%idnode == 0) Then
          Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Status='replace', Recl=recsz)
          Write(Unit=nhist, Fmt='(a52,a1)',      Rec=Int(1,li)) cfgname(1:52),lf
          Write(Unit=nhist, Fmt='(i10,2i21,a1)', Rec=Int(2,li)) megatm,frm,rec,lf
          Close(Unit=nhist)
        End If
        rec=Int(2,li)
        frm=Int(0,li)

! Get some sense of it

    Else

        safe=.true.
        If (comm%idnode == 0) Then

          Open(Unit=nhist, File=fname, Form='formatted')

          Do

              record=' '

! Assume new style of MSDTMP with bookkeeping.

              If (fast) Then

                Read(Unit=nhist, Fmt=*, End=10)            ! title record
                rec=rec+Int(1,li)
                record=' '
                Read(Unit=nhist, Fmt='(a)', End=10) record ! bookkeeping record
                Call tabs_2_blanks(record)
                rec=rec+Int(1,li)

                Call get_word(record,word)
                If (word(1:Len_Trim(word)) /= 'timestep') Then
                    Call get_word(record,word) ; Call get_word(record,word)
                    Call get_word(record,word) ; frm=Nint(word_2_real(word,0.0_wp),li)
                    Call get_word(record,word) ; rec=Nint(word_2_real(word,0.0_wp),li)
                    If (frm /= Int(0,li) .and. rec > Int(2,li)) Then
                      Go To 20 ! New style
                    Else
                      fast=.false. ! TOUGH, old style
                      rec=Int(2,li)
                      frm=Int(0,li)
                    End If
                Else
                    safe=.false. ! Overwrite the file, it's junk to me
                    Go To 20
                End If

! TOUGH, it needs scanning through

              Else

                Read(Unit=nhist, Fmt='(a)', End=20) record              ! timestep record
                Call tabs_2_blanks(record)
                rec=rec+Int(1,li)

                Call get_word(record,word)
                Call get_word(record,word) ; jj=Nint(word_2_real(word)) ! total number of lines to read

                word=' '
                Write(word,'( "(", i0, "( / ) )" )') jj-1
                Read(Unit=nhist, Fmt=word, End=20)
                rec=rec+Int(jj,li)
                frm=frm+Int(1,li)

              End If

          End Do

20         Continue
          Close(Unit=nhist)

        End If

        Call gcheck(comm,safe,"enforce")
        If (.not.safe) Then
          lexist=.false.

          rec=Int(0,li)
          frm=Int(0,li)

          Go To 10
        Else If (comm%mxnode > 1) Then
          buffer(1)=Real(frm,wp)
          buffer(2)=Real(rec,wp)

          Call gbcast(comm,buffer(1:2),9)

          frm=Nint(buffer(1),li)
          rec=Nint(buffer(2),li)
        End If

    End If
  End If

! Get offsets and define batch

  fail=0
  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
    Allocate (n_atm(0:comm%mxnode),        Stat=fail(1))
    Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
    If (Any(fail > 0)) Then
        Write(message,'(a,i0)') 'msd_write allocation failure 0'
        Call error(0,message)
    End If

    chbat=' '
    n_atm=0 ; n_atm(comm%idnode+1)=natms
    Call gsum(comm,n_atm)
    n_atm(0)=Sum(n_atm(0:comm%idnode))
  End If

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = offset_kind)

! Update frame

  frm=frm+Int(1,li)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
          io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header information, where just one node is needed
! Start of file

    rec_mpi_io=Int(rec,offset_kind)
    If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_open( io_write, comm_self, fname, mode_wronly, fh )

        Write(record, Fmt='(a8,2i10,2f12.6,a1)') 'timestep',nstep,megatm,tstep,time,lf

! Dump header information

        Call io_write_record( fh, rec_mpi_io, record )

        Call io_close( fh )
        Call io_finalize

    End If

! Start of file

    rec=rec+Int(1,li)
    rec_mpi_io=Int(rec,offset_kind)+Int(n_atm(0),offset_kind)
    jj=0

    Call io_set_parameters( user_comm = comm%comm )
    Call io_init( recsz )
    Call io_open( io_write, comm%comm, fname, mode_wronly, fh )

    Do i=1,natms
        k=2*i

        If (Abs(dof_site(lsite(i))) > zero_plus) tmp=weight(i)*stpval(27+k)/(boltz*3.0_wp)

        Write(record, Fmt='(a8,i10,1p,2e13.4,a8,a1)') atmnam(i),ltg(i),Sqrt(stpval(27+k-1)),tmp,Repeat(' ',8),lf
        jj=jj+1
        Do k=1,recsz
          chbat(k,jj) = record(k:k)
        End Do

! Dump batch and update start of file

        If (jj >= batsz .or. i == natms) Then
          Call io_write_batch( fh, rec_mpi_io, jj, chbat )
          rec_mpi_io=rec_mpi_io+Int(jj,offset_kind)
          jj=0
        End If
    End Do

! Update and save offset pointer

    rec=rec+Int(megatm,li)
    If (comm%idnode == 0) Then
        Write(record, Fmt='(i10,2i21,a1)') megatm,frm,rec,lf
        Call io_write_record( fh, Int(1,offset_kind), record )
    End If

    Call io_close( fh )
    Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

    Allocate (chbuf(1:mxatms),iwrk(1:mxatms),ddd(1:mxatms),eee(1:mxatms), Stat=fail(1))
    If (fail(1) > 0) Then
        Write(message,'(a)') 'msd_write allocation failure'
        Call error(0,message)
    End If

! node 0 handles I/O

    If (comm%idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,li)
        Write(Unit=nhist, Fmt='(a8,2i10,2f12.6,a1)', Rec=rec) 'timestep',nstep,megatm,tstep,time,lf

        Do i=1,natms
          k=2*i

          iwrk(i)=ltg(i)
          chbuf(i)=atmnam(i)

          ddd(i)=Sqrt(stpval(27+k-1))

          If (Abs(dof_site(lsite(i))) > zero_plus) tmp=weight(i)*stpval(27+k)/(boltz*3.0_wp)
          eee(i)=tmp
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
          If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,MsdWrite_tag)

              Call grecv(comm,jatms,jdnode,MsdWrite_tag)
              If (jatms > 0) Then
                Call grecv(comm,chbuf(1:jatms),jdnode,MsdWrite_tag)
                Call grecv(comm,iwrk(1:jatms),jdnode,MsdWrite_tag)

                Call grecv(comm,ddd(1:jatms),jdnode,MsdWrite_tag)
                Call grecv(comm,eee(1:jatms),jdnode,MsdWrite_tag)
              End If
          End If

          jj=0
          Do i=1,jatms
              Write(record, Fmt='(a8,i10,1p,2e13.4,a8,a1)') chbuf(i),iwrk(i),ddd(i),eee(i),Repeat(' ',8),lf
              jj=jj+1
              Do k=1,recsz
                chbat(k,jj) = record(k:k)
              End Do

! Dump batch and update start of file

              If (jj >= batsz .or. i == jatms) Then
                Write(Unit=nhist, Fmt='(53a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
                rec=rec+Int(jj,li)
                jj=0
              End If
          End Do
        End Do

! Update main header

        Write(Unit=nhist, Fmt='(i10,2i21,a1)', Rec=Int(2,li)) megatm,frm,rec,lf

        Close(Unit=nhist)

    Else

        Do i=1,natms
          k=2*i

          ddd(i)=Sqrt(stpval(27+k-1))

          If (Abs(dof_site(lsite(i))) > zero_plus) tmp=weight(i)*stpval(27+k)/(boltz*3.0_wp)
          eee(i)=tmp
        End Do

        Call grecv(comm,ready,0,MsdWrite_tag)

        Call gsend(comm,natms,0,MsdWrite_tag)
        If (natms > 0) Then
          Call gsend(comm,atmnam(1:natms),0,MsdWrite_tag)
          Call gsend(comm,ltg(1:natms),0,MsdWrite_tag)

          Call gsend(comm,ddd(1:natms),0,MsdWrite_tag)
          Call gsend(comm,eee(1:natms),0,MsdWrite_tag)
        End If

! Save offset pointer

        rec=rec+Int(megatm+1,li)

    End If

    Deallocate (chbuf,iwrk,ddd,eee, Stat=fail(1))
    If (fail(1) > 0) Then
        Write(message,'(a)') 'msd_write deallocation failure'
        Call error(0,message)
    End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
          io_write == IO_WRITE_SORTED_DIRECT) Then

    rec_mpi_io=Int(rec,offset_kind)
    jj=0
    If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = comm_self )
        Call io_init( recsz )
        Call io_open( io_write, comm_self, fname, mode_wronly, fh )

! Write header information

        Write(record, Fmt='(a8,2i10,2f12.6,a1)') 'timestep',nstep,megatm,tstep,time,lf
        Call io_write_record( fh, rec_mpi_io, record )

        Call io_close( fh )
        Call io_finalize

    End If
    Call gsync(comm)

! Start of file

    rec_mpi_io=Int(rec,offset_kind)+Int(1,li)

! Write the rest

    Call io_set_parameters( user_comm = comm%comm )
    Call io_init( recsz )
    Call io_open( io_write, comm%comm, fname,  mode_wronly, fh )

    Allocate (ddd(1:mxatms),eee(1:mxatms), Stat=fail(1))
    If (fail(1) > 0) Then
        Write(message,'(a)') 'msd_write allocation failure'
        Call error(0,message)
    End If

    Do i=1,natms
        k=2*i

        ddd(i)=Sqrt(stpval(27+k-1))

        If (Abs(dof_site(lsite(i))) > zero_plus) tmp=weight(i)*stpval(27+k)/(boltz*3.0_wp)
        eee(i)=tmp
    End Do

    Call io_write_sorted_file( fh, 0, IO_MSDTMP, rec_mpi_io, natms, &
          ltg, atmnam, ddd, eee, (/ 0.0_wp /),                       &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),  ierr )

    If ( ierr /= 0 ) Then
        Select Case( ierr )
        Case( IO_BASE_COMM_NOT_SET )
          Call error( 1050 )
        Case( IO_ALLOCATION_ERROR )
          Call error( 1053 )
        Case( IO_UNKNOWN_WRITE_OPTION )
          Call error( 1056 )
        Case( IO_UNKNOWN_WRITE_LEVEL )
          Call error( 1059 )
        End Select
    End If

! Update and save offset pointer

    rec=rec+Int(1,li)+Int(megatm,li)
    If (comm%idnode == 0) Then
        Write(record, Fmt='(i10,2i21,a1)') megatm,frm,rec,lf
        Call io_write_record( fh, Int(1,offset_kind), record )
    End If

    Call io_close( fh )
    Call io_finalize

    Deallocate (ddd,eee, Stat=fail(1))
    If (fail(1) > 0) Then
        Write(message,'(a)') 'msd_write deallocation failure'
        Call error(0,message)
    End If

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

    Allocate (chbuf(1:mxatms),iwrk(1:mxatms),ddd(1:mxatms),eee(1:mxatms), Stat=fail(1))
    If (fail(1) > 0) Then
        Write(message,'(a)') 'msd_write allocation failure'
        Call error(0,message)
    End If

! node 0 handles I/O

    If (comm%idnode == 0) Then

        Open(Unit=nhist, File='fname', Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,li)
        Write(Unit=nhist, Fmt='(a8,2i10,2f12.6,a1)', Rec=rec) 'timestep',nstep,megatm,tstep,time,lf

        Do i=1,natms
          k=2*i

          iwrk(i)=ltg(i)
          chbuf(i)=atmnam(i)

          ddd(i)=Sqrt(stpval(27+k-1))

          If (Abs(dof_site(lsite(i))) > zero_plus) tmp=weight(i)*stpval(27+k)/(boltz*3.0_wp)
          eee(i)=tmp
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,comm%mxnode-1
          If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,MsdWrite_tag)

              Call grecv(comm,jatms,jdnode,MsdWrite_tag)
              If (jatms > 0) Then
                Call grecv(comm,chbuf(1:jatms),jdnode,MsdWrite_tag)
                Call grecv(comm,iwrk(1:jatms),jdnode,MsdWrite_tag)

                Call grecv(comm,ddd(1:jatms),jdnode,MsdWrite_tag)
                Call grecv(comm,eee(1:jatms),jdnode,MsdWrite_tag)
              End If
          End If

          Do i=1,jatms
              Write(Unit=nhist, Fmt='(a8,i10,1p,2e13.4,a8,a1)', Rec=rec+Int(iwrk(i),li)) &
                  chbuf(i),iwrk(i),ddd(i),eee(i),Repeat(' ',8),lf
          End Do
        End Do

! Update and save offset pointer

        rec=rec+Int(megatm,li)
        Write(Unit=nhist, Fmt='(i10,2i21,a1)', Rec=Int(2,li)) megatm,frm,rec,lf

        Close(Unit=nhist)

    Else

        Do i=1,natms
          k=2*i

          ddd(i)=Sqrt(stpval(27+k-1))

          If (Abs(dof_site(lsite(i))) > zero_plus) tmp=weight(i)*stpval(27+k)/(boltz*3.0_wp)
          eee(i)=tmp
        End Do

        Call grecv(comm,ready,0,MsdWrite_tag)

        Call gsend(comm,natms,0,MsdWrite_tag)
        If (natms > 0) Then
          Call gsend(comm,atmnam(1:natms),0,MsdWrite_tag)
          Call gsend(comm,ltg(1:natms),0,MsdWrite_tag)

          Call gsend(comm,ddd(1:natms),0,MsdWrite_tag)
          Call gsend(comm,eee(1:natms),0,MsdWrite_tag)
        End If

        rec=rec+Int(megatm+1,li)

    End If

    Deallocate (chbuf,iwrk,ddd,eee, Stat=fail(1))
    If (fail(1) > 0) Then
        Write(message,'(a)') 'msd_write deallocation failure'
        Call error(0,message)
    End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
    Deallocate (n_atm, Stat=fail(1))
    Deallocate (chbat, Stat=fail(2))
    If (Any(fail > 0)) Then
        Write(message,'(a)') 'msd_write deallocation failure 0'
        Call error(0,message)
    End If
  End If

  Call gsync(comm)

End Subroutine msd_write
End Module msd
