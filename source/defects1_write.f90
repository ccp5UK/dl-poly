Subroutine defects1_write &
           (imcon,rcut,keyres,keyens,nsdef,isdef,rdef,nstep,tstep,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing DEFECTS1 file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2011
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use site_module,     Only : ntpshl,unqshl
  Use config_module,   Only : cfgname,cell,natms,nlast, &
                              atmnam,ltg,lfrzn,xxx,yyy,zzz
  Use defects1_module
  Use parse_module,    Only : tabs_2_blanks, get_word, word_2_real
  Use io_module,       Only : io_set_parameters,        &
                              io_get_parameters,        &
                              io_init, io_open,         &
                              io_write_record,          &
                              io_write_batch,           &
                              io_close, io_finalize,    &
                              IO_WRITE_UNSORTED_MPIIO,  &
                              IO_WRITE_UNSORTED_DIRECT, &
                              IO_WRITE_UNSORTED_MASTER, &
                              IO_WRITE_SORTED_MPIIO,    &
                              IO_WRITE_SORTED_DIRECT,   &
                              IO_WRITE_SORTED_NETCDF,   &
                              IO_WRITE_SORTED_MASTER

  Implicit None

  Integer,           Intent( In    ) :: imcon,keyres,keyens, &
                                        nsdef,isdef,nstep
  Real( Kind = wp ), Intent( In    ) :: rcut,rdef,tstep,time

  Integer, Parameter :: recsz = 73 ! default record size

  Logical,           Save :: newjob = .true.
  Integer(Kind=ip),  Save :: rec    = 0_ip , &
                             frm    = 0_ip
  Integer,           Save :: mxlcdef
  Real( Kind = wp ), Save :: rdefsq,rcell(1:9),cwx,cwy,cwz, &
                             dxl,dxr,dyl,dyr,dzl,dzr,cutdef

  Logical                 :: safe,lexist,l_tmp,ready
  Character( Len = 40 )   :: word
  Integer                 :: fail(1:7),i,j,k,nlx,nly,nlz,   &
                             ic,jc,kk,ix,iy,iz,jx,jy,jz,    &
                             taken,ni,megni,nv,megnv,jdnode,jatms
  Real( Kind = wp )       :: cut,x,y,z,xs,ys,zs,buffer(1:2)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: fh, io_write, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

! Number of neighbouring cells to look around for counting defects

  Integer, Parameter :: nsbcll = 27

! Direction arrays for jumping around in link-cell space

  Integer, Dimension( 1:nsbcll ), Parameter :: &
  nix = (/ 0,  -1,-1,-1, 0, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 1, 1, 1, -1,-1,-1, 0, 0, 0, 1, 1, 1 /) , &
  niy = (/ 0,  -1, 0, 1,-1, 0, 1,-1, 0, 1, -1, 0, 1,-1, 1,-1, 0, 1, -1, 0, 1,-1, 0, 1,-1, 0, 1 /) , &
  niz = (/ 0,  -1,-1,-1,-1,-1,-1,-1,-1,-1,  0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1, 1 /)

  Real( Kind = wp ),    Dimension( : ),    Allocatable :: dr
  Character( Len = 8 ), Dimension( : ),    Allocatable :: namv,nami
  Integer,              Dimension( : ),    Allocatable :: indv,indi,interstitial,occupies
  Integer,              Dimension( : ),    Allocatable :: linkr,lctr,link,lct

  Integer,              Dimension( : ),    Allocatable :: ni_n,nv_n

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: cxx,cyy,czz


  If (.not.(nstep >= nsdef .and. Mod(nstep-nsdef,isdef) == 0)) Return

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! netCDF not implemented for DEFECTS1.  Switch to DEFAULT temporarily.

  If (io_write == IO_WRITE_SORTED_NETCDF) io_write = IO_WRITE_SORTED_MPIIO

  If (newjob) Then
     newjob = .false.

! rdef squared

     rdefsq=rdef**2

! Build lattice sites list from REFERENCE

     Call allocate_defects1_arrays()
     Call defects_reference_read &
           ('REFERENCE1',imcon,nstep,celr1,nrefs1,namr1,indr1,xr1,yr1,zr1)

! Assume that the MD cell will not change much in size and shape from
! the one provided in REFERENCE, a smaller halo(cutoff(rdef)) is to be set

     cut=rdef+0.15_wp
     Call defects_reference_set_halo   &
           (imcon,cut,cwx,cwy,cwz,dxl, &
           dxr,dyl,dyr,dzl,dzr,        &
           nrefs1,nlrefs1,namr1,lri1,lra1,indr1,xr1,yr1,zr1)

! If the keyres=1, is DEFECTS1 old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (idnode == 0) Inquire(File='DEFECTS1', Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist)
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (idnode == 0) Then
           Open(Unit=ndefdt, File='DEFECTS1', Form='formatted', Access='direct', Status='replace', Recl=recsz)
           Write(Unit=ndefdt, Fmt='(a72,a1)',           Rec=1) cfgname(1:72),lf
           Write(Unit=ndefdt, Fmt='(f7.3,a23,2i21,a1)', Rec=2) rdef,Repeat(' ',23),frm,rec,lf
           Close(Unit=ndefdt)
        End If
        rec=Int(2,ip)
        frm=Int(0,ip)

! Get some sense of it

     Else

        safe=.true.
        If (idnode == 0) Then

           Open(Unit=ndefdt, File='DEFECTS1', Form='formatted')

           Do While (.true.)

              record=' '
              If (l_tmp) Then

                 Read(Unit=ndefdt, Fmt=*, End=20)            ! title record
                 rec=rec+Int(1,ip)
                 Read(Unit=ndefdt, Fmt='(a)', End=20) record ! bookkeeping record
                 rec=rec+Int(1,ip)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 If (word(1:Len_Trim(word)) /= 'timestep') Then
                    Call get_word(record,word) ; Call get_word(record,word)
                    Call get_word(record,word) ; frm=Nint(word_2_real(word,0.0_wp),ip)
                    Call get_word(record,word) ; rec=Nint(word_2_real(word,0.0_wp),ip)
                    If (frm /= Int(0,ip) .and. rec > Int(2,ip)) Then
                       Go To 20 ! New style
                    Else
                       l_tmp=.false. ! TOUGH, old style
                       rec=Int(2,ip)
                       frm=Int(0,ip)
                    End If
                 Else
                    safe=.false. ! Overwrite the file, it's junk to me
                    Go To 20
                 End If

              Else

                 Read(Unit=ndefdt, Fmt=*, End=20)            ! timestep record
                 rec=rec+Int(1,ip)

                 Read(Unit=ndefdt, Fmt='(a)', End=20) record ! defects record
                 rec=rec+Int(1,ip)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 Call get_word(record,word) ; j=Nint(word_2_real(word))

                 Do i=1,3+2*j ! 3 lines for cell parameters and 2*j entries for defects
                    Read(Unit=ndefdt, Fmt=*, End=20)
                    rec=rec+Int(1,ip)
                 End Do
                 frm=frm+Int(1,ip)

              End If

           End Do

20         Continue
           Close(Unit=ndefdt)

        End If

        If (mxnode > 1) Call gcheck(safe)
        If (.not.safe) Then
           lexist=.false.

           rec=Int(0,ip)
           frm=Int(0,ip)

           Go To 10
        Else If (mxnode > 1) Then
           buffer(1)=Real(frm,wp)
           buffer(2)=Real(rec,wp)

           Call MPI_BCAST(buffer(1:2), 2, wp_mpi, 0, dlp_comm_world, ierr)

           frm=Nint(buffer(1),ip)
           rec=Nint(buffer(2),ip)
        End If

     End If

! Get rcell

     Call invert(cell,rcell,cut)

! New real space cutoff and expansion for defects link-cells

     cutdef=Min(rcut/3.0_wp,2.0_wp*rdef)
     mxlcdef=Nint(((rcut/cutdef)**3+0.15_wp)*mxcell)
  End If

! Update rcell

  If (keyens >= 20) Call invert(cell,rcell,cut)

  fail=0
  Allocate (dr(1:mxatms),                                                Stat=fail(1))
  Allocate (namv(1:mxatms),indv(1:mxatms),nami(1:mxatms),indi(1:mxatms), Stat=fail(2))
  Allocate (interstitial(1:mxatms),occupies(1:mxatms),                   Stat=fail(3))
  Allocate (linkr(1:mxatms),link(1:mxatms),                              Stat=fail(4))
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms),                   Stat=fail(5))
  Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms),                   Stat=fail(6))
  Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms),                   Stat=fail(7))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects1_write allocation failure, node: ', idnode
     Call error(0)
  End If

! Build bookkeeping lists: interstitial, occupies

  Do i=1,nlast

! Consider all particles at this point:j=1 - consider in,j=0 - consider out

     j=1

! Get all domain+halo particles' coordinates in MD cell centred reduced space

     cxx(i)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     cyy(i)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     czz(i)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! Exclude particles in the domain's halo farther than cutoff(rdef)
! smaller halo with a width defined by cutoff(rdef)

     If ( i > natms .and.                        &
          ((cxx(i) < dxl .or. cxx(i) > dxr) .or. &
           (cyy(i) < dyl .or. cyy(i) > dyr) .or. &
           (czz(i) < dzl .or. czz(i) > dzr)) ) j=0

! Exclude frozen and shell particles from consideration

     If ( j == 1 .and. (lfrzn(i) /= 0 .or. Any(unqshl(1:ntpshl) == atmnam(i))) ) j=0

! Assume that every considered particles (1) is an interstitial and
! and (2) does not occupy a site yet

     If (j == 1) Then
        interstitial(i)=-1 ! considered particle is assumed to be an interstitial
        occupies(i)    = 0 ! and to not occupy a site yet (but may occupy one later)
     Else
        interstitial(i)= 0 ! excluded particle can neither ever be an interstitial
        occupies(i)    =-1 ! nor can it ever occupy a site
     End If
  End Do

! Partition sites and atoms in link-celled space with same imcon, cell and cutoff!!!

  Allocate (lctr(1:mxlcdef),lct(1:mxlcdef), Stat=fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects1_write allocation failure 1, node: ', idnode
     Call error(0)
  End If
  Call defects_link_cells &
           (imcon,cell,cutdef,mxlcdef,nrefs1,nlrefs1,xr1,yr1,zr1,nlx,nly,nlz,linkr,lctr)
  Call defects_link_cells &
           (imcon,cell,cutdef,mxlcdef,natms,nlast,cxx,cyy,czz,nlx,nly,nlz,link,lct)

  safe = .true.          ! Initialise safety flag to all safe
  nv = 0                 ! Assume no vacancies (all sites are occupied)
  ni = 0                 ! Assume no interstitials (all particles are occupying a site)
  dr = (rdef+0.15_wp)**2 ! Assume all particles are > rdef distance away from the nearest to them site

! Actual defect detection *********************************************

! primary loop over domain subcells of sites

  Do iz=0,nlz+1
     Do iy=0,nly+1
        Do ix=0,nlx+1

! index of primary reference link-cell

           ic=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! bypass reference subcell if empty

           If (lctr(ic) > 0) Then

! head of chain of i-th subcell (reference)

              i=lctr(ic)

! loop over primary cell contents

100           Continue

! Bypass if the site is a shell

              If (Any(unqshl(1:ntpshl) == namr1(i))) Go To 400

! Assume the site is vacant

              taken=0 ! Assume the site is vacant

! secondary loop over subcells of paticles

              Do kk=1,nsbcll

                 jx=ix+nix(kk)
                 jy=iy+niy(kk)
                 jz=iz+niz(kk)

! SHIFT BACK to the LEFT
! Disregard cells outside the look-up scope

  If ( (jx >= 0) .and. (jx <= nlx+1) .and. &
       (jy >= 0) .and. (jy <= nly+1) .and. &
       (jz >= 0) .and. (jz <= nlz+1) ) Then

! index of neighbouring cell

     jc=1+jx+(nlx+2)*(jy+(nly+2)*jz)

! bypass real subcell if empty

     If (lct(jc) > 0) Then

! head of chain of j-th subcell (real)

        j=lct(jc)

! loop over secondary cell contents

200     Continue

! Bypass if the atom is excluded or has already claimed (occupied) a site

        If (occupies(j) /= 0) Go To 300

! Get vectors particle-site (in reduced space)

        xs=xr1(i)-cxx(j) ; If (Abs(xs) > cwx) Go To 300
        ys=yr1(i)-cyy(j) ; If (Abs(ys) > cwy) Go To 300
        zs=zr1(i)-czz(j) ; If (Abs(zs) > cwz) Go To 300

! Get in real space

        x=cell(1)*xs+cell(4)*ys+cell(7)*zs
        y=cell(2)*xs+cell(5)*ys+cell(8)*zs
        z=cell(3)*xs+cell(6)*ys+cell(9)*zs

! Get particle-site squared distance

        cut=x**2+y**2+z**2

! An atom qualifies to claim a site (occupies > 0)
! when the site is vacant - not occupied/taken (taken=0).
! If the site is not claimed yet then the atom is marked
! as a claimee of this site (occupies > 0) and not an
! interstitial (interstitial = 0).  If the site is
! already taken by another atom then the atom is marked
! as an interstitial of this site (interstitial > 0).

        If (cut < rdefsq) Then
           dr(j)=cut             ! Save the distance (site-atom)

! Only general interstitial atoms are allowed to claim.  Otherwise, we are in a big trouble

           If (interstitial(j) >= 0) safe=.false.

           If (taken == 0) Then
              taken=1            ! This site is taken (no longer vacant)
              occupies(j)=i      ! This atom claims a site (this site)
              interstitial(j)=0  ! This atom is not an interstitial any longer
           Else
              interstitial(j)=i  ! This atom is not a general interstitial any longer
           End If
        End If

300     Continue

        j=link(j)
        If (j /= 0) Go To 200

! end of loop over the real subcell contents jc

     End If

! End If of disregard cell outside the look-up scope

  End If

! SHIFT BACK to the RIGHT
! End Do secondary loop around real subcells

              End Do

! Vacancies: taken=0 means we've found a vacancy

              If (taken == 0 .and. i <= nrefs1) Then
                 nv=nv+1
                 namv(nv)=namr1(i)
                 indv(nv)=indr1(i)
                 axx(nv)=cell(1)*xr1(i)+cell(4)*yr1(i)+cell(7)*zr1(i)
                 ayy(nv)=cell(2)*xr1(i)+cell(5)*yr1(i)+cell(8)*zr1(i)
                 azz(nv)=cell(3)*xr1(i)+cell(6)*yr1(i)+cell(9)*zr1(i)
              End If

400           Continue
              i=linkr(i)
              If (i /= 0) Go To 100

! end of bypass of empty subcell ic

           End If
        End Do
     End Do
  End Do

! Check safety on rdef length

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(560)

! Define real interstitials and real claimees of a site by swapping
! interstitial/occupation statuses between particles based on
! site-atom distances of these.  Amongst all atoms close to sites -
! interstitials of the site and the site claimee -
! the real claimee is the atom closest to the site.

! primary loop over domain subcells of particles

  Do iz=0,nlz+1
     Do iy=0,nly+1
        Do ix=0,nlx+1

! index of primary link-cell

           ic=1+ix+(nlx+2)*(iy+(nly+2)*iz)

! bypass subcell if empty

           If (lct(ic) > 0) Then

! head of chain of i-th subcell

              i=lct(ic)

! loop over primary cell contents

500           Continue

! Bypass if the atom is excluded or is not a claimee

              If (occupies(i) <= 0) Go To 700

! Circle around domain+halo particles, k is the index of the
! closest to a site atom (real claimee) assuming that the first
! claimee is the true one.  Check for an interstitial of the
! same site, that is closer to it than the claimee.

              k=i
              cut=dr(i)

! secondary loop over subcells of particles

              Do kk=1,nsbcll

                 jx=ix+nix(kk)
                 jy=iy+niy(kk)
                 jz=iz+niz(kk)

! SHIFT BACK to the LEFT
! Disregard cells outside the look-up scope

  If ( (jx >= 0) .and. (jx <= nlx+1) .and. &
       (jy >= 0) .and. (jy <= nly+1) .and. &
       (jz >= 0) .and. (jz <= nlz+1) ) Then

! index of neighbouring cell

     jc=1+jx+(nlx+2)*(jy+(nly+2)*jz)

! bypass real subcell if empty

     If (lct(jc) > 0) Then

! head of chain of j-th subcell (real)

        j=lct(jc)

! loop over secondary cell contents

600     Continue

        If (interstitial(j) == occupies(i) .and. dr(j) < cut) Then
           cut=dr(j)
           k=j
        End If

        j=link(j)
        If (j /= 0) Go To 600

! end of loop over the subcell contents jc

     End If

! End If of disregard cell outside the look-up scope

  End If

! SHIFT BACK to the RIGHT
! End Do secondary loop around subcells

              End Do

! If there is a new true claimee change statuses

              If (k /= i) Then
                 interstitial(i)=interstitial(k) ; interstitial(k)=0
                 occupies(k)=occupies(i)         ; occupies(i)=0
              End If

700           Continue

              i=link(i)
              If (i /= 0) Go To 500

! end of bypass of empty subcell ic

           End If
        End Do
     End Do
  End Do

! Include the following if you want occupant/site type mismatches printed out
! printing from many nodes can be messy
!
!     taken = 0
!     Do i=1,natms
!        If (occupies(i) > 0 .and. atmnam(i) /= namr1(occupies(i))) Then
!!           Write(Unit=*,'(3(1x,a))') 'occupant/site type mismatch:', atmnam(i), namr1(occupies(i))
!           taken = taken + 1
!        End If
!     End Do
!     If (mxnode > 1) Call gsum(taken)
!     If (idnode == 0) Write(nrite,'(3(1x,a,i10))') 'occupant/site type mismatches', taken

! Interstitials: i <= natms & interstitial(i) /= 0 means we've found an interstitial

  Do i=1,natms
     If (interstitial(i) /= 0) Then
        ni=ni+1
        nami(ni)=atmnam(i)
        indi(ni)=ltg(i)
        bxx(ni)=xxx(i)
        byy(ni)=yyy(i)
        bzz(ni)=zzz(i)
     End If
  End Do

! Sum global defects values

  megni = ni
  If (mxnode > 1) Call gsum(megni)

  megnv = nv
  If (mxnode > 1) Call gsum(megnv)

! Get relative offsets (int,vac) for parallel printing

  Allocate (ni_n(0:mxnode),nv_n(0:mxnode), Stat=fail(1))
  Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects1_write allocation failure 2, node: ', idnode
     Call error(0)
  End If

  ni_n=0 ; ni_n(idnode+1)=ni
  If (mxnode > 1) Call gsum(ni_n)
  ni_n(0)=2*Sum(ni_n(0:idnode)) ! 2 lines per record

  nv_n=0 ; nv_n(idnode+1)=nv
  If (mxnode > 1) Call gsum(nv_n)
  nv_n(0)=2*Sum(nv_n(0:idnode)) ! 2 lines per record

  chbat=' '

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(kind = MPI_OFFSET_KIND)

! Update frame

  frm=frm+Int(1,ip)

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     j=0
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_open( io_write, MPI_COMM_SELF, 'DEFECTS1', MPI_MODE_WRONLY, fh )

        Write(record, Fmt='(a8,i10,2f12.6,i5,f7.3,a18,a1)') &
           'timestep',nstep,tstep,time,imcon,rdef,Repeat(' ',18),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(a8,i10,a15,i10,a11,i10,a8,a1)') &
           'defects ',megni+megnv, ' interstitials ',megni, ' vacancies ',megnv,Repeat(' ',8),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Do i = 0, 2
           Write( record, '( 3f20.10, a12, a1 )' ) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           j=j+1
           Do k=1,recsz
              chbat(k,j) = record(k:k)
           End Do
        End Do

! Dump header and cell information

        Call io_write_batch( fh, rec_mpi_io, j, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        j=j+5

     End If
     Call gsync()

! Start of file

     rec=rec+Int(j,ip)

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, 'DEFECTS1', MPI_MODE_WRONLY, fh )

! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(ni_n(0),MPI_OFFSET_KIND)
     j=0
     Do i=1,ni
        Write(record, Fmt='(a2,a8,i10,a52,a1)') 'i_',nami(i),indi(i),Repeat(' ',52),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

! Dump batch and update start of file

        If (j + 2 >= batsz .or. i == ni) Then
           Call io_write_batch( fh, rec_mpi_io, j, chbat )
           rec_mpi_io=rec_mpi_io+Int(j,MPI_OFFSET_KIND)
           j=0
        End If
     End Do

! Start of file

     rec=rec+Int(2*megni,ip)

! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(nv_n(0),MPI_OFFSET_KIND)
     Do i=1,nv
        Write(record, Fmt='(a2,a8,i10,a52,a1)') 'v_',namv(i),indv(i),Repeat(' ',52),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

! Dump batch and update start of file

        If (j + 2 >= batsz .or. i == nv) Then
           Call io_write_batch( fh, rec_mpi_io, j, chbat )
           rec_mpi_io=rec_mpi_io+Int(j,MPI_OFFSET_KIND)
           j=0
        End If
     End Do
     rec=rec+Int(2*megnv,ip)

! Update main header

     If (idnode == 0) Then
        Write(record, Fmt='(f7.3,a23,2i21,a1)') rdef,Repeat(' ',23),frm,rec,lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record )
     End If

     Call io_close( fh )
     Call io_finalize

  Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
           io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'defects1_write allocation failure 3, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O
! Start of file

     j=0
     If (idnode == 0) Then
        Open(Unit=ndefdt, File='DEFECTS1', Form='formatted', Access='direct', Recl=recsz)

! Accumulate header

        Write(record, Fmt='(a8,i10,2f12.6,i5,f7.3,a18,a1)') &
           'timestep',nstep,tstep,time,imcon,rdef,Repeat(' ',18),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(a8,i10,a15,i10,a11,i10,a8,a1)') &
           'defects ',megni+megnv, ' interstitials ',megni, ' vacancies ',megnv,Repeat(' ',8),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record, Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           j=j+1
           Do k=1,recsz
              chbat(k,j) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=ndefdt, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,j)
        rec=rec+Int(5,ip)
        j=0

        Do i=1,ni
           chbuf(i)=nami(i)
           iwrk(i)=indi(i)

           cxx(i)=bxx(i)
           cyy(i)=byy(i)
           czz(i)=bzz(i)
        End Do

        jatms=ni
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,DefWrite_tag,dlp_comm_world,ierr)
              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)

              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(czz,jatms,wp_mpi,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
              End If
           End If

           Do i=1,jatms
              Write(record, Fmt='(a2,a8,i10,a52,a1)') 'i_',chbuf(i),iwrk(i),Repeat(' ',52),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

! Dump batch and update start of file

              If (j + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=ndefdt, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,j)
                 rec=rec+Int(j,ip)
                 j=0
              End If
           End Do
        End Do
     Else
        Call MPI_RECV(ready,1,MPI_LOGICAL,0,DefWrite_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(ni,1,MPI_INTEGER,0,DefWrite_tag,dlp_comm_world,ierr)

        If (ni > 0) Then
           Call MPI_SEND(nami,8*ni,MPI_CHARACTER,0,DefWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(indi,ni,MPI_INTEGER,0,DefWrite_tag,dlp_comm_world,ierr)

           Call MPI_SEND(bxx,ni,wp_mpi,0,DefWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(byy,ni,wp_mpi,0,DefWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(bzz,ni,wp_mpi,0,DefWrite_tag,dlp_comm_world,ierr)
        End If
     End If

     If (idnode == 0) Then
        Do i=1,nv
           chbuf(i)=namv(i)
           iwrk(i)=indv(i)

           cxx(i)=axx(i)
           cyy(i)=ayy(i)
           czz(i)=azz(i)
        End Do

        jatms=nv
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,DefWrite_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)

              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(czz,jatms,wp_mpi,jdnode,DefWrite_tag,dlp_comm_world,status,ierr)
              End If
           End If

           Do i=1,jatms
              Write(record, Fmt='(a2,a8,i10,a52,a1)') 'v_',chbuf(i),iwrk(i),Repeat(' ',52),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

! Dump batch and update start of file

              If (j + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=ndefdt, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,j)
                 rec=rec+Int(j,ip)
                 j=0
              End If
           End Do
        End Do
     Else
        Call MPI_RECV(ready,1,MPI_LOGICAL,0,DefWrite_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(nv,1,MPI_INTEGER,0,DefWrite_tag,dlp_comm_world,ierr)

        If (nv > 0) Then
           Call MPI_SEND(namv,8*nv,MPI_CHARACTER,0,DefWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(indv,nv,MPI_INTEGER,0,DefWrite_tag,dlp_comm_world,ierr)

           Call MPI_SEND(axx,nv,wp_mpi,0,DefWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ayy,nv,wp_mpi,0,DefWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(azz,nv,wp_mpi,0,DefWrite_tag,dlp_comm_world,ierr)
        End If
     End If

! Update main header

     If (idnode == 0) Then
        Write(Unit=ndefdt, Fmt='(f7.3,a23,2i21,a1)', Rec=2) rdef,Repeat(' ',23),frm,rec,lf
        Close(Unit=ndefdt)
     Else
        rec=rec+Int(5+2*(megni+megnv),ip)
     End If

     Deallocate (chbuf,iwrk, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'defects1_write deallocation failure 3, node: ', idnode
        Call error(0)
     End If

  End If

  If (mxnode > 1) Call gsync()

  Deallocate (ni_n,nv_n, Stat=fail(1))
  Deallocate (chbat,     Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects1_write deallocation failure 2, node: ', idnode
     Call error(0)
  End If

  Deallocate (dr,                    Stat=fail(1))
  Deallocate (namv,indv,nami,indi,   Stat=fail(2))
  Deallocate (interstitial,occupies, Stat=fail(3))
  Deallocate (linkr,lctr,link,lct,   Stat=fail(4))
  Deallocate (axx,ayy,azz,           Stat=fail(5))
  Deallocate (bxx,byy,bzz,           Stat=fail(6))
  Deallocate (cxx,cyy,czz,           Stat=fail(7))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'defects1_write deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine defects1_write
