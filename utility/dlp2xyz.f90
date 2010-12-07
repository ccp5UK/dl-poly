Program dlp2xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly utility to convert DL_POLY CONFIG format to XYZ format
! (0) standard dl_poly CONFIG file
!     title                                                   LINE
!     levcfg & imcon                                          LINE
!     lattice matrix - vector a components (Angstroms)        LINE
!     lattice matrix - vector b components (Angstroms)        LINE
!     lattice matrix - vector c components (Angstroms)        LINE
!     particle name (and global index)                        LINE
!     particle absolute coordinates (centred MD box origin)   LINE
!     optional particle velocity components                   LINE
!     optional particle force components                      LINE
! (1) xyz cell file
!     total number of particles                               LINE
!     title                                                   LINE
!     lattice matrix - vector a components (Angstroms)        LINE
!     lattice matrix - vector b components (Angstroms)        LINE
!     lattice matrix - vector c components (Angstroms)        LINE
!     particle name and absolute coordinates                  LINE
!
! copyright - daresbury laboratory
! authors   - i.t.todorov january 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use parse_module

  Implicit None

  Character         :: line*62
  Data                 line / '                                                              ' /
  Character         :: cfgname*72,record*100,word*40,fin*40,fout*40,name*8,lf*1
  Integer           :: nx,ny,nz,nread,nrite,levcfg,imcon,n,natms,i,j,k,l,rec
  Real( Kind = wp ) :: fx,fy,fz,t,cell(1:9),conv,cosa,cosb,cosg,sing,z2,z3,        &
                       f1,f2,f3,f4,f5,f6,f7,f8,f9,xx0,yy0,zz0,xx,yy,zz,celprp(1:10)

! define newline character Char(10) or Char(Iachar('\n'))

  lf = Char(10)

! Print header

  Write(*,'(18(1x,a,/))')                                                                &
     "**************  dl_poly utility to expand an MD simulation cell from a format   ", &
     "                (0) standard dl_poly CONFIG file                                ", &
     "                    title                                                   LINE", &
     "                    levcfg & imcon                                          LINE", &
     "                    lattice matrix - vector a components (Angstroms)        LINE", &
     "                    lattice matrix - vector b components (Angstroms)        LINE", &
     "                    lattice matrix - vector c components (Angstroms)        LINE", &
     "                    particle name (and global index)                        LINE", &
     "                    particle absolute coordinates (centred MD box origin)   LINE", &
     "                    optional particle velocity components                   LINE", &
     "                    optional particle force components                      LINE", &
     "                (1) xyz cell file                                               ", &
     "                    total number of particles                               LINE", &
     "                    title                                                   LINE", &
     "                    lattice matrix - vector a components (Angstroms)        LINE", &
     "                    lattice matrix - vector b components (Angstroms)        LINE", &
     "                    lattice matrix - vector c components (Angstroms)        LINE", &
     "                    particle name and absolute coordinates                  LINE", &
     "**************  authors: i.t.todorov                                            "

! Read command line

  Write(*,*)
  Write(*,'(1x,a)',Advance='No') 'Enter configuration filename: '
  record=' '
  Read(*,'(a)') record
  Call get_word(record,fin)

! Set time up

  Call gtime(t)

! proceeding

  Write(*,*)
  Write(*,'(1x,3a,i2)') 'You specified file: ', fin(1:Len_Trim(fin)) 
  Write(*,*)

! Create name for the XYZ configuration

  fout=' '
  fout=fin(1:Len_Trim(fin)) // ".XYZ"
  Write(*,'(1x,2a)') 'DLP2XYZ proceeding with the generation of ', fout(1:Len_Trim(fout))
  Write(*,*)

! Open channels

  nread=10
  Open(Unit=nread, File=fin, Form='formatted', Access='sequential')
  nrite=20
  Open(Unit=nrite, File=fout, Form='formatted', Access='direct', Recl=73)

! Read title

  cfgname=' '
  Read(Unit=nread, Fmt='(a)', End=20) cfgname

! Make sure cfgname is 72 characters long

  k=Len(cfgname)
  Do i=k,72
     cfgname(i:i)=' '
  End Do

! Define levcfg & imcon

  record=' '
  Read(Unit=nread, Fmt='(a)', End=20) record

  Call get_word(record,word)
  levcfg=Nint(word_2_real(word))

  Call get_word(record,word)
  imcon=Nint(word_2_real(word))

! Holt execution if levcfg or imcon is unsupported

  If ((levcfg < 0 .or. levcfg > 2) .or. (imcon < 0 .or. imcon > 7)) Go To 30

! Read lattice parameters

  record=' '
  Read(Unit=nread, Fmt='(a)', End=20) record
  Call get_word(record,word)
  cell(1)=word_2_real(word)
  Call get_word(record,word)
  cell(2)=word_2_real(word)
  Call get_word(record,word)
  cell(3)=word_2_real(word)

  record=' '
  Read(Unit=nread, Fmt='(a)', End=20) record
  Call get_word(record,word)
  cell(4)=word_2_real(word)
  Call get_word(record,word)
  cell(5)=word_2_real(word)
  Call get_word(record,word)
  cell(6)=word_2_real(word)

  record=' '
  Read(Unit=nread, Fmt='(a)', End=20) record
  Call get_word(record,word)
  cell(7)=word_2_real(word)
  Call get_word(record,word)
  cell(8)=word_2_real(word)
  Call get_word(record,word)
  cell(9)=word_2_real(word)

! lines to the first particle record

  n=2

! Set particle and line counters

  natms=0
  rec=1

! Write configuration file headers

  Write(*,'(1x,2a)')'File header: ',cfgname
  Write(*,*)
  rec=rec+1
  Write(Unit=nrite, Fmt='(72a,1a)',         Rec=rec) cfgname,lf

  Do While ( .true. )
     record=' '
     Read(Unit=nread, Fmt='(a)', End=10) record
     natms=natms+1
     Call get_word(record,name) ! atom label
     record=' '
     Read(Unit=nread, Fmt='(a)', End=20) record
     Call get_word(record,word) ! x coordinate
     xx0=word_2_real(word)
     Call get_word(record,word) ! y coordinate
     yy0=word_2_real(word)
     Call get_word(record,word) ! z coordinate
     zz0=word_2_real(word)
     If (levcfg > 0) Read(Unit=nread, Fmt=*, End=20)
     If (levcfg > 1) Read(Unit=nread, Fmt=*, End=20)

     rec=rec+1
     Write(Unit=nrite, Fmt='(a8,3f20.12,a4,1a)', Rec=rec) name,xx0,yy0,zz0,line(1:4),lf
  End Do

10 Continue

  Write(Unit=nrite, Fmt='(i10,a62,1a)', Rec=1) natms,line(1:62),lf
  Close (nrite)

! Time up

  Call gtime(t)
  Write(*,'(1x,a,f10.2,a)') 'Elapsed execution time', t, ' sec.'
  STOP

20 Continue

  Write(*,'(1x,a)') 'DLP2XYZ ABORTING... due to semantics failure'
  STOP

30 Continue

  Write(*,'(1x,a,i1,a,i1,a)') 'DLP2XYZ ABORTING... due to unacceptable levcfg/imcon values (',levcfg,'/',imcon,')'
  STOP

End Program dlp2xyz
