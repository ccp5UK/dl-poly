Program nfold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly utility to expand a triclinic MD simulation cell from format
!
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
!
! (1) general unit cell file
!     title                                                   LINE
!     lattice matrix - vector a components (Angstroms)        LINE
!     lattice matrix - vector b components (Angstroms)        LINE
!     lattice matrix - vector c components (Angstroms)        LINE
!     particle name and fractional coordinates                LINE
!
! (2) standard unit cell file
!     title                                                   LINE
!     lattice vectors' lengths (Angstroms) & angles (degrees) LINE
!     particle name and fractional coordinates                LINE
!
! (3) xyz cell file
!     total number of particles                               LINE
!     title                                                   LINE
!     lattice matrix - vector a components (Angstroms)        LINE
!     lattice matrix - vector b components (Angstroms)        LINE
!     lattice matrix - vector c components (Angstroms)        LINE
!     particle name and absolute coordinates                  LINE
!
! (type+10) indicates a hexagonal (type) cell file as INPUT ONLY!

! to format (0,1,2,3) by a nx*ny*nz volumetric multiplication of
! the MD cell contents along the MD cell vectors.  Hexagonal shaped
! cells (V) are only expanded to the nearest orthorhombic (2V)!!!
!
! f95 -o nfold.exe kinds_f90.f90 parse_module.f90 numeric_container.f90 nfold.f90
!
! copyright - daresbury laboratory
! authors   - i.t.todorov november 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use parse_module

  Implicit None

  Logical           :: exists,l_hex=.false.
  Character         :: cfgname*72,record*100,word*40,f_in*40,f_out*40,name*8,lf*1
  Integer           :: type_in,type_out,nx,ny,nz,nread,nrite,levcfg,imcon,n,natms,i,j,k,l,rec
  Real( Kind = wp ) :: half_minus,t,conv,xl,yl,zl,alpha,beta,gamma,cosg,sing,cosa,cosb,z2,z3, &
                       cell(1:9),cel0(1:9),celh(1:9),celprp(1:10),rcell(1:9),fx,fy,fz, &
                       wx,wy,wz,f1,f2,f3,f4,f5,f6,f7,f8,f9,xx0,yy0,zz0,xx,yy,zz

! define half-

  half_minus = Nearest(0.5_wp,-1.0_wp)

! define newline character Char(10) or Char(Iachar('\n'))

  lf = Char(10)

! Print header

  Write(*,'(38(1x,a,/))')                                                                &
     "**************  dl_poly utility to expand a triclinic atom cell from format     ", &
     "                                                                                ", &
     "                (0) standard dl_poly CONFIG file                                ", &
     "                    title                                                   LINE", &
     "                    levcfg & imcon                                          LINE", &
     "                    lattice matrix - vector a components (Angstroms)        LINE", &
     "                    lattice matrix - vector b components (Angstroms)        LINE", &
     "                    lattice matrix - vector c components (Angstroms)        LINE", &
     "      NN            particle name (and global index)                        LINE", &
     "      NN            particle absolute coordinates (centred MD box origin)   LINE", &
     "      ||            optional particle velocity components                   LINE", &
     "      ||            optional particle force components                      LINE", &
     "      FF                                                                        ", &
     "      FF        (1) general unit cell file                                      ", &
     "      ||            title                                                   LINE", &
     "      ||            lattice matrix - vector a components (Angstroms)        LINE", &
     "      OO            lattice matrix - vector b components (Angstroms)        LINE", &
     "      OO            lattice matrix - vector c components (Angstroms)        LINE", &
     "      ||            particle name and fractional coordinates                LINE", &
     "      ||                                                                        ", &
     "      LL        (2) standard unit cell file                                     ", &
     "      LL            title                                                   LINE", &
     "      ||            lattice vectors' lengths (Angstroms) & angles (degrees) LINE", &
     "      ||            particle name and fractional coordinates                LINE", &
     "      DD                                                                        ", &
     "      DD        (3) xyz cell file                                               ", &
     "                    total number of particles                               LINE", &
     "                    title                                                   LINE", &
     "                    lattice matrix - vector a components (Angstroms)        LINE", &
     "                    lattice matrix - vector b components (Angstroms)        LINE", &
     "                    lattice matrix - vector c components (Angstroms)        LINE", &
     "                    particle name and absolute coordinates                  LINE", &
     "                                                                                ", &
     "                (type+10) indicates a hexagonal (type) cell file as INPUT ONLY! ", &
     "                                                                                ", &
     "                to format (0,1,2,3) by a nx*ny*nz volumetric multiplication of  ", &
     "                the MD cell contents along the MD cell vectors. Hexagonal shaped", &
     "                cells (V) are only expanded to the nearest orthorhombic (2V)!!! ", &
     "**************  author: i.t.todorov                                             "

! Read command line

  Write(*,*)
  Write(*,'(1x,a,/,1x,a)') 'Enter space separated input configuration filename and', &
                           'format (0,1,2,3,10,11,12,13): '
  record=' '
  Do While (Len_Trim(record) == 0)
     Write(*,'(1x,a)',Advance='No')
     Read(*,'(a)') record
  End Do
  Call get_word(record,f_in)
  Call get_word(record,word)
  Do While (Len_Trim(word) == 0)
     Write(*,'(3a)',Advance='No') ' ',f_in(1:Len_Trim(f_in)),' '
     record=' '
     Read(*,'(a)') record
     Call get_word(record,word)
  End Do
  Write(*,*)

! Check file format

  Write(*,'(1x,3a)') 'NFOLD will abort if ', word(1:Len_Trim(word)), ' is NaN!!!'
  type_in=Nint(word_2_real(word))
  If (.not.((type_in > 0 .and. type_in < 4) .or. (type_in > 9 .and. type_in < 14))) Then
     Write(*,'(1x,3a)') 'NFOLD ABORTING... format ', word(1:Len_Trim(word)), ' is not implemented by NFOLD'
     STOP
  Else
     If (type_in > 9) Then
        l_hex=.true.

        type_in=type_in-10
     End If
     If (type_in == 3) Then
        Write(*,'(1x,3a)') 'NFOLD expects bounding box lattice vectors for XYZ format!!!'
        Write(*,*)
     End If
  End If
  Write(*,*)

! Check file existence

  Inquire(File=f_in(1:Len_Trim(f_in)),Exist=exists)
  If (exists) Then
     Write(*,'(1x,4a)') 'The program will abort if ', f_in(1:Len_Trim(f_in)), &
     ' format does not match the specified ', word(1:Len_Trim(word))
  Else
     Write(*,'(1x,3a)') 'NFOLD ABORTING... ', f_in(1:Len_Trim(f_in)), ' does not exist in the current directory!!!'
     STOP
  End If
  Write(*,*)

! Read command line

  Write(*,'(1x,a,/,1x,a)') 'Enter space separated output format (0,1,2,3) and ', &
                           'optional configuration filename: '

  record=' '
  Do While (Len_Trim(record) == 0)
     Write(*,'(1x,a)',Advance='No')
     Read(*,'(a)') record
  End Do
  Call get_word(record,word)

! Check file format

  Write(*,'(1x,3a)') 'NFOLD will abort if ', word(1:Len_Trim(word)), ' is NaN!!!'
  type_out=Nint(word_2_real(word))
  If (type_out < 0 .or. type_out > 3) Then
     Write(*,'(1x,3a)') 'NFOLD ABORTING... format ', word(1:Len_Trim(word)), ' is not implemented by NFOLD!!!'
     STOP
  End If
  Write(*,*)

  Call get_word(record,f_out)
  If (Len_Trim(f_out) == 0) Write(*,'(1x,a)') 'NFOLD will generate an aoutomatic file name!!!'
  Write(*,*)

! Read command line

  Write(*,'(1x,a)',Advance='No') 'Enter three space separated positive integers (nx,ny,nz): '
  record=' '
  Read(*,'(a)') record
  Write(*,*)
  Write(*,'(1x,a)') 'Specifying no integers defaults to 1x1x1 multiplication.'

  Call get_word(record,word)
  nx=Max(1,Abs(Nint(word_2_real(word))))
  fx=Real(nx,wp)
  Call get_word(record,word)
  ny=Max(1,Abs(Nint(word_2_real(word))))
  fy=Real(ny,wp)
  Call get_word(record,word)
  nz=Max(1,Abs(Nint(word_2_real(word))))
  fz=Real(nz,wp)

! Open input channel

  nread=222
  Open(Unit=nread, File=f_in, Form='formatted', Access='sequential')

  If (type_in == 3) Then
     record=' '
     Read(Unit=nread, Fmt='(a)', End=50) record
     n=Abs(Nint(word_2_real(record)))
  End if

! Read title

  cfgname=' '
  Read(Unit=nread, Fmt='(a)', End=50) cfgname

! Make sure cfgname is 72 characters long

  k=Len(cfgname)
  Do i=k,72
     cfgname(i:i)=' '
  End Do

! Define levcfg & imcon

  levcfg=0
  imcon=3
  If (type_in == 0) Then
     record=' '
     Read(Unit=nread, Fmt='(a)', End=50) record

     Call get_word(record,word)
     levcfg=Nint(word_2_real(word))

     Call get_word(record,word)
     imcon=Nint(word_2_real(word))

! Holt execution if levcfg or imcon is unsupported

     If ((levcfg < 0 .or. levcfg > 2) .or. &
         (imcon <= 0 .or. imcon == 4 .or. imcon ==5 .or. imcon > 7)) Then
        Write(*,'(1x,2(a,i1),a)') 'NFOLD ABORTING... due to unacceptable levcfg/imcon values (',levcfg,'/',imcon,')'
        STOP
     End If
     If (imcon == 6 .and. nz > 1) Then
        nz=1
        Write(*,'(1x,a)') 'Multiplication nz defaults to 1 for slabs !!!'
        Write(*,*)
     End if
  End If

  If (l_hex .or. (type_in == 0 .and. imcon == 7)) Then
     cel0=0.0_wp
     celh=0.0_wp

     Write(*,*)
     Write(*,'(1x,a)') 'Hexagonal to orthorhombic cell transformation will take priority!!!'
     Write(*,'(1x,a)') 'V(orthorhombic) = 2*V(hexagonal) and particle count will double!!!'
     Write(*,*)

     nx=1 ; ny=1 ; nz=1
     fx=1.0_wp ; fy=1.0_wp ; fz=1.0_wp
  End If

  conv=Atan(1.0_wp)/45.0_wp  ! degree to radians conversion

! Read lattice parameters

  If (type_in == 2) Then
     record=' '
     Read(Unit=nread, Fmt='(a)', End=50) record
     Call get_word(record,word)
     xl=word_2_real(word)  ! length a
     Call get_word(record,word)
     yl=word_2_real(word)  ! length b
     Call get_word(record,word)
     zl=word_2_real(word)  ! length c
     Call get_word(record,word)
     alpha=word_2_real(word)  ! angle bc
     Call get_word(record,word)
     beta =word_2_real(word)  ! angle ac
     Call get_word(record,word)
     gamma=word_2_real(word)  ! angle ab

     cosg=Cos(gamma*conv)
     sing=Sin(gamma*conv)
     cosb=Cos(beta*conv)
     cosa=Cos(alpha*conv)

     z2=(cosa-cosb*cosg)/sing
     z3=Sqrt(1.0_wp-cosb**2-z2**2)
!    z3=Sqrt(1.0_wp-cosa**2-cosb**2-cosg**2+2.0_wp*cosa*cosb*cosg)/sing

     cell(1)=xl
     cell(2)=0.0_wp
     cell(3)=0.0_wp

     cell(4)=yl*cosg
     cell(5)=yl*sing
     cell(6)=0.0_wp

     cell(7)=zl*cosb
     cell(8)=zl*z2
     cell(9)=zl*z3

! Catch hex2orth conversion

     If ( (.not.l_hex) .and. ( (Nint(gamma) == 120) .and.           &
           (Nint(alpha) == 90) .and. (alpha-beta < 5.0e-2_wp) .and. &
           (Abs(4.0_wp*alpha-3.0_wp*gamma) < 5.0e-2_wp) ) ) Then
        l_hex=.true.

        cel0=0.0_wp
        celh=0.0_wp

        Write(*,*)
        Write(*,'(1x,a)') 'Hexagonal to orthorhombic cell transformation will take priority!!!'
        Write(*,'(1x,a)') 'V(orthorhombic) = 2*V(hexagonal) and particle count will double!!!'
        Write(*,*)

        nx=1 ; ny=1 ; nz=1
        fx=1.0_wp ; fy=1.0_wp ; fz=1.0_wp
     End If

     If (l_hex) Then
        cel0(1)=3.0_wp*xl
        cel0(5)=Sqrt(3.0_wp)*xl
        cel0(9)=zl

        celh=cell
     End If
  Else
     record=' '
     Read(Unit=nread, Fmt='(a)', End=50) record
     Call get_word(record,word)
     cell(1)=word_2_real(word)
     Call get_word(record,word)
     cell(2)=word_2_real(word)
     Call get_word(record,word)
     cell(3)=word_2_real(word)

     record=' '
     Read(Unit=nread, Fmt='(a)', End=50) record
     Call get_word(record,word)
     cell(4)=word_2_real(word)
     Call get_word(record,word)
     cell(5)=word_2_real(word)
     Call get_word(record,word)
     cell(6)=word_2_real(word)

     record=' '
     Read(Unit=nread, Fmt='(a)', End=50) record
     Call get_word(record,word)
     cell(7)=word_2_real(word)
     Call get_word(record,word)
     cell(8)=word_2_real(word)
     Call get_word(record,word)
     cell(9)=word_2_real(word)

! Get cell descriptors

     Call dcell(cell,celprp)
     xl=celprp(1) ! length a
     yl=celprp(2) ! length b
     zl=celprp(3) ! length c
     alpha=Acos(celprp(6))/conv  ! angle bc
     beta =Acos(celprp(5))/conv  ! angle ac
     gamma=Acos(celprp(4))/conv  ! angle ab
     wx=celprp(7) ! width along x axis
     wy=celprp(8) ! width along y axis
     wz=celprp(9) ! width along z axis

     If ( (.not.l_hex) .and. ( (Nint(gamma) == 120) .and.           &
           (Nint(alpha) == 90) .and. (alpha-beta < 5.0e-2_wp) .and. &
           (Abs(4.0_wp*alpha-3.0_wp*gamma) < 5.0e-2_wp) ) ) Then
        l_hex=.true.

        cel0=0.0_wp
        celh=0.0_wp

        Write(*,*)
        Write(*,'(1x,a)') 'Hexagonal to orthorhombic cell transformation will take priority!!!'
        Write(*,'(1x,a)') 'V(orthorhombic) = 2*V(hexagonal) and particle count will double!!!'
        Write(*,*)

        nx=1 ; ny=1 ; nz=1
        fx=1.0_wp ; fy=1.0_wp ; fz=1.0_wp
     End If

     If (l_hex) Then
        If (type_in == 0) Then
           cel0=cell
           celh=cell

           cell(1)=cell(1)/3.0_wp
           cell(2)=0.0_wp
           cell(3)=0.0_wp

           cell(4)=cell(1)*Cos(conv*120.0_wp)
           cell(5)=cell(1)*Sin(conv*120.0_wp)
           cell(6)=0.0_wp

           cell(7)=0.0_wp
           cell(8)=0.0_wp
           cell(9)=cell(9)
        Else
           cel0(1)=3.0_wp*cell(1)
           cel0(5)=Sqrt(3.0_wp)*cell(1)
           cel0(9)=cell(9)

           celh=cell
        End If

! Correct cell descriptors for reporting

        xl=cell(1)        ! length a
        yl=xl             ! length b
        zl=cell(9)        ! length c
        alpha=90.0_wp     ! angle bc
        beta =90.0_wp     ! angle ac
        gamma=120.0_wp    ! angle ab
        wx=0.5_wp*cell(5) ! width along x axis
        wy=yl             ! width along y axis
        wz=zl             ! width along z axis
     End If
  End If

  If (l_hex) Then

! Create name for the expanded configuration

     If (Len_Trim(f_out) == 0) f_out=f_in(1:Len_Trim(f_in)) // '_hex2orth'

     Call invert(cel0(1:9),rcell,celprp(10))

! Set time up

     Call gtime(t)

! Open output channel

     nrite=555
     Open(Unit=nrite, File=f_out, Form='formatted', Access='direct', Recl=73)

! Set particle and line counters

     natms=0
     rec=0

! Write configuration file headers

     Write(*,'(1x,2a)') 'File header: ',cfgname
     Write(*,*)
     If (type_out == 3) rec=rec+1
     rec=rec+1
     Write(Unit=nrite, Fmt='(72a,1a)', Rec=rec) cfgname,lf
     If (type_out == 0) Then
        rec=rec+1
        Write(Unit=nrite, Fmt='(2i10,a52,1a)', Rec=rec) 0,3,Repeat(' ',52),lf
     End If
     If (type_out == 2) Then
        rec=rec+1
        Write(Unit=nrite, Fmt='(6f12.7,1a)', Rec=rec) cel0(1),cel0(5),cel0(9),90.0_wp,90.0_wp,90.0_wp,lf
     Else
        rec=rec+1
        Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) cel0(1),cel0(2),cel0(3),Repeat(' ',12),lf
        rec=rec+1
        Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) cel0(4),cel0(5),cel0(6),Repeat(' ',12),lf
        rec=rec+1
        Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) cel0(7),cel0(8),cel0(9),Repeat(' ',12),lf
     End If

     Do
        record=' '
        Read(Unit=nread, Fmt='(a)', End=20) record
        natms=natms+1
        Call get_word(record,name) ! atom label
        If (type_in == 0) Then
           record=' '
           Read(Unit=nread, Fmt='(a)', End=50) record
        End If
        Call get_word(record,word) ! x coordinate
        xx0=word_2_real(word)
        Call get_word(record,word) ! y coordinate
        yy0=word_2_real(word)
        Call get_word(record,word) ! z coordinate
        zz0=word_2_real(word)
        If (type_in == 0 .and. levcfg > 0) Read(Unit=nread, Fmt=*, End=50)
        If (type_in == 0 .and. levcfg > 1) Read(Unit=nread, Fmt=*, End=50)

! Convert fractional to absolute

        If (Mod(type_in,3) /= 0) Then
           xx0=xx0-Anint(xx0) ; If (xx0 >= half_minus) xx0=-xx0
           yy0=yy0-Anint(yy0) ; If (yy0 >= half_minus) yy0=-yy0
           zz0=zz0-Anint(zz0) ; If (zz0 >= half_minus) zz0=-zz0

           xx=celh(1)*xx0+celh(4)*yy0+celh(7)*zz0
           yy=celh(2)*xx0+celh(5)*yy0+celh(8)*zz0
           zz=celh(3)*xx0+celh(6)*yy0+celh(9)*zz0

           xx0=xx
           yy0=yy
           zz0=zz
        End If

        xx=xx0
        yy=yy0
        zz=zz0

! Convert absolute to fractional

        If (Mod(type_out,3) /= 0) Then
           xx0=xx
           yy0=yy
           zz0=zz

           xx=rcell(1)*xx0+rcell(4)*yy0+rcell(7)*zz0
           yy=rcell(2)*xx0+rcell(5)*yy0+rcell(8)*zz0
           zz=rcell(3)*xx0+rcell(6)*yy0+rcell(9)*zz0

           xx=xx-Anint(xx) ; If (xx >= half_minus) xx=-xx
           yy=yy-Anint(yy) ; If (yy >= half_minus) yy=-yy
           zz=zz-Anint(zz) ; If (zz >= half_minus) zz=-zz
        End If

        If (type_out == 0) Then
           rec=rec+1
           Write(Unit=nrite, Fmt='(a8,i10,a54,a1)', Rec=rec) name,natms,Repeat(' ',54),lf
           rec=rec+1
           Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) xx,yy,zz,Repeat(' ',12),lf
        Else
           rec=rec+1
           Write(Unit=nrite, Fmt='(a8,3f20.12,a4,a1)', Rec=rec) name,xx,yy,zz,Repeat(' ',4),lf
        End If
     End Do

20   Continue
     Rewind(Unit=nread, Err=50)
     If      (type_in == 0 .or. type_in == 3) Then
        Do l=1,5
           Read(Unit=nread, Fmt=*, End=50)
        End Do
     Else If (type_in == 1) Then
        Do l=1,4
           Read(Unit=nread, Fmt=*, End=50)
        End Do
     Else If (type_in == 2) Then
        Do l=1,2
           Read(Unit=nread, Fmt=*, End=50)
        End Do
     End If

! Second pass

     Do
        record=' '
        Read(Unit=nread, Fmt='(a)', End=30) record
        natms=natms+1
        Call get_word(record,name) ! atom label
        If (type_in == 0) Then
           record=' '
           Read(Unit=nread, Fmt='(a)', End=50) record
        End If
        Call get_word(record,word) ! x coordinate
        xx0=word_2_real(word)
        Call get_word(record,word) ! y coordinate
        yy0=word_2_real(word)
        Call get_word(record,word) ! z coordinate
        zz0=word_2_real(word)
        If (type_in == 0 .and. levcfg > 0) Read(Unit=nread, Fmt=*, End=50)
        If (type_in == 0 .and. levcfg > 1) Read(Unit=nread, Fmt=*, End=50)

! Convert fractional to absolute

        If (Mod(type_in,3) /= 0) Then
           xx0=xx0-Anint(xx0) ; If (xx0 >= half_minus) xx0=-xx0
           yy0=yy0-Anint(yy0) ; If (yy0 >= half_minus) yy0=-yy0
           zz0=zz0-Anint(zz0) ; If (zz0 >= half_minus) zz0=-zz0

           xx=celh(1)*xx0+celh(4)*yy0+celh(7)*zz0
           yy=celh(2)*xx0+celh(5)*yy0+celh(8)*zz0
           zz=celh(3)*xx0+celh(6)*yy0+celh(9)*zz0

           xx0=xx
           yy0=yy
           zz0=zz
        End If

! Add displacements

        xx0=xx0+cel0(1)/2
        yy0=yy0+cel0(5)/2
        zz0=zz0

! Wrap up

        xx=rcell(1)*xx0+rcell(4)*yy0+rcell(7)*zz0
        yy=rcell(2)*xx0+rcell(5)*yy0+rcell(8)*zz0
        zz=rcell(3)*xx0+rcell(6)*yy0+rcell(9)*zz0

        xx=xx-Anint(xx) ; If (xx >= half_minus) xx=-xx
        yy=yy-Anint(yy) ; If (yy >= half_minus) yy=-yy
        zz=zz-Anint(zz) ; If (zz >= half_minus) zz=-zz

! Convert fractional to absolute

        xx0=xx
        yy0=yy
        zz0=zz

        xx=cel0(1)*xx0+cel0(4)*yy+cel0(7)*zz
        yy=cel0(2)*xx0+cel0(5)*yy+cel0(8)*zz
        zz=cel0(3)*xx0+cel0(6)*yy+cel0(9)*zz

! Convert absolute to fractional

        If (Mod(type_out,3) /= 0) Then
           xx0=xx
           yy0=yy
           zz0=zz

           xx=rcell(1)*xx0+rcell(4)*yy0+rcell(7)*zz0
           yy=rcell(2)*xx0+rcell(5)*yy0+rcell(8)*zz0
           zz=rcell(3)*xx0+rcell(6)*yy0+rcell(9)*zz0

           xx=xx-Anint(xx) ; If (xx >= half_minus) xx=-xx
           yy=yy-Anint(yy) ; If (yy >= half_minus) yy=-yy
           zz=zz-Anint(zz) ; If (zz >= half_minus) zz=-zz
        End If

        If (type_out == 0) Then
           rec=rec+1
           Write(Unit=nrite, Fmt='(a8,i10,a54,a1)', Rec=rec) name,natms,Repeat(' ',54),lf
           rec=rec+1
           Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) xx,yy,zz,Repeat(' ',12),lf
        Else
           rec=rec+1
           Write(Unit=nrite, Fmt='(a8,3f20.12,a4,a1)', Rec=rec) name,xx,yy,zz,Repeat(' ',4),lf
        End If
     End Do

30   Continue

     Close(Unit=nread)

     If (type_out == 0) Write(Unit=nrite, Fmt='(3i10,a42,1a)', Rec=2) 0,imcon,natms,Repeat(' ',42),lf
     If (type_out == 3) Write(Unit=nrite, Fmt='(i10,a62,1a)', Rec=1) natms,Repeat(' ',62),lf
     Close (Unit=nrite)

     If (type_in == 3 .and. natms /= 2*n) Go To 50

! Write summary data

     Write(*,'(1x,2a,i3)') 'INPUT NAME & FORMAT: ', f_in(1:Len_Trim(f_in)),type_in
     Write(*,'(1x,a,6f8.3)') "DIMENSIONS':",xl,yl,zl,alpha,beta,gamma
     Write(*,'(1x,a,i15,f8.3)') 'SIZE & MAXIMIM CUTOFF RADIUS:',natms/2,0.5_wp*Min(wx,wy,wz)
     Write(*,*)
     Write(*,'(1x,2a,i3)') 'OUTPUT NAME & FORMAT: ', f_out(1:Len_Trim(f_out)),type_out
     Write(*,'(1x,a,6f8.3)') 'DIMENSIONS:',cel0(1),cel0(5),cel0(9),90.0_wp,90.0_wp,90.0_wp
     Write(*,'(1x,a,i15,f8.3)') 'SIZE & MAXIMIM CUTOFF RADIUS:',natms,0.5_wp*Min(cel0(1),cel0(5),cel0(9))
     Write(*,*)

  Else

     If (Mod(type_out,3) /= 0) Then
        celprp(1)=fx*cell(1)
        celprp(2)=fx*cell(2)
        celprp(3)=fx*cell(3)

        celprp(4)=fy*cell(4)
        celprp(5)=fy*cell(5)
        celprp(6)=fy*cell(6)

        celprp(7)=fz*cell(7)
        celprp(8)=fz*cell(8)
        celprp(9)=fz*cell(9)

        Call invert(celprp(1:9),rcell,celprp(10))
     End If

! Create name for the expanded configuration

     If (Len_Trim(f_out) == 0) Then
        record=' '
        word=' ' ; Write(word,'(i5)') nx ; Call strip_blanks(word)
        record=record(1:Len_Trim(record)) // '_' // word(1:Len_Trim(word))

        word=' ' ; Write(word,'(i5)') ny ; Call strip_blanks(word)
        record=record(1:Len_Trim(record)) // '_' // word(1:Len_Trim(word))

        word=' ' ; Write(word,'(i5)') nz ; Call strip_blanks(word)
        record=record(1:Len_Trim(record)) // '_' // word(1:Len_Trim(word))

        f_out=' '
        f_out=f_in(1:Len_Trim(f_in)) // record(1:Len_Trim(record))
     End If
     Write(*,'(1x,2a)') 'NFOLD proceeding with the generation of ', f_out(1:Len_Trim(f_out))
     Write(*,*)

! Set time up

     Call gtime(t)

! Open output channel

     nrite=20
     Open(Unit=nrite, File=f_out, Form='formatted', Access='direct', Recl=73)

! Set particle and line counters

     natms=0
     rec=0

     Do k=1,nz

! Define cell vector displacement in z direction

        f7=cell(7)*Real(2*k-nz-1,wp)/2.0_wp
        f8=cell(8)*Real(2*k-nz-1,wp)/2.0_wp
        f9=cell(9)*Real(2*k-nz-1,wp)/2.0_wp

        Do j=1,ny

! Define cell vector displacement in y direction

           f4=cell(4)*Real(2*j-ny-1,wp)/2.0_wp
           f5=cell(5)*Real(2*j-ny-1,wp)/2.0_wp
           f6=cell(6)*Real(2*j-ny-1,wp)/2.0_wp

           Do i=1,nx

! Define cell vector displacement in x direction

              f1=cell(1)*Real(2*i-nx-1,wp)/2.0_wp
              f2=cell(2)*Real(2*i-nx-1,wp)/2.0_wp
              f3=cell(3)*Real(2*i-nx-1,wp)/2.0_wp

! Write configuration file headers

              If (i == 1 .and. j == 1 .and. k == 1) Then
                 Write(*,'(1x,2a)') 'File header: ',cfgname
                 Write(*,*)
                 If (type_out == 3) rec=rec+1
                 rec=rec+1
                 Write(Unit=nrite, Fmt='(72a,1a)', Rec=rec) cfgname,lf
                 If (type_out == 0) Then
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(2i10,a52,1a)', Rec=rec) 0,imcon,Repeat(' ',52),lf
                 End If
                 If (type_out == 2) Then
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(6f12.7,1a)', Rec=rec) xl,yl,zl,alpha,beta,gamma,lf
                 Else
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) fx*cell(1),fx*cell(2),fx*cell(3),Repeat(' ',12),lf
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) fy*cell(4),fy*cell(5),fy*cell(6),Repeat(' ',12),lf
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) fz*cell(7),fz*cell(8),fz*cell(9),Repeat(' ',12),lf
                 End If
              End If

              Do
                 record=' '
                 Read(Unit=nread, Fmt='(a)', End=40) record
                 natms=natms+1
                 Call get_word(record,name) ! atom label
                 If (type_in == 0) Then
                    record=' '
                    Read(Unit=nread, Fmt='(a)', End=50) record
                 End If
                 Call get_word(record,word) ! x coordinate
                 xx0=word_2_real(word)
                 Call get_word(record,word) ! y coordinate
                 yy0=word_2_real(word)
                 Call get_word(record,word) ! z coordinate
                 zz0=word_2_real(word)
                 If (type_in == 0 .and. levcfg > 0) Read(Unit=nread, Fmt=*, End=50)
                 If (type_in == 0 .and. levcfg > 1) Read(Unit=nread, Fmt=*, End=50)

! Convert fractional to absolute

                 If (Mod(type_in,3) /= 0) Then
                    xx0=xx0-Anint(xx0) ; If (xx0 >= half_minus) xx0=-xx0
                    yy0=yy0-Anint(yy0) ; If (yy0 >= half_minus) yy0=-yy0
                    zz0=zz0-Anint(zz0) ; If (zz0 >= half_minus) zz0=-zz0

                    xx=cell(1)*xx0+cell(4)*yy0+cell(7)*zz0
                    yy=cell(2)*xx0+cell(5)*yy0+cell(8)*zz0
                    zz=cell(3)*xx0+cell(6)*yy0+cell(9)*zz0

                    xx0=xx
                    yy0=yy
                    zz0=zz
                 End If

! Create replicas of each atomic coordinate

                 xx=xx0+f1+f4+f7
                 yy=yy0+f2+f5+f8
                 zz=zz0+f3+f6+f9

! Convert absolute to fractional

                 If (Mod(type_out,3) /= 0) Then
                    xx0=xx
                    yy0=yy
                    zz0=zz

                    xx=rcell(1)*xx0+rcell(4)*yy0+rcell(7)*zz0
                    yy=rcell(2)*xx0+rcell(5)*yy0+rcell(8)*zz0
                    zz=rcell(3)*xx0+rcell(6)*yy0+rcell(9)*zz0

                    xx=xx-Anint(xx) ; If (xx >= half_minus) xx=-xx
                    yy=yy-Anint(yy) ; If (yy >= half_minus) yy=-yy
                    zz=zz-Anint(zz) ; If (zz >= half_minus) zz=-zz
                 End If

                 If (type_out == 0) Then
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(a8,i10,a54,a1)', Rec=rec) name,natms,Repeat(' ',54),lf
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(3f20.12,a12,1a)', Rec=rec) xx,yy,zz,Repeat(' ',12),lf
                 Else
                    rec=rec+1
                    Write(Unit=nrite, Fmt='(a8,3f20.12,a4,a1)', Rec=rec) name,xx,yy,zz,Repeat(' ',4),lf
                 End If
              End Do

40            Continue
              Rewind(Unit=nread, Err=50)
              If      (type_in == 0 .or. type_in == 3) Then
                 Do l=1,5
                    Read(Unit=nread, Fmt=*, End=50)
                 End Do
              Else If (type_in == 1) Then
                 Do l=1,4
                    Read(Unit=nread, Fmt=*, End=50)
                 End Do
              Else If (type_in == 2) Then
                 Do l=1,2
                    Read(Unit=nread, Fmt=*, End=50)
                 End Do
              End If

           End Do

        End Do

     End Do
     Close(Unit=nread)

     If (type_out == 0) Write(Unit=nrite, Fmt='(3i10,a42,1a)', Rec=2) 0,imcon,natms,Repeat(' ',42),lf
     If (type_out == 3) Write(Unit=nrite, Fmt='(i10,a62,1a)', Rec=1) natms,Repeat(' ',62),lf
     Close (Unit=nrite)

     If (type_in == 3 .and. natms /= n) Go To 50

! Write summary data

     Write(*,'(1x,2a,i3)') 'INPUT NAME & FORMAT: ', f_in(1:Len_Trim(f_in)),type_in
     Write(*,'(1x,a,6f8.3)') 'DIMENSIONS:',xl,yl,zl,alpha,beta,gamma
     Write(*,'(1x,a,i15,f8.3)') 'SIZE & MAXIMIM CUTOFF RADIUS:',natms/(nx*ny*nz),0.5_wp*Min(wx,wy,wz)
     Write(*,*)
     Write(*,'(1x,2a,i3)') 'OUTPUT NAME & FORMAT: ', f_out(1:Len_Trim(f_out)),type_out
     Write(*,'(1x,a,6f8.3)') 'DIMENSIONS:',Real(nx,wp)*xl,Real(ny,wp)*yl,Real(nz,wp)*zl,alpha,beta,gamma
     Write(*,'(1x,a,i15,f8.3)') 'SIZE & MAXIMIM CUTOFF RADIUS:',natms,0.5_wp*Min(Real(nx,wp)*wx,Real(ny,wp)*wy,Real(nz,wp)*wz)
     Write(*,*)

  End If

  If (type_out == 3) Then
     Write(*,'(1x,3a)') 'NFOLD produces a bounding box lattice vectors for XYZ format!!!'
     Write(*,*)
  End If

! Time up

  Call gtime(t)
  Write(*,'(1x,a,f10.2,a)') 'Elapsed execution time', t, ' sec.'
  STOP

50 Continue

  Write(*,'(1x,a)') 'NFOLD ABORTING... due to semantics failure'
  STOP

End Program nfold
