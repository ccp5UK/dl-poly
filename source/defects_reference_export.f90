Subroutine defects_reference_export(mdir,ixyz,nlrefs,namr,indr,xr,yr,zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic REFERENCE data in domain boundary
! regions for halo formation
!
! all particle coordinates are in reduced space with origin localised
! onto this node (idnode)
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,  Only : nrite,mxatms,mxbfxp
  Use domains_module

  Implicit None

  Integer,              Intent( In    ) :: mdir
  Integer,              Intent( InOut ) :: ixyz(1:mxatms),nlrefs
  Character( Len = 8 ), Intent( InOut ) :: namr(1:mxatms)
  Integer,              Intent( InOut ) :: indr(1:mxatms)
  Real( Kind = wp ),    Intent( InOut ) :: xr(1:mxatms),yr(1:mxatms),zr(1:mxatms)

  Logical           :: safe,lsx,lsy,lsz,lex,ley,lez
  Integer           :: fail,iadd,limit,iblock,          &
                       i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

! Number of transported quantities per particle

  iadd=13

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_export allocation failure, node: ', idnode
     Call error(0)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/Merge(2,1,mxnode > 1)

! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
! kxyz - corrected halo reduction factor particles haloing both +&- sides
! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
! jdnode - destination (send to), knode - source (receive from)

  kx = 0 ; ky = 0 ; kz = 0
  lsx = .false. ; lex = .false.
  lsy = .false. ; ley = .false.
  lsz = .false. ; lez = .false.
  If      (mdir == -1) Then ! Direction -x
     kx  = 1
     jxyz= 1
     kxyz= 3
     lsx = (idx == 0)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= 3
     lex = (idx == nprx-1)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= 30
     lsy = (idy == 0)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= 30
     ley = (idy == npry-1)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= 300
     lsz = (idz == 0)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= 300
     lez = (idz == nprz-1)

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(557)
  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,nlrefs

! If the particle is within the remaining 'inverted halo' of this domain

     If (ixyz(i) > 0) Then

! Get the necessary halo indices

        ix=Mod(ixyz(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

! If the particle is halo to both +&- sides
! use the corrected halo reduction factor

        If (j > jxyz .and. Mod(j,3) == 0) jxyz=kxyz

! If the particle is within the correct halo for the selected direction

        If (j == jxyz) Then

! If safe to proceed

           If ((imove+13) <= iblock) Then

! pack positions and apply possible wrap-around corrections for receiver

              buffer(imove+ 1)=xr(i)
              If (lsx) buffer(imove+1)=buffer(imove+1)+1.0_wp
              If (lex) buffer(imove+1)=buffer(imove+1)-1.0_wp
              buffer(imove+ 2)=yr(i)
              If (lsy) buffer(imove+2)=buffer(imove+2)+1.0_wp
              If (ley) buffer(imove+2)=buffer(imove+2)-1.0_wp
              buffer(imove+ 3)=zr(i)
              If (lsz) buffer(imove+3)=buffer(imove+3)+1.0_wp
              If (lez) buffer(imove+3)=buffer(imove+3)-1.0_wp

! pack config indexing, site name and remaining halo indexing arrays

              buffer(imove+ 4)=Real(indr(i),wp)
              buffer(imove+ 5)=Real(Ichar(namr(i)(1:1)),wp)
              buffer(imove+ 6)=Real(Ichar(namr(i)(2:2)),wp)
              buffer(imove+ 7)=Real(Ichar(namr(i)(3:3)),wp)
              buffer(imove+ 8)=Real(Ichar(namr(i)(4:4)),wp)
              buffer(imove+ 9)=Real(Ichar(namr(i)(5:5)),wp)
              buffer(imove+10)=Real(Ichar(namr(i)(6:6)),wp)
              buffer(imove+11)=Real(Ichar(namr(i)(7:7)),wp)
              buffer(imove+12)=Real(Ichar(namr(i)(8:8)),wp)
              buffer(imove+13)=Real(ixyz(i)-jxyz,wp)

           Else

              safe=.false.

           End If
           imove=imove+iadd

        End If

     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=Merge(2,1,mxnode > 1)*imove
     If (mxnode > 1) Call gmax(itmp)
     Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
     Call error(558)
  End If

! exchange information on buffer sizes

  If (mxnode > 1) Then
     Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,DefExport_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(imove,1,MPI_INTEGER,jdnode,DefExport_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  Else
     jmove=imove
  End If

! Check for array bound overflow (can arrays cope with incoming data)

  safe=((nlrefs+jmove/iadd) <= mxatms)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=nlrefs+jmove/iadd
     If (mxnode > 1) Call gmax(itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
     Call error(559)
  End If

! exchange buffers between nodes (this is a MUST)

  If (mxnode > 1) Then
     Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,DefExport_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(buffer(1),imove,wp_mpi,jdnode,DefExport_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  End If

! load transferred data

  j=Merge(iblock,0,mxnode > 1)
  Do i=1,jmove/13
     nlrefs=nlrefs+1

! unpack positions

     xr(nlrefs)=buffer(j+1)
     yr(nlrefs)=buffer(j+2)
     zr(nlrefs)=buffer(j+3)

! unpack config indexing, site name halo indexing arrays

     indr(nlrefs)=Nint(buffer(j+4))
     namr(nlrefs)(1:1)=Char(Nint(buffer(j+ 5)))
     namr(nlrefs)(2:2)=Char(Nint(buffer(j+ 6)))
     namr(nlrefs)(3:3)=Char(Nint(buffer(j+ 7)))
     namr(nlrefs)(4:4)=Char(Nint(buffer(j+ 8)))
     namr(nlrefs)(5:5)=Char(Nint(buffer(j+ 9)))
     namr(nlrefs)(6:6)=Char(Nint(buffer(j+10)))
     namr(nlrefs)(7:7)=Char(Nint(buffer(j+11)))
     namr(nlrefs)(8:8)=Char(Nint(buffer(j+12)))
     ixyz(nlrefs)=Nint(buffer(j+13))

     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_export deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine defects_reference_export
