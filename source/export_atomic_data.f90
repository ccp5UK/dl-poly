Subroutine export_atomic_data(mdir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic data in domain boundary regions
! for halo formation
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,  Only : nrite,mxatms,mxbfxp
  Use domains_module
  Use config_module, Only : nlast,ltg,lsite,ixyz,cell,xxx,yyy,zzz

  Implicit None

  Integer,           Intent( In    ) :: mdir

  Logical           :: safe,lsx,lsy,lsz,lex,ley,lez
  Integer           :: fail,iadd,limit,iblock,          &
                       i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       jdnode,kdnode,imove,jmove,itmp
  Real( Kind = wp ) :: uuu,vvv,www

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

! Number of transported quantities per particle

  iadd=6

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'export_atomic_data allocation failure, node: ', idnode
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
     Call error(46)
  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,nlast

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

           If ((imove+6) <= iblock) Then

! pack positions and apply possible wrap-around corrections for the receiver

              uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
              vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
              www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

              buffer(imove+1)=xxx(i)+cell(1)*uuu+cell(4)*vvv+cell(7)*www
              buffer(imove+2)=yyy(i)+cell(2)*uuu+cell(5)*vvv+cell(8)*www
              buffer(imove+3)=zzz(i)+cell(3)*uuu+cell(6)*vvv+cell(9)*www

! pack config indexing, site and remaining halo indexing arrays

              buffer(imove+4)=Real(ltg(i),wp)
              buffer(imove+5)=Real(lsite(i),wp)
              buffer(imove+6)=Real(ixyz(i)-jxyz,wp)

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
     Call error(54)
  End If

! exchange information on buffer sizes

  If (mxnode > 1) Then
     Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,Export_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(imove,1,MPI_INTEGER,jdnode,Export_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  Else
     jmove=imove
  End If

! Check for array bound overflow (can arrays cope with incoming data)

  safe=((nlast+jmove/iadd) <= mxatms)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=nlast+jmove/iadd
     If (mxnode > 1) Call gmax(itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
     Call error(56)
  End If

! exchange buffers between nodes (this is a MUST)

  If (mxnode > 1) Then
     Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,Export_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(buffer(1),imove,wp_mpi,jdnode,Export_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  End If

! load transferred data

  j=Merge(iblock,0,mxnode > 1)
  Do i=1,jmove/iadd
     nlast=nlast+1

! unpack positions

     xxx(nlast)=buffer(j+1)
     yyy(nlast)=buffer(j+2)
     zzz(nlast)=buffer(j+3)

! unpack config indexing, site and halo indexing arrays

     ltg(nlast)  =Nint(buffer(j+4))
     lsite(nlast)=Nint(buffer(j+5))
     ixyz(nlast) =Nint(buffer(j+6))

     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'export_atomic_data deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine export_atomic_data
