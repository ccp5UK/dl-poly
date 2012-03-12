Subroutine export_atomic_data(mdir,xxt,yyt,zzt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic data in domain boundary regions
! for halo formation
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module
  Use config_module, Only : nlast,ltg,lsite,ixyz

  Implicit None

  Integer,           Intent( In    ) :: mdir
  Real( Kind = wp ), Intent( InOut ) :: xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms)

  Logical           :: safe,lsx,lsy,lsz,lex,ley,lez
  Integer           :: fail,i,j,iblock,jxyz,ix,iy,iz,kx,ky,kz, &
                       jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  Allocate (buffer(1:mxbuff), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'export_atomic_data allocation failure, node: ', idnode
     Call error(0)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  If (mxnode > 1) Then
     iblock=mxbuff/2
  Else
     iblock=mxbuff
  End If

! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
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
     lsx = (idx == 0)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     lex = (idx == nprx-1)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     lsy = (idy == 0)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     ley = (idy == npry-1)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     lsz = (idz == 0)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
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

! If the particle is within the correct halo for the selected direction

        j=ix*kx+iy*ky+iz*kz
        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

! Reduce halo index (new remaining halo)
!
!           ixyz(i)=ixyz(i)-jxyz

! If safe to proceed

           If ((imove+6) <= iblock) Then

! pack positions and apply possible wrap-around corrections for receiver

              buffer(imove+1)=xxt(i)
              If (lsx) buffer(imove+1)=buffer(imove+1)+1.0_wp
              If (lex) buffer(imove+1)=buffer(imove+1)-1.0_wp
              buffer(imove+2)=yyt(i)
              If (lsy) buffer(imove+2)=buffer(imove+2)+1.0_wp
              If (ley) buffer(imove+2)=buffer(imove+2)-1.0_wp
              buffer(imove+3)=zzt(i)
              If (lsz) buffer(imove+3)=buffer(imove+3)+1.0_wp
              If (lez) buffer(imove+3)=buffer(imove+3)-1.0_wp

! pack config indexing, site and halo indexing arrays

              buffer(imove+4)=Real(ltg(i),wp)
              buffer(imove+5)=Real(lsite(i),wp)
              buffer(imove+6)=Real(ixyz(i),wp)

              imove=imove+6

           Else

              safe=.false.

           End If

        End If

     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     If (mxnode > 1) Then
        itmp=2*imove
     Else
        itmp=imove
     End If
     If (mxnode > 1) Call gmax(itmp)
     Call warning(150,Real(itmp,wp),Real(mxbuff,wp),0.0_wp)
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

  safe=((nlast+jmove/6) <= mxatms)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=nlast+jmove/6
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

  If (mxnode > 1) Then
     j=iblock
  Else
     j=0
  End If

  Do i=1,jmove/6
     nlast=nlast+1

! unpack positions

     xxt(nlast)=buffer(j+1)
     yyt(nlast)=buffer(j+2)
     zzt(nlast)=buffer(j+3)

! unpack config indexing, site and halo indexing arrays

     ltg(nlast)  =Nint(buffer(j+4))
     lsite(nlast)=Nint(buffer(j+5))
     ixyz(nlast) =Nint(buffer(j+6))

     j=j+6
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'export_atomic_data deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine export_atomic_data
