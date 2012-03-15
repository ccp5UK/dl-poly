Subroutine metal_ld_export(mdir,mlast,rho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export metal density data in domain boundary
! regions for halo formation
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,  Only : nrite,mxatms,mxbuff
  Use domains_module
  Use config_module, Only : ixyz

  Implicit None

  Integer,           Intent( In    ) :: mdir
  Integer,           Intent( InOut ) :: mlast
  Real( Kind = wp ), Intent( InOut ) :: rho(1:mxatms)

  Logical           :: safe
  Integer           :: fail,i,j,iblock,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  Allocate (buffer(1:mxbuff), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_export allocation failure, node: ', idnode
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
! kxyz - corrected halo reduction factor particles haloing both +&- sides
! jdnode - destination (send to), knode - source (receive from)

  kx = 0 ; ky = 0 ; kz = 0
  If      (mdir == -1) Then ! Direction -x
     kx  = 1
     jxyz= 1
     kxyz= 3

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= 3

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= 30

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= 30

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= 300

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= 300

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(47)
  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PARTICLES ON THIS NODE


  Do i=1,mlast

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

           If ((imove+2) <= iblock) Then

! pack particle density and halo indexing

              buffer(imove+1)=rho(i)
              buffer(imove+2)=Real(ixyz(i)-jxyz,wp)

              imove=imove+2

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
     Call warning(170,Real(itmp,wp),Real(mxbuff,wp),0.0_wp)
     Call error(38)
  End If

! exchange information on buffer sizes

  If (mxnode > 1) Then
     Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,Metldexp_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(imove,1,MPI_INTEGER,jdnode,Metldexp_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  Else
     jmove=imove
  End If

! Check for array bound overflow (can arrays cope with incoming data)

  safe=((mlast+jmove/2) <= mxatms)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=mlast+jmove/5
     If (mxnode > 1) Call gmax(itmp)
     Call warning(180,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
     Call error(39)
  End If

! exchange buffers between nodes (this is a MUST)

  If (mxnode > 1) Then
     Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,Metldexp_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(buffer(1),imove,wp_mpi,jdnode,Metldexp_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  End If

! load transferred data

  If (mxnode > 1) Then
     j=iblock
  Else
     j=0
  End If

  Do i=1,jmove/2
     mlast=mlast+1

! unpack particle density and remaining halo indexing

     rho(mlast) =buffer(j+1)
     ixyz(mlast)=Nint(buffer(j+2))

     j=j+2
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_export deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine metal_ld_export
