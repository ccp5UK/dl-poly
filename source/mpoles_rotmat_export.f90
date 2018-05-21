Subroutine mpoles_rotmat_export(mdir,mlast,ixyz0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export multipoles rotation and infinitesimally
! rotated matrices in the halo
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,  Only : nrite,mximpl,mxatms,mxbfxp
  Use domains_module
  Use mpoles_module, Only : mplgfr,mprotx,mproty,mprotz

  Implicit None

  Integer, Intent( In    ) :: mdir
  Integer, Intent( InOut ) :: mlast,ixyz0(1:mxatms)


  Logical           :: safe
  Integer           :: fail,iadd,limit,iblock,          &
                       i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       idl1,idl2,idl3,idl4,             &
                       jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

! Number of transported quantities per particle

  iadd=4*mximpl+1
  idl1=mximpl ; idl2=2*mximpl ; idl3=3*mximpl ; idl4=4*mximpl

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
  Allocate (buffer(1:limit), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'export_atomic_positions allocation failure, node: ', idnode
     Call error(0)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=limit/Merge(2,1,mxnode > 1)

! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
! jdnode - destination (send to), kdnode - source (receive from)

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
     Call error(176)
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

     If (ixyz0(i) > 0) Then

! Get the necessary halo indices

        ix=Mod(ixyz0(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz0(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz0(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

! If safe to proceed

           If ((imove+iadd) <= iblock) Then

! pack rotation matrices and infinitely rotated matrices

              buffer(imove+1        : imove + idl1) = mplgfr(:,i)
              buffer(imove+1 + idl1 : imove + idl2) = mprotx(:,i)
              buffer(imove+1 + idl2 : imove + idl3) = mproty(:,i)
              buffer(imove+1 + idl3 : imove + idl4) = mprotz(:,i)

! Use the corrected halo reduction factor when the particle is halo to both +&- sides

              buffer(imove+iadd)=Real(ixyz0(i)-Merge(jxyz,kxyz,j == jxyz),wp)

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
     Call error(178)
  End If

! exchange information on buffer sizes

  If (mxnode > 1) Then
     Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,ExpMplRM_tag,dlp_comm_world,request,ierr)
     Call MPI_SEND(imove,1,MPI_INTEGER,jdnode,ExpMplRM_tag,dlp_comm_world,ierr)
     Call MPI_WAIT(request,status,ierr)
  Else
     jmove=imove
  End If

! Check for array bound overflow (can arrays cope with incoming data)

  safe=((mlast+jmove/iadd) <= mxatms)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=mlast+jmove/iadd
     If (mxnode > 1) Call gmax(itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
     Call error(180)
  End If

! exchange buffers between nodes (this is a MUST)

  If (mxnode > 1) Then
     If (jmove > 0) Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,ExpMplRM_tag,dlp_comm_world,request,ierr)
     If (imove > 0) Call MPI_SEND(buffer(1),imove,wp_mpi,jdnode,ExpMplRM_tag,dlp_comm_world,ierr)
     If (jmove > 0) Call MPI_WAIT(request,status,ierr)
  End If

! load transferred data

  j=Merge(iblock,0,mxnode > 1)
  Do i=1,jmove/iadd
     mlast=mlast+1

! unpack rotation matrices and infinitesimal rotation matrices

     mplgfr(:,mlast) = buffer(j+1        : j + idl1)
     mprotx(:,mlast) = buffer(j+1 + idl1 : j + idl2)
     mproty(:,mlast) = buffer(j+1 + idl2 : j + idl3)
     mprotz(:,mlast) = buffer(j+1 + idl3 : j + idl4)

     ixyz0(mlast)=Nint(buffer(j+iadd))
     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'mpoles_rotmat_export deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine mpoles_rotmat_export