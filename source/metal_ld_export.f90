Subroutine metal_ld_export(mdir,sidex,sidey,sidez,cwx,cwy,cwz,mlast,iwrk,rho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic density data in domain boundary
! regions
!
! all particle coordinates are in reduced space with origin localised
! onto this node (idnode)
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module
  Use config_module, Only : xxx,yyy,zzz

  Implicit None

  Integer,           Intent( In    ) :: mdir
  Real( Kind = wp ), Intent( In    ) :: sidex,sidey,sidez,cwx,cwy,cwz
  Integer,           Intent( InOut ) :: mlast,iwrk(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: rho(1:mxatms)

  Logical           :: safe
  Integer           :: fail,i,j,iblock,jdnode,kdnode,imove,jmove,itmp
  Real( Kind = wp ) :: shovex,shovey,shovez,begin,final,xyz

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

! define the relative displacement between the coordinate systems'
! origins of neighbouring domains with respect to this domain
! (idnode) in the direction of mdir (shovex,shovey,shovez)
! in order to push all particles of this domain (idnode) in the
! direction of mdir

! define 'minus halo' limits in the direction of mdir (begin,final)
! |begin-final| is the link cell width in reduced space (cwx,cwy,cwz)

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! jdnode - destination (send to), knode - source (receive from)

! Direction -x

  If      (mdir == -1) Then

     shovex=sidex
     shovey=0.0_wp
     shovez=0.0_wp

     begin=-0.5_wp*sidex
     final=begin+cwx

     jdnode=map(1)
     kdnode=map(2)

! Direction +x

  Else If (mdir ==  1) Then

     shovex=-sidex
     shovey=0.0_wp
     shovez=0.0_wp

     final=0.5_wp*sidex
     begin=final-cwx

     jdnode=map(2)
     kdnode=map(1)

! Direction -y

  Else If (mdir == -2) Then

     shovex=0.0_wp
     shovey=sidey
     shovez=0.0_wp

     begin=-0.5_wp*sidey
     final=begin+cwy

     jdnode=map(3)
     kdnode=map(4)

! Direction +y

  Else If (mdir ==  2) Then

     shovex=0.0_wp
     shovey=-sidey
     shovez=0.0_wp

     final=0.5_wp*sidey
     begin=final-cwy

     jdnode=map(4)
     kdnode=map(3)

! Direction -z

  Else If (mdir == -3) Then

     shovex=0.0_wp
     shovey=0.0_wp
     shovez=sidez

     begin=-0.5_wp*sidez
     final=begin+cwz

     jdnode=map(5)
     kdnode=map(6)

! Direction +z

  Else If (mdir ==  3) Then

     shovex=0.0_wp
     shovey=0.0_wp
     shovez=-sidez

     final=0.5_wp*sidez
     begin=final-cwz

     jdnode=map(6)
     kdnode=map(5)

  Else

     Call error(47)

  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! Find whether a particle that belongs to this domain (idnode) falls
! into the halo of the neighbouring domain in the direction of mdir,
! i.e. the particle is within the 'minus halo' of this domain and its
! local density (rho) and index (ltg) exported into the neighbouring
! domain in the direction of mdir.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,mlast

     If      (mdir == -1 .or. mdir == 1) Then

        xyz=xxx(i)

     Else If (mdir == -2 .or. mdir == 2) Then

        xyz=yyy(i)

     Else If (mdir == -3 .or. mdir == 3) Then

        xyz=zzz(i)

     End If

! Is this particle from this domain in the halo of the neighbouring
! domain in the direction of mdir?

     If (xyz >= begin .and. xyz < final) Then

! Is it safe to proceed?

        If ((imove+5) > iblock) Then

           safe=.false.

        Else

! pack positions

           buffer(imove+1)=xxx(i)+shovex
           buffer(imove+2)=yyy(i)+shovey
           buffer(imove+3)=zzz(i)+shovez

! pack density and ltg number

           buffer(imove+4)=rho(i)
           buffer(imove+5)=Real(iwrk(i),wp)

        End If

        imove=imove+5

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

  safe=((mlast+jmove/5) <= mxatms)
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

  Do i=1,jmove/5
     mlast=mlast+1

! unpack positions

     xxx(mlast)=buffer(j+1)
     yyy(mlast)=buffer(j+2)
     zzz(mlast)=buffer(j+3)

! unpack density and ltg number

     rho(mlast)=buffer(j+4)
     iwrk(mlast)=Nint(buffer(j+5))

     j=j+5
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_export deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine metal_ld_export
