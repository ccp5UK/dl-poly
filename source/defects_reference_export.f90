Subroutine defects_reference_export &
           (mdir,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic REFERENCE data in domain boundary
! regions for halo formation
!
! all particle coordinates are in reduced space with origin localised
! onto this node (idnode)
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module

  Implicit None

  Integer,              Intent( In    ) :: mdir
  Real( Kind = wp ),    Intent( In    ) :: sidex,sidey,sidez,cwx,cwy,cwz
  Integer,              Intent( InOut ) :: nlrefs
  Character( Len = 8 ), Intent( InOut ) :: namr(1:mxatms)
  Integer,              Intent( InOut ) :: indr(1:mxatms)
  Real( Kind = wp ),    Intent( InOut ) :: xr(1:mxatms),yr(1:mxatms),zr(1:mxatms)


  Logical           :: safe
  Integer           :: fail,i,j,iblock,jdnode,kdnode,imove,jmove,itmp
  Real( Kind = wp ) :: shovex,shovey,shovez,begin,final,xyz

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer


  fail=0
  Allocate (buffer(1:mxbuff), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_export allocation failure, node: ', idnode
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

     Call error(557)

  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! Find whether a particle that belongs to this domain (idnode) falls
! into the halo of the nighbouring domain in the direction of mdir, i.e.
! the particle is is within the 'minus halo' of this domain and has
! to be exported into the neighbouring domain in the direction of mdir.
! If a particle is to be exported across domains then all config
! properties of the particle (i.e. the particle itself) also have to
! be exported across to the receiving domain.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,nlrefs

     If      (mdir == -1 .or. mdir == 1) Then

        xyz=xr(i)

     Else If (mdir == -2 .or. mdir == 2) Then

        xyz=yr(i)

     Else If (mdir == -3 .or. mdir == 3) Then

        xyz=zr(i)

     End If

! Is this particle from this domain in the halo of the neighbouring
! domain in the direction of mdir?

     If (xyz >= begin .and. xyz < final) Then

! Is it safe to proceed?

        If ((imove+12) > iblock) Then

           safe=.false.

        Else

! pack positions

           buffer(imove+ 1)=xr(i)+shovex
           buffer(imove+ 2)=yr(i)+shovey
           buffer(imove+ 3)=zr(i)+shovez

! pack config indexing and site name arrays

           buffer(imove+ 4)=Real(indr(i),wp)
           buffer(imove+ 5)=Real(Ichar(namr(i)(1:1)),wp)
           buffer(imove+ 6)=Real(Ichar(namr(i)(2:2)),wp)
           buffer(imove+ 7)=Real(Ichar(namr(i)(3:3)),wp)
           buffer(imove+ 8)=Real(Ichar(namr(i)(4:4)),wp)
           buffer(imove+ 9)=Real(Ichar(namr(i)(5:5)),wp)
           buffer(imove+10)=Real(Ichar(namr(i)(6:6)),wp)
           buffer(imove+11)=Real(Ichar(namr(i)(7:7)),wp)
           buffer(imove+12)=Real(Ichar(namr(i)(8:8)),wp)

        End If

        imove=imove+12

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

  safe=((nlrefs+jmove/12) <= mxatms)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     itmp=nlrefs+jmove/12
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

  If (mxnode > 1) Then
     j=iblock
  Else
     j=0
  End If

  Do i=1,jmove/12
     nlrefs=nlrefs+1

! unpack positions

     xr(nlrefs)=buffer(j+1)
     yr(nlrefs)=buffer(j+2)
     zr(nlrefs)=buffer(j+3)

! unpack config indexing and site arrays

     indr(nlrefs)=Nint(buffer(j+4))
     namr(nlrefs)(1:1)=Char(Nint(buffer(j+ 5)))
     namr(nlrefs)(2:2)=Char(Nint(buffer(j+ 6)))
     namr(nlrefs)(3:3)=Char(Nint(buffer(j+ 7)))
     namr(nlrefs)(4:4)=Char(Nint(buffer(j+ 8)))
     namr(nlrefs)(5:5)=Char(Nint(buffer(j+ 9)))
     namr(nlrefs)(6:6)=Char(Nint(buffer(j+10)))
     namr(nlrefs)(7:7)=Char(Nint(buffer(j+11)))
     namr(nlrefs)(8:8)=Char(Nint(buffer(j+12)))

     j=j+12
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_export deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine defects_reference_export
