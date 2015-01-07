Subroutine statistics_connect_spread(mdir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to spread atomic and topological data of particles
! leaving this domain
!
! NOTE: When executing on one node we need not get here at all!
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,  Only : nrite,mxstak,mxatdm,mxbfss
  Use domains_module
  Use config_module, Only : ixyz

  Use msd_module,    Only : l_msd

  Use statistics_module


  Implicit None

  Integer, Intent( In    ) :: mdir

  Logical           :: safe,stay,move
  Integer           :: fail,iblock,jdnode,kdnode,   &
                       imove,jmove,kmove,keep,send, &
                       i,j,l,jj,kk,jxyz,kxyz,lxyz,  &
                       ix,iy,iz,kx,ky,kz,newatm

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  Allocate (buffer(1:mxbfss), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'statistics_connect_spread allocation failure, node: ', idnode
     Call error(0)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=mxbfss/2

! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
! kxyz - corrected halo reduction factor particles haloing both +&- sides
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
     Call error(160)
  End If

! Initialise counters for length of sending and receiving buffers
! buffer(1) and buffer(iblock+1) contain the actual number of
! particles to get transfered, imove and jmove are the lengths of
! the buffers

  imove=1
  jmove=1

! Initialise how many particles are to be kept and sent

  keep=0
  send=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL NON-RESIDENT PARTICLES ON THIS NODE

  Do i=1,natms0

! particle designated directions

     ix=Mod(ixyz(i),10)           ! [0,1,2,3]
     iy=Mod(ixyz(i)-ix,100)       ! [0,10,20,30]
     iz=Mod(ixyz(i)-(ix+iy),1000) ! [0,100,200,300]

! Filter the move index for the selected direction

     j=ix*kx+iy*ky+iz*kz

! If the particle is scheduled to be copied in the selected
! direction then indicate it (move)

     If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then
        move=.true.

! Use the corrected halo reduction factor when the particle is halo to both +&- sides

        lxyz=ixyz(i)-Merge(jxyz,kxyz,j == jxyz)

! reduce particle move index (ixyz) and decide on keeping it (stay)

        ixyz(i)=ixyz(i)-jxyz
     End If
     stay = (ixyz(i) /= 0)

     If (stay) Then ! keep it
        keep=keep+1

! retain config indexing and move indexing arrays

        ltg0(keep)=ltg0(i)
        ixyz(keep)=ixyz(i)

! retain initial positions

        xin0(keep)=xin0(i)
        yin0(keep)=yin0(i)
        zin0(keep)=zin0(i)

! retain final displacements

        xto0(keep)=xto0(i)
        yto0(keep)=yto0(i)
        zto0(keep)=zto0(i)

        If (l_msd) Then
           jj=2*i
           j =2*keep
           stpvl00(j-1)=stpvl00(jj-1)
           stpvl00(j  )=stpvl00(jj  )
           stpval0(j-1)=stpval0(jj-1)
           stpval0(j  )=stpval0(jj  )
           zumval0(j-1)=zumval0(jj-1)
           zumval0(j  )=zumval0(jj  )
           ravval0(j-1)=ravval0(jj-1)
           ravval0(j  )=ravval0(jj  )
           ssqval0(j-1)=ssqval0(jj-1)
           ssqval0(j  )=ssqval0(jj  )
           sumval0(j-1)=sumval0(jj-1)
           sumval0(j  )=sumval0(jj  )
           Do kk=1,mxstak
              stkval0(kk,j-1)=stkval0(kk,jj-1)
              stkval0(kk,j  )=stkval0(kk,jj  )
           End Do
        End If
     End If

     If (move) Then ! copy it
        send=send+1
        If (imove+8 <= iblock) Then ! If safe to proceed

! pack config indexing and move indexing arrays

           buffer(imove+1)=Real(ltg0(i),wp)
           buffer(imove+2)=Real(lxyz,wp)

! pack initial positions

           buffer(imove+3)=xin0(i)
           buffer(imove+4)=yin0(i)
           buffer(imove+5)=zin0(i)

! pack final displacements

           buffer(imove+6)=xto0(i)
           buffer(imove+7)=yto0(i)
           buffer(imove+8)=zto0(i)
        Else
           safe=.false.
        End If
        imove=imove+8

! pack MSD arrays

        If (l_msd) Then
           If (imove+2*(6+mxstak) <= iblock) Then
              jj=2*i
              buffer(imove+ 1)=stpvl0(jj-1)
              buffer(imove+ 2)=stpvl0(jj  )
              buffer(imove+ 3)=stpval(jj-1)
              buffer(imove+ 4)=stpval(jj  )
              buffer(imove+ 5)=zumval(jj-1)
              buffer(imove+ 6)=zumval(jj  )
              buffer(imove+ 7)=ravval(jj-1)
              buffer(imove+ 8)=ravval(jj  )
              buffer(imove+ 9)=ssqval(jj-1)
              buffer(imove+10)=ssqval(jj  )
              buffer(imove+11)=sumval(jj-1)
              buffer(imove+12)=sumval(jj  )
              Do kk=1,mxstak
                 l=2*kk
                 buffer(imove+12+l-1)=stkval(kk,jj-1)
                 buffer(imove+12+l  )=stkval(kk,jj  )
              End Do
           Else
              safe=.false.
           End If
           imove=imove+2*(6+mxstak)
        End If
     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  Call gcheck(safe)
  If (.not.safe) Call error(163)

! record of number of atoms for transfer

  buffer(1)=Real(send,wp)

! exchange information on buffer sizes

  Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,Spread_tag,dlp_comm_world,request,ierr)
  Call MPI_SEND(imove,1,MPI_INTEGER,jdnode,Spread_tag,dlp_comm_world,ierr)
  Call MPI_WAIT(request,status,ierr)

! exchange buffers between nodes (this is a MUST)

  Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,Spread_tag,dlp_comm_world,request,ierr)
  Call MPI_SEND(buffer(1),imove,wp_mpi,jdnode,Spread_tag,dlp_comm_world,ierr)
  Call MPI_WAIT(request,status,ierr)

! check arrays can cope with incoming atom numbers

  kmove=iblock+1
  jmove=Nint(buffer(kmove))

  natms0=keep+jmove

! Check for array bound overflow (can arrays cope with incoming data)

  safe=(natms0 <= mxatdm)
  Call gcheck(safe)
  If (.not.safe) Call error(164)

! load transferred data

  Do i=1,jmove
     newatm=i+keep

! unpack config indexing, site and move indexing arrays

     ltg0(newatm)=Nint(buffer(kmove+1))
     ixyz(newatm)=Nint(buffer(kmove+2))

! unpack initial positions arrays

     xin0(newatm)=buffer(kmove+3)
     yin0(newatm)=buffer(kmove+4)
     zin0(newatm)=buffer(kmove+5)

! unpack initial positions arrays

     xto0(newatm)=buffer(kmove+6)
     yto0(newatm)=buffer(kmove+7)
     zto0(newatm)=buffer(kmove+8)

     kmove=kmove+8

! unpack MSD arrays

     If (l_msd) Then
        jj=2*newatm
        stpvl00(jj-1)=buffer(kmove+1 )
        stpvl00(jj  )=buffer(kmove+2 )
        stpval0(jj-1)=buffer(kmove+3 )
        stpval0(jj  )=buffer(kmove+4 )
        zumval0(jj-1)=buffer(kmove+5 )
        zumval0(jj  )=buffer(kmove+6 )
        ravval0(jj-1)=buffer(kmove+7 )
        ravval0(jj  )=buffer(kmove+8 )
        ssqval0(jj-1)=buffer(kmove+9 )
        ssqval0(jj  )=buffer(kmove+10)
        sumval0(jj-1)=buffer(kmove+11)
        sumval0(jj  )=buffer(kmove+12)
        Do kk=1,mxstak
           l=2*kk
           stkval0(kk,jj-1)=buffer(kmove+12+l-1)
           stkval0(kk,jj  )=buffer(kmove+12+l  )
        End Do

        kmove=kmove+2*(6+mxstak)
     End If
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'statistics_connect_spread deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine statistics_connect_spread
