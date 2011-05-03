Subroutine deport_atomic_data(mdir,lbook,sidex,sidey,sidez,cwx,cwy,cwz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to deport atomic data in domain boundary regions
!
! all particle coordinates are in reduced space with origin localised
! onto this node (idnode)
!
! NOTE: When executing on one node we need not get here at all!
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov april 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use domains_module

  Use config_module

  Use core_shell_module,   Only : ntshl,listshl,legshl

  Use constraints_module,  Only : ntcons,listcon,legcon
  Use pmf_module,          Only : ntpmf,listpmf,legpmf

  Use rigid_bodies_module, Only : ntrgd,listrgd,legrgd, &
                                  q0,q1,q2,q3,          &
                                  rgdvxx,rgdvyy,rgdvzz, &
                                  rgdoxx,rgdoyy,rgdozz

  Use tethers_module,      Only : ntteth,listtet,legtet

  Use bonds_module,        Only : ntbond,listbnd,legbnd
  Use angles_module,       Only : ntangl,listang,legang
  Use dihedrals_module,    Only : ntdihd,listdih,legdih
  Use inversions_module,   Only : ntinv,listinv,leginv

  Use statistics_module

  Use langevin_module,     Only : l_lan,fxl,fyl,fzl
  Use minimise_module,     Only : l_x,oxx,oyy,ozz
  Use ewald_module,        Only : l_fce,fcx,fcy,fcz
  Use msd_module

  Implicit None

  Logical,           Intent( In    ) :: lbook
  Integer,           Intent( In    ) :: mdir
  Real( Kind = wp ), Intent( In    ) :: sidex,sidey,sidez,cwx,cwy,cwz

  Logical           :: safe,safe1,check
  Integer           :: fail(1:3),iblock,jdnode,kdnode,         &
                       imove,jmove,kmove,keep,                 &
                       i,j,k,l,jj,kk,                          &
                       newatm,iatm,jatm,katm,latm,             &
                       jshels,kshels,jconst,kconst,jpmf,kpmf,  &
                       jrigid,krigid,jteths,kteths,            &
                       jbonds,kbonds,jangle,kangle,            &
                       jdihed,kdihed,jinver,kinver
  Real( Kind = wp ) :: shovex,shovey,shovez,begin,final,xyz

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Integer,           Dimension( : ), Allocatable :: i1pmf,i2pmf
  Integer,           Dimension( : ), Allocatable :: lrgd

  fail=0
  Allocate (buffer(1:mxbuff),                      Stat=fail(1))
  Allocate (i1pmf(1:mxtpmf(1)),i2pmf(1:mxtpmf(2)), Stat=fail(2))
  Allocate (lrgd(-1:Max(mxlrgd,mxrgd)),            Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'deport_atomic_data allocation failure, node: ', idnode
     Call error(0)
  End If


! Set buffer limit (half for outgoing data - half for incoming)

  iblock=mxbuff/2

! DIRECTION SETTINGS INITIALISATION

! define the relative displacement between the coordinate systems'
! origins of neighbouring domains with respect to this domain
! (idnode) in the direction of mdir (shovex,shovey,shovez)
! in order to push all particles of this domain (idnode) in the
! direction of mdir

! define halo limits in the direction of mdir (begin,final)
! |begin-final| is the link cell width in reduced space (cwx,cwy,cwz)

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! jdnode - destination (send to), knode - source (receive from)

! Direction -x

  If      (mdir == -1) Then

     shovex=sidex
     shovey=0.0_wp
     shovez=0.0_wp

     final=-0.5_wp*sidex
     begin=final-cwx

     jdnode=map(1)
     kdnode=map(2)

! Direction +x

  Else If (mdir ==  1) Then

     shovex=-sidex
     shovey=0.0_wp
     shovez=0.0_wp

     begin=0.5_wp*sidex
     final=begin+cwx

     jdnode=map(2)
     kdnode=map(1)

! Direction -y

  Else If (mdir == -2) Then

     shovex=0.0_wp
     shovey=sidey
     shovez=0.0_wp

     final=-0.5_wp*sidey
     begin=final-cwy

     jdnode=map(3)
     kdnode=map(4)

! Direction +y

  Else If (mdir ==  2) Then

     shovex=0.0_wp
     shovey=-sidey
     shovez=0.0_wp

     begin=0.5_wp*sidey
     final=begin+cwy

     jdnode=map(4)
     kdnode=map(3)

! Direction -z

  Else If (mdir == -3) Then

     shovex=0.0_wp
     shovey=0.0_wp
     shovez=sidez

     final=-0.5_wp*sidez
     begin=final-cwz

     jdnode=map(5)
     kdnode=map(6)

! Direction +z

  Else If (mdir ==  3) Then

     shovex=0.0_wp
     shovey=0.0_wp
     shovez=-sidez

     begin=0.5_wp*sidez
     final=begin+cwz

     jdnode=map(6)
     kdnode=map(5)

  Else

     Call error(42)

  End If


! Initialise counters for length of sending and receiving buffers
! buffer(1) and buffer(iblock+1) contain the actual number of
! particles to get transfered, imove and jmove are the lengths of
! the buffers

  imove=1
  jmove=1

! Initialise how many particles are to be kept

  keep=0

! Initialise array overflow flags

  safe=.true.
  safe1=.true.

! Initialise RB data set for preventing duplication

  lrgd(1:ntrgd)=0


! Find whether a particle that belongs to this domain (idnode) has
! moved into the halo of the domain in the direction of mdir, i.e.
! the particle is not within the boundaries of this domain and has
! to be exported into the neighbouring domain in the direction of mdir.
! If a particle is to be moved across domains then all config and intra
! properties of the particle (i.e. the particle itself) also have to
! be moved across to the receiving domain.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,natms

     If      (mdir == -1 .or. mdir == 1) Then

        xyz=xxx(i)

     Else If (mdir == -2 .or. mdir == 2) Then

        xyz=yyy(i)

     Else If (mdir == -3 .or. mdir == 3) Then

        xyz=zzz(i)

     End If

! Has this particle moved from this domain into the halo of this domain?

     If (xyz >= begin .and. xyz < final) Then

! Is it safe to proceed?

        If (imove+17 <= iblock) Then

! pack positions

           buffer(imove+1)=xxx(i)+shovex
           buffer(imove+2)=yyy(i)+shovey
           buffer(imove+3)=zzz(i)+shovez

! pack velocities

           buffer(imove+4)=vxx(i)
           buffer(imove+5)=vyy(i)
           buffer(imove+6)=vzz(i)

! pack forces

           buffer(imove+7)=fxx(i)
           buffer(imove+8)=fyy(i)
           buffer(imove+9)=fzz(i)

! pack config bookkeeping arrays

           buffer(imove+10)=Real(ltg(i),wp)
           buffer(imove+11)=Real(lsite(i),wp)

! pack initial positions

           buffer(imove+12)=xin(i)
           buffer(imove+13)=yin(i)
           buffer(imove+14)=zin(i)

! pack final displacements

           buffer(imove+15)=xto(i)
           buffer(imove+16)=yto(i)
           buffer(imove+17)=zto(i)

           imove=imove+17

        Else

           safe=.false.

        End If

! pack Langevin force arrays

        If (l_lan) Then
           If (imove+3 <= iblock) Then
              buffer(imove+1)=fxl(i)
              buffer(imove+2)=fyl(i)
              buffer(imove+3)=fzl(i)

              imove=imove+3
           Else
              safe=.false.
           End If
        End If

! pack minimisation arrays

        If (l_x) Then
           If (imove+3 <= iblock) Then
              buffer(imove+1)=oxx(i)
              buffer(imove+2)=oyy(i)
              buffer(imove+3)=ozz(i)

              imove=imove+3
           Else
              safe=.false.
           End If
        End If

! pack k-space SPME force arrays

        If (.not.l_fce) Then
           If (imove+3 <= iblock) Then
              buffer(imove+1)=fcx(i)
              buffer(imove+2)=fcy(i)
              buffer(imove+3)=fcz(i)

              imove=imove+3
           Else
              safe=.false.
           End If
        End If

! pack MSD arrays

        If (l_msd) Then
           If (imove+2*(6+mxstak) <= iblock) Then
              jj=27+2*i
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

              imove=imove+2*(6+mxstak)
           Else
              safe=.false.
           End If
        End If

! If intra-molecular entities exist in the system

        If (lbook .and. safe) Then

! pack the exclusion list

           kk=lexatm(0,i)
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=Real(kk,wp)
           Else
              safe=.false.
           End If
           If (imove+kk <= iblock .and. safe) Then
              Do k=1,kk
                 imove=imove+1
                 buffer(imove)=Real(lexatm(k,i),wp)
              End Do
           Else
              safe=.false.
           End If

! the order of packing up intra bookkeeping arays must be the same as
! their scanning in build_book_intra

! pack core-shell details

           jj=1
           Do While (legshl(jj,i) > 0 .and. safe)
              If (imove+3 <= iblock) Then
                 kk=legshl(jj,i)

                 Do k=0,2
                    imove=imove+1
                    buffer(imove)=Real(listshl(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack contraint details

           jj=1
           Do While (legcon(jj,i) > 0 .and. safe)
              If (imove+3 <= iblock) Then
                 kk=legcon(jj,i)

                 Do k=0,2
                    imove=imove+1
                    buffer(imove)=Real(listcon(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack PMF details

           jj=1
           Do While (legpmf(jj,i) > 0 .and. safe)
              If (imove+mxtpmf(1)+mxtpmf(2)+2 <= iblock) Then
                 kk=legpmf(jj,i)

                 Do k=0,mxtpmf(1)
                    imove=imove+1
                    buffer(imove)=Real(listpmf(k,1,kk),wp)
                 End Do

                 Do k=0,mxtpmf(2)
                    imove=imove+1
                    buffer(imove)=Real(listpmf(k,2,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack RB details

!           jj=1
!           Do While (legrgd(jj,i) > 0 .and. safe)
!              l=12
!              If (imove+mxlrgd+l <= iblock) Then
!                 kk=legrgd(jj,i)
!
!                 Do k=-1,listrgd(-1,kk)
!                    imove=imove+1
!                    buffer(imove)=Real(listrgd(k,kk),wp)
!                 End Do
!
!                 buffer(imove+1)=q0(kk)
!                 buffer(imove+2)=q1(kk)
!                 buffer(imove+3)=q2(kk)
!                 buffer(imove+4)=q3(kk)
!                 imove=imove+4
!
!                 buffer(imove+1)=rgdvxx(kk)
!                 buffer(imove+2)=rgdvyy(kk)
!                 buffer(imove+3)=rgdvzz(kk)
!                 imove=imove+3
!
!                 buffer(imove+1)=rgdoxx(kk)
!                 buffer(imove+2)=rgdoyy(kk)
!                 buffer(imove+3)=rgdozz(kk)
!                 imove=imove+3
!
!                 jj=jj+1
!              Else
!                 safe=.false.
!              End If
!           End Do
!           If (imove+1 <= iblock) Then
!              imove=imove+1
!              buffer(imove)=0.0_wp
!           Else
!              safe=.false.
!           End If

           jj=1
           Do While (legrgd(jj,i) > 0 .and. safe)
              l=12
              If (imove+mxlrgd+l <= iblock) Then
                 kk=legrgd(jj,i)

                 Do k=-1,1
                    imove=imove+1
                    buffer(imove)=Real(listrgd(k,kk),wp)
                 End Do

! Bypass if the data has already been sent

                 If (lrgd(kk) == 0) Then
                    Do k=2,listrgd(-1,kk)
                       imove=imove+1
                       buffer(imove)=Real(listrgd(k,kk),wp)
                    End Do

                    buffer(imove+1)=q0(kk)
                    buffer(imove+2)=q1(kk)
                    buffer(imove+3)=q2(kk)
                    buffer(imove+4)=q3(kk)
                    imove=imove+4

                    buffer(imove+1)=rgdvxx(kk)
                    buffer(imove+2)=rgdvyy(kk)
                    buffer(imove+3)=rgdvzz(kk)
                    imove=imove+3

                    buffer(imove+1)=rgdoxx(kk)
                    buffer(imove+2)=rgdoyy(kk)
                    buffer(imove+3)=rgdozz(kk)
                    imove=imove+3

                    lrgd(kk)=1
                 End If

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack tether details

           jj=1
           Do While (legtet(jj,i) > 0 .and. safe)
              If (imove+2 <= iblock) Then
                 kk=legtet(jj,i)

                 Do k=0,1
                    imove=imove+1
                    buffer(imove)=Real(listtet(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack bond details

           jj=1
           Do While (legbnd(jj,i) > 0 .and. safe)
              If (imove+3 <= iblock) Then
                 kk=legbnd(jj,i)

                 Do k=0,2
                    imove=imove+1
                    buffer(imove)=Real(listbnd(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack valence angle details

           jj=1
           Do While (legang(jj,i) > 0 .and. safe)
              If (imove+4 <= iblock) Then
                 kk=legang(jj,i)

                 Do k=0,3
                    imove=imove+1
                    buffer(imove)=Real(listang(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack dihedral angle details

           jj=1
           Do While (legdih(jj,i) > 0 .and. safe)
              If (imove+5 <= iblock) Then
                 kk=legdih(jj,i)

                 Do k=0,4
                    imove=imove+1
                    buffer(imove)=Real(listdih(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

! pack inversion angle details

           jj=1
           Do While (leginv(jj,i) > 0 .and. safe)
              If (imove+5 <= iblock) Then
                 kk=leginv(jj,i)

                 Do k=0,4
                    imove=imove+1
                    buffer(imove)=Real(listinv(k,kk),wp)
                 End Do

                 jj=jj+1
              Else
                 safe=.false.
              End If
           End Do
           If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=0.0_wp
           Else
              safe=.false.
           End If

        End If

     Else

        keep=keep+1

        xxx(keep)=xxx(i)
        yyy(keep)=yyy(i)
        zzz(keep)=zzz(i)

        vxx(keep)=vxx(i)
        vyy(keep)=vyy(i)
        vzz(keep)=vzz(i)

        fxx(keep)=fxx(i)
        fyy(keep)=fyy(i)
        fzz(keep)=fzz(i)

        ltg(keep)=ltg(i)
        lsite(keep)=lsite(i)

        xin(keep)=xin(i)
        yin(keep)=yin(i)
        zin(keep)=zin(i)

        xto(keep)=xto(i)
        yto(keep)=yto(i)
        zto(keep)=zto(i)

        If (l_lan) Then
           fxl(keep)=fxl(i)
           fyl(keep)=fyl(i)
           fzl(keep)=fzl(i)
        End If

        If (l_x) Then
           oxx(keep)=oxx(i)
           oyy(keep)=oyy(i)
           ozz(keep)=ozz(i)
        End If

        If (.not.l_fce) Then
           fcx(keep)=fcx(i)
           fcy(keep)=fcy(i)
           fcz(keep)=fcz(i)
        End If

        If (l_msd) Then
           jj=27+2*i
           j =27+2*keep
           stpvl0(j-1)=stpvl0(jj-1)
           stpvl0(j  )=stpvl0(jj  )
           stpval(j-1)=stpval(jj-1)
           stpval(j  )=stpval(jj  )
           zumval(j-1)=zumval(jj-1)
           zumval(j  )=zumval(jj  )
           ravval(j-1)=ravval(jj-1)
           ravval(j  )=ravval(jj  )
           ssqval(j-1)=ssqval(jj-1)
           ssqval(j  )=ssqval(jj  )
           sumval(j-1)=sumval(jj-1)
           sumval(j  )=sumval(jj  )
           Do kk=1,mxstak
              stkval(kk,j-1)=stkval(kk,jj-1)
              stkval(kk,j  )=stkval(kk,jj  )
           End Do
        End If

        If (lbook) Then
           lexatm(:,keep)=lexatm(:,i)

           legshl(:,keep)=legshl(:,i)

           legcon(:,keep)=legcon(:,i)
           legpmf(:,keep)=legpmf(:,i)

           legrgd(:,keep)=legrgd(:,i)

           legtet(:,keep)=legtet(:,i)

           legbnd(:,keep)=legbnd(:,i)
           legang(:,keep)=legang(:,i)
           legdih(:,keep)=legdih(:,i)
           leginv(:,keep)=leginv(:,i)
        End If

     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  Call gcheck(safe)
  If (.not.safe) Call error(43)

! record of number of atoms for transfer

  buffer(1)=Real(natms-keep,wp)

! exchange information on buffer sizes

  Call MPI_IRECV(jmove,1,MPI_INTEGER,kdnode,Deport_tag,dlp_comm_world,request,ierr)
  Call MPI_SEND(imove,1,MPI_INTEGER,jdnode,Deport_tag,dlp_comm_world,ierr)
  Call MPI_WAIT(request,status,ierr)

! exchange buffers between nodes (this is a MUST)

  Call MPI_IRECV(buffer(iblock+1),jmove,wp_mpi,kdnode,Deport_tag,dlp_comm_world,request,ierr)
  Call MPI_SEND(buffer(1),imove,wp_mpi,jdnode,Deport_tag,dlp_comm_world,ierr)
  Call MPI_WAIT(request,status,ierr)

! check arrays can cope with incoming atom numbers

  kmove=iblock+1
  jmove=Nint(buffer(kmove))

  natms=keep+jmove

! Check for array bound overflow (can arrays cope with incoming data)

  safe=(natms <= mxatms)
  Call gcheck(safe)
  If (.not.safe) Call error(44)

! load transferred data

  Do i=1,jmove
     newatm=i+keep

! unpack positions

     xxx(newatm)=buffer(kmove+1)
     yyy(newatm)=buffer(kmove+2)
     zzz(newatm)=buffer(kmove+3)

! unpack velocities

     vxx(newatm)=buffer(kmove+4)
     vyy(newatm)=buffer(kmove+5)
     vzz(newatm)=buffer(kmove+6)

! unpack forces

     fxx(newatm)=buffer(kmove+7)
     fyy(newatm)=buffer(kmove+8)
     fzz(newatm)=buffer(kmove+9)

! unpack config bookkeeping arrays

     ltg(newatm)=Nint(buffer(kmove+10))
     lsite(newatm)=Nint(buffer(kmove+11))

! unpack initial positions arrays

     xin(newatm)=buffer(kmove+12)
     yin(newatm)=buffer(kmove+13)
     zin(newatm)=buffer(kmove+14)

! unpack initial positions arrays

     xto(newatm)=buffer(kmove+15)
     yto(newatm)=buffer(kmove+16)
     zto(newatm)=buffer(kmove+17)

     kmove=kmove+17

! unpack Langevin force arrays

     If (l_lan) Then
        fxl(newatm)=buffer(kmove+1)
        fyl(newatm)=buffer(kmove+2)
        fzl(newatm)=buffer(kmove+3)

        kmove=kmove+3
     End If

! unpack minimisation arrays

     If (l_x) Then
        oxx(newatm)=buffer(kmove+1)
        oyy(newatm)=buffer(kmove+2)
        ozz(newatm)=buffer(kmove+3)

        kmove=kmove+3
     End If

! unpack k-space SPME force arrays

     If (.not.l_fce) Then
        fcx(newatm)=buffer(kmove+1)
        fcy(newatm)=buffer(kmove+2)
        fcz(newatm)=buffer(kmove+3)

        kmove=kmove+3
     End If

! unpack MSD arrays

     If (l_msd) Then
        jj=27+2*newatm
        stpvl0(jj-1)=buffer(kmove+1 )
        stpvl0(jj  )=buffer(kmove+2 )
        stpval(jj-1)=buffer(kmove+3 )
        stpval(jj  )=buffer(kmove+4 )
        zumval(jj-1)=buffer(kmove+5 )
        zumval(jj  )=buffer(kmove+6 )
        ravval(jj-1)=buffer(kmove+7 )
        ravval(jj  )=buffer(kmove+8 )
        ssqval(jj-1)=buffer(kmove+9 )
        ssqval(jj  )=buffer(kmove+10)
        sumval(jj-1)=buffer(kmove+11)
        sumval(jj  )=buffer(kmove+12)
        Do kk=1,mxstak
           l=2*kk
           stkval(kk,jj-1)=buffer(kmove+12+l-1)
           stkval(kk,jj  )=buffer(kmove+12+l  )
        End Do

        kmove=kmove+2*(6+mxstak)
     End If

     If (lbook) Then

! unpack the exclusion list

        kmove=kmove+1
        kk=Nint(buffer(kmove))
        lexatm(0,newatm)=kk
        Do k=1,kk
           kmove=kmove+1
           lexatm(k,newatm)=Nint(buffer(kmove))
        End Do
        lexatm(kk+1:mxexcl,newatm)=0

! the order of unpacking intra bookkeeping arays must be the same as
! the order of their scanning in build_book_intra in order to rebuild
! correctly the new list arrayes and create new legend arrays

! set initial intra counters

        jshels=ntshl

        jconst=ntcons
        jpmf  =ntpmf

        jrigid=ntrgd

        jteths=ntteth

        jbonds=ntbond
        jangle=ntangl
        jdihed=ntdihd
        jinver=ntinv

! unpack core-shell details

        legshl(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           kmove=kmove+3

! check if core-shell unit already specified

           kshels=0
           check=.true.
           Do While (check .and. kshels < Min(jshels,mxshl))
              kshels=kshels+1
              check=.not.( jj   == listshl(0,kshels) .and. & ! core-shell units don't intersect
                           iatm == listshl(1,kshels))        ! .and. jatm == listshl(2,kshels) )
           End Do

! add new core-shell unit

           If (check) Then
              jshels=jshels+1

              If (jshels <= mxshl) Then
                 listshl(0,jshels)=jj
                 listshl(1,jshels)=iatm
                 listshl(2,jshels)=jatm

                 Call tag_legend(safe1,newatm,jshels,legshl,mxfshl)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kshels,legshl,mxfshl)
           End If
        End Do
        kmove=kmove+1

! unpack bond contraint details

        legcon(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           kmove=kmove+3

! check if constraint already specified

           kconst=0
           check=.true.
           Do While (check .and. kconst < Min(jconst,mxcons))
              kconst=kconst+1
              check=.not.( jj   == listcon(0,kconst) .and. &
                           iatm == listcon(1,kconst) .and. &
                           jatm == listcon(2,kconst) )
           End Do

! insert new constraint unit

           If (check) Then
              jconst=jconst+1

              If (jconst <= mxcons) Then
                 listcon(0,jconst)=jj
                 listcon(1,jconst)=iatm
                 listcon(2,jconst)=jatm

                 Call tag_legend(safe1,newatm,jconst,legcon,mxfcon)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kconst,legcon,mxfcon)
           End If
        End Do
        kmove=kmove+1

! unpack PMF details

        legpmf(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1)) ! PMF global indentifier
           kmove=kmove+1
           Do k=1,mxtpmf(1)
              kmove=kmove+1
              i1pmf(k)=Nint(buffer(kmove))
           End Do

           kmove=kmove+1            ! omit PMF units presence indentifier
           Do k=1,mxtpmf(2)         ! and deal with it in pmf_units_set
              kmove=kmove+1
              i2pmf(k)=Nint(buffer(kmove))
           End Do

! check if PMF already specified

           kpmf=0
           check=.true.
           Do While (check .and. kpmf < Min(jpmf,mxpmf))
              kpmf=kpmf+1
              check=.not.(jj == listpmf(0,1,kpmf))
           End Do

! insert new PMF constraint

           If (check) Then
              jpmf=jpmf+1

              If (jpmf <= mxpmf) Then
                 listpmf(0,1,jpmf)=jj
                 Do k=1,mxtpmf(1)
                    listpmf(k,1,jpmf)=i1pmf(k)
                 End Do

! PMF units presence indentifier holds zero temporarely
! it's dealt with in pmf_units_set

                 listpmf(0,2,jpmf)=0
                 Do k=1,mxtpmf(2)
                    listpmf(k,2,jpmf)=i2pmf(k)
                 End Do

                 Call tag_legend(safe1,newatm,jpmf,legpmf,mxfpmf)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kpmf,legpmf,mxfpmf)
           End If
        End Do
        kmove=kmove+1

! unpack RB details

!        legrgd(:,newatm) = 0
!        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
!           lrgd(-1)=Nint(buffer(kmove+1))
!           kmove=kmove+1
!           Do k=0,lrgd(-1)
!              kmove=kmove+1
!              lrgd(k)=Nint(buffer(kmove))
!           End Do
!
!! check if RB already specified
!
!           krigid=0
!           check=.true.
!           Do While (check .and. krigid < Min(jrigid,mxrgd))
!              krigid=krigid+1
!              check=( lrgd( 0) == listrgd( 0,krigid) .and. & ! Type
!                      lrgd(-1) == listrgd(-1,krigid) .and. & ! Size
!                      lrgd( 1) == listrgd( 1,krigid) )       ! Global ID of just the first
!!                                                            ! member as RBs don't intersect
!              check=.not.check
!           End Do
!
!! insert new RB unit
!
!           If (check) Then
!              jrigid=jrigid+1
!
!              If (jrigid <= mxrgd) Then
!                 Do k=-1,lrgd(-1)
!                    listrgd(k,jrigid)=lrgd(k)
!                 End Do
!
!                 Call tag_legend(safe1,newatm,jrigid,legrgd,mxfrgd)
!
!                 q0(jrigid)=buffer(kmove+1)
!                 q1(jrigid)=buffer(kmove+2)
!                 q2(jrigid)=buffer(kmove+3)
!                 q3(jrigid)=buffer(kmove+4)
!                 kmove=kmove+4
!
!                 rgdvxx(jrigid)=buffer(kmove+1)
!                 rgdvyy(jrigid)=buffer(kmove+2)
!                 rgdvzz(jrigid)=buffer(kmove+3)
!                 kmove=kmove+3
!
!                 rgdoxx(jrigid)=buffer(kmove+1)
!                 rgdoyy(jrigid)=buffer(kmove+2)
!                 rgdozz(jrigid)=buffer(kmove+3)
!                 kmove=kmove+3
!              Else
!                 safe=.false.
!              End If
!           Else
!              l=10
!              kmove=kmove+l ! Compensate for the 'l' extra unread buffers
!              Call tag_legend(safe1,newatm,krigid,legrgd,mxfrgd)
!           End If
!        End Do
!        kmove=kmove+1

        legrgd(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           lrgd(-1)=Nint(buffer(kmove+1))
           kmove=kmove+1
           Do k=0,2
              kmove=kmove+1
              lrgd(k)=Nint(buffer(kmove))
           End Do

! check if RB already specified

           krigid=0
           check=.true.
           Do While (check .and. krigid < Min(jrigid,mxrgd))
              krigid=krigid+1
              check=( lrgd( 0) == listrgd( 0,krigid) .and. & ! Type
                      lrgd(-1) == listrgd(-1,krigid) .and. & ! Size
                      lrgd( 1) == listrgd( 1,krigid) )       ! Global ID of just the first
!                                                            ! member as RBs don't intersect
              check=.not.check
           End Do

! insert new RB unit

           If (check) Then
              Do k=3,lrgd(-1)
                 kmove=kmove+1
                 lrgd(k)=Nint(buffer(kmove))
              End Do

              jrigid=jrigid+1

              If (jrigid <= mxrgd) Then
                 Do k=-1,lrgd(-1)
                    listrgd(k,jrigid)=lrgd(k)
                 End Do

                 Call tag_legend(safe1,newatm,jrigid,legrgd,mxfrgd)

                 q0(jrigid)=buffer(kmove+1)
                 q1(jrigid)=buffer(kmove+2)
                 q2(jrigid)=buffer(kmove+3)
                 q3(jrigid)=buffer(kmove+4)
                 kmove=kmove+4

                 rgdvxx(jrigid)=buffer(kmove+1)
                 rgdvyy(jrigid)=buffer(kmove+2)
                 rgdvzz(jrigid)=buffer(kmove+3)
                 kmove=kmove+3

                 rgdoxx(jrigid)=buffer(kmove+1)
                 rgdoyy(jrigid)=buffer(kmove+2)
                 rgdozz(jrigid)=buffer(kmove+3)
                 kmove=kmove+3
              Else
                 safe=.false.
              End If
           Else
              If (lrgd(2) == 0) Then  ! Unduplication: Details have already been sent -
                 kmove=kmove-1        ! back up once for the zero read at end of datapack
              Else                    ! Data already persent - jump over it as it's not needed
                 l=10                 ! Compensate for the 'l' extra unread buffers
                 kmove=kmove+l+lrgd(-1)-2
              End If
              Call tag_legend(safe1,newatm,krigid,legrgd,mxfrgd)
           End If
        End Do
        kmove=kmove+1

! unpack tether details

        legtet(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           kmove=kmove+2

! check if tether unit already specified

           kteths=0
           check=.true.
           Do While (check .and. kteths < Min(jteths,mxteth))
              kteths=kteths+1
              check=.not.( jj   == listtet(0,kteths) .and. &
                           iatm == listtet(1,kteths) )
           End Do

! add new tether unit

           If (check) Then
              jteths=jteths+1

              If (jteths <= mxteth) Then
                 listtet(0,jteths)=jj
                 listtet(1,jteths)=iatm

                 Call tag_legend(safe1,newatm,jteths,legtet,mxftet)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kteths,legtet,mxftet)
           End If
        End Do
        kmove=kmove+1

! unpack bond details

        legbnd(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           kmove=kmove+3

! check if bond already specified

           kbonds=0
           check=.true.
           Do While (check .and. kbonds < Min(jbonds,mxbond))
              kbonds=kbonds+1
              check=.not.( jj   == listbnd(0,kbonds) .and. &
                           iatm == listbnd(1,kbonds) .and. &
                           jatm == listbnd(2,kbonds) )
           End Do

! insert new bond details

           If (check) Then
              jbonds=jbonds+1

              If (jbonds <= mxbond) Then
                 listbnd(0,jbonds)=jj
                 listbnd(1,jbonds)=iatm
                 listbnd(2,jbonds)=jatm

                 Call tag_legend(safe1,newatm,jbonds,legbnd,mxfbnd)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kbonds,legbnd,mxfbnd)
           End If
        End Do
        kmove=kmove+1

! unpack valence angle details

        legang(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           katm=Nint(buffer(kmove+4))
           kmove=kmove+4

! check if angle already specified

           kangle=0
           check=.true.
           Do While (check .and. kangle < Min(jangle,mxangl))
              kangle=kangle+1
              check=.not.( jj   == listang(0,kangle) .and. &
                           iatm == listang(1,kangle) .and. &
                           jatm == listang(2,kangle) .and. &
                           katm == listang(3,kangle) )
           End Do

! insert new angle details

           If (check) Then
              jangle=jangle+1

              If (jangle <= mxangl) Then
                 listang(0,jangle)=jj
                 listang(1,jangle)=iatm
                 listang(2,jangle)=jatm
                 listang(3,jangle)=katm

                 Call tag_legend(safe1,newatm,jangle,legang,mxfang)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kangle,legang,mxfang)
           End If
        End Do
        kmove=kmove+1

! unpack dihedral angle details

        legdih(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           katm=Nint(buffer(kmove+4))
           latm=Nint(buffer(kmove+5))
           kmove=kmove+5

! check if dihedral already specified

           kdihed=0
           check=.true.
           Do While (check .and. kdihed < Min(jdihed,mxdihd))
              kdihed=kdihed+1
              check=.not.( jj   == listdih(0,kdihed) .and. &
                           iatm == listdih(1,kdihed) .and. &
                           jatm == listdih(2,kdihed) .and. &
                           katm == listdih(3,kdihed) .and. &
                           latm == listdih(4,kdihed) )
           End Do

! add new dihedral details

           If (check) Then
              jdihed=jdihed+1

              If (jdihed <= mxdihd) Then
                 listdih(0,jdihed)=jj
                 listdih(1,jdihed)=iatm
                 listdih(2,jdihed)=jatm
                 listdih(3,jdihed)=katm
                 listdih(4,jdihed)=latm

                 Call tag_legend(safe1,newatm,jdihed,legdih,mxfdih)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kdihed,legdih,mxfdih)
           End If
        End Do
        kmove=kmove+1

! unpack inversion angle details

        leginv(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           katm=Nint(buffer(kmove+4))
           latm=Nint(buffer(kmove+5))
           kmove=kmove+5

! check if inversion already specified

           kinver=0
           check=.true.
           Do While (check .and. kinver < Min(jinver,mxinv))
              kinver=kinver+1
              check=.not.( jj   == listinv(0,kinver) .and. &
                           iatm == listinv(1,kinver) .and. &
                           jatm == listinv(2,kinver) .and. &
                           katm == listinv(3,kinver) .and. &
                           latm == listinv(4,kinver) )
           End Do

! add new inversion details

           If (check) Then
              jinver=jinver+1

              If (jinver <= mxinv) Then
                 listinv(0,jinver)=jj
                 listinv(1,jinver)=iatm
                 listinv(2,jinver)=jatm
                 listinv(3,jinver)=katm
                 listinv(4,jinver)=latm

                 Call tag_legend(safe1,newatm,jinver,leginv,mxfinv)
              Else
                 safe=.false.
              End If
           Else
              Call tag_legend(safe1,newatm,kinver,leginv,mxfinv)
           End If
        End Do
        kmove=kmove+1

! redefine bond counters

        ntshl =jshels

        ntcons=jconst
        ntpmf =jpmf

        ntrgd =jrigid

        ntteth=jteths

        ntbond=jbonds
        ntangl=jangle
        ntdihd=jdihed
        ntinv =jinver

     End If
  End Do

! check error flags

  Call gcheck(safe)
  If (.not.safe) Call error(113)
  Call gcheck(safe1)
  If (.not.safe1) Call error(114)

  Deallocate (buffer,      Stat=fail(1))
  Deallocate (i1pmf,i2pmf, Stat=fail(2))
  Deallocate (lrgd,        Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'deport_atomic_data deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine deport_atomic_data
