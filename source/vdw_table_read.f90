Subroutine vdw_table_read(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABLE file (for van der waals forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 1994
! amended   - i.t.todorov may 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : delr_max,ntable,nrite, &
                           mxgvdw,zero_plus,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use vdw_module,   Only : ntpvdw,ls_vdw,lstvdw,ltpvdw,prmvdw,gvdw,vvdw,sigeps
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail,ngrid,katom1,katom2,ivdw,jtpatm,keyvdw,i,j,l
  Real( Kind = wp )      :: delpot,cutpot,dlrpot,rdr,rrr,ppp,vk,vk1,vk2,t,t1,t2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer


  If (idnode == 0) Open(Unit=ntable, File='TABLE')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

! read mesh resolution

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  Call get_word(record,word)
  delpot = word_2_real(word)

  Call get_word(record,word)
  cutpot = word_2_real(word)

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  dlrpot = rvdw/Real(mxgvdw-4,wp)

! check grid spacing

  safe=.false.
  If (Abs(delpot-dlrpot) <= 1.0e-8_wp) Then
     safe=.true.
     delpot=dlrpot
  End If
  If (delpot > delr_max .and. (.not.safe)) Then
     If (idnode == 0) Then
        Write(nrite,"(/,                                             &
             & ' expected (maximum) radial increment : ',1p,e15.7,/, &
             & ' TABLE  file actual radial increment : ',1p,e15.7)") &
             delr_max, delpot
        Write(nrite,"(/,                                                &
             & ' expected (minimum) number of grid points : ',0p,i10,/, &
             & ' TABLE  file actual number of grid points : ',0p,i10)") &
             mxgvdw, ngrid
     End If
     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     If (idnode == 0) Write(nrite,"(/,' TABLE arrays resized for mxgrid = ',i10)") mxgvdw-4
  End If

! compare grids dimensions

  If (ngrid < mxgvdw-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgvdw-4,wp),0.0_wp)
     Call error(48)
  End If

  If (cutpot < rvdw) Call error(504)

  fail=0
  Allocate (buffer(0:ngrid), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_table_read allocation failure, node: ', idnode
     Call error(0)
  End If

! read potential arrays for all pairs

  Do ivdw=1,ntpvdw

! read potential arrays if potential not already defined

     If (ltpvdw(ivdw) == 0) Then

! read pair potential labels and long-range corrections

        Call get_line(safe,ntable,record)
        If (.not.safe) Go To 100

        Call get_word(record,atom1)
        Call get_word(record,atom2)

        Call get_word(record,word)
        prmvdw(1,ivdw)=word_2_real(word)

        Call get_word(record,word)
        prmvdw(2,ivdw)=word_2_real(word)

        katom1=0
        katom2=0

        Do jtpatm=1,ntpatm
           If (atom1 == unqatm(jtpatm)) katom1=jtpatm
           If (atom2 == unqatm(jtpatm)) katom2=jtpatm
        End Do

        If (katom1 == 0 .or. katom2 == 0) Then
           If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'**** entry in TABLE'
           Call error(81)
        End If

        keyvdw=(Max(katom1,katom2)*(Max(katom1,katom2)-1))/2 + Min(katom1,katom2)

! Only one vdw potential per pair is allowed
! (FIELD AND TABLE potentials overlapping)

        If (lstvdw(keyvdw) /= ivdw) Call error(23)

! read in potential arrays

        Do i=1,(ngrid+3)/4
           l=Min(4,ngrid-(i-1)*4)
           If (idnode == 0) Then
              Read(Unit=ntable, Fmt=*, End=100) (buffer((i-1)*4+j),j=1,l)
           Else
              buffer((i-1)*4+1:(i-1)*4+l)=0.0_wp
           End If
        End Do
        If (mxnode > 1) Call MPI_BCAST(buffer(1:ngrid), ngrid, wp_mpi, 0, dlp_comm_world, ierr)

! reconstruct arrays using 3pt interpolation

        If (remake) Then
           Do i=1,mxgvdw-3
              rrr = Real(i,wp)*dlrpot
              l   = Int(rrr*rdr)
              ppp=rrr*rdr-Real(l,wp)

              vk  = buffer(l)

! linear extrapolation for the grid points just beyond the cutoff

              If (l+2 > ngrid) Then
                 If (l+1 > ngrid) Then
                    vk1 = 2.0_wp*buffer(l)-buffer(l-1)
                    vk2 = 2.0_wp*vk1-buffer(l)
                 Else
                    vk1 = buffer(l+1)
                    vk2 = 2.0_wp*buffer(l+1)-buffer(l)
                 End If
              Else
                 vk1 = buffer(l+1)
                 vk2 = buffer(l+2)
              End If

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
              vvdw(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,mxgvdw-4
              vvdw(i,ivdw) = buffer(i)
           End Do

! linear extrapolation for the grid point just beyond the cutoff

           vvdw(mxgvdw-3,ivdw) = 2.0_wp*vvdw(mxgvdw-4,ivdw) - vvdw(mxgvdw-5,ivdw)
        End If

! linear extrapolation for the grid points at 0 and at mxgvdw-2

        vvdw(0,ivdw) = 2.0_wp*buffer(1)-buffer(2)
        vvdw(mxgvdw-2,ivdw) = 2.0_wp*vvdw(mxgvdw-3,ivdw) - vvdw(mxgvdw-4,ivdw)

! read in force arrays

        Do i=1,(ngrid+3)/4
           l=Min(4,ngrid-(i-1)*4)
           If (idnode == 0) Then
              Read(Unit=ntable, Fmt=*, End=100) (buffer((i-1)*4+j),j=1,l)
           Else
              buffer((i-1)*4+1:(i-1)*4+l)=0.0_wp
           End If
        End Do
        If (mxnode > 1) Call MPI_BCAST(buffer(1:ngrid), ngrid, wp_mpi, 0, dlp_comm_world, ierr)

! reconstruct arrays using 3pt interpolation

        If (remake) Then
           Do i=1,mxgvdw-3
              rrr = Real(i,wp)*dlrpot
              l   = Int(rrr*rdr)
              ppp=rrr*rdr-Real(l,wp)

              vk  = buffer(l)

! linear extrapolation for the grid points just beyond the cutoff

              If (l+2 > ngrid) Then
                 If (l+1 > ngrid) Then
                    vk1 = 2.0_wp*buffer(l)-buffer(l-1)
                    vk2 = 2.0_wp*vk1-buffer(l)
                 Else
                    vk1 = buffer(l+1)
                    vk2 = 2.0_wp*buffer(l+1)-buffer(l)
                 End If
              Else
                 vk1 = buffer(l+1)
                 vk2 = buffer(l+2)
              End If

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              gvdw(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,mxgvdw-4
              gvdw(i,ivdw) = buffer(i)
           End Do

! linear extrapolation for the grid point just beyond the cutoff

           gvdw(mxgvdw-3,ivdw) = 2.0_wp*gvdw(mxgvdw-4,ivdw) - gvdw(mxgvdw-5,ivdw)
        End If

! linear extrapolation for the grid points at 0 and at mxgvdw-2

        gvdw(0,ivdw) = (2.0_wp*buffer(1)-0.5_wp*buffer(2))*rdr
        gvdw(mxgvdw-2,ivdw) = 2.0_wp*gvdw(mxgvdw-3,ivdw) - gvdw(mxgvdw-4,ivdw)

! We must distinguish that something has been defined

        If (Abs(vvdw(0,ivdw)) <= zero_plus) vvdw(0,ivdw) = Sign(Tiny(vvdw(0,ivdw)),vvdw(0,ivdw))

     End If

  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABLE file'
  End If

! convert to internal units

  Do ivdw=1,ntpvdw
     If (ltpvdw(ivdw) == 0) Then

! Sigma-epsilon initialisation

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        Do i=0,mxgvdw
           vvdw(i,ivdw)=vvdw(i,ivdw)*engunit
           gvdw(i,ivdw)=gvdw(i,ivdw)*engunit

! Sigma-epsilon search

           If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                       vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                      (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                       vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                       vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                    sigeps(2,ivdw)=-vvdw(i-1,ivdw)
              End If
           End If
        End Do
     End If
  End Do

  If (ls_vdw) Then
     Do ivdw=1,ntpvdw
        If (ltpvdw(ivdw) == 0) Then

! Sigma-epsilon initialisation

           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp

! Sigma-epsilon search

           Do i=1,mxgvdw-4
              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 t  = vvdw(i  ,ivdw) + gvdw(mxgvdw-4,ivdw) * &
                      (Real(i  ,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
                 t1 = vvdw(i-1,ivdw) + gvdw(mxgvdw-4,ivdw) * &
                      (Real(i-1,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,t1)) == -Nint(Sign(1.0_wp,t))) &
                       sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    t2 = vvdw(i-2,ivdw) + gvdw(mxgvdw-4,ivdw) * &
                         (Real(i-2,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)

                    If ( (t2 >= t1 .and. t1 <= t) .and.         &
                         (t2 /= t1 .or. t2 /= t .or. t1 /= t) ) &
                       sigeps(2,ivdw)=-t1
                 End If
              End If
           End Do
           vvdw(mxgvdw-3,ivdw) = 0.0_wp ; vvdw(mxgvdw-2,ivdw) = 0.0_wp
           gvdw(mxgvdw-3,ivdw) = 0.0_wp ; gvdw(mxgvdw-2,ivdw) = 0.0_wp
        End If
     End Do
  End If

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_table_read deallocation failure, node: ', idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine vdw_table_read
