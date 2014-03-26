Subroutine vdw_table_read(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABLE file (for van der waals forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 1994
! amended   - i.t.todorov december 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : ntable,nrite,mxgrid,engunit
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

  fail=0
  Allocate (buffer(1:mxgrid), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_table_read allocation failure, node: ', idnode
     Call error(0)
  End If

  remake=.false.

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

  dlrpot=rvdw/Real(mxgrid-4,wp)

  safe=.false.
  If (Abs(delpot-dlrpot) <= 1.0e-8_wp) Then
     safe=.true.
     delpot=dlrpot
  End If
  If ( (delpot > dlrpot .and. (.not.safe)) .or. &
       (ngrid-4 /= Nint(cutpot/delpot)) ) Then
     If (idnode == 0) Then
        Write(nrite,"(/,                                          &
           & 'expected (minimum) radial increment : ',1p,e15.7,/, &
           & 'TABLE file         radial increment : ',1p,e15.7)") &
           dlrpot, delpot
        Write(nrite,"(/,                                             &
           & 'expected (minimum) number of grid points : ',0p,i10,/, &
           & 'TABLE file stated  number of grid points : ',0p,i10,/, &
           & 'TABLE file derived number of grid points : ',0p,i10)") &
           mxgrid, ngrid, Nint(cutpot/delpot)+4
     End If
     Call error(22)
  End If
  safe=.true.

  If (cutpot < rvdw) Call error(504)
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     If (idnode == 0) Write(nrite,"(/,' TABLE arrays resized for mxgrid = ',i10)") mxgrid
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

! check array dimensions

        If (ngrid > mxgrid) Then
           Call warning(270,Real(ngrid,wp),Real(mxgrid,wp),0.0_wp)
           Call error(48)
        End If

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
           rdr=1.0_wp/delpot
           Do i=1,mxgrid
              rrr = Real(i,wp)*dlrpot
              l   = Int(rrr*rdr)

              ppp=rrr*rdr-Real(l,wp)
              vk  = buffer(l)
              vk1 = buffer(l+1)
              vk2 = buffer(l+2)

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
              vvdw(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,mxgrid
              vvdw(i,ivdw) = buffer(i)
           End Do
        End If

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
           Do i=1,mxgrid
              rrr = Real(i,wp)*dlrpot
              l   = Int(rrr*rdr)

              ppp=rrr*rdr-Real(l,wp)
              vk  = buffer(l)
              vk1 = buffer(l+1)
              vk2 = buffer(l+2)

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              gvdw(i,ivdw) = t1 + (t2-t1)*ppp*0.5_wp
           End Do
        Else
           Do i=1,mxgrid
              gvdw(i,ivdw) = buffer(i)
           End Do
        End If

     End If

  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABLE file'
  End If

! convert to internal units

  Do l=1,ntpvdw
     If (ltpvdw(l) == 0) Then

! Sigma-epsilon initialisation

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        Do i=1,mxgrid
           vvdw(i,l)=vvdw(i,l)*engunit
           gvdw(i,l)=gvdw(i,l)*engunit

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
     Do l=1,ntpvdw
        If (ltpvdw(l) == 0) Then

! Sigma-epsilon initialisation

           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp

! Sigma-epsilon search

           Do i=1,mxgrid
              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 t  = gvdw(mxgrid,ivdw)*(Real(i  ,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgrid,ivdw)
                 t1 = gvdw(mxgrid,ivdw)*(Real(i-1,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgrid,ivdw)

                 vk  = vvdw(i  ,ivdw) + t
                 vk1 = vvdw(i-1,ivdw) + t1

                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Sign(1.0_wp,vk1) == -Sign(1.0_wp,vk)) &
                          sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    t2 = gvdw(mxgrid,ivdw)*(Real(i-2,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgrid,ivdw)
                    vk2 = vvdw(i-2,ivdw) + ppp

                    If ( (vk2 >= vk1 .and. vk1 <= vk) .and.           &
                         (vk2 /= vk1 .or. vk2 /= vk .or. vk1 /= vk) ) &
                       sigeps(2,ivdw)=-vk1
                 End If
              End If
           End Do
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
