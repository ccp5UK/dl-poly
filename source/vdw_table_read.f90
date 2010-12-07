Subroutine vdw_table_read(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABLE file (for van der waals forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 1994
! amended   - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode,mxnode,gsum
  Use setup_module, Only : ntable,nrite,mxgrid,mxbuff,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use vdw_module,   Only : ntpvdw,lstvdw,ltpvdw,prmvdw,gvdw,vvdw
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail,ngrid,katom1,katom2,ivdw,jtpatm,keyvdw,i,j,l
  Real( Kind = wp )      :: delpot,cutpot,dlrpot,rdr,rrr,ppp,vk,vk1,vk2,t1,t2

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  Allocate (buffer(1:mxbuff), Stat=fail)
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

  If (Abs(delpot-dlrpot) <= 1.0e-8_wp) delpot=dlrpot
  If ( (delpot > dlrpot) .or. (ngrid-4 /= Nint(cutpot/delpot)) ) Then
     If (idnode == 0) Write(nrite,"(                    &
        & 'expected radial increment : ',1p,e15.7,/,    &
        & 'TABLE    radial increment : ',1p,e15.7,/,/,  &
        & 'expected number of grid points : ',0p,i10,/, &
        & 'grid points in TABLE           : ',i10)") dlrpot, delpot, mxgrid, ngrid

     Call error(22)
  End If

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
           If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'****'
           Call  error(81)
        End If

        keyvdw=(Max(katom1,katom2)*(Max(katom1,katom2)-1))/2 + Min(katom1,katom2)

! Only one vdw potential per pair is allowed
! (FIELD AND TABLE potentials overlapping)

        If (lstvdw(keyvdw) /= ivdw) Call error(23)

! check array dimensions

        If (ngrid > mxbuff) Then
           Call warning(270,Real(ngrid,wp),Real(mxbuff,wp),0.0_wp)
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
        If (mxnode > 1) Call gsum(buffer(1:ngrid))

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
        If (mxnode > 1) Call gsum(buffer(1:ngrid))

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
        Do i=1,mxgrid
           vvdw(i,l)=vvdw(i,l)*engunit
           gvdw(i,l)=gvdw(i,l)*engunit
        End Do
     End If
  End Do

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
