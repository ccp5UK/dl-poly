Subroutine scan_field                                 &
           (l_n_e,mxsite,mxatyp,megatm,mxtmls,mxexcl, &
           mxtshl,mxshl,mxfshl,mxtcon,mxcons,mxfcon,  &
           mxtpmf,mxpmf,mxfpmf,                       &
           mxtrgd,mxrgd,mxlrgd,mxfrgd,                &
           mxtteth,mxteth,mxftet,                     &
           mxtbnd,mxbond,mxfbnd,mxtang,mxangl,mxfang, &
           mxtdih,mxdihd,mxfdih,mxtinv,mxinv,mxfinv,  &
           mxrdf,mxgrid,mxvdw,rvdw,mxmet,rmet,        &
           mxter,rcter,mxtbp,rctbp,mxfbp,rcfbp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of the filed file
!
! copyright - daresbury laboratory
! author    - w.smith november 1994
! amended   - i.t.todorov april 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode,mxnode,gcheck
  Use setup_module, Only : nfield,ntable
  Use parse_module, Only : get_line,get_word,lower_case,word_2_real

  Implicit None

! Max number of different atom types

  Integer, Parameter :: mmk = 1000

! Average maximum number of intra-like bonds per atom

  Integer, Parameter :: mxb = 6

  Character( Len = 8   ), Dimension( 1:mmk ) :: chr

  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: name

  Logical           :: l_n_e,check,ltable,safe
  Integer           :: mxtmls,itmols,nummols,numsit,isite,ksite,nrept,         &
                       mxsite,mxatyp,megatm,i,j,k,nfld,mxexcl,                 &
                       numshl,mtshl,mxtshl,mxshl,ishls,mxfshl,                 &
                       numcon,mtcons,mxtcon,mxcons,icon,mxfcon,                &
                       mxtpmf(1:2),mxpmf,ipmf,jpmf,mxfpmf,                     &
                       numrgd,mtrgd,mxtrgd,mxlrgd,mxrgd,irgd,jrgd,lrgd,mxfrgd, &
                       numteth,mtteth,mxtteth,mxteth,iteth,mxftet,             &
                       numbonds,mtbond,mxtbnd,mxbond,ibonds,mxfbnd,            &
                       numang,mtangl,mxtang,mxangl,iang,mxfang,                &
                       numdih,mtdihd,mxtdih,mxdihd,idih,mxfdih,                &
                       numinv,mtinv,mxtinv,mxinv,iinv,mxfinv,                  &
                       mxrdf,itprdf,mxvdw,itpvdw,mxmet,itpmet,mxgrid,          &
                       mxter,itpter,mxtbp,itptbp,mxfbp,itpfbp,                 &
                       mxt(1:9),mxf(1:9)
  Real( Kind = wp ) :: rvdw,rmet,rcter,rctbp,rcfbp,rct

  l_n_e=.true.

  nummols=0

  numsit=0
  mxsite=0
  mxatyp=0
  megatm=0
  mxtmls=0

  numshl=0
  mtshl =0
  mxshl =0
  mxtshl=0
  mxfshl=1

  numcon=0
  mtcons=0
  mxcons=0
  mxtcon=0
  mxfcon=1

  mxpmf =0
  mxtpmf=0
  mxfpmf=1

  numrgd=0
  mtrgd =0
  mxrgd =0
  mxtrgd=0
  mxlrgd=0
  mxfrgd=1

  numteth=0
  mtteth =0
  mxteth =0
  mxtteth=0
  mxftet =1

  numbonds=0
  mtbond=0
  mxbond=0
  mxtbnd=0
  mxfbnd=1

  numang=0
  mtangl=0
  mxangl=0
  mxtang=0
  mxfang=1

  numdih=0
  mtdihd=0
  mxdihd=0
  mxtdih=0
  mxfdih=1

  numinv=0
  mtinv =0
  mxinv =0
  mxtinv=0
  mxfinv=1

  mxrdf =0

  mxgrid=0

  mxvdw=0
  rvdw =0.0_wp

  mxmet=0
  rmet =0.0_wp

  mxter=0
  rcter=0.0_wp

  mxtbp=0
  rctbp=0.0_wp

  mxfbp=0
  rcfbp=0.0_wp

  mxexcl=0

! Set safe flag

  safe=.true.

! Open the interactions input file

  If (idnode == 0) Inquire(File='FIELD', Exist=safe)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     Go To 20
  Else
     If (idnode == 0) Open(Unit=nfield, File='FIELD', Status='old')
  End If

  Call get_line(safe,nfield,record)
  If (.not.safe) Go To 30

  Do While (.true.)

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nfield,record)
        If (.not.safe) Go To 30
        Call lower_case(record)
        Call get_word(record,word)
     End Do

     If (word(1:7) == 'molecul') Then

        Call get_word(record,word)
        If (word(1:4) == 'type') Call get_word(record,word)
        mxtmls=Nint(word_2_real(word))

        Do itmols=1,mxtmls

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call lower_case(record)
              Call get_word(record,word)
           End Do

           Do While (.true.)

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record)
                 If (.not.safe) Go To 30
                 Call lower_case(record)
                 Call get_word(record,word)
              End Do

              If (word(1:6) == 'nummol') Then

                 Call get_word(record,word)
                 nummols=Nint(word_2_real(word))

              Else If (word(1:5) == 'atoms') Then

                 Call get_word(record,word)
                 numsit=Nint(word_2_real(word))
                 megatm=megatm+numsit*nummols
                 mxsite=mxsite+numsit

                 ksite=0

                 Do isite=1,numsit

                    If (ksite < numsit) Then

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record)
                          If (.not.safe) Go To 30
                          Call get_word(record,word)
                       End Do

                       name=word(1:8)

                       Call get_word(record,word)
                       Call get_word(record,word)
                       l_n_e=(l_n_e.and.(Abs(word_2_real(word)) < 1.0e-5_wp))
                       Call get_word(record,word)
                       nrept=Nint(word_2_real(word))
                       If (nrept == 0) nrept=1

                       If (mxatyp == 0) Then
                          mxatyp=1
                          chr(1)=name
                       Else
                          check=.true.
                          Do j=1,mxatyp
                             If (name == chr(j)) check=.false.
                          End Do

                          If (check) Then
                             mxatyp=mxatyp+1
                             If (mxatyp <= mmk) chr(mxatyp)=name
                          End If
                       End If

                       ksite=ksite+nrept

                    End If

                 End Do

                 If (mmk < mxatyp) Call error(2)

              Else If (word(1:5) == 'shell') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numshl=Nint(word_2_real(word))
                 mtshl=Max(mtshl,numshl)
                 mxtshl=mxtshl+numshl
                 mxshl=mxshl+nummols*numshl

                 Do ishls=1,numshl
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:6) == 'constr') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numcon=Nint(word_2_real(word))
                 mtcons=Max(mtcons,numcon)
                 mxtcon=mxtcon+numcon
                 mxcons=mxcons+nummols*numcon

                 Do icon=1,numcon
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:3) == 'pmf') Then

                 mxpmf=mxpmf+nummols

                 Do ipmf=1,2
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call lower_case(record)
                       Call get_word(record,word)
                    End Do

                    If (word(1:3) == 'pmf') Then
                       Call get_word(record,word)
                    Else
                       Go To 30
                    End If
                    If (word(1:4) == 'unit')  Then
                       Call get_word(record,word)
                    Else
                       Go To 30
                    End If
                    mxtpmf(ipmf)=Nint(word_2_real(word))

                    Do jpmf=1,mxtpmf(ipmf)
                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record)
                          If (.not.safe) Go To 30
                          Call get_word(record,word)
                       End Do
                    End Do
                 End Do

              Else If (word(1:5) == 'rigid') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units' .or. word(1:3) == 'bod') Call get_word(record,word)
                 numrgd=Nint(word_2_real(word))
                 mtrgd=Max(mtrgd,numrgd)
                 mxtrgd=mxtrgd+numrgd
                 mxrgd=mxrgd+nummols*numrgd

                 Do irgd=1,numrgd
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                    jrgd=Nint(word_2_real(word))
                    mxlrgd=Max(mxlrgd,jrgd)

                    Do lrgd=1,jrgd
                       If (Mod(lrgd+1,12) == 0) Then
                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nfield,record)
                             If (.not.safe) Go To 30
                          End Do
                       End If
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:4) == 'teth') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numteth=Nint(word_2_real(word))
                 mtteth=Max(mtteth,numteth)
                 mxtteth=mxtteth+numteth
                 mxteth=mxteth+numteth*nummols

                 Do iteth=1,numteth
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:5) == 'bonds') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numbonds=Nint(word_2_real(word))
                 mtbond=Max(mtbond,numbonds)
                 mxtbnd=mxtbnd+numbonds
                 mxbond=mxbond+nummols*numbonds

                 Do ibonds=1,numbonds
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:6) == 'angles') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numang=Nint(word_2_real(word))
                 mtangl=Max(mtangl,numang)
                 mxtang=mxtang+numang
                 mxangl=mxangl+nummols*numang

                 Do iang=1,numang
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:6) == 'dihedr') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numdih=Nint(word_2_real(word))
                 mtdihd=Max(mtdihd,numdih)
                 mxtdih=mxtdih+numdih
                 mxdihd=mxdihd+nummols*numdih

                 Do idih=1,numdih
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:6) == 'invers') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numinv=Nint(word_2_real(word))
                 mtinv=Max(mtinv,numinv)
                 mxtinv=mxtinv+numinv
                 mxinv=mxinv+nummols*numinv

                 Do iinv=1,numinv
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:6) == 'finish') Then

                 Go To 1000

              End If

           End Do

 1000      Continue

        End Do

     Else If (word(1:3) == 'rdf') Then

        Call get_word(record,word)
        mxrdf=Nint(word_2_real(word))

        Do itprdf=1,mxrdf
            word(1:1)='#'
            Do While (word(1:1) == '#' .or. word(1:1) == ' ')
               Call get_line(safe,nfield,record)
               If (.not.safe) Go To 30
               Call get_word(record,word)
            End Do
        End Do

        mxrdf=Max(mxrdf,(mxatyp*(mxatyp+1))/2)

     Else If (word(1:3) == 'vdw') Then

        ltable=.false.
        Call get_word(record,word)
        If (word(1:3) == 'tab') Then
           ltable=.true.
        Else
           mxvdw=Nint(word_2_real(word))
        End If

        Do itpvdw=1,mxvdw
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call lower_case(record)
              Call get_word(record,word)
           End Do

           Call get_word(record,word)
           Call get_word(record,word)
           If (word(1:3) == 'tab') ltable=.true.
        End Do

        mxvdw=Max(mxvdw,(mxatyp*(mxatyp+1))/2)

        If (ltable) Then
           If (idnode == 0) Open(Unit=ntable, File='TABLE')

           Call get_line(safe,ntable,record)
           If (.not.safe) Go To 40

           Call get_line(safe,ntable,record)
           If (.not.safe) Go To 40
           Call get_word(record,word)

           Call get_word(record,word)
           rvdw=Max(rvdw,word_2_real(word))

           Call get_word(record,word)
           mxgrid=Nint(word_2_real(word))

           If (idnode == 0) Close(Unit=ntable)
        End If

     Else If (word(1:3) == 'met') Then

        ltable=.false.

        Call get_word(record,word)
        mxmet=Nint(word_2_real(word))

        Do itpmet=1,mxmet
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call lower_case(record)
              Call get_word(record,word)
           End Do

           Call get_word(record,word)
           Call get_word(record,word)
           If (word(1:3) == 'eam') ltable=.true.
        End Do

        mxmet=Max(mxmet,(mxatyp*(mxatyp+1))/2)

        If (ltable) Then
           If (idnode == 0) Open(Unit=ntable, File='TABEAM')

           Call get_line(safe,ntable,record)
           If (.not.safe) Go To 40
           Call get_line(safe,ntable,record)
           If (.not.safe) Go To 40
           Call get_word(record,word)

           Do i=1,Nint(word_2_real(word))
              Call get_line(safe,ntable,record)
              If (.not.safe) Go To 40
              Call get_word(record,word)
              Call lower_case(word)
              j=0
              If (word(1:4) == 'embe') j=1
              If (word(1:4) == 'pair') Call get_word(record,word)
              Call get_word(record,word)

              Call get_word(record,word)
              k=Nint(word_2_real(word))
              mxgrid=Max(mxgrid,k+4)

              Call get_word(record,word)
              Call get_word(record,word)
              If (j == 0) rmet=Max(rmet,word_2_real(word))

              Do j=1,(k+3)/4
                 Call get_line(safe,ntable,record)
                 If (.not.safe) Go To 40
              End Do
           End Do

           If (idnode == 0) Close(Unit=ntable)
        End If

     Else If (word(1:7) == 'tersoff') Then

        Call get_word(record,word)
        mxter=Nint(word_2_real(word))

        Do itpter=1,mxter
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           rct=word_2_real(word)
           rcter=Max(rcter,rct)
        End Do

        Do itpter=1,(mxter*(mxter+1)) / 2
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do
        End Do

        mxter=Max(mxter,(mxatyp*(mxatyp+1))/2)

     Else If (word(1:3) == 'tbp') Then

        Call get_word(record,word)
        mxtbp=Nint(word_2_real(word))

        Do itptbp=1,mxtbp
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           Do i=1,8
              Call get_word(record,word)
           End Do

           rct=word_2_real(word)
           rctbp=Max(rctbp,rct)
        End Do

     Else If (word(1:3) == 'fbp') Then

        Call get_word(record,word)
        mxfbp=Nint(word_2_real(word))

        Do itpfbp=1,mxfbp
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           Do i=1,7
              Call get_word(record,word)
           End Do

           rct=word_2_real(word)
           rcfbp=Max(rcfbp,rct)
        End Do

     Else If (word(1:6) == 'extern') Then

        Call get_word(record,word)
        nfld=0

        If (word(1:1) /= '#' .and. word(1:1) /= ' ') nfld=Nint(word_2_real(word))
        If (nfld <= 0) nfld=5

        word(1:1)='#'
        Do While (word(1:1) == '#' .or. word(1:1) == ' ')
           Call get_line(safe,nfield,record)
           If (.not.safe) Go To 30
           Call get_word(record,word)
        End Do

     Else If (word(1:5) == 'close') Then

        Go To 10

     End If

  End Do

10 Continue
  If (idnode == 0) Close(Unit=nfield)

! Estimate per-node maxima and define legend arrays length,
! where length is length+1 for the final zero used used as
! end of read in deport_particles

  mxshl=Max(mxshl,mxnode*mtshl)
  If (mxshl >  0) mxfshl=mxfshl+1
  mxf(1)=mxfshl

  mxcons=Max(mxcons,mxnode*mtcons)
  If (mxcons > 0) mxfcon=mxfcon+mxb
  mxf(2)=mxfcon

  If (mxpmf  > 0) mxfpmf=mxfpmf+1 ! PMFs are global
  mxf(3)=mxfpmf

  mxrgd=Max(mxrgd,mxnode*mtrgd)
  If (mxrgd  > 0) mxfrgd=mxfrgd+1
  mxf(4)=mxlrgd

  mxteth=Max(mxteth,mxnode*mtteth)
  If (mxteth > 0) mxftet=mxftet+1
  mxf(5)=mxftet

  mxbond=Max(mxbond,mxnode*mtbond)
  If (mxbond > 0) mxfbnd=mxfbnd+(mxb*(mxb+1))/2
  mxf(6)=mxfbnd

  mxangl=Max(mxangl,mxnode*mtangl)
  If (mxangl > 0) mxfang=mxfang+(mxb*(mxb+1))/2
  mxf(7)=mxfang

  mxdihd=Max(mxdihd,mxnode*mtdihd)
  If (mxdihd > 0) mxfdih=mxfdih+((mxb-1)*mxb*(mxb+1))/2
  mxf(8)=mxfdih

  mxinv=Max(mxinv,mxnode*mtinv)
  If (mxinv  > 0) mxfinv=mxfinv+(mxb*(mxb+1))/4
  mxf(9)=mxfinv

  Do i=1,9
     mxt(i)=Min(1,mxf(i))
  End Do
  Call shellsort(9,mxf)
  mxexcl=(mxf(9)+Sum(mxf)/Max(1,Sum(mxt)))*(Min(1,mxshl)+1)+1

! (vdw,met) = rdf scanning

  If (mxrdf == 0 .and. (mxvdw > 0 .or. mxmet > 0)) mxrdf = Max(mxvdw,mxmet)

  Return

! FILED file does not exist

20 Continue
  Call error(122)
  Return

30 Continue
  If (idnode == 0) Close(Unit=nfield)
  Call error(52)
  Return

40 Continue
  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine scan_field
