Subroutine scan_field                                &
           (l_n_e,mxompl,mximpl,                     &
           mxsite,mxatyp,megatm,mxtmls,mxexcl,       &
           mtshl,mxtshl,mxshl,mxfshl,                &
           mtcons,mxtcon,mxcons,mxfcon,              &
           mxtpmf,mxpmf,mxfpmf,                      &
           mtrgd,mxtrgd,mxrgd,mxlrgd,mxfrgd,         &
           mtteth,mxtteth,mxteth,mxftet,             &
           mtbond,mxtbnd,mxbond,mxfbnd,rcbnd,mxgbnd, &
           mtangl,mxtang,mxangl,mxfang,mxgang,       &
           mtdihd,mxtdih,mxdihd,mxfdih,mxgdih,       &
           mtinv,mxtinv,mxinv,mxfinv,mxginv,         &
           mxrdf,mxvdw,rvdw,mxgvdw,                  &
           mxmet,mxmed,mxmds,rmet,mxgmet,            &
           mxter,rcter,mxtbp,rctbp,mxfbp,rcfbp,lext)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of the FIELD file
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith february 2016
! contrib   - b.palmer (2band) may 2013
! contrib   - a.v.brukhno & i.t.todorov march 2014 (itramolecular TPs)
! contrib   - h.a.boateng february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gcheck
  Use setup_module,      Only : nrite,nfield,ntable
  Use parse_module,      Only : get_line,strip_blanks, &
                                get_word,lower_case,word_2_real
  Use bonds_module,      Only : lt_bnd
  Use angles_module,     Only : lt_ang
  Use dihedrals_module,  Only : lt_dih
  Use inversions_module, Only : lt_inv
  Use vdw_module,        Only : lt_vdw
  Use metal_module,      Only : tabmet
  Use tersoff_module,    Only : potter
  Use kim_module,        Only : kim,rkim,kim_cutoff

  Implicit None

! Max number of different atom types

  Integer, Parameter :: mmk = 1000

! Average maximum number of intra-like bonds per atom

  Integer, Parameter :: mxb = 6

  Character( Len = 8   ), Dimension( 1:mmk ) :: chr

  Character( Len = 200 ) :: record,record_raw
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: name

  Logical           :: l_n_e,check,safe,lext
  Integer           :: mxtmls,itmols,nummols,numsit,mxnmst,isite,ksite,nrept,  &
                       mxompl,mximpl,mxsite,mxatyp,megatm,i,j,k,mxexcl,        &
                       numshl,mtshl,mxtshl,mxshl,ishls,mxfshl,                 &
                       numcon,mtcons,mxtcon,mxcons,icon,mxfcon,                &
                       mxtpmf(1:2),mxpmf,ipmf,jpmf,mxfpmf,                     &
                       numrgd,mtrgd,mxtrgd,mxlrgd,mxrgd,irgd,jrgd,lrgd,mxfrgd, &
                       numteth,mtteth,mxtteth,mxteth,iteth,mxftet,             &
                       numbonds,mtbond,mxtbnd,mxbond,ibonds,mxfbnd,mxgbnd,     &
                       numang,mtangl,mxtang,mxangl,iang,mxfang,mxgang,         &
                       numdih,mtdihd,mxtdih,mxdihd,idih,mxfdih,mxgdih,         &
                       numinv,mtinv,mxtinv,mxinv,iinv,mxfinv,mxginv,           &
                       mxrdf,itprdf,mxvdw,itpvdw,mxgvdw,                       &
                       mxmet,mxmed,mxmds,itpmet,mxgmet,                        &
                       mxter,itpter,mxtbp,itptbp,mxfbp,itpfbp,                 &
                       mxt(1:9),mxf(1:9)
  Real( Kind = wp ) :: rcbnd,rvdw,rmet,rcter,rctbp,rcfbp,rct,tmp,tmp1,tmp2

  l_n_e=.true.  ! no electrostatics opted
  mxompl=0      ! default of maximum order of poles (charges)
  mximpl=0      ! default maximum number of independent poles values
                ! it initialises to 0 if no MULT directive exists in FIELD

  nummols=0

  numsit=0
  mxnmst=0
  mxsite=0
  mxatyp=0
  megatm=0
  mxtmls=0

  numshl=0
  mtshl =0
  mxshl =0
  mxtshl=0
  mxfshl=0

  numcon=0
  mtcons=0
  mxcons=0
  mxtcon=0
  mxfcon=0

  mxpmf =0
  mxtpmf=0
  mxfpmf=0

  numrgd=0
  mtrgd =0
  mxrgd =0
  mxtrgd=0
  mxlrgd=0
  mxfrgd=0

  numteth=0
  mtteth =0
  mxteth =0
  mxtteth=0
  mxftet =0

  numbonds=0
  mtbond=0
  mxbond=0
  mxtbnd=0
  mxfbnd=0
  rcbnd =0.0_wp
  mxgbnd=-2

  numang=0
  mtangl=0
  mxangl=0
  mxtang=0
  mxfang=0
  mxgang=-2

  numdih=0
  mtdihd=0
  mxdihd=0
  mxtdih=0
  mxfdih=0
  mxgdih=-2

  numinv=0
  mtinv =0
  mxinv =0
  mxtinv=0
  mxfinv=0
  mxginv=-2

  mxrdf =0

  mxvdw =0
  rvdw  =0.0_wp
  mxgvdw=0

  mxmet =0
  mxmed =0
  mxmds =0
  rmet  =0.0_wp
  mxgmet=0

  mxter=0
  rcter=0.0_wp

  mxtbp=0
  rctbp=0.0_wp

  mxfbp=0
  rcfbp=0.0_wp

  mxexcl=0

  lext=.false.

! Set safe flag

  safe=.true.

! Open the interactions input file

  If (idnode == 0) Inquire(File='FIELD', Exist=safe)
  If (mxnode > 1) Call gcheck(safe,"enforce")
  If (.not.safe) Then
     Go To 20
  Else
     If (idnode == 0) Open(Unit=nfield, File='FIELD', Status='old')
  End If

  Call get_line(safe,nfield,record)
  If (.not.safe) Go To 30

  Do

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nfield,record)
        record_raw = record ! KIM needs case sensitive IM name
        If (.not.safe) Go To 30
        Call lower_case(record)
        Call get_word(record,word)
     End Do

! multipoles container detection

     If (word(1:9) == 'multipole') Then

        l_n_e=.false.

        Call get_word(record,word)
        mxompl = Max(0,Nint(word_2_real(word)))

        mximpl = (mxompl+3)*(mxompl+2)*(mxompl+1)/6

        If (idnode == 0 .and. mxompl > 4) &
  Write(nrite,'(1x,a,i0)') "Extending electrostatics to multipolar interactions of order ", mxompl

     Else If (word(1:7) == 'molecul') Then

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

           Do

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
                 mxnmst=Max(mxnmst,numsit)
                 megatm=megatm+nummols*numsit
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
                       If (Mod(lrgd,16) == 0) Then
                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nfield,record)
                             If (.not.safe) Go To 30
                             Call get_word(record,word)
                          End Do
                       Else
                          Call get_word(record,word)
                       End If
                    End Do
                 End Do

              Else If (word(1:4) == 'teth') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numteth=Nint(word_2_real(word))
                 mtteth=Max(mtteth,numteth)
                 mxtteth=mxtteth+numteth
                 mxteth=mxteth+nummols*numteth

                 Do iteth=1,numteth
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:5) == 'bonds') Then

!                 lt_bnd=.false. ! initialised in bonds_module.f90

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

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) lt_bnd=.true.
                 End Do

                 If (lt_bnd) Then
                    If (idnode == 0) Open(Unit=ntable, File='TABBND')

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABBND if it's in .xvg format

                    Call get_word(record,word)
                    rcbnd=Max(rcbnd,word_2_real(word))

                    Call get_word(record,word)
                    k=Nint(word_2_real(word))
                    mxgbnd=Max(mxgbnd,k+4)

                    If (idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'angles') Then

!                 lt_ang=.false. ! initialised in angles_module.f90

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

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) lt_ang=.true.
                 End Do

                 If (lt_ang) Then
                    If (idnode == 0) Open(Unit=ntable, File='TABANG')

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABANG if it's in .xvg format

                    Call get_word(record,word) ! no need for cutoff in angles (max is always 180 degrees)
                    k=Nint(word_2_real(word))
                    mxgang=Max(mxgang,k+4)

                    If (idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'dihedr') Then

!                 lt_dih=.false. ! initialised in dihedrals_module.f90

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

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) lt_dih=.true.
                 End Do

                 If (lt_dih) Then
                    If (idnode == 0) Open(Unit=ntable, File='TABDIH')

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABDIH if it's in .xvg format

                    Call get_word(record,word) ! no need for cutoff in angles (max is always 360 degrees
                    k=Nint(word_2_real(word))  ! from -180 to 180)
                    mxgdih=Max(mxgdih,k+4)

                    If (idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'invers') Then

!                 lt_dih=.false. ! initialised in dihedrals_module.f90

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

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) lt_inv=.true.
                 End Do

                 If (lt_inv) Then
                    If (idnode == 0) Open(Unit=ntable, File='TABINV')

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

                    Call get_word(record,word) ! no need for cutoff in angles (max is always 180 degrees)
                    k=Nint(word_2_real(word))
                    mxginv=Max(mxginv,k+4)

                    If (idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'finish') Then

                 Go To 1000

              End If

           End Do

1000       Continue

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

        If (mxrdf > 0) mxrdf=Max(mxrdf,(mxatyp*(mxatyp+1))/2)

     Else If (word(1:3) == 'vdw') Then

!        lt_vdw=.false. ! initialised in vdw_module

        Call get_word(record,word)
        If (word(1:3) == 'tab') Then
           lt_vdw=.true.
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

           If      (word(1:3) == 'tab') Then
              lt_vdw=.true.
           Else If (word(1:3) == 'snm') Then
              Call get_word(record,word)
              Call get_word(record,word)
              Call get_word(record,word)
              Call get_word(record,word)
              rvdw=Max(rvdw,word_2_real(word))
           Else If (word(1:3) == 'wca') Then
              Call get_word(record,word)
              Call get_word(record,word) ; tmp1=word_2_real(word)
              Call get_word(record,word) ; tmp2=word_2_real(word)
              tmp=tmp2+tmp1*2.0_wp**(1.0_wp/6.0_wp)
              rvdw=Max(rvdw,tmp)
           Else If (word(1:3) == 'dpd') Then
              Call get_word(record,word)
              rvdw=Max(rvdw,word_2_real(word))
           End If
        End Do

        If (mxvdw > 0) Then
           mxvdw=Max(mxvdw,(mxatyp*(mxatyp+1))/2)

           If (lt_vdw) Then
              If (idnode == 0) Open(Unit=ntable, File='TABLE')

              Call get_line(safe,ntable,record)
              If (.not.safe) Go To 40

              Call get_line(safe,ntable,record)
              If (.not.safe) Go To 40
              Call get_word(record,word)

              Call get_word(record,word)
              rvdw=Max(rvdw,word_2_real(word))

              Call get_word(record,word)
              k=Nint(word_2_real(word))
              mxgvdw=Max(mxgvdw,k)

              If (idnode == 0) Close(Unit=ntable)
           End If
        End If

     Else If (word(1:3) == 'met') Then

!        tabmet=-1 ! initialised in metal_module

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
           tabmet=0 ! for FST metal potentials
           If      (word(1:3) ==  'eam') Then
              tabmet=1
           Else If (word(1:4) == 'eeam') Then
              tabmet=2
           Else If (word(1:4) == '2bea') Then
              tabmet=3
           Else If (word(1:4) == '2bee') Then
              tabmet=4
           Else If (word(1:4) == 'fnsc') Then
              Call get_word(record,word) ; Call get_word(record,word)
              Call get_word(record,word) ; Call get_word(record,word)
              rmet=Max(rmet,word_2_real(word))
              Call get_word(record,word) ; Call get_word(record,word)
              rmet=Max(rmet,word_2_real(word))
           Else If (word(1:4) == 'exfs') Then
              Call get_word(record,word) ; Call get_word(record,word)
              Call get_word(record,word) ; Call get_word(record,word)
              Call get_word(record,word) ; Call get_word(record,word)
              rmet=Max(rmet,word_2_real(word))
              Call get_word(record,word) ; Call get_word(record,word)
              rmet=Max(rmet,word_2_real(word))
           End If
        End Do

        If (mxmet > 0) Then
           mxmet=Max(mxmet,(mxatyp*(mxatyp+1))/2)

           If      (tabmet == 0) Then
              mxmed=mxmet
           Else If (tabmet == 1) Then
              mxmed=mxatyp
           Else If (tabmet == 2) Then
              mxmed=mxatyp**2
           Else If (tabmet == 3) Then
              mxmed=mxatyp
              mxmds=mxatyp*(mxatyp+1)/2
           Else If (tabmet == 4) Then
              mxmed=mxatyp**2
              mxmds=mxatyp**2
           End If

           If (tabmet > 0) Then
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
                 j=0 ! assume rmet is specified
                 If (word(1:4) == 'embe' .or. & ! 2-band embedding functionals
                     word(1:4) == 'demb' .or. word(1:4) == 'semb') j=1 ! no rmet is specified
                 If ((word(1:4) == 'dens' .and. tabmet == 2) .or. & ! EEAM
                     (word(2:4) == 'den' .and. (tabmet == 3 .or. tabmet == 4)) .or. & ! sden & dden for 2B extensions
                     word(1:4) == 'pair') Call get_word(record,word) ! skip over one species
                 Call get_word(record,word)                          ! skip over one species

                 Call get_word(record,word)
                 k=Nint(word_2_real(word))
                 mxgmet=Max(mxgmet,k+4)

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
        End If

     Else If (word(1:7) == 'tersoff') Then

        Call get_word(record,word)
        mxter=Nint(word_2_real(word))

        Do itpter=1,mxter
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call lower_case(record)
              Call get_word(record,word)
           End Do

           Call get_word(record,word)
           If      (word(1:4) == 'ters') Then
              potter=1
           Else If (word(1:4) == 'kihs') Then
              potter=2
           End If

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           rct=word_2_real(word)
           rcter=Max(rcter,rct)

           If (potter == 2) Then
              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record)
                 If (.not.safe) Go To 30
                 Call get_word(record,word)
              End Do
           End If
        End Do

        If (mxter > 0) Then
           If (potter == 1) Then
              Do itpter=1,(mxter*(mxter+1))/2
                 word(1:1)='#'
                 Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                    Call get_line(safe,nfield,record)
                    If (.not.safe) Go To 30
                    Call get_word(record,word)
                 End Do
              End Do
           End If

           mxter=Max(mxter,(mxatyp*(mxatyp+1))/2)
        End If

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

     Else If (word(1:7) == 'kim') Then

! Get KIM's IM name and cutoff

        Call get_word(record_raw,word)
        Call strip_blanks(record_raw)
        kim=record(1:Len_Trim(record_raw))
        Call kim_cutoff(mxatyp,chr,kim,rkim)

     Else If (word(1:6) == 'extern') Then

        lext=.true.

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

! Define legend arrays lengths.  If length > 0 then
! length=Max(length)+1 for the violation excess element

  If (mxshl >  0) mxfshl=1+1 ! One shell per core
  mxf(1)=mxfshl

  If (mxcons > 0) mxfcon=mxb+1
  mxf(2)=mxfcon

  If (mxpmf  > 0) mxfpmf=1+1 ! PMFs are global
  mxf(3)=mxfpmf

  If (mxrgd  > 0) mxfrgd=1+1 ! One RB per particle
  mxf(4)=mxlrgd

  If (mxteth > 0) mxftet=1+1 ! One tether per particle
  mxf(5)=mxftet

  If (mxbond > 0) mxfbnd=(mxb*(mxb+1))+1
  mxf(6)=mxfbnd

  If (mxangl > 0) mxfang=(mxb+1)**2/2+1
  mxf(7)=mxfang

  If (mxdihd > 0) mxfdih=((mxb-1)*mxb*(mxb+1))/2+2*mxb+1
  mxf(8)=mxfdih

  If (mxinv  > 0) mxfinv=(mxb*(mxb+1))/4+1
  mxf(9)=mxfinv

  Do i=1,9
     mxt(i)=Min(1,mxf(i))
  End Do
  mxexcl = Min( mxnmst , Max( mxfrgd , Sum(mxf)/Max(1,Sum(mxt)) ) * (Min(1,mxshl)+1) )
  If (mxexcl > 0) mxexcl=mxexcl+1 ! violation excess element

  Return

! FIELD file does not exist

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
