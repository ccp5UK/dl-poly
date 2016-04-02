Subroutine read_mpoles(l_top,sumchg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the molecular mulitpole
! specifications of the system to be simulated
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module,  Only : nrite,nmpldt,mxompl

! SITE & MPOLES MODULE

  Use site_module
  Use mpoles_module, Only : mpllfr

! PARSE MODULE

  Use parse_module

  Implicit None

  Logical,           Intent ( In    ) :: l_top
  Real( Kind = wp ), Intent ( InOut ) :: sumchg

  Logical                :: safe, l_ord=.false.

  Character( Len = 200 ) :: record,record1,record2
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: atom

  Integer                :: itmols,nrept,i,j,k,l,                 &
                            isite,jsite,ksite,lsite,nsite,sitmpl, &
                            ordmpl,ordmpl_start,ordmpl_next,      &
                            ordmpl_min,ordmpl_max,                &
                            indmpl,indmpl_start,indmpl_final

  Real( Kind = wp )      :: Factorial,charge,scl

! open MPOLES data file

  If (idnode == 0) Then
     Open(Unit=nmpldt, File = 'MPOLES', Status = 'old')
     Write(nrite,"(/,/,1x,'ELECTROSTATICS MULTIPOLES SPECIFICATION')")
     If (.not.l_top) Write(nrite,"(/,1x,'detailed specification opted out')")
  End If

  Call get_line(safe,nmpldt,record)
  If (.not.safe) Go To 2000

! omit first line

  nsite  = 0
  sumchg = 0.0_wp

  ordmpl_min = 4
  ordmpl_max = 0

! read and process directives from mpols/field file

  Do

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nmpldt,record)
        If (.not.safe) Go To 2000
        Call get_word(record,word) ; Call lower_case(word)
     End Do

! specify molecular species

     If (word(1:7) == 'molecul') Then

        Call get_word(record,word) ; Call lower_case(word)
        If (word(1:4) == 'type') Call get_word(record,word)

        If (ntpmls == Nint(word_2_real(word))) Then
           If (idnode == 0) Write(nrite,"(/,/,1x,'number of molecular types',6x,i10)") ntpmls
        Else
           If (idnode == 0) Write(nrite,'(/,1x,a,2(/,1x,a,i0))')                        &
  "*** warning - number of molecular types mistmatch between FIELD and MPOLES !!! ***", &
  "***           FIELD  reports: ", ntpmls,                                             &
  "***           MPOLES reports: ", Nint(word_2_real(word))

           Call error(623)
        End If

! read in molecular characteristics for every molecule

        Do itmols=1,ntpmls

           If (idnode == 0 .and. l_top) Write(nrite,"(/,/,1x,'molecular species type',9x,i10)") itmols

! name of molecular species

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nmpldt,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do
           Call strip_blanks(record)
           record1=word(1:Len_Trim(word)+1)//record ; Call lower_case(record1)
           record2=molnam(itmols) ;                   Call lower_case(record2)

           If (record1 == record2) Then
              If (idnode == 0 .and. l_top) Write(nrite,"(/,1x,'name of species:',13x,a40)") molnam(itmols)
           Else
              If (idnode == 0) Write(nrite,'(/,1x,a,i0)') &
  "*** warning - molecular names mistmatch between FIELD and MPOLES for type !!! *** ", itmols

              Call error(623)
           End If

! read molecular data

           Do

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nmpldt,record)
                 If (.not.safe) Go To 2000
                 Call get_word(record,word) ; Call lower_case(word)
              End Do

! number of molecules of this type

              If (word(1:6) == 'nummol') Then

                 Call get_word(record,word)

                 If (nummols(itmols) == Nint(word_2_real(word))) Then
                    If (idnode == 0 .and. l_top) Write(nrite,"(/,1x,'number of molecules  ',10x,i10)") nummols(itmols)
                 Else
                    If (idnode == 0) Write(nrite,'(/,1x,a,2(/,1x,a,i0))')         &
  "*** warning - number of molecules mistmatch between FIELD and MPOLES !!! ***", &
  "***           FIELD  reports: ", nummols(itmols),                              &
  "***           MPOLES reports: ", Nint(word_2_real(word))

                    Call error(623)
                 End If

! read in atomic details

              Else If (word(1:5) == 'atoms') Then

                 Call get_word(record,word)

                 If (numsit(itmols) == Nint(word_2_real(word))) Then
                    If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,1x,'number of atoms/sites',10x,i10)") numsit(itmols)
  Write(nrite,"(/,1x,'atomic characteristics:', &
       & /,/,15x,'site',4x,'name',2x,'multipolar order',2x,'repeat'/)")
                    End If
                 Else
                    If (idnode == 0) Write(nrite,'(/,1x,a,2(/,1x,a,i0))')                        &
  "*** warning - number of atoms/sites per molecule mistmatch between FIELD and MPOLES !!! ***", &
  "***           FIELD  reports: ", numsit(itmols),                                              &
  "***           MPOLES reports: ", Nint(word_2_real(word))

                    Call error(623)
                 End If

! for every molecule of this type

! get site and atom description

! reference point

                 ksite=0
                 Do isite=1,numsit(itmols)
                    If (ksite < numsit(itmols)) Then

! read atom name, highest pole order supplied, repeat

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nmpldt,record)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       atom=word(1:8)

! read supplied pole order

                       Call get_word(record,word)
                       ordmpl=Abs(Nint(word_2_real(word)))
                       indmpl=(ordmpl+3)*(ordmpl+2)*(ordmpl+1)/6

! get the min and max order defined

                       ordmpl_min=Min(ordmpl_min,ordmpl)
                       ordmpl_max=Max(ordmpl_max,ordmpl)

                       Call get_word(record,word)
                       nrept=Abs(Nint(word_2_real(word)))
                       If (nrept == 0) nrept=1

                       jsite=nsite+1
                       lsite=jsite+nrept-1

                       Do i=jsite,lsite
                          If (sitnam(i) /= atom) Then
                             If (idnode == 0) Write(nrite,'(/,1x,a,i0,a)') &
  "*** warning - site names mistmatch between FIELD and MPOLES for site ", ksite+1+i-jsite, " !!! ***"

                             Call error(623)
                          End If
                       End Do

                       If (idnode == 0 .and. l_top) &
  Write(nrite,"(9x,i10,4x,a8,4x,i2,5x,i10)") ksite+1,atom,ordmpl,nrept

! monopole=charge

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nmpldt,record)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       sitmpl = 1

                       charge=word_2_real(word)

                       chgsit(jsite:lsite)=charge
                       mpllfr(sitmpl,jsite:lsite)=charge

! sum absolute charges

                       sumchg=sumchg+Abs(charge)

! report

                       If (idnode == 0 .and. l_top) &
  Write(nrite,"(3x,a12,3x,f10.5)") 'charge',charge

! higher poles counters

                       ordmpl_start = 0
                       ordmpl_next  = ordmpl_start+1
                       indmpl_start = sitmpl+1
                       indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                       Do While (ordmpl_next <= ordmpl)

! read line per pole order

                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nmpldt,record)
                             If (.not.safe) Go To 2000
                             Call get_word(record,word)
                          End Do

! Only assign what FIELD says is needed

                          If (ordmpl_next <= mxompl) Then

                             Do i=indmpl_start,indmpl_final
                                sitmpl = sitmpl+1
                                mpllfr(sitmpl,jsite:lsite)=word_2_real(word)
                                Call get_word(record,word)
                             End Do

! report

                             If (idnode == 0 .and. l_top) Then
                                If      (ordmpl_next == 1) Then
  Write(nrite,"(3x,a12,3x, 3f10.5)") 'dipole',       mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 2) Then
  Write(nrite,"(3x,a12,3x, 6f10.5)") 'quadrupole',   mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 3) Then
  Write(nrite,"(3x,a12,3x,10f10.5)") 'octupole',     mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 4) Then
  Write(nrite,"(3x,a12,3x,15f10.5)") 'hexadecapole', mpllfr(indmpl_start:indmpl_final,jsite)
                                End If
                             End If

! rescale poles values by their degeneracy

                             If (ordmpl_next > 1) Then
                                sitmpl=sitmpl-(indmpl_final-indmpl_start+1) ! rewind
                                Do i=ordmpl_next,0,-1
                                   l=ordmpl_next-i
                                   Do j=l,0,-1
                                      k=l-j

                                      scl=Exp(Factorial(ordmpl_next)-Factorial(k)-Factorial(j)-Factorial(i))
!                                      Write(*,*) i,j,k,Nint(scl)

                                      sitmpl = sitmpl+1 ! forward and apply scaling if degeneracy exists
                                      If (Nint(scl) /= 1) mpllfr(sitmpl,jsite:lsite)=mpllfr(sitmpl,jsite:lsite)/scl
                                   End Do
                                End Do
                             End If

                          Else

                             l_ord=.true.

! update actual order marker

                             sitmpl=indmpl_final

! report

                             If (idnode == 0 .and. l_top) Then
                                If      (ordmpl_next == 1) Then
  Write(nrite,"(3x,a12,1x,a)") 'dipole',                 '     *** supplied but not required ***'
                                Else If (ordmpl_next == 2) Then
  Write(nrite,"(3x,a12,1x,a)") 'quadrupole',             '     *** supplied but not required ***'
                                Else If (ordmpl_next == 3) Then
  Write(nrite,"(3x,a12,1x,a)") 'octupole',               '     *** supplied but not required ***'
                                Else If (ordmpl_next == 4) Then
  Write(nrite,"(3x,a12,1x,a)") 'hexadecapole',           '     *** supplied but not required ***'
                                Else
  Write(nrite,"(3x,a12,i0,a)") 'pole order ',ordmpl_next,'     *** supplied but not required ***'
                                End If
                             End If

                          End If

! update poles counters

                          ordmpl_next  = ordmpl_next+1
                          indmpl_start = sitmpl+1
                          indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                       End Do

                       nsite=nsite+nrept
                       ksite=ksite+nrept

                    End If
                 End Do

! finish of data for one molecular type

              Else If (word(1:6) == 'finish') Then

                 If (idnode == 0) Then
                    Write(nrite,'(/,1x,3(a,i0),a)') &
  "*** warning - multipolar electrostatics requested up to order ", &
  ordmpl, " with specified interactions up order ",                 &
  ordmpl_max," and least order ", ordmpl_min," !!! ***"
                    If (ordmpl_max*ordmpl == 0) Write(nrite,'(1x,2a)') &
  "*** warning - multipolar electrostatics machinery to be used for ", &
  "monopoles only electrostatic interactions (point charges only) !!! ***"
                    If (ordmpl_max > 4) Write(nrite,'(1x,2a)')     &
  "*** warning - electrostatic interactions beyond hexadecapole ", &
  "order can not be considered and are thus ignored !!! ***"
                 End If

                 Go To 1000

              Else

! error exit for unidentified directive in molecular data

                 Call strip_blanks(record)
                 If (idnode == 0) Write(nrite,'(/,1x,2a)') word(1:Len_Trim(word)+1),record
                 Call error(12)

              End If

           End Do

! just finished with this type molecule data

1000       Continue

        End Do

! close MPOLES data file

     Else If (word(1:5) == 'close') Then

        If (idnode == 0) Close(Unit=nmpldt)

! EXIT IF ALL IS OK

        Return

     Else

! error exit for unidentified directive

        If (idnode == 0) Write(nrite,'(/,1x,a)') word(1:Len_Trim(word))
        Call error(4)

     End If

  End Do

  Return

2000 Continue

  If (idnode == 0) Close(Unit=nmpldt)
  Call error(52)

End Subroutine read_mpoles
