Subroutine read_field                      &
           (imcon,l_n_v,l_str,l_top,       &
           rcut,rvdw,rmet,width,keyfce,    &
           lbook,lexcl,keyshl,             &
           rcter,rctbp,rcfbp,              &
           atmfre,atmfrz,megatm,megfrz,    &
           megshl,megcon,megpmf,megrgd,    &
           megtet,megbnd,megang,megdih,meginv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the molecular specifications
! of the system to be simulated
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES

  Use kinds_f90
  Use comms_module, Only : idnode
  Use setup_module

! SITE MODULE

  Use site_module

! CONFIG MODULE
! Fuch's correction of charge non-neutral systems
! Global_To_Local variables

  Use config_module, Only : sumchg,gtl_b

! INTERACTION MODULES

  Use core_shell_module

  Use constraints_module
  Use pmf_module

  Use rigid_bodies_module

  Use tethers_module

  Use bonds_module
  Use angles_module
  Use dihedrals_module
  Use inversions_module

  Use vdw_module
  Use metal_module
  Use tersoff_module
  Use three_body_module
  Use four_body_module

  Use external_field_module

! STATTISTICS MODULE

  Use statistics_module       ! for rdf

! PARSE MODULE

  Use parse_module

  Implicit None

  Logical,           Intent( In    ) :: l_n_v,l_str,l_top
  Integer,           Intent( In    ) :: imcon
  Integer,           Intent( InOut ) :: keyfce
  Real( Kind = wp ), Intent( In    ) :: rcut,rvdw,rmet,width

  Logical,           Intent(   Out ) :: lbook,lexcl
  Real( Kind = wp ), Intent(   Out ) :: rcter,rctbp,rcfbp
  Integer,           Intent(   Out ) :: keyshl,                             &
                                        atmfre,atmfrz,megatm,megfrz,        &
                                        megshl,megcon,megpmf,megrgd,        &
                                        megtet,megbnd,megang,megdih,meginv

  Logical                :: safe,lunits,lmols,atmchk,                        &
                            l_shl,l_con,l_rgd,l_tet,l_bnd,l_ang,l_dih,l_inv, &
                            lshl_one,lshl_all,lpmf,ltable,lmet_safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 4   ) :: keyword
  Character( Len = 8   ) :: atom0,atom1,atom2,atom3
  Integer                :: itmols,isite,jsite,ksite,msite,nsite,                &
                            isite1,isite2,isite3,isite4,is(1:4),js(1:4),         &
                            irept,nrept,ifrz,ntmp,i,j,ia,ja,jtpatm,ntab,keypot,  &
                            iatm1,iatm2,iatm3,iatm4,katom0,katom1,katom2,katom3, &
                            ishls,nshels,icnst,nconst,frzcon,ipmf,jpmf,kpmf,     &
                            nrigid,irgd,jrgd,krgd,lrgd,frzrgd,iteth,nteth,       &
                            ibond,nbonds,iang,nangle,idih,ndihed,iinv,ninver,    &
                            itprdf,keyrdf,itpvdw,keyvdw,itpmet,keymet,           &
                            itpter,keyter,icross,                                &
                            itbp,itptbp,keytbp,ktbp,                             &
                            ifbp,itpfbp,keyfbp,ka1,ka2,ka3,kfbp,nfld
  Real( Kind = wp )      :: weight,charge,pmf_tmp(1:2),parpot(1:30),tmp



! Initialise number of unique atom and shell types and of different types of molecules

  ntpatm = 0
  ntpshl = 0
  ntpmls = 0

! Default flag for existence of molecules

  lmols =.false.

! Initialise energy units for interactions flag

  lunits  = .false.

! Initialise global atomic site index and running counters of
! shells, constraints, PMF, RBs,
! tethers, bonds, angles, dihedrals and inversions in the system

  nsite  = 0

  nshels = 0
  nconst = 0
  lpmf   = .false. ! no PMF defined yet
  nrigid = 0

  nteth  = 0

  nbonds = 0
  nangle = 0
  ndihed = 0
  ninver = 0

! Initialise total number of particles(shells are particles in MD),
! frozen particles, shells, constraints, PMFs, RBs,
! tethers, bonds, angles, dihedrals and inversions in the system.
! Some are used as switches for various force calculations.

  megatm = 0
  megfrz = 0

  megshl = 0

  megcon = 0
  megpmf = 0

  megrgd = 0

  megtet = 0

  megbnd = 0
  megang = 0
  megdih = 0
  meginv = 0

! Default for existence of intra-like interactions (including
! shells and tethers)

  lbook = .false.

! Default flag for existence of excluded intra-interaction

  lexcl = .false.

! default shell model - none

  keyshl   = 1
  lshl_one = .false. ! A massless shell existence indicator
  lshl_all = .true.  ! All shells are massless indicator

! Default flags for existence a table for vdw interactions

  ltable=.false.

! Default potential cutoffs

  rcter=0.0_wp
  rctbp=0.0_wp
  rcfbp=0.0_wp


! open force field data file

  If (idnode == 0) Then
     Open(Unit=nfield, File = 'FIELD', Status = 'old')
     Write(nrite,"(/,/,1x,'SYSTEM SPECIFICATION')")
     If (.not. l_top) Write(nrite,"(/,1x,'detailed tiopology opted out')")
  End If

! omit first line

  Call get_line(safe,nfield,record)
  If (.not.safe) Go To 2000

! read and process directives from field file

  Do While (.true.)

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nfield,record)
        If (.not.safe) Go To 2000
        Call lower_case(record)
        Call get_word(record,word)
     End Do

! energy unit for input/output

     If (word(1:5) == 'units') Then

        lunits=.true.
        Call get_word(record,word)

        If (word(1:2) == 'ev') Then

           engunit = 9648.530821_wp
           If (idnode == 0) Write(nrite,"(/,1x,'energy units = eV')")

        Else If (word(1:4) == 'kcal') Then

           engunit = 418.4_wp
           If (idnode == 0) Write(nrite,"(/,1x,'energy units = kcal/mol')")

        Else If (word(1:2) == 'kj') Then

           engunit = 100.0_wp
           If (idnode == 0) Write(nrite,"(/,1x,'energy units = kJ/mol')")

        Else If (word(1:8) == 'internal') Then

           If (idnode == 0) Write(nrite,"(/,1x,'energy units = dl_poly internal units (10 J/mol)')")

        Else If (word(1:1) == 'k') Then

           engunit = boltz
           If (idnode == 0) Write(nrite,"(/,1x,'energy units = Kelvin/Boltzmann')")

        Else If (word(1:1) == ' ') Then

           If (idnode == 0) Write(nrite,"(/,1x,'energy units = dl_poly internal units (10 J/mol)')")

        Else

           If (idnode == 0) Write(nrite,'(/,1x,a)') word(1:Len_Trim(word)+1)
           Call error(5)

        End If

! neutral group control option

     Else If (word(1:4) == 'neut') Then

        Call error(26)

! specify molecular species

     Else If (word(1:7) == 'molecul') Then

        Call get_word(record,word)
        If (word(1:4) == 'type') Call get_word(record,word)

! number of molecular types

        If (lmols) Call error(11)
        lmols=.true.
        ntpmls=Nint(word_2_real(word))

        If (idnode == 0) Write(nrite,"(/,/,1x,'number of molecular types',6x,i10)") ntpmls

        If (ntpmls > mxtmls) Call error(10)

! read in molecular characteristics for every molecule

        Do itmols=1,ntpmls

! expectation values for once off definitions of bonded quantities
! PMFs make an exception (as defined once and only once per system)

           l_shl=.true. ; l_con=.true. ; l_rgd=.true. ; l_tet=.true.
           l_bnd=.true. ; l_ang=.true. ; l_dih=.true. ; l_inv=.true.

           If (idnode == 0 .and. l_top) Write(nrite,"(/,/,1x,'molecular species type',9x,i10)") itmols

! name of molecular species

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do
           Call strip_blanks(record)
           molnam(itmols)=word(1:Len_Trim(word)+1)//record

           If (idnode == 0 .and. l_top) Write(nrite,"(/,1x,'name of species:',13x,a40)") molnam(itmols)

! stop processing if energy unit has not been specified

           If (.not.lunits) Call error(6)

! read molecular data

           Do While (.true.)

! initialise frozen constraints & RBs counters

              frzcon=0
              frzrgd=0

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record)
                 If (.not.safe) Go To 2000
                 Call lower_case(record)
                 Call get_word(record,word)
              End Do

! number of molecules of this type

              If (word(1:6) == 'nummol') Then

                 Call get_word(record,word)
                 nummols(itmols)=Nint(word_2_real(word))

                 If (idnode == 0 .and. l_top) Write(nrite,"(/,1x,'number of molecules  ',10x,i10)") nummols(itmols)

! read in atomic details

              Else If (word(1:5) == 'atoms') Then

                 Call get_word(record,word)
                 numsit(itmols)=Nint(word_2_real(word))

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,1x,'number of atoms/sites',10x,i10)") numsit(itmols)
  Write(nrite,"(/,1x,'atomic characteristics:', &
       & /,/,15x,'site',4x,'name',13x,'mass',9x,'charge',6x,'repeat',6x,'freeze',/)")
                 End If

! for every molecule of this type

! get site and atom description

! reference point

                 ksite=0

                 Do isite=1,numsit(itmols)

                    If (ksite < numsit(itmols)) Then

! read atom name, mass, charge, freeze option

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       atom1=word(1:8)

                       Call get_word(record,word)
                       weight=word_2_real(word)

                       Call get_word(record,word)
                       charge=word_2_real(word)

                       sumchg=sumchg+Abs(charge)

                       Call get_word(record,word)
                       nrept=Nint(word_2_real(word))
                       If (nrept == 0) nrept=1

                       Call get_word(record,word)
                       ifrz=Nint(word_2_real(word))
                       If (ifrz /= 0) ifrz=1

                       numfrz(itmols)=numfrz(itmols)+ifrz*nrept

                       If (idnode == 0 .and. l_top) &
  Write(nrite,"(9x,i10,4x,a8,2f15.6,2i10)") ksite+1,atom1,weight,charge,nrept,ifrz

                       Do irept=1,nrept
                          nsite=nsite+1
                          If (nsite > mxsite) Call error(20)

                          sitnam(nsite)=atom1
                          wgtsit(nsite)=weight
                          chgsit(nsite)=charge
                          frzsit(nsite)=ifrz
                          If (wgtsit(nsite) > 1.0E-6_wp) dofsit(nsite)=3.0_wp*Real(Abs(1-ifrz),wp)
                       End Do
                       ksite=ksite+nrept

! establish list of unique atom types

                       atmchk=.true.
                       Do jsite=1,ntpatm
                          If (atom1 == unqatm(jsite)) Then
                             atmchk=.false.

                             Do irept=nsite,nsite-nrept+1,-1
                                typsit(irept)=jsite
                             End Do
                          End If
                       End Do

                       If (atmchk) Then
                          ntpatm=ntpatm+1
                          If (ntpatm > mxatyp) Call error(14)

                          unqatm(ntpatm)=atom1

                          Do irept=nsite,nsite-nrept+1,-1
                             typsit(irept)=ntpatm
                          End Do
                       End If
                    End If
                 End Do

! read interaction/field units

! read core - shell spring parameters

              Else If (word(1:5) == 'shell') Then

                 If (.not.l_shl) Call error(477)
                 l_shl=.false.

                 lbook=.true.
                 lexcl=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numshl(itmols)=numshl(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of core-shell units',5x,i10)") ntmp
  Write(nrite,"(/,1x,'core-shell details:', &
       & /,/,18x,'unit',5x,'index',5x,'index',6x,'parameter',/)")
                 End If

                 Do ishls=1,numshl(itmols)
                    nshels=nshels+1
                    If (nshels > mxtshl) Call error(57)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read core & shell atom indices

                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))

                    lstshl(1,nshels)=iatm1
                    lstshl(2,nshels)=iatm2

                    Call get_word(record,word)
                    prmshl(nshels)=word_2_real(word)

! test for frozen core-shell unit

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1)*frzsit(isite2) /= 0) Then
  Write(nrite,"(4x,a8,3i10,2i10,f15.6)") &
       '*frozen*',ishls,lstshl(1,nshels),lstshl(2,nshels),prmshl(nshels)
                       Else
  Write(nrite,"(12x,3i10,f15.6)") &
                   ishls,lstshl(1,nshels),lstshl(2,nshels),prmshl(nshels)
                       End If
                    End If

! abort if a shell is frozen

                    If (frzsit(isite2) /= 0) Call error(49)

! establish list of unique shell types (most certainly ntpshl < ntpatm < mxatyp)

                    If (.not.Any(unqshl(1:ntpshl) == sitnam(isite2))) Then
                       ntpshl=ntpshl+1
                       unqshl(ntpshl)=sitnam(isite2)

                       If (ntpshl > mxatyp) Call error(14)
                    End If

! There is a massless shell, all shells are massless

                    lshl_one=lshl_one .or.  (wgtsit(isite2) < 1.0e-6_wp)
                    lshl_all=lshl_all .and. (wgtsit(isite2) < 1.0e-6_wp)

! test for mistyped core-shell unit (core must be /= shell)

                    If (iatm1 == iatm2) Call error(32)

! convert energy units to internal units

                    prmshl(nshels)=prmshl(nshels)*engunit
                    smax=Max(smax,prmshl(nshels))
                 End Do

! Check for mixed or multiple core-shell entries

                 Do i=nshels-numshl(itmols)+1,nshels
                    is(1)=lstshl(1,i)
                    is(2)=lstshl(2,i)

                    Do j=i+1,nshels
                       js(1)=lstshl(1,j)
                       js(2)=lstshl(2,j)

                       If (js(1) == is(1) .or. js(2) == is(2) .or. &
                           js(1) == is(2) .or. js(2) == is(1)) Then
                          Call warning(390,Real(i,wp),Real(j,wp),0.0_wp)
                          Call error(620)
                       End If
                    End Do
                 End Do

! read constraint bond atom indices and bondlength

              Else If (word(1:6) == 'constr') Then

                 If (.not.l_con) Call error(112)
                 l_con=.false.
                 m_con=1

                 lbook=.true.
                 lexcl=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numcon(itmols)=numcon(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of bond constraints',5x,i10)") ntmp
  Write(nrite,"(/,1x,'constraint bond details:', &
       & /,/,18x,'unit',5x,'index',5x,'index',7x,'bondlength',/)")
                 End If

                 Do icnst=1,numcon(itmols)
                    nconst=nconst+1
                    If (nconst > mxtcon) Call error(40)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read constraint bond atom indices

                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))

                    lstcon(1,nconst)=iatm1
                    lstcon(2,nconst)=iatm2

                    Call get_word(record,word)
                    prmcon(nconst)=word_2_real(word)

! test for frozen atom pairs

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2

! number of completely frozen constraints

                    If (frzsit(isite1)+frzsit(isite2) == 2) Then
                       frzcon=frzcon+1
                    Else If (frzsit(isite1) == 1) Then
                       dofsit(isite2)=dofsit(isite2)-1.0_wp
                    Else If (frzsit(isite2) == 1) Then
                       dofsit(isite1)=dofsit(isite1)-1.0_wp
                    Else
                       dofsit(isite2)=dofsit(isite2)-0.5_wp
                       dofsit(isite1)=dofsit(isite1)-0.5_wp
                    End If

                    If (dofsit(isite1) < -zero_plus) Then
                       Call warning(308,Real(isite1,wp),Real(icnst,wp),Real(itmols,wp))
                       Call error(646)
                    End If

                    If (dofsit(isite2) < -zero_plus) Then
                       Call warning(308,Real(isite2,wp),Real(icnst,wp),Real(itmols,wp))
                       Call error(646)
                    End If

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1)*frzsit(isite2) /= 0) Then
  Write(nrite,"(4x,a8,3i10,f15.6)") &
       '*frozen*',icnst,lstcon(1,nconst),lstcon(2,nconst),prmcon(nconst)
                       Else
  Write(nrite,"(12x,3i10,f15.6)") &
                  icnst,lstcon(1,nconst),lstcon(2,nconst),prmcon(nconst)
                       End If
                    End If

! test for mistyped constraint bond unit

                    If (iatm1 == iatm2) Call error(33)

! test for length of constraint bond unit > rcut

                    If (prmcon(nconst) >= rcut) Call error(34)

                 End Do

! Check for multiple constraint bond entries

                 Do i=nconst-numcon(itmols)+1,nconst
                    is(1)=Min(lstcon(1,i),lstcon(2,i))
                    is(2)=Max(lstcon(1,i),lstcon(2,i))

                    Do j=i+1,nconst
                       js(1)=Min(lstcon(1,j),lstcon(2,j))
                       js(2)=Max(lstcon(1,j),lstcon(2,j))

                       If (js(1) == is(1) .and. js(2) == is(2)) Then
                          Call warning(400,Real(i,wp),Real(j,wp),0.0_wp)
                          Call error(620)
                       End If
                    End Do
                 End Do

! read PMF bond atom indices, weights and bondlength

              Else If (word(1:3) == 'pmf') Then

                 lbook=.true.

                 If (lpmf) Call error(484)
                 lpmf=.true.               ! Only one PMF type per
                 numpmf(itmols)=1          ! MD system is allowed

! read PMF bondlength

                 Call get_word(record,word)
                 prmpmf=word_2_real(word)

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'PMF constraint details')")
  Write(nrite,"(/,5x,'bondlength:',5x,f15.6)") prmpmf
                 End If

                 If (prmpmf > width/2.0_wp) Call error(480)

                 Do ipmf=1,2
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call lower_case(record)
                       Call get_word(record,word)
                    End Do

! read PMF indices and weights

                    pmffrz(ipmf) = 0
                    pmf_tmp(ipmf) = 0.0_wp

! test for zero length units

                    If (mxtpmf(ipmf) == 0) Then
                       Call strip_blanks(record)
                       If (idnode == 0) Write(nrite,"(/,/,2a)") word(1:Len_Trim(word)+1),record
                       Call error(500)
                    End If

                    Do jpmf=1,mxtpmf(ipmf)
                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       iatm1=Nint(word_2_real(word))
                       lstpmf(jpmf,ipmf)=iatm1
                       isite1 = nsite - numsit(itmols) + iatm1

! test for frozen units

                       pmffrz(ipmf)=pmffrz(ipmf)+frzsit(isite1)

                       Call get_word(record,word)
                       weight=word_2_real(word)
                       pmfwgt(jpmf,ipmf)=weight

! test for weightless units

                       pmf_tmp(ipmf)=pmf_tmp(ipmf)+weight
                    End Do
                 End Do

! get real masses
! if a PMF unit is massless supply weights from atoms' masses

                 Do ipmf=1,2
                    Do jpmf=1,mxtpmf(ipmf)
                       isite1 = nsite - numsit(itmols) + lstpmf(jpmf,ipmf)

                       pmfwg1(jpmf,ipmf)=wgtsit(isite1)
                       If (pmf_tmp(ipmf) < 1.0e-6_wp) pmfwgt(jpmf,ipmf)=wgtsit(isite1)
                    End Do

! if a PMF unit is still weightless set all members' masses to 1

                    If (pmf_tmp(ipmf) < 1.0e-6_wp) Then
                       pmf_tmp(ipmf)=Sum(pmfwgt(1:mxtpmf(ipmf),ipmf))

                       If (pmf_tmp(ipmf) < 1.0e-6_wp) Then
                          pmfwgt(:,ipmf)=1.0_wp
                          pmfwg1(:,ipmf)=1.0_wp
                          Call warning(230,Real(ipmf,wp),0.0_wp,0.0_wp)
                       End If
                    End If
                 End Do

! test for frozen atoms and print units

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,18x,'unit',5x,'index',8x,'weight',/)")

                    Do ipmf=1,2
                       Do jpmf=1,mxtpmf(ipmf)
                          isite1 = nsite - numsit(itmols) + lstpmf(jpmf,ipmf)
                          If (frzsit(isite1) /= 0) Then
  Write(nrite,"(4x,a8,2i10,f15.6)") &
       '*frozen*',ipmf,lstpmf(jpmf,ipmf),pmfwgt(jpmf,ipmf)
                          Else
  Write(nrite,"(12x,2i10,f15.6)") &
                  ipmf,lstpmf(jpmf,ipmf),pmfwgt(jpmf,ipmf)
                          End If
                       End Do
                       If (ipmf == 1) Write(nrite, Fmt=*)
                    End Do
                 End If

! PMF reciprocal total unit masses

                 Do ipmf=1,2
                    pmfwgt(0,ipmf)=1.0_wp/Sum(pmfwgt(1:mxtpmf(ipmf),ipmf))
                    pmfwg1(0,ipmf)=1.0_wp/Sum(pmfwg1(1:mxtpmf(ipmf),ipmf))
                 End Do

! abort if there are frozen atoms on both PMF units

                 If (pmffrz(1)*pmffrz(2) > 0) Then
                    Call error(486)
                 Else
                    If (pmffrz(1) == mxtpmf(1) .or. pmffrz(2) == mxtpmf(2)) Then
                       charge=1.0_wp
                    Else
                       charge=0.5_wp
                    End If

                    Do ipmf=1,2
                       ntmp=mxtpmf(ipmf)-pmffrz(ipmf)
                       If (ntmp > 0) Then
                          tmp=charge/Real(ntmp,wp)
                          Do jpmf=1,mxtpmf(ipmf)
                             isite1 = nsite - numsit(itmols) + lstpmf(jpmf,ipmf)
                             If (frzsit(isite1) == 0) Then
                                dofsit(isite1)=dofsit(isite1)-tmp

                                If (dofsit(isite1) < -zero_plus) Then
                                   Call warning(309,Real(isite1,wp),Real(jpmf,wp),Real(ipmf,wp))
                                   Call error(949)
                                End If
                             End If
                          End Do
                       End If
                    End Do
                 End If

! Check for mistyped PMF unit's entries and joined PMF units

                 Do ipmf=1,2
                    Do jpmf=1,mxtpmf(ipmf)
                       Do kpmf=jpmf+1,mxtpmf(ipmf)
                          If (lstpmf(kpmf,ipmf) == lstpmf(jpmf,ipmf)) Call error(501)
                       End Do
                    End Do
                 End Do
                 Do jpmf=1,mxtpmf(2)
                    If (Any(lstpmf(1:mxtpmf(1),1) == lstpmf(jpmf,2))) Call error(502)
                 End Do

! read RBs: number, size, indices

              Else If (word(1:5) == 'rigid') Then

                 If (.not.l_rgd) Call error(625)
                 l_rgd=.false.
                 m_rgd=1

                 lbook=.true.
                 lexcl=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units' .or. word(1:3) == 'bod') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numrgd(itmols)=numrgd(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of rigid bodies',9x,i10)") ntmp
  Write(nrite,"(/,1x,'rigid body details:', &
       & /,/,18x,'unit',6x,'size',12x,'indices',/)")
                 End If

                 Do irgd=1,numrgd(itmols)
                    nrigid=nrigid+1
                    If (nrigid > mxtrgd) Call error(630)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read RB size and indices

                    lrgd=Nint(word_2_real(word))
                    If (lrgd < 2 .or. lrgd > mxlrgd) Call error(632)

                    Do jrgd=1,lrgd
                       If (Mod(jrgd+1,16) == 0) Then
                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nfield,record)
                             If (.not.safe) Go To 2000
                          End Do
                       End If
                       Call get_word(record,word)

                       iatm1=Nint(word_2_real(word))
                       lstrgd(jrgd,nrigid)=iatm1

                       isite1 = nsite - numsit(itmols) + iatm1

! test for frozen and weightless atoms

                       rgdfrz(jrgd,nrigid)=frzsit(isite1)
                       rgdwgt(jrgd,nrigid)=wgtsit(isite1)
                    End Do
                    lstrgd(0,nrigid)=lrgd
                    rgdfrz(0,nrigid)=Sum(rgdfrz(1:lrgd,nrigid))
                    If (rgdfrz(0,nrigid) /= 0) Then
                       If (rgdfrz(0,nrigid) == lrgd) frzrgd=frzrgd+1
                       Do jrgd=1,lrgd
                          rgdwg1(jrgd,nrigid)=Real(rgdfrz(jrgd,nrigid),wp)
                       End Do
                       Do jrgd=1,lrgd
                          rgdwgt(0,nrigid)=rgdwgt(0,nrigid) + &
                                           Real(rgdfrz(jrgd,nrigid)-1,wp)*rgdwgt(jrgd,nrigid)
                       End Do
                       rgdwg1(0,nrigid)=Sum(rgdwg1(1:lrgd,nrigid))
                    Else
                       rgdwgt(0,nrigid)=Sum(rgdwgt(1:lrgd,nrigid))
                       rgdwg1(0:lrgd,nrigid)=rgdwgt(0:lrgd,nrigid)
                    End If


! print RB

                    If (idnode == 0 .and. l_top) Then
                       If (rgdfrz(0,nrigid) /= 0) Then
  Write(nrite,"(4x,a8,12i10,100(/,36x,10i10))") &
       '*frozen*',irgd,(lstrgd(i,nrigid),i=0,lrgd)
                       Else
  Write(nrite,"(12x,12i10,100(/,36x,10i10))") &
                  irgd,(lstrgd(i,nrigid),i=0,lrgd)
                       End If

! test for weightless RB

                       If (rgdwg1(0,nrigid) < 1.0e-6_wp) Call error(634)

! test for frozen RB
!
!                       If (rgdfrz(0,nrigid) > 0) Call error(636)
                    End If

! test for mistyped RB unit

                    Do jrgd=1,lrgd
                       Do krgd=jrgd+1,lrgd
                          If (lstrgd(krgd,nrigid) == lstrgd(jrgd,nrigid)) Call error(638)
                       End Do
                    End Do
                 End Do

! Check for duplicate or joined RB entries

                 Do irgd=nrigid-numrgd(itmols)+1,nrigid
                    lrgd=lstrgd(0,irgd)

                    Do jrgd=irgd+1,nrigid
                       krgd=lstrgd(0,jrgd)
                       If ((Any(lstrgd(1:krgd,jrgd) == lstrgd(1:lrgd,irgd)))) Then
                           Call warning(400,Real(irgd,wp),Real(jrgd,wp),0.0_wp)
                           Call error(620)
                       End If
                    End Do
                 End Do

! read tethered atom indices and tethering parameters

              Else If (word(1:4) == 'teth') Then

                 If (.not.l_tet) Call error(240)
                 l_tet=.false.

                 lbook=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numteth(itmols)=numteth(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of tethered sites',7x,i10)") ntmp
  Write(nrite,"(/,1x,'tethered site details:', &
       & /,/,18x,'unit',5x,'key',6x,'site',19x,'parameters',/) ")
                 End If

                 Do iteth=1,numteth(itmols)
                    nteth=nteth+1
                    If (nteth > mxteth) Call error(62)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of tethering

                    Call lower_case(word)
                    keyword=word(1:4)

                    If      (keyword == 'harm') Then
                       keytet(nteth)=1
                    Else If (keyword == 'rhrm') Then
                       keytet(nteth)=2
                    Else If (keyword == 'quar') Then
                       keytet(nteth)=3
                    Else

                       If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
                       Call error(450)

                    End If

! read tethered atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))

                    lsttet(nteth)=iatm1

                    Call get_word(record,word)
                    prmtet(1,nteth)=word_2_real(word)
                    Call get_word(record,word)
                    prmtet(2,nteth)=word_2_real(word)
                    Call get_word(record,word)
                    prmtet(3,nteth)=word_2_real(word)

! test for frozen atom

                    isite1 = nsite - numsit(itmols) + iatm1

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1) /= 0) Then
  Write(nrite,"(4x,a8,i10,a8,i10,2x,10f15.6)") &
       '*frozen*',iteth,keyword,lsttet(nteth),(prmtet(ja,nteth),ja=1,mxpteth)
                       Else
  Write(nrite,"(12x,i10,a8,i10,2x,10f15.6)") &
                  iteth,keyword,lsttet(nteth),(prmtet(ja,nteth),ja=1,mxpteth)
                       End If
                    End If

! convert energy units to internal units

                    prmtet(:,nteth)=prmtet(:,nteth)*engunit

                 End Do

! Check for multiple tether entries

                 Do i=nteth-numteth(itmols)+1,nteth
                    is(1)=lsttet(i)

                    Do j=i+1,nteth
                       js(1)=lsttet(j)

                       If (js(1) == is(1)) Then
                          If (l_str) Call warning(410,Real(i,wp),Real(j,wp),0.0_wp)
!                          Call error(620)
                       End If
                    End Do
                 End Do

! read chemical bond potential parameters

              Else If (word(1:5) == 'bonds') Then

                 If (.not.l_bnd) Call error(36)
                 l_bnd=.false.

                 lbook=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numbonds(itmols)=numbonds(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of chemical bonds',7x,i10)") ntmp
  Write(nrite,"(/,1x,'chemical bond details:', &
       & /,/,18x,'unit',5x,'key',5x,'index',5x,'index',28x,'parameters',/)")
                 End If

                 Do ibond=1,numbonds(itmols)
                    nbonds=nbonds+1
                    If (nbonds > mxtbnd) Call error(30)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of chemical bond

                    Call lower_case(word)
                    keyword=word(1:4)

                    If (keyword(1:1) /= '-') lexcl=.true.

                    If      (keyword == 'harm') Then
                       keybnd(nbonds)=1
                    Else If (keyword == '-hrm') Then
                       keybnd(nbonds)=-1
                    Else If (keyword == 'mors') Then
                       keybnd(nbonds)=2
                    Else If (keyword == '-mrs') Then
                       keybnd(nbonds)=-2
                    Else If (keyword == '12-6') Then
                       keybnd(nbonds)=3
                    Else If (keyword == '-126') Then
                       keybnd(nbonds)=-3
                    Else If (keyword == 'lj'  ) Then
                       keybnd(nbonds)=4
                    Else If (keyword == '-lj' ) Then
                       keybnd(nbonds)=-4
                    Else If (keyword == 'rhrm') Then
                       keybnd(nbonds)=5
                    Else If (keyword == '-rhm') Then
                       keybnd(nbonds)=-5
                    Else If (keyword == 'quar') Then
                       keybnd(nbonds)=6
                    Else If (keyword == '-qur') Then
                       keybnd(nbonds)=-6
                    Else If (keyword == 'buck') Then
                       keybnd(nbonds)=7
                    Else If (keyword == '-bck') Then
                       keybnd(nbonds)=-7
                    Else If (keyword == 'coul') Then
                       keybnd(nbonds)=8
                    Else If (keyword == '-cul') Then
                       keybnd(nbonds)=-8
                    Else If (keyword == 'fene') Then
                       keybnd(nbonds)=9
                    Else If (keyword == '-fne') Then
                       keybnd(nbonds)=-9
                    Else

                       If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
                       Call error(444)

                    End If

! read bond atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))

                    lstbnd(1,nbonds)=iatm1
                    lstbnd(2,nbonds)=iatm2

                    Call get_word(record,word)
                    prmbnd(1,nbonds)=word_2_real(word)
                    Call get_word(record,word)
                    prmbnd(2,nbonds)=word_2_real(word)
                    Call get_word(record,word)
                    prmbnd(3,nbonds)=word_2_real(word)
                    Call get_word(record,word)
                    prmbnd(4,nbonds)=word_2_real(word)

                    If (Abs(keybnd(nbonds)) == 9) Then
                       prmbnd(2,nbonds)=Abs(prmbnd(2,nbonds))
                       If (Abs(prmbnd(3,nbonds)) > prmbnd(2,nbonds)/2.0_wp) &
                          prmbnd(3,nbonds)=Sign(1.0_wp,prmbnd(3,nbonds))*prmbnd(2,nbonds)/2.0_wp
                    End If

! test for frozen atom pairs

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1)*frzsit(isite2) /= 0) Then
  Write(nrite,"(4x,a8,i10,a8,2i10,2x,10f15.6)") &
       '*frozen*',ibond,keyword,lstbnd(1,nbonds),lstbnd(2,nbonds),(prmbnd(ja,nbonds),ja=1,mxpbnd)
                       Else
  Write(nrite,"(12x,i10,a8,2i10,2x,10f15.6)") &
                  ibond,keyword,lstbnd(1,nbonds),lstbnd(2,nbonds),(prmbnd(ja,nbonds),ja=1,mxpbnd)
                       End If
                    End If

! test for mistyped chemical bond unit

                    If (iatm1 == iatm2) Call error(33)

! convert energy units to internal units

                    prmbnd(1,nbonds)=prmbnd(1,nbonds)*engunit

                    If (Abs(keybnd(nbonds)) == 3) Then
                       prmbnd(2,nbonds)=prmbnd(2,nbonds)*engunit
                    End If

                    If (Abs(keybnd(nbonds)) == 6) Then
                       prmbnd(3,nbonds)=prmbnd(3,nbonds)*engunit
                       prmbnd(4,nbonds)=prmbnd(4,nbonds)*engunit
                    End If

                    If (Abs(keybnd(nbonds)) == 7) Then
                       prmbnd(3,nbonds)=prmbnd(3,nbonds)*engunit
                    End If

                 End Do

! Check for multiple chemical bond entries

                 Do i=nbonds-numbonds(itmols)+1,nbonds
                    is(1)=Min(lstbnd(1,i),lstbnd(2,i))
                    is(2)=Max(lstbnd(1,i),lstbnd(2,i))

                    Do j=i+1,nbonds
                       js(1)=Min(lstbnd(1,j),lstbnd(2,j))
                       js(2)=Max(lstbnd(1,j),lstbnd(2,j))

                       If (js(1) == is(1) .and. js(2) == is(2)) Then
                          If (l_str) Call warning(420,Real(i,wp),Real(j,wp),0.0_wp)
!                          Call error(620)
                       End If
                    End Do
                 End Do

! read intramolecular angular potential parameters

              Else If (word(1:6) == 'angles') Then

                 If (.not.l_ang) Call error(210)
                 l_ang=.false.

                 lbook=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numang(itmols)=numang(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of bond angles',10x,i10)") ntmp
  Write(nrite,"(/,1x,'bond angle details:', &
       & /,/,18x,'unit',5x,'key',5x,'index',5x,'index',5x,'index',7x,'f-const',8x,'angle',/)")
                 End If

                 Do iang=1,numang(itmols)
                    nangle=nangle+1
                    If (nangle > mxtang) Call error(50)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of bond angle

                    Call lower_case(word)
                    keyword=word(1:4)

                    If (keyword(1:1) /= '-') lexcl=.true.

                    If      (keyword == 'harm') Then
                       keyang(nangle)=1
                    Else If (keyword == '-hrm') Then
                       keyang(nangle)=-1
                    Else If (keyword == 'quar') Then
                       keyang(nangle)=2
                    Else If (keyword == '-qur') Then
                       keyang(nangle)=-2
                    Else If (keyword == 'thrm') Then
                       keyang(nangle)=3
                    Else If (keyword == '-thm') Then
                       keyang(nangle)=-3
                    Else If (keyword == 'shrm') Then
                       keyang(nangle)=4
                    Else If (keyword == '-shm') Then
                       keyang(nangle)=-4
                    Else If (keyword == 'bvs1') Then
                       keyang(nangle)=5
                    Else If (keyword == '-bv1') Then
                       keyang(nangle)=-5
                    Else If (keyword == 'bvs2') Then
                       keyang(nangle)=6
                    Else If (keyword == '-bv2') Then
                       keyang(nangle)=-6
                    Else If (keyword == 'hcos') Then
                       keyang(nangle)=7
                    Else If (keyword == '-hcs') Then
                       keyang(nangle)=-7
                    Else If (keyword == 'cos' ) Then
                       keyang(nangle)=8
                    Else If (keyword == '-cos') Then
                       keyang(nangle)=-8
                    Else If (keyword == 'mmsb') Then
                       keyang(nangle)=9
                    Else If (keyword == '-msb') Then
                       keyang(nangle)=-9
                    Else If (keyword == 'stst') Then
                       keyang(nangle)=10
                    Else If (keyword == '-sts') Then
                       keyang(nangle)=-10
                    Else If (keyword == 'stbe') Then
                       keyang(nangle)=11
                    Else If (keyword == '-stb') Then
                       keyang(nangle)=-11
                    Else If (keyword == 'cmps') Then
                       keyang(nangle)=12
                    Else If (keyword == '-cmp') Then
                       keyang(nangle)=-12
                    Else

                       If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
                       Call error(440)

                    End If

! read angle atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm3=Nint(word_2_real(word))

                    lstang(1,nangle)=iatm1
                    lstang(2,nangle)=iatm2
                    lstang(3,nangle)=iatm3

                    Call get_word(record,word)
                    prmang(1,nangle)=word_2_real(word)
                    Call get_word(record,word)
                    prmang(2,nangle)=word_2_real(word)
                    Call get_word(record,word)
                    prmang(3,nangle)=word_2_real(word)
                    Call get_word(record,word)
                    prmang(4,nangle)=word_2_real(word)

! test for frozen atom pairs

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2
                    isite3 = nsite - numsit(itmols) + iatm3

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3) /= 0) Then
  Write(nrite,"(4x,a8,i10,a8,3i10,10f15.6)") &
       '*frozen*',iang,keyword,(lstang(ia,nangle),ia=1,3),(prmang(ja,nangle),ja=1,mxpang)
                       Else
  Write(nrite,"(12x,i10,a8,3i10,10f15.6)") &
                  iang,keyword,(lstang(ia,nangle),ia=1,3),(prmang(ja,nangle),ja=1,mxpang)
                       End If
                    End If

! test for mistyped bond angle unit

                    If (iatm1 == iatm2 .or. iatm1 == iatm3 .or. iatm2 == iatm3) Call error(66)

! convert energies to internal units

                    prmang(1,nangle) = prmang(1,nangle)*engunit

                    If      (Abs(keyang(nangle)) == 2) Then
                       prmang(3,nangle) = prmang(3,nangle)*engunit
                       prmang(4,nangle) = prmang(4,nangle)*engunit
                    Else If (Abs(keyang(nangle)) == 12) Then
                       prmang(2,nangle) = prmang(2,nangle)*engunit
                       prmang(3,nangle) = prmang(3,nangle)*engunit
                    End If

! convert angles to radians

                    If      (Abs(keyang(nangle)) == 12) Then
                       prmang(4,nangle) = prmang(4,nangle)*(pi/180.0_wp)
                    Else If (Abs(keyang(nangle)) /= 10) Then
                       prmang(2,nangle) = prmang(2,nangle)*(pi/180.0_wp)
                    End If

                 End Do

! Check for multiple bond angle entries

                 Do i=nangle-numang(itmols)+1,nangle
                    is(1)=Min(lstang(1,i),lstang(3,i))
                    is(2)=lstang(2,i)
                    is(3)=Max(lstang(1,i),lstang(3,i))

                    Do j=i+1,nangle
                       js(1)=Min(lstang(1,j),lstang(3,j))
                       js(2)=lstang(2,j)
                       js(3)=Max(lstang(1,j),lstang(3,j))

                       If (js(1) == is(1) .and. js(2) == is(2) .and. &
                           js(3) == is(3)) Then
                          If (l_str) Call warning(430,Real(i,wp),Real(j,wp),0.0_wp)
!                          Call error(620)
                       End If
                    End Do
                 End Do

! read intramolecular dihedral potential parameters

              Else If (word(1:6) == 'dihedr') Then

                 If (.not.l_dih) Call error(220)
                 l_dih=.false.

                 lbook=.true.
                 lexcl=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numdih(itmols)=numdih(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of dihedral angles',6x,i10)") ntmp
  Write(nrite,"(/,1x,'dihedral angle details:', &
       & /,/,18x,'unit',5x,'key',5x,'index',5x,'index',5x,'index',5x,'index', &
       & 7x,'f-const',8x,'angle',9x,'trig',11x,'1-4 elec',7x,'1-4 vdw',/)")
                 End If

                 Do idih=1,numdih(itmols)
                    ndihed=ndihed+1
                    If (ndihed > mxtdih) Call error(60)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of dihedral bond angle

                    Call lower_case(word)
                    keyword=word(1:4)

                    If      (keyword == 'cos' ) Then
                       keydih(ndihed)=1
                    Else If (keyword == 'harm') Then
                       keydih(ndihed)=2
                    Else If (keyword == 'hcos') Then
                       keydih(ndihed)=3
                    Else If (keyword == 'cos3') Then
                       keydih(ndihed)=4
                    Else If (keyword == 'ryck') Then
                       keydih(ndihed)=5
                    Else If (keyword == 'rbf' ) Then
                       keydih(ndihed)=6
                    Else If (keyword == 'opls') Then
                       keydih(ndihed)=7
                    Else

                       If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
                       Call error(448)

                    End If

! read dihedral atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm3=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm4=Nint(word_2_real(word))

                    lstdih(1,ndihed)=iatm1
                    lstdih(2,ndihed)=iatm2
                    lstdih(3,ndihed)=iatm3
                    lstdih(4,ndihed)=iatm4

                    Call get_word(record,word)
                    prmdih(1,ndihed)=word_2_real(word)
                    Call get_word(record,word)
                    prmdih(2,ndihed)=word_2_real(word)
                    Call get_word(record,word)
                    prmdih(3,ndihed)=word_2_real(word)
                    Call get_word(record,word)
                    prmdih(4,ndihed)=word_2_real(word)
                    Call get_word(record,word)
                    prmdih(5,ndihed)=word_2_real(word)
                    Call get_word(record,word)
                    prmdih(6,ndihed)=word_2_real(word)
                    Call get_word(record,word)
                    prmdih(7,ndihed)=word_2_real(word)

! test for frozen atom pairs

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2
                    isite3 = nsite - numsit(itmols) + iatm3
                    isite4 = nsite - numsit(itmols) + iatm4

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) /= 0) Then
  Write(nrite,"(4x,a8,i10,a8,4i10,10f15.6)") &
       '*frozen*',idih,keyword,(lstdih(ia,ndihed),ia=1,4),(prmdih(ja,ndihed),ja=1,mxpdih)
                       Else
  Write(nrite,"(12x,i10,a8,4i10,10f15.6)") &
                  idih,keyword,(lstdih(ia,ndihed),ia=1,4),(prmdih(ja,ndihed),ja=1,mxpdih)
                       End If
                    End If

! test for mistyped dihedral unit

                    If ( iatm1 == iatm2 .or. iatm1 == iatm3 .or. &
                         iatm2 == iatm3 .or. iatm1 == iatm4 .or. &
                         iatm2 == iatm4 .or. iatm3 == iatm4 ) Call error(67)

! convert energies to internal units and angles to radians

                    prmdih(1,ndihed)=prmdih(1,ndihed)*engunit

                    If (keydih(ndihed) == 4) Then
                       prmdih(2,ndihed)=prmdih(2,ndihed)*engunit
                       prmdih(3,ndihed)=prmdih(3,ndihed)*engunit
                    Else If (keydih(ndihed) == 7) Then
                       prmdih(2,ndihed)=prmdih(2,ndihed)*engunit
                       prmdih(3,ndihed)=prmdih(3,ndihed)*engunit
                       prmdih(6,ndihed)=prmdih(6,ndihed)*engunit
                       prmdih(7,ndihed)=prmdih(7,ndihed)*(pi/180.0_wp)
                    Else
                       prmdih(2,ndihed)=prmdih(2,ndihed)*(pi/180.0_wp)
                    End If
                 End Do

! Check for multiple dihedral angle entries

                 Do i=ndihed-numdih(itmols)+1,ndihed
                    is(1)=lstdih(1,i)
                    is(2)=lstdih(2,i)
                    is(3)=lstdih(3,i)
                    is(4)=lstdih(4,i)

                    Do j=i+1,ndihed
                       js(1)=lstdih(1,j)
                       js(2)=lstdih(2,j)
                       js(3)=lstdih(3,j)
                       js(4)=lstdih(4,j)

                       If ((js(1) == is(1) .and. js(2) == is(2) .and. &
                            js(3) == is(3) .and. js(4) == is(4)) .or. &
                           (js(1) == is(4) .and. js(2) == is(3) .and. &
                            js(3) == is(2) .and. js(4) == is(1))) Then
                          If (l_str) Call warning(440,Real(i,wp),Real(j,wp),0.0_wp)
!                          Call error(620)
                       End If
                    End Do
                 End Do

! read intramolecular inversion potential parameters

              Else If (word(1:6) == 'invers') Then

                 If (.not.l_inv) Call error(230)
                 l_inv=.false.

                 lbook=.true.
                 lexcl=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 numinv(itmols)=numinv(itmols)+ntmp

                 If (idnode == 0 .and. l_top) Then
  Write(nrite,"(/,/,1x,'number of inversion angles',5x,i10)") ntmp
  Write(nrite,"(/,1x,'inversion angle details:', &
       & /,/,18x,'unit',5x,'key',5x,'index',5x,'index',5x,'index',5x,'index', &
       & 7x,'f-const',8x,'angle',8x,'factor',/)")
                 End If

                 Do iinv=1,numinv(itmols)
                    ninver=ninver+1
                    If (ninver > mxtinv) Call error(73)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of inversion potential

                    Call lower_case(word)
                    keyword=word(1:4)

                    If      (keyword == 'harm') Then
                       keyinv(ninver)=1
                    Else If (keyword == 'hcos') Then
                       keyinv(ninver)=2
                    Else If (keyword == 'plan') Then
                       keyinv(ninver)=3
                    Else If (keyword == 'xpln') Then
                       keyinv(ninver)=4
                    Else If (keyword == 'calc') Then
                       keyinv(ninver)=5
                    Else

                       If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
                       Call error(449)

                    End If

! read inversion atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm3=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm4=Nint(word_2_real(word))

                    lstinv(1,ninver)=iatm1
                    lstinv(2,ninver)=iatm2
                    lstinv(3,ninver)=iatm3
                    lstinv(4,ninver)=iatm4

                    Call get_word(record,word)
                    prminv(1,ninver)=word_2_real(word)
                    Call get_word(record,word)
                    prminv(2,ninver)=word_2_real(word)
                    Call get_word(record,word)
                    prminv(3,ninver)=word_2_real(word)

! test for frozen atom pairs

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2
                    isite3 = nsite - numsit(itmols) + iatm3
                    isite4 = nsite - numsit(itmols) + iatm4

                    If (idnode == 0 .and. l_top) Then
                       If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) /= 0) Then
  Write(nrite,"(4x,a8,i10,a8,4i10,10f15.6)") &
       '*frozen*',iinv,keyword,(lstinv(ia,ninver),ia=1,4),(prminv(ja,ninver),ja=1,mxpinv)
                       Else
  Write(nrite,"(12x,i10,a8,4i10,10f15.6)") &
                  iinv,keyword,(lstinv(ia,ninver),ia=1,4),(prminv(ja,ninver),ja=1,mxpinv)
                       End If
                    End If

! test for mistyped inversion unit

                    If ( iatm1 == iatm2 .or. iatm1 == iatm3 .or. &
                         iatm2 == iatm3 .or. iatm1 == iatm4 .or. &
                         iatm2 == iatm4 .or. iatm3 == iatm4 ) Call error(68)

! convert energies to internal units and angles to radians

                    prminv(1,ninver)=prminv(1,ninver)*engunit

                    If (keyinv(ninver) == 5) Then
                       prminv(2,ninver)=prminv(2,ninver)*engunit
                    Else
                       prminv(2,ninver)=prminv(2,ninver)*(pi/180.0_wp)
                       If (keyinv(ninver) == 2) &
                          prminv(2,ninver)=Cos(prminv(2,ninver))
                    End If
                 End Do

! Check for multiple inversion angle entries

                 Do i=ninver-numinv(itmols)+1,ninver
                    is(1)=lstinv(1,i)
                    is(2)=lstinv(2,i)
                    is(3)=lstinv(3,i)
                    is(4)=lstinv(4,i)
                    Call shellsort(3,is(2:4))

                    Do j=i+1,ninver
                       js(1)=lstinv(1,j)
                       js(2)=lstinv(2,j)
                       js(3)=lstinv(3,j)
                       js(4)=lstinv(4,j)
                       Call shellsort(3,js(2:4))

                       If (js(1) == is(1) .and. js(2) == is(2) .and. &
                           js(3) == is(3) .and. js(4) == is(4)) Then
                          If (l_str) Call warning(450,Real(i,wp),Real(j,wp),0.0_wp)
!                          Call error(620)
                       End If
                    End Do
                 End Do

! finish of data for one molecular type

              Else If (word(1:6) == 'finish') Then

! running totals of number of atoms and frozen atoms, and general types of
! intra-like interactions in system

                 megatm=megatm+nummols(itmols)*numsit(itmols)
                 megfrz=megfrz+nummols(itmols)*numfrz(itmols)

                 megshl=megshl+nummols(itmols)*numshl(itmols)

                 megcon=megcon+nummols(itmols)*(numcon(itmols)-frzcon)
                 megpmf=megpmf+nummols(itmols)*numpmf(itmols)

                 megrgd=megrgd+nummols(itmols)*(numrgd(itmols)-frzrgd)

                 megtet=megtet+nummols(itmols)*numteth(itmols)

                 megbnd=megbnd+nummols(itmols)*numbonds(itmols)
                 megang=megang+nummols(itmols)*numang(itmols)
                 megdih=megdih+nummols(itmols)*numdih(itmols)
                 meginv=meginv+nummols(itmols)*numinv(itmols)

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

! report total molecules and sites

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'total number of molecules',6x,i10)") Sum(nummols(1:ntpmls))
  Write(nrite,"(/,1x,'total number of sites',10x,i10)") nsite
        End If

! just finished with molecular data

! Initialise number of free (of RB structures)
! and free frozen atoms/particles

        atmfre=megatm
        atmfrz=megfrz

! test shells masses and define model

        If (megshl > 0) Then
           If (.not.lshl_one) Then
              keyshl = 1
  If (idnode == 0) Write(nrite,"(/,/,1x,'adiabatic shell model in operation')")
           Else
              If (lshl_all) Then
                 keyshl = 2
  If (idnode == 0) Write(nrite,"(/,/,1x,'relaxed shell model in operation')")
              Else
                 Call error(476)
              End If
           End If
        End If

! Check for charges in the system

        If (keyfce /= 0 .and. Abs(sumchg) <= zero_plus) Then
           If (idnode == 0) Call warning(4,sumchg,0.0_wp,0.0_wp)
           If (l_str) Then
              keyfce=0
              If (idnode == 0) Write(nrite,"(1x,'Electrostatics switched off!!!')")
           End If
        End If

! calculate total system charge

        sumchg=0.0_wp
        jsite=0
        Do itmols=1,ntpmls
           Do msite=1,numsit(itmols)
              jsite=jsite+1
              sumchg=sumchg+Real(nummols(itmols),wp)*chgsit(jsite)
           End Do
        End Do

        If (Abs(sumchg) > 1.0e-6_wp .and. idnode == 0) Call warning(5,sumchg,0.0_wp,0.0_wp)

! test for constrained, RBed and tethered shells in core-shell units

        If (megshl > 0) Then
           nsite =0
           nshels=0
           nconst=0
           nrigid=0
           nteth =0
           Do itmols=1,ntpmls
              Do ishls=1,numshl(itmols)
                 nshels=nshels+1
                 iatm1=lstshl(2,nshels)

                 dofsit(nsite+iatm1)=-3.0_wp ! shells have no DoF

                 Do icnst=1,numcon(itmols)
                    nconst=nconst+1
                    iatm2=lstcon(1,nconst)
                    iatm3=lstcon(2,nconst)

                    If (iatm2 == iatm1 .or. iatm3 == iatm1) Then
                       Call warning(301,Real(ishls,wp),Real(icnst,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nconst=nconst-numcon(itmols)

                 Do irgd=1,numrgd(itmols)
                    nrigid=nrigid+1

                    lrgd=lstrgd(0,nrigid)
                    If (Any(lstrgd(1:lrgd,nrigid) == iatm1)) Then
                       Call warning(302,Real(ishls,wp),Real(irgd,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nrigid=nrigid-numrgd(itmols)

                 Do iteth=1,numteth(itmols)
                    nteth=nteth+1
                    iatm4=lsttet(nteth)

                    If (iatm4 == iatm1) Then
                       Call warning(303,Real(ishls,wp),Real(iteth,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nteth=nteth-numteth(itmols)
              End Do

              nsite=nsite+numsit(itmols)
           End Do
        End If

! RB particulars

        If (m_rgd > 0) Then

! test for constraint units on RB units

           If (m_con > 0) Then
              nsite =0
              nconst=0
              nrigid=0
              Do itmols=1,ntpmls
                 Do icnst=1,numcon(itmols)
                    nconst=nconst+1
                    iatm1=lstcon(1,nconst)
                    iatm2=lstcon(2,nconst)

                    Do irgd=1,numrgd(itmols)
                       nrigid=nrigid+1

                       lrgd=lstrgd(0,nrigid)
                       If (Any(lstrgd(1:lrgd,nrigid) == iatm1) .and. Any(lstrgd(1:lrgd,nrigid) == iatm2)) Then
                          Call warning(304,Real(icnst,wp),Real(irgd,wp),Real(itmols,wp))

                          Call error(97)
                       End If
                    End Do
                 End Do
                 nsite=nsite+numsit(itmols)
              End Do
           End If

! test for PMF units on RB units

           If (megpmf > 0) Then
              nsite =0
              nrigid=0
              Do itmols=1,ntpmls
                 Do i=1,numpmf(itmols)
                    Do ipmf=1,2
                       Do jpmf=1,mxtpmf(ipmf)
                          iatm1=lstpmf(jpmf,ipmf)
                          isite1=nsite+iatm1

                          Do irgd=1,numrgd(itmols)
                             nrigid=nrigid+1

                             lrgd=lstrgd(0,nrigid)
                             If (Any(lstrgd(1:lrgd,nrigid) == iatm1)) Then
                                Call warning(295,Real(ipmf,wp),Real(irgd,wp),Real(itmols,wp))

                                Call error(93)
                             End If
                          End Do
                       End Do
                    End Do
                 End Do
                 nsite=nsite+numsit(itmols)
              End Do
           End If

! Index RBs' sites (fresit=1), correct atmfre & atmfrz
! and test for unfrozen weightless members of a RB unit type
! (partly frozen RB but with unfrozen members being weightless)
! and correct frzsit,dofsit,rgdwg1,frzrgd,megfrz,megrgd if needed

           nsite =0
           nrigid=0
           Do itmols=1,ntpmls
              ntmp=0
              ntab=0


              ifrz=0

              frzrgd=0

              Do irgd=1,numrgd(itmols)
                 nrigid=nrigid+1

                 lrgd=lstrgd(0,nrigid)

                 ntmp=ntmp+lrgd

                 krgd=0
                 If (rgdwgt(0,nrigid) < 1.0e-6_wp .and. rgdfrz(0,nrigid) < lrgd) Then
                    krgd=1

                    rgdfrz(0,nrigid)=lrgd
                    rgdwg1(0,nrigid)=Real(lrgd,wp)

                    Call warning(305,Real(irgd,wp),Real(itmols,wp),0.0_wp)

                    frzrgd=frzrgd+1
                 End If

                 Do jrgd=1,lrgd
                    iatm1=lstrgd(jrgd,nrigid)
                    isite1=nsite+iatm1

                    fresit(isite1)=1

                    If (frzsit(isite1) == 1) Then
                       ntab=ntab+1
                    Else
                       If (krgd == 1) Then
                          ifrz=ifrz+1

                          frzsit(isite1)=1
                          dofsit(isite1)=0.0_wp

                          rgdfrz(jrgd,nrigid)=1
                       End If
                    End If
                 End Do
              End Do

              atmfre=atmfre-ntmp  *nummols(itmols)
              atmfrz=atmfrz-ntab  *nummols(itmols)

              megfrz=megfrz+ifrz  *nummols(itmols)

              megrgd=megrgd-frzrgd*nummols(itmols)

              nsite=nsite+numsit(itmols)
           End Do

! Globalise megrgd,imcon,rcut for RBs

           rgdmeg=megrgd
           rgdimc=imcon
           rgdrct=rcut

        End If

! read in rdf pairs

     Else If (word(1:3) == 'rdf') Then

        Call get_word(record,word)
        ntprdf=Nint(word_2_real(word))

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'number of specified rdf look up pairs    ',i10)") ntprdf
           If (l_top) &
  Write(nrite,"(/,10x,'pair',2x,'atom 1',2x,'atom 2',/)")
        End If

        If (ntprdf > mxrdf) Call error(107)

        Do itprdf=1,ntprdf

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom1=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)

           If (idnode == 0 .and. l_top) &
  Write(nrite,"(4x,i10,2a8)") itprdf,atom1,atom2

           katom1=0
           katom2=0

           Do jtpatm=1,ntpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
           End Do

           If (katom1 == 0 .or. katom2 == 0) Call error(108)

           ka1=Max(katom1,katom2)
           ka2=Min(katom1,katom2)

           keyrdf=(ka1*(ka1-1))/2+ka2

           If (keyrdf > mxrdf) Call error(109)

           If (lstrdf(keyrdf) /= 0) Call error(110)

           lstrdf(keyrdf)=itprdf

        End Do

! read in the vdw potential energy parameters

     Else If (word(1:3) == 'vdw') Then

        Call get_word(record,word)
        ntpvdw=Nint(word_2_real(word))
        Call get_word(record,word)
        ltable=(word(1:5) == 'table')

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'number of specified vdw potentials       ',i10)") ntpvdw
           If (l_top) &
  Write(nrite,"(/,7x,'pair',5x,'atom 1',2x,'atom 2',5x,'key',30x,'parameters',/)")
        End If

        If (ntpvdw > mxvdw) Call error(80)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

        Do itpvdw=1,ntpvdw

           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom1=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)

           Call get_word(record,word)
           Call lower_case(word)
           keyword=word(1:4)

           If      (keyword == 'tab' ) Then
              keypot=0
           Else If (keyword == '12-6') Then
              keypot=1
           Else If (keyword == 'lj'  ) Then
              keypot=2
           Else If (keyword == 'nm'  ) Then
              keypot=3
           Else If (keyword == 'buck') Then
              keypot=4
           Else If (keyword == 'bhm' ) Then
              keypot=5
           Else If (keyword == 'hbnd') Then
              keypot=6
           Else If (keyword == 'snm' ) Then
              keypot=7
           Else If (keyword == 'hcnm') Then
              keypot=7
           Else If (keyword == 'mors') Then
              keypot=8
           Else If (keyword == 'wca' ) Then
              keypot=9
           Else

              If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
              Call error(452)

           End If

           ltable=(ltable .or. (keypot == 0))

           If (keypot == 0) Then
              If (idnode == 0 .and. l_top) &
  Write(nrite,"(1x,i10,5x,2a8,30x,a9)") itpvdw,atom1,atom2,"tabulated"
           Else
              Call get_word(record,word)
              parpot(1)=word_2_real(word)
              Call get_word(record,word)
              parpot(2)=word_2_real(word)
              Call get_word(record,word)
              parpot(3)=word_2_real(word)
              Call get_word(record,word)
              parpot(4)=word_2_real(word)
              Call get_word(record,word)
              parpot(5)=word_2_real(word)

              If (idnode == 0 .and. l_top) &
  Write(nrite,"(1x,i10,5x,2a8,3x,a4,1x,10f20.6)") itpvdw,atom1,atom2,keyword,(parpot(j),j=1,mxpvdw)

! convert energies to internal unit

              parpot(1) = parpot(1)*engunit

              If (keypot == 1) Then
                 parpot(2)=parpot(2)*engunit
              Else If (keypot == 4) Then
                 parpot(3)=parpot(3)*engunit
              Else If (keypot == 5) Then
                 parpot(4)=parpot(4)*engunit
                 parpot(5)=parpot(5)*engunit
              Else If (keypot == 6) Then
                 parpot(2)=parpot(2)*engunit
              Else If (keypot == 9) Then
                 parpot(2)=Abs(parpot(2))
                 If (parpot(3) > parpot(2)/2.0_wp) &
                    parpot(3)=Sign(1.0_wp,parpot(3))*parpot(2)/2.0_wp
                 parpot(4)=2.0_wp**(1.0_wp/6.0_wp)*parpot(2)+parpot(3)
              End If
           End If

           katom1=0
           katom2=0

           Do jtpatm=1,ntpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
           End Do

           If (katom1 == 0 .or. katom2 == 0) Call error(81)

           ka1=Max(katom1,katom2)
           ka2=Min(katom1,katom2)

           keyvdw=(ka1*(ka1-1))/2+ka2

           If (keyvdw > mxvdw) Call error(82)

           If (lstvdw(keyvdw) /= 0) Call error(15)

           lstvdw(keyvdw)=itpvdw
           ltpvdw(itpvdw)=keypot

           Do i=1,mxpvdw
              prmvdw(i,itpvdw)=parpot(i)
           End Do
        End Do

        If (ntpvdw > 0) Then

! test for unspecified atom-atom potentials

           ntab =(ntpatm*(ntpatm+1))/2
           If (ntpvdw < ntab) Then
              Call warning(120,0.0_wp,0.0_wp,0.0_wp)

              If (ntpvdw > mxvdw) Call error(80)

! put undefined potentials outside range

              Do i=1,ntab
                 If (lstvdw(i) == 0) lstvdw(i) = ntpvdw+1
              End Do

              Do i=ntpvdw+1,ntab
                 ltpvdw(i) = -1
              End Do
           End If

! generate nonbonded force arrays

           If (.not.l_n_v) Then
              Call vdw_generate(ltable,rvdw)

              If (ltable) Call vdw_table_read(rvdw)
           End If

        End If

! read in the metal potential energy parameters

     Else If (word(1:3) == 'met') Then

        Call get_word(record,word)
        ntpmet=Nint(word_2_real(word))

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'number of specified metal potentials     ',i10)") ntpmet
           If (l_top) &
  Write(nrite,"(/,7x,'pair',5x,'atom 1',2x,'atom 2',5x,'key',30x,'parameters',/)")
        End If

        If (ntpmet > mxmet) Call error(71)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

        lmet_safe=.true.
        Do itpmet=1,ntpmet

           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom1=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)

           Call get_word(record,word)
           Call lower_case(word)
           keyword=word(1:4)

           If      (keyword(1:3) == 'eam' ) Then
              keypot=0
           Else If (keyword(1:4) == 'fnsc') Then
              keypot=1
           Else If (keyword(1:4) == 'exfs') Then
              keypot=2
           Else If (keyword(1:4) == 'stch') Then
              keypot=3
           Else If (keyword(1:4) == 'gupt') Then
              keypot=4
           Else

              If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
              Call error(461)

           End If

           If (keypot > 0) Then
              Call get_word(record,word)
              parpot(1)=word_2_real(word)
              Call get_word(record,word)
              parpot(2)=word_2_real(word)
              Call get_word(record,word)
              parpot(3)=word_2_real(word)
              Call get_word(record,word)
              parpot(4)=word_2_real(word)
              Call get_word(record,word)
              parpot(5)=word_2_real(word)
              Call get_word(record,word)
              parpot(6)=word_2_real(word)
              Call get_word(record,word)
              parpot(7)=word_2_real(word)
              Call get_word(record,word)
              parpot(8)=word_2_real(word)
              Call get_word(record,word)
              parpot(9)=word_2_real(word)

              If (idnode == 0 .and. l_top) &
  Write(nrite,"(1x,i10,5x,2a8,3x,a4,1x,10f15.6)") itpmet,atom1,atom2,keyword,(parpot(j),j=1,mxpmet)
           End If

           katom1=0
           katom2=0

           Do jtpatm=1,ntpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
           End Do

           If (katom1 == 0 .or. katom2 == 0) Call error(81)

           ka1=Max(katom1,katom2)
           ka2=Min(katom1,katom2)

           keymet=(ka1*(ka1-1))/2+ka2

           If (keymet > mxmet) Call error(82)
           If (lstmet(keymet) /= 0) Call error(141)

           lstmet(keymet)=itpmet
           ltpmet(itpmet)=keypot

           If (itpmet > 1) Then
              If (keypot /= ltpmet(itpmet-1)) lmet_safe=.false.
           End If

! convert energies to internal unit

           If (keypot > 0) Then
              parpot(1)=parpot(1)*engunit

              If      (keypot == 1) Then
                 parpot(2)=parpot(2)*engunit
                 parpot(3)=parpot(3)*engunit
                 parpot(5)=parpot(5)*engunit
              Else If (keypot == 2) Then
                 parpot(2)=parpot(2)*engunit
                 parpot(3)=parpot(3)*engunit
                 parpot(4)=parpot(4)*engunit
                 parpot(5)=parpot(5)*engunit
                 parpot(7)=parpot(7)*engunit
              Else If (keypot == 4) Then
                 parpot(4)=parpot(4)*engunit
              End If

              Do i=1,mxpmet
                 prmmet(i,itpmet)=parpot(i)
              End Do
           End If

        End Do

        If (ntpmet /= 0) Then

! test metal potentials mix-up

           If (.not.lmet_safe) Call error(92)

! test for unspecified atom-atom potentials

           ntab=(ntpatm*(ntpatm+1))/2
           If (ntpmet < ntab) Then
              Call warning(120,0.0_wp,0.0_wp,0.0_wp)

              If (ntpmet > mxmet) Call error(80)

! put undefined potentials outside range

              Do i=1,ntab
                 If (lstmet(i) == 0) lstmet(i)=ntpmet+1
              End Do

              Do i=ntpmet+1,ntab
                 ltpmet(i) = -1
              End Do
           End If

! generate metal force arrays

           If (keypot == 0) Then
              Call metal_table_read()
           Else
              Call metal_generate(rmet)
           End If

        End If

! read in the tersoff potential energy parameters

     Else If (word(1:7) == 'tersoff') Then

        Call get_word(record,word)
        ntpter=Nint(word_2_real(word))

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'number of specified tersoff potentials   ',i10)") ntpter
           If (l_top) &
  Write(nrite,"(/,5x,'number',5x,'atom',7x,'key',37x,'parameters'/)")
        End If

        If (ntpter > mxter) Call error(72)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

        Do itpter=1,ntpter
           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom0=word(1:8)

           Call get_word(record,word)
           Call lower_case(word)
           keyword=word(1:4)

           If      (keyword == 'ters') Then
              keypot=1
           Else

              If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
              Call error(432)

           End If

           Call get_word(record,word)
           parpot(1 )=word_2_real(word)      ! A_i
           Call get_word(record,word)
           parpot(2 )=word_2_real(word)      ! a_i
           Call get_word(record,word)
           parpot(3 )=word_2_real(word)      ! B_i
           Call get_word(record,word)
           parpot(4 )=word_2_real(word)      ! b_i
           Call get_word(record,word)
           parpot(5 )=word_2_real(word)      ! R_i

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           parpot(6 )=word_2_real(word)      ! S_i
           Call get_word(record,word)
           parpot(7 )=word_2_real(word)      ! beta_i
           Call get_word(record,word)
           parpot(8 )=word_2_real(word)      ! eta_i
           Call get_word(record,word)
           parpot(9 )=word_2_real(word)      ! c_i
           Call get_word(record,word)
           parpot(10)=word_2_real(word)      ! d_i
           Call get_word(record,word)
           parpot(11)=word_2_real(word)      ! h_i

           If (idnode == 0 .and. l_top) Then
  Write(nrite,"(1x,i10,5x,a8,3x,a4,1x,10(1p,e13.4))") itpter,atom0,keyword,(parpot(j),j=1,5)
  Write(nrite,"(32x,10(1p,e13.4))") (parpot(j),j=6,mxpter)
           End If

           katom0=0

           Do jtpatm=1,ntpatm
              If (atom0 == unqatm(jtpatm)) katom0=jtpatm
           End Do

           If (katom0 == 0) Call error(74)

! convert parameters to internal units

           If (keypot == 1) Then
              parpot(1)=parpot(1)*engunit
              parpot(3)=parpot(3)*engunit
           End If

           If (lstter(katom0) > 0) Call error(76)

           lfrter(katom0)=.true.

           lstter(katom0)=itpter
           ltpter(itpter)=keypot

! calculate max tersoff cutoff

           rcter=Max(rcter,parpot(6))

! store tersoff single atom potential parameters

           Do i=1,mxpter
              prmter(i,itpter)=parpot(i)
           End Do

        End Do

        If (ntpter > 0) Then
           If (rcter < 1.0e-6_wp) Call error(79)
           If (rcut < 2.0_wp*rcter) Call error(102)
        End If

! start processing cross atom potential parameters

        If (idnode == 0) Then
  Write(nrite,"(/,1x,'number of tersoff cross terms            ',i10)") (ntpter*(ntpter+1))/2
           If (l_top) &
  Write(nrite,"(/,7x,'pair',5x,'atom 1',2x,'atom 2',14x,'parameters',/)")
        End If

        Do icross=1,(ntpter*(ntpter+1))/2

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom1=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)

           Call get_word(record,word)
           parpot(1)=word_2_real(word,1.0_wp)                                       ! chi_ij
           Call get_word(record,word)
           parpot(2)=word_2_real(word,1.0_wp)                                       ! omega_ij
           Call get_word(record,word)
           parpot(3)=Merge(1.0_wp,0.0_wp,Abs(word_2_real(word,0.0_wp)) > zero_plus) ! delta_ij

           katom1=0
           katom2=0

           Do jtpatm=1,ntpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
           End Do

           If (katom1 == 0 .or. katom2 == 0) Call error(74)

           ka1=Max(lstter(katom1),lstter(katom2))
           ka2=Min(lstter(katom1),lstter(katom2))

           keyter=(ka1*(ka1-1))/2+ka2

           prmter2(keyter,1)=parpot(1)
           prmter2(keyter,2)=parpot(2)
           prmter2(keyter,3)=parpot(3)

           If (idnode == 0 .and. l_top) &
  Write(nrite,"(1x,i10,5x,2a8,1x,2f15.6,i6)") icross,atom1,atom2,(parpot(j),j=1,2),Nint(parpot(3))

        End Do

! generate tersoff interpolation arrays

        If (ntpter > 0) Call tersoff_generate(rcter)

! read in the three-body potential energy parameters

     Else If (word(1:3) == 'tbp') Then

        Call get_word(record,word)
        ntptbp=Nint(word_2_real(word))

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'number of specified three-body potentials',i10)") ntptbp
           If (l_top) &
  Write(nrite,"(/,4x,'triplet',5x,'atom 1',2x,'atom 2',2x,'atom 3',5x,'key',22x,'parameters',/)")
        End If

        If (ntptbp > mxtbp) Call error(83)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

! Initialise head of group atom to not participating in an interaction

        Do itbp=1,mxtbp,mx2tbp
           lsttbp(itbp)=-1
        End Do

        Do itptbp=1,ntptbp
           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

! Note the order!! atom0 is the central atom

           atom1=word(1:8)
           Call get_word(record,word)
           atom0=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)

           Call get_word(record,word)
           Call lower_case(word)
           keyword=word(1:4)

           If      (keyword == 'harm') Then
              keypot=1
           Else If (keyword == 'thrm') Then
              keypot=2
           Else If (keyword == 'shrm') Then
              keypot=3
           Else If (keyword == 'bvs1') Then
              keypot=4
           Else If (keyword == 'bvs2') Then
              keypot=5
           Else If (keyword == 'hbnd') Then
              keypot=6
           Else

              If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
              Call error(442)

           End If

           Call get_word(record,word)
           parpot(1)=word_2_real(word)
           Call get_word(record,word)
           parpot(2)=word_2_real(word)
           Call get_word(record,word)
           parpot(3)=word_2_real(word)
           Call get_word(record,word)
           parpot(4)=word_2_real(word)
           Call get_word(record,word)
           parpot(5)=word_2_real(word)

           If (idnode == 0 .and. l_top) &
  Write(nrite,"(1x,i10,5x,3a8,3x,a4,1x,10f15.6)") itptbp,atom1,atom0,atom2,keyword,(parpot(j),j=1,mxptbp)

           katom0=0
           katom1=0
           katom2=0

           Do jtpatm=1,ntpatm
              If (atom0 == unqatm(jtpatm)) katom0=jtpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
           End Do

           If (katom0 == 0 .or. katom1 == 0 .or. katom2 == 0) Call error(84)

           lfrtbp(katom0)=.true.
           lfrtbp(katom1)=.true.
           lfrtbp(katom2)=.true.

           ka1=Max(katom1,katom2)
           ka2=Min(katom1,katom2)

           keytbp=(ka1*(ka1-1))/2+ka2+(katom0-1)*mx2tbp

           If (keytbp > mxtbp) Call error(86)

! convert parameters to internal units and angles to radians

           parpot(1) = parpot(1)*engunit
           If (keypot /= 6) parpot(2) = parpot(2)*(pi/180.0_wp)

           If (lsttbp(keytbp) > 0) Call error(18)

           lsttbp(keytbp)=itptbp
           ktbp=mx2tbp*((keytbp-1)/mx2tbp)+1 ! == mx2tbp*(katom0-1)+1
           If (lsttbp(ktbp) < 0) lsttbp(ktbp)=0

           ltptbp(itptbp)=keypot

! calculate max three-body cutoff

           rctbp=Max(rctbp,parpot(5))
           If (parpot(5) < 1.0e-6_wp) Then
              rcttbp(itptbp)=rctbp
              parpot(5)=rctbp
           Else
              rcttbp(itptbp)=parpot(5)
           End If

! store three-body potential parameters

           Do i=1,mxptbp
              prmtbp(i,itptbp)=parpot(i)
           End Do
        End Do

        If (ntptbp > 0) Then
           If (rctbp < 1.0e-6_wp) Call error(451)
           If (rcut < 2.0_wp*rctbp) Call error(471)
        End If

! read in the four-body potential energy parameters

     Else If (word(1:3) == 'fbp') Then

        Call get_word(record,word)
        ntpfbp=Nint(word_2_real(word))

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'number of specified four-body potentials ',i10)") ntpfbp
           If (l_top) &
  Write(nrite,"(/,4x,'quartet',5x,'atom 1',2x,'atom 2',2x,'atom 3',2x,'atom 4',5x,'key',14x,'parameters',/)")
        End If

        If (ntpfbp > mxfbp) Call error(89)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

! Initialise head of group atom to not participating in an interaction

        Do ifbp=1,mxfbp,mx3fbp
           lstfbp(ifbp)=-1
        End Do

        Do itpfbp=1,ntpfbp

           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

! Note the order!! atom0 is the central atom

           atom0=word(1:8)
           Call get_word(record,word)
           atom1=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)
           Call get_word(record,word)
           atom3=word(1:8)

           Call get_word(record,word)
           Call lower_case(word)
           keyword=word(1:4)

           If      (keyword == 'harm') Then
              keypot=1
           Else If (keyword == 'hcos') Then
              keypot=2
           Else If (keyword == 'plan') Then
              keypot=3
           Else

              If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
              Call error(443)

           End If

           Call get_word(record,word)
           parpot(1)=word_2_real(word)
           Call get_word(record,word)
           parpot(2)=word_2_real(word)
           Call get_word(record,word)
           parpot(3)=word_2_real(word)

           If (idnode == 0 .and. l_top) &
  Write(nrite,"(1x,i10,3x,4a8,3x,a4,2x,10f15.6)") itpfbp,atom0,atom1,atom2,atom3,keyword,(parpot(j),j=1,mxpfbp)

           katom0=0
           katom1=0
           katom2=0
           katom3=0

           Do jtpatm=1,ntpatm
              If (atom0 == unqatm(jtpatm)) katom0=jtpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
              If (atom3 == unqatm(jtpatm)) katom3=jtpatm
           End Do

           If (katom0 == 0 .or. katom1 == 0 .or. katom2 == 0 .or. katom3 == 0) Call error(91)

           lfrfbp(katom0)=.true.
           lfrfbp(katom1)=.true.
           lfrfbp(katom2)=.true.
           lfrfbp(katom3)=.true.

           ka1=Max(katom1,katom2,katom3)
           ka3=Min(katom1,katom2,katom3)
           ka2=katom1+katom2+katom3-ka1-ka3

           keyfbp=ka3+(ka2*(ka2-1))/2+(ka1*(ka1**2-1))/6+(katom0-1)*mx3fbp

           If (keyfbp > mxfbp) Call error(101)

! convert parameters to internal units and angles to radians

           parpot(1) = parpot(1)*engunit
           parpot(2) = parpot(2)*(pi/180.0_wp)

           If (keypot == 2) Then
              parpot(2)=Cos(parpot(2))
           End If

           If (lstfbp(keyfbp) > 0) Call error(19)

           lstfbp(keyfbp)=itpfbp
           kfbp=mx3fbp*((keyfbp-1)/mx3fbp)+1 ! == mx3fbp*(katom0-1)+1
           If (lstfbp(kfbp) < 0) lstfbp(kfbp)=0

           ltpfbp(itpfbp)=keypot

! calculate max four-body cutoff

           rcfbp=Max(rcfbp,parpot(3))
           If (parpot(3) < 1.0e-6_wp) Then
              rctfbp(itpfbp)=rcfbp
              parpot(3)=rcfbp
           Else
              rctfbp(itpfbp)=parpot(3)
           End If

! store four-body potential parameters

           Do i=1,mxpfbp
              prmfbp(i,itpfbp)=parpot(i)
           End Do
        End Do

        If (ntpfbp > 0) Then
           If (rcfbp < 1.0e-6_wp) Call error(453)
           If (rcut < 2.0_wp*rcfbp) Call error(472)
        End If

! read external field data

     Else If (word(1:6) == 'extern') Then

        Call get_word(record,word)
        nfld=0
        If (word(1:1) /= '#' .and. word(1:1) /= ' ') nfld=Nint(word_2_real(word))
        If (nfld <= 0) nfld=5

        word(1:1)='#'
        Do While (word(1:1) == '#' .or. word(1:1) == ' ')
           Call get_line(safe,nfield,record)
           If (.not.safe) Go To 2000
           Call get_word(record,word)
        End Do

        Call lower_case(word)
        keyword=word(1:4)

        If      (keyword == 'elec') Then
           keyfld=1
        Else If (keyword == 'oshr') Then
           keyfld=2
        Else If (keyword == 'shrx') Then
           keyfld=3
        Else If (keyword == 'grav') Then
           keyfld=4
        Else If (keyword == 'magn') Then
           keyfld=5
        Else If (keyword == 'sphr') Then
           keyfld=6
        Else If (keyword == 'zbnd') Then
           keyfld=7
        Else

           If (idnode == 0) Write(nrite,'(/,1x,a)') keyword
           Call error(454)

        End If

        Do i=1,nfld
           Call get_word(record,word)
           prmfld(i)=word_2_real(word)
        End Do

        If (idnode == 0) Then
  Write(nrite,"(/,/,1x,'external field key ',13x,a4)") keyword
           If (l_top) Then
  Write(nrite,"(/,30x,'parameters')")
  Write(nrite,"(/,1x,10f15.6)") prmfld
           End If
        End If

! convert to internal units

        If (keyfld == 1 .or. keyfld == 4 .or. keyfld == 5) Then
           If (.not.lunits) Call error(6)

           Do i=1,3
              prmfld(i) = prmfld(i)*engunit
           End Do
        Else If (keyfld == 2 .or. keyfld == 6 .or. keyfld == 7) Then
           prmfld(1) = prmfld(1)*engunit
        End If

        If (idnode == 0) Then
           If (keyfld == 2 .and. (imcon /= 1 .and. imcon /= 2)) &
  Write(nrite,"(/,1x,a)") '*** warning - external field only applicable for imcon=1,2 (orthorhombic geometry)!!!'
             If (keyfld == 3 .and. imcon /= 6) &
  Write(nrite,"(/,/,1x,a)") '*** warning - external field only applicable for imcon=6 (SLAB geometry)!!!'
        End If

! close force field file

     Else If (word(1:5) == 'close') Then

        If (idnode == 0) Close(Unit=nfield)

! Precautions: (vdw,met) = rdf scanning

        If (ntprdf == 0 .and. (ntpvdw > 0 .or. ntpmet > 0)) Then
           Do ia=1,ntpatm
              Do ja=ia,ntpatm
                 keyrdf=(ja*(ja-1))/2+ia
                 i=0
                 If (ntpvdw > 0) i=Max(i,lstvdw(keyrdf))
                 If (ntpmet > 0) i=Max(i,lstmet(keyrdf))
                 If (i > 0) Then
                    ntprdf = ntprdf+1
                    lstrdf(keyrdf) = ntprdf
                 End If
              End Do
           End Do
        End If

! Precautions: if vdw are cancelled, nullify ntpvdw as
! it is a switch for vdw_lrc and vdw_forces

        If (l_n_v) ntpvdw = 0

! check and resolve any conflicting 14 dihedral specifications

        Call dihedrals_14_check(l_str,ntpmls,nummols,numang,numdih,lstang,lstdih,prmdih)

! test for existence/appliance of any two-body or tersoff interactions!!!

        If ( keyfce == 0 .and. ntpvdw == 0 .and. &
             ntpmet == 0 .and. ntpter == 0 ) Call error(145)

! check for use of fast Global_To_Local search for special systems
!
!        If (megatm <= 50000000 .and.               &
!            (Int(megatm*1,ip)+Int(megfrz*1,ip)   + &
!             Int(megshl*2,ip)+Int(megcon*2,ip)   + &
!             Int(megpmf*(mxtpmf(1)+mxtpmf(2)),ip)+ &
!             Int(megrgd*mxlrgd,ip)+                &
!             Int(megtet*1,ip)+                     &
!             Int(megbnd*2,ip)+Int(megang*3,ip)   + &
!             Int(megdih*4,ip)+Int(meginv*4,ip))  / &
!            Int(megatm,ip) > Int(2,ip)) gtl_b=megatm

! EXIT IF ALL IS OK

        Return

     Else

! error exit for unidentified directive

        If (idnode == 0) Write(nrite,'(/,1x,a)') word(1:Len_Trim(word)+1)
        Call error(4)

     End If
  End Do

! uncontrolled error exit from field file processing

  If (idnode == 0) Close(Unit=nfield)
  Call error(16)

! end of field file error exit

2000 Continue

  If (idnode == 0) Close(Unit=nfield)
  Call error(52)

End Subroutine read_field
