Module ffield
  Use kinds, Only : wp
  Use comms, Only : comms_type
  Use setup
  Use kim,   Only : kimim,rkim,kim_cutoff

! SITE MODULE

  Use site

! CONFIG MODULE
! Fuchs correction of charge non-neutral systems
! Global_To_Local variables

  Use configuration, Only : imcon,cell,sumchg

! DPD module

  Use dpd,    Only : sigdpd

! INTERACTION MODULES

  Use core_shell
  Use mpole, Only : keyind,thole,mpllfr,plrsit,dmpsit

  Use pmf

  Use rigid_bodies

  Use tethers, Only : tethers_type

  Use bonds, Only : bonds_type,bonds_table_read,allocate_bond_dst_arrays
  Use angles, Only : angles_type,angles_table_read,allocate_angl_dst_arrays
  Use dihedrals, Only : dihedrals_type,dihedrals_table_read,allocate_dihd_dst_arrays, &
                        dihedrals_14_check
  Use inversions, Only : inversions_type,inversions_table_read,allocate_invr_dst_arrays
  Use three_body, Only : threebody_type
  Use vdw
  Use metal, Only : metal_type,allocate_metal_arrays,allocate_metal_table_arrays, &
                    metal_generate_erf,metal_table_read,metal_generate
  Use tersoff
  Use four_body

  Use external_field

  Use kinetics, Only : l_vom

! RDF MODULE

  Use rdfs

! PARSE MODULE

  Use parse

  Use numerics, Only : factorial, shellsort
  Use errors_warnings, Only : error,warning,info
  Use thermostat, Only : thermostat_type,ENS_NVE
  Use constraints, Only : constraints_type
  Use pmf, Only : pmf_type

  Implicit None

  Private
  Public :: read_field
  Public :: read_mpoles
  Public :: report_topology
  Public :: scan_field
Contains
Subroutine read_field                      &
           (l_str,l_top,l_n_v,             &
           rcut,rvdw,width,epsq, &
           keyfce,keyshl,           &
           lecx,lbook,lexcl,               &
           rcter,rcfbp,              &
           atmfre,atmfrz,megatm,megfrz,    &
           megshl,megrgd,    &
           megtet,    &
           pmf,cons,  &
           thermo,met,bond,angle,dihedral,inversion,tether,threebody,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading in the molecular specifications
! of the system to be simulated
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2017
! contrib   - r.davidchak (eeam) july 2012
! contrib   - b.palmer (2band) may 2013
! contrib   - a.v.brukhno & i.t.todorov march 2014 (itramolecular TPs & PDFs)
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena february 2017
! contrib   - a.m.elena august 2017 (mstw)
! contrib   - a.m.elena september 2017 (ryd)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (blb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES


  Logical,           Intent( In    ) :: l_str,l_top,l_n_v
  Integer,           Intent( InOut ) :: keyfce
  Real( Kind = wp ), Intent( In    ) :: rcut,rvdw,width,epsq
  Logical,           Intent( InOut ) :: lecx

  Logical,           Intent(   Out ) :: lbook,lexcl
  Real( Kind = wp ), Intent(   Out ) :: rcter,rcfbp
  Integer,           Intent(   Out ) :: keyshl,                             &
                                        atmfre,atmfrz,megatm,megfrz,        &
                                        megshl,megrgd,        &
                                        megtet
  Type( constraints_type ), Intent( InOut ) :: cons
  Type( pmf_type ), Intent( InOut ) :: pmf
  Type( thermostat_type ), Intent( InOut ) :: thermo
  Type( metal_type ), Intent( InOut ) :: met
  Type( bonds_type ), Intent( InOut ) :: bond
  Type( angles_type ), Intent( InOut ) :: angle
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( inversions_type ), Intent( InOut ) :: inversion
  Type( tethers_type), Intent( InOut ) :: tether
  Type( threebody_type), Intent( InOut ) :: threebody
  Type( comms_type), Intent( InOut ) :: comm

  Logical                :: safe,lunits,lmols,atmchk,                        &
                            l_shl,l_con,l_rgd,l_tet,l_bnd,l_ang,l_dih,l_inv, &
                            lshl_one,lshl_all,lshl_abort,lpmf,               &
                            lmet_safe,lter_safe,ldpd_safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 32  ) :: iddihd,idinvr
  Character( Len = 24  ) :: idangl
  Character( Len = 16  ) :: idbond
  Character( Len = 8   ) :: atom0,atom1,atom2,atom3
  Character( Len = 4   ) :: keyword
  Integer                :: fail(1:4),itmols,isite,jsite,ksite,msite,nsite,     &
                            isite1,isite2,isite3,isite4,is(0:4),js(0:4),        &
                            irept,nrept,ifrz,ntmp,i,j,ia,ja,jtpatm,ntab,keypot, &
                            iatm1,iatm2,iatm3,iatm4,                            &
                            katom0,katom1,katom2,katom3,katom4,                 &
                            ishls,nshels,icnst,nconst,frzcon,ipmf,jpmf,kpmf,    &
                            nrigid,irgd,jrgd,krgd,lrgd,frzrgd,iteth,nteth,      &
                            ibond,nbonds,ntpbnd,iang,nangle,ntpang,             &
                            idih,ndihed,ntpdih,iinv,ninver,ntpinv,              &
                            itprdf,keyrdf,itpvdw,keyvdw,itpmet,keymet,          &
                            itpter,keyter,icross,                               &
                            itbp,itptbp,keytbp,ktbp,                            &
                            ifbp,itpfbp,keyfbp,ka1,ka2,ka3,kfbp,nfld,itmp
  Real( Kind = wp )      :: weight,charge,pmf_tmp(1:2),parpot(1:30),tmp,        &
                            sig(0:2),eps(0:2),del(0:2),                         &
                            q_core_p, q_core_s , q_core,                        &
                            q_shel_p, q_shel_s , q_shel,                        &
                            k_crsh_p, k_crsh_s , k_crsh,                        &
                            p_core_p, p_core_s , p_core,                        &
                            d_core_p, d_core_s , d_core

  Character ( Len = 32 ), Allocatable :: dihd_name(:),invr_name(:)
  Character ( Len = 24 ), Allocatable :: angl_name(:)
  Character ( Len = 16 ), Allocatable :: bond_name(:)
  Character ( Len = 256 )             :: message,messages(3)
  Integer                             :: iter,rwidth
  Character ( Len = 100 )             :: rfmt
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

! Allocate auxiliary identity arrays for intramolecular TPs and PDFs

  fail = 0
  If (bond%l_tab .or. bond%bin_pdf > 0) Allocate (bond_name(1:bond%max_types), Stat=fail(1))
  If (angle%l_tab .or. angle%bin_adf > 0) Allocate (angl_name(1:angle%max_types), Stat=fail(2))
  If (dihedral%l_tab .or. dihedral%bin_adf > 0) Allocate (dihd_name(1:dihedral%max_types), Stat=fail(3))
  If (inversion%l_tab .or. inversion%bin_adf > 0) Allocate (invr_name(1:inversion%max_types), Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'read_field allocation failure'
     Call error(0,message)
  End If

! Initialise number of selected unique intramolecular TPs
! and corresponding identity arrays

  If (bond%l_tab) Then
     ntpbnd = 0 ! TABBND
     bond_name = ' '
  End If
  If (angle%l_tab) Then
     ntpang = 0 ! TABANG
     angl_name = ' '
  End If
  If (dihedral%l_tab) Then
     ntpdih = 0 ! TABDIH
     dihd_name = ' '
  End If
  If (inversion%l_tab) Then
     ntpinv = 0 ! TABINV
     invr_name = ' '
  End If

! Initialise total number of particles(shells are particles in MD),
! frozen particles, shells, constraints, PMFs, RBs,
! tethers, bonds, angles, dihedrals and inversions in the system.
! Some are used as switches for various force calculations.

  megatm = 0
  megfrz = 0

  megshl = 0

  cons%megcon = 0
  pmf%megpmf = 0

  megrgd = 0

  megtet = 0

  bond%total = 0
  angle%total = 0
  dihedral%total = 0
  inversion%total = 0

! Default for existence of intra-like interactions (including
! shells and tethers)

  lbook = .false.

! Default flag for existence of excluded intra-interaction

  lexcl = .false.

! default shell model - none

  keyshl   = 0
  lshl_one = .false. ! A massless shell existence indicator
  lshl_all = .true.  ! All shells are massless indicator

! Default potential cutoffs

  rcter=0.0_wp
  threebody%cutoff=0.0_wp
  rcfbp=0.0_wp


! open force field data file

  If (comm%idnode == 0) Then
     Open(Unit=nfield, File = Trim(field), Status = 'old')
  End If
  Call info(' ',.true.)
  Call info('system specification',.true.)
  If (.not.l_top) Then
    Call info('detailed topology opted out',.true.)
  End If

! omit first line

  Call get_line(safe,nfield,record,comm)
  If (.not.safe) Go To 2000

! read and process directives from field file

  Do

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nfield,record,comm)
        If (.not.safe) Go To 2000
        Call get_word(record,word) ; Call lower_case(word)
     End Do

! energy unit for input/output

     If (word(1:5) == 'units') Then
        lunits=.true.
        Call get_word(record,word) ; Call lower_case(word)
        If (word(1:2) == 'ev') Then
           engunit = eu_ev
           Call info('energy units = eV',.true.)
        Else If (word(1:4) == 'kcal') Then
           engunit = eu_kcpm
           Call info('energy units = kcal/mol',.true.)
        Else If (word(1:2) == 'kj') Then
           engunit = eu_kjpm
           Call info('energy units = kJ/mol',.true.)
        Else If (word(1:8) == 'internal') Then
           Call info('energy units = dl_poly internal units (10 J/mol)',.true.)
        Else If (word(1:1) == 'k') Then
           engunit = boltz
           Call info('energy units = Kelvin/Boltzmann',.true.)
        Else If (word(1:1) == ' ') Then
           Call info('energy units = dl_poly internal units (10 J/mol)',.true.)
        Else
           Call info(word(1:Len_Trim(word)),.true.)
           Call error(5)
        End If

! multipolar electrostatics control option

     Else If (word(1:5) == 'multi') Then
        Call get_word(record,word)
        Call lower_case(word)
        If (word(1:5) == 'order') Then
          Call get_word(record,word)
        End If
        Write(message,'(a,i0)') &
           "Multipolar electrostatics opted with poles up to order ", Nint(word_2_real(word))
        Call info(message,.true.)

        If (Nint(word_2_real(word)) > mxompl) Then
           Call warning("supplied multipolar expansion order reduced to the maximum allowed - 4",.true.)
        End If

        Call info("MPOLES file scheduled for reading after reading all intramolecular topology in FIELD",.true.)

        If (mximpl == 0) Then
          Call info('MPOLES file scheduling abandoned due to the "no elec" option in CONTROL',.true.)
        End If

! neutral group control option

     Else If (word(1:4) == 'neut') Then
        Call error(26)

! specify molecular species

     Else If (word(1:7) == 'molecul') Then
        Call get_word(record,word)
        Call lower_case(word)
        If (word(1:4) == 'type') Then
          Call get_word(record,word)
        End If

! number of molecular types

        If (lmols) Then
          Call error(11)
        End If
        lmols=.true.
        ntpmls=Nint(word_2_real(word))

        Write(message,'(a,6x,i10)') 'number of molecular types', ntpmls
        Call info(message,.true.)

        If (ntpmls > mxtmls) Then
          Call error(10)
        End If

! read in molecular characteristics for every molecule

        Do itmols=1,ntpmls

! initialise frozen constraints & RBs counters

           frzcon=0
           frzrgd=0

! expectation values for once off definitions of bonded quantities
! PMFs make an exception (as defined once and only once per system)

           l_shl=.true. ; l_con=.true. ; l_rgd=.true. ; l_tet=.true.
           l_bnd=.true. ; l_ang=.true. ; l_dih=.true. ; l_inv=.true.

           Write(message,'(a,9x,i10)') 'molecular species type', itmols
           Call info(message,.true.)

! name of molecular species

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do
           Call strip_blanks(record)
           molnam(itmols)=word(1:Len_Trim(word)+1)//record

           Write(message,'(a,13x,a40)') 'name of species:', molnam(itmols)

! stop processing if energy unit has not been specified

           If (.not.lunits) Call error(6)

! read molecular data

           Do

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record,comm)
                 If (.not.safe) Go To 2000
                 Call lower_case(record)
                 Call get_word(record,word)
              End Do

! number of molecules of this type

              If (word(1:6) == 'nummol') Then

                 Call get_word(record,word)
                 nummols(itmols)=Nint(word_2_real(word))

                 Write(message,'(a,10x,i10)') 'number of molecules  ', nummols(itmols)

! read in atomic details

              Else If (word(1:5) == 'atoms') Then

                 Call get_word(record,word)
                 numsit(itmols)=Nint(word_2_real(word))

                 Write(message,'(a,10x,i10)') 'number of atoms/sites', numsit(itmols)
                 Call info(message,.true.)

                 Write(messages(1),'(a)') 'atomic characteristics:'
                 Write(messages(2),'(8x,a4,4x,a4,13x,a4,9x,a6,6x,a6,6x,a6)') &
                  'site','name','mass','charge','repeat','freeze'
                 Call info(messages,2,.true.)

! for every molecule of this type get site and atom description

! reference point

                 ksite=0
                 Do isite=1,numsit(itmols)
                    If (ksite < numsit(itmols)) Then

! read atom name, mass, charge, repeat, freeze option

                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record,comm)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       atom1=word(1:8)

                       Call get_word(record,word)
                       weight=Abs(word_2_real(word))

                       Call get_word(record,word)
                       charge=word_2_real(word)

                       Call get_word(record,word)
                       nrept=Nint(word_2_real(word))
                       If (nrept == 0) nrept=1

! sum absolute charges

                       sumchg=sumchg+Abs(charge)

                       Call get_word(record,word)
                       ifrz=Nint(word_2_real(word))
                       If (ifrz /= 0) ifrz=1

                       numfrz(itmols)=numfrz(itmols)+ifrz*nrept

                       Write(message,'(2x,i10,4x,a8,2f15.6,2i10)') &
                         ksite+1,atom1,weight,charge,nrept,ifrz
                       Call info(message,.true.)

                       Do irept=1,nrept
                          ksite=ksite+1
                          If (ksite > numsit(itmols)) Call error(21)

                          nsite=nsite+1
                          If (nsite > mxsite) Call error(20)

                          sitnam(nsite)=atom1
                          wgtsit(nsite)=weight
                          chgsit(nsite)=charge
                          frzsit(nsite)=ifrz
                          If (wgtsit(nsite) > 1.0e-6_wp) dofsit(nsite)=3.0_wp*Real(Abs(1-ifrz),wp)
                       End Do

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

                 Write(message,'(a,5x,i10)') 'number of core-shell units', ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Call info('core-shell details:',.true.)
                   Write(message,'(8x,a4,5x,a5,5x,a5,5x,a10)') &
                     'unit','index','index','parameters'
                   Call info(message,.true.)
                 End If

                 Do ishls=1,numshl(itmols)
                    nshels=nshels+1
                    If (nshels > mxtshl) Call error(57)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
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
                    prmshl(1,nshels)=word_2_real(word)

                    Call get_word(record,word)
                    prmshl(2,nshels)=word_2_real(word)

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2

! test for frozen core-shell unit and print unit

                    If (l_top) Then
                      If (frzsit(isite1)*frzsit(isite2) /= 0) Then
                        Write(message,'(2x,3i10,2f15.6,1x,a8)') &
                          ishls,lstshl(1,nshels),lstshl(2,nshels), &
                          prmshl(1,nshels),prmshl(2,nshels),'*frozen*'
                      Else
                        Write(message,'(2x,3i10,2f15.6)') &
                          ishls,lstshl(1,nshels),lstshl(2,nshels), &
                          prmshl(1,nshels),prmshl(2,nshels)
                      End If
                      Call info(message,.true.)
                    End If

! catch unidentified entry

                    If (Any(lstshl(1:2,nshels) < 1) .or. Any(lstshl(1:2,nshels) > numsit(itmols))) Call error(27)

! abort if a shell is frozen

                    If (frzsit(isite2) /= 0) Call error(49)

! establish list of unique shell types (most certainly ntpshl <= ntpatm <= mxatyp)

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

                    prmshl(1:2,nshels)=prmshl(1:2,nshels)*engunit
                    smax=Max(smax,prmshl(1,nshels)+6.0_wp*prmshl(2,nshels))
                 End Do

! Check for mixed or multiple core-shell entries (no inter units linkage!)

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
                 cons%m_con=1

                 lbook=.true.
                 lexcl=.true.

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 ntmp=Nint(word_2_real(word))
                 cons%numcon(itmols)=cons%numcon(itmols)+ntmp

                 Write(message,'(a,5x,i10)') 'number of bond constraints', ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Call info('constraint bond details:',.true.)
                   Write(message,'(7x,a5,5x,a5,5x,a5,5x,a10)') &
                     'unit','index','index','bondlength'
                   Call info(message,.true.)
                 End If


                 Do icnst=1,cons%numcon(itmols)
                    nconst=nconst+1
                    If (nconst > cons%mxtcon) Call error(40)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read constraint bond atom indices

                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))

                    cons%lstcon(1,nconst)=iatm1
                    cons%lstcon(2,nconst)=iatm2

                    Call get_word(record,word)
                    cons%prmcon(nconst)=word_2_real(word)

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

! test for frozen atoms and print unit

                    If (l_top) Then
                      If (frzsit(isite1)*frzsit(isite2) /= 0) Then
                        Write(message,'(2x,3i10,f15.6,1x,a8)') &
                          icnst,cons%lstcon(1,nconst),cons%lstcon(2,nconst), &
                          cons%prmcon(nconst),'*frozen*'
                      Else
                        Write(message,'(2x,3i10,f15.6)') &
                          icnst,cons%lstcon(1,nconst),cons%lstcon(2,nconst), &
                          cons%prmcon(nconst)
                      End If
                      Call info(message,.true.)
                    End If

! catch unidentified entry

                    If (Any(cons%lstcon(1:2,nconst) < 1) .or. Any(cons%lstcon(1:2,nconst) > numsit(itmols))) Call error(27)

! test for mistyped constraint bond unit

                    If (iatm1 == iatm2) Call error(33)

! test for length of constraint bond unit > rcut

                    If (cons%prmcon(nconst) >= rcut) Call error(34)

                 End Do

! Check for multiple constraint bond entries

                 Do i=nconst-cons%numcon(itmols)+1,nconst
                    is(1)=Min(cons%lstcon(1,i),cons%lstcon(2,i))
                    is(2)=Max(cons%lstcon(1,i),cons%lstcon(2,i))

                    Do j=i+1,nconst
                       js(1)=Min(cons%lstcon(1,j),cons%lstcon(2,j))
                       js(2)=Max(cons%lstcon(1,j),cons%lstcon(2,j))

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
                 pmf%numpmf(itmols)=1          ! MD system is allowed

! read PMF bondlength

                 Call get_word(record,word)
                 pmf%prmpmf=word_2_real(word)

                 Call info('PMF constraint details',.true.)
                 Write(message,'(a,5x,f15.6)') 'bondlength:', pmf%prmpmf
                 Call info(message,.true.)

                 If (pmf%prmpmf > width/2.0_wp) Call error(480)

                 Do ipmf=1,2
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call lower_case(record)
                       Call get_word(record,word)
                    End Do

! read PMF indices and weights

                    pmf%pmffrz(ipmf) = 0
                    pmf_tmp(ipmf) = 0.0_wp

! test for zero length units

                    If (pmf%mxtpmf(ipmf) == 0) Then
                       Call strip_blanks(record)
                       Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                       Call info(message,.true.)
                       Call error(500)
                    End If

                    Do jpmf=1,pmf%mxtpmf(ipmf)
                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record,comm)
                          If (.not.safe) Go To 2000
                          Call get_word(record,word)
                       End Do

                       iatm1=Nint(word_2_real(word))
                       pmf%lstpmf(jpmf,ipmf)=iatm1
                       isite1 = nsite - numsit(itmols) + iatm1

! test for frozen units

                       pmf%pmffrz(ipmf)=pmf%pmffrz(ipmf)+frzsit(isite1)

                       Call get_word(record,word)
                       weight=word_2_real(word)
                       pmf%pmfwgt(jpmf,ipmf)=weight

! test for weightless units

                       pmf_tmp(ipmf)=pmf_tmp(ipmf)+weight
                    End Do
                 End Do

! get real masses
! if a PMF unit is massless supply weights from atoms' masses

                 Do ipmf=1,2
                    Do jpmf=1,pmf%mxtpmf(ipmf)
                       isite1 = nsite - numsit(itmols) + pmf%lstpmf(jpmf,ipmf)

                       pmf%pmfwg1(jpmf,ipmf)=wgtsit(isite1)
                       If (pmf_tmp(ipmf) < 1.0e-6_wp) pmf%pmfwgt(jpmf,ipmf)=wgtsit(isite1)
                    End Do

! if a PMF unit is still weightless set all members' masses to 1

                    If (pmf_tmp(ipmf) < 1.0e-6_wp) Then
                       pmf_tmp(ipmf)=Sum(pmf%pmfwgt(1:pmf%mxtpmf(ipmf),ipmf))

                       If (pmf_tmp(ipmf) < 1.0e-6_wp) Then
                          pmf%pmfwgt(:,ipmf)=1.0_wp
                          pmf%pmfwg1(:,ipmf)=1.0_wp
                          Call warning(230,Real(ipmf,wp),0.0_wp,0.0_wp)
                       End If
                    End If
                 End Do

! test for frozen atoms and print units

                 If (l_top) Then
                   Write(message,'(8x,a4,5x,a5,8x,a6)') &
                     'unit','index','weight'
                   Call info(message,.true.)

                    Do ipmf=1,2
                       Do jpmf=1,pmf%mxtpmf(ipmf)
                          isite1 = nsite - numsit(itmols) + pmf%lstpmf(jpmf,ipmf)
                          If (frzsit(isite1) /= 0) Then
                            Write(message,'(2x,2i10,f15.6,1x,a8)') &
                              ipmf,pmf%lstpmf(jpmf,ipmf),pmf%pmfwgt(jpmf,ipmf),'*frozen*'
                          Else
                            Write(message,'(2x,2i10,f15.6)') &
                              ipmf,pmf%lstpmf(jpmf,ipmf),pmf%pmfwgt(jpmf,ipmf)
                          End If
                          Call info(message,.true.)
                       End Do
                       If (ipmf == 1) Call info(' ',.true.)
                    End Do
                 End If

! catch unidentified entry

                 If ( (Any(pmf%lstpmf(1:pmf%mxtpmf(1),1) < 1) .or. Any(pmf%lstpmf(1:pmf%mxtpmf(1),1) > numsit(itmols))) .or. &
                      (Any(pmf%lstpmf(1:pmf%mxtpmf(2),2) < 1) .or. Any(pmf%lstpmf(1:pmf%mxtpmf(2),2) > numsit(itmols))) ) &
                         Call error(27)

! PMF reciprocal total unit masses

                 Do ipmf=1,2
                    pmf%pmfwgt(0,ipmf)=1.0_wp/Sum(pmf%pmfwgt(1:pmf%mxtpmf(ipmf),ipmf))
                    pmf%pmfwg1(0,ipmf)=1.0_wp/Sum(pmf%pmfwg1(1:pmf%mxtpmf(ipmf),ipmf))
                 End Do

! abort if there are frozen atoms on both PMF units

                 If (pmf%pmffrz(1)*pmf%pmffrz(2) > 0) Then
                    Call error(486)
                 Else
                    If (pmf%pmffrz(1) == pmf%mxtpmf(1) .or. pmf%pmffrz(2) == pmf%mxtpmf(2)) Then
                       charge=1.0_wp
                    Else
                       charge=0.5_wp
                    End If

                    Do ipmf=1,2
                       ntmp=pmf%mxtpmf(ipmf)-pmf%pmffrz(ipmf)
                       If (ntmp > 0) Then
                          tmp=charge/Real(ntmp,wp)
                          Do jpmf=1,pmf%mxtpmf(ipmf)
                             isite1 = nsite - numsit(itmols) + pmf%lstpmf(jpmf,ipmf)
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
                    Do jpmf=1,pmf%mxtpmf(ipmf)
                       Do kpmf=jpmf+1,pmf%mxtpmf(ipmf)
                          If (pmf%lstpmf(kpmf,ipmf) == pmf%lstpmf(jpmf,ipmf)) Call error(501)
                       End Do
                    End Do
                 End Do
                 Do jpmf=1,pmf%mxtpmf(2)
                    If (Any(pmf%lstpmf(1:pmf%mxtpmf(1),1) == pmf%lstpmf(jpmf,2))) Call error(502)
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

                 Write(message,'(a,9x,i10)') 'number of rigid bodies', ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Call info('rigid body details:',.true.)
                   Write(message,'(8x,a4,6x,a4,12x,a7)') &
                     'unit','size','indices'
                   Call info(message,.true.)
                 End If


                 Do irgd=1,numrgd(itmols)
                    nrigid=nrigid+1
                    If (nrigid > mxtrgd) Call error(630)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read RB size and indices

                    lrgd=Nint(word_2_real(word))
                    If (lrgd < 2 .or. lrgd > mxlrgd) Call error(632)

                    Do jrgd=1,lrgd
                       If (Mod(jrgd,16) == 0) Then
                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nfield,record,comm)
                             If (.not.safe) Go To 2000
                             Call get_word(record,word)
                          End Do
                       Else
                          Call get_word(record,word)
                       End If

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
                                           Real(1-rgdfrz(jrgd,nrigid),wp)*rgdwgt(jrgd,nrigid)
                       End Do
                       rgdwg1(0,nrigid)=Sum(rgdwg1(1:lrgd,nrigid))
                    Else
                       rgdwgt(0,nrigid)=Sum(rgdwgt(1:lrgd,nrigid))
                       rgdwg1(0:lrgd,nrigid)=rgdwgt(0:lrgd,nrigid)
                    End If

! print RB unit

                    If (comm%idnode == 0 .and. l_top) Then
                      rwidth = lrgd+1
                      Write(rfmt,'(a,i0,a)') '(26x,',rwidth,'i10)'

                      If (rgdfrz(0,nrigid) /= 0) Then
                        Write(message,'(2x,i10,2x,a)') irgd, '*frozen*'
                      Else
                        Write(message,'(2x,i10)') irgd
                      End If
                      Call info(message,.true.)

                      Do iter = 1, nrigid
                        Write(message,rfmt) lstrgd(0:lrgd,iter)
                        Call info(message,.true.)
                      End Do

! test for weightless RB

                       If (rgdwg1(0,nrigid) < 1.0e-6_wp) Call error(634)

! catch unidentified entry

                       If (Any(lstrgd(1:lrgd,nrigid) < 1) .or. Any(lstrgd(1:lrgd,nrigid) > numsit(itmols))) Call error(27)

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
                       Do jsite=1,krgd
                          If (Any(lstrgd(1:lrgd,irgd) == lstrgd(jsite,jrgd))) Then
                              Call warning(400,Real(irgd,wp),Real(jrgd,wp),0.0_wp)
                              Call error(620)
                          End If
                       End Do
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
                 tether%numteth(itmols)=tether%numteth(itmols)+ntmp

                 Write(message,'(a,7x,i10)') 'number of tethered sites', ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Call info('tethered site details:',.true.)
                   Write(message,'(8x,a4,5x,a3,6x,a4,5x,a10)') &
                     'unit','key','site','parameters'
                   Call info(message,.true.)
                 End If

                 Do iteth=1,tether%numteth(itmols)
                    nteth=nteth+1
                    If (nteth > tether%mxteth) Call error(62)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of tethering

                    Call lower_case(word)
                    keyword=word(1:4)

                    If      (keyword == 'harm') Then
                       tether%keytet(nteth)=1
                    Else If (keyword == 'rhrm') Then
                       tether%keytet(nteth)=2
                    Else If (keyword == 'quar') Then
                       tether%keytet(nteth)=3
                    Else

                       Call info(keyword,.true.)
                       Call error(450)

                    End If

! read tethered atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))

                    tether%lsttet(nteth)=iatm1

                    Call get_word(record,word)
                    tether%prmtet(1,nteth)=word_2_real(word)
                    Call get_word(record,word)
                    tether%prmtet(2,nteth)=word_2_real(word)
                    Call get_word(record,word)
                    tether%prmtet(3,nteth)=word_2_real(word)

                    isite1 = nsite - numsit(itmols) + iatm1

! test for frozen atom and print unit

                    If (l_top) Then
                       If (frzsit(isite1) /= 0) Then
                         Write(rfmt,'(a,i0,a)') '(2x,i10,a8,i10,2x,',tether%mxpteth,'f15.6,2x,a8)'
                         Write(message,rfmt) iteth,keyword,tether%lsttet(nteth),tether%prmtet(1:tether%mxpteth,nteth),'*frozen*'
                       Else
                         Write(rfmt,'(a,i0,a)') '(2x,i10,a8,i10,2x,',tether%mxpteth,'f15.6)'
                         Write(message,rfmt) iteth,keyword,tether%lsttet(nteth),tether%prmtet(1:tether%mxpteth,nteth)
                       End If
                       Call info(message,.true.)
                    End If

! catch unidentified entry

                    If (tether%lsttet(nteth) < 1 .or. tether%lsttet(nteth) > numsit(itmols)) Call error(27)

! convert energy units to internal units

                    tether%prmtet(:,nteth)=tether%prmtet(:,nteth)*engunit

                 End Do

! Check for multiple tether entries

                 Do i=nteth-tether%numteth(itmols)+1,nteth
                    is(0)=tether%keytet(i)
                    is(1)=tether%lsttet(i)

                    Do j=i+1,nteth
                       js(0)=tether%keytet(j)
                       js(1)=tether%lsttet(j)

                       If (js(1) == is(1)) Then
                          If (l_str .and. l_top) Call warning(410,Real(i,wp),Real(j,wp),0.0_wp)
                          If (is(0) == js(0)) Call error(620)
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
                 bond%num(itmols)=bond%num(itmols)+ntmp

                 Write(message,'(a,i10)') 'number of chemical bonds ',ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Call info('chemical bond details:',.true.)
                   Write(message,'(8x,a4,5x,a3,5x,a5,5x,a5,5x,a10)') &
                     'unit','key','index','index','parameters'
                   Call info(message,.true.)
                 End If

                 Do ibond=1,bond%num(itmols)
                    nbonds=nbonds+1
                    If (nbonds > bond%max_types) Call error(30)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of chemical bond

                    Call lower_case(word)
                    keyword=word(1:4)

                    If (keyword(1:1) /= '-') lexcl=.true.

                    If      (keyword == 'tab' ) Then
                       bond%key(nbonds)=20
                    Else If (keyword == '-tab') Then
                       bond%key(nbonds)=-20
                    Else If (keyword == 'harm') Then
                       bond%key(nbonds)=1
                    Else If (keyword == '-hrm') Then
                       bond%key(nbonds)=-1
                    Else If (keyword == 'mors') Then
                       bond%key(nbonds)=2
                    Else If (keyword == '-mrs') Then
                       bond%key(nbonds)=-2
                    Else If (keyword == '12-6') Then
                       bond%key(nbonds)=3
                    Else If (keyword == '-126') Then
                       bond%key(nbonds)=-3
                    Else If (keyword == 'lj'  ) Then
                       bond%key(nbonds)=4
                    Else If (keyword == '-lj' ) Then
                       bond%key(nbonds)=-4
                    Else If (keyword == 'rhrm') Then
                       bond%key(nbonds)=5
                    Else If (keyword == '-rhm') Then
                       bond%key(nbonds)=-5
                    Else If (keyword == 'quar') Then
                       bond%key(nbonds)=6
                    Else If (keyword == '-qur') Then
                       bond%key(nbonds)=-6
                    Else If (keyword == 'buck') Then
                       bond%key(nbonds)=7
                    Else If (keyword == '-bck') Then
                       bond%key(nbonds)=-7
                    Else If (keyword == 'coul') Then
                       bond%key(nbonds)=8
                    Else If (keyword == '-cul') Then
                       bond%key(nbonds)=-8
                    Else If (keyword == 'fene') Then
                       bond%key(nbonds)=9
                    Else If (keyword == '-fne') Then
                       bond%key(nbonds)=-9
                    Else If (keyword == 'mmst') Then
                       bond%key(nbonds)=10
                    Else If (keyword == '-mst') Then
                       bond%key(nbonds)=-10
                    Else
                       Call info(keyword,.true.)
                       Call error(444)
                    End If

! read bond atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))

                    bond%lst(1,nbonds)=iatm1
                    bond%lst(2,nbonds)=iatm2

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2

                    If (Abs(bond%key(nbonds)) /= 20) Then

                       Call get_word(record,word)
                       bond%param(1,nbonds)=word_2_real(word)
                       Call get_word(record,word)
                       bond%param(2,nbonds)=word_2_real(word)
                       Call get_word(record,word)
                       bond%param(3,nbonds)=word_2_real(word)
                       Call get_word(record,word)
                       bond%param(4,nbonds)=word_2_real(word)

                       If (Abs(bond%key(nbonds)) == 9) Then
                          bond%param(2,nbonds)=Abs(bond%param(2,nbonds))
                          If (Abs(bond%param(3,nbonds)) > bond%param(2,nbonds)/2.0_wp) &
                             bond%param(3,nbonds)=Sign(1.0_wp,bond%param(3,nbonds))*bond%param(2,nbonds)/2.0_wp
                       End If

! test for frozen atoms and print unit

                       If (l_top) Then
                          If (frzsit(isite1)*frzsit(isite2) /= 0) Then
                            Write(rfmt,'(a,i0,a)') '(2x,i10,a8,2i10,',bond%max_param,'f15.6,2x,a8)'
                            Write(message,rfmt) ibond,keyword,bond%lst(1,nbonds), &
                              bond%lst(2,nbonds),bond%param(1:bond%max_param,nbonds),'*frozen*'
                            Call info(message,.true.)
                          Else
                            Write(rfmt,'(a,i0,a)') '(2x,i10,a8,2i10,',bond%max_param,'f15.6)'
                            Write(message,rfmt) ibond,keyword,bond%lst(1,nbonds), &
                              bond%lst(2,nbonds),bond%param(1:bond%max_param,nbonds)
                            Call info(message,.true.)
                          End If
                       End If

! convert energy units to internal units

                       If (Abs(bond%key(nbonds)) /= 8) Then
                          bond%param(1,nbonds)=bond%param(1,nbonds)*engunit
                       End If

                       If (Abs(bond%key(nbonds)) == 3) Then
                          bond%param(2,nbonds)=bond%param(2,nbonds)*engunit
                       End If

                       If (Abs(bond%key(nbonds)) == 6) Then
                          bond%param(3,nbonds)=bond%param(3,nbonds)*engunit
                          bond%param(4,nbonds)=bond%param(4,nbonds)*engunit
                       End If

                       If (Abs(bond%key(nbonds)) == 7) Then
                          bond%param(3,nbonds)=bond%param(3,nbonds)*engunit
                       End If

                    Else ! TABBND to read

! Construct unique name for the tabulated bond

                       Do jsite=1,ntpatm
                          If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                          If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                       End Do

                       If (katom1 <= katom2) Then
                          idbond = sitnam(iatm1)//sitnam(iatm2)
                       Else
                          idbond = sitnam(iatm2)//sitnam(iatm1)
                       End If

! ntpbnd total number of unique table potentials to read from TABBND

                       Do i=1,ntpbnd
                          If (bond_name(i) == idbond) Then
                             bond%ltp(nbonds)=i ! Re-point from zero to type
                             Exit
                          End If
                       End Do

                       If (bond%ltp(nbonds) == 0) Then
                          ntpbnd=ntpbnd+1
                          bond_name(ntpbnd)=idbond

                          bond%ltp(0)=ntpbnd      ! NUTBP
                          bond%ltp(nbonds)=ntpbnd ! Re-point from zero to type
                       End If

                       If (l_top) Then
                          If (frzsit(isite1)*frzsit(isite2) /= 0) Then
                            Write(message,'(2x,i10,a8,2i10,2x,a9,2x,a8)') &
                              ibond,keyword,bond%lst(1,nbonds),bond%lst(2,nbonds), &
                              "tabulated",'*frozen*'
                            Call info(message,.true.)
                          Else
                            Write(message,'(2x,i10,a8,2i10,2x,a9)') &
                              ibond,keyword,bond%lst(1,nbonds),bond%lst(2,nbonds), &
                              "tabulated"
                            Call info(message,.true.)
                          End If
                       End If

                    End If

! catch unidentified entry

                    If (Any(bond%lst(1:2,nbonds) < 1) .or. Any(bond%lst(1:2,nbonds) > numsit(itmols))) Call error(27)

! test for mistyped chemical bond unit

                    If (iatm1 == iatm2) Call error(33)
                 End Do

! Check for multiple chemical bond entries

                 Do i=nbonds-bond%num(itmols)+1,nbonds
                    is(0)=bond%key(i)
                    is(1)=Min(bond%lst(1,i),bond%lst(2,i))
                    is(2)=Max(bond%lst(1,i),bond%lst(2,i))

                    Do j=i+1,nbonds
                       js(0)=bond%key(j)
                       js(1)=Min(bond%lst(1,j),bond%lst(2,j))
                       js(2)=Max(bond%lst(1,j),bond%lst(2,j))

                       If (js(1) == is(1) .and. js(2) == is(2)) Then
                          If (l_str .and. l_top) Call warning(420,Real(i,wp),Real(j,wp),0.0_wp)
                          If (is(0) == js(0)) Call error(620)
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
                 angle%num(itmols)=angle%num(itmols)+ntmp

                 Write(message,'(a,i10)') 'number of bond angles ',ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Write(messages(1),'(a)') 'bond angle details:'
                   Write(messages(2),'(8x,a4,5x,a3,3(5x,a5),8x,a7,8x,a5)') &
                     'unit','key','index','index','index','f-const','angle'
                   Call info(messages,2,.true.)
                 End If

                 Do iang=1,angle%num(itmols)
                    nangle=nangle+1
                    If (nangle > angle%max_types) Call error(50)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of bond angle

                    Call lower_case(word)
                    keyword=word(1:4)

                    If (keyword(1:1) /= '-') lexcl=.true.

                    If      (keyword == 'tab' ) Then
                       angle%key(nangle)=20
                    Else If (keyword == '-tab') Then
                       angle%key(nangle)=-20
                    Else If (keyword == 'harm') Then
                       angle%key(nangle)=1
                    Else If (keyword == '-hrm') Then
                       angle%key(nangle)=-1
                    Else If (keyword == 'quar') Then
                       angle%key(nangle)=2
                    Else If (keyword == '-qur') Then
                       angle%key(nangle)=-2
                    Else If (keyword == 'thrm') Then
                       angle%key(nangle)=3
                    Else If (keyword == '-thm') Then
                       angle%key(nangle)=-3
                    Else If (keyword == 'shrm') Then
                       angle%key(nangle)=4
                    Else If (keyword == '-shm') Then
                       angle%key(nangle)=-4
                    Else If (keyword == 'bvs1') Then
                       angle%key(nangle)=5
                    Else If (keyword == '-bv1') Then
                       angle%key(nangle)=-5
                    Else If (keyword == 'bvs2') Then
                       angle%key(nangle)=6
                    Else If (keyword == '-bv2') Then
                       angle%key(nangle)=-6
                    Else If (keyword == 'hcos') Then
                       angle%key(nangle)=7
                    Else If (keyword == '-hcs') Then
                       angle%key(nangle)=-7
                    Else If (keyword == 'cos' ) Then
                       angle%key(nangle)=8
                    Else If (keyword == '-cos') Then
                       angle%key(nangle)=-8
                    Else If (keyword == 'mmsb') Then
                       angle%key(nangle)=9
                    Else If (keyword == '-msb') Then
                       angle%key(nangle)=-9
                    Else If (keyword == 'stst') Then
                       angle%key(nangle)=10
                    Else If (keyword == '-sts') Then
                       angle%key(nangle)=-10
                    Else If (keyword == 'stbe') Then
                       angle%key(nangle)=11
                    Else If (keyword == '-stb') Then
                       angle%key(nangle)=-11
                    Else If (keyword == 'cmps') Then
                       angle%key(nangle)=12
                    Else If (keyword == '-cmp') Then
                       angle%key(nangle)=-12
                    Else If (keyword == 'mmbd') Then
                       angle%key(nangle)=13
                    Else If (keyword == '-mbd') Then
                       angle%key(nangle)=-13
                    Else If (keyword == 'kky' ) Then
                       angle%key(nangle)=14
                    Else If (keyword == '-kky') Then
                       angle%key(nangle)=-14
                    Else
                       Call info(keyword,.true.)
                       Call error(440)
                    End If

! read angle atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word)) ! central atom
                    Call get_word(record,word)
                    iatm3=Nint(word_2_real(word))

                    angle%lst(1,nangle)=iatm1
                    angle%lst(2,nangle)=iatm2
                    angle%lst(3,nangle)=iatm3

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2
                    isite3 = nsite - numsit(itmols) + iatm3

                    If (Abs(angle%key(nangle)) /= 20) Then

                       Call get_word(record,word)
                       angle%param(1,nangle)=word_2_real(word)
                       Call get_word(record,word)
                       angle%param(2,nangle)=word_2_real(word)
                       Call get_word(record,word)
                       angle%param(3,nangle)=word_2_real(word)
                       Call get_word(record,word)
                       angle%param(4,nangle)=word_2_real(word)
                       Call get_word(record,word)
                       angle%param(5,nangle)=word_2_real(word)
                       Call get_word(record,word)
                       angle%param(6,nangle)=word_2_real(word)

! test for frozen atoms and print unit

                       If (l_top) Then
                         If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3) /= 0) Then
                           write(rfmt,'(a,i0,a)') '(2x,i10,a8,3i10,',angle%max_param,'f15.6,2x,a8)'
                           write(message,rfmt) iang,keyword,angle%lst(1:3,nangle),angle%param(1:angle%max_param,nangle),'*frozen*'
                         Else
                           write(rfmt,'(a,i0,a)') '(2x,i10,a8,3i10,',angle%max_param,'f15.6)'
                           write(message,rfmt) iang,keyword,angle%lst(1:3,nangle),angle%param(1:angle%max_param,nangle)
                         End If
                         Call info(message,.true.)
                       End If

! convert energies to internal units

                       angle%param(1,nangle) = angle%param(1,nangle)*engunit

                       If      (Abs(angle%key(nangle)) == 2) Then
                          angle%param(3,nangle) = angle%param(3,nangle)*engunit
                          angle%param(4,nangle) = angle%param(4,nangle)*engunit
                       Else If (Abs(angle%key(nangle)) == 12) Then
                          angle%param(2,nangle) = angle%param(2,nangle)*engunit
                          angle%param(3,nangle) = angle%param(3,nangle)*engunit
                       End If

! convert angles to radians

                       If      (Abs(angle%key(nangle)) == 12) Then
                          angle%param(4,nangle) = angle%param(4,nangle)*(pi/180.0_wp)
                       Else If (Abs(angle%key(nangle)) /= 10) Then
                          angle%param(2,nangle) = angle%param(2,nangle)*(pi/180.0_wp)
                       End If

                    Else ! TABANG to read

! Construct unique name for the tabulated angle

                       Do jsite=1,ntpatm
                          If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                          If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                       End Do

                       If (katom1 <= katom3) Then
                          idangl = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)
                       Else
                          idangl = sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm1)
                       End If

! ntpang total number of unique table potentials to read from TABANG

                       Do i=1,ntpang
                          If (angl_name(i) == idangl) Then
                             angle%ltp(nangle)=i ! Re-point from zero to type
                             Exit
                          End If
                       End Do

                       If (angle%ltp(nangle) == 0) Then
                          ntpang=ntpang+1
                          angl_name(ntpang)=idangl

                          angle%ltp(0)=ntpang      ! NUTAP
                          angle%ltp(nangle)=ntpang ! Re-point from zero to type
                       End If

! test for frozen atoms and print unit

                       If (l_top) Then
                          If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3) /= 0) Then
                            Write(message,'(2x,i10,a8,3i10,2x,a9,2x,a8)') &
                              iang,keyword,angle%lst(1:3,nangle),'tabulated','*frozen*'
                          Else
                            Write(message,'(2x,i10,a8,3i10,2x,a9)') &
                              iang,keyword,angle%lst(1:3,nangle),'tabulated'
                          End If
                          Call info(message,.true.)
                       End If

                    End If

! catch unidentified entry

                    If (Any(angle%lst(1:3,nangle) < 1) .or. Any(angle%lst(1:3,nangle) > numsit(itmols))) Call error(27)

! test for mistyped bond angle unit

                    If (iatm1 == iatm2 .or. iatm1 == iatm3 .or. iatm2 == iatm3) Call error(66)
                 End Do

! Check for multiple bond angle entries

                 Do i=nangle-angle%num(itmols)+1,nangle
                    is(0)=angle%key(i)
                    is(1)=Min(angle%lst(1,i),angle%lst(3,i))
                    is(2)=angle%lst(2,i)
                    is(3)=Max(angle%lst(1,i),angle%lst(3,i))

                    Do j=i+1,nangle
                       js(0)=angle%key(j)
                       js(1)=Min(angle%lst(1,j),angle%lst(3,j))
                       js(2)=angle%lst(2,j)
                       js(3)=Max(angle%lst(1,j),angle%lst(3,j))

                       If (js(1) == is(1) .and. js(2) == is(2) .and. &
                           js(3) == is(3)) Then
                          If (l_str .and. l_top) Call warning(430,Real(i,wp),Real(j,wp),0.0_wp)
                          If (is(0) == js(0)) Call error(620)
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
                 dihedral%num(itmols)=dihedral%num(itmols)+ntmp

                 Write(message,'(a,i10)') 'number of dihedral angles ',ntmp
                 Call info(message,.true.)

                 If (l_top) Then
                   Write(messages(1),'(a)') 'dihedral angle details:'
                   Write(messages(2),'(8x,a4,5x,a3,4(5x,a5),8x,a7,10x,a5,11x,a4,7x,a8,8x,a7)') &
                     'unit','key','index','index','index','index','f-const', &
                     'angle','trig','1-4 elec','1-4 vdw'
                   Call info(messages,2,.true.)
                 End If

                 Do idih=1,dihedral%num(itmols)
                    ndihed=ndihed+1
                    If (ndihed > dihedral%max_types) Call error(60)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of dihedral bond angle

                    Call lower_case(word)
                    keyword=word(1:4)
                    If      (keyword == 'tab' ) Then
                       dihedral%key(ndihed)=20
                    Else If (keyword == '-tab') Then
                       dihedral%key(ndihed)=-20
                    Else If (keyword == 'cos' ) Then
                       dihedral%key(ndihed)=1
                    Else If (keyword == '-cos') Then
                       dihedral%key(ndihed)=-1
                    Else If (keyword == 'harm') Then
                       dihedral%key(ndihed)=2
                    Else If (keyword == '-hrm') Then
                       dihedral%key(ndihed)=-2
                    Else If (keyword == 'hcos') Then
                       dihedral%key(ndihed)=3
                    Else If (keyword == '-hcs') Then
                       dihedral%key(ndihed)=-3
                    Else If (keyword == 'cos3') Then
                       dihedral%key(ndihed)=4
                    Else If (keyword == '-cs3') Then
                       dihedral%key(ndihed)=-4
                    Else If (keyword == 'ryck') Then
                       dihedral%key(ndihed)=5
                    Else If (keyword == '-rck') Then
                       dihedral%key(ndihed)=-5
                    Else If (keyword == 'rbf' ) Then
                       dihedral%key(ndihed)=6
                    Else If (keyword == '-rbf') Then
                       dihedral%key(ndihed)=-6
                    Else If (keyword == 'opls') Then
                       dihedral%key(ndihed)=7
                    Else If (keyword == '-opl') Then
                       dihedral%key(ndihed)=-7
                    Else
                       Call info(keyword,.true.)
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

                    dihedral%lst(1,ndihed)=iatm1
                    dihedral%lst(2,ndihed)=iatm2
                    dihedral%lst(3,ndihed)=iatm3
                    dihedral%lst(4,ndihed)=iatm4

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2
                    isite3 = nsite - numsit(itmols) + iatm3
                    isite4 = nsite - numsit(itmols) + iatm4

                    If (dihedral%key(ndihed) /= 20) Then

                       Call get_word(record,word)
                       dihedral%param(1,ndihed)=word_2_real(word)
                       Call get_word(record,word)
                       dihedral%param(2,ndihed)=word_2_real(word)
                       Call get_word(record,word)
                       dihedral%param(3,ndihed)=word_2_real(word)
                       Call get_word(record,word)
                       dihedral%param(4,ndihed)=word_2_real(word)
                       Call get_word(record,word)
                       dihedral%param(5,ndihed)=word_2_real(word)
                       Call get_word(record,word)
                       dihedral%param(6,ndihed)=word_2_real(word)
                       Call get_word(record,word)
                       dihedral%param(7,ndihed)=word_2_real(word)

! test for frozen atoms and print unit

                       If (l_top) Then
                          If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) /= 0) Then
                            Write(rfmt,'(a,i0,a)') '(2x,i10,a8,4i10,',dihedral%max_param,'f15.6,2x,a8)'
                            Write(message,rfmt) idih,keyword,dihedral%lst(1:4,ndihed), &
                              dihedral%param(1:dihedral%max_param,ndihed),'*frozen*'
                          Else
                            Write(rfmt,'(a,i0,a)') '(2x,i10,a8,4i10,',dihedral%max_param,'f15.6)'
                            Write(message,rfmt) idih,keyword,dihedral%lst(1:4,ndihed), &
                              dihedral%param(1:dihedral%max_param,ndihed)
                          End If
                          Call info(message,.true.)
                       End If

! convert energies to internal units and angles to radians

                       dihedral%param(1,ndihed)=dihedral%param(1,ndihed)*engunit

                       If (dihedral%key(ndihed) == 4) Then
                          dihedral%param(2,ndihed)=dihedral%param(2,ndihed)*engunit
                          dihedral%param(3,ndihed)=dihedral%param(3,ndihed)*engunit
                       Else If (dihedral%key(ndihed) == 7) Then
                          dihedral%param(2,ndihed)=dihedral%param(2,ndihed)*engunit
                          dihedral%param(3,ndihed)=dihedral%param(3,ndihed)*engunit
                          dihedral%param(6,ndihed)=dihedral%param(6,ndihed)*engunit
                          dihedral%param(7,ndihed)=dihedral%param(7,ndihed)*(pi/180.0_wp)
                       Else
                          dihedral%param(2,ndihed)=dihedral%param(2,ndihed)*(pi/180.0_wp)
                       End If

                    Else ! TABDIH to read

! Construct unique name for the tabulated dihedral

                       Do jsite=1,ntpatm
                          If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                          If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                          If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                          If (sitnam(isite4) == unqatm(jsite)) katom4=jsite
                       End Do

                       If      (katom1 == katom4) Then
                          If (katom2 <= katom3) Then
                             iddihd = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)//sitnam(iatm4)
                          Else
                             iddihd = sitnam(iatm1)//sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm4)
                          End If
                       Else If (katom1 <  katom4) Then
                          iddihd = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)//sitnam(iatm4)
                       Else
                          iddihd = sitnam(iatm4)//sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm1)
                       End If

! ntpdih total number of unique table potentials to read from TABDIH

                       Do i=1,ntpdih
                          If (dihd_name(i) == iddihd) Then
                             dihedral%ltp(ndihed)=i ! Re-point from zero to type
                             Exit
                          End If
                       End Do

                       If (dihedral%ltp(ndihed) == 0) Then
                          ntpdih=ntpdih+1
                          dihd_name(ntpdih)=iddihd

                          dihedral%ltp(0)=ntpdih      ! NUTDP
                          dihedral%ltp(ndihed)=ntpdih ! Re-point from zero to type
                       End If

! test for frozen atoms and print unit

                       If (l_top) Then
                          If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) /= 0) Then
                            Write(message,'(2x,i10,a8,4i10,2x,a9,2x,a8)') &
                              idih,keyword,dihedral%lst(1:4,ndihed),'tabulated','*frozen*'
                          Else
                            Write(message,'(2x,i10,a8,4i10,2x,a9)') &
                              idih,keyword,dihedral%lst(1:4,ndihed),'tabulated'
                          End If
                          Call info(message,.true.)
                       End If

                    End If

! catch unidentified entry

                    If (Any(dihedral%lst(1:4,ndihed) < 1) .or. Any(dihedral%lst(1:4,ndihed) > numsit(itmols))) Call error(27)

! test for mistyped dihedral unit

                    If ( iatm1 == iatm2 .or. iatm1 == iatm3 .or. &
                         iatm2 == iatm3 .or. iatm1 == iatm4 .or. &
                         iatm2 == iatm4 .or. iatm3 == iatm4 ) Call error(67)
                 End Do

! Check for multiple dihedral angle entries

                 Do i=ndihed-dihedral%num(itmols)+1,ndihed
                    is(0)=dihedral%key(i)
                    is(1)=dihedral%lst(1,i)
                    is(2)=dihedral%lst(2,i)
                    is(3)=dihedral%lst(3,i)
                    is(4)=dihedral%lst(4,i)

                    Do j=i+1,ndihed
                       js(0)=dihedral%key(j)
                       js(1)=dihedral%lst(1,j)
                       js(2)=dihedral%lst(2,j)
                       js(3)=dihedral%lst(3,j)
                       js(4)=dihedral%lst(4,j)

                       If ((js(1) == is(1) .and. js(2) == is(2) .and. &
                            js(3) == is(3) .and. js(4) == is(4)) .or. &
                           (js(1) == is(4) .and. js(2) == is(3) .and. &
                            js(3) == is(2) .and. js(4) == is(1))) Then
                          If (l_str .and. l_top) Call warning(440,Real(i,wp),Real(j,wp),0.0_wp)
!                          If (is(0) == js(0)) Call error(620)
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
                 inversion%num(itmols)=inversion%num(itmols)+ntmp

                 Write(message,'(a,i10)') 'number of inversion angles ',ntmp
                 Call info(message,.true.)
                 If (l_top) Then
                   Write(messages(1),'(a)') 'inversion angle details:'
                   Write(messages(2),'(8x,a4,5x,a3,4(5x,a5),7x,a7,8x,a5,8x,a6)') &
                    'unit','key','index','index','index','index','f-const', &
                    'angle','factor'
                  Call info(messages,2,.true.)
                End If

                 Do iinv=1,inversion%num(itmols)
                    ninver=ninver+1
                    If (ninver > inversion%max_types) Call error(73)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 2000
                       Call get_word(record,word)
                    End Do

! read type of inversion potential

                    Call lower_case(word)
                    keyword=word(1:4)

                    If      (keyword == 'tab' ) Then
                       inversion%key(ninver)=20
                    Else If (keyword == '-tab') Then
                       inversion%key(ninver)=-20
                    Else If (keyword == 'harm') Then
                       inversion%key(ninver)=1
                    Else If (keyword == '-hrm') Then
                       inversion%key(ninver)=-1
                    Else If (keyword == 'hcos') Then
                       inversion%key(ninver)=2
                    Else If (keyword == '-hcs') Then
                       inversion%key(ninver)=-2
                    Else If (keyword == 'plan') Then
                       inversion%key(ninver)=3
                    Else If (keyword == '-pln') Then
                       inversion%key(ninver)=-3
                    Else If (keyword == 'xpln') Then
                       inversion%key(ninver)=4
                    Else If (keyword == '-xpl') Then
                       inversion%key(ninver)=-4
                    Else If (keyword == 'calc') Then
                       inversion%key(ninver)=5
                    Else If (keyword == '-clc') Then
                       inversion%key(ninver)=-5
                    Else
                       Call info(keyword,.true.)
                       Call error(449)
                    End If

! read inversion atom indices

                    Call get_word(record,word)
                    iatm1=Nint(word_2_real(word)) ! central atom
                    Call get_word(record,word)
                    iatm2=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm3=Nint(word_2_real(word))
                    Call get_word(record,word)
                    iatm4=Nint(word_2_real(word))

                    inversion%lst(1,ninver)=iatm1
                    inversion%lst(2,ninver)=iatm2
                    inversion%lst(3,ninver)=iatm3
                    inversion%lst(4,ninver)=iatm4

                    isite1 = nsite - numsit(itmols) + iatm1
                    isite2 = nsite - numsit(itmols) + iatm2
                    isite3 = nsite - numsit(itmols) + iatm3
                    isite4 = nsite - numsit(itmols) + iatm4

                    If (inversion%key(ninver) /= 20) Then

                       Call get_word(record,word)
                       inversion%param(1,ninver)=word_2_real(word)
                       Call get_word(record,word)
                       inversion%param(2,ninver)=word_2_real(word)
                       Call get_word(record,word)
                       inversion%param(3,ninver)=word_2_real(word)

! test for frozen atoms and print unit

                       If (comm%idnode == 0 .and. l_top) Then
                          If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) /= 0) Then
                            Write(rfmt,'(a,i0,a)') '(2x,i10,a8,4i10,',inversion%max_param,'f15.6,2x,a8)'
                            Write(message,rfmt) iinv,keyword,inversion%lst(1:4,ninver), &
                              inversion%param(1:inversion%max_param,ninver),'*frozen*'
                          Else
                            Write(rfmt,'(a,i0,a)') '(2x,i10,a8,4i10,',inversion%max_param,'f15.6)'
                            Write(message,rfmt) iinv,keyword,inversion%lst(1:4,ninver), &
                              inversion%param(1:inversion%max_param,ninver)
                          End If
                          Call info(message,.true.)
                       End If

! convert energies to internal units and angles to radians

                       inversion%param(1,ninver)=inversion%param(1,ninver)*engunit

                       If (inversion%key(ninver) == 5) Then
                          inversion%param(2,ninver)=inversion%param(2,ninver)*engunit
                       Else
                          inversion%param(2,ninver)=inversion%param(2,ninver)*(pi/180.0_wp)
                          If (inversion%key(ninver) == 2) &
                             inversion%param(2,ninver)=Cos(inversion%param(2,ninver))
                       End If

                    Else ! TABINV to read

! Construct unique name for the tabulated inversions

                       Do jsite=1,ntpatm
                          If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                          If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                          If (sitnam(isite4) == unqatm(jsite)) katom4=jsite
                       End Do

                       If      (Min(katom2,katom3,katom4) == katom2) Then
                          If (katom3 <= katom4) Then
                             idinvr = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)//sitnam(iatm4)
                          Else
                             idinvr = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm4)//sitnam(iatm3)
                          End If
                       Else If (Min(katom2,katom3,katom4) == katom3) Then
                          If (katom2 <= katom4) Then
                             idinvr = sitnam(iatm1)//sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm4)
                          Else
                             idinvr = sitnam(iatm1)//sitnam(iatm3)//sitnam(iatm4)//sitnam(iatm2)
                          End If
                       Else
                          If (katom2 <= katom3) Then
                             idinvr = sitnam(iatm1)//sitnam(iatm4)//sitnam(iatm2)//sitnam(iatm3)
                          Else
                             idinvr = sitnam(iatm1)//sitnam(iatm4)//sitnam(iatm3)//sitnam(iatm2)
                          End If
                       End If

! ntpinv total number of unique table potentials to read from TABINV

                       Do i=1,ntpinv
                          If (invr_name(i) == idinvr) Then
                             inversion%ltp(ninver)=i ! Re-point from zero to type
                             Exit
                          End If
                       End Do

                       If (inversion%ltp(ninver) == 0) Then
                          ntpinv=ntpinv+1
                          invr_name(ntpinv)=idinvr

                          inversion%ltp(0)=ntpinv      ! NUTIP
                          inversion%ltp(ninver)=ntpinv ! Re-point from zero to type
                       End If

! test for frozen atoms and print unit

                       If (comm%idnode == 0 .and. l_top) Then
                          If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) /= 0) Then
                            Write(message,'(2x,i10,a8,4i10,2x,a9,2x,a8)') &
                              iinv,keyword,inversion%lst(1:4,ninver),'tabulated','*frozen*'
                          Else
                            Write(message,'(2x,i10,a8,4i10,2x,a9)') &
                              iinv,keyword,inversion%lst(1:4,ninver),'tabulated'
                          End If
                          Call info(message,.true.)
                       End If

                    End If

! catch unidentified entry

                    If (Any(inversion%lst(1:4,ninver) < 1) .or. Any(inversion%lst(1:4,ninver) > numsit(itmols))) Call error(27)

! test for mistyped inversion unit

                    If ( iatm1 == iatm2 .or. iatm1 == iatm3 .or. &
                         iatm2 == iatm3 .or. iatm1 == iatm4 .or. &
                         iatm2 == iatm4 .or. iatm3 == iatm4 ) Call error(68)
                 End Do

! Check for multiple inversion angle entries

                 Do i=ninver-inversion%num(itmols)+1,ninver
                    is(0)=inversion%key(i)
                    is(1)=inversion%lst(1,i)
                    is(2)=inversion%lst(2,i)
                    is(3)=inversion%lst(3,i)
                    is(4)=inversion%lst(4,i)
                    Call shellsort(3,is(2:4))

                    Do j=i+1,ninver
                       js(0)=inversion%key(j)
                       js(1)=inversion%lst(1,j)
                       js(2)=inversion%lst(2,j)
                       js(3)=inversion%lst(3,j)
                       js(4)=inversion%lst(4,j)
                       Call shellsort(3,js(2:4))

                       If (js(1) == is(1) .and. js(2) == is(2) .and. &
                           js(3) == is(3) .and. js(4) == is(4)) Then
                          If (l_str .and. l_top) Call warning(450,Real(i,wp),Real(j,wp),0.0_wp)
                          If (is(0) == js(0)) Call error(620)
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

                 cons%megcon=cons%megcon+nummols(itmols)*(cons%numcon(itmols)-frzcon)
                 pmf%megpmf=pmf%megpmf+nummols(itmols)*pmf%numpmf(itmols)

                 megrgd=megrgd+nummols(itmols)*(numrgd(itmols)-frzrgd)

                 megtet=megtet+nummols(itmols)*tether%numteth(itmols)

                 bond%total=bond%total+nummols(itmols)*bond%num(itmols)
                 angle%total=angle%total+nummols(itmols)*angle%num(itmols)
                 dihedral%total=dihedral%total+nummols(itmols)*dihedral%num(itmols)
                 inversion%total=inversion%total+nummols(itmols)*inversion%num(itmols)

                 Go To 1000

              Else

! error exit for unidentified directive in molecular data

                 Call strip_blanks(record)
                 Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                 Call info(message,.true.)
                 Call error(12)

              End If

           End Do

! just finished with this type molecule data

1000       Continue

        End Do

! report total molecules and sites

        Write(messages(1),'(a,i10)') 'total number of molecules ',Sum(nummols(1:ntpmls))
        Write(messages(2),'(a,i10)') 'total number of sites ',nsite
        Call info(messages,2,.true.)

! Deal with intarmolecular potential tables:
! read & generate intramolecular potential & virial arrays

        If (bond%l_tab) Call bonds_table_read(bond_name,bond,comm)
        If (angle%l_tab) Call angles_table_read(angl_name,angle,comm)
        If (dihedral%l_tab) Call dihedrals_table_read(dihd_name,dihedral,comm)
        If (inversion%l_tab) Call inversions_table_read(invr_name,inversion,comm)

! If some intramolecular PDFs analysis is opted for

        If (mxgana > 0) Then

! Only for the requested types of PDFs (re)initialise:
! number of unique intramolecular PDFs and auxiliary identity arrays

           If (bond%bin_pdf > 0) Then
              ntpbnd = 0 ! for bonds
              bond_name = ' '
           End If
           If (angle%bin_adf > 0) Then
              ntpang = 0 ! for angles
              angl_name = ' '
           End If
           If (dihedral%bin_adf > 0) Then
              ntpdih = 0 ! for dihedrals
              dihd_name = ' '
           End If
           If (inversion%bin_adf > 0) Then
              ntpinv = 0 ! for inversions
              invr_name = ' '
           End If

           nsite =0
           nbonds=0
           nangle=0
           ndihed=0
           ninver=0
           Do itmols=1,ntpmls
              Do ibond=1,bond%num(itmols)*Merge(1,0,bond%bin_pdf > 0)
                 nbonds=nbonds+1

                 iatm1=bond%lst(1,nbonds)
                 iatm2=bond%lst(2,nbonds)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2

! Construct unique name for the bond

                 Do jsite=1,ntpatm
                    If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                    If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                 End Do

                 If (katom1 <= katom2) Then
                    idbond = sitnam(iatm1)//sitnam(iatm2)
                 Else
                    idbond = sitnam(iatm2)//sitnam(iatm1)
                 End If

! ntpbnd total number of unique BPDFs

                 Do i=1,ntpbnd
                    If (bond_name(i) == idbond) Then
                       bond%ldf(nbonds)=i ! Re-point from zero to type
                       Exit
                    End If
                 End Do

                 If (bond%ldf(nbonds) == 0) Then
                    ntpbnd=ntpbnd+1

                    bond_name(ntpbnd)=idbond
                    bond%ldf(0)=ntpbnd      ! NUTBPDF
                    bond%ldf(nbonds)=ntpbnd ! Re-point from zero to type
                 End If
              End Do

              Do iang=1,angle%num(itmols)*Merge(1,0,angle%bin_adf > 0)
                 nangle=nangle+1

                 iatm1=angle%lst(1,nangle)
                 iatm2=angle%lst(2,nangle)
                 iatm3=angle%lst(3,nangle)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2
                 isite3 = nsite + iatm3

! Construct unique name for the angle

                 Do jsite=1,ntpatm
                    If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                    If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                 End Do

                 If (katom1 <= katom3) Then
                    idangl = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)
                 Else
                    idangl = sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm1)
                 End If

! ntpang total number of unique APDFs

                 Do i=1,ntpang
                    If (angl_name(i) == idangl) Then
                       angle%ldf(nangle)=i ! Re-point from zero to type
                       Exit
                    End If
                 End Do

                 If (angle%ldf(nangle) == 0) Then
                    ntpang=ntpang+1
                    angl_name(ntpang)=idangl

                    angle%ldf(0)=ntpang      ! NUTAPDF
                    angle%ldf(nangle)=ntpang ! Re-point from zero to type
                 End If
              End Do

              Do idih=1,dihedral%num(itmols)*Merge(1,0,dihedral%bin_adf > 0)
                 ndihed=ndihed+1

                 iatm1=dihedral%lst(1,ndihed)
                 iatm2=dihedral%lst(2,ndihed)
                 iatm3=dihedral%lst(3,ndihed)
                 iatm4=dihedral%lst(4,ndihed)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2
                 isite3 = nsite + iatm3
                 isite4 = nsite + iatm4

! Construct unique name for the dihedral

                 Do jsite=1,ntpatm
                    If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                    If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                    If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                    If (sitnam(isite4) == unqatm(jsite)) katom4=jsite
                 End Do

                 If      (katom1 == katom4) Then
                    If (katom2 <= katom3) Then
                       iddihd = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)//sitnam(iatm4)
                    Else
                       iddihd = sitnam(iatm1)//sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm4)
                    End If
                 Else If (katom1 <  katom4) Then
                    iddihd = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)//sitnam(iatm4)
                 Else
                    iddihd = sitnam(iatm4)//sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm1)
                 End If

! ntpdih total number of unique DPDFs

                 Do i=1,ntpdih
                    If (dihd_name(i) == iddihd) Then
                       dihedral%ldf(ndihed)=i ! Re-point from zero to type
                       Exit
                    End If
                 End Do

                 If (dihedral%ldf(ndihed) == 0) Then
                    ntpdih=ntpdih+1
                    dihd_name(ntpdih)=iddihd

                    dihedral%ldf(0)=ntpdih      ! NUTDPDF
                    dihedral%ldf(ndihed)=ntpdih ! Re-point from zero to type
                 End If
              End Do

              Do iinv=1,inversion%num(itmols)*Merge(1,0,inversion%bin_adf > 0)
                 ninver=ninver+1

                 If (inversion%key(ninver) /= 5) Cycle ! avoid the calcite OoP potential

                 iatm1=inversion%lst(1,ninver)
                 iatm2=inversion%lst(2,ninver)
                 iatm3=inversion%lst(3,ninver)
                 iatm4=inversion%lst(4,ninver)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2
                 isite3 = nsite + iatm3
                 isite4 = nsite + iatm4

! Construct unique name for the tabulated inversions

                 Do jsite=1,ntpatm
                    If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                    If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                    If (sitnam(isite4) == unqatm(jsite)) katom4=jsite
                 End Do

                 If      (Min(katom2,katom3,katom4) == katom2) Then
                    If (katom3 <= katom4) Then
                       idinvr = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm3)//sitnam(iatm4)
                    Else
                       idinvr = sitnam(iatm1)//sitnam(iatm2)//sitnam(iatm4)//sitnam(iatm3)
                    End If
                 Else If (Min(katom2,katom3,katom4) == katom3) Then
                    If (katom2 <= katom4) Then
                       idinvr = sitnam(iatm1)//sitnam(iatm3)//sitnam(iatm2)//sitnam(iatm4)
                    Else
                       idinvr = sitnam(iatm1)//sitnam(iatm3)//sitnam(iatm4)//sitnam(iatm2)
                    End If
                 Else
                    If (katom2 <= katom3) Then
                       idinvr = sitnam(iatm1)//sitnam(iatm4)//sitnam(iatm2)//sitnam(iatm3)
                    Else
                       idinvr = sitnam(iatm1)//sitnam(iatm4)//sitnam(iatm3)//sitnam(iatm2)
                    End If
                 End If

! ntpinv total number of unique IPDFs

                 Do i=1,ntpinv
                    If (invr_name(i) == idinvr) Then
                       inversion%ldf(ninver)=i ! Re-point from zero to type
                       Exit
                    End If
                 End Do

                 If (inversion%ldf(ninver) == 0) Then
                    ntpinv=ntpinv+1
                    invr_name(ntpinv)=idinvr

                    inversion%ldf(0)=ntpinv      ! NUTIPDFs
                    inversion%ldf(ninver)=ntpinv ! Re-point from zero to type
                 End If
              End Do

              nsite=nsite+numsit(itmols)
           End Do

! Only for the requested types of PDFs (re)initialise number of unique intramolecular PDFs
! and allocate PDFs arrays and record species and presence(frozen and non-frozen)

           If (bond%bin_pdf > 0) Then
              ntpbnd = 0 ! for bonds
              Call allocate_bond_dst_arrays(bond) ! as it depends on bond%ldf(0)
!             bond%typ = 0 ! initialised in bonds_module
           End If
           If (angle%bin_adf > 0) Then
              ntpang = 0 ! for angles
              Call allocate_angl_dst_arrays(angle) ! as it depends on angle%ldf(0)
!             angle%typ = 0 ! initialised in angles_module
           End If
           If (dihedral%bin_adf > 0) Then
              ntpdih = 0 ! for dihedrals
              Call allocate_dihd_dst_arrays(dihedral) ! as it depends on dihedral%ldf(0)
!             dihedral%typ = 0 ! initialised in dihedrals
           End If
           If (inversion%bin_adf > 0) Then
              ntpinv = 0 ! for inversions
              Call allocate_invr_dst_arrays(inversion) ! as it depends on inversion%ldf(0)
!             inversion%typ = 0 ! initialised in inversions
           End If

           nsite =0
           nbonds=0
           nangle=0
           ndihed=0
           ninver=0
           Do itmols=1,ntpmls
              Do ibond=1,bond%num(itmols)*Merge(1,0,bond%bin_pdf > 0)
                 nbonds=nbonds+1

                 iatm1=bond%lst(1,nbonds)
                 iatm2=bond%lst(2,nbonds)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2

! ntpbnd total number of unique BPDFs

                 j=bond%ldf(nbonds)
                 If (j > ntpbnd) Then

! record species and presence(frozen and non-frozen)

                    ntpbnd=ntpbnd+1

                    Do jsite=1,ntpatm
                       If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                       If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                    End Do

                    If (katom1 <= katom2) Then
                       bond%typ(1,ntpbnd)=katom1
                       bond%typ(2,ntpbnd)=katom2
                    Else
                       bond%typ(1,ntpbnd)=katom2
                       bond%typ(2,ntpbnd)=katom1
                    End If

                    If (frzsit(isite1)*frzsit(isite2) == 0) Then
                       bond%typ(0,ntpbnd)=bond%typ(0,ntpbnd)+nummols(itmols)
                    Else
                       bond%typ(-1,ntpbnd)=bond%typ(-1,ntpbnd)+nummols(itmols)
                    End If

                 Else If (j > 0) Then

! accumulate the existing type and presence(frozen and non-frozen)

                    If (frzsit(isite1)*frzsit(isite2) == 0) Then
                       bond%typ(0,j)=bond%typ(0,j)+nummols(itmols)
                    Else
                       bond%typ(-1,j)=bond%typ(-1,j)+nummols(itmols)
                    End If

                 End If
              End Do

              Do iang=1,angle%num(itmols)*Merge(1,0,angle%bin_adf > 0)
                 nangle=nangle+1

                 iatm1=angle%lst(1,nangle)
                 iatm2=angle%lst(2,nangle)
                 iatm3=angle%lst(3,nangle)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2
                 isite3 = nsite + iatm3

                 j=angle%ldf(nangle)
                 If (j > ntpang) Then

! record species and presence(frozen and non-frozen)

                    ntpang=ntpang+1

                    Do jsite=1,ntpatm
                       If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                       If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                       If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                    End Do

                    angle%typ(2,ntpang)=katom2
                    If (katom1 <= katom3) Then
                       angle%typ(1,ntpang)=katom1
                       angle%typ(3,ntpang)=katom3
                    Else
                       angle%typ(1,ntpang)=katom3
                       angle%typ(3,ntpang)=katom1
                    End If

                    If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3) == 0) Then
                       angle%typ(0,ntpang)=angle%typ(0,ntpang)+nummols(itmols)
                    Else
                       angle%typ(-1,ntpang)=angle%typ(-1,ntpang)+nummols(itmols)
                    End If

                 Else If (j > 0) Then

! accumulate the existing type and presence(frozen and non-frozen)

                    If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3) == 0) Then
                       angle%typ(0,j)=angle%typ(0,j)+nummols(itmols)
                    Else
                       angle%typ(-1,j)=angle%typ(-1,j)+nummols(itmols)
                    End If

                 End If
              End Do

              Do idih=1,dihedral%num(itmols)*Merge(1,0,dihedral%bin_adf > 0)
                 ndihed=ndihed+1

                 iatm1=dihedral%lst(1,ndihed)
                 iatm2=dihedral%lst(2,ndihed)
                 iatm3=dihedral%lst(3,ndihed)
                 iatm4=dihedral%lst(4,ndihed)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2
                 isite3 = nsite + iatm3
                 isite4 = nsite + iatm4

                 j=dihedral%ldf(ndihed)
                 If (j > ntpdih) Then

! record species and presence(frozen and non-frozen)

                    ntpdih=ntpdih+1

                    Do jsite=1,ntpatm
                       If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                       If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                       If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                       If (sitnam(isite4) == unqatm(jsite)) katom4=jsite
                    End Do

                    If      (katom1 == katom4) Then
                       dihedral%typ(1,ntpdih)=katom1
                       dihedral%typ(4,ntpdih)=katom4
                       If (katom2 <= katom3) Then
                          dihedral%typ(2,ntpdih)=katom2
                          dihedral%typ(3,ntpdih)=katom3
                       Else
                          dihedral%typ(2,ntpdih)=katom3
                          dihedral%typ(3,ntpdih)=katom2
                       End If
                    Else If (katom1 <  katom4) Then
                       dihedral%typ(1,ntpdih)=katom1
                       dihedral%typ(2,ntpdih)=katom2
                       dihedral%typ(3,ntpdih)=katom3
                       dihedral%typ(4,ntpdih)=katom4
                    Else
                       dihedral%typ(1,ntpdih)=katom4
                       dihedral%typ(2,ntpdih)=katom3
                       dihedral%typ(3,ntpdih)=katom2
                       dihedral%typ(4,ntpdih)=katom1
                    End If

                    If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) == 0) Then
                       dihedral%typ(0,ntpdih)=dihedral%typ(0,ntpdih)+nummols(itmols)
                    Else
                       dihedral%typ(-1,ntpdih)=dihedral%typ(-1,ntpdih)+nummols(itmols)
                    End If

                 Else If (j > 0) Then

! accumulate the existing type and presence(frozen and non-frozen)

                    If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) == 0) Then
                       dihedral%typ(0,j)=dihedral%typ(0,j)+nummols(itmols)
                    Else
                       dihedral%typ(-1,j)=dihedral%typ(-1,j)+nummols(itmols)
                    End If
                 End If
              End Do

              Do iinv=1,inversion%num(itmols)*Merge(1,0,inversion%bin_adf > 0)
                 ninver=ninver+1

                 If (inversion%key(ninver) /= 5) Cycle ! avoid the calcite OoP potential

                 iatm1=inversion%lst(1,ninver)
                 iatm2=inversion%lst(2,ninver)
                 iatm3=inversion%lst(3,ninver)
                 iatm4=inversion%lst(4,ninver)

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2
                 isite3 = nsite + iatm3
                 isite4 = nsite + iatm4

                 j=inversion%ldf(ninver)
                 If (j > ntpinv) Then

! record species and presence(frozen and non-frozen)

                    ntpinv=ntpinv+1

                    Do jsite=1,ntpatm
                       If (sitnam(isite1) == unqatm(jsite)) katom1=jsite
                       If (sitnam(isite2) == unqatm(jsite)) katom2=jsite
                       If (sitnam(isite3) == unqatm(jsite)) katom3=jsite
                       If (sitnam(isite4) == unqatm(jsite)) katom4=jsite
                    End Do

                    inversion%typ(1,ntpinv)=katom1
                    If      (Min(katom2,katom3,katom4) == katom2) Then
                       If (katom3 <= katom4) Then
                          inversion%typ(2,ntpinv)=katom2
                          inversion%typ(3,ntpinv)=katom3
                          inversion%typ(4,ntpinv)=katom4
                       Else
                          inversion%typ(2,ntpinv)=katom2
                          inversion%typ(3,ntpinv)=katom4
                          inversion%typ(4,ntpinv)=katom3
                       End If
                    Else If (Min(katom2,katom3,katom4) == katom3) Then
                       If (katom2 <= katom4) Then
                          inversion%typ(2,ntpinv)=katom3
                          inversion%typ(3,ntpinv)=katom2
                          inversion%typ(4,ntpinv)=katom4
                       Else
                          inversion%typ(2,ntpinv)=katom3
                          inversion%typ(3,ntpinv)=katom4
                          inversion%typ(4,ntpinv)=katom2
                       End If
                    Else
                       If (katom2 <= katom3) Then
                          inversion%typ(2,ntpinv)=katom4
                          inversion%typ(3,ntpinv)=katom2
                          inversion%typ(4,ntpinv)=katom3
                       Else
                          inversion%typ(2,ntpinv)=katom4
                          inversion%typ(3,ntpinv)=katom3
                          inversion%typ(4,ntpinv)=katom2
                       End If
                    End If

                    If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) == 0) Then
                       inversion%typ(0,ntpinv)=inversion%typ(0,ntpinv)+nummols(itmols)
                    Else
                       inversion%typ(-1,ntpinv)=inversion%typ(-1,ntpinv)+nummols(itmols)
                    End If

                 Else If (j > 0) Then

! accumulate the existing type and presence(frozen and non-frozen)

                    If (frzsit(isite1)*frzsit(isite2)*frzsit(isite3)*frzsit(isite4) == 0) Then
                       inversion%typ(0,j)=inversion%typ(0,j)+nummols(itmols)
                    Else
                       inversion%typ(-1,j)=inversion%typ(-1,j)+nummols(itmols)
                    End If
                 End If
              End Do

              nsite=nsite+numsit(itmols)
           End Do

           mxtana = Max(ntpbnd*Merge(1,0,bond%bin_pdf > 0), &
                        ntpang*Merge(1,0,angle%bin_adf > 0), &
                        ntpdih*Merge(1,0,dihedral%bin_adf > 0), &
                        ntpinv*Merge(1,0,inversion%bin_adf > 0))
        End If

! Deallocate possibly allocated auxiliary intramolecular TPs/PDFs arrays

        If (bond%l_tab .or. bond%bin_pdf > 0) Deallocate (bond_name, Stat=fail(1))
        If (angle%l_tab .or. angle%bin_adf > 0) Deallocate (angl_name, Stat=fail(2))
        If (dihedral%l_tab .or. dihedral%bin_adf > 0) Deallocate (dihd_name, Stat=fail(3))
        If (inversion%l_tab .or. inversion%bin_adf > 0) Deallocate (invr_name, Stat=fail(4))
        If (Any(fail > 0)) Then
           Write(message,'(a)') 'read_field deallocation failure'
           Call error(0,message)
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
              Call info('adiabatic shell model in operation',.true.)
           Else
              If (lshl_all) Then
                 keyshl = 2
                 Call info('relaxed shell model in operation',.true.)
              Else
                 Call error(476)
              End If
           End If
        Else If (lecx) Then ! previously selected option for accounting for
           lecx = .false.   ! extended coulombic exclusion is abandoned
        End If

! Process MPOLES

        If (mximpl > 0) Call read_mpoles(l_top,sumchg,comm)

! check charmming shells (megshl) globalisation

        If (megshl > 0) Then
           nsite =0
           nshels=0

!          product             sum
           q_core_p = 1.0_wp ; q_core_s = 0.0_wp
           q_shel_p = 1.0_wp ; q_shel_s = 0.0_wp
           k_crsh_p = 1.0_wp ; k_crsh_s = 0.0_wp
           p_core_p = 1.0_wp ; p_core_s = 0.0_wp
           d_core_p = 1.0_wp ; d_core_s = 0.0_wp
           Do itmols=1,ntpmls
              Do ishls=1,numshl(itmols)
                 nshels=nshels+1
                 iatm1=lstshl(1,nshels) ! core
                 iatm2=lstshl(2,nshels) ! shell

                 isite1 = nsite + iatm1
                 isite2 = nsite + iatm2

                 q_core_p=q_core_p*chgsit(isite1)
                 q_core_s=q_core_s+Abs(chgsit(isite1))

                 q_shel_p=q_shel_p*chgsit(isite2)
                 q_shel_s=q_shel_s+Abs(chgsit(isite2))

                 k_crsh_p=k_crsh_p*prmshl(1,nshels)
                 k_crsh_s=k_crsh_s+prmshl(1,nshels)

                 If (mximpl > 0) Then
                    p_core_p=p_core_p*plrsit(isite1)
                    p_core_s=p_core_s+plrsit(isite1)

                    d_core_p=d_core_p*dmpsit(isite1)
                    d_core_s=d_core_s+dmpsit(isite1)
                 End If
              End Do
              nsite=nsite+numsit(itmols)
           End Do

! Checks for ABORTS

           lshl_abort=.false. ! no aborting

           If (Abs(q_core_p) <= zero_plus) Then
              lshl_abort=.true.
              Call warning('a core of a core-shell unit bears a zero charge',.true.)
           End If

           If (mximpl > 0) Then ! sort k_charm
              If (p_core_s <= zero_plus .and. &
                  k_crsh_s <= zero_plus .and. &
                  q_shel_s <= zero_plus) Then
                 lshl_abort=.true.
                 Call warning('core-shell units polarisability is compromised',.true.)
              Else
                 If (k_crsh_s <= zero_plus .and. q_shel_s <= zero_plus) Then
                    If (keyind > 0) Then
                       k_crsh_p = 1000 * eu_kcpm ! reset to k_charmm = 1000 kcal*mol^−1*Å^−2
                       smax=k_crsh_p             ! set smax
                       Call warning('all core-shell force constants and shell ' &
                         //'charges are zero force constants to use 1000' &
                         //'kcal*mol^−1*Å^−2 CHARMM default',.true.)
                    Else
                       lshl_abort=.true.
                       Call warning('all core-shell force constants and charges are zero',.true.)
                    End If
                 End If
              End If
           Else ! particles' polarisabilities not really needed
              If (k_crsh_p <= zero_plus) Then
                 lshl_abort=.true.
                 Call warning('a core-shell force constant is undefined/zero',.true.)
              End If
              If (Abs(q_core_p*q_shel_p) <= zero_plus) Then
                 lshl_abort=.true.
                 Call warning('core-shell units polarisability is compromised',.true.)
              End If
           End If

           If (mximpl > 0) Then ! Time for possible resets
              d_core_p=0.0_wp ! the new d_core_s check for switching charming off!!!

              nsite =0
              nshels=0

              Do itmols=1,ntpmls
                 Do ishls=1,numshl(itmols)
                    nshels=nshels+1
                    iatm1=lstshl(1,nshels) ! core
                    iatm2=lstshl(2,nshels) ! shell

                    isite1 = nsite + iatm1
                    isite2 = nsite + iatm2

                    q_core=chgsit(isite1)
                    p_core=plrsit(isite1)

                    q_shel=chgsit(isite2)
                    If (k_crsh_s <= zero_plus .and. q_shel_s <= zero_plus .and. &
                        k_crsh_p > zero_plus) Then
                       prmshl(1,nshels)=k_crsh_p
                       prmshl(2,nshels)=0.0_wp
                    End If
                    k_crsh=prmshl(1,nshels)

                    If (d_core_s <= zero_plus .and. thole >= -zero_plus) dmpsit(isite1) = thole
                    d_core=dmpsit(isite1)
                    d_core_p=d_core_p+d_core

                    If      (Abs(q_shel) <= zero_plus) Then ! set shell's charge
                       charge=-Sign(1.0_wp,q_core)*Sqrt(p_core*k_crsh/(r4pie0/epsq))
                       If (Abs(charge) <= zero_plus) Then
                          lshl_abort=.true.
                          Call warning(296,Real(ishls,wp),Real(itmols,wp),0.0_wp)
                       Else
                          chgsit(isite2)=charge
                          mpllfr(1,isite2)=charge
                          sumchg=sumchg+Abs(charge)
                       End If
                    Else If (p_core <= zero_plus) Then ! set drude force constants
                       If (k_crsh <= zero_plus) Then
                          lshl_abort=.true.
                          Call warning(296,Real(ishls,wp),Real(itmols,wp),0.0_wp)
                       Else
                          plrsit(isite1)=(r4pie0/epsq)*q_shel**2/k_crsh
                       End If
                    Else If (k_crsh <= zero_plus) Then ! set polarisability
                       If (p_core <= zero_plus) Then
                          lshl_abort=.true.
                          Call warning(296,Real(ishls,wp),Real(itmols,wp),0.0_wp)
                       Else
                          prmshl(1,nshels)=(r4pie0/epsq)*q_shel**2/p_core
                          prmshl(2,nshels)=0.0_wp
                       End If
                    End If

! Redefine plrsit as reciprocal of polarisability

                    plrsit(isite1)=1.0_wp/plrsit(isite1)

! copy polarisability and dumping from cores to their shells

                    plrsit(isite2)=plrsit(isite1)
                    dmpsit(isite2)=dmpsit(isite1)
                 End Do
                 nsite=nsite+numsit(itmols)
              End Do

              If (keyind == 1 .and. d_core_p <= zero_plus) Then
                 keyind = 0
                 Call warning('CHARMM polarisation scheme deselected due to zero dumping factor',.true.)
              End If
           Else
              nsite =0
              nshels=0

              Do itmols=1,ntpmls
                 Do ishls=1,numshl(itmols)
                    nshels=nshels+1
                    iatm1=lstshl(1,nshels) ! core
                    iatm2=lstshl(2,nshels) ! shell

                    isite1 = nsite + iatm1
                    isite2 = nsite + iatm2

                    q_core=chgsit(isite1)
                    q_shel=chgsit(isite2)
                    k_crsh=prmshl(1,nshels)

                    If (Abs(q_core*q_shel*k_crsh) <= zero_plus) Then
                       lshl_abort=.true.
                       Call warning(296,Real(ishls,wp),Real(itmols,wp),0.0_wp)
                    End If
                 End Do
                 nsite=nsite+numsit(itmols)
              End Do
           End If

           If (lshl_abort) Call error(615)
        End If

! check (keyfce) for charges in the system

        If (keyfce /= 0 .and. Abs(sumchg) <= zero_plus) Then
           If (comm%idnode == 0) Call warning(4,sumchg,0.0_wp,0.0_wp)
           If (l_str) Then
              keyfce=0
              Call info('Electrostatics switched off!!!',.true.)
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

        If (Abs(sumchg) > 1.0e-6_wp .and. comm%idnode == 0) Call warning(5,sumchg,0.0_wp,0.0_wp)

! check (megshl) shell topological irregularities

        If (megshl > 0) Then
           nsite =0
           nshels=0
           nconst=0
           nrigid=0
           nteth =0
! bonds are allowed
           nangle=0
           ndihed=0
           ninver=0

           Do itmols=1,ntpmls
              Do ishls=1,numshl(itmols)
                 nshels=nshels+1
                 ia=lstshl(1,nshels) ! core
                 ja=lstshl(2,nshels) ! shell

! shells have no DoF, even if they are moved dynamically
! their DoFs don't contribute towards any dynamical properties

                 dofsit(nsite+ja)=-3.0_wp

! test for constrained, RBed and tethered shells

                 Do icnst=1,cons%numcon(itmols)
                    nconst=nconst+1

                    If (Any(cons%lstcon(1:2,nconst) == ja)) Then
                       Call warning(301,Real(ishls,wp),Real(icnst,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nconst=nconst-cons%numcon(itmols)

                 Do i=1,pmf%numpmf(itmols)
                    If (Any(pmf%lstpmf(1:pmf%mxtpmf(1),1) == ja) .or. Any(pmf%lstpmf(1:pmf%mxtpmf(2),2) == ja)) Then
                       Call warning(300,Real(ishls,wp),Real(i,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do

                 Do irgd=1,numrgd(itmols)
                    nrigid=nrigid+1

                    lrgd=lstrgd(0,nrigid)
                    If (Any(lstrgd(1:lrgd,nrigid) == ja)) Then
                       Call warning(302,Real(ishls,wp),Real(irgd,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nrigid=nrigid-numrgd(itmols)

                 Do iteth=1,tether%numteth(itmols)
                    nteth=nteth+1

                    If (tether%lsttet(nteth) == ja) Then
                       Call warning(303,Real(ishls,wp),Real(iteth,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nteth=nteth-tether%numteth(itmols)

! test for core-shell units fully overlapped on angles, dihedrals and inversions

                 Do iang=1,angle%num(itmols)
                    nangle=nangle+1

                    If (Any(angle%lst(1:3,nangle) == ia) .and. Any(angle%lst(1:3,nangle) == ja)) Then
                       Call warning(297,Real(ishls,wp),Real(iang,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) nangle=nangle-angle%num(itmols)

                 Do idih=1,dihedral%num(itmols)
                    ndihed=ndihed+1

                    If (Any(dihedral%lst(1:4,ndihed) == ia) .and. Any(dihedral%lst(1:4,ndihed) == ja)) Then
                       Call warning(298,Real(ishls,wp),Real(idih,wp),Real(itmols,wp))
                       Call error(99)
                    End If

! core-shell up the 1 and 4 members

                    If (lecx) Then ! dihedral%l_core_shell=.false. is the default in dihedrals
                       If (dihedral%lst(1,ndihed) == ia) Then
                          dihedral%l_core_shell=.true.
                          dihedral%lst(5,ndihed)=ja
                       End If
                       If (dihedral%lst(1,ndihed) == ja) Then
                          dihedral%l_core_shell=.true.
                          dihedral%lst(5,ndihed)=ia
                       End If

                       If (dihedral%lst(4,ndihed) == ia) Then
                          dihedral%l_core_shell=.true.
                          dihedral%lst(6,ndihed)=ja
                       End If
                       If (dihedral%lst(4,ndihed) == ja) Then
                          dihedral%l_core_shell=.true.
                          dihedral%lst(6,ndihed)=ia
                       End If
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) ndihed=ndihed-dihedral%num(itmols)

                 Do iinv=1,inversion%num(itmols)
                    ninver=ninver+1

                    If (Any(inversion%lst(1:4,ninver) == ia) .and. Any(inversion%lst(1:4,ninver) == ja)) Then
                       Call warning(299,Real(ishls,wp),Real(iinv,wp),Real(itmols,wp))
                       Call error(99)
                    End If
                 End Do
                 If (ishls /= numshl(itmols)) ninver=ninver-inversion%num(itmols)
              End Do
              nsite=nsite+numsit(itmols)
           End Do

! if core-shelling up has occurred to 1 or/and 4 members then
! default the unshelled cores of 1 or/and 4 to the corresponding 5 & 6

           If (dihedral%l_core_shell) Then
              ndihed=0 ! initialise unshelled units
              Do itmols=1,ntpmls
                 Do idih=1,dihedral%num(itmols)
                    ndihed=ndihed+1

                    If (dihedral%lst(5,ndihed) == 0) dihedral%lst(5,ndihed)=dihedral%lst(1,ndihed)
                    If (dihedral%lst(6,ndihed) == 0) dihedral%lst(6,ndihed)=dihedral%lst(4,ndihed)
                 End Do
              End Do
           End If
        End If

! RB particulars

        If (m_rgd > 0) Then

! test for constraint units on RB units

           If (cons%m_con > 0) Then
              nsite =0
              nconst=0
              nrigid=0
              Do itmols=1,ntpmls
                 Do icnst=1,cons%numcon(itmols)
                    nconst=nconst+1
                    iatm1=cons%lstcon(1,nconst)
                    iatm2=cons%lstcon(2,nconst)

                    Do irgd=1,numrgd(itmols)
                       nrigid=nrigid+1
                       lrgd=lstrgd(0,nrigid)
                       If (Any(lstrgd(1:lrgd,nrigid) == iatm1) .and. Any(lstrgd(1:lrgd,nrigid) == iatm2)) Then
                          Call warning(304,Real(icnst,wp),Real(irgd,wp),Real(itmols,wp))
                          Call error(97)
                       End If
                    End Do
                    If (icnst /= cons%numcon(itmols)) nrigid=nrigid-numrgd(itmols)
                 End Do
                 nsite=nsite+numsit(itmols)
              End Do
           End If

! test for PMF units on RB units

           If (pmf%megpmf > 0) Then
              nsite =0
              nrigid=0
              Do itmols=1,ntpmls
                 Do i=1,pmf%numpmf(itmols)
                    Do ipmf=1,2
                       Do jpmf=1,pmf%mxtpmf(ipmf)
                          iatm1=pmf%lstpmf(jpmf,ipmf)
                          isite1=nsite+iatm1

                          Do irgd=1,numrgd(itmols)
                             nrigid=nrigid+1

                             lrgd=lstrgd(0,nrigid)
                             If (Any(lstrgd(1:lrgd,nrigid) == iatm1)) Then
                                Call warning(295,Real(ipmf,wp),Real(irgd,wp),Real(itmols,wp))

                                Call error(93)
                             End If
                          End Do
                          If (i /= pmf%numpmf(itmols) .and. ipmf /= 2) nrigid=nrigid-numrgd(itmols)
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

! Globalise megrgd & rcut for RBs

           rgdmeg=megrgd
           rgdrct=rcut

        End If

! read in rdf pairs

     Else If (word(1:3) == 'rdf') Then

        Call get_word(record,word)
        ntprdf=Nint(word_2_real(word))

        Write(message,'(a,i10)') 'number of specified rdf look up pairs ',ntprdf
        Call info(message,.true.)
        If (l_top) Then
          write(message,'(8x,a4,2(2x,a6))') 'pair','atom 1','atom 2'
          Call info(message,.true.)
        End If

        If (ntprdf > mxrdf) Call error(107)

        Do itprdf=1,ntprdf

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom1=word(1:8)
           Call get_word(record,word)
           atom2=word(1:8)

           If (l_top) Then
             Write(message,"(2x,i10,2a8)") itprdf,atom1,atom2
             Call info(message,.true.)
           End If

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

        Write(message,'(a,i10)') 'number of specified vdw potentials ',ntpvdw
        Call info(message,.true.)
        If (l_top) Then
          Write(message,'(8x,a4,5x,a6,2x,a6,5x,a3,7x,a10)') &
            'pair','atom 1','atom 2','key','parameters'
          Call info(message,.true.)
        End If

        If (ntpvdw > mxvdw) Call error(80)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

        Do itpvdw=1,ntpvdw

           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
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
           Else If (keyword == 'mors') Then
              keypot=8
           Else If (keyword == 'wca' ) Then
              keypot=9
           Else If (keyword == 'dpd' ) Then
              keypot=10
           Else If (keyword == '14-7') Then
              keypot=11
           Else If (keyword == 'ljc') Then
              keypot=12
           Else If (keyword == 'mstw') Then
              keypot=13
           Else If (keyword == 'ryd') Then
              keypot=14
           Else If (keyword == 'zbl') Then
              keypot=15
           Else If (keyword == 'zbls') Then
              keypot=16
           Else If (keyword == 'zblb') Then
              keypot=17
           Else
              Call info(keyword,.true.)
              Call error(452)
           End If

           If (keypot == 0) Then
              If (thermo%key_dpd > 0) Then ! make sure thermo%gamdpd is read and reported for DPD
                 Call get_word(record,word)
                 parpot(1)=word_2_real(word)
                 If (l_top) Then
                   Write(message,'(2x,i10,5x,2a8,8x,f20.6,1x,a9)') &
                     itpvdw,atom1,atom2,parpot(1),'tabulated'
                   Call info(message,.true.)
                 End If
              Else
                 If (l_top) Then
                   Write(message,'(2x,i10,5x,2a8,1x,a9)') &
                     itpvdw,atom1,atom2,'tabulated'
                   Call info(message,.true.)
                 End If
              End If
           Else
              itmp=Merge(mxpvdw+1,mxpvdw,thermo%key_dpd > 0)
              Do i=1,itmp ! make sure thermo%gamdpd is read and reported for DPD
                 Call get_word(record,word)
                 parpot(i)=word_2_real(word)
              End Do
              If (l_top) Then
                Write(rfmt,'(a,i0,a)') '(2x,i10,5x,2a8,3x,a4,1x,',itmp,'f15.6)'
                Write(message,rfmt) itpvdw,atom1,atom2,keyword,parpot(1:itmp)
                Call info(message,.true.)
              End If

! convert energies to internal unit

              parpot(1) = parpot(1)*engunit

              If      (keypot == 1) Then
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
              Else If (keypot == 13) Then
                 parpot(4)=parpot(4)*engunit
              Else If (keypot == 14) Then
                 parpot(2)=parpot(2)*engunit
              Else If (keypot == 15) Then
                 parpot(1)=parpot(1)/engunit
              Else If (keypot == 16) Then
                 parpot(1)=parpot(1)/engunit
                 parpot(5)=parpot(5)*engunit
              Else If (keypot == 17) Then
                 parpot(1)=parpot(1)/engunit
                 parpot(5)=parpot(5)*engunit
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

           If (thermo%key_dpd > 0) Then ! store possible specification of DPD's gamma_ij
              If      (keypot ==  0) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(1))
              Else If (keypot ==  1) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(3))
              Else If (keypot ==  2) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(3))
              Else If (keypot ==  3) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(5))
              Else If (keypot ==  4) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(4))
              Else If (keypot ==  5) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(6))
              Else If (keypot ==  6) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(3))
              Else If (keypot ==  7) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(6))
              Else If (keypot ==  8) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(4))
              Else If (keypot ==  9) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(4))
              Else If (keypot == 10) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(3))
              Else If (keypot == 11) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(3))
              Else If (keypot == 12) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(4))
              Else If (keypot == 13) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(5))
              Else If (keypot == 14) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(4))
              Else If (keypot == 15) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(5))
              Else If (keypot == 16) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(8))
              Else If (keypot == 17) Then
                 thermo%gamdpd(keyvdw)=Abs(parpot(8))
              End If
              If (thermo%gamdpd(0) > zero_plus) thermo%gamdpd(keyvdw)=thermo%gamdpd(0) ! override
           End If
        End Do

        If (ntpvdw > 0) Then

! test for unspecified atom-atom potentials

           ntab=(ntpatm*(ntpatm+1))/2
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

              If (thermo%key_dpd > 0) Then
                 If (All(thermo%gamdpd(1:mxvdw) <= zero_plus)) Then ! So thermo%gamdpd(0) <= zero_plus too
                    thermo%key_dpd = 0
                    Call info( &
                      'Ensemble NVT dpd defaulting to NVE (Microcanonical) ' &
                      //'due to all drag coefficients equal to zero',.true.)

                 Else
                    If (thermo%gamdpd(0) > zero_plus) Then
                      Call warning('all defined interactions have their drag coefficient overridden',.true.)
                    End If

                    If (mxtvdw == 0) Then
                      Call info('vdw/dpd cross terms mixing (for undefined mixed potentials) may be required',.true.)

                       If (thermo%gamdpd(0) > zero_plus .and. (.not.l_str)) Then
                          mxtvdw = 1
                          Write(messages(1),'(a)') &
                            'type of mixing defaulted - Lorentz–Berthelot :: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2'
                          Write(messages(2),'(a)') &
                            'mixing is limited to potentials of the same type only'
                          Write(messages(3),'(a)') &
                            'mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,14-7,LJC)'
                          Call info(messages,3,.true.)
                       End If
                    End If
                 End If
              End If

! If the user opted for possible vdw potential mixing or when required by DPD thermostat

              If (mxtvdw > 0) Then

                If (l_top .or. thermo%key_dpd > 0) Then
                  If (thermo%key_dpd > 0) Then
                    Call info('vdw potential mixing under testing...',.true.)
                  Else
                    Call info('dpd potential mixing under testing...',.true.)
                  End If
                End If

! Detect if there are qualifying candidates

                 ldpd_safe = .true. ! safe flag for DPD thermostats
                 nsite=0 ! number of new cross pair potentials
                 Do i=1,ntpatm
                    isite=(i*(i-1))/2+i
                    If (lstvdw(isite) <= ntpvdw) Then ! if it exists
                       ia=ltpvdw(lstvdw(isite))
                       Do j=i+1,ntpatm
                          jsite=(j*(j-1))/2+j
                          If (lstvdw(jsite) <= ntpvdw) Then ! if it exists
                             ja=ltpvdw(lstvdw(jsite))
                             If (ia == ja .and.               & ! only if of the same type
                                 (ia == 1  .or. ia == 2  .or. & ! and the type is allowed mixing
                                  ia == 9  .or. ia == 10 .or. & ! LJ, 12-6, WCA, DPD, 14-7, LJC
                                  ia == 11 .or. ia == 12)) Then
                                ksite=isite+j-i
                                If (lstvdw(ksite) > ntpvdw) Then ! if it does not exist - no overriding
                                   nsite=nsite+1
                                   lstvdw(ksite)=-1 ! set a temporary qualifier flag
                                End If
                             Else
                                If (lstvdw(ksite) > ntpvdw) Then
                                  If (thermo%key_dpd > 0) Then
                                    If (thermo%gamdpd(0) <= zero_plus) Then
                                      If (l_str) Then
                                        ldpd_safe = .false. ! test for non-definable interactions
                                        Call warning('the interaction between bead types: ' &
                                          //unqatm(i)//' & '//unqatm(j) &
                                          //' is unresolved and thus thermostating' &
                                          //' is ill defined in a DPD context', &
                                          .true.)
                                      Else
                                        Call warning('the interaction between bead types: ' &
                                          //unqatm(i)//' & '//unqatm(j) &
                                          //' is unresolved and thus thermostating is ill ' &
                                          //'defined in a DPD context but may be OK in CG MD', &
                                          .true.)
                                      End If
                                    End If
                                  Else
                                    Call warning('the interaction between atom types: ' &
                                      //unqatm(i)//' & '//unqatm(j) &
                                      //' is unresolved and thus this ' &
                                      //'cross-interaction will be left undefined', &
                                      .true.)
                                  End If
                                End If
                             End If
                          End If
                       End Do
                    End If
                 End Do
                 If (.not.ldpd_safe) Call error(512)

! Qualification has happened

                 If (nsite > 0) Then

                   If (l_top .or. thermo%key_dpd > 0) Then
                     If (thermo%key_dpd > 0) Then
                       Call info('vdw potential mixing underway...',.true.)
                     Else
                       Call info('dpd potential mixing underway...',.true.)
                     End If
                   End If

! As the range of defined potentials must extend
! put undefined potentials outside the new range

                    Do i=1,ntab
                       If (lstvdw(i) == ntpvdw+1) lstvdw(i) = lstvdw(i) + nsite
                    End Do

! Apply mixing

                    Do i=1,ntpatm
                       isite=(i*(i-1))/2+i
                       ia=lstvdw(isite)
                       keypot=ltpvdw(ia)
                       Do j=i+1,ntpatm
                          jsite=(j*(j-1))/2+j
                          ja=lstvdw(jsite)
                          ksite=isite+j-i
                          If (lstvdw(ksite) == -1) Then ! filter for action
                             ntpvdw = ntpvdw + 1        ! extend range
                             lstvdw(ksite)=ntpvdw       ! connect
                             ltpvdw(ntpvdw)=keypot

! Get mixing in LJ's style characteristic energy(EPSILON) & distance(SIGMA) terms

                             eps = 0.0_wp ; sig = 0.0_wp ; del = 0.0_wp
                             If      (keypot == 1)  Then ! 12-6
                                keyword='12-6'

                                eps(1)=prmvdw(2,ia)**2/(4.0_wp*prmvdw(1,ia))
                                sig(1)=(prmvdw(1,ia)/prmvdw(2,ia))**(1.0_wp/6.0_wp)

                                eps(2)=prmvdw(2,ja)**2/(4.0_wp*prmvdw(1,ja))
                                sig(2)=(prmvdw(1,ja)/prmvdw(2,ja))**(1.0_wp/6.0_wp)
                             Else If (keypot == 2  .or. keypot == 10 .or. &
                                      keypot == 11 .or. keypot == 12) Then ! LJ, DPD, 14-7, LJC
                                If (keypot == 2 ) keyword='lj  '
                                If (keypot == 10) keyword='dpd '
                                If (keypot == 11) keyword='14-7'
                                If (keypot == 12) keyword='ljc '

                                eps(1)=prmvdw(1,ia)
                                sig(1)=prmvdw(2,ia)

                                eps(2)=prmvdw(1,ja)
                                sig(2)=prmvdw(2,ja)
                             Else If (keypot == 9) Then ! WCA
                                keyword='wca '

                                eps(1)=prmvdw(1,ia)
                                sig(1)=prmvdw(2,ia)
                                del(1)=prmvdw(3,ia)

                                eps(2)=prmvdw(1,ja)
                                sig(2)=prmvdw(2,ja)
                                del(2)=prmvdw(3,ja)
                             End If

                             If      (mxtvdw == 1) Then

! Lorentz–Berthelot: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2

                                eps(0) = Sqrt(eps(1)*eps(2))

                                sig(0) = 0.5_wp*(sig(1)+sig(2))

                                If (Any(del > zero_plus)) &
                                del(0) = 0.5_wp*(del(1)+del(2))

                                If (thermo%key_dpd > 0) &
                                thermo%gamdpd(ksite) = Sqrt(thermo%gamdpd(isite)*thermo%gamdpd(jsite))

                             Else If (mxtvdw == 2) Then

! Fender-Halsey : e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2

                                eps(0) = 2.0_wp*eps(1)*eps(2) / (eps(1)+eps(2))

                                sig(0) = 0.5_wp*(sig(1)+sig(2))

                                If (Any(del > zero_plus)) &
                                del(0) = 0.5_wp*(del(1)+del(2))

                              If (thermo%key_dpd > 0) Then
                                If (thermo%gamdpd(isite)+thermo%gamdpd(jsite) > zero_plus) &
                                  thermo%gamdpd(ksite) = &
                                  2.0_wp*thermo%gamdpd(isite)*thermo%gamdpd(jsite) / &
                                  (thermo%gamdpd(isite)+thermo%gamdpd(jsite))
                              End If

                             Else If (mxtvdw == 3) Then

! Hogervorst good hope : e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)

                                eps(0) = Sqrt(eps(1)*eps(2))

                                sig(0) = Sqrt(sig(1)*sig(2))

                                If (Any(del > zero_plus)) &
                                del(0) = Sqrt(del(1)*del(2))

                                If (thermo%key_dpd > 0) &
                                thermo%gamdpd(ksite) = Sqrt(thermo%gamdpd(isite)*thermo%gamdpd(jsite))

                             Else If (mxtvdw == 4) Then

! Halgren HHG: e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)

                                eps(0) = 4.0_wp*eps(1)*eps(2) / (Sqrt(eps(1))+Sqrt(eps(2)))**2

                                sig(0) = (sig(1)**3+sig(2)**3) / (sig(1)**2+sig(2)**2)

                                If (Any(del > zero_plus)) &
                                del(0) = (del(1)**3+del(2)**3) / (del(1)**2+del(2)**2)

                              If (thermo%key_dpd > 0) Then
                                If (thermo%gamdpd(isite) >= zero_plus .and. thermo%gamdpd(jsite) >= zero_plus) Then
                                  If (Sqrt(thermo%gamdpd(isite))+Sqrt(thermo%gamdpd(jsite)) > zero_plus) Then
                                    thermo%gamdpd(ksite) = &
                                      4.0_wp*thermo%gamdpd(isite)*thermo%gamdpd(jsite) / &
                                      (Sqrt(thermo%gamdpd(isite))+Sqrt(thermo%gamdpd(jsite)))**2
                                  End If
                                End If
                              End If

                             Else If (mxtvdw == 5) Then

! Waldman–Hagler : e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ; s_ij=[(s_i^6+s_j^6)/2]^(1/6)

                                tmp    = 0.5_wp*(sig(1)**6+sig(2)**6)

                                eps(0) = Sqrt(eps(1)*eps(2)) * ((sig(1)*sig(2))**3) / tmp

                                sig(0) = tmp**(1.0_wp/6.0_wp)

                                If (Any(del > zero_plus)) &
                                del(0) = (0.5_wp*(del(1)**6+del(2)**6))**(1.0_wp/6.0_wp)

                                If (thermo%key_dpd > 0) &
                                thermo%gamdpd(ksite) = Sqrt(thermo%gamdpd(isite)*thermo%gamdpd(jsite)) * ((sig(1)*sig(2))**3) / tmp

                             Else If (mxtvdw == 6) Then

! Tang-Toennies : e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13 ;
!                 s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)

                                tmp    = (eps(1)*sig(1)**6) * (eps(2)*sig(2)**6)

                                eps(0) = tmp / ( ((eps(1)*sig(1)**12)**(1.0_wp/13.0_wp) +             &
                                                  (eps(2)*sig(2)**12)**(1.0_wp/13.0_wp)) * 0.5_wp )**13

                                sig(0) = (Sqrt(tmp)/eps(0))**(1.0_wp/6.0_wp)

                                If (Any(del > zero_plus)) &
                                del(0) = 0.5_wp*sig(0)*(del(1)/sig(1) + del(2)/sig(2))

                                If (thermo%key_dpd > 0) &
                                thermo%gamdpd(ksite) = 0.5_wp*eps(0)*(thermo%gamdpd(isite)/eps(1) + thermo%gamdpd(jsite)/eps(2))

                             Else If (mxtvdw == 7) Then

! Functional : e_ij=3 * (e_i*e_j)^(1/2) * (s_i*s_j)^3 / SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(6/(6-2L))} ;
!              s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}

                                eps(0) = 0.0_wp ; sig(0) = 0.0_wp
                                Do itmp=0,2
                                   tmp = (sig(1)**3+sig(2)**3)**2 / (4.0_wp*(sig(1)*sig(2))**itmp)

                                   eps(0) = eps(0) + tmp**(Real(6,wp)/Real(6-2*itmp,wp))

                                   sig(0) = sig(0) + tmp**(Real(1,wp)/Real(6-2*itmp,wp))
                                End Do
                                tmp = 1.0_wp / eps(0)

                                eps(0) = 3.0_wp * Sqrt(eps(1)*eps(2)) * (sig(1)*sig(2))**3 * tmp
                                sig(0)=sig(0) / 3.0_wp

                                If (Any(del > zero_plus)) &
                                del(0) = 0.5_wp*sig(0)*(del(1)/sig(1) + del(2)/sig(2))

                                If (thermo%key_dpd > 0) &
                                thermo%gamdpd(ksite) = 0.5_wp*eps(0)*(thermo%gamdpd(isite)/eps(1) + thermo%gamdpd(jsite)/eps(2))

                             End If

! Recover and/or paste in the vdw parameter array

                             If      (keypot == 1)  Then ! 12-6
                                prmvdw(1,ntpvdw)=4.0_wp*eps(0)*(sig(0)**12)
                                prmvdw(2,ntpvdw)=4.0_wp*eps(0)*(sig(0)**6)
                             Else If (keypot == 2  .or. keypot == 10 .or. &
                                      keypot == 11 .or. keypot == 12) Then ! LJ, DPD, 14-7, LJC
                                prmvdw(1,ntpvdw)=eps(0)
                                prmvdw(2,ntpvdw)=sig(0)
                             Else If (keypot == 9) Then ! WCA
                                prmvdw(1,ntpvdw)=eps(0)
                                prmvdw(2,ntpvdw)=sig(0)
                                prmvdw(3,ntpvdw)=del(0)
                             End If

                             If (l_top) Then
                               If (thermo%key_dpd > 0) Then
                                 Write(rfmt,'(a,i0,a)') '(2x,i10,5x,2a8,3x,a4,1x,',mxpvdw+1,'f20.6)'
                                 Write(message,rfmt) ntpvdw,unqatm(i),unqatm(j),keyword,parpot(1:mxpvdw+1)
                               Else
                                 Write(rfmt,'(a,i0,a)') '(2x,i10,5x,2a8,3x,a4,1x,',mxpvdw,'f20.6)'
                                 Write(message,rfmt) ntpvdw,unqatm(i),unqatm(j),keyword,parpot(1:mxpvdw)
                               End If
                               Call info(message,.true.)
                             End If

                          End If
                       End Do
                    End Do
                  Else
                    If (l_top) Then
                      Call info('vdw potential mixing unsuccessful or abandoned',.true.)
                    End If
                  End If

              End If

           End If

           If (thermo%key_dpd > 0) Then
             If      (All(thermo%gamdpd(1:mxvdw) <= zero_plus)) Then
               thermo%key_dpd = 0

               Call info('Ensemble NVT dpd defaulting to NVE (Microcanonical)' &
                 //'due to all drag coefficients equal to zero',.true.)
             Else If (Any(thermo%gamdpd(1:mxvdw) <= zero_plus)) Then ! in principle we should come up with the error before here
               Call warning('there is a two-body interaction with a' &
                 //'non-zero mutual drag coefficient',.true.)
               sigdpd(1:mxvdw) = Sqrt(2.0_wp*boltz*thermo%temp*thermo%gamdpd(1:mxvdw)) ! define sigdpd
             Else
               sigdpd(1:mxvdw) = Sqrt(2.0_wp*boltz*thermo%temp*thermo%gamdpd(1:mxvdw)) ! define sigdpd
             End If
           End If

! generate vdw force arrays

           If (.not.l_n_v) Then
              If ((.not. ld_vdw) .or. lt_vdw) Call vdw_generate(rvdw)
              If (lt_vdw) Call vdw_table_read(rvdw,comm)
              If (ld_vdw .and. Any(ltpvdw(1:ntpvdw) > 0)) Call vdw_direct_fs_generate(rvdw)
           End If

        End If

! read in the metal potential energy parameters

     Else If (word(1:3) == 'met') Then

        Call get_word(record,word)
        met%n_potentials=Nint(word_2_real(word))

        Write(message,'(a,i10)') 'number of specified metal potentials ',met%n_potentials
        Call info(message,.true.)
        If (l_top) Then
          Write(message,'(8x,a4,5x,a6,2x,a6,5x,a3,5x,a10)') &
            'pair','atom 1','atom 2','key','parameters'
          Call info(message,.true.)
        End If

        If (met%n_potentials > mxmet) Call error(71)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

        lmet_safe=.true.
        Do itpmet=1,met%n_potentials

           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
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
              keypot=0 ! met%tab=1 set in scan_field
              lmet_safe=(met%tab == 1)
           Else If (keyword(1:4) == 'eeam' ) Then
              keypot=0 ! met%tab=2 set in scan_field
              lmet_safe=(met%tab == 2)
           Else If (keyword(1:4) == '2bea' ) Then
              keypot=0 ! met%tab=3 set in scan_field
              lmet_safe=(met%tab == 3)
           Else If (keyword(1:4) == '2bee' ) Then
              keypot=0 ! met%tab=4 set in scan_field
              lmet_safe=(met%tab == 4)
           Else If (keyword(1:4) == 'fnsc') Then
              keypot=1
           Else If (keyword(1:4) == 'exfs') Then
              keypot=2
           Else If (keyword(1:4) == 'stch') Then
              keypot=3
           Else If (keyword(1:4) == 'gupt') Then
              keypot=4
           Else If (keyword(1:4) == 'mbpc') Then
              keypot=5
           Else
              Call info(keyword,.true.)
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

              If (l_top) Then
                Write(rfmt,'(a,i0,a)') '(2x,i10,5x,2a8,3x,a4,1x,',mxpmet,'f15.6)'
                Write(message,rfmt) itpmet,atom1,atom2,keyword,parpot(1:mxpmet)
                Call info(message,.true.)
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

           keymet=(ka1*(ka1-1))/2+ka2

           If (keymet > mxmet) Call error(82)
           If (met%list(keymet) /= 0) Call error(141)

           met%list(keymet)=itpmet
           met%ltp(itpmet)=keypot

           If (itpmet > 1) Then
              If (keypot /= met%ltp(itpmet-1)) lmet_safe=.false.
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
                 met%prm(i,itpmet)=parpot(i)
              End Do
           End If

        End Do

        If (met%n_potentials /= 0) Then

! test metal potentials mix-up

           If (.not.lmet_safe) Call error(92)

! test for unspecified atom-atom potentials

           ntab=(ntpatm*(ntpatm+1))/2
           If (met%n_potentials < ntab) Then
              Call warning(120,0.0_wp,0.0_wp,0.0_wp)

              If (met%n_potentials > mxmet) Call error(71)

! put undefined potentials outside range

              Do i=1,ntab
                 If (met%list(i) == 0) met%list(i)=met%n_potentials+1
              End Do

              Do i=met%n_potentials+1,ntab
                 met%ltp(i) = -1
              End Do
           End If

! generate metal force arrays

           Call metal_generate_erf(met)
           If (.not.met%l_direct) Then
              Call allocate_metal_table_arrays(met)
              If (met%tab > 0) Then ! keypot == 0
                 Call metal_table_read(l_top,met,comm)
              Else ! If (met%tab == 0) Then
                 Call metal_generate(met)
              End If
           End If

        End If

! read in the tersoff potential energy parameters

     Else If (word(1:7) == 'tersoff') Then

        Call get_word(record,word)
        ntpter=Nint(word_2_real(word))

        Write(message,'(a,i10)') 'number of specified tersoff potentials ',ntpter
        Call info(message,.true.)
        If (l_top) Then
          Write(message,'(6x,a6,5x,a5,7x,a3,5x,a10)') &
            'number','atom','key','parameters'
          Call info(message,.true.)
        End If

        If (ntpter > mxter) Call error(72)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

        lter_safe=.true.
        Do itpter=1,ntpter
           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           atom0=word(1:8)

           Call get_word(record,word)
           Call lower_case(word)
           keyword=word(1:4)

           If      (keyword == 'ters') Then
              keypot=1
           Else If (keyword == 'kihs') Then
              keypot=2
           Else
              Call info(keyword,.true.)
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
           parpot(5 )=Abs(word_2_real(word)) ! R_i

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 2000
              Call get_word(record,word)
           End Do

           parpot(6 )=Abs(word_2_real(word))    ! S_i
           Call get_word(record,word)

           If      (keypot == 1) Then

              parpot(7 )=word_2_real(word)      ! beta_i
              Call get_word(record,word)
              parpot(8 )=word_2_real(word)      ! eta_i
              Call get_word(record,word)
              parpot(9 )=word_2_real(word)      ! c_i
              Call get_word(record,word)
              parpot(10)=word_2_real(word)      ! d_i
              Call get_word(record,word)
              parpot(11)=word_2_real(word)      ! h_i

              If (l_top) Then
                Write(messages(1),'(2x,i10,5x,a8,3x,a4,1x,5(1p,e13.4))') &
                  itpter,atom0,keyword,parpot(1: 5)
                Write(messages(2),'(33x,6(1p,e13.4))') parpot(6:11)
                Call info(messages,2,.true.)
              End If

           Else If (keypot == 2) Then

              parpot(7 )=word_2_real(word)      ! eta_i
              Call get_word(record,word)
              parpot(8 )=word_2_real(word)      ! delta_i
              Call get_word(record,word)
              parpot(9 )=word_2_real(word)      ! c1_i
              Call get_word(record,word)
              parpot(10)=word_2_real(word)      ! c2_i
              Call get_word(record,word)
              parpot(11)=word_2_real(word)      ! c3_i

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record,comm)
                 If (.not.safe) Go To 2000
                 Call get_word(record,word)
              End Do

              parpot(12)=word_2_real(word)      ! c4_i
              Call get_word(record,word)
              parpot(13)=word_2_real(word)      ! c5_i
              Call get_word(record,word)
              parpot(14)=word_2_real(word)      ! h_i
              Call get_word(record,word)
              parpot(15)=word_2_real(word)      ! alpha_i
              Call get_word(record,word)
              parpot(16)=word_2_real(word)      ! beta_i

              If (l_top) Then
                Write(messages(1),'(2x,i10,5x,a8,3x,a4,1x,5(1p,e13.4))') &
                  itpter,atom0,keyword,parpot(1:5)
                Write(messages(2),'(33x,6(1p,e13.4))') parpot(6:11)
                Write(rfmt,'(a,i0,a)') '(33x,',mxpter-11,'(1p,e13.4))'
                Write(messages(3),rfmt) parpot(12:mxpter)
                Call info(messages,3,.true.)
              End If
           End If

           katom0=0

           Do jtpatm=1,ntpatm
              If (atom0 == unqatm(jtpatm)) katom0=jtpatm
           End Do

           If (katom0 == 0) Call error(74)

! convert parameters to internal units

           parpot(1)=parpot(1)*engunit
           parpot(3)=parpot(3)*engunit

           If (lstter(katom0) > 0) Call error(76)

           lfrter(katom0)=.true.

           lstter(katom0)=itpter
           ltpter(itpter)=keypot

           If (itpter > 1) Then
              If (keypot /= ltpter(itpter-1)) lter_safe=.false.
           End If

! calculate max tersoff cutoff

           rcter=Max(rcter,parpot(6))

! store tersoff single atom potential parameters

           Do i=1,mxpter
              prmter(i,itpter)=parpot(i)
           End Do

        End Do

! test tersoff potentials mix-up and cutoff conditions

        If (ntpter > 0) Then
           If (.not.lter_safe) Call error(90) ! Now potter holds keypot globally
           If (rcter < 1.0e-6_wp) Call error(79)
           If (rcut < 2.0_wp*rcter) Call error(102)
        End If

! generate tersoff interpolation arrays

        If (ntpter > 0) Call tersoff_generate(rcter)

! start processing cross atom potential parameters

        If (keypot == 1) Then

          Write(message,'(a,i10)') 'number of tersoff cross terms ', (ntpter*(ntpter+1))/2
          Call info(message,.true.)
          If (l_top) Then
            Write(message,'(8x,a4,5x,a6,2x,a6,5x,a10)') &
              'pair','atom 1','atom 2','paramters'
            Call info(message,.true.)
          End If

           Do icross=1,(ntpter*(ntpter+1))/2

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record,comm)
                 If (.not.safe) Go To 2000
                 Call get_word(record,word)
              End Do

              atom1=word(1:8)
              Call get_word(record,word)
              atom2=word(1:8)

              Call get_word(record,word)
              parpot(1)=word_2_real(word,1.0_wp) ! chi_ij
              Call get_word(record,word)
              parpot(2)=word_2_real(word,1.0_wp) ! omega_ij

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

              If (l_top) Then
                Write(message,'(2x,i10,5x,2a8,1x,2f15.6)') &
                  icross,atom1,atom2,parpot(1:2)
                Call info(message,.true.)
              End If
           End Do
        End If

! read in the three-body potential energy parameters

     Else If (word(1:3) == 'tbp') Then

        Call get_word(record,word)
        threebody%ntptbp=Nint(word_2_real(word))

        Write(message,'(a,i10)') 'number of specified three-body potentials ',threebody%ntptbp
        Call info(message,.true.)
        If (l_top) Then
          Write(message,'(5x,a7,5x,a6,2(2x,a6),5x,a3,5x,a10)') &
            'triplet','atom 1','atom 2','atom 3','key','parameters'
          Call info(message,.true.)
        End If

        If (threebody%ntptbp > threebody%mxtbp) Call error(83)
        If (.not.lunits) Call error(6)
        If (.not.lmols) Call error(13)

! Initialise head of group atom to not participating in an interaction

        Do itbp=1,threebody%mxtbp,threebody%mx2tbp
           threebody%lsttbp(itbp)=-1
        End Do

        Do itptbp=1,threebody%ntptbp
           parpot=0.0_wp

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
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
              Call info(keyword,.true.)
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

           If (l_top) Then
             Write(rfmt,'(a,i0,a)') '(2x,i10,5x,3a8,3x,a4,1x,',threebody%mxptbp,'f15.6)'
             Write(message,rfmt) itptbp,atom1,atom0,atom2,keyword,parpot(1:threebody%mxptbp)
             Call info(message,.true.)
           End If

           katom0=0
           katom1=0
           katom2=0

           Do jtpatm=1,ntpatm
              If (atom0 == unqatm(jtpatm)) katom0=jtpatm
              If (atom1 == unqatm(jtpatm)) katom1=jtpatm
              If (atom2 == unqatm(jtpatm)) katom2=jtpatm
           End Do

           If (katom0 == 0 .or. katom1 == 0 .or. katom2 == 0) Call error(84)

           threebody%lfrtbp(katom0)=.true.
           threebody%lfrtbp(katom1)=.true.
           threebody%lfrtbp(katom2)=.true.

           ka1=Max(katom1,katom2)
           ka2=Min(katom1,katom2)

           keytbp=(ka1*(ka1-1))/2+ka2+(katom0-1)*threebody%mx2tbp

           If (keytbp > threebody%mxtbp) Call error(86)

! convert parameters to internal units and angles to radians

           parpot(1) = parpot(1)*engunit
           If (keypot /= 6) parpot(2) = parpot(2)*(pi/180.0_wp)

           If (threebody%lsttbp(keytbp) > 0) Call error(18)

           threebody%lsttbp(keytbp)=itptbp
           ktbp=threebody%mx2tbp*((keytbp-1)/threebody%mx2tbp)+1 ! == threebody%mx2tbp*(katom0-1)+1
           If (threebody%lsttbp(ktbp) < 0) threebody%lsttbp(ktbp)=0

           threebody%ltptbp(itptbp)=keypot

! calculate max three-body cutoff

           threebody%cutoff=Max(threebody%cutoff,parpot(5))
           If (parpot(5) < 1.0e-6_wp) Then
              threebody%rcttbp(itptbp)=threebody%cutoff
              parpot(5)=threebody%cutoff
           Else
              threebody%rcttbp(itptbp)=parpot(5)
           End If

! store three-body potential parameters

           Do i=1,threebody%mxptbp
              threebody%prmtbp(i,itptbp)=parpot(i)
           End Do
        End Do

        If (threebody%ntptbp > 0) Then
           If (threebody%cutoff < 1.0e-6_wp) Call error(451)
           If (rcut < 2.0_wp*threebody%cutoff) Call error(471)
        End If

! read in the four-body potential energy parameters

     Else If (word(1:3) == 'fbp') Then

        Call get_word(record,word)
        ntpfbp=Nint(word_2_real(word))

        Write(message,'(a,i10)') 'number of specified four-body potentials ',ntpfbp
        Call info(message,.true.)
        If (l_top) Then
          Write(message,'(5x,a7,5x,a6,3(2x,a6),5x,a3,5x,a10)') &
            'quartet','atom 1','atom 2','atom 3','atom 4','key','parameters'
          Call info(message,.true.)
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
              Call get_line(safe,nfield,record,comm)
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
              Call info(keyword,.true.)
              Call error(443)
           End If

           Call get_word(record,word)
           parpot(1)=word_2_real(word)
           Call get_word(record,word)
           parpot(2)=word_2_real(word)
           Call get_word(record,word)
           parpot(3)=word_2_real(word)

           If (l_top) Then
             Write(rfmt,'(a,i0,a)') '(2x,i10,3x,4a8,3x,a4,2x,',mxpfbp,'f15.6)'
             Write(message,rfmt) itpfbp,atom0,atom1,atom2,atom3,keyword,parpot(1:mxpfbp)
             Call info(message,.true.)
           End If

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

! read kim interaction data - kim and rkim set in scan_field

     Else If (word(1:3) == 'kim') Then

       Write(message,'(2a)') 'using open KIM interaction model: ',kimim
       Call info(message,.true.)

! read external field data

     Else If (word(1:6) == 'extern') Then

        Call get_word(record,word)
        nfld=0
        If (word(1:1) /= '#' .and. word(1:1) /= ' ') nfld=Nint(word_2_real(word))
        If (nfld <= 0) nfld=6

        word(1:1)='#'
        Do While (word(1:1) == '#' .or. word(1:1) == ' ')
           Call get_line(safe,nfield,record,comm)
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
        Else If (keyword == 'xpis') Then
           keyfld=8
           If (l_vom) Then
             Call info('"no vom" option auto-switched on - COM momentum removal will be abandoned',.true.)
             Call warning('this may lead to a build up of the COM momentum' &
               //'and a manifestation of the "flying ice-cube" effect',.true.)
              l_vom=.false. ! exclude COM momentum rescaling by default
           End If
        Else If (keyword == 'zres') Then
           keyfld=9
        Else If (keyword == 'zrs-') Then
           keyfld=10
        Else If (keyword == 'zrs+') Then
           keyfld=11
        Else If (keyword == 'osel') Then
           keyfld=12
        Else If (keyword == 'ushr') Then
           keyfld=13
        Else
           Call info(keyword,.true.)
           Call error(454)
        End If

        Do i=1,nfld
           Call get_word(record,word)
           prmfld(i)=word_2_real(word)
        End Do

        Write(message,'(2a)') 'external field key ',keyword
        Call info(message,.true.)
        If (l_top) Then
          Write(messages(1),'(2x,a)') 'parameters'
          Write(rfmt,'(a,i0,a)') '(2x,',mxpfld,'f15.6)'
          Write(messages(2),rfmt) prmfld(1:mxpfld)
          Call info(messages,2,.true.)
        End If

! convert to internal units

        If      (keyfld == 1 .or. keyfld == 4 .or. keyfld == 5 .or. keyfld == 12) Then

           If (.not.lunits) Call error(6)

           Do i=1,3
              prmfld(i) = prmfld(i)*engunit
           End Do

        Else If (keyfld == 2 .or. keyfld == 6 .or. keyfld == 7) Then

           If (.not.lunits) Call error(6)

           prmfld(1) = prmfld(1)*engunit

        Else If (keyfld == 8) Then

           prmfld(3) = prmfld(3)/prsunt ! piston pressure specified in k-atm
           prmfld(3) = prmfld(3)*cell(5)*cell(9) ! convert to force

        Else If (keyfld == 9 .or. keyfld == 10 .or. keyfld == 11) Then

           If (.not.lunits) Call error(6)

           prmfld(3) = prmfld(3)*engunit

        Else If (keyfld == 13) Then

           If (.not.lunits) Call error(6)

           prmfld(5) = prmfld(5)*engunit

        End If

        If ((keyfld == 2 .or. keyfld == 8) .and. (imcon /= 1 .and. imcon /= 2)) Then
          Call warning('external field is ignored as only applicable for imcon=1,2 (orthorhombic geometry)',.true.)
        End If
        If (keyfld == 3 .and. imcon /= 6) Then
          Call warning('external field is ignored as only applicable for imcon=6 (SLAB geometry)',.true.)
        End If

        If (keyfld == 8 .and. thermo%ensemble /= ENS_NVE) Call error(7)

! close force field file

     Else If (word(1:5) == 'close') Then

        If (comm%idnode == 0) Close(Unit=nfield)

! Precautions: (vdw,met) may have led to rdf scanning (mxrdf > 0), see set_bounds

        If (ntprdf == 0 .and. mxrdf > 0) Then
           Do ia=1,ntpatm
              Do ja=ia,ntpatm
                 keyrdf=(ja*(ja-1))/2+ia
                 i=0
                 If (ntpvdw > 0) i=Max(i,lstvdw(keyrdf))
                 If (met%n_potentials > 0) i=Max(i,met%list(keyrdf))
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

        Call dihedrals_14_check &
           (l_str,l_top,ntpmls,nummols,angle,dihedral,comm)

! test for existence/appliance of any two-body or tersoff or KIM model defined interactions!!!

        If ( keyfce == 0 .and. ntpvdw == 0 .and. &
             met%n_potentials == 0 .and. ntpter == 0 .and. kimim == ' ') Call error(145)

! test for mixing KIM model with external interactions

        If ( (keyfce /= 0 .or. ntpvdw /= 0 .or. &
              met%n_potentials /= 0 .or. ntpter /= 0) .and. kimim /= ' ') Then
          Call warning('open KIM model in use together with extra intermolecular interactions',.true.)
        End If

! EXIT IF ALL IS OK

        Return

     Else

! error exit for unidentified directive

        Call info(word(1:Len_Trim(word)),.true.)
        Call error(4)

     End If

  End Do

! uncontrolled error exit from field file processing

  If (comm%idnode == 0) Close(Unit=nfield)
  Call error(16)

! end of field file error exit

2000 Continue

  If (comm%idnode == 0) Close(Unit=nfield)
  Call error(52)

End Subroutine read_field

Subroutine report_topology               &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megrgd,  &
           megtet,  &
           cons,pmf, &
           bond,angle,dihedral,inversion,tether,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reporting FIELD topology
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer, Intent( In    ) :: megatm,megfrz,atmfre,atmfrz, &
                              megshl,megrgd, &
                              megtet
  Type( constraints_type ), Intent( In    ) :: cons
  Type( pmf_type ), Intent( In    ) :: pmf
  Type( bonds_type ), Intent( In    ) :: bond
  Type( angles_type ), Intent( In    ) :: angle
  Type( dihedrals_type ), Intent( In    ) :: dihedral
  Type( inversions_type ), Intent( In    ) :: inversion
  Type(tethers_type), Intent( InOut ) :: tether
  Type(comms_type), Intent( InOut ) :: comm

  Integer :: itmols,nsite,                &
             isite1,isite2,isite3,isite4, &
             iatm1,iatm2,iatm3,iatm4,     &
             ishls,nshels,frzshl,mgfrsh,  &
             mgcon,                       &
             lpmf,ipmf,jpmf,              &
             mgrgd,                       &
             iteth,nteth,frztet,mgfrtt,   &
             ibond,nbonds,frzbnd,mgfrbn,  &
             iang,nangle,frzang,mgfran,   &
             idih,ndihed,frzdih,mgfrdh,   &
             iinv,ninver,frzinv,mgfrin

  Character( Len = 4) :: frzpmf

  Character( Len = : ),Allocatable  :: fmt1,fmt2,fmt3
  Character( Len = 256) :: banner(18)

  nshels = 0
  mgfrsh = 0

  mgcon  = 0

  frzpmf='NONE'

  mgrgd  = 0

  nteth  = 0
  mgfrtt = 0

  nbonds = 0
  mgfrbn = 0

  nangle = 0
  mgfran = 0

  ndihed = 0
  mgfrdh = 0

  ninver = 0
  mgfrin = 0

  nsite  = 0
  Do itmols=1,ntpmls
     frzshl=0
     Do ishls=1,numshl(itmols)
        nshels=nshels+1

        iatm1=lstshl(1,nshels)

        isite1 = nsite + iatm1

        If (frzsit(isite1) == 1) frzshl=frzshl+1
     End Do

     Do lpmf=1,pmf%numpmf(itmols) ! pmf%numpmf can only be 1 or 0, so the 'Do' loop is used as an 'If' condition
        Do ipmf=1,2
           Do jpmf=1,pmf%mxtpmf(ipmf)
              iatm1=pmf%lstpmf(jpmf,ipmf)

              isite1 = nsite + iatm1

              If (frzsit(isite1) == 1) frzpmf='ALL'
           End Do
        End Do
     End Do

     frztet=0
     Do iteth=1,tether%numteth(itmols)
        nteth=nteth+1

        iatm1=tether%lsttet(nteth)

        isite1 = nsite + iatm1

        If (frzsit(isite1) == 1) frztet=frztet+1
     End Do

     frzbnd=0
     Do ibond=1,bond%num(itmols)
        nbonds=nbonds+1

        iatm1=bond%lst(1,nbonds)
        iatm2=bond%lst(2,nbonds)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2

        If (frzsit(isite1)+frzsit(isite2) == 2) frzbnd=frzbnd+1
     End Do

     frzang=0
     Do iang=1,angle%num(itmols)
        nangle=nangle+1

        iatm1=angle%lst(1,nangle)
        iatm2=angle%lst(2,nangle)
        iatm3=angle%lst(3,nangle)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2
        isite3 = nsite + iatm3

        If (frzsit(isite1)+frzsit(isite2)+frzsit(isite3) == 3) frzang=frzang+1
     End Do

     frzdih=0
     Do idih=1,dihedral%num(itmols)
        ndihed=ndihed+1

        iatm1=dihedral%lst(1,ndihed)
        iatm2=dihedral%lst(2,ndihed)
        iatm3=dihedral%lst(3,ndihed)
        iatm4=dihedral%lst(4,ndihed)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2
        isite3 = nsite + iatm3
        isite4 = nsite + iatm4

        If (frzsit(isite1)+frzsit(isite2)+frzsit(isite3)+frzsit(isite4) == 4) frzdih=frzdih+1
     End Do

     frzinv=0
     Do iinv=1,inversion%num(itmols)
        ninver=ninver+1

        iatm1=inversion%lst(1,ninver)
        iatm2=inversion%lst(2,ninver)
        iatm3=inversion%lst(3,ninver)
        iatm4=inversion%lst(4,ninver)

        isite1 = nsite + iatm1
        isite2 = nsite + iatm2
        isite3 = nsite + iatm3
        isite4 = nsite + iatm4

        If (frzsit(isite1)+frzsit(isite2)+frzsit(isite3)+frzsit(isite4) == 4) frzinv=frzinv+1
     End Do

     mgcon=mgcon+nummols(itmols)*cons%numcon(itmols)
     mgrgd=mgrgd+nummols(itmols)*numrgd(itmols)

     mgfrsh=mgfrsh+nummols(itmols)*frzshl
     mgfrtt=mgfrtt+nummols(itmols)*frztet
     mgfrbn=mgfrbn+nummols(itmols)*frzbnd
     mgfran=mgfran+nummols(itmols)*frzang
     mgfrdh=mgfrdh+nummols(itmols)*frzdih
     mgfrin=mgfrin+nummols(itmols)*frzinv

     nsite=nsite+numsit(itmols)
  End Do

  fmt1 = '(a66)'
  fmt2 = '(a29,i11,a8,i11,a7)'
  fmt3 = '(a29,i11,a15,a4,a7)'
  Write(banner(1),fmt1) '//'//Repeat('=',62)//'\\'
  Write(banner(2),fmt1) '||'//Repeat(' ',62)//'||'
  Write(banner(3),fmt1) '||           SUMMARY  OF  TOPOLOGICAL  DECOMPOSITION            ||'
  Write(banner(4),fmt1) '||'//Repeat(' ',62)//'||'
  Write(banner(5),fmt1) '||'//Repeat('-',62)//'||'
  Write(banner(6),fmt1) '||  INTERACTION  OR  TYPE  |  GRAND TOTAL | Fully/Partly FROZEN ||'
  Write(banner(6),fmt1) '||-------------------------+--------------+---------------------||'
  Write(banner(7),fmt2) '||  all particles/sites    | ',megatm,'  |  F  ',megfrz,'     ||'
  Write(banner(8),fmt2) '||  free particles         | ',atmfre,'  |  F  ',atmfrz,'     ||'
  Write(banner(9),fmt2) '||  core-shell units       | ',megshl,'  |  P  ',mgfrsh,'     ||'
  Write(banner(10),fmt2) '||  constraint bond units  | ',mgcon,'  |  F  ',mgcon-cons%megcon,'     ||'
  Write(banner(11),fmt3) '||  PMF units              | ',pmf%megpmf,'  |  P         ',frzpmf,'     ||'
  Write(banner(12),fmt2) '||  rigid body units       | ',mgrgd,'  |  F  ',mgrgd-megrgd,'     ||'
  Write(banner(13),fmt2) '||  tethered atom units    | ',megtet,'  |  F  ',mgfrtt,'     ||'
  Write(banner(14),fmt2) '||  chemical bond units    | ',bond%total,'  |  F  ',mgfrbn,'     ||'
  Write(banner(15),fmt2) '||  bond angle units       | ',angle%total,'  |  F  ',mgfran,'     ||'
  Write(banner(16),fmt2) '||  dihedral angle units   | ',dihedral%total,'  |  F  ',mgfrdh,'     ||'
  Write(banner(17),fmt2) '||  inversion angle units  | ',inversion%total,'  |  F  ',mgfrin,'     ||'
  Write(banner(18),fmt1) '\\'//Repeat('=',62)//'//'
  Call info(banner,18,.true.)
End Subroutine report_topology

Subroutine scan_field                                &
           (l_n_e,mxompl,mximpl,                     &
           mxsite,mxatyp,megatm,mxtmls,max_exclude,       &
           mtshl,mxtshl,mxshl,mxfshl,                &
           mtcons,              &
           l_usr,                &
           mtrgd,mxtrgd,mxrgd,mxlrgd,mxfrgd,         &
           mtteth,             &
           mtbond, &
           mtangl,       &
           mtdihd,       &
           mtinv,         &
           mxrdf,mxvdw,rvdw,mxgvdw,                  &
           mxmet,mxmed,mxmds,            &
           mxter,rcter,mxfbp,rcfbp,lext,cons,pmf,met,&
           bond,angle,dihedral,inversion,tether,threebody,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of the FIELD file
!
! copyright - daresbury laboratory
! author    - i.t.todorov & w.smith november 2016
! contrib   - b.palmer (2band) may 2013
! contrib   - a.v.brukhno & i.t.todorov march 2014 (itramolecular TPs)
! contrib   - h.a.boateng february 2015
! contrib   - a.m.elena february 2017
! contrib   - i.t.todorov may 2017
! contrib   - v.sokhan may 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( pmf_type ), Intent( InOut ) :: pmf
  Type( constraints_type ), Intent( InOut ) :: cons
  Type( metal_type ), Intent( InOut ) :: met
  Type( bonds_type ), Intent( InOut ) :: bond
  Type( angles_type ), Intent( InOut ) :: angle
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( inversions_type ), Intent( InOut ) :: inversion
  Type( tethers_type ), Intent( InOut ) :: tether
  Type( threebody_type ), Intent( InOut ) :: threebody
  Type( comms_type ), Intent( InOut ) :: comm
  Integer( Kind = wi ), Intent(   Out ) :: max_exclude
! Max number of different atom types

  Integer, Parameter :: mmk = 1000

! Average maximum number of intra-like bonds per atom

  Integer, Parameter :: mxb = 6

  Character( Len = 8   ), Dimension( 1:mmk ) :: chr

  Character( Len = 200 ) :: record,record_raw
  Character( Len = 40  ) :: word
  Character( Len = 8   ) :: name

  Logical           :: l_n_e,check,safe,l_usr,lext
  Integer           :: mxtmls,itmols,nummols,numsit,mxnmst,ksite,nrept,        &
                       mxompl,mximpl,mxsite,mxatyp,megatm,i,j,k,        &
                       numshl,mtshl,mxtshl,mxshl,ishls,mxfshl,                 &
                       numcon,mtcons,icon,                &
                       ipmf,jpmf,                     &
                       numrgd,mtrgd,mxtrgd,mxlrgd,mxrgd,irgd,jrgd,lrgd,mxfrgd, &
                       inumteth,mtteth,iteth,  &
                       numbonds,mtbond,ibonds, &
                       numang,mtangl,iang,         &
                       numdih,mtdihd,idih,         &
                       numinv,mtinv,iinv,           &
                       mxrdf,itprdf,mxvdw,itpvdw,mxgvdw,                       &
                       mxmet,mxmed,mxmds,itpmet,                        &
                       mxter,itpter,itptbp,mxfbp,itpfbp,                 &
                       mxt(1:9),mxf(1:9)
  Real( Kind = wp ) :: rvdw,rcter,rcfbp,rct,tmp,tmp1,tmp2

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
  cons%mxcons=0
  cons%mxtcon=0
  cons%mxfcon=0

  pmf%mxpmf =0
  pmf%mxtpmf=0
  pmf%mxfpmf=0

  numrgd=0
  mtrgd =0
  mxrgd =0
  mxtrgd=0
  mxlrgd=0
  mxfrgd=0

  inumteth=0
  mtteth =0
  tether%mxteth =0
  tether%mxtteth=0
  tether%mxftet =0

  numbonds=0
  mtbond=0
  bond%max_bonds=0
  bond%max_types=0
  bond%max_legend=0
  bond%rcut =0.0_wp
  bond%bin_tab=-2

  numang=0
  mtangl=0
  angle%max_angles=0
  angle%max_types=0
  angle%max_legend=0
  angle%bin_tab=-2

  numdih=0
  mtdihd=0
  dihedral%max_angles=0
  dihedral%max_types=0
  dihedral%max_legend=0
  dihedral%bin_tab=-2

  numinv=0
  mtinv =0
  inversion%max_angles =0
  inversion%max_types=0
  inversion%max_legend=0
  inversion%bin_tab=-2

  mxrdf =0

  mxvdw =0
  rvdw  =0.0_wp
  mxgvdw=0

  mxmet =0
  mxmed =0
  mxmds =0
  met%rcut  =0.0_wp
  met%maxgrid=0

  mxter=0
  rcter=0.0_wp

  threebody%mxtbp=0
  threebody%cutoff=0.0_wp

  mxfbp=0
  rcfbp=0.0_wp

  max_exclude=0

  lext =.false.
  l_usr=.false.

! Set safe flag

  safe=.true.

! Open the interactions input file

  If (comm%idnode == 0) Inquire(File=Trim(field), Exist=safe)
  Call gcheck(comm,safe,"enforce")
  If (.not.safe) Then
     Go To 20
  Else
     If (comm%idnode == 0) Open(Unit=nfield, File=Trim(field), Status='old')
  End If

  Call get_line(safe,nfield,record,comm)
  If (.not.safe) Go To 30

  Do

     word(1:1)='#'
     Do While (word(1:1) == '#' .or. word(1:1) == ' ')
        Call get_line(safe,nfield,record,comm)
        record_raw = record ! KIM needs case sensitive IM name
        If (.not.safe) Go To 30
        Call lower_case(record)
        Call get_word(record,word)
     End Do

! multipoles container detection

     If (word(1:5) == 'multi') Then

        Call get_word(record,word) ; Call lower_case(word)
        If (word(1:5) == 'order') Call get_word(record,word)

        l_n_e=.false. ! abandon assumptions

        mxompl = Min(Max(0,Nint(word_2_real(word))),4)

        mximpl = (mxompl+3)*(mxompl+2)*(mxompl+1)/6

     Else If (word(1:7) == 'molecul') Then

        Call get_word(record,word)
        If (word(1:4) == 'type') Call get_word(record,word)
        mxtmls=Nint(word_2_real(word))

        Do itmols=1,mxtmls

           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 30
              Call lower_case(record)
              Call get_word(record,word)
           End Do

           Do

              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record,comm)
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
                 Do While (ksite < numsit)

                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
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
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:6) == 'constr') Then

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numcon=Nint(word_2_real(word))
                 mtcons=Max(mtcons,numcon)
                 cons%mxtcon=cons%mxtcon+numcon
                 cons%mxcons=cons%mxcons+nummols*numcon

                 Do icon=1,numcon
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:3) == 'pmf') Then

                 pmf%mxpmf=pmf%mxpmf+nummols

                 Do ipmf=1,2
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
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
                    pmf%mxtpmf(ipmf)=Nint(word_2_real(word))

                    Do jpmf=1,pmf%mxtpmf(ipmf)
                       word(1:1)='#'
                       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                          Call get_line(safe,nfield,record,comm)
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
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                    jrgd=Nint(word_2_real(word))
                    mxlrgd=Max(mxlrgd,jrgd)

                    Do lrgd=1,jrgd
                       If (Mod(lrgd,16) == 0) Then
                          word(1:1)='#'
                          Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                             Call get_line(safe,nfield,record,comm)
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
                 inumteth=Nint(word_2_real(word))
                 mtteth=Max(mtteth,inumteth)
                 tether%mxtteth=tether%mxtteth+inumteth
                 tether%mxteth=tether%mxteth+nummols*inumteth

                 Do iteth=1,inumteth
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do
                 End Do

              Else If (word(1:5) == 'bonds') Then

!                 bond%l_tab=.false. ! initialised in bonds_module.f90

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numbonds=Nint(word_2_real(word))
                 mtbond=Max(mtbond,numbonds)
                 bond%max_types=bond%max_types+numbonds
                 bond%max_bonds=bond%max_bonds+nummols*numbonds

                 Do ibonds=1,numbonds
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) bond%l_tab=.true.
                 End Do

                 If (bond%l_tab) Then
                    If (comm%idnode == 0) Open(Unit=ntable, File='TABBND')

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABBND if it's in .xvg format

                    Call get_word(record,word)
                    bond%rcut=Max(bond%rcut,word_2_real(word))

                    Call get_word(record,word)
                    k=Nint(word_2_real(word))
                    bond%bin_tab=Max(bond%bin_tab,k+4)

                    If (comm%idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'angles') Then

!                 angle%l_tab=.false. ! initialised in angles_module.f90

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numang=Nint(word_2_real(word))
                 mtangl=Max(mtangl,numang)
                 angle%max_types=angle%max_types+numang
                 angle%max_angles=angle%max_angles+nummols*numang

                 Do iang=1,numang
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) angle%l_tab=.true.
                 End Do

                 If (angle%l_tab) Then
                    If (comm%idnode == 0) Open(Unit=ntable, File='TABANG')

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABANG if it's in .xvg format

                    Call get_word(record,word) ! no need for cutoff in angles (max is always 180 degrees)
                    k=Nint(word_2_real(word))
                    angle%bin_tab=Max(angle%bin_tab,k+4)

                    If (comm%idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'dihedr') Then

!                 dihedral%l_tab=.false. ! initialised in dihedrals.f90

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numdih=Nint(word_2_real(word))
                 mtdihd=Max(mtdihd,numdih)
                 dihedral%max_types=dihedral%max_types+numdih
                 dihedral%max_angles=dihedral%max_angles+nummols*numdih

                 Do idih=1,numdih
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) dihedral%l_tab=.true.
                 End Do

                 If (dihedral%l_tab) Then
                    If (comm%idnode == 0) Open(Unit=ntable, File='TABDIH')

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABDIH if it's in .xvg format

                    Call get_word(record,word) ! no need for cutoff in angles (max is always 360 degrees
                    k=Nint(word_2_real(word))  ! from -180 to 180)
                    dihedral%bin_tab=Max(dihedral%bin_tab,k+4)

                    If (comm%idnode == 0) Close(Unit=ntable)
                 End If

              Else If (word(1:6) == 'invers') Then

!                 dihedral%l_tab=.false. ! initialised in dihedrals.f90

                 Call get_word(record,word)
                 If (word(1:5) == 'units') Call get_word(record,word)
                 numinv=Nint(word_2_real(word))
                 mtinv=Max(mtinv,numinv)
                 inversion%max_types=inversion%max_types+numinv
                 inversion%max_angles=inversion%max_angles+nummols*numinv

                 Do iinv=1,numinv
                    word(1:1)='#'
                    Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                       Call get_line(safe,nfield,record,comm)
                       If (.not.safe) Go To 30
                       Call get_word(record,word)
                    End Do

                    If (word(1:3) == 'tab' .or. word(1:4)=='-tab' ) inversion%l_tab=.true.
                 End Do

                 If (inversion%l_tab) Then
                    If (comm%idnode == 0) Open(Unit=ntable, File='TABINV')

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40

                    i = Index(record,'#')      ! replace hash as it may occur in
                    If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

                    Call get_word(record,word) ! no need for cutoff in angles (max is always 180 degrees)
                    k=Nint(word_2_real(word))
                    inversion%bin_tab=Max(inversion%bin_tab,k+4)

                    If (comm%idnode == 0) Close(Unit=ntable)
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
               Call get_line(safe,nfield,record,comm)
               If (.not.safe) Go To 30
               Call get_word(record,word)
            End Do
        End Do

        If (mxrdf > 0) mxrdf=Max(mxrdf,(mxatyp*(mxatyp+1))/2)

     Else If (word(1:3) == 'vdw') Then

!        lt_vdw=.false. ! initialised in vdw

        Call get_word(record,word)
        If (word(1:3) == 'tab') Then
           lt_vdw=.true.
        Else
           mxvdw=Nint(word_2_real(word))
        End If

        Do itpvdw=1,mxvdw
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
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
              If (comm%idnode == 0) Open(Unit=ntable, File='TABLE')

              Call get_line(safe,ntable,record,comm)
              If (.not.safe) Go To 40

              Call get_line(safe,ntable,record,comm)
              If (.not.safe) Go To 40
              Call get_word(record,word)

              Call get_word(record,word)
              rvdw=Max(rvdw,word_2_real(word))

              Call get_word(record,word)
              k=Nint(word_2_real(word))
              mxgvdw=Max(mxgvdw,k)

              If (comm%idnode == 0) Close(Unit=ntable)
           End If
        End If

     Else If (word(1:3) == 'met') Then

!        met%tab=-1 ! initialised in metal_module

        Call get_word(record,word)
        mxmet=Nint(word_2_real(word))

        Do itpmet=1,mxmet
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 30
              Call lower_case(record)
              Call get_word(record,word)
           End Do

           Call get_word(record,word)
           Call get_word(record,word)
           met%tab=0 ! for FST metal potentials
           If      (word(1:3) ==  'eam') Then
              met%tab=1
           Else If (word(1:4) == 'eeam') Then
              met%tab=2
           Else If (word(1:4) == '2bea') Then
              met%tab=3
           Else If (word(1:4) == '2bee') Then
              met%tab=4
           Else If (word(1:4) == 'fnsc') Then
              Call get_word(record,word) ; Call get_word(record,word)
              Call get_word(record,word) ; Call get_word(record,word)
              met%rcut=Max(met%rcut,word_2_real(word))
              Call get_word(record,word) ; Call get_word(record,word)
              met%rcut=Max(met%rcut,word_2_real(word))
           Else If (word(1:4) == 'exfs') Then
              Call get_word(record,word) ; Call get_word(record,word)
              Call get_word(record,word) ; Call get_word(record,word)
              Call get_word(record,word) ; Call get_word(record,word)
              met%rcut=Max(met%rcut,word_2_real(word))
              Call get_word(record,word) ; Call get_word(record,word)
              met%rcut=Max(met%rcut,word_2_real(word))
           End If
        End Do

        If (mxmet > 0) Then
           mxmet=Max(mxmet,(mxatyp*(mxatyp+1))/2)

           If      (met%tab == 0) Then
              mxmed=mxmet
           Else If (met%tab == 1) Then
              mxmed=mxatyp
           Else If (met%tab == 2) Then
              mxmed=mxatyp**2
           Else If (met%tab == 3) Then
              mxmed=mxatyp
              mxmds=mxatyp*(mxatyp+1)/2
           Else If (met%tab == 4) Then
              mxmed=mxatyp**2
              mxmds=mxatyp**2
           End If

           If (met%tab > 0) Then
              If (comm%idnode == 0) Open(Unit=ntable, File='TABEAM')

              Call get_line(safe,ntable,record,comm)
              If (.not.safe) Go To 40
              Call get_line(safe,ntable,record,comm)
              If (.not.safe) Go To 40
              Call get_word(record,word)

              Do i=1,Nint(word_2_real(word))
                 Call get_line(safe,ntable,record,comm)
                 If (.not.safe) Go To 40
                 Call get_word(record,word)
                 Call lower_case(word)
                 j=0 ! assume met%rcut is specified
                 If (word(1:4) == 'embe' .or. & ! 2-band embedding functionals
                     word(1:4) == 'demb' .or. word(1:4) == 'semb') j=1 ! no met%rcut is specified
                 If ((word(1:4) == 'dens' .and. met%tab == 2) .or. & ! EEAM
                     (word(2:4) == 'den' .and. (met%tab == 3 .or. met%tab == 4)) .or. & ! sden & dden for 2B extensions
                     word(1:4) == 'pair') Call get_word(record,word) ! skip over one species
                 Call get_word(record,word)                          ! skip over one species

                 Call get_word(record,word)
                 k=Nint(word_2_real(word))
                 met%maxgrid=Max(met%maxgrid,k+4)

                 Call get_word(record,word)
                 Call get_word(record,word)
                 If (j == 0) met%rcut=Max(met%rcut,word_2_real(word))

                 Do j=1,(k+3)/4
                    Call get_line(safe,ntable,record,comm)
                    If (.not.safe) Go To 40
                 End Do
              End Do

              If (comm%idnode == 0) Close(Unit=ntable)
           End If
        End If

     Else If (word(1:7) == 'tersoff') Then

        Call get_word(record,word)
        mxter=Nint(word_2_real(word))

        Do itpter=1,mxter
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
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
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           rct=word_2_real(word)
           rcter=Max(rcter,rct)

           If (potter == 2) Then
              word(1:1)='#'
              Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                 Call get_line(safe,nfield,record,comm)
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
                    Call get_line(safe,nfield,record,comm)
                    If (.not.safe) Go To 30
                    Call get_word(record,word)
                 End Do
              End Do
           End If

           mxter=Max(mxter,(mxatyp*(mxatyp+1))/2)
        End If

     Else If (word(1:3) == 'tbp') Then

        Call get_word(record,word)
        threebody%mxtbp=Nint(word_2_real(word))

        Do itptbp=1,threebody%mxtbp
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
              If (.not.safe) Go To 30
              Call get_word(record,word)
           End Do

           Do i=1,8
              Call get_word(record,word)
           End Do

           rct=word_2_real(word)
           threebody%cutoff=Max(threebody%cutoff,rct)
        End Do

     Else If (word(1:3) == 'fbp') Then

        Call get_word(record,word)
        mxfbp=Nint(word_2_real(word))

        Do itpfbp=1,mxfbp
           word(1:1)='#'
           Do While (word(1:1) == '#' .or. word(1:1) == ' ')
              Call get_line(safe,nfield,record,comm)
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
        kimim=record(1:Len_Trim(record_raw))
        Call kim_cutoff(mxatyp,chr,kimim,rkim,comm)

     Else If (word(1:6) == 'extern') Then

        lext=.true.

        word(1:1)='#'
        Do While (word(1:1) == '#' .or. word(1:1) == ' ')
           Call get_line(safe,nfield,record,comm)
           If (.not.safe) Go To 30
           Call get_word(record,word)
        End Do

        l_usr = (word(1:4) == 'ushr')

     Else If (word(1:5) == 'close') Then

        Go To 10

     End If

  End Do

10 Continue
  If (comm%idnode == 0) Close(Unit=nfield)

! Define legend arrays lengths.  If length > 0 then
! length=Max(length)+1 for the violation excess element

  If (mxshl >  0) mxfshl=1+1 ! One shell per core
  mxf(1)=mxfshl

  If (cons%mxcons > 0) cons%mxfcon=mxb+1
  mxf(2)=cons%mxfcon

  If (pmf%mxpmf  > 0) pmf%mxfpmf=1+1 ! PMFs are global
  mxf(3)=pmf%mxfpmf

  If (mxrgd  > 0) mxfrgd=1+1 ! One RB per particle
  mxf(4)=mxlrgd

  If (tether%mxteth > 0) tether%mxftet=1+1 ! One tether per particle
  mxf(5)=tether%mxftet

  If (bond%max_bonds > 0) bond%max_legend=(mxb*(mxb+1))+1
  mxf(6)=bond%max_legend

  If (angle%max_angles > 0) angle%max_legend=(mxb+1)**2/2+1
  mxf(7)=angle%max_legend

  If (dihedral%max_angles > 0) dihedral%max_legend=((mxb-1)*mxb*(mxb+1))/2+2*mxb+1
  mxf(8)=dihedral%max_legend

  If (inversion%max_angles  > 0) inversion%max_legend=(mxb*(mxb+1))/4+1
  mxf(9)=inversion%max_legend

  Do i=1,9
     mxt(i)=Min(1,mxf(i))
  End Do
  max_exclude = Min( mxnmst , Max( mxfrgd , Sum(mxf)/Max(1,Sum(mxt)) ) * (Max(1,mxshl)+1) )
  If (max_exclude > 0) max_exclude=max_exclude+1 ! violation excess element

  Return

! FIELD file does not exist

20 Continue
  Call error(122)
  Return

30 Continue
  If (comm%idnode == 0) Close(Unit=nfield)
  Call error(52)
  Return

40 Continue
  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine scan_field

  Subroutine read_mpoles(l_top,sumchg,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for reading in the molecular mulitpole
  ! specifications of the system to be simulated
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( In    ) :: l_top
    Real( Kind = wp ),  Intent( InOut ) :: sumchg
    Type( comms_type ), Intent( InOut ) :: comm

    Logical                :: safe,l_rsh,l_ord=.false.

    Character( Len = 200 ) :: record,record1,record2
    Character( Len = 40  ) :: word
    Character( Len = 8   ) :: atom

    Integer                :: itmols,nrept,i,j,k,l,                 &
                              isite,jsite,ksite,lsite,nsite,sitmpl, &
                              ishls,nshels,kshels,isite2,           &
                              ordmpl,ordmpl_start,ordmpl_next,      &
                              ordmpl_min,ordmpl_max,                &
                              indmpl,indmpl_start,indmpl_final

    Real( Kind = wp )      :: charge,scl,polarity,dumping

    Character( Len = 256 ) :: message,messages(3)

  ! open MPOLES data file

    If (comm%idnode == 0) Then
       Open(Unit=nmpldt, File = 'MPOLES', Status = 'old')
    End If
    Call info('electrostatics multipoles specification',.true.)
    If (.not.l_top) Then
      Call info('detailed specification opted out',.true.)
    End If


    Call get_line(safe,nmpldt,record,comm)
    If (.not.safe) Go To 2000

  ! omit first line

    nsite  = 0
    nshels = 0
    sumchg = 0.0_wp

    ordmpl_min = 4
    ordmpl_max = 0

  ! read and process directives from mpols/field file

    Do

       word(1:1)='#'
       Do While (word(1:1) == '#' .or. word(1:1) == ' ')
          Call get_line(safe,nmpldt,record,comm)
          If (.not.safe) Go To 2000
          Call get_word(record,word) ; Call lower_case(word)
       End Do

  ! specify molecular species

       If (word(1:7) == 'molecul') Then

          Call get_word(record,word) ; Call lower_case(word)
          If (word(1:4) == 'type') Call get_word(record,word)

          If (ntpmls == Nint(word_2_real(word))) Then
            Write(message,'(a,i10)') 'number of molecular types ',ntpmls
            Call info(message,.true.)
          Else
            Write(message,'(2(a,i0),a)') &
              'number of molecular types mismatch between FIELD(',ntpmls, &
              ') and MPOLES(',Nint(word_2_real(word)),')'
            Call warning(message,.true.)
            Call error(623)
          End If

  ! read in molecular characteristics for every molecule

          Do itmols=1,ntpmls

             If (l_top) Then
               Write(message,'(a,i10)') 'molecular species type ',itmols
               Call info(message,.true.)
             End If

  ! name of molecular species

             word(1:1)='#'
             Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                Call get_line(safe,nmpldt,record,comm)
                If (.not.safe) Go To 2000
                Call get_word(record,word)
             End Do
             Call strip_blanks(record)
             record1=word(1:Len_Trim(word)+1)//record ; Call lower_case(record1)
             record2=molnam(itmols) ;                   Call lower_case(record2)

             If (record1 == record2) Then
               If (l_top) Then
                 Write(message,'(2a)') 'name of species: ',Trim(molnam(itmols))
                 Call info(message,.true.)
               End If
             Else
               Call warning('molecular names mismatch between FIELD and MPOLES for type',.true.)
               Call error(623)
             End If

  ! read molecular data

             Do

                word(1:1)='#'
                Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                   Call get_line(safe,nmpldt,record,comm)
                   If (.not.safe) Go To 2000
                   Call get_word(record,word) ; Call lower_case(word)
                End Do

  ! number of molecules of this type

                If (word(1:6) == 'nummol') Then

                   Call get_word(record,word)

                   If (nummols(itmols) == Nint(word_2_real(word))) Then
                     If (l_top) Then
                       Write(message,'(a,i10)') 'number of molecules ',nummols(itmols)
                       Call info(message,.true.)
                     End If
                   Else
                     Write(message,'(2(a,i0),a)') &
                       'number of molecular types mismatch between FIELD(',ntpmls, &
                       ') and MPOLES(',Nint(word_2_real(word)),')'
                     Call warning(message,.true.)
                     Call error(623)
                   End If

  ! read in atomic details

                Else If (word(1:5) == 'atoms') Then

                   Call get_word(record,word)

                   If (numsit(itmols) == Nint(word_2_real(word))) Then
                     If (l_top) Then
                       Write(messages(1),'(a,i10)') 'number of atoms/sites ',numsit(itmols)
                       Write(messages(2),'(a)') 'atomic characteristics:'
                       Write(messages(3),'(8x,a4,4x,a4,2x,a16,2x,a6)') &
                         'site','name','multipolar order','repeat'
                       Call info(messages,3,.true.)
                     End If
                   Else
                     Write(message,'(2(a,i0),a)') &
                       'number of molecular types mismatch between FIELD(',ntpmls, &
                       ') and MPOLES(',Nint(word_2_real(word)),')'
                     Call warning(message,.true.)
                     Call error(623)
                   End If

  ! for every molecule of this type get site and atom description

                   ksite=0 ! reference point
                   Do isite=1,numsit(itmols)
                      If (ksite < numsit(itmols)) Then

  ! read atom name, highest pole order supplied, repeat

                         word(1:1)='#'
                         Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                            Call get_line(safe,nmpldt,record,comm)
                            If (.not.safe) Go To 2000
                            Call get_word(record,word)
                         End Do

                         atom=word(1:8)

  ! read supplied pole order

                         Call get_word(record,word)
                         ordmpl=Abs(Nint(word_2_real(word)))
                         indmpl=(ordmpl+3)*(ordmpl+2)*(ordmpl+1)/6

  ! read supplied repetition

                         Call get_word(record,word)
                         nrept=Abs(Nint(word_2_real(word)))
                         If (nrept == 0) nrept=1

                         jsite=nsite+1
                         lsite=jsite+nrept-1

                         Do i=jsite,lsite
                           If (sitnam(i) /= atom) Then ! detect mish-mash
                             Write(message,'(a,i0)') &
                               'site names mismatch between FIELD and MPOLES for site ',ksite+1+i-jsite
                             Call warning(message,.true.)
                             Call error(623)
                           End If
                         End Do

  ! read supplied site polarisation and dumping factor

                         Call get_word(record,word) ; polarity=Abs(word_2_real(word,0.0_wp))
                         Call get_word(record,word) ; dumping =Abs(word_2_real(word,0.0_wp))

                         l_rsh=.true. ! regular or no shelling (Drude)
                         kshels=nshels
                         Do ishls=1,numshl(itmols) ! detect beyond charge shelling
                            kshels=kshels+1

                            isite2=nsite+lstshl(2,kshels)
                            If ((isite2 >= jsite .and. isite2 <= lsite)) Then
                              l_rsh=.false.
                              If (comm%idnode == 0) Then
                                If (ordmpl > 0) Then
                                  Call warning( &
                                    'a shell (of a polarisable multipolar ion)' &
                                    //'can only bear a charge to emulate a' &
                                    //'self-iduced dipole',.true.)
                                End If
                                If (polarity > zero_plus) Then
                                  Call warning( &
                                    'a shell (of a polarisable multipolar ion)' &
                                    //'cannot have its own associated polarisability', &
                                    .true.)
                                End If
                                If (dumping  > zero_plus) Then
                                  Call warning( &
                                    'a shell (of a polarisable multipolar ion)' &
                                    //'cannot have its own associated dumping factor', &
                                    .true.)
                                End IF
                              End If
                            End If
                         End Do

  ! get the min and max order defined for cores/nucleus, ignore irregular shells

                         If (l_rsh) Then
                            ordmpl_min=Min(ordmpl_min,ordmpl)
                            ordmpl_max=Max(ordmpl_max,ordmpl)
                         End If

                         If (l_top) Then
                           If (l_rsh) Then
                             Write(message,'(2x,i10,4x,a8,4x,i2,5x,i10,2f7.3)') &
                               ksite+1,atom,ordmpl,nrept,polarity,dumping
                           Else
                             Write(message,'(2x,i10,4x,a8,4x,i2,5x,i10,2a)') &
                               ksite+1,atom,ordmpl,nrept,polarity,dumping, &
                               Merge(' *ignored* ','           ',polarity>zero_plus), &
                               Merge(' *ignored* ','           ',dumping >zero_plus)
                           End If
                         End If

  ! monopole=charge

                         word(1:1)='#'
                         Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                            Call get_line(safe,nmpldt,record,comm)
                            If (.not.safe) Go To 2000
                            Call get_word(record,word)
                         End Do

                         sitmpl = 1

                         charge=word_2_real(word)

                         chgsit(jsite:lsite)=charge
                         mpllfr(sitmpl,jsite:lsite)=charge
                         If (l_rsh) Then
                            plrsit(jsite:lsite)=polarity
                            dmpsit(jsite:lsite)=dumping
  !                      Else ! initilised to zero in mpoles_module
                         End If

  ! sum absolute charges

                         sumchg=sumchg+Abs(charge)

  ! report

                         If (l_top) Then
                           Write(message,'(2x,a,f10.5)') 'charge ',charge
                           Call info(message,.true.)
                         End If

  ! higher poles counters

                         ordmpl_start = 0
                         ordmpl_next  = ordmpl_start+1
                         indmpl_start = sitmpl+1
                         indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                         Do While (ordmpl_next <= ordmpl)

  ! read line per pole order

                            word(1:1)='#'
                            Do While (word(1:1) == '#' .or. word(1:1) == ' ')
                               Call get_line(safe,nmpldt,record,comm)
                               If (.not.safe) Go To 2000
                               Call get_word(record,word)
                            End Do

  ! Only assign what FIELD says is needed or it is a shell

                            If (ordmpl_next <= Merge(mxompl,1,l_rsh)) Then

                              Do i=indmpl_start,indmpl_final
                                sitmpl = sitmpl+1
                                mpllfr(sitmpl,jsite:lsite)=word_2_real(word)
                                Call get_word(record,word)
                              End Do

  ! report

                              If (l_top) Then
                                If      (ordmpl_next == 1) Then
                                  Write(message,'(2x,a12,1x,3f10.5)') 'dipole', &
                                    mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 2) Then
                                  Write(message,'(2x,a12,1x,6f10.5)') 'quadrupole', &
                                    mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 3) Then
                                  Write(message,'(2x,a12,1x,10f10.5)') 'octupole', &
                                    mpllfr(indmpl_start:indmpl_final,jsite)
                                Else If (ordmpl_next == 4) Then
                                  Write(message,'(2x,a12,1x,15f10.5)') 'hexadecapole', &
                                    mpllfr(indmpl_start:indmpl_final,jsite)
                                End If
                                Call info(message,.true.)
                              End If

  ! rescale poles values by their degeneracy

                               If (ordmpl_next > 1) Then
                                  sitmpl=sitmpl-(indmpl_final-indmpl_start+1) ! rewind
                                  Do i=ordmpl_next,0,-1
                                     l=ordmpl_next-i
                                     Do j=l,0,-1
                                        k=l-j

                                        scl=Exp(factorial(ordmpl_next)-factorial(k)-factorial(j)-factorial(i))
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

                               If (l_top) Then
                                 If (l_rsh) Then
                                   If      (ordmpl_next == 1) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'dipole',' supplied but not required'
                                   Else If (ordmpl_next == 2) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'quadrupole',' supplied but not required'
                                   Else If (ordmpl_next == 3) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'octupole',' supplied but not required'
                                   Else If (ordmpl_next == 4) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'hexadecapole',' supplied but not required'
                                   Else
                                     Write(message,'(2x,a12,i0,a)') &
                                       'pole order ',ordmpl_next,' supplied but not required'
                                   End If
                                 Else
                                   If      (ordmpl_next == 1) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'dipole',' supplied but ignored as invalid'
                                   Else If (ordmpl_next == 2) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'quadrupole',' supplied but ignored as invalid'
                                   Else If (ordmpl_next == 3) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'octupole',' supplied but ignored as invalid'
                                   Else If (ordmpl_next == 4) Then
                                     Write(message,'(2x,a12,1x,a)') &
                                       'hexadecapole',' supplied but ignored as invalid'
                                   Else
                                     Write(message,'(2x,a11,i0,a)') &
                                       'pole order ',ordmpl_next,' supplied but ignored as invalid'
                                   End If
                                 End If
                                 Call info(message,.true.)
                               End If

                            End If

  ! update poles counters

                            ordmpl_next  = ordmpl_next+1
                            indmpl_start = sitmpl+1
                            indmpl_final = (ordmpl_next+3)*(ordmpl_next+2)*(ordmpl_next+1)/6

                         End Do

                         nsite=nsite+nrept
                         ksite=ksite+nrept

                         If (ksite == numsit(itmols)) nshels=kshels

                      End If
                   End Do

  ! finish of data for one molecular type

                Else If (word(1:6) == 'finish') Then

                  Write(message,'(3(a,i0))') &
                    'multipolar electrostatics requested up to order ', &
                    mxompl, ' with specified interactions up order ',  &
                    ordmpl_max,' and least order ', ordmpl_min
                  Call warning(message,.true.)

                  If (ordmpl_max*mxompl == 0) Then
                    Call warning( &
                      'multipolar electrostatics machinery to be used for monompoles ' &
                      //'only electrostatic interactions (point charges only)', &
                      .true.)
                  End If
                  If (ordmpl_max > 4) Then
                    Call warning( &
                      'electrostatic interactions beyond hexadecapole order can ' &
                      //'not be considered and are thus ignored',.true.)
                  End If

                  Go To 1000

                Else

  ! error exit for unidentified directive in molecular data

                   Call strip_blanks(record)
                   Write(message,'(2a)') word(1:Len_Trim(word)+1),record
                   Call info(message,.true.)
                   Call error(12)

                End If

             End Do

  ! just finished with this type molecule data

  1000       Continue

          End Do

  ! close MPOLES data file

       Else If (word(1:5) == 'close') Then

          If (comm%idnode == 0) Close(Unit=nmpldt)

  ! EXIT IF ALL IS OK

          Return

       Else

  ! error exit for unidentified directive

          Call info(word(1:Len_Trim(word)),.true.)
          Call error(4)

       End If

    End Do

    Return

  2000 Continue

    If (comm%idnode == 0) Close(Unit=nmpldt)
    Call error(52)

  End Subroutine read_mpoles
End Module ffield
