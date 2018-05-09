Module bonds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global bond interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
! contrib   - a.v.brukhno april 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use setup,           Only : mxtmls,mxtbnd,mxbond,mxfbnd,mxpbnd,mxgbnd1,mxatdm, &
                           mxgbnd, fourpi,boltz,delr_max,nrite,npdfdt,npdgdt, &
                           engunit,zero_plus,r4pie0,mximpl,ntable,delr_max,nrite, &
                           mxtbnd,mxgbnd,zero_plus,engunit
  Use comms,           Only : comms_type,gsum, gsync, gcheck, gbcast
  Use configuration,   Only : imcon,cell,natms,nlast,lsi,lsa,lfrzn, &
                            chge,xxx,yyy,zzz,fxx,fyy,fzz, cfgname
  Use site,            Only : ntpatm,unqatm
  Use parse,           Only : get_line,get_word,word_2_real
  Use errors_warnings, Only : error, warning, info
  Use numerics,        Only : images, local_index 
  Use coul_mpole,     Only : intra_mcoul 
  Use coul_spole,      Only : intra_coul 

  Implicit None

  Logical,                        Save :: lt_bnd = .false. ! no tabulated potentials opted

  Integer,                        Save :: ntbond  = 0 , &
                                          ntbond1 = 0 , &
                                          ncfbnd  = 0

  Real( Kind = wp ),              Save :: rcbnd = 0.0_wp


  Integer,           Allocatable, Save :: numbonds(:),keybnd(:)
  Integer,           Allocatable, Save :: lstbnd(:,:),listbnd(:,:),legbnd(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmbnd(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpbnd(:)
  Real( Kind = wp ), Allocatable, Save :: vbnd(:,:),gbnd(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfbnd(:),typbnd(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstbnd(:,:)

  Public :: allocate_bonds_arrays , deallocate_bonds_arrays , &
            allocate_bond_pot_arrays , allocate_bond_dst_arrays, &
            bonds_compute, bonds_forces, bonds_table_read

Contains

  Subroutine allocate_bonds_arrays()


    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numbonds(1:mxtmls),        Stat = fail(1))
    Allocate (keybnd(1:mxtbnd),          Stat = fail(2))
    Allocate (lstbnd(1:2,1:mxtbnd),      Stat = fail(3))
    Allocate (listbnd(0:2,1:mxbond),     Stat = fail(4))
    Allocate (legbnd(0:mxfbnd,1:mxatdm), Stat = fail(5))
    Allocate (prmbnd(1:mxpbnd,1:mxtbnd), Stat = fail(6))
    If (lt_bnd) &
    Allocate (ltpbnd(0:mxtbnd),          Stat = fail(7))
    If (mxgbnd1 > 0) &
    Allocate (ldfbnd(0:mxtbnd),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1014)

    numbonds = 0
    keybnd   = 0
    lstbnd   = 0
    listbnd  = 0
    legbnd   = 0

    prmbnd   = 0.0_wp

    If (lt_bnd) &
    ltpbnd   = 0

    If (mxgbnd1 > 0) &
    ldfbnd   = 0

  End Subroutine allocate_bonds_arrays

  Subroutine deallocate_bonds_arrays()

    Integer :: fail

    fail = 0

    Deallocate (numbonds,lstbnd, Stat = fail)

    If (fail > 0) Call error(1029)

  End Subroutine deallocate_bonds_arrays

  Subroutine allocate_bond_pot_arrays()


    Integer :: fail(1:2)

    fail = 0

    Allocate (vbnd(-1:mxgbnd,1:ltpbnd(0)), Stat = fail(1))
    Allocate (gbnd(-1:mxgbnd,1:ltpbnd(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1072)

    vbnd = 0.0_wp
    gbnd = 0.0_wp

  End Subroutine allocate_bond_pot_arrays

  Subroutine allocate_bond_dst_arrays()

    Integer :: fail

    fail = 0

    Allocate (typbnd(-1:2,1:ldfbnd(0)),dstbnd(1:mxgbnd1,1:ldfbnd(0)), Stat = fail)

    If (fail > 0) Call error(1073)

    typbnd = 0
    dstbnd = 0.0_wp

  End Subroutine allocate_bond_dst_arrays


  Subroutine bonds_compute(temp,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bonds distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ), Intent( In    ) :: temp
  Type( comms_type), Intent( InOut ) :: comm

  Logical           :: zero
  Integer           :: fail,ngrid,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,delr,rdlr,dgrid,factor,pdfzero,   &
                       factor1,rrr,dvol,pdfbnd,sum,pdfbnd1,sum1, &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstdbnd(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)
  Character ( Len = 256 )  :: message
  Character ( Len = 256 )  :: messages(2)

  fail = 0
  Allocate (dstdbnd(0:mxgbnd1,1:ldfbnd(0)),pmf(0:mxgbnd1+2),vir(0:mxgbnd1+2), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'bonds_compute - allocation failure'
     Call error(0,message)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! grid interval for pdf tables

  delr = rcbnd/Real(mxgbnd1,wp)
  rdlr = 1.0_wp/delr

! resampling grid and grid interval for pmf tables

  ngrid = Max(Nint(rcbnd/delr_max),mxgbnd1,mxgbnd-4)
  dgrid = rcbnd/Real(ngrid,wp)

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        kk=kk+1
        ll=ll+typbnd(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(ncfbnd,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-5_wp

  Call info('',.true.)
  Call info('BONDS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)',.true.)
  Write(message,'(a,i10,1x,f8.3,3(1x,i10))') '# bins, cutoff, frames, types: ',mxgbnd1,rcbnd,ncfbnd,kk,ll
  Call info(message,.true.)

! open RDF file and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdfdt, File='BNDDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a)') '# BONDS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dr_bin'
     Write(npdfdt,'(a,i10,1x,f8.3,2(1x,i10))') '# bins, cutoff, frames, types: ',mxgbnd1,rcbnd,ncfbnd,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# r(Angstroms)  PDF_norm(r)  PDF_norm(r)/dVol(r)   @   dr_bin = ',delr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        j=j+1

        Write(messages(1),'(a,2(a8,1x),2(i10,1x))') 'type, index, instances: ', &
          unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)
        Write(messages(2),'(a,f8.5)') 'r(Angstroms)  P_bond(r)  Sum_P_bond(r)   @   dr_bin = ',delr
        Call info(messages,2,.true.)
        If (comm%idnode == 0) Then
           Write(npdfdt,'(/,a,2(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i)
        End If

! global sum of data on all nodes

        Call gsum(comm,dstbnd(1:mxgbnd1,i))

! factor in instances (first, pdfbnd is normalised to unity)

        factor1=factor/Real(typbnd(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxgbnd1
           If (zero .and. ig < (mxgbnd1-3)) zero=(dstbnd(ig+2,i) <= 0.0_wp)

           pdfbnd= dstbnd(ig,i)*factor1
           sum = sum + pdfbnd

! null it if < pdfzero

           If (pdfbnd < pdfzero) Then
              pdfbnd1 = 0.0_wp
           Else
              pdfbnd1 = pdfbnd
           End If

           If (sum < pdfzero) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           rrr = (Real(ig,wp)-0.5_wp)*delr
           dvol= fourpi*delr*(rrr*rrr+delr*delr/12.0_wp)

! now pdfbnd is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfbnd= pdfbnd*rdlr

! print out information

           If (.not.zero) Then
              Write(message,'(f11.5,1p,2e14.6)') rrr,pdfbnd1,sum1
              Call info(message,.true.)
           End If
           If (comm%idnode == 0) Then
              Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,pdfbnd,pdfbnd*delr/dvol
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

           dstdbnd(ig,i) = pdfbnd1/dvol ! PDFs density
        End Do
     Else
        dstdbnd(:,i) = 0.0_wp ! PDFs density
     End If
  End Do

  If (comm%idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdgdt, File='BNDPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',delr*Real(mxgbnd1,wp),mxgbnd1,delr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo

     Open(Unit=npdfdt, File='BNDTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',dgrid*Real(ngrid,wp),ngrid,dgrid,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfbnd(0)
     If (typbnd(0,i) > 0) Then
        j=j+1

        If (comm%idnode == 0)  Then
           Write(npdgdt,'(/,a,2(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i), &
                ' (type, index, instances)'
           Write(npdfdt,'(/,a,2(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typbnd(1,i)),unqatm(typbnd(2,i)),j,typbnd(0,i), &
                ' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxgbnd1
           tmp = Real(ig,wp)-0.5_wp
           rrr = tmp*delr

           If (dstdbnd(ig,i) > zero_plus) Then
              fed = -Log(dstdbnd(ig,i))-fed0
              If (fed0 <= zero_plus) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < mxgbnd1-1) Then
                 If (dstdbnd(ig+1,i) <= zero_plus .and. dstdbnd(ig+2,i) > zero_plus) &
                    dstdbnd(ig+1,i) = 0.5_wp*(dstdbnd(ig,i)+dstdbnd(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstdbnd(ig,i) > zero_plus .and. dstdbnd(ig+1,i) > zero_plus) Then
                 dfed = Log(dstdbnd(ig+1,i)/dstdbnd(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == mxgbnd1) Then
              If      (dstdbnd(ig,i) > zero_plus .and. dstdbnd(ig-1,i) > zero_plus) Then
                 dfed = Log(dstdbnd(ig,i)/dstdbnd(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstdbnd(ig-1,i) > zero_plus) Then
              If (dstdbnd(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstdbnd(ig+1,i)/dstdbnd(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstdbnd(ig-1,i))
              End If
           Else If (dstdbnd(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstdbnd(ig+1,i))
           Else If (dfed > 0.0_wp) Then
              dfed = dfed0
           Else
              dfed =-dfed0
           End If

           pmf(ig) = fed
           vir(ig) = dfed

! Print

           If (comm%idnode == 0) &
              Write(npdgdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*tmp
        End Do

! Define edges

        pmf(0)         = 2.0_wp*pmf(1)        -pmf(2)
        vir(0)         = 2.0_wp*vir(1)        -vir(2)
        pmf(mxgbnd1+1) = 2.0_wp*pmf(mxgbnd1)  -pmf(mxgbnd1-1)
        vir(mxgbnd1+1) = 2.0_wp*vir(mxgbnd1)  -vir(mxgbnd1-1)
        pmf(mxgbnd1+2) = 2.0_wp*pmf(mxgbnd1+1)-pmf(mxgbnd1)
        vir(mxgbnd1+2) = 2.0_wp*vir(mxgbnd1+1)-vir(mxgbnd1)

! resample using 3pt interpolation

        Do ig=1,ngrid
           rrr = Real(ig,wp)*dgrid
           ll = Int(rrr*rdlr)

! +0.5_wp due to half-a-bin shift in the original data

           coef = rrr*rdlr-Real(ll,wp)+0.5_wp

           fed0 = pmf(ll)
           fed1 = pmf(ll+1)
           fed2 = pmf(ll+2)

           t1 = fed0 + (fed1 - fed0)*coef
           t2 = fed1 + (fed2 - fed1)*(coef - 1.0_wp)

           fed = t1 + (t2-t1)*coef*0.5_wp

           dfed0 = vir(ll)
           dfed1 = vir(ll+1)
           dfed2 = vir(ll+2)

           t1 = dfed0 + (dfed1 - dfed0)*coef
           t2 = dfed1 + (dfed2 - dfed1)*(coef - 1.0_wp)

           dfed = t1 + (t2-t1)*coef*0.5_wp

           If (comm%idnode == 0) &
              Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*rrr*rdlr
        End Do
     End If
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstdbnd,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'bonds_compute - deallocation failure'
     Call error(0,message)
  End If

End Subroutine bonds_compute

Subroutine bonds_forces(isw,engbnd,virbnd,stress,rcut,keyfce,alpha,epsq,engcpe,vircpe,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating chemical bond energy and force
! terms
!
! isw = 0 - collect statistics
! isw = 1 - calculate forces
! isw = 2 - do both
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov september 2016
! contrib   - a.v.brukhno & i.t.todorov april 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,                             Intent( In    ) :: isw
  Real( Kind = wp ),                   Intent(   Out ) :: engbnd,virbnd
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
  Real( Kind = wp ),                   Intent( In    ) :: rcut,alpha,epsq
  Integer,                             Intent( In    ) :: keyfce
  Real( Kind = wp ),                   Intent( InOut ) :: engcpe,vircpe
  Type( comms_type),                   Intent( InOut ) :: comm

  Logical           :: safe(1:3)
  Integer           :: fail(1:2),i,j,l,ia,ib,keyb,kk
  Real( Kind = wp ) :: rab,rab2,fx,fy,fz,gamma,omega,    &
                       term,term1,term2,eps,sig,         &
                       k,k2,k3,k4,r0,dr,dra,dr2,         &
                       e0,rc,a,b,c,rho,delta,chgprd,     &
                       rdelr,rdr,ppp,vk,vk1,vk2,t1,t2,   &
                       viracc,engc12,virc12,buffer(1:4), &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Character ( Len = 256 )        :: message

  fail=0
  Allocate (lunsafe(1:mxbond),lstopt(0:2,1:mxbond),       Stat=fail(1))
  Allocate (xdab(1:mxbond),ydab(1:mxbond),zdab(1:mxbond), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'bond_forces allocation failure'
     Call error(0,message)
  End If


! calculate atom separation vectors

  Do i=1,ntbond
     lunsafe(i)=.false.

! indices of bonded atoms

     ia=local_index(listbnd(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listbnd(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0) Then ! Tag
        If (lfrzn(ia)*lfrzn(ib) == 0) Then
           If (ia <= natms .or. ib <= natms) Then
              lstopt(0,i)=1
           End If
        End If
     Else                          ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
              (ib > 0 .and. ib <= natms)) .and. &
             (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
     End If

! components of bond vector

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
     End If
  End Do

! Check for uncompressed units

  safe(1) = .not. Any(lunsafe(1:ntbond))
  Call gcheck(comm,safe(1))
  If (.not.safe(1)) Then
     Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
           Do i=1,ntbond
              If (lunsafe(i)) Then
                Write(message,'(a,2(i10,a))') 'global unit number', listbnd(0,1), &
                 ' , with a head particle number', listbnd(1,i), &
                 ' contributes towards next error'
                Call warning(message)
              End If
           End Do
        End If
        Call gsync(comm)
     End Do
     Call error(128)
  End If

! periodic boundary condition

  Call images(imcon,cell,ntbond,xdab,ydab,zdab)

! Initialise safety flag

  safe=.true.

  If (Mod(isw,3) > 0) Then

! zero bond energy and virial accumulators

     engbnd=0.0_wp
     virbnd=0.0_wp

! zero scaled 1-2 electrostatic potential accumulators

     engc12=0.0_wp
     virc12=0.0_wp

     viracc=0.0_wp

! initialise stress tensor accumulators

     strs1=0.0_wp
     strs2=0.0_wp
     strs3=0.0_wp
     strs5=0.0_wp
     strs6=0.0_wp
     strs9=0.0_wp

  End If

! Recover bin size and increment counter

  If (Mod(isw,2) == 0) Then
     rdelr  = Real(mxgbnd1,wp)/rcbnd
     ncfbnd = ncfbnd + 1
  End If

! loop over all specified chemical bond potentials

  Do i=1,ntbond
     If (lstopt(0,i) > 0) Then

! indices of bonded atoms

        ia=lstopt(1,i)
        ib=lstopt(2,i)

! define components of bond vector

        rab2 = xdab(i)**2+ydab(i)**2+zdab(i)**2
        rab  = Sqrt(rab2)

! index of potential function parameters

        kk=listbnd(0,i)
        keyb = Abs(keybnd(kk))

! accumulate the histogram (distribution)

        If (Mod(isw,2) == 0 .and. ia <= natms) Then
           j = ldfbnd(kk)
           l = Min(1+Int(rab*rdelr),mxgbnd1)

           dstbnd(l,j) = dstbnd(l,j) + 1.0_wp

           If (rab > rcbnd) safe(3)=.false. ! catch bondbreaking
        End If
        If (isw == 0) Cycle

! calculate scalar constant terms

        If      (keyb == 0) Then

! null interaction

           omega=0.0_wp
           gamma=0.0_wp

        Else If (keyb == 1) Then

! harmonic potential

           k =prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0

           term=k*dr

           omega=term*0.5_wp*dr
           gamma=-term/rab

        Else If (keyb == 2) Then

! Morse potential

           e0=prmbnd(1,kk)
           r0=prmbnd(2,kk)
           k =prmbnd(3,kk)

           term=Exp(-k*(rab-r0))

           omega=e0*term*(term-2.0_wp)
           gamma=-2.0_wp*e0*k*term*(1.0_wp-term)/rab

        Else If (keyb == 3) Then

! 12-6 potential

           a=prmbnd(1,kk)
           b=prmbnd(2,kk)

           term=rab**(-6)

           omega=term*(a*term-b)
           gamma=6.0_wp*term*(2.0_wp*a*term-b)/rab2

        Else If (keyb == 4) Then

! Lennard-Jones potential

           eps=prmbnd(1,kk)
           sig=prmbnd(2,kk)

           term=(sig/rab)**6

           omega=4.0_wp*eps*term*(term-1.0_wp)
           gamma=24.0_wp*eps*term*(2.0_wp*term-1.0_wp)/rab2

        Else If (keyb == 5) Then

! restrained harmonic

           k =prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0
           dra=Abs(dr)
           rc=prmbnd(3,kk)

           omega=k*(0.5_wp*Min(dra,rc)**2 + rc*Max(dra-rc,0.0_wp))
           gamma=-k*Sign(Min(dra,rc),dr)/rab

        Else If (keyb == 6) Then

! quartic potential

           k2=prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0
           k3=prmbnd(3,kk)
           k4=prmbnd(4,kk)

           dr2=dr**2

           omega=dr2 * (0.5_wp*k2+k3*dr/3.0_wp+0.25_wp*k4*dr2)
           gamma=-dr*(k2+k3*dr+k4*dr2)/rab

        Else If (keyb == 7) Then

! Buckingham exp-6 potential

           a  =prmbnd(1,kk)
           rho=prmbnd(2,kk)
           c  =prmbnd(3,kk)

           term1=a*Exp(-rab/rho)
           term2=-c/rab**6

           omega=term1+term2
           gamma=(term1/rho+6.0_wp*term2/rab)/rab

        Else If (keyb == 8) Then

           omega=0.0_wp
           gamma=0.0_wp

! scaled charge product times dielectric constants

           chgprd=prmbnd(1,kk)*chge(ia)*chge(ib)*r4pie0/epsq
           If ((Abs(chgprd) > zero_plus .or. mximpl > 0) .and. keyfce > 0) Then
              If (mximpl > 0) Then
                 Call intra_mcoul(keyfce,rcut,alpha,epsq,ia,ib,chgprd, &
                      rab,xdab(i),ydab(i),zdab(i),omega,viracc,fx,fy,fz,safe(1))
              Else
                 Call intra_coul(keyfce,rcut,alpha,epsq,chgprd,rab,rab2,omega,gamma,safe(1))

                 fx = gamma*xdab(i)
                 fy = gamma*ydab(i)
                 fz = gamma*zdab(i)

                 viracc = -gamma*rab2
              End If

! correct electrostatic energy and virial

              If (ia <= natms) Then
                 engc12 = engc12 + omega
                 virc12 = virc12 + viracc
              End If
           End If

        Else If (keyb == 9) Then

! extended FENE (Finite Extensive Non-linear Elastic) potential

           k    =prmbnd(1,kk)
           r0   =prmbnd(2,kk)
           delta=prmbnd(3,kk)
           dr   =rab-delta

           term =1.0_wp-(dr/r0)**2
           If (term > 0.0_wp) Then
              omega=-0.5_wp*k*r0**2 * Log(term)
              gamma=-k*dr/term/rab
           Else
              safe(2)=.false.
              omega=0.0_wp
              gamma=0.0_wp
           End If

        Else If (keyb == 10) Then

! MM3-bond-stretch potential

           k =prmbnd(1,kk)
           r0=prmbnd(2,kk)
           dr=rab-r0

           e0=2.55_wp*dr

           omega=k*dr**2*(1.0_wp-e0*(1.0_wp-(7.0_wp/12.0_wp)*e0))
           gamma=2.0_wp*k*dr*(1.0_wp-e0*dr*(1.5_wp+2.0_wp*(7.0_wp/12.0_wp)*e0))/rab

        Else If (keyb == 20) Then

! TABBND potential

           j = ltpbnd(kk)
           If (rab <= vbnd(-1,j)) Then ! rab <= cutpot
              rdr = gbnd(-1,j) ! 1.0_wp/delpot

              l   = Int(rab*rdr)
              ppp = rab*rdr - Real(l,wp)

              vk  = vbnd(l,j)
              vk1 = vbnd(l+1,j)
              vk2 = vbnd(l+2,j)

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              omega = t1 + (t2-t1)*ppp*0.5_wp

              vk  = gbnd(l,j) ; If (l == 0) vk = vk*rab
              vk1 = gbnd(l+1,j)
              vk2 = gbnd(l+2,j)

              t1 = vk  + (vk1 - vk)*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rab
           Else ! bond breaking
              safe(3)=.false.
              omega=0.0_wp
              gamma=0.0_wp
           End If

        Else

! undefined potential

           safe(1)=.false.
           omega=0.0_wp
           gamma=0.0_wp

        End If

! calculate forces and virial additions

        If (keyb /= 8) Then
           fx = gamma*xdab(i)
           fy = gamma*ydab(i)
           fz = gamma*zdab(i)

           viracc=-gamma*rab2
        End If

! add forces

        If (ia <= natms) Then

           fxx(ia)=fxx(ia)+fx
           fyy(ia)=fyy(ia)+fy
           fzz(ia)=fzz(ia)+fz

! calculate bond energy and virial

           engbnd=engbnd+omega
           virbnd=virbnd+viracc

! calculate stress tensor

           strs1 = strs1 + xdab(i)*fx
           strs2 = strs2 + xdab(i)*fy
           strs3 = strs3 + xdab(i)*fz
           strs5 = strs5 + ydab(i)*fy
           strs6 = strs6 + ydab(i)*fz
           strs9 = strs9 + zdab(i)*fz

        End If

        If (ib <= natms) Then

           fxx(ib)=fxx(ib)-fx
           fyy(ib)=fyy(ib)-fy
           fzz(ib)=fzz(ib)-fz

        End If

     End If
  End Do

  If (Mod(isw,3) > 0) Then

! sum contributions to potential and virial

     If (comm%mxnode > 1) Then
        buffer(1) = engbnd
        buffer(2) = virbnd
        buffer(3) = engc12
        buffer(4) = virc12
        Call gsum(comm,buffer(1:4))
        engbnd = buffer(1)
        virbnd = buffer(2)
        engc12 = buffer(3)
        virc12 = buffer(4)
     End If

     engbnd = engbnd - engc12
     virbnd = virbnd - virc12

     engcpe = engcpe + engc12
     vircpe = vircpe + virc12

! complete stress tensor

     stress(1) = stress(1) + strs1
     stress(2) = stress(2) + strs2
     stress(3) = stress(3) + strs3
     stress(4) = stress(4) + strs2
     stress(5) = stress(5) + strs5
     stress(6) = stress(6) + strs6
     stress(7) = stress(7) + strs3
     stress(8) = stress(8) + strs6
     stress(9) = stress(9) + strs9

  End If

! check for undefined potentials + extras

  Call gcheck(comm,safe)
  If (.not.safe(1)) Call error(444)
  If (.not.safe(2)) Call error(655)
  If (.not.safe(3)) Call error(660)

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'bond_forces deallocation failure, node'
     Call error(0,message)
  End If

End Subroutine bonds_forces

Subroutine bonds_table_read(bond_name,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABBND file (for bond potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Character( Len = 16 ), Intent( In    ) :: bond_name(1:mxtbnd)
  Type( comms_type),     Intent( InOut ) :: comm
  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 16 )  :: idbond
  Character( Len = 8   ) :: atom1,atom2

  Integer                :: fail(1:2),ngrid,rtbnd,itbnd,jtbnd,katom1,katom2,jtpatm,i,l
  Real( Kind = wp )      :: cutpot,delpot,dlrpot,rdr,rrr,rrr0, &
                            ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

  Integer,           Allocatable :: read_type(:)
  Real( Kind = wp ), Allocatable :: bufpot(:),bufvir(:)

  Character (Len = 256 )        :: message
  Character (Len = 256 )        :: messages(4)

  If (comm%idnode == 0) Open(Unit=ntable, File='TABBND')

! skip header record

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

! read mesh resolution

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABBND if it's in .xvg format

  Call get_word(record,word)
  cutpot = word_2_real(word)

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  delpot = cutpot/Real(ngrid,wp)

  dlrpot = rcbnd/Real(mxgbnd-4,wp)

! check grid spacing

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > delr_max .and. (.not.safe)) Then
    Write(messages(1),'(a,1p,e15.7)') 'expected (maximum) radial increment : ', delr_max
    Write(messages(2),'(a,1p,e15.7)') 'TABBND file actual radial increment : ', delpot
    Write(messages(3),'(a,0p,i10)') ' expected (minimum) number of grid points : ', mxgbnd-4
    Write(messages(4),'(a,0p,i10)') ' TABBND file actual number of grid points : ', ngrid
    Call info(messages,4,.true.)

     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     Write(message,'(a,i10)') ' TABBND arrays resized for mxgrid = ', mxgbnd-4
     Call info(message,.true.)
  End If

! compare grids dimensions

  If (ngrid < mxgbnd-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgbnd-4,wp),0.0_wp)
     Call error(48)
  End If

  fail=0
  Allocate (read_type(1:ltpbnd(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'error - bonds_table_read allocation failure'
     Call error(0,message)
  End If
  Call allocate_bond_pot_arrays()

  read_type=0 ! initialise read_type
  Do rtbnd=1,ltpbnd(0)
     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABBND if it's in .xvg format

     Call get_word(record,atom1)
     Call get_word(record,atom2)

     katom1=0
     katom2=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0) Then
        Write(message,'(a)') '****',atom1,'***',atom2,'**** entry in TABBND'
        Call info(message,.true.)
        Call error(81)
     End If

! Construct unique name for the tabulated bond

     If (katom1 <= katom2) Then
        idbond = atom1//atom2
     Else
        idbond = atom2//atom1
     End If

! read potential arrays if potential is defined

     itbnd=0
     Do jtbnd=1,ltpbnd(0)
        If (bond_name(jtbnd) == idbond) Then
           Do itbnd=1,mxtbnd
              If (ltpbnd(itbnd) == jtbnd) Exit
           End Do
           Exit
        End If
     End Do

     If (itbnd == 0) Then ! All(bond_name /= idbond)
        Write(message,'(a)') '****',atom1,'***',atom2,'**** entry in TABBND'
        Call info(message,.true.)
        Call error(80)
     End If
     If (Any(read_type == jtbnd)) Then
        Write(message,'(a)') '****',atom1,'***',atom2,'**** entry in TABBND'
        Call info(message,.true.)
        Call error(172)
     Else
        read_type(jtbnd)=jtbnd
     End If

! read in potential & force arrays

     Do i=0,2
        bufpot(i) = 0.0_wp
        bufvir(i) = 0.0_wp
     End Do

! read in the zero and/or first & second data elements (potential & virial)

     If (comm%idnode == 0) Then
        rrr=0.0_wp
        Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

        If (rrr > zero_plus) Then ! no zero element data => extrapolate to zero
           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              Write(messages(1),'(a,1p,e15.7)') ' TABBND stated  radial increment : ', delpot
              Write(messages(2),'(a,1p,e15.7)') ' TABBND read-in radial increment : ', rrr
              Call info(messages,2,.true.)
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0
           rrr0      = rrr

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              Write(messages(1),'(a,1p,e15.7)') ' TABBND stated  radial increment : ', delpot
              Write(messages(2),'(a,1p,e15.7)') ' TABBND read-in radial increment : ', rrr-rrr0
              Call info(messages,2,.true.)
           End If

           bufpot(2) = bufp0
           bufvir(2) = bufv0

! linear extrapolation for grid point 0 at distances close to 0

           bufpot(0) = 2.0_wp*bufpot(1)-bufpot(2)
           bufvir(0) = (2.0_wp*bufvir(1)-0.5_wp*bufvir(2))/delpot
        Else ! zero element data found => read in the first element for checking delr
           bufpot(0) = bufp0
           bufvir(0) = bufv0 ! virial/distance - finite force!

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              Write(messages(1),'(a,1p,e15.7)') ' TABBND stated  radial increment : ', delpot
              Write(messages(2),'(a,1p,e15.7)') ' TABBND read-in radial increment : ', rrr
              Call info(messages,2,.true.)
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           bufpot(2) = bufp0
           bufvir(2) = bufv0
        End If
     End If

     Do i=3,ngrid
        If (comm%idnode == 0) Then
           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufpot(i),bufvir(i)
        Else
           bufpot(i) = 0.0_wp
           bufvir(i) = 0.0_wp
        End If
     End Do

     Call gbcast(comm,bufpot,0)
     Call gbcast(comm,bufvir,0)

! reconstruct arrays using 3pt interpolation

     If (remake) Then
        Do i=1,mxgbnd-4
           rrr = Real(i,wp)*dlrpot
           l   = Int(rrr*rdr)
           ppp = rrr*rdr-Real(l,wp)

           vk  = bufpot(l)

! linear extrapolation for the grid points just beyond the cutoff

           If (l+2 > ngrid) Then
              If (l+1 > ngrid) Then
                 vk1 = 2.0_wp*bufpot(l)-bufpot(l-1)
                 vk2 = 2.0_wp*vk1-bufpot(l)
              Else
                 vk1 = bufpot(l+1)
                 vk2 = 2.0_wp*bufpot(l+1)-bufpot(l)
              End If
           Else
              vk1 = bufpot(l+1)
              vk2 = bufpot(l+2)
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           vbnd(i,jtbnd) = t1 + (t2-t1)*ppp*0.5_wp
           vbnd(i,jtbnd) = vbnd(i,jtbnd)*engunit ! convert to internal units

           vk  = bufvir(l)

! linear extrapolation for the grid points just beyond the cutoff

           If (l+2 > ngrid) Then
              If (l+1 > ngrid) Then
                 vk1 = 2.0_wp*bufvir(l)-bufvir(l-1)
                 vk2 = 2.0_wp*vk1-bufvir(l)
              Else
                 vk1 = bufpot(l+1)
                 vk2 = 2.0_wp*bufvir(l+1)-bufvir(l)
              End If
           Else
              vk1 = bufvir(l+1)
              vk2 = bufvir(l+2)
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           gbnd(i,jtbnd) = t1 + (t2-t1)*ppp*0.5_wp
           gbnd(i,jtbnd) = gbnd(i,jtbnd)*engunit ! convert to internal units
        End Do

        vbnd(-1,jtbnd) = cutpot
        gbnd(-1,jtbnd) = 1.0_wp/delpot
     Else
        Do i=1,mxgbnd-4
           vbnd(i,jtbnd) = bufpot(i)*engunit ! convert to internal units
           gbnd(i,jtbnd) = bufvir(i)*engunit ! convert to internal units
        End Do

! linear extrapolation for the grid point just beyond the cutoff

        vbnd(mxgbnd-3,jtbnd) = 2.0_wp*vbnd(mxgbnd-4,jtbnd) - vbnd(mxgbnd-5,jtbnd)
        gbnd(mxgbnd-3,jtbnd) = 2.0_wp*gbnd(mxgbnd-4,jtbnd) - gbnd(mxgbnd-5,jtbnd)

        vbnd(-1,jtbnd) = cutpot
        gbnd(-1,jtbnd) = 1.0_wp/delpot
     End If

! grid point at 0 and linear extrapolation for the grid point at mxgbnd-2

     vbnd(0,jtbnd) = bufpot(0)
     gbnd(0,jtbnd) = bufvir(0)

     vbnd(mxgbnd-2,jtbnd) = 2.0_wp*vbnd(mxgbnd-3,jtbnd) - vbnd(mxgbnd-4,jtbnd)
     gbnd(mxgbnd-2,jtbnd) = 2.0_wp*gbnd(mxgbnd-3,jtbnd) - gbnd(mxgbnd-4,jtbnd)
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=ntable)
  End If
  Call info('',.true.)
  Call info('potential tables read from TABBND file')

! Break if not safe

  Call gcheck(comm,safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'error - bonds_table_read deallocation failure'
     Call error(0,message)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine bonds_table_read


End Module bonds
