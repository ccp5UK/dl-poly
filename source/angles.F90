Module angles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global valence angle interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2014
! contrib   - a.v.brukhno march 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,  Only : comms_type,gsum,gbcast,gsync,gcheck
  Use setup,  Only : pi,boltz,delth_max,nrite,npdfdt,npdgdt, &
                            mxgang,mxgang1,engunit,zero_plus, mxangl, twopi, &
                            delth_max,ntable,mxtang,mxatdm,mxfang,mxpang,mxtmls
  Use site,   Only : unqatm,ntpatm
  Use configuration, Only : imcon,cell,natms,nlast,lsi,lsa,lfrzn, &
                            xxx,yyy,zzz,fxx,fyy,fzz,cfgname
  Use parse, Only : get_line,get_word,word_2_real
  Use errors_warnings, Only : error,warning
  Use numerics, Only : local_index,images

  Implicit None

  Logical,                        Save :: lt_ang = .false. ! no tabulated potentials opted

  Integer,                        Save :: ntangl  = 0 , &
                                          ntangl1 = 0 , &
                                          ncfang  = 0


  Integer,           Allocatable, Save :: numang(:),keyang(:)
  Integer,           Allocatable, Save :: lstang(:,:),listang(:,:),legang(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmang(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpang(:)
  Real( Kind = wp ), Allocatable, Save :: vang(:,:),gang(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfang(:),typang(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstang(:,:)

  Public :: allocate_angles_arrays , deallocate_angles_arrays , &
            allocate_angl_pot_arrays , allocate_angl_dst_arrays, angles_compute, &
            angles_forces, angles_table_read

Contains

  Subroutine allocate_angles_arrays()

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numang(1:mxtmls),          Stat = fail(1))
    Allocate (keyang(1:mxtang),          Stat = fail(2))
    Allocate (lstang(1:3,1:mxtang),      Stat = fail(3))
    Allocate (listang(0:3,1:mxangl),     Stat = fail(4))
    Allocate (legang(0:mxfang,1:mxatdm), Stat = fail(5))
    Allocate (prmang(1:mxpang,1:mxtang), Stat = fail(6))
    If (lt_ang) &
    Allocate (ltpang(0:mxtang),          Stat = fail(7))
    If (mxgang1 > 0) &
    Allocate (ldfang(0:mxtang),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1013)

    numang  = 0
    keyang  = 0
    lstang  = 0
    listang = 0
    legang  = 0

    prmang  = 0.0_wp

    If (lt_ang) &
    ltpang  = 0

    If (mxgang1 > 0) &
    ldfang  = 0

  End Subroutine allocate_angles_arrays

  Subroutine deallocate_angles_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numang,lstang, Stat = fail)

    If (fail > 0) Call error(1028)

  End Subroutine deallocate_angles_arrays

  Subroutine allocate_angl_pot_arrays()

    Integer :: fail(1:2)

    fail = 0

    Allocate (vang(-1:mxgang,1:ltpang(0)), Stat = fail(1))
    Allocate (gang(-1:mxgang,1:ltpang(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1074)

    vang = 0.0_wp
    gang = 0.0_wp

  End Subroutine allocate_angl_pot_arrays

  Subroutine allocate_angl_dst_arrays()

    Integer :: fail

    fail = 0

    Allocate (typang(-1:3,1:ldfang(0)),dstang(1:mxgang1,1:ldfang(0)), Stat = fail)

    If (fail > 0) Call error(1075)

    typang = 0
    dstang = 0.0_wp

  End Subroutine allocate_angl_dst_arrays
  
  Subroutine angles_compute(temp,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating angles distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ),  Intent( In    ) :: temp
  Type( comms_type ), Intent( InOut ) :: comm

  Logical           :: zero
  Integer           :: fail,ngrid,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,rad2dgr,dgr2rad,delth,rdlth,dgrid,factor,pdfzero, &
                       factor1,theta,sinth,rsint,pdfang,sum,pdfang1,sum1,        &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstdang(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  Character( Len = 256 ) :: message

  fail = 0
  Allocate (dstdang(0:mxgang1,1:ldfang(0)),pmf(0:mxgang1+2),vir(0:mxgang1+2), Stat = fail)
  If (fail > 0) Then
     Write(message,'(/,1x,a)') 'angles_compute - allocation failure, node'
     Call error(0,message)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! conversion: radians <-> degrees (affects not only angle units but also force units!)

  rad2dgr = 180.0_wp/pi
  dgr2rad = pi/180.0_wp

! grid interval for pdf/pmf tables

  delth = pi/Real(mxgang1,wp)
  rdlth = Real(mxgang1,wp)/180.0_wp

! resampling grid and grid interval for pmf tables

  ngrid = Max(Nint(180.0_wp/delth_max),mxgang1,mxgang-4)
  dgrid = pi/Real(ngrid,wp)

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,ldfang(0)
     If (typang(0,i) > 0) Then
        kk=kk+1
        ll=ll+typang(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(ncfang,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-5_wp

  If (comm%idnode == 0) Then
     Write(nrite,'(/,/,12x,a)') 'ANGLES : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)'
     Write(nrite,'(/,1x,a,5(1x,i10))') &
           '# bins, cutoff, frames, types: ',mxgang1,180,ncfang,kk,ll
  End If

! open RDF file and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdfdt, File='ANGDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a)') '# ANGLES: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
     Write(npdfdt,'(a,4(1x,i10))') '# bins, cutoff, frames, types: ',mxgang1,180,ncfang,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)  PDF_norm(Theta)/Sin(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfang(0)
     If (typang(0,i) > 0) Then
        j=j+1

        If (comm%idnode == 0) Then
           Write(nrite,'(/,1x,a,3(a8,1x),2(i10,1x))') 'type, index, instances: ', &
                unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i)
           Write(nrite,'(/,1x,a,f8.5,/)') &
                'Theta(degrees)  PDF_ang(Theta)  Sum_PDF_ang(Theta)   @   dTheta_bin = ',delth*rad2dgr

           Write(npdfdt,'(/,a,3(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i)
        End If

! global sum of data on all nodes

        Call gsum(comm,dstang(1:mxgang1,i))

! factor in instances (first, pdfang is normalised to unity)

        factor1=factor/Real(typang(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,mxgang1
           If (zero .and. ig < (mxgang1-3)) zero=(dstang(ig+2,i) <= 0.0_wp)

           pdfang = dstang(ig,i)*factor1
           sum = sum + pdfang

! null it if < pdfzero

           If (pdfang < pdfzero) Then
              pdfang1 = 0.0_wp
           Else
              pdfang1 = pdfang
           End If

           If (sum < pdfzero) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           theta = (Real(ig,wp)-0.5_wp)*delth

! Jacobian rsint = 1.0_wp/sinth must only be used for the PMF calculations below,
! which require renormalising: dstdang(ig,i) = pdfang1*rsint

           sinth = Max(1.0e-10_wp,Sin(theta))
           rsint = 1.0_wp/sinth

! now pdfang is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfang = pdfang*rdlth

! print out information

           theta  = theta*rad2dgr
           If (comm%idnode == 0) Then
              If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") theta,pdfang1,sum1
              Write(npdfdt,"(f11.5,1p,2e14.6)") theta,pdfang,pdfang*rsint
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF/sinth)
! would have poorly-defined noisy "borders/walls"

           dstdang(ig,i) = pdfang1*rsint ! PDFs density
        End Do
     Else
        dstdang(:,i) = 0.0_wp ! PDFs density
     End If
  End Do

  If (comm%idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdgdt, File='ANGPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,2i6,f12.5,i10,a,e15.7)') '# ',mxgang1,180,delth*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) = ', kT2engo

     Open(Unit=npdfdt, File='ANGTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,2i6,f12.5,i10,a,e15.7)') '# ',ngrid,180,dgrid*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) = ', kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,ldfang(0)
     If (typang(0,i) > 0) Then
        j=j+1

        If (comm%idnode == 0) Then
           Write(npdgdt,'(/,a,3(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i), &
                ' (type, index, instances)'
           Write(npdfdt,'(/,a,3(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(typang(1,i)),unqatm(typang(2,i)),unqatm(typang(3,i)),j,typang(0,i), &
                ' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,mxgang1
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth

           If (dstdang(ig,i) > zero_plus) Then
              fed = -Log(dstdang(ig,i))-fed0
              If (fed0 <= zero_plus) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < mxgang1-1) Then
                 If (dstdang(ig+1,i) <= zero_plus .and. dstdang(ig+2,i) > zero_plus) &
                    dstdang(ig+1,i) = 0.5_wp*(dstdang(ig,i)+dstdang(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstdang(ig,i) > zero_plus .and. dstdang(ig+1,i) > zero_plus) Then
                 dfed = Log(dstdang(ig+1,i)/dstdang(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == mxgang1) Then
              If      (dstdang(ig,i) > zero_plus .and. dstdang(ig-1,i) > zero_plus) Then
                 dfed = Log(dstdang(ig,i)/dstdang(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstdang(ig-1,i) > zero_plus) Then
              If (dstdang(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstdang(ig+1,i)/dstdang(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstdang(ig-1,i))
              End If
           Else If (dstdang(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstdang(ig+1,i))
           Else If (dfed > 0.0_wp) Then
              dfed = dfed0
           Else
              dfed =-dfed0
           End If

           pmf(ig) = fed
           vir(ig) = dfed

! Print

           If (comm%idnode == 0) &
              Write(npdgdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do

! Define edges

        pmf(0)         = 2.0_wp*pmf(1)        -pmf(2)
        vir(0)         = 2.0_wp*vir(1)        -vir(2)
        pmf(mxgang1+1) = 2.0_wp*pmf(mxgang1)  -pmf(mxgang1-1)
        vir(mxgang1+1) = 2.0_wp*vir(mxgang1)  -vir(mxgang1-1)
        pmf(mxgang1+2) = 2.0_wp*pmf(mxgang1+1)-pmf(mxgang1)
        vir(mxgang1+2) = 2.0_wp*vir(mxgang1+1)-vir(mxgang1)

! resample using 3pt interpolation

        Do ig=1,ngrid
           theta = Real(ig,wp)*dgrid
           ll = Int(theta/delth)

! +0.5_wp due to half-a-bin shift in the original data

           coef = theta/delth-Real(ll,wp)+0.5_wp

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
              Write(npdfdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do
     End If
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstdang,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(message,'(/,1x,a)') 'angles_compute - deallocation failure'
     Call error(0,message)
  End If

End Subroutine angles_compute

Subroutine angles_forces(isw,engang,virang,stress,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bond angle energy and force terms
!
! isw = 0 - collect statistics
! isw = 1 - calculate forces
! isw = 2 - do both
!
! copyright - daresbury laboratory
! author    - w.smith may 1992
! amended   - i.t.todorov march 2016
! contrib   - a.v.brukhno & i.t.todorov april 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,                             Intent( In    ) :: isw
  Real( Kind = wp ),                   Intent(   Out ) :: engang,virang
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
  Type( comms_type),                   Intent( InOut ) :: comm

  Logical           :: safe
  Integer           :: fail(1:3),i,j,l,ia,ib,ic,keya,kk
  Real( Kind = wp ) :: xab,yab,zab,rab,rrab, xbc,ybc,zbc,rbc,rrbc, &
                       theta,cost,sint,rsint,                      &
                       fxa,fxc,fya, fyc,fza,fzc,                   &
                       rdelth,rdr,ppp,vk,vk1,vk2,t1,t2,            &
                       k,k2,k3,k4,theta0,dtheta,dthpi,dth0pi,dth,  &
                       rho,rho1,rho2,switch,a,b,c,delta,m,dr1,dr2, &
                       gr,rm,tmp,pterm,gamma,gamsa,gamsc,vterm,    &
                       strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Real( Kind = wp ), Allocatable :: xdbc(:),ydbc(:),zdbc(:)

  Character( Len = 256 ) :: message

  fail=0
  Allocate (lunsafe(1:mxangl),lstopt(0:3,1:mxangl),       Stat=fail(1))
  Allocate (xdab(1:mxangl),ydab(1:mxangl),zdab(1:mxangl), Stat=fail(2))
  Allocate (xdbc(1:mxangl),ydbc(1:mxangl),zdbc(1:mxangl), Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(/,1x,a)') 'angles_forces allocation failure'
     Call error(0,message)
  End If


! calculate atom separation vectors

  Do i=1,ntangl
     lunsafe(i)=.false.

! indices of angle bonded atoms

     ia=local_index(listang(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listang(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib
     ic=local_index(listang(3,i),nlast,lsi,lsa) ; lstopt(3,i)=ic

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0 .and. ic > 0) Then ! Tag
        If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic) == 0) Then
           If (ia <= natms .or. ib <= natms .or. ic <= natms) Then
              lstopt(0,i)=1
            End If
        End If
     Else                                       ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
              (ib > 0 .and. ib <= natms) .or.   &
              (ic > 0 .and. ic <= natms)) .and. &
             (ia == 0 .or. ib == 0 .or. ic == 0) ) lunsafe(i)=.true.
     End If

! define components of bond vectors

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)

        xdbc(i)=xxx(ic)-xxx(ib)
        ydbc(i)=yyy(ic)-yyy(ib)
        zdbc(i)=zzz(ic)-zzz(ib)
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
!
!        xdbc(i)=0.0_wp
!        ydbc(i)=0.0_wp
!        zdbc(i)=0.0_wp
     End If
  End Do

! Check for uncompressed units

  safe = .not. Any(lunsafe(1:ntangl))
  Call gcheck(comm,safe)
  If (.not.safe) Then
     Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
           Do i=1,ntangl
             If (lunsafe(i)) Then
               Write(message,'(/,1x,a,2(i10,a))')     &
                 'global unit number', listang(0,i), &
                 ' , with a head particle number', listang(1,i),   &
                 ' contributes towards next error'
               Call warning(message)
             End If
           End Do
        End If
        Call gsync(comm)
     End Do
     Call error(130)
  End If

! periodic boundary condition

  Call images(imcon,cell,ntangl,xdab,ydab,zdab)
  Call images(imcon,cell,ntangl,xdbc,ydbc,zdbc)

  If (Mod(isw,3) > 0) Then

! Initialise safety flag

     safe=.true.

! zero angle energy accumulator

     engang=0.0_wp
     virang=0.0_wp

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
     rdelth = Real(mxgang1,wp)/pi
     ncfang = ncfang + 1
  End If

! loop over all specified angle potentials

  Do i=1,ntangl
     If (lstopt(0,i) > 0) Then

! indices of bonded atoms

        ia=lstopt(1,i)
        ib=lstopt(2,i)
        ic=lstopt(3,i)

! define components of first bond vector

        rab = Sqrt(xdab(i)**2+ydab(i)**2+zdab(i)**2)
        rrab = 1.0_wp/rab

        xab=xdab(i)*rrab
        yab=ydab(i)*rrab
        zab=zdab(i)*rrab

! define components of second bond vector

        rbc = Sqrt(xdbc(i)**2+ydbc(i)**2+zdbc(i)**2)
        rrbc = 1.0_wp/rbc

        xbc=xdbc(i)*rrbc
        ybc=ydbc(i)*rrbc
        zbc=zdbc(i)*rrbc

! determine bond angle and calculate potential energy

        cost=(xab*xbc+yab*ybc+zab*zbc)
        If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)
        theta=Acos(cost)
        sint=Max(1.0e-10_wp,Sqrt(1.0_wp-cost**2))
        rsint=1.0_wp/sint

! index of potential function parameters

        kk=listang(0,i)
        keya = Abs(keyang(kk))

! accumulate the histogram (distribution)

        If (Mod(isw,2) == 0 .and. ib <= natms) Then
           j = ldfang(kk)
           l = Min(1+Int(theta*rdelth),mxgang1)

           dstang(l,j) = dstang(l,j) + 1.0_wp
        End If
        If (isw == 0) Cycle

        If      (keya == 1) Then

! harmonic potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0

           tmp   =k*dtheta

           pterm=tmp*0.5_wp*dtheta
           gamma=tmp*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 2) Then

! quartic potential

           k2    =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           k3    =prmang(3,kk)
           k4    =prmang(4,kk)

           pterm=0.5_wp*k2*dtheta**2+(k3/3.0_wp)*dtheta**3+0.25*k4*dtheta**4
           gamma=dtheta*(k2+k3*dtheta+k4*dtheta**2)*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 3) Then

! truncated Harmonic potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho   =prmang(3,kk)
           switch=-(rab**8+rbc**8)/rho**8

           tmp   =k*dtheta*Exp(switch)

           pterm=tmp*0.5_wp*dtheta
           gamma=tmp*rsint
           gamsa=pterm*8.0_wp*rab**7/rho**8
           gamsc=pterm*8.0_wp*rbc**7/rho**8
           vterm=pterm*8.0_wp*switch

        Else If (keya == 4) Then

! screened Harmonic potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho1  =prmang(3,kk)
           rho2  =prmang(4,kk)
           switch=-(rab/rho1+rbc/rho2)

           tmp   =k*dtheta*Exp(switch)

           pterm=tmp*0.5_wp*dtheta
           gamma=tmp*rsint
           gamsa=pterm/rho1
           gamsc=pterm/rho2
           vterm=pterm*switch

        Else If (keya == 5) Then

! screened Vessal potential (type 1)

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dth0pi=theta0-pi
           dthpi =theta -pi
           dth   =dth0pi**2-dthpi**2
           rho1  =prmang(3,kk)
           rho2  =prmang(4,kk)
           switch=-(rab/rho1+rbc/rho2)

           tmp   =(k*dth/(2.0_wp*dth0pi**2)) * Exp(switch)

           pterm=tmp*0.25_wp*dth
           gamma=tmp*dthpi*rsint
           gamsa=pterm/rho1
           gamsc=pterm/rho2
           vterm=pterm*switch

        Else If (keya == 6) Then

! truncated Vessal potential (type 2)

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           dth0pi=theta0-pi
           dthpi =theta -pi
           a     =prmang(3,kk)
           rho   =prmang(4,kk)
           switch=-(rab**8+rbc**8)/rho**8

           tmp   =k*dtheta*Exp(switch)

           pterm=tmp * dtheta * (theta**a * (dthpi+dth0pi)**2 + 0.5_wp*a * pi**(a-1.0_wp) * dth0pi**3)
           gamma=tmp * (theta**(a-1.0_wp) * (dthpi+dth0pi) *                                 &
                 ((a+4.0_wp)*theta**2 - twopi*(a+2.0_wp)*theta - a*theta0*(dth0pi-pi)) + &
                 a*pi**(a-1.0_wp) * dth0pi**3) * rsint
           gamsa=pterm*8.0_wp*rab**7/rho**8
           gamsc=pterm*8.0_wp*rbc**7/rho**8
           vterm=pterm*8.0_wp*switch

        Else If (keya == 7) Then

! harmonic cosine potential (note cancellation of sint in gamma)

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=Cos(theta)-Cos(theta0)

           tmp   =k*dtheta

           pterm=tmp*0.5_wp*dtheta
           gamma=-tmp
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 8) Then

! ordinary cosine potential

           k    =prmang(1,kk)
           delta=prmang(2,kk)
           m    =prmang(3,kk)
           a    =m*theta-delta

           pterm=k*(1.0_wp+Cos(a))
           gamma=-k*m*Sin(a)*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 9) Then

! MM3-stretch-bend potential

           a     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho1  =prmang(3,kk)
           dr1   =rab-rho1
           rho2  =prmang(4,kk)
           dr2   =rbc-rho2

           tmp   =a*dr1*dr2

           pterm=tmp*dtheta
           gamma=tmp*rsint
           gamsa=-pterm/dr1
           gamsc=-pterm/dr2
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 10) Then

! compass stretch-stretch potential

           a     =prmang(1,kk)
           rho1  =prmang(2,kk)
           dr1   =rab-rho1
           rho2  =prmang(3,kk)
           dr2   =rbc-rho2

           pterm=a*dr1*dr2
           gamma=0.0_wp
           gamsa=-a*dr2
           gamsc=-a*dr1
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 11) Then

! compass stretch-bend potential

           a     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           rho1  =prmang(3,kk)
           dr1   =rab-rho1

           tmp   =a*dr1

           pterm=tmp*dtheta
           gamma=tmp*rsint
           gamsa=-a*dtheta
           gamsc=0.0_wp
           vterm=-gamsa*rab

        Else If (keya == 12) Then

! combined compass angle potential with 3 coupling terms

           a     =prmang(1,kk)
           b     =prmang(2,kk)
           c     =prmang(3,kk)
           theta0=prmang(4,kk)
           dtheta=theta-theta0
           rho1  =prmang(5,kk)
           dr1   =rab-rho1
           rho2  =prmang(6,kk)
           dr2   =rbc-rho2

           tmp   =b*dr1+c*dr2

           pterm=a*dr1*dr2 + dtheta*tmp
           gamma=tmp*rsint
           gamsa=-a*dr2-b*dtheta
           gamsc=-a*dr1-c*dtheta
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 13) Then

! MM3-angle-bend potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0

           pterm=k*dtheta**2 * (1.0_wp-1.4e-2_wp*dtheta+5.60e-5_wp*dtheta**2 - &
                                7.0e-7_wp*dtheta**3+2.20e-8_wp*dtheta**4)
           gamma=k*dtheta    * (2.0_wp-4.2e-2_wp*dtheta+2.24e-4_wp*dtheta**2 - &
                                3.5e-6_wp*dtheta**3+1.32e-7_wp*dtheta**4)*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 14) Then

! KKY potential

           k     =prmang(1,kk)
           theta0=prmang(2,kk)
           dtheta=theta-theta0
           gr    =prmang(3,kk)
           rm    =prmang(4,kk)
           dr1   =rab-rm
           rho1  =Exp(-0.5_wp*gr*dr1)
           dr2   =rbc-rm
           rho2  =Exp(-0.5_wp*gr*dr2)
           rho   =rho1*rho2

           pterm=k*(Cos(2.0_wp*dtheta)-1.0_wp)*rho
           gamma=-k*2.0_wp*Sin(2.0_wp*dtheta)/sint
           gamsa=-0.5_wp*gr*pterm
           gamsc=gamsa
           vterm=-gamsa*(rab+rbc)

        Else If (keya == 20) Then

! TABANG potential

           j = ltpang(kk)
           rdr = gang(-1,j) ! 1.0_wp/delpot (in rad^-1)

           l   = Int(theta*rdr)
           ppp = theta*rdr - Real(l,wp)

           vk  = vang(l,j)
           vk1 = vang(l+1,j)
           vk2 = vang(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           pterm = t1 + (t2-t1)*ppp*0.5_wp

           vk  = gang(l,j) ; If (l == 0) vk = vk*theta
           vk1 = gang(l+1,j)
           vk2 = gang(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           gamma =-(t1 + (t2-t1)*ppp*0.5_wp)*rsint

           vterm=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp

        Else

! undefined potential

           safe=.false.
           pterm=0.0_wp
           gamma=0.0_wp
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        End If

! calculate atomic forces

        fxa = gamma*(xbc-xab*cost)*rrab+gamsa*xab
        fya = gamma*(ybc-yab*cost)*rrab+gamsa*yab
        fza = gamma*(zbc-zab*cost)*rrab+gamsa*zab

        fxc = gamma*(xab-xbc*cost)*rrbc+gamsc*xbc
        fyc = gamma*(yab-ybc*cost)*rrbc+gamsc*ybc
        fzc = gamma*(zab-zbc*cost)*rrbc+gamsc*zbc

        If (ia <= natms) Then

           fxx(ia)=fxx(ia)+fxa
           fyy(ia)=fyy(ia)+fya
           fzz(ia)=fzz(ia)+fza

        End If

        If (ib <= natms) Then

! energy and virial (associated to the head atom)

           engang=engang+pterm
           virang=virang+vterm

! calculate stress tensor (associated to the head atom)

           strs1 = strs1 + rab*xab*fxa + rbc*xbc*fxc
           strs2 = strs2 + rab*xab*fya + rbc*xbc*fyc
           strs3 = strs3 + rab*xab*fza + rbc*xbc*fzc
           strs5 = strs5 + rab*yab*fya + rbc*ybc*fyc
           strs6 = strs6 + rab*yab*fza + rbc*ybc*fzc
           strs9 = strs9 + rab*zab*fza + rbc*zbc*fzc

           fxx(ib)=fxx(ib)-fxa-fxc
           fyy(ib)=fyy(ib)-fya-fyc
           fzz(ib)=fzz(ib)-fza-fzc

        End If

        If (ic <= natms) Then

           fxx(ic)=fxx(ic)+fxc
           fyy(ic)=fyy(ic)+fyc
           fzz(ic)=fzz(ic)+fzc

        End If

     End If
  End Do

  If (Mod(isw,3) > 0) Then

! global sum of angular potential and virial

     If (comm%mxnode > 1) Then
        buffer(1)=engang
        buffer(2)=virang
        Call gsum(comm,buffer(1:2))
        engang=buffer(1)
        virang=buffer(2)
     End If

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

! check for undefined potentials

     Call gcheck(comm,safe)
     If (.not.safe) Call error(440)

  End If

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  Deallocate (xdbc,ydbc,zdbc, Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(/,1x,a)') 'angles_forces deallocation failure'
     Call error(0,message)
  End If

End Subroutine angles_forces

Subroutine angles_table_read(angl_name,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABANG file (for angle potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Character( Len = 24 ), Intent( In    ) :: angl_name(1:mxtang)
  Type(comms_type),      Intent( InOut ) :: comm

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 24 )  :: idangl
  Character( Len = 8   ) :: atom1,atom2,atom3

  Integer                :: fail(1:2),ngrid,rtang,itang,jtang,katom1,katom2,katom3,jtpatm,i,l
  Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,rrr0, &
                            ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

  Integer,           Allocatable :: read_type(:)
  Real( Kind = wp ), Allocatable :: bufpot(:),bufvir(:)

  Character( Len = 256 ) :: message

  If (comm%idnode == 0) Open(Unit=ntable, File='TABANG')

! skip header record

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

! read mesh resolution not needed for angle dependent
! potentials/forces as delpot=180/ngrid running from 0 to 180

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABANG if it's in .xvg format

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  delpot = 180.0_wp/Real(ngrid,wp)

  dlrpot = 180.0_wp/Real(mxgang-4,wp)

! check grid spacing

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > delth_max .and. (.not.safe)) Then
     If (comm%idnode == 0) Then
        Write(nrite,"(/,                                              &
             & ' expected (maximum) angular increment : ',1p,e15.7,/, &
             & ' TABANG file actual angular increment : ',1p,e15.7)") &
             delth_max, delpot
        Write(nrite,"(/,                                                &
             & ' expected (minimum) number of grid points : ',0p,i10,/, &
             & ' TABANG file actual number of grid points : ',0p,i10)") &
             mxgang-4, ngrid
     End If

     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     If (comm%idnode == 0) Write(nrite,"(/,' TABANG arrays resized for mxgrid = ',i10)") mxgang-4
  End If

! compare grids dimensions

  If (ngrid < mxgang-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgang-4,wp),0.0_wp)
     Call error(48)
  End If

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  fail=0
  Allocate (read_type(1:ltpang(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(/,1x,a)') 'error - angles_table_read allocation failure'
     Call error(0,message)
  End If
  Call allocate_angl_pot_arrays()

  read_type=0 ! initialise read_type
  Do rtang=1,ltpang(0)
     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABANG if it's in .xvg format

     Call get_word(record,atom1)
     Call get_word(record,atom2)
     Call get_word(record,atom3)

     katom1=0
     katom2=0
     katom3=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
        If (atom3 == unqatm(jtpatm)) katom3=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0 .or. katom3 == 0) Then
        If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call error(84)
     End If

! Construct unique name for the tabulated angle

     If (katom1 <= katom3) Then
        idangl = atom1//atom2//atom3
     Else
        idangl = atom3//atom2//atom1
     End If

! read potential arrays if potential is defined

     itang=0
     Do jtang=1,ltpang(0)
        If (angl_name(jtang) == idangl) Then
           Do itang=1,mxtang
              If (ltpang(itang) == jtang) Exit
           End Do
           Exit
        End If
     End Do

     If (itang == 0) Then ! All(angl_name /= idangl)
        If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call error(83)
     End If
     If (Any(read_type == jtang)) Then
        If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call error(172)
     Else
        read_type(jtang)=jtang
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
              If (comm%idnode == 0) Write(nrite,"(/,                       &
                 & ' TABANG stated  angular increment : ',1p,e15.7,/, &
                 & ' TABANG read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0
           rrr0      = rrr

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (comm%idnode == 0) Write(nrite,"(/,                       &
                 & ' TABANG stated  angular increment : ',1p,e15.7,/, &
                 & ' TABANG read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr-rrr0
           End If

           bufpot(2) = bufp0
           bufvir(2) = bufv0

! linear extrapolation for grid point 0 at distances close to 0

           bufpot(0) = 2.0_wp*bufpot(1)-bufpot(2)
           bufvir(0) = (2.0_wp*bufvir(1)-0.5_wp*bufvir(2))/dlrpot
        Else ! zero element data found => read in the first element for checking delr
           bufpot(0) = bufp0
           bufvir(0) = bufv0 ! virial/angle - finite force!

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (comm%idnode == 0) Write(nrite,"(/,                       &
                 & ' TABANG stated  angular increment : ',1p,e15.7,/, &
                 & ' TABANG read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
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
        Do i=1,mxgang-4
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
           vang(i,jtang) = t1 + (t2-t1)*ppp*0.5_wp
           vang(i,jtang) = vang(i,jtang)*engunit ! convert to internal units

           vk  = bufvir(l)

! linear extrapolation for the grid points just beyond the cutoff

           If (l+2 > ngrid) Then
              If (l+1 > ngrid) Then
                 vk1 = 2.0_wp*bufvir(l)-bufvir(l-1)
                 vk2 = 2.0_wp*vk1-bufvir(l)
              Else
                 vk1 = bufvir(l+1)
                 vk2 = 2.0_wp*bufvir(l+1)-bufvir(l)
              End If
           Else
              vk1 = bufvir(l+1)
              vk2 = bufvir(l+2)
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           gang(i,jtang) = t1 + (t2-t1)*ppp*0.5_wp
           gang(i,jtang) = gang(i,jtang)*engunit*rad2dgr ! convert to internal units
        End Do

        gang(-1,jtang) = rad2dgr/dlrpot
     Else
        Do i=1,mxgang-4
           vang(i,jtang) = bufpot(i)*engunit         ! convert to internal units
           gang(i,jtang) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do

! linear extrapolation for the grid point just beyond the cutoff

        vang(mxgang-3,jtang) = 2.0_wp*vang(mxgang-4,jtang) - vang(mxgang-5,jtang)
        gang(mxgang-3,jtang) = 2.0_wp*gang(mxgang-4,jtang) - gang(mxgang-5,jtang)

        gang(-1,jtang) = rad2dgr/delpot
     End If

! grid point at 0 and linear extrapolation for the grid point at mxgang-2

     vang(0,jtang) = bufpot(0)
     gang(0,jtang) = bufvir(0)

     vang(mxgang-2,jtang) = 2.0_wp*vang(mxgang-3,jtang) - vang(mxgang-4,jtang)
     gang(mxgang-2,jtang) = 2.0_wp*gang(mxgang-3,jtang) - gang(mxgang-4,jtang)
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABANG file'
  End If

! Break if not safe

  Call gcheck(comm,safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(/,1x,a)') 'error - angles_table_read deallocation failure'
     Call error(0,message)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine angles_table_read



End Module angles
