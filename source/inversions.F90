Module inversions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global valence invle interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms, Only : comms_type,gsum,gsync,gcheck,gbcast
  Use setup_module, Only : mxtmls,mxtinv,mxinv,mxfinv,mxpinv,mxginv,mxginv1, &
                           mxatdm,pi,boltz,delth_max,nrite,npdfdt,npdgdt, &
                           engunit,zero_plus,ntable

  Implicit None

  Logical,                        Save :: lt_inv=.false. ! no tabulated potentials opted

  Integer,                        Save :: ntinv  = 0 , &
                                          ntinv1 = 0 , &
                                          ncfinv = 0


  Integer,           Allocatable, Save :: numinv(:),keyinv(:)
  Integer,           Allocatable, Save :: lstinv(:,:),listinv(:,:),leginv(:,:)

  Real( Kind = wp ), Allocatable, Save :: prminv(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpinv(:)
  Real( Kind = wp ), Allocatable, Save :: vinv(:,:),ginv(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfinv(:),typinv(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstinv(:,:)

  Public :: allocate_inversions_arrays , deallocate_inversions_arrays , &
            allocate_invr_pot_arrays , allocate_invr_dst_arrays, &
            inversions_compute, inversions_forces

Contains

  Subroutine allocate_inversions_arrays()

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numinv(1:mxtmls),          Stat = fail(1))
    Allocate (keyinv(1:mxtinv),          Stat = fail(2))
    Allocate (lstinv(1:4,1:mxtinv),      Stat = fail(3))
    Allocate (listinv(0:4,1:mxinv),      Stat = fail(4))
    Allocate (leginv(0:mxfinv,1:mxatdm), Stat = fail(5))
    Allocate (prminv(1:mxpinv,1:mxtinv), Stat = fail(6))
    If (lt_inv) &
    Allocate (ltpinv(0:mxtinv),          Stat = fail(7))
    If (mxginv1 > 0) &
    Allocate (ldfinv(0:mxtinv),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1021)

    numinv  = 0
    keyinv  = 0
    lstinv  = 0
    listinv = 0
    leginv  = 0

    prminv  = 0.0_wp

    If (lt_inv) &
    ltpinv  = 0

    If (mxginv1 > 0) &
    ldfinv  = 0

  End Subroutine allocate_inversions_arrays

  Subroutine deallocate_inversions_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numinv,lstinv, Stat = fail)

    If (fail > 0) Call error(1034)

  End Subroutine deallocate_inversions_arrays

  Subroutine allocate_invr_pot_arrays()

    Integer :: fail(1:2)

    fail = 0

    Allocate (vinv(-1:mxginv,1:ltpinv(0)), Stat = fail(1))
    Allocate (ginv(-1:mxginv,1:ltpinv(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1078)

    vinv = 0.0_wp
    ginv = 0.0_wp

  End Subroutine allocate_invr_pot_arrays

  Subroutine allocate_invr_dst_arrays()

    Integer :: fail

    fail = 0

    Allocate (typinv(-1:4,1:ldfinv(0)),dstinv(1:mxginv1,1:ldfinv(0)), Stat = fail)

    If (fail > 0) Call error(1079)

    typinv = 0
    dstinv = 0.0_wp

  End Subroutine allocate_invr_dst_arrays

  Subroutine inversions_compute(temp, comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating inversions distribution functions
  ! from accumulated data
  !
  ! copyright - daresbury laboratory
  ! author    - a.v.brukhno & i.t.todorov march 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use site_module,   Only : unqatm
    Use configuration, Only : cfgname

    Real( Kind = wp ),  Intent( In    ) :: temp
    Type( comms_type ), Intent( InOut ) :: comm

    Logical           :: zero
    Integer           :: fail,ngrid,i,j,ig,kk,ll
    Real( Kind = wp ) :: kT2engo,rad2dgr,dgr2rad,delth,rdlth,dgrid,factor,pdfzero, &
                         factor1,theta,pdfinv,sum,pdfinv1,sum1,                    &
                         fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

    Real( Kind = wp ), Allocatable :: dstdinv(:,:)
    Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

    fail = 0
    Allocate (dstdinv(0:mxginv1,1:ldfinv(0)),pmf(0:mxginv1+2),vir(0:mxginv1+2), Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'inversions_compute - allocation failure, node: ', comm%idnode
       Call error(0)
    End If

  ! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

    kT2engo = boltz*temp/engunit

  ! conversion: radians <-> degrees (affects not only angle units but also force units!)

    rad2dgr = 180.0_wp/pi
    dgr2rad = pi/180.0_wp

  ! grid interval for pdf tables

    delth = pi/Real(mxginv1,wp)
    rdlth = Real(mxginv1,wp)/180.0_wp

  ! resampling grid and grid interval for pmf tables

    ngrid = Max(Nint(360.0_wp/delth_max),mxginv1,mxginv-4)
    dgrid = pi/Real(ngrid,wp)

  ! loop over all valid PDFs to get valid totals

    kk=0
    ll=0
    Do i=1,ldfinv(0)
       If (typinv(0,i) > 0) Then
          kk=kk+1
          ll=ll+typinv(0,i)
       End If
    End Do

  ! normalisation factor

    factor = 1.0_wp/Real(ncfinv,wp)

  ! the lower bound to nullify the nearly-zero histogram (PDF) values

    pdfzero = 1.0e-5_wp

    If (comm%idnode == 0) Then
       Write(nrite,'(/,/,12x,a)') 'INVERSIONS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)'
       Write(nrite,'(/,1x,a,5(1x,i10))') &
              '# bins, cutoff, frames, types: ',mxginv1,180,ncfinv,kk,ll
    End If

  ! open RDF file and write headers

    If (comm%idnode == 0) Then
       Open(Unit=npdfdt, File='INVDAT', Status='replace')
       Write(npdfdt,'(a)') '# '//cfgname
       Write(npdfdt,'(a)') '# INVERSIONS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
       Write(npdfdt,'(a,4(1x,i10))') '# bins, cutoff, frames, types: ',mxginv1,180,ncfinv,kk
       Write(npdfdt,'(a)') '#'
       Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)   @   dTheta_bin = ',delth*rad2dgr
       Write(npdfdt,'(a)') '#'
    End If

  ! loop over all valid PDFs

    j=0
    Do i=1,ldfinv(0)
       If (typinv(0,i) > 0) Then
          j=j+1

          If (comm%idnode == 0) Then
             Write(nrite,'(/,1x,a,4(a8,1x),2(i10,1x))') 'type, index, instances: ', &
                  unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)
             Write(nrite,'(/,1x,a,f8.5)') 'Theta(degrees)  P_inv(Theta)  Sum_P_inv(Theta)   @   dTheta_bin = ',delth*rad2dgr

             Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                  unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)),unqatm(typinv(4,i)),j,typinv(0,i)
          End If

  ! global sum of data on all nodes

          Call gsum(comm,dstinv(1:mxginv1,i))

  ! factor in instances (first, pdfinv is normalised to unity)

          factor1=factor/Real(typinv(0,i),wp)

  ! running integration of pdf

          sum=0.0_wp

  ! loop over distances

          zero=.true.
          Do ig=1,mxginv1
             If (zero .and. ig < (mxginv1-3)) zero=(dstinv(ig+2,i) <= 0.0_wp)

             pdfinv = dstinv(ig,i)*factor1
             sum = sum + pdfinv

  ! null it if < pdfzero

             If (pdfinv < pdfzero) Then
                pdfinv1 = 0.0_wp
             Else
                pdfinv1 = pdfinv
             End If

             If (sum < pdfzero) Then
                sum1 = 0.0_wp
             Else
                sum1 = sum
             End If

             theta = (Real(ig,wp)-0.5_wp)*delth

  ! now pdfinv is normalised by the volume element (as to go to unity at infinity in gases and liquids)

             pdfinv = pdfinv*rdlth

  ! print out information

             theta  = theta*rad2dgr
             If (comm%idnode == 0) Then
                If (.not.zero) Write(nrite,"(f11.5,1p,2e14.6)") theta,pdfinv1,sum1
                Write(npdfdt,"(f11.5,1p,e14.6)") theta,pdfinv
             End If

  ! We use the non-normalised tail-truncated PDF version,
  ! pdf...1 (not pdf...) in order to exclude the nearly-zero
  ! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
  ! would have poorly-defined noisy "borders/walls"

             dstdinv(ig,i) = pdfinv1 ! PDFs density
          End Do
       Else
          dstdinv(:,i) = 0.0_wp ! PDFs density
       End If
    End Do

    If (comm%idnode == 0) Close(Unit=npdfdt)

  ! open PDF files and write headers

    If (comm%idnode == 0) Then
       Open(Unit=npdgdt, File='INVPMF', Status='replace')
       Write(npdgdt,'(a)') '# '//cfgname
       Write(npdgdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',mxginv1,delth*Real(mxginv1,wp)*rad2dgr,delth*rad2dgr,kk, &
            '   conversion factor(kT -> energy units) =',kT2engo

       Open(Unit=npdfdt, File='INVPMF', Status='replace')
       Write(npdfdt,'(a)') '# '//cfgname
       Write(npdfdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',ngrid,dgrid*Real(ngrid,wp)*rad2dgr,dgrid*rad2dgr,kk, &
            '   conversion factor(kT -> energy units) =',kT2engo
    End If

  ! loop over all valid PDFs

    j=0
    Do i=1,ldfinv(0)
       If (typinv(0,i) > 0) Then
          j=j+1

          If (comm%idnode == 0) Then
             Write(npdgdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                  unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)), &
                  unqatm(typinv(4,i)),j,typinv(0,i),' (type, index, instances)'
             Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                  unqatm(typinv(1,i)),unqatm(typinv(2,i)),unqatm(typinv(3,i)), &
                  unqatm(typinv(4,i)),j,typinv(0,i),' (type, index, instances)'
          End If

  ! Smoothen and get derivatives

          fed0  = 0.0_wp
          dfed0 = 10.0_wp
          dfed  = 10.0_wp

          Do ig=1,mxginv1
             tmp = Real(ig,wp)-0.5_wp
             theta = tmp*delth

             If (dstdinv(ig,i) > zero_plus) Then
                fed = -Log(dstdinv(ig,i))-fed0
                If (fed0 <= zero_plus) Then
                   fed0 = fed
                   fed  = 0.0_wp
                End If

                If (ig < mxginv1-1) Then
                   If (dstdinv(ig+1,i) <= zero_plus .and. dstdinv(ig+2,i) > zero_plus) &
                      dstdinv(ig+1,i) = 0.5_wp*(dstdinv(ig,i)+dstdinv(ig+2,i))
                End If
             Else
                fed = 0.0_wp
             End If

             If      (ig == 1) Then
                If      (dstdinv(ig,i) > zero_plus .and. dstdinv(ig+1,i) > zero_plus) Then
                   dfed = Log(dstdinv(ig+1,i)/dstdinv(ig,i))
                Else If (dfed > 0.0_wp) Then
                   dfed = dfed0
                Else
                   dfed =-dfed0
                End If
             Else If (ig == mxginv1) Then
                If      (dstdinv(ig,i) > zero_plus .and. dstdinv(ig-1,i) > zero_plus) Then
                   dfed = Log(dstdinv(ig,i)/dstdinv(ig-1,i))
                Else If (dfed > 0.0_wp) Then
                   dfed = dfed0
                Else
                   dfed =-dfed0
                End If
             Else If (dstdinv(ig-1,i) > zero_plus) Then
                If (dstdinv(ig+1,i) > zero_plus) Then
                   dfed = 0.5_wp*(Log(dstdinv(ig+1,i)/dstdinv(ig-1,i)))
                Else
                   dfed = 0.5_wp*Log(dstdinv(ig-1,i))
                End If
             Else If (dstdinv(ig+1,i) > zero_plus) Then
                dfed =-0.5_wp*Log(dstdinv(ig+1,i))
             Else If (dfed > 0.0_wp) Then
                dfed = dfed0
             Else
                dfed =-dfed0
             End If

             pmf(ig) = fed
             vir(ig) = dfed

  ! Print

             If (comm%idnode == 0) Write(npdgdt,"(f11.5,1p,2e14.6)") theta*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
          End Do

  ! Define edges

          pmf(0)         = 2.0_wp*pmf(1)        -pmf(2)
          vir(0)         = 2.0_wp*vir(1)        -vir(2)
          pmf(mxginv1+1) = 2.0_wp*pmf(mxginv1)  -pmf(mxginv1-1)
          vir(mxginv1+1) = 2.0_wp*vir(mxginv1)  -vir(mxginv1-1)
          pmf(mxginv1+2) = 2.0_wp*pmf(mxginv1+1)-pmf(mxginv1)
          vir(mxginv1+2) = 2.0_wp*vir(mxginv1+1)-vir(mxginv1)

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

    Deallocate (dstdinv,pmf,vir, Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'inversions_compute - deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

  End Subroutine inversions_compute

  Subroutine inversions_forces(isw,enginv,virinv,stress,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating inversion energy and force terms
  !
  ! isw = 0 - collect statistics
  ! isw = 1 - calculate forces
  ! isw = 2 - do both
  !
  ! copyright - daresbury laboratory
  ! author    - w.smith may 1996
  ! amended   - i.t.todorov february 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use configuration,     Only : imcon,cell,natms,nlast,lsi,lsa,lfrzn, &
                                  xxx,yyy,zzz,fxx,fyy,fzz

    Integer,                             Intent( In    ) :: isw
    Real( Kind = wp ),                   Intent(   Out ) :: enginv,virinv
    Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
    Type( comms_type ),                  Intent( InOut ) :: comm

    Logical           :: safe
    Integer           :: fail(1:4),i,j,l,ia,ib,ic,id,kk,keyi,local_index
    Real( Kind = wp ) :: xab,yab,zab,rab2,rrab, xac,yac,zac,rac2,rrac,   &
                         xad,yad,zad,rad2,rrad, rbc,rcd,rdb,             &
                         ubx,uby,ubz,ubn,rub, vbx,vby,vbz,vbn,rvb,wwb,   &
                         ucx,ucy,ucz,ucn,ruc, vcx,vcy,vcz,vcn,rvc,wwc,   &
                         udx,udy,udz,udn,rud, vdx,vdy,vdz,vdn,rvd,wwd,   &
                         cosb,cosc,cosd, rubc,rubd,rucd,rucb,rudb,rudc,  &
                         rvbc,rvbd,rvcd,rvcb,rvdb,rvdc,                  &
                         thb,thc,thd,  k,th0,cos0,a,b,m,                 &
                         uuu,uu2,uun,uux,uuy,uuz,                        &
                         pterm,vterm,gamma,gamb,gamc,gamd,               &
                         rdelth,rdr,ppp,vk,vk1,vk2,t1,t2,                &
                         fax,fay,faz, fbx,fby,fbz,                       &
                         fcx,fcy,fcz, fdx,fdy,fdz,                       &
                         strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

    Logical,           Allocatable :: lunsafe(:)
    Integer,           Allocatable :: lstopt(:,:)
    Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
    Real( Kind = wp ), Allocatable :: xdac(:),ydac(:),zdac(:)
    Real( Kind = wp ), Allocatable :: xdad(:),ydad(:),zdad(:)

    fail=0
    Allocate (lunsafe(1:mxinv),lstopt(0:4,1:mxinv),      Stat=fail(1))
    Allocate (xdab(1:mxinv),ydab(1:mxinv),zdab(1:mxinv), Stat=fail(2))
    Allocate (xdac(1:mxinv),ydac(1:mxinv),zdac(1:mxinv), Stat=fail(3))
    Allocate (xdad(1:mxinv),ydad(1:mxinv),zdad(1:mxinv), Stat=fail(4))
    If (Any(fail > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'inversions_forces allocation failure, node: ', comm%idnode
       Call error(0)
    End If


  ! calculate atom separation vectors

    Do i=1,ntinv
       lunsafe(i)=.false.

  ! indices of atoms involved

       ia=local_index(listinv(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
       ib=local_index(listinv(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib
       ic=local_index(listinv(3,i),nlast,lsi,lsa) ; lstopt(3,i)=ic
       id=local_index(listinv(4,i),nlast,lsi,lsa) ; lstopt(4,i)=id

       lstopt(0,i)=0
       If (ia > 0 .and. ib > 0 .and. ic > 0 .and. id > 0) Then !Tag
          If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic)*lfrzn(id) == 0) Then
             If (ia <= natms .or. ib <= natms .or. ic <= natms .or. id <= natms) Then
                lstopt(0,i)=1
              End If
          End If
       Else                                                    ! Detect uncompressed unit
          If ( ((ia > 0 .and. ia <= natms) .or.   &
                (ib > 0 .and. ib <= natms) .or.   &
                (ic > 0 .and. ic <= natms) .or.   &
                (id > 0 .and. id <= natms)) .and. &
               (ia == 0 .or. ib == 0 .or. ic == 0 .or. id == 0) ) lunsafe(i)=.true.
       End If

  ! define components of bond vectors

       If (lstopt(0,i) > 0) Then
          xdab(i)=xxx(ib)-xxx(ia)
          ydab(i)=yyy(ib)-yyy(ia)
          zdab(i)=zzz(ib)-zzz(ia)

  ! select potential energy function type

          kk=listinv(0,i)
          keyi = Abs(keyinv(kk))

          If (keyi == 5) Then
             xdac(i)=xxx(ic)-xxx(ib)
             ydac(i)=yyy(ic)-yyy(ib)
             zdac(i)=zzz(ic)-zzz(ib)

             xdad(i)=xxx(id)-xxx(ib)
             ydad(i)=yyy(id)-yyy(ib)
             zdad(i)=zzz(id)-zzz(ib)
          Else
             xdac(i)=xxx(ic)-xxx(ia)
             ydac(i)=yyy(ic)-yyy(ia)
             zdac(i)=zzz(ic)-zzz(ia)

             xdad(i)=xxx(id)-xxx(ia)
             ydad(i)=yyy(id)-yyy(ia)
             zdad(i)=zzz(id)-zzz(ia)
          End If
  !     Else ! (DEBUG)
  !        xdab(i)=0.0_wp
  !        ydab(i)=0.0_wp
  !        zdab(i)=0.0_wp
  !
  !        xdac(i)=0.0_wp
  !        ydac(i)=0.0_wp
  !        zdac(i)=0.0_wp
  !
  !        xdad(i)=0.0_wp
  !        ydad(i)=0.0_wp
  !        zdad(i)=0.0_wp
       End If
    End Do

  ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:ntinv))
    Call gcheck(comm,safe)
    If (.not.safe) Then
       Do j=0,comm%mxnode-1
          If (comm%idnode == j) Then
             Do i=1,ntinv
                If (lunsafe(i)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                   '*** warning - global unit number', listinv(0,i), &
                   ' , with a head particle number', listinv(1,i),   &
                   ' contributes towards next error !!! ***'
             End Do
          End If
          Call gsync(comm)
       End Do
       Call error(134)
    End If

  ! periodic boundary condition

    Call images(imcon,cell,ntinv,xdab,ydab,zdab)
    Call images(imcon,cell,ntinv,xdac,ydac,zdac)
    Call images(imcon,cell,ntinv,xdad,ydad,zdad)

    If (Mod(isw,3) > 0) Then

  ! Initialise safety flag

       safe=.true.

  ! zero inversion energy accumulator

       enginv=0.0_wp
       virinv=0.0_wp

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
       rdelth = Real(mxginv1,wp)/pi
       ncfinv = ncfinv + 1
    End If

  ! loop over all specified inversions

    Do i=1,ntinv
       If (lstopt(0,i) > 0) Then

  ! indices of atoms involved

          ia=lstopt(1,i)
          ib=lstopt(2,i)
          ic=lstopt(3,i)
          id=lstopt(4,i)

  ! define components of bond vectors

          xab=xdab(i)
          yab=ydab(i)
          zab=zdab(i)
          rab2=xab*xab+yab*yab+zab*zab
          rrab=1.0_wp/Sqrt(rab2)

          xac=xdac(i)
          yac=ydac(i)
          zac=zdac(i)
          rac2=xac*xac+yac*yac+zac*zac
          rrac=1.0_wp/Sqrt(rac2)

          xad=xdad(i)
          yad=ydad(i)
          zad=zdad(i)
          rad2=xad*xad+yad*yad+zad*zad
          rrad=1.0_wp/Sqrt(rad2)

  ! select potential energy function type

          kk=listinv(0,i)
          keyi=keyinv(kk)

          If (keyi == 5) Then

  ! calculate vector normal to plane

             uux=yac*zad-zac*yad
             uuy=zac*xad-xac*zad
             uuz=xac*yad-yac*xad
             uun=1.0_wp/Sqrt(uux**2+uuy**2+uuz**2)
             uux=uun*uux
             uuy=uun*uuy
             uuz=uun*uuz
             uuu=xab*uux+yab*uuy+zab*uuz

          Else

  ! scalar products of bond vectors

             rbc=xab*xac+yab*yac+zab*zac
             rcd=xac*xad+yac*yad+zac*zad
             rdb=xad*xab+yad*yab+zad*zab

  ! calculate bond-angle-plane vectors

             ubx=xac*rrac+xad*rrad
             uby=yac*rrac+yad*rrad
             ubz=zac*rrac+zad*rrad
             ubn=1.0_wp/Sqrt(ubx**2+uby**2+ubz**2)
             ubx=ubn*ubx
             uby=ubn*uby
             ubz=ubn*ubz
             rub=xab*ubx+yab*uby+zab*ubz

             vbx=xac*rrac-xad*rrad
             vby=yac*rrac-yad*rrad
             vbz=zac*rrac-zad*rrad
             vbn=1.0_wp/Sqrt(vbx**2+vby**2+vbz**2)
             vbx=vbn*vbx
             vby=vbn*vby
             vbz=vbn*vbz
             rvb=xab*vbx+yab*vby+zab*vbz
             wwb=Sqrt(rub**2+rvb**2)

             ucx=xad*rrad+xab*rrab
             ucy=yad*rrad+yab*rrab
             ucz=zad*rrad+zab*rrab
             ucn=1.0_wp/Sqrt(ucx**2+ucy**2+ucz**2)
             ucx=ucn*ucx
             ucy=ucn*ucy
             ucz=ucn*ucz
             ruc=xac*ucx+yac*ucy+zac*ucz

             vcx=xad*rrad-xab*rrab
             vcy=yad*rrad-yab*rrab
             vcz=zad*rrad-zab*rrab
             vcn=1.0_wp/Sqrt(vcx**2+vcy**2+vcz**2)
             vcx=vcn*vcx
             vcy=vcn*vcy
             vcz=vcn*vcz
             rvc=xac*vcx+yac*vcy+zac*vcz
             wwc=Sqrt(ruc**2+rvc**2)

             udx=xab*rrab+xac*rrac
             udy=yab*rrab+yac*rrac
             udz=zab*rrab+zac*rrac
             udn=1.0_wp/Sqrt(udx**2+udy**2+udz**2)
             udx=udn*udx
             udy=udn*udy
             udz=udn*udz
             rud=xad*udx+yad*udy+zad*udz

             vdx=xab*rrab-xac*rrac
             vdy=yab*rrab-yac*rrac
             vdz=zab*rrab-zac*rrac
             vdn=1.0_wp/Sqrt(vdx**2+vdy**2+vdz**2)
             vdx=vdn*vdx
             vdy=vdn*vdy
             vdz=vdn*vdz
             rvd=xad*vdx+yad*vdy+zad*vdz
             wwd=Sqrt(rud**2+rvd**2)

  ! calculate inversion angle cosines

             cosb=wwb*rrab ; If (Abs(cosb) > 1.0_wp) cosb=Sign(1.0_wp,cosb)
             cosc=wwc*rrac ; If (Abs(cosc) > 1.0_wp) cosc=Sign(1.0_wp,cosb)
             cosd=wwd*rrad ; If (Abs(cosd) > 1.0_wp) cosd=Sign(1.0_wp,cosb)

  ! accumulate the histogram (distribution)

             If (Mod(isw,2) == 0 .and. ib <= natms) Then
                j = ldfinv(kk)

                thb=Acos(cosb)
                l = Min(1+Int(thb*rdelth),mxginv1)
                dstinv(l,j) = dstinv(l,j) + 1.0_wp/3.0_wp

                thc=Acos(cosc)
                l = Min(1+Int(thc*rdelth),mxginv1)
                dstinv(l,j) = dstinv(l,j) + 1.0_wp/3.0_wp

                thd=Acos(cosd)
                l = Min(1+Int(thd*rdelth),mxginv1)
                dstinv(l,j) = dstinv(l,j) + 1.0_wp/3.0_wp
             End If

          End If
          If (isw == 0) Cycle

  ! calculate potential energy and scalar force term

          If      (keyi == 1) Then

  ! harmonic inversion potential

             k  =prminv(1,kk)/6.0_wp
             th0=prminv(2,kk)

             thb=Acos(cosb)
             thc=Acos(cosc)
             thd=Acos(cosd)

             pterm=k*((thb-th0)**2+(thc-th0)**2+(thd-th0)**2)
             vterm=0.0_wp
             gamma=0.0_wp
             gamb=0.0_wp
             If (Abs(thb) > 1.0e-10_wp) gamb=2.0_wp*k*(thb-th0)/Sin(thb)
             gamc=0.0_wp
             If (Abs(thc) > 1.0e-10_wp) gamc=2.0_wp*k*(thc-th0)/Sin(thc)
             gamd=0.0_wp
             If (Abs(thd) > 1.0e-10_wp) gamd=2.0_wp*k*(thd-th0)/Sin(thd)

          Else If (keyi == 2) Then

  ! harmonic cosine inversion potential

             k   =prminv(1,kk)/6.0_wp
             cos0=prminv(2,kk)

             pterm=k*((cosb-cos0)**2+(cosc-cos0)**2+(cosd-cos0)**2)
             vterm=0.0_wp
             gamma=0.0_wp
             gamb=-2.0_wp*k*(cosb-cos0)
             gamc=-2.0_wp*k*(cosc-cos0)
             gamd=-2.0_wp*k*(cosd-cos0)

          Else If (keyi == 3) Then

  ! planar inversion potentials

             a=prminv(1,kk)

             pterm=a*(1.0_wp-(cosb+cosc+cosd)/3.0_wp)
             vterm=0.0_wp
             gamma=0.0_wp
             gamb=a/3.0_wp
             gamc=a/3.0_wp
             gamd=a/3.0_wp

          Else If (keyi == 4) Then

  ! extended planar inversion potentials

             k  =prminv(1,kk)/6.0_wp
             th0=prminv(2,kk)
             m  =prminv(3,kk)

             thb=Acos(cosb)
             thc=Acos(cosc)
             thd=Acos(cosd)

             pterm=3.0_wp*k*(1.0_wp-(Cos(m*thb-th0)+Cos(m*thc-th0)+Cos(m*thd-th0))/3.0_wp)
             vterm=0.0_wp
             gamma=0.0_wp
             gamb=0.0_wp
             If (Abs(thb) > 1.0e-10_wp) gamb=k*Sin(m*thb-th0)/Sin(thb)
             gamc=0.0_wp
             If (Abs(thc) > 1.0e-10_wp) gamc=k*Sin(m*thc-th0)/Sin(thc)
             gamd=0.0_wp
             If (Abs(thd) > 1.0e-10_wp) gamd=k*Sin(m*thd-th0)/Sin(thd)

          Else If (keyi == 5) Then

  ! planar calcite potential

             a=prminv(1,kk)
             b=prminv(2,kk)

             uu2=uuu*uuu
             m=2.0_wp*a+4.0_wp*b*uu2

             pterm=uu2*(a+b*uu2)
             vterm=uu2*m
             gamma=-uuu*m
             gamb=0.0_wp
             gamc=0.0_wp
             gamd=0.0_wp

          Else If (keyi == 20) Then

  ! TABINV potential

             pterm=0.0_wp

             j = ltpinv(kk)
             rdr = ginv(-1,j) ! 1.0_wp/delpot (in rad^-1)

             thb=Acos(cosb)
             thc=Acos(cosc)
             thd=Acos(cosd)

             l   = Int(thb*rdr)
             ppp = thb*rdr - Real(l,wp)

             vk  = vinv(l,j)
             vk1 = vinv(l+1,j)
             vk2 = vinv(l+2,j)

             t1 = vk  + (vk1 - vk)*ppp
             t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

             pterm = pterm + t1 + (t2-t1)*ppp*0.5_wp

             vk  = ginv(l,j) ; If (l == 0) vk = vk*thb
             vk1 = ginv(l+1,j)
             vk2 = ginv(l+2,j)

             t1 = vk  + (vk1 - vk)*ppp
             t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

             gamb = -t1 + (t2-t1)*ppp*0.5_wp

             l   = Int(thc*rdr)
             ppp = thc*rdr - Real(l,wp)

             vk  = vinv(l,j)
             vk1 = vinv(l+1,j)
             vk2 = vinv(l+2,j)

             t1 = vk  + (vk1 - vk)*ppp
             t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

             pterm = pterm + t1 + (t2-t1)*ppp*0.5_wp

             vk  = ginv(l,j) ; If (l == 0) vk = vk*thc
             vk1 = ginv(l+1,j)
             vk2 = ginv(l+2,j)

             t1 = vk  + (vk1 - vk)*ppp
             t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

             gamc = -t1 + (t2-t1)*ppp*0.5_wp

             l   = Int(thd*rdr)
             ppp = thd*rdr - Real(l,wp)

             vk  = Merge(vinv(l,j), 0.0_wp, l > 0)
             vk1 = vinv(l+1,j)
             vk2 = vinv(l+2,j)

             t1 = vk  + (vk1 - vk)*ppp
             t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

             pterm = pterm + t1 + (t2-t1)*ppp*0.5_wp

             vk  = ginv(l,j) ; If (l == 0) vk = vk*thd
             vk1 = ginv(l+1,j)
             vk2 = ginv(l+2,j)

             t1 = vk  + (vk1 - vk)*ppp
             t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

             gamd = -t1 + (t2-t1)*ppp*0.5_wp

             vterm=0.0_wp
             gamma=0.0_wp

          Else

  ! undefined potential

             safe=.false.
             pterm=0.0_wp
             vterm=0.0_wp
             gamb=0.0_wp
             gamc=0.0_wp
             gamd=0.0_wp

          End If

          If (keyinv(kk) == 5) Then

  ! calculate atomic forces

             fax=-gamma*uux
             fay=-gamma*uuy
             faz=-gamma*uuz

             fcx=gamma*uun*((yad*zab-zad*yab)-uuu*(yad*uuz-zad*uuy))
             fcy=gamma*uun*((zad*xab-xad*zab)-uuu*(zad*uux-xad*uuz))
             fcz=gamma*uun*((xad*yab-yad*xab)-uuu*(xad*uuy-yad*uux))

             fdx=gamma*uun*((yab*zac-zab*yac)-uuu*(zac*uuy-yac*uuz))
             fdy=gamma*uun*((zab*xac-xab*zac)-uuu*(xac*uuz-zac*uux))
             fdz=gamma*uun*((xab*yac-yab*xac)-uuu*(yac*uux-xac*uuy))

             fbx=-(fax+fcx+fdx)
             fby=-(fay+fcy+fdy)
             fbz=-(faz+fcz+fdz)

          Else

  ! calculate bond and u,v scalar products

             rubc=xab*ucx+yab*ucy+zab*ucz
             rubd=xab*udx+yab*udy+zab*udz
             rucd=xac*udx+yac*udy+zac*udz
             rucb=xac*ubx+yac*uby+zac*ubz
             rudb=xad*ubx+yad*uby+zad*ubz
             rudc=xad*ucx+yad*ucy+zad*ucz

             rvbc=xab*vcx+yab*vcy+zab*vcz
             rvbd=xab*vdx+yab*vdy+zab*vdz
             rvcd=xac*vdx+yac*vdy+zac*vdz
             rvcb=xac*vbx+yac*vby+zac*vbz
             rvdb=xad*vbx+yad*vby+zad*vbz
             rvdc=xad*vcx+yad*vcy+zad*vcz

  ! calculate atomic forces

             fbx = gamb*(-cosb*xab*rrab**2+rrab*(rub*ubx+rvb*vbx)/wwb) +       &
                   ( ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2) -   &
                     rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2) ) * &
                   gamc*rrac/wwc +                                             &
                   ( rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2) +   &
                     rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2) ) * &
                   gamd*rrad/wwd

             fby = gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb) +       &
                   ( ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2) -   &
                     rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2) ) * &
                   gamc*rrac/wwc +                                             &
                   ( rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2) +   &
                     rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2) ) * &
                   gamd*rrad/wwd

             fbz = gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb) +       &
                   ( ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2) -   &
                     rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2) ) * &
                   gamc*rrac/wwc +                                             &
                   ( rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2) +   &
                     rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2) ) * &
                   gamd*rrad/wwd

             fcx = gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc) +       &
                   ( rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2) -   &
                     rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2) ) * &
                   gamd*rrad/wwd +                                             &
                   ( rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2) +   &
                     rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2) ) * &
                   gamb*rrab/wwb

             fcy = gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc) +       &
                   ( rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2) -   &
                     rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2) ) * &
                   gamd*rrad/wwd +                                             &
                   ( rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2) +   &
                     rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2) ) * &
                   gamb*rrab/wwb

             fcz = gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc) +       &
                   ( rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2) -   &
                     rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2) ) * &
                   gamd*rrad/wwd +                                             &
                   ( rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2) +   &
                     rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2) ) * &
                   gamb*rrab/wwb

             fdx = gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd) +       &
                   ( rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2) -   &
                     rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2) ) * &
                   gamb*rrab/wwb +                                             &
                   ( ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2) +   &
                     rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2) ) * &
                   gamc*rrac/wwc

             fdy = gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd) +       &
                   ( rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2) -   &
                     rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2) ) * &
                   gamb*rrab/wwb +                                             &
                   ( ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2) +   &
                     rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2) ) * &
                   gamc*rrac/wwc

             fdz = gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd) +       &
                   ( rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2) -   &
                     rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2) ) * &
                   gamb*rrab/wwb +                                             &
                   ( ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2) +   &
                     rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2) ) * &
                   gamc*rrac/wwc

             fax = -(fbx+fcx+fdx)
             fay = -(fby+fcy+fdy)
             faz = -(fbz+fcz+fdz)

          End If

          If (ia <= natms) Then

  ! inversion energy and virial (associated to the head atom)

             enginv=enginv+pterm
             virinv=virinv+vterm

  ! stress tensor calculation for inversion terms

             If (keyinv(kk) == 5) Then
                strs1 = strs1 + uuu*gamma*uux*uux
                strs2 = strs2 + uuu*gamma*uux*uuy
                strs3 = strs3 + uuu*gamma*uux*uuz
                strs5 = strs5 + uuu*gamma*uuy*uuy
                strs6 = strs6 + uuu*gamma*uuy*uuz
                strs9 = strs9 + uuu*gamma*uuz*uuz
             Else
                strs1 = strs1 + xab*fbx + xac*fcx + xad*fdx
                strs2 = strs2 + yab*fbx + yac*fcx + yad*fdx
                strs3 = strs3 + zab*fbx + zac*fcx + zad*fdx
                strs5 = strs5 + yab*fby + yac*fcy + yad*fdy
                strs6 = strs6 + yab*fbz + yac*fcz + yad*fdz
                strs9 = strs9 + zab*fbz + zac*fcz + zad*fdz
             End If

             fxx(ia)=fxx(ia)+fax
             fyy(ia)=fyy(ia)+fay
             fzz(ia)=fzz(ia)+faz

          End If

          If (ib <= natms) Then

             fxx(ib)=fxx(ib)+fbx
             fyy(ib)=fyy(ib)+fby
             fzz(ib)=fzz(ib)+fbz

          End If

          If (ic <= natms) Then

             fxx(ic)=fxx(ic)+fcx
             fyy(ic)=fyy(ic)+fcy
             fzz(ic)=fzz(ic)+fcz

          End If

          If (id <= natms) Then

             fxx(id)=fxx(id)+fdx
             fyy(id)=fyy(id)+fdy
             fzz(id)=fzz(id)+fdz

          End If

       End If
    End Do

    If (Mod(isw,3) > 0) Then

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
       If (.not.safe) Call error(449)

  ! global sum of inversion potential and virial

       buffer(1)=enginv
       buffer(2)=virinv
       Call gsum(comm,buffer(1:2))
       enginv=buffer(1)
       virinv=buffer(2)

    End If

    Deallocate (lunsafe,lstopt, Stat=fail(1))
    Deallocate (xdab,ydab,zdab, Stat=fail(2))
    Deallocate (xdac,ydac,zdac, Stat=fail(3))
    Deallocate (xdad,ydad,zdad, Stat=fail(4))
    If (Any(fail > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'inversions_forces deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

  End Subroutine inversions_forces

  Subroutine inversions_table_read(invr_name,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for reading potential energy and force arrays
  ! from TABINV file (for inversion potentials & forces only)
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov january 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use site_module,  Only : ntpatm,unqatm
    Use parse_module, Only : get_line,get_word,word_2_real

    Implicit None

    Character( Len = 32 ), Intent( In    ) :: invr_name(1:mxtinv)
    Type( comms_type ),    Intent( InOut ) :: comm

    Logical                :: safe,remake
    Character( Len = 200 ) :: record
    Character( Len = 40  ) :: word
    Character( Len = 32 )  :: idinvr
    Character( Len = 8   ) :: atom1,atom2,atom3,atom4

    Integer                :: fail(1:2),ngrid,rtinv,itinv,jtinv,katom1,katom2,katom3,katom4,jtpatm,i,l
    Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,rrr0, &
                              ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

    Integer,           Allocatable :: read_type(:)
    Real( Kind = wp ), Allocatable :: bufpot(:),bufvir(:)


    If (comm%idnode == 0) Open(Unit=ntable, File='TABINV')

  ! skip header record

    Call get_line(safe,ntable,record,comm)
    If (.not.safe) Go To 100

  ! read mesh resolution not needed for inversion angle dependent
  ! potentials/forces as delpot=180/ngrid running from 0 to 180

    Call get_line(safe,ntable,record,comm)
    If (.not.safe) Go To 100

    i = Index(record,'#')      ! replace hash as it may occur in
    If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

    Call get_word(record,word)
    ngrid = Nint(word_2_real(word,comm))

    delpot = 180.0_wp/Real(ngrid,wp)

    dlrpot = 180.0_wp/Real(mxginv-4,wp)

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
               & ' TABINV file actual angular increment : ',1p,e15.7)") &
               delth_max, delpot
          Write(nrite,"(/,                                                &
               & ' expected (minimum) number of grid points : ',0p,i10,/, &
               & ' TABINV file actual number of grid points : ',0p,i10)") &
               mxginv-4, ngrid
       End If
       Call error(22)
    End If
    safe=.true.

    remake=.false.
    If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
       remake=.true.
       rdr=1.0_wp/delpot
       If (comm%idnode == 0) Write(nrite,"(/,' TABINV arrays resized for mxgrid = ',i10)") mxginv-4
    End If

  ! compare grids dimensions

    If (ngrid < mxginv-4) Then
       Call warning(270,Real(ngrid,wp),Real(mxginv-4,wp),0.0_wp)
       Call error(48)
    End If

    rad2dgr= 180.0_wp/pi
    dgr2rad= pi/180.0_wp

    fail=0
    Allocate (read_type(1:ltpinv(0)),          Stat=fail(1))
    Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
    If (Any(fail > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'error - inversions_table_read allocation failure, node: ', comm%idnode
       Call error(0)
    End If
    Call allocate_invr_pot_arrays()

    read_type=0 ! initialise read_type
    Do rtinv=1,ltpinv(0)
       Call get_line(safe,ntable,record,comm)
       If (.not.safe) Go To 100

       Call get_line(safe,ntable,record,comm)
       If (.not.safe) Go To 100

       i = Index(record,'#')      ! replace hash as it may occur in
       If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

       Call get_word(record,atom1)
       Call get_word(record,atom2)
       Call get_word(record,atom3)
       Call get_word(record,atom4)

       katom1=0
       katom2=0
       katom3=0
       katom4=0

       Do jtpatm=1,ntpatm
          If (atom1 == unqatm(jtpatm)) katom1=jtpatm
          If (atom2 == unqatm(jtpatm)) katom2=jtpatm
          If (atom3 == unqatm(jtpatm)) katom3=jtpatm
          If (atom4 == unqatm(jtpatm)) katom4=jtpatm
       End Do

       If (katom1 == 0 .or. katom2 == 0 .or. katom3 == 0 .or. katom4 == 0) Then
          If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABINV'
          Call error(91)
       End If

  ! Construct unique name for the tabulated inversion

       If      (Min(katom2,katom3,katom4) == katom2) Then
          If (katom3 <= katom4) Then
             idinvr = atom1//atom2//atom3//atom4
          Else
             idinvr = atom1//atom2//atom4//atom3
          End If
       Else If (Min(katom2,katom3,katom4) == katom3) Then
          If (katom2 <= katom4) Then
             idinvr = atom1//atom3//atom2//atom4
          Else
             idinvr = atom1//atom3//atom4//atom2
          End If
       Else
          If (katom2 <= katom3) Then
             idinvr = atom1//atom4//atom2//atom3
          Else
             idinvr = atom1//atom4//atom3//atom2
          End If
       End If

  ! read potential arrays if potential is defined

       itinv=0
       Do jtinv=1,ltpinv(0)
          If (invr_name(jtinv) == idinvr) Then
             Do itinv=1,mxtinv
                If (ltpinv(itinv) == jtinv) Exit
             End Do
             Exit
          End If
       End Do

       If (itinv == 0) Then ! All(invr_name /= idinvr)
          If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABINV'
          Call error(89)
       End If
       If (Any(read_type == jtinv)) Then
          If (comm%idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABINV'
          Call error(172)
       Else
          read_type(jtinv)=jtinv
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
                   & ' TABINV stated  angular increment : ',1p,e15.7,/, &
                   & ' TABINV read-in angular increment : ',1p,e15.7)") &
                   delpot,rrr
             End If

             bufpot(1) = bufp0
             bufvir(1) = bufv0
             rrr0      = rrr

             Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

             If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
                safe=.false.
                If (comm%idnode == 0) Write(nrite,"(/,                       &
                   & ' TABINV stated  angular increment : ',1p,e15.7,/, &
                   & ' TABINV read-in angular increment : ',1p,e15.7)") &
                   delpot,rrr-rrr0
             End If

             bufpot(2) = bufp0
             bufvir(2) = bufv0

  ! linear extrapolation for grid point 0 at distances close to 0

             bufpot(0) = 2.0_wp*bufpot(1)-bufpot(2)
             bufvir(0) = (2.0_wp*bufvir(1)-0.5_wp*bufvir(2))/dlrpot
          Else ! zero element data found => read in the first element for checking delr
             bufpot(0) = bufp0
             bufvir(0) = bufv0

             Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

             If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
                safe=.false.
                If (comm%idnode == 0) Write(nrite,"(/,                       &
                   & ' TABINV stated  angular increment : ',1p,e15.7,/, &
                   & ' TABINV read-in angular increment : ',1p,e15.7)") &
                   delpot,rrr
             End If

             bufpot(1) = bufp0
             bufvir(1) = bufv0

             Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

             bufpot(2) = bufp0
             bufvir(2) = bufv0
          End If
       End If

  ! read in potential & force arrays

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
          Do i=1,mxginv-4
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
             vinv(i,jtinv) = t1 + (t2-t1)*ppp*0.5_wp
             vinv(i,jtinv) = vinv(i,jtinv)*engunit ! convert to internal units

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
             ginv(i,jtinv) = t1 + (t2-t1)*ppp*0.5_wp
             ginv(i,jtinv) = ginv(i,jtinv)*engunit*rad2dgr ! convert to internal units
          End Do

          ginv(-1,jtinv) = rad2dgr/dlrpot
       Else
          Do i=1,mxginv-4
             vinv(i,jtinv) = bufpot(i)*engunit ! convert to internal units
             ginv(i,jtinv) = bufvir(i)*engunit*rad2dgr ! convert to internal units
          End Do

  ! linear extrapolation for the grid point just beyond the cutoff

          vinv(mxginv-3,jtinv) = 2.0_wp*vinv(mxginv-4,jtinv) - vinv(mxginv-5,jtinv)
          ginv(mxginv-3,jtinv) = 2.0_wp*ginv(mxginv-4,jtinv) - ginv(mxginv-5,jtinv)

          ginv(-1,jtinv) = rad2dgr/delpot
       End If

  ! grid point at 0 and linear extrapolation for the grid point at mxginv-2


       vinv(0,jtinv) = bufpot(0)
       ginv(0,jtinv) = bufvir(0)

       vinv(mxginv-2,jtinv) = 2.0_wp*vinv(mxginv-3,jtinv) - vinv(mxginv-4,jtinv)
       ginv(mxginv-2,jtinv) = 2.0_wp*ginv(mxginv-3,jtinv) - ginv(mxginv-4,jtinv)
    End Do

    If (comm%idnode == 0) Then
       Close(Unit=ntable)
       Write(nrite,'(/,1x,a)') 'potential tables read from TABINV file'
    End If

  ! Break if not safe

    Call gcheck(comm,safe)
    If (.not.safe) Call error(22)

    Deallocate (read_type,     Stat=fail(1))
    Deallocate (bufpot,bufvir, Stat=fail(2))
    If (Any(fail > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'error - inversions_table_read deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

    Return

  ! end of file error exit

  100 Continue

    If (comm%idnode == 0) Close(Unit=ntable)
    Call error(24)

  End Subroutine inversions_table_read
End Module inversions
