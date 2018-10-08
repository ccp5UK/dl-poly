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

  Use kinds, Only : wp,wi
  Use comms,  Only : comms_type,gsum,gbcast,gsync,gcheck
  Use constants,  Only : pi,boltz,delth_max,npdfdt,npdgdt, &
                     engunit,zero_plus,twopi, &
    delth_max,ntable
  Use site, Only : site_type
  Use configuration, Only : configuration_type
  Use parse, Only : get_line,get_word,word_2_real
  Use errors_warnings, Only : error,warning,info
  Use numerics, Only : local_index,images
  Use particle, Only : corePart

  Implicit None

  !> Type containing angles variables
  Type, Public :: angles_type
    Private

    !> Tabulated potential switch
    Logical, Public :: l_tab = .false.

    !> Number of angle types (potentials)
    Integer( Kind = wi ), Public :: n_types  = 0
    Integer( Kind = wi ), Public :: n_types1 = 0
    !> Number of frames
    Integer( Kind = wi ), Public :: n_frames  = 0
    !> Total number of angles (all nodes)
    Integer( Kind = wi ), Public :: total

    Integer( Kind = wi ), Allocatable, Public :: num(:),key(:)

    !> Atom indices (local)
    Integer( Kind = wi ), Allocatable, Public :: lst(:,:)
    !> Atom indices
    Integer( Kind = wi ), Allocatable, Public :: list(:,:)
    !> Legend
    Integer( Kind = wi ), Allocatable, Public :: legend(:,:)

    !> Angle parameters (force constant, etc.)
    Real( Kind = wp ), Allocatable, Public :: param(:,:)

    ! Possible tabulated calculation arrays
    Integer,           Allocatable, Public :: ltp(:)
    !> Tabulated potential
    Real( Kind = wp ), Allocatable, Public :: tab_potential(:,:)
    !> Tabulated force
    Real( Kind = wp ), Allocatable, Public :: tab_force(:,:)

    ! Possible distribution arrays
    Integer,           Allocatable, Public :: ldf(:),typ(:,:)
    Real( Kind = wp ), Allocatable, Public :: dst(:,:)

    ! Maximums
    !> Maximum number of angle types
    Integer( Kind = wi ), Public :: max_types
    !> Maximum number of angles per node
    Integer( Kind = wi ), Public :: max_angles
    !> Length of legend array
    Integer( Kind = wi ), Public :: max_legend
    !> Maximum number of angle parameters
    Integer( Kind = wi ), Public :: max_param

    ! Number of bins
    !> Angular distribution function bins
    Integer( Kind = wi ), Public :: bin_adf
    !> Tabulated potential bins
    Integer( Kind = wi ), Public :: bin_tab

  Contains
    Private

    Final :: cleanup
  End Type angles_type

  Public :: allocate_angles_arrays , &
            allocate_angl_pot_arrays , allocate_angl_dst_arrays, angles_compute, &
            angles_forces, angles_table_read

Contains

  Subroutine allocate_angles_arrays(angle,mxatdm,mxtmls)
    Type( angles_type ), Intent( InOut ) :: angle
    Integer( Kind =wi ), Intent( In ) :: mxatdm,mxtmls

    Integer, Dimension(8) :: fail

    fail = 0

    Allocate (angle%num(1:mxtmls),          Stat = fail(1))
    Allocate (angle%key(1:angle%max_types),          Stat = fail(2))
    Allocate (angle%lst(1:3,1:angle%max_types),      Stat = fail(3))
    Allocate (angle%list(0:3,1:angle%max_angles),     Stat = fail(4))
    Allocate (angle%legend(0:angle%max_legend,1:mxatdm), Stat = fail(5))
    Allocate (angle%param(1:angle%max_param,1:angle%max_types), Stat = fail(6))
    If (angle%l_tab) &
    Allocate (angle%ltp(0:angle%max_types),          Stat = fail(7))
    If (angle%bin_adf > 0) &
    Allocate (angle%ldf(0:angle%max_types),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1013)

    angle%num  = 0
    angle%key  = 0
    angle%lst  = 0
    angle%list = 0
    angle%legend  = 0

    angle%param  = 0.0_wp

    If (angle%l_tab) &
    angle%ltp  = 0

    If (angle%bin_adf > 0) &
    angle%ldf  = 0

  End Subroutine allocate_angles_arrays

  Subroutine allocate_angl_pot_arrays(angle)
    Type( angles_type ), Intent( InOut ) :: angle

    Integer :: fail(2)

    fail = 0

    Allocate (angle%tab_potential(-1:angle%bin_tab,1:angle%ltp(0)), Stat = fail(1))
    Allocate (angle%tab_force(-1:angle%bin_tab,1:angle%ltp(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1074)

    angle%tab_potential = 0.0_wp
    angle%tab_force = 0.0_wp

  End Subroutine allocate_angl_pot_arrays

  Subroutine allocate_angl_dst_arrays(angle)
    Type( angles_type ), Intent( InOut ) :: angle

    Integer :: fail

    fail = 0

    Allocate (angle%typ(-1:3,1:angle%ldf(0)),angle%dst(1:angle%bin_adf,1:angle%ldf(0)), Stat = fail)

    If (fail > 0) Call error(1075)

    angle%typ = 0
    angle%dst = 0.0_wp

  End Subroutine allocate_angl_dst_arrays
  
  Subroutine angles_compute(temp,unique_atom,angle,config,comm)

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
  Character( Len = 8 ), Dimension(:), Intent( In    ) :: unique_atom
  Type( angles_type ), Intent( InOut ) :: angle
  Type( configuration_type ), Intent( InOut ) :: config
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
  Allocate (dstdang(0:angle%bin_adf,1:angle%ldf(0)),pmf(0:angle%bin_adf+2),vir(0:angle%bin_adf+2), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'angles_compute - allocation failure, node'
     Call error(0,message)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! conversion: radians <-> degrees (affects not only angle units but also force units!)

  rad2dgr = 180.0_wp/pi
  dgr2rad = pi/180.0_wp

! grid interval for pdf/pmf tables

  delth = pi/Real(angle%bin_adf,wp)
  rdlth = Real(angle%bin_adf,wp)/180.0_wp

! resampling grid and grid interval for pmf tables

  ngrid = Max(Nint(180.0_wp/delth_max),angle%bin_adf,angle%bin_tab-4)
  dgrid = pi/Real(ngrid,wp)

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,angle%ldf(0)
     If (angle%typ(0,i) > 0) Then
        kk=kk+1
        ll=ll+angle%typ(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(angle%n_frames,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-5_wp

  Call info('',.true.)
  Call info('ANGLES : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)',.true.)
  Write(message,'(a,5(1x,i10))') '# bins, cutoff, frames, types: ',angle%bin_adf,180,angle%n_frames,kk,ll
  Call info(message,.true.)

! open RDF file and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdfdt, File='ANGDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//config%cfgname
     Write(npdfdt,'(a)') '# ANGLES: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
     Write(npdfdt,'(a,4(1x,i10))') '# bins, cutoff, frames, types: ',angle%bin_adf,180,angle%n_frames,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)  PDF_norm(Theta)/Sin(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,angle%ldf(0)
     If (angle%typ(0,i) > 0) Then
        j=j+1

        Write(message,'(a,3(a8,1x),2(i10,1x))') 'type, index, instances: ', &
          unique_atom(angle%typ(1,i)),unique_atom(angle%typ(2,i)), &
          unique_atom(angle%typ(3,i)),j,angle%typ(0,i)
        Call info(message,.true.)
        Write(message,'(a,f8.5)') &
         'Theta(degrees)  PDF_ang(Theta)  Sum_PDF_ang(Theta)   @   dTheta_bin = ',delth*rad2dgr
        Call info(message,.true.)

        If (comm%idnode == 0) Then
           Write(npdfdt,'(/,a,3(a8,1x),2(i10,1x))') '# type, index, instances: ', &
                unique_atom(angle%typ(1,i)),unique_atom(angle%typ(2,i)), &
                unique_atom(angle%typ(3,i)),j,angle%typ(0,i)
        End If

! global sum of data on all nodes

        Call gsum(comm,angle%dst(1:angle%bin_adf,i))

! factor in instances (first, pdfang is normalised to unity)

        factor1=factor/Real(angle%typ(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,angle%bin_adf
           If (zero .and. ig < (angle%bin_adf-3)) zero=(angle%dst(ig+2,i) <= 0.0_wp)

           pdfang = angle%dst(ig,i)*factor1
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
           If (.not.zero) Then
             Write(message,"(f11.5,1p,2e14.6)") theta,pdfang1,sum1
             Call info(message,.true.)
           End If
           If (comm%idnode == 0) Then
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
     Write(npdgdt,'(a)') '# '//config%cfgname
     Write(npdgdt,'(a,2i6,f12.5,i10,a,e15.7)') '# ',angle%bin_adf,180,delth*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) = ', kT2engo

     Open(Unit=npdfdt, File='ANGTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//config%cfgname
     Write(npdfdt,'(a,2i6,f12.5,i10,a,e15.7)') '# ',ngrid,180,dgrid*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) = ', kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,angle%ldf(0)
     If (angle%typ(0,i) > 0) Then
        j=j+1

        If (comm%idnode == 0) Then
           Write(npdgdt,'(/,a,3(a8,1x),2(i10,1x),a)') '# ', &
             unique_atom(angle%typ(1,i)),unique_atom(angle%typ(2,i)), &
             unique_atom(angle%typ(3,i)),j,angle%typ(0,i), &
             ' (type, index, instances)'
           Write(npdfdt,'(/,a,3(a8,1x),2(i10,1x),a)') '# ', &
             unique_atom(angle%typ(1,i)),unique_atom(angle%typ(2,i)), &
             unique_atom(angle%typ(3,i)),j,angle%typ(0,i), &
             ' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,angle%bin_adf
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth

           If (dstdang(ig,i) > zero_plus) Then
              fed = -Log(dstdang(ig,i))-fed0
              If (fed0 <= zero_plus) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < angle%bin_adf-1) Then
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
           Else If (ig == angle%bin_adf) Then
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
        pmf(angle%bin_adf+1) = 2.0_wp*pmf(angle%bin_adf)  -pmf(angle%bin_adf-1)
        vir(angle%bin_adf+1) = 2.0_wp*vir(angle%bin_adf)  -vir(angle%bin_adf-1)
        pmf(angle%bin_adf+2) = 2.0_wp*pmf(angle%bin_adf+1)-pmf(angle%bin_adf)
        vir(angle%bin_adf+2) = 2.0_wp*vir(angle%bin_adf+1)-vir(angle%bin_adf)

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
     Write(message,'(a)') 'angles_compute - deallocation failure'
     Call error(0,message)
  End If

End Subroutine angles_compute

Subroutine angles_forces(isw,engang,virang,stress,angle,config,comm)

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
  Type( configuration_type ),          Intent( InOut ) :: config
  Type( angles_type ), Intent( InOut ) :: angle
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
  Allocate (lunsafe(1:angle%max_angles),lstopt(0:3,1:angle%max_angles),       Stat=fail(1))
  Allocate (xdab(1:angle%max_angles),ydab(1:angle%max_angles),zdab(1:angle%max_angles), Stat=fail(2))
  Allocate (xdbc(1:angle%max_angles),ydbc(1:angle%max_angles),zdbc(1:angle%max_angles), Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'angles_forces allocation failure'
     Call error(0,message)
  End If


! calculate atom separation vectors

  Do i=1,angle%n_types
     lunsafe(i)=.false.

! indices of angle bonded atoms

     ia=local_index(angle%list(1,i),config%nlast,config%lsi,config%lsa) ; lstopt(1,i)=ia
     ib=local_index(angle%list(2,i),config%nlast,config%lsi,config%lsa) ; lstopt(2,i)=ib
     ic=local_index(angle%list(3,i),config%nlast,config%lsi,config%lsa) ; lstopt(3,i)=ic

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0 .and. ic > 0) Then ! Tag
        If (config%lfrzn(ia)*config%lfrzn(ib)*config%lfrzn(ic) == 0) Then
           If (ia <= config%natms .or. ib <= config%natms .or. ic <= config%natms) Then
              lstopt(0,i)=1
            End If
        End If
     Else                                       ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= config%natms) .or.   &
              (ib > 0 .and. ib <= config%natms) .or.   &
              (ic > 0 .and. ic <= config%natms)) .and. &
             (ia == 0 .or. ib == 0 .or. ic == 0) ) lunsafe(i)=.true.
     End If

! define components of bond vectors

     If (lstopt(0,i) > 0) Then
        xdab(i)=config%parts(ia)%xxx-config%parts(ib)%xxx
        ydab(i)=config%parts(ia)%yyy-config%parts(ib)%yyy
        zdab(i)=config%parts(ia)%zzz-config%parts(ib)%zzz

        xdbc(i)=config%parts(ic)%xxx-config%parts(ib)%xxx
        ydbc(i)=config%parts(ic)%yyy-config%parts(ib)%yyy
        zdbc(i)=config%parts(ic)%zzz-config%parts(ib)%zzz
     Else ! (DEBUG)
        xdab(i)=0.0_wp
        ydab(i)=0.0_wp
        zdab(i)=0.0_wp

        xdbc(i)=0.0_wp
        ydbc(i)=0.0_wp
        zdbc(i)=0.0_wp
     End If
  End Do

! Check for uncompressed units

  safe = .not. Any(lunsafe(1:angle%n_types))
  Call gcheck(comm,safe)
  If (.not.safe) Then
     Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
           Do i=1,angle%n_types
             If (lunsafe(i)) Then
               Write(message,'(a,2(i10,a))')     &
                 'global unit number', angle%list(0,i), &
                 ' , with a head particle number', angle%list(1,i),   &
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

  Call images(config%imcon,config%cell,angle%n_types,xdab,ydab,zdab)
  Call images(config%imcon,config%cell,angle%n_types,xdbc,ydbc,zdbc)

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
     rdelth = Real(angle%bin_adf,wp)/pi
     angle%n_frames = angle%n_frames + 1
  End If

! loop over all specified angle potentials

  Do i=1,angle%n_types
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

        kk=angle%list(0,i)
        keya = Abs(angle%key(kk))

! accumulate the histogram (distribution)

        If (Mod(isw,2) == 0 .and. ib <= config%natms) Then
           j = angle%ldf(kk)
           l = Min(1+Int(theta*rdelth),angle%bin_adf)

           angle%dst(l,j) = angle%dst(l,j) + 1.0_wp
        End If
        If (isw == 0) Cycle

        If      (keya == 1) Then

! harmonic potential

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0

           tmp   =k*dtheta

           pterm=tmp*0.5_wp*dtheta
           gamma=tmp*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 2) Then

! quartic potential

           k2    =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           k3    =angle%param(3,kk)
           k4    =angle%param(4,kk)

           pterm=0.5_wp*k2*dtheta**2+(k3/3.0_wp)*dtheta**3+0.25*k4*dtheta**4
           gamma=dtheta*(k2+k3*dtheta+k4*dtheta**2)*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 3) Then

! truncated Harmonic potential

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           rho   =angle%param(3,kk)
           switch=-(rab**8+rbc**8)/rho**8

           tmp   =k*dtheta*Exp(switch)

           pterm=tmp*0.5_wp*dtheta
           gamma=tmp*rsint
           gamsa=pterm*8.0_wp*rab**7/rho**8
           gamsc=pterm*8.0_wp*rbc**7/rho**8
           vterm=pterm*8.0_wp*switch

        Else If (keya == 4) Then

! screened Harmonic potential

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           rho1  =angle%param(3,kk)
           rho2  =angle%param(4,kk)
           switch=-(rab/rho1+rbc/rho2)

           tmp   =k*dtheta*Exp(switch)

           pterm=tmp*0.5_wp*dtheta
           gamma=tmp*rsint
           gamsa=pterm/rho1
           gamsc=pterm/rho2
           vterm=pterm*switch

        Else If (keya == 5) Then

! screened Vessal potential (type 1)

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dth0pi=theta0-pi
           dthpi =theta -pi
           dth   =dth0pi**2-dthpi**2
           rho1  =angle%param(3,kk)
           rho2  =angle%param(4,kk)
           switch=-(rab/rho1+rbc/rho2)

           tmp   =(k*dth/(2.0_wp*dth0pi**2)) * Exp(switch)

           pterm=tmp*0.25_wp*dth
           gamma=tmp*dthpi*rsint
           gamsa=pterm/rho1
           gamsc=pterm/rho2
           vterm=pterm*switch

        Else If (keya == 6) Then

! truncated Vessal potential (type 2)

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           dth0pi=theta0-pi
           dthpi =theta -pi
           a     =angle%param(3,kk)
           rho   =angle%param(4,kk)
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

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=Cos(theta)-Cos(theta0)

           tmp   =k*dtheta

           pterm=tmp*0.5_wp*dtheta
           gamma=-tmp
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 8) Then

! ordinary cosine potential

           k    =angle%param(1,kk)
           delta=angle%param(2,kk)
           m    =angle%param(3,kk)
           a    =m*theta-delta

           pterm=k*(1.0_wp+Cos(a))
           gamma=-k*m*Sin(a)*rsint
           gamsa=0.0_wp
           gamsc=0.0_wp
           vterm=0.0_wp

        Else If (keya == 9) Then

! MM3-stretch-bend potential

           a     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           rho1  =angle%param(3,kk)
           dr1   =rab-rho1
           rho2  =angle%param(4,kk)
           dr2   =rbc-rho2

           tmp   =a*dr1*dr2

           pterm=tmp*dtheta
           gamma=tmp*rsint
           gamsa=-pterm/dr1
           gamsc=-pterm/dr2
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 10) Then

! compass stretch-stretch potential

           a     =angle%param(1,kk)
           rho1  =angle%param(2,kk)
           dr1   =rab-rho1
           rho2  =angle%param(3,kk)
           dr2   =rbc-rho2

           pterm=a*dr1*dr2
           gamma=0.0_wp
           gamsa=-a*dr2
           gamsc=-a*dr1
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 11) Then

! compass stretch-bend potential

           a     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           rho1  =angle%param(3,kk)
           dr1   =rab-rho1

           tmp   =a*dr1

           pterm=tmp*dtheta
           gamma=tmp*rsint
           gamsa=-a*dtheta
           gamsc=0.0_wp
           vterm=-gamsa*rab

        Else If (keya == 12) Then

! combined compass angle potential with 3 coupling terms

           a     =angle%param(1,kk)
           b     =angle%param(2,kk)
           c     =angle%param(3,kk)
           theta0=angle%param(4,kk)
           dtheta=theta-theta0
           rho1  =angle%param(5,kk)
           dr1   =rab-rho1
           rho2  =angle%param(6,kk)
           dr2   =rbc-rho2

           tmp   =b*dr1+c*dr2

           pterm=a*dr1*dr2 + dtheta*tmp
           gamma=tmp*rsint
           gamsa=-a*dr2-b*dtheta
           gamsc=-a*dr1-c*dtheta
           vterm=-(gamsa*rab+gamsc*rbc)

        Else If (keya == 13) Then

! MM3-angle-bend potential

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
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

           k     =angle%param(1,kk)
           theta0=angle%param(2,kk)
           dtheta=theta-theta0
           gr    =angle%param(3,kk)
           rm    =angle%param(4,kk)
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

           j = angle%ltp(kk)
           rdr = angle%tab_force(-1,j) ! 1.0_wp/delpot (in rad^-1)

           l   = Int(theta*rdr)
           ppp = theta*rdr - Real(l,wp)

           vk  = angle%tab_potential(l,j)
           vk1 = angle%tab_potential(l+1,j)
           vk2 = angle%tab_potential(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           pterm = t1 + (t2-t1)*ppp*0.5_wp

           vk  = angle%tab_force(l,j) ; If (l == 0) vk = vk*theta
           vk1 = angle%tab_force(l+1,j)
           vk2 = angle%tab_force(l+2,j)

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

        If (ia <= config%natms) Then

           config%parts(ia)%fxx=config%parts(ia)%fxx+fxa
           config%parts(ia)%fyy=config%parts(ia)%fyy+fya
           config%parts(ia)%fzz=config%parts(ia)%fzz+fza

        End If

        If (ib <= config%natms) Then

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

           config%parts(ib)%fxx=config%parts(ib)%fxx-fxa-fxc
           config%parts(ib)%fyy=config%parts(ib)%fyy-fya-fyc
           config%parts(ib)%fzz=config%parts(ib)%fzz-fza-fzc

        End If

        If (ic <= config%natms) Then

           config%parts(ic)%fxx=config%parts(ic)%fxx+fxc
           config%parts(ic)%fyy=config%parts(ic)%fyy+fyc
           config%parts(ic)%fzz=config%parts(ic)%fzz+fzc

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
     Write(message,'(a)') 'angles_forces deallocation failure'
     Call error(0,message)
  End If

End Subroutine angles_forces

Subroutine angles_table_read(angl_name,angle,sites,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABANG file (for angle potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( angles_type ), Intent( InOut ) :: angle
  Character( Len = 24 ), Intent( In    ) :: angl_name(1:angle%max_types)
  Type( site_type ), Intent( In    ) :: sites
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
  Character( Len = 256 ) :: messages(4)

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

  dlrpot = 180.0_wp/Real(angle%bin_tab-4,wp)

! check grid spacing

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > delth_max .and. (.not.safe)) Then
    Write(messages(1),'(a,1p,e15.7)') 'expected (maximum) angular increment : ', delth_max
    Write(messages(2),'(a,1p,e15.7)') 'TABANG file actual angular increment : ', delpot
    Write(messages(3),'(a,0p,i10)') ' expected (minimum) number of grid points : ', angle%bin_tab-4
    Write(messages(4),'(a,0p,i10)') ' TABANG file actual number of grid points : ', ngrid
    Call info(messages,4,.true.)
    Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     Write(message,'(a,i10)') 'TABANG arrays resized for mxgrid = ', angle%bin_tab-4
     Call info(message,.true.)
  End If

! compare grids dimensions

  If (ngrid < angle%bin_tab-4) Then
     Call warning(270,Real(ngrid,wp),Real(angle%bin_tab-4,wp),0.0_wp)
     Call error(48)
  End If

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  fail=0
  Allocate (read_type(1:angle%ltp(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'error - angles_table_read allocation failure'
     Call error(0,message)
  End If
  Call allocate_angl_pot_arrays(angle)

  read_type=0 ! initialise read_type
  Do rtang=1,angle%ltp(0)
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

     Do jtpatm=1,sites%ntype_atom
        If (atom1 == sites%unique_atom(jtpatm)) katom1=jtpatm
        If (atom2 == sites%unique_atom(jtpatm)) katom2=jtpatm
        If (atom3 == sites%unique_atom(jtpatm)) katom3=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0 .or. katom3 == 0) Then
        Write(message, '(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call info(message,.true.)
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
     Do jtang=1,angle%ltp(0)
        If (angl_name(jtang) == idangl) Then
           Do itang=1,angle%max_types
              If (angle%ltp(itang) == jtang) Exit
           End Do
           Exit
        End If
     End Do

     If (itang == 0) Then ! All(angl_name /= idangl)
        Write(message, '(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call info(message,.true.)
        Call error(83)
     End If
     If (Any(read_type == jtang)) Then
        Write(message, '(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call info(message,.true.)
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
              Write(messages(1),'(a,1p,e15.7)') 'TABANG stated angular increment : ', delpot
              Write(messages(2),'(a,1p,e15.7)') 'TABANG read-in angular increment : ', rrr
              Call info(messages,2,.true.)
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0
           rrr0      = rrr

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              Write(messages(1),'(a,1p,e15.7)') 'TABANG stated angular increment : ', delpot
              Write(messages(2),'(a,1p,e15.7)') 'TABANG read-in angular increment : ', rrr-rrr0
              Call info(messages,2,.true.)
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
              Write(messages(1),'(a,1p,e15.7)') 'TABANG stated angular increment : ', delpot
              Write(messages(2),'(a,1p,e15.7)') 'TABANG read-in angular increment : ', rrr
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
        Do i=1,angle%bin_tab-4
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
           angle%tab_potential(i,jtang) = t1 + (t2-t1)*ppp*0.5_wp
           angle%tab_potential(i,jtang) = angle%tab_potential(i,jtang)*engunit ! convert to internal units

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
           angle%tab_force(i,jtang) = t1 + (t2-t1)*ppp*0.5_wp
           angle%tab_force(i,jtang) = angle%tab_force(i,jtang)*engunit*rad2dgr ! convert to internal units
        End Do

        angle%tab_force(-1,jtang) = rad2dgr/dlrpot
     Else
        Do i=1,angle%bin_tab-4
           angle%tab_potential(i,jtang) = bufpot(i)*engunit         ! convert to internal units
           angle%tab_force(i,jtang) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do

! linear extrapolation for the grid point just beyond the cutoff

        angle%tab_potential(angle%bin_tab-3,jtang) = &
          2.0_wp*angle%tab_potential(angle%bin_tab-4,jtang) - angle%tab_potential(angle%bin_tab-5,jtang)
        angle%tab_force(angle%bin_tab-3,jtang) = &
          2.0_wp*angle%tab_force(angle%bin_tab-4,jtang) - angle%tab_force(angle%bin_tab-5,jtang)

        angle%tab_force(-1,jtang) = rad2dgr/delpot
     End If

! grid point at 0 and linear extrapolation for the grid point at angle%bin_tab-2

     angle%tab_potential(0,jtang) = bufpot(0)
     angle%tab_force(0,jtang) = bufvir(0)

     angle%tab_potential(angle%bin_tab-2,jtang) = &
       2.0_wp*angle%tab_potential(angle%bin_tab-3,jtang) - angle%tab_potential(angle%bin_tab-4,jtang)
     angle%tab_force(angle%bin_tab-2,jtang) = &
       2.0_wp*angle%tab_force(angle%bin_tab-3,jtang) - angle%tab_force(angle%bin_tab-4,jtang)
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=ntable)
  End If
  Call info('',.true.)
  Call info('potential tables read from TABANG file',.true.)

! Break if not safe

  Call gcheck(comm,safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'error - angles_table_read deallocation failure'
     Call error(0,message)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine angles_table_read

  Subroutine cleanup(angle)
    Type(angles_type) :: angle

    If (Allocated(angle%num)) Then
      Deallocate(angle%num)
    End If
    If (Allocated(angle%key)) Then
      Deallocate(angle%key)
    End If

    If (Allocated(angle%lst)) Then
      Deallocate(angle%lst)
    End If
    If (Allocated(angle%list)) Then
      Deallocate(angle%list)
    End If
    If (Allocated(angle%legend)) Then
      Deallocate(angle%legend)
    End If

    If (Allocated(angle%param)) Then
      Deallocate(angle%param)
    End If

    If (Allocated(angle%ltp)) Then
      Deallocate(angle%ltp)
    End If
    If (Allocated(angle%tab_potential)) Then
      Deallocate(angle%tab_potential)
    End If
    If (Allocated(angle%tab_force)) Then
      Deallocate(angle%tab_force)
    End If
  End Subroutine cleanup
End Module angles
