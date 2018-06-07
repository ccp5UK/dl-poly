Module dihedrals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global dihedral interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
! contrib   - a.v.brukhno march 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp,wi
  Use comms, Only : comms_type,gcheck,gsum,gsync,gbcast
  Use setup,  Only : pi,twopi,boltz,delth_max,npdfdt,npdgdt, &
                            engunit,zero_plus, mxtmls,     &
                            rtwopi,r4pie0,    &
                            mximpl, ntable,mxgvdw,mxatdm
  Use site,   Only : unqatm,ntpatm
  Use configuration, Only : imcon,cell,natms,nlast,lsi,lsa,ltg,lfrzn,ltype, &
                                chge,xxx,yyy,zzz,fxx,fyy,fzz,cfgname
  Use vdw,    Only : ntpvdw,gvdw,vvdw,afs,prmvdw,bfs,ls_vdw,ld_vdw,  &
                               lstvdw,ltpvdw
  Use parse,  Only : get_line,get_word,word_2_real
  Use errors_warnings, Only : error,warning,info
  Use numerics, Only : local_index,images
  Use coul_spole, Only : intra_coul
  Use coul_mpole, Only : intra_mcoul
  Use angles, Only : angles_type
  Implicit None

  Private

  Type, Public :: dihedrals_type
    Private

    !> Tabulated potential switch
    Logical, Public :: l_tab = .false.
    !> Core shell switch
    Logical, Public :: l_core_shell = .false.

    !> Number of dihedral angle types (potentials)
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
    !> Maximum number of dihedral angle types
    Integer( Kind = wi ), Public :: max_types
    !> Maximum number of dihedral angles per node
    Integer( Kind = wi ), Public :: max_angles
    !> Length of legend array
    Integer( Kind = wi ), Public :: max_legend
    !> Maximum number of dihedral parameters
    Integer( Kind = wi ), Public :: max_param

    ! Number of bins
    !> Angular distribution function bins
    Integer( Kind = wi ), Public :: bin_adf
    !> Tabulated potential bins
    Integer( Kind = wi ), Public :: bin_tab
  Contains
    Private

    Final :: cleanup
  End Type dihedrals_type


  Public :: allocate_dihedrals_arrays, allocate_dihd_pot_arrays, &
            allocate_dihd_dst_arrays, dihedrals_compute, dihedrals_table_read, &
            dihedrals_14_check, dihedrals_forces

Contains

  Subroutine allocate_dihedrals_arrays(dihedral)
    Type( dihedrals_type ), Intent( InOut ) :: dihedral

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (dihedral%num(1:mxtmls),          Stat = fail(1))
    Allocate (dihedral%key(1:dihedral%max_types),          Stat = fail(2))
    Allocate (dihedral%lst(1:6,1:dihedral%max_types),      Stat = fail(3))
    Allocate (dihedral%list(0:6,1:dihedral%max_angles),     Stat = fail(4))
    Allocate (dihedral%legend(0:dihedral%max_legend,1:mxatdm), Stat = fail(5))
    Allocate (dihedral%param(1:dihedral%max_param,1:dihedral%max_types), Stat = fail(6))
    If (dihedral%l_tab) &
    Allocate (dihedral%ltp(0:dihedral%max_types),          Stat = fail(7))
    If (dihedral%bin_adf > 0) &
    Allocate (dihedral%ldf(0:dihedral%max_types),          Stat = fail(8))


    If (Any(fail > 0)) Call error(1020)

    dihedral%num  = 0
    dihedral%key  = 0
    dihedral%lst  = 0
    dihedral%list = 0
    dihedral%legend  = 0

    dihedral%param  = 0.0_wp

    If (dihedral%l_tab) &
    dihedral%ltp  = 0

    If (dihedral%bin_adf > 0) &
    dihedral%ldf  = 0

  End Subroutine allocate_dihedrals_arrays

  Subroutine allocate_dihd_pot_arrays(dihedral)
    Type( dihedrals_type ), Intent( InOut ) :: dihedral

    Integer :: fail(1:2)

    fail = 0

    Allocate (dihedral%tab_potential(-1:dihedral%bin_tab,1:dihedral%ltp(0)), Stat = fail(1))
    Allocate (dihedral%tab_force(-1:dihedral%bin_tab,1:dihedral%ltp(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1076)

    dihedral%tab_potential = 0.0_wp
    dihedral%tab_force = 0.0_wp

  End Subroutine allocate_dihd_pot_arrays

  Subroutine allocate_dihd_dst_arrays(dihedral)
    Type( dihedrals_type ), Intent( InOut ) :: dihedral

    Integer :: fail

    fail = 0

    Allocate (dihedral%typ(-1:4,1:dihedral%ldf(0)),dihedral%dst(1:dihedral%bin_adf,1:dihedral%ldf(0)), Stat = fail)

    If (fail > 0) Call error(1077)

    dihedral%typ = 0
    dihedral%dst = 0.0_wp

  End Subroutine allocate_dihd_dst_arrays
  
  Subroutine dihedrals_14_check &
           (l_str,l_top,ntpmls,nummols,angle,dihedral,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to check and resolve conflicting molecular
! forcefield specification for 1-4 interactions
!
! copyright - daresbury laboratory
! author    - w.smith march 1999
! amended   - i.t.todorov september 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Logical,           Intent( In    ) :: l_str,l_top
  Integer,           Intent( In    ) :: ntpmls,nummols(1:mxtmls)
  Type( angles_type ), Intent( In    ) :: angle
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( comms_type), Intent( InOut ) :: comm

  Logical :: l_print,l_reset,l_reset_l
  Integer :: kangle,kdihed,itmols,imols,langle,ldihed,mdihed, &
             iang,jang,idih,jdih,kdih,ldih,mdih,ndih,odih,pdih

  l_print = (l_str .and. l_top)
  l_reset = .false.

! Initialise angle and dihedral interaction counters

  kangle=0
  kdihed=0

! loop over molecular types

  Do itmols=1,ntpmls

! loop over molecules in system

     Do imols=1,nummols(itmols)

! check for valence angle on dihedral angle conflicts

        Do langle=1,angle%num(itmols)

           If (angle%key(langle+kangle) > 0) Then

              iang=angle%lst(1,langle+kangle)
              jang=angle%lst(3,langle+kangle)

              Do ldihed=1,dihedral%num(itmols)

                 idih=dihedral%lst(1,ldihed+kdihed)
                 jdih=dihedral%lst(4,ldihed+kdihed)

                 If (Min(iang,jang) == Min(idih,jdih) .and. Max(iang,jang) == Max(idih,jdih)) Then
                    dihedral%param(4,ldihed+kdihed)=0.0_wp
                    dihedral%param(5,ldihed+kdihed)=0.0_wp

                    l_reset = .true.
                    If (l_print) Call warning(20,Real(itmols,wp),Real(idih,wp),Real(jdih,wp))
                 End If

              End Do

           End If

        End Do

! check for double dihedral angle conflicts

        Do ldihed=1,dihedral%num(itmols)-1

           idih=dihedral%lst(1,ldihed+kdihed)
           jdih=dihedral%lst(4,ldihed+kdihed)
           If (dihedral%l_core_shell) Then
              mdih=dihedral%lst(5,ldihed+kdihed)
              ndih=dihedral%lst(6,ldihed+kdihed)
           End If

           Do mdihed=ldihed+1,dihedral%num(itmols)

              kdih=dihedral%lst(1,mdihed+kdihed)
              ldih=dihedral%lst(4,mdihed+kdihed)
              If (dihedral%l_core_shell) Then
                 odih=dihedral%lst(5,mdihed+kdihed)
                 pdih=dihedral%lst(6,mdihed+kdihed)
              End If

              l_reset_l=.false.
              If (dihedral%l_core_shell) Then
                 If (Min(kdih,ldih,odih,pdih) == Min(idih,jdih,mdih,ndih) .and. &
                     Max(kdih,ldih,odih,pdih) == Max(idih,jdih,mdih,ndih)) Then
                    If (dihedral%param(4,ldihed+kdihed)*dihedral%param(4,mdihed+kdihed) > 1.0e-10_wp) Then
                       dihedral%param(4,mdihed+kdihed) = 0.0_wp
                       l_reset_l = .true.
                    End If

                    If (dihedral%param(5,ldihed+kdihed)*dihedral%param(5,mdihed+kdihed) > 1.0e-10_wp) Then
                       dihedral%param(5,mdihed+kdihed) = 0.0_wp
                       l_reset_l = .true.
                    End If

                    If (l_reset_l .and. l_print) Call warning(20,Real(itmols,wp),Real(kdih,wp),Real(ldih,wp))
                 End If
              Else
                 If (Min(kdih,ldih) == Min(idih,jdih) .and. Max(kdih,ldih) == Max(idih,jdih)) Then
                    If (dihedral%param(4,ldihed+kdihed)*dihedral%param(4,mdihed+kdihed) > 1.0e-10_wp) Then
                       dihedral%param(4,mdihed+kdihed) = 0.0_wp
                       l_reset_l = .true.
                    End If

                    If (dihedral%param(5,ldihed+kdihed)*dihedral%param(5,mdihed+kdihed) > 1.0e-10_wp) Then
                       dihedral%param(5,mdihed+kdihed) = 0.0_wp
                       l_reset_l = .true.
                    End If

                    If (l_reset_l .and. l_print) Call warning(20,Real(itmols,wp),Real(kdih,wp),Real(ldih,wp))
                 End If
              End If
              l_reset = (l_reset .or. l_reset_l)

           End Do

        End Do

     End Do

! Update counters

     kangle=kangle+angle%num(itmols)
     kdihed=kdihed+dihedral%num(itmols)

  End Do

  If (.not.l_print) Then
     Call gcheck(comm,l_reset,"enforce")
     If (l_reset) Call warning(22,0.0_wp,0.0_wp,0.0_wp)
  End If

End Subroutine dihedrals_14_check

Subroutine dihedrals_compute(temp,dihedral,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating dihedrals distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ), Intent( In    ) :: temp
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( comms_type ), Intent( InOut ) :: comm

  Logical           :: zero
  Integer           :: fail,ngrid,i,j,ig,kk,ll
  Real( Kind = wp ) :: kT2engo,rad2dgr,dgr2rad,delth,rdlth,dgrid,factor,pdfzero, &
                       factor1,theta,pdfdih,sum,pdfdih1,sum1,                    &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstddih(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  Character( Len = 256 ) :: message,messages(3)

  fail = 0
  Allocate (dstddih(0:dihedral%bin_adf,1:dihedral%ldf(0)),pmf(0:dihedral%bin_adf+2),vir(0:dihedral%bin_adf+2), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'dihedrals_compute - allocation failure'
     Call error(0,message)
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! conversion: radians <-> degrees (affects not only dihedral units but also force units!)

  rad2dgr = 180.0_wp/pi
  dgr2rad = pi/180.0_wp

! grid interval for pdf/pmf tables

  delth = twopi/Real(dihedral%bin_adf,wp)
  rdlth = Real(dihedral%bin_adf,wp)/360.0_wp

! resampling grid and grid interval for pmf tables

  ngrid = Max(Nint(360.0_wp/delth_max),dihedral%bin_adf,dihedral%bin_tab-4)
  dgrid = twopi/Real(ngrid,wp)

! loop over all valid PDFs to get valid totals

  kk=0
  ll=0
  Do i=1,dihedral%ldf(0)
     If (dihedral%typ(0,i) > 0) Then
        kk=kk+1
        ll=ll+dihedral%typ(0,i)
     End If
  End Do

! normalisation factor

  factor = 1.0_wp/Real(dihedral%n_frames,wp)

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-5_wp

  Write(messages(1),*) ''
  Write(messages(2),'(a)') &
     'DIHEDRALS : Probability Distribution Functions (PDF) := histogram(bin)/hist_sum(bins)'
  Write(messages(3),'(a,i10,1x,a,2(i0,a),3(1x,i10))') &
     '# bins, range, frames, types: ',dihedral%bin_adf,'[',-180,',',180,']',dihedral%n_frames,kk,ll
  Call info(messages,3,.true.)

! open RDF file and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdfdt, File='DIHDAT', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a)') '# DIHEDRALS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin'
     Write(npdfdt,'(a,4(1x,i10))') '# bins, cutoff, frames, types: ',dihedral%bin_adf,360,dihedral%n_frames,kk
     Write(npdfdt,'(a)') '#'
     Write(npdfdt,'(a,f8.5)') '# Theta(degrees)  PDF_norm(Theta)   @   dTheta_bin = ',delth*rad2dgr
     Write(npdfdt,'(a)') '#'
  End If

! loop over all valid PDFs

  j=0
  Do i=1,dihedral%ldf(0)
     If (dihedral%typ(0,i) > 0) Then
        j=j+1

        Write(messages(1),*) ''
        Write(messages(2),'(a,4(a8,1x),2(i10,1x))') 'type, index, instances: ', &
           unqatm(dihedral%typ(1,i)),unqatm(dihedral%typ(2,i)),unqatm(dihedral%typ(3,i)), &
           unqatm(dihedral%typ(4,i)),j,dihedral%typ(0,i)
        Write(messages(3),'(a,f8.5)') &
           'Theta(degrees)  P_dih(Theta)  Sum_P_dih(Theta)   @   dTheta_bin = ',delth*rad2dgr
        Call info(messages,3,.true.)
        If (comm%idnode == 0) Then
           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x))') '# type, index, instances: ', &
              unqatm(dihedral%typ(1,i)),unqatm(dihedral%typ(2,i)),unqatm(dihedral%typ(3,i)), &
              unqatm(dihedral%typ(4,i)),j,dihedral%typ(0,i)
        End If

! global sum of data on all nodes

        Call gsum(comm,dihedral%dst(1:dihedral%bin_adf,i))

! factor in instances (first, pdfdih is normalised to unity)

        factor1=factor/Real(dihedral%typ(0,i),wp)

! running integration of pdf

        sum=0.0_wp

! loop over distances

        zero=.true.
        Do ig=1,dihedral%bin_adf
           If (zero .and. ig < (dihedral%bin_adf-3)) zero=(dihedral%dst(ig+2,i) <= 0.0_wp)

           pdfdih = dihedral%dst(ig,i)*factor1
           sum = sum + pdfdih

! null it if < pdfzero

           If (pdfdih < pdfzero) Then
              pdfdih1 = 0.0_wp
           Else
              pdfdih1 = pdfdih
           End If

           If (sum < pdfzero) Then
              sum1 = 0.0_wp
           Else
              sum1 = sum
           End If

           theta = (Real(ig,wp)-0.5_wp)*delth-pi

! now pdfdih is normalised by the volume element (as to go to unity at infinity in gases and liquids)

           pdfdih = pdfdih*rdlth

! print out information

           theta  = theta*rad2dgr
           If (.not.zero) Then
             Write(message,'(f11.5,1p,2e14.6)') theta,pdfdih1,sum1
             Call info(message,.true.)
           End If
           If (comm%idnode == 0) Then
              Write(npdfdt,"(f11.5,1p,e14.6)") theta,pdfdih
           End If

! We use the non-normalised tail-truncated PDF version,
! pdf...1 (not pdf...) in order to exclude the nearly-zero
! pdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

           dstddih(ig,i) = pdfdih1 ! PDFs density
        End Do
     Else
        dstddih(:,i) = 0.0_wp ! PDFs density
     End If
  End Do

  If (comm%idnode == 0) Close(Unit=npdfdt)

! open PDF files and write headers

  If (comm%idnode == 0) Then
     Open(Unit=npdgdt, File='DIHPMF', Status='replace')
     Write(npdgdt,'(a)') '# '//cfgname
     Write(npdgdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',dihedral%bin_adf,delth*Real(dihedral%bin_adf,wp)*rad2dgr,delth*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo

     Open(Unit=npdfdt, File='DIHTAB', Status='replace')
     Write(npdfdt,'(a)') '# '//cfgname
     Write(npdfdt,'(a,i10,2f12.5,i10,a,e15.7)') '# ',ngrid,dgrid*Real(ngrid,wp)*rad2dgr,dgrid*rad2dgr,kk, &
          '   conversion factor(kT -> energy units) =',kT2engo
  End If

! loop over all valid PDFs

  j=0
  Do i=1,dihedral%ldf(0)
     If (dihedral%typ(0,i) > 0) Then
        j=j+1

        If (comm%idnode == 0) Then
           Write(npdgdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(dihedral%typ(1,i)),unqatm(dihedral%typ(2,i)),unqatm(dihedral%typ(3,i)), &
                unqatm(dihedral%typ(4,i)),j,dihedral%typ(0,i),' (type, index, instances)'
           Write(npdfdt,'(/,a,4(a8,1x),2(i10,1x),a)') '# ', &
                unqatm(dihedral%typ(1,i)),unqatm(dihedral%typ(2,i)),unqatm(dihedral%typ(3,i)), &
                unqatm(dihedral%typ(4,i)),j,dihedral%typ(0,i),' (type, index, instances)'
        End If

! Smoothen and get derivatives

        fed0  = 0.0_wp
        dfed0 = 10.0_wp
        dfed  = 10.0_wp

        Do ig=1,dihedral%bin_adf
           tmp = Real(ig,wp)-0.5_wp
           theta = tmp*delth-pi

           If (dstddih(ig,i) > zero_plus) Then
              fed = -Log(dstddih(ig,i))-fed0
              If (fed0 <= zero_plus) Then
                 fed0 = fed
                 fed  = 0.0_wp
              End If

              If (ig < dihedral%bin_adf-1) Then
                 If (dstddih(ig+1,i) <= zero_plus .and. dstddih(ig+2,i) > zero_plus) &
                    dstddih(ig+1,i) = 0.5_wp*(dstddih(ig,i)+dstddih(ig+2,i))
              End If
           Else
              fed = 0.0_wp
           End If

           If      (ig == 1) Then
              If      (dstddih(ig,i) > zero_plus .and. dstddih(ig+1,i) > zero_plus) Then
                 dfed = Log(dstddih(ig+1,i)/dstddih(ig,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (ig == dihedral%bin_adf) Then
              If      (dstddih(ig,i) > zero_plus .and. dstddih(ig-1,i) > zero_plus) Then
                 dfed = Log(dstddih(ig,i)/dstddih(ig-1,i))
              Else If (dfed > 0.0_wp) Then
                 dfed = dfed0
              Else
                 dfed =-dfed0
              End If
           Else If (dstddih(ig-1,i) > zero_plus) Then
              If (dstddih(ig+1,i) > zero_plus) Then
                 dfed = 0.5_wp*(Log(dstddih(ig+1,i)/dstddih(ig-1,i)))
              Else
                 dfed = 0.5_wp*Log(dstddih(ig-1,i))
              End If
           Else If (dstddih(ig+1,i) > zero_plus) Then
              dfed =-0.5_wp*Log(dstddih(ig+1,i))
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

! Cyclic grid

        pmf(0)         = 0.5_wp*(pmf(1)-pmf(dihedral%bin_adf))
        vir(0)         = 0.5_wp*(vir(1)-vir(dihedral%bin_adf))
        pmf(dihedral%bin_adf+1) = pmf(0)
        vir(dihedral%bin_adf+1) = vir(0)
        pmf(dihedral%bin_adf+2) = pmf(1)
        vir(dihedral%bin_adf+2) = vir(1)

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
              Write(npdfdt,"(f11.5,1p,2e14.6)") (theta-pi)*rad2dgr,fed*kT2engo,dfed*kT2engo*dgr2rad/delth
        End Do
     End If
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=npdgdt)
     Close(Unit=npdfdt)
  End If

  Deallocate (dstddih,pmf,vir, Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'dihedrals_compute - deallocation failure'
     Call error(0,message)
  End If

End Subroutine dihedrals_compute

Subroutine dihedrals_forces &
           (isw,engdih,virdih,stress,rcut,rvdw,keyfce,alpha,epsq,engcpe,vircpe,engsrp,virsrp,dihedral,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating dihedral energy and force terms
!
! isw = 0 - collect statistics
! isw = 1 - calculate forces
! isw = 2 - do both
!
! Note: scale factors for reduces electrostatic and vdw 1-4 interactions
!       assumes 1-4 interactions are in the exclude list
!
! copyright - daresbury laboratory
! author    - w.smith march 1992
! amended   - i.t.todorov march 2016
! contrib   - a.v.brukhno & i.t.todorov april 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,                             Intent( In    ) :: isw
  Real( Kind = wp ),                   Intent(   Out ) :: engdih,virdih
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
  Real( Kind = wp ),                   Intent( In    ) :: rcut,rvdw, &
                                                          alpha,epsq
  Integer,                             Intent( In    ) :: keyfce
  Real( Kind = wp ),                   Intent( InOut ) :: engcpe,vircpe, &
                                                          engsrp,virsrp
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( comms_type),                   Intent( InOut ) :: comm

  Logical                 :: safe(1:3),csa,csd
  Integer                 :: fail(1:5),i,j,l,ia,ib,ic,id,kk,keyd, &
                             ai,aj,ia0,id0
  Real( Kind = wp )       :: xab,yab,zab, xac,yac,zac,                 &
                             xad,yad,zad,rad(0:3),rad2(0:3),           &
                             xbc,ybc,zbc,rrbc,                         &
                             xcd,ycd,zcd,                              &
                             fax,fay,faz, fb1x,fb1y,fb1z,              &
                             fcx,fcy,fcz, fd1x,fd1y,fd1z,              &
                             fx,fy,fz,                                 &
                             rdelth,rdr,ppp,vk,vk1,vk2,t1,t2,          &
                             pbx,pby,pbz,pb2,rpb1,rpb2,                &
                             pcx,pcy,pcz,pc2,rpc1,rpc2,                &
                             pbpc,cost,sint,rsint,theta,theta0,dtheta, &
                             a,a0,a1,a2,a3,d,m,term,pterm,gamma,       &
                             scale,chgprd,coul,fcoul,eng,virele,       &
                             engc14,virc14,engs14,virs14,              &
                             strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:5)

  Logical,           Allocatable :: lunsafe(:),lad(:,:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Real( Kind = wp ), Allocatable :: xdbc(:),ydbc(:),zdbc(:)
  Real( Kind = wp ), Allocatable :: xdcd(:),ydcd(:),zdcd(:)
  Real( Kind = wp ), Allocatable :: xdad(:,:),ydad(:,:),zdad(:,:)
  Character( Len = 256 ) :: message,messages(7)

  fail=0
  Allocate (lunsafe(1:dihedral%max_angles),lstopt(0:6,1:dihedral%max_angles),lad(1:3,1:dihedral%max_angles), Stat=fail(1))
  Allocate (xdab(1:dihedral%max_angles),ydab(1:dihedral%max_angles),zdab(1:dihedral%max_angles),             Stat=fail(2))
  Allocate (xdbc(1:dihedral%max_angles),ydbc(1:dihedral%max_angles),zdbc(1:dihedral%max_angles),             Stat=fail(3))
  Allocate (xdcd(1:dihedral%max_angles),ydcd(1:dihedral%max_angles),zdcd(1:dihedral%max_angles),             Stat=fail(4))
  Allocate (xdad(1:3,1:dihedral%max_angles),ydad(1:3,1:dihedral%max_angles),zdad(1:3,1:dihedral%max_angles), Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'dihedrals_forces allocation failure'
     Call error(0,message)
  End If


! calculate atom separation vectors

  Do i=1,dihedral%n_types
     lunsafe(i)=.false.

! indices of dihedral atoms

     ia=local_index(dihedral%list(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(dihedral%list(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib
     ic=local_index(dihedral%list(3,i),nlast,lsi,lsa) ; lstopt(3,i)=ic
     id=local_index(dihedral%list(4,i),nlast,lsi,lsa) ; lstopt(4,i)=id
     If (dihedral%l_core_shell) Then
        ia0=local_index(dihedral%list(5,i),nlast,lsi,lsa) ; lstopt(5,i)=ia0
        id0=local_index(dihedral%list(6,i),nlast,lsi,lsa) ; lstopt(6,i)=id0
     End If

     lstopt(0,i)=0
     If (dihedral%l_core_shell) Then
        If (ia > 0 .and. ib > 0 .and. ic > 0 .and. id > 0 .and. &
            ia0 > 0 .and. id0 > 0) Then !Tag
           If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic)*lfrzn(id) == 0) Then
              If (ia <= natms .or. ib <= natms .or. ic <= natms .or. id <= natms .or. &
                  ia0 <= natms .or. id0 <= natms) Then
                 lstopt(0,i)=1
              End If
           End If
        Else                            ! Detect uncompressed unit
           If ( ((ia > 0  .and. ia <= natms)  .or.   &
                 (ib > 0  .and. ib <= natms)  .or.   &
                 (ic > 0  .and. ic <= natms)  .or.   &
                 (id > 0  .and. id <= natms)  .or.   &
                 (ia0 > 0 .and. ia0 <= natms) .or.   &
                 (id0 > 0 .and. id0 <= natms)) .and. &
                (ia == 0 .or. ib == 0 .or. ic == 0 .or. id == 0 .or. &
                 ia0 == 0 .or. id0 == 0) ) lunsafe(i)=.true.
        End If
     Else
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
     End If

! define components of bond vectors

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)

        xdbc(i)=xxx(ib)-xxx(ic)
        ydbc(i)=yyy(ib)-yyy(ic)
        zdbc(i)=zzz(ib)-zzz(ic)

        xdcd(i)=xxx(ic)-xxx(id)
        ydcd(i)=yyy(ic)-yyy(id)
        zdcd(i)=zzz(ic)-zzz(id)

        If (dihedral%l_core_shell) Then
           csa=(dihedral%list(1,i) /= dihedral%list(5,i))
           csd=(dihedral%list(4,i) /= dihedral%list(6,i))

           lad(:,i)=.false.
           If (csa .or. csd) Then
              If (csa .and. csd) Then
                 lad(1,i)=.true.
                 xdad(1,i)=xxx(ia0)-xxx(id)
                 ydad(1,i)=yyy(ia0)-yyy(id)
                 zdad(1,i)=zzz(ia0)-zzz(id)

                 lad(2,i)=.true.
                 xdad(2,i)=xxx(ia)-xxx(id0)
                 ydad(2,i)=yyy(ia)-yyy(id0)
                 zdad(2,i)=zzz(ia)-zzz(id0)

                 lad(3,i)=.true.
                 xdad(3,i)=xxx(ia0)-xxx(id0)
                 ydad(3,i)=yyy(ia0)-yyy(id0)
                 zdad(3,i)=zzz(ia0)-zzz(id0)
              Else If (csa) Then
                 lad(1,i)=.true.
                 xdad(1,i)=xxx(ia0)-xxx(id)
                 ydad(1,i)=yyy(ia0)-yyy(id)
                 zdad(1,i)=zzz(ia0)-zzz(id)
              Else If (csd) Then
                 lad(2,i)=.true.
                 xdad(2,i)=xxx(ia)-xxx(id0)
                 ydad(2,i)=yyy(ia)-yyy(id0)
                 zdad(2,i)=zzz(ia)-zzz(id0)
              End If
           End If
        End If
     Else ! (DEBUG)
        xdab(i)=0.0_wp
        ydab(i)=0.0_wp
        zdab(i)=0.0_wp

        xdbc(i)=0.0_wp
        ydbc(i)=0.0_wp
        zdbc(i)=0.0_wp

        xdcd(i)=0.0_wp
        ydcd(i)=0.0_wp
        zdcd(i)=0.0_wp

        If (dihedral%l_core_shell) Then
           lad(i,:)=.false.
           xdad(i,:)=0.0_wp
           ydad(i,:)=0.0_wp
           zdad(i,:)=0.0_wp
        End If
     End If
  End Do

! Check for uncompressed units

  safe(1) = .not. Any(lunsafe(1:dihedral%n_types))
  Call gcheck(comm,safe(1))
  If (.not.safe(1)) Then
     Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
           Do i=1,dihedral%n_types
             If (lunsafe(i)) Then
               Write(message,'(2(a,i10))') &
                 'global unit number', dihedral%list(0,i), &
                 ' , with a head particle number', dihedral%list(1,i)
               Call info(message)
               Call warning('contributes towards next error')
             End If
           End Do
        End If
        Call gsync(comm)
     End Do
     Call error(132)
  End If

! periodic boundary condition

  Call images(imcon,cell,dihedral%n_types,xdab,ydab,zdab)
  Call images(imcon,cell,dihedral%n_types,xdbc,ydbc,zdbc)
  Call images(imcon,cell,dihedral%n_types,xdcd,ydcd,zdcd)

  If (Mod(isw,3) > 0) Then

    If (dihedral%l_core_shell) Then
      If (Any(lad(1,1:dihedral%n_types))) Then
        Call images(imcon,cell,dihedral%n_types,xdad(1,1:dihedral%n_types), &
          ydad(1,1:dihedral%n_types),zdad(1,1:dihedral%n_types))
      End If
      If (Any(lad(2,1:dihedral%n_types))) Then
        Call images(imcon,cell,dihedral%n_types,xdad(2,1:dihedral%n_types), &
          ydad(2,1:dihedral%n_types),zdad(2,1:dihedral%n_types))
      End If
      If (Any(lad(3,1:dihedral%n_types))) Then
        Call images(imcon,cell,dihedral%n_types,xdad(3,1:dihedral%n_types), &
          ydad(3,1:dihedral%n_types),zdad(3,1:dihedral%n_types))
      End If
    End If

! Initialise safety flags

     safe=.true.

! zero dihedral energy accumulator

     engdih=0.0_wp
     virdih=0.0_wp

! zero scaled 1-4 electrostatic and short-range potential accumulators

     engc14=0.0_wp
     virc14=0.0_wp
     engs14=0.0_wp
     virs14=0.0_wp

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
     rdelth = Real(dihedral%bin_adf,wp)*rtwopi
     dihedral%n_frames = dihedral%n_frames + 1
  End If

! loop over all specified dihedrals

  Do i=1,dihedral%n_types
     If (lstopt(0,i) > 0) Then

! indices of dihedral atoms

        ia=lstopt(1,i)
        ib=lstopt(2,i)
        ic=lstopt(3,i)
        id=lstopt(4,i)

! indices of 1-4 shelled dihedral atoms

        If (dihedral%l_core_shell) Then
           ia0=lstopt(5,i)
           id0=lstopt(6,i)
        End If

! define components of bond vectors

        xab=xdab(i)
        yab=ydab(i)
        zab=zdab(i)

        xbc=xdbc(i)
        ybc=ydbc(i)
        zbc=zdbc(i)
        rrbc=1.0_wp/Sqrt(xbc*xbc+ybc*ybc+zbc*zbc)

        xcd=xdcd(i)
        ycd=ydcd(i)
        zcd=zdcd(i)

        xac=xab+xbc
        yac=yab+ybc
        zac=zab+zbc

! construct first dihedral vector

        pbx=yab*zbc-zab*ybc
        pby=zab*xbc-xab*zbc
        pbz=xab*ybc-yab*xbc

        pb2=pbx*pbx+pby*pby+pbz*pbz

        rpb1=1.0_wp/Sqrt(pb2)
        rpb2=rpb1*rpb1

! construct second dihedral vector

        pcx=ybc*zcd-zbc*ycd
        pcy=zbc*xcd-xbc*zcd
        pcz=xbc*ycd-ybc*xcd

        pc2=pcx*pcx+pcy*pcy+pcz*pcz

        rpc1=1.0_wp/Sqrt(pc2)
        rpc2=rpc1*rpc1

! determine dihedral angle

        pbpc=pbx*pcx+pby*pcy+pbz*pcz
        cost=pbpc*rpb1*rpc1
        If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)
        sint=( xbc*(pcy*pbz-pcz*pby) + ybc*(pbx*pcz-pbz*pcx) + &
               zbc*(pcx*pby-pcy*pbx) ) * (rpb1*rpc1*rrbc)

        theta=Atan2(sint,cost)

! avoid singularity in sint

        sint=Sign( Max(1.0e-10_wp,Abs(sint)) , sint )
        rsint=1.0_wp/sint

! selection of potential energy function type

        kk=dihedral%list(0,i)
        keyd = Abs(dihedral%key(kk))

! accumulate the histogram (distribution)

        If (Mod(isw,2) == 0 .and. ia <= natms) Then
           j = dihedral%ldf(kk)
           l = Min(1+Int((theta+pi)*rdelth),dihedral%bin_adf)

           dihedral%dst(l,j) = dihedral%dst(l,j) + 1.0_wp
        End If
        If (isw == 0) Cycle

! calculate potential energy and scalar force term

        If      (keyd == 1) Then

! torsion dihedral potential

           a=dihedral%param(1,kk)
           d=dihedral%param(2,kk)
           m=dihedral%param(3,kk)

           term=m*theta-d

           pterm=a*(1.0_wp+Cos(term))
           gamma=-a*m*Sin(term)*rsint * rpb1*rpc1

        Else If (keyd == 2) Then

! harmonic improper dihedral

           a     =dihedral%param(1,kk)
           theta0=dihedral%param(2,kk)
           dtheta=theta-theta0
           dtheta=dtheta-Real(Nint(dtheta*rtwopi),wp)*twopi

           term  =a*dtheta

           pterm=0.5_wp*term*dtheta
           gamma=term*rsint * rpb1*rpc1

        Else If (keyd == 3) Then

! harmonic cosine dihedral (note sint is cancelled)

           a     =dihedral%param(1,kk)
           theta0=dihedral%param(2,kk)
           dtheta=Cos(theta)-Cos(theta0)

           term  =a*dtheta

           pterm=0.5_wp*term*dtheta
           gamma=-term * rpb1*rpc1

        Else If (keyd == 4) Then

! 3-term cosine dihedral

           a1=dihedral%param(1,kk)
           a2=dihedral%param(2,kk)
           a3=dihedral%param(3,kk)

           pterm=0.5_wp*(a1*(1.0_wp+Cos(theta)) +        &
                         a2*(1.0_wp-Cos(2.0_wp*theta)) + &
                         a3*(1.0_wp+Cos(3.0_wp*theta)) )
           gamma=-0.5_wp*(a1*Sin(theta) -               &
                          2.0_wp*a2*Sin(2.0_wp*theta) + &
                          3.0_wp*a3*Sin(3.0_wp*theta) )*rsint * rpb1*rpc1

        Else If (keyd == 5) Then

! ryckaert-bellemans potential
!
! reference: chem. phys. lett., vol. 30, p. 123 (1975)
! ATTENTION: Modified to have the transition configuration correspond
!            to theta=180 rather than theta=0 as in original form

           a=dihedral%param(1,kk)
           m=Cos(theta)

           pterm=a*( 1.116_wp      - 1.462_wp*m    - 1.578_wp*m**2 + &
                     0.368_wp*m**3 + 3.156_wp*m**4 + 3.788_wp*m**5 )
           gamma=a*( 1.462_wp       + 3.156_wp*m   - 1.104_wp*m**2 - &
                     12.624_wp*m**3 - 18.94_wp*m**4 ) * rpb1*rpc1

        Else If (keyd == 6) Then

! fluorinated ryckaert-bellemans potential
! reference: Rice at al., JCP 104, p. 2101 (1996)

           a=dihedral%param(1,kk)
           m=Cos(theta)
           d=Exp(-56.0_wp*(theta-pi)**2)
           term=-1083.04_wp*(theta-pi)*d

           pterm=a*( 3.55_wp      - 2.78_wp*m    - 3.56_wp*m**2  - &
                     1.64_wp*m**3 + 7.13_wp*m**4 + 12.84_wp*m**5 + &
                     9.67_wp*d )
           gamma=( a*(2.78_wp       + 7.12_wp*m   + 4.92_wp*m**2 - &
                      28.52_wp*m**3 - 64.2_wp*m**4) + term*rsint ) * rpb1*rpc1

        Else If (keyd == 7) Then

! opls cosine dihedral

           a0=dihedral%param(1,kk)
           a1=dihedral%param(2,kk)
           a2=dihedral%param(3,kk)
           a3=dihedral%param(6,kk)
           theta0=dihedral%param(7,kk)
           dtheta=theta-theta0

           pterm=a0 + 0.5_wp*( a1*(1.0_wp+Cos(dtheta)) +        &
                               a2*(1.0_wp-Cos(2.0_wp*dtheta)) + &
                               a3*(1.0_wp+Cos(3.0_wp*dtheta)) )
           gamma=-0.5_wp*(       a1*Sin(dtheta) -               &
                          2.0_wp*a2*Sin(2.0_wp*dtheta)  + &
                          3.0_wp*a3*Sin(3.0_wp*dtheta))*rsint * rpb1*rpc1

        Else If (keyd == 20) Then

! TABDIH potential

           j = dihedral%ltp(kk)
           rdr = dihedral%tab_force(-1,j) ! 1.0_wp/delpot (in rad^-1)

           l   = Int((theta+pi)*rdr)         ! theta (-pi,+pi) is shifted
           ppp = (theta+pi)*rdr - Real(l,wp) ! by +pi so l is [1,ngrid]

           vk  = dihedral%tab_potential(l,j)
           vk1 = dihedral%tab_potential(l+1,j)
           vk2 = dihedral%tab_potential(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           pterm = t1 + (t2-t1)*ppp*0.5_wp

           vk  = dihedral%tab_force(l,j)
           vk1 = dihedral%tab_force(l+1,j)
           vk2 = dihedral%tab_force(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           gamma =-(t1 + (t2-t1)*ppp*0.5_wp)*rsint * rpb1*rpc1

        Else

! flag undefined potential

           safe(1)=.false.
           pterm=0.0_wp
           gamma=0.0_wp

        End If

! calculate atomic forces

        fax = gamma*((-pcy*zbc+pcz*ybc)-pbpc*rpb2*(-pby*zbc+pbz*ybc))
        fay = gamma*(( pcx*zbc-pcz*xbc)-pbpc*rpb2*( pbx*zbc-pbz*xbc))
        faz = gamma*((-pcx*ybc+pcy*xbc)-pbpc*rpb2*(-pbx*ybc+pby*xbc))

        fcx = gamma*((-pcy*zab+pcz*yab)-pbpc*rpb2*(-pby*zab+pbz*yab))
        fcy = gamma*(( pcx*zab-pcz*xab)-pbpc*rpb2*( pbx*zab-pbz*xab))
        fcz = gamma*((-pcx*yab+pcy*xab)-pbpc*rpb2*(-pbx*yab+pby*xab))

        fb1x= gamma*((-pby*zcd+pbz*ycd)-pbpc*rpc2*(-pcy*zcd+pcz*ycd))
        fb1y= gamma*(( pbx*zcd-pbz*xcd)-pbpc*rpc2*( pcx*zcd-pcz*xcd))
        fb1z= gamma*((-pbx*ycd+pby*xcd)-pbpc*rpc2*(-pcx*ycd+pcy*xcd))

        fd1x= gamma*((-pby*zbc+pbz*ybc)-pbpc*rpc2*(-pcy*zbc+pcz*ybc))
        fd1y= gamma*(( pbx*zbc-pbz*xbc)-pbpc*rpc2*( pcx*zbc-pcz*xbc))
        fd1z= gamma*((-pbx*ybc+pby*xbc)-pbpc*rpc2*(-pcx*ybc+pcy*xbc))

        If (ia <= natms) Then

! sum of dihedral energy (dihedral virial is zero!!!)

           engdih=engdih+pterm

! stress tensor calculation for dihedral terms

           strs1 = strs1 + xab*fax + xbc*(fb1x-fcx) - xcd*fd1x
           strs2 = strs2 + yab*fax + ybc*(fb1x-fcx) - ycd*fd1x
           strs3 = strs3 + zab*fax + zbc*(fb1x-fcx) - zcd*fd1x
           strs5 = strs5 + yab*fay + ybc*(fb1y-fcy) - ycd*fd1y
           strs6 = strs6 + yab*faz + ybc*(fb1z-fcz) - ycd*fd1z
           strs9 = strs9 + zab*faz + zbc*(fb1z-fcz) - zcd*fd1z

           fxx(ia)=fxx(ia)+fax
           fyy(ia)=fyy(ia)+fay
           fzz(ia)=fzz(ia)+faz

        End If

        If (ib <= natms) Then

           fxx(ib)=fxx(ib)-fax-fcx+fb1x
           fyy(ib)=fyy(ib)-fay-fcy+fb1y
           fzz(ib)=fzz(ib)-faz-fcz+fb1z

        End If

        If (ic <= natms) Then

           fxx(ic)=fxx(ic)+fcx-fb1x-fd1x
           fyy(ic)=fyy(ic)+fcy-fb1y-fd1y
           fzz(ic)=fzz(ic)+fcz-fb1z-fd1z

        End If

        If (id <= natms) Then

           fxx(id)=fxx(id)+fd1x
           fyy(id)=fyy(id)+fd1y
           fzz(id)=fzz(id)+fd1z

        End If

        xad=xac+xcd
        yad=yac+ycd
        zad=zac+zcd

        rad2=0.0_wp ; rad=0.0_wp

        rad2(0)=xad**2+yad**2+zad**2
        rad(0) =Sqrt(rad2(0))

        If (dihedral%l_core_shell) Then
           If (lad(1,i)) Then
              rad2(1) = xdad(1,i)**2+ydad(1,i)**2+zdad(1,i)**2
              rad(1)  = Sqrt(rad2(1))
           End If
           If (lad(2,i)) Then
              rad2(2) = xdad(2,i)**2+ydad(2,i)**2+zdad(2,i)**2
              rad(2)  = Sqrt(rad2(2))
           End If
           If (lad(3,i)) Then
              rad2(3) = xdad(3,i)**2+ydad(3,i)**2+zdad(3,i)**2
              rad(3)  = Sqrt(rad2(3))
           End If
        End If

! flag error if rad > cutoff

        If (Any(rad > rcut)) Then
           Write(messages(1),*) 'AB',xab,yab,zab
           Write(messages(2),*) 'BC',xbc,ybc,zbc,xac,yac,zac
           Write(messages(3),*) 'CD',xcd,ycd,zcd,xad,yad,zad
           Write(messages(4),*) 'A',xxx(ia),yyy(ia),zzz(ia)
           Write(messages(5),*) 'B',xxx(ib),yyy(ib),zzz(ib)
           Write(messages(6),*) 'C',xxx(ic),yyy(ic),zzz(ic)
           Write(messages(7),*) 'D',xxx(id),yyy(id),zzz(id)
           Call info(messages,7)
           If (dihedral%l_core_shell) Then
              If (lad(1,i)) Then
                Write(message,*) 'A0',xxx(ia0),yyy(ia0),zzz(ia0)
                Call info(message)
              End If
              If (lad(2,i)) Then
                Write(message,*) 'D0',xxx(id0),yyy(id0),zzz(id0)
                Call info(message)
              End If
           End If
           Write(message,*) i,ltg(ia),ltg(ib),ltg(ic),ltg(id),rcut,rad(0)
           Call info(message)
           If (dihedral%l_core_shell) Then
              If (lad(1,i)) Then
                Write(message,*) i,ltg(ia0),ltg(id),rad(1)
                Call info(message)
              End If
              If (lad(2,i)) Then
                Write(message,*) i,ltg(ia),ltg(id0),rad(2)
                Call info(message)
              End If
              If (lad(3,i)) Then
                Write(message,*) i,ltg(ia0),ltg(id0),rad(3)
                Call info(message)
              End If
           End If
           safe(2) = .false.
        End If

! 1-4 electrostatics: adjust by weighting factor
! assumes 1-4 interactions are in the exclude list and Rad < rcut

        scale=dihedral%param(4,kk)

! scaled charge product times dielectric constants

        chgprd=scale*chge(ia)*chge(id)*r4pie0/epsq
        If ((Abs(chgprd) > zero_plus .or. mximpl > 0) .and. keyfce > 0) Then

           If (mximpl > 0) Then
              Call intra_mcoul(keyfce,rcut,alpha,epsq,ia,id,scale, &
                      rad(0),xad,yad,zad,coul,virele,fx,fy,fz,safe(3))
           Else
              Call intra_coul(keyfce,rcut,alpha,epsq,chgprd,rad(0),rad2(0),coul,fcoul,safe(3))

              fx = fcoul*xad
              fy = fcoul*yad
              fz = fcoul*zad

              virele = -fcoul*rad2(0)
           End If

           If (ia <= natms) Then

! correction to electrostatic energy and virial

              engc14 = engc14 + coul
              virc14 = virc14 + virele

! calculate stress tensor

              strs1 = strs1 + xad*fx
              strs2 = strs2 + xad*fy
              strs3 = strs3 + xad*fz
              strs5 = strs5 + yad*fy
              strs6 = strs6 + yad*fz
              strs9 = strs9 + zad*fz

              fxx(ia) = fxx(ia) + fx
              fyy(ia) = fyy(ia) + fy
              fzz(ia) = fzz(ia) + fz

           End If

           If (id <= natms) Then

              fxx(id) = fxx(id) - fx
              fyy(id) = fyy(id) - fy
              fzz(id) = fzz(id) - fz

           End If

        End If

        If (dihedral%l_core_shell) Then
           If (lad(1,i)) Then
              chgprd=scale*chge(ia0)*chge(id)*r4pie0/epsq
              If ((Abs(chgprd) > zero_plus .or. mximpl > 0) .and. keyfce > 0) Then
                 If (mximpl > 0) Then
                    Call intra_mcoul(keyfce,rcut,alpha,epsq,ia0,id,scale, &
                      rad(1),xdad(1,i),ydad(1,i),zdad(1,i),coul,virele,fx,fy,fz,safe(3))
                 Else
                    Call intra_coul(keyfce,rcut,alpha,epsq,chgprd,rad(1),rad2(1),coul,fcoul,safe(3))

                    fx = fcoul*xdad(1,i)
                    fy = fcoul*ydad(1,i)
                    fz = fcoul*zdad(1,i)

                    virele = -fcoul*rad2(1)
                 End If

                 If (ia0 <= natms) Then

! correction to electrostatic energy and virial

                    engc14 = engc14 + coul
                    virc14 = virc14 + virele

! calculate stress tensor

                    strs1 = strs1 + xdad(1,i)*fx
                    strs2 = strs2 + xdad(1,i)*fy
                    strs3 = strs3 + xdad(1,i)*fz
                    strs5 = strs5 + ydad(1,i)*fy
                    strs6 = strs6 + ydad(1,i)*fz
                    strs9 = strs9 + zdad(1,i)*fz

                    fxx(ia0) = fxx(ia0) + fx
                    fyy(ia0) = fyy(ia0) + fy
                    fzz(ia0) = fzz(ia0) + fz

                 End If

                 If (id <= natms) Then

                    fxx(id) = fxx(id) - fx
                    fyy(id) = fyy(id) - fy
                    fzz(id) = fzz(id) - fz

                 End If

              End If
           End If

           If (lad(2,i)) Then
              chgprd=scale*chge(ia)*chge(id0)*r4pie0/epsq
              If ((Abs(chgprd) > zero_plus .or. mximpl > 0) .and. keyfce > 0) Then
                 If (mximpl > 0) Then
                    Call intra_mcoul(keyfce,rcut,alpha,epsq,ia,id0,scale, &
                      rad(2),xdad(2,i),ydad(2,i),zdad(2,i),coul,virele,fx,fy,fz,safe(3))
                 Else
                    Call intra_coul(keyfce,rcut,alpha,epsq,chgprd,rad(2),rad2(2),coul,fcoul,safe(3))

                    fx = fcoul*xdad(2,i)
                    fy = fcoul*ydad(2,i)
                    fz = fcoul*zdad(2,i)

                    virele = -fcoul*rad2(2)
                 End If

                 If (ia <= natms) Then

! correction to electrostatic energy and virial

                    engc14 = engc14 + coul
                    virc14 = virc14 + virele

! calculate stress tensor

                    strs1 = strs1 + xdad(2,i)*fx
                    strs2 = strs2 + xdad(2,i)*fy
                    strs3 = strs3 + xdad(2,i)*fz
                    strs5 = strs5 + ydad(2,i)*fy
                    strs6 = strs6 + ydad(2,i)*fz
                    strs9 = strs9 + zdad(2,i)*fz

                    fxx(ia) = fxx(ia) + fx
                    fyy(ia) = fyy(ia) + fy
                    fzz(ia) = fzz(ia) + fz

                 End If

                 If (id0 <= natms) Then

                    fxx(id0) = fxx(id0) - fx
                    fyy(id0) = fyy(id0) - fy
                    fzz(id0) = fzz(id0) - fz

                 End If
              End If
           End If

           If (lad(3,i)) Then
              chgprd=scale*chge(ia0)*chge(id0)*r4pie0/epsq
              If ((Abs(chgprd) > zero_plus .or. mximpl > 0) .and. keyfce > 0) Then
                 If (mximpl > 0) Then
                    Call intra_mcoul(keyfce,rcut,alpha,epsq,ia0,id0,scale, &
                      rad(3),xdad(3,i),ydad(3,i),zdad(3,i),coul,virele,fx,fy,fz,safe(3))
                 Else
                    Call intra_coul(keyfce,rcut,alpha,epsq,chgprd,rad(3),rad2(3),coul,fcoul,safe(3))

                    fx = fcoul*xdad(3,i)
                    fy = fcoul*ydad(3,i)
                    fz = fcoul*zdad(3,i)

                    virele = -fcoul*rad2(3)
                 End If

                 If (ia0 <= natms) Then

! correction to electrostatic energy and virial

                    engc14 = engc14 + coul
                    virc14 = virc14 + virele

! calculate stress tensor

                    strs1 = strs1 + xdad(3,i)*fx
                    strs2 = strs2 + xdad(3,i)*fy
                    strs3 = strs3 + xdad(3,i)*fz
                    strs5 = strs5 + ydad(3,i)*fy
                    strs6 = strs6 + ydad(3,i)*fz
                    strs9 = strs9 + zdad(3,i)*fz

                    fxx(ia0) = fxx(ia0) + fx
                    fyy(ia0) = fyy(ia0) + fy
                    fzz(ia0) = fzz(ia0) + fz

                 End If

                 If (id0 <= natms) Then

                    fxx(id0) = fxx(id0) - fx
                    fyy(id0) = fyy(id0) - fy
                    fzz(id0) = fzz(id0) - fz

                 End If
              End If
           End If
        End If

! 1-4 short-range (vdw) interactions: adjust by weighting factor
! assumes 1-4 interactions are in the exclude list and Rad < rvdw

        scale=dihedral%param(5,kk)
        If (Abs(scale) > zero_plus .and. ntpvdw > 0) Then

! atomic type indices

           ai=ltype(ia)
           aj=ltype(id)

           Call dihedrals_14_vdw(rvdw,ai,aj,rad(0),rad2(0),eng,gamma)

           gamma = scale*gamma
           eng   = scale*eng

           fx = gamma*xad
           fy = gamma*yad
           fz = gamma*zad

           If (ia <= natms) Then

! add scaled 1-4 short-range potential energy and virial

              engs14=engs14+eng
              virs14=virs14-gamma*rad2(0)

! calculate stress tensor

              strs1 = strs1 + xad*fx
              strs2 = strs2 + xad*fy
              strs3 = strs3 + xad*fz
              strs5 = strs5 + yad*fy
              strs6 = strs6 + yad*fz
              strs9 = strs9 + zad*fz

              fxx(ia) = fxx(ia) + fx
              fyy(ia) = fyy(ia) + fy
              fzz(ia) = fzz(ia) + fz

           End If

           If (id <= natms) Then

              fxx(id) = fxx(id) - fx
              fyy(id) = fyy(id) - fy
              fzz(id) = fzz(id) - fz

           End If

           If (dihedral%l_core_shell) Then
              If (lad(1,i)) Then
                 ai=ltype(ia0)
                 aj=ltype(id)

                 Call dihedrals_14_vdw(rvdw,ai,aj,rad(1),rad2(1),eng,gamma)

                 gamma = scale*gamma
                 eng   = scale*eng

                 fx = gamma*xdad(1,i)
                 fy = gamma*ydad(1,i)
                 fz = gamma*zdad(1,i)

                 If (ia0 <= natms) Then

! add scaled 1-4 short-range potential energy and virial

                    engs14=engs14+eng
                    virs14=virs14-gamma*rad2(1)

! calculate stress tensor

                    strs1 = strs1 + xdad(1,i)*fx
                    strs2 = strs2 + xdad(1,i)*fy
                    strs3 = strs3 + xdad(1,i)*fz
                    strs5 = strs5 + ydad(1,i)*fy
                    strs6 = strs6 + ydad(1,i)*fz
                    strs9 = strs9 + zdad(1,i)*fz

                    fxx(ia0) = fxx(ia0) + fx
                    fyy(ia0) = fyy(ia0) + fy
                    fzz(ia0) = fzz(ia0) + fz

                 End If

                 If (id <= natms) Then

                    fxx(id) = fxx(id) - fx
                    fyy(id) = fyy(id) - fy
                    fzz(id) = fzz(id) - fz

                 End If
              End If

              If (lad(2,i)) Then
                 ai=ltype(ia)
                 aj=ltype(id0)

                 Call dihedrals_14_vdw(rvdw,ai,aj,rad(2),rad2(2),eng,gamma)

                 gamma = scale*gamma
                 eng   = scale*eng

                 fx = gamma*xdad(2,i)
                 fy = gamma*ydad(2,i)
                 fz = gamma*zdad(2,i)

                 If (ia <= natms) Then

! add scaled 1-4 short-range potential energy and virial

                    engs14=engs14+eng
                    virs14=virs14-gamma*rad2(2)

! calculate stress tensor

                    strs1 = strs1 + xdad(2,i)*fx
                    strs2 = strs2 + xdad(2,i)*fy
                    strs3 = strs3 + xdad(2,i)*fz
                    strs5 = strs5 + ydad(2,i)*fy
                    strs6 = strs6 + ydad(2,i)*fz
                    strs9 = strs9 + zdad(2,i)*fz

                    fxx(ia) = fxx(ia) + fx
                    fyy(ia) = fyy(ia) + fy
                    fzz(ia) = fzz(ia) + fz

                 End If

                 If (id0 <= natms) Then

                    fxx(id0) = fxx(id0) - fx
                    fyy(id0) = fyy(id0) - fy
                    fzz(id0) = fzz(id0) - fz

                 End If
              End If

              If (lad(3,i)) Then
                 ai=ltype(ia0)
                 aj=ltype(id0)

                 Call dihedrals_14_vdw(rvdw,ai,aj,rad(3),rad2(3),eng,gamma)

                 gamma = scale*gamma
                 eng   = scale*eng

                 fx = gamma*xdad(3,i)
                 fy = gamma*ydad(3,i)
                 fz = gamma*zdad(3,i)

                 If (ia0 <= natms) Then

! add scaled 1-4 short-range potential energy and virial

                    engs14=engs14+eng
                    virs14=virs14-gamma*rad2(3)

! calculate stress tensor

                    strs1 = strs1 + xdad(3,i)*fx
                    strs2 = strs2 + xdad(3,i)*fy
                    strs3 = strs3 + xdad(3,i)*fz
                    strs5 = strs5 + ydad(3,i)*fy
                    strs6 = strs6 + ydad(3,i)*fz
                    strs9 = strs9 + zdad(3,i)*fz

                    fxx(ia0) = fxx(ia0) + fx
                    fyy(ia0) = fyy(ia0) + fy
                    fzz(ia0) = fzz(ia0) + fz

                 End If

                 If (id0 <= natms) Then

                    fxx(id0) = fxx(id0) - fx
                    fyy(id0) = fyy(id0) - fy
                    fzz(id0) = fzz(id0) - fz

                 End If
              End If
           End If

        End If

     End If
  End Do

  If (Mod(isw,3) > 0) Then

! sum contributions to potentials


        buffer(1) = engdih
        buffer(2) = engc14
        buffer(3) = virc14
        buffer(4) = engs14
        buffer(5) = virs14

        Call gsum(comm,buffer(1:5))

        engdih = buffer(1)
        engc14 = buffer(2)
        virc14 = buffer(3)
        engs14 = buffer(4)
        virs14 = buffer(5)


     engcpe = engcpe + engc14
     vircpe = vircpe + virc14
     engsrp = engsrp + engs14
     virsrp = virsrp + virs14

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

! check safety to continue

     Call gcheck(comm,safe)
     If (.not.safe(1)) Call error(448)
     If (.not.safe(2)) Call error(445)
     If (.not.safe(3)) Call error(446)

  End If

  Deallocate (lunsafe,lstopt,lad, Stat=fail(1))
  Deallocate (xdab,ydab,zdab,     Stat=fail(2))
  Deallocate (xdbc,ydbc,zdbc,     Stat=fail(3))
  Deallocate (xdcd,ydcd,zdcd,     Stat=fail(4))
  Deallocate (xdad,ydad,zdad,     Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'dihedrals_forces deallocation failure'
     Call error(0,message)
  End If

End Subroutine dihedrals_forces

Subroutine dihedrals_table_read(dihd_name,dihedral,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABDIH file (for dihedral potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov january 2017
! contrib   - a.v.brukhno, february 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Character( Len = 32 ), Intent( In    ) :: dihd_name(1:dihedral%max_types)
  Type( comms_type), Intent( InOut ) :: comm

  Logical                :: safe,remake,zero
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 32 )  :: iddihd
  Character( Len = 8   ) :: atom1,atom2,atom3,atom4

  Integer                :: fail(1:2),ngrid,rtdih,itdih,jtdih,katom1,katom2,katom3,katom4,jtpatm,i,l
  Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,rrr0, &
                            ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

  Integer,           Allocatable :: read_type(:)
  Real( Kind = wp ), Allocatable :: bufpot(:),bufvir(:)
  Character( Len = 256 ) :: message,messages(5)

  If (comm%idnode == 0) Open(Unit=ntable, File='TABDIH')

! skip header record

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

! read mesh resolution not needed for dihedral angle dependent
! potentials/forces as delpot=360/ngrid running from -180 to 180

  Call get_line(safe,ntable,record,comm)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABDIH if it's in .xvg format

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  delpot = 360.0_wp/Real(ngrid,wp)

  dlrpot = 360.0_wp/Real(dihedral%bin_tab-4,wp)

! check grid spacing

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > delth_max .and. (.not.safe)) Then
     Write(messages(1),*) ''
     Write(messages(2),'(a,1p,e15.7)') 'expected (maximum) angular increment : ', delth_max
     Write(messages(3),'(a,1p,e15.7)') 'TABDIH file actual angular increment : ', delpot
     Write(messages(4),'(a,0p,i10)') 'expected (minimum) number of grid points : ', dihedral%bin_tab-4
     Write(messages(5),'(a,0p,i10)') 'TABDIH file actual number of grid points : ', ngrid
     Call info(messages,5,.true.)
     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     Write(message,'(a,i10)') 'TABDIH arrays resized for mxgrid = ', dihedral%bin_tab-4
     Call info('',.true.)
     Call info(message,.true.)
  End If

! compare grids dimensions

  If (ngrid < dihedral%bin_tab-4) Then
     Call warning(270,Real(ngrid,wp),Real(dihedral%bin_tab-4,wp),0.0_wp)
     Call error(48)
  End If

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  fail=0
  Allocate (read_type(1:dihedral%ltp(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'error - dihedrals_table_read allocation failure'
     Call error(0,message)
  End If
  Call allocate_dihd_pot_arrays(dihedral)

  read_type=0 ! initialise read_type
  Do rtdih=1,dihedral%ltp(0)
     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record,comm)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABDIH if it's in .xvg format

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
        Write(message,'(a,i0,a,i0,a,i0,a,i0,a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABDIH'
        Call error(91,message,.true.)
     End If

! Construct unique name for the tabulated dihedral

     If (katom1 <= katom4) Then
        iddihd = atom1//atom2//atom3//atom4
     Else
        iddihd = atom4//atom3//atom2//atom1
     End If

! read potential arrays if potential is defined

     itdih=0
     Do jtdih=1,dihedral%ltp(0)
        If (dihd_name(jtdih) == iddihd) Then
           Do itdih=1,dihedral%max_types
              If (dihedral%ltp(itdih) == jtdih) Exit
           End Do
           Exit
        End If
     End Do

     If (itdih == 0) Then ! All(dihd_name /= iddihd)
        Write(message,'(a,i0,a,i0,a,i0,a,i0,a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABDIH'
        Call error(89,message,.true.)
     End If
     If (Any(read_type == jtdih)) Then
        Write(message,'(a,i0,a,i0,a,i0,a,i0,a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABDIH'
        Call error(172,message,.true.)
     Else
        read_type(jtdih)=jtdih
     End If

! read in potential & force arrays

     Do i=0,2
        bufpot(i) = 0.0_wp
        bufvir(i) = 0.0_wp
     End Do

! read in the zero and/or first & second data elements (potential & virial)

     zero=.false.
     If (comm%idnode == 0) Then
        rrr=0.0_wp
        Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

        If (Abs(rrr) > zero_plus) Then ! no zero element data => extrapolate to zero
           bufpot(1) = bufp0
           bufvir(1) = bufv0
           rrr0      = rrr

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              Write(messages(1),*) ''
              Write(messages(2),'(a,1p,e15.7)') 'TABDIH stated  angular increment : ', delpot
              Write(messages(3),'(a,1p,e15.7)') 'TABDIH read-in angular increment : ', rrr-rrr0
              Call info(messages,3,.true.)
           End If

           bufpot(2) = bufp0
           bufvir(2) = bufv0
        Else ! zero element data found => read in the first element for checking delr
           zero=.true.
           bufpot(0) = bufp0
           bufvir(0) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              Write(messages(1),*) ''
              Write(messages(2),'(a,1p,e15.7)') 'TABDIH stated  angular increment : ', delpot
              Write(messages(3),'(a,1p,e15.7)') 'TABDIH read-in angular increment : ',rrr
              Call info(messages,3,.true.)
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

! linear extrapolation for grid point 0 at distances close to 0 -
! midpoint for a periodic function

     If (.not.zero) Then
        bufpot(0) = 0.5_wp*(bufpot(1)-bufpot(ngrid))
        bufvir(0) = 0.5_wp*(bufvir(1)-bufvir(ngrid))
     End If

! reconstruct arrays using 3pt interpolation

     If (remake) Then
        Do i=0,dihedral%bin_tab-4
           rrr = Real(i,wp)*delth_max
           l   = Int(rrr*rdr)
           ppp = rrr*rdr-Real(l,wp)

! cyclic grid

           If (l  <= ngrid) Then
              vk  = bufpot(l)
           Else
              vk  = bufpot(Mod(l,ngrid+1))
           End If
           If (l+1 <= ngrid) Then
              vk1 = bufpot(l+1)
           Else
              vk1 = bufpot(Mod(l+1,ngrid+1))
           End If
           If (l+2 <= ngrid) Then
              vk2 = bufpot(l+2)
           Else
              vk2 = bufpot(Mod(l+2,ngrid+1))
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           dihedral%tab_potential(i,jtdih) = t1 + (t2-t1)*ppp*0.5_wp
           dihedral%tab_potential(i,jtdih) = dihedral%tab_potential(i,jtdih)*engunit ! convert to internal units

! cyclic grid

           If (l  <= ngrid) Then
              vk  = bufvir(l)
           Else
              vk  = bufvir(Mod(l,ngrid+1))
           End If
           If (l+1 <= ngrid) Then
              vk1 = bufvir(l+1)
           Else
              vk1 = bufvir(Mod(l+1,ngrid+1))
           End If
           If (l+2 <= ngrid) Then
              vk2 = bufvir(l+2)
           Else
              vk2 = bufvir(Mod(l+2,ngrid+1))
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           dihedral%tab_force(i,jtdih) = t1 + (t2-t1)*ppp*0.5_wp
           dihedral%tab_force(i,jtdih) = dihedral%tab_force(i,jtdih)*engunit*rad2dgr ! convert to internal units
        End Do

        dihedral%tab_force(-1,jtdih) = rad2dgr/dlrpot
     Else
        Do i=0,dihedral%bin_tab-4
           dihedral%tab_potential(i,jtdih) = bufpot(i)*engunit ! convert to internal units
           dihedral%tab_force(i,jtdih) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do

! cyclic grid

        dihedral%tab_potential(dihedral%bin_tab-3,jtdih) = dihedral%tab_potential(0,jtdih)
        dihedral%tab_potential(dihedral%bin_tab-2,jtdih) = dihedral%tab_potential(1,jtdih)
        dihedral%tab_force(dihedral%bin_tab-3,jtdih) = dihedral%tab_force(0,jtdih)
        dihedral%tab_force(dihedral%bin_tab-2,jtdih) = dihedral%tab_force(1,jtdih)

        dihedral%tab_force(-1,jtdih) = rad2dgr/delpot
     End If
  End Do

  If (comm%idnode == 0) Then
     Close(Unit=ntable)
  End If
  Call info('',.true.)
  Call info('potential tables read from TABDIH file',.true.)

! Break if not safe

  Call gcheck(comm,safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'error - dihedrals_table_read deallocation failure'
     Call error(0,message)
  End If

  Return

! end of file error exit

100 Continue

  If (comm%idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine dihedrals_table_read

Subroutine dihedrals_14_vdw(rvdw,ai,aj,rad,rad2,eng,gamma)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating dihedrals 1-4 vdw interaction:
! adjust by weighting factor
!
! copyright - daresbury laboratory
! amended   - i.t.todorov january 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: ai,aj
  Real( Kind = wp ), Intent( In    ) :: rvdw,rad,rad2
  Real( Kind = wp ), Intent(   Out ) :: eng,gamma

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: dlrpot,rdr

  Integer           :: key,k,l,ityp,n,m
  Real( Kind = wp ) :: rsq,rrr,ppp,                 &
                       r0,r0rn,r0rm,r_6,sor6,       &
                       rho,a,b,c,d,e0,kk,           &
                       nr,mr,rc,sig,eps,alpha,beta, &
                       gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3

  If (newjob) Then
     newjob = .false.

! define grid resolution for potential arrays and interpolation spacing

     dlrpot = rvdw/Real(mxgvdw-4,wp)
     rdr    = 1.0_wp/dlrpot
  End If

! Zero energy and force components

  eng   = 0.0_wp
  gamma = 0.0_wp

! potential function indices

  If (ai > aj) Then
     key=ai*(ai-1)/2 + aj
  Else
     key=aj*(aj-1)/2 + ai
  End If

  k=lstvdw(key)

! validity of potential

  ityp=ltpvdw(k)
  If (ityp >= 0) Then

! Get separation distance

     rrr = rad
     rsq = rad2

     If (ld_vdw) Then ! direct calculation

        If      (ityp == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

           a=prmvdw(1,k)
           b=prmvdw(2,k)

           r_6=rrr**(-6)

           eng   = r_6*(a*r_6-b)
           gamma = 6.0_wp*r_6*(2.0_wp*a*r_6-b)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

           eps=prmvdw(1,k)
           sig=prmvdw(2,k)

           sor6=(sig/rrr)**6

           eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)
           gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

           e0=prmvdw(1,k)
           n =Nint(prmvdw(2,k)) ; nr=Real(n,wp)
           m =Nint(prmvdw(3,k)) ; mr=Real(m,wp)
           r0=prmvdw(4,k)

           a=r0/rrr
           b=1.0_wp/(nr-mr)
           r0rn=a**n
           r0rm=a**m

           eng   = e0*(mr*r0rn-nr*r0rm)*b
           gamma = e0*mr*nr*(r0rn-r0rm)*b/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

           a  =prmvdw(1,k)
           rho=prmvdw(2,k)
           c  =prmvdw(3,k)

           If (Abs(rho) <= zero_plus) Then
              If (Abs(a) <= zero_plus) Then
                 rho=1.0_wp
              Else
                 Call error(467)
              End If
           End If

           b=rrr/rho
           t1=a*Exp(-b)
           t2=-c/rrr**6

           eng   = t1+t2
           gamma = (t1*b+6.0_wp*t2)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

           a  =prmvdw(1,k)
           b  =prmvdw(2,k)
           sig=prmvdw(3,k)
           c  =prmvdw(4,k)
           d  =prmvdw(5,k)

           t1=a*Exp(b*(sig-rrr))
           t2=-c/rrr**6
           t3=-d/rrr**8

           eng   = t1+t2+t3
           gamma = (t1*rrr*b+6.0_wp*t2+8.0_wp*t3)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

           a=prmvdw(1,k)
           b=prmvdw(2,k)

           t1=a/rrr**12
           t2=-b/rrr**10

           eng   = t1+t2
           gamma = (12.0_wp*t1+10.0_wp*t2)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

           e0=prmvdw(1,k)
           n =Nint(prmvdw(2,k)) ; nr=Real(n,wp)
           m =Nint(prmvdw(3,k)) ; mr=Real(m,wp)
           r0=prmvdw(4,k)
           rc=prmvdw(5,k) ; If (rc < 1.0e-6_wp) rc=rvdw

           If (n <= m) Call error(470)

           b=1.0_wp/(nr-mr)
           c=rc/r0 ; If (c < 1.0_wp) Call error(468)

           beta = c*( (c**(m+1)-1.0_wp) / (c**(n+1)-1.0_wp) )**b
           alpha= -(nr-mr) / (  mr*(beta**n)*(1.0_wp+(nr/c-nr-1.0_wp)/c**n) &
                               -nr*(beta**m)*(1.0_wp+(mr/c-mr-1.0_wp)/c**m) )
           e0 = e0*alpha

           If (rrr <= rc) Then
              a=r0/rrr

              eng   = e0*(  mr*(beta**n)*(a**n-(1.0_wp/c)**n) &
                           -nr*(beta**m)*(a**m-(1.0_wp/c)**m) &
                           +nr*mr*((rrr/rc-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
              gamma = e0*mr*nr*( (beta**n)*a**n-(beta**m)*a**m &
                               -rrr/rc*((beta/c)**n-(beta/c)**m) )*b/rsq
           End If

        Else If (ityp == 8) Then

! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}

           e0=prmvdw(1,k)
           r0=prmvdw(2,k)
           kk=prmvdw(3,k)

           t1=Exp(-kk*(rrr-r0))

           eng   = e0*t1*(t1-2.0_wp)
           gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)/rrr

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (ityp == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

           eps=prmvdw(1,k)
           sig=prmvdw(2,k)
           d  =prmvdw(3,k)

           If (rrr < prmvdw(4,k) .or. Abs(rrr-d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6=(sig/(rrr-d))**6

              eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rrr*(rrr-d))

              If (ls_vdw) Then ! force-shifting
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
              End If
           End If

        Else If (ityp == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

           a =prmvdw(1,k)
           rc=prmvdw(2,k)

           If (rrr < rc) Then ! Else leave them zeros
              t2=rrr/rc
              t1=0.5_wp*a*rrr*(1.0_wp-t2)

              eng   = t1*(1.0_wp-t2)
              gamma = t1*(3.0_wp*t2-1.0_wp)/rsq
           End If

        Else If (ityp == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

           eps=prmvdw(1,k)
           sig=prmvdw(2,k)

           rho=sig/rrr
           t1=1.0_wp/(0.07_wp+rho)
           t2=1.0_wp/(0.12_wp+rho**7)
           t3=eps*(1.07_wp/t1**7)

           eng   = t3*((1.12_wp/t2)-2.0_wp)
           gamma =-7.0_wp*t3*rho*(((1.12_wp/t2)-2.0_wp)/t1 + (1.12_wp/t2**2)*rho**6)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma - afs(k)/rrr
           End If

        Else If (Abs(vvdw(0,k)) > zero_plus) Then ! potential read from TABLE - (ityp == 0)

           l   = Int(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

           vk  = vvdw(l,k)
           vk1 = vvdw(l+1,k)
           vk2 = vvdw(l+2,k)

           t1 = vk  + (vk1 - vk )*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           eng = t1 + (t2-t1)*ppp*0.5_wp
           If (ls_vdw) eng = eng + gvdw(mxgvdw-4,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgvdw-4,k) ! force-shifting

! calculate forces using 3-point interpolation

           gk  = gvdw(l,k) ; If (l == 0) gk = gk*rrr
           gk1 = gvdw(l+1,k)
           gk2 = gvdw(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rsq
           If (ls_vdw) gamma = gamma - gvdw(mxgvdw-4,k)/(rrr*rvdw) ! force-shifting

        End If

     Else If (Abs(vvdw(0,k)) > zero_plus) Then ! no direct = fully tabulated calculation

        l   = Int(rrr*rdr)
        ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

        vk  = vvdw(l,k)
        vk1 = vvdw(l+1,k)
        vk2 = vvdw(l+2,k)

        t1 = vk  + (vk1 - vk )*ppp
        t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

        eng = t1 + (t2-t1)*ppp*0.5_wp
        If (ls_vdw) eng = eng + gvdw(mxgvdw-4,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgvdw-4,k) ! force-shifting

! calculate forces using 3-point interpolation

        gk  = gvdw(l,k) ; If (l == 0) gk = gk*rrr
        gk1 = gvdw(l+1,k)
        gk2 = gvdw(l+2,k)

        t1 = gk  + (gk1 - gk )*ppp
        t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

        gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rsq
        If (ls_vdw) gamma = gamma - gvdw(mxgvdw-4,k)/(rrr*rvdw) ! force-shifting

     End If

  End If

End Subroutine dihedrals_14_vdw

  Subroutine cleanup(dihedral)
    Type(dihedrals_type) :: dihedral

    If (Allocated(dihedral%num)) Then
      Deallocate(dihedral%num)
    End If
    If (Allocated(dihedral%key)) Then
      Deallocate(dihedral%key)
    End If

    If (Allocated(dihedral%lst)) Then
      Deallocate(dihedral%lst)
    End If
    If (Allocated(dihedral%list)) Then
      Deallocate(dihedral%list)
    End If
    If (Allocated(dihedral%legend)) Then
      Deallocate(dihedral%legend)
    End If

    If (Allocated(dihedral%param)) Then
      Deallocate(dihedral%param)
    End If

    If (Allocated(dihedral%ltp)) Then
      Deallocate(dihedral%ltp)
    End If
    If (Allocated(dihedral%tab_potential)) Then
      Deallocate(dihedral%tab_potential)
    End If
    If (Allocated(dihedral%tab_force)) Then
      Deallocate(dihedral%tab_force)
    End If
  End Subroutine cleanup
End Module dihedrals
