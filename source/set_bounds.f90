Subroutine set_bounds                                 &
           (levcfg,l_str,lsim,l_vv,l_n_e,l_n_v,l_ind, &
           dvar,rcut,rpad,rlnk,rvdw,rmet,rbin,nstfce,alpha,width)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to determine various limits as array bounds,
! grid sizes, paddings, iterations, etc. as specified in setup_module
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2016
! contrib   - i.j.bush february 2014
! contrib   - m.a.seaton june 2014 (VAF)
! contrib   - m.a.seaton march 2017 (TTM)
! contrib   - i.t.todorov march 2018 (rpad reset if 'strict' & rpad undefined)
! contrib   - i.t.todorov june 2018 (spme suggestions)
! contrib   - i.t.todorov june 2018 (fdens & mxatms fixes)
! contrib   - i.t.todorov august 2018 (rpad refinements)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode
  Use setup_module
  Use domains_module,     Only : map_domains,nprx,npry,nprz,r_nprx,r_npry,r_nprz
  Use config_module,      Only : imcon,imc_n,cfgname,cell,volm
  Use vnl_module,         Only : llvnl ! Depends on l_str,lsim & rpad
  Use msd_module
  Use rdf_module,         Only : rusr
  Use kim_module,         Only : kim
  Use bonds_module,       Only : rcbnd
  Use tersoff_module,     Only : potter
  Use development_module, Only : l_trm
  Use greenkubo_module,   Only : vafsamp
  Use mpoles_module,      Only : keyind,induce
  Use ttm_module,         Only : delx,dely,delz,volume,rvolume,ntsys,eltsys,redistribute,sysrho,l_ttm

  Implicit None

  Logical,           Intent(   Out ) :: l_str,lsim,l_vv,l_n_e,l_n_v,l_ind
  Integer,           Intent(   Out ) :: levcfg,nstfce
  Real( Kind = wp ), Intent(   Out ) :: dvar,rcut,rpad,rlnk
  Real( Kind = wp ), Intent(   Out ) :: rvdw,rmet,rbin,alpha,width

  Logical           :: l_usr,l_n_r,lzdn,lext,lrpad,lrpad0=.false.
  Integer           :: megatm,i,ilx,ily,ilz,qlx,qly,qlz, &
                       mtshl,mtcons,mtrgd,mtteth,mtbond,mtangl,mtdihd,mtinv
  Real( Kind = wp ) :: ats,celprp(1:10),cut,   &
                       dens0,dens,fdens,fdvar, &
                       test,vcell,tol,         &
                       rpad1,rpad2,            &
                       rcter,rctbp,rcfbp,      &
                       xhi,yhi,zhi

! define zero+ and half+/- (setup_module)

  zero_plus  = Nearest( 0.0_wp , +1.0_wp)
  half_plus  = Nearest( 0.5_wp , +1.0_wp)
  half_minus = Nearest( 0.5_wp , -1.0_wp)

! scan the FIELD file data

  Call scan_field                                    &
           (l_n_e,mxompl,mximpl,                     &
           mxsite,mxatyp,megatm,mxtmls,mxexcl,       &
           mtshl,mxtshl,mxshl,mxfshl,                &
           mtcons,mxtcon,mxcons,mxfcon,              &
           mxtpmf,mxpmf,mxfpmf,l_usr,                &
           mtrgd,mxtrgd,mxrgd,mxlrgd,mxfrgd,         &
           mtteth,mxtteth,mxteth,mxftet,             &
           mtbond,mxtbnd,mxbond,mxfbnd,rcbnd,mxgbnd, &
           mtangl,mxtang,mxangl,mxfang,mxgang,       &
           mtdihd,mxtdih,mxdihd,mxfdih,mxgdih,       &
           mtinv,mxtinv,mxinv,mxfinv,mxginv,         &
           mxrdf,mxvdw,rvdw,mxgvdw,                  &
           mxmet,mxmed,mxmds,rmet,mxgmet,            &
           mxter,rcter,mxtbp,rctbp,mxfbp,rcfbp,lext)

! Get imc_r & set dvar

  Call scan_control_pre(imc_n,dvar)

! scan CONFIG file data

  Call scan_config(megatm,imc_n,dvar,cfgname,levcfg,imcon,cell,xhi,yhi,zhi)

! halt execution for unsupported image conditions in DD
! checks for some inherited from DL_POLY_2 are though kept

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(514)

! scan CONTROL file data

  Call scan_control                                        &
           (rcbnd,mxrdf,mxvdw,rvdw,mxmet,rmet,mxter,rcter, &
           mxrgd,imcon,imc_n,cell,xhi,yhi,zhi,             &
           mxgana,mxgbnd1,mxgang1,mxgdih1,mxginv1,         &
           l_str,lsim,l_vv,l_n_e,l_n_r,lzdn,l_n_v,l_ind,   &
           rcut,lrpad,rpad,rbin,mxstak,                    &
           mxshl,mxompl,mximpl,keyind,                     &
           nstfce,mxspl,alpha,kmaxa1,kmaxb1,kmaxc1)

! check integrity of cell vectors: for cubic, TO and RD cases
! i.e. cell(1)=cell(5)=cell(9) (or cell(9)/Sqrt(2) for RD)

  If (imcon == 1 .or. imcon == 4 .or. imcon == 5) Then

     ats = (Abs(cell(1))+Abs(cell(5)))/2.0_wp
     test = 1.0e-10_wp*ats
     If (Abs(cell(1)-ats) > test) Call error(410)
     If (Abs(cell(5)-ats) > test) Call error(410)
     If (imcon == 5) Then
        If (Abs(cell(9)-ats*rt2) > test) Call error(410)
     Else
        If (Abs(cell(9)-ats) > test) Call error(410)
     End If
  End If

! check integrity of hexagonal prism cell vectors

  If (imcon == 7) Then
     If (Abs(cell(1)-rt3*cell(5)) > 1.0e-6_wp) Call error(410)
  End If

! check for diagonal cell matrix if appropriate: imcon=1,2,4,5,7

  If (imcon /= 0 .and. imcon /= 3 .and. imcon /= 6) Then
     If (Abs(cell(2)) > zero_plus) Call error(410)
     If (Abs(cell(3)) > zero_plus) Call error(410)
     If (Abs(cell(4)) > zero_plus) Call error(410)
     If (Abs(cell(6)) > zero_plus) Call error(410)
     If (Abs(cell(7)) > zero_plus) Call error(410)
     If (Abs(cell(8)) > zero_plus) Call error(410)
  End If

! calculate dimensional properties of simulation cell
! (for use in link-cells) and volume and define min cell width

  Call dcell(cell,celprp)
  width=Min(celprp(7),celprp(8),celprp(9))

  volm = celprp(10)

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) volm=0.5_wp*volm

! check value of cutoff and reset if necessary

  If (imcon > 0) Then
     If (imcon == 4) width=rt3*cell(1)/2.0_wp
     If (imcon == 5) width=cell(1)
     If (imcon == 6) width=Min(celprp(7),celprp(8))

! halt program if potential cutoff exceeds the minimum half-cell width

     If (rcut > width/2.0_wp) Then
        Call warning(3,rcut,width/2.0_wp,0.0_wp)
        Call error(95)
     End If
  End If


! dvar function

  fdvar = dvar**1.7_wp

! dvar push of mxfdih and mxexcl ranges as the usual suspects

  mxfdih = Nint(fdvar * Real(mxfdih,wp))
  mxexcl = Nint(fdvar * Real(mxexcl,wp))

!!! INTRA-LIKE POTENTIAL PARAMETERS !!!

! maximum number of core-shell units per node

  If (mxshl > 0 .and. mxnode > 1) Then
     mxshl = Max(mxshl,mxnode*mtshl)
     mxshl = (3*(Nint(fdvar*Real(mxshl,wp))+mxnode-1))/mxnode
  End If


! maximum number of constraints per node

  If (mxcons > 0 .and. mxnode > 1) Then
     mxcons = Max(mxcons,mxnode*mtcons)
     mxcons = (3*(Nint(fdvar*Real(mxcons,wp))+mxnode-1))/mxnode
  End If


! maximum number of PMF constraints per MD cell - mxpmf
! (only one type of PMF in only one type of molecule !!!)
! maximum number of atoms per PMF unit - mxtpmf(1:2) (only two units per pmf)


! maximum number of RBs per node

  If (mxrgd > 0 .and. mxnode > 1) Then
     mxrgd = Max(mxrgd,mxnode*mtrgd)
     mxrgd = (3*(Nint(fdvar*Real(mxrgd,wp))+mxnode-1))/mxnode
  End If


! dimension of shared atoms arrays for core-shell, constraint and RB units
! Max=Max#(members-per-unit)*Max#(units-per-domain)/2
! and maximum number of neighbouring domains/nodes in 3D DD (3^3 - 1)

  If (mxnode > 1) Then
     mxlshp = Max((2*mxshl)/2,(2*mxcons)/2,(mxlrgd*mxrgd)/2)
     mxproc = 26
  Else ! nothing is to be shared on one node
     mxlshp = 0
     mxproc = 0
  End If


! maximum number of tethered atoms per node and tether potential parameters

  If (mxteth > 0) Then
     If (mxnode > 1) Then
        mxteth = Max(mxteth,mxnode*mtteth)
        mxteth = (3*(Nint(fdvar*Real(mxteth,wp))+mxnode-1))/mxnode
     End If
     mxpteth = 3
  Else
     mxpteth = 0
  End If


! maximum number of chemical bonds per node and bond potential parameters

  If (mxbond > 0) Then
     If (mxnode > 1) Then
        mxbond = Max(mxbond,mxnode*mtbond)
        mxbond = (3*(Nint(fdvar*Real(mxbond,wp))+mxnode-1))/mxnode
     End If
     mxpbnd = 4
  Else
     mxpbnd = 0
  End If


! maximum number of bond angles per node and angular potential parameters

  If (mxangl > 0) Then
     If (mxnode > 1) Then
        mxangl = Max(mxangl,mxnode*mtangl)
        mxangl = (3*(Nint(fdvar*Real(mxangl,wp))+mxnode-1))/mxnode
     End If
     mxpang = 6
  Else
     mxpang = 0
  End If


! maximum number of torsion angles per node and dihedral potential parameters

  If (mxdihd > 0) Then
     If (mxnode > 1) Then
        mxdihd = Max(mxdihd,mxnode*mtdihd)
        mxdihd = (3*(Nint(fdvar*Real(mxdihd,wp))+mxnode-1))/mxnode
        mxdihd = mxdihd + (mxdihd+4)/5 ! allow for 25% higher density
     End If
     mxpdih = 7
  Else
     mxpdih = 0
  End If


! maximum number of inversions per node and inversion potential parameters
! allow for 20% higher density

  If (mxinv > 0) Then
     If (mxnode > 1) Then
        mxinv = Max(mxinv,mxnode*mtinv)
        mxinv = (3*(Nint(fdvar*Real(mxinv,wp))+mxnode-1))/mxnode
        mxinv = mxinv + (mxinv+4)/5 ! allow for 25% higher density
     End If
     mxpinv = 3
  Else
     mxpinv = 0
  End If



!!! GRIDDING PARAMETERS !!!

! Set grids for opted intramolecular distribution analysis if unset
! SO THEY ARE SWITCHES FOR EXISTENCE TOO

  If (mxgana > 0) Then
     If (mxgbnd1 == -1) Then
        If (mxgbnd > 0) Then
           mxgbnd1 = mxgbnd-4
        Else
           mxgbnd1 = Nint(rcbnd/delr_max)
        End If
     End If
     If (mxgang1 == -1) Then
        If (mxgang > 0) Then
           mxgang1 = mxgang-4
        Else
           mxgang1 = Nint(180.0_wp/delth_max)
        End If
     End If
     If (mxgdih1 == -1) Then
        If (mxgdih > 0) Then
           mxgdih1 = mxgdih-4
        Else
           mxgdih1 = Nint(360.0_wp/delth_max)
        End If
     End If
     If (mxginv1 == -1) Then
        If (mxginv > 0) Then
           mxginv1 = mxginv-4
        Else
           mxgdih1 = Nint(180.0_wp/delth_max)
        End If
     End If
     mxgana = Max(mxgbnd1,mxgang1,mxgdih1,mxginv1)
  End If
  mxtana = 0 ! initialise for buffer size purposes, set in read_field

! maximum number of rdf potentials (mxrdf = mxrdf)
! mxgrdf - maximum dimension of rdf and z-density arrays

  If ((.not. l_n_r) .or. lzdn) Then
     If (((.not. l_n_r) .and. mxrdf == 0) .and. (mxvdw > 0 .or. mxmet > 0)) &
        mxrdf = Max(mxvdw,mxmet) ! (vdw,met) == rdf scanning
     mxgrdf = Nint(rcut/rbin)
  Else
     mxgrdf = 0 ! RDF and Z-density function MUST NOT get called!!!
  End If

! RDFs particulars for USR (umbrella sampling restraints)

  If (l_usr) Then
     rusr   = 0.45_wp*width
     mxgusr = Nint(rusr/rbin)      ! allows for up to ~75% system volume shrinkage
     rusr   = Real(mxgusr,wp)*rbin ! round up and beautify for Andrey Brukhno's sake
  Else
     rusr   = 0.0_wp
     mxgusr = 0 ! decider on calling USR RDF
  End If

! maximum of all maximum numbers of grid points for all grids - used for mxbuff

  mxgrid = Max(mxgana,mxgvdw,mxgmet,mxgrdf,mxgusr,1004,Nint(rcut/delr_max)+4)

! grids setting and overrides

! maximum number of grid points for bonds

  mxgbnd = Merge(mxgbnd,Min(mxgbnd,Max(1004,Nint(rcbnd/delr_max)+4)),mxgbnd < 0)

! maximum number of grid points for angles

  mxgang = Merge(mxgang,Min(mxgang,Nint(180.0_wp/delth_max)+4),mxgang < 0)

! maximum number of grid points for dihedrals

  mxgdih = Merge(mxgdih,Min(mxgdih,Nint(360.0_wp/delth_max)+4),mxgdih < 0)

! maximum number of grid points for inversions

  mxginv = Merge(mxginv,Min(mxginv,Nint(180.0_wp/delth_max)+4),mxginv < 0)

! maximum number of grid points for electrostatics

  mxgele = Merge(-1,Max(1004,Nint(rcut/delr_max)+4),l_n_e)

! maximum number of grid points for vdw interactions - overwritten

  mxgvdw = Merge(-1,Max(1004,Nint(rvdw/delr_max)+4),l_n_v)

! maximum number of grid points for metal interactions

  mxgmet = Max(mxgmet,1004,Nint(rmet/delr_max)+4)

! maximum number of grid points for tersoff interaction arrays

  mxgter = Merge(-1,Max(1004,Nint(rcter/delr_max)+4),mxter <= 0)

! maximum of all maximum numbers of grid points for all grids - used for mxbuff

  mxgrid = Max(mxgrid,mxgbnd,mxgang,mxgdih,mxginv,mxgele,mxgvdw,mxgmet,mxgter)



!!! INTER-LIKE POTENTIAL PARAMETERS !!!

! maximum number of vdw potentials and parameters

  If (mxvdw > 0) Then
     mxvdw = mxvdw+1
     mxpvdw = 7
  Else
     mxpvdw = 0
  End If


! maximum number of metal potentials and parameters

  If (mxmet > 0) Then
     mxmet = mxmet+1
     mxpmet = 9
  Else
     mxpmet = 0
  End If


! maximum number of tersoff potentials (mxter = mxter) and parameters

  If (mxter > 0) Then
     If      (potter == 1) Then
        mxpter = 11
     Else If (potter == 2) Then
        mxpter = 16
     End If
  Else
     mxpter = 0
  End If


! maximum number of three-body potentials and parameters

  If (mxtbp > 0) Then
     mx2tbp = (mxatyp*(mxatyp+1))/2
     mxtbp  = mx2tbp*mxatyp
     If (rctbp < 1.0e-6_wp) rctbp=0.5_wp*rcut

     mxptbp = 5
  Else
     mx2tbp = 0
     mxtbp  = 0

     mxptbp = 0
  End If


! maximum number of four-body potentials and parameters

  If (mxfbp > 0) Then
     mx3fbp = (mxatyp*(mxatyp+1)*(mxatyp+2))/6
     mxfbp  = mx3fbp*mxatyp
     If (rcfbp < 1.0e-6_wp) rcfbp=0.5_wp*rcut

     mxpfbp = 3
  Else
     mx3fbp = 0
     mxfbp  = 0

     mxpfbp = 0
  End If


! maximum number of external field parameters

  If (lext) Then
     mxpfld = 6
  Else
     mxpfld = 0
  End If



! DD PARAMETERS - by hypercube mapping of MD cell onto machine resources
! Dependences: MD cell widths (explicit) and machine resources (implicit)

  Call map_domains(imc_n,celprp(7),celprp(8),celprp(9))

  If (idnode == 0) Write(nrite,'(/,/,1x,a,3i6)') 'node/domain decomposition (x,y,z): ', nprx,npry,nprz



! TTM matters

  rpad2=0.0_wp
  If (l_ttm) Then

! two-temperature model: determine number of CITs
! in x- and y-directions based on number in z-direction
! and system size

     delz     = cell(9)/Real(ntsys(3),wp)
     ntsys(1) = Nint(cell(1)/delz)
     ntsys(2) = Nint(cell(5)/delz)
     delx     = cell(1)/Real(ntsys(1),wp)
     dely     = cell(5)/Real(ntsys(2),wp)
     volume   = delx*dely*delz
     rvolume  = 1.0_wp/volume

! Check number of electronic temperature cells is
! >= to number of ionic temperature cells

     If (Any(eltsys < ntsys)) Call error(670)

! Check max rpad does not go too far for determining ionic temperatures

     tol=Min(delx,dely,delz)
     Do i=1,nprx-1
        test=Real(i,wp)*cell(1)*r_nprx
        test=test-delx*Floor(test/delx)
        If (test > zero_plus) tol=Min(tol,test)
     End Do
     Do i=1,npry-1
        test=Real(i,wp)*cell(5)*r_npry
        test=test-dely*Floor(test/dely)
        If (test > zero_plus) tol=Min(tol,test)
     End Do
     Do i=1,nprz-1
        test=Real(i,wp)*cell(9)*r_nprz
        test=test-delz*Floor(test/delz)
        If (test > zero_plus) tol=Min(tol,test)
     End Do
     rpad2=tol

  End If

! If redistribute option selected, check for sufficient electronic temperature
! cells to redistribute energy when ionic temperature cells are switched off:
! if not available, switch off this option

  If (redistribute .and. (eltsys(1) < ntsys(1)+2 .or. &
                          eltsys(2) < ntsys(2)+2 .or. &
                          eltsys(3)<ntsys(3)+2)) Then
    Call warning(500,0.0_wp,0.0_wp,0.0_wp)
    redistribute = .false.
  End If

! Calculate average atomic density: if not overridden by
! 'ttm atomdens' directive in CONTROL file, will be used
! to convert specific heat capacities to volumetric
! heat capacity etc.

  sysrho = Real(megatm,Kind=wp)/(cell(1)*cell(5)*cell(9))



! LC and VNL matters

5 Continue

  If (rpad > zero_plus) Then

! define cut

     cut=rcut+1.0e-6_wp

! Provide advise on decomposition

     qlx=Int(celprp(7)/cut)
     qly=Int(celprp(8)/cut)
     qlz=Int(celprp(9)/cut)

     If (idnode == 0) Write(nrite,'(/,1x,a,i6,a,3(i0,a))')                           &
        'pure cutoff driven limit on largest possible decomposition:', qlx*qly*qlz , &
        ' nodes/domains (', qlx,',',qly,',',qlz,')'

     qlx=Max(1,qlx/2)
     qly=Max(1,qly/2)
     qlz=Max(1,qlz/2)

     If (idnode == 0) Write(nrite,'(/,1x,a,i6,a,3(i0,a))')                           &
        'pure cutoff driven limit on largest balanced decomposition:', qlx*qly*qlz , &
        ' nodes/domains (', qlx,',',qly,',',qlz,')'

  Else

     If (lrpad) lrpad0=.true.

  End If

10 Continue ! possible rcut redefinition...

! Define link-cell cutoff (minimum width)

  rlnk = rcut + rpad

! define cut

  cut=rlnk+1.0e-6_wp

! Provide advise on decomposition

  qlx=Int(celprp(7)/cut)
  qly=Int(celprp(8)/cut)
  qlz=Int(celprp(9)/cut)

  If (idnode == 0) Write(nrite,'(/,1x,a,i6,a,3(i0,a))')                       &
     'cutoffs driven limit on largest possible decomposition:', qlx*qly*qlz , &
     ' nodes/domains (', qlx,',',qly,',',qlz,')'

  qlx=Max(1,qlx/2)
  qly=Max(1,qly/2)
  qlz=Max(1,qlz/2)

  If (idnode == 0) Write(nrite,'(/,1x,a,i6,a,3(i0,a))')                       &
     'cutoffs driven limit on largest balanced decomposition:', qlx*qly*qlz , &
     ' nodes/domains (', qlx,',',qly,',',qlz,')'

! calculate link cell dimensions per node

  ilx=Int(r_nprx*celprp(7)/cut)
  ily=Int(r_npry*celprp(8)/cut)
  ilz=Int(r_nprz*celprp(9)/cut)

! print link cell algorithm and check for violations or...

  If (idnode == 0) Write(nrite,'(/,1x,a,3i6)') "link-cell decomposition 1 (x,y,z): ",ilx,ily,ilz

  tol=Min(0.05_wp,0.005_wp*rcut)                                        ! tolerance
  test = 0.02_wp * Merge( 1.0_wp, 2.0_wp, mxspl > 0)                    ! 2% (w/ SPME/PS) or 4% (w/o SPME/PS)
  cut=Min(r_nprx*celprp(7),r_npry*celprp(8),r_nprz*celprp(9))-1.0e-6_wp ! domain size

  If (ilx*ily*ilz == 0) Then
     If (l_trm) Then ! we are prepared to exit gracefully(-:
        rcut = cut   ! - rpad (was zeroed in scan_control)
        If (idnode == 0) Write(nrite,'(/,1x,a)') &
           "*** warning - real space cutoff reset has occurred, early run termination is due !!! ***"
        Go To 10
     Else
        If (cut < rcut) Then
           If (idnode == 0) Write(nrite,*) '*** warning - rcut <= Min(domain width) < rlnk = rcut + rpad !!! ***'
           Call error(307)
        Else ! rpad is defined & in 'no strict' mode
           If (rpad > zero_plus .and. (.not.l_str)) Then ! Re-set rpad with some slack
              rpad = Min( 0.95_wp * (cut - rcut) , test * rcut)
              If (rpad2 > zero_plus) rpad=Min(rpad,rpad2)
              rpad = Real( Int( 100.0_wp * rpad ) , wp ) / 100.0_wp
              If (rpad < tol) rpad = 0.0_wp ! Don't bother
              Go To 10
           Else
              If (idnode == 0) Write(nrite,*) '*** warning - rcut <= Min(domain width) < rlnk = rcut + rpad !!! ***'
              Call error(307)
           End If
        End If
     End If
  Else ! push/reset the limits in 'no strict' mode and be somewhat hopeful for undefined rpad in 'strict' when rpad=0
     If (.not.(mxmet == 0 .and. l_n_e .and. l_n_v .and. mxrdf == 0 .and. kim == ' ')) Then ! 2b link-cells are needed
        rpad1=0.0_wp
        If (mxnode == 1 .and. Min(ilx,ily,ilz) < 2) Then ! catch & handle exception
           rpad1 = 0.95_wp * (0.5_wp*width - rcut - 1.0e-6_wp)
        End If

        If (rpad <= zero_plus) Then ! When rpad is undefined give it some value
           If (Int(Real(Min(ilx,ily,ilz),wp)/(1.0_wp+test)) >= 2) Then ! good non-exception
              rpad = test * rcut
              If (rpad1 > zero_plus) rpad=Min(rpad,rpad1)
              If (rpad2 > zero_plus) rpad=Min(rpad,rpad2)
              rpad = Real( Int( 100.0_wp * rpad ) , wp ) / 100.0_wp
              If (rpad > tol) Go To 10
           Else ! not so good non-exception
              rpad = Min( 0.95_wp * ( Min ( r_nprx * celprp(7) / Real(ilx,wp) , &
                                            r_npry * celprp(8) / Real(ily,wp) , &
                                            r_nprz * celprp(9) / Real(ilz,wp) ) &
                                      - rcut - 1.0e-6_wp ) , test * rcut )
              If (rpad1 > zero_plus) rpad=Min(rpad,rpad1)
              If (rpad2 > zero_plus) rpad=Min(rpad,rpad2)
              rpad = Real( Int( 100.0_wp * rpad ) , wp ) / 100.0_wp ! round up
           End If
        End If

        If (rpad > zero_plus) Then
           If (rpad < tol) rpad = 0.0_wp ! Don't bother
        End If
     Else
        If (rpad >= zero_plus) rpad = 0.0_wp ! Don't bother
     End If

     If (l_str) Then
        If (lrpad0) Then ! forget about it
           rpad = 0.0_wp
        Else ! less hopeful for undefined rpad
           rpad = 0.95_wp * rpad
           rpad = Real( Int( 100.0_wp * rpad ) , wp ) / 100.0_wp ! round up
           If (rpad > zero_plus) Then
              If (rpad < tol) rpad = 0.0_wp ! Don't bother
           End If
        End If
     Else ! 'no strict'
        If (lrpad0) Then ! forget about it
           rpad = 0.0_wp
        End If
     End If

     rlnk = rcut + rpad ! recalculate rlnk respectively
  End If
  llvnl = (rpad > zero_plus) ! Determine/Detect conditional VNL updating at start

  If (ilx < 3 .or. ily < 3 .or. ilz < 3) Call warning(100,0.0_wp,0.0_wp,0.0_wp)

! get total link cells per domain (boundary padding included)
! total link-cells per node/domain is ncells = (ilx+2)*(ily+2)*(ilz+2)
! allow for more (possible b-spline SPME triggered increase in nlast),
! total link-cells per node/domain is ncells = (ilx+3)*(ily+3)*(ilz+3)
! allow for thermal expansion of unsettled systems
! total link-cells per node/domain is ncells = (ilx+4)*(ily+4)*(ilz+4)

  If      (imcon == 0                ) Then
     mxcell = Nint((fdvar**4) * Real((ilx+4)*(ily+4)*(ilz+4),wp))
  Else If (imcon == 6 .or. imc_n == 6) Then
     mxcell = Nint((fdvar**3) * Real((ilx+4)*(ily+4)*(ilz+4),wp))
  Else
     mxcell = Nint((fdvar**2) * Real((ilx+4)*(ily+4)*(ilz+4),wp))
  End If



! SPME electrostatics particularities

! qlx,qly,qlz - SPME fictional link-cell dimensions postulating that:
! nprx <= kmaxa/mxspl, npry <= kmaxb/mxspl, nprz <= kmaxc/mxspl.
! Otherwise, this node's b-splines in SPME will need 'positive halo'
! that is not on the immediate neighbouring nodes in negative
! directions but beyond them (which may mean self-halo in some cases)

  kmaxa = kmaxa1
  kmaxb = kmaxb1
  kmaxc = kmaxc1

  qlx = ilx
  qly = ily
  qlz = ilz

! mxspl = 0 is an indicator for no SPME or Poisson Solver electrostatics in CONTROL

  If (mxspl /= 0) Then

! ensure (kmaxa,kmaxb,kmaxc) consistency between the DD
! processor grid (map_domains is already called) and the grid
! method or comment out adjustments if using ewald_spme_force~

     Call adjust_kmax( kmaxa, nprx )
     Call adjust_kmax( kmaxb, npry )
     Call adjust_kmax( kmaxc, nprz )

! Calculate and check ql.

     qlx = Min(qlx , kmaxa/(mxspl*nprx))
     qly = Min(qly , kmaxb/(mxspl*npry))
     qlz = Min(qlz , kmaxc/(mxspl*nprz))

     If (.not.llvnl) Then
        mxspl1=mxspl
     Else
        mxspl1=mxspl+Ceiling((rpad*Real(mxspl,wp))/rcut)

! Redifine ql.

        qlx = Min(ilx , kmaxa/(mxspl1*nprx))
        qly = Min(ily , kmaxb/(mxspl1*npry))
        qlz = Min(ilz , kmaxc/(mxspl1*nprz))
     End If

! Hard luck, giving up after trying once more

     If (qlx*qly*qlz == 0) Then
        If ((lrpad0 .eqv. lrpad) .and. rpad > zero_plus) Then ! defaulted padding must be removed
           rpad=0.0_wp
           lrpad0=.true.
           Go To 5
        Else
           test = Min( Real(kmaxa1,wp)/Real(nprx,wp), &
                       Real(kmaxb1,wp)/Real(npry,wp), &
                       Real(kmaxc1,wp)/Real(nprz,wp) ) / Real(mxspl,wp)
           tol  = Min( Real(kmaxa,wp)/Real(nprx,wp), &
                       Real(kmaxb,wp)/Real(npry,wp), &
                       Real(kmaxc,wp)/Real(nprz,wp) ) / Real(mxspl1,wp)
           If (idnode == 0) Then
              Write(nrite,'(/,1x,a,i6,a,3(i0,a))')                               &
                 'SPME driven limit on largest possible domain decomposition: ', &
                 (kmaxa/mxspl1)*(kmaxb/mxspl1)*(kmaxc/mxspl1),                   &
                 ' nodes/domains (',                                             &
                 kmaxa/mxspl1,',',kmaxb/mxspl1,',',kmaxc/mxspl1,                 &
                 ') for currently specified cutoff (with padding) & Ewald precision (sum parameters)'
              Write(nrite,'(2(/,1x,a,f6.2,a))')                                          &
                 'SPME suggested factor to decrease currently specified cutoff (with padding) by: ', &
                 1.0_wp/test, ' for currently speccified Ewald precision & domain decomposition',    &
                 'SPME suggested factor to increase current Ewald precision by: ',                   &
                 (1.0_wp-tol)*100.0_wp, ' for currently specified cutoff (with padding) & domain decomposition'
           End If
        End If
        Call error(308)
     End If

  End If



! decide on MXATMS while reading CONFIG and scan particle density

  Call read_config(megatm,levcfg,l_ind,l_str,rcut,dvar,xhi,yhi,zhi,dens0,dens)

! Create f(fdvar,dens0,dens)

  If (mxnode == 1 .or. (imcon == 0 .or. imcon == 6 .or. imc_n == 6)) Then
     fdens = fdvar * (0.65_wp*dens0 + 0.35_wp*dens)
  Else If (Min(ilx,ily,ilz) == 1) Then
     fdens = fdvar * (0.50_wp*dens0 + 0.50_wp*dens)
  Else
     fdens = fdvar * (0.35_wp*dens0 + 0.65_wp*dens)
  End If

! Get reasonable - all particles in one link-cell

  tol   = Real(megatm,wp) / (Real(ilx*ily*ilz,wp) * Real(mxnode,wp))
  fdens = Min(fdens,tol)

! density variation affects the link-cell arrays' dimension
! more than domains(+halo) arrays' dimensions, in case of
! events of extreme collapse in atomic systems (aggregation)

! mxlist is the maximum length of link-cell list (dens * 4/3 pi rlnk^3)
! + 75% extra tolerance - i.e f(dens0,dens)*(7.5/3)*pi*rlnk^3

  mxlist = Nint(fdens*2.5_wp*pi*rlnk**3)
  mxlist = Min(mxlist,megatm-1) ! mxexcl

  If (mxlist < mxexcl-1) Then
     Call warning(6,Real(mxlist,wp),Real(mxexcl,wp),0.0_wp)
     mxlist=mxexcl-1
  End If

! get link-cell volume

  vcell = volm / (Real(ilx*ily*ilz,wp) * Real(mxnode,wp))

! get averaged link-cell particle number, boosted by fdens
! + 25% extra tolerance

  test = fdens*vcell * 1.25_wp

! set dimension of working coordinate arrays

  mxatms = Max(1 , Nint(test * Real((ilx+3)*(ily+3)*(ilz+3),wp)))
  If       (mxnode == 1 .or. imcon == 0) Then ! ilx >= 2 && ily >= 2 && ilz >= 2
     mxatms = Nint(Min(Real(mxatms,wp),9.0_wp*fdvar*Real(megatm,wp)))
  Else If (imcon == 6 .or. imc_n == 6)   Then ! mxnode >= 4 .or. (ilx >= 2 && ily >= 2)
     mxatms = Nint(Min(Real(mxatms,wp),6.0_wp*fdvar*Real(megatm,wp)))
  Else If (ilx*ily*ilz < 2) Then ! mxnode >= 8
     mxatms = Nint(Min(Real(mxatms,wp),6.0_wp*fdvar*Real(megatm,wp)))
  End If

! maximum number of particles per domain (no halo)

  mxatdm = Max(1 , Nint(test * Real((ilx+1)*(ily+1)*(ilz+1),wp)))
  mxatdm = Min(mxatdm,megatm)

! maximum number of timesteps in stack arrays

  mxstak = Max(100,mxstak)

! maximum number of variables in stack arrays
! update number if the MSD option is used, 51=1+27+...+9+9+1+2+2
! consult statistic_collect for more information

  mxnstk = 51 + mxatyp + Merge(2*mxatdm,0,l_msd)

! maximum dimension of principal transfer buffer

! deport_atomic_data & export_atomic_data (+ metal_ld_export.f90
! defects_reference_export & statistics_connect_deport if used)
! are supposed to be the most MPI/buffer consuming routines

! deporting total per atom

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) / Real(ilx*ily*ilz,wp)
  dens0 = dens0/Max(rlnk/0.2_wp,1.0_wp)
  mxbfdp = Merge( 2, 0, mxnode > 1) * Nint( Real(                          &
           mxatdm*(18+12 + Merge(3,0,llvnl) + (mxexcl+1)                 + &
           Merge(mxexcl+1 + Merge(mxexcl+1,0,keyind == 1),0,mximpl > 0)  + &
           Merge(2*(6+mxstak), 0, l_msd)) + 3*vafsamp                    + &
           4*mxshl+4*mxcons+(Sum(mxtpmf(1:2)+3))*mxpmf+(mxlrgd+13)*mxrgd + &
           3*mxteth+4*mxbond+5*mxangl+8*mxdihd+6*mxinv,wp) * dens0)

! statistics connect deporting total per atom

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) / Real(ilx*ily*ilz,wp)
  dens0 = dens0/Max(rlnk/0.2_wp,1.0_wp)
  mxbfss = Merge( 4, 0, mxnode > 1) * Nint( Real(mxatdm*(8 + Merge(2*(6+mxstak), 0, l_msd)),wp) * dens0)

! exporting single per atom (times 13 up to 35)

!  dens  = Real(((ilx+3)*(ily+3)*(ilz+3))/Min(ilx,ily,ilz)+3,wp) / Real(ilx*ily*ilz,wp)
  dens  = Real(((qlx+2)*(qly+2)*(qlz+2))/Min(qlx,qly,qlz)+2,wp) / Real(qlx*qly*qlz,wp)
  mxbfxp = 2 * Nint(Real(mxatdm,wp) * dens) ! included induced dipoles

! shared units single per atom of all shared unit

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) - 1.0_wp
  dens0 = dens0/Max(rlnk/2.0_wp,1.0_wp)
  mxbfsh = Merge( 1, 0, mxnode > 1) * Nint(Real(Max(2*mxshl,2*mxcons,mxlrgd*mxrgd),wp) * dens0)

  mxbuff = Max( mxbfdp , 35*mxbfxp , 4*mxbfsh , 2*(kmaxa/nprx)*(kmaxb/npry)*(kmaxc/nprz)+10 , &
                mxnstk*mxstak , mxgrid , mxgrdf , mxlrgd*Max(mxrgd,mxtrgd), mxtrgd*(4+3*mxlrgd), 10000 )

! reset (increase) link-cell maximum (mxcell)
! if tersoff or three- or four-body potentials exist

  If (mxter > 0 .or. mxtbp > 0 .or. mxfbp > 0) Then
     cut=rcut+1.0e-6_wp ! reduce cut
     If (mxter > 0) cut = Min(cut,rcter+1.0e-6_wp)
     If (mxtbp > 0) cut = Min(cut,rctbp+1.0e-6_wp)
     If (mxfbp > 0) cut = Min(cut,rcfbp+1.0e-6_wp)

     ilx=Int(r_nprx*celprp(7)/cut)
     ily=Int(r_npry*celprp(8)/cut)
     ilz=Int(r_nprz*celprp(9)/cut)

     If (idnode == 0) Write(nrite,'(/,1x,a,3i6)') "link-cell decomposition 2 (x,y,z): ",ilx,ily,ilz

     If (ilx < 3 .or. ily < 3 .or. ilz < 3) Call error(305)

     If      (imcon == 0                ) Then
        mxcell = Max(mxcell,Nint((fdvar**4) * Real((ilx+5)*(ily+5)*(ilz+5),wp)))
     Else If (imcon == 6 .or. imc_n == 6) Then
        mxcell = Max(mxcell,Nint((fdvar**3) * Real((ilx+5)*(ily+5)*(ilz+5),wp)))
     Else
        mxcell = Max(mxcell,Nint((fdvar**2) * Real((ilx+5)*(ily+5)*(ilz+5),wp)))
     End If
  End If

End Subroutine set_bounds
