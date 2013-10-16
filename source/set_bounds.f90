Subroutine set_bounds                                  &
           (levcfg,imcon,l_vv,l_str,l_n_e,l_n_v,l_ind, &
           rcut,rvdw,rmet,rbin,nstfce,alpha,width)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to determine various limits: array bounds,
! iteration and others as specified in setup_module
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode
  Use setup_module
  Use domains_module,     Only : map_domains,nprx,npry,nprz
  Use config_module,      Only : cfgname,imc_n,cell,volm
  Use msd_module
  Use tersoff_module,     Only : potter
  Use development_module, Only : l_trm

  Implicit None

  Logical,           Intent(   Out ) :: l_vv,l_str,l_n_e,l_n_v,l_ind
  Integer,           Intent(   Out ) :: levcfg,imcon,nstfce
  Real( Kind = wp ), Intent(   Out ) :: rvdw,rmet,rbin,alpha,width
  Real( Kind = wp ), Intent(   Out ) :: rcut

  Logical           :: l_n_r,lzdn,lter,ltbp,lfbp,lext
  Integer           :: ilx,ily,ilz,qlx,qly,qlz,megatm, &
                       mtshl,mtcons,mtrgd,mtteth,      &
                       mtbond,mtangl,mtdihd,mtinv
  Real( Kind = wp ) :: ats,celprp(1:10),cut,         &
                       dens0,dens,dvar,fdens,        &
                       test,vcell,rcter,rctbp,rcfbp, &
                       sidex,sidey,sidez,xhi,yhi,zhi

! define zero+ and half+/- (setup_module)

  zero_plus  = Nearest( 0.0_wp , +1.0_wp)
  half_plus  = Nearest( 0.5_wp , +1.0_wp)
  half_minus = Nearest( 0.5_wp , -1.0_wp)

! scan the FIELD file data

  Call scan_field                                     &
           (l_n_e,mxsite,mxatyp,megatm,mxtmls,mxexcl, &
           mtshl,mxtshl,mxshl,mxfshl,                 &
           mtcons,mxtcon,mxcons,mxfcon,               &
           mxtpmf,mxpmf,mxfpmf,                       &
           mtrgd,mxtrgd,mxrgd,mxlrgd,mxfrgd,          &
           mtteth,mxtteth,mxteth,mxftet,              &
           mtbond,mxtbnd,mxbond,mxfbnd,               &
           mtangl,mxtang,mxangl,mxfang,               &
           mtdihd,mxtdih,mxdihd,mxfdih,               &
           mtinv,mxtinv,mxinv,mxfinv,                 &
           mxrdf,mxgrid,mxvdw,rvdw,                   &
           mxmet,mxmed,mxmds,rmet,                    &
           mxter,rcter,mxtbp,rctbp,mxfbp,rcfbp,lext)

! scan CONFIG file data

  Call scan_config(megatm,imc_n,cfgname,levcfg,imcon,cell,xhi,yhi,zhi)

! halt execution for unsupported image conditions in DD
! checks for some inherited from DL_POLY_2 are though kept

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(514)

! scan CONTROL file data

  Call scan_control                                  &
           (mxrdf,mxvdw,rvdw,mxmet,rmet,mxter,rcter, &
           imcon,imc_n,cell,xhi,yhi,zhi,             &
           l_str,l_vv,l_n_e,l_n_r,lzdn,l_n_v,l_ind,  &
           dvar,rcut,rbin,mxstak,                    &
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
! (for use in link-cells) and volume

  Call dcell(cell,celprp)

  volm = celprp(10)

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) volm=0.5_wp*volm

! check value of cutoff and reset if necessary

  If (imcon > 0) Then
     width=Min(celprp(7),celprp(8),celprp(9))
     If (imcon == 4) width=rt3*cell(1)/2.0_wp
     If (imcon == 5) width=cell(1)
     If (imcon == 6) width=Min(celprp(7),celprp(8))

! halt program if potential cutoff exceeds the minimum half-cell width

     If (rcut > width/2.0_wp) Then
        Call warning(3,rcut,width/2.0_wp,0.0_wp)
        Call error(95)
     End If
  End If



!!! INTRA-LIKE POTENTIAL PARAMETERS !!!

! maximum number of core-shell units per node

  If (mxshl > 0 .and. mxnode > 1) Then
     mxshl = Max(mxshl,mxnode*mtshl)
     mxshl = (3*(Nint((dvar**1.7_wp)*Real(mxshl,wp))+mxnode-1))/mxnode
  End If


! maximum number of constraints per node

  If (mxcons > 0 .and. mxnode > 1) Then
     mxcons = Max(mxcons,mxnode*mtcons)
     mxcons = (3*(Nint((dvar**1.7_wp)*Real(mxcons,wp))+mxnode-1))/mxnode
  End If


! maximum number of PMF constraints per MD cell - mxpmf
! (only one type of PMF in only one type of molecule !!!)
! maximum number of atoms per PMF unit - mxtpmf(1:2) (only two units per pmf)


! maximum number of RBs per node

  If (mxrgd > 0 .and. mxnode > 1) Then
     mxrgd = Max(mxrgd,mxnode*mtrgd)
     mxrgd = (3*(Nint((dvar**1.7_wp)*Real(mxrgd,wp))+mxnode-1))/mxnode
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
        mxteth = (3*(Nint((dvar**1.7_wp)*Real(mxteth,wp))+mxnode-1))/mxnode
     End If
     mxpteth = 3
  Else
     mxpteth = 0
  End If


! maximum number of chemical bonds per node and bond potential parameters

  If (mxbond > 0) Then
     If (mxnode > 1) Then
        mxbond = Max(mxbond,mxnode*mtbond)
        mxbond = (3*(Nint((dvar**1.7_wp)*Real(mxbond,wp))+mxnode-1))/mxnode
     End If
     mxpbnd = 4
  Else
     mxpbnd = 0
  End If


! maximum number of bond angles per node and angular potential parameters

  If (mxangl > 0) Then
     If (mxnode > 1) Then
        mxangl = Max(mxangl,mxnode*mtangl)
        mxangl = (3*(Nint((dvar**1.7_wp)*Real(mxangl,wp))+mxnode-1))/mxnode
     End If
     mxpang = 6
  Else
     mxpang = 0
  End If


! maximum number of torsion angles per node and dihedral potential parameters

  If (mxdihd > 0) Then
     If (mxnode > 1) Then
        mxdihd = Max(mxdihd,mxnode*mtdihd)
        mxdihd = (3*(Nint((dvar**1.7_wp)*Real(mxdihd,wp))+mxnode-1))/mxnode
     End If
     mxpdih = 7
  Else
     mxpdih = 0
  End If


! maximum number of inversions per node and inversion potential parameters

  If (mxinv > 0) Then
     If (mxnode > 1) Then
        mxinv = Max(mxinv,mxnode*mtinv)
        mxinv = (3*(Nint((dvar**1.7_wp)*Real(mxinv,wp))+mxnode-1))/mxnode
     End If
     mxpinv = 3
  Else
     mxpinv = 0
  End If



!!! INTER-LIKE POTENTIAL PARAMETERS !!!

! maximum number of grid points in potentials arrays

  mxgrid = Max(mxgrid,1000,Int(rcut/0.01_wp+0.5_wp)+4)


! maximum number of rdf potentials (mxrdf = mxrdf)
! mxgrdf - maximum dimension of rdf and z-density arrays

  If ((.not. l_n_r) .or. lzdn) Then
     If (((.not.l_n_r) .and. mxrdf == 0) .and. (mxvdw > 0 .or. mxmet > 0)) &
        mxrdf = Max(mxvdw,mxmet) ! (vdw,met) == rdf scanning
     mxgrdf = Nint(rcut/rbin)
  Else
     mxgrdf = 0 ! RDF and Z-density function MUST NOT get called!!!
  End If


! maximum number of vdw potentials and parameters

  If (mxvdw > 0) Then
     mxvdw = mxvdw+1
     mxpvdw = 5
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
     lter = .true.

     If      (potter == 1) Then
        mxpter = 11
     Else If (potter == 2) Then
        mxpter = 16
     End If
  Else
     lter = .false.

     mxpter = 0
  End If


! maximum number of three-body potentials and parameters

  If (mxtbp > 0) Then
     ltbp   = .true.
     mx2tbp = (mxatyp*(mxatyp+1))/2
     mxtbp  = mx2tbp*mxatyp
     If (rctbp < 1.0e-6_wp) rctbp=0.5_wp*rcut

     mxptbp = 5
  Else
     ltbp   = .false.
     mx2tbp = 0
     mxtbp  = 0

     mxptbp = 0
  End If


! maximum number of four-body potentials and parameters

  If (mxfbp > 0) Then
     lfbp   = .true.
     mx3fbp = (mxatyp*(mxatyp+1)*(mxatyp+2))/6
     mxfbp  = mx3fbp*mxatyp
     If (rcfbp < 1.0e-6_wp) rcfbp=0.5_wp*rcut

     mxpfbp = 3
  Else
     lfbp   = .false.
     mx3fbp = 0
     mxfbp  = 0

     mxpfbp = 0
  End If


! maximum number of external field parameters

  If (lext) Then
     mxpfld = 4
  Else
     mxpfld = 0
  End If



! DD PARAMETERS - by hypercube mapping of MD cell onto machine resources
! Dependences: MD cell widths (explicit) and machine resources (implicit)

  Call map_domains(imc_n,celprp(7),celprp(8),celprp(9))

  If (idnode == 0) Write(nrite,'(/,/,1x,a,3i6)') 'node/domain decomposition (x,y,z): ', nprx,npry,nprz



10 Continue ! possible rcut redefinition if l_trm=.true.

! define cut

  cut=rcut+1.0e-6_wp

! calculate link cell dimensions per node

  sidex=1.0_wp/Real(nprx,wp)
  sidey=1.0_wp/Real(npry,wp)
  sidez=1.0_wp/Real(nprz,wp)

  ilx=Int(sidex*celprp(7)/cut)
  ily=Int(sidey*celprp(8)/cut)
  ilz=Int(sidez*celprp(9)/cut)

! print link cell algorithm and check for violations

  If (idnode == 0) Write(nrite,'(/,1x,a,3i6)') "link-cell decomposition 1 (x,y,z): ",ilx,ily,ilz

  If (ilx*ily*ilz == 0) Then
     If (.not.l_trm) Then
        Call error(307)
     Else ! we are prepared to exit gracefully(-:
        rcut=Min(sidex*celprp(7),sidey*celprp(8),sidez*celprp(9))-1.0e-6_wp
        If (idnode == 0) Write(nrite,'(/,1x,a)') &
           "*** warning - real space cutoff reset has occured, early run termination is due !!! ***"
        Go To 10
     End If
  End If

  If (ilx < 4 .or. ily < 4 .or. ilz < 4) Call warning(100,0.0_wp,0.0_wp,0.0_wp)

! get total link cells per domain (boundary padding included)
! total link-cells per node/domain is ncells = (ilx+2)*(ily+2)*(ilz+2)
! allow for more (possible b-spline SPME triggered increase in nlast),
! total link-cells per node/domain is ncells = (ilx+3)*(ily+3)*(ilz+3)
! allow for thermal expansion of unsettled systems
! total link-cells per node/domain is ncells = (ilx+4)*(ily+4)*(ilz+4)

  mxcell = (ilx+4)*(ily+4)*(ilz+4)



! SPME electrostatics particularities

! qlx,qly,qlz - SPME fictional link-cell dimensions postulating that:
! nprx <= kmaxa/mxspl, npry <= kmaxb/mxspl, nprz <= kmaxc/mxspl.
! Otherwise, this node's b-splines in SPME will need 'positive halo'
! that is not on the immediate neighbouring nodes in negative
! directions but beyond them (which may mean self-halo in some cases)

  kmaxa = kmaxa1
  kmaxb = kmaxb1
  kmaxc = kmaxc1

! mxspl = 0 is an indicator for no SPME electrostatics in CONTROL

  If (mxspl == 0) Then

     qlx = ilx
     qly = ily
     qlz = ilz

  Else

! ensure (kmaxa,kmaxb,kmaxc) consistency with what DD
! (map_domains is already called) and DaFT are capable of
! or comment out adjustments if using ewald_spme_force~

     Call adjust_kmax( kmaxa, nprx )
     Call adjust_kmax( kmaxb, npry )
     Call adjust_kmax( kmaxc, nprz )

! Calculate and check ql.

     qlx = Min(ilx , kmaxa/(mxspl*nprx))
     qly = Min(ily , kmaxb/(mxspl*npry))
     qlz = Min(ilz , kmaxc/(mxspl*nprz))

     If (qlx*qly*qlz == 0) Call error(308)

  End If



! decide on MXATMS while reading CONFIG and scan particle density

  Call read_config &
           (megatm,levcfg,imcon,l_ind,l_str,rcut,dvar,xhi,yhi,zhi,dens0,dens)

! Create f(dvar,dens0,dens)

  If (mxnode == 1 .or. (imcon == 0 .or. imcon == 6 .or. imc_n == 6)) Then
     fdens = (dvar**1.7_wp) * (0.65_wp*dens0 + 0.35_wp*dens)
  Else If (Min(ilx,ily,ilz) == 1) Then
     fdens = (dvar**1.7_wp) * (0.50_wp*dens0 + 0.50_wp*dens)
  Else
     fdens = (dvar**1.7_wp) * (0.35_wp*dens0 + 0.65_wp*dens)
  End If

! density variation affects the link-cell arrays' dimension
! more than domains(+halo) arrays' dimensions, in case of
! events of extreme collapse in atomic systems (aggregation)

! mxlist is the maximum length of link-cell list (dens * 4/3 pi rcut^3)
! + 75% extra tolerance - i.e f(dens0,dens)*(7.5/3)*pi*rcut^3

  mxlist = Nint(fdens*2.5_wp*pi*rcut**3)
  mxlist = Min(mxlist,megatm-1) ! mxexcl

  If (mxlist+1 < mxexcl) Then
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
  If (mxnode == 1 .or. (imcon == 0 .or. imcon == 6 .or. imc_n == 6)) Then
    mxatms = Nint(Min(Real(mxatms,wp),Real(27.0_wp,wp)*Real(megatm,wp)))
  Else
    mxatms = Nint(Min(Real(mxatms,wp),Real(13.5_wp,wp)*Real(megatm,wp)))
  End If

! maximum number of particles per domain (no halo)

  mxatdm = Max(1 , Nint(test * Real((ilx+1)*(ily+1)*(ilz+1),wp)))
  mxatdm = Min(mxatdm,megatm)

! maximum number of timesteps in stack arrays

  mxstak = Max(100,mxstak)

! maximum number of variables in stack arrays

  mxnstk = 50+mxatyp

! update maximum number of variables in stack arrays if MSD option is used

  If (l_msd) mxnstk = mxnstk+2*mxatdm

! maximum dimension of principal transfer buffer

! deport_atomic_data & export_atomic_data (and
! and defects_reference_export & metal_ld_export.f90 if used)
! are supposed to be the most MPI/buffer consuming routines

! deporting total per atom

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) / Real(ilx*ily*ilz,wp)
  dens0 = dens0/Max(rcut/0.2_wp,1.0_wp)
  mxbfdp = Merge( 2, 0, mxnode > 1) * Nint( Real(                          &
          mxatdm*(18+12+mxexcl + Merge(2*6+mxstak, 0, l_msd))           + &
          4*mxshl+4*mxcons+(Sum(mxtpmf(1:2)+3))*mxpmf+(mxlrgd+13)*mxrgd + &
          3*mxteth+4*mxbond+5*mxangl+8*mxdihd+6*mxinv,wp) * dens0)

! exporting single per atom

  dens  = Real(((qlx+2)*(qly+2)*(qlz+2))/Min(qlx,qly,qlz)+2,wp) / Real(qlx*qly*qlz,wp)
  mxbfxp = 2 * Nint(Real(mxatdm,wp) * dens)

! shared units single per atom

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) - 1.0_wp
  dens0 = dens0/Max(rcut/2.0_wp,1.0_wp)
  mxbfsh = Merge( 1, 0, mxnode > 1) * Nint(Real(Max(2*mxshl,2*mxcons,mxlrgd*mxrgd),wp) * dens0)

  mxbuff = Max( mxbfdp , 13*mxbfxp , 4*mxbfsh , 2*(kmaxa/nprx)*(kmaxb/npry)*(kmaxc/nprz)+10 , &
                mxnstk*mxstak , mxgrid , mxgrdf , mxlrgd*Max(mxrgd,mxtrgd), mxtrgd*(4+3*mxlrgd), 10000 )

! reset (increase) link-cell maximum (mxcell)
! if tersoff or three- or four-body potentials exist

  If (lter .or. ltbp .or. lfbp) Then
     If (lter) cut = Min(cut,rcter+1.0e-6_wp)
     If (ltbp) cut = Min(cut,rctbp+1.0e-6_wp)
     If (lfbp) cut = Min(cut,rcfbp+1.0e-6_wp)

     ilx=Int(sidex*celprp(7)/cut)
     ily=Int(sidey*celprp(8)/cut)
     ilz=Int(sidez*celprp(9)/cut)

     If (idnode == 0) Write(nrite,'(/,1x,a,3i6)') "link-cell decomposition 2 (x,y,z): ",ilx,ily,ilz

     If (ilx < 3 .or. ily < 3 .or. ilz < 3) Call error(305)

     mxcell = Max(mxcell,(ilx+5)*(ily+5)*(ilz+5))
  End If

End Subroutine set_bounds
