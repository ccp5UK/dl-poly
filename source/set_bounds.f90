Subroutine set_bounds                                        &
           (levcfg,imcon,l_vv,l_str,l_n_e,l_n_r,l_n_v,l_ind, &
           rcut,rvdw,rmet,rbin,nstfce,alpha,width)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to determine various limits: array bounds,
! iteration and others as specified in setup_module
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode
  Use setup_module
  Use domains_module, Only : map_domains,nprx,npry,nprz
  Use config_module,  Only : cfgname,imc_n,cell,volm
  Use msd_module

  Implicit None

  Logical,           Intent(   Out ) :: l_vv,l_str,l_n_e,l_n_r,l_n_v,l_ind
  Integer,           Intent(   Out ) :: levcfg,imcon,nstfce
  Real( Kind = wp ), Intent(   Out ) :: rcut,rvdw,rmet,rbin,alpha,width

  Logical           :: lter,ltbp,lfbp
  Integer           :: ilx,ily,ilz,qlx,qly,qlz,megatm
  Real( Kind = wp ) :: ats,celprp(1:10),cut,dens0,dens,dvar, &
                       test,vcell,rcter,rctbp,rcfbp,         &
                       sidex,sidey,sidez,xhi,yhi,zhi

! define zero+ and half- (setup_module)

  zero_plus  = Nearest(0.0_wp,1.0_wp)
  half_minus = Nearest(0.5_wp,-1.0_wp)

! scan the FIELD file data

  Call scan_field                                     &
           (l_n_e,mxsite,mxatyp,megatm,mxtmls,mxexcl, &
           mxtshl,mxshl,mxfshl,mxtcon,mxcons,mxfcon,  &
           mxtpmf,mxpmf,mxfpmf,                       &
           mxtrgd,mxrgd,mxlrgd,mxfrgd,                &
           mxtteth,mxteth,mxftet,                     &
           mxtbnd,mxbond,mxfbnd,mxtang,mxangl,mxfang, &
           mxtdih,mxdihd,mxfdih,mxtinv,mxinv,mxfinv,  &
           mxrdf,mxgrid,mxvdw,rvdw,mxmet,rmet,        &
           mxter,rcter,mxtbp,rctbp,mxfbp,rcfbp)

! scan CONFIG file data

  Call scan_config(megatm,imc_n,cfgname,levcfg,imcon,cell,xhi,yhi,zhi)

! halt execution for unsupported image conditions in DD
! checks for some inherited from DL_POLY_2 are though kept

  If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(514)

! scan CONTROL file data

  Call scan_control                                  &
           (mxrdf,mxvdw,rvdw,mxmet,rmet,mxter,rcter, &
           imcon,imc_n,cell,xhi,yhi,zhi,             &
           l_str,l_vv,l_n_e,l_n_r,l_n_v,l_ind,       &
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



!!! SITE PARAMETERS !!!

! maximum number of different sites in system

  mxsite = Max(1,mxsite)

! set maximum number of unique atom types

  mxatyp = Max(1,mxatyp)

! maximum number of molecule types

  mxtmls = Max(1,mxtmls)

! maximum number of excluded atoms per atom

  mxexcl = Max(1,mxexcl)



!!! INTRA-LIKE POTENTIAL PARAMETERS !!!

! maximum number of core-shell units per node

  If (mxshl > 0 .and. mxnode > 1) &
     mxshl = (3*(Nint((dvar**1.7_wp)*Real(mxshl,wp))+mxnode-1))/mxnode


! maximum number of constraints per node

  If (mxcons > 0 .and. mxnode > 1) &
     mxcons = (3*(Nint((dvar**1.7_wp)*Real(mxcons,wp))+mxnode-1))/mxnode


! maximum number of PMF constraints per MD cell - mxpmf
! (only one type of PMF in only one type of molecule !!!)
! maximum number of atoms per PMF unit - mxtpmf(1:2) (only two units per pmf)


! maximum number of RBs per node

  If (mxrgd > 0 .and. mxnode > 1) &
     mxrgd = (3*(Nint((dvar**1.7_wp)*Real(mxrgd,wp))+mxnode-1))/mxnode



! dimension of shared atoms arrays for core-shell, constraint and RB units
! Max=Max#(members-per-unit)*Max#(units-per-domain)/2
! and maximum number of neighbouring domains/nodes in 3D DD (3^3 - 1)

  mxlshp = Max((2*mxshl)/2,(2*mxcons)/2,(mxlrgd*mxrgd)/2)
  mxproc = 26

! nothing is to be shared on one node

  If (mxnode == 1) Then
     mxlshp = 0
     mxproc = 0
  End If


! maximum number of tethered atoms per node

  If (mxteth > 0 .and. mxnode > 1) &
     mxteth = (3*(Nint((dvar**1.7_wp)*Real(mxteth,wp))+mxnode-1))/mxnode

! maximum number of parameters for tethers

  mxpteth = 3


! maximum number of chemical bonds per node

  If (mxbond > 0 .and. mxnode > 1) &
     mxbond = (3*(Nint((dvar**1.7_wp)*Real(mxbond,wp))+mxnode-1))/mxnode

! maximum number of parameters for bond potentials

  mxpbnd = 4


! maximum number of bond angles per node

  If (mxangl > 0 .and. mxnode > 1) &
     mxangl = (3*(Nint((dvar**1.7_wp)*Real(mxangl,wp))+mxnode-1))/mxnode

! maximum number of parameters for angular potentials

  mxpang = 6


! maximum number of torsion angles per node

  If (mxdihd > 0 .and. mxnode > 1) &
     mxdihd = (3*(Nint((dvar**1.7_wp)*Real(mxdihd,wp))+mxnode-1))/mxnode

! maximum number of parameters for dihedral (torsional) potentials

  mxpdih = 7


! maximum number of inversion potentials per node

  If (mxinv > 0 .and. mxnode > 1) &
     mxinv = (3*(Nint((dvar**1.7_wp)*Real(mxinv,wp))+mxnode-1))/mxnode

! maximum number of parameters for inversion potentials

  mxpinv = 3



!!! INTER-LIKE POTENTIAL PARAMETERS !!!

! maximum number of grid points in potentials arrays

  mxgrid = Max(mxgrid,1000,Int(rcut/0.01_wp+0.5_wp)+4)



! maximum number of rdf potentials (mxrdf = mxrdf)
! mxgrdf - maximum dimension of rdf and z-density arrays

  mxgrdf = Nint(rcut/rbin)



! maximum number of vdw potentials

  If (mxvdw > 0) mxvdw = mxvdw+1

! maximum number of parameters for vdw potentials

  mxpvdw = 5



! maximum number of metal potentials

  If (mxmet > 0) mxmet = mxmet+1

! maximum number of parameters for metal potentials

  mxpmet = 9



! maximum number of tersoff potentials (mxter = mxter)

  If (mxter == 0) Then
     lter = .false.
  Else
     lter = .true.
  End If

! maximum number of parameters for tersoff potentials

  mxpter = 11



! maximum number of three-body potentials

  If (mxtbp == 0) Then
     ltbp   = .false.
     mx2tbp = 0
  Else
     ltbp   = .true.
     mx2tbp = (mxatyp*(mxatyp+1))/2
     mxtbp  = mx2tbp*mxatyp
     If (rctbp < 1.0e-6_wp) rctbp=0.5_wp*rcut
  End If

! maximum number of three-body potential parameters

  mxptbp = 5



! maximum number of four-body potentials

  If (mxfbp == 0) Then
     lfbp   = .false.
     mx3fbp = 0
  Else
     lfbp   = .true.
     mx3fbp = (mxatyp*(mxatyp+1)*(mxatyp+2))/6
     mxfbp  = mx3fbp*mxatyp
     If (rcfbp < 1.0e-6_wp) rcfbp=0.5_wp*rcut
  End If

! maximum number of four-body potential parameters

  mxpfbp = 3



! maximum number of external field parameters

  mxpfld = 5



! DD PARAMETERS - by hypercube mapping of MD cell onto machine resources
! Dependences: MD cell widths (explicit) and machine resources (implicit)

  Call map_domains(imc_n,celprp(7),celprp(8),celprp(9))

  If (idnode == 0) Write(nrite,'(/,/,1x,a,3i6)') 'node/domain decomposition (x,y,z): ', nprx,npry,nprz



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

  If (ilx*ily*ilz == 0) Call error(307)

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
     qlx = Min(ilx , kmaxa1/(mxspl*nprx))
     qly = Min(ily , kmaxb1/(mxspl*npry))
     qlz = Min(ilz , kmaxc1/(mxspl*nprz))

     If (qlx*qly*qlz == 0) Call error(308)

! ensure (kmaxa,kmaxb,kmaxc) consistency with what DD
! (map_domains is already called) and DaFT are capble of
! or comment out adjustments if using ewald_spme_force~

     Call adjust_kmax( kmaxa, nprx )
     Call adjust_kmax( kmaxb, npry )
     Call adjust_kmax( kmaxc, nprz )
  End If



! decide on MXATMS while reading CONFIG and scan particle density

  Call read_config &
           (megatm,levcfg,imcon,l_ind,l_str,rcut,dvar,xhi,yhi,zhi,dens0,dens)

! density variation affects the link-cell arrays' dimension
! more than domains(+halo) arrays' dimensions, in case of
! events of extreme collapse in atomic systems (aggregation)

! mxlist is the maximum length of link-cell list (dens * 4/3 pi rcut^3)
! + 75% extra tolerance (f(dens0,dens) * 7.5/3 pi rcut^3)

  If (mxnode == 1 .or. (imcon == 0 .or. imcon == 6 .or. imc_n == 6)) Then
     mxlist = Nint( (dvar**1.7_wp) * (0.75_wp*dens0+0.25_wp*dens)*2.5_wp*pi*rcut**3)
  Else
     mxlist = Nint( (dvar**1.7_wp) * (0.25_wp*dens0+0.75_wp*dens)*2.5_wp*pi*rcut**3)
  End If
  mxlist = Min(mxlist,megatm-1)

! get link-cell volume

  vcell = volm / (Real(ilx*ily*ilz,wp) * Real(mxnode,wp))

! get averaged link-cell particle number, boosted by densvar
! + 25% extra tolerance

  test = (dvar**1.7_wp) * dens*vcell * 1.25_wp

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

  mxnstk = 45+mxatyp

! update maximum number of variables in stack arrays if MSD option is used

  If (l_msd) mxnstk = mxnstk+2*mxatdm

! maximum dimension of principal transfer buffer

! deport_atomic_data & export_atomic_data (and ewald_spme_force~,
! and defects_reference_export & metal_ld_export.f90 if used)
! are supposed to be the most MPI/buffer consuming routines

! REMOVE 0*(...) in mxbuff estimation (penultimate line) if using ewald_spme_force~

  mxbuff = Max( (Merge( 2, 0, mxnode > 1)                   * &
                  ((17+9+mxexcl+6+mxstak)*mxatdm            + &
                   3*mxshl + 3*mxcons                       + &
                   (mxtpmf(1)+mxtpmf(2)+2)*mxpmf            + &
                   (mxlrgd+19)*Max(mxtrgd,mxrgd)            + &
                   2*mxteth + 3*mxbond + 4*mxangl           + &
                   5*mxdihd + 5*mxinv))                     / &
                  (Min(ilx,ily,ilz)*Max(Nint(rcut),1))      , &
                ((Mod(mxnode-1,1)+1)*11*mxatms)             / &
                 Min(qlx,qly,qlz)                           , &
                2*(kmaxa/nprx)*(kmaxb/npry)*(kmaxc/nprz)+10 , &
                0*(2*kmaxa*kmaxb*kmaxc+10)                  , &
                mxnstk*mxstak , mxgrid , mxgrdf , mxtrgd*mxlrgd , 10000 )



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
