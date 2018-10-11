Module bounds
  Use kinds,           Only : wp,wi
  Use comms,           Only : comms_type
  Use constants,       Only : rt2, rt3, zero_plus, delr_max, delth_max, pi 
  Use domains,         Only : domains_type,map_domains
  Use configuration,   Only : configuration_type, scan_config,read_config
  Use neighbours,      Only : neighbours_type
  Use msd,             Only : msd_type
  Use rdfs,            Only : rdf_type
  Use msd,             Only : msd_type
  Use z_density,       Only : z_density_type
  Use kim,             Only : kim_type
  Use bonds,           Only : bonds_type
  Use angles,          Only : angles_type
  Use dihedrals,       Only : dihedrals_type
  Use inversions,      Only : inversions_type
  Use tersoff,         Only : tersoff_type
  Use development,     Only : development_type
  Use greenkubo,       Only : greenkubo_type
  Use mpole,           Only : mpole_type,POLARISATION_CHARMM
  Use ttm,             Only : ttm_type
  Use numerics,        Only : dcell
  Use control,         Only : scan_control, scan_control_pre
  Use ffield,          Only : scan_field
  Use errors_warnings, Only : error,warning,info
  Use parallel_fft,    Only : adjust_kmax
  Use thermostat,      Only : thermostat_type
  Use statistics,      Only : stats_type
  Use metal,           Only : metal_type
  Use poisson,         Only : poisson_type
  Use tethers,         Only : tethers_type
  Use constraints, Only : constraints_type
  Use pmf, Only : pmf_type
  Use core_shell, Only : core_shell_type
  Use three_body,      Only : threebody_type
  Use vdw,             Only : vdw_type
  Use four_body, Only : four_body_type
  Use rdfs,            Only : rdf_type
  Use external_field, Only : external_field_type
  Use rigid_bodies, Only : rigid_bodies_type
  Use electrostatic, Only : electrostatic_type
  Use ewald, Only : ewald_type
  Use io, Only : io_type
  Use filename, Only : file_type
  Use site, Only : site_type
  Use flow_control, Only : flow_type
  Implicit None

  Private

  Public :: set_bounds

Contains

  Subroutine set_bounds(site,ttm,io,cshell,cons,pmf,stats,thermo,green,devel, &
    msd_data,met,pois,bond,angle,dihedral,inversion,tether,threebody,zdensity, &
    neigh,vdws,tersoffs,fourbody,rdf,mpoles,ext_field,rigid,electro,domain, &
    config,ewld,kim_data,files,flow,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to determine various limits as array bounds,
! grid sizes, paddings, iterations, etc. as specified in setup
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2016
! contrib   - i.j.bush february 2014
! contrib   - m.a.seaton june 2014 (VAF)
! contrib   - m.a.seaton march 2017 (TTM)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( site_type), Intent( InOut ) :: site
  Type( ttm_type ), Intent( InOut ) :: ttm
  Type( io_type ), Intent( InOut ) :: io
  Type( pmf_type ), Intent( InOut ) :: pmf
  Type( core_shell_type ), Intent( InOut ) :: cshell
  Type( constraints_type ), Intent( InOut ) :: cons
  Type( stats_type ), Intent( InOut ) :: stats
  Type( thermostat_type ), Intent( InOut ) :: thermo
  Type( development_type ), Intent( InOut ) :: devel
  Type( greenkubo_type ), Intent( InOut ) :: green
  Type( msd_type ), Intent( InOut ) :: msd_data
  Type( metal_type ), Intent( InOut ) :: met
  Type( poisson_type ), Intent( InOut ) :: pois
  Type( bonds_type ), Intent( InOut ) :: bond
  Type( angles_type ), Intent( InOut ) :: angle
  Type( dihedrals_type ), Intent( InOut ) :: dihedral
  Type( inversions_type ), Intent( InOut ) :: inversion
  Type( tethers_type ), Intent( InOut ) :: tether
  Type( threebody_type ), Intent( InOut ) :: threebody
  Type( z_density_type ), Intent( InOut ) :: zdensity
  Type( neighbours_type ), Intent( InOut ) :: neigh
  Type( vdw_type ), Intent( InOut ) :: vdws
  Type( tersoff_type ), Intent( InOut )  :: tersoffs
  Type( four_body_type ), Intent( InOut ) :: fourbody
  Type( rdf_type ), Intent( InOut ) :: rdf
  Type( mpole_type ), Intent( InOut ) :: mpoles
  Type( external_field_type ), Intent( InOut ) :: ext_field
  Type( rigid_bodies_type ), Intent( InOut ) :: rigid
  Type( electrostatic_type ), Intent( InOut ) :: electro
  Type( configuration_type ), Intent( InOut ) :: config
  Type( domains_type ), Intent( InOut ) :: domain
  Type( ewald_type ), Intent( InOut ) :: ewld
  Type( kim_type ), Intent( InOut ) :: kim_data
  Type( file_type ), Intent( InOut ) :: files(:)
  Type( flow_type ), Intent( InOut ) :: flow
  Type( comms_type ), Intent( InOut ) :: comm

  Logical           :: l_usr,l_n_r,lzdn,lext
  Integer           :: megatm,ilx,ily,ilz,qlx,qly,qlz, &
    mtshl,mtcons,mtrgd,mtteth,mtbond,mtangl,mtdihd,mtinv
  Real( Kind = wp ) :: ats,celprp(1:10),cut,    &
    dens0,dens,fdens,fdvar,  &
    test,vcell,tol,          &
    rcter,rctbp,rcfbp,       &
    xhi,yhi,zhi
  Integer( Kind = wi ) :: mxgrid
  Character( Len = 256 ) :: message

! scan the FIELD file data

  Call scan_field(megatm,site,neigh%max_exclude,mtshl, &
    mtcons,l_usr,mtrgd,mtteth,mtbond,mtangl,mtdihd,mtinv,rcter,rctbp,rcfbp, &
    lext,cshell,cons,pmf,met,bond,angle,dihedral,inversion,tether,threebody, &
    vdws,tersoffs,fourbody,rdf,mpoles,rigid,kim_data,files,electro,comm)

! Get imc_r & set config%dvar

  Call scan_control_pre(config%imc_n,config%dvar,files,comm)

! scan CONFIG file data

  Call scan_config(config,megatm,config%imc_n,config%dvar,config%levcfg,xhi,yhi,zhi,io,domain,files,comm)

! halt execution for unsupported image conditions in DD
! checks for some inherited from DL_POLY_2 are though kept

  If (config%imcon == 4 .or. config%imcon == 5 .or. config%imcon == 7) Call error(514)

! scan CONTROL file data

  Call scan_control(rcter,rigid%max_rigid,config%imcon,config%imc_n,config%cell, &
    xhi,yhi,zhi,config%mxgana,l_n_r,lzdn,config%l_ind,electro%nstfce, &
    ttm,cshell,stats,thermo,green,devel,msd_data,met,pois,bond,angle,dihedral, &
    inversion,zdensity,neigh,vdws,tersoffs,rdf,mpoles,electro,ewld,kim_data, &
    files,flow,comm)

! check integrity of cell vectors: for cubic, TO and RD cases
! i.e. cell(1)=cell(5)=cell(9) (or cell(9)/Sqrt(2) for RD)

  If (config%imcon == 1 .or. config%imcon == 4 .or. config%imcon == 5) Then

     ats = (Abs(config%cell(1))+Abs(config%cell(5)))/2.0_wp
     test = 1.0e-10_wp*ats
     If (Abs(config%cell(1)-ats) > test) Then
       Call error(410)
     End If
     If (Abs(config%cell(5)-ats) > test) Then
       Call error(410)
     End If
     If (config%imcon == 5) Then
        If (Abs(config%cell(9)-ats*rt2) > test) Then
          Call error(410)
        End If
     Else
        If (Abs(config%cell(9)-ats) > test) Then
          Call error(410)
        End If
     End If
  End If

! check integrity of hexagonal prism cell vectors

  If (config%imcon == 7) Then
     If (Abs(config%cell(1)-rt3*config%cell(5)) > 1.0e-6_wp) Then
       Call error(410)
     End If
  End If

! check for diagonal cell matrix if appropriate: imcon=1,2,4,5,7

  If (config%imcon /= 0 .and. config%imcon /= 3 .and. config%imcon /= 6) Then
    If (Any(Abs(config%cell(2:4)) > zero_plus)) Then
      Call error(410)
    End If
    !If (Abs(config%cell(2)) > zero_plus) Call error(410)
    !If (Abs(config%cell(3)) > zero_plus) Call error(410)
    !If (Abs(config%cell(4)) > zero_plus) Call error(410)
    If (Any(Abs(config%cell(6:8)) > zero_plus)) Then
      Call error(410)
    End If
    !If (Abs(config%cell(6)) > zero_plus) Call error(410)
    !If (Abs(config%cell(7)) > zero_plus) Call error(410)
    !If (Abs(config%cell(8)) > zero_plus) Call error(410)
  End If

! calculate dimensional properties of simulation cell
! (for use in link-cells) and ttm%volume and define min cell config%width

  Call dcell(config%cell,celprp)
  config%width=Min(celprp(7),celprp(8),celprp(9))

  config%volm = celprp(10)

  If (config%imcon == 4 .or. config%imcon == 5 .or. config%imcon == 7) config%volm=0.5_wp*config%volm

! check value of cutoff and reset if necessary

  If (config%imcon > 0) Then
     If (config%imcon == 4) config%width=rt3*config%cell(1)/2.0_wp
     If (config%imcon == 5) config%width=config%cell(1)
     If (config%imcon == 6) config%width=Min(celprp(7),celprp(8))

! halt program if potential cutoff exceeds the minimum half-cell config%width

     If (neigh%cutoff > config%width/2.0_wp) Then
        Call warning(3,neigh%cutoff,config%width/2.0_wp,0.0_wp)
        Call error(95)
     End If
  End If


! config%dvar function

  fdvar = config%dvar**1.7_wp

! config%dvar push of dihedral%max_legend and neigh%max_exclude ranges as the usual suspects

  dihedral%max_legend = Nint(fdvar * Real(dihedral%max_legend,wp))
  neigh%max_exclude = Nint(fdvar * Real(neigh%max_exclude,wp))

!!! INTRA-LIKE POTENTIAL PARAMETERS !!!

! maximum number of core-shell units per node

  If (cshell%mxshl > 0 .and. comm%mxnode > 1) Then
     cshell%mxshl = Max(cshell%mxshl,comm%mxnode*mtshl)
     cshell%mxshl = (3*(Nint(fdvar*Real(cshell%mxshl,wp))+comm%mxnode-1))/comm%mxnode
  End If


! maximum number of constraints per node

  If (cons%mxcons > 0 .and. comm%mxnode > 1) Then
     cons%mxcons = Max(cons%mxcons,comm%mxnode*mtcons)
     cons%mxcons = (3*(Nint(fdvar*Real(cons%mxcons,wp))+comm%mxnode-1))/comm%mxnode
  End If


! maximum number of PMF constraints per MD cell - pmf%mxpmf
! (only one type of PMF in only one type of molecule !!!)
! maximum number of atoms per PMF unit - pmf%mxtpmf(1:2) (only two units per pmf)


! maximum number of RBs per node

  If (rigid%max_rigid > 0 .and. comm%mxnode > 1) Then
     rigid%max_rigid = Max(rigid%max_rigid,comm%mxnode*mtrgd)
     rigid%max_rigid = (3*(Nint(fdvar*Real(rigid%max_rigid,wp))+comm%mxnode-1))/comm%mxnode
  End If


! dimension of shared atoms arrays for core-shell, constraint and RB units
! Max=Max#(members-per-unit)*Max#(units-per-domain)/2
! and maximum number of neighbouring domains/nodes in 3D DD (3^3 - 1)

  If (comm%mxnode > 1) Then
     config%mxlshp = Max((2*cshell%mxshl)/2,(2*cons%mxcons)/2,(rigid%max_list*rigid%max_rigid)/2)
     domain%neighbours = 26
  Else ! nothing is to be shared on one node
     config%mxlshp = 0
     domain%neighbours = 0
  End If


! maximum number of tethered atoms per node and tether potential parameters

  If (tether%mxteth > 0) Then
     If (comm%mxnode > 1) Then
        tether%mxteth = Max(tether%mxteth,comm%mxnode*mtteth)
        tether%mxteth = (3*(Nint(fdvar*Real(tether%mxteth,wp))+comm%mxnode-1))/comm%mxnode
     End If
     tether%mxpteth = 3
  Else
     tether%mxpteth = 0
  End If


! maximum number of chemical bonds per node and bond potential parameters

  If (bond%max_bonds > 0) Then
     If (comm%mxnode > 1) Then
        bond%max_bonds = Max(bond%max_bonds,comm%mxnode*mtbond)
        bond%max_bonds = (3*(Nint(fdvar*Real(bond%max_bonds,wp))+comm%mxnode-1))/comm%mxnode
     End If
     bond%max_param = 4
  Else
     bond%max_param = 0
  End If


! maximum number of bond angles per node and angular potential parameters

  If (angle%max_angles > 0) Then
     If (comm%mxnode > 1) Then
        angle%max_angles = Max(angle%max_angles,comm%mxnode*mtangl)
        angle%max_angles = (3*(Nint(fdvar*Real(angle%max_angles,wp))+comm%mxnode-1))/comm%mxnode
     End If
     angle%max_param = 6
  Else
     angle%max_param = 0
  End If


! maximum number of torsion angles per node and dihedral potential parameters

  If (dihedral%max_angles > 0) Then
     If (comm%mxnode > 1) Then
        dihedral%max_angles = Max(dihedral%max_angles,comm%mxnode*mtdihd)
        dihedral%max_angles = (3*(Nint(fdvar*Real(dihedral%max_angles,wp))+comm%mxnode-1))/comm%mxnode
        dihedral%max_angles = dihedral%max_angles + (dihedral%max_angles+4)/5 ! allow for 25% higher density
     End If
     dihedral%max_param = 7
  Else
     dihedral%max_param = 0
  End If


! maximum number of inversions per node and inversion potential parameters
! allow for 20% higher density

  If (inversion%max_angles > 0) Then
     If (comm%mxnode > 1) Then
        inversion%max_angles = Max(inversion%max_angles,comm%mxnode*mtinv)
        inversion%max_angles = (3*(Nint(fdvar*Real(inversion%max_angles,wp))+comm%mxnode-1))/comm%mxnode
        inversion%max_angles = inversion%max_angles + (inversion%max_angles+4)/5 ! allow for 25% higher density
     End If
     inversion%max_param = 3
  Else
     inversion%max_param = 0
  End If



!!! GRIDDING PARAMETERS !!!

! Set grids for opted intramolecular distribution analysis if unset
! SO THEY ARE SWITCHES FOR EXISTENCE TOO

  If (config%mxgana > 0) Then
     If (bond%bin_pdf == -1) Then
        If (bond%bin_tab > 0) Then
           bond%bin_pdf = bond%bin_tab-4
        Else
           bond%bin_pdf = Nint(bond%rcut/delr_max)
        End If
     End If
     If (angle%bin_adf == -1) Then
        If (angle%bin_tab > 0) Then
           angle%bin_adf = angle%bin_tab-4
        Else
           angle%bin_adf = Nint(180.0_wp/delth_max)
        End If
     End If
     If (dihedral%bin_adf == -1) Then
        If (dihedral%bin_tab > 0) Then
           dihedral%bin_adf = dihedral%bin_tab-4
        Else
           dihedral%bin_adf = Nint(360.0_wp/delth_max)
        End If
     End If
     If (inversion%bin_adf == -1) Then
        If (inversion%bin_tab > 0) Then
           inversion%bin_adf = inversion%bin_tab-4
        Else
           dihedral%bin_adf = Nint(180.0_wp/delth_max)
        End If
     End If
     config%mxgana = Max(bond%bin_pdf,angle%bin_adf,dihedral%bin_adf,inversion%bin_adf)
   End If
   config%mxtana = 0 ! initialise for buffer size purposes, set in read_field

! maximum number of rdf potentials (rdf%max_rdf = rdf%max_rdf)
! rdf%max_grid - maximum dimension of rdf%rdf and z-density arrays

  If ((.not. l_n_r) .or. lzdn) Then
     If (((.not. l_n_r) .and. rdf%max_rdf == 0) .and. (vdws%max_vdw > 0 .or. met%max_metal > 0)) &
        rdf%max_rdf = Max(vdws%max_vdw,met%max_metal) ! (vdws,met) == rdf scanning
     rdf%max_grid = Nint(neigh%cutoff/rdf%rbin)
  Else
     rdf%max_grid = 0 ! RDF and Z-density function MUST NOT get called!!!
  End If

! RDFs particulars for USR (umbrella sampling restraints)

  If (l_usr) Then
     rdf%cutoff_usr   = 0.45_wp*config%width
     rdf%max_grid_usr = Nint(rdf%cutoff_usr/rdf%rbin)      ! allows for up to ~75% system ttm%volume shrinkage
     rdf%cutoff_usr   = Real(rdf%max_grid_usr,wp)*rdf%rbin ! round up and beautify for Andrey Brukhno's sake
  Else
     rdf%cutoff_usr   = 0.0_wp
     rdf%max_grid_usr = 0 ! decider on calling USR RDF
  End If

! maximum of all maximum numbers of grid points for all grids - used for mxbuff

  mxgrid = Max(config%mxgana,vdws%max_grid,met%maxgrid,rdf%max_grid,rdf%max_grid_usr,1004,Nint(neigh%cutoff/delr_max)+4)

! grids setting and overrides

! maximum number of grid points for bonds

  bond%bin_tab = Merge(bond%bin_tab,Min(bond%bin_tab,Max(1004,Nint(bond%rcut/delr_max)+4)),bond%bin_tab < 0)

! maximum number of grid points for angles

  angle%bin_tab = Merge(angle%bin_tab,Min(angle%bin_tab,Nint(180.0_wp/delth_max)+4),angle%bin_tab < 0)

! maximum number of grid points for dihedrals

  dihedral%bin_tab = Merge(dihedral%bin_tab,Min(dihedral%bin_tab,Nint(360.0_wp/delth_max)+4),dihedral%bin_tab < 0)

! maximum number of grid points for inversions

  inversion%bin_tab = Merge(inversion%bin_tab,Min(inversion%bin_tab,Nint(180.0_wp/delth_max)+4),inversion%bin_tab < 0)

! maximum number of grid points for electrostatics

  electro%ewald_exclusion_grid = Merge(-1,Max(1004,Nint(neigh%cutoff/delr_max)+4),electro%no_elec)

! maximum number of grid points for vdw interactions - overwritten

  vdws%max_grid = Merge(-1,Max(1004,Nint(vdws%cutoff/delr_max)+4),vdws%no_vdw)

! maximum number of grid points for metal interactions

  met%maxgrid = Max(met%maxgrid,1004,Nint(met%rcut/delr_max)+4)

! maximum number of grid points for tersoff interaction arrays

  tersoffs%max_grid = Merge(-1,Max(1004,Nint(rcter/delr_max)+4),tersoffs%max_ter <= 0)

! maximum of all maximum numbers of grid points for all grids - used for mxbuff

  mxgrid = Max(mxgrid,bond%bin_tab,angle%bin_tab,dihedral%bin_tab, &
    inversion%bin_tab,electro%ewald_exclusion_grid,vdws%max_grid,met%maxgrid,tersoffs%max_grid)



!!! INTER-LIKE POTENTIAL PARAMETERS !!!

! maximum number of vdw potentials and parameters

  If (vdws%max_vdw > 0) Then
     vdws%max_vdw = vdws%max_vdw+1
     vdws%max_param = 7
  Else
     vdws%max_param = 0
  End If


! maximum number of metal potentials and parameters

  If (met%max_metal > 0) Then
     met%max_metal = met%max_metal+1
     met%max_param = 9
  Else
     met%max_param = 0
  End If


! maximum number of tersoff potentials (tersoffs%max_ter = tersoffs%max_ter) and parameters

  If (tersoffs%max_ter > 0) Then
     If      (tersoffs%key_pot == 1) Then
        tersoffs%max_param = 11
     Else If (tersoffs%key_pot == 2) Then
        tersoffs%max_param = 16
     End If
  Else
     tersoffs%max_param = 0
  End If


! maximum number of three-body potentials and parameters

  If (threebody%mxtbp > 0) Then
    threebody%mx2tbp = (site%mxatyp*(site%mxatyp+1))/2
    threebody%mxtbp  = threebody%mx2tbp*site%mxatyp
     If (rctbp < 1.0e-6_wp) rctbp=0.5_wp*neigh%cutoff

     threebody%mxptbp = 5
  Else
     threebody%mx2tbp = 0
     threebody%mxtbp  = 0

     threebody%mxptbp = 0
  End If


! maximum number of four-body potentials and parameters

  If (fourbody%max_four_body > 0) Then
    fourbody%mx3fbp = (site%mxatyp*(site%mxatyp+1)*(site%mxatyp+2))/6
    fourbody%max_four_body  = fourbody%mx3fbp*site%mxatyp
     If (rcfbp < 1.0e-6_wp) rcfbp=0.5_wp*neigh%cutoff

     fourbody%max_param = 3
  Else
     fourbody%mx3fbp = 0
     fourbody%max_four_body  = 0

     fourbody%max_param = 0
  End If


! maximum number of external field parameters

  If (lext) Then
     ext_field%max_param = 6
  Else
     ext_field%max_param = 0
  End If



! DD PARAMETERS - by hypercube mapping of MD cell onto machine resources
! Dependences: MD cell config%widths (explicit) and machine resources (implicit)

  Call map_domains(config%imc_n,celprp(7),celprp(8),celprp(9),domain,comm)

  Call info(' ',.true.)
  Write(message,'(a,3(i6,1x))') 'node/domain decomposition (x,y,z): ', &
    domain%nx,domain%ny,domain%nz
  Call info(message,.true.)

  If (neigh%padding > zero_plus) Then

! define cut

     cut=neigh%cutoff+1.0e-6_wp

! Provide advise on decomposition

     qlx=Int(celprp(7)/cut)
     qly=Int(celprp(8)/cut)
     qlz=Int(celprp(9)/cut)

     Write(message,'(a,i6,a,3(i0,a))')                           &
       'pure cutoff driven limit on largest possible decomposition:', qlx*qly*qlz , &
       ' nodes/domains (', qlx,',',qly,',',qlz,')'
     Call info(message,.true.)

     qlx=Max(1,qlx/2)
     qly=Max(1,qly/2)
     qlz=Max(1,qlz/2)

     Write(message,'(a,i6,a,3(i0,a))')                           &
       'pure cutoff driven limit on largest balanced decomposition:', qlx*qly*qlz , &
       ' nodes/domains (', qlx,',',qly,',',qlz,')'
     Call info(message,.true.)

  End If

10 Continue ! possible neigh%cutoff redefinition...

! Define link-cell cutoff (minimum config%width)

  neigh%cutoff_extended = neigh%cutoff + neigh%padding

! define cut

  cut=neigh%cutoff_extended+1.0e-6_wp

! Provide advise on decomposition

  qlx=Int(celprp(7)/cut)
  qly=Int(celprp(8)/cut)
  qlz=Int(celprp(9)/cut)

  Write(message,'(a,i6,a,3(i0,a))')                       &
    'cutoffs driven limit on largest possible decomposition:', qlx*qly*qlz , &
    ' nodes/domains (', qlx,',',qly,',',qlz,')'
  Call info(message,.true.)

  qlx=Max(1,qlx/2)
  qly=Max(1,qly/2)
  qlz=Max(1,qlz/2)

  Write(message,'(a,i6,a,3(i0,a))')                       &
    'cutoffs driven limit on largest balanced decomposition:', qlx*qly*qlz , &
    ' nodes/domains (', qlx,',',qly,',',qlz,')'
  Call info(message,.true.)

! calculate link cell dimensions per node

  ilx=Int(domain%nx_recip*celprp(7)/cut)
  ily=Int(domain%ny_recip*celprp(8)/cut)
  ilz=Int(domain%nz_recip*celprp(9)/cut)

! print link cell algorithm and check for violations or...

  Write(message,'(a,3i6)') "link-config%cell decomposition 1 (x,y,z): ",ilx,ily,ilz
  Call info(message,.true.)

  tol=Min(0.05_wp,0.005_wp*neigh%cutoff)                                        ! tolerance
  test = 0.02_wp * Merge( 1.0_wp, 2.0_wp, ewld%bspline > 0)                    ! 2% (w/ SPME/PS) or 4% (w/o SPME/PS)
  cut=Min(domain%nx_recip*celprp(7),domain%ny_recip*celprp(8),domain%nz_recip*celprp(9))-1.0e-6_wp ! domain size

  If (ilx*ily*ilz == 0) Then
     If (devel%l_trm) Then ! we are prepared to exit gracefully(-:
        neigh%cutoff = cut   ! - neigh%padding (was zeroed in scan_control)
        Write(message,'(a)') &
          "real space cutoff reset has occurred, early run termination is due"
        Call warning(message,.true.)
        Go To 10
     Else
        If (cut < neigh%cutoff) Then
           Write(message,'(a)') 'neigh%cutoff <= Min(domain config%width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
           Call warning(message,.true.)
           Call error(307)
        Else ! neigh%padding is defined & in 'no strict' mode
           If (neigh%padding > zero_plus .and. (.not.flow%strict)) Then ! Re-set neigh%padding with some slack
              neigh%padding = Min( 0.95_wp * (cut - neigh%cutoff) , test * neigh%cutoff)
              neigh%padding = Real( Int( 100.0_wp * neigh%padding ) , wp ) / 100.0_wp
              If (neigh%padding < tol) neigh%padding = 0.0_wp ! Don't bother
              Go To 10
           Else
              Write(message,'(a)') 'neigh%cutoff <= Min(domain config%width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
              Call warning(message,.true.)
              Call error(307)
           End If
        End If
     End If
  Else ! push/reset the limits in 'no strict' mode
     If (.not.flow%strict) Then
        If (.not.(met%max_metal == 0 .and. electro%no_elec .and. vdws%no_vdw .and. &
          rdf%max_rdf == 0 .and. kim_data%active)) Then
           ! 2b link-cells are needed
           If (comm%mxnode == 1 .and. Min(ilx,ily,ilz) < 2) Then
              ! catch & handle exception
              neigh%padding = 0.95_wp * (0.5_wp*config%width - neigh%cutoff - 1.0e-6_wp)
              ! round up
              neigh%padding = Real( Int( 100.0_wp * neigh%padding ) , wp ) / 100.0_wp
           End If

           If (neigh%padding <= zero_plus) Then ! When neigh%padding is undefined give it some value
              If (Int(Real(Min(ilx,ily,ilz),wp)/(1.0_wp+test)) >= 2) Then ! good non-exception
                 neigh%padding = test * neigh%cutoff
                 neigh%padding = Real( Int( 100.0_wp * neigh%padding ) , wp ) / 100.0_wp
                 If (neigh%padding > tol) Go To 10
              Else ! not so good non-exception
                 neigh%padding = Min( 0.95_wp * ( Min ( domain%nx_recip * celprp(7) / Real(ilx,wp) , &
                                               domain%ny_recip * celprp(8) / Real(ily,wp) , &
                                               domain%nz_recip * celprp(9) / Real(ilz,wp) ) &
                                         - neigh%cutoff - 1.0e-6_wp ) , test * neigh%cutoff )
              neigh%padding = Real( Int( 100.0_wp * neigh%padding ) , wp ) / 100.0_wp ! round up
              End If
           End If

           If (neigh%padding > zero_plus) Then
              If (neigh%padding < tol) neigh%padding = 0.0_wp ! Don't bother
           End If
        Else
           If (neigh%padding >= zero_plus) neigh%padding = 0.0_wp ! Don't bother
        End If

        neigh%cutoff_extended = neigh%cutoff + neigh%padding ! recalculate neigh%cutoff_extended respectively
     End If
  End If

  ! Ensure padding is large enough for KIM model
  If (kim_data%padding_neighbours_required) Then
    If (neigh%padding < kim_data%influence_distance) Then
      neigh%padding = kim_data%influence_distance
        neigh%cutoff_extended = neigh%cutoff + neigh%padding
    End If
  End If

  neigh%unconditional_update = (neigh%padding > zero_plus) ! Determine/Detect conditional VNL updating at start

  If (ilx < 3 .or. ily < 3 .or. ilz < 3) Call warning(100,0.0_wp,0.0_wp,0.0_wp)

! get total link cells per domain (boundary padding included)
! total link-cells per node/domain is ncells = (ilx+2)*(ily+2)*(ilz+2)
! allow for more (possible b-spline SPME triggered increase in nlast),
! total link-cells per node/domain is ncells = (ilx+3)*(ily+3)*(ilz+3)
! allow for thermal expansion of unsettled systems
! total link-cells per node/domain is ncells = (ilx+4)*(ily+4)*(ilz+4)

  If      (config%imcon == 0                ) Then
     neigh%max_cell = Nint((fdvar**4) * Real((ilx+4)*(ily+4)*(ilz+4),wp))
  Else If (config%imcon == 6 .or. config%imc_n == 6) Then
     neigh%max_cell = Nint((fdvar**3) * Real((ilx+4)*(ily+4)*(ilz+4),wp))
  Else
     neigh%max_cell = Nint((fdvar**2) * Real((ilx+4)*(ily+4)*(ilz+4),wp))
  End If



! SPME electrostatics particularities

! qlx,qly,qlz - SPME fictional link-cell dimensions postulating that:
! domain%nx <= ewld%fft_dim_a/ewld%bspline, domain%ny <= ewld%fft_dim_b/ewld%bspline, domain%nz <= ewld%fft_dim_c/ewld%bspline.
! Otherwise, this node's b-splines in SPME will need 'positive halo'
! that is not on the immediate neighbouring nodes in negative
! directions but beyond them (which may mean self-halo in some cases)

  ewld%fft_dim_a = ewld%fft_dim_a1
  ewld%fft_dim_b = ewld%fft_dim_b1
  ewld%fft_dim_c = ewld%fft_dim_c1

  qlx = ilx
  qly = ily
  qlz = ilz

! ewld%bspline = 0 is an indicator for no SPME or Poisson Solver electrostatics in CONTROL

  If (ewld%bspline /= 0) Then

! ensure (ewld%fft_dim_a,ewld%fft_dim_b,ewld%fft_dim_c) consistency between the DD
! processor grid (map_domains is already called) and the grid
! method or comment out adjustments if using ewald_spme_force~

     Call adjust_kmax( ewld%fft_dim_a, domain%nx )
     Call adjust_kmax( ewld%fft_dim_b, domain%ny )
     Call adjust_kmax( ewld%fft_dim_c, domain%nz )

! Calculate and check ql.

     qlx = Min(qlx , ewld%fft_dim_a/(ewld%bspline*domain%nx))
     qly = Min(qly , ewld%fft_dim_b/(ewld%bspline*domain%ny))
     qlz = Min(qlz , ewld%fft_dim_c/(ewld%bspline*domain%nz))

     If (.not.neigh%unconditional_update) Then
        ewld%bspline1=ewld%bspline
     Else
        ewld%bspline1=ewld%bspline+Ceiling((neigh%padding*Real(ewld%bspline,wp))/neigh%cutoff)

! Redifine ql.

        qlx = Min(ilx , ewld%fft_dim_a/(ewld%bspline1*domain%nx))
        qly = Min(ily , ewld%fft_dim_b/(ewld%bspline1*domain%ny))
        qlz = Min(ilz , ewld%fft_dim_c/(ewld%bspline1*domain%nz))
     End If

! Hard luck, giving up

    If (qlx*qly*qlz == 0) Then
      Write(message,'(a,i6,a,3(i0,a))') &
        'SPME driven limit on largest possible decomposition:',  &
        (ewld%fft_dim_a/ewld%bspline1)*(ewld%fft_dim_b/ewld%bspline1)*(ewld%fft_dim_c/ewld%bspline1) ,           &
        ' nodes/domains (', ewld%fft_dim_a/ewld%bspline1,',',ewld%fft_dim_b/ewld%bspline1,',',ewld%fft_dim_c/ewld%bspline1,')'
      Call info(message)
      Call error(308)
    End If

  End If



! decide on MXATMS while reading CONFIG and scan particle density

  Call read_config(config,megatm,config%levcfg,config%l_ind,flow%strict,neigh%cutoff,config%dvar,xhi,yhi, &
    zhi,dens0,dens,io,domain,files,comm)

! Create f(fdvar,dens0,dens)

  If (comm%mxnode == 1 .or. (config%imcon == 0 .or. config%imcon == 6 .or. config%imc_n == 6)) Then
     fdens = fdvar * (0.65_wp*dens0 + 0.35_wp*dens)
  Else If (Min(ilx,ily,ilz) == 1) Then
     fdens = fdvar * (0.50_wp*dens0 + 0.50_wp*dens)
  Else
     fdens = fdvar * (0.35_wp*dens0 + 0.65_wp*dens)
  End If

! density variation affects the link-cell arrays' dimension
! more than domains(+halo) arrays' dimensions, in case of
! events of extreme collapse in atomic systems (aggregation)

! neigh%max_list is the maximum length of link-cell neigh%list (dens * 4/3 pi neigh%cutoff_extended^3)
! + 75% extra tolerance - i.e f(dens0,dens)*(7.5/3)*pi*neigh%cutoff_extended^3

  neigh%max_list = Nint(fdens*2.5_wp*pi*neigh%cutoff_extended**3)
  neigh%max_list = Min(neigh%max_list,megatm-1) ! neigh%max_exclude

  If (neigh%max_list < neigh%max_exclude-1) Then
     Call warning(6,Real(neigh%max_list,wp),Real(neigh%max_exclude,wp),0.0_wp)
     neigh%max_list=neigh%max_exclude-1
  End If

! get link-cell volume

  vcell = config%volm / (Real(ilx*ily*ilz,wp) * Real(comm%mxnode,wp))

! get averaged link-cell particle number, boosted by fdens
! + 25% extra tolerance

  test = fdens*vcell * 1.25_wp

! set dimension of working coordinate arrays

  config%mxatms = Max(1 , Nint(test * Real((ilx+3)*(ily+3)*(ilz+3),wp)))
  If (comm%mxnode == 1 .or. (config%imcon == 0 .or. config%imcon == 6 .or. config%imc_n == 6)) Then
    config%mxatms = Nint(Min(Real(config%mxatms,wp),Real(27.00_wp,wp)*Real(megatm,wp)))
!  Else If (Min(ilx,ily,ilz) == 1) Then
!    config%mxatms = Nint(Min(Real(config%mxatms,wp),Real(20.25_wp,wp)*Real(megatm,wp)))
  Else
    config%mxatms = Nint(Min(Real(config%mxatms,wp),Real(13.50_wp,wp)*Real(megatm,wp)))
  End If

! maximum number of particles per domain (no halo)

  config%mxatdm = Max(1 , Nint(test * Real((ilx+1)*(ily+1)*(ilz+1),wp)))
  config%mxatdm = Min(config%mxatdm,megatm)

! maximum number of timesteps in stack arrays

  stats%mxstak = Max(100,stats%mxstak)

! maximum number of variables in stack arrays
! update number if the MSD option is used, 51=1+27+...+9+9+1+2+2
! consult statistic_collect for more information

  stats%mxnstk = 51 + site%mxatyp + Merge(2*config%mxatdm,0,msd_data%l_msd)

! maximum dimension of principal transfer buffer

! deport_atomic_data & export_atomic_data (+ metal_ld_export.f90
! defects_reference_export & statistics_connect_deport if used)
! are supposed to be the most MPI/buffer consuming routines

! deporting total per atom

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) / Real(ilx*ily*ilz,wp)
  dens0 = dens0/Max(neigh%cutoff_extended/0.2_wp,1.0_wp)
  domain%mxbfdp = Merge( 2, 0, comm%mxnode > 1) * Nint( Real( &
           config%mxatdm*(18+12 + Merge(3,0,neigh%unconditional_update) + (neigh%max_exclude+1) + &
           Merge(neigh%max_exclude+1 + Merge(neigh%max_exclude+1,0,mpoles%key == POLARISATION_CHARMM),0,mpoles%max_mpoles > 0) + &
           Merge(2*(6+stats%mxstak), 0, msd_data%l_msd)) + 3*green%samp  + &
           4*cshell%mxshl+4*cons%mxcons+(Sum(pmf%mxtpmf(1:2)+3))*pmf%mxpmf+(rigid%max_list+13)*rigid%max_rigid + &
           3*tether%mxteth+4*bond%max_bonds+5*angle%max_angles+8*dihedral%max_angles+6*inversion%max_angles,wp) * dens0)

! statistics connect deporting total per atom

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) / Real(ilx*ily*ilz,wp)
  dens0 = dens0/Max(neigh%cutoff_extended/0.2_wp,1.0_wp)
  config%mxbfss = Merge( 4, 0, comm%mxnode > 1) * Nint( Real(config%mxatdm*(8 + Merge(2*(6+stats%mxstak), 0, msd_data%l_msd)),wp)&
    * dens0)

! exporting single per atom (times 13 up to 35)

!  dens  = Real(((ilx+3)*(ily+3)*(ilz+3))/Min(ilx,ily,ilz)+3,wp) / Real(ilx*ily*ilz,wp)
  dens  = Real(((qlx+2)*(qly+2)*(qlz+2))/Min(qlx,qly,qlz)+2,wp) / Real(qlx*qly*qlz,wp)
  domain%mxbfxp = 2 * Nint(Real(config%mxatdm,wp) * dens) ! included induced dipoles

! shared units single per atom of all shared unit

  dens0 = Real(((ilx+2)*(ily+2)*(ilz+2))/Min(ilx,ily,ilz)+2,wp) - 1.0_wp
  dens0 = dens0/Max(neigh%cutoff_extended/2.0_wp,1.0_wp)
  domain%mxbfsh = Merge( 1, 0, comm%mxnode > 1) * &
    Nint(Real(Max(2*cshell%mxshl,2*cons%mxcons,rigid%max_list*rigid%max_rigid),wp) * dens0)

  config%mxbuff = Max(domain%mxbfdp, 35*domain%mxbfxp, 4*domain%mxbfsh, &
    2*(ewld%fft_dim_a/domain%nx)*(ewld%fft_dim_b/domain%ny)*(ewld%fft_dim_c/domain%nz)+10, &
    stats%mxnstk*stats%mxstak, mxgrid, rdf%max_grid,  &
    rigid%max_list*Max(rigid%max_rigid,rigid%max_type),  &
    rigid%max_type*(4+3*rigid%max_list), 10000 )

! reset (increase) link-cell maximum (neigh%max_cell)
! if tersoff or three- or four-body potentials exist

  If (tersoffs%max_ter > 0 .or. threebody%mxtbp > 0 .or. fourbody%max_four_body > 0) Then
     cut=neigh%cutoff+1.0e-6_wp ! reduce cut
     If (tersoffs%max_ter > 0) cut = Min(cut,rcter+1.0e-6_wp)
     If (threebody%mxtbp > 0) cut = Min(cut,rctbp+1.0e-6_wp)
     If (fourbody%max_four_body > 0) cut = Min(cut,rcfbp+1.0e-6_wp)

     ilx=Int(domain%nx_recip*celprp(7)/cut)
     ily=Int(domain%ny_recip*celprp(8)/cut)
     ilz=Int(domain%nz_recip*celprp(9)/cut)

     Write(message,'(a,3i6)') "link-ccell decomposition 2 (x,y,z): ",ilx,ily,ilz
     Call info(message,.true.)

     If (ilx < 3 .or. ily < 3 .or. ilz < 3) Call error(305)

     If      (config%imcon == 0                ) Then
        neigh%max_cell = Max(neigh%max_cell,Nint((fdvar**4) * Real((ilx+5)*(ily+5)*(ilz+5),wp)))
     Else If (config%imcon == 6 .or. config%imc_n == 6) Then
        neigh%max_cell = Max(neigh%max_cell,Nint((fdvar**3) * Real((ilx+5)*(ily+5)*(ilz+5),wp)))
     Else
        neigh%max_cell = Max(neigh%max_cell,Nint((fdvar**2) * Real((ilx+5)*(ily+5)*(ilz+5),wp)))
     End If
  End If

! two-temperature model: determine number of CITs
! in x- and y-directions based on number in z-direction
! and system size

  ttm%delz     = config%cell(9)/Real(ttm%ntsys(3),wp)
  ttm%ntsys(1) = Nint(config%cell(1)/ttm%delz)
  ttm%ntsys(2) = Nint(config%cell(5)/ttm%delz)
  ttm%delx     = config%cell(1)/Real(ttm%ntsys(1),wp)
  ttm%dely     = config%cell(5)/Real(ttm%ntsys(2),wp)
  ttm%volume   = ttm%delx*ttm%dely*ttm%delz
  ttm%rvolume  = 1.0_wp/ttm%volume

! Check number of electronic temperature cells is greater than/
! equal to number of ionic temperature cells

  If (Any(ttm%eltsys<ttm%ntsys)) Call error(670)

! If ttm%redistribute option selected, check for sufficient electronic temperature
! cells to redistribute energy when ionic tmeperature cells are switched off:
! if not available, switch off this option

  If (ttm%redistribute .and. (ttm%eltsys(1)<ttm%ntsys(1)+2 &
    .or. ttm%eltsys(2)<ttm%ntsys(2)+2 .or. ttm%eltsys(3)<ttm%ntsys(3)+2)) Then
    Call warning(500,0.0_wp,0.0_wp,0.0_wp)
    ttm%redistribute = .false.
  End If

! Calculate average atomic density: if not overridden by
! 'ttm atomdens' directive in CONTROL file, will be used
! to convert specific heat capacities to ttm%volumetric
! heat capacity etc.

  ttm%sysrho = Real(megatm,Kind=wp)/(config%cell(1)*config%cell(5)*config%cell(9))

End Subroutine set_bounds

End Module bounds
