      program parset

c
c***********************************************************************
c
c     dl_poly utility program to prepare the dl_params.inc file
c     for specific dl_poly applications
c
c     copyright daresbury laboratory 1995
c     author - w.smith, t.forester jan 1995
c
c***********************************************************************
c

      implicit real*8(a-h,o-z)
      logical loglnk
      dimension cell(9),celprp(10)
      pi = 3.1415926536d0

      write(*,'(a)')'enter the minimum number of nodes for job'
      read(*,*)minnode
      write(*,'(a)')'enter the maximum number of nodes for job'
      read(*,*)mxproc

c
c     open main printing file

      open(7,file ='new_params.inc')

c
c     write main file header

      write(7,'(a,a)')'c*********************************************'
     x   ,'************************'
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     dl_poly insert file specifying array sizes for the'
      write(7,'(a)')
     x   'c     entire package'
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     copyright - daresbury laboratory 1994'
      write(7,'(a)')
     x   'c     authors - w. smith & t. forester november  1994.'
      write(7,'(a)')'c'
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     note the following internal units apply everywhere'
      write(7,'(a)')'c'
      write(7,'(a,a)')
     x   'c     unit of time      (to)    =          1 x 10**(-12)',
     x   ' seconds'
      write(7,'(a,a)')
     x   'c     unit of length    (lo)    =          1 x 10**(-10)',
     x   ' metres'
      write(7,'(a,a)')
     x   'c     unit of mass      (mo)    = 1.6605402  x 10**(-27)',
     x   ' kilograms'
      write(7,'(a,a)')
     x   'c     unit of charge    (qo)    = 1.60217733 x 10**(-19)',
     x   ' coulombs'
      write(7,'(a,a)')
     x   'c     unit of energy    (eo)    = 1.6605402  x 10**(-23)',
     x   ' joules'
      write(7,'(a,a)')
     x   'c     unit of pressure  (po)    = 1.6605402  x 10**(  7)',
     x   ' pascals'
      write(7,'(a,a)')
     x   'c                               = 163.842151            ',
     x   ' atmospheres'
      write(7,'(a)')'c'
      write(7,'(a,a)')'c*********************************************'
     x   ,'************************'
      write(7,'(a)')'c'

c
c     double precision specification

      write(7,'(a)')'      implicit real*8 (a-h,o-z)'

c
c     numerical and universal constants

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     standard pi values'
      write(7,'(a)')' '
      write(7,'(a,a)')
     x  '      parameter (pi =3.141592653589793d0,',
     x  'sqrpi =1.7724538509055159d0)'

      write(7,'(a)')'c'
      write(7,'(a,a)')'c     conversion factor for coulombic terms in',
     x   ' internal units'
      write(7,'(a,a)')'c     i.e. (unit(charge)**2/(4 pi eps0 unit',
     x   '(length))/unit(energy)'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (r4pie0 =138935.4835d0)'
      write(7,'(a)')'c'
      write(7,'(a)')'c     boltzmann constant in internal units'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (boltz =8.31451115d-1)'
      write(7,'(a)')'c'
      write(7,'(a,a)')'c     conversion factor for pressure from ',
     x   'internal units to kbar'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (prsunt = 1.63842151d-1)'

c
c     minimum number of nodes specification

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     minimum number of nodes for program execution'
      write(7,'(a)')' '
      write(7,'(a,i4,a)')'      parameter (minnode =',minnode,')'

c
c     maximum number of nodes (used in parallel constraint algorithm)

      write(7,'(a)')'c'
      write(7,'(a,a)')
     x   'c     maximum number of nodes (used in parallel',
     x   ' constraint algorithm)'
      write(7,'(a)')' '
      write(7,'(a,i4,a)')'      parameter (mxproc =',mxproc,')'

c
c     specify default i/o channels

      write(7,'(a)')'c'
      write(7,'(a)')'c     main input channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nread = 7)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     main output channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nrite = 8)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     force field input channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nfield = 9)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     configuration file input channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nconf = 10)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     statistical data file output channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nstats = 20)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     trajectory history file channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nhist = 21)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     acummulators restart dump file'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nrest = 22)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     tabulated potential file channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (ntable = 23)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     rdf file channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nrdfdt = 24)'

      write(7,'(a)')'c'
      write(7,'(a)')'c     z density file channel'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (nzdndt = 25)'

      write(7,'(a)')'c'
      write(7,'(a,a)')'c     data dumping interval in event of system',
     x   ' crash'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (ndump = 1000)'

c
c     input the system force field

      call fldscan
     x   (mxtmls,mxtcon,mxtbnd,mxtang,mxtdih,mxungp,mxngp,
     x   mxgrp,mxatms,mxgatm,mxneut,mxteth,mxtet1,mxsvdw,mxvdw,mxn1,
     x   mxtbp,mxexcl,mxsite,mxbond,mxcons,mxangl,mxdihd,mxtshl,
     x   mxshl,mxpmf,mxspmf,mxtinv,mxinv,mxfbp,ngrid,rctbp,rcfbp)

c
c     input  system size

      call cfgscan(imcon,volm,xhi,yhi,zhi,cell)
c
c     scan control file for required parameters

      call conscan(imcon,mxstak,kmaxa,kmaxb,kmaxc,rcut,rvdw,delr,cell)

c
c     set maximum number of atoms

      mxatms=max(1,mxatms)
      write(7,'(a)')'c'
      write(7,'(a)')'c     maximum number of atoms'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxatms =',mxatms,')'

c
c     set dimension of verlet test arrays

      mslst = max(1,(mxatms+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of atoms in verlet test arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mslst =',mslst,')'

c
c     set dimension of working coordinate arrays

      msatms=max(1,(mxatms+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of atoms in working coordinate arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (msatms =',msatms,')'

c
c     maximum number of neutral groups

      mxneut1=max(mxneut+1,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of neutral groups'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxneut =',mxneut1,')'

c
c     maximum number of molecule types

      mxtmls=max(mxtmls,1)
      write(7,'(a)')'c'
      write(7,'(a)')'c     maximum number of molecule types'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtmls =',mxtmls,')'

c
c     maximum number of specified bondlength constraints

      mxtcon=max(mxtcon,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of specified bondlength constraints'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtcon =',mxtcon,')'

c
c     maximum number of chemical bond potentials

      mxtbnd=max(mxtbnd,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of chemical bond potentials'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtbnd =',mxtbnd,')'

c
c     maximum number of different bond angle potentials

      mxtang=max(mxtang,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of different bond angle potentials'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtang =',mxtang,')'

c
c     maximum number of different torsional potentials

      mxtdih=max(mxtdih,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of different torsional potentials'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtdih =',mxtdih,')'

c
c     maximum number of different inversion potentials

      mxtinv=max(mxtinv,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of different inversion potentials'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtinv =',mxtinv,')'

c
c     maximum number of unique rigid body units

      mxungp=max(mxungp,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of unique rigid body units'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxungp =',mxungp,')'

c
c     maximum number of tethered atom potentials

      mxteth=max(mxteth,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of tethered atom potentials'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxteth =',mxteth,')'

c
c     set maximum number of unique atom types

      write(7,'(a)')'c'
      write(7,'(a)')'c     maximum number of unique atom types'
      write(7,'(a)')' '
      write(7,'(a,i4,a)')'      parameter (mxsvdw =',mxsvdw,')'

c
c     maximum number of different pair potentials

      mxvdw=max(mxvdw,1)+1
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of different pair potentials'
      write(7,'(a)')' '
      write(7,'(a,i4,a)')'      parameter (mxvdw  =',mxvdw,')'

c
c     maximum number of three body potentials
      if(mxtbp.eq.0)then

         mxtbp1=1
         mx2tbp=1

      else

         mx2tbp=mxvdw
         mxtbp1=mxvdw*mxsvdw

      endif
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of three body potentials'
      write(7,'(a)')' '
      write(7,'(a,i4,a)')'      parameter (mxtbp =',mxtbp1,')'
      write(7,'(a,i4,a)')'      parameter (mx2tbp =',mx2tbp,')'

c
c     maximum number of four body potentials

      if(mxfbp.eq.0)then

         mxfbp1=1
         mx3fbp=1

      else

         mx3fbp=(mxsvdw*(mxsvdw+1)*(mxsvdw+2))/6
         mxfbp1=mxsvdw*mx3fbp

      endif

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of four body potentials'
      write(7,'(a)')' '
      write(7,'(a,i4,a)')'      parameter (mxfbp =',mxfbp1,')'
      write(7,'(a,i4,a)')'      parameter (mx3fbp =',mx3fbp,')'

c
c     maximum number of angular potential parameters

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of angular potential parameters'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxpang = 4)'

c
c     maximum number of three body potential parameters

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of three body potential parameters'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxptbp = mxpang+1)'

c
c     maximum number of four body potential parameters

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of four body potential parameters'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxpfbp = 3)'

c
c     maximum number of parameters for dihedrals
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of parameters for dihedrals'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxpdih = 5)'

c
c     maximum number of parameters for inversion potentials
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of parameters for inversion potentials'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxpinv = 2)'

c
c     maximum number of parameters for bond potentials
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of parameters for bond potentials'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxpbnd = 4)'

c
c     maximum number of parameters for vdw potentials
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of parameters for vdw potentials'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxpvdw = 5)'

c
c     maximum number of external field parameters

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of external field parameters'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxfld = 10)'

c
c     maximum number of excluded atoms per atom

      mxexcl=max(mxexcl,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of excluded atoms per atom'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxexcl =',mxexcl,')'

c
c     maximum number of different sites in system

      mxsite=max(mxsite,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of different sites in system'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxsite =',mxsite,')'

c
c     maximum number of chemical bonds per node

      mxbond =max(1,(mxbond+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of chemical bonds per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxbond =',mxbond,')'

c
c     maximum number of bond angles per node

      mxangl =max(1,(mxangl+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of bond angles per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxangl =',mxangl,')'

c
c     maximum number of torsion angles per node

      mxdihd =max(1,(mxdihd+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of torsion angles per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxdihd =',mxdihd,')'

c
c     maximum number of inversion potentials per node

      mxinv =max(1,(mxinv+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of inversion potentials per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxinv =',mxinv,')'

c
c     maximum number of constraints per node

      mxcons = max(1,2*((mxcons+minnode-1)/minnode))
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of constraints per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxcons =',mxcons,')'

c
c     maximum number of tethered atoms per node

      msteth = max(1,(mxtet1+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of tethered atoms per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (msteth =',msteth,')'

c
c     maximum size for working arrays for bonds, angles, dihedrals
c     inversion potentials, tethers and core-shell units

      msbad = max(mxbond,mxangl,mxdihd,mxinv,msteth,mxshl)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum size for bond, angle, dihedral, etc.  arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (msbad =',msbad,')'

c
c     maximum number of grid points in potentials arrays

      if(ngrid.eq.0)ngrid = max(500,int(rvdw/0.01d0+0.5d0)+4)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of grid points in potentials arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxgrid = ',ngrid,')'

c
c     maximum dimension of rdf arrays

      nrdf =max(128,int(rvdw/0.05d0+0.5d0))
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum dimension of rdf arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxrdf = ',nrdf,')'

c
c     maximum number of rigid groups in system

      mxgrp=max(mxgrp,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of rigid groups in system'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxgrp =',mxgrp,')'

c
c     maximum number of rigid groups per node

      msgrp=max(1,(mxgrp+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of rigid groups per node'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (msgrp =',msgrp,')'

c
c     maximum number of sites per rigid unit

      mxngp=max(mxngp,3)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of sites per rigid unit'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxngp =',mxngp,')'

c
c     maximum number of sites in rigid units

      mxgatm = max(1,mxgatm)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of sites in rigid units'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxgatm =',mxgatm,')'

c
c     maximum iterations in quaternion integration

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum iterations in quaternion integration'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxquat = 100)'

c
c     maximum number of timesteps in stack arrays

      mxstak=max(100,mxstak)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of timesteps in stack arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxstak =',mxstak,')'

c
c     maximum number of variables in stack arrays

      mxnstk =45+mxsvdw
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of variables in stack arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxnstk =',mxnstk,')'

c
c     maximum dimension of atomic transfer buffer

      mxbuff =max(6*mxatms,8*(mxcons+1),8*(mxgrp+1),mxnstk*mxstak,
     x  mxebuf,ngrid)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum dimension of atomic transfer buffer'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxbuff =',mxbuff,')'

c
c     maximum number of cycles in shake algorithm

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of shake cycles'
      write(7,'(a)')' '
      write(7,'(a)')'      parameter (mxshak = 100)'

c
c     dimension of shake shared atoms array

      mxlshp = max(mxcons*2,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     dimension of shake shared atoms array'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxlshp = ',mxlshp,')'

c
c     set dimension of working coordinate arrays in ewald sum

      mxewld =1
      if(kmaxa.gt.0) mxewld = msatms
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum dimension in ewald working arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxewld =',mxewld,')'

c
c     maximum number of k vectors in ewald sum (kmaxb and kmaxc)

      write(7,'(a)')'c'
      write(7,'(a,a)')
     x   'c     maximum reciprocal space vectors index in b and',
     x   ' c direction'
      write(7,'(a)')' '
      write(7,'(a,i3,a)')'      parameter (kmaxb =',kmaxb,')'
      write(7,'(a,i3,a)')'      parameter (kmaxc =',kmaxc,')'

c
c     maximum size of verlet neighbour/link cell list for each atom

c     decide if link-cells in use or not for van der waals interactions

      loglnk=.true.
      cut=rcut+delr
      dens=dble(mxatms)/volm
      ratio=1.5d0*dens*(4.d0*pi/3.d0)*cut**3
      mxlist=min(nint(ratio),(mxatms+1)/2)
      if(imcon.eq.0) then

        cell(1) = max(2.d0*xhi+cut,3.d0*cut)
        cell(5) = max(2.d0*yhi+cut,3.d0*cut)
        cell(9) = max(2.d0*zhi+cut,3.d0*cut)

      endif
      if(imcon.eq.6)then

        cell(9) = max(2.d0*zhi+cut,3.d0*cut,cell(9))

      endif

      call dcell(cell,celprp)

      ilx = max(3,int(celprp(7)/cut))
      ily = max(3,int(celprp(8)/cut))
      ilz = max(3,int(celprp(9)/cut))
      ncells = ilx*ily*ilz

      if(ncells.eq.27) loglnk=.false.
      if(imcon.eq.4.or.imcon.eq.5) loglnk = .false.
      if(mxneut.gt.0.and.ncells.le.36) loglnk=.false.

      mxcell=1
      if(loglnk)then

        mxcell=ncells
        mxlist=14*nint(1.5d0*dens*celprp(10)/dble(ncells))

      endif

      if(mxneut.gt.0) mxlist = (mxneut1+1)/2

c     set link cells for three and four body forces

      if(mxtbp.gt.0.or.mxfbp.gt.0)then

        if(mxtbp.gt.0)cut=min(cut,rctbp)
        if(mxfbp.gt.0)cut=min(cut,rcfbp)
        ilx = max(3,int(celprp(7)/cut))
        ily = max(3,int(celprp(8)/cut))
        ilz = max(3,int(celprp(9)/cut))
        ncells = ilx*ily*ilz
        mxcell=max(mxcell,ncells)

      endif

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum size of verlet neighbour list for each atom'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxlist =',mxlist,')'

      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of link cells'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxcell =',mxcell,')'

c
c     maximum size for coordinate difference arrays

      mxxdf = max(mxlist,mxatms,mxcons,mxn1*(mxneut1+1)/2)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum size for coordinate difference arrays'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxxdf =',mxxdf,')'

c
c     maximum number of core-shell unit types

      mxtshl=max(mxtshl,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of core-shell unit types'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxtshl =',mxtshl,')'

c
c     maximum number of core-shell units

      mxshl=max(mxshl,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of core-shell units'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxshl =',mxshl,')'

c
c     potential of mean force array parameter

      mxpmf=max(1,mxpmf)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     potential of mean force array parameter'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxpmf =',mxpmf,')'

c
c     number of pmf constraints on a processor

      mspmf=max(1,(mxpmf+minnode-1)/minnode)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     number of pmf constraints on a processor'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mspmf =',mspmf,')'

c
c     maximum number of sites to define pmf units

      mxspmf=max(mxspmf,1)
      write(7,'(a)')'c'
      write(7,'(a)')
     x   'c     maximum number of sites to define pmf units'
      write(7,'(a)')' '
      write(7,'(a,i10,a)')'      parameter (mxspmf =',mxspmf,')'

      close (7)
      stop

      end
      subroutine conscan
     x  (imcon,mxstak,kmaxa,kmaxb,kmaxc,rcut,rvdw,delr,cell)
c
c***********************************************************************
c
c     dl_poly subroutine for scanning the contents of the control file
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith  november   1994
c
c     based on modified simdef routine:
c     author   - t. forester       may 1993
c     keywords - t.forester        april 1994
c
c***********************************************************************
c

      implicit real*8(a-h,o-z)
      parameter (mega=1000000,pi=3.1415926536d0)

      character*100 record
      dimension celprp(10),cell(9)

      mxstak=1
      kmaxa=0
      kmaxb=1
      kmaxc=1
      rcut=0.0
      rvdw=0.0
      delr=0.0
c
c     open the simulation input file

      open (12,file='CONTROL')

      read(12,'(a100)',end=100) record
      write(6,'(a,a)')'CONTROL file header: ',record

      do nrecs= 1,mega

         read(12,'(a100)',end=100) record
         call lowcase(record,100)
         call strip(record,100)

         if(record(1:5).eq.'stack') then

            mxstak= intstr(record,100,idum)

         elseif(record(1:5).eq.'ewald') then

           if(record(7:15).eq.'precision') then

             call dcell(cell,celprp)
             if(rcut.lt.1.d-6)then
               write(*,*)'# using working rcut value of 10 A',
     x           ' to estimate Ewald parameters'
               rcut=10.d0
             endif
             eps = dblstr(record(16:16),64,idum)
             eps = min(abs(eps),0.5d0)
             tol = sqrt(abs(log(eps*rcut)))
             alpha = sqrt(abs(log(eps*rcut*tol)))/rcut
             tol1 = sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
             fac = 1.d0
             if(imcon.eq.4.or.imcon.eq.5) fac = 2.d0**(1.d0/3.d0)
             kmaxa = nint(0.25d0 + fac*celprp(1)*alpha*tol1/pi)
             kmaxb = nint(0.25d0 + fac*celprp(2)*alpha*tol1/pi)
             kmaxc = nint(0.25d0 + fac*celprp(3)*alpha*tol1/pi)

           else

             alpha= dblstr(record,100,idum)
             ilen= 100-idum
             record(1:ilen)= record(idum:100)

             kmaxa=intstr(record,ilen,idum)
             ilen1= ilen-idum
             record(1:ilen1)= record(idum:ilen)

             kmaxb=intstr(record,ilen1,idum)
             ilen= ilen1-idum
             record(1:ilen)= record(idum:ilen1)

            kmaxc=intstr(record,ilen,idum)

          endif

         elseif(record(1:3).eq.'cut') then

            rcut= dblstr(record,100,idum)

         elseif(record(1:4).eq.'rvdw') then

           rvdw= dblstr(record,100,idum)

         elseif(record(1:4).eq.'delr')then

            delr= dblstr(record,100,idum)

         elseif(record(1:6).eq.'finish')then

            go to 100

         endif

      enddo

  100 close (12)

      if(rvdw.eq.0.0) rvdw=rcut

      return
      end
      subroutine cfgscan(imcon,volm,xhi,yhi,zhi,cell)

c***********************************************************************
c
c     dl_poly subroutine for scanning the initial configuration
c     file to determine the number of atoms present
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith  november   1994
c
c***********************************************************************
c

      implicit real*8(a-h,o-z)

      parameter (mega=1000000)

      character*80 header
      character*8 name
      logical lvolm
      dimension cell(9),celprp(10)

      xhi=0.d0
      yhi=0.d0
      zhi=0.d0
      volm=0.d0
      open (10,file='CONFIG')

c
c     read the CONFIG file header

      read(10,'(a80)',end=100)header
      write(6,'(a,a)')'CONFIG file header: ',header
      read(10,'(2i10)',end=100) levcfg,imcon
      lvolm=(imcon.eq.0.or.imcon.eq.6)

c
c     specify molecular dynamics simulation cell

      if(imcon.gt.0)then

        read(10,*,end=100)cell(1),cell(2),cell(3)
        read(10,*,end=100)cell(4),cell(5),cell(6)
        read(10,*,end=100)cell(7),cell(8),cell(9)
        call dcell(cell,celprp)

      endif

      if(.not.lvolm)volm=celprp(10)

      do i=1,mega

        if(levcfg.eq.0)then

          read(10,'(a8)',end=100) name
          read(10,*)xxx,yyy,zzz

        else if(levcfg.eq.1)then

          read(10,'(a8)',end=100) name
          read(10,*)xxx,yyy,zzz
          read(10,*)uuu,vvv,www

        else

          read(10,'(a8)',end=100) name
          read(10,*)xxx,yyy,zzz
          read(10,*)uuu,vvv,www
          read(10,*)uuu,vvv,www

        endif

        if(lvolm)then

          if(i.eq.1)then

            xhi=abs(xxx)
            yhi=abs(yyy)
            zhi=abs(zzz)

          else

            xhi=max(xhi,abs(xxx))
            yhi=max(yhi,abs(yyy))
            zhi=max(zhi,abs(zzz))

          endif

        endif

      enddo

  100 continue

      if(imcon.eq.0)then

        volm=8.d0*xhi*yhi*zhi

      else if(imcon.eq.6)then

        coz=(cell(1)*cell(4)*cell(2)*cell(5)+cell(3)*cell(6))/
     x    (celprp(1)*celprp(2))
        volm=2.d0*zhi*celprp(1)*celprp(2)*sqrt(1.d0-coz**2)

      endif

      close (10)

      return

      end
      subroutine fldscan
     x   (mxtmls,mxtcon,mxtbnd,mxtang,mxtdih,mxungp,mxngp,
     x   mxgrp,mxatms,mxgatm,mxneut,mxteth,mxtet1,mxsvdw,mxvdw,mxn1,
     x   mxtbp,mxexcl,mxsite,mxbond,mxcons,mxangl,mxdihd,mxtshl,
     x   mxshl,mxpmf,mxspmf,mxtinv,mxinv,mxfbp,ngrid,rctbp,rcfbp)

c
c***********************************************************************
c
c     dl_poly routine for scanning the field file to determine the
c     required parameters
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith  november   1994
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)

      parameter (mega=10000)

      character*80 header
      character*40 record
      character*8 name,unique(1000)
      logical check,ltable,lneut

      mxtmls=0
      mxatms=0
      mxgrp=0
      mxtcon=0
      mxtbnd=0
      mxtang=0
      mxtdih=0
      mxtinv=0
      mxpmf=0
      mxspmf=0
      mxungp=0
      mxngp=0
      mxneut=0
      mxn1=0
      nxn1=0
      nold = -1
      mxgatm=0
      mxteth=0
      mxtet1=0
      mxsvdw=0
      mxvdw=0
      mxtbp=0
      mxexcl=0
      mxsite=0
      mxbond=0
      mxcons=0
      mxangl=0
      mxdihd=0
      mxinv=0
      mxshl=0
      mxtshl=0
      mxfbp=0
      ngrid=0
      rctbp=0.d0
      rcfbp=0.d0
      lneut=.false.
      ltable=.false.

c
c     open force field data file

      open (11,file='FIELD')

      read(11,'(a)',end=200)header
      write(6,'(a,a)')'FIELD file header: ',header

c
c     read and process directives from field file

      do nrecs=1,mega

         read(11,'(a40)',end=200) record
         call lowcase(record,40)
         call strip(record,40)
         if(record(1:4).eq.'neut')then

           lneut=.true.

         elseif(record(1:6).eq.'molecu')then

            mxtmls=intstr(record,40,idum)

            do itmols=1,mxtmls

               read(11,*,end=200)

               do itrec=1,mega

                  read(11,'(a40)',end=200) record
                  call lowcase(record,40)
                  call strip(record,40)

                  ksite=0

                  if(record(1:6).eq.'nummol')then

                     nummols=intstr(record,40,idum)

                   elseif(record(1:5).eq.'atoms')then

                     numsit=intstr(record,40,idum)
                     mxatms = mxatms+numsit*nummols
                     mxsite=mxsite+numsit
                     ksite=0
                     do isite=1,numsit

                        if(ksite.lt.numsit)then

                           read(11,'(a8,24x,4i5)',end=200)name,nrept,
     x                       ifrz,nneu

                           if(nrept.eq.0)nrept=1
                           if(lneut)then
                             if(nneu.ne.nold)nxn1=0
                             nxn1 = nxn1+nrept
                             mxn1 = max(mxn1,nxn1)
                             nold = nneu
                           endif

                           if(mxsvdw.eq.0)then

                              mxsvdw=1
                              unique(1)=name

                           else

                              check=.true.
                              do j=1,mxsvdw

                                 if(name.eq.unique(j))check=.false.

                              enddo
                              if(check)then

                                 mxsvdw=mxsvdw+1
                                 unique(mxsvdw)=name

                              endif

                           endif

                           ksite=ksite+nrept

                        endif

                     enddo

                     if(lneut)mxneut = mxneut+nneu*nummols

                  elseif(record(1:5).eq.'shell')then

                     numshl=intstr(record,40,idum)
                     mxtshl=mxtshl+numshl
                     mxshl=mxshl+nummols*numshl

                     do ishls=1,numshl

                        read(11,*,end=100)

                     enddo

                  elseif(record(1:5).eq.'bonds')then

                     numbonds=intstr(record,40,idum)
                     mxtbnd=mxtbnd+numbonds
                     mxbond=mxbond+nummols*numbonds

                     do ibonds=1,numbonds

                        read(11,*,end=100)

                     enddo

                  elseif(record(1:6).eq.'constr')then

                     numcon=intstr(record,40,idum)
                     mxtcon=mxtcon+numcon
                     mxcons=mxcons+nummols*numcon

                     do icon=1,numcon

                        read(11,*,end=100)

                     enddo

                  elseif(record(1:6).eq.'angles')then

                     numang=intstr(record,40,idum)
                     mxtang=mxtang+numang
                     mxangl=mxangl+nummols*numang

                     do iang=1,numang

                        read(11,*,end=100)

                     enddo

                  elseif(record(1:6).eq.'dihedr')then

                     numdih=intstr(record,40,idum)
                     mxtdih=mxtdih+numdih
                     mxdihd=mxdihd+nummols*numdih

                     do idih=1,numdih

                        read(11,*,end=100)

                     enddo

                  elseif(record(1:6).eq.'invers')then

                     numinv=intstr(record,40,idum)
                     mxtinv=mxtinv+numinv
                     mxinv=mxinv+nummols*numinv

                     do iinv=1,numinv

                        read(11,*,end=100)

                     enddo

                  elseif(record(1:5).eq.'rigid')then

                     numgrp=intstr(record,40,idum)
                     mxungp=mxungp+numgrp
                     mxgrp = mxgrp+numgrp*nummols

                     do kgrp=1,numgrp

                        read(11,*,end=200) numgsit
                        mxgatm = mxgatm+numgsit*nummols
                        mxngp=max(mxngp,numgsit)
                        do j=16,numgsit,16

                          read(11,*,end=200)

                        enddo

                     enddo

                  elseif(record(1:4).eq.'teth')then

                     numteth=intstr(record,40,idum)
                     mxteth=mxteth+numteth
                     mxtet1 = mxtet1+numteth*nummols

                     do iteth=1,numteth

                        read(11,*,end=100)

                     enddo

                   elseif(record(1:3).eq.'pmf')then

                     do ipmf = 1,2

                       read(11,'(a40)',end=200) record
                       call strip(record,40)
                       call lowcase(record,40)
                       npmf=intstr(record,40,idum)
                       mxspmf=mxspmf+npmf

                       do jpmf=1,npmf

                         read(11,*,end=200)

                       enddo

                     enddo

                     mxpmf=mxpmf+nummols

                  elseif(record(1:6).eq.'finish')then

                     go to 100

                  endif

               enddo

  100          continue

            enddo

         elseif(record(1:3).eq.'vdw') then

            call strip(record(4:4),37)
            if(record(4:6).eq.'tab')ltable=.true.
            ntpvdw=intstr(record,40,idum)
            mxvdw=max(ntpvdw,(mxsvdw*(mxsvdw+1))/2)
            do itpvdw=1,ntpvdw

               read(11,'(a)',end=200)record
               call lowcase(record,40)
               if(record(17:20).eq.'stch')mxvdw=2*mxvdw
               if(record(17:19).eq.'tab')ltable=.true.
               if(record(17:19).eq.'   ')ltable=.true.

            enddo
            if(ltable)then

              open(13,file='TABLE')
              read(13,*)
              read(13,*)a,b,ngrid

              close (13)
            endif


         elseif(record(1:3).eq.'tbp') then

            ntptbp=intstr(record,40,idum)
            mxtbp=ntptbp

            do itptbp=1,ntptbp

               read(11,'(76x,f12.0)',end=200)rct
               rctbp=max(rctbp,rct)

            enddo

         elseif(record(1:3).eq.'fbp') then

            ntpfbp=intstr(record,40,idum)
            mxfbp=ntpfbp
            do itpfbp=1,ntpfbp

               read(11,'(60x,f12.0)',end=200)rct
               rcfbp=max(rcfbp,rct)

            enddo

         elseif(record(1:6).eq.'extern') then

            read(11,*,end=200)

         elseif(record(1:5).eq.'close')then

            go to 200

         endif

      enddo

  200 close (11)

      if(mxpmf.gt.0)mxpmf=mxatms
      if(mxtcon.gt.0)mxexcl=max(mxexcl,6)
      if(mxtbnd.gt.0)mxexcl=max(mxexcl,6)
      if(mxtang.gt.0)mxexcl=max(mxexcl,16)
      if(mxtdih.gt.0)mxexcl=max(mxexcl,50)
      if(mxtinv.gt.0)mxexcl=max(mxexcl,50)
      if(mxneut.gt.0)mxexcl=max(mxexcl,10*mxn1*mxn1)
      if(mxgrp.gt.0)mxexcl=max(mxexcl,mxngp)

      return
      end
      subroutine dcell(aaa,bbb)

c
c***********************************************************************
c
c     dl_poly subroutine to calculate the dimensional properies of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c
c     copyright daresbury laboratory 1992
c     author - w. smith         july 1992
c
c***********************************************************************
c

      implicit real*8 (a-h,o-z)

      dimension aaa(9),bbb(10)
c
c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))
c
c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))
c
c     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)
c
c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)
c
c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end
      function dblstr(word,len,lst)

c
c***********************************************************************
c
c     dl_poly function for extracting double precisions from a
c     character string.
c     modified from dl_poly function intstr
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     modified  - t. forester april 1994
c
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c              double precision string
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)

      character*1 n(0:9),word(len),ksn,dot,d,e,work(255)
      logical flag,final,ldot,ldd,start

      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/

      sn=1.d0
      ksn='+'
      ten = 10.d0
      one = 1.d0

      call lowcase(word,len)

      dblstr=0.d0
      iexp = 0
      idum =0
      final=.false.
      start=.false.
      ldot = .false.
      ldd = .false.

      do lst=1,len

         flag=.false.

         do j=0,9

            if(n(j).eq.word(lst))then

               dblstr= ten*dblstr+one*dble(j)
               flag=.true.
               start=.true.

            endif

         enddo

         if(flag.and.ksn.eq.'-') sn=-1.d0

          ksn=word(lst)

         if(dot.eq.ksn) then

            flag=.true.
            ten = one
            ldot = .true.

         endif

         if(ldot) then

           one = one/10.d0

        endif

        if(start) then

            if(d.eq.ksn.or.e.eq.ksn) then

               flag=.true.
               idd = lst
               ldd = .true.

            endif

         endif

         if(ldd) goto 100

         if(flag)then

            final=.true.

         else

            if(final) goto 100

         endif

      enddo

  100 dblstr = dblstr*sn

      if(ldd) then

         do i = 1,len-idd
            work(i) = word(i+idd)
         enddo

         len = len-idd
         iexp = intstr(work,len,idum)

      endif

      dblstr = dblstr*(10.d0**iexp)
      lst = lst +idum

      return
      end


      function intstr(word,len,lst)

c
c***********************************************************************
c
c     dl_poly function for extracting integers from a
c     character string
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c
c***********************************************************************
c

      character*1 n(0:9),word(len),ksn
      logical flag,final

      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      ksn='+'
      intstr=0
      final=.false.

      do lst=1,len

         flag=.false.

         do j=0,9

            if(n(j).eq.word(lst))then

               intstr=10*intstr+j
               flag=.true.

            endif

         enddo

         if(flag.and.ksn.eq.'-')isn=-1
         ksn=word(lst)

         if(flag)then

            final=.true.

         else

            if(final)then

               intstr=isn*intstr
               return

            endif

         endif

      enddo

      intstr=isn*intstr

      return

      end
      subroutine strip(string,length)

c***********************************************************************
c
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c
c***********************************************************************

      character*(*) string

      imax = min(length,255)
      do i = 1,imax

         if(string(1:1).eq.' ') then

            do j = 1,imax-1

               string(j:j) = string(j+1:j+1)

            enddo

            string(imax:imax) = ' '

         endif

      enddo

      return
      end
      subroutine lowcase(string,length)

c***********************************************************************
c
c     DL_POLY routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      character*1 string(*)
      character*1 letter

      do i = 1,min(255,length)

         letter = string(i)

         if(letter.eq.'A') letter = 'a'
         if(letter.eq.'B') letter = 'b'
         if(letter.eq.'C') letter = 'c'
         if(letter.eq.'D') letter = 'd'
         if(letter.eq.'E') letter = 'e'
         if(letter.eq.'F') letter = 'f'
         if(letter.eq.'G') letter = 'g'
         if(letter.eq.'H') letter = 'h'
         if(letter.eq.'I') letter = 'i'
         if(letter.eq.'J') letter = 'j'
         if(letter.eq.'K') letter = 'k'
         if(letter.eq.'L') letter = 'l'
         if(letter.eq.'M') letter = 'm'
         if(letter.eq.'N') letter = 'n'
         if(letter.eq.'O') letter = 'o'
         if(letter.eq.'P') letter = 'p'
         if(letter.eq.'Q') letter = 'q'
         if(letter.eq.'R') letter = 'r'
         if(letter.eq.'S') letter = 's'
         if(letter.eq.'T') letter = 't'
         if(letter.eq.'U') letter = 'u'
         if(letter.eq.'V') letter = 'v'
         if(letter.eq.'W') letter = 'w'
         if(letter.eq.'X') letter = 'x'
         if(letter.eq.'Y') letter = 'y'
         if(letter.eq.'Z') letter = 'z'

         string(i) = letter

      enddo

      return
      end
