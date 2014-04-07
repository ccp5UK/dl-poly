Subroutine scan_control                                    &
           (rcbnd,mxrdf,mxvdw,rvdw,mxmet,rmet,mxter,rcter, &
           imcon,imc_n,cell,xhi,yhi,zhi,                   &
           mxgana,mxgbnd1,mxgang1,mxgdih1,mxginv1,         &
           l_str,lsim,l_vv,l_n_e,l_n_r,lzdn,l_n_v,l_ind,   &
           rcut,rpad,rbin,mxstak,                          &
           nstfce,mxspl,alpha,kmaxa1,kmaxb1,kmaxc1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of the control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2014
! contrib   - i.j.bush february 2014
! contrib   - a.v.brukhno march 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gcheck
  Use setup_module,       Only : nread,nrite,pi,zero_plus
  Use parse_module,       Only : get_line,get_word,lower_case,word_2_real
  Use msd_module
  Use development_module, Only : l_trm
  Use kim_module,         Only : l_kim,rkim


  Implicit None

  Logical,           Intent( InOut ) :: l_n_e
  Logical,           Intent(   Out ) :: l_str,lsim,l_vv,l_n_r,lzdn,l_n_v,l_ind
  Integer,           Intent( In    ) :: mxrdf,mxvdw,mxmet,mxter,imcon
  Integer,           Intent( InOut ) :: imc_n
  Integer,           Intent(   Out ) :: mxgana,mxgbnd1,mxgang1,mxgdih1,mxginv1, &
                                        mxstak,nstfce,mxspl,kmaxa1,kmaxb1,kmaxc1
  Real( Kind = wp ), Intent( In    ) :: xhi,yhi,zhi,rcter
  Real( Kind = wp ), Intent( InOut ) :: rvdw,rmet,rcbnd,cell(1:9)
  Real( Kind = wp ), Intent(   Out ) :: rcut,rpad,rbin,alpha

  Logical                :: carry,safe,la_ana,la_bnd,la_ang,la_dih,la_inv, &
                            lrcut,lrpad,lrvdw,lrmet,lelec,lrdf,lvdw,lmet,l_n_m,lter,l_exp
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word,akey
  Integer                :: itmp,nstrun
  Real( Kind = wp )      :: celprp(1:10),cut,eps,fac,tol,tol1

  Integer,           Parameter :: mxspl_def = 8           ! default spline for SPME
  Real( Kind = wp ), Parameter :: rcut_def  = 1.0_wp  , & ! minimum real space cutoff
                                  rbin_def  = 0.05_wp , & ! default bin size
                                  rcbnd_def = 2.0_wp      ! minimum bond length for bond analysis

! default reading indices options

  l_ind=.true.

! strict flag

  l_str = .true.

! replay history option (real dynamics = no replay)

  lsim = .true. ! don't replay history

! slab option default

  imc_n = imcon

! integration flavour - velocity verlet assumed

  l_vv = .true.

! default switches for intramolecular analysis grids

  la_ana = .false. ; mxgana  = 0
  la_bnd = .false. ; mxgbnd1 = 0
  la_ang = .false. ; mxgang1 = 0
  la_dih = .false. ; mxgdih1 = 0
  la_inv = .false. ; mxginv1 = 0

! electrostatics and no elctrostatics, rdf and no rdf, vdw and no vdw,
! metal and no metal, tersoff and no tersoff interactions,
! real and binsize defaults

  lelec = .false.
! l_n_e is now first determined in scan_field l_n_e = (.false.)

  lrdf  = (mxrdf > 0)
  l_n_r = .not.lrdf

  lvdw  = (mxvdw > 0)
  l_n_v = .false.
  lrvdw = .false. ! Even though it rvdw may have been read from TABLE

  lmet  = (mxmet > 0)
  l_n_m = .not.lmet
  lrmet = (rmet > 1.0e-6_wp)

  lter  = (mxter > 0)

  lrcut = .false.
  rcut  = 0.0_wp

  lrpad = .false.
  rpad  = 0.0_wp

  rbin  = rbin_def

! Frequency of the SPME k-space evaluation

  nstfce = -1 ! None defined

! Ewald sum parameters defaults

  mxspl = 0
  alpha = 0.0_wp
  kmaxa1 = 0
  kmaxb1 = 0
  kmaxc1 = 0

! default number of steps and expansion option

  nstrun = 0
  l_exp = .false.

! default stack size

  mxstak = 1

! Set safe flag

  safe=.true.

! Open the simulation input file

  If (idnode == 0) Inquire(File='CONTROL', Exist=safe)
  If (mxnode > 1) Call gcheck(safe,"enforce")
  If (.not.safe) Then
     Go To 10
  Else
     If (idnode == 0) Open(Unit=nread, File='CONTROL', Status='old')
  End If

! First Pass.  Get cutoff distances, stacksize and density variation.

  Call get_line(safe,nread,record)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

! read slab option (limiting DD slicing in z direction to 2)

     If      (word(1:4) == 'slab') Then

        If (imcon /= 0 .and. imcon /= 6) imc_n=6

! read real space cut off

     Else If (word(1:3) == 'cut' .or. word(1:4) == 'rcut') Then

        lrcut = .true.
        Call get_word(record,word)
        rcut = Abs(word_2_real(word))
        lrcut = (rcut > zero_plus) ! if zero or nothing is entered

! read real space cut off

     Else If (word(1:3) == 'pad' .or. word(1:4) == 'rpad') Then

        lrpad = .true.
        Call get_word(record,word)
        rpad = Abs(word_2_real(word))
        lrpad = (rpad > zero_plus) ! if zero or nothing is entered

! read vdw cutoff

     Else If (word(1:4) == 'rvdw') Then

        lrvdw=.true.
        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        If (rvdw > 1.0e-6_wp) Then
           rvdw = Min(rvdw,word_2_real(word))
        Else
           rvdw = Abs(word_2_real(word))
        End If
        lrvdw = (rvdw > zero_plus) ! if zero or nothing is entered

! read binsize option

     Else If (word(1:7) == 'binsize') Then

        Call get_word(record,word)
        rbin = Abs(word_2_real(word))

! read replay history option

     Else If (word(1:6) == 'replay') Then

        lsim = .false.

! read number of timesteps

     Else If (word(1:5) == 'steps') Then

        Call get_word(record,word)
        nstrun = Nint(Abs(word_2_real(word)))

! read expansion option

     Else If (word(1:5) == 'nfold') Then

        l_exp = .true.

! read stack size

     Else If (word(1:5) == 'stack') Then

        Call get_word(record,word)
        If (word(1:4) == 'size') Call get_word(record,word)

        mxstak = Nint(Abs(word_2_real(word)))

! read MSD option

     Else If (word(1:6) == 'msdtmp') Then

        l_msd = .true.

! read DL_POLY_2 multiple timestep option (compatibility)
! as DL_POLY_4 infrequent k-space SPME evaluation option

     Else If (word(1:4) == 'mult') Then

        Call get_word(record,word)
        If (word(1:5) == 'times' .or. word(1:4) == 'step') Call get_word(record,word)
        nstfce=Max(nstfce,Nint(Abs(word_2_real(word))))

! read electrostatics

     Else If (word(1:5) == 'ewald' .or. word(1:4) == 'spme') Then

        Call get_word(record,word)

        If (word(1:5) == 'evalu') Then

! infrequent k-space SPME evaluation

           Call get_word(record,word)
           If (word(1:5) == 'every') Call get_word(record,word)
           nstfce=Max(nstfce,Nint(Abs(word_2_real(word))))

        Else

           lelec = .true.

        End If

     Else If (word(1:6) == 'distan') Then

        lelec = .true.

     Else If (word(1:4) == 'coul') Then

        lelec = .true.

     Else If (word(1:5) == 'shift') Then

        lelec = .true.

     Else If (word(1:8) == 'reaction') Then

        lelec = .true.

! read "no vdw", "no elec" and "no str" options

     Else If (word(1:2) == 'no') Then

        Call get_word(record,word)

        If (word(1:3) == 'vdw') Then

           l_n_v = .true.

        Else If (word(1:4) == 'elec') Then

           l_n_e = .true.

        Else If (word(1:3) == 'ind' ) Then

           l_ind=.false.
           If (idnode == 0) Write(nrite,"(/,1x,a)") "no index (reading in CONFIG) option on"

        Else If (word(1:3) == 'str' ) Then

           l_str=.false.

        End If

! read integration flavour

     Else If (word(1:8) == 'integrat') Then

        Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:4) == 'type' .or. word(1:6) == 'verlet') Call get_word(record,word)
        If (word(1:8) == 'leapfrog') l_vv=.false.

! read analysis (bonded distributions calculation) option

     Else If (word(1:3) == 'ana') Then

        la_ana = .true.

        Call get_word(record,word)
        akey = word(1:3)

        Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)
        If (word(1:7) == 'collect' .or. word(1:5) == 'sampl' .or. word(1:5) == 'every') Call get_word(record,word)

        Call get_word(record,word)
        If (word(1:5) == 'nbins' .or. word(1:5) == 'ngrid' .or. word(1:4) == 'grid') Call get_word(record,word)

        If      (akey == 'all') Then
           la_bnd = .true.
           la_ang = .true.
           la_dih = .true.
           la_inv = .true.

           mxgana = Abs(Nint(word_2_real(word)))
           mxgbnd1 = mxgana
           mxgang1 = mxgana
           mxgdih1 = mxgana
           mxginv1 = mxgana

           Call get_word(record,word) ! AB: for "rbnd"/"rmax"/"max"/figure
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           rcbnd=Max(rcbnd,word_2_real(word,0.0_wp))
        Else If (akey == 'bon') Then
           la_bnd = .true.

           mxgbnd1 = Max(mxgbnd1,Abs(Nint(word_2_real(word))))
           mxgana  = Max(mxgana,mxgbnd1)

           Call get_word(record,word) ! AB: for "rbnd"/"rmax"/"max"/figure
           If (word(1:4) == 'rbnd' .or. word(1:4) == 'rmax' .or. word(1:3) == 'max') Call get_word(record,word)
           rcbnd=Max(rcbnd,word_2_real(word,0.0_wp))
        Else If (akey == 'ang') Then
           la_ang = .true.

           mxgang1 = Max(mxgang1,Abs(Nint(word_2_real(word))))
           mxgana  = Max(mxgana,mxgang1)
        Else If (akey == 'dih') Then
           la_dih = .true.

           mxgdih1 = Max(mxgdih1,Abs(Nint(word_2_real(word))))
           mxgana  = Max(mxgana,mxgdih1)
        Else If (akey == 'inv') Then
           la_inv = .true.

           mxginv1 = Max(mxginv1,Abs(Nint(word_2_real(word))))
           mxgana  = Max(mxgana,mxginv1)
        End If

! read rdf calculation option

     Else If (word(1:3) == 'rdf') Then

        lrdf = .true.
        l_n_r = .not.lrdf

! read z-density profile option

     Else If (word(1:4) == 'zden') Then

        lzdn = .true.

! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     End If

  End Do

! in case of bonded interactions analysis

  If (la_ana) Then
     If (mxgana == 0) mxgana = -1 ! switch indicator for set_bounds

! mxgana by construction (above) equals the largest grid
! deal with overloaded grid sizes

     If (la_bnd) Then
        mxgbnd1 = mxgana
        rcbnd=Max(rcbnd,rcbnd_def)
     End If
     If (la_ang) mxgang1 = mxgana
     If (la_dih) mxgdih1 = mxgana
     If (la_inv) mxginv1 = mxgana
  End If

! Sort electrostatics

  If (lelec) Then
     If (l_n_e) lelec = .not.l_n_e
  Else
     l_n_e = .true.
  End If

! Sort vdw

  If (lvdw) Then
     If (.not.lrvdw) Then
        lrvdw = (rvdw > 1.0e-6_wp)
        rvdw = Min(rvdw,Max(rcut,rcut_def))
     End If

     If (l_n_v) lvdw = .not.l_n_v
  Else
     l_n_v = .true.
  End If

! Sort rcut as the maximum of all valid cutoffs

  rcut=Max(rcut,rcbnd,rvdw,rmet,rkim,2.0_wp*rcter+1.0e-6_wp)

  If (idnode == 0) Rewind(nread)

! Second Pass.  Sort out cutoffs, cell parameters and Ewald precision.

  Call get_line(safe,nread,record)
  If (.not.safe) Go To 20

  carry = .true.
  Do While (carry)

     Call get_line(safe,nread,record)
     If (.not.safe) Go To 20

     Call lower_case(record)
     Call get_word(record,word)

     If (lelec .and. (word(1:5) == 'ewald' .or. word(1:4) == 'spme')) Then

! Double the kmax size if specified "ewald sum"

        If (word(1:5) == 'ewald') Then
           itmp=2
        Else
           itmp=1
        End If

        Call get_word(record,word)

        If      (word(1:5) == 'evalu')     Then

        Else

! rcut MUST be >= rcut_def

           If (rcut < rcut_def) rcut=rcut_def

! define cut

           cut=rcut+1.0e-6_wp

! fix cell vectors for image conditions with discontinuities

           If (imcon == 0) Then

              cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
              cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

              cell(2) = 0.0_wp
              cell(3) = 0.0_wp
              cell(4) = 0.0_wp
              cell(6) = 0.0_wp
              cell(7) = 0.0_wp
              cell(8) = 0.0_wp

           Else If (imcon == 6) Then

              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

           End If

           If (word(1:9) == 'precision') Then

              Call dcell(cell,celprp)

              Call get_word(record,word)
              eps = Abs(word_2_real(word))
              eps = Max(Min(eps,0.5_wp),1.0e-20_wp)

              Call get_word(record,word)
              mxspl = Abs(Nint(word_2_real(word)))

              tol = Sqrt(Abs(Log(eps*rcut)))
              alpha = Sqrt(Abs(Log(eps*rcut*tol)))/rcut
              tol1 = Sqrt(-Log(eps*rcut*(2.0_wp*tol*alpha)**2))

              fac = 1.0_wp
              If (imcon == 4 .or. imcon == 5 .or. imcon == 7) fac = 2.0_wp**(1.0_wp/3.0_wp)

              kmaxa1 = 2*Nint(0.25_wp + fac*celprp(1)*alpha*tol1/pi)
              kmaxb1 = 2*Nint(0.25_wp + fac*celprp(2)*alpha*tol1/pi)
              kmaxc1 = 2*Nint(0.25_wp + fac*celprp(3)*alpha*tol1/pi)

! rcut is needed directly for the SPME and it MUST exist

              If (.not.lrcut) Call error(433)

           Else

              If (word(1:3) == 'sum') Call get_word(record,word)
              alpha = Abs(word_2_real(word))

              Call get_word(record,word)
              kmaxa1 = itmp*Nint(Abs(word_2_real(word)))

              Call get_word(record,word)
              kmaxb1 = itmp*Nint(Abs(word_2_real(word)))

              Call get_word(record,word)
              kmaxc1 = itmp*Nint(Abs(word_2_real(word)))

              Call get_word(record,word)
              mxspl = Nint(Abs(word_2_real(word)))

! rcut is not needed directly for the SPME but it's needed
! for the link-cell division of the domains
! let's not fail here if no cutoff is specified

           End If

! Get default spline order if none is specified

           If (mxspl == 0) mxspl = mxspl_def

        End If

     Else If (word(1:6) == 'finish') Then

! Sort rvdw

        If ((.not.lrvdw) .and. lvdw) Then
           If (lrcut) Then
              lrvdw=.true.
              rvdw=rcut
           Else
              Call error(402)
           End If
        End If

! Sort rmet

        If ((.not.lrmet) .and. lmet) Then
           If (lrcut .or. lrvdw) Then
              lrmet=.true.
              rmet=Max(rcut,rvdw)
           Else
              Call error(382)
           End If
        End If

! Sort rcut by a reset sequence
! rcut may be >= rcut_def but lrcut may still be .false.
! mxspl = 0 is an indicator for no SPME electrostatics in CONTROL

        If (mxspl /= 0) Then

! (1) to Max(rcut,Max(cell_width*mxspl/kmax)) satisfying SPME b-splines
! propagation width

           If (.not.lrcut) Then
              lrcut=.true.
              Call dcell(cell,celprp)
              rcut=Max( rcut                                         , &
                        Max(celprp(7)*Real(mxspl,wp)/Real(kmaxa1,wp) , &
                            celprp(8)*Real(mxspl,wp)/Real(kmaxb1,wp) , &
                            celprp(9)*Real(mxspl,wp)/Real(kmaxc1,wp)) )
           End If

! Reset rvdw, rmet and rcut when only tersoff potentials are opted for

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r) Then
              rvdw=0.0_wp
              rmet=0.0_wp
              If (.not.l_str) Then
                 rcut=Max(rcbnd,2.0_wp*rcter)+1.0e-6_wp
              Else
                 rcut=Max(rcut,rcbnd,2.0_wp*rcter+1.0e-6_wp)
              End If
           End If

        Else

! no SPME electrostatics is specified but rcut is still needed for
! domain decompositioning and link-celling
! It is needed for the rest of the types of electrostatics

           If ((.not.lrcut) .and. lelec) Call error(382)

! So there is rcut and some kind of electrostatics(-: or neither

! Reset rcut to something sensible if sensible is an option

           If ( ((.not.lrcut) .or. (.not.l_str)) .and. &
                (lrvdw .or. lrmet .or. lter .or. l_kim) ) Then
              lrcut=.true.
              rcut=Max(rcbnd,rvdw,rmet,rkim,2.0_wp*rcter+1.0e-6_wp)
           End If

! Reset rvdw and rmet when only tersoff potentials are opted for and
! possibly reset rcut to 2.0_wp*rcter+1.0e-6_wp (leaving room for failure)

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r .and. (.not.l_kim)) Then
              rvdw=0.0_wp
              rmet=0.0_wp
              If (.not.l_str) Then
                 lrcut=.true.
                 rcut=Max(rcbnd,2.0_wp*rcter)+1.0e-6_wp
              End If
           End If

! rcut must exist

           If (.not.lrcut) Call error(382)

! define cut

           cut=rcut+1.0e-6_wp

! fix cell vectors for image conditions with discontinuities

           If (imcon == 0) Then

              cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
              cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

              cell(2) = 0.0_wp
              cell(3) = 0.0_wp
              cell(4) = 0.0_wp
              cell(6) = 0.0_wp
              cell(7) = 0.0_wp
              cell(8) = 0.0_wp

           Else If (imcon == 6) Then

              cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

           End If

        End If

! Sort rmet=rcut if metal interactions are in play, even if
! they are defined by EAM since rmet can be /= rcut in such
! instances, this can break the NLAST check in metal_ld_set_halo

        If (lmet) rmet = rcut

! Sort rvdw=rcut if VDW interactions are in play

        If (lvdw .and. rvdw > rcut) rvdw = rcut

! Sort rbin as now rcut is already pinned down

        If (rbin < 1.0e-05_wp .or. rbin > rcut/4.0_wp) rbin = Min(rbin_def,rcut/4.0_wp)

        carry=.false.

     End If

  End Do

  If (idnode == 0) Close(Unit=nread)

! When not having dynamics or prepared to terminate
! expanding and not running the small system prepare to exit gracefully

  l_trm = (l_exp .and. nstrun == 0)
  If (((.not.lsim) .or. l_trm) .and. lrpad) rpad=0.0_wp

  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control
