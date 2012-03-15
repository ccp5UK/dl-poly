Subroutine scan_control                              &
           (mxrdf,mxvdw,rvdw,mxmet,rmet,mxter,rcter, &
           imcon,imc_n,cell,xhi,yhi,zhi,             &
           l_str,l_vv,l_n_e,l_n_r,l_n_v,l_ind,       &
           dvar,rcut,rbin,mxstak,                    &
           nstfce,mxspl,alpha,kmaxa1,kmaxb1,kmaxc1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for raw scanning the contents of the control file
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode,mxnode,gcheck
  Use setup_module, Only : nread,nrite,pi,zero_plus
  Use parse_module, Only : get_line,get_word,lower_case,word_2_real
  Use msd_module

  Implicit None

  Logical,           Intent( InOut ) :: l_n_e
  Logical,           Intent(   Out ) :: l_str,l_vv,l_n_r,l_n_v,l_ind
  Integer,           Intent( In    ) :: mxrdf,mxvdw,mxmet,mxter,imcon
  Integer,           Intent(   Out ) :: imc_n,mxstak, &
                                        nstfce,mxspl,kmaxa1,kmaxb1,kmaxc1
  Real( Kind = wp ), Intent( In    ) :: xhi,yhi,zhi,rcter
  Real( Kind = wp ), Intent( InOut ) :: rvdw,rmet,cell(1:9)
  Real( Kind = wp ), Intent(   Out ) :: dvar,rcut,rbin,alpha

  Logical                :: carry,safe,lrcut,lrvdw,lrmet, &
                            lelec,lrdf,lvdw,lmet,l_n_m,lter
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Integer                :: itmp
  Real( Kind = wp )      :: celprp(1:10),cut,eps,fac,tol,tol1,rbin1

! default spline for SPME and minimum real space cutoff

  Integer,           Parameter :: mxspl_def = 8
  Real( Kind = wp ), Parameter :: rcut_def  = 1.5_wp

! default reading indices options

  l_ind=.true.

! density variation parameter default

  dvar = 1.0_wp

! strict flag

  l_str = .true.

! slab option default

  imc_n = imcon

! integration flavour - velocity verlet assumed

  l_vv = .true.

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

  rbin1 = 0.05_wp
  rbin  = rbin1

! Frequency of the SPME k-space evaluation

  nstfce = -1 ! None defined

! Ewald sum parameters defaults

  mxspl = 0
  alpha = 0.0_wp
  kmaxa1 = 0
  kmaxb1 = 0
  kmaxc1 = 0

! default stack size

  mxstak = 1

! Set safe flag

  safe=.true.

! Open the simulation input file

  If (idnode == 0) Inquire(File='CONTROL', Exist=safe)
  If (mxnode > 1) Call gcheck(safe)
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

! read density variation option

     Else If (word(1:7) == 'densvar') Then

        Call get_word(record,word)
        dvar = word_2_real(word)
        dvar = 1.0_wp + Abs(dvar)/100.0_wp

! read real space cut off

     Else If (word(1:3) == 'cut' .or. word(1:4) == 'rcut') Then

        lrcut = .true.
        Call get_word(record,word)
        rcut = word_2_real(word)
        lrcut = (rcut > zero_plus) ! if zero or nothing is entered

! read vdw cutoff

     Else If (word(1:4) == 'rvdw') Then

        lrvdw=.true.
        Call get_word(record,word)
        If (word(1:3) == 'cut') Call get_word(record,word)
        If (rvdw > 1.0e-6_wp) Then
           rvdw = Min(rvdw,word_2_real(word))
        Else
           rvdw = word_2_real(word)
        End If
        lrvdw = (rvdw > zero_plus) ! if zero or nothing is entered

! read binsize option

     Else If (word(1:7) == 'binsize') Then

        Call get_word(record,word)
        rbin = Abs(word_2_real(word))

! read stack size

     Else If (word(1:5) == 'stack') Then

        Call get_word(record,word)
        If (word(1:4) == 'size') Call get_word(record,word)

        mxstak = Nint(word_2_real(word))

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

! read finish

     Else If (word(1:6) == 'finish') Then

        carry=.false.

     End If

  End Do

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

  rcut=Max(rcut,rvdw,rmet,2.0_wp*rcter+1.0e-6_wp)

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

! fix cell vectors for image conditions with discontinuties

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
              eps = word_2_real(word)
              eps = Max(Min(Abs(eps),0.5_wp),1.0e-20_wp)

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
              alpha = word_2_real(word)

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

! (1) to Max(rcut,Max(cellwidth*mxspl/kmax)) satisfying SPME b-splines
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
                 rcut=2.0_wp*rcter+1.0e-6_wp
              Else
                 rcut=Max(rcut,2.0_wp*rcter+1.0e-6_wp)
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
                (lrvdw .or. lrmet .or. lter) ) Then
              lrcut=.true.
              rcut=Max(rvdw,rmet,2.0_wp*rcter+1.0e-6_wp)
           End If

! Reset rvdw and rmet when only tersoff potentials are opted for and
! possibly reset rcut to 2.0_wp*rcter+1.0e-6_wp (leaving room for failure)

           If (lter .and. l_n_e .and. l_n_v .and. l_n_m .and. l_n_r) Then
              rvdw=0.0_wp
              rmet=0.0_wp
              If (.not.l_str) Then
                 lrcut=.true.
                 rcut=2.0_wp*rcter+1.0e-6_wp
              End If
           End If

! rcut must exist

           If (.not.lrcut) Call error(382)

! define cut

           cut=rcut+1.0e-6_wp

! fix cell vectors for image conditions with discontinuties

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

        If (rbin < 1.0e-05_wp .or. rbin > rcut/4.0_wp) rbin = Min(rbin1,rcut/4.0_wp)

        carry=.false.

     End If

  End Do

  If (idnode == 0) Close(Unit=nread)

  Return

! CONTROL file does not exist

10 Continue
  Call error(126)
  Return

! No finish in CONTROL file or unexpected break

20 Continue
  Call error(17)

End Subroutine scan_control
