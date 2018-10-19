Module ewald_spole
  Use kinds,           Only : wp
  Use comms,           Only : comms_type, gcheck, gsum
  Use constants,           Only : r4pie0,sqrpi,twopi,zero_plus
  Use configuration,   Only : configuration_type
  Use particle,        Only : corePart
  Use numerics,        Only : erfcgen, invert, dcell
  Use errors_warnings, Only : error
  Use ewald,           Only : ewald_type,spl_cexp, bspcoe, bspgen, exchange_grid
  Use domains, Only : domains_type
  Use parallel_fft, Only : initialize_fft, pfft, pfft_indices
  Use neighbours, Only : neighbours_type
  Use electrostatic, Only : electrostatic_type
  Implicit None

  Private

  Public :: ewald_real_forces, ewald_spme_forces, ewald_excl_forces, ewald_frzn_forces
Contains

  Subroutine ewald_real_forces(iatm,xxt,yyt,zzt,rrt,engcpe_rl,vircpe_rl,stress, &
      neigh,electro,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating coulombic energy and force terms
    ! in a periodic system using ewald's method
    !
    ! Note: real space terms
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1998
    ! amended   - i.t.todorov april 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_rl,vircpe_rl
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( comms_type),                        Intent( In    ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config


    Integer           :: fail,m,idi,jatm,k
    Real( Kind = wp ) :: chgea,chgprd,rrr,ppp,egamma,   &
      fix,fiy,fiz,fx,fy,fz,          &
      vk0,vk1,vk2,gk0,gk1,gk2,t1,t2, &
      strs1,strs2,strs3,strs5,strs6,strs9

    Character ( Len = 256 )                              :: message


    If (electro%newjob_spole) Then
      electro%newjob_spole = .false.

      fail=0
      Allocate (electro%erc_spole(0:electro%ewald_exclusion_grid),electro%fer_spole(0:electro%ewald_exclusion_grid), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_real_forces allocation failure'
        Call error(0,message)
      End If

      ! interpolation interval

      electro%drewd_spole = neigh%cutoff/Real(electro%ewald_exclusion_grid-4,wp)

      ! reciprocal of interpolation interval

      electro%drewd_spole = 1.0_wp/electro%drewd_spole

      ! generate error function complement tables for ewald sum

      Call erfcgen(neigh%cutoff,electro%alpha,electro%ewald_exclusion_grid,electro%erc_spole,electro%fer_spole)
    End If

    ! initialise potential energy and virial

    engcpe_rl=0.0_wp
    vircpe_rl=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/electro%eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index and charge

        jatm=neigh%list(m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! calculate forces

          k   = Int(rrr*electro%drewd_spole)
          ppp = rrr*electro%drewd_spole - Real(k,wp)

          ! calculate forces using 3pt interpolation

          gk0 = electro%fer_spole(k) ; If (k == 0) gk0 = gk0*rrr
          gk1 = electro%fer_spole(k+1)
          gk2 = electro%fer_spole(k+2)

          t1 = gk0 + (gk1 - gk0)*ppp
          t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

          egamma = (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

          ! calculate forces

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! calculate interaction energy using 3-point interpolation

            vk0 = electro%erc_spole(k)
            vk1 = electro%erc_spole(k+1)
            vk2 = electro%erc_spole(k+2)

            t1 = vk0 + (vk1 - vk0)*ppp
            t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

            engcpe_rl = engcpe_rl + (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

            ! calculate virial

            vircpe_rl = vircpe_rl - egamma*rrr**2

            ! calculate stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

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

    End If

  End Subroutine ewald_real_forces

  Subroutine ewald_spme_forces(engcpe_rc,vircpe_rc,stress,ewld,electro,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating coulombic energy and force terms
    ! in a periodic system using the smooth particle mesh ewald method
    ! by Essmann et al. J. Chem. Phys. 103 (1995) 8577
    !
    ! Note: (fourier) reciprocal space terms
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & w.smith & i.j.bush february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent(   Out ) :: engcpe_rc,vircpe_rc
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Type( ewald_type), Intent( InOut ) :: ewld
    Type( electrostatic_type ), Intent( InOut    ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type), Intent( InOut ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config


    Logical              :: llspl=.true.
    Integer              :: fail(1:4), i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt

    Real( Kind = wp )    :: det,rcell(1:9),celprp(1:10),ralph,rvolm,scale,   &
      rcpcut,rcpct2,strs(1:9),eng,akv,tmp,bb1,bb2,bb3, &
      rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3

    Complex( Kind = wp ) :: vterm

    ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter    :: &
      uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)



    ! B-spline coefficients

    Complex( Kind = wp ), Dimension( : ),   Allocatable       :: ww1,ww2,ww3

    Real( Kind = wp ),    Dimension( : ),   Allocatable       :: csp
    Real( Kind = wp ),    Dimension( : ),   Allocatable       :: txx,tyy,tzz
    Integer,              Dimension( : ),   Allocatable       :: ixx,iyy,izz,it
    Real( Kind = wp ),    Dimension( :,: ), Allocatable       :: bsdx,bsdy,bsdz
    Real( Kind = wp ),    Dimension( :,: ), Allocatable       :: bspx,bspy,bspz



    ! temporary qqc

    Real( Kind = wp )    :: qqc_tmp


    ! DaFT arrays local indices

    Integer              :: j_local, k_local, l_local

    ! message for error
    Character ( Len = 256 ) :: message


    fail=0
    If (electro%newjob_sspme) Then
      electro%newjob_sspme = .false.

      !!! BEGIN DD SPME VARIABLES
      ! 3D charge array construction (bottom and top) indices

      electro%ixb=domain%idx*(ewld%fft_dim_a/domain%nx)+1
      electro%ixt=(domain%idx+1)*(ewld%fft_dim_a/domain%nx)
      electro%iyb=domain%idy*(ewld%fft_dim_b/domain%ny)+1
      electro%iyt=(domain%idy+1)*(ewld%fft_dim_b/domain%ny)
      electro%izb=domain%idz*(ewld%fft_dim_c/domain%nz)+1
      electro%izt=(domain%idz+1)*(ewld%fft_dim_c/domain%nz)

      electro%ixbm1_r=Real(electro%ixb-1,wp)
      electro%ixtm0_r=Nearest( Real(electro%ixt,wp) , -1.0_wp )
      electro%iybm1_r=Real(electro%iyb-1,wp)
      electro%iytm0_r=Nearest( Real(electro%iyt,wp) , -1.0_wp )
      electro%izbm1_r=Real(electro%izb-1,wp)
      electro%iztm0_r=Nearest( Real(electro%izt,wp) , -1.0_wp )

      ! Real values of kmax vectors

      electro%kmaxa_r=Real(ewld%fft_dim_a,wp)
      electro%kmaxb_r=Real(ewld%fft_dim_b,wp)
      electro%kmaxc_r=Real(ewld%fft_dim_c,wp)

      !!! END DD SPME VARIABLES

      !!! BEGIN CARDINAL B-SPLINES SET-UP
      ! allocate the complex exponential arrays

      Allocate (ww1(1:ewld%fft_dim_a),ww2(1:ewld%fft_dim_b),ww3(1:ewld%fft_dim_c), Stat = fail(1))
      If (fail(1) > 0) Then
        Write(message,'(a)') 'ww arrays allocation failure'
        Call error(0,message)
      End If

      ! initialise the complex exponential arrays

      Call spl_cexp(ewld%fft_dim_a,ewld%fft_dim_b,ewld%fft_dim_c,ww1,ww2,ww3)

      ! allocate the global B-spline coefficients and the helper array

      Allocate (electro%bscx(1:ewld%fft_dim_a),electro%bscy(1:ewld%fft_dim_b),electro%bscz(1:ewld%fft_dim_c), Stat = fail(1))
      Allocate (csp(1:ewld%bspline),                              Stat = fail(2))
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'bsc and cse arrays allocation failure'
        Call error(0,message)
      End If

      ! calculate the global B-spline coefficients

      Call bspcoe(ewld,csp,electro%bscx,electro%bscy,electro%bscz,ww1,ww2,ww3)

      ! deallocate the helper array and complex exponential arrays

      Deallocate (csp,         Stat = fail(1))
      Deallocate (ww1,ww2,ww3, Stat = fail(2))
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'cse and ww arrays deallocation failure'
        Call error(0,message)
      End If

      !!! END CARDINAL B-SPLINES SET-UP

      !!! BEGIN DAFT SET-UP
      ! domain local block limits of kmax space

      electro%block_x = ewld%fft_dim_a / domain%nx
      electro%block_y = ewld%fft_dim_b / domain%ny
      electro%block_z = ewld%fft_dim_c / domain%nz

      ! set up the parallel fft and useful related quantities

      Call initialize_fft( 3, (/ ewld%fft_dim_a, ewld%fft_dim_b, ewld%fft_dim_c /), &
        (/ domain%nx, domain%ny, domain%nz /), (/ domain%idx, domain%idy, domain%idz /),   &
        (/ electro%block_x, electro%block_y, electro%block_z /),               &
        comm%comm, electro%context )

      ! set up the indexing arrays for each dimension (NOT deallocated manually)

      Allocate ( electro%index_x( 1:electro%block_x ), Stat = fail(1) )
      Allocate ( electro%index_y( 1:electro%block_y ), Stat = fail(2) )
      Allocate ( electro%index_z( 1:electro%block_z ), Stat = fail(3) )
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'SPME index arrays allocation failure'
        Call error(0,message)
      End If

      Call pfft_indices( ewld%fft_dim_a, electro%block_x, domain%idx, domain%nx, electro%index_x )
      Call pfft_indices( ewld%fft_dim_b, electro%block_y, domain%idy, domain%ny, electro%index_y )
      Call pfft_indices( ewld%fft_dim_c, electro%block_z, domain%idz, domain%nz, electro%index_z )

      ! workspace arrays for DaFT

      Allocate ( electro%qqc_local( 1:electro%block_x, 1:electro%block_y, 1:electro%block_z ), Stat = fail(1) )
      Allocate ( electro%qqq_local( 1:electro%block_x, 1:electro%block_y, 1:electro%block_z ), Stat = fail(2) )
      Allocate ( electro%pfft_work( 1:electro%block_x, 1:electro%block_y, 1:electro%block_z ), Stat = fail(3) )
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'SPME DaFT workspace arrays allocation failure'
        Call error(0,message)
      End If

      !!! END DAFT SET-UP

      ! calculate self-interaction correction

      ewld%engsic=0.0_wp
      Do i=1,config%natms
        ewld%engsic=ewld%engsic+config%parts(i)%chge**2
      End Do
      Call gsum(comm,ewld%engsic)
      ewld%engsic=-r4pie0/electro%eps * electro%alpha*ewld%engsic/sqrpi
    End If

    Allocate (txx(1:config%mxatms),tyy(1:config%mxatms),tzz(1:config%mxatms),                            Stat = fail(1))
    Allocate (ixx(1:config%mxatms),iyy(1:config%mxatms),izz(1:config%mxatms),it(1:config%mxatms),               Stat = fail(2))
    Allocate (bsdx(1:ewld%bspline,1:config%mxatms),bsdy(1:ewld%bspline,1:config%mxatms),bsdz(1:ewld%bspline,1:config%mxatms), &
      Stat = fail(3))
    Allocate (bspx(1:ewld%bspline,1:config%mxatms),bspy(1:ewld%bspline,1:config%mxatms),bspz(1:ewld%bspline,1:config%mxatms), &
      Stat = fail(4))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_forces allocation failure'
      Call error(0,message)
    End If

    ! initialise coulombic potential energy and virial

    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    ! set working parameters

    rvolm=twopi/config%volm
    ralph=-0.25_wp/electro%alpha**2

    ! set scaling constant

    scale=rvolm*r4pie0/electro%eps

    ! Convert cell coordinates to fractional coordinates intervalled [0,1)
    ! (bottom left corner of MD cell) and stretch over kmaxs in different
    ! directions.  Only the halo (config%natms,config%nlast] has fractional coordinates
    ! outside the [0,1) interval.  In the worst case scenario of one
    ! "effective" link-cell per domain and one domain in the MD cell only,
    ! the halo will have fractional coordinates intervalled as
    ! [n,0)u[1,2), where -1 <= n < 0.  Only the positive halo is needed by
    ! the B-splines since they distribute/spread charge density in
    ! negative direction with length the length of the spline.
    !
    ! The story has become more complicated with cutoff padding and the
    ! conditional updates of the VNL and thus the halo as now a domain
    ! (1:config%natms) particle can enter the halo and vice versa.  So DD
    ! bounding is unsafe!!!

    Call invert(config%cell,rcell,det)
    If (Abs(det) < 1.0e-6_wp) Call error(120)

    Do i=1,config%nlast
      txx(i)=electro%kmaxa_r*(rcell(1)*config%parts(i)%xxx+rcell(4)*config%parts(i)%yyy+&
        rcell(7)*config%parts(i)%zzz+0.5_wp)
      tyy(i)=electro%kmaxb_r*(rcell(2)*config%parts(i)%xxx+rcell(5)*config%parts(i)%yyy+&
        rcell(8)*config%parts(i)%zzz+0.5_wp)
      tzz(i)=electro%kmaxc_r*(rcell(3)*config%parts(i)%xxx+rcell(6)*config%parts(i)%yyy+&
        rcell(9)*config%parts(i)%zzz+0.5_wp)

      ! If not DD bound in kmax grid space when .not.neigh%unconditional_update = (ewld%bspline1 == ewld%bspline)

      If (ewld%bspline1 == ewld%bspline .and. i <= config%natms) Then
        If (txx(i) < electro%ixbm1_r .or. txx(i) > electro%ixtm0_r .or. &
          tyy(i) < electro%iybm1_r .or. tyy(i) > electro%iytm0_r .or. &
          tzz(i) < electro%izbm1_r .or. tzz(i) > electro%iztm0_r) llspl=.false.
      End If

      ixx(i)=Int(txx(i))
      iyy(i)=Int(tyy(i))
      izz(i)=Int(tzz(i))

      ! Detect if a particle is charged and in the MD config%cell or in its positive halo
      ! (t(i) >= -zero_plus) as the B-splines are negative directionally by propagation

      If (tzz(i) >= -zero_plus .and. &
        tyy(i) >= -zero_plus .and. &
        txx(i) >= -zero_plus .and. &
        Abs(config%parts(i)%chge) > zero_plus) Then
        it(i)=1
      Else
        it(i)=0
      End If
    End Do

    ! Check for breakage of llspl when .not.neigh%unconditional_update = (ewld%bspline1 == ewld%bspline)

    ewld%bspline2=ewld%bspline1
    If (ewld%bspline1 == ewld%bspline) Then
      Call gcheck(comm,llspl)
      If (.not.llspl) ewld%bspline2=ewld%bspline+1
    End If

    ! construct B-splines for atoms

    Call bspgen(config,config%nlast,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,ewld,comm)

    Deallocate (txx,tyy,tzz, Stat = fail(1))
    If (fail(1) > 0) Then
      Write(message,'(a)') 'ewald_spme_forces allocation failure'
      Call error(0,message)
    End If

    ! zero 3D charge array
    ! DaFT version - only need set local bit to zero

    electro%qqc_local = 0.0_wp

    ! construct 3D charge array
    ! DaFT version - use array that holds only the local data

    Do i=1,config%nlast

      ! If a particle is charged and in the MD config%cell or in its positive halo
      ! (t(i) >= 0) as the B-splines are negative directionally by propagation

      If (it(i) == 1) Then
        bb3=config%parts(i)%chge

        !        Do l=1,ewld%bspline
        !           ll=izz(i)-l+2
        !
        !! If a particle's B-spline is entering this domain (originating from its
        !! positive halo), i.e. <= i.t, and not just to start exiting it, i.e. >= i.b
        !! In the limit of one domain in the MD config%cell (npr.=1, id.=0) i.t=kmax. and i.b=1
        !
        !           If (ll >= electro%izb .and. ll <= electro%izt) Then
        !              l_local = ll - electro%izb + 1
        !              bb2=bb3*bspz(l,i)
        !
        !              Do k=1,ewld%bspline
        !                 kk=iyy(i)-k+2
        !
        !                 If (kk >= electro%iyb .and. kk <= electro%iyt) Then
        !                    k_local = kk - electro%iyb + 1
        !                    bb1=bb2*bspy(k,i)
        !
        !                    Do j=1,ewld%bspline
        !                       jj=ixx(i)-j+2
        !
        !                       If (jj >= electro%ixb .and. jj <= electro%ixt) Then
        !                          j_local = jj - electro%ixb + 1
        !
        !                          det=bb1*bspx(j,i)
        !
        !                          electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
        !                       End If
        !                    End Do
        !                 End If
        !              End Do
        !           End If
        !        End Do

        llb = Max( electro%izb, izz(i) - ewld%bspline + 2 )
        llt = Min( electro%izt, izz(i) + 1 )

        kkb = Max( electro%iyb, iyy(i) - ewld%bspline + 2 )
        kkt = Min( electro%iyt, iyy(i) + 1 )

        jjb = Max( electro%ixb, ixx(i) - ewld%bspline + 2 )
        jjt = Min( electro%ixt, ixx(i) + 1 )

        Select Case( jjt - jjb + 1 )

        Case Default

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              Do jj = jjb, jjt
                j = ixx(i) - jj + 2

                j_local = jj - electro%ixb + 1
                det=bb1*bspx(j,i)

                electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
              End Do
            End Do
          End Do

        Case( 0 )

        Case( 1 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 2 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !2
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 3 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !2
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !3
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 4 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !2
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !3
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !4
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 5 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !2
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !3
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !4
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !5
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 6 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !2
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !3
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !4
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !5
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !6
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 7 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              !1
              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !2
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !3
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !4
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !5
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !6
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              !7
              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        Case( 8 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb + 1
            bb2=bb3*bspz(l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb + 1
              bb1=bb2*bspy(k,i)

              jj = jjb

              j = ixx(i) - jj + 2
              j_local = jj - electro%ixb + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det

              j = j - 1
              j_local = j_local + 1
              det=bb1*bspx(j,i)

              electro%qqc_local(j_local,k_local,l_local)=electro%qqc_local(j_local,k_local,l_local)+det
            End Do
          End Do

        End Select
      End If

    End Do

    ! load charge array into complex array for FFT

    electro%qqq_local=Cmplx(electro%qqc_local , Kind = wp)

    ! calculate inverse 3D FFT of charge array (in place)

    Call pfft(electro%qqq_local,electro%pfft_work,electro%context,1)

    ! set reciprocal space cutoff

    Call dcell(rcell,celprp)

    rcpcut=0.5_wp*Min(electro%kmaxa_r*celprp(7),electro%kmaxb_r*celprp(8),electro%kmaxc_r*celprp(9))
    rcpcut=rcpcut*1.05_wp*twopi
    rcpct2=rcpcut**2

    ! initialise temporary stress tensor

    strs = 0.0_wp

    ! calculate convolution of charge array with gaussian function
    ! DaFT Version - only loop over the local stuff

    Do l_local=1,electro%block_z
      l=electro%index_z(l_local)

      ll=l-1
      If (l > ewld%fft_dim_c/2) ll=ll-ewld%fft_dim_c
      tmp=twopi*Real(ll,wp)

      rkx1=tmp*rcell(3)
      rky1=tmp*rcell(6)
      rkz1=tmp*rcell(9)

      bb3=Real( electro%bscz(l)*Conjg(electro%bscz(l)),wp )

      Do k_local=1,electro%block_y
        k=electro%index_y(k_local)

        kk=k-1
        If (k > ewld%fft_dim_b/2) kk=kk-ewld%fft_dim_b
        tmp=twopi*Real(kk,wp)

        rkx2=rkx1+tmp*rcell(2)
        rky2=rky1+tmp*rcell(5)
        rkz2=rkz1+tmp*rcell(8)

        bb2=bb3*Real( electro%bscy(k)*Conjg(electro%bscy(k)),wp )

        Do j_local=1,electro%block_x
          j=electro%index_x(j_local)

          jj=j-1
          If (j > ewld%fft_dim_a/2) jj=jj-ewld%fft_dim_a
          tmp=twopi*Real(jj,wp)

          rkx3=rkx2+tmp*rcell(1)
          rky3=rky2+tmp*rcell(4)
          rkz3=rkz2+tmp*rcell(7)

          bb1=bb2*Real( electro%bscx(j)*Conjg(electro%bscx(j)),wp )

          rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

          If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

            vterm=bb1*Exp(ralph*rksq)/rksq*electro%qqq_local(j_local,k_local,l_local)
            akv=2.0_wp*(1.0_wp/rksq-ralph)*Real( vterm*Conjg(electro%qqq_local(j_local,k_local,l_local)),wp )

            strs(1)=strs(1)-rkx3*rkx3*akv
            strs(5)=strs(5)-rky3*rky3*akv
            strs(9)=strs(9)-rkz3*rkz3*akv
            strs(2)=strs(2)-rkx3*rky3*akv
            strs(3)=strs(3)-rkx3*rkz3*akv
            strs(6)=strs(6)-rky3*rkz3*akv

            electro%qqq_local(j_local,k_local,l_local)=vterm

          Else

            electro%qqq_local(j_local,k_local,l_local)=(0.0_wp,0.0_wp)

          End If
        End Do
      End Do
    End Do

    ! complete strs

    strs(4) = strs(2)
    strs(7) = strs(3)
    strs(8) = strs(6)

    ! as only looped over local stuff, we need to gsum strs

    Call gsum(comm,strs)

    ! scale strs and distribute per node

    strs = strs * scale / Real(comm%mxnode,wp)

    ! calculate atomic energy

    Call pfft(electro%qqq_local,electro%pfft_work,electro%context,-1)

    eng = 0.0_wp
    Do l=1,electro%block_z
      Do k=1,electro%block_y
        Do j=1,electro%block_x
          qqc_tmp=Real(electro%qqq_local(j,k,l),wp)
          eng=eng+electro%qqc_local(j,k,l)*qqc_tmp
          electro%qqc_local(j,k,l)=qqc_tmp
        End Do
      End Do
    End Do

    ! as only looped over local stuff, we need to gsum the eng

    Call gsum(comm,eng)

    ! scale eng and distribute per node

    eng = eng * scale / Real(comm%mxnode,wp)

    ! Second part of the monopole contribution to the stress tensor
    ! calculate stress tensor (symmetrical, per node)

    strs   = strs + eng*uni
    stress = stress + strs

    ! distribute energy and virial terms (per node)

    engcpe_rc = eng + ewld%engsic / Real(comm%mxnode,wp)
    vircpe_rc = -(strs(1)+strs(5)+strs(9))

    ! infrequent calculations copying

    If (ewld%l_cp) Then
      ewld%e_rc=engcpe_rc
      ewld%v_rc=vircpe_rc
      ewld%s_rc=strs
    End If

    ! calculate atomic forces

    Call spme_forces(rcell,scale, ixx,iyy,izz, bspx,bspy,bspz, bsdx,bsdy,bsdz, electro%qqc_local, &
      electro%ixb,electro%ixt, electro%iyb,electro%iyt, electro%izb,electro%izt)

    Deallocate (ixx,iyy,izz,it, Stat = fail(1))
    Deallocate (bsdx,bsdy,bsdz, Stat = fail(2))
    Deallocate (bspx,bspy,bspz, Stat = fail(3))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_forces deallocation failure'
      Call error(0,message)
    End If

  Contains

    Subroutine spme_forces(rcell,scale, ixx,iyy,izz, bspx,bspy,bspz, bsdx,bsdy,bsdz, qqc_local, ixb,ixt, iyb,iyt, izb,izt)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! dl_poly_4 subroutine for calculating coulombic forces in a periodic
      ! system using smooth particle mesh ewald method (fourier part)
      !
      ! Note: electro%qqc_local is shifted from its definition from above
      !       and therefore there is no need for periodic images (!!)
      !
      ! copyright - daresbury laboratory
      ! author    - w.smith & i.t.todorov june 2014
      ! refactoring:
      !           - a.m.elena march-october 2018
      !           - j.madge march-october 2018
      !           - a.b.g.chalk march-october 2018
      !           - i.scivetti march-october 2018
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      Integer,           Intent( In    ) :: ixx(1:config%mxatms),iyy(1:config%mxatms),izz(1:config%mxatms), &
        ixb,ixt, iyb,iyt, izb,izt
      Real( Kind = wp ), Intent( In    ) :: scale,rcell(1:9),                   &
        bsdx(1:ewld%bspline,1:config%mxatms),bsdy(1:ewld%bspline,1:config%mxatms),bsdz(1:ewld%bspline,1:config%mxatms), &
        bspx(1:ewld%bspline,1:config%mxatms),bspy(1:ewld%bspline,1:config%mxatms),bspz(1:ewld%bspline,1:config%mxatms), &
        qqc_local( ixb:ixt, iyb:iyt, izb:izt )

      Integer           :: delspl, ixdb,ixdt,iydb,iydt,izdb,izdt, fail, i,j,k,l, jj,kk,ll
      Real( Kind = wp ) :: tmp,facx,facy,facz,fff(0:3),fx,fy,fz,fix,fiy,fiz,qsum, &
        bdxl,bdyl,bdzl,bdxk,bdyk,bdzk,bdxj,bdyj,bdzj

      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: qqc_domain

      Character ( Len = 256 ) :: message

      ! Define extended ranges for the domain = local + halo slice and allocate

      ixdb = ixb - ewld%bspline2
      iydb = iyb - ewld%bspline2
      izdb = izb - ewld%bspline2

      delspl = ewld%bspline2 - ewld%bspline

      ixdt = ixt + delspl
      iydt = iyt + delspl
      izdt = izt + delspl

      fail=0
      Allocate (qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt ), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'spme_forces allocation failure'
        call error(0,message) 
      End If

      Call exchange_grid(ixb , ixt , iyb , iyt , izb , izt , electro%qqc_local, &
        ixdb, iydb, izdb, ixdt, iydt, izdt, qqc_domain, domain, ewld, comm)

      ! Real values of kmax vectors

      If (electro%newjob_fspme) Then
        electro%newjob_fspme = .false.

        electro%kmaxa_lr=Real(ewld%fft_dim_a,wp)
        electro%kmaxb_lr=Real(ewld%fft_dim_b,wp)
        electro%kmaxc_lr=Real(ewld%fft_dim_c,wp)
      End If

      tmp=-2.0_wp*scale
      facx=tmp*electro%kmaxa_lr
      facy=tmp*electro%kmaxb_lr
      facz=tmp*electro%kmaxc_lr

      fff=0.0_wp
      Do i=1,config%natms
        tmp=config%parts(i)%chge

        If (Abs(tmp) > zero_plus) Then

          ! initialise forces

          fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp

          Do l=1,ewld%bspline
            ll=izz(i)-l+2

            bdxl=tmp*facx*bspz(l,i)
            bdyl=tmp*facy*bspz(l,i)
            bdzl=tmp*facz*bsdz(l,i)

            Do k=1,ewld%bspline
              kk=iyy(i)-k+2

              bdxk=bdxl*bspy(k,i)
              bdyk=bdyl*bsdy(k,i)
              bdzk=bdzl*bspy(k,i)

              Do j=1,ewld%bspline
                jj=ixx(i)-j+2

                qsum=qqc_domain(jj,kk,ll)

                bdxj=qsum*bdxk*bsdx(j,i)
                bdyj=qsum*bdyk*bspx(j,i)
                bdzj=qsum*bdzk*bspx(j,i)

                fix=fix+bdxj
                fiy=fiy+bdyj
                fiz=fiz+bdzj
              End Do
            End Do
          End Do

          fx=fix*rcell(1)+fiy*rcell(2)+fiz*rcell(3)
          fy=fix*rcell(4)+fiy*rcell(5)+fiz*rcell(6)
          fz=fix*rcell(7)+fiy*rcell(8)+fiz*rcell(9)

          ! accumulate forces

          fff(0)=fff(0)+1.0_wp
          fff(1)=fff(1)+fx
          fff(2)=fff(2)+fy
          fff(3)=fff(3)+fz

          ! load forces

          config%parts(i)%fxx=config%parts(i)%fxx+fx
          config%parts(i)%fyy=config%parts(i)%fyy+fy
          config%parts(i)%fzz=config%parts(i)%fzz+fz

          ! infrequent calculations copying

          If (ewld%l_cp) Then
            ewld%fcx(i)=ewld%fcx(i)+fx
            ewld%fcy(i)=ewld%fcy(i)+fy
            ewld%fcz(i)=ewld%fcz(i)+fz
          End If

        End If
      End Do

      ! remove COM drift arising from SPME approximations

      Call gsum(comm,fff)
      If (fff(0) > zero_plus) Then
        fff(1:3)=fff(1:3)/fff(0)

        Do i=1,config%natms
          If (Abs(config%parts(i)%chge) > zero_plus) Then

            config%parts(i)%fxx=config%parts(i)%fxx-fff(1)
            config%parts(i)%fyy=config%parts(i)%fyy-fff(2)
            config%parts(i)%fzz=config%parts(i)%fzz-fff(3)

            ! infrequent calculations copying

            If (ewld%l_cp) Then
              ewld%fcx(i)=ewld%fcx(i)-fff(1)
              ewld%fcy(i)=ewld%fcy(i)-fff(2)
              ewld%fcz(i)=ewld%fcz(i)-fff(3)
            End If
          End If
        End Do
      End If

      Deallocate (qqc_domain, Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'spme_forces dealocation failure'
        call error(0,message)
      End If

    End Subroutine spme_forces

  End Subroutine ewald_spme_forces

  Subroutine ewald_excl_forces(iatm,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex,stress, &
      neigh,electro,config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating coulombic energy and force terms
    ! in a periodic system using ewald's method
    !
    ! Note: exclusion correction terms
    !       frozen pairs are ignored by default, they are not dealt with here
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( configuration_type ),               Intent( InOut ) :: config

    Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp =  0.3275911_wp
    Real( Kind = wp ), Parameter :: rr3  = 1.0_wp/3.0_wp
    Real( Kind = wp ), Parameter :: r10  = 0.1_wp
    Real( Kind = wp ), Parameter :: r42  = 1.0_wp/42.0_wp
    Real( Kind = wp ), Parameter :: r216 = 1.0_wp/216.0_wp

    Integer           :: limit,idi,jatm,m
    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,alpr,alpr2, &
      erfr,egamma,exp1,tt,             &
      fix,fiy,fiz,fx,fy,fz,            &
      strs1,strs2,strs3,strs5,strs6,strs9

    ! initialise potential energy and virial

    engcpe_ex=0.0_wp
    vircpe_ex=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/electro%eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! Get neigh%list limit

      limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

      ! start of primary loop for forces evaluation

      Do m=1,limit

        ! atomic index and charge

        jatm=neigh%list(neigh%list(0,iatm)+m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! Squared distance

          rsq=rrr**2

          ! calculate forces

          alpr =rrr*electro%alpha
          alpr2=alpr*alpr

          ! calculate error function and derivative

          If (alpr < 1.0e-2_wp) Then

            ! close particles (core-shell units) - small distances limit

            erfr=2.0_wp*chgprd*(electro%alpha/sqrpi) * &
              (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

            egamma=-4.0_wp*chgprd*(electro%alpha**3/sqrpi) * &
              (rr3+alpr2*(-2.0_wp*r10+alpr2*(3.0_wp*r42-4.0_wp*alpr2*r216)))

          Else

            ! distant particles - traditional

            exp1=Exp(-(electro%alpha*rrr)**2)
            tt  =1.0_wp/(1.0_wp+pp*electro%alpha*rrr)

            erfr=chgprd * &
              (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

            egamma=-(erfr-2.0_wp*chgprd*(electro%alpha/sqrpi)*exp1)/rsq

          End If

          ! calculate forces

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! add potential energy and virial

            engcpe_ex = engcpe_ex - erfr
            vircpe_ex = vircpe_ex - egamma*rsq

            ! add stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

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

    End If

  End Subroutine ewald_excl_forces

  Subroutine ewald_frzn_forces(engcpe_fr,vircpe_fr,stress,ewld,neigh,electro,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating corrections to coulombic energy
    ! and forces in a periodic system arising from frozen pairs
    !
    ! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
    !       end (and any COM drift removed) but corrections to the stress
    !       and the virial are important as they feed into the system
    !       pressure response.  Constant volume ensembles (ensemble < 20)
    !       need this calculation just once! - controlled by ewld%lf_fce in
    !       ewald_check<-two_body_forces
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp  ),                   Intent(   Out ) :: engcpe_fr,vircpe_fr
    Real( Kind = wp  ), Dimension( 1:9 ), Intent( InOut ) :: stress
    Type( ewald_type ),                   Intent( InOut ) :: ewld
    Type( neighbours_type ),              Intent( In    ) :: neigh
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( comms_type ),                   Intent( InOut ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config

    Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

    Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
    Real( Kind = wp ) :: scl,det,rcell(1:9),xrr,yrr,zrr,rrr,rsq, &
      chgprd,erfr,egamma,exp1,tt,             &
      fx,fy,fz,xss,yss,zss,                   &
      strs1,strs2,strs3,strs5,strs6,strs9

    Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr
    Real( Kind = wp ), Dimension( : ), Allocatable :: cfr,xfr,yfr,zfr
    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt
    Character( Len = 256 ) :: message

    If (.not.ewld%lf_fce) Then ! All's been done but needs copying
      Do i=1,config%natms
        config%parts(i)%fxx=config%parts(i)%fxx+ewld%ffx(i)
        config%parts(i)%fyy=config%parts(i)%fyy+ewld%ffy(i)
        config%parts(i)%fzz=config%parts(i)%fzz+ewld%ffz(i)
      End Do

      engcpe_fr=ewld%ef_fr
      vircpe_fr=ewld%vf_fr
      stress=stress+ewld%sf_fr

      If (ewld%l_cp) Then
        Do i=1,config%natms
          ewld%fcx(i)=ewld%fcx(i)+ewld%ffx(i)
          ewld%fcy(i)=ewld%fcy(i)+ewld%ffy(i)
          ewld%fcz(i)=ewld%fcz(i)+ewld%ffz(i)
        End Do

        ewld%e_fr=ewld%ef_fr
        ewld%v_fr=ewld%vf_fr
        ewld%s_fr=ewld%sf_fr
      End If

      Return
    End If

    fail=0
    Allocate (l_ind(1:config%mxatdm),nz_fr(0:comm%mxnode), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'ewald_frzn_forces allocation failure'
      Call error(0,message)
    End If

    Call invert(config%cell,rcell,det)

    ! Initialise contributions

    engcpe_fr=0.0_wp
    vircpe_fr=0.0_wp

    strs1 = 0.0_wp
    strs2 = 0.0_wp
    strs3 = 0.0_wp
    strs5 = 0.0_wp
    strs6 = 0.0_wp
    strs9 = 0.0_wp

    l_ind=0 ; nz_fr=0
    Do i=1,config%natms
      If (config%lfrzn(i) > 0 .and. Abs(config%parts(i)%chge) > zero_plus) Then
        nz_fr(comm%idnode+1)=nz_fr(comm%idnode+1)+1
        l_ind(nz_fr(comm%idnode+1))=i
      End If
    End Do
    Call gsum(comm, nz_fr)
    nz_fr(0) = Sum(nz_fr(0:comm%idnode)) ! Offset

    scl=r4pie0/electro%eps
    nzfr = Sum(nz_fr(1:comm%mxnode))     ! Total
    If (nzfr <= 10*config%mxatms) Then

      Allocate (cfr(1:nzfr),xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a,i0)') 'ewald_frzn_forces allocation failure 1'
        Call error(0,message)
      End If

      cfr=0.0_wp
      xfr=0.0_wp
      yfr=0.0_wp
      zfr=0.0_wp
      Do i=1,nz_fr(comm%idnode+1)
        ii=nz_fr(0)+i

        cfr(ii)=config%parts(l_ind(i))%chge
        xfr(ii)=config%parts(l_ind(i))%xxx
        yfr(ii)=config%parts(l_ind(i))%yyy
        zfr(ii)=config%parts(l_ind(i))%zzz
      End Do
      Call gsum(comm, cfr)
      Call gsum(comm, xfr)
      Call gsum(comm, yfr)
      Call gsum(comm, zfr)

      Do i=1,nz_fr(comm%idnode+1)
        ii=nz_fr(0)+i

        Do jj=1,nz_fr(0) ! -, on nodes<comm%idnode
          xrr=xfr(ii)-xfr(jj)
          yrr=yfr(ii)-yfr(jj)
          zrr=zfr(ii)-zfr(jj)

          xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
          yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
          zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

          xss=xss-Anint(xss)
          yss=yss-Anint(yss)
          zss=zss-Anint(zss)

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          rsq=xrr**2+yrr**2+zrr**2

          rrr=Sqrt(rsq)
          chgprd=cfr(ii)*cfr(jj)*scl

          ! calculate error function and derivative

          exp1 =Exp(-(electro%alpha*rrr)**2)
          tt   =1.0_wp/(1.0_wp+pp*electro%alpha*rrr)

          erfr=chgprd * &
            (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

          egamma=-(erfr-2.0_wp*chgprd*(electro%alpha/sqrpi)*exp1)/rsq

          fx = egamma*xrr
          fy = egamma*yrr
          fz = egamma*zrr

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-fx
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-fy
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-fz

          ! redundant calculations copying

          If (ewld%lf_cp) Then
            ewld%ffx(l_ind(i))=ewld%ffx(l_ind(i))-fx
            ewld%ffy(l_ind(i))=ewld%ffy(l_ind(i))-fy
            ewld%ffz(l_ind(i))=ewld%ffz(l_ind(i))-fz
          End If

          ! infrequent calculations copying

          If (ewld%l_cp) Then
            ewld%fcx(l_ind(i))=ewld%fcx(l_ind(i))-fx
            ewld%fcy(l_ind(i))=ewld%fcy(l_ind(i))-fy
            ewld%fcz(l_ind(i))=ewld%fcz(l_ind(i))-fz
          End If
        End Do

        Do j=i+1,nz_fr(comm%idnode+1) ! =, node=comm%idnode (OVERLAP but no SELF)!
          jj=nz_fr(0)+j

          xrr=xfr(ii)-xfr(jj)
          yrr=yfr(ii)-yfr(jj)
          zrr=zfr(ii)-zfr(jj)

          xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
          yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
          zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

          xss=xss-Anint(xss)
          yss=yss-Anint(yss)
          zss=zss-Anint(zss)

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          rsq=xrr**2+yrr**2+zrr**2

          rrr=Sqrt(rsq)
          chgprd=cfr(ii)*cfr(jj)*scl

          ! calculate error function and derivative

          exp1 =Exp(-(electro%alpha*rrr)**2)
          tt   =1.0_wp/(1.0_wp+pp*electro%alpha*rrr)

          erfr=chgprd * &
            (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

          egamma=-(erfr-2.0_wp*chgprd*(electro%alpha/sqrpi)*exp1)/rsq

          fx = egamma*xrr
          fy = egamma*yrr
          fz = egamma*zrr

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-fx
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-fy
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-fz

          config%parts(l_ind(j))%fxx=config%parts(l_ind(j))%fxx+fx
          config%parts(l_ind(j))%fyy=config%parts(l_ind(j))%fyy+fy
          config%parts(l_ind(j))%fzz=config%parts(l_ind(j))%fzz+fz

          ! redundant calculations copying

          If (ewld%lf_cp) Then
            ewld%ffx(l_ind(i))=ewld%ffx(l_ind(i))-fx
            ewld%ffy(l_ind(i))=ewld%ffy(l_ind(i))-fy
            ewld%ffz(l_ind(i))=ewld%ffz(l_ind(i))-fz

            ewld%ffx(l_ind(j))=ewld%ffx(l_ind(j))+fx
            ewld%ffy(l_ind(j))=ewld%ffy(l_ind(j))+fy
            ewld%ffz(l_ind(j))=ewld%ffz(l_ind(j))+fz
          End If

          ! infrequent calculations copying

          If (ewld%l_cp) Then
            ewld%fcx(l_ind(i))=ewld%fcx(l_ind(i))-fx
            ewld%fcy(l_ind(i))=ewld%fcy(l_ind(i))-fy
            ewld%fcz(l_ind(i))=ewld%fcz(l_ind(i))-fz

            ewld%fcx(l_ind(j))=ewld%fcx(l_ind(j))+fx
            ewld%fcy(l_ind(j))=ewld%fcy(l_ind(j))+fy
            ewld%fcz(l_ind(j))=ewld%fcz(l_ind(j))+fz
          End If

          ! calculate potential energy and virial

          engcpe_fr = engcpe_fr - erfr
          vircpe_fr = vircpe_fr - egamma*rsq

          ! calculate stress tensor

          strs1 = strs1 + xrr*fx
          strs2 = strs2 + xrr*fy
          strs3 = strs3 + xrr*fz
          strs5 = strs5 + yrr*fy
          strs6 = strs6 + yrr*fz
          strs9 = strs9 + zrr*fz
        End Do

        Do jj=nz_fr(0)+nz_fr(comm%idnode+1)+1,nzfr ! +, on nodes>comm%idnode
          xrr=xfr(ii)-xfr(jj)
          yrr=yfr(ii)-yfr(jj)
          zrr=zfr(ii)-zfr(jj)

          xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
          yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
          zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

          xss=xss-Anint(xss)
          yss=yss-Anint(yss)
          zss=zss-Anint(zss)

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          rsq=xrr**2+yrr**2+zrr**2

          rrr=Sqrt(rsq)
          chgprd=cfr(ii)*cfr(jj)*scl

          ! calculate error function and derivative

          exp1 =Exp(-(electro%alpha*rrr)**2)
          tt   =1.0_wp/(1.0_wp+pp*electro%alpha*rrr)

          erfr=chgprd * &
            (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

          egamma=-(erfr-2.0_wp*chgprd*(electro%alpha/sqrpi)*exp1)/rsq

          fx = egamma*xrr
          fy = egamma*yrr
          fz = egamma*zrr

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-fx
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-fy
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-fz

          ! redundant calculations copying

          If (ewld%lf_cp) Then
            ewld%ffx(l_ind(i))=ewld%ffx(l_ind(i))-fx
            ewld%ffy(l_ind(i))=ewld%ffy(l_ind(i))-fy
            ewld%ffz(l_ind(i))=ewld%ffz(l_ind(i))-fz
          End If

          ! infrequent calculations copying

          If (ewld%l_cp) Then
            ewld%fcx(l_ind(i))=ewld%fcx(l_ind(i))-fx
            ewld%fcy(l_ind(i))=ewld%fcy(l_ind(i))-fy
            ewld%fcz(l_ind(i))=ewld%fcz(l_ind(i))-fz
          End If

          ! calculate potential energy and virial

          engcpe_fr = engcpe_fr - erfr
          vircpe_fr = vircpe_fr - egamma*rsq

          ! calculate stress tensor

          strs1 = strs1 + xrr*fx
          strs2 = strs2 + xrr*fy
          strs3 = strs3 + xrr*fz
          strs5 = strs5 + yrr*fy
          strs6 = strs6 + yrr*fz
          strs9 = strs9 + zrr*fz
        End Do
      End Do

      Deallocate (cfr,xfr,yfr,zfr, Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_forces deallocation failure 1'
        Call error(0,message)
      End If

    Else

      ! We resort to approximating N*(N-1)/2 interactions
      ! with the short-range one from the two body linked config%cell neigh%list

      Allocate (xxt(1:neigh%max_list),yyt(1:neigh%max_list),zzt(1:neigh%max_list),rrt(1:neigh%max_list), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_forces allocation failure 2'
        Call error(0,message)
      End If

      Do ii=1,nz_fr(comm%idnode+1)
        i=l_ind(nz_fr(comm%idnode+1))
        idi=config%ltg(ii)

        ! Get neigh%list limit

        limit=neigh%list(-2,i)-neigh%list(-1,i)
        If (limit > 0) Then

          ! calculate interatomic distances

          Do k=1,limit
            j=neigh%list(neigh%list(-1,i)+k,i)

            xxt(k)=config%parts(i)%xxx-config%parts(j)%xxx
            yyt(k)=config%parts(i)%yyy-config%parts(j)%yyy
            zzt(k)=config%parts(i)%zzz-config%parts(j)%zzz
          End Do

          ! periodic boundary conditions not needed by LC construction
          !
          !           Call images(config%imcon,config%cell,limit,xxt,yyt,zzt)

          ! square of distances

          Do k=1,limit
            rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
          End Do

          Do k=1,limit
            j=neigh%list(neigh%list(-1,i)+k,i)

            rrr=rrt(k)
            If (Abs(config%parts(j)%chge) > zero_plus .and. rrr < neigh%cutoff) Then
              chgprd=config%parts(i)%chge*config%parts(j)%chge*scl
              rsq=rrr**2

              ! calculate error function and derivative

              exp1 =Exp(-(electro%alpha*rrr)**2)
              tt   =1.0_wp/(1.0_wp+pp*electro%alpha*rrr)

              erfr=chgprd * &
                (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

              egamma=-(erfr-2.0_wp*chgprd*(electro%alpha/sqrpi)*exp1)/rsq

              fx = egamma*xxt(k)
              fy = egamma*yyt(k)
              fz = egamma*zzt(k)

              ! calculate forces

              config%parts(i)%fxx=config%parts(i)%fxx-fx
              config%parts(i)%fyy=config%parts(i)%fyy-fy
              config%parts(i)%fzz=config%parts(i)%fzz-fz

              ! redundant calculations copying

              If (ewld%lf_cp) Then
                ewld%ffx(i)=ewld%ffx(i)-fx
                ewld%ffy(i)=ewld%ffy(i)-fy
                ewld%ffz(i)=ewld%ffz(i)-fz
              End If

              ! infrequent calculations copying

              If (ewld%l_cp) Then
                ewld%fcx(i)=ewld%fcx(i)-fx
                ewld%fcy(i)=ewld%fcy(i)-fy
                ewld%fcz(i)=ewld%fcz(i)-fz
              End If

              If (j <= config%natms) Then

                config%parts(j)%fxx=config%parts(j)%fxx+fx
                config%parts(j)%fyy=config%parts(j)%fyy+fy
                config%parts(j)%fzz=config%parts(j)%fzz+fz

                ! redundant calculations copying

                If (ewld%lf_cp) Then
                  ewld%ffx(j)=ewld%ffx(j)+fx
                  ewld%ffy(j)=ewld%ffy(j)+fy
                  ewld%ffz(j)=ewld%ffz(j)+fz
                End If

                ! infrequent calculations copying

                If (ewld%l_cp) Then
                  ewld%fcx(j)=ewld%fcx(j)+fx
                  ewld%fcy(j)=ewld%fcy(j)+fy
                  ewld%fcz(j)=ewld%fcz(j)+fz
                End If

              End If

              If (j <= config%natms .or. idi < config%ltg(j)) Then

                ! calculate potential energy and virial

                engcpe_fr = engcpe_fr - erfr
                vircpe_fr = vircpe_fr - egamma*rsq

                ! calculate stress tensor

                strs1 = strs1 + xxt(k)*fx
                strs2 = strs2 + xxt(k)*fy
                strs3 = strs3 + xxt(k)*fz
                strs5 = strs5 + yyt(k)*fy
                strs6 = strs6 + yyt(k)*fz
                strs9 = strs9 + zzt(k)*fz

              End If
            End If
          End Do

        End If
      End Do

      Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_forces deallocation failure 2'
        Call error(0,message)
      End If

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

    ! redundant calculations copying

    If (ewld%lf_cp) Then
      ewld%ef_fr=engcpe_fr
      ewld%vf_fr=vircpe_fr

      ewld%sf_fr(1) = strs1
      ewld%sf_fr(2) = strs2
      ewld%sf_fr(3) = strs3
      ewld%sf_fr(4) = strs2
      ewld%sf_fr(5) = strs5
      ewld%sf_fr(6) = strs6
      ewld%sf_fr(7) = strs3
      ewld%sf_fr(8) = strs6
      ewld%sf_fr(9) = strs9
    End If

    ! infrequent calculations copying

    If (ewld%l_cp) Then
      ewld%e_fr=engcpe_fr
      ewld%v_fr=vircpe_fr

      ewld%s_fr(1) = strs1
      ewld%s_fr(2) = strs2
      ewld%s_fr(3) = strs3
      ewld%s_fr(4) = strs2
      ewld%s_fr(5) = strs5
      ewld%s_fr(6) = strs6
      ewld%s_fr(7) = strs3
      ewld%s_fr(8) = strs6
      ewld%s_fr(9) = strs9
    End If

    Deallocate (l_ind,nz_fr, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'ewald_frzn_forces deallocation failure'
      Call error(0,message)
    End If
  End Subroutine ewald_frzn_forces
End Module ewald_spole
