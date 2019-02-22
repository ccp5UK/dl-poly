Module ewald_mpole
  !! This module has no header !
  use bspline,          Only : bspline_type, bspline_coeffs_gen, bspline_splines_gen
  Use comms,            Only : comms_type, gcheck, gsum
  Use configuration,    Only : configuration_type
  Use domains,          Only : domains_type,exchange_grid
  Use electrostatic,    Only : electrostatic_type
  Use errors_warnings,  Only : error, error_alloc, error_dealloc
  Use ewald,            Only : ewald_type,dtpbsp
  Use kinds,            Only : wp, sp
  Use mpole,            Only : mpole_type
  Use mpoles_container, Only : ewald_deriv,explicit_ewald_real_loops,explicit_spme_loops,limit_erfr_deriv
  Use neighbours,       Only : neighbours_type
  Use numerics,         Only : invert,dcell,erfcgen, erfc, erfc_deriv
  Use parallel_fft,     Only : initialize_fft, pfft, pfft_indices
  Use particle,         Only : corePart
  Use constants,        Only : sqrpi,zero_plus,twopi, r4pie0
  Implicit None

  Private

  Public :: ewald_real_mforces, ewald_real_mforces_d, ewald_spme_mforces, ewald_spme_mforces_d, &
    ewald_excl_mforces, ewald_excl_mforces_d,ewald_frzn_mforces

Contains

  Subroutine ewald_real_mforces(iatm,xxt,yyt,zzt,rrt,engcpe_rl,vircpe_rl,stress, &
    ewld,neigh,mpoles,electro,domain,config,comm)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using multipoles with the ewald real space
    !! kernel
    !!
    !! copyright - daresbury laboratory
    !! author    - h.a.boateng & i.t.todorov february 2016
    !!
    !!----------------------------------------------------------------------!

    Integer,                                          Intent( In    ) :: iatm
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe_rl,vircpe_rl
    Real( Kind = wp ), Dimension( 1:9 ),              Intent( InOut ) :: stress
    Type( mpole_type ),                               Intent( InOut ) :: mpoles
    Type( electrostatic_type ),                       Intent( InOut ) :: electro
    Type( ewald_type ),                          Intent( In    ) :: ewld
    Type( domains_type ),                             Intent( In    ) :: domain
    Type( comms_type ),                               Intent( In    ) :: comm
    Type( configuration_type ),                       Intent( InOut ) :: config


    Integer           :: fail,idi,jatm,k1,k2,k3,s1,s2,s3,m,n, &
      k,ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,rrr,engmpl,fix,fiy,fiz,fx,fy,fz, &
      strs1,strs2,strs3,strs5,strs6,strs9, &
      ppp,vk0,vk1,vk2,t1,t2,kx,ky,kz,      &
      txyz,erfcr,tmp,tmpi,tmpj,tix,        &
      alphan,tiy,tiz,tjx,tjy,tjz,sx,sy,sz


    Real( Kind = wp ) :: d1(0:2*mpoles%max_order+1,0:2*mpoles%max_order+1,0:2*mpoles%max_order+1)
    Real( Kind = wp ) :: imp(1:mpoles%max_mpoles),jmp(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: jmpx(1:mpoles%max_mpoles),jmpy(1:mpoles%max_mpoles),jmpz(1:mpoles%max_mpoles)
    Logical, save :: newjob
    Character( Len = 256 ) :: message

    If (newjob) Then
      newjob = .false.

      Call erfcgen(neigh%cutoff,ewld%alpha,erfc,erfc_deriv)
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

    ! get the multipoles for site i

    imp=mpoles%global_frame(:,iatm)

    ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

      ! get the components for site i infinitesimal rotations

      impx=mpoles%rotation_x(:,iatm)
      impy=mpoles%rotation_y(:,iatm)
      impz=mpoles%rotation_z(:,iatm)

      ! multipole scaler

      scl=2.0_wp*ewld%alpha*r4pie0/(sqrpi*electro%eps)

      ! scale imp multipoles

      imp=imp*scl

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! initialize torques for atom i (temporary)

      tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index

        jatm=neigh%list(m,iatm)

        ! get the multipoles for site j

        jmp=mpoles%global_frame(:,jatm)

        ! interatomic distance

        rrr = rrt(m)

        ! truncation of potential

        If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

          ! get the components for site j infinitesimal rotations

          jmpx=mpoles%rotation_x(:,jatm)
          jmpy=mpoles%rotation_y(:,jatm)
          jmpz=mpoles%rotation_z(:,jatm)

          ! get the value of the kernel using 3pt interpolation

          ! erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/ewld%alpha
          erfcr = three_p_interp(erfc,rrr)/ewld%alpha

          ! compute derivatives of kernel

          Call ewald_deriv(0,2*mpoles%max_order+1,1,erfcr,ewld%alpha*xxt(m), &
            ewld%alpha*yyt(m),ewld%alpha*zzt(m),ewld%alpha*rrr,mpoles%max_order,d1)

          ! calculate forces

          engmpl = 0.0_wp
          fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
          tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

          If (mpoles%max_order < 5) Then

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  jj = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(jj)) > zero_plus) Call explicit_ewald_real_loops &
                    ( 0,2*mpoles%max_order+1, k1,k2,k3, ewld%alpha, d1,               &
                    imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                    kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                    engmpl,fx,fy,fz,mpoles)

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          Else

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  jj=mpoles%map(k1,k2,k3)

                  If (Abs(jmp(jj)) > zero_plus) Then

                    txyz=kx*jmp(jj)

                    sz = 1.0_wp
                    Do s3=0,mpoles%max_order
                      ks3=k3+s3; ks31=ks3+1

                      sy = sz
                      Do s2=0,mpoles%max_order-s3
                        ks2=k2+s2; ks21=ks2+1

                        sx = sy
                        Do s1=0,mpoles%max_order-s3-s2
                          ks1=k1+s1; ks11=ks1+1

                          n      = ks1+ks2+ks3
                          alphan = ewld%alpha**n

                          ii     = mpoles%map(s1,s2,s3)

                          tmp    = alphan*d1(ks1,ks2,ks3)

                          tmpi   = txyz       * tmp
                          tmpj   = sx*imp(ii) * tmp

                          t1     = alphan     * txyz*imp(ii)

                          ! energy

                          engmpl  = engmpl + t1*d1(ks1,ks2,ks3)

                          ! force

                          t1      = t1*ewld%alpha

                          fx      = fx      - t1*d1(ks11,ks2,ks3)
                          fy      = fy      - t1*d1(ks1,ks21,ks3)
                          fz      = fz      - t1*d1(ks1,ks2,ks31)

                          ! torque on iatm

                          tix     = tix     + impx(ii)*tmpi
                          tiy     = tiy     + impy(ii)*tmpi
                          tiz     = tiz     + impz(ii)*tmpi

                          ! torque on jatm

                          tjx     = tjx     + jmpx(jj)*tmpj
                          tjy     = tjy     + jmpy(jj)*tmpj
                          tjz     = tjz     + jmpz(jj)*tmpj

                          sx = -sx
                        End Do

                        sy = -sy
                      End Do

                      sz = -sz
                    End Do

                  End If

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          End If

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

            mpoles%torque_x(jatm)=mpoles%torque_x(jatm)+tjx
            mpoles%torque_y(jatm)=mpoles%torque_y(jatm)+tjy
            mpoles%torque_z(jatm)=mpoles%torque_z(jatm)+tjz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! accumulate potential energy

            engcpe_rl = engcpe_rl + engmpl

            ! calculate virial

            vircpe_rl = vircpe_rl - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

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

      ! and torques due to multipoles

      mpoles%torque_x(iatm)=mpoles%torque_x(iatm)+scl*tix
      mpoles%torque_y(iatm)=mpoles%torque_y(iatm)+scl*tiy
      mpoles%torque_z(iatm)=mpoles%torque_z(iatm)+scl*tiz

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

  End Subroutine ewald_real_mforces

  Subroutine ewald_real_mforces_d(iatm,xxt,yyt,zzt,rrt,engcpe_rl,vircpe_rl, &
    stress,ewld,neigh,mpoles,electro,config,comm)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using multipoles with the ewald real space
    !! kernel
    !!
    !! copyright - daresbury laboratory
    !! author    - h.a.boateng february 2014
    !! amended   - i.t.todorov february 2016
    !!
    !!----------------------------------------------------------------------!

    Integer,                                          Intent( In    ) :: iatm
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe_rl,vircpe_rl
    Real( Kind = wp ), Dimension( 1:9 ),              Intent( InOut ) :: stress
    Type( ewald_type ),                          Intent( InOut ) :: ewld
    Type( mpole_type ),                               Intent( InOut ) :: mpoles
    Type( electrostatic_type ), Intent( InOut    ) :: electro
    Type( comms_type ),                               Intent( In    ) :: comm
    Type( configuration_type ),                       Intent( InOut ) :: config

    Integer           :: fail,idi,jatm,k,m

    Real( Kind = wp ) :: scl,rsq,rrr,engmpl,fix,fiy,fiz,fx,fy,fz,  &
      strs1,strs2,strs3,strs5,strs6,strs9,      &
      ppp,vk0,vk1,vk2,t1,t2,tx0,ty0,tz0,        &
      siceng,tmpi,tmpj,tix,talp2,               &
      tiy,tiz,tjx,tjy,tjz,b0,b1,b2,b3,b4,b5,    &
      alpsqrpi,exparr,ecc,ecd,ecq,edd,          &
      edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,      &
      zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,      &
      tm9,tm10,tm11,tm12,tm13,tm14,tm15,        &
      dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,     &
      dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,       &
      dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,            &
      dpxzzz,dpxxxx,dpyyyy,dpzzzz,              &
      dpxxxy,dpxxxz,dpxxyy,dpxxzz,dpxxyz,       &
      dpxyyy,dpyzzz,dpyyzz,dpxyyz,dpxyzz,       &
      dpyyyz,dpxxxxx,dpyyyyy,dpzzzzz,dpxxxxy,   &
      dpxxxxz,dpyyyyz,dpxzzzz,dpyzzzz,dpxyyyy,  &
      dpxxxyy,dpxxxzz,dpyyyzz,dpyyzzz,dpxxyyy,  &
      dpxxzzz,dpxyyzz,dpxxyyz,dpxxyzz,dpxyyyz,  &
      dpxyzzz,dpxxxyz,xyz,                      &
      thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,      &
      imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,  &
      imp9,imp10,jmp1,jmp2,jmp3,jmp4,jmp5,      &
      jmp6,jmp7,jmp8,jmp9,jmp10,impx1,impx2,    &
      impx3,impx4,impx5,impx6,impx7,impx8,      &
      impx9,impx10,jmpx1,jmpx2,jmpx3,jmpx4,     &
      jmpx5,jmpx6,jmpx7,jmpx8,jmpx9,jmpx10,     &
      impy1,impy2,impy3,impy4,impy5,impy6,      &
      impy7,impy8,impy9,impy10,jmpy1,jmpy2,     &
      jmpy3,jmpy4,jmpy5,jmpy6,jmpy7,jmpy8,      &
      jmpy9,jmpy10,impz1,impz2,impz3,impz4,     &
      impz5,impz6,impz7,impz8,impz9,impz10,     &
      jmpz1,jmpz2,jmpz3,jmpz4,jmpz5,jmpz6,      &
      jmpz7,jmpz8,jmpz9,jmpz10,                 &
      tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10


    Real( Kind = wp ) :: imp(1:mpoles%max_mpoles),jmp(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: jmpx(1:mpoles%max_mpoles),jmpy(1:mpoles%max_mpoles),jmpz(1:mpoles%max_mpoles)
    Logical, save :: newjob
    Character( Len = 256 ) :: message

    If (newjob) Then
      newjob = .false.

      ! generate error function complement tables for ewald sum

      Call erfcgen(neigh%cutoff,ewld%alpha,erfc,erfc_deriv)

      ! coefficients for exponential in recurrence relation

      talp2 = 2.0_wp*ewld%alpha*ewld%alpha
      alpsqrpi = 1.0_wp/(ewld%alpha*sqrpi)

      electro%co1 = talp2*alpsqrpi
      electro%co2 = talp2*electro%co1
      electro%co3 = talp2*electro%co2
      electro%co4 = talp2*electro%co3
      electro%co5 = talp2*electro%co4

      electro%alp2 = ewld%alpha*ewld%alpha

      electro%exclcoef = r4pie0*ewld%alpha /sqrpi/electro%eps

      electro%twzz=-2.0_wp*ewld%alpha**3 *r4pie0/(3.0_wp*sqrpi*electro%eps)
      electro%twtwz=4.0_wp*ewld%alpha**5 *r4pie0/(5.0_wp*sqrpi*electro%eps)
      electro%fozz=12.0_wp*ewld%alpha**5 *r4pie0/(5.0_wp*sqrpi*electro%eps)
    End If

    ! initialise potential energy and virial

    engcpe_rl=0.0_wp
    vircpe_rl=0.0_wp
    siceng   =0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! get the multipoles for site i

    imp=mpoles%global_frame(:,iatm)

    tx0 = 0.0_wp ; ty0 = 0.0_wp ; tz0 = 0.0_wp

    ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

      ! get the components for site i infinitesimal rotations

      impx=mpoles%rotation_x(:,iatm)
      impy=mpoles%rotation_y(:,iatm)
      impz=mpoles%rotation_z(:,iatm)

      imp1=imp(1); impx1=impx(1); impy1=impy(1); impz1=impz(1)

      ! self-interaction

      siceng = siceng+electro%exclcoef*imp1*imp1

      If (mpoles%max_order >= 1) Then

        imp2=imp(2); imp3=imp(3); imp4=imp(4)

        impx2=impx(2); impx3=impx(3); impx4=impx(4)
        impy2=impy(2); impy3=impy(3); impy4=impy(4)
        impz2=impz(2); impz3=impz(3); impz4=impz(4)

        siceng = siceng - electro%twzz*(imp2*imp2+imp3*imp3+imp4*imp4)
        !        tx0    = tx0    - electro%twzz*(impx2*imp2+impx3*imp3+impx4*imp4)
        !        ty0    = ty0    - electro%twzz*(impy2*imp2+impy3*imp3+impy4*imp4)
        !        tz0    = tz0    - electro%twzz*(impz2*imp2+impz3*imp3+impz4*imp4)

      End If

      If (mpoles%max_order == 2) Then

        imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

        impx5=impx(5); impx6=impx(6); impx7=impx(7); impx8=impx(8); impx9=impx(9); impx10=impx(10)
        impy5=impy(5); impy6=impy(6); impy7=impy(7); impy8=impy(8); impy9=impy(9); impy10=impy(10)
        impz5=impz(5); impz6=impz(6); impz7=impz(7); impz8=impz(8); impz9=impz(9); impz10=impz(10)

        siceng = siceng + 2.0_wp*electro%twzz*(imp1*(imp5+imp8+imp10)) + electro%twtwz*(2.0_wp*imp5*(imp8+imp10)+ &
          2.0_wp*imp8*imp10+imp6*imp6+imp7*imp7+imp9*imp9) + electro%fozz*(imp5*imp5+imp8*imp8+    &
          imp10*imp10)
        !        tx0    = tx0 + electro%twzz*imp1*(impx5+impx8+impx10)+electro%twtwz*(impx5*(imp8+imp10)+impx6*imp6 +     &
        !                 impx7*imp7+impx8*(imp5+imp10)+impx9*imp9+impx10*(imp5+imp8))+electro%fozz*(impx5*imp5 + &
        !                 impx8*imp8+impx10*imp10)
        !        ty0    = ty0 + electro%twzz*imp1*(impy5+impy8+impy10)+electro%twtwz*(impy5*(imp8+imp10)+impy6*imp6 +     &
        !                 impy7*imp7+impy8*(imp5+imp10)+impy9*imp9+impy10*(imp5+imp8))+electro%fozz*(impy5*imp5 + &
        !                 impy8*imp8+impy10*imp10)
        !        tz0    = tz0 + electro%twzz*imp1*(impz5+impz8+impz10)+electro%twtwz*(impz5*(imp8+imp10)+impz6*imp6 +     &
        !                 impz7*imp7+impz8*(imp5+imp10)+impz9*imp9+impz10*(imp5+imp8))+electro%fozz*(impz5*imp5 + &
        !                 impz8*imp8+impz10*imp10)

      End If

      ! multipole scaler

      scl=r4pie0/electro%eps

      ! rescale multipoles

      imp=imp*scl

      imp1=imp(1)

      If (mpoles%max_order >= 1) Then

        imp2=imp(2); imp3=imp(3); imp4=imp(4)

      End If

      If (mpoles%max_order == 2) Then

        imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

      End If

      ! add self-interaction energy

      engcpe_rl = engcpe_rl - siceng
      ewld%engsic    = ewld%engsic    - siceng

      ! initialize torques for atom i (temporary)

      tix = .0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

      !     tix = tix-tx0
      !     tiy = tiy-ty0
      !     tiz = tiz-tz0

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! start of primary loop for forces evaluation

      Do m=1,neigh%list(0,iatm)

        ! atomic index

        jatm=neigh%list(m,iatm)

        ! get the multipoles for site j and the components for its infinitesimal rotations

        jmp=mpoles%global_frame(:,jatm)

        jmpx=mpoles%rotation_x(:,jatm)
        jmpy=mpoles%rotation_y(:,jatm)
        jmpz=mpoles%rotation_z(:,jatm)

        jmp1=jmp(1); jmpx1=jmpx(1); jmpy1=jmpy(1); jmpz1=jmpz(1)

        If (mpoles%max_order >= 1) Then

          jmp2=jmp(2); jmp3=jmp(3); jmp4=jmp(4)

          jmpx2=jmpx(2); jmpx3=jmpx(3); jmpx4=jmpx(4)
          jmpy2=jmpy(2); jmpy3=jmpy(3); jmpy4=jmpy(4)
          jmpz2=jmpz(2); jmpz3=jmpz(3); jmpz4=jmpz(4)

        End If

        If (mpoles%max_order == 2) Then

          jmp5=jmp(5); jmp6=jmp(6); jmp7=jmp(7); jmp8=jmp(8); jmp9=jmp(9); jmp10=jmp(10)

          jmpx5=jmpx(5); jmpx6=jmpx(6); jmpx7=jmpx(7); jmpx8=jmpx(8); jmpx9=jmpx(9); jmpx10=jmpx(10)
          jmpy5=jmpy(5); jmpy6=jmpy(6); jmpy7=jmpy(7); jmpy8=jmpy(8); jmpy9=jmpy(9); jmpy10=jmpy(10)
          jmpz5=jmpz(5); jmpz6=jmpz(6); jmpz7=jmpz(7); jmpz8=jmpz(8); jmpz9=jmpz(9); jmpz10=jmpz(10)

        End If

        ! interatomic distance

        rrr = rrt(m)

        ! truncation of potential

        If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

          ! Squared distance

          rsq=rrr**2

          engmpl = 0.0_wp
          ecc = 0.0_wp; ecd = 0.0_wp; ecq = 0.0_wp
          edd = 0.0_wp; edq = 0.0_wp; eqq = 0.0_wp
          fx = 0.0_wp; fy = 0.0_wp; fz = 0.0_wp
          tjx = 0.0_wp; tjy = 0.0_wp; tjz = 0.0_wp

          xx=xxt(m); yy=yyt(m); zz=zzt(m)
          xx2=xx*xx; yy2=yy*yy; zz2=zz*zz

          ! evaluate the exponential

          exparr = Exp(-electro%alp2*rsq)

          ! get the value of the kernel using 3pt interpolation

          ! compute recurrence terms

          ! b0 = (t1 + (t2-t1)*ppp*0.5_wp)
          b0 = three_p_interp(erfc, rrr)
          b1 = (b0        + electro%co1*exparr)/rsq
          b2 = (3.0_wp*b1 + electro%co2*exparr)/rsq
          b3 = (5.0_wp*b2 + electro%co3*exparr)/rsq
          b4 = (7.0_wp*b3 + electro%co4*exparr)/rsq
          b5 = (9.0_wp*b4 + electro%co5*exparr)/rsq

          ! charge-charge interaction

          ijmp= imp(1)*jmp(1)

          ecc = ijmp*b0

          bijmp = ijmp*b1
          fx  = bijmp*xx
          fy  = bijmp*yy
          fz  = bijmp*zz

          ! There is no torque for charges

          If (mpoles%max_order == 1) Then

            ! charge-dipole interactions

            tm1 = imp2*jmp1-imp1*jmp2
            tm2 = imp3*jmp1-imp1*jmp3
            tm3 = imp4*jmp1-imp1*jmp4

            tmpi= jmp1*b1
            tmpj= imp1*b1

            dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
            dpxx= xx2*b2-b1
            dpyy= yy2*b2-b1
            dpzz= zz2*b2-b1
            dpxy= xx*yy*b2
            dpxz= xx*zz*b2
            dpyz= yy*zz*b2

            ecd = tm1*dpx + tm2*dpy + tm3*dpz

            fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
            fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
            fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

            ! dipole-dipole interactions

            dpxxx = xx*(3.0_wp*b2-xx2*b3)
            dpyyy = yy*(3.0_wp*b2-yy2*b3)
            dpzzz = zz*(3.0_wp*b2-zz2*b3)
            dpxxy = yy*(b2-xx2*b3)
            dpxxz = zz*(b2-xx2*b3)
            dpxyy = xx*(b2-yy2*b3)
            dpyyz = zz*(b2-yy2*b3)
            dpxzz = xx*(b2-zz2*b3)
            dpyzz = yy*(b2-zz2*b3)
            dpxyz = -xx*yy*zz*b3

            tm1 = imp2*jmp2
            tm2 = imp2*jmp3+imp3*jmp2
            tm3 = imp2*jmp4+imp4*jmp2
            tm4 = imp3*jmp3
            tm5 = imp3*jmp4+imp4*jmp3
            tm6 = imp4*jmp4

            edd = -(tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz)

            fx = fx + (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
            fy = fy + (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
            fz = fz + (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

            ! Torques

            tmp2 = jmp1*dpx-jmp2*dpxx-jmp3*dpxy-jmp4*dpxz
            tmp3 = jmp1*dpy-jmp2*dpxy-jmp3*dpyy-jmp4*dpyz
            tmp4 = jmp1*dpz-jmp2*dpxz-jmp3*dpyz-jmp4*dpzz

            tix  = tix + impx2*tmp2 + impx3*tmp3 + impx4*tmp4
            tiy  = tiy + impy2*tmp2 + impy3*tmp3 + impy4*tmp4
            tiz  = tiz + impz2*tmp2 + impz3*tmp3 + impz4*tmp4

            tmp2 = imp1*dpx-imp2*dpxx-imp3*dpxy-imp4*dpxz
            tmp3 = imp1*dpy-imp2*dpxy-imp3*dpyy-imp4*dpyz
            tmp4 = imp1*dpz-imp2*dpxz-imp3*dpyz-imp4*dpzz

            tjx  = tjx + jmpx2*tmp2 + jmpx3*tmp3 + jmpx4*tmp4
            tjy  = tjy + jmpy2*tmp2 + jmpy3*tmp3 + jmpy4*tmp4
            tjz  = tjz + jmpz2*tmp2 + jmpz3*tmp3 + jmpz4*tmp4

          End If

          If (mpoles%max_order == 2) Then

            ! charge-dipole interactions

            tm1 = imp2*jmp1-imp1*jmp2
            tm2 = imp3*jmp1-imp1*jmp3
            tm3 = imp4*jmp1-imp1*jmp4

            tmpi= jmp1*b1
            tmpj= imp1*b1

            dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
            dpxx= xx2*b2-b1
            dpyy= yy2*b2-b1
            dpzz= zz2*b2-b1
            dpxy= xx*yy*b2
            dpxz= xx*zz*b2
            dpyz= yy*zz*b2

            ecd = tm1*dpx + tm2*dpy + tm3*dpz

            fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
            fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
            fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

            ! dipole-dipole interactions

            dpxxx = xx*(3.0_wp*b2-xx2*b3)
            dpyyy = yy*(3.0_wp*b2-yy2*b3)
            dpzzz = zz*(3.0_wp*b2-zz2*b3)
            dpxxy = yy*(b2-xx2*b3)
            dpxxz = zz*(b2-xx2*b3)
            dpxyy = xx*(b2-yy2*b3)
            dpyyz = zz*(b2-yy2*b3)
            dpxzz = xx*(b2-zz2*b3)
            dpyzz = yy*(b2-zz2*b3)
            dpxyz = -xx*yy*zz*b3

            tm1 = imp2*jmp2
            tm2 = imp2*jmp3+imp3*jmp2
            tm3 = imp2*jmp4+imp4*jmp2
            tm4 = imp3*jmp3
            tm5 = imp3*jmp4+imp4*jmp3
            tm6 = imp4*jmp4

            edd = -(tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz)

            ! The forces for dipole-dipole interactions will be computed as part of the charge-quadrupole interactions

            ! charge-quadrupole interactions

            tm1 = imp1*jmp5+imp5*jmp1-tm1; tm2=imp1*jmp6+imp6*jmp1-tm2
            tm3 = imp1*jmp7+imp7*jmp1-tm3; tm4=imp1*jmp8+imp8*jmp1-tm4
            tm5 = imp1*jmp9+imp9*jmp1-tm5; tm6=imp1*jmp10+imp10*jmp1-tm6

            ecq = tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz-edd

            fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
            fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
            fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

            ! dipole-quadrupole interactions

            thrb2  = 3.0_wp*b2; thrb3 = 3.0_wp*b3

            dpxxxx = thrb2 - xx2*(6.0_wp*b3 - xx2*b4)
            dpyyyy = thrb2 - yy2*(6.0_wp*b3 - yy2*b4)
            dpzzzz = thrb2 - zz2*(6.0_wp*b3 - zz2*b4)
            dpxxxy = yy*xx*(xx2*b4-thrb3)
            dpxxxz = zz*xx*(xx2*b4-thrb3)
            dpxyyy = xx*yy*(yy2*b4-thrb3)
            dpyyyz = yy*zz*(yy2*b4-thrb3)
            dpyzzz = yy*zz*(zz2*b4-thrb3)
            dpxzzz = xx*zz*(zz2*b4-thrb3)
            dpxxyy = b2-(xx2+yy2)*b3+xx2*yy2*b4
            dpxxzz = b2-(xx2+zz2)*b3+xx2*zz2*b4
            dpyyzz = b2-(yy2+zz2)*b3+yy2*zz2*b4
            dpxxyz = yy*zz*(xx2*b4-b3)
            dpxyyz = xx*zz*(yy2*b4-b3)
            dpxyzz = xx*yy*(zz2*b4-b3)

            tm1 = imp2*jmp5-imp5*jmp2
            tm2 = imp2*jmp6-imp6*jmp2+imp3*jmp5-imp5*jmp3
            tm3 = imp2*jmp7-imp7*jmp2+imp4*jmp5-imp5*jmp4
            tm4 = imp2*jmp8-imp8*jmp2+imp3*jmp6-imp6*jmp3
            tm5 = imp2*jmp9-imp9*jmp2+imp3*jmp7-imp7*jmp3+imp4*jmp6-imp6*jmp4
            tm6 = imp2*jmp10-imp10*jmp2+imp4*jmp7-imp7*jmp4
            tm7 = imp3*jmp8-imp8*jmp3
            tm8 = imp3*jmp9-imp9*jmp3+imp4*jmp8-imp8*jmp4
            tm9 = imp3*jmp10-imp10*jmp3+imp4*jmp9-imp9*jmp4
            tm10= imp4*jmp10-imp10*jmp4

            edq = tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz + &
              tm6*dpxzz+tm7*dpyyy+tm8*dpyyz+tm9*dpyzz+tm10*dpzzz

            fx = fx - (tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz + &
              tm6*dpxxzz+tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz)

            fy = fy - (tm1*dpxxxy+tm2*dpxxyy+tm3*dpxxyz+tm4*dpxyyy+tm5*dpxyyz + &
              tm6*dpxyzz+tm7*dpyyyy+tm8*dpyyyz+tm9*dpyyzz+tm10*dpyzzz)

            fz = fz - (tm1*dpxxxz+tm2*dpxxyz+tm3*dpxxzz+tm4*dpxyyz+tm5*dpxyzz + &
              tm6*dpxzzz+tm7*dpyyyz+tm8*dpyyzz+tm9*dpyzzz+tm10*dpzzzz)

            ! quadrupole-quadrupole interactions

            fifb3 = 15.0_wp*b3; thrb4 = 3.0_wp*b4; sixb4 = 6.0_wp*b4; tenb4 = 10.0_wp*b4
            xyz   = xx*yy*zz

            dpxxxxx = xx*(xx2*(tenb4-xx2*b5)-fifb3)
            dpyyyyy = yy*(yy2*(tenb4-yy2*b5)-fifb3)
            dpzzzzz = zz*(zz2*(tenb4-zz2*b5)-fifb3)
            dpxxxxy = yy*(xx2*(sixb4-xx2*b5)-thrb3)
            dpxxxxz = zz*(xx2*(sixb4-xx2*b5)-thrb3)
            dpyyyyz = zz*(yy2*(sixb4-yy2*b5)-thrb3)
            dpyzzzz = yy*(zz2*(sixb4-zz2*b5)-thrb3)
            dpxzzzz = xx*(zz2*(sixb4-zz2*b5)-thrb3)
            dpxyyyy = xx*(yy2*(sixb4-yy2*b5)-thrb3)
            dpxxxyy = xx*((xx2+3.0_wp*yy2)*b4-xx2*yy2*b5-thrb3)
            dpxxxzz = xx*((xx2+3.0_wp*zz2)*b4-xx2*zz2*b5-thrb3)
            dpyyyzz = yy*((yy2+3.0_wp*zz2)*b4-yy2*zz2*b5-thrb3)
            dpyyzzz = zz*((3.0_wp*yy2+zz2)*b4-yy2*zz2*b5-thrb3)
            dpxxyyy = yy*((3.0_wp*xx2+yy2)*b4-xx2*yy2*b5-thrb3)
            dpxxzzz = zz*((3.0_wp*xx2+zz2)*b4-xx2*zz2*b5-thrb3)
            dpxyyzz = xx*((yy2+zz2)*b4-yy2*zz2*b5-b3)
            dpxxyyz = zz*((xx2+yy2)*b4-xx2*yy2*b5-b3)
            dpxxyzz = yy*((xx2+zz2)*b4-xx2*zz2*b5-b3)
            dpxyyyz = xyz*(thrb4-yy2*b5)
            dpxyzzz = xyz*(thrb4-zz2*b5)
            dpxxxyz = xyz*(thrb4-xx2*b5)

            tm1 = imp5*jmp5
            tm2 = imp5*jmp6+imp6*jmp5
            tm3 = imp5*jmp7+imp7*jmp5
            tm4 = imp5*jmp8+imp8*jmp5+imp6*jmp6
            tm5 = imp5*jmp9+imp9*jmp5+imp6*jmp7+imp7*jmp6
            tm6 = imp5*jmp10+imp10*jmp5+imp7*jmp7
            tm7 = imp6*jmp8+imp8*jmp6
            tm8 = imp6*jmp9+imp9*jmp6+imp7*jmp8+imp8*jmp7
            tm9 = imp6*jmp10+imp10*jmp6+imp7*jmp9+imp9*jmp7
            tm10= imp7*jmp10+imp10*jmp7
            tm11= imp8*jmp8
            tm12= imp8*jmp9+imp9*jmp8
            tm13= imp8*jmp10+imp10*jmp8+imp9*jmp9
            tm14= imp9*jmp10+imp10*jmp9
            tm15= imp10*jmp10

            eqq = tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz+tm6*dpxxzz +    &
              tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz+tm11*dpyyyy+tm12*dpyyyz + &
              tm13*dpyyzz+tm14*dpyzzz+tm15*dpzzzz

            fx = fx - (tm1*dpxxxxx+tm2*dpxxxxy+tm3*dpxxxxz+tm4*dpxxxyy+tm5*dpxxxyz+tm6*dpxxxzz +    &
              tm7*dpxxyyy+tm8*dpxxyyz+tm9*dpxxyzz+tm10*dpxxzzz+tm11*dpxyyyy+tm12*dpxyyyz + &
              tm13*dpxyyzz+tm14*dpxyzzz+tm15*dpxzzzz)

            fy = fy - (tm1*dpxxxxy+tm2*dpxxxyy+tm3*dpxxxyz+tm4*dpxxyyy+tm5*dpxxyyz+tm6*dpxxyzz +    &
              tm7*dpxyyyy+tm8*dpxyyyz+tm9*dpxyyzz+tm10*dpxyzzz+tm11*dpyyyyy+tm12*dpyyyyz + &
              tm13*dpyyyzz+tm14*dpyyzzz+tm15*dpyzzzz)

            fz = fz - (tm1*dpxxxxz+tm2*dpxxxyz+tm3*dpxxxzz+tm4*dpxxyyz+tm5*dpxxyzz+tm6*dpxxzzz +    &
              tm7*dpxyyyz+tm8*dpxyyzz+tm9*dpxyzzz+tm10*dpxzzzz+tm11*dpyyyyz+tm12*dpyyyzz + &
              tm13*dpyyzzz+tm14*dpyzzzz+tm15*dpzzzzz)

            ! second-order torques

            tmp2 = jmp1*dpx-jmp2*dpxx-jmp3*dpxy-jmp4*dpxz+jmp5*dpxxx+jmp6*dpxxy+jmp7*dpxxz+jmp8*dpxyy + &
              jmp9*dpxyz+jmp10*dpxzz
            tmp3 = jmp1*dpy-jmp2*dpxy-jmp3*dpyy-jmp4*dpyz+jmp5*dpxxy+jmp6*dpxyy+jmp7*dpxyz+jmp8*dpyyy + &
              jmp9*dpyyz+jmp10*dpyzz
            tmp4 = jmp1*dpz-jmp2*dpxz-jmp3*dpyz-jmp4*dpzz+jmp5*dpxxz+jmp6*dpxyz+jmp7*dpxzz+jmp8*dpyyz + &
              jmp9*dpyzz+jmp10*dpzzz

            tmp5 = jmp1*dpxx-jmp2*dpxxx-jmp3*dpxxy-jmp4*dpxxz+jmp5*dpxxxx+jmp6*dpxxxy+jmp7*dpxxxz + &
              jmp8*dpxxyy+jmp9*dpxxyz+jmp10*dpxxzz
            tmp6 = jmp1*dpxy-jmp2*dpxxy-jmp3*dpxyy-jmp4*dpxyz+jmp5*dpxxxy+jmp6*dpxxyy+jmp7*dpxxyz + &
              jmp8*dpxyyy+jmp9*dpxyyz+jmp10*dpxyzz
            tmp7 = jmp1*dpxz-jmp2*dpxxz-jmp3*dpxyz-jmp4*dpxzz+jmp5*dpxxxz+jmp6*dpxxyz+jmp7*dpxxzz + &
              jmp8*dpxyyz+jmp9*dpxyzz+jmp10*dpxzzz
            tmp8 = jmp1*dpyy-jmp2*dpxyy-jmp3*dpyyy-jmp4*dpyyz+jmp5*dpxxyy+jmp6*dpxyyy+jmp7*dpxyyz + &
              jmp8*dpyyyy+jmp9*dpyyyz+jmp10*dpyyzz
            tmp9 = jmp1*dpyz-jmp2*dpxyz-jmp3*dpyyz-jmp4*dpyzz+jmp5*dpxxyz+jmp6*dpxyyz+jmp7*dpxyzz + &
              jmp8*dpyyyz+jmp9*dpyyzz+jmp10*dpyzzz
            tmp10= jmp1*dpzz-jmp2*dpxzz-jmp3*dpyzz-jmp4*dpzzz+jmp5*dpxxzz+jmp6*dpxyzz+jmp7*dpxzzz + &
              jmp8*dpyyzz+jmp9*dpyzzz+jmp10*dpzzzz

            tix = tix + impx2*tmp2+impx3*tmp3+impx4*tmp4+impx5*tmp5+impx6*tmp6+impx7*tmp7+impx8*tmp8 + &
              impx9*tmp9+impx10*tmp10
            tiy = tiy + impy2*tmp2+impy3*tmp3+impy4*tmp4+impy5*tmp5+impy6*tmp6+impy7*tmp7+impy8*tmp8 + &
              impy9*tmp9+impy10*tmp10
            tiz = tiz + impz2*tmp2+impz3*tmp3+impz4*tmp4+impz5*tmp5+impz6*tmp6+impz7*tmp7+impz8*tmp8 + &
              impz9*tmp9+impz10*tmp10

            tmp2 = imp1*dpx-imp2*dpxx-imp3*dpxy-imp4*dpxz+imp5*dpxxx+imp6*dpxxy+imp7*dpxxz+imp8*dpxyy + &
              imp9*dpxyz+imp10*dpxzz
            tmp3 = imp1*dpy-imp2*dpxy-imp3*dpyy-imp4*dpyz+imp5*dpxxy+imp6*dpxyy+imp7*dpxyz+imp8*dpyyy + &
              imp9*dpyyz+imp10*dpyzz
            tmp4 = imp1*dpz-imp2*dpxz-imp3*dpyz-imp4*dpzz+imp5*dpxxz+imp6*dpxyz+imp7*dpxzz+imp8*dpyyz + &
              imp9*dpyzz+imp10*dpzzz

            tmp5 = imp1*dpxx-imp2*dpxxx-imp3*dpxxy-imp4*dpxxz+imp5*dpxxxx+imp6*dpxxxy+imp7*dpxxxz + &
              imp8*dpxxyy+imp9*dpxxyz+imp10*dpxxzz
            tmp6 = imp1*dpxy-imp2*dpxxy-imp3*dpxyy-imp4*dpxyz+imp5*dpxxxy+imp6*dpxxyy+imp7*dpxxyz + &
              imp8*dpxyyy+imp9*dpxyyz+imp10*dpxyzz
            tmp7 = imp1*dpxz-imp2*dpxxz-imp3*dpxyz-imp4*dpxzz+imp5*dpxxxz+imp6*dpxxyz+imp7*dpxxzz + &
              imp8*dpxyyz+imp9*dpxyzz+imp10*dpxzzz
            tmp8 = imp1*dpyy-imp2*dpxyy-imp3*dpyyy-imp4*dpyyz+imp5*dpxxyy+imp6*dpxyyy+imp7*dpxyyz + &
              imp8*dpyyyy+imp9*dpyyyz+imp10*dpyyzz
            tmp9 = imp1*dpyz-imp2*dpxyz-imp3*dpyyz-imp4*dpyzz+imp5*dpxxyz+imp6*dpxyyz+imp7*dpxyzz + &
              imp8*dpyyyz+imp9*dpyyzz+imp10*dpyzzz
            tmp10= imp1*dpzz-imp2*dpxzz-imp3*dpyzz-imp4*dpzzz+imp5*dpxxzz+imp6*dpxyzz+imp7*dpxzzz + &
              imp8*dpyyzz+imp9*dpyzzz+imp10*dpzzzz

            tjx = tjx + jmpx2*tmp2+jmpx3*tmp3+jmpx4*tmp4+jmpx5*tmp5+jmpx6*tmp6+jmpx7*tmp7+jmpx8*tmp8 + &
              jmpx9*tmp9+jmpx10*tmp10
            tjy = tjy + jmpy2*tmp2+jmpy3*tmp3+jmpy4*tmp4+jmpy5*tmp5+jmpy6*tmp6+jmpy7*tmp7+jmpy8*tmp8 + &
              jmpy9*tmp9+jmpy10*tmp10
            tjz = tjz + jmpz2*tmp2+jmpz3*tmp3+jmpz4*tmp4+jmpz5*tmp5+jmpz6*tmp6+jmpz7*tmp7+jmpz8*tmp8 + &
              jmpz9*tmp9+jmpz10*tmp10

          End If

          engmpl = ecc+ecd+ecq+edd+edq+eqq

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

            mpoles%torque_x(jatm)=mpoles%torque_x(jatm)+tjx
            mpoles%torque_y(jatm)=mpoles%torque_y(jatm)+tjy
            mpoles%torque_z(jatm)=mpoles%torque_z(jatm)+tjz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! accumulate potential energy

            engcpe_rl = engcpe_rl + engmpl

            ! calculate virial

            vircpe_rl = vircpe_rl - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

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

      ! and torques due to multipoles

      mpoles%torque_x(iatm)=mpoles%torque_x(iatm)+scl*tix
      mpoles%torque_y(iatm)=mpoles%torque_y(iatm)+scl*tiy
      mpoles%torque_z(iatm)=mpoles%torque_z(iatm)+scl*tiz

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

  End Subroutine ewald_real_mforces_d

  Subroutine ewald_spme_mforces(engcpe_rc,vircpe_rc,stress,ewld,mpoles,electro, &
    domain,config,comm)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! due to multipolar interactions in a periodic system using the smooth
    !! particle mesh ewald method for multipoles
    !!
    !! This version allows for extension to arbitrary order
    !!
    !! Note: (fourier) reciprocal space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - h.a.boateng february 2016
    !! amended   - i.t.todorov march 2016
    !!
    !!----------------------------------------------------------------------!

    Real( Kind = wp ), Intent(   Out ) :: engcpe_rc,vircpe_rc
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Type( ewald_type ), Intent( InOut ) :: ewld
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm
    Type( configuration_type ),                       Intent( InOut ) :: config


    Logical              :: llspl
    Integer              :: fail(1:4), i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt, &
      jjtjjb,counter, k1,k2,k3,s1,s2,s3,ks1,ks2,ks3, mm,nn

    Real( Kind = wp )    :: det,rcell(1:9),celprp(1:10),ralph,rvolm,scale,exclcoef, &
      rcpcut,rcpct2,strs(1:9),eng,akv,tmp,bb1,bb2,bb3,mptmp,         &
      rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3,           &
      rrkxy,rrkxz,rrkyx,rrkyz,rrkzx,rrkzy, sq1,sq2,sq3, tx,ty,tz,    &
      tmpi,tix,tiy,tiz,timp,dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3

    Real( Kind = wp )    :: imp(1:mpoles%max_mpoles)
    Real( Kind = wp )    :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)

    Complex( Kind = wp ) :: vterm,pterm

    ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter    :: &
      uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


    ! B-spline coefficients

    Complex( Kind = wp ), Dimension( : ),       Allocatable       :: ww1,ww2,ww3

    Real( Kind = wp ),    Dimension( : ),       Allocatable       :: csp
    Real( Kind = wp ),    Dimension( : ),       Allocatable       :: txx,tyy,tzz
    Integer,              Dimension( : ),       Allocatable       :: ixx,iyy,izz,it
    Real( Kind = wp ),    Dimension( :,:,: ),   Allocatable       :: bsddx,bsddy,bsddz
    Real( Kind = wp ),    Dimension( :,: ),     Allocatable       :: bspx,bspy,bspz
    Real( Kind = wp ),    Dimension( : ),       Allocatable       :: bdx,bdy,bdz


    ! temporary qqc

    Real( Kind = wp )    :: qqc_tmp

    ! DaFT arrays local indices

    Integer              :: j_local, k_local, l_local
    Logical,  save :: newjob
    Character( Len = 256 ) :: message

    llspl=.true.


    fail=0
    If (newjob) Then
      newjob = .false.

!!! BEGIN DD SPME VARIABLES
      ! 3D charge array construction (bottom and top) indices

      electro%ixb_mf=domain%idx*(ewld%kspace%k_vec_dim(1)/domain%nx)+1
      electro%ixt_mf=(domain%idx+1)*(ewld%kspace%k_vec_dim(1)/domain%nx)
      electro%iyb_mf=domain%idy*(ewld%kspace%k_vec_dim(2)/domain%ny)+1
      electro%iyt_mf=(domain%idy+1)*(ewld%kspace%k_vec_dim(2)/domain%ny)
      electro%izb_mf=domain%idz*(ewld%kspace%k_vec_dim(3)/domain%nz)+1
      electro%izt_mf=(domain%idz+1)*(ewld%kspace%k_vec_dim(3)/domain%nz)

      electro%ixbm1_r=Real(electro%ixb_mf-1,wp)
      electro%ixtm0_r_mf=Nearest( Real(electro%ixt_mf,wp) , -1.0_wp )
      electro%iybm1_r_mf=Real(electro%iyb_mf-1,wp)
      electro%iytm0_r_mf=Nearest( Real(electro%iyt_mf,wp) , -1.0_wp )
      electro%izbm1_r_mf=Real(electro%izb_mf-1,wp)
      electro%iztm0_r_mf=Nearest( Real(electro%izt_mf,wp) , -1.0_wp )

      ! Real values of kmax vectors

      electro%kmaxa_r_mf=Real(ewld%kspace%k_vec_dim(1),wp)
      electro%kmaxb_r_mf=Real(ewld%kspace%k_vec_dim(2),wp)
      electro%kmaxc_r_mf=Real(ewld%kspace%k_vec_dim(3),wp)

      !!! END DD SPME VARIABLES

      !!! BEGIN CARDINAL B-SPLINES SET-UP
      ! allocate the complex exponential arrays

      call bspline_coeffs_gen(ewld%kspace,ewld%bspline)

!!! END CARDINAL B-SPLINES SET-UP

      !!! BEGIN DAFT SET-UP
      ! domain local block limits of kmax space

      electro%block_x_mf = ewld%kspace%k_vec_dim(1) / domain%nx
      electro%block_y_mf = ewld%kspace%k_vec_dim(2) / domain%ny
      electro%block_z_mf = ewld%kspace%k_vec_dim(3) / domain%nz

      ! set up the parallel fft and useful related quantities

      Call initialize_fft( 3, (/ ewld%kspace%k_vec_dim(1), ewld%kspace%k_vec_dim(2), ewld%kspace%k_vec_dim(3) /), &
        (/ domain%nx, domain%ny, domain%nz /), (/ domain%idx, domain%idy, domain%idz /),   &
        (/ electro%block_x_mf, electro%block_y_mf, electro%block_z_mf /),               &
        comm%comm, electro%context_mf )

      ! set up the indexing arrays for each dimension (NOT deallocated manually)

      Allocate ( electro%index_x_mf( 1:electro%block_x_mf ), Stat = fail(1) )
      Allocate ( electro%index_y_mf( 1:electro%block_y_mf ), Stat = fail(2) )
      Allocate ( electro%index_z_mf( 1:electro%block_z_mf ), Stat = fail(3) )
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'SPME index arrays allocation failure'
        Call error(0,message)
      End If

      Call pfft_indices( ewld%kspace%k_vec_dim(1), electro%block_x_mf, domain%idx, domain%nx, electro%index_x_mf )
      Call pfft_indices( ewld%kspace%k_vec_dim(2), electro%block_y_mf, domain%idy, domain%ny, electro%index_y_mf )
      Call pfft_indices( ewld%kspace%k_vec_dim(3), electro%block_z_mf, domain%idz, domain%nz, electro%index_z_mf )

      ! workspace arrays for DaFT

      Allocate ( electro%qqc_local_mf( 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ), Stat = fail(1) )
      Allocate ( electro%qqq_local_mf( 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ), Stat = fail(2) )
      Allocate ( electro%qtc_local_mf( 1:3, 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ), &
        electro%qt1_local_mf( 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ),      &
        electro%qt2_local_mf( 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ),      &
        electro%qt3_local_mf( 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ), Stat = fail(3) )
      Allocate ( electro%pfft_work_mf( 1:electro%block_x_mf, 1:electro%block_y_mf, 1:electro%block_z_mf ), Stat = fail(4) )
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'SPME DaFT workspace arrays allocation failure'
        Call error(0,message)
      End If

      !!! END DAFT SET-UP

      ! compute derivatives of kernel

      Call limit_erfr_deriv(8,ewld%alpha,electro%d1_mf)
    End If

    Allocate (txx(1:config%mxatms),tyy(1:config%mxatms),tzz(1:config%mxatms),                            Stat = fail(1))
    Allocate (ixx(1:config%mxatms),iyy(1:config%mxatms),izz(1:config%mxatms),it(1:config%mxatms),               Stat = fail(2))
    Allocate ( &
      & bdx(0:ewld%bspline%num_splines),bsddx(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
      & bdy(0:ewld%bspline%num_splines),bsddy(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
      & bdz(0:ewld%bspline%num_splines),bsddz(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
      & Stat = fail(3))
    Allocate ( &
      & bspx(1:ewld%bspline%num_splines,1:config%mxatms), &
      & bspy(1:ewld%bspline%num_splines,1:config%mxatms), &
      & bspz(1:ewld%bspline%num_splines,1:config%mxatms), Stat = fail(4))
    ! Allocate ( ewld%bspline%derivs(3, 0:2, ewld%bspline%num_splines, config%mxatms, stat = fail(4))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_mforces allocation failure'
      Call error(0,message)
    End If

    ! compute self-interaction energy (per node) and torques

    ewld%engsic=0.0_wp; exclcoef = 0.5_wp*r4pie0/electro%eps
    Do i=1,config%natms

      ! get the multipoles for site i

      imp=mpoles%global_frame(:,i)

      ! ignore interaction if the charge is zero

      If (Maxval(Abs(imp)) > zero_plus) Then

        ! get the components for site i infinitesimal rotations

        impx=mpoles%rotation_x(:,i)
        impy=mpoles%rotation_y(:,i)
        impz=mpoles%rotation_z(:,i)

        ! initialize torques for atom i (temporary)

        tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

        tz=1.0_wp
        Do k3=0,mpoles%max_order

          ty=tz
          Do k2=0,mpoles%max_order-k3

            tx=ty
            Do k1=0,mpoles%max_order-k3-k2
              nn = mpoles%map(k1,k2,k3)

              If (Abs(imp(nn)) > zero_plus) Then
                timp=imp(nn)

                Do s3=0,mpoles%max_order
                  ks3=k3+s3

                  If (Mod(ks3,2) == 0) Then
                    Do s2=0,mpoles%max_order-s3
                      ks2=k2+s2

                      If (Mod(ks2,2) == 0) Then
                        Do s1=0,mpoles%max_order-s3-s2
                          ks1=k1+s1

                          mm = mpoles%map(s1,s2,s3)

                          If (Mod(ks1,2) == 0) Then
                            tmpi   = tx*timp*electro%d1_mf(ks1,ks2,ks3)

                            ! energy

                            ewld%engsic = ewld%engsic + imp(mm)*tmpi

                            ! torque

                            !                                      tix    = tix + impx(mm)*tmpi
                            !                                      tiy    = tiy + impy(mm)*tmpi
                            !                                      tiz    = tiz + impz(mm)*tmpi
                          End If
                        End Do
                      End If
                    End Do
                  End If
                End Do
              End If

              tx=-tx
            End Do

            ty=-ty
          End Do

          tz=-tz
        End Do

        ! collect torques due to multipoles selfinteraction

        !        mpoles%torque_x(i)=mpoles%torque_x(i)-exclcoef*tix
        !        mpoles%torque_y(i)=mpoles%torque_y(i)-exclcoef*tiy
        !        mpoles%torque_z(i)=mpoles%torque_z(i)-exclcoef*tiz

      End If

    End Do
    Call gsum(comm,ewld%engsic)
    ewld%engsic = -exclcoef * ewld%engsic / Real(comm%mxnode,wp)

    ! initialise coulombic potential energy and virial

    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    ! set working parameters

    rvolm=twopi/config%volm
    ralph=-0.25_wp/ewld%alpha**2

    ! set scaling constant

    scale=rvolm*r4pie0/electro%eps

    ! Convert config%cell coordinates to fractional coordinates intervalled [0,1)
    ! (bottom left corner of MD config%cell) and stretch over kmaxs in different
    ! directions.  Only the halo (config%natms,config%nlast] has fractional coordinates
    ! outside the [0,1) interval.  In the worst case scenario of one
    ! "effective" link-config%cell per domain and one domain in the MD config%cell only,
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
      txx(i)=electro%kmaxa_r_mf*(rcell(1)*config%parts(i)%xxx+rcell(4)*config%parts(i)%yyy+&
        rcell(7)*config%parts(i)%zzz+0.5_wp)
      tyy(i)=electro%kmaxb_r_mf*(rcell(2)*config%parts(i)%xxx+rcell(5)*config%parts(i)%yyy+&
        rcell(8)*config%parts(i)%zzz+0.5_wp)
      tzz(i)=electro%kmaxc_r_mf*(rcell(3)*config%parts(i)%xxx+rcell(6)*config%parts(i)%yyy+&
        rcell(9)*config%parts(i)%zzz+0.5_wp)

      ! If not DD bound in kmax grid space when .not.neigh%unconditional_update = (ewld%bspline%num_spline_pad == ewld%bspline%num_splines)

      If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines .and. i <= config%natms) Then
        If (txx(i) < electro%ixbm1_r .or. txx(i) > electro%ixtm0_r_mf .or. &
          tyy(i) < electro%iybm1_r_mf .or. tyy(i) > electro%iytm0_r_mf .or. &
          tzz(i) < electro%izbm1_r_mf .or. tzz(i) > electro%iztm0_r_mf) llspl=.false.
      End If

      ixx(i)=Int(txx(i))
      iyy(i)=Int(tyy(i))
      izz(i)=Int(tzz(i))

      ! Detect if a particle is charged and in the MD config%cell or in its positive halo
      ! (t(i) >= -zero_plus) as the B-splines are negative directionally by propagation

      If (tzz(i) >= -zero_plus .and. &
        tyy(i) >= -zero_plus .and. &
        txx(i) >= -zero_plus .and. &
        Maxval(Abs(mpoles%global_frame(:,i))) > zero_plus)  Then
        it(i)=1
      Else
        it(i)=0
      End If
    End Do

    ! Check for breakage of llspl when .not.neigh%unconditional_update = (ewld%bspline1 == ewld%bspline)

    ewld%bspline%num_spline_padded=ewld%bspline%num_spline_pad
    If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines) Then
      Call gcheck(comm,llspl)
      If (.not.llspl) ewld%bspline%num_spline_padded=ewld%bspline%num_splines+1
    End If

    ! construct B-splines for atoms

    ! Call bspgen_mpoles(config%nlast,txx,tyy,tzz,bspx,bspy,bspz,bsddx,bsddy,bsddz,mpoles%n_choose_k,config,ewld,comm)

    Deallocate (txx,tyy,tzz,    Stat = fail(1))
    Deallocate (bspx,bspy,bspz, Stat = fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_mforces allocation failure'
      Call error(0,message)
    End If

    ! zero 3D charge array
    ! DaFT version - only need set local bit to zero

    electro%qqc_local_mf = 0.0_wp
    electro%qtc_local_mf = 0.0_wp

    ! construct 3D charge array
    ! DaFT version - use array that holds only the local data

    If (mpoles%max_order < 5) Then

      Do i=1,config%nlast

        ! If a particle is charged and in the MD config%cell or in its positive halo
        ! (t(i) >= 0) as the B-splines are negative directionally by propagation

        If (it(i) == 1) Then

          ! get the multipoles for site i

          imp=mpoles%global_frame(:,i)

          llb = Max( electro%izb_mf, izz(i) - ewld%bspline%num_splines + 2 )
          llt = Min( electro%izt_mf, izz(i) + 1 )

          kkb = Max( electro%iyb_mf, iyy(i) - ewld%bspline%num_splines + 2 )
          kkt = Min( electro%iyt_mf, iyy(i) + 1 )

          jjb = Max( electro%ixb_mf, ixx(i) - ewld%bspline%num_splines + 2 )
          jjt = Min( electro%ixt_mf, ixx(i) + 1 )

          jjtjjb = jjt - jjb + 1

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mf + 1

            bdz=bsddz(:,l,i)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mf + 1

              bdy=bsddy(:,k,i)

              If (jjtjjb > 0 .and. jjtjjb <= 8) Then

                jj = jjb

                j = ixx(i) - jj + 2

                j_local = jj - electro%ixb_mf + 1

                bdx=bsddx(:,j,i)

                Call explicit_spme_loops        &
                  (0,rcell,bdx,bdy,bdz,imp,impx,impy,impz, &
                  dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3,mpoles,config,ewld)

                electro%qqc_local_mf(j_local,k_local,l_local)   = electro%qqc_local_mf(j_local,k_local,l_local)   + dtp

                electro%qtc_local_mf(1,j_local,k_local,l_local) = electro%qtc_local_mf(1,j_local,k_local,l_local) + tq1
                electro%qtc_local_mf(2,j_local,k_local,l_local) = electro%qtc_local_mf(2,j_local,k_local,l_local) + tq2
                electro%qtc_local_mf(3,j_local,k_local,l_local) = electro%qtc_local_mf(3,j_local,k_local,l_local) + tq3

                counter = 1

                Do While (counter < jjtjjb)
                  j = j - 1

                  j_local = j_local + 1

                  bdx=bsddx(:,j,i)

                  Call explicit_spme_loops     &
                    (0,rcell,bdx,bdy,bdz,imp,impx,impy,impz, &
                    dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3,mpoles,config,ewld)

                  electro%qqc_local_mf(j_local,k_local,l_local)   = electro%qqc_local_mf(j_local,k_local,l_local)   + dtp

                  electro%qtc_local_mf(1,j_local,k_local,l_local) = electro%qtc_local_mf(1,j_local,k_local,l_local) + tq1
                  electro%qtc_local_mf(2,j_local,k_local,l_local) = electro%qtc_local_mf(2,j_local,k_local,l_local) + tq2
                  electro%qtc_local_mf(3,j_local,k_local,l_local) = electro%qtc_local_mf(3,j_local,k_local,l_local) + tq3

                  counter = counter + 1
                End Do

              Else

                Do jj = jjb, jjt
                  j = ixx(i) - jj + 2

                  j_local = jj - electro%ixb_mf + 1

                  bdx=bsddx(:,j,i)

                  Call explicit_spme_loops     &
                    (0,rcell,bdx,bdy,bdz,imp,impx,impy,impz, &
                    dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3,mpoles,config,ewld)

                  electro%qqc_local_mf(j_local,k_local,l_local)   = electro%qqc_local_mf(j_local,k_local,l_local)   + dtp

                  electro%qtc_local_mf(1,j_local,k_local,l_local) = electro%qtc_local_mf(1,j_local,k_local,l_local) + tq1
                  electro%qtc_local_mf(2,j_local,k_local,l_local) = electro%qtc_local_mf(2,j_local,k_local,l_local) + tq2
                  electro%qtc_local_mf(3,j_local,k_local,l_local) = electro%qtc_local_mf(3,j_local,k_local,l_local) + tq3
                End Do

              End If
            End Do
          End Do

        End If

      End Do

    Else

      Do i=1,config%nlast

        ! If a particle is charged and in the MD config%cell or in its positive halo
        ! (t(i) >= 0) as the B-splines are negative directionally by propagation

        If (it(i) == 1) Then

          ! get the multipoles for site i

          imp=mpoles%global_frame(:,i)

          llb = Max( electro%izb_mf, izz(i) - ewld%bspline%num_splines + 2 )
          llt = Min( electro%izt_mf, izz(i) + 1 )

          kkb = Max( electro%iyb_mf, iyy(i) - ewld%bspline%num_splines + 2 )
          kkt = Min( electro%iyt_mf, iyy(i) + 1 )

          jjb = Max( electro%ixb_mf, ixx(i) - ewld%bspline%num_splines + 2 )
          jjt = Min( electro%ixt_mf, ixx(i) + 1 )

          jjtjjb = jjt - jjb + 1

          Do s3 = 0, mpoles%max_order
            Do s2 = 0, mpoles%max_order - s3
              Do s1 = 0, mpoles%max_order - s3 - s2

                mptmp=imp(mpoles%map(s1,s2,s3))

                Do ll = llb, llt
                  l = izz(i) - ll + 2

                  l_local = ll - electro%izb_mf + 1

                  bdz=bsddz(:,l,i)

                  Do kk = kkb, kkt
                    k = iyy(i) - kk + 2

                    k_local = kk - electro%iyb_mf + 1

                    bdy=bsddy(:,k,i)

                    If (jjtjjb > 0 .and. jjtjjb <= 8) Then

                      jj = jjb

                      j = ixx(i) - jj + 2

                      j_local = jj - electro%ixb_mf + 1

                      bdx=bsddx(:,j,i)

                      tmp = mptmp*dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)

                      electro%qqc_local_mf(j_local,k_local,l_local)   = electro%qqc_local_mf(j_local,k_local,l_local)   &
                        + tmp

                      electro%qtc_local_mf(1,j_local,k_local,l_local) = electro%qtc_local_mf(1,j_local,k_local,l_local) &
                        + Real(s1,wp)*tmp
                      electro%qtc_local_mf(2,j_local,k_local,l_local) = electro%qtc_local_mf(2,j_local,k_local,l_local) &
                        + Real(s2,wp)*tmp
                      electro%qtc_local_mf(3,j_local,k_local,l_local) = electro%qtc_local_mf(3,j_local,k_local,l_local) &
                        + Real(s3,wp)*tmp

                      counter = 1

                      Do While (counter < jjtjjb)
                        j = j - 1

                        j_local = j_local + 1

                        bdx=bsddx(:,j,i)

                        tmp = mptmp*dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)

                        electro%qqc_local_mf(j_local,k_local,l_local)   = electro%qqc_local_mf(j_local,k_local,l_local)&
                          + tmp

                        electro%qtc_local_mf(1,j_local,k_local,l_local) = electro%qtc_local_mf(1,j_local,k_local,l_local)&
                          + Real(s1,wp)*tmp
                        electro%qtc_local_mf(2,j_local,k_local,l_local) = electro%qtc_local_mf(2,j_local,k_local,l_local)&
                          + Real(s2,wp)*tmp
                        electro%qtc_local_mf(3,j_local,k_local,l_local) = electro%qtc_local_mf(3,j_local,k_local,l_local)&
                          + Real(s3,wp)*tmp

                        counter = counter + 1
                      End Do

                    Else

                      Do jj = jjb, jjt
                        j = ixx(i) - jj + 2

                        j_local = jj - electro%ixb_mf + 1

                        bdx=bsddx(:,j,i)
                        tmp = mptmp*dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)

                        electro%qqc_local_mf(j_local,k_local,l_local)   = electro%qqc_local_mf(j_local,k_local,l_local)&
                          + tmp

                        electro%qtc_local_mf(1,j_local,k_local,l_local) = electro%qtc_local_mf(1,j_local,k_local,l_local)&
                          + Real(s1,sp)*tmp
                        electro%qtc_local_mf(2,j_local,k_local,l_local) = electro%qtc_local_mf(2,j_local,k_local,l_local)&
                          + Real(s2,sp)*tmp
                        electro%qtc_local_mf(3,j_local,k_local,l_local) = electro%qtc_local_mf(3,j_local,k_local,l_local)&
                          + Real(s3,sp)*tmp
                      End Do

                    End If
                  End Do
                End Do

              End Do
            End Do
          End Do

        End If

      End Do

    End If

    ! load charge array into complex array for FFT

    electro%qqq_local_mf=Cmplx(electro%qqc_local_mf , Kind = wp)

    electro%qt1_local_mf=Cmplx(electro%qtc_local_mf(1,:,:,:), Kind = wp)
    electro%qt2_local_mf=Cmplx(electro%qtc_local_mf(2,:,:,:), Kind = wp)
    electro%qt3_local_mf=Cmplx(electro%qtc_local_mf(3,:,:,:), Kind = wp)

    ! calculating inverse 3D FFT of generalized multipolar array (in place)

    Call pfft(electro%qqq_local_mf,electro%pfft_work_mf,electro%context_mf,1)

    Call pfft(electro%qt1_local_mf,electro%pfft_work_mf,electro%context_mf,1)
    Call pfft(electro%qt2_local_mf,electro%pfft_work_mf,electro%context_mf,1)
    Call pfft(electro%qt3_local_mf,electro%pfft_work_mf,electro%context_mf,1)

    ! set reciprocal space cutoff

    Call dcell(rcell,celprp)

    rcpcut=0.5_wp*Min(electro%kmaxa_r_mf*celprp(7),electro%kmaxb_r_mf*celprp(8),electro%kmaxc_r_mf*celprp(9))
    rcpcut=rcpcut*1.05_wp*twopi
    rcpct2=rcpcut**2

    ! initialise temporary stress tensor

    strs = 0.0_wp

    ! calculate convolution of charge array with gaussian function
    ! DaFT Version - only loop over the local stuff

    Do l_local=1,electro%block_z_mf
      l=electro%index_z_mf(l_local)

      ll=l-1
      If (l > ewld%kspace%k_vec_dim(3)/2) ll=ll-ewld%kspace%k_vec_dim(3)
      tmp=twopi*Real(ll,wp)

      rkx1=tmp*rcell(3)
      rky1=tmp*rcell(6)
      rkz1=tmp*rcell(9)

      bb3=Real( electro%bscz_mf(l)*Conjg(electro%bscz_mf(l)),wp )

      Do k_local=1,electro%block_y_mf
        k=electro%index_y_mf(k_local)

        kk=k-1
        If (k > ewld%kspace%k_vec_dim(2)/2) kk=kk-ewld%kspace%k_vec_dim(2)
        tmp=twopi*Real(kk,wp)

        rkx2=rkx1+tmp*rcell(2)
        rky2=rky1+tmp*rcell(5)
        rkz2=rkz1+tmp*rcell(8)

        bb2=bb3*Real( electro%bscy_mf(k)*Conjg(electro%bscy_mf(k)),wp )

        Do j_local=1,electro%block_x_mf
          j=electro%index_x_mf(j_local)

          jj=j-1
          If (j > ewld%kspace%k_vec_dim(1)/2) jj=jj-ewld%kspace%k_vec_dim(1)
          tmp=twopi*Real(jj,wp)

          rkx3=rkx2+tmp*rcell(1)
          rky3=rky2+tmp*rcell(4)
          rkz3=rkz2+tmp*rcell(7)

          bb1=bb2*Real( electro%bscx_mf(j)*Conjg(electro%bscx_mf(j)),wp )

          rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

          !=================================================================
          ! For higher order contributions to the stress tensor

          rrkxy=0.0_wp ; rrkxz=0.0_wp
          rrkyx=0.0_wp ; rrkyz=0.0_wp
          rrkzx=0.0_wp ; rrkzy=0.0_wp

          If (rkx3 > 0.0_wp) Then

            rrkyx=rky3/rkx3; rrkzx=rkz3/rkx3

          End If

          If (rky3 > 0.0_wp) Then

            rrkxy=rkx3/rky3; rrkzy=rkz3/rky3

          End If

          If (rkz3 > 0.0_wp) Then

            rrkxz=rkx3/rkz3; rrkyz=rky3/rkz3

          End If

          !====================================================================

          If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

            ! For  monopole contribution to the stress tensor

            vterm=bb1*Exp(ralph*rksq)/rksq*electro%qqq_local_mf(j_local,k_local,l_local)
            akv=2.0_wp*(1.0_wp/rksq-ralph)*Real( vterm*Conjg(electro%qqq_local_mf(j_local,k_local,l_local)),wp )

            !=======================================================================
            ! For higher order contributions to the stress tensor

            pterm=2.0_wp*bb1*Exp(ralph*rksq)/rksq*Conjg(electro%qqq_local_mf(j_local,k_local,l_local))
            sq1=Real( pterm*electro%qt1_local_mf(j_local,k_local,l_local),wp )
            sq2=Real( pterm*electro%qt2_local_mf(j_local,k_local,l_local),wp )
            sq3=Real( pterm*electro%qt3_local_mf(j_local,k_local,l_local),wp )

            !=======================================================================

            strs(1)=strs(1)-rkx3*rkx3*akv+sq1
            strs(5)=strs(5)-rky3*rky3*akv+sq2
            strs(9)=strs(9)-rkz3*rkz3*akv+sq3
            strs(2)=strs(2)-rkx3*rky3*akv+rrkxy*sq2
            strs(3)=strs(3)-rkx3*rkz3*akv+rrkxz*sq3
            strs(4)=strs(4)-rkx3*rky3*akv+rrkyx*sq1
            strs(6)=strs(6)-rky3*rkz3*akv+rrkyz*sq3
            strs(7)=strs(7)-rkx3*rkz3*akv+rrkzx*sq1
            strs(8)=strs(8)-rky3*rkz3*akv+rrkzy*sq2

            electro%qqq_local_mf(j_local,k_local,l_local)=vterm

          Else

            electro%qqq_local_mf(j_local,k_local,l_local)=(0.0_wp,0.0_wp)

          End If

        End Do

      End Do

    End Do

    ! as only looped over local stuff, we need to gsum strs

    Call gsum(comm,strs)

    ! scale strs and distribute per node

    strs = strs * scale / Real(comm%mxnode,wp)

    ! calculate atomic energy

    Call pfft(electro%qqq_local_mf,electro%pfft_work_mf,electro%context_mf,-1)

    eng = 0.0_wp
    Do l=1,electro%block_z_mf
      Do k=1,electro%block_y_mf
        Do j=1,electro%block_x_mf
          qqc_tmp=Real(electro%qqq_local_mf(j,k,l),wp)
          eng=eng+electro%qqc_local_mf(j,k,l)*qqc_tmp
          electro%qqc_local_mf(j,k,l)=qqc_tmp
        End Do
      End Do
    End Do

    ! as only looped over local stuff, we need to gsum the eng

    Call gsum(comm,eng)

    ! scale eng and distribute per node

    eng = eng * scale / Real(comm%mxnode,wp)

    ! Second part of the monopole contribution to the stress tensor
    ! calculate stress tensor (symmetrical, per node)

    strs      = strs      + eng*uni
    stress(1) = stress(1) + strs(1)
    stress(2) = stress(2) + 0.5_wp * (strs(2) + strs(4))
    stress(3) = stress(3) + 0.5_wp * (strs(3) + strs(7))
    stress(4) = stress(4) + 0.5_wp * (strs(2) + strs(4))
    stress(5) = stress(5) + strs(5)
    stress(6) = stress(6) + 0.5_wp * (strs(6) + strs(8))
    stress(7) = stress(7) + 0.5_wp * (strs(3) + strs(7))
    stress(8) = stress(8) + 0.5_wp * (strs(6) + strs(8))
    stress(9) = stress(9) + strs(9)

    ! distribute energy and virial terms (per node)

    engcpe_rc = eng + ewld%engsic
    vircpe_rc = -(strs(1)+strs(5)+strs(9))

    ! infrequent calculations copying

    If (ewld%l_cp) Then
      ewld%e_rc=engcpe_rc
      ewld%v_rc=vircpe_rc
      ewld%s_rc=strs
    End If

    ! calculate atomic forces

    Call spme_mforces(rcell,scale,ixx,iyy,izz,bsddx,bsddy,bsddz,electro%qqc_local_mf, &
      electro%ixb_mf,electro%ixt_mf,electro%iyb_mf,electro%iyt_mf,electro%izb_mf,electro%izt_mf,mpoles,domain)

    Deallocate (ixx,iyy,izz,it,    Stat = fail(1))
    Deallocate (bdx,bdy,bdz,       Stat = fail(2))
    Deallocate (bsddx,bsddy,bsddz, Stat = fail(3))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_mforces deallocation failure'
      Call error(0,message)
    End If

  Contains

    Subroutine spme_mforces(rcell,scale,ixx,iyy,izz,bsddx,bsddy,bsddz, &
      qqc_local,ixb,ixt,iyb,iyt,izb,izt,mpoles,domain)

      !!----------------------------------------------------------------------!
      !!
      !! dl_poly_4 subroutine for calculating coulombic forces due to
      !! multipolar interactions in a periodic system using smooth particle
      !! mesh ewald method (fourier part)
      !!
      !! Note: qqc_local is shifted from its definition from above
      !!       and therefore there is no need for periodic images (!!)
      !!
      !! copyright - daresbury laboratory
      !! author    - w.smith & i.t.todorov february 2016
      !! amended   - h.a.boateng may 2014
      !!
      !!----------------------------------------------------------------------!

      Integer,           Intent( In    ) :: ixx(1:config%mxatms),iyy(1:config%mxatms),izz(1:config%mxatms), &
        ixb,ixt, iyb,iyt, izb,izt
      Real( Kind = wp ), Intent( In    ) :: scale,rcell(1:9),                &
        bsddx(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
        bsddy(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
        bsddz(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
        qqc_local( ixb:ixt, iyb:iyt, izb:izt )
      Type( mpole_type ), Intent( InOut ) :: mpoles
      Type( domains_type ), Intent( In    ) :: domain

      Integer           :: fail(1:2), delspl, ixdb,iydb,izdb,ixdt,iydt,izdt, &
        i,j,k,l, jj,kk,ll, s1,s2,s3, mm
      Real( Kind = wp ) :: tmp,gmp,fff(0:3),fix,fiy,fiz,qsum, &
        tmpi,tix,tiy,tiz,dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3

      Real( Kind = wp ) :: imp(1:mpoles%max_mpoles)
      Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)

      Real( Kind = wp ), Dimension( : ),       Allocatable :: bdx,bdy,bdz
      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: qqc_domain

      ! Define extended ranges for the domain = local + halo slice and allocate

      ixdb = ixb - ewld%bspline%num_spline_padded
      iydb = iyb - ewld%bspline%num_spline_padded
      izdb = izb - ewld%bspline%num_spline_padded

      delspl = ewld%bspline%num_spline_padded - ewld%bspline%num_splines

      ixdt = ixt + delspl
      iydt = iyt + delspl
      izdt = izt + delspl

      fail=0
      Allocate (bdx(0:ewld%bspline%num_splines),bdy(0:ewld%bspline%num_splines),&
        & bdz(0:ewld%bspline%num_splines), Stat = fail(1))
      Allocate (qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt ), Stat = fail(2))
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'spme_mforces allocation failure'
        Call error(0,message)
      End If

      Call exchange_grid(ixb , ixt , iyb , iyt , izb , izt , qqc_local, &
        ixdb, iydb, izdb, ixdt, iydt, izdt, qqc_domain, &
        domain, comm) !, ewld

      tmp=-2.0_wp*scale

      fff=0.0_wp

      If (mpoles%max_order < 5) Then

        Do i=1,config%natms

          ! get the multipoles for site i

          imp=mpoles%global_frame(:,i)

          If (Maxval(Abs(imp)) > zero_plus) Then

            ! scale imp multipoles

            imp=tmp*imp

            ! get the components for site i infinitesimal rotations

            impx=scale*mpoles%rotation_x(:,i)
            impy=scale*mpoles%rotation_y(:,i)
            impz=scale*mpoles%rotation_z(:,i)

            ! initialise forces & torques

            fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp
            tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp

            Do l=1,ewld%bspline%num_splines
              ll=izz(i)-l+2

              bdz=bsddz(:,l,i)

              Do k=1,ewld%bspline%num_splines
                kk=iyy(i)-k+2

                bdy=bsddy(:,k,i)

                Do j=1,ewld%bspline%num_splines
                  jj=ixx(i)-j+2

                  bdx=bsddx(:,j,i)

                  qsum=qqc_domain(jj,kk,ll)

                  Call explicit_spme_loops      &
                    (1,rcell,bdx,bdy,bdz,imp,impx,impy,impz, &
                    dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3,mpoles,config,ewld)

                  ! force

                  fix = fix + qsum*dt1
                  fiy = fiy + qsum*dt2
                  fiz = fiz + qsum*dt3

                  ! torque

                  tix = tix + qsum*td1
                  tiy = tiy + qsum*td2
                  tiz = tiz + qsum*td3
                End Do
              End Do
            End Do

            ! accumulate forces

            fff(0)=fff(0)+1.0_wp
            fff(1)=fff(1)+fix
            fff(2)=fff(2)+fiy
            fff(3)=fff(3)+fiz

            ! load forces

            config%parts(i)%fxx=config%parts(i)%fxx+fix
            config%parts(i)%fyy=config%parts(i)%fyy+fiy
            config%parts(i)%fzz=config%parts(i)%fzz+fiz

            ! and torque (ITT - is sum of all torques zero? I.e. does SPME machinery generate non-0 torque)

            mpoles%torque_x(i)=mpoles%torque_x(i)+0.5_wp*tix
            mpoles%torque_y(i)=mpoles%torque_y(i)+0.5_wp*tiy
            mpoles%torque_z(i)=mpoles%torque_z(i)+0.5_wp*tiz

            ! infrequent calculations copying

            If (ewld%l_cp) Then
              ewld%fcx(i)=ewld%fcx(i)+fix
              ewld%fcy(i)=ewld%fcy(i)+fiy
              ewld%fcz(i)=ewld%fcz(i)+fiz
            End If

          End If

        End Do

      Else

        Do i=1,config%natms

          ! get the multipoles for site i

          imp=mpoles%global_frame(:,i)

          If (Maxval(Abs(imp)) > zero_plus) Then

            ! scale imp multipoles

            imp=tmp*imp

            ! get the components for site i infinitesimal rotations

            impx=scale*mpoles%rotation_x(:,i)
            impy=scale*mpoles%rotation_y(:,i)
            impz=scale*mpoles%rotation_z(:,i)

            ! initialise forces & torques

            fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp
            tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp

            Do l=1,ewld%bspline%num_splines
              ll=izz(i)-l+2

              bdz=bsddz(:,l,i)

              Do k=1,ewld%bspline%num_splines
                kk=iyy(i)-k+2

                bdy=bsddy(:,k,i)

                Do j=1,ewld%bspline%num_splines
                  jj=ixx(i)-j+2

                  bdx=bsddx(:,j,i)

                  qsum=qqc_domain(jj,kk,ll)

                  Do s3 = 0, mpoles%max_order
                    Do s2 = 0, mpoles%max_order - s3
                      Do s1 = 0, mpoles%max_order - s3 - s2
                        mm = mpoles%map(s1,s2,s3)

                        ! force

                        gmp = qsum*imp(mm)

                        fix = fix + gmp*dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)
                        fiy = fiy + gmp*dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)
                        fiz = fiz + gmp*dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)

                        ! torque

                        tmpi = qsum*dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz,mpoles%n_choose_k,config,ewld)

                        tix = tix + impx(mm)*tmpi
                        tiy = tiy + impy(mm)*tmpi
                        tiz = tiz + impz(mm)*tmpi
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End Do

            ! accumulate forces

            fff(0)=fff(0)+1.0_wp
            fff(1)=fff(1)+fix
            fff(2)=fff(2)+fiy
            fff(3)=fff(3)+fiz

            ! load forces

            config%parts(i)%fxx=config%parts(i)%fxx+fix
            config%parts(i)%fyy=config%parts(i)%fyy+fiy
            config%parts(i)%fzz=config%parts(i)%fzz+fiz

            ! and torque

            mpoles%torque_x(i)=mpoles%torque_x(i)+0.5_wp*tix
            mpoles%torque_y(i)=mpoles%torque_y(i)+0.5_wp*tiy
            mpoles%torque_z(i)=mpoles%torque_z(i)+0.5_wp*tiz

            ! infrequent calculations copying

            If (ewld%l_cp) Then
              ewld%fcx(i)=ewld%fcx(i)+fix
              ewld%fcy(i)=ewld%fcy(i)+fiy
              ewld%fcz(i)=ewld%fcz(i)+fiz
            End If

          End If

        End Do

      End If

      ! remove COM drift arising from SPME approximations

      Call gsum(comm,fff)
      If (fff(0) > zero_plus) Then
        fff(1:3)=fff(1:3)/fff(0)

        Do i=1,config%natms
          imp=mpoles%global_frame(:,i)
          If (Maxval(Abs(imp)) > zero_plus) Then

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

      Deallocate (bdx,bdy,bdz, Stat = fail(1))
      Deallocate (qqc_domain,  Stat = fail(2))
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'spme_mforces dealocation failure'
        Call error(0,message)
      End If

    End Subroutine spme_mforces

  End Subroutine ewald_spme_mforces

  Subroutine ewald_spme_mforces_d(engcpe_rc,vircpe_rc,stress,ewld,mpoles,electro, &
    domain,config,comm)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! due to multipolar interactions in a periodic system using the smooth
    !! particle mesh ewald method for multipoles
    !!
    !! This version allows for extension to arbitrary order
    !!
    !! Note: (fourier) reciprocal space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - h.a.boateng february 2016
    !! amended   - i.t.todorov february 2015
    !!
    !!----------------------------------------------------------------------!

    Real( Kind = wp ), Intent(   Out ) :: engcpe_rc,vircpe_rc
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Type( ewald_type ), Intent( InOut ) :: ewld
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( electrostatic_type ), Intent( InOut    ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm
    Type( configuration_type ),                       Intent( InOut ) :: config


    Logical              :: llspl=.true.
    Integer              :: fail(1:4), i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt

    Real( Kind = wp )    :: det,rcell(1:9),celprp(1:10),ralph,rvolm,scale,       &
      rcpcut,rcpct2,strs(1:9),eng,akv,tmp,bb1,bb2,bb3,     &
      rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3, &
      rrkxy,rrkxz,rrkyx,rrkyz,rrkzx,rrkzy, sq1,sq2,sq3,    &
      dtp,tq1,tq2,tq3,                                     &
      imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10,  &
      bdy0,bdy1,bdy2,bdz0,bdz1,bdz2,                       &
      ka11,kb22,kc33,ka11sq,kb22sq,kc33sq,kakb,kakc,kbkc

    Real( Kind = wp )    :: imp(1:mpoles%max_mpoles)

    Complex( Kind = wp ) :: vterm,pterm

    ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter   :: &
      uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

    ! B-spline coefficients

    Complex( Kind = wp ), Dimension( : ),       Allocatable       :: ww1,ww2,ww3

    Real( Kind = wp ),    Dimension( : ),       Allocatable       :: csp
    Real( Kind = wp ),    Dimension( : ),       Allocatable       :: txx,tyy,tzz
    Integer,              Dimension( : ),       Allocatable       :: ixx,iyy,izz,it
    Real( Kind = wp ),    Dimension( :,:,: ),   Allocatable       :: bsddx,bsddy,bsddz
    Real( Kind = wp ),    Dimension( :,: ),     Allocatable       :: bspx,bspy,bspz
    Real( Kind = wp ),    Dimension( : ),       Allocatable       :: bdx,bdy,bdz


    ! temporary qqc

    Real( Kind = wp )    :: qqc_tmp


    ! DaFT arrays local indices

    Integer              :: j_local, k_local, l_local
    Character( Len = 256 ) :: message


    fail=0
    If (electro%newjob_mfd) Then
      electro%newjob_mfd = .false.

!!! BEGIN DD SPME VARIABLES
      ! 3D charge array construction (bottom and top) indices

      electro%ixb_mfd=domain%idx*(ewld%kspace%k_vec_dim(1)/domain%nx)+1
      electro%ixt_mfd=(domain%idx+1)*(ewld%kspace%k_vec_dim(1)/domain%nx)
      electro%iyb_mfd=domain%idy*(ewld%kspace%k_vec_dim(2)/domain%ny)+1
      electro%iyt_mfd=(domain%idy+1)*(ewld%kspace%k_vec_dim(2)/domain%ny)
      electro%izb_mfd=domain%idz*(ewld%kspace%k_vec_dim(3)/domain%nz)+1
      electro%izt_mfd=(domain%idz+1)*(ewld%kspace%k_vec_dim(3)/domain%nz)

      electro%ixbm1_r_mfd=Real(electro%ixb_mfd-1,wp)
      electro%ixtm0_r_mfd=Nearest( Real(electro%ixt_mfd,wp) , -1.0_wp )
      electro%iybm1_r_mfd=Real(electro%iyb_mfd-1,wp)
      electro%iytm0_r_mfd=Nearest( Real(electro%iyt_mfd,wp) , -1.0_wp )
      electro%izbm1_r_mfd=Real(electro%izb_mfd-1,wp)
      electro%iztm0_r_mfd=Nearest( Real(electro%izt_mfd,wp) , -1.0_wp )

      ! Real values of kmax vectors

      electro%kmaxa_r_mfd=Real(ewld%kspace%k_vec_dim(1),wp)
      electro%kmaxb_r_mfd=Real(ewld%kspace%k_vec_dim(2),wp)
      electro%kmaxc_r_mfd=Real(ewld%kspace%k_vec_dim(3),wp)

!!! END DD SPME VARIABLES

!!! BEGIN CARDINAL B-SPLINES SET-UP
      call bspline_coeffs_gen(ewld%kspace, ewld%bsplne)
!!! END CARDINAL B-SPLINES SET-UP

      ! set up the parallel fft and useful related quantities

      Call initialize_fft( 3, (/ ewld%kspace%k_vec_dim(1), ewld%kspace%k_vec_dim(2), ewld%kspace%k_vec_dim(3) /), &
        (/ domain%nx, domain%ny, domain%nz /), (/ domain%idx, domain%idy, domain%idz /),   &
        (/ electro%block_x_mfd, electro%block_y_mfd, electro%block_z_mfd /),               &
        comm%comm, electro%context_mfd )

      ! set up the indexing arrays for each dimension (NOT deallocated manually)

      Allocate ( electro%index_x_mfd( 1:electro%block_x_mfd ), Stat = fail(1) )
      Allocate ( electro%index_y_mfd( 1:electro%block_y_mfd ), Stat = fail(2) )
      Allocate ( electro%index_z_mfd( 1:electro%block_z_mfd ), Stat = fail(3) )
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'SPME index arrays allocation failure'
        Call error(0,message)
      End If

      Call pfft_indices( ewld%kspace%k_vec_dim(1), electro%block_x_mfd, domain%idx, domain%nx, electro%index_x_mfd )
      Call pfft_indices( ewld%kspace%k_vec_dim(2), electro%block_y_mfd, domain%idy, domain%ny, electro%index_y_mfd )
      Call pfft_indices( ewld%kspace%k_vec_dim(3), electro%block_z_mfd, domain%idz, domain%nz, electro%index_z_mfd )

      ! workspace arrays for DaFT

      Allocate ( electro%qqc_local_mfd(1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ), &
        & Stat = fail(1) )
      Allocate ( electro%qqq_local_mfd(1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ), &
        & Stat = fail(2) )
      Allocate ( electro%qtc_local_mfd(1:3, 1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ), &
        & electro%qt1_local_mfd( 1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ),      &
        & electro%qt2_local_mfd( 1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ),      &
        & electro%qt3_local_mfd( 1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ), Stat = fail(3) )
      Allocate ( electro%pfft_work_mfd( 1:electro%block_x_mfd, 1:electro%block_y_mfd, 1:electro%block_z_mfd ), Stat = fail(4) )
      If (Any(fail > 0)) call error_alloc('SPME DaFT workspace arrays','ewald_spme_mforces_d')

      !!! END DAFT SET-UP
    End If

    Allocate (txx(1:config%mxatms),tyy(1:config%mxatms),tzz(1:config%mxatms),                            Stat = fail(1))
    Allocate (ixx(1:config%mxatms),iyy(1:config%mxatms),izz(1:config%mxatms),it(1:config%mxatms),               Stat = fail(2))
    Allocate ( &
      & bdx(0:ewld%bspline%num_splines),&
      & bsddx(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
      & bdy(0:ewld%bspline%num_splines),&
      & bsddy(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
      & bdz(0:ewld%bspline%num_splines),&
      & bsddz(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
      & Stat = fail(3))
    Allocate ( &
      & bspx(1:ewld%bspline%num_splines,1:config%mxatms), &
      & bspy(1:ewld%bspline%num_splines,1:config%mxatms), &
      & bspz(1:ewld%bspline%num_splines,1:config%mxatms), &
      & Stat = fail(4))
    If (Any(fail > 0)) call error_alloc('SPME arrays','ewald_spme_mforces_d')
    ! initialise coulombic potential energy and virial

    ewld%engsic    = 0.0_wp
    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    ! set working parameters

    rvolm=twopi/config%volm
    ralph=-0.25_wp/ewld%alpha**2

    ! set scaling constant

    scale=rvolm*r4pie0/electro%eps

    ! Convert config%cell coordinates to fractional coordinates intervalled [0,1)
    ! (bottom left corner of MD config%cell) and stretch over kmaxs in different
    ! directions.  Only the halo (config%natms,config%nlast] has fractional coordinates
    ! outside the [0,1) interval.  In the worst case scenario of one
    ! "effective" link-config%cell per domain and one domain in the MD config%cell only,
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
      txx(i)=electro%kmaxa_r_mfd*(rcell(1)*config%parts(i)%xxx+rcell(4)*config%parts(i)%yyy+&
        rcell(7)*config%parts(i)%zzz+0.5_wp)
      tyy(i)=electro%kmaxb_r_mfd*(rcell(2)*config%parts(i)%xxx+rcell(5)*config%parts(i)%yyy+&
        rcell(8)*config%parts(i)%zzz+0.5_wp)
      tzz(i)=electro%kmaxc_r_mfd*(rcell(3)*config%parts(i)%xxx+rcell(6)*config%parts(i)%yyy+&
        rcell(9)*config%parts(i)%zzz+0.5_wp)

      ! If not DD bound in kmax grid space when .not.neigh%unconditional_update = (ewld%bspline%num_spline_pad == ewld%bspline)

      If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines .and. i <= config%natms) Then
        If (txx(i) < electro%ixbm1_r_mfd .or. txx(i) > electro%ixtm0_r_mfd .or. &
          tyy(i) < electro%iybm1_r_mfd .or. tyy(i) > electro%iytm0_r_mfd .or. &
          tzz(i) < electro%izbm1_r_mfd .or. tzz(i) > electro%iztm0_r_mfd) llspl=.false.
      End If

      ixx(i)=Int(txx(i))
      iyy(i)=Int(tyy(i))
      izz(i)=Int(tzz(i))

      ! Detect if a particle is charged and in the MD config%cell or in its positive halo
      ! (t(i) >= -zero_plus) as the B-splines are negative directionally by propagation

      If (tzz(i) >= -zero_plus .and. &
        tyy(i) >= -zero_plus .and. &
        txx(i) >= -zero_plus .and. &
        Maxval(Abs(mpoles%global_frame(:,i))) > zero_plus)  Then
        it(i)=1
      Else
        it(i)=0
      End If
    End Do

    ! Check for breakage of llspl when .not.neigh%unconditional_update = (ewld%bspline%num_spline_pad == ewld%bspline)

    ewld%bspline%num_spline_padded=ewld%bspline%num_spline_pad
    If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines) Then
      Call gcheck(comm,llspl)
      If (.not.llspl) ewld%bspline%num_spline_padded=ewld%bspline%num_splines+1
    End If

    ! construct B-splines for atoms

    ! Call bspgen_mpoles(config%nlast,txx,tyy,tzz,bspx,bspy,bspz,bsddx,bsddy,bsddz,mpoles%n_choose_k,config,ewld,comm)

    Deallocate (txx,tyy,tzz,    Stat = fail(1))
    Deallocate (bspx,bspy,bspz, Stat = fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_mforces allocation failure'
      Call error(0,message)
    End If

    ! zero 3D charge array
    ! DaFT version - only need set local bit to zero

    electro%qqc_local_mfd = 0.0_wp
    electro%qtc_local_mfd = 0.0_wp

    ! construct 3D charge array

    ka11 = electro%kmaxa_r_mfd*rcell(1)
    kb22 = electro%kmaxb_r_mfd*rcell(5)
    kc33 = electro%kmaxc_r_mfd*rcell(9)

    ka11sq = ka11*ka11; kb22sq = kb22*kb22; kc33sq = kc33*kc33
    kakb = ka11*kb22; kakc = ka11*kc33; kbkc = kb22*kc33

    ! DaFT version - use array that holds only the local data

    Do i=1,config%nlast

      ! If a particle is charged and in the MD config%cell or in its positive halo
      ! (t(i) >= 0) as the B-splines are negative directionally by propagation

      If (it(i) == 1) Then

        ! get the multipoles for site i

        imp=mpoles%global_frame(:,i)

        imp1=imp(1)

        If (mpoles%max_order >= 1) Then
          imp2=imp(2) ; imp3=imp(3) ; imp4=imp(4)
        End If

        If (mpoles%max_order == 2) Then
          imp5=imp(5) ; imp6=imp(6) ; imp7=imp(7)
          imp8=imp(8) ; imp9=imp(9) ; imp10=imp(10)
        End If

        llb = Max( electro%izb_mfd, izz(i) - ewld%bspline%num_splines + 2 )
        llt = Min( electro%izt_mfd, izz(i) + 1 )

        kkb = Max( electro%iyb_mfd, iyy(i) - ewld%bspline%num_splines + 2 )
        kkt = Min( electro%iyt_mfd, iyy(i) + 1 )

        jjb = Max( electro%ixb_mfd, ixx(i) - ewld%bspline%num_splines + 2 )
        jjt = Min( electro%ixt_mfd, ixx(i) + 1 )

        Select Case( jjt - jjb + 1 )

        Case Default

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              Do jj = jjb, jjt
                j = ixx(i) - jj + 2

                j_local = jj - electro%ixb_mfd + 1

                bdx=bsddx(:,j,i)

                dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

                dtp = imp1*bdx(0)*bdy0*bdz0

                If (mpoles%max_order >= 1 ) Then
                  tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                  tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                  tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
                End If

                If (mpoles%max_order == 2) Then
                  tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                  tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                  tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                  tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                  tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                  tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
                End If

                electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

                electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
                electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
                electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
              End Do
            End Do
          End Do

        Case( 0 )

        Case( 1 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb

              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 2 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 3 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !3
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 4 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !3
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !4
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 5 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !3
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !4
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !5
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 6 )


          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !3
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !4
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !5
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !6
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 7 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !3
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !4
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !5
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !6
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !7
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        Case( 8 )

          Do ll = llb, llt
            l = izz(i) - ll + 2

            l_local = ll - electro%izb_mfd + 1

            bdz=bsddz(:,l,i)
            bdz0=bdz(0); bdz1=bdz(1); bdz2=bdz(2)

            Do kk = kkb, kkt
              k = iyy(i) - kk + 2

              k_local = kk - electro%iyb_mfd + 1

              bdy=bsddy(:,k,i)
              bdy0=bdy(0); bdy1=bdy(1); bdy2=bdy(2)

              jj = jjb
              !1
              j = ixx(i) - jj + 2

              j_local = jj - electro%ixb_mfd + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !2
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !3
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !4
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !5
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !6
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !7
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3

              !8
              j = j - 1

              j_local = j_local + 1

              bdx=bsddx(:,j,i)

              dtp=0.0_wp; tq1=0.0_wp; tq2=0.0_wp; tq3=0.0_wp

              dtp = imp1*bdx(0)*bdy0*bdz0

              If (mpoles%max_order >= 1 ) Then
                tmp = imp2*ka11*bdx(1)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+tmp
                tmp = imp3*kb22*bdx(0)*bdy1*bdz0; dtp=dtp+tmp; tq2=tq2+tmp
                tmp = imp4*kc33*bdx(0)*bdy0*bdz1; dtp=dtp+tmp; tq3=tq3+tmp
              End If

              If (mpoles%max_order == 2) Then
                tmp = imp5*ka11sq  *bdx(2)*bdy0*bdz0; dtp=dtp+tmp; tq1=tq1+2.0_wp*tmp
                tmp = imp8*kb22sq  *bdx(0)*bdy2*bdz0; dtp=dtp+tmp; tq2=tq2+2.0_wp*tmp
                tmp = imp10*kc33sq *bdx(0)*bdy0*bdz2; dtp=dtp+tmp; tq3=tq3+2.0_wp*tmp

                tmp = imp6*kakb*bdx(1)*bdy1*bdz0; dtp=dtp+tmp; tq1=tq1+tmp; tq2=tq2+tmp
                tmp = imp7*kakc*bdx(1)*bdy0*bdz1; dtp=dtp+tmp; tq1=tq1+tmp; tq3=tq3+tmp
                tmp = imp9*kbkc*bdx(0)*bdy1*bdz1; dtp=dtp+tmp; tq2=tq2+tmp; tq3=tq3+tmp
              End If

              electro%qqc_local_mfd(j_local,k_local,l_local)   = electro%qqc_local_mfd(j_local,k_local,l_local)   + dtp

              electro%qtc_local_mfd(1,j_local,k_local,l_local) = electro%qtc_local_mfd(1,j_local,k_local,l_local) + tq1
              electro%qtc_local_mfd(2,j_local,k_local,l_local) = electro%qtc_local_mfd(2,j_local,k_local,l_local) + tq2
              electro%qtc_local_mfd(3,j_local,k_local,l_local) = electro%qtc_local_mfd(3,j_local,k_local,l_local) + tq3
            End Do
          End Do

        End Select

      End If

    End Do

    ! load charge array into complex array for FFT

    electro%qqq_local_mfd=Cmplx(electro%qqc_local_mfd , Kind = wp)

    electro%qt1_local_mfd=Cmplx(electro%qtc_local_mfd(1,:,:,:), Kind = wp)
    electro%qt2_local_mfd=Cmplx(electro%qtc_local_mfd(2,:,:,:), Kind = wp)
    electro%qt3_local_mfd=Cmplx(electro%qtc_local_mfd(3,:,:,:), Kind = wp)

    ! calculate inverse 3D FFT of charge array (in place)

    Call pfft(electro%qqq_local_mfd,electro%pfft_work_mfd,electro%context_mfd,1)

    Call pfft(electro%qt1_local_mfd,electro%pfft_work_mfd,electro%context_mfd,1)
    Call pfft(electro%qt2_local_mfd,electro%pfft_work_mfd,electro%context_mfd,1)
    Call pfft(electro%qt3_local_mfd,electro%pfft_work_mfd,electro%context_mfd,1)

    ! set reciprocal space cutoff

    Call dcell(rcell,celprp)

    rcpcut=0.5_wp*Min(electro%kmaxa_r_mfd*celprp(7),electro%kmaxb_r_mfd*celprp(8),electro%kmaxc_r_mfd*celprp(9))
    rcpcut=rcpcut*1.05_wp*twopi
    rcpct2=rcpcut**2

    ! initialise temporary stress tensor

    strs = 0.0_wp

    ! calculate convolution of charge array with gaussian function
    ! DaFT Version - only loop over the local stuff

    Do l_local=1,electro%block_z_mfd
      l=electro%index_z_mfd(l_local)

      ll=l-1
      If (l > ewld%kspace%k_vec_dim(3)/2) ll=ll-ewld%kspace%k_vec_dim(3)
      tmp=twopi*Real(ll,wp)

      rkx1=tmp*rcell(3)
      rky1=tmp*rcell(6)
      rkz1=tmp*rcell(9)

      bb3=Real( electro%bscz_mfd(l)*Conjg(electro%bscz_mfd(l)),wp )

      Do k_local=1,electro%block_y_mfd
        k=electro%index_y_mfd(k_local)

        kk=k-1
        If (k > ewld%kspace%k_vec_dim(2)/2) kk=kk-ewld%kspace%k_vec_dim(2)
        tmp=twopi*Real(kk,wp)

        rkx2=rkx1+tmp*rcell(2)
        rky2=rky1+tmp*rcell(5)
        rkz2=rkz1+tmp*rcell(8)

        bb2=bb3*Real( electro%bscy_mfd(k)*Conjg(electro%bscy_mfd(k)),wp )

        Do j_local=1,electro%block_x_mfd
          j=electro%index_x_mfd(j_local)

          jj=j-1
          If (j > ewld%kspace%k_vec_dim(1)/2) jj=jj-ewld%kspace%k_vec_dim(1)
          tmp=twopi*Real(jj,wp)

          rkx3=rkx2+tmp*rcell(1)
          rky3=rky2+tmp*rcell(4)
          rkz3=rkz2+tmp*rcell(7)

          bb1=bb2*Real( electro%bscx_mfd(j)*Conjg(electro%bscx_mfd(j)),wp )

          rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

          !=================================================================
          ! For higher order contributions to the stress tensor

          rrkxy=0.0_wp ; rrkxz=0.0_wp
          rrkyx=0.0_wp ; rrkyz=0.0_wp
          rrkzx=0.0_wp ; rrkzy=0.0_wp

          If (rkx3 > 0.0_wp) Then

            rrkyx=rky3/rkx3; rrkzx=rkz3/rkx3

          End If

          If (rky3 > 0.0_wp) Then

            rrkxy=rkx3/rky3; rrkzy=rkz3/rky3

          End If

          If (rkz3 > 0.0_wp) Then

            rrkxz=rkx3/rkz3; rrkyz=rky3/rkz3

          End If

          !====================================================================

          If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

            ! For  monopole contribution to the stress tensor

            vterm=bb1*Exp(ralph*rksq)/rksq*electro%qqq_local_mfd(j_local,k_local,l_local)
            akv=2.0_wp*(1.0_wp/rksq-ralph)*Real( vterm*Conjg(electro%qqq_local_mfd(j_local,k_local,l_local)),wp )

            !=======================================================================
            ! For higher order contributions to the stress tensor

            pterm=2.0_wp*bb1*Exp(ralph*rksq)/rksq*Conjg(electro%qqq_local_mfd(j_local,k_local,l_local))
            sq1=Real(pterm*electro%qt1_local_mfd(j_local,k_local,l_local),wp)
            sq2=Real(pterm*electro%qt2_local_mfd(j_local,k_local,l_local),wp)
            sq3=Real(pterm*electro%qt3_local_mfd(j_local,k_local,l_local),wp)

            !=======================================================================

            strs(1)=strs(1)-rkx3*rkx3*akv+sq1
            strs(5)=strs(5)-rky3*rky3*akv+sq2
            strs(9)=strs(9)-rkz3*rkz3*akv+sq3
            strs(2)=strs(2)-rkx3*rky3*akv+rrkxy*sq2
            strs(3)=strs(3)-rkx3*rkz3*akv+rrkxz*sq3
            strs(4)=strs(4)-rkx3*rky3*akv+rrkyx*sq1
            strs(6)=strs(6)-rky3*rkz3*akv+rrkyz*sq3
            strs(7)=strs(7)-rkx3*rkz3*akv+rrkzx*sq1
            strs(8)=strs(8)-rky3*rkz3*akv+rrkzy*sq2

            electro%qqq_local_mfd(j_local,k_local,l_local)=vterm

          Else

            electro%qqq_local_mfd(j_local,k_local,l_local)=(0.0_wp,0.0_wp)

          End If

        End Do

      End Do

    End Do

    ! as only looped over local stuff, we need to gsum strs

    Call gsum(comm,strs)

    ! scale strs and distribute per node

    strs = strs * scale / Real(comm%mxnode,wp)

    ! calculate atomic energy

    Call pfft(electro%qqq_local_mfd,electro%pfft_work_mfd,electro%context_mfd,-1)

    eng = 0.0_wp
    Do l=1,electro%block_z_mfd
      Do k=1,electro%block_y_mfd
        Do j=1,electro%block_x_mfd
          qqc_tmp=Real(electro%qqq_local_mfd(j,k,l),wp)
          eng=eng+electro%qqc_local_mfd(j,k,l)*qqc_tmp
          electro%qqc_local_mfd(j,k,l)=qqc_tmp
        End Do
      End Do
    End Do

    ! as only looped over local stuff, we need to gsum the eng

    Call gsum(comm,eng)

    ! scale eng and distribute per node

    eng = eng * scale / Real(comm%mxnode,wp)

    ! Second part of the monopole contribution to the stress tensor
    ! calculate stress tensor (symmetrical, per node)

    strs      = strs      + eng*uni
    stress(1) = stress(1) + strs(1)
    stress(2) = stress(2) + 0.5_wp * (strs(2) + strs(4))
    stress(3) = stress(3) + 0.5_wp * (strs(3) + strs(7))
    stress(4) = stress(4) + 0.5_wp * (strs(2) + strs(4))
    stress(5) = stress(5) + strs(5)
    stress(6) = stress(6) + 0.5_wp * (strs(6) + strs(8))
    stress(7) = stress(7) + 0.5_wp * (strs(3) + strs(7))
    stress(8) = stress(8) + 0.5_wp * (strs(6) + strs(8))
    stress(9) = stress(9) + strs(9)

    ! distribute energy and virial terms (per node)

    engcpe_rc = eng
    vircpe_rc = -(strs(1)+strs(5)+strs(9))

    ! infrequent calculations copying

    If (ewld%l_cp) Then
      ewld%e_rc=engcpe_rc
      ewld%v_rc=vircpe_rc
      ewld%s_rc=strs
    End If

    ! calculate atomic forces

    Call spme_mforces(rcell,scale,ixx,iyy,izz,bsddx,bsddy,bsddz,electro%qqc_local_mfd,&
      electro%ixb_mfd,electro%ixt_mfd,electro%iyb_mfd,electro%iyt_mfd,electro%izb_mfd,electro%izt_mfd,mpoles)

    Deallocate (ixx,iyy,izz,it,    Stat = fail(1))
    Deallocate (bdx,bdy,bdz,       Stat = fail(2))
    Deallocate (bsddx,bsddy,bsddz, Stat = fail(3))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'ewald_spme_mforces deallocation failure'
      Call error(0,message)
    End If

  Contains

    Subroutine spme_mforces(rcell,scale,ixx,iyy,izz,bsddx,bsddy,bsddz, &
      qqc_local,ixb,ixt,iyb,iyt,izb,izt,mpoles)

      !!----------------------------------------------------------------------!
      !!
      !! dl_poly_4 subroutine for calculating coulombic forces due to
      !! multipolar interactions in a periodic system using smooth particle
      !! mesh ewald method (fourier part)
      !!
      !! Note: qqc_local is shifted from its definition from above
      !!       and therefore there is no need for periodic images (!!)
      !!
      !! copyright - daresbury laboratory
      !! author    - w.smith & i.t.todorov february 2016
      !! amended   - h.a.boateng may 2014
      !!
      !!----------------------------------------------------------------------!

      Integer,           Intent( In    ) :: ixx(1:config%mxatms),iyy(1:config%mxatms),izz(1:config%mxatms), &
        ixb,ixt, iyb,iyt, izb,izt
      Real( Kind = wp ), Intent( In    ) :: scale,rcell(1:9),                &
        bsddx(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
        bsddy(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
        bsddz(0:ewld%bspline%num_splines,1:ewld%bspline%num_splines,1:config%mxatms), &
        qqc_local( ixb:ixt, iyb:iyt, izb:izt )
      Type( mpole_type ), Intent( InOut ) :: mpoles

      Integer           :: fail(1:2), delspl, ixdb,iydb,izdb,ixdt,iydt,izdt, &
        i,j,k,l, jj,kk,ll
      Real( Kind = wp ) :: ka11,kb22,kc33,ka11sq,kb22sq,kc33sq,                          &
        ka11cu,kb22cu,kc33cu,kakb,kakc,kbkc,                          &
        kasqkb,kasqkc,kakbsq,kakcsq,kbsqkc,kbkcsq,kakbkc,             &
        tmp,fff(0:3),fix,fiy,fiz,qsum,                                &
        tix,tiy,tiz,dtp,tq2,tq3,tq4,tq5,tq6,tq7,tq8,tq9,tq10,         &
        imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10,           &
        impx1,impx2,impx3,impx4,impx5,impx6,impx7,impx8,impx9,impx10, &
        impy1,impy2,impy3,impy4,impy5,impy6,impy7,impy8,impy9,impy10, &
        impz1,impz2,impz3,impz4,impz5,impz6,impz7,impz8,impz9,impz10

      Real( Kind = wp ) :: imp(1:mpoles%max_mpoles)
      Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)

      Real( Kind = wp ), Dimension( : ),       Allocatable :: bdx,bdy,bdz
      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: qqc_domain

      ! Define extended ranges for the domain = local + halo slice and allocate

      ixdb = ixb - ewld%bspline%num_spline_padded
      iydb = iyb - ewld%bspline%num_spline_padded
      izdb = izb - ewld%bspline%num_spline_padded

      delspl = ewld%bspline%num_spline_padded - ewld%bspline%num_splines

      ixdt = ixt + delspl
      iydt = iyt + delspl
      izdt = izt + delspl

      fail=0
      Allocate ( &
        & bdx(0:ewld%bspline%num_splines), &
        & bdy(0:ewld%bspline%num_splines), &
        & bdz(0:ewld%bspline%num_splines), Stat = fail(1))
      Allocate (qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt ), Stat = fail(2))
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'spme_mforces allocation failure'
        Call error(0,message)
      End If

      Call exchange_grid( ixb , ixt , iyb , iyt , izb , izt , qqc_local, &
        ixdb, iydb, izdb, ixdt, iydt, izdt, qqc_domain, &
        domain, comm)! , ewld

      ! Real values of kmax vectors

      ka11 = Real(ewld%kspace%k_vec_dim(1),wp)*rcell(1)
      kb22 = Real(ewld%kspace%k_vec_dim(2),wp)*rcell(5)
      kc33 = Real(ewld%kspace%k_vec_dim(3),wp)*rcell(9)

      ka11sq=ka11*ka11; ka11cu=ka11sq*ka11; kb22sq=kb22*kb22; kb22cu=kb22sq*kb22
      kc33sq=kc33*kc33; kc33cu=kc33sq*kc33; kakb=ka11*kb22; kakc=ka11*kc33
      kbkc=kb22*kc33; kasqkb=ka11sq*kb22; kasqkc=ka11sq*kc33; kakbsq=ka11*kb22sq
      kakcsq=ka11*kc33sq; kbsqkc=kb22sq*kc33; kbkcsq=kb22*kc33sq; kakbkc=ka11*kb22*kc33

      tmp=-2.0_wp*scale

      fff=0.0_wp

      Do i=1,config%natms

        ! get the multipoles for site i

        imp=mpoles%global_frame(:,i)

        If (Maxval(Abs(imp)) > zero_plus) Then

          ! scale imp multipoles

          imp=tmp*imp

          ! get the components for site i infinitesimal rotations

          impx=scale*mpoles%rotation_x(:,i)
          impy=scale*mpoles%rotation_y(:,i)
          impz=scale*mpoles%rotation_z(:,i)

          imp1=imp(1); impx1=impx(1); impy1=impy(1); impz1=impz(1)

          If (mpoles%max_order >= 1) Then

            imp2=imp(2); imp3=imp(3); imp4=imp(4)

            impx2=impx(2); impx3=impx(3); impx4=impx(4)
            impy2=impy(2); impy3=impy(3); impy4=impy(4)
            impz2=impz(2); impz3=impz(3); impz4=impz(4)

          End If

          If (mpoles%max_order == 2) Then

            imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

            impx5=impx(5); impx6=impx(6); impx7=impx(7); impx8=impx(8); impx9=impx(9); impx10=impx(10)
            impy5=impy(5); impy6=impy(6); impy7=impy(7); impy8=impy(8); impy9=impy(9); impy10=impy(10)
            impz5=impz(5); impz6=impz(6); impz7=impz(7); impz8=impz(8); impz9=impz(9); impz10=impz(10)

          End If

          ! initialise forces & torques

          fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp
          tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp

          Do l=1,ewld%bspline%num_splines
            ll=izz(i)-l+2

            bdz=bsddz(:,l,i)

            Do k=1,ewld%bspline%num_splines
              kk=iyy(i)-k+2

              bdy=bsddy(:,k,i)

              Do j=1,ewld%bspline%num_splines
                jj=ixx(i)-j+2

                bdx=bsddx(:,j,i)

                qsum=qqc_domain(jj,kk,ll)

                dtp = imp1*bdx(0)*bdy(0)*bdz(0)

                tq2 = qsum*ka11*bdx(1)*bdy(0)*bdz(0)
                tq3 = qsum*kb22*bdx(0)*bdy(1)*bdz(0)
                tq4 = qsum*kc33*bdx(0)*bdy(0)*bdz(1)

                fix = fix+imp1*tq2
                fiy = fiy+imp1*tq3
                fiz = fiz+imp1*tq4

                If (mpoles%max_order >= 1) Then

                  ! torques

                  tix = tix + impx2*tq2 + impx3*tq3 + impx4*tq4
                  tiy = tiy + impy2*tq2 + impy3*tq3 + impx4*tq4
                  tiz = tiz + impz2*tq2 + impz3*tq3 + impz4*tq4

                  ! forces

                  tq5 = qsum*ka11sq*bdx(2)*bdy(0)*bdz(0); tq6 = qsum*kakb*bdx(1)*bdy(1)*bdz(0)
                  tq7 = qsum*kakc*bdx(1)*bdy(0)*bdz(1);   tq8 = qsum*kb22sq*bdx(0)*bdy(2)*bdz(0)
                  tq9 = qsum*kbkc*bdx(0)*bdy(1)*bdz(1);   tq10= qsum*kc33sq*bdx(0)*bdy(0)*bdz(2)

                  fix = fix+imp2*tq5+imp3*tq6+imp4*tq7
                  fiy = fiy+imp2*tq6+imp3*tq8+imp4*tq9
                  fiz = fiz+imp2*tq7+imp3*tq9+imp4*tq10

                End If

                If (mpoles%max_order==2) Then

                  ! torques

                  tix = tix + impx5*tq5+impx6*tq6+impx7*tq7+impx8*tq8+impx9*tq9+impx10*tq10
                  tiy = tiy + impy5*tq5+impy6*tq6+impy7*tq7+impy8*tq8+impy9*tq9+impy10*tq10
                  tiz = tiz + impz5*tq5+impz6*tq6+impz7*tq7+impz8*tq8+impz9*tq9+impz10*tq10

                  ! forces

                  fix = fix+qsum*(imp5*ka11cu*bdx(3)*bdy(0)*bdz(0)+imp6*kasqkb*bdx(2)*bdy(1)*bdz(0)+ &
                    imp7*kasqkc*bdx(2)*bdy(0)*bdz(1)+imp8*kakbsq*bdx(1)*bdy(2)*bdz(0)+ &
                    imp9*kakbkc*bdx(1)*bdy(1)*bdz(1)+imp10*kakcsq*bdx(1)*bdy(0)*bdz(2))

                  fiy = fiy+qsum*(imp5*kasqkb*bdx(2)*bdy(1)*bdz(0)+imp6*kakbsq*bdx(1)*bdy(2)*bdz(0)+ &
                    imp7*kakbkc*bdx(1)*bdy(1)*bdz(1)+imp8*kb22cu*bdx(0)*bdy(3)*bdz(0)+ &
                    imp9*kbsqkc*bdx(0)*bdy(2)*bdz(1)+imp10*kbkcsq*bdx(0)*bdy(1)*bdz(2))

                  fiz = fiz+qsum*(imp5*kasqkc*bdx(2)*bdy(0)*bdz(1)+imp6*kakbkc*bdx(1)*bdy(1)*bdz(1)+ &
                    imp7*kakcsq*bdx(1)*bdy(0)*bdz(2)+imp8*kbsqkc*bdx(0)*bdy(2)*bdz(1)+ &
                    imp9*kbkcsq*bdx(0)*bdy(1)*bdz(2)+imp10*kc33cu*bdx(0)*bdy(0)*bdz(3))

                End If

              End Do

            End Do

          End Do

          ! accumulate forces

          fff(0)=fff(0)+1.0_wp
          fff(1)=fff(1)+fix
          fff(2)=fff(2)+fiy
          fff(3)=fff(3)+fiz

          ! load forces

          config%parts(i)%fxx=config%parts(i)%fxx+fix
          config%parts(i)%fyy=config%parts(i)%fyy+fiy
          config%parts(i)%fzz=config%parts(i)%fzz+fiz

          ! and torque

          mpoles%torque_x(i)=mpoles%torque_x(i)+0.5_wp*tix
          mpoles%torque_y(i)=mpoles%torque_y(i)+0.5_wp*tiy
          mpoles%torque_z(i)=mpoles%torque_z(i)+0.5_wp*tiz

          ! infrequent calculations copying

          If (ewld%l_cp) Then
            ewld%fcx(i)=ewld%fcx(i)+fix
            ewld%fcy(i)=ewld%fcy(i)+fiy
            ewld%fcz(i)=ewld%fcz(i)+fiz
          End If

        End If

      End Do

      ! remove COM drift arising from SPME approximations

      Call gsum(comm,fff)
      If (fff(0) > zero_plus) Then
        fff(1:3)=fff(1:3)/fff(0)

        Do i=1,config%natms
          imp=mpoles%global_frame(:,i)
          If (Maxval(Abs(imp)) > zero_plus) Then

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

      Deallocate (bdx,bdy,bdz, Stat = fail(1))
      Deallocate (qqc_domain,  Stat = fail(2))
      If (Any(fail > 0)) Then
        Write(message,'(a)') 'spme_mforces dealocation failure'
        Call error(0,message)
      End If

    End Subroutine spme_mforces

  End Subroutine ewald_spme_mforces_d

  Subroutine ewald_excl_mforces(iatm,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex,stress, &
    ewld,neigh,mpoles,electro,domain,config)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using multipoles with the ewald real space
    !! kernel
    !!
    !! Note: exclusion correction term
    !!
    !! copyright - daresbury laboratory
    !! author    - h.a.boateng june 2016
    !! amended   - i.t.todorov february 2016
    !!
    !!----------------------------------------------------------------------!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( ewald_type ), Intent( In    ) :: ewld
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ),                       Intent( InOut ) :: config

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
    Real( Kind = wp ), Parameter :: rreg = 0.001_wp

    Integer           :: limit,idi,jatm,k1,k2,k3,s1,s2,s3,m,n
    Integer           :: ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,rrr,engmpl,fix,fiy,fiz,fx,fy,fz,    &
      strs1,strs2,strs3,strs5,strs6,strs9,    &
      alpr,alpr2,t1,txyz,erfr,exp1,kx,ky,kz,  &
      tt,tmp,tmpi,tmpj,tix,tiy,tiz, &
      alphan,tjx,tjy,tjz,sx,sy,sz

    Real( Kind = wp ) :: d1(-2:2*mpoles%max_order+1,-2:2*mpoles%max_order+1,-2:2*mpoles%max_order+1)
    Real( Kind = wp ) :: imp(1:mpoles%max_mpoles),jmp(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: jmpx(1:mpoles%max_mpoles),jmpy(1:mpoles%max_mpoles),jmpz(1:mpoles%max_mpoles)

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

    ! get the multipoles for site i

    imp=mpoles%global_frame(:,iatm)

    ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

      ! get the components for site i infinitesimal rotations

      impx=mpoles%rotation_x(:,iatm)
      impy=mpoles%rotation_y(:,iatm)
      impz=mpoles%rotation_z(:,iatm)

      ! multipole scaler

      scl=2.0_wp*ewld%alpha*r4pie0/(sqrpi*electro%eps)

      ! scale imp multipoles

      imp=imp*scl

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! initialize torques for atom i (temporary)

      tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

      ! Get neigh%list limit

      limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

      ! start of primary loop for forces evaluation

      Do m=1,limit

        ! atomic index

        jatm=neigh%list(neigh%list(0,iatm)+m,iatm)

        ! get the multipoles for site j

        jmp=mpoles%global_frame(:,jatm)

        ! interatomic distance

        rrr = rrt(m)

        ! truncation of potential

        If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

          ! get the components for site j infinitesimal rotations

          jmpx=mpoles%rotation_x(:,jatm)
          jmpy=mpoles%rotation_y(:,jatm)
          jmpz=mpoles%rotation_z(:,jatm)

          ! get the value of the kernel using 3pt interpolation

          alpr =rrr*ewld%alpha
          alpr2=alpr*alpr

          ! calculate error function and derivative

          If (alpr < 1.0e-2_wp) Then

            ! close particles (core-shell units) - small distances limit

            erfr=2.0_wp/sqrpi * &
              (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

            ! compute derivatives of kernel using a regularization

            If (rrr < rreg) Then
              Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr, &
                ewld%alpha*xxt(m),ewld%alpha*yyt(m),ewld%alpha*zzt(m), &
                ewld%alpha*sqrt(rrr**2+rreg**2),mpoles%max_order,d1)
            Else
              Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr, &
                ewld%alpha*xxt(m),ewld%alpha*yyt(m),ewld%alpha*zzt(m), &
                ewld%alpha*rrr,mpoles%max_order,d1)
            End If

          Else

            ! distant particles - traditional

            exp1=Exp(-(ewld%alpha*rrr)**2)
            tt  =1.0_wp/(1.0_wp+pp*ewld%alpha*rrr)

            erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(ewld%alpha*rrr)

            ! compute derivatives of kernel

            Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr, &
              ewld%alpha*xxt(m),ewld%alpha*yyt(m),ewld%alpha*zzt(m), &
              ewld%alpha*rrr,mpoles%max_order,d1)

          End If

          ! calculate forces

          engmpl = 0.0_wp
          fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
          tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

          If (mpoles%max_order < 5) Then

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  jj = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(jj)) > zero_plus) Call explicit_ewald_real_loops &
                    (-2,2*mpoles%max_order+1, k1,k2,k3, ewld%alpha, d1,               &
                    imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                    kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                    engmpl,fx,fy,fz,mpoles)

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          Else

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  jj = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(jj)) > zero_plus) Then

                    txyz=kx*jmp(jj)

                    sz = 1.0_wp
                    Do s3=0,mpoles%max_order
                      ks3=k3+s3; ks31=ks3+1

                      sy = sz
                      Do s2=0,mpoles%max_order-s3
                        ks2=k2+s2; ks21=ks2+1

                        sx = sy
                        Do s1=0,mpoles%max_order-s3-s2
                          ks1=k1+s1; ks11=ks1+1

                          n      = ks1+ks2+ks3
                          alphan = ewld%alpha**n

                          ii     = mpoles%map(s1,s2,s3)

                          tmp    = alphan*d1(ks1,ks2,ks3)

                          tmpi   = txyz       * tmp
                          tmpj   = sx*imp(ii) * tmp

                          t1     = alphan     * txyz*imp(ii)

                          ! energy

                          engmpl = engmpl   + t1*d1(ks1,ks2,ks3)

                          ! force

                          t1      = t1*ewld%alpha

                          fx      = fx      - t1*d1(ks11,ks2,ks3)
                          fy      = fy      - t1*d1(ks1,ks21,ks3)
                          fz      = fz      - t1*d1(ks1,ks2,ks31)

                          ! torque on iatm

                          tix     = tix     + impx(ii)*tmpi
                          tiy     = tiy     + impy(ii)*tmpi
                          tiz     = tiz     + impz(ii)*tmpi

                          ! torque on jatm

                          tjx     = tjx     + jmpx(jj)*tmpj
                          tjy     = tjy     + jmpy(jj)*tmpj
                          tjz     = tjz     + jmpz(jj)*tmpj

                          sx = -sx
                        End Do

                        sy = -sy
                      End Do

                      sz = -sz
                    End Do

                  End If

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          End If

          fix=fix-fx
          fiy=fiy-fy
          fiz=fiz-fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx+fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy+fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz+fz

            mpoles%torque_x(jatm)=mpoles%torque_x(jatm)-tjx
            mpoles%torque_y(jatm)=mpoles%torque_y(jatm)-tjy
            mpoles%torque_z(jatm)=mpoles%torque_z(jatm)-tjz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! accumulate potential energy

            engcpe_ex = engcpe_ex - engmpl

            ! calculate virial

            vircpe_ex = vircpe_ex + (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

            ! calculate stress tensor

            strs1 = strs1 - xxt(m)*fx
            strs2 = strs2 - xxt(m)*fy
            strs3 = strs3 - xxt(m)*fz
            strs5 = strs5 - yyt(m)*fy
            strs6 = strs6 - yyt(m)*fz
            strs9 = strs9 - zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! and torques due to multipoles

      mpoles%torque_x(iatm)=mpoles%torque_x(iatm)-scl*tix
      mpoles%torque_y(iatm)=mpoles%torque_y(iatm)-scl*tiy
      mpoles%torque_z(iatm)=mpoles%torque_z(iatm)-scl*tiz

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

  End Subroutine ewald_excl_mforces

  Subroutine ewald_excl_mforces_d(iatm,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex, &
    stress,ewld,neigh,mpoles,electro,config)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using multipoles with the ewald real space
    !! kernel
    !!
    !! Note: exclusion correction term
    !!
    !! copyright - daresbury laboratory
    !! author    - h.a.boateng june 2016
    !! amended   - i.t.todorov march 2016
    !!
    !!----------------------------------------------------------------------!

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( electrostatic_type ), Intent( InOut    ) :: electro
    Type( ewald_type ), Intent ( In    ) :: ewld
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( configuration_type ),                       Intent( InOut ) :: config

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


    Integer           :: idi,jatm,limit,m

    Real( Kind = wp ) :: scl,rsq,rrr,engmpl,fix,fiy,fiz,fx,fy,fz, &
      strs1,strs2,strs3,strs5,strs6,strs9,     &
      alpr,alpr2,tt,tmpi,tmpj,tix,talp2,       &
      tiy,tiz,tjx,tjy,tjz,b0,b1,b2,b3,b4,b5,   &
      alpsqrpi,exparr,ecc,ecd,ecq,edd,         &
      edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,     &
      zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,     &
      tm9,tm10,tm11,tm12,tm13,tm14,tm15,       &
      dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,    &
      dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,      &
      dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,           &
      dpxzzz,dpxxxx,dpyyyy,dpzzzz,             &
      dpxxxy,dpxxxz,dpxxyy,dpxxzz,dpxxyz,      &
      dpxyyy,dpyzzz,dpyyzz,dpxyyz,dpxyzz,      &
      dpyyyz,dpxxxxx,dpyyyyy,dpzzzzz,dpxxxxy,  &
      dpxxxxz,dpyyyyz,dpxzzzz,dpyzzzz,dpxyyyy, &
      dpxxxyy,dpxxxzz,dpyyyzz,dpyyzzz,dpxxyyy, &
      dpxxzzz,dpxyyzz,dpxxyyz,dpxxyzz,dpxyyyz, &
      dpxyzzz,dpxxxyz,xyz,                     &
      thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,     &
      imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8, &
      imp9,imp10,jmp1,jmp2,jmp3,jmp4,jmp5,     &
      jmp6,jmp7,jmp8,jmp9,jmp10,impx1,impx2,   &
      impx3,impx4,impx5,impx6,impx7,impx8,     &
      impx9,impx10,jmpx1,jmpx2,jmpx3,jmpx4,    &
      jmpx5,jmpx6,jmpx7,jmpx8,jmpx9,jmpx10,    &
      impy1,impy2,impy3,impy4,impy5,impy6,     &
      impy7,impy8,impy9,impy10,jmpy1,jmpy2,    &
      jmpy3,jmpy4,jmpy5,jmpy6,jmpy7,jmpy8,     &
      jmpy9,jmpy10,impz1,impz2,impz3,impz4,    &
      impz5,impz6,impz7,impz8,impz9,impz10,    &
      jmpz1,jmpz2,jmpz3,jmpz4,jmpz5,jmpz6,     &
      jmpz7,jmpz8,jmpz9,jmpz10,                &
      tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10

    Real( Kind = wp ) :: imp(1:mpoles%max_mpoles),jmp(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: jmpx(1:mpoles%max_mpoles),jmpy(1:mpoles%max_mpoles),jmpz(1:mpoles%max_mpoles)

    If (electro%newjob_emf) Then
      electro%newjob_emf = .false.

      ! coefficients for exponential in recurrence relation

      talp2    = 2.0_wp*ewld%alpha*ewld%alpha
      alpsqrpi = 1.0_wp/(ewld%alpha*sqrpi)

      electro%co1_emf = talp2*alpsqrpi
      electro%co2_emf = talp2*electro%co1_emf
      electro%co3_emf = talp2*electro%co2_emf
      electro%co4_emf = talp2*electro%co3_emf
      electro%co5_emf = talp2*electro%co4_emf
      electro%co6_emf = talp2*electro%co5_emf

      electro%alp2_emf= ewld%alpha*ewld%alpha
    End If

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

    ! get the multipoles for site i

    imp=mpoles%global_frame(:,iatm)

    ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

      ! get the components for site i infinitesimal rotations

      impx=mpoles%rotation_x(:,iatm)
      impy=mpoles%rotation_y(:,iatm)
      impz=mpoles%rotation_z(:,iatm)

      ! multipole scaler

      scl=r4pie0/electro%eps

      ! rescale multipoles

      imp=imp*scl

      imp1=imp(1); impx1=impx(1); impy1=impy(1); impz1=impz(1)

      If (mpoles%max_order >= 1) Then

        imp2=imp(2); imp3=imp(3); imp4=imp(4)

        impx2=impx(2); impx3=impx(3); impx4=impx(4)
        impy2=impy(2); impy3=impy(3); impy4=impy(4)
        impz2=impz(2); impz3=impz(3); impz4=impz(4)

      End If

      If (mpoles%max_order == 2) Then

        imp5=imp(5); imp6=imp(6); imp7=imp(7); imp8=imp(8); imp9=imp(9); imp10=imp(10)

        impx5=impx(5); impx6=impx(6); impx7=impx(7); impx8=impx(8); impx9=impx(9); impx10=impx(10)
        impy5=impy(5); impy6=impy(6); impy7=impy(7); impy8=impy(8); impy9=impy(9); impy10=impy(10)
        impz5=impz(5); impz6=impz(6); impz7=impz(7); impz8=impz(8); impz9=impz(9); impz10=impz(10)

      End If

      ! initialize torques for atom i (temporary)

      tix = 0.0_wp; tiy = 0.0_wp; tiz = 0.0_wp

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! Get neigh%list limit

      limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

      ! start of primary loop for forces evaluation

      Do m=1,limit

        ! atomic index

        jatm=neigh%list(neigh%list(0,iatm)+m,iatm)

        ! get the multipoles for site j and the components for its infinitesimal rotations

        jmp=mpoles%global_frame(:,jatm)

        jmpx=mpoles%rotation_x(:,jatm)
        jmpy=mpoles%rotation_y(:,jatm)
        jmpz=mpoles%rotation_z(:,jatm)

        jmp1=jmp(1); jmpx1=jmpx(1); jmpy1=jmpy(1); jmpz1=jmpz(1)

        If (mpoles%max_order >= 1) Then

          jmp2=jmp(2); jmp3=jmp(3); jmp4=jmp(4)

          jmpx2=jmpx(2); jmpx3=jmpx(3); jmpx4=jmpx(4)
          jmpy2=jmpy(2); jmpy3=jmpy(3); jmpy4=jmpy(4)
          jmpz2=jmpz(2); jmpz3=jmpz(3); jmpz4=jmpz(4)

        End If

        If (mpoles%max_order == 2) Then

          jmp5=jmp(5); jmp6=jmp(6); jmp7=jmp(7); jmp8=jmp(8); jmp9=jmp(9); jmp10=jmp(10)

          jmpx5=jmpx(5); jmpx6=jmpx(6); jmpx7=jmpx(7); jmpx8=jmpx(8); jmpx9=jmpx(9); jmpx10=jmpx(10)
          jmpy5=jmpy(5); jmpy6=jmpy(6); jmpy7=jmpy(7); jmpy8=jmpy(8); jmpy9=jmpy(9); jmpy10=jmpy(10)
          jmpz5=jmpz(5); jmpz6=jmpz(6); jmpz7=jmpz(7); jmpz8=jmpz(8); jmpz9=jmpz(9); jmpz10=jmpz(10)

        End If

        ! interatomic distance

        rrr = rrt(m)

        ! truncation of potential

        If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

          ! Squared distance

          rsq=rrr**2

          engmpl = 0.0_wp
          ecc = 0.0_wp; ecd = 0.0_wp; ecq = 0.0_wp
          edd = 0.0_wp; edq = 0.0_wp; eqq = 0.0_wp
          fx = 0.0_wp; fy = 0.0_wp; fz = 0.0_wp
          tjx = 0.0_wp; tjy = 0.0_wp; tjz = 0.0_wp

          xx=xxt(m); yy=yyt(m); zz=zzt(m)
          xx2=xx*xx; yy2=yy*yy; zz2=zz*zz

          ! evaluate the exponential

          exparr = Exp(-electro%alp2_emf*rsq)

          ! get the value of the kernel-erf(ewld%alpha*r)/r

          alpr =rrr*ewld%alpha
          alpr2=alpr*alpr

          ! calculate error function and derivative

          If (alpr < 1.0e-2_wp) Then

            ! close particles (core-shell units) - small distances limit

            b0 = 2.0_wp*(ewld%alpha/sqrpi) * &
              (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

            b1 = electro%co2_emf*(1.0_wp/3.0_wp-alpr2*(1.0_wp/5.0_wp-alpr2 * &
              (1.0_wp/14.0_wp-alpr2*(1.0_wp/54.0_wp-alpr2   * &
              (1.0_wp/264.0_wp-alpr2/2340.0_wp)))))

            b2 = electro%co3_emf*(1.0_wp/5.0_wp-alpr2*(1.0_wp/7.0_wp-alpr2 * &
              (1.0_wp/18.0_wp-alpr2*(1.0_wp/66.0_wp-alpr2/468.0_wp))))

            b3 = electro%co4_emf*(1.0_wp/7.0_wp-alpr2*(1.0_wp/9.0_wp-alpr2 * &
              (1.0_wp/22.0_wp-alpr2/117.0_wp)))

            b4 = electro%co5_emf*(1.0_wp/9.0_wp-alpr2*(1.0_wp/11.0_wp-alpr2/39.0_wp))

            b5 = electro%co6_emf*(1.0_wp/11.0_wp-2.0_wp*alpr2/39.0_wp)

          Else

            ! distant particles - traditional

            tt = 1.0_wp/(1.0_wp+pp*alpr)

            b0 = (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exparr)/rrr

            ! compute other recurrence terms

            b1 = (b0        - electro%co1_emf*exparr)/rsq
            b2 = (3.0_wp*b1 - electro%co2_emf*exparr)/rsq
            b3 = (5.0_wp*b2 - electro%co3_emf*exparr)/rsq
            b4 = (7.0_wp*b3 - electro%co4_emf*exparr)/rsq
            b5 = (9.0_wp*b4 - electro%co5_emf*exparr)/rsq

          End If

          ! charge-charge interaction

          ijmp= imp(1)*jmp(1)

          ecc = ijmp*b0

          bijmp = ijmp*b1
          fx  = bijmp*xx
          fy  = bijmp*yy
          fz  = bijmp*zz

          ! There is no torque for charges

          If (mpoles%max_order == 1) Then

            ! charge-dipole interactions

            tm1 = imp2*jmp1-imp1*jmp2
            tm2 = imp3*jmp1-imp1*jmp3
            tm3 = imp4*jmp1-imp1*jmp4

            tmpi= jmp1*b1
            tmpj= imp1*b1

            dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
            dpxx= xx2*b2-b1
            dpyy= yy2*b2-b1
            dpzz= zz2*b2-b1
            dpxy= xx*yy*b2
            dpxz= xx*zz*b2
            dpyz= yy*zz*b2

            ecd = tm1*dpx + tm2*dpy + tm3*dpz

            fx = fx - ( dpxx*tm1 + dpxy*tm2+dpxz*tm3 )
            fy = fy - ( dpxy*tm1 + dpyy*tm2+dpyz*tm3 )
            fz = fz - ( dpxz*tm1 + dpyz*tm2+dpzz*tm3 )

            ! dipole-dipole interactions

            dpxxx = xx*(3.0_wp*b2-xx2*b3)
            dpyyy = yy*(3.0_wp*b2-yy2*b3)
            dpzzz = zz*(3.0_wp*b2-zz2*b3)
            dpxxy = yy*(b2-xx2*b3)
            dpxxz = zz*(b2-xx2*b3)
            dpxyy = xx*(b2-yy2*b3)
            dpyyz = zz*(b2-yy2*b3)
            dpxzz = xx*(b2-zz2*b3)
            dpyzz = yy*(b2-zz2*b3)
            dpxyz = -xx*yy*zz*b3

            tm1 = imp2*jmp2
            tm2 = imp2*jmp3+imp3*jmp2
            tm3 = imp2*jmp4+imp4*jmp2
            tm4 = imp3*jmp3
            tm5 = imp3*jmp4+imp4*jmp3
            tm6 = imp4*jmp4

            edd = -(tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz)

            fx = fx + (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
            fy = fy + (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
            fz = fz + (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

            ! Torques

            tmp2 = jmp1*dpx-jmp2*dpxx-jmp3*dpxy-jmp4*dpxz
            tmp3 = jmp1*dpy-jmp2*dpxy-jmp3*dpyy-jmp4*dpyz
            tmp4 = jmp1*dpz-jmp2*dpxz-jmp3*dpyz-jmp4*dpzz

            tix  = tix + impx2*tmp2 + impx3*tmp3 + impx4*tmp4
            tiy  = tiy + impy2*tmp2 + impy3*tmp3 + impy4*tmp4
            tiz  = tiz + impz2*tmp2 + impz3*tmp3 + impz4*tmp4

            tmp2 = imp1*dpx-imp2*dpxx-imp3*dpxy-imp4*dpxz
            tmp3 = imp1*dpy-imp2*dpxy-imp3*dpyy-imp4*dpyz
            tmp4 = imp1*dpz-imp2*dpxz-imp3*dpyz-imp4*dpzz

            tjx  = tjx + jmpx2*tmp2 + jmpx3*tmp3 + jmpx4*tmp4
            tjy  = tjy + jmpy2*tmp2 + jmpy3*tmp3 + jmpy4*tmp4
            tjz  = tjz + jmpz2*tmp2 + jmpz3*tmp3 + jmpz4*tmp4

          End If

          If (mpoles%max_order == 2) Then

            ! charge-dipole interactions

            tm1 = imp2*jmp1-imp1*jmp2
            tm2 = imp3*jmp1-imp1*jmp3
            tm3 = imp4*jmp1-imp1*jmp4

            tmpi= jmp1*b1
            tmpj= imp1*b1

            dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
            dpxx= xx2*b2-b1
            dpyy= yy2*b2-b1
            dpzz= zz2*b2-b1
            dpxy= xx*yy*b2
            dpxz= xx*zz*b2
            dpyz= yy*zz*b2

            ecd = tm1*dpx + tm2*dpy + tm3*dpz

            fx = fx - ( dpxx*tm1 + dpxy*tm2 + dpxz*tm3 )
            fy = fy - ( dpxy*tm1 + dpyy*tm2 + dpyz*tm3 )
            fz = fz - ( dpxz*tm1 + dpyz*tm2 + dpzz*tm3 )

            ! dipole-dipole interactions

            dpxxx = xx*(3.0_wp*b2-xx2*b3)
            dpyyy = yy*(3.0_wp*b2-yy2*b3)
            dpzzz = zz*(3.0_wp*b2-zz2*b3)
            dpxxy = yy*(b2-xx2*b3)
            dpxxz = zz*(b2-xx2*b3)
            dpxyy = xx*(b2-yy2*b3)
            dpyyz = zz*(b2-yy2*b3)
            dpxzz = xx*(b2-zz2*b3)
            dpyzz = yy*(b2-zz2*b3)
            dpxyz = -xx*yy*zz*b3

            tm1 = imp2*jmp2
            tm2 = imp2*jmp3+imp3*jmp2
            tm3 = imp2*jmp4+imp4*jmp2
            tm4 = imp3*jmp3
            tm5 = imp3*jmp4+imp4*jmp3
            tm6 = imp4*jmp4

            edd = -(tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz)

            ! The forces for dipole-dipole interactions will be computed as part of the charge-quadrupole interactions

            ! charge-quadrupole interactions

            tm1 = imp1*jmp5+imp5*jmp1-tm1; tm2=imp1*jmp6+imp6*jmp1-tm2
            tm3 = imp1*jmp7+imp7*jmp1-tm3; tm4=imp1*jmp8+imp8*jmp1-tm4
            tm5 = imp1*jmp9+imp9*jmp1-tm5; tm6=imp1*jmp10+imp10*jmp1-tm6

            ecq = tm1*dpxx+tm2*dpxy+tm3*dpxz+tm4*dpyy+tm5*dpyz+tm6*dpzz-edd

            fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
            fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
            fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

            ! dipole-quadrupole interactions

            thrb2  = 3.0_wp*b2; thrb3 = 3.0_wp*b3

            dpxxxx = thrb2 - xx2*(6.0_wp*b3 - xx2*b4)
            dpyyyy = thrb2 - yy2*(6.0_wp*b3 - yy2*b4)
            dpzzzz = thrb2 - zz2*(6.0_wp*b3 - zz2*b4)
            dpxxxy = yy*xx*(xx2*b4-thrb3)
            dpxxxz = zz*xx*(xx2*b4-thrb3)
            dpxyyy = xx*yy*(yy2*b4-thrb3)
            dpyyyz = yy*zz*(yy2*b4-thrb3)
            dpyzzz = yy*zz*(zz2*b4-thrb3)
            dpxzzz = xx*zz*(zz2*b4-thrb3)
            dpxxyy = b2-(xx2+yy2)*b3+xx2*yy2*b4
            dpxxzz = b2-(xx2+zz2)*b3+xx2*zz2*b4
            dpyyzz = b2-(yy2+zz2)*b3+yy2*zz2*b4
            dpxxyz = yy*zz*(xx2*b4-b3)
            dpxyyz = xx*zz*(yy2*b4-b3)
            dpxyzz = xx*yy*(zz2*b4-b3)

            tm1 = imp2*jmp5-imp5*jmp2
            tm2 = imp2*jmp6-imp6*jmp2+imp3*jmp5-imp5*jmp3
            tm3 = imp2*jmp7-imp7*jmp2+imp4*jmp5-imp5*jmp4
            tm4 = imp2*jmp8-imp8*jmp2+imp3*jmp6-imp6*jmp3
            tm5 = imp2*jmp9-imp9*jmp2+imp3*jmp7-imp7*jmp3+imp4*jmp6-imp6*jmp4
            tm6 = imp2*jmp10-imp10*jmp2+imp4*jmp7-imp7*jmp4
            tm7 = imp3*jmp8-imp8*jmp3
            tm8 = imp3*jmp9-imp9*jmp3+imp4*jmp8-imp8*jmp4
            tm9 = imp3*jmp10-imp10*jmp3+imp4*jmp9-imp9*jmp4
            tm10= imp4*jmp10-imp10*jmp4

            edq = tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz + &
              tm6*dpxzz+tm7*dpyyy+tm8*dpyyz+tm9*dpyzz+tm10*dpzzz

            fx = fx - (tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz + &
              tm6*dpxxzz+tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz)

            fy = fy - (tm1*dpxxxy+tm2*dpxxyy+tm3*dpxxyz+tm4*dpxyyy+tm5*dpxyyz + &
              tm6*dpxyzz+tm7*dpyyyy+tm8*dpyyyz+tm9*dpyyzz+tm10*dpyzzz)

            fz = fz - (tm1*dpxxxz+tm2*dpxxyz+tm3*dpxxzz+tm4*dpxyyz+tm5*dpxyzz + &
              tm6*dpxzzz+tm7*dpyyyz+tm8*dpyyzz+tm9*dpyzzz+tm10*dpzzzz)

            ! quadrupole-quadrupole interactions

            fifb3 = 15.0_wp*b3; thrb4 = 3.0_wp*b4; sixb4 = 6.0_wp*b4; tenb4 = 10.0_wp*b4
            xyz   = xx*yy*zz

            dpxxxxx = xx*(xx2*(tenb4-xx2*b5)-fifb3)
            dpyyyyy = yy*(yy2*(tenb4-yy2*b5)-fifb3)
            dpzzzzz = zz*(zz2*(tenb4-zz2*b5)-fifb3)
            dpxxxxy = yy*(xx2*(sixb4-xx2*b5)-thrb3)
            dpxxxxz = zz*(xx2*(sixb4-xx2*b5)-thrb3)
            dpyyyyz = zz*(yy2*(sixb4-yy2*b5)-thrb3)
            dpyzzzz = yy*(zz2*(sixb4-zz2*b5)-thrb3)
            dpxzzzz = xx*(zz2*(sixb4-zz2*b5)-thrb3)
            dpxyyyy = xx*(yy2*(sixb4-yy2*b5)-thrb3)
            dpxxxyy = xx*((xx2+3.0_wp*yy2)*b4-xx2*yy2*b5-thrb3)
            dpxxxzz = xx*((xx2+3.0_wp*zz2)*b4-xx2*zz2*b5-thrb3)
            dpyyyzz = yy*((yy2+3.0_wp*zz2)*b4-yy2*zz2*b5-thrb3)
            dpyyzzz = zz*((3.0_wp*yy2+zz2)*b4-yy2*zz2*b5-thrb3)
            dpxxyyy = yy*((3.0_wp*xx2+yy2)*b4-xx2*yy2*b5-thrb3)
            dpxxzzz = zz*((3.0_wp*xx2+zz2)*b4-xx2*zz2*b5-thrb3)
            dpxyyzz = xx*((yy2+zz2)*b4-yy2*zz2*b5-b3)
            dpxxyyz = zz*((xx2+yy2)*b4-xx2*yy2*b5-b3)
            dpxxyzz = yy*((xx2+zz2)*b4-xx2*zz2*b5-b3)
            dpxyyyz = xyz*(thrb4-yy2*b5)
            dpxyzzz = xyz*(thrb4-zz2*b5)
            dpxxxyz = xyz*(thrb4-xx2*b5)

            tm1 = imp5*jmp5
            tm2 = imp5*jmp6+imp6*jmp5
            tm3 = imp5*jmp7+imp7*jmp5
            tm4 = imp5*jmp8+imp8*jmp5+imp6*jmp6
            tm5 = imp5*jmp9+imp9*jmp5+imp6*jmp7+imp7*jmp6
            tm6 = imp5*jmp10+imp10*jmp5+imp7*jmp7
            tm7 = imp6*jmp8+imp8*jmp6
            tm8 = imp6*jmp9+imp9*jmp6+imp7*jmp8+imp8*jmp7
            tm9 = imp6*jmp10+imp10*jmp6+imp7*jmp9+imp9*jmp7
            tm10= imp7*jmp10+imp10*jmp7
            tm11= imp8*jmp8
            tm12= imp8*jmp9+imp9*jmp8
            tm13= imp8*jmp10+imp10*jmp8+imp9*jmp9
            tm14= imp9*jmp10+imp10*jmp9
            tm15= imp10*jmp10

            eqq = tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz+tm6*dpxxzz +    &
              tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz+tm11*dpyyyy+tm12*dpyyyz + &
              tm13*dpyyzz+tm14*dpyzzz+tm15*dpzzzz

            fx = fx - (tm1*dpxxxxx+tm2*dpxxxxy+tm3*dpxxxxz+tm4*dpxxxyy+tm5*dpxxxyz+tm6*dpxxxzz +    &
              tm7*dpxxyyy+tm8*dpxxyyz+tm9*dpxxyzz+tm10*dpxxzzz+tm11*dpxyyyy+tm12*dpxyyyz + &
              tm13*dpxyyzz+tm14*dpxyzzz+tm15*dpxzzzz)

            fy = fy - (tm1*dpxxxxy+tm2*dpxxxyy+tm3*dpxxxyz+tm4*dpxxyyy+tm5*dpxxyyz+tm6*dpxxyzz +    &
              tm7*dpxyyyy+tm8*dpxyyyz+tm9*dpxyyzz+tm10*dpxyzzz+tm11*dpyyyyy+tm12*dpyyyyz + &
              tm13*dpyyyzz+tm14*dpyyzzz+tm15*dpyzzzz)

            fz = fz - (tm1*dpxxxxz+tm2*dpxxxyz+tm3*dpxxxzz+tm4*dpxxyyz+tm5*dpxxyzz+tm6*dpxxzzz +    &
              tm7*dpxyyyz+tm8*dpxyyzz+tm9*dpxyzzz+tm10*dpxzzzz+tm11*dpyyyyz+tm12*dpyyyzz + &
              tm13*dpyyzzz+tm14*dpyzzzz+tm15*dpzzzzz)

            ! second-order torques

            tmp2 = jmp1*dpx-jmp2*dpxx-jmp3*dpxy-jmp4*dpxz+jmp5*dpxxx+jmp6*dpxxy+jmp7*dpxxz+jmp8*dpxyy + &
              jmp9*dpxyz+jmp10*dpxzz
            tmp3 = jmp1*dpy-jmp2*dpxy-jmp3*dpyy-jmp4*dpyz+jmp5*dpxxy+jmp6*dpxyy+jmp7*dpxyz+jmp8*dpyyy + &
              jmp9*dpyyz+jmp10*dpyzz
            tmp4 = jmp1*dpz-jmp2*dpxz-jmp3*dpyz-jmp4*dpzz+jmp5*dpxxz+jmp6*dpxyz+jmp7*dpxzz+jmp8*dpyyz + &
              jmp9*dpyzz+jmp10*dpzzz

            tmp5 = jmp1*dpxx-jmp2*dpxxx-jmp3*dpxxy-jmp4*dpxxz+jmp5*dpxxxx+jmp6*dpxxxy+jmp7*dpxxxz + &
              jmp8*dpxxyy+jmp9*dpxxyz+jmp10*dpxxzz
            tmp6 = jmp1*dpxy-jmp2*dpxxy-jmp3*dpxyy-jmp4*dpxyz+jmp5*dpxxxy+jmp6*dpxxyy+jmp7*dpxxyz + &
              jmp8*dpxyyy+jmp9*dpxyyz+jmp10*dpxyzz
            tmp7 = jmp1*dpxz-jmp2*dpxxz-jmp3*dpxyz-jmp4*dpxzz+jmp5*dpxxxz+jmp6*dpxxyz+jmp7*dpxxzz + &
              jmp8*dpxyyz+jmp9*dpxyzz+jmp10*dpxzzz
            tmp8 = jmp1*dpyy-jmp2*dpxyy-jmp3*dpyyy-jmp4*dpyyz+jmp5*dpxxyy+jmp6*dpxyyy+jmp7*dpxyyz + &
              jmp8*dpyyyy+jmp9*dpyyyz+jmp10*dpyyzz
            tmp9 = jmp1*dpyz-jmp2*dpxyz-jmp3*dpyyz-jmp4*dpyzz+jmp5*dpxxyz+jmp6*dpxyyz+jmp7*dpxyzz + &
              jmp8*dpyyyz+jmp9*dpyyzz+jmp10*dpyzzz
            tmp10= jmp1*dpzz-jmp2*dpxzz-jmp3*dpyzz-jmp4*dpzzz+jmp5*dpxxzz+jmp6*dpxyzz+jmp7*dpxzzz + &
              jmp8*dpyyzz+jmp9*dpyzzz+jmp10*dpzzzz

            tix = tix + impx2*tmp2+impx3*tmp3+impx4*tmp4+impx5*tmp5+impx6*tmp6+impx7*tmp7+impx8*tmp8 + &
              impx9*tmp9+impx10*tmp10
            tiy = tiy + impy2*tmp2+impy3*tmp3+impy4*tmp4+impy5*tmp5+impy6*tmp6+impy7*tmp7+impy8*tmp8 + &
              impy9*tmp9+impy10*tmp10
            tiz = tiz + impz2*tmp2+impz3*tmp3+impz4*tmp4+impz5*tmp5+impz6*tmp6+impz7*tmp7+impz8*tmp8 + &
              impz9*tmp9+impz10*tmp10

            tmp2 = imp1*dpx-imp2*dpxx-imp3*dpxy-imp4*dpxz+imp5*dpxxx+imp6*dpxxy+imp7*dpxxz+imp8*dpxyy + &
              imp9*dpxyz+imp10*dpxzz
            tmp3 = imp1*dpy-imp2*dpxy-imp3*dpyy-imp4*dpyz+imp5*dpxxy+imp6*dpxyy+imp7*dpxyz+imp8*dpyyy + &
              imp9*dpyyz+imp10*dpyzz
            tmp4 = imp1*dpz-imp2*dpxz-imp3*dpyz-imp4*dpzz+imp5*dpxxz+imp6*dpxyz+imp7*dpxzz+imp8*dpyyz + &
              imp9*dpyzz+imp10*dpzzz

            tmp5 = imp1*dpxx-imp2*dpxxx-imp3*dpxxy-imp4*dpxxz+imp5*dpxxxx+imp6*dpxxxy+imp7*dpxxxz + &
              imp8*dpxxyy+imp9*dpxxyz+imp10*dpxxzz
            tmp6 = imp1*dpxy-imp2*dpxxy-imp3*dpxyy-imp4*dpxyz+imp5*dpxxxy+imp6*dpxxyy+imp7*dpxxyz + &
              imp8*dpxyyy+imp9*dpxyyz+imp10*dpxyzz
            tmp7 = imp1*dpxz-imp2*dpxxz-imp3*dpxyz-imp4*dpxzz+imp5*dpxxxz+imp6*dpxxyz+imp7*dpxxzz + &
              imp8*dpxyyz+imp9*dpxyzz+imp10*dpxzzz
            tmp8 = imp1*dpyy-imp2*dpxyy-imp3*dpyyy-imp4*dpyyz+imp5*dpxxyy+imp6*dpxyyy+imp7*dpxyyz + &
              imp8*dpyyyy+imp9*dpyyyz+imp10*dpyyzz
            tmp9 = imp1*dpyz-imp2*dpxyz-imp3*dpyyz-imp4*dpyzz+imp5*dpxxyz+imp6*dpxyyz+imp7*dpxyzz + &
              imp8*dpyyyz+imp9*dpyyzz+imp10*dpyzzz
            tmp10= imp1*dpzz-imp2*dpxzz-imp3*dpyzz-imp4*dpzzz+imp5*dpxxzz+imp6*dpxyzz+imp7*dpxzzz + &
              imp8*dpyyzz+imp9*dpyzzz+imp10*dpzzzz

            tjx = tjx + jmpx2*tmp2+jmpx3*tmp3+jmpx4*tmp4+jmpx5*tmp5+jmpx6*tmp6+jmpx7*tmp7+jmpx8*tmp8 + &
              jmpx9*tmp9+jmpx10*tmp10
            tjy = tjy + jmpy2*tmp2+jmpy3*tmp3+jmpy4*tmp4+jmpy5*tmp5+jmpy6*tmp6+jmpy7*tmp7+jmpy8*tmp8 + &
              jmpy9*tmp9+jmpy10*tmp10
            tjz = tjz + jmpz2*tmp2+jmpz3*tmp3+jmpz4*tmp4+jmpz5*tmp5+jmpz6*tmp6+jmpz7*tmp7+jmpz8*tmp8 + &
              jmpz9*tmp9+jmpz10*tmp10

          End If

          engmpl = ecc+ecd+ecq+edd+edq+eqq

          fix=fix-fx
          fiy=fiy-fy
          fiz=fiz-fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx+fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy+fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz+fz

            mpoles%torque_x(jatm)=mpoles%torque_x(jatm)-tjx
            mpoles%torque_y(jatm)=mpoles%torque_y(jatm)-tjy
            mpoles%torque_z(jatm)=mpoles%torque_z(jatm)-tjz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! accumulate potential energy

            engcpe_ex = engcpe_ex - engmpl

            ! calculate virial

            vircpe_ex = vircpe_ex + (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

            ! calculate stress tensor

            strs1 = strs1 - xxt(m)*fx
            strs2 = strs2 - xxt(m)*fy
            strs3 = strs3 - xxt(m)*fz
            strs5 = strs5 - yyt(m)*fy
            strs6 = strs6 - yyt(m)*fz
            strs9 = strs9 - zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! and torques due to multipoles

      mpoles%torque_x(iatm)=mpoles%torque_x(iatm)-scl*tix
      mpoles%torque_y(iatm)=mpoles%torque_y(iatm)-scl*tiy
      mpoles%torque_z(iatm)=mpoles%torque_z(iatm)-scl*tiz

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

  End Subroutine ewald_excl_mforces_d

  Subroutine ewald_frzn_mforces(engcpe_fr,vircpe_fr,stress,ewld,neigh,mpoles, &
    electro,domain,config,comm)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating corrections to coulombic forces
    !! in a periodic system arising from multipoles on frozen pairs
    !!
    !! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
    !!       end (and any COM drift removed) but corrections to the stress
    !!       and the virial are important as they feed into the system
    !!       pressure response.  Constant volume ensembles (ensemble < 20)
    !!       need this calculation just once! - controlled by lf_fce in
    !!       ewald_check<-two_body_forces
    !!
    !! copyright - daresbury laboratory
    !! author    - i.t.todorov & h.a.boateng february 2016
    !!
    !!----------------------------------------------------------------------!

    Real( Kind = wp ),                     Intent(   Out ) :: engcpe_fr,vircpe_fr
    Real( Kind = wp ), Dimension( 1:9 ),   Intent( InOut ) :: stress
    Type( ewald_type ),                    Intent( InOut ) :: ewld
    Type( neighbours_type ),               Intent( In    ) :: neigh
    Type( mpole_type ),                    Intent( InOut ) :: mpoles
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ),                    Intent( InOut ) :: comm
    Type( configuration_type ),                       Intent( InOut ) :: config

    Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

    Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
    Integer           :: k1,k2,k3,s1,s2,s3,n
    Integer           :: ks1,ks2,ks3,ks11,ks21,ks31,nn,mm

    Real( Kind = wp ) :: scl,det,rcell(1:9),xrr,yrr,zrr,rrr,      &
      engmpl,erfr,exp1,tt,t1,kx,ky,kz,         &
      txyz,fx,fy,fz,xss,yss,zss,               &
      strs1,strs2,strs3,strs5,strs6,strs9,tmp, &
      alphan,tmpi,tmpj,tix,tiy,tiz,tjx,tjy,tjz,sx,sy,sz

    Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr
    Real( Kind = wp ), Dimension( : ), Allocatable :: xfr,yfr,zfr
    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

    Real( Kind = wp ), Dimension( : ), Allocatable :: mmp(:,:)
    Real( Kind = wp ), Dimension( : ), Allocatable :: mmpx(:,:)
    Real( Kind = wp ), Dimension( : ), Allocatable :: mmpy(:,:)
    Real( Kind = wp ), Dimension( : ), Allocatable :: mmpz(:,:)

    Real( Kind = wp ) :: d1(-2:2*mpoles%max_order+1,-2:2*mpoles%max_order+1,-2:2*mpoles%max_order+1)
    Real( Kind = wp ) :: imp(1:mpoles%max_mpoles),jmp(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: impx(1:mpoles%max_mpoles),impy(1:mpoles%max_mpoles),impz(1:mpoles%max_mpoles)
    Real( Kind = wp ) :: jmpx(1:mpoles%max_mpoles),jmpy(1:mpoles%max_mpoles),jmpz(1:mpoles%max_mpoles)
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
      Write(message,'(a)') 'ewald_frzn_mforces allocation failure'
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

    ! initialize torques for atom i (temporary)

    tix = 0.0_wp; tiy = 0.0_wp; tiz = 0.0_wp

    l_ind=0 ; nz_fr=0
    Do i=1,config%natms

      ! If using multipoles then get multipoles in global frame for atom i

      imp=mpoles%global_frame(:,i)

      If (config%lfrzn(i) > 0 .and. Maxval(Abs(imp)) > zero_plus) Then
        nz_fr(comm%idnode+1)=nz_fr(comm%idnode+1)+1
        l_ind(nz_fr(comm%idnode+1))=i
      End If
    End Do
    Call gsum(comm, nz_fr)
    nz_fr(0) = Sum(nz_fr(0:comm%idnode)) ! Offset

    scl=2.0_wp*ewld%alpha*r4pie0/(sqrpi*electro%eps)
    nzfr = Sum(nz_fr(1:comm%mxnode))     ! Total
    If (nzfr <= 10*config%mxatms) Then

      Allocate (mmp(1:mpoles%max_mpoles,1:nzfr),                                              &
        mmpx(1:mpoles%max_mpoles,1:nzfr),mmpy(1:mpoles%max_mpoles,1:nzfr),mmpz(1:mpoles%max_mpoles,1:nzfr), &
        xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_mforces allocation failure 1'
        Call error(0,message)
      End If

      mmp=0.0_wp
      mmpx=0.0_wp ; mmpz=0.0_wp ; mmpz=0.0_wp
      xfr=0.0_wp ; yfr=0.0_wp; zfr=0.0_wp
      Do i=1,nz_fr(comm%idnode+1)
        ii=nz_fr(0)+i

        ! if using multipoles then get them in global frame for atom ii

        mmp(:,ii)=mpoles%global_frame(:,l_ind(i))

        ! get the components for site ii infinitesimal rotations

        mmpx(:,ii)=mpoles%rotation_x(:,l_ind(i))
        mmpy(:,ii)=mpoles%rotation_y(:,l_ind(i))
        mmpz(:,ii)=mpoles%rotation_z(:,l_ind(i))

        xfr(ii)=config%parts(l_ind(i))%xxx
        yfr(ii)=config%parts(l_ind(i))%yyy
        zfr(ii)=config%parts(l_ind(i))%zzz
      End Do
      Call gsum(comm,mmp)
      Call gsum(comm,mmpx)
      Call gsum(comm,mmpy)
      Call gsum(comm,mmpz)

      Call gsum(comm,xfr)
      Call gsum(comm,yfr)
      Call gsum(comm,zfr)

      Do i=1,nz_fr(comm%idnode+1)
        ii=nz_fr(0)+i

        ! get the multipoles for site ii

        imp=mmp(:,ii)*scl

        ! get the components for site ii infinitesimal rotations

        impx=mmpx(:,ii)
        impy=mmpy(:,ii)
        impz=mmpz(:,ii)

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

          rrr=Sqrt(xrr**2+yrr**2+zrr**2)

          ! get the multipoles for site jj

          jmp=mmp(:,jj)

          ! get the components for site jj infinitesimal rotations

          jmpx=mmpx(:,jj)
          jmpy=mmpy(:,jj)
          jmpz=mmpz(:,jj)

          ! calculate error function and derivative

          exp1=Exp(-(ewld%alpha*rrr)**2)
          tt  =1.0_wp/(1.0_wp+pp*ewld%alpha*rrr)

          erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(ewld%alpha*rrr)

          ! compute derivatives of kernel

          Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr,ewld%alpha*xrr,ewld%alpha*yrr, &
            ewld%alpha*zrr,ewld%alpha*rrr,mpoles%max_order,d1)

          ! calculate forces

          engmpl = 0.0_wp
          fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
          tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

          If (mpoles%max_order < 5) Then

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  nn = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(nn)) > zero_plus) Then
                    Call explicit_ewald_real_loops &
                      (-2,2*mpoles%max_order+1, k1,k2,k3, ewld%alpha, d1, &
                      imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                      kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
                      engmpl,fx,fy,fz,mpoles)
                  End If

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          Else

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  nn = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(nn)) > zero_plus) Then

                    txyz=kx*jmp(nn)

                    sz = 1.0_wp
                    Do s3=0,mpoles%max_order
                      ks3=k3+s3; ks31=ks3+1

                      sy = sz
                      Do s2=0,mpoles%max_order-s3
                        ks2=k2+s2; ks21=ks2+1

                        sx = sy
                        Do s1=0,mpoles%max_order-s3-s2
                          ks1=k1+s1; ks11=ks1+1

                          n       = ks1+ks2+ks3
                          alphan  = ewld%alpha**n

                          mm      = mpoles%map(s1,s2,s3)

                          tmp     = alphan*d1(ks1,ks2,ks3)

                          tmpi    = txyz       * tmp
                          tmpj    = sx*imp(mm) * tmp

                          t1      = alphan     * txyz*imp(mm)

                          ! energy
                          engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                          t1      = t1*ewld%alpha

                          ! force

                          fx      = fx      - t1*d1(ks11,ks2,ks3)
                          fy      = fy      - t1*d1(ks1,ks21,ks3)
                          fz      = fz      - t1*d1(ks1,ks2,ks31)

                          ! torque on iatm

                          tix     = tix     + impx(mm)*tmpi
                          tiy     = tiy     + impy(mm)*tmpi
                          tiz     = tiz     + impz(mm)*tmpi

                          ! torque on jatm

                          tjx     = tjx     + jmpx(nn)*tmpj
                          tjy     = tjy     + jmpy(nn)*tmpj
                          tjz     = tjz     + jmpz(nn)*tmpj

                          sx = -sx
                        End Do

                        sy = -sy
                      End Do

                      sz = -sz
                    End Do

                  End If

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          End If

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

          ! interatomic distance

          rrr=Sqrt(xrr**2+yrr**2+zrr**2)

          ! get the multipoles for site jj

          jmp=mmp(:,jj)

          ! get the components for site jj infinitesimal rotations

          jmpx=mmpx(:,jj)
          jmpy=mmpy(:,jj)
          jmpz=mmpz(:,jj)

          ! calculate error function and derivative

          exp1=Exp(-(ewld%alpha*rrr)**2)
          tt  =1.0_wp/(1.0_wp+pp*ewld%alpha*rrr)

          erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(ewld%alpha*rrr)

          ! compute derivatives of kernel

          Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr,ewld%alpha*xrr, &
            ewld%alpha*yrr,ewld%alpha*zrr,ewld%alpha*rrr,mpoles%max_order,d1)

          ! calculate forces

          engmpl = 0.0_wp
          fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
          tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

          If (mpoles%max_order < 5) Then

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  nn = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(nn)) > zero_plus) Call explicit_ewald_real_loops &
                    (-2,2*mpoles%max_order+1, k1,k2,k3, ewld%alpha, d1, &
                    imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                    kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
                    engmpl,fx,fy,fz,mpoles)

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          Else

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  nn = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(nn)) > zero_plus) Then

                    txyz=kx*jmp(nn)

                    sz = 1.0_wp
                    Do s3=0,mpoles%max_order
                      ks3=k3+s3; ks31=ks3+1

                      sy = sz
                      Do s2=0,mpoles%max_order-s3
                        ks2=k2+s2; ks21=ks2+1

                        sx = sy
                        Do s1=0,mpoles%max_order-s3-s2
                          ks1=k1+s1; ks11=ks1+1

                          n       = ks1+ks2+ks3
                          alphan  = ewld%alpha**n

                          mm      = mpoles%map(s1,s2,s3)

                          tmp     = alphan*d1(ks1,ks2,ks3)

                          tmpi    = txyz       * tmp
                          tmpj    = sx*imp(mm) * tmp

                          t1      = alphan     * txyz*imp(mm)

                          ! energy
                          engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                          t1      = t1*ewld%alpha

                          ! force

                          fx      = fx      - t1*d1(ks11,ks2,ks3)
                          fy      = fy      - t1*d1(ks1,ks21,ks3)
                          fz      = fz      - t1*d1(ks1,ks2,ks31)

                          ! torque on iatm

                          tix     = tix     + impx(mm)*tmpi
                          tiy     = tiy     + impy(mm)*tmpi
                          tiz     = tiz     + impz(mm)*tmpi

                          ! torque on jatm

                          tjx     = tjx     + jmpx(nn)*tmpj
                          tjy     = tjy     + jmpy(nn)*tmpj
                          tjz     = tjz     + jmpz(nn)*tmpj

                          sx = -sx
                        End Do

                        sy = -sy
                      End Do

                      sz = -sz
                    End Do

                  End If

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          End If

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-fx
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-fy
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-fz

          config%parts(l_ind(j))%fxx=config%parts(l_ind(j))%fxx+fx
          config%parts(l_ind(j))%fyy=config%parts(l_ind(j))%fyy+fy
          config%parts(l_ind(j))%fzz=config%parts(l_ind(j))%fzz+fz

          mpoles%torque_x(l_ind(j))=mpoles%torque_x(l_ind(j))+tjx
          mpoles%torque_y(l_ind(j))=mpoles%torque_y(l_ind(j))+tjy
          mpoles%torque_z(l_ind(j))=mpoles%torque_z(l_ind(j))+tjz

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

          engcpe_fr = engcpe_fr - engmpl
          vircpe_fr = vircpe_fr - (fx*xrr + fy*yrr + fz*zrr)

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

          ! interatomic distance

          rrr=Sqrt(xrr**2+yrr**2+zrr**2)

          ! get the multipoles for site jj

          jmp=mmp(:,jj)

          ! get the components for site jj infinitesimal rotations

          jmpx=mmpx(:,jj)
          jmpy=mmpy(:,jj)
          jmpz=mmpz(:,jj)

          ! calculate error function and derivative

          exp1=Exp(-(ewld%alpha*rrr)**2)
          tt  =1.0_wp/(1.0_wp+pp*ewld%alpha*rrr)

          erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(ewld%alpha*rrr)

          ! compute derivatives of kernel

          Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr,ewld%alpha*xrr, &
            ewld%alpha*yrr,ewld%alpha*zrr,ewld%alpha*rrr,mpoles%max_order,d1)

          ! calculate forces

          engmpl = 0.0_wp
          fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
          tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

          If (mpoles%max_order < 5) Then

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  nn = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(nn)) > zero_plus) Call explicit_ewald_real_loops &
                    (-2,2*mpoles%max_order+1, k1,k2,k3, ewld%alpha, d1,               &
                    imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                    kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
                    engmpl,fx,fy,fz,mpoles)

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          Else

            kz = 1.0_wp
            Do k3=0,mpoles%max_order

              ky = kz
              Do k2=0,mpoles%max_order-k3

                kx = ky
                Do k1=0,mpoles%max_order-k3-k2

                  nn = mpoles%map(k1,k2,k3)

                  If (Abs(jmp(nn)) > zero_plus) Then

                    txyz=kx*jmp(nn)

                    sz = 1.0_wp
                    Do s3=0,mpoles%max_order
                      ks3=k3+s3; ks31=ks3+1

                      sy = sz
                      Do s2=0,mpoles%max_order-s3
                        ks2=k2+s2; ks21=ks2+1

                        sx = sy
                        Do s1=0,mpoles%max_order-s3-s2
                          ks1=k1+s1; ks11=ks1+1

                          n       = ks1+ks2+ks3
                          alphan  = ewld%alpha**n

                          mm      = mpoles%map(s1,s2,s3)

                          tmp     = alphan*d1(ks1,ks2,ks3)

                          tmpi    = txyz       * tmp
                          tmpj    = sx*imp(mm) * tmp

                          t1      = alphan     * txyz*imp(mm)

                          ! energy
                          engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                          t1      = t1*ewld%alpha

                          ! force

                          fx      = fx      - t1*d1(ks11,ks2,ks3)
                          fy      = fy      - t1*d1(ks1,ks21,ks3)
                          fz      = fz      - t1*d1(ks1,ks2,ks31)

                          ! torque on iatm

                          tix     = tix     + impx(mm)*tmpi
                          tiy     = tiy     + impy(mm)*tmpi
                          tiz     = tiz     + impz(mm)*tmpi

                          ! torque on jatm

                          tjx     = tjx     + jmpx(nn)*tmpj
                          tjy     = tjy     + jmpy(nn)*tmpj
                          tjz     = tjz     + jmpz(nn)*tmpj

                          sx = -sx
                        End Do

                        sy = -sy
                      End Do

                      sz = -sz
                    End Do

                  End If

                  kx = -kx

                End Do

                ky = -ky

              End Do

              kz = -kz

            End Do

          End If

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

          engcpe_fr = engcpe_fr - engmpl
          vircpe_fr = vircpe_fr - (fx*xrr + fy*yrr + fz*zrr)

          ! calculate stress tensor

          strs1 = strs1 + xrr*fx
          strs2 = strs2 + xrr*fy
          strs3 = strs3 + xrr*fz
          strs5 = strs5 + yrr*fy
          strs6 = strs6 + yrr*fz
          strs9 = strs9 + zrr*fz
        End Do
      End Do

      Deallocate (mmp,mmpx,mmpy,mmpz,xfr,yfr,zfr, Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_mforces deallocation failure 1'
        Call error(0,message)
      End If

    Else

      ! We resort to approximating N*(N-1)/2 interactions
      ! with the short-range one from the two body linked config%cell neigh%list

      Allocate (xxt(1:neigh%max_list),yyt(1:neigh%max_list),zzt(1:neigh%max_list),rrt(1:neigh%max_list), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_mforces allocation failure 2'
        Call error(0,message)
      End If

      Do ii=1,nz_fr(comm%idnode+1)
        i=l_ind(nz_fr(comm%idnode+1))
        idi=config%ltg(ii)

        ! get the multipoles for site i

        imp=mpoles%global_frame(:,i)*scl

        ! get the components for site i infinitesimal rotations

        impx=mpoles%rotation_x(:,i)
        impy=mpoles%rotation_y(:,i)
        impz=mpoles%rotation_z(:,i)

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

          ! periodic boundary conditions
          !
          !           Call images(config%imcon,config%cell,limit,xxt,yyt,zzt)

          ! get distances

          Do k=1,limit
            rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
          End Do

          Do k=1,limit
            j=neigh%list(neigh%list(-1,i)+k,i)

            ! get the multipoles for site j

            jmp=mpoles%global_frame(:,j)

            ! get the components for site j infinitesimal rotations

            jmpx=mpoles%rotation_x(:,j)
            jmpy=mpoles%rotation_y(:,j)
            jmpz=mpoles%rotation_z(:,j)

            ! interatomic distance

            rrr=rrt(k)

            ! truncation of potential

            If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

              ! calculate error function and derivative

              exp1=Exp(-(ewld%alpha*rrr)**2)
              tt  =1.0_wp/(1.0_wp+pp*ewld%alpha*rrr)

              erfr=(1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/(ewld%alpha*rrr)

              ! compute derivatives of kernel

              Call ewald_deriv(-2,2*mpoles%max_order+1,2,erfr,ewld%alpha*xxt(k), &
                ewld%alpha*yyt(k),ewld%alpha*zzt(k),ewld%alpha*rrr,mpoles%max_order,d1)

              ! calculate forces

              engmpl = 0.0_wp
              fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
              tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

              If (mpoles%max_order < 5) Then ! fully rewritten by ITT

                kz = 1.0_wp
                Do k3=0,mpoles%max_order

                  ky = kz
                  Do k2=0,mpoles%max_order-k3

                    kx = ky
                    Do k1=0,mpoles%max_order-k3-k2

                      nn = mpoles%map(k1,k2,k3)

                      If (Abs(jmp(nn)) > zero_plus) Then
                        Call explicit_ewald_real_loops &
                          (-2,2*mpoles%max_order+1, k1,k2,k3, ewld%alpha, d1, &
                          imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                          kx*jmp(nn),jmpx(nn),jmpy(nn),jmpz(nn),tjx,tjy,tjz, &
                          engmpl,fx,fy,fz,mpoles)
                      End If

                      kx = -kx

                    End Do

                    ky = -ky

                  End Do

                  kz = -kz

                End Do

              Else

                kz = 1.0_wp
                Do k3=0,mpoles%max_order

                  ky = kz
                  Do k2=0,mpoles%max_order-k3

                    kx = ky
                    Do k1=0,mpoles%max_order-k3-k2

                      nn = mpoles%map(k1,k2,k3)

                      If (Abs(jmp(nn)) > zero_plus) Then

                        txyz=kx*jmp(nn)

                        sz = 1.0_wp
                        Do s3=0,mpoles%max_order
                          ks3=k3+s3; ks31=ks3+1

                          sy = sz
                          Do s2=0,mpoles%max_order-s3
                            ks2=k2+s2; ks21=ks2+1

                            sx = sy
                            Do s1=0,mpoles%max_order-s3-s2
                              ks1=k1+s1; ks11=ks1+1

                              n      = ks1+ks2+ks3
                              alphan = ewld%alpha**n

                              mm     = mpoles%map(s1,s2,s3)

                              tmp    = alphan*d1(ks1,ks2,ks3)

                              tmpi   = txyz       * tmp
                              tmpj   = sx*imp(mm) * tmp

                              t1     = alphan     * txyz*imp(mm)

                              ! energy

                              engmpl  = engmpl  + t1*d1(ks1,ks2,ks3)

                              ! force

                              t1      = t1*ewld%alpha

                              fx      = fx      - t1*d1(ks11,ks2,ks3)
                              fy      = fy      - t1*d1(ks1,ks21,ks3)
                              fz      = fz      - t1*d1(ks1,ks2,ks31)

                              ! torque on iatm

                              tix     = tix     + impx(mm)*tmpi
                              tiy     = tiy     + impy(mm)*tmpi
                              tiz     = tiz     + impz(mm)*tmpi

                              ! torque on jatm

                              tjx     = tjx     + jmpx(nn)*tmpj
                              tjy     = tjy     + jmpy(nn)*tmpj
                              tjz     = tjz     + jmpz(nn)*tmpj

                              sx = -sx
                            End Do

                            sy = -sy
                          End Do

                          sz = -sz
                        End Do

                      End If

                      kx = -kx

                    End Do

                    ky = -ky

                  End Do

                  kz = -kz

                End Do

              End If

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

                mpoles%torque_x(j)=mpoles%torque_x(j)+tjx
                mpoles%torque_y(j)=mpoles%torque_y(j)+tjy
                mpoles%torque_z(j)=mpoles%torque_z(j)+tjz

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

                engcpe_fr = engcpe_fr - engmpl
                vircpe_fr = vircpe_fr - (fx*xxt(k) + fy*yyt(k) + fz*zzt(k))

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

        !  update torques due to multipoles

        mpoles%torque_x(i)=mpoles%torque_x(i)+scl*tix
        mpoles%torque_y(i)=mpoles%torque_y(i)+scl*tiy
        mpoles%torque_z(i)=mpoles%torque_z(i)+scl*tiz

      End Do

      Deallocate (xxt,yyt,zzt,rrt, Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_mforces deallocation failure 2'
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
      Write(message,'(a)') 'ewald_frzn_mforces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine ewald_frzn_mforces
End Module ewald_mpole
