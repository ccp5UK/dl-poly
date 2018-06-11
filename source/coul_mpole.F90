Module coul_mpole

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring configuration variables and arrays for
! multipoles
!
! copyright - daresbury laboratory
! author    - h.a.boateng & i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,  Only : comms_type
  Use setup, Only : mxsite,mxexcl,mxspl,mxompl,mximpl,mxatdm,mxatms, &
                            r4pie0, zero_plus, mxgele, nrite, sqrpi
  Use configuration, Only : imcon,natms,ltg,fxx,fyy,fzz,xxx,yyy,zzz,cell, &
                            chge
  Use mpole
  Use mpoles_container, Only : coul_deriv, ewald_deriv, &
                               explicit_fscp_rfp_loops, explicit_ewald_real_loops, &
                               explicit_ewald_real_loops
  Use numerics,         Only : erfcgen, images_s
  Use neighbours,       Only : neighbours_type

  Implicit None

  Private

  Public :: intra_mcoul, coul_fscp_mforces, coul_rfp_mforces, coul_cp_mforces, coul_dddp_mforces, &
            coul_chrm_forces, d_ene_trq_mpoles

Contains

  Subroutine intra_mcoul(keyfce,rcut,alpha,epsq,iatm,jatm,scale, &
      rrr,xdf,ydf,zdf,coul,virele,fx,fy,fz,safe)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating bond's or 1-4 dihedral
    ! electrostatics: adjusted by a weighting factor
    !
    ! copyright - daresbury laboratory
    ! amended   - i.t.todorov & h.a.boateng february 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: keyfce,iatm,jatm
    Real( Kind = wp ), Intent( In    ) :: rcut,alpha,epsq,scale
    Real( Kind = wp ), Intent( In    ) :: xdf,ydf,zdf,rrr
    Real( Kind = wp ), Intent(   Out ) :: coul,virele,fx,fy,fz
    Logical,           Intent( InOut ) :: safe

    Integer                 :: k1,k2,k3,s1,s2,s3,n
    Integer                 :: ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: aa     = 0.0_wp , &
      bb     = 0.0_wp , &
      rfld0  = 0.0_wp , &
      rfld1  = 0.0_wp , &
      rfld2  = 0.0_wp

    Real( Kind = wp ) :: exp1,tt,erc,fer,b0,      &
      tix,tiy,tiz,tjx,tjy,tjz, &
      talpha,alphan,tmp,tmpi,tmpj,t1,t2,sx,sy,sz,kx,ky,kz,txyz

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: b1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

    Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

    If (newjob) Then
      newjob = .false.

      ! Check for damped force-shifted coulombic and reaction field interactions
      ! and set force and potential shifting parameters dependingly

      damp=.false.
      If (alpha > zero_plus) Then
        damp=.true.

        exp1= Exp(-(alpha*rcut)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rcut)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rcut
        fer = (erc + 2.0_wp*(alpha/sqrpi)*exp1)/rcut**2

        aa  = fer*rcut
        bb  = -(erc + aa*rcut)
      Else If (keyfce == 8) Then
        aa =  1.0_wp/rcut**2
        bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)
      End If

      ! set reaction field terms for RFC

      If (keyfce == 10) Then
        b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
        rfld0 = b0/rcut**3
        rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
        rfld2 = 0.5_wp*rfld0
      End If
    End If

    ! initialise defaults for coulombic energy and force contributions

    coul =0.0_wp ; virele=0.0_wp
    fx =0.0_wp ; fy =0.0_wp ; fz =0.0_wp
    tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp
    tjx=0.0_wp ; tjy=0.0_wp ; tjz=0.0_wp

    ! get the multipoles for sites i and j

    imp=mplgfr(:,iatm)
    jmp=mplgfr(:,jatm)

    ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) <= zero_plus .or. Maxval(Abs(jmp)) <= zero_plus) Return

    If(mxompl > 0 .and. induce) Then

      imp(2)=imp(2)+indipx(iatm)
      imp(3)=imp(3)+indipy(iatm)
      imp(4)=imp(4)+indipz(iatm)

      jmp(2)=jmp(2)+indipx(jatm)
      jmp(3)=jmp(3)+indipy(jatm)
      jmp(4)=jmp(4)+indipz(jatm)

    End If

    ! get the components for site i and j infinitesimal rotations

    impx=mprotx(:,iatm)
    impy=mproty(:,iatm)
    impz=mprotz(:,iatm)

    jmpx=mprotx(:,jatm)
    jmpy=mproty(:,jatm)
    jmpz=mprotz(:,jatm)

    ! default convergence factor and derivative of 'r'

    talpha = 1.0_wp
    a1     = 0.0_wp

    ! scale imp multipoles

    imp = scale*imp*r4pie0/epsq

    ! Electrostatics by ewald sum = direct coulombic

    If      (keyfce ==  2 .or. keyfce ==  6) Then

      ! compute derivatives of 1/r kernel

      Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

      ! distance dependent dielectric

    Else If (keyfce ==  4) Then

      ! Compute derivatives of 1/r^2 kernel

      Call coul_deriv(2,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

      ! force shifted coulombic and reaction field

    Else If (keyfce ==  8 .or. keyfce == 10) Then

      If (damp) Then ! calculate damping contributions

        ! compute derivatives of 'r'

        Call coul_deriv(-1,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

        ! scale the derivatives of 'r'

        a1 = aa*a1

        exp1= Exp(-(alpha*rrr)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rrr)

        fer = erc/alpha

        ! compute derivatives of the ewald real space kernel

        Call ewald_deriv(-2,2*mxompl+1,1,fer,alpha*xdf,alpha*ydf,alpha*zdf,alpha*rrr,d1)

        ! scale the derivatives into the right form

        d1 = 2.0_wp*alpha*d1/sqrpi

      End If

      If      (keyfce ==  8) Then ! force shifted coulombic
        If (.not.damp) Then ! pure

          ! compute derivatives of '1/r'

          Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

          ! compute derivatives of 'r'

          Call coul_deriv(-1,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

          ! scale the derivatives of 'r' and add to d1

          d1 = d1 + aa*a1
          a1 = 0.0_wp

        Else                ! damped

          talpha = alpha

        End If
      Else If (keyfce == 10) Then ! reaction field
        If (.not.damp) Then ! pure

          ! compute derivatives of '1/r'

          Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

          ! compute derivatives of 'r^2'

          Call coul_deriv(-2,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

          ! scale the derivatives of 'r' and add to d1

          d1 = d1 + rfld2*a1
          a1 = 0.0_wp

        Else

          ! compute derivatives of 'r^2'

          Call coul_deriv(-2,2*mxompl+1,xdf,ydf,zdf,rrr,b1)

          ! scale the derivatives of 'r^2' and add to scaled derivatives of 'r'

          a1 = a1 + rfld2*b1

          talpha = alpha

        End If
      End If

    Else

      safe = .false.

    End If

    If (safe) Then

      kz = 1.0_wp
      Do k3=0,mxompl

        ky = kz
        Do k2=0,mxompl-k3

          kx = ky
          Do k1=0,mxompl-k3-k2

            jj = mplmap(k1,k2,k3)

            If (Abs(jmp(jj)) > zero_plus) Then

              txyz=kx*jmp(jj)

              sz = 1.0_wp
              Do s3=0,mxompl
                ks3=k3+s3; ks31=ks3+1

                sy = sz
                Do s2=0,mxompl-s3
                  ks2=k2+s2; ks21=ks2+1

                  sx = sy
                  Do s1=0,mxompl-s3-s2
                    ks1=k1+s1; ks11=ks1+1

                    n      = ks1+ks2+ks3
                    alphan = talpha**n

                    ii     = mplmap(s1,s2,s3)

                    tmp    = alphan * d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                    tmpi   = txyz       * tmp
                    tmpj   = sx*imp(ii) * tmp

                    t2     = txyz*imp(ii)
                    t1     = alphan*t2

                    !  energy

                    coul    = coul + t1*d1(ks1,ks2,ks3)  + t2*a1(ks1,ks2,ks3)

                    !  force
                    t1      = t1*talpha

                    fx      = fx   - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                    fy      = fy   - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                    fz      = fz   - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

                    !  torque on iatm

                    tix     = tix  + impx(ii)*tmpi
                    tiy     = tiy  + impy(ii)*tmpi
                    tiz     = tiz  + impz(ii)*tmpi

                    !  torque on jatm

                    tjx     = tjx  + jmpx(jj)*tmpj
                    tjy     = tjy  + jmpy(jj)*tmpj
                    tjz     = tjz  + jmpz(jj)*tmpj

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

      If (keyfce == 2 .or. keyfce == 4 .or. keyfce == 6) Then

        virele = -coul

      Else If (keyfce ==  8) Then ! force shifted coulombic

        virele = -(fx*xdf + fy*ydf + fz*zdf)

        ! shift potential

        tmp    = aa*rrr + bb
        coul   = coul   + tmp*imp(1)*jmp(1)

        ! shift torque

        tmpi   = tmp*jmp(1)
        tix    = tix    + impx(1)*tmpi
        tiy    = tiy    + impy(1)*tmpi
        tiz    = tiz    + impz(1)*tmpi

        tmpj   = tmp*imp(1)
        tjx    = tjx    + jmpx(1)*tmpj
        tjy    = tjy    + jmpy(1)*tmpj
        tjz    = tjz    + jmpz(1)*tmpj

      Else If (keyfce == 10) Then ! reaction field

        virele = -(fx*xdf + fy*ydf + fz*zdf)

        ! shift potential

        coul   = coul   - rfld1*imp(1)*jmp(1)

        ! shift torque

        tmpi   = -rfld1*jmp(1)
        tix    = tix    + impx(1)*tmpi
        tiy    = tiy    + impy(1)*tmpi
        tiz    = tiz    + impz(1)*tmpi

        tmpj   = -rfld1*imp(1)
        tjx    = tjx    + jmpx(1)*tmpj
        tjy    = tjy    + jmpy(1)*tmpj
        tjz    = tjz    + jmpz(1)*tmpj

      End If

      tix = tix*r4pie0/epsq
      tiy = tiy*r4pie0/epsq
      tiz = tiz*r4pie0/epsq

      If (iatm <= natms) Then

        mptrqx(iatm)=mptrqx(iatm)+tix
        mptrqy(iatm)=mptrqy(iatm)+tiy
        mptrqz(iatm)=mptrqz(iatm)+tiz

      End If

      If (jatm <= natms) Then

        mptrqx(jatm)=mptrqx(jatm)+tjx
        mptrqy(jatm)=mptrqy(jatm)+tjy
        mptrqz(jatm)=mptrqz(jatm)+tjz

      End If

    End If
  End Subroutine intra_mcoul

  Subroutine coul_fscp_mforces &
             (iatm,alpha,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles assuming a force shifted
  ! coulomb potential kernel
  !
  ! U is proportional to ( 1/r + aa*r  + bb ) such that dU(neigh%cutoff)/dr = 0
  ! therefore aa = 1/(neigh%cutoff)**2 and U(neigh%cutoff) = 0 therefore bb = -2/(neigh%cutoff)
  !
  ! Note: FS potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: alpha,epsq
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( comms_type ),                       Intent( In    ) :: comm

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                               rdrewd = 0.0_wp , &
                               aa     = 0.0_wp , &
                               bb     = 0.0_wp

    Integer           :: fail,idi,jatm,k1,k2,k3,s1,s2,s3,m,n, &
                         k,ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,rrr,alphan,engmpl,fix,fiy,fiz,fx,fy,fz, &
                         strs1,strs2,strs3,strs5,strs6,strs9,        &
                         ppp,vk0,vk1,vk2,t1,t2,kx,ky,kz,             &
                         txyz,erfcr,tix,tiy,tiz,tjx,tjy,             &
                         tjz,tmp,tmpi,tmpj,sx,sy,sz

    Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

    Character ( Len = 256 )  ::  message

    If (newjob) Then
       newjob = .false.

       If (alpha > zero_plus) Then
          damp = .true.
       Else
          damp = .false.
       End If

       If (damp) Then

  ! interpolation interval

          drewd = neigh%cutoff/Real(mxgele-4,wp)

  ! reciprocal of interpolation interval

          rdrewd = 1.0_wp/drewd

          fail=0
          Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
          If (fail > 0) Then
             Write(message,'(a)') 'coul_fscp_mforces allocation failure'
             Call error(0,message)
          End If

  ! generate error function complement tables for ewald sum

          Call erfcgen(neigh%cutoff,alpha,mxgele,erc,fer)

  ! set force and potential shifting parameters (screened terms)

          aa =   fer(mxgele-4)*neigh%cutoff
          bb = -(erc(mxgele-4)+aa*neigh%cutoff)

       Else

  ! set force and potential shifting parameters (screened terms)

          aa =  1.0_wp/neigh%cutoff**2
          bb = -2.0_wp/neigh%cutoff ! = -(1.0_wp/neigh%cutoff+aa*neigh%cutoff)

       End If
    End If

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

  ! start of primary loop for forces evaluation

       Do m=1,neigh%list(0,iatm)

  ! atomic index

          jatm=neigh%list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! interatomic distance

          rrr = rrt(m)

  ! truncation of potential

          If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             If (damp) Then

  ! compute derivatives of 'r'

                Call coul_deriv(-1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r'

                a1 = aa*a1

  ! get the value of the ewald real space kernel using 3pt interpolation

                k   = Int(rrr*rdrewd)
                ppp = rrr*rdrewd - Real(k,wp)

                vk0 = erc(k)
                vk1 = erc(k+1)
                vk2 = erc(k+2)

                t1 = vk0 + (vk1 - vk0)*ppp
                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/alpha

  ! compute derivatives of the ewald real space kernel

                Call ewald_deriv(-2,2*mxompl+1,1,erfcr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

  ! scale the derivatives into the right form

                d1 = 2.0_wp*alpha*d1/sqrpi

  ! calculate forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then
                              Call explicit_fscp_rfp_loops &
                               (2*mxompl+1, k1,k2,k3, alpha, d1,a1,               &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n      = ks1+ks2+ks3
                                        alphan = alpha**n

                                        ii     = mplmap(s1,s2,s3)

                                        tmp    = alphan*d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                                        tmpi   = txyz       * tmp
                                        tmpj   = sx*imp(ii) * tmp

                                        t2     = txyz*imp(ii)
                                        t1     = alphan*t2

  ! energy

                                        engmpl = engmpl + t1*d1(ks1,ks2,ks3) + t2*a1(ks1,ks2,ks3)

  ! force
                                        t1     = t1*alpha

                                        fx     = fx     - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

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

             Else

  ! compute derivatives of '1/r'

                Call coul_deriv( 1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

  ! compute derivatives of 'r'

                Call coul_deriv(-1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r' and add to d1

                d1 = d1 + a1/neigh%cutoff**2

  ! calculate potential forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus)  Then
                              Call explicit_ewald_real_loops &
                               (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj=mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n    = ks1+ks2+ks3

                                        ii   = mplmap(s1,s2,s3)

                                        tmp  = d1(ks1,ks2,ks3)

                                        tmpi = txyz       * tmp
                                        tmpj = sx*imp(ii) * tmp

                                        t1   = txyz*imp(ii)

  ! energy

                                        engmpl = engmpl + t1*tmp

  ! force

                                        fx     = fx     - t1*d1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

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

             End If

  ! shift potential

             tmp     = aa*rrr + bb
             engmpl  = engmpl + tmp*imp(1)*jmp(1)

  ! shift torque

             tmpi    = tmp*jmp(1)
             tix     = tix    + impx(1)*tmpi
             tiy     = tiy    + impy(1)*tmpi
             tiz     = tiz    + impz(1)*tmpi

             tmpj    = tmp*imp(1)
             tjx     = tjx    + jmpx(1)*tmpj
             tjy     = tjy    + jmpy(1)*tmpj
             tjz     = tjz    + jmpz(1)*tmpj

             fix=fix+fx
             fiy=fiy+fy
             fiz=fiz+fz

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

  ! calculate virial

                vircpe = vircpe - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

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

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

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

  End Subroutine coul_fscp_mforces

  Subroutine coul_rfp_mforces &
             (iatm,alpha,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles assuming a reaction field
  ! potential (corrected for the existence of a dipole moment outside neigh%cutoff)
  !
  ! Note: RF potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  ! R2: M Neumann, J Chem Phys, 82 (12), 5663, (1985)
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: alpha,epsq
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( comms_type ),                       Intent( In    ) :: comm

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                               rdrewd = 0.0_wp , &
                               aa     = 0.0_wp , &
                               bb     = 0.0_wp , &
                               b0     = 0.0_wp , &
                               rfld0  = 0.0_wp , &
                               rfld1  = 0.0_wp , &
                               rfld2  = 0.0_wp

    Integer           :: fail,idi,jatm,k1,k2,k3,s1,s2,s3,m,n, &
                         k,ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,rrr,alphan,engmpl,fix,fiy,fiz,fx,fy,fz, &
                         strs1,strs2,strs3,strs5,strs6,strs9,        &
                         ppp,vk0,vk1,vk2,t1,t2,kx,ky,kz,             &
                         txyz,erfcr,tmp,tmpi,tmpj,tix,tiy,tiz,       &
                         tjx,tjy,tjz,sx,sy,sz

    Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

    Character ( Len = 256 )   ::  message

    If (newjob) Then
       newjob = .false.

       If (alpha > zero_plus) Then
          damp = .true.
       Else
          damp = .false.
       End If

  ! reaction field terms

       b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
       rfld0 = b0/neigh%cutoff**3
       rfld1 = (1.0_wp + 0.5_wp*b0)/neigh%cutoff
       rfld2 = 0.5_wp*rfld0

       If (damp) Then

  ! interpolation interval

          drewd = neigh%cutoff/Real(mxgele-4,wp)

  ! reciprocal of interpolation interval

          rdrewd = 1.0_wp/drewd

          fail=0
          Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
          If (fail > 0) Then
             Write(message,'(a)') 'coul_rfp_mforces allocation failure'
             Call error(0,message)
          End If

  ! generate error function complement tables for ewald sum

          Call erfcgen(neigh%cutoff,alpha,mxgele,erc,fer)

  ! set force and potential shifting parameters (screened terms)

          aa =   fer(mxgele-4)*neigh%cutoff
          bb = -(erc(mxgele-4)+aa*neigh%cutoff)

       Else

  ! set force and potential shifting parameters (screened terms)

          aa =  1.0_wp/neigh%cutoff**2
          bb = -2.0_wp/neigh%cutoff ! = -(1.0_wp/neigh%cutoff+aa*neigh%cutoff)

       End If
    End If

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

       Do m=1,neigh%list(0,iatm)

  ! atomic index

          jatm=neigh%list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! interatomic distance

          rrr = rrt(m)

  ! truncation of potential

          If (Maxval(Abs(jmp)) > zero_plus .and. rrr < neigh%cutoff) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             If (damp) Then

  ! compute derivatives of 'r^2'

                Call coul_deriv(-2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r^2'

                a1 = rfld2*a1

  ! compute derivatives of 'r'

                Call coul_deriv(-1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

  ! scale the derivatives of 'r' and add to a1

                a1 = a1 + aa*d1

  ! get the value of the ewald real space kernel using 3pt interpolation

                k   = Int(rrr*rdrewd)
                ppp = rrr*rdrewd - Real(k,wp)

                vk0 = erc(k)
                vk1 = erc(k+1)
                vk2 = erc(k+2)

                t1 = vk0 + (vk1 - vk0)*ppp
                t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                erfcr = (t1 + (t2-t1)*ppp*0.5_wp)/alpha

  ! compute derivatives of the ewald real space kernel

                Call ewald_deriv(-2,2*mxompl+1,1,erfcr,alpha*xxt(m),alpha*yyt(m),alpha*zzt(m),alpha*rrr,d1)

  ! scale the derivatives into the right form

                d1 = 2.0_wp*alpha*d1/sqrpi

  ! calculate forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then
                              Call explicit_fscp_rfp_loops &
                               (2*mxompl+1, k1,k2,k3, alpha, d1,a1,               &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n      = ks1+ks2+ks3
                                        alphan = alpha**n

                                        ii     = mplmap(s1,s2,s3)

                                        tmp    = alphan*d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                                        tmpi   = txyz       * tmp
                                        tmpj   = sx*imp(ii) * tmp

                                        t2     = txyz*imp(ii)
                                        t1     = alphan*t2

  ! energy

                                        engmpl = engmpl + t1*d1(ks1,ks2,ks3) + t2*a1(ks1,ks2,ks3)

  ! force
                                        t1     = t1*alpha

                                        fx     = fx     - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

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

  ! shift potential

                tmp    = bb-0.5_wp*b0/neigh%cutoff
                engmpl = engmpl + tmp*imp(1)*jmp(1)

  ! shift torque

                tmpi   = tmp*jmp(1)
                tix    = tix    + impx(1)*tmpi
                tiy    = tiy    + impy(1)*tmpi
                tiz    = tiz    + impz(1)*tmpi

                tmpj   = tmp*imp(1)
                tjx    = tjx    + jmpx(1)*tmpj
                tjy    = tjy    + jmpy(1)*tmpj
                tjz    = tjz    + jmpz(1)*tmpj

             Else

  ! compute derivatives of '1/r'

                Call coul_deriv( 1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,d1)

  ! compute derivatives of 'r^2'

                Call coul_deriv(-2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrr,a1)

  ! scale the derivatives of 'r' and add to d1

                d1 = d1 + rfld2*a1

  ! calculate potential forces

                engmpl = 0.0_wp
                fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
                tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

                If (mxompl < 5) Then

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then
                              Call explicit_ewald_real_loops &
                               (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                               imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                               kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                               engmpl,fx,fy,fz)
                           End If

                            kx = -kx

                         End Do

                         ky = -ky

                      End Do

                      kz = -kz

                   End Do

                Else

                   kz = 1.0_wp
                   Do k3=0,mxompl

                      ky = kz
                      Do k2=0,mxompl-k3

                         kx = ky
                         Do k1=0,mxompl-k3-k2

                            jj = mplmap(k1,k2,k3)

                            If (Abs(jmp(jj)) > zero_plus) Then

                               txyz=kx*jmp(jj)

                               sz = 1.0_wp
                               Do s3=0,mxompl
                                  ks3=k3+s3; ks31=ks3+1

                                  sy = sz
                                  Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                     sx = sy
                                     Do s1=0,mxompl-s3-s2
                                        ks1=k1+s1; ks11=ks1+1

                                        n    = ks1+ks2+ks3

                                        ii   = mplmap(s1,s2,s3)

                                        tmp  = d1(ks1,ks2,ks3)

                                        tmpi = txyz       * tmp
                                        tmpj = sx*imp(ii) * tmp

                                        t1   = txyz*imp(ii)

  ! energy

                                        engmpl = engmpl + t1*tmp

  ! force

                                        fx     = fx     - t1*d1(ks11,ks2,ks3)
                                        fy     = fy     - t1*d1(ks1,ks21,ks3)
                                        fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                        tix    = tix    + impx(ii)*tmpi
                                        tiy    = tiy    + impy(ii)*tmpi
                                        tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                        tjx    = tjx    + jmpx(jj)*tmpj
                                        tjy    = tjy    + jmpy(jj)*tmpj
                                        tjz    = tjz    + jmpz(jj)*tmpj

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

  ! shift potential

                engmpl = engmpl - rfld1*imp(1)*jmp(1)

  ! shift torque

                tmpi   = -rfld1*jmp(1)
                tix    = tix    + impx(1)*tmpi
                tiy    = tiy    + impy(1)*tmpi
                tiz    = tiz    + impz(1)*tmpi

                tmpj   = -rfld1*imp(1)
                tjx    = tjx    + jmpx(1)*tmpj
                tjy    = tjy    + jmpy(1)*tmpj
                tjz    = tjz    + jmpz(1)*tmpj

             End If

             fix=fix+fx
             fiy=fiy+fy
             fiz=fiz+fz

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

  ! calculate virial

                vircpe = vircpe - (fx*xxt(m) + fy*yyt(m) + fz*zzt(m))

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

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

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

  End Subroutine coul_rfp_mforces

  Subroutine coul_cp_mforces &
             (iatm,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles with 1/r kernel with no
  ! truncation or damping
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: epsq
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: idi,jatm,k1,k2,k3,s1,s2,s3,m, &
                         ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,engmpl,fix,fiy,fiz,fx,fy,fz,     &
                         strs1,strs2,strs3,strs5,strs6,strs9, &
                         t1,kx,ky,kz,txyz,tix,tiy,tiz,tjx,    &
                         tjy,tjz,tmp,tmpi,tmpj,sx,sy,sz

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

  ! start of primary loop for forces evaluation

       Do m=1,neigh%list(0,iatm)

  ! atomic index

          jatm=neigh%list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! truncation of potential - rrt(m) is the interatomic distance

          If (Maxval(Abs(jmp)) > zero_plus .and. rrt(m) < neigh%cutoff) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             Call coul_deriv(1,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrt(m),d1)

  ! calculate forces

             engmpl = 0.0_wp
             fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
             tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

             If (mxompl < 5) Then

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then
                           Call explicit_ewald_real_loops &
                             (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                             imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                             kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                             engmpl,fx,fy,fz)
                         End If

                         kx = -kx

                      End Do

                      ky = -ky

                   End Do

                   kz = -kz

                End Do

             Else

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then

                            txyz=kx*jmp(jj)

                            sz = 1.0_wp
                            Do s3=0,mxompl
                               ks3=k3+s3; ks31=ks3+1

                               sy = sz
                               Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                  sx = sy
                                  Do s1=0,mxompl-s3-s2
                                     ks1=k1+s1; ks11=ks1+1

                                     ii   = mplmap(s1,s2,s3)

                                     tmp  = d1(ks1,ks2,ks3)

                                     tmpi = txyz       * tmp
                                     tmpj = sx*imp(ii) * tmp

                                     t1   = txyz*imp(ii)

  ! energy

                                     engmpl = engmpl + t1*tmp

  ! force

                                     fx     = fx     - t1*d1(ks11,ks2,ks3)
                                     fy     = fy     - t1*d1(ks1,ks21,ks3)
                                     fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                     tix    = tix    + impx(ii)*tmpi
                                     tiy    = tiy    + impy(ii)*tmpi
                                     tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                     tjx    = tjx    + jmpx(jj)*tmpj
                                     tjy    = tjy    + jmpy(jj)*tmpj
                                     tjz    = tjz    + jmpz(jj)*tmpj

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

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

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

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

  ! virial

       vircpe = -engcpe

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

  End Subroutine coul_cp_mforces

  Subroutine coul_dddp_mforces &
             (iatm,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,neigh)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using multipoles with 1/r kernel assuming a
  ! distance dependent dielectric 'constant'
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: epsq
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: idi,jatm,k1,k2,k3,s1,s2,s3,m, &
                         ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

    Real( Kind = wp ) :: scl,engmpl,fix,fiy,fiz,fx,fy,fz,     &
                         strs1,strs2,strs3,strs5,strs6,strs9, &
                         t1,kx,ky,kz,txyz,tix,tiy,tiz,tjx,    &
                         tjy,tjz,tmp,tmpi,tmpj,sx,sy,sz

    Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
    Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
    Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  ! initialise potential energy and virial

    engcpe=0.0_wp
    vircpe=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! get the multipoles for site i

    imp=mplgfr(:,iatm)

    If (mxompl > 0 .and. induce) Then

       imp(2)=imp(2)+indipx(iatm)
       imp(3)=imp(3)+indipy(iatm)
       imp(4)=imp(4)+indipz(iatm)

    End If

  ! ignore interaction if the charge is zero

    If (Maxval(Abs(imp)) > zero_plus) Then

  ! get the components for site i infinitesimal rotations

       impx=mprotx(:,iatm)
       impy=mproty(:,iatm)
       impz=mprotz(:,iatm)

  ! multipole scaler

       scl=r4pie0/epsq

  ! scale imp multipoles

       imp=imp*scl

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! initialize torques for atom i (temporary)

       tix = 0.0_wp ; tiy = 0.0_wp ; tiz = 0.0_wp

  ! start of primary loop for forces evaluation

       Do m=1,neigh%list(0,iatm)

  ! atomic index

          jatm=neigh%list(m,iatm)

  ! get the multipoles for site j

          jmp=mplgfr(:,jatm)

          If (mxompl > 0 .and. induce) Then

             jmp(2)=jmp(2)+indipx(jatm)
             jmp(3)=jmp(3)+indipy(jatm)
             jmp(4)=jmp(4)+indipz(jatm)

          End If

  ! truncation of potential - rrt(m) is the interatomic distance

          If (Maxval(Abs(jmp)) > zero_plus .and. rrt(m) < neigh%cutoff) Then

  ! get the components for site j infinitesimal rotations

             jmpx=mprotx(:,jatm)
             jmpy=mproty(:,jatm)
             jmpz=mprotz(:,jatm)

  ! compute derivatives of kernel

             Call coul_deriv(2,2*mxompl+1,xxt(m),yyt(m),zzt(m),rrt(m),d1)

  ! calculate forces

             engmpl = 0.0_wp
             fx  = 0.0_wp ; fy  = 0.0_wp ; fz  = 0.0_wp
             tjx = 0.0_wp ; tjy = 0.0_wp ; tjz = 0.0_wp

             If (mxompl < 5) Then

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then
                           Call explicit_ewald_real_loops &
                             (-2,2*mxompl+1, k1,k2,k3, 1.0_wp, d1,              &
                             imp,       impx,    impy,    impz,    tix,tiy,tiz, &
                             kx*jmp(jj),jmpx(jj),jmpy(jj),jmpz(jj),tjx,tjy,tjz, &
                             engmpl,fx,fy,fz)
                         End If

                         kx = -kx

                      End Do

                      ky = -ky

                   End Do

                   kz = -kz

                End Do

             Else

                kz = 1.0_wp
                Do k3=0,mxompl

                   ky = kz
                   Do k2=0,mxompl-k3

                      kx = ky
                      Do k1=0,mxompl-k3-k2

                         jj = mplmap(k1,k2,k3)

                         If (Abs(jmp(jj)) > zero_plus) Then

                            txyz=kx*jmp(jj)

                            sz = 1.0_wp
                            Do s3=0,mxompl
                               ks3=k3+s3; ks31=ks3+1

                               sy = sz
                               Do s2=0,mxompl-s3
                                  ks2=k2+s2; ks21=ks2+1

                                  sx = sy
                                  Do s1=0,mxompl-s3-s2
                                     ks1=k1+s1; ks11=ks1+1

                                     ii   = mplmap(s1,s2,s3)

                                     tmp  = d1(ks1,ks2,ks3)

                                     tmpi = txyz       * tmp
                                     tmpj = sx*imp(ii) * tmp

                                     t1   = txyz*imp(ii)

  ! energy

                                     engmpl = engmpl  + t1*tmp

  ! force

                                     fx     = fx     - t1*d1(ks11,ks2,ks3)
                                     fy     = fy     - t1*d1(ks1,ks21,ks3)
                                     fz     = fz     - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

                                     tix    = tix    + impx(ii)*tmpi
                                     tiy    = tiy    + impy(ii)*tmpi
                                     tiz    = tiz    + impz(ii)*tmpi

  ! torque on jatm

                                     tjx    = tjx    + jmpx(jj)*tmpj
                                     tjy    = tjy    + jmpy(jj)*tmpj
                                     tjz    = tjz    + jmpz(jj)*tmpj

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

             If (jatm <= natms) Then

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz

                mptrqx(jatm)=mptrqx(jatm)+tjx
                mptrqy(jatm)=mptrqy(jatm)+tjy
                mptrqz(jatm)=mptrqz(jatm)+tjz

             End If

             If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! accumulate potential energy

                engcpe = engcpe + engmpl

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

       fxx(iatm)=fix
       fyy(iatm)=fiy
       fzz(iatm)=fiz

  ! and torques due to multipoles

       mptrqx(iatm)=mptrqx(iatm)+scl*tix
       mptrqy(iatm)=mptrqy(iatm)+scl*tiy
       mptrqz(iatm)=mptrqz(iatm)+scl*tiz

  ! virial

       vircpe = -engcpe

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

  End Subroutine coul_dddp_mforces

  Subroutine coul_chrm_forces(iatm,epsq,xxt,yyt,zzt,rrt,engcpe_ch,vircpe_ch,stress,neigh)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! for CHARMM model intra-core-shell interactions
  !
  ! S(r)=1-(1+u/2).exp(-u) ; u=u(r)=r_ij.(a_i+a_j)/(p_i.p_j)^(1/6)
  !
  ! S'(r)=(1+u).u'.exp(-u)/2 ; |r|'=r/|r|
  !
  ! Uchrm(r_ij) =  S(r_ij)*U(r_ij) ; F(r_ij) = -U'(r_ij)
  ! Fchrm(r_ij) = -Uchrm'(r_ij) = S(r_ij)*F(r_ij) - S'(r_ij)*U(r_ij)
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: epsq
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ch,vircpe_ch
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: limit,idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,r_r,coul,fcoul,tmp, &
                         plra,plrprd,dmpa,dmpsum,u,scr,spr,   &
                         fix,fiy,fiz,fx,fy,fz,                &
                         strs1,strs2,strs3,strs5,strs6,strs9

  ! initialise potential energy and virial

    engcpe_ch=0.0_wp
    vircpe_ch=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

  ! global identity of iatm

    idi=ltg(iatm)

  ! charge, inverse polarisability and dumping

    chgea = chge(iatm)
    plra  = plratm(iatm)
    dmpa  = dmpatm(iatm)

  ! scale main charge

    chgea = chgea*r4pie0/epsq

  ! load forces

    fix=fxx(iatm)
    fiy=fyy(iatm)
    fiz=fzz(iatm)

  ! start of primary loop for forces evaluation

  ! Get neigh%list limit

    limit=neigh%list(-3,iatm)-neigh%list(0,iatm)

    Do m=1,limit

  ! interatomic distance and derivatives

       rrr=rrt(m)
       r_r=1.0_wp/rrr

  ! atomic index, charge & inverse polarisability products
  ! and total inter-atomic summed dumping

       jatm=neigh%list(neigh%list(0,iatm)+m,iatm)
       chgprd=chgea*chge(jatm)
       plrprd=plra*plratm(jatm)
       dmpsum=dmpa+dmpatm(jatm)

       u = (dmpsum*plrprd**6) * rrr         ! dimensionless

       tmp = Exp(-u)
       scr = 1.0_wp-(1.0_wp+0.5_wp*u) * tmp ! S(r)
       spr = (1.0_wp+u) * tmp * u           ! S'(r).r

  ! calculate forces

       coul  = scr*chgprd*r_r
       tmp   = (scr-spr)*chgprd             ! used later for the virial
       fcoul = tmp*r_r**3

       fx = fcoul*xxt(m)
       fy = fcoul*yyt(m)
       fz = fcoul*zzt(m)

       fix=fix+fx
       fiy=fiy+fy
       fiz=fiz+fz

       If (jatm <= natms) Then

          fxx(jatm)=fxx(jatm)-fx
          fyy(jatm)=fyy(jatm)-fy
          fzz(jatm)=fzz(jatm)-fz

       End If

       If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! calculate potential energy and virial

          engcpe_ch = engcpe_ch + coul
          vircpe_ch = vircpe_ch - tmp*r_r

  ! calculate stress tensor

          strs1 = strs1 + xxt(m)*fx
          strs2 = strs2 + xxt(m)*fy
          strs3 = strs3 + xxt(m)*fz
          strs5 = strs5 + yyt(m)*fy
          strs6 = strs6 + yyt(m)*fz
          strs9 = strs9 + zzt(m)*fz

       End If

    End Do

  ! load back forces

    fxx(iatm)=fix
    fyy(iatm)=fiy
    fzz(iatm)=fiz

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

  End Subroutine coul_chrm_forces


  Subroutine d_ene_trq_mpoles(vircpe_dt,stress)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating the change in energy produced by
  ! an infinitesimal rotation of multipoles
  !
  ! Reference : Sagui, Pedersen, Darden, J. Chem. Phys. 120, 73 (2004)
  !             doi: 10.1063/1.1630791
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng april 2015
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Real( Kind = wp ),                   Intent(   Out ) :: vircpe_dt
    Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

    Integer           :: idi,j,iatm,jatm

    Real( Kind = wp ) :: rsq,rrr,p1(1:3),p2(1:3),u(1:3),v(1:3),w(1:3),magp1,magp2, &
                         p2u(1:3),p1perp2(1:3),p2perp1(1:3),magp1perp2,magp2perp1, &
                         tx,ty,tz,dedv,dedw,stx(1:2),sty(1:2),stz(1:2),            &
                         dux_dx,dux_dy,dux_dz,duy_dy,duy_dz,duz_dz,rmag,rmag3,     &
                         fx,fy,fz,fix,fiy,fiz,dedux,deduy,deduz,                   &
                         tmptx,tmpty,tmptz,xdf,ydf,zdf,                            &
                         strs1,strs2,strs3,strs4,strs5,strs6,strs7,strs8,strs9

  ! initialise virial

    vircpe_dt=0.0_wp

  ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs4=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs7=0.0_wp
    strs8=0.0_wp
    strs9=0.0_wp

    Do iatm = 1, natms

  ! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

  ! global identity of iatm

       idi=ltg(iatm)

       If (mprotm(iatm)%flag == 1) Then

  ! p1 and p2 define the local frame

          p1 = mprotm(iatm)%p1
          p2 = mprotm(iatm)%p2

          magp1 = Sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
          magp2 = Sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))

          p2u = p2/magp2

  ! standard basis for coordinate system

          u(1)  = mprotm(iatm)%mtrxa(1)
          u(2)  = mprotm(iatm)%mtrxa(4)
          u(3)  = mprotm(iatm)%mtrxa(7)

          v(1)  = mprotm(iatm)%mtrxa(2)
          v(2)  = mprotm(iatm)%mtrxa(5)
          v(3)  = mprotm(iatm)%mtrxa(8)

          w(1)  = mprotm(iatm)%mtrxa(3)
          w(2)  = mprotm(iatm)%mtrxa(6)
          w(3)  = mprotm(iatm)%mtrxa(9)

  ! change in energy (E) due to infinitesimal rotation (torque => \tau) dE_{\omega} = -\tau * d\omega
  ! we omit the negative here and introduce it in the final force computation, i.e., because force is
  ! the negative of the change in energy with respect to energy, we'll add the magnitude of the force
  ! computed instead of subtracting the magnitude

          tmptx = mptrqx(iatm)
          tmpty = mptrqy(iatm)
          tmptz = mptrqz(iatm)

          tx = tmptx*u(1) + tmpty*u(2) + tmptz*u(3)
          ty = (tmptx*p2(1) + tmpty*p2(2) + tmptz*p2(3)) / magp2
          tz = tmptx*w(1) + tmpty*w(2) + tmptz*w(3)

  ! component of p2 perpendicular to p1

          p2perp1    = p2 - (u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3))*u
          magp2perp1 = Sqrt(p2perp1(1)*p2perp1(1)+p2perp1(2)*p2perp1(2)+p2perp1(3)*p2perp1(3))

  ! component of p1 perpendicular to p2

          p1perp2    = p1 - (p2u(1)*p1(1) + p2u(2)*p1(2) + p2u(3)*p1(3))*p2u
          magp1perp2 = Sqrt(p1perp2(1)*p1perp2(1)+p1perp2(2)*p1perp2(2)+p1perp2(3)*p1perp2(3))

  ! For p1

  ! dedu = 0.0 since movement du along the u axis does not rotate the frame

          dedv   = tz/magp1
          dedw   =-ty/magp1perp2

          stx(1) = dedv*v(1) + dedw*w(1)
          sty(1) = dedv*v(2) + dedw*w(2)
          stz(1) = dedv*v(3) + dedw*w(3)

  ! for p2

          dedw   = tx/magp2perp1

          stx(2) = dedw*w(1)
          sty(2) = dedw*w(2)
          stz(2) = dedw*w(3)

  ! now compute forces and virial

          fx = 0.0_wp ; fy = 0.0_wp ; fz = 0.0_wp

          Do j = 1, 2

             jatm = mprotm(iatm)%mbnd(j)

             If (jatm > 0) Then
                xdf = xxx(jatm) - xxx(iatm)
                ydf = yyy(jatm) - yyy(iatm)
                zdf = zzz(jatm) - zzz(iatm)

                Call images_s(imcon,cell,xdf,ydf,zdf)

                rsq   = xdf*xdf + ydf*ydf + zdf*zdf
                rrr   = sqrt(rsq)
                rmag  = 1.0_wp/rrr
                rmag3 = 1.0_wp/(rsq*rrr)

                dux_dx = rmag - xdf * xdf * rmag3
                dux_dy =      - xdf * ydf * rmag3       ! duy_dx = dux_dy
                dux_dz =      - xdf * zdf * rmag3       ! duz_dx = dux_dz
                duy_dy = rmag - ydf * ydf * rmag3
                duy_dz =      - ydf * zdf * rmag3       ! duz_dy = duy_dz
                duz_dz = rmag - zdf * zdf * rmag3

                dedux  = stx(j)
                deduy  = sty(j)
                deduz  = stz(j)

  ! Now to find the forces (derivatives of energy with respect to cartesian positions),
  ! i.e. fx=dedx, fy=dedy, fz=dedz

                fx = dedux * dux_dx + deduy * dux_dy + deduz * dux_dz
                fy = dedux * dux_dy + deduy * duy_dy + deduz * duy_dz
                fz = dedux * dux_dz + deduy * duy_dz + deduz * duz_dz

                fix = fix + fx
                fiy = fiy + fy
                fiz = fiz + fz

                If (jatm <= natms) Then

                   fxx(jatm)=fxx(jatm)-fx
                   fyy(jatm)=fyy(jatm)-fy
                   fzz(jatm)=fzz(jatm)-fz

                End If

                If (jatm <= natms .or. idi < ltg(jatm)) Then

  ! calculate virial

                   vircpe_dt = vircpe_dt - (fx*xdf + fy*ydf + fz*zdf)

  ! calculate stress tensor

                   strs1 = strs1 + xdf*fx
                   strs2 = strs2 + xdf*fy
                   strs3 = strs3 + xdf*fz
                   strs4 = strs4 + ydf*fx
                   strs5 = strs5 + ydf*fy
                   strs6 = strs6 + ydf*fz
                   strs7 = strs7 + zdf*fx
                   strs8 = strs8 + zdf*fy
                   strs9 = strs9 + zdf*fz

                End If

             End If

          End Do

  ! load back forces

          fxx(iatm)=fix
          fyy(iatm)=fiy
          fzz(iatm)=fiz

       End If

  ! complete stress tensor (and symmetrise)

       stress(1) = stress(1) + strs1
       stress(2) = stress(2) + 0.5_wp * (strs2 + strs4)
       stress(3) = stress(3) + 0.5_wp * (strs3 + strs7)
       stress(4) = stress(4) + 0.5_wp * (strs2 + strs4)
       stress(5) = stress(5) + strs5
       stress(6) = stress(6) + 0.5_wp * (strs6 + strs8)
       stress(7) = stress(7) + 0.5_wp * (strs3 + strs7)
       stress(8) = stress(8) + 0.5_wp * (strs6 + strs8)
       stress(9) = stress(9) + strs9

    End Do

  End Subroutine d_ene_trq_mpoles

End Module coul_mpole
