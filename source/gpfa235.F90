Module gpfa235

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module containing GPFA routines and file also containing
  ! GPFA_WRAP  routine as required by parallel_fft (DaFT)
  !
  ! copyright - daresbury laboratory
  ! author    - i.j.bush august 2010
  ! amended   - i.t.todorov september 2010
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use errors_warnings, Only : error
  Implicit None

  Private

  Interface gpfa_set
    Module Procedure setgpfa
  End Interface
  Public :: gpfa, gpfa_set

Contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! I.T.Todorov march 2004 - F77_2_F90 transformation via SPAG (plusFORT)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine GPFA(A,B,Trigs,Inc,Jump,N,Lot,Isign)

    ! ----------------------------------------------------------------------
    !       self-sorting in-place generalized prime factor (complex) fft
    !
    !       **! this is the all-FORTRAN version ***
    !           -------------------------------
    !
    !       Call gpfa(a,b,trigs,inc,jump,n,lot,isign)
    !
    !       a is first real input/output vector
    !       b is first imaginary input/output vector
    !       trigs is a table of twiddle factors, precalculated
    !             by calling subroutine 'setgpfa'
    !       inc is the increment within each data vector
    !       jump is the increment between data vectors
    !       n is the length of the transforms:
    !         -----------------------------------
    !           n = (2**li) ! (3**iq) * (5**ir)
    !         -----------------------------------
    !       lot is the number of transforms
    !       isign = +1 for forward transform
    !             = -1 for inverse transform
    !
    !       written by clive temperton
    !       recherche en prevision numerique
    !       atmospheric environment service, canada
    !
    ! ----------------------------------------------------------------------
    !
    !       definition of transform
    !       -----------------------
    !
    !       x(j) = sum(k=0,...,n-1)(c(k)*Exp(isign*2*i*j*k*pi/n))
    !
    ! _---------------------------------------------------------------------
    !
    !       for a mathematical development of the algorithm used,
    !       see:
    !
    !       c temperton : "a generalized prime factor fft algorithm
    !         for any n = (2**p)(3**q)(5**r)",
    !         siam j. sci. stat. comp., may 1992.
    ! ----------------------------------------------------------------------

    !*** Start of declarations inserted by SPAG (plusFORT)
    Real( Kind = wp ) :: A(*),B(*),Trigs(*)

    Integer           :: i,ifac,Inc,li,iq,ir,Isign,Jump,kk,ll,Lot,N,nj(3),nn
    !*** End of declarations inserted by SPAG (plusFORT)

    ! decompose n into factors 2,3,5
    nn = N
    ifac = 2
    Do ll = 1,3
      kk = 0
      50      Continue
      If ( Mod(nn,ifac) /= 0 ) Then
        nj(ll) = kk
        ifac = ifac + ll
      Else
        kk = kk + 1
        nn = nn/ifac
        Go To 50
      End If
    End Do

    If ( nn /= 1 ) Then
      Write(Unit=*, Fmt=100) N
      100     Format(' *** WARNING!!!',i10,' IS NOT A LEGAL VALUE OF N *** !!!')
      Return
    End If

    li = nj(1)
    iq = nj(2)
    ir = nj(3)

    ! compute the transform
    ! ---------------------

    i = 1
    If ( li > 0 ) Then
      Call GPFA2F(A,B,Trigs,Inc,Jump,N,li,Lot,Isign)
      i = i + 2*(2**li)
    End If
    If ( iq > 0 ) Then
      Call GPFA3F(A,B,Trigs(i),Inc,Jump,N,iq,Lot,Isign)
      i = i + 2*(3**iq)
    End If
    If ( ir > 0 ) Call GPFA5F(A,B,Trigs(i),Inc,Jump,N,ir,Lot,Isign)

  End Subroutine GPFA



  Subroutine GPFA2F(A,B,Trigs,Inc,Jump,N,Mm,Lot,Isign)

    ! ----------------------------------------------------------------------
    !
    ! FORTRAN version of *gpfa2* -
    ! radix-2 section of self-sorting, in-place, generalized pfa
    ! central radix-2 and radix-8 passes included
    !     so that transform length can be any power of 2
    !
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    ! *                                                                *
    ! *     n.b. lvr = length of vector registers, set to 128 for c90. *
    ! *     reset to 64 for other cray machines, or to any large value *
    ! *     (greater than or equal to lot) for a scalar computer.      *
    ! *                                                                *
    ! ----------------------------------------------------------------------

    !*** Start of declarations inserted by SPAG (plusFORT)
    Real( Kind = wp ) :: A(*),B(*),Trigs(*),                      &
      aja,ajb,ajc,ajd,aje,ajf,ajg,ajh,aji,     &
      ajj,ajk,ajl,ajm,ajn,ajo,ajp,bja,bjb,     &
      bjc,bjd,bje,bjf,bjg,bjh,bji,bjj,bjk,bjl, &
      bjm,bjn,bjo,bjp,c1,c2,c3,co1,co2,co3,    &
      co4,co5,co6,co7,s,si1,si2,si3,si4,si5,   &
      si6,si7,ss,t0,t1,t2,t3,u0,u1,u2,u3

    Integer           :: Inc,ink,inq,ipass,Isign,istart,j,ja,jb,  &
      jc,jd,je,jf,jg,jh,ji,jj,jjj,jk,jl,       &
      jm,jn,jo,jp,jstep,jstepl,jstepx,Jump,k,  &
      kk,l,la,laincl,left,ll,Lot,lvr,m,m2,m8,  &
      mh,Mm,mu,N,n2,nb,nblox,ninc,nu,nvex

    Integer, Parameter :: ncache = 65536
    !*** End of declarations inserted by SPAG (plusFORT)

    ! set lvr dynamically - ncache = # complex elements that fit in cache *!

    If ( ncache < N ) Then
      Write(Unit=*, Fmt=*) ' ***** gpfa ***** ncache too small '
      Call ABORT_IT()
    End If
    lvr = ncache/N

    n2 = 2**Mm
    inq = N/n2
    jstepx = (n2-N)*Inc
    ninc = N*Inc
    ink = Inc*inq

    m2 = 0
    m8 = 0
    If ( Mod(Mm,2) == 0 ) Then
      m = Mm/2
    Else If ( Mod(Mm,4) == 1 ) Then
      m = (Mm-1)/2
      m2 = 1
    Else If ( Mod(Mm,4) == 3 ) Then
      m = (Mm-3)/2
      m8 = 1
    End If
    mh = (m+1)/2

    nblox = 1 + (Lot-1)/lvr
    left = Lot
    s = Real(Isign,wp)
    istart = 1

    !     loop on blocks of lvr transforms
    !     --------------------------------
    Do nb = 1,nblox

      If ( left <= lvr ) Then
        nvex = left
      Else If ( left < (2*lvr) ) Then
        nvex = left/2
        nvex = nvex + Mod(nvex,2)
      Else
        nvex = lvr
      End If
      left = left - nvex

      la = 1

      !     loop on type I radix-4 passes
      !     -----------------------------
      mu = Mod(inq,4)
      If ( Isign == -1 ) mu = 4 - mu
      ss = 1.0_wp
      If ( mu == 3 ) ss = -1.0_wp

      If ( mh /= 0 ) Then

        Do ipass = 1,mh
          jstep = (N*Inc)/(4*la)
          jstepl = jstep - ninc

          !     k = 0 loop (no twiddle factors)
          !     -------------------------------
          Do jjj = 0,(N-1)*Inc,4*jstep
            ja = istart + jjj

            ! "transverse" loop
            ! -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              jd = jc + jstepl
              If ( jd < istart ) jd = jd + ninc
              j = 0

              !     loop across transforms
              !     ----------------------
              !cdir$ ivdep, shortloop
              Do l = 1,nvex
                aja = A(ja+j)
                ajc = A(jc+j)
                t0 = aja + ajc
                t2 = aja - ajc
                ajb = A(jb+j)
                ajd = A(jd+j)
                t1 = ajb + ajd
                t3 = ss*(ajb-ajd)
                bja = B(ja+j)
                bjc = B(jc+j)
                u0 = bja + bjc
                u2 = bja - bjc
                bjb = B(jb+j)
                bjd = B(jd+j)
                u1 = bjb + bjd
                u3 = ss*(bjb-bjd)
                A(ja+j) = t0 + t1
                A(jc+j) = t0 - t1
                B(ja+j) = u0 + u1
                B(jc+j) = u0 - u1
                A(jb+j) = t2 - u3
                A(jd+j) = t2 + u3
                B(jb+j) = u2 + t3
                B(jd+j) = u2 - t3
                j = j + Jump
              End Do
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do

          !     finished if n2 = 4
          !     ------------------
          If ( n2 == 4 ) Go To 50
          kk = 2*la

          !     loop on non-zero k
          !     -----------------
          Do k = ink,jstep - ink,ink
            co1 = Trigs(kk+1)
            si1 = s*Trigs(kk+2)
            co2 = Trigs(2*kk+1)
            si2 = s*Trigs(2*kk+2)
            co3 = Trigs(3*kk+1)
            si3 = s*Trigs(3*kk+2)

            !     loop along transform
            !     --------------------
            Do jjj = k,(N-1)*Inc,4*jstep
              ja = istart + jjj

              ! "transverse" loop
              ! -----------------
              Do nu = 1,inq
                jb = ja + jstepl
                If ( jb < istart ) jb = jb + ninc
                jc = jb + jstepl
                If ( jc < istart ) jc = jc + ninc
                jd = jc + jstepl
                If ( jd < istart ) jd = jd + ninc
                j = 0

                !     loop across transforms
                !     ----------------------
                !cdir$ ivdep,shortloop
                Do l = 1,nvex
                  aja = A(ja+j)
                  ajc = A(jc+j)
                  t0 = aja + ajc
                  t2 = aja - ajc
                  ajb = A(jb+j)
                  ajd = A(jd+j)
                  t1 = ajb + ajd
                  t3 = ss*(ajb-ajd)
                  bja = B(ja+j)
                  bjc = B(jc+j)
                  u0 = bja + bjc
                  u2 = bja - bjc
                  bjb = B(jb+j)
                  bjd = B(jd+j)
                  u1 = bjb + bjd
                  u3 = ss*(bjb-bjd)
                  A(ja+j) = t0 + t1
                  B(ja+j) = u0 + u1
                  A(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                  A(jc+j) = co2*(t0-t1) - si2*(u0-u1)
                  B(jc+j) = si2*(t0-t1) + co2*(u0-u1)
                  A(jd+j) = co3*(t2+u3) - si3*(u2-t3)
                  B(jd+j) = si3*(t2+u3) + co3*(u2-t3)
                  j = j + Jump
                End Do
                ! -----( end of loop across transforms )
                ja = ja + jstepx
                If ( ja < istart ) ja = ja + ninc
              End Do
            End Do
            ! -----( end of loop along transforms )
            kk = kk + 2*la
          End Do
          ! -----( end of loop on non-zero k )
          la = 4*la
        End Do
      End If
      ! -----( end of loop on type I radix-4 passes)

      !     central radix-2 pass
      !     --------------------
      If ( m2 /= 0 ) Then

        jstep = (N*Inc)/(2*la)
        jstepl = jstep - ninc

        !     k=0 loop (no twiddle factors)
        !     -----------------------------
        Do jjj = 0,(N-1)*Inc,2*jstep
          ja = istart + jjj

          ! "transverse" loop
          ! -----------------
          Do nu = 1,inq
            jb = ja + jstepl
            If ( jb < istart ) jb = jb + ninc
            j = 0

            !     loop across transforms
            !     ----------------------
            !cdir$ ivdep, shortloop
            Do l = 1,nvex
              aja = A(ja+j)
              ajb = A(jb+j)
              t0 = aja - ajb
              A(ja+j) = aja + ajb
              A(jb+j) = t0
              bja = B(ja+j)
              bjb = B(jb+j)
              u0 = bja - bjb
              B(ja+j) = bja + bjb
              B(jb+j) = u0
              j = j + Jump
            End Do
            ! -----(end of loop across transforms)
            ja = ja + jstepx
            If ( ja < istart ) ja = ja + ninc
          End Do
        End Do

        !     finished if n2=2
        !     ----------------
        If ( n2 == 2 ) Go To 50

        kk = 2*la

        !     loop on non-zero k
        !     -----------------
        Do k = ink,jstep - ink,ink
          co1 = Trigs(kk+1)
          si1 = s*Trigs(kk+2)

          !     loop along transforms
          !     ---------------------
          Do jjj = k,(N-1)*Inc,2*jstep
            ja = istart + jjj

            ! "transverse" loop
            ! -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              j = 0

              !     loop across transforms
              !     ----------------------
              If ( kk == n2/2 ) Then
                !cdir$ ivdep, shortloop
                Do l = 1,nvex
                  aja = A(ja+j)
                  ajb = A(jb+j)
                  t0 = ss*(aja-ajb)
                  A(ja+j) = aja + ajb
                  bjb = B(jb+j)
                  bja = B(ja+j)
                  A(jb+j) = ss*(bjb-bja)
                  B(ja+j) = bja + bjb
                  B(jb+j) = t0
                  j = j + Jump
                End Do

              Else

                !cdir$ ivdep, shortloop
                Do l = 1,nvex
                  aja = A(ja+j)
                  ajb = A(jb+j)
                  t0 = aja - ajb
                  A(ja+j) = aja + ajb
                  bja = B(ja+j)
                  bjb = B(jb+j)
                  u0 = bja - bjb
                  B(ja+j) = bja + bjb
                  A(jb+j) = co1*t0 - si1*u0
                  B(jb+j) = si1*t0 + co1*u0
                  j = j + Jump
                End Do

              End If

              ! -----(end of loop across transforms)
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
          ! -----(end of loop along transforms)
          kk = kk + 2*la
        End Do
        ! -----(end of loop on non-zero k)
        ! -----(end of radix-2 pass)

        la = 2*la

        !     central radix-8 pass
        !     --------------------
      Else If ( m8 /= 0 ) Then
        jstep = (N*Inc)/(8*la)
        jstepl = jstep - ninc
        mu = Mod(inq,8)
        If ( Isign == -1 ) mu = 8 - mu
        c1 = 1.0_wp
        If ( mu == 3 .or. mu.EQ.7 ) c1 = -1.0_wp
        c2 = SQRT(0.5_wp)
        If ( mu == 3 .or. mu.EQ.5 ) c2 = -c2
        c3 = c1*c2

        !     stage 1
        !     -------
        Do k = 0,jstep - ink,ink
          Do jjj = k,(N-1)*Inc,8*jstep
            ja = istart + jjj

            ! "transverse" loop
            ! -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              jd = jc + jstepl
              If ( jd < istart ) jd = jd + ninc
              je = jd + jstepl
              If ( je < istart ) je = je + ninc
              jf = je + jstepl
              If ( jf < istart ) jf = jf + ninc
              jg = jf + jstepl
              If ( jg < istart ) jg = jg + ninc
              jh = jg + jstepl
              If ( jh < istart ) jh = jh + ninc
              j = 0
              !cdir$ ivdep, shortloop
              Do l = 1,nvex
                aja = A(ja+j)
                aje = A(je+j)
                t0 = aja - aje
                A(ja+j) = aja + aje
                ajc = A(jc+j)
                ajg = A(jg+j)
                t1 = c1*(ajc-ajg)
                A(je+j) = ajc + ajg
                ajb = A(jb+j)
                ajf = A(jf+j)
                t2 = ajb - ajf
                A(jc+j) = ajb + ajf
                ajd = A(jd+j)
                ajh = A(jh+j)
                t3 = ajd - ajh
                A(jg+j) = ajd + ajh
                A(jb+j) = t0
                A(jf+j) = t1
                A(jd+j) = c2*(t2-t3)
                A(jh+j) = c3*(t2+t3)
                bja = B(ja+j)
                bje = B(je+j)
                u0 = bja - bje
                B(ja+j) = bja + bje
                bjc = B(jc+j)
                bjg = B(jg+j)
                u1 = c1*(bjc-bjg)
                B(je+j) = bjc + bjg
                bjb = B(jb+j)
                bjf = B(jf+j)
                u2 = bjb - bjf
                B(jc+j) = bjb + bjf
                bjd = B(jd+j)
                bjh = B(jh+j)
                u3 = bjd - bjh
                B(jg+j) = bjd + bjh
                B(jb+j) = u0
                B(jf+j) = u1
                B(jd+j) = c2*(u2-u3)
                B(jh+j) = c3*(u2+u3)
                j = j + Jump
              End Do
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
        End Do

        !     stage 2
        !     -------

        !     k=0 (no twiddle factors)
        !     ------------------------
        Do jjj = 0,(N-1)*Inc,8*jstep
          ja = istart + jjj

          ! "transverse" loop
          ! -----------------
          Do nu = 1,inq
            jb = ja + jstepl
            If ( jb < istart ) jb = jb + ninc
            jc = jb + jstepl
            If ( jc < istart ) jc = jc + ninc
            jd = jc + jstepl
            If ( jd < istart ) jd = jd + ninc
            je = jd + jstepl
            If ( je < istart ) je = je + ninc
            jf = je + jstepl
            If ( jf < istart ) jf = jf + ninc
            jg = jf + jstepl
            If ( jg < istart ) jg = jg + ninc
            jh = jg + jstepl
            If ( jh < istart ) jh = jh + ninc
            j = 0
            !cdir$ ivdep, shortloop
            Do l = 1,nvex
              aja = A(ja+j)
              aje = A(je+j)
              t0 = aja + aje
              t2 = aja - aje
              ajc = A(jc+j)
              ajg = A(jg+j)
              t1 = ajc + ajg
              t3 = c1*(ajc-ajg)
              bja = B(ja+j)
              bje = B(je+j)
              u0 = bja + bje
              u2 = bja - bje
              bjc = B(jc+j)
              bjg = B(jg+j)
              u1 = bjc + bjg
              u3 = c1*(bjc-bjg)
              A(ja+j) = t0 + t1
              A(je+j) = t0 - t1
              B(ja+j) = u0 + u1
              B(je+j) = u0 - u1
              A(jc+j) = t2 - u3
              A(jg+j) = t2 + u3
              B(jc+j) = u2 + t3
              B(jg+j) = u2 - t3
              ajb = A(jb+j)
              ajd = A(jd+j)
              t0 = ajb + ajd
              t2 = ajb - ajd
              ajf = A(jf+j)
              ajh = A(jh+j)
              t1 = ajf - ajh
              t3 = ajf + ajh
              bjb = B(jb+j)
              bjd = B(jd+j)
              u0 = bjb + bjd
              u2 = bjb - bjd
              bjf = B(jf+j)
              bjh = B(jh+j)
              u1 = bjf - bjh
              u3 = bjf + bjh
              A(jb+j) = t0 - u3
              A(jh+j) = t0 + u3
              B(jb+j) = u0 + t3
              B(jh+j) = u0 - t3
              A(jd+j) = t2 + u1
              A(jf+j) = t2 - u1
              B(jd+j) = u2 - t1
              B(jf+j) = u2 + t1
              j = j + Jump
            End Do
            ja = ja + jstepx
            If ( ja < istart ) ja = ja + ninc
          End Do
        End Do

        If ( n2 == 8 ) Go To 50

        !     loop on non-zero k
        !     -----------------
        kk = 2*la

        Do k = ink,jstep - ink,ink

          co1 = Trigs(kk+1)
          si1 = s*Trigs(kk+2)
          co2 = Trigs(2*kk+1)
          si2 = s*Trigs(2*kk+2)
          co3 = Trigs(3*kk+1)
          si3 = s*Trigs(3*kk+2)
          co4 = Trigs(4*kk+1)
          si4 = s*Trigs(4*kk+2)
          co5 = Trigs(5*kk+1)
          si5 = s*Trigs(5*kk+2)
          co6 = Trigs(6*kk+1)
          si6 = s*Trigs(6*kk+2)
          co7 = Trigs(7*kk+1)
          si7 = s*Trigs(7*kk+2)

          Do jjj = k,(N-1)*Inc,8*jstep
            ja = istart + jjj

            ! "transverse" loop
            ! -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              jd = jc + jstepl
              If ( jd < istart ) jd = jd + ninc
              je = jd + jstepl
              If ( je < istart ) je = je + ninc
              jf = je + jstepl
              If ( jf < istart ) jf = jf + ninc
              jg = jf + jstepl
              If ( jg < istart ) jg = jg + ninc
              jh = jg + jstepl
              If ( jh < istart ) jh = jh + ninc
              j = 0
              !cdir$ ivdep, shortloop
              Do l = 1,nvex
                aja = A(ja+j)
                aje = A(je+j)
                t0 = aja + aje
                t2 = aja - aje
                ajc = A(jc+j)
                ajg = A(jg+j)
                t1 = ajc + ajg
                t3 = c1*(ajc-ajg)
                bja = B(ja+j)
                bje = B(je+j)
                u0 = bja + bje
                u2 = bja - bje
                bjc = B(jc+j)
                bjg = B(jg+j)
                u1 = bjc + bjg
                u3 = c1*(bjc-bjg)
                A(ja+j) = t0 + t1
                B(ja+j) = u0 + u1
                A(je+j) = co4*(t0-t1) - si4*(u0-u1)
                B(je+j) = si4*(t0-t1) + co4*(u0-u1)
                A(jc+j) = co2*(t2-u3) - si2*(u2+t3)
                B(jc+j) = si2*(t2-u3) + co2*(u2+t3)
                A(jg+j) = co6*(t2+u3) - si6*(u2-t3)
                B(jg+j) = si6*(t2+u3) + co6*(u2-t3)
                ajb = A(jb+j)
                ajd = A(jd+j)
                t0 = ajb + ajd
                t2 = ajb - ajd
                ajf = A(jf+j)
                ajh = A(jh+j)
                t1 = ajf - ajh
                t3 = ajf + ajh
                bjb = B(jb+j)
                bjd = B(jd+j)
                u0 = bjb + bjd
                u2 = bjb - bjd
                bjf = B(jf+j)
                bjh = B(jh+j)
                u1 = bjf - bjh
                u3 = bjf + bjh
                A(jb+j) = co1*(t0-u3) - si1*(u0+t3)
                B(jb+j) = si1*(t0-u3) + co1*(u0+t3)
                A(jh+j) = co7*(t0+u3) - si7*(u0-t3)
                B(jh+j) = si7*(t0+u3) + co7*(u0-t3)
                A(jd+j) = co3*(t2+u1) - si3*(u2-t1)
                B(jd+j) = si3*(t2+u1) + co3*(u2-t1)
                A(jf+j) = co5*(t2-u1) - si5*(u2+t1)
                B(jf+j) = si5*(t2-u1) + co5*(u2+t1)
                j = j + Jump
              End Do
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
          kk = kk + 2*la
        End Do

        la = 8*la
      End If

      !     loop on type II radix-4 passes
      !     ------------------------------
      mu = Mod(inq,4)
      If ( Isign == -1 ) mu = 4 - mu
      ss = 1.0_wp
      If ( mu == 3 ) ss = -1.0_wp

      Do ipass = mh + 1,m
        jstep = (N*Inc)/(4*la)
        jstepl = jstep - ninc
        laincl = la*ink - ninc

        !     k=0 loop (no twiddle factors)
        !     -----------------------------
        Do ll = 0,(la-1)*ink,4*jstep

          Do jjj = ll,(N-1)*Inc,4*la*ink
            ja = istart + jjj

            ! "transverse" loop
            ! -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              jd = jc + jstepl
              If ( jd < istart ) jd = jd + ninc
              je = ja + laincl
              If ( je < istart ) je = je + ninc
              jf = je + jstepl
              If ( jf < istart ) jf = jf + ninc
              jg = jf + jstepl
              If ( jg < istart ) jg = jg + ninc
              jh = jg + jstepl
              If ( jh < istart ) jh = jh + ninc
              ji = je + laincl
              If ( ji < istart ) ji = ji + ninc
              jj = ji + jstepl
              If ( jj < istart ) jj = jj + ninc
              jk = jj + jstepl
              If ( jk < istart ) jk = jk + ninc
              jl = jk + jstepl
              If ( jl < istart ) jl = jl + ninc
              jm = ji + laincl
              If ( jm < istart ) jm = jm + ninc
              jn = jm + jstepl
              If ( jn < istart ) jn = jn + ninc
              jo = jn + jstepl
              If ( jo < istart ) jo = jo + ninc
              jp = jo + jstepl
              If ( jp < istart ) jp = jp + ninc
              j = 0

              !     loop across transforms
              !     ----------------------
              !cdir$ ivdep, shortloop
              Do l = 1,nvex
                aja = A(ja+j)
                ajc = A(jc+j)
                t0 = aja + ajc
                t2 = aja - ajc
                ajb = A(jb+j)
                ajd = A(jd+j)
                t1 = ajb + ajd
                t3 = ss*(ajb-ajd)
                aji = A(ji+j)
                ajc = aji
                bja = B(ja+j)
                bjc = B(jc+j)
                u0 = bja + bjc
                u2 = bja - bjc
                bjb = B(jb+j)
                bjd = B(jd+j)
                u1 = bjb + bjd
                u3 = ss*(bjb-bjd)
                aje = A(je+j)
                ajb = aje
                A(ja+j) = t0 + t1
                A(ji+j) = t0 - t1
                B(ja+j) = u0 + u1
                bjc = u0 - u1
                bjm = B(jm+j)
                bjd = bjm
                A(je+j) = t2 - u3
                ajd = t2 + u3
                bjb = u2 + t3
                B(jm+j) = u2 - t3
                ! ----------------------
                ajg = A(jg+j)
                t0 = ajb + ajg
                t2 = ajb - ajg
                ajf = A(jf+j)
                ajh = A(jh+j)
                t1 = ajf + ajh
                t3 = ss*(ajf-ajh)
                ajj = A(jj+j)
                ajg = ajj
                bje = B(je+j)
                bjg = B(jg+j)
                u0 = bje + bjg
                u2 = bje - bjg
                bjf = B(jf+j)
                bjh = B(jh+j)
                u1 = bjf + bjh
                u3 = ss*(bjf-bjh)
                B(je+j) = bjb
                A(jb+j) = t0 + t1
                A(jj+j) = t0 - t1
                bjj = B(jj+j)
                bjg = bjj
                B(jb+j) = u0 + u1
                B(jj+j) = u0 - u1
                A(jf+j) = t2 - u3
                ajh = t2 + u3
                B(jf+j) = u2 + t3
                bjh = u2 - t3
                ! ----------------------
                ajk = A(jk+j)
                t0 = ajc + ajk
                t2 = ajc - ajk
                ajl = A(jl+j)
                t1 = ajg + ajl
                t3 = ss*(ajg-ajl)
                bji = B(ji+j)
                bjk = B(jk+j)
                u0 = bji + bjk
                u2 = bji - bjk
                ajo = A(jo+j)
                ajl = ajo
                bjl = B(jl+j)
                u1 = bjg + bjl
                u3 = ss*(bjg-bjl)
                B(ji+j) = bjc
                A(jc+j) = t0 + t1
                A(jk+j) = t0 - t1
                bjo = B(jo+j)
                bjl = bjo
                B(jc+j) = u0 + u1
                B(jk+j) = u0 - u1
                A(jg+j) = t2 - u3
                A(jo+j) = t2 + u3
                B(jg+j) = u2 + t3
                B(jo+j) = u2 - t3
                ! ----------------------
                ajm = A(jm+j)
                t0 = ajm + ajl
                t2 = ajm - ajl
                ajn = A(jn+j)
                ajp = A(jp+j)
                t1 = ajn + ajp
                t3 = ss*(ajn-ajp)
                A(jm+j) = ajd
                u0 = bjd + bjl
                u2 = bjd - bjl
                bjn = B(jn+j)
                bjp = B(jp+j)
                u1 = bjn + bjp
                u3 = ss*(bjn-bjp)
                A(jn+j) = ajh
                A(jd+j) = t0 + t1
                A(jl+j) = t0 - t1
                B(jd+j) = u0 + u1
                B(jl+j) = u0 - u1
                B(jn+j) = bjh
                A(jh+j) = t2 - u3
                A(jp+j) = t2 + u3
                B(jh+j) = u2 + t3
                B(jp+j) = u2 - t3
                j = j + Jump
              End Do
              ! -----( end of loop across transforms )
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
        End Do
        ! -----( end of double loop for k=0 )

        !     finished if last pass
        !     ---------------------
        If ( ipass == m ) Go To 50

        kk = 2*la

        ! loop on non-zero k
        ! -----------------
        Do k = ink,jstep - ink,ink
          co1 = Trigs(kk+1)
          si1 = s*Trigs(kk+2)
          co2 = Trigs(2*kk+1)
          si2 = s*Trigs(2*kk+2)
          co3 = Trigs(3*kk+1)
          si3 = s*Trigs(3*kk+2)

          !     double loop along first transform in block
          !     ------------------------------------------
          Do ll = k,(la-1)*ink,4*jstep

            Do jjj = ll,(N-1)*Inc,4*la*ink
              ja = istart + jjj

              ! "transverse" loop
              ! -----------------
              Do nu = 1,inq
                jb = ja + jstepl
                If ( jb < istart ) jb = jb + ninc
                jc = jb + jstepl
                If ( jc < istart ) jc = jc + ninc
                jd = jc + jstepl
                If ( jd < istart ) jd = jd + ninc
                je = ja + laincl
                If ( je < istart ) je = je + ninc
                jf = je + jstepl
                If ( jf < istart ) jf = jf + ninc
                jg = jf + jstepl
                If ( jg < istart ) jg = jg + ninc
                jh = jg + jstepl
                If ( jh < istart ) jh = jh + ninc
                ji = je + laincl
                If ( ji < istart ) ji = ji + ninc
                jj = ji + jstepl
                If ( jj < istart ) jj = jj + ninc
                jk = jj + jstepl
                If ( jk < istart ) jk = jk + ninc
                jl = jk + jstepl
                If ( jl < istart ) jl = jl + ninc
                jm = ji + laincl
                If ( jm < istart ) jm = jm + ninc
                jn = jm + jstepl
                If ( jn < istart ) jn = jn + ninc
                jo = jn + jstepl
                If ( jo < istart ) jo = jo + ninc
                jp = jo + jstepl
                If ( jp < istart ) jp = jp + ninc
                j = 0

                !     loop across transforms
                !     ----------------------
                !cdir$ ivdep, shortloop
                Do l = 1,nvex
                  aja = A(ja+j)
                  ajc = A(jc+j)
                  t0 = aja + ajc
                  t2 = aja - ajc
                  ajb = A(jb+j)
                  ajd = A(jd+j)
                  t1 = ajb + ajd
                  t3 = ss*(ajb-ajd)
                  aji = A(ji+j)
                  ajc = aji
                  bja = B(ja+j)
                  bjc = B(jc+j)
                  u0 = bja + bjc
                  u2 = bja - bjc
                  bjb = B(jb+j)
                  bjd = B(jd+j)
                  u1 = bjb + bjd
                  u3 = ss*(bjb-bjd)
                  aje = A(je+j)
                  ajb = aje
                  A(ja+j) = t0 + t1
                  B(ja+j) = u0 + u1
                  A(je+j) = co1*(t2-u3) - si1*(u2+t3)
                  bjb = si1*(t2-u3) + co1*(u2+t3)
                  bjm = B(jm+j)
                  bjd = bjm
                  A(ji+j) = co2*(t0-t1) - si2*(u0-u1)
                  bjc = si2*(t0-t1) + co2*(u0-u1)
                  ajd = co3*(t2+u3) - si3*(u2-t3)
                  B(jm+j) = si3*(t2+u3) + co3*(u2-t3)
                  ! ----------------------------------------
                  ajg = A(jg+j)
                  t0 = ajb + ajg
                  t2 = ajb - ajg
                  ajf = A(jf+j)
                  ajh = A(jh+j)
                  t1 = ajf + ajh
                  t3 = ss*(ajf-ajh)
                  ajj = A(jj+j)
                  ajg = ajj
                  bje = B(je+j)
                  bjg = B(jg+j)
                  u0 = bje + bjg
                  u2 = bje - bjg
                  bjf = B(jf+j)
                  bjh = B(jh+j)
                  u1 = bjf + bjh
                  u3 = ss*(bjf-bjh)
                  B(je+j) = bjb
                  A(jb+j) = t0 + t1
                  B(jb+j) = u0 + u1
                  bjj = B(jj+j)
                  bjg = bjj
                  A(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                  A(jj+j) = co2*(t0-t1) - si2*(u0-u1)
                  B(jj+j) = si2*(t0-t1) + co2*(u0-u1)
                  ajh = co3*(t2+u3) - si3*(u2-t3)
                  bjh = si3*(t2+u3) + co3*(u2-t3)
                  ! ----------------------------------------
                  ajk = A(jk+j)
                  t0 = ajc + ajk
                  t2 = ajc - ajk
                  ajl = A(jl+j)
                  t1 = ajg + ajl
                  t3 = ss*(ajg-ajl)
                  bji = B(ji+j)
                  bjk = B(jk+j)
                  u0 = bji + bjk
                  u2 = bji - bjk
                  ajo = A(jo+j)
                  ajl = ajo
                  bjl = B(jl+j)
                  u1 = bjg + bjl
                  u3 = ss*(bjg-bjl)
                  B(ji+j) = bjc
                  A(jc+j) = t0 + t1
                  B(jc+j) = u0 + u1
                  bjo = B(jo+j)
                  bjl = bjo
                  A(jg+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(jg+j) = si1*(t2-u3) + co1*(u2+t3)
                  A(jk+j) = co2*(t0-t1) - si2*(u0-u1)
                  B(jk+j) = si2*(t0-t1) + co2*(u0-u1)
                  A(jo+j) = co3*(t2+u3) - si3*(u2-t3)
                  B(jo+j) = si3*(t2+u3) + co3*(u2-t3)
                  ! ----------------------------------------
                  ajm = A(jm+j)
                  t0 = ajm + ajl
                  t2 = ajm - ajl
                  ajn = A(jn+j)
                  ajp = A(jp+j)
                  t1 = ajn + ajp
                  t3 = ss*(ajn-ajp)
                  A(jm+j) = ajd
                  u0 = bjd + bjl
                  u2 = bjd - bjl
                  A(jn+j) = ajh
                  bjn = B(jn+j)
                  bjp = B(jp+j)
                  u1 = bjn + bjp
                  u3 = ss*(bjn-bjp)
                  B(jn+j) = bjh
                  A(jd+j) = t0 + t1
                  B(jd+j) = u0 + u1
                  A(jh+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(jh+j) = si1*(t2-u3) + co1*(u2+t3)
                  A(jl+j) = co2*(t0-t1) - si2*(u0-u1)
                  B(jl+j) = si2*(t0-t1) + co2*(u0-u1)
                  A(jp+j) = co3*(t2+u3) - si3*(u2-t3)
                  B(jp+j) = si3*(t2+u3) + co3*(u2-t3)
                  j = j + Jump
                End Do
                ! -----(end of loop across transforms)
                ja = ja + jstepx
                If ( ja < istart ) ja = ja + ninc
              End Do
            End Do
          End Do
          ! -----( end of double loop for this k )
          kk = kk + 2*la
        End Do
        ! -----( end of loop over values of k )
        la = 4*la
      End Do
      ! -----( end of loop on type II radix-4 passes )
      ! -----( nvex transforms completed)
      50      istart = istart + nvex*Jump
    End Do
    ! -----( end of loop on blocks of transforms )

  End Subroutine GPFA2F



  Subroutine GPFA3F(A,B,Trigs,Inc,Jump,N,Mm,Lot,Isign)

    ! -------------------------------------------------------------------
    !
    ! FORTRAN version of *gpfa3* -
    ! radix-3 section of self-sorting, in-place
    !       generalized PFA
    !
    ! -------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    ! *                                                                *
    ! *     N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
    ! *     RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
    ! *     (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
    ! *                                                                *
    ! ----------------------------------------------------------------------

    !*** Start of declarations inserted by SPAG (plusFORT)
    Real( Kind = wp ) :: A(*),B(*),Trigs(*),                  &
      aja,ajb,ajc,ajd,aje,ajf,ajg,ajh,aji, &
      bja,bjb,bjc,bjd,bje,bjf,bjg,bjh,bji, &
      c1,co1,co2,s,si1,si2,t1,t2,t3,u1,u2,u3

    Integer           :: Inc,ink,inq,ipass,Isign,istart,j,ja, &
      jb,jc,jd,je,jf,jg,jh,ji,jjj,jstep,   &
      jstepl,jstepx,Jump,k,kk,l,la,laincl, &
      left,ll,Lot,lvr,m,mh,Mm,mu,N,n3,nb,  &
      nblox,ninc,nu,nvex

    Real( Kind = wp ), Parameter :: sin60  = 0.866025403784437_wp
    Integer,           Parameter :: ncache = 65536
    !*** End of declarations inserted by SPAG (plusFORT)

    ! set lvr dynamically - ncache = # complex elements that fit in cache *!

    If ( ncache < N ) Then
      Write(Unit=*, Fmt=*) ' ***** gpfa ***** ncache too small '
      Call ABORT_IT()
    End If
    lvr = ncache/N

    n3 = 3**Mm
    inq = N/n3
    jstepx = (n3-N)*Inc
    ninc = N*Inc
    ink = Inc*inq
    mu = Mod(inq,3)
    If ( Isign == -1 ) mu = 3 - mu
    m = Mm
    mh = (m+1)/2
    s = Real(Isign,wp)
    c1 = sin60
    If ( mu == 2 ) c1 = -c1

    nblox = 1 + (Lot-1)/lvr
    left = Lot
    s = Real(Isign,wp)
    istart = 1

    !     loop on blocks of lvr transforms
    !     --------------------------------
    Do nb = 1,nblox

      If ( left <= lvr ) Then
        nvex = left
      Else If ( left < (2*lvr) ) Then
        nvex = left/2
        nvex = nvex + Mod(nvex,2)
      Else
        nvex = lvr
      End If
      left = left - nvex

      la = 1
      !
      !     loop on type I radix-3 passes
      !     -----------------------------
      Do ipass = 1,mh
        jstep = (N*Inc)/(3*la)
        jstepl = jstep - ninc

        !     k = 0 loop (no twiddle factors)
        !     -------------------------------
        Do jjj = 0,(N-1)*Inc,3*jstep
          ja = istart + jjj

          !     "transverse" loop
          !     -----------------
          Do nu = 1,inq
            jb = ja + jstepl
            If ( jb < istart ) jb = jb + ninc
            jc = jb + jstepl
            If ( jc < istart ) jc = jc + ninc
            j = 0

            !     loop across transforms
            !     ----------------------
            !cdir$ ivdep, shortloop
            Do l = 1,nvex
              ajb = A(jb+j)
              ajc = A(jc+j)
              t1 = ajb + ajc
              aja = A(ja+j)
              t2 = aja - 0.5_wp*t1
              t3 = c1*(ajb-ajc)
              bjb = B(jb+j)
              bjc = B(jc+j)
              u1 = bjb + bjc
              bja = B(ja+j)
              u2 = bja - 0.5_wp*u1
              u3 = c1*(bjb-bjc)
              A(ja+j) = aja + t1
              B(ja+j) = bja + u1
              A(jb+j) = t2 - u3
              B(jb+j) = u2 + t3
              A(jc+j) = t2 + u3
              B(jc+j) = u2 - t3
              j = j + Jump
            End Do
            ja = ja + jstepx
            If ( ja < istart ) ja = ja + ninc
          End Do
        End Do

        !     finished if n3 = 3
        !     ------------------
        If ( n3 == 3 ) Go To 50
        kk = 2*la

        !     loop on non-zero k
        !     -----------------
        Do k = ink,jstep - ink,ink
          co1 = Trigs(kk+1)
          si1 = s*Trigs(kk+2)
          co2 = Trigs(2*kk+1)
          si2 = s*Trigs(2*kk+2)

          !     loop along transform
          !     --------------------
          Do jjj = k,(N-1)*Inc,3*jstep
            ja = istart + jjj

            !     "transverse" loop
            !     -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              j = 0

              !     loop across transforms
              !     ----------------------
              !cdir$ ivdep,shortloop
              Do l = 1,nvex
                ajb = A(jb+j)
                ajc = A(jc+j)
                t1 = ajb + ajc
                aja = A(ja+j)
                t2 = aja - 0.5_wp*t1
                t3 = c1*(ajb-ajc)
                bjb = B(jb+j)
                bjc = B(jc+j)
                u1 = bjb + bjc
                bja = B(ja+j)
                u2 = bja - 0.5_wp*u1
                u3 = c1*(bjb-bjc)
                A(ja+j) = aja + t1
                B(ja+j) = bja + u1
                A(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                B(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                A(jc+j) = co2*(t2+u3) - si2*(u2-t3)
                B(jc+j) = si2*(t2+u3) + co2*(u2-t3)
                j = j + Jump
              End Do
              ! -----( end of loop across transforms )
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
          ! -----( end of loop along transforms )
          kk = kk + 2*la
        End Do
        ! -----( end of loop on non-zero k )
        la = 3*la
      End Do
      ! -----( end of loop on type I radix-3 passes)

      !     loop on type II radix-3 passes
      !     ------------------------------

      Do ipass = mh + 1,m
        jstep = (N*Inc)/(3*la)
        jstepl = jstep - ninc
        laincl = la*ink - ninc

        !     k=0 loop (no twiddle factors)
        !     -----------------------------
        Do ll = 0,(la-1)*ink,3*jstep

          Do jjj = ll,(N-1)*Inc,3*la*ink
            ja = istart + jjj

            !     "transverse" loop
            !     -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              jd = ja + laincl
              If ( jd < istart ) jd = jd + ninc
              je = jd + jstepl
              If ( je < istart ) je = je + ninc
              jf = je + jstepl
              If ( jf < istart ) jf = jf + ninc
              jg = jd + laincl
              If ( jg < istart ) jg = jg + ninc
              jh = jg + jstepl
              If ( jh < istart ) jh = jh + ninc
              ji = jh + jstepl
              If ( ji < istart ) ji = ji + ninc
              j = 0

              !     loop across transforms
              !     ----------------------
              !cdir$ ivdep, shortloop
              Do l = 1,nvex
                ajb = A(jb+j)
                ajc = A(jc+j)
                t1 = ajb + ajc
                aja = A(ja+j)
                t2 = aja - 0.5_wp*t1
                t3 = c1*(ajb-ajc)
                ajd = A(jd+j)
                ajb = ajd
                bjb = B(jb+j)
                bjc = B(jc+j)
                u1 = bjb + bjc
                bja = B(ja+j)
                u2 = bja - 0.5_wp*u1
                u3 = c1*(bjb-bjc)
                bjd = B(jd+j)
                bjb = bjd
                A(ja+j) = aja + t1
                B(ja+j) = bja + u1
                A(jd+j) = t2 - u3
                B(jd+j) = u2 + t3
                ajc = t2 + u3
                bjc = u2 - t3
                ! ----------------------
                aje = A(je+j)
                ajf = A(jf+j)
                t1 = aje + ajf
                t2 = ajb - 0.5_wp*t1
                t3 = c1*(aje-ajf)
                ajh = A(jh+j)
                ajf = ajh
                bje = B(je+j)
                bjf = B(jf+j)
                u1 = bje + bjf
                u2 = bjb - 0.5_wp*u1
                u3 = c1*(bje-bjf)
                bjh = B(jh+j)
                bjf = bjh
                A(jb+j) = ajb + t1
                B(jb+j) = bjb + u1
                A(je+j) = t2 - u3
                B(je+j) = u2 + t3
                A(jh+j) = t2 + u3
                B(jh+j) = u2 - t3
                ! ----------------------
                aji = A(ji+j)
                t1 = ajf + aji
                ajg = A(jg+j)
                t2 = ajg - 0.5_wp*t1
                t3 = c1*(ajf-aji)
                t1 = ajg + t1
                A(jg+j) = ajc
                bji = B(ji+j)
                u1 = bjf + bji
                bjg = B(jg+j)
                u2 = bjg - 0.5_wp*u1
                u3 = c1*(bjf-bji)
                u1 = bjg + u1
                B(jg+j) = bjc
                A(jc+j) = t1
                B(jc+j) = u1
                A(jf+j) = t2 - u3
                B(jf+j) = u2 + t3
                A(ji+j) = t2 + u3
                B(ji+j) = u2 - t3
                j = j + Jump
              End Do
              ! -----( end of loop across transforms )
              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
        End Do
        ! -----( end of double loop for k=0 )

        !     finished if last pass
        !     ---------------------
        If ( ipass == m ) Go To 50

        kk = 2*la

        ! loop on non-zero k
        ! -----------------
        Do k = ink,jstep - ink,ink
          co1 = Trigs(kk+1)
          si1 = s*Trigs(kk+2)
          co2 = Trigs(2*kk+1)
          si2 = s*Trigs(2*kk+2)

          !     double loop along first transform in block
          !     ------------------------------------------
          Do ll = k,(la-1)*ink,3*jstep

            Do jjj = ll,(N-1)*Inc,3*la*ink
              ja = istart + jjj

              !     "transverse" loop
              !     -----------------
              Do nu = 1,inq
                jb = ja + jstepl
                If ( jb < istart ) jb = jb + ninc
                jc = jb + jstepl
                If ( jc < istart ) jc = jc + ninc
                jd = ja + laincl
                If ( jd < istart ) jd = jd + ninc
                je = jd + jstepl
                If ( je < istart ) je = je + ninc
                jf = je + jstepl
                If ( jf < istart ) jf = jf + ninc
                jg = jd + laincl
                If ( jg < istart ) jg = jg + ninc
                jh = jg + jstepl
                If ( jh < istart ) jh = jh + ninc
                ji = jh + jstepl
                If ( ji < istart ) ji = ji + ninc
                j = 0

                !     loop across transforms
                !     ----------------------
                !cdir$ ivdep, shortloop
                Do l = 1,nvex
                  ajb = A(jb+j)
                  ajc = A(jc+j)
                  t1 = ajb + ajc
                  aja = A(ja+j)
                  t2 = aja - 0.5_wp*t1
                  t3 = c1*(ajb-ajc)
                  ajd = A(jd+j)
                  ajb = ajd
                  bjb = B(jb+j)
                  bjc = B(jc+j)
                  u1 = bjb + bjc
                  bja = B(ja+j)
                  u2 = bja - 0.5_wp*u1
                  u3 = c1*(bjb-bjc)
                  bjd = B(jd+j)
                  bjb = bjd
                  A(ja+j) = aja + t1
                  B(ja+j) = bja + u1
                  A(jd+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(jd+j) = si1*(t2-u3) + co1*(u2+t3)
                  ajc = co2*(t2+u3) - si2*(u2-t3)
                  bjc = si2*(t2+u3) + co2*(u2-t3)
                  ! ----------------------
                  aje = A(je+j)
                  ajf = A(jf+j)
                  t1 = aje + ajf
                  t2 = ajb - 0.5_wp*t1
                  t3 = c1*(aje-ajf)
                  ajh = A(jh+j)
                  ajf = ajh
                  bje = B(je+j)
                  bjf = B(jf+j)
                  u1 = bje + bjf
                  u2 = bjb - 0.5_wp*u1
                  u3 = c1*(bje-bjf)
                  bjh = B(jh+j)
                  bjf = bjh
                  A(jb+j) = ajb + t1
                  B(jb+j) = bjb + u1
                  A(je+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(je+j) = si1*(t2-u3) + co1*(u2+t3)
                  A(jh+j) = co2*(t2+u3) - si2*(u2-t3)
                  B(jh+j) = si2*(t2+u3) + co2*(u2-t3)
                  ! ----------------------
                  aji = A(ji+j)
                  t1 = ajf + aji
                  ajg = A(jg+j)
                  t2 = ajg - 0.5_wp*t1
                  t3 = c1*(ajf-aji)
                  t1 = ajg + t1
                  A(jg+j) = ajc
                  bji = B(ji+j)
                  u1 = bjf + bji
                  bjg = B(jg+j)
                  u2 = bjg - 0.5_wp*u1
                  u3 = c1*(bjf-bji)
                  u1 = bjg + u1
                  B(jg+j) = bjc
                  A(jc+j) = t1
                  B(jc+j) = u1
                  A(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                  B(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                  A(ji+j) = co2*(t2+u3) - si2*(u2-t3)
                  B(ji+j) = si2*(t2+u3) + co2*(u2-t3)
                  j = j + Jump
                End Do
                ! -----(end of loop across transforms)
                ja = ja + jstepx
                If ( ja < istart ) ja = ja + ninc
              End Do
            End Do
          End Do
          ! -----( end of double loop for this k )
          kk = kk + 2*la
        End Do
        ! -----( end of loop over values of k )
        la = 3*la
      End Do
      ! -----( end of loop on type II radix-3 passes )
      ! -----( nvex transforms completed)
      50     istart = istart + nvex*Jump
    End Do
    ! -----( end of loop on blocks of transforms )

  End Subroutine GPFA3F



  Subroutine GPFA5F(A,B,Trigs,Inc,Jump,N,Mm,Lot,Isign)

    ! -------------------------------------------------------------------
    !
    ! FORTRAN version of *gpfa5* -
    ! radix-5 section of self-sorting, in-place,
    !       generalized pfa
    !
    ! -------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    ! *                                                                *
    ! *     N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
    ! *     RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
    ! *     (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
    ! *                                                                *
    ! ----------------------------------------------------------------------

    !*** Start of declarations inserted by SPAG (plusFORT)
    Real( Kind = wp ) :: A(*),B(*),Trigs(*),                      &
      aja,ajb,ajc,ajd,aje,ajf,ajg,ajh,aji,     &
      ajj,ajk,ajl,ajm,ajn,ajo,ajp,ajq,ajr,ajs, &
      ajt,aju,ajv,ajw,ajx,ajy,ax,bja,bjb,      &
      bjc,bjd,bje,bjf,bjg,bjh,bji,bjj,bjk,bjl, &
      bjm,bjn,bjo,bjp,bjq,bjr,bjs,bjt,bju,bjv, &
      bjw,bjx,bjy,bx,c1,c2,c3,co1,co2,co3,     &
      co4,s,si1,si2,si3,si4,                   &
      t1,t10,t11,t2,t3,t4,t5,t6,t7,t8,t9,      &
      u1,u10,u11,u2,u3,u4,u5,u6,u7,u8,u9

    Integer           :: Inc,ink,inq,ipass,Isign,istart,j,ja,jb,  &
      jc,jd,je,jf,jg,jh,ji,jj,jjj,jk,jl,       &
      jm,jn,jo,jp,jq,jr,js,jstep,jstepl,       &
      jstepx,jt,ju,Jump,jv,jw,jx,jy,k,kk,l,    &
      la,laincl,left,ll,Lot,lvr,m,mh,Mm,mu,    &
      N,n5,nb,nblox,ninc,nu,nvex

    Real( Kind = wp ), Parameter :: sin36  = 0.587785252292473_wp
    Real( Kind = wp ), Parameter :: sin72  = 0.951056516295154_wp
    Real( Kind = wp ), Parameter :: qrt5   = 0.559016994374947_wp
    Integer,           Parameter :: ncache = 65536

    !*** End of declarations inserted by SPAG (plusFORT)

    ! set lvr dynamically - ncache = # complex elements that fit in cache *!

    If ( ncache < N ) Then
      Write(Unit=*, Fmt=*) ' ***** gpfa ***** ncache too small '
      Call ABORT_IT()
    End If
    lvr = ncache/N

    n5 = 5**Mm
    inq = N/n5
    jstepx = (n5-N)*Inc
    ninc = N*Inc
    ink = Inc*inq
    mu = Mod(inq,5)
    If ( Isign == -1 ) mu = 5 - mu

    m = Mm
    mh = (m+1)/2
    s = Real(Isign,wp)
    c1 = qrt5
    c2 = sin72
    c3 = sin36
    If ( mu == 2 .or. mu.EQ.3 ) Then
      c1 = -c1
      c2 = sin36
      c3 = sin72
    End If
    If ( mu == 3 .or. mu.EQ.4 ) c2 = -c2
    If ( mu == 2 .or. mu.EQ.4 ) c3 = -c3

    nblox = 1 + (Lot-1)/lvr
    left = Lot
    s = Real(Isign,wp)
    istart = 1

    !     loop on blocks of lvr transforms
    !     --------------------------------
    Do nb = 1,nblox

      If ( left <= lvr ) Then
        nvex = left
      Else If ( left < (2*lvr) ) Then
        nvex = left/2
        nvex = nvex + Mod(nvex,2)
      Else
        nvex = lvr
      End If
      left = left - nvex

      la = 1

      !     loop on type I radix-5 passes
      !     -----------------------------
      Do ipass = 1,mh
        jstep = (N*Inc)/(5*la)
        jstepl = jstep - ninc
        kk = 0

        !     loop on k
        !     ---------
        Do k = 0,jstep - ink,ink

          If ( k > 0 ) Then
            co1 = Trigs(kk+1)
            si1 = s*Trigs(kk+2)
            co2 = Trigs(2*kk+1)
            si2 = s*Trigs(2*kk+2)
            co3 = Trigs(3*kk+1)
            si3 = s*Trigs(3*kk+2)
            co4 = Trigs(4*kk+1)
            si4 = s*Trigs(4*kk+2)
          End If

          !     loop along transform
          !     --------------------
          Do jjj = k,(N-1)*Inc,5*jstep
            ja = istart + jjj

            ! "transverse" loop
            ! -----------------
            Do nu = 1,inq
              jb = ja + jstepl
              If ( jb < istart ) jb = jb + ninc
              jc = jb + jstepl
              If ( jc < istart ) jc = jc + ninc
              jd = jc + jstepl
              If ( jd < istart ) jd = jd + ninc
              je = jd + jstepl
              If ( je < istart ) je = je + ninc
              j = 0

              !     loop across transforms
              !     ----------------------
              If ( k == 0 ) Then

                !cdir$ ivdep, shortloop
                Do l = 1,nvex
                  ajb = A(jb+j)
                  aje = A(je+j)
                  t1 = ajb + aje
                  ajc = A(jc+j)
                  ajd = A(jd+j)
                  t2 = ajc + ajd
                  t3 = ajb - aje
                  t4 = ajc - ajd
                  t5 = t1 + t2
                  t6 = c1*(t1-t2)
                  aja = A(ja+j)
                  t7 = aja - 0.25_wp*t5
                  A(ja+j) = aja + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  t10 = c3*t3 - c2*t4
                  t11 = c2*t3 + c3*t4
                  bjb = B(jb+j)
                  bje = B(je+j)
                  u1 = bjb + bje
                  bjc = B(jc+j)
                  bjd = B(jd+j)
                  u2 = bjc + bjd
                  u3 = bjb - bje
                  u4 = bjc - bjd
                  u5 = u1 + u2
                  u6 = c1*(u1-u2)
                  bja = B(ja+j)
                  u7 = bja - 0.25_wp*u5
                  B(ja+j) = bja + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  u10 = c3*u3 - c2*u4
                  u11 = c2*u3 + c3*u4
                  A(jb+j) = t8 - u11
                  B(jb+j) = u8 + t11
                  A(je+j) = t8 + u11
                  B(je+j) = u8 - t11
                  A(jc+j) = t9 - u10
                  B(jc+j) = u9 + t10
                  A(jd+j) = t9 + u10
                  B(jd+j) = u9 - t10
                  j = j + Jump
                End Do

              Else

                !cdir$ ivdep,shortloop
                Do l = 1,nvex
                  ajb = A(jb+j)
                  aje = A(je+j)
                  t1 = ajb + aje
                  ajc = A(jc+j)
                  ajd = A(jd+j)
                  t2 = ajc + ajd
                  t3 = ajb - aje
                  t4 = ajc - ajd
                  t5 = t1 + t2
                  t6 = c1*(t1-t2)
                  aja = A(ja+j)
                  t7 = aja - 0.25_wp*t5
                  A(ja+j) = aja + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  t10 = c3*t3 - c2*t4
                  t11 = c2*t3 + c3*t4
                  bjb = B(jb+j)
                  bje = B(je+j)
                  u1 = bjb + bje
                  bjc = B(jc+j)
                  bjd = B(jd+j)
                  u2 = bjc + bjd
                  u3 = bjb - bje
                  u4 = bjc - bjd
                  u5 = u1 + u2
                  u6 = c1*(u1-u2)
                  bja = B(ja+j)
                  u7 = bja - 0.25_wp*u5
                  B(ja+j) = bja + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  u10 = c3*u3 - c2*u4
                  u11 = c2*u3 + c3*u4
                  A(jb+j) = co1*(t8-u11) - si1*(u8+t11)
                  B(jb+j) = si1*(t8-u11) + co1*(u8+t11)
                  A(je+j) = co4*(t8+u11) - si4*(u8-t11)
                  B(je+j) = si4*(t8+u11) + co4*(u8-t11)
                  A(jc+j) = co2*(t9-u10) - si2*(u9+t10)
                  B(jc+j) = si2*(t9-u10) + co2*(u9+t10)
                  A(jd+j) = co3*(t9+u10) - si3*(u9-t10)
                  B(jd+j) = si3*(t9+u10) + co3*(u9-t10)
                  j = j + Jump
                End Do

              End If

              ! -----( end of loop across transforms )

              ja = ja + jstepx
              If ( ja < istart ) ja = ja + ninc
            End Do
          End Do
          ! -----( end of loop along transforms )
          kk = kk + 2*la
        End Do
        ! -----( end of loop on non-zero k )
        la = 5*la
      End Do
      ! -----( end of loop on type I radix-5 passes)

      If ( N /= 5 ) Then

        !     loop on type II radix-5 passes
        !     ------------------------------

        Do ipass = mh + 1,m
          jstep = (N*Inc)/(5*la)
          jstepl = jstep - ninc
          laincl = la*ink - ninc
          kk = 0

          ! loop on k
          ! ---------
          Do k = 0,jstep - ink,ink

            If ( k > 0 ) Then
              co1 = Trigs(kk+1)
              si1 = s*Trigs(kk+2)
              co2 = Trigs(2*kk+1)
              si2 = s*Trigs(2*kk+2)
              co3 = Trigs(3*kk+1)
              si3 = s*Trigs(3*kk+2)
              co4 = Trigs(4*kk+1)
              si4 = s*Trigs(4*kk+2)
            End If

            !     double loop along first transform in block
            !     ------------------------------------------
            Do ll = k,(la-1)*ink,5*jstep

              Do jjj = ll,(N-1)*Inc,5*la*ink
                ja = istart + jjj

                ! "transverse" loop
                ! -----------------
                Do nu = 1,inq
                  jb = ja + jstepl
                  If ( jb < istart ) jb = jb + ninc
                  jc = jb + jstepl
                  If ( jc < istart ) jc = jc + ninc
                  jd = jc + jstepl
                  If ( jd < istart ) jd = jd + ninc
                  je = jd + jstepl
                  If ( je < istart ) je = je + ninc
                  jf = ja + laincl
                  If ( jf < istart ) jf = jf + ninc
                  jg = jf + jstepl
                  If ( jg < istart ) jg = jg + ninc
                  jh = jg + jstepl
                  If ( jh < istart ) jh = jh + ninc
                  ji = jh + jstepl
                  If ( ji < istart ) ji = ji + ninc
                  jj = ji + jstepl
                  If ( jj < istart ) jj = jj + ninc
                  jk = jf + laincl
                  If ( jk < istart ) jk = jk + ninc
                  jl = jk + jstepl
                  If ( jl < istart ) jl = jl + ninc
                  jm = jl + jstepl
                  If ( jm < istart ) jm = jm + ninc
                  jn = jm + jstepl
                  If ( jn < istart ) jn = jn + ninc
                  jo = jn + jstepl
                  If ( jo < istart ) jo = jo + ninc
                  jp = jk + laincl
                  If ( jp < istart ) jp = jp + ninc
                  jq = jp + jstepl
                  If ( jq < istart ) jq = jq + ninc
                  jr = jq + jstepl
                  If ( jr < istart ) jr = jr + ninc
                  js = jr + jstepl
                  If ( js < istart ) js = js + ninc
                  jt = js + jstepl
                  If ( jt < istart ) jt = jt + ninc
                  ju = jp + laincl
                  If ( ju < istart ) ju = ju + ninc
                  jv = ju + jstepl
                  If ( jv < istart ) jv = jv + ninc
                  jw = jv + jstepl
                  If ( jw < istart ) jw = jw + ninc
                  jx = jw + jstepl
                  If ( jx < istart ) jx = jx + ninc
                  jy = jx + jstepl
                  If ( jy < istart ) jy = jy + ninc
                  j = 0

                  !     loop across transforms
                  !     ----------------------
                  If ( k == 0 ) Then

                    !cdir$ ivdep, shortloop
                    Do l = 1,nvex
                      ajb = A(jb+j)
                      aje = A(je+j)
                      t1 = ajb + aje
                      ajc = A(jc+j)
                      ajd = A(jd+j)
                      t2 = ajc + ajd
                      t3 = ajb - aje
                      t4 = ajc - ajd
                      ajf = A(jf+j)
                      ajb = ajf
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      aja = A(ja+j)
                      t7 = aja - 0.25_wp*t5
                      A(ja+j) = aja + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      ajk = A(jk+j)
                      ajc = ajk
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      bjb = B(jb+j)
                      bje = B(je+j)
                      u1 = bjb + bje
                      bjc = B(jc+j)
                      bjd = B(jd+j)
                      u2 = bjc + bjd
                      u3 = bjb - bje
                      u4 = bjc - bjd
                      bjf = B(jf+j)
                      bjb = bjf
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      bja = B(ja+j)
                      u7 = bja - 0.25_wp*u5
                      B(ja+j) = bja + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      bjk = B(jk+j)
                      bjc = bjk
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      A(jf+j) = t8 - u11
                      B(jf+j) = u8 + t11
                      aje = t8 + u11
                      bje = u8 - t11
                      A(jk+j) = t9 - u10
                      B(jk+j) = u9 + t10
                      ajd = t9 + u10
                      bjd = u9 - t10
                      ! ----------------------
                      ajg = A(jg+j)
                      ajj = A(jj+j)
                      t1 = ajg + ajj
                      ajh = A(jh+j)
                      aji = A(ji+j)
                      t2 = ajh + aji
                      t3 = ajg - ajj
                      t4 = ajh - aji
                      ajl = A(jl+j)
                      ajh = ajl
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      t7 = ajb - 0.25_wp*t5
                      A(jb+j) = ajb + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      ajq = A(jq+j)
                      aji = ajq
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      bjg = B(jg+j)
                      bjj = B(jj+j)
                      u1 = bjg + bjj
                      bjh = B(jh+j)
                      bji = B(ji+j)
                      u2 = bjh + bji
                      u3 = bjg - bjj
                      u4 = bjh - bji
                      bjl = B(jl+j)
                      bjh = bjl
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      u7 = bjb - 0.25_wp*u5
                      B(jb+j) = bjb + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      bjq = B(jq+j)
                      bji = bjq
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      A(jg+j) = t8 - u11
                      B(jg+j) = u8 + t11
                      ajj = t8 + u11
                      bjj = u8 - t11
                      A(jl+j) = t9 - u10
                      B(jl+j) = u9 + t10
                      A(jq+j) = t9 + u10
                      B(jq+j) = u9 - t10
                      ! ----------------------
                      ajo = A(jo+j)
                      t1 = ajh + ajo
                      ajm = A(jm+j)
                      ajn = A(jn+j)
                      t2 = ajm + ajn
                      t3 = ajh - ajo
                      t4 = ajm - ajn
                      ajr = A(jr+j)
                      ajn = ajr
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      t7 = ajc - 0.25_wp*t5
                      A(jc+j) = ajc + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      ajw = A(jw+j)
                      ajo = ajw
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      bjo = B(jo+j)
                      u1 = bjh + bjo
                      bjm = B(jm+j)
                      bjn = B(jn+j)
                      u2 = bjm + bjn
                      u3 = bjh - bjo
                      u4 = bjm - bjn
                      bjr = B(jr+j)
                      bjn = bjr
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      u7 = bjc - 0.25_wp*u5
                      B(jc+j) = bjc + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      bjw = B(jw+j)
                      bjo = bjw
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      A(jh+j) = t8 - u11
                      B(jh+j) = u8 + t11
                      A(jw+j) = t8 + u11
                      B(jw+j) = u8 - t11
                      A(jm+j) = t9 - u10
                      B(jm+j) = u9 + t10
                      A(jr+j) = t9 + u10
                      B(jr+j) = u9 - t10
                      ! ----------------------
                      ajt = A(jt+j)
                      t1 = aji + ajt
                      ajs = A(js+j)
                      t2 = ajn + ajs
                      t3 = aji - ajt
                      t4 = ajn - ajs
                      ajx = A(jx+j)
                      ajt = ajx
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      ajp = A(jp+j)
                      t7 = ajp - 0.25_wp*t5
                      ax = ajp + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      A(jp+j) = ajd
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      A(jd+j) = ax
                      bjt = B(jt+j)
                      u1 = bji + bjt
                      bjs = B(js+j)
                      u2 = bjn + bjs
                      u3 = bji - bjt
                      u4 = bjn - bjs
                      bjx = B(jx+j)
                      bjt = bjx
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      bjp = B(jp+j)
                      u7 = bjp - 0.25_wp*u5
                      bx = bjp + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      B(jp+j) = bjd
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      B(jd+j) = bx
                      A(ji+j) = t8 - u11
                      B(ji+j) = u8 + t11
                      A(jx+j) = t8 + u11
                      B(jx+j) = u8 - t11
                      A(jn+j) = t9 - u10
                      B(jn+j) = u9 + t10
                      A(js+j) = t9 + u10
                      B(js+j) = u9 - t10
                      ! ----------------------
                      ajv = A(jv+j)
                      ajy = A(jy+j)
                      t1 = ajv + ajy
                      t2 = ajo + ajt
                      t3 = ajv - ajy
                      t4 = ajo - ajt
                      A(jv+j) = ajj
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      aju = A(ju+j)
                      t7 = aju - 0.25_wp*t5
                      ax = aju + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      A(ju+j) = aje
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      A(je+j) = ax
                      bjv = B(jv+j)
                      bjy = B(jy+j)
                      u1 = bjv + bjy
                      u2 = bjo + bjt
                      u3 = bjv - bjy
                      u4 = bjo - bjt
                      B(jv+j) = bjj
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      bju = B(ju+j)
                      u7 = bju - 0.25_wp*u5
                      bx = bju + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      B(ju+j) = bje
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      B(je+j) = bx
                      A(jj+j) = t8 - u11
                      B(jj+j) = u8 + t11
                      A(jy+j) = t8 + u11
                      B(jy+j) = u8 - t11
                      A(jo+j) = t9 - u10
                      B(jo+j) = u9 + t10
                      A(jt+j) = t9 + u10
                      B(jt+j) = u9 - t10
                      j = j + Jump
                    End Do

                  Else

                    !cdir$ ivdep, shortloop
                    Do l = 1,nvex
                      ajb = A(jb+j)
                      aje = A(je+j)
                      t1 = ajb + aje
                      ajc = A(jc+j)
                      ajd = A(jd+j)
                      t2 = ajc + ajd
                      t3 = ajb - aje
                      t4 = ajc - ajd
                      ajf = A(jf+j)
                      ajb = ajf
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      aja = A(ja+j)
                      t7 = aja - 0.25_wp*t5
                      A(ja+j) = aja + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      ajk = A(jk+j)
                      ajc = ajk
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      bjb = B(jb+j)
                      bje = B(je+j)
                      u1 = bjb + bje
                      bjc = B(jc+j)
                      bjd = B(jd+j)
                      u2 = bjc + bjd
                      u3 = bjb - bje
                      u4 = bjc - bjd
                      bjf = B(jf+j)
                      bjb = bjf
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      bja = B(ja+j)
                      u7 = bja - 0.25_wp*u5
                      B(ja+j) = bja + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      bjk = B(jk+j)
                      bjc = bjk
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      A(jf+j) = co1*(t8-u11) - si1*(u8+t11)
                      B(jf+j) = si1*(t8-u11) + co1*(u8+t11)
                      aje = co4*(t8+u11) - si4*(u8-t11)
                      bje = si4*(t8+u11) + co4*(u8-t11)
                      A(jk+j) = co2*(t9-u10) - si2*(u9+t10)
                      B(jk+j) = si2*(t9-u10) + co2*(u9+t10)
                      ajd = co3*(t9+u10) - si3*(u9-t10)
                      bjd = si3*(t9+u10) + co3*(u9-t10)
                      ! ----------------------
                      ajg = A(jg+j)
                      ajj = A(jj+j)
                      t1 = ajg + ajj
                      ajh = A(jh+j)
                      aji = A(ji+j)
                      t2 = ajh + aji
                      t3 = ajg - ajj
                      t4 = ajh - aji
                      ajl = A(jl+j)
                      ajh = ajl
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      t7 = ajb - 0.25_wp*t5
                      A(jb+j) = ajb + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      ajq = A(jq+j)
                      aji = ajq
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      bjg = B(jg+j)
                      bjj = B(jj+j)
                      u1 = bjg + bjj
                      bjh = B(jh+j)
                      bji = B(ji+j)
                      u2 = bjh + bji
                      u3 = bjg - bjj
                      u4 = bjh - bji
                      bjl = B(jl+j)
                      bjh = bjl
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      u7 = bjb - 0.25_wp*u5
                      B(jb+j) = bjb + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      bjq = B(jq+j)
                      bji = bjq
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      A(jg+j) = co1*(t8-u11) - si1*(u8+t11)
                      B(jg+j) = si1*(t8-u11) + co1*(u8+t11)
                      ajj = co4*(t8+u11) - si4*(u8-t11)
                      bjj = si4*(t8+u11) + co4*(u8-t11)
                      A(jl+j) = co2*(t9-u10) - si2*(u9+t10)
                      B(jl+j) = si2*(t9-u10) + co2*(u9+t10)
                      A(jq+j) = co3*(t9+u10) - si3*(u9-t10)
                      B(jq+j) = si3*(t9+u10) + co3*(u9-t10)
                      ! ----------------------
                      ajo = A(jo+j)
                      t1 = ajh + ajo
                      ajm = A(jm+j)
                      ajn = A(jn+j)
                      t2 = ajm + ajn
                      t3 = ajh - ajo
                      t4 = ajm - ajn
                      ajr = A(jr+j)
                      ajn = ajr
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      t7 = ajc - 0.25_wp*t5
                      A(jc+j) = ajc + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      ajw = A(jw+j)
                      ajo = ajw
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      bjo = B(jo+j)
                      u1 = bjh + bjo
                      bjm = B(jm+j)
                      bjn = B(jn+j)
                      u2 = bjm + bjn
                      u3 = bjh - bjo
                      u4 = bjm - bjn
                      bjr = B(jr+j)
                      bjn = bjr
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      u7 = bjc - 0.25_wp*u5
                      B(jc+j) = bjc + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      bjw = B(jw+j)
                      bjo = bjw
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      A(jh+j) = co1*(t8-u11) - si1*(u8+t11)
                      B(jh+j) = si1*(t8-u11) + co1*(u8+t11)
                      A(jw+j) = co4*(t8+u11) - si4*(u8-t11)
                      B(jw+j) = si4*(t8+u11) + co4*(u8-t11)
                      A(jm+j) = co2*(t9-u10) - si2*(u9+t10)
                      B(jm+j) = si2*(t9-u10) + co2*(u9+t10)
                      A(jr+j) = co3*(t9+u10) - si3*(u9-t10)
                      B(jr+j) = si3*(t9+u10) + co3*(u9-t10)
                      ! ----------------------
                      ajt = A(jt+j)
                      t1 = aji + ajt
                      ajs = A(js+j)
                      t2 = ajn + ajs
                      t3 = aji - ajt
                      t4 = ajn - ajs
                      ajx = A(jx+j)
                      ajt = ajx
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      ajp = A(jp+j)
                      t7 = ajp - 0.25_wp*t5
                      ax = ajp + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      A(jp+j) = ajd
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      A(jd+j) = ax
                      bjt = B(jt+j)
                      u1 = bji + bjt
                      bjs = B(js+j)
                      u2 = bjn + bjs
                      u3 = bji - bjt
                      u4 = bjn - bjs
                      bjx = B(jx+j)
                      bjt = bjx
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      bjp = B(jp+j)
                      u7 = bjp - 0.25_wp*u5
                      bx = bjp + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      B(jp+j) = bjd
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      B(jd+j) = bx
                      A(ji+j) = co1*(t8-u11) - si1*(u8+t11)
                      B(ji+j) = si1*(t8-u11) + co1*(u8+t11)
                      A(jx+j) = co4*(t8+u11) - si4*(u8-t11)
                      B(jx+j) = si4*(t8+u11) + co4*(u8-t11)
                      A(jn+j) = co2*(t9-u10) - si2*(u9+t10)
                      B(jn+j) = si2*(t9-u10) + co2*(u9+t10)
                      A(js+j) = co3*(t9+u10) - si3*(u9-t10)
                      B(js+j) = si3*(t9+u10) + co3*(u9-t10)
                      ! ----------------------
                      ajv = A(jv+j)
                      ajy = A(jy+j)
                      t1 = ajv + ajy
                      t2 = ajo + ajt
                      t3 = ajv - ajy
                      t4 = ajo - ajt
                      A(jv+j) = ajj
                      t5 = t1 + t2
                      t6 = c1*(t1-t2)
                      aju = A(ju+j)
                      t7 = aju - 0.25_wp*t5
                      ax = aju + t5
                      t8 = t7 + t6
                      t9 = t7 - t6
                      A(ju+j) = aje
                      t10 = c3*t3 - c2*t4
                      t11 = c2*t3 + c3*t4
                      A(je+j) = ax
                      bjv = B(jv+j)
                      bjy = B(jy+j)
                      u1 = bjv + bjy
                      u2 = bjo + bjt
                      u3 = bjv - bjy
                      u4 = bjo - bjt
                      B(jv+j) = bjj
                      u5 = u1 + u2
                      u6 = c1*(u1-u2)
                      bju = B(ju+j)
                      u7 = bju - 0.25_wp*u5
                      bx = bju + u5
                      u8 = u7 + u6
                      u9 = u7 - u6
                      B(ju+j) = bje
                      u10 = c3*u3 - c2*u4
                      u11 = c2*u3 + c3*u4
                      B(je+j) = bx
                      A(jj+j) = co1*(t8-u11) - si1*(u8+t11)
                      B(jj+j) = si1*(t8-u11) + co1*(u8+t11)
                      A(jy+j) = co4*(t8+u11) - si4*(u8-t11)
                      B(jy+j) = si4*(t8+u11) + co4*(u8-t11)
                      A(jo+j) = co2*(t9-u10) - si2*(u9+t10)
                      B(jo+j) = si2*(t9-u10) + co2*(u9+t10)
                      A(jt+j) = co3*(t9+u10) - si3*(u9-t10)
                      B(jt+j) = si3*(t9+u10) + co3*(u9-t10)
                      j = j + Jump
                    End Do

                  End If

                  ! -----(end of loop across transforms)

                  ja = ja + jstepx
                  If ( ja < istart ) ja = ja + ninc
                End Do
              End Do
            End Do
            ! -----( end of double loop for this k )
            kk = kk + 2*la
          End Do
          ! -----( end of loop over values of k )
          la = 5*la
        End Do
      End If
      ! -----( end of loop on type II radix-5 passes )
      ! -----( nvex transforms completed)
      istart = istart + nvex*Jump
    End Do
    ! -----( end of loop on blocks of transforms )

  End Subroutine GPFA5F



  Subroutine SETGPFA(Trigs,N)

    ! ----------------------------------------------------------------------
    !       SETUP ROUTINE FOR SELF-SORTING IN-PLACE
    !           GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
    !
    !       Call SETGPFA(TRIGS,N)
    !
    !       INPUT :
    !       -----
    !       N IS THE LENGTH OF THE TRANSFORMS. N MUST BE OF THE FORM:
    !         -----------------------------------
    !           N = (2**IP) * (3**IQ) * (5**IR)
    !         -----------------------------------
    !
    !       OUTPUT:
    !       ------
    !       TRIGS IS A TABLE OF TWIDDLE FACTORS,
    !         OF LENGTH 2*IPQR (REAL) WORDS, WHERE:
    !         --------------------------------------
    !           IPQR = (2**IP) + (3**IQ) + (5**IR)
    !         --------------------------------------
    !
    !       WRITTEN BY CLIVE TEMPERTON 1990
    ! ----------------------------------------------------------------------

    !*** Start of declarations inserted by SPAG (plusFORT)
    !     Real( Kind = wp ) :: Trigs(*),angle,del,twopi
    Complex( Kind = wp ) :: Trigs(*)
    Real( Kind = wp )    :: angle,del,twopi

    Integer              :: i,ifac,li,iq,ir,irot,k,kink,kk,ll,N,ni,nj(3),nn
    !*** End of declarations inserted by SPAG (plusFORT)

    ! DECOMPOSE N INTO FACTORS 2,3,5
    nn = N
    ifac = 2
    Do ll = 1,3
      kk = 0
      50      Continue
      If ( Mod(nn,ifac) /= 0 ) Then
        nj(ll) = kk
        ifac = ifac + ll
      Else
        kk = kk + 1
        nn = nn/ifac
        Go To 50
      End If
    End Do

    If ( nn /= 1 ) Then
      Write(Unit=*, Fmt=100) N
      100     Format(' *** WARNING!!!',i10,' IS NOT A LEGAL VALUE OF N *** !!!')
      Return
    End If

    li = nj(1)
    iq = nj(2)
    ir = nj(3)

    ! COMPUTE LIST OF ROTATED TWIDDLE FACTORS
    ! ---------------------------------------

    nj(1) = 2**li
    nj(2) = 3**iq
    nj(3) = 5**ir

    twopi = 4.0_wp*Asin(1.0_wp)
    i = 1

    Do ll = 1,3
      ni = nj(ll)
      If ( ni /= 1 ) Then

        del = twopi/Real(ni,wp)
        irot = N/ni
        kink = Mod(irot,ni)
        kk = 0

        Do k = 1,ni
          angle = Real(kk,wp)*del
          !              Trigs(i) = Cos(angle)
          !              Trigs(i+1) = Sin(angle)
          !              i = i + 2
          Trigs( i ) = Exp( Cmplx( 0.0_wp, angle, Kind = wp ) )
          i = i + 1
          kk = kk + kink
          If ( kk > ni ) kk = kk - ni
        End Do
      End If
    End Do

  End Subroutine SETGPFA

  Subroutine ABORT_IT()

    Character( Len = 256 ) :: message
    Write(message, Fmt=100)
    100  Format(/'************************************'/ &
      & '*    CODE is aborting from GPFA    *'/ &
      & '************************************'/)
    Call error(0,message)

  End Subroutine ABORT_IT


End Module gpfa235


!!!!!!!!!!!!!! This is horrible... relies on implicit conversion from complex to real...
Subroutine gpfa_wrap( a, trigs, stride, jump, n, n_ffts, direction, start )
  Use kinds, Only : wp
  Use gpfa235, Only : gpfa
  Implicit None

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine wrapper for GPFA going into DaFT package
  !
  ! copyright - daresbury laboratory
  ! author    - i.j.bush march 1999
  ! amended   - i.t.todorov november 2009
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Real( Kind = wp ), Dimension( * ), Intent( InOut ) :: a
  Real( Kind = wp ), Dimension( * ), Intent( In    ) :: trigs
  Integer,                           Intent( In    ) :: stride
  Integer,                           Intent( In    ) :: jump
  Integer,                           Intent( In    ) :: n
  Integer,                           Intent( In    ) :: n_ffts
  Integer,                           Intent( In    ) :: direction
  Integer,                           Intent( In    ) :: start

  Call gpfa( a( start ), a( start + 1 ), trigs, stride, jump, n, n_ffts, direction )

End Subroutine gpfa_wrap

