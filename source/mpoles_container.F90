Module mpoles_container
  Use kinds, Only : wp
  Use comms, Only : comms_type
  Use setup, Only : mxompl, sqrpi, mximpl, mxspl, mxompl, nrite
  Use configuration, Only : imcon,cell,lsite,nlast,lsi,lsa,ltg,xxx,yyy,zzz
  Use mpole, Only : ltpatm,mpllfr,mprotm,mplltg,mplflg, &
                     mplgfr,mprotx,mproty,mprotz,mplmap


  Use ewald, Only : dtpbsp
  Implicit None

  Private

  Public :: coul_deriv, ewald_deriv, limit_erfr_deriv, explicit_ewald_real_loops, &
            explicit_fscp_rfp_loops, explicit_spme_loops, explicit_spme_loop_s, &
            rotate_mpoles, infinitesimal_rotation, rotate_mpoles_d
  !!!!!!!!!!!!!!!!!!!!!!!! THIS IS MPOLES_CONTAINER !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Subroutine coul_deriv - computes the derivatives of the 1/r kernel
  !                         for multipolar interactions
  !
  ! Subroutine ewald_deriv - computes the derivatives for a modification
  !                          of the Ewald sum for the real space kernel -
  !                          sqrt(pi)/2*erfc(alpha*r)/(alpha*r)
  !                          and the excluded reciprocal space kernel -
  !                          sqrt(pi)/2*erf(alpha*r)/(alpha*r)
  !
  ! Subroutine limit_erfr_deriv - computes the limit of the derivatives of
  !                               erf(alpha*r)/r as r->0
  !
  ! Subroutine explicit_ewald_real_loops - computes energy, forces and torques
  !                                        for multipolar interactions in real
  !                                        space and in excluded interactions
  !                                        reciprocal space without loops
  !
  ! Subroutine explicit_fscp_rfp_loops - computes energy, forces and torques
  !                                      for multipolar interactions for
  !                                      coul_fscp_mpol and coul_rfp_mpol
  !                                      without loops
  !
  ! Subroutine explicit_spme_loops - computes energy, forces and torques
  !                                  for multipolar interactions in reciprocal
  !                                  space without loops
  !
  ! Subroutine explicit_spme_loop_s - computes field components
  !                                   for multipolar interactions in reciprocal
  !                                   space without loops
  !
  ! Subroutine rotate_mpoles - rotates a multipole in the local molecular
  !                            frame into the global frame
  !
  ! Subroutine infinitesimal_rotation - performs infinitesimal counter-clockwise
  !                                     rotations of multipoles around x, y, & z
  !                                     axes.
  !
  ! Subroutine rotate_mpoles_d - rotates a multipole in the local molecular frame
  !                              into the global frame for multipole order <= 2
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Contains

  Subroutine coul_deriv(nu,torderlim,tx,ty,tz,tr,d1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for computing the derivatives of the 1/r kernel for
  ! multipolar interactions.  If the order of the multipole is 'p' then
  ! derivatives are computed up to 2p+1 in order to be able to compute the
  ! force arising from the p-p multipole interactions.
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng march 2014
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent( In    ) :: nu,torderlim
    Real( Kind = wp), Intent( In    ) :: tx,ty,tz,tr
    Real( Kind = wp), Intent(   Out ) :: d1(-2:torderlim,-2:torderlim,-2:torderlim)

  ! Local variables

    Integer          :: k1,k2,k3,k21,k22,k31,k32,nun,n
    Real( Kind = wp) :: rnu2,fac,rk2,rk221,rk3,rk331
    Real( Kind = wp) :: cf1(1:torderlim),cf2(1:torderlim)

  ! setup variables

    nun  = -nu ; rnu2 = Real(2-nu,Kind=wp)
    fac = 1.0_wp/tr**2

  ! First initialize all the derivatives to zero

    d1 = 0.0_wp

  ! Coefficients in recurrence relation

    Do n=1,torderlim
       cf1(n) = fac*rnu2/Real(n,Kind=wp) - 2.0_wp*fac
       cf2(n) = fac*rnu2/Real(n,Kind=wp) - fac
    End Do

    d1(0,0,0)=tr**nun

    Do k1=1,torderlim
       d1(k1,0,0) = cf1(k1) * Real(k1,Kind=wp)*tx*d1(k1-1,0,0) + &
                    cf2(k1) * Real(k1*(k1-1),Kind=wp)*d1(k1-2,0,0)
    End Do

    k3 = 0
    Do k2=1,torderlim
       k21 = k2-1
       k22 = k2-2

       rk2   = Real(k2,Kind=wp)*ty
       rk221 = Real(k2*k21,Kind=wp)

       Do k1=0,torderlim-k2
          n=k1+k2

          d1(k1,k2,0) = cf1(n) * (Real(k1,Kind=wp)*tx*d1(k1-1,k2,0) + rk2*d1(k1,k21,0)) + &
                        cf2(n) * (Real(k1*(k1-1),Kind=wp)*d1(k1-2,k2,0) + rk221*d1(k1,k22,0))
       End Do
    End Do

    Do k3=1,torderlim
       k31 = k3-1
       k32 = k3-2

       rk3   = Real(k3,Kind=wp)*tz
       rk331 = Real(k3*k31,Kind=wp)

       Do k2=0,torderlim-k3
          k21 = k2-1
          k22 = k2-2

          rk2   = Real(k2,Kind=wp)*ty
          rk221 = Real(k2*k21,Kind=wp)

          Do k1=0,torderlim-k3-k2
             n = k1+k2+k3

             d1(k1,k2,k3) = cf1(n) * (Real(k1,Kind=wp)*tx*d1(k1-1,k2,k3) + rk2*d1(k1,k21,k3) + rk3*d1(k1,k2,k31)) + &
                            cf2(n) * (Real(k1*(k1-1),Kind=wp)*d1(k1-2,k2,k3) + rk221*d1(k1,k22,k3) + rk331*d1(k1,k2,k32))
          End Do
       End Do
    End Do

  End Subroutine coul_deriv

  Subroutine ewald_deriv(lbl,torderlim,itype,funct,tx,ty,tz,tr,a1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for computing the derivatives of the ewald
  ! real space kernel - sqrt(pi)/2*erfc(alpha*r)/(alpha*r) -
  ! and the excluded interactions reciprocal space kernel -
  ! sqrt(pi)/2*erf(alpha*r)/(alpha*r) -
  ! for multipolar interactions.  If the order of the multipole is 'p'
  ! then derivatives are computed up to 2p+1 in order to be able to compute
  ! the force arising from the p-p multipole interactions.
  !
  ! Note: itype = 1 for the complementary error function
  !       itype = 2 for the error function
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov april 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent( In    ) :: lbl,torderlim,itype
    Real( Kind = wp), Intent( In    ) :: funct,tx,ty,tz,tr
    Real( Kind = wp), Intent(   Out ) :: a1(lbl:torderlim,lbl:torderlim,lbl:torderlim)

  ! Local variables

    Integer          :: k1,k2,k3,k11,k12,k21,k22,k31,k32,n
    Real( Kind = wp) :: rsq,fac,rk1,rk111,rk2,rk221,rk3,rk331,tt
    Real( Kind = wp) :: tdx,tdy,tdz,tc,tc1,tc2,tc3,ta
    Real( Kind = wp) :: c1(-2:torderlim,-2:torderlim,-2:torderlim)
    Real( Kind = wp) :: d1(-2:torderlim,-2:torderlim,-2:torderlim)
    Real( Kind = wp) :: cf1(1:torderlim),cf2(1:torderlim),cf3(1:torderlim)

    If (lbl /= 0 .and. lbl /= -2) STOP ! Allowed Lower_Bound_Limit
    If (itype /= 1 .and. itype /= 2) STOP ! erfcr and erfr only

  ! setup variables

    rsq = tr**2
    fac = 1.0_wp/rsq

    If (mxompl == 0) Then

       tt  =-2.0_wp*fac
       tdx = tt*tx ; tdy = tt*ty ; tdz=tt*tz

       tc  = 0.5_wp*Exp(-rsq)
       tc1 = tdx*tc ; tc2=tdy*tc ; tc3=tdz*tc

       a1(0,0,0) = 0.5_wp*sqrpi*funct

       ta = -fac*a1(0,0,0)

       If (itype == 1) Then
          a1(1,0,0) = tx*ta+tc1
          a1(0,1,0) = ty*ta+tc2
          a1(0,0,1) = tz*ta+tc3
       Else
          a1(1,0,0) = tx*ta-tc1
          a1(0,1,0) = ty*ta-tc2
          a1(0,0,1) = tz*ta-tc3
       End If

    Else

  ! First initialize all the derivatives to zero

       c1 = 0.0_wp
       d1 = 0.0_wp

  ! Coefficients in recurrence relation

       Do n=1,torderlim
          tt = fac/Real(n,Kind=wp)
          cf1(n) = tt - 2.0_wp*fac
          cf2(n) = tt - fac
          cf3(n) =-2.0_wp/Real(n,Kind=wp)
       End Do

       If (itype == 2) fac=-fac ! Only change of sign for the last term in the summations

       c1(0,0,0) = 0.5_wp*Exp(-rsq)
       d1(0,0,0) = 0.5_wp*sqrpi*funct

       Do k1=1,torderlim
          k11 = k1-1
          k12 = k1-2

          rk1   = Real(k1,Kind=wp)*tx
          rk111 = Real(k1*k11,Kind=wp)

          c1(k1,0,0) = cf3(k1) * (rk1*c1(k11,0,0) + rk111*c1(k12,0,0))

          d1(k1,0,0) = cf1(k1) * rk1*d1(k11,0,0)   + &
                       cf2(k1) * rk111*d1(k12,0,0) + &
                       fac     * c1(k1,0,0)
       End Do

       k3 = 0
       Do k2=1,torderlim
          k21 = k2-1
          k22 = k2-2

          rk2   = Real(k2,Kind=wp)*ty
          rk221 = Real(k2*k21,Kind=wp)

          Do k1=0,torderlim-k2
             n = k1+k2

             k11 = k1-1
             k12 = k1-2

             rk1   = Real(k1,Kind=wp)*tx
             rk111 = Real(k1*k11,Kind=wp)

             c1(k1,k2,0) = cf3(n) * (rk1*c1(k11,k2,0)   + rk2*c1(k1,k21,0) + &
                                     rk111*c1(k12,k2,0) + rk221*c1(k1,k22,0))

             d1(k1,k2,0) = cf1(n) * (rk1*d1(k11,k2,0)   + rk2*d1(k1,k21,0))   + &
                           cf2(n) * (rk111*d1(k12,k2,0) + rk221*d1(k1,k22,0)) + &
                           fac    * c1(k1,k2,0)
          End Do
       End Do

       Do k3=1,torderlim
          k31 = k3-1
          k32 = k3-2

          rk3   = Real(k3,Kind=wp)*tz
          rk331 = Real(k3*k31,Kind=wp)

          Do k2=0,torderlim-k3
             k21 = k2-1
             k22 = k2-2

             rk2   = Real(k2,Kind=wp)*ty
             rk221 = Real(k2*k21,Kind=wp)

             Do k1=0,torderlim-k3-k2
                n = k1+k2+k3

                k11 = k1-1
                k12 = k1-2

                rk1   = Real(k1,Kind=wp)*tx
                rk111 = Real(k1*k11,Kind=wp)

                c1(k1,k2,k3) = cf3(n) * (rk1*c1(k11,k2,k3)   + rk2*c1(k1,k21,k3)   + rk3*c1(k1,k2,k31) + &
                                         rk111*c1(k12,k2,k3) + rk221*c1(k1,k22,k3) + rk331*c1(k1,k2,k32))

                d1(k1,k2,k3) = cf1(n) * (rk1*d1(k11,k2,k3)    + rk2*d1(k1,k21,k3)   + rk3*d1(k1,k2,k31))   + &
                               cf2(n) * (rk111*d1(k1-2,k2,k3) + rk221*d1(k1,k22,k3) + rk331*d1(k1,k2,k32)) + &
                               fac    * c1(k1,k2,k3)
             End Do
          End Do
       End Do

       If (lbl == 0) Then
          a1=d1(0:torderlim,0:torderlim,0:torderlim)
       Else
          a1=d1
       End If

    End If

  End Subroutine ewald_deriv

  Subroutine limit_erfr_deriv(torderlim,alpha,d1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for computing the limit of the derivatives of
  ! erf(alpha*r)/r as r->0
  !
  ! Note: This is for the self-interaction term
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng october 2014
  ! amended   - i.t.todorov march 2015
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent( In    ) :: torderlim
    Real( Kind = wp), Intent( In    ) :: alpha
    Real( Kind = wp), Intent(   Out ) :: d1(0:torderlim,0:torderlim,0:torderlim)

  ! Local variables

    Real( Kind = wp) :: taspi,twzz,twtwz,fozz,sixzz,fotwz,twtwtw,etzz,sixtwz,&
                        fofoz,fotwtw

  ! setup variables

    taspi  =   2.0_wp * alpha / sqrpi ;             twzz  =  -2.0_wp * alpha**2 * taspi / 3.0_wp
    twtwz  =   4.0_wp * alpha**4 * taspi / 5.0_wp ; fozz  =  12.0_wp * alpha**4 * taspi / 5.0_wp
    sixzz  =-120.0_wp * alpha**6 * taspi / 7.0_wp ; fotwz = -24.0_wp * alpha**6 * taspi / 7.0_wp
    twtwtw =  -8.0_wp * alpha**6 * taspi / 7.0_wp ; etzz  = 560.0_wp * alpha**8 * taspi / 3.0_wp
    sixtwz =  80.0_wp * alpha**8 * taspi / 3.0_wp ; fofoz =  16.0_wp * alpha**9 * taspi
    fotwtw =  16.0_wp * alpha**9 * taspi / 3.0_wp

    d1 = 0.0_wp

    d1(0,0,0) = taspi

    d1(2,0,0) = twzz ;   d1(0,2,0) = twzz;    d1(0,0,2)=twzz

    d1(4,0,0) = fozz ;   d1(2,2,0) = twtwz;   d1(2,0,2)=twtwz ;    d1(0,4,0) = fozz
    d1(0,2,2) = twtwz ;  d1(0,0,4) = fozz

    d1(6,0,0) = sixzz ;  d1(4,2,0) = fotwz
    d1(4,0,2) = fotwz ;  d1(2,4,0) = fotwz ;  d1(2,2,2) = twtwtw ; d1(2,0,4) = fotwz
    d1(0,6,0) = sixzz ;  d1(0,4,2) = fotwz ;  d1(0,2,4) = fotwz ;  d1(0,0,6) = sixzz

    d1(8,0,0) = etzz ;   d1(6,2,0) = sixtwz ; d1(6,0,2) = sixtwz ; d1(4,4,0) = fofoz
    d1(4,2,2) = fotwtw ; d1(4,0,4) = fofoz ;  d1(2,6,0) = sixtwz ; d1(2,0,6) = sixtwz
    d1(2,4,2) = fotwtw ; d1(2,2,4) = fotwtw ; d1(0,8,0) = etzz;    d1(0,6,2) = sixtwz
    d1(0,4,4) = fofoz ;  d1(0,2,6) = sixtwz ; d1(0,0,8) = etzz

  End subroutine limit_erfr_deriv

  Subroutine explicit_ewald_real_loops         &
             (lbl,torderlim,k1,k2,k3,alpha,d1, &
             imp,impx,impy,impz,tix,tiy,tiz,   &
             jmp,jmpx,jmpy,jmpz,tjx,tjy,tjz,   &
             enempol,fx,fy,fz)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for computing the energy, forces and torques
  ! for multipolar interactions.  An explicit formula is used for multipolar
  ! orders up to 4 to avoid the 6-deep do loop.
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent( In    ) :: lbl,torderlim,k1,k2,k3

    Real( Kind = wp), Intent( In    ) :: alpha,d1(lbl:torderlim,lbl:torderlim,lbl:torderlim)

    Real( Kind = wp), Intent( In    ) :: imp(1:mximpl)
    Real( Kind = wp), Intent( In    ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp), Intent( InOut ) :: tix,tiy,tiz

    Real( Kind = wp), Intent( In    ) :: jmp
    Real( Kind = wp), Intent( In    ) :: jmpx,jmpy,jmpz
    Real( Kind = wp), Intent( InOut ) :: tjx,tjy,tjz

    Real( Kind = wp), Intent( InOut ) :: enempol,fx,fy,fz

  ! Local variables

    Integer          :: ii,n,s1,s2,s3,ks1,ks2,ks3,ks11,ks21,ks31
    Real( Kind = wp) :: tmp,tmpi,tmpj,t1,alphan,td

    If (lbl /= 0 .and. lbl /= -2) STOP ! Allowed Lower_Bound_Limit

    If (mxompl >= 0) Then

  !=======================================
       s1=0; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan*jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl >= 1) Then

  !=======================================
       s1=1; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl >= 2) Then

  !=======================================
       s1=2; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=2; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl >= 3) Then

  !=======================================
       s1=3; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=3; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force
       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=3
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=2; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=0; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=1; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=2; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl == 4) Then

  !=======================================
       s1=4; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=4; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=4
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=3; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=3; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=2; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=0; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=1; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=3; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=2; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=1; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=0; s3=3
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=2; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=3; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=3
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)

       tmp=alphan*td

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t1 = alphan * jmp*imp(ii)

  ! energy

       enempol = enempol + t1*td

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

  End Subroutine explicit_ewald_real_loops

  Subroutine explicit_fscp_rfp_loops          &
             (torderlim,k1,k2,k3,alpha,d1,a1, &
             imp,impx,impy,impz,tix,tiy,tiz,  &
             jmp,jmpx,jmpy,jmpz,tjx,tjy,tjz,  &
             enempol,fx,fy,fz)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for computing the energy, forces and torques
  ! for multipolar interactions.  An explicit formula is used for multipolar
  ! orders up to 4 to avoid the 6-deep do loop.
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng july 2014
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent( In    ) :: torderlim,k1,k2,k3

    Real( Kind = wp), Intent( In    ) :: alpha, d1(-2:torderlim,-2:torderlim,-2:torderlim), &
                                                a1(-2:torderlim,-2:torderlim,-2:torderlim)

    Real( Kind = wp), Intent( In    ) :: imp(1:mximpl)
    Real( Kind = wp), Intent( In    ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp), Intent( InOut ) :: tix,tiy,tiz

    Real( Kind = wp), Intent( In    ) :: jmp
    Real( Kind = wp), Intent( In    ) :: jmpx,jmpy,jmpz
    Real( Kind = wp), Intent( InOut ) :: tjx,tjy,tjz

    Real( Kind = wp), Intent( InOut ) :: enempol,fx,fy,fz

  ! Local variables

    Integer          :: ii,n,s1,s2,s3,ks1,ks2,ks3,ks11,ks21,ks31
    Real( Kind = wp) :: tmp,tmpi,tmpj,t1,t2,alphan,td,ta

    If (mxompl >= 0) Then

  !=======================================
       s1=0; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td+t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl >= 1) Then

  !=======================================
       s1=1; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl >= 2) Then

  !=======================================
       s1=2; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=2; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl >= 3) Then

  !=======================================
       s1=3; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=3; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=3
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=2; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=0; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=1; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=2; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj =-imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

    If (mxompl == 4) Then

  !=======================================
       s1=4; s2=0; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=4; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=0; s3=4
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=3; s2=1; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=3; s2=0; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=2; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=0; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=2; s2=1; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=3; s3=0
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=2; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=1; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=1; s2=0; s3=3
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=2; s3=2
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=3; s3=1
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

  !=======================================
       s1=0; s2=1; s3=3
  !=======================================

       ks3=k3+s3; ks31=ks3+1
       ks2=k2+s2; ks21=ks2+1
       ks1=k1+s1; ks11=ks1+1

       n=ks1+ks2+ks3
       alphan=alpha**n

       ii=mplmap(s1,s2,s3)
       td=d1(ks1,ks2,ks3)
       ta=a1(ks1,ks2,ks3)

       tmp=alphan*td+a1(ks1,ks2,ks3)

       tmpi = jmp     * tmp
       tmpj = imp(ii) * tmp

       t2 = jmp*imp(ii)
       t1 = alphan*t1

  ! energy

       enempol = enempol + t1*td + t2*ta

  ! force

       t1 = t1*alpha

       fx = fx - t1*d1(ks11,ks2,ks3) - t2*a1(ks11,ks2,ks3)
       fy = fy - t1*d1(ks1,ks21,ks3) - t2*a1(ks1,ks21,ks3)
       fz = fz - t1*d1(ks1,ks2,ks31) - t2*a1(ks1,ks2,ks31)

  ! torque on iatm

       tix = tix + impx(ii)*tmpi
       tiy = tiy + impy(ii)*tmpi
       tiz = tiz + impz(ii)*tmpi

  ! torque on jatm

       tjx = tjx + jmpx*tmpj
       tjy = tjy + jmpy*tmpj
       tjz = tjz + jmpz*tmpj

    End If

  End Subroutine explicit_fscp_rfp_loops

  Subroutine explicit_spme_loops &
             (flag,rcell,bdx,bdy,bdz,imp,impx,impy,impz, &
             dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to speed up inner loop of spme force computations
  ! for multipolar interactions for orders less than 5
  !
  ! Note: reference to an spme_container routine
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng july 2014
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,          Intent( In    ) :: flag
    Real( Kind = wp), Intent( In    ) :: rcell(1:9)
    Real( Kind = wp), Intent( In    ) :: imp(1:mximpl)
    Real( Kind = wp), Intent( In    ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
    Real( Kind = wp), Intent( In    ) :: bdx(0:mxspl),bdy(0:mxspl),bdz(0:mxspl)
    Real( Kind = wp), Intent(   Out ) :: dtp,tq1,tq2,tq3,dt1,dt2,dt3,td1,td2,td3

  ! Local variables

    Integer          :: mm,s1,s2,s3
    Real( Kind = wp) :: Dtpbsp,tmp,impt,tid

  ! Initialise per loop contributions

    dtp=0.0_wp
    tq1=0.0_wp ; tq2=0.0_wp ; tq3=0.0_wp ! pole
    dt1=0.0_wp ; dt2=0.0_wp ; dt3=0.0_wp ! force
    td1=0.0_wp ; td2=0.0_wp ; td3=0.0_wp ! torque

    If (mxompl >= 0) Then

  !===================================================
       s1=0; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

    End If

    If (mxompl >= 1) Then

  !===================================================
       s1=1; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq1 = tq1 + tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq2 = tq2 + tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq3 = tq3 + tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

    End If

    If (mxompl >= 2) Then

  !===================================================
       s1=2; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq1 = tq1 + 2.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=2; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq2 = tq2 + 2.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=0; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq3 = tq3 + 2.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq2 = tq2 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=1; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq2 = tq2 + tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

    End If

    If (mxompl >= 3) Then

  !===================================================
       s1=3; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq1 = tq1 + 3.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=3; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq2 = tq2 + 3.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=0; s3=3
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq3 = tq3 + 3.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=2; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 2.0_wp*tid
          tq2 = tq2 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=2; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 2.0_wp*tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=2; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq2 = tq2 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=0; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq3 = tq3 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=1; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq2 = tq2 + tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=2; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq2 = tq2 + 2.0_wp*tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=1; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq2 = tq2 + tid
          tq3 = tq3 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

    End If

    If (mxompl == 4) Then

  !===================================================
       s1=4; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq1 = tq1 + 4.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=4; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq2 = tq2 + 4.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=0; s3=4
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) tq3 = tq3 + 4.0_wp*tid

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=3; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 3.0_wp*tid
          tq2 = tq2 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=3; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 3.0_wp*tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=2; s2=2; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 2.0_wp*tid
          tq2 = tq2 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=2; s2=0; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 2.0_wp*tid
          tq3 = tq3 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=2; s2=1; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + 2.0_wp*tid
          tq2 = tq2 + tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=3; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq2 = tq2 + 3.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=2; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq2 = tq2 + 2.0_wp*tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=1; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq2 = tq2 + tid
          tq3 = tq3 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=1; s2=0; s3=3
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq1 = tq1 + tid
          tq3 = tq3 + 3.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=2; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq2 = tq2 + 2.0_wp*tid
          tq3 = tq3 + 2.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=3; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq2 = tq2 + 3.0_wp*tid
          tq3 = tq3 + tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

  !===================================================
       s1=0; s2=1; s3=3
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

       If (flag == 0) Then

          tq2 = tq2 + tid
          tq3 = tq3 + 3.0_wp*tid

       End If

       If (flag == 1) Then

          dt1 = dt1 + impt*Dtpbsp(s1+1,s2,s3,rcell,bdx,bdy,bdz)
          dt2 = dt2 + impt*Dtpbsp(s1,s2+1,s3,rcell,bdx,bdy,bdz)
          dt3 = dt3 + impt*Dtpbsp(s1,s2,s3+1,rcell,bdx,bdy,bdz)

          td1 = td1 + impx(mm)*tmp
          td2 = td2 + impy(mm)*tmp
          td3 = td3 + impz(mm)*tmp

       End If

    End If

  End Subroutine explicit_spme_loops

  Subroutine explicit_spme_loop_s(rcell,bdx,bdy,bdz,imp,dtp)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine to speed up inner loop of spme mfield computations
  ! for multipolar interactions for orders less than 5
  !
  ! Note: reference to an spme_container routine
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp), Intent( In    ) :: rcell(1:9)
    Real( Kind = wp), Intent( In    ) :: bdx(0:mxspl),bdy(0:mxspl),bdz(0:mxspl)
    Real( Kind = wp), Intent( In    ) :: imp(1:mximpl)
    Real( Kind = wp), Intent(   Out ) :: dtp

  ! Local variables

    Integer          :: mm,s1,s2,s3
    Real( Kind = wp) :: Dtpbsp,tmp,impt,tid

    If (mxompl >= 0) Then

  !===================================================
       s1=0; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

    End If

    If (mxompl >= 1) Then

  !===================================================
       s1=1; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

    End If

    If (mxompl >= 2) Then

  !===================================================
       s1=2; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=2; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=0; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=1; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

    End If

    If (mxompl >= 3) Then

  !===================================================
       s1=3; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=3; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=0; s3=3
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=2; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=2; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=2; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=0; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=1; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=2; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=1; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

    End If

    If (mxompl == 4) Then

  !===================================================
       s1=4; s2=0; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=4; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=0; s3=4
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=3; s2=1; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=3; s2=0; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=2; s2=2; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=2; s2=0; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=2; s2=1; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=3; s3=0
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=2; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=1; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=1; s2=0; s3=3
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=2; s3=2
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=3; s3=1
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

  !===================================================
       s1=0; s2=1; s3=3
  !===================================================

       mm  = mplmap(s1,s2,s3)
       impt= imp(mm)

       tmp = Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

       tid = impt*tmp
       dtp = dtp + tid

    End If

  End Subroutine explicit_spme_loop_s

  Subroutine rotate_mpoles(iatm,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for rotating multipoles from local frame to global
  ! frame
  !
  ! Note: reference to numeric_container routines
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng april 2014
  ! amended   - i.t.todorov december 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent( In    )  :: iatm
    Type( comms_type ), Intent( In    ) :: comm

  ! Local variables

    Integer           :: i,j,k,m,n,numnbh,idi,fail,local_index
    Real( Kind = wp ) :: a(9),temp(1:mximpl),mpole(1:mximpl) !(1:(3**(mxompl+1)-1)/2),temp(1:mximpl)
    Real( Kind = wp ) :: ai,ai3,ai6,aj3,aj6,ak3,ak6,am3,am6,                &
                         tmp,t1,t2,t3,t4,t5,t6,p1x,p1y,p1z,p2x,p2y,p2z,dpu, &
                         rrp1,rrp2,rrp3,rrp12,rrp22,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10

    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rsqt

  ! Now on to permanent multipoles

    If (mplflg(iatm) == 1) Return ! multipole is already rotated

    mplflg(iatm) = 1

    idi = ltg(iatm)

    If (mxompl == 0) Then

       mplgfr(:,iatm) = mpllfr(:,lsite(iatm))

       Return

    End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find the rotation matrix a(3,3) = a(9)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    If (mprotm(iatm)%flag == 0) Then

       numnbh = ltpatm(0,iatm)

       If (numnbh == 0) Then

          mplgfr(:,iatm) = mpllfr(:,lsite(iatm))

          Return ! Atom has no bonds => global frame = local frame

       End If

       mprotm(iatm)%flag = 1

       If (numnbh == 1) Then

          j    = local_index(ltpatm(1,iatm),nlast,lsi,lsa)

          p1x  = xxx(j) - xxx(iatm)
          p1y  = yyy(j) - yyy(iatm)
          p1z  = zzz(j) - zzz(iatm)

          Call images_s(imcon,cell,p1x,p1y,p1z)

          rrp1 = Sqrt(p1x**2+p1y**2+p1z**2)

          a(1) = p1x/rrp1
          a(4) = p1y/rrp1
          a(7) = p1z/rrp1
          a(2) =-a(4)
          a(5) = a(1)
          a(8) = a(7)
          a(3) = a(4)*a(8)-a(5)*a(7)
          a(6) = a(2)*a(7)-a(1)*a(8)
          a(9) = a(1)*a(5)-a(5)*a(4)

          rrp3 = Sqrt(a(3)*a(3)+a(6)*a(6)+a(9)*a(9))

          a(3) = a(3)/rrp3
          a(6) = a(6)/rrp3
          a(9) = a(9)/rrp3

          mprotm(iatm)%p1 = (/ p1x, p1y, p1z /)
          mprotm(iatm)%p2 = (/-p1y, p1x, p1z /)

          mprotm(iatm)%mbnd(1) = j
          mprotm(iatm)%mbnd(2) = 0

  ! Is there no torque when atom has only one bond?

          mprotm(iatm)%flag = 1

       End If

       If (numnbh == 2) Then

          j    = local_index(ltpatm(1,iatm),nlast,lsi,lsa)
          k    = local_index(ltpatm(2,iatm),nlast,lsi,lsa)

          p1x  = xxx(j) - xxx(iatm)
          p1y  = yyy(j) - yyy(iatm)
          p1z  = zzz(j) - zzz(iatm)

          Call images_s(imcon,cell,p1x,p1y,p1z)

          p2x  = xxx(k) - xxx(iatm)
          p2y  = yyy(k) - yyy(iatm)
          p2z  = zzz(k) - zzz(iatm)

          Call images_s(imcon,cell,p2x,p2y,p2z)

          rrp12 = p1x**2+p1y**2+p1z**2
          rrp22 = p2x**2+p2y**2+p2z**2

          If ((p1x*p2x+p1y*p2y+p1z*p2z)**2 < rrp12*rrp22 - 1.0e-6_wp) Then

  ! find the bisector of the angle between the two local vectors

             p1x  = p1x + p2x
             p1y  = p1y + p2y
             p1z  = p1z + p2z

             rrp1 = Sqrt(p1x**2+p1y**2+p1z**2)

             mprotm(iatm)%p1 = (/ p1x, p1y, p1z /)
             mprotm(iatm)%p2 = (/ p2x, p2y, p2z /)

             mprotm(iatm)%mbnd(1) = j
             mprotm(iatm)%mbnd(2) = k

             a(1) = p1x/rrp1
             a(4) = p1y/rrp1
             a(7) = p1z/rrp1

             dpu  = p2x*a(1)+p2y*a(4)+p2z*a(7)

             p2x  = p2x - dpu*a(1)
             p2y  = p2y - dpu*a(4)
             p2z  = p2z - dpu*a(7)

             rrp2 = Sqrt(p2x**2+p2y**2+p2z**2)

             a(2) = p2x/rrp2
             a(5) = p2y/rrp2
             a(8) = p2z/rrp2

             a(3) = a(4)*a(8)-a(5)*a(7)
             a(6) = a(2)*a(7)-a(1)*a(8)
             a(9) = a(1)*a(5)-a(5)*a(4)

             rrp3 = Sqrt(a(3)**2+a(6)**2+a(9)**2)

             a(3) = a(3)/rrp3
             a(6) = a(6)/rrp3
             a(9) = a(9)/rrp3

          Else

             rrp1 = Sqrt(rrp12)

             a(1) = p1x/rrp1
             a(4) = p1y/rrp1
             a(7) = p1z/rrp1
             a(2) = -a(4)
             a(5) = a(1)
             a(8) = a(7)
             a(3) = a(4)*a(8)-a(5)*a(7)
             a(6) = a(2)*a(7)-a(1)*a(8)
             a(9) = a(1)*a(5)-a(5)*a(4)

             rrp3 = Sqrt(a(3)**2+a(6)**2+a(9)**2)

             a(3) = a(3)/rrp3
             a(6) = a(6)/rrp3
             a(9) = a(9)/rrp3

             mprotm(iatm)%p1 = (/ p1x, p1y, p1z /)
             mprotm(iatm)%p2 = (/-p1y, p1x, p1z /)

             mprotm(iatm)%mbnd(1) = j
             mprotm(iatm)%mbnd(2) = k

  ! Is there no torque when the bond angle is pi?

             mprotm(iatm)%flag = 1

          End If

       End If

       If (numnbh > 2) Then

          fail = 0
          Allocate (xxt(1:numnbh),yyt(1:numnbh),zzt(1:numnbh),rsqt(1:numnbh), Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'allocation failure in rotate_mpoles, node: ', comm%idnode
             Call error(0)
          End If

          Do i = 1,numnbh

             j       = local_index(ltpatm(i,iatm),nlast,lsi,lsa)

             xxt(i)  = xxx(j) - xxx(iatm)
             yyt(i)  = yyy(j) - yyy(iatm)
             zzt(i)  = zzz(j) - zzz(iatm)

             Call images_s(imcon,cell,xxt(i),yyt(i),zzt(i))

             rsqt(i) = xxt(i)**2+yyt(i)**2+zzt(i)**2

          End Do

          k = Minloc(rsqt,Dim=1)
          rrp12 = rsqt(k)
          rsqt(k) = 0.0_wp

          m = Minloc(rsqt,Dim=1,Mask=(rsqt > 0.0_wp))
          rrp22 = rsqt(m)

          Do While ((xxt(k)*xxt(m)+yyt(k)*yyt(m)+zzt(k)*zzt(m))**2 >= rrp12*rrp22 - 1.0e-6_wp)

             rsqt(m) = 0.0_wp
             m = Minloc(rsqt,Dim=1,Mask=(rsqt > 0.0_wp))
             rrp22 = rsqt(m)

          End Do

          mprotm(iatm)%mbnd(1) = k
          mprotm(iatm)%mbnd(2) = m

          p1x  = xxt(k) + xxt(m)
          p1y  = yyt(k) + yyt(m)
          p1z  = zzt(k) + zzt(m)

          rrp1 = Sqrt(p1x**2+p1y**2+p1z**2)

          mprotm(iatm)%p1 = (/ p1x,    p1y,    p1z    /)
          mprotm(iatm)%p2 = (/ xxt(m), yyt(m), zzt(m) /)

          a(1) = p1x/rrp1
          a(4) = p1y/rrp1
          a(7) = p1z/rrp1

          dpu  = xxt(m)*a(1)+yyt(m)*a(4)+zzt(m)*a(7)

          p2x  = xxt(m) - dpu*a(1)
          p2y  = yyt(m) - dpu*a(4)
          p2z  = zzt(m) - dpu*a(7)

          rrp2 = Sqrt(p2x**2+p2y**2+p2z**2)

          a(2) = p2x/rrp2
          a(5) = p2y/rrp2
          a(8) = p2z/rrp2

          a(3) = a(4)*a(8)-a(5)*a(7)
          a(6) = a(2)*a(7)-a(1)*a(8)
          a(9) = a(1)*a(5)-a(5)*a(4)

          rrp3 = Sqrt(a(3)**2+a(6)**2+a(9)**2)

          a(3) = a(3)/rrp3
          a(6) = a(6)/rrp3
          a(9) = a(9)/rrp3

          Deallocate (xxt,yyt,zzt,rsqt, Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'deallocation failure in rotate_mpoles, node: ', comm%idnode
             Call error(0)
          End If

       End If

       mprotm(iatm)%mtrxa = a

    Else

       a = mprotm(iatm)%mtrxa

    End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now on to rotate multipoles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! get the local multipoles based on atom type

    mpole(:)=mpllfr(:,lsite(iatm))

    temp=0.0_wp

  ! rotate monopole

    temp(1)=mpole(1)

    If (mxompl >= 1) Then

  ! rotate dipoles

       temp(2)=a(1)*mpole(2) + a(2)*mpole(3) + a(3)*mpole(4)
       temp(3)=a(4)*mpole(2) + a(5)*mpole(3) + a(6)*mpole(4)
       temp(4)=a(7)*mpole(2) + a(8)*mpole(3) + a(9)*mpole(4)

    End If

    If (mxompl >= 2) Then

  ! rotate quadrupoles

       k=4

       Do i = 1, 3
          ai = a(i) ; ai3=a(i+3) ; ai6=a(i+6)

          Do j = 1, 3
             k = k + 1 ; aj3=a(j+3 ); aj6=a(j+6)

             tmp = mpole(mplltg(k))

             temp(5)  = temp(5)  + ai  * a(j) * tmp !Q_xx
             temp(6)  = temp(6)  + ai  * aj3  * tmp !Q_xy
             temp(7)  = temp(7)  + ai  * aj6  * tmp !Q_xz

             temp(8)  = temp(8)  + ai3 * aj3  * tmp !Q_yy
             temp(9)  = temp(9)  + ai3 * aj6  * tmp !Q_yz

             temp(10) = temp(10) + ai6 * aj6  * tmp !Q_zz
          End Do
       End Do

  ! Q_xy, Q_xz and Q_yz have to be doubled to account for
  ! Q_yx, Q_zx and Q_yz

       temp(6) = 2.0_wp*temp(6)
       temp(7) = 2.0_wp*temp(7)
       temp(9) = 2.0_wp*temp(9)

    End If

    If (mxompl >= 3) Then

  ! rotate octupoles

       m=k

       Do i = 1, 3
          ai = a(i) ; ai3 = a(i+3) ; ai6 = a(i+6)

          Do j = 1, 3
             aj3 = a(j+3) ; aj6 = a(j+6)

             t1 = ai  * a(j) ; t2 = ai  * aj3 ; t3 = ai  * aj6
             t4 = ai3 * aj3  ; t5 = ai3 * aj6
             t6 = ai6 * aj6

             Do k = 1, 3
                m = m + 1; ak3=a(k+3); ak6=a(k+6)

                tmp = mpole(mplltg(m))

                temp(11) = temp(11) + t1 * a(k) * tmp !O_xxx
                temp(12) = temp(12) + t1 * ak3  * tmp !O_xxy
                temp(13) = temp(13) + t1 * ak6  * tmp !O_xxz

                temp(14) = temp(14) + t2 * ak3  * tmp !O_xyy
                temp(15) = temp(15) + t2 * ak6  * tmp !O_xyz

                temp(16) = temp(16) + t3 * ak6  * tmp !O_xzz

                temp(17) = temp(17) + t4 * ak3  * tmp !O_yyy
                temp(18) = temp(18) + t4 * ak6  * tmp !O_yyz

                temp(19) = temp(19) + t5 * ak6  * tmp !O_yzz

                temp(20) = temp(20) + t6 * ak6  * tmp !O_zzz
             End Do
          End Do
       End Do

  ! O_xxy,O_xxz,O_xyy,O_xzz,O_yyz,O_yzz have to be tripled
  ! to account for permutations.  O_xyz has to be multiplied by six

       temp(12) = 3.0_wp*temp(12) ; temp(13) = 3.0_wp*temp(13)
       temp(14) = 3.0_wp*temp(14) ; temp(15) = 6.0_wp*temp(15)
       temp(16) = 3.0_wp*temp(16) ; temp(18) = 3.0_wp*temp(18)
       temp(19) = 3.0_wp*temp(19)

    End If

    If (mxompl >= 4) Then

  ! rotate hexadecapoles

       n = m

       Do i = 1, 3
          ai = a(i) ; ai3=a(i+3) ; ai6=a(i+6)

          Do j = 1, 3
             aj3 = a(j+3) ; aj6 = a(j+6)

             t1 = ai  * a(j) ; t2 = ai  * aj3 ; t3 = ai  * aj6
             t4 = ai3 * aj3  ; t5 = ai3 * aj6
             t6 = ai6 * aj6

             Do k = 1, 3
                ak3=a(k+3) ; ak6= a(k+6)
                s1  = t1*a(k); s2 = t1*ak3; s3 = t1*ak6
                s4  = t2*ak3 ; s5 = t2*ak6
                s6  = t3*ak6
                s7  = t4*ak3 ; s8 = t4*ak6
                s9  = t5*ak6
                s10 = t6*ak6

                Do m = 1, 3
                   n = n + 1 ; am3 = a(m+3); am6 = a(m+6)

                   tmp = mpole(mplltg(n))

                   temp(21) = temp(21) + s1  * a(m) * tmp !H_xxxx
                   temp(22) = temp(22) + s1  * am3  * tmp !H_xxxy
                   temp(23) = temp(23) + s1  * am6  * tmp !H_xxxz

                   temp(24) = temp(24) + s2  * am3  * tmp !H_xxyy
                   temp(25) = temp(25) + s2  * am6  * tmp !H_xxyz

                   temp(26) = temp(26) + s3  * am6  * tmp !H_xxzz

                   temp(27) = temp(27) + s4  * am3  * tmp !H_xyyy
                   temp(28) = temp(28) + s4  * am6  * tmp !H_xyyz

                   temp(29) = temp(29) + s5  * am6  * tmp !H_xyzz

                   temp(30) = temp(30) + s6  * am6  * tmp !H_xzzz

                   temp(31) = temp(31) + s7  * am3  * tmp !H_yyyy
                   temp(32) = temp(32) + s7  * am6  * tmp !H_yyyz

                   temp(33) = temp(33) + s8  * am6  * tmp !H_yyzz

                   temp(34) = temp(34) + s9  * am6  * tmp !H_yzzz

                   temp(35) = temp(35) + s10 * am6  * tmp !H_zzzz
                End Do
             End Do
          End Do
       End Do

  ! Account for permutations of H_xxxy, H_xxxz, H_xxyy, H_xxyz,
  ! H_xxzz, H_xyyy, H_xyyz, H_xyzz, H_xzzz, H_yyyz, H_yyzz, H_yzzz

       temp(22) =  4.0_wp * temp(22) ; temp(23) =  4.0_wp * temp(23)
       temp(24) =  6.0_wp * temp(24) ; temp(25) = 12.0_wp * temp(25)
       temp(26) =  6.0_wp * temp(26) ; temp(27) =  4.0_wp * temp(27)
       temp(28) = 12.0_wp * temp(28) ; temp(29) = 12.0_wp * temp(29)
       temp(30) =  4.0_wp * temp(30) ; temp(32) =  4.0_wp * temp(32)
       temp(33) =  6.0_wp * temp(33) ; temp(34) =  4.0_wp * temp(34)

    End If

    mplgfr(:,iatm) = temp

  ! Find infinitesimal rotation of global multipoles
  !==================================================

    Call infinitesimal_rotation(iatm)

  !===================================================

  End Subroutine rotate_mpoles

  Subroutine infinitesimal_rotation(iatm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for infinitesimal counter-clockwise rotation of
  ! global multipoles needed to find the torque at each atomic site
  !
  ! Reference : Sagui, Pedersen, Darden, J. Chem. Phys. 120, 73 (2004)
  !             doi: 10.1063/1.1630791
  !
  ! Note : The infinitesimal rotations given in the above paper are actually
  !        clockwise rotations.  The torque, T = -M * Phi, i.e. the negative
  !        of the scalar product of the multipole and the potential, where M
  !        is the counter-clockwise infinitesimal rotation of the multipole.
  !        However, M' = -M, where M' is the clockwise rotation of the
  !        multipoles, hence T = M' * Phi
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov february 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent( In    )  :: iatm

  ! Local variables

    Integer          :: i,j,k,m,n,ll
    Real( Kind = wp) :: uni(1:9),irotdir(1:9,1:3),irot(1:9), &
                        mpole(1:mximpl),temp(1:mximpl)
    Real( Kind = wp) :: tmp,ai,ai3,ai6,aj3,aj6,ak3,ak6,am3,am6,           &
                        bi,bi3,bi6,bj3,bj6,bk3,bk6,bm3,bm6,               &
                        t1,t2,t3,t4,t5,t6,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, &
                        r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,p1,p2,p3,p4,p5,p6

  ! Identity matrix

    uni(1:9)       = (/ 1.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 1.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 1.0_wp /)

  ! Infinitesimal rotation matrices

    irotdir(1:9,1) = (/ 0.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 1.0_wp , 0.0_wp ,-1.0_wp , 0.0_wp /)
    irotdir(1:9,2) = (/ 0.0_wp , 0.0_wp ,-1.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 1.0_wp , 0.0_wp , 0.0_wp /)
    irotdir(1:9,3) = (/ 0.0_wp , 1.0_wp , 0.0_wp ,-1.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 0.0_wp /)

    mpole = mplgfr(:,iatm)

    Do ll=1,3

       irot = irotdir(:,ll)

       temp=0.0_wp

  ! infinitesimal rotation has no effect on charge/monopole => temp(1)=0.0_wp

       If (mxompl >= 1) Then

  ! rotate dipoles

          temp(2) = irot(1)*mpole(2) + irot(2)*mpole(3) + irot(3)*mpole(4)
          temp(3) = irot(4)*mpole(2) + irot(5)*mpole(3) + irot(6)*mpole(4)
          temp(4) = irot(7)*mpole(2) + irot(8)*mpole(3) + irot(9)*mpole(4)

       End If

       If (mxompl >= 2) Then

  ! rotate quadrupoles

          k=4

          Do i = 1, 3
             ai =  uni(i) ; ai3 =  uni(i+3) ; ai6 =  uni(i+6)
             bi = irot(i) ; bi3 = irot(i+3) ; bi6 = irot(i+6)

             Do j = 1, 3
                aj3 = irot(j+3) ; aj6 = irot(j+6)
                bj3 =  uni(j+3) ; bj6 =  uni(j+6)

                k = k + 1

                tmp = mpole(mplltg(k))

                temp(5)  = temp(5)  + (ai  * irot(j) + bi  * uni(j)) * tmp !Q_xx
                temp(6)  = temp(6)  + (ai  * aj3     + bi  * bj3   ) * tmp !Q_xy
                temp(7)  = temp(7)  + (ai  * aj6     + bi  * bj6   ) * tmp !Q_xz

                temp(8)  = temp(8)  + (ai3 * aj3     + bi3 * bj3   ) * tmp !Q_yy
                temp(9)  = temp(9)  + (ai3 * aj6     + bi3 * bj6   ) * tmp !Q_yz

                temp(10) = temp(10) + (ai6 * aj6     + bi6 * bj6   ) * tmp !Q_zz
             End Do
          End Do

       End If

       If (mxompl >= 3) Then

  ! rotate octupoles

          m=k

          Do i = 1, 3
             ai =  uni(i) ; ai3 =  uni(i+3) ; ai6 =  uni(i+6)
             bi = irot(i) ; bi3 = irot(i+3) ; bi6 = irot(i+6)

             Do j = 1, 3
                aj3 =  uni(j+3) ; aj6 =  uni(j+6)
                bj3 = irot(j+3) ; bj6 = irot(j+6)

                t1 = ai  * uni(j) ; t2 = ai  * aj3 ; t3 = ai  * aj6
                t4 = ai3 * aj3    ; t5 = ai3 * aj6
                t6 = ai6 * aj6

                s1 = ai  * irot(j) + bi  * uni(j) ; s2 = ai  * bj3 + bi  * aj3 ; s3 = ai * bj6 + bi * aj6
                s4 = ai3 * bj3     + bi3 * aj3    ; s5 = ai3 * bj6 + bi3 * aj6
                s6 = ai6 * bj6     + bi6 * aj6

                Do k = 1, 3
                   ak3 = irot(k+3) ; ak6 = irot(k+6)
                   bk3 = uni(k+3)    ; bk6 = uni(k+6)

                   m = m + 1

                   tmp = mpole(mplltg(m))

                   temp(11) = temp(11) + (t1 * irot(k) + s1 * uni(k)) * tmp !O_xxx
                   temp(12) = temp(12) + (t1 * ak3     + s1 * bk3   ) * tmp !O_xxy
                   temp(13) = temp(13) + (t1 * ak6     + s1 * bk6   ) * tmp !O_xxz
                   temp(14) = temp(14) + (t2 * ak3     + s2 * bk3   ) * tmp !O_xyy
                   temp(15) = temp(15) + (t2 * ak6     + s2 * bk6   ) * tmp !O_xyz

                   temp(16) = temp(16) + (t3 * ak6     + s3 * bk6   ) * tmp !O_xzz

                   temp(17) = temp(17) + (t4 * ak3     + s4 * bk3   ) * tmp !O_yyy
                   temp(18) = temp(18) + (t4 * ak6     + s4 * bk6   ) * tmp !O_yyz

                   temp(19) = temp(19) + (t5 * ak6     + s5 * bk6   ) * tmp !O_yzz

                   temp(20) = temp(20) + (t6 * ak6     + s6 * bk6   ) * tmp !O_zzz
                End Do
             End Do
          End Do

       End If

       If (mxompl >= 4) Then

  ! rotate hexadecapoles

          n = m

          Do i = 1, 3
             ai =  uni(i) ; ai3 =  uni(i+3) ; ai6 =  uni(i+6)
             bi = irot(i) ; bi3 = irot(i+3) ; bi6 = irot(i+6)

             Do j = 1, 3
                aj3 =  uni(j+3) ; aj6 =  uni(j+6)
                bj3 = irot(j+3) ; bj6 = irot(j+6)

                t1 = ai  * uni(j) ; t2 = ai  * aj3 ; t3 = ai * aj6
                t4 = ai3 * aj3    ; t5 = ai3 * aj6
                t6 = ai6 * aj6

                p1 = ai  * irot(j) + bi  * uni(j) ; p2 = ai  * bj3 + bi  * aj3 ; p3 = ai * bj6 + bi * aj6
                p4 = ai3 * bj3     + bi3 * aj3    ; p5 = ai3 * bj6 + bi3 * aj6
                p6 = ai6 * bj6     + bi6 * aj6

                Do k = 1, 3
                   ak3 =  uni(k+3) ; ak6 =  uni(k+6)
                   bk3 = irot(k+3) ; bk6 = irot(k+6)

                   s1  = t1 * uni(k) ; s2 = t1 * ak3 ; s3 = t1 * ak6
                   s4  = t2 * ak3    ; s5 = t2 * ak6
                   s6  = t3 * ak6
                   s7  = t4 * ak3    ; s8 = t4 * ak6
                   s9  = t5 * ak6
                   s10 = t6 * ak6

                   r1  = t1 * irot(k) + p1 * uni(k) ; r2 = t1 * bk3 + p1 * ak3 ; r3 = t1 * bk6 + p1 * ak6
                   r4  = t2 * bk3     + p2 * ak3    ; r5 = t2 * bk6 + p2 * ak6
                   r6  = t3 * bk6     + p3 * ak6
                   r7  = t4 * bk3     + p4 * ak3    ; r8 = t4 * bk6 + p4 * ak6
                   r9  = t5 * bk6     + p5 * ak6
                   r10 = t6 * bk6     + p6 * ak6

                   Do m = 1, 3
                      am3 = irot(m+3) ; am6 = irot(m+6)
                      bm3 =  uni(m+3) ; bm6 =  uni(m+6)

                      n = n + 1

                      tmp = mpole(mplltg(n))

                      temp(21) = temp(21) + (s1  * irot(m) + r1  * uni(m)) * tmp !H_xxxx
                      temp(22) = temp(22) + (s1  * am3     + r1  * bm3   ) * tmp !H_xxxy
                      temp(23) = temp(23) + (s1  * am6     + r1  * bm6   ) * tmp !H_xxxz

                      temp(24) = temp(24) + (s2  * am3     + r2  * bm3   ) * tmp !H_xxyy
                      temp(25) = temp(25) + (s2  * am6     + r2  * bm6   ) * tmp !H_xxyz

                      temp(26) = temp(26) + (s3  * am6     + r3  * bm6   ) * tmp !H_xxzz

                      temp(27) = temp(27) + (s4  * am3     + r4  * bm3   ) * tmp !H_xyyy
                      temp(28) = temp(28) + (s4  * am6     + r4  * bm6   ) * tmp !H_xyyz

                      temp(29) = temp(29) + (s5  * am6     + r5  * bm6   ) * tmp !H_xyzz

                      temp(30) = temp(30) + (s6  * am6     + r6  * bm6   ) * tmp !H_xzzz

                      temp(31) = temp(31) + (s7  * am3     + r7  * bm3   ) * tmp !H_yyyy
                      temp(32) = temp(32) + (s7  * am6     + r7  * bm6   ) * tmp !H_yyyz

                      temp(33) = temp(33) + (s8  * am6     + r8  * bm6   ) * tmp !H_yyzz

                      temp(34) = temp(34) + (s9  * am6     + r9  * bm6   ) * tmp !H_yzzz

                      temp(35) = temp(35) + (s10 * am6     + r10 * bm6   ) * tmp !H_zzzz
                   End Do
                End Do
             End Do
          End Do

       End If

       If (ll == 1) mprotx(:,iatm) = temp
       If (ll == 2) mproty(:,iatm) = temp
       If (ll == 3) mprotz(:,iatm) = temp

    End Do

  End Subroutine infinitesimal_rotation

  Subroutine rotate_mpoles_d(iatm,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for rotating multipoles from local frame to global
  ! frame for multipole order <= 2
  !
  ! Note: reference to numeric_container routines
  !
  ! copyright - daresbury laboratory
  ! author    - h.a.boateng february 2016
  ! amended   - i.t.todorov december 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent( In    )  :: iatm
    Type( comms_type ), Intent( In    ) :: comm

  ! Local variables

    Integer          :: i,j,k,m,numnbh,idi,fail,local_index
    Real( Kind = wp) :: temp1,temp2,temp3,temp4,temp5,temp6,mpole(1:mximpl)
    Real( Kind = wp) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,           &
                        tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10, &
                        p1x,p1y,p1z,p2x,p2y,p2z,dpu,          &
                        rrp1,rrp2,rrp3,rrp12,rrp22

    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rsqt

  ! Now on to permanent multipoles

    If (mplflg(iatm) == 1) Return ! multipole is already rotated

    mplflg(iatm) = 1

    idi = ltg(iatm)

    If (mxompl == 0) Then

       mplgfr(:,iatm) = mpllfr(:,lsite(iatm))

       Return

    End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find the rotation matrix a(3,3) = a(9)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    If (mprotm(iatm)%flag == 0) Then

       numnbh = ltpatm(0,iatm)

       If (numnbh == 0) Then

          mplgfr(:,iatm) = mpllfr(:,lsite(iatm))

          Return ! Atom has no bonds => global frame = local frame

       End If

       mprotm(iatm)%flag = 1

       If (numnbh == 1) Then

          j    = local_index(ltpatm(1,iatm),nlast,lsi,lsa)

          p1x  = xxx(j) - xxx(iatm)
          p1y  = yyy(j) - yyy(iatm)
          p1z  = zzz(j) - zzz(iatm)

          Call images_s(imcon,cell,p1x,p1y,p1z)

          rrp1 = Sqrt(p1x**2+p1y**2+p1z**2)

          a1   = p1x/rrp1
          a4   = p1y/rrp1
          a7   = p1z/rrp1
          a2   =-a4
          a5   = a1
          a8   = a7
          a3   = a4*a8-a5*a7
          a6   = a2*a7-a1*a8
          a9   = a1*a5-a5*a4

          rrp3 = Sqrt(a3**2+a6**2+a9**2)

          a3   = a3/rrp3
          a6   = a6/rrp3
          a9   = a9/rrp3

          mprotm(iatm)%p1 = (/ p1x, p1y, p1z /)
          mprotm(iatm)%p2 = (/-p1y, p1x, p1z /)

          mprotm(iatm)%mbnd(1) = j
          mprotm(iatm)%mbnd(2) = 0

  ! Is there no torque when atom has only one bond?

          mprotm(iatm)%flag = 1

       End If

       If (numnbh == 2) Then

          j    = local_index(ltpatm(1,iatm),nlast,lsi,lsa)
          k    = local_index(ltpatm(2,iatm),nlast,lsi,lsa)

          p1x  = xxx(j) - xxx(iatm)
          p1y  = yyy(j) - yyy(iatm)
          p1z  = zzz(j) - zzz(iatm)

          Call images_s(imcon,cell,p1x,p1y,p1z)

          p2x  = xxx(k) - xxx(iatm)
          p2y  = yyy(k) - yyy(iatm)
          p2z  = zzz(k) - zzz(iatm)

          Call images_s(imcon,cell,p2x,p2y,p2z)

          rrp12 = p1x**2+p1y**2+p1z**2
          rrp22 = p2x**2+p2y**2+p2z**2

          If ((p1x*p2x+p1y*p2y+p1z*p2z)**2 < rrp12*rrp22 - 1.0e-6_wp) Then

  ! find the bisector of the angle between the two local vectors

             p1x  = p1x + p2x
             p1y  = p1y + p2y
             p1z  = p1z + p2z

             rrp1 = Sqrt(p1x**2+p1y**2+p1z**2)

             mprotm(iatm)%p1 = (/ p1x, p1y, p1z /)
             mprotm(iatm)%p2 = (/ p2x, p2y, p2z /)

             mprotm(iatm)%mbnd(1) = j
             mprotm(iatm)%mbnd(2) = k

             a1   = p1x/rrp1
             a4   = p1y/rrp1
             a7   = p1z/rrp1

             dpu  = p2x*a1+p2y*a4+p2z*a7

             p2x  = p2x - dpu*a1
             p2y  = p2y - dpu*a4
             p2z  = p2z - dpu*a7

             rrp2 = Sqrt(p2x**2+p2y**2+p2z**2)

             a2   = p2x/rrp2
             a5   = p2y/rrp2
             a8   = p2z/rrp2

             a3   = a4*a8-a5*a7
             a6   = a2*a7-a1*a8
             a9   = a1*a5-a5*a4

             rrp3 = Sqrt(a3**2+a6**2+a9**2)

             a3   = a3/rrp3
             a6   = a6/rrp3
             a9   = a9/rrp3

          Else

             rrp1 = Sqrt(rrp12)

             a1   = p1x/rrp1
             a4   = p1y/rrp1
             a7   = p1z/rrp1
             a2   =-a4
             a5   = a1
             a8   = a7
             a3   = a4*a8-a5*a7
             a6   = a2*a7-a1*a8
             a9   = a1*a5-a5*a4

             rrp3 = Sqrt(a3*a3+a6*a6+a9*a9)

             a3   = a3/rrp3
             a6   = a6/rrp3
             a9   = a9/rrp3

             mprotm(iatm)%p1 = (/ p1x, p1y, p1z /)
             mprotm(iatm)%p2 = (/-p1y, p1x, p1z /)

             mprotm(iatm)%mbnd(1) = j
             mprotm(iatm)%mbnd(2) = k

  ! Is there no torque when the bond angle is pi?

             mprotm(iatm)%flag = 1

          End If

       End If

       If (numnbh > 2) Then

          fail = 0
          Allocate (xxt(1:numnbh),yyt(1:numnbh),zzt(1:numnbh),rsqt(1:numnbh), Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'allocation failure in rotate_mpoles_d, node: ', comm%idnode
             Call error(0)
          End If

          Do i = 1,numnbh

             j       = local_index(ltpatm(i,iatm),nlast,lsi,lsa)

             xxt(i)  = xxx(j) - xxx(iatm)
             yyt(i)  = yyy(j) - yyy(iatm)
             zzt(i)  = zzz(j) - zzz(iatm)

             Call images_s(imcon,cell,xxt(i),yyt(i),zzt(i))

             rsqt(i) = xxt(i)**2 + yyt(i)**2 + zzt(i)**2

          End Do

          k = Minloc(rsqt,Dim=1)
          rrp12 = rsqt(k)
          rsqt(k) = 0.0_wp

          m = Minloc(rsqt,Dim=1,Mask=(rsqt > 0.0_wp))
          rrp22 = rsqt(m)

          Do While ((xxt(k)*xxt(m)+yyt(k)*yyt(m)+zzt(k)*zzt(m))**2 >= rrp12*rrp22 - 1.0e-6_wp)

             rsqt(m) = 0.0_wp
             m = Minloc(rsqt,Dim=1,Mask=(rsqt > 0.0_wp))
             rrp22 = rsqt(m)

          End Do

          mprotm(iatm)%mbnd(1) = k
          mprotm(iatm)%mbnd(2) = m

          p1x  = xxt(k) + xxt(m)
          p1y  = yyt(k) + yyt(m)
          p1z  = zzt(k) + zzt(m)

          rrp1 = Sqrt(p1x**2+p1y**2+p1z**2)

          mprotm(iatm)%p1 = (/ p1x   , p1y   , p1z    /)
          mprotm(iatm)%p2 = (/ xxt(m), yyt(m), zzt(m) /)

          a1   = p1x/rrp1
          a4   = p1y/rrp1
          a7   = p1z/rrp1

          dpu  = xxt(m)*a1+yyt(m)*a4+zzt(m)*a7

          p2x  = xxt(m) - dpu*a1
          p2y  = yyt(m) - dpu*a4
          p2z  = zzt(m) - dpu*a7

          rrp2 = Sqrt(p2x**2+p2y**2+p2z**2)

          a2   = p2x/rrp2
          a5   = p2y/rrp2
          a8   = p2z/rrp2

          a3   = a4*a8-a5*a7
          a6   = a2*a7-a1*a8
          a9   = a1*a5-a5*a4

          rrp3 = Sqrt(a3**2+a6**2+a9**2)

          a3   = a3/rrp3
          a6   = a6/rrp3
          a9   = a9/rrp3

          Deallocate (xxt,yyt,zzt,rsqt, Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'deallocation failure in rotate_mpoles_d, node: ', comm%idnode
             Call error(0)
          End If

       End If

       mprotm(iatm)%mtrxa(1) = a1 ; mprotm(iatm)%mtrxa(2) = a2 ; mprotm(iatm)%mtrxa(3) = a3
       mprotm(iatm)%mtrxa(4) = a4 ; mprotm(iatm)%mtrxa(5) = a5 ; mprotm(iatm)%mtrxa(6) = a6
       mprotm(iatm)%mtrxa(7) = a7 ; mprotm(iatm)%mtrxa(8) = a2 ; mprotm(iatm)%mtrxa(9) = a9

    Else

       a1 = mprotm(iatm)%mtrxa(1) ; a2 = mprotm(iatm)%mtrxa(2) ; a3 = mprotm(iatm)%mtrxa(3)
       a4 = mprotm(iatm)%mtrxa(4) ; a5 = mprotm(iatm)%mtrxa(5) ; a6 = mprotm(iatm)%mtrxa(6)
       a7 = mprotm(iatm)%mtrxa(7) ; a8 = mprotm(iatm)%mtrxa(8) ; a9 = mprotm(iatm)%mtrxa(9)

    End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now on to rotate multipoles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! get the local multipoles based on atom type

    mpole(:) = mpllfr(:,lsite(iatm))

  ! rotate monopole

    mplgfr(1,iatm) = mpole(1)

    If (mxompl >= 1) Then

  ! rotate dipoles

       tt2=mpole(2) ; tt3=mpole(3) ; tt4=mpole(4)

       mplgfr(2,iatm)  = a1*tt2 + a2*tt3 + a3*tt4
       mplgfr(3,iatm)  = a4*tt2 + a5*tt3 + a6*tt4
       mplgfr(4,iatm)  = a7*tt2 + a8*tt3 + a9*tt4

    End If

    If (mxompl == 2) Then

  ! rotate quadrupoles

       tt5=mpole(5) ; tt6=mpole(6) ; tt7=mpole(7)
       tt8=mpole(8) ; tt9=mpole(9) ; tt10=mpole(10)

       temp1 = a4*tt5+a5*tt6+a6*tt7 ; temp2=a4*tt6+a5*tt8+a6*tt9 ; temp3 = a4*tt7+a5*tt9+a6*tt10
       temp4 = a7*tt5+a8*tt6+a9*tt7 ; temp5=a7*tt6+a8*tt8+a9*tt9 ; temp6 = a7*tt7+a8*tt9+a9*tt10

       mplgfr(5,iatm)  =         a1*(a1*tt5+a2*tt6+a3*tt7) + &
                                 a2*(a1*tt6+a2*tt8+a3*tt9) + &
                                 a3*(a1*tt7+a2*tt9+a3*tt10)
       mplgfr(6,iatm)  = 2.0_wp*(a1*temp1+a2*temp2+a3*temp3)
       mplgfr(7,iatm)  = 2.0_wp*(a1*temp4+a2*temp5+a3*temp6)
       mplgfr(8,iatm)  =         a4*temp1+a5*temp2+a6*temp3
       mplgfr(9,iatm)  = 2.0_wp*(a4*temp4+a5*temp5+a6*temp6)
       mplgfr(10,iatm) =         a7*temp4+a8*temp5+a9*temp6

    End If

    mpole(:) = mplgfr(:,iatm)

  ! Find infinitesimal rotation of global multipoles
  !==================================================

    mprotx(1,iatm) = 0.0_wp
    mproty(1,iatm) = 0.0_wp
    mprotz(1,iatm) = 0.0_wp

    If (mxompl >= 1) Then

  ! x-coordinate

       mprotx(2,iatm)  = 0.0_wp;
       mprotx(3,iatm)  = mpole(4)
       mprotx(4,iatm)  =-mpole(3)

  ! y-coordinate

       mproty(2,iatm)  =-mpole(4)
       mproty(3,iatm)  = 0.0_wp
       mproty(4,iatm)  = mpole(2)

  ! z-coordinate

       mprotz(2,iatm)  = mpole(3)
       mprotz(3,iatm)  =-mpole(2)
       mprotz(4,iatm)  = 0.0_wp

    End If

    If (mxompl == 2) Then

  ! x-coordinate

       mprotx(5,iatm)  = 0.0_wp
       mprotx(6,iatm)  = mpole(7)
       mprotx(7,iatm)  =-mpole(6)
       mprotx(8,iatm)  = 2.0_wp*mpole(9)
       mprotx(9,iatm)  = mpole(10)-mpole(8)
       mprotx(10,iatm) =-2.0_wp*mpole(9)

  ! y-coordinate

       mproty(5,iatm)  =-2.0_wp*mpole(7)
       mproty(6,iatm)  =-mpole(9)
       mproty(7,iatm)  = mpole(5)-mpole(10)
       mproty(8,iatm)  = 0.0_wp
       mproty(9,iatm)  = mpole(6)
       mproty(10,iatm) = 2.0_wp*mpole(7)

  ! z-coordinate

       mprotz(5,iatm)  = 2.0_wp*mpole(6)
       mprotz(6,iatm)  = mpole(8)-mpole(5)
       mprotz(7,iatm)  = mpole(9)
       mprotz(8,iatm)  =-2.0_wp*mpole(6)
       mprotz(9,iatm)  =-mpole(7)
       mprotz(10,iatm) = 0.0_wp

    End If

  End Subroutine rotate_mpoles_d

End Module mpoles_container
