Subroutine ewald_spme_mfield(alpha,epsq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field
! due to multipolar interactions in a periodic system using the smooth
! particle mesh ewald method for multipoles
!
! This version allows for extension to arbitrary order
!
! Note: (fourier) reciprocal space terms
!
! copyright - daresbury laboratory
! author    - h.a.boateng december 2014
! amended   - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,   Only : idnode,mxnode,gcheck,gsum,dlp_comm_world
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz,idx,idy,idz
  Use configuration,  Only : cell,natms,nlast,xxx,yyy,zzz
  Use mpoles_module,  Only : mplmap,mplgfr
  Use parallel_fft

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: alpha,epsq

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ixb,iyb,izb, ixt,iyt,izt
  Real( Kind = wp ), Save :: ixbm1_r,iybm1_r,izbm1_r, &
                             ixtm0_r,iytm0_r,iztm0_r, &
                             kmaxa_r,kmaxb_r,kmaxc_r

  Logical              :: llspl=.true.
  Integer              :: fail(1:4), i,j,k,l, jj,kk,ll, jjb,jjt, kkb,kkt, llb,llt, &
                          jjtjjb,counter,s1,s2,s3

  Real( Kind = wp )    :: Dtpbsp,det,rcell(1:9),celprp(1:10),ralph,             &
                          rcpcut,rcpct2,tmp,bb1,bb2,bb3,mptmp,                  &
                          rksq,rkx1,rkx2,rkx3, rky1,rky2,rky3, rkz1,rkz2,rkz3,  &
                          dtp

  Real( Kind = wp )    :: imp(1:mximpl)

  Complex( Kind = wp ) :: vterm

! blocking factors for fft

  Integer,           Save :: block_x,block_y,block_z

! B-spline coefficients

  Complex( Kind = wp ), Dimension( : ),     Allocatable, Save :: bscx,bscy,bscz
  Complex( Kind = wp ), Dimension( : ),     Allocatable       :: ww1,ww2,ww3

  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: csp
  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: txx,tyy,tzz
  Integer,              Dimension( : ),     Allocatable       :: ixx,iyy,izz,it
  Real( Kind = wp ),    Dimension( :,:,: ), Allocatable       :: bsddx,bsddy,bsddz
  Real( Kind = wp ),    Dimension( :,: ),   Allocatable       :: bspx,bspy,bspz
  Real( Kind = wp ),    Dimension( : ),     Allocatable       :: bdx,bdy,bdz

! context for parallel fft

  Integer,           Save :: context

! indexing arrays for x, y & z as used in parallel fft

  Integer,              Dimension( : ),     Allocatable, Save :: index_x,index_y,index_z

! temporary qqc

  Real( Kind = wp )    :: qqc_tmp

! temporary workspace for parallel fft

  Real( Kind = wp ),    Dimension( :,:,: ), Allocatable, Save :: qqc_local
  Complex( Kind = wp ), Dimension( :,:,: ), Allocatable, Save :: qqq_local
  Complex( Kind = wp ), Dimension( :,:,: ), Allocatable, Save :: pfft_work

! DaFT arrays local indices

  Integer              :: j_local, k_local, l_local


  fail=0
  If (newjob) Then
     newjob = .false.

!!! BEGIN DD SPME VARIABLES
! 3D charge array construction (bottom and top) indices

     ixb=idx*(kmaxa/nprx)+1
     ixt=(idx+1)*(kmaxa/nprx)
     iyb=idy*(kmaxb/npry)+1
     iyt=(idy+1)*(kmaxb/npry)
     izb=idz*(kmaxc/nprz)+1
     izt=(idz+1)*(kmaxc/nprz)

     ixbm1_r=Real(ixb-1,wp)
     ixtm0_r=Nearest( Real(ixt,wp) , -1.0_wp )
     iybm1_r=Real(iyb-1,wp)
     iytm0_r=Nearest( Real(iyt,wp) , -1.0_wp )
     izbm1_r=Real(izb-1,wp)
     iztm0_r=Nearest( Real(izt,wp) , -1.0_wp )

! Real values of kmax vectors

     kmaxa_r=Real(kmaxa,wp)
     kmaxb_r=Real(kmaxb,wp)
     kmaxc_r=Real(kmaxc,wp)

!!! END DD SPME VARIABLES

!!! BEGIN CARDINAL B-SPLINES SET-UP
! allocate the complex exponential arrays

     Allocate (ww1(1:kmaxa),ww2(1:kmaxb),ww3(1:kmaxc), Stat = fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ww arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! initialise the complex exponential arrays

     Call spl_cexp(kmaxa,kmaxb,kmaxc,ww1,ww2,ww3)

! allocate the global B-spline coefficients and the helper array

     Allocate (bscx(1:kmaxa),bscy(1:kmaxb),bscz(1:kmaxc), Stat = fail(1))
     Allocate (csp(1:mxspl),                              Stat = fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'bsc and cse arrays allocation failure, node: ', idnode
        Call error(0)
     End If

! calculate the global B-spline coefficients

     Call bspcoe(mxspl,kmaxa,kmaxb,kmaxc,csp,bscx,bscy,bscz,ww1,ww2,ww3)

! deallocate the helper array and complex exponential arrays

     Deallocate (csp,         Stat = fail(1))
     Deallocate (ww1,ww2,ww3, Stat = fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'cse and ww arrays deallocation failure, node: ', idnode
        Call error(0)
     End If

!!! END CARDINAL B-SPLINES SET-UP

!!! BEGIN DAFT SET-UP
! domain local block limits of kmax space

     block_x = kmaxa / nprx
     block_y = kmaxb / npry
     block_z = kmaxc / nprz

! set up the parallel fft and useful related quantities

     Call initialize_fft( 3, (/ kmaxa, kmaxb, kmaxc /), &
         (/ nprx, npry, nprz /), (/ idx, idy, idz /),   &
         (/ block_x, block_y, block_z /),               &
         dlp_comm_world, context )

! set up the indexing arrays for each dimension (NOT deallocated manually)

     Allocate ( index_x( 1:block_x ), Stat = fail(1) )
     Allocate ( index_y( 1:block_y ), Stat = fail(2) )
     Allocate ( index_z( 1:block_z ), Stat = fail(3) )
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'SPME index arrays allocation failure, node: ', idnode
        Call error(0)
     End If

     Call pfft_indices( kmaxa, block_x, idx, nprx, index_x )
     Call pfft_indices( kmaxb, block_y, idy, npry, index_y )
     Call pfft_indices( kmaxc, block_z, idz, nprz, index_z )

! workspace arrays for DaFT

     Allocate ( qqc_local( 1:block_x, 1:block_y, 1:block_z ), Stat = fail(1) )
     Allocate ( qqq_local( 1:block_x, 1:block_y, 1:block_z ), Stat = fail(2) )
     Allocate ( pfft_work( 1:block_x, 1:block_y, 1:block_z ), Stat = fail(3) )
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'SPME DaFT workspace arrays allocation failure, node: ', idnode
        Call error(0)
     End If

!!! END DAFT SET-UP
  End If

  Allocate (txx(1:mxatms),tyy(1:mxatms),tzz(1:mxatms),                            Stat = fail(1))
  Allocate (ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms),it(1:mxatms),               Stat = fail(2))
  Allocate (bdx(0:mxspl),bsddx(0:mxspl,1:mxspl,1:mxatms), &
            bdy(0:mxspl),bsddy(0:mxspl,1:mxspl,1:mxatms), &
            bdz(0:mxspl),bsddz(0:mxspl,1:mxspl,1:mxatms),                         Stat = fail(3))
  Allocate (bspx(1:mxspl,1:mxatms),bspy(1:mxspl,1:mxatms),bspz(1:mxspl,1:mxatms), Stat = fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mfield allocation failure, node: ', idnode
     Call error(0)
  End If

! set working parameters

  ralph=-0.25_wp/alpha**2

! Convert cell coordinates to fractional coordinates intervalled [0,1)
! (bottom left corner of MD cell) and stretch over kmaxs in different
! directions.  Only the halo (natms,nlast] has fractional coordinates
! outside the [0,1) interval.  In the worst case scenario of one
! "effective" link-cell per domain and one domain in the MD cell only,
! the halo will have fractional coordinates intervalled as
! [n,0)u[1,2), where -1 <= n < 0.  Only the positive halo is needed by
! the B-splines since they distribute/spread charge density in
! negative direction with length the length of the spline.
!
! The story has become more complicated with cutoff padding and the
! conditional updates of the VNL and thus the halo as now a domain
! (1:natms) particle can enter the halo and vice versa.  So DD
! bounding is unsafe!!!

  Call invert(cell,rcell,det)
  If (Abs(det) < 1.0e-6_wp) Call error(120)

  Do i=1,nlast
     txx(i)=kmaxa_r*(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)+0.5_wp)
     tyy(i)=kmaxb_r*(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)+0.5_wp)
     tzz(i)=kmaxc_r*(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)+0.5_wp)

! If not DD bound in kmax grid space when .not.llvnl = (mxspl1 == mxspl)

     If (mxspl1 == mxspl .and. i <= natms) Then
        If (txx(i) < ixbm1_r .or. txx(i) > ixtm0_r .or. &
            tyy(i) < iybm1_r .or. tyy(i) > iytm0_r .or. &
            tzz(i) < izbm1_r .or. tzz(i) > iztm0_r) llspl=.false.
     End If

     ixx(i)=Int(txx(i))
     iyy(i)=Int(tyy(i))
     izz(i)=Int(tzz(i))

! Detect if a particle is charged and in the MD cell or in its positive halo
! (t(i) >= -zero_plus) as the B-splines are negative directionally by propagation

     If (tzz(i) >= -zero_plus .and. &
         tyy(i) >= -zero_plus .and. &
         txx(i) >= -zero_plus .and. &
         Maxval(Abs(mplgfr(:,i))) > zero_plus)  Then
        it(i)=1
     Else
        it(i)=0
     End If
  End Do

! Check for breakage of llspl when .not.llvnl = (mxspl1 == mxspl)

  mxspl2=mxspl1
  If (mxspl1 == mxspl) Then
     If (mxnode > 1) Call gcheck(llspl)
     If (.not.llspl) mxspl2=mxspl+1
  End If

! construct B-splines for atoms

  Call bspgen_mpoles(nlast,mxspl,txx,tyy,tzz,bspx,bspy,bspz,bsddx,bsddy,bsddz)

  Deallocate (txx,tyy,tzz, Stat = fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mfield allocation failure, node: ', idnode
     Call error(0)
  End If

! zero 3D charge array
! DaFT version - only need set local bit to zero

  qqc_local = 0.0_wp

! construct 3D charge array
! DaFT version - use array that holds only the local data

  If (mxompl < 5) Then

     Do i=1,nlast

! If a particle is charged and in the MD cell or in its positive halo
! (t(i) >= 0) as the B-splines are negative directionally by propagation

        If (it(i) == 1) Then

! get the multipoles for site i

           imp=mplgfr(:,i)

           llb = Max( izb, izz(i) - mxspl + 2 )
           llt = Min( izt, izz(i) + 1 )

           kkb = Max( iyb, iyy(i) - mxspl + 2 )
           kkt = Min( iyt, iyy(i) + 1 )

           jjb = Max( ixb, ixx(i) - mxspl + 2 )
           jjt = Min( ixt, ixx(i) + 1 )

           jjtjjb = jjt - jjb + 1

           Do ll = llb, llt
              l = izz(i) - ll + 2

              l_local = ll - izb + 1

              bdz=bsddz(:,l,i)

              Do kk = kkb, kkt
                 k = iyy(i) - kk + 2

                 k_local = kk - iyb + 1

                 bdy=bsddy(:,k,i)

                 If (jjtjjb > 0 .and. jjtjjb <= 8) Then

                    jj = jjb

                    j = ixx(i) - jj + 2

                    j_local = jj - ixb + 1

                    bdx=bsddx(:,j,i)

                    Call explicit_spme_loop_s(rcell,bdx,bdy,bdz,imp,dtp)

                    qqc_local(j_local,k_local,l_local) = qqc_local(j_local,k_local,l_local) + dtp

                    counter = 1

                    Do While (counter < jjtjjb)
                       j = j - 1

                       j_local = j_local + 1

                       bdx=bsddx(:,j,i)

                       Call explicit_spme_loop_s(rcell,bdx,bdy,bdz,imp,dtp)

                       qqc_local(j_local,k_local,l_local) = qqc_local(j_local,k_local,l_local) + dtp

                       counter = counter + 1
                    End Do

                  Else

                    Do jj = jjb, jjt
                       j = ixx(i) - jj + 2

                       j_local = jj - ixb + 1

                       bdx=bsddx(:,j,i)

                       Call explicit_spme_loop_s(rcell,bdx,bdy,bdz,imp,dtp)

                       qqc_local(j_local,k_local,l_local) = qqc_local(j_local,k_local,l_local) + dtp
                    End Do

                  End If
               End Do
           End Do

        End If

     End Do

  Else

     Do i=1,nlast

! If a particle is charged and in the MD cell or in its positive halo
! (t(i) >= 0) as the B-splines are negative directionally by propagation

        If (it(i) == 1) Then

! get the multipoles for site i

           imp=mplgfr(:,i)

           llb = Max( izb, izz(i) - mxspl + 2 )
           llt = Min( izt, izz(i) + 1 )

           kkb = Max( iyb, iyy(i) - mxspl + 2 )
           kkt = Min( iyt, iyy(i) + 1 )

           jjb = Max( ixb, ixx(i) - mxspl + 2 )
           jjt = Min( ixt, ixx(i) + 1 )

           jjtjjb = jjt - jjb + 1

           Do s3 = 0, mxompl

              Do s2 = 0, mxompl - s3

                 Do s1 = 0, mxompl - s3 - s2

                    mptmp = imp(mplmap(s1,s2,s3))

                    Do ll = llb, llt
                       l = izz(i) - ll + 2

                       l_local = ll - izb + 1

                       bdz=bsddz(:,l,i)

                       Do kk = kkb, kkt
                          k = iyy(i) - kk + 2

                          k_local = kk - iyb + 1

                          bdy=bsddy(:,l,i)

                          If (jjtjjb > 0 .and. jjtjjb <= 8) Then

                             jj = jjb

                             j = ixx(i) - jj +2

                             j_local = jj - ixb + 1

                             bdx=bsddx(:,j,i)

                             tmp = mptmp*Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

                             qqc_local(j_local,k_local,l_local) = qqc_local(j_local,k_local,l_local) + tmp

                             counter = 1

                             Do While (counter < jjtjjb)
                                j = j - 1

                                j_local = j_local + 1

                                bdx=bsddx(:,j,i)

                                tmp = mptmp*Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

                                qqc_local(j_local,k_local,l_local) = qqc_local(j_local,k_local,l_local) + tmp

                                counter = counter + 1
                             End Do

                          Else

                             Do jj = jjb, jjt
                                j = ixx(i) - jj + 2

                                j_local = jj - ixb + 1

                                bdx=bsddx(:,j,i)

                                tmp = mptmp*Dtpbsp(s1,s2,s3,rcell,bdx,bdy,bdz)

                                qqc_local(j_local,k_local,l_local) = qqc_local(j_local,k_local,l_local) + tmp

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

  qqq_local=Cmplx(qqc_local , Kind = wp)

! calculate inverse 3D FFT of charge array (in place)

  Call pfft(qqq_local,pfft_work,context,1)

! set reciprocal space cutoff

  Call dcell(rcell,celprp)

  rcpcut=0.5_wp*Min(kmaxa_r*celprp(7),kmaxb_r*celprp(8),kmaxc_r*celprp(9))
  rcpcut=rcpcut*1.05_wp*twopi
  rcpct2=rcpcut**2

! calculate convolution of charge array with gaussian function
! DaFT Version - only loop over the local stuff

  Do l_local=1,block_z
     l=index_z(l_local)

     ll=l-1
     If (l > kmaxc/2) ll=ll-kmaxc
     tmp=twopi*Real(ll,wp)

     rkx1=tmp*rcell(3)
     rky1=tmp*rcell(6)
     rkz1=tmp*rcell(9)

     bb3=Real( bscz(l)*Conjg(bscz(l)),wp )

     Do k_local=1,block_y
        k=index_y(k_local)

        kk=k-1
        If (k > kmaxb/2) kk=kk-kmaxb
        tmp=twopi*Real(kk,wp)

        rkx2=rkx1+tmp*rcell(2)
        rky2=rky1+tmp*rcell(5)
        rkz2=rkz1+tmp*rcell(8)

        bb2=bb3*Real( bscy(k)*Conjg(bscy(k)),wp )

        Do j_local=1,block_x
           j=index_x(j_local)

           jj=j-1
           If (j > kmaxa/2) jj=jj-kmaxa
           tmp=twopi*Real(jj,wp)

           rkx3=rkx2+tmp*rcell(1)
           rky3=rky2+tmp*rcell(4)
           rkz3=rkz2+tmp*rcell(7)

           bb1=bb2*Real( bscx(j)*Conjg(bscx(j)),wp )

           rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

           If (rksq > 1.0e-6_wp .and. rksq <= rcpct2) Then

              vterm=bb1*Exp(ralph*rksq)/rksq*qqq_local(j_local,k_local,l_local)

              qqq_local(j_local,k_local,l_local)=vterm

           Else

              qqq_local(j_local,k_local,l_local)=(0.0_wp,0.0_wp)

           End If
        End Do
     End Do
  End Do

  Call pfft(qqq_local,pfft_work,context,-1)

! calculate electrostatic field

  Call spme_mfield(rcell, ixx,iyy,izz, bsddx,bsddy,bsddz, qqc_local,ixb,ixt, iyb,iyt, izb,izt)

  Deallocate (ixx,iyy,izz,it,    Stat = fail(1))
  Deallocate (bdx,bdy,bdz,       Stat = fail(2))
  Deallocate (bsddx,bsddy,bsddz, Stat = fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'ewald_spme_mfield deallocation failure, node: ', idnode
     Call error(0)
  End If

Contains

  Subroutine spme_mfield(rcell, ixx,iyy,izz, bsddx,bsddy,bsddz, qqc_local, ixb,ixt, iyb,iyt, izb,izt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the electrostatic field due to
! multipolar interactions in periodic system using the smooth particle
! mesh ewald method (fourier part)
!
! Note: qqc_local is shifted from its definition from above
!       and therefore there is no need for periodic images (!!)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov december 2016
! amended   - h.a.boateng december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use mpoles_module,  Only : mpfldx,mpfldy,mpfldz

    Implicit None

    Integer,           Intent( In    ) :: ixx(1:mxatms),iyy(1:mxatms),izz(1:mxatms), &
                                          ixb,ixt, iyb,iyt, izb,izt
    Real( Kind = wp ), Intent( In    ) :: rcell(1:9),                      &
                                          bsddx(0:mxspl,1:mxspl,1:mxatms), &
                                          bsddy(0:mxspl,1:mxspl,1:mxatms), &
                                          bsddz(0:mxspl,1:mxspl,1:mxatms), &
                                          qqc_local( ixb:ixt, iyb:iyt, izb:izt )

    Integer           :: fail(1:2), delspl, ixdb,iydb,izdb,ixdt,iydt,izdt, &
                         i,j,k,l, jj,kk,ll
    Real( Kind = wp ) :: Dtpbsp,fix,fiy,fiz,qsum

    Real( Kind = wp ), Dimension( : ),     Allocatable :: bdx,bdy,bdz
    Real( Kind = wp ), Dimension( :,:,: ), Allocatable :: qqc_domain

    ixdb = ixb - mxspl
    iydb = iyb - mxspl
    izdb = izb - mxspl

    delspl = mxspl2 - mxspl

    ixdt = ixt + delspl
    iydt = iyt + delspl
    izdt = izt + delspl

    fail=0
    Allocate (bdx(0:mxspl),bdy(0:mxspl),bdz(0:mxspl),        Stat = fail(1))
    Allocate (qqc_domain( ixdb:ixdt, iydb:iydt, izdb:izdt ), Stat = fail(2))
    If (Any(fail > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_mfield allocation failure, node: ', idnode
    End If

    Call exchange_grid( ixb , ixt , iyb , iyt , izb , izt , qqc_local , &
                        ixdb, iydb, izdb, ixdt, iydt, izdt, qqc_domain  )

    Do i=1,natms

! initialise field

       fix=0.0_wp ; fiy=0.0_wp ; fiz=0.0_wp

       Do l=1,mxspl
          ll=izz(i)-l+2

          bdz=bsddz(:,l,i)

          Do k=1,mxspl
             kk=iyy(i)-k+2

             bdy=bsddy(:,k,i)

             Do j=1,mxspl
                jj=ixx(i)-j+2

                bdx=bsddx(:,j,i)

                qsum=qqc_domain(jj,kk,ll)

! field
                fix=fix+qsum*Dtpbsp(1,0,0,rcell,bdx,bdy,bdz)
                fiy=fiy+qsum*Dtpbsp(0,1,0,rcell,bdx,bdy,bdz)
                fiz=fiz+qsum*Dtpbsp(0,0,1,rcell,bdx,bdy,bdz)
             End Do
          End Do
       End Do

! load field

       mpfldx(i)=mpfldx(i)+fix
       mpfldy(i)=mpfldy(i)+fiy
       mpfldz(i)=mpfldz(i)+fiz

    End Do

    Deallocate (bdx,bdy,bdz, Stat = fail(1))
    Deallocate (qqc_domain,  Stat = fail(2))
    If (Any(fail > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'spme_mfield dealocation failure, node: ', idnode
    End If

  End Subroutine spme_mfield

End Subroutine ewald_spme_mfield
