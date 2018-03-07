Module poisson_module

  Use kinds, Only : wp
  Use comms, Only : gsum, comms_type,wp_mpi,ExchgGrid_tag
  Use domains_module
  Use setup_module,  Only : fourpi,r4pie0,nrite,            &
                            kmaxa,kmaxb,kmaxc,mxspl,mxspl1, &
                            mxlist,mxatms,mxatdm,half_minus,zero_plus
  Use configuration, Only : imcon,cell,natms,nlast,list,ltg,lfrzn, &
                            chge,xxx,yyy,zzz,fxx,fyy,fzz
  Use ewald_module
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None

  Logical :: debug =.true., pinitialized=.false.
  Logical :: pconverged

  Integer :: lnxl,lnxu, lnyl,lnyu, lnzl,lnzu, maxsteps=50, maxbicgst, icut
  Integer :: ixb,iyb,izb, ixt,iyt,izt, block_x,block_y,block_z
  Integer :: lmap(1:26) ! local copy of map to be modified
  Integer :: mxitjb, mxitcg, xhalo, fail(1:10)=0

  Real( Kind = wp ) :: eps,alpha,delta,celprp(1:10),epsilon,  &
                       bchange=1.0e-1_wp,deltaB_norm,         &
                       normb, normphi0, normphi1,             &
                       ccxxx,ccyyy,cczzz, ccsum,              &
                       kmaxa_r,kmaxb_r,kmaxc_r

  Real( Kind = wp ), Dimension(:,:,:), Allocatable :: phi ! electrostatic potential
  Real( Kind = wp ), Dimension(:,:,:), Allocatable :: b   ! charge density
  Real( Kind = wp ), Dimension(:,:,:), Allocatable :: F, s, r, p, r0, v, t ! biCGStab variables
  Real( Kind = wp ), Dimension(:,:,:), Allocatable :: pattern, phi0

Contains

  Subroutine poisson_forces(alphain,epsq,engcpe,vircpe,stress,comm)

    Real( Kind = wp ), Intent( In    ) :: epsq,alphain
    Real( Kind = wp ), Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Type(comms_type), Intent( InOut)   :: comm

    Real( Kind = wp ) :: eng,virr
    Real( Kind = wp ) :: strs(1:9)

    If (.not.pinitialized) Then

       pconverged=.false.
       maxbicgst =0 ! the number of steps to solve eq. without a guess (phi=0)
       maxsteps  =mxitcg
       epsilon=epsq
       alpha=alphain
       delta=1.0_wp/alpha
       Call biCGStab_init(comm)

       Call biCGStab_charge_density(comm)

       If (normb > zero_plus) Then
          Call biCGStab_solver_omp(eps,comm)
       Else
          STOP "normb too small"
       End If

       pinitialized=.true.

    Else

       Call biCGStab_charge_density(comm)

    End If

    If (normb > zero_plus) Then
       Call P_solver_omp(eps,comm)
       Call biCGStab_calc_forces(eng,virr,strs,comm)
    Else
       STOP "normb too small"
    End If

    engcpe=engcpe+eng
    vircpe=vircpe+virr
    stress=stress+strs

  End Subroutine poisson_forces

  Subroutine biCGStab_init(comm)

! calculates preambles
   Type(comms_type), Intent( InOut ) :: comm
! copy DD mapping

    lmap=map

    If (imcon == 0 .or. imcon == 6) Then
       If      (imcon == 0) Then
          If (idx == 0     ) lmap(1)=-1
          If (idy == 0     ) lmap(3)=-1
          If (idz == 0     ) lmap(5)=-1
          If (idx == nprx-1) lmap(2)=-1
          If (idy == npry-1) lmap(4)=-1
          If (idz == nprz-1) lmap(6)=-1
       Else If (imcon == 6) Then
          If (idz == 0     ) lmap(5)=-1
          If (idz == nprz-1) lmap(6)=-1
       End If
    End If

! (re)set grid size along x, y and z direction

    Call dcell(cell,celprp)

    kmaxa=Nint(celprp(7)/delta)
    kmaxb=Nint(celprp(8)/delta)
    kmaxc=Nint(celprp(9)/delta)

! adjust accordingly to processor grid restrictions

    Call adjust_kmax(kmaxa, nprx)
    Call adjust_kmax(kmaxb, npry)
    Call adjust_kmax(kmaxc, nprz)

! 3D charge array construction (bottom and top) indices and block size

    ixb=idx*(kmaxa/nprx)+1
    ixt=(idx+1)*(kmaxa/nprx)
    iyb=idy*(kmaxb/npry)+1
    iyt=(idy+1)*(kmaxb/npry)
    izb=idz*(kmaxc/nprz)+1
    izt=(idz+1)*(kmaxc/nprz)

    kmaxa_r=Real(kmaxa, wp)
    kmaxb_r=Real(kmaxb, wp)
    kmaxc_r=Real(kmaxc, wp)

    block_x=kmaxa/nprx
    block_y=kmaxb/npry
    block_z=kmaxc/nprz

! new 1 grid link halo inclusive object per domain distributed sizes

    lnxl=1-1
    lnyl=1-1
    lnzl=1-1
    lnxu=block_x+1
    lnyu=block_y+1
    lnzu=block_z+1

! extra halo on top of the one grid cell halo, according to size of differentiation stencil

    xhalo=2+mxspl1-mxspl

! allocate vectors for biCGStab

    Allocate ( t(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) ,  Stat = fail( 1) )
    Allocate ( b(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) ,  Stat = fail( 2) )
    Allocate ( F(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) ,  Stat = fail( 3) )
    Allocate ( s(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) ,  Stat = fail( 4) )
    Allocate ( r(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) ,  Stat = fail( 5) )
    Allocate ( p(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) ,  Stat = fail( 6) )
    Allocate (  phi(lnxl-xhalo:lnxu+xhalo, &
                    lnyl-xhalo:lnyu+xhalo, &
                    lnzl-xhalo:lnzu+xhalo) ,       Stat = fail( 7) )
    Allocate ( phi0(lnxl-xhalo:lnxu+xhalo, &
                    lnyl-xhalo:lnyu+xhalo, &
                    lnzl-xhalo:lnzu+xhalo) ,       Stat = fail( 8) )
    Allocate ( r0(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) , Stat = fail( 9) )
    Allocate (  v(lnxl:lnxu,lnyl:lnyu,lnzl:lnzu) , Stat = fail(10) )
    If (Any(fail(1:10) > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'biCGStab_init allocation failure, node: ', comm%idnode
       Call error(0)
    End If

! Initialise

    F   = 0.0_wp
    s   = 0.0_wp
    r   = 0.0_wp
    p   = 0.0_wp
    phi = 0.0_wp
    r0  = 0.0_wp
    v   = 0.0_wp
    t   = 0.0_wp

  End Subroutine biCGStab_init

  Subroutine biCGStab_charge_density(comm)

! calculates charge dansity at 0th order

    Type(comms_type), Intent( InOut ) :: comm 
    Integer           :: i,j,k,n
    Real( Kind = wp ) :: reps0dv, txx,tyy,tzz, det,rcell(9)

    reps0dv=fourpi*r4pie0*alpha/epsilon ! dv collapsed to dr=delta=1/alpha

    ! get reciprocal cell
    Call invert(cell,rcell,det)
    If (Abs(det) < 1.0e-6_wp) Call error(120)

    b=0.0_wp

    ccsum=0.0_wp
    ccxxx=0.0_wp
    ccyyy=0.0_wp
    cczzz=0.0_wp

    normb=0.0_wp

    !$omp paralleldo default(shared) private(n) reduction(+:normb,b,ccsum,ccxxx,ccyyy,cczzz)
    Do n=1,natms !nlast
       txx=kmaxa_r*(rcell(1)*xxx(n)+rcell(4)*yyy(n)+rcell(7)*zzz(n)+0.5_wp)
       tyy=kmaxb_r*(rcell(2)*xxx(n)+rcell(5)*yyy(n)+rcell(8)*zzz(n)+0.5_wp)
       tzz=kmaxc_r*(rcell(3)*xxx(n)+rcell(6)*yyy(n)+rcell(9)*zzz(n)+0.5_wp)

! global indeces

       i=Int(txx)
       j=Int(tyy)
       k=Int(tzz)

! get local indces

       i = i - ixb + 2
       j = j - iyb + 2
       k = k - izb + 2

       If ( (i < 1 .or. i > block_x) .or. &
            (j < 1 .or. j > block_y) .or. &
            (k < 1 .or. k > block_z) .or. &
            (Abs(chge(n)) <= zero_plus) ) Cycle

       b(i,j,k)=b(i,j,k)+chge(n)*reps0dv

       ccsum=ccsum+Abs(chge(n))
       ccxxx=ccxxx+Abs(chge(n))*xxx(n)
       ccyyy=ccyyy+Abs(chge(n))*yyy(n)
       cczzz=cczzz+Abs(chge(n))*zzz(n)
    End Do
    !$omp End paralleldo

    normb=Sum(b**2)
    Call gsum(comm,normb)

    Call gsum(comm,ccsum)
    Call gsum(comm,ccxxx)
    Call gsum(comm,ccyyy)
    Call gsum(comm,cczzz)

  End Subroutine biCGStab_charge_density

  Recursive Subroutine P_solver_omp(occ, comm)


    Real( Kind = wp ) :: Totstart, occ, dphi
    Type(comms_type), Intent( InOut)   :: comm

    Integer           :: mmm, i,j,k
    Real( Kind = wp ) :: element,  sm1,sm2,sm3,sm4 ! SM stands for Stoyan Markov (long live!!!)

    Totstart=MPI_WTIME()
    pconverged=.false.
    normphi0=0.0_wp

    !$omp paralleldo default(shared) private(k) reduction(+:normphi0)
    Do k=1,block_z
       normphi0 = normphi0 + Sum(phi(1:block_x,1:block_y,k)**2)
    End Do
    !$omp End paralleldo

    Call gsum(comm,normphi0)

    Do mmm=1,mxitjb+3 ! Za seki sluchaj, proverka sa stabilno shozhdane

       Call biCGStab_exchange_halo(phi,xhalo,comm)

       sm1 = -600.0_wp/144.0_wp
       sm2 =   60.0_wp/144.0_wp
       sm3 =   18.0_wp/144.0_wp
       sm4 =    3.0_wp/144.0_wp

       !$omp paralleldo default(shared) private(element,k,j,i)
       Do k=lnzl+1,lnzu-1
          Do j=lnyl+1,lnyu-1
             Do i=lnxl+1,lnxu-1
                element=0.0_wp
                element=element+sm4*phi(i - 1 ,j - 1 ,k - 1 )
                element=element+sm3*phi(i     ,j - 1 ,k - 1 )
                element=element+sm4*phi(i + 1 ,j - 1 ,k - 1 )
                element=element+sm3*phi(i - 1 ,j     ,k - 1 )
                element=element+sm2*phi(i     ,j     ,k - 1 )
                element=element+sm3*phi(i + 1 ,j     ,k - 1 )
                element=element+sm4*phi(i - 1 ,j + 1 ,k - 1 )
                element=element+sm3*phi(i     ,j + 1 ,k - 1 )
                element=element+sm4*phi(i + 1 ,j + 1 ,k - 1 )
                element=element+sm3*phi(i - 1 ,j - 1 ,k     )
                element=element+sm2*phi(i     ,j - 1 ,k     )
                element=element+sm3*phi(i + 1 ,j - 1 ,k     )
                element=element+sm2*phi(i - 1 ,j     ,k     )
                element=element+sm2*phi(i + 1 ,j     ,k     )
                element=element+sm3*phi(i - 1 ,j + 1 ,k     )
                element=element+sm2*phi(i     ,j + 1 ,k     )
                element=element+sm3*phi(i + 1 ,j + 1 ,k     )
                element=element+sm4*phi(i - 1 ,j - 1 ,k + 1 )
                element=element+sm3*phi(i     ,j - 1 ,k + 1 )
                element=element+sm4*phi(i + 1 ,j - 1 ,k + 1 )
                element=element+sm3*phi(i - 1 ,j     ,k + 1 )
                element=element+sm2*phi(i     ,j     ,k + 1 )
                element=element+sm3*phi(i + 1 ,j     ,k + 1 )
                element=element+sm4*phi(i - 1 ,j + 1 ,k + 1 )
                element=element+sm3*phi(i     ,j + 1 ,k + 1 )
                element=element+sm4*phi(i + 1 ,j + 1 ,k + 1 )
                phi0(i,j,k)=(-b(i,j,k)+element)/(-sm1)
             End Do
          End Do
       End Do
       !$omp End paralleldo

       normphi1=0.0_wp
       !$omp paralleldo default(shared) private(k) reduction(+:normphi1)
       Do k=1,block_z
          phi(:,:,k)=phi0(:,:,k)
          normphi1=normphi1+Sum(phi(1:block_x,1:block_y,k)**2)
       End Do
       !$omp End paralleldo

       Call gsum(comm, normphi1)

       dphi=Abs(normphi1-normphi0)
       If (dphi <= occ*normb) Then
          Call biCGStab_exchange_halo(phi,xhalo,comm)
          If (comm%idnode==0) Write (*,*)  "Jacobi it = ", mmm ,&
               &"d|Coulomb potential|/NormB >>>",&
               &Abs(normphi1 - normphi0)/normb, MPI_WTIME()-Totstart
          pconverged=.true.
          Return
       End If

       If (mmm > mxitjb) Call biCGStab_solver_omp(occ,comm)
       normphi0=normphi1

    End Do

    Print*,"poisson solver not converged in ",maxsteps,"steps..."
    pconverged=.false.

  End Subroutine P_solver_omp

  Recursive Subroutine biCGStab_solver_omp(occ,comm)

    Type(comms_type), Intent( InOut )  :: comm
    Real( Kind = wp ) :: alfa, beta, omega, rho1,rho0,rv,tt,ts
    Real( Kind = wp ) :: Totstart,occ

    Integer :: mmm, kk

    Totstart=MPI_WTIME()

    alfa=1.0_wp
    omega=1.0_wp
    rho1=1.0_wp
    rho0=1.0_wp

    Call biCGStab_exchange_halo(phi,xhalo,comm)
    Call Adot_omp(phi,xhalo,F,0)

    !$omp parallel default(shared) private(kk)
    !$omp paralleldo
    Do kk=1,block_z
       r0(1:block_x,1:block_y,kk) = b(1:block_x,1:block_y,kk) - f(1:block_x,1:block_y,kk)

       r(1:block_x,1:block_y,kk)  = r0(1:block_x,1:block_y,kk)
    End Do
    !$omp End paralleldo

    !$omp paralleldo
    Do kk=lnzl,lnzu
       p(:,:,kk) = 0.0_wp
       v(:,:,kk) = 0.0_wp
    End Do
    !$omp End paralleldo

    normphi0=0.0_wp
    !$omp paralleldo reduction(+:normphi0)
    Do kk=1,block_z
       normphi0 = normphi0 + Sum(phi(1:block_x,1:block_y,kk)**2)
    End Do
    !$omp End parallel

    Call gsum(comm,normphi0)
    Do mmm=1,mxitcg
       rho0 = rho1
       rho1 = 0.0_wp
       !$omp paralleldo reduction(+:rho1)
       Do kk=1,block_z
          rho1 = rho1 + Sum(r0(1:block_x,1:block_y,kk)*r(1:block_x,1:block_y,kk))
       End Do
       !$omp End paralleldo
       Call gsum(comm,rho1)

       beta = (rho1/rho0) * (alfa/omega)

       !$omp paralleldo default(shared) private(kk)
       Do kk=1,block_z
          p(1:block_x,1:block_y,kk) = r(1:block_x,1:block_y,kk) + &
                                      beta * (p(1:block_x,1:block_y,kk) - omega * v(1:block_x,1:block_y,kk))
       End Do
       !$omp End paralleldo

       Call biCGStab_exchange_halo(p,0,comm)
       Call Adot_omp(p,0,v,0)

       rv=0.0_wp
       !$omp parallel Do default(shared) reduction(+:rv) private(kk)
       Do kk=1,block_z
          rv = rv + Sum(r0(1:block_x,1:block_y,kk) * v(1:block_x,1:block_y,kk))
       End Do
       !$omp End paralleldo
       Call gsum(comm,rv)

       alfa=rho1/rv

       !$omp paralleldo default(shared) private(kk)
       Do kk=1,block_z
          s(1:block_x,1:block_y,kk) = r(1:block_x,1:block_y,kk) - alfa * v(1:block_x,1:block_y,kk)
       End Do
       !$omp End paralleldo

       Call biCGStab_exchange_halo(s,0,comm)
       Call Adot_omp(s,0,t,0)

       tt=0.0_wp
       !$omp paralleldo default(shared) private(kk) reduction(+:tt)
       Do kk=1,block_z
          tt = tt + Sum(t(1:block_x,1:block_y,kk)**2)
       End Do
       !$omp End paralleldo
       Call gsum(comm,tt)

       ts=0.0_wp
       !$omp paralleldo default(shared) private(kk) reduction(+:ts)
       Do kk=1,block_z
          ts = ts + Sum(t(1:block_x,1:block_y,kk)*s(1:block_x,1:block_y,kk))
       End Do
       !$omp End paralleldo
       Call gsum(comm,ts)

       omega=ts/tt

       !$omp paralleldo default(shared) private(kk)
       Do kk=1,block_z
          phi(1:block_x,1:block_y,kk) = phi(1:block_x,1:block_y,kk) + &
               alfa * p(1:block_x,1:block_y,kk) + omega * s(1:block_x,1:block_y,kk)
       End Do
       !$omp End paralleldo

       normphi1=0.0_wp
       !$omp paralleldo default(shared) private(kk) reduction(+:normphi1)
       Do kk=1,block_z
          normphi1 = normphi1 + Sum(phi(1:block_x,1:block_y,kk)**2)
       End Do
       !$omp End paralleldo
       Call gsum(comm,normphi1)

       If (maxbicgst > 0 .and. mmm > maxbicgst) Then
          !$omp paralleldo default(shared) private(kk)
          Do kk=lnzl,lnzu
             phi(:,:,kk) = 0.0_wp
          End Do
          !$omp End paralleldo

          maxbicgst=0
          pconverged=.false.
          If (comm%idnode==0) &
             Print*, "bicgstab exceeded *maxbicgst*... reset potential... pconverged=.false."
          Return
       End If
       Write (*,*)  "biCGStab it = ", mmm ,&
            &"d|Coulomb potential|/NormB >>>",(Abs(normphi1 - normphi0)/normb),&
            &MPI_WTIME()-Totstart

       If (Abs(normphi1 - normphi0) <= occ*normb) Then
          Call biCGStab_exchange_halo(phi,xhalo,comm)

         If (comm%idnode == 0) Write (*,*)  "biCGStab it = ", mmm ,&
              &"d|Coulomb potential|/NormB >>>",(Abs(normphi1 - normphi0)/normb),&
              &MPI_WTIME()-Totstart
          If (maxbicgst == 0) Then
             maxbicgst=mmm
          End If
          pconverged=.true.
          If (comm%idnode == 0) Print*, "*maxbicgst* set to",maxbicgst,"pconverged=.true."
          Return
       End If

       !$omp paralleldo default(shared) private(kk)
       Do kk=1,block_z
          r(1:block_x,1:block_y,kk) = s(1:block_x,1:block_y,kk) - omega * t(1:block_x,1:block_y,kk)
       End Do
       !$omp End paralleldo

       normphi0=normphi1
    End Do
    pconverged=.false.
    If (comm%idnode == 0) Print*, "maxit reached... pconverged=.false."

  End Subroutine biCGStab_solver_omp

  Subroutine biCGStab_exchange_halo(vec,xtra,comm)

    Integer,           Intent( In    ) :: xtra

    Real( Kind = wp ), Intent( InOut ) :: vec( lnxl-xtra:lnxu+xtra, &
                                               lnyl-xtra:lnyu+xtra, &
                                               lnzl-xtra:lnzu+xtra )
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: me, lx,ly,lz

! What's my name?

    me = comm%idnode

! Find length of sides of the domain

    lx = ixt - ixb + 1
    ly = iyt - iyb + 1
    lz = izt - izb + 1

! +X direction face - negative halo

    Call exchange_grid_halo( lmap(1),                     lmap(2), &
         xtra+1,                      ly,                          lz, &
         lx-(xtra+1)+1, lx  ,         1,              ly,            1,              lz, &
         1-(xtra+1),    1-1,          1,              ly,            1,              lz,comm)

! -X direction face - positive halo

    Call exchange_grid_halo( lmap(2),                     lmap(1), &
         xtra+1,                      ly,                          lz, &
         1,             1+(xtra+1)-1, 1,              ly,            1,              lz, &
         lx+1,          lx+(xtra+1),  1,              ly,            1,              lz,comm)

! +Y direction face (including the +&-X faces extensions) - negative halo

    Call exchange_grid_halo( lmap(3),                     lmap(4), &
         lx+2*(xtra+1),              xtra+1,                       lz, &
         1-(xtra+1),    lx+(xtra+1),  ly-(xtra+1)+1, ly,             1,              lz, &
         1-(xtra+1),    lx+(xtra+1),  1-(xtra+1)  ,  1-1,            1,              lz,comm)

! -Y direction face (including the +&-X faces extensions) - positive halo

    Call exchange_grid_halo( lmap(4),                     lmap(3), &
         lx+2*(xtra+1),              xtra+1,                       lz, &
         1-(xtra+1),    lx+(xtra+1),  1,              1+(xtra+1)-1,  1, lz, &
         1-(xtra+1),    lx+(xtra+1),  ly+1,           ly+(xtra+1),   1, lz,comm)

! +Z direction face (including the +&-Y+&-X faces extensions) - negative halo

    Call exchange_grid_halo( lmap(5),                     lmap(6), &
         lx+2*(xtra+1),              ly+2*(xtra+1),              xtra+1, &
         1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  lz-(xtra+1)+1, lz, &
         1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  1-(xtra+1)  ,  1-1,comm)

! -Z direction face (including the +&-Y+&-X faces extensions) - positive halo

    Call exchange_grid_halo( lmap(6),                     lmap(5), &
         lx+2*(xtra+1),              ly+2*(xtra+1),              xtra+1, &
         1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  1,              1+(xtra+1)-1, &
         1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  lz+1,           lz+(xtra+1),comm)

  Contains

    Subroutine exchange_grid_halo(     from,       to,           &
                                   lx,       ly,       lz,       &
                                   xlb, xlt, ylb, ylt, zlb, zlt, &
                                   xdb, xdt, ydb, ydt, zdb, zdt,comm )

      Integer, Intent( In    ) :: from, to
      Integer, Intent( In    ) :: lx, ly, lz
      Integer, Intent( In    ) :: xlb, ylb, zlb
      Integer, Intent( In    ) :: xlt, ylt, zlt
      Integer, Intent( In    ) :: xdb, ydb, zdb
      Integer, Intent( In    ) :: xdt, ydt, zdt
      Type( comms_type), Intent( InOut ) :: comm

      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: send_buffer
      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: recv_buffer

      Integer :: length

! If the processor to receive FROM is actually ME it means there is
! only one processor along this axis (so the processor to send TO is
! also ME) and so no message passing need be done.  However, there
! is a catch.  The domain decomposed spme_forces routine expects data
! that would be `seen' through the periodic boundary conditions to be
! actually copied from the high positive indices to negative ones and
! vice-versa.  The Else clause catches this.

      If ( from /= me ) Then

! Length of message to send is the same as that to receive

         length = lx * ly * lz

! Allocate send and receive buffers (of the same size!!!)
! so all can be sent and received as one message!!!

         If (from > -1) Then
            Allocate ( recv_buffer( xdb:xdt, ydb:ydt, zdb:zdt ) , Stat = fail(1) )
            If (fail(1) > 0) Then
               Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo receive allocation failure, node: ', comm%idnode
               Call error(0)
            End If

            Call MPI_IRECV( recv_buffer, length, wp_mpi, from, ExchgGrid_tag, comm%comm, comm%request, comm%ierr )
         End If

         If (to   > -1) Then
            Allocate ( send_buffer( xlb:xlt, ylb:ylt, zlb:zlt ) , Stat = fail(1) )
            If (fail(1) > 0) Then
               Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo send allocation failure, node: ', comm%idnode
               Call error(0)
            End If

! Copy the data to be sent

            send_buffer = vec( xlb:xlt, ylb:ylt, zlb:zlt )

            Call MPI_SEND(  send_buffer, length, wp_mpi, to  , ExchgGrid_tag, comm%comm, comm%ierr )
         End If

! Exchange the data

         If (from > -1) Then
            Call MPI_WAIT(  comm%request, comm%status, comm%ierr )

! Copy the received data into the domain halo

            vec( xdb:xdt, ydb:ydt, zdb:zdt ) = recv_buffer

! And, as my mum told me, leave things as we found them

            Deallocate ( recv_buffer , Stat = fail(1) )
            If (fail(1) > 0) Then
               Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo receive deallocation failure, node: ', comm%idnode
               Call error(0)
            End If
         End If

         If (to   > -1) Then
            Deallocate ( send_buffer , Stat = fail(1) )
            If (fail(1) > 0) Then
               Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo send deallocation failure, node: ', comm%idnode
               Call error(0)
            End If
         End If

      Else

! Simple on node copy - as sizes are the same as 3D shapes

         vec( xdb:xdt, ydb:ydt, zdb:zdt ) = vec( xlb:xlt, ylb:ylt, zlb:zlt )

      End If

    End Subroutine exchange_grid_halo

  End Subroutine biCGStab_exchange_halo

  Subroutine biCGStab_Deallocate_grids(comm)

    Type( comms_type ), Intent( InOut ) :: comm
    Deallocate(t ,   b ,  v , Stat = fail(1) )
    Deallocate(F ,   s ,  r , Stat = fail(2) )
    Deallocate(p , phi , r0 , Stat = fail(3) )
    If (Any(fail(1:3) > 0)) Then
       Write(nrite,'(/,1x,a,i0)') 'exchange_grid_halo send deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

  End Subroutine biCGStab_Deallocate_grids

  Subroutine Adot_omp(vec,vx,res,rx)

    Integer          , Intent( In    ) :: vx,rx
    Real( Kind = wp ), Intent( In    ) :: vec(lnxl-vx:lnxu+vx,lnyl-vx:lnyu+vx,lnzl-vx:lnzu+vx)
    Real( Kind = wp ), Intent(   Out ) :: res(lnxl-rx:lnxu+rx,lnyl-rx:lnyu+rx,lnzl-rx:lnzu+rx)

    Integer          :: i,j,k,ioff,joff,koff
    Real( Kind = wp) :: sm(4), element

    sm(1) = -600.0_wp/144.0_wp
    sm(2) =   60.0_wp/144.0_wp
    sm(3) =   18.0_wp/144.0_wp
    sm(4) =    3.0_wp/144.0_wp

    !$omp paralleldo default(shared) private(element,k,j,i,ioff,joff,koff)
    Do k=lnzl+1,lnzu-1
       Do j=lnyl+1,lnyu-1
          Do i=lnxl+1,lnxu-1
             element=0.0
             Do koff=-1,1
                Do joff=-1,1
                   Do ioff=-1,1
                      element=element+&
                           sm(ioff**2+joff**2+koff**2+1)*&
                           vec(i+ioff,j+joff,k+koff)
                   End Do
                End Do
             End Do
             res(i,j,k)=element
          End Do
       End Do
    End Do
    !$omp End paralleldo

  End Subroutine Adot_omp

  Subroutine biCGStab_calc_forces(cenergy,vir,stress,comm)

    Real( Kind = wp ), Intent( InOut ) :: cenergy, vir, stress(1:9)
    Type( comms_type), Intent( InOut ) :: comm

    Integer       :: i,j,k,n, ii,jj,kk
    Real(Kind=wp) :: reps0dv,txx,tyy,tzz, det,rcell(1:9),        &
                     cfxx, cfyy, cfzz,                           &
                     uenergy,r8veps0

    If(.Not.pconverged) Then
       if (comm%idnode == 0) Print*, "poisson solver not converged"
       Return
    End If

    reps0dv=fourpi*r4pie0/delta
    r8veps0=fourpi*r4pie0/(8.0_wp*delta**3)

! initialization

    cenergy=0.0_wp
    vir=0.0_wp
    stress=0.0_wp

    Call invert(cell,rcell,det)
    If (Abs(det) < 1.0e-6_wp) Call error(120)

    !$omp paralleldo default(shared) private(n,txx,tyy,tzz,ii,jj,kk,i,j,k,cfxx,cfyy,cfzz) &
    !$omp reduction(+:stress,vir,cenergy)
    Do n=1,natms
       If (Abs(chge(n)) < Abs(half_minus)) Cycle

       txx=Real(kmaxa,wp)*(rcell(1)*xxx(n)+rcell(4)*yyy(n)+rcell(7)*zzz(n)+half_minus)
       tyy=Real(kmaxb,wp)*(rcell(2)*xxx(n)+rcell(5)*yyy(n)+rcell(8)*zzz(n)+half_minus)
       tzz=Real(kmaxc,wp)*(rcell(3)*xxx(n)+rcell(6)*yyy(n)+rcell(9)*zzz(n)+half_minus)

! global indeces

       ii=Int(txx)
       jj=Int(tyy)
       kk=Int(tzz)

! get local indces

       i = ii - ixb + 2
       j = jj - iyb + 2
       k = kk - izb + 2

       cfxx=chge(n)*dPhidX(i,j,k)
       cfyy=chge(n)*dPhidY(i,j,k)
       cfzz=chge(n)*dPhidZ(i,j,k)

! calculate atomic energy

       uenergy=chge(n)*phi(i,j,k)

!caclulate atomic contribution to the stress tensor

       stress(1)=stress(1)+xxx(n)*cfxx
       stress(2)=stress(2)+xxx(n)*cfyy
       stress(3)=stress(3)+xxx(n)*cfzz
       stress(4)=stress(4)+yyy(n)*cfxx
       stress(5)=stress(5)+yyy(n)*cfyy
       stress(6)=stress(6)+yyy(n)*cfzz
       stress(7)=stress(7)+zzz(n)*cfxx
       stress(8)=stress(8)+zzz(n)*cfyy
       stress(9)=stress(9)+zzz(n)*cfzz

       fxx(n)=fxx(n)+cfxx
       fyy(n)=fyy(n)+cfyy
       fzz(n)=fzz(n)+cfzz

!add enegry to the total energy

       cenergy=cenergy+uenergy

       vir=vir-uenergy
    End Do
    !$omp End paralleldo

  End Subroutine biCGStab_calc_forces

  Subroutine Write_potential
    Integer :: unit, kpkp,kk,kpk,tt
    character(len=128) :: filename, line
    unit=39
    line="pot"
    Write(filename,'(a,i1,i1,i1,a)') Trim(Adjustl(line)),idx,idy,idz,".dx"
    Open(Unit=unit, File=Trim(Adjustl(filename)), Status='replace')
    Write(unit,'(a,1x,I7,1x,I7,1x,I7)') "object 1 class gridpositions counts ", block_x,block_y,block_z
    Write(unit,'(a,1x,f9.3,1x,f9.3,1x,f9.3)') "origin ", &
         Real(idx*(kmaxa/nprx),wp)*delta,Real(idy*(kmaxb/npry),wp)*delta,Real(idz*(kmaxc/nprz),wp)*delta
    Write(unit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "delta ", (/delta, 0.0_wp,0.0_wp/)
    Write(unit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "delta ", (/0.0_wp,delta, 0.0_wp/)
    Write(unit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "delta ", (/0.0_wp,0.0_wp, delta/)
    Write(unit,'(a,I5,1x,I5,1x,I5)') "object 2 class gridconnections counts ",block_x,block_y,block_z
    Write(unit,'(a,1x,I20,1x,a)')'object 3 class array type Double rank 0 items ', &
         &block_x*block_y*block_z , 'data follows'

    tt=0
    line=""
    Do kpkp=1,block_x
       Do kk=1,block_y
          Do kpk=1,block_z
             Write(line,"(a,3E15.5)") Trim(line), Real(phi(kpkp,kk,kpk),4)
             tt=tt+1
             If ( mod(tt,3) == 0 ) Then
                Write(unit,'(a)') Trim(Adjustl(line))
                line=""
                tt=0
             End If
          End Do
       End Do
    End Do

    Write(unit,'(a)') Trim(Adjustl(line))
    Write(unit,'(a)') 'attribute "dep" string "positions"'
    Write(unit,'(a)') 'object "PME potential (kT/e, T=300K)" class field'
    Write(unit,'(a)') 'component "positions" value 1'
    Write(unit,'(a)') 'component "connections" value 2'
    Write(unit,'(a)') 'component "data" value 3'

    close(unit=unit)

  End Subroutine Write_potential

  Subroutine Write_b
    Integer :: unit, kpkp,kk,kpk,tt
    character(len=128) :: filename, line
    unit=39
    line="b"
    Write(filename,'(a,i1,i1,i1,a)') Trim(Adjustl(line)),idx,idy,idz,".dx"
    Open(Unit=unit, File=Trim(Adjustl(filename)), Status='replace')
    Write(unit,'(a,1x,I7,1x,I7,1x,I7)') "object 1 class gridpositions counts ", block_x,block_y,block_z
    Write(unit,'(a,1x,f9.3,1x,f9.3,1x,f9.3)') "origin ", &
         Real(idx*(kmaxa/nprx),wp)*delta,Real(idy*(kmaxb/npry),wp)*delta,Real(idz*(kmaxc/nprz),wp)*delta
    Write(unit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "delta ", (/Real(delta), 0.0,0.0/)
    Write(unit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "delta ", (/ 0.0,Real(delta), 0.0/)
    Write(unit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "delta ", (/0.0,0.0, Real(delta)/)
    Write(unit,'(a,I5,1x,I5,1x,I5)') "object 2 class gridconnections counts ",block_x,block_y,block_z
    Write(unit,'(a,1x,I20,1x,a)')'object 3 class array type Double rank 0 items ', &
         &block_x*block_y*block_z , 'data follows'

    tt=0
    line=""
    Do kpkp=1,block_x
       Do kk=1,block_y
          Do kpk=1,block_z
             Write(line,"(a,3E15.5)") Trim(line), Real(b(kpkp,kk,kpk),4)
             tt=tt+1
             If ( mod(tt,3) == 0 ) Then
                Write(unit,'(a)') Trim(Adjustl(line))
                line=""
                tt=0
             End If
          End Do
       End Do
    End Do

    Write(unit,'(a)') Trim(Adjustl(line))
    Write(unit,'(a)') 'attribute "dep" string "positions"'
    Write(unit,'(a)') 'object "PME potential (kT/e, T=300K)" class field'
    Write(unit,'(a)') 'component "positions" value 1'
    Write(unit,'(a)') 'component "connections" value 2'
    Write(unit,'(a)') 'component "data" value 3'

    close(unit=unit)

  End Subroutine Write_b

  Function dPhidX(i,j,k) Result(dphi)

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Integer           :: kk
    Real( Kind = wp ) :: par(3), par4, dphi

    par(1)=45.0_wp
    par(2)=-9.0_wp
    par(3)=1.0_wp
    par4=60.0_wp

    ! sm  du/dx = [45(ui+1 – ui-1 ) -9(ui+2 – ui-2) + (ui+3 –ui-3)] / 60h
    dphi=0.0_wp
    Do kk=1,3
       dphi = dphi + par(kk)*(phi(i+kk,j,k)-phi(i-kk,j,k))
    End Do
    dphi = dphi/(par4*delta)

  End Function dPhidX

  Function dPhidY(i,j,k) Result(dphi)

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Integer           :: kk
    Real( Kind = wp ) :: par(3), par4, dphi

    par(1)=45.0_wp
    par(2)=-9.0_wp
    par(3)=1.0_wp
    par4=60.0_wp

    ! sm  du/dx = [45(ui+1 – ui-1 ) -9(ui+2 – ui-2) + (ui+3 –ui-3)] / 60h
    dphi=0.0_wp
    Do kk=1,3
       dphi = dphi + par(kk)*(phi(i,j+kk,k)-phi(i,j-kk,k))
    End Do
    dphi = dphi/(par4*delta)

  End Function dPhidY

  Function dPhidZ(i,j,k) Result(dphi)

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Integer           :: kk
    Real( Kind = wp ) :: par(3), par4, dphi

    par(1)=45.0_wp
    par(2)=-9.0_wp
    par(3)=1.0_wp
    par4=60.0_wp

    ! sm  du/dx = [45(ui+1 – ui-1 ) -9(ui+2 – ui-2) + (ui+3 –ui-3)] / 60h
    dphi=0.0_wp
    Do kk=1,3
       dphi = dphi + par(kk)*(phi(i,j,k+kk)-phi(i,j,k-kk))
    End Do
    dphi = dphi/(par4*delta)

  End Function dPhidZ

  Function d2PhidX2(i,j,k) Result(dphi)

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Real( Kind = wp ) :: dphi

    dphi=( phi(i+2,j,k) + phi(i-1,j,k) - phi(i+1,j,k) - phi(i-2,j,k) ) / (3.0_wp*delta**2)

  End Function d2PhidX2

  Function d2PhidY2(i,j,k) Result(dphi)

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Real( Kind = wp ) :: dphi

    dphi=( phi(i,j+2,k) + phi(i,j-1,k) - phi(i,j+1,k) - phi(i,j-2,k) ) / (3.0_wp*delta**2)

  End Function d2PhidY2

  Function d2PhidZ2(i,j,k) Result(dphi)

    Implicit None

    Integer, Intent( In    ) :: i,j,k

    Real( Kind = wp ) :: dphi

    dphi=( phi(i,j,k+2) + phi(i,j,k-1) - phi(i,j,k+1) - phi(i,j,k-2) ) / (3.0_wp*delta**2)

  End Function d2PhidZ2

  Subroutine poisson_excl_forces &
             (iatm,rcut,epsq,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the exclusion coulombic energy
! and! force terms in a periodic system using 1/r potential as required
! for a direct space poisson solver
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: limit,idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,coul,fcoul, &
                         fix,fiy,fiz,fx,fy,fz,        &
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

! global identity and type of of iatm

    idi=ltg(iatm)
    chgea = chge(iatm)

! ignore interaction if the charge is zero

    If (Abs(chgea) > zero_plus) Then

       chgea = chgea*r4pie0/epsq

! load forces

       fix=fxx(iatm)
       fiy=fyy(iatm)
       fiz=fzz(iatm)

! Get limit

       limit=list(-1,iatm)-list(0,iatm)

! start of primary loop for forces evaluation

       Do m=1,limit

! atomic index and charge

          jatm=list(m,iatm)
          chgprd=chge(jatm)

! interatomic distance

          rrr=rrt(m)

! interaction validity and truncation of potential

          If (Abs(chgprd) > zero_plus .and. rrr < rcut) Then

! charge product

             chgprd=-chgprd*chgea ! the MINUS signifies the exclusion!!!

! calculate forces

             coul = chgprd/rrr
             fcoul = coul/rrr**2

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

! calculate potential energy

                engcpe_ex = engcpe_ex + coul

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

! virial

       vircpe_ex = -engcpe_ex

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

  End Subroutine poisson_excl_forces

  Subroutine poisson_frzn_forces(rcut,epsq,engcpe_fr,vircpe_fr,stress,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating corrections to coulombic energy
! and forces in a periodic system arising from frozen pairs
! for a direct space poisson solver
!
! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
!       end (and any COM drift removed) but corrections to the stress
!       and the virial are important as they feed into the system
!       pressure response.  Constant volume ensembles (keyens < 20)
!       need this calculation just once (NOT! - controlled by lf_fce in
!       ewald_check<-two_body_forces
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),                   Intent( In    ) :: rcut,epsq
    Real( Kind = wp ),                   Intent(   Out ) :: engcpe_fr,vircpe_fr
    Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
    Type(comms_type), Intent( InOut)                     :: comm

    Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
    Real( Kind = wp ) :: det,rcell(1:9),xrr,yrr,zrr,rrr,rsq, &
                         chgprd,coul,fcoul,                  &
                         fx,fy,fz,xss,yss,zss,               &
                         strs1,strs2,strs3,strs5,strs6,strs9

    Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr
    Real( Kind = wp ), Dimension( : ), Allocatable :: cfr,xfr,yfr,zfr
    Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt,rrt

    If (.not.lf_fce) Then ! All's been done but needs copying
       Do i=1,natms
          fxx(i)=fxx(i)+ffx(i)
          fyy(i)=fyy(i)+ffy(i)
          fzz(i)=fzz(i)+ffz(i)
       End Do

       engcpe_fr=ef_fr
       vircpe_fr=vf_fr
       stress=stress+sf_fr

       If (l_cp) Then
          Do i=1,natms
             fcx(i)=fcx(i)+ffx(i)
             fcy(i)=fcy(i)+ffy(i)
             fcz(i)=fcz(i)+ffz(i)
          End Do

          e_fr=ef_fr
          v_fr=vf_fr
          s_fr=sf_fr
       End If

       Return
    End If


    fail=0
    Allocate (l_ind(1:mxatdm),nz_fr(0:comm%mxnode), Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'poisson_frzn_forces allocation failure, node: ', comm%idnode
       Call error(0)
    End If

    Call invert(cell,rcell,det)

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
    Do i=1,natms
       If (lfrzn(i) > 0 .and. Abs(chge(i)) > zero_plus) Then
          nz_fr(comm%idnode+1)=nz_fr(comm%idnode+1)+1
          l_ind(nz_fr(comm%idnode+1))=i
       End If
    End Do
    Call gsum(comm,nz_fr)
    nz_fr(0) = Sum(nz_fr(0:comm%idnode)) ! Offset

    nzfr = Sum(nz_fr(1:comm%mxnode))     ! Total
    If (nzfr <= 10*mxatms) Then

       Allocate (cfr(1:nzfr),xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
       If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'poisson_frzn_forces allocation failure 1, node: ', comm%idnode
          Call error(0)
       End If

       cfr=0.0_wp
       xfr=0.0_wp
       yfr=0.0_wp
       zfr=0.0_wp
       Do i=1,nz_fr(comm%idnode+1)
          ii=nz_fr(0)+i

          cfr(ii)=chge(l_ind(i))
          xfr(ii)=xxx(l_ind(i))
          yfr(ii)=yyy(l_ind(i))
          zfr(ii)=zzz(l_ind(i))
       End Do
       If (comm%mxnode > 1) Then
          Call gsum(comm,cfr)
          Call gsum(comm,xfr)
          Call gsum(comm,yfr)
          Call gsum(comm,zfr)
       End If

       Do i=1,nz_fr(comm%idnode+1)
          ii=nz_fr(0)+i

          Do jj=1,nz_fr(0) ! -, on nodes<idnode
             xrr=xfr(ii)-xfr(jj)
             yrr=yfr(ii)-yfr(jj)
             zrr=zfr(ii)-zfr(jj)

             xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
             yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
             zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

             xss=xss-Anint(xss)
             yss=yss-Anint(yss)
             zss=zss-Anint(zss)

             xrr=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
             yrr=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
             zrr=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

! calculate interatomic distance

             rsq=xrr**2+yrr**2+zrr**2

             rrr=Sqrt(rsq)
             chgprd=-cfr(ii)*cfr(jj)/epsq*r4pie0 ! the MINUS signifies the exclusion!!!

! calculate forces

             coul = chgprd/rrr
             fcoul = coul/rsq

             fx = fcoul*xrr
             fy = fcoul*yrr
             fz = fcoul*zrr

! calculate forces

             fxx(l_ind(i))=fxx(l_ind(i))-fx
             fyy(l_ind(i))=fyy(l_ind(i))-fy
             fzz(l_ind(i))=fzz(l_ind(i))-fz

! redundant calculations copying

             If (lf_cp) Then
                ffx(l_ind(i))=ffx(l_ind(i))-fx
                ffy(l_ind(i))=ffy(l_ind(i))-fy
                ffz(l_ind(i))=ffz(l_ind(i))-fz
             End If

! infrequent calculations copying

             If (l_cp) Then
                fcx(l_ind(i))=fcx(l_ind(i))-fx
                fcy(l_ind(i))=fcy(l_ind(i))-fy
                fcz(l_ind(i))=fcz(l_ind(i))-fz
             End If
          End Do

          Do j=i+1,nz_fr(comm%idnode+1) ! =, node=idnode (OVERLAP but no SELF)!
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

             xrr=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
             yrr=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
             zrr=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

! calculate interatomic distance

             rsq=xrr**2+yrr**2+zrr**2

             rrr=Sqrt(rsq)
             chgprd=-cfr(ii)*cfr(jj)/epsq*r4pie0 ! the MINUS signifies the exclusion!!!

! calculate forces

             coul = chgprd/rrr
             fcoul = coul/rsq

             fx = fcoul*xrr
             fy = fcoul*yrr
             fz = fcoul*zrr

! calculate forces

             fxx(l_ind(i))=fxx(l_ind(i))-fx
             fyy(l_ind(i))=fyy(l_ind(i))-fy
             fzz(l_ind(i))=fzz(l_ind(i))-fz

             fxx(l_ind(j))=fxx(l_ind(j))+fx
             fyy(l_ind(j))=fyy(l_ind(j))+fy
             fzz(l_ind(j))=fzz(l_ind(j))+fz

! redundant calculations copying

             If (lf_cp) Then
                ffx(l_ind(i))=ffx(l_ind(i))-fx
                ffy(l_ind(i))=ffy(l_ind(i))-fy
                ffz(l_ind(i))=ffz(l_ind(i))-fz

                ffx(l_ind(j))=ffx(l_ind(j))+fx
                ffy(l_ind(j))=ffy(l_ind(j))+fy
                ffz(l_ind(j))=ffz(l_ind(j))+fz
             End If

! infrequent calculations copying

             If (l_cp) Then
                fcx(l_ind(i))=fcx(l_ind(i))-fx
                fcy(l_ind(i))=fcy(l_ind(i))-fy
                fcz(l_ind(i))=fcz(l_ind(i))-fz

                fcx(l_ind(j))=fcx(l_ind(j))+fx
                fcy(l_ind(j))=fcy(l_ind(j))+fy
                fcz(l_ind(j))=fcz(l_ind(j))+fz
             End If

! calculate potential energy

             engcpe_fr = engcpe_fr + coul

! calculate stress tensor

             strs1 = strs1 + xrr*fx
             strs2 = strs2 + xrr*fy
             strs3 = strs3 + xrr*fz
             strs5 = strs5 + yrr*fy
             strs6 = strs6 + yrr*fz
             strs9 = strs9 + zrr*fz
          End Do

          Do jj=nz_fr(0)+nz_fr(comm%idnode+1)+1,nzfr ! +, on nodes>idnode
             xrr=xfr(ii)-xfr(jj)
             yrr=yfr(ii)-yfr(jj)
             zrr=zfr(ii)-zfr(jj)

             xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
             yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
             zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

             xss=xss-Anint(xss)
             yss=yss-Anint(yss)
             zss=zss-Anint(zss)

             xrr=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
             yrr=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
             zrr=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

! calculate interatomic distance

             rsq=xrr**2+yrr**2+zrr**2

             rrr=Sqrt(rsq)
             chgprd=-cfr(ii)*cfr(jj)/epsq*r4pie0 ! the MINUS signifies the exclusion!!!

! calculate forces

             coul = chgprd/rrr
             fcoul = coul/rsq

             fx = fcoul*xrr
             fy = fcoul*yrr
             fz = fcoul*zrr

             fxx(l_ind(i))=fxx(l_ind(i))-fx
             fyy(l_ind(i))=fyy(l_ind(i))-fy
             fzz(l_ind(i))=fzz(l_ind(i))-fz

! redundant calculations copying

             If (lf_cp) Then
                ffx(l_ind(i))=ffx(l_ind(i))-fx
                ffy(l_ind(i))=ffy(l_ind(i))-fy
                ffz(l_ind(i))=ffz(l_ind(i))-fz
             End If

! infrequent calculations copying

             If (l_cp) Then
                fcx(l_ind(i))=fcx(l_ind(i))-fx
                fcy(l_ind(i))=fcy(l_ind(i))-fy
                fcz(l_ind(i))=fcz(l_ind(i))-fz
             End If

! calculate potential energy

             engcpe_fr = engcpe_fr + coul

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
          Write(nrite,'(/,1x,a,i0)') 'poisson_frzn_forces deallocation failure 1, node: ', comm%idnode
          Call error(0)
       End If

    Else

! We resort to approximating N*(N-1)/2 interactions
! with the short-range one from the two body linked cell list

       Allocate (xxt(1:mxlist),yyt(1:mxlist),zzt(1:mxlist),rrt(1:mxlist), Stat=fail)
       If (fail > 0) Then
          Write(nrite,'(/,1x,a,i0)') 'poisson_frzn_forces allocation failure 2, node: ', comm%idnode
          Call error(0)
       End If

       Do ii=1,nz_fr(comm%idnode+1)
          i=l_ind(nz_fr(comm%idnode+1))
          idi=ltg(ii)

! Get list limit

          limit=list(-2,i)-list(-1,i)
          If (limit > 0) Then

! calculate interatomic distances

             Do k=1,limit
                j=list(list(-1,i)+k,i)

                xxt(k)=xxx(i)-xxx(j)
                yyt(k)=yyy(i)-yyy(j)
                zzt(k)=zzz(i)-zzz(j)
             End Do

! periodic boundary conditions not needed by LC construction
!
!           Call images(imcon,cell,limit,xxt,yyt,zzt)

! square of distances

             Do k=1,limit
                rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
             End Do

             Do k=1,limit
                j=list(list(-1,i)+k,i)

                rrr=rrt(k)
                If (Abs(chge(j)) > zero_plus .and. rrr < rcut) Then
                   chgprd=-chge(i)*chge(j)/epsq*r4pie0 ! the MINUS signifies the exclusion!!!
                   rsq=rrr**2

! calculate forces

                   coul = chgprd/rrr
                   fcoul = coul/rsq

                   fx = fcoul*xrr
                   fy = fcoul*yrr
                   fz = fcoul*zrr


                   fxx(i)=fxx(i)-fx
                   fyy(i)=fyy(i)-fy
                   fzz(i)=fzz(i)-fz

! redundant calculations copying

                   If (lf_cp) Then
                      ffx(i)=ffx(i)-fx
                      ffy(i)=ffy(i)-fy
                      ffz(i)=ffz(i)-fz
                   End If

! infrequent calculations copying

                   If (l_cp) Then
                      fcx(i)=fcx(i)-fx
                      fcy(i)=fcy(i)-fy
                      fcz(i)=fcz(i)-fz
                   End If

                   If (j <= natms) Then

                      fxx(j)=fxx(j)+fx
                      fyy(j)=fyy(j)+fy
                      fzz(j)=fzz(j)+fz

! redundant calculations copying

                      If (lf_cp) Then
                         ffx(j)=ffx(j)+fx
                         ffy(j)=ffy(j)+fy
                         ffz(j)=ffz(j)+fz
                      End If

! infrequent calculations copying

                      If (l_cp) Then
                         fcx(j)=fcx(j)+fx
                         fcy(j)=fcy(j)+fy
                         fcz(j)=fcz(j)+fz
                      End If

                   End If

                   If (j <= natms .or. idi < ltg(j)) Then

! calculate potential energy

                     engcpe_fr = engcpe_fr + coul

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
          Write(nrite,'(/,1x,a,i0)') 'poisson_frzn_forces deallocation failure 2, node: ', comm%idnode
          Call error(0)
       End If

    End If

! virial

    vircpe_fr = -engcpe_fr

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

    If (lf_cp) Then
       ef_fr=engcpe_fr
       vf_fr=vircpe_fr

       sf_fr(1) = strs1
       sf_fr(2) = strs2
       sf_fr(3) = strs3
       sf_fr(4) = strs2
       sf_fr(5) = strs5
       sf_fr(6) = strs6
       sf_fr(7) = strs3
       sf_fr(8) = strs6
       sf_fr(9) = strs9
    End If

! infrequent calculations copying

    If (l_cp) Then
       e_fr=engcpe_fr
       v_fr=vircpe_fr

       s_fr(1) = strs1
       s_fr(2) = strs2
       s_fr(3) = strs3
       s_fr(4) = strs2
       s_fr(5) = strs5
       s_fr(6) = strs6
       s_fr(7) = strs3
       s_fr(8) = strs6
       s_fr(9) = strs9
    End If

    Deallocate (l_ind,nz_fr, Stat=fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'poisson_frzn_forces deallocation failure, node: ', comm%idnode
       Call error(0)
    End If

  End Subroutine poisson_frzn_forces

End Module poisson_module
