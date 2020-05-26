Module poisson
  Use kinds,           Only : wp,wi
  Use comms,           Only : gsum,comms_type,ExchgGrid_tag,gsend, &
    gwait,girecv,gtime
  Use domains,         Only : domains_type
  Use constants,           Only : fourpi,r4pie0,half_minus,zero_plus
  Use configuration,   Only : configuration_type,&
                             IMCON_NOPBC,&
                             IMCON_CUBIC,&
                             IMCON_ORTHORHOMBIC,&
                             IMCON_PARALLELOPIPED,&
                             IMCON_SLAB,&
                             IMCON_TRUNC_OCTO,&
                             IMCON_RHOMBIC_DODEC,&
                             IMCON_HEXAGONAL
  Use ewald,           Only : ewald_type
  Use errors_warnings, Only : error,info
  Use numerics,        Only : dcell,invert
  Use parallel_fft,    Only : adjust_kmax
  Use neighbours,      Only : neighbours_type
  Use electrostatic,   Only : electrostatic_type
  Implicit None

  Private

  !> Type containing poisson data
  Type, Public :: poisson_type
    Private

    Logical :: debug =.true.
    Logical :: initialized=.false.
    Logical :: converged

    Integer( Kind = wi ) :: lnxl,lnyl,lnzl
    Integer( Kind = wi ) :: lnxu,lnyu,lnzu
    Integer( Kind = wi ) :: maxbicgst

    Integer( Kind = wi ) :: ixb,iyb,izb
    Integer( Kind = wi ) :: ixt,iyt,izt
    Integer( Kind = wi ) :: block_x,block_y,block_z

    !> Local copy of map to be modified
    Integer( Kind = wi ) :: lmap(1:26)

    Integer( Kind = wi ), Public :: mxitjb
    !> Conjugate gradients maximum iterations
    Integer( Kind = wi ), Public :: mxitcg
    Integer( Kind = wi ) :: xhalo
    !Integer( Kind = wi ) :: fail(1:10)=0

    !> Solver tolerance
    Real( Kind = wp ), Public :: eps
    !> Solver spacing
    Real( Kind = wp ) :: delta
    Real( Kind = wp ) :: celprp(1:10)
    Real( Kind = wp ) :: normb
    Real( Kind = wp ) :: normphi0,normphi1
    Real( Kind = wp ) :: ccxxx,ccyyy,cczzz
    Real( Kind = wp ) :: ccsum
    Real( Kind = wp ) :: kmaxa_r,kmaxb_r,kmaxc_r

    !> electrostatic potential
    Real( Kind = wp ), Dimension(:,:,:), Allocatable :: phi
    !> charge density
    Real( Kind = wp ), Dimension(:,:,:), Allocatable :: b
    !> biCGStab variables
    Real( Kind = wp ), Dimension(:,:,:), Allocatable :: F, s, r, p, r0, v, t
    Real( Kind = wp ), Dimension(:,:,:), Allocatable :: pattern, phi0
  Contains
    Private

    Final :: cleanup
  End Type poisson_type

  Public :: poisson_forces,poisson_excl_forces,poisson_frzn_forces

Contains

  Subroutine poisson_forces(engcpe,vircpe,stress,pois,electro,domain,config,ewld,comm)

    Real( Kind = wp ), Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Intent( InOut ) :: stress(1:9)
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ), Intent( InOut )     :: config
    Type( ewald_type ), Intent( InOut ) :: ewld
    Type(comms_type), Intent( InOut )     :: comm

    Real( Kind = wp ) :: eng,virr
    Real( Kind = wp ) :: strs(1:9)

    If (.not.pois%initialized) Then

      pois%converged=.false.
      pois%maxbicgst =0 ! the number of steps to solve eq. without a guess (pois%phi=0)
      pois%delta=1.0_wp/electro%alpha
      Call biCGStab_init(pois,config,domain,ewld)

      Call biCGStab_charge_density(pois,electro,config,comm)

      If (pois%normb > zero_plus) Then
        Call biCGStab_solver_omp(pois,comm)
      Else
        Call error(0,'pois%normb too small')
      End If

      pois%initialized=.true.

    Else

      Call biCGStab_charge_density(pois,electro,config,comm)

    End If

    If (pois%normb > zero_plus) Then
      Call P_solver_omp(pois,comm)
      Call biCGStab_calc_forces(eng,virr,strs,pois,config,ewld,comm)
    Else
      Call error(0,'pois%normb too small')
    End If

    engcpe=engcpe+eng
    vircpe=vircpe+virr
    stress=stress+strs

  End Subroutine poisson_forces

  Subroutine biCGStab_init(pois,config,domain,ewld)
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ),     Intent( InOut ) :: config
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( ewald_type ), Intent( InOut ) :: ewld

    ! calculates preambles
    ! copy DD mapping
    Integer :: fail(10)
    Character( Len = 256 ) :: message
    pois%lmap=domain%map

    If (config%imcon == IMCON_NOPBC .or. config%imcon == IMCON_SLAB) Then
      If      (config%imcon == IMCON_NOPBC) Then
        If (domain%idx == 0     ) pois%lmap(1)=-1
        If (domain%idy == 0     ) pois%lmap(3)=-1
        If (domain%idz == 0     ) pois%lmap(5)=-1
        If (domain%idx == domain%nx-1) pois%lmap(2)=-1
        If (domain%idy == domain%ny-1) pois%lmap(4)=-1
        If (domain%idz == domain%nz-1) pois%lmap(6)=-1
      Else If (config%imcon == IMCON_SLAB) Then
        If (domain%idz == 0     ) pois%lmap(5)=-1
        If (domain%idz == domain%nz-1) pois%lmap(6)=-1
      End If
    End If

    ! (re)set grid size along x, y and z direction

    Call dcell(config%cell,pois%celprp)

    ewld%fft_dim_a=Nint(pois%celprp(7)/pois%delta)
    ewld%fft_dim_b=Nint(pois%celprp(8)/pois%delta)
    ewld%fft_dim_c=Nint(pois%celprp(9)/pois%delta)

    ! adjust accordingly to processor grid restrictions

    Call adjust_kmax(ewld%fft_dim_a, domain%nx)
    Call adjust_kmax(ewld%fft_dim_b, domain%ny)
    Call adjust_kmax(ewld%fft_dim_c, domain%nz)

    ! 3D charge array construction (bottom and top) indices and block size

    pois%ixb=domain%idx*(ewld%fft_dim_a/domain%nx)+1
    pois%ixt=(domain%idx+1)*(ewld%fft_dim_a/domain%nx)
    pois%iyb=domain%idy*(ewld%fft_dim_b/domain%ny)+1
    pois%iyt=(domain%idy+1)*(ewld%fft_dim_b/domain%ny)
    pois%izb=domain%idz*(ewld%fft_dim_c/domain%nz)+1
    pois%izt=(domain%idz+1)*(ewld%fft_dim_c/domain%nz)

    pois%kmaxa_r=Real(ewld%fft_dim_a, wp)
    pois%kmaxb_r=Real(ewld%fft_dim_b, wp)
    pois%kmaxc_r=Real(ewld%fft_dim_c, wp)

    pois%block_x=ewld%fft_dim_a/domain%nx
    pois%block_y=ewld%fft_dim_b/domain%ny
    pois%block_z=ewld%fft_dim_c/domain%nz

    ! new 1 grid link halo inclusive object per domain distributed sizes

    pois%lnxl=1-1
    pois%lnyl=1-1
    pois%lnzl=1-1
    pois%lnxu=pois%block_x+1
    pois%lnyu=pois%block_y+1
    pois%lnzu=pois%block_z+1

    ! extra halo on top of the one grid cell halo, according to size of differentiation stencil

    pois%xhalo=2+ewld%bspline1-ewld%bspline

    ! allocate vectors for biCGStab

    Allocate ( pois%t(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) ,  Stat = fail( 1) )
    Allocate ( pois%b(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) ,  Stat = fail( 2) )
    Allocate ( pois%F(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) ,  Stat = fail( 3) )
    Allocate ( pois%s(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) ,  Stat = fail( 4) )
    Allocate ( pois%r(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) ,  Stat = fail( 5) )
    Allocate ( pois%p(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) ,  Stat = fail( 6) )
    Allocate (  pois%phi(pois%lnxl-pois%xhalo:pois%lnxu+pois%xhalo, &
      pois%lnyl-pois%xhalo:pois%lnyu+pois%xhalo, &
      pois%lnzl-pois%xhalo:pois%lnzu+pois%xhalo) ,       Stat = fail( 7) )
    Allocate ( pois%phi0(pois%lnxl-pois%xhalo:pois%lnxu+pois%xhalo, &
      pois%lnyl-pois%xhalo:pois%lnyu+pois%xhalo, &
      pois%lnzl-pois%xhalo:pois%lnzu+pois%xhalo) ,       Stat = fail( 8) )
    Allocate ( pois%r0(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) , Stat = fail( 9) )
    Allocate (  pois%v(pois%lnxl:pois%lnxu,pois%lnyl:pois%lnyu,pois%lnzl:pois%lnzu) , Stat = fail(10) )
    If (Any(fail(1:10) > 0)) Then
      Write(message,'(a)') 'biCGStab_init allocation failure'
      Call error(0,message)
    End If

    ! Initialise

    pois%F   = 0.0_wp
    pois%s   = 0.0_wp
    pois%r   = 0.0_wp
    pois%p   = 0.0_wp
    pois%phi = 0.0_wp
    pois%r0  = 0.0_wp
    pois%v   = 0.0_wp
    pois%t   = 0.0_wp

  End Subroutine biCGStab_init

  Subroutine biCGStab_charge_density(pois,electro,config,comm)

    ! calculates charge dansity at 0th order

    Type( poisson_type ), Intent( InOut ) :: pois
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( configuration_type ), Intent( InOut ) :: config
    Type(comms_type), Intent( InOut ) :: comm

    Integer           :: i,j,k,n
    Real( Kind = wp ) :: reps0dv, txx,tyy,tzz, det,rcell(9)

    reps0dv=fourpi*r4pie0*electro%alpha/electro%eps ! dv collapsed to dr=pois%delta=1/electro%alpha

    ! get reciprocal config%cell
    Call invert(config%cell,rcell,det)
    If (Abs(det) < 1.0e-6_wp) Call error(120)

#ifdef __OPENMP
    Block
      Real( Kind = wp ), Dimension(:,:,:), Allocatable :: b
      Real( Kind = wp ) :: ccsum
      Real( Kind = wp ) :: ccxxx,ccyyy,cczzz

      Allocate(b(Lbound(pois%b,1):Ubound(pois%b,1), &
        Lbound(pois%b,2):Ubound(pois%b,2), &
        Lbound(pois%b,3):Ubound(pois%b,3)))
      b=0.0_wp

      ccsum=0.0_wp
      ccxxx=0.0_wp
      ccyyy=0.0_wp
      cczzz=0.0_wp

      !$omp paralleldo default(shared) private(n) reduction(+:b,ccsum,ccxxx,ccyyy,cczzz)
      Do n=1,config%natms !config%nlast
        txx=pois%kmaxa_r*(rcell(1)*config%parts(n)%xxx+rcell(4)*config%parts(n)%yyy+&
          rcell(7)*config%parts(n)%zzz+0.5_wp)
        tyy=pois%kmaxb_r*(rcell(2)*config%parts(n)%xxx+rcell(5)*config%parts(n)%yyy+&
          rcell(8)*config%parts(n)%zzz+0.5_wp)
        tzz=pois%kmaxc_r*(rcell(3)*config%parts(n)%xxx+rcell(6)*config%parts(n)%yyy+&
          rcell(9)*config%parts(n)%zzz+0.5_wp)

        ! global indices

        i=Int(txx)
        j=Int(tyy)
        k=Int(tzz)

        ! get local indces

        i = i - pois%ixb + 2
        j = j - pois%iyb + 2
        k = k - pois%izb + 2

        If ( (i < 1 .or. i > pois%block_x) .or. &
          (j < 1 .or. j > pois%block_y) .or. &
          (k < 1 .or. k > pois%block_z) .or. &
          (Abs(config%parts(n)%chge) <= zero_plus) ) Cycle

        b(i,j,k)=b(i,j,k)+config%parts(n)%chge*reps0dv

        ccsum=ccsum+Abs(config%parts(n)%chge)
        ccxxx=ccxxx+Abs(config%parts(n)%chge)*config%parts(n)%xxx
        ccyyy=ccyyy+Abs(config%parts(n)%chge)*config%parts(n)%yyy
        cczzz=cczzz+Abs(config%parts(n)%chge)*config%parts(n)%zzz
      End Do
      !$omp End paralleldo

      pois%b=b

      pois%ccsum=ccsum
      pois%ccxxx=ccxxx
      pois%ccyyy=ccyyy
      pois%cczzz=cczzz
    End Block
#else
    Do n=1,config%natms !config%nlast
      txx=pois%kmaxa_r*(rcell(1)*config%parts(n)%xxx+rcell(4)*config%parts(n)%yyy+&
        rcell(7)*config%parts(n)%zzz+0.5_wp)
      tyy=pois%kmaxb_r*(rcell(2)*config%parts(n)%xxx+rcell(5)*config%parts(n)%yyy+&
        rcell(8)*config%parts(n)%zzz+0.5_wp)
      tzz=pois%kmaxc_r*(rcell(3)*config%parts(n)%xxx+rcell(6)*config%parts(n)%yyy+&
        rcell(9)*config%parts(n)%zzz+0.5_wp)

      ! global indices

      i=Int(txx)
      j=Int(tyy)
      k=Int(tzz)

      ! get local indces

      i = i - pois%ixb + 2
      j = j - pois%iyb + 2
      k = k - pois%izb + 2

      If ( (i < 1 .or. i > pois%block_x) .or. &
        (j < 1 .or. j > pois%block_y) .or. &
        (k < 1 .or. k > pois%block_z) .or. &
        (Abs(config%parts(n)%chge) <= zero_plus) ) Cycle

      pois%b(i,j,k)=pois%b(i,j,k)+config%parts(n)%chge*reps0dv

      pois%ccsum=pois%ccsum+Abs(config%parts(n)%chge)
      pois%ccxxx=pois%ccxxx+Abs(config%parts(n)%chge)*config%parts(n)%xxx
      pois%ccyyy=pois%ccyyy+Abs(config%parts(n)%chge)*config%parts(n)%yyy
      pois%cczzz=pois%cczzz+Abs(config%parts(n)%chge)*config%parts(n)%zzz
    End Do

#endif

    pois%normb=Sum(pois%b**2)
    Call gsum(comm,pois%normb)

    Call gsum(comm,pois%ccsum)
    Call gsum(comm,pois%ccxxx)
    Call gsum(comm,pois%ccyyy)
    Call gsum(comm,pois%cczzz)

  End Subroutine biCGStab_charge_density

  Recursive Subroutine P_solver_omp(pois,comm)

    Type( poisson_type ), Intent( InOut ) :: pois
    Type(comms_type), Intent( InOut )   :: comm

    Real( Kind = wp ) :: normphi0_local, normphi1_local
    Real( Kind = wp ) :: Totstart, Totend, dphi
    Integer           :: mmm, i,j,k
    Real( Kind = wp ) :: element,  sm1,sm2,sm3,sm4 ! SM stands for Stoyan Markov (long live!!!)
    Character ( Len = 80 ) :: message

    Call gtime(Totstart)
    pois%converged=.false.
    normphi0_local=0.0_wp

    !$omp paralleldo default(shared) private(k) reduction(+:normphi0_local)
    Do k=1,pois%block_z
      normphi0_local = normphi0_local + Sum(pois%phi(1:pois%block_x,1:pois%block_y,k)**2)
    End Do
    !$omp End paralleldo

    pois%normphi0=normphi0_local
    Call gsum(comm,pois%normphi0)

    Do mmm=1,pois%mxitjb+3 ! Za seki sluchaj, proverka sa stabilno shozhdane

      Call biCGStab_exchange_halo(pois%phi,pois%xhalo,pois,comm)

      sm1 = -600.0_wp/144.0_wp
      sm2 =   60.0_wp/144.0_wp
      sm3 =   18.0_wp/144.0_wp
      sm4 =    3.0_wp/144.0_wp

      !$omp paralleldo default(shared) private(element,k,j,i)
      Do k=pois%lnzl+1,pois%lnzu-1
        Do j=pois%lnyl+1,pois%lnyu-1
          Do i=pois%lnxl+1,pois%lnxu-1
            element=0.0_wp
            element=element+sm4*pois%phi(i - 1 ,j - 1 ,k - 1 )
            element=element+sm3*pois%phi(i     ,j - 1 ,k - 1 )
            element=element+sm4*pois%phi(i + 1 ,j - 1 ,k - 1 )
            element=element+sm3*pois%phi(i - 1 ,j     ,k - 1 )
            element=element+sm2*pois%phi(i     ,j     ,k - 1 )
            element=element+sm3*pois%phi(i + 1 ,j     ,k - 1 )
            element=element+sm4*pois%phi(i - 1 ,j + 1 ,k - 1 )
            element=element+sm3*pois%phi(i     ,j + 1 ,k - 1 )
            element=element+sm4*pois%phi(i + 1 ,j + 1 ,k - 1 )
            element=element+sm3*pois%phi(i - 1 ,j - 1 ,k     )
            element=element+sm2*pois%phi(i     ,j - 1 ,k     )
            element=element+sm3*pois%phi(i + 1 ,j - 1 ,k     )
            element=element+sm2*pois%phi(i - 1 ,j     ,k     )
            element=element+sm2*pois%phi(i + 1 ,j     ,k     )
            element=element+sm3*pois%phi(i - 1 ,j + 1 ,k     )
            element=element+sm2*pois%phi(i     ,j + 1 ,k     )
            element=element+sm3*pois%phi(i + 1 ,j + 1 ,k     )
            element=element+sm4*pois%phi(i - 1 ,j - 1 ,k + 1 )
            element=element+sm3*pois%phi(i     ,j - 1 ,k + 1 )
            element=element+sm4*pois%phi(i + 1 ,j - 1 ,k + 1 )
            element=element+sm3*pois%phi(i - 1 ,j     ,k + 1 )
            element=element+sm2*pois%phi(i     ,j     ,k + 1 )
            element=element+sm3*pois%phi(i + 1 ,j     ,k + 1 )
            element=element+sm4*pois%phi(i - 1 ,j + 1 ,k + 1 )
            element=element+sm3*pois%phi(i     ,j + 1 ,k + 1 )
            element=element+sm4*pois%phi(i + 1 ,j + 1 ,k + 1 )
            pois%phi0(i,j,k)=(-pois%b(i,j,k)+element)/(-sm1)
          End Do
        End Do
      End Do
      !$omp End paralleldo

      normphi1_local=0.0_wp
      !$omp paralleldo default(shared) private(k) reduction(+:normphi1_local)
      Do k=1,pois%block_z
        pois%phi(:,:,k)=pois%phi0(:,:,k)
        normphi1_local=normphi1_local+Sum(pois%phi(1:pois%block_x,1:pois%block_y,k)**2)
      End Do
      !$omp End paralleldo

      pois%normphi1=normphi1_local
      Call gsum(comm, pois%normphi1)

      dphi=Abs(pois%normphi1-pois%normphi0)
      If (dphi <= pois%eps*pois%normb) Then
        Call biCGStab_exchange_halo(pois%phi,pois%xhalo,pois,comm)
        Call gtime(Totend)
        Write (message,'(a,i0,a,g0,g0)')  "Jacobi it = ", mmm ,&
          &"d|Coulomb potential|/NormB >>>",&
          &Abs(pois%normphi1 - pois%normphi0)/pois%normb, Totend-Totstart
        Call info(message,.true.)
        pois%converged=.true.
        Return
      End If

      If (mmm > pois%mxitjb) Call biCGStab_solver_omp(pois,comm)
      pois%normphi0=pois%normphi1

    End Do

    Write(message,'(a,i0,a)') "poisson solver not converged in ",pois%mxitcg,"steps..."
    Call info(message)
    pois%converged=.false.

  End Subroutine P_solver_omp

  Recursive Subroutine biCGStab_solver_omp(pois,comm)

    Type( poisson_type ), Intent( InOut ) :: pois
    Type(comms_type), Intent( InOut )  :: comm

    Real( Kind = wp ) :: normphi0_local,normphi1_local
    Real( Kind = wp ) :: alfa, beta, omega, rho1,rho0,rv,tt,ts
    Real( Kind = wp ) :: Totstart,Totend

    Integer :: mmm, kk
    Character( Len = 80 ) :: message

    Call gtime(Totstart)

    alfa=1.0_wp
    omega=1.0_wp
    rho1=1.0_wp
    rho0=1.0_wp

    Call biCGStab_exchange_halo(pois%phi,pois%xhalo,pois,comm)
    Call Adot_omp(pois%phi,pois%xhalo,pois%F,0,pois)

    !$omp parallel default(shared) private(kk)
    !$omp paralleldo
    Do kk=1,pois%block_z
      pois%r0(1:pois%block_x,1:pois%block_y,kk) = &
        pois%b(1:pois%block_x,1:pois%block_y,kk) - pois%F(1:pois%block_x,1:pois%block_y,kk)

      pois%r(1:pois%block_x,1:pois%block_y,kk)  = pois%r0(1:pois%block_x,1:pois%block_y,kk)
    End Do
    !$omp End paralleldo

    !$omp paralleldo
    Do kk=pois%lnzl,pois%lnzu
      pois%p(:,:,kk) = 0.0_wp
      pois%v(:,:,kk) = 0.0_wp
    End Do
    !$omp End paralleldo

    normphi0_local=0.0_wp
    !$omp paralleldo reduction(+:normphi0_local)
    Do kk=1,pois%block_z
      normphi0_local = normphi0_local + Sum(pois%phi(1:pois%block_x,1:pois%block_y,kk)**2)
    End Do
    !$omp End parallel

    pois%normphi0=normphi0_local
    Call gsum(comm,pois%normphi0)
    Do mmm=1,pois%mxitcg
      rho0 = rho1
      rho1 = 0.0_wp
      !$omp paralleldo reduction(+:rho1)
      Do kk=1,pois%block_z
        rho1 = rho1 + &
          Sum(pois%r0(1:pois%block_x,1:pois%block_y,kk)*pois%r(1:pois%block_x,1:pois%block_y,kk))
      End Do
      !$omp End paralleldo
      Call gsum(comm,rho1)

      beta = (rho1/rho0) * (alfa/omega)

      !$omp paralleldo default(shared) private(kk)
      Do kk=1,pois%block_z
        pois%p(1:pois%block_x,1:pois%block_y,kk) = pois%r(1:pois%block_x,1:pois%block_y,kk) + &
          beta * (pois%p(1:pois%block_x,1:pois%block_y,kk) - omega * pois%v(1:pois%block_x,1:pois%block_y,kk))
      End Do
      !$omp End paralleldo

      Call biCGStab_exchange_halo(pois%p,0,pois,comm)
      Call Adot_omp(pois%p,0,pois%v,0,pois)

      rv=0.0_wp
      !$omp parallel Do default(shared) reduction(+:rv) private(kk)
      Do kk=1,pois%block_z
        rv = rv + Sum(pois%r0(1:pois%block_x,1:pois%block_y,kk) * pois%v(1:pois%block_x,1:pois%block_y,kk))
      End Do
      !$omp End paralleldo
      Call gsum(comm,rv)

      alfa=rho1/rv

      !$omp paralleldo default(shared) private(kk)
      Do kk=1,pois%block_z
        pois%s(1:pois%block_x,1:pois%block_y,kk) = pois%r(1:pois%block_x,1:pois%block_y,kk) - &
          alfa * pois%v(1:pois%block_x,1:pois%block_y,kk)
      End Do
      !$omp End paralleldo

      Call biCGStab_exchange_halo(pois%s,0,pois,comm)
      Call Adot_omp(pois%s,0,pois%t,0,pois)

      tt=0.0_wp
      !$omp paralleldo default(shared) private(kk) reduction(+:tt)
      Do kk=1,pois%block_z
        tt = tt + Sum(pois%t(1:pois%block_x,1:pois%block_y,kk)**2)
      End Do
      !$omp End paralleldo
      Call gsum(comm,tt)

      ts=0.0_wp
      !$omp paralleldo default(shared) private(kk) reduction(+:ts)
      Do kk=1,pois%block_z
        ts = ts + Sum(pois%t(1:pois%block_x,1:pois%block_y,kk)*pois%s(1:pois%block_x,1:pois%block_y,kk))
      End Do
      !$omp End paralleldo
      Call gsum(comm,ts)

      omega=ts/tt

      !$omp paralleldo default(shared) private(kk)
      Do kk=1,pois%block_z
        pois%phi(1:pois%block_x,1:pois%block_y,kk) = pois%phi(1:pois%block_x,1:pois%block_y,kk) + &
          alfa * pois%p(1:pois%block_x,1:pois%block_y,kk) + omega * pois%s(1:pois%block_x,1:pois%block_y,kk)
      End Do
      !$omp End paralleldo

      normphi1_local=0.0_wp
      !$omp paralleldo default(shared) private(kk) reduction(+:normphi1_local)
      Do kk=1,pois%block_z
        normphi1_local = normphi1_local + Sum(pois%phi(1:pois%block_x,1:pois%block_y,kk)**2)
      End Do
      !$omp End paralleldo
      pois%normphi1=normphi1_local
      Call gsum(comm,pois%normphi1)

      If (pois%maxbicgst > 0 .and. mmm > pois%maxbicgst) Then
        !$omp paralleldo default(shared) private(kk)
        Do kk=pois%lnzl,pois%lnzu
          pois%phi(:,:,kk) = 0.0_wp
        End Do
        !$omp End paralleldo

        pois%maxbicgst=0
        pois%converged=.false.
        Call info("bicgstab exceeded *pois%maxbicgst*... reset potential... pois%converged=.false.",.true.)
        Return
      End If
      Call gtime(Totend)
      Write (message,'(a,i0,a,g0,g0)')  "biCGStab it = ", mmm ,&
        &"d|Coulomb potential|/NormB >>>",(Abs(pois%normphi1 - pois%normphi0)/pois%normb),&
        &Totend-Totstart
      Call info(message)

      If (Abs(pois%normphi1 - pois%normphi0) <= pois%eps*pois%normb) Then
        Call biCGStab_exchange_halo(pois%phi,pois%xhalo,pois,comm)

        Call gtime(Totend)
        Write (message,'(a,i0,a,g0,g0)')  "biCGStab it = ", mmm ,&
          &"d|Coulomb potential|/NormB >>>",(Abs(pois%normphi1 - pois%normphi0)/pois%normb),&
          &Totend-Totstart
        Call info(message,.true.)
        If (pois%maxbicgst == 0) Then
          pois%maxbicgst=mmm
        End If
        pois%converged=.true.
        Write(message,'(a,i0,a)') "*pois%maxbicgst* set to ",pois%maxbicgst," pois%converged=.true."
        Call info(message,.true.)
        Return
      End If

      !$omp paralleldo default(shared) private(kk)
      Do kk=1,pois%block_z
        pois%r(1:pois%block_x,1:pois%block_y,kk) = pois%s(1:pois%block_x,1:pois%block_y,kk) - &
          omega * pois%t(1:pois%block_x,1:pois%block_y,kk)
      End Do
      !$omp End paralleldo

      pois%normphi0=pois%normphi1
    End Do
    pois%converged=.false.
    Call info("maxit reached... pois%converged=.false.",.true.)

  End Subroutine biCGStab_solver_omp

  Subroutine biCGStab_exchange_halo(vec,xtra,pois,comm)

    Integer,           Intent( In    ) :: xtra
    Type( poisson_type ), Intent( InOut ) :: pois
    Real( Kind = wp ), Intent( InOut ) :: vec( pois%lnxl-xtra:pois%lnxu+xtra, &
      pois%lnyl-xtra:pois%lnyu+xtra, &
      pois%lnzl-xtra:pois%lnzu+xtra )
    Type( comms_type ), Intent( InOut ) :: comm

    Integer :: me, lx,ly,lz

    ! What's my name?

    me = comm%idnode

    ! Find length of sides of the domain

    lx = pois%ixt - pois%ixb + 1
    ly = pois%iyt - pois%iyb + 1
    lz = pois%izt - pois%izb + 1

    ! +X direction face - negative halo

    Call exchange_grid_halo( pois%lmap(1),                     pois%lmap(2), &
      lx-(xtra+1)+1, lx  ,         1,              ly,            1,              lz, &
      1-(xtra+1),    1-1,          1,              ly,            1,              lz,comm)

    ! -X direction face - positive halo

    Call exchange_grid_halo( pois%lmap(2),                     pois%lmap(1), &
      1,             1+(xtra+1)-1, 1,              ly,            1,              lz, &
      lx+1,          lx+(xtra+1),  1,              ly,            1,              lz,comm)

    ! +Y direction face (including the +&-X faces extensions) - negative halo

    Call exchange_grid_halo( pois%lmap(3),                     pois%lmap(4), &
      1-(xtra+1),    lx+(xtra+1),  ly-(xtra+1)+1, ly,             1,              lz, &
      1-(xtra+1),    lx+(xtra+1),  1-(xtra+1)  ,  1-1,            1,              lz,comm)

    ! -Y direction face (including the +&-X faces extensions) - positive halo

    Call exchange_grid_halo( pois%lmap(4),                     pois%lmap(3), &
      1-(xtra+1),    lx+(xtra+1),  1,              1+(xtra+1)-1,  1, lz, &
      1-(xtra+1),    lx+(xtra+1),  ly+1,           ly+(xtra+1),   1, lz,comm)

    ! +Z direction face (including the +&-Y+&-X faces extensions) - negative halo

    Call exchange_grid_halo( pois%lmap(5),                     pois%lmap(6), &
      1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  lz-(xtra+1)+1, lz, &
      1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  1-(xtra+1)  ,  1-1,comm)

    ! -Z direction face (including the +&-Y+&-X faces extensions) - positive halo

    Call exchange_grid_halo( pois%lmap(6),                     pois%lmap(5), &
      1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  1,              1+(xtra+1)-1, &
      1-(xtra+1),    lx+(xtra+1),  1-(xtra+1),    ly+(xtra+1),  lz+1,           lz+(xtra+1),comm)

  Contains

    Subroutine exchange_grid_halo(     from,       to,           &
        xlb, xlt, ylb, ylt, zlb, zlt, &
        xdb, xdt, ydb, ydt, zdb, zdt,comm )

      Integer, Intent( In    ) :: from, to
      Integer, Intent( In    ) :: xlb, ylb, zlb
      Integer, Intent( In    ) :: xlt, ylt, zlt
      Integer, Intent( In    ) :: xdb, ydb, zdb
      Integer, Intent( In    ) :: xdt, ydt, zdt
      Type( comms_type), Intent( InOut ) :: comm

      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: send_buffer
      Real( Kind = wp ), Dimension( :, :, : ), Allocatable :: recv_buffer

      Integer :: fail
      Integer :: length
      Character( Len = 256 ) :: message
      ! If the processor to receive FROM is actually ME it means there is
      ! only one processor along this axis (so the processor to send TO is
      ! also ME) and so no message passing need be done.  However, there
      ! is a catch.  The domain decomposed spme_forces routine expects data
      ! that would be `seen' through the periodic boundary conditions to be
      ! actually copied from the high positive indices to negative ones and
      ! vice-versa.  The Else clause catches this.

      If ( from /= me ) Then

        ! Allocate send and receive buffers (of the same size!!!)
        ! so all can be sent and received as one message!!!

        If (from > -1) Then
          Allocate (recv_buffer(xdb:xdt,ydb:ydt,zdb:zdt), Stat = fail)
          If (fail > 0) Then
            Write(message,'(a)') 'exchange_grid_halo receive allocation failure'
            Call error(0,message)
          End If

          Call girecv(comm,recv_buffer(:,:,:),from,ExchgGrid_tag)
        End If

        If (to   > -1) Then
          Allocate (send_buffer(xlb:xlt,ylb:ylt,zlb:zlt), Stat = fail)
          If (fail > 0) Then
            Write(message,'(a)') 'exchange_grid_halo send allocation failure'
            Call error(0,message)
          End If

          ! Copy the data to be sent

          send_buffer = vec( xlb:xlt, ylb:ylt, zlb:zlt )

          Call gsend(comm,send_buffer(:,:,:),to,ExchgGrid_tag)
        End If

        ! Exchange the data

        If (from > -1) Then
          Call gwait(comm)

          ! Copy the received data into the domain halo

          vec(xdb:xdt, ydb:ydt, zdb:zdt) = recv_buffer

          ! And, as my mum told me, leave things as we found them

          Deallocate(recv_buffer, Stat = fail)
          If (fail > 0) Then
            Write(message,'(a,i0)') 'exchange_grid_halo receive deallocation failure'
            Call error(0,message)
          End If
        End If

        If (to   > -1) Then
          Deallocate(send_buffer, Stat = fail)
          If (fail > 0) Then
            Write(message,'(a,i0)') 'exchange_grid_halo send deallocation failure'
            Call error(0,message)
          End If
        End If

      Else

        ! Simple on node copy - as sizes are the same as 3D shapes

        vec(xdb:xdt, ydb:ydt, zdb:zdt) = vec(xlb:xlt, ylb:ylt, zlb:zlt)

      End If

    End Subroutine exchange_grid_halo

  End Subroutine biCGStab_exchange_halo

  Subroutine Adot_omp(vec,vx,res,rx,pois)

    Integer          , Intent( In    ) :: vx,rx
    Type( poisson_type ), Intent( InOut ) :: pois
    Real( Kind = wp ), Intent( In    ) :: vec(pois%lnxl-vx:pois%lnxu+vx, &
      pois%lnyl-vx:pois%lnyu+vx, &
      pois%lnzl-vx:pois%lnzu+vx)
    Real( Kind = wp ), Intent(   Out ) :: res(pois%lnxl-rx:pois%lnxu+rx, &
      pois%lnyl-rx:pois%lnyu+rx, &
      pois%lnzl-rx:pois%lnzu+rx)

    Integer          :: i,j,k,ioff,joff,koff
    Real( Kind = wp) :: sm(4), element

    sm(1) = -600.0_wp/144.0_wp
    sm(2) =   60.0_wp/144.0_wp
    sm(3) =   18.0_wp/144.0_wp
    sm(4) =    3.0_wp/144.0_wp

    !$omp paralleldo default(shared) private(element,k,j,i,ioff,joff,koff)
    Do k=pois%lnzl+1,pois%lnzu-1
      Do j=pois%lnyl+1,pois%lnyu-1
        Do i=pois%lnxl+1,pois%lnxu-1
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

  Subroutine biCGStab_calc_forces(cenergy,vir,stress,pois,config,ewld,comm)

    Real( Kind = wp ), Intent( InOut ) :: cenergy, vir, stress(1:9)
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( configuration_type ),  Intent( InOut ) :: config
    Type( ewald_type ), Intent( In    ) :: ewld
    Type( comms_type), Intent( InOut ) :: comm

    Integer       :: i,j,k,n, ii,jj,kk
    Real(Kind=wp) :: reps0dv,txx,tyy,tzz, det,rcell(1:9),        &
      cfxx, cfyy, cfzz,                           &
      uenergy,r8veps0

    If(.Not.pois%converged) Then
      Call info("poisson solver not converged",.true.)
      Return
    End If

    reps0dv=fourpi*r4pie0/pois%delta
    r8veps0=fourpi*r4pie0/(8.0_wp*pois%delta**3)

    ! initialization

    cenergy=0.0_wp
    vir=0.0_wp
    stress=0.0_wp

    Call invert(config%cell,rcell,det)
    If (Abs(det) < 1.0e-6_wp) Call error(120)

    !$omp paralleldo default(shared) private(n,txx,tyy,tzz,ii,jj,kk,i,j,k,cfxx,cfyy,cfzz) &
    !$omp reduction(+:stress,vir,cenergy)
    Do n=1,config%natms
      If (Abs(config%parts(n)%chge) < Abs(half_minus)) Cycle

      txx=Real(ewld%fft_dim_a,wp)*(rcell(1)*config%parts(n)%xxx+rcell(4)*config%parts(n)%yyy+&
        rcell(7)*config%parts(n)%zzz+half_minus)
      tyy=Real(ewld%fft_dim_b,wp)*(rcell(2)*config%parts(n)%xxx+rcell(5)*config%parts(n)%yyy+&
        rcell(8)*config%parts(n)%zzz+half_minus)
      tzz=Real(ewld%fft_dim_c,wp)*(rcell(3)*config%parts(n)%xxx+rcell(6)*config%parts(n)%yyy+&
        rcell(9)*config%parts(n)%zzz+half_minus)

      ! global indices

      ii=Int(txx)
      jj=Int(tyy)
      kk=Int(tzz)

      ! get local indces

      i = ii - pois%ixb + 2
      j = jj - pois%iyb + 2
      k = kk - pois%izb + 2

      cfxx=config%parts(n)%chge*dPhidX(i,j,k,pois)
      cfyy=config%parts(n)%chge*dPhidY(i,j,k,pois)
      cfzz=config%parts(n)%chge*dPhidZ(i,j,k,pois)

      ! calculate atomic energy

      uenergy=config%parts(n)%chge*pois%phi(i,j,k)

      !caclulate atomic contribution to the stress tensor

      stress(1)=stress(1)+config%parts(n)%xxx*cfxx
      stress(2)=stress(2)+config%parts(n)%xxx*cfyy
      stress(3)=stress(3)+config%parts(n)%xxx*cfzz
      stress(4)=stress(4)+config%parts(n)%yyy*cfxx
      stress(5)=stress(5)+config%parts(n)%yyy*cfyy
      stress(6)=stress(6)+config%parts(n)%yyy*cfzz
      stress(7)=stress(7)+config%parts(n)%zzz*cfxx
      stress(8)=stress(8)+config%parts(n)%zzz*cfyy
      stress(9)=stress(9)+config%parts(n)%zzz*cfzz

      config%parts(n)%fxx=config%parts(n)%fxx+cfxx
      config%parts(n)%fyy=config%parts(n)%fyy+cfyy
      config%parts(n)%fzz=config%parts(n)%fzz+cfzz

      !add enegry to the total energy

      cenergy=cenergy+uenergy

      vir=vir-uenergy
    End Do
    !$omp End paralleldo

  End Subroutine biCGStab_calc_forces

  Subroutine Write_potential(pois,domain,ewld)
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( domains_type ), Intent( In    ) :: domain
    Type( ewald_type ), Intent( In    ) :: ewld

    Integer :: ounit, kpkp,kk,kpk,tt
    character(len=128) :: filename, line

    ounit=39
    line="pot"
    Write(filename,'(a,i1,i1,i1,a)') Trim(Adjustl(line)),domain%idx,domain%idy,domain%idz,".dx"
    Open(Unit=ounit, File=Trim(Adjustl(filename)), Status='replace')
    Write(ounit,'(a,1x,I7,1x,I7,1x,I7)') "object 1 class gridpositions counts ", pois%block_x,pois%block_y,pois%block_z
    Write(ounit,'(a,1x,f9.3,1x,f9.3,1x,f9.3)') "origin ", &
      Real(domain%idx*(ewld%fft_dim_a/domain%nx),wp)*pois%delta, &
      Real(domain%idy*(ewld%fft_dim_b/domain%ny),wp)*pois%delta, &
      Real(domain%idz*(ewld%fft_dim_c/domain%nz),wp)*pois%delta
    Write(ounit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "pois%delta ", (/pois%delta, 0.0_wp,0.0_wp/)
    Write(ounit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "pois%delta ", (/0.0_wp,pois%delta, 0.0_wp/)
    Write(ounit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "pois%delta ", (/0.0_wp,0.0_wp, pois%delta/)
    Write(ounit,'(a,I5,1x,I5,1x,I5)') "object 2 class gridconnections counts ",pois%block_x,pois%block_y,pois%block_z
    Write(ounit,'(a,1x,I20,1x,a)')'object 3 class array type Double rank 0 items ', &
      &pois%block_x*pois%block_y*pois%block_z , 'data follows'

    tt=0
    line=''
    Do kpkp=1,pois%block_x
      Do kk=1,pois%block_y
        Do kpk=1,pois%block_z
          Write(line,"(a,3E15.5)") Trim(line), Real(pois%phi(kpkp,kk,kpk),4)
          tt=tt+1
          If ( mod(tt,3) == 0 ) Then
            Write(ounit,'(a)') Trim(Adjustl(line))
            line=''
            tt=0
          End If
        End Do
      End Do
    End Do

    Write(ounit,'(a)') Trim(Adjustl(line))
    Write(ounit,'(a)') 'attribute "dep" string "positions"'
    Write(ounit,'(a)') 'object "PME potential (kT/e, T=300K)" class field'
    Write(ounit,'(a)') 'component "positions" value 1'
    Write(ounit,'(a)') 'component "connections" value 2'
    Write(ounit,'(a)') 'component "data" value 3'

    close(Unit=ounit)

  End Subroutine Write_potential

  Subroutine Write_b(pois,domain,ewld)
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( domains_type ), Intent( In    ) :: domain
    Type( ewald_type ), Intent( In    ) :: ewld

    Integer :: ounit, kpkp,kk,kpk,tt
    character(len=128) :: filename, line
    ounit=39
    line="pois%b"
    Write(filename,'(a,i1,i1,i1,a)') Trim(Adjustl(line)),domain%idx,domain%idy,domain%idz,".dx"
    Open(Unit=ounit, File=Trim(Adjustl(filename)), Status='replace')
    Write(ounit,'(a,1x,I7,1x,I7,1x,I7)') "object 1 class gridpositions counts ", pois%block_x,pois%block_y,pois%block_z
    Write(ounit,'(a,1x,f9.3,1x,f9.3,1x,f9.3)') "origin ", &
      Real(domain%idx*(ewld%fft_dim_a/domain%nx),wp)*pois%delta, &
      Real(domain%idy*(ewld%fft_dim_b/domain%ny),wp)*pois%delta, &
      Real(domain%idz*(ewld%fft_dim_c/domain%nz),wp)*pois%delta
    Write(ounit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "pois%delta ", (/Real(pois%delta), 0.0,0.0/)
    Write(ounit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "pois%delta ", (/ 0.0,Real(pois%delta), 0.0/)
    Write(ounit,'(a,1x,f9.6,1x,f9.6,1x,f9.6)') "pois%delta ", (/0.0,0.0, Real(pois%delta)/)
    Write(ounit,'(a,I5,1x,I5,1x,I5)') "object 2 class gridconnections counts ",pois%block_x,pois%block_y,pois%block_z
    Write(ounit,'(a,1x,I20,1x,a)')'object 3 class array type Double rank 0 items ', &
      &pois%block_x*pois%block_y*pois%block_z , 'data follows'

    tt=0
    line=''
    Do kpkp=1,pois%block_x
      Do kk=1,pois%block_y
        Do kpk=1,pois%block_z
          Write(line,"(a,3E15.5)") Trim(line), Real(pois%b(kpkp,kk,kpk),4)
          tt=tt+1
          If ( mod(tt,3) == 0 ) Then
            Write(ounit,'(a)') Trim(Adjustl(line))
            line=''
            tt=0
          End If
        End Do
      End Do
    End Do

    Write(ounit,'(a)') Trim(Adjustl(line))
    Write(ounit,'(a)') 'attribute "dep" string "positions"'
    Write(ounit,'(a)') 'object "PME potential (kT/e, T=300K)" class field'
    Write(ounit,'(a)') 'component "positions" value 1'
    Write(ounit,'(a)') 'component "connections" value 2'
    Write(ounit,'(a)') 'component "data" value 3'

    close(Unit=ounit)

  End Subroutine Write_b

  Function dPhidX(i,j,k,pois) Result(dphi)

    Integer, Intent( In    ) :: i,j,k
    Type( poisson_type ), Intent( InOut ) :: pois

    Integer           :: kk
    Real( Kind = wp ) :: par(3), par4, dphi

    par(1)=45.0_wp
    par(2)=-9.0_wp
    par(3)=1.0_wp
    par4=60.0_wp

    ! sm  du/dx = [45(ui+1 – ui-1 ) -9(ui+2 – ui-2) + (ui+3 –ui-3)] / 60h
    dphi=0.0_wp
    Do kk=1,3
      dphi = dphi + par(kk)*(pois%phi(i+kk,j,k)-pois%phi(i-kk,j,k))
    End Do
    dphi = dphi/(par4*pois%delta)

  End Function dPhidX

  Function dPhidY(i,j,k,pois) Result(dphi)

    Integer, Intent( In    ) :: i,j,k
    Type( poisson_type ), Intent( InOut ) :: pois

    Integer           :: kk
    Real( Kind = wp ) :: par(3), par4, dphi

    par(1)=45.0_wp
    par(2)=-9.0_wp
    par(3)=1.0_wp
    par4=60.0_wp

    ! sm  du/dx = [45(ui+1 – ui-1 ) -9(ui+2 – ui-2) + (ui+3 –ui-3)] / 60h
    dphi=0.0_wp
    Do kk=1,3
      dphi = dphi + par(kk)*(pois%phi(i,j+kk,k)-pois%phi(i,j-kk,k))
    End Do
    dphi = dphi/(par4*pois%delta)

  End Function dPhidY

  Function dPhidZ(i,j,k,pois) Result(dphi)

    Integer, Intent( In    ) :: i,j,k
    Type( poisson_type ), Intent( InOut ) :: pois

    Integer           :: kk
    Real( Kind = wp ) :: par(3), par4, dphi

    par(1)=45.0_wp
    par(2)=-9.0_wp
    par(3)=1.0_wp
    par4=60.0_wp

    ! sm  du/dx = [45(ui+1 – ui-1 ) -9(ui+2 – ui-2) + (ui+3 –ui-3)] / 60h
    dphi=0.0_wp
    Do kk=1,3
      dphi = dphi + par(kk)*(pois%phi(i,j,k+kk)-pois%phi(i,j,k-kk))
    End Do
    dphi = dphi/(par4*pois%delta)

  End Function dPhidZ

  Function d2PhidX2(i,j,k,pois) Result(dphi)

    Integer, Intent( In    ) :: i,j,k
    Type( poisson_type ), Intent( InOut ) :: pois

    Real( Kind = wp ) :: dphi

    dphi=( pois%phi(i+2,j,k) + pois%phi(i-1,j,k) - pois%phi(i+1,j,k) - pois%phi(i-2,j,k) ) / (3.0_wp*pois%delta**2)

  End Function d2PhidX2

  Function d2PhidY2(i,j,k,pois) Result(dphi)

    Integer, Intent( In    ) :: i,j,k
    Type( poisson_type ), Intent( InOut ) :: pois

    Real( Kind = wp ) :: dphi

    dphi=( pois%phi(i,j+2,k) + pois%phi(i,j-1,k) - pois%phi(i,j+1,k) - pois%phi(i,j-2,k) ) / (3.0_wp*pois%delta**2)

  End Function d2PhidY2

  Function d2PhidZ2(i,j,k,pois) Result(dphi)

    Integer, Intent( In    ) :: i,j,k
    Type( poisson_type ), Intent( InOut ) :: pois

    Real( Kind = wp ) :: dphi

    dphi=( pois%phi(i,j,k+2) + pois%phi(i,j,k-1) - pois%phi(i,j,k+1) - pois%phi(i,j,k-2) ) / (3.0_wp*pois%delta**2)

  End Function d2PhidZ2

  Subroutine poisson_excl_forces(iatm,rcut,eps,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex,stress,neigh,config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the exclusion coulombic energy
    ! and! force terms in a periodic system using 1/pois%r potential as required
    ! for a direct space poisson solver
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,eps
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( configuration_type ),               Intent( InOut ) :: config

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

    idi=config%ltg(iatm)
    chgea = config%parts(iatm)%chge

    ! ignore interaction if the charge is zero

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*r4pie0/eps

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! Get limit

      limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

      ! start of primary loop for forces evaluation

      Do m=1,limit

        ! atomic index and charge

        jatm=neigh%list(m,iatm)
        chgprd=config%parts(jatm)%chge

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

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

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

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

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

  Subroutine poisson_frzn_forces(eps,engcpe_fr,vircpe_fr,stress,ewld,neigh,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating corrections to coulombic energy
    ! and forces in a periodic system arising from frozen pairs
    ! for a direct space poisson solver
    !
    ! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
    !       end (and any COM drift removed) but corrections to the stress
    !       and the virial are important as they feed into the system
    !       pressure response.  Constant volume ensembles (ensemble < 20)
    !       need this calculation just once (NOT! - controlled by ewld%lf_fce in
    !       ewald_check<-two_body_forces
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

    Real( Kind = wp ),                   Intent( In    ) :: eps
    Real( Kind = wp ),                   Intent(   Out ) :: engcpe_fr,vircpe_fr
    Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
    Type(ewald_type), Intent( InOut )                    :: ewld
    Type( neighbours_type ), Intent( In    )             :: neigh
    Type( configuration_type ),          Intent( InOut ) :: config
    Type(comms_type), Intent( InOut )                    :: comm

    Integer           :: fail,i,j,k,ii,jj,idi,nzfr,limit
    Real( Kind = wp ) :: det,rcell(1:9),xrr,yrr,zrr,rrr,rsq, &
      chgprd,coul,fcoul,                  &
      fx,fy,fz,xss,yss,zss,               &
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
      Write(message,'(a)') 'poisson_frzn_forces allocation failure'
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
    Call gsum(comm,nz_fr)
    nz_fr(0) = Sum(nz_fr(0:comm%idnode)) ! Offset

    nzfr = Sum(nz_fr(1:comm%mxnode))     ! Total
    If (nzfr <= 10*config%mxatms) Then

      Allocate (cfr(1:nzfr),xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'poisson_frzn_forces allocation failure 1'
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
      Call gsum(comm,cfr)
      Call gsum(comm,xfr)
      Call gsum(comm,yfr)
      Call gsum(comm,zfr)

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

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          rsq=xrr**2+yrr**2+zrr**2

          rrr=Sqrt(rsq)
          chgprd=-cfr(ii)*cfr(jj)/eps*r4pie0 ! the MINUS signifies the exclusion!!!

          ! calculate forces

          coul = chgprd/rrr
          fcoul = coul/rsq

          fx = fcoul*xrr
          fy = fcoul*yrr
          fz = fcoul*zrr

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

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          rsq=xrr**2+yrr**2+zrr**2

          rrr=Sqrt(rsq)
          chgprd=-cfr(ii)*cfr(jj)/eps*r4pie0 ! the MINUS signifies the exclusion!!!

          ! calculate forces

          coul = chgprd/rrr
          fcoul = coul/rsq

          fx = fcoul*xrr
          fy = fcoul*yrr
          fz = fcoul*zrr

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

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          rsq=xrr**2+yrr**2+zrr**2

          rrr=Sqrt(rsq)
          chgprd=-cfr(ii)*cfr(jj)/eps*r4pie0 ! the MINUS signifies the exclusion!!!

          ! calculate forces

          coul = chgprd/rrr
          fcoul = coul/rsq

          fx = fcoul*xrr
          fy = fcoul*yrr
          fz = fcoul*zrr

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
        Write(message,'(a,i0)') 'poisson_frzn_forces deallocation failure 1'
        Call error(0,message)
      End If

    Else

      ! We resort to approximating N*(N-1)/2 interactions
      ! with the short-range one from the two body linked cell neigh%list

      Allocate (xxt(1:neigh%max_list),yyt(1:neigh%max_list),zzt(1:neigh%max_list),rrt(1:neigh%max_list), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a,i0)') 'poisson_frzn_forces allocation failure 2'
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
          !           Call images(imcon,cell,limit,xxt,yyt,zzt)

          ! square of distances

          Do k=1,limit
            rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
          End Do

          Do k=1,limit
            j=neigh%list(neigh%list(-1,i)+k,i)

            rrr=rrt(k)
            If (Abs(config%parts(j)%chge) > zero_plus .and. rrr < neigh%cutoff) Then
              chgprd=-config%parts(i)%chge*config%parts(j)%chge/eps*r4pie0 ! the MINUS signifies the exclusion!!!
              rsq=rrr**2

              ! calculate forces

              coul = chgprd/rrr
              fcoul = coul/rsq

              fx = fcoul*xrr
              fy = fcoul*yrr
              fz = fcoul*zrr


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
        Write(message,'(a)') 'poisson_frzn_forces deallocation failure 2'
        Call error(0,message)
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
      Write(message,'(a)') 'poisson_frzn_forces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine poisson_frzn_forces

  Subroutine cleanup(pois)
    Type(poisson_type) :: pois

    If (Allocated(pois%phi)) Then
      Deallocate(pois%phi)
    End If

    If (Allocated(pois%b)) Then
      Deallocate(pois%b)
    End If

    If (Allocated(pois%F)) Then
      Deallocate(pois%F)
    End If

    If (Allocated(pois%s)) Then
      Deallocate(pois%s)
    End If

    If (Allocated(pois%r)) Then
      Deallocate(pois%r)
    End If

    If (Allocated(pois%p)) Then
      Deallocate(pois%p)
    End If

    If (Allocated(pois%r0)) Then
      Deallocate(pois%r0)
    End If

    If (Allocated(pois%v)) Then
      Deallocate(pois%v)
    End If

    If (Allocated(pois%t)) Then
      Deallocate(pois%t)
    End If

    If (Allocated(pois%pattern)) Then
      Deallocate(pois%pattern)
    End If

    If (Allocated(pois%phi0)) Then
      Deallocate(pois%phi0)
    End If
  End Subroutine cleanup
End Module poisson
