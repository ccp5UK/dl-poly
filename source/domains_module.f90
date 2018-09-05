Module domains_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring fundamental domain decomposition variables
! and arrays for the entire package
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
! contrib   - i.j.bush august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer          , Public, Save :: nprx,npry,nprz, &   ! dimensions of the 3D processor/domain grid
                                     idx,idy,idz,    &   ! this domain's coordinates on the grid
                                     map(1:26),mop(1:26) ! neighbourhood coordinates map of the grid
  Real( Kind = wp ), Public, Save :: nprx_r,npry_r,nprz_r, & ! RS = reduced space [-0.5:0.5)^3
                                     r_nprx,r_npry,r_nprz    ! domain's length in RS

  Public  :: map_domains
  Private

! Size of array for factorization. It allows for the first max_factor-1 prime numbers in the
! factorization, so 10 allows for the factors 2, 3, 5, 7, 11, 13, 17, 19, 23

  Integer, Parameter :: max_factor = 10

Contains

  Subroutine map_domains(imcon,wx,wy,wz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for obtaining characteristics of a parallel
! computer and constructing a domain mapping
!
! copyright - daresbury laboratory
! author    - i.t.todorov & i.j.bush september 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use comms_module, Only : idnode,mxnode,gsync

    Implicit None

    Integer,           Intent( In    ) :: imcon
    Real( Kind = wp ), Intent( In    ) :: wx,wy,wz ! MD cell Cartesian widths

    Integer                            :: i,j, jdx,jdy,jdz

! Limits on the (x,y,z) sizes of the processor grid

    Integer                            :: limx,limy,limz

! processor factorisations of and factorised grid vectors over the x and y
! sides of the Cartesian rectangular parallelepiped approximating the MD cell

    Integer                            :: nfacsx, nfacsy
    Integer, Dimension( 1:max_factor ) :: pfacsx, pfacsy

! Target grid product and running value of the (y,z) processor grid product

    Integer                            :: P,Pyz

! Running values of the (x,y,z) sizes of the processor producing P and Pyz

    Integer                            :: nx,ny,nz

! Running values of the (x,y,z) sizes of the grid cell volume and surface
! Minimum surface

    Real( Kind = wp )                  :: dx,dy,dz,S,min_S

! Search tolerance

    Real( Kind = wp ), Parameter       :: tol = 1.0e-6_wp


! DD SEARCH

    If (mxnode == 1) Then

       nprx = 1
       npry = 1
       nprz = 1

    Else ! mxnode > 1

! Impose imcon driven limits on decomposition

       limx = Merge( Huge( 1 ), 2, imcon /= 0 )
       limy = Merge( Huge( 1 ), 2, imcon /= 0 )
       limz = Merge( Huge( 1 ), 2, imcon /= 0 .and. imcon /= 6 )

! Minimum surface (the largest possible) and decomposition (none)

       min_S = Huge( 1.0_wp )
       nprx = -1
       npry = -1
       nprz = -1

       P = mxnode
       Call factor( P, pfacsx )
       nfacsx = get_n_factors( pfacsx )
       Do i = 1, nfacsx
          nx = get_nth_factor( pfacsx, i )
          If ( nx > limx ) Cycle
          dx = wx / Real(nx,wp)

          Pyz = P / nx
          Call factor( Pyz, pfacsy )
          nfacsy = get_n_factors( pfacsy )
          Do j = 1, nfacsy
             ny = get_nth_factor( pfacsy, j )
             If ( ny > limy ) Cycle
             dy = wy / Real(ny,wp)

             nz = Pyz / ny
             If ( nz > limz ) Cycle
             dz = wz / Real(nz,wp)

             S = 2.0_wp * ( dx*dy + dy*dz + dz*dx )

             If      ( min_S - S         > tol ) Then

! Qualifier

                min_S = S
                nprx  = nx
                npry  = ny
                nprz  = nz

             Else If ( Abs( min_S - S ) < tol ) Then

! For degenerate cases choose case where have least procs down any dimension (it should help FFT)

                If ( Max( nx, ny, nz ) < Max( nprx, npry, nprz ) ) Then

                   min_S = S
                   nprx  = nx
                   npry  = ny
                   nprz  = nz

                Else If ( Max( nx, ny, nz ) == Max( nprx, npry, nprz ) ) Then

! The the new case has the same number of procs

                   If ( nx < nprx ) Then

! Choose first the case which has the least down x

                      min_S = S
                      nprx  = nx
                      npry  = ny
                      nprz  = nz

                   Else If ( nx == nprx .and. ny < npry ) Then

! If max the same AND x the same choose it if y is less

                      min_S = S
                      nprx  = nx
                      npry  = ny
                      nprz  = nz

                   End If
                End If
             End If
          End Do
       End Do

       Call gsync()
       If ( nprx == - 1 .or. npry == -1 .or. nprz == -1 ) Call error(520)

    End If

! DD MAPPING

    map=0
    mop=0

! Define npr._r & r_npr.

    nprx_r = Real(nprx,wp) ; r_nprx = 1.0_wp/nprx_r
    npry_r = Real(npry,wp) ; r_npry = 1.0_wp/npry_r
    nprz_r = Real(nprz,wp) ; r_nprz = 1.0_wp/nprz_r

! construct map of neighbouring nodes and domains

    idz=idnode/(nprx*npry)
    idy=idnode/nprx-idz*npry
    idx=Mod(idnode,nprx)

    jdz=nprz+idz
    jdy=npry+idy
    jdx=nprx+idx

    map(1)=idcube(Mod(jdx-1,nprx),idy,idz)
    map(2)=idcube(Mod(idx+1,nprx),idy,idz)
    map(3)=idcube(idx,Mod(jdy-1,npry),idz)
    map(4)=idcube(idx,Mod(idy+1,npry),idz)
    map(5)=idcube(idx,idy,Mod(jdz-1,nprz))
    map(6)=idcube(idx,idy,Mod(idz+1,nprz))

! map(1) points to the node number that is in -x direction to me (idnode)
! map(2) points to the node number that is in +x direction to me (idnode)
! map(3) points to the node number that is in -y direction to me (idnode)
! map(4) points to the node number that is in +y direction to me (idnode)
! map(5) points to the node number that is in -z direction to me (idnode)
! map(6) points to the node number that is in +z direction to me (idnode)

    map(7)=idcube(Mod(jdx-1,nprx),Mod(idy+1,npry),idz)
    map(8)=idcube(Mod(idx+1,nprx),Mod(jdy-1,npry),idz)
    map(9)=idcube(Mod(jdx-1,nprx),Mod(jdy-1,npry),idz)
    map(10)=idcube(Mod(idx+1,nprx),Mod(idy+1,npry),idz)

    map(11)=idcube(Mod(jdx-1,nprx),idy,Mod(idz+1,nprz))
    map(12)=idcube(Mod(idx+1,nprx),idy,Mod(jdz-1,nprz))
    map(13)=idcube(Mod(jdx-1,nprx),idy,Mod(jdz-1,nprz))
    map(14)=idcube(Mod(idx+1,nprx),idy,Mod(idz+1,nprz))

    map(15)=idcube(idx,Mod(jdy-1,npry),Mod(idz+1,nprz))
    map(16)=idcube(idx,Mod(idy+1,npry),Mod(jdz-1,nprz))
    map(17)=idcube(idx,Mod(jdy-1,npry),Mod(jdz-1,nprz))
    map(18)=idcube(idx,Mod(idy+1,npry),Mod(idz+1,nprz))

    map(19)=idcube(Mod(jdx-1,nprx),Mod(jdy-1,npry),Mod(jdz-1,nprz))
    map(20)=idcube(Mod(idx+1,nprx),Mod(idy+1,npry),Mod(idz+1,nprz))
    map(21)=idcube(Mod(jdx-1,nprx),Mod(jdy-1,npry),Mod(idz+1,nprz))
    map(22)=idcube(Mod(idx+1,nprx),Mod(idy+1,npry),Mod(jdz-1,nprz))

    map(23)=idcube(Mod(jdx-1,nprx),Mod(idy+1,npry),Mod(jdz-1,nprz))
    map(24)=idcube(Mod(idx+1,nprx),Mod(jdy-1,npry),Mod(idz+1,nprz))
    map(25)=idcube(Mod(jdx-1,nprx),Mod(idy+1,npry),Mod(idz+1,nprz))
    map(26)=idcube(Mod(idx+1,nprx),Mod(jdy-1,npry),Mod(jdz-1,nprz))

! Determine which processors appear more than once
! (mop(the first unique node)=0 if repeated then =1)
! Up to Min(mxnode-1,26) elements of mop can be zero
! the remainder is images of idnode
!
! NEEDED FOR CATCHING SELF-HALOING

    Do i=1,26
       If (idnode == map(i)) mop(i)=1
       Do j=i+1,26
          If (map(i) == map(j)) mop(j)=1
       End Do
    End Do

    Contains

    Subroutine factor( n, facs )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 prime factorisability decomposition function
! (as in parallel_fft)
!
! copyright - daresbury laboratory
! author    - i.j.bush august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Implicit None

      Integer                , Intent( In    ) :: n
      Integer, Dimension( : ), Intent(   Out ) :: facs

      Integer :: left
      Integer :: p
      Integer :: i

      facs = 0

      left = n
      Do i = 1, Size( facs ) - 1
         p = get_nth_prime( i )

         If ( p <= 0 ) Exit

         Do While ( p * ( left / p ) == left )
            left = left / p

            facs( i ) = facs( i ) + 1
         End Do
      End Do

      facs( Size( facs ) ) = left

    End Subroutine factor

    Function get_nth_prime( n )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 return n-th prime function (as in parallel_fft)
!
! copyright - daresbury laboratory
! author    - i.j.bush august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Implicit None

      Integer                  :: get_nth_prime

      Integer, Intent( In    ) :: n

      Integer, Dimension( 1:170 ), Parameter :: primes = (/                             &
             2,      3,      5,      7,     11,     13,     17,     19,     23,     29, &
            31,     37,     41,     43,     47,     53,     59,     61,     67,     71, &
            73,     79,     83,     89,     97,    101,    103,    107,    109,    113, &
           127,    131,    137,    139,    149,    151,    157,    163,    167,    173, &
           179,    181,    191,    193,    197,    199,    211,    223,    227,    229, &
           233,    239,    241,    251,    257,    263,    269,    271,    277,    281, &
           283,    293,    307,    311,    313,    317,    331,    337,    347,    349, &
           353,    359,    367,    373,    379,    383,    389,    397,    401,    409, &
           419,    421,    431,    433,    439,    443,    449,    457,    461,    463, &
           467,    479,    487,    491,    499,    503,    509,    521,    523,    541, &
           547,    557,    563,    569,    571,    577,    587,    593,    599,    601, &
           607,    613,    617,    619,    631,    641,    643,    647,    653,    659, &
           661,    673,    677,    683,    691,    701,    709,    719,    727,    733, &
           739,    743,    751,    757,    761,    769,    773,    787,    797,    809, &
           811,    821,    823,    827,    829,    839,    853,    857,    859,    863, &
           877,    881,    883,    887,    907,    911,    919,    929,    937,    941, &
           947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013 /)

      If ( n <= Size( primes ) ) Then
         get_nth_prime = primes( n )
      Else
         get_nth_prime = -1
      End If

    End Function get_nth_prime

    Function get_n_factors( factors ) Result( nfacs )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 function to return the total number of
! possible integer factorisations
!
! copyright - daresbury laboratory
! author    - i.j.bush august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Implicit None

      Integer                                  :: nfacs

      Integer, Dimension( : ), Intent( In    ) :: factors

      nfacs = Product( factors( 1:Size( factors ) - 1 ) + 1 )

      If ( factors( Size( factors ) ) /= 1 ) nfacs = nfacs * 2

    End function get_n_factors

    Function get_nth_factor( factors, n )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 function to return the n-th factoriser from the list of
! all possible integer factorisers
!
! copyright - daresbury laboratory
! author    - i.j.bush august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Implicit None

      Integer                                  :: get_nth_factor

      Integer, Dimension( : ), Intent( In    ) :: factors
      Integer                , Intent( In    ) :: n

      Integer, Dimension( 1:Size( factors ) ) :: fac_counts

      Integer :: nfacs, nt
      Integer :: dim_prod
      Integer :: i

      nfacs = get_n_factors( factors )
      If ( n > nfacs ) Then

         get_nth_factor = -1

      Else

         If ( factors( Size( factors ) ) /= 1 .and. n > nfacs / 2 ) Then
            nt = n - nfacs / 2
         Else
            nt = n
         End If
         nt = nt - 1

         dim_prod = Product( factors( 1:Size( factors ) - 2 ) + 1 )
         Do i = Size( factors ) - 1, 2, -1
            fac_counts( i ) = nt / dim_prod
            nt = nt - fac_counts( i ) * dim_prod
            dim_prod = dim_prod / ( factors( i - 1 ) + 1 )
         End Do
         fac_counts( 1 ) = nt

         get_nth_factor = 1
         Do i = 1, Size( factors ) - 1
            get_nth_factor = get_nth_factor * ( get_nth_prime( i ) ** fac_counts( i ) )
         End Do

         If ( factors( Size( factors ) ) /= 1 .and. n > nfacs / 2 ) &
            get_nth_factor = get_nth_factor * factors( Size( factors ) )

      End If

    End Function get_nth_factor

    Function idcube(i,j,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 hypercube mapping function
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2006
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Implicit None

      Integer                  :: idcube

      Integer, Intent( In    ) :: i,j,k

      idcube = i + nprx * ( j + npry * k )

    End Function idcube

  End Subroutine map_domains

End Module domains_module
