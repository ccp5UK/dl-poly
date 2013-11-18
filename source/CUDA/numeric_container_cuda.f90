!!!!!!!!!!!!!!!!!!!!!!!! THIS IS NUMERIC_CONTAINER !!!!!!!!!!!!!!!!!!!!!
!
! Function uni - generates a random number
!
! Subroutine box_mueller - generates gaussian random numbers of unit
!                          variance (with zero mean and standard
!                          variation of 1)
!
! Subroutine gauss_old - constructs velocity arrays with a gaussian
!                        distribution of unit variance (zero mean) by
!                        an approximation of the Central Limit Theorem
!
! Subroutine gauss - constructs velocity arrays with a gaussian
!                    distribution of unit variance (zero mean) using
!                    the box-mueller method
!
! Subroutine erfcgen - generates interpolation tables for erfc and its
!                      derivative
!
! Function match - determines a match between integer value 'n' and an
!                  array of integers in ascending order
!
! Subroutine shellsort - sorts an integer array in ascending order
!
! Subroutine shellsort2 - sorts an integer array in ascending order,
!                         keeping the original ranking of the array
!
! Function local_index - finds the local atom number given the global
!                        atom number
!
! Subroutine dcell - calculates the dimensional properies of a
!                    simulation cell
!
! Subroutine invert - calculates the invert of a 3x3 matrix using
!                     cofactors
!
! Subroutine images - calculates the minimum image distance of
!                     atom pairs within a specified MD cell
!
! Subroutine pbcshift - calculates the minimum image of atoms within
!                       a specified MD cell in accordance with the DD
!                       boundary convention
!
! Subroutine jacobi - diagonalises real symmetric matrices by the
!                     Jacobi method
!
! Subroutine mat_mul - calculates product of two 3x3 matrices written
!                     in a DL_POLY format as vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function uni()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 random number generator based on the universal random number
! generator of marsaglia, zaman and tsang
! (stats and prob. lett. 8 (1990) 35-39.)
!
! This random number generator originally appeared in "Toward a
! Universal Random Number Generator" by George Marsaglia, Arif Zaman and
! W.W. Tsang in Florida State University Report: FSU-SCRI-87-50 (1987).
! It was later modified by F. James and published in "A Review of
! Pseudo-random Number Generators".
! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
! It passes ALL of the tests for random number generators and has a
! period of 2^144, is completely portable (gives bit identical results
! on all machines with at least 24-bit mantissas in the floating point
! representation).
! The algorithm is a combination of a Fibonacci sequence (with lags of
! 97 and 33, and operation "subtraction plus one, modulo one") and an
! "arithmetic sequence" (using subtraction).
! Use IJ = 1802 & KL = 9373 (idnode=0) to test the random number
! generator. The subroutine RANMAR should be used to generate 20000
! random numbers.  Then display the next six random numbers generated
! multiplied by 4096*4096.  If the random number generator is working
! properly, the random numbers should be:
!         6533892.0  14220222.0  7275067.0
!         6172232.0  8354498.0   10633180.0
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov april 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode
  Use setup_module, Only : lseed,seed

  Implicit None

  Logical,           Save :: newjob = .true.
  Integer,           Save :: ir,jr
  Integer                 :: i,ii,ij,j,jj,k,kl,l,m
  Real( Kind = wp ), Save :: c,cd,cm,u(1:97)
  Real( Kind = wp )       :: s,t,uni

! initialise parameters u,c,cd,cm

  If (newjob .or. lseed) Then
     newjob = .false.

! If no seeding is specified then default to DL_POLY scheme

     If (lseed) Then

        lseed=.false.

! First random number seed must be between 0 and 31328
! Second seed must have a value between 0 and 30081

        ij=Mod(Abs(seed(1)+idnode),31328)
        i = Mod(ij/177,177) + 2;
        j = Mod(ij,177)     + 2;

        kl=Mod(Abs(seed(2)+idnode),30081)
        k = Mod(kl/169,178) + 1
        l = Mod(kl,169)

     Else

! initial values of i,j,k must be in range 1 to 178 (not all 1)
! initial value of l must be in range 0 to 168

        i = Mod(idnode,166) + 12
        j = Mod(idnode,144) + 34
        k = Mod(idnode,122) + 56
        l = Mod(idnode,90)  + 78

     End If

     ir = 97
     jr = 33

     Do ii=1,97

        s = 0.0_wp
        t = 0.5_wp

        Do jj=1,24

           m = Mod(Mod(i*j,179)*k,179)
           i = j
           j = k
           k = m
           l = Mod(53*l+1,169)
           If (Mod(l*m,64) >= 32) s = s+t
           t = 0.5_wp*t

        End Do

        u(ii)=s

     End Do

     c  =   362436.0_wp/16777216.0_wp
     cd =  7654321.0_wp/16777216.0_wp
     cm = 16777213.0_wp/16777216.0_wp

  End If

! calculate random number

  uni=u(ir)-u(jr)
  If (uni < 0.0_wp) uni = uni + 1.0_wp

  u(ir)=uni

  ir=ir-1
  If (ir == 0) ir = 97

  jr=jr-1
  If (jr == 0) jr = 97

  c = c-cd
  If (c < 0.0_wp) c = c+cm

  uni = uni-c
  If (uni < 0.0_wp) uni = uni + 1.0_wp

End Function uni

Subroutine box_mueller(gauss1,gauss2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine using the box-mueller method for generating
! gaussian random numbers of unit variance (with zero mean and standard
! variation of 1).  Otherwise, an approximation of the Central Limit
! Theorem must be used: G = (1/A)*[Sum_i=1,N(Ri) - AN/2]*(12/N)^(1/2),
! where A is the number of outcomes from the random throw Ri and N is
! the number of tries.
!
! dependent on uni
!
! copyright - daresbury laboratory
! author    - w.smith may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Intent(   Out ) :: gauss1,gauss2

  Logical           :: newjob = .true.
  Real( Kind = wp ) :: uni,ran0,ran1,ran2

! make sure uni is initialised

  If (newjob) Then
     newjob = .false.
     ran0=uni()
  End If

  ran0=1.0_wp

! generate uniform random numbers on [-1, 1)

  Do While (ran0 >= 1.0_wp)
     ran1=2.0_wp*uni()-1.0_wp
     ran2=2.0_wp*uni()-1.0_wp
     ran0=ran1**2+ran2**2
  End Do

! calculate gaussian random numbers

  ran0=Sqrt(-2.0_wp*Log(ran0)/ran0)
  gauss1=ran0*ran1
  gauss2=ran0*ran2

End Subroutine box_mueller

Subroutine gauss_old(natms,vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for constructing velocity arrays with a gaussian
! distribution of unit variance (zero mean), based on the method
! described by Allen and Tildesley in "Computer Simulation of Liquids",
! Clarendon Press 1987, P347.  It is based on an approximation of the
! Central Limit Theorem : G = (1/A)*[Sum_i=1,N(Ri) - AN/2]*(12/N)^(1/2),
! where A is the number of outcomes from the random throw Ri and N is
! the number of tries.
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov july 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Parameter :: a1 = 3.949846138_wp
  Real( Kind = wp ), Parameter :: a3 = 0.252408784_wp
  Real( Kind = wp ), Parameter :: a5 = 0.076542912_wp
  Real( Kind = wp ), Parameter :: a7 = 0.008355968_wp
  Real( Kind = wp ), Parameter :: a9 = 0.029899776_wp

  Integer,                             Intent( In    ) :: natms
  Real( Kind = wp ), Dimension( 1:* ), Intent(   Out ) :: vxx,vyy,vzz

  Integer           :: i,j
  Real( Kind = wp ) :: uni,rrr,rr2

  Do i=1,natms
     rrr=0.0_wp
     Do j=1,12
        rrr=rrr+uni()
     End Do
     rrr=(rrr-6.0_wp)/4.0_wp
     rr2=rrr*rrr
     vxx(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))

     rrr=0.0_wp
     Do j=1,12
        rrr=rrr+uni()
     End Do
     rrr=(rrr-6.0_wp)/4.0_wp
     rr2=rrr*rrr
     vyy(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))

     rrr=0.0_wp
     Do j=1,12
        rrr=rrr+uni()
     End Do
     rrr=(rrr-6.0_wp)/4.0_wp
     rr2=rrr*rrr
     vzz(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))
  End Do

End Subroutine gauss_old

Subroutine gauss(natms,vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for constructing velocity arrays with a gaussian
! distribution of unit variance (zero mean), based on the box-mueller
! method
!
! copyright - daresbury laboratory
! author    - w.smith july 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Integer,                             Intent( In    ) :: natms
  Real( Kind = wp ), Dimension( 1:* ), Intent(   Out ) :: vxx,vyy,vzz

  Integer           :: i,j
  Real( Kind = wp ) :: gauss1,gauss2

  Do i=1,(natms+1)/2
     j=natms+1-i

     Call box_mueller(gauss1,gauss2)
     vxx(i)=gauss1
     vxx(j)=gauss2

     Call box_mueller(gauss1,gauss2)
     vyy(i)=gauss1
     vyy(j)=gauss2

     Call box_mueller(gauss1,gauss2)
     vzz(i)=gauss1
     vzz(j)=gauss2
  End Do

End Subroutine gauss


Subroutine erfcgen_helper(rcut,alpha,mxgrid,erc,fer)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for generating interpolation tables for erfc and its
! derivative - for use with Ewald sum
!
! copyright - daresbury laboratory
! author    - t.forester december 1994
! amended   - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : sqrpi

#ifdef _OPENMP
  Use omp_lib
#endif

  Implicit None

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

  Integer,                                  Intent( In    ) :: mxgrid
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha
  Real( Kind = wp ), Dimension( 1:mxgrid ), Intent(   Out ) :: erc,fer

  Integer           :: i,j,ib,ie
  Real( Kind = wp ) :: drewd,exp1,rrr,rsq,tt

! look-up tables for real space part of ewald sum

  drewd = rcut/Real(mxgrid-4,wp)

#ifdef _OPENMP
  i = 1 + omp_get_thread_num()
#else
  i = 1
#endif

  Do While (i<=mxgrid)
     rrr = Real(i,wp)*drewd
     rsq = rrr*rrr

     tt = 1.0_wp/(1.0_wp + pp*alpha*rrr)
     exp1 = Exp(-(alpha*rrr)**2)

     erc(i) = tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rrr
     fer(i) = (erc(i) + 2.0_wp*(alpha/sqrpi)*exp1)/rsq

#ifdef _OPENMP
     i = i + omp_get_num_threads()
#else
     i = i + 1
#endif

  End Do
  !$OMP BARRIER
End Subroutine erfcgen_helper


Subroutine erfcgen(rcut,alpha,mxgrid,erc,fer)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for generating interpolation tables for erfc and its
! derivative - for use with Ewald sum
!
! copyright - daresbury laboratory
! author    - t.forester december 1994
! amended   - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : sqrpi

  Implicit None

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

  Integer,                                  Intent( In    ) :: mxgrid
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha
  Real( Kind = wp ), Dimension( 1:mxgrid ), Intent(   Out ) :: erc,fer

  Integer           :: i
  Real( Kind = wp ) :: drewd,exp1,rrr,rsq,tt

! look-up tables for real space part of ewald sum

  drewd = rcut/Real(mxgrid-4,wp)

  Do i=1,mxgrid
     rrr = Real(i,wp)*drewd
     rsq = rrr*rrr

     tt = 1.0_wp/(1.0_wp + pp*alpha*rrr)
     exp1 = Exp(-(alpha*rrr)**2)

     erc(i) = tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rrr
     fer(i) = (erc(i) + 2.0_wp*(alpha/sqrpi)*exp1)/rsq
  End Do

End Subroutine erfcgen

Function match(n,ind_top,list)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 function to determine a match between a positive integer
! 'n' and an array of positive integer 'list(1:ind_top)' sorted in
! ascending order
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2006
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer, Intent( In    ) :: n,ind_top,list(1:*)

  Logical :: match
  Integer :: ind_old,ind_now

  If (n < 1) Call error(0)

  match = .false.

  If (ind_top < 1) Return

  ind_old = 1
  ind_now = 1

  Do
     If      (n == list(ind_now)) Then
        match=.true.
        Return
     Else If (n >  list(ind_now)) Then
        If (ind_old == ind_top) Return
        ind_old = ind_now
        ind_now = (ind_old+ind_top+1)/2
     Else If (n <  list(ind_now)) Then
        ind_now = (ind_old+ind_now)/2
        If (ind_now == ind_old) Return
     End If
  End Do

End Function match

Subroutine shellsort(n,list)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 shell sort routine.  Sorts an array of integers into
! ascending order.
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer,                   Intent( In    ) :: n
  Integer, Dimension( 1:* ), Intent( InOut ) :: list

  Logical :: go
  Integer :: nl,index,value,i

! set up sort

  If (n > 1) Then

! number of lists

     nl = n/2

! iterate shell sort until there is a list

     Do While (nl > 0)

! for all lists from next-to-ground-level up to their end

        Do i=nl+1,n

           value = list(i)
           index = i

! Antibubble down between levels of the same list

           go = .true.
           Do While (index > nl .and. go)
              go = (list(index-nl) > value)

              If (go) Then
                 list(index) = list(index-nl)
                 index = index-nl
              End If
           End Do

! Last insertion as close to the ground as it gets

           list(index) = value

        End Do

! Decrease the number of lists

        nl = nl/2

     End Do

  End If

End Subroutine shellsort

Subroutine shellsort2(n,rank,list)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 shell sort routine.  Sorts an array of integers (list) into
! ascending order.  The original rank of array list is kept in rank.
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer,                   Intent( In    ) :: n
  Integer, Dimension( 1:* ), Intent( InOut ) :: list,rank

  Logical :: go
  Integer :: nl,index,value,i,rang

! set up sort

  If (n > 1) Then

! number of lists

     nl = n/2

! iterate shell sort until there is a list

     Do While (nl > 0)

! for all lists from next-to-ground-level up to their end

        Do i=nl+1,n

           value = list(i)
           rang  = rank(i)
           index = i

! Antibubble down between levels of the same list

           go = .true.
           Do While (index > nl .and. go)
              go = (list(index-nl) > value)

              If (go) Then
                 list(index) = list(index-nl)
                 rank(index) = rank(index-nl)
                 index = index-nl
              End If
           End Do

! Last insertion as close to the ground as it gets

           list(index) = value
           rank(index) = rang

        End Do

! Decrease the number of lists

        nl = nl/2

     End Do

  End If

End Subroutine shellsort2

Function local_index(global_index,search_limit,rank,list)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 function to find the local atom number given the
! global atom number.  (For use with DD codes only)
!
! If multiple copies present it returns the lowest local atom number
! If no copy is present it returns zero
!
! rank(1,*) - array of local atom indices, ranking list
! list(1,*) - array of sorted global atom indices
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer,                   Intent( In    ) :: global_index,search_limit
  Integer, Dimension( 1:* ), Intent( In    ) :: rank,list

  Integer local_index,point,lower_bound,upper_bound,down

! Initialise

  local_index = 0

! Limits for the search

  lower_bound=1
  upper_bound=search_limit

! Get smart, check for exceptions (whether it's a false pass)
! and check for a match on bounds

  If      (global_index <= 0) Then

     local_index = 0

     Return

  Else If (global_index == list(lower_bound)) Then

     local_index = rank(lower_bound)

     Return

  Else If (global_index == list(search_limit)) Then

     down = search_limit

     Do While (global_index == list(down) .and. down >= lower_bound)

        local_index = rank(down)
        down = down - 1

     End Do

     Return

  End If

! Carry on then

  Do While (upper_bound-lower_bound > 1)

     point=(lower_bound+upper_bound)/2

     If (global_index < list(point)) Then

        upper_bound = point

     Else If (global_index > list(point)) Then

        lower_bound = point

     Else

        down=point

        Do While (global_index == list(down))

           local_index = rank(down)
           down = down - 1

        End Do

        Return

     End If

  End Do

End Function local_index

Subroutine dcell(aaa,bbb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to calculate the dimensional properties of a
! simulation cell specified by the input 3x3 matrix aaa (cell vectors in
! rows, the matrix is in the form of one dimensional reading
! (row1,row2,row3).
!
! The results are returned in the array bbb, with:
!
! bbb(1 to 3) - lengths of cell vectors: a(x,y,z) , b(x,y,z) , c(x,y,z)
! bbb(4 to 6) - cosines of cell angles: gamma(a,c) , alpha(b,c) , beta(a,c)
! bbb(7 to 9) - perpendicular cell widths : wx(y,z) , wy(x,z) , wz(x,y)
! bbb(10)     - cell volume
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Dimension( 1:9 ),  Intent( In    ) :: aaa
  Real( Kind = wp ), Dimension( 1:10 ), Intent(   Out ) :: bbb

  Real( Kind = wp ) :: axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3, &
                       x(1:3),y(1:3),z(1:3),d(1:3)

! calculate lengths of cell vectors

  bbb(1)=Sqrt(aaa(1)**2+aaa(2)**2+aaa(3)**2)
  bbb(2)=Sqrt(aaa(4)**2+aaa(5)**2+aaa(6)**2)
  bbb(3)=Sqrt(aaa(7)**2+aaa(8)**2+aaa(9)**2)

! calculate cosines of cell angles

  bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
  bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
  bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))

! calculate vector products of cell vectors

  axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
  axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
  axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)

  bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
  bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
  bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)

  cxa1=aaa(8)*aaa(3)-aaa(9)*aaa(2)
  cxa2=aaa(9)*aaa(1)-aaa(7)*aaa(3)
  cxa3=aaa(7)*aaa(2)-aaa(8)*aaa(1)

! calculate volume of cell

  bbb(10)=Abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)

! calculate cell perpendicular widths

  d(1)=bbb(10)/Sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
  d(2)=bbb(10)/Sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
  d(3)=bbb(10)/Sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

  x(1)=Abs(aaa(1))/bbb(1) ; y(1)=Abs(aaa(2))/bbb(1) ; z(1)=Abs(aaa(3))/bbb(1)
  x(2)=Abs(aaa(4))/bbb(2) ; y(2)=Abs(aaa(5))/bbb(2) ; z(2)=Abs(aaa(6))/bbb(2)
  x(3)=Abs(aaa(7))/bbb(3) ; y(3)=Abs(aaa(8))/bbb(3) ; z(3)=Abs(aaa(9))/bbb(3)

! distribute widths

  If      (x(1) >= x(2) .and. x(1) >= x(3)) Then
     bbb(7)=d(1)
     If (y(2) >= y(3)) Then
        bbb(8)=d(2)
        bbb(9)=d(3)
     Else
        bbb(8)=d(3)
        bbb(9)=d(2)
     End If
  Else If (x(2) >= x(1) .and. x(2) >= x(3)) Then
     bbb(7)=d(2)
     If (y(1) >= y(3)) Then
        bbb(8)=d(1)
        bbb(9)=d(3)
     Else
        bbb(8)=d(3)
        bbb(9)=d(1)
     End If
  Else
     bbb(7)=d(3)
     If (y(1) >= y(2)) Then
        bbb(8)=d(1)
        bbb(9)=d(2)
     Else
        bbb(8)=d(2)
        bbb(9)=d(1)
     End If
  End If

End Subroutine dcell

Subroutine invert(a,b,d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to invert a 3x3 matrix using cofactors
! matrices are in the form of one dimensional array reading
! (row1,row2,row3)
!
! a - input matrix
! b - inverted matrix
! d - determinant
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Dimension( 1:9 ), Intent( In    ) :: a
  Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: b
  Real( Kind = wp ),                   Intent(   Out ) :: d

  Real( Kind = wp ) :: r

! calculate adjoint matrix

  b(1)=a(5)*a(9)-a(6)*a(8)
  b(2)=a(3)*a(8)-a(2)*a(9)
  b(3)=a(2)*a(6)-a(3)*a(5)
  b(4)=a(6)*a(7)-a(4)*a(9)
  b(5)=a(1)*a(9)-a(3)*a(7)
  b(6)=a(3)*a(4)-a(1)*a(6)
  b(7)=a(4)*a(8)-a(5)*a(7)
  b(8)=a(2)*a(7)-a(1)*a(8)
  b(9)=a(1)*a(5)-a(2)*a(4)

! calculate determinant

  d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
  r=0.0_wp
  If (Abs(d) > 0.0_wp) r=1.0_wp/d

! complete inverse matrix

  b(1)=r*b(1)
  b(2)=r*b(2)
  b(3)=r*b(3)
  b(4)=r*b(4)
  b(5)=r*b(5)
  b(6)=r*b(6)
  b(7)=r*b(7)
  b(8)=r*b(8)
  b(9)=r*b(9)

End Subroutine invert


! Based on 'images' and to be used by two_body_forces_helper
! for OpenMP-based parallelism in an orphaned manner.
Subroutine images_helper_orphan(imcon,cell,pairs,xxx,yyy,zzz)
  Use kinds_f90
  Use setup_module, Only : rt2,rt3

  Implicit None

  Integer,                              Intent( In    ) :: imcon,pairs
  Real( Kind = wp ), Dimension( 1:9 ),  Intent( In    ) :: cell
  Real( Kind = wp ), Dimension( 1:* ),  Intent( InOut ) :: xxx,yyy,zzz

  Integer           :: i
  Real( Kind = wp ) :: aaa,bbb,ccc,ddd,det,rcell(1:9), &
                       xss,yss,zss

  If (imcon == 2) Then

     ! rectangular (slab) boundary conditions

     aaa=1.0_wp/cell(1)
     bbb=1.0_wp/cell(5)
     ccc=1.0_wp/cell(9)
!$OMP DO
     Do i=1,pairs
        xxx(i)=xxx(i)-cell(1)*Anint(aaa*xxx(i))
        yyy(i)=yyy(i)-cell(5)*Anint(bbb*yyy(i))
        zzz(i)=zzz(i)-cell(9)*Anint(ccc*zzz(i))
     End Do
!$OMP END DO NOWAIT
  Else If (imcon == 3) Then

! parallelepiped boundary conditions

     Call invert(cell,rcell,det)
!$OMP DO
     Do i=1,pairs
        xss=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        yss=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        zss=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

        xss=xss-Anint(xss)
        yss=yss-Anint(yss)
        zss=zss-Anint(zss)

        xxx(i)=cell(1)*xss+cell(4)*yss+cell(7)*zss
        yyy(i)=cell(2)*xss+cell(5)*yss+cell(8)*zss
        zzz(i)=cell(3)*xss+cell(6)*yss+cell(9)*zss
     End Do
!$OMP END DO NOWAIT
  End If
End Subroutine images_helper_orphan

Subroutine images(imcon,cell,pairs,xxx,yyy,zzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the minimum image vector of
! atom pairs within a specified MD cell.  The cell matrix is in the form
! of one dimensional array reading (row1,row2,row3).
!
! Image conditions
!
! imcon=0 no boundary conditions apply
! imcon=1 standard cubic boundaries apply
! imcon=2 orthorhombic boundaries apply
! imcon=3 parallelepiped boundaries apply
! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
!
! Note: in all cases the centre of the MD cell is at (0,0,0)
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : rt2,rt3

  Implicit None

  Integer,                              Intent( In    ) :: imcon,pairs
  Real( Kind = wp ), Dimension( 1:9 ),  Intent( In    ) :: cell
  Real( Kind = wp ), Dimension( 1:* ),  Intent( InOut ) :: xxx,yyy,zzz

  Integer           :: i
  Real( Kind = wp ) :: aaa,bbb,ccc,ddd,det,rcell(1:9), &
                       xss,yss,zss

  If (imcon == 1) Then

! standard cubic boundary conditions

     aaa=1.0_wp/cell(1)

     Do i=1,pairs
        xxx(i)=xxx(i)-cell(1)*Anint(aaa*xxx(i))
        yyy(i)=yyy(i)-cell(1)*Anint(aaa*yyy(i))
        zzz(i)=zzz(i)-cell(1)*Anint(aaa*zzz(i))
     End Do

  Else If (imcon == 2) Then

! rectangular (slab) boundary conditions

     aaa=1.0_wp/cell(1)
     bbb=1.0_wp/cell(5)
     ccc=1.0_wp/cell(9)

     Do i=1,pairs
        xxx(i)=xxx(i)-cell(1)*Anint(aaa*xxx(i))
        yyy(i)=yyy(i)-cell(5)*Anint(bbb*yyy(i))
        zzz(i)=zzz(i)-cell(9)*Anint(ccc*zzz(i))
     End Do

  Else If (imcon == 3) Then

! parallelepiped boundary conditions

     Call invert(cell,rcell,det)

     Do i=1,pairs
        xss=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        yss=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        zss=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

        xss=xss-Anint(xss)
        yss=yss-Anint(yss)
        zss=zss-Anint(zss)

        xxx(i)=cell(1)*xss+cell(4)*yss+cell(7)*zss
        yyy(i)=cell(2)*xss+cell(5)*yss+cell(8)*zss
        zzz(i)=cell(3)*xss+cell(6)*yss+cell(9)*zss
     End Do

  Else If (imcon == 4) Then

! truncated octahedral boundary conditions

     If (.not.(Abs(cell(1)-cell(5)) < 1.0e-6_wp .and. Abs(cell(5)-cell(9)) < 1.0e-6_wp)) Call error(130)

     aaa=1.0_wp/cell(1)

     Do i=1,pairs
        xxx(i)=xxx(i)-cell(1)*Anint(aaa*xxx(i))
        yyy(i)=yyy(i)-cell(1)*Anint(aaa*yyy(i))
        zzz(i)=zzz(i)-cell(1)*Anint(aaa*zzz(i))

        If ((Abs(xxx(i))+Abs(yyy(i))+Abs(zzz(i))) >= 0.75_wp*cell(1)) Then
           xxx(i)=xxx(i)-0.5_wp*Sign(cell(1),xxx(i))
           yyy(i)=yyy(i)-0.5_wp*Sign(cell(1),yyy(i))
           zzz(i)=zzz(i)-0.5_wp*Sign(cell(1),zzz(i))
        End If
     End Do

  Else If (imcon == 5) Then

! rhombic Dodecahedral boundary conditions

     If (.not.(Abs(cell(1)-cell(5)) < 1.0e-6_wp .and. Abs(cell(9)-cell(1)*rt2) < 1.0e-6_wp)) Call error(140)

     aaa=1.0_wp/cell(1)
     bbb=1.0_wp/cell(9)

     Do i=1,pairs
        xxx(i)=xxx(i)-cell(1)*Anint(aaa*xxx(i))
        yyy(i)=yyy(i)-cell(1)*Anint(aaa*yyy(i))
        zzz(i)=zzz(i)-cell(9)*Anint(bbb*zzz(i))

        If ((Abs(xxx(i))+Abs(yyy(i))+Abs(rt2*zzz(i))) >= cell(1)) Then
           xxx(i)=xxx(i)-0.5_wp*Sign(cell(1),xxx(i))
           yyy(i)=yyy(i)-0.5_wp*Sign(cell(1),yyy(i))
           zzz(i)=zzz(i)-0.5_wp*Sign(cell(9),zzz(i))
        End If
     End Do

  Else If (imcon == 6) Then

! x-y boundary conditions

     det=cell(1)*cell(5)-cell(2)*cell(4)

     If (Abs(det) < 1.0e-6_wp) Call error(120)

     det=1.0_wp/det

     rcell(1) =  det*cell(5)
     rcell(2) = -det*cell(2)
     rcell(4) = -det*cell(4)
     rcell(5) =  det*cell(1)

     Do i=1,pairs
        xss=rcell(1)*xxx(i)+rcell(4)*yyy(i)
        yss=rcell(2)*xxx(i)+rcell(5)*yyy(i)

        xss=xss-Anint(xss)
        yss=yss-Anint(yss)

        xxx(i)=cell(1)*xss+cell(4)*yss
        yyy(i)=cell(2)*xss+cell(5)*yss
     End Do

  Else If (imcon == 7) Then

! hexagonal prism boundary conditions

     If (Abs(cell(1)-rt3*cell(5)) > 1.0e-6_wp) Call error(135)

     aaa=cell(1)/(rt3*2.0_wp)
     bbb=cell(1)/rt3
     ccc=rt3/cell(1)
     ddd=1.0_wp/cell(9)

     Do i=1,pairs
        yyy(i)=yyy(i)-bbb*Anint(ccc*yyy(i))
        zzz(i)=zzz(i)-cell(9)*Anint(ddd*zzz(i))

        If ((Abs(yyy(i))+Abs(rt3*xxx(i))) >= bbb) Then
           xxx(i)=xxx(i)-rt3*Sign(aaa,xxx(i))
           yyy(i)=yyy(i)-Sign(aaa,yyy(i))
        End If
    End Do

  End If

End Subroutine images

Subroutine pbcshift(imcon,cell,natms,xxx,yyy,zzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the minimum image of atoms within
! a specified MD cell in accordance with the domain decomposition
! boundary convention for fractional coordinates: every coordinate must
! be intervalled as [-0.5,+0.5)
!
! Note: in all cases the centre of the MD cell is at (0,0,0)
!
! imcon=0 no boundary conditions apply
! imcon=1 standard cubic boundaries apply
! imcon=2 orthorhombic boundaries apply
! imcon=3 parallelepiped boundaries apply
! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2003
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : rt2,rt3,half_minus

  Implicit None

  Integer,                              Intent( In    ) :: imcon,natms
  Real( Kind = wp ), Dimension( 1:9 ),  Intent( In    ) :: cell
  Real( Kind = wp ), Dimension( 1:* ),  Intent( InOut ) :: xxx,yyy,zzz

  Integer           :: i
  Real( Kind = wp ) :: aaa,bbb,ccc,ddd,det,rcell(1:9), &
                       xss,yss,zss

  If (imcon == 1) Then

! standard cubic boundary conditions

     aaa=1.0_wp/cell(1)

     Do i=1,natms
        xss=aaa*xxx(i)
        yss=aaa*yyy(i)
        zss=aaa*zzz(i)

        xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
        zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

        xxx(i)=cell(1)*xss
        yyy(i)=cell(1)*yss
        zzz(i)=cell(1)*zss
     End Do

  Else If (imcon == 2) Then

! rectangular boundary conditions

     aaa=1.0_wp/cell(1)
     bbb=1.0_wp/cell(5)
     ccc=1.0_wp/cell(9)

!$OMP PARALLEL DO PRIVATE(i,xss,yss,zss)
     Do i=1,natms
        xss=aaa*xxx(i)
        yss=bbb*yyy(i)
        zss=ccc*zzz(i)

        xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
        zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

        xxx(i)=cell(1)*xss
        yyy(i)=cell(5)*yss
        zzz(i)=cell(9)*zss
     End Do
!$OMP END PARALLEL DO
  Else If (imcon == 3) Then

! parallelepiped boundary conditions

     Call invert(cell,rcell,det)

     Do i=1,natms
        xss=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        yss=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        zss=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

        xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
        zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

        xxx(i)=cell(1)*xss+cell(4)*yss+cell(7)*zss
        yyy(i)=cell(2)*xss+cell(5)*yss+cell(8)*zss
        zzz(i)=cell(3)*xss+cell(6)*yss+cell(9)*zss
     End Do

  Else If (imcon == 4) Then

! truncated octahedral boundary conditions

     If (.not.(Abs(cell(1)-cell(5)) < 1.0e-6_wp .and. Abs(cell(5)-cell(9)) < 1.0e-6_wp)) Call error(130)

     aaa=1.0_wp/cell(1)

     Do i=1,natms
        xss=aaa*xxx(i)
        yss=aaa*yyy(i)
        zss=aaa*zzz(i)

        xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
        zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

        xxx(i)=cell(1)*xss
        yyy(i)=cell(1)*yss
        zzz(i)=cell(1)*zss

        If ((Abs(xxx(i))+Abs(yyy(i))+Abs(zzz(i))) >= 0.75_wp*cell(1)) Then
           xxx(i)=xxx(i)-0.5_wp*Sign(cell(1),xxx(i))
           yyy(i)=yyy(i)-0.5_wp*Sign(cell(1),yyy(i))
           zzz(i)=zzz(i)-0.5_wp*Sign(cell(1),zzz(i))

           xss=aaa*xxx(i)
           yss=aaa*yyy(i)
           zss=aaa*zzz(i)

           xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
           yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
           zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

           xxx(i)=cell(1)*xss
           yyy(i)=cell(1)*yss
           zzz(i)=cell(1)*zss
        End If
     End Do

  Else If (imcon == 5) Then

! rhombic Dodecahedral boundary conditions

     If (.not.(Abs(cell(1)-cell(5)) < 1.0e-6_wp .and. Abs(cell(9)-cell(1)*rt2) < 1.0e-6_wp)) Call error(140)

     aaa=1.0_wp/cell(1)
     bbb=1.0_wp/cell(9)

     Do i=1,natms
        xss=aaa*xxx(i)
        yss=aaa*yyy(i)
        zss=bbb*zzz(i)

        xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
        zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

        xxx(i)=cell(1)*xss
        yyy(i)=cell(1)*yss
        zzz(i)=cell(9)*zss

        If ((Abs(xxx(i))+Abs(yyy(i))+Abs(rt2*zzz(i))) >= cell(1)) Then
           xxx(i)=xxx(i)-0.5_wp*Sign(cell(1),xxx(i))
           yyy(i)=yyy(i)-0.5_wp*Sign(cell(1),yyy(i))
           zzz(i)=zzz(i)-0.5_wp*Sign(cell(9),zzz(i))

           xss=aaa*xxx(i)
           yss=aaa*yyy(i)
           zss=bbb*zzz(i)

           xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
           yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
           zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss

           xxx(i)=cell(1)*xss
           yyy(i)=cell(1)*yss
           zzz(i)=cell(9)*zss
        End If
     End Do

  Else If (imcon == 6) Then

! x-y boundary conditions (SLAB)

     det=cell(1)*cell(5)-cell(2)*cell(4)

     If (Abs(det) < 1.0e-6_wp) Call error(120)

     det=1.0_wp/det

     rcell(1) =  det*cell(5)
     rcell(2) = -det*cell(2)
     rcell(4) = -det*cell(4)
     rcell(5) =  det*cell(1)

     Do i=1,natms
        xss=rcell(1)*xxx(i)+rcell(4)*yyy(i)
        yss=rcell(2)*xxx(i)+rcell(5)*yyy(i)

        xss=xss-Anint(xss) ; If (xss >= half_minus) xss=-xss
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss

        xxx(i)=cell(1)*xss+cell(4)*yss
        yyy(i)=cell(2)*xss+cell(5)*yss
     End Do

  Else If (imcon == 7) Then

! hexagonal prism boundary conditions

     If (Abs(cell(1)-rt3*cell(5)) > 1.0e-6_wp) Call error(135)

     aaa=cell(1)/(rt3*2.0_wp)
     bbb=cell(1)/rt3
     ccc=rt3/cell(1)
     ddd=1.0_wp/cell(9)

     Do i=1,natms
        zss=ddd*zzz(i)
        zss=zss-Anint(zss) ; If (zss >= half_minus) zss=-zss
        zzz(i)=cell(9)*zss

        yss=ccc*yyy(i)
        yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
        yyy(i)=bbb*yss

        If ((Abs(yyy(i))+Abs(rt3*xxx(i))) >= bbb) Then
           xxx(i)=xxx(i)-rt3*Sign(aaa,xxx(i))
           yyy(i)=yyy(i)-Sign(aaa,yyy(i))

           yss=ccc*yyy(i)
           yss=yss-Anint(yss) ; If (yss >= half_minus) yss=-yss
           yyy(i)=bbb*yss
        End If
     End Do

  End If

End Subroutine pbcshift

Subroutine jacobi(n,aaa,vvv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Diagonalisation of real square symmetric matrices by the Jacobi method:
! a sequence of Jacobi rotations
!
! Users must ensure the symmetry of the input matrix
!
! input parameters: n   - matrix dimension
!                   aaa - the matrix to be diagonalised
!                   vvv - the (diagonalised) eigenvector matrix
!
! Jacobi processes lower triangle only - strictly upper triangle
!                                        remains unchanged
!
! Variable rho sets absolute tolerance on convergence
! Variable test is a moving tolerance that diminishes on each pass
! until true convergence test<rho
! ZERO matrices are accepted and returned
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : zero_plus

  Integer,                                    Intent( In    ) :: n
  Real( Kind = wp ), Dimension ( 1:n , 1:n ), Intent( InOut ) :: aaa,vvv

  Real( Kind = wp ), Parameter :: rho = 1.0e-20_wp
  Logical           :: pass
  Integer           :: i,j,k ! ,l
  Real( Kind = wp ) :: scale,test,                      &
                       v_d_hor,v_d_off,v_d_ver,v_d_mid, &
                       omg,s_s,s,c_c,c,s_c,tmp

!  l=0 ! Iteration counter

! Rescale (lower triangle) matrix for optimal accuracy
! by the largest by magnitude diagonal element

  scale=0.0_wp
  Do i=1,n
     If (Abs(aaa(i,i)) > scale) scale=Abs(aaa(i,i))
  End Do
  If (scale <= zero_plus) Then
     vvv=aaa
     Return ! Accept & Return zero matrices
  Else
     Do i=1,n
        Do j=1,i
           aaa(j,i)=aaa(j,i)/scale
        End Do
     End Do
  End If

! Set initial value of moving tolerance
! Sum of all off-diagonal elements (strictly lower triangle)

  test=0.0_wp
  Do i=2,n
     Do j=1,i-1
        test=test+aaa(j,i)**2
     End Do
  End Do
  test=Sqrt(2.0_wp*test)

! Initialise eigenvectors

  vvv=0.0_wp
  Do i=1,n
     vvv(i,i)=1.0_wp
  End Do

! Accept & Return already diagonalised matrices
! (as well as zero matrices)

  If (test < rho) Return

! Recycle until absolute tolerance satisfied

  Do While (test > rho)
     test=test/Real(n,wp)
     If (test < rho) test=rho

! Jacobi diagonalisation

     pass=.true.

! Recycle until moving tolerance satisfied

     Do While (pass)
        pass=.false.

! Loop around the strictly lower triangle matrix

        Do i=2,n
           Do j=1,i-1
              If (Abs(aaa(j,i)) >= test) Then
!                 l=l+1
                 pass=.true.

                 v_d_hor=aaa(i,i)
                 v_d_ver=aaa(j,j)
                 v_d_off=aaa(j,i)
                 v_d_mid=0.5_wp*(v_d_ver-v_d_hor)
                 If (Abs(v_d_mid) < rho) Then
                    omg=-1.0_wp
                 Else
                    omg=-v_d_off/Sqrt(v_d_off**2+v_d_mid**2)
                    If (v_d_mid < 0.0_wp) omg=-omg
                 End If
                 s=omg/Sqrt(2.0_wp*(1.0_wp+Sqrt(1.0_wp-omg**2)))
                 s_s=s*s ; c_c=1.0_wp-s_s ; c=Sqrt(c_c) ; s_c=s*c

                 Do k=1,n
                    If      (k <= j) Then
                       tmp     =aaa(k,j)*c-aaa(k,i)*s
                       aaa(k,i)=aaa(k,j)*s+aaa(k,i)*c
                       aaa(k,j)=tmp
                    Else If (k >  i) Then
                       tmp     =aaa(j,k)*c-aaa(i,k)*s
                       aaa(i,k)=aaa(j,k)*s+aaa(i,k)*c
                       aaa(j,k)=tmp
                    Else
                       tmp     =aaa(j,k)*c-aaa(k,i)*s
                       aaa(k,i)=aaa(j,k)*s+aaa(k,i)*c
                       aaa(j,k)=tmp
                    End If

                    tmp     =vvv(k,j)*c-vvv(k,i)*s
                    vvv(k,i)=vvv(k,j)*s+vvv(k,i)*c
                    vvv(k,j)=tmp
                 End Do

                 aaa(i,i)=v_d_hor*c_c+v_d_ver*s_s+2.0_wp*v_d_off*s_c
                 aaa(j,j)=v_d_hor*s_s+v_d_ver*c_c-2.0_wp*v_d_off*s_c
                 aaa(j,i)=(v_d_ver-v_d_hor)*s_c+v_d_off*(c_c-s_s)
              End If
           End Do
        End Do
     End Do
  End Do

! Rescale back the lower triangle matrix

  Do i=1,n
     Do j=1,i
        aaa(j,i)=aaa(j,i)*scale
     End Do
  End Do

End Subroutine jacobi

Subroutine mat_mul(aaa,bbb,ccc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! matrix multiply routine: A*B=C (note order!)
!
! Note: A, B and C are 3x3 matrices in linear arrays as used in dl_poly
!
! copyright - daresbury laboratory
! author    - w.smith april 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Real( Kind = wp ), Intent( In    ) :: aaa(1:9),bbb(1:9)
  Real( Kind = wp ), Intent(   Out ) :: ccc(1:9)

  ccc(1)=aaa(1)*bbb(1)+aaa(4)*bbb(2)+aaa(7)*bbb(3)
  ccc(2)=aaa(2)*bbb(1)+aaa(5)*bbb(2)+aaa(8)*bbb(3)
  ccc(3)=aaa(3)*bbb(1)+aaa(6)*bbb(2)+aaa(9)*bbb(3)

  ccc(4)=aaa(1)*bbb(4)+aaa(4)*bbb(5)+aaa(7)*bbb(6)
  ccc(5)=aaa(2)*bbb(4)+aaa(5)*bbb(5)+aaa(8)*bbb(6)
  ccc(6)=aaa(3)*bbb(4)+aaa(6)*bbb(5)+aaa(9)*bbb(6)

  ccc(7)=aaa(1)*bbb(7)+aaa(4)*bbb(8)+aaa(7)*bbb(9)
  ccc(8)=aaa(2)*bbb(7)+aaa(5)*bbb(8)+aaa(8)*bbb(9)
  ccc(9)=aaa(3)*bbb(7)+aaa(6)*bbb(8)+aaa(9)*bbb(9)

End Subroutine mat_mul
