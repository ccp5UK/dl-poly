Module zbl_pots
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring constants and ZBL related routines
  !
  ! copyright - daresbury laboratory
  ! author    - a.m.elena october 2017
  ! contrib   - a.m.elena december 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use Kinds, Only : wp
  Implicit None
  Private 

  Real(wp), Parameter, Dimension(4) :: b=[0.18175_wp,0.50986_wp,0.28022_wp,0.02817_wp], &
    c=[3.1998_wp,0.94229_wp,0.40290_wp,0.20162_wp]
  Real(wp), Parameter, Public :: ab=0.52917721067_wp 
  Public :: zbl
  Public :: zbls
  Public :: zblb
  Public :: intRadZBL
  Public :: intdRadZBL

Contains

  Pure Subroutine zbl(r,k,ia,phi,dphi)
    Real(wp), intent( In    )  :: r,k,ia
    Real(wp), intent(   Out ) :: phi,dphi

    Integer :: i
    Real(wp) :: x,t1,t2,ir

    phi = 0.0_wp
    dphi = 0.0_wp
    x = r*ia
    ir = 1.0_wp/r
    do i = 1,4
      t1 = b(i)*exp(-x*c(i))
      phi = phi+t1
      dphi = dphi-c(i)*t1
    end do
    phi = k*phi*ir
    ! -r∂U/δr
    dphi = phi - ia*k*dphi

  End Subroutine 

  Pure Subroutine fm(r,rm,ic,f,df)
    Real(wp), Intent( In    )  :: r,rm,ic
    Real(wp), Intent(   Out ) :: f,df

    Real(wp) :: t


    If (r<rm) Then
      t = Exp(-(rm-r)*ic)*0.5_wp
      f = 1.0_wp-t
      df = r*ic*t
    Else
      t = Exp(-(r-rm)*ic)*0.5_wp
      f = t
      df = r*ic*t
    End If
  End Subroutine fm

  Pure Subroutine morse(r,d,k,r0,m,dm)
    Real(wp), Intent( In    )  :: r,d,k,r0
    Real(wp), Intent(   Out ) ::m,dm

    Real(wp) :: t

    t = Exp(-k*(r-r0))

    m = d*((1.0_wp-t)**2-1.0_wp)
    dm=-2.0_wp*r*d*k*(1.0_wp-t)*t

  End Subroutine morse

  Pure Subroutine buckingham(r,A,r0,C,b,db)
    Real(wp), Intent( In    )  :: r,A,r0,C
    Real(wp), Intent(   Out ) ::b,db

    Real(wp) :: t1,t2,ir0

    ir0=1.0_wp/r0
    t1 = A*Exp(-r*ir0)
    t2 = -C*(1.0_wp/r)**6

    b = t1+t2
    db=r*ir0*t1+6.0_wp*t2
  End Subroutine buckingham

  Pure Subroutine zbls(r,kk,ia,rm,ic,d,k,r0,V,dV)
    Real(wp), Intent( In    )  :: r,kk,ia,ic,rm,d,k,r0
    Real(wp), Intent(   Out ) :: V,dV

    Real(wp) :: z,dz,f,df,m,dm

    call zbl(r,kk,ia,z,dz)
    call fm(r,rm,ic,f,df)
    call morse(r,d,k,r0,m,dm)
    V = f*z+(1.0_wp-f)*m
    dV = f*dz+df*z+(1.0_wp-f)*dm-df*m
  End Subroutine zbls

  Pure Subroutine zblb(r,kk,ia,rm,ic,A,r0,C,V,dV)
    Real(wp), Intent( In    )  :: r,kk,ia,ic,rm,A,C,r0
    Real(wp), Intent(   Out ) :: V,dV

    Real(wp) :: z,dz,f,df,b,db

    call zbl(r,kk,ia,z,dz)
    call fm(r,rm,ic,f,df)
    call buckingham(r,A,r0,C,b,db)
    V = f*z+(1.0_wp-f)*b
    dV = f*dz+df*z+(1.0_wp-f)*db-df*b

  End Subroutine zblb

  Pure Real(wp) Function intRadZBL(kk,a,rw,prec)
    Real(wp), Intent( In    ) :: kk,a,rw,prec

    Real(wp) :: is,ie,h,s,x,sold,f1,f0,f2,df0

    Integer           :: i,n,j

    n=10000
    is=rw
    ie=2*is
    h=(ie-is)/Real(n,wp)
    s=0.0_wp
    sold=Huge(1.0_wp)
    j=1

    Do While (Abs(s-sold)*h/3.0_wp > prec)
       sold = s

       Do i = (j-1)*n, j*n, 2

          x = is + i*h
          Call zbl(x,kk,a,f0,df0)
          Call zbl(x+h,kk,a,f1,df0)
          Call zbl(x+2.0_wp*h,kk,a,f2,df0)
          s = s  +             x*x * f0 + &
                   4.0_wp*(x+h)**2 * f1 + &
                   (x+2.0_wp*h)**2 * f2
       End Do

       j=j+1
    End Do
    intRadZBL=s*h/3.0_wp   
  End Function intRadZBL

  Pure Real(wp) Function intdRadZBL(kk,a,rw,prec)
    Real(wp), Intent( In    ) :: kk,a,rw,prec

    Real(wp) :: is,ie,h,s,x,sold,df1,df0,df2,f0

    Integer           :: i,n,j

    n=10000
    is=rw
    ie=2*is
    h=(ie-is)/Real(n,wp)
    s=0.0_wp
    sold=Huge(1.0_wp)
    j=1

    Do While (Abs(s-sold)*h/3.0_wp > prec)
       sold = s

       Do i = (j-1)*n, j*n, 2

          x = is + i*h
          Call zbl(x,kk,a,f0,df0)
          Call zbl(x+h,kk,a,f0,df1)
          Call zbl(x+2.0_wp*h,kk,a,f0,df1)
          s = s  +             x*x * df0 + &
                   4.0_wp*(x+h)**2 * df1 + &
                   (x+2.0_wp*h)**2 * df2
       End Do

       j=j+1
    End Do
    intdRadZBL=s*h/3.0_wp
  End Function intdRadZBL

End Module zbl_pots
