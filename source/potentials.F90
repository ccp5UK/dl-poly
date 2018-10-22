Module potentials
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring potentials and long range correctins for them
  !
  ! copyright - daresbury laboratory
  ! author    - a.m.elena october 2017 (zbl)
  ! contrib   - a.m.elena december 2017 (zbl extended)
  ! contrib   - a.m.elena may 2018 (mdf)
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
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
  Public :: mlj126
  Public :: mbuck
  Public :: mlj
  Public :: intdRadMDF
  Public :: intRadMDF
  Public :: intRadMM3
  Public :: intRaddMM3

Contains

  Pure Subroutine zbl(r,k,ia,phi,dphi)
    Real(wp), intent( In    )  :: r,k,ia
    Real(wp), intent(   Out ) :: phi,dphi

    Integer :: i
    Real(wp) :: x,t1,ir

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
        Call zbl(x+2.0_wp*h,kk,a,f0,df2)
        s = s  +             x*x * df0 + &
          4.0_wp*(x+h)**2 * df1 + &
          (x+2.0_wp*h)**2 * df2
      End Do

      j=j+1
    End Do
    intdRadZBL=s*h/3.0_wp
  End Function intdRadZBL


  Pure Subroutine LJ(r,eps,sig,e,v)
    Real(wp), Intent( In    ) :: r,eps,sig
    Real(wp), Intent(   Out ) :: e,v

    Real(wp) :: sor6

    sor6=(sig/r)**6
    e = 4.0_wp*eps*sor6*(sor6-1.0_wp)
    v = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)

  End Subroutine LJ

  Pure Subroutine LJ126(r,A,B,e,v)
    Real(wp), Intent( In    ) :: r,A,B
    Real(wp), Intent(   Out ) :: e,v
    Real(wp) :: sor6

    sor6=(1.0_wp/r)**6
    e = (A*sor6-B)*sor6
    v = 6.0_wp*sor6*(2.0_wp*a*sor6-b)

  End Subroutine LJ126

  Pure Subroutine MDF(r,ri,rc,e,v)
    Real(wp), Intent( In    ) :: r,ri,rc
    Real(wp), Intent(   Out ) :: e,v

    Real(wp) :: rci

    If (r<ri) Then
      e = 1.0_wp
      v = 0.0_wp
    Else IF(r>rc) Then
      e = 0.0_wp
      v = 0.0_wp
    Else
      rci=(rc-ri)**5
      e = (rc-r)**3*(10*ri**2-5*rc*ri-15*r*ri+rc**2+3*r*rc+6*r**2)/rci
      v = 30.0_wp*r*(r-rc)**2*(r-ri)**2/rci
    End If

  End Subroutine MDF

  Pure Subroutine mlj(r,eps,sig,ri,rc,e,v)
    Real(wp), Intent( In    ) :: r,eps,sig,ri,rc
    Real(wp), Intent(   Out ) :: e,v

    Real(wp) :: el,vl,em,vm

    Call LJ(r,eps,sig,el,vl)
    Call MDF(r,ri,rc,em,vm)
    e = el*em
    v = vl*em+vm*el

  End Subroutine mlj

  Pure Subroutine mbuck(r,A,r0,c,ri,rc,e,v)
    Real(wp), Intent( In    ) :: r,A,r0,c,ri,rc
    Real(wp), Intent(   Out ) :: e,v

    Real(wp) :: eb,vb,em,vm

    Call buckingham(r,A,r0,c,eb,vb)
    Call MDF(r,ri,rc,em,vm)
    e = eb*em
    v = vb*em+vm*eb

  End Subroutine mbuck

  Pure Subroutine mlj126(r,A,B,ri,rc,e,v)
    Real(wp), Intent( In    ) :: r,A,B,ri,rc
    Real(wp), Intent(   Out ) :: e,v

    Real(wp) :: el,vl,em,vm

    Call LJ126(r,A,B,el,vl)
    Call MDF(r,ri,rc,em,vm)
    e = el*em
    v = vl*em+vm*el
  End Subroutine mlj126

  Pure Real(wp) Function intRadMDF(pot,a,b,c,ri,rw,prec)
    Real(wp), Intent( In    )             :: a,b,c,ri,rw,prec
    Character( Len = * ), Intent( In    ) :: pot

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
        ! this is stupid by remember we do it only once
        If (pot == 'm126' ) Then
          Call mlj126(x,a,b,ri,rw,f0,df0)
          Call mlj126(x+h,a,b,ri,rw,f1,df0)
          Call mlj126(x+2.0_wp*h,a,b,ri,rw,f2,df0)
        Else If (pot == 'mbuc' ) Then
          Call mbuck(x,a,b,c,ri,rw,f0,df0)
          Call mbuck(x+h,a,b,c,ri,rw,f1,df0)
          Call mbuck(x+2.0_wp*h,a,b,c,ri,rw,f2,df0)
        Else If (pot == 'mlj' ) Then
          Call mlj(x,a,b,ri,rw,f0,df0)
          Call mlj(x+h,a,b,ri,rw,f1,df0)
          Call mlj(x+2.0_wp*h,a,b,ri,rw,f2,df0)
        End If

        s = s  +             x*x * f0 + &
          4.0_wp*(x+h)**2 * f1 + &
          (x+2.0_wp*h)**2 * f2
      End Do

      j=j+1
    End Do
    intRadMDF=s*h/3.0_wp
  End Function intRadMDF

  Pure Real(wp) Function intdRadMDF(pot,a,b,c,ri,rw,prec)
    Real(wp), Intent( In    )             :: a,b,c,ri,rw,prec
    Character( Len = * ), Intent( In    ) :: pot

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
        If (pot == 'm126' ) Then
          Call mlj126(x,a,b,ri,rw,f0,df0)
          Call mlj126(x+h,a,b,ri,rw,f0,df1)
          Call mlj126(x+2.0_wp*h,a,b,ri,rw,f0,df2)
        Else If (pot == 'mbuc' ) Then
          Call mbuck(x,a,b,c,ri,rw,f0,df0)
          Call mbuck(x+h,a,b,c,ri,rw,f0,df1)
          Call mbuck(x+2.0_wp*h,a,b,c,ri,rw,f0,df2)
        Else If (pot == 'mlj' ) Then
          Call mlj(x,a,b,ri,rw,f0,df0)
          Call mlj(x+h,a,b,ri,rw,f0,df1)
          Call mlj(x+2.0_wp*h,a,b,ri,rw,f0,df2)
        End If
        s = s  +             x*x * df0 + &
          4.0_wp*(x+h)**2 * df1 + &
          (x+2.0_wp*h)**2 * df2
      End Do

      j=j+1
    End Do
    intdRadMDF=s*h/3.0_wp
  End Function intdRadMDF

  Function mm3(x,A,B,eps)

    Real( Kind = wp ) :: mm3
    Real( Kind = wp ), Intent( In    ) :: x,A,B,eps

    Real( Kind = wp ) :: ia,ib

    ia  = 1.0_wp/(A+x   )
    ib  = 1.0_wp/(B+x**7)

    mm3 = eps*((1.0_wp+A)*ia)**7 * ((1.0_wp+B)*ib - 2.0_wp)

  End Function mm3


  Function intRadMM3(r0,A,B,eps,rw,prec)

    Real( Kind = wp ) :: intRadMM3
    Real( Kind = wp ), Intent( In    ) :: r0,A,B,eps,rw,prec

    Real( Kind = wp ) :: is,ie,h,s,x,sold

    Integer           :: i,n,j

    n=10000
    is=rw/r0
    ie=2*is
    h=(ie-is)/Real(n,wp)
    s=0.0_wp
    sold=Huge(1.0_wp)
    j=1

    Do While (Abs(s-sold)*h/3.0_wp > prec)
      sold = s

      Do i = (j-1)*n, j*n, 2
        x = is + i*h
        s = s  +             x*x * mm3(x,    A,B,eps) + &
          4.0_wp*(x+h)**2 * mm3(x+h,  A,B,eps) + &
          (x+2.0_wp*h)**2 * mm3(x+2*h,A,B,eps)
      End Do

      j=j+1
    End Do

    intRadMM3=s*h*r0**3/3.0_wp

  End Function intRadMM3

  Function dmm3(x,r0,A,B,eps)

    Real( Kind = wp ) :: dmm3
    Real( kind = wp ), Intent( In    ) :: x,r0,A,B,eps

    Real( kind = wp ) :: ia,ib

    ia   = 1.0_wp/(A+x   )
    ib   = 1.0_wp/(B+x**7)

    dmm3 = -7.0_wp*(mm3(x,A,B,eps)*ia - ((1.0_wp+A)*ia)**7 * ((1.0_wp+B)*ib**2) * x**6)/r0

  End Function dmm3

  Function intRaddMM3(r0,A,B,eps,rw,prec)

    Real( Kind = wp ) :: intRaddMM3
    Real( kind = wp ), Intent( In    ) :: r0,A,B,eps,rw,prec

    Real( kind = wp ) :: is,ie,h,s,x,sold
    Integer           :: i,n,j

    n=10000
    is=rw/r0
    ie=2*is
    h=(ie-is)/Real(n,wp)
    s=0.0_wp
    sold=Huge(1.0_wp)
    j=1

    Do While (Abs(s-sold)*h/3.0_wp > prec)
      sold=s

      Do i=(j-1)*n,j*n,2
        x = is+i*h

        s = s + x**3            * dmm3(x,    r0,A,B,eps) + &
          4.0_wp*(x+h)**3 * dmm3(x+h,  r0,A,B,eps) + &
          (x+2.0_wp*h)**3 * dmm3(x+2*h,r0,A,B,eps)
      End Do

      j=j+1
    End Do

    intRaddMM3=s*h*r0**4/3.0_wp

  End Function intRaddMM3

End Module potentials
