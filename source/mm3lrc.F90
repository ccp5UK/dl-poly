Module mm3lrc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to with functions to evaluate numerically the LRC of AMOEBA
! 14-7 buffered vdw potential
!
! copyright - daresbury laboratory
! author    - a.m.elena march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp

  Implicit None

  Private

  Public :: intRadMM3
  Public :: intRaddMM3

Contains

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

End Module mm3lrc

