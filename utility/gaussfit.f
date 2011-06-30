      program gaussfit

c***********************************************************************
c
c     program to fit a lennard jones potential with a sum of gaussians
c
c     author k. singer (modified by w. smith)
c
c***********************************************************************

      implicit none

      integer i,k,l,m
      real(8) varm,r,rr,f1,f2,dr1,dr2,a11,dr3,a12,a22,a21,a31,a32,a33
      real(8) a13,a23,b1,b2,b3,det,absd,var,ca1,ca2,ca3,da1,da2,da3
      real(8) c1,c2,c3
      real(8) f(200),g1(200),g2(200),g3(200),gg(200),dd(3),dx(3),d(3)

      varm=1.0d+25
      write(6,"(1x,'gaussian potential fitting program',//)")

c     set up array for standard potential

      r=0.91d0
      do i=1,163

        r=r+0.01d0
        rr=r*r
        f1=1.d0/rr
        f2=f1*f1*f1
        f1=f2*f2
        f(i)=4.d0*(f1-f2)

      enddo

c     set initial guesses for exponents

      dx(1)=0.443d0
      dx(2)=1.544d0
      dx(3)=8.502d0
      dd(1)=3.889d-3
      dd(2)=3.889d-3
      dd(3)=7.779d-3

c     adjustment of first exponent

      do k=0,5

        d(1)=dx(1)+dble(k-3)*dd(1)

c     adjustment of second exponent

        do l=0,5

          d(2)=dx(2)+dble(l-3)*dd(2)

c     adjustment of third exponent

          do m=0,5

            d(3)=dx(3)+float(m-3)*dd(3)

c     store current exponents

            dr1=d(1)
            dr2=d(2)
            dr3=d(3)

c     initialise least squares parameters

            a11=0.d0
            a12=0.d0
            a22=0.d0
            a21=0.d0
            a31=0.d0
            a32=0.d0
            a33=0.d0
            a13=0.d0
            a23=0.d0
            b1=0.d0
            b2=0.d0
            b3=0.d0

c     calculate least squares parameters

            r=0.91d0
            do i=1,163

              r=r+0.01d0
              rr=r*r
              g1(i)=exp(-rr*dr1)
              g2(i)=exp(-rr*dr2)
              g3(i)=exp(-rr*dr3)
              a11=a11+g1(i)**2
              a12=a12+g1(i)*g2(i)
              a22=a22+g2(i)**2
              a13=a13+g1(i)*g3(i)
              a23=a23+g2(i)*g3(i)
              a33=a33+g3(i)**2
              b1=b1+g1(i)*f(i)
              b2=b2+g2(i)*f(i)
              b3=b3+g3(i)*f(i)

            enddo

            det=a11*(a22*a33-a23*a23)
     x        +a12*(a23*a13-a12*a33)
     x        +a13*(a12*a23-a22*a13)
            absd=abs(det)
            if(absd.le.2.0d-20)then
              write(6,"(1x,'error - determinant too small')")
            else

c     calculate coefficients

              c1=(b1*(a22*a33-a23*a23)+b2*(a13*a23-a12*a33)+
     x          b3*(a12*a23-a13*a22))/det
              c2=(b1*(a23*a13-a12*a33)+b2*(a11*a33-a13*a13)+
     x          b3*(a12*a13-a11*a23))/det
              c3=(b1*(a12*a23-a22*a13)+b2*(a12*a13-a11*a23)+
     x          b3*(a11*a22-a12*a12))/det

            endif

c     construct approximating potential array

            var=0.d0
            r=0.91d0
            do i=1,163

              r=r+0.01d0
              gg(i)=c1*g1(i)+c2*g2(i)+c3*g3(i)
              var=var+r*(f(i)-gg(i))**2

            enddo

            if(var.le.varm)then

c     update potential parameters

              write(6,'(1x,3i5)')k,l,m
              varm=var
              ca1=c1
              ca2=c2
              ca3=c3
              da1=d(1)
              da2=d(2)
              da3=d(3)

c     print out current values

              write(6,"(1x,'var = ',e15.6,/,1x,'c(l)= ',3e15.6,/,
     x          1x,'d(l)= ',3e15.6)")var,ca1,ca2,ca3,da1,da2,da3

            endif

          enddo

          dx(3)=da3

        enddo

        dx(2)=da2

      enddo

      dx(1)=da1

c     print out best parameters

      write(6,"(1x,'best coefficients',5x,3e15.6,/,
     x  1x,'best exponents   ',5x,3e15.6,/,/,
     x  1x,10x,'rad',11x,'f(i)',10x,'gg(i)',/)")
     x  ca1,ca2,ca3,da1,da2,da3

c     print out potential arrays for comparison

      r=0.91d0
      do i=1,163

        r=r+0.01d0
        write(6,'(3e15.6)')r,f(i),gg(i)

      enddo

      end

