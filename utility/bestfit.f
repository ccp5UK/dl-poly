      program bestfit

c***********************************************************************
c
c     Least squares fitting of arbitary functions to data
c     uses Gram-Schmidt Orthogonalisation procedure.
c
c     copyright daresbury laboratory 1993
c     author t. forester may 1993.
c
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c
c***********************************************************************

      implicit real*8(a-h,o-z)


c     maximum number of functions to fit

      Parameter (M1 = 20)

c     maximum number of data points (> M)

      Parameter (N1 = 1000)

      Dimension A(m1,n1), B(m1,m1),C(m1,n1),x(n1),y(n1),d(m1)


      write(*,*) 'Gram Schmidt Orthogonalisation with a maximum of '
     x     ,m1,' functions and ',n1, ' data points'

      write(*,'(//,a30)')'number of functions to use ?'
      read(*,*) m
      if(m.gt.m1) call error(1)

      write(*,*)'enter points x, y {type" -999 -999" to terminate list}'

      n = 0
   10 read(*,*) x1,y1
      if (x1.eq.-999.d0.and.y1.eq.-999.d0) goto 20
      if(n.ge.n1) call error(2)
      n = n+1
      x(n) = x1
      y(n) = y1
      goto 10
      
   20 continue

      if(n.lt.m) call error(4)
c
c     create a matrix - change function definitions here - only the first m
c     functions are used.

      do i = 1,n
c
c     define functions as power series
c     but in general may be any function you like

         if(m.ge.1) a(1,i) = 1
         if(m.ge.2) a(2,i) = x(i)
         if(m.ge.3) a(3,i) = x(i)*a(2,i)
         if(m.ge.4) a(4,i) = x(i)*a(3,i)
         if(m.ge.5) a(5,i) = x(i)*a(4,i)
         if(m.ge.6) a(6,i) = x(i)*a(5,i)
         if(m.ge.7) a(7,i) = x(i)*a(6,i)
         if(m.ge.8) a(8,i) = x(i)*a(7,i)
         if(m.ge.9) a(9,i) = x(i)*a(8,i)
         if(m.ge.10) a(10,i) = x(i)*a(9,i)
         if(m.ge.11) a(11,i) = x(i)*a(10,i)
         if(m.ge.12) a(12,i) = x(i)*a(11,i)
         if(m.ge.13) a(13,i) = x(i)*a(12,i)
         if(m.ge.14) a(14,i) = x(i)*a(13,i)
         if(m.ge.15) a(15,i) = x(i)*a(14,i)
         if(m.ge.16) a(16,i) = x(i)*a(15,i)
         if(m.ge.17) a(17,i) = x(i)*a(16,i)
         if(m.ge.18) a(18,i) = x(i)*a(17,i)
         if(m.ge.19) a(19,i) = x(i)*a(18,i)
         if(m.ge.20) a(20,i) = x(i)*a(19,i)

      enddo

c
c     construct b = a.(a transpose)  (m x m)

      do i = 1,m

         do j = i,m

            b(i,j) = 0.d0

            do k = 1,n

               b(i,j) = b(i,j) + a(i,k)*a(j,k)

            enddo

            b(j,i) = b(i,j)

         enddo

      enddo
      
c
c     find inverse of b

      call invert(m,m1,b,b)

c
c     construct (b^-1)a  (m x n)

      do i = 1,m

         do j = 1,n

            c(i,j) = 0.0d0

            do k = 1,m

               c(i,j) = c(i,j) + b(i,k)*a(k,j)

            enddo

         enddo

      enddo

c
c     find coefficients 

      do i = 1,m

         d(i) = 0.d0

         do j = 1,n

            d(i) = d(i) + c(i,j)*y(j)

         enddo

         write(*,198) i,d(i)
  198    format (' coefficient',i5,4x,1p,e16.8)

      enddo


c
c     calculate best fit of data points
      write(*,*)
      write(*,'(3(5x,a5,5x))') 'x','y','fit'
      do i = 1,n

         ax = 0.d0

         do j = 1,m

            ax = ax + a(j,i)*d(j)

         enddo

         write(*,95) x(i), y(i), ax

   95    format(1p,3e16.8)

      enddo

      end
      subroutine error(i)

    1 format(/,' error - too many functions specified')
    2 format(/,' error - too many data points specified')
    3 format(/,' error - matrix too big for "invert"')
    4 format(/,' error - too few data points specified')

      if(i.eq.1) write(*,1)
      if(i.eq.2) write(*,2)
      if(i.eq.3) write(*,3)
      if(i.eq.4) write(*,4)

      call exit()
      return
      end

      subroutine invert(nnn,nmax,aaa,bbb)
c     
c***********************************************************************
c     
c     routine to invert a real symmetric matrix (reference:
c     computing methods v.ii, i.s. berezin and n.p.zhidkov,
c     pergamon press 1965). note that the matrices are
c     packed in the minimum storage mode, (i.e. a(k)=a(i,j),
c     where i>j and k=(i*(i-1))/2+j ).
c     the matrices aaa and bbb may be equivalenced though
c     this will destroy the contents of the original
c     array.
c     
c     general version for all real symmetric matrices
c     
c     copyright daresbury laboratory 1993
c     author w.smith (added to dl_poly july 1993)
c     
c***********************************************************************
c     
      implicit real*8(a-h,o-z)
c
c     maximum matrix (mbig*mbig) to invert

      parameter (mbig =100, bmax = mbig*(mbig+1)/2)
      dimension aaa(nmax*nmax),bbb(nmax*nmax),bbi(bmax)

      ind(i,j)=(max0(i,j)*(max0(i,j)-1))/2+min0(i,j)

      if(nnn.gt.mbig) call error(3)
c
c     pack matrix aaa into a triangular form

      do i = 1,nnn
         do j = 1,i
            k = ind(i,j)
            k1 = (i-1)*nmax + j
            bbi(k) = aaa(k1)
         enddo
      enddo
      do i = 1,nnn
         do j = 1,i
            k = ind(i,j)
            aaa(k) = bbi(k)
         enddo
      enddo

         
c     
c     factorize matrix aaa into upper and lower forms
      do l=1,nnn
         do m=1,l
            k=ind(l,m)
            bbi(k)=0.d0
            bbb(k)=aaa(k)
            if(m.ne.1)then
               lm=m-1
               do n=1,lm
                  ln=ind(l,n)
                  mn=ind(m,n)
                  temv  =bbb(k)-bbb(ln)*bbb(mn)+bbi(ln)*bbi(mn)
                  bbi(k)=bbi(k)-bbi(ln)*bbb(mn)-bbb(ln)*bbi(mn)
                  bbb(k)=temv
               enddo
            endif
            if(l.eq.m)then
               if(bbb(k).lt.0.d0)then
                  bbi(k)=sqrt(abs(bbb(k)))
                  bbb(k)=0.d0
               else
                  bbb(k)=sqrt(abs(bbb(k)))
                  bbi(k)=0.d0
               endif
            else
               mm=ind(m,m)
               den=1.d0/(bbb(mm)**2+bbi(mm)**2)
               temv  =den*(bbb(k)*bbb(mm)+bbi(k)*bbi(mm))
               bbi(k)=den*(bbi(k)*bbb(mm)-bbb(k)*bbi(mm))
               bbb(k)=temv
            endif
         enddo
      enddo
c     
c     invert lower triangular matrix
      do l=1,nnn
         n=ind(l,l)
         den=1.d0/(bbb(n)*bbb(n)+bbi(n)*bbi(n))
         do m=1,l
            k=ind(l,m)
            if(l.eq.m)then
               bbb(k)= den*bbb(k)
               bbi(k)=-den*bbi(k)
            else
               lm=l-1
               suma=0.d0
               sumb=0.d0
               do j=m,lm
                  lj=ind(l,j)
                  mj=ind(m,j)
                  suma=suma-bbb(lj)*bbb(mj)+bbi(lj)*bbi(mj)
                  sumb=sumb-bbb(lj)*bbi(mj)-bbi(lj)*bbb(mj)
               enddo
               temv  =den*(suma*bbb(n)+sumb*bbi(n))
               bbi(k)=den*(sumb*bbb(n)-suma*bbi(n))
               bbb(k)=temv
            endif
         enddo
      enddo
c     
c     form product of upper and lower inverse triangular matrices
      do l=1,nnn
         do m=1,l
            sum=0.d0
            k=ind(l,m)
            do j=l,nnn
               lj=ind(l,j)
               mj=ind(m,j)
               sum=sum+bbb(lj)*bbb(mj)-bbi(lj)*bbi(mj)
            enddo
            bbb(k)=sum
         enddo
      enddo

c
c     unpack inverse matrix bbb

      do i = 1,nnn
         do j = 1,i
            k=ind(i,j)
            bbi(k) = bbb(k)
         enddo
      enddo

      do i = 1,nnn
         do j = 1,nnn
            k1 = (i-1)*nmax + j
            k = ind(i,j)
            bbb(k1) = bbi(k)
         enddo
      enddo

      return
      end

