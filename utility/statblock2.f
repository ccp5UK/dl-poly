      program statblock2
c
c**********************************************************************
c
c     blocking method for estimating standard error of mean
c
c     (reference - tildesley and allen, p191  (1987)
c
c     copyright daresbury laboratory 1996
c     author - t. forester sept 1996
c
c**********************************************************************
c
      parameter (ndeg=25,mxrows=10000,mxcols=20)
      implicit real*8 (a-h,o-z)
      dimension list(mxcols),record(mxcols),d(mxrows,mxcols)
      dimension aver(mxcols),var(ndeg,mxcols)
      dimension x(ndeg),y(ndeg)
      dimension err(mxcols),seff(mxcols)
      character*40 statfile

      write(*,*)'error estimation using blocking method'
      write(*,*)'maximum number of time points allowed =',mxrows
      write(*,*)'name of statistical data file (40 chars max)?'
      read(*,'(a40)')statfile
      write(*,*)'enter number of columns in data file (max 20)'
      read(*,*)mcols
      mcols=min(mcols,mxcols)
      write(*,*)'enter number of columns to be processed (max 20)'
      read(*,*)ncols
      ncols=min(ncols,mxcols)
      write(*,*)'identify columns to be processed'
      read(*,*)(list(j),j=1,ncols)
c
c     open the statistical data file

      open(7,file=statfile)

c
c     read in statistical data
      do i=1,mxrows
        read(7,*,end=100)(record(j),j=1,mcols)
        do k=1,ncols
          d(i,k)=record(list(k))
        enddo
      enddo
  100 m=min(i-1,mxrows)
      write(*,*)'number of data points read',m
      if(m.le.2) then
        write(*,*) 'too few data points to process'
        call exit()
      endif

c
c     calculate errors using blocking method

      niter = min(m,ndeg)
      write(*,*) 'number of blocks used ',niter

c
c     calculate averages
      write(*,*)'calculated averages'
      do l=1,ncols
        aver(l)=0.d0
        do i=1,m
          aver(l)=aver(l)+d(i,l)
        enddo
        aver(l)=aver(l)/dble(m)
      enddo
      write(*,'(3x,1p,20e12.4)')(aver(l),l=1,ncols)

      do l=1,ncols

        var(1,l) =0.d0

        do  i=1,m
          var(1,l) = var(1,l) + (d(i,l)- aver(l))**2
        enddo

        var(1,l)= var(1,l)/dble(m)

      enddo

      do iter = 2,niter
c
c     block length

        tb = dble(m)/dble(iter)

        do l=1,ncols

          var(iter,l)=0.d0

          do  j=1,iter

            dd = 0.d0
            j1 = 1 + nint(dble(j-1)*tb)
            j2 = nint(dble(j)*tb)
            nb = j2 -j1 + 1

            do i = j1,j2

              dd = dd + d(i,l)

            enddo
c
c     dd is mean of block value

            dd = dd/dble(nb)
            var(iter,l)=var(iter,l)+((dd-aver(l))**2)

          enddo
          var(iter,l)=var(iter,l)/dble(iter)

        enddo

      enddo
c
c     linear regression of s = nb* var()^2/ var(1)^2
c     to find statistical uncertainty

      n0 = 2
      n1 = n0 -1

      do l = 1,ncols

        do j = n0,niter

          x(j-n1) = dble(j)/dble(m)
          if(var(1,l).gt.1d-8) then
            y(j-n1) = (var(j,l)/var(1,l))/x(j-n1)
          else
            y(j-n1) = 0.d0
          endif
        enddo
        nn = niter-n1
        call bestfit(nn,x,y,a1)
        a1 = max(a1,1.d0)
        err(l) = sqrt(var(1,l)*a1/dble(m))
        seff(l)= a1
      enddo
c$$$      write(*,'(/,a)')'# of blocks      r.m.s. fluctuations'
c$$$      do j = 1,niter
c$$$        write(*,'(i3,1p,20e12.4)')j,(sqrt(var(j,l)),l=1,ncols)
c$$$      enddo

      write(*,'(/,a)')' calculated uncertainty in averages'
      write(*,'(3x,1p,20e12.4)')(err(j),j=1,ncols)

      write(*,'(/,a)')' calculated statistical inefficiencies'
      write(*,'(3x,1p,20e12.4)')(seff(j),j=1,ncols)

      end

      subroutine bestfit(n,x,y,a1)

c***********************************************************************
c
c     least squares fitting of arbitary functions to data
c     uses gram-schmidt orthogonalisation procedure.
c
c     copyright daresbury laboratory 1993
c     author t. forester may 1993.
c
c     $author: wl $
c     $date: 1996/02/15 14:33:28 $
c     $revision: 1.1.1.1 $
c     $state: exp $
c
c***********************************************************************

      implicit real*8(a-h,o-z)


c     maximum number of functions to fit

      parameter (m1 = 5)

c     maximum number of data points (> m)

      parameter (n1 = 25)

      dimension a(m1,n1), b(m1,m1),c(m1,n1),x(n1),y(n1),d(m1)

c
c     linear fit

      m = 2
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

c     write(*,198) i,d(i)
  198    format (' coefficient',i5,4x,1p,e16.8)

      enddo
      a1 = d(1)


c
c     calculate best fit of data points
c      write(*,*)
c      write(*,'(3(5x,a5,5x))') 'x','y','fit'
      do i = 1,n

         ax = 0.d0

         do j = 1,m

            ax = ax + a(j,i)*d(j)

         enddo

c         write(*,95) x(i), y(i), ax

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




