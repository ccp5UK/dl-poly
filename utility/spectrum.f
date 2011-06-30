      program spectrum
c
c***********************************************************************
c
c     dl_poly utility program for spectral analysis using the fast
c     fourier transform
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith june 1992.
c
c     note - this program includes a naff fft routine purely for
c     portability and verification reasons. users are strongly
c     recommended to substitute a machine specific fft, which
c     should vastly improve performance.
c
c     note - you will probably need to change some of the array
c     dimensions using the following parameters.
c
c     note - ndiv must be a power of 2. ndiv2 will change according to
c     the fft routine you use.
c
c***********************************************************************
c

      implicit real*8(a-h,o-z)

      parameter (ndiv=8192,ndiv2=ndiv)
      parameter (mxcols=15,mxspec=10)

      character*2 name

      complex*16 g(ndiv),w(ndiv2),cmplx,sum
      dimension sample(mxcols),u(ndiv)
      dimension icol(mxspec),key(ndiv)

      data pi/3.1415926536d0/

c
c     define the window function (blackman-harris 4 term)

      win(x)=0.35875d0-0.48829d0*cos(x)+0.14128d0*cos(2.d0*x)-
     x       0.01168d0*cos(3.d0*x)


c
c     declare how many fourier transforms are required
c     and how many columns of data there in the input data
c     (excluding the time column)

      write(*,'(a)')'enter number of required transforms'
      read(*,*)nspecs
      write(*,'(a)')'enter number of data columns in input'
      read(*,*)ncols
      write(*,'(a)')'enter data stride length'
      read(*,*)nstrid

      ncols=min(mxcols,ncols)
      nspecs=min(mxspec,nspecs)

c
c     specify the columns to be fourier transformed

      write(*,'(a)')'specify columns to be transformed'
      read(*,*)(icol(i),i=1,nspecs)

c
c     initialise the fft work arrays

      call fft(1,1,ndiv,key,g,w,g)

c
c     now process each fourier transform in turn

      do ispecs=1,nspecs

c
c     zero the complex data array

         do i=1,ndiv
            g(i)=(0.d0,0.d0)
         enddo

c
c     open the file containing the columns of data

         open(7,file='input.data')

c
c     now read the required columns (note first column is time)

         k=0
         sum=(0.d0,0.d0)
         do i=1,ndiv
            read(7,*,end=100)t,(sample(j),j=1,ncols)
            if(i.eq.1)t0=t
            if(mod(i-1,nstrid).eq.0)then
               k=k+1
               g(k)=cmplx(sample(icol(ispecs)),0.d0)
               sum=sum+g(k)
            endif
         enddo

  100    m=k
         sum=sum/dble(m)

c
c     close data file

         close(7)

         dtim=(t-t0)/dble(m-1)
         period=dble(m)*dtim

c
c     apply window function

         do i=1,m
            g(i)=(g(i)-sum)*win((2.d0*pi*dtim/period)*dble(i-1))
         enddo

c
c     fourier transform the data

         call fft(0,-1,ndiv,key,g,w,g)

c
c     normalise the fourier transform

         do i=1,ndiv
            g(i)=dtim*g(i)
            u(i)=sqrt(dreal(g(i))**2+dimag(g(i))**2)
         enddo

c
c     open output file to receive fourier transform

         write(name,'(i2.2)')icol(ispecs)
         open(8,file='ft'//name)

c
c     write out fourier transform

         do i=1,ndiv/2
            fre=dble(i-1)/(dtim*dble(ndiv))
            write(8,'(1p,4e12.4)')fre,u(i),g(i)
         enddo

c
c     close fourier transform file

         close(8)

      enddo

      stop
      end

      subroutine fft(ind,isw,ndiv,key,aaa,wfft,bbb)
c***********************************************************************
c
c     fast fourier transform routine
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith june 1992.
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      logical check
      complex*16 aaa(ndiv),bbb(ndiv),wfft(ndiv),ttt
      dimension key(ndiv)
      data tpi/6.2831853072d0/
   10 format(1h0,'error - number of points not a power of two')

c
c     check that array is of suitable length
      nt=1
      check=.true.
      do i=1,20
         nt=2*nt
         if(nt.eq.ndiv)then
            check=.false.
            nu=i
         endif
      enddo
      if(check)then
         write(*,10)
         stop
      endif

      if(ind.gt.0)then
c
c     set reverse bit address array

         do kkk=1,ndiv
            iii=0
            jjj=kkk-1
            do j=1,nu
               jj2=jjj/2
               iii=2*(iii-jj2)+jjj
               jjj=jj2
            enddo
            key(kkk)=iii+1
         enddo
c
c     initialise complex exponential factors

         tpn=tpi/dble(ndiv)
         arg=0.d0
         np1=ndiv+1
         np2=ndiv/2
         wfft(1)=(1.d0,0.d0)
         do i=1,np2
            arg=tpn*dble(i)
            wfft(i+1)=cmplx(cos(arg),sin(arg))
            wfft(np1-i)=conjg(wfft(i+1))
         enddo

         return
      endif

c
c     take conjugate of exponentials if required

      if(isw.lt.0)then

         do i=1,ndiv
            wfft(i)=conjg(wfft(i))
         enddo

      endif

c
c     take copy input array

      do i=1,ndiv
         bbb(i)=aaa(i)
      enddo
c
c     perform fourier transform

      kkk=0
      nu1=nu-1
      np2=ndiv/2
      do l=1,nu

  100    do i=1,np2
            iii=key(kkk/2**nu1+1)
            kk1=kkk+1
            k12=kk1+np2
            ttt=bbb(k12)*wfft(iii)
            bbb(k12)=bbb(kk1)-ttt
            bbb(kk1)=bbb(kk1)+ttt
            kkk=kkk+1
         enddo
         kkk=kkk+np2
         if(kkk.lt.ndiv)go to 100
         kkk=0
         nu1=nu1-1
         np2=np2/2

      enddo
c
c     unscramble the fft using bit address array

      do kkk=1,ndiv
         iii=key(kkk)
         if(iii.gt.kkk)then
            ttt=bbb(kkk)
            bbb(kkk)=bbb(iii)
            bbb(iii)=ttt
         endif
      enddo
c
c     restore exponentials to unconjugated values if necessary

      if(isw.lt.0)then

         do i=1,ndiv
            wfft(i)=conjg(wfft(i))
         enddo

      endif

      return
      end
