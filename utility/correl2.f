      program correl2
c     
c***********************************************************************
c     
c     dl_poly utility program for calculating correlation functions
c     using the fast fourier transform strategy
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     
c     note - this program includes a naff fft routine purely for
c     portability and verification reasons. users are strongly
c     recommended to substitute a machine specific fft, which
c     should vastly improve performance.
c     
c     note - you will probably need to change some of the array
c     dimensions using the following parameters. 
c     
c     note - ndiv must be a power of 2. ndiv3 will change according to
c     the fft routine you use. leave ndiv2 as it is.
c     
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c
c***********************************************************************
c     

      implicit real*8(a-h,o-z)

      parameter (ndiv=8192,ndiv2=2*ndiv,ndiv3=2*ndiv)
      parameter (mxcols=15,mxcorr=10)
      
      character*4 name
      
      complex*16 g(ndiv2),c(ndiv2),w(ndiv3),ave,cmplx
      dimension sample(mxcols),u(ndiv),t(ndiv)
      dimension icol(mxcorr),jcol(mxcorr),key(ndiv2)
      
c     
c     declare how many correlation functions are required
c     and how many columns of data there in the input data
c     (excluding the time column)
      
      write(*,'(a)')'enter number of required correlation functions'
      read(*,*)npairs
      write(*,'(a)')'enter number of data columns in input file'
      read(*,*)ncols
      ncols=min(mxcols,ncols)
      npairs=min(mxcorr,npairs)
      
c     
c     specify which columns in the input data are to be correlated
      
      write(*,'(a)')'enter ids of column pairs to be correlated'
      read(*,*)(icol(i),jcol(i),i=1,npairs)
      
c     
c     initialise the fft work arrays
      
      call fft(1,1,ndiv2,key,g,w,g)
      
c     
c     now process each correlation function in turn
      
      do ipairs=1,npairs
         
c     
c     zero the complex data array
         
         do i=1,ndiv2
            g(i)=(0.d0,0.d0)
         enddo
         
c     
c     open the file containing the columns of data
         
         open(7,file='input.data')
         
c     
c     now read the required columns (note first column is time)
         
         do i=1,ndiv
            read(7,*,end=100)t(i),(sample(j),j=1,ncols)
            g(i)=cmplx(sample(icol(ipairs)),sample(jcol(ipairs)))
         enddo
         
  100    m=i-1
         
c     
c     close data file
         
         close(7)
         
         dtim=t(2)-t(1)
         
c     
c     calculate column  average
         
         ave=(0.d0,0.d0)
         
         do i=1,m
            ave=ave+g(i)
         enddo
         
         ave=ave/dble(m)
         
c     
c     subtract average from columns
         
         do i=1,m
            g(i)=g(i)-ave
         enddo
         
c     
c     fourier transform the data
         
         call fft(0,-1,ndiv2,key,g,w,g)
         
c     
c     calculate correlation as product in frequency domain
         
         ur=real(g(1))
         ui=0.d0
         vr=aimag(g(1))
         vi=0.d0
         c(1)=cmplx(ur*vr+ui*vi,ur*vi-ui*vr)
         
         do i=2,ndiv2
            ur=0.5d0*(real(g(i))+real(g(ndiv2+2-i)))
            ui=0.5d0*(aimag(g(i))-aimag(g(ndiv2+2-i)))
            vr=0.5d0*(aimag(g(i))+aimag(g(ndiv2+2-i)))
            vi=0.5d0*(-real(g(i))+real(g(ndiv2+2-i)))
            c(i)=cmplx(ur*vr+ui*vi,ur*vi-ui*vr)
         enddo
         
c     
c     inverse fourier transform
         
         call fft(0,1,ndiv2,key,c,w,c)
         
c     
c     normalise the correlation function
         
         rnorm=1.d0/dble(ndiv2)
         do i=1,m
            u(i)=rnorm*real(c(i))/dble(m+1-i)
         enddo
         
c     
c     open output file to receive correlation function
         
         namfil=100*icol(ipairs)+jcol(ipairs)
         write(name,'(i4.4)')namfil
         open(8,file='c'//name)
         
c     
c     write out correlation function
         
         do i=1,m
            tim=dble(i-1)*dtim
            write(8,'(2f10.6)')tim,u(i)
         enddo
         
c     
c     close correlation function file
         
         close(8)
         
      enddo
      
      stop
      end
      
      subroutine fft(ind,isw,ndiv,key,aaa,wfft,bbb)
c***********************************************************************
c     
c     fast fourier transform routine
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
