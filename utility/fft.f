c**********************************************************************
c
c     test program for fft routine 
c     fourier transforms a cosine function
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c
c**********************************************************************

      implicit real*8 (a-h,o-z)

      parameter (ndiv=64)

      complex*16 aaa(ndiv),bbb(ndiv),wfft(ndiv)
      dimension key(ndiv)
      data tpi/6.2831853072d0/

c     set up the cosine function (with period 1/4 of window width)

      tpn=tpi/dble(ndiv)
      np1=ndiv+1
      np2=ndiv/2
      aaa(1)=(1.d0,0.d0)
      do i=1,np2
        arg=tpn*dble(4*i)
        aaa(i+1)=cmplx(cos(arg),0.d0)
        aaa(np1-i)=conjg(aaa(i+1))
      enddo
      
c     initialise the fft routine

      call fft(1,1,ndiv,key,aaa,wfft,bbb)

c     perform FT of function

      call fft(0,1,ndiv,key,aaa,wfft,bbb)

c     print result

      do i=1,ndiv

        write(*,'(i5,1p,2e12.4)')i,real(aaa(i)),real(bbb(i))

      enddo

      end
      subroutine fft(ind,isw,ndiv,key,aaa,wfft,bbb)
c***********************************************************************
c     
c     fast fourier transform routine
c     
c     copyright daresbury laboratory 1994
c
c     author w smith
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
         
         np1=ndiv+1
         np2=ndiv/2
         tpn=tpi/dble(ndiv)
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
      np2=ndiv/2
      do l=1,nu

  100    do i=1,np2
            iii=key(kkk/np2+1)
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
