      program vacf

c*********************************************************************
c
c     dl_poly routine to calculate velocity autocorrelation function 
c     for selected atoms
c
c     copyright daresbury laboratory 1995
c     author  w.smith jan 1995
c
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c
c*********************************************************************

      implicit real*8(a-h,o-z)
      
      parameter (mxatms=1080,npts=128,nstr=250)
      character*8 atname
      character*8 atmnam
      character*40 fname
      character*80 title
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension vacf(npts,nstr),vcr0(npts,nstr,3),vcf(npts)
      dimension idv(npts),nvm(npts),key(npts)
      complex*16 aaa(npts),www(npts)
      
      write(*,'(a)')'Velocity Autocorrelation Program'
      
c     read the control variables

c     number of atoms
      read(*,*)natms
c     working length of vacf array
      read(*,*)nvacf
c     interval between origins
      read(*,*)iovacf
c     time interval between successive configurations
      read(*,*)tstep
c     name of species for autocorrelation
      read(*,*)atname
c     name of selected HISTORY file
      read(*,*)fname
c     number of configurations in HISTORY file
      read(*,*)nconf

c     check on specified control variables

      if(nvacf.gt.npts)then

         write(*,'(a,2i5)')'error - vacf array too large',nvacf,npts
         stop

      endif
      if(mod(nvacf,iovacf).ne.0)then

         nvacf=iovacf*(nvacf/iovacf)
         write(*,'(a,i5)')'warning - vacf array reset to ',nvacf

      endif

      novacf=nvacf/iovacf

c     initialise velocity autocorrelation arrays

      lor=0
      mor=0
      nsvacf=0
      do j=1,npts
         nvm(j)=0
         do i=1,nstr
            vacf(j,i)=0.d0
            vcr0(j,i,1)=0.d0
            vcr0(j,i,2)=0.d0
            vcr0(j,i,3)=0.d0
         enddo
      enddo
            
c     open HISTORY file and read headers

      open (10,file=fname,form='formatted')

      read(10,'(a80)')title
      write(*,'(a80)')title
      read(10,*)levcfg,imcon

      if(levcfg.eq.0)then

         write(*,'(a)')'error - HISTORY file contains no velocities'
         stop

      endif

      do iconf=1,nconf

         read(10,*,end=100)

         if(imcon.gt.0)then

            read(10,*,end=100)
            read(10,*,end=100)
            read(10,*,end=100)

         endif

         nsampl=0
         do i=1,natms
               
            read(10,*,end=100) atmnam
            read(10,*) 
            read(10,*) xxx,yyy,zzz
            if(levcfg.gt.1)read(10,*)
            if(atmnam.eq.atname)then

               nsampl=nsampl+1

               if(nsampl.le.nstr)then

                  vxx(nsampl)=xxx
                  vyy(nsampl)=yyy
                  vzz(nsampl)=zzz

               endif

            endif
         
         enddo

         if(nsampl.gt.nstr)then

            write(*,'(a)')'error - too many atoms of required type'
            stop

         endif

c     calculate velocity autocorrelation function
         
         if(mod(nsvacf,iovacf).eq.0)then

            lor=min(lor+1,novacf)
            mor=mod(mor,novacf)+1
            idv(mor)=1

            do i=1,nsampl

               vcr0(mor,i,1)=vxx(i)
               vcr0(mor,i,2)=vyy(i)
               vcr0(mor,i,3)=vzz(i)

            enddo

         endif

         nsvacf=nsvacf+1

         do l=1,lor

            m=idv(l)
            idv(l)=m+1
            nvm(m)=nvm(m)+1

            do i=1,nsampl

               vacf(m,i)=vacf(m,i)+vxx(i)*vcr0(l,i,1)+
     x              vyy(i)*vcr0(l,i,2)+vzz(i)*vcr0(l,i,3)

            enddo

         enddo

      enddo
         
  100 continue

      last=min(nsvacf,nvacf)
         
c     normalise velocity autocorrelation functions
         
      do i=1,nsampl

         rnorm=dble(nvm(1))/vacf(1,i)

         do j=1,last

            vacf(j,i)=rnorm*vacf(j,i)/dble(nvm(j))

         enddo

      enddo
      
c     calculate normalised vacf

      vint=0.0d0
      do j=1,last

         vsum=0.0d0
         do i=1,nsampl
            
            vsum=vsum+vacf(j,i)
            
         enddo

         vcf(j)=vsum/dble(nsampl)
         vint=vint+vcf(j)
         time=tstep*dble(j-1)
         write(*,'(1p,3e15.7)')time,vcf(j),vint

      enddo
      write(*,'(a)')'&'

c     initialise FFT arrays

      do i=1,npts

         if(i.le.nvacf)then

            aaa(i)=cmplx(vcf(i),0.d0)

         else

            aaa(i)=(0.d0,0.d0)

         endif

      enddo
      aaa(1)=0.5d0*aaa(1)

c     perform power spectrum

      call fft(1,-1,npts,key,aaa,www,aaa)
      call fft(0,-1,npts,key,aaa,www,aaa)

      do i=1,npts

         freq=dble(i-1)/(dble(npts)*tstep)
         write(*,'(1p,2e15.7)')freq,real(aaa(i))
         
      enddo
      write(*,'(a)')'&'
      
      stop
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


