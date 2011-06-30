      program slowmo
c
c**********************************************************************
c
c     program to process dl_poly trajectories, stripping out
c     selected frequencies
c
c     copyright daresbury laboratory 1993
c     author    w.smith   october    1993
c
c**********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (mxatms=200,mxconf=4096)
      character*80 source,cfgname
      character*8 name(mxatms)

      dimension kkk(mxconf)
      complex*16 ggg(mxconf),www(mxconf)
      dimension aaa(mxconf),bbb(mxconf),ccc(mxconf)
      dimension xyz(3,mxatms),wrk(3*mxatms),cell(9)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      equivalence (xyz,wrk)

c     read in name of trajectory file

      write(*,*)'enter name of trajectory file'
      read(*,'(a)')source
      open(7,file=source)

c     read number of atoms in each configuration

      write(*,*)'enter number of atoms in each configuration'
      read(*,*)natms
      if(natms.gt.mxatms)then

         write(*,*)'error- too many atoms. recompile program'
         stop

      endif

c     read the zapp number

      write(*,'(a,i6)')'enter the zapp number: between 1 and ',mxconf/2
      read(*,*)nzapp
      if(mxconf.lt.2*nzapp)then
        write(*,*)'error nzapp > mxconf/2'
        stop
      endif

c     read first two lines in trajectory file

      read(7,'(a80)')cfgname
      read(7,'(2i10)')levtrj,imcon

c     open temporary storage file

      open(8,file='pivot.01',form='unformatted')

      do k=1,mxconf

c     read configuration data

         read(7,*,end=100)

         if(imcon.gt.0) read(7,'(3g12.4)') cell

         do i=1,natms

            read(7,'(a8)') name(i)
            read(7,'(3e12.4)') xyz(1,i),xyz(2,i),xyz(3,i)
            if(levtrj.eq.1)then

               read(7,*)

            else if(levtrj.eq.2)then

               read(7,*)
               read(7,*)

            endif

         enddo

c     write data to temporary storage

         write(8)xyz

      enddo

  100 continue
      nconf=k-1
      close(7)
      close(8)

c     number of transforms required

      ntrans=3*natms

c     open temporary storage files

      open(9,file='pivot.02',form='unformatted')

c     initialise FFT routine

      call fft(1,-1,mxconf,kkk,ggg,www,ggg)

      do k=1,ntrans,2

         open(8,file='pivot.01',form='unformatted')

         do m=1,nconf

            read(8)wrk

            if(k.lt.ntrans)then

               ggg(m)=cmplx(wrk(k),wrk(k+1))

            else

               ggg(m)=cmplx(wrk(k),0.d0)

            endif

         enddo

         do m=nconf+1,mxconf

            ggg(m)=cmplx(0.d0,0.d0)

         enddo

c     fourier transform data

         call fft(0,-1,mxconf,kkk,ggg,www,ggg)

c     zapp high frequency modes

         do i=1,mxconf

            if(i.gt.nzapp.and.i.le.mxconf-nzapp)then

               ggg(i)=cmplx(0.d0,0.d0)

            else

               ggg(i)=ggg(i)/dble(mxconf)

            endif

         enddo

c     inverse fourier transform

         call fft(0,1,mxconf,kkk,ggg,www,ggg)

         do i=1,nconf

            aaa(i)=real(ggg(i))
            bbb(i)=imag(ggg(i))

         enddo

         write(9)(aaa(i),i=1,nconf)
         write(9)(bbb(i),i=1,nconf)

         close(8)

      enddo

      close(9)

c     now recreate trajectory file

      open(3,file='New_'//source)
      open(4,file='Dis_'//source)

      write(3,'(a80)')cfgname
      write(3,'(2i10)')0,0

      do k=1,nconf

         open(9,file='pivot.02',form='unformatted')

         write(3,'(a8,i10)')'timestep',k

         do i=1,natms

            read(9)(aaa(j),j=1,nconf)
            read(9)(bbb(j),j=1,nconf)
            read(9)(ccc(j),j=1,nconf)

            write(3,'(a8)') name(i)
            write(3,'(1p,3e12.4)')aaa(k),bbb(k),ccc(k)

            xxx(i)=aaa(k)
            yyy(i)=bbb(k)
            zzz(i)=ccc(k)

         enddo

         close(9)

         call wrdisp(cfgname,natms,name,xxx,yyy,zzz)

      enddo

      close(3)
      close(4)

      end
      subroutine fft(ind,isw,ndiv,key,aaa,wfft,bbb)
c***********************************************************************
c
c     fast fourier transform routine
c
c     copyright daresbury laboratory
c     author w.smith
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
      subroutine wrdisp(cfgname,natms,name,xxx,yyy,zzz)

c
c**********************************************************************
c
c     write out a configuration in xyz format
c
c     copyright daresbury laboratory
c     author w.smith
c
c**********************************************************************
c
      implicit real*8(a-h,o-z)
      character*80 cfgname
      character*8 name(*)
      dimension xxx(*),yyy(*),zzz(*)

      write(4,'(i6)')natms
      write(4,'(a)')cfgname
      do i=1,natms

        write(4,'(a8,3f15.8)')name(i),xxx(i),yyy(i),zzz(i)

      enddo
      return
      end
