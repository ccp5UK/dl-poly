      program statblock
c
c**********************************************************************
c
c     blocking method for estimating standard error of mean
c
c     (reference - Flyvbjerg and Petersen  JCP 91 (1989) 461)
c
c     copyright daresbury laboratory 1993
c     author - w. smith nov 1992
c
c**********************************************************************
c
      parameter (ndeg=13,mxrows=10000,mxcols=20)
      implicit real*8 (a-h,o-z)
      dimension list(mxcols),record(mxcols),d(mxrows,mxcols)
      dimension aver(mxcols),sigma(ndeg,mxcols)
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
         read(7,'(a)',end=100)
         read(7,'(5e14.6)')(record(j),j=1,mcols)
         do k=1,ncols
            d(i,k)=record(list(k))
         enddo
      enddo
  100 m=min(i-1,mxrows)
      write(*,*)'number of data points read',m
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
      write(*,'(1p,20e12.4)')(aver(j),j=1,ncols)
c
c     calculate errors using blocking method
      write(*,*)'calculated standard errors'
      do l=1,ncols
         sigma(1,l)=0.d0
         do  i=1,m
            d(i,l)=d(i,l)-aver(l)
            sigma(1,l)=sigma(1,l)+d(i,l)**2
         enddo
         sigma(1,l)=sqrt(sigma(1,l)/(dble(m)*dble(m-1)))
         mmm=m
         do n=2,ndeg
            j=0
            mmm=mmm/2
            sigma(n,l)=0.d0
            do i=1,mmm
               j=j+2
               d(i,l)=0.5d0*(d(j-1,l)+d(j,l))
               sigma(n,l)=sigma(n,l)+d(i,l)**2
            enddo
            sigma(n,l)=sqrt(sigma(n,l)/(dble(mmm)*dble(mmm-1)))
            if(mmm.le.3)go to 200
         enddo
  200    continue
      enddo
      do  i=1,min(n,ndeg)
         write(*,'(1p,20e12.4)')(sigma(i,j),j=1,ncols)
      enddo
c
c     calculate uncertainties in errors
      write(*,*)'error of standard error'
      do l=1,ncols
         mmm=2*m
         do  n=1,ndeg
            j=0
            mmm=mmm/2
            sigma(n,l)=sigma(n,l)/sqrt(dble(2*(mmm-1)))
            if(mmm.le.3)go to 300
         enddo
  300    continue
      enddo
      do i=1,min(n,ndeg)
         write(*,'(1p,20e12.4)')(sigma(i,j),j=1,ncols)
      enddo
      stop
      end




