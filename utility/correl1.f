      program correl1
c     
c***********************************************************************
c     
c     dl_poly utility program for calculating correlation functions
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)
      
      parameter (ndiv=500,mxcols=15,mxcorr=10)
      
      character*4 name
      
      dimension sample(mxcols),c(ndiv),u(ndiv),v(ndiv),t(ndiv)
      dimension icol(mxcorr),jcol(mxcorr)
      
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
c     now process each correlation function in turn
      
      do ipairs=1,npairs
         
c     
c     open the file containing the columns of data
         
         open(7,file='input.data')
         
c     
c     now read the required columns (note first column is time)
         
         do i=1,ndiv
            read(7,*,end=100)t(i),(sample(j),j=1,ncols)
            u(i)=sample(icol(ipairs))
            v(i)=sample(jcol(ipairs))
         enddo
         
  100    m=i-1
         
c     
c     close data file
         
         close(7)
         
         dtim=t(2)-t(1)
         
c     
c     calculate column  averages
         
         sumu=0.d0
         sumv=0.d0
         
         do i=1,m
            sumu=sumu+u(i)
            sumv=sumv+v(i)
         enddo
         
         sumu=sumu/dble(m)
         sumv=sumv/dble(m)
         
c     
c     subtract averages from columns
         
         do i=1,m
            u(i)=u(i)-sumu
            v(i)=v(i)-sumv
         enddo
         
c     
c     calculate correlation function
         
         do i=1,m
            c(i)=0.d0
            do j=1,m+1-i
               c(i)=c(i)+u(j)*v(i+j-1)
            enddo
         enddo
         
c     
c     normalise the correlation function
         
         do i=1,m
            u(i)=c(i)/dble(m+1-i)
         enddo
         
c     
c     open output file to receive correlation function
         
         namfil=100*icol(ipairs)+jcol(ipairs)
         write(name,'(i4.4)')namfil
         open(8,file='x'//name)
         
c     
c     write out correlation function
         
         do i=1,m
            tim=dble(i-1)*dtim
            write(8,'(f10.6,1p,e15.7)')tim,u(i)
         enddo
         
c     
c     close correlation function file
         
         close(8)
         
      enddo
      
      stop
      end
