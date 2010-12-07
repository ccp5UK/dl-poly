      program cfgcmp
c*********************************************************************
c
c     dl_poly utility program for comparing two CONFIG files
c
c     author w.smith december 1998
c     copyright daresbury laboratory 1998
c
c*********************************************************************

      implicit real*8 (a-h,o-z)
      parameter (mxatms=13824)

      character*8 name,file1,file2
      real*8 cel0(9),cel1(9)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension vx0(mxatms),vy0(mxatms),vz0(mxatms)
      dimension fx0(mxatms),fy0(mxatms),fz0(mxatms)
      dimension xx1(mxatms),yy1(mxatms),zz1(mxatms)
      dimension vx1(mxatms),vy1(mxatms),vz1(mxatms)
      dimension fx1(mxatms),fy1(mxatms),fz1(mxatms)

      write(*,*)'Enter name of first file'
      read(*,*)file1
      write(*,*)'Enter name of second file'
      read(*,*)file2
      open(7,file=file1)
      open(8,file=file2)

      read(7,*)
      read(7,*)
      read(7,*)cel0(1),cel0(2),cel0(3)
      read(7,*)cel0(4),cel0(5),cel0(6)
      read(7,*)cel0(7),cel0(8),cel0(9)

      read(8,*)
      read(8,*)
      read(8,*)cel1(1),cel1(2),cel1(3)
      read(8,*)cel1(4),cel1(5),cel1(6)
      read(8,*)cel1(7),cel1(8),cel1(9)

      celmx=0.d0

      do i=1,9

        celmx=max(celmx,abs(cel0(i)-cel1(i)))

      enddo

      do k=1,mxatms

        read(7,*,end=100)name,k
        read(7,*)xx0(k),yy0(k),zz0(k)
        read(7,*)vx0(k),vy0(k),vz0(k)
        read(7,*)fx0(k),fy0(k),fz0(k)

      enddo
  100 continue
      do k=1,mxatms

        read(8,*,end=200)name,k
        read(8,*)xx1(k),yy1(k),zz1(k)
        read(8,*)vx1(k),vy1(k),vz1(k)
        read(8,*)fx1(k),fy1(k),fz1(k)

      enddo
  200 continue
      erx=0.d0
      ery=0.d0
      erz=0.d0
      evx=0.d0
      evy=0.d0
      evz=0.d0
      efx=0.d0
      efy=0.d0
      efz=0.d0
      
      do i=1,mxatms

        if(erx.lt.abs(xx0(i)-xx1(i)))then
          erx=abs(xx0(i)-xx1(i))
          k1=i
        endif
        if(ery.lt.abs(yy0(i)-yy1(i)))then
          ery=abs(yy0(i)-yy1(i))
          k2=i
        endif
        if(erz.lt.abs(zz0(i)-zz1(i)))then
          erz=abs(zz0(i)-zz1(i))
          k3=i
        endif
        if(evx.lt.abs(vx0(i)-vx1(i)))then
          evx=abs(vx0(i)-vx1(i))
          k4=i
        endif
        if(evy.lt.abs(vy0(i)-vy1(i)))then
          evy=abs(vy0(i)-vy1(i))
          k5=i
        endif
        if(evz.lt.abs(vz0(i)-vz1(i)))then
          evz=abs(vz0(i)-vz1(i))
          k6=i
        endif
        if(efx.lt.abs(fx0(i)-fx1(i)))then
          efx=abs(fx0(i)-fx1(i))
          k7=i
        endif
        if(efy.lt.abs(fy0(i)-fy1(i)))then
          efy=abs(fy0(i)-fy1(i))
          k8=i
        endif
        if(efz.lt.abs(fz0(i)-fz1(i)))then
          efz=abs(fz0(i)-fz1(i))
          k9=i
        endif

      enddo

      write(*,*)'Maximum error in cell vectors'
      write(*,'(1p,e12.4)')celmx
      write(*,*)'Maximum errors in position'
      write(*,'(3i5,1p,3e12.4)')k1,k2,k3,erx,ery,erz
      write(*,*)'Maximum errors in velocity'
      write(*,'(3i5,1p,3e12.4)')k4,k5,k6,evx,evy,evz
      write(*,*)'Maximum errors in force'
      write(*,'(3i5,1p,3e12.4)')k7,k8,k9,efx,efy,efz

      end
