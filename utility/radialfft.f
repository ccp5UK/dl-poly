      program lado
c
c**********************************************************************
c
c     dl_poly routine for 3d radial fourier tranform
c
c     copyright daresbury laboratory 1994
c     author w.smith
c
c**********************************************************************

      implicit real*8(a-h,o-z)
      parameter (pi=3.141592653589793d0)
      parameter (mxd=1000)
      dimension aaa(mxd),bbb(mxd)
      data a0,a1,a2/0.42d0,0.50d0,0.08d0/

c
c     read in data

      do i=1,mxd

         read(*,*,end=100)bbb(i),aaa(i)
         aaa(i)=aaa(i)-1.d0

      enddo
  100 nnn=i-1
      do i=1,nnn-1

         bbb(i)=0.5d0*(bbb(i)+bbb(i+1))
         aaa(i)=0.5d0*(aaa(i)+aaa(i+1))

      enddo
      nnn=nnn-1
c
c     set up window function (blackman function)
c$$$      arg=pi/dble(nnn)
c$$$      do i=1,nnn
c$$$
c$$$         ccc=cos(arg*dble(i+nnn-1))
c$$$         aaa(i)=(a0-a1*ccc+a2*(2.d0*ccc**2-1.d0))*aaa(i)
c$$$
c$$$      enddo
c$$$      aaa(1)=aaa(1)/2.d0

      delr=(bbb(nnn)-bbb(1))/dble(nnn-1)
      write(*,'(a,i6)')'number of data points = ',nnn
      write(*,'(a,1pe12.4)')'radial increment delr = ',delr
      do i=nnn+1,500

         aaa(i)=0.d0

      enddo
      nnn=500

      isw=1
      call radfft(isw,nnn,delr,aaa,bbb)


      delk=pi/(delr*dble(nnn))
      do i=1,nnn

         write(*,'(1p,2e14.6)')delk*dble(i),bbb(i)

      enddo

      stop

      end
      subroutine radfft(isw,nnn,delr,aaa,bbb)
c***********************************************************************
c
c     dl_poly 3D radial fourier transform routine using lado's method
c     reference: j. comput. phys. 8 (1971) 417
c
c     copyright daresbury laboratory 1994
c     author w smith
c
c     note: first data point is i*delr not 0
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      parameter (pi=3.141592653589793d0)

      dimension aaa(nnn),bbb(nnn)

c
c     perform fourier transform

      sw=pi*dble(isw)/dble(nnn)

      do j=1,nnn

         bbb(j)=0.d0

         do i=1,nnn

            bbb(j)=dble(i)*aaa(i)*sin(sw*dble(j)*dble(i))+bbb(j)

         enddo

         bbb(j)=(4.d0*dble(nnn)*delr**3/dble(j))*bbb(j)

      enddo

      return
      end
