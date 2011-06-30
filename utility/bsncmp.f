      program bsncmp
c*********************************************************************
c
c     dl_poly utility program for comparing two BASIN files
c
c     copyright daresbury laboratory
c     author w.smith march 2008
c
c*********************************************************************

      implicit none

      integer, parameter:: mxatms=998

      character*8 name(mxatms)
      character*10 file0,file1
      integer i,lev0,lev1,imcon0,imcon1,natms
      real(8) disp,det,xxx,yyy,zzz,ddx,ddy,ddz,ddd,txx,tyy,tzz
      real(8) ccc(9),cel0(9),cel1(9),rcel0(9),rcel1(9)
      real(8) xx0(mxatms),yy0(mxatms),zz0(mxatms)
      real(8) xx1(mxatms),yy1(mxatms),zz1(mxatms)

      write(*,*)'Enter name of first basin file'
      read(*,*)file0
      write(*,*)'Enter name of second basin file'
      read(*,*)file1
      write(*,*)'Enter max allowed displacement'
      read(*,*)disp

c     open BASIN files

      open(7,file=file0)
      open(8,file=file1)

      read(7,*)
      read(7,*)lev0,imcon0,natms
      if(imcon0.gt.0)then
        read(7,*)cel0(1),cel0(2),cel0(3)
        read(7,*)cel0(4),cel0(5),cel0(6)
        read(7,*)cel0(7),cel0(8),cel0(9)
      endif
      call invert(cel0,rcel0,det)

      read(8,*)
      read(8,*)lev1,imcon1,natms
      if(imcon1.gt.0)then
        read(8,*)cel1(1),cel1(2),cel1(3)
        read(8,*)cel1(4),cel1(5),cel1(6)
        read(8,*)cel1(7),cel1(8),cel1(9)
      endif
      call invert(cel1,rcel1,det)

      do i=1,9
        ccc(i)=0.5d0*(cel0(i)+cel1(i))
      enddo

      do i=1,natms

        read(7,*,end=100)name(i)
        read(7,*)xxx,yyy,zzz
        xx0(i)=rcel0(1)*xxx+rcel0(4)*yyy+rcel0(7)*zzz
        yy0(i)=rcel0(2)*xxx+rcel0(5)*yyy+rcel0(8)*zzz
        zz0(i)=rcel0(3)*xxx+rcel0(6)*yyy+rcel0(9)*zzz

      enddo
  100 continue
      do i=1,natms

        read(8,*,end=200)name(i)
        read(8,*)xxx,yyy,zzz
        xx1(i)=rcel1(1)*xxx+rcel1(4)*yyy+rcel1(7)*zzz
        yy1(i)=rcel1(2)*xxx+rcel1(5)*yyy+rcel1(8)*zzz
        zz1(i)=rcel1(3)*xxx+rcel1(6)*yyy+rcel1(9)*zzz

      enddo
  200 continue

      write(*,*)'structural differences'

      do i=1,natms

        ddx=xx0(i)-xx1(i)
        ddy=yy0(i)-yy1(i)
        ddz=zz0(i)-zz1(i)
        txx=ddx-cel0(1)*nint(ddx/cel0(1))
        tyy=ddy-cel0(5)*nint(ddy/cel0(5))
        tzz=ddz-cel0(9)*nint(ddz/cel0(9))
        ddx=ccc(1)*txx+ccc(4)*tyy+ccc(7)*tzz
        ddy=ccc(2)*txx+ccc(5)*tyy+ccc(8)*tzz
        ddz=ccc(3)*txx+ccc(6)*tyy+ccc(9)*tzz
        ddd=ddx**2+ddy**2+ddz**2
        if(ddd.gt.disp**2)then
          write(*,'(a8,i10,3f16.8)')name(i),i,ddx,ddy,ddz
        endif

      enddo

      write(*,*)'Analysis complete'

      end
      subroutine invert(a,b,d)

c***********************************************************************
c
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert
