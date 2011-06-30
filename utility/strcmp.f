      program strcmp
c*********************************************************************
c
c     dl_poly utility program for comparing two CONFIG files
c     compares configurations only
c
c     author w.smith september 2006
c     copyright daresbury laboratory
c
c*********************************************************************

      implicit real*8 (a-h,o-z)
      parameter (mxatms=13824)

      logical same
      character*8 name,file0,file1
      real*8 cel0(9),cel1(9),rcel0(9),rcel1(9),cell(9)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension xx1(mxatms),yy1(mxatms),zz1(mxatms)

c     read instructions

      write(*,*)'Enter name of first file'
      read(*,*)file0
      write(*,*)'Enter name of second file'
      read(*,*)file1
      write(*,*)'Enter required tolerance'
      read(*,*)rcut
      rct2=rcut**2

c     read first file

      open(7,file=file0)
      read(7,*)
      read(7,*)lev1,imcon
      if(imcon.gt.0)then
        read(7,*)cel0(1),cel0(2),cel0(3)
        read(7,*)cel0(4),cel0(5),cel0(6)
        read(7,*)cel0(7),cel0(8),cel0(9)
        call invert(cel0,rcel0,det)
      endif
      do k=1,mxatms
        read(7,*,end=100)name
        read(7,*)xx0(k),yy0(k),zz0(k)
        if(lev1.gt.0)read(7,*)
        if(lev1.gt.1)read(7,*)
      enddo
  100 continue
      close (7)

c     read second file

      open(8,file=file1)
      read(8,*)
      read(8,*)lev2,imcon
      if(imcon.gt.0)then
        read(8,*)cel1(1),cel1(2),cel1(3)
        read(8,*)cel1(4),cel1(5),cel1(6)
        read(8,*)cel1(7),cel1(8),cel1(9)
        call invert(cel1,rcel1,det)
      endif
      do k=1,mxatms
        read(8,*,end=200)name
        read(8,*)xx1(k),yy1(k),zz1(k)
        if(lev2.gt.0)read(8,*)
        if(lev2.gt.1)read(8,*)
      enddo
  200 continue
      close(8)
      natms=k-1

c     average cell matrix

      do i=1,9
        cell(i)=0.5d0*(cel0(i)+cel1(i))
      enddo

c     compare structures

      ddd=0.d0
      do i=1,natms

        xxa=rcel0(1)*xx0(i)+rcel0(4)*yy0(i)+rcel0(7)*zz0(i)
        yya=rcel0(2)*xx0(i)+rcel0(5)*yy0(i)+rcel0(8)*zz0(i)
        zza=rcel0(3)*xx0(i)+rcel0(6)*yy0(i)+rcel0(9)*zz0(i)

        xxb=rcel1(1)*xx1(i)+rcel1(4)*yy1(i)+rcel1(7)*zz1(i)
        yyb=rcel1(2)*xx1(i)+rcel1(5)*yy1(i)+rcel1(8)*zz1(i)
        zzb=rcel1(3)*xx1(i)+rcel1(6)*yy1(i)+rcel1(9)*zz1(i)

        ssx=xxa-xxb
        ssy=yya-yyb
        ssz=zza-zzb

        xss=ssx-nint(ssx)
        yss=ssy-nint(ssy)
        zss=ssz-nint(ssz)

        dxx=cell(1)*xss+cell(4)*yss+cell(7)*zss
        dyy=cell(2)*xss+cell(5)*yss+cell(8)*zss
        dzz=cell(3)*xss+cell(6)*yss+cell(9)*zss
        ddd=max(ddd,dxx*dxx+dyy*dyy+dzz*dzz)

      enddo

      same=(ddd.le.rct2)

      if(same)then
        write(*,*)file0//' and '//file1//' same'
      else
        write(*,*)file0//' and '//file1//' different'
      endif

      end
      subroutine invert(a,b,d)
c
c***********************************************************************
c
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c
c***********************************************************************
c

      real*8 a,b,d,r

      dimension a(9),b(9)
c
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
c
c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d
c
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
      end

