      program cfgchk

c*********************************************************************
c
c     dl_poly utility program for checking a CONFIG file to ensure
c     an absence of close or zero atomic contacts
c
c     copyright daresbury laboratory 1998
c     author w.smith december 1998
c
c*********************************************************************

      implicit real*8 (a-h,o-z)
      parameter (mxatms=30000)

      character*8 name,filnam
      dimension cell(9),rcell(9)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)

      write(*,*)'Enter name of CONFIG file'
      read(*,*)filnam
      open(7,file=filnam)

      read(7,*)
      read(7,*)levcfg,imcon
      if(imcon.gt.0)then

        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)
        call invert(cell,rcell,det)
        if(abs(det).lt.1.d-8)write(*,'(a)')
     x    'Error - zero determinant for cell vectors'

      endif

      do i=1,mxatms

        read(7,*,end=100)name
        read(7,*)xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)read(7,*)vxx,vyy,vzz
        if(levcfg.gt.1)read(7,*)fxx,fyy,fzz

      enddo

  100 continue

      rchk=1.d10
      natms=i-1
      last=natms
      mlast=natms/2
      nlast=(natms-1)/2

      write(*,'(a)')'Checking pair separations.......'

      do m=1,mlast

        if(m.gt.nlast)last=mlast

        do i=1,last

          j=i+m
          if(j.gt.natms)j=j-natms

          ddx=xxx(i)-xxx(j)
          ddy=yyy(i)-yyy(j)
          ddz=zzz(i)-zzz(j)
          if(imcon.gt.0)then

            sdx=rcell(1)*ddx+rcell(4)*ddy+rcell(7)*ddz
            sdy=rcell(2)*ddx+rcell(5)*ddy+rcell(8)*ddz
            sdz=rcell(3)*ddx+rcell(6)*ddy+rcell(9)*ddz
            sdx=sdx-nint(sdx)
            sdy=sdy-nint(sdy)
            sdz=sdz-nint(sdz)
            ddx=cell(1)*sdx+cell(4)*sdy+cell(7)*sdz
            ddy=cell(2)*sdx+cell(5)*sdy+cell(8)*sdz
            ddz=cell(3)*sdx+cell(6)*sdy+cell(9)*sdz

          endif

          rsq=ddx**2+ddy**2+ddz**2

          rchk=min(rchk,rsq)

          if(rsq.lt.1.d-8)then

            write(*,'(i6,1p,3e12.4)')i,xxx(i),yyy(i),zzz(i)
            write(*,'(i6,1p,3e12.4)')j,xxx(j),yyy(j),zzz(j)
            write(*,'(1p,3e12.4)')rsq,sqrt(rsq)

          endif

        enddo

      enddo

      write(*,'(a)')'Checking complete.'
      write(*,'(a,1p,e12.4)')'Minimum separation found:',sqrt(rchk)


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

      real*8 a,b,d,r

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
      end

