      subroutine rotate(xo,yo,zo,x1,y1,z1,x2,y2,z2,alp
     x     ,bet,gam)

c***********************************************************************
c
c     DL_POLY utility to rotate water molecule through the
c     Euler angles (alp,bet,gam) in lab fixed axes system.
c     The first data point is taken as the origin for the rotation
c
c     copyright daresbury laboratory 1993
c     author  - t. forester     july 1993
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      dimension x(2),y(2),z(2)
c
c     create hydrogen positions - bondlength 1.0d0

      x(1) = cos(109.5d0*3.1415926d0/360.d0)
      y(1) = sin(109.5d0*3.1415926d0/360.d0)
      z(1) = 0.d0

      x(2) = x(1)
      y(2) = -y(1)
      z(2) = z(1)
c
c     rotate about z by gam

      cs = cos(gam)
      sn = sin(gam)

      do i = 1,2

         xa = cs*x(i) - sn*y(i)
         ya = sn*x(i) + cs*y(i)
         x(i) = xa
         y(i) = ya

      enddo
c
c     rotate about x by bet

      cs = cos(bet)
      sn = sin(bet)

      do i = 1,2

         ya = cs*y(i) - sn*z(i)
         za = sn*y(i) + cs*z(i)
         y(i) = ya
         z(i) = za

      enddo
c
c     rotate about z by alp

      cs = cos(alp)
      sn = sin(alp)

      do i = 1,2

         xa = cs*x(i) - sn*y(i)
         ya = sn*x(i) + cs*y(i)
         x(i) = xa
         y(i) = ya

      enddo

      x1 = x(1) + xo
      x2 = x(2) + xo
      y1 = y(1) + yo
      y2 = y(2) + yo
      z1 = z(1) + zo
      z2 = z(2) + zo

      return
      end

      subroutine error(idnode,kode)
      implicit real*8(a-h,o-z)
      write(6,10) kode
   10 format(' error ',i5,' called by invert' )
      return
      end










