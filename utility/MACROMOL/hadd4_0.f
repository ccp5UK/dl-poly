      subroutine hadd4_0(x,y,z,hx,hy,hz,bl)
c
c***********************************************************************
c
c     dl_poly utility to add 4 bonds to an uncoordinate site
c
c              h1
c              :
c         h4...a...h2
c              :
c              h3
c
c     copyright daresbury laboratory 1996
c     author        t forester   feb 1996
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      dimension x(*),y(*),z(*)
      dimension hx(*),hy(*),hz(*)

c
c     construct first a-h vector along x axis

      hx(1) = x(1) + bl
      hy(1) = y(1)
      hz(1) = z(1)
c
c     use hadd3_1 to get the rest!

      x(2) = hx(1)
      y(2) = hy(1)
      z(2) = hz(1)

      angl = 109.5d0
      call hadd3_1(x,y,z,hx(2),hy(2),hz(2),bl,angl)

      return
      end





