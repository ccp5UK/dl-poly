       subroutine hadd1_3(x,y,z,hx,hy,hz,bl)
c                  
c***********************************************************************
c     
c     dlpoly utility to add 1 hydrogen to a tertiary site
c
c     copyright daresbury laboratory 1993
c     author -     t forester    feb 1993
c     modified     t forester    feb 1996
c
c     itt
c     2010-10-30 17:20:53
c     1.3
c     Exp
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      dimension x(*),y(*),z(*)
      dimension hx(*),hy(*),hz(*)

c
c	add one H to x(1),y(1),z(1) with bondlength bl

c
c	define (2),(3),(4) plane

	xab = x(2) - x(3)
	yab = y(2) - y(3)
	zab = z(2) - z(3)

	xbc = x(3) - x(4)
	ybc = y(3) - y(4)
	zbc = z(3) - z(4)

c     plane normal

	rx = yab*zbc - zab*ybc
        ry = zab*xbc - xab*zbc
        rz = xab*ybc - yab*xbc

	rsq = 1.d0/sqrt(rx*rx +ry*ry +rz*rz)

        rx = rx*rsq
        ry = ry*rsq
        rz = rz*rsq
c
c   x(1) to point in plane

	xp = x(1) - x(2)
        yp = y(1) - y(2)
        zp = z(1) - z(2)

	dot = xp*rx + yp*ry + zp*rz

	bl = sign(bl,dot)

	hx(1) = x(1) + bl*rx
        hy(1) = y(1) + bl*ry
        hz(1) = z(1) + bl*rz

       return
       end



