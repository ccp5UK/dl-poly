       subroutine hqadd
     x   (x0,y0,z0,xa,ya,za,xb,yb,zb,xc,yc,zc,hx,hy,hz,bl)
c
c***********************************************************************
c
c     dlpoly utility to add 1 hydrogen to a tertiary amide site
c
c     copyright daresbury laboratory 1993
c     author -     t forester    feb 1993
c
c***********************************************************************


      implicit real*8(a-h,o-z)
c
c	add one H to x0,y0,z0 with bondlength bl

c
c	define a,b,c plane

	xab = xa - xb
	yab = ya - yb
	zab = za - zb

	xbc = xb - xc
	ybc = yb - yc
	zbc = zb - zc

c     plane normal

	rx = yab*zbc - zab*ybc
        ry = zab*xbc - xab*zbc
        rz = xab*ybc - yab*xbc

	rsq = 1.d0/sqrt(rx*rx +ry*ry +rz*rz)

        rx = rx*rsq
        ry = ry*rsq
        rz = rz*rsq
c
c   x0 to point in plane

	xp = x0 - xa
        yp = y0 - ya
        zp = z0 - za

	dot = xp*rx + yp*ry + zp*rz

	bl = sign(bl,dot)

	hx = x0 + bl*rx
        hy = y0 + bl*ry
        hz = z0 + bl*rz

       return
       end



