      subroutine hadd3_1(x,y,z,hx,hy,hz,bl,angd)
c
c***********************************************************************
c
c     dl_poly utility to add 3 hydrogens to a 1 coordinate site
c
c              h1
c              :
c            b-a...h2
c              :
c              h3
c
c     copyright daresbury laboratory 1993
c     author        t forester   feb 1993
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      dimension x(2),y(2),z(2)
      dimension hx(3),hy(3),hz(3)
c
c     convert angle between -180 and 180 degrees

      angd = angd - nint(angd/360.d0)*360.d0

c
c     convert to radians

      ang = angd*3.1415926d0/180.0d0
      csa = cos(ang)
c
c     construct a-b vector

      xab = x(2) - x(1)
      yab = y(2) - y(1)
      zab = z(2) - z(1)

      rab = xab*xab + yab*yab + zab*zab
      rab = 1.0d0/sqrt(rab)
      xab = xab*rab
      yab = yab*rab
      zab = zab*rab

c     find which way bond is oriented
      ip = 0
      if (abs(yab).gt.abs(zab)) then
         if (abs(xab).gt.abs(yab)) then
            ip = 2
         else
            ip = 1
         endif
      else if (abs(xab).gt.abs(zab)) then
         ip = 2
      endif

      do i = 1,ip
         temp = zab
         zab = yab
         yab = xab
         xab = temp
      enddo

c
c     construct a-h vector in transformed x-y plane

      ahx = sqrt(0.5d0)*sign(1.d0,xab)*sign(1.d0,90.d0-abs(angd))
      ahy = sqrt(0.5d0)*sign(1.d0,yab)*sign(1.d0,90.d0-abs(angd))

c
c     solve quadratic for z component

      aa = xab*ahx + yab*ahy
      c = aa*aa - csa*csa
      b = 2.0d0*aa
      a = 1.d0 - (csa/zab)**2
      if (abs(angd).ge.90.d0) then
         ahz = (-b - sqrt(b*b - 4.d0*a*c))/(2.d0*a*zab)
      else
         ahz = (-b + sqrt(b*b - 4.d0*a*c))/(2.d0*a*zab)
      endif

c
c     normalise bond to desired length

      scale = bl/sqrt(1.d0 + ahz*ahz)
      ahx = ahx*scale
      ahy = ahy*scale
      ahz = ahz*scale

      do i = 1,ip
         temp = ahx
         ahx = ahy
         ahy = ahz
         ahz = temp
         temp = xab
         xab = yab
         yab = zab
         zab = temp
      enddo
c
c     first H position

      hx(1) = x(1) + ahx
      hy(1) = y(1) + ahy
      hz(1) = z(1) + ahz

c     generate basis vectors perpendicular to b-a vector

      rd = 1.0d0/sqrt(ahx*ahx + ahy*ahy + ahz*ahz)
      dox = ahx*rd
      doy = ahy*rd
      doz = ahz*rd

      ax = yab*doz - zab*doy
      ay = zab*dox - xab*doz
      az = xab*doy - yab*dox
      ra = 1.0d0/sqrt(ax*ax + ay*ay + az*az)
      ax = ax*ra
      ay = ay*ra
      az = az*ra

      bx = ay*zab - az*yab
      by = az*xab - ax*zab
      bz = ax*yab - ay*xab

c     project first a-h onto b-a unit vector

      dpn = ahx*xab + ahy*yab + ahz*zab

c     project first a-h onto b axis

      dpb = ahx*bx + ahy*by + ahz*bz

c     rotate around b-a bond by 120 degrees

      cs = -0.5d0
      sn = sqrt(3.0d0)*0.5d0

      b1 = dpb*bx
      b2 = dpb*by
      b3 = dpb*bz

      a1 = dpb*ax
      a2 = dpb*ay
      a3 = dpb*az

      rx = cs*b1 + sn*a1
      ry = cs*b2 + sn*a2
      rz = cs*b3 + sn*a3

      ahx = dpn*xab + rx
      ahy = dpn*yab + ry
      ahz = dpn*zab + rz

      hx(2) = x(1) + ahx
      hy(2) = y(1) + ahy
      hz(2) = z(1) + ahz

c     rotate around b-a bond by -120 degrees for third H

      sn = -sn

      rx = cs*b1 + sn*a1
      ry = cs*b2 + sn*a2
      rz = cs*b3 + sn*a3

      ahx = dpn*xab + rx
      ahy = dpn*yab + ry
      ahz = dpn*zab + rz

      hx(3) = x(1) + ahx
      hy(3) = y(1) + ahy
      hz(3) = z(1) + ahz

      return

      end

