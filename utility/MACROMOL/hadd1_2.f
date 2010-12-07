       subroutine hadd1_2(x,y,z,hx,hy,hz,bl,angd)

c***********************************************************************
c
c     subroutine to add 1 hydrogen onto a 2 coordinate group
c
c                b - a ... H
c                    |
c                    c
c
c     copyright daresbury laboratory 1996
c     author  -    t forester    feb 1996
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
c     convert angle between -180 and 180 degrees

      angd = angd - anint(angd/360.0)*360.0
c
c     construct a-b vector

      xab = x(2) - x(1)
      yab = y(2) - y(1)
      zab = z(2) - z(1)

      rab = xab*xab + yab*yab + zab*zab
      rab = 1.0d0/dsqrt(rab)
      xab = xab*rab
      yab = yab*rab
      zab = zab*rab
c
c     construct a-c vector

      xac = x(3) - x(1)
      yac = y(3) - y(1)
      zac = z(3) - z(1)

      rac = xac*xac + yac*yac + zac*zac
      rac = 1.0d0/dsqrt(rac)
      xac = xac*rac
      yac = yac*rac
      zac = zac*rac
c
c     construct 2 basis vectors in plane perpendicular to b-a-c

      ax = xab + xac
      ay = yab + yac
      az = zab + zac
      ra = dsqrt(ax*ax + ay*ay + az*az)
      if(ra.gt.1.d-4) then
        ax = ax/ra
        ay = ay/ra
        az = az/ra
      else
c
c     bonds  b-a-c are linear!!
        ra = dsqrt(xab*xab + yab*yab + zab*zab)
        ax = zab/ra
        ay = xab/ra
        az = yab/ra
      endif

      bx = yab*zac - zab*yac
      by = zab*xac - xab*zac
      bz = xab*yac - zab*xac
      rb = dsqrt(bx*bx + by*by + bz*bz)
      if(rb.gt.1.d-4) then
        bx = bx*rb
        by = by*rb
        bz = bz*rb
      else
        if(abs(ax).gt.0.5d0) then
          bx = az
          by = -ax
          bz = ay
        elseif(abs(ay).gt.0.5d0) then
          bx = az
          by = ax
          bz = -ay
        else
          bx = -az
          by = ax
          bz = ay
        endif

      endif
      dp = ax*bx + ay*by + az*bz

      bx = bx - dp*ax
      by = by - dp*ay
      bz = bz - dp*az
      rb = sqrt(bx*bx + by*by + bz*bz)

      bx = -bx/rb
      by = -by/rb
      bz = -bz/rb

      ax = -ax
      ay = -ay
      az = -az

c
c     construct a-h vector 

      theta = 90.d0 - angd*0.5d0
      theta = theta*(3.141592653589793d0/180.d0)

      c1 = cos(theta)
      s1 = sin(theta)

      ahx = s1*ax + c1*bx
      ahy = s1*ay + c1*by
      ahz = s1*az + c1*bz

c     normalise a-h bond to desired length

      ahx = ahx*bl
      ahy = ahy*bl
      ahz = ahz*bl

c     get H positions

      hx(1) = x(1) + ahx
      hy(1) = y(1) + ahy
      hz(1) = z(1) + ahz

      return

      end
