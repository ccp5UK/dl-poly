      subroutine hadd2_1(x,y,z,hx,hy,hz,bl,angd,angt)

c***********************************************************************
c
c     subroutine to add 2 bonds to a 1 coordinate site
c
c     copyright daresbury laboratory 1996
c     author -     t forester    feb 1996
c
c     itt
c     2010-10-30 17:20:53
c     1.3
c     Exp
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      dimension x(2),y(2),z(2)
      dimension hx(2),hy(2),hz(2)

c
c     convert angle between -180 and 180 degrees

      angd = angd - anint(angd/360.0)*360.0

c
c     convert to radians
      
      ang = angd*(3.141592653589793d0/180.d0)
      ang1 = angt*(3.141592653589793d0/180.d0)
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

c
c     arbitary unit vector

      if(abs(xab).lt.0.5d0) then

         ax = 1.0d0
         ay = 0.d0
         az = 0.d0

      elseif(abs(yab).lt.0.5d0) then

         ax = 0.0d0
         ay = 1.d0
         az = 0.d0

      elseif(abs(zab).lt.0.5d0) then

         ax = 0.0d0
         ay = 0.d0
         az = 1.d0

      endif

c
c     axis vector orthogonal to b-a

      bx = yab*az - zab*ay
      by = zab*ax - xab*az
      bz = xab*ay - yab*ax

c
c     rotate b-a vector onto first h-a vector

      cs = cos(ang)
      sn = sin(ang)

      ahx = xab*cs - bx*sn
      ahy = yab*cs - by*sn
      ahz = zab*cs - bz*sn

c
c     second a-h  : rotate first a-h about b-a bond by angt degrees
c     but first find remaining axis vector:

      ax = yab*bz - zab*by
      ay = zab*bx - xab*bz
      az = xab*by - yab*bx

      cs = cos(ang1)
      sn = sin(ang1)
c
c     project ah onto ab axis

      dpab = ahx*xab + yab*ahy + zab*ahz
c
c     project ah onto b axis

      dpb = bx*ahx + by*ahy + bz*ahz
c
c     rotate in a-b plane

      cs = cs*dpb
      sn = sn*dpb

      vex = bx*cs - ax*sn
      vey = by*cs - ay*sn
      vez = bz*cs - az*sn
c
c     add back in ab component

      vex = vex + dpab*xab
      vey = vey + dpab*yab
      vez = vez + dpab*zab
c
c     hydrogen positions

      hx(1) = x(1) + ahx*bl
      hy(1) = y(1) + ahy*bl
      hz(1) = z(1) + ahz*bl

      hx(2) = x(1) + vex*bl
      hy(2) = y(1) + vey*bl
      hz(2) = z(1) + vez*bl

      return
      
      end

