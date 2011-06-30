      subroutine hadd1_1(x,y,z,hx,hy,hz,bl,angd)

c***********************************************************************
c
c     subroutine to add 1 terminal hydrogen
c              b-a...h
c
c     copyright daresbury laboratory 1996
c     author -      t forester   feb 1996
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      dimension x(*),y(*),z(*)
      dimension hx(*),hy(*),hz(*)
c
c     convert angle between -180 and 180 degrees

      angd = angd - anint(angd/360.d0)*360.d0

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

c     construct arbitary vector

      ohx = 1.0d0
      ohy = 1.d-2
      ohz = 1.d-2

c     make sure it is not too close to c-o vector

      if(abs(xab).gt. 0.9d0) then

         ohx = 1.d-2
         ohy = 1.0d0

      endif

c     construct vector orthogonal to both

      a1 = yab*ohz - zab*ohy
      a2 = zab*ohx - xab*ohz
      a3 = xab*ohy - yab*ohx

c     normalise a

      raa = 1.0d0/sqrt(a1*a1 + a2*a2 + a3*a3)
      a1 = a1*raa
      a2 = a2*raa
      a3 = a3*raa

c     rotate co vector in (co ,a) plane by angd

      ang = angd*3.141592659d0/180.d0

      cs = cos(ang)
      sn = sin(ang)

      ahx = cs*xab - sn*a1
      ahy = cs*yab - sn*a2
      ahz = cs*zab - sn*a3
c
c
      ahx = ahx*bl
      ahy = ahy*bl
      ahz = ahz*bl

c     get H positions

      hx(1) = x(1) + ahx
      hy(1) = y(1) + ahy
      hz(1) = z(1) + ahz

      return

      end



