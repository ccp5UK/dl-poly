      program cfgcmp

c*********************************************************************
c
c     dl_poly utility program for comparing two CONFIG files
c
c     copyright daresbury laboratory 1998
c     author w.smith december 1998
c
c*********************************************************************

      implicit real*8 (a-h,o-z)
      parameter (mxatms=10000)

      character*8 name
      character*40 file1,file2
      real*8 cel0(9),cel1(9)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension vx0(mxatms),vy0(mxatms),vz0(mxatms)
      dimension fx0(mxatms),fy0(mxatms),fz0(mxatms)
      dimension xx1(mxatms),yy1(mxatms),zz1(mxatms)
      dimension vx1(mxatms),vy1(mxatms),vz1(mxatms)
      dimension fx1(mxatms),fy1(mxatms),fz1(mxatms)

      write(*,*)'Enter name of first file'
      read(*,*)file1
      write(*,*)'Enter name of second file'
      read(*,*)file2
      open(7,file=file1)
      open(8,file=file2)

      read(7,*)
      read(7,*)lcfg0,imcn0
      if(imcn0.gt.0)then
        read(7,*)cel0(1),cel0(2),cel0(3)
        read(7,*)cel0(4),cel0(5),cel0(6)
        read(7,*)cel0(7),cel0(8),cel0(9)
      endif
      read(8,*)
      read(8,*)lcfg1,imcn1
      if(imcn1.gt.0)then
        read(8,*)cel1(1),cel1(2),cel1(3)
        read(8,*)cel1(4),cel1(5),cel1(6)
        read(8,*)cel1(7),cel1(8),cel1(9)
      endif
      if(imcn0.ne.imcn1)then
        write(*,*)'error - image convention different in CONFIG files'
        call exit(0)
      endif
      celmx=0.d0

      do i=1,9
        celmx=max(celmx,abs(cel0(i)-cel1(i)))
      enddo

      do i=1,mxatms
        read(7,*,end=100)name
        read(7,*)xx0(i),yy0(i),zz0(i)
        if(lcfg0.gt.0)read(7,*)vx0(i),vy0(i),vz0(i)
        if(lcfg0.gt.1)read(7,*)fx0(i),fy0(i),fz0(i)
      enddo
  100 natms=i-1
      do i=1,mxatms
        read(8,*,end=200)name
        read(8,*)xx1(i),yy1(i),zz1(i)
        if(lcfg1.gt.0)read(8,*)vx1(i),vy1(i),vz1(i)
        if(lcfg1.gt.1)read(8,*)fx1(i),fy1(i),fz1(i)
      enddo
  200 matms=i-1
      if(natms.ne.matms)then
        write(*,*)'error - CONFIG files have different atom numbers'
        call exit(0)
      endif

c     periodic boundary correction

      do i=1,natms
        xx0(i)=xx1(i)-xx0(i)
        yy0(i)=yy1(i)-yy0(i)
        zz0(i)=zz1(i)-zz0(i)
      enddo

      call images(imcn0,0,1,natms,cel0,xx0,yy0,zz0)

      erx=0.d0
      ery=0.d0
      erz=0.d0
      evx=0.d0
      evy=0.d0
      evz=0.d0
      efx=0.d0
      efy=0.d0
      efz=0.d0

      do i=1,mxatms

        if(erx.lt.abs(xx0(i)))then
          erx=abs(xx0(i))
          k1=i
        endif
        if(ery.lt.abs(yy0(i)))then
          ery=abs(yy0(i))
          k2=i
        endif
        if(erz.lt.abs(zz0(i)))then
          erz=abs(zz0(i))
          k3=i
        endif
        if(lcfg0.gt.0.and.lcfg1.gt.0)then
          if(evx.lt.abs(vx0(i)-vx1(i)))then
            evx=abs(vx0(i)-vx1(i))
            k4=i
          endif
          if(evy.lt.abs(vy0(i)-vy1(i)))then
            evy=abs(vy0(i)-vy1(i))
            k5=i
          endif
          if(evz.lt.abs(vz0(i)-vz1(i)))then
            evz=abs(vz0(i)-vz1(i))
            k6=i
          endif
        endif
        if(lcfg0.gt.1.and.lcfg1.gt.1)then
          if(efx.lt.abs(fx0(i)-fx1(i)))then
            efx=abs(fx0(i)-fx1(i))
            k7=i
          endif
          if(efy.lt.abs(fy0(i)-fy1(i)))then
            efy=abs(fy0(i)-fy1(i))
            k8=i
          endif
          if(efz.lt.abs(fz0(i)-fz1(i)))then
            efz=abs(fz0(i)-fz1(i))
            k9=i
          endif
        endif

      enddo
      write(*,*)'Maximum error in cell vectors'
      write(*,'(1p,e12.4)')celmx
      write(*,*)'Maximum errors in position'
      write(*,'(3i5,1p,3e12.4)')k1,k2,k3,erx,ery,erz
      if(lcfg0.gt.0.and.lcfg1.gt.0)then
        write(*,*)'Maximum errors in velocity'
        write(*,'(3i5,1p,3e12.4)')k4,k5,k6,evx,evy,evz
      endif
      if(lcfg0.gt.1.and.lcfg1.gt.1)then
        write(*,*)'Maximum errors in force'
        write(*,'(3i5,1p,3e12.4)')k7,k8,k9,efx,efy,efz
      endif

      end

      subroutine images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)

c***********************************************************************
c
c     dl_poly subroutine for calculating the minimum image
c     of atom pairs within a specified MD cell
c
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     T3D optimised version. t.forester july 1994
c
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge
c     coordinate arrays
c
c***********************************************************************

      implicit none

      integer imcon,idnode,mxnode,natms,iatm1,iatm2,i
      real*8 cell,xxx,yyy,zzz,aaa,bbb,ccc,det,rt2,rt3,ssx
      real*8 ssy,ssz,ddd,xss,yss,zss,rcell

      dimension xxx(*),yyy(*),zzz(*)
      dimension cell(9),rcell(9)

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/

      if(imcon.gt.0) then

c     block indices

        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif

      if(imcon.eq.1)then

c     standard cubic boundary conditions


        aaa=1.d0/cell(1)

        do i=iatm1,iatm2
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
        enddo

      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions

        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)

        do i=iatm1,iatm2

          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))

        enddo

      else if(imcon.eq.3)then

c     parallelepiped boundary conditions

        call invert(cell,rcell,det)

        do i=iatm1,iatm2

          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))

          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)

          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

        enddo

      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions

        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(5)-cell(9)).lt.1.d-6)) then
          write(*,*)'error - TO boundaries incorrectly defined'
          call exit(0)
        endif

        aaa=1.d0/cell(1)

        do i=iatm1,iatm2

          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))

          if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.
     x      (0.75d0*cell(1)))then

            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(1),zzz(i))

          endif

        enddo

      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions

        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) then
          write(*,*)'error - RD boundaries incorrectly defined'
          call exit(0)
        endif


        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)

        do i=iatm1,iatm2

          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(bbb*zzz(i))

          if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge.
     x      cell(1))then

            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(9),zzz(i))

          endif

        enddo

      else if(imcon.eq.6) then

c     x-y boundary conditions

        det = cell(1)*cell(5) - cell(2)*cell(4)

        if(abs(det).lt.1.d-6)then
          write(*,*)'error - slab boundaries incorrectly defined'
          call exit(0)
        endif

        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)

        do i=iatm1,iatm2

          ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
          ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

          xss = ssx - nint(ssx)
          yss = ssy - nint(ssy)

          xxx(i)=cell(1)*xss + cell(4)*yss
          yyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions

        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)then
          write(*,*)'error - hexagonal boundaries incorrectly defined'
          call exit(0)
        endif

        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)

        do i=iatm1,iatm2

          yyy(i)=yyy(i)-bbb*nint(ccc*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ddd*zzz(i))

          if((abs(yyy(i))+abs(rt3*xxx(i))).ge.bbb)then

            xxx(i)=xxx(i)-rt3*sign(aaa,xxx(i))
            yyy(i)=yyy(i)-sign(aaa,yyy(i))

          endif

        enddo

      endif

      return
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
