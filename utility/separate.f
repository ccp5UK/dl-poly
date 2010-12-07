c*********************************************************************
c
c     dl_poly utility to select two atoms in a system and adjust
c     their distances apart
c
c     author w smith june 2001
c     copyright daresbury laboratory
c
c*********************************************************************

      implicit real*8(a-h,o-z)
      parameter (natms=6124,nlist=100)
      character*80 title
      character*8 name(natms)
      dimension lst1(nlist),lst2(nlist)
      dimension x(natms),y(natms),z(natms)

      write(*,*)'Enter required C-C separation (A) (max 10-A)'
      read(*,*)ddd
      dd=min(dd,10.d0)
      write(*,*)'Enter IDs of atoms concerned'
      read(*,*)id1,id2
      write(*,*)'Number of atoms in molecule 1'
      read(*,*)num1
      write(*,*)'list atoms in molecule 1'
      read(*,*)(lst1(i),i=1,num1)
      write(*,*)'Number of atoms in molecule 2'
      read(*,*)num2
      write(*,*)'list atoms in molecule 2'
      read(*,*)(lst2(i),i=1,num2)

c     open config files

      open(7,file='CFG0')
      open(8,file='CFG1')

c     read config file CFG0

      read(7,'(a80)')title
      write(8,'(a80)')title
      read(7,'(2i10)')i,j
      write(8,'(2i10)')i,j
      read(7,'(3f20.0)')ax,ay,az
      write(8,'(3f20.12)')ax,ay,az
      read(7,'(3f20.0)')bx,by,bz
      write(8,'(3f20.12)')bx,by,bz
      read(7,'(3f20.0)')cx,cy,cz
      write(8,'(3f20.12)')cx,cy,cz
      do i=1,natms

         read(7,'(a8)')name(i)
         read(7,'(3f20.0)')x(i),y(i),z(i)

      enddo

c     calculate separation vector

      dx=x(id2)-x(id1)-ax*nint((x(id2)-x(id1))/ax)
      dy=y(id2)-y(id1)-by*nint((y(id2)-y(id1))/by)
      dz=z(id2)-z(id1)-cz*nint((z(id2)-z(id1))/cz)
      rrr=sqrt(dx*dx+dy*dy+dz*dz)
      dx=dx/rrr
      dy=dy/rrr
      dz=dz/rrr
      dis=0.5d0*(ddd-rrr)

c     shift molecule 1

      do i=1,num1

         id=lst1(i)
         x(id)=x(id)-dis*dx
         y(id)=y(id)-dis*dy
         z(id)=z(id)-dis*dz
         x(id)=x(id)-ax*nint(x(id)/ax)
         y(id)=y(id)-by*nint(y(id)/by)
         z(id)=z(id)-cz*nint(z(id)/cz)

      enddo

c     shift molecule 2

      do i=1,num2

         id=lst2(i)
         x(id)=x(id)+dis*dx
         y(id)=y(id)+dis*dy
         z(id)=z(id)+dis*dz
         x(id)=x(id)-ax*nint(x(id)/ax)
         y(id)=y(id)-by*nint(y(id)/by)
         z(id)=z(id)-cz*nint(z(id)/cz)
         
      enddo

c     write new config file

      do i=1,natms

         write(8,'(a8)')name(i)
         write(8,'(3f20.12)')x(i),y(i),z(i)

      enddo

      end
