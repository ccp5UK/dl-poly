      subroutine cdens
     x  (idnode,mxnode,imcon,natms,nclust,weight,csep,ccut,cell,
     x  xxx,yyy,zzz,xdf,ydf,zdf,crho)

c*********************************************************************
c
c     dl_poly subroutine to calculate cylindrical density about
c     the centre-centre axis of two clusters a fixed distance apart
c
c     copyright daresbury laboratory
c     author w.smith jan 2005
c
c*********************************************************************

      include "dl_params.inc"

      integer idnode,mxnode,imcon,natms,i,last
      real*8 csep,ccut,mass1,mass2,vnm,rrr,rzz,drad,dzed
      integer nclust(4)
      real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
      real*8 xdf(mxatms),ydf(mxatms),zdf(mxatms)
      real*8 crho(mxzed,mxrad),cell(9),com(3),com1(3),com2(3),dir(3)

c     density array bin widths

      drad=ccut/mxrad
      dzed=csep/mxzed

c     calculate cluster centres of mass

      do i=1,3

        com1(i)=0.d0
        com2(i)=0.d0

      enddo

      mass1=0.d0
      do i=nclust(1),nclust(2)

        mass1=mass1+weight(i)
        com1(1)=com1(1)+weight(i)*xxx(i)
        com1(2)=com1(2)+weight(i)*yyy(i)
        com1(3)=com1(3)+weight(i)*zzz(i)

      enddo
      com1(1)=com1(1)/mass1
      com1(2)=com1(2)/mass1
      com1(3)=com1(3)/mass1

      mass2=0.d0
      do i=nclust(3),nclust(4)

        mass2=mass2+weight(i)
        com2(1)=com2(1)+weight(i)*xxx(i)
        com2(2)=com2(2)+weight(i)*yyy(i)
        com2(3)=com2(3)+weight(i)*zzz(i)

      enddo
      com2(1)=com2(1)/mass2
      com2(2)=com2(2)/mass2
      com2(3)=com2(3)/mass2

c     cluster pair centre of mass

      com(1)=0.5d0*(com1(1)+com2(1))
      com(2)=0.5d0*(com1(2)+com2(2))
      com(3)=0.5d0*(com1(3)+com2(3))

c     vector between cluster centres of mass

      dir(1)=com1(1)-com2(1)
      dir(2)=com1(2)-com2(2)
      dir(3)=com1(3)-com2(3)

      call images(imcon,0,1,1,cell,dir(1),dir(2)dir(3))

      vnm=1.d0/sqrt(dir(1)**2+dir(2)**2+dir(3)**3)

c     calculate displacement vectors

      j=0
      do i=1,natms

        if(.not.((i.ge.nclust(1).and.i.le.nclust(2)).or.
     x    (i.ge.nclust(3).and.i.le.nclust(4))))then

          j=j+1
          xdf(j)=xxx(i)-com(1)
          ydf(j)=yyy(i)-com(2)
          zdf(j)=zzz(i)-com(3)

        endif

      enddo

      call images(imcon,0,1,j,cell,xdf,ydf,zdf)
      last=j

c     accumulate cylindrical density

      do i=1,last

        rzz=vnm*(dir(1)*xdf(i)+dir(2)*ydf(i)+dir(3)*zdf(i))
        rrr=sqrt(rzz**2-(xdf(i)**2+ydf(i)**2+zdf(i)**2))
        irad=int(rrr/drad)+1
        ized=int(rzz/dzed)+1
        crho(ized,irad)=crho(ized,irad)+1.d0

      enddo

      return
      end

      subroutine ccom(idnode,mxnode,ia,ib,com,weight,xxx,yyy,zzz)

c*********************************************************************
c
c     dl_poly subroutine to calculate centre of mass
c
c     copyright daresbury laboratory
c     author w.smith jan 2005
c
c*********************************************************************

      implict none

      integer idnode,mxnode,ia,ib,i
      real*8 mass
      real*8 weight(*),xxx(*),yyy(*),zzz(*),com(3)

      mass=0.d0
      com(1)=0.d0
      com(2)=0.d0
      com(3)=0.d0

      do i=ia,ib

        mass=mass+weight(i)
        com(1)=com(1)+weight(i)*xxx(i)
        com(2)=com(2)+weight(i)*yyy(i)
        com(3)=com(3)+weight(i)*zzz(i)

      enddo

      com(1)=com(1)/mass
      com(2)=com(2)/mass
      com(3)=com(3)/mass

      return
      end
