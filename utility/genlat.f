      program genlat
c
c***********************************************************************
c
c     program to generate a perfect lattice with general unit cell
c
c     copyright daresbury laboratory
c     author - w.smith oct. 1992
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter (mxatom=1000)
      dimension cell(9),cprp(10),xb(mxatom),yb(mxatom),zb(mxatom)
      character*8 name(mxatom)
      character*80 title
      write(*,*)'enter suitable one-line title for lattice file'
      read(*,'(a80)')title
      write(*,*)'enter number of basis atoms (max 1000)'
      read(*,*)natms
      write(*,*)'define unit cell vector A(1-3)'
      read(*,*)cell(1),cell(2),cell(3)
      write(*,*)'define unit cell vector B(1-3)'
      read(*,*)cell(4),cell(5),cell(6)
      write(*,*)'define unit cell vector C(1-3)'
      read(*,*)cell(7),cell(8),cell(9)
      write(*,*)'define 3 unit cell multiplications LxA,MxB,NxC'
      read(*,*)nx,ny,nz
      write(*,*)'enter ion names (8 chars) and fractional coordinates'
      do i=1,natms
         read(*,*)name(i),xb(i),yb(i),zb(i)
      enddo
      open(10,file='LATTICE')
      write(10,'(a80)')title
      write(10,'(2i10)')0,3
      write(10,'(3f20.8)')dble(nx)*cell(1),dble(nx)*cell(2),
     x                    dble(nx)*cell(3),dble(ny)*cell(4),
     x                    dble(ny)*cell(5),dble(ny)*cell(6),
     x                    dble(nz)*cell(7),dble(nz)*cell(8),
     x                    dble(nz)*cell(9)
c
c     set up lattice
      m=0
      do  k=1,nz
         do  j=1,ny
            do  i=1,nx
               do  n=1,natms
                  m=m+1
                  xs=dble(i-1)+xb(n)-0.5d0*dble(nx)
                  ys=dble(j-1)+yb(n)-0.5d0*dble(ny)
                  zs=dble(k-1)+zb(n)-0.5d0*dble(nz)
                  xx=cell(1)*xs+cell(4)*ys+cell(7)*zs
                  yy=cell(2)*xs+cell(5)*ys+cell(8)*zs
                  zz=cell(3)*xs+cell(6)*ys+cell(9)*zs
                  write(10,'(a8,i10,/,3f20.6)')name(n),m,xx,yy,zz
               enddo
            enddo
         enddo
      enddo
      call dcell(cell,cprp)
      wdth=0.5d0*min(dble(nx)*cprp(7),dble(ny)*cprp(8),dble(nz)*cprp(9))
      write(*,*)'number of ions in system = ',m
      write(*,*)'maximum radius of cutoff = ',wdth
      stop
      end
      subroutine dcell(aaa,bbb)

c
c***********************************************************************
c
c     dl_poly subroutine to calculate the dimensional properies of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c
c     copyright daresbury laboratory
c     author - w. smith july 1992
c
c***********************************************************************
c

      implicit real*8 (a-h,o-z)

      dimension aaa(9),bbb(10)
c
c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))
c
c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))
c
c     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)
c
c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)
c
c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end
