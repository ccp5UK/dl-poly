      program nfold

c**********************************************************************
c
c     dl_poly utility to expand a simulation cell by a
c     nx*ny*nz multiplication
c
c     copyright daresbury laboratory 1997
c     author w. smith june 1997
c
c     Note: under no circumstances should the velocity arrays be
c           treated using this program
c
c**********************************************************************

      implicit real*8(a-h,o-z)

      character*8 name
      character*80 title
      character*40 fname
      dimension cell(9),cprp(10)

c     read control parameters

      write(*,'(a)')'enter name of config file'
      read(*,'(a40)')fname
      write(*,'(a)')'enter nx,ny,nz'
      read(*,*)nx,ny,nz

c     open output file

      open(8,file='CFGBIG')

c     read configuration file headers

      open(7,file=fname)
      read(7,'(a)')title
      write(8,'(a)')title
      write(*,'(a)')'File header: ',title
      read(7,*)levcfg,imcon
      write(8,'(2i10)')0,imcon

      if(imcon.eq.0)then
         write(*,'(a)')'error - nfold is for periodic cells only'
         stop
      endif

      fx=dble(nx)
      fy=dble(ny)
      fz=dble(nz)
      read(7,*)cell(1),cell(2),cell(3)
      read(7,*)cell(4),cell(5),cell(6)
      read(7,*)cell(7),cell(8),cell(9)

      write(8,'(3f20.12)')fx*cell(1),fx*cell(2),fx*cell(3)
      write(8,'(3f20.12)')fy*cell(4),fy*cell(5),fy*cell(6)
      write(8,'(3f20.12)')fz*cell(7),fz*cell(8),fz*cell(9)

c     scan contents of cell

      n=0
      natms=0
      do while(.true.)

        read(7,'(a)',end=100)name
        read(7,*,end=100)xx0,yy0,zz0
        if(levcfg.gt.0)read(7,*,end=100)vxx,vyy,vzz
        if(levcfg.gt.1)read(7,*,end=100)fxx,fyy,fzz

c     create replicas of each atomic coordinate

        do  k=1,nz

          f7=cell(7)*dble(2*k-nz-1)/2.d0
          f8=cell(8)*dble(2*k-nz-1)/2.d0
          f9=cell(9)*dble(2*k-nz-1)/2.d0

          do  j=1,ny

            f4=cell(4)*dble(2*j-ny-1)/2.d0
            f5=cell(5)*dble(2*j-ny-1)/2.d0
            f6=cell(6)*dble(2*j-ny-1)/2.d0

            do  i=1,nx

              f1=cell(1)*dble(2*i-nx-1)/2.d0
              f2=cell(2)*dble(2*i-nx-1)/2.d0
              f3=cell(3)*dble(2*i-nx-1)/2.d0

              n=n+1

              xx=xx0+f1+f4+f7
              yy=yy0+f2+f5+f8
              zz=zz0+f3+f6+f9
              write(8,'(a8,i10,/,3f20.12)')name,n,xx,yy,zz

            enddo
          enddo
        enddo

        natms=natms+1

      enddo

  100 close(7)
      close (8)

c     write summary data

      call dcell(cell,cprp)
      wdth=0.5d0*min(dble(nx)*cprp(7),dble(ny)*cprp(8),dble(nz)*cprp(9))
      write(*,*)'number of ions in system = ',natms*nx*ny*nz
      write(*,*)'maximum radius of cutoff = ',wdth
      write(*,*)'output file name: CFGBIG'

      stop
      end

      subroutine dcell(aaa,bbb)

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


      implicit real*8 (a-h,o-z)

      dimension aaa(9),bbb(10)

c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))

c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))

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

c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)

c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end
