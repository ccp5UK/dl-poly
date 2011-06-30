      program fold8
c
c**********************************************************************
c
c     dl_poly utility to expand a simulation cell by a factor of
c     eight in size using a 2x2x2 multiplication
c
c     copyright daresbury laboratory 1994
c     author w. smith september 1994
c
c     Note: under no circumstances should the velocity arrays be
c           treated using this program
c
c**********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (mega=1000000)
      character*8 name
      character*80 title
      character*40 fname
      dimension cell(9)

c
c     read control parameters

      write(*,'(a)')'enter name of config file'
      read(*,'(a40)')fname

c
c     open files

      open(7,file=fname)
      open(8,file='CONFIG_EIGHT')

c
c     read configuration data
      read(7,'(a)')title
      write(8,'(a)')title
      write(*,'(a)')'File header: ',title
      read(7,*)levcfg,imcon
      write(8,'(2i10)')0,imcon

      if(imcon.eq.0)then
         write(*,'(a)')'error - eightfold is for periodic cells only'
         stop
      endif

      read(7,*)cell(1),cell(2),cell(3)
      read(7,*)cell(4),cell(5),cell(6)
      read(7,*)cell(7),cell(8),cell(9)

      write(8,'(3f20.12)')2.d0*cell(1),2.d0*cell(2),2.d0*cell(3)
      write(8,'(3f20.12)')2.d0*cell(4),2.d0*cell(5),2.d0*cell(6)
      write(8,'(3f20.12)')2.d0*cell(7),2.d0*cell(8),2.d0*cell(9)

      do i=1,9
         cell(i)=0.5d0*cell(i)
      enddo

      do i=1,mega

         read(7,'(a)',end=100)name
         read(7,*,end=100)xx0,yy0,zz0
         if(levcfg.gt.0)read(7,*,end=100)vxx,vyy,vzz
         if(levcfg.gt.1)read(7,*,end=100)fxx,fyy,fzz

         xxx=xx0+cell(1)+cell(4)+cell(7)
         yyy=yy0+cell(2)+cell(5)+cell(8)
         zzz=zz0+cell(3)+cell(8)+cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0-cell(1)+cell(4)+cell(7)
         yyy=yy0-cell(2)+cell(5)+cell(8)
         zzz=zz0-cell(3)+cell(8)+cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0+cell(1)-cell(4)+cell(7)
         yyy=yy0+cell(2)-cell(5)+cell(8)
         zzz=zz0+cell(3)-cell(8)+cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0-cell(1)-cell(4)+cell(7)
         yyy=yy0-cell(2)-cell(5)+cell(8)
         zzz=zz0-cell(3)-cell(8)+cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0+cell(1)+cell(4)-cell(7)
         yyy=yy0+cell(2)+cell(5)-cell(8)
         zzz=zz0+cell(3)+cell(8)-cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0-cell(1)+cell(4)-cell(7)
         yyy=yy0-cell(2)+cell(5)-cell(8)
         zzz=zz0-cell(3)+cell(8)-cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0+cell(1)-cell(4)-cell(7)
         yyy=yy0+cell(2)-cell(5)-cell(8)
         zzz=zz0+cell(3)-cell(8)-cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

         xxx=xx0-cell(1)-cell(4)-cell(7)
         yyy=yy0-cell(2)-cell(5)-cell(8)
         zzz=zz0-cell(3)-cell(8)-cell(9)
         write(8,'(a8)')name
         write(8,'(3f20.12)')xxx,yyy,zzz

      enddo
  100 natm=i-1

      write(*,'(a,i6)')'number of atoms in new config file: ',natm*8

      close (7)
      close (8)

      stop
      end
