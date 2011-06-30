      program resize
c
c***********************************************************************
c
c     dl_poly utility program to resize a config file by linearly
c     scaling all the dimensions
c
c     author - w.smith december 1994
c     copyright daresbury laboratory 1994
c
c***********************************************************************
c
      parameter (mega=1000000)
      implicit real*8(a-h,o-z)

      character*80 header,name,fname
      dimension cell(9)

      write(*,'(a)')'enter the config file name'
      read(*,*)fname
      write(*,'(a)')'enter the scaling factor'
      read(*,*)scale
      write(*,*)'output file will be called NEW_CONFIG'

c
c     open the config files

      open (10,file=fname)
      open (11,file='NEW_CONFIG')

c
c     read the CONFIG file header

      read(10,'(a80)',end=100)header
      write(11,'(a80)')header
      read(10,'(2i10)',end=100) levcfg,imcon
      write(11,'(2i10)') levcfg,imcon

c
c     specify molecular dynamics simulation cell

      if(imcon.gt.0)then

         read(10,*,end=100)cell(1),cell(2),cell(3)
         read(10,*,end=100)cell(4),cell(5),cell(6)
         read(10,*,end=100)cell(7),cell(8),cell(9)
         do i=1,9

            cell(i)=scale*cell(i)

         enddo
         write(11,'(3f20.15)')cell(1),cell(2),cell(3)
         write(11,'(3f20.15)')cell(4),cell(5),cell(6)
         write(11,'(3f20.15)')cell(7),cell(8),cell(9)

      endif

      do i=1,mega

         if(levcfg.eq.0)then

            read(10,'(a80)',end=100) name
            write(11,'(a80)') name
            read(10,*)xxx,yyy,zzz
            write(11,'(3g20.10)')scale*xxx,scale*yyy,scale*zzz

         else if(levcfg.eq.1)then

            read(10,'(a80)',end=100) name
            write(11,'(a80)') name
            read(10,*)xxx,yyy,zzz
            write(11,'(3g20.10)')scale*xxx,scale*yyy,scale*zzz
            read(10,*)xxx,yyy,zzz
            write(11,'(3g20.10)')xxx,yyy,zzz

         else

            read(10,'(a80)',end=100) name
            write(11,'(a80)') name
            read(10,*)xxx,yyy,zzz
            write(11,'(3g20.10)')scale*xxx,scale*yyy,scale*zzz
            read(10,*)xxx,yyy,zzz
            write(11,'(3g20.10)')xxx,yyy,zzz
            read(10,*)xxx,yyy,zzz
            write(11,'(3g20.10)')xxx,yyy,zzz

         endif

      enddo

  100 continue

      write(*,'(a,i6)')'number of atoms processed = ',i-1

      close (10)
      close (11)

      end
