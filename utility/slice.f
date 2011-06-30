c      program slice

c*********************************************************************
c
c     dl_poly program to extract a slice of atoms from the centre
c     in a given direction from a dl_poly CONFIG file
c
c     copyright daresbury laboratory 1997
c     author  w.smith june 1997
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      character*40 fname
      character*80 title
      character*8 atname
      dimension cell(9)

c     open the I/O files

c      open(5,file='slice_input')
c      open(6,file='slice_output')

      write(6,'(a)')'# X-density Program'

c     name of selected CONFIG file
      write(6,'(a)')'Enter name of CONFIG file'
      read(5,*)fname

c     total number of atoms in a CONFIG file configuration
      write(6,'(a)')'Enter number of atoms in MD cell'
      read(5,*)natms

c     bin width for slice
      write(6,'(a)')'Enter slice width (A)'
      read(5,*)delx

c     projection vector
      write(6,'(a)')'Enter direction vector'
      read(5,*)xx0,yy0,zz0
      rrr=sqrt(xx0**2+yy0**2+zz0**2)
      xx0=xx0/rrr
      yy0=yy0/rrr
      zz0=zz0/rrr

c     check on specified control variables

      write(6,'(a,a40)')'# name of target CONFIG file   : ',fname
      write(6,'(a,i8)')'# total no. of atoms in config  : ',natms
      write(6,'(a,1p,e12.4)')'# x-density histogram bin width : ',delx
      write(6,'(a,1p,3e12.4)')'# normalised projection vector  : ',
     x    xx0,yy0,zz0

c     open CONFIG file

      open(7,file=fname)
      open(8,file='SLICE')

      write(8,'(i5)')natms

      read(7,'(a80)')title
      write(6,'(a10,a80)')'# Header: ',title
      write(8,'(a80)')title
      read(7,*)levcfg,imcon

      if(imcon.gt.0)then

        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)

        write(6,'(a)')     '# cell vectors:'
        write(6,'(a,3f12.6)')'# vector A :',cell(1),cell(2),cell(3)
        write(6,'(a,3f12.6)')'# vector B :',cell(4),cell(5),cell(6)
        write(6,'(a,3f12.6)')'# vector C :',cell(7),cell(8),cell(9)

      endif

c     read configuration

      do i=1,natms

        read(7,'(a8)')atname

        read(7,*)xxx,yyy,zzz
        if(levcfg.gt.0)read(7,*)uuu,vvv,www
        if(levcfg.gt.1)read(7,*)uuu,vvv,www

        rrr=xxx*xx0+yyy*yy0+zzz*zz0
        ix=nint(rrr/delx)+ndiv/2
        if(ix.eq.0)then

          write(8,'(a8,1p,3e15.7)')atname,xxx,yyy,zzz

        endif

      enddo

      close (7)
      close (8)

      end
