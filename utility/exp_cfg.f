c      program sclexp

c*********************************************************************
c     
c     dl_poly program to scale the volume of a CONFIG file
c     
c     copyright daresbury laboratory 2000
c     author  w.smith june 2000
c     
c     wl
c     1996/02/15 14:33:26
c     1.1.1.1
c     Exp
c
c*********************************************************************
      
      implicit real*8(a-h,o-z)
      
      parameter (mxatms=999999)

      character*40 fname
      character*80 title
      character*8 atname
      dimension cell(9)
        
c     open the I/O files
      
c      open(5,file='scale_input')
c      open(6,file='scale_output')
      
      write(*,'(a)')'# Scale CONFIG Program'
      
c     name of selected CONFIG file
      write(*,'(a)')'Enter name of CONFIG file'
      read(*,*)fname

c     required volume scale factor

      write(*,'(a)')'Enter required volume scale factor'
      read(*,*)scale

c     check on specified control variables
      
      write(*,'(a,a40)')    '# name of target CONFIG file   : ',fname
      write(*,'(a,1p,e12.4)')'# requred volume scale factor  : ',scale
      scale=scale**(1.d0/3.d0)
        
c     open CONFIG file

      open(7,file=fname)
      open(8,file='CFGVOL')

      read(7,'(a80)')title
      write(*,'(a10,a80)')'# Header: ',title
      write(8,'(a80)')title
      read(7,*)levcfg,imcon
      write(8,'(2i10)')levcfg,imcon

      if(imcon.gt.0)then

        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)

        write(8,'(3f20.10)')scale*cell(1),scale*cell(2),scale*cell(3)
        write(8,'(3f20.10)')scale*cell(4),scale*cell(5),scale*cell(6)
        write(8,'(3f20.10)')scale*cell(7),scale*cell(8),scale*cell(9)

      endif

c     read configuration

      do i=1,mxatms

        read(7,'(a8)',end=100)atname
        write(8,'(a8)')atname

        read(7,*)xxx,yyy,zzz
        write(8,'(3f20.10)')scale*xxx,scale*yyy,scale*zzz
        if(levcfg.gt.0)then

          read(7,*)uuu,vvv,www
          write(8,'(3f20.10)')uuu,vvv,www
          if(levcfg.gt.1)then

            read(7,*)uuu,vvv,www
            write(8,'(3f20.10)')uuu,vvv,www

          endif

        endif

      enddo

  100 continue

      close (7)
      close (8)
      write(*,'(a)')'job done'
      stop
      
      end
