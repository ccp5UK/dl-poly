      program statsum

c---------------------------------------------------------------------
c
c     dl_poly utility program to summarise the contents of a
c     dl_poly STATIS file
c
c     copyright daresbury laboratory
c     author w.smith september 2007
c
c---------------------------------------------------------------------

      implicit none

      logical fixvolume
      character*80 record
      integer i,step,nitems,count,ntypes
      real(8) time,fac1,fac2
      real(8), allocatable :: data(:),average(:),variance(:)

c     open the STATIS file

      open(7,file='STATIS')

c     write headers

      read(7,'(a80)')record
      write(*,'(a80)')record
      read(7,'(a80)')record
      write(*,'(a80)')record

c     read file contents

      count=0
      do while(.true.)

        read(7,*,end=100)step,time,nitems

c     allocate data arrays

        if(count.eq.0)then

          allocate(data(nitems),average(nitems),variance(nitems))
          do i=1,nitems
            average(i)=0.d0
            variance(i)=0.d0
          enddo

        endif

c     update counter

        count=count+1
        fac2=1.d0/dble(count)
        fac1=dble(count-1)/dble(count)

c     read data records

        read(7,'(5e14.6)')data

c     accumulate average and variance data

        do i=1,nitems

          variance(i)=fac1*(variance(i)+fac2*(data(i)-average(i))**2)
          average(i)=fac1*average(i)+fac2*data(i)

        enddo

      enddo

c     close STATIS file

 100  close(7)
      write(*,*)'number of data records read: ',count

c     summarise average data

      do i=1,nitems
        variance(i)=sqrt(variance(i))
      enddo

c     check if system has a constant volume

      fixvolume=(variance(19).lt.1.d-8)
      if(fixvolume)then
        ntypes=nitems-36
      else
        ntypes=nitems-45
      endif

c     write out calculated properties

 10   format(1x,a20,1p,2e12.4)
      write(*,10)'system temperature  ',average(2),variance(2)
      write(*,10)'system pressure     ',average(27),variance(27)
      write(*,10)'system volume       ',average(19),variance(19)
      write(*,10)'system energy       ',average(1),variance(1)
      write(*,10)'configuration energy',average(3),variance(3)
      write(*,10)'short range energy  ',average(4),variance(4)
      write(*,10)'coulombic energy    ',average(5),variance(5)
      write(*,10)'bond energy         ',average(6),variance(6)
      write(*,10)'bond angle energy   ',average(7),variance(7)
      write(*,10)'dihedral energy     ',average(8),variance(8)
      write(*,10)'tethering energy    ',average(9),variance(9)
      write(*,10)'system enthalpy     ',average(10),variance(10)
      write(*,10)'rotational energy   ',average(11),variance(11)
      write(*,10)'system virial       ',average(12),variance(12)
      write(*,10)'short range virial  ',average(13),variance(13)
      write(*,10)'coulombic virial    ',average(14),variance(14)
      write(*,10)'bond virial         ',average(15),variance(15)
      write(*,10)'three body virial   ',average(16),variance(16)
      write(*,10)'constraint virial   ',average(17),variance(17)
      write(*,10)'tethering virial    ',average(18),variance(18)
      write(*,10)'constraint virial   ',average(26),variance(26)
      write(*,10)'shell tmperature    ',average(20),variance(20)
      write(*,10)'shell energy        ',average(21),variance(21)
      write(*,10)'shell virial        ',average(22),variance(22)
      write(*,10)'cell angle (alpha)  ',average(23),variance(23)
      write(*,10)'cell angle (beta)   ',average(24),variance(24)
      write(*,10)'cell angle (gamma)  ',average(25),variance(25)
      i=28+ntypes
      write(*,10)'pressure tensor (xx)',average(i),variance(i)
      write(*,10)'pressure tensor (xy)',average(i+1),variance(i+1)
      write(*,10)'pressure tensor (xz)',average(i+2),variance(i+2)
      write(*,10)'pressure tensor (yy)',average(i+4),variance(i+4)
      write(*,10)'pressure tensor (yz)',average(i+5),variance(i+5)
      write(*,10)'pressure tensor (zz)',average(i+8),variance(i+8)
      if(.not.fixvolume)then
        i=i+9
        write(*,10)'cell vector     (Ax)',average(i),variance(i)
        write(*,10)'cell vector     (Ay)',average(i+1),variance(i+1)
        write(*,10)'cell vector     (Az)',average(i+2),variance(i+2)
        write(*,10)'cell vector     (Bx)',average(i+3),variance(i+3)
        write(*,10)'cell vector     (By)',average(i+4),variance(i+4)
        write(*,10)'cell vector     (Bz)',average(i+5),variance(i+5)
        write(*,10)'cell vector     (Cx)',average(i+6),variance(i+6)
        write(*,10)'cell vector     (Cy)',average(i+7),variance(i+7)
        write(*,10)'cell vector     (Cz)',average(i+8),variance(i+8)
      endif
 20   format(1x,a16,i4,1p,2e12.4)
      do i=28,ntypes+27
        write(*,20)'MSD of Atom Type',i-27,average(i),variance(i)
      enddo

      end


