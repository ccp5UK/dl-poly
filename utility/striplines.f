c*********************************************************************
c
c     dl_poly utility to strip duplicate lines from a file
c
c     copyright daresbury laboratory
c     author w.smith dec 2006
c
c*********************************************************************

      logical same
      character*1 oldline(80),newline(80)
      character*40 file1,file2

      write(*,*)'enter file name:'
      read(*,*)file1
      file2="STRIPPED"
      open(7,file=file1)
      open(8,file=file2)

      read(7,'(80a1)',end=100)oldline
      write(8,'(80a1)')oldline

      do while(.true.)

        read(7,'(80a1)',end=100)newline

        same=.true.
        do i=1,80
          if(oldline(i).ne.newline(i))same=.false.
        enddo
        if(.not.same)write(8,'(80a1)')newline

        do i=1,80
          oldline(i)=newline(i)
        enddo

      enddo

  100 close (7)
      close (8)
      end
