      program tadlook

c***********************************************************************
c
c     dlpoly utility program for inspecting TADREV files
c
c     copyright daresbury laboratory
c     author      w.smith  dec  2007
c
c**********************************************************************

      implicit none

      character*40 filename
      integer natms,nsteql,i,j,netdif,numbsn,numtrk,numpro
      integer numdark,last
      real(8) dnumbsn,dnumtrk,dnumpro,dnetdif,dnumdark,dnsteql
      real(8) timtad,timres,tstop,tkeres
      real(8) strres(9),celres(9)
      real(8), allocatable, buffer(:)

      write(*,*)'enter name of TADREV file'
      read(*,'(a40)')filename
      write(*,*)'enter number of atoms'
      read(*,*)natms

c     open TADREV file

      open(7,file=filename,form="unformatted")

c     write control variables

      read(7)dnumbsn,dnumtrk,dnumpro,dnetdif,dnumdark,
     x  dnsteql,timtad,timres,tstop,tkeres,strres,celres

      numbsn=nint(dnumbsn)
      numtrk=nint(dnumtrk)
      numpro=nint(dnumpro)
      netdif=nint(dnetdif)
      numdark=nint(dnumdark)
      nsteql=nint(dnsteql)

      write(*,*)'number of basins          ',numbsn
      write(*,*)'track file number         ',numtrk
      write(*,*)'number of profiles        ',numpro
      write(*,*)'number of diffs           ',netdif
      write(*,*)'end of dark period        ',numdark
      write(*,*)'equilibration period      ',nsteql
      write(*,*)'hyperdynamics time        ',timtad
      write(*,*)'time of rewind file       ',timres
      write(*,*)'stopping time             ',tstop
      write(*,*)'rewind file KE            ',tkeres

      write(*,*)'rewind cell vectors'
      write(*,'(1p,3e16.8)')celres(1),celres(2),celres(3)
      write(*,'(1p,3e16.8)')celres(4),celres(5),celres(6)
      write(*,'(1p,3e16.8)')celres(7),celres(8),celres(9)

      write(*,*)'rewind stress tensor'
      write(*,'(1p,3e16.8)')strres(1),strres(2),strres(3)
      write(*,'(1p,3e16.8)')strres(4),strres(5),strres(6)
      write(*,'(1p,3e16.8)')strres(7),strres(8),strres(9)

c     write basin difference data

      if(numbsn.gt.1)then

        write(*,*)'structural difference data'

        last=5*netdif
        allocate(buffer(last))
        read(7)buffer

        do i=1,last,5

          write(*,'(2i10,1p,3e16.8)')nint(buffer(i)),nint(buffer(i+1)),
     x      buffer(i+2),buffer(i+3),buffer(i+4)

        enddo

        deallocate(buffer)

      endif

c     write rewind configuration data

      write(*,*)'sample coordinates of rewind structure'
      allocate(buffer(3*natms))
      read(7)buffer

      j=1
      do i=1,3*natms,3

        write(*,'(i6,1p,3e16.8)')j,buffer(i),buffer(i+1),buffer(i+2)
        j=j+1

      enddo

      close(7)
      end
