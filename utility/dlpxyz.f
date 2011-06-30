c      program dlpxyz
c*********************************************************************
c
c     dl_poly program to convert CONFIG and REVCON files to
c     a standard XYZ file for the XMOL package
c
c     copyright daresbury laboratory 1997
c     author  w.smith april 1997
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      parameter (mxatms=2490)
      character*40 fname
      character*80 title
      character*8 atmnam(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)

c     open I/O files

c      open(5,file='dlpxyz_input')
c      open)6,file='dlpxyz_output')

      write(6,'(a)')
     x  'CONFIG/REVCON to XYZ File Converter'

      write(6,'(a)')'Enter name of file to convert'
      read(5,*)fname
      write(6,'(a)')'Output file will be called XYZ.xyz'
      write(6,'(a)')'Enter number of atoms in MD cell'
      read(5,*)natms
      if(natms.gt.mxatms)then
        write(6,*)'Error - too many atoms in cell. Max=',mxatms
        stop
      endif

      open(7,file=fname)

      read(7,'(a80)')title
      read(7,'(2i10)')levcfg,imcon
      if(imcon.gt.0)then
        read(7,'(3f20.0)')a,b,c
        read(7,'(3f20.0)')a,b,c
        read(7,'(3f20.0)')a,b,c
      endif
      do i=1,natms

        read(7,'(a8)',end=100)atmnam(i)
        read(7,'(3f20.0)')xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)read(7,'(3f20.0)')a,b,c
        if(levcfg.gt.1)read(7,'(3f20.0)')a,b,c

      enddo

  100 close (7)
      natms=i-1

      do i=1,natms

        call strip(atmnam(i),8)
        if(atmnam(i)(1:1).eq."K")then
          atmnam(i)="K"
        else if(atmnam(i)(1:2).eq."Ti")then
          atmnam(i)="Ti"
        else if(atmnam(i)(1:2).eq."Cl")then
          atmnam(i)="Cl"
        else if(atmnam(i)(1:2).eq."Ag")then
          atmnam(i)="Ag"
        else if(atmnam(i)(1:1).eq."I")then
          atmnam(i)="I"
        else if(atmnam(i)(1:1).eq."C")then
          atmnam(i)="C"
        else if(atmnam(i)(1:1).eq."N")then
          atmnam(i)="N"
        else if(atmnam(i)(1:1).eq."O")then
          atmnam(i)="O"
        else if(atmnam(i)(1:1).eq."H")then
          atmnam(i)="H"
        endif

      enddo
      open(8,file='XYZ.xyz')
      write(8,'(i5)')natms
      write(8,'(a80)')title
      do i=1,natms

        write(8,'(a2,3f15.8)')atmnam(i)(1:2),xxx(i),yyy(i),zzz(i)

      enddo

      close (8)
      end
      subroutine strip(string,length)

c***********************************************************************
c
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c
c     copyright daresbury laboratory 1993
c
c***********************************************************************

      character*(*) string

      imax = min(length,255)
      do i = 1,imax

        if(string(1:1).eq.' ') then

          do j = 1,imax-1

            string(j:j) = string(j+1:j+1)

          enddo

          string(imax:imax) = ' '

        endif

      enddo

      return
      end
