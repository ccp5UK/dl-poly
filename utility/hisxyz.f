c     program hisxyz
c*********************************************************************
c     
c     dl_poly program to convert HISTORY files to
c     a standard XYZ file for the XMOL package
c     
c     copyright daresbury laboratory 1997
c     author  w.smith june 1997
c     
c     wl
c     1996/02/15 14:33:26
c     1.1.1.1
c     Exp
c     
c*********************************************************************
      
      implicit real*8(a-h,o-z)
      
      parameter (mxatms=2490)
      character*40 fname
      character*80 title
      character*8 atmnam(mxatms),newnam(mxatms),dummy
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      
      write(*,'(a)')
     x  'HISTORY to XYZ File Converter'
      
      write(*,'(a)')'Enter name of file to convert'
      read(5,*)fname
      write(*,'(a)')'Output file will be called HIS.xyz'
      write(*,'(a)')'Enter number of configurations required'
      read(5,*)nconf
      write(*,'(a)')'Enter number of atoms in file'
      read(5,*)natms

      open(7,file=fname)
      open(8,file='HIS.xyz')

      read(7,'(a80)')title
      read(7,'(2i10)')levcfg,imcon

      do iconf=1,nconf

        read(7,'(a)')dummy

        if(imcon.gt.0)then
          read(7,*)a,b,c
          read(7,*)a,b,c
          read(7,*)a,b,c
        endif
        do i=1,natms

          read(7,'(a8)',end=100)atmnam(i)
          read(7,*)xxx(i),yyy(i),zzz(i)
          if(levcfg.gt.0)read(7,*)a,b,c
          if(levcfg.gt.1)read(7,*)a,b,c

        enddo

        if(iconf.eq.1)call whatname(natms,newnam,atmnam)

        write(8,'(i5)')natms
        write(8,'(a80)')title
        do i=1,natms
          
          write(8,'(a2,3f12.4)')newnam(i)(1:2),xxx(i),yyy(i),zzz(i)

        enddo

      enddo

  100 close (7)
      close (8)
      end
      subroutine strip(string,length)

c***********************************************************************
c     
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
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
      subroutine whatname(natms,newnam,atmnam)
c*********************************************************************
c     
c     dl_poly utility program
c     reassign atom names in a known list of atoms
c     copyright daresbury laboratory jun 1997
c     author w.smith jun 1997
c     
c*********************************************************************

      character*8 atmnam(*),newnam(*)

      do i=1,natms

        newnam(i)="X"
        call strip(atmnam(i),8)
        if(atmnam(i)(1:1).eq."K")then
          newnam(i)="K"
        else if(atmnam(i)(1:2).eq."Ti")then
          newnam(i)="Ti"
        else if(atmnam(i)(1:2).eq."Cl")then
          newnam(i)="Cl"
        else if(atmnam(i)(1:2).eq."Ag")then
          newnam(i)="Ag"
        else if(atmnam(i)(1:1).eq."I")then
          newnam(i)="I"
        else if(atmnam(i)(1:1).eq."C")then
          newnam(i)="C"
        else if(atmnam(i)(1:1).eq."N")then
          newnam(i)="N"
        else if(atmnam(i)(1:1).eq."O")then
          newnam(i)="O"
        else if(atmnam(i)(1:1).eq."H")then
          newnam(i)="H"
        endif

      enddo

      return
      end
