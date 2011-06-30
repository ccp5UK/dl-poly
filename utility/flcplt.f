c*********************************************************************
c
c     dl_poly utility program to extract energy fluctuations from
c     STATIS file data for plotting
c
c     copyright daresbury laboratory
c     author w.smith feb 2010
c
c*********************************************************************

      implicit none
      logical findstring
      character record(40)
      integer i,j,k,m,df,idum
      real(8) x,a,b,c,d,e,f,g,h,eu,bb,fac

c     ask number of degrees of freedom

      write(*,*)'DL_POLY energy fluctuation program'
      write(*,*)'Enter number of degrees of freedom in system:'
      read(*,*)df
      write(*,*)'Processing STATIS file...'

c     open STATIS file

      open (77,file='STATIS')

c     determine energy units of STATIS file

      eu=0.d0
      bb=0.831451115d0
      read(77,*)
      read(77,'(40a1)')record
      if(findstring('electron Volts',record,idum))then
        eu=9648.530821d0
      elseif(findstring('kilo electron Volts',record,idum))then
        eu=9648530.821d0
      elseif(findstring('kcal/mol',record,idum))then
        eu=418.4d0
      elseif(findstring('kjoule/mol',record,idum))then
        eu=1.d2
      elseif(findstring('kelvin',record,idum))then
        eu=bb
      elseif(findstring('Internal Units',record,idum))then
        eu=1.d0
      endif
      fac=bb*dble(df)/eu

c     calculate averages from STATIS file contents

      i=0
      d=0.d0
      e=0.d0
      f=0.d0
      do while(.true.)
        read(77,*,end=100)k,x,m
        read(77,*)a,b,c
        do j=1,(m+4)/5-1
          read(77,*)
        enddo
        if(k.gt.0)then
          i=i+1
          d=d+a
          e=e+b
          f=f+c
        endif
      enddo
  100 continue
      close (77)

c     calculate averages

      d=d/dble(i)
      e=e/dble(i)
      f=f/dble(i)

c     calculate fluctuations from STATIS file contents

      open (77,file='STATIS')
      open (88,file='FLCPLT')

      read(77,*)
      read(77,*)

      i=0
      do while(.true.)
        read(77,*,end=200)k,x,m
        read(77,*)a,b,c
        do j=1,(m+4)/5-1
          read(77,*)
        enddo
        if(k.gt.0)then
          i=i+1
          a=a-d
          b=(b-e)*fac
          c=c-f
          g=b+c
          h=a-g
          write(88,'(1p,6e14.6)')x,a,b,c,g,h
        endif
      enddo
  200 continue

      close (77)
      close (88)

      write(*,*)'FLCPLT file written'
      write(*,*)'End of Job'

      end

      logical function findstring(seek,string,here)

c***********************************************************************
c
c     DL_POLY routine to find an explicit string in an input record
c     note: variable `seek' is a character string while variable
c    `string' is a character*1 array i.e. code is application specific
c
c     copyright daresbury laboratory
c     author    w.smith   jan   2004
c
c***********************************************************************

      implicit none

      integer i,n,m,here
      character*(*) seek
      character*1 string(40)

      m=40
      n=len(seek)
      findstring=.false.

      here=0
      do while(here.le.m-n.and.(.not.findstring))

        findstring=.true.

        do i=1,n
          if(seek(i:i).ne.string(here+i))findstring=.false.
        enddo

        here=here+1

      enddo

      return
      end function findstring
