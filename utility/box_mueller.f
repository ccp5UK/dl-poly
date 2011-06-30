      subroutine box_mueller(gran1,gran2)

c*********************************************************************
c
c     box mueller method for generating gaussian random numbers
c      with a zero mean and standard deviation of 1
c
c     duni() is a double precision random number on interval [0,1)
c     gran1 and gran2 are both gaussian random numbers
c
c     copyright daresbury laboratory
c     author w.smith may 2008
c
c*********************************************************************

      implicit none

      logical newjob
      real(8) ran0,ran1,ran2,gran1,gran2
      save newjob
      data newjob/.true./

c     make sure duni is initialised

      if(newjob)then

        newjob=.false.
        ran0=duni()

      endif

      ran0=1.d0

c     generate uniform random numbers on [-1, 1)

      do while(ran0.ge.1.d0)

        ran1=2.d0*duni()-1.d0
        ran2=2.d0*duni()-1.d0
        ran0=ran1**2+ran2**2

      enddo

c     calculate gaussian random numbers

      ran0=sqrt((-2.d0*log(ran0))/ran0)
      gran1=ran0*ran1
      gran2=ran0*ran2

      return
      end

