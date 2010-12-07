      program potent
c     
c**********************************************************************
c     
c     dl_poly utility program
c
c     construct potential arrays for silver iodide
c     based on model by Ray, Rahman and Vashishta:
c     Superionics and Solid Electrolytes, Academic Press 1989
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith july 1994
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c
c**********************************************************************
c     
      implicit real*8(a-h,o-z)
      parameter (mxgrid=2000,ntable=10)
      dimension uuu(mxgrid,3),ggg(mxgrid,3),elrc(3),vlrc(3)
      character*8 silver,iodine
      write(*,'(a)')'dl_poly utility to create potential tables for'
      write(*,'(a)')'Ray-Rahman-Vashista silver iodide potential'
      write(*,*)
      write(*,'(a)')'Please specify potential cutoff range (rcut)'
      read(*,*) rcut
      write(*,'(a)')'Please specify number of grid points (mxgrid)'
      read(*,*) maxgrd
      if(maxgrd.gt.mxgrid)then
         write(*,'(a)')'Error - too many grid points requested'
         stop
      endif
      rpd=rcut/dble(maxgrd-4)
      silver='Ag+     '
      iodine='I-      '

c     
c     silver-iodine potential
      j=1
      do i=1,maxgrd
        rrr=dble(i)*rpd
        rsq=rrr**2
        uuu(i,j)=14.39d0*(114.48d0/rrr**5-1.1736d0)/rsq**2
        ggg(i,j)=14.39d0*(9.d0*114.48d0/rrr**5-4.d0*1.1736d0)/rsq**2
      enddo
      elrc(1)=14.39d0*(114.48d0/(6.d0*rcut**6)-1.1736d0/rcut)
      vlrc(1)=-14.39d0*(9.d0*114.48d0/(6.d0*rcut**6)-4.d0*1.1736d0/rcut)
c     
c     iodine-iodine potential
      j=2
      do i=1,maxgrd
        rrr=dble(i)*rpd
        rsq=rrr**2
        uuu(i,j)=14.39d0*(446.04d0/rrr**3-6.9331d0/rsq-2.3472d0)/rsq**2
        ggg(i,j)=14.39d0*(7.d0*446.04d0/rrr**3-6.d0*6.9331d0/rsq-
     x    4.d0*2.3472d0)/rsq**2
      enddo
      elrc(2)=14.39d0*(446.04d0/(4.d0*rcut**4)-6.9331d0/(3.d0*rcut**3)-
     x  2.3472d0/rcut)
      vlrc(2)=-14.39d0*(7.d0*446.04d0/(4.d0*rcut**4)-6.d0*6.9331d0/
     x  (3.d0*rcut**3)-4.d0*2.3472d0/rcut)
c     
c     silver-silver potential
      j=3
      do i=1,maxgrd
        rrr=dble(i)*rpd
        rsq=rrr**2
        uuu(i,j)=14.39d0*(0.014804d0/rrr**11)
        ggg(i,j)=14.39d0*(11.d0*0.014804d0/rrr**11)
      enddo
      elrc(3)=14.39d0*(0.014804d0/(8.d0*rcut**8))
      vlrc(3)=-11.d0*14.39d0*(0.014804d0/(8.d0*rcut**8))
c     
c     conversions to required units
      do  j=1,3
        elrc(j)=elrc(j)*1.60217733e-19/1.6605402e-23
        vlrc(j)=vlrc(j)*1.60217733e-19/1.6605402e-23
        do  i=1,maxgrd
          uuu(i,j)=uuu(i,j)*1.60217733e-19/1.6605402e-23
          ggg(i,j)=ggg(i,j)*1.60217733e-19/1.6605402e-23
        enddo
      enddo
c     
c     write out potentials
      open(ntable,file='TABLE')
      write(ntable,'(a)')'Silver Iodide Potential - DL_POLY units'
      write(ntable,'(2e15.8,i10)')rpd,rcut,maxgrd
      do k=1,3
        if(k.eq.1)write(ntable,'(2a8,1p,2e15.8)')
     x    silver,iodine,elrc(1),vlrc(1)
        if(k.eq.2)write(ntable,'(2a8,1p,2e15.8)')
     x    iodine,iodine,elrc(2),vlrc(2)
        if(k.eq.3)write(ntable,'(2a8,1p,2e15.8)')
     x    silver,silver,elrc(3),vlrc(3)
        write(ntable,'(4e15.8)')(uuu(i,k),i=1,maxgrd)
        write(ntable,'(4e15.8)')(ggg(i,k),i=1,maxgrd)
      enddo
      close (ntable)
      stop
      end
