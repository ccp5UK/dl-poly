      program potent
c     
c**********************************************************************
c     
c     dl_poly utility program
c
c     construct potential arrays for silica and related systems
c     based on Vessal's model for constant volume systems: 
c     Phil. Mag. B 60 (1989) 753-775
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
      parameter (mxgrid=5000,ntable=10)
      dimension uuu(mxgrid,2),ggg(mxgrid,2),elrc(2),vlrc(2)
      character*8 silicon,oxygen
      write(*,'(a)')'dl_poly utility to create potential tables'
      write(*,'(a)')'for Vessal constant volume silica potential'
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
      silicon='Si4+    '
      oxygen='O2-     '

c     
c     silicon-oxygen potential
      j=1
      do i=1,maxgrd
        rrr=dble(i)*rpd
        rsq=rrr**2
        if(rrr.lt.1.5d0)then
          uuu(i,j)=990.617d0*exp(-rrr/0.3297d0)
          ggg(i,j)=rsq*uuu(i,j)/(rrr*0.3297d0)
        else if(rrr.lt.2.5d0)then
          uuu(i,j)=728.7931332d0+rrr*(
     x         -1623.5166500d0+rrr*(1496.7899792d0+
     x         rrr*(-704.0045699d0+rrr*(167.0744823d0
     x         -rrr*15.8848138d0))))
          ggg(i,j)=-rsq*(-1623.5166500d0+rrr*(1496.7899792d0*2.d0+
     x         rrr*(-704.0045699d0*3.d0+rrr*(167.0744823d0*4.d0-
     x         rrr*15.8848138d0*5.d0))))/rrr
        else if(rrr.lt.3.5d0)then
          uuu(i,j)=0.6187898d0+rrr*(-0.6702747d0+rrr*(
     x         0.2214821d0-rrr*0.0233139d0))
          ggg(i,j)=-rsq*(-0.67027470+rrr*(0.2214821d0*2.d0-
     x         rrr*0.0233139d0*3.d0))/rrr
        else if(i.le.maxgrd)then
          uuu(i,j)=-25.d0/rrr**6
          ggg(i,j)=-6.d0*25.d0*rsq/rrr**8
        endif
      enddo
      elrc(1)=uuu(maxgrd-4,1)*rcut**3/3.d0
      vlrc(1)=-2.d0*uuu(maxgrd-4,1)*rcut**3
c     
c     oxygen-oxygen potential
      j=2
      do i=1,maxgrd
        rrr=dble(i)*rpd
        rsq=rrr**2
        if(rrr.lt.2.9d0)then
          uuu(i,j)=4511887.2d0*exp(-rrr/0.149d0)
          ggg(i,j)=rsq*uuu(i,j)/(rrr*0.149d0)
        else if(rrr.lt.3.6d0)then
          uuu(i,j)=298.0818367d0+rrr*(-453.1251863d0+
     x         rrr*(275.1200462d0+
     x         rrr*(-83.3698778d0+rrr*(12.6055881d0-
     x         rrr*0.7606781d0))))
          ggg(i,j)=-rsq*(-453.1251863d0+rrr*(275.1200462d0*2.d0+
     x         rrr*(-83.3698778d0*3.d0+rrr*(12.6055881d0*4.d0-
     x         rrr*0.7606781d0*5.d0))))/rrr
        else if(rrr.lt.4.2d0)then
          uuu(i,j)=1.5952103d0+rrr*(-1.2208242d0+rrr*
     x         (0.3052061d0-
     x         rrr*0.0251198d0))
          ggg(i,j)=-rsq*(-1.2208242d0+rrr*(0.3052061d0*2.d0-
     x         rrr*0.0251198d0*3.d0))/rrr
        else if(i.le.maxgrd)then
          uuu(i,j)=-52.12d0/rrr**6
          ggg(i,j)=-6.d0*52.12d0*rsq/rrr**8
        endif
      enddo
      elrc(2)=uuu(maxgrd-4,2)*rcut**3/3.d0
      vlrc(2)=-2.d0*uuu(maxgrd-4,2)*rcut**3
c     
c     conversions to required units
      do  j=1,2
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
      write(ntable,'(a)')'Silica Potential - DL_POLY units'
      write(ntable,'(1p,2e15.8,i10)')rpd,rcut,maxgrd
      do k=1,2
        if(k.eq.1)write(ntable,'(2a8,1p,2e15.8)')
     x    silicon,oxygen,elrc(1),vlrc(1)
        if(k.eq.2)write(ntable,'(2a8,1p,2e15.8)')
     x    oxygen,oxygen,elrc(2),vlrc(2)
        write(ntable,'(1p,4e15.8)')(uuu(i,k),i=1,maxgrd)
        write(ntable,'(1p,4e15.8)')(ggg(i,k),i=1,maxgrd)
      enddo
      close (ntable)
      stop
      end

