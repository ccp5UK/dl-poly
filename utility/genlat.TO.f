      program genlat2
c
c***********************************************************************
c
c     program to generate a truncated octahedral cell from a perfect
c     lattice
c
c     copyright daresbury laboratory
c     author - w.smith oct. 1992
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter (mxatom=10)
      dimension xb(mxatom),yb(mxatom),zb(mxatom)
      character*8 name(mxatom)
      character*80 title
      write(*,*)'enter suitable one-line title for lattice file'
      read(*,'(a80)')title
      write(*,*)'enter number of basis atoms (max 10)'
      read(*,*)natms
      write(*,*)'define cubic unit cell width'
      read(*,*)cell_param
      write(*,*)'define required spherical cutoff'
      read(*,*)radius
      write(*,*)'enter ion names (8 chars) and fractional coordinates'
      wdth=radius*4.d0/sqrt(3.d0)
      num=int(wdth/cell_param)
      if(dble(num)*cell_param.lt.wdth)num=num+1
      wdth=dble(num)*cell_param
      test=0.75d0*wdth
      density=dble(natms)/cell_param**3
      do i=1,natms
         read(*,*)name(i),xb(i),yb(i),zb(i)
      enddo
      open(10,file='LATTICE')
      write(10,'(a80)')title
      write(10,'(2i10)')0,4
      write(10,'(3f20.8)')wdth,0.d0,0.d0,0.d0,wdth,0.d0,0.d0,0.d0,wdth
c
c     set up lattice
      matms=0
      do  n=1,natms
         do  k=1,num
            do  j=1,num
               do  i=1,num
                  xs=dble(i-1)+xb(n)-0.5d0*dble(num)
                  xx=cell_param*xs
                  if(abs(xx).lt.0.5d0*wdth)then
                     yy=cell_param*ys
                     ys=dble(j-1)+yb(n)-0.5d0*dble(num)
                     if(abs(yy).lt.0.5d0*wdth)then
                        zs=dble(k-1)+zb(n)-0.5d0*dble(num)
                        zz=cell_param*zs
                        if(abs(zz).lt.0.5d0*wdth)then
                           if(abs(xx)+abs(yy)+abs(zz).lt.test)then
                              matms=matms+1
                              write(10,'(a8,/,3f20.6)')name(n),xx,yy,zz
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
      tdensty=dble(matms)/(0.5d0*wdth**3)
      actcut=wdth*sqrt(3.d0)/4.d0
      write(*,*)'width parameter of TO cell =',wdth
      write(*,*)'number of ions in system = ',matms
      write(*,*)'target density =',density
      write(*,*)'actual density =',tdensty
      write(*,*)'target cutoff  =',radius
      write(*,*)'maximum cutoff =',actcut
      stop
      end
