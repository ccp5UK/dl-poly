c      program xdenc

c*********************************************************************
c     
c     dl_poly program to calculate planar densities of all atom types
c     in a given direction from a dl_poly CONFIG file
c     
c     copyright daresbury laboratory 1996
c     author  w.smith oct 1996
c     
c     wl
c     1996/02/15 14:33:26
c     1.1.1.1
c     Exp
c
c*********************************************************************
      
      implicit real*8(a-h,o-z)
      
      parameter (ndiv=500,ntyp=4)
      character*40 fname
      character*80 title
      character*8 atname,atmnam(ntyp)
      dimension cell(9),xden(ndiv,ntyp)
        
c     open the I/O files
      
      open(5,file='xden_input')
      open(6,file='xden_output')
      
      write(6,'(a)')'# X-density Program'
      
c     name of selected CONFIG file
      read(5,*)fname

c     total number of atoms in a CONFIG file configuration
      read(5,*)natms
      
c     bin width for X density histogram
      read(5,*)delx

c     projection vector
      read(5,*)xx0,yy0,zz0
      rrr=sqrt(xx0**2+yy0**2+zz0**2)
      xx0=xx0/rrr
      yy0=yy0/rrr
      zz0=zz0/rrr

c     check on specified control variables
      
      write(6,'(a,a40)')'# name of target CONFIG file   : ',fname
      write(6,'(a,i8)')'# total no. of atoms in config  : ',natms
      write(6,'(a,1p,e12.4)')'# x-density histogram bin width : ',delx
      write(6,'(a,1p,3e12.4)')'# normalised projection vector  : ',
     x    xx0,yy0,zz0
        
c     open CONFIG file

      open(7,file=fname)

      read(7,'(a80)')title
      write(6,'(a10,a80)')'# Header: ',title
      read(7,*)levcfg,imcon

      if(imcon.gt.0)then

        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)

        write(6,'(a)')     '# cell vectors:'
        write(6,'(a,3f12.6)')'# vector A :',cell(1),cell(2),cell(3)
        write(6,'(a,3f12.6)')'# vector B :',cell(4),cell(5),cell(6)
        write(6,'(a,3f12.6)')'# vector C :',cell(7),cell(8),cell(9)

      endif

c     initialise x-density array

      
      do j=1,ntyp
        do i=1,ndiv
          
          xden(i,j)=0.d0
          
        enddo
      enddo
      
c     read configuration

      ityp=0
      nmin=ndiv/2
      nmax=ndiv/2

      do i=1,natms

        read(7,'(a8)')atname

        read(7,*)xxx,yyy,zzz
        if(levcfg.gt.0)read(7,*)uuu,vvv,www
        if(levcfg.gt.1)read(7,*)uuu,vvv,www

        ktyp=0
        do jtyp=1,ityp

          if(atname.eq.atmnam(jtyp))ktyp=jtyp

        enddo
        if(ktyp.eq.0.and.ityp.lt.ntyp)then

          ityp=ityp+1
          atmnam(ityp)=atname
          ktyp=ityp

        endif

        if(ktyp.le.ntyp)then

          rrr=xxx*xx0+yyy*yy0+zzz*zz0
          ix=nint(rrr/delx)+ndiv/2
          if(ix.gt.0.and.ix.le.ndiv)then

            xden(ix,ktyp)=xden(ix,ktyp)+1.d0
            nmin=min(nmin,ix)
            nmax=max(nmax,ix)

          endif

        endif

      enddo

c     print out final histogram

      do i=nmin,nmax

        xxx=delx*dble(i-ndiv/2)
        write(6,'(1p,5e12.4)')xxx,(xden(i,j),j=1,ityp)

      enddo

      stop
      
      end
