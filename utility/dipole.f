      program dipole

c*********************************************************************
c     
c     dl_poly program to calculate cluster dipoles
c     for selected atom clusters from dl_poly HISTORY file
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
      
      parameter (mxatms=1080,mxcls=2)
      character*1 fmt
      character*40 fname
      character*80 title
      character*8 atmnam(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension weight(mxatms),chge(mxatms),cell(9)
      dimension cmass(3,mxcls),dipole(3,mxcls)
      dimension ibc(mxcls),iec(mxcls)
        
c     open the I/O files
      
      open(5,file='dip_input')
      open(6,file='dip_output')
      
      write(6,'(a)')
     x  '# Cluster Dipole Program'
      
c     read the control variables
      
c     number of clusters
      read(5,*)ncls

c     specify atom clusters

      do i=1,ncls
        read(5,*)ibc(i),iec(i)
      enddo
      
c     name of selected HISTORY file
      read(5,*)fname

c     HISTORY file formatted or unformatted
      read(5,*)fmt
      
c     total number of atoms in a HISTORY file configuration
      read(5,*)natms
      
c     max number of configurations to sample
      read(5,*)nconf
      
c     read time interval between configurations
      read(5,*)tstep

c     check on specified control variables
      
      write(6,'(a,a40)')'# name of target HISTORY file   : ',fname
      write(6,'(a,i8)')'# number of atomic clusters     : ',ncls
      do i=1,ncls
        write(6,'(a,2i8)')'# start and end atoms   : ',ibc(i),iec(i)
      enddo
      write(6,'(a,i8)')'# total no. of atoms in config  : ',natms
      write(6,'(a,i8)')'# number of configurations      : ',nconf
      write(6,'(a,f12.6)')'# time interval between configs : ',tstep
        
      if(natms.gt.mxatms)then
        
        write(6,'(a,2i5)')
     x    '# error - too many atoms in system'
        stop
        
      endif
      
c     set default cell properties
      
      do i=1,9
        cell(i)=0.d0
      enddo
      cell(1)=1.d0
      cell(5)=1.d0
      cell(9)=1.d0
      
      iflg=0
      keytrj=0

      do iconf=0,nconf-1
        
        if(fmt.eq."U".or.fmt.eq."u")then

          call uread
     x      (fname,title,atmnam,iflg,imcon,keytrj,natms,nstep,timstp,
     x       cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

        else

          call hread
     x      (fname,title,atmnam,iflg,imcon,keytrj,natms,nstep,timstp,
     x       cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

        endif

        if(imcon.ne.0)then

          write(6,*)'error - periodic boundary condition in operation'
          stop

        endif

        if(iflg.lt.0)go to 100

c     initialise dipole and centre of mass arrays
      
        do i=1,mxcls
          
          cmass(1,i)=0.d0
          cmass(2,i)=0.d0
          cmass(3,i)=0.d0
          dipole(1,i)=0.d0
          dipole(2,i)=0.d0
          dipole(3,i)=0.d0
          
        enddo
      

        do j=1,ncls

          huge=0.d0
          
          do i=ibc(j),iec(j)
            
            huge=huge+weight(i)
            
            cmass(1,j)=weight(i)*xxx(i)+cmass(1,j)
            cmass(2,j)=weight(i)*yyy(i)+cmass(2,j)
            cmass(3,j)=weight(i)*zzz(i)+cmass(3,j)
            
          enddo
          
          cmass(1,j)=cmass(1,j)/huge
          cmass(2,j)=cmass(2,j)/huge
          cmass(3,j)=cmass(3,j)/huge
          
          do i=ibc(j),iec(j)
          
            dipole(1,j)=chge(i)*(xxx(i)-cmass(1,j))+dipole(1,j)
            dipole(2,j)=chge(i)*(yyy(i)-cmass(2,j))+dipole(2,j)
            dipole(3,j)=chge(i)*(zzz(i)-cmass(3,j))+dipole(3,j)
            
          enddo

        enddo

        time=tstep*dble(iconf)
        write(6,'(1p,10e12.4)')time,dipole
        
      enddo

  100 continue

      stop
      
      end
      subroutine uread
     x  (history,cfgname,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading unformatted history files
c     
c     double precision, single processor version
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c     
c     wl
c     1996/02/15 14:33:26
c     1.1.1.1
c     Exp
c     
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)
      
      logical new
      
      character*80 cfgname
      character*40 history
      character*8 atmnam(*)
      
      dimension cell(9)
      dimension chge(*),weight(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension fxx(*),fyy(*),fzz(*)
      
      save new
      
      data new/.true./,nhist/77/
      
      iflg=0

c     open the history file if new job
      
      if(new)then
        
        open(nhist,file=history,form='unformatted',
     x       status='old',err=100)
        
        read(nhist,err=200) cfgname
        write(*,'(a,a)')'# History file header: ',cfgname
        read(nhist,end=200) datms
        if(natms.ne.nint(datms))then
          
          matms=nint(datms)
          write(*,'(a)')'# error - incorrect number of atoms in file'
          write(*,'(a,i6,a)')'# file contains',matms,' atoms'
          stop
          
        endif
        read(nhist,end=200) (atmnam(i),i=1,natms)
        read(nhist,end=200) (weight(i),i=1,natms)
        read(nhist,end=200) (chge(i),i=1,natms)
        
      endif
      
      read(nhist,end=200)dstep,datms,trjkey,dimcon
      nstep=nint(dstep)
      ktrj=nint(trjkey)
      imcon=nint(dimcon)
      if(keytrj.gt.ktrj)then
        
        if(ktrj.eq.0)write(*,'(a)')'# error - no velocities in file'
        if(keytrj.gt.1)write(*,'(a)')'# error - no forces in file'
        stop

      endif
      
      if(imcon.gt.0) read(nhist,end=200) cell
      
      read(nhist,end=200) (xxx(i),i = 1,natms)
      read(nhist,end=200) (yyy(i),i = 1,natms)
      read(nhist,end=200) (zzz(i),i = 1,natms)
      
      if(keytrj.ge.1)then
        read(nhist,end=200) (vxx(i),i = 1,natms)
        read(nhist,end=200) (vyy(i),i = 1,natms)
        read(nhist,end=200) (vzz(i),i = 1,natms)
      else if(ktrj.ge.1)then
        read(nhist,end=200)
        read(nhist,end=200)
        read(nhist,end=200)
      endif
      if(keytrj.ge.2)then
        read(nhist,end=200) (fxx(i),i = 1,natms)
        read(nhist,end=200) (fyy(i),i = 1,natms)
        read(nhist,end=200) (fzz(i),i = 1,natms)
      else if(ktrj.ge.2)then
        read(nhist,end=200)
        read(nhist,end=200)
        read(nhist,end=200)
      endif

      iflg=1
      
      return

  100 continue

      write(*,'(a)')'# error - History file not found'
      stop

  200 continue
      write(*,'(a)')'# warning - end of History file encountered'
      close (nhist)
      iflg=-1

      return
      end
      subroutine hread
     x  (history,cfgname,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading the formatted history file 
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c
c     single processor version
c
c     wl
c     1996/02/15 14:33:26
c     1.1.1.1
c     Exp
c     
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)

      logical new

      character*80 cfgname
      character*40 history
      character*8 atmnam(*),step
      
      dimension cell(9)
      dimension chge(*),weight(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension fxx(*),fyy(*),fzz(*)

      save new

      data new/.true./,nhist/77/

      iflg=0

c     open history file if new job

      if(new)then
      
        open(nhist,file=history,status='old',err=100)
        
        read(nhist,'(a80)',err=200) cfgname
        write(*,'(a,a)')'# History file header: ',cfgname
        read(nhist,'(2i10)',end=200) ktrj,imcon
        if(keytrj.gt.ktrj)then

          if(ktrj.eq.0)write(*,'(a)')'# error - no velocities in file'
          if(keytrj.gt.1)write(*,'(a)')'# error - no forces in file'
          stop

        endif

        new=.false.

      endif
        
      read(nhist,*,end=200)
     x     step,nstep,matms,ktrj,imcon,tstep
      
      if(natms.ne.matms)then

        write(*,'(a)')'# error - incorrect number of atoms in file'
        write(*,'(a,i6,a)')'# file contains',matms,' atoms'
        stop

      endif
      
      if(imcon.gt.0) read(nhist,'(3g12.4)',end=200) cell
      
      do i = 1,natms
        read(nhist,'(a8,i10,2f12.6)',end=200)
     x    atmnam(i),j,weight(i),chge(i)
        read(nhist,'(1p,3e12.4)',end=200) xxx(i),yyy(i),zzz(i)
        if(keytrj.ge.1)then
          read(nhist,'(1p,3e12.4)',end=200) vxx(i),vyy(i),vzz(i)
        else if(ktrj.ge.1)then
          read(nhist,'(1p,3e12.4)',end=200) vx,vy,vz
        endif
        if(keytrj.ge.2)then
          read(nhist,'(1p,3e12.4)',end=200) fxx(i),fyy(i),fzz(i)
        else if(ktrj.ge.2)then
          read(nhist,'(1p,3e12.4)',end=200) fx,fy,fz
        endif
      enddo
      
      iflg=1

      return

  100 continue

      write(*,'(a)')'# error - History file not found'
      stop

  200 continue
      write(*,'(a)')'# warning - end of History file encountered'
      close (nhist)
      iflg=-1

      return
      end

