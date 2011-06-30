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

        new=.false.

      endif

      read(nhist,end=200)dstep,datms,trjkey,dimcon,tstep
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


