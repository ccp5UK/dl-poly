      program gdiff

c*********************************************************************
c
c     dl_poly program to calculate cross term of Van Hove density
c     correlation function for selected atoms
c
c     copyright daresbury laboratory 1995
c     author  w.smith march 1995
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      parameter (mxatms=1080,nstr=240,mxcrs=64,mxtim=64)
      parameter (pi=3.1415926536d0)

      character*1 fmt
      character*8 atnam1
      character*8 atnam2
      character*40 fname
      character*80 title
      character*8 atmnam(mxatms),name(nstr)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension weight(mxatms),chge(mxatms)
      dimension xx0(nstr),yy0(nstr),zz0(nstr)
      dimension xx1(nstr),yy1(nstr),zz1(nstr)
      dimension xx2(nstr),yy2(nstr),zz2(nstr)
      dimension gcross(mxcrs,mxtim,4),gzero(nstr,mxtim,3)
      dimension avcell(9),cell(9),rcell(9),time(7)
      dimension imd(mxtim),msm(mxtim)

c     open the I/O files
      open(5,file='gdiff_input')
      open(6,file='gdiff_output')

      write(6,'(a)')
     x  '# van Hove Cross Correlation Function Program'

c     read the control variables

c     name of species for autocorrelation
      read(5,*)atnam1
      read(5,*)atnam2

c     name of selected HISTORY file
      read(5,*)fname

c     HISTORY file formatted or unformatted
      read(5,*)fmt

c     total number of atoms in a HISTORY file configuration
      read(5,*)natms

c     time interval between successive configurations
      read(5,*)tstep

c     max number of configurations
      read(5,*)nconf

c     sampling interval
      read(5,*)isampl

c     interval between origins
      read(5,*)iogcrs

c     maximum radius of correlation
      read(5,*)rcut

c     average cell vectors (taken from dl_poly OUTPUT file)
      read(5,*)avcell

c     working length of gcross array
      ngcrs=mxtim

      delr=rcut/dble(mxcrs)
      rcut2=rcut**2

c     check on specified control variables

      if(mod(ngcrs,iogcrs).ne.0)then

        ngcrs=iogcrs*(ngcrs/iogcrs)
        write(6,'(a,i5)')
     x    '# warning - gcross array reset to ',ngcrs

      endif

      nogcrs=ngcrs/iogcrs

      write(6,'(a,a40)')'# label  of target HISTORY file : ',fname
      write(6,'(a,a8)')'# label  of first atom          : ',atnam1
      write(6,'(a,a8)')'# label  of second atom         : ',atnam2
      write(6,'(a,i6)')'# length of correlation arrays  : ',ngcrs
      write(6,'(a,i6)')'# number of requested configs.  : ',nconf
      write(6,'(a,i6)')'# sampling interval             : ',isampl
      write(6,'(a,i6)')'# interval between origins      : ',iogcrs
      write(6,'(a,f6.3)')'# time interval between configs : ',tstep
      write(6,'(a,f6.3)')'# selected correlation radius   : ',rcut
      write(6,'(a,f6.3)')'# bin width in correlation      : ',delr
      write(6,'(a)')     '# average cell vectors:'
      write(6,'(a,3f12.6)')'# vector A :',avcell(1),avcell(2),avcell(3)
      write(6,'(a,3f12.6)')'# vector B :',avcell(4),avcell(5),avcell(6)
      write(6,'(a,3f12.6)')'# vector C :',avcell(7),avcell(8),avcell(9)

      if(natms.gt.mxatms)then

        write(6,'(a,2i5)')
     x    'error - too many atoms in system',natms,mxatms
        stop

      endif

c     initialise gcross arrays

      do j=1,mxtim
        msm(j)=0
        do i=1,mxcrs

          gcross(i,j,1)=0.d0
          gcross(i,j,2)=0.d0
          gcross(i,j,3)=0.d0
          gcross(i,j,4)=0.d0

        enddo
      enddo

c     multiple pass over HISTORY file

      do ipass=0,isampl-1
c      do ipass=0,0

        lsr=0
        msr=0
        nsgcrs=0

c     set default cell vectors

        do i=1,9
          cell(i)=0.d0
        enddo
        cell(1)=1.d5
        cell(5)=1.d5
        cell(9)=1.d5

        iflg=0
        keytrj=0

        do iconf=0,nconf-1

          if(fmt.eq."U".or.fmt.eq."u")then

            call uread
     x        (fname,title,atmnam,iflg,imcon,keytrj,natms,nstep,
     x        cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


          else

            call hread
     x        (fname,title,atmnam,iflg,imcon,keytrj,natms,nstep,
     x        cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


          endif

          if(iflg.lt.0)go to 100

          if(imcon.ge.1.and.imcon.le.3)then

            call invert(cell,rcell,det)

          else

            write(6,'(a)')
     x        '# error - incorrect boundary condition'
            stop

          endif

          j0=0
          do i=1,natms

            if(atmnam(i).eq.atnam1.or.atmnam(i).eq.atnam2)then

              j0=j0+1
              if(j0.gt.nstr)then

                write(6,'(a,2i5)')
     x            'error - too many atoms of selected type',j0,nstr
                stop

              endif

              name(j0)=atmnam(i)
              xx2(j0)=xxx(i)*rcell(1)+yyy(i)*rcell(4)+zzz(i)*rcell(7)
              yy2(j0)=xxx(i)*rcell(2)+yyy(i)*rcell(5)+zzz(i)*rcell(8)
              zz2(j0)=xxx(i)*rcell(3)+yyy(i)*rcell(6)+zzz(i)*rcell(9)

            endif

          enddo

          matms=j0

          if(atnam1.eq.atnam2)then

            ntrm=1
          else

            ntrm=4

          endif

          if(iconf.eq.ipass)then

c     set initial positions of atoms

            do i=1,matms

              xx0(i)=xx2(i)
              yy0(i)=yy2(i)
              zz0(i)=zz2(i)

            enddo

          else if(iconf.gt.ipass)then

c     accumulate incremental distances

            do i=1,matms

              uuu=xx2(i)-xx1(i)
              vvv=yy2(i)-yy1(i)
              www=zz2(i)-zz1(i)

              uuu=uuu-nint(uuu)
              vvv=vvv-nint(vvv)
              www=www-nint(www)

              xx0(i)=xx0(i)+uuu
              yy0(i)=yy0(i)+vvv
              zz0(i)=zz0(i)+www

            enddo

          endif

          do i=1,matms

            xx1(i)=xx2(i)
            yy1(i)=yy2(i)
            zz1(i)=zz2(i)

          enddo

c     calculate cross correlation function

          if(mod(iconf,isampl).eq.ipass)then

            if(mod(nsgcrs,iogcrs).eq.0)then

              lsr=min(lsr+1,nogcrs)
              msr=mod(msr,nogcrs)+1
              imd(msr)=1
              do i=1,matms

                gzero(i,msr,1)=xx0(i)
                gzero(i,msr,2)=yy0(i)
                gzero(i,msr,3)=zz0(i)

              enddo

            endif

            nsgcrs=nsgcrs+1

            do j=1,lsr

              m=imd(j)
              imd(j)=m+1
              msm(m)=msm(m)+1
              do ii=1,matms
                do jj=1,matms

                  if(ii.ne.jj)then

                    uuu=xx0(ii)-gzero(jj,j,1)
                    vvv=yy0(ii)-gzero(jj,j,2)
                    www=zz0(ii)-gzero(jj,j,3)

                    uuu=uuu-nint(uuu)
                    vvv=vvv-nint(vvv)
                    www=www-nint(www)

                    rmsx=uuu*avcell(1)+vvv*avcell(4)+www*avcell(7)
                    rmsy=uuu*avcell(2)+vvv*avcell(5)+www*avcell(8)
                    rmsz=uuu*avcell(3)+vvv*avcell(6)+www*avcell(9)

                    rsq=rmsx**2+rmsy**2+rmsz**2

                    if(rsq.lt.rcut2)then
                      k=int(sqrt(rsq)/delr)+1

                      if(name(ii).eq.atnam1)then

                        if(name(jj).eq.atnam1)then

                          gcross(k,m,1)=gcross(k,m,1)+1.d0

                        else

                          gcross(k,m,2)=gcross(k,m,2)+1.d0

                        endif

                      else

                        if(name(jj).eq.atnam1)then

                          gcross(k,m,3)=gcross(k,m,3)+1.d0

                        else

                          gcross(k,m,4)=gcross(k,m,4)+1.d0

                        endif

                      endif

                    endif

                  endif

                enddo
              enddo

            enddo

          endif

        enddo

  100   continue

      enddo

      if(iflg.lt.0)iconf=iconf-1
      last=min(nsgcrs,ngcrs)
      write(6,'(a,i6)')
     x  '# number of configurations sampled: ',iconf

c     normalise and print correlation function

      tstep=tstep*dble(isampl)

      volm=abs(avcell(1)*(avcell(5)*avcell(9)-avcell(6)*avcell(8))
     x  +avcell(4)*(avcell(3)*avcell(8)-avcell(2)*avcell(9))
     x  +avcell(7)*(avcell(2)*avcell(6)-avcell(3)*avcell(5)))

      nta=0
      ntb=0
      do i=1,matms

        if(name(i).eq.atnam1)nta=nta+1
        if(name(i).eq.atnam2)ntb=ntb+1

      enddo

      do k=1,ntrm

        if(k.eq.1) then

          rnorm=volm/(4.d0*pi*delr**3*dble(nta)*dble(nta))

        else if(k.eq.4)then

          rnorm=volm/(4.d0*pi*delr**3*dble(ntb)*dble(ntb))

        else

          rnorm=volm/(4.d0*pi*delr**3*dble(nta)*dble(ntb))

        endif

        do j=1,last

          do i=1,mxcrs

            gcross(i,j,k)=(rnorm*gcross(i,j,k)/dble(msm(j)))
     x        /((dble(i)-0.5d0)**2+1.d0/12.d0)

          enddo

        enddo

c     print out final correlation function


        do j=1,last,7

          do n=0,6

            time(n+1)=tstep*dble(n+j-1)

          enddo

          write(6,'(a10,7f10.5)')'# rad/time',(time(n),n=1,7)

          do i=1,mxcrs

            rrr=delr*(dble(i)-0.5d0)
            write(6,'(8f10.5)')rrr,(gcross(i,n,k),n=j,min(j+6,last))

          enddo

          write(6,'(a)')'&'

        enddo

      enddo

      stop

      end
      subroutine invert(a,b,d)
c
c***********************************************************************
c
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c
c***********************************************************************
c

      real*8 a,b,d,r

      dimension a(9),b(9)
c
c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)
c
c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d
c
c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)
      return
      end
      subroutine hread
     x  (history,cfgname,atmnam,iflg,imcon,keytrj,natms,nstep,
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
c***********************************************************************
c

      implicit real*8(a-h,o-z)

      character*80 cfgname
      character*40 history
      character*8 atmnam(*),step

      dimension cell(9)
      dimension chge(*),weight(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension fxx(*),fyy(*),fzz(*)

      data nhist/77/

c     open history file if new job

      if(iflg.eq.0)then

        close(nhist)
        open(nhist,file=history,status='old', err=100)

        read(nhist,'(a80)',err=200) cfgname
        write(*,'(a,a)')'# History file header: ',cfgname
        read(nhist,'(2i10)',end=200) ktrj,imcon
        if(keytrj.gt.ktrj)then

          if(ktrj.eq.0)write(*,'(a)')'# error - no velocities in file'
          if(keytrj.gt.1)write(*,'(a)')'# error - no forces in file'
          stop

        endif

      endif

      read(nhist,'(a8,4i10)',end=200)step,nstep,matms,ktrj,imcon

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
      subroutine uread
     x  (history,cfgname,atmnam,iflg,imcon,keytrj,natms,nstep,
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

      read(nhist,end=200)dstep,datms,trjkey,dimcon
      if(natms.ne.nint(datms))then

        matms=nint(datms)
        write(*,'(a)')'# error - incorrect number of atoms in file'
        write(*,'(a,i6,a)')'# file contains',matms,' atoms'
        stop

      endif
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
