      program gself

c*********************************************************************
c
c     dl_poly program to calculate van Hove self correlation
c     function for selected atoms from dl_poly HISTORY file
c     includes gaussian analysis as described by rahman
c     in phys rev 136 (1964) a405
c
c     copyright daresbury laboratory 1996
c     author  w.smith jan 1996
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      parameter (mxatms=1080,nstr=240,mxslf=64,mxtim=64)
      parameter (pi=3.1415926536d0)
      character*1 fmt
      character*8 atname
      character*40 fname
      character*80 title
      character*8 atmnam(mxatms)
      dimension gslf(mxslf,mxtim),gslf0(nstr,mxtim,3)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension weight(mxatms),chge(mxatms)
      dimension avcell(9),cell(9),rcell(9),time(7)
      dimension acx(nstr),acy(nstr),acz(nstr)
      dimension xx1(nstr),yy1(nstr),zz1(nstr)
      dimension xx2(nstr),yy2(nstr),zz2(nstr)
      dimension order(0:10),alpha(0:10),ccn(0:10)
      dimension imd(mxtim),msm(mxtim)

c     open the i/o files

      open(5,file='gself_input')
      open(6,file='gself_output')

      write(6,'(a)')
     x  '# van Hove Self Correlation Function Program'

c     read the control variables

c     name of atomic species for autocorrelation
      read(5,*)atname

c     name of selected HISTORY file
      read(5,*)fname

c     HISTORY file formatted or unformatted
      read(5,*)fmt

c     total number of atoms in a HISTORY file configuration
      read(5,*)natms

c     time interval between successive configurations
      read(5,*)tstep

c     max number of configurations to sample
      read(5,*)nconf

c     sampling interval
      read(5,*)isampl

c     interval between origins
      read(5,*)iogslf

c     maximum correlation radius
      read(5,*)rcut

c     average cell vectors (taken from dl_poly OUTPUT file)
      read(5,*)avcell

      rcut2=rcut**2
      delr=rcut/dble(mxslf)

c     check on specified control variables

      ngslf=mxtim

      if(mod(ngslf,iogslf).ne.0)then

        ngslf=iogslf*(ngslf/iogslf)
        write(6,'(a,i5)')
     x    '# warning - gself array reset to ',ngslf

      endif

      nogslf=ngslf/iogslf

      write(6,'(a,a40)')'# name of target HISTORY file   : ',fname
      write(6,'(a,a8)')'# label  of atom  of interest   : ',atname
      write(6,'(a,i8)')'# total no. of atoms in config  : ',natms
      write(6,'(a,i8)')'# length of correlation arrays  : ',ngslf
      write(6,'(a,i8)')'# number of configurations      : ',nconf
      write(6,'(a,i8)')'# sampling interval             : ',isampl
      write(6,'(a,i8)')'# interval between origins      : ',iogslf
      write(6,'(a,f8.3)')'# time interval between configs : ',tstep
      write(6,'(a,f8.3)')'# selected correlation radius   : ',rcut
      write(6,'(a,f8.3)')'# radial correlation bin width  : ',delr
      write(6,'(a)')     '# average cell vectors:'
      write(6,'(a,3f12.6)')'# vector A :',avcell(1),avcell(2),avcell(3)
      write(6,'(a,3f12.6)')'# vector B :',avcell(4),avcell(5),avcell(6)
      write(6,'(a,3f12.6)')'# vector C :',avcell(7),avcell(8),avcell(9)

      if(natms.gt.mxatms)then

        write(6,'(a,2i5)')
     x    '# error - too many atoms in system'
        stop

      endif

c     initialise gself arrays

      do j=1,mxtim

        msm(j)=0

        do i=1,mxslf

          gslf(i,j)=0.d0

        enddo

      enddo

      do ipass=0,isampl-1

      lsr=0
      msr=0
      nsgslf=0
      do i=1,nstr

        acx(i)=0.d0
        acy(i)=0.d0
        acz(i)=0.d0

      enddo

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
     x      cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


        else

          call hread
     x      (fname,title,atmnam,iflg,imcon,keytrj,natms,nstep,timstp,
     x      cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


        endif

        if(iflg.lt.0)go to 100

        if(imcon.ge.1.and.imcon.le.3)then

          call invert(cell,rcell,det)

        else

          write(6,'(a)')
     x      '# error - incorrect boundary condition'
          stop

        endif

        j=0
        do i=1,natms

          if(atmnam(i).eq.atname)then

            j=j+1

            if(j.gt.nstr)then

              write(6,'(a)')
     x          '# error - too many atoms of specified type'
              stop

            endif
            xx2(j)=xxx(i)*rcell(1)+yyy(i)*rcell(4)+zzz(i)*rcell(7)
            yy2(j)=xxx(i)*rcell(2)+yyy(i)*rcell(5)+zzz(i)*rcell(8)
            zz2(j)=xxx(i)*rcell(3)+yyy(i)*rcell(6)+zzz(i)*rcell(9)

          endif

        enddo
        nat=j

        if(iconf.eq.0.and.ipass.eq.0)write(6,'(a,i8)')
     x    '# no. atoms of selected type    : ',nat

        if(iconf.gt.ipass)then

c     accumulate incremental distances

          do i=1,nat

            uuu=xx2(i)-xx1(i)
            vvv=yy2(i)-yy1(i)
            www=zz2(i)-zz1(i)

            uuu=uuu-nint(uuu)
            vvv=vvv-nint(vvv)
            www=www-nint(www)

            acx(i)=acx(i)+uuu*avcell(1)+vvv*avcell(4)+www*avcell(7)
            acy(i)=acy(i)+uuu*avcell(2)+vvv*avcell(5)+www*avcell(8)
            acz(i)=acz(i)+uuu*avcell(3)+vvv*avcell(6)+www*avcell(9)

          enddo

        endif

        do i=1,nat

          xx1(i)=xx2(i)
          yy1(i)=yy2(i)
          zz1(i)=zz2(i)

        enddo

c     calculate self correlation function

        if(iconf.gt.ipass)then

          if(mod(iconf,isampl).eq.ipass)then

            if(mod(nsgslf,iogslf).eq.0)then

              lsr=min(lsr+1,nogslf)
              msr=mod(msr,nogslf)+1
              imd(msr)=1
              do i=1,nat

                gslf0(i,msr,1)=0.d0
                gslf0(i,msr,2)=0.d0
                gslf0(i,msr,3)=0.d0

              enddo

            endif

            nsgslf=nsgslf+1

            do j=1,lsr

              m=imd(j)
              imd(j)=m+1
              msm(m)=msm(m)+1
              do i=1,nat

                rmsx=gslf0(i,j,1)+acx(i)
                rmsy=gslf0(i,j,2)+acy(i)
                rmsz=gslf0(i,j,3)+acz(i)

                rsq=rmsx**2+rmsy**2+rmsz**2

                gslf0(i,j,1)=rmsx
                gslf0(i,j,2)=rmsy
                gslf0(i,j,3)=rmsz

                if(rsq.lt.rcut2)then

                  k=int(sqrt(rsq)/delr)+1
                  gslf(k,m)=gslf(k,m)+1.d0

                endif

              enddo

            enddo

            do i=1,nat

              acx(i)=0.d0
              acy(i)=0.d0
              acz(i)=0.d0

            enddo

          endif

        endif

      enddo

  100 continue

      enddo

      if(iflg.lt.0)iconf=iconf-1
      last=min(nsgslf,ngslf)
      write(6,'(a,i6)')
     x  '# number of configurations read: ',iconf

c     normalise self correlation function

      do j=1,last

        do i=1,mxslf

          gslf(i,j)=(gslf(i,j)/(dble(msm(j))*dble(nat)))
     x      /(4.d0*pi*delr**3*((dble(i)-0.5d0)**2+1.d0/12.d0))

        enddo

      enddo

c     print out final self correlation function


      tstep=tstep*dble(isampl)

      do j=1,last,7

        do k=1,7

          time(k)=tstep*dble(k+j-1)

        enddo

        write(6,'(a10,7f10.5)')'# rad/time',(time(k),k=1,7)

        do i=1,mxslf

          rrr=delr*(dble(i)-0.5d0)
          write(6,'(8f10.5)')rrr,(gslf(i,k),k=j,min(j+6,last))

        enddo

        write(6,'(a)')'&'

      enddo

c     analysis of gaussian behaviour (after Rahman)

      ccn(0)=1.d0

      do i=1,9

        ccn(i)=ccn(i-1)*dble(2*i+1)/3.d0

      enddo

      do k=1,last

        order(0)=0.d0
        alpha(0)=0.d0
        do i=1,mxslf

          gslf(i,k)=4.d0*pi*gslf(i,k)*(delr*(dble(i)-0.5d0))**2

        enddo
        call simpson(mxslf,delr,order(0),gslf(1,k))
        order(0)=order(0)+0.25d0*delr*gslf(1,k)

        do j=1,9

          order(j)=0.d0

          do i=1,mxslf

            gslf(i,k)=gslf(i,k)*(delr*(dble(i)-0.5d0))**2

          enddo
          call simpson(mxslf,delr,order(j),gslf(1,k))
          order(j)=order(j)+0.25d0*delr*gslf(1,k)
          alpha(j)=order(j)/(ccn(j)*order(1)**j)-1.d0

        enddo

        write(6,'(a,f10.5)')'# moments for time',tstep*dble(k)
        write(6,'(a1,1p,5e12.4)')'#',(order(j),j=0,4)
        write(6,'(a1,1p,5e12.4)')'#',(order(j),j=5,9)
        write(6,'(a)')'# nongaussian terms alpha'
        write(6,'(a1,1p,5e12.4)')'#',(alpha(j),j=0,4)
        write(6,'(a1,1p,5e12.4)')'#',(alpha(j),j=5,9)

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

      read(nhist,'(a8,4i10,f12.6)',end=200)
     x    step,nstep,matms,ktrj,imcon,tstep

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

      character*80 cfgname
      character*40 history
      character*8 atmnam(*)

      dimension cell(9)
      dimension chge(*),weight(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension fxx(*),fyy(*),fzz(*)

      data nhist/77/

c     open the history file if new job

      if(iflg.eq.0)then

        close(nhist)
        open(nhist,file=history,form='unformatted',
     x       status='old',err=100)

        read(nhist,err=200) cfgname
        write(*,'(a,a)')'# History file header: ',cfgname
        read(nhist,end=200) datms
        if(natms.ne.nint(datms))then

          matms=dint(datms)
          write(*,'(a)')'# error - incorrect number of atoms in file'
          write(*,'(a,i6,a)')'# file contains',matms,' atoms'
          stop

        endif
        read(nhist,end=200) (atmnam(i),i=1,natms)
        read(nhist,end=200) (weight(i),i=1,natms)
        read(nhist,end=200) (chge(i),i=1,natms)

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
      subroutine simpson(n,del,sum,aaa)

c*********************************************************************
c
c     dl_poly subroutine for performing integration using
c     Simpson's rule.
c
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      dimension aaa(*)

      j=n
      m=2*(n/2)
      if(m.eq.n)j=n-1

      sum=(aaa(1)+aaa(j))/2.d0

      do i=2,j,2

        sum=sum+2.d0*aaa(i)+aaa(i+1)

      enddo

      sum=2.d0*del*sum/3.d0

      if(m.eq.n)sum=sum+del*(aaa(j)+aaa(n))/2.d0

      return
      end


