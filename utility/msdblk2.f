      program msdblk
c*********************************************************************
c
c     dl_poly program to calculate mean square displacement
c     for selected atoms from dl_poly HISTORY file
c     this version incorporates a blocking average method
c     to calculate the statistical error
c
c     copyright daresbury laboratory 1996
c     author  w.smith jan 1996
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      parameter (mxatms=1080,npts=600,nstr=1080,ndeg=10,mxconf=2500)
      logical safe
      character*1 fmt
      character*8 atname
      character*40 fname
      character*80 title
      character*8 atmnam(mxatms)
      real*8 msd(npts,mxconf),msd0(npts,nstr,3),sig(ndeg)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension weight(mxatms),chge(mxatms)
      dimension avcell(9),cell(9),rcell(9)
      dimension xx0(nstr),yy0(nstr),zz0(nstr)
      dimension xx1(nstr),yy1(nstr),zz1(nstr)
      dimension acx(nstr),acy(nstr),acz(nstr)
      dimension imd(npts),msm(npts)

c     open the I/O files

      open(5,file='msd_input')
c      open(6,file='msd_output')

      write(6,'(a)')
     x  '# Mean Square Displacement Blocking Program'

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

c     working length of msd array
      read(5,*)nmsd

c     sampling interval
      read(5,*)isampl

c     interval between origins
      read(5,*)iomsd

c     average cell vectors (taken from dl_poly OUTPUT file)
      read(5,*)avcell

      nconf=min(nconf,mxconf*isampl)

c     check on specified control variables

      write(6,'(a,a40)')'# name of target HISTORY file   : ',fname
      write(6,'(a,a8)')'# label  of atom  of interest   : ',atname
      write(6,'(a,i8)')'# total no. of atoms in config  : ',natms
      write(6,'(a,i8)')'# length of correlation arrays  : ',nmsd
      write(6,'(a,i8)')'# number of configurations      : ',nconf
      write(6,'(a,i8)')'# sampling interval             : ',isampl
      write(6,'(a,i8)')'# interval between origins      : ',iomsd
      write(6,'(a,f8.3)')'# time interval between configs : ',tstep
      write(6,'(a)')     '# average cell vectors:'
      write(6,'(a,3f12.6)')'# vector A :',avcell(1),avcell(2),avcell(3)
      write(6,'(a,3f12.6)')'# vector B :',avcell(4),avcell(5),avcell(6)
      write(6,'(a,3f12.6)')'# vector C :',avcell(7),avcell(8),avcell(9)

      if(natms.gt.mxatms)then

        write(6,'(a,2i5)')
     x    'error - too many atoms in system'
        stop

      endif
      if(nmsd.gt.npts)then

        write(6,'(a,2i5)')
     x    'error - msd array set too large',nmsd,npts
        stop

      endif
      if(mod(nmsd,iomsd).ne.0)then

        nmsd=iomsd*(nmsd/iomsd)
        write(6,'(a,i5)')
     x    '# warning - msd array reset to ',nmsd

      endif

      nomsd=nmsd/iomsd

c     initialise msd arrays

      lsr=0
      msr=0
      nsmsd=0
      do i=1,nstr

        acx(i)=0.d0
        acy(i)=0.d0
        acz(i)=0.d0

      enddo
      do j=1,npts
        msm(j)=0
        do i=1,nstr

          msd0(j,i,1)=0.d0
          msd0(j,i,2)=0.d0
          msd0(j,i,3)=0.d0

        enddo
        do i=1,mxconf

          msd(j,i)=0.d0

        enddo
      enddo

c     set default cell properties

      do i=1,9
        cell(i)=0.d0
      enddo
      cell(1)=1.d0
      cell(5)=1.d0
      cell(9)=1.d0

      keytrj=0

      do iconf=0,nconf-1

        if(fmt.eq."U".or.fmt.eq."u")then

          call uread
     x      (fname,title,atmnam,imcon,keytrj,natms,nstep,cell,chge,
     x      weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


        else

          call hread
     x      (fname,title,atmnam,imcon,keytrj,natms,nstep,cell,chge,
     x      weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


        endif

        if(imcon.ge.1.and.imcon.le.3)then

          call invert(cell,rcell,det)

        else

          write(6,'(a)')
     x      'error - incorrect boundary condition'
          stop

        endif

        j=0
        do i=1,natms

          if(atmnam(i).eq.atname)then

            j=j+1

            if(j.gt.nstr)then

              write(6,'(a)')
     x          'error - too many atoms of specified type'
              stop

            endif
            xx1(j)=xxx(i)*rcell(1)+yyy(i)*rcell(4)+zzz(i)*rcell(7)
            yy1(j)=xxx(i)*rcell(2)+yyy(i)*rcell(5)+zzz(i)*rcell(8)
            zz1(j)=xxx(i)*rcell(3)+yyy(i)*rcell(6)+zzz(i)*rcell(9)

          endif

        enddo
        nat=j

        if(iconf.gt.0)then

c     accumulate incremental distances

          do i=1,nat

            uuu=xx1(i)-xx0(i)
            vvv=yy1(i)-yy0(i)
            www=zz1(i)-zz0(i)

            uuu=uuu-nint(uuu)
            vvv=vvv-nint(vvv)
            www=www-nint(www)

            acx(i)=acx(i)+uuu*avcell(1)+vvv*avcell(4)+www*avcell(7)
            acy(i)=acy(i)+uuu*avcell(2)+vvv*avcell(5)+www*avcell(8)
            acz(i)=acz(i)+uuu*avcell(3)+vvv*avcell(6)+www*avcell(9)

          enddo

        endif

        do i=1,nat

          xx0(i)=xx1(i)
          yy0(i)=yy1(i)
          zz0(i)=zz1(i)

        enddo

c     calculate mean square displacement

        if(iconf.gt.0)then

          if(mod(iconf,isampl).eq.0)then

            if(mod(nsmsd,iomsd).eq.0)then

              lsr=min(lsr+1,nomsd)
              msr=mod(msr,nomsd)+1
              imd(msr)=1
              do i=1,nat

                msd0(msr,i,1)=0.d0
                msd0(msr,i,2)=0.d0
                msd0(msr,i,3)=0.d0

              enddo

            endif

            nsmsd=nsmsd+1

            do j=1,lsr

              m=imd(j)
              imd(j)=m+1
              msm(m)=msm(m)+1
              k=msm(m)
              do i=1,nat

                rmsx=msd0(j,i,1)+acx(i)
                rmsy=msd0(j,i,2)+acy(i)
                rmsz=msd0(j,i,3)+acz(i)
                msd(m,k)=msd(m,k)+rmsx**2+rmsy**2+rmsz**2
                msd0(j,i,1)=rmsx
                msd0(j,i,2)=rmsy
                msd0(j,i,3)=rmsz

              enddo

            enddo

            do i=1,nstr

              acx(i)=0.d0
              acy(i)=0.d0
              acz(i)=0.d0

            enddo

          endif

        endif

      enddo

      last=min(nsmsd,nmsd)
      write(6,'(a,i6)')
     x  '# number of configurations sampled: ',iconf

c     normalise mean square displacement

      rnat=1.d0/dble(nat)
      do j=1,last

        do i=1,msm(j)

          msd(j,i)=rnat*msd(j,i)

        enddo

      enddo

c     start of blocking statistics

      do j=1,last

c     calculate average value at time t

        avg=0.d0

        do i=1,msm(j)

          avg=avg+msd(j,i)

        enddo

        avg=avg/dble(msm(j))

c     calculate estimates of sigma

        big=0.d0
        mmm=msm(j)

        do n=1,ndeg

          if(mmm.gt.1)then

            sig(n)=0.d0

            do i=1,mmm

              sig(n)=sig(n)+(msd(j,i)-avg)**2

            enddo

            sig(n)=sqrt(sig(n)/(dble(mmm)*dble(mmm-1)))
            big=max(big,sig(n))
            mmm=mmm/2

            k=0
            do i=1,mmm

              k=k+2
              msd(j,i)=0.5d0*(msd(j,k-1)+msd(j,k))

            enddo

          endif
          eno=dble((iconf/isampl)/j)
          bug=avg*sqrt(2.d0*rnat/(3.d0*eno))
        enddo

c     write out results for time t

        time=tstep*dble(isampl)*dble(j)
        write(6,'(1p,4e11.4,i10)')time,avg,big,bug,msm(j)

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
     x  (history,cfgname,atmnam,imcon,keytrj,natms,nstep,cell,chge,
     x   weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

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

c     open history file if new job

      if(new)then

        open(nhist,file=history,status='old', err=100)

        read(nhist,'(a80)',err=200) cfgname
        write(*,'(a,a)')'#History file header: ',cfgname
        read(nhist,'(2i10)',end=200) ktrj,imcon
        if(keytrj.gt.ktrj)then

          if(ktrj.eq.0)write(*,'(a)')'error - no velocities in file'
          if(keytrj.gt.1)write(*,'(a)')'error - no forces in file'
          stop

        endif

        new=.false.

      endif

      read(nhist,'(a8,4i10)',end=200) step,nstep,matms,ktrj,imcon

      if(natms.ne.matms)then

        write(*,'(a)')'error - incorrect number of atoms in file'
        write(*,'(i6)')'file contains',matms,' atoms'
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

      return

  100 continue

      write(*,'(a)')'error - History file not found'
      stop

  200 continue
      write(*,'(a)')'warning - end of History file encountered'
      close (nhist)

      return
      end
      subroutine uread
     x  (history,cfgname,atmnam,imcon,keytrj,natms,nstep,cell,chge,
     x  weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

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

c     open the history file if new job

      if(new)then

        open(nhist,file=history,form='unformatted',
     x       status='old',err=100)

        read(nhist,err=200) cfgname
        write(*,'(a,a)')'#History file header: ',cfgname
        read(nhist,end=200) datms
        if(natms.ne.nint(datms))then

          write(*,'(a)')'error - incorrect number of atoms in file'
          write(*,'(i6)')'file contains',matms,' atoms'
          stop

        endif
        read(nhist,end=200) (atmnam(i),i=1,natms)
        read(nhist,end=200) (weight(i),i=1,natms)
        read(nhist,end=200) (chge(i),i=1,natms)

        new=.false.

      endif

      read(nhist,end=200)dstep,datms,trjkey,dimcon
      nstep=nint(dstep)
      ktrj=nint(trjkey)
      imcon=nint(dimcon)
      if(keytrj.gt.ktrj)then

        if(ktrj.eq.0)write(*,'(a)')'error - no velocities in file'
        if(keytrj.gt.1)write(*,'(a)')'error - no forces in file'
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

      return

  100 continue

      write(*,'(a)')'error - History file not found'
      stop

  200 continue
      write(*,'(a)')'warning - end of History file encountered'
      close (nhist)

      return
      end


