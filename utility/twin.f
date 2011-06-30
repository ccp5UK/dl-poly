c*********************************************************************
c
c     dl_poly utility for twinning metal clusters for potential of mean
c     force calculation.
c
c     author - w.smith july 2000
c     copyright daresbury laboratory 2000
c
c*********************************************************************

      parameter (mxatms=2000)

      implicit real*8(a-h,o-z)

      character*40 fname
      character*80 text

      character*8 name(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)

c     select CONFIG file

      write(*,*)'Enter name of required CONFIG file:'
      read(*,'(a)')fname
      write(*,'(a)')'Output file will be called: CONFIG.TWIN'

c     specify cluster separation

      write(*,*)'Specify cluster separation (A):'
      read(*,*)sep

c     open required files

      open(7,file=fname)
      open(8,file='CONFIG.TWIN')

c     read CONFIG file

      read(7,'(a80)')title
      read(7,*)levcfg,imcon
      if(levcfg.eq.0)then
        write(*,*)'Error - CONFIG file contains no velocities'
        stop
      endif

      if(imcon.gt.0)then
        write(*,*)
     x  'Warning - CONFIG file does not represent a cluster'
        read(7,*)
        read(7,*)
        read(7,*)
      endif
      do i=1,mxatms
        write(*,*)'===>',i
        read(7,'(a8)',end=100)name(i)
        read(7,*)xxx(i),yyy(i),zzz(i)
        read(7,*)vxx(i),vyy(i),vzz(i)
        if(levcfg.gt.1)read(7,*)a,b,c

      enddo

  100 close(7)
      natms=i-1
      write(*,*)'Number of atoms in cluster: ',natms

c     fix centre of mass

      cmx=0.d0
      cmy=0.d0
      cmz=0.d0

      do i=1,natms

        cmx=cmx+xxx(i)
        cmy=cmy+yyy(i)
        cmz=cmz+zzz(i)

      enddo

      cmx=cmx/dble(natms)
      cmy=cmy/dble(natms)
      cmz=cmz/dble(natms)

      m=0
      rmin=1.d6
      do i=1,natms

        xxx(i)=xxx(i)-cmx
        yyy(i)=yyy(i)-cmy
        zzz(i)=zzz(i)-cmz
        rr2=xxx(i)**2+yyy(i)**2+zzz(i)**2
        if(rr2.lt.rmin)then
          m=i
          rmin=rr2
        endif

      enddo

      write(*,*)'Central atom is: ',m,' at '
      write(*,*) xxx(m),yyy(m),zzz(m)

c     double up the cluster

      do i=1,natms

        name(i+natms)=name(i)
        xxx(i)=xxx(i)+0.5d0*sep
        xxx(i+natms)=xxx(i)-sep
        yyy(i+natms)=yyy(i)
        zzz(i+natms)=zzz(i)
        vxx(i+natms)=vxx(i)
        vyy(i+natms)=vyy(i)
        vzz(i+natms)=vzz(i)

      enddo

c     write twinned CONFIG file

      write(8,'(a80)')'TWIN: '//text(1:74)
      write(8,'(2i10)')1,0
      do i=1,2*natms

        write(8,'(a8)')name(i)
        write(8,'(3f20.12)')xxx(i),yyy(i),zzz(i)
        write(8,'(3f20.12)')vxx(i),vyy(i),vzz(i)

      enddo

c     close(8)

      end
