      program curden
c***********************************************************************
c
c     daresbury laboratory dl_poly program for the calculation
c     of the fourier transform of the particle current density
c     in space and time
c
c     original curden program written by w.smith march 1982
c     adapted for dl_poly by w.smith august 1994
c     parallel version
c
c     copyright daresbury laboratory 1994
c     author w smith 1994
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,msp=1080,iread=5,irite=6)
      parameter (pi=3.141592653589793d0,tpi=6.283185307179586d0)
      logical lden,lcor,ltim,lchg
      dimension xxx(msp),yyy(msp),zzz(msp)
      dimension vxx(msp),vyy(msp),vzz(msp)
      dimension man(10),space(mcore)
      data man/10*0/
      call timchk(0,wtime)

c     determine parallel machine characteristics

      open(iread,file='curden_input')
      open(irite,file='curden_output')
      call machine(idnode,mxnode)
      if(idnode.eq.0)then

         write(irite,'(a)')
     x   '# DL_POLY Current Density Correlation Program'

      endif

c     read in control variables

      call start
     x   (idnode,mxnode,natm,ncon,kmax,ntime,ngap,lden,
     x   lcor,ltim,lchg,tstep)

      klim=((2*kmax+1)**3-1)/2

c     calculate fourier transform of current density

      if(lden)then
        iman=(natm+mxnode-1)/mxnode
        man(1)=1
        man(2)=man(1)+2*iman
        man(3)=man(2)+2*iman
        man(4)=man(3)+2*iman
        man(5)=man(4)+4*klim
        man(6)=man(5)+8*klim
        man(7)=man(6)+8*klim
        man(8)=man(7)+2*iman*(kmax+1)
        man(9)=man(8)+2*iman*(kmax+1)
        man(10)=man(9)+2*iman*(kmax+1)
        call fourdc
     x     (idnode,mxnode,natm,ncon,kmax,lchg,xxx,yyy,zzz,vxx,vyy,vzz,
     x      space(man(1)),space(man(2)),space(man(3)),space(man(4)),
     x      space(man(5)),space(man(6)),space(man(7)),space(man(8)),
     X      space(man(9)))
      endif

c     calculate current density correlation function

      if(lcor)then
        nkt=(klim+mxnode-1)/mxnode
        man(1)=1
        man(2)=man(1)+ntime
        man(3)=man(2)+ntime
        man(4)=man(3)+3*klim
        man(5)=man(4)+8*nkt*ntime
        man(6)=man(5)+8*nkt*ntime
        man(7)=man(6)+8*klim
        call correl
     x     (idnode,mxnode,ncon,kmax,ntime,ngap,tstep,space(man(1)),
     x      space(man(2)),space(man(3)),space(man(4)),space(man(5)),
     x      space(man(6)))
      endif

c     calculate correlation function fourier transform

      if(ltim)then
        nkt=(klim+mxnode-1)/mxnode
        man(1)=1
        man(2)=man(1)+ntime
        man(3)=man(2)+3*klim
        man(4)=man(3)+ntime
        man(5)=man(4)+4*nkt*ntime
        man(6)=man(5)+4*ntime
        man(7)=man(6)+4*ntime
        call curfft
     x     (idnode,mxnode,kmax,ntime,ngap,tstep,space(man(1)),
     x      space(man(2)),space(man(3)),space(man(4)),space(man(5)),
     x      space(man(6)))
      endif
      call exit()
      end
      subroutine start
     x   (idnode,mxnode,natm,ncon,kmax,ntime,ngap,lden,lcor,
     x   ltim,lchg,tstep)

c***********************************************************************
c
c     read in control variables and check parameters
c
c     copyright daresbury laboratory 1994
c     author w smith 1994
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,msp=1080)
      parameter (iread=5,irite=6)
      logical lden,ltim,lcor,lchg,kill
      character*80 title
      kill=.false.
      read(iread,'(a80)')title
      read(iread,*)natm
      read(iread,*)ncon
      read(iread,*)kmax
      read(iread,*)lden
      read(iread,*)lcor
      read(iread,*)ltim
      read(iread,*)lchg
      read(iread,*)ntime
      read(iread,*)ngap
      read(iread,*)tstep
      if(idnode.eq.0)then
         write(irite,'(a,i10)')
     x     '# maximum size of dynamic core area         = ',mcore
         write(irite,'(a,i10)')
     x     '# number of atoms in simulation box         = ',natm
         write(irite,'(a,i10)')
     x     '# number of configurations                  = ',ncon
         write(irite,'(a,i10)')
     x     '# largest k vector index                    = ',kmax
         write(irite,'(a,i10)')
     x     '# length of correlation arrays              = ',ntime
         write(irite,'(a,i10)')
     x     '# interval between time origins             = ',ngap
         write(irite,'(a,1p,e12.4)')
     x     '# magnitude of time interval                = ',tstep
      endif
      if(idnode.eq.0)then

         if(lchg)then
            write(irite,'(a)')'# charge current option selected'
         else
            write(irite,'(a)')'# density current option selected'
         endif

      endif

c     check on control parameters

      if(natm.gt.msp)then

         if(idnode.eq.0)write(irite,'(a)')
     x        'error - number of atoms too large'
         kill=.true.

      endif
      if(tstep.le.0.d0) then

         if(idnode.eq.0)write(irite,'(a)')
     x        'error - time step assigned zero value'
         kill=.true.

      endif
      iman=(natm+mxnode-1)/mxnode
      klim=((2*kmax+1)**3-1)/2
      nkt=(klim+mxnode-1)/mxnode
      m1=0
      m2=0
      m3=0
      if(lden)m1=6*iman*(kmax+2)+24*klim
      if(lcor)m2=16*nkt*ntime+7*klim+2*ntime
      if(ltim)m3=ntime*(8*nkt+10)+3*klim
      kcore=max(m1,m2,m3)
      if(mcore.lt.kcore)then

         if(idnode.eq.0)write(irite,'(a,i8,a)')
     x        'error - insufficient core allocated. ',kcore,
     x        ' words required'
         kill=.true.

      endif

      if(nbits(ntime).ne.1)then
         if(idnode.eq.0)write(irite,'(a)')
     x        'error - parameter ntime not of the form 2**n'
         kill=.true.

      endif
      if(nbits(ngap).ne.1)then

         if(idnode.eq.0)write(irite,'(a)')
     x        'error - parameter ngap not of the form 2**m'
         kill=.true.

      endif

      if(kill)call exit()

      return
      end
      subroutine fourdc
     x     (idnode,mxnode,natm,ncon,kmax,lchg,xxx,yyy,zzz,
     x     vxx,vyy,vzz,elc,emc,enc,sxyz,vxyz,buff,ekx,eky,ekz)

c***********************************************************************
c
c     calculate spatial fourier transform of current density
c
c     copyright daresbury laboratory 1994
c     author w smith 1994
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (irite=6,iconf=7,iwork=8)
      parameter (pi=3.141592653589793d0,tpi=6.283185307179586d0)
      character*8 atname
      character*80 header
      logical lchg
      complex*16 ekx,eky,ekz,elc,emc,enc,vxyz,crdc
      complex*16 ckx,cky,ckz,buff
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension ekx(*),eky(*),ekz(*),buff(*)
      dimension elc(*),emc(*),enc(*),vxyz(*),sxyz(*)
      dimension cell(9),rcell(9)
      data cl/1.d0/
      rcl=tpi/cl
      rnatm=sqrt(1.d0/dble(natm))
      klim=((2*kmax+1)**3-1)/2

c     zero mean square density components

      do i=1,4*klim

         sxyz(i)=0.d0

      enddo

c     open history and density files

      open(iconf,file='HISTORY',status='old',err=200)
      open(iwork,file='CUR_SPC',form='unformatted')

      read(iconf,'(a80)')header
      if(idnode.eq.0)then
         write(irite,'(a)')'# HISTORY file header:'
         write(irite,'(a1,a80)')'#',header
      endif
      read(iconf,'(2i10)')levtrj,imcon

c     check file particulars

      if(imcon.eq.0)then
        if(idnode.eq.0)
     x       write(irite,'(a)')
     x       'error - no periodic boundary in HISTORY file'
        call exit()
      endif
      if(levtrj.lt.1)then
        if(idnode.eq.0)
     x       write(irite,'(a)')
     x       'error - no velocity data in HISTORY file'
        call exit()
      endif

      do nstp=1,ncon

        read(iconf,*,end=100)

c     read in configuration data

        read(iconf,*)cell(1),cell(2),cell(3)
        read(iconf,*)cell(4),cell(5),cell(6)
        read(iconf,*)cell(7),cell(8),cell(9)

c     calculate reciprocal lattice basis

        call invert(cell,rcell,det)

c     read atomic data

        do i=1,natm

           if(lchg)then
              read(iconf,*)atname,idx,mass,chg
           else
              read(iconf,*)atname
           endif
          read(iconf,*)xx,yy,zz
          read(iconf,*)vxx(i),vyy(i),vzz(i)
          if(levtrj.eq.2)read(iconf,*)fxx,fyy,fzz
          xxx(i)=xx*rcell(1)+yy*rcell(4)+zz*rcell(7)
          yyy(i)=xx*rcell(2)+yy*rcell(5)+zz*rcell(8)
          zzz(i)=xx*rcell(3)+yy*rcell(6)+zz*rcell(9)
          if(lchg)then
             vxx(i)=chg*vxx(i)
             vyy(i)=chg*vyy(i)
             vzz(i)=chg*vzz(i)
          endif

        enddo

c     calculate fourier exponential terms

        m=0
        do i=idnode+1,natm,mxnode
          m=m+1
          ekx(m)=(1.d0,0.d0)
          eky(m)=(1.d0,0.d0)
          ekz(m)=(1.d0,0.d0)
          elc(m)=cmplx(cos(rcl*xxx(i)),sin(rcl*xxx(i)))
          emc(m)=cmplx(cos(rcl*yyy(i)),sin(rcl*yyy(i)))
          enc(m)=cmplx(cos(rcl*zzz(i)),sin(rcl*zzz(i)))
        enddo
        matm=m
        do l=1,kmax
          do i=1,matm
            ekx(i+l*matm)=ekx(i+matm*(l-1))*elc(i)
            eky(i+l*matm)=eky(i+matm*(l-1))*emc(i)
            ekz(i+l*matm)=ekz(i+matm*(l-1))*enc(i)
          enddo
        enddo

c     start loop over k vectors

        kkk=0
        mmin=0
        lmin=1
        do n=0,kmax
          rn=dble(n)
          do i=1,matm
            enc(i)= ekz(i+n*matm)
          enddo
          do m=mmin,kmax
            rm=dble(m)
            if(m.ge.0)then
              do i=1,matm
                emc(i)= eky(i+m*matm)*enc(i)
              enddo
            else
              do i=1,matm
                emc(i)= conjg(eky(i-m*matm))*enc(i)
              enddo
            endif
            do l=lmin,kmax
              rl=dble(l)
              if(l.ge.0)then
                do i=1,matm
                  elc(i)= ekx(i+l*matm)*emc(i)
                enddo
              else
                do i=1,matm
                  elc(i)= conjg(ekx(i-l*matm))*emc(i)
                enddo
              endif

              rkx=rcell(1)*rl+rcell(2)*rm+rcell(3)*rn
              rky=rcell(4)*rl+rcell(5)*rm+rcell(6)*rn
              rkz=rcell(7)*rl+rcell(8)*rm+rcell(9)*rn
              cnm=rnatm/sqrt(rkx*rkx+rky*rky+rkz*rkz)
              ckx= cnm*conjg(crdc(matm,vxx(idnode+1),mxnode,elc,1))
              cky= cnm*conjg(crdc(matm,vyy(idnode+1),mxnode,elc,1))
              ckz= cnm*conjg(crdc(matm,vzz(idnode+1),mxnode,elc,1))
              vxyz(1+4*kkk)=rkx*ckx+rky*cky+rkz*ckz
              vxyz(2+4*kkk)=rky*ckz-rkz*cky
              vxyz(3+4*kkk)=rkz*ckx-rkx*ckz
              vxyz(4+4*kkk)=rkx*cky-rky*ckx
              kkk=kkk+1

            enddo
            lmin=-kmax
          enddo
          mmin=-kmax
        enddo

c     global sum of fourier transform of current density

        call gdsum(vxyz,8*klim,buff)

c     accumulate mean square values of current density components

        do i=1,4*klim

           sxyz(i)=sxyz(i)+real(vxyz(i)*conjg(vxyz(i)))

        enddo

c     store fourier transform of current density

        if(idnode.eq.0)write(iwork)(vxyz(j),j=1,4*klim)

      enddo

  100 continue

      nstp=nstp-1

      if(idnode.eq.0)then
         write(irite,'(a,i10)')
     x     '# number of timesteps in simulation data    = ',nstp
         close (iwork)
         write(irite,'(a)')
     x     '# CUR_SPC file written and closed'
         write(irite,'(a,a,a)')'#    L    M    N',
     x     '   |Cl(k)|^2','   |Ct(k)|^2'
         k=0
         mmin=0
         lmin=1
         do n=0,kmax
            do m=mmin,kmax
               do l=lmin,kmax
                  slong=sxyz(k+1)/dble(nstp)
                  stran=(sxyz(k+2)+sxyz(k+3)+sxyz(k+4))/dble(nstp)
                  write(irite,'(a1,3i5,1p,2e12.4)')"#",l,m,n,
     x               slong,stran
                  k=k+4
               enddo
               lmin=-kmax
            enddo
            mmin=-kmax
         enddo
         close (iwork)
         write(irite,'(a)')
     x     '# CUR_DEN file written and closed'
      endif

      close (iconf)
      call timchk(1,wtime)

      return

  200 continue

      if(idnode.eq.0)
     x   write(irite,'(a)')'error - HISTORY file not found'
      call exit()

      end
      subroutine correl
     x  (idnode,mxnode,ncon,kmax,ntime,ngap,tstep,ind,num,kvec,cfkt,
     x   ckr0,ckr)

c***********************************************************************
c
c     calculate current density correlation function
c
c     copyright daresbury laboratory 1994
c     author w smith 1994
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      parameter (irite=6,iwork=8,isave1=19)
      complex*16 ckr,ckr0,cfkt
      dimension cfkt(*),ckr0(*),ckr(*),ind(*),num(*),kvec(3,*)

c     open files isave1 and iwork

      open(iwork,file='CUR_SPC',form='unformatted',status='old',err=200)

      if(idnode.eq.0)then

         open(isave1,file='CUR_COR',form='unformatted')

      endif

c     set control parameters

      norg=ntime/ngap
      klim=((2*kmax+1)**3-1)/2
      nblk=mod(klim,mxnode)
      if(idnode.lt.nblk)then
         iblk=4*(klim/mxnode+1)
         ibgn=idnode*iblk
      else
         iblk=4*(klim/mxnode)
         ibgn=idnode*iblk+4*nblk
      endif

      lor=0
      mor=0

c     initialise cfkt array


      do l=1,iblk*ntime
        cfkt(l)=(0.d0,0.d0)
      enddo

      do l=1,ntime
        num(l)=0
      enddo

c     start of loop over time steps

      do n=0,ncon

        read(iwork,end=100)(ckr(j),j=1,4*klim)

        if(mod(n,ngap).le.0)then

          lor=min(lor+1,norg)
          mor=mod(mor,norg)+1
          ind(mor)=1

          do k=1,iblk
            ckr0(mor+ntime*(k-1))=conjg(ckr(k+ibgn))
          enddo

        endif

        do l=1,lor

          m=ind(l)
          ind(l)=m+1
          num(m)=num(m)+1
          do k=1,iblk
            cfkt(m+ntime*(k-1))=cfkt(m+ntime*(k-1))+
     x                          ckr(k+ibgn)*ckr0(l+ntime*(k-1))
          enddo

        enddo

      enddo
  100 continue

      last=min(n-1,ntime)

c     normalise correlation functions, combine transverse components

      i=0
      do k=1,iblk,4

        rnorm1=real(cfkt(1+(k-1)*ntime))
        if(rnorm1.gt.0.d0)rnorm1=dble(num(1))/rnorm1
        rnorm2=real(cfkt(1+k*ntime)+cfkt(1+(k+1)*ntime)+
     x    cfkt(1+(k+2)*ntime))
        if(rnorm2.gt.0.d0)rnorm2=dble(num(1))/rnorm2

        do l=1,last
            cfkt(l+i*ntime)=rnorm1*cfkt(l+(k-1)*ntime)/dble(num(l))
            cfkt(l+(i+1)*ntime)=rnorm2*(cfkt(l+k*ntime)+
     x        cfkt(l+(k+1)*ntime)+cfkt(l+(k+2)*ntime))/dble(num(l))
        enddo
        i=i+2

      enddo

c     calculate k vector indices for printing

      k=0
      mmin=0
      lmin=1
      do n=0,kmax
         do m=mmin,kmax
            do l=lmin,kmax
               k=k+1
               kvec(1,k)=l
               kvec(2,k)=m
               kvec(3,k)=n
            enddo
            lmin=-kmax
         enddo
         mmin=-kmax
      enddo

c     store correlation functions in disc file

      iblk=iblk/2

      if(idnode.eq.0)then

         write(irite,'(a)')
     x   '# Correlation Function of FT of Current Density'
         write(isave1)(cfkt(j),j=1,ntime*iblk)
         call tabulate(idnode,mxnode,kmax,ntime,tstep,kvec,cfkt)

         do kdnode=1,mxnode-1

            call csend(1744549,kdnode,4,kdnode,0)
            call crecv(1744550,cfkt,16*ntime*iblk)
            write(isave1)(cfkt(j),j=1,ntime*iblk)
            call tabulate(kdnode,mxnode,kmax,ntime,tstep,kvec,cfkt)

         enddo

      else

         call crecv(1744549,kdnode,4)
         call csend(1744550,cfkt,16*ntime*iblk,0,0)

      endif

c     close files

      if(idnode.eq.0)then
         close (isave1)
         write(irite,'(a)')
     x      '# CUR_COR file written and closed'
      endif

      close (iwork)
      call timchk(1,wtime)

      return

  200 continue

      if(idnode.eq.0)
     x   write(irite,'(a)')'error - CUR_SPC file not found'
      call exit

      end
      subroutine curfft
     x (idnode,mxnode,kmax,ntime,ngap,tstep,key,kvec,wind,cfkt,work,fta)
c
c***********************************************************************
c
c    calculate fourier transform of current density correlation function
c
c     copyright daresbury laboratory 1994
c     author w smith 1994
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (irite=6,isave1=19,isave2=20)
      parameter (pi=3.141592653589793d0,tpi=6.283185307179586d0)
      complex*16 cfkt(*),work(*),fta(*)
      dimension key(*),wind(*),kvec(3,*)
      data a0,a1,a2/0.42d0,0.50d0,0.08d0/

c     set control parameters

      ntime2=2*ntime
      norg=ntime/ngap
      klim=((2*kmax+1)**3-1)/2
      nblk=mod(klim,mxnode)
      if(idnode.lt.nblk)then
         iblk=2*(klim/mxnode+1)
         ibgn=idnode*iblk
      else
         iblk=2*(klim/mxnode)
         ibgn=idnode*iblk+2*nblk
      endif
      omega=tpi/(dble(ntime2)*tstep)

c     open files isave1 and isave2

      open(isave1,file='CUR_COR',form='unformatted',status='old',
     x    err=200)

      if(idnode.eq.0)then

         open(isave2,file='CUR_FFT',form='unformatted')

      endif

c     initialise complex fast fourier transform routine

      ind=1
      isw=-1
      call fft(ind,isw,ntime2,key,cfkt,work,cfkt)
      ind=0

c     set up window function (blackman function)

      arg=tpi/dble(ntime2)
      do i=1,ntime

         ccc=cos(arg*dble(i+ntime-1))
         wind(i)=a0-a1*ccc+a2*(2.d0*ccc**2-1.d0)

      enddo

c     read the correlation functions from disc

      do i=0,mxnode-1

         if(i.eq.idnode)then

            read(isave1)(cfkt(j),j=1,ntime*iblk)

         else

            read(isave1)

         endif

      enddo

c     loop over correlation functions (longitudinal and transverse)

      do k=1,iblk

c     apply window function

         do j=1,ntime2

            if(j.le.ntime)then

               fta(j)=cfkt(j+ntime*(k-1))*wind(j)

            else

               fta(j)=(0.d0,0.d0)

            endif

         enddo

         fta(1)=fta(1)/2.d0

c     apply complex fourier transform

         call fft(ind,isw,ntime2,key,fta,work,fta)

c     store fourier coefficients

         do j=1,ntime

            cfkt(j+ntime*(k-1))=fta(j)

         enddo

      enddo

c     calculate k vector indices for printing

      k=0
      mmin=0
      lmin=1
      do n=0,kmax
         do m=mmin,kmax
            do l=lmin,kmax
               k=k+1
               kvec(1,k)=l
               kvec(2,k)=m
               kvec(3,k)=n
            enddo
            lmin=-kmax
         enddo
         mmin=-kmax
      enddo

c     save fourier coefficients

      if(idnode.eq.0)then

         write(irite,'(a)')
     x   '# Fourier Transform of Current Density Correlation Functions'
         write(isave2)(cfkt(j),j=1,ntime*iblk)
         call tabulate(idnode,mxnode,kmax,ntime,omega,kvec,cfkt)

         do kdnode=1,mxnode-1

            call csend(1534549,kdnode,4,kdnode,0)
            call crecv(1534550,cfkt,16*ntime*iblk)
            write(isave2)(cfkt(j),j=1,ntime*iblk)
            call tabulate(kdnode,mxnode,kmax,ntime,omega,kvec,cfkt)

         enddo

      else

         call crecv(1534549,kdnode,4)
         call csend(1534550,cfkt,16*ntime*iblk,0,0)

      endif

c     close files

      close (isave1)
      close (isave2)
      if(idnode.eq.0)write(irite,'(a)')
     x   '# CUR_FFT file written and closed'
      call timchk(1,wtime)

      return

  200 continue

      if(idnode.eq.0)
     x   write(irite,'(a)')'error - CUR_COR file not found'
      call exit()

      end
      complex*16 function crdc(n,a,i,b,j)
c
c***********************************************************************
c
c     dl_poly utility routine for dot product of real array A and
c     complex array B
c
c     copyright daresbury laboratory 1994
c     author w smith
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      dimension a(n),b(2*n)
      u=sdot(n,a(1),i,b(1),2*j)
      v=sdot(n,a(1),i,b(2),2*j)
      crdc=cmplx(u,v)
      return
      end
      double precision function sdot(n,a,i,b,j)
c
c***********************************************************************
c
c     dl_poly utility routine for dot product of real array A and
c     real array B
c
c     copyright daresbury laboratory 1994
c     author w smith
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)

      dimension a(*),b(*)

      l=1
      m=1
      sdot=0.d0
      do k=1,n

         sdot=sdot+a(l)*b(m)

         l=l+i
         m=m+j

      enddo

      return
      end
      function nbits(n)
c
c**********************************************************************
c
c     dl_poly utility for counting number of set bits in an integer N
c
c     copyright daresbury laboratory 1994
c     author w smith
c
c**********************************************************************
c
      m=n
      nbits=0

      do i=1,128

         if(m.eq.0)return

         nbits=nbits+m-2*(m/2)
         m=m/2

      enddo

      return
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
      subroutine fft(ind,isw,ndiv,key,aaa,wfft,bbb)
c***********************************************************************
c
c     fast fourier transform routine
c
c     copyright daresbury laboratory 1994
c     author w smith
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      logical check
      complex*16 aaa(ndiv),bbb(ndiv),wfft(ndiv),ttt
      dimension key(ndiv)
      data tpi/6.2831853072d0/
   10 format(1h0,'error - number of points not a power of two')

c
c     check that array is of suitable length
      nt=1
      check=.true.
      do i=1,20
         nt=2*nt
         if(nt.eq.ndiv)then
            check=.false.
            nu=i
         endif
      enddo
      if(check)then
         write(*,10)
         call exit()
      endif

      if(ind.gt.0)then
c
c     set reverse bit address array

         do kkk=1,ndiv
            iii=0
            jjj=kkk-1
            do j=1,nu
               jj2=jjj/2
               iii=2*(iii-jj2)+jjj
               jjj=jj2
            enddo
            key(kkk)=iii+1
         enddo
c
c     initialise complex exponential factors

         tpn=tpi/dble(ndiv)
         arg=0.d0
         np1=ndiv+1
         np2=ndiv/2
         wfft(1)=(1.d0,0.d0)
         do i=1,np2
            arg=tpn*dble(i)
            wfft(i+1)=cmplx(cos(arg),sin(arg))
            wfft(np1-i)=conjg(wfft(i+1))
         enddo

         return
      endif

c
c     take conjugate of exponentials if required

      if(isw.lt.0)then

         do i=1,ndiv
            wfft(i)=conjg(wfft(i))
         enddo

      endif

c
c     take copy input array

      do i=1,ndiv
         bbb(i)=aaa(i)
      enddo
c
c     perform fourier transform

      kkk=0
      nu1=nu-1
      np2=ndiv/2
      do l=1,nu

  100    do i=1,np2
            iii=key(kkk/2**nu1+1)
            kk1=kkk+1
            k12=kk1+np2
            ttt=bbb(k12)*wfft(iii)
            bbb(k12)=bbb(kk1)-ttt
            bbb(kk1)=bbb(kk1)+ttt
            kkk=kkk+1
         enddo
         kkk=kkk+np2
         if(kkk.lt.ndiv)go to 100
         kkk=0
         nu1=nu1-1
         np2=np2/2

      enddo
c
c     unscramble the fft using bit address array

      do kkk=1,ndiv
         iii=key(kkk)
         if(iii.gt.kkk)then
            ttt=bbb(kkk)
            bbb(kkk)=bbb(iii)
            bbb(iii)=ttt
         endif
      enddo
c
c     restore exponentials to unconjugated values if necessary

      if(isw.lt.0)then

         do i=1,ndiv
            wfft(i)=conjg(wfft(i))
         enddo

      endif

      return
      end
      subroutine tabulate(idnode,mxnode,kmax,ntime,tstep,kvec,cfkt)
c
c*********************************************************************
c
c     routine to write out blocked arrays of complex numbers
c
c     copyright daresbury laboratory
c     author w smith 1994
c
c*********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (irite=6)
      complex*16 cfkt(*)
      dimension kvec(3,*)

      klim=((2*kmax+1)**3-1)/2
      nblk=mod(klim,mxnode)
      if(idnode.lt.nblk)then
         iblk=2*(klim/mxnode+1)
         ibgn=idnode*iblk
      else
         iblk=2*(klim/mxnode)
         ibgn=idnode*iblk+2*nblk
      endif

      do k=1,iblk

         kk=(ibgn+k-1)/2+1
         if(mod(k,2).eq.1)then
            write(irite,'(a,3i5)')
     x     '# Longitudinal Component ',kvec(1,kk),kvec(2,kk),kvec(3,kk)
         else
            write(irite,'(a,3i5)')
     x     '# Transverse  Component  ',kvec(1,kk),kvec(2,kk),kvec(3,kk)
         endif
         do j=1,ntime

            ttt=tstep*dble(j-1)
            write(irite,'(1p,3e12.4)')ttt,cfkt(j+(k-1)*ntime)

         enddo

         write(irite,'(a)')'&'

      enddo

      return
      end
      subroutine machine(idnode,mxnode)
c
c*********************************************************************
c
c     dl_poly subroutine for obtaining charcteristics of
c     the computer on which the program is being run
c
c     this version is for the intel hypercube
c
c     copyright daresbury laboratory 1992
c     author - w.smith july 1992
c
c*********************************************************************
c

      implicit real*8(a-h,o-z)

c
c     number of nodes available to program

      mxnode=numnodes()

c
c     identity of executing node

      idnode=mynode()

      return
      end
      subroutine timchk(ktim,time)
c
c***********************************************************************
c
c     timing routine (time elapsed in seconds)
c
c     copyright daresbury laboratory 1994
c     author w. smith
c
c***********************************************************************
c
      parameter (irite=6)
      real*8 time
      save init
      data init/0/
      idnode=mynode()
c
c     synchronise time
      call gsync()
      if(init.eq.0)init=mclock()
      itime=mclock()-init
      time=dble(itime)/1000.d0
      if(ktim.gt.0.and.idnode.eq.0)write(irite,
     x  "('# time elapsed since job start = ',f15.8,' seconds',/)")
     x  time
      return
      end
