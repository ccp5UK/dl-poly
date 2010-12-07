      program acf3d

c***********************************************************************
c
c     DLPOLY routine for auto correlation function of three
c     dimensional vector -eg. velocities, forces.
c
c     Data is assumed to be formatted and compatible with the
c     HISTORY file written by DL_POLY subroutine traject.f
c
c     parallel version.
c     This program contains dummy routines for machine,gisum,and gdsum
c     to allow running on a single processor.
c
c     copyright daresbury laboratory 1993.
c
c     author t. forester  April 1993.
c
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c
c***********************************************************************


      parameter (mxcor = 300, mxatms = 1000, msatms = 500)

      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension ax(msatms,mxcor),ay(msatms,mxcor),az(msatms,mxcor)
      dimension acf(mxcor)
      dimension norm(mxcor)
      dimension buffer(mxcor)
      dimension ibuff(mxcor),list(msatms)

      character*80 dumpfile,outfile,headline
      character*8 acfme,atmnam(mxatms),oneline
c
c     determine number of nodes in use

      call machine(idnode,mxnode)
c
c     zero all arrays

      do i = 1,mxcor

         acf(i) = 0.d0
         norm(i) = 0

         do j = 1,msatms

            ax(j,i) =0.d0
            ay(j,i) =0.d0
            az(j,i) =0.d0

         enddo

      enddo
c
c     counter for correlation function index

      itcor = -1
c
c    interactive data

      write(*,*) 'dumpfile name ?'
      read(*,'(a80)') dumpfile

      write(*,*) 'number of atoms in file ?'
      read(*,*) natms
      if(natms.gt.mxatms) call error(30,idnode)

      write(*,*) 'atom label of interest ?'
      read(*,'(a8)') acfme

      write(*,*) 'time interval (in ps) for dumping ?'
      read(*,*) tstep

      write(*,*) 'number of steps for correlation function ?'
      read(*,*) ncor
      ncor = max(1,min(ncor,mxcor))

      write(*,*) 'output file'
      read(*,'(a80)') outfile

      write(*,*) '"1" for velocities, "2" for forces '
      read(*,*) lvf

c
c     open input and output files

      open(10,file=dumpfile,form='formatted')
      read(10,'(a80)') headline
      read(10,*) levcfg,imcon

      open(11,file=outfile)
      write(11,'(a80)') headline
      if(lvf.eq.1) write(11,*) ' correlating velocities'
      if(lvf.eq.2) write(11,*) ' correlating forces'

      if(levcfg.lt.lvf) call error(10,idnode)

c     start of correlation process loop ********************************

   99 continue
c
c     read in data

      read(10,'(a8)',end=100) oneline

      if(imcon.gt.0) then
         read(10,*) 
         read(10,*) 
         read(10,*) 
      endif

      if(lvf.eq.1) then

         if(levcfg.eq.1) then

            do i = 1,natms
               read(10,'(a8)') atmnam(i)
               read(10,*)
               read(10,*) vxx(i),vyy(i),vzz(i)
            enddo

         elseif(levcfg.eq.2) then

            do i = 1,natms
               read(10,'(a8)') atmnam(i)
               read(10,*)
               read(10,*) vxx(i),vyy(i),vzz(i)
               read(10,*)
            enddo

         endif

      else if(lvf.eq.2) then

         do i = 1,natms
            read(10,'(a8)') atmnam(i)
            read(10,*)
            read(10,*)
            read(10,*) vxx(i),vyy(i),vzz(i)
         enddo
      
      endif
c
c     counter for correlation function index

      itcor = itcor +1
      it = mod(itcor,mxcor) + 1
c
c     generate list of atoms of interest

      if(itcor.eq.0) then

         j = 0
         do i = idnode+1,natms,mxnode

            if(atmnam(i).eq.acfme) then

               j = j+1
               if(j.gt.msatms) call error(20,idnode)
               list(j) = i

            endif

         enddo
         
         jmax = j

         write(*,*) jmax ,' atoms of interest'
      endif
c
c     assign current variables to acf array

      do i = 1,jmax

         ii = list(i)

         ax(i,it) = vxx(ii)
         ay(i,it) = vyy(ii)
         az(i,it) = vzz(ii)

      enddo
c
c     correlate with previous time steps
c     outer loop over sites

      do i = 1,jmax

         ii = list(i)

         a0x = ax(i,it)
         a0y = ay(i,it)
         a0z = az(i,it)

c
c     number of steps to correlate

         itt = min(itcor,ncor)
c
c     inner loop over timesteps

         do j = 1,itt

            is = mod(it-j,ncor) + 1
            dotprd = ax(i,is)*a0x + ay(i,is)*a0y + az(i,is)*a0z
            acf(j) = acf(j) + dotprd
            norm(j) = norm(j) + 1

         enddo

      enddo

      goto 99
c
c     end of correlation process ***************************************

  100 continue
c
c     normalisation starts here

c     number of points in acf

      it = min(itcor,ncor)

c
c     sum arrays across nodes - replicated data strategy
      
      if(mxnode.gt.1) call gdsum(acf(1),ncor,buffer)
      if(mxnode.gt.1) call gisum(norm(1),ncor,ibuff)
c
c     normalise correlation functions
      if(acf(1).le.1d-20) acf(1) = 1.d0
      racf = dble(norm(1))/acf(1)
      if(racf.le.1d-20) racf= 1.d0

      acf(1) = 1.d0

      do i = 2,it
         norm(i) = max(1,norm(i))
         acf(i) = acf(i)/dble(norm(i))*racf

      enddo
c
c     write out results

      if(idnode.eq.0) then

         write(11,10) acfme,itcor+1,it,1.d0/racf
c
c     integral by Trapezium rule

         sum = 0.d0
         do i = 1,it
            
            if(i.gt.1) sum = sum+0.5d0*(acf(i)+acf(i-1))*tstep
            write(11,20) dble(i-1)*tstep,acf(i),sum

         enddo

      endif

   10 format(//,' Autocorrelation function calculated for :',a8,
     x     ' using',/,i20,' timesteps  ',/,i20,' bins ',//,4x,
     x     'normalisation constant : ',1p,e14.6,//,
     x     5x,'time(ps)',10x,'acf(i)',10x,'integral')
   20 format (f10.3,1p,6x,e14.6,6x,e14.6)

      write(*,*) 'done: output in ',outfile 
      end

      subroutine error(ierr,idnode)

   10 format(' error -- required information not in dump file')
   20 format(' error -- parameter msatms too small')
   30 format(' error -- parameter mxatms too small')

      if(idnode.eq.0) then

         if(ierr.eq.10) write(11,10)
         if(ierr.eq.20) write(11,20)
         if(ierr.eq.30) write(11,30)

      endif

      call exit()

      return
      end
      subroutine machine (idnode,mxnode)

      idnode = 0
      mxnode = 1

      return
      end
      subroutine gisum(i,j,k)
      return
      end
      subroutine gdsum(a,j,b)
      return
      end
