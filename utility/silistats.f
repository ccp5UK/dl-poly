      program silistats
c*********************************************************************
c
c     dl_poly program to calculate silicate bond statistics
c     in a HISTORY file
c
c     copyright daresbury laboratory march 1995
c     author w.smith march 1995
c
c*********************************************************************
c
      implicit real*8 (a-h,o-z)

      parameter (mxatms=1080,mxlist=60,mxrdf=128,mxtp=4,mxshl=30)

      character*8 name,atmnam(mxatms)
      character*80 text
      character*40 fname

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),nbr(0:4)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension chge(mxatms),weight(mxatms)
      dimension lentry(mxatms),list(mxatms,mxlist),nqr(0:4),nal(2)
      dimension ltype(mxatms),kode(mxatms,mxtp),ihist(2000,mxtp)
      dimension key(2000,mxtp),num(mxtp),radc(mxtp,mxtp),ind(mxtp)
      dimension cell(9),rcell(9),nspec(mxtp),nshell(mxshl,mxtp)

      open (unit=5,file='bond_input')
      open (unit=6,file='bond_output')

      read(5,*)fname
      read(5,*)natms
      if(natms.gt.mxatms)then

        write(6,'(a)')'error -too many atoms'
        stop

      endif
      read(5,*)ntp
      if(ntp.gt.mxtp)then

        write(6,'(a)')'error -too many species types'
        stop

      endif
      do i=1,ntp

        read(5,*)nspec(i)

      enddo
      read(5,*)rcut
      do j=1,ntp
        do i=j,ntp

          read(5,*)radc(i,j)
          radc(j,i)=radc(i,j)

        enddo
      enddo
      read(5,*)nconf

      write(6,'(a,a40)')'filename=',fname
      write(6,'(a,i12)')'natms   =',natms
      write(6,'(a,10i6)')'nspec   =',(nspec(j),j=1,ntp)
      write(6,'(a, i4)')'ntp     =',ntp
      write(6,'(a, i4)')'nconf   =',nconf
      write(6,'(a,1pe12.4)')'rcut    =',rcut

      do i=1,ntp

        do j=i,ntp

          write(6,'(a,i2,a,i2,a,1p,3e12.4)')
     x      'radc(',i,',',j,') = ',radc(i,j)
        enddo

      enddo

c     square of distances

      do i=1,ntp
        do j=1,ntp

          radc(i,j)=radc(i,j)**2

        enddo
      enddo

c     construct type array and zero accumulators

      i=0

      do k=1,mxtp

        num(k)=0

        do j=1,nspec(k)

          i=i+1
          ltype(i)=k

        enddo

        do j=1,2000

          key(j,k)=0
          ihist(j,k)=0

        enddo

      enddo

c     read history file

      keytrj=0

      do iconf=1,nconf

      call hread(fname,text,atmnam,iflg,keytrj,imcon,natms,
     x    tstep,cell,chge,weight,xxx,yyy,zzz,fxx,fyy,fzz)

      if(imcon.gt.0)call invert(cell,rcell,det)


c     convert coordinates to fractional coordinates

        do i=1,natms

          xx=xxx(i)
          yy=yyy(i)
          zz=zzz(i)
          xxx(i)=rcell(1)*xx+rcell(4)*yy+rcell(7)*zz
          yyy(i)=rcell(2)*xx+rcell(5)*yy+rcell(8)*zz
          zzz(i)=rcell(3)*xx+rcell(6)*yy+rcell(9)*zz

        enddo

        call parlst
     x    (natms,rcut,lentry,list,cell,xxx,yyy,zzz)

        call bndstat
     x    (ntp,natms,rcut,ltype,lentry,list,cell,kode,ihist,key,
     x    num,radc,xxx,yyy,zzz)

        goto 200
  100   continue
        backspace(9)
  200   continue
      enddo

      do k=1,ntp

        write(6,'(a,i2)')'histograms for species ',k

        do j=1,ntp

          do i=1,mxshl

            nshell(i,j)=0

          enddo

        enddo

        lshl=0
        do i=1,num(k)

          n=key(i,k)
          do j=1,ntp

            ind(j)=n-100*(n/100)
            lshl=max(lshl,ind(j))
            if(ind(j).gt.(mxshl-1))then
              write(6,'(a)')'error in nshell construction'
              stop
            endif
            nshell(ind(j)+1,j)=nshell(ind(j)+1,j)+ihist(i,k)
            n=n/100

          enddo

c$$$      write(6,'(10i8)')(ind(j),j=1,ntp),ihist(i,k)

        enddo

        write(6,'(a)')'coordination number frequencies'
        write(6,
     x    "('number','   Spec ',i2,'   Spec ',i2,'   Spec ',i2,
     x    '   Spec ',i2,'   Spec ',i2,'   Spec ',i2)")(j,j=1,ntp)
        do i=1,lshl

          write(6,'(i6,10i10)')i-1,(nshell(i,j),j=1,ntp)

        enddo

      enddo

c     bridging oxygen analysis

      nal(1)=0
      nal(2)=0
      do i=0,4
         nbr(i)=0
      enddo

      do i=1,num(ntp)

        n=key(i,ntp)
        do j=1,ntp
          ind(j)=n-100*(n/100)
          n=n/100
        enddo
        if(ind(ntp-1).eq.0)then
          nbr(0)=nbr(0)+ihist(i,ntp)
        elseif(ind(ntp-1).eq.1)then
          if(ntp.eq.3)then
            nal(1)=nal(1)+ind(1)*ihist(i,ntp)
            if(ind(1).eq.0)then
              nbr(1)=nbr(1)+ihist(i,ntp)
            else
              nbr(2)=nbr(2)+ihist(i,ntp)
            endif
          else
            nal(1)=nal(1)+(ind(1)+ind(2))*ihist(i,ntp)
            if(ind(1).eq.0.and.ind(2).eq.0)then
              nbr(1)=nbr(1)+ihist(i,ntp)
            else
              nbr(2)=nbr(2)+ihist(i,ntp)
            endif
          endif
        elseif(ind(ntp-1).eq.2)then
          if(ntp.eq.3)then
            nal(2)=nal(2)+ind(1)*ihist(i,ntp)
            if(ind(1).eq.0)then
              nbr(3)=nbr(3)+ihist(i,ntp)
            else
              nbr(4)=nbr(4)+ihist(i,ntp)
            endif
          else
            nal(2)=nal(2)+(ind(1)+ind(2))*ihist(i,ntp)
            if(ind(1).eq.0.and.ind(2).eq.0)then
              nbr(3)=nbr(3)+ihist(i,ntp)
            else
              nbr(4)=nbr(4)+ihist(i,ntp)
            endif
          endif
        endif

      enddo

      write(6,'(a)')'Analysis of bridging/non-bridging oxygens'

      write(6,'(a,i8)')'uncoordinated oxygens        :',nbr(0)
      write(6,'(a,i8)')'non bridging oxygens - alkali:',nbr(1)
      write(6,'(a,2i8)')'non bridging oxygens + alkali:',nbr(2),nal(1)
      write(6,'(a,i8)')'bridging oxygens     - alkali:',nbr(3)
      write(6,'(a,2i8)')'bridging oxygens     + alkali:',nbr(4),nal(2)

c     determine Q ratios

      do i=0,4
         nqr(i)=0
      enddo

      do i=1,num(ntp-1)

         n=key(i,ntp-1)
         do j=1,ntp
            ind(j)=n-100*(n/100)
            n=n/100
         enddo
         if((ind(ntp).eq.4).and.(ind(ntp-1).le.4))then
            nqr(ind(ntp-1))=nqr(ind(ntp-1))+ihist(i,ntp-1)
         endif

      enddo

      write(6,'(a)')'Analysis of Q Ratios'

      do i=0,4
         write(6,'(a,i2,i8)')'Q_ratio',i,nqr(i)
      enddo

      stop

  300 continue
      write(6,*)'error - problem with history file'

      stop
      end

      subroutine parlst
     x     (natms,rcut,lentry,list,cell,xxx,yyy,zzz)

c***********************************************************************
c
c     subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs decomposition
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c
c***********************************************************************

      implicit real*8 (a-h,o-z)

      parameter (mxatms=1080,mxlist=60)
      logical lchk

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension lentry(mxatms),list(mxatms,mxlist)
      dimension cell(9)

c     set control variables

      mlist=0
      last=natms
      lchk=.true.
      mpm2=natms/2
      npm2=(natms-1)/2

c
c     set cutoff radius

      rclim=rcut**2

c
c     construct pair force neighbour list

      do i=1,mxatms

         lentry(i)=0

      enddo

c
c     outer loop over atoms

      do m=1,mpm2

         if(m.gt.npm2)last=mpm2

c
c     inner loop over atoms

         do i=1,last

c
c     calculate atom indices

            j=i+m
            if(j.gt.natms)j=j-natms

c
c     calculate interatomic displacements

            xdf=xxx(i)-xxx(j)-nint(xxx(i)-xxx(j))
            ydf=yyy(i)-yyy(j)-nint(yyy(i)-yyy(j))
            zdf=zzz(i)-zzz(j)-nint(zzz(i)-zzz(j))

            xdd=cell(1)*xdf+cell(4)*ydf+cell(7)*zdf
            ydd=cell(2)*xdf+cell(5)*ydf+cell(8)*zdf
            zdd=cell(3)*xdf+cell(6)*ydf+cell(9)*zdf

c
c     calculate interatomic distance

            rsq=xdd*xdd+ydd*ydd+zdd*zdd

c
c     running check of neighbour list array capacity

            if(rsq.lt.rclim)then

               lentry(i)=lentry(i)+1
               lentry(j)=lentry(j)+1

               if(lentry(i).gt.mxlist)then

                  mlist=max(mlist,lentry(i))
                  lchk=.false.

               endif

               if(lentry(j).gt.mxlist)then

                  mlist=max(mlist,lentry(j))
                  lchk=.false.

               endif

c
c     compile neighbour list array

               if(lchk)then

                  list(i,lentry(i))=j
                  list(j,lentry(j))=-i

               endif

            endif

         enddo

      enddo
c
c     terminate job if neighbour list array exceeded

      if(.not.lchk)then
         write(6,"('error - neighbour list exceeded',2i10)")
     x        mlist,mxlist
          stop
      endif


      return
      end

      subroutine bndstat
     x     (ntp,natms,rcut,ltype,lentry,list,cell,kode,ihist,key,num,
     x     radc,xxx,yyy,zzz)
c
c***********************************************************************
c
c     calculate bond  statistics
c
c     copyright daresbury laboratory
c     author w.smith
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter (mxatms=1080,mxlist=60,mxrdf=128,mxtp=4)

      logical lfind

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension lentry(mxatms),list(mxatms,mxlist),num(mxtp)
      dimension ltype(mxatms),kode(mxatms,mxtp)
      dimension ihist(2000,mxtp),key(2000,mxtp)
      dimension cell(9),radc(mxtp,mxtp)

      data pi/3.1415926536d0/

      rcsq=rcut**2

c     initialise kode array

      do j=1,mxatms

        do i=1,mxtp

          kode(j,i)=0

        enddo

      enddo
c
c     determine bonds

      do i=1,natms

         ityp=ltype(i)

         rxi=xxx(i)
         ryi=yyy(i)
         rzi=zzz(i)

         do jj=1,lentry(i)

            j=iabs(list(i,jj))
            jtyp=ltype(j)

            rjix=xxx(j)-rxi-nint(xxx(j)-rxi)
            rjiy=yyy(j)-ryi-nint(yyy(j)-ryi)
            rjiz=zzz(j)-rzi-nint(zzz(j)-rzi)

            xji=cell(1)*rjix+cell(4)*rjiy+cell(7)*rjiz
            yji=cell(2)*rjix+cell(5)*rjiy+cell(8)*rjiz
            zji=cell(3)*rjix+cell(6)*rjiy+cell(9)*rjiz

            rji2=xji*xji+yji*yji+zji*zji

            if(rcsq.gt.rji2)then

               if(radc(ityp,jtyp).gt.rji2)then
                  kode(i,jtyp)=kode(i,jtyp)+1
                  if(kode(i,jtyp).gt.100)then
                     write(6,'(a)')'error - kode > 100'
                     stop
                  endif
               endif

            endif

         enddo

      enddo

      do i=1,natms

         ik=0
         it=ltype(i)
         do j=1,ntp

            ik=ik+kode(i,j)*100**(j-1)

         enddo

         lfind=.false.
         do j=1,num(it)
            if(ik.eq.key(j,it))then
               lfind=.true.
               ihist(j,it)=ihist(j,it)+1
            endif
         enddo
         if(.not.lfind)then
            num(it)=num(it)+1
            if(num(it).gt.2000)then
               write(6,'(a)')'error - num >2000'
               stop
            endif
            key(num(it),it)=ik
            ihist(num(it),it)=1
         endif

      enddo

      return
      end
      subroutine invert(a,b,d)
c
c***********************************************************************
c
c     routine to invert a 3 * 3 matrix using cofactors
c
c     copyright daresbury laboratory
c     author w.smith
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
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

