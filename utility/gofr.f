      subroutine gofr
     x  (dl_fmt,atnam1,atnam2,fname,title,name,mxrad,mxatms,natms,
     x  nconf,ncorr,isampl,numfil,rcut,avcell,rdf,xyz,vel,frc,cell,
     x  bcell,weight,chge,status)

c*********************************************************************
c
c     dl_poly program to calculate radial distribution function
c     for selected atoms
c
c     copyright daresbury laboratory 1999
c     author  w.smith march 1999
c
c*********************************************************************

      implicit none

      real*8 pi
      parameter (pi=3.1415926536d0)

      logical dl_fmt,all,new
      character*2 fnum
      character*40 fname
      character*80 title
      integer i,iconf,iflg,ii,imcon,isampl,j,jj,k,keytrj,last,matm1
      integer matm2,matms,mxatms,mxrad,natms,nconf,nstep,nta,ntb
      integer numfil,ncorr,status,kk
      character*8 atnam1,atnam2,name(mxatms)
      real*8 delr,det,f1,f2,rcut,rcut2,rmsx,rmsy,rmsz,rnorm
      real*8 rrr,rsq,timstp,uuu,volm,vvv,www,xxx,yyy,zzz,tcut
      real*8 rdf(mxrad),avcell(9),cell(9),bcell(10),weight(mxatms)
      real*8 xyz(3,mxatms),vel(3,mxatms),frc(3,mxatms),rcell(9)
      real*8 chge(mxatms)

      status=0

      call strcut(atnam1,8)
      call strcut(atnam2,8)
      call strcut(fname,40)

      all=(atnam1.eq."ALL".or.atnam2.eq."ALL".or.
     x     atnam1.eq."all".or.atnam2.eq."all")

      if(all)then

         atnam1="ALL"
         atnam2="ALL"

      endif

c     open report file

      open(9, file="RDFPRG.rep")

      write(9,'(a)')
     x  'Radial Distribution Function Program'

c     check rdf array dimension

      if(ncorr.gt.mxrad)then

        write(9,'(a,i5)')'Warning - RDF array dimension reset to',ncorr

      endif

      rcut2=rcut**2
      delr=rcut/dble(ncorr)

      write(9,'(a,a40)')'Label  of target HISTORY file : ',fname
      write(9,'(a,a8)')'Label  of first atom          : ',atnam1
      write(9,'(a,a8)')'Label  of second atom         : ',atnam2
      write(9,'(a,i6)')'Length of RDF arrays          : ',ncorr
      write(9,'(a,i6)')'Number of requested configs.  : ',nconf
      write(9,'(a,i6)')'Sampling interval             : ',isampl
      write(9,'(a,f6.3)')'Selected cutoff radius        : ',rcut
      write(9,'(a,f6.3)')'Bin width in RDF              : ',delr

c     initialise rdf array

      do i=1,mxrad

        rdf(i)=0.d0

      enddo

c     initialize average cell vectors

      do i=1,9

        avcell(i)=0.d0

      enddo

c     set default cell vectors

      do i=1,9
        cell(i)=0.d0
      enddo
      cell(1)=1.d0
      cell(5)=1.d0
      cell(9)=1.d0

      iflg=0
      keytrj=0
      new=.true.

      do iconf=0,nconf-1

        if(dl_fmt)then

          call hread
     x      (new,fname,title,name,iflg,imcon,keytrj,natms,nstep,timstp,
     x      cell,chge,weight,xyz,vel,frc,status)

        else

          call uread
     x      (new,fname,title,name,iflg,imcon,keytrj,natms,nstep,timstp,
     x      cell,chge,weight,xyz,vel,frc,status)


        endif

c     check specified cutoff

        if(iconf.eq.0)then

          call dcell(cell,bcell)

          tcut=0.5d0*min(bcell(7),bcell(8),bcell(9))
          if(rcut.gt.tcut)then

            rcut=tcut
            rcut2=rcut**2
            delr=rcut/dble(ncorr)
            write(9,'(a,f10.6)')'Warning - RDF  cutoff  reset to',tcut
            write(9,'(a,f10.6)')'Warning - RDF binwidth reset to',delr

          endif

        endif

        if(iflg.lt.0)go to 100

c     check boundary condition

        if(imcon.ge.1.and.imcon.le.3)then

          call invert(cell,rcell,det)

        else

          write(9,'(a)')
     x      'Error - incorrect boundary condition'
          status=-2
          return

        endif

        j=0
        do i=1,natms

          if(all.or.name(i).eq.atnam1.or.name(i).eq.atnam2)then

            j=j+1
            name(j)=name(i)
            xxx=xyz(1,i)
            yyy=xyz(2,i)
            zzz=xyz(3,i)
            xyz(1,j)=xxx*rcell(1)+yyy*rcell(4)+zzz*rcell(7)
            xyz(2,j)=xxx*rcell(2)+yyy*rcell(5)+zzz*rcell(8)
            xyz(3,j)=xxx*rcell(3)+yyy*rcell(6)+zzz*rcell(9)
          endif

        enddo

        matms=j

c     running average of cell vectors

        f1=(dble(iconf)/dble(iconf+1))
        f2=1.d0/dble(iconf+1)
        do i=1,9

          avcell(i)= f1*avcell(i)+f2*cell(i)

        enddo

        last=matms
        matm1=matms/2
        matm2=(matms-1)/2

        do kk=1,matm1

          if(kk.gt.matm2)last=matm1

          do ii=1,last

            jj=ii+kk
            if(jj.gt.matms)jj=jj-matms

            if(all.or.
     x        (atnam1.eq.name(ii).and.atnam2.eq.name(jj)).or.
     x        (atnam1.eq.name(jj).and.atnam2.eq.name(ii)))then

              uuu=xyz(1,ii)-xyz(1,jj)
              vvv=xyz(2,ii)-xyz(2,jj)
              www=xyz(3,ii)-xyz(3,jj)

              uuu=uuu-nint(uuu)
              vvv=vvv-nint(vvv)
              www=www-nint(www)

              rmsx=uuu*avcell(1)+vvv*avcell(4)+www*avcell(7)
              rmsy=uuu*avcell(2)+vvv*avcell(5)+www*avcell(8)
              rmsz=uuu*avcell(3)+vvv*avcell(6)+www*avcell(9)

              rsq=rmsx**2+rmsy**2+rmsz**2

              if(rsq.lt.rcut2)then

                k=int(sqrt(rsq)/delr)+1
                rdf(k)=rdf(k)+1.d0

              endif

            endif

          enddo

        enddo

      enddo

  100 continue

      if(iflg.lt.0)iconf=iconf-1
      write(9,'(a,i6)')'Number of configurations sampled: ',iconf

c     normalise and print radial distribution function

      volm=abs(avcell(1)*(avcell(5)*avcell(9)-avcell(6)*avcell(8))
     x  +avcell(4)*(avcell(3)*avcell(8)-avcell(2)*avcell(9))
     x  +avcell(7)*(avcell(2)*avcell(6)-avcell(3)*avcell(5)))

      if(all)then

         nta=matms
         ntb=matms

      else

         nta=0
         ntb=0
         do i=1,matms

            if(name(i).eq.atnam1)nta=nta+1
            if(name(i).eq.atnam2)ntb=ntb+1

         enddo

      endif

c     open the correlation data files

      numfil=numfil+1
      if(numfil.gt.99)then
        fnum="XX"
        close (88)
        open(88,file='RDFDAT.'//fnum)
      else
        write(fnum,'(i2.2)')numfil
        close(88)
        open(88,file='RDFDAT.'//fnum)
      endif

      write(88,'(a80)')title
      write(88,'(2i10)')1,ncorr

      if(atnam1.eq.atnam2) then

        rnorm=volm/(2.d0*pi*delr**3*dble(nta)*dble(nta)*dble(iconf))
        write(88,'(2a8,1pe14.6)')atnam1,atnam1,rcut

      else

        rnorm=volm/(4.d0*pi*delr**3*dble(nta)*dble(ntb)*dble(iconf))
        write(88,'(2a8,1pe14.6)')atnam1,atnam2,rcut

      endif

      do i=1,ncorr

        rdf(i)=rnorm*rdf(i)/((dble(i)-0.5d0)**2+1.d0/12.d0)

      enddo

      do i=1,ncorr

        rrr=delr*(dble(i)-0.5d0)
        write(88,'(1p,2e14.6)')rrr,rdf(i)

      enddo

      close(88)
      write(9,'(a)')'RDF file RDFDAT.'//fnum//' created'
      close(9)
      return
      end
      subroutine strcut(a,n)

      character*1 a(n)
      logical y

      y=.false.
      do i=1,n

        if(ichar(a(i)).eq.0)y=.true.
        if(y)a(i)=" "

      enddo
      return
      end
      subroutine hread
     x  (new,history,cfgname,name,iflg,imcon,keytrj,natms,
     x   nstep,tstep,cell,chge,weight,xyz,vel,frc,status)

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

      implicit none

      logical new

      character*80 cfgname
      character*40 history
      character*8 name(*),step

      integer iflg,imcon,keytrj,natms,nstep,status,nhist,ktrj,matms,i,j

      real*8 tstep,vx,vy,vz,fx,fy,fz

      real*8 cell(9)
      real*8 chge(*),weight(*)
      real*8 xyz(3,*),vel(3,*),frc(3,*)

      save ktrj
      data nhist/77/

      status=0

c     open history file if new job

      if(new)then

        open(nhist,file=history,status='old',err=100)

        read(nhist,'(a80)',err=200) cfgname
        write(*,'(a,a)')'History file header: ',cfgname
        read(nhist,'(2i10)',end=200) ktrj,imcon
        if(keytrj.gt.ktrj)then

          if(ktrj.eq.0)then

            write(*,'(a)')'Error - no velocities in file'
            status=-1
            close (nhist)
            new=.true.
            return

          endif
          if(keytrj.gt.1)then

            write(*,'(a)')'Error - no forces in file'
            status=-2
            close (nhist)
            new=.true.
            return

          endif
        endif

        new=.false.

      endif

      read(nhist,'(a8,4i10,f12.6)',end=200)
     x     step,nstep,matms,ktrj,imcon,tstep

      if(matms.gt.natms)then

        write(*,'(a)')'Error - too many atoms in MD cell'
        write(*,'(a,i6,a)')'File contains',matms,' atoms'
        status=-3
        close (nhist)
        new=.true.
        return

      endif

      natms=matms

      if(imcon.gt.0) then

         read(nhist,*,end=200) cell(1),cell(2),cell(3)
         read(nhist,*,end=200) cell(4),cell(5),cell(6)
         read(nhist,*,end=200) cell(7),cell(8),cell(9)

      endif

      do i = 1,natms

        read(nhist,'(a8,i10,2f12.6)',end=200)
     x    name(i),j,weight(i),chge(i)
        read(nhist,*,end=200) xyz(1,i),xyz(2,i),xyz(3,i)
        if(keytrj.ge.1)then
          read(nhist,*,end=200) vel(1,i),vel(2,i),vel(3,i)
        else if(ktrj.ge.1)then
          read(nhist,*,end=200) vx,vy,vz
        endif
        if(keytrj.ge.2)then
          read(nhist,*,end=200) frc(1,i),frc(2,i),frc(3,i)
        else if(ktrj.ge.2)then
          read(nhist,*,end=200) fx,fy,fz
        endif

      enddo

      if(iflg.lt.0)then

        close (nhist)
        new=.true.

      endif

      return

  100 continue

      write(*,'(a)')'Error - History file not found'
      status=-4
      return

  200 continue
      write(*,'(a)')'Warning - end of History file encountered'
      close (nhist)
      iflg=-1
      new=.true.

      return
      end
