      implicit real*8(a-h,o-z)
      parameter(pi=3.1415926536d0)
      parameter(mxtheta=512,mxatms=1080,mxconf=3000,mxcell=2000)
      
      character*40   fname
      character*80   title
      character*8    name(mxatms)
      dimension lst(mxcell),lct(mxcell),ltype(mxatms)   
      dimension nix(27),niy(27),link(mxatms),niz(27)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension cell(9),rcell(9),cprp(10),theta(mxtheta) 
      dimension latinx(mxatms),chge(mxatms),weight(mxatms)
      data nix/0,-1,-1,-1,0,0,-1,-1,-1,0,1,-1,0,1,
     x  1,1,1,0,0,1,-1,1,0,-1,1,0,-1/
      data niy/0,0,-1,1,1,0,0,0,-1,-1,-1,1,1,1,
     x  0,1,-1,-1,0,0,0,1,1,1,-1,-1,-1/
      data niz/0,0,0,0,0,1,1,1,1,1,1,1,1,1,
     x  0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1/ 
c     name of the history
      fname="HISTORY"
c     total number of atoms in the history file configuration
      natms=1080
c     max number of configuration to sample
      nconf=1000
c     the length of theta array
      ncorr=100
      rcut=2.d0
      deltheta=pi/dble(ncorr) 
      iflg=0
      keytrj=0
      do i=1,mxtheta
        theta(i)=0.d0
      enddo
      do iconf=0,nconf-1
        call hread(fname,title,name,iflg,imcon,keytrj,
     x    natms,nstep,timestp,cell,chge,weight,xxx,yyy,
     x    zzz,vxx,vyy,vzz,fxx,fyy,fzz)
        if(iflg.lt.0)go to 100

c     check for appropriate boundary conditions

        if(imcon.ge.1.and.imcon.le.3) then
          
          call invert(cell,rcell,det)
          call dcell(cell,cprp)

        endif

c     calculate link cell numbers

        nbx=int(cprp(7)/(rcut+1.d-6))
        nby=int(cprp(8)/(rcut+1.d-6))
        nbz=int(cprp(9)/(rcut+1.d-6))
        ncells=nbx*nby*nbz

c     transform atomic coodinates and construct link cells
        do i=1,ncells
          lct(i)=0
          lst(i)=0
        enddo
        xdc=dble(nbx)
        ydc=dble(nby)
        zdc=dble(nbz)
        do i=1,natms

          sxx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          syy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          szz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          xxx(i)=sxx
          yyy(i)=syy
          zzz(i)=szz
          ix=int(xdc*(sxx+0.5d0))
          iy=int(ydc*(syy+0.5d0))
          iz=int(zdc*(szz+0.5d0))
          k=1+ix+nbx*(iy+nby*iz)
          lst(k)=lst(k)+1
          link(i)=lct(k)
          lct(k)=i

        enddo

c     loop over central atoms of angles

        ix=0
        iy=1
        iz=1
        do icell=1,ncells
          ix=ix+1
          if(ix.gt.nbx)then
            ix=1
            iy=iy+1
            if(iy.gt.nby)then
              iy=1
              iz=iz+1
            endif
          endif

c     construct mini-list of neighbour cell contents

        k=0
        do kk=1,27

          jx=ix+nix(kk)
          jy=iy+niy(kk)
          jz=iz+niz(kk)
          if(jx.gt.nbx)jx=1
          if(jy.gt.nby)jy=1
          if(jz.gt.nbz)jz=1
          if(jx.lt.1)jx=jx+nbx
          if(jy.lt.1)jy=jy+nby
          if(jz.lt.1)jz=jz+nbz
          jcell=jx+nbx*(jy-1+nby*(jz-1))
          j=lct(jcell)

          do ii=1,lst(jcell)
            k=k+1
            latinx(k)=j
            j=link(j)
          enddo

        enddo
        limit=k

        do ii=1,lst(icell)

          i=latinx(ii)
          if(name(i).eq."Si4+") then 
            last=limit
            do kk=1,limit/2

              if(kk.gt.(limit-1)/2)last=limit/2
              do jj=1,last
                j=latinx(jj)
                jk=jj+kk
                if(jk.gt.limit)jk=jk-limit
                k=latinx(jk)

c     make labels etc consistent with angfrc.f

                ia = j
                ib = i
                ic = k

                if(name(ia).eq."O2-")then          
                  sxab = xxx(ia)-xxx(ib)
                  sxab = sxab-nint(sxab)
                  syab = yyy(ia)-yyy(ib)
                  syab = syab-nint(syab)
                  szab = zzz(ia)-zzz(ib)
                  szab = szab-nint(szab)
                  xab=cell(1)*sxab+cell(4)*syab+cell(7)*szab
                  if(abs(xab).lt.rcut)then
                    yab=cell(2)*sxab+cell(5)*syab+cell(8)*szab
                    if(abs(yab).lt.rcut)then
                      zab=cell(3)*sxab+cell(6)*syab+cell(9)*szab
                      if(name(ic).eq."O2-")then 
                        if(abs(zab).lt.rcut)then
                          sxbc = xxx(ic)-xxx(ib)
                          sxbc = sxbc-nint(sxbc)
                          sybc = yyy(ic)-yyy(ib)
                          sybc = sybc-nint(sybc)
                          szbc = zzz(ic)-zzz(ib)
                          szbc = szbc-nint(szbc)
                          xbc=cell(1)*sxbc+cell(4)*sybc+cell(7)*szbc
                          if(abs(xbc).lt.rcut)then
                            ybc=cell(2)*sxbc+cell(5)*sybc+cell(8)*szbc
                            if(abs(ybc).lt.rcut)then
                              zbc=cell(3)*sxbc+cell(6)*sybc+cell(9)*szbc
                              if(abs(zbc).lt.rcut)then
                                rab=sqrt(xab*xab+yab*yab+zab*zab)
                                rbc=sqrt(xbc*xbc+ybc*ybc+zbc*zbc)
                                if(rcut.ge.max(rab,rbc))then

                                  rrab = 1.d0/rab
                                  rrbc = 1.d0/rbc

c     normalise direction vectors

                                  xab = xab*rrab
                                  yab = yab*rrab
                                  zab = zab*rrab

                                  xbc = xbc*rrbc
                                  ybc = ybc*rrbc
                                  zbc = zbc*rrbc

                                  thet=acos(xab*xbc+yab*ybc+zab*zbc)

                                  kkk=int(thet/deltheta)
                                  theta(kkk)=theta(kkk)+1.d0 

                                endif
                              endif
                            endif
                          endif
                        endif
                      endif
                    endif
                  endif
                endif
              enddo
            enddo
          endif
        enddo
      enddo
      enddo
  100 continue
      nconf=iconf-1

c     normalise the angle

      rnorm=0
      do i=2,ncorr-1,2
        rnorm=rnorm+(theta(i-1)+4.d0*theta(i)+theta(i+1))/3.d0 
      enddo
      do i=1,ncorr
        theta(i)=theta(i)/rnorm 
      enddo
      open(8,file='output')
      do i=1,ncorr
        write(8,'(1p,2e14.6)')deltheta*(dble(i)-0.5d0),theta(i)
      enddo
      close(8) 
      end
      subroutine dcell(aaa,bbb)

c     
c***********************************************************************
c     
c     dl_poly subroutine to calculate the dimensional properties of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c     
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c     
c     copyright daresbury laboratory 1992
c     author - w. smith         july 1992
c     
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c     
c***********************************************************************
c     

      implicit real*8 (a-h,o-z)

      dimension aaa(9),bbb(10)
c     
c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))
c     
c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))
c     
c     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)
c     
c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)
c     
c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

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
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
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
c     itt
c     2010-10-30 17:20:49
c     1.3
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

