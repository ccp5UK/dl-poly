      program solvadd
     
c***********************************************************************
c     
c     DL_POLY utility
c     
c     Program to add single site solvent molecules to a structure to
c     fill out the MD cell.
c     Assumes atomic postions are in a form compatible
c     with the CONFIG file used in DL_POLY with periodic boundary
c     conditions.
c     
c     solvent is added from a f.c.c. lattice
c     
c     input file = CONFIG
c     output file = CONFIG.solv
c     
c     copyright Daresbury laboratory 1993
c     author -     t.  forester  feb 1993
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c
c***********************************************************************
      
      implicit real*8(a-h,o-z)
      parameter(mxatm = 100000)
      parameter(mxcell= 1000000)
      
      character*80 rdline
      character*8 name(mxatm)
      dimension cell(9),rcell(9)
      dimension xxx(mxatm),yyy(mxatm),zzz(mxatm)
      dimension ssx(mxatm),ssy(mxatm),ssz(mxatm)
      dimension ox(4),oy(4),oz(4)
      dimension nix(27),niy(27),niz(27)
      dimension lstrem(mxatm)
      dimension xd(mxatm),yd(mxatm),zd(mxatm)
      
      dimension lct(mxcell),link(mxatm)
      
      data nix/  0, 1, 0, 0,-1, 1, 0,-1, 1, 0,-1, 1,-1, 1,
     x   -1, 0, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1/
      data niy/  0, 0, 1, 0, 1, 1,-1, 0, 0, 1,-1,-1, 1, 1,
     x   0,-1, 0,-1,-1, 1, 0, 0,-1, 1, 1,-1,-1/
      data niz/  0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
     x   0, 0,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1/
      
c     
c     define solvent
      
      write(*,*)' input lattice spacing for solvent (f.c.c.) '
      read (*,*) axcub
      
      write(*,*)' solvent label (a8)'
      read(*,'(a8)') ow
c     
c     tolerance for solvent - solute distance
      
      write(*,*) 'minimum solvent-solute distance ? '
      read(*,*) tolnce
      tolnce= tolnce**2
      
c     
c     tolerance for solvent/solvent distance across periodic boundary
      
      write(*,*)'minimum permitted solvent/solvent distance ? ',
     x   '  (across periodic boundary)'
      read(*,*) tolow
      tolow=tolow**2
c     
c     device identities
      
      ird = 8
      iwr = 9
      
      open(ird,file='CONFIG')
      open(iwr,file='CONFIG.solv')
      
      read(ird,'(a80)') rdline
      write(iwr,'(a80)') rdline
      
      read(ird,'(2i10)') levcfg,imcon
      if(imcon.eq.0) call error(50)

      read(ird,'(3f20.8)') cell

      if(imcon.ne.10) then
        write(iwr,'(2i10)') 0,imcon
        write(iwr,'(3f20.8)') cell
      else
        write(iwr,'(2i10)') 0,0
      endif
      
      do i = 1,mxatm
        
        read(ird,'(a8)',end=100) name(i)
        read(ird,'(3f20.8)',end=100) xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0) read(ird,*,end=100)
        if(levcfg.gt.1) read(ird,*,end=100)
        
      enddo
      
  100 natm = i-1
      
      call invert(cell,rcell,det)
      
      imcon1=imcon
      imcon = min(imcon,3)
c     
c     place solute atoms into link list
      
      call cells
     x   (natm,ncells,imcon,tolnce,xdc,ydc,zdc,
     x   lct,link,cell,rcell,ssx,ssy,ssz,xxx,yyy,zzz)
      
c     
c     create f.c.c. lattice
      
      axh = axcub*0.5d0
      
      ox(1) = axh
      oy(1) = axh
      oz(1) = 0.d0
      
      ox(2) = axh
      oy(2) = 0.d0
      oz(2) = axh
      
      ox(3) = 0.d0
      oy(3) = axh
      oz(3) = axh
      
      ox(4) = 0.d0
      oy(4) = 0.d0
      oz(4) = 0.d0
      
c     
c     define number of solvent blocks to use
      
      imax=int((abs(cell(1))+abs(cell(4))+abs(cell(7)))/axcub+.99d0)/2+1 
      jmax=int((abs(cell(2))+abs(cell(5))+abs(cell(8)))/axcub+.99d0)/2+1
      kmax=int((abs(cell(3))+abs(cell(6))+abs(cell(9)))/axcub+.99d0)/2+1 
      
c     
c     inverse of MD cell
      
      call invert(cell,rcell,det)
      
      write(*,*) 'inserting solvent '
      write(*,*)
      
      iplus = 0
      nsbcll=27
      call solvins
     x   (name,ow,iplus,imax,jmax,kmax,4,natm,nsbcll,mxatm,axcub,
     x   tolnce,xdc,ydc,zdc,lct,link,nix,niy,niz,xxx,yyy,zzz,
     x   ssx,ssy,ssz,cell,rcell,ox,oy,oz)
      
c     
c     check for bad interactions (less than tolow apart)
      
      natm1 = natm+(iplus)
      tolnce = max(tolow,9.d0)
      
c     
c     place all atoms in link cells
      
      call cells
     x   (natm1,ncells,imcon,tolnce,xdc,ydc,zdc,
     x   lct,link,cell,rcell,ssx,ssy,ssz,xxx,yyy,zzz)
      
      tolnce = tolow
c     
c     additional checks on non standard periodic boundaries
      
      if(imcon1.gt.3.and.imcon1.ne.6) then
        
        imcon = imcon1
        call imgchk
     x     (imcon,natm,iplus,tolnce,lstrem,
     x     ssx,ssy,ssz,xxx,yyy,zzz,cell,xd,yd,zd)
        
      else
        
        call solvchk
     x     (iplus,ncells,natm,14,tolnce,xdc,ydc,zdc,lct,
     x     link,lstrem,nix,niy,niz,ssx,ssy,ssz,xxx,yyy,zzz,cell)
        
        
      endif
      
c     
c     write out final structure
      
      do i = 1,natm+iplus
        
        write(iwr,'(a8,i10)') name(i),i
        write(iwr,'(3f20.8)') xxx(i),yyy(i),zzz(i)
      enddo
      
      write(*,*)
      write(*,*) 'all done: added ',iplus,' solvents '
      write(*,*) 'output in CONFIG.solv'
      
      end
      
      subroutine images
     x   (imcon,idnode,mxnode,natm,cell,xxx,yyy,zzz)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating the minimum image
c     of atom pairs within a specified MD cell
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 rectangular (slab) boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram, z non-periodic
c     
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge 
c     coordinate arrays
c     
c***********************************************************************
c     
      
      implicit real*8 (a-h,o-z)
      
      dimension xxx(*),yyy(*),zzz(*)
      dimension cell(9),rcell(9)
      
      if(imcon.eq.1)then
c     
c     standard cubic boundary conditions
        
        aaa=1.d0/cell(1)
        
        do i=idnode+1,natm,mxnode
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
          
        enddo
        
      else if(imcon.eq.2)then
c     
c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=idnode+1,natm,mxnode
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))
          
        enddo
        
      else if(imcon.eq.3)then
c     
c     parallelepiped boundary conditions
        
        call invert(cell,rcell,det)
        if(abs(det).lt.1.d-6)call error(120)
        
        do i=idnode+1,natm,mxnode
          
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
        
      else if(imcon.eq.4)then
c     
c     truncated octahedral boundary conditions
        
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-3.and.
     x     abs(cell(5)-cell(9)).lt.1.d-3)) call error(120)
        
        aaa=1.d0/cell(1)
        
        do i=idnode+1,natm,mxnode
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.
     x       (0.75d0*cell(1)))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(1),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then
c     
c     rhombic dodecahedral boundary conditions
        
        rt2=sqrt(2.d0)
        if(abs(cell(1)-cell(5)).gt.1.d-3.or.
     x     abs(cell(9)-cell(1)*rt2).gt.1.d-3)
     x     call error(120)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=idnode+1,natm,mxnode
          
          ssx= xxx(i)*aaa
          ssy= yyy(i)*aaa
          ssz= zzz(i)*bbb
          
          ssx = ssx-nint(ssx)
          ssy = ssy-nint(ssy)
          ssz = ssz-nint(ssz)
          
          if(abs(ssx)+abs(ssy)+2.d0*abs(ssz).gt.1.d0) then
            
            ssx = ssx - sign(0.5d0,ssx)
            ssy = ssy - sign(0.5d0,ssy)
            ssz = ssz - sign(0.5d0,ssz)
            
          endif
          
          xxx(i) = ssx*cell(1)
          yyy(i) = ssy*cell(1)
          zzz(i) = ssz*cell(9)
          
        enddo
        
      else if(imcon.eq.6)then
c     
c     parallelogram in x-y, z non periodic
        
        call invert(cell,rcell,det)
        if(abs(det).lt.1.d-6)call error(120)
        
        do i=idnode+1,natm,mxnode
          
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          
        enddo
        
      endif
      
      return
      end
      
      subroutine invert(a,b,d)
c     
c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith april 1992.
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
      
      subroutine error(kode)
      implicit real*8(a-h,o-z)
      
      
      write(*,*)
      if(kode.eq.50) then
        write(*,*) ' error - invalid image convention'
      elseif(kode.eq.60) then
        write(*,*) ' error : system volume too small'
      elseif(kode.eq.100) then
        write(*,*) ' error : system too large for link cells'
      elseif(kode.eq.110) then
        write(*,*) ' error : system too small for link cells'
      elseif(kode.eq.120) then
        write(*,*) ' error : invalid cell vectors '
      elseif(kode.eq.200) then
        write(*,*) ' too many atoms in system: increase parameter mxatm'
      else
        write(*,*) ' unnamed error found'
      endif
      
      if(kode.gt.0) stop
      
      return
      end
      
      subroutine cells
     x   (natm,ncells,imcon,tolnce,xdc,ydc,zdc,
     x   lct,link,cell,rcell,ssx,ssy,ssz,xxx,yyy,zzz)
      
      implicit real*8(a-h,o-z)
      parameter(mxcell=1000000)
      
      logical linc
      
      dimension lct(*),link(*)
      dimension ssx(*),ssy(*),ssz(*)
      dimension cell(*),rcell(*),celprp(10)
      dimension xxx(*),yyy(*),zzz(*)
      
      call dcell(cell,celprp)
      if(celprp(10).lt.1.d0) call error(60)

      rlink = sqrt(tolnce)
      ilx = int(celprp(7)/rlink)
      ily = int(celprp(8)/rlink)
      ilz = int(celprp(9)/rlink)
      
      ncells = ilx*ily*ilz
c     
c     check for enough link cells
      
      if(ncells.gt.mxcell) call error(100)
c     
c     check system is big enough
      
      linc =.false.
      if(ilx.lt.3) linc = .true.
      if(ily.lt.3) linc = .true.
      if(ilz.lt.3) linc = .true.
      
      if(linc) call error(110)
      
c     
c     calculate link cell indices
      
      do i = 1,ncells
        
        lct(i)=0
        
      enddo
c     
c     link-cell cutoff for reduced space
      
      xdc = dble(ilx)
      ydc = dble(ily)
      zdc = dble(ilz)
      
c     
c     reduced space coordinates
      
      call images(imcon,0,1,natm,cell,xxx,yyy,zzz)
      
      do i = 1,natm
        
        ssx(i)=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))+0.5d0
        ssy(i)=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))+0.5d0
        ssz(i)=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))+0.5d0
        
      enddo
c     
c     link neighbours 
      
      do i = 1,natm
        
        ix = int(xdc*ssx(i))
        iy = int(ydc*ssy(i))
        iz = int(zdc*ssz(i))
        
        icell = 1+ix+ilx*(iy+ily*iz)
        
        j = lct(icell)
        lct(icell)=i
        link(i)=j
        
      enddo
      
      return 
      end
      
      subroutine lnkimg(jy,ily,cy)
      implicit real*8(a-h,o-z)
      
      if(jy.ge.ily) then
        
        jy = jy-ily
        cy = 1.d0
        
      elseif(jy.lt.0) then
        
        jy = jy+ily
        cy =-1.d0
        
      endif
      
      return
      end
      
      subroutine solvins
     x   (name,ow,iplus,imax,jmax,kmax,naa,natm,nsbcll,mxatm,axcub,
     x   tolnce,xdc,ydc,zdc,lct,link,nix,niy,niz,xxx,yyy,zzz,
     x   ssx,ssy,ssz,cell,rcell,ox,oy,oz)
      
      implicit real*8(a-h,o-z)
      
      logical ladd
      
      character*8 ow,name(*)
      dimension lct(*),link(*)
      dimension ssx(*),ssy(*),ssz(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension ox(*),oy(*),oz(*)
      dimension nix(*),niy(*),niz(*)
      dimension cell(*),rcell(*)
      
      ilx=int(xdc+0.5d0)
      ily=int(ydc+0.5d0)
      ilz=int(zdc+0.5d0)
      
      do i = -imax,imax
        
        xi = dble(i)*axcub
        
        do j = -jmax,jmax
          
          yi = dble(j)*axcub
          
          do k = -kmax,kmax
            
            zi = dble(k)*axcub 
            
            do l = 1,naa
              
              xo = xi + ox(l) 
              yo = yi + oy(l)
              zo = zi + oz(l)
c     
c     check if water is in basic MD cell
              
              sx=(rcell(1)*xo+rcell(4)*yo+rcell(7)*zo)
              sy=(rcell(2)*xo+rcell(5)*yo+rcell(8)*zo)
              sz=(rcell(3)*xo+rcell(6)*yo+rcell(9)*zo)
              
              if(abs(sx).le..5d0) then
                if(abs(sy).le..5d0) then
                  if(abs(sz).le..5d0)then
c     
c     compute link cell index
                    
                    sx=sx+0.5d0
                    sy=sy+0.5d0
                    sz=sz+0.5d0
                    
                    ix = int(xdc*sx)
                    iy = int(ydc*sy)
                    iz = int(zdc*sz)
c     
c     flag to add water
                    
                    ladd = .true.
c     
c     loop over nearby link cells of solute
                    
                    nsbcll = 27
                    do kc = 1,nsbcll
                      
                      cx = 0.d0
                      cy = 0.d0
                      cz = 0.d0
                      
                      jx=ix+nix(kc)
                      jy=iy+niy(kc)
                      jz=iz+niz(kc)
c     
c     minimum image convention for link cells
                      
                      call lnkimg(jx,ilx,cx)
                      call lnkimg(jy,ily,cy)
                      call lnkimg(jz,ilz,cz)
c     
c     index of neighbouring cell
                      
                      jc = 1+jx+ilx*(jy+ily*jz)
                      jj=lct(jc)
c     
c     ignore if empty
                      
                      if(jj.gt.0) then
                        
  300                   continue
                        
c     
c     distance in real space : minimum image applied
                        
                        sxd = ssx(jj)-sx+cx
                        syd = ssy(jj)-sy+cy
                        szd = ssz(jj)-sz+cz
                        
                        xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
                        yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
                        zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                        
                        rsq = xd*xd+yd*yd+zd*zd
c     
c     test of distance
                        if(rsq.lt.tolnce) then
                          ladd =.false.
                        endif
                        
                        jj=link(jj)
                        if(jj.gt.0) goto 300
                        
                      endif
                      
                    enddo
c     
c     add solvent to structure
                    
                    if(ladd) then
                      
                      iplus = iplus + 1
                      k1 = natm+iplus
                      if(k1.gt.mxatm) call error(200)
                      
                      xxx(k1) = xo
                      yyy(k1) = yo
                      zzz(k1) = zo
                      name(k1) = ow
                      
                    endif
                    
                  endif
                endif
              endif
              
            enddo
            
          enddo
          
        enddo
        
      enddo
      
      return 
      end
      
      subroutine solvchk
     x   (iplus,ncells,natm,nsbcll,tolnce,xdc,ydc,zdc,lct,
     x   link,lstrem,nix,niy,niz,ssx,ssy,ssz,xxx,yyy,zzz,cell)
      
      implicit real*8(a-h,o-z)
      
      dimension lct(*),link(*)
      dimension ssx(*),ssy(*),ssz(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension nix(*),niy(*),niz(*)
      dimension cell(*)
      dimension lstrem(*)
      
      ilx = int(xdc+0.5d0)
      ily = int(ydc+0.5d0)
      ilz = int(zdc+0.5d0)
      
      ix=0
      iy=0
      iz=0
      
      irem =0
      
      do ic = 1,ncells
        
        ii=lct(ic)
        if(ii.gt.0) then
c     
c     secondary loop over subcells
          
          do kc = 1,nsbcll
            
            ii = lct(ic)
            
            cx = 0.d0
            cy = 0.d0
            cz = 0.d0
            
            jx=ix+nix(kc)
            jy=iy+niy(kc)
            jz=iz+niz(kc)
            
c     
c     minimum image convention for link cells
            
            call lnkimg(jx,ilx,cx)
            call lnkimg(jy,ily,cy)
            call lnkimg(jz,ilz,cz)
c     
c     index of neighbouring cell
            
            jc =1+jx+ilx*(jy+ily*jz)
            jj=lct(jc)
            
c     
c     ignore if empty of water
            
            if(jj.gt.0) then
              
  200         continue
c     
c     solvent id for site ii
              
              ik = ii-natm
              
              if(ic.eq.jc) jj=link(ii)
              if(jj.gt.0) then
                
  300           continue
                
c     
c     distance in real space : minimum image applied
                
                
                sxd = ssx(jj)-ssx(ii)+cx
                syd = ssy(jj)-ssy(ii)+cy
                szd = ssz(jj)-ssz(ii)+cz
                
                xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
                yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
                zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                
                rsq = xd*xd + yd*yd + zd*zd
                
c     
c     test of distance
                
                if(rsq.lt.tolnce) then
c     
c     make sure sites are not on same water molecule!

                  jk = jj-natm

                  if(ik.gt.0) then
                    if(jk.gt.0) then
                      if(ik.ne.jk) then

                        irem = irem+1
                        lstrem(irem) = max(ik,jk)
                    
                      endif
                    endif
                  endif

                endif
                
                jj=link(jj)
                if(jj.gt.natm) goto 300
                
              endif
              
              jj=lct(jc)
              ii=link(ii)
              
              if(ii.gt.natm) goto 200
              
            endif
            
          enddo
          
        endif
        
        ix=ix+1
        if(ix.ge.ilx) then
          
          ix=0
          iy=iy+1
          
          if(iy.ge.ily) then
            
            iy=0
            iz=iz+1
            
          endif
          
        endif
        
      enddo
c
c     remove unwanted water molecules
 
      call cleanup(natm,iplus,irem,lstrem,xxx,yyy,zzz)

      return 
      end

      subroutine shellsort(n,list)

c***********************************************************************
c
c     dlpoly shell sort routine. 
c     Sorts an array of integers into ascending order
c
c     copyright daresbury laboratory 1993
c     author - t.forester   november 1993
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      dimension list(*)

c
c     set up sort

      if(n.gt.1) then

c     number of lists
         nl = n/2

c     iterate shell sort

   10  do nn = 1,nl
c
c     begin insertion sort on nnth list
            
            do i = nn+nl,n,nl

               imax = list(i)
               ix = i
c
c     find location for insertion
               
               j = i
  100          j = j-nl

               if(j.lt.1) goto 110

               if (list(j).gt.imax) then
                     
                  ix = j
                     
               else
               
                  j =1

               endif
                  
               goto 100
  110           continue
               
c
c     insert in index array

               do j = i,ix+nl,-nl
                  
                  list(j) = list(j-nl)
                  
               enddo
               
               list(ix) = imax

            enddo

         enddo
         
         nl = nl/2
         if(nl.gt.0) goto 10
         
      endif

      return
      end

      subroutine dcell(aaa,bbb)

c
c***********************************************************************
c
c     dl_poly subroutine to calculate the dimensional properies of
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
      
      subroutine  imgchk
     x   (imcon,natm,iplus,tolnce,lstrem,
     x   ssx,ssy,ssz,xxx,yyy,zzz,cell,xd,yd,zd)
      
c     
c**********************************************************************
c     
c     dlpoly utility subroutine to remove additional solvent from
c     unusual (periodic) bondaries: TO, RD, ellipsoidal
c     
c     copyright daresbury laboratory 1994
c     author      t.forester october 1994
c     
c**********************************************************************
c     
      implicit real*8(a-h,o-z)
      
      logical check
      
      dimension ssx(*),ssy(*),ssz(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension cell(*)
      dimension xd(*),yd(*),zd(*)
      dimension lstrem(*)
      
      natm1 = natm+iplus
      irem = 0
c     
c     Truncated octahedral (TO)
      
      if(imcon.eq.4) then
c     
c     check that solvent is in central box
        
        do i = natm+1,natm1
          
          sx = abs(ssx(i) -0.5d0)
          sy = abs(ssy(i) -0.5d0)
          sz = abs(ssz(i) -0.5d0)
          
          if(sx+sy+sz.gt.0.75d0) then
            
            irem = irem+1
            lstrem(irem) = i-natm
            
          endif
          
        enddo
        
      elseif(imcon.eq.5) then
c     
c     Rhombic dodecahedral
        
        do i = natm+1,natm1
          
          sx = abs(ssx(i) - 0.5d0)
          sy = abs(ssy(i) - 0.5d0)
          sz = abs(ssz(i) - 0.5d0)
          
          if(sx+sy+2.d0*sz.gt.1.d0) then
            
            irem = irem+1
            lstrem(irem) = i-natm
            
          endif
          
        enddo
        
      elseif(imcon.eq.10) then
c     
c     ellipsoidal shell
        
        do i = natm+1,natm1
          
          sx = abs(ssx(i) - 0.5d0)
          sy = abs(ssy(i) - 0.5d0)
          sz = abs(ssz(i) - 0.5d0)
          
          if(sx*sx+sy*sy+sz*sz.gt.0.25d0) then
            
            irem = irem+1
            lstrem(irem) = i-natm
            
          endif
          
        enddo
        
      endif
c     
c     remove unwanted solvent
      
      call cleanup(natm,iplus,irem,lstrem,xxx,yyy,zzz)
      
      natm1 = natm+iplus
      irem = 0
      
      do i = natm+1,natm1-1
c     
c     solvent number
        
        ik = i-natm
        
        ii = 0
        
        do j = i+1,natm1
          
          ii=ii+1
          xd(ii) = xxx(i)-xxx(j)
          yd(ii) = yyy(i)-yyy(j)
          zd(ii) = zzz(i)-zzz(j)
          
        enddo
c     
c     periodic image
        
        call images(imcon,0,1,ii,cell,xd,yd,zd)
        
        check = .true.
        
        do j = 1,ii
          
          rsq = xd(j)**2+ yd(j)**2 +zd(j)**2
          if(rsq.lt.tolnce) check = .false.
          
        enddo
        
        if(.not.check) then
          irem = irem+1
          lstrem(irem) = ik
        endif
        
      enddo
c     
c     remove unwanted waters
      
      call cleanup(natm,iplus,irem,lstrem,xxx,yyy,zzz)
      
      natm1 = natm+iplus
      
      return
      end

      subroutine  cleanup(natm,iplus,irem,lstrem,xxx,yyy,zzz)

      implicit real*8(a-h,o-z)
      dimension xxx(*),yyy(*),zzz(*)
      dimension lstrem(*)
c     
c     sort into order
      
      call shellsort(irem,lstrem)
c     
c     remove redundancies from list
      
      ii = 0
      do i = irem,2,-1
        
        if(lstrem(i-1).eq.lstrem(i)) then
          
          ii=ii+1
          
          do jj = i,irem
            lstrem(jj)=lstrem(jj+1)
          enddo
          
        endif
        
      enddo
      
      irem = irem-ii
c     
c     go through list and remove all the ones we don't want
      
      do ik = irem,1,-1
        
        j = natm + lstrem(ik)
        ntot=natm+(iplus)
        
        do i = j,ntot-1
          
          xxx(i) = xxx(i+1)
          yyy(i) = yyy(i+1)
          zzz(i) = zzz(i+1)
          
        enddo

        iplus = iplus-1        

      enddo

      return
      end










