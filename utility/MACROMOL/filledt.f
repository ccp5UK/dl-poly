      program filledt
      
c***********************************************************************
c     
c     DL_POLY utility to add missing atoms to a structure.
c     files are in dl_poly edit.out style
c     
c     author t. forester  feb 1996
c     copyright daresbury laboratory 1996
c     
c     itt
c     2010-10-30 17:20:53
c     1.3
c     Exp
c     
c***********************************************************************
      
      implicit real*8(a-h,o-z)
      parameter (mxatm = 50000, mxtmls=10, mxlist=10)
      parameter (mxres=50,mxsit=30, mxunt=20000, mxbond = 50000)
c     
c     information arrays
      character*3 resnam(mxres),untnam(mxunt),res
      character*4 sitnam(mxsit,mxres),ffnam(mxsit,mxres),h,h1
      dimension join(mxsit,mxres),kkres(mxunt)
      dimension nsite(mxres),count(mxsit),index(mxsit),count1(mxsit)
      dimension rept(mxsit,mxres),qsit(mxsit,mxres)
      dimension nhsit(mxsit,mxres)
c     
c     system arrays
      logical yes ,safe
      character*4 ambnam(mxatm),name(mxatm)
      dimension xx(mxatm),yy(mxatm),zz(mxatm),qq(mxatm)
      dimension x(4),y(4),z(4)
      dimension ibda(mxbond),ibdb(mxbond)
      dimension nbds(mxatm)
      dimension hx(3),hy(3),hz(3)
      character*80 title
      dimension nh(mxsit),nexatm(mxatm,mxlist)
      dimension molatm(mxtmls),nunt(mxtmls),lres(mxunt)
      character*40 input,output
      
      kres = 0
      iunt = 0
      imol = 0
c     
c     read database file
      
      call data_read(sitnam,ffnam,resnam,mxsit,mxres,
     x     nsite,rept,nres,join,qsit,nhsit)
      
c     
c     input data in default file (*)
      
      ird1 = 9
      write(*,'(/,a)') 'name of file to be processed ?'
      read(*,'(a40)') input
      open(ird1,file=input)

      write(*,'(/,a)') 'name of output file ?'
      read(*,'(a40)') output
      open(11,file=output)
c
c     read the header record

      read(ird1,'(a80)') title
      write(11,'(a80)') title

      iatm = 0
      
 202  call bndtitle(ird1,ierr,untnam,iunt)
      
      if (ierr.eq.1) goto 201
c     
c     new molecule
      if(ierr.eq.-999) then 
        
        imol = imol+1
        if(imol.gt.mxtmls) then
          write(*,*) ' too many molecules: ',
     x         ' make mxtmls at least ',imol
          call exit(0)
        endif
        
        if(imol.gt.1) then
          molatm(imol-1) = iatm
          nunt(imol-1) = kres
        endif
      endif
      
      kres = kres+1
      call bnddata(ird1,iatm,ibnd,kres,
     x     name,ambnam,ibda,ibdb,ires,xx,yy,zz,qq)
      
      kkres(kres) = iatm
      
      goto 202
      
 201  continue
      
      if(imol.gt.0) then
        molatm(imol) = iatm
        nunt(imol) = kres
      endif
      
      if(ibnd.gt.mxbond) then
        write(*,*) 'too many bonds: make mxbond at least ',ibnd
      endif
c     
c     construct connectivity table 
      
      do i = 1,mxatm
        
        nbds(i) = 0
        
        do j = 1,mxlist
          
          nexatm(i,j)  = 0
          
        enddo
        
      enddo
c     
      ibig = 0
      do i = 1,ibnd
c     
c     counters for connectivity
        
        ia = ibda(i)
        ib = ibdb(i)
        
        if(ia.gt.0 .and. ib.gt.0) then
          nbds(ia) = nbds(ia) + 1
          nbds(ib) = nbds(ib) + 1
c     
c     connectivity table
          
          if((nbds(ia).le.mxlist).and.(nbds(ib).le.mxlist))then
            nexatm(ia,nbds(ia)) = ib
            nexatm(ib,nbds(ib)) = ia
          else
            ibig = max(ibig,nbds(ia),nbds(ib))
          endif
        endif
        
      enddo
c     
c     improve connectivity ...
      
c     
c     check residues
      
      i = 0
      kk1 = 0
      do ires = 1,iunt
        
        res = untnam(ires)
        
        yes=.false.
        do jres = 1,nres
          if(resnam(jres)(1:3).eq.res) then
            lres(ires) = jres
            lres1 = jres
            yes=.true.
          endif
        enddo
        if(.not.yes) then
          write(*,*) 'residue ',ires,res,' not in database '
          call exit(0)
        endif
      
        kk0 = kk1+1
        if(ires.gt.1) kk0 = kkres(ires-1)+1
        kk1 = kkres(ires)

        iatm1 = iatm
        iloop =1
          
        call concheck
     x       (ambnam,kk0,kk1,iatm1,ibnd,nhsit,lres1,mxsit,mxres,
     x       iloop,mxatm,xx,yy,zz,ibda,ibdb,nbds,nexatm,res,
     x       name,sitnam,nsite)
          
      enddo
      
c     
c     safety check
      if(ibig.gt.0) then 
        write(*,*) ' too many neighbours : make mxlist at least ',ibig
        call exit(0)
      endif

c
c     remove zero entries from nexatms

        do i = 1,iatm
          nn = nbds(i)
          i9 =0
 1000     i9 = i9+1
 1010     if(i9.gt.nn) goto 1020

          if(nexatm(i,i9).le.0) then
            do j9 = i9,nn-1
              nexatm(i,j9) = nexatm(i,j9+1)
            enddo
            nn = nn-1
            goto 1010
          endif
          goto 1000
 1020     nbds(i) = nn
        enddo
c     
c     check residues
      
      i = 0
      kk1 = 0
      do ires = 1,iunt
        
        res = untnam(ires)
        
        yes=.false.
        do jres = 1,nres
          if(resnam(jres)(1:3).eq.res) then
            lres(ires) = jres
            lres1 = jres
            yes=.true.
          endif
        enddo
        if(.not.yes) then
          write(*,*) 'residue ',ires,res,' not in database '
          call exit(0)
        endif

        kk0 = 1
        if(ires.gt.1) kk0 = kkres(ires-1)+1
        kk1 = kkres(ires)
        
        do i9 = 1,mxsit
          count(i9) = 0.d0
          count1(i9) = 0.d0
          index(i9) = 0
        enddo
        
        do i = kk0,kk1
          i9=0
          call data_find
     x         (lres1,i9,name(i),res,sitnam,count,mxres,
     x         mxsit,nsite)
        enddo

        safe = .true.
        do iloop = 1,3

          i0 = kk0-1
          call res_check
     x      (safe,res,sitnam,ires,lres1,count,rept,
     x       nsite,mxres,mxsit,nh,iloop,ffnam,nbds,i0,name,kk0,kk1)
          
c$$$          m1 = min(20,nsite(lres1))
c$$$          write(*,'(a,20i4)')'count',(nint(count(k9)),k9=1,m1)
c$$$          write(*,'(a,20i4)')'rept ',(nint(rept(k9,lres1)),k9=1,m1)
c$$$          write(*,'(a,20i4)')'nh   ',(nh(k9),k9=1,m1)

          if(.not.safe) then
            
            do i9  = 1,nsite(lres1)

              if(nh(i9).gt.0) then
                i = i0
                if(count(i9).gt.1.1d0) then
                  write(*,*) 'error count(',i9,') ',count(i9),i9
                  call exit(0)
                endif

                kk1 = kkres(ires)
                call site_find(i,kk0,kk1,name,sitnam(i9,lres1),res)
c     
c     locate site in backbone

                ia = 0
                ib = 0
                ic = 0
                if(nbds(i).ge.1) ia = nexatm(i,1)
                if(nbds(i).ge.2) ib = nexatm(i,2)
                if(nbds(i).ge.3) ic = nexatm(i,3)
c     
c     x-h bond length : default length is 1.0 angstrom
                bl = 1.d0
                if(sitnam(i9,lres1)(1:1).eq.'C') bl = 1.085d0
                if(sitnam(i9,lres1)(1:1).eq.'N') bl = 1.010d0
                if(sitnam(i9,lres1)(1:1).eq.'O') bl = 0.960d0
                if(sitnam(i9,lres1)(1:1).eq.'S') then
                  bl = 1.336d0
                  if(sitnam(i9+iloop,lres1)(1:2).eq.'LP') bl = 0.679d0
                endif
c     
c     assign force field and site name of atoms to be added

                h = ffnam(i9+iloop,lres1)
                h1= sitnam(i9+iloop,lres1)
                qnew = qsit(i9+iloop,lres1)
                
                write(*,'(a,i3,3a,i6,3a,i3)') 'adding ',nh(i9),
     x               ' ',h,' to site ',
     x               i,' (',name(i),')'
                
                
                do jj = 1,4
                  if(jj.eq.1) kk = i
                  if(jj.eq.2) kk = ia
                  if(jj.eq.3) kk = ib
                  if(jj.eq.4) kk = ic
                  
                  if(kk.gt.0) then
                    x(jj) = xx(kk)
                    y(jj) = yy(kk)
                    z(jj) = zz(kk)
                  endif
                enddo
                
                if(nh(i9).eq.1) then
                  
                  if(ic.gt.0) then

                    call  hadd1_3(x,y,z,hx,hy,hz,bl)
                    
                  elseif(ib.gt.0) then
                    
                    angl = 0.d0
                    call  hadd1_2(x,y,z,hx,hy,hz,bl,angl)
                    
                  else
                    
                    angl = 109.48d0
                    call  hadd1_1(x,y,z,hx,hy,hz,bl,angl)
                    
                  endif
                  
                else if(nh(i9).eq.2) then
                  
                  if(ib.gt.0) then
                    
                    angl = 109.48d0
                    if(bl.lt.0.8d0) angl = 160.d0
                    call  hadd2_2(x,y,z,hx,hy,hz,bl,angl)
                    
                  elseif(ia.gt.0) then

                    angl = 120.d0
                    if(bl.lt.0.8d0) angl = 96.d0
                    call  hadd2_1(x,y,z,hx,hy,hz,bl,angl,180.d0)

                  else
c
c     hydrogens onto water essentially random orientation
                    
                    bl = 1.d0
                    x(2) = 0.d0
                    y(2) = 0.d0
                    z(2) = 0.d0
                    angl = 109.48
                    call hadd2_1(x,y,z,hx,hy,hz,bl,angl,120.d0)
                    
                  endif
                  
                else if(nh(i9).eq.3) then
                  
                  angl = 109.5d0
                  call  hadd3_1(x,y,z,hx,hy,hz,bl,angl)  
                  
                else if(nh(i9).eq.4) then
                  
                  call  hadd4_0(x,y,z,hx,hy,hz,bl)  
                  
                endif
                
c     
c     readjust book keeping to make room for iadd h's
                
                iadd = nh(i9)
                count(i9+iloop) = count(i9+iloop) + dble(iadd)
                
                do ir = ires,kres
                  kkres(ir) = kkres(ir)+iadd
                enddo
                kk1 = kk1+iadd
                
                do j = 1,iatm
                  do jk = 1,nbds(j)
                    if(nexatm(j,jk).gt.i)
     x                   nexatm(j,jk) = nexatm(j,jk)+iadd
                  enddo
                enddo

                do j = iatm,i+1,-1
                  
                  xx(j+iadd) = xx(j)
                  yy(j+iadd) = yy(j)
                  zz(j+iadd) = zz(j)
                  qq(j+iadd) = qq(j)
                  ambnam(j+iadd) = ambnam(j)
                  name(j+iadd) = name(j)
                  nbds(j+iadd) = nbds(j)
                  
                  do jk = 1,nbds(j)
                    nexatm(j+iadd,jk) = nexatm(j,jk)
                  enddo
                  
                enddo
                
                iatm = iatm+iadd

c     
c     place h's in structure
                
                do j = 1,iadd
                  
                  xx(i+j) = hx(j)
                  yy(i+j) = hy(j)
                  zz(i+j) = hz(j)
                  qq(i+j) = qnew
                  
                  name(i+j) = h1
                  ambnam(i+j) = h
                  
                  nexatm(i+j,1) = i
                  nexatm(i+j,2) = 0
                  nexatm(i+j,3) = 0
                  nbds(i+j) = 1

                  nbds(i) = nbds(i) + 1
                  nexatm(i,nbds(i))= i+j
                  
                enddo
                
              endif
              
            enddo
            
          else
            goto 150
          endif
        
        enddo
          
 150    continue
      enddo
c     
c     print out 'connect_dat' style
c     
c     adjust nunt array
      kres = 0
      do imols = 1,imol
        nunt(imols) = nunt(imols) - kres
        kres = kres + nunt(imols)
      enddo
      
      kres = 0
      kk1 = 0
      i = 0
      do imols = 1,imol
        
        write(11,'(//,6x,a,//)')'MOLECULE'
        
        do junt = 1,nunt(imols)
          
          kres = kres + 1 
          lres1 = lres(kres)
          write(11,'(/,5x,a,i5,2a)') 'RESIDUE',kres,' =  ',resnam(lres1)
          write(11,'(/,5x,a,/)') 'BOND ARRAY'
          
          kk0 = kk1+1
          kk1 = kkres(kres)
          
          do i = kk0,kk1
            
            ia = nexatm(i,1)
            if(ia.gt.i.or.ia.le.0)  then 
              ia = -99
              do jk = 1,nbds(i)
                in = nexatm(i,jk)
                if((in.gt.0).and.(in.lt.i)) ia = in
              enddo
            endif
            
            write(11,99) i,ia,name(i),ambnam(i),
     x           xx(i),yy(i),zz(i),qq(i)
            
 99         format(2i5,3x,a4,2x,a4,6x,3f10.4,f8.4)
            
          enddo
          
        enddo
        
      enddo
      
      end

      subroutine bndtitle(ird1,ierr,untnam,iunt)

      implicit real*8(a-h,o-z)
      character*60 oneline
      character*3 untnam(*)

      ierr= 0
   10 read(ird1,'(a60)',err=20,end=20) oneline

      if (oneline(7:14).eq."MOLECULE") ierr=-999
      if (oneline(6:12).eq."RESIDUE") then
        do jk = 19,30
          if(oneline(jk:jk).eq.'=') oneline(jk:jk)=' '
        enddo
        oneline(1:12)=oneline(19:30)
        call strip(oneline,12)
        call lowcase(oneline,3)
        iunt = iunt+1
        untnam(iunt) = oneline(1:3)
      endif
      if(oneline(1:12).eq.'CELL_VECTORS') then
        write(11,'(a)') oneline(1:12)
        do jk = 1,3
          read(ird1,'(a60)',err=20,end=20) oneline
          write(11,*) oneline
        enddo
      endif
      if (oneline(6:15).ne."BOND ARRAY") goto 10

      read(ird1,*)
      return 
c
c     end of file reached

   20 ierr = 1

      return
      end

      subroutine bnddata(ird1,iatm,ibnd,kres,
     x     name,ambnam,ibda,ibdb,ires,xx,yy,zz,qq)

      implicit real*8(a-h,o-z)
      
      character*4 ambnam(*),an,at
      character*4 name(*)
      dimension xx(*),yy(*),zz(*),qq(*)
      dimension ibda(*),ibdb(*)

   10 read(ird1,'(2i5,3x,a4,2x,a4,6x,3f10.4,f8.4)',err=20,end=30)
     x      ia,ib,at,an,x1,y1,z1,q1

      if (ia.gt.0)  then 

        ambnam(ia) = an
        qq(ia) = q1

        if(ia.gt.iatm) then

          name(ia) = at
          iatm = ia
          xx(ia) = x1
          yy(ia) = y1
          zz(ia) = z1
          ires = kres

        endif
         
      endif
      
      if (ib.gt.0) then

         ibnd = ibnd + 1

         ibda(ibnd) = ia
         ibdb(ibnd) = ib

      endif
      
      goto 10

 20   backspace(ird1)
 30   return
      end
      
      subroutine concheck
     x     (ambnam,kk0,kk1,iatm1,ibnd,nhsit,lres1,mxsit,mxres,
     x     iloop,mxatm,x,y,z,ibda,ibdb,nbds,nexatm,res,
     x     name,sitnam,nsite)

      implicit real*8(a-h,o-z)
      
      logical now,found,safe
      dimension x(*),y(*),z(*)
      character*4 ambnam(*),name(*),sitnam(mxsit,mxres)
      dimension ibda(*),ibdb(*),nbds(*),nexatm(mxatm,*)
      dimension rmin(4), imin(4)
      dimension count1(50)
      dimension nhsit(mxsit,mxres),nsite(*)
      character*3 res
c     
c     checks that the correct number of connections have been found for each 
c     atom type
      
c     set up tolerances for bond lengths
      if(iloop.eq.1) then
        rl1 = 1.2d0**2
        rl2 = 1.6d0**2
      elseif(iloop.eq.2) then
        rl1 = 0.5d0**2
        rl2 = 1.6d0**2
      elseif(iloop.eq.3) then
        rl1 = 0.5d0**2
        rl2 = 3.0d0**2
      endif


      do i = kk0,kk1
        
        icn = 0
        call nconnect(icn,ambnam(i))
c     
c     search for missing connections - taken to be the shortest
c     remaining distances to the central atom.
        
        iter = icn - nbds(i)
        do ll = 1,nbds(i)
          l2 = nexatm(i,ll)
          if(ambnam(l2)(1:1).eq.'H') iter = iter+1
          if(ambnam(l2)(1:2).eq.'LP') iter = iter+1
        enddo
c
c     locate site in residue

        call data_find
     x     (lres1,j1,name(i),res,sitnam,count1,mxres,
     x     mxsit,nsite)
        iter = iter - nhsit(j1,lres1)
c
c     apply known patches ...

c     i) N terminal AA has  R-N-H_3
        if(iter.eq.-1.and.ambnam(i).eq.'N   ') then
          if(iloop.eq.3)
     x         write(6,*) ' warning: 4 coordinate Nitrogen ',i
          iter =0
        endif
c     ii) OH in tyo has missing H
        if(iter.eq.1.and.ambnam(i).eq.'OH  ') then
          if(iloop.eq.3) 
     x         write(6,*) ' warning: 1 coordinate Hydroxy Oxygen ',i
          iter =0
        endif
        
        
        if(iter.gt.0) then
          
          do kk = 1,4
            rmin(kk) = rl2
            imin(kk) = 0
          enddo
          
          do j = kk0,iatm1
            
            now = .true.
            if(iloop.ne.2.and.ambnam(i)(1:1).eq.'H') now = .false.
            if(iloop.ne.2.and.ambnam(j)(1:1).eq.'H') now = .false.

            if (i.eq.j) now =.false.
            call nconnect(jcn,ambnam(j))
            if(jcn.le.nbds(j)) now = .false.
            
            if(now) then

              do k = 1,nbds(i)
              
                if(nexatm(i,k).eq.j) now = .false.
              
              enddo

            endif
            if(now) then
              
              r = (x(i)-x(j))**2+(y(i)-y(j))**2 + (z(i)-z(j))**2
              
              safe =.false.
              do k1 = 1,4
                
                if(r.lt.rmin(k1)) then
                  
                  if(r.gt.rl1) then
                    safe = .true.
                    do k2 = 4,k1+1,-1
                    
                      rmin(k2) = rmin(k2-1)
                      imin(k2) = imin(k2-1)
                    
                    enddo
                  
                    imin(k1) = j
                    rmin(k1) = r

                  endif
                  
                endif
                
                if(safe) goto 10
                
              enddo
              
   10       endif
            
          enddo
          
c     
c     update connection list 
          
          do kk = 1,iter
            ij = imin(kk)
            rr = sqrt(rmin(kk))
            if(ij.gt.0) then
c     
c     insert near entry for i or ij
              
              i1 = min(i,ij)
              i3 = ibnd+1
              found = .false.
              do i2 = 1,ibnd
                
                if(ibda(i2).ge.i1.or.ibdb(i2).ge.i1) then
                  i3 = i2
                  found = .true.
                endif
                if(found) goto 20
              enddo
   20         do i2 = ibnd+1,i3+1,-1
                ibda(i2) = ibda(i2-1)
                ibdb(i2) = ibdb(i2-1)
              enddo
              if(ij.gt.0.and.i.gt.0) then
               ibnd = ibnd + 1
               ibda(i3) = i
               ibdb(i3) = ij
               nbds(i) = nbds(i) + 1
               nbds(ij) = nbds(ij) +1
               nexatm(i,nbds(i)) = ij
               nexatm(ij,nbds(ij)) = i
c
c     write out the longer bonds only

               if(rr.gt.1.7d0) write(6,'(a,i6,3a,i6,3a,f9.4)')
     x        'adding neighbour to ',i,
     x        '(',ambnam(i),') : ',ij,'(',ambnam(ij),') : dist= ',rr
              
               endif
             endif
          enddo
          
        elseif (nbds(i).gt.icn) then
c     
c     too many connections for this atom
          
          if(ambnam(i).ne.'OW  ') 
     x         write(6,*)'too many neighbours for ',i,'(',ambnam(i),')',
     x      ': found ',nbds(i),': should be ',icn
          
c$$$  write(6,*) ' neighbours and distances are : '
c$$$  do k = 1,nbds(i)
c$$$  
c$$$  j = nexatm(i,k)
c$$$  r = (x(i)-x(j))**2+(y(i)-y(j))**2 + (z(i)-z(j))**2
c$$$  r = sqrt(r)
c$$$  
c$$$  write(6,*) j,'(',ambnam(j),') ',r
c$$$  
c$$$  enddo
c     
c     remove largest bonds
          
          kremove = nbds(i)-icn
          do krem = 1,kremove
            
            rbig=0.d0
            do k = 1,nbds(i)
              
              j = nexatm(i,k)
              r = (x(i)-x(j))**2+(y(i)-y(j))**2 + (z(i)-z(j))**2
              if(r.gt.rbig) then
                ibig = k
                rbig=r
              endif
              
            enddo
            
            j = nexatm(i,ibig)
c     
c     remove j from entry in i
            
            found = .false.
            do i1 = 1,nbds(i)
              
              if(nexatm(i,i1).eq.j) then
                
                found=.true.
                do i2 = i1,nbds(i)-1
                  nexatm(i,i2) = nexatm(i,i2+1)
                enddo
              endif
            enddo
            if(found) then
              nexatm(i,nbds(i))=0
              nbds(i) = nbds(i)-1
ctrf              write(6,*) 'removed ',j,' from list ',i
            else
              write(6,*) ' failed to find ',j,' in list ',i
              write(6,*) nbds(i), (nexatm(i,i1),i1=1,5)
            endif
c     
c     remove i from entry in j
            
            found = .false.
            do j1 = 1,nbds(j)
              
              if(nexatm(j,j1).eq.i) then
                
                found=.true.
                do j2 = j1,nbds(j)-1
                  nexatm(j,j2) = nexatm(j,j2+1)
                enddo
              endif
            enddo
            if(found) then
              nexatm(j,nbds(j))=0
              nbds(j) = nbds(j)-1
ctrf              write(6,*) 'removed ',i,' from list ',j
            else
              write(6,*) ' failed to find ',i,' in list ',j
              write(6,*) nbds(j), (nexatm(j,j1),j1=1,5)
            endif
            
c     
c     remove i,j from list of bonds in system
            
            found=.false.
            do k = 1,ibnd
              if(ibda(k).eq.i) then
                if(ibdb(k).eq.j) then
                  
                  do k1 = k,ibnd-1
                    ibda(k1) = ibda(k1+1)
                    ibdb(k1) = ibdb(k1+1)
                  enddo
                  
                  ibda(ibnd) = 0
                  ibdb(ibnd) = 0
                  found=.true.
                  
                endif
                
              elseif(ibda(k).eq.j) then
                if(ibdb(k).eq.i) then
                  
                  do k1 = k,ibnd-1
                    ibda(k1) = ibda(k1+1)
                    ibdb(k1) = ibdb(k1+1)
                  enddo
                  
                  ibda(ibnd) = 0
                  ibdb(ibnd) = 0
                  found=.true.
                  
                endif
                
              endif
              
            enddo
            
            if(found) then
              ibnd = ibnd -1
            else
              write(6,*) 'failed to locate bond ',i,j 
            endif
            
          enddo
          
        endif
        
      enddo
      return
      end

      subroutine lowcase(string,length)

c***********************************************************************
c
c     DL_POLY routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c
c***********************************************************************

      character*255 string
      character*1 letter

      do i = 1,min(length,255)

         letter = string(i:i)

         if(letter.eq.'A') letter = 'a'
         if(letter.eq.'B') letter = 'b'
         if(letter.eq.'C') letter = 'c'
         if(letter.eq.'D') letter = 'd'
         if(letter.eq.'E') letter = 'e'
         if(letter.eq.'F') letter = 'f'
         if(letter.eq.'G') letter = 'g'
         if(letter.eq.'H') letter = 'h'
         if(letter.eq.'I') letter = 'i'
         if(letter.eq.'J') letter = 'j'
         if(letter.eq.'K') letter = 'k'
         if(letter.eq.'L') letter = 'l'
         if(letter.eq.'M') letter = 'm'
         if(letter.eq.'N') letter = 'n'
         if(letter.eq.'O') letter = 'o'
         if(letter.eq.'P') letter = 'p'
         if(letter.eq.'Q') letter = 'q'
         if(letter.eq.'R') letter = 'r'
         if(letter.eq.'S') letter = 's'
         if(letter.eq.'T') letter = 't'
         if(letter.eq.'U') letter = 'u'
         if(letter.eq.'V') letter = 'v'
         if(letter.eq.'W') letter = 'w'
         if(letter.eq.'X') letter = 'x'
         if(letter.eq.'Y') letter = 'y'
         if(letter.eq.'Z') letter = 'z'

         string(i:i) = letter

      enddo

      return
      end

      subroutine strip(string,length)

c***********************************************************************
c
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c
c***********************************************************************

      character*255 string
      
      imax = min(length,255)
      do i = 1,imax

         if(string(1:1).eq.' ') then

            do j = 1,imax-1

               string(j:j) = string(j+1:j+1)

            enddo

            string(imax:imax) = ' '

         endif

      enddo

      return
      end
      function dblstr(word,len,lst)
c     
c***********************************************************************
c     
c     dl_poly function for extracting double precisions from a 
c     character string. 
c     modified from dl_poly function intstr
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     modified  - t. forester april 1994
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     double precision string
c
c***********************************************************************
c     
      
      implicit real*8 (a-h,o-z)
      
      character*1 n(0:9),word(len),ksn,dot,d,e,work(100)
      logical flag,ldot,start
      
      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/
      
      sn=1.d0
      ksn='+'
      ten = 10.d0
      one = 1.d0
      
      call lowcase(word,len)
      
      dblstr=0.d0
      iexp = 0
      idum =0
      start=.false.
      ldot = .false.
      
      do lst=1,len
        
        flag=.false.
        
        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            dblstr= ten*dblstr+one*dble(j)
            flag=.true.
            start=.true.
            
          endif
          
        enddo
        
        if(dot.eq.word(lst)) then
          
          flag=.true.
          ten = 1.d0
          ldot =.true.
          start=.true.
          
        endif

        if(flag.and.ksn.eq.'-') sn=-1.d0
        
        ksn=word(lst)
        
        if(ldot) then
          
          one = one/10.d0
          
        endif
        
        if(start) then
          
          if(d.eq.ksn.or.e.eq.ksn) then
            
            do i = 1,len-lst
              work(i) = word(i+lst)
            enddo
            iexp = intstr(work,len-lst,idum)
            go to 100

          endif

          if(.not.flag)go to 100
          
        endif
        
      enddo
      
  100 dblstr = sn*dblstr*(10.d0**iexp)
      lst = lst +idum
      
      return
      end
      function intstr(word,len,lst)

c     
c***********************************************************************
c     
c     dl_poly function for extracting integers from a 
c     character string
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c     
c***********************************************************************
c     
      
      character*1 n(0:9),word(len),ksn
      logical flag,final
      
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      ksn='+'
      intstr=0
      final=.false.
      
      do lst=1,len
        
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            flag=.true.
            
          endif
          
        enddo

        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

        if(flag)then

          final=.true.

        else

          if(final)then

            intstr=isn*intstr
            return

          endif

        endif

      enddo

      intstr=isn*intstr
      
      return
      
      end

      subroutine data_read(sitnam,ffnam,resnam,mxsit,mxres,
     x     nsite,rept,nres,join,qsit,nhsit)

      implicit real*8(a-h,o-z)

      logical dna,rna,yes
      character*3 bac
      dimension nsite(*),join(mxsit,mxres)
      dimension rept(mxsit,mxres),qsit(mxsit,mxres)
      character*4 sitnam(mxsit,mxres),ffnam(mxsit,mxres)
      character*3 resnam(mxres)
      dimension nhsit(mxsit,mxres)
      character*40 oneline

      write(*,*) 'enter name of residue map file (e.g. resmap.amb ',
     x     'or resmap_u.amb) '
      read(*,'(a40)') oneline
      open(20,file=oneline,status='old')

      read(20,*)
      read(20,*) nres

      do i = 1,nres
        read(20,'(10x,a3)') resnam(i)
        call lowcase(resnam(i),3)
        call strip(resnam(i),3)

        read(20,'(7x,i10)') nsite(i)

        do j = 1,nsite(i)
          read(20,'(4x,i4,2a4,2f8.0)') 
     x         join(j,i),sitnam(j,i),ffnam(j,i),qsit(j,i),
     x         rept(j,i)
          rept(j,i) = max(1.d0,rept(j,i))
          nhsit(j,i) = 0
          if((sitnam(j,i)(1:1).eq.'H').or.
     x         (sitnam(j,i)(1:2).eq.'LP')) then

            do j1 = j-1,1,-1
              if((sitnam(j1,i)(1:1).ne.'H').and.
     x          (sitnam(j1,i)(1:2).ne.'LP')) then
                nhsit(j1,i) = nhsit(j1,i) +nint(rept(j,i))
                goto 100
              endif
            enddo
          endif
 100      continue

        enddo
      enddo


 199  yes = .true.
      write(*,'(/,a)')
     x     'do you have RNA units ? (y/n)'
      read(*,'(a40)') oneline
      call lowcase(oneline,40)
      call strip(oneline,40)
      if(oneline(1:1).eq.'n') then
        dna = .true.
        bac = 'dna'
        rna = .false.
      elseif (oneline(1:1).eq.'y') then
        rna =.true.
        bac='rna'
        dna=.false.
      else
        write(6,*) ' you must answer "y" or "n"'
        write(6,*) 
        yes = .false.
      endif
      if(.not.yes) goto 199
c
c     modify dna/rna residues in database to include backbone

      if(dna.or.rna) call res_back
     x   (bac,mxres,mxsit,nres,join,
     x   resnam,sitnam,ffnam,qsit,rept,nsite)

      return 
      end
      subroutine res_check
     x  (safe,res,sitnam,ires,lres,count,rept,nsite,mxres,
     x     mxsit,nh,iloop,ffnam,nbds,i0,name,kk0,kk1)

      implicit real*8(a-h,o-z)

      character*3 res
      character*4 sitnam(mxsit,mxres),ffnam(mxsit,mxres),name(*)
      dimension count(mxsit),rept(mxsit,mxres)
      dimension nsite(mxres)
      dimension nh(*),nbds(*),index(100)
      logical safe

      safe=.true.
      
      do i = 1,nsite(lres)

        index(i) = 0
        nh(i) = 0
        if(abs(count(i)-rept(i,lres)).gt.0.1d0) then
          safe=.false.
        endif
      enddo

      if(.not.safe) then

        if(res.ne.'mtl') then

          if(iloop.eq.1) write(*,*) 'residue', ires, ' (',res,')'
      
          do i = 1,nsite(lres)

            j = nint(count(i)-rept(i,lres))
            if(j.lt.0) then
c
c     test for missing hydrogens 
c     exit if atoms other than H or lone pairs are missing!

            
              if((sitnam(i,lres)(1:1).ne.'H').and.
     x           (sitnam(i,lres)(1:1).ne.'L')) then
                if(iloop.eq.1)write(*,'(4x,a,i3,3a,i3,a)')'missing ',
     x               -j,' atoms of type ',sitnam(i,lres),'(',i,')'
              else
                j1 = i-iloop
                if(j1.gt.0) then
                  if(nint(count(j1)-rept(j1,lres)).eq.0) then
                    call nconnect(jcn,ffnam(j1,lres))
                    
                    call site_find
     x                   (j2,kk0,kk1,name,sitnam(j1,lres),res)
                    jcn = jcn - nbds(j2)
                    if(jcn+j.ge.0) nh(i-iloop) = -j
                  endif
                endif
              endif

            elseif(j.gt.0) then
              write(*,'(4x,i3,3a,i3,a)')j,' too many  atoms of type ',
     x             sitnam(i,lres),'(',i,')'
            endif

          enddo

        endif

      endif
      return
      end

      subroutine data_find
     x     (lres,jj,name,res,sitnam,count,mxres,
     x     mxsit,nsite)
      
      implicit real*8(a-h,o-z)
      character*3 res
      character*4 name
      character*4 sitnam(mxsit,mxres)
      dimension count(*)
      dimension nsite(mxres)
      
c     
c     look at as many characters as necessry to uniquely identify
c     the site
      
      k = 1
      
 10   look = 0
      do i = 1,nsite(lres)
        
        if(sitnam(i,lres)(1:k).eq.name(1:k)) then

          jj = i
          look = look+1
          
        endif
        
      enddo
      
      if(look.eq.0) then
        write(*,*) 'could not locate site "',name,'" in residue ',
     x       res(1:3)
        return
      elseif(look.eq.1) then
        count(jj) = count(jj)+1.d0
        return
      endif
      
      k = k+1
      if(k.gt.4) then
        write(*,*) look,' sites with "',name,'" in residue ',res(1:3)
        call exit(0)
      endif
      
      goto 10
      end

      subroutine site_find(i,kk0,kk1,name,sitnam,res)
      character*4 name(*),sitnam
      character*3 res
c     
c     look at as many characters as necessry to uniquely identify
c     the site
      
      k = 1
      
 10   look = 0
      do j = kk0,kk1
        
        if(sitnam(1:k).eq.name(j)(1:k)) then

          i=j
          look = look+1
          
        endif
        
      enddo
      
      if(look.eq.0) then
        write(*,*) 'could not locate site "',sitnam,'" in residue ',
     x       res
        return
      elseif(look.eq.1) then
        return
      endif
      
      k = k+1
      if(k.gt.4) then
        write(*,*) look,' sites with "',sitnam,'" in residue ',res
        call exit(0)
      endif
      
      goto 10
      end

      subroutine res_back
     x   (bac,mxres,mxsit,nres,join,
     x   resnam,sitnam,ffnam,qsit,rept,nsite)

c     subroutine to append 'backbone' data into 'residue' data
c     
      implicit real*8(a-h,o-z)

      logical yes
      character*3 bac,resnam(mxres)
      character*4 sitnam(mxsit,mxres),ffnam(mxsit,mxres)
      dimension qsit(mxsit,mxres),rept(mxsit,mxres)
      dimension nsite(mxres),join(mxsit,mxres)
c
c     locate which backbone: rna or dna

      yes=.false.
      do jres = 1,nres
        if(resnam(jres)(1:3).eq.bac) then
          lbac = jres
          yes=.true.
        endif
      enddo
      if(.not.yes) then
        write(*,*) bac,' backbone not in database '
        call exit(0)
      endif
c
c     add in backbone to nuc residue data:
c     assume backbone comes before residues...

      nadd = nsite(lbac)
      ibig = 0

      do jres = lbac+1,nres

        nnn = nsite(jres)

        if(nnn+nadd.gt.mxres) then
          ibig = max(ibig,nnn+nadd)
        else

          nsite(jres) = nsite(jres)+nadd

          do k = 1,nadd

            sitnam(k+nnn,jres) = sitnam(k,lbac)
            ffnam(k+nnn,jres) = ffnam(k,lbac)
            qsit(k+nnn,jres) = qsit(k,lbac)
            rept(k+nnn,jres) = rept(k,lbac)
            
            join(k+nnn,jres) = join(k,lbac)+nnn

          enddo

        endif

      enddo

      if(ibig.gt.0) then
        write(*,'(/,a)') 'parameter mxres is too small'
        write(*,*) ' it should be at least ',ibig
        call exit(0)
      endif

      return
      end
