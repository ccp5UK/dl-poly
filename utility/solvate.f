      program solvate
c*********************************************************************
c
c     dl_poly utility to solvate a large molecule in any solvent
c
c1. A file containing your protein in CONFIG format 2. A file containing the solvent in CONFIG format.
c
cIt is assumed that you have simulated the solvent file to equilibrate the system and that the system is big enough to accommodate the protein. In other words, the final simulation box will be the size of your solvent box. It is also assumed the box is cubic.
c
c3. Run the solvate program and answer a few simple questions. The
c   final solvated molecule is written in the CFGSOL file.
c
c     copyright daresbury laboratory
c     author w. smith sep 2002
c
c*********************************************************************

      implicit real*8(a-h,o-z)
      parameter (mega=2000000)

      character*80 title
      character*8 name(mega)
      character*40 molfil,solfil
      dimension xxx(mega),yyy(mega),zzz(mega),cell(9)
      dimension lct(mega),link(mega),nil(mega)
      dimension nix(27),niy(27),niz(27)
      data nix/ 0,-1,-1,-1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1,
     x  1, 1, 1, 0, 0, 1,-1, 1, 0,-1, 1, 0,-1/
      data niy/ 0, 0,-1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1,
     x  0, 1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1/
      data niz/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     x  0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/

c     input parameters

      write(*,'(a)')'enter name of molecule file'
      read(*,*)molfil
      write(*,'(a)')'enter name of solvent file'
      read(*,*)solfil
      write(*,'(a)')'enter solute-solvent separation'
      read(*,*)tolnce
      write(*,'(a)')'enter number of atoms in solvent molecule'
      read(*,*)nsat

      rcsq=tolnce**2

c     real molecule file

      open(7,file=molfil)

      read(7,'(a80)')title
      read(7,*)levcfg,imcon
      if(imcon.gt.0)then
        read(7,*)
        read(7,*)
        read(7,*)
      endif
      do i=1,mega

        read(7,'(a8)',end=100)name(i)
        read(7,*)xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)read(7,*)
        if(levcfg.gt.1)read(7,*)

      enddo

  100 close(7)
      nmol=i-1
      write(*,*)'number of molecule atoms',nmol

c     read solvent file

      open(7,file=solfil)
      read(7,*)
      read(7,*)levcfg,imcon
      if(imcon.gt.0)then
        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)
      endif
      do j=1,mega

        if(i.gt.mega)then
           write(*,*)'error - too many atoms in file'
           stop
        endif
        read(7,'(a8)',end=200)name(i)
        read(7,*)xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)read(7,*)
        if(levcfg.gt.1)read(7,*)
        i=i+1

      enddo

  200 close(7)
      natms=i-1
      write(*,*)'number of solvent molecules ',(natms-nmol)/nsat

c     enforce periodic boundary

      aaa=1.d0/cell(1)
      bbb=1.d0/cell(5)
      ccc=1.d0/cell(9)
      
      do i=1,natms
         
         xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
         yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
         zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))
         
      enddo

c     calculate link cells (for solvent atoms only)

      ilx=int(cell(1)/tolnce)
      ily=int(cell(5)/tolnce)
      ilz=int(cell(9)/tolnce)
      ncells = ilx*ily*ilz
      write(*,*)'number of link cells',ncells
      if(ncells.gt.mega)then
        write(*,*)'error - too many link cells'
        stop
      endif

      do i=1,mega
        lct(i)=0
        nil(i)=0
      enddo

      xdc = aaa*dble(ilx)
      ydc = bbb*dble(ily)
      zdc = ccc*dble(ilz)

      do i = nmol+1,natms
        
        ix = min(int(xdc*(xxx(i)+cell(1)/2.d0)),ilx-1)
        iy = min(int(ydc*(yyy(i)+cell(5)/2.d0)),ily-1)
        iz = min(int(zdc*(zzz(i)+cell(9)/2.d0)),ilz-1)
        icell = 1+ix+ilx*(iy+ily*iz)
        j = lct(icell)
        lct(icell)=i
        link(i)=j
        
      enddo

c     loop over molecule atoms

      do i=1,nmol

        ix = min(int(xdc*(xxx(i)+cell(1)/2.d0)),ilx-1)+1
        iy = min(int(ydc*(yyy(i)+cell(5)/2.d0)),ily-1)+1
        iz = min(int(zdc*(zzz(i)+cell(9)/2.d0)),ilz-1)+1
        
        do kc = 1,27
            
          cx = 0.d0
          cy = 0.d0
          cz = 0.d0
          jx=ix+nix(kc)
          jy=iy+niy(kc)
          jz=iz+niz(kc)
          
c     minimum image convention
            
          if(jx.gt.ilx) then
            
            jx = jx-ilx
            cx = cell(1)
            
          elseif(jx.lt.1) then
            
            jx = jx+ilx
            cx =-cell(1)
            
          endif
          
          if(jy.gt.ily) then
            
            jy = jy-ily
            cy = cell(5)
            
          elseif(jy.lt.1) then
            
            jy = jy+ily
            cy =-cell(5)
            
          endif
          
          if(jz.gt.ilz) then
            
            jz = jz-ilz
            cz = cell(9)
            
          elseif(jz.lt.1) then
            
            jz = jz+ilz
            cz =-cell(9)
            
          endif
          
c     index of neighbouring cell
          
          jc =jx+ilx*((jy-1)+ily*(jz-1))
          j=lct(jc)
          
          do while(j.gt.0)
             
c     distance in real space : minimum image applied
             
             if(nil(j).eq.0)then
                
                xd = xxx(j)-xxx(i)+cx
                yd = yyy(j)-yyy(i)+cy
                zd = zzz(j)-zzz(i)+cz
                rsq=xd*xd+yd*yd+zd*zd
                
c     test of distance
                
                if(rcsq.gt.rsq) then
                   
                   k=nmol+nsat*((j-nmol-1)/nsat)
                   if(nsat.eq.3)then
                     nil(k+1)=1
                     nil(k+2)=1
                     nil(k+3)=1
                   else
                     do m=1,nsat
                       nil(k+m)=1
                     enddo
                   endif
                   
                endif
                
             endif
            
             j=link(j)
              
          enddo
          
        enddo
        
      enddo
      
      open(8,file="CFGSOL")
      write(8,'(a80)')title
      write(8,'(2i10)')0,imcon
      write(8,'(3f20.10)')cell(1),cell(2),cell(3)
      write(8,'(3f20.10)')cell(4),cell(5),cell(6)
      write(8,'(3f20.10)')cell(7),cell(8),cell(9)
      k=0
      do i=1,natms

        if(nil(i).eq.0)then

          k=k+1
          write(8,'(a8,i10)')name(i),k
          write(8,'(3f20.10)')xxx(i),yyy(i),zzz(i)

        endif

      enddo

      write(*,*)'number of atoms in file CFGSOL',k

      close(8)

      end
