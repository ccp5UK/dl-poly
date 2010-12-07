      program jong3
c programme to test for clashes and cut out the membrane molecules c c the logic of version 1 is that a membrane molecule is removed if at c least one of its atoms clashes with the protein c
c copyright P.-L. Chau   2000/01/29
c
c version 1 altered  P.-L. Chau   2006/07/03
c
c the logic of version 2 is that a membrane molecule is removed if and c only if all its atoms clash with the protein c
c version 2 copyright P.-L. Chau   2008/04/28
c
c the logic of version 3 is that (1) if a molecule is DMPC, it is removed if c and only if at least half of its atoms clash with the protein (2) if a c molecule is water, it is removed even if one of its atoms clashes with the c protein c
c version 3 copyright P.-L. Chau   2008/05/16

      implicit none
      integer patms,mxatms,mxjong,mxmol
      parameter (patms = 20000, mxatms = 200000, mxjong = 250000)
      parameter (mxmol = 19000)

      real*8 px(patms),py(patms),pz(patms)
      real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
      real*8 cell(9),cell2(9),cmx,cmy,cmz,tram
      real*8 dxx(1),dyy(1),dzz(1),dist2,djong,djong2,status(mxmol)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,rejfrac

      integer nconfig,nconf2,nconfout,levcfg,imcon
      integer resnum(patms),memnum(mxatms)
      integer i,j,k,l,ic,icount,ix,natms,natms2
      integer levcfg2,imcon2
      integer mark,index(mxatms),ijong
      integer jonglist(mxjong),jonglist2(mxjong)
      integer maxmemnum,natpmol(mxmol),indexsum(mxmol)
      integer itarget,iassess,iclash,inoclash

      character*1 anam,elenam
      character*3 resnam(patms),memnam(mxatms)
      character*4 pdbatm(patms),mematm(mxatms)
      character*4 atmnam0(patms),atmnam1(patms)
      character*8 atmnam2(mxatms)
      character*40 cjunk,title,title2

c define units
      nconfig = 1
      nconf2 = 2
      nconfout = 8
c define clash distance
      djong = 4.0d0
      djong2 = djong * djong
c define input files: CONFIG.shr is the protein, CONFIG.mb4a is the membrane
      open(nconfig,file='CONFIG.shr',status='unknown')
      open(nconf2,file='CONFIG.mb4a',status='unknown')
c the output file is CONFIG.nj3
      open(nconfout,file='CONFIG.nj3',status='unknown')

c define rejection fraction, rejfrac. this number defines the fraction of c clashed atoms in a molecule beyond which the molecule is rejected. it is c ALWAYS a negative number. e.g., if rejfrac = -0.7, then a molecule can c have up to 70% of its atoms clashing with the protein and still be retained c by this programme.
c qualitatively, the more negative the rejection fraction, the more tolerant c the programme is with clashed molecules
      rejfrac = -0.7
c read in CONFIG.gabaar file header
      read(nconfig,'(a40)') title
      read(nconfig,*)levcfg,imcon
      if (imcon.ne.0) then
         read(nconfig,*) cell(1),cell(2),cell(3)
         read(nconfig,*) cell(4),cell(5),cell(6)
         read(nconfig,*) cell(7),cell(8),cell(9)
      endif

c zero variables
      xmin = 1.0d12
      xmax = -1.0d12
      ymin = 1.0d12
      ymax = -1.0d12
      zmin = 1.0d12
      zmax = -1.0d12

      icount = 1
c read in protein coordinates in DL_POLY format  100  read(nconfig,'(2a4,i10,6x,a4,x,a3,i6)',end=140) atmnam1(icount)
     +     ,atmnam0(icount),ic,pdbatm(icount),resnam(icount)
     +     ,resnum(icount)
      read(nconfig,*,end=140) px(icount),py(icount),pz(icount) c test for max and min
      if (px(icount).lt.xmin) xmin = px(icount)
      if (px(icount).gt.xmax) xmax = px(icount)
      if (py(icount).lt.ymin) ymin = py(icount)
      if (py(icount).gt.ymax) ymax = py(icount)
      if (pz(icount).lt.zmin) zmin = pz(icount)
      if (pz(icount).gt.zmax) zmax = pz(icount)

      if (levcfg.eq.1) then
         read(nconfig,'(a40)',end=140) cjunk
      else if (levcfg.eq.2) then
         read(nconfig,'(a40)',end=140) cjunk
         read(nconfig,'(a40)',end=140) cjunk
      endif
      icount = icount + 1
      goto 100
 140  natms = icount - 1
      print*,'limit distance for JONG =',djong
      print*,'number of atoms in GABA_A receptor =',natms
      write(6,'(''    xmin, xmax = '',2f10.5)') xmin,xmax
      write(6,'(''    ymin, ymax = '',2f10.5)') ymin,ymax
      write(6,'(''    zmin, zmax = '',2f10.5)') zmin,zmax
c set trap
      if (natms.gt.patms) then
         write(0,*)'patms too small, set it to ',natms
         write(0,*)'Abort execution.'
         stop
      endif

c read in CONFIG.mb4a file header
      read(nconf2,'(a40)') title2
      read(nconf2,*)levcfg2,imcon2
      if (imcon2.ne.0) then
         read(nconf2,*) cell2(1),cell2(2),cell2(3)
         read(nconf2,*) cell2(4),cell2(5),cell2(6)
         read(nconf2,*) cell2(7),cell2(8),cell2(9)
      endif
      icount = 1
c read in original membrane atom names in DL_POLY format  200  read(nconf2,'(a8,i10,6x,a4,x,a3,i6)',end=220) atmnam2(icount),ic
     +     ,mematm(icount),memnam(icount),memnum(icount)
      read(nconf2,*,end=220) xxx(icount),yyy(icount),zzz(icount)
      if (levcfg2.eq.1) then
         read(nconf2,'(a40)',end=220) cjunk
      else if (levcfg2.eq.2) then
         read(nconf2,'(a40)',end=220) cjunk
         read(nconf2,'(a40)',end=220) cjunk
      endif
      icount = icount + 1
      goto 200
 220  natms2 = icount - 1
      print*,'number of atoms in hydrated DMPC membrane =',natms2

c test for clashes
      mark = 0
      do 300 i = 1,natms
         do 340 j = 1,natms2
            dxx(1) = px(i) - xxx(j)
            dyy(1) = py(i) - yyy(j)
            dzz(1) = pz(i) - zzz(j)
            call images(imcon,0,1,1,cell,dxx,dyy,dzz)
            if ((dxx(1).lt.djong).and.(dyy(1).lt.djong)
     +           .and.(dzz(1).lt.djong)) then
               dist2 = dxx(1)**2 + dyy(1)**2 + dzz(1)**2
               if (dist2.lt.djong2) then c stick the first two sets of clashed data in
                  if (mark.lt.2) then
                     mark = mark + 1
c clash is found; save the clashed atomic number into a list
                     jonglist(mark) = j
c and then save the clashed resiude number into a different list
                     jonglist2(mark) = memnum(j)
                  else
c go through the current list; if this atom is already listed, then there is c no need to re-list it
                     do 360 k = 1,mark
                        if (jonglist(k).eq.j) then
                           goto 320
                        else
                           continue
                        endif
 360                 continue
c if atom not on list already, then go on with distance testing c set trap for size of jonglist array
                     if (mark.gt.mxjong) then
                        write(0,*)'mark > mxjong; abort.'
                        stop
                     endif
c if the list array is large enough, then place it in storage
                     mark = mark + 1
c clash is found; save the clashed atomic number into a list
                     jonglist(mark) = j
c and then save the clashed resiude number into a different list
                     jonglist2(mark) = memnum(j) c once entry is saved, get out of all this mess
                     goto 320
                  endif
               endif
 320           continue
            endif
 340     continue
 300  continue

c determine cell size
      do 380 i= 1,9
         write(0,*)cell(i),cell2(i)
         if (cell(i).lt.cell2(i)) cell(i) = cell2(i)  380  continue

c print header
      levcfg = 0
      write(nconfout,'(a40)') title
      write(nconfout,'(2i10)')levcfg,imcon2
      write(nconfout,'(3f20.12)') cell2(1),cell2(2),cell2(3)
      write(nconfout,'(3f20.12)') cell2(4),cell2(5),cell2(6)
      write(nconfout,'(3f20.12)') cell2(7),cell2(8),cell2(9) c printout only those atoms that do not clash, firstly ALL protein atoms
      do 400 i = 1,natms
         write(nconfout,'(2a4,i10,6x,a4,x,a3,i6)') atmnam1(i)
     +        ,atmnam0(i),i,pdbatm(i),resnam(i),resnum(i)
         write(nconfout,'(2(f16.10,4x),f16.10)')px(i),py(i),pz(i)
 400  continue

c mark out the clashed atoms
      do 440 i = 1,natms2
         index(i) = 0
 440  continue
      l = 0
      do 460 i = 1,natms2
         do 480 j = 1,mark
            if (i.eq.jonglist(j)) then
               index(i) = -1
               l = l + 1
            endif
 480     continue
 460  continue

c determine the highest memnum in the system
      maxmemnum = memnum(natms2)
      do 500 j = 1,maxmemnum
         natpmol(j) = 0
         indexsum(j) = 0
 500  continue
c logic of this programme: (1) if a molecule is DMPC, it is removed if and c only if at least half of its atoms clash with the protein (2) if a molecule c is water, it is removed even if one of its atoms clashes with the protein c go through jong list to pull out all the molecules which clash with the c protein
      do 520 i = 1, natms2
c go through all molecules
         do 540 j = 1,maxmemnum
            if (memnum(i).eq.j) then
               natpmol(j) = natpmol(j) + 1
               if (index(i).eq.-1) indexsum(j) = indexsum(j) - 1
            endif
 540     continue
 520  continue

      do 560 j = 1,maxmemnum
         status(j) = dble(indexsum(j))/dble(natpmol(j))
 560  continue

c re-define index(i)
      do 580 i = 1,natms2
         index(i) = 0
 580  continue
      do 600 i = 1,natms2
         do 620 j = 1,maxmemnum
c this rejfrac business only applies to the phospholipids
            if (memnam(i).ne.'HOH') then
               if ((memnum(i).eq.j).and.(status(j).lt.rejfrac)) then
                  index(i) = -1
               endif
c if the molecule is water, then even if there is one atom which clashes with c the protein, the water molecule is rejected
            else
               if ((memnum(i).eq.j).and.(status(j).lt.-0.1)) then
                  index(i) = -1
               endif
            endif
 620     continue
 600  continue

c write out results
      ix = natms
      do 700 i = 1,natms2
         if (index(i).ne.-1) then
            ix = ix + 1
            write(nconfout,'(a8,i10,6x,a4,x,a3,i6,f6.2)') atmnam2(i),ix
     +           ,mematm(i),memnam(i),memnum(i),status(memnum(i))
            write(nconfout,'(2(f16.10,4x),f16.10)')xxx(i),yyy(i),zzz(i)
         endif
 700  continue

      stop
      end
