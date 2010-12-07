      program pbd2edt

c***********************************************************************
c     
c     dl_poly routine to read a broohaven pdb file
c     (or the AMBER processed equivalent) and to extract ...
c     1)  atom  names and types
c     2)  residue name and numbers
c     3)  x,y,z coordinates
c     4)  site charges
c     and then print in an edit.out compatible form
c
c     author t.forester feb 1996
c     copyright daresbury laboratory 1996
c
c     itt
c     2010-10-30 17:20:53
c     1.3
c     Exp
c
c***********************************************************************
      
      parameter(mxres=50, mxsit=30)
      implicit real*8(a-h,o-z)
      
      character*8 header
      character*80 oneline,pdbfile
      character*4 name,ambnam,oxt
      character*3 res,rold,his,bac
      character*2 strand,strand1
      logical yes, dna,rna,charge,hsafe,safe
      character*4 resnam(mxres),sitnam(mxsit,mxres),ffnam(mxsit,mxres)
      dimension qsit(mxsit,mxres),rept(mxsit,mxres),count(mxsit)
      dimension nsite(mxres),join(mxsit,mxres),index(mxsit)
      dimension cell(9)
c
c     flag for all hydrogens present

      hsafe=.true.
      safe=.true.

      write(6,'(/,a)') 'name of pdb file'
      read(*,'(a80)') pdbfile
      open(10,file=pdbfile,status='old')

      oneline ='EDTOUT'
      open(11,file=oneline(1:6))
      write(*,'(/,2a)') 'enter header record for ',oneline(1:6)
      read(*,'(a80)') oneline
      write(11,'(a)') oneline
      strand1 ='%%'

 199  yes = .true.
      write(*,'(/,a)') 'do you have RNA units ? (y/n)'
      read(*,'(a80)') oneline
      call lowcase(oneline,80)
      call strip(oneline,80)
      if(oneline(1:1).eq.'n') then
        dna=.true.
        bac = 'dna'
        rna=.false.
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
      
 399  write(*,'(/,2a)')'do you want atom charges read from ',
     x     pdbfile(1:20)
      read(*,'(a80)') oneline
      call lowcase(oneline,80)
      call strip(oneline,80)
      if(oneline(1:1).eq.'y') then
        charge=.true.
      elseif (oneline(1:1).eq.'n') then
        charge=.false.
        write(*,'(/,a)') ' charges will be read from residue map file'
      else
        write(6,*) ' you must answer "y" or "n"'
        write(6,*) 
        yes = .false.
      endif
      
      if(.not.yes) goto 399
c     
c     read charge and name database
      
      nres=0
      call data_read
     x  (mxres,mxsit,nres,join,qoxt,oxt,header,
     x   resnam,sitnam,ffnam,qsit,rept,nsite)
c
c     replacement name for HIS residues

      if(header.eq.'amber   ') then
        his = 'HID'
      elseif(header.eq.'gromos  ') then
        his = 'HIA'
      endif

c
c     modify dna/rna residues in database to include backbone

      if(dna.or.rna) call res_back
     x   (bac,mxres,mxsit,nres,join,
     x   resnam,sitnam,ffnam,qsit,rept,nsite)
      
c     
c     process pdb file to make into AMBER style edit.out      
      blank = -9999.99
      irs = 0
      lres = 0
      i = 0
      do imega= 1,1000000
        read(10,'(a80)',end=100,err=100) oneline
        if(oneline(1:6).eq.'CRYST1') then
          backspace(10)
c
c     get cell vectors
          read(10,'(6x,3f9.0,3f7.0)') a1,b1,c1,alp,bet,gam
          gam = gam*(3.141592653589793d0/180.d0)
          bet = bet*(3.141592653589793d0/180.d0)
          alp = alp*(3.141592653589793d0/180.d0)
          cell(1) = a1
          cell(2) = 0.d0
          cell(3) = 0.d0
          cell(4) = b1*cos(gam)
          cell(5) = b1*sin(gam)
          cell(6) = 0.d0
          cell(7) = c1*cos(bet)
          cell(8) = (b1*c1*cos(alp) - cell(4)*cell(7))/cell(5)
          cell(9) = sqrt(c1*c1 - cell(8)**2-cell(7)**2)
          
          write(11,'(/,a,/,3(3f20.10,/))')'CELL_VECTORS',cell
          
        elseif((oneline(1:4).eq.'ATOM').or.
     x     (oneline(1:6).eq.'HETATM')) then
          backspace(10)
          i = i+1
        read(10,'(6x,i5,2x,a4,a3,a2,i4,4x,3f8.3,f7.3)',end=100,err=100)
     x       idum,name,res,strand,ires,x,y,z,q 

        if(strand.ne.strand1) then
          write(11,'(/,6x,a,/)') 'MOLECULE'
          strand1=strand
        endif
            
        if(irs.ne.ires) then
c     
c     a new residue has been found
c     
c     check all atoms found in previous residue
          
          safe=.true.
          
          if(irs.gt.0.and.rold.ne.'mtl') then

            call  res_check
     x         (safe,rold,sitnam,irs,lres,count,rept,nsite,mxres,mxsit)
          endif

          hsafe = (safe.and.hsafe)
          
          irs=ires
          if(res.eq.'HIS'.or.res.eq.'his') res = his
          rold =res
          ll =0
          call strip(rold,3)
          call lowcase(rold,3)
          write(11,*)
          write(11,'(a,i4,a,a3)')'     RESIDUE',ires,' =  ',res
          write(11,*)
          write(11,'(a,i6)')'     BOND ARRAY BEGINS WITH ',i
          write(11,*) 
c     
c     find residue in database
          
          call lowcase(res,3)
          yes=.false.
          do jres = 1,nres
            if(resnam(jres)(1:3).eq.res) then
              lres = jres
              yes=.true.
            endif
          enddo
          if(.not.yes) then
            write(*,*) 'residue ',res,' not in database '
            call exit(0)
          endif
          do i9 = 1,mxsit
            count(i9) = 0.d0
            index(i9) = 0
          enddo
        endif

        call strip(res,3)
        call lowcase(res,3)
        call strip(name,4)
c
c     site counter 
        ll = ll+1
c     
c     assign ff names and check charges
        
        q1 = blank
        j = -99

        call data_find
     x       (lres,jj,name,res,sitnam,count,rept,mxres,
     x       mxsit,nsite)
        if(jj.gt.0) index(jj)=i
c
c     find connection ... 

        if(ll.gt.1.and.jj.gt.0) then
          kk = join(jj,lres)
          if(kk.gt.0) then
            if(count(kk).gt.0.) j = index(kk)
          endif
        endif
        
        if(.not.charge) then
          if(jj.gt.0) q1 = qsit(jj,lres)
        else
          q1 = q
        endif
        if(jj.gt.0) ambnam = ffnam(jj,lres)
        
c     
c     apply patch for terminating oxygen
        
        if(name.eq.'OXT ') then 
            q1 = qoxt
            ambnam=oxt
        endif
        
        if(j.eq.0) j = -99
        
        if(q1.eq.blank) then 
          write(6,*) 'could not find atom ',name,' ', i, 
     x         ' residue ',res,' in database'
          call exit(0)
          
        elseif(q1.ne.q) then
          if(abs(q1-q).le.0.0006) then
            q = q1
          else
            if(charge)
     x           write(6,'(a,i7,a,a4,a,a4,a,f8.4,a,f8.4)')
     x           'warning: non-standard charge on atom ',i,
     x           ' ',name,'(',res,') : ',q, ' vs ',q1
          endif
        endif
c     
c     write out edit out file
     
        if(header.eq.'amber   ') then
          if(dna.or.rna) call name_fix(name,4)
        endif
        
        write(11,'(2i5,3x,a4,2x,a4,6x,3f10.4,f10.6)')
     x       i,j,name,ambnam,x,y,z,q1
        
      endif
      enddo
      
 100  close(10)
      close(11)

      write(*,*)
      write(*,*) 'To process your output file (EDTOUT) into '
      write(*,*) 'DL_POLY FIELD and CONFIG files '
      write(*,*) 'run the utility "ffgen". If you have missing'
      write(*,*) 'atoms you will first need to run the utility'
      write(*,*) '"filledt".   In any event '
      write(*,*) 'you will also need to insert into EDTOUT the lines'
      write(*,'(/a,/)')'      MOLECULE    (with 6 leading spaces)'
      write(*,*) 'between entries for residues of different molecules'
      write(*,*) '(with only one such record for the solvent residues)'
      write(*,*)
      if(hsafe) then
        write(*,*)
        write(*,*) 'missing atoms have been detected'
        write(*,*) 'run the utility "molfill" on the edited EDTOUT'
        write(*,*) 'file **before** you run "ffgen"'
        write(*,*)
      endif

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

      character*1 string(*)
      character*1 letter

      do i = 1,min(255,length)

        letter = string(i)

        if(letter.eq.'A') then
          letter = 'a'
        else if(letter.eq.'B') then
          letter = 'b'
        else if(letter.eq.'C') then
          letter = 'c'
        else if(letter.eq.'D') then
          letter = 'd'
        else if(letter.eq.'E') then
          letter = 'e'
        else if(letter.eq.'F') then
          letter = 'f'
        else if(letter.eq.'G') then
          letter = 'g'
        else if(letter.eq.'H') then
          letter = 'h'
        else if(letter.eq.'I') then
          letter = 'i'
        else if(letter.eq.'J') then
          letter = 'j'
        else if(letter.eq.'K') then
          letter = 'k'
        else if(letter.eq.'L') then
          letter = 'l'
        else if(letter.eq.'M') then
          letter = 'm'
        else if(letter.eq.'N') then
          letter = 'n'
        else if(letter.eq.'O') then
          letter = 'o'
        else if(letter.eq.'P') then
          letter = 'p'
        else if(letter.eq.'Q') then
          letter = 'q'
        else if(letter.eq.'R') then
          letter = 'r'
        else if(letter.eq.'S') then
          letter = 's'
        else if(letter.eq.'T') then
          letter = 't'
        else if(letter.eq.'U') then
          letter = 'u'
        else if(letter.eq.'V') then
          letter = 'v'
        else if(letter.eq.'W') then
          letter = 'w'
        else if(letter.eq.'X') then
          letter = 'x'
        else if(letter.eq.'Y') then
          letter = 'y'
        else if(letter.eq.'Z') then
          letter = 'z'
        endif

        string(i) = letter

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

      character*(*) string
      
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

      subroutine data_read
     x  (mxres,mxsit,nres,join,qoxt,oxt,header,
     x   resnam,sitnam,ffnam,qsit,rept,nsite)

      implicit real*8(a-h,o-z)
      character*40 oneline
      character*4 oxt
      character*8 header
      character*4 resnam(mxres),sitnam(mxsit,mxres),ffnam(mxsit,mxres)
      dimension qsit(mxsit,mxres),rept(mxsit,mxres)
      dimension nsite(mxres),join(mxsit,mxres)
      logical safe
c
c     open data file

      write(*,*) 'enter name of residue map file (e.g. respdb.amb ',
     x     'or respdb_u.amb)'
      read(*,'(a40)') oneline
      open(20,file=oneline,status='old')
      read(20,'(a8)') header
      call strip(header,8)
      call lowcase(header,8)
      safe = .false.
      if(header.eq.'gromos  ') safe = .true.
      if(header.eq.'amber   ') safe = .true.
      if(.not.safe) then
        write(*,*) 'unrecognised forcefield type in ',pdbfile
        call exit(0)
      endif
c
c     qoxt is charge on terminal oxylate group
c     to be found from force field

      qoxt = -9.d50

c     header
      read(20,*) 
c     number of residues
      read(20,*) nres  
      
      do lres = 1,nres
        read(20,'(10x,a4)') resnam(lres)
        call strip (resnam(lres),4)
        call lowcase(resnam(lres),4)

        read(20,'(7x,i10)') nsite(lres)
        do i = 1,nsite(lres)
          read(20,'(2i4,2a4,2f8.0)') idummy,
     x         join(i,lres),sitnam(i,lres),ffnam(i,lres),
     x         qsit(i,lres),rept(i,lres)
          if(join(i,lres).eq.0) join(i,lres)=-99
          rept(i,lres) = max(rept(i,lres),1.d0)
          call strip(sitnam(i,lres),4)
          call strip(ffnam(i,lres),4)
c
c     pick up qoxt charge ...
          if(resnam(lres).eq.'asp') then
            if(sitnam(i,lres)(1:1).eq.'O') then
              if(sitnam(i,lres)(2:2).ne.' ') then
                qoxt = qsit(i,lres)
                oxt = ffnam(i,lres)
              endif
            endif
          endif
        enddo

      enddo

      write(*,'(/,i3,2a)') nres,
     x     ' units read in from database ',oneline
      return
      end
      
      subroutine data_find
     x     (lres,jj,name,res,sitnam,count,rept,mxres,
     x     mxsit,nsite)
      
      implicit real*8(a-h,o-z)
      character*4 name,sitnam(mxsit,mxres),res
      dimension rept(mxsit,mxres),count(mxsit)
      dimension nsite(mxres)
      
c     
c     look at as many characters as necessry to uniquely identify
c     the site
      
      k = 1
      
 10   look = 0
      jj= 0
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

      subroutine res_check
     x  (safe,res,sitnam,ires,lres,count,rept,nsite,mxres,mxsit)

      implicit real*8(a-h,o-z)
      character*3 res
      character*4 sitnam(mxsit,mxres)
      dimension count(mxsit),rept(mxsit,mxres)
      dimension nsite(mxres)
      logical safe

      safe=.true.
      
      do i = 1,nsite(lres)

        if(abs(count(i)-rept(i,lres)).gt.0.5) then
          safe=.false.
        endif
      enddo

      if(.not.safe) then

        do i = 1,nsite(lres)

          j = nint(count(i)-rept(i,lres))
          
          if(j.lt.0) then
            if((sitnam(i,lres)(1:1).ne.'H').and.
     x           (sitnam(i,lres)(1:2).ne.'LP'))
     x        write(*,'(4x,a,i3,3a,i3,3a,i6)')'missing ',-j,
     x        ' atoms of type ',sitnam(i,lres),'(',i,') from ',
     x        ' residue ',res,ires

          elseif(j.gt.0) then
            write(*,'(4x,i3,3a,i3,a)')j,' too many  atoms of type ',
     x           sitnam(i,lres),'(',i,')'

          endif

        enddo

      endif
      return
      end
      subroutine res_back
     x   (bac,mxres,mxsit,nres,join,
     x   resnam,sitnam,ffnam,qsit,rept,nsite)

c     subroutine to append 'backbone' data into 'residue' data
c     
      implicit real*8(a-h,o-z)

      logical yes
      character*3 bac
      character*4 resnam(mxres),sitnam(mxsit,mxres),ffnam(mxsit,mxres)
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

      subroutine name_fix(name,kk)
      character*1 name(kk)
c
c     replaces * (pdb convention) in backbone name with ' : amber style
      do i = 1,kk

        if(name(i).eq.'*') name(i) = "'"
        
      enddo
      return
      end
