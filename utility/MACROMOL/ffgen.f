      program ffgen

c***********************************************************
c
c     program to extract force field for a macromolecule
c     and interface it to DL_POLY.
c     input file: EDTOUT (AMBER 3, 4 edit.out style)
c               : FF.dat  (title plus additional info)
c     output    : CONFIG (coordinates etc)
c               : FIELD  (force-field)
c
c     copyright -  daresbury laboratory 1993
c
c     author    -  t forester    April 1993
c     modified  - t.forester       feb 1996
c
c**********************************************************

      implicit real*8(a-h,o-z)

      parameter (ntyp = 100, mxatm = 30000,
     x     mgb =250, mga = 500,
     x     mgi = 200, mgd = 500, mlj =100,mhb = 80,
     x     mxbnd = 10000, mxtmls=20)
      parameter(mxang=mxbnd*3, mxdih=mxang*3)

c     mxatm = max number of atoms in molecule
c     ntyp = number of atom  types in forcefield
c     mgb  = number of bond  types in forcefield
c     mga  = number of bond angles  types in forcefield
c     mgi  = number of improper dihedrals  types in forcefield
c     mgd  = number of dihedrals  types in forcefield
c     mhb  = number of hygrogen bond types in forcefield
c     mxbnd = max number of bonds in macromolecule
c     mxang = max number of angles in macromolecule
c     mxtmls = max number of molecular types

      logical found,hcons,allcons,ldo
      character*8 rname(mxatm),blank8,ffnam
      character*4 ambnam(mxatm),name(ntyp),water,hwater
      character*4 bndtyp,angtyp,dihtyp,diityp,vdtyp,hbtyp
      character*3 resn(mxatm),res
      character*8 nami,namj
      character*80 title,oneline,unitlin

      dimension molatm(mxtmls)

      dimension cell(9)
      dimension qq(mxatm)
      dimension ityp(mxatm),amass(ntyp)
      dimension ires(mxatm)
      dimension x(mxatm),y(mxatm),z(mxatm)
      dimension nbds(mxatm),nexatm(mxatm,10),iamb(mxatm)
      dimension l14a(mxdih),l14b(mxdih),l14c(mxdih),l14d(mxdih)

      character*4 bnda(mgb),bndb(mgb)
      dimension   cb(mgb),bo(mgb)
      dimension   ibda(mxbnd),ibdb(mxbnd)

      character*4 dia(mgd),dib(mgd),dic(mgd),did(mgd)
      dimension   cq(mgd),qd(mgd),td(mgd)
      dimension   np(mgd)

      character*4 disa(mgd),disb(mgd),disc(mgd),disd(mgd)
      dimension   cqs(mgd),qds(mgd),tds(mgd)
      dimension   nps(mgd)

      character*4 anga(mga),angb(mga),angc(mga)
      dimension   ct(mga),to(mga)
      dimension langa(mxang),langb(mxang),langc(mxang)
      dimension ldihb(mxang),ldihc(mxang)

      character*4 ima(mgi),imb(mgi),imc(mgi),imd(mgi)
      dimension   aimp(mgi),bimp(mgi),cimp(mgi)
      character*4 imsa(mgi),imsb(mgi),imsc(mgi),imsd(mgi)
      dimension   aimps(mgi),bimps(mgi),cimps(mgi)

      character*4 hba(mhb),hbb(mhb)
      dimension h12(mhb),h10(mhb)

      parameter(ntyp2 = (ntyp*(ntyp+1))/2)
      character*4 nlj(ntyp),nlja(ntyp2),nljb(ntyp2),nljh(ntyp)
      dimension   c12(ntyp),c6(ntyp)
      dimension   c12b(ntyp),c6b(ntyp)
      dimension   c12a(ntyp2),c6a(ntyp2)
      dimension   ljeqiv(ntyp),ijchk(ntyp)


      water = "OW  "
      hwater = "HW  "
c****************************************************************************
c
c     initialize atom names as blanks

      blank8 = "        "

      do i = 1,mxatm

        rname(i) = blank8

      enddo
c
c     open output file: 7 = FIELD file, 8 = CONFIG file,
c     20 = WARNING file

      ifield = 7
      conf = 8
      ier = 20

      open(ifield,file='FIELD')
      open(conf,file='CONFIG')
      open(ier,file='WARNING.ff')

c
c     define input streams

      ird1 = 10
      ird2 = 11
      ird3 = 12

c*****************************************************************************
c
c     read forcefield data
c
c***************************************************************************

      write(*,*) 'Enter filename with forcefield parameter database'
      read(*,'(a80)') oneline
      open(ird2,file=oneline,status='old')

      read(ird2,'(a8)') ffnam
      call strip(ffnam,8)
      call lowcase(ffnam,8)
      found = .false.
      if(ffnam.eq.'gromos  ') found = .true.
      if(ffnam.eq.'amber   ') found = .true.
      if(.not.found) then
        write(*,'(/,2a)') ' error - unkown force field ',ffnam
        call exit(0)
      endif

      call block(inum,ird2)
      natm = inum

      do i = 1,natm

        read(ird2,'(a4,f12.0)',err=997) name(i),amass(i)

      enddo

c
c     read bond stretch parameters

 997  call block(inum,ird2)

      nbd = inum
      if(nbd.gt.0) then
        read(ird2,*,err=996) bfac
        read(ird2,'(a4)') bndtyp

        if(bndtyp.ne.'harm') then
          write(*,*) ' unrecognized bond potential form ',bndtyp
          call exit(0)
        endif
      endif
      do i = 1,nbd
        read(ird2,'(2a4,2f12.0)',err=996) bnda(i),bndb(i),cb(i),bo(i)
        call strip(bnda(i),4)
        call strip(bndb(i),4)
      enddo

c
c     read angle bend parameters

 996  call block(inum,ird2)

      nang = inum
      if(nang.gt.0) then
        read(ird2,*,err=996) afac
        read(ird2,'(a4)') angtyp

        if(angtyp.ne.'harm') then
          write(*,*) ' unrecognized angular potential form ',angtyp
          call exit(0)
        endif
      endif
      do i = 1,nang

        read(ird2,'(3a4,2f12.0)',err=995) anga(i),angb(i),angc(i)
     x       ,ct(i),to(i)
        call strip(anga(i),4)
        call strip(angb(i),4)
        call strip(angc(i),4)
      enddo
c
c     read general proper torsion parameters

 995  call block(inum,ird2)

      ndi = inum
      if(ndi.gt.0) then
        read(ird2,*,err=996) dfac
        read(ird2,'(a4)') dihtyp

        if((dihtyp.ne.'cos ').and.(dihtyp.ne.'harm')) then
          write(*,*) ' unrecognized dihedral potential form ',dihtyp
          call exit(0)
        endif

        read(ird2,*) delec
        read(ird2,*) dvdw
      endif
      do i = 1,ndi

        read(ird2,'(4a4,i5,3f12.0)',err=994)
     x       dia(i),dib(i),dic(i),did(i),
     x       np(i),cq(i),qd(i),td(i)

        call strip(dia(i),4)
        call strip(dib(i),4)
        call strip(dic(i),4)
        call strip(did(i),4)

      enddo

c
c     read specific proper torsion parameters

 994  call block(inum,ird2)

      ndis = inum
      do i = 1,ndis

        read(ird2,'(4a4,i5,3f12.0)',err=993)
     x       disa(i),disb(i),disc(i),disd(i),
     x       nps(i),cqs(i),qds(i),tds(i)

        nps(i) = max(nps(i),1)
        call strip(disa(i),4)
        call strip(disb(i),4)
        call strip(disc(i),4)
        call strip(disd(i),4)
      enddo

c
c     read general improper torsion parameters

 993  call block(inum,ird2)

      nim = inum
      if(nim.gt.0) then
        read(ird2,*,err=996) difac
        read(ird2,'(a4)') diityp

        if((diityp.ne.'cos ').and.(diityp.ne.'harm')) then
          write(*,*) ' unrecognized dihedral potential form ',diityp
          call exit(0)
        endif

        read(ird2,*) dielec
        read(ird2,*) divdw
      endif
      do i = 1,nim

        read(ird2,'(4a4,3f12.0)',err=893) ima(i),imb(i),imc(i),
     x       imd(i),aimp(i),bimp(i),cimp(i)

        call strip(ima(i),4)
        call strip(imb(i),4)
        call strip(imc(i),4)
        call strip(imd(i),4)

      enddo

c
c     read specific improper torsion parameters

 893  call block(inum,ird2)

      nims = inum

      do i = 1,nims

        read(ird2,'(4a4,3f12.0)',err=793) imsa(i),imsb(i),imsc(i),
     x       imsd(i),aimps(i),bimps(i),cimps(i)

        call strip(imsa(i),4)
        call strip(imsb(i),4)
        call strip(imsc(i),4)
        call strip(imsd(i),4)

      enddo

c
c     read Hydrogen bond parameters

 793  call block(inum,ird2)

      ihb = inum

      if(ihb.gt.0) then
        read(ird2,*,err=992) hbfac
        read(ird2,'(a4)') hbtyp

        if(hbtyp.ne.'hbnd') then
          write(*,*) ' unrecognized h-bond potential form ',hbtyp
          call exit(0)
        endif

      endif
      do i = 1,ihb

        read(ird2,'(2a4,2f12.0)',err=992) hba(i),hbb(i),h12(i),h10(i)
        call strip(hba(i),4)
        call strip(hbb(i),4)
      enddo

 992  continue

c
c     Lennard Jones equivalents

      do i = 1,natm

        ljeqiv(i) = i

      enddo

      call block(inum,ird2)

      ieqiv = inum

      do i = 1, ieqiv

        inam = 0

        read(ird2,'(a4,a76)')
     x       nami,oneline
        call strip(nami,4)

        do j = 1,73,4

          namj = oneline(j:j+3)
          call strip(namj,4)
c
c     identify types

          jnam = 0

          do k = 1,natm

            if(namj.eq.name(k)) jnam = k
            if(nami.eq.name(k)) inam = k

          enddo

          if(inam.eq.0) call error(8,i)
          if(jnam.ne.0) ljeqiv(jnam) = inam

        enddo

      enddo

c
c     read non-bonded parameters

      call block(inum,ird2)

      ilj = inum
      if(ilj.gt.0) then
        read(ird2,*,err=996) vdfac
        read(ird2,'(a4)') vdtyp

        if((vdtyp.ne.'ljmn').and.(vdtyp.ne.'12-6').and.
     x    (vdtyp.ne.'lj  ')) then
          write(*,*) ' unrecognized VdW potential form ',vdtyp
          call exit(0)
        endif
      endif
      do i = 1,ilj
        read(ird2,'(a4,2f12.0)',err=991) nlj(i),c12(i),c6(i)
        call strip(nlj(i),4)
      enddo

 991  continue
c
c     convert to equivalent LJ c12 and c6 values

      if(vdtyp.eq.'ljmn')then

c     c12 = rmin --> sig,     c6 = eps

        rmin = 1.d0/(2.d0**(1.d0/6.d0))

        do i = 1,ilj
          c12(i) = c12(i)*rmin
        enddo

      endif
c
c     read specific vdw cross term parameters

      call block(inum,ird2)

      ljx = inum

      do i = 1,ljx
        read(ird2,'(2a4,2f12.0)',err=991) nlja(i),nljb(i),c12a(i),c6a(i)
        call strip(nlja(i),4)
        call strip(nljb(i),4)
        if(vdtyp.eq.'ljmn') c12a(i)=c12a(i)*rmin
      enddo
c
c     read 1-4 Vdw parameters

      call block(inum,ird2)

      lj14 = inum

      do i = 1,lj14
        read(ird2,'(a4,2f12.0)',err=991) nljh(i),c12b(i),c6b(i)
        call strip(nljh(i),4)
        if(vdtyp.eq.'ljmn') c12b(i)=c12b(i)*rmin
      enddo

      if(vdtyp.eq.'ljmn') vdtyp='lj  '

c
c*****************************************************************************
c     end of forcefield parameters
c
c     start of MACROMOLECULE file
c
c****************************************************************************
c
c     open protein data file

      open(ird3,file='FF.dat',status='old')
c
c     read sequence

      read(ird3,'(a80)')  title

c
c     check for which bonds should be constrained

      hcons = .false.
      allcons = .false.

      read(ird3,'(a80)',end=145) oneline

      call strip(oneline,80)
      call lowcase(oneline,80)

      if (oneline(1:9).eq."constrain") then

        oneline(1:71) = oneline(10:80)
        call strip(oneline,71)

        if (oneline(1:1).eq.'h') hcons = .true.
        if (oneline(1:3).eq.'all') allcons = .true.

      else

        backspace(ird3)

      endif
  145 continue
c
c     check energy units

      eunit =1.d0
      unitlin ='UNITS kcal                             '
      oneline='        '
      read(ird3,'(a80)',end=345) oneline

      call strip(oneline,80)
      call lowcase(oneline,80)

      if(oneline(1:5).eq.'units') then

        oneline(1:75) = oneline(6:80)
        call strip(oneline,75)

        if(oneline(1:2).eq.'ev') then

          unitlin = 'UNITS eV                                 '
          eunit = 23.060542d0

        elseif (oneline(1:2).eq.'kj') then

          unitlin ='UNITS kJ                                  '
          eunit = 4.184d0

        elseif (oneline(1:4).eq.'kcal') then

          unitlin ='UNITS kcal                                 '
          eunit = 1.d0

        elseif (oneline(1:8).eq.'internal') then

          unitlin ='UNITS internal                            '
          eunit = 4.184d2

        else

          unitlin ='UNITS kcal                                 '
          eunit = 1.d0

        endif

      else

        backspace(ird3)
        unitlin ='UNITS kcal                                 '
        eunit = 1.d0

      endif

  345 continue
c
c     convert all force field tables to new units
c     taking into account prefactors in potential forms

      call units(eunit*bfac,cb,nbd)
      call units(eunit*afac,ct,nang)
      call units(eunit*dfac,cq,ndi)
      call units(eunit*dfac,cqs,ndis)
      call units(eunit*difac,aimp,nim)
      call units(eunit*difac,aimps,nims)
      call units(eunit*hbfac,h12,ihb)
      call units(eunit*hbfac,h10,ihb)
      call units(eunit*vdfac,c6a,ljx)
      if(vdtyp.eq.'12-6') then
        call units(eunit*vdfac,c12a,ljx)
        call units(sqrt(eunit*vdfac),c12,ilj)
        call units(sqrt(eunit*vdfac),c6,ilj)
        call units(sqrt(eunit*vdfac),c12b,lj14)
        call units(sqrt(eunit*vdfac),c6b,lj14)
      else
        call units(eunit*vdfac,c6,ilj)
        call units(eunit*vdfac,c6b,lj14)
      endif

      write(6,*) 'reading EDTOUT ..... '
c****************************************************************************
c
c     read connectivity information - data produced from 'edit'
c     run under AMBER or using DLPOLY utility 'pdb2edt'
c
c****************************************************************************

c
c     number of different species in simulation

      imol = 0
      iatm = 0

      open(ird1,file='EDTOUT',status='old')

c
c     Process bond entries -- harmonic and constraints

      ibnd = 0
      icon = 0
      kres = 0

 202  call bndtitle(ird1,ierr,res)

      if (ierr.eq.1) goto 201
c
c     new molecule
      if(ierr.eq.-999) then

        imol = imol+1
        kres =0
        if(imol.gt.mxtmls) call error(10,imol)

        if(imol.gt.1) then
          molatm(imol-1) = iatm
        endif
      endif

      kres = kres+1
      call bnddata(ird1,ifield,iatm,ibnd,kres,mgb,cb,bo,
     x     bnda,bndb,rname,ambnam,ibda,ibdb,ires,x,y,z,qq,
     x     res,resn)

      goto 202

 201  continue

      if(imol.gt.0) molatm(imol) = iatm
      if(ibnd.gt.mxbnd) call error(4,ibnd)

c*********************************************************************
c     modified atomic positions and charges
c********************************************************************

 110  read(ird3,'(a8,i10)',end = 120) nami,i
      call strip(nami,8)
      rname(i) = nami
      ambnam(i) = nami(1:4)
      read(ird3,*) x(i),y(i),z(i),qq(i)
      iatm = max(iatm,i)
      goto 110
 120  continue

c
c*********************************************************
c     interlude - write out CONFIG now ambnam's are known
c     assign masses and charges
c*******************************************************

ctrf  write(6,*) ' assigning solvent names .....'
c
c     assign force field names to solvent species

c
c     check for water - assumed to be tips3p

      do i = 1,iatm
        if((rname(i).eq."O   ").and.
     x       (rname(i+1).eq."H1  ").and.
     x       (rname(i+2).eq."H2  ")) then

          ambnam(i) = water
          ambnam(i+1) = hwater
          ambnam(i+2) = hwater

        endif
      enddo

      do i = 1,iatm
c
c     identify atom type

        call ambtype(ambnam(i),name,natm,iamb(i))

      enddo

      write(conf,'(a80)') title
      call write1(conf,ird1,iatm,ljeqiv,iamb,x,y,z,ambnam,name,rname,
     x     cell,imcon)

c*****************************************************************
c     end of logical interlude
c***************************************************************
c
c     construct connectivity table

ctrf  write(6,*) ' constructing connectivity .....'

      do i = 1,mxatm

        nbds(i) = 0

        do j = 1,10

          nexatm(i,j)  = 0

        enddo

      enddo

      write(6,*) 'creating connectivity tables .....'
      do i = 1,ibnd
c
c     indices of atoms in ith bond

        ia = ibda(i)
        ib = ibdb(i)
c
c     counters for connectivity

        nbds(ia) = nbds(ia) + 1
        nbds(ib) = nbds(ib) + 1
c
c     connectivity table

        nexatm(ia,nbds(ia)) = ib
        nexatm(ib,nbds(ib)) = ia

      enddo

      write(6,*) 'checking connectivities ..... '

      write(6,*) imol,' Molecules found '
c
c     loop twice - joining up non - hydrogens first

      do iloop = 1,3

        iatm0 = 1
        do imols = 1,imol

          iatm1 = molatm(imols)

          call concheck(ambnam,iatm0,iatm1,ibnd,
     x         iloop,mxatm,x,y,z,ibda,ibdb,nbds,nexatm)

          iatm0 = iatm1+1

        enddo

      enddo

      write(6,*) 'writing the FIELD file .....'

c-----------------------------------------------------------------------
c     now start processing FIELD file
c-----------------------------------------------------------------------

      write(ifield,'(a80)') title

      write(ifield,'(a80)') unitlin
      write(ifield,*) ' neutral  groups'
      write(ifield,*)' molecular types ',imol
      iatm0 = 0
      do imols = 1,imol
        if(imols.gt.1)iatm0 = molatm(imols-1)
        iatm1 = molatm(imols)
        katm = iatm1 - iatm0

        write(ifield,*) 'Molecule ',imols
        write(*,'(/,a,i5,/)') 'Molecule ',imols
c
c     check for tip3p for water - most obvious repeated molecule

        nummols=1
        if(ambnam(iatm0+1).eq.water) then
          nummols = katm/3
          katm=3
        endif

        write(ifield,*) 'nummols ',nummols

        if(nummols.gt.1) iatm1=iatm0+katm

        write(ifield,*) 'atoms ', katm

        do ii = 1,katm
          i = ii + iatm0
          write(ifield,'(a8,2f12.4,4i10,3x,a4,2x,a4,a3)')
     x         name(ljeqiv(iamb(i))),amass(iamb(i)),qq(i),1,0,
     x         ires(i),ii,ambnam(i),rname(i),resn(i)

        enddo
c
c     determine bond  parameters and print out

          ibc = 0
          ibnd1 = 0
          ibnd0 = 1
          ih0 = 0
          found = .false.
          do i = 1,ibnd

            ia = ibda(i)
            if((ia.gt.iatm0).and.(ia.le.iatm1)) then
              ibnd0 = i
              found = .true.
            endif

            if(found) goto 310

          enddo

 310      if(.not.found) then
            write(*,*) " warning - no bonds found"
          endif

            ibnd1 = ibnd
            found = .false.
            do i = ibnd0,ibnd

              ia = ibda(i)
              if(ia.gt.iatm1) then
                ibnd1 = i-1
                found =.true.
              endif

              if(found) goto 320

            enddo

 320        continue
c
c     search for constrained bonds -

            if(hcons)then
              do i = ibnd0,ibnd1
                nami = ambnam(ibda(i))
                namj = ambnam(ibdb(i))

                if(nami(1:1).eq.'H'.or.namj(1:1).eq.'H') ih0=ih0+1

              enddo

            endif

        if(.not.allcons) then

            ibnd2 = 0
            if(ibnd1.gt.0) ibnd2 = ibnd1-ibnd0+1-ih0
            if(ibnd2.gt.0) write(ifield,*) 'bonds ',ibnd2
            do i = ibnd0,ibnd1

              ia = min(ibda(i),ibdb(i))
              ib = max(ibda(i),ibdb(i))
              call bndfind(ifield,hcons,ia,ib,ibc,ibnd,mgb,cb,bo,
     x             bnda,bndb,ambnam,x,y,z,iatm0,bndtyp)

            enddo
            if(.not.hcons) go to 312
          endif
c
c     determine bond constraints and print out

            if(allcons) ih0 = ibnd1-ibnd0+1
            if(ih0.gt.0) then
              write(ifield,*) 'constraints ', ih0

              do i = ibnd0,ibnd1

                ia = min(ibda(i),ibdb(i))
                ib = max(ibda(i),ibdb(i))
                call consfind(ifield,hcons,allcons,ia,ib,i,icon,
     x               mgb,bo,bnda,bndb,ambnam,x,y,z,iatm0)

              enddo

            endif
  312     continue
c
c     Find valence angle arrays from bond angle arrays
c
c     initialise counter and arrays

          iang = 0
c
c     apply a patch for tip3p water.

          found = .true.
          if((ambnam(iatm0+1).eq.water).and.(ibnd1-ibnd0+1.eq.3))
     x         found = .false.

          if(found) then

            do i = iatm0+1,iatm1

              if(ambnam(i).ne.hwater) then

                do j = 1,nbds(i)-1

                  ij = nexatm(i,j)

                  do k = j+1, nbds(i)

                    ik = nexatm(i,k)

                    iang = iang + 1

                  enddo
                enddo
              endif
            enddo
          endif
          if(iang.gt.0)  then
            write(ifield,*) 'angles ',iang

            iang = 0
            do i = iatm0+1,iatm1

              do j = 1,nbds(i)-1

                ij = nexatm(i,j)

                do k = j+1, nbds(i)

                  ik = nexatm(i,k)

                  iang = iang + 1
                  langa(iang) = min(ij,ik)
                  langb(iang) = i
                  langc(iang) = max(ij,ik)

                  call angfind(ifield,mga,iang,langa(iang),i,
     x                 langc(iang),anga,angb,angc,ct,to,ambnam,iatm0,
     x                 angtyp)

                enddo

              enddo

            enddo

          endif

c
c     generate dihedral angle arrays - using angle arrays and bond arrays

          idih = 0
          do i = 1,mxang
            ldihb(i) = 0
            ldihc(i) = 0
          enddo
          do iloop = 1,2

            if(iloop.eq.2.and.idih.gt.0)
     x           write(ifield,*) 'dihedrals ',idih

            idmax = 0
            idih = 0
            do i = 1,iang

              ia = langa(i)
              itpa = ityp(ia)

              ib = langb(i)
              itpb = ityp(ib)

              ic = langc(i)
              itpc = ityp(ic)
c
c     check for dihedral id-- ia -- ib -- ic

              nloop = nbds(ia)

              do j = 1,nloop

                id = nexatm(ia,j)

                found = .false.
                if(ib.eq.id) found = .true.
                if(ic.eq.id) found = .true.
                if(id.le.0) found = .true.

                if(.not.found) then
                  idmax = idmax+1

                  l14a(idmax)=id
                  l14b(idmax)=ia
                  l14c(idmax)=ib
                  l14d(idmax)=ic
                  if(ffnam.eq.'gromos  ') then
c
c     have only 1 dihedral for each central pair
c     unless are dealing with sugar carbons

                    do i9 = 1,idih

                      ia1 = ldihb(i9)
                      ib1 = ldihc(i9)

                      if(ia1.eq.ia) then
                        if(ib1.eq.ib) found = .true.
                      elseif(ia1.eq.ib) then
                        if(ib1.eq.ia) found = .true.
                      endif

                    enddo

                    isug = 0
                    if(found) then
                      if(ambnam(ia)(1:2).eq.'CS'.and.
     x                     ambnam(ib)(1:2).eq.'CS') then
                        found = .false.
                        isug = 1
                      endif
                    endif
                  endif
                endif
                if(.not.found) then
                  do i9 = 1,idmax-1

                    ia1 = l14b(i9)
                    ib1 = l14c(i9)
                    ic1 = l14d(i9)
                    id1 = l14a(i9)

                    if(ia1.eq.ia) then
                      if(ib1.eq.ib) then
                        if(ic1.eq.ic) then
                          if(id1.eq.id) found = .true.
                        endif
                      endif
                    elseif(ia1.eq.ib) then
                      if(ib1.eq.ia) then
                        if(ic1.eq.id) then
                          if(id1.eq.ic) found = .true.
                        endif
                      endif
                    endif
                  enddo
                endif

                if(.not.found) then

                  ifnd = 0
c
c     search for specific interactions

                  call dihfind(ifield,mgd,idih,mxatm,ndis,id,ia,ib,ic,
     x                 cqs,qds,tds,nps,disa,disb,disc,disd,ambnam,ifnd,
     x                 iatm0,iloop,dihtyp,delec,dvdw,isug)
c
c     search for general interactions

                  if(ifnd.eq.0) then

                    call dihfind
     x                   (ifield,mgd,idih,mxatm,ndi,id,ia,ib,ic,cq,qd,
     x                   td,np,dia,dib,dic,did,ambnam,ifnd,iatm0,iloop,
     x                   dihtyp,delec,dvdw,isug)

                  endif
c
c     store central pair of found dihedrals

                  if(ifnd.gt.0) then

                    ldihb(idih)= ia
                    ldihc(idih)= ib
                  endif


                  if(ifnd.eq.0) then
c
c     try improper dihedrals

                    idih1 = idih
                    id1 =id
                    ic1 = ic
                    ib1 =ib
                    ia1 = ia
                    call imdfind(ifield,nims,id1,ia1,ib1,ic1,aimps,
     x                   bimps,cimps,ambnam,imsa,imsb,imsc,imsd,
     x                   idih,found,iatm0,iloop,diityp,dielec,divdw,
     x                   x,y,z,cell,imcon)

                    id1 =id
                    ic1 = ic
                    ib1 =ib
                    ia1 = ia

                    if(idih1.eq.idih)
     x                   call imdfind(ifield,nim,id1,ia1,ib1,ic1,aimp,
     x                   bimp,cimp,ambnam,ima,imb,imc,imd,idih,found,
     x                   iatm0,iloop,diityp,dielec,divdw,x,y,z,
     x                   cell,imcon)

c     print warning for undefined parameters

                    if(iloop.eq.2.and.idih1.eq.idih) then

                      if(ffnam.ne.'gromos  ') then
                        idih = idih+1
                        write(ier,'(a4,4i5,a36,2g12.5,i10)')dihtyp,
     x                       id-iatm0,ia-iatm0,ib-iatm0,ic-iatm0,
     x                       ' ****   unknown   parameters   **** ',
     x                       delec,dvdw,idih

                      endif

                    endif

                  endif

                endif

              enddo
c
c     check for dihedral ia -- ib -- ic -- id

              nloop = nbds(ic)

              do j = 1,nloop

                id = nexatm(ic,j)

                found = .false.
                if(ib.eq.id) found = .true.
                if(ia.eq.id) found = .true.
                if(id.le.0) found = .true.

                if(.not.found) then
                  idmax = idmax+1
                  l14a(idmax)=ia
                  l14b(idmax)=ib
                  l14c(idmax)=ic
                  l14d(idmax)=id

                  if(ffnam.eq.'gromos  ') then

                    do i9 = 1,idih

                      ib1 = ldihb(i9)
                      ic1 = ldihc(i9)

                      if(ib1.eq.ib) then
                        if(ic1.eq.ic) found = .true.
                      elseif(ib1.eq.ic) then
                        if(ic1.eq.ib) found = .true.
                      endif

                    enddo
                    isug = 0
                    if(found) then
                      if(ambnam(ib)(1:2).eq.'CS'.and.
     x                     ambnam(ic)(1:2).eq.'CS') then
                        found = .false.
                        isug = 1
                      endif
                    endif
                  endif
                endif
                if(.not.found) then
                  do i9 = 1,idmax-1

                    ia1 = l14a(i9)
                    ib1 = l14b(i9)
                    ic1 = l14c(i9)
                    id1 = l14d(i9)

                    if(ia1.eq.ia) then
                      if(ib1.eq.ib) then
                        if(ic1.eq.ic) then
                          if(id1.eq.id) found = .true.
                        endif
                      endif
                    elseif(ia1.eq.id) then
                      if(ib1.eq.ic) then
                        if(ic1.eq.ib) then
                          if(id1.eq.ia) found = .true.
                        endif
                      endif
                    endif
                  enddo
                endif

                if(.not.found) then

                  ifnd = 0
c
c     search for specific interactions

                  call dihfind(ifield,mgd,idih,mxatm,ndis,ia,ib,ic,id,
     x                 cqs,qds,tds,nps,disa,disb,disc,disd,ambnam,ifnd,
     x                 iatm0,iloop,dihtyp,delec,dvdw,isug)
c
c     search for general interactions

                  if(ifnd.eq.0) then

                    call dihfind
     x                   (ifield,mgd,idih,mxatm,ndi,ia,ib,ic,id,cq,qd,
     x                   td,np,dia,dib,dic,did,ambnam,ifnd,iatm0,iloop,
     x                   dihtyp,delec,dvdw,isug)

                  endif
c
c     store central pair of found dihedrals

                  if(ifnd.gt.0) then
                    ldihb(idih)= ib
                    ldihc(idih)= ic
                  endif

                  if(ifnd.eq.0) then
c
c     try improper dihedrals

                    idih1 = idih
                    id1 =id
                    ic1 = ic
                    ib1 =ib
                    ia1 = ia
                    call imdfind(ifield,nims,ia1,ib1,ic1,id1,aimps,
     x                   bimps,cimps,ambnam,imsa,imsb,imsc,imsd,
     x                   idih,found,iatm0,iloop,diityp,dielec,divdw,
     x                   x,y,z,cell,imcon)

                    id1 =id
                    ic1 = ic
                    ib1 =ib
                    ia1 = ia

                    if(idih1.eq.idih)
     x                   call imdfind(ifield,nim,ia1,ib1,ic1,id1,aimp,
     x                   bimp,cimp,ambnam,ima,imb,imc,imd,idih,found,
     x                   iatm0,iloop,diityp,dielec,divdw,x,y,z,
     x                   cell,imcon)

c     print warning for undefined parameters

                    if(iloop.eq.2.and.idih1.eq.idih) then

                      if(ffnam.ne.'gromos  ') then
                        idih = idih+1
                        write(ier,'(a4,4i5,a36,2g12.5,i10)')dihtyp,
     x                       ia-iatm0,ib-iatm0,ic-iatm0,id-iatm0,
     x                       ' ****   unknown   parameters   **** ',
     x                       delec,dvdw,idih

                      endif

                    endif

                  endif

                endif

              enddo

            enddo

c
c     Improper dihedral force field

            do i = iatm0+1,iatm1

              ldo = .false.
              if(nbds(i).eq.3) ldo = .true.

              if(ldo) then

                if(ffnam.eq.'gromos  ') then
                  ia = i
                  ib = nexatm(i,1)
                  ic = nexatm(i,2)
                  id = nexatm(i,3)

                elseif(ffnam.eq.'amber   ') then
                  ic = i
                  ia = nexatm(i,1)
                  ib = nexatm(i,2)
                  id = nexatm(i,3)
                endif

                irot=0
                irot1=0

 195            continue
                call imdfind(ifield,nims,ia,ib,ic,id,aimps,bimps,cimps,
     x               ambnam,imsa,imsb,imsc,imsd,idih,found,iatm0,iloop,
     x               diityp,dielec,divdw,x,y,z,cell,imcon)

                if (.not.found) call imdfind
     x               (ifield,nim,ia,ib,ic,id,aimp,bimp,cimp,
     x               ambnam,ima,imb,imc,imd,idih,found,iatm0,iloop,
     x               diityp,dielec,divdw,x,y,z,cell,imcon)

                if(.not.found) then
c
c     rotate indices around central atom if search unsuccessful

                  if(ffnam.eq.'amber   ') then
c
c     amber
                    if(irot.lt.2) then

                      irot = irot + 1
                      ie = ia
                      ia = ib
                      ib = id
                      id = ie
                      is = 0

                    else

                      irot = 0
                      irot1 = irot1 + 1

                      ie = ia
                      ia = ib
                      ib = ie

                    endif
                  elseif(ffnam.eq.'gromos  ') then
c
c     gromos
                    if(irot.lt.2) then

                      irot = irot + 1
                      ie = ib
                      ib = ic
                      ic = id
                      id = ie

                    else

                      irot = 0
                      irot1 = irot1 + 1

                      ie = ic
                      ib = ic
                      ic = ie

                    endif

                  endif
                  if(irot1.le.1) goto 195
                endif

              endif

            enddo
          enddo

c
c     special 1..4 bonds (12-6) type

        if(ffnam.eq.'gromos  ') then
          if(idmax.gt.0) then
            ncc = 0
            do iloop = 1,2
              if(ncc.gt.0) write(ifield,'(a,i10)')'bonds ',ncc

              ncc = 0
              do i = 1,idmax
                ia = min(l14a(i),l14d(i))
                ib = max(l14a(i),l14d(i))
                found = .false.
c
c     check for redundancies
                do j = 1,i-1
                  ja = l14a(j)
                  jb = l14d(j)
                  if(ia.eq.ja) then
                    if(ib.eq.jb) found = .true.
                  elseif(ia.eq.jb) then
                    if(ib.eq.ja) found = .true.
                  endif
                enddo

c
c     check for 4 or 5 membered rings ....
                if(.not.found) then
                  do j = 1,iang
                    if(ia.eq.langa(j)) then
                      if(ib.eq.langb(j)) found = .true.
                      if(ib.eq.langc(j)) found = .true.
                    elseif(ia.eq.langc(j)) then
                      if(ib.eq.langb(j)) found = .true.
                      if(ib.eq.langa(j)) found = .true.
                    endif
                  enddo
                endif
c
c     got one!
                if(.not.found) then

                  nami = ambnam(ia)
                  namj = ambnam(ib)
c
c     set standard lj terms

                  do j = 1,ilj
                    if(nami.eq.nlj(j)) then
                      a12 = c12(j)
                      a6  = c6(j)
                    endif
                    if(namj.eq.nlj(j)) then
                      b12 = c12(j)
                      b6  = c6(j)
                    endif
                  enddo
c
c     set special cross terms

                  do j = 1,ljx
                    if(nami.eq.nlja(j)) then
                      if(namj.eq.nljb(j)) then
                        a12 = sqrt(c12a(j))
                        a6 =  sqrt(c6a(j))
                        b12 = sqrt(c12a(j))
                        b6  = sqrt(c6a(j))
                      endif
                    elseif(namj.eq.nlja(j)) then
                      if(nami.eq.nljb(j)) then
                        a12 = sqrt(c12a(j))
                        a6 =  sqrt(c6a(j))
                        b12 = sqrt(c12a(j))
                        b6  = sqrt(c6a(j))
                      endif
                    endif
                  enddo
c
c     set smaller values for united atoms

                  do j = 1,lj14
                    if(nami.eq.nljh(j)) then
                      a12 = c12b(j)
                      a6  = c6b(j)
                    endif
                    if(namj.eq.nljh(j)) then
                      b12 = c12b(j)
                      b6  = c6b(j)
                    endif
                  enddo

                  if(a12*b12.gt.0.d0) then
                    ncc = ncc+1

                    if(iloop.eq.2)
     x                   write(ifield,'(a4,2i5,f12.3,f12.4,i10)')
     x                   '12-6',ia,ib,a12*b12,a6*b6,ncc
                  endif
                endif
              enddo
            enddo
          endif
        endif
        write(ifield,*) 'finish'
      enddo

c************************************************************************
c
c     generate non-bonded interaction arrays
c
c************************************************************************


      do i = 1,ilj

c     check for appearance of atom (or equivalents) in atom list

        ijchk(i) = 0

      enddo

      do i = 1,iatm

        nami = ambnam(i)

        do j = 1,natm

          lj = ljeqiv(j)

          if (name(j).eq.nami) ijchk(lj) = 1

        enddo

      enddo

      llj = 0

      do i = 1,natm

        if(ijchk(i).eq.1) llj = llj + 1

      enddo
c
c     number of short range forces

      llj = llj*(llj+1)/2
      write(ifield,*) 'vdw ',llj


      llj = 0
      do i = 1,natm

        if (ijchk(i).eq.1) then

          ik = 0

          do k = 1,ilj

            if(name(i).eq.nlj(k)) ik = k

          enddo

          do j = i,natm

            if(ijchk(j).eq.1) then

              jk = 0

              do k = 1,ilj

                if(name(j).eq.nlj(k)) jk = k

              enddo

              llj = llj + 1

c
c     check for Hydrogen bonds also

              if(jk.gt.0.and.ik.gt.0) then
                call hbond (nlj(ik),nlj(jk),hba,hbb,
     x               is,ihb)

                if(is.eq.0) then

                  if(vdtyp(1:2).eq.'lj') then
                    write(ifield,'(2(a4,4x),a4,2g12.5,3x,i6)')
     x                   nlj(ik),nlj(jk),'lj  ',
     x                   sqrt(c6(ik)*c6(jk)),
     x                   c12(ik)+c12(jk),llj

                  else

                    a12 = c12(ik)
                    a6 = c6(ik)
                    b12= c12(jk)
                    b6 = c6(jk)
c
c     set special cross terms

                    do j9 = 1,ljx
                      if(name(i).eq.nlja(j9)) then
                        if(name(j).eq.nljb(j9)) then
                          a12 = sqrt(c12a(j9))
                          a6 =  sqrt(c6a(j9))
                          b12 = sqrt(c12a(j9))
                          b6  = sqrt(c6a(j9))
                        endif
                      elseif(name(j).eq.nlja(j9)) then
                        if(name(i).eq.nljb(j9)) then
                          a12 = sqrt(c12a(j9))
                          a6 =  sqrt(c6a(j9))
                          b12 = sqrt(c12a(j9))
                          b6  = sqrt(c6a(j9))
                        endif
                      endif
                    enddo

                    write(ifield,'(2(a4,4x),a4,f12.3,f12.5,3x,i6)')
     x                   nlj(ik),nlj(jk),'12-6',
     x                   a12*b12,a6*b6,llj

                  endif

                else

                  write(ifield,'(2(a4,4x),a4,2g12.5,3x,i6)')
     x                 nlj(ik),nlj(jk),'hbnd',h12(is),h10(is),llj

                endif

              endif

            endif

          enddo

        endif

      enddo

      write(ifield,*) 'close'
      close(ifield)
      close(ier)

      end


      subroutine error(i,j)

      ifield = 6
      if (i.eq.1) write(ifield,1) j
      if (i.eq.2) write(ifield,2) j
      if (i.eq.3) write(ifield,3) j
      if (i.eq.4) write(ifield,4) j
      if (i.eq.5) write(ifield,5) j
      if (i.eq.6) write(ifield,6) j
      if (i.eq.7) write(ifield,7) j
      if (i.eq.8) write(ifield,8) j
      if (i.eq.9) write(ifield,9) j
      if (i.eq.10) write(ifield,10) j


 1    format(3x,'*** - amino acid unrecognised - ***',i5)
 2    format(3x,'*** - hydrogen miscount on amino acid - ***',i5)
 3    format(3x,'*** - bndfind failed - ***',i6)
 4    format(3x,'*** - too many bonds in system : ',i9,' ***')
 5    format(3x,'*** - angfind failed - ***', i6)
 6    format(3x,'*** - dihfind failed - ***', i6)
 7    format(3x,i6)
 8    format(3x,'*** - LJ atom unrecognised',i6)
 9    format(3x,'*** - LJ map unrecognised ',i6)
 10   format(3x,'*** - too many molecular types in system: found ',i6)

      if(i.eq.2) goto 999
      if(i.eq.3) goto 999
      if(i.eq.5) goto 999
      if(i.eq.6) goto 999
      if(i.eq.8) goto 999
      if(i.eq.9) goto 999
      stop
 999  return
      end

      subroutine block(inum,ird2)
c
c     subroutine to read past comment statements and find out how many
c     entries follow

      character*5 blck

  999 read(ird2,10) blck

      if (blck.ne."BLOCK") goto 999

      read(ird2,*) inum


   10 format(a5)

      return
      end

      subroutine ambtype(at,names,ntyp,i)

c     identifies the type of atom

      character*4 at, names(*)


   20 format(' atom name : ',a4, ' not recognized, type set to ',a4,
     x     ' index = ',i4)


      if(at.eq."NA+ ") at = "IP  "
      i = 0
   10 i = i + 1

      if(at.eq.names(i)) goto 999

      if (i.lt.ntyp) goto 10


c
c     atom name not identified -- signal an error


      i = 1
      write(20,20) at,names(i), i

  999 return
      end

      subroutine bndtitle(ird1,ierr,res)

      implicit real*8(a-h,o-z)
      character*60 oneline
      character*3 res

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
        res = oneline(1:3)
      endif
      if (oneline(6:15).ne."BOND ARRAY") goto 10

      read(ird1,*)
      return
c
c     end of file reached

   20 ierr = 1

      return
      end

      subroutine bnddata(ird1,ifield,iatm,ibnd,kres,
     x     mgb,cb,bo,bnda,bndb,name,ambnam,ibda,ibdb,ires,x,y,z,qq,
     x     res,resn)


      implicit real*8(a-h,o-z)

      character*3 res,resn(*)
      character*4 bnda(*),bndb(*),ambnam(*),an,at
      character*8 name(*)
      dimension cb(*),bo(*)
      dimension x(*),y(*),z(*),qq(*)
      dimension ibda(*),ibdb(*),ires(*)

   10 read(ird1,'(2i5,3x,a4,2x,a4,6x,3f10.4,f8.4)',err=20,end=30)
     x      ia,ib,at,an,x1,y1,z1,q1

      if (ia.gt.0)  then

         ambnam(ia) = an
         qq(ia) = q1

         if(ia.gt.iatm) then

            resn(ia) = res
            name(ia) = at
            iatm = ia
            x(ia) = x1
            y(ia) = y1
            z(ia) = z1
            ires(ia) = kres

         endif


         if (ib.gt.0) then

            ibnd = ibnd + 1

            ibda(ibnd) = ia
            ibdb(ibnd) = ib

         endif

      endif

      goto 10

 20   backspace(ird1)
 30   return
      end

      subroutine bndfind(ifield,hcons,ia,ib,i,ibnd,
     x     mgb,cb,bo,bnda,bndb,name,x,y,z,iatm0,bndtyp)

      implicit real*8(a-h,o-z)
      logical found,hcons,wanted

      character*4 bnda(*),bndb(*),name(*),nami,namj,bndtyp

      dimension cb(*),bo(*)
      dimension x(*),y(*),z(*)

      found = .false.
      wanted = .true.
c
c     check if harmonic bonds required for this bond

      if(hcons) then

         nami = name(ia)
         namj = name(ib)

         if(nami(1:1).eq."H".or.namj(1:1).eq."H") wanted = .false.

      endif

      if (wanted) then

         is = 0
   10    is = is + 1

         if(bnda(is).eq.name(ia)) then

            if (bndb(is).eq.name(ib)) found = .true.

         else

            if(bndb(is).eq.name(ia)) then

               if (bnda(is).eq.name(ib)) found = .true.

            endif

         endif

         if ((is.eq.mgb).and.(.not.found)) then

            write(6,*) name(ia),name(ib)
            write(6,*) ia-iatm0,ib-iatm0
            call error (3,ibnd)
            found = .true.
            is = -1

         endif

         if (.not.found) goto 10

         i = i+1
c
c     write out bond data in FIELD file

         if(is.gt.0) then
           write(ifield,'(a4,2i5,g12.6,f12.5,3x,i6)') bndtyp,
     x      ia-iatm0,ib-iatm0,cb(is),bo(is),i
         else
           write(ifield,'(a4,2i5,2a12,3x,i6)') bndtyp,
     x      ia-iatm0,ib-iatm0,' unknown k  ',' and r0     ',i
         endif
      endif
      return
      end

      subroutine consfind(ifield,hcons,allcons,ia,ib,i,icon,
     x  mgb,bo,bnda,bndb,ambnam,x,y,z,iatm0)

      implicit real*8(a-h,o-z)

      logical hcons,allcons,wanted,found
      dimension x(*),y(*),z(*)
      character*4 bnda(*),bndb(*),ambnam(*),nami,namj
      dimension bo(*)

      found = .false.
      wanted = .false.

      if(allcons) then

         wanted = .true.

      else

         nami = ambnam(ia)
         namj = ambnam(ib)

         if(hcons.and.nami(1:1).eq."H") wanted = .true.
         if(hcons.and.namj(1:1).eq."H") wanted = .true.

      endif
c
c     constraint bonds

      if(wanted) then

         is = 0
   10    is = is + 1

         if(bnda(is).eq.ambnam(ia)) then

            if (bndb(is).eq.ambnam(ib)) found = .true.

         else

            if(bndb(is).eq.ambnam(ia)) then

               if (bnda(is).eq.ambnam(ib)) found = .true.

            endif

         endif

         if ((is.eq.mgb).and.(.not.found)) then

            write(6,*) ambnam(ia),ambnam(ib)
            write(6,*) ia-iatm0,ib-iatm0
            call error (3,ibnd)
            found = .true.
            is = -1

         endif

         if (.not.found) goto 10

         icon = icon+1
c
c     write out bond data in FIELD file

         if(is.gt.0) then
           write(ifield,'(2i5,f12.5,3x,i6)')
     x          ia-iatm0,ib-iatm0,bo(is),icon
         else
           write(ifield,'(2i5,a12,3x,i6)')
     x          ia-iatm0,ib-iatm0,' unknown r0 ',icon
         endif

      endif
      return
      end

      subroutine angfind(ifield,mga,iang,ia,ib,ic,
     x  anga,angb,angc,ct,to,name,iatm0,angtyp)

      implicit real*8 (a-h,o-z)
      character*4 name(*),anga(*),angb(*),angc(*),angtyp
      dimension ct(*),to(*)

      logical found

      found = .false.
c
c     apply patch for HW as central atom

      if(name(ib).eq.'HW  ') goto 20

      is = 0
   10 is = is + 1

      if(angb(is).eq.name(ib)) then

        if (anga(is).eq.name(ia)) then

          if(angc(is).eq.name(ic)) found = .true.

        else if(angc(is).eq.name(ia)) then

          if (anga(is).eq.name(ic)) found = .true.

        endif

      endif

      if ((is.eq.mga).and.(.not.found)) then

        write(6,*) name(ia),name(ib),name(ic)
        write(6,*) ia-iatm0,ib-iatm0,ic-iatm0
        call error (5,iang)
        found = .true.
        is = -1

      endif

      if (.not.found) goto 10

      if(is.gt.0) then
        write(ifield,'(a4,3i5,2g12.5,3x,i6)') angtyp,ia-iatm0,
     x  ib-iatm0,ic-iatm0,ct(is),to(is),iang
      else
        write(ifield,'(a4,3i5,2a12,3x,i6)') angtyp,ia-iatm0,
     x  ib-iatm0,ic-iatm0,' unknown k  ',' and theta  ',iang
      endif

   20 return

      end
      subroutine dihfind(ifield,mgd,idih,mxatm,ndi,ia,ib,ic,id,
     x  cq,qd,td,np,dia,dib,dic,did,name,ifnd,iatm0,iloop,dihtyp,
     x  delec,dvdw,isug)


      implicit real*8(a-h,o-z)
      character*4 name(mxatm),dia(mgd),dib(mgd),dic(mgd),did(mgd)
      character*4 wild,dihtyp
      dimension cq(mgd),qd(mgd),td(mgd)
      dimension np(mgd)

      logical fndnow

      wild = "*   "

      is = 0
   10 is = is + 1

      fndnow = .false.

      if((dia(is).eq.name(ia)).or.(dia(is).eq.wild)) then
        if((dib(is).eq.name(ib)).or.(dib(is).eq.wild)) then
          if((dic(is).eq.name(ic)).or.(dic(is).eq.wild)) then
            if((did(is).eq.name(id)).or.(did(is).eq.wild)) then

              fndnow = .true.

            endif
          endif
        endif
      endif

      if((dia(is).eq.name(id)).or.(dia(is).eq.wild)) then
        if((dib(is).eq.name(ic)).or.(dib(is).eq.wild)) then
          if((dic(is).eq.name(ib)).or.(dic(is).eq.wild)) then
            if((did(is).eq.name(ia)).or.(did(is).eq.wild)) then

              fndnow = .true.

            endif
          endif
        endif
      endif

      if(fndnow) then

        ifnd = ifnd+1
        idih = idih + 1

        ah = cq(is)/dble(np(is))

c
c     force field to include reduced electrostatic term

        vdw = dvdw
        coul = delec
        if(ifnd.gt.1)  then
          vdw = 0.d0
          coul =0.d0
        endif
c
c     patch for gromos sugar units!
        if(isug.gt.0.and.abs(ah-1.4d0).lt.1d-4) then

          idih = idih-1

        else

          if(ah.gt.0.d0.or.ifnd.eq.1) then
            if(iloop.eq.2)
     x           write(ifield,'(a4,4i5,5f12.5,i10)') dihtyp,
     x           ia-iatm0,ib-iatm0,ic-iatm0,id-iatm0,
     x           ah,qd(is),td(is),coul,vdw,idih

          else

            idih = idih-1

          endif

        endif

      endif

      if(is.lt.ndi) goto 10

      return
      end

      subroutine units(e,x,ix)
      implicit real*8 (a-h,o-z)

      dimension x(*)

      do i = 1, ix

        x(i) = x(i)*e

      enddo

      return
      end

      subroutine hbond(nami,namj,hba,hbb,is,ihb)

      character*4 nami,namj, hba(*),hbb(*)
      logical found

      found = .false.
      is = 0
   10 is = is + 1

      if(nami.eq.hba(is).and.namj.eq.hbb(is)) found = .true.
      if(nami.eq.hbb(is).and.namj.eq.hba(is)) found = .true.

      if(found) goto 20

      if(is.lt.ihb.and.(.not.found)) goto 10

c
c     no Hydrogen bond here

      is = 0

   20 return
      end

      subroutine write1(conf,ird1,iatm,ljeqiv,iamb,
     x  x,y,z,ambnam,name,rname,cell,imcon)

      implicit real*8(a-h,o-z)
      character*80 oneline
      character*4 ambnam(*),name(*),rname(*),nami
      dimension x(*),y(*),z(*),cell(9)
      dimension ljeqiv(*),iamb(*)
      logical found

      write(6,*) 'writing the CONFIG file .....'
c
c     search for box vectors in EDTOUT

      imcon = 0
      do i = 1,9
        cell(i) = 0.d0
      enddo
      rewind(ird1)
      found = .false.
      do imega = 1,100000
        read(ird1,'(a80)',end=20,err=20) oneline
        if(oneline(6:9).eq.'BOXX') found = .true.
        if(found) goto 20
        if(oneline(1:12).eq.'CELL_VECTORS') found = .true.
        if(found) goto 25
      enddo

   20 if(found) then
        cell(1) = dblstr(oneline,80,idum)
        cell(5) = dblstr(oneline(idum:idum),80-idum,idum2)
        idum = idum + idum2
        cell(9) = dblstr(oneline(idum:idum),80-idum,idum2)
        imcon = 2

        write(conf,'(2i10)') 0,imcon
        write(conf,'(3f20.10)') cell

        found= .false.
      else
c    can't find cell vectors

        write(conf,'(2i10)') 0,0

      endif
 25   if(found) then
        imcon = 3
        read(ird1,*) cell
        write(conf,'(2i10)') 0,imcon
        write(conf,'(3f20.10)') cell
      endif
c
c     write all names and coordinates in CONFIG file

      do i = 1,iatm

         nami = name(ljeqiv(iamb(i)))

         write(conf,'(4x,a4,i5,3x,a4,3x,a4)')
     x       nami,i,ambnam(i)
         write(conf,'(3f20.8)') x(i),y(i),z(i)

      enddo
      close(conf)
      return
      end
      subroutine imdfind(ifield,nim,ia,ib,ic,id,aimp,bimp,cimp,
     x     name,ima,imb,imc,imd,idih,found,iatm0,iloop,diityp,
     x     dielec,divdw,x,y,z,cell,imcon)

      implicit real*8 (a-h,o-z)
      dimension aimp(*), bimp(*),cimp(*)
      character*4 name(*),ima(*),imb(*),imc(*),imd(*),wild,diityp
      dimension x(*),y(*),z(*),cell(*)

      logical found

      wild = "*   "

      found = .false.

      is = 0
 10   is = is + 1

      if(ima(is).eq.name(ia).or.ima(is).eq.wild) then
        if(imb(is).eq.name(ib).or.imb(is).eq.wild) then
          if(imc(is).eq.name(ic).or.imc(is).eq.wild) then
            if(imd(is).eq.name(id).or.imd(is).eq.wild) then

              found = .true.

            endif
          endif
        endif
      endif

      if(ima(is).eq.name(id).or.ima(is).eq.wild) then
        if(imb(is).eq.name(ic).or.imb(is).eq.wild) then
          if(imc(is).eq.name(ib).or.imc(is).eq.wild) then
            if(imd(is).eq.name(ia).or.imd(is).eq.wild) then

              found = .true.

            endif
          endif
        endif
      endif

      if(found.and.aimp(is).ne.0.d0) then

        idih = idih + 1

        if(iloop.eq.2) then

          bss = bimp(is)
          if(diityp.eq.'harm') then
            call dihang(ia,ib,ic,id,x,y,z,theta,cell,imcon)
            theta= theta*180.d0/3.1415926d0
            if(abs(bss).lt.1.d0) then
              bss = 180.d0*dnint(theta/180.d0)
            else
              bss = sign(bss,theta)
            endif
          endif

          write(ifield,'(a4,4i5,5f12.5,3x,i6)') diityp,
     x         ia-iatm0,ib-iatm0,ic-iatm0,id-iatm0,
     x         aimp(is),bss,cimp(is),dielec,divdw,idih

        endif
      endif


      if(is.lt.nim.and..not.found) goto 10

      return
      end

      subroutine concheck
     x     (ambnam,iatm0,iatm1,ibnd,
     x     iloop,mxatm,x,y,z,ibda,ibdb,nbds,nexatm)

      implicit real*8(a-h,o-z)

      logical now,found,safe
      dimension x(*),y(*),z(*)
      character*4 ambnam(*),nami
      dimension ibda(*),ibdb(*),nbds(*),nexatm(mxatm,10)
      dimension rmin(4), imin(4)
c
c     checks that the correct number of connections have been found for each
c     atom type

c     set up tolerancesfor bond lengths
      if(iloop.eq.1) then
        rl1 = 1.2d0**2
        rl2 = 1.7d0**2
      elseif(iloop.eq.2) then
        rl1 = 1.1d0**2
        rl2 = 1.9d0**2
      elseif(iloop.eq.3) then
        rl1 = 0.5d0**2
        rl2 = 3.0d0**2
      endif

      do i = iatm0,iatm1

        icn = 0
        nami = ambnam(i)

        call nconnect(icn,nami)
c
c     search for missing connections - taken to be the shortest
c     remaining distances to the central atom.

        iter = icn - nbds(i)
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

          do j = iatm0,iatm1

            now = .true.
            if(iloop.eq.1.and.ambnam(i)(1:1).eq.'H') now = .false.
            if(iloop.eq.1.and.ambnam(j)(1:1).eq.'H') now = .false.

            if (i.eq.j) now =.false.
c
c     reject fully bonded atoms

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

          do kk = 1,min(iter,4)
            ij = imin(kk)
            rr = sqrt(rmin(kk))
            if(ij.ne.0) then
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
     x             'adding neighbour to ',i,
     x           '(',nami,') : ',ij,'(',ambnam(ij),') : dist= ',rr

            else

              if(iloop.eq.3) write(6,*) 'missing neighbour for ',i,
     x             ' (',nami,')'

            endif

          enddo

        elseif (iter.lt.0) then
c
c     too many connections for this atom

          if(nami.ne.'OW  ')
     x         write(6,*)'too many neighbours for ',i,'(',nami,')',
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

      subroutine dihang(ia,ib,ic,id,x,y,z,theta,cell,imcon)
      implicit real*8(a-h,o-z)
      dimension x(*),y(*),z(*),cell(*)
c
c     determine dihedral angle ia-ib-ic-id

      xab=x(ia)-x(ib)
      yab=y(ia)-y(ib)
      zab=z(ia)-z(ib)

      xbc=x(ib)-x(ic)
      ybc=y(ib)-y(ic)
      zbc=z(ib)-z(ic)

      xcd=x(ic)-x(id)
      ycd=y(ic)-y(id)
      zcd=z(ic)-z(id)
c
c     periodic boundary condition

      call images(imcon,0,1,1,cell,xab,yab,zab)
      call images(imcon,0,1,1,cell,xbc,ybc,zbc)
      call images(imcon,0,1,1,cell,xcd,ycd,zcd)

      rrbc=1.d0/sqrt(xbc*xbc+ybc*ybc+zbc*zbc)

      xac=xab+xbc
      yac=yab+ybc
      zac=zab+zbc
c
c     construct first dihedral vector

      pbx=yab*zbc-zab*ybc
      pby=zab*xbc-xab*zbc
      pbz=xab*ybc-yab*xbc
      pb2=pbx*pbx+pby*pby+pbz*pbz
      rpb1=1.d0/sqrt(pb2)
      rpb2 =rpb1*rpb1
c
c     construct second dihedral vector

      pcx=ybc*zcd-zbc*ycd
      pcy=zbc*xcd-xbc*zcd
      pcz=xbc*ycd-ybc*xcd
      pc2=pcx*pcx+pcy*pcy+pcz*pcz
      rpc1=1.d0/sqrt(pc2)
      rpc2 = rpc1*rpc1
c
c     determine dihedral angle

      pbpc=pbx*pcx+pby*pcy+pbz*pcz
      cost=pbpc*rpb1*rpc1
      sint=(xbc*(pcy*pbz-pcz*pby)+ybc*(pbx*pcz-pbz*pcx)+
     x     zbc*(pcx*pby-pcy*pbx))*(rpb1*rpc1*rrbc)

      theta=atan2(sint,cost)

      return
      end

      subroutine images
     x  (imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
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
c     T3D optimised version. t.forester july 1994
c
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
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

      if(imcon.gt.0) then

c
c     block indices

        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif

      if(imcon.eq.1)then
c
c     standard cubic boundary conditions


        aaa=1.d0/cell(1)


        do i=iatm1,iatm2
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

        do i=iatm1,iatm2

          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))

        enddo

      else if(imcon.eq.3)then
c
c     parallelepiped boundary conditions

        call invert(cell,rcell,det)

        do i=iatm1,iatm2

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

        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)

        aaa=1.d0/cell(1)

        do i=iatm1,iatm2

          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))

          if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.
     x      (0.75d0*cell(1)))then

            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(1),zzz(i))

          endif

        enddo

      else if(imcon.eq.5)then
c
c     rhombic dodecahedral boundary conditions

        rt2=sqrt(2.d0)
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6))
     x    call error(idnode,140)

        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)

        do i=iatm1,iatm2

          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(bbb*zzz(i))

          if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge.
     x      cell(1))then

            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(9),zzz(i))

          endif

        enddo

      else if(imcon.eq.6) then
c
c     x-y boundary conditions

        det = cell(1)*cell(5) - cell(2)*cell(4)

        if(abs(det).lt.1.d-6)call error(idnode,120)

        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)

        do i=iatm1,iatm2

          ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
          ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

          xss = ssx - nint(ssx)
          yss = ssy - nint(ssy)

          xxx(i)=cell(1)*xss + cell(4)*yss
          yyy(i)=cell(2)*xss + cell(5)*yss

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
