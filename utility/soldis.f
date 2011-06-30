      program soldis

c*********************************************************************
c
c     dl_poly utility to calculate average solvation energy
c     distribution functions
c
c     copyright - daresbury laboratory
c     author    -   w.smith oct 2008
c
c*********************************************************************

      implicit none

      integer, parameter :: dim=128

      character header(80),units(40)
      logical lexcite,lswitch,lcomp(9),jcomp(9)
      integer nstep,natms,mxtmls,mxtmls2,mxtmls3,mxtmls4,numacc
      integer i,j,k
      real(8) tstep,elrc,elrc2,hdim,deld,delv
      real(8), allocatable :: tmpsol(:),tmpbot(:),tmptop(:)
      real(8), allocatable :: bndsol(:),angsol(:),dihsol(:)
      real(8), allocatable :: invsol(:),shlsol(:),cousol(:)
      real(8), allocatable :: vdwsol(:),en3sol(:),en4sol(:)
      real(8), allocatable :: bndexc(:),angexc(:),dihexc(:)
      real(8), allocatable :: invexc(:),shlexc(:),couexc(:)
      real(8), allocatable :: vdwexc(:),en3exc(:),en4exc(:)
      real(8), allocatable :: bndbot(:),angbot(:),dihbot(:)
      real(8), allocatable :: invbot(:),shlbot(:),coubot(:)
      real(8), allocatable :: vdwbot(:),en3bot(:),en4bot(:)
      real(8), allocatable :: bndbet(:),angbet(:),dihbet(:)
      real(8), allocatable :: invbet(:),shlbet(:),coubet(:)
      real(8), allocatable :: vdwbet(:),en3bet(:),en4bet(:)
      real(8), allocatable :: bndtop(:),angtop(:),dihtop(:)
      real(8), allocatable :: invtop(:),shltop(:),coutop(:)
      real(8), allocatable :: vdwtop(:),en3top(:),en4top(:)
      real(8), allocatable :: bndtep(:),angtep(:),dihtep(:)
      real(8), allocatable :: invtep(:),shltep(:),coutep(:)
      real(8), allocatable :: vdwtep(:),en3tep(:),en4tep(:)
      real(8), allocatable :: tmpdis(:,:),bnddis(:,:),angdis(:,:)
      real(8), allocatable :: dihdis(:,:),invdis(:,:),shldis(:,:)
      real(8), allocatable :: coudis(:,:),vdwdis(:,:),en3dis(:,:)
      real(8), allocatable :: en4dis(:,:)
      real(8), allocatable :: tmpdes(:,:),bnddes(:,:),angdes(:,:)
      real(8), allocatable :: dihdes(:,:),invdes(:,:),shldes(:,:)
      real(8), allocatable :: coudes(:,:),vdwdes(:,:),en3des(:,:)
      real(8), allocatable :: en4des(:,:)

      hdim=dble(dim-1)

c     open main input SOLVAT file

      open(7,file="SOLVAT")

c     read header data

      read(7,'(80a1)')header
      write(*,'(80a1)')header
      read(7,'(40a1)')units
      write(*,'(40a1)')units
      read(7,'(2i10)')natms,mxtmls
      write(*,'(2i10)')natms,mxtmls
      read(7,*)
      read(7,'(11l4)')lexcite,lswitch,lcomp

c     set working array dimensions

      jcomp(:)=lexcite.and.lcomp(:)
      mxtmls2=((mxtmls+1)*mxtmls)/2
      mxtmls3=(((mxtmls+3)*mxtmls+2)*mxtmls)/6
      mxtmls4=((((mxtmls+6)*mxtmls+11)*mxtmls+6)*mxtmls)/24

c     create required data arrays

      allocate(tmpsol(mxtmls),tmpbot(mxtmls),tmptop(mxtmls))
      allocate(tmpdis(dim,mxtmls))
      tmpdis(:,:)=0.d0
      if(lcomp(1))then
        allocate(bndsol(mxtmls),bndbot(mxtmls),bndtop(mxtmls))
        allocate(bnddis(dim,mxtmls))
        bnddis(:,:)=0.d0
      endif
      if(lcomp(2))then
        allocate(angsol(mxtmls),angbot(mxtmls),angtop(mxtmls))
        allocate(angdis(dim,mxtmls))
        angdis(:,:)=0.d0
      endif
      if(lcomp(3))then
        allocate(dihsol(mxtmls),dihbot(mxtmls),dihtop(mxtmls))
        allocate(dihdis(dim,mxtmls))
        dihdis(:,:)=0.d0
      endif
      if(lcomp(4))then
        allocate(invsol(mxtmls),invbot(mxtmls),invtop(mxtmls))
        allocate(invdis(dim,mxtmls))
        invdis(:,:)=0.d0
      endif
      if(lcomp(5))then
        allocate(shlsol(mxtmls),shlbot(mxtmls),shltop(mxtmls))
        allocate(shldis(dim,mxtmls))
        shldis(:,:)=0.d0
      endif
      if(lcomp(6))then
        allocate(cousol(mxtmls2),coubot(mxtmls2),coutop(mxtmls2))
        allocate(coudis(dim,mxtmls2))
        coudis(:,:)=0.d0
      endif
      if(lcomp(7))then
        allocate(vdwsol(mxtmls2),vdwbot(mxtmls2),vdwtop(mxtmls2))
        allocate(vdwdis(dim,mxtmls2))
        vdwdis(:,:)=0.d0
      endif
      if(lcomp(8))then
        allocate(en3sol(mxtmls3),en3bot(mxtmls3),en3top(mxtmls3))
        allocate(en3dis(dim,mxtmls3))
        en3dis(:,:)=0.d0
      endif
      if(lcomp(9))then
        allocate(en4sol(mxtmls4),en4bot(mxtmls4),en4top(mxtmls4))
        allocate(en4dis(dim,mxtmls4))
        en4dis(:,:)=0.d0
      endif
      if(jcomp(1))then
        allocate(bndexc(mxtmls),bndbet(mxtmls),bndtep(mxtmls))
        allocate(bnddes(dim,mxtmls))
        bnddes(:,:)=0.d0
      endif
      if(jcomp(2))then
        allocate(angexc(mxtmls),angbet(mxtmls),angtep(mxtmls))
        allocate(angdes(dim,mxtmls))
        angdes(:,:)=0.d0
      endif
      if(jcomp(3))then
        allocate(dihexc(mxtmls),dihbet(mxtmls),dihtep(mxtmls))
        allocate(dihdes(dim,mxtmls))
        dihdes(:,:)=0.d0
      endif
      if(jcomp(4))then
        allocate(invexc(mxtmls),invbet(mxtmls),invtep(mxtmls))
        allocate(invdes(dim,mxtmls))
        invdes(:,:)=0.d0
      endif
      if(jcomp(5))then
        allocate(shlexc(mxtmls),shlbet(mxtmls),shltep(mxtmls))
        allocate(shldes(dim,mxtmls))
        shldes(:,:)=0.d0
      endif
      if(jcomp(6))then
        allocate(couexc(mxtmls2),coubet(mxtmls2),coutep(mxtmls2))
        allocate(coudes(dim,mxtmls2))
        coudes(:,:)=0.d0
      endif
      if(jcomp(7))then
        allocate(vdwexc(mxtmls2),vdwbet(mxtmls2),vdwtep(mxtmls2))
        allocate(vdwdes(dim,mxtmls2))
        vdwdes(:,:)=0.d0
      endif
      if(jcomp(8))then
        allocate(en3exc(mxtmls3),en3bet(mxtmls3),en3tep(mxtmls3))
        allocate(en3des(dim,mxtmls3))
        en3des(:,:)=0.d0
      endif
      if(jcomp(9))then
        allocate(en4exc(mxtmls4),en4bet(mxtmls4),en4tep(mxtmls4))
        allocate(en4des(dim,mxtmls4))
        en4des(:,:)=0.d0
      endif

c     read periodic data

      numacc=0

      do while(.true.)

c     skip first record

        read(7,*,end=100)

c     statistical counters

        numacc=numacc+1

c     read species temperatures and and find max and min

        read(7,*)tmpsol(:)
        if(numacc.eq.1)then
          tmptop(:)=tmpsol(:)
          tmpbot(:)=tmpsol(:)
        endif
        tmptop(:)=max(tmptop(:),tmpsol(:))
        tmpbot(:)=min(tmpbot(:),tmpsol(:))

c     read species bond energies and find max and min

        if(lcomp(1))then
          read(7,*)bndsol(:)
          if(numacc.eq.1)then
            bndtop(:)=bndsol(:)
            bndbot(:)=bndsol(:)
          endif
          bndtop(:)=max(bndtop(:),bndsol(:))
          bndbot(:)=min(bndbot(:),bndsol(:))
        endif
        if(jcomp(1))then
          read(7,*)bndexc(:)
          if(numacc.eq.1)then
            bndtep(:)=bndexc(:)
            bndbet(:)=bndexc(:)
          endif
          bndtep(:)=max(bndtep(:),bndsol(:))
          bndbet(:)=min(bndbet(:),bndsol(:))
        endif

c     read species valence angle energies and find max and min

        if(lcomp(2))then
          read(7,*)angsol(:)
          if(numacc.eq.1)then
            angtop(:)=angsol(:)
            angbot(:)=angsol(:)
          endif
          angtop(:)=max(angtop(:),angsol(:))
          angbot(:)=min(angbot(:),angsol(:))
        endif
        if(jcomp(2))then
          read(7,*)angexc(:)
          if(numacc.eq.1)then
            angtep(:)=angexc(:)
            angbet(:)=angexc(:)
          endif
          angtep(:)=max(angtep(:),angexc(:))
          angbet(:)=min(angbet(:),angexc(:))
        endif

c     read species dihedral energies and find max and min

        if(lcomp(3))then
          read(7,*)dihsol(:)
          if(numacc.eq.1)then
            dihtop(:)=dihsol(:)
            dihbot(:)=dihsol(:)
          endif
          dihtop(:)=max(dihtop(:),dihsol(:))
          dihbot(:)=min(dihbot(:),dihsol(:))
        endif
        if(jcomp(3))then
          read(7,*)dihexc(:)
          if(numacc.eq.1)then
            dihtep(:)=dihexc(:)
            dihbet(:)=dihexc(:)
          endif
          dihtep(:)=max(dihtep(:),dihexc(:))
          dihbet(:)=min(dihbet(:),dihexc(:))
        endif

c     read species inversion energies and find max and min

        if(lcomp(4))then
          read(7,*)invsol(:)
          if(numacc.eq.1)then
            invtop(:)=invsol(:)
            invbot(:)=invsol(:)
          endif
          invtop(:)=max(invtop(:),invsol(:))
          invbot(:)=min(invbot(:),invsol(:))
        endif
        if(jcomp(4))then
          read(7,*)invexc(:)
          if(numacc.eq.1)then
            invtep(:)=invexc(:)
            invbet(:)=invexc(:)
          endif
          invtep(:)=max(invtep(:),invexc(:))
          invbet(:)=min(invbet(:),invexc(:))
        endif

c     read species polarisation energies and find max and min

        if(lcomp(5))then
          read(7,*)shlsol(:)
          if(numacc.eq.1)then
            shltop(:)=shlsol(:)
            shlbot(:)=shlsol(:)
          endif
          shltop(:)=max(shltop(:),shlsol(:))
          shlbot(:)=min(shlbot(:),shlsol(:))
        endif
        if(jcomp(5))then
          read(7,*)shlexc(:)
          if(numacc.eq.1)then
            shltop(:)=shlexc(:)
            shlbot(:)=shlexc(:)
          endif
          shltep(:)=max(shltep(:),shlexc(:))
          shlbet(:)=min(shlbet(:),shlexc(:))
        endif

c     read species electrostatic energies and find max and min

        if(lcomp(6))then
          read(7,*)cousol(:)
          if(numacc.eq.1)then
            coutop(:)=cousol(:)
            coubot(:)=cousol(:)
          endif
          coutop(:)=max(coutop(:),cousol(:))
          coubot(:)=min(coubot(:),cousol(:))
        endif
        if(jcomp(6))then
          read(7,*)couexc(:)
          if(numacc.eq.1)then
            coutep(:)=couexc(:)
            coubet(:)=couexc(:)
          endif
          coutep(:)=max(coutep(:),couexc(:))
          coubet(:)=min(coubet(:),couexc(:))
        endif

c     read species van der waals energies and find max and min

        if(lcomp(7))then
          read(7,*)vdwsol(:)
          if(numacc.eq.1)then
            vdwtop(:)=vdwsol(:)
            vdwbot(:)=vdwsol(:)
          endif
          vdwtop(:)=max(vdwtop(:),vdwsol(:))
          vdwbot(:)=min(vdwbot(:),vdwsol(:))
        endif
        if(jcomp(7))then
          read(7,*)vdwexc(:)
          if(numacc.eq.1)then
            vdwtep(:)=vdwexc(:)
            vdwbet(:)=vdwexc(:)
          endif
          vdwtep(:)=max(vdwtep(:),vdwexc(:))
          vdwbet(:)=min(vdwbet(:),vdwexc(:))
        endif

c     read species 3 body energies and find max and min

        if(lcomp(8))then
          read(7,*)en3sol(:)
          if(numacc.eq.1)then
            en3top(:)=en3sol(:)
            en3bot(:)=en3sol(:)
          endif
          en3top(:)=max(en3top(:),en3sol(:))
          en3bot(:)=min(en3bot(:),en3sol(:))
        endif
        if(jcomp(8))then
          read(7,*)en3exc(:)
          if(numacc.eq.1)then
            en3tep(:)=en3exc(:)
            en3bet(:)=en3exc(:)
          endif
          en3tep(:)=max(en3tep(:),en3exc(:))
          en3bet(:)=min(en3bet(:),en3exc(:))
        endif

c     read species 4 body energies and find max and min

        if(lcomp(9))then
          read(7,*)en4sol(:)
          if(numacc.eq.1)then
            en4top(:)=en4sol(:)
            en4bot(:)=en4sol(:)
          endif
          en4top(:)=max(en4top(:),en4sol(:))
          en4bot(:)=min(en4bot(:),en4sol(:))
        endif
        if(jcomp(9))then
          read(7,*)en4exc(:)
          if(numacc.eq.1)then
            en4tep(:)=en4exc(:)
            en4bet(:)=en4exc(:)
          endif
          en4tep(:)=max(en4tep(:),en4exc(:))
          en4bet(:)=min(en4bet(:),en4exc(:))
        endif

      enddo

c     end of SOLVAT file

  100 close (7)

      write(*,'("SOLVAT file pass 1 completed")')
      write(*,'("Number of data points processed",i6)')numacc

c     reopen main input SOLVAT file and accumulate histograms

      open(7,file="SOLVAT")

      deld=1.d2/dble(numacc)

c     skip header data

      do i=1,5
        read(7,*)
      enddo

c     read periodic data

      do while(.true.)

c     skip first record

        read(7,*,end=200)

c     read species temperatures and calculate histogram

        read(7,*)tmpsol(:)
        do i=1,mxtmls
          delv=(tmptop(i)-tmpbot(i))/hdim
          if(delv.gt.0.d0)then
            k=int((tmpsol(i)-tmpbot(i))/delv+1.5d0)
            tmpdis(k,i)=tmpdis(k,i)+deld
          endif
        enddo

c     read species bond energies and calculate histogram

        if(lcomp(1))then
          read(7,*)bndsol(:)
          do i=1,mxtmls
            delv=(bndtop(i)-bndbot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((bndsol(i)-bndbot(i))/delv+1.5d0)
              bnddis(k,i)=bnddis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(1))then
          read(7,*)bndexc(:)
          do i=1,mxtmls
            delv=(bndtep(i)-bndbet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((bndexc(i)-bndbet(i))/delv+1.5d0)
              bnddes(k,i)=bnddes(k,i)+deld
            endif
          enddo
        endif

c     read species valence angle energies and calculate histogram

        if(lcomp(2))then
          read(7,*)angsol(:)
          do i=1,mxtmls
            delv=(angtop(i)-angbot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((angsol(i)-angbot(i))/delv+1.5d0)
              angdis(k,i)=angdis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(2))then
          read(7,*)angexc(:)
          do i=1,mxtmls
            delv=(angtep(i)-angbet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((angexc(i)-angbet(i))/delv+1.5d0)
              angdes(k,i)=angdes(k,i)+deld
            endif
          enddo
        endif

c     read species dihedral energies and calculate histogram

        if(lcomp(3))then
          read(7,*)dihsol(:)
          do i=1,mxtmls
            delv=(dihtop(i)-dihbot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((dihsol(i)-dihbot(i))/delv+1.5d0)
              dihdis(k,i)=dihdis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(3))then
          read(7,*)dihexc(:)
          do i=1,mxtmls
            delv=(dihtep(i)-dihbet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((dihexc(i)-dihbet(i))/delv+1.5d0)
              dihdes(k,i)=dihdes(k,i)+deld
            endif
          enddo
        endif

c     read species inversion energies and calculate histogram

        if(lcomp(4))then
          read(7,*)invsol(:)
          do i=1,mxtmls
            delv=(invtop(i)-invbot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((invsol(i)-invbot(i))/delv+1.5d0)
              invdis(k,i)=invdis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(4))then
          read(7,*)invexc(:)
          do i=1,mxtmls
            delv=(invtep(i)-invbet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((invexc(i)-invbet(i))/delv+1.5d0)
              invdes(k,i)=invdes(k,i)+deld
            endif
          enddo
        endif

c     read species polarisation energies and calculate histogram

        if(lcomp(5))then
          read(7,*)shlsol(:)
          do i=1,mxtmls
            delv=(shltop(i)-shlbot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((shlsol(i)-shlbot(i))/delv+1.5d0)
              shldis(k,i)=shldis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(5))then
          read(7,*)shlexc(:)
          do i=1,mxtmls
            delv=(shltep(i)-shlbet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((shlexc(i)-shlbet(i))/delv+1.5d0)
              shldes(k,i)=shldes(k,i)+deld
            endif
          enddo
        endif

c     read species electrostatic energies and calculate histogram

        if(lcomp(6))then
          read(7,*)cousol(:)
          do i=1,mxtmls2
            delv=(coutop(i)-coubot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((cousol(i)-coubot(i))/delv+1.5d0)
              coudis(k,i)=coudis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(6))then
          read(7,*)couexc(:)
          do i=1,mxtmls2
            delv=(coutep(i)-coubet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((couexc(i)-coubet(i))/delv+1.5d0)
              coudes(k,i)=coudes(k,i)+deld
            endif
          enddo
        endif

c     read species van der waals energies and calculate histogram

        if(lcomp(7))then
          read(7,*)vdwsol(:)
          do i=1,mxtmls2
            delv=(vdwtop(i)-vdwbot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((vdwsol(i)-vdwbot(i))/delv+1.5d0)
              vdwdis(k,i)=vdwdis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(7))then
          read(7,*)vdwexc(:)
          do i=1,mxtmls2
            delv=(vdwtep(i)-vdwbet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((vdwexc(i)-vdwbet(i))/delv+1.5d0)
              vdwdes(k,i)=vdwdes(k,i)+deld
            endif
          enddo
        endif

c     read species 3 body energies and calculate histogram

        if(lcomp(8))then
          read(7,*)en3sol(:)
          do i=1,mxtmls3
            delv=(en3top(i)-en3bot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((en3sol(i)-en3bot(i))/delv+1.5d0)
              en3dis(k,i)=en3dis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(8))then
          read(7,*)en3exc(:)
          do i=1,mxtmls3
            delv=(en3tep(i)-en3bet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((en3exc(i)-en3bet(i))/delv+1.5d0)
              en3des(k,i)=en3des(k,i)+deld
            endif
          enddo
        endif

c     read species 4 body energies and calculate histogram

        if(lcomp(9))then
          read(7,*)en4sol(:)
          do i=1,mxtmls4
            delv=(en4top(i)-en4bot(i))/hdim
            if(delv.gt.0.d0)then
              k=int((en4sol(i)-en4bot(i))/delv+1.5d0)
              en4dis(k,i)=en4dis(k,i)+deld
            endif
          enddo
        endif

        if(jcomp(9))then
          read(7,*)en4exc(:)
          do i=1,mxtmls4
            delv=(en4tep(i)-en4bet(i))/hdim
            if(delv.gt.0.d0)then
              k=int((en4exc(i)-en4bet(i))/delv+1.5d0)
              en4des(k,i)=en4des(k,i)+deld
            endif
          enddo
        endif

      enddo

c     end of SOLVAT file

  200 close (7)

      write(*,'("SOLVAT file pass 2 completed")')

c     open main output file SOLDENS

      open(8,file="SOLDENS")

      write(8,'(80a1)')header
      write(8,'(40a1)')units
      write(8,'("number of data points sampled",i10)')numacc

c     write species temperature distributions

      do i=1,mxtmls
        delv=(tmptop(i)-tmpbot(i))/hdim
        if(delv.gt.0.d0)then
          write(8,'("#temperature: species",i5)')i
          do j=1,dim
            write(8,'(1p,2e12.4)')dble(j-1)*delv+tmpbot(i),tmpdis(j,i)
          enddo
          write(8,'("&")')
        endif
      enddo

c     write bond energy distributions

      if(lcomp(1))then
        do i=1,mxtmls
          delv=(bndtop(i)-bndbot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#bond energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+bndbot(i),bnddis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(1))then
        do i=1,mxtmls
          delv=(bndtep(i)-bndbet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#bond energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+bndbet(i),bnddes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write valence angle energy distributions

      if(lcomp(2))then
        do i=1,mxtmls
          delv=(angtop(i)-angbot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#valence angle energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+angbot(i),angdis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(2))then
        do i=1,mxtmls
          delv=(angtep(i)-angbet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#valence angle energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+angbet(i),angdes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write dihedral energy distributions

      if(lcomp(3))then
        do i=1,mxtmls
          delv=(dihtop(i)-dihbot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#dihedral angle energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+dihbot(i),dihdis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(3))then
        do i=1,mxtmls
          delv=(dihtep(i)-dihbet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#dihedral angle energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+dihbet(i),dihdes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write inversion energy distributions

      if(lcomp(4))then
        do i=1,mxtmls
          delv=(invtop(i)-invbot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#inversion angle energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+invbot(i),invdis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(4))then
        do i=1,mxtmls
          delv=(invtep(i)-invbet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#inversion angle energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+invbet(i),invdes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write polarisation energy distributions

      if(lcomp(5))then
        do i=1,mxtmls
          delv=(shltop(i)-shlbot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#polarisation energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+shlbot(i),shldis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(5))then
        do i=1,mxtmls
          delv=(shltep(i)-shlbet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#polarisation energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+shlbet(i),shldes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write coulomb energy distributions

      if(lcomp(6))then
        do i=1,mxtmls2
          delv=(coutop(i)-coubot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#coulombic energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+coubot(i),coudis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(6))then
        do i=1,mxtmls2
          delv=(coutep(i)-coubet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#coulombic energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+coubet(i),coudes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write van der waals energy distributions

      if(lcomp(7))then
        do i=1,mxtmls2
          delv=(vdwtop(i)-vdwbot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#vdw energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+vdwbot(i),vdwdis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(7))then
        do i=1,mxtmls2
          delv=(vdwtep(i)-vdwbet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#vdw energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+vdwbet(i),vdwdes(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write 3 body energy dstributions

      if(lcomp(8))then
        do i=1,mxtmls3
          delv=(en3top(i)-en3bot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#3-body energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+en3bot(i),en3dis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(8))then
        do i=1,mxtmls3
          delv=(en3tep(i)-en3bet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#3-body energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+en3bet(i),en3des(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif

c     write 4 body energy distributions

      if(lcomp(9))then
        do i=1,mxtmls4
          delv=(en4top(i)-en4bot(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#4-body energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+en4bot(i),en4dis(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      if(jcomp(9))then
        do i=1,mxtmls4
          delv=(en4tep(i)-en4bet(i))/hdim
          if(delv.gt.0.d0)then
            write(8,'("#4-body energy: species",i5)')i
            do j=1,dim
              write(8,'(1p,2e14.6)')dble(j-1)*delv+en4bet(i),en4des(j,i)
            enddo
            write(8,'("&")')
          endif
        enddo
      endif
      write(8,'("END")')
      close(8)
      write(*,'("SOLDENS data file completed")')
      write(*,'("Job Done")')

      end
