      program solsta

c*********************************************************************
c
c     dl_poly utility to calculate average solvation energies
c
c     copyright - daresbury laboratory
c     author    -   w.smith oct 2008
c
c*********************************************************************

      implicit none

      character header(80),units(40)
      logical lexcite,lswitch,lcomp(9),jcomp(9)
      integer nstep,natms,mxtmls,mxtmls2,mxtmls3,mxtmls4,numacc
      integer i,j,k,m,n
      real(8) tstep,elrc,elrc2,sclnv1,sclnv2
      real(8), allocatable :: tmpsol(:),tmpavg(:),tmpav2(:)
      real(8), allocatable :: bndsol(:),angsol(:),dihsol(:)
      real(8), allocatable :: invsol(:),shlsol(:),cousol(:)
      real(8), allocatable :: vdwsol(:),en3sol(:),en4sol(:)
      real(8), allocatable :: bndexc(:),angexc(:),dihexc(:)
      real(8), allocatable :: invexc(:),shlexc(:),couexc(:)
      real(8), allocatable :: vdwexc(:),en3exc(:),en4exc(:)
      real(8), allocatable :: bndavg(:),angavg(:),dihavg(:)
      real(8), allocatable :: invavg(:),shlavg(:),couavg(:)
      real(8), allocatable :: vdwavg(:),en3avg(:),en4avg(:)
      real(8), allocatable :: bndevg(:),angevg(:),dihevg(:)
      real(8), allocatable :: invevg(:),shlevg(:),couevg(:)
      real(8), allocatable :: vdwevg(:),en3evg(:),en4evg(:)
      real(8), allocatable :: bndav2(:),angav2(:),dihav2(:)
      real(8), allocatable :: invav2(:),shlav2(:),couav2(:)
      real(8), allocatable :: vdwav2(:),en3av2(:),en4av2(:)
      real(8), allocatable :: bndev2(:),angev2(:),dihev2(:)
      real(8), allocatable :: invev2(:),shlev2(:),couev2(:)
      real(8), allocatable :: vdwev2(:),en3ev2(:),en4ev2(:)

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

      allocate(tmpsol(mxtmls),tmpavg(mxtmls),tmpav2(mxtmls))
      if(lcomp(1))
     x  allocate(bndsol(mxtmls),bndavg(mxtmls),bndav2(mxtmls))
      if(lcomp(2))
     x  allocate(angsol(mxtmls),angavg(mxtmls),angav2(mxtmls))
      if(lcomp(3))
     x  allocate(dihsol(mxtmls),dihavg(mxtmls),dihav2(mxtmls))
      if(lcomp(4))
     x  allocate(invsol(mxtmls),invavg(mxtmls),invav2(mxtmls))
      if(lcomp(5))
     x  allocate(shlsol(mxtmls),shlavg(mxtmls),shlav2(mxtmls))
      if(lcomp(6))
     x  allocate(cousol(mxtmls2),couavg(mxtmls2),couav2(mxtmls2))
      if(lcomp(7))
     x  allocate(vdwsol(mxtmls2),vdwavg(mxtmls2),vdwav2(mxtmls2))
      if(lcomp(8))
     x  allocate(en3sol(mxtmls3),en3avg(mxtmls3),en3av2(mxtmls3))
      if(lcomp(9))
     x  allocate(en4sol(mxtmls4),en4avg(mxtmls4),en4av2(mxtmls4))
      if(jcomp(1))
     x  allocate(bndexc(mxtmls),bndevg(mxtmls),bndev2(mxtmls))
      if(jcomp(2))
     x  allocate(angexc(mxtmls),angevg(mxtmls),angev2(mxtmls))
      if(jcomp(3))
     x  allocate(dihexc(mxtmls),dihevg(mxtmls),dihev2(mxtmls))
      if(jcomp(4))
     x  allocate(invexc(mxtmls),invevg(mxtmls),invev2(mxtmls))
      if(jcomp(5))
     x  allocate(shlexc(mxtmls),shlevg(mxtmls),shlev2(mxtmls))
      if(jcomp(6))
     x  allocate(couexc(mxtmls2),couevg(mxtmls2),couev2(mxtmls2))
      if(jcomp(7))
     x  allocate(vdwexc(mxtmls2),vdwevg(mxtmls2),vdwev2(mxtmls2))
      if(jcomp(8))
     x  allocate(en3exc(mxtmls3),en3evg(mxtmls3),en3ev2(mxtmls3))
      if(jcomp(9))
     x  allocate(en4exc(mxtmls4),en4evg(mxtmls4),en4ev2(mxtmls4))

c     set statistical variables

      numacc=0
      tmpavg(:)=0.d0
      if(lcomp(1))bndavg(:)=0.d0
      if(jcomp(1))bndevg(:)=0.d0
      if(lcomp(2))angavg(:)=0.d0
      if(jcomp(2))angevg(:)=0.d0
      if(lcomp(3))dihavg(:)=0.d0
      if(jcomp(3))dihevg(:)=0.d0
      if(lcomp(4))invavg(:)=0.d0
      if(jcomp(4))invevg(:)=0.d0
      if(lcomp(5))shlavg(:)=0.d0
      if(jcomp(5))shlevg(:)=0.d0
      if(lcomp(6))couavg(:)=0.d0
      if(jcomp(6))couevg(:)=0.d0
      if(lcomp(7))vdwavg(:)=0.d0
      if(jcomp(7))vdwevg(:)=0.d0
      if(lcomp(8))en3avg(:)=0.d0
      if(jcomp(8))en3evg(:)=0.d0
      if(lcomp(9))en4avg(:)=0.d0
      if(jcomp(9))en4evg(:)=0.d0

c     read periodic data

      do while(.true.)

        if(lexcite)then
          read(7,'(8x,i10,f12.5,2e14.6)',end=100)nstep,tstep,elrc,elrc2
        else
          read(7,'(8x,i10,f12.5,e14.6)',end=100)nstep,tstep,elrc
        endif

c     statistical counters

        numacc=numacc+1
        sclnv2=1.d0/dble(numacc)
        sclnv1=dble(numacc-1)/dble(numacc)

c     read species temperatures and calculate averages

        read(7,*)tmpsol(:)
        tmpav2(:)=sclnv1*(tmpav2(:)+sclnv2*(tmpsol(:)-tmpavg(:))**2)
        tmpavg(:)=sclnv1*tmpavg(:)+sclnv2*tmpsol(:)

c     read species bond energies calculate averages

        if(lcomp(1))then
          read(7,*)bndsol(:)
          bndav2(:)=sclnv1*(bndav2(:)+sclnv2*(bndsol(:)-bndavg(:))**2)
          bndavg(:)=sclnv1*bndavg(:)+sclnv2*bndsol(:)
        endif
        if(jcomp(1))then
          read(7,*)bndexc(:)
          bndev2(:)=sclnv1*(bndev2(:)+sclnv2*(bndexc(:)-bndevg(:))**2)
          bndevg(:)=sclnv1*bndevg(:)+sclnv2*bndexc(:)
        endif

c     read species valence angle energies calculate averages

        if(lcomp(2))then
          read(7,*)angsol(:)
          angav2(:)=sclnv1*(angav2(:)+sclnv2*(angsol(:)-angavg(:))**2)
          angavg(:)=sclnv1*angavg(:)+sclnv2*angsol(:)
        endif
        if(jcomp(2))then
          read(7,*)angexc(:)
          angev2(:)=sclnv1*(angev2(:)+sclnv2*(angexc(:)-angevg(:))**2)
          angevg(:)=sclnv1*angevg(:)+sclnv2*angexc(:)
        endif

c     read species dihedral energies calculate averages

        if(lcomp(3))then
          read(7,*)dihsol(:)
          dihav2(:)=sclnv1*(dihav2(:)+sclnv2*(dihsol(:)-dihavg(:))**2)
          dihavg(:)=sclnv1*dihavg(:)+sclnv2*dihsol(:)
        endif
        if(jcomp(3))then
          read(7,*)dihexc(:)
          dihev2(:)=sclnv1*(dihev2(:)+sclnv2*(dihexc(:)-dihevg(:))**2)
          dihevg(:)=sclnv1*dihevg(:)+sclnv2*dihexc(:)
        endif

c     read species inversion energies calculate averages

        if(lcomp(4))then
          read(7,*)invsol(:)
          invav2(:)=sclnv1*(invav2(:)+sclnv2*(invsol(:)-invavg(:))**2)
          invavg(:)=sclnv1*invavg(:)+sclnv2*invsol(:)
        endif
        if(jcomp(4))then
          read(7,*)invexc(:)
          invev2(:)=sclnv1*(invev2(:)+sclnv2*(invexc(:)-invevg(:))**2)
          invevg(:)=sclnv1*invevg(:)+sclnv2*invexc(:)
        endif

c     read species polarisation energies calculate averages

        if(lcomp(5))then
          read(7,*)shlsol(:)
          shlav2(:)=sclnv1*(shlav2(:)+sclnv2*(shlsol(:)-shlavg(:))**2)
          shlavg(:)=sclnv1*shlavg(:)+sclnv2*shlsol(:)
        endif
        if(jcomp(5))then
          read(7,*)shlexc(:)
          shlev2(:)=sclnv1*(shlev2(:)+sclnv2*(shlexc(:)-shlevg(:))**2)
          shlevg(:)=sclnv1*shlevg(:)+sclnv2*shlexc(:)
        endif

c     read species electrostatic energies calculate averages

        if(lcomp(6))then
          read(7,*)cousol(:)
          couav2(:)=sclnv1*(couav2(:)+sclnv2*(cousol(:)-couavg(:))**2)
          couavg(:)=sclnv1*couavg(:)+sclnv2*cousol(:)
        endif
        if(jcomp(6))then
          read(7,*)couexc(:)
          couev2(:)=sclnv1*(couev2(:)+sclnv2*(couexc(:)-couevg(:))**2)
          couevg(:)=sclnv1*couevg(:)+sclnv2*couexc(:)
        endif

c     read species van der waals energies calculate averages

        if(lcomp(7))then
          read(7,*)vdwsol(:)
          vdwav2(:)=sclnv1*(vdwav2(:)+sclnv2*(vdwsol(:)-vdwavg(:))**2)
          vdwavg(:)=sclnv1*vdwavg(:)+sclnv2*vdwsol(:)
        endif
        if(jcomp(7))then
          read(7,*)vdwexc(:)
          vdwev2(:)=sclnv1*(vdwev2(:)+sclnv2*(vdwexc(:)-vdwevg(:))**2)
          vdwevg(:)=sclnv1*vdwevg(:)+sclnv2*vdwexc(:)
        endif

c     read species 3 body energies calculate averages

        if(lcomp(8))then
          read(7,*)en3sol(:)
          en3av2(:)=sclnv1*(en3av2(:)+sclnv2*(en3sol(:)-en3avg(:))**2)
          en3avg(:)=sclnv1*en3avg(:)+sclnv2*en3sol(:)
        endif
        if(jcomp(8))then
          read(7,*)en3exc(:)
          en3ev2(:)=sclnv1*(en3ev2(:)+sclnv2*(en3exc(:)-en3evg(:))**2)
          en3evg(:)=sclnv1*en3evg(:)+sclnv2*en3exc(:)
        endif

c     read species 4 body energies calculate averages

        if(lcomp(9))then
          read(7,*)en4sol(:)
          en4av2(:)=sclnv1*(en4av2(:)+sclnv2*(en4exc(:)-en4avg(:))**2)
          en4avg(:)=sclnv1*en4avg(:)+sclnv2*en4exc(:)
        endif
        if(jcomp(9))then
          read(7,*)en4exc(:)
          en4ev2(:)=sclnv1*(en4ev2(:)+sclnv2*(en4sol(:)-en4evg(:))**2)
          en4evg(:)=sclnv1*en4evg(:)+sclnv2*en4sol(:)
        endif

      enddo

c     end of SOLVAT file

  100 close (7)
      write(*,'("Finished reading SOLVAT file")')
      write(*,'("number of data points sampled",i10)')numacc

c     open main output file SOLSTATS

      open(8,file="SOLSTATS")

      write(8,'(80a1)')header
      write(8,'(40a1)')units
      write(8,'("number of data points sampled",i10)')numacc

c     write average species temperatures and rms deviations

      write(8,'("species temperatures and rms deviations")')
      do i=1,mxtmls
        write(8,'(i6,6x,1p,2e14.6)')i,tmpavg(i),sqrt(tmpav2(i))
      enddo

c     write average bond energies and rms deviations

      if(lcomp(1))then
        write(8,'("species bond energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,bndavg(i),sqrt(bndav2(i))
        enddo
      endif
      if(jcomp(1))then
        write(8,'("species bond energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,bndevg(i),sqrt(bndev2(i))
        enddo
      endif

c     write average valence angle energies and rms deviations

      if(lcomp(2))then
        write(8,'("species valence angle energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,angavg(i),sqrt(angav2(i))
        enddo
      endif
      if(jcomp(2))then
        write(8,'("species valence angle energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,angevg(i),sqrt(angev2(i))
        enddo
      endif

c     write average dihedral energies and rms deviations

      if(lcomp(3))then
        write(8,
     x    '("species dihedral angle energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,dihavg(i),sqrt(dihav2(i))
        enddo
      endif
      if(jcomp(3))then
        write(8,
     x    '("species dihedral angle energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,dihevg(i),sqrt(dihev2(i))
        enddo
      endif

c     write average inversion energies and rms deviations

      if(lcomp(4))then
        write(8,
     x    '("species inversion angle energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,invavg(i),sqrt(invav2(i))
        enddo
      endif
      if(jcomp(4))then
        write(8,
     x    '("species inversion angle energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,invevg(i),sqrt(invev2(i))
        enddo
      endif

c     write average polarisation energies and rms deviations

      if(lcomp(5))then
        write(8,'("species polarisation energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,shlavg(i),sqrt(shlav2(i))
        enddo
      endif
      if(jcomp(5))then
        write(8,'("species polarisation energies and rms deviations")')
        do i=1,mxtmls
          write(8,'(i6,6x,1p,2e14.6)')i,shlevg(i),sqrt(shlev2(i))
        enddo
      endif

c     write average coulomb energies and rms deviations

      if(lcomp(6))then
        write(8,'("inter-species coulomb energies and rms deviations")')
        k=0
        do i=1,mxtmls
          do j=1,i
            k=k+1
            write(8,'(2i6,6x,1p,2e14.6)')i,j,couavg(k),sqrt(couav2(k))
          enddo
        enddo
      endif
      if(jcomp(6))then
        write(8,'("inter-species coulomb energies and rms deviations")')
        k=0
        do i=1,mxtmls
          do j=1,i
            k=k+1
            write(8,'(2i6,6x,1p,2e14.6)')i,j,couevg(k),sqrt(couev2(k))
          enddo
        enddo
      endif

c     write average van der waals energies and rms deviations

      if(lcomp(7))then
        write(8,
     x    "('inter-species van der waals energies and rms deviations')")
        k=0
        do i=1,mxtmls
          do j=1,i
            k=k+1
            write(8,'(2i6,6x,1p,2e14.6)')i,j,vdwavg(k),sqrt(vdwav2(k))
          enddo
        enddo
      endif
      if(jcomp(7))then
        write(8,
     x    "('inter-species van der waals energies and rms deviations')")
        k=0
        do i=1,mxtmls
          do j=1,i
            k=k+1
            write(8,'(2i6,6x,1p,2e14.6)')i,j,vdwevg(k),sqrt(vdwev2(k))
          enddo
        enddo
      endif

c     write average 3 body energies and rms deviations

      if(lcomp(8))then
        write(8,'("inter-species 3-body energies and rms deviations")')
        m=0
        do i=1,mxtmls
          do j=1,i
            do k=1,j
              m=m+1
              write(8,'(3i6,6x,1p,2e14.6)')i,j,k,en3avg(m),
     x          sqrt(en3av2(m))
            enddo
          enddo
        enddo
      endif
      if(jcomp(8))then
        write(8,'("inter-species 3-body energies and rms deviations")')
        m=0
        do i=1,mxtmls
          do j=1,i
            do k=1,j
              m=m+1
              write(8,'(3i6,6x,1p,2e14.6)')i,j,k,en3evg(m),
     x          sqrt(en3ev2(m))
            enddo
          enddo
        enddo
      endif

c     write average 4 body energies and rms deviations

      if(lcomp(9))then
        write(8,'("inter-species 4-body energies and rms deviations")')
        n=0
        do i=1,mxtmls
          do j=1,i
            do k=1,j
              do m=1,k
                n=n+1
                write(8,'(4i6,6x,1p,2e14.6)')i,j,k,m,en4avg(n),
     x            sqrt(en4av2(n))
              enddo
            enddo
          enddo
        enddo
      endif
      if(jcomp(9))then
        write(8,'("inter-species 4-body energies and rms deviations")')
        n=0
        do i=1,mxtmls
          do j=1,i
            do k=1,j
              do m=1,k
                n=n+1
                write(8,'(4i6,6x,1p,2e14.6)')i,j,k,m,en4evg(n),
     x            sqrt(en4ev2(n))
              enddo
            enddo
          enddo
        enddo
      endif
      write(8,'("END")')

c     close SOLSTATS file

      close(8)
      write(*,'("Finished writing SOLSTATS file")')
      write(*,'("Job done")')

      end
