      program proseq
c
c***********************************************************************
c
c     dl_poly utility program for reading SEQNET data files and
c     extracting the protein structure in a format compatible
c     with the dl_poly CONFIG file
c
c     copyright - daresbury laboratory 1993
c     author    - w. smith jan 1993
c
c***********************************************************************
c

      implicit double precision (a-h,o-z)

      parameter (mxatm = 100000, mxseq = 5000)
      character*6 label,header,atom,hetatm,seqres
      character*60 title
      character*80 record
      character*4 name,tseq(13),seq(mxseq)

      logical sflag
      data header/'HEADER'/,atom/'ATOM'/,hetatm/'HETATM'/
     $     ,seqres/'SEQRES'/

      sflag = .true.
      kseq = 0
      do i=1,mxatm
         read(*,'(a80)',end=1000)record
         label=record
         if(label.eq.header)then

            title=record(11:70)
            write(*,'(a60)')title

         else if(label.eq.atom.or.label.eq.hetatm)then

            read(record,'(12x,a4,14x,3f8.3)')name,x,y,z
            do k=1,4
               if(name(1:1).eq." ")then
                  do j=1,3
                     name(j:j)=name(j+1:j+1)
                  enddo
                  name(4:4)=" "
               endif
            enddo
            if(label.eq.hetatm.and.name.eq.'O   ') name = 'OW  '
            write(*,'(a4)')name
            write(*,'(3f20.8)')x,y,z

         endif

      enddo

 1000 stop
      end
