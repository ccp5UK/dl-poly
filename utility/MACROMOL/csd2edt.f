      program csd2edt

c***********************************************************************
c
c     DL_POLY utility to convert fractional coordinates into real space
c     coordinates in form suitable for EDTOUT (see ffgen.f). The
c     FF.dat file is also created.
c
c     input a,b,c (cell lengths)
c     input alp,bet,gam  (cell angles)
c     input xi,yi,zi - fractional coordinates
c
c     copyright darsebury laboratory 1993
c     author -  t. forester     june 1993
c
c***********************************************************************

      implicit real*8 (a-h,o-z)
      parameter (mxn= 10000)
      dimension cell(9),xb(mxn),yb(mxn),zb(mxn),q(mxn)
      dimension xx(mxn),yy(mxn),zz(mxn)
      dimension icn(mxn,4),ict(mxn)
      character*6 name(mxn),ambnam(mxn),nami,namj
      character*40 title
      logical print
c
c     read title
      read(*,'(a40,3f8.0)')title,a0,b0,c0
c
c     read cell angles
      read(*,'(22x,3f8.0)') alp,bet,gam

c
c     read number of sites in unit cell
      read(*,*)natms
      read(*,*)
c
c     create cell vectors

      if(a0.le.0.or.b0.le.0.or.c0.le.0) then

         write(6,*) 'error : zero vector length : ',a0,b0,c0
         call exit()

      endif

      pi = 3.141592653589793d0
      alp = alp *pi/180.d0
      bet = bet *pi/180.d0
      gam = gam *pi/180.d0
c
c     a direction is along x

      cell(1) = a0
      cell(2) = 0.d0
      cell(3) = 0.d0

c
c     b direction is in x-y plane

      cell(4) = b0*cos(gam)
      cell(5) = b0*sin(gam)
      cell(6) = 0.d0

c
c     c direction

      cell(7) = c0*cos(bet)
      cell(8) = c0*sin(bet)*cos(alp)
      cell(9) = c0*sin(bet)*sin(alp)
c
c     remove rounding errors for small entries

      do i = 1,9

         if(abs(cell(i)).lt.1d-10) cell(i) = 0.d0

      enddo
c
c     write out cell vectors, and angles

      write(*,*) 'cell vectors are'
      write(*,'(3g14.5)') cell(1),cell(2),cell(3)
      write(*,'(3g14.5)') cell(4),cell(5),cell(6)
      write(*,'(3g14.5)') cell(7),cell(8),cell(9)

c
c     calculate angles

      ga = cell(1)*cell(4) + cell(2)*cell(5) + cell(3)*cell(6)
      ga = acos(ga /a0/b0)*180.d0/pi

      be = cell(1)*cell(7) + cell(2)*cell(8) + cell(3)*cell(9)
      be = acos(be /a0/c0)*180.d0/pi

      al = cell(7)*cell(4) + cell(8)*cell(5) + cell(9)*cell(6)
      al = acos(al /c0/b0)*180.d0/pi

      write(*,*)' cell angles ',al,be,ga

      do i=1,natms
         read(5,'(5x,a5,3f10.5,1x,4i4)') name(i),xb(i),yb(i),zb(i),
     x        (icn(i,k),k=1,4)

      enddo

      open(10,file='EDTOUT')
      write(10,'(a80)')title
      write(10,'(a)')'CELL_VECTORS'
      write(10,'(3f20.10)') cell
      write(10,'(/,6x,a,/,5x,a,/) 'MOLECULE','RESIDUE = ....'
      write(10,'(5x,a./)')'BOND ARRAY'
c
c     fix names and connections

      do i = 1,natms
         k = 0
         do j = 1,4
            if(icn(i,j).gt.0) k = k+1
         enddo
         ict(i) = k
      enddo
c
c     assign amber names and charges (unique for valinomycin molecule!)

      do i = 1,natms
         ambnam(i) = name(i)
         q(i) = 0.d0

         do j = ict(i)+1,4
            icn(i,j) = -99
         enddo

         nami = name(i)
         if(nami(1:1).eq.'N') then
            ambnam(i) = 'N     '

         elseif(nami(1:1).eq.'H') then
            ambnam(i) = 'H     '

         elseif (nami(1:1).eq.'C') then
            ambnam(i) = 'C     '

         elseif (nami(1:1).eq.'O') then
            ambnam(i) = 'O     '

         endif
      enddo

      do i = 1,natms
         nami = ambnam(i)
         namj = ambnam(icn(i,1))

         if (nami(1:1).eq.'H'.and.namj(1:1).eq.'C') then
            ambnam(i)='HC    '

         elseif (nami(1:1).eq.'C'.and.ict(i).eq.4) then
            ambnam(i) = 'CT    '

         elseif (nami(1:1).eq.'O'.and.ict(i).eq.2) then
            ambnam(i) = 'OS    '

         endif

      enddo
c
c     refinement of labels C labels

c     if C not bonded to O => is CT
      do i = 1,natms
         if(ambnam(i).eq.'C     ') then
            ik = 0
            do j = 1,ict(i)
               nami = ambnam(icn(i,j))
               if(nami(1:1).eq.'O') ik = ik+ 1
               if(nami(1:1).eq.'N') ik = ik+ 1
            enddo
            if(ik.le.1) then
               ambnam(i)='CT    '
               iadd(i) = 4-ict(i)
            endif
         endif
      enddo
c
c     set up lattice

      do  i=1,natms

         xx(i)=cell(1)*xb(i)+cell(4)*yb(i)+cell(7)*zb(i)
         yy(i)=cell(2)*xb(i)+cell(5)*yb(i)+cell(8)*zb(i)
         zz(i)=cell(3)*xb(i)+cell(6)*yb(i)+cell(9)*zb(i)

         print = .true.
         do j = 1,4

            if(icn(i,j).gt.i.or.icn(i,j).lt.0) then
               icn(i,j) = -99
               if(j.gt.1) print = .false.
            endif
            if(print) write(10,'(2i5,3x,2a6,4x,3f10.4,f8.3)')
     x           i,icn(i,j),name(i),ambnam(i),xx(i),yy(i),zz(i),q(i)
            if(icn(i,j).lt.0) print = .false.

         enddo

      enddo

c
c     write FF.dat file

      open (11, file='FF.dat')
      write(11,'(a80)') title
      write(11,*)'CONSTRAIN none                     '
      write(11,*) 'UNITS kj'

      stop
      end




