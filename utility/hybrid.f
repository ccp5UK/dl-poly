      program hybrid
c
c**********************************************************************
c
c     dl_poly utility to determine the probable hybridization of atoms
c     in an organic molecule based on a bondlength criterion (warning -
c     not foolproof, but a useful start)
c
c     molecular configuration specified in dl_poly CONFIG file
c
c     for use with Dreiding forcefield
c
c     copyright daresbury laboratory 1999
c     author w. smith november 1999
c
c**********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (mxatm=2500,mega=1000000,mxpair=15)
      parameter (pi=3.141592653589793d0)
      character*8 name(mxatm)
      character*80 title
      character*40 fname
      dimension xxx(mxatm),yyy(mxatm),zzz(mxatm)
      dimension list(mxatm,mxpair),num(mxatm)
      dimension cell(9),rcell(9),celprp(10)
c
c     read control parameters

      write(*,'(a)')'enter the maximum permitted bondlength:'
      read(*,*)bond
      write(*,'(a)')'enter name of config file'
      read(*,'(a40)')fname
c
c     open files

      open(7,file=fname,err=200)
      open(8,file='CFGHYB')
c
c     read configuration data
      read(7,'(a)')title
      write(*,'(a)')'File header: ',title
      write(8,'(a)')title
      read(7,*)levcfg,imcon
      write(8,'(2i10)')0,imcon
      if(imcon.gt.0)then
        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)
        write(8,'(3f20.10)')cell(1),cell(2),cell(3)
        write(8,'(3f20.10)')cell(4),cell(5),cell(6)
        write(8,'(3f20.10)')cell(7),cell(8),cell(9)
        call invert(cell,rcell,det)
        call dcell(cell,celprp)
        rchk=0.5d0*min(celprp(7),celprp(8),celprp(9))
        if(rchk.lt.bond)then
          write(*,'(a)')'error - chosen bondlength exceeds cell width'
          stop
        endif
      endif
      do i=1,mega
        if(i.gt.mxatm)then
          write(*,'(a)')'error - too many atoms in config file'
          stop
        endif
        read(7,'(a)',end=100)name(i)
        read(7,*,end=100)xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)read(7,*,end=100)vxx,vyy,vzz
        if(levcfg.gt.1)read(7,*,end=100)fxx,fyy,fzz
      enddo
  100 natm=i-1
      write(*,'(a,i6)')'number of atoms in config file: ',natm
      close (7)
c
c     initialize num array

      do i=1,natm

        num(i)=0

      enddo
c
c     determine bonds in system
      k=0
      last=natm
      m1=natm/2
      m2=(natm-1)/2
      do m=1,m1
        if(m.gt.m2)last=m1
        do i=1,last
          j=i+m
          if(j.gt.natm)j=j-natm
          xd=xxx(i)-xxx(j)
          yd=yyy(i)-yyy(j)
          zd=zzz(i)-zzz(j)
          if(imcon.gt.0)then
            ssx=(rcell(1)*xd+rcell(4)*yd+rcell(7)*zd)
            ssy=(rcell(2)*xd+rcell(5)*yd+rcell(8)*zd)
            ssz=(rcell(3)*xd+rcell(6)*yd+rcell(9)*zd)
            xss=ssx-nint(ssx)
            yss=ssy-nint(ssy)
            zss=ssz-nint(ssz)
            xd=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
            yd=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
            zd=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          endif
          rsq=xd*xd+yd*yd+zd*zd
          if(rsq.le.bond**2)then

            if(num(i).eq.mxpair.or.num(j).eq.mxpair)then
              write(*,'(a)')'error - too many bonds found'
              stop
            endif
            k=k+1
            num(i)=num(i)+1
            num(j)=num(j)+1
            list(i,num(i))=j
            list(j,num(j))=i

          endif
        enddo
      enddo
      npairs=k
      write(*,'(a,i6)')'number of bonds in config file: ',npairs
c
c     determine probable hybridizations

      do i=1,natm

        if(name(i)(1:1).eq."H")then

          name(i)="H_"

        elseif(name(i)(1:1).eq."P")then

          name(i)="P_3"

        elseif(name(i)(1:1).eq."C")then

          if(num(i).ge.4)then

            name(i)="C_3"

          elseif(num(i).eq.3)then

            name(i)="C_2"

          elseif(num(i).eq.2)then

            name(i)="C_1"

          endif

        elseif(name(i)(1:1).eq."O")then

          if(num(i).ge.2)then

            name(i)="O_3"

          elseif(num(i).eq.1)then

            name(i)="O_2"

          endif

        elseif(name(i)(1:1).eq."N")then

          if(num(i).ge.3)then

            name(i)="N_3"

          elseif(num(i).eq.2)then

            name(i)="N_2"

          elseif(num(i).eq.1)then

            name(i)="N_1"

          endif
        endif
      enddo
c
c     now check for possible resonance condition

      do i=1,natm

        if(name(i).eq."O_3".or.name(i).eq."N_3")then

          do j=1,num(i)

            if(name(j).eq."C_2")name(i)=name(i)(1:1)//"_2"

          enddo

        endif

      enddo
c
c     write out new config file

      do i=1,natm

        write(8,'(a8,i10)')name(i),i
        write(8,'(3f20.10)')xxx(i),yyy(i),zzz(i)

      enddo
      close (8)
      stop

  200 continue
      write(*,*)'Error - file not found'
      stop

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
