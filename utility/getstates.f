      program getstates
c*********************************************************************
c
c     dl_poly utility program for establishing equivalence between
c     states generated by hyperdynamics
c
c     author w.smith september 2006
c     copyright daresbury laboratory
c
c*********************************************************************

      implicit none
      integer, parameter :: mxatms=1080
      logical same,flag
      character*4 tail
      character*8 name,file0,file1
      integer i,j,k,m,n,nlast,imcon,lev1,lev2,natms,nstates
      real*8 rcut,xx0,yy0,zz0,xx1,yy1,zz1,cell,cel0,cel1,rcel0,rcel1
      real*8 rct2,det,xxa,yya,zza,xxb,yyb,zzb,ssx,ssy,ssz,xss,yss,zss
      real*8 dxx,dyy,dzz,ddd

      dimension cel0(9),cel1(9),rcel0(9),rcel1(9),cell(9)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension xx1(mxatms),yy1(mxatms),zz1(mxatms)

      integer, allocatable :: sort(:),state(:)

      data cel0/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data cel1/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/

c     read instructions

      read(*,*)nlast,rcut
      rct2=rcut**2
      nstates=nlast+1

      allocate(sort(0:nlast),state(0:nlast))

c     initialise state array

      do i=0,nlast
        sort(i)=i
      enddo

c     outer loop over all files

      do m=0,nlast-1

c     read first file

        write(tail,'(i4.4)')m
        file0='BSN'//tail
        open(7,file=file0)
        read(7,*)
        read(7,*)lev1,imcon
        if(imcon.gt.0)then
          read(7,*)cel0(1),cel0(2),cel0(3)
          read(7,*)cel0(4),cel0(5),cel0(6)
          read(7,*)cel0(7),cel0(8),cel0(9)
        endif
        call invert(cel0,rcel0,det)
        do i=1,mxatms
          read(7,*,end=100)name
          read(7,*)xxa,yya,zza
          if(lev1.gt.0)read(7,*)
          if(lev1.gt.1)read(7,*)
          xx0(i)=rcel0(1)*xxa+rcel0(4)*yya+rcel0(7)*zza
          yy0(i)=rcel0(2)*xxa+rcel0(5)*yya+rcel0(8)*zza
          zz0(i)=rcel0(3)*xxa+rcel0(6)*yya+rcel0(9)*zza

        enddo
 100    continue
        close (7)
        natms=i-1

c     inner loop over all files

        do n=m+1,nlast

c     read second file

          write(tail,'(i4.4)')n
          file1='BSN'//tail
          open(8,file=file1)
          read(8,*)
          read(8,*)lev2,imcon
          if(imcon.gt.0)then
            read(8,*)cel1(1),cel1(2),cel1(3)
            read(8,*)cel1(4),cel1(5),cel1(6)
            read(8,*)cel1(7),cel1(8),cel1(9)
            call invert(cel1,rcel1,det)
          endif
          do i=1,natms
            read(8,*,end=200)name
            read(8,*)xxb,yyb,zzb
            if(lev2.gt.0)read(8,*)
            if(lev2.gt.1)read(8,*)
            xx1(i)=rcel1(1)*xxb+rcel1(4)*yyb+rcel1(7)*zzb
            yy1(i)=rcel1(2)*xxb+rcel1(5)*yyb+rcel1(8)*zzb
            zz1(i)=rcel1(3)*xxb+rcel1(6)*yyb+rcel1(9)*zzb

          enddo
 200      continue
          close(8)

c     average cell matrix of both structures

          do i=1,9
            cell(i)=0.5d0*(cel0(i)+cel1(i))
          enddo

c     compare structures

          ddd=0.d0
          do i=1,natms

            ssx=xx0(i)-xx1(i)
            ssy=yy0(i)-yy1(i)
            ssz=zz0(i)-zz1(i)

            xss=ssx-nint(ssx)
            yss=ssy-nint(ssy)
            zss=ssz-nint(ssz)

            dxx=cell(1)*xss+cell(4)*yss+cell(7)*zss
            dyy=cell(2)*xss+cell(5)*yss+cell(8)*zss
            dzz=cell(3)*xss+cell(6)*yss+cell(9)*zss
            ddd=max(ddd,dxx*dxx+dyy*dyy+dzz*dzz)

          enddo

          same=(ddd.le.rct2)
          if(same.and.(sort(m).ne.sort(n)))then
            j=min(sort(m),sort(n))
            k=max(sort(m),sort(n))
            sort(m)=j
            sort(n)=j
            nstates=nstates-1
            do i=0,nlast
              if(sort(i).eq.k)sort(i)=j
            enddo
          endif

        enddo
      enddo

c     sort states

      n=0
      state(0)=0
      do i=1,nlast
        if(sort(i).eq.sort(i-1))then
          state(i)=state(i-1)
        else
          j=0
          flag=.true.
          do while(flag.and.j.lt.i)
            if(sort(j).eq.sort(i))then
              flag=.false.
              state(i)=state(j)
            endif
            j=j+1
          enddo
          if(flag)then
            n=n+1
            state(i)=n
          endif
        endif
      enddo

c     write out state table

      write(*,*)'number of states ',nstates
      do i=0,nlast
        write(*,'(2i10)')i,state(i)
      enddo

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

