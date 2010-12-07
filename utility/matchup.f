c*********************************************************************
c
c     molecular matching program
c
c     copyright daresbury laboratory
c     author w.smith march 2004
c
c*********************************************************************
      
      implicit none

      integer, parameter :: npass=10
      integer, parameter :: mxatms=10000

      character*40 file0,file1
      integer i,pass,natms,natm0,natm1,imcon0,imcon1
      real*8 cell0,cell1,quality,compare,fitold,fitnew,offset

      real*8 xx0(mxatms),yy0(mxatms),zz0(mxatms)
      real*8 xx1(mxatms),yy1(mxatms),zz1(mxatms)

      dimension cell0(9),cell1(9)

c     read template stucture

      write(*,*)'Enter name of template CONFIG file'
      read(*,'(a40)')file0
      call cfgread(file0,natm0,imcon0,cell0,xx0,yy0,zz0)

c     read trial stucture

      write(*,*)'Enter name of trial CONFIG file'
      read(*,'(a40)')file1
      call cfgread(file1,natm1,imcon1,cell1,xx1,yy1,zz1)

      if(natm0.ne.natm1)then

        write(*,*)'Error - CONFIG files not equivalent'
        stop

      endif

      natms=natm0
      write(*,*)'Number of atoms in template',natms

c     centre the template

      call centre(natms,xx0,yy0,zz0)

c     optimise match

      quality=1.d-6
      compare=1.d20
      fitnew=compare

      pass=0
      do while(compare.gt.quality.and.pass.lt.npass)

        pass=pass+1

c     slide positions to best match of template

        call slide(natms,xx0,yy0,zz0,xx1,yy1,zz1)

c     rotate positions to best match of template

        call rotate(natms,xx0,yy0,zz0,xx1,yy1,zz1)

c     calculate current mismatch with template

        fitold=fitnew
        fitnew=offset(natms,xx0,yy0,zz0,xx1,yy1,zz1)
        compare=abs(fitnew-fitold)

      enddo

c     print out best fit

      write(*,*)'RMS structural fit (Angstroms) = ',fitnew

      end

      subroutine cfgread(fname,natms,imcon,cell,xxx,yyy,zzz)
      
c*********************************************************************
c
c     read a DL_POLY CONFIG file
c
c     copyright daresbury laboratory
c     author w.smith march 2004
c
c*********************************************************************

      implicit none

      character*8 name
      character*80 head
      character*40 fname

      integer i,imcon,levcfg,natms

      real*8 a,b,c,cell,xxx,yyy,zzz
      dimension cell(*),xxx(*),yyy(*),zzz(*)

c     open CONFIG file

      open(7,file=fname)

c     read file header

      read(7,'(a80)')head
      write(*,*)'File header: ',head
      read(7,'(2i10)')levcfg,imcon

c     set default cell vectors

      do i=1,9
        cell(i)=0.d0
      enddo

c     read cell vectors

      if(imcon.gt.0)then

        read(7,*)cell(1),cell(2),cell(3)
        read(7,*)cell(4),cell(5),cell(6)
        read(7,*)cell(7),cell(8),cell(9)

      endif

      i=0
      do while(.true.)

        read(7,'(a8)',end=100)name
        i=i+1
        read(7,*)xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)read(7,*)a,b,c
        if(levcfg.gt.1)read(7,*)a,b,c

      enddo

  100 natms=i

c     close CONFIG file

      close(7)

      return
      end

      subroutine centre(natms,xx0,yy0,zz0)
      
c*********************************************************************
c     
c     subroutine to centre the template
c     
c     copyright daresbury laboratory
c     author w.smith march 2004
c     
c*********************************************************************
      
      implicit none
      
      integer i,natms
      real*8 rrr,xx0,yy0,zz0
      
      dimension xx0(*),yy0(*),zz0(*)
      dimension rrr(3)
      
      rrr(1)=0.d0
      rrr(2)=0.d0
      rrr(3)=0.d0
      
      do i=1,natms
        
        rrr(1)=rrr(1)+xx0(i)
        rrr(2)=rrr(2)+yy0(i)
        rrr(3)=rrr(3)+zz0(i)
        
      enddo
      
      rrr(1)=rrr(1)/dble(natms)
      rrr(2)=rrr(2)/dble(natms)
      rrr(3)=rrr(3)/dble(natms)
      
      do i=1,natms
        
        xx0(i)=xx0(i)-rrr(1)
        yy0(i)=yy0(i)-rrr(2)
        zz0(i)=zz0(i)-rrr(3)
        
      enddo
      
      return
      end

      subroutine slide(natms,xx0,yy0,zz0,xx1,yy1,zz1)
      
c*********************************************************************
c     
c     slide structure to best match of template
c     
c     copyright daresbury laboratory
c     author w.smith march 2004
c     
c*********************************************************************
      
      implicit none
      
      integer i,natms
      real*8 rrr,xx0,yy0,zz0,xx1,yy1,zz1
      
      dimension xx0(*),yy0(*),zz0(*),xx1(*),yy1(*),zz1(*)
      dimension rrr(3)
      
      rrr(1)=0.d0
      rrr(2)=0.d0
      rrr(3)=0.d0
      
      do i=1,natms
        
        rrr(1)=rrr(1)+(xx1(i)-xx0(i))
        rrr(2)=rrr(2)+(yy1(i)-yy0(i))
        rrr(3)=rrr(3)+(zz1(i)-zz0(i))
        
      enddo
      
      rrr(1)=rrr(1)/dble(natms)
      rrr(2)=rrr(2)/dble(natms)
      rrr(3)=rrr(3)/dble(natms)
      
      do i=1,natms
        
        xx1(i)=xx1(i)-rrr(1)
        yy1(i)=yy1(i)-rrr(2)
        zz1(i)=zz1(i)-rrr(3)
        
      enddo
      
      return
      end

      subroutine rotate(natms,xx0,yy0,zz0,xx1,yy1,zz1)
      
c*********************************************************************
c     
c     subroutine to centre the template
c     
c     copyright daresbury laboratory
c     author w.smith march 2004
c     
c*********************************************************************
      
      implicit none
      
      integer i,j,k,natms
      real*8 xx0,yy0,zz0,xx1,yy1,zz1,mat,vec,aaa,qqq,rot,txx,tyy,tzz
      
      dimension xx0(*),yy0(*),zz0(*),xx1(*),yy1(*),zz1(*)
      dimension aaa(3,3),mat(4,4),vec(4,4),qqq(4),rot(3,3)
      
c     zero work arrays

      do i=1,3
        do j=1,3
          aaa(i,j)=0.d0
        enddo
      enddo

c     calculate optimisation parameters

      do i=1,natms
        
        aaa(1,1)=aaa(1,1)+xx0(i)*xx1(i)
        aaa(2,1)=aaa(2,1)+yy0(i)*xx1(i)
        aaa(3,1)=aaa(3,1)+zz0(i)*xx1(i)
        aaa(1,2)=aaa(1,2)+xx0(i)*yy1(i)
        aaa(2,2)=aaa(2,2)+yy0(i)*yy1(i)
        aaa(3,2)=aaa(3,2)+zz0(i)*yy1(i)
        aaa(1,3)=aaa(1,3)+xx0(i)*zz1(i)
        aaa(2,3)=aaa(2,3)+yy0(i)*zz1(i)
        aaa(3,3)=aaa(3,3)+zz0(i)*zz1(i)
        
      enddo

c     construct optimisation matrix

      mat(1,1)=aaa(1,1)+aaa(2,2)+aaa(3,3)
      mat(2,2)=aaa(1,1)-aaa(2,2)-aaa(3,3)
      mat(3,3)=aaa(2,2)-aaa(1,1)-aaa(3,3)
      mat(4,4)=aaa(3,3)-aaa(2,2)-aaa(1,1)
      mat(1,2)=aaa(2,3)-aaa(3,2)
      mat(2,1)=aaa(2,3)-aaa(3,2)
      mat(1,3)=aaa(3,1)-aaa(1,3)
      mat(3,1)=aaa(3,1)-aaa(1,3)
      mat(1,4)=aaa(1,2)-aaa(2,1)
      mat(4,1)=aaa(1,2)-aaa(2,1)
      mat(2,3)=aaa(1,2)+aaa(2,1)
      mat(3,2)=aaa(1,2)+aaa(2,1)
      mat(2,4)=aaa(1,3)+aaa(3,1)
      mat(4,2)=aaa(1,3)+aaa(3,1)
      mat(3,4)=aaa(2,3)+aaa(3,2)
      mat(4,3)=aaa(2,3)+aaa(3,2)

c     diagonalise optimisation matrix

      call jacobi(4,mat,vec)

c     find largest eigenvalue

      k=1
      do i=2,4

        if(mat(i,i).gt.mat(k,k))k=i

      enddo

c     obtain optimal quaternion

      do i=1,4

        qqq(i)=vec(i,k)

      enddo

c     construct rotation matrix

      rot(1,1)=qqq(1)**2+qqq(2)**2-qqq(3)**2-qqq(4)**2
      rot(1,2)=2.d0*(qqq(2)*qqq(3)+qqq(1)*qqq(4))
      rot(1,3)=2.d0*(qqq(2)*qqq(4)-qqq(1)*qqq(3))
      rot(2,1)=2.d0*(qqq(2)*qqq(3)-qqq(1)*qqq(4))
      rot(2,2)=qqq(1)**2-qqq(2)**2+qqq(3)**2-qqq(4)**2
      rot(2,3)=2.d0*(qqq(3)*qqq(4)+qqq(1)*qqq(2))
      rot(3,1)=2.d0*(qqq(2)*qqq(4)+qqq(1)*qqq(3))
      rot(3,2)=2.d0*(qqq(3)*qqq(4)-qqq(1)*qqq(2))
      rot(3,3)=qqq(1)**2-qqq(2)**2-qqq(3)**2+qqq(4)**2

c     rotate trial structure

      do i=1,natms

        txx=xx1(i)
        tyy=yy1(i)
        tzz=zz1(i)
        xx1(i)=rot(1,1)*txx+rot(1,2)*tyy+rot(1,3)*tzz
        yy1(i)=rot(2,1)*txx+rot(2,2)*tyy+rot(2,3)*tzz
        zz1(i)=rot(3,1)*txx+rot(3,2)*tyy+rot(3,3)*tzz

      enddo

      return
      end

      function offset(natms,xx0,yy0,zz0,xx1,yy1,zz1)

c*********************************************************************
c
c     calculate degree of fit between trial and template structures
c
c     copyright daresbury laboratory
c     author w.smith march 2004
c
c*********************************************************************

      implicit none

      integer i,natms
      real*8 xx0,yy0,zz0,xx1,yy1,zz1,offset

      dimension xx0(*),yy0(*),zz0(*),xx1(*),yy1(*),zz1(*)

      offset=0.d0

      do i=1,natms

        offset=offset+(xx0(i)-xx1(i))**2+(yy0(i)-yy1(i))**2+
     x    (zz0(i)-zz1(i))**2

      enddo

      offset=sqrt(offset/dble(natms))

      return
      end

      subroutine jacobi(n,a,v)

c***********************************************************************
c     
c     diagonalisation of real symmetric matices by jacobi method
c     
c     input parameters:
c     
c     a(n,n) is the matrix to be diagonalised
c     v(n,n) is the eigenvector matrix
c     n   is the dimension of the matrices
c     
c     jacobi processes lower triangle only (upper triangle unchanged)
c     
c     variable rho sets absolute tolerance on convergence
c     variable tes is a moving tolerance that diminishes
c     on each pass until at true convergence tes<rho
c     
c     author w.smith 1993
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      logical pass

      integer n,i,j,k
      real*8 a,v,rho,tes,scl,v1,v2,v3,omg,s,c,u,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        enddo
        v(i,i)=1.0d0
      enddo

c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl)scl=abs(a(i,i))
      enddo
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        enddo
      enddo

c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

c     recycle until absolute tolerance satisfied
        
      do while(tes.gt.rho)
        
        pass=.true.
        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho
        
c     recycle until moving tolerance satisfied
        
        do while(pass)
          
          pass=.false.
          
c     jacobi diagonalisation

          do i=2,n
            
            do j=1,i-1
              
              if(abs(a(i,j)).ge.tes)then
                
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                
                if(abs(u).lt.rho)then
                  
                  omg=-1.0d0
                  
                else
                  
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                  
                endif
                
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                
                do k=1,n
                  
                  if(k.ge.i)then
                    
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                    
                  else if(k.lt.j)then
                    
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                    
                  else
                    
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                    
                  endif
                  
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                  
                enddo
                
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
                
              endif
              
            enddo
            
          enddo
          
        enddo
        
      enddo

c     rescale matrix

      do i=1,n

        do j=1,i

          a(i,j)=scl*a(i,j)

        enddo

      enddo

      return
      end

