      subroutine jacobi(a,v,n)
c     
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
c     copyright daresbury laboratory 1993
c     author w.smith 1993
c
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c
c***********************************************************************
c     
      implicit real*8(a-h,o-z)      
      logical pass
      dimension a(n,n),v(n,n)
      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0
c     
c     initialize eigenvectors
      do i=1,n
         do j=1,n
            v(i,j)=0.0d0
         enddo
         v(i,i)=1.0d0
      enddo
c     
c     rescale matrix for optimal accuracy
      do i=1,n
         if(abs(a(i,i)).gt.scl)scl=abs(a(i,i))
      enddo
      do i=1,n
         do j=1,i
            a(i,j)=a(i,j)/scl
         enddo
      enddo
c     
c     set initial value of moving tolerance
      do i=2,n
         do j=1,i-1
            tes=tes+2.0d0*a(i,j)*a(i,j)
         enddo
      enddo
      tes=sqrt(tes)
  100 tes=tes/dble(n)
      if(tes.lt.rho)tes=rho
c     
c     jacobi diagonalisation
  200 pass=.false.
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
c     
c     recycle until moving tolerance satisfied
      if(pass)go to 200
c     
c     recycle until absolute tolerance satisfied
      if(tes.gt.rho)go to 100
c     
c     rescale matrix
      do i=1,n
         do j=1,i
            a(i,j)=scl*a(i,j)
         enddo
      enddo
      return
      end
