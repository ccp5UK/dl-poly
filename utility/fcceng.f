c*********************************************************************
c
c     dl_poly utility for calculating energy of perfect LJ FCC lattice
c
c     author - w.smith Nov. 1998
c     copyright daresbury laboratory 1998
c
c*********************************************************************

      implicit none
      real*8 rd2,www,eps,sig,rcut,dis,rcut2,sig2,dis2,rr2
      real*8 eng,rho,pi,vir
      integer i,nnn
      dimension rd2(20),www(20)

      data pi/3.1415926536d0/
      data rd2/1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0,11.d0,
     x  12.d0,13.d0,15.d0,16.d0,17.d0,18.d0,19.d0,20.d0,21.d0/
      data www/12.d0,6.d0,24.d0,12.d0,24.d0,8.d0,48.d0,6.d0,36.d0,24.d0,
     x  24.d0,24.d0,72.d0,48.d0,12.d0,48.d0,30.d0,72.d0,24.d0,48.d0/

      write(*,*)'Input epsilon'
      read(*,*)eps
      write(*,*)'Input sigma'
      read(*,*)sig
      write(*,*)'Input cutoff'
      read(*,*)rcut
      write(*,*)'Input n-n distance'
      read(*,*)dis
      write(*,*)'Input no. atoms'
      read(*,*)nnn

      sig2=sig**2
      dis2=dis**2
      rcut2=rcut**2
      rho=4.d0/(2.d0**1.5d0*dis**3)
      eng=16.d0*pi*rho*eps*sig**3*((sig/rcut)**9/9.d0-
     x    (sig/rcut)**3/3.d0)
      vir=16.d0*pi*rho*eps*sig**3*(-4.d0*(sig/rcut)**9/3.d0+
     x    2.d0*(sig/rcut)**3)

      do i=1,20

        rr2=dis2*rd2(i)
        if(rr2.le.rcut2)then

          eng=eng+4.d0*eps*www(i)*((sig2/rr2)**6-(sig2/rr2)**3)
          vir=vir-24.d0*eps*www(i)*(2.d0*(sig2/rr2)**6-(sig2/rr2)**3)

        endif

      enddo

      write(*,'(a,1pe14.6)')'Site   Energy = ',eng
      write(*,'(a,1pe14.6)')'Site   Virial = ',vir
      write(*,'(a,1pe14.6)')'System Energy = ',eng*dble(nnn)/2.d0
      write(*,'(a,1pe14.6)')'System Virial = ',vir*dble(nnn)/2.d0

      end
