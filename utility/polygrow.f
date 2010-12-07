      program polygrow
c     
c*********************************************************************
c     
c     dl_poly utility for generating an amorphous polymer chain
c     
c     copyright daresbury laboratory march 1995
c     
c     author    w.smith march 1995
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c
c*********************************************************************
c     
      implicit real*8(a-h,o-z)
      
      parameter (mxatms=906,nsearch=100,mxtry=25,mxfail=1000)
      
      logical opr,lrej

      dimension xbs(4),ybs(4),zbs(4),rot(9)
      dimension eps(3),sig(3),edih(4)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      
      data opr/.false./,pi/3.141592653589793d0/

c     enter the control variables
      
      write(*,*)'Enter number of monomers in chain'
      read(*,*)nlinks
      write(*,'(i12)')nlinks
      write(*,*)'Enter system volume and temperature'
      read(*,*)volm,temp
      write(*,'(2f12.4)')volm,temp
      write(*,*)'Enter C-C and C-H bondlengths (A)'
      read(*,*)ccbond,chbond
      write(*,'(2f12.4)')ccbond,chbond
      write(*,*)'Enter C-C bond energy (kJ/mol)'
      read(*,*)ebond
      write(*,'(2f12.4)')ebond
      write(*,*)'Enter C-C, C-H and H-H epsilons (kJ/mol)'
      read(*,*)eps(1),eps(2),eps(3)
      write(*,'(3f12.4)')eps(1),eps(2),eps(3)
      write(*,*)'Enter C-C, C-H and H-H sigmas (A)'
      read(*,*)sig(1),sig(2),sig(3)
      write(*,'(3f12.4)')sig(1),sig(2),sig(3)
      write(*,*)'Enter dihedral LJ scaling factor'
      read(*,*)dihscl
      write(*,'(f12.4)')dihscl
      write(*,*)'Enter dihedral energy parameters (kJ/mol)'
      read(*,*)edih(1),edih(2),edih(3),edih(4)
      write(*,'(4f12.4)')edih(1),edih(2),edih(3),edih(4)

c     test atom numbers

      matms=3*nlinks+2
      if(mxatms.lt.matms)then

        write(*,'(a,i6)')
     x    'Error parameter mxatms must be reset to',matms
        stop

      endif

c     calculate mass density (amu/A^3)

      rho=(14.0268d0*dble(nlinks)+2.0158d0)/volm

c     initialise random numbers
      
      start=duni()

c     linear width of cell

      size=volm**(1.d0/3.d0)

c     define tetrahedral groups
      
      xbs(1)=0.d0
      ybs(1)=0.d0
      zbs(1)=1.d0
      
      xbs(2)=2.d0*sqrt(2.d0)/3.d0
      ybs(2)=0.d0
      zbs(2)=-1.d0/3.d0
      
      xbs(3)=-sqrt(2.d0)/3.d0
      ybs(3)=-sqrt(2.d0/3.d0)
      zbs(3)=-1.d0/3.d0
      
      xbs(4)=-sqrt(2.d0)/3.d0
      ybs(4)=sqrt(2.d0/3.d0)
      zbs(4)=-1.d0/3.d0
      
c     position first tetrahedron
      
      natms=1
      xx0=0.d0
      yy0=0.d0
      zz0=0.d0
      base=volm**(1.d0/3.d0)
      alp=2.d0*pi*duni()
      bet=pi*duni()
      gam=2.d0*pi*duni()
      call euler(alp,bet,gam,rot)
      call rotate(opr,xbs(1),ybs(1),zbs(1),xxt,yyt,zzt,rot)
      xxx(1)=chbond*xxt+xx0
      yyy(1)=chbond*yyt+yy0
      zzz(1)=chbond*zzt+zz0
      
      link=0
      nstep=0
      do ntrial=1,nsearch*nlinks
         
         lhist=0
         nfail=0
         link=link+1
         if(link.gt.nlinks)go to 300

  100    continue

         lrej=.false.
         lhist=lhist+1

         if(lhist.gt.mxtry)then
         
            lhist=0
            nfail=nfail+1
            if(nfail.gt.mxfail)go to 200
            natms=max(4,natms-3*nint(3.d0*duni()))
            link=(natms-1)/3+1
            xx0=xxx(natms+1)
            yy0=yyy(natms+1)
            zz0=zzz(natms+1)
            xxt=xxx(natms-2)-xx0
            yyt=yyy(natms-2)-yy0
            zzt=zzz(natms-2)-zz0
            xxt=xxt-size*anint(xxt/size)
            yyt=yyt-size*anint(yyt/size)
            zzt=zzt-size*anint(zzt/size)
            alp=atan2(yyt,xxt)
            bet=acos(zzt/sqrt(xxt**2+yyt**2+zzt**2))
            gam=2.d0*pi*duni()

         endif

         call euler(alp,bet,gam,rot)

         xxx(natms+1)=xx0
         yyy(natms+1)=yy0
         zzz(natms+1)=zz0
         call rotate(opr,xbs(2),ybs(2),zbs(2),xxt,yyt,zzt,rot)
         xxx(natms+2)=chbond*xxt+xx0
         yyy(natms+2)=chbond*yyt+yy0
         zzz(natms+2)=chbond*zzt+zz0
         call rotate(opr,xbs(3),ybs(3),zbs(3),xxt,yyt,zzt,rot)
         xxx(natms+3)=chbond*xxt+xx0
         yyy(natms+3)=chbond*yyt+yy0
         zzz(natms+3)=chbond*zzt+zz0
         call rotate(opr,xbs(4),ybs(4),zbs(4),xxt,yyt,zzt,rot)
         xxx(natms+4)=chbond*xxt+xx0
         yyy(natms+4)=chbond*yyt+yy0
         zzz(natms+4)=chbond*zzt+zz0

c     check if added unit is energetically acceptable

         nstep=nstep+1
         if(link.gt.1)call select
     x     (lrej,natms,size,temp,ebond,dihscl,eps,sig,edih,xxx,yyy,zzz)

         gam=2.d0*pi*duni()

         if(lrej)go to 100

c     growth direction vector and rotation matrix

         xx0=ccbond*xxt+xx0
         yy0=ccbond*yyt+yy0
         zz0=ccbond*zzt+zz0
         alp=atan2(-yyt,-xxt)
         bet=acos(-zzt)
         
         natms=natms+3

      enddo

      write(*,'("error - unable to propagate growth")')
      write(*,'("number of successful links :",i5)')link
      call outcon(natms,rho,temp,xxx,yyy,zzz)
      stop

  200 continue

      write(*,'("error - too many trial failures")')
      write(*,'("number of successful links :",i5)')link
      call outcon(natms,rho,temp,xxx,yyy,zzz)
      stop

  300 continue

      write(*,'("number of trial moves:",i6)')nstep

      natms=natms+1
      call outcon(natms,rho,temp,xxx,yyy,zzz)

      end
      subroutine select
     x (lrej,natms,size,temp,ebond,dihscl,eps,sig,edih,xxx,yyy,zzz)

c*********************************************************************
c
c     dl_poly routine to select or reject a new CH2 monomer
c     added to a chain using the boltzman factor as the 
c     selection criterion
c
c     copyright daresbury laboratory march 1995
c     author w.smith march 1995
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      parameter (mxatms=906)

      logical lrej
      dimension eps(3),sig(3),edih(4)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      save heng
      data rgas/8.31451d-3/,heng/0.d0/

      
      lrej=.false.
      boltz=1.d0/(temp*rgas)
      eng=edih(1)-ebond-heng

c     bond vectors

      if(natms.eq.4) then

         b1x=xxx(2)-xxx(1)
         b1y=yyy(2)-yyy(1)
         b1z=zzz(2)-zzz(1)
         b1x=b1x-size*anint(b1x/size)
         b1y=b1y-size*anint(b1y/size)
         b1z=b1z-size*anint(b1z/size)

      else

         b1x=xxx(natms-2)-xxx(natms-5)
         b1y=yyy(natms-2)-yyy(natms-5)
         b1z=zzz(natms-2)-zzz(natms-5)
         b1x=b1x-size*anint(b1x/size)
         b1y=b1y-size*anint(b1y/size)
         b1z=b1z-size*anint(b1z/size)

      endif

      b2x=xxx(natms+1)-xxx(natms-2)
      b2y=yyy(natms+1)-yyy(natms-2)
      b2z=zzz(natms+1)-zzz(natms-2)
      b2x=b2x-size*anint(b2x/size)
      b2y=b2y-size*anint(b2y/size)
      b2z=b2z-size*anint(b2z/size)

      b3x=xxx(natms+4)-xxx(natms+1)
      b3y=yyy(natms+4)-yyy(natms+1)
      b3z=zzz(natms+4)-zzz(natms+1)
      b3x=b3x-size*anint(b3x/size)
      b3y=b3y-size*anint(b3y/size)
      b3z=b3z-size*anint(b3z/size)

c     calculate dihedral angle energy

      b12x=b1y*b2z-b2y*b1z
      b12y=b1z*b2x-b2z*b1x
      b12z=b1x*b2y-b2x*b1y
      bb12=sqrt(b12x**2+b12y**2+b12z**2)

      b23x=b2y*b3z-b3y*b2z
      b23y=b2z*b3x-b3z*b2x
      b23z=b2x*b3y-b3x*b2y
      bb23=sqrt(b23x**2+b23y**2+b23z**2)

      cost=(b12x*b23x+b12y*b23y+b12z*b23z)/(bb12*bb23)
      eng=eng+edih(1)+cost*(edih(2)+cost*(edih(3)+cost*edih(4)))

c     C-C backbone LJ interactions

      i=natms+1

      do j=2,natms-8,3
         
         ddx=xxx(i)-xxx(j)
         ddy=yyy(i)-yyy(j)
         ddz=zzz(i)-zzz(j)
         ddx=ddx-size*anint(ddx/size)
         ddy=ddy-size*anint(ddy/size)
         ddz=ddz-size*anint(ddz/size)
         rat=(sig(1)**2/(ddx*ddx+ddy*ddy+ddz*ddz))**3
         eng=eng+4.d0*eps(1)*(rat**2-rat)
         
      enddo

c     C-H interactions

      k=2
      do j=1,natms-3

         if(j.eq.k)then

            k=k+3

         else

            fact=1.d0
            if(j.gt.i-6.or.i.eq.8)fact=dihscl
            ddx=xxx(i)-xxx(j)
            ddy=yyy(i)-yyy(j)
            ddz=zzz(i)-zzz(j)
            ddx=ddx-size*anint(ddx/size)
            ddy=ddy-size*anint(ddy/size)
            ddz=ddz-size*anint(ddz/size)
            rat=(sig(2)**2/(ddx*ddx+ddy*ddy+ddz*ddz))**3
            eng=eng+4.d0*fact*eps(2)*(rat**2-rat)

         endif

      enddo
            
c     initialise energy for replacement Hydrogen atom

      teng=0.d0

c     H-H 1-4 LJ interactions

      do i=natms+2,natms+4

         k=2
         do j=1,natms

            if(j.eq.k)then
               
               k=k+3
               
            else
               
               fact=1.d0
               if(j.ge.natms-1.or.natms.eq.4)fact=dihscl
               ddx=xxx(i)-xxx(j)
               ddy=yyy(i)-yyy(j)
               ddz=zzz(i)-zzz(j)
               ddx=ddx-size*anint(ddx/size)
               ddy=ddy-size*anint(ddy/size)
               ddz=ddz-size*anint(ddz/size)
               rat=(sig(3)**2/(ddx*ddx+ddy*ddy+ddz*ddz))**3
               eterm=4.d0*fact*eps(3)*(rat**2-rat)
               eng=eng+eterm
               if(i.eq.natms+4)teng=teng+eterm
            
            endif

         enddo

      enddo
      
c     apply selection criterion

      if(eng.gt.0.d0)then

         if(duni().gt.exp(-eng*boltz))lrej=.true.

      endif
      if(.not.lrej)heng=teng

      return
      end
      function duni()
c     
c*********************************************************************
c     
c     dl_poly random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c     copyright daresbury laboratory 1992
c     author -  w.smith         july 1992
c     
c*********************************************************************
c     
      logical new
      real*8 duni
      real*4 u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./
      if(new)then
c     
c     initial values of i,j,k must be in range 1 to 178 (not all 1)
c     initial value of l must be in range 0 to 168.
         i=12
         j=34
         k=56
         l=78
c     
         ir=97
         jr=33
         new=.false.
         do 200 ii=1,97
            s=0.0
            t=0.5
            do 100 jj=1,24
               m=mod(mod(i*j,179)*k,179)
               i=j
               j=k
               k=m
               l=mod(53*l+1,169)
               if(mod(l*m,64).ge.32)s=s+t
               t=0.5*t
  100       continue
            u(ii)=s
  200    continue
         c =  362436.0/16777216.0
         cd= 7654321.0/16777216.0
         cm=16777213.0/16777216.0
      else
c     
c     calculate random number
         uni=u(ir)-u(jr)
         if(uni.lt.0.0)uni=uni+1.0
         u(ir)=uni
         ir=ir-1
         if(ir.eq.0)ir=97
         jr=jr-1
         if(jr.eq.0)jr=97
         c=c-cd
         if(c.lt.0.0)c=c+cm
         uni=uni-c
         if(uni.lt.0.0)uni=uni+1.0
         duni=dble(uni)
      endif
      return
      end
      subroutine euler(alp,bet,gam,rot)

c*********************************************************************
c
c     dl_poly routine to construct a rotation matrix as defined by the
c     Z-(-)Y'-Z'' euler angle convention and the equation r'=Rr, where
c     r' is a vector in a molecular frame, r the corresponding vector
c     in the laboratory frame and R is the rotation matrix
c
c     copyright daresbury laboratory march 1995
c
c     author w. smith march 1995
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      dimension rot(9)

      ca=cos(alp)
      cb=cos(bet)
      cg=cos(gam)
      sa=sin(alp)
      sb=sin(bet)
      sg=sin(gam)

      rot(1)= ca*cb*cg-sa*sg
      rot(2)=-ca*cb*sg-sa*cg
      rot(3)= ca*sb
      rot(4)= sa*cb*cg+ca*sg
      rot(5)=-sa*cb*sg+ca*cg
      rot(6)= sa*sb
      rot(7)=-sb*cg
      rot(8)= sb*sg
      rot(9)= cb

      return
      end
      subroutine rotate(op,xl,yl,zl,xg,yg,zg,rot)

c*********************************************************************
c
c     dl_poly routine to transform a vector from the local (molecule-
c     fixed) frame of reference (xl,yl,zl) to the global frame of 
c     reference (xg,yg,zg) and vice-versa using the rotation matrix 
c     defined by the matrix equation rl = (Rot) rg.
c
c     copyright daresbury laboratory march 1995
c     author w.smith march 1995
c
c     note: if op=.true.  converts global to local
c           if op=.false. converts local to global
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      logical op
      dimension rot(9)

      if(op)then

         xl=xg*rot(1)+yg*rot(4)+zg*rot(7)
         yl=xg*rot(2)+yg*rot(5)+zg*rot(8)
         zl=xg*rot(3)+yg*rot(6)+zg*rot(9)

      else

         xg=xl*rot(1)+yl*rot(2)+zl*rot(3)
         yg=xl*rot(4)+yl*rot(5)+zl*rot(6)
         zg=xl*rot(7)+yl*rot(8)+zl*rot(9)

      endif

      return
      end
      subroutine outcon(n,rho,temp,xxx,yyy,zzz)
c
c*********************************************************************

c     subroutine to write out configuration files for polygrow
c
c     copyright daresbury laboratory february 1995
c     author - w.smith february 1995
c
c*********************************************************************
c

      implicit real*8(a-h,o-z)

      dimension xxx(*),yyy(*),zzz(*)


c     open configuration file

      open(8,file='POLYCON')

      write(8,'(i6)')n
      write(8,'(a,2f12.4)')"Linear Chain Hydrocarbon (dens,temp)",
     x  rho,temp

c     write out configuration

      write(8,'(a8,3f15.8)')"H       ",xxx(1),yyy(1),zzz(1)
      do i=1,n-2,3

      write(8,'(a8,3f15.8)')"C       ",xxx(i+1),yyy(i+1),zzz(i+1)
      write(8,'(a8,3f15.8)')"H       ",xxx(i+2),yyy(i+2),zzz(i+2)
      write(8,'(a8,3f15.8)')"H       ",xxx(i+3),yyy(i+3),zzz(i+3)

      enddo
      write(8,'(a8,3f15.8)')"H       ",xxx(n),yyy(n),zzz(n)

c     close configuration file

      close (8)

      return
      end
