      subroutine ewald1a
     x  (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x  engcpe,vircpe,alpha,volm,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     
c     parallel replicated data version (part 1)
c     k vectors distributed over processors
c     
c     copyright - daresbury laboratory 1993.
c     author    - t. forester october 1993.
c     
c     version 2
c     author    - t. forester april 1993
c     
c     part 1 - reciprocal space terms (fourier part)
c     
c     note - in loop over all k vectors k=2pi(ll/cl,mm/cl,nn/cl)
c     the values of ll,mm and nn are selected so that the symmetry of
c     reciprocal lattice is taken into account i.e. the following
c     rules apply.
c     
c     ll ranges over the values 0 to kmax1 only.
c     
c     mm ranges over 0 to kmax2 when ll=0 and over
c     -kmax2 to kmax2 otherwise.
c     nn ranges over 1 to kmax3 when ll=mm=0 and over
c     -kmax3 to kmax3 otherwise.
c     
c     hence the result of the summation must be doubled at the end.
c     
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c     
c***********************************************************************
      
      use setup_module
      use config_module
      use ewald_module
      use property_module
      
      implicit none

      logical newjob,lconsw,leven
      integer idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3
      integer i,j,limit,ivec,l,ll,m,mm,n,nn,mmin,nmin,ierr
      real*8 engcpe,vircpe,alpha,volm,epsq,rcell,omg
      real*8 twopi,rvolm,ralph,det,rcpcut,rcpct2,engsic,ssx
      real*8 ssy,ssz,rkx1,rky1,rkz1,rkx2,rky2,rkz2,rkx3,rky3
      real*8 rkz3,cs,rksq,ckcs,ckss,rrksq,akk,akv,bkk,virprs
      real*8 qforce,eng1,scal1,scale

      dimension rcell(9),omg(9)

      save newjob,engsic
      
      data newjob/.true./,lconsw/.true./,leven/.true./

      twopi=2.d0*pi

CDPP ifdef VAMPIR
      call VTBEGIN(98, ierr)
CDPP endif
      if(mxewld.ne.mxatms) call error(idnode,330)

c     initialise coulombic potential energy
      
      engcpe=0.d0

c     initalize stress tensor working arrays

      do i=1,9
        omg(i)=0.d0
      enddo

c     set working parameters

      rvolm=twopi/volm
      ralph=-0.25d0/alpha**2

c     construct reciprocal lattice vectors and set k vector range

      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)
      call dcell(rcell,buffer)

      rcpcut=min(dble(kmax1)*buffer(7),dble(kmax2)*buffer(8),
     x  dble(kmax3)*buffer(9))
      rcpcut=rcpcut*1.05d0*twopi
      rcpct2=rcpcut**2

c     set switch for TO, RD and HP boundary conditions

      if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) then

        lconsw=.false.
        rvolm=0.5d0*rvolm

      endif
      
      if(newjob)then

c     calculate self interaction correction
        
        engsic=0.d0
        
        do i=idnode+1,natms,mxnode
          
          engsic=engsic+chge(i)**2
          
        enddo

        engsic=-r4pie0/epsq*alpha*engsic/sqrpi
        newjob=.false.
        
      endif

c     calculate and store exponential factors
c     convert real to reciprocal space coordinates
      
      i=0
      do j=1,natms
        
        i=i+1
        elc(i,0)=1.d0
        emc(i,0)=1.d0
        enc(i,0)=1.d0
        els(i,0)=0.d0
        ems(i,0)=0.d0
        ens(i,0)=0.d0
        ssx=rcell(1)*xxx(j)+rcell(4)*yyy(j)+rcell(7)*zzz(j)
        ssy=rcell(2)*xxx(j)+rcell(5)*yyy(j)+rcell(8)*zzz(j)
        ssz=rcell(3)*xxx(j)+rcell(6)*yyy(j)+rcell(9)*zzz(j)
        elc(i,1)=cos(twopi*ssx)
        emc(i,1)=cos(twopi*ssy)
        enc(i,1)=cos(twopi*ssz)
        els(i,1)=sin(twopi*ssx)
        ems(i,1)=sin(twopi*ssy)
        ens(i,1)=sin(twopi*ssz)
        
      enddo
      
      limit=i
      ivec=0

      do l=2,kmax2
        
        do i=1,limit
          
          emc(i,l)=emc(i,l-1)*emc(i,1)-ems(i,l-1)*ems(i,1)
          ems(i,l)=ems(i,l-1)*emc(i,1)+emc(i,l-1)*ems(i,1)
          
        enddo
        
      enddo

      do l=2,kmax3
        
        do i=1,limit
          
          enc(i,l)=enc(i,l-1)*enc(i,1)-ens(i,l-1)*ens(i,1)
          ens(i,l)=ens(i,l-1)*enc(i,1)+enc(i,l-1)*ens(i,1)
          
        enddo
        
      enddo

c     start of main loop over k vectors
      
      mmin=0
      nmin=1
      
      do ll=0,kmax1
        
        l=ll
        rkx1=twopi*dble(ll)*rcell(1)
        rky1=twopi*dble(ll)*rcell(4)
        rkz1=twopi*dble(ll)*rcell(7)

        if(l.ge.1) then
          
          do i=1,limit

            cs=elc(i,0)
            elc(i,0)=cs      *elc(i,1)-els(i,0)*els(i,1)
            els(i,0)=els(i,0)*elc(i,1)+cs      *els(i,1)
            
          enddo
          
        endif

        do mm=mmin,kmax2
          
          m=iabs(mm)
          rkx2=rkx1+twopi*dble(mm)*rcell(2)
          rky2=rky1+twopi*dble(mm)*rcell(5)
          rkz2=rkz1+twopi*dble(mm)*rcell(8)

c     set temporary products of exponential terms
          
          if(mm.ge.0)then
            
            do i=1,limit
              
              clm(i)=elc(i,0)*emc(i,m)-els(i,0)*ems(i,m)
              slm(i)=els(i,0)*emc(i,m)+ems(i,m)*elc(i,0)
              
            enddo
            
          else
            
            do i=1,limit
              
              clm(i)=elc(i,0)*emc(i,m)+els(i,0)*ems(i,m)
              slm(i)=els(i,0)*emc(i,m)-ems(i,m)*elc(i,0)
              
            enddo
            
          endif
          
          do nn=nmin,kmax3
            
            n=iabs(nn)

              if(.not.lconsw)then

                if(imcon.eq.4)then

                  leven=(mod(l+m+n,2).eq.0)

                elseif(imcon.eq.5)then

                  leven=(mod(l+m+n,2).eq.0)

                elseif(imcon.eq.7)then

                  leven=(mod(l+m,2).eq.0)

                endif

              endif

            if(lconsw.or.leven)then

              rkx3=rkx2+twopi*dble(nn)*rcell(3)
              rky3=rky2+twopi*dble(nn)*rcell(6)
              rkz3=rkz2+twopi*dble(nn)*rcell(9)
              
c     test on magnitude of k vector
              
              rksq=rkx3**2+rky3**2+rkz3**2

c     calculate exp(ikr) terms and product with charges

              if(rksq.le.rcpct2)then
                
c     test if vector is of interest to this processor

                ivec=ivec+1
                if(mod(ivec-1,mxnode).eq.idnode) then

                  i=0
                  if(nn.ge.0)then
                    
                    do j=1,natms
                      
                      i=i+1
                      ckc(i)=chge(j)*
     x                  (clm(i)*enc(i,n)-slm(i)*ens(i,n))
                      cks(i)=chge(j)*
     x                  (slm(i)*enc(i,n)+clm(i)*ens(i,n))
                      
                    enddo
                    
                  else
                    
                    do j=1,natms
                      
                      i=i+1
                      ckc(i)=chge(j)*
     x                  (clm(i)*enc(i,n)+slm(i)*ens(i,n))
                      cks(i)=chge(j)*
     x                  (slm(i)*enc(i,n)-clm(i)*ens(i,n))
                      
                    enddo
                    
                  endif

c     calculate vector sums
                  
                  ckcs=0.d0
                  ckss=0.d0
                  
                  do i=1,limit
                    
                    ckcs=ckcs+ckc(i)
                    ckss=ckss+cks(i)
                    
                  enddo
                  
c     calculate akk coefficients
                  
                rrksq=1.d0/rksq
                if(lconsw)then
                  akk=exp(ralph*rksq)*rrksq
                  akv=akk*2.d0*(rrksq-ralph)
                else
                  akk=4.0d0*exp(ralph*rksq)*rrksq
                  akv=akk*2.d0*(rrksq-ralph)
                endif
                bkk=akk

c     accumulate potential energy and virial terms
                
                engcpe=engcpe+akk*(ckcs*ckcs+ckss*ckss)
                virprs=akv*(ckcs*ckcs+ckss*ckss)
                omg(1)=omg(1)-virprs*rkx3*rkx3
                omg(5)=omg(5)-virprs*rky3*rky3
                omg(9)=omg(9)-virprs*rkz3*rkz3
CDPP ifdef STRESS
                omg(2)=omg(2)-virprs*rkx3*rky3
                omg(3)=omg(3)-virprs*rkx3*rkz3
                omg(6)=omg(6)-virprs*rky3*rkz3
CDPP endif

c     calculate force on each site
                  
                  i=0
                  
                  do j=1,natms
                    
                    i=i+1
                    qforce=bkk*(cks(i)*ckcs-ckc(i)*ckss)
                    fxx(j)=fxx(j)+rkx3*qforce
                    fyy(j)=fyy(j)+rky3*qforce
                    fzz(j)=fzz(j)+rkz3*qforce
                    
                  enddo
                  
                endif

c     end vector loop
              endif
              
            endif

          enddo
          
          nmin=-kmax3
          
        enddo
        
        mmin=-kmax2
        
      enddo
      
c     add self interaction correction to potential
      
      if(lconsw)then

        eng1=engcpe
        engcpe=2.d0*rvolm*r4pie0*engcpe/epsq+engsic
        scal1=2.d0*rvolm*r4pie0/epsq
        scale=4.d0*rvolm*r4pie0/epsq

      else

        eng1=engcpe
        engcpe=rvolm*r4pie0*engcpe/epsq+engsic
        scal1=rvolm*r4pie0/epsq
        scale=2.d0*rvolm*r4pie0/epsq

      endif

c     calculate final forces
      
      do i=1,natms
        
        fxx(i)=scale*fxx(i)
        fyy(i)=scale*fyy(i)
        fzz(i)=scale*fzz(i)

      enddo

CDPP ifdef STRESS

c     calculate stress tensor (symmetrical)

      stress(1)=stress(1)+scal1*(omg(1)+eng1)
      stress(2)=stress(2)+scal1*omg(2)
      stress(3)=stress(3)+scal1*omg(3)
      stress(4)=stress(4)+scal1*omg(2)
      stress(5)=stress(5)+scal1*(omg(5)+eng1)
      stress(6)=stress(6)+scal1*omg(6)
      stress(7)=stress(7)+scal1*omg(3)
      stress(8)=stress(8)+scal1*omg(6)
      stress(9)=stress(9)+scal1*(omg(9)+eng1)
CDPP endif

c     virial term

      vircpe=-scal1*(omg(1)+omg(5)+omg(9)+3.d0*eng1)

CDPP ifdef VAMPIR
      call VTEND(98, ierr)
CDPP endif

      return
      end

