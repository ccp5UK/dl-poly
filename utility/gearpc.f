      subroutine gearpc
     x  (mode,idnode,mxnode,natms,imcon,engke,weight,cell,
     x  xxx,yyy,zzz,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,
     x  zz4,xx5,yy5,zz5,fxx,fyy,fzz,buffer)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - gear predictor corrector
c     
c     parallel replicated data version / block data
c     
c     copyright - daresbury laboratory 1993
c     author    - w. smith february 1993.
c     amended   - t.forester dec 1994 : block data
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
      include 'dl_params.inc'
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xx1(mxatms),yy1(mxatms),zz1(mxatms)
      dimension xx2(mxatms),yy2(mxatms),zz2(mxatms)
      dimension xx3(mxatms),yy3(mxatms),zz3(mxatms)
      dimension xx4(mxatms),yy4(mxatms),zz4(mxatms)
      dimension xx5(mxatms),yy5(mxatms),zz5(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension weight(mxatms),buffer(mxbuff),cell(9)
      
c     
c     block indices

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      if(mode.eq.1)then
c     
c     predictor step for atomic motion
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)+xx1(i)+xx2(i)+xx3(i)+xx4(i)+xx5(i)
          xx1(i)=xx1(i)+2.d0*xx2(i)+3.d0*xx3(i)+4.d0*xx4(i)+
     x      5.d0*xx5(i)
          xx2(i)=xx2(i)+3.d0*xx3(i)+6.d0*xx4(i)+10.d0*xx5(i)
          xx3(i)=xx3(i)+4.d0*xx4(i)+10.d0*xx5(i)
          xx4(i)=xx4(i)+5.d0*xx5(i)
          
        enddo
        
        do i=iatm1,iatm2
          
          yyy(i)=yyy(i)+yy1(i)+yy2(i)+yy3(i)+yy4(i)+yy5(i)
          yy1(i)=yy1(i)+2.d0*yy2(i)+3.d0*yy3(i)+4.d0*yy4(i)+
     x      5.d0*yy5(i)
          yy2(i)=yy2(i)+3.d0*yy3(i)+6.d0*yy4(i)+10.d0*yy5(i)
          yy3(i)=yy3(i)+4.d0*yy4(i)+10.d0*yy5(i)
          yy4(i)=yy4(i)+5.d0*yy5(i)
          
        enddo
        
        do i=iatm1,iatm2
          
          zzz(i)=zzz(i)+zz1(i)+zz2(i)+zz3(i)+zz4(i)+zz5(i)
          zz1(i)=zz1(i)+2.d0*zz2(i)+3.d0*zz3(i)+4.d0*zz4(i)+
     x      5.d0*zz5(i)
          zz2(i)=zz2(i)+3.d0*zz3(i)+6.d0*zz4(i)+10.d0*zz5(i)
          zz3(i)=zz3(i)+4.d0*zz4(i)+10.d0*zz5(i)
          zz4(i)=zz4(i)+5.d0*zz5(i)
          
        enddo
        
c     
c     calculate predicted value of kinetic energy
        
        engke=0.d0
        
        do i=iatm1,iatm2
          
          engke=engke+0.5d0*weight(i)*(xx1(i)**2+yy1(i)**2+zz1(i)**2)
          
        enddo
        
        if(mxnode.gt.1)call gdsum(engke,1,buffer)
        
        return
        
c     
c     calculate corrected motion of atoms
        
      else if(mode.eq.2)then
        
c     
c     calculate correction vector
        
        do i=iatm1,iatm2
          
          fxx(i)=0.5d0*fxx(i)/weight(i)-xx2(i)
          fyy(i)=0.5d0*fyy(i)/weight(i)-yy2(i)
          fzz(i)=0.5d0*fzz(i)/weight(i)-zz2(i)
          
        enddo
        
        do i=iatm1,iatm2
          
          xx1(i)=xx1(i)+(251.0d0/360.0d0)*fxx(i)
          xx3(i)=xx3(i)+(11.0d0/18.0d0)*fxx(i)
          xxx(i)=xxx(i)+(3.0d0/16.0d0)*fxx(i)
          xx5(i)=xx5(i)+(1.0d0/60.0d0)*fxx(i)
          xx4(i)=xx4(i)+(1.0d0/6.0d0)*fxx(i)
          xx2(i)=xx2(i)+fxx(i)
          
        enddo
        
        do i=iatm1,iatm2
          
          yy1(i)=yy1(i)+(251.0d0/360.0d0)*fyy(i)
          yy3(i)=yy3(i)+(11.0d0/18.0d0)*fyy(i)
          yyy(i)=yyy(i)+(3.0d0/16.0d0)*fyy(i)
          yy5(i)=yy5(i)+(1.0d0/60.0d0)*fyy(i)
          yy4(i)=yy4(i)+(1.0d0/6.0d0)*fyy(i)
          yy2(i)=yy2(i)+fyy(i)
          
        enddo
        
        do i=iatm1,iatm2
          
          zz1(i)=zz1(i)+(251.0d0/360.0d0)*fzz(i)
          zz3(i)=zz3(i)+(11.0d0/18.0d0)*fzz(i)
          zzz(i)=zzz(i)+(3.0d0/16.0d0)*fzz(i)
          zz5(i)=zz5(i)+(1.0d0/60.0d0)*fzz(i)
          zz4(i)=zz4(i)+(1.0d0/6.0d0)*fzz(i)
          zz2(i)=zz2(i)+fzz(i)
          
        enddo
        
c     
c     periodic boundary condition
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        
c     
c     calculate corrected value of kinetic energy
        
        engke=0.d0
        
        do i=iatm1,iatm2
          
          engke=engke+0.5d0*weight(i)*(xx1(i)**2+yy1(i)**2+zz1(i)**2)
          
        enddo
        
c     
c     global exchange of configuration data (note higher derivatives
c     are not exchanged)
        
        if(mxnode.gt.1)then
          
          nbuff=mxbuff
          call gdsum(engke,1,buffer)
          call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,natms,nbuff,xx1,yy1,zz1,buffer)
          
        endif
        
      endif
      
      return
      end
