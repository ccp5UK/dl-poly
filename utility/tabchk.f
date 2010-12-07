      program tabchk
c     
c***********************************************************************
c     
c     dl_poly utility to check accuracy of tabulation procedures
c     
c     copyright daresbury laboratory 1994
c     author -    t.forester october 1994
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c
c***********************************************************************
c     
      implicit real*8(a-h,o-z)
      parameter (mega=10000,mxvdw=1000,mxpvdw=5)
      parameter (nfield=10)
      
      character*8  atom1,atom2
      character*80 record
      
      dimension parvdw(mxpvdw)
      dimension lstvdw(mxvdw),ltpvdw(mxvdw),prmvdw(mxvdw,mxpvdw)
      
      write(*,*) 'trial mxgrid '
      read(*,*) ngrid
      ngrid = max(ngrid,5)

      write(*,*)
      write(*,*) 'trial cutoff ?'
      read(*,*) rcut

      write(*,*)
      write(*,*) ' what is the minimum permitted pair separation ?'
      read(*,*) rmin

      write(*,*)
      write(*,*) 'number of test points ? (eg. 2311) '
      read(*,*) npts
      npts = max(npts,2)
      
c     
c     increment for distance
      
      dr = (rcut-rmin)/dble(npts-1)
c     
c     set up constants
      idnode=0
      mxnode = 1
      keyvdw = 0
      
      open(nfield,file='FIELD')
c     
c     search for vdw tables
      
      do irec = 1,mega
        
        read(nfield,'(a80)') record
        
        call lowcase(record,80)
        call strip(record,80)
        
        if(record(1:3).eq.'vdw') then
          
          ntpvdw=intstr(record,80,idum)
          
          if(ntpvdw.gt.mxvdw) call error(idnode,80)
          
          do ivdw=1,mxvdw
            lstvdw(ivdw)=0
          enddo
          
          do itpvdw=1,ntpvdw
            
            do i=1,mxpvdw
              parvdw(i)=0.d0
            enddo
            
            read(nfield,'(2a8,i5,5f12.0)',end=9999)
     x         atom1,atom2,keypot,parvdw
            
            keyvdw=keyvdw + 1
            
            if(keyvdw.gt.mxvdw) call error(idnode,82)
            
            lstvdw(keyvdw)=itpvdw
            ltpvdw(itpvdw)=keypot
            
            do i=1,5
              
              prmvdw(itpvdw,i)=parvdw(i)
              
            enddo
            
          enddo
          
        elseif(record(1:5).eq.'close')then
          
          close (nfield)
          goto 100
          
        endif
        
      enddo
      
c     
c     uncontrolled error exit from field file procesing
      
      close (nfield)
      call error(idnode,16)
      
c     
c     end of field file error exit
      
 9999 close (nfield)
      call error(idnode,52)
c     
c     evaluate accuracy of tabulation procedures
      
  100 continue
      
c
c     begin testing

      do iter = 1,3
        
        write(*,*)
        write(*,*)

        if(iter.eq.1) then
c     
c     three point interpolation
          
          write(*,*) 'Three point interpolation: '
          
          dlrpot = rcut/dble(ngrid-4)
          
        elseif(iter.eq.2) then
          
          write(*,*) 'four  point interpolation '
          dlrpot = rcut/dble(ngrid-4)
          
        elseif(iter.eq.3) then
          
          write(*,*) 'r-squared 2 point interpolation '
          
          dlrpot = (rcut**2)/dble(ngrid-4)
          
        endif
        
        write(*,*) ' functions : ',ntpvdw
        write(*,'(/,45x,a,60x,a)') 'energies','forces'
        write(*,'(1x,4(22x,a))') 'absolute error','relative error',
     x     'absolute error','relative error'
        write(*,'(a,5x,4(a,3x))') 'ivdw',
     x     'maximum      mean     variance   ',
     x     'maximum      mean     variance   ',
     x     '   maximum      mean     variance   ',
     x     'maximum      mean     variance'
        
        rdlpot = 1.d0/dlrpot
        
        do ivdw = 1,ntpvdw
          
          ltp = ltpvdw(ivdw)
          
          aev=0.d0
          aeg=0.d0
          rev=0.d0
          reg=0.d0
          accv=0.d0
          ac2v=0.d0
          arcv=0.d0
          ar2v=0.d0
          accg=0.d0
          ac2g=0.d0
          arcg=0.d0
          ar2g=0.d0
          
          do ii = 0,npts-1
            
            rrr = rmin + dble(ii)*dr
            
            call lookup(ltp,ivdw,rrr,exactv,exactg,prmvdw)
 
            if(iter.eq.1) then
              
              call three(ltp,ivdw,rrr,dlrpot,engsrp,gamma,prmvdw)
              
            elseif(iter.eq.2) then
              
              call four(ltp,ivdw,rrr,dlrpot,engsrp,gamma,prmvdw)
              
            elseif(iter.eq.3) then
              
              call rsquare(ltp,ivdw,rrr,dlrpot,engsrp,gamma,prmvdw)
              
            endif
            
c     
c     absolute errors
            
            errv = engsrp - exactv
            errg = gamma  - exactg
c     
c     make sure both non-zero!!!

            ev = sign(max(abs(exactv),1.d-6),exactv)
            eg = sign(max(abs(exactg),1.d-6),exactg)
c
c     relative errors

            relv = errv/ev
            relg = errg/eg
c
c     maximum solute error in forces and energy

            aev = max(aev,abs(errv))
            aeg = max(aeg,abs(errg))

            rev= max(rev,abs(relv))
            reg= max(reg,abs(relg))
c
c     accumulators for mean error + sd

            accv = accv + errv
            ac2v = ac2v + errv*errv

            accg = accg + errg
            ac2g = ac2g + errg*errg
c     
c     accumulators for mean relative error + sd

            arcv = arcv + relv
            ar2v = ar2v + relv*relv

            arcg = arcg + relg
            ar2g = ar2g + relg*relg

          enddo
c
c     report on errors

          averv = accv/dble(npts)
          vdv   = abs((ac2v - accv**2)/dble(npts-1))

          averg = accg/dble(npts)
          vdg   = abs((ac2g - accg**2)/dble(npts-1))

          avrv = arcv/dble(npts)
          vrdv   = abs((ar2v - arcv**2)/dble(npts-1))

          avrg = arcg/dble(npts)
          vrdg   = abs((ar2g - arcg**2)/dble(npts-1))

          write(*,'(i4,1p,6e12.4,3x,6e12.4)')
     x       ivdw,aev,averv,vdv,rev,avrv,vrdv,aeg,averg,vdg,
     x       reg,avrg,vrdg

        enddo

      enddo
      
      end
      
      subroutine lookup(ltp,ivdw,rrr,eng,for,prmvdw)
      
      implicit real*8(a-h,o-z)
      parameter (mxvdw=1000,mxpvdw=5)
      dimension prmvdw(mxvdw,mxpvdw)
      
c     
c     12 - 6 potential
      
      vv1(r,a,b)=(a/r**6-b)/r**6      
      gg1(r,a,b)=6.d0*(2.d0*a/r**6-b)/r**6
c     
c     lennard-jones potential
      
      vv2(r,a,b)=4.d0*a*(b/r)**6*((b/r)**6-1.d0)
      gg2(r,a,b)=24.d0*a*(b/r)**6*(2.d0*(b/r)**6-1.d0)
c     
c     n - m potential
      
      vv3(r,a,b,c,d)=a/(b-c)*(c*(d/r)**b-b*(d/r)**c)
      gg3(r,a,b,c,d)=a*c*b/(b-c)*((d/r)**b-(d/r)**c)
c     
c     buckingham exp - 6 potential
      
      vv4(r,a,b,c)=a*exp(-r/b)-c/r**6
      gg4(r,a,b,c)=r*a*exp(-r/b)/b-6.d0*c/r**6
c     
c     born-huggins-meyer exp - 6 - 8 potential
      
      vv5(r,a,b,c,d,e)=a*exp(b*(c-r))-d/r**6-e/r**8
      gg5(r,a,b,c,d,e)=r*a*b*exp(b*(c-r))-6.d0*d/r**6-8.d0*e/r**8
      
c     
c     Hydrogen-bond 12 - 10 potential
      
      vv6(r,a,b) = a/r**12 - b/r**10
      gg6(r,a,b) = 12.0d0*a/r**12 - 10.d0*b/r**10
      
      if(ltp.eq.1) then
        
        eng = vv1(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
        for = gg1(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
        
      elseif(ltp.eq.2) then
        
        eng = vv2(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
        for = gg2(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
        
      elseif(ltp.eq.3) then
        
        eng = vv3(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x     prmvdw(ivdw,3),prmvdw(ivdw,4))
        for = gg3(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x     prmvdw(ivdw,3),prmvdw(ivdw,4))
        
      elseif(ltp.eq.4) then
        
        eng = vv4(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x     prmvdw(ivdw,3))
        for = gg4(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x     prmvdw(ivdw,3))
        
      elseif(ltp.eq.5) then
        
        eng = vv5(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x     prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5))
        for = gg5(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x     prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5))
        
        
      elseif(ltp.eq.6) then
        
        eng = vv6(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
        for = gg6(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
        
      else

        call error(0,150)

      endif

      return
      end
      
      subroutine three(ltp,ivdw,rrr,dlrpot,engsrp,gamma,prmvdw)
c     
c     mimic 3pt interpolation

      implicit real*8(a-h,o-z)
      parameter (mxvdw=1000,mxpvdw=5)
      dimension prmvdw(mxvdw,mxpvdw)
      
      rdlpot = 1.d0/dlrpot
      l=int(rrr*rdlpot)
      ppp=rrr*rdlpot-dble(l)
      
c     
c     calculate interaction energy using 3-point interpolation
      
      r0 = dble(l)*dlrpot
      call lookup(ltp,ivdw,r0,vk,gk,prmvdw)
      
      r1 = dble(l+1)*dlrpot
      call lookup(ltp,ivdw,r1,vk1,gk1,prmvdw)
      
      r2 = dble(l+2)*dlrpot
      call lookup(ltp,ivdw,r2,vk2,gk2,prmvdw)
      
      t1 = vk + (vk1-vk)*ppp
      t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)
      
      engsrp = t1 + (t2-t1)*ppp*0.5d0
c     
c     calculate forces using 3-point interpolation
      
      t1 = gk + (gk1-gk)*ppp
      t2 = gk1 + (gk2-gk1)*(ppp - 1.0d0)
      
      gamma = (t1 +(t2-t1)*ppp*0.5d0)
      
      return
      end
      
      subroutine four(ltp,ivdw,rrr,dlrpot,engsrp,gamma,prmvdw)
c     
c     mimic 4pt interpolation
      
      implicit real*8(a-h,o-z)
      parameter (mxvdw=1000,mxpvdw=5)
      dimension prmvdw(mxvdw,mxpvdw)
      
      rdlpot = 1.d0/dlrpot
      l=int(rrr*rdlpot)
      ppp=rrr*rdlpot-dble(l)
      
c     
c     calculate interaction energy using 3-point interpolation
      
      rm = dble(l-1)*dlrpot
      call lookup(ltp,ivdw,rm,vkm,gkm,prmvdw)
      
      r0 = dble(l)*dlrpot
      call lookup(ltp,ivdw,r0,vk,gk,prmvdw)
      
      r1 = dble(l+1)*dlrpot
      call lookup(ltp,ivdw,r1,vk1,gk1,prmvdw)
      
      r2 = dble(l+2)*dlrpot
      call lookup(ltp,ivdw,r2,vk2,gk2,prmvdw)
      
      engsrp=vk+
     x   ppp*(-2.d0*vkm-3.d0*vk+6.d0*vk1-vk2+
     x   ppp*(3.d0*(vkm-vk-vk+vk1)+
     x   ppp*(-vkm+vk2+3.d0*(vk-vk1))))/6.d0
      
c     
c     calculate forces using 4-point interpolation
      
      gamma=gk+
     x   ppp*(-2.d0*gkm-3.d0*gk+6.d0*gk1-gk2+
     x   ppp*(3.d0*(gkm-gk-gk+gk1)+
     x   ppp*(-gkm+gk2+3.d0*(gk-gk1))))/6.d0
      
      
      return
      end
      
      subroutine rsquare(ltp,ivdw,rrr,dlrpot,engsrp,gamma,prmvdw)
c     
c     mimic r-squared interpolation
      
      implicit real*8(a-h,o-z)
      parameter (mxvdw=1000,mxpvdw=5)
      dimension prmvdw(mxvdw,mxpvdw)
      
      rdlpot = 1.d0/dlrpot
      rsq = rrr*rrr
      l=int(rsq*rdlpot)
      ppp=rsq*rdlpot-dble(l)
      
c     
c     calculate interaction energy using 2-point interpolation
      
      r0 = sqrt(dble(l)*dlrpot)
      call lookup(ltp,ivdw,r0,vk,gk,prmvdw)
      
      r1 = sqrt(dble(l+1)*dlrpot)
      call lookup(ltp,ivdw,r1,vk1,gk1,prmvdw)
      
      engsrp=vk+ (vk1-vk)*ppp
c     
c     calculate forces using 2-point interpolation
      
      gamma=gk+(gk1-gk)*ppp
      
      return
      end

      subroutine strip(string,length)

c***********************************************************************
c
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      character*(*) string
      
      imax = min(length,255)
      do i = 1,imax

         if(string(1:1).eq.' ') then

            do j = 1,imax-1

               string(j:j) = string(j+1:j+1)

            enddo

            string(imax:imax) = ' '

         endif

      enddo

      return
      end

      subroutine lowcase(string,length)

c***********************************************************************
c
c     DL_POLY routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      character*(*) string
      character*1 letter

      do i = 1,min(255,length)

         letter = string(i:i)

         if(letter.eq.'A') letter = 'a'
         if(letter.eq.'B') letter = 'b'
         if(letter.eq.'C') letter = 'c'
         if(letter.eq.'D') letter = 'd'
         if(letter.eq.'E') letter = 'e'
         if(letter.eq.'F') letter = 'f'
         if(letter.eq.'G') letter = 'g'
         if(letter.eq.'H') letter = 'h'
         if(letter.eq.'I') letter = 'i'
         if(letter.eq.'J') letter = 'j'
         if(letter.eq.'K') letter = 'k'
         if(letter.eq.'L') letter = 'l'
         if(letter.eq.'M') letter = 'm'
         if(letter.eq.'N') letter = 'n'
         if(letter.eq.'O') letter = 'o'
         if(letter.eq.'P') letter = 'p'
         if(letter.eq.'Q') letter = 'q'
         if(letter.eq.'R') letter = 'r'
         if(letter.eq.'S') letter = 's'
         if(letter.eq.'T') letter = 't'
         if(letter.eq.'U') letter = 'u'
         if(letter.eq.'V') letter = 'v'
         if(letter.eq.'W') letter = 'w'
         if(letter.eq.'X') letter = 'x'
         if(letter.eq.'Y') letter = 'y'
         if(letter.eq.'Z') letter = 'z'

         string(i:i) = letter

      enddo

      return
      end

      function intstr(word,len,lst)

c     
c***********************************************************************
c     
c     dl_poly function for extracting integers from a 
c     character string
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c     
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)
      character*1 n(0:9),word(len),ksn
      logical flag,final
      
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      ksn='+'
      intstr=0
      final=.false.
      
      do lst=1,len
         
         flag=.false.

         do j=0,9
            
            if(n(j).eq.word(lst))then
               
               intstr=10*intstr+j
               flag=.true.
               
            endif
            
         enddo

         if(flag.and.ksn.eq.'-')isn=-1
         ksn=word(lst)

         if(flag)then

            final=.true.

         else

            if(final)then

               intstr=isn*intstr
               return

            endif

         endif

      enddo
      
      return
      
      end

      subroutine error(idnode,iode)
c     
c***********************************************************************
c     
c     dl_poly subroutine for printing error messages and bringing
c     about a controlled termination of the program
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     amended   - t.forester sept 1994
c     
c     warning - this routine will terminate the job if iode > 0.
c     Users must ensure that all nodes are informed of error condition
c     efore this subroutine is called. e.g. using subroutine gstate().
c     
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)
      
      logical kill
      
      kill=(iode.ge.0)
      kode = abs(iode)
      
      if(idnode.eq.0) then
        
        if(kill) then
          write(nrite,'(//,1x,a,i5)') 
     x       'TABCHK terminated due to error ', kode
          
        else
          
          write(nrite,'(//,1x,a,i5)') 
     x       'TABCHK will terminate due to error ', kode
          
        endif
        
        if (kode.eq. 0) then
c     
c     dummy entry
          
        elseif (kode.eq.16) then
          write(nrite,'(//,1x,a)')
     x       'error - strange exit from FIELD file processing'
        elseif (kode.eq.52) then
          write(nrite,'(//,1x,a)')
     x       'error - end of FIELD file encountered'
        elseif (kode.eq.80) then
          write(nrite,'(//,1x,a)')
     x       'error - too many pair potentials specified'
        elseif (kode.eq.82) then
          write(nrite,'(//,1x,a)')
     x       'error - calculated pair potential index too large'
        elseif (kode.eq.150) then
          write(nrite,'(//,1x,a)')
     x       'error - unknown van der waals potential selected'
        else
          write(nrite,'(//,1x,a)')
     x       'error - unnamed error found'
        endif
        
      endif
      
      if (kill) then

       call exit()
        
      endif
      
      return
      end







