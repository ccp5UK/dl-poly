      program shape

c***********************************************************************
c
c     calculate the shape of a gold cluster with attached chains
c
c     copyright daresbury laboratory
c     author w.smith
c
c***********************************************************************

      implicit none

      integer nread,nprint,nsave,mxatms,mxarms,mxcor
      parameter (nread=5,nprint=6,nsave=7)
      parameter (mxatms=10000,mxarms=5,mxcor=256)

      logical safe
      character*8 word,name

      integer i,imcon,j,k,keytrj,lenchn,length,m,narms,natms,nclust
      integer ngold,nstep,nstr

      real*8 a13,aaa,aat1,aat2,aat3,aln2,amass,arg2,arga,argb,argc
      real*8 bbb,ccc,cell,cmx,cmy,cmz,ddd,det,dmass,eee,espa,fat1
      real*8 fat2,fat3,fff,fln2,frg2,frga,frgb,frgc,ggg,hhh,hln2
      real*8 pat1,pat2,pat3,pi,pln2,prg2,prga,prgb,prgc,range
      real*8 ravstr,rcell,re3,re54,re9,rrr,rxx,rxy,rxz,ryy,ryz,rzz
      real*8 snstr1,snstr2,spha,sqr2,sqrq,sum,sxx,syy,szz,theta
      real*8 tstep,xdf,xln2,xxx,ydf,yyy,zdf,zzz

      dimension cell(9),rcell(9),name(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension amass(mxatms),hln2(mxcor),length(2,mxarms)

   10 format(1h ,'number of configurations processed =',i8)
   20 format(/,1h ,'average structure of clusters',/,
     x       1h ,'chain MS end-to-end    ',1p,e12.4,' +/- ',e12.4,/,
     x       1h ,'inertia tensor 1       ',1p,e12.4,' +/- ',e12.4,/,
     x       1h ,'inertia tensor 2       ',1p,e12.4,' +/- ',e12.4,/,
     x       1h ,'inertia tensor 3       ',1p,e12.4,' +/- ',e12.4,/,
     x       1h ,'moment of inertia      ',1p,e12.4,' +/- ',e12.4,/,
     x       1h ,'sphericity factor  a   ',1p,e12.4,' +/- ',e12.4,/,
     x       1h ,'sphericity factor <a>  ',1p,e12.4,' +/- ',e12.4)

      safe=.true.
      pi=3.1415926536d0

c     enter number of atoms, gold atoms, chains and chain length

      write(nprint,'(a)')'
     x  enter number of atoms, gold atoms, chains and chain length'

      read(nread,*)natms,ngold,narms,lenchn
      if(natms.gt.mxatms)then

        write(nprint,'(a,2i8)')
     x    'error - maximum number of atoms exceeded',mxatms,natms
        stop

      endif
      if(natms.gt.mxarms)then

        write(nprint,'(a,2i8)')
     x    'error - maximum number of chains exceeded',mxarms,narms
        stop

      endif
      nclust=ngold+lenchn*narms
      range=3.5d0*dble(lenchn)

c     zero the shape parameter variables

        aln2=0.d0
        arga=0.d0
        argb=0.d0
        argc=0.d0
        arg2=0.d0
        aat1=0.d0
        aat2=0.d0
        aat3=0.d0
        fln2=0.d0
        frga=0.d0
        frgb=0.d0
        frgc=0.d0
        frg2=0.d0
        fat1=0.d0
        fat2=0.d0
        fat3=0.d0

c     define chain start and end

      do k=1,narms

        length(1,k)=ngold+(k-1)*lenchn+1
        length(2,k)=ngold+k*lenchn

      enddo

c     open HISTORY file

      open(nsave,file="HISTORY")

      read(nsave,*)
      read(nsave,*)

      nstr=0

      do while(.true.)

        read(nsave,'(a8,4i10,f12.6)',end=100)
     x    word,nstep,natms,keytrj,imcon,tstep

        if(imcon.gt.0)then

          read(nsave,*)cell(1),cell(2),cell(3)
          read(nsave,*)cell(4),cell(5),cell(6)
          read(nsave,*)cell(7),cell(8),cell(9)

          call invert(cell,rcell,det)

        endif

        do i=1,natms

          read(nsave,'(a8)')name(i)
          read(nsave,*)xxx(i),yyy(i),zzz(i)
          if(keytrj.gt.0)read(nsave,*)aaa,bbb,ccc
          if(keytrj.gt.1)read(nsave,*)ddd,eee,fff

        enddo

        if(nstr.eq.1)then

          do i=1,nclust

            if(name(i)(1:2).eq.'Au')then

              amass(i)=196.967d0

            else if(name(i)(1:2).eq.'C_')then

              amass(i)=12.011d0

            else if(name(i)(1:2).eq.'S_')then

              amass(i)=32.06d0

            else

              write(nprint,'(a)')
     x          'error - unidentified atom type : '//name(i)
              safe=.false.

            endif

          enddo

        endif

        if(.not.safe)then

          close(nsave)
          stop

        endif

        nstr=nstr+1

        if(imcon.gt.0)then

c     convert to reduced coordinates

          do i=1,nclust

            sxx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
            syy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
            szz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
            xxx(i)=sxx
            yyy(i)=syy
            zzz(i)=szz

          enddo

c     unscramble periodic boundary effects

          do i=2,nclust

            xdf=xxx(i)-xxx(i-1)
            ydf=yyy(i)-yyy(i-1)
            zdf=zzz(i)-zzz(i-1)

            xxx(i)=xxx(i-1)+(xdf-nint(xdf))
            yyy(i)=yyy(i-1)+(ydf-nint(ydf))
            zzz(i)=zzz(i-1)+(zdf-nint(zdf))

          enddo

c     restore to natural coordinates

          do i=1,nclust

            sxx=cell(1)*xxx(i)+cell(4)*yyy(i)+cell(7)*zzz(i)
            syy=cell(2)*xxx(i)+cell(5)*yyy(i)+cell(8)*zzz(i)
            szz=cell(3)*xxx(i)+cell(6)*yyy(i)+cell(9)*zzz(i)
            xxx(i)=sxx
            yyy(i)=syy
            zzz(i)=szz

          enddo

        endif

c     calculate structural properties

        pln2=0.d0
        prga=0.d0
        prgb=0.d0
        prgc=0.d0
        prg2=0.d0
        pat1=0.d0
        pat2=0.d0
        pat3=0.d0

c     calculate square of end-to-end length and histogram

        do k=1,narms

          i=length(1,k)
          j=length(2,k)
          xln2=(xxx(i)-xxx(j))**2+(yyy(i)-yyy(j))**2+(zzz(i)-zzz(j))**2
          pln2=pln2+xln2
          m=xln2*dble(mxcor)/range**2
          if(m.le.mxcor)hln2(m)=hln2(m)+1.d0

        enddo
        pln2=pln2/dble(narms)

c     locate centre of mass of cluster

        cmx=0.d0
        cmy=0.d0
        cmz=0.d0
        dmass=0.d0
        do i=1,nclust

          cmx=cmx+amass(i)*xxx(i)
          cmy=cmy+amass(i)*yyy(i)
          cmz=cmz+amass(i)*zzz(i)
          dmass=dmass+amass(i)

        enddo
        cmx=cmx/dmass
        cmy=cmy/dmass
        cmz=cmz/dmass

c     calculate inertia tensor (exclude gold atoms)

        rxx=0.d0
        ryy=0.d0
        rzz=0.d0
        rxy=0.d0
        rxz=0.d0
        ryz=0.d0
        do i=ngold+1,nclust

          rxx=rxx+amass(i)*(xxx(i)-cmx)*(xxx(i)-cmx)
          ryy=ryy+amass(i)*(yyy(i)-cmy)*(yyy(i)-cmy)
          rzz=rzz+amass(i)*(zzz(i)-cmz)*(zzz(i)-cmz)
          ryz=ryz+amass(i)*(yyy(i)-cmy)*(zzz(i)-cmz)
          rxz=rxz+amass(i)*(xxx(i)-cmx)*(zzz(i)-cmz)
          rxy=rxy+amass(i)*(xxx(i)-cmx)*(yyy(i)-cmy)

        enddo
        rxx=rxx/dble(nclust-ngold)
        ryy=ryy/dble(nclust-ngold)
        rzz=rzz/dble(nclust-ngold)
        rxy=rxy/dble(nclust-ngold)
        rxz=rxz/dble(nclust-ngold)
        ryz=ryz/dble(nclust-ngold)

c     calculate principal moments of inertia

        re3=1.d0/3.d0
        re9=1.d0/9.d0
        re54=1.d0/54.d0
        aaa=-(rxx+ryy+rzz)
        bbb=rxx*ryy+rxx*rzz+ryy*rzz-rxy*rxy-rxz*rxz-ryz*ryz
        ccc=rxx*(ryz*ryz-ryy*rzz)+rxy*(rzz*rxy-2.d0*rxz*ryz)
     x    +ryy*rxz*rxz

c     solve for roots using viete's method

        rrr=re54*(2.d0*aaa**3-9.d0*aaa*bbb+27.d0*ccc)
        sqrq=sqrt(re9*(aaa*aaa-3.d0*bbb))
        theta=acos(rrr/sqrq**3)
        sqr2=-2.d0*sqrq
        a13=re3*aaa
        ddd=sqr2*cos(re3*theta)-a13
        eee=sqr2*cos(re3*(theta+2.d0*pi))-a13
        fff=sqr2*cos(re3*(theta+4.d0*pi))-a13

c     arrange in ascending order

        sum=ddd+eee+fff
        ggg=max(ddd,eee)
        hhh=min(ddd,eee)
        ggg=max(ggg,fff)
        hhh=min(hhh,fff)
        prga=hhh
        prgb=sum-hhh-ggg
        prgc=ggg
        pat1=(ddd-eee)**2+(ddd-fff)**2+(eee-fff)**2
        pat2=sum**2
        pat3=pat1/(2.d0*sum**2)
        prg2=prga+prgb+prgc

c     accumulate statistics

        snstr1=dble(nstr-1)/dble(nstr)
        snstr2=1.d0/dble(nstr)
        fln2=snstr1*(fln2+snstr2*(pln2-aln2)**2)
        frga=snstr1*(frga+snstr2*(prga-arga)**2)
        frgb=snstr1*(frgb+snstr2*(prgb-argb)**2)
        frgc=snstr1*(frgc+snstr2*(prgc-argc)**2)
        frg2=snstr1*(frg2+snstr2*(prg2-arg2)**2)
        fat1=snstr1*(fat1+snstr2*(pat1-aat1)**2)
        fat2=snstr1*(fat2+snstr2*(pat2-aat2)**2)
        fat3=snstr1*(fat3+snstr2*(pat3-aat3)**2)
        aln2=snstr1*aln2+snstr2*pln2
        arga=snstr1*arga+snstr2*prga
        argb=snstr1*argb+snstr2*prgb
        argc=snstr1*argc+snstr2*prgc
        arg2=snstr1*arg2+snstr2*prg2
        aat1=snstr1*aat1+snstr2*pat1
        aat2=snstr1*aat2+snstr2*pat2
        aat3=snstr1*aat3+snstr2*pat3

      enddo

  100 continue
      close(nsave)

c     calculate averages of structural data and standard errors

      write(nprint,10)nstr
      ravstr=1.d0/dble(nstr-1)
      fln2=sqrt(ravstr*fln2)
      frga=sqrt(ravstr*frga)
      frgb=sqrt(ravstr*frgb)
      frgc=sqrt(ravstr*frgc)
      frg2=sqrt(ravstr*frg2)
      fat1=sqrt(ravstr*fat1)
      fat2=sqrt(ravstr*fat2)
      fat3=sqrt(ravstr*fat3)
      spha=0.5d0*aat1/aat2
      espa=spha*(fat1/aat1+fat2/aat2)
      write(nprint,20)aln2,fln2,arga,frga,argb,frgb,argc,frgc,
     x                arg2,frg2,spha,espa,aat3,fat3

      end

