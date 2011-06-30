      subroutine bragg(natms,kmax,cell,xxx,yyy,zzz)

c*********************************************************************
c
c     dl_poly utility routine to calculate the bragg planes in a
c     CONFIG file.
c
c     copyright daresbury laboratory
c     author - w.smith oct 2009
c
c*********************************************************************

      implicit none

      integer i,k,l,m,n,ll,mm,nn,kk,fail,nmin,mmin,natms,kmax,kmax3
      real(8) twopi,peak,stx,sty,stz,det
      real(8) cell(9),xxx(natms),yyy(natms),zzz(natms),rcell(9)

      real(8), allocatable :: ckc(:),cks(:),clm(:),slm(:)
      real(8), allocatable :: emc(:,:),enc(:,:),els(:,:)
      real(8), allocatable :: ems(:,:),ens(:,:),elc(:,:)
      real(8), allocatable :: strk(:)

      data twopi/6.2831853072d0/

      kmax3=((kmax+1)**3)/2
      allocate(ckc(natms),cks(natms),clm(natms),slm(natms),stat=fail)
      allocate(emc(natms,0:kmax),enc(natms,0:kmax),stat=fail)
      allocate(els(natms,0:kmax),ems(natms,0:kmax),stat=fail)
      allocate(ens(natms,0:kmax),elc(natms,0:kmax),stat=fail)
      allocate(strk(kmax3),stat=fail)

      write(*,"(1x,'bragg plane projection',//,
     x  1x,' l m n   str.fac')")

      peak=1.d0
      open(8,file="BRAGG")

c     invert cell matrix

      call invert(cell,rcell,det)

c     calculate and store exponential factors

      do i=1,natms

        elc(i,0)=1.0d0
        emc(i,0)=1.0d0
        enc(i,0)=1.0d0
        els(i,0)=0.0d0
        ems(i,0)=0.0d0
        ens(i,0)=0.0d0

        stx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        sty=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        stz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

        elc(i,1)=cos(twopi*stx)
        emc(i,1)=cos(twopi*sty)
        enc(i,1)=cos(twopi*stz)
        els(i,1)=sin(twopi*stx)
        ems(i,1)=sin(twopi*sty)
        ens(i,1)=sin(twopi*stz)

      enddo

      do k=2,kmax

        do i=1,natms

          elc(i,k)=elc(i,k-1)*elc(i,1)-els(i,k-1)*els(i,1)
          emc(i,k)=emc(i,k-1)*emc(i,1)-ems(i,k-1)*ems(i,1)
          enc(i,k)=enc(i,k-1)*enc(i,1)-ens(i,k-1)*ens(i,1)
          els(i,k)=els(i,k-1)*elc(i,1)+elc(i,k-1)*els(i,1)
          ems(i,k)=ems(i,k-1)*emc(i,1)+emc(i,k-1)*ems(i,1)
          ens(i,k)=ens(i,k-1)*enc(i,1)+enc(i,k-1)*ens(i,1)

        enddo

      enddo

c     loop over all k vectors  k=2pi(ll/cl,mm/cl,nn/cl)

      kk=0
      mmin=0
      nmin=1

      do l=0,klim
        do mm=mmin,kmax
          m=iabs(mm)
          if(mm.ge.0)then
            do i=1,natms
              clm(i)=elc(i,l)*emc(i,m)-els(i,l)*ems(i,m)
              slm(i)=els(i,l)*emc(i,m)+ems(i,m)*elc(i,l)
            enddo
          else
            do i=1,natms
              clm(i)=elc(i,l)*emc(i,m)+els(i,l)*ems(i,m)
              slm(i)=els(i,l)*emc(i,m)-ems(i,m)*elc(i,l)
            enddo
          endif

          do nn=nmin,kmax
            n=iabs(nn)+1

            if(nn.ge.0)then
              do i=1,natms
                ckc(i)=clm(i)*enc(i,n)-slm(i)*ens(i,n)
                cks(i)=slm(i)*enc(i,n)+clm(i)*ens(i,n)
              enddo
            else
              do i=1,natms
                ckc(i)=clm(i)*enc(i,n)+slm(i)*ens(i,n)
                cks(i)=slm(i)*enc(i,n)-clm(i)*ens(i,n)
              enddo
            endif

            ckcs=0.0d0
            ckss=0.0d0
            do i=1,natm
              ckcs=ckcs+ckc(i)
              ckss=ckss+cks(i)
            enddo
            kk=kk+1
            strk(kk)=(ckcs*ckcs+ckss*ckss)/dble(natms)
            if(strk(kk).gt.peak)
     x        write(*,'(1x,3i4,f10.4)')ll,mm,nn,strk(kk)
            write(8,'(1x,3i4,f10.4)')ll,mm,nn,strk(kk)

          enddo

          nmin=-kmax

        enddo

        mmin=-kmax

      enddo

      deallocate(ckc,cks,clm,slm,strk,stat=fail)
      deallocate(emc,enc,els,ems,ens,elc,stat=fail)
      close(8)

      return
      end
      subroutine invert(a,b,d)

c***********************************************************************
c
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

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

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

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
      end subroutine invert

