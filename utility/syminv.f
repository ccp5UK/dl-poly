      subroutine syminv(nnn,aaa,bbb,bbi)
c
c***********************************************************************
c
c     routine to invert a real symmetric matrix (reference:
c     computing methods v.ii, i.s. berezin and n.p.zhidkov,
c     pergamon press 1965). note that the matrices are
c     packed in the minimum storage mode, (i.e. a(k)=a(i,j),
c     where i>j and k=(i*(i-1))/2+j ).
c     the matrices aaa and bbb may be equivalenced though
c     this will destroy the contents of the original
c     array.
c
c     general version for all real symmetric matrices
c
c     copyright - daresbury laboratory 1993
c     author w.smith (added to dl_poly july 1993)
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      dimension aaa(*),bbb(*),bbi(*)
      ind(i,j)=(max0(i,j)*(max0(i,j)-1))/2+min0(i,j)
c
c     factorize matrix aaa into upper and lower forms
      do l=1,nnn
        do m=1,l
          k=ind(l,m)
          bbi(k)=0.d0
          bbb(k)=aaa(k)
          if(m.ne.1)then
            lm=m-1
            do n=1,lm
              ln=ind(l,n)
              mn=ind(m,n)
              temv  =bbb(k)-bbb(ln)*bbb(mn)+bbi(ln)*bbi(mn)
              bbi(k)=bbi(k)-bbi(ln)*bbb(mn)-bbb(ln)*bbi(mn)
              bbb(k)=temv
            enddo
          endif
          if(l.eq.m)then
            if(bbb(k).lt.0.d0)then
              bbi(k)=sqrt(abs(bbb(k)))
              bbb(k)=0.d0
            else
              bbb(k)=sqrt(abs(bbb(k)))
              bbi(k)=0.d0
            endif
          else
            mm=ind(m,m)
            den=1.d0/(bbb(mm)**2+bbi(mm)**2)
            temv  =den*(bbb(k)*bbb(mm)+bbi(k)*bbi(mm))
            bbi(k)=den*(bbi(k)*bbb(mm)-bbb(k)*bbi(mm))
            bbb(k)=temv
          endif
        enddo
      enddo
c
c     invert lower triangular matrix
      do l=1,nnn
        n=ind(l,l)
        den=1.d0/(bbb(n)*bbb(n)+bbi(n)*bbi(n))
        do m=1,l
          k=ind(l,m)
          if(l.eq.m)then
            bbb(k)= den*bbb(k)
            bbi(k)=-den*bbi(k)
          else
            lm=l-1
            suma=0.d0
            sumb=0.d0
            do j=m,lm
              lj=ind(l,j)
              mj=ind(m,j)
              suma=suma-bbb(lj)*bbb(mj)+bbi(lj)*bbi(mj)
              sumb=sumb-bbb(lj)*bbi(mj)-bbi(lj)*bbb(mj)
            enddo
            temv  =den*(suma*bbb(n)+sumb*bbi(n))
            bbi(k)=den*(sumb*bbb(n)-suma*bbi(n))
            bbb(k)=temv
          endif
        enddo
      enddo
c
c     form product of upper and lower inverse triangular matrices
      do l=1,nnn
        do m=1,l
          sum=0.d0
          k=ind(l,m)
          do j=l,nnn
            lj=ind(l,j)
            mj=ind(m,j)
            sum=sum+bbb(lj)*bbb(mj)-bbi(lj)*bbi(mj)
          enddo
          bbb(k)=sum
        enddo
      enddo
      return
      end
