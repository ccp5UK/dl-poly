      subroutine shellsort(n,lst,aaa)

c********************************************************************
c
c     basic shellsort routine
c
c     author w.smith oct 2007
c     daresbury laboratory
c
c*********************************************************************

      implicit none

      real(8) tmp
      integer i,j,k,l,m,n,kk
      real(8) aaa(*)
      integer lst(*)

      m=n

      do while(m.gt.0)

        m=m/2
        k=n-m

        do j=1,k

          i=j

          do while(i.gt.0)

            l=i+m

            if(aaa(l).lt.aaa(i))then

              kk=lst(i)
              lst(i)=lst(l)
              lst(l)=kk
              tmp=aaa(i)
              aaa(i)=aaa(l)
              aaa(l)=tmp
              i=i-m

            else

              i=0

            endif

          enddo

        enddo

      enddo

      return
      end
