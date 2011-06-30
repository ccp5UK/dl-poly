      subroutine sort(n,index,x)

c***********************************************************************
c
c     DL_POLY utility
c     shell sort of index array
c
c     copyright daresbury laboratory 1993
c     author -        t.forester may 1993
c
c*************************************************************************

      implicit real*8 (a-h,o-z)
      dimension x(*)
      dimension index(*)
c
c     set up index array

      do i = 1,n
         index(i) = i
      enddo

c
c     set up sort

      if(n.gt.1) then

c     number of lists
         nl = n/2

c     iterate shell sort

   10  do nn = 1,nl
c
c     begin insertion sort on nnth list

            do i = nn+nl,n,nl

               xmax = x(index(i))
               ix = i
               ixind = index(i)
c
c     find location for insertion

               do j = i-nl,1,-nl

                  if (x(index(j)).gt.xmax) then

                     ix = j

                  else

                     j = 1

                  endif

               enddo

c
c     insert in index array

               do j = i,ix+nl,-nl

                  index(j) = index(j-nl)

               enddo

               index(ix) = ixind

            enddo

         enddo


         nl = nl/2
         if(nl.gt.0) goto 10

      endif
      return
      end
