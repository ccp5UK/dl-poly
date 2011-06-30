      subroutine solve(ndim,lpv,aaa,xxx,bbb)
c
c***********************************************************************
c
c     gaussian elimination routine with pivoting
c
c     copyright daresbury laboratory 1987
c     w. smith aug. 1987
c
c***********************************************************************
c

      implicit real*8(a-h,o-z)

      dimension aaa(ndim,ndim),xxx(ndim),bbb(ndim),lpv(ndim)

c
c     initialise pivot arrays

      do j=1,ndim

         lpv(j)=j

      enddo

c
c     construct upper triangular matrix
c

      do  i=1,ndim-1

c
c     locate pivot elements

         lpc=i
         lpr=i
         piv=aaa(i,i)

         do j=i,ndim

            do k=i,ndim

               if(abs(aaa(j,k)).gt.abs(piv))then
                  piv=aaa(j,k)
                  lpc=k
                  lpr=j
               endif

            enddo

         enddo
c
c     perform row pivots

         saveit=bbb(lpr)
         bbb(lpr)=bbb(i)
         bbb(i)=saveit

         do j=i,ndim

            saveit=aaa(lpr,j)
            aaa(lpr,j)=aaa(i,j)
            aaa(i,j)=saveit

         enddo

c
c     perform column pivots

         do j=1,ndim

            saveit=aaa(j,lpc)
            aaa(j,lpc)=aaa(j,i)
            aaa(j,i)=saveit

         enddo

c
c     remove pivot columns

         do j=i+1,ndim

            bbb(j)=bbb(j)-(aaa(j,i)/piv)*bbb(i)

         enddo

         do j=i+1,ndim

            do k=ndim,i,-1

               aaa(j,k)=aaa(j,k)-(aaa(j,i)/piv)*aaa(i,k)

            enddo

         enddo

c
c     record pivoting data

         keepit=lpv(lpc)
         lpv(lpc)=lpv(i)
         lpv(i)=keepit

c
c     end of upper triangulation

      enddo

c
c     solve for unknowns by back substitution
      do i=ndim,2,-1

         bbb(i)=bbb(i)/aaa(i,i)

         do j=i-1,1,-1

            bbb(j)=bbb(j)-bbb(i)*aaa(j,i)

         enddo

      enddo

      bbb(1)=bbb(1)/aaa(1,1)

c
c     re-order solution vectors according to pivot

      do i=1,ndim

         xxx(lpv(i))=bbb(i)

      enddo

      return
      end
