      program smooth

c***********************************************************************
c
c     dl_poly utility for smoothing data in "xvgr" type file
c     format  x,y data sets separated by &
c     lines beginning "#" = comments
c     lines beginning "@TYPE xy" => x,y data follow on following lines
c
c     smoothing algorithm from Allen and Tildesley, p204
c     3rd degree, five point smooth.
c
c     copyright daresbury laboratory 1996
c     author t.forester     March 1996
c
c***********************************************************************

      implicit real*8(a-h,o-z)
      parameter (mxn = 10000)
      dimension x(mxn),y(mxn),ys(mxn)
      character*80 oneline

      write(*,*) ' name of input file '
      read(*,*) oneline

      open(10,file=oneline,form='formatted',status='old')

      write(*,*) ' name of output file '
      read(*,*) oneline

      open(11,file=oneline,form='formatted')

      write(*,*) 'number of smooth iterations ?'
      read(*,*) nsmooth
      nsmooth = max(1,nsmooth)

  101 read(10,'(a80)',end=100) oneline

      if(oneline(1:8).eq.'@TYPE xy') then

        write(11,'(a8)') oneline

        i = 1
   10   read(10,*,end=20,err=20) x(i),y(i)
        i=i+1
        goto 10

   20   nn = i-1

        do ismooth = 1,nsmooth

        if(nn.gt.5) then

          ys(1) = 1/7.d1*(69.d0*y(1)+4.d0*y(2)-6.d0*y(3)+4.d0*y(4)
     x       -y(5))

          ys(2) = 1/35.d0*(2.d0*y(1)+27.d0*y(2)+12.d0*y(3)-8.d0*y(4)
     x       +2.d0*y(5))

          do i = 3,nn-2

            ys(i) = 1.d0/35.d0*(-3.d0*y(i-2)+12.d0*y(i-1)+17.d0*y(i)
     x         +12.d0*y(i+1)-3.d0*y(i+2))


          enddo

          ys(nn) = 1/7.d1*(69.d0*y(nn)+4.d0*y(nn-1)-6.d0*y(nn-2)
     x       +4.d0*y(nn-3) -y(nn-4))

          ys(nn-1) = 1/35.d0*(2.d0*y(nn)+27.d0*y(nn-1)+12.d0*y(nn-2)
     x       -8.d0*y(nn-3)+2.d0*y(nn-4))

          do ii = 1,nn
            y(ii) = max(ys(ii),0.d0)
          enddo

        endif

        enddo

c
c     write out smoothed results

        do i=1,nn

          write(11,'(f8.4,1p,2e15.6)') x(i),y(i)

        enddo
        write(11,*) '&'

      else

        write(11,'(a80)') oneline
      endif

      goto 101

  100 continue
      end

