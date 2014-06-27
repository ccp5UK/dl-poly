       PROGRAM PMFtoTABLE
c c c  
c c c  The program creates a tabulated force-field file in DL_POLY format - TABLE
c c c  by converting a *PMF or *TAB data file (e.g. BNDPMF/BNDTAB, ANGPMF/ANGTAB etc.)
c c c  
       implicit double precision (a-h,o-z)
c
       PARAMETER( finp=21,fout=31,Mgrid=20000,Msets=100 )
c
       double precision u(Mgrid),f(Mgrid)
       double precision rcut,delr,delta,r,rrr
       integer      iflag,i,k,n,ngrd,nset,nrgd4
       character*80 title,header
       character    fname*6,ch2*2
c c c  
!       write(*,*)
!     * 'The name of force-field file (*PMF/*TAB) to convert [6 chars]:'
!       read(*,*) fname
c
!       open(unit=finp,file=fname,ERR=13)
c c c  
c c c  using standard input stream (terminal=stdin)
c c c  read stdin from a file by supplying the file name preceded by "<"
c c c  
       iflag = 0
       read(*,'(a2,a80)',ERR=13,END=13) ch2,title
       read(*,*,ERR=13,END=13) ch2,rcut,ngrd,delr,nset
       iflag = 1
c
 13    continue
c
       if(iflag.lt.1) then
         write(*,*)'*!* ERROR: reading input - FULL STOP *!*'
         STOP
       elseif( mod(ngrd,4).gt.0 ) then
         write(*,*)'*!* ERROR: input Ngrid must be a multiple of 4',
     *             ' - FULL STOP *!*'
         STOP
       elseif( rcut.le.0.0 ) then
         write(*,*)'*!* ERROR: input Rcut must be greater than 0.0',
     *             ' - FULL STOP *!*'
         STOP
       endif
c
       delta = rcut/float(ngrd)
c
       if( abs(delta-delr).gt.1.0e-6 ) then
         write(*,*)'*!* ERROR: input delr inconsistent with Ngrid',
     *             ' - FULL STOP *!*'
         STOP
       endif
c
       ngrd4 = ngrd/4
c
       iflag = -2
c
       open(unit=fout,file='TABLE',ERR=113)
c
       write(fout,'(a80)',ERR=113) title
       write(fout,'(2f10.5,i10)',ERR=113)delta,rcut,ngrd+4
c
       iflag = 1
c
       DO n=1,nset

         iflag = -1
         read(*,*,ERR=113,END=113) 
         iflag = 0
         read(*,'(a2,a16)',ERR=113,END=113) ch2,header
         iflag = 1

         Do i=1,ngrd
           r = float(i)*delta

           read(*,*,ERR=113,END=113) rrr,u(i),f(i)
c
!         if(u(i).gt.1.d4) then
!            u(i)=1.d4
!         endif
!         if(u(i).lt.-1.d4) then
!            u(i)=-1.d4
!         endif
!         if(f(i).gt.1.d4) then
!            f(i)=1.d4
!         endif
!         if(f(i).lt.-1.d4) then
!            f(i)=-1.d4
!         endif
         EndDo
c
         write(fout,'(a16)') header
c
         do k=1,ngrd4
            write(fout,1001) u(4*k-3),u(4*k-2),u(4*k-1),u(4*k)
         enddo
         write(fout,1001) 0.0,0.0,0.0,0.0
c
         do k=1,ngrd4
            write(fout,1001) f(4*k-3),f(4*k-2),f(4*k-1),f(4*k)
         enddo
         write(fout,1001) f(ngrd),f(ngrd),f(ngrd),f(ngrd)

 1001    format(4(e15.7))

       ENDDO
c
       close(fout)
c
 113   continue
c
       if(iflag.lt.-1) then
         write(*,*)'*!* ERROR: opening output - FULL STOP *!*'
         STOP
       elseif(iflag.lt.1) then
         write(*,*)'*!* ERROR: reading input - FULL STOP *!* ',n,i
         STOP
       elseif( i.lt.ngrd+1 ) then
         write(*,*)'*!* ERROR: input too short - FULL STOP *!* ',n,i
         STOP
       endif
c
       STOP
       END
