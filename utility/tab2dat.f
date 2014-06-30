       PROGRAM TABLEtoDATA
c*********************************************************************
c     
c     dl_poly program to convert a tabulated force-field file in DL_POLY format - TABLE
c     into a three column text-data file VDWTAB (readable by [Xm]Grace and other plotters)
c     
c     copyright daresbury laboratory 2014
c     author    Andrey Brukhno, June 2014
c
c     compilation: f77 -o tab2dat.exe tab2dat.f
c     usage: tab2dat.exe < TABLE
c     
c*********************************************************************
c     
       implicit double precision (a-h,o-z)
c
       PARAMETER( finp=21,fout=31,Mgrid=20000,Msets=1000 )
c
       double precision u(Mgrid),f(Mgrid)
       double precision rcut,delr,delta,r,rrr
       integer      iflag,i,k,n,ngrd,nset,nrgd4
       character*80 title,header
       character    fname*6,ch2*2
c c c  
c c c  using standard input stream (terminal=stdin)
c c c  
       nset  = -1
       ngrd  = 0
       delr  = 0.0
       rcut  = 0.0
c
       iflag = 0
c c c
c c c  NOTE: introduce number of sets to read at the end of the second line!
c c c
       read(*,'(a80)',ERR=13,END=13) title
       read(*,*,ERR=13,END=13) delr,rcut,ngrd,nset
       iflag = 1
c
 13    continue
c
       if(iflag.lt.1) then
           write(*,*)'*!* ERROR: reading input - FULL STOP *!*'
           write(*,*)'*!* Must read: delr, rcut, Ngrid, Nsets *!*'
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
       ngrd4 = ngrd/4
       ngrd  = ngrd-4
       delta = rcut/float(ngrd)
c
       if( abs(delta-delr).gt.1.0e-6 ) then
         write(*,*)'*!* ERROR: input delr inconsistent with Ngrid',
     *             ' - FULL STOP *!*'
         STOP
       endif
c
       iflag = -2
c
       open(unit=fout,file='VDWTAB',ERR=113)
c
       write(fout,'(a2,a80)') '# ',title
       write(fout,'(a2,2(f10.5,i10),a)')'# ',rcut,ngrd,delta,nset,
     * '   (rcut, Ngrid, delr, Nsets)'
c
       iflag = 1
c
       DO n=1,nset

         iflag = -1

         read(*,'(a16)',ERR=113,END=113) header

         iflag = 0

         Do i=1,ngrd4
           read(*,*,ERR=113,END=113) (u((i-1)*4+k),k=1,4)
         EndDo

         Do i=1,ngrd4
           read(*,*,ERR=113,END=113) (f((i-1)*4+k),k=1,4)
         EndDo

         iflag = 1

         write(fout,*)
         write(fout,'(a2,a16)') '# ',header

         Do k=1,ngrd
           rrr = float(k)*delta
           write(fout,'(f10.5,2e15.7)') rrr,u(k),f(k)
         EndDo

       ENDDO
c
       close(fout)
c
 113   continue
c
       if(iflag.lt.-1) then
         write(*,*)'*!* ERROR: opening output - FULL STOP *!*'
         STOP
       elseif(iflag.eq.0) then
         write(*,*)'*!* ERROR: reading input - FULL STOP *!* ',n,i
         STOP
       elseif( i.lt.ngrd4+1 ) then
         write(*,*)'*!* ERROR: input too short - FULL STOP *!* ',n,i
         STOP
       endif
c
       STOP
       END
