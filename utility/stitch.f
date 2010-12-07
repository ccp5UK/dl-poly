      implicit none

      integer, parameter :: lenrec=80
      integer, parameter :: li=9
      integer, parameter :: no=8
      integer, parameter :: in=7
      integer, parameter :: mxatms=100000
      
      logical go,findstring
      character*8 aname
      character*40 CNAME,FNAME
      character*80 text
      integer i,j,m,k,n,nmols,levcfg,imcon,natms,nmat,nlnk,link,ia,ib
      integer nbnds,ncons,intstr,afix,bfix,ibnds,mxbnds
      real*8 cell,xxx,yyy,zzz,sx,sy,sz,rcell,det,rsq,rmx,rmn
      real*8 xdf,ydf,zdf
      
      character*8 name(mxatms)
      dimension cell(9),rcell(9),link(2,mxatms)
      dimension xdf(mxatms),ydf(mxatms),zdf(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension afix(mxatms),bfix(mxatms)

c     declare input files

      write(*,*)'Enter name of FIELD file'
      read(*,*)FNAME

      write(*,*)'Enter name of CONFIG file'
      read(*,*)CNAME

c     search from molecules

      go=.true.
      open(li,file=FNAME)
      do while(go)
         read(li,'(a80)',end=400)text
         if(findstring("numm",text,i).or.
     x      findstring("NUMM",text,i))go=.false.
      enddo
      nmols=intstr(text,lenrec,i)
      write(*,*)'Number of molecules',nmols
      close (li)

c     search for bonds

      nbnds=0
      go=.true.
      open(li,file=FNAME)
      do while(go)
         read(li,'(a80)',end=100)text
         if(findstring("bond",text,i).or.
     x      findstring("BOND",text,i))go=.false.
      enddo
      nbnds=intstr(text,lenrec,i)
 100  write(*,*)'Number of bonds',nbnds
      do i=1,nbnds
         read(li,*)ia,ib
         link(1,i)=min(ia,ib)
         link(2,i)=max(ia,ib)
      enddo
      close (li)

c     search for constraints

      ncons=0
      go=.true.
      open(li,file=FNAME)
      do while(go)
         read(li,'(a80)',end=200)text
         if(findstring("cons",text,i).or.
     x      findstring("CONS",text,i))go=.false.
      enddo
      ncons=intstr(text,lenrec,i)
 200  write(*,*)'Number of constraints',ncons
      do i=1,ncons
         read(li,*)ia,ib
         link(1,i+nbnds)=min(ia,ib)
         link(2,i+nbnds)=max(ia,ib)
      enddo
      close (li)
      nlnk=nbnds+ncons
      if(nlnk.eq.0)then
         write(*,*)'No connection information found'
         stop
      endif
      write(*,*)'Total number of links',nlnk

c     read CONFIG file

      open(in,file=CNAME)
      read(in,'(a80)')text
      write(*,*)text
      read(in,'(2i10)')levcfg,imcon
      if(imcon.gt.0)then
         read(in,'(3f20.12)')cell(1),cell(2),cell(3)
         read(in,'(3f20.12)')cell(4),cell(5),cell(6)
         read(in,'(3f20.12)')cell(7),cell(8),cell(9)
      endif
      do i=1,mxatms
         read(in,'(a8,i10)',end=300)aname,n
         if(n.eq.0)n=i
         name(n)=aname
         read(in,'(3f20.12)')xxx(n),yyy(n),zzz(n)
         if(levcfg.gt.0)read(in,*)
         if(levcfg.gt.1)read(in,*)
      enddo
 300  continue
      natms=i-1
      write(*,*)'atoms in CONFIG file = ',natms
      close(in)

c     construct reduced coordinates

      call invert(cell,rcell,det)
      do m=1,natms
         sx=xxx(m)*rcell(1)+yyy(m)*rcell(4)+zzz(m)*rcell(7)
         sy=xxx(m)*rcell(2)+yyy(m)*rcell(5)+zzz(m)*rcell(8)
         sz=xxx(m)*rcell(3)+yyy(m)*rcell(6)+zzz(m)*rcell(9)
         xxx(m)=sx
         yyy(m)=sy
         zzz(m)=sz
      enddo     

c     set up fix arrays

      do i=1,mxatms
         afix(i)=0
         bfix(i)=0
      enddo

      mxbnds=nmols*nlnk
      if(mxatms.lt.mxbnds)then
         write(*,*)"too many bonds to process"
         stop
      endif

c     define bond vectors

      k=0
      nmat=natms/nmols
      do j=1,nmols
         afix(k+1)=1
         do i=1,nlnk
            m=k+link(1,i)
            n=k+link(2,i)
            xdf(k+i)=xxx(n)-xxx(m)-nint(xxx(n)-xxx(m))
            ydf(k+i)=yyy(n)-yyy(m)-nint(yyy(n)-yyy(m))
            zdf(k+i)=zzz(n)-zzz(m)-nint(zzz(n)-zzz(m))
         enddo
         k=k+nmat
      enddo

c     loop over bond vectors

      ibnds=0
      do while(ibnds.lt.mxbnds)

c     remove periodic boundary condition

         k=0
         do j=1,nmols
            do i=1,nlnk
               if(bfix(k+i).eq.0)then
                  m=k+link(1,i)
                  n=k+link(2,i)
                  if(afix(m).eq.1)then
                     xxx(n)=xxx(m)+xdf(k+i)
                     yyy(n)=yyy(m)+ydf(k+i)
                     zzz(n)=zzz(m)+zdf(k+i)
                     ibnds=ibnds+1
                     bfix(k+i)=1
                     afix(n)=1
                  elseif(afix(n).eq.1)then
                     xxx(m)=xxx(n)-xdf(k+i)
                     yyy(m)=yyy(n)-ydf(k+i)
                     zzz(m)=zzz(n)-zdf(k+i)
                     ibnds=ibnds+1
                     bfix(k+i)=1
                     afix(m)=1
                  endif
               endif
            enddo
            k=k+nmat
         enddo
         write(*,*)'===> ibnds',ibnds
      enddo

c     write new CONFIG file

      open(no,file='CONFIG.NEW')
      write(no,'(a80)')text
      write(no,'(2i10)')0,imcon
      write(no,'(3f20.12)')cell(1),cell(2),cell(3)
      write(no,'(3f20.12)')cell(4),cell(5),cell(6)
      write(no,'(3f20.12)')cell(7),cell(8),cell(9)
      do m=1,natms
         sx=xxx(m)*cell(1)+yyy(m)*cell(4)+zzz(m)*cell(7)
         sy=xxx(m)*cell(2)+yyy(m)*cell(5)+zzz(m)*cell(8)
         sz=xxx(m)*cell(3)+yyy(m)*cell(6)+zzz(m)*cell(9)
         write(no,'(a8)')name(m)
         write(no,'(3f20.12)')sx,sy,sz
         xxx(m)=sx
         yyy(m)=sy
         zzz(m)=sz
      enddo     
      close(no)

c     check for periodic boundary effects

      k=0
      rmx=0.d0
      rmn=1.d9
      do j=1,nmols
         do i=1,nlnk
            m=k+link(1,i)
            n=k+link(2,i)
            rsq=(xxx(m)-xxx(n))**2+(yyy(m)-yyy(n))**2+
     x         (zzz(m)-zzz(n))**2
            rmx=max(rmx,rsq)
            rmn=min(rmn,rsq)
         enddo
         k=k+nmat
      enddo
      write(*,*)'largest  separation ',sqrt(rmx)
      write(*,*)'smallest separation ',sqrt(rmn)

      stop

c     abort file readers

 400  continue

      write(*,*)'Error end of file encountered'
      stop

      end
      subroutine invert(a,b,d)

c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c     wl
c     2002/05/31 14:52:39
c     1.1
c     Exp
c     
c***********************************************************************
      
      implicit none

      real*8 a,b,d,r

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
      end
      function intstr(word,len,lst)

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
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical flag,count,final
      character*1 n,word,ksn
      integer lst,len,intstr,j,isn

      dimension n(0:9),word(len)
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      lst=0
      ksn='+'
      intstr=0
      flag=.false.
      final=.false.
      count=.false.
      
      do while(lst.lt.len.and.(.not.final))

        lst=lst+1
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            count=.true.
            flag=.true.
            
          endif
          
        enddo

        if(count.and.(.not.flag))final=.true.
        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

      enddo

      intstr=isn*intstr

      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      return
      
      end
      logical function findstring(seek,string,here)

c***********************************************************************
c     
c     DL_POLY routine to find an explicit string in an input record
c     note: variable `seek' is a character string while variable
c    `string' is a character*1 array i.e. code is application specific
c
c     copyright daresbury laboratory
c     author    w.smith   jan   2004
c     
c     itt
c     2010-10-30 17:20:50
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      integer, parameter :: lenrec=80

      integer i,n,m,here
      character*(*) seek
      character*1 string(lenrec)
      character*1 record(lenrec)

      m=lenrec
      n=len(seek)
      findstring=.false.

      here=0
      do while(here.le.m-n.and.(.not.findstring))

        findstring=.true.

        do i=1,n
          if(seek(i:i).ne.string(here+i))findstring=.false.
        enddo

        here=here+1

      enddo

      return
      end
