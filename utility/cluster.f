      program cluster

c*********************************************************************
c
c     dl_poly program to analyse clusters of selected atoms in
c     a dl_poly configuration file
c
c     copyright daresbury laboratory may 1995
c     author w.smith may 1995
c
c*********************************************************************

      implicit real*8(a-h,o-z)

      parameter (mxatms=8640)

      character*8 atname,target
      character*40 fname

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9),rcell(9)
      dimension idc(mxatms),ist(mxatms),lst(mxatms)

c     enter control variables

      write(*,'(a)')'Enter name of config file'
      read(*,*)fname
      write(*,'(a40)')fname
      write(*,'(a)')'Enter target atom identity'
      read(*,*)target
      write(*,'(a)')target
      write(*,'(a)')'Enter number of config atoms'
      read(*,*)natms
      write(*,'(i6)')natms
      write(*,'(a)')'Enter cutoff radius'
      read(*,*)rcut
      write(*,'(f10.5)')rcut

c     initialise cluster array

      do i=1,mxatms

         idc(i)=0
         ist(i)=0
         lst(i)=0

      enddo

c     read the config file

      open(7,file=fname)

      read(7,*)
      read(7,*)levcfg,imcon
      write(*,'(2i6)')levcfg,imcon

      if(imcon.gt.0)then

         read(7,*)cell(1),cell(2),cell(3)
         read(7,*)cell(4),cell(5),cell(6)
         read(7,*)cell(7),cell(8),cell(9)
         call invert(cell,rcell,det)
         write(*,'(3f12.6)')cell(1),cell(5),cell(9)

      endif

      j=0
      do i=1,natms

         read(7,*)atname
         read(7,*)xx,yy,zz
         if(atname.eq.target)then
            j=j+1
            if(j.gt.mxatms)then

               write(*,'(a)')'error - too many atoms in config file'
               stop

            endif
            xxx(j)=xx
            yyy(j)=yy
            zzz(j)=zz

         endif

         if(levcfg.gt.0)read(7,*)
         if(levcfg.gt.1)read(7,*)

      enddo
      natms=j
      write(*,'(a,i6)')'number of target atoms found: ',natms

c     cluster determination

      ncl=0

      do i=1,natms

         icl=idc(i)
         if(idc(i).eq.0)then

            icl=i
            idc(i)=i
            ncl=ncl+1

         endif

         do j=i+1,natms

            xxd=xxx(i)-xxx(j)
            yyd=yyy(i)-yyy(j)
            zzd=zzz(i)-zzz(j)

            if(imcon.gt.0)then

               xdd=rcell(1)*xxd+rcell(4)*yyd+rcell(7)*zzd
               ydd=rcell(2)*xxd+rcell(5)*yyd+rcell(8)*zzd
               zdd=rcell(3)*xxd+rcell(6)*yyd+rcell(9)*zzd
               sxd=xdd-nint(xdd)
               syd=ydd-nint(ydd)
               szd=zdd-nint(zdd)
               xxd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
               yyd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
               zzd=cell(3)*sxd+cell(6)*syd+cell(9)*szd

            endif

            rsq=xxd**2+yyd**2+zzd**2

            if(rsq.le.rcut**2)then

               if(idc(j).eq.0)then

                  idc(j)=icl

               elseif (idc(j).ne.icl) then

                  jcl=min(icl,idc(j))
                  kcl=max(icl,idc(j))

                  do k=1,natms

                     if(idc(k).eq.kcl)idc(k)=jcl

                  enddo

                  ncl=ncl-1
                  icl=jcl

               endif

            endif

         enddo

      enddo

      write(*,'(a,i5)')'number of clusters found: ',ncl
      close (7)

c     gather cluster statistics

      nlim=0
      do i=1,natms

c     write(*,'(2i6)')i,idc(i)
         ist(idc(i))=ist(idc(i))+1
         nlim=max(nlim,idc(i))

      enddo

      nmax=0
      do i=1,nlim

c     write(*,'(2i6)')i,ist(i)
         nmax=max(nmax,ist(i))
         lst(ist(i))=lst(ist(i))+1

      enddo
      write(*,'(a,i6)')'largest cluster found: ',nmax

      nsum=0
      do i=1,nmax

         if(lst(i).gt.0)then

            nsum=nsum+i*lst(i)
            write(*,'(2i5)')i,lst(i)

         endif

      enddo
      write(*,'(a,i6)')'total number of atoms: ',nsum

      end

      subroutine invert(a,b,d)
c
c***********************************************************************
c
c     routine to invert a 3 * 3 matrix using cofactors
c
c     copyright daresbury laboratory
c     author w.smith
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension a(9),b(9)
c
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
c
c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d
c
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
