Module adf

!Module to calculate the angular distribution functions during a DLPOLY run
!
! author  Oliver Dicks and Aaron Diver
! date    22/02/2018

  Use kinds_f90, only : wp
  Use comms_module, Only : status,dlp_comm_world,MPI_INTEGER,MPI_CHARACTER
  Implicit None
  type, public :: adf_type
    Integer :: myarray_len,coordinterval,coordstart!,ncoordatoms
    real(wp), allocatable :: arraycuts(:)
    character( Len = 8 ), allocatable :: arraypairs(:,:)!,coordatoms(:)
    integer,allocatable :: ltype(:,:),cstat(:,:)
  contains
    procedure :: init=>init_adf
    final :: clean_adf
  end type adf_type
  type, public :: adf_type2
    Integer :: coordinterval,coordstart
    Logical :: coordon
  end type adf_type2
  private
  public :: adf_hello,init_coord_list,check_coord

contains

  subroutine init_adf(T)
    class(adf_type) :: T
    allocate(T%arraycuts(1:T%myarray_len))
    allocate(T%arraypairs(1:T%myarray_len,1:2))
    allocate(T%ltype(1:T%myarray_len,1:2))
    allocate(T%cstat(-3:1000,1:(2*T%myarray_len)))
 !   allocate(T%coordatoms(0:2*T%myarray_len))
  end subroutine init_adf

  subroutine clean_adf(T)
    type(adf_type) :: T

    if (allocated(T%arraycuts)) then
      deallocate(T%arraycuts)
    end if
  end subroutine clean_adf

  subroutine init_coord_list(x,y,z,n,mxlist,mxatdm,list,ltg,coordlist,adfv,adfw,idnode,mxnode,ltype,unqatm,nstep)
    use config_module, only: cfgname ! these two use we do not party like in 80s anymore
    Use setup_module,  only: nicrdt
    Real(kind=wp) , Intent(In) :: x(:)
    Real(kind=wp) , Intent(In) :: y(:)
    Real(kind=wp) , Intent(In) :: z(:)
    Integer, Intent(Out) :: coordlist(0:mxlist,1:mxatdm)
    Integer, Intent(In) :: list(-3:mxlist,1:mxatdm)
    Integer , Intent(In) :: n,mxlist,mxatdm
    Integer, Intent(In) :: ltg(:),ltype(:)
    Integer, Intent(In) :: idnode,mxnode,nstep
    Type(adf_type), Intent(InOut) :: adfv
    Type(adf_type2),Intent(In) :: adfw
    Character(len=8), Intent(In) :: unqatm(:) 
    integer :: i,ii,j,jj,k,kk, ncoord,en,ierr,m,mcoord,lgcoord,nmax,ncb !ncb size coordbuff
    real :: rcut,rab
    Character(len=8) :: aux 
    Character(len=1000),allocatable :: cbuff(:)
    integer, allocatable :: buff(:),coordbuff(:)
    Logical :: newatom
    If(adfw%coordon .Eqv. .False.)Return
    If(adfv%myarray_len==0)Return
    If(mod(nstep,adfw%coordinterval).NE.0)Return
    coordlist(0,:)=0
    ncb=0
    do j = 1, n
      ncoord = 0
      k=list(0,j)
      Do i=1,k
        kk=list(i,j)
        rab = (x(j)-x(kk))**2+(y(j)-y(kk))**2+(z(j)-z(kk))**2
        rcut=0.00

        Do ii= 1 , adfv%myarray_len
          if (((ltype(j)==adfv%ltype(ii,1)) .and. (ltype(kk)==adfv%ltype(ii,2)))&
            .or.&
            ((ltype(j)==adfv%ltype(ii,2)) .and. (ltype(kk)==adfv%ltype(ii,1)))) Then
            rcut=adfv%arraycuts(ii)*adfv%arraycuts(ii)
          endif
        end Do


        if (rab <= rcut) Then
          coordlist(0,j)=coordlist(0,j)+1
          coordlist(coordlist(0,j),j)=kk
          if (kk<=n) Then
            coordlist(0,kk)=coordlist(0,kk)+1
            coordlist(coordlist(0,kk),kk)=j
          endif
        End if
      End Do
    End Do
    adfv%cstat(-3:,:)=0
    Do i=1,adfv%myarray_len
      adfv%cstat(-3,(2*i)-1)=adfv%ltype(i,1)
      adfv%cstat(-2,(2*i)-1)=adfv%ltype(i,2)
      adfv%cstat(-3,(2*i))=adfv%ltype(i,2)
      adfv%cstat(-2,(2*i))=adfv%ltype(i,1)
    End Do

    
    Do i=1,n
      Do j=1,2*adfv%myarray_len
      mcoord=0
        if ((ltype(i)) == adfv%cstat(-3,j)) then
          Do k=1,coordlist(0,i)
            if (ltype(coordlist(k,i)) == adfv%cstat(-2,j)) then
              mcoord=mcoord+1
            end if
          End do
        adfv%cstat(mcoord,j)=adfv%cstat(mcoord,j)+1
        if (mcoord>adfv%cstat(-1,j)) then
          adfv%cstat(-1,j)=mcoord
        End if
        End if
      End do 
    End do

!Set coordbuff size
    Do i=1,2*adfv%myarray_len
      ncb=ncb+(adfv%cstat(-1,i)+4)
    End do

    write(*,*),"NCB:",ncb

    allocate(buff(2*mxatdm))
    allocate(cbuff(mxatdm))
    If (idnode==0) Then
      Open(Unit=nicrdt, File='ICOORD', Form='formatted')
      Write(Unit=nicrdt, Fmt='(a72)') cfgname(1:72)
      Write(Unit=nicrdt, Fmt='(a72)')'Initial coordination between atoms'
!write(*,*)adfv%arraycuts(1:adfv%myarray_len),adfv%arraypairs(1:adfv%myarray_len,1:2)
      Do i=1,n
        m=coordlist(0,i)
        write (nicrdt,Fmt='(i12,1x,i12,1x)',advance="no") &
          ltg(i),coordlist(0,i)
        do ii=1,m
          write(nicrdt,'(i0,1x)',advance="no") ltg(coordlist(ii,i))
        enddo 
        do ii=1,m
          write(nicrdt,'(a,1x)',advance="no") trim(unqatm(ltype(coordlist(ii,i))))
        enddo 
        write(nicrdt,*)
      enddo


      do j=1,mxnode-1
        Call MPI_RECV(en,1,MPI_INTEGER,j,j,dlp_comm_world,status,ierr)
        if (en>0) Then
          Call MPI_RECV(buff,en,MPI_INTEGER,j,j,dlp_comm_world,status,ierr)
          Call MPI_RECV(cbuff,1000*en/2,MPI_CHARACTER,j,j,dlp_comm_world,status,ierr)
          Do i=1,en/2
            write (nicrdt,Fmt='(i12,1x,i12,a)') &
              buff(2*i-1),buff(2*i),trim(cbuff(i))
          enddo
        endif
      enddo
!      Close(Unit=nicrdt)
    else
      en=2*n
      do i=1,n
        buff(2*i-1) = ltg(i)
        buff(2*i) = coordlist(0,i)
        cbuff(i)=''
        do ii=1,coordlist(0,i)
          write(aux,'(i0)') ltg(coordlist(ii,i))
          cbuff(i)=trim(cbuff(i))//" "//trim(aux)
        enddo 
        do ii=1,coordlist(0,i)
          write(aux,'(a)') trim(unqatm(ltype(coordlist(ii,i))))
          cbuff(i)=trim(cbuff(i))//" "//trim(aux)
        enddo 
      enddo
      Call MPI_SEND(en,1,MPI_INTEGER,0,idnode,dlp_comm_world,ierr)
      if (en>0) then
        Call MPI_SEND(buff,en,MPI_INTEGER,0,idnode,dlp_comm_world,ierr)
        Call MPI_SEND(cbuff,1000*n,MPI_CHARACTER,0,idnode,dlp_comm_world,ierr)
      endif
    endif
    deallocate(buff)
    deallocate(cbuff)


    allocate(coordbuff(2*ncb))

    If (idnode==0) Then
      Open(Unit=nicrdt, File='ICOORD', Form='formatted')

      do j=1,mxnode-1
        Call MPI_RECV(ncb,1,MPI_INTEGER,j,j,dlp_comm_world,status,ierr)
          write(*,*),"NCB:",ncb
        if (ncb>1) then
          Call MPI_RECV(coordbuff,ncb,MPI_INTEGER,j,j,dlp_comm_world,status,ierr)
          write(*,*) coordbuff(1:10)
          jj=1
          do ii=1,2*adfv%myarray_len
            nmax=coordbuff(jj)
            if (nmax>adfv%cstat(-1,ii)) then
              adfv%cstat(-1,ii)=nmax
            end if
            if (adfv%cstat(-3,ii)/=coordbuff(jj+1) .or. adfv%cstat(-2,ii)/=coordbuff(jj+2)) then
              write(*,*), "ERROR: coord pairs do not match in MPI"
            end if
            do kk=0,nmax
              adfv%cstat(kk,ii)=adfv%cstat(kk,ii)+coordbuff(jj+3+kk)
            end do
            jj=jj+nmax+4
          end do
        end if
      enddo

      write(nicrdt,*),"Coordination distribution statistics"
      Do i=1,2*adfv%myarray_len
        Do j=0,adfv%cstat(-1,i)
          write(nicrdt,*),unqatm(adfv%cstat(-3,i)),unqatm(adfv%cstat(-2,i)),j,adfv%cstat(j,i)
        End Do
      End Do
    else
    k=1
      Do i=1,2*adfv%myarray_len
        coordbuff(k)=adfv%cstat(-1,i)
        k=k+1
        coordbuff(k:k+1)=adfv%cstat(-3:-2,i)
        k=k+2
        Do j=0,adfv%cstat(-1,i)
          coordbuff(k)=adfv%cstat(j,i)
          k=k+1
        End Do
      End Do
      Call MPI_SEND(ncb,1,MPI_INTEGER,0,idnode,dlp_comm_world,ierr)
      if (ncb>0) then
        Call MPI_SEND(coordbuff,ncb,MPI_INTEGER,0,idnode,dlp_comm_world,ierr)
      endif
    endif
    deallocate(coordbuff)


  end subroutine init_coord_list


  subroutine check_coord(x,y,z,n,mxlist,mxatdm,list,ltg,coordlist,newcoordlist,defectlist,adfv,nstep,time,adfw)
    use config_module, only: atmnam,cfgname,cell
    Use setup_module,  only: ncrdcdt
    Real(kind=wp) , Intent(In) :: x(:)
    Real(kind=wp) , Intent(In) :: y(:)
    Real(kind=wp) , Intent(In) :: z(:)
    Integer , Intent(Out) :: coordlist(0:mxlist,1:mxatdm)
    Integer , Intent(Out):: newcoordlist(0:mxlist,1:mxatdm),defectlist(0:mxlist,0:mxatdm)
    Integer, Intent(In) :: list(-3:mxlist,1:mxatdm)
    Integer , Intent(In) :: n,mxlist,mxatdm
    Integer, Intent(In) :: ltg(:)
    Type(adf_type), Intent(In) :: adfv
    type(adf_type2), Intent(In) :: adfw
    Integer,           Intent( In    ) :: nstep

    Real( Kind = wp ), Intent( In    ) :: time
    integer :: i,ii,j,jj,k,kk, oldnum, newnum,defectcnt
    real :: rcut,rab
    logical :: coordchange,coordfound
    If(adfw%coordon .Eqv. .False.)Return
    Open(Unit=ncrdcdt, File='CCOORD', Form='formatted')
    If(nstep.eq.1)then
      Write(Unit=ncrdcdt, Fmt='(a72)') cfgname(1:72)
      Write(Unit=ncrdcdt, Fmt='(a72)')'File showing change in coordination between atoms'
    endif
    If(mod(nstep,adfw%coordinterval).NE.0)Return
    write(*,*)"Runs check coord",nstep,adfw%coordinterval,adfv%myarray_len

    newcoordlist(0,:)=0
    do j = 1, n
      k=list(0,j)
      Do i=1,k
        kk=list(i,j)
        rab = (x(j)-x(kk))**2+(y(j)-y(kk))**2+(z(j)-z(kk))**2

        rcut=0.00

        Do ii= 1 , adfv%myarray_len
          if ((atmnam(ltg(j))==adfv%arraypairs(ii,1)) .and. (atmnam(ltg(kk))==adfv%arraypairs(ii,2))) Then
            rcut=adfv%arraycuts(ii)*adfv%arraycuts(ii)
          else if ((atmnam(ltg(j))==adfv%arraypairs(ii,2)) .and. (atmnam(ltg(kk))==adfv%arraypairs(ii,1))) Then
            rcut=adfv%arraycuts(ii)*adfv%arraycuts(ii)
          endif
        end Do

        if (rab <= rcut) Then
          newcoordlist(0,j)=newcoordlist(0,j)+1
          newcoordlist(newcoordlist(0,j),j)=ltg(kk)
          if (kk<=n) Then
            newcoordlist(0,kk)=newcoordlist(0,kk)+1
            newcoordlist(newcoordlist(0,kk),kk)=ltg(j)
          endif

        End if
      End Do
    enddo
    write(*,*) coordlist(0:coordlist(0,1),1)
    write(*,*) newcoordlist(0:newcoordlist(0,1),1)

    defectcnt=0
    defectlist(:,:)=0

    do i = 1, n
      coordchange=.False.
      coordfound=.False.
      newnum = newcoordlist(0,i)
      do j = 1, newnum
        coordfound=.False.
        oldnum = coordlist(0,i)
        do k=1, oldnum
          if (newcoordlist(j,i) .eq. coordlist(k,i)) Then
            coordfound=.True.
          endif
        enddo
        if (coordfound .eqv. .False.) Then
          coordchange = .True.
        endif

      enddo

      if (coordchange .eqv. .True.) Then
        defectcnt=defectcnt+1
        defectlist(0,0)=defectcnt
        defectlist(1,defectcnt)=i
      endif
    enddo
    Write(ncrdcdt,Fmt='(A30,I10,I10,f20.6)')'Number of coordination changes',defectcnt,nstep,time
    Do i = 0, 2
      Write(ncrdcdt, Fmt='(3f20.10)') &
        cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 )
    enddo
    Do i=1, defectcnt
      write(ncrdcdt,Fmt='(a6,I10,3f20.10)') &
        atmnam(defectlist(1,i)),defectlist(1,i),x(defectlist(1,i)),y(defectlist(1,i)),z(defectlist(1,i))
    enddo

    coordlist = newcoordlist

!close(unit=ncrdcdt)
  end subroutine check_coord

  subroutine adf_hello
    write (*,*) "Olly says hi..."

  end subroutine adf_hello
end Module adf
