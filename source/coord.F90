Module coord

!Module to calculate the coordination number distributions during a DLPOLY run
!
! author  Oliver Dicks and Aaron Diver
! date    22/02/2018
  Use site, Only : site_type
  Use configuration, Only : configuration_type
  Use neighbours, Only : neighbours_type 
  Use flow_control, Only : flow_type  
  Use kinds, only : wp
  Use comms, Only : comms_type,grecv,gsend
  Use constants, Only : nicrdt,nccrdt


  Implicit None

  type, public :: coord_type
    Integer :: ncoordpairs,coordinterval,coordstart,coordops
    real(wp), allocatable :: arraycuts(:)
    character( Len = 8 ), allocatable :: arraypairs(:,:)!,coordatoms(:)
    Integer, allocatable :: coordlist(:,:),icoordlist(:,:),defectlist(:)
    integer, allocatable :: ltype(:,:),cstat(:,:)
    Logical :: coordon 
  contains
    procedure :: init=>init_coord
    procedure :: init_coordlist
    final :: clean_coord
  end type coord_type
  private
  public :: init_coord_list,checkcoord

contains

  subroutine init_coord(T)
    class(coord_type) :: T
    allocate(T%arraycuts(1:T%ncoordpairs))
    allocate(T%arraypairs(1:T%ncoordpairs,1:2))
    allocate(T%ltype(1:T%ncoordpairs,1:2))
    allocate(T%cstat(-3:1000,1:(2*T%ncoordpairs)))
!   allocate(T%coordatoms(0:2*T%ncoordpairs))
  end subroutine init_coord
  subroutine init_coordlist(T,n,m)
    class(coord_type) :: T
    Integer, Intent(In) :: n,m
    allocate(T%coordlist(0:n,1:m))
    allocate(T%icoordlist(0:n,1:m))
    allocate(T%defectlist(0:n))
  end subroutine init_coordlist

  subroutine clean_coord(T)
    type(coord_type) :: T

    if (allocated(T%arraycuts)) then
      deallocate(T%arraycuts)
    end if
    if (allocated(T%coordlist)) then
      deallocate(T%coordlist)
    end if
    if (allocated(T%icoordlist)) then
      deallocate(T%icoordlist)
    end if
    if (allocated(T%defectlist)) then
      deallocate(T%defectlist)
    end if
  end subroutine clean_coord

  subroutine init_coord_list(config,neigh,crd,sites,flow,comm)
    Type(neighbours_type), Intent(In) :: neigh
    Type(configuration_type), Intent(In)  :: config
    Type(comms_type), Intent(InOut) :: comm
    Type(flow_type), Intent(In) :: flow
    Type(site_type), Intent(In) :: sites
    Type(coord_type), Intent(InOut) :: crd 
    integer :: i,ii,j,jj,k,kk, ncoord,en,ierr,m,mcoord,lgcoord,nmax,ncb !ncb size coordbuff
    real :: rcut,rab
    Character(len=8) :: aux 
    Character(len=1000),allocatable :: cbuff(:)
    integer, allocatable :: buff(:),coordbuff(:)
    Logical :: newatom,itsopen

    !Check whether option is called

    If(crd%coordon .Eqv. .False.)Return
    If(crd%ncoordpairs==0)Return
    If(crd%coordstart>flow%step)Return
    If(mod(flow%step,crd%coordinterval).NE.0)Return
    crd%coordlist(0,:)=0
    ncb=0
    Do j = 1, config%natms
      ncoord = 0
      k=neigh%list(0,j)
      Do i=1,k
        kk=neigh%list(i,j)
        rab = (config%parts(j)%xxx-config%parts(kk)%xxx)**2+(config%parts(j)%yyy-config%parts(kk)%yyy)**2 &
              + (config%parts(j)%zzz-config%parts(kk)%zzz)**2

        rcut=0.00

        Do ii= 1 , crd%ncoordpairs
          if (((config%ltype(j)==crd%ltype(ii,1)) .and. (config%ltype(kk)==crd%ltype(ii,2)))&
            .or.&
            ((config%ltype(j)==crd%ltype(ii,2)) .and. (config%ltype(kk)==crd%ltype(ii,1)))) Then
            rcut=crd%arraycuts(ii)*crd%arraycuts(ii)
          endif
        end Do

        if (rab <= rcut) Then
          crd%coordlist(0,j)=crd%coordlist(0,j)+1
          crd%coordlist(crd%coordlist(0,j),j)=config%ltg(kk)
          if (kk<=config%natms) Then
            crd%coordlist(0,kk)=crd%coordlist(0,kk)+1
            crd%coordlist(crd%coordlist(0,kk),kk)=config%ltg(j)
          endif
        End if
      End Do
    End Do

     !Create icoordlist
     if (flow%step==crd%coordstart) then
      crd%icoordlist=crd%coordlist
     end if

    !Create coordination statistics array
    crd%cstat(-3:,:)=0
    Do i=1,crd%ncoordpairs
      crd%cstat(-3,(2*i)-1)=crd%ltype(i,1)
      crd%cstat(-2,(2*i)-1)=crd%ltype(i,2)
      crd%cstat(-3,(2*i))=crd%ltype(i,2)
      crd%cstat(-2,(2*i))=crd%ltype(i,1)
    End Do

    !Collect coordination statistics
    Do i=1,config%natms
      Do j=1,2*crd%ncoordpairs
      mcoord=0
        if ((config%ltype(i)) == crd%cstat(-3,j)) then
          Do k=1,crd%coordlist(0,i)
            if (config%ltype(crd%coordlist(k,i)) == crd%cstat(-2,j)) then
              mcoord=mcoord+1
            end if
          End do
        crd%cstat(mcoord,j)=crd%cstat(mcoord,j)+1
        if (mcoord>crd%cstat(-1,j)) then
          crd%cstat(-1,j)=mcoord
        End if
        End if
      End do 
    End do

   
!Set coordbuff size
    Do i=1,2*crd%ncoordpairs
      ncb=ncb+(crd%cstat(-1,i)+4)
    End do


    allocate(buff(2*config%mxatms))
    allocate(cbuff(config%mxatms))
    If (comm%idnode==0) Then
      inquire(Unit=nicrdt, opened=itsopen)
      if ( .not. itsopen ) Then
        Open(Unit=nicrdt, File='ICOORD', Form='formatted')
      end if
      Write(Unit=nicrdt, Fmt='(a72)') config%cfgname(1:72)
      Write(Unit=nicrdt, Fmt='(a60,I10)')'Initial coordination between atoms',flow%step
      If(crd%coordops .eq.0) then
      Do i=1,config%natms
        m=crd%coordlist(0,i)
        write (nicrdt,Fmt='(i12,1x,i12,1x)',advance="no") &
          config%ltg(i),crd%coordlist(0,i)
        do ii=1,m
          write(nicrdt,'(i0,1x)',advance="no") config%ltg(crd%coordlist(ii,i))
        enddo 
        do ii=1,m
          write(nicrdt,'(a,1x)',advance="no") trim(sites%unique_atom(config%ltype(crd%coordlist(ii,i))))
        enddo 
        write(nicrdt,*)
      enddo


      do j=1,comm%mxnode-1
        Call grecv(comm,en,j,j)
        if (en>0) Then
          Call grecv(comm,buff,j,j)
          Call grecv(comm,cbuff,j,j)
          Do i=1,en/2
            write (nicrdt,Fmt='(i12,1x,i12,a)') &
              buff(2*i-1),buff(2*i),trim(cbuff(i))
          enddo
        endif
      enddo
      endif
    else
     If(crd%coordops .eq.0) then
      en=2*config%natms
      do i=1,config%natms
        buff(2*i-1) = config%ltg(i)
        buff(2*i) = crd%coordlist(0,i)
        cbuff(i)=''
        do ii=1,crd%coordlist(0,i)
          write(aux,'(i0)') config%ltg(crd%coordlist(ii,i))
          cbuff(i)=trim(cbuff(i))//" "//trim(aux)
        enddo 
        do ii=1,crd%coordlist(0,i)
          write(aux,'(a)') trim(sites%unique_atom(config%ltype(crd%coordlist(ii,i))))
          cbuff(i)=trim(cbuff(i))//" "//trim(aux)
        enddo 
      enddo
      Call gsend(comm,en,0,comm%idnode)
      if (en>0) then
        Call gsend(comm,buff,0,comm%idnode)
        Call gsend(comm,cbuff,0,comm%idnode)
      endif
    endif
    endif
    deallocate(buff)
    deallocate(cbuff)

    
    allocate(coordbuff(ncb))

    If (comm%idnode==0) Then
      Open(Unit=nicrdt, File='ICOORD', Form='formatted')

      do j=1,comm%mxnode-1
        Call grecv(comm,ncb,j,j)
        if (ncb>1) then
          Call grecv(comm,coordbuff,j,j)
          jj=1
          do ii=1,2*crd%ncoordpairs
            nmax=coordbuff(jj)
            if (nmax>crd%cstat(-1,ii)) then
              crd%cstat(-1,ii)=nmax
            end if
            if (crd%cstat(-3,ii)/=coordbuff(jj+1) .or. crd%cstat(-2,ii)/=coordbuff(jj+2)) then
              write(*,*) "ERROR: coord pairs do not match in MPI"
            end if
            do kk=0,nmax
              crd%cstat(kk,ii)=crd%cstat(kk,ii)+coordbuff(jj+3+kk)
            end do
            jj=jj+nmax+4
          end do
        end if
      enddo

      write(nicrdt,*)"Coordination distribution statistics"
      Do i=1,2*crd%ncoordpairs
        Do j=0,crd%cstat(-1,i)
          write(nicrdt,*)sites%unique_atom(crd%cstat(-3,i)),sites%unique_atom(crd%cstat(-2,i)),j,crd%cstat(j,i)
        End Do
      End Do
    else
    k=1
      Do i=1,2*crd%ncoordpairs
        coordbuff(k)=crd%cstat(-1,i)
        k=k+1
        coordbuff(k:k+1)=crd%cstat(-3:-2,i)
        k=k+2
        Do j=0,crd%cstat(-1,i)
          coordbuff(k)=crd%cstat(j,i)
          k=k+1
        End Do
      End Do
      Call gsend(comm,ncb,0,comm%idnode)
      if (ncb>0) then
        Call gsend(comm,coordbuff,0,comm%idnode)
      endif
    endif
    deallocate(coordbuff)

  end subroutine init_coord_list


  subroutine checkcoord(config,neigh,crd,sites,flow,comm) 
    Type(neighbours_type), Intent(In) :: neigh
    Type(configuration_type), Intent(In)  :: config
    Type(comms_type), Intent(InOut) :: comm
    Type(flow_type), Intent(In) :: flow
    Type(site_type), Intent(In) :: sites
    Type(coord_type), Intent(InOut) :: crd

    integer :: i,ii,j,jj,k,kk,defn,defectcnt,totdefectcnt
    real :: rcut,rab
    integer, allocatable :: buff(:)
    logical :: coordchange,coordfound,thisopen
    
    If(crd%coordon .Eqv. .False.)Return
    If(crd%ncoordpairs==0)Return
    If(crd%coordops .eq.1)Return
    If(crd%coordstart>flow%step)Return    
    If(mod(flow%step,crd%coordinterval).NE.0)Return

    defectcnt=0
    crd%defectlist(:)=0

   do i = 1, config%natms 
     coordchange=.False.
     coordfound=.False.
     if (crd%coordlist(0,i) /= crd%icoordlist(0,i)) then
      coordchange=.true.
     else
       do j = 1, crd%coordlist(0,i)
        coordfound=.False.
         do k=1, crd%icoordlist(0,i)
           if (crd%coordlist(j,i) == crd%icoordlist(k,i)) Then
             coordfound=.True.
            endif
          enddo
          if (coordfound .eqv. .False.) Then
            coordchange = .True.
          endif
        enddo
      end if
      if (coordchange .eqv. .True.) Then
        defectcnt=defectcnt+1
        crd%defectlist(0)=defectcnt
        crd%defectlist(defectcnt)=i
      endif
    enddo
  allocate(buff(2*crd%defectlist(0))) 
    If (comm%idnode==0) Then
      totdefectcnt=0
      do j=1,comm%mxnode-1
        Call grecv(comm,defn,j,j)
        totdefectcnt=totdefectcnt+defn
      end do
      totdefectcnt=totdefectcnt+crd%defectlist(0)

      inquire(Unit=nccrdt, opened=thisopen)
      if ( .not. thisopen ) Then
        Open(Unit=nccrdt, File='CCOORD', Form='formatted')
      end if
      Write(Unit=nccrdt,Fmt='(A30,I10,I10,f20.6)')'Number of coordination changes',totdefectcnt
!    Do i = 0, 2
!      Write(ncrdcdt, Fmt='(3f20.10)') &
!        cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 )
!    enddo
      Do i=1, crd%defectlist(0)
        write(Unit=nccrdt,Fmt='(a6,I10,3f20.10)') &
          trim(sites%unique_atom(config%ltype(crd%defectlist(i)))),config%ltg(crd%defectlist(i)), &
          config%parts(crd%defectlist(i))%xxx,config%parts(crd%defectlist(i))%yyy,config%parts(crd%defectlist(i))%zzz
      enddo

      do j=1,comm%mxnode-1
        Call grecv(comm,defectcnt,j,j)
        if (defectcnt>0) Then
          Call grecv(comm,buff,j,j)
          do i=1,defectcnt
            write(nccrdt,Fmt='(a6,I10)') trim(sites%unique_atom(buff(2*i-1))),buff(2*i)
          end do
        end if
      end do

    else
      defectcnt=crd%defectlist(0)
      defn=crd%defectlist(0)
      Call gsend(comm,defn,0,comm%idnode)
      do i =1, defectcnt
        buff(2*i-1)=config%ltype(crd%defectlist(i))
        buff(2*i)=config%ltg(crd%defectlist(i))
      End do
      Call gsend(comm,defectcnt,0,comm%idnode)
      if (defectcnt>0) then
        Call gsend(comm,buff,0,comm%idnode)
      End if
    End if
    deallocate(buff)

  end subroutine checkcoord

end Module coord
