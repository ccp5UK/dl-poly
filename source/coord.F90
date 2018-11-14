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



  Implicit None

  type, public :: coord_type
    Integer :: ncoordpairs,coordinterval,coordstart
    real(wp), allocatable :: arraycuts(:)
    character( Len = 8 ), allocatable :: arraypairs(:,:)!,coordatoms(:)
    Integer, allocatable :: coordlist(:,:),newcoordlist(:,:),defectlist(:,:)
    integer,allocatable :: ltype(:,:),cstat(:,:)
    Logical :: coordon 
  contains
    procedure :: init=>init_coord
    procedure :: init_coordlist
    final :: clean_coord
  end type coord_type
  private
  public :: coord_hello,init_coord_list

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
  end subroutine init_coordlist

  subroutine clean_coord(T)
    type(coord_type) :: T

    if (allocated(T%arraycuts)) then
      deallocate(T%arraycuts)
    end if
    if (allocated(T%coordlist)) then
      deallocate(T%coordlist)
    end if
  end subroutine clean_coord

  subroutine init_coord_list(config,neigh,crd,sites,flow,comm)
    Type(neighbours_type), Intent(In) :: neigh
    Type(configuration_type), Intent(In)  :: config
    Type(comms_type), Intent(InOut) :: comm
    Type(flow_type), Intent(In) :: flow
    Type(site_type), Intent(In) :: sites
    Type(coord_type), Intent(InOut) :: crd 
    integer :: i,ii,j,jj,k,kk, ncoord,en,ierr,m,mcoord,lgcoord,nmax,ncb,nicrdt !ncb size coordbuff
    real :: rcut,rab
    Character(len=8) :: aux 
    Character(len=1000),allocatable :: cbuff(:)
    integer, allocatable :: buff(:),coordbuff(:)
    Logical :: newatom
    If(crd%coordon .Eqv. .False.)Return
    If(crd%ncoordpairs==0)Return
    If(mod(flow%step,crd%coordinterval).NE.0)Return
    crd%coordlist(0,:)=0
    ncb=0
    nicrdt=2323
    do j = 1, config%natms
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
          crd%coordlist(crd%coordlist(0,j),j)=kk
          if (kk<=config%natms) Then
            crd%coordlist(0,kk)=crd%coordlist(0,kk)+1
            crd%coordlist(crd%coordlist(0,kk),kk)=j
          endif
        End if
      End Do
    End Do
    crd%cstat(-3:,:)=0
    Do i=1,crd%ncoordpairs
      crd%cstat(-3,(2*i)-1)=crd%ltype(i,1)
      crd%cstat(-2,(2*i)-1)=crd%ltype(i,2)
      crd%cstat(-3,(2*i))=crd%ltype(i,2)
      crd%cstat(-2,(2*i))=crd%ltype(i,1)
    End Do

    
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
      Open(Unit=nicrdt, File='ICOORD', Form='formatted')
      Write(Unit=nicrdt, Fmt='(a72)') config%cfgname(1:72)
      Write(Unit=nicrdt, Fmt='(a60,I10)')'Initial coordination between atoms',flow%step
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
!      Close(Unit=nicrdt)
    else
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

  subroutine checkcoord(config,neigh,crd,sites,comm)
!    Integer , Intent(Out):: newcrd%coordlist(0:neigh%max_list,1:config%mxatms),defectlist(0:neigh%max_list,0:config%mxatms)
  
    Type(neighbours_type), Intent(In) :: neigh
    Type(configuration_type), Intent(In)  :: config
    Type(comms_type), Intent(InOut) :: comm
    Type(flow_type), Intent(In) :: flow
    Type(site_type), Intent(In) :: sites
    Type(coord_type), Intent(InOut) :: crd
!
    integer :: i,ii,j,jj,k,kk, oldnum, newnum,defectcnt
    real :: rcut,rab
    logical :: coordchange,coordfound
    If(crd%coordon .Eqv. .False.)Return
    Open(Unit=ncrdcdt, File='CCOORD', Form='formatted')
    If(flow%step.eq.1)then
      Write(Unit=ncrdcdt, Fmt='(a72)') cfgname(1:72)
      Write(Unit=ncrdcdt, Fmt='(a72)')'File showing change in coordination between atoms'
    endif
    If(mod(flow%step,crd%coordinterval).NE.0)Return
    write(*,*)"Runs check coord",flow%step,crd%coordinterval,crd%ncoordpairs
!
    newcrd%coordlist(0,:)=0
    do j = 1, config%natms
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
          coordlist(0,j)=coordlist(0,j)+1
          coordlist(coordlist(0,j),j)=kk
          if (kk<=n) Then
            coordlist(0,kk)=coordlist(0,kk)+1
            coordlist(coordlist(0,kk),kk)=j
          endif
        End if
      End Do
    End Do
    
    write(*,*) crd%coordlist(0:crd%coordlist(0,1),1)
    write(*,*) crd%newcoordlist(0:crd%newcoordlist(0,1),1)

    defectcnt=0
    defectlist(:,:)=0

    do i = 1, config%natms 
      coordchange=.False.
      coordfound=.False.
      newnum = crd%newcoordlist(0,i)
      do j = 1, newnum
       coordfound=.False.
       oldnum = crd%coordlist(0,i)
        do k=1, oldnum
          if (crd%newcoordlist(j,i) .eq. crd%coordlist(k,i)) Then
            coordfound=.True.
          endif
        enddo
        if (coordfound .eqv. .False.) Then
          coordchange = .True.
        endif
!
      enddo
!
      if (coordchange .eqv. .True.) Then
        defectcnt=defectcnt+1
        defectlist(0,0)=defectcnt
        defectlist(1,defectcnt)=i
      endif
    enddo
    Write(ncrdcdt,Fmt='(A30,I10,I10,f20.6)')'Number of coordination changes',defectcnt,nstep,time
!    Do i = 0, 2
!      Write(ncrdcdt, Fmt='(3f20.10)') &
!        cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 )
!    enddo
!    Do i=1, defectcnt
!      write(ncrdcdt,Fmt='(a6,I10,3f20.10)') &
!        atmnam(defectlist(1,i)),defectlist(1,i),x(defectlist(1,i)),y(defectlist(1,i)),z(defectlist(1,i))
!    enddo
!
!    crd%coordlist = newcrd%coordlist
!
!!close(unit=ncrdcdt)
!  end subroutine check_coord

  subroutine coord_hello
    write (*,*) "Olly says hi..."

  end subroutine coord_hello
end Module coord
