module currents
  use kinds, only : wp
  Use configuration, Only : configuration_type
  Use constants, Only : czero
  Use comms, Only : comms_type, gsum
  Use filename, Only : file_type
  implicit none

  type, public :: current_type
    complex(kind=wp),allocatable :: jlk(:,:)
    complex(kind=wp), allocatable :: c(:,:,:) ! lags,kpoints,xyz
    complex(kind=wp), allocatable :: fc(:,:,:) ! lags,kpoints,xyz
    integer :: nkpoints,lag
    integer :: file_handle
    logical :: on = .False.
  contains
    private
    procedure,public :: init
    procedure,public :: compute
    final :: cleanup
  end type

contains

  subroutine init(T,nk,lag,fcurrent,comm)
  class(current_type) :: T
    Type( file_type ), Intent( InOut ) :: fcurrent
    Type(comms_type), intent( In ) :: comm

    integer, intent(in) :: nk
    integer, intent(in) :: lag

    allocate(T%jlk(nk,3))
    T%nkpoints = nk
    T%lag=lag
    T%on = .True.
    If (comm%idnode == 0) Open(Newunit=fcurrent%unit_no, File=fcurrent%filename, Status='unknown',&
      action="write")
    T%file_handle = fcurrent%unit_no
  end subroutine init


  subroutine compute(T,config,comm)
  class(current_type) :: T
    type(configuration_type), intent(in) :: config
    Type(comms_type), intent(inout) :: comm

    integer :: i,k 
    real(kind=wp) :: tmp
    complex(kind=wp) :: h(3)

      do k=1,config%k%n
        T%jlk(k,:)=czero
        h = czero
        do i=1,config%natms
          tmp=dot_product(config%k%r(:,k),[config%parts(i)%xxx,config%parts(i)%yyy,config%parts(i)%zzz])
          h = h + [config%vxx(i),config%vyy(i),config%vzz(i)]*exp(cmplx(0.0_wp,tmp,wp))
        enddo
        !current%jlk(:,k,j)=kp%u(:,k)*dot_product(kp%u(:,k),h)
        call gsum(comm,h)
        T%jlk(k,:) = h
      end do

      if (comm%idnode==0) then
        write(T%file_handle,*) T%jlk(:,:)
      end if  
    end subroutine compute

  subroutine cleanup(T)
    type(current_type) :: T
    if (allocated(T%jlk)) deallocate(T%jlk)
    close(t%file_handle)
  end subroutine cleanup

end module currents

