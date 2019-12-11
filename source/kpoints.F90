module kpoints
  use kinds, only : wp
  Use constants, Only : epsilon_wp
  Use comms, Only : comms_type, gbcast
  implicit none
  private

  type, public :: kpoints_type
    integer :: n
    real(kind = wp),allocatable :: r(:,:)
    real(kind = wp),allocatable :: u(:,:)
    integer :: uf
  contains
    private 
    procedure, public :: init

    final :: cleanup
  end type

contains

  subroutine init(T,  filename,comm)
  class(kpoints_type)               :: T
    character(len=*), intent(in)    :: filename
    type(comms_type), intent(inout) :: comm
    
    integer         :: i
    real(kind = wp) :: h


    if (comm%idnode == 0) then 
      open(file=trim(filename),newunit=T%uf, action="read",status="old")
      read(t%uf,*) T%n
    endif

    call gbcast(comm,T%n,0)
    
    if (T%n>0) allocate(T%r(3,t%n))
    if (T%n>0) allocate(T%u(3,t%n))
    
    if (comm%idnode == 0) then 
      do i = 1,t%n 
        read(t%uf,*) t%r(:,i)
      end do
    endif
    call gbcast(comm,T%r,0)
      do i=1,T%n
        h=norm2(T%r(:,i))
        if (abs(h) > epsilon_wp) & 
          T%u(:,i) = T%r(:,i)/h
      end do
    if (comm%idnode == 0) then 
      close(t%uf)
    endif

  end subroutine init

  subroutine cleanup(T)
    type(kpoints_type) :: T
    if (allocated(t%r)) deallocate(t%r)
    if (allocated(t%u)) deallocate(t%u)
  end subroutine cleanup

end module kpoints
