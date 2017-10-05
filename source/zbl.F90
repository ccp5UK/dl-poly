Module m_zbl
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring constants and ZBL related routines
  !
  ! copyright - daresbury laboratory
  ! author    - a.m.elena october 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use Kinds_f90, only : wp
  Implicit None
  Private 

  Real( Kind = wp ),parameter,dimension(4) :: b=[0.18175_wp,0.50986_wp,0.28022_wp,0.02817_wp], &
    c=[3.1998_wp,0.94229_wp,0.40290_wp,0.20162_wp]
  Real( Kind = wp ),parameter,public :: ab=0.52917721067_wp 
  Public :: zbl

Contains

  Pure Subroutine zbl(r,k,ia,phi,dphi)
    Real(wp), intent(in)  :: r,k,ia
    Real(wp), intent(out) :: phi,dphi

    Integer :: i
    Real(wp) :: x,t1,t2,ir

    phi = 0.0_wp
    dphi = 0.0_wp
    x = r*ia
    ir = 1.0_wp/r
    do i = 1,4
      t1 = b(i)*exp(-x*c(i))
      phi = phi+t1
      dphi = dphi-c(i)*t1
    end do
    phi = k*phi*ir
    ! -r∂U/δr
    dphi = phi - ia*k*dphi

  End Subroutine 

end module m_zbl
