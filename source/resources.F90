!! Author: Alin M Elena
!! Date: 04-07-2020
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0

module resources
  use, intrinsic :: iso_c_binding,   only: C_INT64_T,&
                                           C_LONG,&
                                           c_int
  use, intrinsic :: iso_fortran_env, only: eu => error_unit,&
                                           wp => real64

  implicit none

  !  struct timeval {
  !   time_t      tv_sec;   // Number of whole seconds of elapsed time
  !   long int    tv_usec;  // Number of microseconds of rest of elapsed time minus tv_sec. Always less than one million
  !};
  private
  type, bind(c) :: timeval
    integer(C_INT64_T) ::  tv_sec
    integer(c_long)    :: tv_usec
  end type

  !struct rusage {
  !      struct timeval ru_utime; /* user CPU time used */
  !      struct timeval ru_stime; /* system CPU time used */
  !      long   ru_maxrss;        /* maximum resident set size */
  !      long   ru_ixrss;         /* integral shared memory size */
  !      long   ru_idrss;         /* integral unshared data size */
  !      long   ru_isrss;         /* integral unshared stack size */
  !      long   ru_minflt;        /* page reclaims (soft page faults) */
  !      long   ru_majflt;        /* page faults (hard page faults) */
  !      long   ru_nswap;         /* swaps */
  !      long   ru_inblock;       /* block input operations */
  !      long   ru_oublock;       /* block output operations */
  !      long   ru_msgsnd;        /* IPC messages sent */
  !      long   ru_msgrcv;        /* IPC messages received */
  !      long   ru_nsignals;      /* signals received */
  !      long   ru_nvcsw;         /* voluntary context switches */
  !      long   ru_nivcsw;        /* involuntary context switches */
  !  };

  type, bind(c), public :: rusage_t
    type(timeval) :: ru_utime
    type(timeval) :: ru_stime
    integer(c_long) :: ru_maxrss
    integer(c_long) :: ru_ixrss
    integer(c_long) :: ru_idrss
    integer(c_long) :: ru_isrss
    integer(c_long) :: ru_minflt
    integer(c_long) :: ru_majflt
    integer(c_long) :: ru_nswap
    integer(c_long) :: ru_inblock
    integer(c_long) :: ru_oublock
    integer(c_long) :: ru_msgsnd
    integer(c_long) :: ru_msgrcv
    integer(c_long) :: ru_nsignals
    integer(c_long) :: ru_nvcsw
    integer(c_long) :: ru_nivcsw
  end type

  !int getrusage(int who, struct rusage *usage);
  integer, parameter, public :: ru_maxrss = 1
  integer, parameter, public :: ru_ixrss = 2
  integer, parameter, public :: ru_idrss = 3
  integer, parameter, public :: ru_isrss = 4
  integer, parameter, public :: ru_minflt = 5
  integer, parameter, public :: ru_majflt = 6
  integer, parameter, public :: ru_nswap = 7
  integer, parameter, public :: ru_inblock = 8
  integer, parameter, public :: ru_oublock = 9
  integer, parameter, public :: ru_msgsnd = 10
  integer, parameter, public :: ru_msgrcv = 11
  integer, parameter, public :: ru_nsignals = 12
  integer, parameter, public :: ru_nvcsw = 13
  integer, parameter, public :: ru_nivcsw = 14

  interface

    integer(c_int) function getrusage(who, rusage) bind(C, name="getrusage")
      import c_int, rusage_t
      integer(c_int), intent(in), value :: who
      type(rusage_t) :: rusage
    end function getrusage

  end interface
  public :: getrusage
  public :: print_res
  public :: mem
contains

  subroutine print_res(sa)
    type(rusage_t), intent(in)    :: sa

    write (eu, '(a,1x,i0,a)') "ru_maxrss vmrss", sa%ru_maxrss, "KiB"
    write (eu, '(a,1x,i0)') "ru_minflt minor fault page", sa%ru_minflt
    write (eu, '(a,1x,i0)') "ru_majflt  major fault page", sa%ru_majflt
    write (eu, '(a,1x,i0)') "ru_inblock io i", sa%ru_inblock
    write (eu, '(a,1x,i0)') "ru_oublock io o", sa%ru_oublock
    write (eu, '(a,1x,i0)') "ru_nvcsw voluntary", sa%ru_nvcsw
    write (eu, '(a,1x,i0)') "ru_nivcsw invonlutary", sa%ru_nivcsw

  end subroutine print_res

  character(len=16) function mem(m)
    real(kind=wp), intent(in)    :: m

    integer         :: i
    integer(kind=8) :: mm(4)

    do i = 1, 4
      mm(i) = ishft(1_8, i * 10)
    end do
    mem = ''
    if (m < mm(1)) then
      write (mem, '(f12.4,a4)') m, " B"
    elseif (m < mm(2)) then
      write (mem, '(f12.4,a4)') m / mm(1), " KiB"
    elseif (m < mm(3)) then
      write (mem, '(f12.4,a4)') m / mm(2), " MiB"
    elseif (m < mm(4)) then
      write (mem, '(f12.4,a4)') m / mm(3), " GiB"
    else
      write (mem, '(f16.4,a4)') m / mm(4), " TiB"
    endif
  end function mem

end module resources
