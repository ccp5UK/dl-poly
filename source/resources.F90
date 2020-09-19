!! Author: Alin M Elena
!! Date: 04-07-2020
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0

Module resources
  Use, Intrinsic :: iso_c_binding,   Only: C_INT64_T,&
                                           C_LONG,&
                                           c_int
  Use, Intrinsic :: iso_fortran_env, Only: eu => error_unit,&
                                           wp => real64

  Implicit None

  !  struct timeval {
  !   time_t      tv_sec;   // Number of whole seconds of elapsed time
  !   long int    tv_usec;  // Number of microseconds of rest of elapsed time minus tv_sec. Always less than one million
  !};
  Private
  Type, bind(c) :: timeval
    Integer(C_INT64_T) ::  tv_sec
    Integer(c_long)    :: tv_usec
  End Type

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

  Type, bind(c), Public :: rusage_t
    Type(timeval) :: ru_utime
    Type(timeval) :: ru_stime
    Integer(c_long) :: ru_maxrss
    Integer(c_long) :: ru_ixrss
    Integer(c_long) :: ru_idrss
    Integer(c_long) :: ru_isrss
    Integer(c_long) :: ru_minflt
    Integer(c_long) :: ru_majflt
    Integer(c_long) :: ru_nswap
    Integer(c_long) :: ru_inblock
    Integer(c_long) :: ru_oublock
    Integer(c_long) :: ru_msgsnd
    Integer(c_long) :: ru_msgrcv
    Integer(c_long) :: ru_nsignals
    Integer(c_long) :: ru_nvcsw
    Integer(c_long) :: ru_nivcsw
  End Type

  !int getrusage(int who, struct rusage *usage);
  Integer, Parameter, Public :: ru_maxrss = 1
  Integer, Parameter, Public :: ru_ixrss = 2
  Integer, Parameter, Public :: ru_idrss = 3
  Integer, Parameter, Public :: ru_isrss = 4
  Integer, Parameter, Public :: ru_minflt = 5
  Integer, Parameter, Public :: ru_majflt = 6
  Integer, Parameter, Public :: ru_nswap = 7
  Integer, Parameter, Public :: ru_inblock = 8
  Integer, Parameter, Public :: ru_oublock = 9
  Integer, Parameter, Public :: ru_msgsnd = 10
  Integer, Parameter, Public :: ru_msgrcv = 11
  Integer, Parameter, Public :: ru_nsignals = 12
  Integer, Parameter, Public :: ru_nvcsw = 13
  Integer, Parameter, Public :: ru_nivcsw = 14

  Interface

    Integer(c_int) Function getrusage(who, rusage) bind(C, name="getrusage")
      import c_int, rusage_t
      Integer(c_int), Intent(in), value :: who
      Type(rusage_t) :: rusage
    End Function getrusage

  End Interface
  Public :: getrusage
  Public :: print_res
  Public :: mem
Contains

  Subroutine print_res(sa)
    Type(rusage_t), Intent(In   ) :: sa

    Write (eu, '(a,1x,i0,a)') "ru_maxrss vmrss", sa%ru_maxrss, "KiB"
    Write (eu, '(a,1x,i0)') "ru_minflt minor fault page", sa%ru_minflt
    Write (eu, '(a,1x,i0)') "ru_majflt  major fault page", sa%ru_majflt
    Write (eu, '(a,1x,i0)') "ru_inblock io i", sa%ru_inblock
    Write (eu, '(a,1x,i0)') "ru_oublock io o", sa%ru_oublock
    Write (eu, '(a,1x,i0)') "ru_nvcsw voluntary", sa%ru_nvcsw
    Write (eu, '(a,1x,i0)') "ru_nivcsw invonlutary", sa%ru_nivcsw

  End Subroutine print_res

  Character(len=16) Function mem(m)
    Real(kind=wp), Intent(In   ) :: m

    Integer         :: i
    Integer(kind=8) :: mm(4)

    Do i = 1, 4
      mm(i) = Ishft(1_8, i * 10)
    End Do
    mem = ''
    If (m < mm(1)) Then
      Write (mem, '(f12.4,a4)') m, " B"
    Elseif (m < mm(2)) Then
      Write (mem, '(f12.4,a4)') m / mm(1), " KiB"
    Elseif (m < mm(3)) Then
      Write (mem, '(f12.4,a4)') m / mm(2), " MiB"
    Elseif (m < mm(4)) Then
      Write (mem, '(f12.4,a4)') m / mm(3), " GiB"
    Else
      Write (mem, '(f16.4,a4)') m / mm(4), " TiB"
    Endif
  End Function mem

End Module resources
