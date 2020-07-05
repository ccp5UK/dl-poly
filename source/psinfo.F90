module psInfo
  use, intrinsic :: iso_c_binding,   only: C_LONG
  use, intrinsic :: iso_fortran_env, only: li => int64,&
                                           wp => real64,&
                                           ou => output_unit

  implicit none
  private
  type, public :: psType
    integer(kind=li)   :: pid
    character(len=16) :: comm
    character(len=1)  :: state
    integer(kind=li)   :: ppid
    integer(kind=li)   :: pgrp
    integer(kind=li)   :: session
    integer(kind=li)   :: tty_nr
    integer(kind=li)   :: tpgid
    integer(kind=li)   :: flags
    real(kind=wp)      :: minflt
    real(kind=wp)      :: cminflt
    real(kind=wp)      :: majflt
    real(kind=wp)      :: cmajflt
    real(kind=wp)      :: utime
    real(kind=wp)      :: stime
    real(kind=wp)      :: cutime
    real(kind=wp)      :: cstime
    integer(kind=li)   :: priority
    integer(kind=li)   :: nice
    integer(kind=li)   :: num_threads
    real(kind=wp)      :: itrealvalue
    real(kind=wp)      :: starttime
    integer(kind=c_long)      :: vsize
    real(kind=wp)      :: rss
    real(kind=wp)      :: rsslim
    real(kind=wp)      :: startcode
    real(kind=wp)      :: endcode
    real(kind=wp)      :: startstack
    real(kind=wp)      :: kstkesp
    real(kind=wp)      :: kstkeip
    integer(kind=li)   :: signal
    real(kind=wp)      :: blocked
    integer(kind=li)   :: sigignore
    integer(kind=li)   :: sigcatch
    real(kind=wp)      :: wchan
    real(kind=wp)      :: nswap
    real(kind=wp)      :: cnswap
    integer(kind=li)   :: exit_signal
    integer(kind=li)   :: processor
    integer(kind=li)   :: rt_priority
    integer(kind=li)   :: policy
    real(kind=wp)      :: delayacct_blkio_ticks
    real(kind=wp)      :: guest_time
    real(kind=wp)      :: cguest_time
  contains
    private
    procedure, public :: inspect => myInfoStat
    procedure, public :: print => printMyInfoStat
  end type psType

  integer, parameter, public   :: sz = 1 ! total # of pages of memory
  integer, parameter, public   :: resident = 2 ! number of resident set (non-swapped) pages (4k)
  integer, parameter, public   :: share = 3 ! number of pages of shared (mmap'd) memory
  integer, parameter, public   :: textrs = 4 ! text resident size
  integer, parameter, public   :: librs = 5 ! shared libs resident size
  integer, parameter, public   :: datars = 6 ! data resident size
  integer, parameter, public   :: dirty = 7 ! dirty pages
  type, public ::pgsType
    integer(kind=li) :: data(7)

  contains
    private
    procedure, public :: inspect => myInfoStatm
    procedure, public :: print => printMyInfoStatm
  end type pgsType

  integer(kind=4), parameter, public   :: VmPeak = 1
  integer(kind=4), parameter, public   :: VmSize = 2
  integer(kind=4), parameter, public   :: VmLck = 3
  integer(kind=4), parameter, public   :: VmPin = 4
  integer(kind=4), parameter, public   :: VmHWM = 5
  integer(kind=4), parameter, public   :: VmRSS = 6
  integer(kind=4), parameter, public   :: RssAnon = 7
  integer(kind=4), parameter, public   :: RssFile = 8
  integer(kind=4), parameter, public   :: RssShmem = 9
  integer(kind=4), parameter, public   :: VmData = 10
  integer(kind=4), parameter, public   :: VmStk = 11
  integer(kind=4), parameter, public   :: VmExe = 12
  integer(kind=4), parameter, public   :: VmLib = 13
  integer(kind=4), parameter, public   :: VmPTE = 14
  integer(kind=4), parameter, public   :: VmSwap = 15
  integer(kind=4), parameter, public   :: HugetlbPages = 16
  integer(kind=4), parameter, public   :: Threads = 17
  integer(kind=4), parameter, public   :: voluntary_ctxt_switches = 18
  integer(kind=4), parameter, public   :: nonvoluntary_ctxt_switches = 19

  type, public :: pssType
    integer :: n = 19
    integer(kind=li) :: data(19)
    character(len=3)  :: sunit
  contains
    private
    procedure, public :: inspect => myInfoStatus
    procedure, public :: print => printMyInfoStatus

  end type pssType

contains

  subroutine myInfoStat(ps)
    class(psType) :: ps

    integer :: u

    open (newunit=u, file="/proc/self/stat", form="formatted", status="old", action="read")
    read (u, *) ps%pid, ps%comm, ps%state, ps%ppid, ps%pgrp, ps%session, &
      ps%tty_nr, ps%tpgid, ps%flags, ps%minflt, ps%cminflt, ps%majflt, &
      ps%cmajflt, ps%utime, ps%stime, ps%cutime, ps%cstime, ps%priority, &
      ps%nice, ps%num_threads, ps%itrealvalue, ps%starttime, ps%vsize, ps%rss, &
      ps%rsslim, ps%startcode, ps%endcode, ps%startstack, ps%kstkesp, ps%kstkeip, &
      ps%signal, ps%blocked, ps%sigignore, ps%sigcatch, ps%wchan, ps%nswap, ps%cnswap, &
      ps%exit_signal, ps%processor, ps%rt_priority, ps%policy, ps%delayacct_blkio_ticks, &
      ps%guest_time, ps%cguest_time
    close (u)
  end subroutine myInfoStat

  subroutine printMyInfoStat(ps, uo)
    class(psType)                    :: ps
    integer, intent(in), optional    :: uo

    integer :: u

    if (present(uo)) then
      u = uo
    else
      u = ou
    endif

    write (u, *) 'pid                   : ', ps%pid
    write (u, *) 'comm                  : ', ps%comm
    write (u, *) 'state                 : ', ps%state
    write (u, *) 'ppid                  : ', ps%ppid
    write (u, *) 'pgrp                  : ', ps%pgrp
    write (u, *) 'session               : ', ps%session
    write (u, *) 'tty_nr                : ', ps%tty_nr
    write (u, *) 'tpgid                 : ', ps%tpgid
    write (u, *) 'flags                 : ', ps%flags
    write (u, *) 'minflt                : ', ps%minflt
    write (u, *) 'cminflt               : ', ps%cminflt
    write (u, *) 'majflt                : ', ps%majflt
    write (u, *) 'cmajflt               : ', ps%cmajflt
    write (u, *) 'utime                 : ', ps%utime
    write (u, *) 'stime                 : ', ps%stime
    write (u, *) 'cutime                : ', ps%cutime
    write (u, *) 'cstime                : ', ps%cstime
    write (u, *) 'priority              : ', ps%priority
    write (u, *) 'nice                  : ', ps%nice
    write (u, *) 'num_threads           : ', ps%num_threads
    write (u, *) 'itrealvalue           : ', ps%itrealvalue
    write (u, *) 'starttime             : ', ps%starttime
    write (u, *) 'vsize                 : ', ps%vsize / 1024.0_wp
    write (u, *) 'rss                   : ', ps%rss / 1024.0_wp
    write (u, *) 'rsslim                : ', ps%rsslim
    write (u, *) 'startcode             : ', ps%startcode
    write (u, *) 'endcode               : ', ps%endcode
    write (u, *) 'startstack            : ', ps%startstack
    write (u, *) 'kstkesp               : ', ps%kstkesp
    write (u, *) 'kstkeip               : ', ps%kstkeip
    write (u, *) 'signal                : ', ps%signal
    write (u, *) 'blocked               : ', ps%blocked
    write (u, *) 'sigignore             : ', ps%sigignore
    write (u, *) 'sigcatch              : ', ps%sigcatch
    write (u, *) 'wchan                 : ', ps%wchan
    write (u, *) 'nswap                 : ', ps%nswap
    write (u, *) 'cnswap                : ', ps%cnswap
    write (u, *) 'exit_signal           : ', ps%exit_signal
    write (u, *) 'processor             : ', ps%processor
    write (u, *) 'rt_priority           : ', ps%rt_priority
    write (u, *) 'policy                : ', ps%policy
    write (u, *) 'delayacct_blkio_ticks : ', ps%delayacct_blkio_ticks
    write (u, *) 'guest_time            : ', ps%guest_time
    write (u, *) 'cguest_time           : ', ps%cguest_time

  end subroutine printMyInfoStat

  subroutine myInfoStatm(pgs)
    class(pgsType) :: pgs

    integer :: u

    open (newunit=u, file="/proc/self/statm", form="formatted", status="old", action="read")
    read (u, *) pgs%data
    close (u)
  end subroutine myInfoStatm

  subroutine printMyInfoStatm(pgs, uo)
    class(pgsType)                   :: pgs
    integer, intent(in), optional    :: uo

    integer :: u

    if (present(uo)) then
      u = uo
    else
      u = ou
    endif
    write (u, *) "Total pages           :", pgs%data(sz)
    write (u, *) "Resident pages        :", pgs%data(resident)
    write (u, *) "Shared pages          :", pgs%data(share)
    write (u, *) "Text pages            :", pgs%data(textrs)
    write (u, *) "Shared libs pages     :", pgs%data(librs)
    write (u, *) "Data pages            :", pgs%data(datars)
    write (u, *) "Dirty pages           :", pgs%data(dirty)
  end subroutine printMyInfoStatm

  subroutine myInfoStatus(pss)
    class(pssType) :: pss

    character(len=256) :: aux
    character(len=50)  :: a, label, val
    integer            :: io, u

    open (newunit=u, file="/proc/self/status", form="formatted", status="old", action="read")
    pss%data = 0
    do
      read (u, '(a)', iostat=io) aux
      if (io /= 0) exit
      read (aux, *) label, val
      select case (trim (label))
      case ("VmPeak:")
        read (val, '(i20)') pss%data(VmPeak)
        read (aux, *) a, a, pss%sunit
      case ("VmSize:")
        read (val, '(i20)') pss%data(VmSize)
      case ("VmLck:")
        read (val, '(i20)') pss%data(VmLck)
      case ("VmPin:")
        read (val, '(i20)') pss%data(VmPin)
      case ("VmHWM:")
        read (val, '(i20)') pss%data(VmHWM)
      case ("VmRSS:")
        read (val, '(i20)') pss%data(VmRSS)
      case ("RssAnon:")
        read (val, '(i20)') pss%data(RssAnon)
      case ("RssFile:")
        read (val, '(i20)') pss%data(RssFile)
      case ("RssShmem:")
        read (val, '(i20)') pss%data(RssShmem)
      case ("VmData:")
        read (val, '(i20)') pss%data(VmData)
      case ("VmStk:")
        read (val, '(i20)') pss%data(VmStk)
      case ("VmExe:")
        read (val, '(i20)') pss%data(VmExe)
      case ("VmLib:")
        read (val, '(i20)') pss%data(VmLib)
      case ("VmPTE:")
        read (val, '(i20)') pss%data(VmPTE)
      case ("VmSwap:")
        read (val, '(i20)') pss%data(VmSwap)
      case ("HugetlbPages:")
        read (val, '(i20)') pss%data(HugetlbPages)
      case ("Threads:")
        read (val, '(i20)') pss%data(Threads)
      case ("voluntary_ctxt_switches:")
        read (val, '(i20)') pss%data(voluntary_ctxt_switches)
      case ("nonvoluntary_ctxt_switches:")
        read (val, '(i20)') pss%data(nonvoluntary_ctxt_switches)
      end select
    end do
    close (u)
  end subroutine myInfoStatus

  subroutine printMyInfoStatus(pss, uo)
    class(pssType)                   :: pss
    integer, intent(in), optional    :: uo

    integer :: u

    if (present(uo)) then
      u = uo
    else
      u = ou
    endif
    write (u, *) "unit", pss%sunit
    write (u, *) "VmPeak                     : ", pss%data(VmPeak), pss%sunit
    write (u, *) "VmSize                     : ", pss%data(VmSize), pss%sunit
    write (u, *) "VmLck                      : ", pss%data(VmLck), pss%sunit
    write (u, *) "VmPin                      : ", pss%data(VmPin), pss%sunit
    write (u, *) "VmHWM                      : ", pss%data(VmHWM), pss%sunit
    write (u, *) "VmRSS                      : ", pss%data(VmRSS), pss%sunit
    write (u, *) "RssAnon                    : ", pss%data(RssAnon), pss%sunit
    write (u, *) "RssFile                    : ", pss%data(RssFile), pss%sunit
    write (u, *) "RssShmem                   : ", pss%data(RssShmem), pss%sunit
    write (u, *) "VmData                     : ", pss%data(VmData), pss%sunit
    write (u, *) "VmStk                      : ", pss%data(VmStk), pss%sunit
    write (u, *) "VmExe                      : ", pss%data(VmExe), pss%sunit
    write (u, *) "VmLib                      : ", pss%data(VmLib), pss%sunit
    write (u, *) "VmPTE                      : ", pss%data(VmPTE), pss%sunit
    write (u, *) "VmSwap                     : ", pss%data(VmSwap), pss%sunit
    write (u, *) "HugetlbPages               : ", pss%data(HugetlbPages), pss%sunit
    write (u, *) "Threads                    : ", pss%data(Threads)
    write (u, *) "voluntary_ctxt_switches    : ", pss%data(voluntary_ctxt_switches)
    write (u, *) "nonvoluntary_ctxt_switches : ", pss%data(nonvoluntary_ctxt_switches)
  end subroutine printMyInfoStatus

end module psInfo
