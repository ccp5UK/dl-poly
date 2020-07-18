Module psInfo
  Use, Intrinsic :: iso_c_binding,   Only: C_LONG
  Use, Intrinsic :: iso_fortran_env, Only: li => int64,&
                                           ou => output_unit,&
                                           wp => real64

  Implicit None
  Private
  Type, Public :: psType
    Integer(kind=li)   :: pid
    Character(len=16) :: comm
    Character(len=1)  :: state
    Integer(kind=li)   :: ppid
    Integer(kind=li)   :: pgrp
    Integer(kind=li)   :: session
    Integer(kind=li)   :: tty_nr
    Integer(kind=li)   :: tpgid
    Integer(kind=li)   :: flags
    Real(kind=wp)      :: minflt
    Real(kind=wp)      :: cminflt
    Real(kind=wp)      :: majflt
    Real(kind=wp)      :: cmajflt
    Real(kind=wp)      :: utime
    Real(kind=wp)      :: stime
    Real(kind=wp)      :: cutime
    Real(kind=wp)      :: cstime
    Integer(kind=li)   :: priority
    Integer(kind=li)   :: nice
    Integer(kind=li)   :: num_threads
    Real(kind=wp)      :: itrealvalue
    Real(kind=wp)      :: starttime
    Integer(kind=c_long)      :: vsize
    Real(kind=wp)      :: rss
    Real(kind=wp)      :: rsslim
    Real(kind=wp)      :: startcode
    Real(kind=wp)      :: endcode
    Real(kind=wp)      :: startstack
    Real(kind=wp)      :: kstkesp
    Real(kind=wp)      :: kstkeip
    Integer(kind=li)   :: signal
    Real(kind=wp)      :: blocked
    Integer(kind=li)   :: sigignore
    Integer(kind=li)   :: sigcatch
    Real(kind=wp)      :: wchan
    Real(kind=wp)      :: nswap
    Real(kind=wp)      :: cnswap
    Integer(kind=li)   :: exit_signal
    Integer(kind=li)   :: processor
    Integer(kind=li)   :: rt_priority
    Integer(kind=li)   :: policy
    Real(kind=wp)      :: delayacct_blkio_ticks
    Real(kind=wp)      :: guest_time
    Real(kind=wp)      :: cguest_time
  Contains
    Private
    Procedure, Public :: inspect => myInfoStat
    Procedure, Public :: Print => printMyInfoStat
  End Type psType

  Integer, Parameter, Public   :: sz = 1 ! total # of pages of memory
  Integer, Parameter, Public   :: resident = 2 ! number of resident set (non-swapped) pages (4k)
  Integer, Parameter, Public   :: share = 3 ! number of pages of shared (mmap'd) memory
  Integer, Parameter, Public   :: textrs = 4 ! text resident size
  Integer, Parameter, Public   :: librs = 5 ! shared libs resident size
  Integer, Parameter, Public   :: datars = 6 ! data resident size
  Integer, Parameter, Public   :: dirty = 7 ! dirty pages
  Type, Public ::pgsType
    Integer(kind=li) :: Data(7)

  Contains
    Private
    Procedure, Public :: inspect => myInfoStatm
    Procedure, Public :: Print => printMyInfoStatm
  End Type pgsType

  Integer(kind=4), Parameter, Public   :: VmPeak = 1
  Integer(kind=4), Parameter, Public   :: VmSize = 2
  Integer(kind=4), Parameter, Public   :: VmLck = 3
  Integer(kind=4), Parameter, Public   :: VmPin = 4
  Integer(kind=4), Parameter, Public   :: VmHWM = 5
  Integer(kind=4), Parameter, Public   :: VmRSS = 6
  Integer(kind=4), Parameter, Public   :: RssAnon = 7
  Integer(kind=4), Parameter, Public   :: RssFile = 8
  Integer(kind=4), Parameter, Public   :: RssShmem = 9
  Integer(kind=4), Parameter, Public   :: VmData = 10
  Integer(kind=4), Parameter, Public   :: VmStk = 11
  Integer(kind=4), Parameter, Public   :: VmExe = 12
  Integer(kind=4), Parameter, Public   :: VmLib = 13
  Integer(kind=4), Parameter, Public   :: VmPTE = 14
  Integer(kind=4), Parameter, Public   :: VmSwap = 15
  Integer(kind=4), Parameter, Public   :: HugetlbPages = 16
  Integer(kind=4), Parameter, Public   :: Threads = 17
  Integer(kind=4), Parameter, Public   :: voluntary_ctxt_switches = 18
  Integer(kind=4), Parameter, Public   :: nonvoluntary_ctxt_switches = 19

  Type, Public :: pssType
    Integer :: n = 19
    Integer(kind=li) :: Data(19)
    Character(len=3)  :: sunit
  Contains
    Private
    Procedure, Public :: inspect => myInfoStatus
    Procedure, Public :: Print => printMyInfoStatus

  End Type pssType

Contains

  Subroutine myInfoStat(ps)
    class(psType) :: ps

    Integer :: u

    Open (newunit=u, file="/proc/self/stat", form="formatted", status="old", action="read")
    Read (u, *) ps%pid, ps%comm, ps%state, ps%ppid, ps%pgrp, ps%session, &
      ps%tty_nr, ps%tpgid, ps%flags, ps%minflt, ps%cminflt, ps%majflt, &
      ps%cmajflt, ps%utime, ps%stime, ps%cutime, ps%cstime, ps%priority, &
      ps%nice, ps%num_threads, ps%itrealvalue, ps%starttime, ps%vsize, ps%rss, &
      ps%rsslim, ps%startcode, ps%endcode, ps%startstack, ps%kstkesp, ps%kstkeip, &
      ps%signal, ps%blocked, ps%sigignore, ps%sigcatch, ps%wchan, ps%nswap, ps%cnswap, &
      ps%exit_signal, ps%processor, ps%rt_priority, ps%policy, ps%delayacct_blkio_ticks, &
      ps%guest_time, ps%cguest_time
    Close (u)
  End Subroutine myInfoStat

  Subroutine printMyInfoStat(ps, uo)
    class(psType)                    :: ps
    Integer, Optional, Intent(In   ) :: uo

    Integer :: u

    If (Present(uo)) Then
      u = uo
    Else
      u = ou
    Endif

    Write (u, *) 'pid                   : ', ps%pid
    Write (u, *) 'comm                  : ', ps%comm
    Write (u, *) 'state                 : ', ps%state
    Write (u, *) 'ppid                  : ', ps%ppid
    Write (u, *) 'pgrp                  : ', ps%pgrp
    Write (u, *) 'session               : ', ps%session
    Write (u, *) 'tty_nr                : ', ps%tty_nr
    Write (u, *) 'tpgid                 : ', ps%tpgid
    Write (u, *) 'flags                 : ', ps%flags
    Write (u, *) 'minflt                : ', ps%minflt
    Write (u, *) 'cminflt               : ', ps%cminflt
    Write (u, *) 'majflt                : ', ps%majflt
    Write (u, *) 'cmajflt               : ', ps%cmajflt
    Write (u, *) 'utime                 : ', ps%utime
    Write (u, *) 'stime                 : ', ps%stime
    Write (u, *) 'cutime                : ', ps%cutime
    Write (u, *) 'cstime                : ', ps%cstime
    Write (u, *) 'priority              : ', ps%priority
    Write (u, *) 'nice                  : ', ps%nice
    Write (u, *) 'num_threads           : ', ps%num_threads
    Write (u, *) 'itrealvalue           : ', ps%itrealvalue
    Write (u, *) 'starttime             : ', ps%starttime
    Write (u, *) 'vsize                 : ', ps%vsize / 1024.0_wp
    Write (u, *) 'rss                   : ', ps%rss / 1024.0_wp
    Write (u, *) 'rsslim                : ', ps%rsslim
    Write (u, *) 'startcode             : ', ps%startcode
    Write (u, *) 'endcode               : ', ps%endcode
    Write (u, *) 'startstack            : ', ps%startstack
    Write (u, *) 'kstkesp               : ', ps%kstkesp
    Write (u, *) 'kstkeip               : ', ps%kstkeip
    Write (u, *) 'signal                : ', ps%signal
    Write (u, *) 'blocked               : ', ps%blocked
    Write (u, *) 'sigignore             : ', ps%sigignore
    Write (u, *) 'sigcatch              : ', ps%sigcatch
    Write (u, *) 'wchan                 : ', ps%wchan
    Write (u, *) 'nswap                 : ', ps%nswap
    Write (u, *) 'cnswap                : ', ps%cnswap
    Write (u, *) 'exit_signal           : ', ps%exit_signal
    Write (u, *) 'processor             : ', ps%processor
    Write (u, *) 'rt_priority           : ', ps%rt_priority
    Write (u, *) 'policy                : ', ps%policy
    Write (u, *) 'delayacct_blkio_ticks : ', ps%delayacct_blkio_ticks
    Write (u, *) 'guest_time            : ', ps%guest_time
    Write (u, *) 'cguest_time           : ', ps%cguest_time

  End Subroutine printMyInfoStat

  Subroutine myInfoStatm(pgs)
    class(pgsType) :: pgs

    Integer :: u

    Open (newunit=u, file="/proc/self/statm", form="formatted", status="old", action="read")
    Read (u, *) pgs%data
    Close (u)
  End Subroutine myInfoStatm

  Subroutine printMyInfoStatm(pgs, uo)
    class(pgsType)                   :: pgs
    Integer, Optional, Intent(In   ) :: uo

    Integer :: u

    If (Present(uo)) Then
      u = uo
    Else
      u = ou
    Endif
    Write (u, *) "Total pages           :", pgs%data(sz)
    Write (u, *) "Resident pages        :", pgs%data(resident)
    Write (u, *) "Shared pages          :", pgs%data(share)
    Write (u, *) "Text pages            :", pgs%data(textrs)
    Write (u, *) "Shared libs pages     :", pgs%data(librs)
    Write (u, *) "Data pages            :", pgs%data(datars)
    Write (u, *) "Dirty pages           :", pgs%data(dirty)
  End Subroutine printMyInfoStatm

  Subroutine myInfoStatus(pss)
    class(pssType) :: pss

    Character(len=256) :: aux
    Character(len=50)  :: a, label, val
    Integer            :: io, u

    Open (newunit=u, file="/proc/self/status", form="formatted", status="old", action="read")
    pss%data = 0
    Do
      Read (u, '(a)', iostat=io) aux
      If (io /= 0) Exit
      Read (aux, *) label, val
      Select Case (Trim (label))
      Case ("VmPeak:")
        Read (val, '(i20)') pss%data(VmPeak)
        Read (aux, *) a, a, pss%sunit
      Case ("VmSize:")
        Read (val, '(i20)') pss%data(VmSize)
      Case ("VmLck:")
        Read (val, '(i20)') pss%data(VmLck)
      Case ("VmPin:")
        Read (val, '(i20)') pss%data(VmPin)
      Case ("VmHWM:")
        Read (val, '(i20)') pss%data(VmHWM)
      Case ("VmRSS:")
        Read (val, '(i20)') pss%data(VmRSS)
      Case ("RssAnon:")
        Read (val, '(i20)') pss%data(RssAnon)
      Case ("RssFile:")
        Read (val, '(i20)') pss%data(RssFile)
      Case ("RssShmem:")
        Read (val, '(i20)') pss%data(RssShmem)
      Case ("VmData:")
        Read (val, '(i20)') pss%data(VmData)
      Case ("VmStk:")
        Read (val, '(i20)') pss%data(VmStk)
      Case ("VmExe:")
        Read (val, '(i20)') pss%data(VmExe)
      Case ("VmLib:")
        Read (val, '(i20)') pss%data(VmLib)
      Case ("VmPTE:")
        Read (val, '(i20)') pss%data(VmPTE)
      Case ("VmSwap:")
        Read (val, '(i20)') pss%data(VmSwap)
      Case ("HugetlbPages:")
        Read (val, '(i20)') pss%data(HugetlbPages)
      Case ("Threads:")
        Read (val, '(i20)') pss%data(Threads)
      Case ("voluntary_ctxt_switches:")
        Read (val, '(i20)') pss%data(voluntary_ctxt_switches)
      Case ("nonvoluntary_ctxt_switches:")
        Read (val, '(i20)') pss%data(nonvoluntary_ctxt_switches)
      End Select
    End Do
    Close (u)
  End Subroutine myInfoStatus

  Subroutine printMyInfoStatus(pss, uo)
    class(pssType)                   :: pss
    Integer, Optional, Intent(In   ) :: uo

    Integer :: u

    If (Present(uo)) Then
      u = uo
    Else
      u = ou
    Endif
    Write (u, *) "unit", pss%sunit
    Write (u, *) "VmPeak                     : ", pss%data(VmPeak), pss%sunit
    Write (u, *) "VmSize                     : ", pss%data(VmSize), pss%sunit
    Write (u, *) "VmLck                      : ", pss%data(VmLck), pss%sunit
    Write (u, *) "VmPin                      : ", pss%data(VmPin), pss%sunit
    Write (u, *) "VmHWM                      : ", pss%data(VmHWM), pss%sunit
    Write (u, *) "VmRSS                      : ", pss%data(VmRSS), pss%sunit
    Write (u, *) "RssAnon                    : ", pss%data(RssAnon), pss%sunit
    Write (u, *) "RssFile                    : ", pss%data(RssFile), pss%sunit
    Write (u, *) "RssShmem                   : ", pss%data(RssShmem), pss%sunit
    Write (u, *) "VmData                     : ", pss%data(VmData), pss%sunit
    Write (u, *) "VmStk                      : ", pss%data(VmStk), pss%sunit
    Write (u, *) "VmExe                      : ", pss%data(VmExe), pss%sunit
    Write (u, *) "VmLib                      : ", pss%data(VmLib), pss%sunit
    Write (u, *) "VmPTE                      : ", pss%data(VmPTE), pss%sunit
    Write (u, *) "VmSwap                     : ", pss%data(VmSwap), pss%sunit
    Write (u, *) "HugetlbPages               : ", pss%data(HugetlbPages), pss%sunit
    Write (u, *) "Threads                    : ", pss%data(Threads)
    Write (u, *) "voluntary_ctxt_switches    : ", pss%data(voluntary_ctxt_switches)
    Write (u, *) "nonvoluntary_ctxt_switches : ", pss%data(nonvoluntary_ctxt_switches)
  End Subroutine printMyInfoStatus

End Module psInfo
