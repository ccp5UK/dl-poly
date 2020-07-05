module meminfo
  use, intrinsic :: iso_c_binding,   only: c_int
  use, intrinsic :: iso_fortran_env, only: li => int64,&
                                           wp => real64
  use comms, only : ROOT=>root_id,comms_type
  use mpi
  use psInfo
  use resources
  use table,                         only: print_header,&
                                           print_row,&
                                           hline

  implicit none

private
  interface reduce_sum
    module procedure reduce_sum_ia
  end interface

  interface reduce_min
    module procedure reduce_min_ia
  end interface

  interface reduce_max
    module procedure reduce_max_ia
  end interface


public  :: mem_picture
contains
  subroutine reduce_sum_ia(s, r, c, source, comm, ierr)
    integer(kind=li), intent(in)    :: s(:)
    integer(kind=li), intent(out)   :: r(:)
    integer, intent(in)             :: c, source
    integer, intent(in)      :: comm
    integer, intent(out)            :: ierr

    call MPI_Reduce(s, r, c, MPI_INT64_T, MPI_SUM, source, comm, ierr)
  end subroutine reduce_sum_ia

  subroutine reduce_min_ia(s, r, c, source, comm, ierr)
    integer(kind=li), intent(in)    :: s(:)
    integer(kind=li), intent(out)   :: r(:)
    integer, intent(in)             :: c, source
    integer, intent(in)      :: comm
    integer, intent(out)            :: ierr

    call MPI_Reduce(s, r, c, MPI_INT64_T, MPI_MIN, source, comm, ierr)
  end subroutine reduce_min_ia

  subroutine reduce_max_ia(s, r, c, source, comm, ierr)
    integer(kind=li), intent(in)    :: s(:)
    integer(kind=li), intent(out)   :: r(:)
    integer, intent(in)             :: c, source
    integer, intent(in)      :: comm
    integer, intent(out)            :: ierr

    call MPI_Reduce(s, r, c, MPI_INT64_T, MPI_MAX, source, comm, ierr)
  end subroutine reduce_max_ia

subroutine mem_picture(u,comm)
    integer, intent(in) :: u
    type(comms_type), intent(in) :: comm

  integer  :: ierror
  type(pssType) :: pss
  type(pgsType) :: pgs
  type(rusage_t) :: sa
  integer(c_int) :: x = 0, b

  integer, parameter :: m = 14
  integer(li) :: sbuff(m), smin(m), smax(m), avg(m)
  integer, parameter :: m2 = 7
  integer(li) :: sbuff2(m2), smin2(m2), smax2(m2), avg2(m2)
  integer, parameter :: m3 = 19
  integer(li) :: sbuff3(m3), smin3(m3), smax3(m3), avg3(m3)


  b = getrusage(x, sa)
  call pss%inspect()
  call pgs%inspect()

  !call print_res(sa)
  sbuff = 0_li
  sbuff(ru_maxrss) = sa%ru_maxrss
  sbuff(ru_minflt) = sa%ru_minflt
  sbuff(ru_majflt) = sa%ru_majflt
  sbuff(ru_inblock) = sa%ru_inblock
  sbuff(ru_oublock) = sa%ru_oublock
  sbuff(ru_nvcsw) = sa%ru_nvcsw
  sbuff(ru_nivcsw) = sa%ru_nivcsw
  call reduce_sum(sbuff, avg, m, ROOT, comm%comm, ierror)
  call reduce_min(sbuff, smin, m, ROOT, comm%comm, ierror)
  call reduce_max(sbuff, smax, m, ROOT, comm%comm, ierror)

  call reduce_sum(pgs%data, avg2, m2, ROOT, comm%comm, ierror)
  call reduce_min(pgs%data, smin2, m2, ROOT, comm%comm, ierror)
  call reduce_max(pgs%data, smax2, m2, ROOT, comm%comm, ierror)

  call reduce_sum(pss%data, avg3, m3, ROOT, comm%comm, ierror)
  call reduce_min(pss%data, smin3, m3, ROOT, comm%comm, ierror)
  call reduce_max(pss%data, smax3, m3, ROOT, comm%comm, ierror)
  if (comm%idnode == ROOT) then
    block
      integer :: w(4) = [20, 16, 16, 16]

      call print_header(["   Memory", "Min [kiB]", "Max [kiB]", 'Avg [kiB]'], w,u)
      call print_row(" resident", [smin(ru_maxrss), smax(ru_maxrss)], [avg(ru_maxrss) / real(comm%mxnode, wp)], w,u)
      call print_row("VmPeak", [smin3(VmPeak), smax3(VmPeak)], [avg2(VmPeak) / real(comm%mxnode, wp)], w,u)
      call print_row("VmSize", [smin3(VmSize), smax3(VmSize)], [avg3(VmSize) / real(comm%mxnode, wp)], w,u)
      call print_row("VmLck", [smin3(VmLck), smax3(VmLck)], [avg3(VmLck) / real(comm%mxnode, wp)], w,u)
      call print_row("VmPin", [smin3(VmPin), smax3(VmPin)], [avg3(VmPin) / real(comm%mxnode, wp)], w,u)
      call print_row("VmHWM", [smin3(VmHWM), smax3(VmHWM)], [avg3(VmHWM) / real(comm%mxnode, wp)], w,u)
      call print_row("VmRSS", [smin3(VmRSS), smax3(VmRSS)], [avg3(VmRSS) / real(comm%mxnode, wp)], w,u)
      call print_row("RssAnon", [smin3(RssAnon), smax3(RssAnon)], [avg3(RssAnon) / real(comm%mxnode, wp)], w,u)
      call print_row("RssFile", [smin3(RssFile), smax3(RssFile)], [avg3(RssFile) / real(comm%mxnode, wp)], w,u)
      call print_row("RssShmem", [smin3(RssShmem), smax3(RssShmem)], [avg3(RssShmem) / real(comm%mxnode, wp)], w,u)
      call print_row("VmData", [smin3(VmData), smax3(VmData)], [avg3(VmData) / real(comm%mxnode, wp)], w,u)
      call print_row("VmStk", [smin3(VmStk), smax3(VmStk)], [avg3(VmStk) / real(comm%mxnode, wp)], w,u)
      call print_row("VmExe", [smin3(VmExe), smax3(VmExe)], [avg3(VmExe) / real(comm%mxnode, wp)], w,u)
      call print_row("VmLib", [smin3(VmLib), smax3(VmLib)], [avg3(VmLib) / real(comm%mxnode, wp)], w,u)
      call print_row("VmPTE", [smin3(VmPTE), smax3(VmPTE)], [avg3(VmPTE) / real(comm%mxnode, wp)], w,u)
      call print_row("VmSwap", [smin3(VmSwap), smax3(VmSwap)], [avg3(VmSwap) / real(comm%mxnode, wp)], w,u)
      call print_row("HugetlbPages", [smin3(HugetlbPages), smax3(HugetlbPages)], [avg3(HugetlbPages) / real(comm%mxnode, wp)], w,u)
      write(u,'(a)')trim(hline(w))
      write (u, *)
      call print_header(["  Varia", "Min [#]", "Max [#]", 'Avg [#]'], w,u)
      call print_row("Minor Faults", [smin(ru_minflt), smax(ru_minflt)], [avg(ru_minflt) / real(comm%mxnode, wp)], w,u)
      call print_row("Major Faults", [smin(ru_majflt), smax(ru_majflt)], [avg(ru_majflt) / real(comm%mxnode, wp)], w,u)
      call print_row("input op", [smin(ru_inblock), smax(ru_inblock)], [avg(ru_inblock) / real(comm%mxnode, wp)], w,u)
      call print_row("output op", [smin(ru_oublock), smax(ru_oublock)], [avg(ru_oublock) / real(comm%mxnode, wp)], w,u)
      call print_row("volunt context sw", [smin(ru_nvcsw), smax(ru_nvcsw)], [avg(ru_nvcsw) / real(comm%mxnode, wp)], w,u)
      call print_row("involunt context sw", [smin(ru_nivcsw), smax(ru_nivcsw)], [avg(ru_nivcsw) / real(comm%mxnode, wp)], w,u)
      call print_row("Total pages", [smin2(sz), smax2(sz)], [avg2(sz) / real(comm%mxnode, wp)], w,u)
      call print_row("Resident pages", [smin2(resident), smax2(resident)], [avg2(resident) / real(comm%mxnode, wp)], w,u)
      call print_row("Shared pages", [smin2(share), smax2(share)], [avg2(share) / real(comm%mxnode, wp)], w,u)
      call print_row("Text pages", [smin2(textrs), smax2(textrs)], [avg2(textrs) / real(comm%mxnode, wp)], w,u)
      call print_row("Shared libs pages", [smin2(librs), smax2(librs)], [avg2(librs) / real(comm%mxnode, wp)], w,u)
      call print_row("Data pages", [smin2(datars), smax2(datars)], [avg2(datars) / real(comm%mxnode, wp)], w,u)
      call print_row("Dirty pages", [smin2(dirty), smax2(datars)], [avg2(dirty) / real(comm%mxnode, wp)], w,u)
      write(u,'(a)')trim(hline(w))

    end block
  end if
  end subroutine mem_picture
end module meminfo

