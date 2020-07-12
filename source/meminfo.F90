Module meminfo
  Use comms,                         Only: ROOT => root_id,&
                                           comms_type
  Use, Intrinsic :: iso_c_binding,   Only: c_int
  Use, Intrinsic :: iso_fortran_env, Only: li => int64,&
                                           wp => real64
  Use mpi
  Use psInfo
  Use resources
  Use table,                         Only: hline,&
                                           print_header,&
                                           print_row

  Implicit None

  Private
  Interface reduce_sum
    Module Procedure reduce_sum_ia
  End Interface

  Interface reduce_min
    Module Procedure reduce_min_ia
  End Interface

  Interface reduce_max
    Module Procedure reduce_max_ia
  End Interface

  Public  :: mem_picture
Contains
  Subroutine reduce_sum_ia(s, r, c, source, comm, ierr)
    Integer(kind=li), Intent(In   ) :: s(:)
    Integer(kind=li), Intent(  Out) :: r(:)
    Integer,          Intent(In   ) :: c, source, comm
    Integer,          Intent(  Out) :: ierr

    Call MPI_Reduce(s, r, c, MPI_INT64_T, MPI_SUM, source, comm, ierr)
  End Subroutine reduce_sum_ia

  Subroutine reduce_min_ia(s, r, c, source, comm, ierr)
    Integer(kind=li), Intent(In   ) :: s(:)
    Integer(kind=li), Intent(  Out) :: r(:)
    Integer,          Intent(In   ) :: c, source, comm
    Integer,          Intent(  Out) :: ierr

    Call MPI_Reduce(s, r, c, MPI_INT64_T, MPI_MIN, source, comm, ierr)
  End Subroutine reduce_min_ia

  Subroutine reduce_max_ia(s, r, c, source, comm, ierr)
    Integer(kind=li), Intent(In   ) :: s(:)
    Integer(kind=li), Intent(  Out) :: r(:)
    Integer,          Intent(In   ) :: c, source, comm
    Integer,          Intent(  Out) :: ierr

    Call MPI_Reduce(s, r, c, MPI_INT64_T, MPI_MAX, source, comm, ierr)
  End Subroutine reduce_max_ia

  Subroutine mem_picture(u, comm)
    Integer,          Intent(In   ) :: u
    Type(comms_type), Intent(In   ) :: comm

    Integer, Parameter :: m = 14, m2 = 7, m3 = 19

    Integer        :: ierror
    Integer(c_int) :: b, x = 0
    Integer(li)    :: avg(m), avg2(m2), avg3(m3), sbuff(m), smax(m), smax2(m2), smax3(m3), &
                      smin(m), smin2(m2), smin3(m3)
    Type(pgsType)  :: pgs
    Type(pssType)  :: pss
    Type(rusage_t) :: sa

    b = getrusage(x, sa)
    Call pss%inspect()
    Call pgs%inspect()

    !call print_res(sa)
    sbuff = 0_li
    sbuff(ru_maxrss) = sa%ru_maxrss
    sbuff(ru_minflt) = sa%ru_minflt
    sbuff(ru_majflt) = sa%ru_majflt
    sbuff(ru_inblock) = sa%ru_inblock
    sbuff(ru_oublock) = sa%ru_oublock
    sbuff(ru_nvcsw) = sa%ru_nvcsw
    sbuff(ru_nivcsw) = sa%ru_nivcsw
    Call reduce_sum(sbuff, avg, m, ROOT, comm%comm, ierror)
    Call reduce_min(sbuff, smin, m, ROOT, comm%comm, ierror)
    Call reduce_max(sbuff, smax, m, ROOT, comm%comm, ierror)

    Call reduce_sum(pgs%data, avg2, m2, ROOT, comm%comm, ierror)
    Call reduce_min(pgs%data, smin2, m2, ROOT, comm%comm, ierror)
    Call reduce_max(pgs%data, smax2, m2, ROOT, comm%comm, ierror)

    Call reduce_sum(pss%data, avg3, m3, ROOT, comm%comm, ierror)
    Call reduce_min(pss%data, smin3, m3, ROOT, comm%comm, ierror)
    Call reduce_max(pss%data, smax3, m3, ROOT, comm%comm, ierror)
    If (comm%idnode == ROOT) Then
      block
        Integer :: w(4) = [20, 16, 16, 16]

        Call print_header(["   Memory", "Min [kiB]", "Max [kiB]", 'Avg [kiB]'], w, u)
        Call print_row(" resident", [smin(ru_maxrss), smax(ru_maxrss)], [avg(ru_maxrss) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmPeak", [smin3(VmPeak), smax3(VmPeak)], [avg2(VmPeak) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmSize", [smin3(VmSize), smax3(VmSize)], [avg3(VmSize) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmLck", [smin3(VmLck), smax3(VmLck)], [avg3(VmLck) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmPin", [smin3(VmPin), smax3(VmPin)], [avg3(VmPin) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmHWM", [smin3(VmHWM), smax3(VmHWM)], [avg3(VmHWM) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmRSS", [smin3(VmRSS), smax3(VmRSS)], [avg3(VmRSS) / Real(comm%mxnode, wp)], w, u)
        Call print_row("RssAnon", [smin3(RssAnon), smax3(RssAnon)], [avg3(RssAnon) / Real(comm%mxnode, wp)], w, u)
        Call print_row("RssFile", [smin3(RssFile), smax3(RssFile)], [avg3(RssFile) / Real(comm%mxnode, wp)], w, u)
        Call print_row("RssShmem", [smin3(RssShmem), smax3(RssShmem)], [avg3(RssShmem) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmData", [smin3(VmData), smax3(VmData)], [avg3(VmData) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmStk", [smin3(VmStk), smax3(VmStk)], [avg3(VmStk) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmExe", [smin3(VmExe), smax3(VmExe)], [avg3(VmExe) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmLib", [smin3(VmLib), smax3(VmLib)], [avg3(VmLib) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmPTE", [smin3(VmPTE), smax3(VmPTE)], [avg3(VmPTE) / Real(comm%mxnode, wp)], w, u)
        Call print_row("VmSwap", [smin3(VmSwap), smax3(VmSwap)], [avg3(VmSwap) / Real(comm%mxnode, wp)], w, u)
      Call print_row("HugetlbPages", [smin3(HugetlbPages), smax3(HugetlbPages)], [avg3(HugetlbPages) / Real(comm%mxnode, wp)], w, u)
        Write (u, '(a)') Trim(hline(w))
        Write (u, *)
        Call print_header(["  Varia", "Min [#]", "Max [#]", 'Avg [#]'], w, u)
        Call print_row("Minor Faults", [smin(ru_minflt), smax(ru_minflt)], [avg(ru_minflt) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Major Faults", [smin(ru_majflt), smax(ru_majflt)], [avg(ru_majflt) / Real(comm%mxnode, wp)], w, u)
        Call print_row("input op", [smin(ru_inblock), smax(ru_inblock)], [avg(ru_inblock) / Real(comm%mxnode, wp)], w, u)
        Call print_row("output op", [smin(ru_oublock), smax(ru_oublock)], [avg(ru_oublock) / Real(comm%mxnode, wp)], w, u)
        Call print_row("volunt context sw", [smin(ru_nvcsw), smax(ru_nvcsw)], [avg(ru_nvcsw) / Real(comm%mxnode, wp)], w, u)
        Call print_row("involunt context sw", [smin(ru_nivcsw), smax(ru_nivcsw)], [avg(ru_nivcsw) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Total pages", [smin2(sz), smax2(sz)], [avg2(sz) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Resident pages", [smin2(resident), smax2(resident)], [avg2(resident) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Shared pages", [smin2(share), smax2(share)], [avg2(share) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Text pages", [smin2(textrs), smax2(textrs)], [avg2(textrs) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Shared libs pages", [smin2(librs), smax2(librs)], [avg2(librs) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Data pages", [smin2(datars), smax2(datars)], [avg2(datars) / Real(comm%mxnode, wp)], w, u)
        Call print_row("Dirty pages", [smin2(dirty), smax2(datars)], [avg2(dirty) / Real(comm%mxnode, wp)], w, u)
        Write (u, '(a)') Trim(hline(w))

      End block
    End If
  End Subroutine mem_picture
End Module meminfo

