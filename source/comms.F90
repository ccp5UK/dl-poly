Module comms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for global communication routines and functions
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  ! contrib   - a.m.elena march 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,sp,dp,qp
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif

  Implicit None

  !  Include 'mpif.h'  ! Needed instead "Use mpi" for some machines
  !  Include 'mpiof.h' ! Needed for ScaliMPI

  ! l_fast is controlled via gsync and affects gcheck - global safety checks

  Integer, Save :: wp_mpi = 0 


  ! MPI-I/O representation

  Character( Len = 6 ), Parameter :: datarep = 'native'

  Integer,                                           Public :: mpi_ver     = -1, &
    mpi_subver  = -1
  Character( Len = MPI_MAX_PROCESSOR_NAME ),         Public :: proc_name   = "*"
#ifndef OLDMPI
  Character( Len = MPI_MAX_LIBRARY_VERSION_STRING ), Public :: lib_version = "*"
#endif

  ! Message tags

  Integer, Parameter :: Deport_tag    = 1100, &
    Export_tag    = 1111, &
    Revive_tag    = 1122, &
    PassUnit_tag  = 1133, &
    UpdShUnit_tag = 1144, &
    SysExpand_tag = 1155, &
    WriteConf_tag = 1166, &
    Traject_tag   = 1177, &
    Spread_tag    = 1188, &
    DpdVExp_tag   = 1199, &
    MetLdExp_tag  = 2200, &
    ExpMplRM_tag  = 2211, &
    ExchgGrid_tag = 2222, &
    DefRWrite_tag = 2233, &
    DefExport_tag = 2244, &
    DefWrite_tag  = 2255, &
    RsdWrite_tag  = 2266, &
    MsdWrite_tag  = 2277, &
    Grid1_tag     = 3300, &
    Grid2_tag     = 3311, &
    Grid3_tag     = 3322, &
    Grid4_tag     = 3333

  Type, Public :: comms_type
    Integer               :: ierr 
    Integer               :: request  
    Integer               :: status(1:MPI_STATUS_SIZE) = 0
    Integer               :: comm
    Integer               :: idnode = 0
    Integer               :: mxnode = 1
    Logical               :: l_fast
  End Type 

  Public :: init_comms, exit_comms, abort_comms, &
    gsync, gcheck, gsum, gmax, gtime

  Interface gcheck
    Module Procedure gcheck_vector
    Module Procedure gcheck_scalar
  End Interface !gcheck

  Interface gsum
    Module Procedure gisum_vector
    Module Procedure gisum_scalar

    Module Procedure grsum_matrix
    Module Procedure grsum_vector
    Module Procedure grsum_scalar
  End Interface !gsum

  Interface gmax
    Module Procedure gimax_vector
    Module Procedure gimax_scalar

    Module Procedure grmax_vector
    Module Procedure grmax_scalar
  End Interface !gmax

  Interface gmin
    Module Procedure gimin_vector
    Module Procedure gimin_scalar

    Module Procedure grmin_vector
    Module Procedure grmin_scalar
  End Interface !gmin

  Interface gbcast
    Module procedure gbcast_integer
    Module procedure gbcast_real
  End Interface !gbcast

Contains

  Subroutine init_comms(comm)
    Type(comms_type), Intent (InOut) :: comm

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for initialising the communication harness
    ! determining the MPI precision, and node identification and count
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    ! contrib   - a.m.elena march 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer :: lversion, lname

    Call MPI_INIT(comm%ierr)

    Call MPI_COMM_DUP(MPI_COMM_WORLD,comm%comm,comm%ierr)

    ! use iso_fortran_env
    If      (wp == sp)  Then
      wp_mpi = MPI_REAL
    Else If (wp == dp)  Then
      ! MPI_REAL8 is apparently not in the strict MPI2 standard
      ! It is just an optional data type in the FORTRAN Bindings
      wp_mpi = MPI_DOUBLE_PRECISION
    Else If (wp == qp) Then
      wp_mpi = MPI_REAL16
    Else
      Call error(1000)
    End If

    Call MPI_COMM_RANK(comm%comm, comm%idnode, comm%ierr)
    Call MPI_COMM_SIZE(comm%comm, comm%mxnode, comm%ierr)

#ifndef OLDMPI
    Call MPI_GET_PROCESSOR_NAME(proc_name,lname, comm%ierr)
    Call MPI_GET_VERSION(mpi_ver,mpi_subver, comm%ierr)
    Call MPI_GET_LIBRARY_VERSION(lib_version,lversion, comm%ierr)
#endif

  End Subroutine init_comms

  Subroutine exit_comms(comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for exiting the communication harness
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent (InOut) :: comm

    Call MPI_FINALIZE(comm%ierr)

  End Subroutine exit_comms

  Subroutine abort_comms(comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for aborting the communication harness
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2008
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut) :: comm
    Call MPI_ABORT(comm%comm, comm%idnode, comm%ierr)

  End Subroutine abort_comms

  Subroutine gsync(comm, lfast)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 barrier/synchronization routine
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut) :: comm
    Logical, Optional :: lfast

    If (Present(lfast)) comm%l_fast = lfast

    If (comm%mxnode == 1) Return

    Call MPI_BARRIER(comm%comm,comm%ierr)

  End Subroutine gsync

  Subroutine gcheck_vector(comm, aaa, enforce)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - boolean vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut) :: comm
    Logical, Dimension( : ), Intent( InOut )           :: aaa
    Character( Len = *),     Intent( In    ), Optional :: enforce

    Integer                                  :: n_l,n_u,n_s,fail
    Logical, Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Then
      Return
    Else
      If (.not.Present(enforce)) Then
        If (comm%l_fast) Return
      End If
    End If

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1001)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_LOGICAL,MPI_LAND,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1002)

  End Subroutine gcheck_vector

  Subroutine gcheck_scalar(comm, aaa, enforce)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - boolean scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)               :: comm
    Logical,             Intent( InOut )           :: aaa
    Character( Len = *), Intent( In    ), Optional :: enforce

    Logical :: bbb

    If (comm%mxnode == 1) Then
      Return
    Else
      If (.not.Present(enforce)) Then
        If (comm%l_fast) Return
      End If
    End If

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_LOGICAL,MPI_LAND,comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine gcheck_scalar

  Subroutine gisum_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - integer vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)         :: comm
    Integer, Dimension( : ), Intent( InOut ) :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Integer, Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1003)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_INTEGER,MPI_SUM,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1004)

  End Subroutine gisum_vector

  Subroutine gisum_scalar(comm,aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - integer scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)         :: comm
    Integer, Intent( InOut )                 :: aaa

    Integer                  :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_INTEGER,MPI_SUM,comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine gisum_scalar

  Subroutine grsum_matrix(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - working precision vector
    !                                         version
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush december 2009
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                      :: comm
    Real( Kind = wp ), Dimension( :, : ), Intent( InOut ) :: aaa

    Integer                                               :: n_l1,n_u1,n_s,fail
    Integer                                               :: n_l2,n_u2
    Real( Kind = wp ), Dimension( :, : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l1 = Lbound(aaa, Dim = 1)
    n_u1 = Ubound(aaa, Dim = 1)

    n_l2 = Lbound(aaa, Dim = 2)
    n_u2 = Ubound(aaa, Dim = 2)

    fail = 0
    Allocate (bbb(n_l1:n_u1,n_l2:n_u2), Stat = fail)
    If (fail > 0) Call error(1048)

    n_s = Size(aaa)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_SUM,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1049)

  End Subroutine grsum_matrix

  Subroutine grsum_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - working precision vector
    !                                         version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Real( Kind = wp ), Dimension( : ), Intent( InOut ) :: aaa

    Integer                                            :: n_l,n_u,n_s,fail
    Real( Kind = wp ), Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1005)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_SUM,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1006)

  End Subroutine grsum_vector

  Subroutine grsum_scalar(comm,aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - working precision scalar
    !                                         version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)   :: comm
    Real( Kind = wp ), Intent( InOut ) :: aaa

    Real( Kind = wp )                  :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,wp_mpi,MPI_SUM,comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine grsum_scalar

  Subroutine gimax_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global maximum subroutine - integer vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)         :: comm
    Integer, Dimension( : ), Intent( InOut ) :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Integer, Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1007)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_INTEGER,MPI_MAX,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1008)

  End Subroutine gimax_vector

  Subroutine gimax_scalar(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global maximum subroutine - integer scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type(comms_type), Intent (InOut)         :: comm
    Integer, Intent( InOut )                 :: aaa

    Integer                  :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_INTEGER,MPI_MAX,comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine gimax_scalar

  Subroutine grmax_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global maximum subroutine - working precision vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Real( Kind = wp ), Dimension( : ), Intent( InOut ) :: aaa

    Integer                                            :: n_l,n_u,n_s,fail
    Real( Kind = wp ), Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1009)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_MAX,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1010)

  End Subroutine grmax_vector

  Subroutine grmax_scalar(comm,aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global maximum subroutine - working precision scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Real( Kind = wp ), Intent( InOut )                 :: aaa

    Real( Kind = wp )                  :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,wp_mpi,MPI_MAX, comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine grmax_scalar

  Subroutine gimin_vector(comm,aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global minimum subroutine - integer vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2007
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Integer, Dimension( : ), Intent( InOut )           :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Integer, Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1044)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_INTEGER,MPI_MIN,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1045)

  End Subroutine gimin_vector

  Subroutine gimin_scalar(comm,aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global minimum subroutine - integer scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2007
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Integer, Intent( InOut )                           :: aaa

    Integer                  :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_INTEGER,MPI_MIN,comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine gimin_scalar

  Subroutine grmin_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global minimum subroutine - working precision vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2007
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Real( Kind = wp ), Dimension( : ), Intent( InOut ) :: aaa

    Integer                                            :: n_l,n_u,n_s,fail
    Real( Kind = wp ), Dimension( : ), Allocatable     :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1046)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_MIN,comm%comm,comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1047)

  End Subroutine grmin_vector

  Subroutine grmin_scalar(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global minimum subroutine - working precision scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2007
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent (InOut)                   :: comm
    Real( Kind = wp ), Intent( InOut )                 :: aaa

    Real( Kind = wp )                  :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,wp_mpi,MPI_MIN,comm%comm,comm%ierr)

    aaa = bbb

  End Subroutine grmin_scalar

  Subroutine gbcast_integer(comm,vec,root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast an integer array subroutine 
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena may 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer, Intent( InOut )                           :: vec(:)
    Integer, Intent( In    )                           :: root
    Type(comms_type), Intent (InOut)                   :: comm

    Integer                                            :: n_l,n_u,n_s,fail

    If (comm%mxnode == 1) Return
    n_l = Lbound(vec, Dim = 1)
    n_u = Ubound(vec, Dim = 1)
    n_s = Size(vec, Dim = 1)

    Call MPI_BCAST(vec(n_l:n_u), n_s, MPI_INTEGER, root, comm%comm, comm%ierr)

  End Subroutine gbcast_integer

  Subroutine gbcast_real(comm,vec,root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast a real array subroutine 
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena march 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind= wp ), Intent( InOut )                   :: vec(:)
    Integer, Intent( In    )                           :: root
    Type(comms_type), Intent (InOut)                   :: comm

    Integer                                            :: n_l,n_u,n_s,fail

    If (comm%mxnode == 1) Return
    n_l = Lbound(vec, Dim = 1)
    n_u = Ubound(vec, Dim = 1)
    n_s = Size(vec, Dim = 1)

    Call MPI_BCAST(vec(n_l:n_u), n_s, wp_mpi, root, comm%comm, comm%ierr)
  End Subroutine gbcast_real


  Subroutine gtime(time)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 timing routine
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent(   Out )                 :: time

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: tzero

    If (newjob) Then
      newjob = .false.

      tzero = MPI_WTIME()
      time = 0.0_wp
    Else
      time = MPI_WTIME()-tzero
    End If

  End Subroutine gtime

End Module comms
