Module comms_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for global communication routines and functions
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use mpi!_module

  Implicit None

!  Include 'mpif.h'  ! Needed instead "Use mpi" for some machines
!  Include 'mpiof.h' ! Needed for ScaliMPI

  Integer, Save :: wp_mpi = 0 , &
                   idnode = 0 , &
                   mxnode = 1

  Integer, Save :: ierr                      = 0 , &
                   request                   = 0 , &
                   dlp_comm_world            = 0 , &
                   status(1:MPI_STATUS_SIZE) = 0

! MPI-I/O representation

  Character( Len = 6), Parameter :: datarep = 'native'

! Message tags

  Integer, Parameter :: Deport_tag   = 1100, &
                        Export_tag   = 1111, &
                        Revive_tag   = 1122, &
                        Traject_tag  = 1133, &
                        Fftcomms_tag = 1144, &
                        Passunit_tag = 1155, &
                        Passpmf_tag  = 1166, &
                        Updshun_tag  = 1177, &
                        Updpmfun_tag = 1188, &
                        Metldexp_tag = 1199, &
                        DefExport_tag= 2200, &
                        DefWrite_tag = 2211

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

Contains

  Subroutine init_comms()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for initialising the communication harness
! determining the MPI precision, and node identification and count
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Call MPI_INIT(ierr)

    Call MPI_COMM_DUP(MPI_COMM_WORLD,dlp_comm_world,ierr)

    If      (wp == Selected_Real_Kind(6,30))  Then
       wp_mpi = MPI_REAL
    Else If (wp == Selected_Real_Kind(14,300))  Then
! MPI_REAL8 is apparently not in the strict MPI2 standard
! It is just an optional data type in the FORTRAN Bindings
       wp_mpi = MPI_DOUBLE_PRECISION
    Else If (wp == Selected_Real_Kind(30,300)) Then
       wp_mpi = MPI_REAL16
    Else
       Call error(1000)
    End If

    Call MPI_COMM_RANK(dlp_comm_world, idnode, ierr)
    Call MPI_COMM_SIZE(dlp_comm_world, mxnode, ierr)

  End Subroutine init_comms

  Subroutine exit_comms()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for exiting the communication harness
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Call MPI_FINALIZE(ierr)

    Stop

  End Subroutine exit_comms

  Subroutine abort_comms()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for aborting the communication harness
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Call MPI_ABORT(dlp_comm_world, idnode, ierr)

    Stop

  End Subroutine abort_comms

  Subroutine gsync()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 barrier/synchronization routine
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    If (mxnode == 1) Return

    Call MPI_BARRIER(dlp_comm_world,ierr)

  End Subroutine gsync

  Subroutine gcheck_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - boolean vector version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Logical, Dimension( : ), Intent( InOut ) :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Logical, Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1001)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_LOGICAL,MPI_LAND,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1002)

  End Subroutine gcheck_vector

  Subroutine gcheck_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - boolean scalar version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Logical, Intent( InOut ) :: aaa

    Logical                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_LOGICAL,MPI_LAND,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine gcheck_scalar

  Subroutine gisum_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - integer vector version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Dimension( : ), Intent( InOut ) :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Integer, Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1003)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_INTEGER,MPI_SUM,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1004)

  End Subroutine gisum_vector

  Subroutine gisum_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - integer scalar version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Intent( InOut ) :: aaa

    Integer                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_INTEGER,MPI_SUM,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine gisum_scalar

  Subroutine grsum_matrix(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - working precision vector
!                                         version
!
! copyright - daresbury laboratory
! author    - i.j.bush december 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Dimension( :, : ), Intent( InOut ) :: aaa

    Integer                                               :: n_l1,n_u1,n_s,fail
    Integer                                               :: n_l2,n_u2
    Real( Kind = wp ), Dimension( :, : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l1 = Lbound(aaa, Dim = 1)
    n_u1 = Ubound(aaa, Dim = 1)

    n_l2 = Lbound(aaa, Dim = 2)
    n_u2 = Ubound(aaa, Dim = 2)

    fail = 0
    Allocate (bbb(n_l1:n_u1,n_l2:n_u2), Stat = fail)
    If (fail > 0) Call error(1048)

    n_s = Size(aaa)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_SUM,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1049)

  End Subroutine grsum_matrix

  Subroutine grsum_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - working precision vector
!                                         version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Dimension( : ), Intent( InOut ) :: aaa

    Integer                                            :: n_l,n_u,n_s,fail
    Real( Kind = wp ), Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1005)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_SUM,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1006)

  End Subroutine grsum_vector

  Subroutine grsum_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global summation subroutine - working precision scalar
!                                         version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Intent( InOut ) :: aaa

    Real( Kind = wp )                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,wp_mpi,MPI_SUM,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine grsum_scalar

  Subroutine gimax_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global maximum subroutine - integer vector version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Dimension( : ), Intent( InOut ) :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Integer, Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1007)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_INTEGER,MPI_MAX,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1008)

  End Subroutine gimax_vector

  Subroutine gimax_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global maximum subroutine - integer scalar version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Intent( InOut ) :: aaa

    Integer                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_INTEGER,MPI_MAX,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine gimax_scalar

  Subroutine grmax_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global maximum subroutine - working precision vector version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Dimension( : ), Intent( InOut ) :: aaa

    Integer                                            :: n_l,n_u,n_s,fail
    Real( Kind = wp ), Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1009)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_MAX,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1010)

  End Subroutine grmax_vector

  Subroutine grmax_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global maximum subroutine - working precision scalar version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Intent( InOut ) :: aaa

    Real( Kind = wp )                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,wp_mpi,MPI_MAX,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine grmax_scalar

  Subroutine gimin_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global minimum subroutine - integer vector version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Dimension( : ), Intent( InOut ) :: aaa

    Integer                                  :: n_l,n_u,n_s,fail
    Integer, Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1044)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,MPI_INTEGER,MPI_MIN,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1045)

  End Subroutine gimin_vector

  Subroutine gimin_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global minimum subroutine - integer scalar version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Integer, Intent( InOut ) :: aaa

    Integer                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,MPI_INTEGER,MPI_MIN,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine gimin_scalar

  Subroutine grmin_vector(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global minimum subroutine - working precision vector version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Dimension( : ), Intent( InOut ) :: aaa

    Integer                                            :: n_l,n_u,n_s,fail
    Real( Kind = wp ), Dimension( : ), Allocatable     :: bbb

    If (mxnode == 1) Return

    n_l = Lbound(aaa, Dim = 1)
    n_u = Ubound(aaa, Dim = 1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat = fail)
    If (fail > 0) Call error(1046)

    n_s = Size(aaa, Dim = 1)

    Call MPI_ALLREDUCE(aaa,bbb,n_s,wp_mpi,MPI_MIN,dlp_comm_world,ierr)

    aaa = bbb

    Deallocate (bbb, Stat = fail)
    If (fail > 0) Call error(1047)

  End Subroutine grmin_vector

  Subroutine grmin_scalar(aaa)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 global minimum subroutine - working precision scalar version
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2007
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Intent( InOut ) :: aaa

    Real( Kind = wp )                  :: bbb

    If (mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa,bbb,1,wp_mpi,MPI_MIN,dlp_comm_world,ierr)

    aaa = bbb

  End Subroutine grmin_scalar

  Subroutine gtime(time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 timing routine
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Real( Kind = wp ), Intent(   Out ) :: time

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

End Module comms_module
