Module comms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for global communication routines and functions
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov february 2015
  ! contrib   - a.m.elena march 2016
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only: dp, &
    qp, &
    si, &
    sp, &
    wi, &
    wp
  Use particle, Only: corePart
  Use, Intrinsic :: iso_fortran_env, Only: CHARACTER_STORAGE_SIZE, error_unit
  Use asserts,  Only : assert
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

  Integer, Public, Parameter  :: root_id = 0

  Integer, Public :: mpi_ver = -1, &
                     mpi_subver = -1

  Character(Len=MPI_MAX_PROCESSOR_NAME), Public :: proc_name = "*"
#ifndef OLDMPI
  Character(Len=MPI_MAX_LIBRARY_VERSION_STRING), Public :: lib_version = "*"
#endif

  ! Message tags

  Integer, Parameter, Public :: Deport_tag = 1100, &
                                Export_tag = 1111, &
                                Revive_tag = 1122, &
                                PassUnit_tag = 1133, &
                                UpdShUnit_tag = 1144, &
                                SysExpand_tag = 1155, &
                                WriteConf_tag = 1166, &
                                Traject_tag = 1177, &
                                Spread_tag = 1188, &
                                DpdVExp_tag = 1199, &
                                MetLdExp_tag = 2200, &
                                ExpMplRM_tag = 2211, &
                                ExchgGrid_tag = 2222, &
                                DefRWrite_tag = 2233, &
                                DefExport_tag = 2244, &
                                DefWrite_tag = 2255, &
                                RsdWrite_tag = 2266, &
                                MsdWrite_tag = 2277, &
                                Grid1_tag = 3300, &
                                Grid2_tag = 3311, &
                                Grid3_tag = 3322, &
                                Grid4_tag = 3333, &
                                Timer_tag = 4000

  ! MPI operations
  Integer, Parameter, Public :: op_sum = MPI_SUM, &
                                op_max = MPI_MAX, &
                                op_min = MPI_MIN, &
                                op_prod = MPI_PROD, &
                                op_land = MPI_LAND, &
                                op_band = MPI_BAND, &
                                op_lor = MPI_LOR, &
                                op_bor = MPI_BOR, &
                                op_lxor = MPI_LXOR, &
                                op_bxor = MPI_BXOR, &
                                op_maxloc = MPI_MAXLOC, &
                                op_minloc = MPI_MINLOC

  Integer, Parameter, Public :: offset_kind = MPI_OFFSET_KIND, &
                                status_size = MPI_STATUS_SIZE, &
                                address_kind = MPI_ADDRESS_KIND, &
                                comm_self = MPI_COMM_SELF, &
                                comm_null = MPI_COMM_NULL

  Integer, Parameter, Public :: mode_wronly = MPI_MODE_WRONLY, &
                                mode_rdonly = MPI_MODE_RDONLY, &
                                mode_create = MPI_MODE_CREATE

  Type, Public :: comms_type
    Integer               :: ierr
    Integer               :: request
    Integer               :: status(1:MPI_STATUS_SIZE) = 0
    Integer               :: comm
    Integer               :: idnode = -1
    Integer               :: mxnode = -1
    Logical               :: l_fast = .false.
    Integer               :: ou
    Integer               :: part_type
    Integer               :: part_array_type
    Integer               :: part_type_positions
    Integer               :: part_array_type_positions
    Integer               :: part_type_forces
    Integer               :: part_array_type_forces
  End Type

  !> @brief Type containing mpi distribution information for a 1D or 2D array
  !!
  !! Holds index arrays required by all/gatherv and all/scatterv mpi routines
  !!
  !! @param  counts        Integer array of length group size.
  !!                       For gatherv: Number of elements that are to be received from each process
  !!                       For scatterv: Number of elements to send to each process
  !! @param  displ         Integer array of length group size.
  !!                       For gatherv: Entry i specifies the displacement relative to recv_buffer
  !!                       at which to place the incoming data from process i
  !!                       For scatterv: Entry i specifies the displacement relative to send_buffer
  !!                       from which to take the outgoing data from root to process i
  !
  Type, Public :: mpi_distribution_type
     Private
     Integer, Public, Allocatable :: counts(:)
     Integer, Public, Allocatable :: displ(:)
   Contains
     Procedure :: init => mpi_distribution_init
     Procedure :: finalise => mpi_distribution_finalise
  End Type mpi_distribution_type

  Public :: init_comms, exit_comms, abort_comms, &
            gsync, gwait, gcheck, gsum, gmax, gtime, gsend, grecv, girecv, &
            gscatter, gscatterv, gscatter_columns, ggatherv, gallgather, galltoall, &
            galltoallv, gallreduce, gallgatherv, gatherv_scatterv_index_arrays, mtime

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

    Module Procedure gcsum_vector
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
    Module Procedure gbcast_integer
    Module Procedure gbcast_integer_scalar
    Module Procedure gbcast_integer_scalar_16
    Module Procedure gbcast_real
    Module Procedure gbcast_real2d
    Module Procedure gbcast_real_scalar
    Module Procedure gbcast_char
    Module Procedure gbcast_char_scalar
    Module Procedure gbcast_logical_scalar
  End Interface !gbcast

  Interface gsend
    Module Procedure gsend_integer_scalar
    Module Procedure gsend_integer_vector
    Module Procedure gsend_real_scalar
    Module Procedure gsend_real_vector
    Module Procedure gsend_real_array3
    Module Procedure gsend_logical_scalar
    Module Procedure gsend_logical_vector
    Module Procedure gsend_character_scalar
    Module Procedure gsend_character_vector
    Module Procedure gsend_particle_scalar
    Module Procedure gsend_particle_vector
  End Interface

  Interface grecv
    Module Procedure grecv_integer_scalar
    Module Procedure grecv_integer_vector
    Module Procedure grecv_real_scalar
    Module Procedure grecv_real_vector
    Module Procedure grecv_real_array3
    Module Procedure grecv_logical_scalar
    Module Procedure grecv_logical_vector
    Module Procedure grecv_character_scalar
    Module Procedure grecv_character_vector
    Module Procedure grecv_particle_scalar
    Module Procedure grecv_particle_vector
  End Interface grecv

  Interface girecv
    Module Procedure girecv_integer_scalar
    Module Procedure girecv_integer_vector
    Module Procedure girecv_real_scalar
    Module Procedure girecv_real_vector
    Module Procedure girecv_real_array3
    Module Procedure girecv_logical_scalar
    Module Procedure girecv_logical_vector
    Module Procedure girecv_character_scalar
    Module Procedure girecv_character_vector
  End Interface girecv

  Interface gscatter
    Module Procedure gscatter_integer_to_scalar
    Module Procedure gscatter_integer_to_vector
    Module Procedure gscatter_real_to_scalar
    Module Procedure gscatter_real_to_vector
  End Interface gscatter

  Interface gscatterv
    Module Procedure gscatterv_integer
    Module Procedure gscatterv_real
    Module Procedure gscatterv_real_2d
    Module Procedure gscatterv_character
  End Interface gscatterv

  Interface gscatter_columns
    Module Procedure gscatter_columns_real
  End Interface gscatter_columns

  Interface ggatherv
     Module Procedure ggatherv_real
     Module Procedure ggatherv_real_2d
  End Interface ggatherv

  Interface gallgather
    Module Procedure gallgather_integer_vector_to_vector
    Module Procedure gallgather_integer_scalar_to_vector
  End Interface gallgather

  Interface gallgatherv
     Module Procedure gallgatherv_character
     Module Procedure gallgatherv_real_1d
     Module Procedure gallgatherv_real_2d
  End Interface gallgatherv

  Interface galltoall
    Module Procedure galltoall_integer
  End Interface galltoall

  Interface galltoallv
    Module Procedure galltoallv_integer
  End Interface galltoallv

  Interface gallreduce
    Module Procedure gallreduce_logical_scalar
    Module Procedure gallreduce_logical_vector
    Module Procedure gallreduce_integer_vector
    Module Procedure gallreduce_real_vector
  End Interface gallreduce

Contains

  Subroutine init_comms(comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for initialising the communication harness
! determining the MPI precision, and node identification and count
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
! contrib   - a.m.elena march 2016
! refactoring:
!           - a.m.elena march-october 2018
!           - j.madge march-october 2018
!           - a.b.g.chalk march-october 2018
!           - i.scivetti march-october 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm

    Integer                        :: lname, lversion
    Integer(KIND=MPI_ADDRESS_KIND) :: base, displacements(1:9), extent, lb
    Integer, Dimension(1:9)        :: block_lengths, types
    Type(corePart)                 :: part_array(1:5), part_temp

#ifdef WITH_DFTBP
    Integer, Parameter :: requiredThreading=MPI_THREAD_FUNNELED
    Integer, Parameter :: errorcode = 0
    Integer            :: providedThreading
    Integer            :: idnode

    Call MPI_INIT_THREAD(requiredThreading, providedThreading, comm%ierr)
    If( providedThreading < requiredThreading )Then
       Call MPI_COMM_RANK(MPI_COMM_WORLD, idnode, comm%ierr)
       if(idnode == root_id)then
          Write(error_unit,'(1x,a,I1,a,I1)') 'error - MPI library threading support, ',&
               providedThreading,', is less than that required by DFTB+ v18.2, ',requiredThreading
          Call MPI_ABORT(MPI_COMM_WORLD, errorcode, comm%ierr)
       End if
    End If
#else
    Call MPI_INIT(comm%ierr)
#endif

    Call MPI_COMM_DUP(MPI_COMM_WORLD, comm%comm, comm%ierr)

    ! use iso_fortran_env
    If (wp == sp) Then
      wp_mpi = MPI_REAL
    Else If (wp == dp) Then
      ! MPI_REAL8 is apparently not in the strict MPI2 standard
      ! It is just an optional data type in the FORTRAN Bindings
      wp_mpi = MPI_DOUBLE_PRECISION
    Else If (wp == qp) Then
      wp_mpi = MPI_REAL16
    Else
      Write (0, '(/,1x,a)') 'error - working precision mismatch between FORTRAN90 and MPI implementation'
      Call abort_comms(comm, 1000)
    End If
    Call MPI_COMM_RANK(comm%comm, comm%idnode, comm%ierr)
    Call MPI_COMM_SIZE(comm%comm, comm%mxnode, comm%ierr)

    !Create the transfer type for the corePart type.
    block_lengths(1:9) = 1
    Call MPI_GET_ADDRESS(part_temp%xxx, displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%yyy, displacements(2), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%zzz, displacements(3), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fxx, displacements(4), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fyy, displacements(5), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fzz, displacements(6), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%chge, displacements(7), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%pad1, displacements(8), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%pad2, displacements(9), comm%ierr)
    base = displacements(1)
    displacements(1:9) = displacements(1:9) - base
    types(1:7) = wp_mpi
    types(8:9) = MPI_INTEGER

    Call MPI_TYPE_CREATE_STRUCT(9, block_lengths, displacements, types, comm%part_type, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_type, comm%ierr)

    Call MPI_GET_ADDRESS(part_array(1), displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_array(2), displacements(2), comm%ierr)
    extent = displacements(2) - displacements(1)
    lb = 0
    Call MPI_TYPE_CREATE_RESIZED(comm%part_type, lb, extent, comm%part_array_type, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_array_type, comm%ierr)

    Call MPI_GET_ADDRESS(part_temp%xxx, displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%yyy, displacements(2), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%zzz, displacements(3), comm%ierr)
    base = displacements(1)
    block_lengths(1:3) = 1
    types(1:3) = wp_mpi
    Call MPI_TYPE_CREATE_STRUCT(3, block_lengths, displacements, types, comm%part_type_positions, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_type_positions, comm%ierr)
    Call MPI_GET_ADDRESS(part_array(1), displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_array(2), displacements(2), comm%ierr)
    extent = displacements(2) - displacements(1)
    lb = 0
    Call MPI_TYPE_CREATE_RESIZED(comm%part_type_positions, lb, extent, comm%part_array_type_positions, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_array_type_positions, comm%ierr)

    Call MPI_GET_ADDRESS(part_temp%fxx, displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fyy, displacements(2), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%fzz, displacements(3), comm%ierr)
    Call MPI_GET_ADDRESS(part_temp%xxx, base, comm%ierr)
    displacements(1:3) = displacements(1:3) - base
    Call MPI_TYPE_CREATE_STRUCT(3, block_lengths, displacements, types, comm%part_type_forces, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_type_forces, comm%ierr)
    Call MPI_GET_ADDRESS(part_array(1), displacements(1), comm%ierr)
    Call MPI_GET_ADDRESS(part_array(2), displacements(2), comm%ierr)
    extent = displacements(2) - displacements(1)
    lb = 0
    Call MPI_TYPE_CREATE_RESIZED(comm%part_type_forces, lb, extent, comm%part_array_type_forces, comm%ierr)
    Call MPI_TYPE_COMMIT(comm%part_array_type_forces, comm%ierr)
#ifndef OLDMPI
    Call MPI_GET_PROCESSOR_NAME(proc_name, lname, comm%ierr)
    Call MPI_GET_VERSION(mpi_ver, mpi_subver, comm%ierr)
    Call MPI_GET_LIBRARY_VERSION(lib_version, lversion, comm%ierr)
#endif

  End Subroutine init_comms


  Subroutine exit_comms(comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for exiting the communication harness
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm(:)

    Integer(Kind=wi) :: i

    Do i = 1, Size(comm, dim=1)
      Call MPI_TYPE_FREE(comm(i)%part_array_type, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_type, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_array_type_positions, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_type_positions, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_array_type_forces, comm(i)%ierr)
      Call MPI_TYPE_FREE(comm(i)%part_type_forces, comm(i)%ierr)
      Call MPI_FINALIZE(comm(i)%ierr)
    End Do

  End Subroutine exit_comms

  Subroutine abort_comms(comm, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for aborting the communication harness
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2008
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: ierr

    Call MPI_ABORT(comm%comm, ierr, comm%ierr)

  End Subroutine abort_comms

  Subroutine gsync(comm, lfast)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 barrier/synchronization routine
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical, Optional               :: lfast

    If (Present(lfast)) comm%l_fast = lfast

    If (comm%mxnode == 1) Return

    Call MPI_BARRIER(comm%comm, comm%ierr)

  End Subroutine gsync

  Subroutine gwait(comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 wait for a non-blocking send or receive to finish
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm

    If (comm%mxnode == 1) Then
      Return
    End If

    Call MPI_WAIT(comm%request, comm%status, comm%ierr)
  End Subroutine gwait

  Subroutine gcheck_vector(comm, aaa, enforce)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - boolean vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2013
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type),           Intent(InOut) :: comm
    Logical, Dimension(:),      Intent(InOut) :: aaa
    Character(Len=*), Optional, Intent(In   ) :: enforce

    Integer                            :: fail, n_l, n_s, n_u
    Logical, Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Then
      Return
    Else
      If (.not. Present(enforce)) Then
        If (comm%l_fast) Return
      End If
    End If

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms_module -> gcheck_vector'
      Call abort_comms(comm, 1001)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, MPI_LOGICAL, MPI_LAND, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms_module -> gcheck_vector'
      Call abort_comms(comm, 1002)
    End If

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

    Type(comms_type),           Intent(InOut) :: comm
    Logical,                    Intent(InOut) :: aaa
    Character(Len=*), Optional, Intent(In   ) :: enforce

    Logical :: bbb

    If (comm%mxnode == 1) Then
      Return
    Else
      If (.not. Present(enforce)) Then
        If (comm%l_fast) Return
      End If
    End If

    Call MPI_ALLREDUCE(aaa, bbb, 1, MPI_LOGICAL, MPI_LAND, comm%comm, comm%ierr)

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

    Type(comms_type),      Intent(InOut) :: comm
    Integer, Dimension(:), Intent(InOut) :: aaa

    Integer                            :: fail, n_l, n_s, n_u
    Integer, Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> gisum_vector'
      Call abort_comms(comm, 1003)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, MPI_INTEGER, MPI_SUM, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)

    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> gisum_vector'
      Call abort_comms(comm, 1004)
    End If

  End Subroutine gisum_vector

  Subroutine gisum_scalar(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - integer scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: aaa

    Integer :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa, bbb, 1, MPI_INTEGER, MPI_SUM, comm%comm, comm%ierr)

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

    Type(comms_type),               Intent(InOut) :: comm
    Real(Kind=wp), Dimension(:, :), Intent(InOut) :: aaa

    Integer                                     :: fail, n_l1, n_l2, n_s, n_u1, n_u2
    Real(Kind=wp), Allocatable, Dimension(:, :) :: bbb

    If (comm%mxnode == 1) Return

    n_l1 = Lbound(aaa, Dim=1)
    n_u1 = Ubound(aaa, Dim=1)

    n_l2 = Lbound(aaa, Dim=2)
    n_u2 = Ubound(aaa, Dim=2)

    fail = 0
    Allocate (bbb(n_l1:n_u1, n_l2:n_u2), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> grsum_matrix'
      Call abort_comms(comm, 1048)
    End If

    n_s = Size(aaa)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, wp_mpi, MPI_SUM, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> grsum_matrix'
      Call abort_comms(comm, 1049)
    End If

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

    Type(comms_type),            Intent(InOut) :: comm
    Real(Kind=wp), Dimension(:), Intent(InOut) :: aaa

    Integer                                  :: fail, n_l, n_s, n_u
    Real(Kind=wp), Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> grsum_vector'
      Call abort_comms(comm, 1005)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, wp_mpi, MPI_SUM, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> grsum_vector'
      Call abort_comms(comm, 1006)
    End If

  End Subroutine grsum_vector

  Subroutine gcsum_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - complex vector
    !                                         version
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena december 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type),               Intent(InOut) :: comm
    Complex(Kind=wp), Dimension(:), Intent(InOut) :: aaa

    Complex(Kind=wp), Allocatable, Dimension(:) :: bbb
    Integer                                     :: fail, n_l, n_s, n_u

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> gcsum_vector'
      Call abort_comms(comm, 1005)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, MPI_DOUBLE_COMPLEX, MPI_SUM, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> gcsum_vector'
      Call abort_comms(comm, 1006)
    End If

  End Subroutine gcsum_vector

  Subroutine grsum_scalar(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global summation subroutine - working precision scalar
    !                                         version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: aaa

    Real(Kind=wp) :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa, bbb, 1, wp_mpi, MPI_SUM, comm%comm, comm%ierr)

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

    Type(comms_type),      Intent(InOut) :: comm
    Integer, Dimension(:), Intent(InOut) :: aaa

    Integer                            :: fail, n_l, n_s, n_u
    Integer, Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> gimax_vector'
      Call abort_comms(comm, 1007)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, MPI_INTEGER, MPI_MAX, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> gimax_vector'
      Call abort_comms(comm, 1008)
    End If

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

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: aaa

    Integer :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa, bbb, 1, MPI_INTEGER, MPI_MAX, comm%comm, comm%ierr)

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

    Type(comms_type),            Intent(InOut) :: comm
    Real(Kind=wp), Dimension(:), Intent(InOut) :: aaa

    Integer                                  :: fail, n_l, n_s, n_u
    Real(Kind=wp), Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> grmax_vector'
      Call abort_comms(comm, 1009)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, wp_mpi, MPI_MAX, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> grmax_vector'
      Call abort_comms(comm, 1010)
    End If

  End Subroutine grmax_vector

  Subroutine grmax_scalar(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global maximum subroutine - working precision scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: aaa

    Real(Kind=wp) :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa, bbb, 1, wp_mpi, MPI_MAX, comm%comm, comm%ierr)

    aaa = bbb

  End Subroutine grmax_scalar

  Subroutine gimin_vector(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global minimum subroutine - integer vector version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2007
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type),      Intent(InOut) :: comm
    Integer, Dimension(:), Intent(InOut) :: aaa

    Integer                            :: fail, n_l, n_s, n_u
    Integer, Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> gimin_vector'
      Call abort_comms(comm, 1044)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, MPI_INTEGER, MPI_MIN, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> gimin_vector'
      Call abort_comms(comm, 1045)
    End If

  End Subroutine gimin_vector

  Subroutine gimin_scalar(comm, aaa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 global minimum subroutine - integer scalar version
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2007
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: aaa

    Integer :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa, bbb, 1, MPI_INTEGER, MPI_MIN, comm%comm, comm%ierr)

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

    Type(comms_type),            Intent(InOut) :: comm
    Real(Kind=wp), Dimension(:), Intent(InOut) :: aaa

    Integer                                  :: fail, n_l, n_s, n_u
    Real(Kind=wp), Allocatable, Dimension(:) :: bbb

    If (comm%mxnode == 1) Return

    n_l = Lbound(aaa, Dim=1)
    n_u = Ubound(aaa, Dim=1)

    fail = 0
    Allocate (bbb(n_l:n_u), Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - allocation failure in comms -> grmin_vector'
      Call abort_comms(comm, 1046)
    End If

    n_s = Size(aaa, Dim=1)

    Call MPI_ALLREDUCE(aaa, bbb, n_s, wp_mpi, MPI_MIN, comm%comm, comm%ierr)

    aaa = bbb

    Deallocate (bbb, Stat=fail)
    If (fail > 0) Then
      Write (comm%ou, '(/,1x,a)') 'error - deallocation failure in comms -> grmin_vector'
      Call abort_comms(comm, 1047)
    End If

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

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: aaa

    Real(Kind=wp) :: bbb

    If (comm%mxnode == 1) Return

    Call MPI_ALLREDUCE(aaa, bbb, 1, wp_mpi, MPI_MIN, comm%comm, comm%ierr)

    aaa = bbb

  End Subroutine grmin_scalar

  Subroutine gbcast_integer(comm, vec, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast an integer array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena may 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: root

    Integer :: n_l, n_s, n_u

    If (comm%mxnode == 1) Return
    n_l = Lbound(vec, Dim=1)
    n_u = Ubound(vec, Dim=1)
    n_s = Size(vec, Dim=1)

    Call MPI_BCAST(vec(n_l:n_u), n_s, MPI_INTEGER, root, comm%comm, comm%ierr)

  End Subroutine gbcast_integer

  Subroutine gbcast_logical_scalar(comm, s, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast an integer array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena february 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(InOut) :: s
    Integer,          Intent(In   ) :: root

    If (comm%mxnode == 1) Return

    Call MPI_BCAST(s, 1, MPI_LOGICAL, root, comm%comm, comm%ierr)

  End Subroutine gbcast_logical_scalar


  Subroutine gbcast_integer_scalar(comm, s, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast an integer array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena may 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: s
    Integer,          Intent(In   ) :: root

    If (comm%mxnode == 1) Return

    Call MPI_BCAST(s, 1, MPI_INTEGER, root, comm%comm, comm%ierr)

  End Subroutine gbcast_integer_scalar

  Subroutine gbcast_integer_scalar_16(comm, s, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast an integer array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena may 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Integer(Kind=si), Intent(InOut) :: s
    Integer,          Intent(In   ) :: root

    If (comm%mxnode == 1) Return

    Call MPI_BCAST(s, 1, MPI_INTEGER, root, comm%comm, comm%ierr)

  End Subroutine gbcast_integer_scalar_16

  Subroutine gbcast_real(comm, vec, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast a real array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena march 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: root

    Integer :: n_l, n_s, n_u

    If (comm%mxnode == 1) Return
    n_l = Lbound(vec, Dim=1)
    n_u = Ubound(vec, Dim=1)
    n_s = Size(vec, Dim=1)

    Call MPI_BCAST(vec(n_l:n_u), n_s, wp_mpi, root, comm%comm, comm%ierr)
  End Subroutine gbcast_real

   Subroutine gbcast_real2d(comm, array, root)
    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: array(:,:)
    Integer,          Intent(In   ) :: root
    If (comm%mxnode == 1) Return
    Call MPI_BCAST(array, size(array), wp_mpi, root, comm%comm, comm%ierr)
  End Subroutine gbcast_real2d

  Subroutine gbcast_real_scalar(comm, s, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast a real array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena march 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: s
    Integer,          Intent(In   ) :: root

    If (comm%mxnode == 1) Return

    Call MPI_BCAST(s, 1, wp_mpi, root, comm%comm, comm%ierr)

  End Subroutine gbcast_real_scalar

  Subroutine gbcast_char_scalar(comm, vec, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast a real array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena Nov 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(InOut) :: vec
    Integer,          Intent(In   ) :: root

    Integer :: n_s

    If (comm%mxnode == 1) Return
    n_s = Len(vec) * CHARACTER_STORAGE_SIZE / 8

    Call MPI_BCAST(vec, n_s, MPI_CHARACTER, root, comm%comm, comm%ierr)
  End Subroutine gbcast_char_scalar

  Subroutine gbcast_char(comm, vec, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 broadcast a real array subroutine
    !
    ! copyright - daresbury laboratory
    ! author    - a.m.elena march 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: root

    Integer :: n_l, n_s, n_u

    If (comm%mxnode == 1) Return
    n_l = Lbound(vec, Dim=1)
    n_u = Ubound(vec, Dim=1)
    n_s = Size(vec, Dim=1) * Len(vec(n_l)) * CHARACTER_STORAGE_SIZE / 8

    Call MPI_BCAST(vec(n_l:n_u), n_s, MPI_CHARACTER, root, comm%comm, comm%ierr)
  End Subroutine gbcast_char

  Subroutine mtime(time)
    Real(Kind=wp), Intent(  Out) :: time

    time = MPI_WTIME()
  End Subroutine mtime

  Subroutine gtime(time)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 timing routine
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov may 2004
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(  Out) :: time

    Logical, Save       :: newjob = .true.
    Real(Kind=wp), Save :: tzero

    If (newjob) Then
      newjob = .false.

      tzero = MPI_WTIME()
      time = 0.0_wp
    Else
      time = MPI_WTIME() - tzero
    End If

  End Subroutine gtime

  Subroutine gsend_integer_scalar(comm, s, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send an integer
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: s, dest, tag

    Call MPI_SEND(s, 1, MPI_INTEGER, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_integer_scalar

  Subroutine gsend_integer_vector(comm, vec, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send an integer array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: vec(:), dest, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_SEND(vec(:), n_s, MPI_INTEGER, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_integer_vector

  Subroutine gsend_real_scalar(comm, s, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a real scalar
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: s
    Integer,          Intent(In   ) :: dest, tag

    Call MPI_SEND(s, 1, wp_mpi, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_real_scalar

  Subroutine gsend_real_vector(comm, vec, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a real array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: vec(:)
    Integer,          Intent(In   ) :: dest, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_SEND(vec(:), n_s, wp_mpi, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_real_vector

  Subroutine gsend_real_array3(comm, arr, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a real three dimensional array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: arr(:, :, :)
    Integer,          Intent(In   ) :: dest, tag

    Integer               :: i
    Integer, Dimension(3) :: n_s

    Do i = 1, 3
      n_s(i) = Size(arr, Dim=i)
    End Do

    Call MPI_SEND(arr(:, :, :), Product(n_s(1:3)), &
                  wp_mpi, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_real_array3

  Subroutine gsend_logical_scalar(comm, s, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a logical scalar
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(In   ) :: s
    Integer,          Intent(In   ) :: dest, tag

    Call MPI_SEND(s, 1, MPI_LOGICAL, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_logical_scalar

  Subroutine gsend_logical_vector(comm, vec, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a logical array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(In   ) :: vec(:)
    Integer,          Intent(In   ) :: dest, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_SEND(vec(:), n_s, MPI_LOGICAL, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_logical_vector

  Subroutine gsend_particle_scalar(comm, s, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a particle scalar
    !
    ! copyright - daresbury laboratory
    ! author    - a.chalk june 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Type(corePart),   Intent(In   ) :: s
    Integer,          Intent(In   ) :: dest, tag

    Call MPI_Send(s, 1, comm%part_type, dest, tag, comm%comm, comm%ierr)

  End Subroutine

  Subroutine gsend_particle_vector(comm, vec, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a particle vector
    !
    ! copyright - daresbury laboratory
    ! author    - a.chalk june 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Type(corePart),   Intent(In   ) :: vec(:)
    Integer,          Intent(In   ) :: dest, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_SEND(vec(:), n_s, comm%part_array_type, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_particle_vector

  Subroutine gsend_character_scalar(comm, s, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a character string
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(In   ) :: s
    Integer,          Intent(In   ) :: dest, tag

    Integer :: n_s

    n_s = Len(s)

    Call MPI_SEND(s, n_s, MPI_CHARACTER, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_character_scalar

  Subroutine gsend_character_vector(comm, vec, dest, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 send a character string array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(In   ) :: vec(:)
    Integer,          Intent(In   ) :: dest, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1) * Len(vec(1))

    Call MPI_SEND(vec(:), n_s, MPI_CHARACTER, dest, tag, comm%comm, comm%ierr)
  End Subroutine gsend_character_vector

  Subroutine grecv_integer_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive an integer scalar
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(  Out) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_RECV(s, 1, MPI_INTEGER, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_integer_scalar

  Subroutine grecv_integer_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive an integer vector
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_RECV(vec(:), n_s, MPI_INTEGER, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_integer_vector

  Subroutine grecv_real_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a real scalar
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(  Out) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_RECV(s, 1, wp_mpi, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_real_scalar

  Subroutine grecv_real_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a real vector
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_RECV(vec(:), n_s, wp_mpi, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_real_vector

  Subroutine grecv_real_array3(comm, arr, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a real three dimensional array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: arr(:, :, :)
    Integer,          Intent(In   ) :: source, tag

    Integer               :: i
    Integer, Dimension(3) :: n_s

    Do i = 1, 3
      n_s(i) = Size(arr, Dim=i)
    End Do

    Call MPI_RECV(arr(:, :, :), Product(n_s(1:3)), &
                  wp_mpi, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_real_array3

  Subroutine grecv_logical_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a logical scalar
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(  Out) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_RECV(s, 1, MPI_LOGICAL, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_logical_scalar

  Subroutine grecv_logical_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a logical vector
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_RECV(vec(:), n_s, MPI_LOGICAL, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_logical_vector

  Subroutine grecv_character_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a character string
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(InOut) :: s
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Len(s)

    Call MPI_RECV(s, n_s, MPI_CHARACTER, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_character_scalar

  Subroutine grecv_character_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a character string array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1) * Len(vec(1))

    Call MPI_RECV(vec(:), n_s, MPI_CHARACTER, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_character_vector
  Subroutine grecv_particle_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a corePart
    !
    ! copyright - daresbury laboratory
    ! author    - a.chalk july 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Type(corePart),   Intent(InOut) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_RECV(s, 1, comm%part_type, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_particle_scalar
  Subroutine grecv_particle_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a corePart array
    !
    ! copyright - daresbury laboratory
    ! author    - a.chalk july 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type), Intent(InOut) :: comm
    Type(corePart),   Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)
    Call MPI_RECV(vec(:), n_s, comm%part_array_type, source, tag, comm%comm, comm%status, comm%ierr)
  End Subroutine grecv_particle_vector
  Subroutine girecv_integer_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive an integer scalar (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(  Out) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_IRECV(s, 1, MPI_INTEGER, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_integer_scalar

  Subroutine girecv_integer_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive an integer vector (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_IRECV(vec(:), n_s, MPI_INTEGER, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_integer_vector

  Subroutine girecv_real_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a real scalar (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(  Out) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_IRECV(s, 1, wp_mpi, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_real_scalar

  Subroutine girecv_real_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a real vector (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_IRECV(vec(:), n_s, wp_mpi, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_real_vector

  Subroutine girecv_real_array3(comm, arr, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a real three dimensional array (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(InOut) :: arr(:, :, :)
    Integer,          Intent(In   ) :: source, tag

    Integer               :: i
    Integer, Dimension(3) :: n_s

    Do i = 1, 3
      n_s(i) = Size(arr, Dim=i)
    End Do

    Call MPI_IRECV(arr(:, :, :), Product(n_s(1:3)), &
                   wp_mpi, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_real_array3

  Subroutine girecv_logical_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a logical scalar (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(  Out) :: s
    Integer,          Intent(In   ) :: source, tag

    Call MPI_IRECV(s, 1, MPI_LOGICAL, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_logical_scalar

  Subroutine girecv_logical_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a logical vector (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1)

    Call MPI_IRECV(vec(:), n_s, MPI_LOGICAL, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_logical_vector

  Subroutine girecv_character_scalar(comm, s, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a character string (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(InOut) :: s
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Len(s)

    Call MPI_IRECV(s, n_s, MPI_CHARACTER, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_character_scalar

  Subroutine girecv_character_vector(comm, vec, source, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 receive a character string array (non-blocking)
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Character(Len=*), Intent(InOut) :: vec(:)
    Integer,          Intent(In   ) :: source, tag

    Integer :: n_s

    n_s = Size(vec, Dim=1) * Len(vec(1))

    Call MPI_IRECV(vec(:), n_s, MPI_CHARACTER, source, tag, comm%comm, comm%request, comm%ierr)
  End Subroutine girecv_character_vector

  Subroutine gscatter_integer_to_scalar(comm, sendbuf, recv, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter an integer buffer to scalar variables
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: sendbuf(:)
    Integer,          Intent(  Out) :: recv
    Integer,          Intent(In   ) :: root

    Call MPI_SCATTER(sendbuf(:), 1, MPI_INTEGER, &
                     recv, 1, MPI_INTEGER, root, comm%comm, comm%ierr)
  End Subroutine gscatter_integer_to_scalar

  Subroutine gscatter_integer_to_vector(comm, sendbuf, scount, recvbuf, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter an integer buffer to vector variables
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: sendbuf(:), scount
    Integer,          Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: root

    Integer :: r_s

    r_s = Size(recvbuf, Dim=1)

    Call MPI_SCATTER(sendbuf(:), scount, MPI_INTEGER, &
                     recvbuf(:), r_s, MPI_INTEGER, root, comm%comm, comm%ierr)
  End Subroutine gscatter_integer_to_vector

  Subroutine gscatter_real_to_scalar(comm, sendbuf, recv, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter a real buffer to scalar variables
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: sendbuf(:)
    Real(Kind=wp),    Intent(  Out) :: recv
    Integer,          Intent(In   ) :: root

    Call MPI_SCATTER(sendbuf(:), 1, wp_mpi, &
                     recv, 1, wp_mpi, root, comm%comm, comm%ierr)
  End Subroutine gscatter_real_to_scalar

  Subroutine gscatter_real_to_vector(comm, sendbuf, scount, recvbuf, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter a real buffer to vector variables
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: sendbuf(:)
    Integer,          Intent(In   ) :: scount
    Real(Kind=wp),    Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: root

    Integer :: r_s

    r_s = Size(recvbuf, Dim=1)

    Call MPI_SCATTER(sendbuf(:), scount, wp_mpi, &
                     recvbuf(:), r_s, wp_mpi, root, comm%comm, comm%ierr)
  End Subroutine gscatter_real_to_vector

  Subroutine gscatterv_integer(comm, sendbuf, scounts, disps, recvbuf, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter an integer buffer with a defined displacement per
    ! node
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: sendbuf(:), scounts(:), disps(:)
    Integer,          Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: root

    Integer :: r_s

    r_s = Size(recvbuf, Dim=1)

    Call MPI_SCATTERV(sendbuf(:), scounts(:), disps(:), MPI_INTEGER, &
                      recvbuf(:), r_s, MPI_INTEGER, &
                      root, comm%comm, comm%ierr)
  End Subroutine gscatterv_integer

  Subroutine gscatterv_real(comm, sendbuf, scounts, disps, recvbuf, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter a real buffer with a defined displacement per
    ! node
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: sendbuf(:)
    Integer,          Intent(In   ) :: scounts(:), disps(:)
    Real(Kind=wp),    Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: root

    Integer :: r_s

    r_s = Size(recvbuf, Dim=1)

    Call MPI_SCATTERV(sendbuf(:), scounts(:), disps(:), wp_mpi, &
                      recvbuf(:), r_s, wp_mpi, &
                      root, comm%comm, comm%ierr)
  End Subroutine gscatterv_real

  Subroutine gscatterv_real_2d(comm, send_buf, send_counts, disps, recv_buf, root)
    !> @brief Scatter the columns of a real two dimensional array
    !! Equivalent operation to gscatter_columns, but send_counts and disps
    !! consistent with a single (composite), contiguous index.
    !!
    !! dl_poly_4 Scatter columns from send buffer on root, to each task.
    !! Each task can receive a unique number of columns.
    !! This implimentation relies on arrays being column major as defined in the
    !! Fortran standard.
    !!
    !! send_counts and disps can be generated using 'gatherv_scatterv_index_arrays'
    !!
    !! @param[in]     comm             Object containing MPI settings and comms
    !! @param[in]     send_buf         Send buffer array, defined on root
    !! @param[inout]  send_counts      Holds size of send buffer for each process
    !! @param[inout]  disps            Holds start index of send buffer for
    !!                                 receive buffer (single, contiguous index)
    !! @param[out]    recv_buf         Local receive buffer array. Each receive
    !!                                 buffer must have the same number of rows,
    !!                                 but can have a unique number of columns
    !!
    !! copyright - daresbury laboratory
    !! author    - A. Buccheri June 2019
    !

    Type( comms_type ), Intent( InOut ) :: comm
    Real( Kind = wp ),  Intent( In    ) :: send_buf(:,:)
    Integer,            Intent( In    ) :: send_counts(:)
    Integer,            Intent( In    ) :: disps(:)
    Real( Kind = wp),   Intent( Out   ) :: recv_buf(:,:)
    Integer,            Intent( In    ) :: root
    Character( Len = 120)               :: error_message

    error_message = 'gscatterv: Sum of buffer sizes in send_counts '//  &
         'does not equal total size of send_buffer.'
    Call assert(Sum(send_counts) == Size(send_buf), error_message)

    error_message = 'gscatterv: Number of rows differs between send and '// &
         'receive buffers. Data will not be correctly ordered.'
    Call assert(Size(send_buf, Dim = 1) == Size(recv_buf, Dim = 1), error_message)

    Call MPI_SCATTERV(send_buf, send_counts, disps, wp_mpi, &
                      recv_buf, Size(recv_buf),     wp_mpi, &
                      root,     comm%comm, comm%ierr)

  End Subroutine gscatterv_real_2d

  Subroutine gscatterv_character(comm, send_buf, send_counts, disps, recv_buf, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter a real buffer with a defined displacement per
    ! node
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018, A.Buccheri June 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(comms_type),     Intent( InOut ) :: comm
    Character( Len = * ), Intent( In    ) :: send_buf(:)
    Integer,              Intent( In    ) :: send_counts(:)
    Integer,              Intent( In    ) :: disps(:)
    Character( Len = * ), Intent(   Out ) :: recv_buf(:)
    Integer,              Intent( In    ) :: root

    Integer               :: send_len, recv_len, send_size, recv_size
    Character( Len = 100 ):: error_message

    send_len  = Len(send_buf(1))
    send_size = Size(send_buf, Dim = 1)
    recv_len  = Len(recv_buf(1))
    recv_size = Size(recv_buf, Dim = 1)

    error_message = 'gscatterv: Character length of send and receive buffers differ.'
    Call assert(send_len == recv_len, error_message)

    error_message = 'gscatterv: Sum of counts to send exceeds the size of send buffer.'
    Call assert(Sum(send_counts) <= send_size, error_message)

    Call MPI_SCATTERV(send_buf, send_counts*send_len, disps*send_len, MPI_CHARACTER, &
                      recv_buf, recv_len*recv_size,                   MPI_CHARACTER, &
                      root, comm%comm, comm%ierr)

  End Subroutine gscatterv_character


  Subroutine gscatter_columns_real(comm, sendbuf, scounts, disps, recvbuf, root)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 scatter the columns of a real two dimensional array
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Real(Kind=wp),    Intent(In   ) :: sendbuf(:, :)
    Integer,          Intent(In   ) :: scounts(:), disps(:)
    Real(Kind=wp),    Intent(  Out) :: recvbuf(:, :)
    Integer,          Intent(In   ) :: root

    Integer :: r_c, r_s, s_c

    s_c = Size(sendbuf, Dim=1)
    r_c = Size(recvbuf, Dim=1)
    r_s = Size(recvbuf, Dim=2)

    ! This implimentation relies on arrays being column major as defined in the
    ! Fortran standard
    Call MPI_SCATTERV(sendbuf(:, :), &
                      scounts(:) * s_c, disps(:) * s_c, wp_mpi, &
                      recvbuf(:, :), &
                      r_s * r_c, wp_mpi, &
                      root, comm%comm, comm%ierr)
  End Subroutine gscatter_columns_real

  Subroutine ggatherv_real(comm, send_buf, recv_counts, disps, recv_buf, recv_process)
    !> @brief Gather send buffer from each task and broadcast the combined data
    !! to single process
    !!
    !! dl_poly_4 Gather data of type real(wp), stored in send buffers of
    !! variable size into a receive buffer, and broadcast to specified process
    !! If no process is specified, data is gathered on root_id
    !!
    !! @param[in]     comm             Object containing MPI settings and comms
    !! @param[in]     send_buf         Local send buffer array
    !! @param[inout]  recv_counts      Holds size of send buffer of each process
    !! @param[inout]  disps            Holds offset of send buffer w.r.t. index in
    !!                                 receive buffer (single, contiguous index)
    !! @param[out]    recv_buf         Global receive buffer array
    !! @param[in]     recv_process     Optional. Process id for received data
    !!
    !! copyright - daresbury laboratory
    !! author    - A. Buccheri Jan 2020
    !

    Type( comms_type ), Intent( InOut ) :: comm
    Real( Kind = wp ),  Intent( In    ) :: send_buf(:)
    Integer,            Intent( In    ) :: recv_counts(:)
    Integer,            Intent( In    ) :: disps(:)
    Real( Kind = wp),   Intent( Out   ) :: recv_buf(:)
    Integer, Optional,  Intent( In    ) :: recv_process

    Integer :: process_id
    Character( Len = 100 )              :: error_message

    error_message = 'ggatherv: Sum of send_buffer sizes in receive_counts '//&
         'does not equal total size of receive_buffer.'
    Call assert(Sum(recv_counts) == Size(recv_buf), error_message)

    If(Present(recv_process)) Then
       process_id = recv_process
    Else
       process_id = root_id
    Endif

    Call MPI_GATHERV(send_buf,   Size(send_buf),     wp_mpi, &
                     recv_buf,   recv_counts, disps, wp_mpi, &
                     process_id, comm%comm, comm%ierr)

  End Subroutine ggatherv_real


  Subroutine ggatherv_real_2d(comm, send_buf, recv_counts, disps, recv_buf, recv_process)
    !> @brief Gather send buffer from each task and broadcast the combined data
    !! to single process
    !!
    !! dl_poly_4 Gather data of type real, stored in 2D send buffers of size(n,m),
    !! into a receive buffer and broadcasts to a specified process.
    !! If a process is not specified, master is chosen.
    !! Exploits column-major storage of fortran arrays. For ordering to be
    !! correct in recv_buf, n must be fixed for all send_buf, but m can vary.
    !!
    !! recv_counts and disps can be generated using 'gatherv_scatterv_index_arrays'
    !!
    !! @param[in]     comm             Object containing MPI settings and comms
    !! @param[in]     send_buf         Local send buffer array. Each send buffer must
    !!                                 have the same number of rows
    !! @param[inout]  recv_counts      Holds size of send buffer of each process
    !! @param[inout]  disps            Holds offset of send buffer w.r.t. index in
    !!                                 receive buffer (single, contiguous index)
    !! @param[out]    recv_buf         Global receive buffer array
    !! @param[in]     recv_process     Optional. Process id for received data
    !!
    !! copyright - daresbury laboratory
    !! author    - A. Buccheri Jan 2020
    !

    Type( comms_type ), Intent( InOut ) :: comm
    Real( Kind = wp ),  Intent( In    ) :: send_buf(:,:)
    Integer,            Intent( In    ) :: recv_counts(:)
    Integer,            Intent( In    ) :: disps(:)
    Real( Kind = wp),   Intent( Out   ) :: recv_buf(:,:)
    Integer, Optional,  Intent( In    ) :: recv_process

    Integer :: process_id
    Character( Len = 120 ) :: error_message

    error_message = 'ggatherv: Sum of send_buffer sizes in receive_counts '//&
         'does not equal total size of receive_buffer.'
    Call assert(Sum(recv_counts) == Size(recv_buf), error_message)

    error_message = 'ggatherv: Number of rows differs between send and '//&
         'receive buffers. Data will not be correctly ordered.'
    Call assert(Size(send_buf, Dim = 1) == Size(recv_buf, Dim = 1), error_message)

    If(Present(recv_process)) Then
       process_id = recv_process
    Else
       process_id = root_id
    Endif

    Call MPI_GATHERV(send_buf,   Size(send_buf),     wp_mpi, &
                     recv_buf,   recv_counts, disps, wp_mpi, &
                     process_id, comm%comm, comm%ierr)

  End Subroutine ggatherv_real_2d


  Subroutine gallgather_integer_vector_to_vector(comm, sendbuf, recvbuf, rcount)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 gather rcount integers from all processes then distributes them
    ! to all processes. Each process therefore receives and identical and
    ! complete buffer.
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: sendbuf(:)
    Integer,          Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: rcount

    Integer :: s_s

    Call MPI_ALLGATHER(sendbuf(:), s_s, MPI_INTEGER, &
                       recvbuf(:), rcount, MPI_INTEGER, &
                       comm%comm, comm%ierr)
  End Subroutine gallgather_integer_vector_to_vector

  Subroutine gallgather_integer_scalar_to_vector(comm, s, recvbuf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 gather rcount integers from all processes then distributes them
    ! to all processes. Each process therefore receives and identical and
    ! complete buffer.
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: s
    Integer,          Intent(  Out) :: recvbuf(:)

    Call MPI_ALLGATHER(s, 1, MPI_INTEGER, &
                       recvbuf(:), 1, MPI_INTEGER, &
                       comm%comm, comm%ierr)
  End Subroutine gallgather_integer_scalar_to_vector

 Subroutine gallgatherv_character(comm, send_buf, recv_counts, disps, recv_buf)
    !> @brief Gather data from each task and broadcast the combined data to all tasks.
    !!
    !! dl_poly_4 Gather data of type character, stored in send buffers of variable
    !! size (fixed length), into a receive buffer, and broadcast to all tasks.
    !!
    !! @param[in]     comm             Object containing MPI settings and comms
    !! @param[in]     send_buf         Local send buffer array.
    !!                                 Fixed character length but variable size.
    !! @param[inout]  recv_counts      Holds size of send buffer of each process
    !! @param[inout]  disps            Holds offset of send buffer w.r.t. index in
    !!                                 receive buffer (single, contiguous index)
    !! @param[out]    recv_buf         Global receive buffer array
    !!
    !! copyright - daresbury laboratory
    !! author    - A. Buccheri June 2019
    !
    Type( comms_type ),   Intent( InOut ) :: comm
    Character( Len = * ), Intent( In    ), Contiguous :: send_buf(:)
    Integer,              Intent( In    ) :: recv_counts(:)
    Integer,              Intent( In    ) :: disps(:)
    Character( Len = * ), Intent( Out   ) :: recv_buf(:)

    Integer               :: send_len, recv_len, send_size, recv_size
    Character( Len = 100) :: error_message

    send_len  = Len(send_buf(1))
    send_size = Size(send_buf, Dim = 1)
    recv_len  = Len(recv_buf(1))
    recv_size = Size(recv_buf, Dim = 1)

    error_message = 'gallgatherv: Character length of send and receive buffers differ.'
    Call assert(send_len == recv_len, error_message)

    error_message = 'gallgatherv: Sum of counts to receive exceeds the size of receive buffer.'
    Call assert(Sum(recv_counts) <= recv_size, error_message)

    Call MPI_ALLGATHERV(send_buf,  send_size*send_len,                   MPI_CHARACTER, &
                        recv_buf,  recv_counts*recv_len, disps*recv_len, MPI_CHARACTER, &
                        comm%comm, comm%ierr )

  End Subroutine gallgatherv_character

    Subroutine gallgatherv_real_1d(comm, send_buf, recv_counts, disps, recv_buf)
    !> @brief Gather send buffer from each task and broadcast the combined data
    !! to all tasks.
    !!
    !! dl_poly_4 Gather data of type real(wp), stored in send buffers of
    !! variable size into a receive buffer, and broadcast to all tasks.
    !!
    !! @param[in]     comm             Object containing MPI settings and comms
    !! @param[in]     send_buf         Local send buffer array
    !! @param[inout]  recv_counts      Holds size of send buffer of each process
    !! @param[inout]  disps            Holds offset of send buffer w.r.t. index in
    !!                                 receive buffer (single, contiguous index)
    !! @param[out]    recv_buf         Global receive buffer array
    !!
    !! copyright - daresbury laboratory
    !! author    - A. Buccheri June 2019
    !

    Type( comms_type ), Intent( InOut ) :: comm
    Real( Kind = wp ),  Intent( In    ) :: send_buf(:)
    Integer,            Intent( In    ) :: recv_counts(:)
    Integer,            Intent( In    ) :: disps(:)
    Real( Kind = wp),   Intent( Out   ) :: recv_buf(:)
    Character( Len = 100 )              :: error_message

    error_message = 'gallgatherv: Sum of send_buffer sizes in receive_counts '//&
         'does not equal total size of receive_buffer.'
    Call assert(Sum(recv_counts) == Size(recv_buf), error_message)

    Call MPI_ALLGATHERV(send_buf,  Size(send_buf),     wp_mpi, &
                        recv_buf,  recv_counts, disps, wp_mpi, &
                        comm%comm, comm%ierr )

  End Subroutine gallgatherv_real_1d


  Subroutine gallgatherv_real_2d(comm, send_buf, recv_counts, disps, recv_buf)
    !> @brief Gather send buffer from each task and broadcast the combined data
    !! to all tasks.
    !!
    !! dl_poly_4 Gather data of type real, stored in 2D send buffers of size(n,m),
    !! into a receive buffer and broadcast to all tasks.
    !! Exploits column-major storage of fortran arrays. For ordering to be
    !! correct in recv_buf, n must be fixed for all send_buf, but m can vary.
    !!
    !! recv_counts and disps can be generated using 'gatherv_scatterv_index_arrays'
    !!
    !!
    !! @param[in]     comm             Object containing MPI settings and comms
    !! @param[in]     send_buf         Local send buffer array. Each send buffer must
    !!                                 have the same number of rows
    !! @param[inout]  recv_counts      Holds size of send buffer of each process
    !! @param[inout]  disps            Holds offset of send buffer w.r.t. index in
    !!                                 receive buffer (single, contiguous index)
    !! @param[out]    recv_buf         Global receive buffer array
    !!
    !! copyright - daresbury laboratory
    !! author    - A. Buccheri June 2019
    !

    Type( comms_type ), Intent( InOut ) :: comm
    Real( Kind = wp ),  Intent( In    ) :: send_buf(:,:)
    Integer,            Intent( In    ) :: recv_counts(:)
    Integer,            Intent( In    ) :: disps(:)
    Real( Kind = wp),   Intent( Out   ) :: recv_buf(:,:)
    Character( Len = 150 )              :: error_message

    error_message = 'gallgatherv: Sum of send_buffer sizes in receive_counts '//&
         'does not equal total size of receive_buffer.'
    Call assert(Sum(recv_counts) == Size(recv_buf), error_message)

    error_message = 'gallgatherv: Number of rows differs between send and '//&
         'receive buffers. Data will not be correctly ordered.'
    Call assert(Size(send_buf, Dim = 1) == Size(recv_buf, Dim = 1), error_message)

    Call MPI_ALLGATHERV(send_buf,  Size(send_buf),     wp_mpi, &
                        recv_buf,  recv_counts, disps, wp_mpi, &
                        comm%comm, comm%ierr )

  End Subroutine gallgatherv_real_2d

  Subroutine galltoall_integer(comm, sendbuf, scount, recvbuf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 all processes send scount integers from sendbuf to each
    ! processes
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: sendbuf(:), scount
    Integer,          Intent(  Out) :: recvbuf(:)

    Call MPI_ALLTOALL(sendbuf(:), scount, MPI_INTEGER, &
                      recvbuf(:), scount, MPI_INTEGER, &
                      comm%comm, comm%ierr)
  End Subroutine galltoall_integer

  Subroutine galltoallv_integer(comm, sendbuf, scounts, sdisps, &
                                recvbuf, rcounts, rdisps)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 all processes send and receive different amounts of data from
    ! all other processes. The amount send and recieved to each process is
    ! defined by scounts and rcounts respecctively. Displacements are defined by
    ! sdisps and rdisps.
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Integer,          Intent(In   ) :: sendbuf(:), scounts(:), sdisps(:)
    Integer,          Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: rcounts(:), rdisps(:)

    Call MPI_ALLTOALLV(sendbuf(:), scounts(:), sdisps(:), MPI_INTEGER, &
                       recvbuf(:), rcounts(:), rdisps(:), MPI_INTEGER, &
                       comm%comm, comm%ierr)
  End Subroutine galltoallv_integer

  Subroutine gallreduce_logical_scalar(comm, send, recv, op)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 perform a reduction with operation 'op' and share the result
    ! with all processes
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(In   ) :: send
    Logical,          Intent(  Out) :: recv
    Integer,          Intent(In   ) :: op

    Call MPI_ALLREDUCE(send, recv, 1, MPI_LOGICAL, op, &
                       comm%comm, comm%ierr)
  End Subroutine gallreduce_logical_scalar

  Subroutine gallreduce_logical_vector(comm, sendbuf, recvbuf, op)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 perform a reduction with operation 'op' and share the result
    ! with all processes
    !
    ! copyright - daresbury laboratory
    ! author    - j.madge april 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(In   ) :: sendbuf(:)
    Logical,          Intent(  Out) :: recvbuf(:)
    Integer,          Intent(In   ) :: op

    Integer :: n_s

    n_s = Size(sendbuf(:), Dim=1)

    Call MPI_ALLREDUCE(sendbuf(:), recvbuf(:), n_s, MPI_LOGICAL, op, &
                       comm%comm, comm%ierr)
  End Subroutine gallreduce_logical_vector

  Subroutine gallreduce_integer_vector(comm, send_buf, recv_buf, operation)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 perform a reduction operation and share the result
    ! with all processes
    !
    ! copyright - daresbury laboratory
    ! author    - A. Buccheri June 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( comms_type ), Intent( InOut ) :: comm
    Integer,            Intent( In    ) :: send_buf(:)
    Integer,            Intent(   Out ) :: recv_buf(:)
    Integer,            Intent( In    ) :: operation
    Character( Len = 100 )              :: error_message

    error_message = 'gallreduce: Size of send buffer /= size of receive buffer'
    Call assert(Size(send_buf) == Size(recv_buf), error_message)

    Call MPI_ALLREDUCE(send_buf, recv_buf, size(send_buf), MPI_INTEGER, operation, &
                       comm%comm,comm%ierr)

  End Subroutine gallreduce_integer_vector

  Subroutine gallreduce_real_vector(comm, send_buf, recv_buf, operation)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 perform a reduction operation and share the result
    ! with all processes
    !
    ! copyright - daresbury laboratory
    ! author    - A. Buccheri June 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( comms_type ), Intent( InOut ) :: comm
    Real(wp),           Intent( In    ) :: send_buf(:)
    Real(wp),           Intent(   Out ) :: recv_buf(:)
    Integer,            Intent( In    ) :: operation
    Character( Len = 100 )              :: error_message

    error_message = 'gallreduce: Size of send buffer /= size of receive buffer'
    Call assert(Size(send_buf) == Size(recv_buf), error_message)

    Call MPI_ALLREDUCE(send_buf, recv_buf, size(send_buf), wp_mpi, operation, &
                       comm%comm,comm%ierr)

  End Subroutine gallreduce_real_vector


  Subroutine mpi_distribution_init(this, comms)
    Type(comms_type),             Intent(In)    :: comms
    Class(mpi_distribution_type), Intent(InOut) :: this

    Allocate(this%counts( comms%mxnode ))
    Allocate(this%displ( comms%mxnode ))
    this%counts = 0
    this%displ  = 0

  End Subroutine mpi_distribution_init

  Subroutine mpi_distribution_finalise(this)
    Class(mpi_distribution_type), Intent(InOut) :: this
    If(Allocated(this%counts)) Deallocate(this%counts)
    If(Allocated(this%displ))  Deallocate(this%displ)
  End Subroutine mpi_distribution_finalise


  !> @brief Fill counts and displ index arrays required by mpi all/gatherv & all/scatterv
  !!
  !! @param[in]     comm        Object holding MPI settings and comms
  !! @param[in]     N           Size of local array
  !! @param[inout]  counts      1D integer array, size numprocs
  !!                            Holds size of local array on each process
  !! @param[inout]  displ       1D integer array, size numprocs
  !!                            Holds offset of local data w.r.t. position in
  !!                            global array (single, contiguous index)
  !
  Subroutine gatherv_scatterv_index_arrays(comm, N, counts, displ)
    !---------------
    !Declarations
    !---------------
    !Suboutine Arguments
    Type( comms_type ),   Intent( InOut ) :: comm
    Integer,              Intent( In    ) :: N
    Integer, Allocatable, Intent( InOut ) :: displ(:),counts(:)

    !Local data
    Integer, Allocatable :: counts_local(:)
    Integer              :: rank,numprocs,i

    !---------------
    !Main Routine
    !---------------
    rank           = comm%idnode
    numprocs       = comm%mxnode

    Call assert(Allocated(counts), 'gatherv_scatterv_index_arrays: counts not allocated')
    Call assert(Allocated(displ), 'gatherv_scatterv_index_arrays: displ not allocated')

    !Fill counts
    Allocate(counts_local(numprocs))
    counts_local=0
    counts_local(rank+1)=N

    !Add element count from each process into global array
    counts=0
    Call gallreduce(comm, counts_local, counts, op_sum)

    !Fill displ
    displ=0
    Do i=2,numprocs
       displ(i)=displ(i-1)+counts(i-1)
    End Do

  End Subroutine gatherv_scatterv_index_arrays

End Module comms
