Module gpfa_wrappers

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module to wrap gpfa without requiring compiler flag
  !   -fallow-argument-mismatch
  !
  ! 
  !
  ! author    - h.l.devereux december 2023
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    Use iso_c_binding, Only: c_f_pointer, c_ptr, c_loc
    Use gpfa235,       Only: gpfa
    Use kinds,         Only: wp

    Interface gpfa_wrap
        ! interface over each rank (:), (:,:,:)
        ! as (*) causes issues with the pointer casting
        ! currently on require rank 1 and 3 as used by
        ! parallel_fft.F90
        Procedure :: gpfa_wrap_rank_1
        Procedure :: gpfa_wrap_rank_3
    End Interface

    Contains 

        Subroutine gpfa_wrap_rank_1(a, trigs, stride, jump, n, n_ffts, direction, start)

            Complex(Kind=wp), Target, Dimension(:), Intent(InOut) :: a
            ! NB trigs are currently always rank 1
            Complex(Kind=wp), Target, Dimension(:), Intent(In   ) :: trigs
            Integer,                                Intent(In   ) :: stride, jump, &
                                                                     n, n_ffts, &
                                                                     direction, start

            
            Type(c_ptr)                          :: c_ptr_a, c_ptr_trigs
            Real(Kind=wp), Pointer, Dimension(:) :: f_ptr_a, f_ptr_trigs
            
            ! address the Size(a) Complex data values by a C pointer
            c_ptr_a = c_loc(a(1))
            ! reinterpret the Size(a) Complex data values as
            ! 2*Size(a) Real data values, assigining to a
            ! Real Fortran pointer
            Call c_f_pointer(c_ptr_a, f_ptr_a, shape=[Size(a)*2])

            c_ptr_trigs = c_loc(trigs(1))
            Call c_f_pointer(c_ptr_trigs, f_ptr_trigs, shape=[Size(trigs)*2])
            
            Call gpfa(f_ptr_a(start:), f_ptr_a(start+1:), f_ptr_trigs, stride, jump, n, n_ffts, direction)
        
        End Subroutine

        Subroutine gpfa_wrap_rank_3(a, trigs, stride, jump, n, n_ffts, direction, start)

            Complex(Kind=wp), Target, Dimension(:,:,:), Intent(InOut) :: a
            ! NB trigs are currently always rank 1
            Complex(Kind=wp), Target, Dimension(:),     Intent(In   ) :: trigs
            Integer,                                    Intent(In   ) :: stride, jump, &
                                                                         n, n_ffts, &
                                                                         direction, start

            
            Type(c_ptr)                              :: c_ptr_a, c_ptr_trigs
            Real(Kind=wp), Pointer, Dimension(:)     :: f_ptr_a
            Real(Kind=wp), Pointer, Dimension(:)     :: f_ptr_trigs
            
            ! address the N Complex data values by a C pointer
            c_ptr_a = c_loc(a(1,1,1))
            ! reinterpret the N Complex data values as
            ! 2*N Real data values, assigining to a
            ! Real Fortran pointer
            Call c_f_pointer(c_ptr_a, f_ptr_a, shape=[Size(a,1)*2*Size(a,2)*2*Size(a,3)*2])

            c_ptr_trigs = c_loc(trigs(1))
            Call c_f_pointer(c_ptr_trigs, f_ptr_trigs, shape=[Size(trigs,1)*2])
            
            Call gpfa(f_ptr_a(start:), f_ptr_a(start+1:), f_ptr_trigs, stride, jump, n, n_ffts, direction)
        
        End Subroutine


End Module