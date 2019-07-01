Module evb

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global EVB variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.scietti june 2019
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use comms, Only : comms_type
  Use particle, Only : corePart

  Use errors_warnings, Only : error
  Use flow_control, Only : flow_type
  Use configuration, Only : configuration_type
  Use statistics, Only : stats_type

  Implicit None
  Private

  Type, Public ::  evb_type
! EVB force
   Real( Kind = wp ), Allocatable :: force(:,:)
! Energy shift for fields
   Real( Kind = wp ), Allocatable :: eshift(:)
! EVB matrix
   Real( Kind = wp ), Allocatable :: matrix(:,:)
! EVB eigenvalues
   Real( Kind = wp ), Allocatable :: eigval(:)
! Energy for each force field, including any potential shift 
   Real( Kind = wp ), Allocatable :: eneFF(:)
! EVB eigenvectors 
   Real( Kind = wp ), Allocatable :: psi(:,:)
! Coupling parameters between fiedls
! The functional form for the off-diagonal elements for the EVB matrix (coupling terms) is
! matrix(i,j)=ac(i,j)*exp[-bc(i,j)*(matrix(i,i)-matrix(j,j))^2]=matrix(j,i)
! Thus, a(i,j)=a(j,i) and b(i,j)=b(j,i)
   Real( Kind = wp ), Allocatable :: ac(:,:), bc(:,:)
! Working arrays for diagonalization
  Real( Kind = wp ), Allocatable :: work(:)
  Integer( Kind = wi ), Allocatable :: ifail(:), iwork(:)

! labelling for components of the stress tensor
  Integer( Kind = wi )          :: strtag(1:6)

  Contains
    Private
    Procedure, Public :: init => allocate_evb_arrays
!    Final :: deallocate_evb_arrays
  End Type evb_type

  Public :: evb_pes

Contains

          
  Subroutine allocate_evb_arrays(evb,nff)
  Class(evb_type), Intent(InOut) :: evb
  Integer( Kind = wi ), Intent (In ) :: nff
  
  Integer :: fail(1:11)

  Allocate (evb%eshift(1:nff)   ,  Stat=fail(1))
  Allocate (evb%matrix(nff,nff) ,  Stat=fail(2))
  Allocate (evb%psi(nff,nff)    ,  Stat=fail(3))
  Allocate (evb%ac(nff,nff)     ,  Stat=fail(4))
  Allocate (evb%bc(nff,nff)     ,  Stat=fail(5))
  Allocate (evb%eigval(nff)     ,  Stat=fail(6))
  Allocate (evb%ifail(nff)      ,  Stat=fail(7))
  Allocate (evb%iwork(5*nff)    ,  Stat=fail(8))
  Allocate (evb%work(8*nff)     ,  Stat=fail(9))
  Allocate (evb%eneFF(nff)      ,  Stat=fail(10)) 
  Allocate (evb%force(3,nff)    ,  Stat=fail(11)) 


  If ( Any(fail /= 0 )) Call error(1025)
 
! Initialise coupling terms and energy shifts
  evb%ac=0.0_wp ! 20000.0_wp
  evb%bc=0.0_wp ! 1.0E-14
  evb%eshift(1)=0.0_wp
  evb%eshift(2)=10.0_wp

! Initialise labelling for components of the stress tensor  
  evb%strtag(1)=1; evb%strtag(2)=2; evb%strtag(3)=3 
  evb%strtag(4)=5; evb%strtag(5)=6; evb%strtag(6)=9 


  End Subroutine allocate_evb_arrays


  Subroutine evb_pes(evb,flow,config,stat,comm)
 
  Type(evb_type), Intent(InOut) :: evb
  Type(comms_type), Intent(InOut) ::comm
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(InOut) :: flow


  Integer( Kind = wi ) :: m,k,l    ! Indices for matrix elements
  Integer( Kind = wi ) :: i,j      ! Indices for atoms and/or coordinates

  Integer( Kind = wi ) :: mevb, info, mxatms
 
#ifdef EVB
! Initialise matrx elements
  evb%matrix=0.0_wp

! Maximum number of atoms
  mxatms=config(1)%mxatms 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Temporary checks for consistency between fields. This should be removed in a later stage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   if( abs(stat(1)%stpcfg-stat(2)%stpcfg).gt.1E-12) then
!     print*, 'stpcfg', flow%step, abs(stat(1)%stpcfg-stat(2)%stpcfg)
!     stop
!   elseif(abs(stat(1)%engcpe-stat(2)%engcpe).gt.1E-12) then
!     print*, 'engcpe', flow%step, abs(stat(1)%engcpe-stat(2)%engcpe)
!     stop
!   elseif(abs(stat(1)%engsrp-stat(2)%engsrp).gt.1E-12) then
!     print*, 'engsrp', flow%step, abs(stat(1)%engsrp-stat(2)%engsrp)
!     stop
!   end if   

!!!!!!!!!!!!!!!!!!!!
! Compute EVB Energy 
!!!!!!!!!!!!!!!!!!!!

! Matrix elements
! Diagonal elements
    Do m=1,flow%NUM_FF 
     evb%eneFF(m)   = stat(m)%stpcfg+evb%eshift(m)
     evb%matrix(m,m)= evb%eneFF(m)
    End Do
! Off-diagonal terms
    Do m=1,flow%NUM_FF
      Do k=m+1,flow%NUM_FF
      evb%matrix(m,k) =evb%ac(m,k)*exp(-evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))**2)
        evb%matrix(k,m) =evb%matrix(m,k)
      End Do
    End Do

! Diagonalisation
    call dsyevx( 'V', 'I', 'U', flow%NUM_FF, evb%matrix, flow%NUM_FF, -1., 0., 1, flow%NUM_FF, &
                  1.d-20, mevb, evb%eigval, evb%psi, flow%NUM_FF, evb%work, 8*flow%NUM_FF,     &
                  evb%iwork, evb%ifail, info)

    Do m=1,flow%NUM_FF
     stat(m)%stpcfg=evb%eigval(1)
    End Do

!!!!!!!!!!!!!!!!!!!!
! Compute EVB forces
!!!!!!!!!!!!!!!!!!!!

! Here we use the evb%matrix to compute the x,y and z components of the EVB force over each particle 
   Do i=1, mxatms
     evb%matrix=0.0_wp

     ! Copy forces from config%parts(i) to evb%force. x->1, y->2, y->3. We do this for each foce field.
     Do m=1,flow%NUM_FF
       evb%force(1,m)=config(m)%parts(i)%fxx
       evb%force(2,m)=config(m)%parts(i)%fyy
       evb%force(3,m)=config(m)%parts(i)%fzz
     End Do
! For particle i, loop over the three coordinates
     Do j=1,3
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Define matrix with gradients
       ! Diagonal elements
       Do m=1,flow%NUM_FF
         evb%matrix(m,m)=evb%force(j,m)
       End Do
       ! Off-diagonal elements     
       Do m=1,flow%NUM_FF
         Do k=m+1,flow%NUM_FF
           evb%matrix(m,k)=-2.0*evb%ac(m,k)*evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))*(evb%force(j,m)-evb%force(j,k))* & 
                        exp(-evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))**2)
           evb%matrix(k,m)=evb%matrix(m,k)
         End Do 
       End Do 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Set the first column of evb%force to zero   
       evb%force(j,1)=0.0_wp

       ! Matrix multiplication (evb%psi)^T*evb%matrix*evb%psi to compute evb%force(j,1)
       ! Store EVB force in column 1 of evb%force
       Do m=1,flow%NUM_FF
         Do k=1,flow%NUM_FF
           evb%force(j,1)=evb%force(j,1)+evb%psi(m,1)*evb%matrix(m,k)*evb%psi(k,1)
         End Do
       End Do
! Finish loop over coordinates
     End Do

   ! Copy evb-force to config-force
     Do m=1,flow%NUM_FF
       config(m)%parts(i)%fxx=evb%force(1,1)
       config(m)%parts(i)%fyy=evb%force(2,1)
       config(m)%parts(i)%fzz=evb%force(3,1)
     End Do

   End Do

!!!!!!!!!!!!!!!!!!!!
! Compute EVB stress
!!!!!!!!!!!!!!!!!!!!
! Here we use the evb%matrix to compute the components of the EVB stress tensor 
! Loop over the 6 independent components of the stress tensor
  Do l=1,6
    j=evb%strtag(l)     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define matrix components
    ! Diagonal elements
    Do m=1,flow%NUM_FF
      evb%matrix(m,m)=stat(m)%stress(j)
    End Do
    ! Off-diagonal elements     
    Do m=1,flow%NUM_FF
      Do k=m+1,flow%NUM_FF
        evb%matrix(m,k)=-2.0*evb%ac(m,k)*evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))*(stat(m)%stress(j)-stat(k)%stress(j))* & 
                         exp(-evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))**2)
        evb%matrix(k,m)=evb%matrix(m,k)
      End Do 
    End Do
    
    stat(1)%stress(j)=0.0_wp

    ! Matrix multiplication (evb%psi)^T*evb%matrix*evb%psi to compute stat(1)%stress(j)
    ! Store EVB stress in stress component of force-field 1
    Do m=1,flow%NUM_FF
      Do k=1,flow%NUM_FF
        stat(1)%stress(j)=stat(1)%stress(j)+evb%psi(m,1)*evb%matrix(m,k)*evb%psi(k,1)
      End Do
    End Do

  End Do
! Build a symmetric tenson
  stat(1)%stress(4)=stat(1)%stress(2)
  stat(1)%stress(7)=stat(1)%stress(3)
  stat(1)%stress(8)=stat(1)%stress(6)

! Since the EVB stress was stored in the stress of force-field 1, we now copy to the rest of the fields  
  
   Do m=2,flow%NUM_FF
     stat(m)%stress(:)=stat(1)%stress(:)
   End Do

! In contrast to standard MD, for EVB it is not possible to have a 
! decomposition of the virial into separate contributions for each type of interaction 
! (e.g. angles, bonds, dihedrals, etc) without recomputing each of the interactions, which would
! increase the computational time innecesarily. Here we obtain the total virial from the evb stress tensor
! computed previously.

   Do m=1,flow%NUM_FF
   stat(m)%virtot=-(stat(m)%stress(1)+stat(m)%stress(5)+stat(m)%stress(9))
   End Do
#endif
  End Subroutine evb_pes
End Module evb        
