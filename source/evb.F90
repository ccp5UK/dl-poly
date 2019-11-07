Module evb
!>
!>  Module declaring global EVB variables, arrays and subroutines
!> 
!>  copyright - daresbury laboratory
!>  author    - i.scivetti june 2019
!>  

  Use kinds, Only : wp,wi
  Use comms, Only : comms_type
  Use particle, Only : corePart

  Use comms,      Only : comms_type,gcheck,grecv, gsend, WriteConf_tag 

  Use errors_warnings, Only : error, info
  Use flow_control, Only : flow_type
  Use configuration, Only : configuration_type
  Use statistics, Only : stats_type
  Use numerics, Only : invert

  Use constants, Only : engunit, eu_ev, eu_kcpm, eu_kjpm, boltz

  Use parse, Only : get_line,get_word,lower_case,clean_string, word_2_real

  Use thermostat, Only : thermostat_type
  Use rigid_bodies, Only : rigid_bodies_type
  Use constraints, Only : constraints_type
  Use pmf, only : pmf_type
  Use core_shell, Only : core_shell_type

  Use filename, Only : file_type, FILE_SETEVB 


  Implicit None
  Private

  Type, Public ::  evb_type
!> EVB force
   Real( Kind = wp ), Allocatable :: force(:,:)
!> Energy shift for fields
   Real( Kind = wp ), Allocatable :: eshift(:)
!> EVB matrix
   Real( Kind = wp ), Allocatable :: matrix(:,:)
!> Matrix to compute EVB stress
   Real( Kind = wp ), Allocatable :: strmat(:,:)
!> Matrix with the derivative of the EVB energy with respect to each component of the three lattice vectors
   Real( Kind = wp ), Allocatable :: dE_da (:,:)
!> EVB stress tensor
   Real( Kind = wp ), Allocatable :: stress(:,:)
!> EVB eigenvalues
   Real( Kind = wp ), Allocatable :: eigval(:)
!> Energy for each force field, including any potential shift 
   Real( Kind = wp ), Allocatable :: eneFF(:)
!> EVB eigenvectors 
   Real( Kind = wp ), Allocatable :: psi(:,:)
!> Coupling parameters between fiedls
!> The functional form for the off-diagonal elements for the EVB matrix (coupling terms) is
!> matrix(i,j)=ac(i,j)*exp[-bc(i,j)*(matrix(i,i)-matrix(j,j))^2]=matrix(j,i)
!> Thus, a(i,j)=a(j,i) and b(i,j)=b(j,i)
   Real( Kind = wp ), Allocatable :: ac(:,:), bc(:,:)
!> Working arrays for diagonalization
  Real( Kind = wp ), Allocatable :: work(:)
  Integer( Kind = wi ), Allocatable :: ifail(:), iwork(:)

  Contains
    Private
    Procedure, Public :: init => allocate_evb_arrays
!>    Final :: deallocate_evb_arrays
  End Type evb_type

  Public :: read_evb
  Public :: evb_pes, evb_setzero
  Public :: evb_checkconfig, evb_merge_stochastic

Contains

          
  Subroutine allocate_evb_arrays(evb,nff)
  Class(evb_type), Intent(InOut) :: evb
  Integer( Kind = wi ), Intent (In ) :: nff
  
  Integer :: fail(1:14)

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
  Allocate (evb%strmat(nff,nff) ,  Stat=fail(12))
  Allocate (evb%dE_da(3,3)      ,  Stat=fail(13))
  Allocate (evb%stress(3,3)     ,  Stat=fail(14))

  If ( Any(fail /= 0 )) Call error(1025)
 
!> Initialise coupling terms and energy shifts
  evb%ac= 0.0_wp
  evb%bc= 0.0_wp

  evb%eshift=0.0_wp

  End Subroutine allocate_evb_arrays

  Subroutine read_evb(evb, flow, files,comm) 
!> !>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!> Read EVB settings and parameters from CONTROL
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type(evb_type), Intent(InOut) :: evb
  Type( flow_type ), Intent(InOut) :: flow
  Type( file_type ), Intent( InOut ) :: files(:)
  Type( comms_type ), Intent( InOut ) :: comm
  Logical                :: carry,safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word

  Integer                :: ncoupl, icoupl
  Integer                :: i, j, ieshift

  Logical, allocatable   :: couplflag(:,:) 
  Logical, allocatable   :: eshiftflag(:) 

  Character( Len = 256 ) :: message, messages(5)
  Character( Len = 256 ) :: evbunit

  ! Allocate matrix for checking we consider all couplings between fields
  Allocate(couplflag(flow%NUM_FF,flow%NUM_FF))
  Allocate(eshiftflag(flow%NUM_FF)) 

  couplflag=.False.
  eshiftflag=.False.

  ! Set the total number of coupling elements for later check
  ncoupl=0
  Do i=1,(flow%NUM_FF-1) 
    ncoupl=ncoupl+1
  End Do

  ! initialise counters
  icoupl=0
  ieshift=0

  ! Set safe flag
  safe=.true.

  ! Open the simulation input file

  If (comm%idnode == 0) Inquire(File=files(FILE_SETEVB)%filename, Exist=safe)
     Call gcheck(comm,safe,"enforce")
  If (.not.safe) Then
    Call error(1092)
  Else
    If (comm%idnode == 0) Then
      Open(Newunit=files(FILE_SETEVB)%unit_no, File=files(FILE_SETEVB)%filename,Status='old')
    End If
  End If
  Call get_line(safe,files(FILE_SETEVB)%unit_no,record,comm)  

  If (safe) Then
    carry = .true.
    Do While (carry)
      Call get_line(safe,files(FILE_SETEVB)%unit_no,record,comm)
      If (.not.safe) Exit
      Call lower_case(record)
      Call get_word(record,word)

!> Read all coupling elements

      If (word(1:8) == 'evbcoupl') Then
        Call get_word(record,word)
        If (word(1:2) == '  ') Then
          Call error(1093)        
        Else    
          i = Nint(word_2_real(word))
          If (i< 1 .or. i > flow%NUM_FF) Then
            Write(message,'(a,i2)') ' ERROR: There is no reference to FIELD ', i 
            Call error(1094)      
          End If        
        End If

        Call get_word(record,word)
        If (word(1:2) == '  ') Then
          Call error(1093)        
        Else    
          j = Nint(word_2_real(word))
          If (j< 1 .or. j > flow%NUM_FF) Then
            Write(message,'(a,i2)') ' ERROR: There is no reference to FIELD ', j 
            Call error(1094)      
          End If        
        End If


        If(i==j) Call error(1094)

        If(couplflag(i,j))Then
          Write(messages(1),'(a,i2,a,i2)') ' ERROR: Coupling between field ',i, ' and ', j
          Write(messages(2),'(a)')        ' cannot be defined more than once. Sorry!'
          Call info(messages,2,.true.)                
          Call error(1095)
        Else  
          couplflag(i,j)=.True.
          couplflag(j,i)=.True.
        End If

        ! Read parameter ac
        Call get_word(record,word)
        If (word(1:2) == '  ') Then
          Call info(message,.true.)                
          Call error(1096)
        Else 
          evb%ac(i,j)= Abs(word_2_real(word))
          evb%ac(j,i)= evb%ac(i,j) 
        End If        

        ! Read parameter bc
        Call get_word(record,word)
        If (word(1:2) == '  ') Then
          Call error(1096)
        Else  
          evb%bc(i,j)= Abs(word_2_real(word))
          evb%bc(j,i)= evb%bc(i,j) 
        End If 

        icoupl=icoupl+1

!> Read the energy shift for each force field

      Else If (word(1:8) == 'evbshift') Then
        Call get_word(record,word)
        If (word(1:2) == '  ') Then
          Call error(1097)        
        End If
        
        i = Nint(word_2_real(word))
        If (i< 1 .or. i > flow%NUM_FF) Then
          Write(message,'(a,i2)') ' ERROR: There is no reference to FIELD  ', i 
          Call info(message,.true.)                
          Call error(1098)
        End If

        If(eshiftflag(i))Then
          Write(message,'(a,i2,a)') ' ERROR: Energy shift for field ',i,' cannot be defined more than once. Sorry!'
          Call info(message,.true.)
          Call error(1099)
        Else
          eshiftflag(i)=.True.
        End If        

        ! Read evb%eshift(i) 
        Call get_word(record,word)
        If (word(1:2) == '  ') Then
          Call error(1097)
        Else  
          evb%eshift(i)= Abs(word_2_real(word))
        End If 

        ieshift=ieshift+1

      Else If (word(1:6) == 'finish') Then
        carry=.false.

      End If

    End Do
  End If

  If(icoupl  /= ncoupl)      Call error(1100)
  If(ieshift /= flow%NUM_FF) Call error(1101)

  Deallocate(eshiftflag, couplflag)

  If (comm%idnode == 0) Close(files(FILE_SETEVB)%unit_no)    

  Call get_word(record,word) ; Call lower_case(word)

!> Recover the label consistent witht the engunit factor   

  If (abs(engunit-eu_ev) <= epsilon(eu_ev)) Then
    evbunit = 'eV'
  Else If (abs(engunit-eu_kcpm) <= epsilon(eu_kcpm)) Then
    evbunit = 'kcal/mol'
  Else If (abs(engunit-eu_kjpm) <= epsilon(eu_kjpm)) Then
    evbunit = 'kJ/mol'
  Else If (abs(engunit-1.0_wp) <= epsilon(1.0_wp)) Then
    evbunit = '10 J/mol'
  Else If (abs(engunit-boltz) <= epsilon(boltz)) Then
    evbunit = 'Kelvin/Boltzmann'
  End If
   
!> Print info about coupling  
  Write(messages(1),'(1x,a)') '------------'
  Write(messages(2),'(1x,a)') 'EVB settings'
  Write(messages(3),'(1x,a)') '------------'
  Call info(messages,3,.true.) 
  Write(messages(1),'(1x,a)') '1) Coupling terms between Force Fields.'
  Write(messages(2),'(1x,a)') 'The adopted functional form for the coupling between fields i and j (C_{ij})' 
  Write(messages(3),'(1x,a)') 'has the following Gaussian expression:' 
  Write(messages(4),'(1x,a)') 'C_{ij}(E_i,E_j)= ac_{ij}*exp[-b_{ij}*(E_i-E_j)^2]'
  Write(messages(5),'(1x,a)') 'Parameters are summarised in the following table:'
  Call info(messages,4,.true.)
  Write(messages(1),'(1x,a)') '______________________________________________________________________'
  Write(messages(2),'(1x,a,20x,a)') 'Force Field pair', 'Parameters'
  Write(messages(3),'(5x,a,5x,a,10x,a,5x,a)') 'i','j', 'ac_{ij} ('//trim(evbunit)//')', 'bc_{ij} ('//trim(evbunit)//')^{-2}' 
  Write(messages(4),'(1x,a)') '______________________________________________________________________'
  Call info(messages,4,.true.)
  Do i=1,flow%NUM_FF
   Do j=i+1,flow%NUM_FF
     Write(message,'(4x,i2,4x,i2,12x,E12.5,15x,E12.5)') i , j, evb%ac(i,j),  evb%bc(i,j)
     Call info(message,.true.)
   End Do
  EnD Do
  Write(message,'(1x,a)') '______________________________________________________________________'
  Call info(message,.true.)
  Call info(' ',.true.)
  Write(messages(1),'(1x,a)') '2) Energy shift for each Force Fields.'
  Write(messages(2),'(1x,a)') 'In addition to the interactions described in each FIELD file, one'
  Write(messages(3),'(1x,a)') 'might want to shift the potential by a certain amount depending on the problem.'
  Write(messages(4),'(1x,a)') 'Energy shifts are described in the following table:'
  Call info(messages,4,.true.)
  Write(messages(1),'(1x,a)') '_______________________________________________________'
  Write(messages(2),'(1x,a,10x,a)') 'Force Field', 'Energy shift ('//trim(evbunit)//')'
  Write(messages(3),'(1x,a)') '_______________________________________________________'
  Call info(messages,3,.true.)
  Do i=1,flow%NUM_FF
    Write(message,'(4x,i2,22x,E12.5)') i , evb%eshift(i)
    Call info(message,.true.)
  EnD Do
  Write(message,'(1x,a)')     '_______________________________________________________'
  Call info(message,.true.)
  Call info(' ',.true.)
  
  ! Convert coupling EVB parameters and energy shifts to internal units
  evb%eshift = engunit*evb%eshift
  evb%ac     = engunit*evb%ac  
  evb%bc     = evb%bc/engunit**2  

  End Subroutine read_evb
          

!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Check all CONFIG files have the same coordinates
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine evb_checkconfig(config,flow,comm)
 
  Type(configuration_type), Intent( In  ) :: config(:)
  Type( flow_type )       , Intent( In  ) :: flow
  Type( comms_type), Intent( InOut ) :: comm

  Integer                :: i, j, m
  Character( Len = 256 ) :: messages(3)
  Character( Len = 1   ) :: coord, coord0
  Character( Len = 6   ) :: string1, string2
  Logical                :: carry

  Integer                :: jdnode

  Real( Kind = wp )      :: cell(flow%NUM_FF,3,3)

  ! Comparison between different CONFIG files: check the value for levcfg 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Do m = 1, flow%NUM_FF-1 
    If(config(m)%levcfg /= config(m+1)%levcfg)then
      If(m == 1) Then       
        write(string1,'(a1)') ' ' 
      Else
        write(string1,'(i2)') m      
      End If        
         write(string2,'(i2)') m+1
      ! In case an inconsistency has been found, complain and abort       
      Write(messages(1),'(a)')  ' Problems: Value for levcfg differs between between'
      Write(messages(2),'(4a)') ' CONFIG'//adjustl(trim(string1)), ' and ','CONFIG'//adjustl(trim(string2)), ' files!'     
      Call info(messages,2,.true.)
      Call error(1106) 
    End If 
  End Do       

  ! Comparison between different CONFIG files: check number of atoms
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Arrange config(1)%cell in a more convenient index notation, stored in cell
  Do m = 1, flow%NUM_FF
    Do i=1,3
      Do j=1,3
        cell(m,i,j)=config(m)%cell(3*(i-1)+j)
      End Do
    End Do
  End Do

  Do i=1,3
    Do j=1,3
      Do m = 1, flow%NUM_FF-1 
        If(abs(cell(m,i,j)-cell(m+1,i,j))>= epsilon(cell(m,i,j)))then
         If(m == 1) Then       
           write(string1,'(a1)') ' ' 
         Else
           write(string1,'(i2)') m      
         End If        
         write(string2,'(i2)') m+1
         ! In case an inconsistency has been found, complain and abort       
         Write(messages(1),'(a,2i2,a)')  ' Problems: Component (', i, j,') of the simulation cell differs between'
         Write(messages(2),'(4a)')       ' CONFIG'//adjustl(trim(string1)), ' and ','CONFIG'//adjustl(trim(string2)), ' files!'     
         Call info(messages,2,.true.)
         Call error(1105) 
        End If 
      End Do       
    End Do
  End Do

  ! Comparison between different CONFIG files: loop over each ionic position
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! set initial flags for comparion
  coord0='w'
  coord=coord0
  carry=.True.

  Do i = 1 , config(1)%mxatdm
     Do m = 1, flow%NUM_FF-1
       If(carry)Then
         ! check x-coordinate of atom i      
         If(abs(config(m)%parts(i)%xxx - config(m+1)%parts(i)%xxx) >= epsilon(config(m)%parts(i)%xxx))then
           coord='x'
         End If
         ! check y-coordinate of atom i      
         If(abs(config(m)%parts(i)%yyy - config(m+1)%parts(i)%yyy) >= epsilon(config(m)%parts(i)%yyy))then
           coord='y'
         End If
         ! check z-coordinate of atom i      
         If(abs(config(m)%parts(i)%zzz - config(m+1)%parts(i)%zzz) >= epsilon(config(m)%parts(i)%zzz))then
           coord='z'
         End If
         If(coord /= coord0)Then
           carry=.False.      
           If(m == 1) Then
             write(string1,'(a1)') ' '
           Else
             write(string1,'(i2)') m      
           End If
           write(string2,'(i2)') m+1 
         End If        
       End If
     End Do
   End Do

  ! Since atoms are divided in different domains and processors we need to pass the flags (computed above) to the master node
   If(comm%idnode == 0)Then
     Do jdnode=0,comm%mxnode-1
       If(jdnode>0)Then
       Call grecv(comm,carry,jdnode,WriteConf_tag)
       Call grecv(comm,coord,jdnode,WriteConf_tag)
       Call grecv(comm,string1,jdnode,WriteConf_tag)
       Call grecv(comm,string2,jdnode,WriteConf_tag)
       End If
       If(.Not. carry)Then
         ! In case an inconsistency has been found, complain and abort       
         Write(messages(1),'(a,i2,7a)')  ' Problems: Process', jdnode, ' found that the ', coord, '-coordinate for one'
         Write(messages(2),'(3a)')      ' of the ions differs between ', 'CONFIG'//adjustl(trim(string1)), ' and '   
         Write(messages(3),'(2a)')      ' CONFIG'//adjustl(trim(string2)),'. This EVB simulation cannot continue.'     
         Call info(messages,3,.true.)
         Call error(1104) 
       End If
     End Do
   Else  
       Call gsend(comm,carry,0,WriteConf_tag)     
       Call gsend(comm,coord,0,WriteConf_tag)     
       Call gsend(comm,string1,0,WriteConf_tag)     
       Call gsend(comm,string2,0,WriteConf_tag)     
   End If  





  Call info(' ',.true.)
  Call info(' EVB check for atomic coordinates and supercell dimensions between different CONFIG files: Passed!',.true.) 
  Call info(' ',.true.)

  End Subroutine evb_checkconfig

!> 
!> This subroutine copy fields of variable types only when regauss is applied
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine evb_merge_stochastic(flow,config,stat,rigid,thermo,cshell,cons,pmf)

  Type( flow_type ), Intent(In   ) :: flow
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type( stats_type), Intent(InOut) :: stat(:)
  Type( rigid_bodies_type), Intent( InOut ) :: rigid(:)
  Type( thermostat_type), Intent( InOut ) :: thermo(:)
  Type( core_shell_type ), Intent( InOut ) :: cshell(:)
  Type( constraints_type ), Intent( InOut ) :: cons(:)
  Type( pmf_type ), Intent( InOut ) :: pmf(:)

  Integer( Kind = wi ) :: m

  Do m=2,flow%NUM_FF

    config(m)%vxx=config(1)%vxx
    config(m)%vyy=config(1)%vyy
    config(m)%vzz=config(1)%vzz
    
    config(m)%parts(:)%fxx=config(1)%parts(:)%fxx
    config(m)%parts(:)%fyy=config(1)%parts(:)%fyy
    config(m)%parts(:)%fzz=config(1)%parts(:)%fzz
    
    config(m)%parts(:)%xxx=config(1)%parts(:)%xxx
    config(m)%parts(:)%yyy=config(1)%parts(:)%yyy
    config(m)%parts(:)%zzz=config(1)%parts(:)%zzz
    
    rigid(m)=rigid(1)
    thermo(m)=thermo(1)
    cshell(m)=cshell(1)
    cons(m)=cons(1)
    pmf(m) = pmf(1)
    stat(m)=stat(1)
  
  End Do 

  End subroutine evb_merge_stochastic
  
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Compute EVB energy, forces and stress tensor
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine evb_pes(evb,flow,config,stat)
 
  Type(evb_type), Intent(InOut) :: evb
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(In   ) :: flow



!> Compute EVB Energy 
!> -----------------
  Call evb_energy(evb,flow,stat)

!> Compute EVB forces
!> -----------------
  Call evb_force(evb,flow,config)

!> Compute EVB stress
!> ------------------
  Call evb_stress(evb,flow,config,stat)

!#endif 

  End Subroutine evb_pes


!> Subroutine to compute the total EVB energy 
  Subroutine evb_energy(evb,flow,stat)
 
  Type(evb_type), Intent(InOut) :: evb
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(In ) :: flow

  Integer( Kind = wi ) :: m,k   ! Indices for matrix elements

  Integer( Kind = wi ) :: mevb, evbinfo

#ifdef EVB
! Initialise matrix elements
  evb%matrix=0.0_wp

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
                  1.d-30, mevb, evb%eigval, evb%psi, flow%NUM_FF, evb%work, 8*flow%NUM_FF,     &
                  evb%iwork, evb%ifail, evbinfo)

    Do m=1,flow%NUM_FF
      stat(m)%stpcfg=evb%eigval(1)
    End Do

#endif    

  End subroutine evb_energy

!> Subroutine to compute the total EVB energy 
  Subroutine evb_force(evb,flow,config)
 
  Type(evb_type), Intent(InOut) :: evb
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type( flow_type ), Intent(In   ) :: flow


  Integer( Kind = wi ) :: m,k     ! Indices for matrix elements
  Integer( Kind = wi ) :: i,j     ! Indices for atoms and/or coordinates

  Integer( Kind = wi ) :: mxatms

! Maximum number of atoms
  mxatms=config(1)%mxatms 

! Here we use the evb%matrix to compute the x,y and z components of the EVB force over each particle 

   Do i=1, mxatms

     ! Copy forces from config%parts(i) to evb%force. x->1, y->2, y->3. We do this for each force field.
     Do m=1,flow%NUM_FF
       evb%force(1,m)=config(m)%parts(i)%fxx
       evb%force(2,m)=config(m)%parts(i)%fyy
       evb%force(3,m)=config(m)%parts(i)%fzz
     End Do
! For particle i, loop over the three coordinates
     Do j=1,3
       evb%matrix=0.0_wp
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


   End subroutine evb_force

!> Subroutine to compute the total EVB energy 
  Subroutine evb_stress(evb,flow,config,stat)
 
  Type(evb_type), Intent(InOut) :: evb
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(In   ) :: flow


  Integer( Kind = wi ) :: m,k     ! Indices for matrix elements
  Integer( Kind = wi ) :: i,j,i2  ! Indices for atoms and/or coordinates

  Integer( Kind = wi ) :: mxatms

  Real( Kind = wp )    :: invcell(1:9), det
  Real( Kind = wp )    :: cell(3,3), rcell(3,3)
  Real( Kind = wp )    :: strff(flow%NUM_FF,3,3)
   

! Maximum number of atoms
  mxatms=config(1)%mxatms 

!> Compute EVB stress
!> ------------------

! Arrange stat(m)%stress in a more convenient index notation, stored in strff. We do that for every field. 
   Do m=1,flow%NUM_FF 
     Do i=1,3
       Do j=1,3
         strff(m,i,j)=stat(m)%stress(3*(i-1)+j)
       End Do
     End Do
   End Do

! Arrange config(1)%cell in a more convenient index notation, stored in cell
   Do i=1,3
     Do j=1,3
       cell(i,j)=config(1)%cell(3*(i-1)+j) 
     End Do
   End Do
     
! Compute the inverse of the lattice vectors   
   Call invert(config(1)%cell,invcell,det)

! Arrange invcell in a more convenient index notation, stored in rcell
   Do i=1,3
     Do j=1,3
       rcell(i,j)=invcell(3*(i-1)+j) 
     End Do
   End Do

! Here we build the evb%dE_da matrix to compute the components of the EVB stress tensor 
! Loop over the 6 independent components of the stress tensor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Do i=1,3
    Do j=1,3
      ! Define matrix components
      ! -----------------------
      ! Diagonal elements
      Do m=1,flow%NUM_FF
        evb%strmat(m,m)= 0.0_wp
        Do i2=1,3  
          evb%strmat(m,m)= evb%strmat(m,m)+strff(m,i,i2)*rcell(j,i2)
        End Do
      End Do
      ! Off-diagonal elements     
      Do m=1,flow%NUM_FF
        Do k=m+1,flow%NUM_FF
          evb%strmat(m,k)=0.0_wp
          Do i2=1,3
            evb%strmat(m,k)=evb%strmat(m,k)+(strff(m,i,i2)-strff(k,i,i2))*rcell(j,i2)
          End Do
          evb%strmat(m,k)=-2.0*evb%ac(m,k)*evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))*evb%strmat(m,k)* & 
                           exp(-evb%bc(m,k)*(evb%eneFF(m)-evb%eneFF(k))**2)
          evb%strmat(k,m)=evb%strmat(m,k)
        End Do 
      End Do
      !------------------------
      evb%dE_da(i,j)=0.0_wp
  
      ! Matrix multiplication (evb%psi)^T*evb%matrix*evb%psi to compute stat(1)%stress(j)
      ! Store EVB stress in stress component of force-field 1
      Do m=1,flow%NUM_FF
        Do k=1,flow%NUM_FF
          evb%dE_da(i,j)=evb%dE_da(i,j)+evb%psi(m,1)*evb%strmat(m,k)*evb%psi(k,1)
        End Do
      End Do
    End Do
  End Do 

! Obtain EVB stress tensor
  Do i=1,3
    Do j=i,3
     ! Set stress to zero
      evb%stress(i,j)=0.0_wp
      Do i2=1,3
        evb%stress(i,j)=evb%stress(i,j)+evb%dE_da(i,i2)*cell(j,i2)
      End Do
        evb%stress(j,i)=evb%stress(i,j)
    End Do
  End Do 

! Copy evb%stress to each FIELD
  Do i=1,3
    Do j=1,3
      Do m=1,flow%NUM_FF  
        stat(m)%stress(3*(i-1)+j)=evb%stress(i,j)
      End Do  
    End Do
  End Do

!> For EVB it is not possible to have a decomposition of the virial into separate contributions for each type of interaction 
!> (e.g. angles, bonds, dihedrals, etc) without recomputing each of the interactions, which would increase the 
!> computational time innecesarily. Here we obtain the total virial from the evb stress tensor computed above.

   Do m=1,flow%NUM_FF
     stat(m)%virtot=-(stat(m)%stress(1)+stat(m)%stress(5)+stat(m)%stress(9))
   End Do

  End subroutine evb_stress

 
!> This subroutine sets to zero most of the decomposed terms of the virial and energy
!> By construction, it is not possible to decompose the EVB virial for the different 
!> types of interactions (bond, angle, dihedral, etc). We thus opt to set them to zero
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine evb_setzero(flow,stat)

  Type( flow_type ), Intent(In   ) :: flow
  Type( stats_type), Intent(InOut) :: stat(:)

  Integer( Kind = wi ) :: ff

  Do ff=1,flow%NUM_FF
!  Energy components
!    stat(ff)%engcpe = 0.0_wp
!    stat(ff)%engsrp = 0.0_wp 
!    stat(ff)%engter = 0.0_wp 
!    stat(ff)%engtbp = 0.0_wp 
!    stat(ff)%engfbp = 0.0_wp 
!    stat(ff)%engshl = 0.0_wp 
!    stat(ff)%engtet = 0.0_wp 
!    stat(ff)%engbnd = 0.0_wp 
!    stat(ff)%engang = 0.0_wp 
!    stat(ff)%engdih = 0.0_wp 
!    stat(ff)%enginv = 0.0_wp
!!  Virial components
!    stat(ff)%vircpe = 0.0_wp 
!    stat(ff)%virsrp = 0.0_wp
!    stat(ff)%virter = 0.0_wp
!    stat(ff)%virtbp = 0.0_wp
!    stat(ff)%virfbp = 0.0_wp
!    stat(ff)%virshl = 0.0_wp
!    stat(ff)%virtet = 0.0_wp 
!    stat(ff)%virbnd = 0.0_wp
!    stat(ff)%virang = 0.0_wp
!    stat(ff)%virdih = 0.0_wp 
!    stat(ff)%virinv = 0.0_wp
  End Do 

  End subroutine evb_setzero

End Module evb        
