Module evb
!>
!>  Module declaring global EVB variables, arrays and subroutines
!> 
!>  copyright - daresbury laboratory
!>  author    - i.scivetti june 2019
!>  

  Use kinds, Only : wp,wi
  Use comms, Only : comms_type, gsync, gbcast
  Use particle, Only : corePart

  Use comms,      Only : comms_type,gcheck,grecv, gsend, WriteConf_tag 

  Use errors_warnings, Only : error, info
  Use flow_control, Only : flow_type, RESTART_KEY_OLD
  Use configuration, Only : configuration_type
  Use statistics, Only : stats_type
  Use numerics, Only : invert

  Use constants, Only : engunit, eu_ev, eu_kcpm, eu_kjpm, boltz

  Use parse, Only : get_line,get_word,lower_case,clean_string, word_2_real, strip_blanks

  Use thermostat, Only : thermostat_type
  Use rigid_bodies, Only : rigid_bodies_type
  Use constraints, Only : constraints_type
  Use pmf, only : pmf_type
  Use core_shell, Only : core_shell_type

  Use filename, Only : file_type, FILE_SETEVB, FILE_POPEVB 

  Use site, Only : site_type

  Implicit None
  Private

  Type, Public ::  evb_type
!> Flag for activating printing of EVB population    
   Logical                        :: population = .false.
!> Flag for opening EVB population file     
   Logical                        :: population_file_open = .false.
!> Flag for newjob
   Logical                        :: newjob = .true.
!> EVB force
   Real( Kind = wp ), Allocatable :: force(:,:)
!> Energy shift for fields
   Real( Kind = wp ), Allocatable :: eshift(:)
!> EVB Energy matrix
   Real( Kind = wp ), Allocatable :: ene_matrix(:,:)
!> EVB Force matrix
   Real( Kind = wp ), Allocatable :: force_matrix(:,:)
!> Matrix to compute EVB stress
   Real( Kind = wp ), Allocatable :: stress_matrix(:,:)
!> Matrix with the derivative of the EVB energy with respect to each component of the three lattice vectors
   Real( Kind = wp ), Allocatable :: dE_dh (:,:)
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
   Integer                             :: maxfitparam=7        !  maximum number of fitting parameters
   Real( Kind = wp )     , Allocatable :: couplparam(:,:,:)    !  Array with coupling parameters
   Real( Kind = wp )     , Allocatable :: grad_coupl(:,:)      !  Array with coupling parameters
   Character( Len = 5  ) , Allocatable :: typcoupl(:,:)
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
  Public :: evb_check_configs, evb_merge_stochastic, evb_population, evb_check_topology

Contains

          
  Subroutine allocate_evb_arrays(evb,nff)
  Class(evb_type), Intent(InOut) :: evb
  Integer( Kind = wi ), Intent (In ) :: nff
  
  Integer :: fail(1:16)


  Allocate (evb%eshift(1:nff)            ,  Stat=fail(1))
  Allocate (evb%ene_matrix(nff,nff)      ,  Stat=fail(2))
  Allocate (evb%psi(nff,nff)             ,  Stat=fail(3))
  Allocate (evb%eigval(nff)              ,  Stat=fail(4))
  Allocate (evb%ifail(nff)               ,  Stat=fail(5))
  Allocate (evb%iwork(5*nff)             ,  Stat=fail(6))
  Allocate (evb%work(8*nff)              ,  Stat=fail(7))
  Allocate (evb%eneFF(nff)               ,  Stat=fail(8)) 
  Allocate (evb%force(3,nff)             ,  Stat=fail(9)) 
  Allocate (evb%force_matrix(nff,nff)    ,  Stat=fail(10))
  Allocate (evb%dE_dh(3,3)               ,  Stat=fail(11))
  Allocate (evb%stress(3,3)              ,  Stat=fail(12))
  Allocate (evb%stress_matrix(nff,nff)   ,  Stat=fail(13))
  Allocate (evb%couplparam(nff,nff,evb%maxfitparam)  ,  Stat=fail(14))
  Allocate (evb%typcoupl(nff,nff)                    ,  Stat=fail(15))
  Allocate (evb%grad_coupl(nff,nff)      ,  Stat=fail(16))

  If ( Any(fail /= 0 )) Call error(1025)
 
!> Initialise coupling terms and energy shifts
  evb%couplparam= 0.0_wp
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
  Integer                :: i, j, k, kparam, ieshift

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

      If (word(1:1) == '#' .or. word(1:3) == '   ') Then

!> Read setting to print EVB population
      Else If (word(1:6) == 'evbpop') Then
        evb%population = .True.

!> Read all coupling elements

      Else If(word(1:8) == 'evbcoupl') Then
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

        ! Read type of functional form for the coupling term
        Call get_word(record,evb%typcoupl(i,j))
        If (word(1:2) == '  ') Then
          Call error(1096)
        Else 
          evb%typcoupl(j,i)= evb%typcoupl(i,j) 
        End If        

        ! Read parameters for coupling. The number of parameters depend on the type, as follows 
        If (evb%typcoupl(i,j) == 'const') Then
          kparam=1      
          Do k=1, kparam        
            Call get_word(record,word)
            If (word(1:2) == '  ') Then
              Call error(1096)
            Else
              evb%couplparam(i,j,k)= word_2_real(word)
              evb%couplparam(j,i,k)=evb%couplparam(i,j,k)
            End If
          End Do
        Else If (evb%typcoupl(i,j) == 'gauss') Then
          kparam=4      
          Do k=1, kparam        
            Call get_word(record,word)
            If (word(1:2) == '  ') Then
              Call error(1096)
            Else
              evb%couplparam(i,j,k)= word_2_real(word)
              evb%couplparam(j,i,k)= evb%couplparam(i,j,k)
            End If
          End Do
        Else If (evb%typcoupl(i,j) == 'gdump') Then
          kparam=evb%maxfitparam
          Do k=1, kparam        
            Call get_word(record,word)
            If (word(1:2) == '  ') Then
              Call error(1096)
            Else
              evb%couplparam(i,j,k)= word_2_real(word)
              evb%couplparam(j,i,k)= evb%couplparam(i,j,k)
            End If
          End Do
        Else 
          Call error(1096)  
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
          evb%eshift(i)= word_2_real(word)
        End If 

        ieshift=ieshift+1

      Else If (word(1:6) == 'finish') Then
        carry=.false.

      Else 
        Call strip_blanks(record)
        Write(message,'(2a)') word(1:Len_Trim(word)+1),record
        Call info(message,.true.)
        call error(1102)

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
  Write(messages(3),'(1x,a)') 'can be of three types: constant (const), Gaussian (gauss) or Gauss-dumping (gdump)'
  Write(messages(4),'(1x,a)') 'See manual for details on the functional forms'
  Write(messages(5),'(1x,a)') 'Details for coupling terms are summarised in the following table:'
  Call info(messages,5,.true.)
  Write(messages(1),'(1x,a)')               '______________________________________________________________________'
  Write(messages(2),'(1x,a,5x,a,5x,a)')           'Force Field pair' , 'Type' , 'Parameters ('//trim(evbunit)//')'
  Write(messages(3),'(5x,a,5x,a,10x,a,7(15x,a))')             'i','j', '    '  , 'A1', 'A2','A3','A4', 'A5', 'A6','A7'
  Write(messages(4),'(1x,a)')               '______________________________________________________________________'
  Call info(messages,4,.true.)
  Do i=1,flow%NUM_FF
   Do j=i+1,flow%NUM_FF
     If(evb%typcoupl(i,j)=='const')Then
       Write(message,'(2(4x,i2),10x,a,10x,1(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%couplparam(i,j,k),k=1,1)
       Call info(message,.true.)
     ElseIf(evb%typcoupl(i,j)=='gauss')Then
       Write(message,'(2(4x,i2),10x,a,10x,4(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%couplparam(i,j,k),k=1,4)
       Call info(message,.true.)
     ElseIf(evb%typcoupl(i,j)=='gdump')Then
       Write(message,'(2(4x,i2),10x,a,10x,7(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%couplparam(i,j,k),k=1,7)
       Call info(message,.true.)
     End If  
   End Do
  EnD Do
  Write(message,'(1x,a)') '______________________________________________________________________'
  Call info(message,.true.)
  Call info(' ',.true.)
  Write(messages(1),'(1x,a)') '2) Energy shift for each Force Fields.'
  Write(messages(2),'(1x,a)') 'In addition to the interactions described in each FIELD file, one'
  Write(messages(3),'(1x,a)') 'might want to shift the potential by a certain amount'
  Write(messages(4),'(1x,a)') 'Energy shifts are reported in the following table:'
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
  evb%eshift     = engunit*evb%eshift
  evb%couplparam = engunit*evb%couplparam  

  End Subroutine read_evb
          


!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Check all CONFIG files have the same number of atoms, coordinates and simulation cell
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine evb_check_configs(config,flow,comm)
 
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


  Call info(' ',.true.)
  Call info(' EVB check for consistency of atomic coordinates and supercell dimensions between different CONFIG files....',.true.) 
  Call info(' ',.true.)

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
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, 
  Do m = 1, flow%NUM_FF-1
    If(config(m)%mxatdm /= config(m+1)%mxatdm)Then
      Write(messages(1),'(a,2i2,a)')  ' Different number of atoms in CONFIG files for FFs', m, ' and ', m+1
      Call info(messages,1,.true.)
      Call error(1108)
    EndIf        
  End Do  

  ! Comparison between different CONFIG files: check simulation cell 
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
  
  Do i = 1 , config(1)%natms
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


  Call info(' Check Passed!',.true.) 
  Call info(' ',.true.)

  End Subroutine evb_check_configs


!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Check consistency in the topology for all the FIELD  files, (frozen atoms, rigid bonds, bond constraints, etc) 
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine evb_check_topology(config,cons,cshell,sites,flow,rigid,comm)

  Type( configuration_type ), Intent( In    ) :: config(:)
  Type( constraints_type )  , Intent( In    ) :: cons(:)
  Type( core_shell_type )   , Intent( In    ) :: cshell(:)
  Type( site_type )         , Intent( In    ) :: sites(:)
  Type( flow_type )         , Intent( In    ) :: flow
  Type( rigid_bodies_type)  , Intent( In    ) :: rigid(:)        
  Type( comms_type)         , Intent( InOut ) :: comm

  Integer                :: m, i, j, l, jdnode
  Integer                :: mol1, mol2   

  Character( Len = 256 ) :: messages(4)

  Logical                :: loopi, loopj
  Logical                :: consmiss, consdist, conserror
  Logical                :: cshlmiss, cshldist, cshlerror

  Integer                :: st1(2)
  Integer                :: st2(2)
  Integer                :: mol(2)

  Real( Kind = wp )      :: dcons1, dcons2
  Real( Kind = wp )      :: kshl1, kshl2

  Call info(' ',.true.)
  Call info(' EVB check for consistency in the topology for all the FIELD files (rigid bodies, & 
           &bond constraints, frozen atoms, etc)',.true.) 
  Call info(' ',.true.)

  !Initilization of vectors
  mol=0
  st1=0
  st2=0

  ! Check consistency in the number of rigid bodies
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Do m = 1, flow%NUM_FF-1
    If(rigid(m)%total /= rigid(m+1)%total)Then
      Write(messages(1),'(a,i2,a,i2)')  ' Problems! Different number of Rigid Bodies for FFs', m, ' and ', m+1
      Call info(messages,1,.true.)
      Call error(1108)
    EndIf        
  End Do  

  
  ! Check consistency in the frozen nature of ions 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Do i = 1 , config(1)%natms
    Do m = 1, flow%NUM_FF-1
      If(config(m)%lfrzn(i) /= config(m+1)%lfrzn(i))Then
        Write(messages(1),'(a)')          ' Problems!' 
        Write(messages(2),'(a)')          ' Ions change from frozen to free (or vice versa)' 
        Write(messages(3),'(a,i2,a,i2)')  ' between FFs ', m, ' and ', m+1
        Call info(messages,2,.true.)
        Call error(1109)
      End If
    End Do
  End Do

  ! Check consistency in the definition of bond constraints 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, check that the number of bonds constraints is the same

  Do m = 1, flow%NUM_FF-1
    If( cons(m)%megcon /= cons(m+1)%megcon )Then
      Write(messages(1),'(a,i2)')       ' Problems! Number bond constraints in FIELD file', m
      Write(messages(2),'(a,i2)')       ' **DIFFERS** from the number of constraints set in FIELD file ', m+1
      Write(messages(3),'(a)')          ' Please revise the settings of constraints in the FIELD files'
      Call info(messages,3,.true.)
      Call error(1110)
    End If  
  End Do

  ! Second, check that bond constraint specification does not change beetween files
  ! Perform the checking only if the number of bond contraints is different from zero 
  If(cons(1)%megcon /= 0)Then

    consmiss=.False.
    consdist=.False.
    conserror=.False.
  
    Do m = 1, flow%NUM_FF-1
      loopi=.True.
      i=1
      Do While (i <= cons(m)%ntcons .And. loopi)
        loopj=.True.
        j=1
          Do While (j <= cons(m+1)%ntcons .And. loopj)
              If((cons(m)%listcon(1,i) == cons(m+1)%listcon(1,j)) .or. (cons(m)%listcon(1,i) == cons(m+1)%listcon(2,j)) )Then
                If((cons(m)%listcon(2,i) == cons(m+1)%listcon(1,j)) .or. (cons(m)%listcon(2,i) == cons(m+1)%listcon(2,j)) )Then
                  loopj=.False.
                  dcons1=cons(m)%prmcon(cons(m)%listcon(0,i))
                  dcons2=cons(m+1)%prmcon(cons(m+1)%listcon(0,j))
                  If(abs(dcons1-dcons2) .gt. epsilon(dcons1))Then
                    loopi=.False.
                    consdist =.True. 
                    Call obtain_sites_from_constraints(i,m  ,sites,cons,st1,mol1)
                    Call obtain_sites_from_constraints(j,m+1,sites,cons,st2,mol2)
                    mol(1)=mol1
                    mol(2)=mol2
                  End If  
                End If  
              End If
            j=j+1
          End do
  
          If(loopj)Then
            loopi=.False.
            consmiss=.True.
            Call obtain_sites_from_constraints(i,m,sites,cons,st1,mol1)
            mol(1)=mol1
          End If        
  
          i=i+1
      End Do
  
      If(comm%idnode == 0)Then
  
        Do jdnode=0,comm%mxnode-1
          If(jdnode>0)Then
            Call grecv(comm,consmiss,jdnode,WriteConf_tag)
            Call grecv(comm,consdist,jdnode,WriteConf_tag)
            Call grecv(comm,st1,jdnode,WriteConf_tag)
            Call grecv(comm,st2,jdnode,WriteConf_tag)
            Call grecv(comm,mol,jdnode,WriteConf_tag)
          End If
  
          If(consmiss) Then
            Write(messages(1),'(a)')                     ' Problems! '
            Write(messages(2),'(a,i4,a,i4,a,i2,a,i2,a)') ' Bond constraint between atomic sites', st1(1), ' and ' , st1(2), &
                                                           ' of molecule type ', mol(1), ' (set in  FF ', m, ')'
            Write(messages(3),'(a,i2,a)')                ' could not find its equivalent in FF ', m+1, '.' 
            Write(messages(4),'(a)')                     ' Please revise the settings of constraints in the FIELD files' 
            Call info(messages,4,.True.)
            conserror=.True.
  
          Else If(consdist) Then
            
            Write(messages(1),'(a)')                     ' Problems! '
            Write(messages(2),'(a,i6,a,i6,a,i2,a,i2,a)') ' Constraint distance between atomic site       ', st1(1), & 
                                                           ' and ' , st1(2) ,' of molecule type ', mol(1),        & 
                                                           ' (set in FF ', m,')'
            Write(messages(3),'(a,i6,a,i6,a,i2,a,i2,a)') ' **DIFFERS** from the distance between atomic sites',     & 
                                                             st2(1), ' and ', st2(2),           &
                                                           ' of molecule type ', mol(2), ' (set in FF ', m+1,')'
            Write(messages(4),'(a)')                     ' Please revise the settings of constraints in the FIELD files.' 
            Call info(messages,4,.True.)
  
            conserror=.True.
          End If
          If(conserror)Then
            Call error(1110)
          End if  
        End Do
      Else  
        Call gsend(comm,consmiss,0,WriteConf_tag)     
        Call gsend(comm,consdist,0,WriteConf_tag)     
        Call gsend(comm,st1,0,WriteConf_tag)     
        Call gsend(comm,st2,0,WriteConf_tag)     
        Call gsend(comm,mol,0,WriteConf_tag)     
      End If 
  
    End Do

  End If

  ! Check consistency in the definition of core-shell units 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! First, check that the number of core-shell units is the same for all FIELD files

  Do m = 1, flow%NUM_FF-1
    If( cshell(m)%megshl /= cshell(m+1)%megshl )Then
      Write(messages(1),'(a,i2)')       ' Problems! Total number of core-shell units in FIELD file', m
      Write(messages(2),'(a,i2)')       ' **DIFFERS** from the number of core-shell units in FIELD file ', m+1
      Write(messages(3),'(a)')          ' Please revise the settings core-shell in the FIELD files'
      Call info(messages,3,.true.)
      Call error(1111)
    End If  
  End Do

  ! Second, check that core-shell specification does not change beetween files
  ! Perform the checking only if the number of core-shell units is different from zero 

  If(cshell(1)%megshl /= 0)Then

    cshlmiss=.False.
    cshldist=.False.
    cshlerror=.False.
  
    Do m = 1, flow%NUM_FF-1
      loopi=.True.
      i=1
      Do While (i <= cshell(m)%ntshl .And. loopi)
        loopj=.True.
        j=1
        Do While (j <= cshell(m+1)%ntshl .And. loopj)
          If((cshell(m)%listshl(1,i) == cshell(m+1)%listshl(1,j)) .or. (cshell(m)%listshl(1,i) == cshell(m+1)%listshl(2,j)) )Then
            If((cshell(m)%listshl(2,i) == cshell(m+1)%listshl(1,j)) .or. (cshell(m)%listshl(2,i) == cshell(m+1)%listshl(2,j)) )Then
              loopj=.False.
              Do l=1,2
                kshl1=cshell(m)%prmshl(l,cshell(m)%listshl(0,i))
                kshl2=cshell(m+1)%prmshl(l,cshell(m+1)%listshl(0,j))
                If(abs(kshl1-kshl2) .gt. epsilon(kshl1))Then
                  loopi=.False.
                  cshldist =.True. 
                  Call obtain_sites_from_coreshell(i,m  ,sites,cshell,st1,mol1)
                  Call obtain_sites_from_coreshell(j,m+1,sites,cshell,st2,mol2)
                  mol(1)=mol1
                  mol(2)=mol2
                End If
              End Do
             End If  
          End If
          j=j+1
        End do
  
        If(loopj)Then
          loopi=.False.
          cshlmiss=.True.
          Call obtain_sites_from_coreshell(i,m,sites,cshell,st1,mol1)
          mol(1)=mol1
         End If        
      i=i+1
    End Do
  
      If(comm%idnode == 0)Then
  
        Do jdnode=0,comm%mxnode-1
          If(jdnode>0)Then
            Call grecv(comm,cshlmiss,jdnode,WriteConf_tag)
            Call grecv(comm,cshldist,jdnode,WriteConf_tag)
            Call grecv(comm,st1,jdnode,WriteConf_tag)
            Call grecv(comm,st2,jdnode,WriteConf_tag)
            Call grecv(comm,mol,jdnode,WriteConf_tag)
          End If
  
          If(cshlmiss) Then
            Write(messages(1),'(a)')                     ' Problems! '
            Write(messages(2),'(a,i4,a,i4,a,i2,a,i2,a)') ' Core-Shell unit between atomic sites', st1(1), ' and ' , st1(2), &
                                                         ' of molecule type ', mol(1), ' (set in  FF ', m, ')'
            Write(messages(3),'(a,i2,a)')                ' could not find its equivalent in FF ', m+1, '.' 
            Write(messages(4),'(a)')                     ' Please check the existence of these units in the FIELD files' 
            Call info(messages,4,.True.)
            cshlerror=.True.
  
          Else If(cshldist) Then
            
            Write(messages(1),'(a)')                     ' Problems! '
            Write(messages(2),'(a,i6,a,i6,a,i2,a,i2,a)') ' Core-Shell unit specification between atomic site ', st1(1), & 
                                                         ' and ' , st1(2) ,' of molecule type ', mol(1),        & 
                                                         ' (set in FF ', m,') **DIFFERS** from'
            Write(messages(3),'(a,i6,a,i6,a,i2,a,i2,a)') ' the core-shell especification between atomic sites', & 
                                                           st2(1), ' and ', st2(2),           &
                                                         ' of molecule type ', mol(2), ' (set in FF ', m+1,')'
            Write(messages(4),'(a)')                     ' Please revise the force constants for these units in the FIELD files' 
            Call info(messages,4,.True.)
  
            cshlerror=.True.
          End If
          If(cshlerror)Then
            Call error(1111)
          End if  
        End Do
      Else  
        Call gsend(comm,cshlmiss,0,WriteConf_tag)     
        Call gsend(comm,cshldist,0,WriteConf_tag)     
        Call gsend(comm,st1,0,WriteConf_tag)     
        Call gsend(comm,st2,0,WriteConf_tag)     
        Call gsend(comm,mol,0,WriteConf_tag)     
      End If 
  
    End Do

  End If



  Call info(' Check Passed!',.true.) 
  Call info(' ',.true.)

  End Subroutine evb_check_topology

  Subroutine obtain_sites_from_constraints(ind,field,sites,cons,st,mol) 
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Subroutine to compute the atomic sites of a molecule type
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   Integer, Intent( In    )                    :: ind   
   Integer, Intent( In    )                    :: field   
   Type( site_type )         , Intent( In    ) :: sites(:)
   Type( constraints_type ), Intent( In    )   :: cons(:)
   Integer, Intent(   Out )                    :: st(2) 
   Integer, Intent(   Out )                    :: mol   

   Integer :: i, j, k
   Integer :: atcons, nsite
   Integer :: num, den       

   Logical :: flag


   Do i=1,2
     flag=.True.
     atcons=cons(field)%listcon(i,ind)

     Do k=sites(field)%ntype_mol-1,1,-1
      If(flag)Then
        nsite=0
        Do j=k,1,-1
          nsite=nsite+sites(field)%num_mols(j)*sites(field)%num_site(j)
        End Do
        If(nsite < atcons)Then
          num=atcons-nsite
          den=sites(field)%num_site(k+1)
          st(i) = mod(num,den)
          If(st(i)==0)Then
            st(i)=sites(field)%num_site(k+1)
          End If        
          mol=k+1
          flag=.False.
        End If  
       End If
     End Do

     If(flag)Then
       num=atcons
       den=sites(field)%num_site(1)      
       st(i)=mod(num,den)
       If(st(i)==0)Then
         st(i)=sites(field)%num_site(1)
       End If
       mol=1     
     End If        

   End Do

   End subroutine  obtain_sites_from_constraints

  Subroutine obtain_sites_from_coreshell(ind,field,sites,cshell,st,mol) 
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Subroutine to compute the atomic sites of a molecule type
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   Integer, Intent( In    )                    :: ind   
   Integer, Intent( In    )                    :: field   
   Type( site_type )         , Intent( In    ) :: sites(:)
   Type( core_shell_type )   , Intent( In    ) :: cshell(:)
   Integer, Intent(   Out )                    :: st(2) 
   Integer, Intent(   Out )                    :: mol   

   Integer :: i, j, k
   Integer :: atshell, nsite
   Integer :: num, den       

   Logical :: flag


   Do i=1,2
     flag=.True.
     atshell=cshell(field)%listshl(i,ind)

     Do k=sites(field)%ntype_mol-1,1,-1
      If(flag)Then
        nsite=0
        Do j=k,1,-1
          nsite=nsite+sites(field)%num_mols(j)*sites(field)%num_site(j)
        End Do
        If(nsite < atshell)Then
          num=atshell-nsite
          den=sites(field)%num_site(k+1)
          st(i) = mod(num,den)
          If(st(i)==0)Then
            st(i)=sites(field)%num_site(k+1)
          End If        
          mol=k+1
          flag=.False.
        End If  
       End If
     End Do

     If(flag)Then
       num=atshell
       den=sites(field)%num_site(1)      
       st(i)=mod(num,den)
       If(st(i)==0)Then
         st(i)=sites(field)%num_site(1)
       End If
       mol=1     
     End If        

   End Do

   End subroutine  obtain_sites_from_coreshell

 
   Subroutine evb_pes(evb,flow,config,stat)
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Compute EVB energy, forces and stress tensor
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  Type(evb_type), Intent(InOut) :: evb
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(In   ) :: flow



!> Compute EVB Energy 
!> -----------------
  Call evb_energy(evb,flow,stat)

!> Compute gradients for the coupling terms of the EVB matrix
!> -----------------
  Call grad_evb_couplings(flow,evb)

!> Compute EVB forces
!> -----------------
  Call evb_force(evb,flow,config)

!> Compute EVB stress
!> ------------------
  Call evb_stress(evb,flow,config,stat)

!#endif 

  End Subroutine evb_pes


  Subroutine evb_energy(evb,flow,stat)
!> Subroutine to compute the total EVB energy 
 
  Type(evb_type), Intent(InOut) :: evb
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(In ) :: flow

  Integer( Kind = wi ) :: m,k   ! Indices for matrix elements

  Integer( Kind = wi ) :: mevb, evbinfo

#ifdef EVB
! Initialise matrix elements
  evb%ene_matrix=0.0_wp

! Matrix elements
! Diagonal elements
    Do m=1,flow%NUM_FF 
     evb%eneFF(m)   = stat(m)%stpcfg+evb%eshift(m)
     evb%ene_matrix(m,m)= evb%eneFF(m)
    End Do
! Off-diagonal terms
    Do m=1,flow%NUM_FF
      Do k=m+1,flow%NUM_FF
        ! Compute coupling, element (k,m) of the EVB matrix   
        Call evb_energy_couplings(m,k,evb)
        ! Make EVB matrix symmetric
        evb%ene_matrix(k,m) =evb%ene_matrix(m,k)
      End Do
    End Do
  
! Diagonalisation
    call dsyevx( 'V', 'I', 'U', flow%NUM_FF, evb%ene_matrix, flow%NUM_FF, -1., 0., 1, flow%NUM_FF, &
                  1.d-30, mevb, evb%eigval, evb%psi, flow%NUM_FF, evb%work, 8*flow%NUM_FF,     &
                  evb%iwork, evb%ifail, evbinfo)

    Do m=1,flow%NUM_FF
      stat(m)%stpcfg=evb%eigval(1)
    End Do

#endif    

  End subroutine evb_energy

!> Compute coupling energy for matrix element (m,k) 
  Subroutine evb_energy_couplings(m,k,evb)

  Type(evb_type), Intent(InOut) :: evb
  Integer       , Intent(In   ) :: m,k

  Real( Kind = wp )             :: ediff
  Real( Kind = wp )             :: A(evb%maxfitparam)
  
  ! Copy evb coupling parameters to array A, for the sake of clarity in coding the formulae
  A(:)  = evb%couplparam(m,k,:)
  ! Define energy difference, which is the reaction corrdinate between fields m and k
  ediff = evb%eneFF(m)-evb%eneFF(k)
 
  If(evb%typcoupl(m,k)=='const')Then
    evb%ene_matrix(m,k) = A(1)
  Else If (evb%typcoupl(m,k)=='gauss')Then
    evb%ene_matrix(m,k) =  A(1)*exp(-((ediff-A(2))/A(3))**2)+A(4)
  Else If (evb%typcoupl(m,k)=='gdump')Then
    evb%ene_matrix(m,k) = A(1)*exp(-((ediff-A(2))/A(3))**2)*(1.0d0+tanh((ediff+A(4))/A(5)))*(1.0d0+tanh((-ediff+A(6))/A(7)))      
  End If        


  End Subroutine evb_energy_couplings

!> Compute gradient of coupling terms for the EVB matrix
  Subroutine grad_evb_couplings(flow,evb)

  Type( flow_type ), Intent(In )   :: flow
  Type( evb_type  ), Intent(InOut) :: evb
  Integer                          :: m,k

  Real( Kind = wp )                :: ediff
  Real( Kind = wp )                :: A(evb%maxfitparam)
  Real( Kind = wp )                :: t1, t2, t3, grad_t1, grad_t2, grad_t3

  evb%grad_coupl=0.0d0

  Do m=1,flow%NUM_FF
    Do k=m+1,flow%NUM_FF
    ! Copy evb coupling parameters to array A, for the sake of clarity in coding the formulae
      A(:)  = evb%couplparam(m,k,:)
    ! Define energy difference, which is the reaction corrdinate between fields m and k
      ediff = evb%eneFF(m)-evb%eneFF(k)

      If(evb%typcoupl(m,k)=='const')Then
        evb%grad_coupl(m,k) = 0.0d0 
      Else If (evb%typcoupl(m,k)=='gauss')Then
        evb%grad_coupl(m,k) = -2.0* (A(1)/A(3)**2)*(ediff-A(2))*exp(-((ediff-A(2))/A(3))**2)
      Else If (evb%typcoupl(m,k)=='gdump')Then
        t1  = (A(1)*exp(-((ediff-A(2))/A(3))**2))
        t2  = (1.0d0+tanh(( ediff+A(4))/A(5)))
        t3  = (1.0d0+tanh((-ediff+A(6))/A(7)))
        grad_t1= -2.0* (A(1)/A(3)**2)*(ediff-A(2))*exp(-((ediff-A(2))/A(3))**2)
        grad_t2= 1.0/A(5)/(cosh(( ediff+A(4))/A(5)))**2
        grad_t3=-1.0/A(7)/(cosh((-ediff+A(6))/A(7)))**2 
        evb%grad_coupl(m,k) = grad_t1*t2*t3+t1*grad_t2*t3+t1*t2*grad_t3
      End If       

    End Do
  End Do 


  End Subroutine grad_evb_couplings

  
   Subroutine evb_force(evb,flow,config)
!> Compute EVB forces
 
  Type(evb_type), Intent(InOut) :: evb
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type( flow_type ), Intent(In   ) :: flow


  Integer( Kind = wi ) :: m,k     ! Indices for matrix elements
  Integer( Kind = wi ) :: i,j     ! Indices for atoms and/or coordinates

  Integer( Kind = wi ) :: natms

  Real( Kind = wp )    :: fdiff

! Maximum number of atoms
  natms=config(1)%natms 

   Do i=1, natms
     ! Copy forces from config%parts(i) to evb%force. x->1, y->2, y->3. We do this for each force field.
     Do m=1,flow%NUM_FF
       evb%force(1,m)=config(m)%parts(i)%fxx
       evb%force(2,m)=config(m)%parts(i)%fyy
       evb%force(3,m)=config(m)%parts(i)%fzz
     End Do
! For particle i, loop over the three coordinates
     Do j=1,3
       evb%force_matrix=0.0_wp
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Define matrix with gradients
       ! Diagonal elements
       Do m=1,flow%NUM_FF
         evb%force_matrix(m,m)=evb%force(j,m)
       End Do
       ! Off-diagonal elements     
       Do m=1,flow%NUM_FF
         Do k=m+1,flow%NUM_FF
         fdiff=evb%force(j,m)-evb%force(j,k)       
         ! Compute coupling forces for element (k,m) of the matrix, gradients of coupling terms are computed in 
           evb%force_matrix(m,k)=fdiff*evb%grad_coupl(m,k) 
           ! Make matrix symmetric
           evb%force_matrix(k,m)=evb%force_matrix(m,k)
         End Do
       End Do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Set the first column of evb%force to zero   
       evb%force(j,1)=0.0_wp

       ! Matrix multiplication (evb%psi)^T*evb%force_matrix*(evb%psi) to compute evb%force(j,1)
       ! Store EVB force in column 1 of evb%force
       Do m=1,flow%NUM_FF
         Do k=1,flow%NUM_FF
           evb%force(j,1)=evb%force(j,1)+evb%psi(m,1)*evb%force_matrix(m,k)*evb%psi(k,1)
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


!> Subroutine to compute EVB stress
  Subroutine evb_stress(evb,flow,config,stat)
 
  Type(evb_type), Intent(InOut) :: evb
  Type( configuration_type ), Intent(InOut) :: config(:)
  Type(stats_type), Intent(InOut) :: stat(:)
  Type( flow_type ), Intent(In   ) :: flow


  Integer( Kind = wi ) :: m,k     ! Indices for matrix elements
  Integer( Kind = wi ) :: i,j,i2  ! Indices for atoms and/or coordinates

  Real( Kind = wp )    :: invcell(1:9), det
  Real( Kind = wp )    :: cell(3,3), rcell(3,3)
  Real( Kind = wp )    :: strff(flow%NUM_FF,3,3)
   

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

! Here we build the evb%dE_dh matrix to compute the components of the EVB stress tensor 
! Loop over the 6 independent components of the stress tensor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Do i=1,3
    Do j=1,3
      ! Define matrix components
      ! -----------------------
      ! Diagonal elements
      Do m=1,flow%NUM_FF
        evb%stress_matrix(m,m)= 0.0_wp
        Do i2=1,3  
          evb%stress_matrix(m,m)= evb%stress_matrix(m,m)+strff(m,i,i2)*rcell(j,i2)
        End Do
      End Do
      ! Off-diagonal elements     
      Do m=1,flow%NUM_FF
        Do k=m+1,flow%NUM_FF
          evb%stress_matrix(m,k)=0.0_wp
          Do i2=1,3
            evb%stress_matrix(m,k)=evb%stress_matrix(m,k)+(strff(m,i,i2)-strff(k,i,i2))*rcell(j,i2)
          End Do
           ! Compute coupling stress for element (k,m) of the matrix
           evb%stress_matrix(m,k) = evb%stress_matrix(m,k)*evb%grad_coupl(m,k)
           !  Make matrix symmetric
          evb%stress_matrix(k,m)=evb%stress_matrix(m,k)
        End Do 
      End Do
      !------------------------
      evb%dE_dh(i,j)=0.0_wp
  
      ! Matrix multiplication (evb%psi)^T*(evb%stress_matrix)*(evb%psi) to compute stat(1)%stress(j)
      ! Store EVB stress in stress component of force-field 1
      Do m=1,flow%NUM_FF
        Do k=1,flow%NUM_FF
          evb%dE_dh(i,j)=evb%dE_dh(i,j)+evb%psi(m,1)*evb%stress_matrix(m,k)*evb%psi(k,1)
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
        evb%stress(i,j)=evb%stress(i,j)+evb%dE_dh(i,i2)*cell(j,i2)
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


!> Population analysis
  Subroutine evb_population(evb,flow,files,comm)
  Type(evb_type)   , Intent(InOut)   :: evb
  Type( flow_type ), Intent(In   )   :: flow
  Type( file_type ), Intent( InOut ) :: files(:)
  Type( comms_type), Intent( InOut ) :: comm
  Integer( Kind = wi ) :: ff
  Logical              :: l_tmp

    ! open EVB file file and write header
    If (evb%newjob .and. comm%idnode == 0) Then

      evb%newjob = .false.
      l_tmp=.false.
      
      ! If the flow%restart_key = RESTART_KEY_OLD is the file old (does it exist)?
      If (flow%restart_key == RESTART_KEY_OLD)Then
        Inquire(File=files(FILE_POPEVB)%filename, Exist=l_tmp)
      End If

      If (.not.l_tmp) Then
        Open(Newunit=files(FILE_POPEVB)%unit_no,File=files(FILE_POPEVB)%filename, Status='Replace')
        write(files(FILE_POPEVB)%unit_no,'(a,17x,a)') '# Time (ps)','Weights of each FFs in the total EVB state (FF 1, 2, etc)'
      Else
        Open(Newunit=files(FILE_POPEVB)%unit_no,File=files(FILE_POPEVB)%filename, Position='append')      
      End If  
  
    End If

  ! Only print after equilibration
  If(flow%step>flow%equil_steps+1)Then
     If (comm%idnode == 0)Then
        write(files(FILE_POPEVB)%unit_no,*) flow%time, (evb%psi(ff,1)**2, ff=1,flow%NUM_FF) 
     End If
  End If

  End subroutine evb_population



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

  Subroutine evb_merge_stochastic(flow,config,stat,rigid,thermo,cshell,cons,pmf)
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This subroutine copy fields of variable types only when regauss is applied
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

End Module evb        
