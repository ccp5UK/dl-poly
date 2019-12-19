module angular_distribution
  Use constants, Only : pi,nchadf
  Use configuration, Only : configuration_type
  Use kinds, Only : wp
  Use site, Only : site_type
  Use neighbours, Only : neighbours_type
  Use coord, Only : coord_type
  Use flow_control, Only : flow_type
  Use comms, Only : comms_type,grecv,gsend
  Implicit none
  type, public :: adf_type
        real(wp) :: rij(1:10),rik,rjk
        Integer, allocatable :: astat(:,:),coordlist(:,:)
        Integer :: adfinterval
        Logical :: adfon
!        real(wp), allocatable :: rij(1:10)

  end type adf_type

contains




  subroutine adf_calculate(config,sites,flow,crd,adf,comm)
    Type(configuration_type), Intent(In) :: config
    Type(coord_type), Intent(In) :: crd
    Type(adf_type), Intent(InOUT) :: adf
    Type(site_type), Intent(In) :: sites
    Type(flow_type), Intent(In) :: flow
    Type(comms_type), Intent(InOut) :: comm
    integer :: i,ii,iii,j,jj,jjj,k,kk,kkk,nab
    integer, allocatable :: adfbuff(:)
    real :: costheta,temptheta
    logical :: itsopen
    if (adf%adfon .Eqv. .False.)Return
    If(crd%coordon .Eqv. .False.)Return
    If(mod(flow%step,adf%adfinterval).NE.0)Return
    Open(Unit=nchadf, File='ADFDAT', Form='formatted')
    If(flow%step.eq.0)then
   allocate(adf%astat(-1:180,1:2*crd%ncoordpairs))
    endif
    adf%astat(:,:)=0
    do i=1,crd%ncoordpairs
       adf%astat(-1,(2*i)-1)=crd%ltype(i,1)
       adf%astat(0,(2*i)-1)=crd%ltype(i,2)
       adf%astat(-1,(2*i))=crd%ltype(i,2)
       adf%astat(0,(2*i))=crd%ltype(i,1)
    enddo
    
    Do i= 1,config%natms
     adf%rij(:)=0_wp
     adf%rjk=0_wp

    do j= 1,crd%adfcoordlist(0,i)
      adf%rij(j)=(config%parts(i)%xxx-config%parts(crd%adfcoordlist(j,i))%xxx)**2 &
      +(config%parts(i)%yyy-config%parts(crd%adfcoordlist(j,i))%yyy)**2 &  
      + (config%parts(i)%zzz-config%parts(crd%adfcoordlist(j,i))%zzz)**2
     End Do
 

     Do j=1, crd%adfcoordlist(0,i)-1
      Do jj=1+j, crd%adfcoordlist(0,i)
         if(config%ltype(j).eq.config%ltype(jj))then
          do ii = 1,2*crd%ncoordpairs
         if(adf%astat(-1,ii)==config%ltype(crd%adfcoordlist(j,i)) .and. adf%astat(0,ii)==config%ltype(i))then       
         adf%rjk=(config%parts(crd%adfcoordlist(j,i))%xxx - config%parts(crd%adfcoordlist(jj,i))%xxx)**2 &
                + (config%parts(crd%adfcoordlist(j,i))%yyy - config%parts(crd%adfcoordlist(jj,i))%yyy)**2 &
                + (config%parts(crd%adfcoordlist(j,i))%zzz - config%parts(crd%adfcoordlist(jj,i))%zzz)**2
        costheta=((adf%rij(j) + adf%rij(jj) - adf%rjk)/(2*sqrt(adf%rij(j))*sqrt(adf%rij(jj)) ))
        temptheta=ACOS(costheta)*(180/pi)

            do iii=1,180
             if(temptheta.ge.(iii-1) .and. temptheta.lt.(iii))then

             adf%astat(iii,ii)=adf%astat(iii,ii)+1
             End If
            End Do
    End If            
    End Do
       End If
     End do    
      End do   
   


    End do
    print*,adf%astat(170,1),adf%astat(180,2)
nab=0
   do i=1,2*crd%ncoordpairs
   nab=nab+182
   enddo

   

   If(comm%idnode==0)then
      inquire(unit=nchadf, opened=itsopen)
        if(.not. itsopen)then
         Open(Unit=nchadf, File='ADFDAT', Form='formatted')
        endif
       
        do j=1,comm%mxnode-1
           call grecv(comm,nab,j,j)
            if(nab>0)then
            allocate(adfbuff(nab))
            call grecv(comm,adfbuff,j,j)
            do ii=1,2*crd%ncoordpairs
              if(adf%astat(-1,ii)/=adfbuff(1+182*(ii-1)) .or. adf%astat(0,ii)/=adfbuff(2+182*(ii-1)))then
                 write(*,*) 'ERROR: adf pairs do not match in MPI'
              endif
            do kk=1,180
               adf%astat(kk,ii)=adf%astat(kk,ii)+adfbuff(2+kk+182*(ii-1))
            enddo
            enddo
            deallocate(adfbuff)
            end if
         enddo
         
!        do ii=1,2*crd%ncoordpairs
!        do kk=1,180
!        adf%astatavg(kk,ii)=adf%astatavg(kk,ii)+adf%astat(kk,ii)
!        enddo
!        enddo


    write(nchadf,'(A29,I10,F20.6)')"Angular distribution function",flow%step,flow%time
    Do i=1,2*crd%ncoordpairs
    write(nchadf,*)trim(sites%unique_atom(adf%astat(-1,i))),'-',trim(sites%unique_atom(adf%astat(0,i)))&
            ,'-',trim(sites%unique_atom(adf%astat(-1,i)))
    do ii= 1,180
    write(nchadf,*)real(ii)-0.5,adf%astat(ii,i)
    End Do
    End DO
   



   else
    allocate(adfbuff(nab))
    k=0
    do i=1,2*crd%ncoordpairs
       k=k+1
       adfbuff(k)=adf%astat(-1,i)    
       k=k+1
       adfbuff(k)=adf%astat(0,i)
       do ii=1,180
          k=k+1
       adfbuff(K)=adf%astat(ii,i)
       enddo
       enddo
      Call gsend(comm,nab,0,comm%idnode)
         if(nab>0)then
            call gsend(comm,adfbuff,0,comm%idnode)
         endif
    deallocate(adfbuff)
     endif 
    end subroutine adf_calculate

    end module angular_distribution



































