module angular_distribution
  Use constants, Only : pi,nchadf
  Use configuration, Only : configuration_type
  Use kinds, Only : wp
  Use site, Only : site_type
  Use neighbours, Only : neighbours_type
  Use coord, Only : coord_type
  Implicit none
  type, public :: adf_type
        real(wp) :: rij(1:10),rik,rjk
        Integer, allocatable :: astat(:,:),coordlist(:,:)
!        real(wp), allocatable :: rij(1:10)

  end type adf_type

contains




  subroutine adf_calculate(config,sites,crd,adf)
    Type(configuration_type), Intent(In) :: config
    Type(coord_type), Intent(In) :: crd
    Type(adf_type), Intent(InOUT) :: adf
    Type(site_type), Intent(In) :: sites
    integer :: i,ii,iii,j,jj,jjj,k,kk,kkk
    real :: costheta,temptheta
    
    Open(Unit=nchadf, File='ADFDAT', Form='formatted')
   allocate(adf%astat(-1:180,1:2*crd%ncoordpairs))
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

    do j= 1,crd%coordlist(0,i)
      adf%rij(j)=(config%parts(i)%xxx-config%parts(crd%adfcoordlist(j,i))%xxx)**2 &
      +(config%parts(i)%yyy-config%parts(crd%adfcoordlist(j,i))%yyy)**2 &  
      + (config%parts(i)%zzz-config%parts(crd%adfcoordlist(j,i))%zzz)**2
     End Do
 

     Do j=1, crd%coordlist(0,i)-1
      Do jj=1+j, crd%coordlist(0,i)
         if(config%ltype(j).eq.config%ltype(jj))then
          do ii = 1,2*crd%ncoordpairs
         if(adf%astat(-1,ii)==config%ltype(crd%coordlist(j,i)) .and. adf%astat(0,ii)==config%ltype(i))then       
           
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

    Do i=1,2*crd%ncoordpairs
    write(nchadf,*)sites%unique_atom(adf%astat(-1,i)),'-',sites%unique_atom(adf%astat(0,i))&
            ,'-',sites%unique_atom(adf%astat(-1,i))
    do ii= 1,180
    write(nchadf,*)real(ii)-0.5,adf%astat(ii,i)
    End Do
    End DO
 
    end subroutine adf_calculate

    end module angular_distribution



































