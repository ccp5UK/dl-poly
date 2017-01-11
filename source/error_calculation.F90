module block_averages_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for error computations
!
! copyright - daresbury laboratory
! author   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use kinds_f90
  Use rdf_module, Only: tmp_rdf, block_averages, block_size, block_number, lstrdf, ncfrdf, rdf, ntprdf, tmp_rdf_sync, num_blocks
  Use kinds_f90
  Use config_module
  Use parse_module
  Use io_module, Only : io_close
Use setup_module, Only : engunit, boltz, mxgrdf, delr_max, fourpi, mxlist, mxrdf, nrdfdt,nrite
Use site_module, Only: ntpatm, dens, unqatm, numtyp
Use comms_module, Only: mxnode, gsum, idnode
  Implicit None
  Public :: calculate_errors, calculate_errors_jackknife
  Private :: calculate_block
  
  Contains
  
  Subroutine calculate_block(temp, rcut)
  
  Implicit None
  
  Real( Kind = wp ), Intent(in) :: temp, rcut
  Integer :: i, loopend, j, ia, ib, ngrid, kk, k, limit, jj
Real( Kind = wp ), Dimension( 1:mxlist ) :: rrt, xxt, yyt, zzt
Real( Kind = wp ) :: kT2engo, delr, rdlr, dgrid, pdfzero, factor1, rrr,dvol,gofr,gofr1
Logical :: zero
  
   kT2engo = boltz*temp/engunit
! grid interval for rdf tables
  delr = rcut/Real(mxgrdf,wp)
  rdlr = 1.0_wp/delr
! resampling grid and grid interval for rdf tables
  ngrid = Max(Nint(rcut/delr_max),mxgrdf)
  dgrid = rcut/Real(ngrid,wp)
  pdfzero = 1.0e-9_wp
  Do ia=1,ntpatm
     Do ib=ia,ntpatm
! number of the interaction by its rdf key
        kk=lstrdf(ib*(ib-1)/2+ia)
! only for valid interactions specified for a look up
! global sum of data on all nodes
! normalisation factor
           factor1=volm*dens(ia)*dens(ib)*Real(ncfrdf,wp)
           If (ia == ib) factor1=factor1*0.5_wp*(1.0_wp-1.0_wp/numtyp(ia))
! loop over distances
           zero=.true.
           Do i=1,mxgrdf
              If (zero .and. i < (mxgrdf-3)) zero=(tmp_rdf(i+2,kk, block_number) <= 0.0_wp) 
              gofr= tmp_rdf(i,kk, block_number)/factor1
              rrr = (Real(i,wp)-0.5_wp)*delr
              dvol= fourpi*delr*(rrr**2+delr**2/12.0_wp)
              gofr= gofr/dvol
! zero it if < pdfzero
              If (gofr < pdfzero) Then
                 gofr1 = 0.0_wp
              Else
                 gofr1 = gofr
              End If
! Store information to compute block average
              block_averages(ia, ib,i, block_number) = block_averages(ia,ib,i, block_number) + gofr1 !!* i_blocksize
           End Do
     End Do
  End Do
  
  
  End Subroutine calculate_block
  
  Subroutine calculate_errors(temp, rcut, num_steps)
  
  Implicit None
  
  Real( Kind = wp ), Intent(in) :: temp, rcut
  Integer, Intent(in) :: num_steps
  Integer :: nr_blocks, i, j,k ,l, ierr, ierr2, ierr3, a, b, ia, ib, kk
  Real( Kind = wp ) :: test1, delr
  Real( kind=wp ), Dimension( :, : , :), Allocatable :: averages, errors
  Real( kind = wp) :: i_nr_blocks, s
  
    test1 = 0.0_wp
  block_number = 1
  
  if(mxnode > 1 .and. .not. tmp_rdf_sync) then
    do i=1, num_blocks+1
      CALL gsum(tmp_rdf(:,:,i))
    end do
    tmp_rdf_sync = .TRUE.
  end if
    Allocate(averages(ntpatm,ntpatm, mxgrdf), stat = ierr2)
  Allocate(errors(ntpatm,ntpatm, mxgrdf), stat = ierr3)
  if(ierr > 0 .or. ierr2 > 0 .or. ierr3 > 0) then
    call error(1084)
  end if  
  averages = 0.0_wp
  errors = 0.0_wp
  !!Compute the rdf for each of the blocks
  do block_number=1, num_blocks+1
    call calculate_block(temp, rcut)
  end do 
nr_blocks = num_blocks+1
!Output blocks!
  If(idnode == 0) Then
    Open(Unit=nrdfdt, File='RDFBLOCKS', Status='replace')
 
    do j=1, ntpatm
      do k=j, ntpatm
        kk=lstrdf(k*(k-1)/2+j)
        if(kk > 0 .and. kk < ntprdf) Then
            Write(nrite,"(/,' g(r)  :',2(1x,a8),/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(j),unqatm(k)
            Write(nrdfdt,'(2a8)') unqatm(j),unqatm(k)
          do i = 1,nr_blocks
            do ia = 1,mxgrdf
              write(nrdfdt, "(1p,2e14.6)") ((Real(ia,wp)-0.5_wp)*delr),block_averages(j,k,ia,i)
            end do
              write(nrdfdt, '()')
          end do 
        End If
     end do
    end do
  Close(Unit=nrdfdt)
  End if
  
  i_nr_blocks = 1.0_wp/real(nr_blocks, wp)
  do k=1, nr_blocks
    do l=1, mxgrdf
      do j=1, ntpatm
        do i=1, ntpatm
          averages(i,j,l) = averages(i,j,l) + block_averages(i,j,l,k) * i_nr_blocks
        end do      
      end do
    end do
  end do
  
  i_nr_blocks = 1.0_wp/real(nr_blocks *(nr_blocks-1), wp)
  do i=1, nr_blocks
    do k=1, ntpatm
      do j=1, ntpatm
        do l=1, mxgrdf
        errors(j,k,l) = errors(j,k,l) + ( (block_averages(j,k,l,i) - averages(j,k,l))**2 * i_nr_blocks )
        end do
      end do
    end do
  end do
  
do l=1, mxgrdf
  do j=1, ntpatm
    do i = 1, ntpatm
      averages(i,j,l) = averages(i,j,l) *real(nr_blocks,wp)
    end do
  end do
end do
  !output errors
    If (idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf
  
  delr = rcut/Real(mxgrdf,wp)
  do j =1, ntpatm
    do k = j, ntpatm
        kk=lstrdf(k*(k-1)/2+j)
        If (kk > 0 .and. kk <= ntprdf) Then
          !If (idnode == 0) Then
            Write(nrite,"(/,' g(r)  :',2(1x,a8),/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(j),unqatm(k)
            Write(nrdfdt,'(2a8)') unqatm(j),unqatm(k)
          !End If
          Do i=1,mxgrdf
            Write(nrdfdt,"(1p,2e14.6,2e14.6)") ((Real(i,wp)-0.5_wp)*delr),averages(j,k,i),errors(j,k,i)
          
          End Do
        End If
    End Do
  End Do
  Close(Unit=nrdfdt)
  End If  
  Deallocate(averages, errors)
  End Subroutine calculate_errors
  
  Subroutine calculate_errors_jackknife(temp, rcut, num_steps)
   
   Implicit None
  
  Real( Kind = wp ), Intent(in) :: temp, rcut
  Integer, Intent(in) :: num_steps
  Integer :: nr_blocks, i, j,k ,l, ierr, ierr2, ierr3, a, b, kk
  Real( Kind = wp ) :: test1
  Real( kind=wp ), Dimension( :, : , :), Allocatable :: averages, errors
  Real(Kind = wp) :: i_nr_blocks, delr, s
  
    test1 = 0.0_wp
  block_number = 1
  print *, mxnode
  if(mxnode > 1 .and. .not. tmp_rdf_sync) then
    do i=1, num_blocks+1
      CALL gsum(tmp_rdf(:,:,i))
    end do
    tmp_rdf_sync = .TRUE.
  end if
  
    Allocate(averages(ntpatm,ntpatm, mxgrdf), stat = ierr2)
  Allocate(errors(ntpatm,ntpatm, mxgrdf), stat = ierr3)
  if(ierr > 0 .or. ierr2 > 0 .or. ierr3 > 0) then
    call error(1084)
  end if  
  averages = 0.0_wp
  errors = 0.0_wp
  block_averages =0.0_wp
    !!Compute the rdf for each of the blocks
  do block_number=1,num_blocks+1
    call calculate_block(temp, rcut)
  end do 
nr_blocks = num_blocks+1
    i_nr_blocks = 1.0_wp/real(nr_blocks, wp)
  
  do k=1, nr_blocks
    do l=1, mxgrdf
      do j=1, ntpatm
        do i=1, ntpatm
          averages(i,j,l) = averages(i,j,l) + block_averages(i,j,l,k) !* i_nr_blocks
        end do      
      end do
    end do
  end do
  
  i_nr_blocks = 1.0_wp / real(nr_blocks-1, wp)
  !Create jackknife bins
  do k=1, nr_blocks
    do l=1, mxgrdf
      do j=1, ntpatm
        do i=1, ntpatm
          block_averages(i,j,l,k) = (averages(i,j,l) - block_averages(i,j,l,k)) * i_nr_blocks
        end do      
      end do
    end do
  end do
  
  !Average
  i_nr_blocks = 1.0_wp/real(nr_blocks,wp)
  do l=1, mxgrdf
    do j=1, ntpatm
      do i=1, ntpatm
        averages(i,j,l) = averages(i,j,l) * i_nr_blocks
      end do
    end do
  end do
  
  !Errors
  !Compute the errors
  i_nr_blocks = real((nr_blocks-1), wp) / real(nr_blocks, wp)
  do i=1, nr_blocks
    do k=1, ntpatm
      do j=1, ntpatm
        do l=1, mxgrdf
        errors(j,k,l) = errors(j,k,l) + ( (block_averages(j,k,l,i) - averages(j,k,l))**2 * i_nr_blocks )
        end do
      end do
    end do
  end do
  
do l=1, mxgrdf
  do j=1, ntpatm
    do i = 1, ntpatm
      averages(i,j,l) = averages(i,j,l)*real(nr_blocks,wp)
    end do
  end do
end do
  !output errors
    If (idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf
  
  delr = rcut/Real(mxgrdf,wp)
  do j =1, ntpatm
    do k = j, ntpatm
        kk=lstrdf(k*(k-1)/2+j)
        If (kk > 0 .and. kk <= ntprdf) Then
            Write(nrite,"(/,' g(r)  :',2(1x,a8),/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") unqatm(j),unqatm(k)
            Write(nrdfdt,'(2a8)') unqatm(j),unqatm(k)
          Do i=1,mxgrdf
            Write(nrdfdt,"(1p,2e14.6,2e14.6)") ((Real(i,wp)-0.5_wp)*delr),averages(j,k,i),errors(j,k,i)
          
          End Do
        End If
    End Do
  End Do
  Close(Unit=nrdfdt)
  End If
  End Subroutine calculate_errors_jackknife
end module block_averages_module
