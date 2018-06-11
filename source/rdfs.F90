Module rdfs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring RDF property variables and arrays
! including USR (umbrella sampling restraint) RDF
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2016
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use site, Only: ntpatm
  Use configuration, Only : natms,ltg,ltype
  Use comms,  Only : comms_type,gsum
  Use setup,  Only : fourpi,boltz,delr_max,nrdfdt,npdfdt,npdgdt, &
                            mxgrdf,engunit,zero_plus,mxrdf,mxgusr
  Use site,   Only : ntpatm,unqatm,numtyp,dens
  Use configuration, Only : cfgname,volm
  Use parse
  Use io
  Use errors_warnings, Only : error,info
  Use neighbours, Only : neighbours_type

  Implicit None

  Public  :: rdf_compute, calculate_errors, calculate_errors_jackknife
  Private :: calculate_block

  Integer,                        Save :: ncfrdf = 0 , &
                                          ntprdf = 0

  Integer,                        Save :: ncfusr = 0

  Real( Kind = wp ),              Save :: rusr   = 0.0_wp ! USR RDF cutoff

  Integer,           Allocatable, Save :: lstrdf(:)

  Real( Kind = wp ), Allocatable, Save :: rdf(:,:),usr(:)

  Real( Kind = wp ), Allocatable, Save :: block_averages(:,:,:,:)
  Integer, Parameter                   :: num_blocks = 25
  Integer, Save                        :: block_size
  Integer,                        Save :: block_number = 1
  Real( Kind = wp ), Allocatable, Save :: tmp_rdf(:,:,:)
  Logical,                        Save :: tmp_rdf_sync = .FALSE.
  Logical,                        Save :: l_errors_block = .FALSE., l_errors_jack = .FALSE.

  Public :: allocate_rdf_arrays, allocate_block_average_array

Contains

  Subroutine allocate_rdf_arrays()


    Integer, Dimension( 1:3 ) :: fail

    fail = 0

    Allocate (lstrdf(1:mxrdf),       Stat = fail(1))
    Allocate (rdf(1:mxgrdf,1:mxrdf), Stat = fail(2))
    Allocate (usr(1:mxgusr),         Stat = fail(3))

    If (Any(fail > 0)) Call error(1016)

    lstrdf = 0

    rdf = 0.0_wp ; usr = 0.0_wp

  End Subroutine allocate_rdf_arrays

Subroutine allocate_block_average_array(nstrun)

  Integer, Intent( In ) :: nstrun
  Integer :: temp1, temp2
  
  Integer, Dimension( 1:2 ) :: fail
  block_size = nstrun/(num_blocks-1)
  if(block_size < 2) then
    block_size = 2
  endif

  temp1 = mxrdf + 16-Mod(mxrdf,16)
  temp2 = mxgrdf + 16-Mod(mxgrdf,16)
  Allocate(block_averages(1:ntpatm,1:ntpatm,1:mxgrdf,1:num_blocks+1), Stat = fail(1))
  Allocate(tmp_rdf( 1:temp2,1:temp1, 1:num_blocks+1 ), Stat = fail(2))

  If (Any(fail > 0)) Call error(1016)
  block_averages = 0.0_wp
  tmp_rdf = 0.0_wp

  End Subroutine allocate_block_average_array


  Subroutine rdf_collect(iatm,rrt,neigh)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for radial
! distribution functions
!
! Note: to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - t.forester march 1994
! amended   - i.t.todorov november 2014
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Type( neighbours_type), Intent( In    ) :: neigh
  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: rrt

  Integer                 :: idi,jatm,ai,aj,keyrdf,kk,ll,m
  Real( Kind = wp )       :: rdelr,rrr

! set cutoff condition for pair forces and grid interval for rdf tables

  rdelr= Real(mxgrdf,wp)/neigh%cutoff

! global identity and type of iatm

  idi=ltg(iatm)
  ai=ltype(iatm)

! start of primary loop for rdf accumulation

  Do m=1,neigh%list(0,iatm)

! atomic and type indices

     jatm=neigh%list(m,iatm)
     aj=ltype(jatm)

     If (jatm <= natms .or. idi < ltg(jatm)) Then

! rdf function indices

        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=lstrdf(keyrdf)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then

! apply truncation of potential

           rrr=rrt(m)

           If (rrr < neigh%cutoff) Then
              ll=Min(1+Int(rrr*rdelr),mxgrdf)

! accumulate correlation

              rdf(ll,kk) = rdf(ll,kk) + 1.0_wp
              If(l_errors_block .or. l_errors_jack) tmp_rdf(ll,kk,block_number) = tmp_rdf(ll,kk,block_number) + 1.0_wp

           End If

        End If

     End If

  End Do

End Subroutine rdf_collect

Subroutine rdf_compute(lpana,rcut,temp,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating radial distribution functions
! from accumulated data
!
! copyright - daresbury laboratory
! author    - t.forester & i.t.todorov march 2016
! contrib   - a.v.brukhno january 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Logical          , Intent( In    ) :: lpana
  Real( Kind = wp ), Intent( In    ) :: rcut,temp
  Type(comms_type), Intent( InOut )  :: comm

  Logical           :: zero
  Integer           :: fail,ngrid,i,ia,ib,kk,ig,ll
  Real( Kind = wp ) :: kT2engo,delr,rdlr,dgrid,pdfzero,      &
                       factor1,rrr,dvol,gofr,gofr1,sum,sum1, &
                       fed0,fed,dfed,dfed0,tmp,fed1,fed2,dfed1,dfed2,coef,t1,t2

  Real( Kind = wp ), Allocatable :: dstdrdf(:,:)
  Real( Kind = wp ), Allocatable :: pmf(:),vir(:)

  Character ( Len = 256 )  :: message,messages(2)

  If (lpana) Then
     fail = 0
     Allocate (dstdrdf(0:mxgrdf,1:ntprdf),pmf(0:mxgrdf+2),vir(0:mxgrdf+2), Stat = fail)
     If (fail > 0) Then
        Write(message,'(a)') 'rdf_compute - allocation failure'
        Call error(0,message)
     End If
  End If

! conversion: internal units -> in/out units (kJ/mol, kcal/mol, eV etc)

  kT2engo = boltz*temp/engunit

! grid interval for rdf tables

  delr = rcut/Real(mxgrdf,wp)
  rdlr = 1.0_wp/delr

! resampling grid and grid interval for rdf tables

  ngrid = Max(Nint(rcut/delr_max),mxgrdf)
  dgrid = rcut/Real(ngrid,wp)

  Write(messages(1),'(a)') 'radial distribution functions:'
  Write(messages(2),'(2x,a,i8,a)') 'calculated using ', ncfrdf, ' configurations'
  Call info(messages,2,.true.)

! open RDF file and Write headers

  If (comm%idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf
  End If

! the lower bound to nullify the nearly-zero histogram (PDF) values

  pdfzero = 1.0e-6_wp

! for all possible unique type-to-type pairs

  Do ia=1,ntpatm
     Do ib=ia,ntpatm

! number of the interaction by its rdf key

        kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then
           Write(messages(1),'(2x,a,2(1x,a8))') 'g(r): ',unqatm(ia),unqatm(ib)
           Write(messages(2),'(8x,a1,6x,a4,9x,a4)') 'r','g(r)','n(r)'
           Call info(messages,2,.true.)
           If (comm%idnode == 0) Then
              Write(nrdfdt,'(2a8)') unqatm(ia),unqatm(ib)
           End If

! global sum of data on all nodes

           Call gsum(comm,rdf(1:mxgrdf,kk))

! normalisation factor

           factor1=volm*dens(ia)*dens(ib)*Real(ncfrdf,wp)
           If (ia == ib) factor1=factor1*0.5_wp*(1.0_wp-1.0_wp/numtyp(ia))

! running integration of rdf

           sum=0.0_wp

! loop over distances

           zero=.true.
           Do i=1,mxgrdf
              If (zero .and. i < (mxgrdf-3)) zero=(rdf(i+2,kk) <= 0.0_wp)

              gofr= rdf(i,kk)/factor1
              sum = sum + gofr*dens(ib)

              rrr = (Real(i,wp)-0.5_wp)*delr
              dvol= fourpi*delr*(rrr**2+delr**2/12.0_wp)
              gofr= gofr/dvol

! zero it if < pdfzero

              If (gofr < pdfzero) Then
                 gofr1 = 0.0_wp
              Else
                 gofr1 = gofr
              End If

              If (sum < pdfzero) Then
                 sum1 = 0.0_wp
              Else
                 sum1 = sum
              End If

! print out information

              If (.not.zero) Then
                Write(message,'(f10.4,1p,2e14.6)') rrr,gofr1,sum1
                Call info(message,.true.)
              End If
              If (comm%idnode == 0) Then
                 Write(nrdfdt,"(1p,2e14.6)") rrr,gofr
              End If

! We use the non-normalised tail-truncated RDF version,
! rdf...1 (not pdf...) in order to exclude the nearly-zero
! rdf... noise in PMF, otherwise the PMF = -ln(PDF)
! would have poorly-defined noisy "borders/walls"

                If (lpana) dstdrdf(i,kk) = gofr1 ! RDFs density
           End Do
        Else
           If (lpana) dstdrdf(:,kk) = 0.0_wp ! RDFs density
        End If

     End Do
  End Do

  If (comm%idnode == 0) Close(Unit=nrdfdt)

! Only when PDF analysis is requested

  If (lpana) Then

! open PDF files and write headers

     If (comm%idnode == 0) Then
        Open(Unit=npdgdt, File='VDWPMF', Status='replace')
        Write(npdgdt,'(a)') '# '//cfgname
        Write(npdgdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',delr*Real(mxgrdf,wp),mxgrdf,delr,ntprdf, &
             '   conversion factor(kT -> energy units) =',kT2engo

        Open(Unit=npdfdt, File='VDWTAB', Status='replace')
        Write(npdfdt,'(a)') '# '//cfgname
        Write(npdfdt,'(a,f12.5,i10,f12.5,i10,a,e15.7)') '# ',dgrid*Real(ngrid,wp),ngrid,dgrid,ntprdf, &
          '   conversion factor(kT -> energy units) =',kT2engo
     End If

! loop over all valid RDFs

     Do ia=1,ntpatm
        Do ib=ia,ntpatm

! number of the interaction by its rdf key

           kk=lstrdf(ib*(ib-1)/2+ia)

! only for valid interactions specified for a look up

           If (kk > 0 .and. kk <= ntprdf) Then
              If (comm%idnode == 0) Then
                 Write(npdgdt,'(/,a2,2a8)') '# ',unqatm(ia),unqatm(ib)
                 Write(npdfdt,'(/,a2,2a8)') '# ',unqatm(ia),unqatm(ib)
              End If

! Smoothen and get derivatives

! RDFs -> 1 at long distances, so we do not shift the PMFs
! but instead put a cap over those, -Log(pdfzero) - hence,
! the upper bound for the PMF due to the zero RDF values

              fed0  = -Log(pdfzero)
              dfed0 = 10.0_wp
              dfed  = 10.0_wp

              Do ig=1,mxgrdf
                 tmp = Real(ig,wp)-0.5_wp
                 rrr = tmp*delr

                 If (dstdrdf(ig,kk) > zero_plus) Then
                    fed = -Log(dstdrdf(ig,kk)+pdfzero) !-fed0
                    If (fed0 <= zero_plus) Then
                       fed0 = fed
                       fed  = fed0
!                       fed = 0.0_wp
                    End If

                    If (ig < mxgrdf-1) Then
                       If (dstdrdf(ig+1,kk) <= zero_plus .and. dstdrdf(ig+2,kk) > zero_plus) &
                          dstdrdf(ig+1,kk) = 0.5_wp*(dstdrdf(ig,kk)+dstdrdf(ig+2,kk))
                    End If
                 Else
                    fed = fed0
!                    fed = 0.0_wp
                 End If

                 If      (ig == 1) Then
                    If      (dstdrdf(ig,kk) > zero_plus .and. dstdrdf(ig+1,kk) > zero_plus) Then
                       dfed = Log(dstdrdf(ig+1,kk)/dstdrdf(ig,kk))
                    Else If (dfed > 0.0_wp) Then
                       dfed = dfed0
                    Else
                       dfed =-dfed0
                    End If
                 Else If (ig == mxgrdf) Then
                    If      (dstdrdf(ig,kk) > zero_plus .and. dstdrdf(ig-1,kk) > zero_plus) Then
                       dfed = Log(dstdrdf(ig,kk)/dstdrdf(ig-1,kk))
                    Else If (dfed > 0.0_wp) Then
                       dfed = dfed0
                    Else
                       dfed =-dfed0
                    End If
                 Else If (dstdrdf(ig-1,kk) > zero_plus) Then
                    If (dstdrdf(ig+1,kk) > zero_plus) Then
                       dfed = 0.5_wp*(Log(dstdrdf(ig+1,kk)/dstdrdf(ig-1,kk)))
                    Else
                       dfed = 0.5_wp*Log(dstdrdf(ig-1,kk))
                    End If
                 Else If (dstdrdf(ig+1,kk) > zero_plus) Then
                    dfed =-0.5_wp*Log(dstdrdf(ig+1,kk))
                 Else If (dfed > 0.0_wp) Then
                    dfed = dfed0
                 Else
                    dfed =-dfed0
                 End If

                 pmf(ig) = fed
                 vir(ig) = dfed

! Print

                 If (comm%idnode == 0) &
                      Write(npdgdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*tmp
              End Do

! Define edges

              pmf(0)        = 2.0_wp*pmf(1)       -pmf(2)
              vir(0)        = 2.0_wp*vir(1)       -vir(2)
              pmf(mxgrdf+1) = 2.0_wp*pmf(mxgrdf)  -pmf(mxgrdf-1)
              vir(mxgrdf+1) = 2.0_wp*vir(mxgrdf)  -vir(mxgrdf-1)
              pmf(mxgrdf+2) = 2.0_wp*pmf(mxgrdf+1)-pmf(mxgrdf)
              vir(mxgrdf+2) = 2.0_wp*vir(mxgrdf+1)-vir(mxgrdf)

! resample using 3pt interpolation

              Do ig=1,ngrid
                 rrr = Real(ig,wp)*dgrid
                 ll = Int(rrr*rdlr)

! +0.5_wp due to half-a-bin shift in the original (bin-centered) grid

                 coef = rrr*rdlr-Real(ll,wp)+0.5_wp

                 fed0 = pmf(ll)
                 fed1 = pmf(ll+1)
                 fed2 = pmf(ll+2)

                 t1 = fed0 + (fed1 - fed0)*coef
                 t2 = fed1 + (fed2 - fed1)*(coef - 1.0_wp)

                 fed = t1 + (t2-t1)*coef*0.5_wp

                 dfed0 = vir(ll)
                 dfed1 = vir(ll+1)
                 dfed2 = vir(ll+2)

                 t1 = dfed0 + (dfed1 - dfed0)*coef
                 t2 = dfed1 + (dfed2 - dfed1)*(coef - 1.0_wp)

                 dfed = t1 + (t2-t1)*coef*0.5_wp

                 If (comm%idnode == 0) &
                    Write(npdfdt,"(f11.5,1p,2e14.6)") rrr,fed*kT2engo,dfed*kT2engo*rrr*rdlr
              End Do
           End If
        End Do
     End Do

     If (comm%idnode == 0) Then
        Close(Unit=npdgdt)
        Close(Unit=npdfdt)
     End If

     Deallocate (dstdrdf,pmf,vir, Stat = fail)
     If (fail > 0) Then
        Write(message,'(a)') 'rdf_compute - deallocation failure'
        Call error(0,message)
     End If
  End If

End Subroutine rdf_compute

Subroutine calculate_block(temp, rcut,neigh)

  Type( neighbours_type), Intent( In    ) :: neigh
  Real( Kind = wp ), Intent(in)            :: temp, rcut
  Real( Kind = wp ), Dimension( 1:neigh%max_list ) :: rrt, xxt, yyt, zzt
  Real( Kind = wp )                        :: kT2engo, delr, rdlr, dgrid, pdfzero, factor1, rrr,dvol,gofr,gofr1

  Integer :: i, loopend, j, ia, ib, ngrid, kk, k, limit, jj
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
           block_averages(ia, ib,i, block_number) = block_averages(ia,ib,i, block_number) + gofr1  
        End Do
     End Do
  End Do

End Subroutine calculate_block

Subroutine calculate_errors(temp, rcut, num_steps, neigh, comm)

  Real( Kind = wp ), Intent( In )                      :: temp, rcut
  Type( neighbours_type ), Intent( In    ) :: neigh
  Type(comms_type), Intent( InOut )                    :: comm

  Real( Kind = wp )                                    :: test1, delr
  Real( kind = wp ), Dimension( :, : , :), Allocatable :: averages, errors
  Real( kind = wp)                                     :: i_nr_blocks, s

  Integer, Intent(In) :: num_steps
  Integer             :: nr_blocks, i, j,k ,l, ierr, ierr2, ierr3, a, b, ia, ib, kk

  Character ( Len = 256 )  :: messages(2)

  test1 = 0.0_wp
  block_number = 1

  If(comm%mxnode > 1 .and. .not. tmp_rdf_sync) Then
     Do i=1, num_blocks+1
        Call gsum(comm,tmp_rdf(:,:,i))
     End Do
     tmp_rdf_sync = .TRUE.
  End If

  Allocate(averages(ntpatm,ntpatm, mxgrdf), stat = ierr2)
  Allocate(errors(ntpatm,ntpatm, mxgrdf), stat = ierr3)
  If(ierr > 0 .or. ierr2 > 0 .or. ierr3 > 0) Then
     Call error(1084)
  End If
  averages = 0.0_wp
  errors = 0.0_wp

!Compute the rdf for each of the blocks
  Do block_number=1, num_blocks+1
     Call calculate_block(temp, rcut,neigh)
  End Do
  nr_blocks = num_blocks+1

!Compute the errors.
  i_nr_blocks = 1.0_wp / Real(nr_blocks, wp)
  Do k=1, nr_blocks
     Do l=1, mxgrdf
        Do j=1, ntpatm
           Do i=1, ntpatm
              averages(i,j,l) = averages(i,j,l) + block_averages(i,j,l,k) * i_nr_blocks
           End Do
        End Do
     End Do
  End Do

  i_nr_blocks = 1.0_wp / Real(nr_blocks *(nr_blocks-1), wp)
  Do i=1, nr_blocks
     Do k=1, ntpatm
        Do j=1, ntpatm
           Do l=1, mxgrdf
               errors(j,k,l) = errors(j,k,l) + ( (block_averages(j,k,l,i) - averages(j,k,l))**2 * i_nr_blocks )
           End Do
        End Do
     End Do
  End Do

  Do l=1, mxgrdf
     Do j=1, ntpatm
        Do i = 1, ntpatm
           averages(i,j,l) = averages(i,j,l) * Real(nr_blocks,wp)
        End Do
    End Do
  End Do

!output errors
  If (comm%idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf

     delr = rcut/Real(mxgrdf,wp)
     Do j =1, ntpatm
        Do k = j, ntpatm
           kk=lstrdf(k*(k-1)/2+j)
           If (kk > 0 .and. kk <= ntprdf) Then
              Write(messages(1),'(2x,a,2(1x,a8))') 'g(r): ',unqatm(j),unqatm(k)
              Write(messages(2),'(8x,a1,6x,a4,9x,a4)') 'r','g(r)','n(r)'
              Call info(messages,2,.true.)
              Write(nrdfdt,'(2a8)') unqatm(j),unqatm(k)
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

Subroutine calculate_errors_jackknife(temp, rcut, num_steps,neigh,comm)

  Real( Kind = wp ), Intent(In)                        :: temp, rcut
  Type( neighbours_type ), Intent( In    ) :: neigh
  Type(comms_type), Intent( InOut )                    :: comm

  Real( Kind = wp )                                    :: test1
  Real( Kind = wp ), Dimension( :, : , :), Allocatable :: averages, errors
  Real(Kind = wp)                                      :: i_nr_blocks, delr, s

  Integer, Intent(in) :: num_steps
  Integer             :: nr_blocks, i, j,k ,l, ierr, ierr2, ierr3, a, b, kk

  Character( Len = 256 ) :: messages(2)

  test1 = 0.0_wp
  block_number = 1
  If(comm%mxnode > 1 .and. .not. tmp_rdf_sync) Then
     Do i=1, num_blocks+1
        Call gsum(comm,tmp_rdf(:,:,i))
     End Do
     tmp_rdf_sync = .TRUE.
  End If

  Allocate(averages(ntpatm,ntpatm, mxgrdf), stat = ierr2)
  Allocate(errors(ntpatm,ntpatm, mxgrdf), stat = ierr3)
  if(ierr > 0 .or. ierr2 > 0 .or. ierr3 > 0) then
     Call error(1084)
  end if
  averages = 0.0_wp
  errors = 0.0_wp
  block_averages =0.0_wp

!Compute the rdf for each of the blocks
  Do block_number=1,num_blocks+1
     Call calculate_block(temp, rcut,neigh)
  End Do
  nr_blocks = num_blocks+1
  i_nr_blocks = 1.0_wp / Real(nr_blocks, wp)

  Do k=1, nr_blocks
     Do l=1, mxgrdf
        Do j=1, ntpatm
           Do i=1, ntpatm
              averages(i,j,l) = averages(i,j,l) + block_averages(i,j,l,k) 
           End Do
        End Do
     End Do
  End Do


  i_nr_blocks = 1.0_wp / Real(nr_blocks-1, wp)
!Create jackknife bins
  Do k=1, nr_blocks
     Do l=1, mxgrdf
        Do j=1, ntpatm
           Do i=1, ntpatm
              block_averages(i,j,l,k) = (averages(i,j,l) - block_averages(i,j,l,k)) * i_nr_blocks
           End Do
        End Do
     End Do
  End Do

!Average
  i_nr_blocks = 1.0_wp / Real(nr_blocks,wp)
  Do l=1, mxgrdf
    Do j=1, ntpatm
      Do i=1, ntpatm
        averages(i,j,l) = averages(i,j,l) * i_nr_blocks
      End Do
    End Do
  End Do

!Errors
!Compute the errors
  i_nr_blocks = Real((nr_blocks-1), wp) / Real(nr_blocks, wp)
  Do i=1, nr_blocks
     Do k=1, ntpatm
        Do j=1, ntpatm
           Do l=1, mxgrdf
              errors(j,k,l) = errors(j,k,l) + ( (block_averages(j,k,l,i) - averages(j,k,l))**2 * i_nr_blocks )
           End Do
        End Do
     End Do
  End Do

  Do l=1, mxgrdf
     Do j=1, ntpatm
        Do i = 1, ntpatm
           averages(i,j,l) = averages(i,j,l)*Real(nr_blocks,wp)
        End Do
     End Do
  End Do

!output errors
  If (comm%idnode == 0) Then
     Open(Unit=nrdfdt, File='RDFDAT', Status='replace')
     Write(nrdfdt,'(a)') cfgname
     Write(nrdfdt,'(2i10)') ntprdf,mxgrdf

     delr = rcut/Real(mxgrdf,wp)
     Do j =1, ntpatm
        Do k = j, ntpatm
           kk=lstrdf(k*(k-1)/2+j)
           If (kk > 0 .and. kk <= ntprdf) Then
              Write(messages(1),'(2x,a,2(1x,a8))') 'g(r): ',unqatm(j),unqatm(k)
              Write(messages(2),'(8x,a1,6x,a4,9x,a4)') 'r','g(r)','n(r)'
              Call info(messages,2,.true.)
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

  Subroutine rdf_excl_collect(iatm,rrt,neigh)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for radial
! distribution functions of excluded pairs
!
! Note: to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( neighbours_type), Intent( In    ) :: neigh
  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: rrt

  Integer                 :: limit,idi,jatm,ai,aj,keyrdf,kk,ll,m
  Real( Kind = wp )       :: rdelr,rrr

! set cutoff condition for pair forces and grid interval for rdf tables

  rdelr= Real(mxgrdf,wp)/neigh%cutoff

! global identity and type of iatm

  idi=ltg(iatm)
  ai=ltype(iatm)

! Get neigh%list limit

  limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

! start of primary loop for rdf accumulation

  Do m=1,limit

! atomic and type indices

     jatm=neigh%list(neigh%list(0,iatm)+m,iatm)
     aj=ltype(jatm)

     If (jatm <= natms .or. idi < ltg(jatm)) Then

! rdf function indices

        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=lstrdf(keyrdf)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then

! apply truncation of potential

           rrr=rrt(m)

           If (rrr < neigh%cutoff) Then
              ll=Min(1+Int(rrr*rdelr),mxgrdf)

! accumulate correlation

              rdf(ll,kk) = rdf(ll,kk) + 1.0_wp
              If(l_errors_block .or. l_errors_jack) tmp_rdf(ll,kk,block_number) = tmp_rdf(ll,kk,block_number) + 1.0_wp
           End If

        End If

     End If

  End Do

End Subroutine rdf_excl_collect

Subroutine rdf_frzn_collect(iatm,rrt,neigh)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for radial
! distribution functions of frozen pairs
!
! Note: to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( neighbours_type), Intent( In    ) :: neigh
  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: rrt

  Integer                 :: limit,idi,jatm,ai,aj,keyrdf,kk,ll,m
  Real( Kind = wp )       :: rdelr,rrr

! set cutoff condition for pair forces and grid interval for rdf tables

  rdelr= Real(mxgrdf,wp)/neigh%cutoff

! global identity and type of iatm

  idi=ltg(iatm)
  ai=ltype(iatm)

! Get neigh%list limit

  limit=neigh%list(-2,iatm)-neigh%list(-1,iatm)

! start of primary loop for rdf accumulation

  Do m=1,limit

! atomic and type indices

     jatm=neigh%list(neigh%list(-1,iatm)+m,iatm)
     aj=ltype(jatm)

     If (jatm <= natms .or. idi < ltg(jatm)) Then

! rdf function indices

        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=lstrdf(keyrdf)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then

! apply truncation of potential

           rrr=rrt(m)

           If (rrr < neigh%cutoff) Then
              ll=Min(1+Int(rrr*rdelr),mxgrdf)

! accumulate correlation

              rdf(ll,kk) = rdf(ll,kk) + 1.0_wp
              If(l_errors_block .or. l_errors_jack) tmp_rdf(ll,kk,block_number) = tmp_rdf(ll,kk,block_number) + 1.0_wp
           End If

        End If

     End If

  End Do

End Subroutine rdf_frzn_collect


Subroutine usr_collect(rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for USR RDFs
!
! Note: to be used in external_field_apply
!
! copyright - daresbury laboratory
! author    - i.t.todorov & a.brukhno november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ), Intent( In    ) :: rrt

  Integer           :: ll
  Real( Kind = wp ) :: rdelr

! set cutoff condition for pair forces and grid interval for rdf tables

  rdelr= Real(mxgusr,wp)/rusr

  If (rrt < rusr) Then ! apply truncation of potential
     ll=Min(1+Int(rrt*rdelr),mxgusr)
     usr(ll) = usr(ll) + 1.0_wp ! accumulate correlation
     ncfusr = ncfusr + 1        ! Increment sample
  End If

End Subroutine usr_collect

Subroutine usr_compute(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating radial distribution function
! from accumulated data for umbrella sampled two COMs separation (ushr)
!
! to be used in exernal_field_apply & statistics_result
!
! copyright - daresbury laboratory
! author    - i.t.todorov & a.brukhno november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( comms_type ), Intent( InOut ) :: comm
  Integer           :: i
  Real( Kind = wp ) :: delr,rdlr,factor1,rrr,dvol,gofr,sum0,sum1

! grid interval for rdf tables

  delr = rusr/Real(mxgusr,wp)
  rdlr = 1.0_wp/delr

! open RDF file and Write headers

  If (comm%idnode == 0) Then
     Open(Unit=nrdfdt, File='USRDAT', Status='replace')
     Write(nrdfdt,'(2a)') '# '//cfgname
     Write(nrdfdt,'(a)')  "# RDF for the two fragments' COMs (umbrella sampling)"
     Write(nrdfdt,'(a,i10,f12.6,i10,e15.6,/)') '# bins, cutoff, frames, volume: ',mxgusr,rusr,ncfusr,volm
     Write(nrdfdt,'(a)') '#'
  End If

! global sum of data on all nodes

  Call gsum(comm,usr(1:mxgusr))

! get normalisation factor

  factor1 = Sum(usr(1:mxgusr))

! running integration of rdf

  sum0 = 0.0_wp
  sum1 = 0.0_wp

! loop over distances

  Do i=1,mxgusr
     gofr = usr(i)/factor1
     sum0 = sum0 + gofr

     rrr  = (Real(i,wp)-0.5_wp)*delr
     dvol = fourpi*delr*(rrr**2+delr**2/12.0_wp)
     gofr = gofr*volm/dvol
     sum1 = sum1 + gofr

! print out information

     If (comm%idnode == 0) Write(nrdfdt,"(1p,4e15.6)") rrr,gofr
!     If (idnode == 0) Write(nrdfdt,"(1p,4e15.6)") rrr,gofr,sum0,sum1
  End Do

  If (comm%idnode == 0) Close(Unit=nrdfdt)

! distribute usr between nodes

  usr(:) = usr(:) / Real(comm%mxnode,wp)

End Subroutine usr_compute

End Module rdfs
