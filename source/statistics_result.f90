Subroutine statistics_result                       &
           (rcut,lmin,lpana,lrdf,lprdf,lzdn,lpzdn, &
           nstrun,keyens,keyshl,iso,               &
           temp,press,strext,nstep,tstep,time,tmst)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing simulation summary
!
! copyright - daresbury laboratory
! author    - w.smith december 1992
! amended   - i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,gtime
  Use setup_module
  Use site_module,       Only : ntpatm,unqatm,numtyp,dens
  Use config_module,     Only : cell,volm
  Use vnl_module,        Only : llvnl,skipvnl
  Use core_shell_module, Only : passshl
  Use minimise_module,   Only : passmin
  Use bonds_module,      Only : ncfbnd
  Use angles_module,     Only : ncfang
  Use dihedrals_module,  Only : ncfdih
  Use inversions_module, Only : ncfinv
  Use rdf_module,        Only : ncfrdf
  Use z_density_module,  Only : ncfzdn
  Use statistics_module
  Use msd_module

  Implicit None

  Logical,           Intent( In    ) :: lmin,lpana,lrdf,lprdf,lzdn,lpzdn
  Integer,           Intent( In    ) :: nstrun,keyens,keyshl,iso,nstep
  Real( Kind = wp ), Intent( In    ) :: rcut,temp,press,strext(1:9),tstep,time,tmst

  Logical           :: check
  Integer           :: i,j,iadd
  Real( Kind = wp ) :: avvol,avcel(1:9),dc,srmsd,timelp,tmp,h_z,tx,ty

! VNL skipping statistics

  If (llvnl .and. idnode == 0) Write(nrite,"(//,            &
     & ' VNL skipping run statistics: average skips', f7.2, &
     & ' minimum skips ', i5, ' maximum skips ', i5)")      &
     skipvnl(3),Nint(Merge(skipvnl(4),skipvnl(5),skipvnl(4)<skipvnl(5))),Nint(skipvnl(5))

! minimisation convergence statistics

  If (lmin .and. idnode == 0) Write(nrite,"(//,              &
     & ' minimisation run statistics: average cycles', f7.2, &
     & ' minimum cycles ', i5, ' maximum cycles ', i5)")     &
     passmin(3),Nint(passmin(4)),Nint(passmin(5))

! shell relaxation convergence statistics

  If (keyshl == 2 .and. idnode == 0) Write(nrite,"(//,           &
     & ' shell relaxation run statistics: average cycles', f7.2, &
     & ' minimum cycles ', i5, ' maximum cycles ', i5)")         &
     passshl(3),Nint(passshl(4)),Nint(passshl(5))

! Get elapsed time

  Call gtime(timelp)

! Get simulation time for averages

  If (numacc == 0) Then
     tmp=0.0_wp
  Else
     tmp=time-tmst
  End If

! Report termination

  If (idnode == 0) &
     Write(nrite,"(/,/,1x,'run terminated after',i9,' steps (',f10.3,   &
          & ' ps), final averages calculated over',i9,' steps (',f10.3, &
          & ' ps).',/,/)") nstep,time,numacc,tmp


! safe average volume and cell

  avvol = volm
  avcel = cell

! If dry/static/minimisation run - NO AVERAGES, but possible RDF and Z-Density

  If (nstep == 0 .and. nstrun == 0) Go To 10

! If still running in the pure equilibration regeime - NO AVERAGES

  If (numacc == 0) Go To 20

! shift back statistical averages as from statistics_collect

  Do i=1,mxnstk
     sumval(i)=sumval(i)+stpvl0(i)
  End Do

! calculate final fluctuations

  Do i=1,mxnstk
     ssqval(i)=Sqrt(ssqval(i))
  End Do

! average volume

  avvol = sumval(19)

! final averages and fluctuations

  If (idnode == 0) Then
     Write(nrite,"(1x,130('-'),                                         &
          &  /,/,1x,'   steps',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg', &
          &  5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',5x,    &
          &  'eng_dih',5x,'eng_tet',/,1x,'time(ps)',5x,' eng_pv',4x,    &
          &  'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',5x,      &
          &  'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,        &
          &  1x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',      &
          &  5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',5x,'vir_pmf', &
          &  7x,'press',/,/,1x,130('-'))")

     Write(nrite,'(1x,i8,1p,9e12.4,/,0p,f9.5,1p,9e12.4,/,1x,0p,f8.1,    &
          & 1p,9e12.4)') numacc,(sumval(i),i=1,9),tmp,                  &
                      (sumval(i),i=10,18),timelp,(sumval(i),i=19,27)

     Write(nrite,"(/,1x,' r.m.s. ',1p,9e12.4,/,1x,'fluctn. ',1p,9e12.4, &
          & /,9x,9e12.4)") (ssqval(i),i=1,27)

     Write(nrite,"(1x,130('-'))")
  End If

! Move at the end of the default 27 quantities

  iadd = 27

  If (l_msd) iadd = iadd+2*mxatdm

! Write out estimated diffusion coefficients

  If (idnode == 0) Then
     Write(nrite,"(/,/,12x,a)") 'Approximate 3D Diffusion Coefficients and square root of MSDs'
     Write(nrite,"(/,12x,'atom',9x,'DC (10^-9 m^2 s^-1)',3x,'Sqrt[MSD] (Ang)',/)")
  End If

  Do i=1,ntpatm
     If (numtyp(i) > zero_plus) Then
        dc = 10.0_wp * (ravval(iadd+i)-sumval(iadd+i)) / &
             (3.0_wp*Real(numacc-Min(mxnstk,numacc-1),wp)*tstep)
        If (dc < 1.0e-10_wp) dc = 0.0_wp

        srmsd = Sqrt(ravval(iadd+i))
        If (idnode == 0) Write(nrite,'(12x,a8,1p,2(7x,e13.4))') unqatm(i),dc,srmsd
     End If
  End Do

  iadd = iadd+ntpatm

! print out average pressure tensor

  If (idnode == 0) Then
     Write(nrite,"(/,/,16x,'Average pressure tensor  (katms)',30x,'r.m.s. fluctuations',/)")

     Do i=iadd,iadd+6,3
        Write(nrite,'(9x,1p,3e12.4,24x,3e12.4)') (sumval(i+j),j = 1,3),(ssqval(i+j),j = 1,3)
     End Do

     Write(nrite,'(/,12x,a,1p,e12.4)') 'trace/3  ', (sumval(iadd+1)+sumval(iadd+5)+sumval(iadd+9))/3.0_wp
  End If

  iadd = iadd+9

! Write out mean cell vectors for npt/nst

  If (keyens >= 20) Then

! average cell (again)

     Do i=1,9
        avcel(i) = sumval(iadd+i)
     End Do

     If (idnode == 0) Then
        Write(nrite,"(/,/,16x,'Average cell vectors     (Angs) ',30x,'r.m.s. fluctuations ',/)")

        Do i=iadd,iadd+6,3
           Write(nrite,'(3f20.10,9x,1p,3e12.4)') (sumval(i+j),j = 1,3),(ssqval(i+j),j = 1,3)
        End Do
     End If

     iadd = iadd+9

     If (iso > 0) Then
        h_z=sumval(iadd+1)

        If (idnode == 0) Then
           Write(nrite,"(/,/,16x,'Average surface area, fluctuations & mean estimate (Angs^2)',/)")
           Write(nrite,'(1p,3e12.4)') sumval(iadd+2),ssqval(iadd+2),avvol/h_z
        End If

        iadd = iadd+2

        If (iso > 1) Then
           tx= -h_z * ( sumval(iadd-9-8-2)/prsunt - (press+strext(1)) ) * tenunt
           ty= -h_z * ( sumval(iadd-9-7-2)/prsunt - (press+strext(5)) ) * tenunt
           If (idnode == 0) Then
              Write(nrite,"(/,16x,'Average surface tension, fluctuations & mean estimate in x (dyn/cm)',/)")
              Write(nrite,'(1p,3e12.4)') sumval(iadd+1),ssqval(iadd+1),tx
              Write(nrite,"(/,16x,'Average surface tension, fluctuations & mean estimate in y (dyn/cm)',/)")
              Write(nrite,'(1p,3e12.4)') sumval(iadd+2),ssqval(iadd+2),ty
           End If

           iadd = iadd+2
        End If
     End If

  End If

! Write out remaining registers

  check = .false.
  Do i=iadd+1,mxnstk
     If (Abs(sumval(i)) > zero_plus .or. Abs(ssqval(i)) > zero_plus) check=.true.
  End Do

  If (check .and. idnode == 0) Then
     Write(nrite,"(/,/,12x,'Remaining non-zero statistics registers ', &
          & /,/,12x,'Register',7x,'Average value',8x,'r.m.s. fluc.')")

     Do i=iadd+1,mxnstk
        If (Abs(sumval(i)) > zero_plus .or. Abs(ssqval(i)) > zero_plus) &
           Write(nrite,'(10x,i10,2f20.10)') i,sumval(i),ssqval(i)
     End Do
  End If

10 Continue

! scale densities for average volume and average volume and cell

  Do i=1,ntpatm
     If (numtyp(i) > zero_plus) dens(i)=dens(i)*(volm/avvol)
  End Do

  volm = avvol
  cell = avcel

! calculate and print radial distribution functions

  If (lrdf .and. lprdf .and. ncfrdf > 0) Call rdf_compute(rcut)

! calculate and print z-density profile

  If (lzdn .and. lpzdn .and. ncfzdn > 0) Call z_density_compute()

! Calculate and print PDFs

  If (lpana) Then
     If (mxgbnd1 > 0 .and. ncfbnd > 0) Call bonds_compute(temp)
     If (mxgang1 > 0 .and. ncfang > 0) Call angles_compute(temp)
     If (mxgdih1 > 0 .and. ncfdih > 0) Call dihedrals_compute(temp)
     If (mxginv1 > 0 .and. ncfinv > 0) Call inversions_compute(temp)
  End If

20 Continue

! print final time check

  Call gtime(timelp)

  If (idnode == 0) Write(nrite,'(/,/,/,1x,"time elapsed since job start: ", f12.3, " sec",/)') timelp

End Subroutine statistics_result
