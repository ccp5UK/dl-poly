Subroutine ttm_thermal_diffusion (tstep,time,nstep,nsteql,temp,nstbpo,ndump,nstrun,lines,npage)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for iterating thermal diffusion equation for
! electrons using numerical integration, after a.rutherford and
! d.duffy
!
! copyright - daresbury laboratory
! author    - s.l.darazewicz & m.a.seaton july 2012
! contrib -   g.khara may 2016
! contrib -   g.khara, s.t.murphy & m.a.seaton september 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup
  Use configuration
  Use comms_module
  Use ttm_module
  Use ttm_utils
  Use ttm_track_module

  Implicit None

  Integer, Intent( In ) :: ndump,nstbpo,nsteql,nstep,nstrun,lines,npage
  Real ( Kind = wp ), Intent( In ) :: temp,tstep,time

  Real ( Kind = wp ), Allocatable :: eltemp1(:,:,:,:)
  Real ( Kind = wp ) :: fomAx,fomAy,fomAz,mintstep,maxtstep,opttstep,delx2,dely2,delz2
  Real ( Kind = wp ) :: fopttstep,del2av,eltempmax,eltempmin,eltempmean,eltempKe,eltempmaxKe,eltempminKe
  Real ( Kind = wp ) :: actsite, actxm, actxp, actym, actyp, actzm, actzp, alploc
  Integer :: i,j,k,ii,jj,kk,ijk
  Logical :: safe
  Integer :: fail,redtstepmx,redtstep

! Debugging flag

  Logical :: debug1=.false.

! Initialise eltemp1 (electronic temperature grid for next timestep) and timestep sizes

  Allocate (eltemp1(1:numcell,-eltcell(1):eltcell(1),-eltcell(2):eltcell(2),-eltcell(3):eltcell(3)), Stat = fail)
  If (fail>0) Call error(1087)
  eltemp1 = 0.0_wp
  redtstepmx = 1

! deposition stage 1 (initialization):
! nstep-nsteql offsets equilibration time

  If ((nstep-nsteql)==1 .and. (sdepoType>0 .and. (dEdX>zero_plus .or. fluence>zero_plus))) Call depoinit (time)

! determine timestep reduction factor (chosen empirically, acts beyond minimum stability condition)

  fopttstep = 0.25_wp

! determine maximum/minimum spacings and electronic temperatures

  delx2 = delx*delx
  dely2 = dely*dely
  delz2 = delz*delz
  del2av = (delx2*dely2*delz2)/(dely2*delz2+delx2*delz2+delx2*dely2)
  Call eltemp_max (eltempmax)
  Call eltemp_min (eltempmin)

! This section of the code establishes the optimum size of fourier mesh to ensure stability of the electronic
! temperature finite difference solver

  Select Case (KeType)
  Case (0)
! infinite thermal conductivity
    redtstepmx = 1
    opttstep = tstep
  Case (1)
! constant thermal conductivity and non-metal systems
    mintstep = 0.5_wp*del2av/alp(eltempmax)
    maxtstep = 0.5_wp*del2av/alp(eltempmin)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  Case (2)
! Drude-type thermal conductivity
    mintstep = 0.5_wp*del2av*Ce(eltempmax)/KeD(eltempmax,temp)
    maxtstep = 0.5_wp*del2av*Ce(eltempmin)/KeD(eltempmin,temp)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  Case (3)
  Call eltemp_maxKe (temp, eltempmaxKe)
  Call eltemp_minKe (temp, eltempminKe)
! tabulated thermal conductivity
    mintstep = 0.5_wp*del2av*Ce(eltempmax)/Ke(eltempmaxKe)
    maxtstep = 0.5_wp*del2av*Ce(eltempmin)/Ke(eltempminKe)
    opttstep = fopttstep*Min(mintstep,maxtstep)
    redtstepmx = Max(Ceiling(tstep/opttstep),2)
  End Select

! reduce timestep further for deposition stage

  If (KeType>0 .and. trackInit) Then
    Select Case (tdepoType)
    Case (1)
      redtstepmx = Max(50, redtstepmx)
    Case (2)
      redtstepmx = Max(10000, redtstepmx)
    Case Default
      redtstepmx = Max(2, redtstepmx)
    End Select
  End If

  fomAx = tstep/(delx2*Real(redtstepmx,Kind=wp))
  fomAy = tstep/(dely2*Real(redtstepmx,Kind=wp))
  fomAz = tstep/(delz2*Real(redtstepmx,Kind=wp))

! write information to OUTPUT

  If (mod(nstep,nstbpo) == 0 .or. nstep == 1) Then
    If (idnode == 0) Then
      Write(nrite,'(6x,"ttm thermal diffusion timesteps:",2x,"optimal/ps",3x,"actual/ps",5x,"diff/md")')
      Write(nrite,'(38x,es12.4,es12.4,2x,i10)') opttstep, tstep/Real(redtstepmx,Kind=wp), redtstepmx
      If (ttmdyndens) Then
        Write(nrite,'(6x,"active ion temperature cells:",5x,"atom dens.",5x,"no. of active cells")')
        Write(nrite,'(38x,es12.4,14x,i10)') cellrho,acell
      End If
      If(nstep>1 .and. Mod(lines,npage)/=0) Write(nrite,"(1x,130('-'))")
    End If
  End If

! apply boundary conditions

  Call boundaryHalo ()
  Call boundaryCond (bcTypeE, temp)

! print statistics to files: electronic and ionic temperatures
! (note timestep is subtracted by 1, as these are values at
!  beginning of MD timestep)

  Call printElecLatticeStatsToFile('PEAK_E', time, temp, nstep-1, ttmstats)
  Call peakProfilerElec('LATS_E', nstep-1, ttmtraj)

  Call printLatticeStatsToFile(tempion, 'PEAK_I', time, nstep-1, ttmstats)
  Call peakProfiler(tempion, 'LATS_I', nstep-1, ttmtraj)

! debugging option: print electron-phonon and electronic stopping source terms
!                   (normally switched off)

  If (debug1) Then
    Call printLatticeStatsToFile(gsource, 'PEAK_G', time, nstep-1, ttmstats)
    Call peakProfiler(gsource, 'LATS_G', nstep-1, ttmtraj)
    Call printLatticeStatsToFile(asource, 'PEAK_A', time, nstep-1, ttmstats)
    Call peakProfiler(asource, 'LATS_A', nstep-1, ttmtraj)
  End If

  safe=.true.

! determine energy redistribution from deactivated ionic temperature voxels for slab geometry

  If (redistribute) Call redistribute_Te (temp)

! Adaptive timestep

  Do redtstep = 1, redtstepmx

! deposition stage 2 (with boundary conditions)

    If (trackInit) Then
      Call depoevolve(time, tstep, redtstep, redtstepmx)
      Call boundaryCond(bcTypeE, temp)
    End If

! MAIN LOOP
! this portion of the code is the main electronic temperature solver, solving
! for heat diffusion with two source terms within ionic cells, i.e.
!
!   Ce*d(T_e)/dt = d/dx(Ke * d(T_e)/dx) + volume*(asource + gsource*(T_e-T_i))
!
! the partial differential equation is solved using an explicit finite
! difference solver: care is needed in choosing timestep to ensure
! numerical stability

    Select Case (KeType)
    Case (0)
! infinite thermal conductivity case: set all electronic temperatures
! to mean value in active cells, to system temperature in inactive cells
      Call eltemp_mean(eltempmean)
      eltemp1 = eltempmean
      Do ijk=1,numcell
        If(act_ele_cell(ijk,0,0,0)<=zero_plus) eltemp1(ijk,0,0,0) = temp
      End Do

    Case (1)
! constant thermal conductivity or non-metal case 
      If (redistribute) Then
      ! system with cell deactivation/energy redistribution
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1,numcell
                  If (adjust (ijk,ii,jj,kk)) eltemp(ijk,ii,jj,kk) = eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      actsite = act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*act_ele_cell(ijk-1,ii,jj,kk)
                      actxp = actsite*act_ele_cell(ijk+1,ii,jj,kk)
                      actym = actsite*act_ele_cell(ijk-(ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*act_ele_cell(ijk+(ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*act_ele_cell(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*act_ele_cell(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      alploc = alp(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      alploc = alp(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do
      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)
              Do k=1,ntcell(3)
                Do j=1,ntcell(2)
                  Do i=1,ntcell(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    alploc = alp(eltemp(ijk,ii,jj,kk))
                    eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                      fomAx*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAx*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    Case (2)
! Drude-type thermal conductivity case
      If (redistribute) Then
      ! system with cell deactivation/energy redistribution
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1, numcell
                  If (adjust (ijk,ii,jj,kk)) eltemp(ijk,ii,jj,kk) = eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temeperature sites (and active neighbours)
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      actsite = act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*act_ele_cell (ijk-1,ii,jj,kk)
                      actxp = actsite*act_ele_cell (ijk+1,ii,jj,kk)
                      actym = actsite*act_ele_cell (ijk-(ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*act_ele_cell (ijk+(ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*act_ele_cell (ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*act_ele_cell (ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-1,ii,jj,kk)),temp)*&
                        (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+1,ii,jj,kk)),temp)*&
                        (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-1,ii,jj,kk)),temp)*&
                        (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk)) +&
                        fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+1,ii,jj,kk)),temp)*&
                        (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                        fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                        (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do

      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)
              Do k=1,ntcell(3)
                Do j=1,ntcell(2)
                  Do i=1,ntcell(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                      fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-1,ii,jj,kk)),temp)*&
                      (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAx*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+1,ii,jj,kk)),temp)*&
                      (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAy*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))+&
                      fomAz*KeD(0.5_wp*(eltemp(ijk,ii,jj,kk)+eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)),temp)*&
                      (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))/Ce(eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    Case (3)
! tabulated thermal conductivity: uses local ionic or system temperature to calculate value
      If (redistribute) Then
      ! system with cell deactivation/energy redistribution
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)

              If (ii>-2 .and. ii<2 .and. jj>-2 .and. jj<2 .and. kk>-2 .and. kk<2) Then
              ! replace electronic temperatures with values required for energy redistribution
                Do ijk = 1, numcell
                  If (adjust (ijk,ii,jj,kk)) eltemp(ijk,ii,jj,kk) = eltemp_adj(ijk,ii,jj,kk)
                End Do
              ! calculate thermal diffusion only for active ionic temperature sites (and active neighbours)
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      actsite = act_ele_cell (ijk,ii,jj,kk)
                      actxm = actsite*act_ele_cell(ijk-1,ii,jj,kk)
                      actxp = actsite*act_ele_cell(ijk+1,ii,jj,kk)
                      actym = actsite*act_ele_cell(ijk-(ntcell(1)+2),ii,jj,kk)
                      actyp = actsite*act_ele_cell(ijk+(ntcell(1)+2),ii,jj,kk)
                      actzm = actsite*act_ele_cell(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      actzp = actsite*act_ele_cell(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)
                      eltempKe = Merge(tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                      alploc = Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*actxm*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*actxp*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actym*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*actyp*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzm*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*actzp*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              Else
              ! standard thermal diffusion calculation applies for electronic cells away from ionic cells
                Do k=1,ntcell(3)
                  Do j=1,ntcell(2)
                    Do i=1,ntcell(1)
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      ! note that temperature for thermal conductivity is always system
                      ! temperature for electronic cells away from ionic cells
                      alploc = Ke(temp)/Ce(eltemp(ijk,ii,jj,kk))
                      eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                        fomAx*alploc*(eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAx*alploc*(eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAy*alploc*(eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                        fomAz*alploc*(eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                    End Do
                  End Do
                End Do
              End If

            End Do
          End Do
        End Do
      Else
      ! standard thermal diffusion calculation applies when energy redistribution is not applicable
        Do kk=-eltcell(3),eltcell(3)
          Do jj=-eltcell(2),eltcell(2)
            Do ii=-eltcell(1),eltcell(1)
              Do k=1,ntcell(3)
                Do j=1,ntcell(2)
                  Do i=1,ntcell(1)
                    ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                    eltempKe = Merge(tempion(ijk),temp,(ii==0 .and. jj==0 .and. kk==0))
                    eltemp1(ijk,ii,jj,kk) = eltemp(ijk,ii,jj,kk)+&
                      fomAx*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk-1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAx*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk+1,ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk-(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAy*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk+(ntcell(1)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk-(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))+&
                      fomAz*Ke(eltempKe)/Ce(eltemp(ijk,ii,jj,kk))*&
                      (eltemp(ijk+(ntcell(1)+2)*(ntcell(2)+2),ii,jj,kk)-eltemp(ijk,ii,jj,kk))
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
      End If

    End Select

! electron stopping and electron-phonon couplings

    If (oneway) Then
      If (nstep > nstepcpl) Then
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
              If (act_ele_cell(ijk,0,0,0)>zero_plus) Then
                ! e-s coupling term
                eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)+tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*asource(ijk)
                ! e-p coupling term: only use if electronic temperature 
                ! exceeds ionic temperature
                If (l_epcp .and. eltemp(ijk,0,0,0)>tempion(ijk)) Then
                  Select Case (gvar)
                  Case (0,1)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))
                  Case (2)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))*&
                                                                                   Gep(eltemp(ijk,0,0,0))
                  End Select
                End If
              End If
            End Do
          End Do
        End Do
      End If
    Else
      If (nstep > nstepcpl) Then
        Do k=1,ntcell(3)
          Do j=1,ntcell(2)
            Do i=1,ntcell(1)
              ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
              If (act_ele_cell(ijk,0,0,0)>zero_plus) Then
                ! e-s coupling term
                eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)+tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*asource(ijk)
                ! e-p coupling term
                If (l_epcp) Then
                  Select Case (gvar)
                  Case (0,1)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))
                  Case (2)
                    eltemp1(ijk,0,0,0) = eltemp1(ijk,0,0,0)-&
                    tstep*rvolume/(Ce(eltemp(ijk,0,0,0))*Real(redtstepmx,Kind=wp))*gsource(ijk)*(eltemp(ijk,0,0,0)-tempion(ijk))*&
                                                                                   Gep(eltemp(ijk,0,0,0))
                  End Select
                End If
              End If
            End Do
          End Do
        End Do
      End If
    End If

! update electronic temperatures to adjusted values

    Do kk=-eltcell(3),eltcell(3)
      Do jj=-eltcell(2),eltcell(2)
        Do ii=-eltcell(1),eltcell(1)
          Do ijk=1,numcell
            eltemp(ijk,ii,jj,kk) = eltemp1(ijk,ii,jj,kk)
          End Do
        End Do
      End Do
    End Do

! update boundary halo values and apply boundary conditions

    Call boundaryHalo ()
    Call boundaryCond (bcTypeE, temp)

! simple stability check for simulation

    If (Any(eltemp < 0.0_wp)) safe = .false.
    Call gcheck(safe)
    If (.not. safe) Call error (683)

  End Do

! Dumping Te file every ndump steps
  Call ttm_system_revive ('DUMP_E',nstep,time,ndump,nstrun)

  Deallocate (eltemp1, Stat = fail)
  If (fail>0) Call error(1088)

End Subroutine ttm_thermal_diffusion
