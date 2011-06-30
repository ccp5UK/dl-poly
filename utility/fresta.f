      program fresta

c*********************************************************************
c
c     dl_poly utility to calculate the average energy differences
c     that form part of a thermodynamic integration
c
c     copyright - daresbury laboratory
c     author    -   w.smith nov 2008
c
c*********************************************************************

      implicit none

      character header(80),units(40)
      integer nstep,numacc
      real(8) engcfg,engkin,avgcfg,avgkin,av2cfg,av2kin
      real(8) sclnv1,sclnv2,pfree,lambda1,lambda2,dlambda

c     open main input FREENG file

      open(7,file="FREENG")

c     read header data

      read(7,'(80a1)')header
      write(*,'(80a1)')header
      read(7,'(40a1)')units
      write(*,'(40a1)')units
      read(7,'(4e16.8)')pfree,lambda1,lambda2,dlambda
      write(*,'("pfree,lambda1,lambda2,dlambda:",/,1p,4e16.8)')
     x  pfree,lambda1,lambda2,dlambda

c     set statistical variables

      numacc=0
      avgcfg=0.d0
      avgkin=0.d0
      av2cfg=0.d0
      av2kin=0.d0

c     read periodic data

      do while(.true.)

        read(7,'(i10,2e16.8)',end=100)nstep,engcfg,engkin

c     statistical counters

        numacc=numacc+1
        sclnv2=1.d0/dble(numacc)
        sclnv1=dble(numacc-1)/dble(numacc)

c     calculate averages

        av2cfg=sclnv1*(av2cfg+sclnv2*(engcfg-avgcfg)**2)
        avgcfg=sclnv1*avgcfg+sclnv2*engcfg
        av2kin=sclnv1*(av2kin+sclnv2*(engkin-avgkin)**2)
        avgkin=sclnv1*avgkin+sclnv2*engkin

      enddo

c     end of FREENG file

  100 close (7)

      write(*,'("Finished reading FREENG file")')
      write(*,'("number of data points sampled",i10)')numacc
      write(*,'("dlambda, average energies and RMS deviations:")')
      write(*,'(1p,5e16.8)')dlambda,avgcfg,avgkin,sqrt(av2cfg),
     x  sqrt(av2kin)

c     open main output file FREEDATA

      open(8,file="FREEDATA")

      write(8,'(80a1)')header
      write(8,'(40a1)')units
      write(8,'("number of data points sampled",i10)')numacc

c     write average energies and RMS deviations

      write(8,'("dlambda, average energies and RMS deviations:")')
      write(8,'(1p,5e16.8)')dlambda,avgcfg,avgkin,sqrt(av2cfg),
     x  sqrt(av2kin)

      write(8,'("END")')

c     close FREEDATA file

      close(8)
      write(*,'("Finished writing FREEDATA file")')
      write(*,'("Job done")')

      end
