Module cell_vdw_forces_module2

Use kinds_f90
Use setup_module
Use vdw_module
#ifdef __OPENMP
  Use comms_module, Only: mxthreads
  Use config_module, Only : natms,ltg,ltype,list
#else
  Use config_module, Only : natms,ltg,ltype,list,fxx,fyy,fzz
#endif
Use cell_list_module

  Implicit none
  Private;
  !Externally visible functions are created by preprocessor from cell_vdw_forces.h
  Public :: cell_local_pair_vdw_non_ld_ls_interaction,&
            cell_pair_vdw_non_ld_ls_interaction,&
            cell_self_vdw_non_ld_ls_interaction,&
            cell_local_pair_vdw_non_ld_non_ls_interaction,&
            cell_pair_vdw_non_ld_non_ls_interaction,&
            cell_self_vdw_non_ld_non_ls_interaction,&
            cell_local_pair_vdw_ld_ls_interaction,&
            cell_pair_vdw_ld_ls_interaction,&
            cell_self_vdw_ld_ls_interaction,&
            cell_local_pair_vdw_ld_non_ls_interaction,&
            cell_pair_vdw_ld_non_ls_interaction,&
            cell_self_vdw_ld_non_ls_interaction
  Private :: non_ld_ls_interaction, non_ld_non_ls_interaction, ld_ls_interaction, ld_non_ls_interaction

Contains

Subroutine ld_ls_interaction( irrr, rrr, rdr, eng, gamma, itype, rvdw,k )
Implicit None
Real(Kind = wp), Intent (In)  :: irrr, rrr, rdr, rvdw
Real(Kind = wp), Intent (Out) :: eng, gamma
Integer, Intent(In)           :: itype,k
Real(Kind = wp) :: eps, sig, sor6
  !Interaction to be compute when ld_vdw = true and ls_vdw = true

        If(itype /= 2) then
          !!TODO Implement other interactions
          print *, "NYI:cells_vdw_forces only implements itype=2 interaction atm.", itype
        else
        ! Lennard-Jones potential :: i=4*eps*[(sig/r)^12-(sig/r)^6]
          eps = prmvdw(1,k)
          sig = prmvdw(2,k)
          sor6 = (sig*irrr)**6
  
          eng = 4.0_wp*eps*sor6*(sor6-1.0_wp)
          gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*irrr*irrr
  
          eng = eng + afs(k)*rrr+bfs(k)
          gamma = gamma - afs(k)*irrr
          End If
End Subroutine ld_ls_interaction


Subroutine ld_non_ls_interaction( irrr,rrr, rdr,  eng, gamma, itype, rvdw ,k)
Implicit None
Real(Kind = wp), Intent (In)  :: irrr, rrr, rdr, rvdw
Real(Kind = wp), Intent (Out)  :: eng, gamma
Integer, Intent(In)           :: itype,k
Real(Kind = wp) :: eps, sig, sor6
  !Interaction to be compute when ld_vdw = true and ls_vdw = false

        If(itype /= 2) then
          !!TODO Implement other interactions
          print *, "NYI:cells_vdw_forces only implements itype=2 interaction atm.", itype
        else
        ! Lennard-Jones potential :: i=4*eps*[(sig/r)^12-(sig/r)^6]
          eps = prmvdw(1,k)
          sig = prmvdw(2,k)
          sor6 = (sig*irrr)**6
  
          eng = 4.0_wp*eps*sor6*(sor6-1.0_wp)
          gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*irrr*irrr
  
        End If

End Subroutine ld_non_ls_interaction

Subroutine non_ld_ls_interaction(irrr, rrr, rdr, eng, gamma, itype,rvdw,k )
  Implicit None
  Real(Kind = wp), Intent (In)  :: irrr, rrr, rdr,rvdw
  Real(Kind = wp), Intent(Out)  :: eng, gamma 
  Integer, Intent(In)           :: itype,k

  Integer  :: l 
  Real(Kind = wp) :: t1,t2,vk,vk1,vk2,gk,gk1,gk2,ppp
  !Interaction to be compute when ld_vdw = false and ls_vdw = true

           l   = Int(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)
! calculate interaction energy using 3-point interpolation
              vk  = vvdw(l,k)
              vk1 = vvdw(l+1,k)
              vk2 = vvdw(l+2,k)

              t1 = vk  + (vk1 - vk )*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              eng = t1 + (t2-t1)*ppp*0.5_wp
              eng = eng + gvdw(mxgvdw-4,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgvdw-4,k) ! force-shifting
! calculate forces using 3-point interpolation

           gk  = gvdw(l,k) ; If (l == 0) gk = gk*rrr
           gk1 = gvdw(l+1,k)
           gk2 = gvdw(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)*(irrr*irrr)
           gamma = gamma - gvdw(mxgvdw-4,k)*(irrr*rvdw) !force-shifting

End Subroutine non_ld_ls_interaction


Subroutine non_ld_non_ls_interaction(irrr, rrr, rdr, eng, gamma, itype,rvdw,k )
  Implicit None
  Real(Kind = wp), Intent (In)  :: irrr,rrr, rdr,rvdw
  Real(Kind = wp), Intent(Out)  :: eng, gamma 
  Integer, Intent(In)           :: itype,k

  Integer  :: l 
  Real(Kind = wp) :: t1,t2,vk,vk1,vk2,gk,gk1,gk2,ppp
  !Interaction to be compute when ld_vdw = false and ls_vdw = false

           l   = Int(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)
! calculate interaction energy using 3-point interpolation
              vk  = vvdw(l,k)
              vk1 = vvdw(l+1,k)
              vk2 = vvdw(l+2,k)

              t1 = vk  + (vk1 - vk )*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              eng = t1 + (t2-t1)*ppp*0.5_wp
! calculate forces using 3-point interpolation

           gk  = gvdw(l,k) ; If (l == 0) gk = gk*rrr
           gk1 = gvdw(l+1,k)
           gk2 = gvdw(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)*(irrr*irrr)

End Subroutine non_ld_non_ls_interaction

#define FUNC ld_ls_interaction
#include "cell_vdw_forces.h"
#undef FUNC

#define FUNC non_ld_ls_interaction
#include "cell_vdw_forces.h"
#undef FUNC

#define FUNC ld_non_ls_interaction
#include "cell_vdw_forces.h"
#undef FUNC

#define FUNC non_ld_non_ls_interaction
#include "cell_vdw_forces.h"
#undef FUNC

End Module cell_vdw_forces_module2
