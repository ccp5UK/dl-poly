set(dlpoly_VV_SRCS
constraints_rattle.F90
constraints_shake_vv.F90
npt_l0_vv.F90
npt_l1_vv.F90
npt_m0_vv.F90
npt_m1_vv.F90
pmf_rattle.F90
pmf_shake_vv.F90
pseudo_vv.F90
)
add_library(dlpolyvv ${dlpoly_VV_SRCS})
target_link_libraries(dlpolyvv dlpolycore dlpolythermo)
INSTALL_SUBLIB(dlpolyvv)
