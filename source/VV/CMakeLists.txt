set(dlpoly_VV_SRCS
pseudo_vv.F90
)
add_library(dlpolyvv ${dlpoly_VV_SRCS})
target_link_libraries(dlpolyvv dlpolycore dlpolythermo)
INSTALL_SUBLIB(dlpolyvv)
