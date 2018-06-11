Module link_cells
  Use kinds,         Only : wp, li
  Use comms,         Only : comms_type,gcheck,gmax,gsum
  Use setup
  Use domains,       Only : idx,idy,idz, nprx,npry,nprz, &
                            r_nprx,r_npry,r_nprz, &
                            nprx_r,npry_r,nprz_r
  Use configuration, Only : cell,natms,nlast,ltg,lfrzn, &
                            xxx,yyy,zzz,neigh%list_excl,neigh%list
  Use core_shell,    Only : listshl,legshl
  Use mpole,         Only : keyind,lchatm
  Use development,   Only : development_type
  Use errors_warnings, Only : error,warning,info
  Use numerics, Only : dcell, invert,match
  Use timer,  Only : timer_type,start_timer,stop_timer
  Use neighbours, Only : neighbours_type
  Implicit None

  Private
Contains

End Module link_cells
