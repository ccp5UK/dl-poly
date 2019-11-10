#%Module1.0

proc ModulesHelp { } {
  puts stderr "Simple module for DL_POLY @DLPOLY_VERSION@"
}
module-whatis   "Simple module for DL_POLY @DLPOLY_VERSION@"

prepend-path PATH @CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_BINDIR@
prepend-path LD_LIBRARY_PATH @CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_LIBDIR@

