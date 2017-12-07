#!/usr/bin/env bash

kim_version=1.7.2
kim_compiler=GCC #GCC or INTEL not too much variation 
kim_os=linux # 'linux', 'darwin', or 'freebsd' 
kim_archive=kim-api-v${kim_version}.tgz
kim_folder=${kim_archive/.tgz/}
kim_inst=/opt/kim/${kim_compiler,,}/$kim_version
FFLAGS="-O3 -mtune=native"
CFLAGS="-O3 -mtune=native"
CXXFLAGS="-O3 -mtune=native"
[[ $(uname -m) == "x86_64"  ]] && kim_arch="64bit" || kim_arch=32bit
[[ ! -f ${kim_archive}  ]] && curl -O https://s3.openkim.org/kim-api/${kim_archive}
rm -rf $kim_folder
tar -xf $kim_archive
pushd $kim_folder
  sed -e "s|^KIM_DIR.*$|KIM_DIR = $(pwd)|" \
    -e "s|^#prefix.*$|prefix = ${kim_inst}|" \
    -e "s|^KIM_COMPILERSUITE.*$|KIM_COMPILERSUITE = ${kim_compiler}|" \
    -e "s|^KIM_SYSTEMARCH.*$|KIM_SYSTEMARCH = ${kim_arch}|" \
    -e "s|^KIM_SYSTEMLINKER.*$|KIM_SYSTEMLINKER = ${kim_os}|" \
      Makefile.KIM_Config.example > Makefile.KIM_Config

  make add-Pair_Lennard_Jones_Truncated_Nguyen_Ar__MO_398194508715_000
  make
  make install
  make install-set-default-to-v1
popd
mkdir -p /opt/modules/kim/${kim_compiler,,}
cat > /opt/modules/kim/${kim_compiler,,}/${kim_version} << EOF
#%Module
set PKG kim
set VER ${kim_version}
set COMP ${kim_compiler,,}
proc ModulesHelp { } {
        global dotversion
        puts stderr "\tLoads the \$COMP \$PKG \$VER  Environment"
}

module-whatis  "Loads the \$COMP \$PKG \$VER Environment."
conflict \$PKG
set PREFIX /opt/\$PKG/\$COMP/\$VER
setenv KIM_ROOT \$PREFIX
prepend-path PATH \$PREFIX/bin
prepend-path INCLUDE \$PREFIX/include/kim-api
prepend-path C_INCLUDE_PATH \$PREFIX/include/kim-api
prepend-path CPLUS_INCLUDE_PATH \$PREFIX/include/kim-api
prepend-path CPATH \$PREFIX/include/kim-api
prepend-path MANPATH \$PREFIX/share/man
prepend-path LD_LIBRARY_PATH \$PREFIX/lib
prepend-path LIBRARY_PATH \$PREFIX/lib
EOF
