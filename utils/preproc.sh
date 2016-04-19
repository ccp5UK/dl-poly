#!/usr/bin/env bash

echo "This script is for developers... "
echo "shall be used if any of the input files from here are edited "
echo "The script is the manual"


for i in plumed netcdf kim; do
cat > ${i}_modul~.F90 <<EOF
! This file is generated manually from ${i}_module.F90 by
! gfortran -E -P ${i}_module.F90 > ${i}_modul~.F90
EOF
gfortran -E -P ${i}_module.F90 >> ${i}_modul~.F90
cat > ${i}_module_pre.F90 <<EOF
! This file is generated manually from ${i}_module.F90 by
! gfortran -E -P -D${i^^} ${i}_module.F90 > ${i}_module_pre.F90
EOF
if [ "x$i" != "xkim" ] ; then
  gfortran -E -P -D${i^^} ${i}_module.F90 >> ${i}_module_pre.F90
else
  sed -i "s;#include;!#include;g" ${i}_module.F90
  sed -i "s;#define;!#define;g" ${i}_module.F90
  gfortran -E -P -D${i^^} ${i}_module.F90 >> ${i}_module_pre.F90
  sed -i "s;!#include;#include;g" ${i}_module.F90
  sed -i "s;!#define;#define;g" ${i}_module.F90
  sed -i "s;!#include;#include;g" ${i}_module_pre.F90
  sed -i "s;!#define;#define;g" ${i}_module_pre.F90
fi
done
