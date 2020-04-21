#!/usr/bin/env bash

DLP_ROOT=/home/drFaustroll/playground/dlpoly/dl-poly-devel
data="${DLP_ROOT}/data"

genTest(){
  test=$1
  iterations="$2"
  name=$(grep "${test/TEST/TEST } " $DLP_ROOT/data/README.txt | cut -d "-" -f 2-) 
  echo "processing $test - $name"
  [[ ! -d $test ]] && tar -xf $data/${test}.tar.xz
  pushd ${test} >> /dev/null
  echo "    <testcase name=\"${test}\">" >> $testf
  echo "      <path archive=\"Yes\">${test}</path>" >> $testf
  echo "      <description>${test} - ${name}</description>" >> $testf

  for iteration in $iterations; do
    for prop in $props; do
      val=$($DLP_ROOT/utils/Scripts/${prop}.py $iteration)
      echo "        <${prop} iteration=\"${iteration}\">${val}</${prop}>" >> $testf
    done
  done
  echo "    </testcase>" >> $testf
  popd >> /dev/null
  rm -rf ${test}
}

props=$(grep "property name" properties.xml | sed 's;<property name=";;g' | sed 's;">;;g' | xargs )
testf="$DLP_ROOT/utils/Testing/all.xml"


echo "<beetest>" > $testf
cat properties.xml >> $testf


echo "  <testsuite name=\"all\">" >> $testf
echo "    <description>All DL_POLY Tests</description>" >> $testf
echo "    <date>$(date)</date>" >> $testf
echo "    <author>Alin Marin Elena</author>" >> $testf
echo "    <outdir>TestsOut</outdir>" >> $testf
   

for i in $(cat tests.dat); do
  genTest $i "0 10 20"
done 

echo "  </testsuite>" >> $testf
echo "</beetest>" >> $testf

sed "s;Scripts;@PYTHON_EXECUTABLE@ @CMAKE_SOURCE_DIR@/utils/Scripts;g"  $testf > ${testf}.cmake
rm $testf
cp ${testf}.cmake $DLP_ROOT/cmake
