#!/usr/bin/env bash

DLROOT=/home/alin/dl-poly-internal
data=new-data

genTest(){
  test=$1
  iterations="$2"
  name=$(grep "TEST $test " $DLROOT/data/README.txt | awk -F \- '{print $2}') 
  echo "processing test $test - $name"
  pushd $data/TEST${test} >> /dev/null
  echo "    <testcase name=\"TEST${test}\">" >> $testf
  echo "      <path archive=\"Yes\">TEST${test}</path>" >> $testf
  echo "      <description>TEST ${test} - ${name}</description>" >> $testf

  for iteration in $iterations; do
    for prop in $props; do
      val=$($DLROOT/utils/Scripts/${prop}.py $iteration)
      echo "        <${prop} iteration=\"${iteration}\">${val}</${prop}>" >> $testf
    done 
  done
  rm -rf TEST${test}
  popd >> /dev/null
  echo "    </testcase>" >> $testf
}

props=$(grep "property name" properties.xml | sed 's;<property name=";;g' | sed 's;">;;g' | xargs )
testf="$DLROOT/utils/Testing/all.xml"


echo "<beetest>" > $testf
cat properties.xml >> $testf


echo "  <testsuite name=\"all\">" >> $testf
echo "    <description>All DL_POLY Tests</description>" >> $testf
echo "    <date>$(date)</date>" >> $testf
echo "    <author>Alin Marin Elena</author>" >> $testf
echo "    <outdir>TestsOut</outdir>" >> $testf
   
genTest 01 "0 10 20" 
genTest 02 "0 10 20" 
genTest 03 "0 10 20" 
genTest 04 "0 10 20" 
genTest 05 "0 10 20" 
genTest 06 "0 10 20" 
genTest 07 "0 10 20" 
genTest 08 "0 10 20" 
genTest 09 "0 10 20" 
genTest 10 "0 10 20" 
genTest 11 "0 10 20" 
genTest 12 "0 10 20" 
genTest 13 "0 10 20" 
genTest 14 "0 10 20" 
genTest 15 "0 10 20" 
genTest 16 "0 10 20" 
genTest 17 "0 10 20" 
genTest 18 "0 10 20" 
genTest 19 "0 10 20" 
genTest 20 "0 10 25" 
genTest 21 "0 10 20" 
genTest 22 "0 10 20" 
genTest 23 "0 10 20" 
genTest 24 "0 10 20" 
genTest 25 "0 10 20" 
#genTest 28 "0 10 20" 
echo "  </testsuite>" >> $testf
echo "</beetest>" >> $testf

sed "s;Scripts;@PYTHON_EXECUTABLE@ @CMAKE_SOURCE_DIR@/utils/Scripts;g"  $testf > ${testf}.cmake
rm $testf
cp ${testf}.cmake $DLROOT/cmake
