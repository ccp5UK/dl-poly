#!/usr/bin/env bash
DLP_ROOT=/home/drFaustroll/playground/dlpoly/dl-poly-alin
data="${DLP_ROOT}/data"
genTest(){
  test=$1
  iterations="$2"
  name=$(grep "${test/TEST/TEST } " $DLP_ROOT/data/README.txt | cut -d "-" -f 2-) 
  echo "processing test $test - $name"
  tar -xf $data/${test}.tar.xz 
  pushd ${test}
  echo "    <testcase name=\"${test}\">" >> ../$testf
  echo "      <path archive=\"Yes\">${test}</path>" >> ../$testf
  echo "      <description>TEST ${test} - ${name}</description>" >> ../$testf
  for iteration in $iterations; do
    for prop in $props; do
      val=$(${DLP_ROOT}/utils/Scripts/${prop}.py $iteration)
      echo "        <${prop} iteration=\"${iteration}\">${val}</${prop}>" >> ../$testf
    done
  done
  echo "    </testcase>" >> ../$testf
  rm -rf ${test}
  popd
}

props=$(grep "property name" properties.xml | sed 's;<property name=";;g' | sed 's;">;;g' | xargs )
testf="one.xml"


#echo "<beetest>" > $testf
#cat properties.xml >> $testf


#echo "  <testsuite name=\"all\">" >> $testf
#echo "    <description>All DL_POLY Tests</description>" >> $testf
#echo "    <date>"$(date)"</date>" >> $testf
#echo "    <author>Alin Marin Elena</author>" >> $testf
#echo "    <outdir>TestsOut</outdir>" >> $testf
   
genTest $1 "0 10 20" 


#echo "  </testsuite>" >> $testf
#echo "</beetest>" >> $testf
