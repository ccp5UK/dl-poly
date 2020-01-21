#!/usr/bin/env bash
DLP_ROOT=/home/drFaustroll/playground/dlpoly/dl-poly-alin
data="${DLP_ROOT}/data"
genTest(){
  test=$1
  iterations="$2"
  name=$(grep "TEST $test" README.txt | awk -F \+ '{print $2}') 
  echo "processing test $test - $name"
  tar -xf $data/TEST${test}.tar.xz 
  pushd TEST${test}
  echo "    <testcase name=\"TEST${test}\">" >> ../$testf
  echo "      <path archive=\"Yes\">TEST${test}</path>" >> ../$testf
  echo "      <description>TEST ${test} - ${name}</description>" >> ../$testf
  for iteration in $iterations; do
    for prop in $props; do
      val=$(${DLP_ROOT}/utils/Scripts/${prop}.py $iteration)
      echo "        <${prop} iteration=\"${iteration}\">${val}</${prop}>" >> ../$testf
    done 
  done
  popd
  rm -rf TEST${test}
  echo "    </testcase>" >> $testf
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
   
genTest $1 "0 10 20 30" 


#echo "  </testsuite>" >> $testf
#echo "</beetest>" >> $testf
