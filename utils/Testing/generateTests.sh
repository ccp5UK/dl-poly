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
   

genTest TEST01 "0 10 20"
genTest TEST02 "0 10 20"
genTest TEST03 "0 10 20"
genTest TEST04 "0 10 20"
genTest TEST05 "0 10 20"
genTest TEST06 "0 10 20"
genTest TEST07 "0 10 20"
genTest TEST08 "0 10 20"
genTest TEST09 "0 10 20"
genTest TEST10 "0 10 20"
genTest TEST11 "0 10 20"
genTest TEST12 "0 10 20"
genTest TEST13 "0 10 20"
genTest TEST14 "0 10 20"
genTest TEST15 "0 10 20"
genTest TEST16 "0 10 20"
genTest TEST17 "0 10 20"
genTest TEST18 "0 10 20"
genTest TEST19 "0 10 20"
genTest TEST20 "0 10 20"
genTest TEST21 "0 10 20"
genTest TEST22 "0 10 20"
genTest TEST23 "0 10 20"
genTest TEST24 "0 10 20"
genTest TEST25 "0 10 20"
genTest TEST29 "0 10 20"
genTest TEST30 "0 10 20"
genTest TEST31 "0 10 20"
genTest TEST32 "0 10 20"
genTest TEST33 "0 10 20"
genTest TEST34 "0 10 20"
genTest TEST35 "0 10 20"
genTest TEST36 "0 10 20"
genTest TEST37 "0 10 20"
genTest TEST38 "0 10 20"
genTest TEST39 "0 10 20"
genTest TEST40 "0 10 20"
genTest TEST41 "0 10 20"
genTest TEST42 "0 10 20"
genTest TEST43 "0 10 20"
genTest TEST44 "0 10 20"
genTest TEST45 "0 10 20"
genTest TEST46 "0 10 20"
genTest TEST47 "0 10 20"
genTest TEST48 "0 10 20"
genTest TEST49 "0 10 20"
genTest TEST50 "0 10 20"
genTest TEST51 "0 10 20"
genTest TEST52 "0 10 20"
genTest TEST53 "0 10 20"
genTest TEST54 "0 10 20"
genTest TEST55 "0 10 20"
genTest TEST56 "0 10 20"
genTest TEST57 "0 10 20"
genTest TEST58 "0 10 20"
genTest TEST59 "0 10 20"
genTest TEST60 "0 10 20"
genTest TEST61 "0 10 20"
genTest TEST62 "0 10 20"
genTest TEST63 "0 10 20"
genTest TEST64 "0 10 20"
genTest TEST65 "0 10 20"
genTest TEST66 "0 10 20"
genTest TEST67 "0 10 20"
genTest TEST68 "0 10 20"
genTest TEST69 "0 10 20"
genTest TEST70 "0 10 20"
genTest TEST71 "0 10 20"
genTest TEST72 "0 10 20"
genTest TEST73 "0 10 20"
genTest TEST74 "0 10 20"
genTest TEST75 "0 10 20"
genTest TEST76 "0 10 20"
genTest TEST77 "0 10 20"
genTest TEST78 "0 10 20"
genTest TEST79 "0 10 20"
genTest TEST80 "0 10 20"
genTest TEST81 "0 10 20"
genTest TEST82 "0 10 20"
genTest TEST83 "0 10 20"
genTest TEST84 "0 10 20"
genTest TEST85 "0 10 20"
genTest TEST86 "0 10 20"
genTest TEST87 "0 10 20"
genTest TEST88 "0 10 20"
genTest TEST89 "0 10 20"
genTest TEST90 "0 10 20"
genTest TEST91 "20"
genTest TEST92 "0 10 20"
genTest TEST93 "0 10 20"
genTest TEST94 "0 10 20"
genTest TEST95 "0 10 20"
genTest TEST96 "0 10 20"
genTest TEST97 "0 10 20"
genTest TEST98 "0 10 20"
genTest TEST99 "0 10 20"
genTest TEST100 "0 10 20"
genTest TEST101 "0 10 20"
genTest TEST102 "0 10 20"
genTest TEST103 "0 10 20"
genTest TEST104 "0 10 20"
genTest TEST105 "0 10 20"
genTest TEST106 "0 10 20"
genTest TEST107 "0 10 20"
genTest TEST108 "0 10 20"
genTest TEST109 "0 10 20"
genTest TEST110 "0 10 20"
genTest TEST111 "0 10 20"
genTest TEST112 "0 10 20"
genTest TEST113 "0 10 20"
genTest TEST114 "0 10 20"
genTest TEST115 "0 10 20"
genTest TEST116 "0 10 20"
genTest TEST117 "0 10 20"
genTest TEST118 "0 10 20"
genTest TEST119 "0 10 20"
genTest TEST120 "0 10 20"
genTest TEST121 "0 10 20"
genTest TEST122 "0 10 20"
genTest TEST123 "0 10 20"
genTest TEST124 "0 10 20"
genTest TEST125 "0 10 20"
genTest TEST126 "0 10 20"
genTest TEST127 "0 10 20"
genTest TEST128 "0 10 20"
genTest TEST129 "0 10 20"
genTest TEST130 "0 10 20"
genTest TEST131 "0 10 20"
genTest TEST132 "0 10 20"
genTest TEST133 "0 10 20"
genTest TEST134 "0 10 20"
genTest TEST135 "0 10 20"
genTest TEST136 "0 10 20"
genTest TEST137 "0 10 20"
genTest TEST138 "0 10 20"
genTest TEST139 "0 10 20"
genTest TEST140 "0 10 20"
genTest TEST141 "0 10 20"
genTest TEST144 "0 10 20"
genTest TEST145 "0 10 20"
genTest TEST146 "0 10 20"
genTest TEST147 "0 10 20"
genTest TEST148 "0 10 20"
genTest TEST149 "0 10 20"
genTest TEST150 "0 10 20"
genTest TEST151 "0 10 20"
genTest TEST152 "0 10 20"
genTest TEST153 "0 10 20"
genTest TEST154 "0 10 20"
genTest TEST155 "0 10 20"
genTest TEST156 "0 10 20"
genTest TEST157 "0 10 20"
genTest TEST158 "0 10 20"
genTest TEST159 "0 10 20"
genTest TEST160 "0 10 20"
genTest TEST161 "0 10 20"
genTest TEST162 "0 10 20"
genTest TEST163 "0 10 20"
genTest TEST169 "0 10 20"
genTest TEST170 "0 10 20"
genTest TEST171 "0 10 20"
genTest TEST173 "0 50 100"
genTest TEST174 "0 50 100"
genTest TEST175 "0 50 100"
genTest TEST176 "0 10 20"
genTest TEST177 "0 10 20"
genTest TEST178 "0 10 20"
genTest TEST179 "0 10 20"
genTest TEST181 "0 10 20"


echo "  </testsuite>" >> $testf
echo "</beetest>" >> $testf

sed "s;Scripts;@PYTHON_EXECUTABLE@ @CMAKE_SOURCE_DIR@/utils/Scripts;g"  $testf > ${testf}.cmake
rm $testf
cp ${testf}.cmake $DLP_ROOT/cmake
