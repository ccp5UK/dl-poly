#!/usr/bin/env bash

if [ "$#" -eq "1" ] ; then
  it=$1
else
  echo "specify the iteration"
  exit
fi
echo $(grep -A 8 -w " $it " STATIS | tail -n 8 | tr -d '\n' | awk '{print $15}')
