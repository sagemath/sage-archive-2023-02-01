#!/usr/bin/env bash
type -atp java
OUT=$?
#echo $OUT
if [ $OUT -eq 0 ]; then
#found java now check 1.5<=version<=1.7
#version >1.7 may be OK just not checked yet.
  java -version 2>&1|grep version.*[1]\.[567]
else
  exit 1
fi
#OUT=$?
#echo $OUT
