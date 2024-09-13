#!/bin/bash

dev=$1

for lib in QDUtilLib AD_dnSVM EVRT_dnSVM nDindex 
do
  rm -f $lib # remove the link
  rm -rf $lib* # remove the directory
  get_lib.sh $lib $dev
done
