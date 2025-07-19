#!/bin/bash

dev=$1

for lib in QDUtilLib AD_dnSVM ConstPhys QuantumModelLib nDindex EVRT_dnSVM FOR_EVRT
do
  rm -fr $lib # remove the link
  rm -rf $lib* # remove the directory
  get_lib.sh $lib $dev
done
