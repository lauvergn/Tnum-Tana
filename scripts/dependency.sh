#!/bin/bash

name_dep=scripts/dependencies.mk
SRCFile=scripts/fortranlist.mk

list=`find SRC  -name "*.f90"`
ExcludeList=' read_para.f90 sub_system.f90 calc_f2_f1Q.f90 Calc_Tab_dnQflex.f90 Sub_X_TO_Q_ana.f90 OneDTransfo_m.f90 '

echo "#===============================================" > $name_dep
echo "#===============================================" > $SRCFile
echo "SRCFILES := \\" >> $SRCFile

for ff90 in $list
do
   ff=`basename $ff90`
   #echo $ff
   if grep -vq $ff <<< $ExcludeList;  then
     echo $ff " \\" >> $SRCFile
     awk -f scripts/mod2file.awk $ff90 >> $name_dep
   fi

done

echo "#===============================================" >> $name_dep
for ff90 in $list
do
   ff=`basename $ff90`
   #echo $ff
   #echo '# '$ff >> $name_dep
   if grep -vq $ff <<< $ExcludeList;  then
     mname=`awk -f scripts/get_modname.awk $ff90`
     #echo "#mname: '"$mname"'"
     #echo "#mname: " $mname >> $name_dep
     test -z $mname || awk -v mod_name=$mname -f scripts/dep2.awk $ff90 >> $name_dep
   fi
done