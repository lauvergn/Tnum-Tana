#!/bin/bash

name_dep=dependencies.mk
SRCFile=fortranlist.mk

list=`ls SRC/Source_PrimOperator/*.f90 SRC/Source_TnumTana_Coord/*.f90 \
     SRC/Source_TnumTana_Coord/Qtransfo*/*.f90 SRC/Source_TnumTana_Coord/Tana/*.f90 SRC/Source_TnumTana_Coord/Tnum/*.f90`
#ExcludeList='Main_TnumTana_FDriver.f90 TEST_TnumTana.f90 Tnum90_MCTDH.f90 Tnum90_MidasCpp.f90 Tnum90.f90 Tnum_OOP.f90'
ExcludeList=' OneDTransfo_m.f90 '

echo "#===============================================" > $name_dep
echo "#===============================================" > $SRCFile
echo "SRCFILES= \\" >> $SRCFile

for ff90 in $list
do
   ff=`awk '{name=$1 ; n=split(name,tab,"/") ; if (n > 0) {l=length(tab[n]) ; print tab[n]}}' <<< $ff90`
   #echo $ff
   if grep -vq $ff <<< $ExcludeList;  then
     echo $ff " \\" >> $SRCFile
     awk -f scripts/mod2file.awk $ff90 >> $name_dep
   fi
done
echo "#===============================================" >> $name_dep
for ff90 in $list
do
   ff=`awk '{name=$1 ; n=split(name,tab,"/") ; if (n > 0) {l=length(tab[n]) ; print tab[n]}}' <<< $ff90`
   if grep -vq $ff <<< $ExcludeList;  then
     awk -f scripts/dep2.awk $ff90 >> $name_dep
   fi
done