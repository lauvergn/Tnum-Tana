#!/bin/bash

EXTLIB_TYPE=$1

BaseName=nDindex

echo "In get_"$BaseName".sh"


SAVE_version="Save_"$BaseName"_devloc"
LOC_version=$BaseName

if (test -d $LOC_version) then
  echo $LOC_version directory exist
  exit 0
fi

rm -rf $BaseName* #always remove the link


#latest release
#latest HEAD version (dev version)
#version=https://github.com/lauvergn/nDindex/archive/refs/tags/v4.0.zip


test -z $EXTLIB_TYPE       &&    curl -LJ $version --output $LOC_version.zip
test $EXTLIB_TYPE != 'loc' &&    curl -LJ $version --output $LOC_version.zip

test -e $LOC_version.zip && echo $LOC_version.zip file exist || cp $SAVE_version.zip $LOC_version.zip

unzip $LOC_version.zip
rm -f $LOC_version.zip


LIBDIR=`ls -d $BaseName*`
#echo $LIBDIR

ln -s $LIBDIR $LOC_version

echo "End get_"$BaseName".sh"