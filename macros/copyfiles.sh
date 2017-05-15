#!/bin/bash
sdir="/net/data2/paschkelab1/pass5c_reduced_slugfiles/run1/"
ddir="/net/data3/paschkelab4/pass5c_reduced_slugfiles/run1/"

cd $sdir

for i in $(ls red*.root);
  do
    if [ -e "${ddir}/${i}" ];then
	echo "File $i already exists in destination directory."
    else
	echo "Copying $i to new location."
	cp $i $ddir
    fi
  done
