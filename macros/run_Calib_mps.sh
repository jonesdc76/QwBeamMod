#!/bin/bash
stem=$3

cd $MPS_ONLY_ROOTFILES
for i in $( ls mps_only_f*.root | cut -d'_' -f 4 | cut -d'.' -f 1);
  do
    if [ "$i" -ge $1 ] && [ "$i" -le $2 ];then
	echo "Analyzing run $i"
	root -b -q $MPS_ONLY_ROOTFILES/runCalib_mps.C+\(${i},\"${3}\"\)
	mv Calib_mps.root Calib_mps_${i}${stem}.root 
    fi
  done

