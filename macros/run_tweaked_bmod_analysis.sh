#!/bin/bash
#This script goes through a given run range and first runs the QwMpsOnly code
#and then executes finBModResiduals2.C to produce the residual cosine and
#sine amplitudes.
#The run range is determined by the arguments passed. This program will only
#analyze runs between 10000 and 19999.
START_RUN=$1
LAST_RUN=$2
POST=$3
TWEAK_FACTOR=1.0
COIL_TO_TWEAK=0
FIELD=3
n=0
prev_run=0
cd $MPS_ONLY_ROOTFILES
for i in $( ls mps_only_1*.root | cut -d'_' -f $FIELD );
#for i in $( find ${OUTPUT}/macrocycle_slopes/coil_coeffs_${1}*_ChiSqMin.set0.dat -mtime +3 -mtime -90 -ls| cut -d'_' -f ${FIELD} );
  do
  if [ "$i" -ne $prev_run ] && [ "$i" -ge $START_RUN ] && [ "$i" -le $LAST_RUN ];then
      echo "$n $i"
#      copy_coeff_files $i $BMOD_SRC 
      nice root -b -q ${BMOD_SRC}/macros/GetSlopesAndResidualsFast.C+\(${i},\"${POST}\"\)
#      nice root -b -q ${BMOD_SRC}/macros/GetSlopesAndResiduals.C+\(${i},\"${POST}\",1\)
#      nice root -b -q ${BMOD_SRC}/macros/SixMonitorSlopesAndResiduals.C+\(${i},\"${POST}\",1\)
      ((++n))
  fi
  prev_run=$i
done
echo "$n runs."

