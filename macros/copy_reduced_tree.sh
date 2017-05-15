#!/bin/bash
#stem can be "", "set11_reg_" or "std_reg_"
slug_start=$1
slug_end=$2
stem=$3
reduced=""
reduced_tree=""

cd /net/data2/paschkelab2/missing_quartets/

for i in $( ls | cut -d'g' -f 2 | cut -d'_' -f 1);
  do
    if [ "$3" != "" ]; then
      reduced="reduced_slug${i}.root"
    fi
    cd ${BMOD_SRC}/macros
    if [ $i -ge $slug_start ] && [ $i -le $slug_end ];then
	echo "$i"
	if [ $i -ne 100 ] && [ $i -ne 101 ];then
	    nice root -b -q ${BMOD_SRC}/macros/copyReducedTree.C+\(\"${stem}reduced_slug${i}.root\",\"slug\",\"${reduced}\",${i}\)
	fi
    fi
  done
