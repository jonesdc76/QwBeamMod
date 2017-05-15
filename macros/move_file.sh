#!/bin/bash


cd /net/data2/paschkelab1/pass5c_reduced_slugfiles/run1/
dir="/net/data3/paschkelab4/pass5c_reduced_slugfiles/run1"
n=0
for i in $(ls red*)
do
    file="${dir}/${i}" 
    if [ -f $file ]; then
	echo "$file already exists."
    else
	echo "Moving ${i} to ${dir}"
	mv $i $dir
	((n++))
    fi
done
echo "$n moved"
